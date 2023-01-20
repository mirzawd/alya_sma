!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Algebraic_Solver
!> @{
!> @file    presol.f90
!> @date    22/10/2013
!> @author  Guillaume Houzeaux
!> @brief   Operations on matrix and solvers
!> @details Perform some modifications on the matrix and
!>          or RHS to account for:
!>          \verbatim
!>
!>          ---------------------------------------------------
!>          ITASK     Periodicity             Dirichlet
!>          ---------------------------------------------------
!>                    amatr                   amatr
!>                      +        rhsid          +        rhsid
!>                    rhsid                   rhsid
!>          ---------------------------------------------------
!>            1        X           X            X          X
!>            2        -           X            -          -
!>            3        X           X            -          -
!>            4        -           X            -          X
!>          ---------------------------------------------------
!>
!>         \endverbatim
!> @}
!-----------------------------------------------------------------------
subroutine presol(itask,ndof1,ndof2,rhsid,amatr,unkno)
  use def_kintyp,         only : ip,rp,lg
  use def_master,         only : IMASTER
  use def_master,         only : zeror
!  use def_master,         only : current_zone
  use def_solver,         only : solve_sol
  use def_domain,         only : r_dom,c_dom,npoin
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_memory,         only : memory_alloca,memory_deallo

  use def_master, only : ID_TEMPER
  implicit none
  integer(ip), intent(in)    :: itask
  integer(ip), intent(in)    :: ndof1
  integer(ip), intent(in)    :: ndof2
  real(rp),    intent(inout) :: rhsid(ndof1,*)
  real(rp),    intent(inout) :: amatr(ndof2,ndof1,*)
  real(rp),    intent(inout) :: unkno(ndof1,*)
  integer(ip)                :: izdom,jzdom
  integer(ip)                :: idofn,ipoin,jpoin,kpoin
  integer(ip)                :: jdofn,izdod,ndofn
  integer(ip)                :: num_blocks
!  integer(ip)                :: lpoin
  real(rp)                   :: amatrd(ndof1)
  integer(ip), pointer       :: kfl_fixno(:,:)
  real(rp),    pointer       :: bvess(:,:)
  logical(lg), pointer       :: in_my_zone(:)
  real(rp),    pointer       :: reaction(:,:)

  nullify( kfl_fixno )
  nullify( bvess )
  nullify( in_my_zone )
  nullify( reaction )

  if( IMASTER ) return

  !----------------------------------------------------------------------
  !
  ! Reaction on Dirichlet nodes: save matrix and RHS before imposing Dirichlet b.c.s
  !
  !----------------------------------------------------------------------

  if( itask == 1 .or. itask == 4 ) then

     if( associated(solve_sol(1) % lpoin_reaction) .and. solve_sol(1) % block_num == 1 ) then
        !
        ! Allocate
        !
        num_blocks = solve_sol(1) % num_blocks

        if( num_blocks == 1 ) then
           !
           ! Monolithic system
           !
           ndofn = solve_sol(1) % ndofn
           do ipoin = 1,npoin
              if( solve_sol(1) % lpoin_reaction(ipoin) ) then
                 solve_sol(1) % lpoin_block(ipoin) % block1_num(1) % rhs(1:ndofn) = rhsid(1:ndofn,ipoin)
                 jzdom = 0
                 do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                    jzdom = jzdom + 1
                    solve_sol(1) % lpoin_block(ipoin) % block2_num(1,1) % matrix(1:ndofn,1:ndofn,jzdom) = amatr(1:ndofn,1:ndofn,izdom)
                 end do
              end if
           end do

        else

           call runend('PRESOL: NOT CODED FOR THIS NUMBER OF BLOCKS')

        end if
     end if

  end if

  !----------------------------------------------------------------------
  !
  ! Impose Neumann b.c.
  !
  !----------------------------------------------------------------------

  if( itask == 1 ) then
     if( solve_sol(1) % kfl_bvnat == 1 ) then
        if( solve_sol(1) % kfl_iffix == 1 ) then
           kfl_fixno => solve_sol(1) % kfl_fixno
           do ipoin = 1,npoin
              do idofn = 1,ndof1
                 if( solve_sol(1) % kfl_fixno(idofn,ipoin) <= 0 ) then
                    rhsid(idofn,ipoin) = rhsid(idofn,ipoin) + solve_sol(1) % bvnat(idofn,ipoin)
                 end if
              end do
           end do
        else
           do ipoin = 1,npoin
              rhsid(1:ndof1,ipoin) = rhsid(1:ndof1,ipoin) + solve_sol(1) % bvnat(1:ndof1,ipoin)
           end do
        end if
     end if
  end if

  !----------------------------------------------------------------------
  !
  ! Impose Dirichlet b.c.
  !
  !----------------------------------------------------------------------

  if( itask == 4 ) then
     !
     ! Richardson residual-based solvers: put RHS=0 on Dirichlet nodes
     !
     kfl_fixno => solve_sol(1) % kfl_fixno

     if( solve_sol(1) % kfl_iffix == 1 ) then
        do ipoin = 1,npoin
           do idofn = 1,ndof1
              if( kfl_fixno(idofn,ipoin) > 0 ) then
                 rhsid(idofn,ipoin) = 0.0_rp
              end if
           end do
        end do
     end if

  else if( itask == 1 ) then

     kfl_fixno => solve_sol(1) % kfl_fixno
     bvess     => solve_sol(1) % bvess

     if( solve_sol(1) % kfl_iffix == 1 .and. associated(solve_sol(1) % bvess) ) then
        !
        ! Dirichlet value is not given in BVESS
        !
        do ipoin = 1,npoin

           do idofn = 1,ndof1

              if( kfl_fixno(idofn,ipoin) > 0 ) then
                 !
                 ! Eliminate dof of IPOIN from other equations (JPOIN)
                 ! Keep rows unchanged in order to compute the reaction force
                 !
                 unkno(idofn,ipoin) = bvess(idofn,ipoin)
                 do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                    jpoin = c_dom(izdom)
                    if( ipoin /= jpoin ) then
                       do jzdom = r_dom(jpoin),r_dom(jpoin+1) - 1
                          kpoin = c_dom(jzdom)
                          if( kpoin == ipoin ) then
                             do jdofn = 1,ndof1
                                rhsid(jdofn,jpoin)       = rhsid(jdofn,jpoin) - amatr(idofn,jdofn,jzdom) * bvess(idofn,ipoin)
                                amatr(idofn,jdofn,jzdom) = 0.0_rp
                             end do
                          end if
                       end do
                    end if
                 end do
                 !
                 ! => Uncomment to cancel out rows with Dirichlet b.c.
                 !
                 !
                 ! IZDOD: Diagonal
                 !
                 izdod = r_dom(ipoin) - 1
                 jpoin = 0
                 do while( jpoin /= ipoin )
                    izdod = izdod + 1
                    jpoin = c_dom(izdod)
                 end do
                 amatrd(idofn) = amatr(idofn,idofn,izdod)
                 if( abs(amatrd(idofn)) < zeror ) amatrd(idofn) = 1.0_rp
                 !
                 ! Set line to zero
                 !
                 do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                    do jdofn = 1,ndof1
                       amatr(jdofn,idofn,izdom) = 0.0_rp
                    end do
                 end do
                 !
                 ! Presrcibe value
                 !
                 amatr(idofn,idofn,izdod) = amatrd(idofn)
                 rhsid(idofn,ipoin)       = bvess(idofn,ipoin) * amatrd(idofn)
                 !
                 ! <= Uncomment
                 !
              end if

           end do

        end do

     else if( solve_sol(1) % kfl_iffix == 2 .or. ( solve_sol(1) % kfl_iffix == 1 .and. .not. associated(solve_sol(1) % bvess) ) ) then
        !
        ! Dirichlet value is not given: Impose zero
        !
        do ipoin = 1,npoin

           do idofn = 1,ndof1

              if( kfl_fixno(idofn,ipoin) > 0 ) then
                 !
                 ! Eliminate dof of IPOIN from other equations (JPOIN)
                 !
                 do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                    jpoin = c_dom(izdom)
                    if( ipoin /= jpoin ) then
                       do jzdom = r_dom(jpoin),r_dom(jpoin+1) - 1
                          kpoin = c_dom(jzdom)
                          if( kpoin == ipoin ) then
                             do jdofn = 1,ndof1
                                amatr(idofn,jdofn,jzdom) = 0.0_rp
                             end do
                          end if
                       end do

                    end if
                 end do
                 !
                 ! IZDOD: Diagonal
                 !
                 izdod = r_dom(ipoin) - 1
                 jpoin = 0
                 do while( jpoin /= ipoin )
                    izdod = izdod + 1
                    jpoin = c_dom(izdod)
                 end do
                 amatrd(idofn) = amatr(idofn,idofn,izdod)
                 if( abs(amatrd(idofn)) < zeror ) amatrd(idofn) = 1.0_rp
                 !
                 ! Set line to zero
                 !
                 do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                    do jdofn = 1,ndof1
                       amatr(jdofn,idofn,izdom) = 0.0_rp
                    end do
                 end do
                 !
                 ! Prescribe value
                 !
                 amatr(idofn,idofn,izdod) = amatrd(idofn)
                 if( solve_sol(1) % kfl_iffix /= 2 ) &
                      rhsid(idofn,ipoin)       = 0.0_rp

              end if

           end do

        end do

     end if

  end if

  !----------------------------------------------------------------------
  !
  ! Recover Dirichlet condition
  !
  !----------------------------------------------------------------------

!!$  if( itask == 5 ) then
!!$
!!$     if( solve_sol(1) % kfl_iffix == 1 ) then
!!$        if( associated(solve_sol(1) % bvess) ) then
!!$           do lpoin = 1,npoiz(current_zone)
!!$              ipoin = lpoiz(current_zone) % l(lpoin)
!!$              do idofn = 1,ndof1
!!$                 if( solve_sol(1) % kfl_fixno(idofn,ipoin) > 0 ) unkno(idofn,ipoin) = solve_sol(1) % bvess(idofn,ipoin)
!!$              end do
!!$           end do
!!$        else
!!$           do lpoin = 1,npoiz(current_zone)
!!$              ipoin = lpoiz(current_zone) % l(lpoin)
!!$              do idofn = 1,ndof1
!!$                 if( solve_sol(1) % kfl_fixno(idofn,ipoin) > 0 ) unkno(idofn,ipoin) = 0.0_rp
!!$              end do
!!$           end do
!!$        end if
!!$     end if
!!$
!!$  end if

  !----------------------------------------------------------------------
  !
  ! Reaction on Dirichlet nodes: save matrix and RHS before imposing Dirichlet b.c.s
  !
  !----------------------------------------------------------------------

  if( itask == 5 ) then

     if( solve_sol(1) % kfl_react == 1 .and. associated(solve_sol(1) % lpoin_reaction) ) then

        num_blocks = solve_sol(1) % num_blocks

        if( num_blocks == 1 ) then
           !
           ! Monolithic system
           !
           ndofn = solve_sol(1) % ndofn
           do ipoin = 1,npoin
              solve_sol(1) % reaction(1:ndofn,ipoin) = 0.0_rp
              if( solve_sol(1) % lpoin_reaction(ipoin) ) then
                 jzdom = 0
                 solve_sol(1) % reaction(1:ndofn,ipoin) = solve_sol(1) % lpoin_block(ipoin) % block1_num(1) % rhs(1:ndofn)
                 do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                    jpoin = c_dom(izdom)
                    jzdom = jzdom + 1
                    do idofn = 1,ndofn
                       do jdofn = 1,ndofn
                          solve_sol(1) % reaction(idofn,ipoin) = solve_sol(1) % reaction(idofn,ipoin) &
                               & - solve_sol(1) % lpoin_block(ipoin) % block2_num(1,1) % matrix(jdofn,idofn,jzdom) &
                               & * unkno(jdofn,jpoin)
                       end do
                    end do
                 end do
              end if
           end do
           call PAR_INTERFACE_NODE_EXCHANGE(solve_sol(1) % reaction,'SUM','IN MY ZONE')
        else
           call runend('PRESOL: NOT CODED')
        end if

     end if

  end if

end subroutine presol
