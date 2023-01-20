!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ker_walnor(itask)
  !-----------------------------------------------------------------------
  !****f* domain/ker_walnor
  ! NAME
  !    ker_walnor
  ! DESCRIPTION
  !    Compute the generalized distance to the wall via a
  !    Poisson equation:
  !    1. Solve Lapl(f)=-1, with f=0 on wall
  !    2. d=sqrt[ grad(f)^2 +2*f ] - sqrt[grad(f)^2]
  !    See the following references:
  !    P.G. Tucker, Differential equation-based wall distance computation for
  !         DES and RANS, J. Comp. Phys. 190 (2003) 229-248.
  !    P.G. Tucker, Int. J. Numer. Fluids 33 (2000) 869.
  !    P.G. Tucker, Appl. Math. Model. 22 (1998) 293.
  ! USES
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_elmtyp
  use def_master
  use def_kermod
  use def_domain
  use def_solver
  use mod_ADR,            only : ADR_assemble_laplacian
  use mod_solver,         only : solver_solve
  use mod_messages,       only : messages_live
  use mod_moduls_conf,    only : moduls_set_current_solver
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: idime,ipoin,ibopo
  real(rp)                :: xnorm

  if( kfl_walln /= 0 ) then

     call moduls_set_current_solver(ID_KERMOD,'WALL_NORMAL')

     if( solve_sol(1) % kfl_algso /= -999 ) then
        !
        ! Initialize solver
        !
        call messages_live('KERMOD: EXTEND NORMAL FROM THE WALL')        
        solve_sol(1) % kfl_iffix = 1

        do idime = 1,ndime

           if( INOTMASTER ) then

              call inisol()

              solve_sol(1) % bvess     => bvess_walln_ker
              solve_sol(1) % kfl_fixno => kfl_fixno_walln_ker

              do ipoin = 1,npoin
                 ibopo = lpoty(ipoin)
                 if( ibopo > 0 .and. kfl_fixno_walln_ker(1,ipoin) == 1 ) then
                    bvess_walln_ker(1,ipoin) = exnor(idime,1,ibopo)
                    unkno(ipoin) = bvess_walln_ker(1,ipoin)
                 else
                    bvess_walln_ker(1,ipoin) = 0.0_rp
                    unkno(ipoin) = 0.0_rp
                 end if
              end do

              call ADR_assemble_laplacian(meshe(ndivi),elmar,amatr)
           end if
           !
           ! Solve system
           !
           call solver_solve(solve_sol,amatr,rhsid,unkno,pmatr)
           !call solver(rhsid,unkno,amatr,pmatr)

           if( INOTMASTER ) then
              do ipoin = 1,npoin
                 walln(idime,ipoin) = unkno(ipoin)
              end do
           end if

        end do
        !
        ! Normalize normal
        !
        if( INOTMASTER ) then
           do ipoin = 1,npoin
              xnorm = sqrt( dot_product(walln(1:ndime,ipoin),walln(1:ndime,ipoin)) + zeror )
              walln(1:ndime,ipoin) = walln(1:ndime,ipoin) / xnorm
           end do
        end if
     end if
  end if

end subroutine ker_walnor
