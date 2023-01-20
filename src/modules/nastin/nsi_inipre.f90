!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



 subroutine nsi_inipre()
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_inipre
  ! NAME
  !    nsi_inipre
  ! DESCRIPTION
  !    This routine solves the initial pressure
  ! USES
  !    nsi_ifconf
  !    nsi_solmon
  !    nsi_solbgs
  !    nsi_rotunk
  ! USED BY
  !    nsi_doiter
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  use def_solver, only : SOL_MATRIX_HAS_CHANGED
  use mod_messages, only : livinf
  use mod_communications
  use mod_matrices,  only : matrices_laplacian
  use mod_memory,    only : memory_alloca
  use mod_memory,    only : memory_deallo
  implicit none
  integer(ip)         :: izdom,ipoin,iboun,inodb,ibopo
  integer(ip)         :: jbopo,izdod,jpoin,jzdom
!  real(rp)            :: dtinv_nsi_save
  real(rp)            :: p0
  real(rp)            :: Qd
  real(rp),   pointer :: lapl(:)
  
  call livinf(59_ip,'INITIAL PRESSURE USING LAPLACIAN',0_ip)

  nullify(lapl)
  ivari_nsi =  ivari_nsi_cont
  solve_sol => solve(2:)
  solve_sol %  kfl_assem = SOL_MATRIX_HAS_CHANGED

  if( INOTMASTER ) then
     !
     ! Initialize matrix
     !
     do ipoin = 1,npoin
        rhsid(ipoin) = 0.0_rp
        unkno(ipoin) = 0.0_rp
     end do
     !
     ! Boundary values
     !
     call memgen(zero,npoin,zero)

     do iboun = 1,nboun
        if( kfl_fixbo_nsi(iboun) == 2 .or. kfl_fixbo_nsi(iboun) == 6 .or. kfl_fixbo_nsi(iboun) == 20 ) then
           do inodb = 1,nnode(ltypb(iboun))
              ipoin = lnodb(inodb,iboun)
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) then
                 if( kfl_fixpr_nsi(1,ipoin) > 0 ) then
                    p0 = bvnat_nsi(1,iboun,1)
                    ! Set pressure to non-zero value if initial value is zero
                    if ( abs(p0) < 1e-10_rp ) p0 = 1e-10_rp
                    gesca(ipoin) = p0
                 end if
              end if
           end do
        end if
     end do
     call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MAX')

     call matrices_laplacian(Lapl)

!!$     !
!!$     ! Define some variables
!!$     !
!!$     dtinv_nsi_save = dtinv_nsi
!!$     !
!!$     ! Assemble equation
!!$     !
!!$     if( kfl_vector_nsi == 0 ) then
!!$        !
!!$        ! Classical version
!!$        !
!!$        call nsi_elmope_omp(6_ip)
!!$     else
!!$        !
!!$        ! Vectorized version
!!$        !
!!$        call nsi_elmope_all(6_ip)
!!$     end if
!!$     !
!!$     ! Recover original variables
!!$     !
!!$     dtinv_nsi = dtinv_nsi_save 
     !
     ! Impose b.c.
     !
     if( solve(2)%kfl_symme == 1 ) then

        call runend('NIS_INIPRE: NOT CODED: CHECK IT')
        do ipoin = 1,npoin
           do izdom = r_sym(ipoin),r_sym(ipoin+1) - 2
              jpoin = c_sym(izdom)
              jbopo = lpoty(jpoin)
              if( jbopo /= 0 ) then
                 if( kfl_fixpr_nsi(1,jpoin) > 0 ) then
                    rhsid(ipoin) = rhsid(ipoin) - lapl(izdom) * gesca(jpoin)
                    lapl(izdom) = 0.0_rp
                 end if
              end if
           end do
        end do

        do ipoin = 1,npoin

           ibopo = lpoty(ipoin)
           if( ibopo > 0 ) then

              if( kfl_fixpr_nsi(1,ipoin) > 0 ) then
                 !
                 ! IZDOD: Diagonal
                 !
                 izdod = r_sym(ipoin+1) - 1
                 Qd = lapl(izdod)
                 if( abs(Qd) < zeror ) Qd = 1.0_rp
                 !
                 ! Set line to zero
                 !
                 do izdom = r_sym(ipoin),r_sym(ipoin+1) - 1
                    lapl(izdom) = 0.0_rp
                 end do
                 !
                 ! Prescribe value
                 !
                 lapl(izdod)  = Qd
                 rhsid(ipoin) = Qd * gesca(ipoin)
                 unkno(ipoin) = gesca(ipoin)

              end if

           end if

        end do

     else

        do jpoin = 1,npoin
           do jzdom = r_dom(jpoin),r_dom(jpoin+1) - 1
              ipoin = c_dom(jzdom)
              ibopo = lpoty(ipoin)
              if ( ibopo /= 0 .and. ipoin /= jpoin ) then
                 if( kfl_fixpr_nsi(1,ipoin) > 0 ) then
                    rhsid(jpoin) = rhsid(jpoin) - lapl(jzdom) * gesca(ipoin)
                    lapl(jzdom) = 0.0_rp
                 end if
              end if
           end do
        end do

        do ipoin = 1,npoin

           ibopo = lpoty(ipoin)
           if( ibopo > 0 ) then

              if( kfl_fixpr_nsi(1,ipoin) > 0 ) then
                 !
                 ! IZDOD: Diagonal
                 !
                 izdod = r_dom(ipoin) - 1
                 jpoin = 0
                 do while( jpoin /= ipoin )
                    izdod = izdod + 1
                    jpoin = c_dom(izdod)
                 end do
                 Qd = lapl(izdod)
                 if( abs(Qd) < zeror ) Qd = 1.0_rp
                 !
                 ! Set line to zero
                 !
                 do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                    lapl(izdom) = 0.0_rp
                 end do
                 !
                 ! Prescribe value
                 !
                 lapl(izdod)  = Qd
                 rhsid(ipoin) = Qd * gesca(ipoin)
                 unkno(ipoin) = gesca(ipoin)

              end if

           end if
        end do

     end if
     call memgen(two,npoin,zero)

  else

     call memory_alloca(memor_dom,'Lapl','matrices_gradient',Lapl,1_ip)
     
  end if
  !
  ! Solve system
  !
  call solver(rhsid,unkno,lapl,pmatr)
  !
  ! Deallocate
  !
  call memory_deallo(memor_dom,'Lapl','matrices_gradient',Lapl)
  !
  ! Tell solver that matrix will be changed
  !
  solve_sol %  kfl_assem = SOL_MATRIX_HAS_CHANGED
  !
  ! Update pressure
  !
  if( INOTMASTER ) then
     do ipoin = 1,npoin
        press(ipoin,1) = unkno(ipoin)
        press(ipoin,2) = unkno(ipoin)
        press(ipoin,3) = unkno(ipoin)
     end do
  end if

end subroutine nsi_inipre
