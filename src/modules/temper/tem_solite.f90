!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_solite()
  !-----------------------------------------------------------------------
  !****f* Temper/tem_solite
  ! NAME 
  !    tem_solite
  ! DESCRIPTION
  !    This routine solves an iteration of the temperature equations.
  ! USES
  !    tem_matrix
  !    Soldir
  !    Solite
  ! USED BY
  !    tem_doiter
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use mod_gradie
  use mod_solver, only : solver_solve

  implicit none
  integer(ip) :: kfl_advec_old,kfl_timei_old
  real(rp)    :: dtinv_tmp2,dtinv_tmp
  !
  ! Update inner iteration counter
  !
  itinn(modul) = itinn(modul) + 1
  ittot_tem    = ittot_tem + 1
  !
  ! Update boundary conditions
  !
  call tem_updbcs(ITASK_INNITE)
  !
  ! Compute temperature gradients
  !
  if( kfl_ellen_tem == -1 .and. INOTEMPTY ) call gradie(tempe(:,1),grtem_tem)
  !
  ! If initial solution is Stokes: save original values
  !
  if( kfl_inidi_tem == 1 ) then
     kfl_advec_old = kfl_advec_tem
     kfl_timei_old = kfl_timei_tem
     dtinv_tmp2    = dtinv_tem
     kfl_advec_tem = 0
     kfl_timei_tem = 0
     dtinv_tem     = 0.0_rp
  end if
  !
  ! Construct the system matrix and right-hand-side
  !
  if( solve(1) % kfl_algso == -2 ) then
     dtinv_tmp     = dtinv_tem
     dtinv_tem     = 0.0_rp
     kfl_timei_tem = 0
  end if

  call tem_matrix()

  if( solve(1) % kfl_algso == -2 ) then
     dtinv_tem     = dtinv_tmp
     kfl_timei_tem = 1
  end if 
  !
  ! Solve the algebraic system
  !
  call solver_solve(momod(modul) % solve,amatr,rhsid,unkno,pmatr)
  !
  ! If initial solution is Stokes: recover original values
  !
  if( kfl_inidi_tem == 1 ) then
     kfl_inidi_tem = 0
     kfl_advec_tem = kfl_advec_old 
     kfl_timei_tem = kfl_timei_old 
     dtinv_tem     = dtinv_tmp2
  end if
  !
  ! Actualize subgrid scale
  !
  !if( INOTMASTER .and. kfl_sgsti_tem == 1 ) call tem_elmope(2_ip)
  
end subroutine tem_solite
