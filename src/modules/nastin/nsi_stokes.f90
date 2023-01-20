!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_stokes()
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_stokes
  ! NAME 
  !    nsi_stokes
  ! DESCRIPTION
  !    This routine solves an iteration for the incompressible
  !    Navier-Stokes equations, using:
  !    - A Monolithic scheme
  !    - A block Gauss-Seidel scheme
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
  implicit none
  integer(ip) :: kfl_advec_old,kfl_timei_old,kfl_cotem_old
  integer(ip) :: kfl_tiacc_old,kfl_sgsco_old
  integer(ip) :: kfl_sgsti_old
  real(rp)    :: dtinv_tmp
  !
  ! Initial solution is Stokes: save original values
  !
  kfl_advec_old = kfl_advec_nsi
  kfl_sgsco_old = kfl_sgsco_nsi
  kfl_timei_old = kfl_timei_nsi
  kfl_sgsti_old = kfl_sgsti_nsi
  kfl_cotem_old = kfl_cotem_nsi
  kfl_tiacc_old = kfl_tiacc_nsi
  dtinv_tmp     = dtinv_nsi
  kfl_advec_nsi = 0
  kfl_sgsco_nsi = 0
  kfl_sgsti_nsi = 0
  kfl_timei_nsi = 0
  kfl_cotem_nsi = 0
  dtinv_nsi     = 0.0_rp
  !
  ! Set up the solver parameters for the NS equations
  !
  call nsi_inisol(one)
  !
  ! Initial solution: UNKNO <= VELOC(:,:,1)
  !
  call nsi_updunk(2200_ip) 
  !
  ! Solve Stokes flow
  !
  call nsi_solite()
  !
  ! If initial solution is Stokes: recover original values
  !
  kfl_advec_nsi = kfl_advec_old 
  kfl_sgsco_nsi = kfl_sgsco_old 
  kfl_timei_nsi = kfl_timei_old 
  kfl_sgsti_nsi = kfl_sgsti_old 
  kfl_cotem_nsi = kfl_cotem_old 
  kfl_tiacc_nsi = kfl_tiacc_old
  dtinv_nsi     = dtinv_tmp
  !
  ! Update velocity and pressure: VELOC(:,:,1) <= UNKNO
  !
  call nsi_updunk(2300_ip)

end subroutine nsi_stokes
