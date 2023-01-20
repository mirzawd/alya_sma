!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_turnon()
  !-----------------------------------------------------------------------
  !****f* magnet/mag_turnon
  ! NAME
  !    mag_turnon
  ! DESCRIPTION
  !    This routine starts the run of the 'magnet' module
  ! USES
  !    mag_reaphy
  !    mag_parall
  !    mag_inivar
  !    mag_reanut
  !    mag_reaous
  !    mag_reabcs
  !    mag_inibcs
  !    mag_memall
  ! USED BY
  !    Magnet
  !***
  !-----------------------------------------------------------------------

  use def_kintyp, only : ip

  implicit none
  !
  ! Initial variables
  !
  call mag_inivar(0_ip)
  !
  ! Read the physical problem data
  !
  call mag_reaphy()
  !
  ! Read the numerical parameters
  !
  call mag_reanut()
  !
  ! Read the output strategy
  !
  call mag_reaous()
  !
  ! Service: Parall
  !
  call mag_parall(1_ip)
  !
  ! Read the boundary conditions
  !
  call mag_reabcs()
  !
  ! Service: Parall
  !
  call mag_parall(2_ip)
  !
  ! Axisymmetric option
  !
  call mag_inivar(1_ip)
  !
  ! Modify the boundary conditions
  !
  call mag_inibcs(1_ip)
  !
  ! Memory allocation: solver and variables
  !
  call mag_memall()
  !
  ! Modify the boundary conditions
  !
  call mag_inibcs(2_ip)
  !
  ! Compute additional geometry data
  !
  call mag_inivar(2_ip)
  !
  ! Initialize time variables
  !
  call mag_inivar(3_ip)
  !
end subroutine mag_turnon
