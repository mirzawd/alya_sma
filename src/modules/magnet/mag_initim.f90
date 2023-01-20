!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_initim()

  use def_master, only: timei, timef, kfl_gotim
  use def_magnet

  implicit none
  !
  ! Initial time step
  !
  dt_mag = dtmin_mag
  !
  ! Initialize my own time variable
  !
  myClockTime_mag = timei
  !
  ! Initialize time step
  !
  timeStep_mag = 1_ip
  !
  if (myClockTime_mag < timef) then
    kfl_gotim = 1_ip
  else
    kfl_gotim = 0_ip
  end if
  !
end subroutine mag_initim
