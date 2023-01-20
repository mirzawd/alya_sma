!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_endste()
  !-----------------------------------------------------------------------
  !****f* Magnet/mag_endste
  ! NAME 
  !    mag_endste
  ! DESCRIPTION
  !    This routine ends a time step of the H-formulation of Maxwell's
  !    equations
  ! USES
  ! USED BY
  !    Magnet
  !***
  !-----------------------------------------------------------------------

  use def_master, only: INOTSLAVE, kfl_gotim, timef !, ittim, cutim
  use def_magnet
!  use mod_communications, only: PAR_SUM

  implicit none
  
!  integer(ip) :: pmate
  !
  ! To get 'further' information check 'master/Endste.f90'
  ! 
  ! if (cutim >= timef - epsilon(1.0_rp)) kfl_gotim = 0
  ! If current time is greater than final time stop time loop
  !
  ! if (ittim >= mitim) kfl_gotim = 0
  ! If time iteration is greater than maximum number of steps stop time loop
  !
  ! Tell Alya to continue till it reached the final time or maximum number of time steps
  ! If set to 0, Alya will perform just one time step (nice for debuging)
  !
  kfl_gotim = 1_ip
  !
  !
  ! Time step converged:
  !
  ! 0. Compute energies and write to file
  !
  call mag_wridat(1_ip, myClockTime_mag + dt_mag * theta_mag)
  !
!!  do pmate = 1, maxmat_mag
!!    call PAR_SUM(magnen_mag(pmate))
!!    call PAR_SUM(joulen_mag(pmate))
!!  end do
!  !
!!  if (INOTSLAVE) then
!!    open(unit=ioun1_mag, file='magnet.nrj', position='append')
!!    write(ioun1_mag, 80) myClockTime_mag + dt_mag * theta_mag, magnen_mag, joulen_mag
!!    80 format (e15.8, e15.8, e15.8)
!!    close(ioun1_mag)
!!  end if
  !
  ! 1. Update my own time variable...
  !
  myClockTime_mag = myClockTime_mag + dt_mag
  !
  ! 2. Update time step
  !
  timeStep_mag = timeStep_mag + 1_ip
  !
  ! 3. Compute new time step
  !
  dt_mag = min( max( dtmin_mag, min( dtmax_mag, real(nlide_mag, rp) * dt_mag / real(nonlCount_mag, rp), 3.0_rp * dt_mag ) ), timef - myclocktime_mag)
  !
  ! 4. Show some info to the user...
  !
  if (INOTSLAVE) then
    write(*, 37) timeStep_mag, myClockTime_mag + dt_mag, dt_mag
    37 format ('    Next Time Step = ', i5, ' t = ', e11.4, ', dt = ', e11.4)
  end if
  !
end subroutine mag_endste
