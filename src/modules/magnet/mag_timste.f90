!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_timste()
  !-----------------------------------------------------------------------
  !****f* Magnet/mag_timste
  ! NAME
  !    mag_timste
  ! DESCRIPTION
  !    This routine computes the time step
  ! USES
  ! USED BY
  !    Magnet
  !***
  !-----------------------------------------------------------------------

  use def_kintyp, only: rp
  use def_master, only: dtime, dtinv
  use def_magnet, only: dt_mag

  implicit none
  
  dtime = dt_mag
  dtinv = 1.0_rp / dt_mag

end subroutine mag_timste
