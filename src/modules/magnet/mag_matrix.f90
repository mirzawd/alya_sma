!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_matrix()
  !-----------------------------------------------------------------------
  !****f* magnet/mag_matrix
  ! NAME
  !    mag_matrix
  ! DESCRIPTION
  !    This routine constructs the system matrix, RHS, and constraints if any
  ! USES
  !    mag_elmope
  !    mag_asscon
  ! USED BY
  !    mag_solite
  !-----------------------------------------------------------------------

  use def_kintyp, only: rp
  use def_master, only: INOTMASTER
  use def_magnet, only: kfl_lagr_mag, magnen_mag, joulen_mag, magtiz_mag, cursum_mag, magsum_mag, volume_mag

  implicit none

  magnen_mag = 0.0_rp
  joulen_mag = 0.0_rp

  magtiz_mag = 0.0_rp
  cursum_mag = 0.0_rp
  magsum_mag = 0.0_rp
  volume_mag = 0.0_rp

  if (INOTMASTER) then
    if (kfl_lagr_mag) call mag_asscon()
    call mag_elmope()
  end if

end subroutine mag_matrix
