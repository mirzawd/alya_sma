!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_chkaxs()

  use def_kintyp, only: ip
  use def_domain, only: ndime
  use def_magnet, only: kfl_axsym_mag

  implicit none

  if (ndime /= 2_ip) then
    !
    ! Turn off axisymmetric option if problem dimension is not 2-D
    !
    kfl_axsym_mag = .false.
    !
  end if

end subroutine mag_chkaxs
