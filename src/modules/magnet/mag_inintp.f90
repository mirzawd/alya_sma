!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_inintp()
  !-----------------------------------------------------------------------
  !****f* magnet/mag_inintp
  ! NAME
  !    mag_turnon
  ! DESCRIPTION
  !    This routine initializes the interpolation data
  ! USES
  ! USED BY
  !    mag_inivar
  !***
  !-----------------------------------------------------------------------

  use def_magnet

  implicit none

  integer(ip) :: intp

  character(20), allocatable :: filename_mag(:)

  allocate(intp1_mag(nintp1_mag), filename_mag(nintp1_mag))

  do intp = 1_ip, nintp1_mag
    write(filename_mag(intp), '(A,I3.3,A)') 'interp1-', intp, '.inp'
    intp1_mag(intp) = interpReadat(filename_mag(intp), ioun5_mag)
  end do

end subroutine mag_inintp

