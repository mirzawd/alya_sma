!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_mag_elmgeo

  use def_kintyp, only: ip, rp
  use mod_mag_linalg, only: matdet

  implicit none

contains

  !##############################################################
  function mag_trivol(n0, n1, n2) result(vol)
    !
    ! Area of a triangle
    !
    implicit none

    real(rp), intent(in), dimension(2) :: n0, n1, n2

    real(rp), dimension(2,2) :: matrix

    real(rp) :: vol

    matrix(:, 1) = n1 - n0
    matrix(:, 2) = n2 - n0

    vol = abs(matdet(matrix)) / 2.0_rp

  end function mag_trivol
  !##############################################################


  !##############################################################
  function mag_quavol(n0, n1, n2, n3) result(vol)
    !
    ! Area of a quadrangle
    !
    implicit none

    real(rp), intent(in), dimension(2) :: n0, n1, n2, n3

    real(rp) :: vol

    vol = mag_trivol(n0, n1, n2) + mag_trivol(n0, n2, n3)

  end function mag_quavol
  !##############################################################


  !##############################################################
  function mag_hexvol(n1, n2, n3, n4, n5, n6, n7, n8) result(vol)
    !
    ! Volume of a heaxaedron from sum of volumes of 3 pyramides
    !
    implicit none

    real(rp), intent(in), dimension(3) :: n1, n2, n3, n4, n5, n6, n7, n8

    real(rp) :: vol

    vol = mag_pyrvol(n2, n1, n5, n4, n8) + mag_pyrvol(n2, n5, n6, n7, n8) + mag_pyrvol(n2, n4, n3, n7, n8)

  end function mag_hexvol
  !##############################################################


  !##############################################################
  function mag_pyrvol(n0, n1, n2, n3, n4) result(vol)
    !
    ! Volume of a pyramide from sum of volumes of 2 tetrahedra
    !
    implicit none

    real(rp), intent(in), dimension(3) :: n0, n1, n2, n3, n4

    real(rp) :: vol

    vol = mag_tetvol(n0, n1, n2, n3) + mag_tetvol(n0, n1, n2, n4)

  end function mag_pyrvol
  !##############################################################


  !##############################################################
  function mag_tetvol(n0, n1, n2, n3) result(vol)
    !
    ! Volume of a tetrahedron
    !
    implicit none

    real(rp), intent(in), dimension(3) :: n0, n1, n2, n3

    real(rp), dimension(3,3) :: matrix

    real(rp) :: vol

    matrix(:, 1) = n1 - n0
    matrix(:, 2) = n2 - n0
    matrix(:, 3) = n3 - n0

    vol = abs(matdet(matrix)) / 6.0_rp

  end function mag_tetvol
  !##############################################################

end module mod_mag_elmgeo
