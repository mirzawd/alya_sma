!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_iniqua()
  !-----------------------------------------------------------------------
  !****f* magnet/mag_quadra
  ! NAME
  !    mag_turnon
  ! DESCRIPTION
  !    This routine initializes the quadrature data
  ! USES
  ! USED BY
  !    mag_inivar
  !***
  !-----------------------------------------------------------------------

  use def_magnet
  use def_domain, only: ndime

  implicit none

  if (ndime == 1_ip) then
     quadLin = quadraLine(gslin_mag)
     maxgau_mag = size(quadLin % wq)
  elseif (ndime == 2_ip) then
     quadTri = quadraTriangle(gstri_mag)
     quadQua = quadraQuadrangle(gsqua_mag)
     maxgau_mag = max(size(quadTri % wq), size(quadQua % wq))
  elseif (ndime == 3_ip) then
     quadTet = quadraTetrahedron(gstet_mag)
     quadHex = quadraHexahedron(gshex_mag)
     maxgau_mag = max(size(quadTet % wq), size(quadHex % wq))
  else
     stop "mag_iniqua: dimension must be less or equal than 3"
  end if

end subroutine mag_iniqua
