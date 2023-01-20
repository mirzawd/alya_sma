!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_dirich()
  !-----------------------------------------------------------------------
  !****f* Magnet/mag_dirich
  ! NAME
  !    mag_dirich
  ! DESCRIPTION
  !    This routine sets Dirichlet boundary conditions
  ! USES
  !    dirichletField
  ! USED BY
  !    mag_begste
  !***
  !-----------------------------------------------------------------------
  use def_magnet, only: ip, rp, diredg_mag, locboucen_mag, locboutan_mag, myclocktime_mag, dt_mag, Hxtl_mag
  use def_domain, only: ndime
  use mod_memory, only: memory_size
  use mod_mag_inpdat, only: mag_dirfie

  implicit none

  integer(ip) :: &
    iedge

  do iedge = 1, memory_size(diredg_mag)
    !
    ! Scalar projection of Dirichlet vector field onto edge unit tangent vector: He_mag = H_dir'*t
    !
    Hxtl_mag(iedge) = dot_product(mag_dirfie(myclocktime_mag + dt_mag, locboucen_mag(1:ndime, iedge)), locboutan_mag(1:ndime, iedge))
    !
  end do

end subroutine mag_dirich
