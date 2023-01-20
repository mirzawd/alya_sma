!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




subroutine mag_upddir()
  !-----------------------------------------------------------------------
  !****f* Magnet/mag_inidir
  ! NAME
  !    mag_upddir
  ! DESCRIPTION
  !    This routine sets Dirichlet boundary conditions to zero
  ! USES
  !    dirichletField
  ! USED BY
  !    mag_begste
  !***
  !-----------------------------------------------------------------------
  use def_magnet, only: ip, He_mag, diredg_mag, Hxtl_mag, Hsfl_mag, kfl_self_mag
  use mod_memory, only: memory_size

  implicit none

  integer(ip) :: &
    iedge

  if (kfl_self_mag) then
    !
    do iedge = 1, memory_size(diredg_mag)
      !
      He_mag(diredg_mag(iedge)) = Hxtl_mag(iedge) + Hsfl_mag(iedge)
      !
    end do
    !
  else
    !
    do iedge = 1, memory_size(diredg_mag)
      !
      He_mag(diredg_mag(iedge)) = Hxtl_mag(iedge)
      !
    end do
    !
  end if
  !
end subroutine mag_upddir
