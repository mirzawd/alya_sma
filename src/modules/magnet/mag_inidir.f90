!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_inidir()
  !-----------------------------------------------------------------------
  !****f* Magnet/mag_inidir
  ! NAME
  !    mag_inidir
  ! DESCRIPTION
  !    This routine sets Dirichlet boundary conditions to zero
  ! USES
  !    dirichletField
  ! USED BY
  !    mag_begste
  !***
  !-----------------------------------------------------------------------
  use def_magnet, only: rp, He_mag, diredg_mag
!  use mod_memory, only: memory_size

  implicit none

!  integer(ip) :: &
!    iedge

  if( associated(He_mag) .and. associated(diredg_mag) ) then
     if( diredg_mag(1) /= 0 ) He_mag(diredg_mag) = 0.0_rp
  end if

!  do iedge = 1, memory_size(diredg_mag)
!    !
!    He_mag(diredg_mag(iedge)) = 0.0_rp
!    !
!  end do
 
end subroutine mag_inidir
