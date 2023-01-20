!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_reaous()
  !--------------------------------------
  !****f* Magnet/mag_reaous
  ! NAME
  !    mag_reaous
  ! DESCRIPTION
  !    This routine reads data
  ! OUTPUT
  ! USES
  ! USED BY
  !    Reapro
  !***
  !--------------------------------------
  use def_inpout
  use def_master
  use def_magnet
  use mod_ecoute, only : ecoute
  use mod_output_postprocess, only : output_postprocess_read

  implicit none

  if (INOTSLAVE) then
    !
    ! Reach the section
    !
    call ecoute('mag_reaous')
    do while (words(1) /= 'OUTPU')
      call ecoute('mag_reaous')
    end do
    call ecoute('mag_reaous')
    !
    ! Begin to read data
    !
    do while (words(1) /= 'ENDOU')
      !
      ! Read Post-Process
      !
       call output_postprocess_read()
      !
!! WIP
      if (words(1) == 'NRJOU') then
        kfl_nrj_mag = option('NRJOU')
      else if (words(1) == 'MTZOU') then
        kfl_mtz_mag = option('MTZOU')
      else if (words(1) == 'CRNOU') then
        kfl_crn_mag = option('CRNOU')
      else if (words(1) == 'VLMOU') then 
        kfl_vlm_mag = option('VLMOU')
      end if
!! WIP
      call ecoute('mag_reaous')
      !
    end do
    !
  end if

end subroutine mag_reaous
