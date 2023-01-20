!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine reaous()
  !------------------------------------------------------------------------
  !****f* kernel/reamod
  ! NAME 
  !    reamod
  ! DESCRIPTION
  !    This routine reads postprocess
  ! USES
  ! USED BY
  !    ***_reapro
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use mod_ecoute, only :  ecoute
  use mod_output_postprocess, only : output_postprocess_read
  implicit none

  if( INOTSLAVE ) then
     !
     ! Reach the section
     !
     call ecoute('reaous')
     do while(words(1)/='OUTPU')
        call ecoute('reaous')
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDOU')
        call ecoute('reaous')
        call output_postprocess_read()
     end do
  end if

end subroutine reaous
