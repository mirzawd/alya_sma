!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup NeutroInput
!> @{
!> @file    neu_reaous.f90
!> @author  Guillaume Houzeaux
!> @date    29/03/2016
!> @brief   Read postprocess data
!> @details Read postprocess data
!> @} 
!-----------------------------------------------------------------------

subroutine neu_reaous()
  use def_parame
  use def_inpout
  use def_master
  use def_neutro
  use def_domain
  use mod_ecoute, only :  ecoute
  use mod_output_postprocess, only : output_postprocess_read
  implicit none

  if( INOTSLAVE ) then
     !
     ! Initializations
     !
     !
     ! Reach the section
     !
     call ecoute('neu_reaous')
     do while( words(1) /= 'OUTPU' )
        call ecoute('neu_reaous')
     end do
     call ecoute('neu_reaous')
     !
     ! Begin to read data
     !
     do while( words(1) /= 'ENDOU' )
        call output_postprocess_read()
        call ecoute('neu_reaous')
     end do

  end if

end subroutine neu_reaous

