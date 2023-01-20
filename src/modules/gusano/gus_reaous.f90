!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup GusanoInput
!> @{
!> @file    gus_reaous.f90
!> @author  Guillaume Houzeaux
!> @brief   Read postprocess data
!> @details Read postprocess data
!!          - Array postprocess
!!          - Witness points
!!          - Element, boundary and node sets
!> @} 
!-----------------------------------------------------------------------
subroutine gus_reaous
!-----------------------------------------------------------------------
 !.md<module>gusano
 !.md<input>case.gus.dat
 !.md<pos>2
 !.md<sec>
  
  use def_parame
  use def_inpout
  use def_master
  use def_gusano
  use def_domain
  use mod_ecoute,             only : ecoute
  use mod_ecoute,             only : ecoute_reach_section
  use mod_output_postprocess, only : output_postprocess_read
  implicit none

  if( INOTSLAVE ) then
     !
     ! Initializations
     !
     call ecoute_reach_section('OUTPU')
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDOU')
        call ecoute('gus_reaous')
        call output_postprocess_read()
     end do
  
  end if

end subroutine gus_reaous
    
 
