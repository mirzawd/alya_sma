!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Turnon
!> @{
!> @file    Reamod.f90
!> @author  houzeaux
!> @date    2019-06-12
!> @brief   Read module data
!> @details Read module data
!> @} 
!-----------------------------------------------------------------------

subroutine Reamod()
  
  use def_kintyp,   only : ip
  use def_master,   only : ITASK_TURNON
  use mod_messages, only : messages_live
  use mod_moduls,   only : moduls

  implicit none
  !
  ! Read modules
  !
  call messages_live('MODULE DATA','START SECTION')
  call moduls(ITASK_TURNON)
  call messages_live('MODULE DATA','END SECTION')
  !
  ! Close module data files and open occasional output files
  !  
  call moddef(2_ip)

end subroutine Reamod
