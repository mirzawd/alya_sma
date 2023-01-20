!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Interp
!> @{
!> @file    Interp.f90
!> @author  houzeaux
!> @date    2019-06-17
!> @brief   Interpribution
!> @details interpribution of data
!> @} 
!-----------------------------------------------------------------------

subroutine Interp()
  
  use def_kintyp,   only : ip
  use def_master,   only : ITASK_INTERP
  use def_master,   only : iblok
  use def_master,   only : nblok
  use mod_moduls,   only : moduls
  use mod_messages, only : messages_live
  implicit none
  
  call messages_live('INTERPOLATION','START SECTION')
  do iblok = 1_ip, nblok
     call moduls(ITASK_INTERP)
  end do
  call Kermod(ITASK_INTERP)
  call messages_live('INTERPOLATION','END SECTION')
    
end subroutine Interp
