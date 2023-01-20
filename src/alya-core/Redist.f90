!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Redist
!> @{
!> @file    Redist.f90
!> @author  houzeaux
!> @date    2019-06-17
!> @brief   Redistribution
!> @details redistribution of data
!> @} 
!-----------------------------------------------------------------------

subroutine Redist()
  
  use def_kintyp,   only : ip
  use def_master,   only : ITASK_REDIST
  use def_master,   only : iblok
  use def_master,   only : nblok
  use mod_moduls,   only : moduls
  use mod_messages, only : messages_live
  implicit none
  
  call messages_live('REDISTRIBUTION','START SECTION')
  do iblok = 1_ip, nblok
     call moduls(ITASK_REDIST)
  end do
  call Kermod(ITASK_REDIST)
  call messages_live('REDISTRIBUTION','END SECTION')
    
end subroutine Redist
