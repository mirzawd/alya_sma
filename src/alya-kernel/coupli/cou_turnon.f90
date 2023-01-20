!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Coupling 
!> @{
!> @file    cou_turnon.f90
!> @date    03/03/2014
!> @author  Guillaume Houzeaux
!> @brief   Turn on coupling 
!> @details Read data and allocate memory
!> @} 
!------------------------------------------------------------------------
subroutine cou_turnon()
  use def_kintyp, only : ip
  implicit none
  !
  ! Open files
  !  
  call cou_openfi()
  !
  ! Read data
  !  
  call cou_readat()
  !
  ! Broadcast data
  !  
  call cou_parall()

end subroutine cou_turnon
