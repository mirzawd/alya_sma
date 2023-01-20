!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Turbul
!> @{
!> @file    tur_restar.f90
!> @author  houzeaux
!> @date    2022-09-22
!> @brief   Restart
!> @details Turbul restart variables
!> @} 
!-----------------------------------------------------------------------

subroutine tur_restar(itask)

  use def_kintyp_basic, only : ip
  use def_master,       only : ITASK_READ_RESTART
  use def_master,       only : ITASK_WRITE_RESTART
  use mod_tur_arrays,   only : tur_arrays
  implicit none
  
  integer(ip), intent(in) :: itask
  
  !----------------------------------------------------------------------
  !
  ! Primary arrays
  !
  !----------------------------------------------------------------------
  
  if( itask == ITASK_READ_RESTART ) then
     call tur_arrays('READ RESTART')
  else if( itask == ITASK_WRITE_RESTART ) then
     call tur_arrays('WRITE RESTART')
  end if
  
end subroutine tur_restar
 
 
