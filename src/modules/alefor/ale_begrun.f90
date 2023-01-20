!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Alefor
!> @{
!> @file    alebegrun.f90
!> @date    14/06/2019
!> @author  Guillaume Houzeaux
!> @brief   Beginning the run... 
!> @details Beginning the run... we can compute matrices!
!> @}
!-----------------------------------------------------------------------

subroutine ale_begrun()

  use def_master
  use def_alefor
  implicit none
  
  !----------------------------------------------------------------------
  !
  ! Define body fitted 
  !
  !---------------------------------------------------------------------- 

  if( kfl_rigid_ale == 1 ) call ale_bodfit()

end subroutine ale_begrun
 
