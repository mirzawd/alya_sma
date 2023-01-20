!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    pts_begrun.f90
!> @date    14/06/2019
!> @author  Guillaume Houzeaux
!> @brief   Beginning the run... 
!> @details Beginning the run... we can compute matrices!
!> @}
!-----------------------------------------------------------------------

subroutine pts_begrun()

  use def_master,                only : velom
  use def_partis,                only : if_moving_mesh_pts
  use mod_communications_global, only : PAR_OR
  implicit none
  !
  ! Wall distance 
  !
  call pts_waldis()
  !
  ! Check if mesh is moving
  !
  if_moving_mesh_pts = associated(velom)
  call PAR_OR(if_moving_mesh_pts)

end subroutine pts_begrun
 
