!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    ker_doiter.f90
!> @author  houzeaux
!> @date    2020-09-16
!> @brief   Solve problems
!> @details Solve some physical problems to be used by all modules
!> @} 
!-----------------------------------------------------------------------

subroutine ker_doiter()

  use mod_tubes, only : tubes_doiter
  implicit none
  !
  ! Check if tubes should be solved for a module
  !
  call tubes_doiter()
  
end subroutine ker_doiter
