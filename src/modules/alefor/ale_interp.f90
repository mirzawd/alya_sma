!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    ale_interp.f90
!> @author  houzeaux
!> @date    2019-06-17
!> @brief   Redistribution of arrays
!> @details Redistribution of arrays. mainly allocated ale_memall
!> @} 
!-----------------------------------------------------------------------

subroutine ale_interp()

  use mod_output_postprocess, only : output_postprocess_check_variable_postprocess
  use mod_ale_arrays,         only : ale_arrays
  implicit none
  !
  ! Recompute boundary conditions
  !
  call ale_inibcs()
  !
  ! Variables in memall
  !
  call ale_arrays('INTERPOLATE')
  !
  ! Update boundary conditions
  !
  call ale_updbcs()

end subroutine ale_interp
