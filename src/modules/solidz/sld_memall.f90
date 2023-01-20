!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_memall.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   General memory allocation
!> @details General memory allocation
!> @}
!-----------------------------------------------------------------------

subroutine sld_memall()

  use mod_sld_arrays, only : sld_arrays
  use mod_sld_cardiac_cycle, only : cardiac_cycle_redist

  implicit none

  !----------------------------------------------------------------------
  !
  ! Primary variables
  !
  !----------------------------------------------------------------------
  call sld_arrays('ALLOCATE')

  !Entering here after reading .dat
  call cardiac_cycle_redist()

end subroutine sld_memall

