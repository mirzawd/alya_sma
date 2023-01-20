!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Alefor
!> @{
!> @file    mod_ale_arrays.f90
!> @author  houzeaux
!> @date    2019-11-16
!> @brief   Alefor arrays
!> @details Alefor arrays
!-----------------------------------------------------------------------

module mod_ale_arrays

  use def_master
  use def_domain 
  use def_alefor
  use def_kermod
  use def_solver
  use mod_output_postprocess,  only : output_postprocess_check_variable_postprocess
  use mod_arrays,              only : arrays
  use mod_arrays,              only : arrays_number
  use mod_memory,              only : memory_alloca
  use mod_ADR,                 only : ADR_initialize_type
  use mod_ADR,                 only : ADR_check_and_compute_data
  use mod_ADR,                 only : ADR_allocate_projections_bubble_sgs
  use mod_ADR,                 only : ADR_arrays

  implicit none

  private

  public :: ale_arrays
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-16
  !> @brief   Alefor arrays
  !> @details Do what you have to do with alefor arrays
  !> 
  !-----------------------------------------------------------------------
  
  subroutine ale_arrays(wtask)

    character(len=*), intent(in) :: wtask

    call arrays(arrays_number('DISPM'),wtask,dispm    ,ndime,npoin,3_ip)
    call arrays(arrays_number('VELOM'),wtask,velom    ,ndime,npoin)
    call arrays(arrays_number('COALE'),wtask,coord_ale,ndime,npoin,3_ip)
    call arrays(arrays_number('COORI'),wtask,coord_ori,ndime,npoin)

  end subroutine ale_arrays
   
end module mod_ale_arrays
!> @}
