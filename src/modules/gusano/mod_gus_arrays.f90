!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    mod_gus_arrays.f90
!> @author  houzeaux
!> @date    2019-11-16
!> @brief   Gusano arrays
!> @details Gusano arrays
!-----------------------------------------------------------------------

module mod_gus_arrays

  use def_master
  use def_domain 
  use def_gusano
  use def_kermod
  use mod_output_postprocess,  only : output_postprocess_check_variable_postprocess
  use mod_arrays,              only : arrays
  use mod_arrays,              only : arrays_number
  use mod_memory,              only : memory_alloca
  use mod_gusano,              only : gusano_memory_allocate
  
  implicit none

  private

  public :: gus_arrays
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-16
  !> @brief   Gusano arrays
  !> @details Do what you have to do with gusano arrays
  !> 
  !-----------------------------------------------------------------------
  
  subroutine gus_arrays(wtask)
    
    character(len=*), intent(in) :: wtask
    !
    ! Unknowns
    !
    call arrays(arrays_number('PRESS'),trim(wtask),press,npoin,3_ip)
    call arrays(arrays_number('FLOWR'),trim(wtask),flowr,npoin,3_ip)
    
  end subroutine gus_arrays
   
end module mod_gus_arrays
!> @}
