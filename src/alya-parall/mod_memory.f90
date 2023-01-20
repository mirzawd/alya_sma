!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!>
!> @defgroup Memory_Toolbox
!> @{
!> @name    ToolBox for memory management
!> @file    mod_memory.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for memory management
!> @details ToolBox for memory management
!>
!------------------------------------------------------------------------

module mod_memory

  use mod_memory_tools                        ! Tools for memory management
  use mod_memory_basic                        ! Memory for basic types defined in def_kintyp_basic.f90
#ifndef I_AM_NOT_ALYA
  use mod_memory_physics                      ! Memory for physics types defined in def_kintyp_physics.f90
  use mod_memory_parall                       ! Memory for parallelization arrays
  use mod_memory_spmat                        ! Memory for matrices
#endif

  implicit none

  private

  public :: memory_initialization             ! Initialization of the module
  public :: memory_initia                     ! Allocate memory
  public :: memory_alloca                     ! Allocate memory
  public :: memory_size                       ! Gives the size of a pointer (=0 if not associated)
  public :: memory_alloca_min                 ! Allocate a minimum memory for arrays
  public :: memory_deallo                     ! Deallocate memory
  public :: memory_copy                       ! Copy an array
  public :: memory_renumber                   ! Renumber an array
  public :: memory_resize                     ! Resize an array
  public :: memory_unit                       ! Memory scaling and units
  public :: memory_output_info                ! Info about memory
  public :: lbytm                             ! Allocated bytes
  public :: mem_curre                         ! Current memory allocation in bytes
  public :: mem_maxim                         ! Current memory allocation in bytes
  public :: lun_memor                         ! Output unit
  public :: lun_varcount                      ! Variable memory counter unit
  public :: kfl_alloc                         ! Allocation mode
  public :: mem_alloc                         ! Number of allocation
  public :: Kbytes                            ! Memory units
  public :: Mbytes                            ! Memory units
  public :: Gbytes                            ! Memory units
  public :: memory_add_to_memory_counter      ! Add bytes to memory counter
  public :: memory_remove_from_memory_counter ! Add bytes to memory counter
  public :: memory_allocation_mode            ! Set allocation mode
  public :: memory_output_variable_counter    ! Output variable counter
  
end module mod_memory
!> @}

