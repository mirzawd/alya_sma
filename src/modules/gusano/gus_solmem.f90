!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    gus_solver.f90
!> @date    14/06/2019
!> @author  Guillaume Houzeaux
!> @brief   Allocate memory
!> @details Allocate memory for solver and variables that should be
!>          reallocated and not redistributed in case of repartitioning
!> @}
!-----------------------------------------------------------------------

subroutine gus_solmem()
  
  use def_master,             only : solve
  use def_kintyp_basic,       only : ip
  use def_kintyp_dims,        only : ndime
  use def_kintyp_dims,        only : npoin
  use def_master,             only : veloc
  use def_master,             only : vel1d
  use def_master,             only : mem_modul
  use def_master,             only : modul
  use mod_gus_arrays,         only : gus_arrays
  use mod_arrays,             only : arrays
  use mod_arrays,             only : arrays_number
  use mod_memory_basic,       only : memory_alloca
  use mod_output_postprocess, only : output_postprocess_check_variable_postprocess
  use def_gusano
  
  implicit none
  
  !----------------------------------------------------------------------
  !
  ! Solver
  !
  !----------------------------------------------------------------------

  if( solve(1) % kfl_block == 0 ) then
     solve(1) % num_blocks          = 1
     solve(1) % block_dimensions(1) = 2
  else 
     solve(1) % num_blocks          = 2
     solve(1) % block_dimensions(1) = 1
     solve(1) % block_dimensions(2) = 1
  end if

  call soldef(4_ip)
  !
  ! Boundary conditions
  !
  solve(1) % bvess                      => bvess_gus(:,:,1)   ! Momentum
  solve(1) % kfl_fixno                  => kfl_fixno_gus
  !
  ! In case a block solver is used: If the SIZE(BVESS_GUS,1) is
  ! not that of the first block, then it is understodd by the
  ! solver that solve(1) % block_array(2) is useless
  !
  solve(1) % block_array(1) % bvess     => bvess_gus(:,:,1)
  solve(1) % block_array(1) % kfl_fixno => kfl_fixno_gus
  solve(1) % block_array(2) % bvess     => bvess_gus(:,:,1)
  solve(1) % block_array(2) % kfl_fixno => kfl_fixno_gus
  solve(1) % block_solve                => solve(2:)

  if( kfl_algor_gus == GUS_SCHUR_COMPLEMENT ) then
     call memory_alloca(mem_modul(1:2,modul),'SCHUR_GUS'    ,'gus_memall',schur_gus,solve(3) % nzmat)
     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXSC_GUS','gus_memall',kfl_fixsc_gus,1_ip,max(npoin,1_ip))
     solve(3) % kfl_iffix =  2
     solve(3) % kfl_fixno => kfl_fixsc_gus
  end if
  
  !----------------------------------------------------------------------
  !
  ! Secondary arrays
  !
  !----------------------------------------------------------------------

  call arrays(arrays_number('PROJM'),'ALLOCATE',projm_gus,npoin,2_ip)
  call arrays(arrays_number('PROJC'),'ALLOCATE',projc_gus,npoin,2_ip)
  call arrays(arrays_number('VELOC'),'ALLOCATE',veloc    ,ndime,npoin,1_ip)
  call arrays(arrays_number('ANGLE'),'ALLOCATE',angle_gus,npoin)
  call arrays(arrays_number('VEL1D'),'ALLOCATE',vel1d,npoin,3_ip)

  call memory_alloca(mem_modul(1:2,modul),'NEUMAN_GUS','gus_memall',neuman_gus,npoin)
  call memory_alloca(mem_modul(1:2,modul),'DIRICH_GUS','gus_memall',dirich_gus,npoin)
  call memory_alloca(mem_modul(1:2,modul),'BENDI_GUS' ,'gus_memall',bendi_gus,2_ip,npoin)
  
  if( kfl_bendi_gus /= GUS_BEND_OFF ) then 
     call memory_alloca(mem_modul(1:2,modul),'DENSI_GUS','gus_memall',densi_gus,npoin)
     call memory_alloca(mem_modul(1:2,modul),'VISCO_GUS','gus_memall',visco_gus,npoin)
  end if


end subroutine gus_solmem
