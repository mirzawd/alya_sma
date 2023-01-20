!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_initialize_nullify.f90
!> @author  Guillaume Houzeaux
!> @date    01/02/2014
!> @brief   Intialize variables
!> @details Initialize and nullify variables of mod_parall
!> @} 
!----------------------------------------------------------------------

subroutine par_initialize_nullify()
  use def_kintyp,            only : ip,rp
  use def_parall
  use mod_parall
  use def_mpio
  use mod_par_bin_structure, only : par_bin_structure_initialization
  implicit none
  !
  ! Bin structure
  !
  call par_bin_structure_initialization()
  !
  ! Initialize variables
  !
  mcolo         = 0  
  mcode         = 0  
  mzone         = 0  
  msubd         = 0
  mcode         = 1
  ncolo         = 0                         
  par_memor     = 0
  par_bin_boxes = 1
  par_bin_comin = 0.0_rp
  par_bin_comax = 0.0_rp
  PAR_CODE_SIZE = 1
  PAR_MY_CODE_RANK = 1
  num_repart_par   = 0
  !
  ! Nullify global pointers 
  !
  nullify(PAR_COMM_COLOR)           
  nullify(PAR_COMM_COLOR_PERM)    
  nullify(I_AM_IN_COLOR)              

  nullify(PAR_COMM_COLOR_ARRAY)       
  nullify(PAR_COMM_MY_CODE_ARRAY)

  nullify(PAR_CPU_TO_COLOR)           
  nullify(PAR_COLOR_TO_CPU)           
  nullify(PAR_COLOR_BIN)              

  nullify(par_bin_part)
  nullify(par_bin_size)
  nullify(par_part_comin)
  nullify(par_part_comax)
  !
  ! OpenMP/OmpSs/vectorization stuffs
  !
  par_hybrid                 = PAR_HYBRID_OFF
  par_omp_num_blocks         = 1
  par_omp_granularity        = 10
  par_omp_coloring_alg       = 1
  par_omp_num_threads        = 0
  
  par_omp_num_colors         = 0
  par_omp_nboun_num_colors   = 0
  
  par_omp_nelem_chunk        = 0  
  par_omp_nboun_chunk        = 0
  par_omp_npoin_chunk        = 0
  num_subd_par               = 0       
  num_subd_norace_par        = 0   
  num_subd_cpu               = 0       
  num_subd_norace_cpu        = 0   
  num_subd_nboun_par         = 0       
  num_subd_norace_nboun_par  = 0   

  nullify(par_omp_list_colors) 
  nullify(par_omp_ia_colors)
  nullify(par_omp_ja_colors)
  
  nullify(par_omp_nboun_list_colors) 
  nullify(par_omp_nboun_ia_colors)
  nullify(par_omp_nboun_ja_colors)
  
  nullify(num_pack_par)
  nullify(list_elements_par)

  nullify(num_pack_norace_par)
  nullify(list_elements_norace_par)
  
  nullify(num_pack_cpu)
  nullify(list_elements_cpu)

  nullify(num_pack_norace_cpu)
  nullify(list_elements_norace_cpu)
  
  nullify(num_pack_nboun_par)
  nullify(list_boundaries_par)

  nullify(num_pack_norace_nboun_par)
  nullify(list_boundaries_norace_par)

  nullify(sendbuff_rp)
  nullify(recvbuff_rp)
  nullify(sendbuff_ip)
  nullify(recvbuff_ip)
  !
  ! Partition type
  !
  kfl_partition_par   = PAR_METIS4                  ! Partition method
  boxes_coarse_par    = 0                           ! Number of boxes coarse bin
  boxes_fine_par      = 0                           ! Number of boxes fine bin
  vect_partition_par  = (/ 1.0_rp,0.0_rp,0.0_rp /)  ! Direction of partition for oriented bin
  vect_partition_par  = 0.0_rp 

end subroutine par_initialize_nullify
