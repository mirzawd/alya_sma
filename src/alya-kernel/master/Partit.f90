!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Partit
!> @{
!> @file    Partit.f90
!> @author  houzeaux
!> @date    2019-06-11
!> @brief   Partition domain
!> @details Partition domain
!> @} 
!-----------------------------------------------------------------------

subroutine Partit()

  use def_master
  use def_domain
  use mod_ker_timeline,     only : ker_timeline
  use def_parall,           only : kfl_partition_par
  use def_parall,           only : kfl_parseq_par
  use def_parall,           only : kfl_virfi_par
  use mod_parall,           only : PAR_USING_RANK
  use mod_parall,           only : PAR_PARALLEL_PARTITION
  use mod_domain,           only : domain_memory_deallocate 
  use mod_par_partitioning, only : par_partitioning
  use mod_par_virfil,       only : par_dumbuf
  use mod_performance,      only : performance_outcpu
  use mod_performance,      only : performance_outmem
  use mod_mpio_config,      only : mpio_config
  use mod_run_config,       only : run_config
  use mod_one_sided_communications, only : par_one_sided_allocate

  implicit none
  real(rp) :: time1,time2

  if( IPARALL ) then
     !
     ! When exporting mesh in parallel, do not partition neither redistribute
     !
     if( mpio_config%output%post_process%export_only ) then
        kfl_partition_par = PAR_USING_RANK
        kfl_parseq_par    = PAR_PARALLEL_PARTITION
     end if
     call cputim(time1)
     if( run_config%timeline ) call ker_timeline(0_ip,'INI_PARTITION_MESH')
     !
     ! Partition mesh
     !
     call par_partitioning()

     call par_errors(2_ip)
     call par_openfi(2_ip)
     !
     ! Output partition info
     !
     if( kfl_virfi_par == 1 ) then                       ! Virtual file
        call par_dumbuf(-1_ip)
        kfl_virfi_par = 0
     end if
     !
     ! Info and possibly end the run
     !
     if( PART_AND_WRITE() ) then                         ! Partition only: end of the run
        call performance_outmem()
        call performance_outcpu() 
        call par_turnof()
        call runend('O.K.!')
     end if
     if( IMASTER .and. .not. READ_AND_RUN() ) then
        call domain_memory_deallocate('ALL MESH')        ! Deallocate Master geometry memory
        call domain_memory_deallocate('LESET')           ! Sets
        call domain_memory_deallocate('LBSET')           ! Sets
        call domain_memory_deallocate('LNSET')           ! Sets
     end if
     if( IMASTER ) then
        call par_memory(4_ip)                            ! Deallocate memory of partition arrays
     end if
     kfl_ptask = 1                                       ! Switch to normal execution
     call vocabu(-1_ip,0_ip,0_ip)
     if( run_config%timeline ) call ker_timeline(0_ip,'END_PARTITION_MESH')
     call cputim(time2)
     cpu_start(CPU_MESH_PARTITION) = time2 - time1
     !
     ! One-sided communications
     !
     call par_one_sided_allocate()

  end if

end subroutine Partit
