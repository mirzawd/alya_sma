!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kernel
!> @{
!> @file    Turnof.f90
!> @author  houzeaux
!> @date    2019-02-27
!> @brief   Turnof 
!> @details Finish Alya... smoothly...
!> @} 
!-----------------------------------------------------------------------

subroutine Turnof
  
  use def_kintyp,               only : ip
  use def_master,               only : IPARALL
  use def_master,               only : ITASK_TURNOF
  use def_master,               only : nblok
  use def_master,               only : iblok
  use mod_par_output_partition, only : par_output_partition
  use mod_par_output_partition, only : par_output_solvers
  use mod_cou_output,           only : cou_output_timings
  use mod_bourgogne_pinotnoir,  only : bourgogne
  use mod_mpio_par_configure,   only : par_mpio_finalize
  use mod_messages,             only : messages_live
  use mod_moduls,               only : moduls
  use mod_repartitioning,       only : repartitioning_timings
  use mod_performance,          only : performance_outcpu
  use mod_performance,          only : performance_outmem
  use mod_perf_csv,             only : performance_csv
  use mod_performance,          only : performance_ann
  use mod_parall_destructor,    only : parall_destructor
  use mod_mpio_par_postpr,      only : posmpio_destructor
  use mod_mass_matrix,          only : mass_matrix_destructor
  use mod_materials,            only : materials_destructor
  use mod_exterior_normal,      only : exterior_normal_destructor

#ifdef ALYA_FTI
  use mod_alya2fti,             only : alya2fti_finalization
#endif

  implicit none
  !
  ! NINJAAAAAAAA
  !
#ifdef NINJA
  call lasterrorninja()
#endif
  call bourgogne(1_ip)
  !
  ! Write info about direct solvers
  !
  call out_direct_solvers()
  !
  ! Writes memory used
  !
  call performance_outmem()
  !
  ! Repartitioning
  !
  call repartitioning_timings()     
  !
  ! Write info about timings in partition file
  !
  call par_output_partition(2_ip)
  !
  ! Postprocess solver information
  !
  call par_output_solvers()
  !
  ! Turn off modules
  !
  do iblok = 1,nblok
     call moduls(ITASK_TURNOF)
  end do
  !
  ! Coupling CPU
  !
  call cou_output_timings()
  !
  ! Write CPU time heading and master's CPU time
  !
  call performance_outcpu()
  !
  ! Write performance file in csv format
  !
  call performance_csv()
  !
  ! Write performance of neural networks 
  !
  call performance_ann()
  !
  ! Close module files
  !
  call moddef(4_ip)
  
  if (IPARALL) then
    call par_mpio_finalize()
  end if
  call par_turnof()    
  !
  !call finalize coprocessing
  !
#ifdef CATA
  call messages_live('CALL FINALIZE COPROCESSING')
  call coprocessorfinalize()
#endif

  call bourgogne(2_ip)
  !
  ! Destroy
  !
  call coupli_destructor()
  call cou_memory(0_ip)
  call cshder(2_ip)
  call solver_destructor() 
  call domain_destructor()
  call geometry_destructor()
  call mass_matrix_destructor()
  call exterior_normal_destructor()
  call materials_destructor()
  call parall_destructor()
  call posmpio_destructor()
  !
  ! Destroy PETSc library
  !
#if PETSC
  block
  use mod_alya2petsc, only : alya2petsc_finalise
  call alya2petsc_finalise()
  endblock
#endif
  !
  ! Uncomment to check state of memory
  !
  !call performance_outmem()
  !
#if ALYA_FTI
  call alya2fti_finalization()
#endif

  ! Stop the run
  !
  call runend('O.K.!')
  
end subroutine Turnof
