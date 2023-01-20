!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Initia
!> @{
!> @file    Turnon.f90
!> @author  Guillaume Houzeaux
!> @brief   Turnon the run
!> @details Initial operatons: read general data, mesh and modules data
!>          Construct the mesh dependent arrays.
!>
!> @}
!-----------------------------------------------------------------------
subroutine Initia

  use def_parame
  use def_elmtyp
  use def_master
  use def_inpout
  use def_domain
  use def_kermod
  use mod_bourgogne_pinotnoir
  use mod_unity_tests
  use mod_finite_volume,             only : finite_volume_arrays
  use mod_auto_tuning,               only : auto_tuning_SpMV_OpenMP
  use mod_messages,                  only : messages_report
  use mod_communications,            only : PAR_BARRIER
  use mod_messages,                  only : messages_live
  use mod_outfor,                    only : outfor
  use mod_messages,                  only : livinf
  use mod_parall_openmp,             only : parall_openmp_adjacency_ompss_unity_test
  use mod_alya2dlb,                  only : alya2dlb_initialization
#ifdef ALYA_FTI
  use mod_alya2fti,                  only : alya2fti_initialization
#endif
  use mod_alya2talp,                 only : alya2talp_register
  use mod_alya2talp,                 only : alya2talp_register_module
  use mod_alya2signal,               only : alya2signal
  use mod_optimum_partition,         only : optimum_partition_setup
  use mod_state,                     only : state_output
  use mod_alya2extrae,               only : alya2extrae_initialization
  implicit none

  integer(8) :: count_rate8
 
  !---------------------------------------------------------------------------------
  !
  ! To force Nanos to work with any MPI versions... here just in case!
  !
  !---------------------------------------------------------------------------------

#ifdef ALYA_OMPSS
  integer :: i
  !$omp do schedule(static)
  do i=1,2
  end do
  !$omp do schedule(dynamic)
  do i=1,2
  end do
#endif  

  !---------------------------------------------------------------------------------
  !
  ! Start run, creates code communicators in case of coupling, and read problem data
  !
  !---------------------------------------------------------------------------------

  call bourgogne(1_ip)
  !
  ! Compute time rate
  !
  call system_clock(count_rate=count_rate8)
  rate_time = 1.0_rp / max(real(count_rate8,rp),zeror)
  !
  ! Initialize file units, par_code_split_world needs to access them
  !
  call units_set()
  !
  ! Splits the MPI_COMM_WORLD for coupling with other codes. This defines the Alya world
  !
  call par_code_split_universe()
  !
  ! Splits the MPI_COMM_WORLD for coupling
  !
  call par_code_split_world()
  !
  ! Initialize DLB, just after initializing MPI
  !
  call alya2dlb_initialization()  
  call alya2talp_register()
  !
  ! Initialize FTI
  !
#ifdef ALYA_FTI
  call alya2fti_initialization()
#endif
  !
  ! Initializations
  !
  call inirun()
  !
  ! Read problem data
  !
  call Reapro()
  !
  ! Register TALP who needs to know used modules
  !
  call alya2talp_register_module()
  !
  ! Define tags for extrae for used modules
  !
  call alya2extrae_initialization()
  !
  ! Initialise PETSc library
  !
#ifdef PETSC
  block
  use mod_alya2petsc, only: alya2petsc_initialise
  call alya2petsc_initialise()
  endblock
#endif
  !
  ! Initialize time counters
  !
  call setgts(ITASK_INITIA)
  !
  ! Initializations
  !
  call inidom()
  !
  ! Define element types
  !
  call elmtyp()                            ! Element lists: NNODE, LTOPO, LDIME, LLAPL...
  !
  ! Open module data files
  !
  call moddef(1_ip)
  !
  ! Signal handling
  !
  call alya2signal()
  !
  ! Automatic partitioning with Pycomms
  !
  call optimum_partition_setup()
  !
  ! State output
  !
  call state_output('INITIALIZING')
  
end subroutine Initia
