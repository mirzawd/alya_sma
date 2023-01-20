!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> Brige to TALP
!> @{
!> @file    mod_alya2talp.f90
!> @author  houzeaux
!> @date    2019-01-07
!> @brief   Bridge to TALP
!> @details Interfaces with TALP
!>          Data for global timings:
!>          DLB_MONITOR_GLOBAL
!>          DLB_HANDLE_GLOBAL
!>          Data for modules:
!>          DLB_MONITOR_MODULE
!>          DLB_HANDLE_MODULE
!>
!-----------------------------------------------------------------------

module mod_alya2talp

  use def_kintyp,         only : ip,rp,lg
  use def_master,         only : mmodu,modul
  use def_master,         only : namod
  use def_master,         only : namda
  use def_master,         only : kfl_modul
  use def_master,         only : intost
  use def_master,         only : zeror
  use def_master,         only : npart
  use mod_messages,       only : messages_live
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_AVERAGE

  use, intrinsic :: ISO_C_BINDING

  implicit none
     
#if defined ALYA_TALP
  include 'dlbf_talp.h'

  type talp_module
     type(dlb_monitor_t), pointer :: monitor
     type(c_ptr)                  :: handle
  end type talp_module
      
  type(talp_module)               :: dlb_module(0:mmodu)
  type(talp_module)               :: dlb_global

#endif  

  private

  public :: alya2talp_register
  public :: alya2talp_register_module
  public :: alya2talp_MonitoringRegionStart
  public :: alya2talp_MonitoringRegionStop
  public :: alya2talp_MonitoringRegionreset
  public :: alya2talp_parallel_efficiency
  public :: alya2talp_initialization

contains

  subroutine alya2talp_initialization()

#if defined ALYA_TALP
    integer(ip) :: imodu

    nullify(dlb_global % monitor)
    do imodu = 0,mmodu
       nullify(dlb_module(imodu) % monitor)
    end do
#endif  

  end subroutine alya2talp_initialization

  !-----------------------------------------------------------------------
  !> 
  !> @author  bsc21943
  !> @date    2019-07-19
  !> @brief   Initialization
  !> @details Initialization of TALP whenevr TALP is used for load balance
  !>          or TALP_Barrier is used
  !> 
  !-----------------------------------------------------------------------

  subroutine alya2talp_register()

#if defined ALYA_TALP
    dlb_global % handle = DLB_MonitoringRegionRegister(c_char_"region global"//C_NULL_CHAR)
    if (.not. c_associated(dlb_global % handle)) call runend('MOD_ALYA2TALP: ERROR WHEN ASSOCIATING HANDLE_GLOBAL')
    call alya2talp_MonitoringRegionStart(GLOBAL_REGION=.true.)
#endif
    
  end subroutine alya2talp_register

  subroutine alya2talp_register_module()

#if defined ALYA_TALP
    integer(ip) :: imodu

    do imodu = 1,mmodu-1
       if( kfl_modul(imodu) == 1 ) then
            dlb_module(imodu) % handle = DLB_MonitoringRegionRegister(c_char_"region module "//trim(namod(imodu))//": "//trim(namda)//C_NULL_CHAR)
          if (.not. c_associated(dlb_module(imodu) % handle) ) call runend('MOD_ALYA2TALP: ERROR WHEN ASSOCIATING HANDLE_MODULE')
       end if
    end do

#endif
    
  end subroutine alya2talp_register_module
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  bsc21943
  !> @date    2020-03-06
  !> @brief   Start monitoring
  !> @details Start monitoring a region
  !> 
  !-----------------------------------------------------------------------

  subroutine alya2talp_MonitoringRegionStart(&
       DYNAMIC_ALLOCATION_REGION,&
       GLOBAL_REGION,&
       MODULE_REGION,&
       CURRENT_MODULE)

    logical(lg), optional, intent(in) :: DYNAMIC_ALLOCATION_REGION
    logical(lg), optional, intent(in) :: GLOBAL_REGION
    logical(lg), optional, intent(in) :: MODULE_REGION
    integer(ip), optional, intent(in) :: CURRENT_MODULE
#if defined ALYA_TALP
    integer(4)                        :: ierr
    logical(lg)                       :: if_global_region
    logical(lg)                       :: if_module_region
    integer(ip)                       :: imodu
#endif

#if defined ALYA_TALP
    if( present(GLOBAL_REGION) ) then
       if_global_region = GLOBAL_REGION
    else
       if_global_region = .false.
    end if
    if( present(MODULE_REGION) ) then
       if_module_region = MODULE_REGION
    else
       if_module_region = .false.
    end if
    if( present(CURRENT_MODULE) ) then
       imodu = CURRENT_MODULE
    else
       imodu = modul
    end if    

    if( if_global_region ) then
       if(c_associated(dlb_global % handle)) &
            ierr = DLB_MonitoringRegionStart(dlb_global % handle)
    end if
    if( if_module_region ) then
       if(c_associated(dlb_module(imodu) % handle)) &
            ierr = DLB_MonitoringRegionStart(dlb_module(imodu) % handle)
    end if
#endif

  end subroutine alya2talp_MonitoringRegionStart

  !-----------------------------------------------------------------------
  !> 
  !> @author  bsc21943
  !> @date    2020-03-06
  !> @brief   Stop monitoring
  !> @details Stop monitoring a region
  !> 
  !-----------------------------------------------------------------------

  subroutine alya2talp_MonitoringRegionStop(&
       DYNAMIC_ALLOCATION_REGION,&
       GLOBAL_REGION,&
       MODULE_REGION,&
       CURRENT_MODULE)

    logical(lg), optional, intent(in) :: DYNAMIC_ALLOCATION_REGION
    logical(lg), optional, intent(in) :: GLOBAL_REGION
    logical(lg), optional, intent(in) :: MODULE_REGION
    integer(ip), optional, intent(in) :: CURRENT_MODULE
#if defined ALYA_TALP
    integer(4)                        :: ierr
    logical(lg)                       :: if_global_region
    logical(lg)                       :: if_module_region
    integer(ip)                       :: imodu
#endif

#if defined ALYA_TALP
     if( present(GLOBAL_REGION) ) then
       if_global_region = GLOBAL_REGION
    else
       if_global_region = .false.
    end if
   if( present(MODULE_REGION) ) then
       if_module_region = MODULE_REGION
    else
       if_module_region = .false.
    end if
    if( present(CURRENT_MODULE) ) then
       imodu =  CURRENT_MODULE
    else
       imodu = modul
    end if    

    if( if_global_region ) then
       if(c_associated(dlb_global % handle)) &
            ierr = DLB_MonitoringRegionStop(dlb_global % handle)
    end if
    if( if_module_region ) then
       if(c_associated(dlb_module(imodu) % handle)) &
            ierr = DLB_MonitoringRegionStop(dlb_module(imodu) % handle)
    end if
#endif

  end subroutine alya2talp_MonitoringRegionStop

  !-----------------------------------------------------------------------
  !> 
  !> @author  bsc21943
  !> @date    2020-03-06
  !> @brief   Reset a region
  !> @details Reset a region: put counter to zero
  !> 
  !-----------------------------------------------------------------------

  subroutine alya2talp_MonitoringRegionReset(&
       DYNAMIC_ALLOCATION_REGION,&
       GLOBAL_REGION,&
       MODULE_REGION,&
       CURRENT_MODULE)

    logical(lg), optional, intent(in) :: DYNAMIC_ALLOCATION_REGION
    logical(lg), optional, intent(in) :: GLOBAL_REGION
    logical(lg), optional, intent(in) :: MODULE_REGION
    integer(ip), optional, intent(in) :: CURRENT_MODULE
#if defined ALYA_TALP
    integer(4)                        :: ierr
    logical(lg)                       :: if_global_region
    logical(lg)                       :: if_module_region
    integer(ip)                       :: imodu
#endif

#if defined ALYA_TALP
    if( present(GLOBAL_REGION) ) then
       if_global_region = GLOBAL_REGION
    else
       if_global_region = .false.
    end if
    if( present(MODULE_REGION) ) then
       if_module_region = MODULE_REGION
    else
       if_module_region = .false.
    end if

    if( if_global_region ) then
       ierr = DLB_MonitoringRegionreset(dlb_global % handle)
    end if
    if( if_module_region ) then
       do imodu = 0,mmodu
          if( kfl_modul(imodu) == 1 ) &
               ierr = DLB_MonitoringRegionreset(dlb_module(imodu) % handle)
       end do
    end if
#endif

  end subroutine alya2talp_MonitoringRegionReset

  !-----------------------------------------------------------------------
  !> 
  !> @author  bsc21943
  !> @date    2020-03-06
  !> @brief   Get values from a region
  !> @details Get TALP counter from a region. Time is given in ns
  !>          by TALP.
  !> 
  !-----------------------------------------------------------------------

#if defined ALYA_TALP
  subroutine alya2talp_time(elapsed_time,accumulated_MPI_time,accumulated_computation_time,handle,monitor)
    
    real(rp),            optional, intent(out) :: elapsed_time
    real(rp),            optional, intent(out) :: accumulated_MPI_time
    real(rp),            optional, intent(out) :: accumulated_computation_time
    type(c_ptr),         optional              :: handle
    type(dlb_monitor_t), optional, pointer     :: monitor

    if( present(handle) .and. present(monitor) ) then
       call c_f_pointer(handle, monitor)
       if( present(elapsed_time) )                 elapsed_time                 = real(monitor % elapsed_time,rp)*1.0e-9_rp
       if( present(accumulated_MPI_time) )         accumulated_MPI_time         = real(monitor % accumulated_MPI_time,rp)*1.0e-9_rp
       if( present(accumulated_computation_time) ) accumulated_computation_time = real(monitor % accumulated_computation_time,rp)*1.0e-9_rp
    end if

  end subroutine alya2talp_time
#endif

  !-----------------------------------------------------------------------
  !> 
  !> @author  bsc21943
  !> @date    2020-03-06
  !> @brief   Performance indicators
  !> @details Parallel peformance indicators
  !>
  !>            Definitions taken from PoP project: https://pop-coe.eu/node/69 
  !>                                                                           
  !>                t_i^w      t_i^MPI                                         
  !>            <------------><------->                                        
  !>                                                                           
  !>            +---------------------+  -+                                    
  !>            |\\\\\\\\\\\\\        +   |                                    
  !>            +---------------------+   |                                    
  !>            |\\\\\\\\\\\\\\\\\    +   | P processors                       
  !>            +---------------------+   |                                    
  !>            |\\\\\\\              +   |                                    
  !>            +---------------------+  -+                                    
  !>                                                                           
  !>            <---------------------> te                                     
  !>                                                                           
  !>            t_i^w ..... working time of process i                          
  !>            t_i^MPI ... MPI communication time of process i                
  !>            t_e ....... Elpased time (same for all)                        
  !>                                                                           
  !>            t_ave^w = sum_i t_i^w / P (average working time)               
  !>            t_max^w = max_i t_i^w     (max working time)                   
  !>                                                                           
  !>                                                                           
  !>            Load balance:             LB = t_ave^w /t_max^w                
  !>            Communication efficiency: CE = t_max^w /t_e                    
  !>            Parallel efficiency:      PE = LB * CE = t_ave^w / te          
  !>
  !-----------------------------------------------------------------------
  
  subroutine alya2talp_parallel_efficiency(&
       PE,&
       LB,&
       CE,&
       time_comp_ave,&
       time_mpi_ave,&
       time_comp_max,&
       time_mpi_max,&
       time_total,&
       GLOBAL_REGION,&
       DYNAMIC_ALLOCATION_REGION,&
       MODULE_REGION,&
       CURRENT_MODULE)
    
    real(rp),              intent(out) :: PE                         !< Parallel efficiency
    real(rp),              intent(out) :: LB                         !< Load balance
    real(rp),              intent(out) :: CE                         !< Communication efficiency
    real(rp),    optional, intent(out) :: time_comp_ave
    real(rp),    optional, intent(out) :: time_mpi_ave
    real(rp),    optional, intent(out) :: time_comp_max
    real(rp),    optional, intent(out) :: time_mpi_max
    real(rp),    optional, intent(out) :: time_total
    logical(lg), optional, intent(in)  :: GLOBAL_REGION
    logical(lg), optional, intent(in)  :: DYNAMIC_ALLOCATION_REGION
    logical(lg), optional, intent(in)  :: MODULE_REGION
    integer(ip), optional, intent(in)  :: CURRENT_MODULE
    integer(ip)                        :: imodu
    logical(lg)                        :: if_global_region
    logical(lg)                        :: if_dynamic_allocation_region
    logical(lg)                        :: if_module_region
#if defined ALYA_TALP
    real(rp)                           :: elapsed_time
    real(rp)                           :: accumulated_MPI_time
    real(rp)                           :: accumulated_computation_time
    real(rp)                           :: ave_elapsed_time
    real(rp)                           :: max_elapsed_time
    real(rp)                           :: max_accumulated_MPI_time
    real(rp)                           :: ave_accumulated_MPI_time
    real(rp)                           :: ave_accumulated_computation_time
    real(rp)                           :: max_accumulated_computation_time
    real(rp)                           :: P
#endif

    if( present(GLOBAL_REGION) ) then
       if_global_region = GLOBAL_REGION
    else
       if_global_region = .false.
    end if
    if( present(DYNAMIC_ALLOCATION_REGION) ) then
       if_dynamic_allocation_region = DYNAMIC_ALLOCATION_REGION
    else
       if_dynamic_allocation_region = .false.
    end if
    if( present(MODULE_REGION) ) then
       if_module_region = MODULE_REGION
    else
       if_module_region = .false.
    end if
    if( present(CURRENT_MODULE) ) then
       imodu = CURRENT_MODULE
    else
       imodu = modul
    end if    

#if defined ALYA_TALP

    if( if_global_region ) then
       call alya2talp_time(elapsed_time,accumulated_MPI_time,accumulated_computation_time,&
            dlb_global % handle,&
            dlb_global % monitor)
    else if( if_module_region ) then
       call alya2talp_time(elapsed_time,accumulated_MPI_time,accumulated_computation_time,&
            dlb_module(imodu) % handle,&
            dlb_module(imodu) % monitor)
    else
       call runend('MOD_ALYA2TAP: DO NOT KNOW WHAT TO DO')
    end if

     ave_accumulated_MPI_time         = accumulated_MPI_time
     max_accumulated_MPI_time         = accumulated_MPI_time
     ave_accumulated_computation_time = accumulated_computation_time
     max_accumulated_computation_time = accumulated_computation_time
     ave_elapsed_time                 = elapsed_time
     max_elapsed_time                 = elapsed_time

     call PAR_AVERAGE(ave_accumulated_MPI_time)
     call PAR_MAX    (max_accumulated_MPI_time)
     call PAR_AVERAGE(ave_accumulated_computation_time)
     call PAR_MAX    (max_accumulated_computation_time)
     call PAR_AVERAGE(ave_elapsed_time)
     call PAR_MAX    (max_elapsed_time)
     
     LB = ave_accumulated_computation_time / (max_accumulated_computation_time+zeror) ! Load balance
     CE = max_accumulated_computation_time / (max_elapsed_time+zeror)                 ! Communication efficiency
     PE = ave_accumulated_computation_time / (max_elapsed_time+zeror)                 ! Parallel efficiency

     if( present(time_comp_ave)  ) time_comp_ave  = ave_accumulated_computation_time
     if( present(time_comp_max)  ) time_comp_max  = max_accumulated_computation_time
     if( present(time_mpi_ave)   ) time_mpi_ave   = ave_accumulated_MPI_time
     if( present(time_mpi_max)   ) time_mpi_max   = max_accumulated_MPI_time
     if( present(time_total)     ) time_total     = max_elapsed_time

#else

     PE = -1.0_rp
     LB = -1.0_rp
     CE = -1.0_rp
     if( present(time_comp_ave)  ) time_comp_ave  = -1.0_rp
     if( present(time_comp_max)  ) time_comp_max  = -1.0_rp
     if( present(time_mpi_ave)   ) time_mpi_ave   = -1.0_rp
     if( present(time_mpi_max)   ) time_mpi_max   = -1.0_rp
     if( present(time_total)     ) time_total     = -1.0_rp

#endif

  end subroutine alya2talp_parallel_efficiency

end module mod_alya2talp
!> @}
