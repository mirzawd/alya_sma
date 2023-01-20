!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    mod_pts_injection.f90
!> @author  bsc21240
!> @date    2018-03-08
!> @brief   Module for injection
!> @details Inject particles
!> For the circle injector you can add two more options:
!> DISTRIBUTION: UNICA|UNIPO  (uniform spatially distrib. particles in cartesian coordinates 
!> or polar coordinates[this one is default previous implementation])
!> NPARTICLES: ASIS -- by default the number of particles provided in the injector parameters is squared, cubed or 
!> transformed in some other way. This parameter disables these transformations (only for CIRCLE)
!-----------------------------------------------------------------------

module mod_pts_injection
    use def_parame
    use def_master
    use def_kermod
    use def_partis
    use def_domain
    use mod_memory,           only : memory_alloca
    use mod_memory,           only : memory_deallo
    use mod_memory,           only : memory_copy
    use mod_communications,   only : PAR_BROADCAST
    use mod_parall,           only : PAR_COMM_MY_CODE_WM
    use mod_messages,         only : livinf
    use mod_messages,         only : messages_live
    use mod_strings,          only : integer_to_string
    use mod_pts_host_element, only : pts_host_element
    use mod_pts_host_element, only : pts_host_element_parall
    implicit none

    !*******************!
    ! PARAMETERS        !
    !*******************!
    
    integer(ip),   parameter :: &
        PTS_INJ_SIZEDIST_CONST        = 0_ip,   &
        PTS_INJ_SIZEDIST_UNIFORM      = 1_ip,   & 
        PTS_INJ_SIZEDIST_ROSINRAMMLER = 2_ip,   &
        PTS_INJ_SIZEDIST_TABULATED    = 3_ip,   &
        PTS_INJ_SIZEDIST_FILED        = 4_ip,   &
        PTS_INJ_SIZEDIST_FILEM        = 5_ip     

    integer(ip),   parameter :: &
        PTS_INJ_TEMPE_CONST           = 0_ip,   &
        PTS_INJ_TEMPE_FILE            = 1_ip

    integer(ip),   parameter :: &
        PTS_INJ_VELOC_ZERO            =-1_ip,   &
        PTS_INJ_VELOC_FLUID           = 0_ip,   &
        PTS_INJ_VELOC_NORMAL          = 1_ip,   &
        PTS_INJ_VELOC_GAUSSIAN        = 2_ip,   &
        PTS_INJ_VELOC_CONIC           = 3_ip,   &
        PTS_INJ_VELOC_SPRAY           = 4_ip,   &
        PTS_INJ_VELOC_CONST           = 5_ip,   &
        PTS_INJ_VELOC_AXSPRAY         = 6_ip,   &
        PTS_INJ_VELOC_FILE            = 7_ip,   &
        PTS_INJ_VELOC_SIZE_DEP_SPRAY  = 8_ip        ! Conditional droplet injection model of Ma (2017)

    integer(ip),   parameter :: &
        PTS_INJ_FLUCT_VELOC_CONST     = 0_ip,   &
        PTS_INJ_FLUCT_VELOC_UNIFORM   = 1_ip,   &
        PTS_INJ_FLUCT_VELOC_NORMAL    = 2_ip

    integer(ip),   parameter :: &
        PTS_INJ_FLOW_NMODIFIED        = 0_ip,   &
        PTS_INJ_FLOW_NASIS            = 1_ip,   &
        PTS_INJ_FLOW_MASSFLOW         = 2_ip,   &
        PTS_INJ_FLOW_VOLUMEFLOW       = 3_ip     

    integer(ip),   parameter :: &
        PTS_INJ_SPATDIST_UNICARTESIAN = 1_ip,   &
        PTS_INJ_SPATDIST_UNIPOLAR     = 2_ip     

    integer(ip),   parameter :: &
        PTS_INJ_GEO_SQUARE            = 1_ip,   &
        PTS_INJ_GEO_SPHERE            = 2_ip,   &
        PTS_INJ_GEO_SEMISPHERE        = 3_ip,   &
        PTS_INJ_GEO_CIRCLE            = 4_ip,   &
        PTS_INJ_GEO_RECTANGLE         = 5_ip,   &
        PTS_INJ_GEO_POINT             = 6_ip,   &
        PTS_INJ_GEO_CONE              = 7_ip,   &
        PTS_INJ_GEO_RANDOM            = 8_ip,   &
        PTS_INJ_GEO_SEGMENT           = 9_ip,   &
        PTS_INJ_GEO_FILE              =10_ip,   &
        PTS_INJ_GEO_ANNULUS           =11_ip

    private

    public :: PTS_INJ_SIZEDIST_CONST        
    public :: PTS_INJ_SIZEDIST_UNIFORM      
    public :: PTS_INJ_SIZEDIST_ROSINRAMMLER 
    public :: PTS_INJ_SIZEDIST_TABULATED 
    public :: PTS_INJ_SIZEDIST_FILED 
    public :: PTS_INJ_SIZEDIST_FILEM 

    public :: PTS_INJ_TEMPE_CONST 
    public :: PTS_INJ_TEMPE_FILE 

    public :: PTS_INJ_VELOC_ZERO      
    public :: PTS_INJ_VELOC_FLUID     
    public :: PTS_INJ_VELOC_NORMAL    
    public :: PTS_INJ_VELOC_GAUSSIAN  
    public :: PTS_INJ_VELOC_CONIC     
    public :: PTS_INJ_VELOC_SPRAY     
    public :: PTS_INJ_VELOC_CONST     
    public :: PTS_INJ_VELOC_AXSPRAY     
    public :: PTS_INJ_VELOC_FILE 
    public :: PTS_INJ_VELOC_SIZE_DEP_SPRAY 

    public :: PTS_INJ_FLUCT_VELOC_CONST     
    public :: PTS_INJ_FLUCT_VELOC_UNIFORM   
    public :: PTS_INJ_FLUCT_VELOC_NORMAL    

    public :: PTS_INJ_FLOW_NMODIFIED 
    public :: PTS_INJ_FLOW_NASIS     
    public :: PTS_INJ_FLOW_MASSFLOW  
    public :: PTS_INJ_FLOW_VOLUMEFLOW

    public :: PTS_INJ_SPATDIST_UNICARTESIAN
    public :: PTS_INJ_SPATDIST_UNIPOLAR    

    public :: PTS_INJ_GEO_SQUARE    
    public :: PTS_INJ_GEO_SPHERE    
    public :: PTS_INJ_GEO_SEMISPHERE
    public :: PTS_INJ_GEO_CIRCLE    
    public :: PTS_INJ_GEO_RECTANGLE 
    public :: PTS_INJ_GEO_POINT     
    public :: PTS_INJ_GEO_CONE      
    public :: PTS_INJ_GEO_RANDOM    
    public :: PTS_INJ_GEO_SEGMENT   
    public :: PTS_INJ_GEO_FILE      
    public :: PTS_INJ_GEO_ANNULUS   

    public :: pts_injection_initialization
    public :: pts_injection_parallelization
    public :: pts_injection_memory
    public :: pts_injection_injectors

contains

    !-----------------------------------------------------------------------
    !> 
    !> @author  bsc21240
    !> @date    2018-03-08
    !> @brief   Initialize injectors
    !> @details Initialize injectors
    !> 
    !-----------------------------------------------------------------------

    subroutine pts_injection_parallelization()
       use def_master,         only : parii, npari, nparr, &
                                      nparc, parin, parre, &
                                      nparc, mem_modul, &
!                                      parch, &
                                      modul, ISLAVE, IMASTER

       use mod_memory,         only : memory_alloca
       use mod_memory,         only : memory_deallo
       use mod_communications, only : PAR_BROADCAST
       use mod_communications, only : PAR_EXCHANGE
       implicit none
       integer(ip) :: iinj, ii
       character(20), parameter :: vacal = 'pts_injection_parallelization'

       !
       ! Exchange non-allocated variables
       !                                 
       do parii=1,2                      
          npari=0                        
          nparr=0                        
          nparc=0                        

          
          do iinj = 1, kfl_imax_injector_pts
             call PAR_EXCHANGE(injection_pts(iinj) % kfl_particle_type,     parin,npari,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % kfl_size_dist,         parin,npari,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % kfl_size_lookupfw,     parin,npari,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % size_diame,            parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % size_dmin,             parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % size_dmax,             parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % size_rr_dbar,          parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % size_rr_n,             parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % size_rr_K,             parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % kfl_tempe,             parin,npari,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % tempe,                 parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % kfl_veloc,             parin,npari,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % kfl_velfun,            parin,npari,parii) 
             call PAR_EXCHANGE(ndime, injection_pts(iinj) % vel_vec,        parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % vel_mag,               parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % vel_ax,                parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % vel_sigma,             parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % vel_angle,             parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % vel_mag_s,             parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % vel_diam_L,            parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % vel_diam_s,            parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % vel_ang_max_L,         parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % vel_ang_min_L,         parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % vel_ang_max_s,         parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % vel_ang_min_s,         parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % kfl_fluct_veloc,       parin,npari,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % fluct_vel_std,         parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % kfl_flow,              parin,npari,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % kfl_flowfun,           parin,npari,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % flow_rate,             parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % num_part,              parin,npari,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % kfl_geometry,          parin,npari,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % kfl_geo_spatial_dist,  parin,npari,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % kfl_random          ,  parin,npari,parii) 
             call PAR_EXCHANGE(ndime, injection_pts(iinj) % geo_coord_min,  parre,nparr,parii) 
             call PAR_EXCHANGE(ndime, injection_pts(iinj) % geo_coord_max,  parre,nparr,parii) 
             call PAR_EXCHANGE(ndime, injection_pts(iinj) % geo_coord    ,  parre,nparr,parii) 
             call PAR_EXCHANGE(ndime, injection_pts(iinj) % geo_normal   ,  parre,nparr,parii) 
             call PAR_EXCHANGE(ndime, injection_pts(iinj) % geo_basis(:,1), parre,nparr,parii) 
             call PAR_EXCHANGE(ndime, injection_pts(iinj) % geo_basis(:,2), parre,nparr,parii) 
             call PAR_EXCHANGE(ndime, injection_pts(iinj) % geo_basis(:,3), parre,nparr,parii) 
             call PAR_EXCHANGE(ndime, injection_pts(iinj) % geo_coord1   ,  parre,nparr,parii) 
             call PAR_EXCHANGE(ndime, injection_pts(iinj) % geo_coord2   ,  parre,nparr,parii) 
             call PAR_EXCHANGE(ndime, injection_pts(iinj) % geo_coord3   ,  parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % geo_rad,               parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % geo_radmin,            parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % geo_height,            parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % time_initial,          parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % time_period,           parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % time_final,            parre,nparr,parii) 
             call PAR_EXCHANGE(injection_pts(iinj) % time_cumulative,       parre,nparr,parii) 
          enddo

          if( parii == 1 ) then
             call memory_alloca(mem_modul(1:2,modul),'PARIN',vacal,parin,npari)
             call memory_alloca(mem_modul(1:2,modul),'PARRE',vacal,parre,nparr)
             if( ISLAVE  ) call PAR_BROADCAST(parin,      'IN MY CODE')
             if( ISLAVE  ) call PAR_BROADCAST(parre,      'IN MY CODE')
          else
             if( IMASTER ) call PAR_BROADCAST(parin,      'IN MY CODE')
             if( IMASTER ) call PAR_BROADCAST(parre,      'IN MY CODE')
          end if
       enddo
       call memory_deallo(mem_modul(1:2,modul),'PARIN',vacal,parin)
       call memory_deallo(mem_modul(1:2,modul),'PARRE',vacal,parre)

       !
       ! Allocate memory for allocated variables
       !
       if ( ISLAVE ) then
          do iinj = 1, kfl_imax_injector_pts
             call pts_injection_memory(iinj)
          enddo
       endif
       
       !
       ! Exchange allocated variables
       !
       do parii=1,2 
          npari=0
          nparr=0
          nparc=0

          
          do iinj = 1, kfl_imax_injector_pts
             if ( associated(injection_pts(iinj) % size_list) ) then
                call PAR_EXCHANGE(injection_pts(iinj) % num_part,injection_pts(iinj) % size_list, parre,nparr,parii)
             endif
             if ( associated(injection_pts(iinj) % tempe_list) ) then
                call PAR_EXCHANGE(injection_pts(iinj) % num_part,injection_pts(iinj) % tempe_list, parre,nparr,parii)
             endif
             if ( associated(injection_pts(iinj) % vel_list) ) then
                do ii = 1, injection_pts(iinj) % num_part
                   call PAR_EXCHANGE(ndime,injection_pts(iinj) % vel_list(1:ndime,ii), parre,nparr,parii)
                enddo
             endif
             if ( associated(injection_pts(iinj) % coord_list) ) then
                do ii = 1, injection_pts(iinj) % num_part
                   call PAR_EXCHANGE(ndime,injection_pts(iinj) % coord_list(1:ndime,ii), parre,nparr,parii)
                enddo
             endif
          enddo
          if( parii == 1 ) then
             call memory_alloca(mem_modul(1:2,modul),'PARRE',vacal,parre,nparr)
             if( ISLAVE  ) call PAR_BROADCAST(parre,      'IN MY CODE')
          else
             if( IMASTER ) call PAR_BROADCAST(parre,      'IN MY CODE')
          end if
       enddo
       call memory_deallo(mem_modul(1:2,modul),'PARRE',vacal,parre)

       !
       ! Deallocate memory for master
       !
       if (IMASTER) then
          do iinj = 1, kfl_imax_injector_pts
             call memory_deallo(mem_modul(1:2, modul), 'INJECTION_PTS % VEL_LIST',   vacal, injection_pts(iinj) % vel_list)
             call memory_deallo(mem_modul(1:2, modul), 'INJECTION_PTS % COORD_LIST', vacal, injection_pts(iinj) % coord_list)
          enddo
       endif

    end subroutine pts_injection_parallelization

    !-----------------------------------------------------------------------
    !> 
    !> @author  bsc21240
    !> @date    2018-03-08
    !> @brief   Allocate memory
    !> @details Allocate memory
    !> 
    !-----------------------------------------------------------------------

    subroutine pts_injection_memory(iinj,itask)
       integer(ip), intent(in)           :: iinj !< Injector number
       integer(ip), intent(in), optional :: itask  !< Task
       integer(ip)                       :: itask_loc

       itask_loc = 0_ip
       if (present(itask)) itask_loc = itask

       !
       ! Allocate memory if necessary
       !
       if (injection_pts(iinj) % num_part > 0 .and. &
           & (injection_pts(iinj) % kfl_size_dist == PTS_INJ_SIZEDIST_FILED .or. &
           &  injection_pts(iinj) % kfl_size_dist == PTS_INJ_SIZEDIST_FILEM) .and. &
           &  (itask_loc == 0 .or. itask_loc == 1 )) then
           call memory_alloca(mem_modul(1:2, modul), 'INJECTION_PTS % SIZE_LIST', 'pts_injection_memory', injection_pts(iinj) % size_list, injection_pts(iinj) % num_part)
       end if

       if (injection_pts(iinj) % num_part > 0 .and. &
           & injection_pts(iinj) % kfl_tempe == PTS_INJ_TEMPE_FILE .and. &
           & (itask_loc == 0 .or. itask_loc == 2 )) then
           call memory_alloca(mem_modul(1:2, modul), 'INJECTION_PTS % TEMPE_LIST', 'pts_injection_memory', injection_pts(iinj) % tempe_list, injection_pts(iinj) % num_part)
       end if

       if (injection_pts(iinj) % num_part > 0 .and. &
           & injection_pts(iinj) % kfl_veloc == PTS_INJ_VELOC_FILE .and. &
           & (itask_loc == 0 .or. itask_loc == 3 )) then
           call memory_alloca(mem_modul(1:2, modul), 'INJECTION_PTS % VEL_LIST', 'pts_injection_memory', injection_pts(iinj) % vel_list, ndime, injection_pts(iinj) % num_part)
       end if

       if (injection_pts(iinj) % num_part > 0 .and. &
           & injection_pts(iinj) % kfl_geometry == PTS_INJ_GEO_FILE .and. &
           & (itask_loc == 0 .or. itask_loc == 4 )) then
           call memory_alloca(mem_modul(1:2, modul), 'INJECTION_PTS % COORD_LIST', 'pts_injection_memory', injection_pts(iinj) % coord_list, ndime, injection_pts(iinj) % num_part)
       end if

    end subroutine pts_injection_memory

    !-----------------------------------------------------------------------
    !> 
    !> @author  bsc21240
    !> @date    2018-03-08
    !> @brief   Initialize injectors
    !> @details Initialize injectors
    !> 
    !-----------------------------------------------------------------------

    subroutine pts_injection_initialization()
       implicit none
       integer(ip) :: ii

       injection_pts(:) % kfl_particle_type    = 0_ip

       injection_pts(:) % kfl_size_dist        = PTS_INJ_SIZEDIST_CONST
       injection_pts(:) % kfl_size_lookupfw    = 0_ip
       injection_pts(:) % size_diame           = 0.0_rp
       injection_pts(:) % size_dmin            = 0.0_rp
       injection_pts(:) % size_dmax            = 0.0_rp
       injection_pts(:) % size_rr_dbar         = 0.0_rp
       injection_pts(:) % size_rr_n            = 0.0_rp
       injection_pts(:) % size_rr_K            = 0.0_rp

       injection_pts(:) % kfl_tempe            = PTS_INJ_TEMPE_CONST
       injection_pts(:) % tempe                = 200.0_rp

       injection_pts(:) % kfl_veloc            = PTS_INJ_VELOC_FLUID
       injection_pts(:) % kfl_velfun           = 0_ip 
       injection_pts(:) % vel_vec(1)           = 0.0_rp
       injection_pts(:) % vel_vec(2)           = 0.0_rp
       injection_pts(:) % vel_vec(3)           = 0.0_rp
       injection_pts(:) % vel_mag              = 0.0_rp
       injection_pts(:) % vel_ax               = 0.0_rp
       injection_pts(:) % vel_sigma            = 0.0_rp
       injection_pts(:) % vel_angle            = 0.0_rp
       injection_pts(:) % vel_mag_s            = 0.0_rp
       injection_pts(:) % vel_diam_L           = 0.0_rp
       injection_pts(:) % vel_diam_s           = 0.0_rp
       injection_pts(:) % vel_ang_max_L        = 0.0_rp
       injection_pts(:) % vel_ang_min_L        = 0.0_rp
       injection_pts(:) % vel_ang_max_s        = 0.0_rp
       injection_pts(:) % vel_ang_min_s        = 0.0_rp

       injection_pts(:) % kfl_fluct_veloc      = PTS_INJ_FLUCT_VELOC_CONST
       injection_pts(:) % fluct_vel_std        = 0.0_rp

       injection_pts(:) % kfl_flow             = PTS_INJ_FLOW_NMODIFIED
       injection_pts(:) % kfl_flowfun          = 0_ip 
       injection_pts(:) % flow_rate            = 0.0_rp
       injection_pts(:) % num_part             = 0_ip

       injection_pts(:) % kfl_geometry         = 0_ip
       injection_pts(:) % kfl_geo_spatial_dist = PTS_INJ_SPATDIST_UNIPOLAR 
       injection_pts(:) % kfl_random           = 0_ip
       injection_pts(:) % geo_coord_min(1)     = 0.0_rp
       injection_pts(:) % geo_coord_max(1)     = 0.0_rp
       injection_pts(:) % geo_coord(1)         = 0.0_rp
       injection_pts(:) % geo_normal(1)        = 0.0_rp
       injection_pts(:) % geo_coord1(1)        = 0.0_rp
       injection_pts(:) % geo_coord2(1)        = 0.0_rp
       injection_pts(:) % geo_coord3(1)        = 0.0_rp
       injection_pts(:) % geo_coord_min(2)     = 0.0_rp
       injection_pts(:) % geo_coord_max(2)     = 0.0_rp
       injection_pts(:) % geo_coord(2)         = 0.0_rp
       injection_pts(:) % geo_normal(2)        = 0.0_rp
       injection_pts(:) % geo_coord1(2)        = 0.0_rp
       injection_pts(:) % geo_coord2(2)        = 0.0_rp
       injection_pts(:) % geo_coord3(2)        = 0.0_rp
       injection_pts(:) % geo_coord_min(3)     = 0.0_rp
       injection_pts(:) % geo_coord_max(3)     = 0.0_rp
       injection_pts(:) % geo_coord(3)         = 0.0_rp
       injection_pts(:) % geo_normal(3)        = 0.0_rp
       injection_pts(:) % geo_coord1(3)        = 0.0_rp
       injection_pts(:) % geo_coord2(3)        = 0.0_rp
       injection_pts(:) % geo_coord3(3)        = 0.0_rp
       injection_pts(:) % geo_rad              = 0.0_rp
       injection_pts(:) % geo_radmin           = 0.0_rp
       injection_pts(:) % geo_height           = 0.0_rp
       injection_pts(:) % geo_basis(1,1)       = 0.0_rp
       injection_pts(:) % geo_basis(2,1)       = 0.0_rp
       injection_pts(:) % geo_basis(3,1)       = 0.0_rp
       injection_pts(:) % geo_basis(1,2)       = 0.0_rp
       injection_pts(:) % geo_basis(2,2)       = 0.0_rp
       injection_pts(:) % geo_basis(3,2)       = 0.0_rp
       injection_pts(:) % geo_basis(1,3)       = 0.0_rp
       injection_pts(:) % geo_basis(2,3)       = 0.0_rp
       injection_pts(:) % geo_basis(3,3)       = 0.0_rp
       injection_pts(:) % time_initial         = -1.0_rp
       injection_pts(:) % time_period          = -1.0_rp
       injection_pts(:) % time_final           = -1.0_rp
       injection_pts(:) % time_cumulative      = 1.0e12_rp

       do ii = 1, mpala
           nullify(injection_pts(ii) % size_list)
           nullify(injection_pts(ii) % tempe_list)
           nullify(injection_pts(ii) % vel_list)
           nullify(injection_pts(ii) % coord_list)
       end do

    end subroutine pts_injection_initialization




    function pts_injection_get_next_diame(iinj) result(diame)
       use mod_random,            only : random_generate_number
       use mod_interp_tab,        only : fw_lookup
       implicit none
       integer(ip),       intent(in)  :: iinj
       real(rp)                       :: diame
       real(rp)                       :: U1
       real(rp)                       :: control(1), scale_control(1), retva(1)
       integer(ip)                    :: ind(1)
            
       select case(injection_pts(iinj) % kfl_size_dist) 
         case(PTS_INJ_SIZEDIST_CONST)
            !
            ! CONST
            !
            diame =  injection_pts(iinj) % size_diame
         case(PTS_INJ_SIZEDIST_UNIFORM)
            !
            ! UNIFORM
            !
            U1      = random_generate_number(broadcast_seed=.true.)
            diame   = injection_pts(iinj) % size_dmin + (injection_pts(iinj) % size_dmax - injection_pts(iinj) % size_dmin) * U1
         case(PTS_INJ_SIZEDIST_ROSINRAMMLER)
            !
            ! ROSIN_RAMMLER
            !
            U1      = random_generate_number(broadcast_seed=.true.)
            diame   = injection_pts(iinj) % size_dmin + injection_pts(iinj) % size_rr_dbar * (-1.0_rp*log(1.0_rp - U1*injection_pts(iinj) % size_rr_K))**(1.0_rp/injection_pts(iinj) % size_rr_n)
            diame   = min(injection_pts(iinj) % size_dmax,max(injection_pts(iinj) % size_dmin,diame))
         case(PTS_INJ_SIZEDIST_TABULATED)
            !
            ! Lookup framework based on tabulating the inverse of the CDF:
            ! D = f(U1)
            !
            control(1) = random_generate_number(broadcast_seed=.true.)
            ind        = 1_ip
            call fw_lookup( control, scale_control, lookup_fw( injection_pts(iinj) % kfl_size_lookupfw ), retva, ind )
            diame      = retva(1) 
       end select
    end function pts_injection_get_next_diame



    !-----------------------------------------------------------------------
    !> 
    !> @author  bsc21304
    !> @date    2019-11-08
    !> @brief   Get particle size distribution
    !> @details Based on the mass to be injected, get a particle size
    !>          dictribution, that fits a given cummulative density function 
    !> 
    !-----------------------------------------------------------------------

    subroutine pts_injection_rate_based_distribution(iinj, mass, denpa, len_list, d_list, mass_rem)
      use mod_physics,         only : physics_sphere_mass
      use def_master,          only : parii, npari, nparr, &
           parin, parre, mem_modul, &
           modul, ISLAVE, IMASTER

      use mod_memory,          only : memory_alloca
      use mod_memory,          only : memory_deallo
      use mod_memory,          only : memory_resize
      use mod_communications,  only : PAR_BROADCAST
      use mod_communications,  only : PAR_EXCHANGE
      implicit none

      integer(ip),       intent(in)  :: iinj
      real(rp),          intent(in)  :: mass
      real(rp),          intent(in)  :: denpa
      integer(ip),       intent(out) :: len_list
      real(rp), pointer, intent(inout) :: d_list(:)
      real(rp),          intent(out) :: mass_rem
      real(rp)                       :: diame, mass_i, mass_all
      integer(ip)                    :: ii, jj, itype, isafety
      logical(lg)                    :: list_full

      real(rp)                       :: d_avg, d_std
      real(rp), save                 :: d_next(pts_minj) = 0.0_rp
      character(20), parameter       :: vacal = 'pts_injection_rate_based_distribution'

      nullify(d_list)

      if( INOTSLAVE ) then

         itype = max(1_ip,injection_pts(iinj) % kfl_particle_type)
         !
         ! Initialize
         !
         len_list = 100_ip
         call memory_alloca(mem_modul(1:2, modul), 'D_LIST % A', vacal, d_list, len_list)

         ii             = 0_ip
         isafety        = 0_ip
         mass_all       = 0.0_rp
         mass_rem       = mass
         list_full= .false.
         if (mass == 0.0_rp) list_full = .true.

         fill_list: do while( .not. list_full)
            !
            ! Break if for some reason too many particles would be injected
            !
            isafety = isafety + 1_ip
            if ( isafety > 99999 ) then
               list_full = .true.
               call messages_live('NUMBER OF INJECTED PARTICLES REACHED LIMIT.', 'WARNING')
            endif
            !
            ! SELECT DIAMETER 
            !
            if (d_next(iinj) > zeror) then
               !
               ! Last selection was stored in d_next
               !
               diame        = d_next(iinj)
               d_next(iinj) = 0.0_rp
            else
               diame        = pts_injection_get_next_diame(iinj)
            endif

            mass_i  = parttyp(itype) % n_drop * physics_sphere_mass(diame, denpa)

            if (mass <= (mass_all + mass_i + 1.0e-20_rp) ) then
               !
               ! EXCEEDED MASS TO BE INJECTED
               !
               mass_rem      = mass - mass_all
               d_next(iinj)  = diame
               diame         = 0.0_rp 
               list_full     = .true.
            endif

            !
            ! Increment counters
            !
            if (diame > 0.0_rp) then
               mass_all = mass_all + mass_i
               ii = ii + 1_ip

               if (ii > size(d_list,KIND=ip)) then
                  !
                  ! Need to reallocate list, because it's too short
                  !
                  len_list = len_list + 100_ip
                  call memory_resize(mem_modul(1:2, modul), 'D_LIST % A', vacal, d_list, len_list)
               endif

               d_list(ii) = diame 
            endif
         end do fill_list

         !
         ! Update real length:
         !
         len_list = ii
         !
         ! Extract statistics 
         !
         d_avg = 0.0_rp
         d_std = 0.0_rp
         do jj = 1,len_list
            d_avg = d_avg + d_list(jj)
            d_std = d_std + d_list(jj)**2
         enddo
         d_avg = d_avg / max(real(len_list,rp),1.0_rp) 
         d_std = sqrt( max(0.0_rp,d_std / max(real(len_list-1,rp),1.0_rp) - d_avg**2) )

         if (len_list > 0_ip) call messages_live('INJECTOR '//trim(intost(iinj))//': '//trim(intost(len_list))// &
              & ' PARTICLES TO INJECT, AVG (STD) DIAMETER: '//trim(retost(d_avg*1.0e6_rp))//' um ('//             &
              & trim(retost(d_std*1.0e6_rp))//' um)', 'MODULE')
      end if
      call PAR_BROADCAST(len_list)

      if( len_list > 0 ) then
         if ( ISLAVE ) then 
            call memory_alloca(mem_modul(1:2, modul), 'D_LIST % A', vacal, d_list, max(1_ip,len_list))
         end if
         call PAR_BROADCAST(len_list,d_list)
      end if
      
    end subroutine pts_injection_rate_based_distribution

    subroutine pts_injection_number_based_distribution(iinj, len_list, d_list)
      use mod_physics,         only : physics_sphere_diameter
      use mod_physics,         only : physics_set_liquid_temperature
      use def_master,          only : parii, nparr, parre, mem_modul, &
           modul, ISLAVE, IMASTER

      use mod_memory,          only : memory_alloca
      use mod_memory,          only : memory_deallo
      use mod_memory,          only : memory_resize
      use mod_communications,  only : PAR_BROADCAST
      use mod_communications,  only : PAR_EXCHANGE
      implicit none

      integer(ip),       intent(in)  :: iinj
      integer(ip),       intent(in)  :: len_list
      real(rp), pointer, intent(inout) :: d_list(:)
      integer(ip)                    :: ii, jj, itype
      real(rp)                       :: tempe_loc 

      real(rp)                       :: d_avg, d_std
      character(20), parameter       :: vacal = 'pts_injection_number_based_distribution'

      nullify(d_list)

      call memory_alloca(mem_modul(1:2, modul), 'D_LIST % A', vacal, d_list, len_list)

      if (injection_pts(iinj) % kfl_size_dist == PTS_INJ_SIZEDIST_CONST) then
         !
         ! Constant
         !
         d_list = injection_pts(iinj) % size_diame
      elseif (injection_pts(iinj) % kfl_size_dist == PTS_INJ_SIZEDIST_FILED) then
         !
         ! Listed diameters 
         !
         d_list(1:len_list) = injection_pts(iinj) % size_list(1:len_list)
      elseif (injection_pts(iinj) % kfl_size_dist == PTS_INJ_SIZEDIST_FILEM) then
         !
         ! Listed masses
         !
         itype = injection_pts(iinj) % kfl_particle_type
         if ( injection_pts(iinj) % kfl_tempe == PTS_INJ_TEMPE_CONST ) then
            !
            ! Constant temperature
            !
            tempe_loc = injection_pts(iinj) % tempe
            call physics_set_liquid_temperature( parttyp(itype) % liq , tempe_loc )
         endif
         do jj = 1,len_list
            !
            ! Needs to know temperature for density  
            !
            if ( injection_pts(iinj) % kfl_tempe == PTS_INJ_TEMPE_FILE ) then
               !
               ! Temperature from input file
               !
               tempe_loc = injection_pts(iinj) % tempe_list(jj)
               call physics_set_liquid_temperature( parttyp(itype) % liq , tempe_loc )
            end if
            d_list(jj) = physics_sphere_diameter(injection_pts(iinj) % size_list(jj),&
                 & parttyp(itype) % liq % rho)
         enddo
      else
         if( INOTSLAVE ) then
            fill_list: do ii = 1,len_list
               !
               ! SELECT DIAMETER 
               !
               d_list(ii)   = pts_injection_get_next_diame(iinj)
            end do fill_list

            !
            ! Extract statistics 
            !
            d_avg = 0.0_rp
            d_std = 0.0_rp
            do jj = 1,len_list
               d_avg = d_avg + d_list(jj)
               d_std = d_std + d_list(jj)**2
            enddo
            d_avg = d_avg / max(real(len_list,rp),1.0_rp) 
            d_std = sqrt( d_std / max(real(len_list-1,rp),1.0_rp) - d_avg**2)

            if (len_list > 0_ip) call messages_live('INJECTOR '//trim(intost(iinj))//': '//trim(intost(len_list))// &
                 & ' PARTICLES TO INJECT, AVG (STD) DIAMETER: '//trim(retost(d_avg*1.0e6_rp))//' um ('//             &
                 & trim(retost(d_std*1.0e6_rp))//' um)', 'MODULE')
         end if
         call PAR_BROADCAST(d_list)

      endif

    end subroutine pts_injection_number_based_distribution



    !-----------------------------------------------------------------------
    !> 
    !> @author  houzeaux
    !> @date    2020-07-06
    !> @brief   Find host elements
    !> @details Find the host elements for a list of particles
    !> 
    !-----------------------------------------------------------------------
    
    subroutine pts_injection_ini_loc_and_cond(&
         nlagr_pos,nlagr_new,nlagr_inj,particle_position,particle_injector,&
         particle_type,particle_diameter)

      implicit none
      integer(ip),                     intent(in)    :: nlagr_pos                   !< Number of injected particles
      integer(ip),                     intent(inout) :: nlagr_new                   !< Number of new particles owned by myself
      integer(ip),                     intent(inout) :: nlagr_inj                   !< Number of injected particles
      real(rp),               pointer, intent(in)    :: particle_position(:,:)      !< (ndime,nlagr_pos)
      integer(ip),            pointer, intent(in)    :: particle_injector(:)        !< (nlagr_pos)
      integer(ip),            pointer, intent(in)    :: particle_type(:)            !< (nlagr_pos)
      real(rp),     optional, pointer, intent(in)    :: particle_diameter(:)        !< (nlagr_pos)
      real(rp),               pointer                :: host_shapf(:,:)          
      integer(ip),            pointer                :: host_element(:)
      integer(ip),            pointer                :: particle_place(:) 
  
      character(20), parameter :: vacal = 'pts_injection_ini_loc_and_cond'
  

      nlagr_new = 0
      nullify(host_shapf)
      nullify(host_element)
      nullify(particle_place)
      !
      ! Alllocate
      !
      call memory_alloca(mem_modul(1:2,modul),'HOST_ELEMENT'  ,vacal,host_element  ,nlagr_pos)
      call memory_alloca(mem_modul(1:2,modul),'HOST_SHAPF'    ,vacal,host_shapf    ,mnode,nlagr_pos)
      call memory_alloca(mem_modul(1:2,modul),'PARTICLE_PLACE',vacal,particle_place,max(1_ip,nlagr_inj))

      if( INOTEMPTY .and. nlagr_pos > 0 ) then
         !
         ! Find host elements
         !
         call pts_host_element(&
              nlagr_pos,nlagr_new,particle_position,&
              host_shapf,host_element)
         !
         ! Initialize particles
         !
         call pts_initial_condition(&
              nlagr_pos,nlagr_new,particle_position,particle_injector,&
              particle_place,particle_type,host_shapf,host_element,&
              particle_diameter)
      end if
      !
      ! Parallelization, treat particles with multiple hosts
      !
      call pts_host_element_parall(&
           nlagr_new,nlagr_inj,particle_place)
      !
      ! Deallocate
      !
      call memory_deallo(mem_modul(1:2,modul),'HOST_ELEMENT'  ,vacal,host_element)
      call memory_deallo(mem_modul(1:2,modul),'HOST_SHAPF'    ,vacal,host_shapf)
      call memory_deallo(mem_modul(1:2,modul),'PARTICLE_PLACE',vacal,particle_place)

    end subroutine pts_injection_ini_loc_and_cond
      
    !-----------------------------------------------------------------------
    !> 
    !> @author  houzeaux
    !> @date    2020-07-06
    !> @brief   Find host elements
    !> @details Find the host elements for a list of particles
    !> 
    !-----------------------------------------------------------------------
    
    subroutine pts_initial_condition(&
         nlagr_pos,nlagr_new,particle_position,particle_injector,&
         particle_place,particle_type,host_shapf,host_element,&
         particle_diameter)

      use def_parame,                  only : xmaxint4,pi
      use def_master,                  only : zeror
      use def_master,                  only : kfl_paral
      use def_master,                  only : dtime
      use def_master,                  only : cutim
      use def_master,                  only : advec
      use def_domain,                  only : ndime
      use def_domain,                  only : mnode
      use def_domain,                  only : lnods
      use def_domain,                  only : lnnod
      use def_kermod,                  only : grnor,gravi
      use mod_ker_proper,              only : ker_proper
      use mod_physics,                 only : physics_sphere_mass
      use mod_pts_thermodynamic,       only : pts_thermodynamic_properties
      use mod_ker_space_time_function, only : ker_space_time_function
      use mod_random,                  only : random_generate_number
      implicit none

      integer(ip),                      intent(in)    :: nlagr_pos                   !< Number of injected particles
      integer(ip),                      intent(in)    :: nlagr_new                   !< Number of new particles owned by myself
      real(rp),               pointer,  intent(in)    :: particle_position(:,:)      !< (ndime,nlagr_pos)
      integer(ip),            pointer,  intent(in)    :: particle_injector(:)        !< (nlagr_pos)
      integer(ip),            pointer,  intent(inout) :: particle_place(:)           !< (nlagr_pos)
      integer(ip),            pointer,  intent(in)    :: particle_type(:)            !< (nlagr_pos)
      real(rp),               pointer,  intent(in)    :: host_shapf(:,:)             !< Shape functions
      integer(ip),            pointer,  intent(in)    :: host_element(:)             !< Host elements
      real(rp),               pointer,  intent(in)    :: particle_diameter(:)        !< (nlagr_pos)

      integer(ip)                                     :: new_size
      integer(ip)                                     :: ielem,inode,ilagr,klagr_last
      integer(ip)                                     :: ipoin,iinj,itype
      integer(ip)                                     :: pnode,dumm0,klagr,nlagr_free
      real(rp)                                        :: xx(3),nn(3),xa(3),xc(3)
      real(rp)                                        :: v_tmp(3),v_scaling
      real(rp)                                        :: rvec(3),ang_min,ang_max,vel_mag,weight
      real(rp)                                        :: dummr(3),shapf(mnode),rad
      real(rp)                                        :: grafo,buofo,denpa,dista
      real(rp)                                        :: confl,sphfl,denfl,visfl,Dvg_m,Yv_surf,Yv_fluid_k
      real(rp)                                        :: Therm_fluid_k,T_fluid_k,conce_seen(nclas_pts)
      real(rp)                                        :: Pr_m, Sc_m, LK_m, xvap, w_nonf_seen
      real(rp)                                        :: tau_p,diame
      real(rp)                                        :: sigma,angle
      real(rp)                                        :: U1,U2,U3 

      xa  = 0.0_rp
      xc  = 0.0_rp
      xx  = 0.0_rp
      
      !----------------------------------------------------------------------
      !
      ! Reallocate LAGRTYP if necessary
      !
      !----------------------------------------------------------------------
      !
      ! Counter number of free places in LAGRTYP
      !
      nlagr_free = 0
      do klagr = 1,mlagr
         if( lagrtyp(klagr) % kfl_exist == 0 ) then
            nlagr_free = nlagr_free + 1 
         end if
      end do
      !
      ! Reallocate LAGRTYP if necessary
      ! I currently have MLAGR places
      ! NLAGR_FREE are available
      ! I need NLAGR_NEW new positions
      !
      if( nlagr_new > nlagr_free ) then
         dista = 1.2_rp*real(mlagr+nlagr_new-nlagr_free,rp)
         if( ip == 4 .and. dista > xmaxint4 ) then
            call runend('PTS_HOST_ELEMENT: GOT TO 8 BYTES INTEGERS')
         else
            new_size = int(dista,ip)
            call pts_reallocate(new_size)
         end if
      end if

      !----------------------------------------------------------------------
      !
      ! Fill in LAGRTYP with new particles with host element
      !
      !----------------------------------------------------------------------

      klagr_last = 0

      !-$OMP PARALLEL DO SCHEDULE (DYNAMIC,500)                                          & 
      !-$OMP DEFAULT       (NONE)                                                        &    
      !-$OMP FIRSTPRIVATE  (klagr_last)                                                  &
      !-$OMP PRIVATE       (ilagr,ielem,xx,shapf,dumm0,klagr,itype,inode,ipoin,iinj,     &
      !-$OMP                nn,sigma,xc,pnode,dummr,denfl,denpa,grafo,buofo,v_scaling,   &   
      !-$OMP                angle)                                                       &   
      !-$OMP SHARED        (nlagr_pos,host_element,number_types_pts,particle_position,   &
      !-$OMP                host_shapf,mlagr,lagrtyp,particle_place,parttyp,gravi,       &
      !-$OMP                lnnod,lnods,injection_pts,particle_injector,ndime,grnor,     &
      !-$OMP                kfl_exacs_pts,mnode,nlacc_pts,kfl_adapt_pts,dtime,dtmin_pts, &
      !-$OMP                advec)                                                       &
      do ilagr = 1,nlagr_pos

         ielem = host_element(ilagr)

         if( ielem > 0 ) then

            xx(1:ndime)    = particle_position(1:ndime,ilagr)
            itype          = particle_type(ilagr)
            iinj           = particle_injector(ilagr)
            shapf(1:mnode) = host_shapf(1:mnode,ilagr)
            pnode          = lnnod(ielem)
            !
            ! Look for new available particle position: 
            !
            ! OPTIMIZE: klagr could be 0 at the begining of injection
            ! so that we do not start from 0 each time!!!!!
            !
            klagr = klagr_last
            loop_klagr: do while( klagr < mlagr )
               klagr = klagr + 1
               if( lagrtyp(klagr) % kfl_exist == 0 ) then
                  exit loop_klagr
               end if
            end do loop_klagr
            klagr_last            = klagr
            particle_place(ilagr) = klagr                         ! To compare with my neighbor
            !
            ! Initialize particle
            !
            call lagrtyp(klagr) % init()                          ! Initial value
            lagrtyp(klagr) % ilagr     =  ilagr + nlacc_pts       ! Particle absolute ID
            lagrtyp(klagr) % itype     =  itype                   ! Type
            lagrtyp(klagr) % ielem     =  ielem                   ! Host element
            lagrtyp(klagr) % ittim     =  0_ip                    ! Time step
            lagrtyp(klagr) % kfl_exist = -5_ip                    ! Just been injected
            lagrtyp(klagr) % coord     =  xx                      ! Initial coordinates
            lagrtyp(klagr) % t         =  cutim - dtime           ! Current time
            lagrtyp(klagr) % t_inject  =  cutim - dtime           ! Injection time
            lagrtyp(klagr) % mpi_rank  =  kfl_paral               ! My rank                

            !
            ! Initialize diameter for all models
            ! It will be preserved across MPI migrations 
            ! if it is used in the calculations or requested as an output variable
            !
            diame                      =  particle_diameter(ilagr)! Particle diameter
            lagrtyp(klagr) % diam_0    =  diame                   
            lagrtyp(klagr) % diam_k    =  diame



            !----------------------------------------------------------------
            !
            ! Thermodynamic model
            !
            !---------------------------------------------------------------- 
            if( parttyp(itype) % kfl_therm /= 0 ) then
               !
               ! Temperature
               !
               if ( injection_pts(iinj) % kfl_tempe == PTS_INJ_TEMPE_CONST ) then
                  !
                  ! Constant temperature
                  !
                  lagrtyp(klagr) % tempe_k = injection_pts(iinj) % tempe
               elseif ( injection_pts(iinj) % kfl_tempe == PTS_INJ_TEMPE_FILE ) then
                  !
                  ! Temperature from input file
                  !
                  lagrtyp(klagr) % tempe_k = injection_pts(iinj) % tempe_list(ilagr-sum(injection_pts(1:iinj-1) % num_part))
               end if
               !
               ! Set liquid state to get density
               ! And calculate fluid density and viscosity according to the 1/3 law
               !
               call pts_thermodynamic_properties(itype, lagrtyp(klagr) % tempe_k, ielem, pnode, lnods(1:pnode,ielem), shapf, &
                    confl, sphfl, denfl, visfl, Dvg_m, Pr_m, Sc_m, LK_m, Yv_surf, Yv_fluid_k, &
                    Therm_fluid_k, T_fluid_k, xvap, w_nonf_seen, conce_seen)
               denpa     = parttyp(itype) % liq % rho               ! Particle density 
               !
               ! Previous solutions
               !
               lagrtyp(klagr) % tempe_km1 = lagrtyp(klagr) % tempe_k
               lagrtyp(klagr) % tempe_km2 = lagrtyp(klagr) % tempe_k
            else
               !
               ! Constant temperature model
               !
               call ker_proper('VISCO','IGAUS',dumm0,ielem,dummr,pnode,1_ip,shapf)   ! Fluid viscosity
               visfl = max(zeror,dummr(1))                                           ! Fluid viscosity
               call ker_proper('DENSI','IGAUS',dumm0,ielem,dummr,pnode,1_ip,shapf)   ! Fluid density
               denfl = dummr(1)           
               denpa = parttyp(itype) % denpa                                        ! Particle density
            endif

            !
            ! Initialize mass for all models
            ! It will be preserved across MPI migrations 
            ! if it is used in the calculations or requested as an output variable
            !
            lagrtyp(klagr) % mass_k = physics_sphere_mass(diame,denpa)
            lagrtyp(klagr) % mass_km1  = lagrtyp(klagr) % mass_k
            lagrtyp(klagr) % mass_km2  = lagrtyp(klagr) % mass_k
            lagrtyp(klagr) % mass_0    = lagrtyp(klagr) % mass_k


            !
            ! Initial local & effective Stokes number
            !
            tau_p  = (denpa * diame * diame)/(18.0_rp*visfl)
            call pts_stk_local(ielem,tau_p, lagrtyp(klagr) % t_inject, lagrtyp(klagr) % Stk(1), lagrtyp(klagr) % Stk(2))    
            !
            ! Initial time step
            !
            if( parttyp(itype) % kfl_modla == 2 ) then
               if(      parttyp(itype) % kfl_tstep == 0 ) then
                  lagrtyp(klagr) % dt_k = dtime
               else if( parttyp(itype) % kfl_tstep <  0 ) then
                  lagrtyp(klagr) % dt_k = parttyp(itype) % dtime
               else
                  lagrtyp(klagr) % dt_k = dtmin_pts
               end if
            else
               lagrtyp(klagr) % dt_k =  dtime
            end if
            lagrtyp(klagr) % dt_km1 = lagrtyp(klagr) % dt_k
            lagrtyp(klagr) % dt_km2 = lagrtyp(klagr) % dt_k
            !
            ! Fluid velocity
            !
            lagrtyp(klagr) % v_fluid_k =  0.0_rp
            if( associated(advec) ) then
               do inode = 1,lnnod(ielem)
                  ipoin = lnods(inode,ielem)  
                  lagrtyp(klagr) % v_fluid_k(1:ndime) = lagrtyp(klagr) % v_fluid_k(1:ndime) &
                       + shapf(inode) * advec(1:ndime,ipoin,3)
               end do
            end if
            lagrtyp(klagr) % v_fluid_km1 = lagrtyp(klagr) % v_fluid_k
            lagrtyp(klagr) % v_fluid_km2 = lagrtyp(klagr) % v_fluid_k

            !----------------------------------------------------------------
            !
            ! Initial particle velocity
            !
            !----------------------------------------------------------------        
            
            v_scaling = 1.0_rp
            if ( injection_pts(iinj) % kfl_velfun /= 0 ) then
               call ker_space_time_function(&
                    injection_pts(iinj) % kfl_velfun,xx(1),xx(2),xx(ndime),cutim,v_scaling)
            endif
            !
            ! Fluctuations
            !
            if ( injection_pts(iinj) % kfl_fluct_veloc /= PTS_INJ_FLUCT_VELOC_CONST ) then
               if ( injection_pts(iinj) % kfl_fluct_veloc == PTS_INJ_FLUCT_VELOC_UNIFORM ) then
                  !
                  ! Uniform distribution with given standard deviation
                  ! var = 1/12 * (max-min)^2             delta
                  ! std = 1/(2*sqrt(3)) * 2 * delta     |------|------|
                  ! std = 1/sqrt(3) * delta                  mean
                  ! delta = sqrt(3) * std
                  !
                  U1      = random_generate_number(broadcast_seed=.true.)
                  v_scaling = v_scaling * ( 1.0_rp + 2.0_rp * ( U1-0.5_rp ) * sqrt(3.0_rp) * injection_pts(iinj) % fluct_vel_std ) 
               elseif ( injection_pts(iinj) % kfl_fluct_veloc == PTS_INJ_FLUCT_VELOC_NORMAL ) then
                  !
                  ! Normal distribution with given standard deviation
                  ! Following the Box-Muller method
                  !
                  U1      = random_generate_number(broadcast_seed=.true.)
                  U2      = random_generate_number(broadcast_seed=.true.)
                  !
                  ! This will be a normally disributed variable:
                  !
                  U3      = sqrt(-2.0_rp * log(U1)) * cos(2.0_rp * pi * U2)
                  v_scaling = v_scaling * ( 1.0_rp + U3 * injection_pts(iinj) % fluct_vel_std ) 
               endif
            endif


            !
            ! VELOCITY STRATEGY:
            !
            if (     injection_pts(iinj) % kfl_veloc == PTS_INJ_VELOC_ZERO ) then
               !
               ! Zero velocity
               !
               lagrtyp(klagr) % veloc(1:ndime) = 0.0_rp 

            else if ( injection_pts(iinj) % kfl_veloc == PTS_INJ_VELOC_FLUID ) then
               !
               ! Fluid velocity
               !
               lagrtyp(klagr) % veloc(1:ndime) = v_scaling * lagrtyp(klagr) % v_fluid_k(1:ndime)

            else if ( injection_pts(iinj) % kfl_veloc == PTS_INJ_VELOC_NORMAL ) then
               !
               ! With respect to normal
               !
               lagrtyp(klagr) % veloc(1:ndime) = v_scaling * injection_pts(iinj) % vel_mag * injection_pts(iinj) % geo_normal(1:ndime)

            else if ( injection_pts(iinj) % kfl_veloc  == PTS_INJ_VELOC_GAUSSIAN ) then
               !
               ! Gaussian velocity injector
               ! f(x) = exp[(x-x0)^2/(2*sigma^2)]
               !
               xc    = injection_pts(iinj) % geo_coord 
               nn    = injection_pts(iinj) % geo_normal
               sigma = injection_pts(iinj) % vel_sigma
               lagrtyp(klagr) % veloc(1:ndime) = v_scaling * injection_pts(iinj) % vel_mag * (-1.0_rp) * nn(1:ndime) * exp(-(dot_product(xx-xc,xx-xc)) / (2.0_rp*sigma**2)) 

            else if ( injection_pts(iinj) % kfl_veloc == PTS_INJ_VELOC_CONIC ) then
               !
               ! Conic velocity 
               !
               angle    = injection_pts(iinj) % vel_angle*(pi/180.0_rp)
               rad      = injection_pts(iinj) % geo_rad      
               nn       = injection_pts(iinj) % geo_normal

               if( ndime == 3 ) then    
                  
                  xa       = injection_pts(iinj) % geo_coord + rad*cos(angle)/sin(angle)*nn
                  v_tmp    = xa - xx
                  v_tmp    = v_tmp / (sqrt(dot_product(v_tmp,v_tmp))+zeror)

                  lagrtyp(klagr) % veloc(1:3) = - v_scaling * injection_pts(iinj) % vel_mag * v_tmp(1:3)
               end if

            else if ( injection_pts(iinj) % kfl_veloc == PTS_INJ_VELOC_SPRAY .or. &
                &     injection_pts(iinj) % kfl_veloc == PTS_INJ_VELOC_AXSPRAY ) then
               !
               ! Spray velocity
               !
               angle     = injection_pts(iinj) % vel_angle*(pi/180.0_rp)  ! angle at outer edge of surface
               rad       = injection_pts(iinj) % geo_rad                  ! radius at outer edge
               nn        = injection_pts(iinj) % geo_normal               ! normal of injection surface

               !
               ! Calculate origin point behind injector plane, 
               ! such that the specified angle is recovered on the edge of the plane
               !                                                       
               !                       | alpha                             
               !                       |     /                           
               !                       | r  /                           
               !                   ----|----                           
               !                    \  |  /   ^                        
               !                     \ | /    | r * ctg(alpha)                       
               !                      \|/     |                        
               !                       o      v                        
               !
               xa        = injection_pts(iinj) % geo_coord + rad*cos(angle)/sin(angle)*nn
               !
               ! Create vector pointing from origin point to the injection point
               !
               v_tmp     = xx - xa

               if ( injection_pts(iinj) % kfl_veloc == PTS_INJ_VELOC_SPRAY ) then
                  !
                  ! MAGNITUDE is specified:
                  !
                  v_tmp                       = v_tmp / (sqrt(dot_product(v_tmp,v_tmp))+zeror) ! normalize vector
                  lagrtyp(klagr) % veloc(1:3) = v_scaling * injection_pts(iinj) % vel_mag * v_tmp(1:3)

               else if ( injection_pts(iinj) % kfl_veloc == PTS_INJ_VELOC_AXSPRAY ) then 
                  !
                  ! AXIAL component is specified
                  ! divide v_tmp r ctg(alpha) so it's unity in the center
                  ! and it increases towards the edges such that the axial 
                  ! component is constant
                  !
                  v_tmp                       = v_tmp * sin(angle) / (rad * cos(angle)+zeror) 
                  lagrtyp(klagr) % veloc(1:3) = v_scaling * injection_pts(iinj) % vel_ax * v_tmp(1:3)
                  
               endif

            else if ( injection_pts(iinj) % kfl_veloc == PTS_INJ_VELOC_SIZE_DEP_SPRAY  ) then
                !
                ! Hollow cone spray with defined spray angles depending on the droplet size
                !
                nn =   injection_pts(iinj) % geo_normal     ! normal of injection surface
                rvec = xx - injection_pts(iinj) % geo_coord ! radial vector of injection location 
                rad = sqrt(dot_product(rvec,rvec))
                rvec = rvec / rad                           ! normalized radial vector

                !
                ! Decide regime
                !
                if (diame <= injection_pts(iinj) % vel_diam_s) then
                    vel_mag = injection_pts(iinj) % vel_mag_s     
                    ang_max = injection_pts(iinj) % vel_ang_max_s 
                    ang_min = injection_pts(iinj) % vel_ang_min_s 
                elseif (diame >= injection_pts(iinj) % vel_diam_L) then
                    vel_mag = injection_pts(iinj) % vel_mag 
                    ang_max = injection_pts(iinj) % vel_ang_max_L 
                    ang_min = injection_pts(iinj) % vel_ang_min_L 
                else
                    weight = (diame - injection_pts(iinj) % vel_diam_s) / (injection_pts(iinj) % vel_diam_L - injection_pts(iinj) % vel_diam_s)
                    vel_mag = injection_pts(iinj) % vel_mag       * weight + injection_pts(iinj) % vel_mag_s     * (1.0_rp - weight) 
                    ang_max = injection_pts(iinj) % vel_ang_max_L * weight + injection_pts(iinj) % vel_ang_max_s * (1.0_rp - weight)
                    ang_min = injection_pts(iinj) % vel_ang_min_L * weight + injection_pts(iinj) % vel_ang_min_s * (1.0_rp - weight) 
                endif
                
                !
                ! Determine angle
                !
                U1 = random_generate_number(broadcast_seed=.true.)
                angle = (ang_min + U1 * (ang_max - ang_min))*(pi/180.0_rp)
                v_tmp = -1.0_rp * cos(angle) * nn + sin(angle) * rvec

                lagrtyp(klagr) % veloc(1:3) = v_scaling * vel_mag * v_tmp(1:3)

            else if ( injection_pts(iinj) % kfl_veloc == PTS_INJ_VELOC_CONST ) then
               !
               ! Constant velocity 
               !
               lagrtyp(klagr) % veloc(1:ndime) = v_scaling * injection_pts(iinj) % vel_vec(1:ndime)        

            else if ( injection_pts(iinj) % kfl_veloc == PTS_INJ_VELOC_FILE ) then
               !
               ! Velocity from input file
               !
               lagrtyp(klagr) % veloc(1:ndime) = v_scaling * injection_pts(iinj) % vel_list(1:ndime,ilagr-sum(injection_pts(1:iinj-1) % num_part))        
            end if

            !----------------------------------------------------------------
            !
            ! Initial acceleration
            !
            !----------------------------------------------------------------

            if( parttyp(itype) % kfl_modla == 2 ) then
               grafo = real( parttyp(itype) % kfl_grafo, rp )                        ! Gravity  force = 1.0
               buofo = real( parttyp(itype) % kfl_buofo, rp )                        ! Buoyancy force = 1.0  
               lagrtyp(klagr) % accel(1:ndime) = grnor * gravi(1:ndime) * ( grafo - buofo * denfl / denpa )
            end if

            !
            ! Exact solution
            !
            if( kfl_exacs_pts /= 0 ) then
               call pts_exacso(2_ip,0.0_rp,lagrtyp(klagr) % accel,lagrtyp(klagr) % veloc,lagrtyp(klagr) % coord)
            end if
         end if
      end do
      !-$OMP END PARALLEL DO

    end subroutine pts_initial_condition



    !-----------------------------------------------------------------------
    !> 
    !> @author  bsc21240
    !> @date    2018-03-08
    !> @brief   Inject particles
    !> @details Inject particles according to the different injectors
    !> 
    !-----------------------------------------------------------------------

    subroutine pts_injection_injectors()
      use mod_maths,                   only : maths_vector_from_new_basis
      use mod_maths,                   only : maths_local_orthonormal_basis
      use mod_ker_timeline,            only : ker_timeline
      use mod_random,                  only : random_generate_number
      use mod_random,                  only : random_initialization
      use mod_random,                  only : random_end
      use mod_physics,                 only : physics_set_liquid_temperature
      use mod_ker_space_time_function, only : ker_space_time_function

      implicit none
      integer(ip)          :: ii, jj, kk, iinj
      integer(ip)          :: i1, i2, i3, j1, j2, j3
      real(rp)             :: zfact, tfact, rfact
      real(rp)             :: xx(3), xmini, xmaxi, ymini, ymaxi, zmini, zmaxi, nn(3)
      real(rp)             :: x1, y1, z1, x2, y2, z2, x3, y3, z3, p1(3), p2(3), p3(3)
      real(rp)             :: ux, uy, uz, wx, wy, wz, d1, d2, d3, U1, U2, U3
      real(rp)             :: t_aux, t, z, r, radius, rad, rad_min, rr, hh
      real(rp)             :: injected_quantity, denpa 
      integer(ip)          :: nside, nlagr_new, nlagr_pos, nlagr_pos2, itype,dummi
      integer(ip)          :: number_injected_pts, num_types_injected, num_particles_per_injection
      
      real(rp),   external :: funcre
      real(rp),    pointer :: particle_position(:,:)
      integer(ip), pointer :: particle_injector(:)
      integer(ip), pointer :: particle_type(:)
      integer(ip), pointer :: particle_place(:)

      type(r1p),   pointer :: particle_diameter_inj(:)
      real(rp),    pointer :: particle_diameter(:)
      integer(ip), pointer :: n_inj(:)

      real(rp),    save    :: quantity_remaining(pts_minj) = 0.0_rp
      real(rp)             :: cutla_loc_pts(pts_minj)         
      real(rp)             :: flow_rate, flow_rate_scale 
      !
      ! Initialize random generator. If kfl_randseed == 0, then a fixed
      ! seed is set to enable reproducibility over several runs
      !
      call random_initialization(broadcast_seed=.true.)
      !
      ! By default NO particle is injected
      !
      kfl_injec           = 0
      nlagr_pos           = 0
      nlagr_pos2          = 0
      number_injected_pts = 0
      cutla_loc_pts       = 0.0_rp
      !
      ! Decide if any particles should be injected:
      !
      do iinj = 1, kfl_imax_injector_pts
         
         if ( cutim >= injection_pts(iinj) % time_initial .and. &
           &  cutim <  injection_pts(iinj) % time_final) then
            !
            ! Accumulate time since last injection 
            !
            injection_pts(iinj) % time_cumulative = injection_pts(iinj) % time_cumulative + dtime

            !
            ! If accumulated time surpasses period time, allow injection
            !
            if (injection_pts(iinj) % time_cumulative >= injection_pts(iinj) % time_period - zeror) then
               if (injection_pts(iinj) % time_cumulative > 1.0e11_rp) then
                  cutla_loc_pts(iinj) = dtime
               else
                  cutla_loc_pts(iinj) = injection_pts(iinj) % time_cumulative
               endif
               kfl_injec = 1
               injection_pts(iinj) % time_cumulative = 0.0_rp
            end if
         end if
      enddo
      if (cutim >= tfila_pts) then
         kfl_injec = 0
      end if

      !
      ! Break execution if none of the injectors are active.
      !
      if (kfl_injec == 0) return


      call ker_timeline('INI_ASSEMBLY')
      !
      ! Nullify pointers
      !
      nullify(particle_diameter_inj)
      nullify(n_inj)
      call memory_alloca(mem_modul(1:2,modul),'D_LIST' ,'pts_injection_injectors',particle_diameter_inj ,kfl_imax_injector_pts)
      call memory_alloca(mem_modul(1:2,modul),'N_INJ'  ,'pts_injection_injectors',n_inj ,kfl_imax_injector_pts)
      !
      ! Number of injected particles when injection
      !
      do iinj = 1, kfl_imax_injector_pts
         n_inj(iinj) = 0
         !
         ! Only inject if given injector is active
         !
         if (cutla_loc_pts(iinj) > 0.0_rp) then
            num_types_injected = 0
            if (injection_pts(iinj) % kfl_geometry /= 0) then
               if (injection_pts(iinj) % kfl_particle_type == 0) then
                  num_types_injected = number_types_pts
               else
                  num_types_injected = 1
               end if
            end if

            !
            ! Evaluate space-time functions if necessary
            ! 
            flow_rate_scale = 1.0_rp
            if ( injection_pts(iinj) % kfl_flowfun /= 0 ) then
               !
               ! Flow rate imposed by space & time function 
               !
               if ( injection_pts(iinj) % kfl_flowfun > 0 ) then
                  xx = 0.0_rp    !! Time dependent function only
                  call ker_space_time_function(&
                             injection_pts(iinj) % kfl_flowfun,xx(1),xx(2),xx(ndime),cutim,flow_rate_scale)
               end if

               !
               ! Flow rate imposed by time function 
               !
               if ( injection_pts(iinj) % kfl_flowfun < 0 ) then
                  flow_rate_scale =   funcre(                                                      &
                         time_function(-1_ip * injection_pts(iinj) % kfl_flowfun) % parameters,    &
                         time_function(-1_ip * injection_pts(iinj) % kfl_flowfun) % npara,         &
                         time_function(-1_ip * injection_pts(iinj) % kfl_flowfun) % kfl_type,      &
                         cutim)
               end if
            endif


            !
            ! NUMBER OF INJECTED PARTICLES N_INJ(IINJ) of type IINJ
            !
            select case (injection_pts(iinj) % kfl_flow)
               
            case(PTS_INJ_FLOW_NMODIFIED)
               !
               ! The number of injected particles has to be modified to accommodate
               ! some specific deterministic arrangements.
               ! I do NOT like this practice, but I will implement it for compatbility.
               !
               nside = injection_pts(iinj) % num_part
               select case ( injection_pts(iinj) % kfl_geometry )
               case ( 1_ip) ; n_inj(iinj) =  nside * nside       ! SQUARE 
               case ( 2_ip) ; n_inj(iinj) =  nside ** ndime      ! SPHERE
               case ( 3_ip) ; n_inj(iinj) =  2 * nside ** ndime  ! SEMI-SPHERE
               case ( 4_ip) ; n_inj(iinj) =  nside * nside       ! CIRCLE
               case ( 5_ip) ; n_inj(iinj) =  nside * nside       ! RECTANGLE
               case ( 6_ip) ; n_inj(iinj) =  nside               ! POINTWISE
               case ( 7_ip) ; n_inj(iinj) =  nside ** 3          ! CONE 
               case ( 8_ip) ; n_inj(iinj) =  nside               ! RANDOM IN CUBE/SQUARE
               case ( 9_ip) ; n_inj(iinj) =  nside               ! SEGMENT
               case (10_ip) ; n_inj(iinj) =  nside               ! FILE
               case (11_ip) ; n_inj(iinj) =  nside * nside       ! ANNULUS
               end select

            case(PTS_INJ_FLOW_NASIS)
               !
               ! Number of injected particles is defined in geometry 
               !
               n_inj(iinj) = injection_pts(iinj) % num_part

            case(PTS_INJ_FLOW_MASSFLOW)
               !
               ! PARTICLE DENSITY
               !
               itype = max(1_ip,injection_pts(iinj) % kfl_particle_type)
               if( parttyp(itype) % kfl_therm /= 0 ) then
                  call physics_set_liquid_temperature( parttyp(itype) % liq , injection_pts(iinj) % tempe )
                  denpa = parttyp(itype) % liq % rho
               else
                  denpa = parttyp(itype) % denpa  
               endif
               flow_rate = injection_pts(iinj) % flow_rate * flow_rate_scale
                
            case(PTS_INJ_FLOW_VOLUMEFLOW)
               denpa = 1.0_rp
               flow_rate = injection_pts(iinj) % flow_rate * flow_rate_scale
            end select

            !
            ! DROPLET SIZES, (and number of particles for flow-rate based methods)
            !
            if ( (injection_pts(iinj) % kfl_flow == PTS_INJ_FLOW_MASSFLOW) .or. &
                &(injection_pts(iinj) % kfl_flow == PTS_INJ_FLOW_VOLUMEFLOW) ) then
               !
               ! Rate based injection methods
               !
               injected_quantity        = cutla_loc_pts(iinj) * flow_rate + quantity_remaining(iinj)
               quantity_remaining(iinj) = 0.0_rp

               call pts_injection_rate_based_distribution(iinj, injected_quantity, denpa, n_inj(iinj), &
                   &                           particle_diameter_inj(iinj)%a, quantity_remaining(iinj) )
            else
               !
               ! Number based injection 
               !
               call pts_injection_number_based_distribution(iinj,n_inj(iinj)* num_types_injected,particle_diameter_inj(iinj)%a)

               !
               ! Overwrite diameter for legacy method of injecting multiple types with one injector 
               !
               if ( num_types_injected > 1 .and. INOTMASTER ) then
                  do jj = 1, n_inj(iinj) * num_types_injected
                     itype = mod(jj,num_types_injected) + 1_ip
                     particle_diameter_inj(iinj) % a(jj) = parttyp(itype) % diame
                  enddo
               endif
            endif

            !
            ! Total number of injected particles
            !
            number_injected_pts = number_injected_pts + n_inj(iinj) * num_types_injected
         endif
      end do

      !
      ! Nullify pointers
      !
      nullify(particle_position)
      nullify(particle_injector)
      nullify(particle_type)
      nullify(particle_place)
      nullify(particle_diameter)

      if (INOTMASTER) then
         ! 
         ! Every subdomain does the same
         ! Later it is decided to which subdomain do they belong.
         !
         call memory_alloca(mem_modul(1:2, modul), 'PARTICLE_POSITION', 'pts_injection_injectors', particle_position, ndime, max(1_ip,number_injected_pts))
         call memory_alloca(mem_modul(1:2, modul), 'PARTICLE_INJECTOR', 'pts_injection_injectors', particle_injector, max(1_ip,number_injected_pts))
         call memory_alloca(mem_modul(1:2, modul), 'PARTICLE_TYPE',     'pts_injection_injectors', particle_type,     max(1_ip,number_injected_pts))
         call memory_alloca(mem_modul(1:2, modul), 'PARTICLE_PLACE'   , 'pts_injection_injectors', particle_place,    max(1_ip,number_injected_pts))
         call memory_alloca(mem_modul(1:2, modul), 'PARTICLE_DIAMETER', 'pts_injection_injectors', particle_diameter, max(1_ip,number_injected_pts))

         do iinj = 1, kfl_imax_injector_pts
            
            !
            ! Only execute injection operations for injectors
            ! where n_inj is non-zero.
            !
            if (n_inj(iinj) > 0_ip) then
               if (injection_pts(iinj) % kfl_particle_type == 0) then
                  num_types_injected = number_types_pts
               else
                  num_types_injected = 1
               end if

               !
               ! Particle diameter from individual injectors
               !
               do jj = 1, n_inj(iinj) * num_types_injected
                  nlagr_pos2 = nlagr_pos2 + 1
                  particle_diameter(nlagr_pos2) = particle_diameter_inj(iinj)%a(jj)
               enddo


               if (injection_pts(iinj) % kfl_geometry == PTS_INJ_GEO_SQUARE) then

                  !----------------------------------------------------------
                  !
                  ! Square injector given the corners
                  !
                  !----------------------------------------------------------

                  nside = injection_pts(iinj) % num_part 
                  xmini = injection_pts(iinj) % geo_coord_min(1)
                  ymini = injection_pts(iinj) % geo_coord_min(2)
                  zmini = injection_pts(iinj) % geo_coord_min(3)
                  xmaxi = injection_pts(iinj) % geo_coord_max(1)
                  ymaxi = injection_pts(iinj) % geo_coord_max(2)
                  zmaxi = injection_pts(iinj) % geo_coord_max(3)

                  if (     zmini == zmaxi) then
                     i1 = 1
                     i2 = 2
                     i3 = 3
                  else if (xmini == xmaxi) then
                     i1 = 3
                     i2 = 2
                     i3 = 1
                  else if (ymini == ymaxi) then
                     i1 = 1
                     i2 = 3
                     i3 = 2
                  else
                     call runend('LAGRAN: WRONG BOX')
                  end if
                  j1 = i1
                  j2 = i2
                  j3 = i3

                  if (xmini == xmaxi .or. ymini == ymaxi .or. zmini == zmaxi) then

                     xx(j3) = injection_pts(iinj) % geo_coord_min(i3)

                     do jj = 1, nside

                        if (nside == 1) then
                           xx(j2) = injection_pts(iinj) % geo_coord_min(i2)
                        else
                           xx(j2) = real(jj - 1, rp)/real(nside - 1, rp)*(injection_pts(iinj) % geo_coord_max(i2) &
                                - injection_pts(iinj) % geo_coord_min(i2)) + injection_pts(iinj) % geo_coord_min(i2)
                        end if

                        do ii = 1, nside

                           if (nside == 1) then
                              xx(j1) = injection_pts(iinj) % geo_coord_min(i1)
                           else
                              xx(j1) = real(ii - 1, rp)/real(nside - 1, rp)*(injection_pts(iinj) % geo_coord_max(i1) &
                                   - injection_pts(iinj) % geo_coord_min(i1)) + injection_pts(iinj) % geo_coord_min(i1)
                           end if

                           if (injection_pts(iinj) % kfl_particle_type > 0) then
                              nlagr_pos = nlagr_pos + 1
                              particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                              particle_injector(nlagr_pos) = iinj
                              particle_type(nlagr_pos) = injection_pts(iinj) % kfl_particle_type

                           else
                              do itype = 1, ntyla_pts
                                 if (parttyp(itype) % kfl_exist /= 0) then
                                    nlagr_pos = nlagr_pos + 1
                                    particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                    particle_injector(nlagr_pos) = iinj
                                    particle_type(nlagr_pos) = itype
                                 end if
                              end do
                           end if

                        end do
                     end do

                  end if

               else if (injection_pts(iinj) % kfl_geometry == PTS_INJ_GEO_RANDOM) then

                  !----------------------------------------------------------
                  !
                  ! Random injector in a square/cubic box
                  !
                  !----------------------------------------------------------

                  nside = injection_pts(iinj) % num_part
                  xmini = injection_pts(iinj) % geo_coord_min(1)
                  ymini = injection_pts(iinj) % geo_coord_min(2)
                  zmini = injection_pts(iinj) % geo_coord_min(3)
                  xmaxi = injection_pts(iinj) % geo_coord_max(1)
                  ymaxi = injection_pts(iinj) % geo_coord_max(2)
                  zmaxi = injection_pts(iinj) % geo_coord_max(3)

                  do ii = 1, n_inj(iinj)
                     U1 = random_generate_number(broadcast_seed=.true.)
                     U2 = random_generate_number(broadcast_seed=.true.)
                     U3 = random_generate_number(broadcast_seed=.true.)
                     xx(1) = xmini + (xmaxi - xmini) * U1
                     xx(2) = ymini + (ymaxi - ymini) * U2
                     if (ndime == 3) xx(3) = zmini + (zmaxi - zmini) * U3

                     if (injection_pts(iinj) % kfl_particle_type > 0) then
                        nlagr_pos = nlagr_pos + 1
                        particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                        particle_injector(nlagr_pos) = iinj
                        particle_type(nlagr_pos) = injection_pts(iinj) % kfl_particle_type
                     else
                        do itype = 1, ntyla_pts
                           if (parttyp(itype) % kfl_exist /= 0) then
                              nlagr_pos = nlagr_pos + 1
                              particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                              particle_injector(nlagr_pos) = iinj
                              particle_type(nlagr_pos) = itype
                           end if
                        end do
                     end if

                  end do

               else if (injection_pts(iinj) % kfl_geometry == PTS_INJ_GEO_RECTANGLE ) then

                  !----------------------------------------------------------
                  !
                  ! Rectangle injector
                  !
                  !----------------------------------------------------------

                  nside = injection_pts(iinj) % num_part
                  if (ndime == 3) then

                     x1 = injection_pts(iinj) % geo_coord1(1)
                     y1 = injection_pts(iinj) % geo_coord1(2)
                     z1 = injection_pts(iinj) % geo_coord1(3)
                     x2 = injection_pts(iinj) % geo_coord2(1)
                     y2 = injection_pts(iinj) % geo_coord2(2)
                     z2 = injection_pts(iinj) % geo_coord2(3)
                     x3 = injection_pts(iinj) % geo_coord3(1)
                     y3 = injection_pts(iinj) % geo_coord3(2)
                     z3 = injection_pts(iinj) % geo_coord3(3)

                     d1 = sqrt((x1 - x2)**2_ip + (y1 - y2)**2_ip + (z1 - z2)**2_ip)
                     d2 = sqrt((x3 - x2)**2_ip + (y3 - y2)**2_ip + (z3 - z2)**2_ip)
                     d3 = sqrt((x3 - x1)**2_ip + (y3 - y1)**2_ip + (z3 - z1)**2_ip)
                     if (d1 - d2 >= 0 .and. d1 - d3 >= 0)then
                        p1(1) = x1
                        p1(2) = y1
                        p1(3) = z1
                        p2(1) = x3
                        p2(2) = y3
                        p2(3) = z3
                        p3(1) = x2
                        p3(2) = y2
                        p3(3) = z2
                     else if (d2 - d1 >= 0 .and. d2 - d3 >= 0)then
                        p1(1) = x2
                        p1(2) = y2
                        p1(3) = z2
                        p2(1) = x1
                        p2(2) = y1
                        p2(3) = z1
                        p3(1) = x3
                        p3(2) = y3
                        p3(3) = z3
                     else if (d3 - d1 >= 0 .and. d3 - d2 >= 0)then
                        p1(1) = x3
                        p1(2) = y3
                        p1(3) = z3
                        p2(1) = x2
                        p2(2) = y2
                        p2(3) = z2
                        p3(1) = x1
                        p3(2) = y1
                        p3(3) = z1
                     end if
                     ux = x3 - x1
                     uy = y3 - y1
                     uz = z3 - z1
                     wx = x2 - x1
                     wy = y2 - y1
                     wz = z2 - z1
                     nn(1) = uy * wz - uz * wy
                     nn(2) = wx * uz - wz * ux
                     nn(3) = ux * wy - wx * uy
                     if (nn(1) == 0.0_rp)then
                        if (nn(3) == 0.0_rp)then
                           i1 = 1
                           i2 = 3
                           i3 = 2
                        else
                           i1 = 1
                           i2 = 2
                           i3 = 3
                           if (p1(i1) == p2(i1) .or. p2(i2) == p3(i2))then
                              i1 = 2
                              i2 = 1
                              i3 = 3
                           end if
                        end if
                     else
                        i1 = 3
                        i2 = 2
                        i3 = 1
                        if (p1(i1) == p2(i1) .or. p2(i2) == p3(i2))then
                           i1 = 2
                           i2 = 3
                           i3 = 1
                        end if

                     end if

                     do ii = 1, nside
                        if (nside == 1) then
                           xx(i1) = p2(i1)
                        else
                           if ((p2(i1) - p1(i1)) >= 0)then
                              xx(i1) = real(ii - 1, rp)/real(nside - 1, rp)*(p2(i1) - p1(i1)) + p1(i1)
                           else
                              xx(i1) = real(ii - 1, rp)/real(nside - 1, rp)*(p1(i1) - p2(i1)) + p2(i1)
                           end if
                        end if
                        do jj = 1, nside
                           if (nside == 1) then
                              xx(i2) = p2(i2)
                           else
                              if ((p3(i2) - p2(i2)) >= 0)then
                                 xx(i2) = real(jj - 1, rp)/real(nside - 1, rp)*(p3(i2) - p2(i2)) + p2(i2)
                              else
                                 xx(i2) = real(jj - 1, rp)/real(nside - 1, rp)*(p2(i2) - p3(i2)) + p3(i2)
                              end if
                           end if
                           xx(i3) = (-(xx(i1) - p1(i1)) * nn(i1) - (xx(i2) - p1(i2)) * nn(i2))/nn(i3) + p1(i3)

                           if (injection_pts(iinj) % kfl_particle_type > 0) then
                              nlagr_pos = nlagr_pos + 1
                              particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                              particle_injector(nlagr_pos) = iinj
                              particle_type(nlagr_pos) = injection_pts(iinj) % kfl_particle_type
                           else
                              do itype = 1, ntyla_pts
                                 if (parttyp(itype) % kfl_exist /= 0) then
                                    nlagr_pos = nlagr_pos + 1
                                    particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                    particle_injector(nlagr_pos) = iinj
                                    particle_type(nlagr_pos) = itype
                                 end if
                              end do
                           end if

                        end do
                     end do
                  end if

               else if ((injection_pts(iinj) % kfl_geometry == PTS_INJ_GEO_CIRCLE) .or. (injection_pts(iinj) % kfl_geometry == PTS_INJ_GEO_ANNULUS)) then
                  !----------------------------------------------------------
                  !
                  ! Circle or Annular injector
                  !
                  !----------------------------------------------------------

                  !
                  ! INPUT
                  !
                  nside     = injection_pts(iinj) % num_part
                  rad       = injection_pts(iinj) % geo_rad      ! outter radius
                  rad_min   = injection_pts(iinj) % geo_radmin   ! inner  radius


                  if (ndime == 3) then !3d mesh
                     if (injection_pts(iinj) % kfl_random == 1) then

                        if (injection_pts(iinj) % kfl_geo_spatial_dist == PTS_INJ_SPATDIST_UNIPOLAR) then !uniform points in polar coordinates
                           call livinf(-7_ip, "Injecting particles unformly distributed in polar coordinates for injector ", iinj)
                           do jj = 1, nside
                              do kk = 1, nside
                                 U1 = random_generate_number(broadcast_seed=.true.)
                                 U2 = random_generate_number(broadcast_seed=.true.)
                                 t = 2.0_rp * pi * U2
                                 radius = rad_min + (rad-rad_min) * U1   ! Don't worry, rad_min is 0.0 for circle

                                 xx(1) = 0.0_rp
                                 xx(2) = radius * cos(t)
                                 xx(3) = radius * sin(t)

                                 !
                                 ! Transforrm local coordinates back to global basis
                                 !
                                 !call maths_vector_to_new_basis(ndime, injection_pts(iinj) % geo_basis, xx)
                                 call maths_vector_from_new_basis(ndime, injection_pts(iinj) % geo_basis, xx)
                                 
                                 !
                                 ! Shift by center point
                                 !
                                 xx = xx + injection_pts(iinj) % geo_coord

                                 if (injection_pts(iinj) % kfl_particle_type > 0) then
                                    nlagr_pos                             = nlagr_pos + 1
                                    particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                    particle_injector(nlagr_pos)          = iinj
                                    particle_type(nlagr_pos)              = injection_pts(iinj) % kfl_particle_type
                                 else
                                    do itype = 1, ntyla_pts
                                       if (parttyp(itype) % kfl_exist /= 0) then
                                          nlagr_pos                             = nlagr_pos + 1
                                          particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                          particle_injector(nlagr_pos)          = iinj
                                          particle_type(nlagr_pos)              = itype
                                       end if
                                    end do
                                 end if

                              end do
                           end do

                        else if (injection_pts(iinj) % kfl_geo_spatial_dist == PTS_INJ_SPATDIST_UNICARTESIAN) then !uniform points in cartesian coordinates

                           do jj = 1, n_inj(iinj) 

                              U1     = random_generate_number(broadcast_seed=.true.) 
                              U2     = random_generate_number(broadcast_seed=.true.) 
                              
                              !
                              ! Take a random angle
                              !
                              t      = 2.0_rp * pi * U2

                              !
                              ! The idea is to give the outter radii higher
                              ! probability, becuse they have a bigger area.
                              ! The area in a drdt segment is exactly r*dr*dt.
                              ! U1 is the probability that the particle is found
                              ! within an area: 
                              ! A(radius) / Atot = (radius**2-rad_min**2) / (rad**2-rad_min**2) = U1
                              !                     -
                              !             -       |                       
                              !      -      |        |                      
                              ! o    |       |       |                      
                              !      -      |        |                      
                              !             -       |                       
                              !                     -                       
                              ! For just the circle, see:  
                              ! !https://stackoverflow.com/questions/5837572/generate-a-random-point-within-a-circle-uniformly
                              
                              radius = sqrt( rad_min**2 + (rad**2-rad_min**2) * U1 )

                              xx(1) = 0.0_rp
                              xx(2) = radius * cos(t)
                              xx(3) = radius * sin(t)

                              !
                              ! Transforrm local coordinates back to global basis
                              !
                              !!!call maths_vector_to_new_basis(ndime,injection_pts(iinj) % geo_basis,xx)
                              call maths_vector_from_new_basis(ndime,injection_pts(iinj) % geo_basis,xx)

                              !
                              ! Shift by center point
                              !
                              xx = xx + injection_pts(iinj) % geo_coord

                              if( injection_pts(iinj) % kfl_particle_type > 0 ) then
                                 nlagr_pos                             = nlagr_pos + 1
                                 particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                 particle_injector(nlagr_pos)          = iinj
                                 particle_type(nlagr_pos)              = injection_pts(iinj) % kfl_particle_type
                              else
                                 do itype = 1, ntyla_pts
                                    if (parttyp(itype) % kfl_exist /= 0) then
                                       nlagr_pos                             = nlagr_pos + 1
                                       particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                       particle_injector(nlagr_pos)          = iinj
                                       particle_type(nlagr_pos)              = itype
                                    end if
                                 end do
                              end if

                           end do


                        end if

                     else  !not injection_pts(iinj) % kfl_random == 1

                        tfact = 2.0_rp * pi / real(nside, rp)
                        rfact = (rad-rad_min) / real(nside, rp)

                        if (injection_pts(iinj) % kfl_geometry == PTS_INJ_GEO_ANNULUS) then
                           call runend('mod_pts_injection: deterministic particle palcement is not implemented for annulus.')
                        endif


                        if(injection_pts(iinj) % kfl_random==2 )then  !PIPE=Injected following mass flow rate
                           tfact = 2.0_rp * pi / real(10, rp)
                           rfact = rad / real(10, rp)
                           num_particles_per_injection =  nside / 10  !!!nside siempre es divisible entre 10!!!!
                           do kk = 1,10
                              do ii = 1, num_particles_per_injection
                                 rr = (real(kk, rp) - 1.0_rp) * rfact
                                 t = (real(kk, rp) - 1.0_rp) * tfact
                                 radius = (rr + rfact * (real(kk, rp) - 1.0_rp)/real(num_particles_per_injection,rp)) * rad

                                 xx(1) = 0.0_rp
                                 xx(2) = radius**0.5_rp * cos(t) 
                                 xx(3) = radius**0.5_rp * sin(t) 
                                !!! call maths_vector_to_new_basis(ndime, injection_pts(iinj) % geo_basis, xx)
                                 call maths_vector_from_new_basis(ndime,injection_pts(iinj) % geo_basis,xx)

                                 xx = xx + injection_pts(iinj) % geo_coord

                                 if (injection_pts(iinj) % kfl_particle_type > 0) then
                                    nlagr_pos = nlagr_pos + 1
                                    particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                    particle_injector(nlagr_pos) = iinj
                                    particle_type(nlagr_pos) = injection_pts(iinj) % kfl_particle_type
                                 else
                                    do itype = 1, ntyla_pts
                                       if (parttyp(itype) % kfl_exist /= 0) then
                                          nlagr_pos = nlagr_pos + 1
                                          particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                          particle_injector(nlagr_pos) = iinj
                                          particle_type(nlagr_pos) = itype
                                       end if
                                    end do
                                 end if
                              end do
                           end do
                        else
                           do jj = 1, nside
                              rr = (real(jj, rp) - 1.0_rp) * rfact
                              do kk = 1, nside
                                 t = (real(kk, rp) - 1.0_rp) * tfact
                                 radius = (rr + rfact * (real(kk, rp) - 1.0_rp)/real(nside,rp)) * rad

                                 xx(1) = 0.0_rp
                                 xx(2) = radius**0.5_rp * cos(t)
                                 xx(3) = radius**0.5_rp * sin(t)

                                 !!!!call maths_vector_to_new_basis(ndime, injection_pts(iinj) % geo_basis, xx)
                                 call maths_vector_from_new_basis(ndime,injection_pts(iinj) % geo_basis,xx)

                                 xx = xx + injection_pts(iinj) % geo_coord

                                 if (injection_pts(iinj) % kfl_particle_type > 0) then
                                    nlagr_pos = nlagr_pos + 1
                                    particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                    particle_injector(nlagr_pos) = iinj
                                    particle_type(nlagr_pos) = injection_pts(iinj) % kfl_particle_type
                                 else
                                    do itype = 1, ntyla_pts
                                       if (parttyp(itype) % kfl_exist /= 0) then
                                          nlagr_pos = nlagr_pos + 1
                                          particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                          particle_injector(nlagr_pos) = iinj
                                          particle_type(nlagr_pos) = itype
                                       end if
                                    end do
                                 end if

                              end do
                           end do
                        end if
                     end if
                  else !not 3d mesh
                     zfact = 2.0_rp / real(nside + 1, rp)
                     tfact = 2.0_rp * pi / real(nside + 1, rp)
                     rfact = rad / real(nside, rp)

                     if (injection_pts(iinj) % kfl_geometry == PTS_INJ_GEO_ANNULUS) then
                        call runend('mod_pts_injection: 2D model is not implemented for annulus.')
                     endif

                     if (injection_pts(iinj) % kfl_random == 1) then
                        do jj = 1, nside
                           t_aux = real(jj, rp) * tfact
                           do kk = 1, nside
                              U1 = random_generate_number(broadcast_seed=.true.)
                              U2 = random_generate_number(broadcast_seed=.true.)
                              t = t_aux * U1
                              radius = real(kk, rp) * rfact * U2
                              xx(1) = radius * cos(t) + injection_pts(iinj) % geo_coord(1)
                              xx(2) = radius * sin(t) + injection_pts(iinj) % geo_coord(2)

                              if (injection_pts(iinj) % kfl_particle_type > 0) then
                                 nlagr_pos = nlagr_pos + 1
                                 particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                 particle_injector(nlagr_pos) = iinj
                                 particle_type(nlagr_pos) = injection_pts(iinj) % kfl_particle_type
                              else
                                 do itype = 1, ntyla_pts
                                    if (parttyp(itype) % kfl_exist /= 0) then
                                       nlagr_pos = nlagr_pos + 1
                                       particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                       particle_injector(nlagr_pos) = iinj
                                       particle_type(nlagr_pos) = itype
                                    end if
                                 end do
                              end if

                           end do
                        end do
                     else if(injection_pts(iinj) % kfl_random==2 )then  !PIPE=Injected following mass flow rate
                        tfact = 2.0_rp * pi / real(10, rp)
                        rfact = rad / real(10, rp)
                        num_particles_per_injection =  nside / 10  !!!nside siempre es divisible entre 10!!!!
                        do kk = 1,10
                           do ii = 1, num_particles_per_injection
                              rr = real(kk, rp) * rfact
                              t =  real(kk, rp) * tfact
                              radius = rr
                              xx(1) = radius * cos(t) + injection_pts(iinj) % geo_coord(1)
                              xx(2) = radius * sin(t) + injection_pts(iinj) % geo_coord(2)
                             

                              if (injection_pts(iinj) % kfl_particle_type > 0) then
                                 nlagr_pos = nlagr_pos + 1
                                 particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                 particle_injector(nlagr_pos) = iinj
                                 particle_type(nlagr_pos) = injection_pts(iinj) % kfl_particle_type
                              else
                                 do itype = 1, ntyla_pts
                                    if (parttyp(itype) % kfl_exist /= 0) then
                                       nlagr_pos = nlagr_pos + 1
                                       particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                       particle_injector(nlagr_pos) = iinj
                                       particle_type(nlagr_pos) = itype
                                    end if
                                 end do
                              end if
                           end do
                        end do
                     else
                        do jj = 1, nside
                           t = real(jj, rp) * tfact
                           do kk = 1, nside
                              radius = real(kk, rp) * rfact
                              xx(1) = radius * cos(t) + injection_pts(iinj) % geo_coord(1)
                              xx(2) = radius * sin(t) + injection_pts(iinj) % geo_coord(2)

                              if (injection_pts(iinj) % kfl_particle_type > 0) then
                                 nlagr_pos = nlagr_pos + 1
                                 particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                 particle_injector(nlagr_pos) = iinj
                                 particle_type(nlagr_pos) = injection_pts(iinj) % kfl_particle_type 
                              else
                                 do itype = 1, ntyla_pts
                                    if (parttyp(itype) % kfl_exist /= 0) then
                                       nlagr_pos = nlagr_pos + 1
                                       particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                       particle_injector(nlagr_pos) = iinj
                                       particle_type(nlagr_pos) = itype
                                    end if
                                 end do
                              end if

                           end do
                        end do

                     end if
                  end if

               else if (injection_pts(iinj) % kfl_geometry == PTS_INJ_GEO_SPHERE) then

                  !----------------------------------------------------------
                  !
                  ! Sphere injector
                  !
                  !----------------------------------------------------------

                  nside     = injection_pts(iinj) % num_part
                  rad       = injection_pts(iinj) % geo_rad      ! outter radius

                  if (ndime == 3) then

                     zfact = 2.0_rp / real(nside + 1, rp)
                     tfact = 2.0_rp * pi / real(nside + 1, rp)
                     rfact = rad / real(nside, rp)

                     do ii = 1, nside

                        z = -1.0_rp + real(ii, rp) * zfact
                        r = sqrt(1.0_rp - z * z)

                        do jj = 1, nside

                           t = real(jj, rp) * tfact

                           do kk = 1, nside

                              radius = real(kk, rp) * rfact

                              xx(1) = radius * r * cos(t) 
                              xx(2) = radius * r * sin(t)
                              xx(3) = radius * z 
                              xx    = xx + injection_pts(iinj) % geo_coord

                              if (injection_pts(iinj) % kfl_particle_type > 0) then
                                 nlagr_pos = nlagr_pos + 1
                                 particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                 particle_injector(nlagr_pos) = iinj
                                 particle_type(nlagr_pos) = injection_pts(iinj) % kfl_particle_type
                              else
                                 do itype = 1, ntyla_pts
                                    if (parttyp(itype) % kfl_exist /= 0) then
                                       nlagr_pos = nlagr_pos + 1
                                       particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                       particle_injector(nlagr_pos) = iinj
                                       particle_type(nlagr_pos) = itype
                                    end if
                                 end do
                              end if

                           end do

                        end do
                     end do

                  else

                     zfact = 2.0_rp / real(nside + 1, rp)
                     tfact = 2.0_rp * pi / real(nside + 1, rp)
                     rfact = rad / real(nside, rp)

                     do jj = 1, nside

                        t = real(jj, rp) * tfact

                        do kk = 1, nside

                           radius = real(kk, rp) * rfact
                           xx(1) = radius * cos(t) + injection_pts(iinj) % geo_coord(1)
                           xx(2) = radius * sin(t) + injection_pts(iinj) % geo_coord(2)

                           if (injection_pts(iinj) % kfl_particle_type > 0) then
                              nlagr_pos = nlagr_pos + 1
                              particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                              particle_injector(nlagr_pos) = iinj
                              particle_type(nlagr_pos) = injection_pts(iinj) % kfl_particle_type
                           else
                              do itype = 1, ntyla_pts
                                 if (parttyp(itype) % kfl_exist /= 0) then
                                    nlagr_pos = nlagr_pos + 1
                                    particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                    particle_injector(nlagr_pos) = iinj
                                    particle_type(nlagr_pos) = itype
                                 end if
                              end do
                           end if

                        end do

                     end do

                  end if

               else if (injection_pts(iinj) % kfl_geometry == PTS_INJ_GEO_SEMISPHERE) then

                  !----------------------------------------------------------
                  !
                  ! Semi-Sphere injector
                  !
                  !----------------------------------------------------------

                  nside     = injection_pts(iinj) % num_part
                  rad       = injection_pts(iinj) % geo_rad      

                  if (ndime == 3) then


                     zfact = 2.0_rp / real(nside + 1, rp)
                     tfact = 2.0_rp * pi / real(nside + 1, rp)
                     rfact = rad / real(nside, rp) / 2.0_rp

                     do ii = 1, nside

                        z = -1.0_rp + real(ii, rp) * zfact
                        r = sqrt(1.0_rp - z * z)

                        do jj = 1, nside

                           t = real(jj, rp) * tfact

                           do kk = 1, nside * 2

                              radius = real(kk, rp) * rfact

                              xx(1) = radius * r * cos(t) 
                              xx(2) = radius * r * sin(t) 
                              xx(3) = radius * z 
                              xx    = xx + injection_pts(iinj) % geo_coord

                              if (dot_product(injection_pts(iinj) % geo_normal,xx-injection_pts(iinj) % geo_coord) >= 0.0_rp) then
                                 if (injection_pts(iinj) % kfl_particle_type > 0) then
                                    nlagr_pos = nlagr_pos + 1
                                    particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                    particle_injector(nlagr_pos) = iinj
                                    particle_type(nlagr_pos) = injection_pts(iinj) % kfl_particle_type
                                 else
                                    do itype = 1, ntyla_pts
                                       if (parttyp(itype) % kfl_exist /= 0) then
                                          nlagr_pos = nlagr_pos + 1
                                          particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                          particle_injector(nlagr_pos) = iinj
                                          particle_type(nlagr_pos) = itype
                                       end if
                                    end do
                                 end if

                              end if

                           end do

                        end do
                     end do

                  else

                     zfact = 2.0_rp / real(nside + 1, rp)
                     tfact = 2.0_rp * pi / real(nside + 1, rp)
                     rfact = rad / real(nside, rp) / 2.0_rp

                     do jj = 1, nside

                        t = real(jj, rp) * tfact

                        do kk = 1, nside * 2

                           radius = real(kk, rp) * rfact

                           xx(1) = radius * cos(t) + injection_pts(iinj) % geo_coord(1)
                           xx(2) = radius * sin(t) + injection_pts(iinj) % geo_coord(2)

                           if (dot_product(injection_pts(iinj) % geo_normal(1:ndime),xx(1:ndime)-injection_pts(iinj) % geo_coord(1:ndime)) >= 0.0_rp) then
                              if (injection_pts(iinj) % kfl_particle_type > 0) then
                                 nlagr_pos = nlagr_pos + 1
                                 particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                 particle_injector(nlagr_pos) = iinj
                                 particle_type(nlagr_pos) = injection_pts(iinj) % kfl_particle_type
                              else
                                 do itype = 1, ntyla_pts
                                    if (parttyp(itype) % kfl_exist /= 0) then
                                       nlagr_pos = nlagr_pos + 1
                                       particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                       particle_injector(nlagr_pos) = iinj
                                       particle_type(nlagr_pos) = itype
                                    end if
                                 end do
                              end if

                           end if

                        end do

                     end do

                  end if


               else if (injection_pts(iinj) % kfl_geometry == PTS_INJ_GEO_CONE) then

                  !----------------------------------------------------------
                  !
                  ! Cone injector
                  !
                  !           _ _     _ _
                  !    /\      |       |  z
                  !   /r \     |  hh  _|_     zfact
                  !  /<-  \    |       |  z
                  ! /______\  _|_     _|_
                  ! rad       radius = [0,rad]
                  ! <--       t      = [0,2pi]
                  !----------------------------------------------------------

                  nside     = injection_pts(iinj) % num_part
                  rad       = injection_pts(iinj) % geo_rad      
                  hh        = injection_pts(iinj) % geo_height   
                  nn        = injection_pts(iinj) % geo_normal

                  if (ndime == 3) then

                     zfact = hh/real(nside, rp)

                     do ii = 1, nside

                        z = (-1.0_rp + real(ii, rp)) * zfact
                        r = (hh - z) / (hh/rad)
                        tfact = 2 * pi /real(nside, rp) + 1
                        rfact = r / real(nside, rp)

                        do jj = 1, nside

                           t = real(jj, rp) * tfact

                           do kk = 1, nside
                              radius = real(kk, rp) * rfact

                              xx(1) = radius * cos(t)
                              xx(2) = radius * sin(t)
                              xx(3) = z
                              ! Rotate axis in function of given normal (input)
                              !!!call maths_vector_to_new_basis(ndime, injection_pts(iinj) % geo_basis, xx)
                              call maths_vector_from_new_basis(ndime,injection_pts(iinj) % geo_basis,xx)

                              xx    = xx + injection_pts(iinj) % geo_coord

                              if (injection_pts(iinj) % kfl_particle_type > 0) then
                                 nlagr_pos = nlagr_pos + 1
                                 particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                 particle_injector(nlagr_pos) = iinj
                                 particle_type(nlagr_pos) = injection_pts(iinj) % kfl_particle_type
                              else
                                 do itype = 1, ntyla_pts
                                    if (parttyp(itype) % kfl_exist /= 0) then
                                       nlagr_pos = nlagr_pos + 1
                                       particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                                       particle_injector(nlagr_pos) = iinj
                                       particle_type(nlagr_pos) = itype
                                    end if
                                 end do
                              end if

                           end do

                        end do
                     end do
                  else
                     call runend('PTS_INJECT: CONE INJECTOR ONLY WORKS IF 3D')
                  end if

               else if (injection_pts(iinj) % kfl_geometry == PTS_INJ_GEO_POINT) then

                  !----------------------------------------------------------
                  !
                  ! Pointwise injector
                  !
                  !----------------------------------------------------------

                  xx(1:3) = injection_pts(iinj) % geo_coord

                  do ii = 1, n_inj(iinj)
                     if (injection_pts(iinj) % kfl_particle_type > 0) then
                        nlagr_pos = nlagr_pos + 1
                        particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                        particle_injector(nlagr_pos) = iinj
                        particle_type(nlagr_pos) = injection_pts(iinj) % kfl_particle_type
                     else
                        do itype = 1, ntyla_pts
                           if (parttyp(itype) % kfl_exist /= 0) then
                              nlagr_pos = nlagr_pos + 1
                              particle_position(1:ndime, nlagr_pos) = xx(1:ndime)
                              particle_injector(nlagr_pos) = iinj
                              particle_type(nlagr_pos) = itype
                           end if
                        end do
                     end if
                  end do

               else if (injection_pts(iinj) % kfl_geometry == PTS_INJ_GEO_SEGMENT) then

                  !----------------------------------------------------------
                  !
                  ! Segment injector
                  !
                  !----------------------------------------------------------

                  nside = injection_pts(iinj) % num_part
                  xmini = injection_pts(iinj) % geo_coord_min(1)
                  ymini = injection_pts(iinj) % geo_coord_min(2)
                  zmini = injection_pts(iinj) % geo_coord_min(3)
                  xmaxi = injection_pts(iinj) % geo_coord_max(1)
                  ymaxi = injection_pts(iinj) % geo_coord_max(2)
                  zmaxi = injection_pts(iinj) % geo_coord_max(3)

                  if( injection_pts(iinj) % kfl_particle_type > 0 ) then
                     do ii = 1,nside
                        nlagr_pos = nlagr_pos + 1
                        particle_position(1, nlagr_pos) = xmini + real(ii-1,rp)/real(nside-1,rp)*(xmaxi - xmini)
                        particle_position(2, nlagr_pos) = ymini + real(ii-1,rp)/real(nside-1,rp)*(ymaxi - ymini)
                        if (ndime == 3) particle_position(3, nlagr_pos) = zmini + real(ii-1,rp)/real(nside-1,rp)*(zmaxi - zmini)
                        particle_injector(nlagr_pos) = iinj
                        particle_type(nlagr_pos) = injection_pts(iinj) % kfl_particle_type
                     end do
                  else
                     do itype = 1, ntyla_pts
                        if (parttyp(itype) % kfl_exist /= 0) then
                           do ii = 1, nside
                              nlagr_pos = nlagr_pos + 1
                              particle_position(1, nlagr_pos) = xmini + real(ii-1,rp)/real(nside-1,rp) * (xmaxi - xmini)
                              particle_position(2, nlagr_pos) = ymini + real(ii-1,rp)/real(nside-1,rp) * (ymaxi - ymini)
                              if (ndime == 3) particle_position(3, nlagr_pos) = zmini + real(ii-1,rp)/real(nside-1,rp) * (zmaxi - zmini)
                              particle_injector(nlagr_pos) = iinj
                              particle_type(nlagr_pos) = itype
                           end do
                        end if
                     end do
                  end if

               else if (injection_pts(iinj) % kfl_geometry == PTS_INJ_GEO_FILE) then

                  !----------------------------------------------------------
                  !
                  ! Read from file
                  !
                  !----------------------------------------------------------
                  do ii = 1, n_inj(iinj)
                     if (injection_pts(iinj) % kfl_particle_type > 0) then
                        nlagr_pos = nlagr_pos + 1
                        particle_position(1:ndime, nlagr_pos) = injection_pts(iinj) % coord_list(1:ndime,ii)
                        particle_injector(nlagr_pos) = iinj
                        particle_type(nlagr_pos) = injection_pts(iinj) % kfl_particle_type
                     else
                        do itype = 1, ntyla_pts
                           if (parttyp(itype) % kfl_exist /= 0) then
                              nlagr_pos = nlagr_pos + 1
                              particle_position(1:ndime, nlagr_pos )= injection_pts(iinj) % coord_list(1:ndime,ii)
                              particle_injector(nlagr_pos) = iinj
                              particle_type(nlagr_pos) = itype
                           end if
                        end do
                     end if
                  end do
               end if
            end if
         end do

      end if
      !
      ! Loop over particles to find host elements
      ! NLAGR_POS = Total number of injected particles
      ! NUMBER_INJECTED_PTS = Total number of particles (NLAGR_POS)
      !
      call pts_injection_ini_loc_and_cond(&
           nlagr_pos,nlagr_new,number_injected_pts,particle_position,particle_injector,&
           particle_type,particle_diameter)
      !
      ! Output injection info
      !
      call messages_live(&
           'INJECT '//integer_to_string(nlagr_new)//&
           ' LAGRANGIAN PARTICLES OUT OF '//integer_to_string(number_injected_pts),&
           'MODULE')
      if( nlagr_new < number_injected_pts ) then
         dummi = number_injected_pts - nlagr_new
         call messages_live(&
              'PARTICLES INJECTED OUT OF THE DOMAIN: '//integer_to_string(dummi),&
              'MODULE')
      end if
      !
      ! Actual total number of particles (over all partitions)
      !
      nlacc_pts = nlacc_pts + nlagr_new
      !
      ! Deallocate memory
      !
      call memory_deallo(mem_modul(1:2, modul), 'PARTICLE_POSITION', 'pts_injection_injectors', particle_position)
      call memory_deallo(mem_modul(1:2, modul), 'PARTICLE_INJECTOR', 'pts_injection_injectors', particle_injector)
      call memory_deallo(mem_modul(1:2, modul), 'PARTICLE_TYPE',     'pts_injection_injectors', particle_type)
      call memory_deallo(mem_modul(1:2, modul), 'PARTICLE_PLACE',    'pts_injection_injectors', particle_place)

      call memory_deallo(mem_modul(1:2, modul), 'PARTICLE_DIAMETER', 'pts_injection_injectors', particle_diameter)
      call memory_deallo(mem_modul(1:2,modul),  'D_LIST',            'pts_injection_injectors', particle_diameter_inj)
      call memory_deallo(mem_modul(1:2,modul),  'N_INJ'  ,           'pts_injection_injectors', n_inj)
      !
      ! Output injected particles and the recover normal flag kfl_exist
      !
      ittyp = ITASK_BEGSTE  ! AB: Manually set this, because Filter(itask) is not called before Begste
#ifdef DBPARTICLES
      call pts_output_parall_db()
      call pts_output()
#else
      call pts_output()
#endif          
      if (INOTMASTER) where( lagrtyp(:) % kfl_exist == -5) lagrtyp(:) % kfl_exist = -1
      kfl_injec = 0

      call ker_timeline('END_ASSEMBLY', nlagr_new)

      call random_end()

    end subroutine pts_injection_injectors

end module mod_pts_injection
!> @}
