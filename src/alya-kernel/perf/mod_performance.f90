!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    mod_performance.f90
!> @author  houzeaux
!> @date    2020-07-01
!> @brief   Performance
!> @details Memory and CPU performance
!-----------------------------------------------------------------------

module mod_performance

  use def_master
  use def_kermod
  use def_domain
  use def_solver
  use mod_communications,    only : PAR_AVERAGE
  use mod_communications,    only : PAR_MAX
  use mod_communications,    only : PAR_GATHER
  use mod_parall,            only : PAR_CODE_SIZE
  use mod_messages,          only : messages_live
  use mod_outfor,            only : outfor
  use def_master,            only : intost
  use mod_maths,             only : maths_solve_overdetermined_system
  use mod_timings,           only : timings_output
  use mod_alya2talp,         only : alya2talp_parallel_efficiency
  use mod_alya2talp,         only : alya2talp_MonitoringRegionStop
  use mod_iofile,            only : iofile_open_unit
  use mod_outfor,            only : outfor
  use mod_memory,            only : lun_varcount
  use mod_memory,            only : memory_output_variable_counter
  use mod_memory,            only : mem_maxim
  use def_performance
  use mod_perf_csv_writer,   only : init_perf
  use mod_perf_csv_writer,   only : write_perf_line
  use mod_perf_csv_writer,   only : write_perf_header
  use mod_strings,           only : integer_to_string
  use mod_optional_argument, only : optional_argument
  use mod_ker_proper,        only : ker_proper_updpro_time
  use mod_ann_framework,     only : ANN_STAGE_READING
  use mod_ann_framework,     only : ANN_STAGE_INPUTSCALE
  use mod_ann_framework,     only : ANN_STAGE_FORWARDPASS
  use mod_ann_framework,     only : ANN_STAGE_OUTPUTSCALE
  use mod_memory_config,     only : memory_config
  use mod_std

  implicit none

  public :: performance_outcpu
  public :: performance_outmem
  public :: performance_summary
  public :: performance_ann
  
contains

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    15/07/2015
  !> @brief   Output general CPU time info
  !> @details Output general CPU time info about starting operation and
  !>          modules
  !>
  !-----------------------------------------------------------------------

  subroutine performance_outcpu()

    real(rp)             :: error
    real(rp)             :: dummr
    integer(ip)          :: imodu,kelty,ielty,jelty,nn,mm
    real(rp)             :: PE,LB
    real(rp)             :: CE
    integer(ip), pointer :: lperm(:)
    integer(ip), pointer :: numel(:)
    integer(ip), pointer :: numel_gat(:,:)
    real(rp),    pointer :: time_gat(:)
    real(rp),    pointer :: aa(:,:)
    real(rp),    pointer :: bb(:)
    real(rp),    pointer :: xx(:)

    !----------------------------------------------------------------------
    !
    ! Compute relative weights of elements during the assembly... this
    ! is just an approximation, as it uses the different timings coming from
    ! the partitions to compute a least square problem:
    !
    ! a11 * TET04 + a12 * PYR05 + a13 * PEN06 = t1
    ! a21 * TET04 + a22 * PYR05 + a23 * PEN06 = t2
    ! ...
    ! an1 * TET04 + an2 * PYR05 + an3 * PEN06 = tn
    !
    ! Note that it does not make sense in sequential
    !
    !----------------------------------------------------------------------

    if( IPARALL ) then

       nullify(lperm,numel,numel_gat,time_gat,aa,xx,bb)
       allocate( lperm(nelty) )
       kelty = 0
       nn    = PAR_CODE_SIZE-1
       do ielty = iesta_dom,iesto_dom
          if( lexis(ielty) /= 0 ) then
             kelty = kelty + 1
             lperm(kelty) = ielty
          end if
       end do
       if( kelty > 1 ) then
          allocate( numel(kelty)          )
          allocate( numel_gat(kelty,0:nn) )
          allocate( time_gat(0:nn)        )
          numel = 0
          if( INOTMASTER .and. nelem > 0 ) then
             do jelty = 1,kelty
                ielty = lperm(jelty)
                numel(jelty) = count(ltype(1:nelem)==ielty,KIND=ip)
             end do
          end if
          call PAR_GATHER(numel,numel_gat,'IN MY CODE')
          mm = kelty
          if( INOTSLAVE ) then
             allocate( aa(nn,mm) )
             allocate( bb(nn)    )
             allocate( xx(mm)    )
          end if
          call outfor(96_ip)
          do imodu = 1,mmodu-1
             if( kfl_modul(imodu) /= 0 ) then
                dummr = cpu_modul(CPU_ASSEMBLY,imodu)
                call PAR_GATHER(dummr,time_gat,'IN MY CODE')
                if( INOTSLAVE ) then
                   bb(1:nn) = time_gat(1:nn)
                   do jelty = 1,kelty
                      aa(1:nn,jelty) = real(numel_gat(jelty,1:nn),rp)
                   end do
                   call maths_solve_overdetermined_system(nn,mm,aa,bb,xx,error)
                   if( error >= 0.0_rp ) then
                      xx            = xx/(xx(1)+zeror)
                      ioutp(1)      = imodu
                      ioutp(2)      = kelty
                      routp(1)      = error
                      routp(2:mm+1) = xx(1:mm)
                      call outfor(97_ip,INT_LIST=lperm)
                   end if
                end if
             end if
          end do
          deallocate( numel     )
          deallocate( numel_gat )
          deallocate( time_gat  )
          if( INOTSLAVE ) then
             deallocate( aa )
             deallocate( bb )
             deallocate( xx )
          end if
       end if
       deallocate( lperm     )
    end if

    !----------------------------------------------------------------------
    !
    ! Parallel efficiency
    !
    !----------------------------------------------------------------------

    if( IPARALL ) then

       call alya2talp_MonitoringRegionStop(GLOBAL_REGION=.true.)
       call alya2talp_parallel_efficiency(&
            PE,LB,CE,&
            time_comp_ave(0),time_mpi_ave(0),&
            time_comp_max(0),time_mpi_max(0),&
            GLOBAL_REGION=.true.)
       routp(1)    = LB
       routp(2)    = CE
       routp(3)    = PE
       routp(4)    = time_comp_ave(0)
       routp(5)    = time_comp_max(0)
       routp(6)    = time_mpi_ave(0)
       routp(7)    = time_mpi_max(0)
       lb_eff(1,0) = LB
       lb_eff(2,0) = CE
       lb_eff(3,0) = PE
       call outfor(104_ip)
       do imodu = 1,mmodu-1
          if( kfl_modul(imodu) /= 0 ) then
             call alya2talp_parallel_efficiency(&
                  PE,LB,CE,&
                  time_comp_ave(imodu),time_mpi_ave(imodu),&
                  time_comp_max(imodu),time_mpi_max(imodu),&
                  MODULE_REGION=.true.,CURRENT_MODULE=imodu)
             ioutp(1)        = imodu
             routp(1)        = LB
             routp(2)        = CE
             routp(3)        = PE
             routp(4)        = time_comp_ave(imodu)
             routp(5)        = time_comp_max(imodu)
             routp(6)        = time_mpi_ave(imodu)
             routp(7)        = time_mpi_max(imodu)
             lb_eff(1,imodu) = LB
             lb_eff(2,imodu) = CE
             lb_eff(3,imodu) = PE
             call outfor(105_ip)
          end if
       end do

    end if

    !----------------------------------------------------------------------
    !
    ! Output computing time summary
    !
    !----------------------------------------------------------------------

    call performance_summary()

    !----------------------------------------------------------------------
    !
    ! Solver statistics
    !
    !----------------------------------------------------------------------

    call performance_solver()
    
  end subroutine performance_outcpu

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-03-11
  !> @brief   Output a computing time summary
  !> @details Output a computing time summary
  !> 
  !-----------------------------------------------------------------------
  
  subroutine performance_summary()

    integer(ip) :: imodu,ivari,ii,jj
    real(rp)    :: cpu_minim,cpu_denom
    real(rp)    :: cpu_start_loc(size(cpu_start))
    real(rp)    :: cpu_refer,dummr
    real(rp)    :: cpu_modut
    
    !----------------------------------------------------------------------
    !
    ! Module counters
    !
    !----------------------------------------------------------------------    
    !
    ! Compute timings for kernel timings
    !
    if( associated(times) ) then
       do ii = 1,size(times)
          if( times(ii) % parent == 'kernel' ) then
             times(ii) % used = .true.
             do jj = 1,size(times)
                if( times(jj) % parent == times(ii) % name ) then
                   times(ii) % time = times(ii) % time + times(jj) % time
                end if
             end do
          end if
       end do
    end if

    do imodu = 0,mmodu
       if( kfl_modul(imodu) /= 0 ) then
          if( associated(momod(imodu) % times) ) then
             do ii = 1,size(momod(imodu) % times)
                if( momod(imodu) % times(ii) % used ) then
                   momod(imodu) % times(ii) % time_ave = momod(imodu) % times(ii) % time
                   momod(imodu) % times(ii) % time_max = momod(imodu) % times(ii) % time
                   call PAR_MAX    (momod(imodu) % times(ii) % time_max)
                   call PAR_AVERAGE(momod(imodu) % times(ii) % time_ave)
                end if
             end do
          end if
       end if
    end do 
    
    !----------------------------------------------------------------------
    !
    ! Compute maximum times over slaves for assembly, solver and output
    !
    !----------------------------------------------------------------------

    do imodu = 1,mmodu
       if( kfl_modul(imodu) /= 0 ) then
          call timings_output(&
               imodu,&
               cpu_element(imodu),&
               cpu_max_element(imodu),&
               cpu_ave_element(imodu),&
               lb_ave_element(imodu),&
               cpu_ave_per_element(imodu), &
               cpu_max_min_element(imodu),&
               cpu_ave_min_element(imodu),&
               lb_ave_min_element(imodu),&
               cpu_delta_per(imodu),&
               cpu_boundary(imodu),&
               cpu_max_boundary(imodu),&
               cpu_ave_boundary(imodu),&
               lb_ave_boundary(imodu),&
               cpu_ave_per_boundary(imodu), &
               cpu_node(imodu),&
               cpu_max_node(imodu),&
               cpu_ave_node(imodu),&
               lb_ave_node(imodu),&
               cpu_ave_per_node(imodu),&       
               cpu_max_min_node(imodu),&
               cpu_ave_min_node(imodu),&
               lb_ave_min_node(imodu), &
               cpu_delta_node_per(imodu),&
               cpu_particle(imodu),&
               cpu_max_particle(imodu),&
               cpu_ave_particle(imodu),&
               lb_ave_particle(imodu),&
               cpu_max_solver(imodu),&
               cpu_ave_solver(imodu),&
               lb_ave_solver(imodu),&
               cpu_max_post(imodu),&
               cpu_ave_post(imodu),&
               lb_ave_post(imodu))

       end if
    end do
    !
    ! Starting operations
    !
    cpu_start_loc = cpu_start
    call PAR_MAX(9_ip,cpu_start_loc,'IN MY CODE',INCLUDE_ROOT=.true.)
    call cputim(cpu_refer)
    cpu_total =   cpu_refer - cpu_initi
    call PAR_MAX(cpu_total,'IN MY CODE',INCLUDE_ROOT=.true.)

    !----------------------------------------------------------------------
    !
    ! Starting operations
    !
    !----------------------------------------------------------------------
    !
    ! Initializations
    !
    cpu_modut = 0.0_rp
    !
    ! Total CPU and CPU for starting operations
    !
    cpu_total = cpu_total + zeror
    routp( 1) = cpu_total

    routp( 2) =   cpu_start_loc(CPU_READ_GEO)            + cpu_start_loc(CPU_READ_SETS)           &
         &      + cpu_start_loc(CPU_READ_BCS)            + cpu_start_loc(CPU_READ_FIELDS)         &
         &      + cpu_start_loc(CPU_MESH_PARTITION)      + cpu_start_loc(CPU_MESH_MULTIPLICATION) &
         &      + cpu_start_loc(CPU_CONSTRUCT_DOMAIN)    + cpu_start_loc(CPU_ADDTIONAL_ARRAYS)
    routp( 3) = 100.0_rp * routp(2) / cpu_total

    routp( 4) = cpu_start_loc(CPU_READ_GEO)
    routp( 5) = 100.0_rp * routp( 4) / cpu_total

    routp( 6) = cpu_start_loc(CPU_READ_SETS)
    routp( 7) = 100.0_rp * routp( 6) / cpu_total

    routp( 8) = cpu_start_loc(CPU_READ_BCS)
    routp( 9) = 100.0_rp * routp( 8) / cpu_total

    routp(20) = cpu_start_loc(CPU_READ_FIELDS)
    routp(21) = 100.0_rp * routp(20) / cpu_total

    routp(10) = cpu_start_loc(CPU_MESH_PARTITION)
    routp(11) = 100.0_rp * routp(10) / cpu_total

    routp(14) = cpu_start_loc(CPU_MESH_MULTIPLICATION)
    routp(15) = 100.0_rp * routp(14) / cpu_total

    routp(16) = cpu_start_loc(CPU_CONSTRUCT_DOMAIN)
    routp(17) = 100.0_rp * routp(16) / cpu_total

    routp(18) = cpu_start_loc(CPU_ADDTIONAL_ARRAYS)
    routp(19) = 100.0_rp * routp(18) / cpu_total
    !
    ! Domain
    !
    routp(22) = cpu_domain(CPU_GROUPS)
    routp(23) = 100.0_rp * routp(22) / cpu_start_loc(CPU_CONSTRUCT_DOMAIN)

    routp(24) = cpu_domain(CPU_HALOS)
    routp(25) = 100.0_rp * routp(24) / cpu_start_loc(CPU_CONSTRUCT_DOMAIN)

    routp(26) = cpu_domain(CPU_ELSEST)
    routp(27) = 100.0_rp * routp(26) / cpu_start_loc(CPU_CONSTRUCT_DOMAIN)

    routp(28) = cpu_domain(CPU_COUPLING)
    routp(29) = 100.0_rp * routp(28) / cpu_start_loc(CPU_CONSTRUCT_DOMAIN)

    routp(30) = cpu_domain(CPU_OUTPUT_DOMAIN)
    routp(31) = 100.0_rp * routp(30) / cpu_start_loc(CPU_CONSTRUCT_DOMAIN)

    call outfor( 18_ip,lun_outpu,' ')
    !
    ! Property calculations
    !
    cpu_max_prope = ker_proper_updpro_time()
    call PAR_MAX(cpu_max_prope)
    cpu_ave_prope = cpu_max_prope !TODO compute the real average
    routp(1)  = cpu_max_prope
    routp(2)  = 100.0_rp*routp(1)/cpu_total
    call outfor(108_ip)     
    !
    ! Module times
    !
    do imodu = 1,mmodu
       if( kfl_modul(imodu) /= 0 ) then
          if(  cpu_modul(CPU_COUNT_ASSEMBLY,imodu) > 0.5_rp .or. &
               cpu_modul(CPU_COUNT_BOUNDARY,imodu) > 0.5_rp .or. &
               cpu_modul(CPU_COUNT_NODE,imodu)     > 0.5_rp .or. &
               cpu_modul(CPU_COUNT_PARTICLE,imodu) > 0.5_rp ) then

             cpu_modut  = cpu_modul(CPU_TOTAL_MODULE,imodu)
             call PAR_MAX(cpu_modut)
             cpu_minim  = 1.0e-6_rp
             cpu_denom  = max(cpu_modut,cpu_minim)
             coutp(1)   = namod(imodu)

             routp( 1)  = cpu_modut                               ! Total time
             routp( 2)  = 100.0_rp*routp(1)/cpu_total             
             call outfor(19_ip)
             !
             ! Compute max and averages
             !
             cpu_max_begrun(imodu)    = cpu_modul(CPU_BEGRUN,imodu)
             cpu_max_rst_read(imodu)  = cpu_modul(CPU_READ_RESTART,imodu)
             cpu_max_iniunk(imodu)    = cpu_modul(CPU_INIUNK,imodu)
             cpu_max_begste(imodu)    = cpu_modul(CPU_BEGSTE,imodu)
             cpu_max_begite(imodu)    = cpu_modul(CPU_BEGITE,imodu)
             cpu_max_endite(imodu)    = cpu_modul(CPU_ENDITE,imodu)
             cpu_max_doiter(imodu)    = cpu_modul(CPU_DOITER,imodu)
             cpu_max_endste(imodu)    = cpu_modul(CPU_ENDSTE,imodu)
             cpu_max_rst_write(imodu) = cpu_modul(CPU_WRITE_RESTART,imodu)

             call PAR_MAX(cpu_max_begrun(imodu))
             call PAR_MAX(cpu_max_rst_read(imodu))
             call PAR_MAX(cpu_max_iniunk(imodu))
             call PAR_MAX(cpu_max_begste(imodu))             
             call PAR_MAX(cpu_max_begite(imodu))
             call PAR_MAX(cpu_max_endite(imodu))
             call PAR_MAX(cpu_max_doiter(imodu))
             call PAR_MAX(cpu_max_endste(imodu))
             call PAR_MAX(cpu_max_rst_write(imodu))

             cpu_ave_begrun(imodu)    = cpu_max_begrun(imodu)
             cpu_ave_rst_read(imodu)  = cpu_max_rst_read(imodu)
             cpu_ave_iniunk(imodu)    = cpu_max_iniunk(imodu)
             cpu_ave_begste(imodu)    = cpu_max_begste(imodu)
             cpu_ave_begite(imodu)    = cpu_max_begite(imodu)
             cpu_ave_endite(imodu)    = cpu_max_endite(imodu)
             cpu_ave_doiter(imodu)    = cpu_max_doiter(imodu)
             cpu_ave_endste(imodu)    = cpu_max_endste(imodu)
             cpu_ave_rst_write(imodu) = cpu_max_rst_write(imodu)
             !
             ! Iteration other
             !
             dummr =                                 &
                  &  cpu_modul(CPU_DOITER,imodu)     &
                  & -cpu_modul(CPU_BEGITE,imodu)     &
                  & -cpu_modul(CPU_ENDITE,imodu)     &
                  & -cpu_modul(CPU_SOLVER,imodu)     &
                  & -cpu_element(imodu)              &
                  & -cpu_boundary(imodu)             &
                  & -cpu_node(imodu)                 &
                  & -cpu_particle(imodu) 

             cpu_max_othite(imodu) = dummr
             call PAR_MAX(cpu_max_othite(imodu))
             cpu_max_othite(imodu) = max(0.0_rp,cpu_max_othite(imodu))

             cpu_ave_othite(imodu) = dummr
             call PAR_AVERAGE(cpu_ave_othite(imodu))
             !
             ! Begrun
             !
             ioutp(1)   = 2 ; 
             coutp(1)   = 'Begin run'
             routp(1)   = cpu_max_begrun(imodu)
             routp(2)   = 100.0_rp*routp(1)/cpu_denom
             call outfor(100_ip)
             !
             ! Read restart
             !
             ioutp(1)   = 2
             coutp(1)   = 'Read restart'
             routp(1)   = cpu_max_rst_read(imodu)
             routp(2)   = 100.0_rp*routp(1)/cpu_denom
             call outfor(100_ip)
             !
             ! Iniunk
             !
             ioutp(1)   = 2
             coutp(1)   = 'Initial solution'
             routp(1)   = cpu_max_iniunk(imodu)
             routp(2)   = 100.0_rp*routp(1)/cpu_denom
             call outfor(100_ip)
             !
             ! Begste
             !
             ioutp(1)   = 2
             coutp(1)   = 'Begin time step'
             routp(1)   = cpu_max_begste(imodu)
             routp(2)   = 100.0_rp*routp(1)/cpu_denom
             call outfor(100_ip)
             !
             ! Begite
             !
             ioutp(1)   = 2
             coutp(1)   = 'Begin iteration'
             routp(1)   = cpu_max_begite(imodu)
             routp(2)   = 100.0_rp*routp(1)/cpu_denom
             call outfor(100_ip)
             !
             ! Doiter (rest)
             !
             ioutp(1)   = 2
             coutp(1)   = 'Do iteration'
             routp(1)   = cpu_max_doiter(imodu)-cpu_max_begite(imodu)-cpu_max_endite(imodu)
             routp(2)   = 100.0_rp*routp(1)/cpu_denom
             call outfor(100_ip)
             !
             ! Endite
             !
             ioutp(1)   = 2
             coutp(1)   = 'End iteration'
             routp(1)   = cpu_max_endite(imodu)
             routp(2)   = 100.0_rp*routp(1)/cpu_denom
             call outfor(100_ip)
             !
             ! Endste
             !
             ioutp(1)   = 2
             coutp(1)   = 'End time step'
             routp(1)   = cpu_max_endste(imodu)
             routp(2)   = 100.0_rp*routp(1)/cpu_denom
             call outfor(100_ip)
             !
             ! Output
             !
             ioutp(1)   = 2
             coutp(1)   = 'Output'
             routp(1)   = cpu_max_post(imodu)
             routp(2)   = 100.0_rp*routp(1)/cpu_denom
             call outfor(100_ip)
             !
             ! Write restart
             !
             ioutp(1)   = 2
             coutp(1)   = 'Write restart'
             routp(1)   = cpu_max_rst_write(imodu)
             routp(2)   = 100.0_rp*routp(1)/cpu_denom
             call outfor(100_ip)
             !
             ! Details of doiter
             !
             call outfor(109_ip)
             !
             ! Assemblies
             !
             if( cpu_modul(CPU_COUNT_ASSEMBLY,imodu) > 0.5_rp ) then
                coutp(1)   = 'Matrix element assembly'
                routp(1)   = cpu_ave_element(imodu)
                routp(2)   = 100.0_rp*routp(1)/cpu_denom
                routp(3)   = cpu_max_element(imodu)
                routp(4)   = 100.0_rp*routp(3)/cpu_denom
                routp(5)   = lb_ave_element(imodu)
                routp(10)  = cpu_ave_per_element(imodu)

                if( cpu_modul(CPU_COUNT_ASSEMBLY,imodu) > 1.5_rp ) then
                   routp(6)   = cpu_ave_min_element(imodu)
                   routp(7)   = cpu_max_min_element(imodu)
                   routp(8)   = lb_ave_min_element(imodu)
                   routp(9)   = cpu_delta_per(imodu)                   ! Max variability percentage
                   call outfor(99_ip)
                else
                   ioutp(1) = 1
                   call outfor(100_ip)
                end if
             end if

             if( cpu_modul(CPU_COUNT_BOUNDARY,imodu) > 0.5_rp ) then
                ioutp(1)   = 1
                coutp(1)   = 'Matrix boundary assembly'
                routp(1)   = cpu_ave_boundary(imodu) 
                routp(2)   = 100.0_rp*routp(1)/cpu_denom
                routp(3)   = cpu_max_boundary(imodu)
                routp(4)   = 100.0_rp*routp(3)/cpu_denom
                routp(5)   = lb_ave_boundary(imodu)
                routp(10)  = cpu_ave_per_boundary(imodu)
                call outfor(100_ip)
             end if

             if( cpu_modul(CPU_COUNT_NODE,imodu) > 0.5_rp ) then
                coutp(1)   = 'Node assembly'
                routp(1)   = cpu_ave_node(imodu)
                routp(2)   = 100.0_rp*routp(1)/cpu_denom
                routp(3)   = cpu_max_node(imodu)
                routp(4)   = 100.0_rp*routp(3)/cpu_denom
                routp(5)   = lb_ave_node(imodu)
                routp(10)  = cpu_ave_per_node(imodu)

                if( cpu_modul(CPU_COUNT_NODE,imodu) > 1.5_rp ) then 
                   routp(6)   = cpu_ave_min_node(imodu)
                   routp(7)   = cpu_max_min_node(imodu)
                   routp(8)   = lb_ave_min_node(imodu)
                   routp(9)   = cpu_delta_node_per(imodu)              ! Max variability percentage
                   call outfor(99_ip)
                else
                   ioutp(1) = 1
                   call outfor(100_ip)                  
                end if
             end if

             if( cpu_modul(CPU_COUNT_PARTICLE,imodu) > 0.5_rp ) then
                ioutp(1)   = 0
                coutp(1)   = 'Particle'
                routp(1)   = cpu_ave_particle(imodu)
                routp(2)   = 100.0_rp*routp(1)/cpu_denom
                routp(3)   = cpu_max_particle(imodu)
                routp(4)   = 100.0_rp*routp(3)/cpu_denom
                routp(5)   = lb_ave_particle(imodu)
                call outfor(100_ip)
             end if
             !
             ! Solver
             !
             ioutp(1)   = 2
             coutp(1)   = 'Solver'
             routp(1)   = cpu_max_solver(imodu)
             routp(2)   = 100.0_rp*routp(1)/cpu_denom
             call outfor(100_ip)             
             !
             ! Write some reporting info
             !
             if( cpu_modul(CPU_COUNT_ASSEMBLY,imodu) > 0.5_rp ) then
                if( lb_ave_element(imodu) < 0.7_rp ) & 
                     call messages_live('YOU HAVE A HIGH LOAD IMBALANCE IN YOUR MATRIX CONSTRUCTION IN MODULE '&
                     //trim(namod(imodu)),'REPORT')
             end if

          end if
       end if
    end do

    !----------------------------------------------------------------------
    !
    ! Module customized and generic counters for module log file
    !
    !----------------------------------------------------------------------

    do imodu = 1,mmodu
       if( kfl_modul(imodu) /= 0 ) then
          if( associated(momod(imodu) % times) ) then
             cpu_modut = cpu_modul(CPU_TOTAL_MODULE,imodu)
             call outfor(110_ip,momod(imodu)%lun_outpu,' ')
             cpu_denom = 100.0_rp / cpu_modut
             if( cpu_modul(CPU_COUNT_ASSEMBLY,imodu) > 0.5_rp ) then
                coutp(1)   = 'element assembly'
                routp(1)   = cpu_max_element(imodu)
                routp(2)   = routp(1) * cpu_denom
                call outfor(111_ip,momod(imodu)%lun_outpu,' ')
             end if
             if( cpu_modul(CPU_COUNT_BOUNDARY,imodu) > 0.5_rp ) then
                ioutp(1)   = 1
                coutp(1)   = 'boundary assembly'
                routp(1)   = cpu_ave_boundary(imodu) 
                routp(2)   = routp(1) * cpu_denom
                call outfor(111_ip,momod(imodu)%lun_outpu,' ')
             end if
             if( cpu_modul(CPU_COUNT_NODE,imodu) > 0.5_rp ) then
                coutp(1)   = 'node assembly'
                routp(1)   = cpu_ave_node(imodu)
                routp(2)   = routp(1) * cpu_denom
                call outfor(111_ip,momod(imodu)%lun_outpu,' ')
             end if
             if( cpu_modul(CPU_COUNT_PARTICLE,imodu) > 0.5_rp ) then
                coutp(1)   = 'particle'
                routp(1)   = cpu_ave_particle(imodu)
                routp(2)   = rp*routp(1) * cpu_denom
                call outfor(111_ip,momod(imodu)%lun_outpu,' ')
             end if
             if( associated(momod(imodu) % solve) ) then
                solve_sol => momod(imodu) % solve
                do ivari = 1,size(solve_sol)
                   if( solve_sol(ivari) % kfl_algso/=-999 .and. solve_sol(ivari) % nsolv > 0 ) then
                      coutp(1) = 'solver '//trim(solve_sol(ivari) % wprob)
                      routp(1) = solve_sol(ivari) % cputi(1)                      
                      routp(2) = routp(1) * cpu_denom
                      call outfor(111_ip,momod(imodu)%lun_outpu,' ')                       
                   end if
                end do
             end if
             !coutp(1)   = 'output'
             !routp(1)   = cpu_max_post(imodu)
             !routp(2)   = routp(1) * cpu_denom
             !call outfor(111_ip,momod(imodu)%lun_outpu,' ')            
             do ii = 1,size(momod(imodu) % times)
                if( momod(imodu) % times(ii) % used ) then
                   routp(1) = momod(imodu) % times(ii) % time_max
                   routp(2) = routp(1) * cpu_denom
                   coutp(1) = momod(imodu) % times(ii) % name
                   call outfor(111_ip,momod(imodu)%lun_outpu,' ')
                end if
             end do
          end if
       end if
    end do
    
  end subroutine performance_summary

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-12-30
  !> @brief   Output solver
  !> @details Output solver information, SpMV and DOT statistics
  !>
  !-----------------------------------------------------------------------

  subroutine performance_solver()

    do modul = 1,mmodu
       
        if ( kfl_modul(modul) /= 0 ) then

           solve => momod(modul) % solve

           if( associated(momod(modul) % solve) ) then

              solve_sol => momod(modul) % solve
              !
              ! Output solver statistics
              !
              call outfor(-37_ip,momod(modul) % lun_outpu,' ')
              !
              ! Write tail for formatted files
              !
              call outfor(6_ip,momod(modul) % lun_outpu,' ')
           end if

        end if
     end do

  end subroutine performance_solver
  
  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-12-30
  !> @brief   Output memory
  !> @details Output information on memory required
  !>
  !-----------------------------------------------------------------------

  subroutine performance_outmem()

    real(rp)     :: rgiga,rmega,rkilo,rbyte
    integer(8)   :: imodu
    character(6) :: lbyte
    integer(ip)  :: ipass,number_passes
    real(rp)     :: r_memor_dom,r_memor_sol,r_memor_els

    if( npart > 1 ) then
       number_passes = 2
    else
       number_passes = 1
    end if
    !
    ! First pass is for Master's max memory
    ! Second pass computes the max memory over the slaves
    !
    do ipass = 1,number_passes

       ioutp(50) = ipass
       !
       ! Main memory: domain+master+solver
       !
       r_memor_dom = real(memor_dom(2),rp)
       r_memor_sol = real(memma(2) + memdi(2) + memit(2),rp)
       r_memor_els = relse(3)
       r_tomax     = r_memor_dom + r_memor_sol + r_memor_els
       max_mem     = real(mem_maxim,rp)
       ave_mem     = real(mem_maxim,rp)
       do imodu = 1,mmodu
          if( kfl_modul(imodu) /= 0 ) then
             r_mem_max_modul(imodu) = real(mem_modul(2,imodu),rp)
             r_mem_ave_modul(imodu) = real(mem_modul(2,imodu),rp)
             r_tomax                = r_tomax + r_mem_max_modul(imodu)
          end if
       end do
       !
       ! Max values
       !
       if( ipass == 2 ) then
          call PAR_MAX    (r_tomax     , 'IN MY CODE')
          call PAR_MAX    (r_memor_dom , 'IN MY CODE')
          call PAR_MAX    (r_memor_sol , 'IN MY CODE')
          call PAR_MAX    (r_memor_els , 'IN MY CODE')
          call PAR_MAX    (max_mem     , 'IN MY CODE')
          call PAR_AVERAGE(ave_mem     , 'IN MY CODE')
          do imodu = 1,mmodu
             if( kfl_modul(imodu) /= 0 ) then
                call PAR_MAX    (r_mem_max_modul(imodu),'IN MY CODE')
                call PAR_AVERAGE(r_mem_ave_modul(imodu),'IN MY CODE')
             end if
          end do
       end if
       !
       ! Gbutes, Mbytes or bytes?
       !
       rgiga = 1024_rp*1024_rp*1024_rp
       rmega = 1024_rp*1024_rp
       rkilo = 1024_rp
       if( r_tomax >= rgiga ) then
          rbyte = rgiga
          lbyte = 'Gbytes'
       else if( r_tomax >= rmega ) then
          rbyte = rmega
          lbyte = 'Mbytes'
       else if( r_tomax >= rkilo ) then
          rbyte = rkilo
          lbyte = 'kbytes'
       else
          rbyte = 1.0_rp
          lbyte = ' bytes'
       end if

       routp(1) = r_tomax     / rbyte
       coutp(1) = lbyte
       routp(2) = r_memor_dom / rbyte
       coutp(2) = lbyte
       routp(3) = r_memor_els / rbyte
       coutp(3) = lbyte

       if( INOTSLAVE ) call outfor(21_ip,lun_outpu,' ')
       !
       ! Memory depending on the module
       !
       do imodu = 1,mmodu
          if( kfl_modul(imodu) /= 0 ) then
             coutp(1) = trim(namod(imodu))
             routp(1) = r_mem_max_modul(imodu) / rbyte
             coutp(2) = lbyte
             if( INOTSLAVE ) call outfor(22_ip,lun_outpu,' ')
          end if
       end do
       !
       ! Solver and maximum memory
       !
       routp(1) = r_memor_sol / rbyte
       coutp(1) = lbyte
       if( INOTSLAVE ) call outfor(24_ip,lun_outpu,' ')

    end do

    !----------------------------------------------------------------------
    !
    ! Variable memory counter
    !
    !----------------------------------------------------------------------

    if( memory_config%varcount ) then
       call memory_output_variable_counter(lun_varcount,OUTPUT_FORMAT="('VARIABLE= ',a,'CALL= ',a,' MEMORY= ',e13.6)")
    end if

  end subroutine performance_outmem


  !
  ! Output on artficial neural network timings
  !
  subroutine performance_ann()
    use mod_ann_framework,  only : max_ann_stage
    use def_kermod,         only : max_ann_fw, ann_fw 

    integer(ip)          :: ifw
    real(rp)             :: maximum(max_ann_stage) 
    real(rp)             :: average(max_ann_stage) 
    real(rp)             :: total(max_ann_stage)

    do ifw = 1, max_ann_fw
       if (ann_fw(ifw) % index /= 0) then
          !
          ! ANN exists, so statistics can be calculated
          !
          call ann_fw(ifw) % tim_statistics(maximum, average, total)

          
          !
          ! Write ANN index
          !
          coutp(1)   = 'ANN '//trim(intost(ifw))
          ioutp(1)   = ifw
          call outfor(61_ip)

          !
          ! Reading
          !
          ioutp(1)   = 2
          coutp(1)   = 'Reading ANN file (AVG)'
          routp(1)   = average(ANN_STAGE_READING)
          routp(2)   = average(ANN_STAGE_READING) / cpu_total * 100.0_rp
          call outfor(100_ip)

          !
          ! Input scaling
          !
          ioutp(1)   = 2
          coutp(1)   = 'Scaling input    (AVG)'
          routp(1)   = average(ANN_STAGE_INPUTSCALE)
          routp(2)   = average(ANN_STAGE_INPUTSCALE) / cpu_total * 100.0_rp
          call outfor(100_ip)

          !
          ! Forward pass 
          !
          ioutp(1)   = 2
          coutp(1)   = 'Forward pass     (AVG)'
          routp(1)   = average(ANN_STAGE_FORWARDPASS)
          routp(2)   = average(ANN_STAGE_FORWARDPASS) / cpu_total * 100.0_rp
          call outfor(100_ip)

          ioutp(1)   = 2
          coutp(1)   = 'Forward pass     (MAX)'
          routp(1)   = maximum(ANN_STAGE_FORWARDPASS)
          routp(2)   = maximum(ANN_STAGE_FORWARDPASS) / cpu_total * 100.0_rp
          call outfor(100_ip)

          !
          ! Output scaling
          !
          ioutp(1)   = 2
          coutp(1)   = 'Scaling output   (AVG)'
          routp(1)   = average(ANN_STAGE_OUTPUTSCALE)
          routp(2)   = average(ANN_STAGE_OUTPUTSCALE) / cpu_total * 100.0_rp
          call outfor(100_ip)

       endif
    enddo
    
  end subroutine performance_ann

end module mod_performance
!> @}
