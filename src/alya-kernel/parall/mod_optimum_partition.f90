!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    mod_alya2pycomms.f90
!> @author  houzeaux
!> @date    2020-07-15
!> @brief   Automatic paritioning
!> @details Automatic decision ion the number of subdomains to be used.
!>          Should be compiled with -DALYA_TALP and -DALYA_PYCOMMS
!-----------------------------------------------------------------------

module mod_optimum_partition

  use def_kintyp_basic
  use def_master
  use def_inpout
  use def_parall
  use mod_iofile
  use mod_strings,            only : integer_to_string
  use mod_messages,           only : messages_live
  use mod_alya2talp,          only : alya2talp_parallel_efficiency, alya2talp_MonitoringRegionStart, alya2talp_MonitoringRegionStop
  use mod_ecoute,             only : ecoute
  use mod_communications,     only : PAR_BARRIER
  implicit none
  private

  integer(ip)             :: num_iter
  character(8), parameter :: dirname = 'PYCOMPSS'
  character(66)           :: case_dir

  public :: optimum_partition_initialization
  public :: optimum_partition_read
  public :: optimum_partition_check
  public :: optimum_partition_parall
  public :: optimum_partition_setup
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-07-15
  !> @brief   Broadcast input data
  !> @details Broadcast input data
  !> 
  !-----------------------------------------------------------------------
  
  subroutine optimum_partition_parall()

    call iexcha(optimum_partition % kfl_method)
    call iexcha(optimum_partition % frequency)
    call iexcha(optimum_partition % modul)
    call iexcha(optimum_partition % min_cores)
    call iexcha(optimum_partition % max_cores)
    call rexcha(optimum_partition % min_criterion)
    call rexcha(optimum_partition % max_criterion)
    call rexcha(optimum_partition % change_rate)

  end subroutine optimum_partition_parall
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-07-15
  !> @brief   Initialization
  !> @details Initialization of the module
  !> 
  !-----------------------------------------------------------------------
  
  subroutine optimum_partition_initialization()

    num_iter      = 0
    
    optimum_partition % kfl_method     = 0
    optimum_partition % frequency      = 10
    optimum_partition % modul          = 0      ! Global
    optimum_partition % min_criterion = 0.5_rp
    optimum_partition % max_criterion = 0.8_rp
    optimum_partition % change_rate    = 1.0_rp
    optimum_partition % min_cores      = 3
    optimum_partition % max_cores      = huge(1_ip)

  end subroutine optimum_partition_initialization

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-07-15
  !> @brief   Initialization
  !> @details Initialization of the module
  !> 
  !-----------------------------------------------------------------------
  
  subroutine optimum_partition_setup()
    use mod_iofile_basic, only : iofile_create_directory 

    if( optimum_partition % kfl_method /= 0 ) then
#if defined ALYA_SIGNAL || ALYA_TALP
       continue
#else
       call runend('COMPILE ALYA WITH -DALYA_SIGNAL -DALYA_TALP')
#endif
       case_dir = namda(1:scan(trim(namda),"/", BACK= .true.)) !TODO: this is not crossplatform
       if( IMASTER ) then
          !call execute_command_line ('mkdir -p '//trim(case_dir)//trim(dirname))
          call iofile_create_directory( trim(case_dir)//trim(dirname) )
       end if
    end if
    
  end subroutine optimum_partition_setup
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-07-15
  !> @brief   Read data
  !> @details Read data
  !> 
  !-----------------------------------------------------------------------
  
  subroutine optimum_partition_read()

    if( words(2) /= 'OFF  ' ) then
       optimum_partition % kfl_method = 1
       do while( words(1)/='ENDOP')

          if( words(1) == 'OBJEC' ) then
             !
             ! Objective
             !
             if( words(2) == 'PARAL' ) then
                optimum_partition % kfl_method = 1 ! Parallel efficiency
             else if( words(2) == 'TIME ' ) then
                optimum_partition % kfl_method = 2 ! Time
             end if

          else if( words(1) == 'MINCO' ) then
             !
             ! Minimum number of cores
             !
             optimum_partition % min_cores = getint('MINCO',3_ip,'#Minimum number of cores')

          else if( words(1) == 'MAXCO' ) then
             !
             ! Maximum number of cores
             !
             optimum_partition % max_cores = getint('MAXCO',1000000_ip,'#Maximum number of cores')

          else if( words(1) == 'MINIM' ) then
             !
             ! Minimum efficiency/time
             !
             optimum_partition % min_criterion = getrea('MINIM',0.5_rp,'#Minimum efficiency')

          else if( words(1) == 'MAXIM' ) then
             !
             ! Maximum efficiency/time
             !
             optimum_partition % max_criterion = getrea('MAXIM',0.8_rp,'#Maximum efficiency')

          else if( words(1) == 'CHANG' ) then
             !
             ! Maximum efficiency
             !
             optimum_partition % change_rate = getrea('CHANG',1.0_rp,'#Change rate')

          else if( words(1) == 'FREQU' ) then
             !
             ! Frequency
             !
             optimum_partition % frequency = getint('FREQU',10_ip,'#Frequency')
             
          else if( words(1) == 'MODUL' ) then
             !
             ! Module
             !
             optimum_partition % modul = idmod(words(2))
             
          end if
          call ecoute('PAR_REAPRO')
       end do
    else
      do while( words(1)/='ENDOP')
         call ecoute('PAR_REAPRO')
       end do        
    end if
    
  end subroutine optimum_partition_read
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-07-15
  !> @brief   Check
  !> @details Check if an automatic partitioning should be carried out
  !> 
  !-----------------------------------------------------------------------
  
  subroutine optimum_partition_check()

    character(200) :: filename_pycompss
    character(200) :: filename_report
    real(rp)       :: PE,LB,CE
!    real(rp)       :: minmax
    real(rp)       :: time_total
    real(rp)       :: time_comp
    real(rp)       :: criterion,criterion_target
    integer(ip)    :: n,imodu,ifile
    logical(lg)    :: ifoun
    
    if( optimum_partition % kfl_method /= 0 ) then
       num_iter = num_iter + 1
       if( num_iter == optimum_partition % frequency ) then
          num_iter = 0
          !
          ! Calling TALP
          !
          imodu = optimum_partition % modul
          if( imodu == 0 ) then
             call alya2talp_MonitoringRegionStop (GLOBAL_REGION=.true.)
             call alya2talp_parallel_efficiency  (PE,LB,CE,time_comp_ave=time_comp,time_total=time_total,GLOBAL_REGION=.true.)  
             call alya2talp_MonitoringRegionStart(GLOBAL_REGION=.true.)
          else
             call alya2talp_MonitoringRegionStop (MODULE_REGION=.true.,CURRENT_MODULE=imodu)
             call alya2talp_parallel_efficiency  (PE,LB,CE,time_comp_ave=time_comp,time_total=time_total,MODULE_REGION=.true.,CURRENT_MODULE=imodu)
             call alya2talp_MonitoringRegionStart(MODULE_REGION=.true.,CURRENT_MODULE=imodu)
          end if
          call messages_live('OPTIMUM PARTITION','START SECTION')
          call messages_live('MEASURED PARALLEL EFFICIENCY= '//integer_to_string(int(Pe*100.0_rp,ip))//' %')
          if( optimum_partition % kfl_method == 1 ) then
             criterion = PE
          else if( optimum_partition % kfl_method == 2 ) then
             criterion = time_total
          end if
          criterion_target = 0.5_rp*( optimum_partition % min_criterion + optimum_partition % max_criterion )
          !
          ! # cores has to be changed
          !
          if( criterion < optimum_partition % min_criterion .or. criterion > optimum_partition % max_criterion ) then
             ifile = 0
             ifoun = .true.
             do while( ifoun )
                ifile = ifile + 1             
                filename_pycompss = trim(case_dir)//trim(dirname)//'/pycompss-'//integer_to_string(ifile)//'.txt'
                if( .not. iofile_file_exists(filename_pycompss) ) ifoun = .false.
             end do
             filename_report = trim(case_dir)//trim(dirname)//'/automatic-paritioning.res'
             call PAR_BARRIER()
             !
             ! Compute number of subdomains to achieve target criterion
             !
             if( optimum_partition % kfl_method == 1 ) then
                !if ( Pe < optimum_partition % min_criterion ) then
                !   minmax = optimum_partition % min_criterion
                !else
                !   minmax = optimum_partition % max_criterion
                !end if
                !n = nint((npart/(1.0_rp+(minmax-Pe))-npart) * optimum_partition % change_rate,ip)
               n = int(real(npart,rp)*(1.0_rp-criterion_target)/(1.0_rp-PE)*PE,ip)
             else if( optimum_partition % kfl_method == 2 ) then
                n = int(real(npart,rp)*time_comp/(criterion_target-time_total+time_comp),ip)
             end if

             n = min(optimum_partition % max_cores,max(optimum_partition % min_cores,n))
             n = n - npart 

             if( n == 0 ) then
                call messages_live('PARALLEL EFFICIENCY IS ALMOST BETWEEN REQUESTED MINIMUM AND MAXIMUM')
             else
                call iofile_open_unit(lun_autom_par,trim(filename_pycompss),'PYCOMPSS DIRECTIVE')                
                optimum_partition % kfl_method = 0
                if(      n > 0 ) then
                   !
                   ! Require more cores
                   !
                   write(lun_autom_par,'(a1,i6)') '+',n
                   call messages_live('REQUIRING '//integer_to_string(n)//' MORE CORES TO PYCOMPSS')
                   
                else
                   !
                   ! Require less cores
                   !
                   write(lun_autom_par,'(a1,i6)') '-',abs(n)
                   call messages_live('REQUIRING '//integer_to_string(abs(n))//' LESS CORES TO PYCOMPSS')
                end if                
                call iofile_close_unit(lun_autom_par,trim(filename_pycompss),'PYCOMPSS DIRECTIVE')
                !
                ! Report 
                !
                call iofile_open_unit(lun_autor_par,trim(filename_report),'PYCOMPSS CAMPAIGN REPORT','unknown','formatted','append')
                if( optimum_partition % kfl_method == 1 ) then
                   write(lun_autor_par,*) npart,PE,n+npart,criterion_target
                else if( optimum_partition % kfl_method == 2 ) then
                   write(lun_autor_par,*) npart,time_total,n+npart,criterion_target
                end if
                call iofile_close_unit(lun_autor_par,trim(filename_report),'PYCOMPSS CAMPAIGN REPORT')
             end if
          else
             call messages_live('PARALLEL EFFICIENCY IS BETWEEN REQUESTED MINIMUM AND MAXIMUM')
          end if 
          call messages_live('OPTIMUM PARTITION','END SECTION')
       end if
    end if
    
  end subroutine optimum_partition_check
  
end module mod_optimum_partition
!> @}

