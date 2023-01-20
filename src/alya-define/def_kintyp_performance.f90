!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kernel
!> @{
!> @file    def_kintyp_performance.f90
!> @author  guillaume
!> @date    2021-02-15
!> @brief   Performance
!> @details Performance type
!>          To add a counter:
!>
!>          1. Register it in xxx_inivar:
!>
!>          call moduls_allocate_timing(10_ip)
!>          call times(1) % reg('schur complement','iteration other','io')
!>          call times(2) % reg('wall exchange')
!>
!>          The arguments are:
!>             1. The name of the counter
!>             2. The parent name
!>             3. The category: io, computation, communication, etc.
!>
!>          The different possible parent and categories can be found
!>          in mod_perf_csv.f90 
!>
!>          2. In your favorite region, add:
!>
!>          call times(1) % ini()
!>          ... things 1 ...
!>          call times(1) % add()
!>          ...
!>          call times(1) % ini()
!>          ... things 2
!>          call times(1) % add()
!>
!>          !!! WARNING !!!
!>          You timings should NOT include assemblies nor solvers
!>
!>
!-----------------------------------------------------------------------

module def_kintyp_performance
  
  use def_kintyp_basic, only : ip,rp,lg

  type perf
     character(len=:), allocatable :: name
     character(len=:), allocatable :: parent
     character(len=:), allocatable :: category
     real(rp)                      :: time_ini
     real(rp)                      :: time
     real(rp)                      :: time_ave
     real(rp)                      :: time_max
     logical(lg)                   :: used
   contains
     procedure,    pass            :: ini       => perf_ini      ! Start a timing
     procedure,    pass            :: end       => perf_end      ! End a timing
     procedure,    pass            :: add       => perf_add      ! Add a timing
     procedure,    pass            :: register  => perf_register ! Register a timing
     procedure,    pass            :: init      => perf_init     ! Called internally to initialize after allocation
  end type perf

  private
  
  public :: perf
  public :: perf_alloca
  public :: perf_deallo
  
contains

  subroutine perf_register(self,wname,wpare,wcate)
    class(perf),                intent(inout) :: self
    character(len=*),           intent(in)    :: wname
    character(len=*), optional, intent(in)    :: wpare
    character(len=*), optional, intent(in)    :: wcate
    
    if( len(wname) > 20 ) call runend('DEF_KINTYP_PERFORMANCE: NAME TOO LONG, 20 CHARACTERS AT MOST')

    allocate( character(len(wname)) :: self % name )
    self % name = wname

    if( present(wpare) ) then
       allocate( character(len(wpare)) :: self % parent )
       self % parent = wpare
    else
       allocate( character(15) :: self % parent )
       self % parent = 'iteration other'
    end if
    
    if( present(wcate) ) then
       allocate( character(len(wcate)) :: self % category )
       self % category = wcate
    else
       allocate( character(2) :: self % category )
       self % category = 'cc'
    end if
    
  end subroutine perf_register

  subroutine perf_init(self)
    class(perf), intent(inout) :: self
    self % time_ini = 0.0_rp
    self % time     = 0.0_rp
    self % time_ave = 0.0_rp
    self % time_max = 0.0_rp
    self % used     = .false.
  end subroutine perf_init

  subroutine perf_ini(self)
    class(perf), intent(inout) :: self
    call cputim(self % time_ini)   
  end subroutine perf_ini

  subroutine perf_add(self)
    class(perf), intent(inout) :: self
    real(rp)                   :: time_end
    call cputim(time_end)
    self % time     = self % time + (time_end - self % time_ini)
    self % time_ini = time_end
    self % used     = .true.
  end subroutine perf_add

  subroutine perf_end(self)
    class(perf), intent(inout) :: self
    real(rp)                   :: time_end
    call cputim(time_end)
    self % time     = (time_end - self % time_ini)
    self % time_ini = self % time
    self % used     = .true.
  end subroutine perf_end

  subroutine perf_alloca(self,nn)
    type(perf),  pointer, intent(inout) :: self(:)
    integer(ip),          intent(in)    :: nn
    integer(ip)                         :: ii
    
    if( .not. associated(self) ) then       
       allocate(self(nn))
       do ii = 1,nn
          call self(ii) % init()
       end do
    end if
    
  end subroutine perf_alloca
  
  subroutine perf_deallo(self)
    type(perf),  pointer, intent(inout) :: self(:)
        
    if( associated(self) ) then
       deallocate(self)
    end if
  end subroutine perf_deallo

end module def_kintyp_performance
!> @}

