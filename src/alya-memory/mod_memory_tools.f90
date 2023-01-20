!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!>
!> @defgroup Memory_Toolbox
!> @{
!> @name    ToolBox for memory management
!> @file    mod_memory_tools.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for memory management
!> @details Tools for mod_memory_module
!>
!------------------------------------------------------------------------

module mod_memory_tools

  use def_kintyp_basic, only : ip,rp
  use mod_memory_config, only : memory_config 

  implicit none

  private

  real(rp),      parameter :: Kbytes           = 1024.0_rp
  real(rp),      parameter :: Mbytes           = 1024.0_rp*1024.0_rp
  real(rp),      parameter :: Gbytes           = 1024.0_rp*1024.0_rp*1024.0_rp

  integer(8),    parameter :: Kbytes_i8        = 1024_8
  integer(8),    parameter :: Mbytes_i8        = 1024_8*1024_8
  integer(8),    parameter :: Gbytes_i8        = 1024_8*1024_8*1024_8

  integer(8),    parameter :: size_i1p_type    = 8
  integer(8),    parameter :: size_i1pp_type   = 8
  integer(8),    parameter :: size_i2p_type    = 8
  integer(8),    parameter :: size_i3p_type    = 8
  integer(8),    parameter :: size_r1p_type    = 8
  integer(8),    parameter :: size_r2p_type    = 8
  integer(8),    parameter :: size_r3p_type    = 8
  integer(8),    parameter :: size_r4p_type    = 8

  character(26), parameter :: cap              = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  character(26), parameter :: low              = 'abcdefghijklmnopqrstuvwxyz'
  
  integer(ip),   parameter :: varcount_max     = 1000                          ! Number of variables to track memory 

  integer(8)               :: lbytm                                            ! Allocated bytes
  integer(8),    protected :: mem_curre                                        ! Current memory allocation in bytes
  integer(8),    protected :: mem_maxim                                        ! Maximum memory allocation in bytes
  integer(ip),   protected :: mem_alloc                                        ! Number of allocations
  integer(ip),   protected :: kfl_alloc                                        ! Deallocate before allocating
  integer(ip),   protected :: varcount_number                                  ! Number of memory-tracked variables
  integer(8),    protected :: varcount_memory(varcount_max)                    ! Memory counter of variables
  character(50), protected :: varcount_name(varcount_max)                      ! Variable names
  character(50), protected :: varcount_call(varcount_max)                      ! Variable last call
  integer(ip)              :: lun_memor                                        ! Memory file unit
  integer(ip)              :: lun_varcount                                     ! Variable memory counter file unit

  public                   :: memory_initialization                            ! Initialization
  public                   :: memory_unit                                      ! Memory scaling and units
  public                   :: memory_output_info                               ! Info about memory
  public                   :: memory_runend                                    ! End Alya
  public                   :: lbytm                                            ! Allocated bytes
  public                   :: mem_curre                                        ! Current memory allocation in bytes
  public                   :: mem_maxim                                        ! Current memory allocation in bytes
  public                   :: lun_memor                                        ! Output unit
  public                   :: lun_varcount                                     ! Variable memory counter unit
  public                   :: mem_alloc                                        ! Number of allocation
  public                   :: Kbytes                                           ! Memory units
  public                   :: Mbytes                                           ! Memory units
  public                   :: Gbytes                                           ! Memory units
  public                   :: memory_add_to_memory_counter                     ! Add bytes to memory counter
  public                   :: memory_remove_from_memory_counter                ! Add bytes to memory counter
  public                   :: memory_allocation_mode                           ! Set allocation mode
  public                   :: kfl_alloc                                        ! Allocation mode
  public                   :: memory_output_variable_counter                   ! Output variable counter
  public                   :: memory_already_associated                        ! Variable already associated
  public                   :: memory_info                                      ! Information about memory
  public                   :: memory_error                                     ! Error treatment
  public                   :: memory_counter_ini                               ! Start a memory counter
  public                   :: memory_counter_end                               ! End a memory counter
  public                   :: memory_get_variable_counter                      ! Get a variable memory counter                    
  public                   :: size_i1p_type  
  public                   :: size_i1pp_type 
  public                   :: size_i2p_type  
  public                   :: size_i3p_type  
  public                   :: size_r1p_type  
  public                   :: size_r2p_type  
  public                   :: size_r3p_type  
  public                   :: size_r4p_type  
  
contains

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Initialize module variables
  !> @details Initialize all the vairables of this module
  !>
  !----------------------------------------------------------------------

  subroutine memory_initialization()

    lbytm           = 0          ! Allocated bytes
    mem_curre       = 0          ! Current memory allocation in bytes
    mem_maxim       = 0          ! Maximum memory allocation in bytes
    mem_alloc       = 0          ! Number of allocations
    !lun_memor       = 0          ! Memory file unit
    !lun_varcount    = 0          ! Variable memory counter unit
    varcount_number = 0          ! Number of memory-t<racked variables
    varcount_memory = 0          ! Memory counter of variables
    varcount_name   = ''         ! Memory counter of variable name
    varcount_call   = ''         ! Memory counter of variable calling subroutine
    kfl_alloc       = 0          ! Do not deallocate before allocating
    
  end subroutine memory_initialization

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Start memory
  !> @details End memory
  !> 
  !-----------------------------------------------------------------------

  subroutine memory_counter_ini(memor,memor_default,MEMORY_COUNTER) 

    integer(8),                     intent(out) :: memor(2)
    integer(8),                     intent(in)  :: memor_default(2)
    integer(8),           optional, intent(in)  :: MEMORY_COUNTER(2)

    if( present(MEMORY_COUNTER) ) then
       memor = MEMORY_COUNTER
    else
       memor = memor_default
    end if

  end subroutine memory_counter_ini
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Start memory
  !> @details End memory
  !> 
  !-----------------------------------------------------------------------

  subroutine memory_counter_end(memor,memor_default,MEMORY_COUNTER)

    integer(8),                     intent(in)    :: memor(2)
    integer(8),                     intent(out)   :: memor_default(2)
    integer(8),           optional, intent(out)   :: MEMORY_COUNTER(2)
    
    if( present(MEMORY_COUNTER) ) then
       MEMORY_COUNTER = memor
    else
       memor_default  = memor
    end if

  end subroutine memory_counter_end
  
  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Variables counter
  !> @details Enable variable counter to track memory
  !>
  !----------------------------------------------------------------------

  integer(ip) function memory_enable_variable_counter()
    
    memory_config%varcount = .true.
    memory_enable_variable_counter = 0
    
  end function memory_enable_variable_counter

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Get memory counter
  !> @details Get memory counter of a given variable
  !>
  !----------------------------------------------------------------------

  function memory_get_variable_counter(varia) result(res)
    
    character(len=*),           intent(in) :: varia
    integer(8)                             :: res
    integer(ip)                            :: ii

    res = 0_8
    if( memory_config%varcount ) then       
       do ii = 1,varcount_number
          if( trim(varcount_name(ii)) == trim( memory_to_upper_case(varia)) ) then
             res = varcount_memory(ii)
             return
          end if
       end do
    end if
    
  end function memory_get_variable_counter
  
  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Output counter
  !> @details Output variable memory counter
  !>
  !----------------------------------------------------------------------

  subroutine memory_output_variable_counter(nunit,OUTPUT_FORMAT)
    
    character(len=*), optional, intent(in) :: OUTPUT_FORMAT
    integer(ip),      optional, intent(in) :: nunit
    integer(4)                             :: nunit4
    integer(ip)                            :: ii
    integer(8)                             :: memor_value
    character(6)                           :: memor_char
    real(rp)                               :: memor_factor
    real(rp)                               :: total_memory

    if( memory_config%varcount ) then

       if( present(nunit) ) then
          nunit4 = int(nunit,4)
       else
          nunit4 = int(lun_varcount,4)
       end if

       memor_value = 0_8
       
       if( present(OUTPUT_FORMAT) ) then
          do ii = 1,varcount_number
             write(nunit4,OUTPUT_FORMAT) varcount_name(ii),varcount_call(ii),real(varcount_memory(ii),rp)
             memor_value = memor_value + varcount_memory(ii)
          end do
       else
          do ii = 1,varcount_number
             write(nunit4,*) varcount_name(ii),varcount_call(ii),varcount_memory(ii)
             memor_value = memor_value + varcount_memory(ii)
          end do
       end if
       
       call memory_unit(memor_value,memor_char,memor_factor)
       total_memory = real(memor_value,rp) * memor_factor 
       write(nunit4,'(a)')         '-----------------------------'
       write(nunit4,'(a,f7.2,a)') 'TOTAL MEMORY= ',total_memory,' '//trim(memor_char)
       
       call memory_unit(mem_maxim,memor_char,memor_factor)
       total_memory = real(mem_maxim,rp) * memor_factor 
       write(nunit4,'(a)')         '-----------------------------'
       write(nunit4,'(a,f7.2,a)') 'MAXIM MEMORY= ',total_memory,' '//trim(memor_char)
       
    end if 

  end subroutine memory_output_variable_counter
  
  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Variables counter
  !> @details Enable variable counter to track memory
  !>
  !----------------------------------------------------------------------

  subroutine memory_variable_counter(vanam,vacal,lbyts)

    character(len=*), intent(in) :: vanam         !< Variable name
    character(*),     intent(in) :: vacal         !< Calling subroutine name
    integer(8),       intent(in) :: lbyts         !< Number of bytes
    character(len_trim(vanam))   :: vanam_cap
    character(len_trim(vacal))   :: vacal_low
    integer(ip)                  :: ii,jj
    
    vanam_cap = memory_to_upper_case(trim(vanam))
    vacal_low = memory_to_lower_case(trim(vacal))

    ii = 0
    loop_jj: do jj = 1,varcount_number
       if( trim(varcount_name(jj)) == trim(vanam_cap) ) then
          ii = jj
          exit loop_jj
       end if
    end do loop_jj

    if( ii <= 0 ) then
       varcount_number = varcount_number + 1
       if( varcount_number > varcount_max ) varcount_number = varcount_max
       ii = varcount_number
       varcount_name(varcount_number) = trim(vanam_cap)
       varcount_call(varcount_number) = trim(vacal_low)
    end if

    varcount_memory(ii) = varcount_memory(ii) + lbyts
    
  end subroutine memory_variable_counter

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-29
  !> @brief   Convert to upper case
  !> @details Convert to upper case
  !> 
  !-----------------------------------------------------------------------

  function memory_to_upper_case (str) Result (string)

    character(*),       intent(In) :: str
    character(len(str))            :: string
    integer(ip)                    :: ic, i

!   Capitalize each letter if it is lowecase
    string = str
    do i = 1, len_trim(str)
        ic = index(low, str(i:i))
        if (ic > 0) string(i:i) = cap(ic:ic)
    end do

  end Function memory_to_upper_case

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-29
  !> @brief   Convert to lower case
  !> @details Convert to lower case
  !> 
  !-----------------------------------------------------------------------

  function memory_to_lower_case (str) Result (string)

    character(*),       intent(In) :: str
    character(len(str))            :: string
    integer(ip)                    :: ic, i

!   Capitalize each letter if it is lowecase
    string = str
    do i = 1, len_trim(str)
        ic = index(low, str(i:i))
        if (ic > 0) string(i:i) = low(ic:ic)
    end do

  end Function memory_to_lower_case

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-10
  !> @brief   Allocation mode
  !> @details Select the allocaiton mode
  !> 
  !-----------------------------------------------------------------------

  subroutine memory_allocation_mode(message)

    character(*), intent(in) :: message

    if( trim(message) == 'DEALLOCATE BEFORE ALLOCATING' ) then
       kfl_alloc = 1
    else
       kfl_alloc = 0
    end if
    
  end subroutine memory_allocation_mode  
  
  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   This routine writes some info on memory
  !> @details This routine writes some info on memory
  !>
  !----------------------------------------------------------------------

  subroutine memory_output_info(memor,vanam,vacal,vatyp)

    character(*), intent(in)    :: vanam         !< Variable name
    character(*), intent(in)    :: vacal         !< Variable name
    integer(8),   intent(inout) :: memor(2)      !< Memory counter
    character(*), intent(in)    :: vatyp
    integer(4)                  :: lun_memor4
    real(rp)                    :: rbyte,rbyt2

    if( memory_config%output ) then

       rbyte = 1.0_rp
       rbyt2 = 1.0_rp

       lun_memor4 = int(lun_memor,4)
       mem_alloc  = mem_alloc + 1
       if( mem_alloc == 1 ) write(lun_memor4,1)
       !
       ! Write info
       !
       write(lun_memor4,2) &
            real(mem_curre,rp),real(lbytm,rp),trim(vanam),&
            trim(vatyp),trim(vacal)
       flush(lun_memor4)

    end if

1   format('# Information on memory in bytes',/,&
         & '# --|  Columns displayed:',/,&
         & '#      1. Current memory      2. Variable memory    3. Variable name                  ',/,&
         & '#      4. Variable type       5. Calling subroutine ',/,&
         & '# ')
2   format(e16.8E3,1x,e16.8E3,1x,a20,1x,a10,1x,a)

  end subroutine memory_output_info

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    19/11/2015
  !> @brief   Errro message
  !> @details Write an error message when memory could not be allocated
  !>          or deallocated
  !>
  !-----------------------------------------------------------------------

  subroutine memory_error(itask,vanam,vacal,istat)
    implicit none
    integer(ip),   intent(in) :: itask
    integer(4),    intent(in) :: istat
    integer(ip)               :: ibyte
    real(rp)                  :: rbyte
    character(6)              :: lbyte
    integer(ip)               :: ibyte_tot
    real(rp)                  :: rbyte_tot
    character(6)              :: lbyte_tot
    character*(*), intent(in) :: vanam         !< Variable name
    character*(*), intent(in) :: vacal         !< Variable name
    character(200)            :: wmess
    character(20)             :: wmes2,wmes3,wmes4

    if( itask == 0 ) then
       !
       ! Allocation
       !
       if(      lbytm >= Gbytes_i8 ) then
          rbyte = Gbytes
          lbyte = 'GBYTES'
       else if( lbytm >= Mbytes_i8 ) then
          rbyte = Mbytes
          lbyte = 'MBYTES'
       else if( lbytm >= Kbytes_i8 ) then
          rbyte = Kbytes
          lbyte = 'KBYTES'
       else
          rbyte = 1.0_rp
          lbyte = 'BYTES'
       end if
       ibyte = int(real(lbytm,rp)/rbyte,KIND=ip)
       wmes2 = memory_intost(ibyte)
       wmes3 = memory_intost(int(istat,ip))

       if(      mem_curre >= Gbytes_i8 ) then
          rbyte_tot = Gbytes
          lbyte_tot = 'GBYTES'
       else if( mem_curre >= Mbytes_i8 ) then
          rbyte_tot = Mbytes
          lbyte_tot = 'MBYTES'
       else if( mem_curre >= Kbytes_i8 ) then
          rbyte_tot = Kbytes
          lbyte_tot = 'KBYTES'
       else
          rbyte_tot = 1.0_rp
          lbyte_tot = 'BYTES'
       end if
       ibyte_tot = int(real(mem_curre,rp)/rbyte_tot)
       wmes4     = memory_intost(ibyte_tot)

       wmess = trim(vacal)&
            //': MEMORY FOR '//trim(vanam)&
            //' COULD NOT BE ALLOCATED.'&
            //' RUN TIME ERROR= '//trim(wmes3)&
            //' WHEN ALLOCATING '//trim(wmes2)//' '//trim(lbyte)&
            //' OUT OF '//trim(wmes4)//' '//trim(lbyte_tot)
       call memory_runend(trim(wmess))

    else if( itask == 1 ) then
       !
       ! Reallocation
       !
       call memory_runend(trim(vacal)//': MEMORY FOR '//trim(vanam)//' COULD NOT BE REALLOCATED')

    else
       !
       ! Deallocation
       !
       call memory_runend(trim(vacal)//': MEMORY FOR '//trim(vanam)//' COULD NOT BE DEALLOCATED')

    end if

  end subroutine memory_error

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    03/03/2016
  !> @brief   Memory info
  !> @details Returns memory scaling info
  !>
  !-----------------------------------------------------------------------

  subroutine memory_unit(memor_value,memor_char,memor_factor)
    implicit none
    integer(8),    intent(in)  :: memor_value  !< Memory in bytes
    character(6),  intent(out) :: memor_char   !< Memory unit character
    real(rp),      intent(out) :: memor_factor !< Memory scaling

    if( memor_value >= Gbytes_i8 ) then
       memor_factor = 1.0_rp / Gbytes
       memor_char   = 'Gbytes'
    else if( memor_value >= Mbytes_i8 ) then
       memor_factor = 1.0_rp / Mbytes
       memor_char   = 'Mbytes'
    else if( memor_value >= Kbytes_i8 ) then
       memor_factor = 1.0_rp / Kbytes
       memor_char   = 'kbytes'
    else
       memor_factor = 1.0_rp
       memor_char   = 'bytes '
    end if

  end subroutine memory_unit

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    19/11/2015
  !> @brief   Integer to string
  !> @details Convert an integer(ip) to a string
  !>
  !-----------------------------------------------------------------------

  function memory_intost(integ)

    integer(ip)   :: integ
    integer(4)    :: integ4
    character(20) :: memory_intost
    character(20) :: intaux

    integ4 = int(integ,4)
    write(intaux,*) integ4
    memory_intost = adjustl(intaux)

  end function memory_intost

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    19/11/2015
  !> @brief   Memory info
  !> @details This routine computes some info on memory
  !>
  !-----------------------------------------------------------------------

  subroutine memory_info(memor,vanam,vacal,vatyp)

    character(*), intent(in)    :: vanam         !< Variable name
    character(*), intent(in)    :: vacal         !< Calling subroutine name
    integer(8),   intent(inout) :: memor(2)      !< Memory counter
    character(*), intent(in)    :: vatyp
    !
    ! Update input memory counter
    !
    memor(1) = memor(1)+lbytm
    memor(2) = max(memor(2),memor(1))
    !
    ! Update this module memory counters
    !
    mem_curre = mem_curre+lbytm
    mem_maxim = max(mem_curre,mem_maxim)
    !
    ! Write memory info
    !
    call memory_output_info(memor,vanam,vacal,vatyp)
    !
    ! Memory counter of variables
    !
    if( memory_config%varcount ) call memory_variable_counter(vanam,vacal,lbytm)
    
  end subroutine memory_info
  
  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-09-19
  !> @brief   Add memory to memory counter
  !> @details Add memory to the memory counter of this module
  !>
  !-----------------------------------------------------------------------

  subroutine memory_add_to_memory_counter(lbytm_in,MEMORY_COUNTER,VARIABLE_NAME,CALLING_SUBROUTINE)

    integer(8),   intent(in)              :: lbytm_in           !< Memory in bytes
    integer(8),   intent(inout), optional :: MEMORY_COUNTER(2)  !< Memory counter
    character(*), intent(in),    optional :: VARIABLE_NAME      !< Name of variable
    character(*), intent(in),    optional :: CALLING_SUBROUTINE !< The calling subroutine
    integer(8)                            :: memor_loc(2)
    character(200)                        :: vanam_loc
    character(200)                        :: vacal_loc

    lbytm     = lbytm_in
    memor_loc = 0_8
    vanam_loc = 'UNKNOWN'
    vacal_loc = 'unknown'
    if( present(MEMORY_COUNTER)     ) memor_loc = MEMORY_COUNTER
    if( present(VARIABLE_NAME)      ) vanam_loc = trim(VARIABLE_NAME)
    if( present(CALLING_SUBROUTINE) ) vacal_loc = trim(CALLING_SUBROUTINE)

    call memory_info(memor_loc,trim(vanam_loc),trim(vacal_loc),'unknown')

  end subroutine memory_add_to_memory_counter

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-09-19
  !> @brief   Add memory to memory counter
  !> @details Add memory to the memory counter of this module
  !>
  !-----------------------------------------------------------------------

  subroutine memory_remove_from_memory_counter(lbytm_in,MEMORY_COUNTER,VARIABLE_NAME,CALLING_SUBROUTINE)

    integer(8),   intent(in)              :: lbytm_in           !< Memory in bytes
    integer(8),   intent(inout), optional :: MEMORY_COUNTER(2)  !< Memory counter
    character(*), intent(in),    optional :: VARIABLE_NAME      !< Name of variable
    character(*), intent(in),    optional :: CALLING_SUBROUTINE !< The calling subroutine
    integer(8)                            :: memor_loc(2)
    character(200)                        :: vanam_loc
    character(200)                        :: vacal_loc

    lbytm     = -lbytm_in
    memor_loc = 0_8
    vanam_loc = 'UNKNOWN'
    vacal_loc = 'unknown'
    if( present(MEMORY_COUNTER)     ) memor_loc = MEMORY_COUNTER
    if( present(VARIABLE_NAME)      ) vanam_loc = trim(VARIABLE_NAME)
    if( present(CALLING_SUBROUTINE) ) vacal_loc = trim(CALLING_SUBROUTINE)

    call memory_info(memor_loc,trim(vanam_loc),trim(vacal_loc),'unknown')

  end subroutine memory_remove_from_memory_counter

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-08-09
  !> @brief   Error message
  !> @details Error message if variable already associated
  !> 
  !-----------------------------------------------------------------------

  subroutine memory_already_associated(vanam,vacal)

    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    
    write(*,'(a,a)') 'POINTER ALREADY ASSOCIATED: ',trim(vanam),', ',trim(vacal)
    call memory_runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
    
  end subroutine memory_already_associated

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-08-09
  !> @brief   Error message
  !> @details Specific runend
  !> 
  !-----------------------------------------------------------------------

  subroutine memory_runend(messa)

#ifndef MPI_OFF
#ifdef USEMPIF08
  use mpi_f08
#else
  include 'mpif.h'
#endif
#endif

    character(*), intent(in) :: messa
    integer(4)               :: istat,ierro

    write(*,'(a)') ''
    write(*,'(a)') trim(messa)

    ierro = 1 

#ifndef MPI_OFF
    call MPI_Abort(MPI_COMM_WORLD,ierro,istat)
#endif
    stop 2

  end subroutine memory_runend

end module mod_memory_tools
!> @}

