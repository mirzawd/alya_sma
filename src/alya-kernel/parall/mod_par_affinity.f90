!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    mod_par_affinity.f90
!> @author  houzeaux
!> @date    2018-06-23
!> @brief   Affinity
!> @details Displays affinity of Alya
!-----------------------------------------------------------------------

module mod_par_affinity
  
  use def_master
  use mod_parall
  use mod_memory
  use mod_communications_global, only : PAR_ALLGATHER
  use mod_communications_tools,  only : PAR_GET_PROCESSOR_NAME
  use mod_outfor,                only : outfor
  use mod_parall,                only : PAR_WORLD_SIZE
  use mod_strings,               only : string_sort
  use, intrinsic :: ISO_C_BINDING
  
  implicit none
  
  private

  interface
     subroutine par_core_num(core_num) bind(C,name="par_core_num")
       use, intrinsic :: ISO_C_BINDING
       implicit none
       integer(C_INT), intent(out), dimension(*) :: core_num
     end subroutine par_core_num
  end interface
  
  public :: par_affinity
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-10-23
  !> @brief   Process affinity
  !> @details Compute:
  !>          CORE_NUM ... Number of the core
  !>          
  !> 
  !-----------------------------------------------------------------------

  subroutine par_affinity(core_num,host_num,num_hosts)
    
    integer(ip),         optional, pointer,   intent(inout) :: core_num(:)
    integer(ip),         optional,            intent(out)   :: host_num
    integer(ip),         optional,            intent(out)   :: num_hosts
    integer(ip),         parameter                          :: LEN_CHAR=30
    character(LEN_CHAR),           pointer                  :: hostname_gat(:)
    integer(ip)                                             :: ihost,nn
    character(len=:),              allocatable              :: hostname
    character(LEN_CHAR)                                     :: hostname_fix
    integer(C_INT),                pointer                  :: core_num4(:)
    !
    ! CORE_NUM: Get core number
    !
    if( present(core_num) ) then
       nullify(core_num4)
       if( .not. associated(core_num) ) & 
            call memory_alloca(par_memor,'core_num' ,'par_affinity',core_num ,max(1_ip,par_omp_num_threads))
       call memory_alloca(par_memor,'core_num4','par_affinity',core_num4,int(max(1_ip,par_omp_num_threads),4_ip))
       call par_core_num(core_num4)   
       core_num = int(core_num4,KIND=ip)
       call memory_deallo(par_memor,'core_num4','par_affinity',core_num4)
    end if

    if( present(host_num) .and. present(num_hosts) ) then
       !
       ! HOSTNAME: Get hostname and avoid null name
       !
       nullify(hostname_gat)
       allocate(hostname_gat(0:PAR_WORLD_SIZE-1))
       hostname = PAR_GET_PROCESSOR_NAME() 
       if(trim(hostname)=='') hostname = 'unknown'
       nn                 = len_trim(hostname)
       nn                 = min(int(len(hostname_fix),KIND=ip),nn)
       hostname_fix       = ' '
       hostname_fix(1:nn) = hostname(1:nn)       
       call PAR_ALLGATHER(hostname_fix,hostname_gat)
       !
       ! HOST_NUM: Look for my host number
       !
       call string_sort(hostname_gat,num_hosts,.false.)
       host_num = 0
       loop_ihost: do ihost = 1,num_hosts
          if( trim(hostname) == trim(hostname_gat(ihost-1)) ) then
             host_num = ihost-1
             exit loop_ihost
          end if
       end do loop_ihost
       !
       ! Deallocate
       !
       deallocate(hostname_gat)
    end if

  end subroutine par_affinity
  
end module mod_par_affinity
  !> @}
  
