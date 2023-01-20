!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> Brige to DLB
!> @{
!> @file    mod_alya2dlb.f90
!> @author  houzeaux
!> @date    2019-01-07
!> @brief   Bridge to DLB
!> @details Interfaces with DLB
!-----------------------------------------------------------------------

module mod_alya2dlb

  use def_kintyp,   only :  ip
  use def_master,   only :  IMASTER,INOTMASTER
  use mod_messages, only :  messages_live
  use, intrinsic         :: ISO_C_BINDING

  implicit none
#if defined ALYA_DLB_BARRIER || defined ALYA_DLB
  include 'dlbf.h'
#endif  

  private

  public :: alya2dlb_DLB_Init
  public :: alya2dlb_DLB_Enable
  public :: alya2dlb_DLB_Disable
  public :: alya2dlb_DLB_Barrier
  public :: alya2dlb_initialization
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  bsc21943
  !> @date    2019-07-19
  !> @brief   Initialization
  !> @details Initialization of DLB whenevr DLB is used for load balance
  !>          or DLB_Barrier is used
  !> 
  !-----------------------------------------------------------------------

  subroutine alya2dlb_initialization(DLB_Init_Error,DLB_Disable_Error)

    integer(ip), intent(out), optional :: DLB_Init_Error
    integer(ip), intent(out), optional :: DLB_Disable_Error
#ifdef ALYA_DLB_BARRIER
    integer(kind=c_int)                :: ierr1
#endif
#ifdef ALYA_DLB
    integer(kind=c_int)                :: ierr2
#endif
    !
    ! Intialize DLB when using barrier
    !
#ifdef ALYA_DLB_BARRIER 
    ierr1 = alya2dlb_DLB_Init()
    if( present(DLB_Init_Error) ) DLB_Init_Error = int(ierr1,ip)
#endif
    !
    ! DLB should be disabled as we only wabnt to activate it for particular loops
    ! Master does not disble to lend its resources automatically
    !
#ifdef ALYA_DLB
    if( INOTMASTER ) ierr2 = alya2dlb_DLB_Disable()
    if( present(DLB_Disable_Error) ) DLB_Disable_Error = int(ierr2,ip)
#endif

  end subroutine alya2dlb_initialization

  !-----------------------------------------------------------------------
  !> 
  !> @author  bsc21943
  !> @date    2019-07-19
  !> @brief   DLB_Init
  !> @details Interface for DLB_Init
  !> 
  !-----------------------------------------------------------------------

  integer(ip) function alya2dlb_DLB_Init()

    integer(kind=c_int)           :: ierr
    
#ifdef ALYA_DLB_BARRIER 

    integer(kind=c_int)           :: ncpus
    type(c_ptr)              :: mask
    character(kind=c_char,len=11) :: dlb_args

    dlb_args = " --barrier"//C_NULL_CHAR
    ncpus    = 0
    mask     = C_NULL_PTR
    ierr     = DLB_Init(ncpus,mask,dlb_args) 

    if( IMASTER ) then
       if( ierr < DLB_SUCCESS ) then
          call messages_live('DLB COULD NOT BE INITIALIZED','WARNING') 
       else       
          call messages_live('DLB HAS BEEN INITIALIZED SUCCESSFULLY') 
       end if
    end if
#else
    ierr = 0
#endif 
    alya2dlb_DLB_Init = int(ierr,ip)

  end function alya2dlb_DLB_Init
 
  !-----------------------------------------------------------------------
  !> 
  !> @author  bsc21943
  !> @date    2019-07-19
  !> @brief   DLB_Enable
  !> @details Interface for DLB_Enable
  !> 
  !-----------------------------------------------------------------------

  integer(ip) function alya2dlb_DLB_Enable()

   integer(c_int) :: ierr

#ifdef ALYA_DLB
   ierr = dlb_enable()
   if( ierr < DLB_SUCCESS ) call messages_live('DLB COULD NOT BE ENABLED','WARNING') 
#else
   ierr = 0
#endif
    alya2dlb_DLB_Enable = int(ierr,ip)

  end function alya2dlb_DLB_Enable

  !-----------------------------------------------------------------------
  !> 
  !> @author  bsc21943
  !> @date    2019-07-19
  !> @brief   DLB_Disable
  !> @details Interface for DLB_Disable
  !> 
  !-----------------------------------------------------------------------

  integer(ip) function alya2dlb_DLB_Disable()

   integer(c_int) :: ierr

#ifdef ALYA_DLB
   ierr = dlb_disable()
   if( ierr < DLB_SUCCESS ) call messages_live('DLB COULD NOT BE ENABLED','WARNING') 
#else
   ierr = 0
#endif
    alya2dlb_DLB_Disable = int(ierr,ip)

  end function alya2dlb_DLB_Disable

  !-----------------------------------------------------------------------
  !> 
  !> @author  bsc21943
  !> @date    2019-07-19
  !> @brief   DLB_Disable
  !> @details Interface for DLB_Disable
  !> 
  !-----------------------------------------------------------------------

  integer(ip) function alya2dlb_DLB_Barrier()

  integer(c_int) :: ierr

#ifdef ALYA_DLB_BARRIER
    ierr = dlb_barrier()
    if( ierr < DLB_SUCCESS ) call messages_live('DLB BARRIER COULD NOT BE ENABLED','WARNING')      
#else
    ierr = 0
#endif
    alya2dlb_DLB_Barrier = int(ierr,ip)

  end function alya2dlb_DLB_Barrier

end module mod_alya2dlb
!> @}
