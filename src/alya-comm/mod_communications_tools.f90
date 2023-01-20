!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @defgroup Communication_Toolbox
!> Toolbox for MPI communication, bridge to MPI
!> @{
!> @name    Parallelization toolbox
!> @file    mod_communications_tools.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!> @brief   ToolBox for parallel communications
!> @details ToolBox for parallel communications
!------------------------------------------------------------------------

module mod_communications_tools

  use def_communications
  use def_mpi
#include "def_mpi.inc"  
#ifdef EXTRAE
  use extrae_module
#endif
  
  implicit none

  private
  !
  ! Communicator operations
  !
  interface PAR_COMM_RANK_AND_SIZE
     module procedure PAR_COMM_RANK_AND_SIZE_MPI_4,  &
          &           PAR_COMM_RANK_AND_SIZE_MPI_8,  &
          &           PAR_COMM_RANK_AND_SIZE_COMM_4, &
          &           PAR_COMM_RANK_AND_SIZE_COMM_8, &
          &           PAR_COMM_RANK_AND_SIZE_CH_4,   &
          &           PAR_COMM_RANK_AND_SIZE_CH_8
  end interface PAR_COMM_RANK_AND_SIZE
  !
  ! Split
  !
  interface PAR_COMM_SPLIT
     module procedure PAR_COMM_SPLIT_4, &
          &           PAR_COMM_SPLIT_8, &
          &           PAR_COMM_SPLIT_IP
  end interface PAR_COMM_SPLIT
  !
  ! PAR_COMM_RANK
  !
  interface PAR_COMM_RANK
     module procedure PAR_COMM_RANK_MPI,  &
          &           PAR_COMM_RANK_COMM, &
          &           PAR_COMM_RANK_CH
  end interface PAR_COMM_RANK
  !
  ! PAR_BARRIER
  !
  interface PAR_BARRIER
     module procedure PAR_BARRIER_MPI, &
          &           PAR_BARRIER_CH
  end interface PAR_BARRIER

  public :: PAR_DEFINE_COMMUNICATOR            ! Define the communicator according to some keywords
  public :: PAR_COMM_RANK_AND_SIZE             ! Give rank (and size) of a communicator
  public :: PAR_INIT                           ! Initialize MPI
  public :: PAR_COMM_SPLIT                     ! Split a communicator
  public :: PAR_COMM_FREE                      ! Free MPI communicator
  public :: PAR_BARRIER                        ! Barrier
  public :: PAR_MPI_ERROR_TO_MESSAGE           ! Transform an MPI error code into a string
  public :: PAR_MPI_RUNEND                     ! End with an MPI message
  public :: PAR_COMM_SET_ERRHANDLER            ! Set error handler
  public :: PAR_GET_PROCESSOR_NAME             ! Get hostname
  public :: PAR_COMM_RANK                      ! Get the rank and communicator
  public :: PAR_COMM_TO_INT
  public :: PAR_WAITALL                        ! Waitall for non-blocking communications
  
contains
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2022-03-18
  !> @brief   Determine rank and comunicator
  !> @details Determine rank and comunicator
  !
  !-----------------------------------------------------------------------

  subroutine PAR_COMM_RANK_MPI(PAR_COMM_TO_USE,MY_RANK,COMM)
    
    MY_MPI_COMM,          intent(out)   :: PAR_COMM_TO_USE
    integer(4),           intent(out)   :: MY_RANK
    MY_MPI_COMM,          intent(in)    :: COMM

    PAR_COMM_TO_USE = COMM
    call PAR_COMM_RANK_AND_SIZE(COMM,MY_RANK)

  end subroutine PAR_COMM_RANK_MPI

  subroutine PAR_COMM_RANK_COMM(PAR_COMM_TO_USE,MY_RANK,COMM)
    
    MY_MPI_COMM,               intent(out)   :: PAR_COMM_TO_USE
    integer(4),                intent(out)   :: MY_RANK
    type(comm_data_par_basic), intent(in)    :: COMM

    PAR_COMM_TO_USE = COMM % PAR_COMM_WORLD
    MY_RANK         = COMM % RANK4

  end subroutine PAR_COMM_RANK_COMM

  subroutine PAR_COMM_RANK_CH(PAR_COMM_TO_USE,MY_RANK,COMM)
    
    MY_MPI_COMM,                 intent(out)   :: PAR_COMM_TO_USE
    integer(4),        optional, intent(out)   :: MY_RANK
    character(LEN=*),  optional, intent(in)    :: COMM

    if( present(COMM) ) then
       
       call PAR_DEFINE_COMMUNICATOR(COMM,PAR_COMM_TO_USE,MY_RANK=MY_RANK)

    else 
       
       PAR_COMM_TO_USE = PAR_COMM_MY_CODE
       MY_RANK         = PAR_MY_CODE_RANK
       
    end if
    
  end subroutine PAR_COMM_RANK_CH

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-12-29
  !> @brief   Define communcator and communication arrays
  !> @details Define the communicator according to a keyword
  !>
  !-----------------------------------------------------------------------

  subroutine PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu,MY_RANK)

    character(*),        optional,          intent(in)    :: wherein
    MY_MPI_COMM,                            intent(out)   :: PAR_COMM_TO_USE
    type(comm_data_par), optional, pointer, intent(inout) :: commu
    integer(4),          optional,          intent(out)   :: MY_RANK
    integer(ip)                                           :: icolo,jcolo
    integer(4)                                            :: istat4

    if( present(wherein) ) then
       if( trim(wherein) == 'IN THE UNIVERSE' ) then
          !
          ! In the universe
          !
#ifndef MPI_OFF
          PAR_COMM_TO_USE = MPI_COMM_WORLD
#endif
       else if( trim(wherein) == 'IN THE WORLD' ) then
          !
          ! In the world
          !
          PAR_COMM_TO_USE = PAR_COMM_WORLD   ! Alya world          
          if( present(commu) .and. associated(PAR_COMM_MY_CODE_ARRAY) )   commu => commd
          if( present(MY_RANK) ) MY_RANK = PAR_MY_WORLD_RANK

       else if( trim(wherein) == 'IN MY CODE' ) then
          !
          ! In my code
          !
          PAR_COMM_TO_USE = PAR_COMM_MY_CODE
          if( present(commu) .and. associated(PAR_COMM_MY_CODE_ARRAY) ) commu => PAR_COMM_MY_CODE_ARRAY(1)
          if( present(MY_RANK) ) MY_RANK = PAR_MY_CODE_RANK

       else if( trim(wherein) == 'IN MY CODE WITHOUT MASTER' ) then
          !
          ! In my code without master
          !
          PAR_COMM_TO_USE = PAR_COMM_MY_CODE_WM
          if( present(commu) .and. associated(PAR_COMM_MY_CODE_ARRAY) ) commu => PAR_COMM_MY_CODE_ARRAY(1)
          if( present(MY_RANK) ) MY_RANK = PAR_MY_CODE_RANK_WM

       else if( trim(wherein) == 'IN MY ZONE' .or. trim(wherein) == 'IN CURRENT ZONE' ) then
          !
          ! In my current zone
          !
          icolo           = par_code_zone_subd_to_color(current_code,current_zone,0_ip)
          PAR_COMM_TO_USE = PAR_COMM_COLOR(icolo,icolo)
          if( present(commu) .and. associated(PAR_COMM_MY_CODE_ARRAY) ) commu => PAR_COMM_COLOR_ARRAY(icolo)
          if( present(MY_RANK) ) MY_RANK = PAR_COMM_COLOR_PERM(icolo,icolo,PAR_MY_WORLD_RANK)

       else if( trim(wherein) == 'IN MY SUBD' ) then
          !
          ! In my current subd
          !
          icolo           = par_code_zone_subd_to_color(current_code,0_ip,current_subd)
          PAR_COMM_TO_USE = PAR_COMM_COLOR(icolo,icolo)
          if( present(commu) .and. associated(PAR_COMM_MY_CODE_ARRAY) ) commu => PAR_COMM_COLOR_ARRAY(icolo)
          if( present(MY_RANK) ) MY_RANK = PAR_COMM_COLOR_PERM(icolo,icolo,PAR_MY_WORLD_RANK)

       else if( trim(wherein) == 'IN CURRENT COUPLING' ) then
          !
          ! In my current coupling
          !
          icolo = color_target
          jcolo = color_source
          PAR_COMM_TO_USE = PAR_COMM_COLOR(icolo,jcolo)
          if( present(commu) ) then
             call runend('PAR_DEFINE_COMMUNICATOR: WRONG OPTION 1')
          end if
          if( present(MY_RANK) ) MY_RANK = PAR_COMM_COLOR_PERM(icolo,jcolo,PAR_MY_WORLD_RANK)

       else if( trim(wherein) == 'IN CURRENT COLOR' .or. trim(wherein) == 'IN CURRENT TARGET COLOR' ) then
          !
          ! In my current target color
          !
          icolo           = color_target
          PAR_COMM_TO_USE = PAR_COMM_COLOR(icolo,icolo)
          if( present(commu) .and. associated(PAR_COMM_MY_CODE_ARRAY) ) commu => PAR_COMM_COLOR_ARRAY(icolo)
          if( present(MY_RANK) ) MY_RANK = PAR_COMM_COLOR_PERM(icolo,icolo,PAR_MY_WORLD_RANK)

       else if( trim(wherein) == 'IN CURRENT SOURCE COLOR' ) then
          !
          ! In my current source color
          !
          icolo           = color_source
          PAR_COMM_TO_USE = PAR_COMM_COLOR(icolo,icolo)
          if( present(commu) .and. associated(PAR_COMM_MY_CODE_ARRAY) ) commu => PAR_COMM_COLOR_ARRAY(icolo)
          if( present(MY_RANK) ) MY_RANK = PAR_COMM_COLOR_PERM(icolo,icolo,PAR_MY_WORLD_RANK)

       else if( trim(wherein) == 'IN CURRENT' ) then
          !
          ! Uses current communicator
          !
          PAR_COMM_TO_USE = PAR_COMM_CURRENT
          if( present(commu) ) then
             call runend('PAR_DEFINE_COMMUNICATOR: WRONG OPTION 2')
          end if
#ifndef MPI_OFF
          if( present(MY_RANK) ) call MPI_Comm_rank(PAR_COMM_TO_USE,MY_RANK,istat4)
#endif

       else if( trim(wherein) == 'IN SFC PARTITION' ) then
          !
          ! In SFC partition
          !
          PAR_COMM_TO_USE = PAR_COMM_SFC_WM
#ifndef MPI_OFF
          if( present(MY_RANK) ) call MPI_Comm_rank(PAR_COMM_TO_USE,MY_RANK,istat4)
#endif
          
       else if( trim(wherein) == 'IN SFC PARTITION WITH MASTER' ) then
          !
          ! In SFC partition 
          !
          PAR_COMM_TO_USE = PAR_COMM_SFC
#ifndef MPI_OFF
          if( present(MY_RANK) ) call MPI_Comm_rank(PAR_COMM_TO_USE,MY_RANK,istat4)
#endif

       else if( trim(wherein) == 'IN MPIO' ) then
          !
          ! In MPIO
          !
          PAR_COMM_TO_USE = PAR_COMM_MPIO_WM
#ifndef MPI_OFF
          if( present(MY_RANK) ) call MPI_Comm_rank(PAR_COMM_TO_USE,MY_RANK,istat4)
#endif

       else if( trim(wherein) == 'IN MPIO WITH MASTER' ) then
          !
          ! In MPIO
          !
          PAR_COMM_TO_USE = PAR_COMM_MPIO
#ifndef MPI_OFF
          if( present(MY_RANK) ) call MPI_Comm_rank(PAR_COMM_TO_USE,MY_RANK,istat4)
#endif

       else if( trim(wherein) == 'IN SHARED MEMORY WITHOUT MASTER' ) then
          !
          ! In shared memory without master
          !
          PAR_COMM_TO_USE = PAR_COMM_SHARED_WM
          if( present(commu)) then
            call runend('PAR_DEFINE_COMMUNICATOR: WRONG OPTION 3')
          end if
          if( present(MY_RANK) ) MY_RANK = PAR_SHARED_RANK_WM

       else

          call runend('PAR DEFINE COMMUNICATOR: INVALID COMMUNICATOR OPTION: '//trim(wherein))

       end if

    else

       PAR_COMM_TO_USE =  PAR_COMM_WORLD
       commu           => commd
       if( present(MY_RANK) ) MY_RANK = PAR_MY_WORLD_RANK

    end if

  end subroutine PAR_DEFINE_COMMUNICATOR
  
  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   Initialize MPI
  !> @details Initialize MPI and define communication datatypes
  !
  !-----------------------------------------------------------------------

  subroutine PAR_INIT()

    integer(4) :: istat4

#ifndef MPI_OFF
    istat4 = 0
    call MPI_Init(istat4)
    if( istat4 /= MPI_SUCCESS ) call runend('COULD NOT INITIALIZE MPI')
#endif
#ifdef EXTRAE
    call extrae_shutdown()
#endif
    !
    ! Communication datatypes for ip and rp
    !
#ifdef I8  
    PAR_INTEGER = MPI_INTEGER8
#else
    PAR_INTEGER = MPI_INTEGER4
#endif

#ifdef R4  
    PAR_REAL = MPI_REAL4
#else
    PAR_REAL = MPI_REAL8
#endif
    !
    ! NULL communicator
    !
    PAR_COMM_NULL = MPI_COMM_NULL

  end subroutine PAR_INIT

  !----------------------------------------------------------------------
  !
  ! SPLIT COMMUNICATOR
  ! IKEY should have the rank
  !
  !----------------------------------------------------------------------

  subroutine PAR_COMM_SPLIT_4(icolor,PAR_COMM_FINAL,my_new_rank,wherein)
    
    integer(4),     intent(in)           :: icolor
    MY_MPI_COMM,    intent(out)          :: PAR_COMM_FINAL
    integer(4),     intent(out)          :: my_new_rank
    character(*),   intent(in), optional :: wherein
    integer(4)                           :: ikey4,istat4,jcolor4
    MY_MPI_COMM                          :: PAR_COMM_INITIAL

#ifndef MPI_OFF
    istat4 = 0_4
    if( icolor == 0 ) then
       jcolor4 = MPI_UNDEFINED
    else
       jcolor4 = int(icolor,4)
    end if
    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_INITIAL)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_INITIAL)
    end if
    call MPI_COMM_RANK(PAR_COMM_INITIAL,ikey4,istat4)

    if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_SPLIT4: ERROR 1')
    call MPI_COMM_SPLIT(PAR_COMM_INITIAL,jcolor4,ikey4,PAR_COMM_FINAL,istat4)
    if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_SPLIT4: ERROR 2')
    if( PAR_COMM_FINAL /= MPI_COMM_NULL ) then
       call MPI_COMM_RANK (PAR_COMM_FINAL,my_new_rank,istat4)
       if( istat4 /= MPI_SUCCESS ) call runend('PAR_COMM_SPLIT4: MPI ERROR 3')
    else
       my_new_rank = -1
    end if
#else
    PAR_COMM_FINAL = PAR_COMM_NULL
#endif

  end subroutine PAR_COMM_SPLIT_4

  subroutine PAR_COMM_SPLIT_8(icolor,PAR_COMM_FINAL,my_new_rank,wherein)
    
    integer(8),     intent(in)           :: icolor
    MY_MPI_COMM,    intent(out)          :: PAR_COMM_FINAL
    integer(8),     intent(out)          :: my_new_rank
    character(*),   intent(in), optional :: wherein
    integer(4)                           :: my_new_rank4
    integer(4)                           :: ikey4,istat4,jcolor4
    MY_MPI_COMM                       :: PAR_COMM_INITIAL4
    MY_MPI_COMM                       :: PAR_COMM_FINAL4

#ifndef MPI_OFF
    istat4 = 0_4
    if( icolor == 0 ) then
       jcolor4 = MPI_UNDEFINED
    else
       jcolor4 = int(icolor,4)
    end if

    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_INITIAL4)
    call MPI_COMM_RANK (PAR_COMM_INITIAL4,ikey4,istat4)

    if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_SPLIT8: ERROR 1')
    call MPI_COMM_SPLIT(PAR_COMM_INITIAL4,jcolor4,ikey4,PAR_COMM_FINAL4,istat4)
    if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_SPLIT8: ERROR 2')
    if( PAR_COMM_FINAL4 /= MPI_COMM_NULL ) then
       call MPI_COMM_RANK (PAR_COMM_FINAL4,my_new_rank4,istat4)
       if( istat4 /= MPI_SUCCESS ) call runend('PAR_COMM_SPLIT8: MPI ERROR 3')
       PAR_COMM_FINAL = PAR_COMM_FINAL4
       my_new_rank    = int(my_new_rank4,8)
    else
       my_new_rank    = -1
    end if
#else
    PAR_COMM_FINAL = PAR_COMM_NULL
#endif

  end subroutine PAR_COMM_SPLIT_8
  
  subroutine PAR_COMM_SPLIT_IP(icolor,PAR_COMM_FINAL,RANK_FINAL,PAR_COMM_INITIAL,RANK_INITIAL)
   
    integer(ip),    intent(in)           :: icolor
    MY_MPI_COMM,    intent(out)          :: PAR_COMM_FINAL
    integer(ip),    intent(out)          :: RANK_FINAL
    MY_MPI_COMM,    intent(in)           :: PAR_COMM_INITIAL
    integer(ip),    intent(in)           :: RANK_INITIAL
    integer(4)                           :: ikey4,istat4,jcolor4

    MY_MPI_COMM                       :: PAR_COMM_FINAL4
    integer(4)                           :: RANK_FINAL4
    MY_MPI_COMM                       :: PAR_COMM_INITIAL4

#ifndef MPI_OFF
    PAR_COMM_INITIAL4 = PAR_COMM_INITIAL
    istat4 = 0_4
    if( icolor == 0 ) then
       jcolor4 = MPI_UNDEFINED
    else
       jcolor4 = int(icolor,4)
    end if
    ikey4 = int(RANK_INITIAL,4)
    call MPI_COMM_SPLIT(PAR_COMM_INITIAL4,jcolor4,ikey4,PAR_COMM_FINAL4,istat4)
    if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_SPLIT_IP: ERROR 1')
    if( PAR_COMM_FINAL4 /= MPI_COMM_NULL ) then
       call MPI_COMM_RANK (PAR_COMM_FINAL4,RANK_FINAL4,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_SPLIT_IP: ERROR 2')
    else
       RANK_FINAL4 = -1
    end if
#else
    PAR_COMM_FINAL4 = PAR_COMM_NULL
#endif

    PAR_COMM_FINAL = PAR_COMM_FINAL4
    RANK_FINAL     = int(RANK_FINAL4,ip)

  end subroutine PAR_COMM_SPLIT_IP

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-06-11
  !> @brief   Free MPI communicator
  !> @details Free MPI communicator
  !> 
  !-----------------------------------------------------------------------
  
  subroutine PAR_COMM_FREE(PAR_COMM_FINAL)
    MY_MPI_COMM,    intent(inout) :: PAR_COMM_FINAL
    integer(4)                    :: istat4
#ifndef MPI_OFF
    if( PAR_COMM_FINAL /= MPI_COMM_NULL ) call MPI_COMM_FREE(PAR_COMM_FINAL,istat4)
#endif    
  end subroutine PAR_COMM_FREE

  !----------------------------------------------------------------------
  !
  ! RANK and SIZE of a communicator
  !
  !----------------------------------------------------------------------

  subroutine PAR_COMM_RANK_AND_SIZE_MPI_4(PAR_COMM_ORIGINAL,my_rank,comm_size)

    MY_MPI_COMM,           intent(in)  :: PAR_COMM_ORIGINAL  !< Communicator
    integer(4),            intent(out) :: my_rank            !< Rank in this communicator
    integer(4), optional,  intent(out) :: comm_size          !< Size of this communicator
    integer(4)                         :: istat4

#ifndef MPI_OFF
    istat4 = 0_4
    if( present(comm_size) ) then
       call MPI_Comm_size(PAR_COMM_ORIGINAL,comm_size,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_RANK_AND_SIZE_MPI_4')
    end if
    call MPI_Comm_rank(PAR_COMM_ORIGINAL,my_rank,istat4)
    if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_RANK_AND_SIZE_MPI_4')
       
#else
    my_rank = 0
    if( present(comm_size) ) comm_size = 0
#endif

  end subroutine PAR_COMM_RANK_AND_SIZE_MPI_4

  subroutine PAR_COMM_RANK_AND_SIZE_MPI_8(PAR_COMM_ORIGINAL,my_rank,comm_size)

    MY_MPI_COMM,           intent(in)  :: PAR_COMM_ORIGINAL  !< Communicator
    integer(8),            intent(out) :: my_rank            !< Rank in this communicator
    integer(8), optional,  intent(out) :: comm_size          !< Size of this communicator
    integer(4)                         :: istat4
    integer(4)                         :: my_rank4,comm_size4

#ifndef MPI_OFF
    istat4 = 0_4
    if( present(comm_size) ) then
       call MPI_Comm_size(PAR_COMM_ORIGINAL,comm_size4,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_RANK_AND_SIZE_MPI_8')
       comm_size = int(comm_size4,8)
    end if
    call MPI_Comm_rank(PAR_COMM_ORIGINAL,my_rank4,istat4)
    if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_RANK_AND_SIZE_MPI_8')
    my_rank = int(my_rank4,8)
     
#else
    my_rank = 0
    if( present(comm_size) ) comm_size = 0
#endif

  end subroutine PAR_COMM_RANK_AND_SIZE_MPI_8

  subroutine PAR_COMM_RANK_AND_SIZE_COMM_4(PAR_COMM_ORIGINAL,my_rank,comm_size)

    class(comm_data_par_basic),     intent(in)  :: PAR_COMM_ORIGINAL  !< Communicator
    integer(4),                     intent(out) :: my_rank            !< Rank in this communicator
    integer(4), optional,           intent(out) :: comm_size          !< Size of this communicator
    integer(4)                                  :: istat4
    MY_MPI_COMM                                 :: PAR_COMM_TO_USE

#ifndef MPI_OFF
    istat4 = 0_4
    PAR_COMM_TO_USE = PAR_COMM_ORIGINAL % PAR_COMM_WORLD
    if( present(comm_size) ) then
       call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_RANK_AND_SIZE_COMM_4')
    end if
    call MPI_Comm_rank(PAR_COMM_TO_USE,my_rank,istat4)
    if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_RANK_AND_SIZE_COMM_4')
#else
    my_rank = 0
    if( present(comm_size) ) comm_size = 0
#endif

  end subroutine PAR_COMM_RANK_AND_SIZE_COMM_4

  subroutine PAR_COMM_RANK_AND_SIZE_COMM_8(PAR_COMM_ORIGINAL,my_rank,comm_size)

    class(comm_data_par_basic),     intent(in)  :: PAR_COMM_ORIGINAL  !< Communicator
    integer(8),                     intent(out) :: my_rank            !< Rank in this communicator
    integer(8), optional,           intent(out) :: comm_size          !< Size of this communicator
    integer(4)                                  :: istat4
    MY_MPI_COMM                                 :: PAR_COMM_TO_USE
    integer(4)                                  :: my_rank4,comm_size4
#ifndef MPI_OFF
    istat4 = 0_4
    PAR_COMM_TO_USE = PAR_COMM_ORIGINAL % PAR_COMM_WORLD
    if( present(comm_size) ) then
       call MPI_Comm_size(PAR_COMM_TO_USE,comm_size4,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_RANK_AND_SIZE_COMM_8')
       comm_size = int(comm_size4,8)
    end if
    call MPI_Comm_rank(PAR_COMM_TO_USE,my_rank4,istat4)
    if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_RANK_AND_SIZE_COMM_8')
    my_rank = int(my_rank4,8)
#else
    my_rank = 0
    if( present(comm_size) ) comm_size = 0
#endif

  end subroutine PAR_COMM_RANK_AND_SIZE_COMM_8

  subroutine PAR_COMM_RANK_AND_SIZE_CH_4(my_rank,comm_size,wherein)
    implicit none
    integer(4),             intent(out) :: my_rank           !< Rank in this communicator
    integer(4),   optional, intent(out) :: comm_size         !< Size of this communicator
    character(*),           intent(in)  :: wherein             !< Wherein
    integer(4)                          :: istat4
    MY_MPI_COMM                         :: PAR_COMM_ORIGINAL

#ifndef MPI_OFF
    istat4 = 0_4
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_ORIGINAL)
    if( present(comm_size) ) then
       call MPI_Comm_size(PAR_COMM_ORIGINAL,comm_size,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_RANK_AND_SIZE_CH_4')
    end if
    call MPI_Comm_rank(PAR_COMM_ORIGINAL,my_rank,istat4)
    if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_RANK_AND_SIZE_CH_4')
#else
    my_rank   = 0
    if( present(comm_size) ) comm_size = 0
#endif

  end subroutine PAR_COMM_RANK_AND_SIZE_CH_4

  subroutine PAR_COMM_RANK_AND_SIZE_CH_8(my_rank,comm_size,wherein)
    implicit none
    integer(8),             intent(out) :: my_rank           !< Rank in this communicator
    integer(8),   optional, intent(out) :: comm_size         !< Size of this communicator
    character(*),           intent(in)  :: wherein           !< Wherein
    integer(4)                          :: istat4
    MY_MPI_COMM                         :: PAR_COMM_ORIGINAL
    integer(4)                          :: my_rank4
    integer(4)                          :: comm_size4

#ifndef MPI_OFF
    istat4 = 0_4
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_ORIGINAL)
    if( present(comm_size) ) then
       call MPI_Comm_size(PAR_COMM_ORIGINAL,comm_size4,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_RANK_AND_SIZE_CH_8')
       comm_size = int(comm_size4,8)
    end if
    call MPI_Comm_rank(PAR_COMM_ORIGINAL,my_rank4,istat4)
    if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_RANK_AND_SIZE_CH_8')
    my_rank = int(my_rank4,8)
#else
    my_rank   = 0
    if( present(comm_size) ) comm_size = 0
#endif

  end subroutine PAR_COMM_RANK_AND_SIZE_CH_8

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    16/05/2014
  !> @brief   MPI Barrier
  !> @details MPI Barrier
  !
  !----------------------------------------------------------------------

  subroutine PAR_BARRIER_CH(COMM)
    character(LEN=*), optional, intent(in) :: COMM
    MY_MPI_COMM                            :: PAR_COMM_TO_USE
    integer(4)                             :: istat4

#ifndef MPI_OFF
    if( IPARALL ) then       
       if( present(COMM) ) then         
          call PAR_DEFINE_COMMUNICATOR(COMM,PAR_COMM_TO_USE)             
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if      
       call MPI_Barrier( PAR_COMM_TO_USE, istat4 )
    end if
#endif

  end subroutine PAR_BARRIER_CH
  
  subroutine PAR_BARRIER_MPI(COMM)
    MY_MPI_COMM,     intent(in) :: COMM
    integer(4)                  :: istat4

#ifndef MPI_OFF
    if( IPARALL ) then       
       call MPI_Barrier( COMM, istat4 )
    end if
#endif

  end subroutine PAR_BARRIER_MPI
  

!!$  subroutine PAR_BARRIER(wherein,PAR_COMM_IN)
!!$    character(*),   optional, intent(in) :: wherein
!!$    MY_MPI_COMM,    optional, intent(in) :: PAR_COMM_IN
!!$    MY_MPI_COMM                       :: PAR_COMM_TO_USE
!!$    integer(4)                           :: istat4
!!$
!!$#ifndef MPI_OFF
!!$    if( IPARALL ) then
!!$       if( present(PAR_COMM_IN) ) then
!!$          PAR_COMM_TO_USE=PAR_COMM_IN
!!$       else
!!$          if( present(wherein) ) then
!!$             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
!!$          else
!!$             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
!!$          end if
!!$       endif
!!$       call MPI_Barrier( PAR_COMM_TO_USE, istat4 )
!!$    end if
!!$#endif
!!$
!!$  end subroutine PAR_BARRIER

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2019-01-02
  !> @brief   Transform an MPI error code into a string
  !> @details Transform an MPI error code into a string
  !>
  !-----------------------------------------------------------------------

  function PAR_MPI_ERROR_TO_MESSAGE(istat4) result(message)

    integer(4),                    intent(in)  :: istat4   !< Error code
    character(len=:), allocatable              :: message   !< Message
    integer(4)                                 :: length4
    integer(4)                                 :: temp4

#ifndef MPI_OFF
    character(len=MPI_MAX_ERROR_STRING) :: message_mpi
    if( istat4 /= MPI_SUCCESS ) then
       call MPI_Error_string(istat4,message_mpi,length4,temp4)
       message = message_mpi(1:length4)
    end if
#endif

  end function PAR_MPI_ERROR_TO_MESSAGE

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2019-01-02
  !> @brief   End
  !> @details End with an MPI error message
  !>
  !-----------------------------------------------------------------------

  subroutine PAR_MPI_RUNEND(istat ,vacal) 

    class(*),                   intent(in) :: istat     !< Error code
    character(len=*), optional, intent(in) :: vacal     !< Caller
    integer(4)                             :: istat4
    
#ifndef MPI_OFF
    select type ( istat )
    type is ( integer(kind=4) ) ; istat4 = istat
    type is ( integer(kind=8) ) ; istat4 = int(istat,4)
    end select
       
    if( istat4 /= MPI_SUCCESS ) then
       if( present(vacal) ) then
          call runend(trim(vacal)//': '//PAR_MPI_ERROR_TO_MESSAGE(istat4))
       else
          call runend('MPI ERROR: '//PAR_MPI_ERROR_TO_MESSAGE(istat4))
       end if
    end if
#endif

  end subroutine PAR_MPI_RUNEND

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2019-01-02
  !> @brief   End
  !> @details End with an MPI error message
  !>
  !-----------------------------------------------------------------------

  subroutine PAR_COMM_SET_ERRHANDLER()

    integer(4) :: istat4

#ifndef MPI_OFF
    call MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_RETURN,istat4)
#endif

  end subroutine PAR_COMM_SET_ERRHANDLER

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2019-01-02
  !> @brief   Name of processor
  !> @details Get host name
  !>
  !-----------------------------------------------------------------------

  function PAR_GET_PROCESSOR_NAME() result(wname)

    character(len=:), allocatable         :: wname
    integer                               :: resultlen
    integer                               :: ierror
#ifndef MPI_OFF
    character(LEN=MPI_MAX_PROCESSOR_NAME) :: wname_loc
    call mpi_get_processor_name(wname_loc,resultlen,ierror)
    wname = wname_loc(1:resultlen)
#else
    wname = 'unknown'
#endif

  end function PAR_GET_PROCESSOR_NAME
  
  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2019-01-02
  !> @brief   Communicator as an integer
  !> @details Communicator as an integer
  !>
  !-----------------------------------------------------------------------

  function PAR_COMM_TO_INT(COMM) result(COMM_INT)
    MY_MPI_COMM :: COMM
    integer(4)  :: COMM_INT

#ifndef MPI_OFF
#ifdef USEMPIF08
    COMM_INT = COMM % MPI_VAL
#else
    COMM_INT = COMM
#endif
#else
     COMM_INT = PAR_COMM_NULL
#endif
     
  end function PAR_COMM_TO_INT 

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2019-01-02
  !> @brief   WAITALL Communicator as an integer
  !> @details WAITALL for non-blocking comunications
  !>
  !-----------------------------------------------------------------------

  subroutine PAR_WAITALL(count4,ireq4,ipass)

    integer(4),               intent(in)    :: count4
    MY_MPI_REQUEST,           intent(inout) :: ireq4(:)
    integer(ip),    optional, intent(inout) :: ipass
    integer(4)                              :: istat4
    MY_MPI_STATUS,  allocatable             :: status4(:)
    
#ifndef MPI_OFF
    if( IPARALL .and. count4 > 0 ) then
       allocate(status4(MPI_STATUS_SIZE*count4))
       call MPI_WAITALL(count4,ireq4,status4,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_WAITALL')
       deallocate(status4)
       if( present(ipass) ) ipass = 0
    end if
#endif
    
  end subroutine PAR_WAITALL
 
end module mod_communications_tools
!> @}
