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
!> @file    mod_communications_point_to_point.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!> @brief   ToolBox for parallel communications
!> @details ToolBox for parallel communications
!------------------------------------------------------------------------

module mod_communications_point_to_point
  use def_kintyp_dims,                         only : npoin
  use def_kintyp_dims,                         only : npoin_2
  use def_kintyp_dims,                         only : nelem_2
  use def_kintyp_dims,                         only : nboun_2
  use def_kintyp_dims,                         only : npoi1
  use mod_maths_sort,                          only : maths_heap_sort
  use def_parall,                              only : kfl_order_exchange_par
  use def_kintyp_comm,                         only : COMM_ONE_SIDED
  use def_kintyp_comm,                         only : COMM_SEND_RECV
  use def_communications 
  use mod_communications_tools
  use mod_communications_global
  use mod_communications_point_to_point_basic
  use def_mpi
#include "def_mpi.inc"
  implicit none
  
  private
  
  integer(ip),            parameter   :: ON_NODES = 1
  integer(ip),            parameter   :: ON_EDGES = 2
  
  logical(lg),            allocatable :: tmp_lsend(:)
  logical(lg),            allocatable :: tmp_lrecv(:)
  integer(ip),            allocatable :: tmp_isend(:)
  integer(ip),            allocatable :: tmp_irecv(:)
  real(rp),               allocatable :: tmp_rsend(:)
  real(rp),               allocatable :: tmp_rrecv(:)  
  MY_MPI_REQUEST,         allocatable :: ireq4(:)
  real(rp),               pointer     :: sendbuff_rp(:)
  real(rp),               pointer     :: recvbuff_rp(:)
  type non_blocking_typ
     MY_MPI_REQUEST,     pointer      :: request4(:)
     integer(4)                       :: count4
  end type non_blocking_typ
  type(non_blocking_typ), pointer     :: non_blocking(:)
  integer(ip)                         :: inonblocking
  !
  ! Symmetric exchange
  ! 
  interface PAR_SYMMETRIC_EXCHANGE
     module procedure PAR_SYMMETRIC_EXCHANGE_RP_1,&
          &           PAR_SYMMETRIC_EXCHANGE_RP_2
  end interface PAR_SYMMETRIC_EXCHANGE
  !
  ! Exchange between slaves:
  ! type:  1. nodes  2. faces    3. fringe nodes
  ! kind:  1. real   2. integer  3. Complex
  ! dim:   1. x(:)   2. x(:,:)   3. x(:,:,:)   4. n,x(n) 
  ! what:  1. SUM    2. MAX      3. MIN
  ! where: in my code, in my zone
  !
  interface PAR_INTERFACE_NODE_EXCHANGE
     module procedure PAR_INTERFACE_NODE_EXCHANGE_IP_0,  &
          &           PAR_INTERFACE_NODE_EXCHANGE_IP_00, &
          &           PAR_INTERFACE_NODE_EXCHANGE_IP_1,  &
          &           PAR_INTERFACE_NODE_EXCHANGE_IP_2,  &
          &           PAR_INTERFACE_NODE_EXCHANGE_IP_3,  &
          &           PAR_INTERFACE_NODE_EXCHANGE_IP_2b, &
          &           PAR_INTERFACE_NODE_EXCHANGE_RP_0,  &
          &           PAR_INTERFACE_NODE_EXCHANGE_RP_00, &
          &           PAR_INTERFACE_NODE_EXCHANGE_RP_1,  &
          &           PAR_INTERFACE_NODE_EXCHANGE_RP_2,  &
          &           PAR_INTERFACE_NODE_EXCHANGE_RP_3,  &
          &           PAR_INTERFACE_NODE_EXCHANGE_RP_2b, &
          &           PAR_INTERFACE_NODE_EXCHANGE_LG_1
  end interface PAR_INTERFACE_NODE_EXCHANGE
  interface PAR_INTERFACE_EDGE_EXCHANGE
     module procedure PAR_INTERFACE_EDGE_EXCHANGE_IP_0,  &
          &           PAR_INTERFACE_EDGE_EXCHANGE_IP_1,  &
          &           PAR_INTERFACE_EDGE_EXCHANGE_IP_2,  &
          &           PAR_INTERFACE_EDGE_EXCHANGE_IP_3,  &
          &           PAR_INTERFACE_EDGE_EXCHANGE_IP_2b, &
          &           PAR_INTERFACE_EDGE_EXCHANGE_RP_00, &
          &           PAR_INTERFACE_EDGE_EXCHANGE_RP_0,  &
          &           PAR_INTERFACE_EDGE_EXCHANGE_RP_1,  &
          &           PAR_INTERFACE_EDGE_EXCHANGE_RP_2,  &
          &           PAR_INTERFACE_EDGE_EXCHANGE_RP_3,  &
          &           PAR_INTERFACE_EDGE_EXCHANGE_RP_2b, &
          &           PAR_INTERFACE_EDGE_EXCHANGE_LG_1
  end interface PAR_INTERFACE_EDGE_EXCHANGE
  !
  ! Exchange on own nodes
  !
  interface PAR_INTERFACE_OWN_NODE_EXCHANGE
     module procedure PAR_INTERFACE_OWN_NODE_EXCHANGE_IP_0,  &
          &           PAR_INTERFACE_OWN_NODE_EXCHANGE_IP_00, &
          &           PAR_INTERFACE_OWN_NODE_EXCHANGE_IP_1,  &
          &           PAR_INTERFACE_OWN_NODE_EXCHANGE_IP_2,  &
          &           PAR_INTERFACE_OWN_NODE_EXCHANGE_RP_0,  &
          &           PAR_INTERFACE_OWN_NODE_EXCHANGE_RP_00, &
          &           PAR_INTERFACE_OWN_NODE_EXCHANGE_RP_1,  &
          &           PAR_INTERFACE_OWN_NODE_EXCHANGE_RP_2
  end interface PAR_INTERFACE_OWN_NODE_EXCHANGE
  !
  ! Interface node matrix exchange
  !
  interface PAR_INTERFACE_MATRIX_EXCHANGE
     module procedure PAR_INTERFACE_MATRIX_EXCHANGE_WHERE,&
          &           PAR_INTERFACE_MATRIX_EXCHANGE_COMM, &
          &           PAR_INTERFACE_MATRIX_EXCHANGE_COMM_SHAPE
  end interface PAR_INTERFACE_MATRIX_EXCHANGE
  !
  ! Interface node matrix with halos exchange
  !
  interface PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE
     module procedure PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP_1,&
          &           PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP_2
  end interface PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE
  !
  ! Exchange between assymetric
  !
  interface PAR_COUPLING_NODE_EXCHANGE
     module procedure PAR_COUPLING_NODE_EXCHANGE_RP_0,  &
          &           PAR_COUPLING_NODE_EXCHANGE_RP_1,  &
          &           PAR_COUPLING_NODE_EXCHANGE_RP_1b, &
          &           PAR_COUPLING_NODE_EXCHANGE_RP_2,  &
          &           PAR_COUPLING_NODE_EXCHANGE_RP_3,  &
          &           PAR_COUPLING_NODE_EXCHANGE_RP_2b
  end interface PAR_COUPLING_NODE_EXCHANGE
  !
  ! Exchange between slaves on ghost elements/boundaries/nodes
  !
  ! +----+----+----+....o
  !             ||   /|
  !             \/   ||
  !           o....+----+----+----+
  !
  !
  interface PAR_GHOST_ELEMENT_EXCHANGE
     module procedure PAR_GHOST_ELEMENT_EXCHANGE_IP_0,  &
          &           PAR_GHOST_ELEMENT_EXCHANGE_IP_1,  &
          &           PAR_GHOST_ELEMENT_EXCHANGE_IP_2,  &
          &           PAR_GHOST_ELEMENT_EXCHANGE_IP_3,  &
          &           PAR_GHOST_ELEMENT_EXCHANGE_IP_2b, &
          &           PAR_GHOST_ELEMENT_EXCHANGE_RP_00, &
          &           PAR_GHOST_ELEMENT_EXCHANGE_RP_0,  &
          &           PAR_GHOST_ELEMENT_EXCHANGE_RP_1,  &
          &           PAR_GHOST_ELEMENT_EXCHANGE_RP_2,  &
          &           PAR_GHOST_ELEMENT_EXCHANGE_RP_3,  &
          &           PAR_GHOST_ELEMENT_EXCHANGE_RP_2b
  end interface PAR_GHOST_ELEMENT_EXCHANGE
  interface PAR_GHOST_BOUNDARY_EXCHANGE
     module procedure PAR_GHOST_BOUNDARY_EXCHANGE_IP_0,  &
          &           PAR_GHOST_BOUNDARY_EXCHANGE_IP_1,  &
          &           PAR_GHOST_BOUNDARY_EXCHANGE_IP_2,  &
          &           PAR_GHOST_BOUNDARY_EXCHANGE_IP_3,  &
          &           PAR_GHOST_BOUNDARY_EXCHANGE_IP_2b, &
          &           PAR_GHOST_BOUNDARY_EXCHANGE_RP_00, &
          &           PAR_GHOST_BOUNDARY_EXCHANGE_RP_0,  &
          &           PAR_GHOST_BOUNDARY_EXCHANGE_RP_1,  &
          &           PAR_GHOST_BOUNDARY_EXCHANGE_RP_2,  &
          &           PAR_GHOST_BOUNDARY_EXCHANGE_RP_3,  &
          &           PAR_GHOST_BOUNDARY_EXCHANGE_RP_2b
  end interface PAR_GHOST_BOUNDARY_EXCHANGE
  interface PAR_GHOST_NODE_EXCHANGE
     module procedure PAR_GHOST_NODE_EXCHANGE_IP_0,  &
          &           PAR_GHOST_NODE_EXCHANGE_IP_1,  &
          &           PAR_GHOST_NODE_EXCHANGE_IP_2,  &
          &           PAR_GHOST_NODE_EXCHANGE_IP_3,  &
          &           PAR_GHOST_NODE_EXCHANGE_IP_2b, &
          &           PAR_GHOST_NODE_EXCHANGE_RP_00, &
          &           PAR_GHOST_NODE_EXCHANGE_RP_0,  &
          &           PAR_GHOST_NODE_EXCHANGE_RP_1,  &
          &           PAR_GHOST_NODE_EXCHANGE_RP_2,  &
          &           PAR_GHOST_NODE_EXCHANGE_RP_3,  &
          &           PAR_GHOST_NODE_EXCHANGE_RP_2b
  end interface PAR_GHOST_NODE_EXCHANGE
  !
  ! Exchange between slaves on ghost elements/boundaries/nodes
  !
  ! +----+----+----+....o
  !             /|   ||
  !             ||   \/
  !           o....+----+----+----+
  !
  !
  interface PAR_FROM_GHOST_ELEMENT_EXCHANGE
     module procedure PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_00, &
          &           PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_0,  &
          &           PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_1,  &
          &           PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_2,  &
          &           PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_3,  &
          &           PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_2b, &
          &           PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_00, &
          &           PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_0,  &
          &           PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_1,  &
          &           PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_2,  &
          &           PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_3,  &
          &           PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_2b
  end interface PAR_FROM_GHOST_ELEMENT_EXCHANGE
  interface PAR_FROM_GHOST_BOUNDARY_EXCHANGE
     module procedure PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_00, &
          &           PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_0,  &
          &           PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_1,  &
          &           PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_2,  &
          &           PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_3,  &
          &           PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_2b, &
          &           PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_00, &
          &           PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_0,  &
          &           PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_1,  &
          &           PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_2,  &
          &           PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_3,  &
          &           PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_2b
  end interface PAR_FROM_GHOST_BOUNDARY_EXCHANGE
  interface PAR_FROM_GHOST_NODE_EXCHANGE
     module procedure PAR_FROM_GHOST_NODE_EXCHANGE_IP_00, &
          &           PAR_FROM_GHOST_NODE_EXCHANGE_IP_0,  &
          &           PAR_FROM_GHOST_NODE_EXCHANGE_IP_1,  &
          &           PAR_FROM_GHOST_NODE_EXCHANGE_IP_2,  &
          &           PAR_FROM_GHOST_NODE_EXCHANGE_IP_3,  &
          &           PAR_FROM_GHOST_NODE_EXCHANGE_IP_2b, &
          &           PAR_FROM_GHOST_NODE_EXCHANGE_RP_1,  &
          &           PAR_FROM_GHOST_NODE_EXCHANGE_RP_2          
  end interface PAR_FROM_GHOST_NODE_EXCHANGE
  !
  ! Send
  !
  interface PAR_SEND
     module procedure PAR_SEND_IP_s,&
          &           PAR_SEND_IP_0,&
          &           PAR_SEND_IP_1,&
          &           PAR_SEND_IP_2,&
          &           PAR_SEND_RP_s,&
          &           PAR_SEND_RP_0,&
          &           PAR_SEND_RP_1,&
          &           PAR_SEND_RP_2
  end interface PAR_SEND
  !
  ! Receive
  !
  interface PAR_RECEIVE
     module procedure PAR_RECEIVE_IP_s,&
          &           PAR_RECEIVE_IP_0,&
          &           PAR_RECEIVE_IP_1,&
          &           PAR_RECEIVE_IP_2,&
          &           PAR_RECEIVE_RP_s,&
          &           PAR_RECEIVE_RP_0,&
          &           PAR_RECEIVE_RP_1,&
          &           PAR_RECEIVE_RP_2
  end interface PAR_RECEIVE
  !
  ! Send/receive
  !
  interface PAR_SEND_RECEIVE
     module procedure PAR_SEND_RECEIVE_IP_s,&
          &           PAR_SEND_RECEIVE_IP_0,&
          &           PAR_SEND_RECEIVE_IP_1,&
          &           PAR_SEND_RECEIVE_IP_2,&
          &           PAR_SEND_RECEIVE_IP_3,&
          &           PAR_SEND_RECEIVE_RP_s,&
          &           PAR_SEND_RECEIVE_RP_0,&
          &           PAR_SEND_RECEIVE_RP_1,&
          &           PAR_SEND_RECEIVE_RP_2,&
          &           PAR_SEND_RECEIVE_RP_3
  end interface PAR_SEND_RECEIVE
  interface PAR_SEND_RECEIVE_TO_ALL
     module procedure PAR_SEND_RECEIVE_TO_ALL_RP_s,&
          &           PAR_SEND_RECEIVE_TO_ALL_RP_0,&
          &           PAR_SEND_RECEIVE_TO_ALL_RP_1,&
          &           PAR_SEND_RECEIVE_TO_ALL_RP_1c,&
          &           PAR_SEND_RECEIVE_TO_ALL_RP_2,&
          &           PAR_SEND_RECEIVE_TO_ALL_RP_3,&
          &           PAR_SEND_RECEIVE_TO_ALL_RP_4,&
          !&           PAR_SEND_RECEIVE_TO_ALL_IP_1,&
          &           PAR_SEND_RECEIVE_TO_ALL_IP_1c,&
          &           PAR_SEND_RECEIVE_TO_ALL_IP_2,&
          &           PAR_SEND_RECEIVE_TO_ALL_IP_0
  end interface PAR_SEND_RECEIVE_TO_ALL
  !
  ! Operations with neighbors
  !
  interface PAR_POINT_TO_POINT_ARRAY_OPERATION
     module procedure PAR_POINT_TO_POINT_ARRAY_OPERATION_IP_1
  end interface PAR_POINT_TO_POINT_ARRAY_OPERATION
  !
  ! Bridge to broadcast data from master to slaves
  !
  interface PAR_EXCHANGE
     module procedure PAR_EXCHANGE_IP_s,PAR_EXCHANGE_IP_0,PAR_EXCHANGE_IP_02,PAR_EXCHANGE_IP_1,PAR_EXCHANGE_IP_2,&
          &           PAR_EXCHANGE_RP_s,PAR_EXCHANGE_RP_0,PAR_EXCHANGE_RP_02,PAR_EXCHANGE_RP_1,PAR_EXCHANGE_RP_2,&
          &           PAR_EXCHANGE_CH,PAR_EXCHANGE_CH_1,&
          &           PAR_EXCHANGE_LG_s,PAR_EXCHANGE_LG_0,PAR_EXCHANGE_LG_1
  end interface PAR_EXCHANGE
  !
  ! Start non-blocking communications
  !
  interface PAR_START_NON_BLOCKING_COMM
     module procedure PAR_START_NON_BLOCKING_COMM_4,&
          &           PAR_START_NON_BLOCKING_COMM_4_OPT,&
          &           PAR_START_NON_BLOCKING_COMM_8
  end interface PAR_START_NON_BLOCKING_COMM
  !
  ! Public functions
  !
  public :: PAR_SYMMETRIC_EXCHANGE             ! Symmetric exchange 
  public :: PAR_INTERFACE_NODE_EXCHANGE        ! Interface nodal exchange (send/receive)
  public :: PAR_INTERFACE_OWN_NODE_EXCHANGE    ! Interface OWN node exchange (send/receive)
  public :: PAR_INTERFACE_EDGE_EXCHANGE        ! Interface edge exchange (send/receive)
  public :: PAR_INTERFACE_MATRIX_EXCHANGE      ! Matrix exchange on interface nodes
  public :: PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE ! Matrix with halos exchange on interface nodes
  public :: PAR_COUPLING_NODE_EXCHANGE         ! Coupling nodal exchange (send/receive)
  public :: PAR_GHOST_ELEMENT_EXCHANGE         ! Ghost element exchange (send/receive): from element to ghost
  public :: PAR_GHOST_NODE_EXCHANGE            ! Ghost nodal exchange (send/receive)
  public :: PAR_FROM_GHOST_ELEMENT_EXCHANGE    ! Ghost element exchange (send/receive): from ghost to element
  public :: PAR_FROM_GHOST_BOUNDARY_EXCHANGE   ! Ghost element exchange (send/receive): from ghost to boundary
  public :: PAR_FROM_GHOST_NODE_EXCHANGE       ! Ghost node exchange (send/receive): from ghost to node
  public :: PAR_GHOST_BOUNDARY_EXCHANGE        ! Ghost boundary exchange (send/receive)
  public :: PAR_SEND                           ! Send arrays to a specific partition
  public :: PAR_RECEIVE                        ! Receive arrays to a specific partition
  public :: PAR_SEND_RECEIVE                   ! Send and receive arrays to a specific partition
  public :: PAR_SEND_RECEIVE_IP                ! Send and receive arrays to a specific partition
  public :: PAR_SEND_RECEIVE_RP                ! Send and receive arrays to a specific partition
  public :: PAR_SEND_RECEIVE_RP_2              ! Send and receive arrays to a specific partition
  public :: PAR_SEND_RECEIVE_TO_ALL            ! Send and receive arrays to all
  public :: PAR_EXCHANGE                       ! Bridge to broadcast data from master to slave
  public :: PAR_INITIALIZE_NON_BLOCKING_COMM   ! Initialize non-blocking communications variables
  public :: PAR_START_NON_BLOCKING_COMM        ! Start non-blocking communications
  public :: PAR_END_NON_BLOCKING_COMM          ! End non-blocking communications
  public :: PAR_SET_NON_BLOCKING_COMM_NUMBER   ! Set the non-blocking communicator number
  public :: PAR_WAITALL                        ! Waitall
  public :: PAR_POINT_TO_POINT_ARRAY_OPERATION ! Array operations between all partitions of communicators
  public :: PAR_COMM_NULL                      ! NULL MPI communicator
  public :: PAR_MIN_RANK_OWNER                 ! Get the minimum owner rank of an interface node
  public :: PAR_INTERFACE_EXCHANGE             ! Generic exchange
  
  public :: communications_point_to_point_initialization

  public :: PAR_INTERFACE_NODE_EXCHANGE_RP_1
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-13
  !> @brief   Initialization
  !> @details Initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine communications_point_to_point_initialization()
    !
    ! Nullify and initialize
    !
    nullify(non_blocking)
    nullify(sendbuff_rp)
    nullify(recvbuff_rp)
    !
    ! Initialize non-blocking communcation arrays
    !
    call PAR_INITIALIZE_NON_BLOCKING_COMM()
     
  end subroutine communications_point_to_point_initialization
  
  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_INTERFACE_NODE_EXCHANGE_IP
  !
  !----------------------------------------------------------------------

  subroutine PAR_INTERFACE_NODE_EXCHANGE_IP_0(n,xx,what,wherein,wsynch,dom_k,EXCLUDE_MYSELF)
    
    integer(ip),                    intent(in)    :: n
    integer(ip),                    intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),         optional, intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    logical(lg),          optional, intent(in)    :: EXCLUDE_MYSELF
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE .or. npoin == 0 ) return
    ndofn = n
    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
    end if
    call PAR_INTERFACE_NODE_EXCHANGE_IP(ON_NODES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k,EXCLUDE_MYSELF)

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_IP_0

  subroutine PAR_INTERFACE_NODE_EXCHANGE_IP_00(n,xx,what,wherein,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: n
    integer(ip),                    intent(inout) :: xx(*)
    character(*),                   intent(in)    :: what
    character(*),         optional, intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
    end if
    call PAR_INTERFACE_NODE_EXCHANGE_IP(ON_NODES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_IP_00

  subroutine PAR_INTERFACE_NODE_EXCHANGE_IP_1(xx,what,wherein,wsynch,dom_k,COMM,EXCLUDE_MYSELF)
    
    integer(ip),          pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),         optional, intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    type(comm_data_par),  optional, intent(in)    :: COMM
    logical(lg),          optional, intent(in)    :: EXCLUDE_MYSELF
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

   if( INOTSLAVE .or. npoin == 0 .or. .not. associated(xx) ) return
    ndofn = 1
    if( present(COMM) ) then
       PAR_COMM_TO_USE = COMM % PAR_COMM_WORLD
       call PAR_INTERFACE_NODE_EXCHANGE_IP(ON_NODES,ndofn,xx,what,COMM,PAR_COMM_TO_USE,wsynch,dom_k,EXCLUDE_MYSELF)
    else
       !if( size(xx,1) /= npoin .and. size(xx,1) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE IP_1')
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
       end if
       call PAR_INTERFACE_NODE_EXCHANGE_IP(ON_NODES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k,EXCLUDE_MYSELF)
    end if

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_IP_1

  subroutine PAR_INTERFACE_NODE_EXCHANGE_IP_2(xx,what,wherein,wsynch,dom_k,EXCLUDE_MYSELF)
    
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),         optional, intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    logical(lg),          optional, intent(in)    :: EXCLUDE_MYSELF
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE .or. npoin == 0 ) return
    ndofn = size(xx,1)
    if( size(xx,2) /= npoin .and. size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE IP_2')
    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
    end if
    call PAR_INTERFACE_NODE_EXCHANGE_IP(ON_NODES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k,EXCLUDE_MYSELF)

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_IP_2

  subroutine PAR_INTERFACE_NODE_EXCHANGE_IP_3(xx,what,wherein,wsynch,dom_k,EXCLUDE_MYSELF)
    
    integer(ip),          pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),         optional, intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    logical(lg),          optional, intent(in)    :: EXCLUDE_MYSELF
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE .or. npoin == 0 ) return
    if( size(xx,3) <= 2 ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= npoin .and. size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE IP_3')
    else
       ndofn = size(xx,1)*size(xx,2)
       if( size(xx,3) /= npoin .and. size(xx,3) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE IP_3')
    end if
    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
    end if
    call PAR_INTERFACE_NODE_EXCHANGE_IP(ON_NODES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k,EXCLUDE_MYSELF)

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_IP_3

  subroutine PAR_INTERFACE_NODE_EXCHANGE_IP_2b(xx,what,commu,wsynch,dom_k,EXCLUDE_MYSELF)
    
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),  pointer,  intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    logical(lg),          optional, intent(in)    :: EXCLUDE_MYSELF
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE .or. npoin == 0 ) return
    ndofn = size(xx,1)
    if( size(xx,2) /= npoin .and. size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
    PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
    call PAR_INTERFACE_NODE_EXCHANGE_IP(ON_NODES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k,EXCLUDE_MYSELF)

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_IP_2b

  !----------------------------------------------------------------------
  !
  ! PAR_INTERFACE_NODE_EXCHANGE_IP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_INTERFACE_NODE_EXCHANGE_IP(onwhat,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k,EXCLUDE_MYSELF)
    
    integer(ip),                    intent(in)    :: onwhat
    integer(ip),                    intent(in)    :: ndofn
    integer(ip),                    intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    MY_MPI_COMM,                    intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    logical(lg),          optional, intent(in)    :: EXCLUDE_MYSELF
    integer(ip)                                   :: ii,nsize,jj,dom_i,kdofn
    integer(ip)                                   :: ipoin,ini,kk,idofn,jdofn
    integer(4)                                    :: istat4,nsize4,count4
    integer(4)                                    :: dom_i4,my_rank4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    MY_MPI_STATUS,                     pointer    :: status4(:)
    integer(ip)                                   :: dom_j
    integer(ip)                                   :: bound_dim
    integer(ip),                       pointer    :: bound_perm(:)
    integer(ip),                       pointer    :: bound_invp(:)
    integer(ip),                       pointer    :: bound_size(:)
    logical(lg)                                   :: if_exclude_myself
    
#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! On nodes or edges
       !
       if( onwhat == ON_NODES ) then
          bound_dim  =  commu % bound_dim
          bound_perm => commu % bound_perm
          bound_invp => commu % bound_invp
          bound_size => commu % bound_size
       else
          bound_dim  =  commu % bedge_dim
          bound_perm => commu % bedge_perm
          bound_invp => commu % bedge_perm
          bound_size => commu % bedge_size
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_isend(bound_dim * ndofn))
             allocate(tmp_irecv(bound_dim * ndofn))
             !
             ! Save in temp_send
             !
             if( trim(what) == 'DISTRIBUTE' ) then
                do jj = 1,bound_dim
                   ipoin = bound_perm(jj)
                   if( ipoin < commu % npoi2 .or. ipoin > commu % npoi3 ) xx(1:ndofn,ipoin) = 0_ip
                end do
             end if

             kk = 0
             do jj = 1,bound_dim
                ipoin = bound_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_isend(kk) = xx(idofn,ipoin)
                   tmp_irecv(kk) = 0
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             istat4 = 0_4
             kk = 0
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   ini   = ndofn * ( bound_size(ii)   - 1 ) + 1
                   nsize = ndofn * ( bound_size(ii+1) - 1 ) + 1 - ini

                   nsize4 = int(nsize,4)
                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_isend(ini:ini+nsize-1), nsize4, &
                           PAR_INTEGER,  dom_i4, 0_4,          &
                           PAR_COMM_TO_USE, ireq4(kk), istat4 )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_irecv(ini:ini+nsize-1), nsize4, &
                           PAR_INTEGER,  dom_i4, 0_4,          &
                           PAR_COMM_TO_USE, ireq4(kk), istat4 )
                   else
                      call MPI_Sendrecv(                       &
                           tmp_isend(ini:), nsize4,            &
                           PAR_INTEGER, dom_i4, 0_4,           &
                           tmp_irecv(ini:), nsize4,            &
                           PAR_INTEGER, dom_i4, 0_4,           &
                           PAR_COMM_TO_USE, status, istat4    )
                   end if
                   if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_INTERFACE_NODE_EXCHANGE_IP')
                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE*2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' .or. trim(what) == 'DISTRIBUTE' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_invp(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_irecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_invp(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_invp(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'TAKE MIN' ) then
                !
                ! TAKE MIN
                !
                if( present(EXCLUDE_MYSELF) ) then
                   if_exclude_myself = EXCLUDE_MYSELF
                else
                   if_exclude_myself = .false.
                end if

                my_rank4 = commd % RANK4
                kk = 0
                do ii = 1,commd % bound_dim
                   ipoin = bound_invp(ii)
                   do idofn = 1,ndofn
                      xx(idofn,ipoin) = 0
                   end do
                end do
                
                if( if_exclude_myself ) then
                   do ii = 1,commu % nneig
                      dom_i = commu % neights(ii)
                      do jj = bound_size(ii),bound_size(ii+1)-1
                         ipoin = bound_invp(jj)
                         do idofn = 1,ndofn
                            kk = kk + 1
                            if( dom_i /= my_rank4 ) then
                               if( my_rank4 <= dom_i ) then
                                  xx(idofn,ipoin) = tmp_irecv(kk)
                              ! else
                              !    xx(idofn,ipoin) = 0
                               end if
                            !else
                            !   xx(idofn,ipoin) = 0
                            end if
                         end do
                      end do
                   end do
                else
                   do ii = 1,commu % nneig
                      dom_i = commu % neights(ii)
                      do jj = bound_size(ii),bound_size(ii+1)-1
                         ipoin = bound_invp(jj)
                         do idofn = 1,ndofn
                            kk = kk + 1
                            if( my_rank4 <= dom_i ) then
                               xx(idofn,ipoin) = tmp_irecv(kk)
                            else
                               xx(idofn,ipoin) = 0
                            end if
                         end do
                      end do
                   end do
                end if
                
             else if( trim(what) == 'MIN RANK OR NEGATIVE' ) then
                !
                ! MIN RANK OR NEGATIVE
                !
                call runend('MIN RANK OR NEGATIVE NOT CODED')
                !call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank4)
                my_rank4 = commu % RANK4
                kk = 0
                do ii = 1,commu % nneig
                   dom_i = commu % neights(ii)
                   ini   = ndofn * ( bound_size(ii)   - 1 ) + 1
                   nsize = ndofn * ( bound_size(ii+1) - 1 ) + 1 - ini
                   do jj = bound_size(ii),bound_size(ii+1)-1
                      ipoin = bound_invp(jj)
                      do idofn = 1,ndofn
                         kk = kk + 1
                         if( xx(idofn,ipoin) > 0 .and. tmp_irecv(kk) > 0 ) then
                            if( my_rank4 > dom_i ) then
                               !pard1 = 1
                               xx(idofn,ipoin) = -abs(xx(idofn,ipoin))
                            end if
                         end if
                      end do
                   end do
                end do

             else if( trim(what) == 'DIF' .or. trim(what) == 'DIFFERENCE' .or. trim(what) == 'DIFF' ) then
                !
                ! DIF
                !
                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_invp(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      if( abs(xx(idofn,ipoin)-tmp_irecv(kk)) > 0 ) then
                         if( ipoin <= npoi1 .or. ( ipoin >= commu % npoi2 .and. ipoin <= commu % npoi3 ) ) write(6,*) 'DIF=',xx(idofn,ipoin)-tmp_irecv(kk)
                         xx(idofn,ipoin) = huge(1_ip)
                      end if
                   end do
                end do

             else if( trim(what) == 'MERGE' ) then
                !
                ! MERGE
                !
                kk = 0

                do jj = 1,bound_dim
                   ipoin = bound_invp(jj)

                   loop_kdofn1: do kdofn = 1,ndofn
                      if( xx(kdofn,ipoin) == 0 ) exit loop_kdofn1
                   end do loop_kdofn1
                   kdofn = min(kdofn,ndofn)

                   do idofn = 1,ndofn
                      kk = kk + 1
                      jdofn = 0
                      do while( jdofn < ndofn )
                         jdofn = jdofn + 1
                         if( tmp_irecv(kk) == xx(jdofn,ipoin) ) jdofn = 2*ndofn
                      end do
                      if( jdofn /= 2*ndofn ) then
                         xx(kdofn,ipoin) = tmp_irecv(kk)
                      end if
                   end do
                end do

                do jj = 1,bound_dim
                   ipoin = bound_invp(jj)
                   loop_kdofn2: do kdofn = 1,ndofn
                      if( xx(kdofn,ipoin) == 0 ) exit loop_kdofn2
                   end do loop_kdofn2
                   kdofn = min(kdofn,ndofn)
                   if( kdofn > 0 ) call maths_heap_sort(1_ip,kdofn,xx(1:kdofn,ipoin))
                end do

             else
                call runend('UNKNOWN ORDER: '//trim(what))
             end if

             ipass = 0
             deallocate(tmp_irecv)
             deallocate(tmp_isend)

          end if

       end if

    end if
#endif

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_IP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_INTERFACE_NODE_EXCHANGE_LG
  !
  !----------------------------------------------------------------------

  subroutine PAR_INTERFACE_NODE_EXCHANGE_LG_1(xx,what,wherein,wsynch,dom_k)
    
    logical(lg),          pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),         optional, intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

   if( INOTSLAVE ) return
    ndofn = 1
    if( size(xx,1) /= npoin .and. size(xx,1) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
    end if
    call PAR_INTERFACE_NODE_EXCHANGE_LG(ON_NODES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_LG_1

  !----------------------------------------------------------------------
  !
  ! PAR_INTERFACE_NODE_EXCHANGE_IP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_INTERFACE_NODE_EXCHANGE_LG(onwhat,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: onwhat
    integer(ip),                    intent(in)    :: ndofn
    logical(lg),                    intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    MY_MPI_COMM   ,                   intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,nsize,jj,dom_i
    integer(ip)                                   :: ipoin,ini,kk,idofn
    integer(4)                                    :: istat4,nsize4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    MY_MPI_STATUS ,                    pointer    :: status4(:)
    integer(ip)                                   :: dom_j
    integer(ip)                                   :: bound_dim
    integer(ip),                       pointer    :: bound_perm(:)
    integer(ip),                       pointer    :: bound_size(:)

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! On nodes or edges
       !
       if( onwhat == ON_NODES ) then
          bound_dim  =  commu % bound_dim
          bound_perm => commu % bound_perm
          bound_size => commu % bound_size
       else
          bound_dim  =  commu % bedge_dim
          bound_perm => commu % bedge_perm
          bound_size => commu % bedge_size
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_lsend(bound_dim * ndofn))
             allocate(tmp_lrecv(bound_dim * ndofn))
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,bound_dim
                ipoin = bound_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_lsend(kk) = xx(idofn,ipoin)
                   tmp_lrecv(kk) = .false.
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             istat4 = 0_4
             kk = 0
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   ini   = ndofn * ( bound_size(ii)   - 1 ) + 1
                   nsize = ndofn * ( bound_size(ii+1) - 1 ) + 1 - ini

                   nsize4 = int(nsize,4)
                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_lsend(ini:ini+nsize-1), nsize4, &
                           MPI_LOGICAL,  dom_i4, 0_4,          &
                           PAR_COMM_TO_USE, ireq4(kk), istat4 )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_lrecv(ini:ini+nsize-1), nsize4, &
                           MPI_LOGICAL,  dom_i4, 0_4,          &
                           PAR_COMM_TO_USE, ireq4(kk), istat4 )
                   else
                      call MPI_Sendrecv(                       &
                           tmp_lsend(ini:), nsize4,            &
                           MPI_LOGICAL, dom_i4, 0_4,           &
                           tmp_lrecv(ini:), nsize4,            &
                           PAR_INTEGER, dom_i4, 0_4,           &
                           PAR_COMM_TO_USE, status, istat4    )
                   end if
                   if( istat4 /= MPI_SUCCESS )call PAR_MPI_RUNEND(istat4,'PAR_INTERFACE_NODE_EXCHANGE_LG')
                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE*2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'OR' ) then
                !
                ! OR
                !
                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) .or. tmp_lrecv(kk)
                   end do
                end do

             else if( trim(what) == 'AND' ) then
                !
                ! AND
                !
                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) .and. tmp_lrecv(kk)
                   end do
                end do

             else
                call runend('UNKNOWN ORDER: '//trim(what))
             end if

             ipass = 0
             deallocate(tmp_lrecv)
             deallocate(tmp_lsend)

          end if

       end if

    end if
#endif

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_LG

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_INTERFACE_NODE_EXCHANGE_RP
  !
  !----------------------------------------------------------------------

  subroutine PAR_INTERFACE_NODE_EXCHANGE_RP_00(n,xx,what,wherein,wsynch,dom_k,COMM)
    
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(*)
    character(*),                   intent(in)    :: what
    character(*),         optional, intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    type(comm_data_par),  optional, intent(in)    :: COMM
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    if( present(COMM) ) then
       PAR_COMM_TO_USE = COMM % PAR_COMM_WORLD
       call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_NODES,ndofn,xx,what,COMM,PAR_COMM_TO_USE,wsynch,dom_k)  
    else
       
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
       end if
       call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_NODES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if
    
  end subroutine PAR_INTERFACE_NODE_EXCHANGE_RP_00

  subroutine PAR_INTERFACE_NODE_EXCHANGE_RP_0(n,xx,what,wherein,wsynch,dom_k,COMM)
    
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),         optional, intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    type(comm_data_par),  optional, intent(in)    :: COMM
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    if( present(COMM) ) then
       PAR_COMM_TO_USE = COMM % PAR_COMM_WORLD
       call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_NODES,ndofn,xx,what,COMM,PAR_COMM_TO_USE,wsynch,dom_k)  
    else
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
       end if
       call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_NODES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if
    
  end subroutine PAR_INTERFACE_NODE_EXCHANGE_RP_0

  subroutine PAR_INTERFACE_NODE_EXCHANGE_RP_1(xx,what,wherein,wsynch,dom_k,COMM)
    
    real(rp),             pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),         optional, intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    type(comm_data_par),  optional, intent(in)    :: COMM
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE .or. npoin == 0 ) return
    ndofn = 1
    if( what(len(what):len(what)) /= 'I' ) then
       if( size(xx,1) /= npoin .and. size(xx,1) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE RP_1')
    end if

    if( present(COMM) ) then
       PAR_COMM_TO_USE = COMM % PAR_COMM_WORLD
       call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_NODES,ndofn,xx,what,COMM,PAR_COMM_TO_USE,wsynch,dom_k)  
    else
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
       end if
       call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_NODES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if
    
  end subroutine PAR_INTERFACE_NODE_EXCHANGE_RP_1

  subroutine PAR_INTERFACE_NODE_EXCHANGE_RP_2(xx,what,wherein,wsynch,dom_k,COMM)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),         optional, intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    type(comm_data_par),  optional, intent(in)    :: COMM
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE .or. npoin == 0 ) return
    ndofn = size(xx,1)
    if( what(len(what):len(what)) /= 'I' ) then
       if( size(xx,2) /= npoin .and. size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE RP_2')
    end if
    if( present(COMM) ) then
       PAR_COMM_TO_USE = COMM % PAR_COMM_WORLD
       call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_NODES,ndofn,xx,what,COMM,PAR_COMM_TO_USE,wsynch,dom_k)
    else
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
       end if
       call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_NODES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if
    
  end subroutine PAR_INTERFACE_NODE_EXCHANGE_RP_2

  subroutine PAR_INTERFACE_NODE_EXCHANGE_RP_3(xx,what,wherein,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),         optional, intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE .or. npoin == 0 ) return
    if( what(len(what):len(what)) /= 'I' ) then
       if( size(xx,3) <= 2 ) then
          ndofn = size(xx,1)
          if( size(xx,2) /= npoin .and. size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE IP_3')
       else
          ndofn = size(xx,1)*size(xx,2)
          if( size(xx,3) /= npoin .and. size(xx,3) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE RP_3')
       end if
    end if
    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
    end if
    call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_NODES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_RP_3

  subroutine PAR_INTERFACE_NODE_EXCHANGE_RP_2b(xx,what,commu,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE .or. npoin == 0 ) return
    ndofn = size(xx,1)
    if( what(len(what):len(what)) /= 'I' ) then
       if( size(xx,2) /= npoin .and. size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
    end if
    PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
    call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_NODES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_RP_2b

  !----------------------------------------------------------------------
  !
  ! PAR_INTERFACE_NODE_EXCHANGE_RP: NODE ASSEMBLY FOR INTEGERS
  ! If last letter of what is I (meaning I interface) then it is assumed
  ! that only the interface part of the vector is passed, from npoi1+1
  ! to npoin.
  !
  !              Interior nodes                    Interface nodes
  ! <----------------------------------------><----------------------->
  ! +-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
  ! |  1  |     |     |     |     |     |NPOI1|     |     |     |NPOIN|
  ! +-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
  !
  !----------------------------------------------------------------------

  subroutine PAR_INTERFACE_NODE_EXCHANGE_RP(onwhat,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

    integer(ip),                    intent(in)    :: onwhat
    integer(ip),                    intent(in)    :: ndofn
    real(rp),                       intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    class(comm_data_par_basic),     intent(in)    :: commu
    MY_MPI_COMM,                    intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,nsize,jj,dom_i
    integer(ip)                                   :: ipoin,ini,kk,idofn,jj2
    integer(ip)                                   :: offset,ipoin_offset
    integer(4)                                    :: istat4,nsize4,count4
    integer(4)                                    :: dom_i4,my_rank4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    MY_MPI_STATUS ,                    pointer    :: status4(:)
    integer(ip)                                   :: dom_j
    integer(ip)                                   :: bound_dim
    integer(ip),                       pointer    :: bound_perm(:)
    integer(ip),                       pointer    :: bound_invp(:)
    integer(ip),                       pointer    :: bound_size(:)
    real(rp),                          pointer    :: bound_diff(:,:)
    real(rp),         save,            pointer    :: my_rxx(:,:)
    integer(KIND=MPI_ADDRESS_KIND)                :: target_disp
    integer(4)                                    :: target_rank
    integer(4)                                    :: target_count
    integer(4)                                    :: origin_count
    integer(ip),  parameter                       :: put_get=1 ! 1: put, 2: get
    
#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1

       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! On nodes or edges
       !
       if( onwhat == ON_NODES ) then
          bound_dim  =  commu % bound_dim
          bound_perm => commu % bound_perm
          select type( commu )
          class is ( comm_data_par_basic )
             bound_invp => commu % bound_perm
          class is ( comm_data_par       )
             if( associated(commu % bound_invp) ) then
                bound_invp => commu % bound_invp
             else
                bound_invp => commu % bound_perm
             end if
          end select
          bound_size => commu % bound_size
       else
          select type( commu )
          class is ( comm_data_par )
             bound_dim  =  commu % bedge_dim
             bound_perm => commu % bedge_perm
             bound_invp => commu % bedge_perm
             bound_size => commu % bedge_size
          end select
       end if
       !
       ! Offset
       !
       if( what(len(what):len(what)) /= 'I' ) then
          offset = 0
       else
          offset = -npoi1
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( asynch .and. commu % type == COMM_SEND_RECV ) allocate(ireq4(commu % nneig*2))
             if( memory_size(sendbuff_rp) < bound_dim * ndofn ) then
                call memory_deallo(par_memor,'SENDBUFF_RP','par_interface_node_exchange_rp',sendbuff_rp)
                call memory_alloca(par_memor,'SENDBUFF_RP','par_interface_node_exchange_rp',sendbuff_rp,bound_dim * ndofn,'DO_NOT_INITIALIZE')
             end if
             if( memory_size(recvbuff_rp) < bound_dim * ndofn .and. commu % type == COMM_SEND_RECV ) then
                call memory_deallo(par_memor,'RECVBUFF_RP','par_interface_node_exchange_rp',recvbuff_rp)
                call memory_alloca(par_memor,'RECVBUFF_RP','par_interface_node_exchange_rp',recvbuff_rp,bound_dim * ndofn,'DO_NOT_INITIALIZE')
             else if( commu % type == COMM_ONE_SIDED ) then
                if( put_get == 1 ) then
                   recvbuff_rp => commd % buffe_recv
                else
                   recvbuff_rp => sendbuff_rp
                end if
             end if
             !
             ! Distribute values among neighbors
             !
             if( trim(what) == 'DISTRIBUTE' ) then
                do jj = 1,bound_dim
                   ipoin = bound_perm(jj)
                   if( ipoin < commu % npoi2 .or. ipoin > commu % npoi3 ) xx(1:ndofn,ipoin) = 0.0_rp
                end do
             end if
             kk = 0

             if( commu % type == COMM_ONE_SIDED .and. put_get == 2 ) then
               do jj = 1,bound_dim
                   ipoin = bound_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      commd % buffe_recv(kk) = xx(idofn,ipoin+offset) 
                   end do
                end do
             else
                do jj = 1,bound_dim
                   ipoin = bound_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      sendbuff_rp(kk) = xx(idofn,ipoin+offset)
                   end do
                end do
             end if
             
             if( commu % type == COMM_SEND_RECV ) then
                !
                ! Save in temp_send
                !
                if( kfl_order_exchange_par == 1 .and. commu % npoin-commu % npoi1 > 0 ) then
                   nullify(my_rxx)
                   allocate(my_rxx(ndofn,commu % npoin-commu % npoi1))
                   do ipoin = commu % npoi1+1,commu % npoin
                      my_rxx(1:ndofn,ipoin-commu % npoi1) = xx(1:ndofn,ipoin)
                   end do
                end if
                !
                ! Send    temp_send
                ! Receive temp_recv
                !
                istat4 = 0_4
                kk = 0
                do ii = 1,commu % nneig

                   dom_i  = commu % neights(ii)
                   dom_i4 = int(dom_i,4)

                   if( dom_j == 0 .or. dom_j == dom_i ) then

                      ini   = ndofn * ( bound_size(ii)   - 1 ) + 1
                      nsize = ndofn * ( bound_size(ii+1) - 1 ) + 1 - ini

                      nsize4 = int(nsize,4)
                      if( asynch ) then
                         kk = kk + 1
                         call MPI_Isend(&
                              sendbuff_rp(ini:ini+nsize-1), nsize4, &
                              PAR_REAL,  dom_i4, 0_4,               &
                              PAR_COMM_TO_USE, ireq4(kk), istat4 )
                         kk = kk + 1
                         call MPI_Irecv(&
                              recvbuff_rp(ini:ini+nsize-1), nsize4, &
                              PAR_REAL,  dom_i4, 0_4,               &
                              PAR_COMM_TO_USE, ireq4(kk), istat4 )
                      else
                         call MPI_Sendrecv(                         &
                              sendbuff_rp(ini:), nsize4,            &
                              PAR_REAL, dom_i4, 0_4,                &
                              recvbuff_rp(ini:), nsize4,            &
                              PAR_REAL, dom_i4, 0_4,                &
                              PAR_COMM_TO_USE, status, istat4 )
                      end if
                      if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_INTERFACE_NODE_EXCHANGE_RP')

                   end if

                end do
                
             else
                !
                ! One-sided
                !
                call MPI_Win_fence(0, commu % WIN_WM,istat4)
                !call MPI_Win_lock_all(0, commu %  WIN_WM,istat4)

                do ii = 1,commd % nneig                      
                   target_rank  = commu % neights(ii)-1
                   target_disp  = int(commu % displ_recv(ii)*ndofn,kind=MPI_ADDRESS_KIND)
                   target_count = (commu % bound_size(ii+1)-commu % bound_size(ii))*ndofn
                   kk           = (commu % bound_size(ii)-1)*ndofn+1
                   origin_count = target_count
                   if( put_get == 1 ) then
                      call MPI_Put(sendbuff_rp(kk:), ORIGIN_COUNT, PAR_REAL, TARGET_RANK,&
                           TARGET_DISP, TARGET_COUNT, PAR_REAL, commd % WIN_WM, istat4)
                   else
                      call MPI_Get(sendbuff_rp(kk:), ORIGIN_COUNT, PAR_REAL, TARGET_RANK,&
                           TARGET_DISP, TARGET_COUNT, PAR_REAL, commd % WIN_WM, istat4)
                   end if
                end do
                if( .not. asynch ) call MPI_Win_fence(0, commu % WIN_WM,istat4)
                !if( .not. asynch ) call MPI_Win_unlock_all(commu % WIN_WM,istat4)

             end if

          else if( asynch .and. ipass == 2 ) then
             !
             ! Second pass for one-sided and asynchronous
             !
             if( commu % type == COMM_SEND_RECV ) then
                
                count4 = 2*int(commu % nneig,4)
                allocate( status4(MPI_STATUS_SIZE*2*commu % nneig) )
                CALL MPI_WAITALL(count4,ireq4,status4,istat4)
                deallocate( status4 )
                deallocate(ireq4)
                
             else
                
                call MPI_Win_fence(0, commu % WIN_WM,istat4)
                !call MPI_Win_unlock_all(commu % WIN_WM,istat4)
                
             end if
             
          end if
          !
          ! sum,max,min on temp_recv
          !         
          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if(         trim(what) == 'SUM'  .or. trim(what) == 'ASSEMBLY' .or. trim(what) == 'SUMI'&
                  & .or. trim(what) == 'DISTRIBUTE' ) then
                !
                ! SUM
                !
                if( kfl_order_exchange_par == 1 .and. commu % npoin-commu % npoi1 > 0 ) then

                   do ipoin = commu % npoi1+1,commu % npoin
                      xx(1:ndofn,ipoin) = 0.0_rp 
                   end do

                   do ii = 1,commd % nneig_1
                      dom_i = commd % neights_ordered(ii)
                      kk    = commd % perm_ordered(ii)                
                      do jj = commd % bound_size(kk),commd % bound_size(kk+1)-1
                         ipoin        = commd % bound_invp(jj)
                         ipoin_offset = ipoin + offset
                         jj2           = (jj-1)*ndofn
                         do idofn = 1,ndofn
                            jj2 = jj2 + 1
                            xx(idofn,ipoin_offset) = xx(idofn,ipoin_offset) + recvbuff_rp(jj2)
                         end do
                      end do
                   end do
                   
                   do ipoin = commu % npoi1+1,commu % npoin
                      xx(1:ndofn,ipoin) = xx(1:ndofn,ipoin) + my_rxx(1:ndofn,ipoin-commu % npoi1) 
                   end do

                   do ii = commd % nneig_2,commd % nneig
                      dom_i = commd % neights_ordered(ii)
                      kk    = commd % perm_ordered(ii)
                      do jj = commd % bound_size(kk),commd % bound_size(kk+1)-1
                         ipoin        = commd % bound_invp(jj)
                         ipoin_offset = ipoin + offset
                         jj2          = (jj-1)*ndofn
                         do idofn = 1,ndofn
                            jj2 = jj2 + 1
                            xx(idofn,ipoin_offset) = xx(idofn,ipoin_offset) + recvbuff_rp(jj2)
                         end do
                      end do
                   end do

                   deallocate(my_rxx)

                else

                   kk = 0
                   do jj = 1,bound_dim
                      ipoin = bound_invp(jj)
                      ipoin_offset = ipoin + offset
                      do idofn = 1,ndofn
                         kk = kk + 1
                         xx(idofn,ipoin_offset) = xx(idofn,ipoin_offset) + recvbuff_rp(kk)
                      end do
                   end do
                end if

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' .or. trim(what) == 'MAXI' )&
                  & then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_invp(jj)
                   ipoin_offset = ipoin + offset
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin_offset) = max(xx(idofn,ipoin_offset),recvbuff_rp(kk))
                   end do
                end do

             else if( trim(what) == 'TAKE MIN' ) then
                !
                ! TAKE MIN
                !
                my_rank4 = commu % RANK4
                kk = 0
                do ii = 1,commu % nneig
                   dom_i = commu % neights(ii)
                   if( dom_i < my_rank4 ) then
                      ini   = ndofn * ( bound_size(ii)   - 1 ) + 1
                      nsize = ndofn * ( bound_size(ii+1) - 1 ) + 1 - ini
                      do jj = bound_size(ii),bound_size(ii+1)-1
                         ipoin = bound_perm(jj)
                         do idofn = 1,ndofn
                            kk = kk + 1
                            xx(idofn,ipoin) = recvbuff_rp(kk)
                         end do
                      end do
                   else
                      kk = kk + (bound_size(ii+1)-bound_size(ii))*ndofn
                   end if
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' .or. trim(what) == 'MINI' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_invp(jj)
                   ipoin_offset = ipoin + offset
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin_offset) = min(xx(idofn,ipoin_offset),recvbuff_rp(kk))
                   end do
                end do

             else if( trim(what) == 'DIF' .or. trim(what) == 'DIFFERENCE' .or. trim(what) == 'DIFF' ) then
                !
                ! DIF
                !
                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_invp(jj)
                   ipoin_offset = ipoin + offset
                   do idofn = 1,ndofn
                      kk = kk + 1
                      if( abs(xx(idofn,ipoin_offset)-recvbuff_rp(kk)) / (abs(xx(idofn,ipoin_offset))+epsilon(1.0_rp)) > 1.0e-12_rp ) then
                         if( ipoin <= npoi1 .or. ( ipoin >= commu % npoi2 .and. ipoin <= commu % npoi3 ) ) write(6,*) 'DIF=',xx(idofn,ipoin_offset)-recvbuff_rp(kk)
                         xx(idofn,ipoin_offset) = huge(1.0_rp)
                      end if
                   end do
                end do

             else if( trim(what) == 'MAX_DIF' .or. trim(what) == 'MAX_DIFFERENCE' .or. trim(what) == 'MAX_DIFF' ) then
                !
                ! MAX_DIF
                !
                nullify(bound_diff)
                call memory_alloca(par_memor,'bound_diff','par_send_receive_to_all_rp_0',bound_diff,ndofn,npoin)
                bound_diff(1:ndofn,1:npoin) = xx(1:ndofn,1:npoin)
                xx(1:ndofn,1:npoin) = 0.0_rp

                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_invp(jj)
                   ipoin_offset = ipoin + offset
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin_offset) = max( abs(bound_diff(idofn,ipoin)-recvbuff_rp(kk)),xx(idofn,ipoin_offset))
                   end do
                end do
                call memory_deallo(par_memor,'bound_diff','par_send_receive_to_all_rp_0',bound_diff)

             else
                call runend('UNKNOWN ORDER: '//trim(what))
             end if

             ipass = 0

          end if

       end if

    end if
#endif

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_RP

  subroutine PAR_SYMMETRIC_EXCHANGE_RP_1(xx,what,commu,wsynch,dom_k,NDOF)
   
    real(rp),             pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    class(comm_data_par_basic),     intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip),          optional, intent(in)    :: NDOF
    MY_MPI_COMM                                   :: PAR_COMM_TO_USE
    integer(ip)                                   :: ii,nsize,jj,dom_i
    integer(ip)                                   :: ipoin,ini,kk,idofn
    integer(ip)                                   :: offset,ipoin_offset,ndofn
    integer(4)                                    :: istat4,nsize4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    MY_MPI_STATUS ,                    pointer    :: status4(:)
    integer(ip)                                   :: dom_j
    integer(ip)                                   :: bound_dim
    integer(ip),                       pointer    :: bound_perm(:)
    integer(ip),                       pointer    :: bound_invp(:)
    integer(ip),                       pointer    :: bound_size(:)

#ifndef MPI_OFF
    if( associated(xx) ) then
       !
       ! Passes
       !
       ipass           =  ipass + 1
       dom_j           =  optional_argument(0_ip,dom_k)
       ndofn           =  optional_argument(1_ip,NDOF)
       PAR_COMM_TO_USE =  commu % PAR_COMM_WORLD
       bound_dim       =  commu % bound_dim
       bound_perm      => commu % bound_perm
       bound_size      => commu % bound_size
       select type( commu )
       class is ( comm_data_par_basic )
          bound_invp => commu % bound_perm
       class is ( comm_data_par       )
          if( associated(commu % bound_invp) ) then
             bound_invp => commu % bound_invp
          else
             bound_invp => commu % bound_perm
          end if
       end select
       !
       ! Offset
       !
       if( what(len(what):len(what)) /= 'I' ) then
          offset = 0
       else
          offset = -commu % npoi1
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if(       trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING'     ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( asynch ) allocate(ireq4(commu % nneig*2))
             if( memory_size(sendbuff_rp) < bound_dim * ndofn ) then
                call memory_deallo(par_memor,'SENDBUFF_RP','par_interface_node_exchange_rp',sendbuff_rp)
                call memory_alloca(par_memor,'SENDBUFF_RP','par_interface_node_exchange_rp',sendbuff_rp,bound_dim * ndofn,'DO_NOT_INITIALIZE')
             end if
             if( memory_size(recvbuff_rp) < bound_dim * ndofn ) then
                call memory_deallo(par_memor,'RECVBUFF_RP','par_interface_node_exchange_rp',recvbuff_rp)
                call memory_alloca(par_memor,'RECVBUFF_RP','par_interface_node_exchange_rp',recvbuff_rp,bound_dim * ndofn,'DO_NOT_INITIALIZE')
             end if
             !
             ! Distribute values among neighbors
             !
             if( trim(what) == 'DISTRIBUTE' ) then
                do jj = 1,bound_dim
                   ipoin = bound_perm(jj)
                   if( ipoin < commu % npoi2 .or. ipoin > commu % npoi3 ) then
                      do idofn = 1,ndofn
                         xx((ipoin-1)*ndofn+idofn) = 0.0_rp
                      end do
                   end if
                end do
             end if
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,bound_dim
                ipoin = bound_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   sendbuff_rp(kk) = xx((ipoin+offset-1)*ndofn+idofn)
                   recvbuff_rp(kk) = 0.0_rp
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             istat4 = 0_4
             kk = 0
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)
                
                if( dom_j == 0 .or. dom_j == dom_i ) then

                   ini   = ndofn * ( bound_size(ii)   - 1 ) + 1
                   nsize = ndofn * ( bound_size(ii+1) - 1 ) + 1 - ini

                   nsize4 = int(nsize,4)
                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           sendbuff_rp(ini:ini+nsize-1), nsize4, &
                           PAR_REAL,  dom_i4, 0_4, &
                           PAR_COMM_TO_USE, ireq4(kk), istat4  )
                      kk = kk + 1
                      call MPI_Irecv(&
                           recvbuff_rp(ini:ini+nsize-1), nsize4, &
                           PAR_REAL,  dom_i4, 0_4, &
                           PAR_COMM_TO_USE, ireq4(kk), istat4  )
                   else
                      call MPI_Sendrecv(                       &
                           sendbuff_rp(ini:), nsize4,            &
                           PAR_REAL, dom_i4, 0_4,  &
                           recvbuff_rp(ini:), nsize4, PAR_REAL, dom_i4, 0_4,&
                           PAR_COMM_TO_USE, status, istat4     )
                   end if
                   if( istat4 /= MPI_SUCCESS )call PAR_MPI_RUNEND(istat4,'PAR_INTERFACE_NODE_EXCHANGE_RP')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE*2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             deallocate( status4 )
             deallocate(ireq4)
          end if

         if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if(         trim(what) == 'SUM'  .or. trim(what) == 'ASSEMBLY' .or. trim(what) == 'SUMI'&
                  & .or. trim(what) == 'DISTRIBUTE' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_invp(jj)
                   ipoin_offset = ipoin + offset
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx((ipoin+offset-1)*ndofn+idofn) = xx((ipoin+offset-1)*ndofn+idofn) + recvbuff_rp(kk)
                   end do
                end do
                
             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' .or. trim(what) == 'MAXI' )&
                  & then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_invp(jj)
                   ipoin_offset = ipoin + offset
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx((ipoin+offset-1)*ndofn+idofn) = max(xx((ipoin+offset-1)*ndofn+idofn),recvbuff_rp(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' .or. trim(what) == 'MINI' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_invp(jj)
                   ipoin_offset = ipoin + offset
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx((ipoin+offset-1)*ndofn+idofn) = min(xx((ipoin+offset-1)*ndofn+idofn),recvbuff_rp(kk))
                   end do
                end do

             else
                call runend('UNKNOWN ORDER: '//trim(what))
             end if

             ipass = 0
 
          end if

       end if

    end if
#endif

  end subroutine PAR_SYMMETRIC_EXCHANGE_RP_1

    subroutine PAR_SYMMETRIC_EXCHANGE_RP_2(xx,what,commu,wsynch,dom_k)
   
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    class(comm_data_par_basic),     intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    integer(ip)                                   :: ii,nsize,jj,dom_i
    integer(ip)                                   :: ipoin,ini,kk,idofn
    integer(ip)                                   :: offset,ipoin_offset,ndofn
    integer(4)                                    :: istat4,nsize4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    MY_MPI_STATUS ,                    pointer    :: status4(:)
    integer(ip)                                   :: dom_j
    integer(ip)                                   :: bound_dim
    integer(ip),                       pointer    :: bound_perm(:)
    integer(ip),                       pointer    :: bound_invp(:)
    integer(ip),                       pointer    :: bound_size(:)

#ifndef MPI_OFF
    if( associated(xx) ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Pointers
       !
       ndofn           =  size(xx,1)
       PAR_COMM_TO_USE =  commu % PAR_COMM_WORLD
       bound_dim       =  commu % bound_dim
       bound_perm      => commu % bound_perm
       bound_size      => commu % bound_size
       select type( commu )
       class is ( comm_data_par_basic )
          bound_invp => commu % bound_perm
       class is ( comm_data_par       )
          if( associated(commu % bound_invp) ) then
             bound_invp => commu % bound_invp
          else
             bound_invp => commu % bound_perm
          end if
       end select
       !
       ! Offset
       !
       if( what(len(what):len(what)) /= 'I' ) then
          offset = 0
       else
          offset = -commu % npoi1
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if(       trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING'     ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( asynch ) allocate(ireq4(commu % nneig*2))
             if( memory_size(sendbuff_rp) < bound_dim * ndofn ) then
                call memory_deallo(par_memor,'SENDBUFF_RP','par_interface_node_exchange_rp',sendbuff_rp)
                call memory_alloca(par_memor,'SENDBUFF_RP','par_interface_node_exchange_rp',sendbuff_rp,bound_dim * ndofn,'DO_NOT_INITIALIZE')
             end if
             if( memory_size(recvbuff_rp) < bound_dim * ndofn ) then
                call memory_deallo(par_memor,'RECVBUFF_RP','par_interface_node_exchange_rp',recvbuff_rp)
                call memory_alloca(par_memor,'RECVBUFF_RP','par_interface_node_exchange_rp',recvbuff_rp,bound_dim * ndofn,'DO_NOT_INITIALIZE')
             end if
             !
             ! Distribute values among neighbors
             !
             if( trim(what) == 'DISTRIBUTE' ) then
                do jj = 1,bound_dim
                   ipoin = bound_perm(jj)
                   if( ipoin < commu % npoi2 .or. ipoin > commu % npoi3 ) then
                      do idofn = 1,ndofn
                         xx(idofn,ipoin) = 0.0_rp
                      end do
                   end if
                end do
             end if
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,bound_dim
                ipoin = bound_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   sendbuff_rp(kk) = xx(idofn,ipoin+offset)
                   recvbuff_rp(kk) = 0.0_rp
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             istat4 = 0_4
             kk = 0
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)
                
                if( dom_j == 0 .or. dom_j == dom_i ) then

                   ini   = ndofn * ( bound_size(ii)   - 1 ) + 1
                   nsize = ndofn * ( bound_size(ii+1) - 1 ) + 1 - ini

                   nsize4 = int(nsize,4)
                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           sendbuff_rp(ini:ini+nsize-1), nsize4, &
                           PAR_REAL,  dom_i4, 0_4, &
                           PAR_COMM_TO_USE, ireq4(kk), istat4  )
                      kk = kk + 1
                      call MPI_Irecv(&
                           recvbuff_rp(ini:ini+nsize-1), nsize4, &
                           PAR_REAL,  dom_i4, 0_4, &
                           PAR_COMM_TO_USE, ireq4(kk), istat4  )
                   else
                      call MPI_Sendrecv(                       &
                           sendbuff_rp(ini:), nsize4,            &
                           PAR_REAL, dom_i4, 0_4,  &
                           recvbuff_rp(ini:), nsize4, PAR_REAL, dom_i4, 0_4,&
                           PAR_COMM_TO_USE, status, istat4     )
                   end if
                   if( istat4 /= MPI_SUCCESS )call PAR_MPI_RUNEND(istat4,'PAR_INTERFACE_NODE_EXCHANGE_RP')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE*2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             deallocate( status4 )
             deallocate(ireq4)
          end if

         if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if(         trim(what) == 'SUM'  .or. trim(what) == 'ASSEMBLY' .or. trim(what) == 'SUMI'&
                  & .or. trim(what) == 'DISTRIBUTE' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_invp(jj)
                   ipoin_offset = ipoin + offset
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin+offset) = xx(idofn,ipoin+offset) + recvbuff_rp(kk)
                   end do
                end do
                
             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' .or. trim(what) == 'MAXI' )&
                  & then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_invp(jj)
                   ipoin_offset = ipoin + offset
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin+offset) = max(xx(idofn,ipoin+offset),recvbuff_rp(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' .or. trim(what) == 'MINI' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_invp(jj)
                   ipoin_offset = ipoin + offset
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin+offset) = min(xx(idofn,ipoin+offset),recvbuff_rp(kk))
                   end do
                end do

             else
                call runend('UNKNOWN ORDER: '//trim(what))
             end if

             ipass = 0
 
          end if

       end if

    end if
#endif

  end subroutine PAR_SYMMETRIC_EXCHANGE_RP_2
  
  subroutine PAR_INTERFACE_EXCHANGE(ndofn,xx,commu,wsynch)

    integer(ip),                    intent(in)    :: ndofn
    real(rp),                       intent(inout) :: xx(ndofn,*)
    type(comm_data_par),            intent(in)    :: commu
    character(*),   optional,       intent(in)    :: wsynch
    integer(ip)                                   :: ii,nsize,jj,dom_i
    integer(ip)                                   :: ipoin,ini,kk,idofn
    integer(4)                                    :: istat4,nsize4,count4
    integer(4)                                    :: dom_i4
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    MY_MPI_STATUS ,                    pointer    :: status4(:)
    integer(ip)                                   :: dom_j
    integer(ip)                                   :: bound_dim
    integer(ip),                       pointer    :: bound_perm(:)
    integer(ip),                       pointer    :: bound_invp(:)
    integer(ip),                       pointer    :: bound_size(:)

#ifndef MPI_OFF
    
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       dom_j = 0
       !
       ! Exchange arrays
       !
       bound_dim  =  commu % bound_dim
       bound_perm => commu % bound_perm
       bound_invp => commu % bound_invp
       bound_size => commu % bound_size
       PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( commu % nneig > 0 ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( asynch ) allocate(ireq4(commu % nneig*2))
             if( memory_size(sendbuff_rp) < bound_dim * ndofn ) then
                call memory_deallo(par_memor,'SENDBUFF_RP','par_interface_node_exchange_rp',sendbuff_rp)
                call memory_alloca(par_memor,'SENDBUFF_RP','par_interface_node_exchange_rp',sendbuff_rp,bound_dim * ndofn,'DO_NOT_INITIALIZE')
             end if
             if( memory_size(recvbuff_rp) < bound_dim * ndofn ) then
                call memory_deallo(par_memor,'RECVBUFF_RP','par_interface_node_exchange_rp',recvbuff_rp)
                call memory_alloca(par_memor,'RECVBUFF_RP','par_interface_node_exchange_rp',recvbuff_rp,bound_dim * ndofn,'DO_NOT_INITIALIZE')
             end if
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,bound_dim
                ipoin = bound_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   sendbuff_rp(kk) = xx(idofn,ipoin)
                   recvbuff_rp(kk) = 0.0_rp
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             istat4 = 0_4
             kk = 0
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)
                
                ini   = ndofn * ( bound_size(ii)   - 1 ) + 1
                nsize = ndofn * ( bound_size(ii+1) - 1 ) + 1 - ini
                
                nsize4 = int(nsize,4)
                if( asynch ) then
                   kk = kk + 1
                   call MPI_Isend(&
                        sendbuff_rp(ini:ini+nsize-1), nsize4, &
                        PAR_REAL,  dom_i4, 0_4, &
                        PAR_COMM_TO_USE, ireq4(kk), istat4  )
                   kk = kk + 1
                   call MPI_Irecv(&
                        recvbuff_rp(ini:ini+nsize-1), nsize4, &
                        PAR_REAL,  dom_i4, 0_4, &
                        PAR_COMM_TO_USE, ireq4(kk), istat4  )
                else
                   call MPI_Sendrecv(                       &
                        sendbuff_rp(ini:), nsize4,            &
                        PAR_REAL, dom_i4, 0_4,  &
                        recvbuff_rp(ini:), nsize4, PAR_REAL, dom_i4, 0_4,&
                        PAR_COMM_TO_USE, status, istat4     )
                end if
                if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_INTERFACE_NODE_EXCHANGE_RP')

             end do

          end if
          !
          ! Wait is asynchronous
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE*2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             deallocate( status4 )
             deallocate( ireq4   )
          end if
          !
          ! SUM
          !
          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then
             
            kk = 0
            do jj = 1,bound_dim
               ipoin = bound_invp(jj)
               do idofn = 1,ndofn
                  kk = kk + 1
                  xx(idofn,ipoin) = xx(idofn,ipoin) + recvbuff_rp(kk)
               end do
            end do            

            ipass = 0
            
         end if
      end if

   end if
#endif

  end subroutine PAR_INTERFACE_EXCHANGE

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_INTERFACE_OWN_NODE_EXCHANGE
  !
  ! KPOIN enables to limit the positions in the receiving arrays.
  ! It can be used for example not to write on halos, for example if XX
  ! is only defined until own+boundary nodes
  !
  !----------------------------------------------------------------------

  subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_IP_0(ndofn,xx,message,kpoin)
    integer(ip),                    intent(in)    :: ndofn
    integer(ip),                    intent(inout) :: xx(ndofn,*)
    character(*),         optional, intent(in)    :: message
    integer(ip),          optional, intent(in)    :: kpoin
    call PAR_INTERFACE_OWN_NODE_EXCHANGE_IP(ndofn,xx,message,kpoin)
  end subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_IP_0

  subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_IP_00(ndofn,xx,message,kpoin)
    integer(ip),                    intent(in)    :: ndofn
    integer(ip),                    intent(inout) :: xx(*)
    character(*),         optional, intent(in)    :: message
    integer(ip),          optional, intent(in)    :: kpoin
    call PAR_INTERFACE_OWN_NODE_EXCHANGE_IP(ndofn,xx,message,kpoin)
  end subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_IP_00

  subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_IP_1(xx,message,kpoin)
    integer(ip),          pointer,  intent(inout) :: xx(:)
    character(*),         optional, intent(in)    :: message
    integer(ip),          optional, intent(in)    :: kpoin
    integer(ip)                                   :: ndofn
    if( .not. associated(xx) ) return
    ndofn = 1
    call PAR_INTERFACE_OWN_NODE_EXCHANGE_IP(ndofn,xx,message,kpoin)
  end subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_IP_1

  subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_IP_2(xx,message,kpoin)
    
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),         optional, intent(in)    :: message
    integer(ip),          optional, intent(in)    :: kpoin
    integer(ip)                                   :: ndofn
    if( .not. associated(xx) ) return
    ndofn = size(xx,1)
    call PAR_INTERFACE_OWN_NODE_EXCHANGE_IP(ndofn,xx,message,kpoin)
  end subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_IP_2

  subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_RP_00(ndofn,xx,message,kpoin)
    integer(ip),                    intent(in)    :: ndofn
    real(rp),                       intent(inout) :: xx(*)
    character(*),         optional, intent(in)    :: message
    integer(ip),          optional, intent(in)    :: kpoin
    call PAR_INTERFACE_OWN_NODE_EXCHANGE_RP(ndofn,xx,message,kpoin)
  end subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_RP_00

  subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_RP_0(ndofn,xx,message,kpoin)
    integer(ip),                    intent(in)    :: ndofn
    real(rp),                       intent(inout) :: xx(ndofn,*)
    character(*),         optional, intent(in)    :: message
    integer(ip),          optional, intent(in)    :: kpoin
    call PAR_INTERFACE_OWN_NODE_EXCHANGE_RP(ndofn,xx,message,kpoin)
  end subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_RP_0

  subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_RP_1(xx,message,kpoin)
    real(rp),             pointer,  intent(inout) :: xx(:)
    character(*),         optional, intent(in)    :: message
    integer(ip),          optional, intent(in)    :: kpoin
    integer(ip)                                   :: ndofn
    if( .not. associated(xx) ) return
    ndofn = 1
    call PAR_INTERFACE_OWN_NODE_EXCHANGE_RP(ndofn,xx,message,kpoin)
  end subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_RP_1

  subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_RP_2(xx,message,kpoin)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),         optional, intent(in)    :: message
    integer(ip),          optional, intent(in)    :: kpoin
    integer(ip)                                   :: ndofn
    if( .not. associated(xx) ) return
    ndofn = size(xx,1)
    call PAR_INTERFACE_OWN_NODE_EXCHANGE_RP(ndofn,xx,message,kpoin)
  end subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_RP_2

  subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_RP(ndofn,xx,message,kpoin)
    
    integer(ip),                    intent(in)    :: ndofn
    real(rp),                       intent(inout) :: xx(ndofn,*)
    character(*), optional,         intent(in)    :: message
    integer(ip),  optional,         intent(in)    :: kpoin
    integer(ip)                                   :: ii,dom_i
    integer(ip)                                   :: ipoin,idofn
    integer(ip)                                   :: ii_recv_size,ii_send_size
    integer(ip)                                   :: ii_recv,ii_send,ineig
    integer(ip)                                   :: size_non_blocking
    real(rp),                          target     :: tmp_rsend_ok(2)
    real(rp),                          target     :: tmp_rrecv_ok(2)
    logical(lg)                                   :: do_send_receive
    logical(lg)                                   :: do_wait_and_assemble

    do_send_receive      = .false.
    do_wait_and_assemble = .false.

    if( present(message) ) then
       if( trim(message) == 'SEND RECEIVE'              ) do_send_receive      = .true.
       if( trim(message) == 'WAIT AND ASSEMBLE'         ) do_wait_and_assemble = .true.
       if( trim(message) == 'SEND RECEIVE AND ASSEMBLE' ) then
          do_send_receive      = .true.
          do_wait_and_assemble = .true.
       end if
    else
       do_send_receive      = .true.
       do_wait_and_assemble = .true.
    end if

    if( do_send_receive ) then
       !
       ! Fill-in sending array
       !
       allocate(tmp_rsend(max(commd % full_row_send_dim * ndofn,1_ip)))
       allocate(tmp_rrecv(max(commd % full_row_recv_dim * ndofn,1_ip)))

       do ii = 1,commd % full_row_send_dim
          ipoin = commd % full_row_send_perm(ii)
          do idofn = 1,ndofn
             !if( (ii-1)*ndofn+idofn > size(tmp_rsend) ) call runend('TROUBLE 0')
             tmp_rsend((ii-1)*ndofn+idofn) = xx(idofn,ipoin)
          end do
       end do
       !
       ! Maximum number of non-blocking operations
       !
       size_non_blocking = commd % full_row_send_nneig + commd % full_row_recv_nneig
       !
       ! Send and receive
       !
       call PAR_START_NON_BLOCKING_COMM(1_ip,size_non_blocking)     ! Set number of requests
       call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)                  ! Set non-blocking communication number to 1
      
       do ineig = 1,commd % full_row_send_nneig
          dom_i        =   commd % full_row_send_neights(ineig)
          ii_send_size = ( commd % full_row_send_size(ineig+1)-commd % full_row_send_size(ineig) ) * ndofn
          ii_send      = ( commd % full_row_send_size(ineig)-1)*ndofn+1
          ii_recv_size = 0
          !if(ii_send < 1 .or. ii_send > size(tmp_rsend) )   call runend('TROUBLE 1')
          !if(ii_send + ii_send_size - 1 > size(tmp_rsend) ) call runend('TROUBLE 2')
          call PAR_SEND_RECEIVE(ii_send_size,ii_recv_size,tmp_rsend(ii_send:),tmp_rrecv_ok,'IN MY CODE',dom_i,'ASYNCHRONOUS')
       end do
       do ineig = 1,commd % full_row_recv_nneig
          dom_i        =   commd % full_row_recv_neights(ineig)
          ii_recv_size = ( commd % full_row_recv_size(ineig+1)-commd % full_row_recv_size(ineig) ) * ndofn
          ii_recv      = ( commd % full_row_recv_size(ineig)-1)*ndofn+1
          ii_send_size = 0
          !if(ii_recv < 1 .or. ii_recv > size(tmp_rrecv) )   call runend('TROUBLE 3')
          !if(ii_recv + ii_recv_size - 1 > size(tmp_rrecv) ) call runend('TROUBLE 4')
          call PAR_SEND_RECEIVE(ii_send_size,ii_recv_size,tmp_rsend_ok,tmp_rrecv(ii_recv:),'IN MY CODE',dom_i,'ASYNCHRONOUS')
       end do 

    end if

    if( do_wait_and_assemble ) then
       !
       ! Wait and assemble
       !
       call PAR_END_NON_BLOCKING_COMM(1_ip)

       if( present(kpoin) ) then
          do ii = 1,commd % full_row_recv_dim
             ipoin = commd % full_row_recv_perm(ii)
             if( ipoin <= kpoin ) then
                do idofn = 1,ndofn
                   xx(idofn,ipoin) = tmp_rrecv((ii-1)*ndofn+idofn)
                end do
             end if
          end do
       else
          do ii = 1,commd % full_row_recv_dim
             ipoin = commd % full_row_recv_perm(ii)
             do idofn = 1,ndofn
                !if( (ii-1)*ndofn+idofn > size(tmp_rrecv) ) call runend('TROUBLE 5')
                xx(idofn,ipoin) = tmp_rrecv((ii-1)*ndofn+idofn)
             end do
          end do
       end if

       deallocate(tmp_rsend)
       deallocate(tmp_rrecv)

    end if

  end subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_RP

  subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_IP(ndofn,xx,message,kpoin)
    
    integer(ip),                    intent(in)    :: ndofn
    integer(ip),                    intent(inout) :: xx(ndofn,*)
    character(*), optional,         intent(in)    :: message
    integer(ip),  optional,         intent(in)    :: kpoin
    integer(ip)                                   :: ii,dom_i
    integer(ip)                                   :: ipoin,idofn
    integer(ip)                                   :: ii_recv_size,ii_send_size
    integer(ip)                                   :: ii_recv,ii_send,ineig
    integer(ip)                                   :: size_non_blocking
    integer(ip),                       target     :: tmp_isend_ok(2)
    integer(ip),                       target     :: tmp_irecv_ok(2)
    logical(lg)                                   :: do_send_receive
    logical(lg)                                   :: do_wait_and_assemble

    do_send_receive      = .false.
    do_wait_and_assemble = .false.

    if( present(message) ) then
       if( trim(message) == 'SEND RECEIVE'      ) do_send_receive      = .true.
       if( trim(message) == 'WAIT AND ASSEMBLE' ) do_wait_and_assemble = .true.
    else
       do_send_receive      = .true.
       do_wait_and_assemble = .true.
    end if

    if( do_send_receive ) then
       !
       ! Send receive
       !
       allocate(tmp_isend(max(commd % full_row_send_dim * ndofn,1_ip)))
       allocate(tmp_irecv(max(commd % full_row_recv_dim * ndofn,1_ip)))

       do ii = 1,commd % full_row_send_dim
          ipoin = commd % full_row_send_perm(ii)
          do idofn = 1,ndofn
             tmp_isend((ii-1)*ndofn+idofn) = xx(idofn,ipoin)
          end do
       end do
       !
       ! Maximum number of non-blocking operations
       !
       size_non_blocking = commd % full_row_send_nneig + commd % full_row_recv_nneig
       !
       ! Send and receive
       !
       call PAR_START_NON_BLOCKING_COMM(1_ip,size_non_blocking)     ! Set number of requests
       call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)                  ! Set non-blocking communication number to 1

       do ineig = 1,commd % full_row_send_nneig
          dom_i        =   commd % full_row_send_neights(ineig)
          ii_send_size = ( commd % full_row_send_size(ineig+1)-commd % full_row_send_size(ineig) ) * ndofn
          ii_send      = ( commd % full_row_send_size(ineig)-1)*ndofn+1
          ii_recv_size = 0
          call PAR_SEND_RECEIVE(ii_send_size,ii_recv_size,tmp_isend(ii_send:),tmp_irecv_ok,'IN MY CODE',dom_i,'ASYNCHRONOUS')
       end do
       do ineig = 1,commd % full_row_recv_nneig
          dom_i        =   commd % full_row_recv_neights(ineig)
          ii_recv_size = ( commd % full_row_recv_size(ineig+1)-commd % full_row_recv_size(ineig) ) * ndofn
          ii_recv      = ( commd % full_row_recv_size(ineig)-1)*ndofn+1
          ii_send_size = 0
          call PAR_SEND_RECEIVE(ii_send_size,ii_recv_size,tmp_isend_ok,tmp_irecv(ii_recv:),'IN MY CODE',dom_i,'ASYNCHRONOUS')
       end do

    end if

    if( do_wait_and_assemble ) then
       !
       ! Wait and assemble
       !
       call PAR_END_NON_BLOCKING_COMM(1_ip)

       if( present(kpoin) ) then
          do ii = 1,commd % full_row_recv_dim
             ipoin = commd % full_row_recv_perm(ii)
             if( ipoin <= kpoin ) then
                do idofn = 1,ndofn
                   xx(idofn,ipoin) = tmp_irecv((ii-1)*ndofn+idofn)
                end do
             end if
          end do
       else
          do ii = 1,commd % full_row_recv_dim
             ipoin = commd % full_row_recv_perm(ii)
             do idofn = 1,ndofn
                xx(idofn,ipoin) = tmp_irecv((ii-1)*ndofn+idofn)
             end do
          end do
       end if

       deallocate(tmp_isend)
       deallocate(tmp_irecv)

    end if

  end subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_IP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_INTERFACE_NODE_EXCHANGE_RP
  !
  !----------------------------------------------------------------------

  subroutine PAR_COUPLING_NODE_EXCHANGE_RP_0(n,xx,what,wherein,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_COUPLING_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_COUPLING_NODE_EXCHANGE_RP_0

  subroutine PAR_COUPLING_NODE_EXCHANGE_RP_1(xx,what,wherein,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = 1
    if( size(xx,1) /= npoin .and. size(xx,1) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE RP_1')
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_COUPLING_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_COUPLING_NODE_EXCHANGE_RP_1

  subroutine PAR_COUPLING_NODE_EXCHANGE_RP_1b(xx,what,commu,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = 1
    if( size(xx,1) /= npoin .and. size(xx,1) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE RP_1b')
    PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
    call PAR_COUPLING_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_COUPLING_NODE_EXCHANGE_RP_1b

  subroutine PAR_COUPLING_NODE_EXCHANGE_RP_2(xx,what,wherein,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    if( size(xx,2) /= npoin .and. size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE RP_2')
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_COUPLING_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_COUPLING_NODE_EXCHANGE_RP_2

  subroutine PAR_COUPLING_NODE_EXCHANGE_RP_3(xx,what,wherein,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    if( size(xx,3) <= 2 ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= npoin .and. size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE IP_3')
    else
       ndofn = size(xx,1)*size(xx,2)
       if( size(xx,3) /= npoin .and. size(xx,3) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE RP_3')
    end if
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_COUPLING_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_COUPLING_NODE_EXCHANGE_RP_3

  subroutine PAR_COUPLING_NODE_EXCHANGE_RP_2b(xx,what,commu,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    if( size(xx,2) /= npoin .and. size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
    PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
    call PAR_COUPLING_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_COUPLING_NODE_EXCHANGE_RP_2b

  !----------------------------------------------------------------------
  !
  ! PAR_COUPLING_NODE_EXCHANGE_RP: NODE ASSEMBLY FOR REALS
  !
  ! From the receiver point of view
  ! +---+
  ! | 1 |
  ! +---+                                                    +---+  /|
  ! | 2 |                                                    |   |  ||
  ! +---+   lscat_dim                                        +---+  ||
  ! | 3 |                                                    |   |  ||
  ! +---+      /\  +---+    +---+---+---+---+---+---+---+    +---+  || lrecv_size(2)-lrecv_size(1)
  ! | 4 |      ||  |   |    |   |   |   |   |   |   |   |    |   |  ||
  ! +---+      ||  +---+    +---+---+---+---+---+---+---+    +---+  ||
  ! | 5 |      ||  |   |  = |   |   |   |   |   |   |   |    |   |  ||
  ! +---+      ||  +---+    +---+---+---+---+---+---+---+    +---+  \/
  ! | 6 |      ||  |   |    |   |   |   |   |   |   |   |    +---+  /|
  ! +---+      \/  +---+    +---+---+---+---+---+---+---+    |   |  ||
  ! | 7 |            |             matrix_ia(:)              +---+  ||
  ! +---+  <---------+             matrix_ja(:)              |   |  || lrecv_size(3)-lrecv_size(2)
  ! | 8 |    lscat_perm(:)         matrix_aa(:)              +---+  ||
  ! +---+                          matrix_nzdom              |   |  ||
  ! | 9 |                                                    +---+  \/
  ! +---+
  !
  ! Nodal array
  !
  !----------------------------------------------------------------------

  subroutine PAR_COUPLING_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

    
    integer(ip),                    intent(in)    :: ndofn
    real(rp),                       intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    MY_MPI_COMM,                    intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,nsize_recv,jj,dom_i
    integer(ip)                                   :: ipoin,ini_recv,kk,idofn
    integer(4)                                    :: istat4,nsize4_send,nsize4_recv,count4
    integer(4)                                    :: ini_send,nsize_send
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    MY_MPI_STATUS ,                    pointer    :: status4(:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_rsend(commu % lsend_dim * ndofn))
             allocate(tmp_rrecv(commu % lrecv_dim * ndofn))
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,commu % lsend_dim
                ipoin = commu % lsend_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_rsend(kk) = xx(idofn,ipoin)
                end do
             end do
             do kk = 1,commu % lrecv_dim * ndofn
                tmp_rrecv(kk) = 0.0_rp
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             istat4 = 0_4
             kk = 0
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   ini_send   = ndofn * ( commu % lsend_size(ii)   - 1 ) + 1
                   nsize_send = ndofn * ( commu % lsend_size(ii+1) - 1 ) + 1 - ini_send
                   ini_recv   = ndofn * ( commu % lrecv_size(ii)   - 1 ) + 1
                   nsize_recv = ndofn * ( commu % lrecv_size(ii+1) - 1 ) + 1 - ini_recv

                   nsize4_send = int(nsize_send,4)
                   nsize4_recv = int(nsize_recv,4)

                   if( asynch ) then
                      !if( nsize_send > 0 ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_rsend(ini_send:ini_send+nsize_send-1), nsize4_send, &
                           PAR_REAL,  dom_i4, 0_4, &
                           PAR_COMM_TO_USE, ireq4(kk), istat4 )
                      !end if
                      !if( nsize_recv > 0 ) then
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_rrecv(ini_recv:ini_recv+nsize_recv-1), nsize4_recv, &
                           PAR_REAL,  dom_i4, 0_4, &
                           PAR_COMM_TO_USE, ireq4(kk), istat4 )
                      !end if
                   else
                      call MPI_Sendrecv(                       &
                           tmp_rsend(ini_send:), nsize4_send,  &
                           PAR_REAL, dom_i4, 0_4,  &
                           tmp_rrecv(ini_recv:), nsize4_recv,  &
                           PAR_REAL, dom_i4, 0_4,  &
                           PAR_COMM_TO_USE, status, istat4    )
                   end if
                   if( istat4 /= MPI_SUCCESS )call PAR_MPI_RUNEND(istat4,'PAR_COUPLING_NODE_EXCHANGE_RP')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE*2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'MATRIX' ) then
                !
                ! MATRIX
                !
                kk = 0
                do ii = 1,commu % lscat_dim
                   ipoin = commu % lscat_perm(ii)
                   do kk = commu % matrix_ia(ii),commu % matrix_ia(ii+1)-1
                      jj = commu % matrix_ja(kk)
                      xx(1:ndofn,ipoin) = xx(1:ndofn,ipoin) + commu % matrix_aa(kk) * tmp_rrecv( jj:jj )
                   end do
                end do

             else
                call runend('UNKNOWN ORDER:'//trim(what))
             end if

             ipass = 0
             deallocate(tmp_rrecv)
             deallocate(tmp_rsend)

          end if

       end if

    end if
#endif

  end subroutine PAR_COUPLING_NODE_EXCHANGE_RP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_SEND and PAR_RECEIVE
  !
  !----------------------------------------------------------------------

  subroutine PAR_RECEIVE_IP_s(xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN)
    
    integer(ip),           intent(out)          :: xx_recv
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in), optional :: dom_i
    character(*),          intent(in), optional :: wsynch
    MY_MPI_COMM   ,       intent(in), optional :: PAR_COMM_IN
    integer(ip)                                 :: yy_send(2)
    integer(ip)                                 :: yy_recv(2)

    if( ISEQUEN ) return
    if( present(wsynch) ) then
       call PAR_SEND_RECEIVE_IP(0_ip,1_ip,yy_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
    else
       call PAR_SEND_RECEIVE_IP(0_ip,1_ip,yy_send,yy_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
    end if
    xx_recv    = yy_recv(1)

  end subroutine PAR_RECEIVE_IP_s

  subroutine PAR_RECEIVE_RP_s(xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN)
    
    real(rp),              intent(out)          :: xx_recv
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in), optional :: dom_i
    character(*),          intent(in), optional :: wsynch
    MY_MPI_COMM   ,       intent(in), optional :: PAR_COMM_IN
    real(rp)                                    :: yy_send(2)
    real(rp)                                    :: yy_recv(2)

    if( ISEQUEN ) return
    if( present(wsynch) ) then
       call PAR_SEND_RECEIVE_RP(0_ip,1_ip,yy_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
    else
       call PAR_SEND_RECEIVE_RP(0_ip,1_ip,yy_send,yy_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
    end if
    xx_recv    = yy_recv(1)

  end subroutine PAR_RECEIVE_RP_s

  subroutine PAR_RECEIVE_RP_0(nrecv,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN)

    integer(ip),           intent(in)           :: nrecv
    real(rp),              intent(inout)        :: xx_recv(:)
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in), optional :: dom_i
    character(*),          intent(in), optional :: wsynch
    MY_MPI_COMM   ,        intent(in), optional :: PAR_COMM_IN
    real(rp)                                    :: yy_send(2)
    
    if( ISEQUEN ) return
    if( nrecv > 0 ) then
       if( present(wsynch) ) then
          call PAR_SEND_RECEIVE_RP(0_ip,nrecv,yy_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else
          call PAR_SEND_RECEIVE_RP(0_ip,nrecv,yy_send,xx_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       end if
    end if
    
  end subroutine PAR_RECEIVE_RP_0

  subroutine PAR_RECEIVE_RP_1(xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN)
   
    real(rp),              intent(inout), pointer    :: xx_recv(:)
    character(*),          intent(in),    optional :: wherein
    integer(ip),           intent(in),    optional :: dom_i
    character(*),          intent(in),    optional :: wsynch
    MY_MPI_COMM   ,        intent(in),    optional :: PAR_COMM_IN
    real(rp)                                       :: yy_send(2)
    integer(ip)                                    :: nsize
    
    if( ISEQUEN ) return
    if( associated(xx_recv) ) then
       nsize = size(xx_recv)
       if( present(wsynch) ) then
          call PAR_SEND_RECEIVE_RP(0_ip,nsize,yy_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else
          call PAR_SEND_RECEIVE_RP(0_ip,nsize,yy_send,xx_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       end if
    end if
    
  end subroutine PAR_RECEIVE_RP_1

  subroutine PAR_RECEIVE_RP_2(xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN)
    
    real(rp),              intent(inout), pointer  :: xx_recv(:,:)
    character(*),          intent(in),    optional :: wherein
    integer(ip),           intent(in),    optional :: dom_i
    character(*),          intent(in),    optional :: wsynch
    MY_MPI_COMM   ,        intent(in),    optional :: PAR_COMM_IN
    real(rp)                                       :: yy_send(2,2)
    integer(ip)                                    :: nsize
    
    if( ISEQUEN ) return
    if( associated(xx_recv) ) then
       nsize = size(xx_recv)
       if( present(wsynch) ) then
          call PAR_SEND_RECEIVE_RP(0_ip,nsize,yy_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else
          call PAR_SEND_RECEIVE_RP(0_ip,nsize,yy_send,xx_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       end if
    end if
    
  end subroutine PAR_RECEIVE_RP_2

  subroutine PAR_RECEIVE_IP_0(nrecv,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN)
    
    integer(ip),           intent(in)           :: nrecv
    integer(ip),           intent(inout)        :: xx_recv(*)
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in), optional :: dom_i
    character(*),          intent(in), optional :: wsynch
    MY_MPI_COMM   ,        intent(in), optional :: PAR_COMM_IN
    integer(ip)                                 :: yy_send(2)
    
    if( ISEQUEN ) return
    if( nrecv > 0 ) then
       if( present(wsynch) ) then
          call PAR_SEND_RECEIVE_IP(0_ip,nrecv,yy_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else
          call PAR_SEND_RECEIVE_IP(0_ip,nrecv,yy_send,xx_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       end if
    end if
    
  end subroutine PAR_RECEIVE_IP_0

  subroutine PAR_RECEIVE_IP_1(xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN)
    
    integer(ip),           intent(inout), pointer  :: xx_recv(:)
    character(*),          intent(in),    optional :: wherein
    integer(ip),           intent(in),    optional :: dom_i
    character(*),          intent(in),    optional :: wsynch
    MY_MPI_COMM   ,        intent(in),    optional :: PAR_COMM_IN
    integer(ip)                                    :: yy_send(2)
    integer(ip)                                    :: nsize
    
    if( ISEQUEN ) return
    if( associated(xx_recv) ) then
       nsize = size(xx_recv)
       if( present(wsynch) ) then
          call PAR_SEND_RECEIVE_IP(0_ip,nsize,yy_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else
          call PAR_SEND_RECEIVE_IP(0_ip,nsize,yy_send,xx_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       end if
    end if
    
  end subroutine PAR_RECEIVE_IP_1

  subroutine PAR_RECEIVE_IP_2(xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN)
    
    integer(ip),           intent(inout), pointer :: xx_recv(:,:)
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in), optional :: dom_i
    character(*),          intent(in), optional :: wsynch
    MY_MPI_COMM   ,       intent(in), optional :: PAR_COMM_IN
    integer(ip)                                 :: yy_send(2,2)
    integer(ip)                                 :: nsize
    
    if( ISEQUEN ) return
    if( associated(xx_recv) ) then
       nsize = size(xx_recv)
       if( present(wsynch) ) then
          call PAR_SEND_RECEIVE_IP(0_ip,nsize,yy_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else
          call PAR_SEND_RECEIVE_IP(0_ip,nsize,yy_send,xx_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       end if
    end if
    
  end subroutine PAR_RECEIVE_IP_2

  subroutine PAR_SEND_IP_s(xx_send,wherein,dom_i,wsynch,PAR_COMM_IN)
    
    integer(ip),           intent(in)           :: xx_send
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in), optional :: dom_i
    character(*),          intent(in), optional :: wsynch
    MY_MPI_COMM   ,       intent(in), optional :: PAR_COMM_IN
    integer(ip)                                 :: yy_send(2)
    integer(ip)                                 :: yy_recv(2)

    if( ISEQUEN ) return
    yy_send(1) = xx_send
    if( present(wsynch) ) then
       call PAR_SEND_RECEIVE_IP(1_ip,0_ip,yy_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
    else
       call PAR_SEND_RECEIVE_IP(1_ip,0_ip,yy_send,yy_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
    end if

  end subroutine PAR_SEND_IP_s

  subroutine PAR_SEND_RP_s(xx_send,wherein,dom_i,wsynch,PAR_COMM_IN)
    
    real(rp),              intent(in)           :: xx_send
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in), optional :: dom_i
    character(*),          intent(in), optional :: wsynch
    MY_MPI_COMM   ,       intent(in), optional :: PAR_COMM_IN
    real(rp)                                    :: yy_send(2)
    real(rp)                                    :: yy_recv(2)

    if( ISEQUEN ) return
    yy_send(1) = xx_send
    if( present(wsynch) ) then
       call PAR_SEND_RECEIVE_RP(1_ip,0_ip,yy_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
    else
       call PAR_SEND_RECEIVE_RP(1_ip,0_ip,yy_send,yy_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
    end if

  end subroutine PAR_SEND_RP_s
  
  subroutine PAR_SEND_RP_0(nsend,xx_send,wherein,dom_i,wsynch,PAR_COMM_IN)
   
    integer(ip),           intent(in)           :: nsend
    real(rp),              intent(in)           :: xx_send(*)
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in), optional :: dom_i
    character(*),          intent(in), optional :: wsynch
    MY_MPI_COMM   ,       intent(in), optional :: PAR_COMM_IN
    real(rp)                                    :: yy_recv(2)
    
    if( ISEQUEN ) return
    if( nsend > 0 ) then
       if( present(wsynch) ) then
          call PAR_SEND_RECEIVE_RP(nsend,0_ip,xx_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else
          call PAR_SEND_RECEIVE_RP(nsend,0_ip,xx_send,yy_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       end if
    end if
    
  end subroutine PAR_SEND_RP_0
  
  subroutine PAR_SEND_RP_1(xx_send,wherein,dom_i,wsynch,PAR_COMM_IN)
   
    real(rp),              intent(in), pointer  :: xx_send(:)
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in), optional :: dom_i
    character(*),          intent(in), optional :: wsynch
    MY_MPI_COMM   ,       intent(in), optional :: PAR_COMM_IN
    real(rp)                                    :: yy_recv(2)
    integer(ip)                                 :: nsize
    
    if( ISEQUEN ) return
    if( associated(xx_send) ) then
       nsize = size(xx_send)
       if( present(wsynch) ) then
          call PAR_SEND_RECEIVE_RP(nsize,0_ip,xx_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else
          call PAR_SEND_RECEIVE_RP(nsize,0_ip,xx_send,yy_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       end if
    end if
    
  end subroutine PAR_SEND_RP_1
  
  subroutine PAR_SEND_RP_2(xx_send,wherein,dom_i,wsynch,PAR_COMM_IN)
   
    real(rp),              intent(in), pointer  :: xx_send(:,:)
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in), optional :: dom_i
    character(*),          intent(in), optional :: wsynch
    MY_MPI_COMM   ,       intent(in), optional :: PAR_COMM_IN
    real(rp)                                    :: yy_recv(2,2)
    integer(ip)                                 :: nsize
    
    if( ISEQUEN ) return
    if( associated(xx_send) ) then
       nsize = size(xx_send)
       if( present(wsynch) ) then
          call PAR_SEND_RECEIVE_RP(nsize,0_ip,xx_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else
          call PAR_SEND_RECEIVE_RP(nsize,0_ip,xx_send,yy_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       end if
    end if
    
  end subroutine PAR_SEND_RP_2
  
  subroutine PAR_SEND_IP_0(nsend,xx_send,wherein,dom_i,wsynch,PAR_COMM_IN)
   
    integer(ip),           intent(in)           :: nsend
    integer(ip),           intent(in)           :: xx_send(*)
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in), optional :: dom_i
    character(*),          intent(in), optional :: wsynch
    MY_MPI_COMM   ,       intent(in), optional :: PAR_COMM_IN
    integer(ip)                                 :: yy_recv(2)
    
    if( ISEQUEN ) return
    if( nsend > 0 ) then
       if( present(wsynch) ) then
          call PAR_SEND_RECEIVE_IP(nsend,0_ip,xx_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else
          call PAR_SEND_RECEIVE_IP(nsend,0_ip,xx_send,yy_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       end if
    end if
    
  end subroutine PAR_SEND_IP_0
  
  subroutine PAR_SEND_IP_1(xx_send,wherein,dom_i,wsynch,PAR_COMM_IN)
   
    integer(ip),           intent(in), pointer  :: xx_send(:)
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in), optional :: dom_i
    character(*),          intent(in), optional :: wsynch
    MY_MPI_COMM   ,       intent(in), optional :: PAR_COMM_IN
    integer(ip)                                 :: yy_recv(2)
    integer(ip)                                 :: nsize
    
    if( ISEQUEN ) return
    if( associated(xx_send) ) then
       nsize = size(xx_send)
       if( present(wsynch) ) then
          call PAR_SEND_RECEIVE_IP(nsize,0_ip,xx_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else
          call PAR_SEND_RECEIVE_IP(nsize,0_ip,xx_send,yy_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       end if
    end if
    
  end subroutine PAR_SEND_IP_1
  
  subroutine PAR_SEND_IP_2(xx_send,wherein,dom_i,wsynch,PAR_COMM_IN)
   
    integer(ip),           intent(in), pointer  :: xx_send(:,:)
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in), optional :: dom_i
    character(*),          intent(in), optional :: wsynch
    MY_MPI_COMM   ,       intent(in), optional :: PAR_COMM_IN
    integer(ip)                                 :: yy_recv(2,2)
    integer(ip)                                 :: nsize
    
    if( ISEQUEN ) return
    if( associated(xx_send) ) then
       nsize = size(xx_send)
       if( present(wsynch) ) then
          call PAR_SEND_RECEIVE_IP(nsize,0_ip,xx_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else
          call PAR_SEND_RECEIVE_IP(nsize,0_ip,xx_send,yy_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       end if
    end if
    
  end subroutine PAR_SEND_IP_2
  
  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_SEND_RECEIVE_IP
  !
  !----------------------------------------------------------------------

  subroutine PAR_SEND_RECEIVE_IP_s(xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN)
    
    integer(ip),           intent(in)           :: xx_send
    integer(ip),           intent(out)          :: xx_recv
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in), optional :: dom_i
    character(*),          intent(in), optional :: wsynch
    MY_MPI_COMM   ,       intent(in), optional :: PAR_COMM_IN
    integer(ip)                                 :: yy_send(2)
    integer(ip)                                 :: yy_recv(2)
    if( ISEQUEN ) return
    yy_send(1) = xx_send
    if( present(wsynch) ) then
       call PAR_SEND_RECEIVE_IP(1_ip,1_ip,yy_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
    else
       call PAR_SEND_RECEIVE_IP(1_ip,1_ip,yy_send,yy_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
    end if
    xx_recv    = yy_recv(1)
  end subroutine PAR_SEND_RECEIVE_IP_s

  subroutine PAR_SEND_RECEIVE_IP_0(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN)
    
    integer(ip),           intent(in)           :: nsend
    integer(ip),           intent(in)           :: nrecv
    integer(ip),           intent(in)           :: xx_send(*)
    integer(ip),           intent(out)          :: xx_recv(*)
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in)           :: dom_i
    character(*),          intent(in), optional :: wsynch
    MY_MPI_COMM   ,       intent(in), optional :: PAR_COMM_IN
    integer(ip)                                 :: yy_send(2)
    integer(ip)                                 :: yy_recv(2)
    if( ISEQUEN ) return
    if( present(wsynch) ) then
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else
          return
       end if
    else
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       else
          return
       end if
    end if
  end subroutine PAR_SEND_RECEIVE_IP_0

  subroutine PAR_SEND_RECEIVE_IP_1(xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN)
    
    integer(ip),  pointer, intent(in)           :: xx_send(:)
    integer(ip),  pointer, intent(inout)          :: xx_recv(:)
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in), optional :: dom_i
    character(*),          intent(in), optional :: wsynch
    MY_MPI_COMM   ,       intent(in), optional :: PAR_COMM_IN
    integer(ip)                                 :: nsend,nrecv
    integer(ip)                                 :: yy_send(2)
    integer(ip)                                 :: yy_recv(2)
    if( ISEQUEN ) return
    if( .not. associated(xx_send) ) then
       nsend = 0
    else
       nsend = size(xx_send)
    end if
    if( .not. associated(xx_recv) ) then
       nrecv = 0
    else
       nrecv = size(xx_recv)
    end if
    if( present(wsynch) ) then
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else
          return
       end if
    else
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       else
          return
       end if
    end if
  end subroutine PAR_SEND_RECEIVE_IP_1

  subroutine PAR_SEND_RECEIVE_IP_2(xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN)
    
    integer(ip),  pointer, intent(in)           :: xx_send(:,:)
    integer(ip),  pointer, intent(inout)        :: xx_recv(:,:)
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in), optional :: dom_i
    character(*),          intent(in), optional :: wsynch
    MY_MPI_COMM   ,       intent(in), optional :: PAR_COMM_IN
    integer(ip)                                 :: nsend,nrecv
    integer(ip)                                 :: yy_send(2)
    integer(ip)                                 :: yy_recv(2)

    if( ISEQUEN ) return

    if( .not. associated(xx_send) ) then
       nsend = 0
    else
       nsend = size(xx_send,1)*size(xx_send,2)
    end if
    if( .not. associated(xx_recv) ) then
       nrecv = 0
    else
       nrecv = size(xx_recv,1)*size(xx_recv,2)
    end if
    if( present(wsynch) ) then
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else
          return
       end if
    else
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       else
          return
       end if
    end if
  end subroutine PAR_SEND_RECEIVE_IP_2
  
  subroutine PAR_SEND_RECEIVE_IP_3(xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN)
    
    integer(ip),  pointer, intent(in)           :: xx_send(:,:,:)
    integer(ip),  pointer, intent(inout)          :: xx_recv(:,:,:)
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in), optional :: dom_i
    character(*),          intent(in), optional :: wsynch
    MY_MPI_COMM   ,       intent(in), optional :: PAR_COMM_IN
    integer(ip)                                 :: nsend,nrecv
    integer(ip)                                 :: yy_send(2)
    integer(ip)                                 :: yy_recv(2)
    if( ISEQUEN ) return
    if( .not. associated(xx_send) ) then
       nsend = 0
    else
       nsend = size(xx_send,1)*size(xx_send,2)*size(xx_send,3)
    end if
    if( .not. associated(xx_recv) ) then
       nrecv = 0
    else
       nrecv = size(xx_recv,1)*size(xx_recv,2)*size(xx_recv,3)
    end if
    if( present(wsynch) ) then
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else
          return
       end if
    else
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       else
          return
       end if
    end if
  end subroutine PAR_SEND_RECEIVE_IP_3

  !----------------------------------------------------------------------
  !
  ! PAR_SEND_RECEIVE_IP
  !
  !----------------------------------------------------------------------

  subroutine PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN)
    
    integer(ip),            intent(in)  :: nsend
    integer(ip),            intent(in)  :: nrecv
    integer(ip),            intent(in)  :: xx_send(*)
    integer(ip),            intent(out) :: xx_recv(*)
    character(*), optional, intent(in)  :: wherein
    integer(ip),            intent(in)  :: dom_i
    character(*), optional, intent(in)  :: wsynch
    MY_MPI_COMM   ,  optional, intent(in)  :: PAR_COMM_IN
    integer(ip)                         :: kk,ii
    integer(4)                          :: istat4,nsend4,nrecv4,dom_i4
    integer(4)                          :: my_rank4
    logical(lg)                         :: asynch
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    !
    ! Define communicator
    !
    if( present(PAR_COMM_IN) ) then
       call PAR_COMM_RANK(PAR_COMM_TO_USE,MY_RANK4,PAR_COMM_IN)
    else if( present(wherein) ) then
       call PAR_COMM_RANK(PAR_COMM_TO_USE,MY_RANK4,wherein)
    else
       call PAR_COMM_RANK(PAR_COMM_TO_USE,MY_RANK4,'IN MY CODE')
    end if
    !
    ! BlockinG/non blocking
    !
    if( present(wsynch) ) then
       if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
          asynch = .false.
       else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
          asynch = .true.
          if( .not. associated(non_blocking(inonblocking) % request4) ) then
             call runend('NON-BLOCKING SEND/RECEIVE SHOULD BE STARTED')
          end if
       else
          call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
       end if
    else
       asynch = .false.
    end if
    !
    ! Synchronous Send/receive
    !
    nsend4 = int(nsend,4)
    nrecv4 = int(nrecv,4)
    dom_i4 = int(dom_i,4)
    istat4 = 0

#ifndef MPI_OFF
    if( asynch ) then
       if( nrecv /= 0 ) then
          non_blocking(inonblocking) % count4 = non_blocking(inonblocking) % count4 + 1
          kk = non_blocking(inonblocking) % count4
          call MPI_IRecv(                                 &
               xx_recv(1:nrecv), nrecv4,                  &
               PAR_INTEGER, dom_i4, 0_4,                  &
               PAR_COMM_TO_USE,                           &
               non_blocking(inonblocking) % request4(kk), &
               istat4                                     )
       end if
       if( nsend /= 0 ) then
          non_blocking(inonblocking) % count4 = non_blocking(inonblocking) % count4 + 1
          kk = non_blocking(inonblocking) % count4
          call MPI_ISend(                                 &
               xx_send(1:nsend), nsend4,                  &
               PAR_INTEGER, dom_i4, 0_4,                  &
               PAR_COMM_TO_USE,                           &
               non_blocking(inonblocking) % request4(kk), &
               istat4                                     )
       end if
    else
       if( nrecv /= 0 .and. nsend == 0 ) then
          call MPI_Recv(                          &
               xx_recv(1:nrecv), nrecv4,          &
               PAR_INTEGER, dom_i4, 0_4,          &
               PAR_COMM_TO_USE, status, istat4    )

       else if( nrecv == 0 .and. nsend /= 0 ) then
          call MPI_Send(                          &
               xx_send(1:nsend), nsend4,          &
               PAR_INTEGER, dom_i4, 0_4,          &
               PAR_COMM_TO_USE, istat4            )

       else if( nrecv /= 0 .and. nsend /= 0 ) then
          if( my_rank4 == dom_i4 ) then
             do ii = 1,min(nsend,nrecv)
                xx_recv(ii) = xx_send(ii)
             end do
          else
             call MPI_Sendrecv(                      &
                  xx_send(1:nsend), nsend4,          &
                  PAR_INTEGER, dom_i4, 0_4,          &
                  xx_recv(1:nrecv), nrecv4,          &
                  PAR_INTEGER, dom_i4, 0_4,          &
                  PAR_COMM_TO_USE, status, istat4    )             
          end if
       end if
    end if

    if( istat4 /= MPI_SUCCESS )call PAR_MPI_RUNEND(istat4,'PAR_SEND_RECEIVE_IP')
#endif

  end subroutine PAR_SEND_RECEIVE_IP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_SEND_RECEIVE_RP
  !
  !----------------------------------------------------------------------

  subroutine PAR_SEND_RECEIVE_RP_s(xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN)
    
    real(rp),            intent(in)           :: xx_send
    real(rp),            intent(out)          :: xx_recv
    character(*),        intent(in)           :: wherein
    integer(ip),         intent(in)           :: dom_i
    character(*),        intent(in), optional :: wsynch
    MY_MPI_COMM   ,      intent(in), optional :: PAR_COMM_IN
    real(rp)                                  :: yy_send(2)
    real(rp)                                  :: yy_recv(2)
    if( ISEQUEN ) return
    yy_send(1) = xx_send
    if( present(wsynch) ) then
       call PAR_SEND_RECEIVE_RP(1_ip,1_ip,yy_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
    else
       call PAR_SEND_RECEIVE_RP(1_ip,1_ip,yy_send,yy_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
    end if
    xx_recv    = yy_recv(1)
  end subroutine PAR_SEND_RECEIVE_RP_s
  subroutine PAR_SEND_RECEIVE_RP_0(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN)
    
    integer(ip),         intent(in)           :: nsend
    integer(ip),         intent(in)           :: nrecv
    real(rp),            intent(in)           :: xx_send(*)
    real(rp),            intent(out)          :: xx_recv(*)
    character(*),        intent(in)           :: wherein
    integer(ip),         intent(in)           :: dom_i
    character(*),        intent(in), optional :: wsynch
    MY_MPI_COMM   ,      intent(in), optional :: PAR_COMM_IN
    real(rp)                                  :: yy_send(2)
    real(rp)                                  :: yy_recv(2)
    if( ISEQUEN ) return
    if( present(wsynch) ) then
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else
          return
       end if
    else
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       else
          return
       end if
    end if
  end subroutine PAR_SEND_RECEIVE_RP_0
  subroutine PAR_SEND_RECEIVE_RP_1(xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN)
    
    real(rp),     pointer, intent(in)           :: xx_send(:)
    real(rp),     pointer, intent(inout)        :: xx_recv(:)
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in)           :: dom_i
    character(*),          intent(in), optional :: wsynch
    MY_MPI_COMM   ,        intent(in), optional :: PAR_COMM_IN
    integer(ip)                                 :: nsend,nrecv
    real(rp)                                    :: yy_send(2)
    real(rp)                                    :: yy_recv(2)
    if( ISEQUEN ) return
    if( .not. associated(xx_send) ) then
       nsend = 0
    else
       nsend = size(xx_send)
    end if
    if( .not. associated(xx_recv) ) then
       nrecv = 0
    else
       nrecv = size(xx_recv)
    end if
    if( present(wsynch) ) then
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else
          return
       end if
    else
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       else
          return
       end if
    end if
  end subroutine PAR_SEND_RECEIVE_RP_1
  subroutine PAR_SEND_RECEIVE_RP_2(xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN)
    
    real(rp),     pointer, intent(in)           :: xx_send(:,:)
    real(rp),     pointer, intent(inout)        :: xx_recv(:,:)
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in)           :: dom_i
    character(*),          intent(in), optional :: wsynch
    MY_MPI_COMM   ,        intent(in), optional :: PAR_COMM_IN
    integer(ip)                                 :: nsend,nrecv
    real(rp)                                    :: yy_send(2)
    real(rp)                                    :: yy_recv(2)
    if( ISEQUEN ) return
    if( .not. associated(xx_send) ) then
       nsend = 0
    else
       nsend = size(xx_send,1)*size(xx_send,2)
    end if
    if( .not. associated(xx_recv) ) then
       nrecv = 0
    else
       nrecv = size(xx_recv,1)*size(xx_recv,2)
    end if
    if( present(wsynch) ) then
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else
          return
       end if
    else
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       else
          return
       end if
    end if
  end subroutine PAR_SEND_RECEIVE_RP_2
  subroutine PAR_SEND_RECEIVE_RP_3(xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN)
    
    real(rp),     pointer, intent(in)           :: xx_send(:,:,:)
    real(rp),     pointer, intent(inout)        :: xx_recv(:,:,:)
    character(*),          intent(in)           :: wherein
    integer(ip),           intent(in)           :: dom_i
    character(*),          intent(in), optional :: wsynch
    MY_MPI_COMM   ,        intent(in), optional :: PAR_COMM_IN
    integer(ip)                                 :: nsend,nrecv
    real(rp)                                    :: yy_send(2)
    real(rp)                                    :: yy_recv(2)
    if( ISEQUEN ) return
    if( .not. associated(xx_send) ) then
       nsend = 0
    else
       nsend = size(xx_send,1)*size(xx_send,2)*size(xx_send,3)
    end if
    if( .not. associated(xx_recv) ) then
       nrecv = 0
    else
       nrecv = size(xx_recv,1)*size(xx_recv,2)*size(xx_recv,3)
    end if
    if( present(wsynch) ) then
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN=PAR_COMM_IN)
       else
          return
       end if
    else
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,PAR_COMM_IN=PAR_COMM_IN)
       else
          return
       end if
    end if
  end subroutine PAR_SEND_RECEIVE_RP_3

  !----------------------------------------------------------------------
  !
  ! PAR_SEND_RECEIVE_RP
  !
  !----------------------------------------------------------------------

  subroutine PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN)
    
    integer(ip),              intent(in)  :: nsend
    integer(ip),              intent(in)  :: nrecv
    real(rp),                 intent(in)  :: xx_send(*)
    real(rp),                 intent(out) :: xx_recv(*)
    character(*),   optional, intent(in)  :: wherein
    integer(ip),              intent(in)  :: dom_i
    character(*),   optional, intent(in)  :: wsynch
    MY_MPI_COMM   , optional, intent(in)  :: PAR_COMM_IN
    integer(ip)                           :: kk
    integer(4)                            :: istat4,nsend4,nrecv4,dom_i4,my_rank4
    MY_MPI_COMM                           :: PAR_COMM_TO_USE
    logical(lg)                           :: asynch

    if( IPARALL ) then
       !
       ! Define communicator
       !
       if( present(PAR_COMM_IN) ) then
          call PAR_COMM_RANK(PAR_COMM_TO_USE,MY_RANK4,PAR_COMM_IN)
       else
          call PAR_COMM_RANK(PAR_COMM_TO_USE,MY_RANK4,wherein)
       end if
       !
       ! BlockinG/non blocking
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
             if( .not. associated(non_blocking(inonblocking) % request4) ) then
                call runend('NON-BLOCKING SEND/RECEIVE SHOULD BE STARTED')
             end if
          else
             call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if
       !
       ! Synchronous Send/receive
       !
       nsend4 = int(nsend,4)
       nrecv4 = int(nrecv,4)
       dom_i4 = int(dom_i,4)
       istat4 = 0_4

#ifndef MPI_OFF
       if( asynch ) then
          if( nrecv /= 0 ) then
             non_blocking(inonblocking) % count4 = non_blocking(inonblocking) % count4 + 1
             kk = non_blocking(inonblocking) % count4
             call MPI_IRecv(                                 &
                  xx_recv(1:nrecv), nrecv4,                  &
                  PAR_REAL, dom_i4, 0_4,         &
                  PAR_COMM_TO_USE,                           &
                  non_blocking(inonblocking) % request4(kk), &
                  istat4                                     )
          end if
          if( nsend /= 0 ) then
             non_blocking(inonblocking) % count4 = non_blocking(inonblocking) % count4 + 1
             kk = non_blocking(inonblocking) % count4
             call MPI_ISend(                                 &
                  xx_send(1:nsend), nsend4,                  &
                  PAR_REAL, dom_i4, 0_4,         &
                  PAR_COMM_TO_USE,                           &
                  non_blocking(inonblocking) % request4(kk), &
                  istat4                                     )
          end if
       else
          if( nrecv /= 0 .and. nsend == 0 ) then
             call MPI_Recv(                          &
                  xx_recv(1:nrecv), nrecv4,          &
                  PAR_REAL, dom_i4, 0_4, &
                  PAR_COMM_TO_USE, status, istat4   )
          else if( nrecv == 0 .and. nsend /= 0 ) then
             call MPI_Send(                          &
                  xx_send(1:nsend), nsend4,          &
                  PAR_REAL, dom_i4, 0_4, &
                  PAR_COMM_TO_USE, istat4           )
          else if( nrecv /= 0 .and. nsend /= 0 ) then
             call MPI_Sendrecv(                      &
                  xx_send(1:nsend), nsend4,          &
                  PAR_REAL, dom_i4, 0_4, &
                  xx_recv(1:nrecv), nrecv4,          &
                  PAR_REAL, dom_i4, 0_4, &
                  PAR_COMM_TO_USE, status, istat4   )
          end if
       end if
       if( istat4 /= MPI_SUCCESS )call PAR_MPI_RUNEND(istat4,'PAR_SEND_RECEIVE_RP')
#endif

    end if

  end subroutine PAR_SEND_RECEIVE_RP

  !-----------------------------------------------------------------------
  !
  !> @brief   Send receive
  !> @details Send and receive to all my neghbors within communicator COMMU
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine PAR_SEND_RECEIVE_TO_ALL_IP_1c(xx_send,xx_recv,commu,wsynch,lsend_perm,PERMUTE_SEND,PERMUTE_RECV)

    integer(ip), pointer,           intent(in)    :: xx_send(:)
    integer(ip), pointer,           intent(inout) :: xx_recv(:)
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip), pointer, optional, intent(in)    :: lsend_perm(:)
    logical(lg),          optional, intent(in)    :: PERMUTE_SEND
    logical(lg),          optional, intent(in)    :: PERMUTE_RECV
    integer(ip)                                   :: ndofn
    integer(ip),                       pointer    :: xx_send_perm(:)
    integer(ip),                       pointer    :: xx_recv_perm(:)
    integer(ip)                                   :: icont,dom_i,ineig
    logical(lg)                                   :: if_permute_send
    logical(lg)                                   :: if_permute_recv
    integer(ip)                                   :: xx_send_tmp(2)
    integer(ip)                                   :: xx_recv_tmp(2)
    logical(lg)                                   :: alltoallv
    integer(ip)                                   :: PAR_CURRENT_RANK,PAR_CURRENT_SIZE
    integer(4),                        pointer    :: sendcount4(:)
    integer(4),                        pointer    :: recvcount4(:)

    ndofn = 1_ip
    if_permute_send = .false.
    if_permute_recv = .false.
    alltoallv       = .false.

    if( present(PERMUTE_SEND) .or. present(lsend_perm) ) if_permute_send = PERMUTE_SEND
    if( present(PERMUTE_RECV) ) if_permute_recv = PERMUTE_RECV
    nullify(xx_send_perm)
    nullify(xx_recv_perm)

    if( IPARALL ) then

       if( present(wsynch) ) then
          if( trim(wsynch) == 'ALLTOALLV' ) then
             alltoallv = .true.
          end if
       end if

       if( if_permute_send ) then
          call memory_alloca(par_memor,'xx_send_perm','par_send_receive_to_all_ip_1c',xx_send_perm,max(1_ip,commu % lsend_dim))
          do icont = 1,commu % lsend_dim
             xx_send_perm(icont) = xx_send(commu % lsend_perm(icont))
          end do

          if( alltoallv ) then

             call PAR_COMM_RANK_AND_SIZE(COMMU % PAR_COMM_WORLD,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
             allocate(sendcount4(0:PAR_CURRENT_SIZE-1))
             allocate(recvcount4(0:PAR_CURRENT_SIZE-1))
             sendcount4 = 0_4
             recvcount4 = 0_4
             do ineig = 1,commu % nneig
                dom_i             = commu % neights(ineig)
                sendcount4(dom_i) = commu % lsend_size(ineig+1) - commu % lsend_size(ineig)
                recvcount4(dom_i) = commu % lrecv_size(ineig+1) - commu % lrecv_size(ineig)
             end do
             call PAR_ALLTOALLV_IP_1(xx_send_perm,xx_recv,sendcount4,recvcount4,PAR_COMM_IN=commu % PAR_COMM_WORLD)
             deallocate(sendcount4)
             deallocate(recvcount4)

          else

             if( commu % lsend_dim > 0 ) then
                if( .not. associated(xx_recv) ) then
                   call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send_perm,xx_recv_tmp,commu,wsynch)
                else
                   call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send_perm,xx_recv,commu,wsynch)
                end if
             else
                if( .not. associated(xx_recv) ) then
                   call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send_tmp,xx_recv_tmp,commu,wsynch)
                else
                   call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send_tmp,xx_recv,commu,wsynch)
                end if
             end if
          end if
          call memory_deallo(par_memor,'xx_send_perm','par_send_receive_to_all_ip_1c',xx_send_perm)

       else
          if( .not. associated(xx_send) .and. .not. associated(xx_recv) ) then
             return
          else if( .not. associated(xx_send) ) then
             call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send_tmp,xx_recv,commu,wsynch)
          else if( .not. associated(xx_recv) ) then
             call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send,xx_recv_tmp,commu,wsynch)
          else
             call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send,xx_recv,commu,wsynch)
          end if
       endif

       if( if_permute_recv .and. commu % lrecv_dim > 0  ) then
          call runend('PAR_SEND_RECEIVE_TO_ALL_RP_0: NOT CODED')
       end if

    else

       if( if_permute_send .and. commu % lsend_dim > 0 ) then
          call runend('PAR_SEND_RECEIVE_TO_ALL: OPTION NOT CODED')
       else
          do icont = lbound(xx_recv,1),ubound(xx_recv,1)
             xx_recv(icont) = xx_send(icont)
          end do
       end if

    end if

  end subroutine PAR_SEND_RECEIVE_TO_ALL_IP_1c

  subroutine PAR_SEND_RECEIVE_TO_ALL_IP_2(xx_send,xx_recv,commu,wsynch,lsend_perm,PERMUTE_SEND,PERMUTE_RECV)
    
    integer(ip),          pointer,  intent(in)    :: xx_send(:,:)
    integer(ip),          pointer,  intent(inout) :: xx_recv(:,:)
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip), pointer, optional, intent(in)    :: lsend_perm(:)
    logical(lg),          optional                :: PERMUTE_SEND
    logical(lg),          optional                :: PERMUTE_RECV
    integer(ip)                                   :: yy_recv(2)
    integer(ip)                                   :: yy_send(2)
    integer(ip),             pointer              :: xx_send_perm(:,:)
    integer(ip),             pointer              :: xx_recv_perm(:,:)
    integer(ip)                                   :: ndofn,icont,dom_i,ineig
    logical(lg)                                   :: if_permute_send
    logical(lg)                                   :: if_permute_recv
    logical(lg)                                   :: alltoallv
    integer(ip)                                   :: PAR_CURRENT_RANK,PAR_CURRENT_SIZE
    integer(4),                        pointer    :: sendcount4(:)
    integer(4),                        pointer    :: recvcount4(:)

    if(      associated(xx_send) ) then
       ndofn = size(xx_send,1)
    else if( associated(xx_recv) ) then
       ndofn = size(xx_recv,1)
    else
       return
    end if

    if_permute_send = .false.
    if_permute_recv = .false.
    alltoallv       = .false.

    if( IPARALL ) then

       if( present(PERMUTE_SEND) .or. present(lsend_perm) ) if_permute_send = PERMUTE_SEND
       if( present(PERMUTE_RECV) ) if_permute_recv = PERMUTE_RECV
       nullify(xx_send_perm)
       nullify(xx_recv_perm)

       if( present(wsynch) ) then
          if( trim(wsynch) == 'ALLTOALLV' ) then
             alltoallv = .true.
          end if
       end if

       if( alltoallv ) then

          call PAR_COMM_RANK_AND_SIZE(COMMU % PAR_COMM_WORLD,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
          allocate(sendcount4(0:PAR_CURRENT_SIZE-1))
          allocate(recvcount4(0:PAR_CURRENT_SIZE-1))
          sendcount4 = 0_4
          recvcount4 = 0_4
          do ineig = 1,commu % nneig
             dom_i             = commu % neights(ineig)
             sendcount4(dom_i) = commu % lsend_size(ineig+1) - commu % lsend_size(ineig)
             recvcount4(dom_i) = commu % lrecv_size(ineig+1) - commu % lrecv_size(ineig)
          end do
          sendcount4 = sendcount4 * ndofn
          recvcount4 = recvcount4 * ndofn

       end if

       if( if_permute_recv .and. associated(xx_recv) ) then
          call memory_alloca(par_memor,'xx_recv_perm','par_send_receive_to_all_ip_1',xx_recv_perm,ndofn,max(1_ip,commu % lrecv_dim))
       else
          xx_recv_perm => xx_recv
       endif

       if( if_permute_send .and. associated(xx_send) ) then
          call memory_alloca(par_memor,'xx_send_perm','par_send_receive_to_all_ip_1',xx_send_perm,ndofn,max(1_ip,commu % lsend_dim))
          do icont= 1,commu % lsend_dim
             xx_send_perm(1:ndofn,icont) = xx_send(1:ndofn,commu % lsend_perm(icont))
          end do
          if( alltoallv ) then
             call PAR_ALLTOALLV_IP_2(xx_send_perm,xx_recv_perm,sendcount4,recvcount4,PAR_COMM_IN=commu % PAR_COMM_WORLD)
          else
             if( .not. associated(xx_recv) ) then
                call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send_perm,yy_recv,commu,wsynch)
             else
                call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send_perm,xx_recv_perm,commu,wsynch)
             end if
          end if
          call memory_deallo(par_memor,'xx_send_perm','par_send_receive_to_all_ip_01',xx_send_perm)
       else
          if( alltoallv ) then
             call PAR_ALLTOALLV_IP_2(xx_send,xx_recv_perm,sendcount4,recvcount4,PAR_COMM_IN=commu % PAR_COMM_WORLD)
          else
             if( .not. associated(xx_recv) ) then
                call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send,yy_recv,commu,wsynch)
             else if( .not. associated(xx_send) ) then
                call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,yy_send,xx_recv_perm,commu,wsynch)
             else
                call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send,xx_recv_perm,commu,wsynch)
             end if
          end if
       endif

       if( if_permute_recv .and. associated(xx_recv) ) then
          do icont= 1,commu % lrecv_dim
             xx_recv(1:ndofn,commu % lrecv_perm(icont)) = xx_recv_perm(1:ndofn,icont)
          end do
       end if

       if( alltoallv ) then
          deallocate(sendcount4)
          deallocate(recvcount4)
       end if

    else

       if( if_permute_send .and. associated(xx_send) ) then
          call runend('PAR_SEND_RECEIVE_TO_ALL_IP_2')
       else
          xx_recv = xx_send
       end if

    end if

  end subroutine PAR_SEND_RECEIVE_TO_ALL_IP_2

  subroutine PAR_SEND_RECEIVE_TO_ALL_IP_0(ndofn,xx_send,xx_recv,commu,wsynch,lsend_perm,PERMUTE_SEND,PERMUTE_RECV)
    
    integer(ip),                    intent(in)  :: ndofn
    integer(ip),                    intent(in)  :: xx_send(ndofn,*)
    integer(ip),                    intent(out) :: xx_recv(ndofn,*)
    type(comm_data_par),            intent(in)  :: commu
    character(*),         optional, intent(in)  :: wsynch
    integer(ip), pointer, optional, intent(in)  :: lsend_perm(:)
    logical(lg),          optional              :: PERMUTE_SEND
    logical(lg),          optional              :: PERMUTE_RECV
    integer(ip),                       pointer  :: xx_send_perm(:,:)
    integer(ip),                       pointer  :: xx_recv_perm(:,:)
    integer(ip)                                 :: icont
    logical(lg)                                 :: if_permute_send
    logical(lg)                                 :: if_permute_recv

    if( ISEQUEN ) return

    if_permute_send = .false.
    if_permute_recv = .false.
    if( present(PERMUTE_SEND) .or. present(lsend_perm) ) if_permute_send = PERMUTE_SEND
    if( present(PERMUTE_RECV) ) if_permute_recv = PERMUTE_RECV
    nullify(xx_send_perm)
    nullify(xx_recv_perm)

    if( if_permute_send ) then
       nullify(xx_send_perm)
       call memory_alloca(par_memor,'xx_send_perm','par_send_receive_to_all_ip_0',xx_send_perm,ndofn,max(1_ip,commu % lsend_dim))
       do icont=1, commu % lsend_dim
          xx_send_perm(1:ndofn,icont) = xx_send(1:ndofn,commu % lsend_perm(icont))
       end do
       call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send_perm,xx_recv,commu,wsynch)
       call memory_deallo(par_memor,'xx_send_perm','par_send_receive_to_all_ip_0',xx_send_perm)
    else
       call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send,xx_recv,commu,wsynch)
    end if

  end subroutine PAR_SEND_RECEIVE_TO_ALL_IP_0

  subroutine PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send,xx_recv,commu,wsynch)
    
    integer(ip),                    intent(in)  :: ndofn
    integer(ip),                    intent(in)  :: xx_send(*)
    integer(ip),                    intent(out) :: xx_recv(*)
    type(comm_data_par),            intent(in)  :: commu
    character(*),         optional, intent(in)  :: wsynch
    integer(ip)                                 :: inise,finse,inire,finre,counti
    integer(ip)                                 :: dom_i,ineig,nsend,nrecv,kk
    integer(4)                                  :: istat4,nsend4,nrecv4,dom_i4
    integer(4)                                  :: count4
    MY_MPI_COMM                                 :: PAR_COMM_TO_USE
    logical(lg)                                 :: asynch
    MY_MPI_REQUEST,     allocatable             :: ireqq4(:)
    MY_MPI_STATUS ,     pointer                 :: status4(:)
    !
    ! Define communicator
    !
    PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
    !
    ! Asynchronous
    !
    if( present(wsynch) ) then
       if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
          asynch = .false.
       else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
          asynch = .true.
       else
          call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
       end if
    else
       asynch = .false.
    end if
    if( asynch ) allocate(ireqq4(commu % nneig*2))
    kk = 0
    !
    ! Synchronous Send/receive
    !
    do ineig = 1,commu % nneig

       dom_i  = commu % neights(ineig)

       inise  = ndofn * ( commu % lsend_size(ineig)   - 1 ) + 1
       nsend  = ndofn * ( commu % lsend_size(ineig+1) - 1 ) + 1 - inise
       finse  = inise + nsend - 1
       inire  = ndofn * ( commu % lrecv_size(ineig)   - 1 ) + 1
       nrecv  = ndofn * ( commu % lrecv_size(ineig+1) - 1 ) + 1 - inire
       finre  = inire + nrecv - 1

       nsend4 = int(nsend,4)
       nrecv4 = int(nrecv,4)
       dom_i4 = int(dom_i,4)

#ifndef MPI_OFF
       if( asynch .and. commu % nneig /= 0_ip ) then
          if( nsend > 0_ip ) then
             kk = kk + 1
             call MPI_Isend(                          &
                  xx_send(inise:finse), nsend4,       &
                  PAR_INTEGER,  dom_i4, 0_4,          &
                  PAR_COMM_TO_USE, ireqq4(kk), istat4 )
          end if
          if( nrecv > 0_ip ) then
             kk = kk + 1
             call MPI_Irecv(                          &
                  xx_recv(inire:finre), nrecv4,       &
                  PAR_INTEGER,  dom_i4, 0_4,          &
                  PAR_COMM_TO_USE, ireqq4(kk), istat4 )
          end if
       else
          if( nrecv /= 0 .and. nsend == 0 ) then
             call MPI_Recv(                       &
                  xx_recv(inire:finre), nrecv4,   &
                  PAR_INTEGER, dom_i4, 0_4,       &
                  PAR_COMM_TO_USE, status, istat4 )

          else if( nrecv == 0 .and. nsend /= 0 ) then
             call MPI_Send(                       &
                  xx_send(inise:finse), nsend4,   &
                  PAR_INTEGER, dom_i4, 0_4,       &
                  PAR_COMM_TO_USE, istat4         )

          else if( nrecv /= 0 .and. nsend /= 0 ) then
             call MPI_Sendrecv(                   &
                  xx_send(inise:finse), nsend4,   &
                  PAR_INTEGER, dom_i4, 0_4,       &
                  xx_recv(inire:finre), nrecv4,   &
                  PAR_INTEGER, dom_i4, 0_4,       &
                  PAR_COMM_TO_USE, status, istat4 )
          end if
       end if
#endif

    end do
    !
    ! Wait in case of asynchronous, count is the total of communications actually performe
    ! and is smaller or equal than the asynchronous communications requested (ireqq)
    !
#ifndef MPI_OFF
    if( asynch .and. commu % nneig /= 0 ) then
       counti = kk
       count4 = int(counti,4)
       allocate( status4(MPI_STATUS_SIZE*counti) )
       CALL MPI_WAITALL(count4,ireqq4,status4,istat4)
       deallocate( status4 )
       deallocate( ireqq4  )
    end if
#endif

  end subroutine PAR_SEND_RECEIVE_TO_ALL_IP

  subroutine PAR_SEND_RECEIVE_TO_ALL_RP_s(xx_send,xx_recv,commu,wsynch)
    
    real(rp),                       intent(in)  :: xx_send
    real(rp),                       intent(out) :: xx_recv
    type(comm_data_par),            intent(in)  :: commu
    character(*),         optional, intent(in)  :: wsynch
    real(rp)                                    :: yy_send(2)
    real(rp)                                    :: yy_recv(2)
    if( ISEQUEN ) return
    yy_send(1) = xx_send
    call PAR_SEND_RECEIVE_TO_ALL_RP(1_ip,yy_send,yy_recv,commu,wsynch)
    xx_recv    = yy_recv(1)
  end subroutine PAR_SEND_RECEIVE_TO_ALL_RP_s

 subroutine PAR_SEND_RECEIVE_TO_ALL_RP_0(ndofn,xx_send,xx_recv,commu,wsynch,lsend_perm,PERMUTE_SEND,PERMUTE_RECV)
    
    integer(ip),                    intent(in)  :: ndofn
    real(rp),                       intent(in)  :: xx_send(ndofn,*)
    real(rp),                       intent(out) :: xx_recv(ndofn,*)
    type(comm_data_par),            intent(in)  :: commu
    integer(ip), pointer, optional, intent(in)  :: lsend_perm(:)
    logical(lg),          optional, intent(in)  :: PERMUTE_SEND
    logical(lg),          optional, intent(in)  :: PERMUTE_RECV
    character(*),         optional, intent(in)  :: wsynch
    real(rp),                          pointer  :: xx_send_perm(:,:)
    real(rp),                          pointer  :: xx_recv_perm(:,:)
    integer(ip)                                 :: icont
    logical(lg)                                 :: if_permute_send
    logical(lg)                                 :: if_permute_recv

    if( ISEQUEN ) return

    if_permute_send = .false.
    if_permute_recv = .false.
    if( present(PERMUTE_SEND) .or. present(lsend_perm) ) if_permute_send = PERMUTE_SEND
    if( present(PERMUTE_RECV) ) if_permute_recv = PERMUTE_RECV
    nullify(xx_send_perm)
    nullify(xx_recv_perm)

    if( if_permute_send .and. commu % lsend_dim > 0 ) then
       call memory_alloca(par_memor,'xx_send_perm','par_send_receive_to_all_rp_0',xx_send_perm,ndofn,max(1_ip,commu % lsend_dim))
       do icont = 1,commu % lsend_dim
          xx_send_perm(1:ndofn,icont) = xx_send(1:ndofn,commu % lsend_perm(icont))
       end do
       call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send_perm,xx_recv,commu,wsynch)
    else
       call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send,xx_recv,commu,wsynch)
    endif

    if( if_permute_send .and. commu % lsend_dim > 0 ) then
       call memory_deallo(par_memor,'xx_send_perm','par_send_receive_to_all_rp_0',xx_send_perm)
    end if

    if( if_permute_recv .and. commu % lrecv_dim > 0  ) then
       call runend('PAR_SEND_RECEIVE_TO_ALL_RP_0: NOT CODED')
    end if

 end subroutine PAR_SEND_RECEIVE_TO_ALL_RP_0

  subroutine PAR_SEND_RECEIVE_TO_ALL_RP_1(xx_send,xx_recv,wherein,wsynch,lsend_perm,PERMUTE_SEND,PERMUTE_RECV)
    
    real(rp),    pointer,           intent(in)    :: xx_send(:)
    real(rp),    pointer,           intent(inout) :: xx_recv(:)
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip), pointer, optional, intent(in)    :: lsend_perm(:)
    logical(lg),          optional, intent(in)    :: PERMUTE_SEND
    logical(lg),          optional, intent(in)    :: PERMUTE_RECV

    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    real(rp),                          pointer    :: xx_send_perm(:)
    real(rp),                          pointer    :: xx_recv_perm(:)
    integer(ip)                                   :: icont
    logical(lg)                                   :: if_permute_send
    logical(lg)                                   :: if_permute_recv
    real(rp)                                      :: xx_send_tmp(2)
    real(rp)                                      :: xx_recv_tmp(2)

    if( ISEQUEN ) return
    ndofn = 1_ip
    if_permute_send = .false.
    if_permute_recv = .false.
    if( present(PERMUTE_SEND) .or. present(lsend_perm) ) if_permute_send = PERMUTE_SEND
    if( present(PERMUTE_RECV) ) if_permute_recv = PERMUTE_RECV
    nullify(xx_send_perm)
    nullify(xx_recv_perm)

    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)

    if( if_permute_send .and. commu % lsend_dim > 0 ) then
       call memory_alloca(par_memor,'xx_send_perm','par_send_receive_to_all_rp_1',xx_send_perm,max(1_ip,commu % lsend_dim))
       do icont = 1,commu % lsend_dim
          xx_send_perm(icont) = xx_send(commu % lsend_perm(icont))
       end do
       if( .not. associated(xx_recv) ) then
          call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send_perm,xx_recv_tmp,commu,wsynch)
       else
          call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send_perm,xx_recv,commu,wsynch)
       end if
    else
      if( .not. associated(xx_send) .and. .not. associated(xx_recv) ) then
          return
       else if( .not. associated(xx_send) ) then
          call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send_tmp,xx_recv,commu,wsynch)
       else if( .not. associated(xx_recv) ) then
          call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send,xx_recv_tmp,commu,wsynch)
       else
          call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send,xx_recv,commu,wsynch)
       end if
    endif

    if( if_permute_send .and. commu % lsend_dim > 0 ) then
       call memory_deallo(par_memor,'xx_send_perm','par_send_receive_to_all_rp_1',xx_send_perm)
    end if

    if( if_permute_recv .and. commu % lrecv_dim > 0  ) then
       call runend('PAR_SEND_RECEIVE_TO_ALL_RP_0: NOT CODED')
    end if

  end subroutine PAR_SEND_RECEIVE_TO_ALL_RP_1

  subroutine PAR_SEND_RECEIVE_TO_ALL_RP_1c(xx_send,xx_recv,commu,wsynch,lsend_perm,PERMUTE_SEND,PERMUTE_RECV)

    real(rp),    pointer,           intent(in)    :: xx_send(:)
    real(rp),    pointer,           intent(inout) :: xx_recv(:)
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip), pointer, optional, intent(in)    :: lsend_perm(:)
    logical(lg),          optional, intent(in)    :: PERMUTE_SEND
    logical(lg),          optional, intent(in)    :: PERMUTE_RECV
    integer(ip)                                   :: ndofn
    real(rp),                          pointer    :: xx_send_perm(:)
    real(rp),                          pointer    :: xx_recv_perm(:)
    integer(ip)                                   :: icont,dom_i,ineig
    logical(lg)                                   :: if_permute_send
    logical(lg)                                   :: if_permute_recv
    real(rp)                                      :: xx_send_tmp(2)
    real(rp)                                      :: xx_recv_tmp(2)
    logical(lg)                                   :: alltoallv
    integer(ip)                                   :: PAR_CURRENT_RANK,PAR_CURRENT_SIZE
    integer(4),                        pointer    :: sendcount4(:)
    integer(4),                        pointer    :: recvcount4(:)

    if( ISEQUEN ) return
    ndofn = 1_ip
    if_permute_send = .false.
    if_permute_recv = .false.
    alltoallv       = .false.

    if( present(PERMUTE_SEND) .or. present(lsend_perm) ) if_permute_send = PERMUTE_SEND
    if( present(PERMUTE_RECV) ) if_permute_recv = PERMUTE_RECV
    nullify(xx_send_perm)
    nullify(xx_recv_perm)

    if( present(wsynch) ) then
       if( trim(wsynch) == 'ALLTOALLV' ) then
          alltoallv = .true.
       end if
    end if

    if( alltoallv ) then

       call PAR_COMM_RANK_AND_SIZE(COMMU % PAR_COMM_WORLD,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
       allocate(sendcount4(0:PAR_CURRENT_SIZE-1))
       allocate(recvcount4(0:PAR_CURRENT_SIZE-1))
       sendcount4 = 0_4
       recvcount4 = 0_4
       do ineig = 1,commu % nneig
          dom_i             = commu % neights(ineig)
          sendcount4(dom_i) = commu % lsend_size(ineig+1) - commu % lsend_size(ineig)
          recvcount4(dom_i) = commu % lrecv_size(ineig+1) - commu % lrecv_size(ineig)
       end do

    end if

    if( if_permute_send .and. commu % lsend_dim > 0 ) then
       call memory_alloca(par_memor,'xx_send_perm','par_send_receive_to_all_rp_1c',xx_send_perm,max(1_ip,commu % lsend_dim))
       do icont = 1,commu % lsend_dim
          xx_send_perm(icont) = xx_send(commu % lsend_perm(icont))
       end do
       if( alltoallv ) then
          call PAR_ALLTOALLV_RP_1(xx_send_perm,xx_recv,sendcount4,recvcount4,PAR_COMM_IN=commu % PAR_COMM_WORLD)
       else
          if( .not. associated(xx_recv) ) then
             call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send_perm,xx_recv_tmp,commu,wsynch)
          else
             call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send_perm,xx_recv,commu,wsynch)
          end if
       end if
    else
       if( .not. associated(xx_send) .and. .not. associated(xx_recv) ) then
          return
       else if( .not. associated(xx_send) ) then
          call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send_tmp,xx_recv,commu,wsynch)
       else if( .not. associated(xx_recv) ) then
          call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send,xx_recv_tmp,commu,wsynch)
       else
          call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send,xx_recv,commu,wsynch)
       end if
    endif

    if( if_permute_send .and. commu % lsend_dim > 0 ) then
       call memory_deallo(par_memor,'xx_send_perm','par_send_receive_to_all_rp_1c',xx_send_perm)
    end if

    if( if_permute_recv .and. commu % lrecv_dim > 0  ) then
       call runend('PAR_SEND_RECEIVE_TO_ALL_RP_0: NOT CODED')
    end if

    if( alltoallv ) then
       deallocate(sendcount4)
       deallocate(recvcount4)
    end if

  end subroutine PAR_SEND_RECEIVE_TO_ALL_RP_1c

  subroutine PAR_SEND_RECEIVE_TO_ALL_RP_2(xx_send,xx_recv,commu,wsynch,lsend_perm,PERMUTE_SEND,PERMUTE_RECV)
    
    real(rp),             pointer,  intent(in)    :: xx_send(:,:)
    real(rp),             pointer,  intent(inout) :: xx_recv(:,:)
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip), pointer, optional, intent(in)    :: lsend_perm(:)
    logical(lg),          optional                :: PERMUTE_SEND
    logical(lg),          optional                :: PERMUTE_RECV
    integer(ip)                                   :: ndofn,icont,dom_i,ineig
    real(rp)                                      :: yy_send(2)
    real(rp)                                      :: yy_recv(2)
    real(rp),             pointer                 :: xx_send_perm(:,:)
    real(rp),             pointer                 :: xx_recv_perm(:,:)
    logical(lg)                                   :: if_permute_send
    logical(lg)                                   :: if_permute_recv
    logical(lg)                                   :: alltoallv
    integer(ip)                                   :: PAR_CURRENT_RANK,PAR_CURRENT_SIZE
    integer(4),                        pointer    :: sendcount4(:)
    integer(4),                        pointer    :: recvcount4(:)

    if_permute_send = .false.
    if_permute_recv = .false.
    alltoallv       = .false.

    if( present(PERMUTE_SEND) .or. present(lsend_perm) ) if_permute_send = PERMUTE_SEND
    if( present(PERMUTE_RECV) ) if_permute_recv = PERMUTE_RECV
    nullify(xx_send_perm)
    nullify(xx_recv_perm)

    if( present(wsynch) ) then
       if( trim(wsynch) == 'ALLTOALLV' ) then
          alltoallv = .true.
       end if
    end if

    if( if_permute_send .and. associated(xx_send) ) then
       ndofn = size(xx_send,1)
       call memory_alloca(par_memor,'xx_send_perm','par_send_receive_to_all_rp_2',xx_send_perm,ndofn,max(1_ip,commu % lsend_dim))
       do icont = 1,commu % lsend_dim
          xx_send_perm(1:ndofn,icont) = xx_send(1:ndofn,commu % lsend_perm(icont))
       end do
    else
       xx_send_perm => xx_send
    end if

    if( if_permute_recv ) then
       call memory_alloca(par_memor,'xx_recv_perm','par_send_receive_to_all_rp_0',xx_recv_perm,ndofn,max(1_ip,commu % lrecv_dim))
   else
       xx_recv_perm => xx_recv
    end if

    if( .not. if_permute_send ) then
       if( associated(xx_send_perm) ) then
          if( size(xx_send_perm,2) < npoin ) call runend('WRONG SIZE: CANNOT INTERPOLATE')
       end if
    end if

    if( alltoallv ) then

       call PAR_COMM_RANK_AND_SIZE(COMMU % PAR_COMM_WORLD,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
       allocate(sendcount4(0:PAR_CURRENT_SIZE-1))
       allocate(recvcount4(0:PAR_CURRENT_SIZE-1))
       sendcount4 = 0_4
       recvcount4 = 0_4
       do ineig = 1,commu % nneig
          dom_i             = commu % neights(ineig)
          sendcount4(dom_i) = commu % lsend_size(ineig+1) - commu % lsend_size(ineig)
          recvcount4(dom_i) = commu % lrecv_size(ineig+1) - commu % lrecv_size(ineig)
       end do

    end if

    if( .not. associated(xx_send_perm) .and. associated(xx_recv_perm) ) then
       ndofn = size(xx_recv_perm,1)
       if( alltoallv ) then
          sendcount4 = sendcount4 * ndofn
          recvcount4 = recvcount4 * ndofn
          call PAR_ALLTOALLV_RP_2(xx_send_perm,xx_recv_perm,sendcount4,recvcount4,PAR_COMM_IN=commu % PAR_COMM_WORLD)
       else
          call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,yy_send,xx_recv_perm,commu,wsynch)
       end if

    else if( .not. associated(xx_send_perm) .and. .not. associated(xx_recv_perm) ) then
       continue

    else if( associated(xx_send_perm) .and. .not. associated(xx_recv_perm) ) then
       ndofn = size(xx_send_perm,1)
       if( alltoallv ) then
          sendcount4 = sendcount4 * ndofn
          recvcount4 = recvcount4 * ndofn
          call PAR_ALLTOALLV_RP_2(xx_send_perm,xx_recv_perm,sendcount4,recvcount4,PAR_COMM_IN=commu % PAR_COMM_WORLD)
       else
          call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send_perm,yy_recv,commu,wsynch)
       end if

    else
       ndofn = size(xx_send_perm,1)
       if( ndofn /= size(xx_recv_perm,1) ) call runend('WRONG SIZE: CANNOT INTERPOLATE')
       if( alltoallv ) then
          sendcount4 = sendcount4 * ndofn
          recvcount4 = recvcount4 * ndofn
          call PAR_ALLTOALLV_RP_2(xx_send_perm,xx_recv_perm,sendcount4,recvcount4,PAR_COMM_IN=commu % PAR_COMM_WORLD)
       else
          call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send_perm,xx_recv_perm,commu,wsynch)
       end if

    end if

    call memory_deallo(par_memor,'xx_send_perm','par_send_receive_to_all_rp_2',xx_send_perm)

   if( if_permute_recv .and. associated(xx_recv) ) then
      do icont = 1,commu % lrecv_dim
         xx_recv(1:ndofn,commu % lrecv_perm(icont)) = xx_recv_perm(1:ndofn,icont)
      end do
      call memory_deallo(par_memor,'xx_recv_perm','par_send_receive_to_all_rp_2',xx_recv_perm)
   end if

   if( alltoallv ) then
      deallocate(sendcount4)
      deallocate(recvcount4)
   end if

  end subroutine PAR_SEND_RECEIVE_TO_ALL_RP_2

  subroutine PAR_SEND_RECEIVE_TO_ALL_RP_3(xx_send,xx_recv,commu,wsynch,PERMUTE_SEND)
    
    real(rp),             pointer,  intent(in)    :: xx_send(:,:,:)
    real(rp),             pointer,  intent(inout) :: xx_recv(:,:,:)
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    logical(lg),          optional                :: PERMUTE_SEND
    integer(ip)                                   :: ndofn,ndof1,ndof2,icont,dom_i,ineig
    real(rp)                                      :: yy_send(2)
    real(rp)                                      :: yy_recv(2)
    real(rp),             pointer                 :: xx_send_perm(:,:,:)
    logical(lg)                                   :: if_permute_send
    logical(lg)                                   :: alltoallv
    integer(ip)                                   :: PAR_CURRENT_RANK,PAR_CURRENT_SIZE
    integer(4),                        pointer    :: sendcount4(:)
    integer(4),                        pointer    :: recvcount4(:)

    if_permute_send = .false.
    if( present(PERMUTE_SEND) ) if_permute_send = PERMUTE_SEND
    alltoallv = .false.

    if( present(wsynch) ) then
       if( trim(wsynch) == 'ALLTOALLV' ) then
          alltoallv = .true.
       end if
    end if

    if( alltoallv ) then

       call PAR_COMM_RANK_AND_SIZE(COMMU % PAR_COMM_WORLD,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
       allocate(sendcount4(0:PAR_CURRENT_SIZE-1))
       allocate(recvcount4(0:PAR_CURRENT_SIZE-1))
       sendcount4 = 0_4
       recvcount4 = 0_4
       do ineig = 1,commu % nneig
          dom_i             = commu % neights(ineig)
          sendcount4(dom_i) = commu % lsend_size(ineig+1) - commu % lsend_size(ineig)
          recvcount4(dom_i) = commu % lrecv_size(ineig+1) - commu % lrecv_size(ineig)
       end do

    end if


    if( if_permute_send ) then

       nullify(xx_send_perm)
       if(      associated(xx_send) ) then
          ndof1 = size(xx_send,1)
          ndof2 = size(xx_send,2)
       else if( associated(xx_recv) ) then
          ndof1 = size(xx_recv,1)
          ndof2 = size(xx_recv,2)
       else
          return
       end if
       ndofn = ndof1*ndof2

       call memory_alloca(par_memor,'xx_send_perm','par_send_receive_to_all_rp_3',xx_send_perm,ndof1,ndof2,max(1_ip,commu % lsend_dim))
        do icont = 1, commu % lsend_dim
          xx_send_perm(1:ndof1,1:ndof2,icont) = xx_send(1:ndof1,1:ndof2,commu % lsend_perm(icont))
       end do

       if( alltoallv ) then
          sendcount4 = sendcount4 * ndofn
          recvcount4 = recvcount4 * ndofn
          call PAR_ALLTOALLV_RP_3(xx_send_perm,xx_recv,sendcount4,recvcount4,PAR_COMM_IN=commu % PAR_COMM_WORLD)
       else
          if( .not. associated(xx_recv) ) then
             call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send_perm,yy_recv,commu,wsynch)
          else
             call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send_perm,xx_recv,commu,wsynch)
          end if
       end if
       call memory_deallo(par_memor,'xx_send_perm','par_send_receive_to_all_rp_3',xx_send_perm)

    else

       if( .not. associated(xx_send) .and. associated(xx_recv) ) then
          ndofn = size(xx_recv,1)*size(xx_recv,2)
          if( alltoallv ) then
             sendcount4 = sendcount4 * ndofn
             recvcount4 = recvcount4 * ndofn
             call PAR_ALLTOALLV_RP_3(xx_send,xx_recv,sendcount4,recvcount4,PAR_COMM_IN=commu % PAR_COMM_WORLD)
          else
             call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,yy_send,xx_recv,commu,wsynch)
          end if
       else if( .not. associated(xx_send) .and. .not. associated(xx_recv) ) then
          ndofn = 1
          if( alltoallv ) then
             sendcount4 = sendcount4
             recvcount4 = recvcount4
             call PAR_ALLTOALLV_RP_3(xx_send,xx_recv,sendcount4,recvcount4,PAR_COMM_IN=commu % PAR_COMM_WORLD)
          else
             call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,yy_send,yy_recv,commu,wsynch)
          end if
       else if( associated(xx_send) .and. .not. associated(xx_recv) ) then
          ndofn = size(xx_send,1)* size(xx_send,2)
          if( alltoallv ) then
             sendcount4 = sendcount4 * ndofn
             recvcount4 = recvcount4 * ndofn
             call PAR_ALLTOALLV_RP_3(xx_send,xx_recv,sendcount4,recvcount4,PAR_COMM_IN=commu % PAR_COMM_WORLD)
          else
             call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send,yy_recv,commu,wsynch)
          end if
       else
          ndofn = size(xx_send,1)*size(xx_send,2)
          if( ndofn /= size(xx_recv,1)*size(xx_recv,2) ) call runend('WRONG SIZE: CANNOT INTERPOLATE')
          if( alltoallv ) then
             sendcount4 = sendcount4 * ndofn
             recvcount4 = recvcount4 * ndofn
             call PAR_ALLTOALLV_RP_3(xx_send,xx_recv,sendcount4,recvcount4,PAR_COMM_IN=commu % PAR_COMM_WORLD)
          else
             call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send,xx_recv,commu,wsynch)
          end if
       end if
    end if

  end subroutine PAR_SEND_RECEIVE_TO_ALL_RP_3

  subroutine PAR_SEND_RECEIVE_TO_ALL_RP_4(xx_send,xx_recv,commu,wsynch,PERMUTE_SEND)
    
    real(rp),             pointer,  intent(in)    :: xx_send(:,:,:,:)
    real(rp),             pointer,  intent(inout) :: xx_recv(:,:,:,:)
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    logical(lg),          optional                :: PERMUTE_SEND
    integer(ip)                                   :: ndofn,ndof1,ndof2,ndof3,icont,dom_i,ineig
    real(rp)                                      :: yy_send(2)
    real(rp)                                      :: yy_recv(2)
    real(rp),             pointer                 :: xx_send_perm(:,:,:,:)
    logical(lg)                                   :: if_permute_send
    logical(lg)                                   :: alltoallv
    integer(ip)                                   :: PAR_CURRENT_RANK,PAR_CURRENT_SIZE
    integer(4),                        pointer    :: sendcount4(:)
    integer(4),                        pointer    :: recvcount4(:)

    if_permute_send = .false.
    if( present(PERMUTE_SEND) ) if_permute_send = PERMUTE_SEND
    alltoallv = .false.

    if( present(wsynch) ) then
       if( trim(wsynch) == 'ALLTOALLV' ) then
          alltoallv = .true.
       end if
    end if

    if( alltoallv ) then

       call PAR_COMM_RANK_AND_SIZE(COMMU % PAR_COMM_WORLD,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
       allocate(sendcount4(0:PAR_CURRENT_SIZE-1))
       allocate(recvcount4(0:PAR_CURRENT_SIZE-1))
       sendcount4 = 0_4
       recvcount4 = 0_4
       do ineig = 1,commu % nneig
          dom_i             = commu % neights(ineig)
          sendcount4(dom_i) = commu % lsend_size(ineig+1) - commu % lsend_size(ineig)
          recvcount4(dom_i) = commu % lrecv_size(ineig+1) - commu % lrecv_size(ineig)
       end do

    end if


    if( if_permute_send ) then

       nullify(xx_send_perm)
       if(      associated(xx_send) ) then
          ndof1 = size(xx_send,1)
          ndof2 = size(xx_send,2)
          ndof3 = size(xx_send,3)
       else if( associated(xx_recv) ) then
          ndof1 = size(xx_recv,1)
          ndof2 = size(xx_recv,2)
          ndof3 = size(xx_recv,3)
       else
          return
       end if
       ndofn = ndof1*ndof2*ndof3

       call memory_alloca(par_memor,'xx_send_perm','par_send_receive_to_all_rp_4',xx_send_perm,ndof1,ndof2,ndof3,max(1_ip,commu % lsend_dim))
        do icont = 1, commu % lsend_dim
          xx_send_perm(1:ndof1,1:ndof2,1:ndof3,icont) = xx_send(1:ndof1,1:ndof2,1:ndof3,commu % lsend_perm(icont))
       end do

       if( alltoallv ) then
          sendcount4 = sendcount4 * ndofn
          recvcount4 = recvcount4 * ndofn
          call PAR_ALLTOALLV_RP_4(xx_send_perm,xx_recv,sendcount4,recvcount4,PAR_COMM_IN=commu % PAR_COMM_WORLD)
       else
          if( .not. associated(xx_recv) ) then
             call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send_perm,yy_recv,commu,wsynch)
          else
             call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send_perm,xx_recv,commu,wsynch)
          end if
       end if
       call memory_deallo(par_memor,'xx_send_perm','par_send_receive_to_all_rp_4',xx_send_perm)

    else

       if( .not. associated(xx_send) .and. associated(xx_recv) ) then
          ndofn = size(xx_recv,1)*size(xx_recv,2)*size(xx_recv,3)
          if( alltoallv ) then
             sendcount4 = sendcount4 * ndofn
             recvcount4 = recvcount4 * ndofn
             call PAR_ALLTOALLV_RP_4(xx_send,xx_recv,sendcount4,recvcount4,PAR_COMM_IN=commu % PAR_COMM_WORLD)
          else
             call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,yy_send,xx_recv,commu,wsynch)
          end if
       else if( .not. associated(xx_send) .and. .not. associated(xx_recv) ) then
          ndofn = 1
          if( alltoallv ) then
             sendcount4 = sendcount4
             recvcount4 = recvcount4
             call PAR_ALLTOALLV_RP_4(xx_send,xx_recv,sendcount4,recvcount4,PAR_COMM_IN=commu % PAR_COMM_WORLD)
          else
             call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,yy_send,yy_recv,commu,wsynch)
          end if
       else if( associated(xx_send) .and. .not. associated(xx_recv) ) then
          ndofn = size(xx_send,1)*size(xx_send,2)*size(xx_send,3)
          if( alltoallv ) then
             sendcount4 = sendcount4 * ndofn
             recvcount4 = recvcount4 * ndofn
             call PAR_ALLTOALLV_RP_4(xx_send,xx_recv,sendcount4,recvcount4,PAR_COMM_IN=commu % PAR_COMM_WORLD)
          else
             call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send,yy_recv,commu,wsynch)
          end if
       else
          ndofn = size(xx_send,1)*size(xx_send,2)*size(xx_send,3)
          if( ndofn /= size(xx_recv,1)*size(xx_recv,2)*size(xx_recv,3) ) call runend('WRONG SIZE: CANNOT INTERPOLATE')
          if( alltoallv ) then
             sendcount4 = sendcount4 * ndofn
             recvcount4 = recvcount4 * ndofn
             call PAR_ALLTOALLV_RP_4(xx_send,xx_recv,sendcount4,recvcount4,PAR_COMM_IN=commu % PAR_COMM_WORLD)
          else
             call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send,xx_recv,commu,wsynch)
          end if
       end if
    end if

  end subroutine PAR_SEND_RECEIVE_TO_ALL_RP_4

  subroutine PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send,xx_recv,commu,wsynch)
    
    integer(ip),                    intent(in)  :: ndofn
    real(rp),                       intent(in)  :: xx_send(*)
    real(rp),                       intent(out) :: xx_recv(*)
    type(comm_data_par),            intent(in)  :: commu
    character(*),         optional, intent(in)  :: wsynch
    integer(ip)                                 :: inise,finse,inire,finre,counti
    integer(ip)                                 :: dom_i,ineig,nsend,nrecv,kk
    integer(4)                                  :: istat4,nsend4,nrecv4,dom_i4
    integer(4)                                  :: count4
    MY_MPI_COMM                                 :: PAR_COMM_TO_USE
    logical(lg)                                 :: asynch
    MY_MPI_REQUEST,     allocatable             :: ireqq4(:)
    MY_MPI_STATUS ,     pointer                 :: status4(:)
    !
    ! Define communicator
    !
    PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
    !
    ! Asynchronous
    !
    if( present(wsynch) ) then
       if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
          asynch = .false.
       else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
          asynch = .true.
       else
          call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
       end if
    else
       asynch = .false.
    end if
    if( asynch .and. commu % nneig /= 0_ip ) allocate(ireqq4(commu % nneig*2_ip))

    kk = 0
    istat4 = 0_4
    !
    ! Synchronous Send/receive
    !
    do ineig = 1,commu % nneig

       dom_i  = commu % neights(ineig)
       inise  = ndofn * ( commu % lsend_size(ineig)   - 1 ) + 1
       nsend  = ndofn * ( commu % lsend_size(ineig+1) - 1 ) + 1 - inise
       finse  = inise + nsend - 1
       inire  = ndofn * ( commu % lrecv_size(ineig)   - 1 ) + 1
       nrecv  = ndofn * ( commu % lrecv_size(ineig+1) - 1 ) + 1 - inire
       finre  = inire + nrecv - 1
       nsend4 = int(nsend,4)
       nrecv4 = int(nrecv,4)
       dom_i4 = int(dom_i,4)

#ifndef MPI_OFF
       if( asynch .and. commu % nneig /= 0_ip ) then

          if( nsend > 0 ) then
             kk = kk + 1
             call MPI_Isend(                            &
                  xx_send(inise:finse), nsend4,         &
                  PAR_REAL,  dom_i4, 0_4,   &
                  PAR_COMM_TO_USE, ireqq4(kk), istat4   )
          end if
          if( nrecv > 0 ) then
             kk = kk + 1
             call MPI_Irecv(                            &
                  xx_recv(inire:finre), nrecv4,         &
                  PAR_REAL,  dom_i4, 0_4,   &
                  PAR_COMM_TO_USE, ireqq4(kk), istat4   )
          end if
       else
          if( nrecv /= 0 .and. nsend == 0 ) then
             call MPI_Recv(                          &
                  xx_recv(inire:finre), nrecv4,      &
                  PAR_REAL, dom_i4, 0_4, &
                  PAR_COMM_TO_USE, status, istat4    )
          else if( nrecv == 0 .and. nsend /= 0 ) then
             call MPI_Send(                          &
                  xx_send(inise:finse), nsend4,      &
                  PAR_REAL, dom_i4, 0_4, &
                  PAR_COMM_TO_USE, istat4            )
          else if( nrecv /= 0 .and. nsend /= 0 ) then
             call MPI_Sendrecv(                      &
                  xx_send(inise:finse), nsend4,      &
                  PAR_REAL, dom_i4, 0_4, &
                  xx_recv(inire:finre), nrecv4,      &
                  PAR_REAL, dom_i4, 0_4, &
                  PAR_COMM_TO_USE, status, istat4    )
          end if
       end if
#endif
    end do
    !
    ! Wait in case of asynchronous, count is the total of communications actually performe
    ! and is smaller or equal than the asynchronous communications requested (ireqq)
    !
#ifndef MPI_OFF
    if( asynch .and. commu % nneig /= 0_ip ) then
       counti = kk
       count4 = int(counti,4)
       allocate( status4(MPI_STATUS_SIZE*counti) )
       CALL MPI_WAITALL(count4,ireqq4,status4,istat4)
       deallocate( status4 )
       deallocate( ireqq4  )
    end if
#endif

  end subroutine PAR_SEND_RECEIVE_TO_ALL_RP

  !-----------------------------------------------------------------------
  !
  !> @brief   Bridge to broadcast data
  !> @details Bridge to broadcast data
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine PAR_EXCHANGE_IP_s(xx,xarra,icoun,ipass,READ_DATA,WRITE_DATA)
    
    integer(ip), intent(inout)           :: xx
    integer(ip), intent(inout), pointer  :: xarra(:)
    integer(ip), intent(inout)           :: icoun
    integer(ip), intent(in)              :: ipass
    logical(lg), intent(in),    optional :: READ_DATA,WRITE_DATA
    logical(lg)                          :: if_read,if_write

    if_read  = .false.
    if_write = .false.
    if( present(READ_DATA ) ) if_read  = READ_DATA
    if( present(WRITE_DATA) ) if_write = WRITE_DATA

    icoun = icoun + 1
    if( ipass == 2 ) then
       if( if_read ) then
          xx = xarra(icoun)
       else if( if_write ) then
          xarra(icoun) = xx          
       else
          if( IMASTER ) xarra(icoun) = xx
          if( ISLAVE  ) xx           = xarra(icoun)
       end if
    end if
  end subroutine PAR_EXCHANGE_IP_s
  subroutine PAR_EXCHANGE_IP_0(n,xx,xarra,icoun,ipass,READ_DATA,WRITE_DATA)
    integer(ip), intent(in)             :: n
    integer(ip), intent(inout)          :: xx(:)
    integer(ip), intent(inout), pointer :: xarra(:)
    integer(ip), intent(inout)          :: icoun
    integer(ip), intent(in)             :: ipass
    integer(ip)                         :: ii
    logical(lg), intent(in),    optional :: READ_DATA,WRITE_DATA
    logical(lg)                          :: if_read,if_write

    if_read  = .false.
    if_write = .false.
    if( present(READ_DATA ) ) if_read  = READ_DATA
    if( present(WRITE_DATA) ) if_write = WRITE_DATA

    do ii = 1,n
       icoun = icoun + 1
       if( ipass == 2 ) then
          if( if_read ) then
             xx(ii) = xarra(icoun)
          else if( if_write ) then
             xarra(icoun) = xx(ii)
          else
             if( IMASTER ) xarra(icoun) = xx(ii)
             if( ISLAVE  ) xx(ii)       = xarra(icoun)
          end if
       end if
    end do
  end subroutine PAR_EXCHANGE_IP_0
  subroutine PAR_EXCHANGE_IP_1(xx,xarra,icoun,ipass,READ_DATA,WRITE_DATA)
    
    integer(ip), intent(inout), pointer :: xx(:)
    integer(ip), intent(inout), pointer :: xarra(:)
    integer(ip), intent(inout)          :: icoun
    integer(ip), intent(in)             :: ipass
    integer(ip)                         :: ii
    logical(lg), intent(in),    optional :: READ_DATA,WRITE_DATA
    logical(lg)                          :: if_read,if_write
    if_read  = .false.
    if_write = .false.
    if( present(READ_DATA ) ) if_read  = READ_DATA
    if( present(WRITE_DATA) ) if_write = WRITE_DATA

    do ii = 1,memory_size(xx)
       icoun = icoun + 1
       if( ipass == 2 ) then
          if( if_read ) then
             xx(ii) = xarra(icoun)
          else if( if_write ) then
             xarra(icoun) = xx(ii)
          else
             if( IMASTER ) xarra(icoun) = xx(ii)
             if( ISLAVE  ) xx(ii)       = xarra(icoun)
          end if
       end if
    end do
  end subroutine PAR_EXCHANGE_IP_1
  subroutine PAR_EXCHANGE_IP_2(xx,xarra,icoun,ipass,READ_DATA,WRITE_DATA)
    
    integer(ip), intent(inout), pointer :: xx(:,:)
    integer(ip), intent(inout), pointer :: xarra(:)
    integer(ip), intent(inout)          :: icoun
    integer(ip), intent(in)             :: ipass
    integer(ip)                         :: ii,jj
    logical(lg), intent(in),    optional :: READ_DATA,WRITE_DATA
    logical(lg)                          :: if_read,if_write
    if_read  = .false.
    if_write = .false.
    if( present(READ_DATA ) ) if_read  = READ_DATA
    if( present(WRITE_DATA) ) if_write = WRITE_DATA

    do jj = 1,memory_size(xx,2_ip)
       do ii = 1,memory_size(xx,1_ip)
          icoun = icoun + 1
          if( ipass == 2 ) then
             if( if_read ) then
                xx(ii,jj) = xarra(icoun)
             else if( if_write ) then
                xarra(icoun) = xx(ii,jj)
             else
                if( IMASTER ) xarra(icoun) = xx(ii,jj)
                if( ISLAVE  ) xx(ii,jj)    = xarra(icoun)
             end if
          end if
       end do
    end do
  end subroutine PAR_EXCHANGE_IP_2
  subroutine PAR_EXCHANGE_IP_02(n1,n2,xx,xarra,icoun,ipass,READ_DATA,WRITE_DATA)
    integer(ip), intent(in)             :: n1
    integer(ip), intent(in)             :: n2
    integer(ip), intent(inout)          :: xx(:,:)
    integer(ip), intent(inout), pointer :: xarra(:)
    integer(ip), intent(inout)          :: icoun
    integer(ip), intent(in)             :: ipass
    integer(ip)                         :: i1,i2
    logical(lg), intent(in),    optional :: READ_DATA,WRITE_DATA
    logical(lg)                          :: if_read,if_write
    if_read  = .false.
    if_write = .false.
    if( present(READ_DATA ) ) if_read  = READ_DATA
    if( present(WRITE_DATA) ) if_write = WRITE_DATA

    do i2 = 1,n2
       do i1 = 1,n1
          icoun = icoun + 1
          if( ipass == 2 ) then
             if( if_read ) then
                xx(i1,i2) = xarra(icoun)
             else if( if_write ) then
                xarra(icoun) = xx(i1,i2)
             else
                if( IMASTER ) xarra(icoun) = xx(i1,i2)
                if( ISLAVE  ) xx(i1,i2)    = xarra(icoun)
             end if
          end if
       end do
    end do
  end subroutine PAR_EXCHANGE_IP_02

  subroutine PAR_EXCHANGE_RP_s(xx,xarra,icoun,ipass,READ_DATA,WRITE_DATA)
    
    real(rp),    intent(inout)          :: xx
    real(rp),    intent(inout), pointer :: xarra(:)
    integer(ip), intent(inout)          :: icoun
    integer(ip), intent(in)             :: ipass
    logical(lg), intent(in),    optional :: READ_DATA,WRITE_DATA
    logical(lg)                          :: if_read,if_write
    if_read  = .false.
    if_write = .false.
    if( present(READ_DATA ) ) if_read  = READ_DATA
    if( present(WRITE_DATA) ) if_write = WRITE_DATA

    icoun = icoun + 1
    if( ipass == 2 ) then
       if( if_read ) then
          xx = xarra(icoun)
       else if( if_write ) then
          xarra(icoun) = xx
       else
          if( IMASTER ) xarra(icoun) = xx
          if( ISLAVE  ) xx           = xarra(icoun)
       end if
    end if
  end subroutine PAR_EXCHANGE_RP_s
  subroutine PAR_EXCHANGE_RP_0(n,xx,xarra,icoun,ipass,READ_DATA,WRITE_DATA)
    integer(ip), intent(in)             :: n
    real(rp),    intent(inout)          :: xx(:)
    real(rp),    intent(inout), pointer :: xarra(:)
    integer(ip), intent(inout)          :: icoun
    integer(ip), intent(in)             :: ipass
    integer(ip)                         :: ii
    logical(lg), intent(in),    optional :: READ_DATA,WRITE_DATA
    logical(lg)                          :: if_read,if_write
    if_read  = .false.
    if_write = .false.
    if( present(READ_DATA ) ) if_read  = READ_DATA
    if( present(WRITE_DATA) ) if_write = WRITE_DATA

    do ii = 1,n
       icoun = icoun + 1
       if( ipass == 2 ) then
          if( if_read ) then
             xx(ii) = xarra(icoun)
          else if( if_write ) then
             xarra(icoun) = xx(ii)
          else
             if( IMASTER ) xarra(icoun) = xx(ii)
             if( ISLAVE  ) xx(ii)       = xarra(icoun)
          end if
       end if
    end do
  end subroutine PAR_EXCHANGE_RP_0
  subroutine PAR_EXCHANGE_RP_02(n1,n2,xx,xarra,icoun,ipass,READ_DATA,WRITE_DATA)
    integer(ip), intent(in)             :: n1
    integer(ip), intent(in)             :: n2
    real(rp),    intent(inout)          :: xx(:,:)
    real(rp),    intent(inout), pointer :: xarra(:)
    integer(ip), intent(inout)          :: icoun
    integer(ip), intent(in)             :: ipass
    integer(ip)                         :: i1,i2
    logical(lg), intent(in),    optional :: READ_DATA,WRITE_DATA
    logical(lg)                          :: if_read,if_write
    if_read  = .false.
    if_write = .false.
    if( present(READ_DATA ) ) if_read  = READ_DATA
    if( present(WRITE_DATA) ) if_write = WRITE_DATA

    do i2 = 1,n2
       do i1 = 1,n1
          icoun = icoun + 1
          if( ipass == 2 ) then
             if( if_read ) then
                xx(i1,i2) = xarra(icoun)
             else if( if_write ) then
                xarra(icoun) = xx(i1,i2)
             else
                if( IMASTER ) xarra(icoun) = xx(i1,i2)
                if( ISLAVE  ) xx(i1,i2)    = xarra(icoun)
             end if
          end if
       end do
    end do
  end subroutine PAR_EXCHANGE_RP_02
  subroutine PAR_EXCHANGE_RP_1(xx,xarra,icoun,ipass,READ_DATA,WRITE_DATA)
    
    real(rp),    intent(inout), pointer :: xx(:)
    real(rp),    intent(inout), pointer :: xarra(:)
    integer(ip), intent(inout)          :: icoun
    integer(ip), intent(in)             :: ipass
    integer(ip)                         :: ii
    logical(lg), intent(in),    optional :: READ_DATA,WRITE_DATA
    logical(lg)                          :: if_read,if_write
    if_read  = .false.
    if_write = .false.
    if( present(READ_DATA ) ) if_read  = READ_DATA
    if( present(WRITE_DATA) ) if_write = WRITE_DATA

    do ii = 1,memory_size(xx)
       icoun = icoun + 1
       if( ipass == 2 ) then
          if( if_read ) then
             xx(ii) = xarra(icoun)
          else if( if_write ) then
             xarra(icoun) = xx(ii)
          else
             if( IMASTER ) xarra(icoun) = xx(ii)
             if( ISLAVE  ) xx(ii)       = xarra(icoun)
          end if
       end if
    end do
  end subroutine PAR_EXCHANGE_RP_1
  subroutine PAR_EXCHANGE_RP_2(xx,xarra,icoun,ipass,READ_DATA,WRITE_DATA)
    real(rp),    intent(inout), pointer :: xx(:,:)
    real(rp),    intent(inout), pointer :: xarra(:)
    integer(ip), intent(inout)          :: icoun
    integer(ip), intent(in)             :: ipass
    integer(ip)                         :: ii,jj
    logical(lg), intent(in),    optional :: READ_DATA,WRITE_DATA
    logical(lg)                          :: if_read,if_write
    if_read  = .false.
    if_write = .false.
    if( present(READ_DATA ) ) if_read  = READ_DATA
    if( present(WRITE_DATA) ) if_write = WRITE_DATA

    do jj = 1,memory_size(xx,2_ip)
       do ii = 1,memory_size(xx,1_ip)
          icoun = icoun + 1
          if( ipass == 2 ) then
             if( if_read ) then
                xx(ii,jj) = xarra(icoun)
             else if( if_write ) then
                xarra(icoun) = xx(ii,jj)
             else
                if( IMASTER ) xarra(icoun) = xx(ii,jj)
                if( ISLAVE  ) xx(ii,jj)    = xarra(icoun)
             end if
          end if
       end do
    end do
  end subroutine PAR_EXCHANGE_RP_2

  subroutine PAR_EXCHANGE_CH(xx,xarra,icoun,ipass,READ_DATA,WRITE_DATA)
    character(*), intent(inout) :: xx
    character(*), intent(inout) :: xarra
    integer(ip),  intent(inout) :: icoun
    integer(ip),  intent(in)    :: ipass
    integer(ip)                 :: iposi,len_xarra,len_xx
    logical(lg), intent(in),    optional :: READ_DATA,WRITE_DATA
    logical(lg)                          :: if_read,if_write
    if_read  = .false.
    if_write = .false.
    if( present(READ_DATA ) ) if_read  = READ_DATA
    if( present(WRITE_DATA) ) if_write = WRITE_DATA

    len_xx    = len(xx)
    len_xarra = len(xarra)
    iposi     = icoun
    icoun     = icoun + len_xx
    if( ipass == 2 ) then
       if( iposi+len_xx > len_xarra ) call runend('CEXCHA: MAXIMUM PARCH LENGTH EXCEEDED')
       if( if_read ) then
          xx(1:len_xx) = xarra(iposi+1:iposi+len_xx)
       else if( if_write ) then
          xarra(iposi+1:iposi+len_xx) = xx(1:len_xx)
       else
          if( IMASTER ) xarra(iposi+1:iposi+len_xx) = xx(1:len_xx)
          if( ISLAVE  ) xx(1:len_xx) = xarra(iposi+1:iposi+len_xx)
       end if
    end if
  end subroutine PAR_EXCHANGE_CH

  subroutine PAR_EXCHANGE_CH_1(ncoun,xx,xarra,icoun,ipass,READ_DATA,WRITE_DATA)

    integer(ip),      intent(in)             :: ncoun
    character(len=*), intent(inout)          :: xx
    character(len=:), intent(inout), pointer :: xarra
    integer(ip),      intent(inout)          :: icoun
    integer(ip),      intent(in)             :: ipass
    integer(ip)                              :: iposi,len_xarra,len_xx
    logical(lg), intent(in),    optional :: READ_DATA,WRITE_DATA
    logical(lg)                          :: if_read,if_write
    if_read  = .false.
    if_write = .false.
    if( present(READ_DATA ) ) if_read  = READ_DATA
    if( present(WRITE_DATA) ) if_write = WRITE_DATA

    len_xx    = ncoun !len(xx)
    len_xarra = len(xarra)
    iposi     = icoun
    icoun     = icoun + len_xx
    if( ipass == 2 ) then
       if( iposi+len_xx > len_xarra ) call runend('CEXCHA: MAXIMUM PARCH LENGTH EXCEEDED')
       if( if_read ) then
          xx(1:len_xx) = xarra(iposi+1:iposi+len_xx)
       else if( if_write ) then
          xarra(iposi+1:iposi+len_xx) = xx(1:len_xx)
       else
          if( IMASTER ) xarra(iposi+1:iposi+len_xx) = xx(1:len_xx)
          if( ISLAVE  ) xx(1:len_xx) = xarra(iposi+1:iposi+len_xx)
       end if
    end if
    
  end subroutine PAR_EXCHANGE_CH_1

  subroutine PAR_EXCHANGE_LG_s(xx,xarra,icoun,ipass,READ_DATA,WRITE_DATA)
    logical(lg), intent(inout)          :: xx
    logical(lg), intent(inout), pointer :: xarra(:)
    integer(ip), intent(inout)          :: icoun
    integer(ip), intent(in)             :: ipass
    logical(lg), intent(in),    optional :: READ_DATA,WRITE_DATA
    logical(lg)                          :: if_read,if_write
    if_read  = .false.
    if_write = .false.
    if( present(READ_DATA ) ) if_read  = READ_DATA
    if( present(WRITE_DATA) ) if_write = WRITE_DATA

    icoun = icoun + 1
    if( ipass == 2 ) then
       if( if_read ) then
          xx           = xarra(icoun)
       else if( if_write ) then
          xarra(icoun) = xx
       else
          if( IMASTER ) xarra(icoun) = xx
          if( ISLAVE  ) xx           = xarra(icoun)
       end if
    end if
  end subroutine PAR_EXCHANGE_LG_s
  subroutine PAR_EXCHANGE_LG_0(n,xx,xarra,icoun,ipass,READ_DATA,WRITE_DATA)
    integer(ip), intent(in)             :: n
    logical(lg), intent(inout)          :: xx(:)
    logical(lg), intent(inout), pointer :: xarra(:)
    integer(ip), intent(inout)          :: icoun
    integer(ip), intent(in)             :: ipass
    integer(ip)                         :: ii
    logical(lg), intent(in),    optional :: READ_DATA,WRITE_DATA
    logical(lg)                          :: if_read,if_write
    if_read  = .false.
    if_write = .false.
    if( present(READ_DATA ) ) if_read  = READ_DATA
    if( present(WRITE_DATA) ) if_write = WRITE_DATA

    do ii = 1,n
       icoun = icoun + 1
       if( ipass == 2 ) then
          if( if_read ) then
             xx(ii)       = xarra(icoun)
          else if( if_write ) then
             xarra(icoun) = xx(ii)
          else
             if( IMASTER ) xarra(icoun) = xx(ii)
             if( ISLAVE  ) xx(ii)       = xarra(icoun)
          end if
       end if
    end do
  end subroutine PAR_EXCHANGE_LG_0
  subroutine PAR_EXCHANGE_LG_1(xx,xarra,icoun,ipass,READ_DATA,WRITE_DATA)
    logical(lg), intent(inout), pointer :: xx(:)
    logical(lg), intent(inout), pointer :: xarra(:)
    integer(ip), intent(inout)          :: icoun
    integer(ip), intent(in)             :: ipass
    integer(ip)                         :: ii
    logical(lg), intent(in),    optional :: READ_DATA,WRITE_DATA
    logical(lg)                          :: if_read,if_write
    if_read  = .false.
    if_write = .false.
    if( present(READ_DATA ) ) if_read  = READ_DATA
    if( present(WRITE_DATA) ) if_write = WRITE_DATA

    do ii = 1,memory_size(xx)
       icoun = icoun + 1
       if( ipass == 2 ) then
          if( if_read ) then
             xx(ii)       = xarra(icoun)
          else if( if_write ) then
             xarra(icoun) = xx(ii)
          else
             if( IMASTER ) xarra(icoun) = xx(ii)
             if( ISLAVE  ) xx(ii)       = xarra(icoun)
          end if
       end if
    end do
  end subroutine PAR_EXCHANGE_LG_1

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    16/05/2014
  !> @brief   Initialize non-blocking communications
  !> @details Initialize non-blocking communications
  !
  !----------------------------------------------------------------------

  subroutine PAR_INITIALIZE_NON_BLOCKING_COMM()
    integer(ip) :: ii
    
    allocate(non_blocking(10))
    do ii = 1,size(non_blocking)
       non_blocking(ii) % count4 = 0
       nullify(non_blocking(ii) % request4)
    end do
  end subroutine PAR_INITIALIZE_NON_BLOCKING_COMM

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    16/05/2014
  !> @brief   Start non-blocking communications
  !> @details Start non-blocking communications
  !>
  !>          KNONBLOCKING ............. communication number
  !>          MAX_NUMBER_NEIGHBORS*2 ... maximum number of requests
  !
  !----------------------------------------------------------------------

  subroutine PAR_START_NON_BLOCKING_COMM_8(knonblocking,max_number_neighbors)
    integer(ip), intent(in)           :: knonblocking
    integer(8),  intent(in)           :: max_number_neighbors

    if( knonblocking <= 0 .or. knonblocking > size(non_blocking) ) then
       call runend('PAR_START_NON_BLOCKING_COMM: WRONG NON-BLOCKING STRUCTURE NUMBER')
    else
       if( non_blocking(knonblocking) % count4 /= 0 ) then
          call runend('PAR_START_NON_BLOCKING_COMM: AN ASYNCHRONOUS SEND/RECEIVE IS ALREADY STARTED')
       else
          if (associated(non_blocking(knonblocking) % request4) ) then
             call runend('PAR_START_NON_BLOCKING_COMM: A NON-BLOCKING COMMUNICATOR NUMBER IS ALREADY ASSOCIATED')
          else
             if( max_number_neighbors > 0 ) &
                  allocate(non_blocking(knonblocking) % request4(max_number_neighbors*2))
             non_blocking(knonblocking) % count4 = 0
          end if
       end if
    end if
  end subroutine PAR_START_NON_BLOCKING_COMM_8

  subroutine PAR_START_NON_BLOCKING_COMM_4(knonblocking,max_number_neighbors)
    integer(ip), intent(in)           :: knonblocking
    integer(4),  intent(in)           :: max_number_neighbors

    if( knonblocking <= 0 .or. knonblocking > size(non_blocking) ) then
       write(*,*)'knonblocking,size(non_blocking)',knonblocking,size(non_blocking)
       call runend('PAR_START_NON_BLOCKING_COMM: WRONG NON-BLOCKING STRUCTURE NUMBER')
    else
       if( non_blocking(knonblocking) % count4 /= 0 ) then
          call runend('PAR_START_NON_BLOCKING_COMM: AN ASYNCHRONOUS SEND/RECEIVE IS ALREADY STARTED')
       else
          if (associated(non_blocking(knonblocking) % request4) ) then
             call runend('PAR_START_NON_BLOCKING_COMM: A NON-BLOCKING COMMUNICATOR NUMBER IS ALREADY ASSOCIATED')
          else
             allocate(non_blocking(knonblocking) % request4(max_number_neighbors*2))
             non_blocking(knonblocking) % count4 = 0
          end if
       end if
    end if
  end subroutine PAR_START_NON_BLOCKING_COMM_4

  subroutine PAR_START_NON_BLOCKING_COMM_4_OPT(knonblocking)
    integer(ip), intent(in)           :: knonblocking
    integer(ip)                       :: max_number_neighbors
!    integer(ip)                       :: PAR_CURRENT_SIZE

    !call PAR_COMM_RANK_AND_SIZE(commd % PAR_COMM_WORLD,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
    
    max_number_neighbors = commd % SIZE4 

    if( knonblocking <= 0 .or. knonblocking > size(non_blocking) ) then
       write(*,*)'knonblocking,size(non_blocking)',knonblocking,size(non_blocking)
       call runend('PAR_START_NON_BLOCKING_COMM: WRONG NON-BLOCKING STRUCTURE NUMBER')
    else
       if( non_blocking(knonblocking) % count4 /= 0 ) then
          call runend('PAR_START_NON_BLOCKING_COMM: AN ASYNCHRONOUS SEND/RECEIVE IS ALREADY STARTED')
       else
          if (associated(non_blocking(knonblocking) % request4) ) then
             call runend('PAR_START_NON_BLOCKING_COMM: A NON-BLOCKING COMMUNICATOR NUMBER IS ALREADY ASSOCIATED')
          else
             allocate(non_blocking(knonblocking) % request4(max_number_neighbors*2))
             non_blocking(knonblocking) % count4 = 0
          end if
       end if
    end if
  end subroutine PAR_START_NON_BLOCKING_COMM_4_OPT

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    16/05/2014
  !> @brief   End non-blocking communications
  !> @details End non-blocking communications waiting all messages have
  !>          arrived
  !>
  !>          KNONBLOCKING ............. communication number
  !
  !----------------------------------------------------------------------

  subroutine PAR_END_NON_BLOCKING_COMM(knonblocking)
    integer(ip), intent(in)    :: knonblocking
    MY_MPI_STATUS ,  pointer :: status4(:)
    integer(4)                 :: istat4

#ifndef MPI_OFF
    if( non_blocking(knonblocking) % count4 > 0 ) then
       allocate( status4(MPI_STATUS_SIZE*non_blocking(knonblocking) % count4) )
       CALL MPI_WAITALL(non_blocking(knonblocking) % count4,non_blocking(knonblocking) % request4,status4,istat4)
       if( istat4 /= MPI_SUCCESS ) call runend('NON BLOCKING SEND/RECEIVE COULD NOT BE COMPLETED')
       deallocate( status4 )
       deallocate( non_blocking(knonblocking) % request4 )
       non_blocking(knonblocking) % count4 = 0
    end if
    if( associated(non_blocking(knonblocking) % request4) ) then
       deallocate( non_blocking(knonblocking) % request4  )
    end if
#endif

  end subroutine PAR_END_NON_BLOCKING_COMM

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    16/05/2014
  !> @brief   Set non-blocking communication number
  !> @details Set non-blocking communication number
  !>
  !>          KNONBLOCKING ............. current communication number
  !
  !----------------------------------------------------------------------

  subroutine PAR_SET_NON_BLOCKING_COMM_NUMBER(knonblocking)
    integer(ip), intent(in) :: knonblocking

    if( knonblocking < 0 .or. knonblocking > size(non_blocking) ) then
       call runend('PAR_START_NON_BLOCKING_COMM: WRONG NON-BLOCKING STRUCTURE NUMBER')
    else
       inonblocking = knonblocking
    end if

  end subroutine PAR_SET_NON_BLOCKING_COMM_NUMBER

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_GHOST_NODE_EXCHANGE_RP
  !
  !----------------------------------------------------------------------

  subroutine PAR_GHOST_NODE_EXCHANGE_RP_00(n,xx,what,wherein,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_GHOST_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_GHOST_NODE_EXCHANGE_RP_00

  subroutine PAR_GHOST_NODE_EXCHANGE_RP_0(n,xx,what,wherein,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_GHOST_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_GHOST_NODE_EXCHANGE_RP_0

  subroutine PAR_GHOST_NODE_EXCHANGE_RP_1(xx,what,wherein,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = 1
    if( associated(xx) ) then
       if( size(xx,1) /= npoin_2 ) call runend('PAR_GHOST_NODE_EXCHANGE_RP_1: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_NODE_EXCHANGE_RP_1

  subroutine PAR_GHOST_NODE_EXCHANGE_RP_2(xx,what,wherein,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    if( associated(xx) ) then
       if( size(xx,2) /= npoin_2 ) call runend('PAR_GHOST_NODE_EXCHANGE_RP_2: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_NODE_EXCHANGE_RP_2

  subroutine PAR_GHOST_NODE_EXCHANGE_RP_3(xx,what,wherein,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)*size(xx,2)
       if( size(xx,3) /= npoin_2 ) call runend('PAR_GHOST_NODE_EXCHANGE_RP_3: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_NODE_EXCHANGE_RP_3

  subroutine PAR_GHOST_NODE_EXCHANGE_RP_2b(xx,what,commu,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= npoin_2 ) call runend('PAR_GHOST_NODE_EXCHANGE_RP_2b: WRONG SIZE')
       PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
       call PAR_GHOST_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_NODE_EXCHANGE_RP_2b

  !----------------------------------------------------------------------
  !
  ! PAR_GHOST_NODE_EXCHANGE_RP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_GHOST_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: ndofn
    real(rp),                       intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    MY_MPI_COMM   ,                   intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,jj,dom_i
    integer(ip)                                   :: nsize_send,nsize_recv
    integer(ip)                                   :: ipoin,ini_send,ini_recv,kk,idofn
    integer(4)                                    :: nsize_send4,nsize_recv4
    integer(4)                                    :: istat4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    MY_MPI_STATUS ,                    pointer    :: status4(:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( commu % ghost_send_node_dim == -1 .or. commu % ghost_recv_node_dim == -1 ) &
                  call runend('FRINGE NODE EXCHANGE NOT COMPUTED')
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_rsend(commu % ghost_send_node_dim * ndofn))
             allocate(tmp_rrecv(commu % ghost_recv_node_dim * ndofn))
             istat4 = 0_4
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,commu % ghost_send_node_dim
                ipoin = commu % ghost_send_node_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_rsend(kk) = xx(idofn,ipoin)
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             kk = 0
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   dom_i      = commd % neights(ii)
                   ini_send   = ndofn * ( commd % ghost_send_node_size(ii)   -1 ) + 1
                   nsize_send = ndofn * ( commd % ghost_send_node_size(ii+1) -1 ) + 1 - ini_send
                   ini_recv   = ndofn * ( commd % ghost_recv_node_size(ii)   -1 ) + 1
                   nsize_recv = ndofn * ( commd % ghost_recv_node_size(ii+1) -1 ) + 1 - ini_recv

                   nsize_send4 = int(nsize_send,4)
                   nsize_recv4 = int(nsize_recv,4)

                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_rsend(ini_send:ini_send+nsize_send-1), nsize_send4, &
                           PAR_REAL,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_rrecv(ini_recv:ini_recv+nsize_recv-1), nsize_recv4, &
                           PAR_REAL,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                   else
                      if( nsize_recv /= 0 .and. nsize_send == 0 ) then
                         call MPI_Recv(                          &
                              tmp_rrecv(ini_recv:), nsize_recv4, &
                              PAR_REAL, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, status, istat4    )
                      else if( nsize_recv == 0 .and. nsize_send /= 0 ) then
                         call MPI_Send(                          &
                              tmp_rsend(ini_send:), nsize_send4, &
                              PAR_REAL, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, istat4            )
                      else if( nsize_recv /= 0 .and. nsize_send /= 0 ) then
                         call MPI_Sendrecv(                         &
                              tmp_rsend(ini_send:), nsize_send4,    &
                              PAR_REAL, dom_i4, 0_4,    &
                              tmp_rrecv(ini_recv:), nsize_recv4,    &
                              PAR_REAL, dom_i4, 0_4,    &
                              PAR_COMM_TO_USE, status, istat4       )
                      end if
                   end if
                   if( istat4 /= MPI_SUCCESS )call PAR_MPI_RUNEND(istat4,'PAR_GHOST_NODE_EXCHANGE_RP')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE*2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             if( istat4 /= MPI_SUCCESS ) call runend('WRONG SEND/RECEIVE')
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,commu % ghost_recv_node_dim
                   ipoin = commu % ghost_recv_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_rrecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,commu % ghost_recv_node_dim
                   ipoin = commu % ghost_recv_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_rrecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,commu % ghost_recv_node_dim
                   ipoin = commu % ghost_recv_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_rrecv(kk))
                   end do
                end do

             else if( trim(what) == 'REPLACE' .or. trim(what) == 'SUBSTITUTE' ) then
                !
                ! Replace value on my fringe nodes according to what I have received
                !
                kk = 0
                do jj = 1,commu % ghost_recv_node_dim
                   ipoin = commu % ghost_recv_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = tmp_rrecv(kk)
                   end do
                end do

             else
                call runend('UNKNOWN ORDER:'//trim(what))

              end if

             ipass = 0
             deallocate(tmp_rrecv)
             deallocate(tmp_rsend)

          end if

       end if

    end if
#endif

  end subroutine PAR_GHOST_NODE_EXCHANGE_RP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_GHOST_NODE_EXCHANGE_IP
  !
  !----------------------------------------------------------------------

  subroutine PAR_GHOST_NODE_EXCHANGE_IP_0(n,xx,what,wherein,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: n
    integer(ip),                    intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_GHOST_NODE_EXCHANGE_IP_0

  subroutine PAR_GHOST_NODE_EXCHANGE_IP_1(xx,what,wherein,wsynch,dom_k)
    
    integer(ip),          pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

   if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = 1
       if( size(xx,1) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_NODE_EXCHANGE_IP_1

  subroutine PAR_GHOST_NODE_EXCHANGE_IP_2(xx,what,wherein,wsynch,dom_k)
    
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_NODE_EXCHANGE_IP_2

  subroutine PAR_GHOST_NODE_EXCHANGE_IP_3(xx,what,wherein,wsynch,dom_k)
    
    integer(ip),          pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)*size(xx,2)
       if( size(xx,3) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_NODE_EXCHANGE_IP_3

  subroutine PAR_GHOST_NODE_EXCHANGE_IP_2b(xx,what,commu,wsynch,dom_k)
    
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),  pointer,  intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= npoin .and. size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
       call PAR_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_NODE_EXCHANGE_IP_2b

 !----------------------------------------------------------------------
  !
  ! PAR_GHOST_NODE_EXCHANGE_RP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: ndofn
    integer(ip),                    intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    MY_MPI_COMM   ,                   intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,jj,dom_i
    integer(ip)                                   :: nsize_send,nsize_recv
    integer(ip)                                   :: ipoin,ini_send,ini_recv,kk,idofn
    integer(4)                                    :: nsize_send4,nsize_recv4
    integer(4)                                    :: istat4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    MY_MPI_STATUS ,                    pointer    :: status4(:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( commu % ghost_send_node_dim == -1 .or. commu % ghost_recv_node_dim == -1 ) &
                  call runend('FRINGE NODE EXCHANGE NOT COMPUTED')
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_isend(commu % ghost_send_node_dim * ndofn))
             allocate(tmp_irecv(commu % ghost_recv_node_dim * ndofn))
             istat4 = 0_4
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,commu % ghost_send_node_dim
                ipoin = commu % ghost_send_node_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_isend(kk) = xx(idofn,ipoin)
                end do
             end do
             istat4 = 0_4
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             kk = 0
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   dom_i      = commd % neights(ii)
                   ini_send   = ndofn * ( commd % ghost_send_node_size(ii)   -1 ) + 1
                   nsize_send = ndofn * ( commd % ghost_send_node_size(ii+1) -1 ) + 1 - ini_send
                   ini_recv   = ndofn * ( commd % ghost_recv_node_size(ii)   -1 ) + 1
                   nsize_recv = ndofn * ( commd % ghost_recv_node_size(ii+1) -1 ) + 1 - ini_recv

                   nsize_send4 = int(nsize_send,4)
                   nsize_recv4 = int(nsize_recv,4)

                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_isend(ini_send:ini_send+nsize_send-1), nsize_send4, &
                           PAR_INTEGER,  dom_i4, 0_4,                              &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_irecv(ini_recv:ini_recv+nsize_recv-1), nsize_recv4, &
                           PAR_INTEGER,  dom_i4, 0_4,                              &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                   else
                      if( nsize_recv /= 0 .and. nsize_send == 0 ) then
                         call MPI_Recv(                          &
                              tmp_irecv(ini_recv:), nsize_recv4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              PAR_COMM_TO_USE, status, istat4    )
                      else if( nsize_recv == 0 .and. nsize_send /= 0 ) then
                         call MPI_Send(                          &
                              tmp_isend(ini_send:), nsize_send4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              PAR_COMM_TO_USE, istat4            )
                      else if( nsize_recv /= 0 .and. nsize_send /= 0 ) then
                         call MPI_Sendrecv(                      &
                              tmp_isend(ini_send:), nsize_send4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              tmp_irecv(ini_recv:), nsize_recv4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              PAR_COMM_TO_USE, status, istat4    )
                      end if
                   end if
                   if( istat4 /= MPI_SUCCESS )call PAR_MPI_RUNEND(istat4,'PAR_GHOST_NODE_EXCHANGE_IP')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE*2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             if( istat4 /= MPI_SUCCESS ) call runend('WRONG SEND/RECEIVE')
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,commu % ghost_recv_node_dim
                   ipoin = commu % ghost_recv_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_irecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,commu % ghost_recv_node_dim
                   ipoin = commu % ghost_recv_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,commu % ghost_recv_node_dim
                   ipoin = commu % ghost_recv_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'REPLACE' .or. trim(what) == 'SUBSTITUTE' ) then
                !
                ! Replace value on my fringe nodes according to what I have received
                !
                kk = 0
                do jj = 1,commu % ghost_recv_node_dim
                   ipoin = commu % ghost_recv_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = tmp_irecv(kk)
                   end do
                end do

             else
                call runend('UNKNOWN ORDER:'//trim(what))
             end if

             ipass = 0
             deallocate(tmp_irecv)
             deallocate(tmp_isend)

          end if

       end if

    end if
#endif

  end subroutine PAR_GHOST_NODE_EXCHANGE_IP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_GHOST_ELEMENT_EXCHANGE_IP
  !
  !----------------------------------------------------------------------

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_IP_0(n,xx,what,wherein,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: n
    integer(ip),                    intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_IP_0

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_IP_1(xx,what,wherein,wsynch,dom_k)
    
    integer(ip),          pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = 1
       if( size(xx,1) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_IP_1

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_IP_2(xx,what,wherein,wsynch,dom_k)
    
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_IP_2

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_IP_3(xx,what,wherein,wsynch,dom_k)
    
    integer(ip),          pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)*size(xx,2)
       if( size(xx,3) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_IP_3

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_IP_2b(xx,what,commu,wsynch,dom_k)
    
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),  pointer,  intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
       call PAR_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_IP_2b

 !----------------------------------------------------------------------
  !
  ! PAR_GHOST_ELEMENT_EXCHANGE_RP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: ndofn
    integer(ip),                    intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    MY_MPI_COMM   ,                   intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,jj,dom_i
    integer(ip)                                   :: nsize_send,nsize_recv
    integer(ip)                                   :: ipoin,ini_send,ini_recv,kk,idofn
    integer(4)                                    :: nsize_send4,nsize_recv4
    integer(4)                                    :: istat4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    MY_MPI_STATUS ,                    pointer    :: status4(:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_ELEM_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( commu % ghost_send_elem_dim == -1 .or. commu % ghost_recv_elem_dim == -1 ) &
                  call runend('FRINGE NODE EXCHANGE NOT COMPUTED')
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_isend(commu % ghost_send_elem_dim * ndofn))
             allocate(tmp_irecv(commu % ghost_recv_elem_dim * ndofn))
             istat4 = 0_4
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,commu % ghost_send_elem_dim
                ipoin = commu % ghost_send_elem_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_isend(kk) = xx(idofn,ipoin)
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             kk = 0
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   dom_i      = commd % neights(ii)
                   ini_send   = ndofn * ( commd % ghost_send_elem_size(ii)   -1 ) + 1
                   nsize_send = ndofn * ( commd % ghost_send_elem_size(ii+1) -1 ) + 1 - ini_send
                   ini_recv   = ndofn * ( commd % ghost_recv_elem_size(ii)   -1 ) + 1
                   nsize_recv = ndofn * ( commd % ghost_recv_elem_size(ii+1) -1 ) + 1 - ini_recv

                   nsize_send4 = int(nsize_send,4)
                   nsize_recv4 = int(nsize_recv,4)

                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_isend(ini_send:ini_send+nsize_send-1), nsize_send4, &
                           PAR_INTEGER,  dom_i4, 0_4,                              &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_irecv(ini_recv:ini_recv+nsize_recv-1), nsize_recv4, &
                           PAR_INTEGER,  dom_i4, 0_4,                              &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                   else
                      if( nsize_recv /= 0 .and. nsize_send == 0 ) then
                         call MPI_Recv(                          &
                              tmp_irecv(ini_recv:), nsize_recv4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              PAR_COMM_TO_USE, status, istat4    )
                      else if( nsize_recv == 0 .and. nsize_send /= 0 ) then
                         call MPI_Send(                          &
                              tmp_isend(ini_send:), nsize_send4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              PAR_COMM_TO_USE, istat4            )
                      else if( nsize_recv /= 0 .and. nsize_send /= 0 ) then
                         call MPI_Sendrecv(                      &
                              tmp_isend(ini_send:), nsize_send4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              tmp_irecv(ini_recv:), nsize_recv4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              PAR_COMM_TO_USE, status, istat4    )
                      end if
                   end if
                   if( istat4 /= MPI_SUCCESS )call PAR_MPI_RUNEND(istat4,'PAR_GHOST_ELEMENT_EXCHANGE_RP')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE*2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             if( istat4 /= MPI_SUCCESS ) call runend('WRONG SEND/RECEIVE')
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,commu % ghost_recv_elem_dim
                   ipoin = commu % ghost_recv_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_irecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,commu % ghost_recv_elem_dim
                   ipoin = commu % ghost_recv_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,commu % ghost_recv_elem_dim
                   ipoin = commu % ghost_recv_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'REPLACE' .or. trim(what) == 'SUBSTITUTE' ) then
                !
                ! Replace value on my fringe nodes according to what I have received
                !
                kk = 0
                do jj = 1,commu % ghost_recv_elem_dim
                   ipoin = commu % ghost_recv_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = tmp_irecv(kk)
                   end do
                end do

             else
                call runend('UNKNOWN ORDER:'//trim(what))
             end if

             ipass = 0
             deallocate(tmp_irecv)
             deallocate(tmp_isend)

          end if

       end if

    end if
#endif

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_IP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_GHOST_ELEMENT_EXCHANGE_RP
  !
  !----------------------------------------------------------------------

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP_00(n,xx,what,wherein,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP_00

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP_0(n,xx,what,wherein,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP_0

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP_1(xx,what,wherein,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = 1
       if( size(xx,1) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP_1

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP_2(xx,what,wherein,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    if( size(xx,2) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP_2

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP_3(xx,what,wherein,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)*size(xx,2)
       if( size(xx,3) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP_3

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP_2b(xx,what,commu,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
       call PAR_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP_2b

  !----------------------------------------------------------------------
  !
  ! PAR_GHOST_ELEMENT_EXCHANGE_RP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: ndofn
    real(rp),                       intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    MY_MPI_COMM   ,                   intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,jj,dom_i
    integer(ip)                                   :: nsize_send,nsize_recv
    integer(ip)                                   :: ipoin,ini_send,ini_recv,kk,idofn
    integer(4)                                    :: nsize_send4,nsize_recv4
    integer(4)                                    :: istat4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    MY_MPI_STATUS ,                    pointer    :: status4(:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_ELEMENT_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( commu % ghost_send_elem_dim == -1 .or. commu % ghost_recv_elem_dim == -1 ) &
                  call runend('FRINGE NODE EXCHANGE NOT COMPUTED')
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_rsend(commu % ghost_send_elem_dim * ndofn))
             allocate(tmp_rrecv(commu % ghost_recv_elem_dim * ndofn))
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,commu % ghost_send_elem_dim
                ipoin = commu % ghost_send_elem_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_rsend(kk) = xx(idofn,ipoin)
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             kk = 0
             istat4 = 0_4
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   dom_i      = commd % neights(ii)
                   ini_send   = ndofn * ( commd % ghost_send_elem_size(ii)   -1 ) + 1
                   nsize_send = ndofn * ( commd % ghost_send_elem_size(ii+1) -1 ) + 1 - ini_send
                   ini_recv   = ndofn * ( commd % ghost_recv_elem_size(ii)   -1 ) + 1
                   nsize_recv = ndofn * ( commd % ghost_recv_elem_size(ii+1) -1 ) + 1 - ini_recv

                   nsize_send4 = int(nsize_send,4)
                   nsize_recv4 = int(nsize_recv,4)

                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_rsend(ini_send:ini_send+nsize_send-1), nsize_send4, &
                           PAR_REAL,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_rrecv(ini_recv:ini_recv+nsize_recv-1), nsize_recv4, &
                           PAR_REAL,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                   else
                      if( nsize_recv /= 0 .and. nsize_send == 0 ) then
                         call MPI_Recv(                          &
                              tmp_rrecv(ini_recv:), nsize_recv4, &
                              PAR_REAL, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, status, istat4    )
                      else if( nsize_recv == 0 .and. nsize_send /= 0 ) then
                         call MPI_Send(                          &
                              tmp_rsend(ini_send:), nsize_send4, &
                              PAR_REAL, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, istat4            )
                      else if( nsize_recv /= 0 .and. nsize_send /= 0 ) then
                         call MPI_Sendrecv(                         &
                              tmp_rsend(ini_send:), nsize_send4,    &
                              PAR_REAL, dom_i4, 0_4,    &
                              tmp_rrecv(ini_recv:), nsize_recv4,    &
                              PAR_REAL, dom_i4, 0_4,    &
                              PAR_COMM_TO_USE, status, istat4       )
                      end if
                   end if
                   if( istat4 /= MPI_SUCCESS )call PAR_MPI_RUNEND(istat4,'PAR_GHOST_ELEMENT_EXCHANGE_RP')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE*2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             if( istat4 /= MPI_SUCCESS ) call runend('WRONG SEND/RECEIVE')
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,commu % ghost_recv_elem_dim
                   ipoin = commu % ghost_recv_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_rrecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,commu % ghost_recv_elem_dim
                   ipoin = commu % ghost_recv_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_rrecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,commu % ghost_recv_elem_dim
                   ipoin = commu % ghost_recv_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_rrecv(kk))
                   end do
                end do

             else if( trim(what) == 'REPLACE' .or. trim(what) == 'SUBSTITUTE' ) then
                !
                ! Replace value on my fringe nodes according to what I have received
                !
                kk = 0
                do jj = 1,commu % ghost_recv_elem_dim
                   ipoin = commu % ghost_recv_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = tmp_rrecv(kk)
                   end do
                end do

             else
                call runend('UNKNOWN ORDER:'//trim(what))
             end if

             ipass = 0
             deallocate(tmp_rrecv)
             deallocate(tmp_rsend)

          end if

       end if

    end if
#endif

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_GHOST_BOUNDARY_EXCHANGE_IP
  !
  !----------------------------------------------------------------------

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_IP_0(n,xx,what,wherein,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: n
    integer(ip),                    intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_IP_0

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_IP_1(xx,what,wherein,wsynch,dom_k)
    
    integer(ip),          pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = 1
    if( associated(xx) ) then
       !
       ! Boundary arrays are generally allocated to a minimum of 1
       ! for robustness reasons, therefore we have not SIZE(XX) = NBOUN_2
       ! when NBOUN_2 = 0
       !
       if( size(xx,1) /= nboun_2 .and. nboun_2 /= 0 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_IP_1

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_IP_2(xx,what,wherein,wsynch,dom_k)
    
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    if( associated(xx) ) then
       !
       ! Boundary arrays are generally allocated to a minimum of 1
       ! for robustness reasons, therefore we have not SIZE(XX) = NBOUN_2
       ! when NBOUN_2 = 0
       !
       if( size(xx,2) /= nboun_2 .and. nboun_2 /= 0 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_IP_2

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_IP_3(xx,what,wherein,wsynch,dom_k)
    
    integer(ip),          pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = size(xx,1)*size(xx,2)
    if( associated(xx) ) then
       !
       ! Boundary arrays are generally allocated to a minimum of 1
       ! for robustness reasons, therefore we have not SIZE(XX) = NBOUN_2
       ! when NBOUN_2 = 0
       !
       if( size(xx,3) /= nboun_2 .and. nboun_2 /= 0 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_IP_3

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_IP_2b(xx,what,commu,wsynch,dom_k)
    
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),  pointer,  intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    if( associated(xx) ) then
       !
       ! Boundary arrays are generally allocated to a minimum of 1
       ! for robustness reasons, therefore we have not SIZE(XX) = NBOUN_2
       ! when NBOUN_2 = 0
       !
       if( size(xx,2) /= nboun_2 .and. nboun_2 /= 0 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
       call PAR_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_IP_2b

 !----------------------------------------------------------------------
  !
  ! PAR_GHOST_BOUNDARY_EXCHANGE_RP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: ndofn
    integer(ip),                    intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    MY_MPI_COMM   ,                   intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,jj,dom_i,dummi
    integer(ip)                                   :: nsize_send,nsize_recv
    integer(ip)                                   :: ipoin,ini_send,ini_recv,kk,idofn
    integer(4)                                    :: nsize_send4,nsize_recv4
    integer(4)                                    :: istat4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    MY_MPI_STATUS ,                    pointer    :: status4(:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_BOUN_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( commu % ghost_send_boun_dim == -1 .or. commu % ghost_recv_boun_dim == -1 ) &
                  call runend('FRINGE NODE EXCHANGE NOT COMPUTED')
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_isend(commu % ghost_send_boun_dim * ndofn))
             allocate(tmp_irecv(commu % ghost_recv_boun_dim * ndofn))
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,commu % ghost_send_boun_dim
                ipoin = commu % ghost_send_boun_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_isend(kk) = xx(idofn,ipoin)
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             kk = 0
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   dom_i      = commd % neights(ii)
                   ini_send   = ndofn * ( commd % ghost_send_boun_size(ii)   -1 ) + 1
                   nsize_send = ndofn * ( commd % ghost_send_boun_size(ii+1) -1 ) + 1 - ini_send
                   ini_recv   = ndofn * ( commd % ghost_recv_boun_size(ii)   -1 ) + 1
                   nsize_recv = ndofn * ( commd % ghost_recv_boun_size(ii+1) -1 ) + 1 - ini_recv

                   nsize_send4 = int(nsize_send,4)
                   nsize_recv4 = int(nsize_recv,4)

                   if( asynch ) then
                      kk = kk + 1

                      if( nsize_send > 0 ) then
                         call MPI_Isend(&
                              dummi, nsize_send4,                                     &
                              PAR_INTEGER,  dom_i4, 0_4,                              &
                              PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      else
                         call MPI_Isend(&
                              tmp_isend(ini_send:ini_send+nsize_send-1), nsize_send4, &
                              PAR_INTEGER,  dom_i4, 0_4,                              &
                              PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      end if

                      kk = kk + 1

                      if( nsize_recv > 0 ) then
                         call MPI_Irecv(&
                              tmp_irecv(ini_recv:ini_recv+nsize_recv-1), nsize_recv4, &
                              PAR_INTEGER,  dom_i4, 0_4,                              &
                              PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      else
                         call MPI_Irecv(&
                              dummi, nsize_recv4,                                     &
                              PAR_INTEGER,  dom_i4, 0_4,                              &
                              PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      end if

                   else
                      if( nsize_recv > 0 .and. nsize_send <= 0 ) then
                         call MPI_Recv(                          &
                              tmp_irecv(ini_recv:), nsize_recv4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              PAR_COMM_TO_USE, status, istat4    )
                      else if( nsize_recv <= 0 .and. nsize_send > 0 ) then
                         call MPI_Send(                          &
                              tmp_isend(ini_send:), nsize_send4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              PAR_COMM_TO_USE, istat4            )
                      else if( nsize_recv > 0 .and. nsize_send > 0 ) then
                         call MPI_Sendrecv(                      &
                              tmp_isend(ini_send:), nsize_send4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              tmp_irecv(ini_recv:), nsize_recv4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              PAR_COMM_TO_USE, status, istat4    )
                      end if
                   end if
                   if( istat4 /= MPI_SUCCESS )call PAR_MPI_RUNEND(istat4,'PAR_GHOST_BOUNDARY_EXCHANGE_IP')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE*2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             if( istat4 /= MPI_SUCCESS ) call runend('WRONG SEND/RECEIVE')
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,commu % ghost_recv_boun_dim
                   ipoin = commu % ghost_recv_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_irecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,commu % ghost_recv_boun_dim
                   ipoin = commu % ghost_recv_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,commu % ghost_recv_boun_dim
                   ipoin = commu % ghost_recv_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'REPLACE' .or. trim(what) == 'SUBSTITUTE' ) then
                !
                ! Replace value on my fringe nodes according to what I have received
                !
                kk = 0
                do jj = 1,commu % ghost_recv_boun_dim
                   ipoin = commu % ghost_recv_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = tmp_irecv(kk)
                   end do
                end do

             else
                call runend('UNKNOWN ORDER:'//trim(what))
             end if

             ipass = 0
             deallocate(tmp_irecv)
             deallocate(tmp_isend)

          end if

       end if

    end if
#endif

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_IP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_GHOST_BOUNDARY_EXCHANGE_RP
  !
  !----------------------------------------------------------------------

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP_00(n,xx,what,wherein,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP_00

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP_0(n,xx,what,wherein,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP_0

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP_1(xx,what,wherein,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = 1
    if( associated(xx) ) then
       if( size(xx,1) /= nboun_2 .and. nboun_2 /= 0 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP_1

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP_2(xx,what,wherein,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    if( associated(xx) ) then
       if( size(xx,2) /= nboun_2 .and. nboun_2 /= 0 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP_2

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP_3(xx,what,wherein,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = size(xx,1)*size(xx,2)
    if( associated(xx) ) then
       if( size(xx,3) /= nboun_2 .and. nboun_2 /= 0 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP_3

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP_2b(xx,what,commu,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    if( associated(xx) ) then
       if( size(xx,2) /= nboun_2 .and. nboun_2 /= 0 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
       call PAR_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP_2b

  !----------------------------------------------------------------------
  !
  ! PAR_GHOST_BOUNDARY_EXCHANGE_RP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: ndofn
    real(rp),                       intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    MY_MPI_COMM   ,                   intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,jj,dom_i
    integer(ip)                                   :: nsize_send,nsize_recv
    integer(ip)                                   :: ipoin,ini_send,ini_recv,kk,idofn
    integer(4)                                    :: nsize_send4,nsize_recv4
    integer(4)                                    :: istat4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    MY_MPI_STATUS ,                    pointer    :: status4(:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_BOUNDARY_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( commu % ghost_send_boun_dim == -1 .or. commu % ghost_recv_boun_dim == -1 ) &
                  call runend('FRINGE NODE EXCHANGE NOT COMPUTED')
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_rsend(commu % ghost_send_boun_dim * ndofn))
             allocate(tmp_rrecv(commu % ghost_recv_boun_dim * ndofn))
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,commu % ghost_send_boun_dim
                ipoin = commu % ghost_send_boun_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_rsend(kk) = xx(idofn,ipoin)
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             istat4 = 0_4
             kk = 0
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   dom_i      = commd % neights(ii)
                   ini_send   = ndofn * ( commd % ghost_send_boun_size(ii)   -1 ) + 1
                   nsize_send = ndofn * ( commd % ghost_send_boun_size(ii+1) -1 ) + 1 - ini_send
                   ini_recv   = ndofn * ( commd % ghost_recv_boun_size(ii)   -1 ) + 1
                   nsize_recv = ndofn * ( commd % ghost_recv_boun_size(ii+1) -1 ) + 1 - ini_recv

                   nsize_send4 = int(nsize_send,4)
                   nsize_recv4 = int(nsize_recv,4)

                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_rsend(ini_send:ini_send+nsize_send-1), nsize_send4, &
                           PAR_REAL,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_rrecv(ini_recv:ini_recv+nsize_recv-1), nsize_recv4, &
                           PAR_REAL,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                   else
                      if( nsize_recv /= 0 .and. nsize_send == 0 ) then
                         call MPI_Recv(                          &
                              tmp_rrecv(ini_recv:), nsize_recv4, &
                              PAR_REAL, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, status, istat4    )
                      else if( nsize_recv == 0 .and. nsize_send /= 0 ) then
                         call MPI_Send(                          &
                              tmp_rsend(ini_send:), nsize_send4, &
                              PAR_REAL, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, istat4            )
                      else if( nsize_recv /= 0 .and. nsize_send /= 0 ) then
                         call MPI_Sendrecv(                         &
                              tmp_rsend(ini_send:), nsize_send4,    &
                              PAR_REAL, dom_i4, 0_4,    &
                              tmp_rrecv(ini_recv:), nsize_recv4,    &
                              PAR_REAL, dom_i4, 0_4,    &
                              PAR_COMM_TO_USE, status, istat4       )
                      end if
                   end if
                   if( istat4 /= MPI_SUCCESS )call runend('PAR_GHOST_BOUNDARY_EXCHANGE_RP')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE*2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             if( istat4 /= MPI_SUCCESS ) call runend('WRONG SEND/RECEIVE')
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,commu % ghost_recv_boun_dim
                   ipoin = commu % ghost_recv_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_rrecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,commu % ghost_recv_boun_dim
                   ipoin = commu % ghost_recv_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_rrecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,commu % ghost_recv_boun_dim
                   ipoin = commu % ghost_recv_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_rrecv(kk))
                   end do
                end do

             else if( trim(what) == 'REPLACE' .or. trim(what) == 'SUBSTITUTE' ) then
                !
                ! Replace value on my fringe nodes according to what I have received
                !
                kk = 0
                do jj = 1,commu % ghost_recv_boun_dim
                   ipoin = commu % ghost_recv_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = tmp_rrecv(kk)
                   end do
                end do

             else
                call runend('UNKNOWN ORDER:'//trim(what))
             end if

             ipass = 0
             deallocate(tmp_rrecv)
             deallocate(tmp_rsend)

          end if

       end if

    end if
#endif

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP

  !----------------------------------------------------------------------
  !
  ! PAR_POINT_TO_POINT_ARRAY_OPERATION_IP: FOR INTEGERS
  ! INPUT/OUTPUT:  XX_SEND_RECV(NDOFN) for all slaves
  ! Operate on arrays between neighbors in the communicator.
  !
  !----------------------------------------------------------------------

  subroutine PAR_POINT_TO_POINT_ARRAY_OPERATION_IP_1(xx_send_recv,wherein,what)
    
    integer(ip),          pointer,  intent(inout) :: xx_send_recv(:)
    character(*),                   intent(in)    :: wherein
    character(*),                   intent(in)    :: what
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( ISEQUEN ) return

    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)

    if( .not. associated(xx_send_recv) ) then
       return
    else
       ndofn = size(xx_send_recv)
       if( ndofn > 0 ) call PAR_POINT_TO_POINT_ARRAY_OPERATION_IP(ndofn,xx_send_recv,PAR_COMM_TO_USE,commu,what)
    end if

  end subroutine PAR_POINT_TO_POINT_ARRAY_OPERATION_IP_1

  subroutine PAR_POINT_TO_POINT_ARRAY_OPERATION_IP(ndofn,xx_send_recv,PAR_COMM_TO_USE,commu,what)
    
    integer(ip),                    intent(in)    :: ndofn
    integer(ip),                    intent(inout) :: xx_send_recv(*)
    MY_MPI_COMM,                    intent(in)    :: PAR_COMM_TO_USE
    type(comm_data_par),            intent(in)    :: commu
    character(*),                   intent(in)    :: what
    integer(ip)                                   :: dom_i,ineig,ii
    integer(ip),          pointer                 :: xx_recv(:)
    integer(4)                                    :: istat4,dom_i4,my_rank4
    integer(4)                                    :: ndofn4
    !
    ! Define communicator
    !
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank4)
    !
    ! Allocate
    !
    nullify(xx_recv)
    allocate(xx_recv(ndofn))
    ndofn4 = int(ndofn,4)
    !
    ! Send/receive and operate
    !
    do ineig = 1,commu % nneig

       dom_i  = commu % neights(ineig)
       dom_i4 = int(dom_i,4)

#ifndef MPI_OFF
       call MPI_Sendrecv(                      &
            xx_send_recv(1:ndofn), ndofn4,     &
            PAR_INTEGER, dom_i4, 0_4,          &
            xx_recv(1:ndofn), ndofn4,          &
            PAR_INTEGER, dom_i4, 0_4,          &
            PAR_COMM_TO_USE, status, istat4    )
#endif
       !
       ! Operate on array
       !
       if(      trim(what) == 'MAX' ) then
          !
          ! MAX
          !
          do ii = 1,ndofn
             xx_send_recv(ii) = max(xx_send_recv(ii),xx_recv(ii))
          end do
       else if( trim(what) == 'MIN' ) then
          !
          ! MIN
          !
          do ii = 1,ndofn
             xx_send_recv(ii) = min(xx_send_recv(ii),xx_recv(ii))
          end do
       else if( trim(what) == 'SUM' ) then
          !
          ! SUM
          !
          do ii = 1,ndofn
             xx_send_recv(ii) = xx_send_recv(ii) + xx_recv(ii)
          end do
       else if( trim(what) == 'MIN RANK OR NEGATIVE' ) then
          !
          ! At the end, only one partition will have the positive sign
          ! if the value is positive
          !
          if( my_rank4 > dom_i ) then
             do ii = 1,ndofn
                if( xx_send_recv(ii) > 0 .and. xx_recv(ii) > 0 ) then
                   xx_send_recv(ii) = -abs(xx_send_recv(ii))
                end if
             end do
          end if
       end if

    end do
    !
    ! Deallocate
    !
    deallocate(xx_recv)

  end subroutine PAR_POINT_TO_POINT_ARRAY_OPERATION_IP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP
  !
  !----------------------------------------------------------------------

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_00(n,xx,what,wherein,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_00

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_0(n,xx,what,wherein,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_0

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_1(xx,what,wherein,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = 1
       if( size(xx,1) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_1

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_2(xx,what,wherein,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    if( size(xx,2) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_2

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_3(xx,what,wherein,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)*size(xx,2)
       if( size(xx,3) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_3

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_2b(xx,what,commu,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
       call PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_2b

  !----------------------------------------------------------------------
  !
  ! PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: ndofn
    real(rp),                       intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    MY_MPI_COMM   ,                   intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,jj,dom_i
    integer(ip)                                   :: nsize_send,nsize_recv
    integer(ip)                                   :: ipoin,ini_send,ini_recv,kk,idofn
    integer(4)                                    :: nsize_send4,nsize_recv4
    integer(4)                                    :: istat4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    MY_MPI_STATUS ,                    pointer    :: status4(:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_ELEMENT_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( commu % ghost_send_elem_dim == -1 .or. commu % ghost_recv_elem_dim == -1 ) &
                  call runend('FRINGE NODE EXCHANGE NOT COMPUTED')
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_rsend(commu % ghost_recv_elem_dim * ndofn))
             allocate(tmp_rrecv(commu % ghost_send_elem_dim * ndofn))
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,commu % ghost_recv_elem_dim
                ipoin = commu % ghost_recv_elem_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_rsend(kk) = xx(idofn,ipoin)
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             kk = 0
             istat4 = 0_4
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   dom_i      = commd % neights(ii)
                   ini_send   = ndofn * ( commd % ghost_recv_elem_size(ii)   -1 ) + 1
                   nsize_send = ndofn * ( commd % ghost_recv_elem_size(ii+1) -1 ) + 1 - ini_send
                   ini_recv   = ndofn * ( commd % ghost_send_elem_size(ii)   -1 ) + 1
                   nsize_recv = ndofn * ( commd % ghost_send_elem_size(ii+1) -1 ) + 1 - ini_recv

                   nsize_send4 = int(nsize_send,4)
                   nsize_recv4 = int(nsize_recv,4)

                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_rsend(ini_send:ini_send+nsize_send-1), nsize_send4, &
                           PAR_REAL,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_rrecv(ini_recv:ini_recv+nsize_recv-1), nsize_recv4, &
                           PAR_REAL,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                   else
                      if( nsize_recv /= 0 .and. nsize_send == 0 ) then
                         call MPI_Recv(                          &
                              tmp_rrecv(ini_recv:), nsize_recv4, &
                              PAR_REAL, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, status, istat4    )
                      else if( nsize_recv == 0 .and. nsize_send /= 0 ) then
                         call MPI_Send(                          &
                              tmp_rsend(ini_send:), nsize_send4, &
                              PAR_REAL, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, istat4            )
                      else if( nsize_recv /= 0 .and. nsize_send /= 0 ) then
                         call MPI_Sendrecv(                         &
                              tmp_rsend(ini_send:), nsize_send4,    &
                              PAR_REAL, dom_i4, 0_4,    &
                              tmp_rrecv(ini_recv:), nsize_recv4,    &
                              PAR_REAL, dom_i4, 0_4,    &
                              PAR_COMM_TO_USE, status, istat4       )
                      end if
                   end if
                   if( istat4 /= MPI_SUCCESS )call PAR_MPI_RUNEND(istat4,'PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE*2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             if( istat4 /= MPI_SUCCESS ) call runend('WRONG SEND/RECEIVE')
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,commu % ghost_send_elem_dim
                   ipoin = commu % ghost_send_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_rrecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,commu % ghost_send_elem_dim
                   ipoin = commu % ghost_send_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_rrecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,commu % ghost_send_elem_dim
                   ipoin = commu % ghost_send_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_rrecv(kk))
                   end do
                end do

             else if( trim(what) == 'REPLACE' .or. trim(what) == 'SUBSTITUTE' ) then
                !
                ! Replace value on my fringe nodes according to what I have received
                !
                kk = 0
                do jj = 1,commu % ghost_send_elem_dim
                   ipoin = commu % ghost_send_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = tmp_rrecv(kk)
                   end do
                end do

             else
                call runend('UNKNOWN ORDER:'//trim(what))
             end if

             ipass = 0
             deallocate(tmp_rrecv)
             deallocate(tmp_rsend)

          end if

       end if

    end if
#endif

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_INTERFACE_EDGE_EXCHANGE_IP
  !
  !----------------------------------------------------------------------

  subroutine PAR_INTERFACE_EDGE_EXCHANGE_IP_0(n,xx,what,wherein,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: n
    integer(ip),                    intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_INTERFACE_NODE_EXCHANGE_IP(ON_EDGES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_EDGE_EXCHANGE_IP_0

  subroutine PAR_INTERFACE_EDGE_EXCHANGE_IP_1(xx,what,wherein,wsynch,dom_k,COMM)
    
    integer(ip),          pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    type(comm_data_par),  optional, intent(in)    :: COMM
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

   if( INOTSLAVE ) return
   ndofn = 1
   if( present(COMM) ) then
      PAR_COMM_TO_USE = COMM % PAR_COMM_WORLD
      call PAR_INTERFACE_NODE_EXCHANGE_IP(ON_EDGES,ndofn,xx,what,COMM,PAR_COMM_TO_USE,wsynch,dom_k)
   else
      call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
      call PAR_INTERFACE_NODE_EXCHANGE_IP(ON_EDGES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
   end if
   
  end subroutine PAR_INTERFACE_EDGE_EXCHANGE_IP_1

  subroutine PAR_INTERFACE_EDGE_EXCHANGE_IP_2(xx,what,wherein,wsynch,dom_k)
    
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),         optional, intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
    end if
    call PAR_INTERFACE_NODE_EXCHANGE_IP(ON_EDGES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_EDGE_EXCHANGE_IP_2

  subroutine PAR_INTERFACE_EDGE_EXCHANGE_IP_3(xx,what,wherein,wsynch,dom_k)
    
    integer(ip),          pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    if( size(xx,3) <= 2 ) then
       ndofn = size(xx,1)
    else
       ndofn = size(xx,1)*size(xx,2)
    end if
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_INTERFACE_NODE_EXCHANGE_IP(ON_EDGES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_EDGE_EXCHANGE_IP_3

  subroutine PAR_INTERFACE_EDGE_EXCHANGE_IP_2b(xx,what,commu,wsynch,dom_k)
    
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),  pointer,  intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
    call PAR_INTERFACE_NODE_EXCHANGE_IP(ON_EDGES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_EDGE_EXCHANGE_IP_2b

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_INTERFACE_EDGE_EXCHANGE_RP
  !
  !----------------------------------------------------------------------

  subroutine PAR_INTERFACE_EDGE_EXCHANGE_RP_00(n,xx,what,wherein,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_EDGES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_EDGE_EXCHANGE_RP_00

  subroutine PAR_INTERFACE_EDGE_EXCHANGE_RP_0(n,xx,what,wherein,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_EDGES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_EDGE_EXCHANGE_RP_0

  subroutine PAR_INTERFACE_EDGE_EXCHANGE_RP_1(xx,what,wherein,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = 1
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_EDGES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_EDGE_EXCHANGE_RP_1

  subroutine PAR_INTERFACE_EDGE_EXCHANGE_RP_2(xx,what,wherein,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_EDGES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_EDGE_EXCHANGE_RP_2

  subroutine PAR_INTERFACE_EDGE_EXCHANGE_RP_3(xx,what,wherein,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    if( what(len(what):len(what)) /= 'I' ) then
       if( size(xx,3) <= 2 ) then
          ndofn = size(xx,1)
       else
          ndofn = size(xx,1)*size(xx,2)
       end if
    end if
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_EDGES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_EDGE_EXCHANGE_RP_3

  subroutine PAR_INTERFACE_EDGE_EXCHANGE_RP_2b(xx,what,commu,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
    call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_EDGES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_EDGE_EXCHANGE_RP_2b

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_INTERFACE_EDGE_EXCHANGE_LG
  !
  !----------------------------------------------------------------------

  subroutine PAR_INTERFACE_EDGE_EXCHANGE_LG_1(xx,what,wherein,wsynch,dom_k)
    
    logical(lg),          pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = 1
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_INTERFACE_NODE_EXCHANGE_LG(ON_EDGES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_EDGE_EXCHANGE_LG_1

  !----------------------------------------------------------------------
  !
  ! PAR_INTERFACE_MATRIX_EXCHANGE: EXCHANGE MATRIX ON INTERFACE NODES
  !
  !----------------------------------------------------------------------

  subroutine PAR_INTERFACE_MATRIX_EXCHANGE_WHERE(ndofn,aa,wherein)
    
    integer(ip),          intent(in)    :: ndofn
    real(rp),             intent(inout) :: aa(ndofn,ndofn,*)
    character(*),         intent(in)    :: wherein
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),  pointer       :: commu
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_INTERFACE_MATRIX_EXCHANGE_RP(ndofn,aa,commu)
  end subroutine PAR_INTERFACE_MATRIX_EXCHANGE_WHERE

  subroutine PAR_INTERFACE_MATRIX_EXCHANGE_COMM_SHAPE(ndofn,aa,commu)
    
    integer(ip),          intent(in)    :: ndofn
    real(rp),             intent(inout) :: aa(*)
    type(comm_data_par),  intent(in)    :: commu
    call PAR_INTERFACE_MATRIX_EXCHANGE_RP(ndofn,aa,commu)
  end subroutine PAR_INTERFACE_MATRIX_EXCHANGE_COMM_SHAPE

  subroutine PAR_INTERFACE_MATRIX_EXCHANGE_COMM(ndofn,aa,commu)
    
    integer(ip),          intent(in)    :: ndofn
    real(rp),             intent(inout) :: aa(ndofn,ndofn,*)
    type(comm_data_par),  intent(in)    :: commu
    call PAR_INTERFACE_MATRIX_EXCHANGE_RP(ndofn,aa,commu)
  end subroutine PAR_INTERFACE_MATRIX_EXCHANGE_COMM

  subroutine PAR_INTERFACE_MATRIX_EXCHANGE_RP(ndofn,aa,commu)
    integer(ip),          intent(in)    :: ndofn
    real(rp),             intent(inout) :: aa(ndofn,ndofn,*)
    type(comm_data_par),  intent(in)    :: commu
    integer(ip)                         :: ineig,dom_i,bsize,jj,kk,ll
    integer(ip)                         :: idofn,jdofn,ii,iz
    integer(ip)                         :: PAR_CURRENT_RANK
    integer(ip)                         :: PAR_CURRENT_SIZE
    type(r1p),            pointer       :: loc_sparr1(:),loc_rparr1(:)
    !
    ! Allocate memory
    !
    nullify(loc_sparr1)
    nullify(loc_rparr1)
    call memory_alloca(par_memor,'LOC_SPARR1','PAR_INTERFACE_MATRIX_EXCHANGE',loc_sparr1,commu % nneig)
    call memory_alloca(par_memor,'LOC_RPARR1','PAR_INTERFACE_MATRIX_EXCHANGE',loc_rparr1,commu % nneig)
    do ineig = 1,commu % nneig
       bsize = ndofn * ndofn *  commu % bound_matrix(ineig) % nzdom
       call memory_alloca(par_memor,'LOC_SPARR1','PAR_INTERFACE_MATRIX_EXCHANGE',loc_sparr1(ineig) % a,bsize)
       call memory_alloca(par_memor,'LOC_RPARR1','PAR_INTERFACE_MATRIX_EXCHANGE',loc_rparr1(ineig) % a,bsize)
    end do
    !
    ! Exchange matrix coefficients
    !
    call PAR_COMM_RANK_AND_SIZE(COMMU % PAR_COMM_WORLD,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
    call PAR_START_NON_BLOCKING_COMM(1_ip,PAR_CURRENT_SIZE)
    call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
    do ineig = 1,commu % nneig
       dom_i = commu % neights(ineig)
       bsize = ndofn * ndofn *  commu % bound_matrix(ineig) % nzdom
       jj = 0
       kk = 0
       do ii = commu % bound_size(ineig),commu % bound_size(ineig+1)-1
          jj = jj + 1
          do ll = 1,commu % bound_matrix(ineig) % nzdom_ii(jj)
             iz = commu % bound_matrix(ineig) % ja(jj) % l(ll)
             do idofn = 1,ndofn
                do jdofn = 1,ndofn
                   kk = kk + 1
                   loc_sparr1(ineig) % a(kk) = aa(idofn,jdofn,iz)
                end do
             end do
          end do
       end do
       if( bsize > 0 ) call PAR_SEND_RECEIVE(bsize,bsize,loc_sparr1(ineig)%a,loc_rparr1(ineig)%a,'IN MY CODE',dom_i,'NON BLOCKING')
       !if( bsize > 0 ) call PAR_SEND_RECEIVE(bsize,bsize,loc_sparr1(ineig)%a,loc_rparr1(ineig)%a,'IN MY CODE',dom_i)
    end do
    call PAR_END_NON_BLOCKING_COMM(1_ip)
    !
    ! Assemble the contribution of my neighbor
    !
    do ineig = 1,commu % nneig
       jj = 0
       kk = 0
       do ii = commu % bound_size(ineig),commu % bound_size(ineig+1)-1
          jj = jj + 1
          do ll = 1,commu % bound_matrix(ineig) % nzdom_ii(jj)
             iz = commu % bound_matrix(ineig) % ja(jj) % l(ll)
             do idofn = 1,ndofn
                do jdofn = 1,ndofn
                   kk = kk + 1
                   aa(idofn,jdofn,iz) = aa(idofn,jdofn,iz) + loc_rparr1(ineig) % a(kk)
                end do
             end do
          end do
       end do

    end do

    call memory_deallo( par_memor,'LOC_SPARR1','PAR_INTERFACE_MATRIX_EXCHANGE',loc_sparr1)
    call memory_deallo( par_memor,'LOC_RPARR1','PAR_INTERFACE_MATRIX_EXCHANGE',loc_rparr1)

  end subroutine PAR_INTERFACE_MATRIX_EXCHANGE_RP

  !----------------------------------------------------------------------
  !
  ! PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE: EXCHANGE MATRIX WITH HALOS ON INTERFACE NODES
  !
  !----------------------------------------------------------------------
  ! ojo aa lo tengo que partir en 2 un aasend y un aarecv  el primero es sin halos y el segundo con !!

  subroutine PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP_2(ndofn,aa_send,aa_recv,commu)
    
    integer(ip),          intent(in)             :: ndofn
    real(rp),             intent(in)             :: aa_send(ndofn,ndofn,*)
    real(rp),             intent(inout)          :: aa_recv(ndofn,ndofn,*)
    type(comm_data_par),  intent(in),   optional :: commu

    if( present(commu) ) then
       call PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP(ndofn,aa_send,aa_recv,commu)
    else
       call PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP(ndofn,aa_send,aa_recv,PAR_COMM_MY_CODE_ARRAY(1))
    end if

  end subroutine PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP_2

  subroutine PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP_1(ndofn,aa_send,aa_recv,commu)
    
    integer(ip),          intent(in)             :: ndofn
    real(rp),             intent(in)             :: aa_send(*)
    real(rp),             intent(inout)          :: aa_recv(*)
    type(comm_data_par),  intent(in),   optional :: commu

    if( present(commu) ) then
       call PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP(ndofn,aa_send,aa_recv,commu)
    else
       call PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP(ndofn,aa_send,aa_recv,PAR_COMM_MY_CODE_ARRAY(1))
    end if

  end subroutine PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP_1


  subroutine PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP(ndofn,aa_send,aa_recv,commu)

    integer(ip),          intent(in)    :: ndofn
    real(rp),             intent(in)    :: aa_send(ndofn,ndofn,*)
    real(rp),             intent(inout) :: aa_recv(ndofn,ndofn,*)
    type(comm_data_par),  intent(in)    :: commu
    integer(ip)                         :: ineig,dom_i,jj,kk,ll
    integer(ip)                         :: idofn,jdofn,ii,iz
    integer(ip)                         :: PAR_CURRENT_RANK
    integer(ip)                         :: PAR_CURRENT_SIZE
    type(r1p),            pointer       :: loc_sparr1(:),loc_rparr1(:)
    integer(ip)                         :: ssize,rsize

    !
    ! Allocate memory
    !
    nullify(loc_sparr1)
    nullify(loc_rparr1)
    call memory_alloca(par_memor,'LOC_SPARR1','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',loc_sparr1,commu % nneig)
    call memory_alloca(par_memor,'LOC_RPARR1','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',loc_rparr1,commu % nneig)
    do ineig = 1,commu % nneig
       ssize = ndofn * ndofn * commu % bound_mat_halo_send(ineig) % nzdom
       rsize = ndofn * ndofn * commu % bound_mat_halo_recv(ineig) % nzdom
       call memory_alloca(par_memor,'LOC_SPARR1','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',loc_sparr1(ineig) % a,max(1_ip,ssize))
       call memory_alloca(par_memor,'LOC_RPARR1','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',loc_rparr1(ineig) % a,max(1_ip,rsize))
    end do
    !
    ! Exchange matrix coefficients
    !
    call PAR_COMM_RANK_AND_SIZE(COMMU % PAR_COMM_WORLD,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
    call PAR_START_NON_BLOCKING_COMM(1_ip,PAR_CURRENT_SIZE)
    call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
    do ineig = 1,commu % nneig
       dom_i = commu % neights(ineig)
       ssize = ndofn * ndofn *  commu % bound_mat_halo_send(ineig) % nzdom
       rsize = ndofn * ndofn *  commu % bound_mat_halo_recv(ineig) % nzdom
       jj = 0
       kk = 0
       do ii = commu % bound_size(ineig),commu % bound_size(ineig+1)-1     ! These are all the nodes that are in the boundary with my neighbour ineig
          jj = jj + 1
          do ll = 1,commu % bound_mat_halo_send(ineig) % nzdom_ii(jj)      ! These are the number of conections to a certain node
             iz = commu % bound_mat_halo_send(ineig) % ja(jj) % l(ll)
             do idofn = 1,ndofn
                do jdofn = 1,ndofn
                   kk = kk + 1
                   loc_sparr1(ineig) % a(kk) = aa_send(idofn,jdofn,iz)
                end do
             end do
          end do
       end do
       if( ssize > 0 .or. rsize > 0 ) call PAR_SEND_RECEIVE(ssize,rsize,loc_sparr1(ineig)%a,loc_rparr1(ineig)%a,'IN MY CODE',dom_i,'NON BLOCKING')
    end do
    call PAR_END_NON_BLOCKING_COMM(1_ip)
    !
    ! Assemble the contribution of my neighbor
    !
    do ineig = 1,commu % nneig
       jj = 0
       kk = 0
       do ii = commu % bound_size(ineig),commu % bound_size(ineig+1)-1
          jj = jj + 1
          do ll = 1,commu % bound_mat_halo_recv(ineig) % nzdom_ii(jj)
             iz = commu % bound_mat_halo_recv(ineig) % ja(jj) % l(ll)
             do idofn = 1,ndofn
                do jdofn = 1,ndofn
                   kk = kk + 1
                   aa_recv(idofn,jdofn,iz) = aa_recv(idofn,jdofn,iz) + loc_rparr1(ineig) % a(kk)
                end do
             end do
          end do
       end do

    end do

    call memory_deallo( par_memor,'LOC_SPARR1','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',loc_sparr1)
    call memory_deallo( par_memor,'LOC_RPARR1','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',loc_rparr1)

  end subroutine PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_FROM_GHOST_NODE_EXCHANGE_RP
  !
  !----------------------------------------------------------------------

  subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP_00(n,xx,what,wherein,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: n
    integer(ip),                    intent(inout) :: xx(*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) then
       return
    else if( n > 0 ) then
       ndofn = n
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP_00

  subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP_0(n,xx,what,wherein,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: n
    integer(ip),                    intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_FROM_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP_0

  subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP_1(xx,what,wherein,wsynch,dom_k)
    
    integer(ip),          pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) then
       return
    else if( associated(xx) ) then
       ndofn = 1
       if( size(xx,1) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP_1

  subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP_2(xx,what,wherein,wsynch,dom_k)
    
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) then
       return
    else if( associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP_2

  subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP_3(xx,what,wherein,wsynch,dom_k)
    
    integer(ip),          pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) then
       return
    else if( associated(xx) ) then
       ndofn = size(xx,1)*size(xx,2)
       if( size(xx,3) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP_3

  subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP_2b(xx,what,commu,wsynch,dom_k)
    
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( INOTSLAVE ) then
       return
    else if( associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
       call PAR_FROM_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP_2b

  !----------------------------------------------------------------------
  !
  ! PAR_FROM_GHOST_NODE_EXCHANGE_IP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: ndofn
    integer(ip),                    intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    MY_MPI_COMM   ,                   intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,jj,dom_i
    integer(ip)                                   :: nsize_send,nsize_recv
    integer(ip)                                   :: ipoin,ini_send,ini_recv,kk,idofn
    integer(4)                                    :: nsize_send4,nsize_recv4
    integer(4)                                    :: istat4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    MY_MPI_STATUS ,                    pointer    :: status4(:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_BOUNDARY_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( commu % ghost_send_node_dim == -1 .or. commu % ghost_recv_node_dim == -1 ) &
                  call runend('FRINGE NODE EXCHANGE NOT COMPUTED')
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_isend(commu % ghost_recv_node_dim * ndofn))
             allocate(tmp_irecv(commu % ghost_send_node_dim * ndofn))
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,commu % ghost_recv_node_dim
                ipoin = commu % ghost_recv_node_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_isend(kk) = xx(idofn,ipoin)
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             kk = 0
             istat4 = 0_4
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   dom_i      = commd % neights(ii)
                   ini_send   = ndofn * ( commd % ghost_recv_node_size(ii)   -1 ) + 1
                   nsize_send = ndofn * ( commd % ghost_recv_node_size(ii+1) -1 ) + 1 - ini_send
                   ini_recv   = ndofn * ( commd % ghost_send_node_size(ii)   -1 ) + 1
                   nsize_recv = ndofn * ( commd % ghost_send_node_size(ii+1) -1 ) + 1 - ini_recv

                   nsize_send4 = int(nsize_send,4)
                   nsize_recv4 = int(nsize_recv,4)

                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_isend(ini_send:ini_send+nsize_send-1), nsize_send4, &
                           PAR_INTEGER,  dom_i4, 0_4,                              &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_irecv(ini_recv:ini_recv+nsize_recv-1), nsize_recv4, &
                           PAR_INTEGER,  dom_i4, 0_4,                              &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                   else
                      if( nsize_recv /= 0 .and. nsize_send == 0 ) then
                         call MPI_Recv(                          &
                              tmp_irecv(ini_recv:), nsize_recv4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              PAR_COMM_TO_USE, status, istat4    )
                      else if( nsize_recv == 0 .and. nsize_send /= 0 ) then
                         call MPI_Send(                          &
                              tmp_isend(ini_send:), nsize_send4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              PAR_COMM_TO_USE, istat4            )
                      else if( nsize_recv /= 0 .and. nsize_send /= 0 ) then
                         call MPI_Sendrecv(                         &
                              tmp_isend(ini_send:), nsize_send4,    &
                              PAR_INTEGER, dom_i4, 0_4,             &
                              tmp_irecv(ini_recv:), nsize_recv4,    &
                              PAR_INTEGER, dom_i4, 0_4,             &
                              PAR_COMM_TO_USE, status, istat4       )
                      end if
                   end if
                   if( istat4 /= MPI_SUCCESS )call PAR_MPI_RUNEND(istat4,'PAR_FROM_GHOST_NODE_EXCHANGE_IP')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE*2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             if( istat4 /= MPI_SUCCESS ) call runend('WRONG SEND/RECEIVE')
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,commu % ghost_send_node_dim
                   ipoin = commu % ghost_send_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_irecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,commu % ghost_send_node_dim
                   ipoin = commu % ghost_send_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,commu % ghost_send_node_dim
                   ipoin = commu % ghost_send_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'REPLACE' .or. trim(what) == 'SUBSTITUTE' ) then
                !
                ! Replace value on my fringe nodes according to what I have received
                !
                kk = 0
                do jj = 1,commu % ghost_send_node_dim
                   ipoin = commu % ghost_send_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = tmp_irecv(kk)
                   end do
                end do

             else
                call runend('UNKNOWN ORDER:'//trim(what))
             end if

             ipass = 0
             deallocate(tmp_irecv)
             deallocate(tmp_isend)

          end if

       end if

    end if
#endif

  end subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP

  !----------------------------------------------------------------------
  !
  ! PAR_FROM_GHOST_NODE_EXCHANGE_IP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_FROM_GHOST_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: ndofn
    real(rp),                       intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    MY_MPI_COMM   ,                   intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,jj,dom_i
    integer(ip)                                   :: nsize_send,nsize_recv
    integer(ip)                                   :: ipoin,ini_send,ini_recv,kk,idofn
    integer(4)                                    :: nsize_send4,nsize_recv4
    integer(4)                                    :: istat4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    MY_MPI_STATUS ,                    pointer    :: status4(:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_BOUNDARY_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( commu % ghost_send_node_dim == -1 .or. commu % ghost_recv_node_dim == -1 ) &
                  call runend('FRINGE NODE EXCHANGE NOT COMPUTED')
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_rsend(commu % ghost_recv_node_dim * ndofn))
             allocate(tmp_rrecv(commu % ghost_send_node_dim * ndofn))
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,commu % ghost_recv_node_dim
                ipoin = commu % ghost_recv_node_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_rsend(kk) = xx(idofn,ipoin)
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             kk = 0
             istat4 = 0_4
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   dom_i      = commd % neights(ii)
                   ini_send   = ndofn * ( commd % ghost_recv_node_size(ii)   -1 ) + 1
                   nsize_send = ndofn * ( commd % ghost_recv_node_size(ii+1) -1 ) + 1 - ini_send
                   ini_recv   = ndofn * ( commd % ghost_send_node_size(ii)   -1 ) + 1
                   nsize_recv = ndofn * ( commd % ghost_send_node_size(ii+1) -1 ) + 1 - ini_recv

                   nsize_send4 = int(nsize_send,4)
                   nsize_recv4 = int(nsize_recv,4)

                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_rsend(ini_send:ini_send+nsize_send-1), nsize_send4, &
                           PAR_REAL,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_rrecv(ini_recv:ini_recv+nsize_recv-1), nsize_recv4, &
                           PAR_REAL,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                   else
                      if( nsize_recv /= 0 .and. nsize_send == 0 ) then
                         call MPI_Recv(                          &
                              tmp_rrecv(ini_recv:), nsize_recv4, &
                              PAR_REAL, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, status, istat4    )
                      else if( nsize_recv == 0 .and. nsize_send /= 0 ) then
                         call MPI_Send(                          &
                              tmp_rsend(ini_send:), nsize_send4, &
                              PAR_REAL, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, istat4            )
                      else if( nsize_recv /= 0 .and. nsize_send /= 0 ) then
                         call MPI_Sendrecv(                         &
                              tmp_rsend(ini_send:), nsize_send4,    &
                              PAR_REAL, dom_i4, 0_4,    &
                              tmp_rrecv(ini_recv:), nsize_recv4,    &
                              PAR_REAL, dom_i4, 0_4,    &
                              PAR_COMM_TO_USE, status, istat4       )
                      end if
                   end if
                   if( istat4 /= MPI_SUCCESS )call PAR_MPI_RUNEND(istat4,'PAR_FROM_GHOST_NODE_EXCHANGE_RP')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE*2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             if( istat4 /= MPI_SUCCESS ) call runend('WRONG SEND/RECEIVE')
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,commu % ghost_send_node_dim
                   ipoin = commu % ghost_send_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_rrecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,commu % ghost_send_node_dim
                   ipoin = commu % ghost_send_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_rrecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,commu % ghost_send_node_dim
                   ipoin = commu % ghost_send_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_rrecv(kk))
                   end do
                end do

             else if( trim(what) == 'REPLACE' .or. trim(what) == 'SUBSTITUTE' ) then
                !
                ! Replace value on my fringe nodes according to what I have received
                !
                kk = 0
                do jj = 1,commu % ghost_send_node_dim
                   ipoin = commu % ghost_send_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = tmp_rrecv(kk)
                   end do
                end do

             else
                call runend('UNKNOWN ORDER:'//trim(what))
             end if

             ipass = 0
             deallocate(tmp_rrecv)
             deallocate(tmp_rsend)

          end if

       end if

    end if
#endif

  end subroutine PAR_FROM_GHOST_NODE_EXCHANGE_RP
  
  subroutine PAR_FROM_GHOST_NODE_EXCHANGE_RP_1(xx,what,wherein,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) then
       return
    else if( associated(xx) ) then
       ndofn = 1
       if( size(xx,1) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_NODE_EXCHANGE_RP_1

  subroutine PAR_FROM_GHOST_NODE_EXCHANGE_RP_2(xx,what,wherein,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) then
       return
    else if( associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_NODE_EXCHANGE_RP_2
   
  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP
  !
  !----------------------------------------------------------------------

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_00(n,xx,what,wherein,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( ISLAVE ) then
       ndofn = n
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_00

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_0(n,xx,what,wherein,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( ISLAVE ) then
       ndofn = n
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_0

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_1(xx,what,wherein,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( ISLAVE .and. associated(xx) ) then
       ndofn = 1
       if( size(xx,1) /= nboun_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_1

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_2(xx,what,wherein,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( ISLAVE .and. associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= nboun_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_2

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_3(xx,what,wherein,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( ISLAVE .and. associated(xx) ) then
       ndofn = size(xx,1)*size(xx,2)
       if( size(xx,3) /= nboun_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_3

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_2b(xx,what,commu,wsynch,dom_k)
    
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( ISLAVE .and. associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= nboun_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
       call PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_2b

  !----------------------------------------------------------------------
  !
  ! PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: ndofn
    real(rp),                       intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    MY_MPI_COMM   ,                   intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,jj,dom_i
    integer(ip)                                   :: nsize_send,nsize_recv
    integer(ip)                                   :: ipoin,ini_send,ini_recv,kk,idofn
    integer(4)                                    :: nsize_send4,nsize_recv4
    integer(4)                                    :: istat4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    MY_MPI_STATUS ,                    pointer    :: status4(:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_BOUNDARY_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( commu % ghost_send_boun_dim == -1 .or. commu % ghost_recv_boun_dim == -1 ) &
                  call runend('FRINGE NODE EXCHANGE NOT COMPUTED')
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_rsend(commu % ghost_recv_boun_dim * ndofn))
             allocate(tmp_rrecv(commu % ghost_send_boun_dim * ndofn))
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,commu % ghost_recv_boun_dim
                ipoin = commu % ghost_recv_boun_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_rsend(kk) = xx(idofn,ipoin)
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             kk = 0
             istat4 = 0_4
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   dom_i      = commd % neights(ii)
                   ini_send   = ndofn * ( commd % ghost_recv_boun_size(ii)   -1 ) + 1
                   nsize_send = ndofn * ( commd % ghost_recv_boun_size(ii+1) -1 ) + 1 - ini_send
                   ini_recv   = ndofn * ( commd % ghost_send_boun_size(ii)   -1 ) + 1
                   nsize_recv = ndofn * ( commd % ghost_send_boun_size(ii+1) -1 ) + 1 - ini_recv

                   nsize_send4 = int(nsize_send,4)
                   nsize_recv4 = int(nsize_recv,4)

                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_rsend(ini_send:ini_send+nsize_send-1), nsize_send4, &
                           PAR_REAL,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_rrecv(ini_recv:ini_recv+nsize_recv-1), nsize_recv4, &
                           PAR_REAL,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                   else
                      if( nsize_recv /= 0 .and. nsize_send == 0 ) then
                         call MPI_Recv(                          &
                              tmp_rrecv(ini_recv:), nsize_recv4, &
                              PAR_REAL, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, status, istat4    )
                      else if( nsize_recv == 0 .and. nsize_send /= 0 ) then
                         call MPI_Send(                          &
                              tmp_rsend(ini_send:), nsize_send4, &
                              PAR_REAL, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, istat4            )
                      else if( nsize_recv /= 0 .and. nsize_send /= 0 ) then
                         call MPI_Sendrecv(                         &
                              tmp_rsend(ini_send:), nsize_send4,    &
                              PAR_REAL, dom_i4, 0_4,    &
                              tmp_rrecv(ini_recv:), nsize_recv4,    &
                              PAR_REAL, dom_i4, 0_4,    &
                              PAR_COMM_TO_USE, status, istat4       )
                      end if
                   end if
                   if( istat4 /= MPI_SUCCESS )call PAR_MPI_RUNEND(istat4,'PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE*2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             if( istat4 /= MPI_SUCCESS ) call runend('WRONG SEND/RECEIVE')
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,commu % ghost_send_boun_dim
                   ipoin = commu % ghost_send_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_rrecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,commu % ghost_send_boun_dim
                   ipoin = commu % ghost_send_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_rrecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,commu % ghost_send_boun_dim
                   ipoin = commu % ghost_send_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_rrecv(kk))
                   end do
                end do

             else if( trim(what) == 'REPLACE' .or. trim(what) == 'SUBSTITUTE' ) then
                !
                ! Replace value on my fringe nodes according to what I have received
                !
                kk = 0
                do jj = 1,commu % ghost_send_boun_dim
                   ipoin = commu % ghost_send_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = tmp_rrecv(kk)
                   end do
                end do

             else
                call runend('UNKNOWN ORDER:'//trim(what))
             end if

             ipass = 0
             deallocate(tmp_rrecv)
             deallocate(tmp_rsend)

          end if

       end if

    end if
#endif

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP
  !
  !----------------------------------------------------------------------

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_00(n,xx,what,wherein,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: n
    integer(ip),                    intent(inout) :: xx(*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( ISLAVE  ) then
       ndofn = n
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_00

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_0(n,xx,what,wherein,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: n
    integer(ip),                    intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( ISLAVE ) then
       ndofn = n
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_0

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_1(xx,what,wherein,wsynch,dom_k)
    
    integer(ip),          pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( ISLAVE .and. associated(xx) ) then
       ndofn = 1
       if( size(xx,1) /= nboun_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_1

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_2(xx,what,wherein,wsynch,dom_k)
    
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( ISLAVE .and. associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= nboun_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_2

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_3(xx,what,wherein,wsynch,dom_k)
    
    integer(ip),          pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( ISLAVE .and. associated(xx) ) then
       ndofn = size(xx,1)*size(xx,2)
       if( size(xx,3) /= nboun_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_3

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_2b(xx,what,commu,wsynch,dom_k)
    
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE

    if( ISLAVE .and. associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= nboun_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
       call PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_2b

  !----------------------------------------------------------------------
  !
  ! PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: ndofn
    integer(ip),                    intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    MY_MPI_COMM   ,                   intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,jj,dom_i
    integer(ip)                                   :: nsize_send,nsize_recv
    integer(ip)                                   :: ipoin,ini_send,ini_recv,kk,idofn
    integer(4)                                    :: nsize_send4,nsize_recv4
    integer(4)                                    :: istat4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    MY_MPI_STATUS ,                    pointer    :: status4(:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_BOUNDARY_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( commu % ghost_send_boun_dim == -1 .or. commu % ghost_recv_boun_dim == -1 ) &
                  call runend('FRINGE NODE EXCHANGE NOT COMPUTED')
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_isend(commu % ghost_recv_boun_dim * ndofn))
             allocate(tmp_irecv(commu % ghost_send_boun_dim * ndofn))
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,commu % ghost_recv_boun_dim
                ipoin = commu % ghost_recv_boun_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_isend(kk) = xx(idofn,ipoin)
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             kk = 0
             istat4 = 0_4
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   dom_i      = commd % neights(ii)
                   ini_send   = ndofn * ( commd % ghost_recv_boun_size(ii)   -1 ) + 1
                   nsize_send = ndofn * ( commd % ghost_recv_boun_size(ii+1) -1 ) + 1 - ini_send
                   ini_recv   = ndofn * ( commd % ghost_send_boun_size(ii)   -1 ) + 1
                   nsize_recv = ndofn * ( commd % ghost_send_boun_size(ii+1) -1 ) + 1 - ini_recv

                   nsize_send4 = int(nsize_send,4)
                   nsize_recv4 = int(nsize_recv,4)

                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_isend(ini_send:ini_send+nsize_send-1), nsize_send4, &
                           PAR_INTEGER,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_irecv(ini_recv:ini_recv+nsize_recv-1), nsize_recv4, &
                           PAR_INTEGER,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                   else
                      if( nsize_recv /= 0 .and. nsize_send == 0 ) then
                         call MPI_Recv(                          &
                              tmp_irecv(ini_recv:), nsize_recv4, &
                              PAR_INTEGER, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, status, istat4    )
                      else if( nsize_recv == 0 .and. nsize_send /= 0 ) then
                         call MPI_Send(                          &
                              tmp_isend(ini_send:), nsize_send4, &
                              PAR_INTEGER, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, istat4            )
                      else if( nsize_recv /= 0 .and. nsize_send /= 0 ) then
                         call MPI_Sendrecv(                         &
                              tmp_isend(ini_send:), nsize_send4,    &
                              PAR_INTEGER, dom_i4, 0_4,    &
                              tmp_irecv(ini_recv:), nsize_recv4,    &
                              PAR_INTEGER, dom_i4, 0_4,    &
                              PAR_COMM_TO_USE, status, istat4       )
                      end if
                   end if
                   if( istat4 /= MPI_SUCCESS )call PAR_MPI_RUNEND(istat4,'PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE*2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             if( istat4 /= MPI_SUCCESS ) call runend('WRONG SEND/RECEIVE')
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,commu % ghost_send_boun_dim
                   ipoin = commu % ghost_send_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_irecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,commu % ghost_send_boun_dim
                   ipoin = commu % ghost_send_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,commu % ghost_send_boun_dim
                   ipoin = commu % ghost_send_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'REPLACE' .or. trim(what) == 'SUBSTITUTE' ) then
                !
                ! Replace value on my fringe nodes according to what I have received
                !
                kk = 0
                do jj = 1,commu % ghost_send_boun_dim
                   ipoin = commu % ghost_send_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = tmp_irecv(kk)
                   end do
                end do

             else
                call runend('UNKNOWN ORDER:'//trim(what))
             end if

             ipass = 0
             deallocate(tmp_irecv)
             deallocate(tmp_isend)

          end if

       end if

    end if
#endif

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP
  !
  !----------------------------------------------------------------------

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_00(n,xx,what,wherein,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: n
    integer(ip),                       intent(inout) :: xx(*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_00

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_0(n,xx,what,wherein,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: n
    integer(ip),                       intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_0

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_1(xx,what,wherein,wsynch,dom_k)
    
    integer(ip),             pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE .or. .not. associated(xx) ) return
    if( associated(xx) ) then
       ndofn = 1
       if( size(xx,1) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_1

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_2(xx,what,wherein,wsynch,dom_k)
    
    integer(ip),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE .or. .not. associated(xx) ) return
    ndofn = size(xx,1)
    if( size(xx,2) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_2

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_3(xx,what,wherein,wsynch,dom_k)
    
    integer(ip),             pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                         :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE .or. .not. associated(xx) ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)*size(xx,2)
       if( size(xx,3) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_3

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_2b(xx,what,commu,wsynch,dom_k)
    
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    MY_MPI_COMM                                   :: PAR_COMM_TO_USE

    if( INOTSLAVE .or. .not. associated(xx) ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
       call PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_2b

  !----------------------------------------------------------------------
  !
  ! PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    
    integer(ip),                    intent(in)    :: ndofn
    integer(ip),                    intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    MY_MPI_COMM,                    intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,jj,dom_i
    integer(ip)                                   :: nsize_send,nsize_recv
    integer(ip)                                   :: ipoin,ini_send,ini_recv,kk,idofn
    integer(4)                                    :: nsize_send4,nsize_recv4
    integer(4)                                    :: istat4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    MY_MPI_STATUS,                    pointer    :: status4(:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_ELEMENT_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( commu % ghost_send_elem_dim == -1 .or. commu % ghost_recv_elem_dim == -1 ) &
                  call runend('FRINGE NODE EXCHANGE NOT COMPUTED')
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_isend(commu % ghost_recv_elem_dim * ndofn))
             allocate(tmp_irecv(commu % ghost_send_elem_dim * ndofn))
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,commu % ghost_recv_elem_dim
                ipoin = commu % ghost_recv_elem_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_isend(kk) = xx(idofn,ipoin)
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             kk = 0
             istat4 = 0_4
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   dom_i      = commd % neights(ii)
                   ini_send   = ndofn * ( commd % ghost_recv_elem_size(ii)   -1 ) + 1
                   nsize_send = ndofn * ( commd % ghost_recv_elem_size(ii+1) -1 ) + 1 - ini_send
                   ini_recv   = ndofn * ( commd % ghost_send_elem_size(ii)   -1 ) + 1
                   nsize_recv = ndofn * ( commd % ghost_send_elem_size(ii+1) -1 ) + 1 - ini_recv

                   nsize_send4 = int(nsize_send,4)
                   nsize_recv4 = int(nsize_recv,4)

                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_isend(ini_send:ini_send+nsize_send-1), nsize_send4, &
                           PAR_INTEGER,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_irecv(ini_recv:ini_recv+nsize_recv-1), nsize_recv4, &
                           PAR_INTEGER,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                   else
                      if( nsize_recv /= 0 .and. nsize_send == 0 ) then
                         call MPI_Recv(                          &
                              tmp_irecv(ini_recv:), nsize_recv4, &
                              PAR_INTEGER, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, status, istat4    )
                      else if( nsize_recv == 0 .and. nsize_send /= 0 ) then
                         call MPI_Send(                          &
                              tmp_isend(ini_send:), nsize_send4, &
                              PAR_INTEGER, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, istat4            )
                      else if( nsize_recv /= 0 .and. nsize_send /= 0 ) then
                         call MPI_Sendrecv(                         &
                              tmp_isend(ini_send:), nsize_send4,    &
                              PAR_INTEGER, dom_i4, 0_4,    &
                              tmp_irecv(ini_recv:), nsize_recv4,    &
                              PAR_INTEGER, dom_i4, 0_4,    &
                              PAR_COMM_TO_USE, status, istat4       )
                      end if
                   end if
                   if( istat4 /= MPI_SUCCESS )call PAR_MPI_RUNEND(istat4,' PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE*2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             if( istat4 /= MPI_SUCCESS ) call runend('WRONG SEND/RECEIVE')
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,commu % ghost_send_elem_dim
                   ipoin = commu % ghost_send_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_irecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,commu % ghost_send_elem_dim
                   ipoin = commu % ghost_send_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,commu % ghost_send_elem_dim
                   ipoin = commu % ghost_send_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'REPLACE' .or. trim(what) == 'SUBSTITUTE' ) then
                !
                ! Replace value on my fringe nodes according to what I have received
                !
                kk = 0
                do jj = 1,commu % ghost_send_elem_dim
                   ipoin = commu % ghost_send_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = tmp_irecv(kk)
                   end do
                end do

             else
                call runend('UNKNOWN ORDER: '//trim(what))
             end if

             ipass = 0
             deallocate(tmp_irecv)
             deallocate(tmp_isend)

          end if

       end if

    end if
#endif

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-12-09
  !> @brief   Get mon rank
  !> @details Min rank of a node owner
  !> 
  !-----------------------------------------------------------------------

  integer(ip) function PAR_MIN_RANK_OWNER(COMM,ipoin) 
    
    type(comm_data_par), intent(inout) :: COMM            !< Input communicator
    integer(ip),         intent(in)    :: ipoin
    integer(ip)                        :: ii,ineig,dom_i
    integer(ip)                        :: jpoin

    PAR_MIN_RANK_OWNER = huge(1_ip)
    do ineig = 1,COMM % NNEIG
       dom_i = COMM % NEIGHTS(ineig)
       do ii = COMM % bound_size(ineig),COMM % bound_size(ineig+1)-1
          jpoin = COMM % bound_perm(ii)
          if( ipoin == jpoin ) PAR_MIN_RANK_OWNER = min(PAR_MIN_RANK_OWNER,dom_i)
       end do
    end do
    
  end function PAR_MIN_RANK_OWNER
    
end module mod_communications_point_to_point
!> @}


