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
!> @file    mod_communications.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!> @brief   ToolBox for parallel communications
!> @details ToolBox for parallel communications
!------------------------------------------------------------------------

module mod_communications

  use def_communications
  use mod_communications_tools
  use mod_communications_global
  use mod_communications_point_to_point
  implicit none  

  private
  !
  ! From mod_communications_tools
  !
  public :: PAR_DEFINE_COMMUNICATOR               ! Define the communicator according to some keywords
  public :: PAR_COMM_RANK_AND_SIZE                ! Give rank (and size) of a communicator
  public :: PAR_INIT                              ! Initialize MPI
  public :: PAR_COMM_SPLIT                        ! Split a communicator
  public :: PAR_COMM_FREE                         ! Free MPI communicator
  public :: PAR_BARRIER                           ! Barrier
  public :: PAR_MPI_ERROR_TO_MESSAGE              ! Transform an MPI error code into a string
  !
  ! From mod_communications_global
  !
  public :: PAR_SUM                               ! AllReduce SUM
  public :: PAR_MAX                               ! AllReduce MAX
  public :: PAR_MIN                               ! AllReduce MIN
  public :: PAR_AVERAGE                           ! Average value
  public :: PAR_LOAD_BALANCE                      ! Load balance: average / max
  public :: PAR_SUM_ALL                           ! AllReduce SUM in the Alya world (PAR_COMM_WORLD)
  public :: PAR_MAX_ALL                           ! AllReduce MAX in the Alya world (PAR_COMM_WORLD)
  public :: PAR_ALLTOALL                          ! Alltoall
  public :: PAR_ALLTOALLV                         ! All to all V
  public :: PAR_ALLTOALLV_RP_1
  public :: PAR_ALLTOALLV_RP_2
  public :: PAR_ALLTOALLV_RP_4
  public :: PAR_ALLTOALLV_RP_3
  public :: PAR_ALLTOALLV_IP_1
  public :: PAR_ALLTOALLV_IP_2
  public :: PAR_BROADCAST                         ! Broadcast
  public :: PAR_ALL_TO_ALL_ARRAY_OPERATION        ! Array operations between all partitions of communicators
  public :: PAR_GATHER                            ! Gather of integers, reals, characters
  public :: PAR_GATHERV                           ! Gather of integers, reals, characters
  public :: PAR_ALLGATHERV                        ! All Gatherv
  public :: PAR_ALLGATHER                         ! All Gather
  public :: PAR_SCATTER                           ! Scatter of integers, reals, characters
  public :: PAR_SCATTERV                          ! Scatter of integers, reals, characters
  !
  ! From mod_communications_point_to_point
  !
  public :: PAR_INTERFACE_NODE_EXCHANGE           ! Interface nodal exchange (send/receive)
  public :: PAR_INTERFACE_OWN_NODE_EXCHANGE       ! Interface OWN node exchange (send/receive)
  public :: PAR_INTERFACE_EDGE_EXCHANGE           ! Interface edge exchange (send/receive)
  public :: PAR_INTERFACE_MATRIX_EXCHANGE         ! Matrix exchange on interface nodes
  public :: PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE ! Matrix with halos exchange on interface nodes
  public :: PAR_COUPLING_NODE_EXCHANGE            ! Coupling nodal exchange (send/receive)
  public :: PAR_GHOST_ELEMENT_EXCHANGE            ! Ghost element exchange (send/receive): from element to ghost
  public :: PAR_GHOST_NODE_EXCHANGE               ! Ghost nodal exchange (send/receive)
  public :: PAR_FROM_GHOST_ELEMENT_EXCHANGE       ! Ghost element exchange (send/receive): from ghost to element
  public :: PAR_FROM_GHOST_BOUNDARY_EXCHANGE      ! Ghost element exchange (send/receive): from ghost to boundary
  public :: PAR_FROM_GHOST_NODE_EXCHANGE          ! Ghost node exchange (send/receive): from ghost to node
  public :: PAR_GHOST_BOUNDARY_EXCHANGE           ! Ghost boundary exchange (send/receive)
  public :: PAR_SEND                              ! Send arrays to a specific partition
  public :: PAR_RECEIVE                           ! Receive arrays to a specific partition
  public :: PAR_SEND_RECEIVE                      ! Send and receive arrays to a specific partition
  public :: PAR_SEND_RECEIVE_IP                   ! Send and receive arrays to a specific partition
  public :: PAR_SEND_RECEIVE_RP                   ! Send and receive arrays to a specific partition
  public :: PAR_SEND_RECEIVE_TO_ALL               ! Send and receive arrays to all
  public :: PAR_EXCHANGE                          ! Bridge to broadcast data from master to slave
  public :: PAR_START_NON_BLOCKING_COMM           ! Start non-blocking communications
  public :: PAR_END_NON_BLOCKING_COMM             ! End non-blocking communications
  public :: PAR_SET_NON_BLOCKING_COMM_NUMBER      ! Set the non-blocking communicator number
  public :: PAR_WAITALL                           ! Waitall
  public :: PAR_POINT_TO_POINT_ARRAY_OPERATION    ! Array operations between all partitions of communicators
  public :: PAR_COMM_NULL                         ! NULL MPI communicator
  public :: PAR_MIN_RANK_OWNER                    ! Get the minimum owner rank of an interface node

  public :: communications_initialization         ! Initialization

contains

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-12-29
  !> @brief   Communication initialization
  !> @details Communication initialization
  !>
  !-----------------------------------------------------------------------

  subroutine communications_initialization()
    !
    ! Local arrays
    !
    call communications_point_to_point_initialization()
    
  end subroutine communications_initialization

end module mod_communications
!> @}


