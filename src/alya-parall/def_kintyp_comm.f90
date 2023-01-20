!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kinds_and_types
!> @{
!> @file    def_kintyp_comm.g90
!> @author  houzeaux
!> @date    2020-04-04
!> @brief   Functions
!> @details Communications
!-----------------------------------------------------------------------

module def_kintyp_comm

  use def_kintyp_basic,      only : ip,rp,lg,i1p
  use mod_memory_basic,      only : memory_alloca
  use mod_memory_basic,      only : memory_deallo
  use mod_memory_basic,      only : memory_copy
  use mod_memory_basic,      only : memory_append
  use mod_memory_basic,      only : memory_resize
  use mod_memory_tools,      only : memory_counter_ini
  use mod_memory_tools,      only : memory_counter_end
  use mod_optional_argument, only : optional_argument
  use mod_std 
  use def_mpi
#include "def_mpi.inc"
 
  private

  character(15), parameter :: vacal = 'def_kintyp_comm'
  integer(ip),   parameter :: COMM_SEND_RECV = 0
  integer(ip),   parameter :: COMM_ONE_SIDED = 1
  
  type comm_bound_matrix
     type(i1p),   pointer :: ja(:)
     integer(ip), pointer :: nzdom_ii(:)
     integer(ip)          :: nzdom
  end type comm_bound_matrix
  !
  ! Automatic partitioning based on parallel efficiency
  !
  type typ_optimum_npart
     integer(ip) :: kfl_method
     integer(ip) :: frequency
     integer(ip) :: modul
     integer(ip) :: min_cores
     integer(ip) :: max_cores
     real(rp)    :: min_criterion     
     real(rp)    :: max_criterion
     real(rp)    :: change_rate
  end type typ_optimum_npart
  
  type comm_data_par_basic
     ! Name
     character(30)                    :: name                   ! Name of the coomunicator
     ! Memory
     integer(8)                       :: memor(2)               ! Memory counter
     ! Neighbors
     integer(ip)                      :: nneig                  ! Number of neigbors
     integer(ip), pointer             :: neights(:)             ! List of neighbors
     ! Interface node communication
     integer(ip)                      :: bound_dim              ! Size of interface
     integer(ip), pointer             :: bound_size(:)          ! Size of interface node esxchange
     integer(ip), pointer             :: bound_perm(:)          ! List of subdomain interface nodes
     ! Unsymmetric send-receive
     integer(ip), pointer             :: lsend_size(:)
     integer(ip), pointer             :: lsend_perm(:)
     integer(ip)                      :: lsend_dim
     integer(ip), pointer             :: lrecv_size(:)
     integer(ip), pointer             :: lrecv_perm(:)
     integer(ip)                      :: lrecv_dim
     ! Offsets
     integer(ip)                      :: offset_npoin           ! Offset of nodes
     integer(ip)                      :: offset_nelem           ! Offset of elements
     ! Communicator
     MY_MPI_COMM                      :: PAR_COMM_WORLD         ! Communicator
     integer(4)                       :: SIZE4                  ! Size 
     integer(4)                       :: RANK4                  ! Rank
      ! Interior, own and other boundary nodes
     integer(ip)                      :: npoi1
     integer(ip)                      :: npoi2
     integer(ip)                      :: npoi3
     integer(ip)                      :: npoin
     ! Type of communications
     integer(ip)                      :: type                   ! send/recv or one-sided 
     ! One-sided communications
     MY_MPI_WIN                       :: WIN                    ! Window with master
     MY_MPI_WIN                       :: WIN_WM                 ! Window without master
     integer(ip), pointer             :: displ_recv(:)          ! Displacement for receiver
     real(rp),    pointer, contiguous :: buffe_recv(:)          ! Buffer for receiver
  contains
     procedure,                  pass :: init   => init_basic   ! Initialize
     procedure,                  pass :: alloca => alloca_basic ! Allocate a communicator
     procedure,                  pass :: deallo => deallo_basic ! Deallocate
     procedure,                  pass :: copy   => copy_basic   ! Copy a communicator
     procedure,                  pass :: move                   ! Move a communicator
     procedure,                  pass :: merge                  ! Merge two communicators
     procedure,                  pass :: null_arrays            ! Nullify arrays
     procedure,                  pass :: permute                ! Permute arrays
  end type comm_data_par_basic
  
  type, extends(comm_data_par_basic)  :: comm_data_par
     ! Ordered neighbors
     integer(ip)                      :: nneig_1
     integer(ip)                      :: nneig_2
     integer(ip), pointer             :: neights_ordered(:)
     ! Interior, own and other boundary edges
     integer(ip)                      :: nedg1
     integer(ip)                      :: nedg2
     integer(ip)                      :: nedg3
     ! Offsets
     integer(ip)                      :: offset_nboun           ! Offset of boundaries     
     ! Interface node communication
     integer(ip), pointer             :: bound_invp(:)
     integer(ip), pointer             :: bound_scal(:)
     integer(ip), pointer             :: bound_multiplicity(:)
     integer(ip), pointer             :: bound_owner_rank(:)
     integer(ip), pointer             :: node_number_in_owner(:)
     integer(ip), pointer             :: perm_ordered(:)
     type(comm_bound_matrix), pointer :: bound_matrix(:)
     type(comm_bound_matrix), pointer :: bound_mat_halo_send(:)
     type(comm_bound_matrix), pointer :: bound_mat_halo_recv(:)
     ! Interface edge communication
     integer(ip), pointer             :: bedge_size(:)
     integer(ip), pointer             :: bedge_perm(:)
     integer(ip), pointer             :: bedge_adja(:)
     integer(ip), pointer             :: bedge_scal(:)
     integer(ip)                      :: bedge_dim
     integer(ip), pointer             :: bedge_multiplicity(:)
     integer(ip), pointer             :: bedge_owner_rank(:)
     ! Non-symmetric send receive
     integer(ip), pointer             :: lscat_perm(:)
     integer(ip)                      :: lscat_dim
     integer(ip)                      :: matrix_nzdom
     integer(ip), pointer             :: matrix_ia(:)
     integer(ip), pointer             :: matrix_ja(:)
     real(rp),    pointer             :: matrix_aa(:)
     ! Face communication
     integer(ip), pointer             :: bface_size(:)
     integer(ip), pointer             :: bface_perm(:)
     integer(ip)                      :: bface_dim
     ! Full row matrix communicator
     integer(ip)                      :: full_row_send_nneig
     integer(ip), pointer             :: full_row_send_neights(:)
     integer(ip), pointer             :: full_row_send_size(:)
     integer(ip), pointer             :: full_row_send_perm(:)
     integer(ip)                      :: full_row_send_dim
     integer(ip)                      :: full_row_recv_nneig
     integer(ip), pointer             :: full_row_recv_neights(:)
     integer(ip), pointer             :: full_row_recv_size(:)
     integer(ip), pointer             :: full_row_recv_perm(:)
     integer(ip)                      :: full_row_recv_dim
     ! Ghost node communication
     integer(ip), pointer             :: ghost_send_node_size(:)
     integer(ip), pointer             :: ghost_send_node_perm(:)
     integer(ip)                      :: ghost_send_node_dim
     integer(ip), pointer             :: ghost_recv_node_size(:)
     integer(ip), pointer             :: ghost_recv_node_perm(:)
     integer(ip)                      :: ghost_recv_node_dim
     ! Ghost element communication
     integer(ip), pointer             :: ghost_send_elem_size(:)
     integer(ip), pointer             :: ghost_send_elem_perm(:)
     integer(ip)                      :: ghost_send_elem_dim
     integer(ip), pointer             :: ghost_recv_elem_size(:)
     integer(ip), pointer             :: ghost_recv_elem_perm(:)
     integer(ip)                      :: ghost_recv_elem_dim
     ! Ghost boundary communication
     integer(ip), pointer             :: ghost_send_boun_size(:)
     integer(ip), pointer             :: ghost_send_boun_perm(:)
     integer(ip)                      :: ghost_send_boun_dim
     integer(ip), pointer             :: ghost_recv_boun_size(:)
     integer(ip), pointer             :: ghost_recv_boun_perm(:)
     integer(ip)                      :: ghost_recv_boun_dim
   contains
     procedure,                  pass :: init   => init_full    ! Initialize
     procedure,                  pass :: alloca => alloca_full  ! Allocate
     procedure,                  pass :: deallo => deallo_full  ! Deallocate
     procedure,                  pass :: copy   => copy_full    ! Copy a communicator
     procedure,                  pass :: sub                    ! Compute a sub communicator using a msk
     procedure,                  pass :: collapse               ! Collapse list of neighbors
  end type comm_data_par
  
  type :: tAdj_par
     integer(ip)            :: node1
     integer(ip)            :: node2
  end type tAdj_par

  public :: typ_optimum_npart
  public :: comm_data_par_basic
  public :: comm_data_par
  public :: tAdj_par
  public :: COMM_SEND_RECV 
  public :: COMM_ONE_SIDED 
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-28
  !> @brief   Initialization
  !> @details Initialzation of a communicator
  !> 
  !-----------------------------------------------------------------------

  subroutine init_basic(comm,PAR_COMM_OPT,COMM_NAME)

    class(comm_data_par_basic), intent(inout)        :: comm
    logical(lg),                intent(in), optional :: PAR_COMM_OPT
    character(*),               intent(in), optional :: COMM_NAME
    logical(lg)                                      :: PAR_comm
    
    PAR_comm              = optional_argument(.true.,PAR_comm_OPT)
    
    comm % name           = optional_argument('COMM',COMM_NAME)
    comm % memor          = 0_8
    comm % nneig          = 0
    comm % bound_dim      = 0
    comm % lsend_dim      = 0
    comm % lrecv_dim      = 0
    comm % offset_npoin   = 0
    comm % offset_nelem   = 0
    if(  PAR_comm ) comm % PAR_COMM_WORLD = PAR_COMM_NULL
    comm % SIZE4          = 0_4
    comm % RANK4          = 0_4
    comm % npoi1          = 0
    comm % npoi2          = 0
    comm % npoi3          = 0
    comm % npoin          = 0
    comm % type           = COMM_SEND_RECV
    comm % WIN            = MPI_WIN_NULL
    comm % WIN_WM         = MPI_WIN_NULL
    
    nullify(comm % neights   )
    nullify(comm % bound_size)
    nullify(comm % bound_perm)
    nullify(comm % lsend_size)
    nullify(comm % lsend_perm)
    nullify(comm % lrecv_size)
    nullify(comm % lrecv_perm)
    nullify(comm % displ_recv)
    nullify(comm % buffe_recv)

  end subroutine init_basic

  subroutine init_full(comm,PAR_COMM_OPT,COMM_NAME)

    class(comm_data_par),      intent(inout)        :: comm
    logical(lg),               intent(in), optional :: PAR_COMM_OPT
    character(*),              intent(in), optional :: COMM_NAME

    call init_basic(comm,PAR_COMM_OPT,COMM_NAME)
    
    comm % nneig_1             =  0
    comm % nneig_2             =  0
    comm % nedg1               =  0
    comm % nedg2               =  0
    comm % nedg3               =  0
    
    comm % offset_nboun        =  0
    
    comm % bedge_dim           =  0
    comm % lscat_dim           =  0
    comm % matrix_nzdom        =  0
    comm % bface_dim           =  0

    comm % full_row_send_nneig =  0
    comm % full_row_send_dim   = -1
    comm % full_row_recv_nneig =  0
    comm % full_row_recv_dim   = -1

    comm % ghost_send_elem_dim = -1
    comm % ghost_recv_elem_dim = -1
    comm % ghost_send_node_dim = -1
    comm % ghost_recv_node_dim = -1
    comm % ghost_send_boun_dim = -1
    comm % ghost_recv_boun_dim = -1

    nullify(comm % neights_ordered)
    nullify(comm % perm_ordered)
    nullify(comm % bound_invp)
    nullify(comm % bound_scal)
    nullify(comm % bound_multiplicity)
    nullify(comm % bound_owner_rank)
    nullify(comm % node_number_in_owner)
    nullify(comm % bound_matrix)
    nullify(comm % bound_mat_halo_send)
    nullify(comm % bound_mat_halo_recv)

    nullify(comm % bedge_size)
    nullify(comm % bedge_perm)
    nullify(comm % bedge_adja)
    nullify(comm % bedge_scal)
    nullify(comm % bedge_multiplicity)
    nullify(comm % bedge_owner_rank)

    nullify(comm % lscat_perm)
    nullify(comm % matrix_ia)
    nullify(comm % matrix_ja)
    nullify(comm % matrix_aa)
    nullify(comm % bface_size)
    nullify(comm % bface_perm)

    nullify(comm % full_row_send_neights)
    nullify(comm % full_row_send_size)
    nullify(comm % full_row_send_perm)
    nullify(comm % full_row_recv_neights)
    nullify(comm % full_row_recv_size)
    nullify(comm % full_row_recv_perm)

    nullify(comm % ghost_send_elem_size)
    nullify(comm % ghost_send_elem_perm)
    nullify(comm % ghost_recv_elem_size)
    nullify(comm % ghost_recv_elem_perm)
    nullify(comm % ghost_send_node_size)
    nullify(comm % ghost_send_node_perm)
    nullify(comm % ghost_recv_node_size)
    nullify(comm % ghost_recv_node_perm)
    nullify(comm % ghost_send_boun_size)
    nullify(comm % ghost_send_boun_perm)
    nullify(comm % ghost_recv_boun_size)
    nullify(comm % ghost_recv_boun_perm)
    
  end subroutine init_full

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-28
  !> @brief   Nullify
  !> @details Nullify arrays
  !> 
  !-----------------------------------------------------------------------

  subroutine null_arrays(comm)

    class(comm_data_par_basic) :: comm

    nullify(comm % neights   )
    nullify(comm % bound_size)
    nullify(comm % bound_perm)
    nullify(comm % lsend_size)
    nullify(comm % lsend_perm)
    nullify(comm % lrecv_size)
    nullify(comm % lrecv_perm)
    nullify(comm % displ_recv)
    nullify(comm % buffe_recv)

  end subroutine null_arrays

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-28
  !> @brief   Deallocate
  !> @details Deallocate arrays
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo_basic(comm,MEMORY_COUNTER,PAR_comm_OPT,comm_NAME,INITIALIZE)

    class(comm_data_par_basic), intent(inout)           :: comm          !< Input communicator
    integer(8),                 intent(inout), optional :: MEMORY_COUNTER(2)  !< Memory counter
    logical(lg),                intent(in),    optional :: PAR_COMM_OPT
    character(*),               intent(in),    optional :: COMM_NAME
    logical(lg),                intent(in),    optional :: INITIALIZE
    integer(8)                                          :: memor(2)
    logical(lg)                                         :: PAR_comm
    character(LEN=:), allocatable                       :: my_comm_name
    logical(lg)                                         :: if_initialize

    memor         = optional_argument((/0_8,0_8/),MEMORY_COUNTER)

    if_initialize = optional_argument(.false.,INITIALIZE)
    memor         = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
    my_comm_name  = optional_argument(trim(comm % name),COMM_NAME)
    PAR_comm      = optional_argument(.true.,PAR_comm_OPT)

    if( if_initialize ) then
       comm % nneig               = 0
       comm % bound_dim           = 0
       comm % npoi1               = 0
       comm % npoi2               = 0
       comm % npoi3               = 0
       comm % npoin               = 0
      
       comm % offset_npoin        = 0
       comm % offset_nelem        = 0
 
       if(  PAR_comm ) comm % PAR_comm_WORLD = PAR_COMM_NULL

       comm % SIZE4               = -1_4
       comm % RANK4               = -1_4

    end if

    call memory_deallo(memor,trim(my_comm_name)//' % NEIGHTS'                ,vacal,comm % neights    )
    call memory_deallo(memor,trim(my_comm_name)//' % BOUND_SIZE'             ,vacal,comm % bound_size )
    call memory_deallo(memor,trim(my_comm_name)//' % BOUND_PERM'             ,vacal,comm % bound_perm )
    call memory_deallo(memor,trim(my_comm_name)//' % LSEND_SIZE'             ,vacal,comm % lsend_size )
    call memory_deallo(memor,trim(my_comm_name)//' % LSEND_PERM'             ,vacal,comm % lsend_perm )
    call memory_deallo(memor,trim(my_comm_name)//' % LRECV_SIZE'             ,vacal,comm % lrecv_size )
    call memory_deallo(memor,trim(my_comm_name)//' % LRECV_PERM'             ,vacal,comm % lrecv_perm )    
    call memory_deallo(memor,trim(my_comm_name)//' % DISPL_RECV'             ,vacal,comm % displ_recv )    

    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor
    if( allocated(my_comm_name) ) deallocate(my_comm_name)

  end subroutine deallo_basic

  subroutine deallo_full(comm,MEMORY_COUNTER,PAR_comm_OPT,comm_NAME,INITIALIZE)

    class(comm_data_par),       intent(inout)           :: comm          !< Input communicator
    integer(8),                 intent(inout), optional :: MEMORY_COUNTER(2)  !< Memory counter
    logical(lg),                intent(in),    optional :: PAR_comm_OPT
    character(*),               intent(in),    optional :: COMM_NAME
    logical(lg),                intent(in),    optional :: INITIALIZE
    integer(8)                                          :: memor(2)
    logical(lg)                                         :: PAR_comm
    character(LEN=:), allocatable                       :: my_comm_name
    logical(lg)                                         :: if_initialize
    integer(ip)                                         :: ii
    
    
    if_initialize = optional_argument(.false.,INITIALIZE)
    memor         = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
    my_comm_name  = optional_argument(trim(comm % name),comm_NAME)
    PAR_comm      = optional_argument(.true.,PAR_comm_OPT)

    call deallo_basic(comm,memor,PAR_comm_OPT,comm_NAME,INITIALIZE)
    
    if( if_initialize ) then
  
       comm % nedg1               = 0
       comm % nedg2               = 0
       comm % nedg3               = 0

       comm % bedge_dim           = 0
       comm % lsend_dim           = 0
       comm % lrecv_dim           = 0
       comm % lscat_dim           = 0
       comm % matrix_nzdom        = 0
       comm % bface_dim           = 0
       comm % full_row_send_dim   = 0
       comm % full_row_recv_nneig = 0
       comm % full_row_recv_dim   = 0
       comm % ghost_send_node_dim = 0
       comm % ghost_recv_node_dim = 0
       comm % ghost_send_elem_dim = 0
       comm % ghost_recv_elem_dim = 0
       comm % ghost_send_boun_dim = 0
       
    end if

    call memory_deallo(memor,trim(my_comm_name)//' % NEIGHTS_ORDERED'        ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % neights_ordered      )
    call memory_deallo(memor,trim(my_comm_name)//' % BOUND_INVP'             ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % bound_invp           )
    call memory_deallo(memor,trim(my_comm_name)//' % PERM_ORDERED'           ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % perm_ordered         )
    call memory_deallo(memor,trim(my_comm_name)//' % FULL_ROW_SEND_NEIGHTS'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % full_row_send_neights)
    call memory_deallo(memor,trim(my_comm_name)//' % FULL_ROW_SEND_SIZE'     ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % full_row_send_size   )
    call memory_deallo(memor,trim(my_comm_name)//' % FULL_ROW_SEND_PERM'     ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % full_row_send_perm   )
    call memory_deallo(memor,trim(my_comm_name)//' % FULL_ROW_RECV_NEIGHTS'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % full_row_recv_neights)
    call memory_deallo(memor,trim(my_comm_name)//' % FULL_ROW_RECV_SIZE'     ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % full_row_recv_size   )
    call memory_deallo(memor,trim(my_comm_name)//' % FULL_ROW_RECV_PERM'     ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % full_row_recv_perm   )

    call memory_deallo(memor,trim(my_comm_name)//' % BOUND_SCAL'             ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % bound_scal         )
    call memory_deallo(memor,trim(my_comm_name)//' % BOUND_MULTIPLICITY'     ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % bound_multiplicity )
    call memory_deallo(memor,trim(my_comm_name)//' % BOUND_OWNER_RANK'       ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % bound_owner_rank   )

    call memory_deallo(memor,trim(my_comm_name)//' % BFACE_SIZE'             ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % bface_size         )
    call memory_deallo(memor,trim(my_comm_name)//' % BFACE_PERM'             ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % bface_perm         )

    call memory_deallo(memor,trim(my_comm_name)//' % NODE_NUMBER_IN_OWNER'   ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % node_number_in_owner)

    if( associated(comm % bound_matrix) ) then
       do ii = 1,size(comm % bound_matrix)
          call memory_deallo(memor,trim(my_comm_name)//' % BOUND_MATRIX % JA'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % bound_matrix(ii) % ja)
          call memory_deallo(memor,trim(my_comm_name)//' % BOUND_MATRIX % NZDOM_ii'    ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % bound_matrix(ii) % nzdom_ii)
          comm % bound_matrix(ii) % nzdom = 0
       end do
       deallocate(comm % bound_matrix)
    end if
    if( associated(comm % bound_mat_halo_send) ) then
       do ii = 1,size(comm % bound_mat_halo_send)
          call memory_deallo(memor,trim(my_comm_name)//' % BOUND_MAT_HALO_SEND % JA'      ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % bound_mat_halo_send(ii) % ja)
          call memory_deallo(memor,trim(my_comm_name)//' % BOUND_MAT_HALO_SEND % NZDOM_II','PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % bound_mat_halo_send(ii) % nzdom_ii)
          comm % bound_mat_halo_send(ii) % nzdom = 0
       end do
       deallocate(comm % bound_mat_halo_send)
    end if
    if( associated(comm % bound_mat_halo_recv) ) then
       do ii = 1,size(comm % bound_mat_halo_recv)
          call memory_deallo(memor,trim(my_comm_name)//' % BOUND_MAT_HALO_RECV % JA'       ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % bound_mat_halo_recv(ii) % ja)
          call memory_deallo(memor,trim(my_comm_name)//' % BOUND_MAT_HALO_RECV % NZDOM_II' ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % bound_mat_halo_recv(ii) % nzdom_ii)
          comm % bound_mat_halo_recv(ii) % nzdom = 0
       end do
       deallocate(comm % bound_mat_halo_recv)
    end if

    call memory_deallo(memor,trim(my_comm_name)//' % BEDGE_SIZE'             ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % bedge_size           )
    call memory_deallo(memor,trim(my_comm_name)//' % BEDGE_PERM'             ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % bedge_perm           )
    call memory_deallo(memor,trim(my_comm_name)//' % BEDGE_ADJA'             ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % bedge_adja           )
    call memory_deallo(memor,trim(my_comm_name)//' % BEDGE_SCAL'             ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % bedge_scal           )
    call memory_deallo(memor,trim(my_comm_name)//' % BEDGE_MULTIPLICITY'     ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % bedge_multiplicity   )
    call memory_deallo(memor,trim(my_comm_name)//' % BEDGE_OWNER_RANK'       ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % bedge_owner_rank     )
    call memory_deallo(memor,trim(my_comm_name)//' % LSCAT_PERM'             ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % lscat_perm           )
    call memory_deallo(memor,trim(my_comm_name)//' % MATRIX_IA'              ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % matrix_ia            )
    call memory_deallo(memor,trim(my_comm_name)//' % MATRIX_JA'              ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % matrix_ja            )
    call memory_deallo(memor,trim(my_comm_name)//' % MATRIX_AA'              ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % matrix_aa            )

    call memory_deallo(memor,trim(my_comm_name)//' % GHOST_SEND_NODE_SIZE'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % ghost_send_node_size   )
    call memory_deallo(memor,trim(my_comm_name)//' % GHOST_SEND_NODE_PERM'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % ghost_send_node_perm   )
    call memory_deallo(memor,trim(my_comm_name)//' % GHOST_RECV_NODE_SIZE'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % ghost_recv_node_size   )
    call memory_deallo(memor,trim(my_comm_name)//' % GHOST_RECV_NODE_PERM'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % ghost_recv_node_perm   )
    call memory_deallo(memor,trim(my_comm_name)//' % GHOST_SEND_ELEM_SIZE'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % ghost_send_elem_size   )
    call memory_deallo(memor,trim(my_comm_name)//' % GHOST_SEND_ELEM_PERM'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % ghost_send_elem_perm   )
    call memory_deallo(memor,trim(my_comm_name)//' % GHOST_RECV_ELEM_SIZE'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % ghost_recv_elem_size   )
    call memory_deallo(memor,trim(my_comm_name)//' % GHOST_RECV_ELEM_PERM'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % ghost_recv_elem_perm   )
    call memory_deallo(memor,trim(my_comm_name)//' % GHOST_SEND_BOUN_SIZE'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % ghost_send_boun_size   )
    call memory_deallo(memor,trim(my_comm_name)//' % GHOST_SEND_BOUN_PERM'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % ghost_send_boun_perm   )
    call memory_deallo(memor,trim(my_comm_name)//' % GHOST_RECV_BOUN_SIZE'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % ghost_recv_boun_size   )
    call memory_deallo(memor,trim(my_comm_name)//' % GHOST_RECV_BOUN_PERM'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',comm % ghost_recv_boun_perm   )
  
    
    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor
    if( allocated(my_comm_name) ) deallocate(my_comm_name)

  end subroutine deallo_full

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-28
  !> @brief   Allocate
  !> @details Allocate
  !> 
  !-----------------------------------------------------------------------

  subroutine alloca_basic(comm,MEMORY_COUNTER,comm_NAME)

    class(comm_data_par_basic),           intent(inout) :: comm
    integer(8),                 optional, intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    character(*),               optional, intent(in)    :: comm_NAME
    integer(8)                                          :: memor(2)
    character(len=:), allocatable                       :: my_comm_name

    my_comm_name = optional_argument('comm',comm_NAME)
    memor        = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
        
    call memory_alloca(memor,trim(my_comm_name)//' % NEIGHTS'    ,vacal,comm % neights    , comm % nneig    )
    call memory_alloca(memor,trim(my_comm_name)//' % BOUND_SIZE' ,vacal,comm % bound_size , comm % nneig+1  )
    call memory_alloca(memor,trim(my_comm_name)//' % BOUND_PERM' ,vacal,comm % bound_perm , comm % bound_dim)
    call memory_alloca(memor,trim(my_comm_name)//' % LSEND_SIZE' ,vacal,comm % lsend_size , comm % nneig+1  )
    call memory_alloca(memor,trim(my_comm_name)//' % LSEND_PERM' ,vacal,comm % lsend_perm , comm % lsend_dim)
    call memory_alloca(memor,trim(my_comm_name)//' % LRECV_SIZE' ,vacal,comm % lrecv_size , comm % nneig+1  )
    call memory_alloca(memor,trim(my_comm_name)//' % LRECV_PERM' ,vacal,comm % lrecv_perm , comm % lrecv_dim)
     
    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor
    if( allocated(my_comm_name) ) deallocate(my_comm_name)

  end subroutine alloca_basic

  subroutine alloca_full(comm,MEMORY_COUNTER,comm_NAME)
    
    class(comm_data_par),          intent(inout) :: comm       !< Input communicator
    integer(8),          optional, intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    character(len=*),    optional, intent(in)    :: COMM_NAME
    integer(8)                                   :: memor(2)
    character(len=:), allocatable                :: my_comm_name

    my_comm_name = optional_argument('comm',comm_NAME)
    memor        = optional_argument((/0_8,0_8/),MEMORY_COUNTER)

    call alloca_basic(comm,MEMORY_COUNTER,comm_NAME)
    
    call memory_alloca(memor,trim(my_comm_name)//' % BOUND_INVP',vacal,comm % bound_invp,comm % bound_dim)

    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor
    if( allocated(my_comm_name) ) deallocate(my_comm_name)

  end subroutine alloca_full
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-28
  !> @brief   Copy
  !> @details Copy COMM = COMM2
  !> 
  !-----------------------------------------------------------------------

  subroutine copy_basic(comm,comm2,MEMORY_COUNTER,COMM_NAME,COPY_NAME)

    class(comm_data_par_basic),       intent(inout) :: comm
    class(*),                         intent(inout) :: comm2
    integer(8),             optional, intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    character(*),           optional, intent(in)    :: COMM_NAME
    character(*),           optional, intent(in)    :: COPY_NAME
    character(len=:),       allocatable             :: my_comm_name
    character(len=:),       allocatable             :: my_copy_name
    integer(8)                                      :: memor(2)
   
    select type ( comm2 )
    class is ( comm_data_par_basic )
       
       memor           = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
       my_comm_name    = optional_argument(trim(comm2 % name),COMM_NAME)
       my_copy_name    = optional_argument(trim(comm  % name),COPY_NAME)
       
       call comm % deallo(MEMORY_COUNTER=memor)
       call comm % init  ()

       comm % nneig          = comm2 % nneig
       comm % bound_dim      = comm2 % bound_dim    
       comm % lsend_dim      = comm2 % lsend_dim    
       comm % lrecv_dim      = comm2 % lrecv_dim    
       comm % offset_npoin   = comm2 % offset_npoin   
       comm % offset_nelem   = comm2 % offset_nelem   
       comm % PAR_COMM_WORLD = comm2 % PAR_COMM_WORLD 
       comm % SIZE4          = comm2 % SIZE4          
       comm % RANK4          = comm2 % RANK4          
       comm % npoin          = comm2 % npoin          
       comm % npoi1          = comm2 % npoi1          
       comm % npoi2          = comm2 % npoi2          
       comm % npoi3          = comm2 % npoi3          

       call comm % alloca(MEMORY_COUNTER)

       if( associated(comm2 % neights)    ) comm % neights    = comm2 % neights
       if( associated(comm2 % bound_size) ) comm % bound_size = comm2 % bound_size
       if( associated(comm2 % bound_perm) ) comm % bound_perm = comm2 % bound_perm
       if( associated(comm2 % lsend_size) ) comm % lsend_size = comm2 % lsend_size
       if( associated(comm2 % lsend_perm) ) comm % lsend_perm = comm2 % lsend_perm
       if( associated(comm2 % lrecv_size) ) comm % lrecv_size = comm2 % lrecv_size
       if( associated(comm2 % lrecv_perm) ) comm % lrecv_perm = comm2 % lrecv_perm

    end select

  end subroutine copy_basic

  subroutine copy_full(comm,comm2,MEMORY_COUNTER,COMM_NAME,COPY_NAME)

    class(comm_data_par),             intent(inout) :: comm
    class(*),                         intent(inout) :: comm2
    integer(8),             optional, intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    character(*),           optional, intent(in)    :: COMM_NAME
    character(*),           optional, intent(in)    :: COPY_NAME
    character(len=:),       allocatable             :: my_comm_name
    character(len=:),       allocatable             :: my_copy_name
    integer(ip)                                     :: ii,isize
    integer(8)                                      :: memor(2)
    
    select type ( comm2 )
       
    class is ( comm_data_par_basic )

       call copy_basic(comm,comm2,MEMORY_COUNTER,COMM_NAME,COPY_NAME)
       
    class is ( comm_data_par )
       
       memor           = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
       my_comm_name    = optional_argument(trim(comm2 % name),COMM_NAME)
       my_copy_name    = optional_argument(trim(comm  % name),COPY_NAME)
       
       call comm % deallo(MEMORY_COUNTER=memor)
       call comm % init  ()

       comm % nneig                 =  comm2 % nneig
       comm % nneig_1               =  comm2 % nneig_1
       comm % nneig_2               =  comm2 % nneig_2
       comm % npoi1                 =  comm2 % npoi1
       comm % npoi2                 =  comm2 % npoi2
       comm % npoi3                 =  comm2 % npoi3
       comm % npoin                 =  comm2 % npoin
       comm % bound_dim             =  comm2 % bound_dim

       comm % nedg1                 =  comm2 % nedg1
       comm % nedg2                 =  comm2 % nedg2
       comm % nedg3                 =  comm2 % nedg3

       comm % offset_npoin          =  comm2 % offset_npoin
       comm % offset_nelem          =  comm2 % offset_nelem
       comm % offset_nboun          =  comm2 % offset_nboun

       comm % lsend_dim             =  comm2 % lsend_dim
       comm % lrecv_dim             =  comm2 % lrecv_dim
       comm % lscat_dim             =  comm2 % lscat_dim
       comm % matrix_nzdom          =  comm2 % matrix_nzdom
       comm % bface_dim             =  comm2 % bface_dim
       comm % ghost_send_elem_dim   =  comm2 % ghost_send_elem_dim
       comm % ghost_recv_elem_dim   =  comm2 % ghost_recv_elem_dim
       comm % ghost_send_node_dim   =  comm2 % ghost_send_node_dim
       comm % ghost_recv_node_dim   =  comm2 % ghost_recv_node_dim
       comm % ghost_send_boun_dim   =  comm2 % ghost_send_boun_dim
       comm % ghost_recv_boun_dim   =  comm2 % ghost_recv_boun_dim

       comm % PAR_COMM_WORLD        =  comm2 % PAR_COMM_WORLD
       comm % SIZE4                 =  comm2 % SIZE4
       comm % RANK4                 =  comm2 % RANK4

       comm % bedge_dim             =  comm2 % bedge_dim

       call memory_copy(memor,trim(my_comm_name)//' % NEIGHTS',         vacal,comm2 % neights         , comm % neights,         'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % NEIGHTS'   )
       call memory_copy(memor,trim(my_comm_name)//' % NEIGHTS_ORDERED', vacal,comm2 % neights_ordered , comm % neights_ordered, 'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % NEIGHTS_ORDERED' )
       call memory_copy(memor,trim(my_comm_name)//' % PERM_ORDERED',    vacal,comm2 % perm_ordered    , comm % perm_ordered,    'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % PERM_ORDERED'  )

       call memory_copy(memor,trim(my_comm_name)//' % BOUND_SIZE',      vacal,comm2 % bound_size      , comm % bound_size,      'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_SIZE'   )
       call memory_copy(memor,trim(my_comm_name)//' % BOUND_PERM',      vacal,comm2 % bound_perm      , comm % bound_perm,      'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_PERM'   )
       call memory_copy(memor,trim(my_comm_name)//' % BOUND_INVP',      vacal,comm2 % bound_invp      , comm % bound_invp,      'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_INVP'   )

       call memory_copy(memor,trim(my_comm_name)//' % BOUND_SCAL',          vacal,comm2 % bound_scal            , comm % bound_scal,          'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_SCAL'          )
       call memory_copy(memor,trim(my_comm_name)//' % BOUND_MULTIPLICITY',  vacal,comm2 % bound_multiplicity    , comm % bound_multiplicity,  'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_MULTIPLICITY'  )
       call memory_copy(memor,trim(my_comm_name)//' % BOUND_OWNER_RANK',    vacal,comm2 % bound_owner_rank      , comm % bound_owner_rank,    'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_OWNER_RANK'    )
       call memory_copy(memor,trim(my_comm_name)//' % NODE_NUMBER_IN_OWNER',vacal,comm2 % node_number_in_owner  , comm % node_number_in_owner,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % NODE_NUMBER_IN_OWNER')

       if( associated(comm2 % bound_matrix) ) then
          isize = size(comm2 % bound_matrix)
          allocate(comm % bound_matrix(isize))
          do ii = 1,isize
             call memory_copy(memor,trim(my_comm_name)//' % BOUND_MATRIX % JA',      vacal,&
                  & comm2 % bound_matrix(ii) % ja,             comm % bound_matrix(ii) % ja,              'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_MATRIX % JA')
             call memory_copy(memor,trim(my_comm_name)//' % BOUND_MATRIX % NZDOM_II',vacal,&
                  & comm2 % bound_matrix(ii) % nzdom_ii,       comm % bound_matrix(ii) % nzdom_ii,        'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_MATRIX % NZDOM_II')
             comm % bound_matrix(ii) % nzdom  =  comm2 % bound_matrix(ii) % nzdom
          end do
       end if

       if( associated(comm2 % bound_mat_halo_send) ) then
          isize = size(comm2 % bound_mat_halo_send)
          allocate(comm % bound_mat_halo_send(isize))
          do ii = 1,isize
             call memory_copy(memor,trim(my_comm_name)//' % BOUND_MAT_HALO_SEND % JA',vacal,&
                  & comm2 % bound_mat_halo_send(ii) % ja,             comm % bound_mat_halo_send(ii) % ja,              'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_MAT_HALO_SEND % JA')
             call memory_copy(memor,trim(my_comm_name)//' % BOUND_MAT_HALO_SEND % NZDOM_II',        vacal,&
                  & comm2 % bound_mat_halo_send(ii) % nzdom_ii,       comm % bound_mat_halo_send(ii) % nzdom_ii,        'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_MAT_HALO_SEND % NZDOM_II')
             comm % bound_mat_halo_send(ii) % nzdom  =  comm2 % bound_mat_halo_send(ii) % nzdom
          end do
       end if

       if( associated(comm2 % bound_mat_halo_recv) ) then
          isize = size(comm2 % bound_mat_halo_recv)
          allocate(comm % bound_mat_halo_recv(isize))
          do ii = 1,isize
             call memory_copy(memor,trim(my_comm_name)//' % BOUND_MAT_HALO_RECV % JA',vacal,&
                  & comm2 % bound_mat_halo_recv(ii) % ja,             comm % bound_mat_halo_recv(ii) % ja,              'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_MAT_HALO_RECV % JA')
             call memory_copy(memor,trim(my_comm_name)//' % BOUND_MAT_HALO_RECV % NZDOM_II',        vacal,&
                  & comm2 % bound_mat_halo_recv(ii) % nzdom_ii,       comm % bound_mat_halo_recv(ii) % nzdom_ii,        'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_MAT_HALO_RECV % NZDOM_II')
             comm % bound_mat_halo_recv(ii) % nzdom  =  comm2 % bound_mat_halo_recv(ii) % nzdom
          end do
       end if

       call memory_copy(memor,trim(my_comm_name)//' % BEDGE_SIZE',          vacal,comm2 % bedge_size            , comm % bedge_size,          'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BEDGE_SIZE'          )
       call memory_copy(memor,trim(my_comm_name)//' % BEDGE_PERM',          vacal,comm2 % bedge_perm            , comm % bedge_perm,          'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BEDGE_PERM'          )
       call memory_copy(memor,trim(my_comm_name)//' % BEDGE_ADJA',          vacal,comm2 % bedge_adja            , comm % bedge_adja,          'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BEDGE_ADJA'          )
       call memory_copy(memor,trim(my_comm_name)//' % BEDGE_SCAL',          vacal,comm2 % bedge_scal            , comm % bedge_scal,          'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BEDGE_SCAL'          )
       call memory_copy(memor,trim(my_comm_name)//' % BEDGE_MULTIPLICITY',  vacal,comm2 % bedge_multiplicity    , comm % bedge_multiplicity,  'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BEDGE_MULTIPLICITY'  )
       call memory_copy(memor,trim(my_comm_name)//' % BEDGE_OWNER_RANK',    vacal,comm2 % bedge_owner_rank      , comm % bedge_owner_rank,    'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BEDGE_OWNER_RANK'    )

       call memory_copy(memor,trim(my_comm_name)//' % LSEND_SIZE',          vacal,comm2 % lsend_size            , comm % lsend_size,          'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % LSEND_SIZE'          )
       call memory_copy(memor,trim(my_comm_name)//' % LRECV_SIZE',          vacal,comm2 % lrecv_size            , comm % lrecv_size,          'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % LRECV_SIZE'          )
       call memory_copy(memor,trim(my_comm_name)//' % LSEND_PERM',          vacal,comm2 % lsend_perm            , comm % lsend_perm,          'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % LSEND_PERM'          )
       call memory_copy(memor,trim(my_comm_name)//' % LRECV_PERM',          vacal,comm2 % lrecv_perm            , comm % lrecv_perm,          'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % LRECV_PERM'          )
       call memory_copy(memor,trim(my_comm_name)//' % LSCAT_PERM',          vacal,comm2 % lscat_perm            , comm % lscat_perm,          'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % LSCAT_PERM'          )

       call memory_copy(memor,trim(my_comm_name)//' % MATRIX_IA',           vacal,comm2 % matrix_ia             , comm % matrix_ia,           'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % MATRIX_IA'           )
       call memory_copy(memor,trim(my_comm_name)//' % MATRIX_JA',           vacal,comm2 % matrix_ja             , comm % matrix_ja,           'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % MATRIX_JA'           )
       call memory_copy(memor,trim(my_comm_name)//' % MATRIX_AA',           vacal,comm2 % matrix_aa             , comm % matrix_aa,           'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % MATRIX_AA'           )

       call memory_copy(memor,trim(my_comm_name)//' % FULL_ROW_SEND_NEIGHTS',vacal,comm2 % full_row_send_neights,comm % full_row_send_neights,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % FULL_ROW_SEND_NEIGHTS')
       call memory_copy(memor,trim(my_comm_name)//' % FULL_ROW_SEND_SIZE'   ,vacal,comm2 % full_row_send_size   ,comm % full_row_send_size   ,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % FULL_ROW_SEND_SIZE'   )
       call memory_copy(memor,trim(my_comm_name)//' % FULL_ROW_SEND_PERM'   ,vacal,comm2 % full_row_send_perm   ,comm % full_row_send_perm   ,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % FULL_ROW_SEND_PERM'   )
       call memory_copy(memor,trim(my_comm_name)//' % FULL_ROW_RECV_NEIGHTS',vacal,comm2 % full_row_recv_neights,comm % full_row_recv_neights,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % FULL_ROW_RECV_NEIGHTS')
       call memory_copy(memor,trim(my_comm_name)//' % FULL_ROW_RECV_SIZE'   ,vacal,comm2 % full_row_recv_size   ,comm % full_row_recv_size   ,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % FULL_ROW_RECV_SIZE'   )
       call memory_copy(memor,trim(my_comm_name)//' % FULL_ROW_RECV_PERM'   ,vacal,comm2 % full_row_recv_perm   ,comm % full_row_recv_perm   ,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % FULL_ROW_RECV_PERM'   )

       call memory_copy(memor,trim(my_comm_name)//' % BFACE_SIZE',          vacal,comm2 % bface_size            , comm % bface_size,          'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BFACE_SIZE'          )
       call memory_copy(memor,trim(my_comm_name)//' % BFACE_PERM',          vacal,comm2 % bface_perm            , comm % bface_perm,          'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BFACE_PERM'          )

       call memory_copy(memor,trim(my_comm_name)//' % GHOST_SEND_ELEM_SIZE',vacal,comm2 % ghost_send_elem_size  , comm % ghost_send_elem_size,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % GHOST_SEND_ELEM_SIZE')
       call memory_copy(memor,trim(my_comm_name)//' % GHOST_SEND_ELEM_PERM',vacal,comm2 % ghost_send_elem_perm  , comm % ghost_send_elem_perm,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % GHOST_SEND_ELEM_PERM')
       call memory_copy(memor,trim(my_comm_name)//' % GHOST_RECV_ELEM_SIZE',vacal,comm2 % ghost_recv_elem_size  , comm % ghost_recv_elem_size,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % GHOST_RECV_ELEM_SIZE')
       call memory_copy(memor,trim(my_comm_name)//' % GHOST_RECV_ELEM_PERM',vacal,comm2 % ghost_recv_elem_perm  , comm % ghost_recv_elem_perm,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % GHOST_RECV_ELEM_PERM')
       call memory_copy(memor,trim(my_comm_name)//' % GHOST_SEND_NODE_SIZE',vacal,comm2 % ghost_send_node_size  , comm % ghost_send_node_size,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % GHOST_SEND_NODE_SIZE')
       call memory_copy(memor,trim(my_comm_name)//' % GHOST_SEND_NODE_PERM',vacal,comm2 % ghost_send_node_perm  , comm % ghost_send_node_perm,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % GHOST_SEND_NODE_PERM')
       call memory_copy(memor,trim(my_comm_name)//' % GHOST_RECV_NODE_SIZE',vacal,comm2 % ghost_recv_node_size  , comm % ghost_recv_node_size,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % GHOST_RECV_NODE_SIZE')
       call memory_copy(memor,trim(my_comm_name)//' % GHOST_RECV_NODE_PERM',vacal,comm2 % ghost_recv_node_perm  , comm % ghost_recv_node_perm,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % GHOST_RECV_NODE_PERM')

       call memory_copy(memor,trim(my_comm_name)//' % GHOST_SEND_BOUN_SIZE',vacal,comm2 % ghost_send_boun_size  , comm % ghost_send_boun_size,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % GHOST_SEND_BOUN_SIZE')
       call memory_copy(memor,trim(my_comm_name)//' % GHOST_SEND_BOUN_PERM',vacal,comm2 % ghost_send_boun_perm  , comm % ghost_send_boun_perm,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % GHOST_SEND_BOUN_PERM')
       call memory_copy(memor,trim(my_comm_name)//' % GHOST_RECV_BOUN_SIZE',vacal,comm2 % ghost_recv_boun_size  , comm % ghost_recv_boun_size,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % GHOST_RECV_BOUN_SIZE')
       call memory_copy(memor,trim(my_comm_name)//' % GHOST_RECV_BOUN_PERM',vacal,comm2 % ghost_recv_boun_perm  , comm % ghost_recv_boun_perm,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % GHOST_RECV_BOUN_PERM')

       if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor
       if( allocated(my_comm_name) ) deallocate(my_comm_name)
       if( allocated(my_copy_name) ) deallocate(my_copy_name)
       
    end select

  end subroutine copy_full

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-28
  !> @brief   Copy
  !> @details Copy COMM = COMM2
  !> 
  !-----------------------------------------------------------------------

  subroutine move(comm,comm2,MEMORY_COUNTER)

    class(comm_data_par_basic),           intent(inout) :: comm
    class(comm_data_par_basic),           intent(inout) :: comm2
    integer(8),                 optional, intent(inout) :: MEMORY_COUNTER(2) !< Memory counters

    call comm  % copy  (comm2,MEMORY_COUNTER)
    call comm2 % deallo(MEMORY_COUNTER)
    
  end subroutine move

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-28
  !> @brief   Merge and/or compact a communicator
  !> @details Merge and compact COMM = COMM U COMM2
  !> 
  !-----------------------------------------------------------------------

  subroutine merge(comm,comm2,MEMORY_COUNTER)

    class(comm_data_par_basic),           intent(inout) :: comm
    class(comm_data_par_basic), optional, intent(in)    :: comm2
    integer(8),                 optional, intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(8)                                          :: memor_loc(2)
    integer(ip)                                         :: ineig,ii,kk,dom
    integer(ip)                                         :: nneig1,kneig,jneig
    logical(lg)                                         :: ifoun
    logical(lg)                                         :: do_something
    type(comm_data_par_basic)                           :: comm1
    integer(ip),                pointer                 :: bound_size(:)
    integer(ip),                pointer                 :: lsend_size(:)
    integer(ip),                pointer                 :: lrecv_size(:)
    integer(ip),                pointer                 :: neights(:)
    integer(ip),                pointer                 :: perm_cpy(:)
    type(i1p),                  pointer                 :: bound_perm(:)
    type(i1p),                  pointer                 :: lsend_perm(:)
    type(i1p),                  pointer                 :: lrecv_perm(:)
    integer(ip),                pointer                 :: p(:)
    !
    ! Something to do?
    !
    do_something = .false.
    if( comm % nneig > 0 ) do_something = .true.
    if( present(comm2) ) then
       if( comm2 % nneig > 0 ) do_something = .true.
    end if
    ! 
    ! Some checks
    !
#ifndef MPI_OFF
    if( present(comm2) ) then
       if( comm % SIZE4          /= comm2 % SIZE4 )          call runend('DEF_KINTYP_COMM: SIZES OF COMMUNICATORS SHOULD BE EQUAL')
       if( comm % RANK4          /= comm2 % RANK4 )          call runend('DEF_KINTYP_COMM: RANKS OF COMMUNICATORS SHOULD BE EQUAL')       
       !if( comm % PAR_COMM_WORLD /= comm2 % PAR_COMM_WORLD ) call runend('DEF_KINTYP_COMM: COMMUNICATORS SHOULD BE EQUAL')       
    end if
#endif
    
    if( do_something ) then
       
       memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)

       nullify(bound_size)
       nullify(lsend_size)
       nullify(lrecv_size)
       nullify(neights)
       nullify(bound_perm)
       nullify(lsend_perm)
       nullify(lrecv_perm)
       nullify(perm_cpy)

       call comm1 % init()
       call comm1 % move(comm,MEMORY_COUNTER=memor_loc)
       call comm  % init()

       comm % offset_npoin   = comm1 % offset_npoin    
       comm % offset_nelem   = comm1 % offset_nelem   
       comm % PAR_COMM_WORLD = comm1 % PAR_COMM_WORLD 
       comm % SIZE4          = comm1 % SIZE4          
       comm % RANK4          = comm1 % RANK4
       !
       ! First communicator COMM1
       !
       call memory_alloca(memor_loc,'NEIGHTS'   ,vacal,neights,   comm1 % nneig)
       call memory_alloca(memor_loc,'BOUND_SIZE',vacal,bound_size,comm1 % nneig)
       call memory_alloca(memor_loc,'LSEND_SIZE',vacal,lsend_size,comm1 % nneig)
       call memory_alloca(memor_loc,'LRECV_SIZE',vacal,lrecv_size,comm1 % nneig)
       call memory_alloca(memor_loc,'BOUND_PERM',vacal,bound_perm,comm1 % nneig)
       call memory_alloca(memor_loc,'LSEND_PERM',vacal,lsend_perm,comm1 % nneig)
       call memory_alloca(memor_loc,'LRECV_PERM',vacal,lrecv_perm,comm1 % nneig)
       
       jneig = 0
       if( comm1 % nneig > 0 ) then
          do ineig = 1,comm1 % nneig
             ifoun = .false.
             if( associated(comm1 % bound_size) ) then
                if( comm1 % bound_size(ineig+1)-comm1 % bound_size(ineig) > 0 ) ifoun= .true.
             end if
             if( associated(comm1 % lsend_size) ) then
                if( comm1 % lsend_size(ineig+1)-comm1 % lsend_size(ineig) > 0 ) ifoun= .true.
             end if
             if( associated(comm1 % lrecv_size) ) then
                if( comm1 % lrecv_size(ineig+1)-comm1 % lrecv_size(ineig) > 0 ) ifoun= .true.
             end if
             if( ifoun ) then
                jneig          = jneig + 1
                neights(jneig) = comm1 % neights(ineig)
                if( associated(comm1 % bound_size) ) then
                   bound_size(jneig) = comm1 % bound_size(ineig+1)-comm1 % bound_size(ineig)
                   if( bound_size(jneig) > 0 ) then
                      call memory_alloca(memor_loc,'BOUND_PERM % L',vacal,bound_perm(jneig) % l,bound_size(jneig))
                      bound_perm(jneig) % l = comm1 % bound_perm(comm1 % bound_size(ineig):comm1 % bound_size(ineig+1)-1)
                   end if
                end if
                if( associated(comm1 % lsend_size) ) then
                   lsend_size(jneig) = comm1 % lsend_size(ineig+1)-comm1 % lsend_size(ineig)
                   if( lsend_size(jneig) > 0 ) then
                      call memory_alloca(memor_loc,'LSEND_PERM % L',vacal,lsend_perm(jneig) % l,lsend_size(jneig))
                      lsend_perm(jneig) % l = comm1 % lsend_perm(comm1 % lsend_size(ineig):comm1 % lsend_size(ineig+1)-1)
                   end if
                end if
                if( associated(comm1 % lrecv_size) ) then
                   lrecv_size(jneig) = comm1 % lrecv_size(ineig+1)-comm1 % lrecv_size(ineig)
                   if( lrecv_size(jneig) > 0 ) then
                      call memory_alloca(memor_loc,'LRECV_PERM % L',vacal,lrecv_perm(jneig) % l,lrecv_size(jneig))
                      lrecv_perm(jneig) % l = comm1 % lrecv_perm(comm1 % lrecv_size(ineig):comm1 % lrecv_size(ineig+1)-1)
                   end if
                end if
             end if
          end do
       end if
       !
       ! Check second communicator COMM2
       !
       if( present(comm2) ) then
          nneig1 = jneig 
          do ineig = 1,comm2 % nneig
             dom  = comm2 % neights(ineig)
             ii   = comm_findloc(neights(1:nneig1),VAL=dom)
             if( ii == 0 ) then
                ifoun = .false.
                if( associated(comm2 % bound_size) ) then
                   if( comm2 % bound_size(ineig+1)-comm2 % bound_size(ineig) > 0 ) ifoun= .true.
                end if
                if( associated(comm2 % lsend_size) ) then
                   if( comm2 % lsend_size(ineig+1)-comm2 % lsend_size(ineig) > 0 ) ifoun= .true.
                end if
                if( associated(comm2 % lrecv_size) ) then
                   if( comm2 % lrecv_size(ineig+1)-comm2 % lrecv_size(ineig) > 0 ) ifoun= .true.
                end if
                if( ifoun ) then
                   jneig          = jneig + 1
                   neights(jneig) = comm2 % neights(ineig)
                   if( associated(comm2 % bound_size) ) then
                      bound_size(jneig) = comm2 % bound_size(ineig+1)-comm2 % bound_size(ineig)                   
                      if( bound_size(jneig) > 0 ) then
                         call memory_alloca(memor_loc,'BOUND_PERM % L',vacal,bound_perm(jneig) % l,bound_size(jneig))
                         bound_perm(jneig) % l = comm2 % bound_perm(comm2 % bound_size(ineig):comm2 % bound_size(ineig+1)-1)
                         ifoun = .true.
                      end if
                   end if
                   if( associated(comm2 % lsend_size) ) then
                      lsend_size(jneig) = comm2 % lsend_size(ineig+1)-comm2 % lsend_size(ineig)
                      if( lsend_size(jneig) > 0 ) then
                         call memory_alloca(memor_loc,'LSEND_PERM % L',vacal,lsend_perm(jneig) % l,lsend_size(jneig))
                         lsend_perm(jneig) % l = comm2 % lsend_perm(comm2 % lsend_size(ineig):comm2 % lsend_size(ineig+1)-1)
                         ifoun = .true.
                      end if
                   end if
                   if( associated(comm2 % lrecv_size) ) then
                      lrecv_size(jneig) = comm2 % lrecv_size(ineig+1)-comm2 % lrecv_size(ineig)
                      if( lrecv_size(jneig) > 0 ) then
                         call memory_alloca(memor_loc,'LRECV_PERM % L',vacal,lrecv_perm(jneig) % l,lrecv_size(jneig))
                         lrecv_perm(jneig) % l = comm2 % lrecv_perm(comm2 % lrecv_size(ineig):comm2 % lrecv_size(ineig+1)-1)
                         ifoun = .true.
                      end if
                   end if
                end if
             else
                kneig          = ii 
                  if( associated(comm2 % bound_size) ) then
                     kk                = comm2 % bound_size(ineig+1)-comm2 % bound_size(ineig)
                     bound_size(kneig) = bound_size(kneig) + kk
                     if( kk > 0 ) then
                        p => comm2 % bound_perm(comm2 % bound_size(ineig):comm2 % bound_size(ineig+1)-1)
                        call memory_append(memor_loc,'BOUND_PERM % L',vacal,bound_perm(kneig) % l,p)
                     end if
                  end if
                  if( associated(comm2 % lsend_size) ) then
                     kk                = comm2 % lsend_size(ineig+1)-comm2 % lsend_size(ineig)
                     lsend_size(kneig) = lsend_size(kneig) + kk
                     if( kk > 0 ) then
                        p => comm2 % lsend_perm(comm2 % lsend_size(ineig):comm2 % lsend_size(ineig+1)-1)
                        call memory_append(memor_loc,'LSEND_PERM % L',vacal,lsend_perm(kneig) % l,p)
                     end if
                  end if
                  if( associated(comm2 % lrecv_size) ) then
                     kk                = comm2 % lrecv_size(ineig+1)-comm2 % lrecv_size(ineig)
                     lrecv_size(kneig) = lrecv_size(kneig) + kk
                     if( kk > 0 ) then
                        p => comm2 % lrecv_perm(comm2 % lrecv_size(ineig):comm2 % lrecv_size(ineig+1)-1)
                        call memory_append(memor_loc,'LRECV_PERM % L',vacal,lrecv_perm(kneig) % l,p)
                     end if
                  end if
               end if
            end do
         end if
!!$         block
!!$           use def_master, only : kfl_paral
!!$           if( kfl_paral == 3 ) then
!!$              do ineig = 1,comm1 % nneig
!!$                 if( comm1 % lsend_dim > 0 ) print*,'a1=',comm1 % neights(ineig),': ',comm1 % lsend_size(ineig+1)-comm1 % lsend_size(ineig)
!!$                 if( comm1 % lrecv_dim > 0 ) print*,'b1=',comm1 % neights(ineig),': ',comm1 % lrecv_size(ineig+1)-comm1 % lrecv_size(ineig)
!!$              end do
!!$              do ineig = 1,comm2 % nneig
!!$                 if( comm2 % lsend_dim > 0 ) print*,'a2=',comm2 % neights(ineig),': ',comm2 % lsend_size(ineig+1)-comm2 % lsend_size(ineig)
!!$                 if( comm2 % lrecv_dim > 0 ) print*,'b2=',comm2 % neights(ineig),': ',comm2 % lrecv_size(ineig+1)-comm2 % lrecv_size(ineig)
!!$              end do
!!$           end if
!!$         end block
         !
         ! Merge
         !
         comm % nneig = jneig
         do ineig = 1,comm % nneig
            comm % bound_dim      = comm % bound_dim + bound_size(ineig)
            comm % lsend_dim      = comm % lsend_dim + lsend_size(ineig)
            comm % lrecv_dim      = comm % lrecv_dim + lrecv_size(ineig)
         end do
         call comm % alloca(MEMORY_COUNTER=memor_loc)          
         do ineig = 1,comm % nneig
            comm % bound_size(ineig) = bound_size(ineig)
            comm % lsend_size(ineig) = lsend_size(ineig)
            comm % lrecv_size(ineig) = lrecv_size(ineig)
            comm % neights(ineig)    = neights   (ineig)
         end do
         if(comm % bound_dim > 0 ) call number_to_linked_list(comm % nneig,comm % bound_size)
         if(comm % lsend_dim > 0 ) call number_to_linked_list(comm % nneig,comm % lsend_size)
         if(comm % lrecv_dim > 0 ) call number_to_linked_list(comm % nneig,comm % lrecv_size)
         do ineig = 1,comm % nneig
            kk = 0
            do ii = comm % bound_size(ineig),comm % bound_size(ineig+1)-1
               kk = kk + 1
               comm % bound_perm(ii) = bound_perm(ineig) % l(kk)
            end do
            kk = 0
            do ii = comm % lsend_size(ineig),comm % lsend_size(ineig+1)-1
               kk = kk + 1
               comm % lsend_perm(ii) = lsend_perm(ineig) % l(kk)
            end do
            kk = 0
            do ii = comm % lrecv_size(ineig),comm % lrecv_size(ineig+1)-1
               kk = kk + 1
               comm % lrecv_perm(ii) = lrecv_perm(ineig) % l(kk)
            end do
         end do
!!$         block
!!$           use def_master, only : kfl_paral
!!$           if( kfl_paral == 3 ) then
!!$              do ineig = 1,comm % nneig
!!$                 if( comm % lsend_dim > 0 ) print*,'a3=',comm % neights(ineig),': ',comm % lsend_size(ineig+1)-comm % lsend_size(ineig)
!!$                 if( comm % lrecv_dim > 0 ) print*,'b3=',comm % neights(ineig),': ',comm % lrecv_size(ineig+1)-comm % lrecv_size(ineig)
!!$              end do
!!$           end if
!!$         end block
         !
         ! Deallocate
         !
         call memory_deallo(memor_loc,'NEIGHTS'   ,vacal,neights   )
         call memory_deallo(memor_loc,'BOUND_SIZE',vacal,bound_size)
         call memory_deallo(memor_loc,'LSEND_SIZE',vacal,lsend_size)
         call memory_deallo(memor_loc,'LRECV_SIZE',vacal,lrecv_size)
         call memory_deallo(memor_loc,'BOUND_PERM',vacal,bound_perm)
         call memory_deallo(memor_loc,'LSEND_PERM',vacal,lsend_perm)
         call memory_deallo(memor_loc,'LRECV_PERM',vacal,lrecv_perm)

         call comm1 % deallo(MEMORY_COUNTER=memor_loc)

         if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc

      end if

    end subroutine merge

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-05-19
  !> @brief   Number to linked list
  !> @details Number to linked list
  !> 
  !-----------------------------------------------------------------------

  pure subroutine number_to_linked_list(nn,ia)
    
    integer(ip), intent(in)    :: nn
    integer(ip), intent(inout) :: ia(nn+1)
    integer(ip)                :: ii,kk,ll

    kk    = ia(1)
    ia(1) = 1 
    do ii = 2,nn+1
       ll     = ia(ii)
       ia(ii) = ia(ii-1) + kk
       kk     = ll
    end do

  end subroutine number_to_linked_list
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-05-19
  !> @brief   Find location
  !> @details Find location. findloc not supported by old versions of gcc
  !> 
  !-----------------------------------------------------------------------

  pure function comm_findloc(xx,val) result(pos)

    integer(ip),          intent(in) :: xx(:)
    integer(ip),          intent(in) :: val
    integer(ip)                      :: kk,pos

    pos = 0
    if( size(xx) > 0 ) then
       kk = 0
       do while( kk < size(xx,KIND=ip) )
          kk = kk + 1
          if( xx(kk) == val ) then
             pos = kk
             return
          end if
       end do
    end if
 
  end function comm_findloc

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-28
  !> @brief   Permute
  !> @details Permute
  !> 
  !-----------------------------------------------------------------------

  subroutine permute(comm,permu)

    class(comm_data_par_basic),           intent(inout) :: comm
    integer(ip),                 pointer, intent(in)    :: permu(:)

    if( .not. associated(permu) ) then
       call runend('PERMUTE: PERMU NOT ASSOCIATED')
    end if
    
    if( associated(comm % bound_perm) ) comm % bound_perm = permu(comm % bound_perm)
    if( associated(comm % lsend_perm) ) comm % lsend_perm = permu(comm % lsend_perm)
    if( associated(comm % lrecv_perm) ) comm % lrecv_perm = permu(comm % lrecv_perm)

  end subroutine permute

    !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    06/03/2014
  !> @brief   Create a sub communication array
  !> @details Given a communicaiton array structure, create a sub
  !           communication structure for a subset of nodes
  !
  !----------------------------------------------------------------------

  subroutine sub(comm,comm2,mask,MEMORY_COUNTER,COMM_NAME)
    
    class(comm_data_par),                    intent(inout) :: comm           !< Output communicator
    class(comm_data_par),                    intent(in)    :: comm2          !< Input communicator
    integer(ip),                   pointer,  intent(in)    :: mask(:)        !< mask (=1 to consider node)
    integer(8),          optional,           intent(inout) :: MEMORY_COUNTER(2)  !< Memory counter
    character(*),        optional,           intent(in)    :: COMM_NAME
    integer(ip),                   pointer                 :: bound_perm(:)
    integer(ip),                   pointer                 :: bound_invp(:)
    integer(ip),                   pointer                 :: bound_size(:)
    integer(ip),                   pointer                 :: neights(:)
    integer(ip)                                            :: nneig,bound_dim
    integer(ip)                                            :: nneig_1,nneig_2
    integer(ip)                                            :: ineig,jj,jneig
    integer(ip)                                            :: ipoin
    character(50)                                          :: my_comm_name
    integer(8)                                             :: memor(2)

    memor                 = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
    my_comm_name          = optional_argument('COMM',COMM_NAME)
    
    comm % PAR_COMM_WORLD = comm2 % PAR_COMM_WORLD
    comm % SIZE4          = comm2 % SIZE4
    comm % RANK4          = comm2 % RANK4
    
    if( comm2 % bound_dim > 0 ) then

       nullify(bound_perm)
       nullify(bound_invp)
       nullify(bound_size)
       nullify(neights)

       allocate( bound_perm(comm2 % bound_dim) )
       allocate( bound_invp(comm2 % bound_dim) )
       allocate( bound_size(comm2 % nneig)     )
       allocate( neights   (comm2 % nneig)     )

       nneig     = 0
       bound_dim = 0

       do jj = 1,comm2 % bound_dim
          bound_perm(jj) = comm2 % bound_perm(jj)
          bound_invp(jj) = comm2 % bound_invp(jj)
       end do
       do ineig = 1,comm2 % nneig
          bound_size(ineig) = 0
       end do

       do ineig = 1,comm2 % nneig
          do jj = comm2 % bound_size(ineig),comm2 % bound_size(ineig+1)-1
             ipoin = comm2 % bound_perm(jj)
             if( mask(ipoin) > 0 ) then
                bound_size(ineig) = bound_size(ineig) + 1
             else
                bound_perm(jj) = 0
                bound_invp(jj) = 0
             end if
          end do
          bound_dim = bound_dim + bound_size(ineig)

          if( bound_size(ineig) > 0 ) nneig = nneig + 1
       end do
       !
       ! Allocate interzone communication array
       !
       comm % nneig     = nneig
       comm % bound_dim = bound_dim
       comm % npoi1     = comm2 % npoi1
       comm % npoi2     = comm2 % npoi2
       comm % npoi3     = comm2 % npoi3
       comm % npoin     = comm2 % npoin
       call memory_alloca(memor,trim(my_comm_name)//' % NEIGHTS'        ,vacal,comm % neights,nneig)
       call memory_alloca(memor,trim(my_comm_name)//' % NEIGHTS_ORDERED',vacal,comm % neights_ordered,nneig)
       call memory_alloca(memor,trim(my_comm_name)//' % PERM_ORDERED'   ,vacal,comm % perm_ordered,nneig)
       call memory_alloca(memor,trim(my_comm_name)//' % BOUND_PERM'     ,vacal,comm % bound_perm,bound_dim)
       call memory_alloca(memor,trim(my_comm_name)//' % BOUND_INVP'     ,vacal,comm % bound_invp,bound_dim)
       call memory_alloca(memor,trim(my_comm_name)//' % BOUND_SIZE'     ,vacal,comm % bound_size,nneig+1)

       !allocate( comm % neights(nneig) )
       !allocate( comm % bound_perm(bound_dim) )
       !allocate( comm % bound_size(nneig+1) )
       !
       ! Permutation array BOUND_PERM(1:BOUND_DIM)
       !
       jneig = 0
       bound_dim = 0
       do ineig = 1,comm2 % nneig
          if( bound_size(ineig) > 0 ) then

             jneig = jneig + 1
             comm % neights(jneig)         = comm2 % neights(ineig)
             comm % neights_ordered(jneig) = comm2 % neights_ordered(ineig)
             comm % perm_ordered(jneig)    = comm2 % perm_ordered(ineig)

             do jj = comm2 % bound_size(ineig),comm2 % bound_size(ineig+1)-1
                ipoin = bound_perm(jj)
                if( ipoin > 0 ) then
                   bound_dim = bound_dim + 1
                   comm % bound_perm(bound_dim) = ipoin
                   comm % bound_invp(bound_dim) = ipoin
                end if
             end do

          end if
       end do
       ineig_1: do ineig = 1,comm2 % nneig
          if( ineig > int(comm2 % RANK4,ip) ) then
             nneig_1 = ineig-1
             nneig_2 = ineig
             exit ineig_1
          end if
       end do ineig_1
       !
       ! Construct linked list BOUND_SIZE(1:NNEIG+1)
       !
       comm % bound_size(1) = 1
       jneig = 0
       do ineig = 1,comm2 % nneig
          if( bound_size(ineig) > 0 ) then
             jneig = jneig + 1
             comm % bound_size(jneig+1) = &
                  comm % bound_size(jneig) + bound_size(ineig)
          end if
       end do

       if( associated(bound_perm) ) deallocate( bound_perm )
       if( associated(bound_invp) ) deallocate( bound_invp )
       if( associated(bound_size) ) deallocate( bound_size )
       if( associated(neights)    ) deallocate( neights    )
       
       if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor

    end if

  end subroutine sub

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-12-09
  !> @brief   Collapse a communicator
  !> @details Collapse a communicator in case a neighbor appears
  !>          more than once 
  !> 
  !-----------------------------------------------------------------------

  subroutine collapse(comm,MEMORY_COUNTER)
    
    class(comm_data_par),          intent(inout) :: comm            !< Input communicator
    integer(8),          optional, intent(inout) :: MEMORY_COUNTER(2)  !< Memory counter
    integer(ip)                                  :: ineig,jneig
    integer(ip)                                  :: nneig_new
    integer(ip)                                  :: dom_i,dom_j,ii,kk,ii1,ii2
    integer(ip), pointer                         :: bound_size(:)
    integer(ip), pointer                         :: neights(:)
    type(i1p),   pointer                         :: bound_perm(:)
    type(i1p),   pointer                         :: bound_invp(:)
    integer(8)                                   :: memor(2)
    
    memor = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
    
    nullify(bound_size)
    nullify(bound_perm)
    nullify(bound_invp)
    nullify(neights)

    call memory_alloca(memor,'BOUND_SIZE',vacal,bound_size,comm % nneig+1)
    call memory_alloca(memor,'BOUND_PERM',vacal,bound_perm,comm % bound_dim)
    call memory_alloca(memor,'BOUND_INVP',vacal,bound_invp,comm % bound_dim)
    call memory_alloca(memor,'NEIGHTS'   ,vacal,neights   ,comm % nneig)
    !
    ! Merge repeated neighbors
    !
    nneig_new = 0
    do ineig = 1,comm % nneig
       dom_i = comm % neights(ineig)
       if( dom_i /= 0 ) then
          kk = 0
          do jneig = ineig,comm % nneig
             dom_j = comm % neights(jneig)
             if( dom_j == dom_i ) then
                if( jneig == ineig ) then
                   nneig_new = nneig_new + 1
                   neights(nneig_new) = dom_i
                end if
                ii1                   = comm % bound_size(jneig)   
                ii2                   = comm % bound_size(jneig+1)-1
                ii                    = ii2-ii1+1
                bound_size(nneig_new) = bound_size(nneig_new) + ii
                if( bound_size(nneig_new) > 0 ) then
                   call memory_resize(memor,'BOUND_PERM % L',vacal,bound_perm(nneig_new) % l,bound_size(nneig_new))
                   call memory_resize(memor,'BOUND_INVP % L',vacal,bound_invp(nneig_new) % l,bound_size(nneig_new))
                   bound_perm(nneig_new) % l(kk+1:kk+ii) = comm % bound_perm(ii1:ii2)
                   bound_invp(nneig_new) % l(kk+1:kk+ii) = comm % bound_invp(ii1:ii2)
                end if
                kk                    = kk + ii
                comm % neights(jneig) = 0
             end if
          end do
       end if
    end do
    comm % nneig = nneig_new
    
    call memory_deallo(memor,'BOUND_SIZE',vacal,comm % bound_size)
    call memory_deallo(memor,'BOUND_PERM',vacal,comm % neights   )
    call memory_alloca(memor,'BOUND_SIZE',vacal,comm % bound_size,comm % nneig+1)
    call memory_alloca(memor,'BOUND_PERM',vacal,comm % neights   ,comm % nneig)
    !
    ! Reconstruct communicator
    !
    if( comm % nneig > 0 ) then
       comm % bound_size(1) = 1
       kk = 0
       do ineig = 1,comm % nneig
          comm % neights(ineig) = neights(ineig) 
          do ii = 1,bound_size(ineig)
             kk = kk + 1
             comm % bound_perm(kk) = bound_perm(ineig) % l(ii)
             comm % bound_invp(kk) = bound_invp(ineig) % l(ii)
          end do
          comm % bound_size(ineig+1) = comm % bound_size(ineig) + bound_size(ineig)
       end do
    end if

    call memory_deallo(memor,'BOUND_SIZE',vacal,bound_size)
    call memory_deallo(memor,'BOUND_PERM',vacal,bound_perm)
    call memory_deallo(memor,'BOUND_INVP',vacal,bound_invp)
    
    call memory_deallo(memor,'BOUND_PERM',vacal,neights)
    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor
    
  end subroutine collapse
  
end module def_kintyp_comm
!> @}
