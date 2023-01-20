!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @defgroup Parall_Toolbox
!> Toolbox for parallelization tools
!> @{
!> @name    Parallelization toolbox
!> @file    mod_parall.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!> @brief   Communication arrays
!> @details Tools for communciation arrays
!>
!------------------------------------------------------------------------

module mod_communication_arrays

  use def_kintyp_basic,      only : ip,rp,lg
  use def_kintyp_comm,       only : comm_data_par
  use def_kintyp_comm,       only : comm_data_par_basic
  use def_parall,            only : par_memor
  use mod_memory_basic,      only : memory_alloca
  use mod_memory_basic,      only : memory_deallo
  use mod_memory_basic,      only : memory_copy
  use mod_optional_argument, only : optional_argument
  use def_mpi
  implicit none
  private
  !
  ! Interfaces
  !
  interface PAR_INITIALIZE_COMMUNICATION_ARRAY
     module procedure PAR_INITIALIZE_COMMUNICATION_ARRAY_s,     &
          &           PAR_INITIALIZE_COMMUNICATION_ARRAY_BASIC, &
          &           PAR_INITIALIZE_COMMUNICATION_ARRAY_1
  end interface PAR_INITIALIZE_COMMUNICATION_ARRAY

  interface PAR_DEALLOCATE_COMMUNICATION_ARRAY
     module procedure PAR_DEALLOCATE_COMMUNICATION_ARRAY_s,     &
          &           PAR_DEALLOCATE_COMMUNICATION_ARRAY_BASIC, &
          &           PAR_DEALLOCATE_COMMUNICATION_ARRAY_1
  end interface PAR_DEALLOCATE_COMMUNICATION_ARRAY

  interface PAR_ALLOCATE_COMMUNICATION_ARRAY
     module procedure PAR_ALLOCATE_COMMUNICATION_ARRAY_s,     &
          &           PAR_ALLOCATE_COMMUNICATION_ARRAY_BASIC
  end interface PAR_ALLOCATE_COMMUNICATION_ARRAY

  interface PAR_COPY_COMMUNICATION_ARRAY
     module procedure PAR_COPY_COMMUNICATION_ARRAY_s,     &
          &           PAR_COPY_COMMUNICATION_ARRAY_BASIC
  end interface PAR_COPY_COMMUNICATION_ARRAY

  public :: PAR_INITIALIZE_COMMUNICATION_ARRAY
  public :: PAR_DEALLOCATE_COMMUNICATION_ARRAY
  public :: PAR_ALLOCATE_COMMUNICATION_ARRAY
  public :: PAR_COPY_COMMUNICATION_ARRAY
  public :: PAR_SUB_COMMUNICATION_ARRAY
  public :: PAR_POINT_COMMUNICATION_ARRAY
  
contains

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    27/03/2018
  !> @brief   Deallocate communication array
  !> @details Deallocate communication array
  !
  !----------------------------------------------------------------------

  subroutine PAR_DEALLOCATE_COMMUNICATION_ARRAY_1(COMM,memor_opt,PAR_COMM_OPT,COMM_NAME,INITIALIZE)

    type(comm_data_par), intent(inout), pointer  :: COMM(:)          !< Input communicator
    integer(8),          intent(inout), optional :: memor_opt(2)  !< Memory counter
    logical(lg),         intent(in),    optional :: PAR_COMM_OPT
    character(*),        intent(in),    optional :: COMM_NAME
    logical(lg),         intent(in),    optional :: INITIALIZE
    integer(ip)                                  :: ii

    if( associated(COMM) ) then
       do ii = lbound(COMM,1),ubound(COMM,1)
          call PAR_DEALLOCATE_COMMUNICATION_ARRAY_s(COMM(ii),memor_opt,PAR_COMM_OPT,COMM_NAME,INITIALIZE)
       end do
       deallocate(COMM)
    end if

  end subroutine PAR_DEALLOCATE_COMMUNICATION_ARRAY_1
  
  subroutine PAR_DEALLOCATE_COMMUNICATION_ARRAY_s(COMM,memor_opt,PAR_COMM_OPT,COMM_NAME,INITIALIZE)

    type(comm_data_par), intent(inout)           :: COMM          !< Input communicator
    integer(8),          intent(inout), optional :: memor_opt(2)  !< Memory counter
    logical(lg),         intent(in),    optional :: PAR_COMM_OPT
    character(*),        intent(in),    optional :: COMM_NAME
    logical(lg),         intent(in),    optional :: INITIALIZE
    integer(8)                                   :: memor(2)
    logical(lg)                                  :: PAR_COMM
    integer(ip)                                  :: ii
    character(20)                                :: my_comm_name
    logical(lg)                                  :: if_initialize

    if_initialize = optional_argument(.false.,INITIALIZE)
    my_comm_name  = optional_argument(trim(COMM % NAME),COMM_NAME)

    PAR_COMM   = .true.

    if( present(PAR_COMM_OPT) ) PAR_COMM = PAR_COMM_OPT

    if( present(memor_opt) ) then
       memor = memor_opt
    else
       memor = 0_8
    end if

    if( if_initialize ) then
       COMM % nneig               = 0
       COMM % npoi1               = 0
       COMM % npoi2               = 0
       COMM % npoi3               = 0
       COMM % npoin               = 0
       COMM % bound_dim           = 0
       COMM % nedg1               = 0
       COMM % nedg2               = 0
       COMM % nedg3               = 0

       COMM % offset_npoin        = 0
       COMM % offset_nelem        = 0
       COMM % offset_nboun        = 0

       COMM % bedge_dim           = 0
       COMM % lsend_dim           = 0
       COMM % lrecv_dim           = 0
       COMM % lscat_dim           = 0
       COMM % matrix_nzdom        = 0
       COMM % bface_dim           = 0
       COMM % full_row_send_dim   = 0
       COMM % full_row_recv_nneig = 0
       COMM % full_row_recv_dim   = 0
       COMM % ghost_send_node_dim = 0
       COMM % ghost_recv_node_dim = 0
       COMM % ghost_send_elem_dim = 0
       COMM % ghost_recv_elem_dim = 0
       COMM % ghost_send_boun_dim = 0
       COMM % ghost_recv_boun_dim = 0
       if(  PAR_COMM ) then
          ! OJO NULL
          !#ifndef MPI_OFF
          !         COMM % PAR_COMM_WORLD            = MPI_COMM_NULL
          !#else
          !          COMM % PAR_COMM_WORLD % MPI_VAL  = 0
          !#endif
          COMM % PAR_COMM_WORLD            = PAR_COMM_NULL
       end if

       COMM % SIZE4               = -1_4
       COMM % RANK4               = -1_4

    end if

    call memory_deallo(memor,trim(my_comm_name)//' % NEIGHTS'                ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % neights              )
    call memory_deallo(memor,trim(my_comm_name)//' % NEIGHTS_ORDERED'        ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % neights_ordered      )
    call memory_deallo(memor,trim(my_comm_name)//' % BOUND_SIZE'             ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bound_size           )
    call memory_deallo(memor,trim(my_comm_name)//' % BOUND_PERM'             ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bound_perm           )
    call memory_deallo(memor,trim(my_comm_name)//' % BOUND_INVP'             ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bound_invp           )
    call memory_deallo(memor,trim(my_comm_name)//' % PERM_ORDERED'           ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % perm_ordered         )
    call memory_deallo(memor,trim(my_comm_name)//' % FULL_ROW_SEND_NEIGHTS'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % full_row_send_neights)
    call memory_deallo(memor,trim(my_comm_name)//' % FULL_ROW_SEND_SIZE'     ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % full_row_send_size   )
    call memory_deallo(memor,trim(my_comm_name)//' % FULL_ROW_SEND_PERM'     ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % full_row_send_perm   )
    call memory_deallo(memor,trim(my_comm_name)//' % FULL_ROW_RECV_NEIGHTS'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % full_row_recv_neights)
    call memory_deallo(memor,trim(my_comm_name)//' % FULL_ROW_RECV_SIZE'     ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % full_row_recv_size   )
    call memory_deallo(memor,trim(my_comm_name)//' % FULL_ROW_RECV_PERM'     ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % full_row_recv_perm   )

    call memory_deallo(memor,trim(my_comm_name)//' % BOUND_SCAL'             ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bound_scal         )
    call memory_deallo(memor,trim(my_comm_name)//' % BOUND_MULTIPLICITY'     ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bound_multiplicity )
    call memory_deallo(memor,trim(my_comm_name)//' % BOUND_OWNER_RANK'       ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bound_owner_rank   )

    call memory_deallo(memor,trim(my_comm_name)//' % BFACE_SIZE'             ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bface_size         )
    call memory_deallo(memor,trim(my_comm_name)//' % BFACE_PERM'             ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bface_perm         )

    call memory_deallo(memor,trim(my_comm_name)//' % NODE_NUMBER_IN_OWNER'   ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % node_number_in_owner)

    if( associated(COMM % bound_matrix) ) then
       do ii = 1,size(COMM % bound_matrix)
          call memory_deallo(memor,trim(my_comm_name)//' % BOUND_MATRIX % JA'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bound_matrix(ii) % ja)
          call memory_deallo(memor,trim(my_comm_name)//' % BOUND_MATRIX % NZDOM_ii'    ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bound_matrix(ii) % nzdom_ii)
          COMM % bound_matrix(ii) % nzdom = 0
       end do
       deallocate(COMM % bound_matrix)
    end if
    if( associated(COMM % bound_mat_halo_send) ) then
       do ii = 1,size(COMM % bound_mat_halo_send)
          call memory_deallo(memor,trim(my_comm_name)//' % BOUND_MAT_HALO_SEND % JA'      ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bound_mat_halo_send(ii) % ja)
          call memory_deallo(memor,trim(my_comm_name)//' % BOUND_MAT_HALO_SEND % NZDOM_II','PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bound_mat_halo_send(ii) % nzdom_ii)
          COMM % bound_mat_halo_send(ii) % nzdom = 0
       end do
       deallocate(COMM % bound_mat_halo_send)
    end if
    if( associated(COMM % bound_mat_halo_recv) ) then
       do ii = 1,size(COMM % bound_mat_halo_recv)
          call memory_deallo(memor,trim(my_comm_name)//' % BOUND_MAT_HALO_RECV % JA'       ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bound_mat_halo_recv(ii) % ja)
          call memory_deallo(memor,trim(my_comm_name)//' % BOUND_MAT_HALO_RECV % NZDOM_II' ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bound_mat_halo_recv(ii) % nzdom_ii)
          COMM % bound_mat_halo_recv(ii) % nzdom = 0
       end do
       deallocate(COMM % bound_mat_halo_recv)
    end if

    call memory_deallo(memor,trim(my_comm_name)//' % BEDGE_SIZE'             ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bedge_size           )
    call memory_deallo(memor,trim(my_comm_name)//' % BEDGE_PERM'             ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bedge_perm           )
    call memory_deallo(memor,trim(my_comm_name)//' % BEDGE_ADJA'             ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bedge_adja           )
    call memory_deallo(memor,trim(my_comm_name)//' % BEDGE_SCAL'             ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bedge_scal           )
    call memory_deallo(memor,trim(my_comm_name)//' % BEDGE_MULTIPLICITY'     ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bedge_multiplicity   )
    call memory_deallo(memor,trim(my_comm_name)//' % BEDGE_OWNER_RANK'       ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bedge_owner_rank     )
    call memory_deallo(memor,trim(my_comm_name)//' % LSEND_SIZE'             ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % lsend_size           )
    call memory_deallo(memor,trim(my_comm_name)//' % LSEND_PERM'             ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % lsend_perm           )
    call memory_deallo(memor,trim(my_comm_name)//' % LRECV_SIZE'             ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % lrecv_size           )
    call memory_deallo(memor,trim(my_comm_name)//' % LRECV_PERM'             ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % lrecv_perm           )
    call memory_deallo(memor,trim(my_comm_name)//' % LSCAT_PERM'             ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % lscat_perm           )
    call memory_deallo(memor,trim(my_comm_name)//' % MATRIX_IA'              ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % matrix_ia            )
    call memory_deallo(memor,trim(my_comm_name)//' % MATRIX_JA'              ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % matrix_ja            )
    call memory_deallo(memor,trim(my_comm_name)//' % MATRIX_AA'              ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % matrix_aa            )

    call memory_deallo(memor,trim(my_comm_name)//' % GHOST_SEND_NODE_SIZE'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % ghost_send_node_size   )
    call memory_deallo(memor,trim(my_comm_name)//' % GHOST_SEND_NODE_PERM'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % ghost_send_node_perm   )
    call memory_deallo(memor,trim(my_comm_name)//' % GHOST_RECV_NODE_SIZE'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % ghost_recv_node_size   )
    call memory_deallo(memor,trim(my_comm_name)//' % GHOST_RECV_NODE_PERM'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % ghost_recv_node_perm   )
    call memory_deallo(memor,trim(my_comm_name)//' % GHOST_SEND_ELEM_SIZE'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % ghost_send_elem_size   )
    call memory_deallo(memor,trim(my_comm_name)//' % GHOST_SEND_ELEM_PERM'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % ghost_send_elem_perm   )
    call memory_deallo(memor,trim(my_comm_name)//' % GHOST_RECV_ELEM_SIZE'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % ghost_recv_elem_size   )
    call memory_deallo(memor,trim(my_comm_name)//' % GHOST_RECV_ELEM_PERM'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % ghost_recv_elem_perm   )
    call memory_deallo(memor,trim(my_comm_name)//' % GHOST_SEND_BOUN_SIZE'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % ghost_send_boun_size   )
    call memory_deallo(memor,trim(my_comm_name)//' % GHOST_SEND_BOUN_PERM'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % ghost_send_boun_perm   )
    call memory_deallo(memor,trim(my_comm_name)//' % GHOST_RECV_BOUN_SIZE'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % ghost_recv_boun_size   )
    call memory_deallo(memor,trim(my_comm_name)//' % GHOST_RECV_BOUN_PERM'  ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % ghost_recv_boun_perm   )

    if( present(memor_opt) ) then
       memor_opt = memor
    end if

  end subroutine PAR_DEALLOCATE_COMMUNICATION_ARRAY_s

  subroutine PAR_DEALLOCATE_COMMUNICATION_ARRAY_BASIC(COMM,memor_opt,PAR_COMM_OPT,COMM_NAME,INITIALIZE)

    type(comm_data_par_basic), intent(inout)           :: COMM          !< Input communicator
    integer(8),                intent(inout), optional :: memor_opt(2)  !< Memory counter
    logical(lg),               intent(in),    optional :: PAR_COMM_OPT
    character(*),              intent(in),    optional :: COMM_NAME
    logical(lg),               intent(in),    optional :: INITIALIZE

    call COMM % deallo(memor_opt,PAR_COMM_OPT,COMM_NAME,INITIALIZE)

  end subroutine PAR_DEALLOCATE_COMMUNICATION_ARRAY_BASIC

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    06/03/2014
  !> @brief   Initialize communication array
  !> @details Initialize communication array
  !
  !----------------------------------------------------------------------

  subroutine PAR_INITIALIZE_COMMUNICATION_ARRAY_s(COMM,PAR_COMM_OPT,COMM_NAME)
    type(comm_data_par), intent(inout)        :: COMM     !< Input communicator
    logical(lg),         intent(in), optional :: PAR_COMM_OPT
    character(*),        intent(in), optional :: COMM_NAME
    logical(lg)                               :: PAR_COMM
    
    PAR_COMM = .true.
    if( present(PAR_COMM_OPT) ) PAR_COMM = PAR_COMM_OPT

    COMM % name                =  optional_argument('COMM',COMM_NAME)
    COMM % nneig               =  0
    COMM % nneig_1             =  0
    COMM % nneig_2             =  0
    COMM % npoi1               =  0
    COMM % npoi2               =  0
    COMM % npoi3               =  0
    COMM % npoin               =  0
    COMM % nedg1               =  0
    COMM % nedg2               =  0
    COMM % nedg3               =  0
    
    COMM % offset_npoin        =  0
    COMM % offset_nelem        =  0
    COMM % offset_nboun        =  0
    
    COMM % bound_dim           =  0
    COMM % bedge_dim           =  0
    COMM % lsend_dim           =  0
    COMM % lrecv_dim           =  0
    COMM % lscat_dim           =  0
    COMM % matrix_nzdom        =  0
    COMM % bface_dim           =  0

    COMM % full_row_send_nneig =  0
    COMM % full_row_send_dim   = -1
    COMM % full_row_recv_nneig =  0
    COMM % full_row_recv_dim   = -1

    COMM % ghost_send_elem_dim = -1
    COMM % ghost_recv_elem_dim = -1
    COMM % ghost_send_node_dim = -1
    COMM % ghost_recv_node_dim = -1
    COMM % ghost_send_boun_dim = -1
    COMM % ghost_recv_boun_dim = -1

    COMM % SIZE4               = -1_4
    COMM % RANK4               = -1_4

    if( PAR_COMM ) then
       ! OJO NULL
!#ifndef MPI_OFF
!       COMM % PAR_COMM_WORLD            = MPI_COMM_NULL
!#else
!       COMM % PAR_COMM_WORLD % MPI_VAL  = -1
       !#endif
       COMM % PAR_COMM_WORLD            = PAR_COMM_NULL
    end if
    nullify(COMM % neights)
    nullify(COMM % neights_ordered)
    nullify(COMM % perm_ordered)
    nullify(COMM % bound_size)
    nullify(COMM % bound_perm)
    nullify(COMM % bound_invp)
    nullify(COMM % bound_scal)
    nullify(COMM % bound_multiplicity)
    nullify(COMM % bound_owner_rank)
    nullify(COMM % node_number_in_owner)
    nullify(COMM % bound_matrix)
    nullify(COMM % bound_mat_halo_send)
    nullify(COMM % bound_mat_halo_recv)

    nullify(COMM % bedge_size)
    nullify(COMM % bedge_perm)
    nullify(COMM % bedge_adja)
    nullify(COMM % bedge_scal)
    nullify(COMM % bedge_multiplicity)
    nullify(COMM % bedge_owner_rank)

    nullify(COMM % lsend_size)
    nullify(COMM % lrecv_size)
    nullify(COMM % lsend_perm)
    nullify(COMM % lrecv_perm)
    nullify(COMM % lscat_perm)
    nullify(COMM % matrix_ia)
    nullify(COMM % matrix_ja)
    nullify(COMM % matrix_aa)
    nullify(COMM % bface_size)
    nullify(COMM % bface_perm)

    nullify(COMM % full_row_send_neights)
    nullify(COMM % full_row_send_size)
    nullify(COMM % full_row_send_perm)
    nullify(COMM % full_row_recv_neights)
    nullify(COMM % full_row_recv_size)
    nullify(COMM % full_row_recv_perm)

    nullify(COMM % ghost_send_elem_size)
    nullify(COMM % ghost_send_elem_perm)
    nullify(COMM % ghost_recv_elem_size)
    nullify(COMM % ghost_recv_elem_perm)
    nullify(COMM % ghost_send_node_size)
    nullify(COMM % ghost_send_node_perm)
    nullify(COMM % ghost_recv_node_size)
    nullify(COMM % ghost_recv_node_perm)
    nullify(COMM % ghost_send_boun_size)
    nullify(COMM % ghost_send_boun_perm)
    nullify(COMM % ghost_recv_boun_size)
    nullify(COMM % ghost_recv_boun_perm)

  end subroutine PAR_INITIALIZE_COMMUNICATION_ARRAY_s

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    06/03/2014
  !> @brief   Initialize communication array
  !> @details Initialize communication array
  !
  !----------------------------------------------------------------------

  subroutine PAR_INITIALIZE_COMMUNICATION_ARRAY_BASIC(COMM,PAR_COMM_OPT,COMM_NAME)
    type(comm_data_par_basic), intent(inout)        :: COMM     !< Input communicator
    logical(lg),               intent(in), optional :: PAR_COMM_OPT
    character(*),              intent(in), optional :: COMM_NAME
    logical(lg)                                     :: PAR_COMM

    call COMM % init()
    
    PAR_COMM = .true.
    if( present(PAR_COMM_OPT) ) PAR_COMM = PAR_COMM_OPT

    COMM % SIZE4               = -1_4
    COMM % RANK4               = -1_4

    if( PAR_COMM ) then
          ! OJO NULL
!#ifndef MPI_OFF
!       COMM % PAR_COMM_WORLD            = MPI_COMM_NULL
!#else
!       COMM % PAR_COMM_WORLD % MPI_VAL  = -1
!#endif
       COMM % PAR_COMM_WORLD            = PAR_COMM_NULL
    end if
    
  end subroutine PAR_INITIALIZE_COMMUNICATION_ARRAY_BASIC

  subroutine PAR_INITIALIZE_COMMUNICATION_ARRAY_1(COMM,PAR_COMM_OPT,COMM_NAME)
    type(comm_data_par), pointer, intent(inout)          :: COMM(:)     !< Input communicator
    logical(lg),                  intent(in),   optional :: PAR_COMM_OPT
    character(*),                   intent(in), optional :: COMM_NAME
    integer(ip)                                          :: icomm

    if( associated(COMM) ) then
       do icomm = lbound(COMM,1),ubound(COMM,1)
          call PAR_INITIALIZE_COMMUNICATION_ARRAY_s(COMM(icomm),PAR_COMM_OPT,COMM_NAME)
       end do
    end if

  end subroutine PAR_INITIALIZE_COMMUNICATION_ARRAY_1

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    06/03/2014
  !> @brief   Copy a communication array
  !> @details Copy a communication array to another
  !
  !----------------------------------------------------------------------

  subroutine PAR_POINT_COMMUNICATION_ARRAY(COMM_IN,COMM_OUT)
    type(comm_data_par), intent(in)    :: COMM_IN(:)     !< Input communicator
    type(comm_data_par), intent(inout) :: COMM_OUT(:)    !< Output communicator

    COMM_OUT(1) % nneig                 =  COMM_IN(1) % nneig
    COMM_OUT(1) % nneig_1               =  COMM_IN(1) % nneig_1
    COMM_OUT(1) % nneig_2               =  COMM_IN(1) % nneig_2
    COMM_OUT(1) % npoi1                 =  COMM_IN(1) % npoi1
    COMM_OUT(1) % npoi2                 =  COMM_IN(1) % npoi2
    COMM_OUT(1) % npoi3                 =  COMM_IN(1) % npoi3
    COMM_OUT(1) % npoin                 =  COMM_IN(1) % npoin
    COMM_OUT(1) % bound_dim             =  COMM_IN(1) % bound_dim
    COMM_OUT(1) % nedg1                 =  COMM_IN(1) % nedg1
    COMM_OUT(1) % nedg2                 =  COMM_IN(1) % nedg2
    COMM_OUT(1) % nedg3                 =  COMM_IN(1) % nedg3

    COMM_OUT(1) % offset_npoin          =  COMM_IN(1) % offset_npoin
    COMM_OUT(1) % offset_nelem          =  COMM_IN(1) % offset_nelem
    COMM_OUT(1) % offset_nboun          =  COMM_IN(1) % offset_nboun
    
    COMM_OUT(1) % lsend_dim             =  COMM_IN(1) % lsend_dim
    COMM_OUT(1) % lrecv_dim             =  COMM_IN(1) % lrecv_dim
    COMM_OUT(1) % lscat_dim             =  COMM_IN(1) % lscat_dim
    COMM_OUT(1) % matrix_nzdom          =  COMM_IN(1) % matrix_nzdom
    COMM_OUT(1) % bface_dim             =  COMM_IN(1) % bface_dim
    COMM_OUT(1) % ghost_send_elem_dim   =  COMM_IN(1) % ghost_send_elem_dim
    COMM_OUT(1) % ghost_recv_elem_dim   =  COMM_IN(1) % ghost_recv_elem_dim
    COMM_OUT(1) % ghost_send_node_dim   =  COMM_IN(1) % ghost_send_node_dim
    COMM_OUT(1) % ghost_recv_node_dim   =  COMM_IN(1) % ghost_recv_node_dim
    COMM_OUT(1) % ghost_send_boun_dim   =  COMM_IN(1) % ghost_send_boun_dim
    COMM_OUT(1) % ghost_recv_boun_dim   =  COMM_IN(1) % ghost_recv_boun_dim

    COMM_OUT(1) % PAR_COMM_WORLD        =  COMM_IN(1) % PAR_COMM_WORLD
    COMM_OUT(1) % SIZE4                 =  COMM_IN(1) % SIZE4
    COMM_OUT(1) % RANK4                 =  COMM_IN(1) % RANK4
    COMM_OUT(1) % neights               => COMM_IN(1) % neights
    COMM_OUT(1) % neights_ordered       => COMM_IN(1) % neights_ordered
    COMM_OUT(1) % perm_ordered          => COMM_IN(1) % perm_ordered

    COMM_OUT(1) % bound_size            => COMM_IN(1) % bound_size
    COMM_OUT(1) % bound_perm            => COMM_IN(1) % bound_perm
    COMM_OUT(1) % bound_invp            => COMM_IN(1) % bound_invp
    COMM_OUT(1) % bound_scal            => COMM_IN(1) % bound_scal
    COMM_OUT(1) % bound_multiplicity    => COMM_IN(1) % bound_multiplicity
    COMM_OUT(1) % bound_owner_rank      => COMM_IN(1) % bound_owner_rank
    COMM_OUT(1) % node_number_in_owner  => COMM_IN(1) % node_number_in_owner
    COMM_OUT(1) % bound_matrix          => COMM_IN(1) % bound_matrix
    COMM_OUT(1) % bound_mat_halo_send   => COMM_IN(1) % bound_mat_halo_send
    COMM_OUT(1) % bound_mat_halo_recv   => COMM_IN(1) % bound_mat_halo_recv

    COMM_OUT(1) % bedge_size            => COMM_IN(1) % bedge_size
    COMM_OUT(1) % bedge_perm            => COMM_IN(1) % bedge_perm
    COMM_OUT(1) % bedge_adja            => COMM_IN(1) % bedge_adja
    COMM_OUT(1) % bedge_scal            => COMM_IN(1) % bedge_scal
    COMM_OUT(1) % bedge_dim             =  COMM_IN(1) % bedge_dim
    COMM_OUT(1) % bedge_multiplicity    => COMM_IN(1) % bedge_multiplicity
    COMM_OUT(1) % bedge_owner_rank      => COMM_IN(1) % bedge_owner_rank

    COMM_OUT(1) % lsend_size            => COMM_IN(1) % lsend_size
    COMM_OUT(1) % lrecv_size            => COMM_IN(1) % lrecv_size
    COMM_OUT(1) % lsend_perm            => COMM_IN(1) % lsend_perm
    COMM_OUT(1) % lrecv_perm            => COMM_IN(1) % lrecv_perm
    COMM_OUT(1) % lscat_perm            => COMM_IN(1) % lscat_perm

    COMM_OUT(1) % matrix_ia             => COMM_IN(1) % matrix_ia
    COMM_OUT(1) % matrix_ja             => COMM_IN(1) % matrix_ja
    COMM_OUT(1) % matrix_aa             => COMM_IN(1) % matrix_aa

    COMM_OUT(1) % full_row_send_neights => COMM_IN(1) % full_row_send_neights
    COMM_OUT(1) % full_row_send_size    => COMM_IN(1) % full_row_send_size
    COMM_OUT(1) % full_row_send_perm    => COMM_IN(1) % full_row_send_perm
    COMM_OUT(1) % full_row_recv_neights => COMM_IN(1) % full_row_recv_neights
    COMM_OUT(1) % full_row_recv_size    => COMM_IN(1) % full_row_recv_size
    COMM_OUT(1) % full_row_recv_perm    => COMM_IN(1) % full_row_recv_perm

    COMM_OUT(1) % bface_size            => COMM_IN(1) % bface_size
    COMM_OUT(1) % bface_perm            => COMM_IN(1) % bface_perm

    COMM_OUT(1) % ghost_send_elem_size  => COMM_IN(1) % ghost_send_elem_size
    COMM_OUT(1) % ghost_send_elem_perm  => COMM_IN(1) % ghost_send_elem_perm
    COMM_OUT(1) % ghost_recv_elem_size  => COMM_IN(1) % ghost_recv_elem_size
    COMM_OUT(1) % ghost_recv_elem_perm  => COMM_IN(1) % ghost_recv_elem_perm
    COMM_OUT(1) % ghost_send_node_size  => COMM_IN(1) % ghost_send_node_size
    COMM_OUT(1) % ghost_send_node_perm  => COMM_IN(1) % ghost_send_node_perm
    COMM_OUT(1) % ghost_recv_node_size  => COMM_IN(1) % ghost_recv_node_size
    COMM_OUT(1) % ghost_recv_node_perm  => COMM_IN(1) % ghost_recv_node_perm

    COMM_OUT(1) % ghost_send_boun_size  => COMM_IN(1) % ghost_send_boun_size
    COMM_OUT(1) % ghost_send_boun_perm  => COMM_IN(1) % ghost_send_boun_perm
    COMM_OUT(1) % ghost_recv_boun_size  => COMM_IN(1) % ghost_recv_boun_size
    COMM_OUT(1) % ghost_recv_boun_perm  => COMM_IN(1) % ghost_recv_boun_perm

  end subroutine PAR_POINT_COMMUNICATION_ARRAY

  subroutine PAR_ALLOCATE_COMMUNICATION_ARRAY_s(COMM_IN,memor_opt,COMM_NAME)

    type(comm_data_par), intent(inout)           :: COMM_IN       !< Input communicator
    integer(8),          intent(inout), optional :: memor_opt(2)  !< Memory counter
    character(len=*),    intent(in),    optional :: COMM_NAME
    integer(8)                                   :: memor(2)
    character(30)                                :: my_comm_name

    my_comm_name  = optional_argument(trim(COMM_IN % NAME),COMM_NAME)
    
    if( present(memor_opt) ) then
       memor = memor_opt
    else
       memor = 0_8
    end if

    call memory_alloca(memor,trim(my_comm_name)//' % NEIGHTS',             'par_allocate_communication_array',COMM_IN % neights               , COMM_IN % nneig )
    call memory_alloca(memor,trim(my_comm_name)//' % BOUND_SIZE',          'par_allocate_communication_array',COMM_IN % bound_size            , COMM_IN % nneig+1 )
    call memory_alloca(memor,trim(my_comm_name)//' % BOUND_PERM',          'par_allocate_communication_array',COMM_IN % bound_perm            , COMM_IN % bound_dim )
    call memory_alloca(memor,trim(my_comm_name)//' % BOUND_INVP',          'par_allocate_communication_array',COMM_IN % bound_invp            , COMM_IN % bound_dim )

!!$    call memory_alloca(memor,'COMM_IN % NEIGHTS_ORDERED',     'par_allocate_communication_array',COMM_IN % neights_ordered       , COMM_IN % nneig )
!!$    call memory_alloca(memor,'COMM_IN % PERM_ORDERED',        'par_allocate_communication_array',COMM_IN % perm_ordered          , COMM_IN % nneig )

!!$    call memory_alloca(memor,'COMM_IN % BOUND_SCAL',          'par_allocate_communication_array',COMM_IN % bound_scal            , COMM_IN %
!!$    call memory_alloca(memor,'COMM_IN % BOUND_MULTIPLICITY',  'par_allocate_communication_array',COMM_IN % bound_owner_rank      , COMM_IN %
!!$    call memory_alloca(memor,'COMM_IN % NODE_NUMBER_IN_OWNER','par_allocate_communication_array',COMM_IN % node_number_in_owner  , COMM_IN %
!!$
!!$    if( associated(COMM_IN % bound_matrix) )        call runend('YOU SHOULD CODE THESE LINES...')
!!$    if( associated(COMM_IN % bound_mat_halo_send) ) call runend('YOU SHOULD CODE THESE LINES...')
!!$    if( associated(COMM_IN % bound_mat_halo_recv) ) call runend('YOU SHOULD CODE THESE LINES...')
!!$    !call memory_alloca(memor,'COMM % LSEND_SIZE','par_allocate_communication_array',COMM_IN % bound_matrix          , COMM_OUT % bound_matrix,'DO_NOT_DEALLOCATE')
!!$    !call memory_alloca(memor,'COMM % LSEND_SIZE','par_allocate_communication_array',COMM_IN % bound_mat_halo_send   , COMM_OUT % bound_mat_halo_send,'DO_NOT_DEALLOCATE')
!!$    !call memory_alloca(memor,'COMM % LSEND_SIZE','par_allocate_communication_array',COMM_IN % bound_mat_halo_recv   , COMM_OUT % bound_mat_halo_recv,'DO_NOT_DEALLOCATE')
!!$
!!$    call memory_alloca(memor,'COMM_IN % BEDGE_SIZE'          ,'par_allocate_communication_array',COMM_IN % bedge_size            ,
!!$    call memory_alloca(memor,'COMM_IN % BEDGE_PERM',          'par_allocate_communication_array',COMM_IN % bedge_perm            ,
!!$    call memory_alloca(memor,'COMM_IN % BEDGE_ADJA',          'par_allocate_communication_array',COMM_IN % bedge_adja            ,
!!$    call memory_alloca(memor,'COMM_IN % BEDGE_SCAL',          'par_allocate_communication_array',COMM_IN % bedge_scal            ,
!!$    call memory_alloca(memor,'COMM_IN % BEDGE_MULTIPLICITY',  'par_allocate_communication_array',COMM_IN % bedge_multiplicity    ,
!!$    call memory_alloca(memor,'COMM_IN % BEDGE_OWNER_RANK',    'par_allocate_communication_array',COMM_IN % bedge_owner_rank      ,
!!$
!!$    call memory_alloca(memor,'COMM_IN % LSEND_SIZE',          'par_allocate_communication_array',COMM_IN % lsend_size            ,
!!$    call memory_alloca(memor,'COMM_IN % LRECV_SIZE',          'par_allocate_communication_array',COMM_IN % lrecv_size            ,
!!$    call memory_alloca(memor,'COMM_IN % LSEND_PERM',          'par_allocate_communication_array',COMM_IN % lsend_perm            ,
!!$    call memory_alloca(memor,'COMM_IN % LRECV_PERM',          'par_allocate_communication_array',COMM_IN % lrecv_perm            ,
!!$    call memory_alloca(memor,'COMM_IN % LSCAT_PERM',          'par_allocate_communication_array',COMM_IN % lscat_perm            ,
!!$
!!$    call memory_alloca(memor,'COMM_IN % MATRIX_IA',           'par_allocate_communication_array',COMM_IN % matrix_ia             ,
!!$    call memory_alloca(memor,'COMM_IN % MATRIX_JA',           'par_allocate_communication_array',COMM_IN % matrix_ja             ,
!!$    call memory_alloca(memor,'COMM_IN % MATRIX_AA',           'par_allocate_communication_array',COMM_IN % matrix_aa             ,
!!$
!!$    call memory_alloca(memor,'COMM_IN %% FULL_ROW_SEND_NEIGHTS','par_allocate_communication_array',COMM_IN % full_row_send_neights,)
!!$    call memory_alloca(memor,'COMM_IN %% FULL_ROW_SEND_SIZE'   ,'par_allocate_communication_array',COMM_IN % full_row_send_size   ,)
!!$    call memory_alloca(memor,'COMM_IN %% FULL_ROW_SEND_PERM'   ,'par_allocate_communication_array',COMM_IN % full_row_send_perm   ,)
!!$    call memory_alloca(memor,'COMM_IN %% FULL_ROW_RECV_NEIGHTS','par_allocate_communication_array',COMM_IN % full_row_recv_neights,)
!!$    call memory_alloca(memor,'COMM_IN %% FULL_ROW_RECV_SIZE'   ,'par_allocate_communication_array',COMM_IN % full_row_recv_size   ,)
!!$    call memory_alloca(memor,'COMM_IN %% FULL_ROW_RECV_PERM'   ,'par_allocate_communication_array',COMM_IN % full_row_recv_perm   ,)
!!$
!!$    call memory_alloca(memor,'COMM_IN % BFACE_SIZE',          'par_allocate_communication_array',COMM_IN % bface_size            ,
!!$    call memory_alloca(memor,'COMM_IN % BFACE_PERM',          'par_allocate_communication_array',COMM_IN % bface_perm            ,
!!$
!!$    call memory_alloca(memor,'COMM_IN % GHOST_SEND_ELEM_SIZE','par_allocate_communication_array',COMM_IN % ghost_send_elem_size  ,
!!$    call memory_alloca(memor,'COMM_IN % GHOST_SEND_ELEM_PERM','par_allocate_communication_array',COMM_IN % ghost_send_elem_perm  ,
!!$    call memory_alloca(memor,'COMM_IN % GHOST_RECV_ELEM_SIZE','par_allocate_communication_array',COMM_IN % ghost_recv_elem_size  ,
!!$    call memory_alloca(memor,'COMM_IN % GHOST_RECV_ELEM_PERM','par_allocate_communication_array',COMM_IN % ghost_recv_elem_perm  ,
!!$    call memory_alloca(memor,'COMM_IN % GHOST_SEND_NODE_SIZE','par_allocate_communication_array',COMM_IN % ghost_send_node_size  ,
!!$    call memory_alloca(memor,'COMM_IN % GHOST_SEND_NODE_PERM','par_allocate_communication_array',COMM_IN % ghost_send_node_perm  ,
!!$    call memory_alloca(memor,'COMM_IN % GHOST_RECV_NODE_SIZE','par_allocate_communication_array',COMM_IN % ghost_recv_node_size  ,
!!$    call memory_alloca(memor,'COMM_IN % GHOST_RECV_NODE_PERM','par_allocate_communication_array',COMM_IN % ghost_recv_node_perm  ,
!!$
!!$    call memory_alloca(memor,'COMM_IN % GHOST_SEND_BOUN_SIZE','par_allocate_communication_array',COMM_IN % ghost_send_boun_size  , )
!!$    call memory_alloca(memor,'COMM_IN % GHOST_SEND_BOUN_PERM','par_allocate_communication_array',COMM_IN % ghost_send_boun_perm  , )
!!$    call memory_alloca(memor,'COMM_IN % GHOST_RECV_BOUN_SIZE','par_allocate_communication_array',COMM_IN % ghost_recv_boun_size  , )
!!$    call memory_alloca(memor,'COMM_IN % GHOST_RECV_BOUN_PERM','par_allocate_communication_array',COMM_IN % ghost_recv_boun_perm  , )

    if( present(memor_opt) ) memor_opt = memor

  end subroutine PAR_ALLOCATE_COMMUNICATION_ARRAY_s

  subroutine PAR_ALLOCATE_COMMUNICATION_ARRAY_BASIC(COMM_IN,memor_opt,COMM_NAME)

    type(comm_data_par_basic), intent(inout)           :: COMM_IN       !< Input communicator
    integer(8),                intent(inout), optional :: memor_opt(2)  !< Memory counter
    character(len=*),          intent(in),    optional :: COMM_NAME
    integer(8)                                         :: memor(2)
    character(30)                                      :: my_comm_name

    memor         = optional_argument((/0_8,0_8/),memor_opt)
    my_comm_name  = optional_argument(trim(COMM_IN % NAME),COMM_NAME)
    
    call memory_alloca(memor,trim(my_comm_name)//' % NEIGHTS',             'par_allocate_communication_array_basic',COMM_IN % neights               , COMM_IN % nneig )
    call memory_alloca(memor,trim(my_comm_name)//' % BOUND_SIZE',          'par_allocate_communication_array_basic',COMM_IN % bound_size            , COMM_IN % nneig+1 )
    call memory_alloca(memor,trim(my_comm_name)//' % BOUND_PERM',          'par_allocate_communication_array_basic',COMM_IN % bound_perm            , COMM_IN % bound_dim )

    if( present(memor_opt) ) memor_opt = memor

  end subroutine PAR_ALLOCATE_COMMUNICATION_ARRAY_BASIC

  subroutine PAR_COPY_COMMUNICATION_ARRAY_BASIC(COMM_IN,COMM_OUT,memor_opt,COMM_NAME,COPY_NAME)

    type(comm_data_par_basic), intent(inout)           :: COMM_IN       !< Input communicator
    type(comm_data_par),       intent(inout)           :: COMM_OUT      !< Output communicator
    character(*),              intent(in),    optional :: COMM_NAME
    character(*),              intent(in),    optional :: COPY_NAME
    integer(8),                intent(inout), optional :: memor_opt(2)  !< Memory counter
    integer(8)                                         :: memor(2)
    character(len=:), allocatable                      :: my_comm_name
    character(len=:), allocatable                      :: my_copy_name

    memor        = optional_argument((/0_8,0_8/),memor_opt)
    my_comm_name = optional_argument(trim(COMM_IN  % NAME),COMM_NAME)
    my_copy_name = optional_argument(trim(COMM_OUT % NAME),COPY_NAME)
    
    COMM_OUT % nneig                 =  COMM_IN % nneig
    COMM_OUT % bound_dim             =  COMM_IN % bound_dim

    COMM_OUT % offset_npoin          =  COMM_IN % offset_npoin
    COMM_OUT % offset_nelem          =  COMM_IN % offset_nelem    

    COMM_OUT % PAR_COMM_WORLD        =  COMM_IN % PAR_COMM_WORLD
    COMM_OUT % SIZE4                 =  COMM_IN % SIZE4
    COMM_OUT % RANK4                 =  COMM_IN % RANK4

    call memory_copy(memor,trim(my_comm_name)//' % NEIGHTS',   'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % neights    , COMM_OUT % neights,   'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % NEIGHTS'   )
    call memory_copy(memor,trim(my_comm_name)//' % BOUND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bound_size , COMM_OUT % bound_size,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_SIZE')
    call memory_copy(memor,trim(my_comm_name)//' % BOUND_PERM','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bound_perm , COMM_OUT % bound_perm,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_PERM')

    if( present(memor_opt) ) memor_opt = memor

    if( allocated(my_comm_name) ) deallocate(my_comm_name)
    if( allocated(my_copy_name) ) deallocate(my_copy_name)

  end subroutine PAR_COPY_COMMUNICATION_ARRAY_BASIC
  
  subroutine PAR_COPY_COMMUNICATION_ARRAY_s(COMM_IN,COMM_OUT,memor_opt,MINIMUM_DATA,COMM_NAME,COPY_NAME)

    type(comm_data_par), intent(inout)           :: COMM_IN       !< Input communicator
    type(comm_data_par), intent(inout)           :: COMM_OUT      !< Output communicator
    character(*),        intent(in),    optional :: COMM_NAME
    character(*),        intent(in),    optional :: COPY_NAME
    integer(8),          intent(inout), optional :: memor_opt(2)  !< Memory counter
    logical(lg),         intent(in),    optional :: MINIMUM_DATA
    integer(8)                                   :: memor(2)
    integer(ip)                                  :: ii,isize
    logical(lg)                                  :: if_minimum_data
    character(len=:), allocatable                :: my_comm_name
    character(len=:), allocatable                :: my_copy_name

    memor           = optional_argument((/0_8,0_8/),memor_opt)
    my_comm_name    = optional_argument(trim(COMM_IN  % NAME),COMM_NAME)
    my_copy_name    = optional_argument(trim(COMM_OUT % NAME),COPY_NAME)
    if_minimum_data = optional_argument(.false.,MINIMUM_DATA)

    COMM_OUT % nneig                 =  COMM_IN % nneig
    COMM_OUT % nneig_1               =  COMM_IN % nneig_1
    COMM_OUT % nneig_2               =  COMM_IN % nneig_2
    COMM_OUT % npoi1                 =  COMM_IN % npoi1
    COMM_OUT % npoi2                 =  COMM_IN % npoi2
    COMM_OUT % npoi3                 =  COMM_IN % npoi3
    COMM_OUT % npoin                 =  COMM_IN % npoin
    COMM_OUT % bound_dim             =  COMM_IN % bound_dim

    COMM_OUT % nedg1                 =  COMM_IN % nedg1
    COMM_OUT % nedg2                 =  COMM_IN % nedg2
    COMM_OUT % nedg3                 =  COMM_IN % nedg3
    
    COMM_OUT % offset_npoin          =  COMM_IN % offset_npoin
    COMM_OUT % offset_nelem          =  COMM_IN % offset_nelem
    COMM_OUT % offset_nboun          =  COMM_IN % offset_nboun
    
    COMM_OUT % lsend_dim             =  COMM_IN % lsend_dim
    COMM_OUT % lrecv_dim             =  COMM_IN % lrecv_dim
    COMM_OUT % lscat_dim             =  COMM_IN % lscat_dim
    COMM_OUT % matrix_nzdom          =  COMM_IN % matrix_nzdom
    COMM_OUT % bface_dim             =  COMM_IN % bface_dim
    COMM_OUT % ghost_send_elem_dim   =  COMM_IN % ghost_send_elem_dim
    COMM_OUT % ghost_recv_elem_dim   =  COMM_IN % ghost_recv_elem_dim
    COMM_OUT % ghost_send_node_dim   =  COMM_IN % ghost_send_node_dim
    COMM_OUT % ghost_recv_node_dim   =  COMM_IN % ghost_recv_node_dim
    COMM_OUT % ghost_send_boun_dim   =  COMM_IN % ghost_send_boun_dim
    COMM_OUT % ghost_recv_boun_dim   =  COMM_IN % ghost_recv_boun_dim

    COMM_OUT % PAR_COMM_WORLD        =  COMM_IN % PAR_COMM_WORLD
    COMM_OUT % SIZE4                 =  COMM_IN % SIZE4
    COMM_OUT % RANK4                 =  COMM_IN % RANK4

    COMM_OUT % bedge_dim             =  COMM_IN % bedge_dim

    call memory_copy(memor,trim(my_comm_name)//' % NEIGHTS',         'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % neights         , COMM_OUT % neights,         'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % NEIGHTS'   )
    call memory_copy(memor,trim(my_comm_name)//' % NEIGHTS_ORDERED', 'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % neights_ordered , COMM_OUT % neights_ordered, 'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % NEIGHTS_ORDERED' )
    call memory_copy(memor,trim(my_comm_name)//' % PERM_ORDERED',    'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % perm_ordered    , COMM_OUT % perm_ordered,    'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % PERM_ORDERED'  )
                                                                                                                                                                                                                            
    call memory_copy(memor,trim(my_comm_name)//' % BOUND_SIZE',      'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bound_size      , COMM_OUT % bound_size,      'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_SIZE'   )
    call memory_copy(memor,trim(my_comm_name)//' % BOUND_PERM',      'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bound_perm      , COMM_OUT % bound_perm,      'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_PERM'   )
    call memory_copy(memor,trim(my_comm_name)//' % BOUND_INVP',      'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bound_invp      , COMM_OUT % bound_invp,      'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_INVP'   )

    if( .not. if_minimum_data ) then 
       call memory_copy(memor,trim(my_comm_name)//' % BOUND_SCAL',          'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bound_scal            , COMM_OUT % bound_scal,          'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_SCAL'          )
       call memory_copy(memor,trim(my_comm_name)//' % BOUND_MULTIPLICITY',  'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bound_multiplicity    , COMM_OUT % bound_multiplicity,  'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_MULTIPLICITY'  )
       call memory_copy(memor,trim(my_comm_name)//' % BOUND_OWNER_RANK',    'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bound_owner_rank      , COMM_OUT % bound_owner_rank,    'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_OWNER_RANK'    )
       call memory_copy(memor,trim(my_comm_name)//' % NODE_NUMBER_IN_OWNER','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % node_number_in_owner  , COMM_OUT % node_number_in_owner,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % NODE_NUMBER_IN_OWNER')

       if( associated(COMM_IN % bound_matrix) ) then
          isize = size(COMM_IN % bound_matrix)
          allocate(COMM_OUT % bound_matrix(isize))
          do ii = 1,isize
             call memory_copy(memor,trim(my_comm_name)//' % BOUND_MATRIX % JA',      'PAR_COPY_COMMUNICATION_ARRAY',&
                  & COMM_IN % bound_matrix(ii) % ja,             COMM_OUT % bound_matrix(ii) % ja,              'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_MATRIX % JA')
             call memory_copy(memor,trim(my_comm_name)//' % BOUND_MATRIX % NZDOM_II','PAR_COPY_COMMUNICATION_ARRAY',&
                  & COMM_IN % bound_matrix(ii) % nzdom_ii,       COMM_OUT % bound_matrix(ii) % nzdom_ii,        'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_MATRIX % NZDOM_II')
             COMM_OUT % bound_matrix(ii) % nzdom  =  COMM_IN % bound_matrix(ii) % nzdom
          end do
       end if

       if( associated(COMM_IN % bound_mat_halo_send) ) then
          isize = size(COMM_IN % bound_mat_halo_send)
          allocate(COMM_OUT % bound_mat_halo_send(isize))
          do ii = 1,isize
             call memory_copy(memor,trim(my_comm_name)//' % BOUND_MAT_HALO_SEND % JA','PAR_COPY_COMMUNICATION_ARRAY',&
                  & COMM_IN % bound_mat_halo_send(ii) % ja,             COMM_OUT % bound_mat_halo_send(ii) % ja,              'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_MAT_HALO_SEND % JA')
             call memory_copy(memor,trim(my_comm_name)//' % BOUND_MAT_HALO_SEND % NZDOM_II',        'PAR_COPY_COMMUNICATION_ARRAY',&
                  & COMM_IN % bound_mat_halo_send(ii) % nzdom_ii,       COMM_OUT % bound_mat_halo_send(ii) % nzdom_ii,        'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_MAT_HALO_SEND % NZDOM_II')
             COMM_OUT % bound_mat_halo_send(ii) % nzdom  =  COMM_IN % bound_mat_halo_send(ii) % nzdom
          end do
       end if

       if( associated(COMM_IN % bound_mat_halo_recv) ) then
          isize = size(COMM_IN % bound_mat_halo_recv)
          allocate(COMM_OUT % bound_mat_halo_recv(isize))
          do ii = 1,isize
             call memory_copy(memor,trim(my_comm_name)//' % BOUND_MAT_HALO_RECV % JA','PAR_COPY_COMMUNICATION_ARRAY',&
                  & COMM_IN % bound_mat_halo_recv(ii) % ja,             COMM_OUT % bound_mat_halo_recv(ii) % ja,              'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_MAT_HALO_RECV % JA')
             call memory_copy(memor,trim(my_comm_name)//' % BOUND_MAT_HALO_RECV % NZDOM_II',        'PAR_COPY_COMMUNICATION_ARRAY',&
                  & COMM_IN % bound_mat_halo_recv(ii) % nzdom_ii,       COMM_OUT % bound_mat_halo_recv(ii) % nzdom_ii,        'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BOUND_MAT_HALO_RECV % NZDOM_II')
             COMM_OUT % bound_mat_halo_recv(ii) % nzdom  =  COMM_IN % bound_mat_halo_recv(ii) % nzdom
          end do
       end if

       call memory_copy(memor,trim(my_comm_name)//' % BEDGE_SIZE'          ,'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bedge_size            , COMM_OUT % bedge_size,          'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BEDGE_SIZE'          )
       call memory_copy(memor,trim(my_comm_name)//' % BEDGE_PERM',          'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bedge_perm            , COMM_OUT % bedge_perm,          'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BEDGE_PERM'          )
       call memory_copy(memor,trim(my_comm_name)//' % BEDGE_ADJA',          'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bedge_adja            , COMM_OUT % bedge_adja,          'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BEDGE_ADJA'          )
       call memory_copy(memor,trim(my_comm_name)//' % BEDGE_SCAL',          'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bedge_scal            , COMM_OUT % bedge_scal,          'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BEDGE_SCAL'          )
       call memory_copy(memor,trim(my_comm_name)//' % BEDGE_MULTIPLICITY',  'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bedge_multiplicity    , COMM_OUT % bedge_multiplicity,  'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BEDGE_MULTIPLICITY'  )
       call memory_copy(memor,trim(my_comm_name)//' % BEDGE_OWNER_RANK',    'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bedge_owner_rank      , COMM_OUT % bedge_owner_rank,    'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BEDGE_OWNER_RANK'    )
                                                                                                                                                                                                                                                  
       call memory_copy(memor,trim(my_comm_name)//' % LSEND_SIZE',          'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % lsend_size            , COMM_OUT % lsend_size,          'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % LSEND_SIZE'          )
       call memory_copy(memor,trim(my_comm_name)//' % LRECV_SIZE',          'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % lrecv_size            , COMM_OUT % lrecv_size,          'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % LRECV_SIZE'          )
       call memory_copy(memor,trim(my_comm_name)//' % LSEND_PERM',          'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % lsend_perm            , COMM_OUT % lsend_perm,          'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % LSEND_PERM'          )
       call memory_copy(memor,trim(my_comm_name)//' % LRECV_PERM',          'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % lrecv_perm            , COMM_OUT % lrecv_perm,          'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % LRECV_PERM'          )
       call memory_copy(memor,trim(my_comm_name)//' % LSCAT_PERM',          'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % lscat_perm            , COMM_OUT % lscat_perm,          'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % LSCAT_PERM'          )
                                                                                                                                                                                                                  
       call memory_copy(memor,trim(my_comm_name)//' % MATRIX_IA',           'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % matrix_ia             , COMM_OUT % matrix_ia,           'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % MATRIX_IA'           )
       call memory_copy(memor,trim(my_comm_name)//' % MATRIX_JA',           'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % matrix_ja             , COMM_OUT % matrix_ja,           'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % MATRIX_JA'           )
       call memory_copy(memor,trim(my_comm_name)//' % MATRIX_AA',           'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % matrix_aa             , COMM_OUT % matrix_aa,           'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % MATRIX_AA'           )
                                                                                                                                                                                                           
       call memory_copy(memor,trim(my_comm_name)//' % FULL_ROW_SEND_NEIGHTS','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % full_row_send_neights,COMM_OUT % full_row_send_neights,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % FULL_ROW_SEND_NEIGHTS')
       call memory_copy(memor,trim(my_comm_name)//' % FULL_ROW_SEND_SIZE'   ,'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % full_row_send_size   ,COMM_OUT % full_row_send_size   ,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % FULL_ROW_SEND_SIZE'   )
       call memory_copy(memor,trim(my_comm_name)//' % FULL_ROW_SEND_PERM'   ,'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % full_row_send_perm   ,COMM_OUT % full_row_send_perm   ,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % FULL_ROW_SEND_PERM'   )
       call memory_copy(memor,trim(my_comm_name)//' % FULL_ROW_RECV_NEIGHTS','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % full_row_recv_neights,COMM_OUT % full_row_recv_neights,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % FULL_ROW_RECV_NEIGHTS')
       call memory_copy(memor,trim(my_comm_name)//' % FULL_ROW_RECV_SIZE'   ,'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % full_row_recv_size   ,COMM_OUT % full_row_recv_size   ,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % FULL_ROW_RECV_SIZE'   )
       call memory_copy(memor,trim(my_comm_name)//' % FULL_ROW_RECV_PERM'   ,'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % full_row_recv_perm   ,COMM_OUT % full_row_recv_perm   ,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % FULL_ROW_RECV_PERM'   )
                                                                                                                                                                                                                   
       call memory_copy(memor,trim(my_comm_name)//' % BFACE_SIZE',          'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bface_size            , COMM_OUT % bface_size,          'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BFACE_SIZE'          )
       call memory_copy(memor,trim(my_comm_name)//' % BFACE_PERM',          'PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bface_perm            , COMM_OUT % bface_perm,          'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % BFACE_PERM'          )
                                                                                                                                                                                                                
       call memory_copy(memor,trim(my_comm_name)//' % GHOST_SEND_ELEM_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % ghost_send_elem_size  , COMM_OUT % ghost_send_elem_size,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % GHOST_SEND_ELEM_SIZE')
       call memory_copy(memor,trim(my_comm_name)//' % GHOST_SEND_ELEM_PERM','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % ghost_send_elem_perm  , COMM_OUT % ghost_send_elem_perm,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % GHOST_SEND_ELEM_PERM')
       call memory_copy(memor,trim(my_comm_name)//' % GHOST_RECV_ELEM_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % ghost_recv_elem_size  , COMM_OUT % ghost_recv_elem_size,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % GHOST_RECV_ELEM_SIZE')
       call memory_copy(memor,trim(my_comm_name)//' % GHOST_RECV_ELEM_PERM','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % ghost_recv_elem_perm  , COMM_OUT % ghost_recv_elem_perm,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % GHOST_RECV_ELEM_PERM')
       call memory_copy(memor,trim(my_comm_name)//' % GHOST_SEND_NODE_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % ghost_send_node_size  , COMM_OUT % ghost_send_node_size,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % GHOST_SEND_NODE_SIZE')
       call memory_copy(memor,trim(my_comm_name)//' % GHOST_SEND_NODE_PERM','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % ghost_send_node_perm  , COMM_OUT % ghost_send_node_perm,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % GHOST_SEND_NODE_PERM')
       call memory_copy(memor,trim(my_comm_name)//' % GHOST_RECV_NODE_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % ghost_recv_node_size  , COMM_OUT % ghost_recv_node_size,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % GHOST_RECV_NODE_SIZE')
       call memory_copy(memor,trim(my_comm_name)//' % GHOST_RECV_NODE_PERM','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % ghost_recv_node_perm  , COMM_OUT % ghost_recv_node_perm,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % GHOST_RECV_NODE_PERM')
                                                                                                                                                                                                            
       call memory_copy(memor,trim(my_comm_name)//' % GHOST_SEND_BOUN_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % ghost_send_boun_size  , COMM_OUT % ghost_send_boun_size,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % GHOST_SEND_BOUN_SIZE')
       call memory_copy(memor,trim(my_comm_name)//' % GHOST_SEND_BOUN_PERM','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % ghost_send_boun_perm  , COMM_OUT % ghost_send_boun_perm,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % GHOST_SEND_BOUN_PERM')
       call memory_copy(memor,trim(my_comm_name)//' % GHOST_RECV_BOUN_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % ghost_recv_boun_size  , COMM_OUT % ghost_recv_boun_size,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % GHOST_RECV_BOUN_SIZE')
       call memory_copy(memor,trim(my_comm_name)//' % GHOST_RECV_BOUN_PERM','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % ghost_recv_boun_perm  , COMM_OUT % ghost_recv_boun_perm,'DO_NOT_DEALLOCATE',COPY_NAME=my_copy_name//' % GHOST_RECV_BOUN_PERM')

    end if

    if( present(memor_opt) ) memor_opt = memor
    if( allocated(my_comm_name) ) deallocate(my_comm_name)
    if( allocated(my_copy_name) ) deallocate(my_copy_name)

  end subroutine PAR_COPY_COMMUNICATION_ARRAY_s

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    06/03/2014
  !> @brief   Create a sub communication array
  !> @details Given a communicaiton array structure, create a sub
  !           communication structure for a subset of nodes
  !
  !----------------------------------------------------------------------

  subroutine PAR_SUB_COMMUNICATION_ARRAY(COMM_IN,COMM_OUT,mask,COMM_NAME)
    type(comm_data_par),          intent(in)    :: COMM_IN        !< Input communicator
    type(comm_data_par),          intent(inout) :: COMM_OUT       !< Output communicator
    integer(ip),         pointer, intent(in)    :: mask(:)        !< mask (=1 to consider node)
    character(*),        optional,intent(in)    :: COMM_NAME
    integer(ip),         pointer                :: bound_perm(:)
    integer(ip),         pointer                :: bound_invp(:)
    integer(ip),         pointer                :: bound_size(:)
    integer(ip),         pointer                :: neights(:)
    integer(ip)                                 :: nneig,bound_dim
    integer(ip)                                 :: nneig_1,nneig_2
    integer(ip)                                 :: ineig,jj,jneig
    integer(ip)                                 :: ipoin
    character(50)                               :: my_comm_name

    COMM_OUT % PAR_COMM_WORLD = COMM_IN % PAR_COMM_WORLD
    COMM_OUT % SIZE4          = COMM_IN % SIZE4
    COMM_OUT % RANK4          = COMM_IN % RANK4
    
    if( present(COMM_NAME) ) then
       my_comm_name = trim(COMM_NAME)
    else
       my_comm_name = 'COMM'
    end if

    if( COMM_IN % bound_dim > 0 ) then

       nullify(bound_perm)
       nullify(bound_invp)
       nullify(bound_size)
       nullify(neights)

       allocate( bound_perm(COMM_IN % bound_dim) )
       allocate( bound_invp(COMM_IN % bound_dim) )
       allocate( bound_size(COMM_IN % nneig)     )
       allocate( neights   (COMM_IN % nneig)     )

       nneig     = 0
       bound_dim = 0

       do jj = 1,COMM_IN % bound_dim
          bound_perm(jj) = COMM_IN % bound_perm(jj)
          bound_invp(jj) = COMM_IN % bound_invp(jj)
       end do
       do ineig = 1,COMM_IN % nneig
          bound_size(ineig) = 0
       end do

       do ineig = 1,COMM_IN % nneig
          do jj = COMM_IN % bound_size(ineig),COMM_IN % bound_size(ineig+1)-1
             ipoin = COMM_IN % bound_perm(jj)
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
       COMM_OUT % nneig     = nneig
       COMM_OUT % bound_dim = bound_dim
       COMM_OUT % npoi1     = COMM_IN % npoi1
       COMM_OUT % npoi2     = COMM_IN % npoi2
       COMM_OUT % npoi3     = COMM_IN % npoi3
       COMM_OUT % npoin     = COMM_IN % npoin
       call memory_alloca(par_memor,trim(my_comm_name)//' % NEIGHTS'        ,'PAR_SUB_COMMUNICATION_ARRAY',COMM_OUT % neights,nneig)
       call memory_alloca(par_memor,trim(my_comm_name)//' % NEIGHTS_ORDERED','PAR_SUB_COMMUNICATION_ARRAY',COMM_OUT % neights_ordered,nneig)
       call memory_alloca(par_memor,trim(my_comm_name)//' % PERM_ORDERED'   ,'PAR_SUB_COMMUNICATION_ARRAY',COMM_OUT % perm_ordered,nneig)
       call memory_alloca(par_memor,trim(my_comm_name)//' % BOUND_PERM'     ,'PAR_SUB_COMMUNICATION_ARRAY',COMM_OUT % bound_perm,bound_dim)
       call memory_alloca(par_memor,trim(my_comm_name)//' % BOUND_INVP'     ,'PAR_SUB_COMMUNICATION_ARRAY',COMM_OUT % bound_invp,bound_dim)
       call memory_alloca(par_memor,trim(my_comm_name)//' % BOUND_SIZE'     ,'PAR_SUB_COMMUNICATION_ARRAY',COMM_OUT % bound_size,nneig+1)

       !allocate( COMM_OUT % neights(nneig) )
       !allocate( COMM_OUT % bound_perm(bound_dim) )
       !allocate( COMM_OUT % bound_size(nneig+1) )
       !
       ! Permutation array BOUND_PERM(1:BOUND_DIM)
       !
       jneig = 0
       bound_dim = 0
       do ineig = 1,COMM_IN % nneig
          if( bound_size(ineig) > 0 ) then

             jneig = jneig + 1
             COMM_OUT % neights(jneig)         = COMM_IN % neights(ineig)
             COMM_OUT % neights_ordered(jneig) = COMM_IN % neights_ordered(ineig)
             COMM_OUT % perm_ordered(jneig)    = COMM_IN % perm_ordered(ineig)

             do jj = COMM_IN % bound_size(ineig),COMM_IN % bound_size(ineig+1)-1
                ipoin = bound_perm(jj)
                if( ipoin > 0 ) then
                   bound_dim = bound_dim + 1
                   COMM_OUT % bound_perm(bound_dim) = ipoin
                   COMM_OUT % bound_invp(bound_dim) = ipoin
                end if
             end do

          end if
       end do
       ineig_1: do ineig = 1,COMM_IN % nneig
          if( ineig > int(COMM_IN % RANK4,ip) ) then
             nneig_1 = ineig-1
             nneig_2 = ineig
             exit ineig_1
          end if
       end do ineig_1
       !
       ! Construct linked list BOUND_SIZE(1:NNEIG+1)
       !
       COMM_OUT % bound_size(1) = 1
       jneig = 0
       do ineig = 1,COMM_IN % nneig
          if( bound_size(ineig) > 0 ) then
             jneig = jneig + 1
             COMM_OUT % bound_size(jneig+1) = &
                  COMM_OUT % bound_size(jneig) + bound_size(ineig)
          end if
       end do

       if( associated(bound_perm) ) deallocate( bound_perm )
       if( associated(bound_invp) ) deallocate( bound_invp )
       if( associated(bound_size) ) deallocate( bound_size )
       if( associated(neights)    ) deallocate( neights    )

    end if

  end subroutine PAR_SUB_COMMUNICATION_ARRAY

end module mod_communication_arrays
!> @}

