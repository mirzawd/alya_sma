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
!> @brief   "The easiest way of making software scalable is to make it sequentially inefficient."
!>          - Bill Gropp, 1999
!>
!>          Michael E. Henderson, Christopher Radcliff Anderson, Stephen L. Lyons
!>          Object Oriented Methods for Interoperable Scientific and Engineering Computing: Proceedings of the 1998 SIAM Workshop
!>          SIAM, Jan 1, 1999 - Technology & Engineering - 321 pages
!>
!> @details ToolBox and definitions for parall.
!>
!>          The variables
!>          -------------
!>
!>          \verbatim
!>
!>          PAR_COMM_UNIVERSE........................................ Communicator of the universe  (Always MPI_COMM_WORLD)
!>          PAR_COMM_MY_CODE ........................................ Communicator of current code
!>          PAR_COMM_MY_CODE_WM ..................................... Communicator of current code without master
!>          PAR_COMM_SHARED_WM ...................................... MPI_COMM_TYPE_SHARED without the master
!>          PAR_COMM_WORLD .......................................... Alya world communicator (MPI_COMM_WORLD or application split of MPI_COMM_WORLD)
!>          PAR_WORLD_SIZE .......................................... World size
!>          PAR_MY_UNIVERSE_RANK..................................... My rank in the universe
!>          PAR_MY_WORLD_RANK_WM .................................... My rank in the world without master
!>          PAR_MY_WORLD_RANK ....................................... My rank in the world
!>          PAR_UNIVERSE_SIZE........................................ Size of the universe
!>          PAR_CODE_SIZE ........................................... Size of my code communicator
!>          PAR_MY_CODE_RANK ........................................ My rank in my code
!>          PAR_MY_PARMETIS_RANK .................................... My rank in the parallel partitioners communicator
!>          PAR_MY_PARMETI2_RANK .................................... My rank in the parallel partitioners communicator + master
!>          PAR_INTEGER ............................................. Length of integer for MPI (ip)
!>
!>          I_AM_IN_COLOR(ICOLO) .................................... if I have color ICOLO (TRUE/FALSE)
!>          PAR_COMM_COLOR(0:MCOLO,0:MCOLO) ......................... Intercolor communicator
!>          PAR_COMM_COLOR_ARRAY(0:MCOLO) ........................... Color communication arrays
!>          PAR_CPU_TO_COLOR(0:PAR_WORLD_SIZE) % L(:) ............... List of colors for each world partition
!>          PAR_COLOR_TO_CPU(0:MCOLO) % L(:) ........................ List of world partitions for each color
!>          PAR_COMM_COLOR_PERM(0:MCOLO,0:MCOLO,0:PAR_WORLD_SIZE) ... Ranks for each communicator
!>          PAR_COMM_WORLD_TO_CODE_PERM(2,0:PAR_WORLD_SIZE) ......... Rank permutation from world to code
!>
!>          MCOLO ................................................... Nb of colors (codes,zones,subds) in the world
!>          MCODE ................................................... Nb of codes in the world
!>          MZONE ................................................... Nb of zones in the world
!>          MSUBD ................................................... Nb of subds in the world
!>          ncolo ................................................... Nb of colors in my code
!>          MAPPS ................................................... Nb of applications in the universe (MPI_COMM_WORLD)
!>
!>          \endverbatim
!>
!>          The functions
!>          -------------
!>
!>          \verbatim
!>
!>          PAR_WORLD_RANK_OF_A_CODE_NEIGHBOR ....................... World rank of a code neighbor
!>          PAR_WORLD_MASTER_RANK_OF_A_CODE ......................... Rank of the master of a code
!>          PAR_COLOR_COUPLING_RANK_TO_WORLD ........................ My rank in the world given my rank in a coupling
!>          PAR_CODE_ZONE_SUBD_TO_COLOR ............................. Mapping (code,zone,subd) => color
!>          PAR_COLOR_TO_CODE ....................................... Mapping color => code
!>          PAR_COLOR_TO_ZONE ....................................... Mapping color => zone
!>          PAR_COLOR_TO_SUBD ....................................... Mapping color => subd
!>          PAR_COLOR_TO_CODE_ZONE_SUBD(3) .......................... Mapping color => (code,zone,subd)
!>          PAR_PART_IN_COLOR ....................................... Check if a partition (rank in world) is in color
!>          PAR_THIS_NODE_IS_MINE ................................... If a node is interior or own boundary
!>          PAR_NODE_NUMBER_OF_PARTITIONS ........................... Number of partitions sharing this node
!>          par_omp_coloring ........................................ Color the elements for openmp
!>
!>          \endverbatim
!>
!>          Communication structure
!>          -----------------------
!>
!>          COMMD % NNEIG .................. Number of neigbors
!>          COMMD % NEIGHTS(:) ............. List of neighbor's ranks
!>          COMMD % BOUND_DIM ............. Total size of interface (repeated nodes)
!>          COMMD % BOUND_SIZE(:) ......... Linked lisT to exchange nodes with my neighbors
!>                                          Number of nodes to exchange with ISUBD: COMMD % BOUND_SIZE(ISUBD+1)-COMMD % BOUND_SIZE(ISUBD)
!>          COMMD % BOUND_PERM(:) ......... Linked list to exchange nodes with my neighbors
!>                                          commd % BOUND_PERM(COMMD % BOUND_SIZE(ISUBD):COMMD % BOUND_SIZE(ISUBD+1)-1) is the
!>                                          list of nodes to exchange with my neighbor ISUBD
!>          COMMD % BOUND_OWNER_RANK(:) ... Rank of the subdomain that owns the node
!>
!>
!------------------------------------------------------------------------

module mod_parall

  use def_kintyp_basic,          only : ip,rp,lg,i1p
  use def_kintyp_domain,         only : ompss_domain
  use def_kintyp_comm,           only : comm_data_par
  use def_kintyp_comm,           only : comm_data_par_basic
  use def_parall,                only : par_memor
  use def_master,                only : ISEQUEN,ISLAVE,IMASTER,INOTMASTER
  use def_master,                only : npoi1,npoi2,npoi3,lninv_loc
  use def_master,                only : intost,gisca,igene,lun_outpu
  use def_master,                only : ioutp,kfl_paral
  use mod_memory,                only : memory_copy
  use mod_memory,                only : memory_alloca
  use mod_memory,                only : memory_deallo
  use mod_memory,                only : memory_size
  use def_mpi
#include "def_mpi.inc"
  
  implicit none 

  private
  
  type comm_redistribution
     type(comm_data_par) :: commd_npoin
     type(comm_data_par) :: commd_nelem
     type(comm_data_par) :: commd_nboun
     type(comm_data_par) :: commd_nbopo
  end type comm_redistribution
  !
  ! MPI communicators
  !
  MY_MPI_COMM                   :: PAR_COMM_CURRENT                  ! Current communicator
  MY_MPI_COMM                   :: PAR_COMM_MY_CODE                  ! Communicator of current code (defined at beginning)
  MY_MPI_COMM                   :: PAR_COMM_MY_CODE_WM               ! Communicator of current code without master
  MY_MPI_COMM                   :: PAR_COMM_SHARED_WM                ! Equivalent to MPI_COMM_TYPE_SHARED without the master
  MY_MPI_COMM                   :: PAR_COMM_VIZ                      ! Communicator for insitu vizulization
  MY_MPI_COMM                   :: PAR_COMM_SIM                      ! Communicator for simulation (same as alya world)
  MY_MPI_COMM                   :: PAR_COMM_SIMVIZ                   ! Communicator combined simulation and viz
  MY_MPI_COMM                   :: PAR_COMM_UNIVERSE                 ! Universe communicator
  MY_MPI_COMM                   :: PAR_COMM_WORLD                    ! Alya World communicator
  MY_MPI_COMM                   :: PAR_MY_PARMETIS_RANK              ! My communicator in the partitioners world
  MY_MPI_COMM                   :: PAR_MY_PARMETI2_RANK              ! My communicator in the partitioners world + master
  MY_MPI_COMM                   :: PAR_COMM_SFC                      ! SFC partition communicator
  MY_MPI_COMM                   :: PAR_COMM_SFC_WM                   ! SFC partition communicator without master
  MY_MPI_COMM                   :: PAR_COMM_MPIO                     ! // READING communicator
  MY_MPI_COMM                   :: PAR_COMM_MPIO_WM                  ! // READING communicator without master
  MY_MPI_COMM,         pointer  :: PAR_COMM_COLOR(:,:)               ! Intercolor MPI communicator
  !
  ! Data types
  !
  MY_MPI_DATATYPE               :: PAR_INTEGER                       ! Integer ip type 
  MY_MPI_DATATYPE               :: PAR_REAL                          ! Real rp type
  !
  ! Rank and size of MPI communicators
  !
  integer(ip)                   :: PAR_MY_UNIVERSE_RANK              ! My rank in the universe
  integer(ip)                   :: PAR_MY_WORLD_RANK                 ! My communicator in the world
  integer(ip)                   :: PAR_MY_CODE_RANK                  ! My communicator in the world
  integer(ip)                   :: PAR_MY_CODE_RANK_WM               ! My communicator in the world
  integer(ip)                   :: PAR_SHARED_RANK_WM                ! My rank in the host
  integer(ip)                   :: PAR_MY_SFC_RANK_WM                ! SFC partition communicator without master
  integer(ip)                   :: PAR_COMM_MPIO_RANK_WM             ! // READING rank without master
  integer(ip)                   :: PAR_UNIVERSE_SIZE                 ! Size of the universe communicator
  integer(ip)                   :: PAR_CODE_SIZE                     ! Size of the code communicator
  integer(ip)                   :: PAR_WORLD_SIZE                    ! Size of the world communicator
  integer(ip)                   :: PAR_COMM_MPIO_WM_SIZE             ! // READING communicator size without master
  integer(ip),         pointer  :: PAR_COMM_COLOR_PERM(:,:,:)        ! Rank communicator permutation array
  !
  ! Color stuffs
  !
  integer(ip),         pointer  :: PAR_COMM_WORLD_TO_CODE_PERM(:,:)  ! Permutation world to code and rank
  logical(lg),         pointer  :: I_AM_IN_COLOR(:)                  ! If I am in color
  !
  ! Alya communicators
  !
  type(comm_data_par), pointer  :: PAR_COMM_COLOR_ARRAY(:)           ! Intercolor communicator
  type(comm_data_par), pointer  :: PAR_COMM_MY_CODE_ARRAY(:)         ! Intercolor communicator
  type(comm_data_par), pointer  :: commd                             ! Generic communication array
  type(comm_redistribution)     :: commd_repartitioning              ! Repartitioning commmunicator
  type(comm_data_par)           :: commd_npoin_from_mpio             ! MPIO npoin communicator
  type(comm_data_par)           :: commd_nelem_from_mpio             ! MPIO nelem communicator
  type(comm_data_par)           :: commd_nboun_from_mpio             ! MPIO nboun communicator  
  type(comm_data_par)           :: commd_host                        ! MPI communicator in host
  !
  ! Color structure
  !
  type(i1p),           pointer  :: PAR_CPU_TO_COLOR(:)               ! Given a CPU, gives the colors and associated rank
  type(i1p),           pointer  :: PAR_COLOR_TO_CPU(:)               ! Given a color, gives the CPUs
  type(i1p),           pointer  :: PAR_COLOR_BIN(:)                  ! Bin structure for colors
  !
  ! Current CODE/ZONE/SUBDOMAIN
  !
  integer(ip)                   :: mcolo                             ! Maximum number of colors
  integer(ip)                   :: mcode                             ! Maximum number of codes (over all processes)
  integer(ip)                   :: mzone                             ! Maximum number of zones (over all processes)
  integer(ip)                   :: msubd                             ! Maximum number of subds (over all processes)
  integer(ip)                   :: ncolo                             ! Number of colors of my code
  integer(ip)                   :: mapps                             ! Number of applications in the universe (MPI_COMM_WORLD)
  !
  ! Bin structure to perform subdomain geometrical queries
  ! and subdomain bounding box
  !
  integer(ip)                   :: par_bin_boxes(3)                  ! # boxes in each direction
  real(rp)                      :: par_bin_comin(3)                  ! Minimum box coordinates
  real(rp)                      :: par_bin_comax(3)                  ! Maximum box coordinates
  integer(ip),         pointer  :: par_bin_size(:,:,:)               ! # partitions per box
  type(i1p),           pointer  :: par_bin_part(:,:,:)               ! Bin structure of world partition
  real(rp),            pointer  :: par_part_comin(:,:)               ! Subdomain minimum coordinates
  real(rp),            pointer  :: par_part_comax(:,:)               ! Subdomain maximum coordinates
  type typ_bin_structure
     integer(ip)                :: boxes(3)                          ! # boxes in each direction
     real(rp)                   :: comin(3)                          ! Minimum box coordinates
     real(rp)                   :: comax(3)                          ! Maximum box coordinates
     integer(ip),      pointer  :: size(:,:,:)                       ! # partitions per box
     type(i1p),        pointer  :: part(:,:,:)                       ! Bin structure of world partition
     real(rp),         pointer  :: part_comin(:,:)                   ! Subdomain minimum coordinates
     real(rp),         pointer  :: part_comax(:,:)                   ! Subdomain maximum coordinates
  end type typ_bin_structure
  !
  ! Others
  !
  integer(ip)                   :: color_target                      ! Target color when using coupling
  integer(ip)                   :: color_source                      ! Source color when using coupling
  integer(ip),        parameter :: lun_outpu_par = 5502              ! Output
  integer(ip),        parameter :: lun_parti_msh = 5507              ! Partition mesh file
  integer(ip),        parameter :: lun_parti_res = 5508              ! Partition result file
  integer(ip),        parameter :: lun_matri_msh = 5509              ! Matrix mesh file
  integer(ip),        parameter :: lun_matri_res = 5510              ! Matrix result file
  integer(ip),        parameter :: lun_conne_par = 5511              ! MPI connectivity
  integer(ip),        parameter :: lun_repar_par = 5512              ! Repartitioning results
  integer(4)                    :: PARMETIS_COMM, PARMETIS_INTERCOMM ! parmetis communicators
  integer(4)                    :: PARMETI2_COMM                     ! parmetis communicators + master
  !
  ! OpenMP
  !
  integer(ip)                   :: par_omp_num_blocks                ! Size of blocks to compute chunks n/par_omp_num_blocks
  integer(ip)                   :: par_omp_granularity               ! Granularity to compute par_omp_num_blocks=par_omp_num_threads*par_omp_granularity
  integer(ip)                   :: par_omp_nelem_chunk   = 0         ! Element loop chunk size
  integer(ip)                   :: par_omp_npoin_chunk               ! Node loop chunk size
  integer(ip)                   :: par_omp_nboun_chunk               ! Boundary loop chunk size
  integer(ip)                   :: par_omp_num_threads               ! Number of openmp threads
  integer(ip),        pointer   :: par_omp_thread_num(:)             ! Thread number
  integer(ip)                   :: par_omp_coloring_alg              ! Coloring algorithm
  integer(ip)                   :: par_omp_partition_alg             ! Partitioning algorithm for OmpSs (use same nomenclature for mesh partitioning)

  integer(ip)                   :: par_omp_num_colors                ! Element: Number of colors
  integer(ip),        pointer   :: par_omp_list_colors(:)            ! Element: Element colors
  integer(ip),        pointer   :: par_omp_ia_colors(:)              ! Element: Linked list IA for colors
  integer(ip),        pointer   :: par_omp_ja_colors(:)              ! Element: Linked list JA for colors

  integer(ip)                   :: par_omp_nboun_num_colors          ! Boundary: Number of colors
  integer(ip),        pointer   :: par_omp_nboun_list_colors(:)      ! Boundary: Element colors
  integer(ip),        pointer   :: par_omp_nboun_ia_colors(:)        ! Boundary: Linked list IA for colors
  integer(ip),        pointer   :: par_omp_nboun_ja_colors(:)        ! Boundary: Linked list JA for colors
  !
  ! One-sided communications
  !
  integer(ip)                          :: par_one_sided
  !
  ! Hybrid parallelization (OpenMP/OmpSs) and vectorization
  !
  integer(ip),        parameter        :: PAR_HYBRID_OFF         = 0
  integer(ip),        parameter        :: PAR_OPENMP_COLORING    = 1
  integer(ip),        parameter        :: PAR_OPENMP_NO_COLORING = 2
  integer(ip),        parameter        :: PAR_OMPSS              = 3
  integer(ip)                          :: par_hybrid
  type typ_list_elements_par
     type(i1p),                pointer :: packs(:)
  end type typ_list_elements_par

  integer(ip)                          :: num_subd_par               ! Element loop: For loops with race condition
  integer(ip),                 pointer :: num_pack_par(:)
  type(typ_list_elements_par), pointer :: list_elements_par(:)

  integer(ip)                          :: num_subd_cpu               ! Element loop: For loops not ported to gpu
  integer(ip),                 pointer :: num_pack_cpu(:)
  type(typ_list_elements_par), pointer :: list_elements_cpu(:)

  integer(ip)                          :: num_subd_norace_par        ! Element loop: For loops without race condition
  integer(ip),                 pointer :: num_pack_norace_par(:)
  type(typ_list_elements_par), pointer :: list_elements_norace_par(:)

  integer(ip)                          :: num_subd_norace_cpu        ! Element loop: For loops without race condition
  integer(ip),                 pointer :: num_pack_norace_cpu(:)
  type(typ_list_elements_par), pointer :: list_elements_norace_cpu(:)

  integer(ip)                          :: num_subd_nboun_par         ! Boundary loop: For loops with race condition
  integer(ip),                 pointer :: num_pack_nboun_par(:)
  type(typ_list_elements_par), pointer :: list_boundaries_par(:)

  integer(ip)                          :: num_subd_norace_nboun_par  ! Boundary loop: For loops without race condition
  integer(ip),                 pointer :: num_pack_norace_nboun_par(:)
  type(typ_list_elements_par), pointer :: list_boundaries_norace_par(:)
  !
  ! Topology
  !
  integer(ip)                          :: par_topo_num_nodes                ! Number of computing nodes
  integer(ip)                          :: par_topo_num_cores_per_node       ! Number of cores per node
  !
  ! Buffers to avoid allocating memory
  !
  real(rp),           pointer          :: sendbuff_rp(:)
  real(rp),           pointer          :: recvbuff_rp(:)
  integer(ip),        pointer          :: sendbuff_ip(:)
  integer(ip),        pointer          :: recvbuff_ip(:)
  !
  ! Partitioning method
  !
  integer(ip),        parameter        :: PAR_METIS4               = 0
  integer(ip),        parameter        :: PAR_SFC                  = 1
  integer(ip),        parameter        :: PAR_ORIENTED_BIN         = 2
  integer(ip),        parameter        :: PAR_FRONTAL              = 3
  integer(ip),        parameter        :: PAR_NUMBERING            = 4
  integer(ip),        parameter        :: PAR_USING_RANK           = 5
  integer(ip),        parameter        :: PAR_RANDOM               = 6
  integer(ip),        parameter        :: PAR_ZOLTAN               = 7

  integer(ip),        parameter        :: PAR_SEQUENTIAL_PARTITION = 0
  integer(ip),        parameter        :: PAR_PARALLEL_PARTITION   = 1

  integer(ip),        parameter        :: PAR_WEIGHT_GAUSS         =  0
  integer(ip),        parameter        :: PAR_WEIGHT_OFF           = -1
  integer(ip),        parameter        :: PAR_WEIGHT_ELEMENT       = -2
  integer(ip),        parameter        :: PAR_WEIGHT_SQUARE        = -3
  integer(ip),        parameter        :: PAR_WEIGHT_MATERIAL      = -4
  integer(ip),        parameter        :: PAR_WEIGHT_STUPID        = -5
  !
  ! Files
  !
  character(150)                       :: fil_parti_msh                     ! Partition mesh
  character(150)                       :: fil_parti_res                     ! Partition result
  !
  ! Interfaces
  !
  interface NODE_IN_NEIGHBOR
     module procedure NODE_IN_NEIGHBOR_P,NODE_IN_NEIGHBOR_s
  end interface NODE_IN_NEIGHBOR

  public :: PAR_TRANSPOSE_COMMUNICATION_ARRAY
  public :: NODE_IN_NEIGHBOR
  !
  ! Subroutines
  !
  public :: par_code_zone_subd_to_color
  public :: par_color_to_code
  public :: par_color_to_zone
  public :: par_color_to_subd
  public :: par_part_in_color
  public :: par_color_coupling_rank_to_world
  public :: par_world_rank_of_a_code_neighbor
  public :: PAR_NODE_NUMBER_OF_PARTITIONS
  !
  ! Types
  !
  public :: typ_bin_structure
  !
  ! Variables
  !
  public :: PAR_COMM_NULL
  public :: PAR_COMM_CURRENT
  public :: PAR_COMM_MY_CODE
  public :: PAR_COMM_MY_CODE_WM
  public :: PAR_COMM_SHARED_WM
  public :: PAR_COMM_VIZ
  public :: PAR_COMM_SIM
  public :: PAR_COMM_SIMVIZ
  public :: PAR_COMM_UNIVERSE
  public :: PAR_UNIVERSE_SIZE
  public :: PAR_MY_UNIVERSE_RANK
  public :: PAR_COMM_WORLD
  public :: PAR_WORLD_SIZE
  public :: PAR_MY_WORLD_RANK
  public :: PAR_CODE_SIZE
  public :: PAR_MY_CODE_RANK
  public :: PAR_MY_CODE_RANK_WM
  public :: PAR_SHARED_RANK_WM
  public :: PAR_MY_PARMETIS_RANK
  public :: PAR_MY_PARMETI2_RANK
  public :: PAR_INTEGER
  public :: PAR_REAL
  public :: PAR_COMM_SFC
  public :: PAR_COMM_SFC_WM
  public :: PAR_MY_SFC_RANK_WM
  public :: PAR_COMM_COLOR
  public :: PAR_COMM_COLOR_PERM
  public :: PAR_COMM_WORLD_TO_CODE_PERM
  public :: I_AM_IN_COLOR
  public :: PAR_COMM_COLOR_ARRAY
  public :: PAR_COMM_MY_CODE_ARRAY
  public :: commd
  public :: commd_npoin_from_mpio
  public :: commd_nelem_from_mpio
  public :: commd_nboun_from_mpio
  public :: commd_host
  
  public :: PAR_CPU_TO_COLOR
  public :: PAR_COLOR_TO_CPU
  public :: PAR_COLOR_BIN
  public :: mcolo
  public :: mcode
  public :: mzone
  public :: msubd
  public :: ncolo
  public :: mapps
  public :: par_bin_boxes
  public :: par_bin_comin
  public :: par_bin_comax
  public :: par_bin_size
  public :: par_bin_part
  public :: par_part_comin
  public :: par_part_comax
  public :: par_memor
  public :: color_target
  public :: color_source
  public :: lun_outpu_par
  public :: lun_parti_msh
  public :: lun_matri_res
  public :: lun_matri_msh
  public :: lun_parti_res
  public :: lun_conne_par
  public :: lun_repar_par
  public :: PARMETIS_COMM
  public :: PARMETI2_COMM
  public :: PARMETIS_INTERCOMM
  public :: par_omp_num_blocks
  public :: par_omp_granularity
  public :: par_omp_nelem_chunk
  public :: par_omp_npoin_chunk
  public :: par_omp_nboun_chunk
  public :: par_omp_num_threads
  public :: par_omp_thread_num

  public :: par_omp_num_colors
  public :: par_omp_list_colors
  public :: par_omp_ia_colors
  public :: par_omp_ja_colors

  public :: par_omp_nboun_num_colors
  public :: par_omp_nboun_list_colors
  public :: par_omp_nboun_ia_colors
  public :: par_omp_nboun_ja_colors

  public :: par_omp_coloring_alg
  public :: par_omp_partition_alg
  public :: ompss_domain
  public :: par_topo_num_nodes                ! Number of computing nodes
  public :: par_topo_num_cores_per_node       ! Number of cores per node

  public :: sendbuff_rp
  public :: recvbuff_rp
  public :: sendbuff_ip
  public :: recvbuff_ip

  public :: PAR_METIS4
  public :: PAR_SFC
  public :: PAR_ZOLTAN
  public :: PAR_ORIENTED_BIN
  public :: PAR_NUMBERING
  public :: PAR_USING_RANK
  public :: PAR_RANDOM
  public :: PAR_SEQUENTIAL_PARTITION
  public :: PAR_PARALLEL_PARTITION
  public :: PAR_WEIGHT_GAUSS
  public :: PAR_WEIGHT_OFF
  public :: PAR_WEIGHT_ELEMENT
  public :: PAR_WEIGHT_SQUARE
  public :: PAR_WEIGHT_MATERIAL
  public :: PAR_WEIGHT_STUPID

  public :: PAR_COMM_MPIO                       ! // READING communicator
  public :: PAR_COMM_MPIO_WM                    ! // READING communicator without master
  public :: PAR_COMM_MPIO_RANK_WM               ! // READING rank without master
  public :: PAR_COMM_MPIO_WM_SIZE

  public :: fil_parti_msh                     ! Partition mesh
  public :: fil_parti_res                     ! Partition result
  !
  ! Communicators
  !
  public :: comm_redistribution
  public :: commd_repartitioning
  !
  ! One-sided
  !
  public :: par_one_sided
  !
  ! Hybrid parallelization
  !
  public :: PAR_OPENMP_COLORING
  public :: PAR_OPENMP_NO_COLORING
  public :: PAR_OMPSS
  public :: PAR_HYBRID_OFF
  public :: par_hybrid

  public :: typ_list_elements_par
  public :: num_subd_par
  public :: num_pack_par
  public :: list_elements_par
  public :: num_subd_cpu
  public :: num_pack_cpu
  public :: list_elements_cpu
  public :: num_subd_norace_par
  public :: num_pack_norace_par
  public :: list_elements_norace_par
  public :: num_subd_norace_cpu
  public :: num_pack_norace_cpu
  public :: list_elements_norace_cpu

  public :: num_subd_nboun_par
  public :: num_pack_nboun_par
  public :: list_boundaries_par
  public :: num_subd_norace_nboun_par
  public :: num_pack_norace_nboun_par
  public :: list_boundaries_norace_par

contains

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/03/2014
  !> @brief   Give the world rank of a code neighbor
  !> @details Given a neighbor INEIG, returns its rank in the MPI_COMM_WORLD
  !
  !----------------------------------------------------------------------

  function par_world_rank_of_a_code_neighbor(ipart,icode)
    integer(ip), intent(in) :: ipart
    integer(ip), intent(in) :: icode
    integer(ip)             :: par_world_rank_of_a_code_neighbor
    integer(ip)             :: ipart_world,ipart_code,jcode

    ipart_world = 0
    par_world_rank_of_a_code_neighbor = -1
    do while( ipart_world <= PAR_WORLD_SIZE-1 )
       jcode      = PAR_COMM_WORLD_TO_CODE_PERM(1,ipart_world)
       ipart_code = PAR_COMM_WORLD_TO_CODE_PERM(2,ipart_world)
       if( ipart_code == ipart .and. icode == jcode ) then
          par_world_rank_of_a_code_neighbor = ipart_world
          ipart_world = PAR_WORLD_SIZE-1
       end if
       ipart_world = ipart_world + 1
    end do
    if( par_world_rank_of_a_code_neighbor == -1 ) &
         call runend('par_world_rank_of_a_code_neighbor: WE ARE IN TROUBLE')

  end function par_world_rank_of_a_code_neighbor

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/03/2014
  !> @brief   Give the rank of the master of a code
  !> @details Given a code ICODE, returns the rank of the master of the
  !>          code in the MPI_COMM_WORLD
  !
  !----------------------------------------------------------------------

  function par_world_master_rank_of_a_code(icode)
    integer(ip), intent(in) :: icode
    integer(ip)             :: par_world_master_rank_of_a_code
    integer(ip)             :: icolo,isize

    icolo = par_code_zone_subd_to_color(icode,0_ip,0_ip)
    isize = 0
    do while( isize <= PAR_WORLD_SIZE-1 )
       par_world_master_rank_of_a_code = PAR_COMM_COLOR_PERM(icolo,icolo,isize)
       if( par_world_master_rank_of_a_code == 0 ) then
          isize = PAR_WORLD_SIZE
       end if
       isize = isize + 1
    end do

  end function par_world_master_rank_of_a_code

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/03/2014
  !> @brief   Give my rank in the world given my rank in a coupling
  !> @details Give my rank in the world given my rank in a coupling
  !
  !----------------------------------------------------------------------

  function par_color_coupling_rank_to_world(my_rank,icolo,jcolo)
    integer(ip), intent(in) :: my_rank
    integer(ip), intent(in) :: icolo
    integer(ip), intent(in) :: jcolo
    integer(ip)             :: par_color_coupling_rank_to_world
    integer(ip)             :: isize

    isize = 0
    do while( isize <= PAR_WORLD_SIZE-1 )
       par_color_coupling_rank_to_world = PAR_COMM_COLOR_PERM(icolo,jcolo,isize)
       if( par_color_coupling_rank_to_world == my_rank ) then
          isize = PAR_WORLD_SIZE
       end if
       isize = isize + 1
    end do

  end function par_color_coupling_rank_to_world

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Mapping (code,zone,subd) => color
  !> @details Given a code, zone and subd, compute the corresponding
  !>          color, where:
  !>          code: 0 to mcode
  !>          zone: 0 to mzone
  !>          subd: 0 to msubd
  !>          3D <=> 1D mapping:
  !>          i = x + Lx * (y+Ly*z)
  !>          x = modulo(i,Lx)
  !>          y = modulo(i/Lx,Ly)
  !>          z = i/(Lx*Ly)
  !
  !----------------------------------------------------------------------

  function par_code_zone_subd_to_color(icode,izone,isubd)
    integer(ip), intent(in) :: icode,izone,isubd
    integer(ip)             :: par_code_zone_subd_to_color

    par_code_zone_subd_to_color = isubd+(msubd+1)*(izone+icode*(mzone+1))

  end function par_code_zone_subd_to_color

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Mapping color => code
  !> @details Given a color, returns the code number
  !
  !----------------------------------------------------------------------

  function par_color_to_code(icolo)
    integer(ip), intent(in) :: icolo
    integer(ip)             :: par_color_to_code

    par_color_to_code = icolo / ((msubd+1)*(mzone+1))

  end function par_color_to_code

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Mapping color => code
  !> @details Given a color, returns the zone number
  !
  !----------------------------------------------------------------------

  function par_color_to_zone(icolo)
    integer(ip), intent(in) :: icolo
    integer(ip)             :: par_color_to_zone

    par_color_to_zone = modulo(icolo/(msubd+1),mzone+1)

  end function par_color_to_zone

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Mapping color => code
  !> @details Given a color, returns the subdomain number
  !
  !----------------------------------------------------------------------

  function par_color_to_subd(icolo)
    integer(ip), intent(in) :: icolo
    integer(ip)             :: par_color_to_subd

    par_color_to_subd = modulo(icolo,msubd+1)

  end function par_color_to_subd

  function par_color_to_code_zone_subd(icolo)
    implicit none
    integer(ip), intent(in) :: icolo
    integer(ip)             :: par_color_to_code_zone_subd(3)

    par_color_to_code_zone_subd(1) = icolo / ((msubd+1)*(mzone+1))
    par_color_to_code_zone_subd(2) = modulo(icolo/(msubd+1),mzone+1)
    par_color_to_code_zone_subd(3) = modulo(icolo,msubd+1)

  end function par_color_to_code_zone_subd

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Check if a partition IPART is in color ICOLO
  !> @details IPART should be given in the world
  !>                                                       ICOLO
  !>                                           +---+---+---+---+---+---+
  !>          PAR_CPU_TO_COLOR(IPART) % L(:) = | 1 | 4 | 7 | 8 | 9 |10 |
  !>                                           +---+---+---+---+---+---+
  !>                                           KCOLO =>
  !>
  !----------------------------------------------------------------------

  function par_part_in_color(ipart,icolo)
    integer(ip), intent(in) :: ipart
    integer(ip), intent(in) :: icolo
    logical(lg)             :: par_part_in_color
    integer(ip)             :: kcolo,jcolo,ksize

    par_part_in_color = .false.
    if( associated(PAR_CPU_TO_COLOR(ipart) % l) ) then
       jcolo = 1
       kcolo = PAR_CPU_TO_COLOR(ipart) % l(jcolo)
       ksize = size(PAR_CPU_TO_COLOR(ipart) % l)
       do while( jcolo <= ksize .and. kcolo < icolo )
          kcolo = PAR_CPU_TO_COLOR(ipart) % l(jcolo)
          jcolo = jcolo + 1
       end do
       if( kcolo == icolo ) par_part_in_color = .true.
    end if

  end function par_part_in_color

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    06/03/2014
  !> @brief   Number of partitions for nodes
  !> @details Gives the number of neighbors for a given node or all
  !>          the nodes. This includes my partition. In next example,
  !>          the middle node (o) has 4 partitions.
  !>
  !>          x----------x----------x
  !>          |          |          |
  !>          |   CPU1   |   CPU2   |
  !>          |          |          |
  !>          x----------o----------x
  !>          |          |          |
  !>          |   CPU3   |   CPU4   |
  !>          |          |          |
  !>          x----------x----------x
  !
  !----------------------------------------------------------------------

  subroutine PAR_NODE_NUMBER_OF_PARTITIONS(npoin_loc,lneig,ipoin)
    integer(ip), intent(in)           :: npoin_loc
    integer(ip), intent(out)          :: lneig(*)
    integer(ip), intent(in), optional :: ipoin
    integer(ip)                       :: kpoin
    integer(ip)                       :: jj

    if( present(ipoin) ) then
       lneig(1) = 1
       do jj = 1,commd % bound_dim
          kpoin = commd % bound_perm(jj)
          if( kpoin == ipoin ) lneig(ipoin) = lneig(ipoin) + 1
       end do
    else
       do kpoin = 1,npoin_loc
          lneig(kpoin) = 1
       end do
       do jj = 1,commd % bound_dim
          kpoin = commd % bound_perm(jj)
          lneig(kpoin) = lneig(kpoin) + 1
       end do
    end if

  end subroutine PAR_NODE_NUMBER_OF_PARTITIONS

  function NODE_IN_NEIGHBOR_P(ipoin,ineig,commu)
    implicit none
    integer(ip),                  intent(in) :: ipoin     !< Node
    integer(ip),                  intent(in) :: ineig     !< Neighbor
    type(comm_data_par), pointer, intent(in) :: commu(:)  !< Generic communication array
    logical(lg)                              :: NODE_IN_NEIGHBOR_P
    integer(ip)                              :: ii

    NODE_IN_NEIGHBOR_P = .false.
    loop_nodes: do ii = commu(1) % bound_size(ineig),commu(1) % bound_size(ineig+1)-1
       if( ipoin == commu(1) % bound_perm(ii) )then
          NODE_IN_NEIGHBOR_P = .true.
          exit loop_nodes
       end if
    end do loop_nodes

  end function NODE_IN_NEIGHBOR_P

  function NODE_IN_NEIGHBOR_s(ipoin,ineig,commu)
    implicit none
    integer(ip),         intent(in) :: ipoin  !< Node
    integer(ip),         intent(in) :: ineig  !< Neighbor
    type(comm_data_par), intent(in) :: commu  !< Generic communication array
    logical(lg)                     :: NODE_IN_NEIGHBOR_s
    integer(ip)                     :: ii

    NODE_IN_NEIGHBOR_s = .false.
    loop_nodes: do ii = commu % bound_size(ineig),commu % bound_size(ineig+1)-1
       if( ipoin == commu % bound_perm(ii) )then
          NODE_IN_NEIGHBOR_s = .true.
          exit loop_nodes
       end if
    end do loop_nodes

  end function NODE_IN_NEIGHBOR_s

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    08/02/2019
  !> @brief   Transpose communication array
  !> @details Transpose communication array
  !
  !----------------------------------------------------------------------

  subroutine PAR_TRANSPOSE_COMMUNICATION_ARRAY(COMM,memor_opt)

    type(comm_data_par), intent(inout)           :: COMM          !< Input communicator
    integer(8),          intent(inout), optional :: memor_opt(2)  !< Memory counter
    integer(8)                                   :: memor(2)
    integer(ip)                                  :: lsend_dim_copy
    integer(ip),  pointer                        :: lsend_size_copy(:)
    integer(ip),  pointer                        :: lsend_perm_copy(:)

    nullify(lsend_size_copy)
    nullify(lsend_perm_copy)

    if( present(memor_opt) ) then
       memor = memor_opt
    else
       memor = 0_8
    end if
    !
    ! SEND => SEND_COPY
    !
    lsend_dim_copy = COMM % lsend_dim
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_TRANSPOSE_COMMUNICATION_ARRAY',COMM % lsend_size,lsend_size_copy,COPY_NAME='LSEND_SIZE_COPY')
    call memory_copy(memor,'COMM % LSEND_PERM','PAR_TRANSPOSE_COMMUNICATION_ARRAY',COMM % lsend_perm,lsend_perm_copy,COPY_NAME='LSEND_PERM_COPY')
    !
    ! RECV => SEND
    !
    COMM % lsend_dim = COMM % lrecv_dim
    call memory_copy(memor,'COMM % LRECV_SIZE','PAR_TRANSPOSE_COMMUNICATION_ARRAY',COMM % lrecv_size,COMM % lsend_size,COPY_NAME='COMM % LSEND_SIZE')
    call memory_copy(memor,'COMM % LRECV_PERM','PAR_TRANSPOSE_COMMUNICATION_ARRAY',COMM % lrecv_perm,COMM % lsend_perm,COPY_NAME='COMM % LSEND_PERM')
    !
    ! SEND_COPY => RECV
    !
    COMM % lrecv_dim = lsend_dim_copy
    call memory_copy(memor,'LSEND_SIZE_COPY','PAR_TRANSPOSE_COMMUNICATION_ARRAY',lsend_size_copy,COMM % lrecv_size,COPY_NAME='COMM % LRECV_SIZE')
    call memory_copy(memor,'LSEND_PERM_COPY','PAR_TRANSPOSE_COMMUNICATION_ARRAY',lsend_perm_copy,COMM % lrecv_perm,COPY_NAME='COMM % LRECV_PERM')

    if( present(memor_opt) ) memor_opt = memor

  end subroutine PAR_TRANSPOSE_COMMUNICATION_ARRAY
 
end module mod_parall
!
!> @}
!-----------------------------------------------------------------------
