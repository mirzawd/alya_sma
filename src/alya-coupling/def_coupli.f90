!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Coupling
!> @{
!> @name    Coupling arrays
!> @file    def_coupli.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!> @brief   ToolBox for coupling
!> @details ToolBox for coupling
!>          \verbatim
!>
!>          o----o----o----o----o----o----o----o Source
!>
!>          x-----------x----------x-----------x Target
!>        
!>          x ................... Wet nodes
!>          x-----x ............. Wet boundaries
!>
!>          NBOUN_WET ........... Number of wet boundaries
!>          LBOUN_WET(:) ........ List of wet boundaries
!>
!>          =>
!>
!>          NPOIN_WET ........... Number of wet nodes
!>          LPOIN_WET(:) ........ List of wet nodes
!>          
!>          NUMBER_WET_POINTS ... Number of wet points (nodes or Gauss points)
!>          COORD_WET(:,:) ...... Coordinates of wet points
!>          
!>          VMASB_WET(:) ........ Boundary mass matrix of wet nodes
!>          WEIGHT_WET(:) ....... Weight of values of wet nodes
!>
!>          \endverbatim
!>
!> @{
!------------------------------------------------------------------------

module def_coupli
  use def_kintyp_basic,         only : ip,rp,lg
  use def_kintyp_comm,          only : comm_data_par
  use def_spmat,                only : spmat
  use mod_kdtree,               only : typ_kdtree
  use def_master,               only : AT_BEGINNING
  use def_domain,               only : nelty
  implicit none  

  !----------------------------------------------------------------------
  !
  ! Types and parameters
  !
  !----------------------------------------------------------------------

  integer(ip),        parameter     :: ON_SET                     = 1    ! Where
  integer(ip),        parameter     :: ON_FIELD                   = 2    ! Where
  integer(ip),        parameter     :: ON_CODE                    = 3    ! Where
  integer(ip),        parameter     :: ON_WHOLE_MESH              = 4    ! Where
  integer(ip),        parameter     :: ON_CHIMERA_MESH            = 5    ! Where
  integer(ip),        parameter     :: ON_MIRROR                  = 6    ! Where
  integer(ip),        parameter     :: ON_IMMERSED_MESH           = 7    ! Where
  integer(ip),        parameter     :: ON_FLOATING_POINTS         = 8    ! Where
  integer(ip),        parameter     :: ON_EMBEDDED_MESH           = 9    ! Where
  integer(ip),        parameter     :: ELEMENT_INTERPOLATION      = 1    ! Type
  integer(ip),        parameter     :: NEAREST_BOUNDARY_NODE      = 2    ! Type
  integer(ip),        parameter     :: BOUNDARY_INTERPOLATION     = 3    ! Type
  integer(ip),        parameter     :: PROJECTION                 = 4    ! Type
  integer(ip),        parameter     :: NEAREST_ELEMENT_NODE       = 5    ! Type
  integer(ip),        parameter     :: BOUNDARY_VECTOR_PROJECTION = 6    ! Type
  integer(ip),        parameter     :: STRESS_PROJECTION          = 7    ! Type
  integer(ip),        parameter     :: TRANSPOSE_MIRROR           = 8    ! Type
  integer(ip),        parameter     :: GLOBAL_NUMBERING           = 9    ! Type
  integer(ip),        parameter     :: SAME_COORDINATE            = 10   ! Type
  integer(ip),        parameter     :: UNKNOWN                    = 1    ! What
  integer(ip),        parameter     :: RESIDUAL                   = 2    ! What
  integer(ip),        parameter     :: DIRICHLET_EXPLICIT         = 3    ! What
  integer(ip),        parameter     :: DIRICHLET_IMPLICIT         = 4    ! What
  integer(ip),        parameter     :: INTRINSIC                  = 5    ! What
  integer(ip),        parameter     :: SCALAR                     = 6    ! What
  integer(ip),        parameter     :: RELAXATION_SCHEME          = 1    ! Scheme (relaxation)
  integer(ip),        parameter     :: AITKEN_SCHEME              = 2    ! Scheme (relaxation)
  integer(ip),        parameter     :: BROYDEN_SCHEME             = 3    ! Scheme (relaxation)
  integer(ip),        parameter     :: IQNLS_SCHEME               = 4    ! Scheme (relaxation)
  integer(ip),        parameter     :: BETWEEN_ZONES              = 1    ! Kind
  integer(ip),        parameter     :: BETWEEN_SUBDOMAINS         = 2    ! Kind
  integer(ip),        parameter     :: I_SOURCE                   = 1    ! Source/target
  integer(ip),        parameter     :: I_TARGET                   = 2    ! Source/target
  integer(ip),        parameter     :: FIXED_UNKNOWN              = 1111 ! Fixity
  integer(ip),        parameter     :: INTERFACE_MASS             = 1    ! Conservation
  integer(ip),        parameter     :: GLOBAL_MASS                = 2    ! Conservation
  integer(ip),        parameter     :: VALUES_ON_NODES            = 1    ! Type
  integer(ip),        parameter     :: VALUES_ON_ELEMENTS         = 2    ! Type
  integer(ip),        parameter     :: VALUES_ON_BOUNDARIES       = 3    ! Type
  
  integer(ip),        parameter     :: IMPOSE_ZERO                = 0    ! Lost wet point strategy
  integer(ip),        parameter     :: DONT_DO_ANYTHING           = 1    ! Lost wet point strategy
  integer(ip),        parameter     :: STOP_ALYA                  = 2    ! Lost wet point strategy
  integer(ip),        parameter     :: STOP_ALYA_WITH_WARNINGS    = 3    ! Lost wet point strategy

  integer(ip),        parameter     :: NODE_WET_POINT             = 0    ! Wet point is a node
  integer(ip),        parameter     :: FLOATING_WET_POINT         = 1    ! Wt point is not a node
  
  integer(ip),        parameter     :: NODE_TARGET_ENTITY         = 0    ! Wet point is a node
  integer(ip),        parameter     :: FLOATING_TARGET_ENTITY     = 1    ! Wt point is not a node
  !
  ! Driver
  !
  integer(ip),        parameter     :: max_block_cou = 10
  integer(ip),        parameter     :: max_coupl_cou = 10
  integer(ip)                       :: coupling_driver_couplings(max_coupl_cou,max_block_cou)  
  integer(ip)                       :: coupling_driver_number_couplings(max_block_cou)  
  integer(ip)                       :: coupling_driver_max_iteration(max_block_cou)       
  integer(ip)                       :: coupling_driver_iteration(max_block_cou)           
  real(rp)                          :: coupling_driver_tolerance(max_block_cou)
  integer(ip)                       :: kfl_gozon
  !
  ! Partition
  ! 
  integer(ip)                       :: kfl_graph_cou = 0_ip
  !
  ! Useful boundary: could be different from original one when
  ! using Chimera. In this case the fringe boundaries are only
  ! modified here
  !
  integer(ip)                       :: number_of_holes                ! Number of holes
  integer(ip)                       :: nboun_cou                      ! Number of boundaries
  integer(ip),      pointer         :: lnodb_cou(:,:)                 ! Boundary connectivity
  integer(ip),      pointer         :: ltypb_cou(:)                   ! Boundary type
  integer(ip),      pointer         :: lboch_cou(:)                   ! Boundary characteristics
  integer(ip),      pointer         :: lnnob_cou(:)                   ! Number of boundary nodes
  integer(ip),      pointer         :: lboel_cou(:,:)                 ! Boundary/element connectivity
  integer(ip),      pointer         :: lelbo_cou(:)                   ! Boundary/element connectivity
  integer(ip)                       :: ncoup_implicit_d               ! Number of implicit coupling
  integer(ip),      pointer         :: lcoup_implicit_d(:)            ! List of implicit couplings
  integer(ip)                       :: ncoup_implicit_n               ! Number of implicit coupling
  integer(ip),      pointer         :: lcoup_implicit_n(:)            ! List of implicit couplings
  integer(ip)                       :: ncoup_implicit                 ! Couplings involving implicit coupling
  real(rp),         pointer         :: mask_cou(:)                    ! Mask for Dirichlet nodes

  ! Cloud of points
  integer(ip)                       :: ngaus_elem(nelty)              ! Number of Gauss points for cloud of points

  ! Immersed method
  logical(lg)                       :: kfl_immer = .false.            ! Any one of couplings is of immersed-type

  ! Immersed boundary method for deformable solid bodies
  logical(lg)                       :: kfl_dimbou = .false.           ! Flag takes into account 2 couplings

  ! EFECT
  logical(lg)                       :: kfl_efect = .false.            ! Embedded Finite Element Coupling Technique

  !
  ! When adding a variables, modify the following subroutine
  ! INITIALIZE_COUPLING_STRUCTURE
  ! DEALLOCATE_COUPLING_STRUCTURE
  !
  !
  ! Projection type to be used on target
  !
  type typ_proje
     integer(ip)          :: pgaub
     integer(ip), pointer :: permu(:)
     real(rp),    pointer :: shapb(:,:)
     real(rp),    pointer :: gbsur(:)
  end type typ_proje

  type typ_coupling_wet
     integer(ip)                    :: number_wet_points              ! Number of wet points (nodes or Gauss points)
     integer(ip)                    :: npoin_wet                      ! Number of wet nodes
     integer(ip)                    :: nboun_wet                      ! Number of wet boundaries
     integer(ip)                    :: nolap_wet                      ! Number of wet nodes in overlapping nodes list
     integer(ip)                    :: point_type                     ! Wet point type
     integer(ip),      pointer      :: lboun_wet(:)                   ! List of wet boundaries
     logical(lg)                    :: kfl_get_fringe                 ! Flag for computing fringe wetnodes
     integer(ip),      pointer      :: lpoin_wet(:)                   ! List of wet nodes
     integer(ip),      pointer      :: lsubd(:)                       ! List of wet points subdomains
     integer(ip),      pointer      :: lolap_wet(:)                   ! List of wet nodes in overlapping nodes list
     integer(ip),      pointer      :: kfl_fringe_wetnodes(:)         ! Flag for fringe wet nodes (=1) interior wet nodes (-1) or dry nodes (0)
     real(rp),         pointer      :: coord_wet(:,:)                 ! Weight nodes coordinates
     real(rp),         pointer      :: coord_wet_ini(:,:)             ! Weight nodes coordinates for moving subdomains
     real(rp),         pointer      :: weight_wet(:)                  ! Weight for force like
     real(rp),         pointer      :: weight_wet_imp(:)              ! Weight for force like used for implicit coupling
     !
     ! Unused yet, DEFINED COU_WET_INTERFACE_MASS_MATRIX
     !
     real(rp),         pointer      :: vmasb_wet(:)                   ! Boundary mass matrix
     !
     ! Unused yet
     !
     real(rp),         pointer      :: mass_ratio_wet(:)              ! Mass matrix around wet nodes
     ! 
     ! Memory allocation or values initialization in COU_CREATE_A_SUBMESH unused yet
     !
     type(typ_proje),  pointer      :: proje_target(:)
     !
     ! Previous wet-nodes
     !
     integer(ip)                    :: npoin_wet_prev                 ! Number of wet nodes at previous time step
     integer(ip),      pointer      :: lpoin_wet_prev(:)              ! List of wet nodes at previous time step
     integer(ip)                    :: npoin_fresh                    ! Number of fresh nodes (convert from wet to regular)
     integer(ip),      pointer      :: lpoin_fresh(:)                 ! List of fresh nodes
   contains
     procedure,   pass              :: init => init_wet               ! Initialize wet geometry     
  end type typ_coupling_wet
  
  type typ_coupling_geometry
     !
     ! Information used as a source
     !
     !
     ! Unused yet
     !
     integer(ip)                    :: ndime                          ! Problem dimension
     !
     ! Memory allocation or integer definition in COU_INIT_INTERPOLATE_POINTS_VALUES
     !
     integer(ip)                    :: nelem_source                   ! Number of source elements
     integer(ip)                    :: nboun_source                   ! Number of source boundaries
     integer(ip)                    :: npoin_source                   ! Number of source nodes
     integer(ip),      pointer      :: lelem_source(:)                ! List of source elements
     integer(ip),      pointer      :: lboun_source(:)                ! List of source boundaries
     integer(ip),      pointer      :: lpoin_source(:)                ! List of source nodes
     integer(ip),      pointer      :: status(:)                      ! Results of the interpolation
     real(rp),         pointer      :: shapf(:,:)                     ! Shape functions for interpolations
     !
     ! Definition in COU_INITIALIZE_COUPLING
     !
     type(typ_kdtree)               :: kdtree                         ! Kd-tree
     !
     ! Definition in COU_SOURCE_INTERFACE_MASS_MATRIX
     !
     real(rp),         pointer      :: vmasb(:)                       ! Mass matrix used for projection
     !
     ! Wet information used as a target. 
     !
     ! Memory allocation or values initialization 
     ! in cou_define_wet_geometry
     !
     integer(ip),      pointer      :: sched_perm(:)                  ! Scheduling permutation of status
   contains
     procedure,   pass              :: init => init_geometry          ! Initialize geometry     
  end type typ_coupling_geometry

  type typ_color_coupling
     integer(ip)                    :: number
     integer(ip)                    :: itype
     integer(ip)                    :: code_source
     integer(ip)                    :: code_target
     integer(ip)                    :: zone_source
     integer(ip)                    :: zone_target
     integer(ip)                    :: color_source
     integer(ip)                    :: color_target
     integer(ip)                    :: module_source
     integer(ip)                    :: module_target
     integer(ip)                    :: subdomain_source
     integer(ip)                    :: subdomain_target
     integer(ip)                    :: where_type
     integer(ip)                    :: where_number
     integer(ip)                    :: where_type_source
     integer(ip)                    :: where_number_source
     integer(ip)                    :: target_entity                  ! Where coupling is eventuall applied, different from wet % point_type
     integer(ip)                    :: kfl_par_transm
     integer(ip)                    :: kfl_check_exha
     integer(ip)                    :: kfl_toda_costa
     integer(ip)                    :: what
     integer(ip)                    :: scheme
     integer(ip)                    :: itera
     integer(ip)                    :: conservation
     integer(ip)                    :: overlap                        ! Number of overlapping elements (0 for disjoint)
     integer(ip)                    :: ngaus                          ! Gauss points strategy when using projection
     integer(ip)                    :: kfl_multi_source               ! Enable multi-source communication
     integer(ip)                    :: kfl_lost_wet_points            ! 0=stop, 1=discard
     integer(ip)                    :: kfl_symmetry                   ! Coupling symmetry (to save time)
     integer(ip)                    :: kfl_source_value               ! Where source value is
     integer(ip)                    :: kfl_fixity                     ! 0=free, 1=coupling imposes fixity
     
     integer(ip)                    :: kind                           ! Kind of coupling: subdomain/zone type
     integer(ip)                    :: mirror_coupling                ! Number of the mirror coupling
     integer(ip)                    :: task_compute_and_send
     integer(ip)                    :: when_compute_and_send
     integer(ip)                    :: task_recv_and_assemble
     integer(ip)                    :: when_recv_and_assemble
     integer(ip)                    :: frequ_send                     ! Frequency of sends in coupling
     integer(ip)                    :: frequ_recv                     ! Frequency of recvs in coupling
     integer(ip)                    :: when_update_wet_geometry
     integer(ip)                    :: when_update_coupling
     integer(ip)                    :: temporal_predictor             ! 0: No temporal prediction will be done, 1: temporal prediction activated
     integer(ip)                    :: temporal_predictor_order       ! -1: No temporal prediction, 0: zeroth order, 1: first order ...
     character(5)                   :: variable
     real(rp)                       :: relax
     real(rp)                       :: aitken                  ! Aitken relaxation factor calculated in cou_update_values
     real(rp)                       :: resid(2)                ! Residual
     real(rp)                       :: cputim(20)              ! cputim
     real(rp),    pointer           :: values(:,:,:)           ! Relaxed exchanged quantity, allocated in mod_couplings
     real(rp),    pointer           :: values_frequ(:,:,:)     ! Values for the exchanged quantities in frequency mode
     real(rp),    pointer           :: values_predicted(:,:)   ! Unrelaxed exchanged quantity, allocated in mod_couplings
     real(rp),    pointer           :: values_converged(:,:,:) ! Converged exchanged quantity (temporal predictor scheme) allocated in cou_check_convergence
     real(rp),    pointer           :: jacobian_inverse(:,:,:) ! inverse of the jacobian for the broyden problem (npoin_wet x npoin_wet)
     real(rp),    pointer           :: dincr_predicted(:)      ! G_ij * v_j for the broyden problem
     integer(ip)                    :: ranku_iqnls             ! Number of past iterations used in INQLS matrix free
     integer(ip)                    :: history_iqnls           ! Number of past TIME iterations used
     real(rp)                       :: scaling_iqnls           ! Scaling for the values in the CIQN algorithm
     real(rp)                       :: efilter_iqnls           ! Epsilon for the CIQN filter
     real(rp),    pointer           :: residues_iqnls(:, :)    ! Residues for the IQNLS for different time steps
     real(rp),    pointer           :: relaxed_iqnls(:, :)     ! Relaxed solutions needed for IQNLS
     real(rp),    pointer           :: unrelaxed_iqnls(:, :)   ! Unrelaxed solutions needed for IQNLS
     real(rp),    pointer           :: valincr_iqnls(:, :)     ! Increment of the values for that iteration 
     real(rp),    pointer           :: residincr_iqnls(:, :)   ! Increment of the residual for that iteration 
     real(rp),    pointer           :: residincr_history_iqnls(:, :, :)   ! Histories of the residue increments
     real(rp),    pointer           :: valincr_history_iqnls(:, :, :)   ! Histories of the value increments
     real(rp),    pointer           :: V_current_history_iqnls(:, :)   ! V current + histories
     real(rp),    pointer           :: W_current_history_iqnls(:, :)   ! V current + histories
     integer(ip), pointer           :: history_tracking_iqnls(:)   ! Increment of the residual for that iteration 
     type(typ_coupling_wet)         :: wet                     ! Wet geometry
     type(typ_coupling_geometry)    :: geome                   ! Geometry
     type(comm_data_par)            :: commd                   ! Communicator
     type(spmat), pointer           :: ltransmat_target(:)     ! Transmission matrix
   contains
     procedure,   pass              :: init => init_coupling   ! Initialize coupling
  end type typ_color_coupling

  !----------------------------------------------------------------------
  !
  ! Read data
  !
  !----------------------------------------------------------------------

  integer(ip)                       :: mcoup                          ! Number of couplings
  integer(ip)                       :: coudt = -1_ip
  integer(ip)                       :: kfl_lost_wet_point_cou         ! Lost wet poiunt strategy
  integer(ip)                       :: kfl_absolute_cou               ! Strategy for absolute tolerance
  real(rp)                          :: toler_absolute_cou             ! Absolute tolerance for partition communciation
  real(rp)                          :: toler_relative_cou             ! Relative tolerance for partition communciation
  integer(ip)                       :: kfl_timco_cou                  ! Time coupling strategy
  
  !----------------------------------------------------------------------
  !
  ! Definitions
  !
  !----------------------------------------------------------------------

  type(typ_color_coupling), pointer :: coupling_type(:)               ! Coupling array
  integer(ip)                       :: lun_coupl_dat                  ! Coupling data file
  integer(ip)                       :: lun_coupl_res                  ! Coupling result file
  integer(ip)                       :: lun_coupl_cvg                  ! Coupling result file
  integer(8)                        :: memor_cou(2)                   ! Memory
  type(typ_kdtree)                  :: kdtree_typ                     ! SKD-Tree structure
  type(comm_data_par),      pointer :: COU_COMM_COUPLING_ARRAY(:)     ! Sub communicator for projection type coupling

  integer(ip)                       :: nscal_cou                      ! Number of scalar couplings
  real(rp),                 pointer :: scala_cou(:,:)                 ! Scalar values to broadcast
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-04
  !> @brief   Coupling initialization
  !> @details Initialization of coupling class
  !> 
  !-----------------------------------------------------------------------

  subroutine init_coupling(coupling,WET,INPUT)

    class(typ_color_coupling), intent(inout)        :: coupling
    logical(lg),               intent(in), optional :: WET
    logical(lg),               intent(in), optional :: INPUT
    logical(lg)                                     :: ifwet
    logical(lg)                                     :: ifinput

    ifwet = .true.
    if( present(WET) ) ifwet = WET

    ifinput = .true.
    if( present(INPUT) ) ifinput = INPUT

    if( ifinput ) then
       coupling % number                   =  0 
       coupling % itype                    =  0 
       coupling % code_source              =  1 
       coupling % code_target              =  1 
       coupling % zone_source              =  0 
       coupling % zone_target              =  0 
       coupling % color_source             =  0 
       coupling % color_target             =  0 
       coupling % module_source            =  0 
       coupling % module_target            =  0 
       coupling % subdomain_source         =  0 
       coupling % subdomain_target         =  0 
       coupling % where_type               =  0 
       coupling % where_number             =  0
       coupling % where_type_source        =  0 
       coupling % where_number_source      =  0
       coupling % target_entity            =  NODE_TARGET_ENTITY 
       coupling % kfl_toda_costa           =  1 
       coupling % kfl_par_transm           =  0 
       coupling % kfl_check_exha           =  0 
       coupling % what                     =  UNKNOWN 
       coupling % scheme                   =  RELAXATION_SCHEME     
       coupling % itera                    =  0 
       coupling % conservation             =  0 
       coupling % overlap                  =  0            
       coupling % ngaus                    =  0
       coupling % kfl_symmetry             =  0
       coupling % kfl_source_value         =  VALUES_ON_NODES
       coupling % kfl_multi_source         =  0
       coupling % kfl_lost_wet_points      =  0
       coupling % kind                     =  0      
       coupling % mirror_coupling          =  0      
       coupling % task_compute_and_send    = -100_ip 
       coupling % when_compute_and_send    = -100_ip 
       coupling % task_recv_and_assemble   = -100_ip 
       coupling % when_recv_and_assemble   = -100_ip 
       coupling % frequ_send               =  0      
       coupling % frequ_recv               =  0      
       coupling % when_update_wet_geometry =  AT_BEGINNING
       coupling % when_update_coupling     =  AT_BEGINNING
       coupling % temporal_predictor       =  0      
       coupling % temporal_predictor_order = -1_ip      
       coupling % variable                 =  '     '
       coupling % relax                    =  1.0_rp
       coupling % resid                    =  0.0_rp
       coupling % cputim                   =  0.0_rp
       coupling % aitken                   =  1.0_rp
       coupling % ranku_iqnls              = -1_ip
       coupling % history_iqnls            =  1_ip
       coupling % scaling_iqnls            = -1.0_rp
       coupling % efilter_iqnls            = -1.0_rp
       coupling % kfl_fixity               =  1
    end if

    nullify(coupling % values)           
    nullify(coupling % values_frequ)     
    nullify(coupling % values_predicted)   
    nullify(coupling % values_converged) 
    nullify(coupling % jacobian_inverse) 
    nullify(coupling % dincr_predicted)      
    nullify(coupling % residues_iqnls)    
    nullify(coupling % relaxed_iqnls)     
    nullify(coupling % unrelaxed_iqnls)   
    nullify(coupling % valincr_iqnls)     
    nullify(coupling % residincr_iqnls)
    nullify(coupling % valincr_history_iqnls)
    nullify(coupling % residincr_history_iqnls)
    nullify(coupling % V_current_history_iqnls)
    nullify(coupling % W_current_history_iqnls)
    nullify(coupling % history_tracking_iqnls)

    nullify(coupling % ltransmat_target)

    call coupling % geome % init()
    if( ifwet ) call coupling % wet % init()
    call coupling % commd % init(PAR_COMM_OPT=.false.,COMM_NAME='COUPLING % COMMD')

  end subroutine init_coupling

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-04
  !> @brief   Geometry initialization
  !> @details Initialization of geometry class
  !> 
  !-----------------------------------------------------------------------
  
  subroutine init_geometry(geome)
    class(typ_coupling_geometry), intent(inout) :: geome

    geome % ndime             = 0     
    geome % nelem_source      = 0      
    geome % nboun_source      = 0         
    geome % npoin_source      = 0 
    
    nullify(geome % lelem_source)            
    nullify(geome % lboun_source)            
    nullify(geome % lpoin_source)            
    nullify(geome % status)                  
    nullify(geome % shapf)               
    nullify(geome % vmasb)
    nullify(geome % sched_perm) 
    
    call geome % kdtree % init()
    
   end subroutine init_geometry

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-04
  !> @brief   Wet geometry initialization
  !> @details Initialization wet of geometry class
  !> 
  !-----------------------------------------------------------------------
  
   subroutine init_wet(wet)

     class(typ_coupling_wet), intent(inout) :: wet

     wet % number_wet_points = 0  
     wet % npoin_wet         = 0   
     wet % nboun_wet         = 0
     wet % point_type        = NODE_WET_POINT
     wet % kfl_get_fringe    = .false.

     nullify(wet % lsubd)               
     nullify(wet % lboun_wet)               
     nullify(wet % lpoin_wet)               
     nullify(wet % lolap_wet)
     nullify(wet % kfl_fringe_wetnodes)
     nullify(wet % coord_wet)             
     nullify(wet % coord_wet_ini)             
     nullify(wet % weight_wet)              
     nullify(wet % weight_wet_imp)              
     nullify(wet % vmasb_wet)               
     nullify(wet % mass_ratio_wet)          
     nullify(wet % proje_target)

     ! Initialization of attributes which save historic values is conditional
     if (.not. associated(wet % lpoin_wet_prev) ) then
        nullify(wet % lpoin_wet_prev)
        wet % npoin_wet_prev = 0
     end if 
     if (.not. associated(wet % lpoin_wet_prev) ) then
        nullify(wet % lpoin_fresh)
        wet % npoin_fresh = 0
     end if 
  
  end subroutine init_wet
  
end module def_coupli
!> @} 
!-----------------------------------------------------------------------
