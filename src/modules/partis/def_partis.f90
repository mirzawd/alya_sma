!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module def_partis
  !------------------------------------------------------------------------
  !****f* Partis/def_partis
  ! NAME 
  !    def_partis
  ! DESCRIPTION
  !    Heading for the Partis routines
  ! USES
  ! USED BY
  !    Almost all
  !***
  !------------------------------------------------------------------------
  use def_kintyp
  use def_inpout,      only : maxwp
  use def_partis_type, only : latyp
  use def_partis_type, only : mlapr
  use def_partis_type, only : mtyla
  use def_partis_type, only : parttyp
  use def_partis_type, only : lagrtyp
  use def_partis_type, only : typ_injector
  use mod_physics,     only : liquid_state
  use mod_interp_tab,  only : typ_lookup_framework
  use mod_interp_tab,  only : typ_tab_coord, typ_lookup_table 

  USE, INTRINSIC :: ISO_C_BINDING

  !------------------------------------------------------------------------
  !
  ! Lagrangian particle, needed by Master for coupling
  !
  !------------------------------------------------------------------------

  integer(ip)               :: mlagr             ! Lagrangian particles: Max. number of particle in each subdomain
  integer(ip), parameter    :: pts_advec_narrays=3 !number of arrays of advec(:,:,:) to use, partis uses only advec(:,:,1) and advec(:,:,3). Needed when nastin is off and timsteps are loaded
  
 
  integer(ip)               :: nlagr

  !------------------------------------------------------------------------
  !
  ! Parameters
  !
  !------------------------------------------------------------------------

  integer(ip), parameter   ::      &
       PTS_OUTFLOW_CONDITION  = 0, &
       PTS_WALL_CONDITION     = 1, &
       PTS_SLIP_CONDITION     = 3, &
       PTS_BOUNCING_CONDITION = 4
  
  integer(ip), parameter   ::      &
       PTS_PARTICLE_EXISTS      = -1,&
       PTS_PARTICLE_HITS_WALL   = -2,&
       PTS_PARTICLE_ZERO_TIME   = -3,&
       PTS_PARTICLE_OUTFLOW     = -4,&
       PTS_PARTICLE_JUST_INJ    = -5,& ! particle has just been injected (beginning of timestep)
       PTS_PARTICLE_EVAPORATED  = -6,& ! particle fully evaporated
       PTS_PARTICLE_MOVING_MESH = -7   ! particle going out of a moving mesh

  integer(ip), parameter   :: &
       mpala     = 20,        & ! Max # parameters for lagrangian particles
       mpara_inj =  4,        & ! Max # parameters for lagrangian particles
       pts_minj  = 100,       & ! Max # injections
       mvarp_pts = 70           ! Max # of postprocess variables
  
  integer(ip), parameter   :: &
       lun_resul_pts = 1725,  & ! Result file unit
       lun_oudep_pts = 1728,  & ! Deposition map file unit
       lun_depsu_pts = 1729     ! Deposition surface file unit

  character(150)           :: &
       fil_resul_pts,         & ! Result file 
       fil_oudep_pts,         & ! Deposition map file
       fil_depsu_pts            ! Deposition surface file

  real(rp),    parameter   :: &
       zepts = epsilon(1.0_rp)   ! Zero

  !------------------------------------------------------------------------
  !
  ! Physical problem: read in pts_reaphy
  !
  !------------------------------------------------------------------------
  !--BEGIN REA GROUP

  real(rp)                   :: &
       dtmin_pts,               &    ! Minimum time step 
       dimin_pts,               &    ! Minimum distance to wall under which particle is deposited
       prthe_pts,               &    ! Thermodynamic pressure (if NASTIN not activated)
       mean_free_path_pts            ! Mean free path of medium
  integer(ip)                :: &
       kfl_momentum_pts,        &    ! Momentum transport  
       kfl_thermo_pts,          &    ! Thermodynamic transport
       kfl_momentum_sink_pts,   &    ! Compute momentum sink for coupling with fluid
       kfl_heat_sink_pts,       &    ! Heat sink
       kfl_mass_sink_pts,       &    ! Mass sink
       nclas_pts                     ! Number of classes in conce, for evaporating droplets
       
  
  !------------------------------------------------------------------------
  !
  ! Numerical problem: read in pts_reanut
  !
  !------------------------------------------------------------------------

  integer(ip)              :: &
       kfl_adapt_pts,         &      ! Adaptive time step
       kfl_usbin_pts,         &      ! Use element bin for element search
       kfl_order_pts,         &      ! Order of integration of AB
       kfl_walld_pts,         &      ! Activate walld calculation for corss wall test
       kfl_vesgs_pts,         &      ! Tracking of sgs velocity
       kfl_thermo_timsch_pts         ! Thermodynamic transport time scheme 
  real(rp)                 :: & 
       gamma_pts,             &      ! Lagrangian particles: Initial time of injection
       beta_pts,              &      ! Lagrangian particles: Time period of injection
       chale_pts,             &      ! Characteristic length for adaptive time step
       safet_pts,             &      ! Safety factor
       dtime_pts                     ! Particle time step
  
  !------------------------------------------------------------------------
  !
  ! Output and Postprocess: read in pts_reaous
  !
  !------------------------------------------------------------------------

  type, BIND(C) :: dbparams
     integer(C_INT) :: kfl_db_deep            !number of cube levels in the daba base
     character  :: kfl_db_url_conn(150)       !url database connection data
     integer(C_INT) :: port                   !port connection
     integer(C_INT) :: maxInserts             !number maximum of parall inserts
     integer(C_INT) :: numParticlesBatch      !number maximum of parall inserts    
  end type dbparams

  type (dbparams) dbSettings

  integer(ip)              :: &     
       kfl_posla_pts,         &              ! Lagrangian particles: Postprocess 
       kfl_oudep_pts,         &              ! Output of deposition map
       kfl_depos_surface_pts, &              ! Postprocess deposition surface
       kfl_oufre_pts,         &              ! Output frequency
       kfl_dbfre_pts,         &              ! Output DB frequency
       kfl_lost_moving_pts,   &              ! Exists particles lost by moving walls
       kfl_exacs_pts                         ! Exact solution

  character(100)           :: wordsdb(maxwp) !database load word

  !------------------------------------------------------------------------
  !
  ! Physical problem: read in pts_reabcs
  !
  !------------------------------------------------------------------------

  integer(ip)                             :: &  
       kfl_imax_injector_pts,                &        ! Highest injector index defined in input
       kfl_boundary_injection ,              &        ! Inject particles from boundary
       codbo_pts(pts_minj)                            ! Code of the boundary to inject particles
  real(rp)                                :: & 
       tinla_pts,                            &        ! Lagrangian particles: Initial time of injection
       tfila_pts,                            &        ! Lagrangian particles: Final time of injection
       tpela_pts                                      ! Lagrangian particles: Time period of injection
  integer(ip),        pointer             :: &
       kfl_fixbo_pts(:)                               ! Element boundary fixity     
  real(rp),           pointer             :: &
       bvnat_pts(:,:)                                 ! Natural bc values
  type(bc_bound),     pointer             :: &     
       tbcod_pts(:)                                   ! Boundary code type  

  type(typ_injector) :: injection_pts(pts_minj)

  !--END REA GROUP

  !------------------------------------------------------------------------
  !
  ! Others
  !
  !------------------------------------------------------------------------

  integer(ip)              :: &
       nlacc_pts,             &        ! Number of accumulated particles
       kfl_depos_pts,         &        ! Postprocess deposition map
       kfl_resid_pts,         &        ! Postprocess residence time
       ntyla_pts,             &        ! # Number of effective particle types (max type used)     
       number_types_pts,      &        ! # Number of particle types
       ninj_pts,              &        !
       kfl_slip_wall_pts,     &        ! If there are slip walls
       kfl_bouncing_wall_pts, &        ! If there are bouncing walls
       kfl_injec,             &        ! If particles have been injected
       kfl_max_out_exist_state_pts, &  ! Maximum existence flag for output
       kfl_min_out_exist_state_pts, &  ! Minimum existence flag for output
       kfl_max_prop_gwit_pts, &        ! Maximum property index used for geometrical witness points
       nlagr_local_pts,       &        ! Number of existing local particles
       nlagr_free_pts,        &        ! Number of free positions in type 
       nlagr_existing_pts,    &        ! # total number of existing particles in all subdomains _1
       nlagr_non_migrating_pts,&       ! # particles going from one subdomain to another        _2
       nlagr_going_out_pts,   &        ! # particles going out of the computational domain      _3
       nlagr_zero_time_pts,   &        ! # particles that disappear because of zero time step   _4
       nlagr_deposited_pts,   &        ! # particles that deposited without nboun               _5
       nlagr_hits_wall_pts,   &        ! # particles that deposited on wall
       nvarp_pts,             &        ! # Number of postprocess variables
       nvard_pts                       ! # Number of deposition variables
  real(rp)                 :: &   
       time_transport_pts              ! CPU time for transport
  integer(ip), pointer     :: &
       lboue_pts(:),          &        ! List of element touching boundaries
       kfl_fixno_walld_slip_pts(:,:),& ! Fixity wall distance to slip boundaries
       kfl_fixno_walld_bouncing_pts(:,:),& ! Fixity wall distance to bouncing boundaries
       permu_nlagr_pts(:)              ! Permutation
  real(rp),    pointer     :: &
       resid_pts(:,:),        &        ! Residence time
       depoe_pts(:,:),        &        ! Deposition map over elements
       depob_pts(:,:),        &        ! Deposition map over boundaries
       defor_pts(:,:),        &        ! Velocity deformation tensor
       hleng_pts(:),          &        ! Element characteristic length
       bouno_pts(:,:),        &        ! Boundary outwards normal
       walld_slip_pts(:),     &        ! Wall distance to slip boundaries
       walld_bouncing_pts(:), &        ! Wall distance to bouncing boundaries
       friction_pts(:),       &        ! Friction coefficient
       avg_mie_pts(:)                  ! Average Mie scattering equivalent 
  type(i1p),   pointer     :: &
       leleboun_pts(:)                 ! List of element boundaries
  type(latyp), pointer     :: &
       lagrtyp_tmp(:)                  ! Copy of particle type
  integer(ip)              :: &
       particles_sent,        &        ! Number of sent particles
       particles_recv,        &        ! Number of received particles
       comm_loops_pts                  ! Number of communicaiton loops
  integer(ip)              :: &
       kfl_rstar_pts                   ! Independent restart of particles
  logical(lg)              :: &
       if_moving_mesh_pts              ! If moving mesh
  !
  ! Postprocess
  !
  logical(lg)              :: &
       postprocess_var_pts(mvarp_pts)  ! Postprocess variables
  character(5)             :: &
       postprocess_name_pts(mvarp_pts) ! Postprocess variable names
  logical(lg)              :: &
       deposition_var_pts(mvarp_pts)   ! Deposition variables
  logical(lg)              :: &
       migrated_variables_pts(mvarp_pts) ! Variables to migrate
  integer(ip)              :: &
       number_migrated_variables_pts   ! Variables to migrate
  integer(ip),  pointer    :: &
       postprocess_list_pts(:)         ! List of postprocess variables
  integer(ip),  pointer    :: &
       deposition_list_pts(:)          ! List of deposition variables

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-10-03
  !> @brief   Find a variable ID
  !> @details According to a string find the variable ID
  !> 
  !-----------------------------------------------------------------------

  integer(ip) function pts_name_to_variable_number(wname)

    implicit none
    character(*), intent(in) :: wname
    integer(ip)              :: ii

    pts_name_to_variable_number = 0
    do ii = 1,mvarp_pts
       if( trim(wname) == trim(postprocess_name_pts(ii)) ) then
          pts_name_to_variable_number = ii
          return
       end if
    end do
    
  end function pts_name_to_variable_number
  
end module def_partis
