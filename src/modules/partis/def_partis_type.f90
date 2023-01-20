!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module def_partis_type

  use def_kintyp_basic, only : ip,rp
  use def_master,       only : dtime
  use mod_physics,      only : liquid_state
  use mod_interp_tab,   only : typ_lookup_framework
  use mod_interp_tab,   only : typ_tab_coord, typ_lookup_table
  implicit none
  private

  integer(ip), parameter    :: mlapr = 10        ! Maximum number of properties in Lagrangian particles
  integer(ip), parameter    :: mtyla = 100       ! Max # of lagrangian particle types
  !
  ! Particle type
  !
  type typlatyp
     integer(ip)            :: kfl_exist         ! If type exists
     integer(ip)            :: kfl_modla         ! Transport model
     integer(ip)            :: kfl_therm         ! Thermodynamic model
     integer(ip)            :: kfl_heattr_corr   ! Correction of heat transfer due to evaporation
     integer(ip)            :: kfl_mass_pot      ! Mass transfer potential formulation
     integer(ip)            :: kfl_grafo         ! Gravity
     integer(ip)            :: kfl_buofo         ! Buoyancy
     integer(ip)            :: kfl_drafo         ! Drag force model
     integer(ip)            :: kfl_brown         ! Brownian force
     integer(ip)            :: kfl_turbu         ! Turbulent diffusion
     integer(ip)            :: kfl_extfo         ! External force
     integer(ip)            :: kfl_saffm         ! Saffman force
     integer(ip)            :: kfl_schem         ! Time integration scheme
     integer(ip)            :: kfl_tstep         ! Time step strategy
     integer(ip)            :: kfl_dmini         ! Minimum diameter strategy >0: eliminate, <0: keep, abs=1: constant, abs=2: mass ratio, abs=3: diameter ratio
     real(rp)               :: denpa             ! Density particle
     real(rp)               :: spher             ! Particle sphericity
     real(rp)               :: diame             ! Particle diameter
     real(rp)               :: calor             ! Calorific capacity
     real(rp)               :: emisi             ! Particle radiation emissivity
     real(rp)               :: scatt             ! Particle radiation scattering factor
     real(rp)               :: diffu             ! Particle diffusion coefficient [m^2/s]
     real(rp)               :: dtime             ! Time step
     real(rp)               :: safet             ! Safety factor
     real(rp)               :: chale             ! Characteristic length
     real(rp)               :: tursc             ! Turbulent schmidt number
     real(rp)               :: prope(mlapr)      ! Default value of properties for each type
     integer(ip)            :: prova(mlapr)      ! ID of variable in properties vector
     !
     ! Thermodynamic
     !
     real(rp)               :: param_dmini       ! Thermodynamic model: minimum diameter
     real(rp)               :: n_drop            ! Number of droplets represented by a parcel
     real(rp)               :: L_vapor           ! Thermodynamic model: enthalpy of vaporization          [ J/kg ] 
     real(rp)               :: cp                ! Thermodynamic model: specific heat capacity            [ J / (kg K) ]        
     real(rp)               :: w                 ! Thermodynamic model: molecular weight                  [ kg / mol ]
     real(rp)               :: T_boiling         ! Thermodynamic model: Boiling temperature               [ K ]
     real(rp)               :: T_crit            ! Thermodynamic model: Critical temperature              [ K ]
     real(rp)               :: P_boiling         ! Thermodynamic model: Boiling pressure                  [ Pa ]
     real(rp)               :: cpcoef_v_chm(6,2) ! Coefficients of NASA polynomial of fuel in the lookup table 
     real(rp)               :: weight_seen       ! Thermodynamic model: weighting factor of seen gas properties for evaluating mean properties, default: 1/3
     type(liquid_state)     :: liq               ! Liquid state
     integer(ip)                        :: kfl_tab_fw              ! index of lookup framework
     type(typ_lookup_framework),pointer :: table_fw                ! lookup framework
     integer(ip)                        :: kfl_spr_fw              ! index of spray terms lookup framework
     type(typ_lookup_framework),pointer :: spr_tab_fw              ! spray terms lookup framework

     integer(ip)                        ::  kfl_W_tab_index,                & ! column index of molecular weight in table 
                                            kfl_k_tab_index,                & ! column index of thermal conductivity in table 
                                            kfl_mu_tab_index,               & ! column index of dynamic viscosity in table 
                                            kfl_cpCoefLT_tab_index,         & ! first column index of low temperature NASA polynomials in table 
                                            kfl_cpCoefHT_tab_index,         & ! first column index of high temperature NASA polynomials in table 
                                            kfl_cpCoefLT_end_tab_index,     & ! last column index of low temperature NASA polynomials in table 
                                            kfl_cpCoefHT_end_tab_index,     & ! last column index of high temperature NASA polynomials in table 
                                            kfl_Yfuel_spr_index,            & ! column index of fuel mass fraction in spray table 
                                            kfl_Dfuel_spr_index,            & ! column index of fuel diffusivity in spray table 
                                            kfl_dkdT_spr_index,             & ! column index of d(k)/d(T) in spray table 
                                            kfl_dmudT_spr_index,            & ! column index of d(mu)/d(T) in spray table 
                                            kfl_dDfueldT_spr_index,         & ! column index of d(Dfuel)/d(T) in spray table 
                                            kfl_T_spr_index                   ! column index of gas temperature in spray table 


  end type typlatyp

  !
  ! Particle
  !
  type latyp
     integer(ip)            :: ilagr             ! Absolute ID
     integer(ip)            :: itype             ! Type of particle
     integer(ip)            :: kfl_exist         ! If I have it
     integer(ip)            :: ielem             ! Last interpolation element
     integer(ip)            :: iboun             ! Boundary element where particle is deposited
     integer(ip)            :: ittim             ! Time iteration
     integer(ip)            :: boundary_set      ! Boundary set where wall intersection
     real(rp)               :: coord(3)          ! Coordinates
     real(rp)               :: veloc(3)          ! Velocity
     real(rp)               :: accel(3)          ! Acceleration
     real(rp)               :: coord_k(3)        ! Coordinates at time k
     real(rp)               :: coord_km1(3)      ! Coordinates at time k-1
     real(rp)               :: v_fluid_k(3)      ! Fluid velocity at time k
     real(rp)               :: v_fluid_km1(3)    ! Fluid velocity at time k-1
     real(rp)               :: v_fluid_km2(3)    ! Fluid velocity at time k-2
     real(rp)               :: acced(3)          ! Drag acceleration
     real(rp)               :: accee(3)          ! External acceleration
     real(rp)               :: acceg(3)          ! Gravity/buoyancy acceleration
     real(rp)               :: stret             ! Stretching factor
     real(rp)               :: t_inject          ! Time particle was injected
     real(rp)               :: t                 ! Time
     real(rp)               :: dt_k              ! Time step: t^k+1-t^k
     real(rp)               :: dt_km1            ! Time step: t^k  -t^k-1
     real(rp)               :: dt_km2            ! Time step: t^k-1-t^k-2
     real(rp)               :: dtg               ! Guessed time step
     real(rp)               :: Cd                ! Drag coefficient
     real(rp)               :: Re                ! Reynolds number
     real(rp)               :: Stk(2)            ! Stokes number (instantaneous & effective) 
     real(rp)               :: dista             ! Distance
     real(rp)               :: coord1d           ! 1D coordinate
     real(rp)               :: sign              ! sign used for 1D coordinate
     !
     ! Thermodynamic model
     !
     real(rp)               :: tempe_k           ! Temperature at time k
     real(rp)               :: tempe_km1         ! Temperature at time k-1
     real(rp)               :: tempe_km2         ! Temperature at time k-2
     real(rp)               :: mass_k            ! Mass at time k
     real(rp)               :: mass_km1          ! Mass at time k-1
     real(rp)               :: mass_km2          ! Mass at time k-2
     real(rp)               :: diam_k            ! Diameter for postprocessing
     real(rp)               :: Temp_fluid_k      ! Seen temperature for postprocessing
     real(rp)               :: Yvap_fluid_k      ! Seen vapour mass fraction for postprocessing
     real(rp)               :: mass_0            ! Initial mass
     real(rp)               :: diam_0            ! Initial diameter
     real(rp)               :: BM                ! Spalding mass transfer number
     !
     ! Redistribution
     !
     integer(ip)            :: mpi_rank          ! MPI rank

   contains

     procedure, pass        :: init              ! Initialize

  end type latyp
  
  !
  ! Injector type
  !
  type typ_injector
     !
     ! Particle properties
     !
     integer(ip)              :: kfl_particle_type

     !
     ! Size distribution
     !
     integer(ip)              :: kfl_size_dist
     integer(ip)              :: kfl_size_lookupfw
     real(rp)                 :: size_diame 
     real(rp)                 :: size_dmin 
     real(rp)                 :: size_dmax 
     real(rp)                 :: size_rr_dbar 
     real(rp)                 :: size_rr_n 
     real(rp)                 :: size_rr_K   
     real(rp),        pointer :: size_list(:)

     !
     ! Temperature
     !
     integer(ip)              :: kfl_tempe
     real(rp)                 :: tempe 
     real(rp),        pointer :: tempe_list(:)

     !
     ! Velocity
     !
     integer(ip)              :: kfl_veloc
     integer(ip)              :: kfl_velfun
     real(rp)                 :: vel_vec(3)
     real(rp)                 :: vel_mag
     real(rp)                 :: vel_ax
     real(rp)                 :: vel_sigma
     real(rp)                 :: vel_angle
     real(rp),        pointer :: vel_list(:,:)
     ! Size dependent velocity:
     real(rp)                 :: vel_mag_s
     real(rp)                 :: vel_diam_L
     real(rp)                 :: vel_diam_s
     real(rp)                 :: vel_ang_max_L
     real(rp)                 :: vel_ang_min_L
     real(rp)                 :: vel_ang_max_s
     real(rp)                 :: vel_ang_min_s

     !
     ! Velocity fluctuations
     !
     integer(ip)              :: kfl_fluct_veloc
     real(rp)                 :: fluct_vel_std

     !
     ! Flow rate of particles
     !
     integer(ip)              :: kfl_flow
     integer(ip)              :: kfl_flowfun
     real(rp)                 :: flow_rate
     integer(ip)              :: num_part

     !
     ! Geometry
     !
     integer(ip)              :: kfl_geometry
     integer(ip)              :: kfl_geo_spatial_dist
     integer(ip)              :: kfl_random 
     real(rp)                 :: geo_coord_min(3)
     real(rp)                 :: geo_coord_max(3)
     real(rp)                 :: geo_coord(3)
     real(rp)                 :: geo_normal(3)
     real(rp)                 :: geo_basis(3,3)
     real(rp)                 :: geo_coord1(3)
     real(rp)                 :: geo_coord2(3)
     real(rp)                 :: geo_coord3(3)
     real(rp)                 :: geo_rad
     real(rp)                 :: geo_radmin
     real(rp)                 :: geo_height
     real(rp),        pointer :: coord_list(:,:)

     !
     ! Injection time management of injector
     !
     real(rp)                 :: time_initial
     real(rp)                 :: time_period
     real(rp)                 :: time_final
     real(rp)                 :: time_cumulative

  end type typ_injector


  type(latyp),    pointer   :: lagrtyp(:)
  type(typlatyp)            :: parttyp(mtyla)    ! Particle types
  
  public :: mlapr
  public :: latyp
  public :: mtyla
  public :: typlatyp
  public :: lagrtyp
  public :: parttyp
  public :: typ_injector
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-07-10
  !> @brief   Initialization
  !> @details Particle initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine init(particle)

    class(latyp) :: particle

    particle % ilagr          = 0
    particle % itype          = 1
    particle % kfl_exist      = 0
    particle % ielem          = 0
    particle % ittim          = 0
    particle % boundary_set   = 0
    particle % iboun          = 0
    particle % coord          = 0.0_rp
    particle % veloc          = 0.0_rp
    particle % accel          = 0.0_rp
    particle % coord_k        = 0.0_rp
    particle % coord_km1      = 0.0_rp
    particle % v_fluid_k      = 0.0_rp
    particle % v_fluid_km1    = 0.0_rp
    particle % v_fluid_km2    = 0.0_rp
    particle % acced          = 0.0_rp
    particle % accee          = 0.0_rp
    particle % acceg          = 0.0_rp 
    particle % stret          = 1.0_rp
    particle % t_inject       = 0.0_rp
    particle % t              = 0.0_rp
    particle % dt_k           = dtime
    particle % dt_km1         = dtime
    particle % dt_km2         = dtime
    particle % dtg            = dtime
    particle % Cd             = 0.0_rp
    particle % Re             = 0.0_rp
    particle % Stk            = 0.0_rp
    particle % dista          = 0.0_rp
    particle % coord1d        = 0.0_rp
    particle % sign           = 1.0_rp
    particle % tempe_k        = 0.0_rp
    particle % tempe_km1      = 0.0_rp
    particle % tempe_km2      = 0.0_rp
    particle % mass_k         = 0.0_rp
    particle % mass_km1       = 0.0_rp
    particle % mass_km2       = 0.0_rp
    particle % diam_k         = 0.0_rp
    particle % Temp_fluid_k   = 0.0_rp
    particle % Yvap_fluid_k   = 0.0_rp
    particle % mass_0         = 0.0_rp
    particle % diam_0         = 0.0_rp
    particle % BM             = 0.0_rp

  end subroutine init

end module def_partis_type
