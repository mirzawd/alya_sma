!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module def_turbul
  !------------------------------------------------------------------------
  !    
  ! Heading for the Turbulence routines
  !
  !------------------------------------------------------------------------
  use def_kintyp
  use mod_ADR,   only : ADR_typ

  !------------------------------------------------------------------------
  ! Parameters
  !------------------------------------------------------------------------

  real(rp), parameter      :: &
       zetur=epsilon(1.0_rp), &
       vonka_tur=0.41_rp
  integer(ip),   parameter :: &
       ncoef_tur=10,          &      ! # coefficient for properties
       nipar_tur=10,          &      ! # integer parameters
       nrpar_tur=20,          &       ! # real parameters
       max_mater_tur = 1000          ! maximum number of materials (allowed by the moment)
  integer(ip),   parameter ::         &
       TUR_BEFORE_TIME_STEP=1,        &
       TUR_BEFORE_GLOBAL_ITERATION=2, &
       TUR_BEFORE_INNER_ITERATION=3

  !------------------------------------------------------------------------
  ! Physical problem: read in tur_reaphy
  !------------------------------------------------------------------------
!--BEGIN REA GROUP
  integer(ip)              :: &
       kfl_model_tur,         &      ! Turbulence model
       kfl_timei_tur,         &      ! Existence of dT/dt 
       kfl_cotem_tur,         &      ! Temperature coupling
       kfl_advec_tur,         &      ! Existence of advection (=1 uses veloc)
       kfl_colev_tur,         &      ! Coupling with LEVELS
       kfl_ddesm_tur,         &      ! Delayed detached eddy simulation model (tur_sstkow)
       kfl_sasim_tur,         &      ! Scale adaptive simulation model (tur_sstkow)
       kfl_rotat_tur,         &      ! Rotating reference frame
       kfl_inifi_tur(3),      &      ! Initial fields
       inits_tur,             &      ! Initial time step
       nturb_tur,             &      ! # turbulence variables
       nfiel_tur(3),          &      ! Fields assignement
       ipara_tur(nipar_tur),  &      ! Turbulence integer parameters
       lawde_tur,             &      ! Law for rho
       lawvi_tur,             &      ! Law for mu   
       kfl_kxmod_tur,         &      ! modification of k eps model (RNG or REALIZABLE)
       kfl_logva_tur,         &      ! works with logarthmic of unknown functions 
       kfl_discd_tur,         &          ! dissipation wake model 
       ldiss_material_tur(max_mater_tur) ! list of materials for extra dissipation wake
     
  real(rp)                 :: &
       boube_tur,             &      ! Boussinesq coupling beta
       grnor_tur,             &      ! Boussinesq coupling |g|
       gravi_tur(3),          &      ! Boussinesq coupling g
       prtur_tur,             &      ! Turbulent Prandtl number  
       param_tur(nrpar_tur),  &      ! Turbulence parameters
       densi_tur(ncoef_tur),  &      ! Density (rho)
       visco_tur(ncoef_tur),  &      ! Viscosity (mu)
       densa_tur,             &      ! Air density
       visca_tur,             &      ! Air viscosity
       cddes_tur,             &      ! CDDES - smagorinsky-like constant for DES
       inv_l_max                     ! Inverse of Max. mixing length (CFDWind2 model)

  !------------------------------------------------------------------------
  ! Numerical problem: read in tur_reanut
  !------------------------------------------------------------------------

  integer(ip)              :: &     
       kfl_algor_tur,         &      ! Algorithm: uncoupled (1) / coupled (2)
       kfl_clipp_tur,         &      ! Clipping strategy
       kfl_ellen_tur,         &      ! =0,1 for min/max element length
       kfl_normc_tur,         &      ! Norm of convergence
       kfl_relax_tur,         &      ! Relaxation strategy
       kfl_repro_tur,         &      ! Stabilization based on residual projection
       kfl_shock_tur,         &      ! Shock capturing type 
       kfl_taust_tur,         &      ! Tau calculation option
       kfl_tiacc_tur,         &      ! Temporal accuracy
       kfl_tisch_tur,         &      ! Time integration scheme
       kfl_walgo_tur,         &      ! Wall distance algorithm
       kfl_weigh_tur,         &      ! Weighting of dT/dt
       kfl_assem_tur,         &      ! Assembly
       kfl_ortho_tur,         &      ! Stabilization strategy
       kfl_limit_tur,         &      ! Limiter
       kfl_produ_tur,         &      ! Production term 2*mut*Sij*Sij or 2*mut*Sij*dui/dxj
       kfl_meshi_tur,         &      ! Mesh interpolator activation flag
       miinn_tur,             &      ! Max # of iterations
       niter_tur,             &      ! Max # of inner iterations
       neule_tur,             &      ! # of Euler time steps
       kfl_sgsti_tur,         &      ! Subscale time tracking
       kfl_sgsno_tur,         &      ! Subscale non-linear tracking
       kfl_tibub_tur,         &      ! Time integration of bubble
       kfl_sgsac_tur,         &      ! SGS time accuracy
       kfl_lmaxi_tur                 ! max mixing length using mellor yamada
  

  real(rp)                 :: &
       staco_tur(3),          &      ! Stability constants
       shock_tur,             &      ! Shock capturing parameter
       sstol_tur,             &      ! Steady state tolerance
       safet_tur,             &      ! Safety factor for time step
       cotol_tur,             &      ! Convergence tolerance
       relax_tur,             &      ! Relaxation factor
       safex_tur,             &      ! Time function parameter for safety factor
       safma_tur,             &      ! Maximum safety factor
       safeo_tur,             &      ! Initial safety factor
       clipfac_tur,           &      ! Clipping factor
       saflo_tur,             &      ! Minimum global safety factor when using local time step
       bemol_tur                     ! Convective term integration

  !------------------------------------------------------------------------
  ! Output and Postprocess: read in tur_reaous
  !------------------------------------------------------------------------

  integer(ip)              :: &
       kfl_exacs_tur

  !------------------------------------------------------------------------
  ! Boundary conditions: read in tur_reabcs
  !------------------------------------------------------------------------

  integer(ip)              :: &
       kfl_inidi_tur,         &      ! Initial diffusion problem
       kfl_wallw_tur,         &      ! w-type of boundary condition (0=classical,1=bredberg)
       kfl_infl1_tur,         &      ! Inflow type for 1st turbulence variable
       kfl_infl2_tur,         &      ! Inflow type for 2nd turbulence variable
       kfl_usrbc_tur,         &      ! User boundary condition
       kfl_initi_tur,         &      ! Initial solution
       kfl_valbc_tur(4),      &      
       npnat_tur=1_ip                           ! # variable natural bc
  real(rp)                 :: &
       delta_tur,             &      ! Distance to the wall for wall functions       
       turin_tur,             &      ! Turbulent intensity    
       hdiam_tur,             &      ! Hydraulic diameter
       turle_tur,             &      ! Turbulent length scale  
       nutnu_tur,             &      ! Ratio nut/nu
       rebcs_tur,             &      ! Relaxation of boundary condition
       xinit_tur(4)                  ! Initial value
  type(bc_nodes), pointer  :: &     
       tncod_tur(:)                  ! Node code type
  type(bc_nodes), pointer  :: &     
       tgcod_tur(:)                  ! Geometrical node code type
  type(bc_bound), pointer  :: &     
       tbcod_tur(:)                  ! Boundary code type
!--END REA GROUP
  !------------------------------------------------------------------------
  ! Others
  !------------------------------------------------------------------------
  !
  ! Boundary conditions
  !
  integer(ip), pointer     :: &
       kfl_fixno_tur(:,:,:),  &      ! Nodal fixity 
       kfl_funno_tur(:,:),    &      ! Functions number for node bc
       kfl_funtn_tur(:,:),    &      ! Functions type for node bc
       kfl_fixbo_tur(:)              ! Element boundary fixity
  real(rp),    pointer     :: &
       bvess_tur(:,:,:),      &      ! Essential bc (or initial) values
       bvnat_tur(:,:,:)              ! Natural bc values
  
  !
  ! Iteration counters
  !
  integer(ip)              :: &     
       itera_tur,             &      ! Inner iteration counter
       ittot_tur,             &      ! Total number of iteration 
       kfl_goite_tur,         &      ! Keep iterating
       kfl_stead_tur,         &      ! Steady-state has been reached 
       kfl_walld_tur,         &      ! Wall distance
       kfl_lwnei_tur,         &      ! First neighbor
       kfl_ustar_tur,         &      ! Friction velocity
       kfl_grve2_tur,         &      ! Velocity 2nd order gradients
       kfl_grsqk_tur,         &      ! Projection of sqrt(k)
       kfl_grono_tur,         &      ! Smoothing
       kfl_greps_tur,         &      ! Smoothing 
       kfl_grk12_tur,         &      ! Smoothing
       kfl_grphi_tur,         &      ! Smoothing grad(phi)
       kfl_vorti_tur,         &      ! Smoothing grad(u)
       kfl_avvel_tur,         &      ! Average velocity needed
       kfl_adapt_tur,         &      ! Adaptive b.c.
       kfl_fixn6_tur,         &      ! Inlet b.c. exists
       kfl_fixn8_tur,         &      ! Freestream b.c. exists
       kfl_tiaor_tur,         &      ! Original time accuracy
       nzmat_tur,             &      ! Matrix size
       nzrhs_tur,             &      ! RHS size
       kfl_exist_fixi7_tur           ! exists fixity of type 7 

      
  real(rp)                 :: &
       dtinv_tur,             &      ! 1/dt
       dtcri_tur,             &      ! Critical time step
       cpuit_tur,             &      ! CPU time inner iterarions
       resid_tur,             &      ! Residual total for outer iterations
       dtmax_tur                     ! Stores maximum time step when using local time step
  !
  ! Dimensions
  !
  integer(ip)              :: &
       ndofn_tur,             &      ! # of d.o.f. of the turbulence unknown
       ndof2_tur,             &      ! ndofn_tur*ndofn_tur
       nunkn_tur,             &      ! # of unknonws ndofn*npoin  
       ncomp_tur,             &      ! Number of components of the turbulence unknown
       nprev_tur,             &      ! Previous time step or iteration
       iunkn_tur                     ! Current turbulent variable
  !
  ! Models
  !
  logical(lg)                 :: &
       TUR_FAMILY_K_EPS,         &
       TUR_FAMILY_K_OMEGA,       &       
       TUR_ONE_EQUATION_MODEL,   &       
       TUR_TWO_EQUATION_MODEL,   &       
       TUR_SPALART_ALLMARAS,     &   
       TUR_K_XU_CHIEN,           &                   
       TUR_K_EPS_STD,            &
       TUR_K_EPS_LAUNDER_SHARMA, &   
       TUR_K_EPS_CHIEN,          &        
       TUR_K_EPS_LAM_BREMHORST,  & 
       TUR_K_EPS_JAW_HWANG,      &
       TUR_K_EPS_NAGANO_TAGAWA,  &
       TUR_K_EPS_V2_F,           &
       TUR_K_EPS_PHI_F,          &
       TUR_TWO_LAYER_RODI,       &        
       TUR_TWO_LAYER_XU_CHEN,    &      
       TUR_K_OMEGA,              &                
       TUR_K_OMEGA_BREDBERG,     &
       TUR_SST_K_OMEGA,          &
       TUR_TKE_SGS
  !  
  ! Others
  !
  integer(ip), allocatable :: & 
       lwnei_tur(:)                  ! List of nearest wall neighbors
  real(rp)                 :: &
       relpa_tur(2),          &      !
       avvel_tur,             &      ! Averaged velocity for b.c.'s
       avkin_tur                     ! Averaged k for b.c.'s
  real(rp),    pointer     :: &
!       walld_tur(:),          &      ! Wall distance
       ustar_tur(:),          &      ! Friction velocity
       grve2_tur(:),          &      ! Velocity 2nd order gradients
       grsqk_tur(:),          &      ! grad(sqrt(k))
       tur_max_mixlen(:)             ! maximum mixing length, (used when cfdw2 model)
 
  !
  ! Internal variables
  !
  real(rp),    pointer     :: &
       dunkn_tur(:),          &      ! DUNKN_TUR: Delta unknown for Aitken relaxation strategy   
       unold_tur(:,:),        &      ! Save previous unknown for postprocess
       turvi_tur(:,:,:)              ! turbulent viscosity at gauss points
  real(rp),    pointer     :: &
       grk12_tur(:),          &      ! Smoothing of near wall values
       grono_tur(:),          &      ! Smoothing of near wall values
       greps_tur(:),          &      ! Smoothing of near wall values
       grphi_tur(:,:),        &      ! grad(phi)
       unpro_tur(:,:),        &      ! Turbulence projection
       detur_tur(:),          &      ! Projected variable density
       vitur_tur(:),          &      ! Projected variable viscosity
       vorti_tur(:),          &      ! Project W (vorticity)
       fddes_tur(:),          &      ! to gather FDDES postprocess information
       gddes_tur(:),          &      ! to gather GDDES (=fd) postprocess information
       sstf1_tur(:),          &      ! to gather SST blending function postprocess information
       sstf2_tur(:),          &      ! to gather SST blending function postprocess information
       sasso_tur(:),          &      ! to gather postprocess information, sas source term
       avkey_tur(:),          &      ! Time-averaged turbulent kinetic energy
       avome_tur(:),          &      ! Time-averaged omega
       avtvi_tur(:),          &      ! Time-averaged turbulent viscosity
       olded_tur(:),          &      ! eddy viscosity of privious time step
       produ_tur(:),          &      ! Smoothed production term
       tupr2_tur(:),          &      ! turbulence projection 
       tuprr_tur(:),          &      ! reactive term projection (split oss)
       tupgr_tur(:,:),        &      ! unknown gradient projection (shock capturing)
       unprr_tur(:,:),        &      ! Turbulence projectio
       unpgr_tur(:,:,:),      &      ! Turbulence gradient projections
       bubbl_tur(:,:)                ! Bubble tracking in time

  real(rp) ::&
       avtim_tur,             &      ! Accumulated time for time-averaging
       pabdf_tur(10)                 ! BDF parameters, also used for CN
  
  integer(ip) ::&
       nbdfp_tur                     ! Number of terms in the temporal derivative

  ! 
  ! Manufactured solution
  !
  real(rp)                 :: &
       err01_tur(2),          & ! L1 error T
       err02_tur(2),          & ! L2 error T
       err0i_tur(2),          & ! Linf error T
       err11_tur(2),          & ! L1 error grad(T)
       err12_tur(2),          & ! L2 error grad(T)
       err1i_tur(2)             ! Linf error grad(T)      
  !------------------------------------------------------------------------
  ! Optimization
  !------------------------------------------------------------------------

  real(rp),   pointer ::                   &
       veloc_fix(:,:)                       ! Partial Derivative of turbulent residuals w.r.t design var
       
  real(rp),   pointer ::                   &
       resdiff_tur(:,:)                       ! Partial Derivative of turbulent residuals w.r.t design var
       
  type(r2p),   pointer                  :: &
       Rhsadjtur_tur(:)                      ! Right hand side of the ith turbulence adjoint equation coming from jth one
  
  real(rp), pointer      ::                &
       sens_wall(:,:)                        ! mesh sensitivity values respecto to the wall node coordinates
  !
  ! ADR type
  !
  type(ADR_typ)            :: &
       ADR_tur(2)               ! ADR type
 !------------------------------------------------------------------------
  ! maximum mixing length
  !------------------------------------------------------------------------
  real(rp),    pointer    :: numer(:), denom(:), z_wall(:)
  

end module def_turbul
