!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    def_nastin.f90
!> @author  houzeaux
!> @date    2020-04-22
!> @brief   Nastin definitions
!> @details All input variables
!-----------------------------------------------------------------------

module def_nastin_input

  use def_kintyp
  
  !--BEGIN REA GROUP
  !------------------------------------------------------------------------
  !
  ! Physical problem: read in nsi_reaphy
  !
  !------------------------------------------------------------------------

  integer(ip)               :: &
       kfl_advec_nsi,          &      ! Existence of (u.grad)u
       kfl_convection_type_nsi,&      ! Type of convection
       kfl_colev_nsi,          &      ! Coupling with LEVELS
       kfl_cotem_nsi,          &      ! Coupling with TEMPER
       kfl_cotur_nsi,          &      ! Coupling with TURBUL
       kfl_timei_nsi,          &      ! Existence of du/dt
       kfl_fvfua_nsi,          &      ! Frame angular velocity function  
       kfl_fvful_nsi,          &      ! Frame linear velocity function  
       kfl_grtur_nsi,          &      ! Add grad(K) to momemtum equations
       kfl_regim_nsi,          &      ! Flow regime: incompressible/compressible 
       kfl_dynco_nsi,          &      ! Dynamical coupling
       kfl_visco_nsi,          &      ! Viscous term 
       kfl_prthe_nsi,          &      ! Thermodynamic pressure calculation
       kfl_lookg_nsi,          &      ! Thermodynamic pressure calculation
       kfl_surte_nsi,          &      ! Include surface tension
       kfl_force_nsi,          &      ! Force term of the momentum equations
       kfl_vegeo_time_nsi,     &      ! Force term of the momentum equations, time dependent force
       kfl_mfrco_nsi,          &      ! Mass flow rate control activation flag
       kfl_bnods_nsi,          &      ! boundary nodes defined
       kfl_hydro_gravity_nsi,  &      ! Add hydrostatic gravity to NS
       kfl_anipo_nsi,          &      ! Existence of anisotropic viscosity
       kfl_fscon_nsi,          &      ! FS consistent algorithm (solve a Poisson each step)
       mfrse_nsi,              &      ! Set from which the mass flow rate is calculated
       nbtim_nod_nsi,          &      ! total number of time instances and nodes for time-space boundary from file
       nfiel_nsi(2),           &      ! Fields assignement 
       nbnod_nsi,              &      ! Total number of nodes on boundary
       ntabf_nsi                      ! Number of tables for spatial and temporal function tables

  real(rp)                  :: &
       boube_nsi    ,          &      ! Boussinesq volume expansion
       bougr_nsi    ,          &      ! Boussinesq gravity
       boutr_nsi    ,          &      ! Boussinesq reference temperature
       lowtr_nsi    ,          &      ! Low Mach reference temperature
       facca_nsi(3) ,          &      ! Frame angular acceleration vector       
       faccl_nsi(3) ,          &      ! Frame linear acceleration vector  
       fadia_nsi(3) ,          &      ! Frame angular acceleration direction 
       fadil_nsi(3) ,          &      ! Frame linear acceleration direction 
       fvdil_nsi(3) ,          &      ! Frame linear velocity direction     
       fanoa_nsi    ,          &      ! Frame angular acceleration norm   
       fanol_nsi    ,          &      ! Frame linear acceleration norm      
       fcons_nsi    ,          &      ! Convection term
       frotc_nsi(3) ,          &      ! Frame rotation center
       centr_nsi    ,          &      ! Centrifugal force
       fvdia_nsi(3) ,          &      ! Frame angular velocity direction 
       fvela_nsi(3) ,          &      ! Frame angular velocity vector    
       fvell_nsi(3) ,          &      ! Frame linear velocity vector        
       fvins_nsi    ,          &      ! Viscous term
       fvnoa_nsi    ,          &      ! Frame angular velocity norm     
       fvnol_nsi    ,          &      ! Frame linear velocity norm 
       fvpaa_nsi(6) ,          &      ! Frame angular velocity parameters  
       fvpal_nsi(6) ,          &      ! Frame linear velocity parameters  
       gravi_nsi(3) ,          &      ! Gravity vector
       gravb_nsi(3) ,          &      ! Gravity vector for Boussinesq coupling
       grnor_nsi    ,          &      ! Gravity norm
       turbu_nsi(2) ,          &      ! Turbulence parameters
       heihy_nsi    ,          &      ! Height for hydrostatic pressure
       surte_nsi    ,          &      ! Surface tension coeficient (sigma)
       mfrub_nsi    ,          &      ! Target bulk velocity when mass flow control activated
       fsrot_nsi,              &      ! FS rotational factor
       ubpre_nsi    ,          &      ! Bulk velocity from previous time-step
       mfccf_nsi                      ! Coefficient for the mass flow control formula

  integer(ip), pointer      :: &
       lforc_material_nsi(:)  ,&      ! List of material force 
       ntabl_nsi(:)           ,&      ! number of tabulated ct and cp parameters (wind turbines)
       ntabr_nsi(:)                   ! number of tabulated rotational parameters (wind turbines)
  real(rp),    pointer      :: &
       xforc_material_nsi(:,:),&        ! Material force parameters
       velta_nsi(:,:)         ,&        ! tabulated input velocity 
       thrta_nsi(:,:)         ,&        ! tabulated thrust coeff
       powta_nsi(:,:)         ,&        ! tabulated power  coeff
       veave_nsi(:,:)         ,&        ! averaged velocity
       radiu_nsi(:,:)         ,&        ! tabulated dimensional radius
       forcn_nsi(:,:)         ,&        ! tabulated normal force distribution
       forct_nsi(:,:)         ,&        ! tabulated tangential force distribution
       nbtdt_nsi(:)                     ! Time-step of the turbulent inlet database

  real(rp),    pointer         :: &
       bntab_nsi(:,:),            &     ! boundary nodes table for time-space boundary from file
       bnval_nsi(:,:)                   ! boundary values table for time-space boundary from file

  integer(ip),    pointer      :: &
       iboun_nsi(:,:),            &     ! boundary correspondence for time-space boundary from file
       nbnod_pos_nsi(:),          &     ! Start of each table for spatial and temporal functions for nodes
       nbtim_pos_nsi(:),          &     ! Start of each table for spatial and temporal functions for number of times
       nbtim_nod_pos_nsi(:)             ! Start of each table for spatial and temporal functions for temporal values
  !
  ! Fluid properties
  !
  real(rp)                  :: &
       sphea_nsi,              &      ! Specific heat (Cp)
       prthe_nsi,              &      ! Thermodynamics pressure (cst or initial)
       tmass_nsi                      ! Initial mean density related to initial mass(low Mach)

  !------------------------------------------------------------------------
  !
  ! Numerical problem: read in nsi_reanut
  !
  !------------------------------------------------------------------------

  integer(ip)              :: &
       kfl_penal_nsi,         &      ! Penalization
       kfl_prepe_nsi,         &      ! Pressure Penalization - Used to avoid temporal pressure oscilations with local dt
       kfl_dttyp_nsi,         &      ! Local strategy of time step
       kfl_ellen_nsi,         &      ! =0,1 for min/max element length
       kfl_relax_nsi,         &      ! Velocity relaxation strategy
       kfl_relap_nsi,         &      ! Pressure relaxation strategy
       kfl_sgsco_nsi,         &      ! Stabilization convection tracking
       kfl_sgsti_nsi,         &      ! Stabilization time tracking
       kfl_sgsac_nsi,         &      ! Stabilization time tracking accuracy
       kfl_sgsli_nsi,         &      ! Stabilization tracking convection linearization
       kfl_sgscp_nsi,         &      ! Coupling of SGS with grid scale
       kfl_shock_nsi,         &      ! Shock capturing type 
       kfl_tiacc_nsi,         &      ! Temporal accuracy
       kfl_normc_nsi,         &      ! Norm of convergence
       kfl_refer_nsi,         &      ! Difference between solutions
       kfl_linea_nsi,         &      ! Linearization (RHS=0, Picard=1, Newton=2, ExactNewton(adj)=3)
       kfl_tisch_nsi,         &      ! Time integration scheme
       kfl_algor_nsi,         &      ! Type of algorithm ( Monolithic, gauss-seidel, etc.)
       kfl_predi_nsi,         &      ! Predictor corrector
       kfl_taush_nsi,         &      ! Schur complement. Tau strategy 
       kfl_ellsh_nsi,         &      ! Schur complement. =0,1 for min/max element length 
       kfl_updpr_nsi,         &      ! Pressure update
       kfl_intpr_nsi,         &      ! Treatment of the pressure term
       kfl_assem_nsi,         &      ! Assembly type (1,2)
       kfl_asbou_nsi,         &      ! Assembly type (1,2)
       kfl_taust_nsi,         &      ! Tau strategy
       kfl_stabi_nsi,         &      ! Orthogonal SGS
       kfl_limit_nsi,         &      ! Limiter for split OSS
       kfl_trres_nsi,         &      ! Transient residual
       kfl_prtre_nsi,         &      ! Pressure treatment (explicit/implicit)
       kfl_matdi_nsi,         &      ! Dirichlet bc on matrix
       kfl_intfo_nsi,         &      ! Internal force calculation (=0: integral, 1=residual based)
       kfl_press_nsi,         &      ! Integrate pressure term in momentum equations
       kfl_bubbl_nsi,         &      ! Pressure bubble
       kfl_fsgrb_nsi,         &      ! Fractional Step with Gravity or Boussinesq correction
       kfl_stab_corio_nsi,      &      ! Coriolis stabilization
       kfl_stain_nsi,         &      ! Time step to start inner iterations
       kfl_immer_nsi,         &      ! Immersed boundary method
       kfl_grad_div_nsi,      &      ! Use grad an div matrices, also Laplacian despite the name does not reflect it
       misgs_nsi,             &      ! Max # of SGS iterations
       npica_nsi,             &      ! Number of Picard iteration (Newton's lin.)
       neule_nsi,             &      ! # of Euler time steps
       kfl_meshi_nsi,         &      ! Mesh interpolator activation flag
       kfl_savco_nsi,         &      ! Save linear matrix
       kfl_corre_nsi,         &      ! Fractional step correction-like 
       kfl_sosch_nsi,         &      ! Schur complement solver
       kfl_modfi_nsi,         &      ! Modify kfl_fixno_nsi
       kfl_expco_nsi,         &      ! Treat the convective term explicitly, that is, assemble the matrix only in the first ortomin iteration
       kfl_addpr_nsi,         &      ! Add contribution due to pressure in matrix side (do nothing BC') on wall law boundaries        
       kfl_grvir_nsi,         &      ! Add  viscous gradient contribution inside the residual
       kfl_wlare_nsi,         &      ! Initialize time-averaged velocity after restart
       kfl_hydro_nsi,         &      ! Hydrostatic initial state
       kfl_update_hydro_nsi,  &      ! When to update hydrostatic pressure
       kfl_hydro_interface_nsi, &    ! Interface height computation
       mitri_nsi,             &      ! Maximum number of Richardson iterations
       kfl_adres_nsi,         &      ! FS, Tau solver: adaptive residual 
       kfl_incre_nsi,         &      ! Solve pressure in incremental form
       kfl_nota1_nsi,         &      ! Does not stabilize of convective reactive and coriolis term in momentum eq.
       kfl_ini_ts_guess_order_nsi,&  ! Order of the initial guess extrapolation at each time step  - for the moment only velocity. 
       kfl_vector_nsi,        &      ! Vectorized version of Alya
       kfl_press_stab_nsi,    &      ! Pressure stabilization: dt(0) or tau(1)
       kfl_stop_by_wit_nsi,   &      ! Stop by convergence of witness points - uses velocity
       kfl_massm_nsi,         &      ! Consistent mass matrix
       kfl_dampi_nsi,         &      ! Damping - Rayleigh and/or viscous
       k_dim_damp_nsi                ! vertical dimension

  
  real(rp)                 :: &
       top_r_damp_nsi,        &      ! Damping
       top_v_damp_nsi,        &      !
       bot_r_damp_nsi,        &      !
       bot_v_damp_nsi,        &      !
       val_r_damp_nsi,        &      !
       mul_v_damp_nsi,        &      !
       v_geo_damp_nsi(3),     &      !
       penal_nsi,             &      ! Penalization factor
       prepe_nsi,             &      ! Pressure penalization factor
       dtcri_nsi,             &      ! Critical time step
       staco_nsi(4),          &      ! Stability constants
       shock_nsi,             &      ! Shock capturing parameter
       safet_nsi,             &      ! Safety factor for time step
       bemol_nsi,             &      ! Integration of convective term by parts
       sstol_nsi,             &      ! Steady state tolerance
       cotol_nsi,             &      ! Convergence tolerance
       resid_nsi,             &      ! Residual for outer iterations (u)
       resip_nsi,             &      ! Residual for outer iterations (p)
       weigh_nsi,             &      ! Weight of dU/dt in the residual
       relax_nsi,             &      ! Relaxation parameter velocity
       relap_nsi,             &      ! Relaxation parameter pressure
       relsg_nsi,             &      ! Relaxation parameter of subgrid scale
       staco_corio_nsi,       &      ! Coriolis stabilization parameter
       tosgs_nsi,             &      ! Subgrid scale tolerance
       strec_nsi,             &      ! Adaptive dt: Stretching factor
       dampi_nsi,             &      ! Adaptive dt: damping
       epsht_nsi,             &      ! Adaptive dt: eps_R
       epstr_nsi,             &      ! Adaptive dt: eps_A
       xfree_nsi,             &      ! X Coordinate of the plane where to free
       safex_nsi,             &      ! Time function parameter for safety factor
       adres_nsi,             &      ! FS, Tau solver: adaptive residual factor
       toler_nsi,             &      ! FS, Tau solver: adaptive residual factor
       safma_nsi,             &      ! Maximum safety factor
       safeo_nsi,             &      ! Initial safety factor
       saflo_nsi,             &      ! Minimum global safety factor for local time steps
       gamma_nsi                     ! Gamma factor for pressure in momentum equation
 
  !------------------------------------------------------------------------
  !
  ! Boundary conditions: read in nsi_reabcs
  !
  !------------------------------------------------------------------------

  integer(ip)              :: &
       kfl_confi_nsi,         &       ! Confined flow
       kfl_local_nsi,         &       ! Local system of reference
       kfl_conbc_nsi,         &       ! Constant b.c.
       kfl_syntu_nsi,         &       ! Synthetic eddy method (SEM)
       kfl_initi_nsi,         &       ! Initial 
       kfl_inico_nsi,         &       ! Solve coarse grid system
       kfl_inipr_nsi,         &       ! Initial pressure
       kfl_nopen_nsi,         &       ! No penetration condition
       kfl_cadan_nsi,         &       ! Coupling with ADAN
       kfl_aiobo_nsi,         &       ! Alya IO boundary to couple with ADAN 
       itebc_nsi,             &       ! Initial step for boundary condition
       nodpr_nsi,             &       ! Node on which pressure is prescribed
       nodpr_global_nsi,      &       ! Node on which pressure is prescribed (global numbering)
       exfpr_nsi,             &       ! Extends fixpr one layer of elements
       neddy_nsi,             &       ! Number of eddies in the inlet box for SEM
       kfl_imppr_nsi,         &       ! Imposses pressure in nodes w/ fixpr>0 
       kfl_divcorrec_nsi,     &       ! Correction div(u)=0 to initial solution
       nflow_rates,           &       ! Number of flow rates
       kfl_press_lev_nsi,     &       ! Level pressure on a set
       kfl_immersed_nsi,      &       ! Immersed boundary strategy
       kfl_velom_nsi                  ! Use mesh velocity in bc
  real(rp) ::&
       delta_nsi,             &      ! Distance to the wall
       relbc_nsi,             &      ! Boundary condition relaxation
       valpr_nsi,             &      ! Pressure value
       mulpr_nsi,             &      ! Multiplicative factor for weak dirichlet pressure
       velin_nsi(3),          &      ! Initial constant velocity
       poise_nsi(6),          &      ! Parameters for Poiseuille law
       divcorrec_nsi,         &      ! div(u)=0 parameter (alpha)
       fact_immersed_nsi,     &      ! Immersed boundary factor
       val_press_lev_nsi,     &      ! Value of pressure level
       inflow_cosine_nsi             ! Cosine to decide if inflow or outflow
  real(rp) ::&
       press_cadan_nsi               ! Pressure coming from ADAN (ADAN COUPLING)
  real(rp), allocatable    :: &
       Q_cadan_nsi(:),        &      ! Flow to send to ADAN
       P_cadan_nsi(:)                ! Pressure to send to ADAN
  type(bc_nodes), pointer  :: &     
       tncod_nsi(:)                  ! Node code type
  type(bc_nodes), pointer  :: &     
       tgcod_nsi(:)                  ! Geometrical node code type
  type(bc_bound), pointer  :: &     
       tbcod_nsi(:)                  ! Boundary code type
  type(bc_nodes), pointer  :: &     
       tscod_nsi(:)                  ! Spare mesh code type

  integer(ip), pointer  ::&
       kfl_flow_rate_codes_nsi(:),&  ! Flow rates codes 
       kfl_flow_rate_set_nsi(:),  &  ! Flow rates set to compute pressure 
       kfl_flow_rate_normal_nsi(:),& ! Flow rates normals 
       kfl_flow_rate_stfn_nsi(:),&   ! Flow rates space time function number       
       kfl_flow_rate_tfn_nsi(:),&    ! Flow rates time function number 
       kfl_flow_rate_pfn_nsi(:),&    ! Flow rates pump function number 
       kfl_flow_rate_alg_nsi(:)      ! Flow rates algorithm

  real(rp), pointer  ::&
       flow_rate_values_nsi(:), &    ! Flow rate values
       flow_rate_press_nsi(:),  &    ! Flow rate impossed pressure values
       flow_rate_relax_nsi(:),  &    ! Flow rate temporal relaxation
       flow_rate_normals_nsi(:,:)    ! Flow rate values

  !------------------------------------------------------------------------
  !
  ! Optimization problem
  !
  !------------------------------------------------------------------------

  real(rp),   pointer ::                   &
       resdiff_nsi(:,:),                   & ! Partial Derivative of residual w.r.t design var
       dcost_dx_nsi(:),                    & ! Partial Derivative of objective function w.r.t coordinates
       costdiff_nsi(:)                       ! Partial Derivative of objective function w.r.t design var
  
  !------------------------------------------------------------------------
  !
  ! Output and Postprocess: read in nsi_reaous
  !
  !------------------------------------------------------------------------

  integer(ip)              :: &
       kfl_exacs_nsi,         &      ! Exact solution for the NS eq.
       kfl_exfix_nsi,         &      ! Fixity imposed for exact solution
       kfl_psmat_nsi,         &      ! Postscript of matrix profile
       kfl_inert_nsi,         &      ! Velocity in inertial frame of ref.
       kfl_miniapp_nsi,       &      ! Nastin miniapp
       avste_nsi                     ! Steps for the element averaging
  real(rp) ::&
       expar_nsi(10),           &    ! Exact solution parameters
       cloth_nsi,               &    ! CLO
       metab_nsi,               &    ! MET
       wetme_nsi,               &    ! WME
       ambie_nsi,               &    ! TA
       radia_nsi,               &    ! TR
       relat_nsi,               &    ! RH
       avtim_nsi,               &    ! Start averaging time
       entim_nsi                     ! End ensemble time
  !--END REA GROUP

end module def_nastin_input
!> @}
