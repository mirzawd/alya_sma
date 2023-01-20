!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module def_chemic
  !------------------------------------------------------------------------
  !****f* Partis/def_chemic
  ! NAME
  !    def_chemic
  ! DESCRIPTION
  !    Heading for the Partis routines
  ! USES
  ! USED BY
  !    Almost all
  !***
  !------------------------------------------------------------------------
  use def_kintyp,       only : ip, rp, r3p, bc_nodes, bc_bound, r2p
  use def_kermod,       only : max_lookup_fw
  use def_kermod,       only : max_ann_fw
  use mod_interp_tab,   only : typ_tab_coord, typ_lookup_table, typ_lookup_framework
  use mod_ADR,          only : ADR_typ
#ifdef CANTERA
  use cantera,          only : phase_t
#endif

  implicit none
  !------------------------------------------------------------------------
  ! Parameters
  !------------------------------------------------------------------------

  integer(ip), parameter ::                  &
       lun_times_chm   = 1911,               &  ! Time step info
       lun_time2_chm   = 1912,               &  ! Time step targets
       lun_resu1_chm   = 1921,               &  ! Result 1
       lun_resu2_chm   = 1922,               &  ! Result 2
       lun_resu3_chm   = 1923,               &  ! Result 3
       lun_resou_chm   = 1931,               &  ! Meteo source file
       lun_spcvg_chm   = 1932,               &  ! Species convergence file
       lun_droplet_chm = 1941                   ! Droplet results unit

  character(150)                        :: &
       fil_times_chm,                      &  ! Time step
       fil_time2_chm,                      &  ! Time step target
       fil_resou_chm,                      &  ! Meteo
       fil_droplet_chm                        ! Droplet results file

  real(rp),      parameter :: &
       zepts = epsilon(1.0_rp)
  integer(ip),   parameter              :: &
       npara_chm=10,                       &   ! # parameters for sets
       npart_chm=10                            ! # temperature parameters
!--BEGIN REA GROUP

  !------------------------------------------------------------------------
  ! TYPES
  !------------------------------------------------------------------------
  type typ_chm_mixedEq_group
      character(5)              :: name               ! Name of group
      integer(ip)               :: kfl_grtype         ! Type of group
      integer(ip)               :: nequa              ! Number of equations in group
      integer(ip)               :: i_start            ! Start index of group
      integer(ip)               :: i_end              ! End index of group
      integer(ip)               :: kfl_therm_phor     ! Group is diffused by thermophoretic velocity
  end type typ_chm_mixedEq_group

  type typ_chm_mixedEq_equa
      character(5)              :: name               ! Name of equation
      integer(ip)               :: kfl_eqtype         ! Type of equation
      integer(ip)               :: kfl_source_type    ! Type of equation volumetric sources
      integer(ip)               :: kfl_source_fw      ! Index of tabulated source term framework
      integer(ip)               :: kfl_diffsource_fw  ! Index of tabulated source term framework (Diffusion)
      integer(ip)               :: kfl_premsource_fw  ! Index of tabulated source term framework (Premixed)
      integer(ip)               :: kfl_source_ann     ! Index of ANN source term framework
      integer(ip)               :: kfl_source_col     ! Column in source term framework
      integer(ip)               :: kfl_diffsource_col ! Column in source term framework (Diffusion)
      integer(ip)               :: kfl_premsource_col ! Column in source term framework (Premixed)
      integer(ip)               :: kfl_consum_col     ! Column in source term framework for consumption rates
      integer(ip)               :: kfl_source_split   ! Group is split by source terms
      integer(ip)               :: kfl_ini_type       ! Type of initialization
      integer(ip)               :: kfl_ini_field      ! Index of initial field
      integer(ip)               :: kfl_fix_diffusion  ! Activation constant diffusion
      real(rp)                  :: ini_value          ! Initial value
      integer(ip)               :: kfl_ieq_mean       ! Index of mean feald for a variance
      integer(ip)               :: kfl_do_post        ! Whether or not postprocess this variable as named unknown
      real(rp)                  :: Lewis              ! Lewis number of transported scalar
      real(rp)                  :: diffusivity        ! Diffusivity
  end type typ_chm_mixedEq_equa

  !------------------------------------------------------------------------
  ! Physical problem: read in chm_reaphy
  !------------------------------------------------------------------------

  integer(ip)                           :: &
       kfl_model_chm,                      & ! PDE-ODE model
       kfl_timei_chm,                      & ! Existence of du/dt
       kfl_advec_chm,                      & ! Existence of (a.grad)u
       kfl_diffu_chm,                      & ! Existence of -div[k. grad(u)]
       kfl_transport_chm,                  & ! Strategy for transport properties
       kfl_norma_chm,                      & ! Normalize mass fractions
       kfl_premix_chm,                     & ! Flag to identify the combustion model: premixed or non-premixed
       kfl_control_gr_chm,                 & ! Index of control group
       kfl_varYc_chm,                      & ! Flag to define Yc variance model
       kfl_varZ_chm,                       & ! Flag to define Z variance model
       kfl_wallc_chm,                      & ! Flag to impose a zero source term at walls for the Flamelet combustion model
       kfl_radia_chm,                      & ! Flag to activate radiation model for Flamelet model
       kfl_pfa_chm,                        & ! Flag to activate PFA reduction model
       kfl_key_chm,                        & ! Flag to determine number of key species
       kfl_cont_chm,                       & ! Number of controlling varibale TDAC
       kfl_z_chm,                          & ! Flag to determine z computation type
       kfl_freq_chm,                       & ! Number of timesteps between DAC reductions
       kfl_tdac_write_chm,                 & ! Write Reduced table
       kfl_spray_chm,                      & ! Flag to activate spray model
       kfl_ufpv_chm,                       & ! Unsteady Flamelet Progress Variable model
       kfl_lookg_chm,                      & ! Lookup on Gauss points
       kfl_min_srcfw_chm,                  & ! lowest framework index involved in source terms
       kfl_max_srcfw_chm,                  & ! highest framework index involved in source terms
       kfl_min_src_annfw_chm,              & ! lowest ANN framework index involved in source terms
       kfl_max_src_annfw_chm,              & ! highest ANN framework index involved in source terms
       kfl_W_tab_index_chm,                & ! column index of molecular weight in table
       kfl_k_tab_index_chm,                & ! column index of thermal conductivity in table
       kfl_mu_tab_index_chm,               & ! column index of dynamic viscosity in table
       kfl_cpCoefLT_tab_index_chm,         & ! first column index of low temperature NASA polynomials in table
       kfl_cpCoefHT_tab_index_chm,         & ! first column index of high temperature NASA polynomials in table
       kfl_cpCoefLT_end_tab_index_chm,     & ! last column index of low temperature NASA polynomials in table
       kfl_cpCoefHT_end_tab_index_chm,     & ! last column index of high temperature NASA polynomials in table
       kfl_T_tab_index_chm,                & ! column index of temperature in table
       kfl_DtRho_tab_index_chm,            & ! column index of Dt*rho in table
       kfl_max_nvar_lookup_chm,            & ! maximum column index used in property, and sourceterm lookup
       kfl_max_nvar_ann_in_chm,            & ! maximum input size with ANNs
       kfl_max_nvar_ann_out_chm,           & ! maximum output size with ANNs
       kfl_int_chm,                        & ! Integration Strategy Finite rate (CVODE/ODEPIM/PYJAC)
       kfl_spec_name_chm,                  & ! Number of Fields (finite rate)
       kfl_entropy_chm,                    & ! Activate entropy stable viscosity stabilization method
       kfl_multimod_chm,                   & ! Activate Multiregime model
       kfl_disable_entpre_chm,             & ! Otion to bypass RK prediction
       kfl_field_chm(4),                   & ! Flag to activate the initialization by fields (1) Species fields,
                                             !  (2) Temperature field, (3) alpha (CMC model), (4) Soot fields
       nclas_chm,                          & ! Number of classes
       ngrou_chm,                          & ! Number of equation groups
       nspec_chm                             ! Total number of species

  !
  ! Flamelet model variables
  !
  type(typ_tab_coord), pointer  :: &
      table_coords(:)              ! control variable discretization for lookup table

  type(typ_lookup_table), pointer :: &
      table_tab                    ! Table for lookup table


  integer(ip), target                   :: &
       kfl_tab_fw_chm,                     & ! Index of table lookup framework
       kfl_post_fw_chm,                    & ! Index of postprocessing table lookup framework
       kfl_hrr_fw_chm,                     & ! Index of heat release rate lookup framework
       kfl_hrr_col_chm,                    & ! Index of HRR column within the framework
       kfl_zg_fw_chm,                      & ! Index of Zgrad lookup framework
       kfl_zg_col_chm,                     & ! Index of Zgrad column within the framework
       kfl_tab_fw_chm_diff,                & ! Diffusion table framework for multi-regime model
       kfl_tab_fw_chm_prem                   ! Premixed table framework for multi-regime model


  type(typ_lookup_framework), pointer :: &
      table_fw,                  & ! Framework for lookup table
      posttable_fw                 ! Framework for lookup table

  integer(ip)                           :: &
      kfl_fw_has_sources(max_lookup_fw),   & ! Which frameworks to be called for sources
      kfl_fw_has_hybrid_sources(max_lookup_fw),   & ! Which frameworks to be called for multi-regime sources
      kfl_annfw_has_sources(max_ann_fw)      ! Which ANN frameworks to be called for sources

  integer(ip), pointer                  :: &
      kfl_fw_src_equa_list(:,:),           & ! Which equations are related to the active frameworks
      kfl_annfw_src_equa_list(:,:),        & ! Which equations are related to the active ANN frameworks
      kfl_table_ind_chm(:,:)                 ! Save starting indexes per element


  ! Flamelet properties on Gauss points
  type(r2p),pointer        :: &
       zgrad_gp(:),           &  ! Zgradient - lookedup from diffusion table            
       flame_index_gp(:)         ! Flame index 


  type(r3p),pointer        :: &
       mass_gp(:),   &             ! mass source term
       DtRho_gp(:),  &             ! mass source term
       hk_gp(:),     &             ! enthalpy at gauss points (sensible + chemical)
       Yk_ssm_gp(:), &             ! Tabulated species for soot model
       massConsumption_gp(:)       ! Tabulated Consumption rates

  !
  ! mixedEquations variables
  !
  type(typ_chm_mixedEq_group), pointer :: &
      mixedEq_groups_chm(:)       ! Groups of equations

  type(typ_chm_mixedEq_equa), pointer  :: &
      mixedEq_eqs_chm(:)          ! Equations

  !
  ! Indeces to set for some compatibility concerns
  !
  integer(ip)              :: &
      kfl_izmean_chm,         &
      kfl_icmean_chm,         &
      kfl_icvar_chm,          &
      kfl_izvar_chm

  !
  ! Finite-rate combustion model
  !
  real(rp),  pointer                    :: &
       grad_Yk(:,:,:),                     & ! Gradient of species Yk
       grad_Yk_hk(:,:,:),                  & ! Gradiente of species Yk times the enthalpy (chem.+sens.) of species k
       grad_T(:,:),                        & ! Gradient of temperature
       grad_h(:,:),                        & ! Gradient of enthalpy (chem.+sens.)
       aux_nodes(:),                       & ! Auxiliary variable at nodes
       coeff_cp_k(:,:),                    & ! NASA polinomial coefficients for individual species
       W_k(:),                             & ! Species molecular weights
       Y_k_n(:,:,:),                       & ! Local mass fractions of all species for PFA reduction methods
       enthalpy_transport_nodes(:,:),      & ! Enthalpy transport by diffusion
       grad_k_single(:,:),                 & ! Gradient related to a given field (e.g. species, etc.). Used as an auxiliary variable
       Le_k(:),                            & ! Lewis number for each species
       hrr_chm(:),                         & ! Heat release field (W/m3)
       hrr_avg_chm(:)                        ! Averaged heat release field (W/m3)

  !
  ! CMC turbulent combustion model
  !

  integer(ip), pointer                  :: &
      kfl_bc_type_spec_CMC_chm(:),         & ! It values 1 if for all the physical points the boundary conditions for all the
                                             !  chemical species are of the same type and 0 otherwise
      write_uncond_spec_CMC_chm(:),        & ! List of the number of species for which <Y> has to be written
      write_cond_spec_CMC_chm(:),          & ! List of the number of species for which <Y|Z> has to be written
      write_iZ_spec_CMC_chm(:),            & ! List of mixture fractions to be written
      chem_int_iZ_CMC_chm(:),              & ! List of mixture fraction were to carry chemical integration
      transf_spec_CMC_chm(:)                 ! Species to be transferred from CMC to CFD


  real(rp),  pointer                    :: & ! The order of the dimensions is (mixt fraction, points, variables)
      enthalp_CMC_chm(:,:),                & ! Conditional enthalpy h at nodes (nZ,npoin)
      temp_CMC_chm(:,:),                   & ! Conditional temperature at nodes (nZ,npoin)
      Yk_CMC_chm(:,:,:),                   & ! Conditional mass fractions for species (and soot) Y_k at nodes (nZ,npoin,nclas)
      src_Yk_CMC_chm(:,:,:),               & ! Conditional chemical source term for species (and soot) (nZ,npoin,nclas)
      Yk_int_CMC_chm(:,:),                 & ! Unconditional mass fractions for species (and soot) Y_k at nodes (npoin,nclas)
      densi_int_CMC_chm(:),                & ! Integrated density from CMC (npoin)
      visco_lam_int_CMC_chm(:),            & ! Integrated laminar dynamic viscosity from CMC (npoin)
      condu_int_CMC_chm(:),                & ! Integrated thermal conductivity from CMC (npoin)
      sphea_int_CMC_chm(:),                & ! Integrated specific heat from CMC (npoin)
      enthalp_int_CMC_chm(:),              & ! Unconditional enthalpy for CMC Y_k at nodes (npoin)
      temp_int_CMC_chm(:),                 & ! Unconditional temperature for CMC Y_k at nodes (npoin)
      src_Yk_int_CMC_chm(:,:),             & ! Unconditional chemical source term dYi/dt for species (and soot) (npoin,nclas)
      hrr_mass_cond_CMC_chm(:,:),          & ! Heat release (nZ,npoin)
      veloc_CMC_chm(:,:),                  & ! Velocity field from CFD calculation (ndime,npoin). It may be unconditional or
                                             !  conditional depending on whether interpolations are done
      Zavg_CMC_chm(:),                     & ! Mixture fraction field from CFD (npoin)
      Zvar_CMC_chm(:),                     & ! Mixture fraction variance field from CFD (npoin)
      Xtot_CMC_chm(:),                     & ! Total scalar dissipation rate from CFD (npoin) (sgs+solved)
      Xtot_whole_mesh_CMC_chm(:,:),        & ! Total scalar dissipation rate from CFD (nZ,npoin) (sgs+solved) (used when doing
                                             !  interpolations)
      grad_Zavg_CMC_chm(:,:),              & ! Gradient of mixture fraction field from CFD (ndime,npoin)
      grad_Zavg_CFD_chm(:,:),              & ! Gradient of mixture fraction obtained in CFD execution (ndime,npoin)
      turb_kin_visc_CMC_chm(:),            & ! Mass turbulent diffusion coefficient from CFD (npoin)
      deriv2_Yk_CMC_chm(:,:,:),            & ! Second derivative in mixture fraction direction for the conditional species
                                             !  (nZ,npoin,nclas)
      deriv2_enthalp_CMC_chm(:,:),         & ! 2nd derivative in mixture fraction direction for the conditional enthalpy (nZ,npoin)
      bvess_CMC_chm(:,:,:),                & ! Essential bc conditional values for CMC (nZ,npoin,nvar_CMC_chm)
      bvess_ufield_CMC_chm(:,:),           & ! Essential bc unconditional fields (npoin,nvar_CMC_chm)
      T_bc_CMC_chm(:),                     & ! Temperature at boundaries if kfl_weigh_in_eq_CMC_chm activated
      react_scal_bc_CMC_chm(:,:),          & ! Reactive scalars (enthalpy and species) at the Z=0 and Zs if
                                             !  kfl_weigh_in_eq_CMC_chm activated
      Z_CMC_chm(:),                        & ! Mixture fraction vector for CMC (nZ)
      diff_Z_CMC_chm(:),                   & ! Vector with the increments in mixture fractions for Z_CMC_chm (nZ-1)
      Xnormalized_prof_CMC_chm(:),         & ! Scalar dissipation rate normalized profile (nZ)
      rscal_inert_CMC_chm(:,:),            & ! Inert profile for reactive scalars (enthalpy and species) (nZ,nclas+1)
      rscal_equil_CMC_chm(:,:),            & ! Equilibrium profile for reactive scalars (enthalpy and species) (nZ,nclas+1)
      temp_inert_CMC_chm(:),               & ! Inert profile for temperature (nZ)
      temp_equil_CMC_chm(:),               & ! Equilibrium profile for temperature (nZ)
      Z_AMC_CMC_chm(:),                    & ! Mixture fraction vector for AMC model (nZ_AMC_CMC_chm)
      S_AMC_CMC_chm(:),                    & ! Segregation factor vector for AMC model (nS_AMC_CMC_chm)
      Xintegrated_table_AMC_CMC_chm(:,:),  & ! Table that contains the integrated profiles for scalar dissipation rate for AMC
                                             !  model from CMC model (nZ_AMC_CMC_chm, nS_AMC_CMC_chm)
      PDF_min_CMC_chm(:),                  & ! PDF minima values for each mixture fraction
      ! Variables for codes communication
      send_CMC_to_CFD_chm(:,:,:),          & ! CMC sends to CFD
      send_CMC_to_CFD_post_chm(:,:),       & ! CMC sends to CFD (post-processing variables only)
      send_CFD_to_CMC_chm(:,:),            & ! CFD sends to CMC
      rec_CMC_from_CFD_chm(:,:),           & ! CMC receives from CFD
      rec_CFD_from_CMC_chm(:,:,:),         & ! CFD receives from CMC
      rec_CFD_old_from_CMC_chm(:,:,:),     & ! CFD receives from CMC from previous time step
      rec_CFD_from_CMC_post_chm(:,:),      & ! CFD receives from CMC (post-processing variables only)
      turb_kin_visco_CFD_chm(:),           & ! Turbulent kinematic viscosity from CFD to be transfered to CMC
      weights_int_RK_CMC_chm(:,:),         & ! Internal weights for RK for mixt. fract. diff.
      weights_ext_RK_CMC_chm(:),           & ! External weights for RK for mixt. fract. diff.
      fact_diag_low_RK_CMC_chm(:),         & ! Mesh factor for lower diagonal for RK for mixt. fract. diff.
      fact_diag_RK_CMC_chm(:),             & ! Mesh factor for diagonal for RK for mixt. fract. diff.
      fact_diag_up_RK_CMC_chm(:),          & ! Mesh factor for upper diagonal for RK for mixt. fract. diff.
      alpha_grid_CMC_chm(:),               & ! Vector of alfa in which save homogeneous reactor solutions (n_alpha_grid_CMC_chm)
      homog_reactors_CMC_chm(:,:,:),       & ! Solutions for homogeneous reactors. Used for initialization and BCs (nZ_CMC_chm,
                                             !  n_alpha_grid_CMC_chm, nclas_chm)
      alpha_val_spec_CMC_chm(:,:,:),       & ! alfa values for each species and point (npoin,nclass_chm). Used for BCs
      Text_cond_CMC_chm(:),                & ! Minimum and maximum conditional temperatures
      sum_Yk_ext_cond_CMC_chm(:),          & ! Minimum and maximum conditional mass fraction sums
      normal_CMC_mesh_CMC_chm(:),          & ! Normalization for the volumetric integrals for CMC mesh for interpolations (npoin)
      aux_cond_fields_CMC_chm(:,:),        & ! Auxiliar conditioned fields for CMC mesh for interpolations
                                             !  (2+2_ndime,gp_total_CMC_chm)
      aux_interp_fields_CMC_chm(:,:),      & ! Auxiliar interpolated fields for CMC mesh for interpolation (npoin,4+2*ndime) or
                                             !  (npoin,2+2*ndime)
      max_diff_Yk_inert_eq_CMC_chm(:),     & ! Maxima difference values for species mass fractions between inert and equilibrium
      aux_val_CMC_chm(:),                  & ! Auxiliar matrix for different purposes (npoin)
      matr_scal_CMC_chm(:),                & ! Used as a pointer to matrices
      ! Temporal averaged fields
      av_Yk_CMC_chm(:,:,:),                & ! Temp. averaged conditional mass fractions for CMC Y_k at nodes (nZ,npoin,nclas)
      av_enthalp_CMC_chm(:,:),             & ! Temp. averaged conditional enthalpy for CMC h at nodes (nZ,npoin)
      av_temp_CMC_chm(:,:),                & ! Temp. averaged conditional temperature for CMC at nodes (nZ,npoin)
      av_Yk_int_CMC_chm(:,:),              & ! Temp. averaged unconditional mass fractions for species (and soot) Y_k at
                                             !  nodes (npoin,nclas)
      av_enthalp_int_CMC_chm(:),           & ! Temp. averaged unconditional enthalpy for CMC Y_k at nodes (npoin)
      av_temp_int_CMC_chm(:),              & ! Temp. averaged unconditional temperature for CMC Y_k at nodes (npoin)
      av_visco_lam_int_CMC_chm(:),         & ! Temp. averaged integrated laminar dynamic viscosity from CMC (npoin)
      t_chem_integ_CMC_chm(:,:)              ! Time spent in chemical integration (nZ,npoin)


  integer(ip)                           :: &
      kfl_weigh_in_eq_CMC_chm,             & ! When starting from scratch activate the option to find the initial solution from
                                             !  weighing the inert and equilibrium solutions
      kfl_split_CFD_CMC,                   & ! Flag to split CFD and CMC into two different executions: 0 CFD and CMC in same
                                             !  execution and 1 if CFD and CMC separated (MODIFY: THIS SHOULD BE IN THE KERNEL).
                                             !  If 1 make turbulent diffusion matrix points to nu_t
      kfl_solve_enth_CMC_chm,              & ! 0 if enthalpy is not solved and 1 if it is transported
      kfl_start_CMC_chm = 1,               & ! 1 for the initial time step; 0 otherwise
      kfl_solve_cond_CMC_chm,              & ! 0: equations for mixture fraction, its variance, etc. are solved (unconditional
                                             !  variables); 1: equations for conditional species (and enthalpy) are solved
      kfl_bc_init_method_CMC_chm,          & ! Method for the calculation of boundary conditions and initial conditions; 1 for
                                             !  linear relationships in alfa, 2 for homogeneous reactors
      kfl_trans_phs_spc_CMC_chm,           & ! 0: do not include transport in physical space in CMC equations, 1: include it
      kfl_trans_mxt_spc_CMC_chm,           & ! 0: do not include transport in mixture fraction space in CMC equations, 1: include it
      kfl_mesh_interp_CMC_chm,             & ! Flag for different used meshes in CFD and CMC. 0: same meshes, 1: different meshes
      kfl_avg_cond_CMC_chm,                & ! Flag to indicate that averages to go from CFD to CMC are made with conditional or
                                             !  unconditional values; 0: unconditional fields, 1: conditional fields
      kfl_post_gp_CMC_chm,                 & ! Where to postprocess variables: 0 on nodes, 1 on Gauss points
      kfl_incl_PDF_trans_CMC_chm,          & ! 0: do not include PDF in CMC transport equations, 1: include it
      kfl_dt_calc_CMC_chm,                 & ! 0: omit the calculation of dt if dt is prescribed, 1: otherwise
      kfl_bc_alpha_CMC_chm,                & ! 0: boundary conditions defined through species; 1: boundary conditions defined
                                             !  through alpha
      kfl_transfer_condField_CMC_chm,      & ! 0: do not transfer conditional fields, 1: transfer conditional fields
      kfl_av_species_CMC_chm,              & ! 0: do not post-process averages for Yk in CFD for CMC model; 1: otherwise
      kfl_av_enthalp_CMC_chm,              & ! 0: do not post-process averages for enthalpy in CFD for CMC model; 1: otherwise
      nZ_CMC_chm,                          & ! Number of slices in mixture fraction space for CMC conditioning. Due to different
                                             !  reasons when using finite rate we take nZ_CMC_chm = 3
      nZ_AMC_CMC_chm,                      & ! Number of mixture fractions for scalar dissipation rate integration (AMC model)
      nS_AMC_CMC_chm,                      & ! Number of segregation factors for scalar dissipation rate integration (AMC model)
      nsize_mf_vec_chm,                    & ! Length of the mixture fraction path
      pos_Zstq_CMC_chm,                    & ! Position for the stoichiometric mixture fraction
      imixf_rk,                            & ! Mixture fraction iterator in RK scheme
      nvar_CMC_chm,                        & ! Number of variables to be solved in CMC: nvar_CMC_chm = nclas_chm(+1) (species +
                                             !  soot + enthalpy)
      nvar_therm_CMC_chm,                  & ! Number of thermochemical variables to be solved in CMC: species (+ enthalpy)
      index_N2,                            & ! Index in the mechanism for N2
      posZ_chem_integr_CMC_chm(2),         & ! Positions in the mixture fraction vector between which do chemical integrations
      order_RK_CMC_chm = 3_ip,             & ! Order for the Runge-Kutta applied to mixture fraction diffusion
      n_alpha_grid_CMC_chm = 51_ip,        & ! Number of values in which to save homogeneous reactor solutions along their chemical
                                             !  evolution
      nspec_uncond_write_CMC_chm,          & ! Number of <Y> to be written
      nspec_cond_write_CMC_chm,            & ! Number of <Y|Z> to be written
      nZ_write_CMC_chm,                    & ! Number of mixture fractions to be written
      nZ_chm_int_CMC_chm,                  & ! Number of mixture fractions to do chemical integration
      nZ_no_chm_int_CMC_chm,               & ! Number of mixture fractions to do not do chemical integration
      size_tncod_CMC_chm,                  & ! Number of fields for boundary conditions
      nspec_transf_CMC_chm,                & ! Number of species to be transferred from CFD to CMC
      transf_entha_CMC_chm,                & ! Transfer enthalpy
      nvar_trans_CMC_chm,                  & ! Total number of variables to be transferred
      nwrite_sol_CFD_CMC_chm,              & ! How often solutions for CFD are written
      freq_CFD_coup_CMC_chm,               & ! Frequency at which CFD transfers information in the coupling
      last_time_step_CFD_CMC_chm,          & ! Last time step from CFD
      last_time_step_CMC_CMC_chm             ! Last time step from CMC


  real(rp)                              :: &
      Zs_CMC_chm,                          & ! Saturation mixture fraction
      Zstq_CMC_chm,                        & ! Stoichiometric mixture fraction
      Smax_AMC_CMC_chm,                    & ! Maximum segregation factor for AMC model
      S_threshold,                         & ! Threshold segregation factor for integrations
      extr_Z_chem_integr_CMC_chm(2),       & ! Mixture fractions limits for which chemical integration is done
      extremes_Z_CMC_chm(2),               & ! Mixture fraction extremes to discern when to apply models close to Z=0 and Z=Zs
      fact_BC_RK_CMC_chm(2),               & ! Mesh factor for BCs for RK for mixt. fract. diff.
      exp_alpha_grid_CMC_chm = 2.0_rp,     & ! Exponent for the distribution of alpha_grid_CMC_chm
      inv_exp_alpha_grid_CMC_chm,          & ! 1/exp_alpha_grid_CMC_chm
      Zavg_const_CMC_chm,                  & ! Constant value for Zavg
      Zvar_const_CMC_chm,                  & ! Constant value for Zvar
      Xtot_const_CMC_chm,                  & ! Constant value for scalar dissipation rate
      Zlim_clipZvar_CMC_chm(2),            & ! Mixture fraction limits for Zvariance clipping
      Slim_clipZvar_CMC_chm(2),            & ! Segregation factor limits for Zvariance clipping
      Text_uncond_CMC_chm(2),              & ! Minimum and maximum unconditional temperature
      dt_CFD_CMC_chm                         ! Time step in the CFD


  type(r3p),pointer                     :: &
      ! Following variables save conditional values
      condu_gp_CMC_chm(:),                & ! Thermal conductivity at Gauss points (nelem)%(pgaus,nZ,1)
      sphec_gp_CMC_chm(:),                & ! Specific heat at Gauss points (nelem)%(pgaus,nZ,1)
      visco_gp_CMC_chm(:),                & ! Viscosity at Gauss points (nelem)%(pgaus,nZ,1)
      spvol_gp_CMC_chm(:),                & ! Specific volume at Gauss points (nelem)%(pgaus,nZ,1)
      react_scalars_gp_CMC_chm(:)           ! Reacting scalars in case we do smoothing (nelem)%(2*nclas_chm+6,pgaus,1)


  character(len=:), allocatable         :: &
       mf_vec_path_CMC_chm                   ! Mixture fraction vector path for CMC

  !
  ! End of CMC turbulent combustion model variables
  !


  character(len=:), allocatable ::         &
       mechanism_path,                     & ! Mechanism name for Finite Rate chemistry
       Red_spec                              ! Species to be reduced

  integer(ip)                           :: &
       nsize_mech_name,                    & ! Number of characters mechanism name
       nsize_red                             ! Number of characters of reduced Species
  !
  ! Global variables
  !
  real(rp),  pointer                    :: &
       rspec_chm(:,:),                     & ! Flamelet model: species mass fraction for radiation model
       zgradmax_chm(:),                    & ! Maximum gradient of the mixture fraction that defines phi_chm=1
       phi_chm(:),                         & ! Weighting factor for the hybrid model
       dt_rho_chm(:),                      & ! Projection of dt/rho
       dt_chm(:),                          & ! Projection of dt for interface tracking in ELSA model
       avden_chm(:),                       & ! Time-averaged density
       av_Z_flux_chm(:,:),                 & ! Time-averaged mixture fraction flux
       av_Named_Unk_chm(:,:),              & ! Time-averaged unknowns mixed equation model
       avY_chm(:),                         & ! Time-averaged reaction progress variable Yc or C
       avY2_chm(:),                        & ! Time-averaged reaction progress variable squared Yc*Yc or C*C
       avZv_chm(:),                        & ! Time-averaged variance of Z
       avZ_chm(:),                         & ! Time-averaged mixture fraction Z
       avYv_chm(:),                        & ! Time-averaged variance reaction progress Yc
       avmsk_chm(:),                       & ! Time-averaged mass source term from spray
       xYr_chm(:),                         & ! Scalar dissipation rate of Yc (resolved part)
       xZr_chm(:),                         & ! Scalar dissipation rate of Z  (resolved part)
       xYs_chm(:),                         & ! Scalar dissipation rate of Yc (subgrid part)
       xZs_chm(:),                         & ! Scalar dissipation rate of Z  (subgrid part)
       avxYr_chm(:),                       & ! Average scalar dissipation rate of Yc (resolved part)
       avxZr_chm(:),                       & ! Average scalar dissipation rate of Z  (resolved part)
       avxYs_chm(:),                       & ! Average scalar dissipation rate of Yc (subgrid part)
       avxZs_chm(:),                       & ! Average scalar dissipation rate of Z  (subgrid part)
       avZ2_chm(:),                        & ! Time-averaged squared of mixture fraction Z*Z
       avposttab_chm(:,:)                    ! Time-averaged tabulated postprocessing quantities

  !
  ! Spray variables
  !
  real(rp),  pointer                    :: &
       avL_chm(:),                         & ! Time-averaged liquid volume fraction phi_L
       avL2_chm(:),                        & ! Time-averaged liquid volume fraction squared phi_L*phi_L
       Sigma_chm(:),                       & ! Interface surface density Sigma
       Sigm0_chm(:),                       & ! Interface surface density Sigma_0 or Sigma_min
       d32_chm(:),                         & ! Sauter mean diameter
       avS_chm(:),                         & ! Time-averaged interface surface density Sigma
       avS0_chm(:),                        & ! Time-averaged interface surface density Sigma_0 or Sigma_min
       avd32_chm(:)                          ! Time-averaged Sauter mean diameter

  type(r3p),  pointer                   :: &
       sigma_gp_chm(:),                    & ! Interface surface density Sigma
       d32_gp_chm(:),                      & ! Sauter mean diameter
       dummy_enthalpy(:),                  & ! enthalpy at gauss points (sensible + chemical)
       sigma0_gp_chm(:)                      ! Interface surface density Sigma_0 or Sigma_min
  !
  ! Level set variables for spray
  !
  real(rp),  pointer                    :: &
       lap_phi_levSet_chm(:),              & ! Laplacian of the gradient of level set function
       grad_phic_levSet_chm(:,:),          & ! Gradient of phi*(1-phi)
       grad_phi_levSet_chm(:,:)              ! Gradient of level set function phi

  !
  ! Sectional soot model
  !
  integer(ip)                           :: &
       kfl_soot_chm,                       & ! Activation soot model
       kfl_yk_fw_ssm                         ! Index of lookup framework for soot precursor species
!!DMM       gasCoupling_ssm,                    & ! Coupling with gas phase
!!DMM       ID_NUCL,                            & ! Index nucleation process
!!DMM       ID_COND,                            & ! Index condensation process
!!DMM       ID_COAG,                            & ! Index coagulation process
!!DMM       ID_SURF,                            & ! Index surface growth process
!!DMM       nsect_chm                             ! Number of sections for soot sectional model

  !
  ! Others
  !
  real(rp)                              :: &
       radwt_chm,                          & ! Wall temprature for radiation model
       dac_crit_chm,                       & ! Flag for crit value in PFA
       dac_cor_chm,                        & ! Corrleation thershold for CODAC
       bf_fuel_chm,                        & ! Bilgers Formula Fuel
       bo_oxy_chm,                         & ! Bilgers Formula Oxy
       sorad_chm,                          & ! Source radius
       socen_chm(3),                       & ! Source center
       prthe_chm,                          & ! Thermodynamic pressure (if NASTIN not activated)
       chm_zmax,                           & ! zmax limit for flame index
       chm_zmin,                           & ! zmin limit for flame index
       surf_tension_chm,                   & ! Surface tension for sprays
       hrr_int_chm                           ! Total heat release in the domain (W)

  integer(ip)                           :: &
       nreac_chm,                          & ! Number of reactions
       stofu_chm(3)                          ! Which is fuel and which is oxygen for equivalence ratio

  !
  ! Variables transmitted to nodes later
  !
  integer(ip), pointer                  :: &
       React_ind(:,:),                     & ! Local reduced mechanims based on PFA
       Field_ind_chm(:),                   & ! Index of field variable
       Corr_chm(:)                           ! Local correlation for dynamic reduction

  real(rp),  pointer                    :: &
       diffu_chm(:,:)                        ! Diffusion constants

  real(rp),  pointer                    :: &
       entha_chm(:,:),                     & ! Enthalpy (sensible + chemical) for each species
       elem_h(:),                          & ! Elemental mass fraction of H
       elem_c(:),                          & ! Elemental mass fraction of C
       elem_o(:),                          & ! Elemental mass fraction of O
       elem_n(:),                          & ! Elemental mass fraction of N
       mixfr_chm(:),                       & ! Mixture Fraction Chemic (finite rate)
       prog_var_chm(:),                    & ! Progress Variable TDAC
       sum_reac_chm(:),                    & ! Sum of Reactions at a node
       src_chm(:,:)                          ! Chemical Source - Finite Rate

  integer,parameter                     :: &
       maxsp_chm = 1000                       ! Max number of species tracked in combustion code

  !
  ! Eulerian droplet identification variables
  !
  integer(ip)                           :: &
       kfl_droplet_id_chm,                 & ! Droplet identification flag
       ndrop_chm,                          & ! Total number of identified Eulerian droplets
       droplet_postprocess_frequency_chm     ! Droplet postprocess frequency

  real(rp)                              :: &
       levelSet_threshold_chm,             & ! Level Set threshold for droplet identification
       droplet_compactness_limit_chm,      & ! Compactness value below which a cluster won't be considered as a droplet
       droplet_max_diameter_chm,           & ! Max. diameter above which a cluster won't be considered as a droplet
       droplet_h_factor_chm                  ! Mesh size factor to define a max. diameter

  real(rp),  pointer                    :: &
       volume_drop_chm(:),                 & ! Droplet volume
       diameter_drop_chm(:),               & ! Droplet diameter
       centroid_drop_chm(:,:),             & ! Droplet centroid
       compactness2_drop_chm(:),           & ! Droplet compactness_2
       volume_cluster_chm(:)                 ! Cluster volume

  !------------------------------------------------------------------------
  ! Numerical problem: read in chm_reanut
  !------------------------------------------------------------------------

  integer(ip)                           :: &
       kfl_ellen_chm,                      & ! =0,1 for min/max element length
       kfl_taust_chm,                      & ! Tau calculation option
       kfl_shock_chm,                      & ! Shock capturing type
       kfl_stabi_chm,                      & ! Stabilization strategy
       kfl_limit_chm,                      & ! Limiter
       kfl_tiacc_chm,                      & ! Temporal accuracy
       kfl_tibub_chm,                      & ! Time integration of bubble
       kfl_split_chm,                      & ! Splitting algorithm
       kfl_tisch_chm,                      & ! Time integration scheme
       kfl_normc_chm,                      & ! Norm of convergence
       kfl_dtcri_chm,                      & ! dt criteria
       kfl_negat_chm,                      & ! Startegy for negative concentrations
       kfl_posit_chm,                      & ! Startegy for too positive concentrations
       kfl_warni_chm,                      & ! Warn about points with zero sum mass
       kfl_temli_chm                         ! Flag to activate a T limiter to compute reaction rates

  real(rp)                              :: &
       staco_chm(3),                       & ! Stability constants
       shock_chm,                          & ! Shock capturing parameter
       bemol_chm,                          & ! Bemol
       temli_chm,                          & ! Temperature limiter to compute reaction rates
       cotol_chm,                          & ! Convergence tolerance
       safet_chm,                          & ! Safety factor for time step
       chemical_time_factor,               & ! Safety factor exclusively for the source term
       cutof_chm,                          & ! Concentration cutoff for critical time computation
       sstol_chm,                          & ! Steady state tolerance
       strec_chm,                          & ! Adaptive dt: Stretching factor
       dampi_chm,                          & ! Adaptive dt: damping
       epsht_chm,                          & ! Adaptive dt: eps_R
       epstr_chm,                          & ! Adaptive dt: eps_A
       dtmin_chm,                          & ! Minimum time step
       dtmax_chm,                          & ! Maximum time step
       relax_chm                             ! Relaxation of update

  !------------------------------------------------------------------------
  ! Output and Postprocess: read in chm_reaous
  !------------------------------------------------------------------------

  integer(ip)                           :: &
       kfl_sized_chm                         ! Size distribution
  integer(ip)                           :: &
       ipara_chm(npara_chm)                  ! Int parameters for sets
  real(rp)                              :: &
       rpara_chm(npara_chm),               & ! Real parameters for sets
       avtim_chm                             ! Accumulated time for time-averaging

  !------------------------------------------------------------------------
  ! Boundary conditions: read in chm_reabcs
  !------------------------------------------------------------------------

  integer(ip)                           :: &
       kfl_conbc_chm,                      & ! Constant boundary conditions
       kfl_allcl_chm,                      & ! Bc on all classes
       kfl_fields_scale_chm                  ! Initialization by unscale (== 1) or scale (== 0) values
  integer(ip),   pointer                :: &
       kfl_initi_chm(:)                      ! Initial condition
  real(rp),      pointer                :: &
       xinit_chm(:,:)                        ! Initial pvalue parameter
  type(bc_nodes), pointer               :: &
       tncod_chm(:)                          ! Node code type
  type(bc_nodes), pointer               :: &
       tgcod_chm(:)                          ! Geometrical node code type
  type(bc_bound), pointer               :: &
       tbcod_chm(:)                          ! Boundary code type
!--END REA GROUP

  !
  ! Others
  !
  integer(ip),   pointer                :: &
       kfl_fixno_chm(:,:),                 & ! Nodal fixity
       kfl_fixbo_chm(:,:),                 & ! Boundary fixity
       kfl_funno_chm(:,:),                 & ! Function # for node BC
       kfl_funtn_chm(:,:)                    ! Function type for node BC

  real(rp),      pointer                :: &
       bvess_chm(:,:)                      ! Essential bc values

  integer(ip)                           :: &
       kfl_robin_chm,                      & ! Robin condition exists
       ncomp_chm,                          & ! Number of components
       kfl_grdif_chm,                      & ! If there are gradients of conductivity
       kfl_tiaor_chm,                      & ! Original time accuracy
       kfl_goite_chm,                      & ! Keep iterating
       iclas_chm,                          & ! Current class being solved
       iclai_chm,                          & ! Initial class
       iclaf_chm,                          & ! Final class
       ittot_chm,                          & ! Total number of iterations
       nskyl_chm,                          & ! Size of the skyline matrix
       kfl_goit2_chm,                      & ! Internal goite
       kfl_under_chm,                      & ! # undershoots
       kfl_overs_chm                         ! # overshoots

  real(rp)                              :: &
       dtinv_chm,                          & ! 1/dt
       dtcri_chm,                          & ! Critical time step
       resid_chm,                          & ! Residual for outer iterations
       pabdf_chm(10),                      & ! BDF factors
       rtpts_chm,                          & ! Global inner residual
       comin_chm,                          & ! Minimum concentration
       comax_chm,                          & ! Maximum concentration
       cputi_chm(10),                      & ! CPU time
       dtmat_chm,                          & ! Matrix time step
       xvael_chm(100)                        ! Values

  real(rp),     pointer                 :: &
       amatr_chm(:),                       & ! Matrix for ODE's
       rhsid_chm(:),                       & ! RHS for ODE's
       ripts_chm(:),                       & ! Class inner residual
       vmass_chm(:),                       & ! Mass matrix
       smatr_chm(:),                       & ! Constant matrix
       shsid_chm(:)                          ! Constant RHS

  integer(ip),  pointer                 :: &
       iarea_chm(:),                       & ! IA CSR format
       jarea_chm(:),                       & ! JA CSR format
       iskyl_chm(:),                       & ! DOE's skyline list
       idiag_chm(:),                       & ! DOE's skyline diagonal
       idima_chm(:)                          ! Position of diagonal in sparse matrix

  type(ADR_typ), allocatable, target    :: &
       ADR_chm(:)                            ! ADR type for nclas_chm

#ifdef CANTERA
  type(phase_t) gas_chm
#endif


end module def_chemic
