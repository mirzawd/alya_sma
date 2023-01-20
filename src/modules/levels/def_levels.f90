!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module def_levels 
  !-----------------------------------------------------------------------
  !    
  ! Heading for the level set advection equation routines
  !
  !-----------------------------------------------------------------------
  use def_kintyp
  !------------------------------------------------------------------------
  ! Parameters
  !------------------------------------------------------------------------
  integer(ip),   parameter       :: &
       lun_volum_lev = 1409,&
       lun_capte_lev = 1410, lun_reave_lev = 1411, lun_inter_lev = 1412,&
       lun_gauge_lev = 1413, lun_cored_lev = 1414
       
  character(150)                 :: &
       fil_rstar_lev
  real(rp),      parameter       :: &
       zelev = epsilon(1.0_rp)
  integer(ip),   parameter       :: &
       ngaug_lev=10,                & ! # number og gauges
       ncoef_lev=10,                & ! # coefficient for properties
       ninit_lev=10                   ! # coefficient for source term

  !------------------------------------------------------------------------
  ! Physical problem
  !------------------------------------------------------------------------
!--BEGIN REA GROUP
  integer(ip)                    :: &
       kfl_advec_lev,               & ! Existence of (u.grad)T
       kfl_reave_lev                  ! Initial velocity read from file 
  !
  ! Physical properties
  !
  integer(ip)                    :: &
       nmate_lev                      ! # of materials
       
  real(rp),     pointer          :: &
       flsex_lev(:)                   ! exact level set solution

  !------------------------------------------------------------------------
  ! Numerical problem
  !------------------------------------------------------------------------
  integer(ip)                    :: &
       kfl_timet_lev,               & ! Time treatment
       kfl_tisch_lev,               & ! Time integration scheme
       kfl_timco_lev,               & ! Temporal accuracy
       kfl_tiacc_lev,               & ! Temporal accuracy
       kfl_normc_lev,               & ! Norm of convergence
       neule_lev,                   & ! Number of Euler time steps 
       miinn_lev,                   & ! Max # of iterations
       inred_lev,                   & ! Initial Redistanciation
       nfred_lev,                   & ! Redistanciation frequency
       tyred_lev,                   & ! Redistanciation type
       kfl_geodi_lev,               & ! Geometrical distance type
       nstre_lev,                   & ! If /=0 indicates the step to stop redistancing, it is used to debug error with redistance 
       kfl_locre_lev,               & ! Local redistanciation, in a small layer (thicl*10). In the rest step over with the value prior to redistance
       nbitr_lev,                   & ! Redistanciation equation iteration number
       kfl_corvo_lev                  ! Volume correction trough Level Set
  real(rp)                       :: &
       safet_lev,                   & ! Safety factor for time step
       sstol_lev,                   & ! Steady state tolerance
       cotol_lev,                   & ! Convergence tolerance
       cpuit_lev,                   & ! CPU time per iteration
       supgp_lev                      ! SUPG stabililization

  !------------------------------------------------------------------------
  ! Boundary conditions
  !------------------------------------------------------------------------
  integer(ip)                    :: &
       kfl_inlev_lev,               & ! Initialisation of Level Set
       kfl_conbc_lev,               &
       kfl_initi_lev

  real(rp)                       :: &
       height_lev                     ! Height of level set
  type(bc_nodes), pointer        :: &     
       tncod_lev(:)                   ! Node code type
  type(bc_nodes), pointer        :: &     
       tgcod_lev(:)                   ! Geometrical node code type
  
  !------------------------------------------------------------------------
  ! Output and Postprocess
  !------------------------------------------------------------------------
  integer(ip)                    :: &
       npp_gauge_lev,               & ! Postprocess interface gauge(s)
       npp_nbgau_lev,               & ! Number of gauges to postprocess
       typga_lev(ngaug_lev),        & ! Type of gauge (direction)
       npp_inter_lev                  ! Interface mesh output
  real(rp)                       :: &
       cogau_lev(3,ngaug_lev)         ! Name and character of the postprocess variables
!--END REA GROUP 
  !------------------------------------------------------------------------
  ! Others
  !------------------------------------------------------------------------

  integer(ip), pointer           :: &
       kfl_fixno_lev(:,:),          &! Nodal fixity 
       kfl_funno_lev(:)
  
  real(rp),    pointer           :: &
       bvess_lev(:,:,:)               ! Essential bc values

  real(rp), pointer             :: &
       bvcod_lev(:)                   ! Node value arrays
       
  integer(ip)                    :: &
       ncomp_lev,                   & ! Number of components of the temperature
       kfl_ellen_lev,               & ! =0,1 for min/max element length
       kfl_zonal_lev,               & ! Redistanzation by zones (default all domain )
       kfl_tiaor_lev,               & ! Original time accuracy
       kfl_stead_lev,               & ! Steady-state has been reached 
       kfl_goite_lev,               & ! Keep iterating
       kfl_alloc_lev                  ! Redistantiation allocated
  real(rp)                       :: &
       dtcri_lev,                   & ! Critical time step
       dtinv_lev,                   & ! 1/dt
       pabdf_lev(10),               & ! BDF factors
       lemin_lev,                   & ! Minimum level set amplitude
       lemax_lev,                   & ! Maximum level set amplitude
       resid_lev                      ! Residual for outer iterations

  !------------------------------------------------------------------------
  ! Interface treatment
  !------------------------------------------------------------------------
  integer(ip), pointer           :: &
       lnodb_lev(:,:),              & ! connectivity 
       lnodp_lev(:,:),              & ! connectivity temporar
       lebsu_lev(:),                & ! element (whole mesh numeration) to which the surface belongs 
       lebsp_lev(:),                & ! element (whole mesh numeration) to which the surface belongs, temporar
       psbel_lev(:),                & ! pointer for lsbel 
       lsbel_lev(:),                & ! list of surfaces that belong to each element
       nredm_lev(:,:),              & ! redistanciation size for the master
       capin_lev(:)                   ! interface transducer
  real(rp), pointer              :: &
       coord_lev(:,:),              & ! interface points coordinates
       coorp_lev(:,:),              & ! interface points coordinates
       dista_lev(:),                & ! vector temporal to compute distance
       flev0_lev(:)                   ! phi0 for sussmann redistanciation
  integer(ip)                    :: &
       nboun_lev,                   & ! number of interface segment
       npoin_lev,                   & ! number of interface point  
       nelcr_lev,                   & ! number of element crossed by the interface
       ncapt_lev                   
  integer(ip), target            :: &
       nredi_lev(2),                & ! redistanciation size  (parallelization)
       nredt_lev(2),                & ! total redistanciation size (parallelisation)  
       findg_lev(ngaug_lev) 
  real(rp)                       :: &
       volrf_lev,                   & ! # Initial reference volume
       volit_lev,                   & ! # volume(t)
       dtred_lev                      ! # Time step for redistanciation
  real(rp),    target      :: &
       l1nor_lev(1),                & ! # L1 error norm
       valga_lev(ngaug_lev)           ! Value of Gauge
  real(rp), pointer              :: &
       norml_lev(:,:),              & ! # grad phi / |grad phi |
       normv_lev(:)                   ! # |v.grad phi / |grad phi ||
  integer(ip), pointer           :: &
       icupt_lev(:),                & ! point in element crossed by interface 
       elcro_lev(:)                   ! elements crossed by interface
  real(rp),    pointer           :: &
       walld_lev(:)

end module def_levels
