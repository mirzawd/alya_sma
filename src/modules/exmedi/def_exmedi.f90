!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    def_exmedi.f90
!> @date    12/03/2013
!> @author  Mariano Vazquez
!> @brief   Def file
!> @details Def file
!> @}
!------------------------------------------------------------------------
module def_exmedi

!-----------------------------------------------------------------------
!    
! Heading for the module subroutines
!
!-----------------------------------------------------------------------
  use def_kintyp
  use def_master
!
! General
!------------------------------------------------------------------------
! Parameters
!------------------------------------------------------------------------
  real(rp), parameter :: &
       zeexm = epsilon(0.0_rp)
  integer(ip), parameter ::&                     ! # postprocess variables
       lun_vinte_exm       = 512, &
       nvarp_exm           = 60_ip, &
       nvars_exm           = 20_ip, &
       nvecp_exm           = 20_ip, &
       nspep_exm           = 35_ip, &
       nscap_exm           = 60_ip, &
       nvars_ttparmate_exm = 17_ip, &
       msets_exm           = 10_ip, &              ! # max sets
       ohara_conductance_ikatp_row = 14_ip       !Row number for Ik_atp in Ohara conductance table TODO: move somewhere else


!------------------------------------------------------------------------
! Physical problem: read, modified or defined in exm_reaphy FIRST CALL
!------------------------------------------------------------------------
!--BEGIN REA GROUP
  integer(ip), pointer :: &
       kfl_voini_exm(:),&          ! Initial voltage flag (per material)
       kfl_user_specified_celltypes_exm(:,:), &          ! Flag to know whether to run initial condition ODEs for different cell types. Each element 0|1, for ENDO, MID, EPI.  1 by default
       kfl_ignore_steadystate_celltypes_exm(:,:), &      ! Flag to know whether to die if different cell types did not reach steady state. Each element 1(ignore)|0(terminate), for ENDO, MID, EPI. 0 by deafult
       kfl_steadystate_variable(:), &     ! Which variable will be used to test for the steady state. See EXM_CELL_STEADY_* variables
       moneclmate_exm(:,:),&        ! Number of beats to steady state and Cycle length
       kfl_hfmodmate_exm(:),&       ! Runs the model on Heart Failure mode
       !kfl_hfmod_exm(:),&           ! Runs the model on MODIFIED mode
       kfl_fract_diffusion_exm(:),& ! Fractional diffusion activated (per material)
       kfl_isac_exm(:)              ! Flag that indicates currents of ISAC (stretch-activated-channel current)  

  logical(lg), pointer :: &
       kfl_active_material(:)       ! Flag to indicate if elements with this material are used in the matrix assembly or not
  logical(lg)          :: &
       kfl_reset_fisoc              ! Flag to check if fisoc was postprocessed and needs reset at the end of the timestep 

  integer(ip) ::&
       !kfl_fento_exm,&             ! Specific sub model of fenton's model
       kfl_appty_exm,&             ! Type of applied currents
       kfl_appva_exm,&             ! What is applied: current, density current, voltage...
       kfl_cemod_exm,&             ! Cell model (propagation model): monodomain o bidomain
       kfl_stree_exm,&             ! Construct a streeter fiber model
       kfl_nodif_exm,&             ! when set to 1, do not compute diffusion terms
       kfl_paced_exm,&             ! Flag for pacing or not the hiPSC-CMs
       kfl_save_convergence_cellmodel,& !for ohara model: 0 - do not save, 1 - save files with sequences of currents and stuff
       kfl_save_init_cellmodel      ! dump initial conditions to hardcode cell model 0 - do not, 1 - dump

  integer(ip) ::         &
       ndofn_exm,      &         ! # of d.o.f. of the problem
       ndof2_exm,      &         ! ndofn_*ndofn_
       nmate_exm,      &         ! # of materials
       nevat_exm,      &         ! Element matrix dim.=(ndime+1)*nnode
       nstim_exm,      &         ! Number of initial stimuli
       nstis_exm,      &         ! If needed, set number of initial stimuli
       modst_exm,      &         ! Stimuli field
       modab_exm,      &         ! Apex to base gradient field
       ncomp_exm,      &         ! Number of components of the fields
       nstrb_exm,      &         ! Number of sets to generate streeter fiber field 
       strbo_exm(3),&              ! Boundaries for the streeter model
       nrootecg_exm,      &         ! Number of electrodes to be considered in the pseudo_ecg
       kcopeecg_exm                 ! Pseudo-ecg counter

  !---------------------------------
  ! ECG arrays          
  !
  type ecg_coords
     real(rp)         :: coords(3)       ! coordinates of the points to measure ECG potential
     character(5)     :: label           ! string label for the columns in .vin, does not have to be unique, but preferable
  end type ecg_coords

  type(ecg_coords), dimension(:), pointer :: ecg_points        => null() ! ECG points
  real(rp),                       pointer :: pseudecg_exm(:)   => null() ! ECG values for the pseudo-ecg
  !
  ! ECG arrays          
  !---------------------------------


  real(rp), pointer :: &
     xmccmmate_exm(:)    ,&        ! Chi membrane capacitance
     gdiff_exm(:,:,:),&            ! Global diffusivities
     ttparmate_exm(:,:,:),&        ! parameters to identify normal vs heart failure 1-cell simulations
     vminimate_exm(:,:), &         ! new voltage(celltype,imate)
     vauxi_exm_initial(:,:,:),&    ! initial values of vauxi for 3D simulation (nvars, ncelltypes_ecc, nmaterials)
     vconc_initial(:,:,:),&        ! initial vaules of vconc for 3D simulation (nvars, ncelltypes_ecc, nmaterials)
     fract_diff_coef_exm(:), &     ! Fractional diffusion coefficient S
     steady_tol_cellm(:), &        ! Tolerance to use to determine the steady state in the 0d cell model. Leave -1 by default, then subroutine will determine the correct one.
     timestep_cellm(:)             ! Timestep for the cellmodel. Default = -1, meaning the cell model decides

  real(rp) ::&
       dtinv_exm    ,&              ! 1/dt
       fiaxe_exm(3),&               ! Ventricular axis
       stran_endo_exm,&             ! Angle for streeter endo
       stran_epi_exm,&              ! Angle for streeter epi
       volcai_exm                   ! Volume integral of calcium
   
  real(rp), dimension(:), pointer :: &
       apval_exm => null(),&         ! Applied current intensity  
       aplap_exm => null(),&         ! Applied current time lapse 
       aprea_exm => null(),&         ! Reach 
       aptim     => null()           ! Activation potential starting times  EXMEDI - SOLIDZ

  real(rp), dimension(:,:), pointer :: &
       apcen_exm => null()           ! Applied current center 


  integer(ip), pointer ::&
       fract_diff_nintp_exm(:)          ! Fractional diffusion number of integration points

  integer(ip) ::&
       nauxi_exm,&                  ! # of auxiliary variables (vauxi_exm_initial)
       nconc_exm,&                  ! # of concentration unknowns (vconc_initial)
       nvint_exm,&                  ! Volume integral postprocess frequency (time steps)
       nicel_exm                    ! # of cell ionic currents

  integer(ip), pointer :: &         ! TODO: number of celltypes should be a parameter in .dat but that will require refactoring ohara
       ncelltypes_per_model(:)      ! # of celltypes (3 for dffect, 9 for courtemanche)

  logical(lg), pointer ::&
       kfl_ortho_diffusion_exm(:)   ! Orthotropic model for the diffusivity tensor
  
!------------------------------------------------------------------------
! Numerical problem: read, modified or defined in exm_reanut
!------------------------------------------------------------------------


  integer(ip) ::&
       kfl_genal_exm,&              ! General algorithm type ---DEPRECATED---
       kfl_timet_exm,&              ! Time treatment
       kfl_tiacc_exm,&              ! Temporal accuracy
       kfl_ticel_exm,&              ! Cell model time advance method
       kfl_tisch_exm,&              ! Temporal accuracy
       kfl_goite_exm,&              ! Keep iterating
       kfl_shock_exm,&              ! Shock capturing type 
       kfl_comat_exm,&              ! Compute amatr
       kfl_weigh_exm,&              ! Weighting of dT/dt
       kfl_adres_exm,&              ! Subiterations adaptive flag.
       kfl_normc_exm,&              ! Norm of convergence
       kfl_algso_exm,&              ! Type of algebraic solver
       kfl_repro_exm,&              ! Stabilization based on residual projection
       kfl_nolim_exm,&              ! Non-linear correction method
       kfl_nolum_exm,&              ! Non-linear terms lumped or not
       kfl_gcoup_exm,&              ! Geometric coupling with SOLIDZ
       miinn_exm,&                  ! Max # of iterations
       msste_exm,&                  ! Time substepping number
       last_iters_exm,&             ! Solver iterations for the current sub-iteration
       mnoli_exm,&                  ! Max # of iterations for the non-linear intracellular problem
       msoit_exm,&                  ! Max # of solver iterations
       nkryd_exm,&                  ! Krylov dimension
       memor_exm(2),&               ! Memory counter
       itera_exm,&                  ! Internal iteration counter
       nunkn_exm

  real(rp) ::&
       dtcri_exm    ,&              ! Critical time step
       shock_exm    ,&              ! Shock capturing parameter
       sstol_exm    ,&              ! Steady state tolerance
       cotol_exm    ,&              ! Convergence tolerance
       corat_exm    ,&              ! Convergence tolerance ratio
       dtext_exm    ,&              ! Externally fixed time step
       safet_exm    ,&              ! Safety factor for time step
       solco_exm    ,&              ! Solver tolerance
       weigh_exm,    &              ! Weight of dU/dt in the residual
       tnoli_exm    ,&              ! Tolerance for the non-linear intracellular problem
       resid_exm(10),&              ! Residual for outer iterations
       resou_exm(10) ,&              ! Reference residual for outer iterations 
       resin_exm(10) ,&              ! Reference residual for inner (sub)iterations 
       resin_first_exm(10) ,&        ! Reference residual for first inner (sub)iterations 
       err01_exm(2),&               ! L1 error T
       err02_exm(2),&               ! L2 error T
       err0i_exm(2),&               ! Linf error T
       err11_exm(2),&               ! L1 error grad(T)
       err12_exm(2),&               ! L2 error grad(T)
       err1i_exm(2),&               ! Linf error grad(T)
       staco_exm(3),&               ! Stability constants
       cpu_exmed(2),&               ! CPU for the EXM problem
       cpu_ass_sol_exm(4)           ! CPU time assembly and solver at each iteration

!------------------------------------------------------------------------
! Physical problem: read, modified or defined in exm_reaphy SECOND CALL
!------------------------------------------------------------------------

!
! Physical properties used in the model
!

  integer(ip),  pointer ::&
       idima_exm(:),      &              ! Diagonal indices for amatr
       kgrfi_exm(:)                      ! Large gradient fibers label vector
       
  real(rp),     pointer ::&
       cedif_exm(:,:,:)  ,&             ! Ex/Intracellular Diffusivity (point-wise)
       grafi_exm(:,:,:)  ,&             ! Fiber orientation gradient (a tensor)
       fiber_exm(:,:,:)  ,&             ! Fiber
       sheet_exm(:,:,:)  ,&             ! ortho fiber 1
       normal_exm(:,:,:) ,&             ! ortho fiber 2
       celty_exm(:)    ,&             ! Cell types 
       atbhe_exm(:,:)    ,&             ! Apex to base Conductance gradient 
       vdiag_exm(:)      ,&             ! Diagonal preconditioner 
!+MRV       
       fibe2_exm(:,:)    

!------------------------------------------------------------------------
! Boundary conditions: read, modified or defined  in exm_reabcs
!------------------------------------------------------------------------
!
! Boundary conditions
! 
  integer(ip)           ::   kfl_exboc_exm    ! Boundary conditions explicitly given
  type(bc_nodes), pointer              :: &     
       tncod_exm(:)                          ! Node code type
  type(bc_bound), pointer              :: &     
       tbcod_exm(:)                          ! Boundary code type

!------------------------------------------------------------------------
! Output and Postprocess: read, modified or defined  in exm_reaous
!------------------------------------------------------------------------

  real(rp)                 :: &
       fisoc_exm(2)              ! Isochrones trigger on depolarisation, >threshold, tuple (0/1, threshold)
              
!--END REA GROUP
!------------------------------------------------------------------------
! Derived variables (slave-local in MPI)
!------------------------------------------------------------------------
  logical                  :: &
       weparal
!
! Physical properties used in the model
!
  real(rp),     pointer ::&
       amatr_auxi_exm(:)            ! Auxiliary amatr

  real(rp),     pointer ::&
       appfi_exm(:) ,&              ! Applied current field 
       eflux_exm(:,:),&             ! Electrical flux 
       bipol_exm(:,:)               ! Bipolar leads electrical fluxes

  !
  ! Cell model types parameters and variables
  !
  real(rp), pointer ::&
       vauxi_exm(:,:,:),&           ! Auxiliary variables
       ticel_exm(:),&               ! Total cell ionic current
       jicel_exm(:),&               ! Jacobian of the total cell ionic current
       vicel_exm(:,:),&             ! Individual Cell ionic currents
       qneto_exm(:)

  integer(ip),    pointer ::&
       kwave_exm(:)                   ! Detect a new depolarization wave. 0 downstroke was detected, 1 upstroke was detected

  integer(ip) ::&
       kfl_refre_exm(10)              ! Compute the reference residual for the exm_cvgunk's itask  

  real(rp), parameter :: &
       sms_conversion_currents_exm = 1000.0_rp

  real(rp) :: &                     ! some local cputime measurement
       cpold_exm = 0.0_rp,& 
       cptot_exm = 0.0_rp


contains
!!!!! NOT USED AND CODE COEVARAGE DIES HERE
! function heavis(a,b)
!    use def_kintyp
!    implicit none
!    real(rp) :: a,b,heavis
!    
!    heavis = 1.0_rp
!    if (a<b) heavis = 0.0_rp
!    
!    
!  end function heavis
 

end module def_exmedi
