!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module def_master

  !-----------------------------------------------------------------------
  !****f* defmod/def_master
  ! NAME
  !   def_master
  ! DESCRIPTION
  !   This module is the header of the master
  !***
  !-----------------------------------------------------------------------
  
  use def_kintyp
  use def_kintyp_postprocess, only : typos
  use def_solver,             only : solve_sol,eigen_sol
  use mod_optional_argument,  only : optional_argument
  use def_kintyp_dims,        only : npoi1
  use def_kintyp_dims,        only : npoi2
  use def_kintyp_dims,        only : npoi3
  use def_kintyp_dims,        only : npoin
  use def_kintyp_dims,        only : nedg1
  use def_kintyp_dims,        only : nedg2
  use def_kintyp_dims,        only : nedg3

  implicit none
  
  !------------------------------------------------------------------------
  !
  ! Parameters
  !
  !------------------------------------------------------------------------

  integer(ip) :: new_periodicity
  
  integer(ip), parameter   :: &
       mmodu=30,              & ! Max. # of modules
       mblok=mmodu              ! Max. # of blocks
  !
  ! Task numbering
  !
  integer(ip), parameter   ::        &
       ITASK_REAPRO        =  1,     &
       ITASK_TURNON        =  2,     &
       ITASK_INIUNK        =  3,     &
       ITASK_TIMSTE        =  4,     &
       ITASK_BEGSTE        =  5,     &
       ITASK_DOITER        =  6,     &
       ITASK_CONCOU        =  7,     &
       ITASK_CONBLK        =  8,     &
       ITASK_ENDSTE        = 10,     &
       ITASK_OUTPUT        = 12,     &
       ITASK_TURNOF        = 13,     &
       ITASK_BEGITE        = 14,     &
       ITASK_ENDITE        = 15,     &
       ITASK_MATRIX        = 16,     &
       ITASK_DOOPTI        = 17,     &
       ITASK_ENDOPT        = 18,     &
       ITASK_BEGZON        = 19,     &
       ITASK_ENDZON        = 20,     &
       ITASK_SOLMEM        = 21,     &
       ITASK_REDIST        = 22,     &
       ITASK_INTERP        = 33,     &
       ITASK_BEGRUN        = 23,     &
       ITASK_READ_RESTART  = 24,     &
       ITASK_WRITE_RESTART = 25,     &
       ITASK_BEGINN        = 27,     &
       ITASK_ENDINN        = 28,     &
       ITASK_INNITE        = 29,     &
       ITASK_INITIA        = 30,     &
       ITASK_ENDTIM        = 31,     &
       ITASK_ENDRUN        = 32,     &
       ITASK_RECOVER_RESET = 33,     &
       ITASK_SAVE_RESET    = 34,     &
       ITASK_BEFORE        =  1,     &
       ITASK_AFTER         =  2
  !
  ! Task names
  !
  character(30)            :: TASK_LONG_NAME(31) = (/ &
       'READ PROBLEM DATA             ', & !  1
       'TURN ON                       ', & !  2
       'INITIAL UNKNOWNS              ', & !  3
       'COMPUTE TIME STEP             ', & !  4
       'BEGIN TIME STEP               ', & !  5
       'DO ITERATIONS                 ', & !  6
       'COUPLING CONVERGENCE          ', & !  7
       'COUPLING BLOCK                ', & !  8
       'NEW MESH                      ', & !  9
       'END TIME STEP                 ', & ! 10
       'FILTERING                     ', & ! 11
       'OUTPUT                        ', & ! 12
       'TURN OFF                      ', & ! 13
       'BEGIN NON-LINEAR ITERATIONS   ', & ! 14
       'END NON-LINEA ITERATIONS      ', & ! 15
       'COMPUTE MATRIX                ', & ! 16
       'BEGIN OPTIMIZATION            ', & ! 17
       'END OPTIMIZAITON              ', & ! 18
       'BEGIN ZONE                    ', & ! 19
       'END ZONE                      ', & ! 20
       'SOLVER MEMORY                 ', & ! 21
       'REDISTRIBUTE                  ', & ! 22
       'BEGIN RUN                     ', & ! 23
       'READ RESTART                  ', & ! 24
       'WRITE RESTART                 ', & ! 25
       'NULL                          ', & ! 26
       'BEGIN INNER ITERATIONS        ', & ! 27
       'END INNTER ITERATIONS         ', & ! 28
       'INSIDE INNTER ITERATIONS      ', & ! 29
       'INITIALIZE                    ', & ! 30
       'END TIME                      '  /)! 31
  !
  ! CPU
  !
  integer(ip), parameter    ::       &
       CPU_SMVP                = 32, &
       CPU_ASSEMBLY            = 31, &
       CPU_ASSEMBLY_BOUNDARY   = 33, &
       CPU_ASSEMBLY_NODE       = 41, &
       CPU_ASSEMBLY_PARTICLE   = 42, &
       CPU_SOLVER              = 28, &
       CPU_EIGEN_SOLVER        = 29, &
       CPU_TOTAL_MODULE        = 30, &
       CPU_MINI_ASSEMBLY       = 34, & 
       CPU_MAXI_ASSEMBLY       = 35, &
       CPU_COUNT_ASSEMBLY      = 36, &
       CPU_COUNT_BOUNDARY      = 38, &
       CPU_COUNT_NODE          = 39, &
       CPU_COUNT_PARTICLE      = 40, &
       CPU_DOITER              = 37, &
       CPU_BEGSTE              = 43, &
       CPU_ENDSTE              = 44, &
       CPU_REPART_LOAD         = 45, &
       CPU_REPART_WLR          = 46, &
       CPU_REPART_OLR          = 47, &
       CPU_BEGITE              = 48, &
       CPU_ENDITE              = 49, &
       CPU_READ_RESTART        = 51, &
       CPU_WRITE_RESTART       = 52, &
       CPU_OUTPUT              = 53, &
       CPU_INIUNK              = 54, &
       CPU_BEGRUN              = 55, &
       CPU_MINI_NODE           = 56, & 
       CPU_MAXI_NODE           = 57


  integer(ip), parameter    ::       &  ! Start operations: cpu_start
       CPU_READ_GEO            =  1, &
       CPU_MESH_PARTITION      =  2, &
       CPU_MESH_MULTIPLICATION =  3, &
       CPU_CONSTRUCT_DOMAIN    =  4, &
       CPU_ADDTIONAL_ARRAYS    =  5, &
       CPU_READ_SETS           =  7, &
       CPU_READ_BCS            =  8, &
       CPU_READ_FIELDS         =  9, &
       CPU_START_END           = 10

  integer(ip), parameter    ::       &  ! Domain operations: cpu_domain
       CPU_GROUPS              =  1, &
       CPU_HALOS               =  2, &
       CPU_OUTPUT_DOMAIN       =  3, &
       CPU_COUPLING            =  4, &
       CPU_ELSEST              =  5, &
       CPU_END_DOMAIN          = 10
  !
  ! When to solve
  !
  integer(ip), parameter   ::  &
       AT_EACH_TIME_STEP  = 0, &
       AT_BEGINNING       = 1, &
       AT_FIRST_TIME_STEP = 2
  !
  ! Module and service identifiers
  !
  integer(ip), parameter   :: &
       ID_KERNEL=  0,         & !
       ID_NASTIN=  1,         & !
       ID_TEMPER=  2,         & !
       ID_CODIRE=  3,         & !
       ID_TURBUL=  4,         & !
       ID_EXMEDI=  5,         & !
       ID_NASTAL=  6,         & !
       ID_ALEFOR=  7,         & !
       ID_LATBOL=  8,         & !
       ID_GUSANO=  9,         & !
       ID_SOLIDZ= 10,         & !
       ID_GOTITA= 11,         & !
       ID_WAVEQU= 12,         & !
       ID_LEVELS= 14,         & !
       ID_QUANTY= 15,         & !
       ID_MAGNET= 16,         & !
       ID_PARTIS= 17,         & !
       ID_NASEDG= 18,         & !
       ID_CHEMIC= 19,         & !
       ID_HELMOZ= 20,         & !
       ID_IMMBOU= 21,         & !
       ID_RADIAT= 22,         & !
       ID_CASIMI= 23,         & !
       ID_POROUS= 24,         & !
       ID_XXXXXX= 25,         & !
       ID_NEUTRO= 26,         & !
       ID_DODEME=  7,         & !
       ID_COUPLI= 25,         & !
       ID_OPTSOL= 13,         & !
       ID_INSITU= 27,         & !
       ID_SOLFE2= 28,         & !
       ID_KERMOD= mmodu
  !
  ! Primary variables identifiers
  !
  !TODO: these numbers do not coincide with numbers in the comments of the variables
  integer(ip), parameter   :: &
       ID_VELOC= 1,           & ! Velocity
       ID_PRESS= 2,           & ! Pressure
       ID_TEMPE= 3,           & ! Temperature
       ID_DENSI= 4,           & ! Density
       ID_ENERG= 5,           & ! Energy
       ID_VISCO= 6,           & ! Viscosity
       ID_UMOME= 7,           & ! Momentum
       ID_FMOME= 8,           & ! Fractional momentum
       ID_UNTUR= 9,           & ! Turbulence unknowns
       ID_UNCDR=10,           & ! CDR unknown
       ID_ELMAG=11,           & ! Electromagnetic
       ID_DISPM=12,           & ! Mesh displacement
       ID_VELOM=13,           & ! Mesh velocity
       ID_DISPL=14,           & ! Displacement
       ID_SPINS=15,           & ! Quantic spins
       ID_DISFU=16,           & ! Distribution funct.
       ID_VDROP=17,           & ! Droplet velocity
       ID_CDROP=18,           & ! Droplet concentration
       ID_VORTI=19,           & ! Vorticity
       ID_CONCE=20,           & ! Concentration
       ID_CDENS=24,           & ! Classical density
       ID_ENTHA=25,           & ! Enthalpy
       ID_PHION=26,           & ! Phi
       ID_FIBER=27,           & ! Fiber fields
       ID_GPFIB=31,           & ! Fiber fields in the deformed configuration (elementary)
       ID_RADSO=32,           & ! Heat radiation source term
       ID_VCONC=33,           & ! Ionic concentrations
       ID_WASAT=34,           & ! Water Saturation for Porous Flow
       ID_VEL1D=35,           & ! 1D velocity for network models
       ID_AREAS=36,           & ! Areas for network models
       ID_FLOWR=37              ! Flow rate for network models

  !
  ! Things to do
  !
  integer(ip), parameter   :: &
       WRITE_RESTART_FILE= 2, & ! Restart
       READ_RESTART_FILE=  1    ! Restart
  !
  ! Time step and iteration labels
  !
  integer(ip), parameter   :: &
       ITERAUX_EQ_TIMEN             =   1, & ! begste
       ITERK_EQ_ITERAUX             =   2, & ! begite
       ITERK_EQ_UNKNO               =   3, & ! endite(1)
       ITERAUX_EQ_ITERK             =   4, & ! endite(2)
       TIMEN_EQ_ITERK               =   5, & ! endste
       DT_PHYSICAL                  =   1, & !
       DT_PSEUDO                    =   2, & !
       ITER_K                       =   1, & ! Current iteration
       ITER_K_STATE                 =   1, & ! Current iteration for state variables
       ITER_AUX                     =   2, & ! useful when coupling modules
       TIME_N                       =   3, & ! Time step (converged)
       TIME_N_MINUS_1               =   4, & !
       TIME_N_MINUS_2               =   5, &
       TIME_N_STATE                 =   2    ! Time step (converged) for state variables
  !
  ! Discretization
  !
  integer(ip), parameter   :: &
       NODAL_SCHEME                 =   0, & !
       CELL_CENTERED_SCHEME         =   1
  !
  ! Boundary conditions
  !
  integer(ip), parameter   :: &
       READ_NODE_CODES              =   1_ip,&
       READ_GEOMETRICAL_CODES       =   4_ip,&
       READ_BOUNDARY_CODES          =   2_ip,&
       READ_EDGE_CODES              =   1_ip,&
       IMPOSE_NODE_CODES            =  10_ip,&
       IMPOSE_BOUNDARY_CODES        =  20_ip,&
       IMPOSE_EDGE_CODES            =  40_ip
  !
  ! Coordinate systems
  !
  integer(ip), parameter   :: &
       LOCAL_BASIS_CARTESIAN        =  0_ip, &
       LOCAL_BASIS_CYLINDRICAL      =  1_ip, &
       LOCAL_BASIS_SPHERICAL        =  2_ip, &
       LOCAL_BASIS_CONSTANT         =  3_ip
  !
  ! Functions
  !
  integer(ip), parameter   :: &
       FUNCTION_OFF                 =   0_ip,&  ! 0 : No special function in boundary
       FUNCTION_SPACE_TIME          =   1_ip,&  ! 1 : Space time function
       FUNCTION_TIME                =   2_ip,&  ! 2 : Time function
       FUNCTION_WINDKESSEL          =   3_ip,&  ! 3 : Windkessel function
       FUNCTION_FIELD               =   4_ip,&  ! 4 : Field function
       FUNCTION_MODULE              =   5_ip,&  ! 5 : Module function
       FUNCTION_DISCRETE            =   6_ip,&  ! 6 : Discrete function
       FUNCTION_PUMP                =   7_ip,&  ! 7 : Pump function
       FUNCTION_TUBES               =   7_ip,&  ! 7 : Tubes thridparty
       FUNCTION_USER                =   0_ip    ! 0 : User function


  !
  ! Input/Output units
  !
  integer(ip)              :: &
       lun_pdata,             & ! Data file unit
       lun_outpu,             & ! Outlout (log) file unit
       lun_conve,             & ! Convergence file unit
!       lun_time,              & ! Convergence at time step file unit
       lun_rstar,             & ! Restart file unit
       lun_rstib,             & ! IB restart file unit
       lun_latex,             & ! Latex file unit: text
       lun_gnupl,             & ! Latex file unit: gnuplot
       lun_syste,             & ! System file unit
       lun_tempo,             & ! Temporary file unit
       lun_commu,             & ! Communication with Alya
       lun_binar,             & ! Geometry binary file
       lun_postp,             & ! Postprocess domain file unit
       lun_posvx,             & ! Postprocess Voxel file unit
       lun_pos00,             & ! Additional output file (VU): Domain
       lun_pos01,             & ! Additional output file (VU): Domain
       lun_pos02,             & ! Additional output file (VU): Domain
       lun_mesh,              & ! Basic mesh output unit
       lun_perf,              & ! Performance
       lun_pos04,             & ! Additional output file (VU): Set
       lun_pos05,             & ! Additional output file (VU): Set
       lun_pos06,             & ! Additional output file (VU): Plane
       lun_pos09,             & ! Additional output file (VU): Filter
       lun_pos10,             & ! Additional output file (VU): Filter
       lun_rstla,             & ! Lagragian particles: restart file
       lun_posla,             & ! Lagragian particles: positions
       lun_cvgla,             & ! Lagragian particles: convergence
       lun_mshib,             & ! IB mesh
       lun_resib,             & ! IB result
       lun_mshi2,             & ! IB mesh (2)
       lun_resi2,             & ! IB result (2)
       lun_detec,             & ! Automatic detection unit
       lun_timeline,          & ! Timeline unit
       lun_memory,            & ! Memory evolution
       lun_state                ! State file
  character(150)           :: &
       fil_pdata,             &
       fil_pos00,             &
       fil_pos01,             &
       fil_pos02,             &
       fil_pos00_save,        &
       fil_pos01_save,        &
       fil_pos02_save,        &
       fil_pos04,             &
       fil_pos05,             &
       fil_pos06,             &
       fil_pos07,             &
       fil_pos08,             &
       fil_pos09,             &
       fil_pos10,             &
       fil_mshib,             & ! IB mesh
       fil_resib,             & ! IB result
       fil_mshi2,             & ! IB mesh (2)
       fil_resi2                ! IB result (2)
  !
  ! Others
  !
  real(rp),   parameter    :: &
       zeror = epsilon(1.0_rp)

  !------------------------------------------------------------------------
  !
  ! Run data: read in rrudat
  !
  !------------------------------------------------------------------------

  integer(ip)              :: &
       current_code,          &      ! Current code treated (my code number!)
       kfl_rstar                ! Problem is a restart
  character(66)            :: &
       title                         ! Problem title

  !------------------------------------------------------------------------
  !
  ! Problem data: read in readat
  !
  !------------------------------------------------------------------------

  integer(ip)              :: &
       micou(mblok),          &      ! Max. # of global iter. for each block
       nblok,                 &      ! Number of blocks
       kfl_timco,             &      ! Time coupling strategy
       kfl_timei,             &      ! If there exist a transient problem
       kfl_timef,             &      ! Time step function
       mitim,                 &      ! Maximum number of steps
       mitsm,                 &      ! Mesh refinement: max. # of smoothing iter.
       mitrf,                 &      ! Mesh refinement: max. # of refinement iter.
       kfl_postp_par,         &      ! Parall postprocess type flag
       kfl_wwork,             &      ! Who works
       kfl_lumped,            &      ! Lumped evolution matrix
       kfl_dtfun,             &      ! Time step defined as a piecewise function
       nfundt                        ! Number of DeltaT in case of piecewise function for Time step

  real(rp)                 :: &
       timei,                 &      ! Initial time
       timef,                 &      ! Final time
       dtime                         ! Time step size dt

  real(rp), pointer       ::  &
       dtfun(:,:)                      ! Vector to allocate the piecewise function for time step

  !
  ! Modules information
  !
  integer(ip)              ::      &
       kfl_modul(0:mmodu),         & ! Existence of module
       kfl_coupl(0:mmodu,0:mmodu), & ! Coupling of modules
       kfl_cowhe(0:mmodu,0:mmodu), & ! Field to couple modules
       lmord(mmodu,mblok),         & ! Order of modules iterations
       kfl_delay(mmodu),           & ! Delay module
       kfl_conve(mmodu),           & ! Convergence required
       kfl_solve(0:mmodu),         & ! When to be solved
       kfl_itask(23,0:mmodu),      & ! Task have been already carried out
       lzone(0:mmodu),             & ! Liste of zones
       ndela(mmodu)                  ! Steps to delay module

  !------------------------------------------------------------------------
  !
  ! Variables read by one module and shared by others
  !
  !------------------------------------------------------------------------

  integer(ip)              :: &
       kfl_coibm,             &      ! Immbou coupling: Read by IMMBOU. Shared by kernel
       kfl_cofor,             &      ! Immbou coupling: Read by IMMBOU. Shared by kernel
       kfl_advec,             &      ! Mesh advection: Read by IMMBOU
       kfl_async,             &      ! (A)synchronous communications: Read by PARALL
       kfl_forca_res,         &      ! Forces in residual formulation for rigid body
       kfl_htran                     ! Flag to specify the way diffusive fluxes are computed in energy eq.

  ! Two terms for chemical heat source term
  type(r3p),pointer        :: &
       div_enthalpy_transport(:),  &
       chemical_heat(:),           &
       radiative_heat(:),          &
       enthalpy_transport(:)

  ! Transport properties at Gauss points for flamelet model
  type(r3p),pointer        :: &
       condu_gp(:),           &     ! thermal conductivity
       sphec_gp(:),           &     ! specific heat
       visco_gp(:),           &     ! viscosity
       densi_gp(:),           &     ! density
       wmean_gp(:),           &     ! mean molar weight at different time steps
       tempe_gp(:),           &     ! mean molar weight at different time steps
       sphec_gp_ht(:),        &     ! high temperature cp coefficients
       sphec_gp_lt(:)               ! low temperature cp coefficients

  ! Quantities interchanged between SOLIDZ-TEMPER-CHEMIC for
  ! Poroelasticity/solute transport coupling, at Gauss points
  type(r1p),pointer        :: &
       donna_gp(:)                  ! Donnan osmosis swelling pressure 

  !------------------------------------------------------------------------
  !
  ! Variables read by one module and shared by others for adjoint
  !
  !------------------------------------------------------------------------

  type(r1p),   pointer                  :: &
       RhsadjTem_chm(:),                   & ! Right hand side of the tempe adjoint equation coming from chemic
       RhsadjTem_nsi(:)                      ! Right hand side of the tempe adjoint equation coming from nastin

  type(r2p),   pointer                  :: &
       ddiv_enthalpy_transport_dtem(:),    & ! Derivative of div_enthalpy_transport respect to temperature
       dchemicalHeat_dtem(:),              & ! Derivative of chemicalHeat respect to temperature
       chemicalHeatdiff(:),                & ! Derivative of chemicalHeat respect to design_vars
       RhsadjNas_tur(:),                   & ! Right hand side of the nastin adjoint equation coming from turbulenece
       RhsadjTur_nsi(:),                   & ! Right hand side of the nastin adjoint equation coming from turbulenece
       RhsadjNas_chm(:),                   & ! Right hand side of the nastin adjoint equation coming from chemic
       RhsadjNas_tem(:),                   & ! Right hand side of the nastin adjoint equation coming from tempe
       RhsadjChe_tem(:)                      ! Right hand side of the chemical adjoint equation coming from tempe

  type(r3p),   pointer                  :: &
       dchemicalHeat(:)                      ! Derivative of chemicalHeat respect to concentrations

  ! Primary variables: module unknowns for adjoint
  real(rp), pointer             :: &
       tempe_forw(:,:),            &      !   Known FORWARD Temperature for the adjoint case
       veloc_forw(:,:,:),          &      !   Known FORWARD Velocity for the adjoint case
       press_forw(:,:),            &      !   Known FORWARD Pressure for the adjoint case
       conce_forw(:,:,:),          &      !   Known FORWARD Concentration for the adjoint case
       untur_forw(:,:,:)                  !   Known FORWARD turbulence values for the adjoint case

  !------------------------------------------------------------------------
  !
  ! Modules and services
  !
  !------------------------------------------------------------------------

  integer(ip)              :: &
       modul,                 &      ! Current module
       nmodu,                 &      ! # modules used
       itinn(0:mmodu)                ! Module inner iteration
  integer(8)               :: &
       mem_modul(2,mmodu)            ! Module memory
  real(rp)                 :: &
       cpu_modul(60,mmodu),   &      ! Module CPU time
       dtc_modul(mmodu),      &      ! Module critical time
       glres(mmodu)                  ! Problem residuals
  character(6)             :: &
       namod(0:mmodu)                ! Module name
  character(3)             :: &
       exmod(0:mmodu)                ! Module extension

#ifdef MYTIMING
  real(rp)    :: time_sld_elmope, time_bcsrax                   ! Counters for timing
#endif

  !------------------------------------------------------------------------
  !
  ! Global variables
  !
  !------------------------------------------------------------------------
  !
  ! File names
  !
  character(150)           :: &
       fil_rstar,             &      ! Restart file
       fil_rstib,             &      ! IB Restart file
       fil_conve,             &      ! Convergence file
       fil_postp,             &      ! Postprocess domain file
       fil_posvx,             &      ! Postprocess voxel
       fil_posbs,             &      ! Postprocess set file
       fil_binar                     ! Geometry binary file
  integer(8)               :: &
       memke(2)                      ! Kernel memory
  character(5)             :: &
       wopos_pos(3)                  ! Postprocess variable name
  integer(ip)              :: &
       kfl_split_plus                ! Whether the symbol + should be used as a splitter in ecoute
  !
  ! If arrays should be computed
  !
  integer(ip)              :: &
       kfl_lelbf,             &      ! If list of element boundary faces is required: LELBF
       kfl_lelp2,             &      ! If extended node-element graph is needed: PELPO_2, LELPO_2
       kfl_lele2,             &      ! If extended node-element graph is needed: PELEL_2, LELEL_2
       kfl_symgr,             &      ! If symmetric graph is needed
       kfl_schur,             &      ! If a Schur solver exists
       kfl_aiipr                     ! If a Aii preconditioner exists
  !
  ! Integer
  !
#ifdef ALYA_FTI
  integer(ip), target      :: &
#else
  integer(ip)              :: &
#endif
       ittim                         ! Current time step

  integer(ip)              :: &
       kfl_check_data_file,   &      ! Check data file
       kfl_goblk,             &      ! Block coupling converged
       kfl_gocou,             &      ! Global problem converged
       kfl_gotim,             &      ! Global problem evolving in time
       kfl_stop,              &      ! Global signal to stop the case after the current interation
       kfl_goopt,             &      ! Global problem optimization
       kfl_naked,             &      ! If there is an argument
       kfl_paral,             &      ! -1: Sequential,0: Master,>=1: Slave
       kfl_ptask,             &      ! Partition type
       kfl_livee,             &      ! Live info each [kfl_livee] steps
       ittim_ini,             &      ! Initial time step
       ittim_rst,             &      ! Last restart time step
       itti2,                 &      ! Current time step (restart not included)
       itcou,                 &      ! Current global iteration
       ittyp,                 &      ! Type of iteration
       iblok,                 &      ! Current block
       iitrf,                 &      ! Current refinement iteration
       nusol,                 &      ! Memory variable
       mxdof,                 &      ! Memory variable
       ioutp(50),             &      ! Integer parameters for output
       nturb,                 &      ! Number of turbulence variables
       ITASK_CURREN                  ! Current task
  !
  ! Current CODE/ZONE/SUBDOMAIN
  !
  integer(ip)              :: &
       current_zone,          &      ! Current zone treated
       current_subd,          &      ! Current subdomain
       current_color                 ! Current color
  !
  ! Save variables
  !
  integer(ip)              :: &
       lun_postp_old                 ! Old postprocess unit
  !
  ! ADR equation
  !
  integer(ip)              :: &
       kfl_timei_adr,         &      ! ADR eqn: Time flag
       kfl_advec_adr,         &      ! ADR eqn: Advection flag
       kfl_diffu_adr,         &      ! ADR eqn: Diffusion flag
       kfl_react_adr,         &      ! ADR eqn: Reaction flag
       kfl_grdif_adr,         &      ! ADR eqn: Diffusion gradient flag
       kfl_tisch_adr,         &      ! ADR eqn: Time scheme flag
       kfl_taust_adr,         &      ! ADR eqn: Stab. strategy
       kfl_sgsti_adr,         &      ! ADR eqn: SGS tracking in time
       kfl_sgsno_adr,         &      ! ADR eqn: SGS tracking in non-linearity
       kfl_tiacc_adr,         &      ! ADR eqn: time scheme accuracy
       kfl_shock_adr,         &      ! ADR eqn: shcok capturing flag
       kfl_ellen_adr,         &      ! ADR eqn: element length strategy
       kfl_stead_adr                 ! ADR eqn: steady state flag
  real(rp)                 :: &
       bemol_adr,             &      ! ADR eqn: Integration by parts convective term
       pabdf_adr(10),         &      ! ADR eqn: BDF coefficients
       staco_adr(4),          &      ! ADR eqn: stability constants
       shock_adr,             &      ! ADR eqn: shock capturing coefficients
       safet_adr                     ! ADR eqn: safety factor
  !
  ! Postprocess
  !
  integer(ip), pointer     :: &
       gefil(:)                      ! Filter list
  real(rp)                 :: &
       pafil(10)                     ! Filter parameters
  integer(ip)              :: &
       ivapo,                 &      ! Postprocess variable
       kfl_filte,             &      ! Filter number in module
       kfilt,                 &      ! Filter number
       isect,                 &      ! Live output sections
       inews,                 &      ! Live output New section
       init_ti,               &      ! Inital time to be used in time maasurament with system_clock
       fini_ti,               &      ! Final time to be used in time measurament with system_clock
       cmax_ti,               &      ! Input for system_cock() function.
       crate_ti,              &      ! Represents the precission for system_clock() function. Integer = milisecnds // Real = microseconds
       kfl_origi,             &
       kfl_permu
  integer(8)               :: &
       vtk_id                        ! vkt_id (double precision)
  logical(lg)              :: &
       file_opened                   ! If an opened file was openend
  !
  ! Real
  !
#ifdef ALYA_FTI
  real(rp),    target      :: &
#else
  real(rp)                 :: &
#endif
       dtinv,                 &      ! 1/dt
       dtold(10),             &      ! dt of previous time steps
       dtinv_old(10),         &      ! 1/dt of all time steps
       cutim                        ! Current time

  real(rp)                 :: &
       cpu_initi,             &      ! Initial CPU time
       cpu_times,             &      ! Initial time integration
       cpu_start(CPU_START_END),   & ! CPU for starting operations
       cpu_domain(CPU_END_DOMAIN), & ! CPU for domain operations
       cpu_outpu,             &      ! CPU for output operations
       oltim,                 &      ! Previous time
       routp(80),             &      ! Real parameters for output
       cpu_other(30),         &      ! CPU time
       ini_tim,               &      ! Initial time to be used in time measurment with cputim
       fin_tim,               &      ! Final time to be used in time measurment with cputim
       rate_time                     ! Time rate computed with system_clock
  !
  ! Character
  !
  character(66)            :: &
       namda                         ! Data file name (naked Z)
  character(400)           :: &
       coutp(20)                     ! Character parameter for output
  !
  ! Primary variables: module unknowns
  !
  real(rp), contiguous, pointer :: &
       veloc(:,:,:),               & ! 1:  Velocity                NASTIN-NASTAL-POROUS-EXMEDI
       press(:,:)                    ! 2:  Pressure                NASTIN-NASTAL-SOLIDZ-POROUS
  
  real(rp),             pointer :: &
       tempe(:,:),            &      ! 3:  Temperature             TEMPER-NASTAL
       densi(:,:),            &      ! 4:  Density                 NASTAL
       energ(:,:),            &      ! 5:  Energy                  NASTAL
       visco(:,:),            &      ! 6:  Viscosity               NASTAL
       umome(:,:,:),          &      ! 7:  Momentum                NASTAL
       untur(:,:,:),          &      ! 9:  Turbulence unknowns     TURBUL
       uncdr(:,:,:),          &      ! 10: CDR unknown             CODIRE
       elmag(:,:),            &      ! 11: Electromagnetic         EXMEDI-MAGNET
       dispm(:,:,:),          &      ! 12: Mesh displacement       ALEFOR, IMMBOU, SOLIDZ
       velom(:,:),            &      ! 13: Mesh velocity           ALEFOR
       displ(:,:,:),          &      ! 14: Displacement            SOLIDZ
       spins(:,:),            &      ! 15: Quantic spins
       disfu(:,:,:),          &      ! 16: Distribution funct.     LATBOL
       vdrop(:,:,:),          &      ! 17: Droplet velocity        GOTITA
       cdrop(:,:),            &      ! 19: Droplet concentration   GOTITA
       wavam(:,:),            &      ! 20: Wave amplitude          WAVEQU
       fleve(:,:),            &      ! 21: Level set               LEVELS
       erres(:),              &      ! 22: Error estimator         ADAPTI-...
       vorti(:,:),            &      ! 23: Vorticity and relatives NASTAL-NASTIN
       conce(:,:,:),          &      ! 20: Concentration           PARTIS-NASTAL-CHEMIC
       cdens(:,:),            &      ! 24: Classical density       NASTAL
       entha(:,:),            &      ! 25: Enthalpy                NASTAL
       therm(:,:),            &      ! 26: Thermal variable        TEMPER
       phion(:,:),            &      ! 27: Wave function           QUANTY
       fiber(:,:),            &      ! 28: Fiber fields for composites  EXMEDI-SOLIDZ
       fisoc(:,:),            &      ! 29: Isochrones nodal field for upstroke.
       gpfib(:,:,:),          &      ! 30: Fiber fields in the deformed configuration EXMEDI-SOLIDZ
       radso(:,:)  ,          &      ! 31: Radiation Heat source term
       vconc(:,:,:)  ,        &      ! 32: Concentrations EXMEDI-SOLIDZ
       taulo(:),              &      ! 33: Local electrical wave passage time  EXMEDI - SOLIDZ
       wasat(:,:),            &      ! 34: Water Saturation        POROUS
       neutr(:,:,:,:),        &      ! 35: Neutron flux            NEUTRO
       gdeinv(:,:,:),         &      ! 36: Inverse of the deformation gradient defined per node     SOLIDZ
       gdedet(:),             &      ! 37: Determinant of the deformation gradient defined per node SOLIDZ
       vel1d(:,:),            &      ! 35: 1D velocity              GUSANO
       areas(:,:),            &      ! 36: 1D velocity              GUSANO
       flowr(:,:),            &      ! 37: Flow rate                GUSANO
       elmag_minmax(:,:)             ! 38: Electromagnetic         EXMEDI-MAGNET
       
   
  integer(ip), pointer     :: &
       kfl_fixno_ale(:,:)            ! ALEFOR module fixity condition
  real(rp),    pointer     :: &
       bvess_ale(:,:,:)              ! ALEFOR fixity value
  !
  ! Generic arrays and scalars
  !
  integer(ip)              :: &
       igene,                 &      ! 1st generic integer scalar
       igen2,                 &      ! 2nd generic integer scalar
       iar3p,                 &      ! If gevec allocated
       iasca,                 &      ! If gesca allocated
       iavec,                 &      ! If gevec allocated
       iaten                         ! If geten allocated

  integer(ip), pointer     :: &
       gisca(:),              &      ! Generic integer array
       givec(:,:),            &      ! Generic integer array
       giscp(:),              &      ! Generic integer array in postprocess
       givep(:,:)                    ! Generic integer array in postprocess
  real(rp)                 :: &
       rgene                         ! Generic real scalar
  real(rp),    pointer     :: &
       geten(:,:,:),          &      ! Generic tensor array
       gevec(:,:),            &      ! Generic vector array
       gesca(:),              &      ! Generic scalar array
       getep(:,:,:),          &      ! Generic tensor array in postprocess
       gevep(:,:),            &      ! Generic vector array in postprocess
       gescp(:)                      ! Generic scalar array in postprocess
  complex(rp), pointer     :: &
       getex(:,:,:),          &      ! Generic complex tensor array
       gevex(:,:),            &      ! Generic complex vector array
       gescx(:)                      ! Generic complex scalar array
  type(r3p),   pointer     :: &
       ger3p(:)                      ! Generic r3p type arrays
  real(rp)                 :: &
       funin(3),              &      ! Function in
       funou(3)                      ! Function out
  !
  ! Secondary variables and auxiliars
  !
  real(rp),        pointer :: &
       turmu(:),              &      ! Turb. viscosity mu_t       TURBUL
       rhoon(:,:),            &      ! Density function           QUANTY
       forcf(:,:),            &      ! Fluid force on solid       NASTIN-SOLIDZ
       forca(:,:,:),          &      ! Fluid force on rigid body  NASTIN-ALEFOR (Residual based)
       wmean(:,:),            &      ! Mean molecular weight
       visck(:,:),            &      ! Species viscosity
       condk(:,:),            &      ! Species heat conductivity
       sphek(:,:),            &      ! Species specific heat
       advec(:,:,:),          &      ! Advection for convection terms
       sphec(:,:,:),          &      ! specific heat coefficients used in combination with the enthalpy equation
       massk(:,:),            &      ! Mass source term for each species
       lescl(:),              &      ! Flamelet model: LES closure term for RPV variance (coming from table)
       gdepo(:,:,:),          &      ! Deformation gradient defined per node (for nodal operations)
       momentum_sink(:,:),    &      ! Momentum sink              PARTIS
       heat_sink(:),          &      ! Heat sink                  PARTIS
       mass_sink(:)                ! Mass sink                  PARTIS

#ifdef ALYA_FTI
  real(rp),    target      :: &
#else
  real(rp)                 :: &
#endif
       dpthe,                 &      ! Low-Mach: dp/dt of therm. pressure
       prthe(4)                      ! Thermodynamics pressure (cst or initial)

  real(rp)                 :: &
       tflux ,                &      ! Low-Mach: Total heat flux=int_S k*grad(T).n ds
       epres ,                &      ! Endo pressure, trigger of initial stimuli EXMEDI - SOLIDZ
       tcardiac_cycle ,       &     ! Cardiac cycle length, read by EXMEDI and used by SOLIDZ
       imass                         ! Total initial mass  for low mach


  integer(ip)              :: &
       nspec,                 &      ! Total number of species around
       kfl_prthe                     ! type of thpressure calculation

  !
  ! Chemical Species with their properties
  !
  type (typ_speci),pointer          :: &
       speci(:)

  !
  ! Subgrid scales of primary variables
  !
  type(r3p),   pointer     :: &
       vesgs(:)                      ! Velocity subgrid scale       NASTIN
  type(r3p),   pointer     :: &
       tesgs(:)                      ! Temperature subgrid scale    TEMPER
  type(r3p),   pointer     :: &
       statt(:)                      ! State variables   SOLIDZ

  !
  ! Sets values
  !
  real(rp),    pointer     :: &
       veset(:,:),            &      ! Set element values
       vbset(:,:),            &      ! Set boundary values
       vnset(:,:),            &      ! Set node values
       witne(:,:),            &      ! Witness point values
       witng(:,:)                    ! Witness geometry values
  !
  ! Null targets
  !
  real(rp),    target      :: &
       nul1r(1),              &      ! Null vector
       nul2r(1,1)                    ! Null vector
  integer(ip), target      :: &
       nul1i(1),              &      ! Null vector
       nul2i(1,1)                    ! Null vector
  type(r3p),   target      :: &
       nur3p(1)                      ! Null type r3p

  !------------------------------------------------------------------------
  !
  ! Postprocess and modules
  !
  !------------------------------------------------------------------------

  integer(ip)               :: kfl_reawr       ! postprocess(0)/read(1)/write(2) mode
  type(tymod),  target      :: momod(0:mmodu)
  type(typos),  pointer     :: postp(:)
  type(soltyp), pointer     :: solve(:)
  type(soltyp), pointer     :: solad(:)
  type(eigtyp), pointer     :: eigeg(:)
  type(restyp), pointer     :: reset 
  type(perf),   pointer     :: times(:)
  type(r3p),    pointer     :: gpgradefo ! Gradient deformation tensor per element per gauss point (computed in SOLIDZ)
       

  !------------------------------------------------------------------------
  !
  ! Algebraic system
  !
  !------------------------------------------------------------------------

  real(rp),    pointer,       &
       &       contiguous  :: &
       unkno(:),              &      ! Working array for inner iter.
       amatr(:),              &      ! System matrix
       rhsid(:)                      ! Right-hand-side

  !$acc declare create(rhsid)  
  
  real(rp),    pointer     :: &
       eigen(:),              &      ! Eigen vectors
       eigva(:),              &      ! Eigen values
       bmatr(:),              &      ! System RHS matrix
       pmatr(:),              &      ! Preconditionner
       pschu(:),              &      ! Schur preconditioner
       aii(:),                &      ! Aii matrix
       aib(:),                &      ! Aib matrix
       abi(:),                &      ! Abi matrix
       abb(:),                &      ! Abb matrix
       xxi(:),                &      ! Ui
       xxb(:),                &      ! Ub
       bbi(:),                &      ! bi
       bbb(:),                &      ! bb
       damatr(:),             &      ! Differentiated System matrix
       drhsid(:),             &      ! Differentiated Right-hand-side
       aunkno(:),             &      ! Adjoint Unknown
       lumma(:)                      ! lumped mass matrix for dual time stepping

  complex(rp), pointer     :: &
       amatx(:),              &      ! System matrix
       rhsix(:),              &      ! RHS
       unknx(:),              &      ! Unknown
       pmatx(:),              &      ! Preconditionner
       damatx(:),             &      ! Differentiated System matrix
       drhsix(:),             &      ! Differentiated Right-hand-side
       aunknx(:)                     ! Adjoint Unknown

  !------------------------------------------------------------------------
  !
  ! Parall service
  !
  !------------------------------------------------------------------------

  integer(ip), parameter    :: &
       NELEM_INTE_1DIM = 111,  &
       NELEM_INTE_2DIM = 112,  &
       NELEM_REAL_1DIM = 121,  &
       NELEM_REAL_2DIM = 122,  &
       NELEM_REAL_3DIM = 123,  &
       NBOUN_INTE_1DIM = 211,  &
       NBOUN_INTE_2DIM = 212,  &
       NBOUN_REAL_2DIM = 222,  &
       NBOUN_REAL_3DIM = 223,  &
       NPOIN_INTE_1DIM = 311,  &
       NPOIN_INTE_2DIM = 312,  &
       NPOIN_REAL_1DIM = 321,  &
       NPOIN_REAL_2DIM = 322,  &
       NPOIN_REAL_3DIM = 323,  &
       NPOIN_REAL_12DI = 351,  &
       NBOPO_INTE_1DIM = 411,  &
       NBOPO_REAL_2DIM = 422,  &
       NBOPO_REAL_3DIM = 423,  &
       NELEM_TYPE      =   1,  &
       NBOUN_TYPE      =   2,  &
       NPOIN_TYPE      =   3,  &
       NBOPO_TYPE      =   4,  &
       NFACE_TYPE      =   5,  &
       NEDGE_TYPE      =   6,  &
       ITASK_BROADCAST =   2,  &
       ITASK_SEND      =   3,  &
       ITASK_RECEIVE   =   4,  &
       ITASK_MINIMUM   =   5,  &
       ITASK_SUM       =   9,  &
       ITASK_MAXIMUM   =   10, &
       ITASK_GATHER    =  300

#ifdef ALYA_FTI
  integer(ip), target      :: &
#else
  integer(ip)              :: &
#endif
       npart                         ! Number of partitions

  integer(ip)              :: &
       npari,                 &      ! parin buffer size
       nparl,                 &      ! parlo buffer size
       npasi,                 &      ! parin buffer size (paris)
       npasr,                 &      ! parin buffer size (parrs)
       nparr,                 &      ! parre buffer size
       nparc,                 &      ! parch buffer size
       nparh,                 &      ! parhh buffer size
       nparx,                 &      ! parcx buffer size
       party,                 &      ! array type (e=1,b=2,n=3)
       pardi,                 &      ! array dimension (1,2,3)
       parki,                 &      ! array kind (i=1,r=2,c=3)
       pard1,                 &      ! Dimension 1
       pard2,                 &      ! Dimension 2
       pard3,                 &      ! Dimension 3
       pard4,                 &      ! Dimension 4
       parii,                 &      ! Counter
       gni,                   &      ! Global number of internal nodes
       gnb,                   &      ! Global number of boundary nodes
       lni,                   &      ! Local number of internal nodes
       lnb,                   &      ! Local number of boundary nodes
       icoml,                 &      ! Communication level
       nedgg,                 &      ! Number of edges (mesh subdivision)
       nfacg,                 &      ! Number of faces (mesh subdivision)
       nedgb,                 &      ! Number of boundary edges (mesh subdivision)
       nfacb,                 &      ! Number of face edges     (mesh subdivision)
       kfl_desti_par                 ! Sender/receiver destination

  integer(ip), pointer     :: &
       parin(:),              &      ! Working integer array
       pari1(:),              &      ! Working integer array
       pari2(:,:),            &      ! Working integer array
       pari3(:,:,:),          &      ! Working integer array
       pari4(:),              &      ! Working integer array
       paris(:),              &      ! Working integer array (to send)
       parig(:),              &      ! Working integer array for gather and gatherv
       lnbin(:),              &      ! Working integer array
       lgpar(:),              &      ! Local node to global group
       lnwit(:),              &      ! Witness nodes
       ledgg(:,:),            &      ! List of edges          (mesh subdivision)
       lfacg(:,:),            &      ! List of edges          (mesh subdivision)
       ledgc(:,:),            &      ! List of boundary edges (mesh subdivision)
       lfacb(:,:),            &      ! List of boundary faces (mesh subdivision)
       nelem_tot(:),          &      ! Number of elements
       npoin_tot(:),          &      ! Number of interior+own nodes
       nboun_tot(:),          &      ! Number of boundaries
       npoia_tot(:),          &      ! Number of all nodes
       npoin_par(:),          &      ! npoin
       nelem_par(:),          &      ! nelem
       nboun_par(:),          &      ! nboun
       lninv_loc(:),          &      ! global node     numbering
       leinv_loc(:),          &      ! global element  numbering
       lbinv_loc(:),          &      ! global boundary numbering
       lginv_loc(:),          &      ! global edge     numbering
       lpoi4(:)                      ! List of node to be added to scalar product
  integer(4)               :: &
       sendcount,             &      ! Arrays to be used by AllGatherv
       recvcount                     ! Arrays to be used by AllGatherv
  integer(ip), pointer     :: &
       recvbuf_ip(:),         &      ! Arrays to be used by AllGatherv
       sendbuf_ip(:)                 ! Arrays to be used by AllGatherv
  real(rp),    pointer     :: &
       recvbuf_rp(:),         &      ! Arrays to be used by AllGatherv
       sendbuf_rp(:)                 ! Arrays to be used by AllGatherv
  integer(4),  pointer     :: &
       displs(:),             &      ! Arrays to be used by AllGatherv
       recvcounts(:)                 ! Arrays to be used by AllGatherv
  logical(lg), pointer     :: &
       parlo(:)                      ! Working logical array
  character(len=:), pointer :: &
       parhh                         ! Working character arrays
  !
  ! Renumbering graphs using METIS
  !
  integer(ip)              :: &
       nnren_par,             &      ! Size of graph
       kfl_sfc_par_part,      &      ! Partition the mesh using sfc (= 0 sequential, =1 parallel)
       kfl_sfc_part                  ! Partition the mesh using sfc
  integer(ip), pointer     :: &
       iaren_par(:),          &      ! Renumbering graph IA
       jaren_par(:),          &      ! Renumbering graph JA
       permr_par(:),          &      ! Renumbering graph permutation array
       invpr_par(:)                  ! Renumbering graph inv. permutation array
  !
  ! Others
  !
  type(i1pp),  pointer     :: &
       lelbf(:)                      ! List of element boundary faces
  type(i1p),   pointer     :: &
       lelfa(:)                      ! Element faces list
  real(rp),    pointer     :: &
       parre(:),              &      ! Working real array
       parr1(:),              &      ! Working real array
       parrs(:),              &      ! Working real array (to send)
       parr2(:,:),            &      ! Working real array
       parr3(:,:,:)                  ! Working real array
  complex(rp), pointer     :: &
       parcx(:),              &      ! Working complex array
       parx1(:),              &      ! Working complex array
       parx2(:,:),            &      ! Working complex array
       parx3(:,:,:)                  ! Working complex array

  type(r3p),   pointer     :: &
       par3p(:)                      ! Working r3p array
  type(i1p),   pointer     :: &
       pai1p(:)                      ! Working i1p array

  character(100000)          :: &
       parch                         ! Working character array

  character(40)            :: &
       strre,                 &      ! Name of real variable
       strin,                 &      ! Name of integer variable
       strch,                 &      ! Name of character variable
       strcx                         ! Name of complex variable

  logical(lg)              :: &
       ISLAVE,                &      ! I am a slave        (kfl_paral> 0)
       IMASTER,               &      ! I am the master     (kfl_paral= 0)
       ISEQUEN,               &      ! I am sequential     (kfl_paral=-1)
       IPARALL,               &      ! I am a parallel, master or slave     (kfl_paral>=0)
       INOTMASTER,            &      ! I am not the master (kfl_paral/=0)
       INOTSLAVE,             &      ! I am not a slave    (kfl_paral<=0)
       IPARSLAVE,             &      ! I am a partition pre process slave ()
       IIOSLAVE,              &      ! I am a IO slave ()
       IEMPTY,                &      ! I have no element
       INOTEMPTY                     ! I have elements
     
  
  !------------------------------------------------------------------------
  !
  ! What is beging done flags
  !
  !------------------------------------------------------------------------

  logical(lg)              :: read_restart
  logical(lg)              :: write_restart
  logical(lg)              :: do_repartitioning
  logical(lg)              :: do_AMR

  !------------------------------------------------------------------------
  !
  ! Dodeme service
  !
  !------------------------------------------------------------------------

  integer(ip)              :: &
       ivari_dod,             &      ! Current interpolation variable
       nzdom_dod,             &      ! Graph numbering
       ninte_dod,             &      ! Number of integral interfaces
       ipoin_dod                     ! IPOIN
  real(rp)                 :: &
       bemol_dod                     ! Robin constant bemol (for Robin bc)
  integer(ip), pointer     :: &
       kfl_fixno_dod(:),      &      ! Fixity
       lnsub_dod(:,:),        &      ! List of node subdomains
       lbsub_dod(:,:)                ! List of boundary subdomains

  !------------------------------------------------------------------------
  !
  ! MPIO post process
  !
  !------------------------------------------------------------------------

  real(rp),      pointer   :: &
        gescar_mpio(:),       &
        gevecr_mpio(:,:)
  integer(ip),   pointer   :: &
        gescai_mpio(:),       &
        geveci_mpio(:,:)
  integer(ip)              :: &
        itste_mpio, pdime_mpio, pleng_mpio, tag1_mpio, tag2_mpio
  real(rp)                 :: &
        ttime_mpio
  character(5)             :: &
        wopos_mpio(3)

  type(ibtyp), pointer    :: imbou(:)      ! IB properties
  type(rbtyp), pointer    :: rbbou(:)      ! Rigid body properties
  type(ibint), pointer    :: lnint(:)      ! Interpolation
  type(ibint), pointer    :: lnin2(:)      ! Interpolation

  integer(ip), pointer    :: lntib(:)      ! Node types
  integer(ip), pointer    :: lnti2(:)      ! Node types (previous time stpe)
  integer(ip), pointer    :: lndom(:)      ! Save graph
  integer(ip), pointer    :: lntra(:)      ! List of travesties nodes
  integer(ip), pointer    :: letib(:)      ! Element types
  real(rp),    pointer    :: massc(:,:)    ! Mass conservation matrices (1st column: diagonal mass matrix, 2nd - 4th column: restriccions )

  !----------------------------------------------------------------------
  !
  ! COUPLING AND PARALLELIZATION STUFFS
  !
  !----------------------------------------------------------------------

  logical(lg),    pointer :: &
       I_AM_IN_CODE(:),      & ! Am I in code
       I_AM_IN_ZONE(:),      & ! Am I in zone of my code
       I_AM_IN_SUBD(:)         ! Am I in subd of my code
  character(128), pointer :: &
       application_names(:)    ! Name of applications
  
  !----------------------------------------------------------------------
  !
  ! OPTIMIZATION
  !
  !----------------------------------------------------------------------


  integer(ip)            :: vector_size_variable = 1_ip    ! Size for vectorization

  !----------------------------------------------------------------------
  !
  ! FUNCTIONS
  !
  !----------------------------------------------------------------------

  interface intost
     module procedure intost_4,intost_8
  end interface intost
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-27
  !> @brief   Integer to string
  !> @details Convert an integer(4,8) to a string. This could use class(*)
  !           but PGI crashes when doing select type ( integ )
  !> 
  !-----------------------------------------------------------------------

  function intost_4(integ)
    
    integer(4),   intent(in)   :: integ
    character(20)              :: intost_4
    character(20)              :: intaux
    
    write(intaux,*) integ
    intost_4 = adjustl(intaux)

  end function intost_4

  function intost_8(integ)
    
    integer(8),   intent(in)   :: integ
    integer(4)                 :: integ4
    character(20)              :: intost_8
    character(20)              :: intaux

    integ4 = int(integ,4)    
    write(intaux,*) integ4
    intost_8 = adjustl(intaux)

  end function intost_8

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-27
  !> @brief   Real to string
  !> @details Convert a real to a string
  !> 
  !-----------------------------------------------------------------------

  function retost(realn,REAL_FORMAT)
    
    real(rp)                               :: realn
    character(len=*), intent(in), optional :: REAL_FORMAT
    character(20)                          :: retost
    character(20)                          :: reaux
    integer(ip)                            :: ierr
    
    if( present(REAL_FORMAT) ) then
       write(reaux,REAL_FORMAT,IOSTAT=ierr) realn
       if( ierr /= 0 ) reaux = '0.0'
    else
       write(reaux,'(e19.12)') realn
    end if
    retost=adjustl(reaux)

  end function retost

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-27
  !> @brief   Check if a is an NaN 
  !> @details Check if a is an NaN 
  !> 
  !-----------------------------------------------------------------------

  function isanan(a)
    
    real(rp),    intent(in) :: a
    logical(lg)             :: isanan
    if (a/=a) then
       isanan = .true.
    else
       isanan = .false.
    end if
    return
  end function isanan

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-27
  !> @brief   Check if a is a +/-Inf
  !> @details Check if a is a +/-Inf 
  !> 
  !-----------------------------------------------------------------------

  function isainf(a)
 
    real(rp),    intent(in) :: a
    logical(lg)             :: isainf

    if ((a*0.0_rp)/=0.0_rp) then
       isainf = .true.
    else
       isainf = .false.
    end if
    return

  end function isainf

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-27
  !> @brief   Check if a is a Nan or +/-Inf
  !> @details Check if a is a Nan or +/-Inf 
  !> 
  !-----------------------------------------------------------------------

  function isnain(a)
    
    real(rp),    intent(in) :: a
    logical(lg)             :: isnain

    if ((a*0.0_rp)/=0.0_rp.or.a/=a) then
       isnain = .true.
    else
       isnain = .false.
    end if
    return

  end function isnain

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-27
  !> @brief   Frequency for live info 
  !> @details Frequency for live info 
  !> 
  !-----------------------------------------------------------------------

  function writeliveinfo()
    
    logical(lg) :: writeliveinfo           ! Check logical for kfl_livee

    writeliveinfo = .false.
    if (mod(ittim,kfl_livee)==0) writeliveinfo = .true.

  end function writeliveinfo

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-27
  !> @brief   Parallelization modes
  !> @details Parallelization modes
  !> 
  !-----------------------------------------------------------------------

  function PART_AND_WRITE()
    
    logical(lg) :: PART_AND_WRITE
    if( kfl_ptask == 0 ) then
       PART_AND_WRITE = .true.
    else
       PART_AND_WRITE = .false.
    end if
    return

  end function PART_AND_WRITE

  function PART_AND_RUN()
    
    logical(lg) :: PART_AND_RUN
    if( kfl_ptask == 1 ) then
       PART_AND_RUN = .true.
    else
       PART_AND_RUN = .false.
    end if
    return

  end function PART_AND_RUN

  function READ_AND_RUN()
    
    logical(lg) :: READ_AND_RUN
    if( kfl_ptask == 2 ) then
       READ_AND_RUN = .true.
    else
       READ_AND_RUN = .false.
    end if
    return

  end function READ_AND_RUN

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-27
  !> @brief   If a node is mine
  !> @details If a node is on my own boundary
  !> 
  !-----------------------------------------------------------------------

  function THIS_NODE_IS_MINE(ipoin)
    
    integer(ip), intent(in) :: ipoin
    logical(lg)             :: THIS_NODE_IS_MINE
    
    THIS_NODE_IS_MINE = .false.
    
    if( ISEQUEN ) then
       if( ipoin <= npoin ) THIS_NODE_IS_MINE = .true.
    else if( INOTMASTER ) then
       if( ipoin <= npoi1 .or. ( ipoin >= npoi2 .and. ipoin <= npoi3 ) ) then
          THIS_NODE_IS_MINE = .true.
       else
          THIS_NODE_IS_MINE = .false.
       end if
    end if

  end function THIS_NODE_IS_MINE

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-27
  !> @brief   Give the coupling function between modules 
  !> @details Give the coupling function between modules 
  !> 
  !-----------------------------------------------------------------------

  function coupling(wnam1,wnam2)

    integer(ip)              :: coupling
    character(*), intent(in) :: wnam1
    character(*), intent(in) :: wnam2
    integer(ip)              :: imodu,jmodu

    imodu = idmod(wnam1)
    jmodu = idmod(wnam2)
    if( jmodu == -1 .or. jmodu == -1 ) call runend('WRONG MODULES FOR COUPLING')
    coupling = kfl_coupl(imodu,jmodu)

  end function coupling

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-27
  !> @brief   Convert a module name into its ID
  !> @details Convert a module name into its ID
  !> 
  !-----------------------------------------------------------------------

  function idmod(wname)
    
    integer(ip)              :: idmod
    character(*), intent(in) :: wname
    integer(ip)              :: kmodu

    kmodu = -1
    idmod = -1
    do while( kmodu < mmodu )
       kmodu = kmodu + 1
       if( wname(1:5) == namod(kmodu)(1:5) ) then
          idmod = kmodu
          kmodu = mmodu
       end if
    end do

  end function idmod

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-21
  !> @brief   Name to id
  !> @details name to id
  !> 
  !-----------------------------------------------------------------------

  integer(ip) function module_name_to_id(wname)

    character(len=*) :: wname
    integer(ip)      :: imodu
        
    do imodu = 1,mmodu
       if( wname(1:5) == namod(imodu)(1:5) .and. trim(namod(imodu)) /= '' ) then
          module_name_to_id = imodu
          return
       end if
    end do
    
    call runend('MODULES: UNKNOWN MODULE')
    
  end function module_name_to_id
  
end module def_master
