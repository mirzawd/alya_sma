!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_reaphy()
  !------------------------------------------------------------------------
  !****f* Chemic/chm_reaphy
  ! NAME
  !    chm_reaphy
  ! DESCRIPTION
  !    This routine reads the physical problem definition
  !    NREAC_CHM ........................ Number of reactions
  !    LREAC(NREAC_CHM)%L(:) ............ List of reactions
  !    REACT_CHM(NCOEF_CHM,NREAC_CHM) ... Reaction coefficients
  !
  ! USES
  ! USED BY
  !    chm_turnon
  !***
  !------------------------------------------------------------------------
  use def_inpout,             only : words, exists, getrea, getint, param, getcha
  use def_master,             only : ID_CHEMIC, ID_NASTIN, INOTSLAVE, momod, mem_modul, modul, kfl_htran, kfl_rstar, prthe, intost,&
                                     speci, kfl_coupl
  use def_kintyp,             only : ip, rp, lg
  use def_chemic,             only : nclas_chm, table_fw, posttable_fw, mixedEq_groups_chm, mixedEq_eqs_chm, bf_fuel_chm,&
                                     bo_oxy_chm, dac_cor_chm, dac_crit_chm, droplet_compactness_limit_chm, droplet_h_factor_chm,&
                                     droplet_max_diameter_chm, Field_ind_chm, index_N2, kfl_advec_chm, kfl_bc_alpha_CMC_chm,&
                                     kfl_bc_init_method_CMC_chm, kfl_cont_chm, kfl_control_gr_chm, kfl_diffu_chm,&
                                     kfl_droplet_id_chm, kfl_field_chm, kfl_freq_chm, kfl_hrr_col_chm, kfl_hrr_fw_chm,&
                                     kfl_incl_PDF_trans_CMC_chm, kfl_key_chm, kfl_lookg_chm, kfl_model_chm, kfl_norma_chm,&
                                     kfl_pfa_chm, kfl_post_fw_chm, kfl_premix_chm, kfl_radia_chm, kfl_solve_cond_CMC_chm,&
                                     kfl_solve_enth_CMC_chm, kfl_soot_chm, kfl_spec_name_chm, kfl_split_CFD_CMC, kfl_spray_chm,&
                                     kfl_tab_fw_chm, kfl_tdac_write_chm, kfl_timei_chm, kfl_trans_mxt_spc_CMC_chm,&
                                     kfl_trans_phs_spc_CMC_chm, kfl_transport_chm, kfl_ufpv_chm, kfl_varYc_chm, kfl_varZ_chm,&
                                     kfl_weigh_in_eq_CMC_chm, kfl_yk_fw_ssm, kfl_z_chm, levelSet_threshold_chm, mechanism_path,&
                                     mf_vec_path_CMC_chm, nclas_chm, ngrou_chm, nreac_chm, nS_AMC_CMC_chm,&
                                     nsize_mech_name, nsize_mf_vec_chm, nsize_red, nspec_chm, nvar_CMC_chm, nZ_AMC_CMC_chm,&
                                     nZ_CMC_chm, posttable_fw, prthe_chm, radwt_chm, Red_spec, S_threshold, Smax_AMC_CMC_chm,&
                                     socen_chm, sorad_chm, surf_tension_chm, table_coords, table_fw, table_tab, Xtot_const_CMC_chm,&
                                     Zavg_const_CMC_chm, Zs_CMC_chm, Zstq_CMC_chm, Zvar_const_CMC_chm, Le_k, Z_CMC_chm,&
                                     T_bc_CMC_chm, react_scal_bc_CMC_chm, diffu_chm, diff_Z_CMC_chm, kfl_multimod_chm,&
                                     chm_zmax, chm_zmin, kfl_tab_fw_chm_diff, kfl_tab_fw_chm_prem, kfl_zg_col_chm, kfl_zg_fw_chm

  use def_kermod,             only : gasco, lookup_fw
  use mod_ecoute,             only : ecoute
  use mod_interp_tab,         only : tab_load_file
#ifdef CANTERA
  use def_chemic,             only : gas_chm
  use cantera,                only : importPhase, nSpecies, nReactions, getSpeciesName, getReactionString
#endif
  use iso_fortran_env,        only : iostat_eor,iostat_end
  use mod_chm_operations_CMC, only : chm_initial_actions_reaphy_CMC
  use mod_chm_sectional_soot_model_fast, only: Vmax_ssm
  use mod_chm_sectional_soot_model_fast, only: RadSoot_ssm
  use mod_chm_sectional_soot_model_fast, only: RhoC_ssm
  use mod_chm_sectional_soot_model_fast, only: nsect_ssm
  use mod_chm_sectional_soot_model_fast, only: nspec_ssm
  use mod_chm_sectional_soot_model_fast, only: nclas_ssm
  use mod_chm_sectional_soot_model_fast, only: gasCoupling_ssm
  use mod_chm_sectional_soot_model_fast, only: ID_NUCL,ID_COND,ID_COAG,ID_SURF

  use mod_chm_mixedEq,       only : chm_mixedEq_setGroupType
  use mod_chm_mixedEq,       only : chm_mixedEq_setEqType
  use mod_chm_mixedEq,       only : chm_mixedEq_setSrcType
  use mod_chm_mixedEq,       only : chm_mixedEq_setIniType
  use mod_chm_mixedEq,       only : CHM_GR_CONTROL
  use mod_chm_mixedEq,       only : CHM_INI_OLDWAY

  use mod_messages,          only : messages_live
  use mod_memory,            only : memory_alloca


  implicit none
  integer(ip)                   :: iclas,kfl_forwa,nsize_spec_name,ipostvar,jpostvar, &
                                   index_zbc,icoun,line_f,state_f
  integer(ip)                   :: iequa, igrou
  character(5)                  :: typ,eqname
  character(5)                  :: stoichio_names(3) = ''
  integer(ip)                   :: nbuff, ioerror, ii, jj
  character(len=:), allocatable :: buffer, spec_name
  integer(ip)                   :: kfl_stoic
  logical                       :: found
  logical(lg)                   :: give_new_name
  real(rp)                      :: Sc_turb_loc

#ifdef CANTERA
  integer(ip)                   :: ireac
  character(80)                 :: eq
#endif

  external                      :: runend
  external                      :: chm_memphy
  external                      :: cantera_initialization
  external                      :: cantera_trim

  if( INOTSLAVE ) then
     !
     ! Initializations (defaults)
     !
     kfl_model_chm    = 0                        ! Defect evolution
                                                 !  (1) Flamelet model / (2) Optimized framework for flamelets /
                                                 !  (3) Finite rate kinetics / (4) CMC
     kfl_timei_chm    = 1                        ! Existence of du/dt
     kfl_advec_chm    = 0                        ! Existence of (a.grad)u
     kfl_diffu_chm    = 1                        ! Existence of -div[k. grad(u)], and its options (include density)
     kfl_field_chm    = 0                        ! Initialization by fields (=0 OFF, /=0 ON)
     kfl_key_chm      = 0                        ! Number of Key species for Dynamic Chemistry
     kfl_multimod_chm = 0                        ! Multi-regime model
     kfl_z_chm        = 0                        ! Flag for z
     kfl_tdac_write_chm    = 0                   ! Flag to write reduced table
     kfl_freq_chm     = 0                        ! Number of interations between DAC reduction
     dac_cor_chm      = 0.0_rp                   ! CODAC Corrleation
     bf_fuel_chm      = 0.0_rp                   ! B_F for fuel
     bo_oxy_chm       = 0.0_rp                   ! B_O for Oxygen
     kfl_radia_chm    = 0                        ! radiation model activated?
     kfl_pfa_chm      = -1                       ! Reduction method based on PFA: =-1 NO REDUCTION, =0 EQUILIBRIUM,
                                                 !   =1 PFA (static), =2 DAC, =3, CODAC
     dac_crit_chm     = 0.0_rp                   ! Critical PFA Value
     kfl_premix_chm   = 0                        ! Combustion model: =0 (NON-PREMIXED) = 1 (PREMIXED)
     kfl_spray_chm    = 0                        ! Activate Spray model
     kfl_ufpv_chm     = 0                        ! Unsteady Flamelet Progress Variable model 0: original, 1: steady FPV,
                                                 !  2: omegaYc, 3: Yc_dot
     kfl_varZ_chm     = 0
     kfl_varYc_chm    = 0
     kfl_control_gr_chm = 0
     kfl_lookg_chm    = 1                        ! Lookup of table properties, 0: on nodes, 1: on gauss points
     kfl_soot_chm     = 0                        ! Activation soot model, 0: OFF, 1: sectional method, -1: validation
     nsect_ssm        = 0                        ! Number of sections soot model
     kfl_yk_fw_ssm    = 0                        ! Index of lookup framework for soot precursor species
     kfl_tab_fw_chm   = -1                       ! Lookup framework
     kfl_post_fw_chm  = -1                       ! Postprocessing lookup framework
     kfl_hrr_fw_chm   = -1                       ! Heat release rate lookup framework
     kfl_hrr_col_chm  = -1                       ! Heat release rate lookup column in framework
     kfl_zg_fw_chm    = -1                       ! Zgrad lookup framework
     kfl_zg_col_chm   = -1                       ! Zgrad lookup column in framework
     kfl_tab_fw_chm_diff  = -1                   ! Properties based on diffusion table
     kfl_tab_fw_chm_prem  = -1                   ! Properties based on premixed table

     kfl_norma_chm = 0                           ! Normalize concentrations
     kfl_forwa     = 0                           ! Forward/backward reaction rate

     prthe_chm     = 0.0_rp                      ! Thermodynamic pressure (if NASTIN is not activated)
     nclas_chm     = 1                           ! Number of classes
     ngrou_chm     = 1                           ! Number of equation groups
     radwt_chm     = 0.0_rp                      ! Wall temperature for radiation model
     sorad_chm     = 0.0_rp                      ! Source radius
     socen_chm     = 0.0_rp                      ! Source center
     nreac_chm     = 0                           ! Number of reactions

     Sc_turb_loc   = 0.0_rp                      ! Turbulent Schmidt number

     kfl_droplet_id_chm            = 0_ip        ! Droplet identification flag
     levelSet_threshold_chm        = 0.0_rp      ! Level Set threshold for droplet identification
     droplet_compactness_limit_chm = 1.0_rp      ! Compactness value below which a cluster won't be considered as a droplet
     droplet_max_diameter_chm      = 0.0_rp      ! Max. diameter above which a cluster won't be considered as a droplet
     droplet_h_factor_chm          = 0.0_rp      ! Mesh size factor to define a max. diameter as:
                                                 !  droplet_h_factor_chm * h_mesh

     kfl_stoic         = 0_ip                    ! Stoichiometric mass ratio defined or not
     stoichio_names    = ''                      !
     mechanism_path    = ''                      ! Path to mechanism: can be global path, relative path
     Red_spec          = ''                      ! Reduced Species
     ID_NUCL           = 1                       ! ID for Nucleation
     ID_COND           = 1                       ! ID for Condensation
     ID_COAG           = 1                       ! ID for Coagulation
     ID_SURF           = 1                       ! ID for Surface growth
     kfl_spec_name_chm = 0_ip                    ! Number of fields read (finite rate)

     !
     ! Variables for CMC model
     !
     kfl_weigh_in_eq_CMC_chm = 0_ip              ! Flag to weigh inert and equilibrium solutions when starting from scratch
     kfl_split_CFD_CMC = 1_ip                    ! Flag to split CFD and CMC into two different executions:
                                                 !  0 CFD and CMC in same execution and 1 if CFD and CMC separated
     kfl_solve_enth_CMC_chm = 0_ip               ! 0 if enthalpy is not solved and 1 if it is transported
     kfl_solve_cond_CMC_chm = -1_ip              ! 0: equations for mixture fraction, its variance, etc. are
                                                 !  solved (unconditional variables);
                                                 !  1: equations for conditional species (default)
     kfl_bc_init_method_CMC_chm = 1_ip           ! 1: linear; 2: homogeneous reactors
     kfl_incl_PDF_trans_CMC_chm = 1_ip           ! 0: do not include PDF in CMC transport equations, 1: include it
     kfl_trans_phs_spc_CMC_chm = 1_ip            ! 0: do not include transport in physical space in CMC equations, 1: include it
     kfl_trans_mxt_spc_CMC_chm = 1_ip            ! 0: do not include transport in mixture fraction space in CMC equations,
                                                 !  1: include it
     kfl_bc_alpha_CMC_chm  = -1_ip               ! 0: boundary conditions defined through species;
                                                 !  1: boundary conditions defined through alpha
     nZ_CMC_chm = 1_ip                           ! Number of slices in mixture fraction space
     nvar_CMC_chm = 0_ip                         ! Number of variables to be solved in CMC:
                                                 !  nvar_CMC_chm = nclas_chm+1 (species + enthalpy)
     nZ_AMC_CMC_chm = 51_ip                      ! Number of mixture fractions for scalar dissipation rate integration (AMC model)
     nS_AMC_CMC_chm = 31_ip                      ! Number of segregation factors for scalar dissipation rate integration (AMC model)
     index_N2 = 0_ip                             ! Index in the mechanism for N2
     Zs_CMC_chm = 1.0_rp                         ! Saturation mixture fraction
     Smax_AMC_CMC_chm = 0.3_rp                   ! Maximum segregation factor for AMC model
     S_threshold = 9.9e-3_rp                     ! Threshold segregation factor for integrations
     Zstq_CMC_chm = -1.0e-4_rp                   ! Stoichiometric mixture fraction
     Zavg_const_CMC_chm = 0.0_rp                 ! Constant value for Zavg
     Zvar_const_CMC_chm = 0.0_rp                 ! Constant value for Zvar
     Xtot_const_CMC_chm = 0.0_rp                 ! Constant value for total scalar dissipation rate

     !
     ! Reach the section
     !
     call ecoute('chm_reaphy')
     do while(words(1)/='PHYSI')
        call ecoute('chm_reaphy')
     end do
     !
     ! Begin to read data
     !
     physical_problem: do while(words(1)/='ENDPH')
        call ecoute('chm_reaphy')

        if( words(1) == 'PROBL' ) then
           !
           ! Problem definition data
           !
           call ecoute('chm_reaphy')

           problem_definition: do while(words(1)/='ENDPR')

              !
              ! Chemistry model:
              !
              if( words(1) == 'MODEL' ) then
                 if( words(2) == 'FINIT' ) then
                    kfl_model_chm   = 3
                    write(momod(modul) % lun_outpu,*) 'Using finite rate model for calculations'
                    write(momod(modul) % lun_outpu,*)

                 else if (words(2) == 'CMC  ') then
                    kfl_model_chm   = 4
                    write(momod(modul) % lun_outpu,*) 'Using CMC model for calculations'
                    write(momod(modul) % lun_outpu,*)
                    if (words(3) == 'CONDI') then
                       kfl_solve_cond_CMC_chm = 1_ip
                       write(momod(modul) % lun_outpu,*) 'Solving species (and enthalpy) conditional CMC transport equations'
                       write(momod(modul) % lun_outpu,*)
                    else if (words(3) == 'UNCON') then
                       kfl_solve_cond_CMC_chm = 0_ip
                       write(momod(modul) % lun_outpu,*) 'Solving mixture fraction, its variance, etc. unconditional transport&
                           & equations'
                       write(momod(modul) % lun_outpu,*)
                       nclas_chm = 4 ! Even if in this case only mixture fraction and its
                                     ! variance transport equations are solved, progress
                                     ! variable and its variance are considered but they
                                     ! are meaningless in this context
                       nspec_chm = nclas_chm
                    end if

                    if (kfl_solve_cond_CMC_chm == -1)   call runend('CHEMIC REAPHY: Type of CMC variables to be solved not&
                        & specified')

                 else if( words(2) == 'FLAME' ) then
                    kfl_model_chm   = 1
                    !
                    ! Type of combustion problem
                    !
                    nclas_chm                   = 4  ! Non-premixed
                    kfl_premix_chm = 0_ip
                    if( .not.exists('NONPR') ) then  ! Premixed
                       nclas_chm                = 2
                       kfl_premix_chm = 1_ip
                    end if

                    !
                    ! Radiation model
                    !
                    if( exists('RADIA') ) then
                       kfl_radia_chm = 1
                       if ( exists('WALLT') ) then
                          radwt_chm = getrea('WALLT',0.0_rp,'#Wall temperature radiation')
                          kfl_radia_chm = 2
                       end if
                    end if

                    !
                    ! Unsteady flamelet progress variable model
                    !
                    if( exists('UFPV ') ) then
                       kfl_ufpv_chm = 3   ! default source is dot{Yc}

                       if( exists('OMEGA') ) kfl_ufpv_chm = 2   ! source is omega_Yc
                       if( exists('STEAD') ) kfl_ufpv_chm = 1   ! steady flamelet model
                    end if

                 else if( words(2) == 'SPRAY' ) then
                    kfl_premix_chm = 1_ip
                    kfl_spray_chm  = 1
                    kfl_model_chm  = 1
                    nclas_chm      = 4

                 else if( words(2) == 'MIXED' ) then
                    !
                    ! Mixed model
                    !
                    kfl_model_chm   = 2
                    nclas_chm       = 0

                 end if

              !
              ! Z variance model for Flamelet models
              !
              else if( words(1) == 'ZVARI' ) then

                 !
                 ! Transport variance
                 !
                 if( words(2) == 'ZV   ' ) then
                    kfl_varZ_chm = 1
                    if( exists('LEA  ') ) then
                       kfl_varZ_chm = -1
                       call runend('CHEMIC REAPHY: LEA model for Z not implemented')
                    end if
                 !
                 ! Transport Z*Z
                 !
                 else if( words(2) == 'ZZ  ') then
                    kfl_varZ_chm = 2
                    if( exists('LEA  ') ) then
                       kfl_varZ_chm = -2
                       call runend('CHEMIC REAPHY: LEA model for Z not implemented')
                    end if
                 !
                 ! Variance is 0
                 !
                 else if( words(2) == 'OFF  ') then
                    kfl_varZ_chm = 0
                 end if

              !
              ! Yc variance model for Flamelet models
              !
              else if( words(1) == 'YCVAR' ) then               ! Yc variance model for flamelet combustion model
                 !
                 ! Transport variance
                 !
                 if( words(2) == 'YCV  ' ) then
                    kfl_varYc_chm = 1
                 !
                 ! Transport Yc*Yc
                 !
                 else if( words(2) == 'YCYC ') then
                    kfl_varYc_chm = 2
                 !
                 ! Variance is 0
                 !
                 else if( words(2) == 'OFF  ') then
                    kfl_varYc_chm = 0
                 end if

              !
              ! LOOKUP LOCATION
              !
              else if( words(1) == 'LOOKU' ) then

                 !
                 ! Transport variance
                 !
                 if( words(2) == 'GAUSS' ) then
                    kfl_lookg_chm = 1
                 elseif ( words(2) == 'NODES' ) then
                    kfl_lookg_chm = 0
                 endif

              else if( words(1) == 'SPRAY' ) then               ! Problem name
                 if( words(2) == 'EULER' ) then
                    kfl_spray_chm = 1
                    if (kfl_premix_chm == 0) nclas_chm = nclas_chm + 2

                 else if( words(2) == 'LEVEL') then
                    kfl_spray_chm = 2
                    !call runend('CHEMIC REAPHY: Level set model for ELSA not available yet, be patient')

                 else if( words(2) == 'ELSA ') then
                    kfl_spray_chm = 3
                    call runend('CHEMIC REAPHY: Full ELSA model not available yet, be patient')
                 end if

              else if( words(1) == 'TEMPO' ) then               ! Temporal evolution
                 if( words(2) == 'ON   ' ) then
                    kfl_timei_chm = 1
                 else
                    kfl_timei_chm = 0
                 end if

              else if( words(1) == 'CONVE' ) then               ! Convective term
                 if( words(2) == 'ON   ' ) then
                    kfl_advec_chm = -2
                    if(words(3)=='VELOC' ) then
                       kfl_advec_chm = -2
                    end if
                 else if( words(2) == 'OFF  ' ) then
                    kfl_advec_chm = 0
                 end if

              else if( words(1) == 'DIFFU' ) then               ! Diffusion
                 if( words(2) == 'ON   ' ) then
                    kfl_diffu_chm = 1
                 else
                    kfl_diffu_chm = 0
                 end if

              else if( words(1) == 'HTRAN' ) then               ! Fluxes for enthalpy transport
                 if( words(2) == 'DETAI' ) then                 ! Detailed flux
                    kfl_htran = 1
                 end if


              !
              ! Fuel and Oxidiser definition (used to calculate mixture Fraction)
              !
              else if( words(1) == 'ZSCAL') then
                 bf_fuel_chm  = getrea('FUEL ',-1.0_rp,'#Fuel - Chemic')
                 bo_oxy_chm   = getrea('OXY  ',-1.0_rp,'#Fuel - Chemic')
                 kfl_z_chm    = -1_ip

                 !
                 ! PFA Reduction - Read Non Key Species
                 !
              else if( words(1) == 'REDUC') then

                 if(words(2)=='PFA  ') then
                    !
                    ! Use equilibrium condition for non-Key species
                    !
                    if (exists('EQUIL')) then
                       kfl_pfa_chm = 0
                       call ecoute('chm_reaphy')
                       backspace(momod(modul) % lun_pdata)
                       nbuff = 1000
                       allocate(character(len=nbuff) :: buffer)
                       read(momod(modul) % lun_pdata,'(a)', advance='no', size=nbuff, iostat=ioerror) buffer
                       if ((ioerror == iostat_eor).or.(ioerror == iostat_end)) ioerror = 0
                       if (ioerror /= 0 ) call runend('CHEMIC REAPHY: Cannot read Reduced Species.')
                       buffer = adjustl(buffer) ! trim leading whitespace
                       buffer = trim(buffer)    ! trim trailing whitespace
                       Red_spec  = buffer
                       nsize_red = len(Red_spec, kind=ip)
                       deallocate(buffer)

                    elseif (exists('STATI')) then
                       kfl_pfa_chm   = 1
                       kfl_key_chm   = getint('KEYSP',0_ip,'#Species starting fields')
                       dac_crit_chm  = getrea('CRITI',0.0_rp,'#Species starting fields')
                       call ecoute('chm_reaphy')
                       backspace(momod(modul) % lun_pdata)
                       nbuff = 1000
                       allocate(character(len=nbuff) :: buffer)
                       read(momod(modul) % lun_pdata,'(a)', advance='no', size=nbuff, iostat=ioerror) buffer
                       if ((ioerror == iostat_eor).or.(ioerror == iostat_end)) ioerror = 0
                       if (ioerror /= 0 ) call runend('CHEMIC REAPHY: Cannot read Reduced Species.')
                       buffer = adjustl(buffer) ! trim leading whitespace
                       buffer = trim(buffer)    ! trim trailing whitespace
                       Red_spec  = buffer
                       nsize_red = len(Red_spec, kind=ip)
                       deallocate(buffer)

                       !
                       ! Dynamic Adaptive Chemistry (DAC)
                       !  ("Red_Spec" holds key species)
                       !
                    else if(exists('DAC  ')) then
                       kfl_pfa_chm = 2
                       kfl_key_chm   = getint('KEYSP',0_ip,'#Species starting fields')
                       dac_crit_chm  = getrea('CRITI',0.0_rp,'#Species starting fields')
                       kfl_freq_chm  = getint('FREQU',0_ip,'#Steps between sucssive DAC reductions')
                       call ecoute('chm_reaphy')
                       backspace(momod(modul) % lun_pdata)
                       nbuff = 1000
                       allocate(character(len=nbuff) :: buffer)
                       read(momod(modul) % lun_pdata,'(a)', advance='no', size=nbuff, iostat=ioerror) buffer
                       if ((ioerror == iostat_eor).or.(ioerror == iostat_end)) ioerror = 0
                       if (ioerror /= 0 ) call runend('CHEMIC REAPHY: Cannot read Reduced Species.')
                       buffer = adjustl(buffer) ! trim leading whitespace
                       buffer = trim(buffer)    ! trim trailing whitespace
                       Red_spec  = buffer
                       nsize_red = len(Red_spec, kind=ip)
                       deallocate(buffer)

                       !
                       ! Correlated Dynamic Adaptive Chemistry (CODAC)
                       !  ("Red_Spec" holds key species)
                       !
                    else if(exists('CODAC')) then
                       kfl_pfa_chm = 3
                       kfl_key_chm  = getint('KEYSP',0_ip,'#Species starting fields')
                       dac_crit_chm = getrea('CRITI',0.0_rp,'#Species starting fields')
                       dac_cor_chm  = getrea('CORRE',0.0_rp,'#Species starting fields')
                       call ecoute('chm_reaphy')
                       backspace(momod(modul) % lun_pdata)
                       nbuff = 1000
                       allocate(character(len=nbuff) :: buffer)
                       read(momod(modul) % lun_pdata,'(a)', advance='no', size=nbuff, iostat=ioerror) buffer
                       if ((ioerror == iostat_eor).or.(ioerror == iostat_end)) ioerror = 0
                       if (ioerror /= 0 ) call runend('CHEMIC REAPHY: Cannot read Reduced Species.')
                       buffer = adjustl(buffer) ! trim leading whitespace
                       buffer = trim(buffer)    ! trim trailing whitespace
                       Red_spec  = buffer
                       nsize_red = len(Red_spec, kind=ip)
                       deallocate(buffer)

                       !
                       ! Tabulated Dynamic Adaptive Chemistry (T-DAC)
                       !  ("Red_Spec" holds defnition of progress variable)
                       !
                    else if(exists('TDAC ')) then
                       kfl_pfa_chm = 4
                       kfl_key_chm  = getint('KEYSP',0_ip,'#Species starting fields')
                       kfl_cont_chm = getint('CONTR',0_ip,'#Controling Varibales TDAC')
                       dac_crit_chm = getrea('CRITI',0.0_rp,'#Species starting fields')
                       if( exists('WRITE') ) kfl_tdac_write_chm = -10_ip
                       if( exists('READ ') ) kfl_tdac_write_chm = 10_ip
                       if( exists('FILL ') ) kfl_tdac_write_chm = 666_ip
                       call ecoute('chm_reaphy')
                       backspace(momod(modul) % lun_pdata)
                       nbuff = 1000
                       allocate(character(len=nbuff) :: buffer)
                       read(momod(modul) % lun_pdata,'(a)', advance='no', size=nbuff, iostat=ioerror) buffer
                       if ((ioerror == iostat_eor).or.(ioerror == iostat_end)) ioerror = 0
                       if (ioerror /= 0 ) call runend('CHEMIC REAPHY: Cannot read Reduced Species.')
                       buffer = adjustl(buffer) ! trim leading whitespace
                       buffer = trim(buffer)    ! trim trailing whitespace
                       Red_spec  = buffer
                       nsize_red = len(Red_spec, kind=ip)
                       deallocate(buffer)
                    end if

                 else if (exists('NONE ') .or. exists('OFF  ')) then
                    !
                    ! Complete Integration of the mechanism
                    !
                    kfl_pfa_chm = -1_ip

                 else
                    call runend('CHEMIC REAPHY: Reduction activated, but no additional information provided')
                 end if
                 !
                 ! Reactions
                 !
              else if( words(1) == 'MECHA' ) then
                 !
                 ! Cantera
                 !
                 if( words(2) == 'CANTE' ) then

                    !
                    ! Read mechanism path from:   momod(modul) % fil_pdata   with ID: momod(modul) % lun_pdata
                    !
                    ii = 0
                    do while(words(1)/='ENDME')
                       ii = ii + 1
                       !
                       ! Call ecoute to look for the END_REACT, and step back
                       !
                       call ecoute('chm_reaphy')
                       backspace(momod(modul) % lun_pdata)

                       !
                       ! Read next line
                       !
                       nbuff = 100
                       allocate(character(len=nbuff) :: buffer)
                       read(momod(modul) % lun_pdata,'(a)', advance='no', size=nbuff, iostat=ioerror) buffer

                       !
                       ! Check read, disregard errors about buffer being larger
                       ! than the line.
                       !
                       if ((ioerror == iostat_eor).or.(ioerror == iostat_end)) ioerror = 0
                       if (ioerror /= 0 ) call runend('CHEMIC REAPHY: Cannot read cantera mechanism path.')

                       !
                       ! Process line
                       !
                       buffer = adjustl(buffer) ! trim leading whitespace
                       buffer = trim(buffer)    ! trim trailing whitespace
                       jj     = 0_ip
                       jj     = index(buffer,'xml', kind=ip)
                       jj     = jj + 2_ip
                       if (mechanism_path=='') then
                          !
                          ! Take first line as path (other lines are ignored
                          ! for now)
                          !
                          mechanism_path  = buffer(1:jj)
                          nsize_mech_name = len(mechanism_path, kind=ip)
                       endif

                       deallocate(buffer)

                    end do
#ifdef CANTERA
                    gas_chm   = importPhase(mechanism_path)
                    nclas_chm = nSpecies(gas_chm)    ! number of species
                    nspec_chm = nclas_chm
                    nreac_chm = nReactions(gas_chm)  ! number of reactions
#else
                    call runend('CHEMIC REAPHY: please compile with Cantera.')
#endif
                 else
                    call runend('CHEMIC REAPHY: Mechanism selection only valid with FINITE RATE model')
                 endif

              end if

              call ecoute('chm_reaphy')
           end do problem_definition

        else if( words(1) == 'PROPE' ) then
           !
           ! Properties
           !
           call chm_memphy(1_ip)

           call ecoute('chm_reaphy')
           !
           write(momod(modul) % lun_outpu,*)'-----------------'
           write(momod(modul) % lun_outpu,*)' FLOW PROPERTIES '
           write(momod(modul) % lun_outpu,*)'-----------------'
           write(momod(modul) % lun_outpu,*) ''
           !
           properties: do while(words(1)/='ENDPR')

              !-------------------------------------------------------
              !
              ! Properties and constants
              !
              !-------------------------------------------------------

              if ( words(1) == 'GASCO' ) then              ! Universal gas constant
                 gasco = getrea('GASCO',8.3144621_rp,'#Gas constant')

              else if ( words(1) == 'SURFA' ) then              ! Surface tension
                 surf_tension_chm = getrea('SURFA',1.0_rp,'#Surface tension liquid')
                 !
                 ! Variable names for ELSA model
                 !
                 if (kfl_spray_chm /= 0) then
                    if (nclas_chm == 4_ip ) then
                       speci(1_ip)%name = 'NONE'
                       speci(2_ip)%name = 'NONE'
                       speci(3_ip)%name = 'PHI'
                       speci(4_ip)%name = 'SIGMA'
                    elseif (nclas_chm == 6_ip ) then
                       speci(5_ip)%name = 'PHI'
                       speci(6_ip)%name = 'SIGMA'
                    end if
                 end if

              else if ( words(1) == 'THERM' ) then              ! Thermodynamic pressure
                 prthe_chm = getrea('THERM',101325.0_rp,'#Thermodynamic pressure')
                 if (kfl_coupl(ID_CHEMIC,ID_NASTIN) == 0 ) prthe = prthe_chm

              else if( words(1) == 'TURBU' ) then
                 Sc_turb_loc = getrea('TURBU',0.9_rp,'#Turbulent Schmidt number')    ! Turbulent Schmidt number
                 if (nclas_chm > 0) diffu_chm(1,1) = Sc_turb_loc
                 write(momod(modul) % lun_outpu,*) 'TURBULENT SCHMIDT NUMBER:',Sc_turb_loc
                 write(momod(modul) % lun_outpu,*) ''

              else if( words(1) == 'KINET' ) then

                 if ( words(2) == 'CANTE' ) then
#ifdef CANTERA
                    write(momod(modul) % lun_outpu,*)'Reading mechanism ',mechanism_path
                    write(momod(modul) % lun_outpu,*)''

                    write(momod(modul) % lun_outpu,*)'# Species=',nclas_chm
                    write(momod(modul) % lun_outpu,*)''

                    write(momod(modul) % lun_outpu,*)'# Removed Species=',Red_spec
                    write(momod(modul) % lun_outpu,*)''

                    do iclas=1,nclas_chm
                       call getSpeciesName(gas_chm,iclas,speci(iclas)%name)
                       write(momod(modul) % lun_outpu,*)'Species ',iclas,'=',speci(iclas)%name
                    end do
                    write(momod(modul) % lun_outpu,*)''
                    write(momod(modul) % lun_outpu,*)'Chemical reactions'
                    write(momod(modul) % lun_outpu,*)''
                    write(momod(modul) % lun_outpu,*)'# Reactions=',nreac_chm
                    write(momod(modul) % lun_outpu,*)''

                    do ireac = 1,nreac_chm
                       call getReactionString(gas_chm, ireac,eq)
                       write(momod(modul) % lun_outpu,*)eq
                    end do
                    write(momod(modul) % lun_outpu,*)''
#endif
                 else if ( words(2) == 'USER' ) then

                 end if

              else if( words(1) == 'TRANS' ) then

                 if( words(2) == 'LEWIS' ) then
                    if( words(3) == 'UNITY' ) then
                       kfl_transport_chm = 1_ip
                       write(momod(modul) % lun_outpu,*)''
                       write(momod(modul) % lun_outpu,*)'Lewis number unity assumption'
                       write(momod(modul) % lun_outpu,*)''
                    elseif( words(3) == 'CONST' ) then
                       kfl_transport_chm = 2_ip
                       write(momod(modul) % lun_outpu,*)''
                       write(momod(modul) % lun_outpu,*)'Lewis number constant assumption'
                       write(momod(modul) % lun_outpu,*)''
                    end if
                 elseif ( words(2) == 'MIXAV' ) then
                    kfl_transport_chm = 3_ip
                    write(momod(modul) % lun_outpu,*)''
                    write(momod(modul) % lun_outpu,*)'Mixtured Average Diffusion'
                    write(momod(modul) % lun_outpu,*)''
                 elseif ( words(2) == 'USER ' ) then
                    kfl_transport_chm = 4_ip
                    call runend('CHEMIC REAPHY: User transport not implemented yet')
                 end if

              !
              ! Reduction Table
              !
              else if( words(1) == 'REDUC' ) then

                 !
                 ! Flamelet model: read thermochemical database
                 !
                 if( words(2) == 'TABLE' ) then
                    if (kfl_pfa_chm /= 4) call runend('CHEMIC REAPHY: Reduction table only for TDAC')
                    !
                    ! TDAC - Dynamic tabulation scheme
                    !
                    if ( kfl_tdac_write_chm > 0_ip ) then
                       call tab_load_file(table_coords, table_tab, dont_read_data=.true.)
                       print *, "Not Reading Table Data"
                    else
                       call tab_load_file(table_coords, table_tab, dont_read_data=.false.)
                       print *, "Reading Table Data"
                    end if
                    table_fw % main_table => table_tab


                    do while(words(1)/='ENDRE')
                       call ecoute('chm_reaphy')
                    end do
                    !
                    ! Print table dimensions
                    !
                    write(momod(modul) % lun_outpu,*)''
                    write(momod(modul) % lun_outpu,*)'---------------------------------------------'
                    write(momod(modul) % lun_outpu,*)'TDAC - REDUCTION:'
                    write(momod(modul) % lun_outpu,*)'---------------------------------------------'
                    write(momod(modul) % lun_outpu,*)''
                    write(momod(modul) % lun_outpu,*)'TABLE DIMENSIONS:'
                    do ii = 1,table_fw % main_table%ndim
                       write(momod(modul) % lun_outpu,'(5X,A5,1X,I5)') table_fw % main_table % coords(ii) % name ,&
                           table_fw % main_table % coords(ii) % leng
                    enddo
                 end if

              !
              ! Reactions
              !
              else if( words(1) == 'TABFR' ) then
                 !
                 ! Table framework
                 !
                 kfl_tab_fw_chm = getint('TABFR',1_ip,'#Index of table framework')
                 table_fw => lookup_fw( kfl_tab_fw_chm )

              else if( words(1) == 'POSTF' ) then
                 !
                 ! Postprocessing table framework
                 !
                 kfl_post_fw_chm = getint('POSTF',1_ip,'#Index of postprocessing table framework')
                 posttable_fw => lookup_fw( kfl_post_fw_chm  )

                 !
                 ! Select names:
                 !
                 do ipostvar=1,posttable_fw % main_table % nvar
                    give_new_name = .false.
                    !
                    ! Check if empty name
                    !
                    if (posttable_fw % main_table % varname(ipostvar) == '') then
                        give_new_name = .true.
                    endif

                    !
                    ! Check if name is repeated
                    !
                    do jpostvar=1,ipostvar-1
                       if (posttable_fw % main_table % varname(ipostvar) == posttable_fw % main_table % varname(jpostvar)) then
                           give_new_name = .true.
                       endif
                    enddo

                    !
                    ! If necessary give new name
                    !
                    if (give_new_name) then
                        if (ipostvar < 10) then
                            posttable_fw % main_table % varname(ipostvar) = 'POS0'//trim(intost(ipostvar))
                        elseif (ipostvar < 100) then
                            posttable_fw % main_table % varname(ipostvar) = 'POS'//trim(intost(ipostvar))
                        else
                            posttable_fw % main_table % varname(ipostvar) = 'PO'//trim(intost(ipostvar))
                        endif
                    endif
                 enddo

              else if( words(1) == 'HRRFR' ) then
                 !
                 ! Framework to look up property that will be volume integrated in the entire domain.
                 ! (Typically heat release rate, but could be arbitrary.)
                 !
                 kfl_hrr_fw_chm = getint('HRRFR',1_ip,'#Index of HRR framework')
                 kfl_hrr_col_chm = getint('HRRCO',1_ip,'#Column in which HRR is tabulated')

              else if( words(1) == 'ZGFRA' ) then
                 !
                 ! Framework to look up property that will be looked up at every
                 ! step.
                 !
                 kfl_zg_fw_chm  = getint('ZGFRA',1_ip,'#Index of ZGRAD framework')
                 kfl_zg_col_chm = getint('ZGCOL',1_ip,'#Column in which ZGRAD is tabulated')

              else if( words(1) == 'ZLIMI' ) then
                 !
                 ! Premixed flamability limits
                 !
                 chm_zmax  = param(2)
                 chm_zmin  = param(3)

              else if( words(1) == 'HYBRI' ) then
                 !
                 ! Activated hybrid model 
                 !
                 kfl_multimod_chm = 1_ip
                 kfl_tab_fw_chm_diff  = getint('DIFFU',1_ip,'#Diffusion Framework')
                 kfl_tab_fw_chm_prem  = getint('PREMI',1_ip,'#Premixed Framwork')
              end if
              call ecoute('chm_reaphy')
           end do properties

           call chm_memphy(10_ip) ! Allocate control variable association
        !
        ! Reading soot model variables
        !
        else if( words(1) == 'SOOTM' ) then

           call ecoute('chm_reaphy')

           !
           !
           !
           write(momod(modul) % lun_outpu,*)'---------------------'
           write(momod(modul) % lun_outpu,*)' SOOT MODEL VARIBLES '
           write(momod(modul) % lun_outpu,*)'---------------------'
           write(momod(modul) % lun_outpu,*) ''

           !
           ! Soot model reading
           !
           soot_model: do while(words(1)/='ENDSO')

              if ( words(1) == 'MODEL' ) then
                 if( words(2) == 'SECTI' ) then
                    kfl_soot_chm = 1
                    write(momod(modul) % lun_outpu,*)''
                    write(momod(modul) % lun_outpu,*)'Soot sectional model activated, model =',kfl_soot_chm
                    write(momod(modul) % lun_outpu,*)''
                 else if( words(2) == 'NONE' ) then
                    kfl_soot_chm = 0
                    write(momod(modul) % lun_outpu,*)''
                    write(momod(modul) % lun_outpu,*)'Soot sectional model not activated, model =',kfl_soot_chm
                    write(momod(modul) % lun_outpu,*)''
                 else
                    call runend('CHM_REAPHY: Only soot sectional method available, be patient!')
                 end if
              else if ( words(1) == 'SECTI' ) then
                 nsect_ssm = getint('SECTI',1_ip,'#Number of sections soot model')
                 write(momod(modul) % lun_outpu,*)'# Number of sections = ',nsect_ssm
                 write(momod(modul) % lun_outpu,*)''

              else if ( words(1) == 'DENSI' ) then
                 RhoC_ssm = getrea('DENSI',1.0_rp,'#Soot particles density')
                 write(momod(modul) % lun_outpu,*)'# Soot particle density = ',RhoC_ssm,'[kg/m3]'
                 write(momod(modul) % lun_outpu,*)''

              else if ( words(1) == 'VMAXS' ) then
                 Vmax_ssm = getrea('VMAXS',1.0_rp,'#Maximum volume soot particles')
                 write(momod(modul) % lun_outpu,*)'# Maximum volume soot particles = ',Vmax_ssm,'[m3]'
                 write(momod(modul) % lun_outpu,*)''

              else if ( words(1) == 'RADPA' ) then
                 RadSoot_ssm = getrea('RADPA',1.0_rp,'#Radiation parameter')
                 write(momod(modul) % lun_outpu,*)'# Radiation parameter = ',RadSoot_ssm
                 write(momod(modul) % lun_outpu,*)''

              else if ( words(1) == 'VALID' ) then
                 if( words(2) == 'CHEM1' ) then
                    kfl_soot_chm = -1
                    write(momod(modul) % lun_outpu,*)'# Validation test with Chem1d '
                    write(momod(modul) % lun_outpu,*)''
                 endif

              else if ( words(1) == 'NUCLE' ) then
                 if( words(2) == 'ON' ) then
                    ID_NUCL = 1
                    write(momod(modul) % lun_outpu,*)'# Nucleation process active '
                    write(momod(modul) % lun_outpu,*)''
                 elseif ( words(2) == 'OFF' ) then
                    ID_NUCL = 0
                 endif

              else if ( words(1) == 'CONDE' ) then
                 if( words(2) == 'ON' ) then
                    ID_COND = 1
                    write(momod(modul) % lun_outpu,*)'# Condensation process active '
                    write(momod(modul) % lun_outpu,*)''
                 elseif ( words(2) == 'OFF' ) then
                    ID_COND = 0
                 endif

              else if ( words(1) == 'COAGU' ) then
                 if( words(2) == 'ON' ) then
                    ID_COAG = 1
                    write(momod(modul) % lun_outpu,*)'# Coagulation process active '
                    write(momod(modul) % lun_outpu,*)''
                 elseif ( words(2) == 'OFF' ) then
                    ID_COAG = 0
                 endif

              else if ( words(1) == 'SURFA' ) then
                 if( words(2) == 'ON' ) then
                    ID_SURF = 1
                    write(momod(modul) % lun_outpu,*)'# Surface growth process active '
                    write(momod(modul) % lun_outpu,*)''
                 elseif ( words(2) == 'OFF' ) then
                    ID_SURF = 0
                 endif

              else if ( words(1) == 'GASCO' ) then
                 if( words(2) == 'ON' ) then
                    gasCoupling_ssm = 1
                 elseif ( words(2) == 'OFF' ) then
                    gasCoupling_ssm = 0
                 else
                    call runend('CHM_REAPHY: Coupling soot with gas phase not specified.')
                 endif

                 write(momod(modul) % lun_outpu,*)'# Coupling with gas phase = ',gasCoupling_ssm
                 write(momod(modul) % lun_outpu,*)''

              else if ( words(1) == 'YKFRA' ) then
                 !
                 ! Index of lookup framework for precursor species
                 !
                 kfl_yk_fw_ssm = getint('YKFRA',1_ip,'#Index of lookup framework for precursor species')

                 write(momod(modul) % lun_outpu,*)'# Index of lookup framework for precursor species = ',kfl_yk_fw_ssm
                 write(momod(modul) % lun_outpu,*)''
              else if ( words(1) == 'FIELD' ) then
                 kfl_field_chm(4) = 1_ip
                 write(momod(modul) % lun_outpu,*)'# Reading Soot Fields'
                 write(momod(modul) % lun_outpu,*)''

              end if

              call ecoute('chm_reaphy')
           end do soot_model

        else if ( words(1) == 'DROPL' ) then
           !
           ! Eulerian droplets identification parameters
           !
           kfl_droplet_id_chm = 1_ip
           call ecoute('chm_reaphy')
           !
           write(momod(modul) % lun_outpu,*) ''
           write(momod(modul) % lun_outpu,*)'---------------------------------'
           write(momod(modul) % lun_outpu,*)' EULERIAN DROPLET IDENTIFICATION '
           write(momod(modul) % lun_outpu,*)'---------------------------------'
           write(momod(modul) % lun_outpu,*) ''
           !
           droplet_id: do while(words(1)/='ENDDR')

              if ( words(1) == 'LEVEL' ) then
                 levelSet_threshold_chm = getrea('LEVEL',0.0_rp,'#Droplet identification parameter')
                 write(momod(modul) % lun_outpu,*) 'LEVEL SET THRESHOLD:',levelSet_threshold_chm
              else if ( words(1) == 'COMPA' ) then
                 droplet_compactness_limit_chm = getrea('COMPA',1.0_rp,'#Droplet identification parameter')
                 write(momod(modul) % lun_outpu,*) 'COMPACTNESS LIMIT:',droplet_compactness_limit_chm
              else if ( words(1) == 'MAXDI' ) then
                 droplet_max_diameter_chm = getrea('MAXDI',0.0_rp,'#Droplet identification parameter')
                 write(momod(modul) % lun_outpu,*) 'MAX DIAMETER:',droplet_max_diameter_chm
              else if ( words(1) == 'HFACT' ) then
                 droplet_h_factor_chm = getrea('HFACT',0.0_rp,'#Droplet identification parameter')
                 write(momod(modul) % lun_outpu,*) 'MESH SIZE FACTOR:',droplet_h_factor_chm
              end if
              call ecoute('chm_reaphy')
           end do droplet_id
        !
        ! Read index of starting fields
        !
        else if (words(1)=='FIELD') then
           if (exists('SPECI')) kfl_field_chm(1) = getint('SPECI',0_ip,'#Species starting fields')
           if (exists('TEMPE')) kfl_field_chm(2) = getint('TEMPE',0_ip,'#Temperature/enthalpy field')
           if (exists('THERM')) kfl_field_chm(2) = getint('THERM',0_ip,'#Temperature/enthalpy field')
           if (exists('ALPHA')) kfl_field_chm(3) = getint('ALPHA',0_ip,'#Alpha field for CMC model')

           if ( kfl_field_chm(2) > 0_ip ) then
              call messages_live('CHEMIC WILL INITIALIZE THERM TO FIELD:'//trim(intost(kfl_field_chm(2))),'WARNING')
           endif

        !
        ! Read species names
        !
        else if (words(1)=='SPECI') then
           kfl_spec_name_chm = getint('SPECI',0_ip,'#Number of species')
           call memory_alloca(mem_modul(1:2,modul),'FIELD_IND_CHM','chm_reaphy',Field_ind_chm,kfl_spec_name_chm)
           call ecoute('chm_reaphy')
           backspace(momod(modul) % lun_pdata)
           nbuff = 1000
           allocate(character(len=nbuff) :: buffer)
           read(momod(modul) % lun_pdata,'(a)', advance='no', size=nbuff, iostat=ioerror) buffer
           if ((ioerror == iostat_eor).or.(ioerror == iostat_end)) ioerror = 0
           if (ioerror /= 0 ) call runend('CHEMIC REAPHY: Cannot read Reduced Species.')
           buffer = adjustl(buffer) ! trim leading whitespace
           buffer = trim(buffer)    ! trim trailing whitespace
           spec_name  = buffer
           nsize_spec_name = len(spec_name, kind=ip)
#ifdef CANTERA
           call cantera_initialization(1_ip,mechanism_path,nsize_mech_name)
           call cantera_trim(spec_name, nsize_spec_name, Field_ind_chm)
#endif
           deallocate(buffer)
        !
        ! Read Lewis Numbers (FORMAT : Species_index Lewis_numner)
        !
        else if (words(1)=='LEWIS') then
           if (kfl_transport_chm /= 2_ip) then
              call runend('CHEMIC REAPHY: Lewis number reading only valid with Constant Lewis transport model')
           else
              call ecoute('chm_reaphy')
              do while( words(1) /= 'ENDLE' )
                 Le_k(int(param(1))) =  param(2)
                 call ecoute('chm_reaphy')
              end do
           end if

        else if( words(1) == 'CMCMO' ) then
           !
           ! Data for CMC model
           !
           call ecoute('chm_reaphy')

           CMC_model: do while(words(1)/='ENDCM')

              !
              ! CMC model: read mixture fraction vector
              !
              ! Read path and length for mixture fraction vector
              if( words(1) == 'MFVEC') then
                 do icoun = 1,2
                    call ecoute('chm_reaphy' )
                    if ( words(1) == 'MFPAT') then   ! Read path to the file
                       !
                       ! Read next line
                       !
                       nbuff = 1000
                       allocate(character(len=nbuff) :: buffer)
                       read(momod(modul) % lun_pdata,'(a)', advance='no', size=nbuff, iostat=ioerror) buffer

                       !
                       ! Check read, disregard errors about buffer being larger
                       ! than the line.
                       !
                       if ((ioerror == iostat_eor) .or. (ioerror == iostat_end)) ioerror = 0
                       if (ioerror /= 0 ) call runend('CHEMIC REAPHY: Cannot read mixture fraction vector path.')

                       !
                       ! Process line
                       !
                       buffer = adjustl(buffer) ! trim leading whitespace
                       buffer = trim(buffer)    ! trim trailing whitespace
                       nsize_mf_vec_chm = len(buffer, kind=ip)
                       mf_vec_path_CMC_chm = buffer(1:nsize_mf_vec_chm)
                       deallocate(buffer)
                    else if ( words(1) == 'MFSIZ') then  ! Read size of the vector
                       nZ_CMC_chm = getint('MFSIZ',1_ip,'#Number of slices in mixture fraction space for CMC conditioning')
                    end if
                 end do

                 if (nsize_mf_vec_chm==0_ip .or. nZ_CMC_chm==1_ip) call runend('CHEMIC REAPHY: bad entries for mixture fraction&
                     & vector definition')


                 ! Read mixture fraction vector

                 write(momod(modul) % lun_outpu,*) 'Reading mixture fraction vector from ', mf_vec_path_CMC_chm
                 call chm_memphy(3_ip)

                 open(unit=1, file=mf_vec_path_CMC_chm, status='old', action='read', iostat=state_f)
                 if (state_f==0) then
                    do line_f = 1,nZ_CMC_chm
                       read(unit=1, fmt=*) Z_CMC_chm(line_f)
                    end do
                 else
                    call runend('CHEMIC REAPHY: Cannot read mixture fraction vector file')
                 end if
                 close(unit=1, iostat=state_f, status='keep')
                 Zs_CMC_chm = Z_CMC_chm(nZ_CMC_chm)

                 ! Do some checks about the vector consistency (Zs<1, monotonically increasing, etc.)
                 if (Z_CMC_chm(1) /= 0.0_rp)        call runend('CHEMIC REAPHY: mixture fraction value for oxidizer different to 0')
                 if (Z_CMC_chm(nZ_CMC_chm) > 1.0_rp)call runend('CHEMIC REAPHY: mixture fraction value for fuel not valid')
                 do line_f = 1,nZ_CMC_chm-1
                    if (Z_CMC_chm(line_f+1) <= Z_CMC_chm(line_f)) call runend('CHEMIC REAPHY: mixture fraction vector not strictly&
                        & monotonically increasing')
                 end do

                 diff_Z_CMC_chm = Z_CMC_chm(2:nZ_CMC_chm) - Z_CMC_chm(1:nZ_CMC_chm-1)

                 write(momod(modul) % lun_outpu,*) 'Values for mixture fraction slices:'
                 do line_f = 1,nZ_CMC_chm
                    write(momod(modul) % lun_outpu,*) line_f, ': ',  Z_CMC_chm(line_f)
                 end do
                 write(momod(modul) % lun_outpu,*)
                 write(momod(modul) % lun_outpu,*) 'Maximum mixture fraction is ', Zs_CMC_chm
                 write(momod(modul) % lun_outpu,*)
                 write(momod(modul) % lun_outpu,*)

                 ! Allocate matrices for boundaries in mixture fraction
                 call chm_memphy(4_ip)

              !
              ! Stoichiometric mixture fraction
              !
              else if( words(1) == 'STOIC' ) then
                 Zstq_CMC_chm = param(1)

              !
              ! Number of mixture fractions and segregation factors for AMC model
              !
              else if( words(1) == 'AMCMF' ) then
                 nZ_AMC_CMC_chm = getint('AMCMF',51_ip,'#Number of mixture fractions for scalar dissipation rate integration&
                     & (AMC model)')

              else if( words(1) == 'AMCSS' ) then
                 nS_AMC_CMC_chm = getint('AMCSS',31_ip,'#Number of segregation factors for scalar dissipation rate integration&
                     & (AMC model)')

              else if( words(1) == 'AMCSM' ) then
                 Smax_AMC_CMC_chm = param(1)

              else if( words(1) == 'STHRE' ) then
                 S_threshold = param(1)

              !
              ! Choose in which dimensions transport will be included in CMC equations
              !
              else if( words(1) == 'TRANS' ) then
                 call ecoute('chm_reaphy')
                 do while (words(1) /= 'ENDTR')
                    if (words(1) == 'PHYSI') then
                       if (words(2) == 'OFF  ')  kfl_trans_phs_spc_CMC_chm = 0_ip
                    end if
                    if (words(1) == 'MIXTU') then
                       if (words(2) == 'OFF  ')  kfl_trans_mxt_spc_CMC_chm = 0_ip
                    end if
                    call ecoute('chm_reaphy')
                 end do

              !
              ! CMC model: define boundary conditions for species and enthalpy
              ! at mixt. frac. = 0 and Zs. This information, if appears, is only
              ! used when starting from scratch.
              !
              else if( words(1) == 'ZBCOX' .or. words(1) == 'ZBCFU' ) then
                 if ( words(1) == 'ZBCOX' ) then
                    index_zbc = 1
                 else
                    index_zbc = 2
                 end if

                 call ecoute('chm_reaphy')
                 if (words(1) == 'TEMPE') then
                    T_bc_CMC_chm(index_zbc) = param(1)
                 else
                    call runend('CHEMIC REAPHY: Temperature not provided at the mixture fraction boundary conditions (CMC model)')
                 end if

                 call ecoute('chm_reaphy')
                 do while(words(1) /= 'ENDZB')
                    iclas = 1
                    found = .false.
                    do while( (iclas <= nclas_chm) .and. (.not. found))
                       ! OJO: Y SI LA ESPECIE TIENE MAS DE 5 LETRAS?
                       if (words(1) == speci(iclas)%name) then
                          react_scal_bc_CMC_chm(iclas+1,index_zbc) = param(1)
                          found = .true.
                       end if
                       iclas = iclas + 1
                    end do
                    if( .not. found ) call runend('CHEMIC REAPHY: Not all the species exist in the mechanism for CMC BC')
                    call ecoute('chm_reaphy')
                 end do

              else if( words(1) == 'CSTFI' ) then
                 call ecoute('chm_reaphy')
                 do while (words(1) /= 'ENDCS')
                    if( words(1) == 'ZAVER' ) then
                       Zavg_const_CMC_chm = param(1)
                       write(momod(modul) % lun_outpu,*)
                       write(momod(modul) % lun_outpu,*) 'Initial value for mixture fraction field equal to ', Zavg_const_CMC_chm
                       write(momod(modul) % lun_outpu,*)
                    else if ( words(1) == 'ZVARI' ) then
                       Zvar_const_CMC_chm = param(1)
                       write(momod(modul) % lun_outpu,*)
                       write(momod(modul) % lun_outpu,*) 'Initial value for mixture fraction variance field equal to ',&
                           Zvar_const_CMC_chm
                       write(momod(modul) % lun_outpu,*)
                    else if ( words(1) == 'XTOTA' ) then
                       Xtot_const_CMC_chm = param(1)
                       write(momod(modul) % lun_outpu,*)
                       write(momod(modul) % lun_outpu,*) 'Initial value for total scalar dissipation rate field equal to ',&
                           Xtot_const_CMC_chm
                       write(momod(modul) % lun_outpu,*)
                    end if
                    call ecoute('chm_reaphy')
                 end do

              else if( words(1) == 'IFWEI' ) then  ! Weigh initial fields
                 if (kfl_rstar == 0 .and. words(2) == 'YES  ') then
                    kfl_weigh_in_eq_CMC_chm = 1_ip
                 else
                    kfl_weigh_in_eq_CMC_chm = 0_ip
                 end if

              else if( words(1) == 'BCINI' ) then ! Method to find boundary and initial conditions
                 if ( words(2) == 'LINEA' ) then  ! Linear relationships
                    kfl_bc_init_method_CMC_chm = 1_ip
                 !!!else if ( words(2) == 'HOMRE' ) then  ! From homogeneous reactors
                 !!!   kfl_bc_init_method_CMC_chm = 2_ip
                 else
                    call runend('CHEMIC REAPHY: only linear relationships are allowed for BCs and initial conditions')
                 end if

              else if( words(1) == 'SOLEN' ) then
                 if ( words(2) == 'NO   ') then
                    kfl_solve_enth_CMC_chm = 0_ip
                 else
                    kfl_solve_enth_CMC_chm = 1_ip
                 end if

              else if( words(1) == 'SPLIT' ) then
                 if ( words(2) == 'NO   ') then
                    kfl_split_CFD_CMC = 0_ip
                 else
                    kfl_split_CFD_CMC = 1_ip
                 end if

              else if( words(1) == 'PDFTR' ) then
                 if  ( words(2) == 'NO   ') kfl_incl_PDF_trans_CMC_chm = 0_ip

              else if( words(1) == 'TYPE ' ) then
                 if (words(2) == 'SPECI') then
                    kfl_bc_alpha_CMC_chm = 0_ip
                 else
                    kfl_bc_alpha_CMC_chm = 1_ip
                 end if

              end if

              call ecoute('chm_reaphy')

           end do CMC_model


        else if ( words(1) == 'EQUAT' ) then
           !
           ! List of equations for equation model
           !
           ngrou_chm   = getint('NGROU',1_ip,'#Number of equation groups')
           nclas_chm   = getint('NEQUA',0_ip,'#Number of equations to solve')

           !
           ! Allocate structures for equations
           !
           call chm_memphy(2_ip)
           call chm_memphy(1_ip)
           call chm_memphy(9_ip)

           if (nclas_chm > 0) diffu_chm(1,1) = Sc_turb_loc

           call ecoute('chm_reaphy')
           !
           !
           equations: do while(words(1)/='ENDEQ')

              if ( words(1) == 'IGROU' ) then
                  !
                  ! GROUP SETTINGS
                  !
                  igrou = getint('IGROU',0_ip,'#Index of group')

                  !
                  ! Name
                  !
                  if (exists('NAME ')) then
                      mixedEq_groups_chm(igrou) % name = getcha('NAME ','NONE ','#Group name')
                  endif

                  !
                  ! Type
                  !
                  if (exists('TYPE ')) then
                      typ = getcha('TYPE ','NONE ','#Group type')
                      call chm_mixedEq_setGroupType(mixedEq_groups_chm(igrou),typ)
                  endif

                  !
                  ! Number of equations
                  !
                  if (exists('NEQUA')) then
                      mixedEq_groups_chm(igrou) % nequa = getint('NEQUA',0_ip,'#Group equation number')
                      !
                      ! Reset indexes whenever a new group is given its number of equations
                      !
                      jj = 0_ip
                      do ii = 1,ngrou_chm
                         mixedEq_groups_chm(ii) % i_start = jj + 1_ip
                         mixedEq_groups_chm(ii) % i_end   = jj + mixedEq_groups_chm(ii) % nequa
                         jj = jj + mixedEq_groups_chm(ii) % nequa
                      enddo
                  endif

                  !
                  ! Turn on thermophoretic diffusion of group
                  !
                  if (exists('TPHDI')) then
                      typ = getcha('TPHDI','ERROR','#Thermophoretic diffusion')
                      if (typ == 'OFF  ') then
                         mixedEq_groups_chm(igrou) % kfl_therm_phor = 0
                      elseif (typ == 'ON   ') then
                         mixedEq_groups_chm(igrou) % kfl_therm_phor = 1
                      endif
                  endif

                  !
                  ! Whether or not to postprocess as named unknown of ALL equations of group
                  !
                  if (exists('DOPOS')) then
                      do iequa = mixedEq_groups_chm(igrou) % i_start, mixedEq_groups_chm(igrou) % i_end
                         mixedEq_eqs_chm(iequa) % kfl_do_post = 1_ip
                      enddo
                  endif

                  !
                  ! Give names to all equations in group
                  !
                  if (exists('EQNAM')) then
                      eqname = getcha('EQNAM','DEFAU','#Equation name')
                      !
                      ! Default is the name of the group
                      !
                      if (eqname=='DEFAU') eqname = mixedEq_groups_chm(igrou) % name
                      !
                      ! The equation name is the given name
                      ! but the last characters are substituated by the equation index WITHIN THE GROUP
                      !
                      jj = mixedEq_groups_chm(igrou) % i_end - mixedEq_groups_chm(igrou) % i_start + 1_ip
                      do iequa = mixedEq_groups_chm(igrou) % i_start, mixedEq_groups_chm(igrou) % i_end
                         ii = iequa - mixedEq_groups_chm(igrou) % i_start + 1_ip
                         if (ii<10) then
                             if (jj < 10) then
                                mixedEq_eqs_chm(iequa) % name = eqname(1:4)//trim(intost(ii))
                             elseif (jj < 100) then
                                mixedEq_eqs_chm(iequa) % name = eqname(1:3)//'0'//trim(intost(ii))
                             elseif (jj < 1000) then
                                mixedEq_eqs_chm(iequa) % name = eqname(1:2)//'00'//trim(intost(ii))
                             endif
                         elseif (ii<100) then
                             if (jj < 100) then
                                mixedEq_eqs_chm(iequa) % name = eqname(1:3)//trim(intost(ii))
                             elseif (jj < 1000) then
                                mixedEq_eqs_chm(iequa) % name = eqname(1:2)//'0'//trim(intost(ii))
                             endif
                         elseif (ii<1000) then
                             mixedEq_eqs_chm(iequa) % name = eqname(1:2)//trim(intost(ii))
                         endif
                      enddo
                  endif

                  !
                  ! Set Lewis number of entire group
                  !
                  if (exists('LEWIS')) then
                      do iequa = mixedEq_groups_chm(igrou) % i_start, mixedEq_groups_chm(igrou) % i_end
                         mixedEq_eqs_chm(iequa) % Lewis = getrea('LEWIS',1.0_rp,'#Lewis number')
                      enddo
                  endif

                  !
                  ! Set constant diffusivity of entire group
                  !
                  if (exists('FIXDI')) then
                      !
                      ! Constant diffusion activated
                      !
                      do iequa = mixedEq_groups_chm(igrou) % i_start, mixedEq_groups_chm(igrou) % i_end
                         mixedEq_eqs_chm(iequa) % kfl_fix_diffusion = 1_ip
                         mixedEq_eqs_chm(iequa) % diffusivity = getrea('FIXDI',0.0_rp,'#Diffusivity')
                      enddo
                  endif


                  !
                  ! Set Equation type to entire group
                  !
                  if (exists('EQTYP')) then
                      typ = getcha('EQTYP','NONE ','#Equation type')
                      do iequa = mixedEq_groups_chm(igrou) % i_start, mixedEq_groups_chm(igrou) % i_end
                         call chm_mixedEq_setEqType(mixedEq_eqs_chm(iequa),typ)
                      enddo
                  endif

                  !
                  ! Set source type to entire group
                  !
                  if (exists('EQSRC')) then
                      typ = getcha('EQSRC','NONE ','#Equation source type')
                      do iequa = mixedEq_groups_chm(igrou) % i_start, mixedEq_groups_chm(igrou) % i_end
                         call chm_mixedEq_setSrcType(mixedEq_eqs_chm(iequa),typ)
                      enddo
                  endif


              else if ( words(1) == 'IEQUA' ) then
                  iequa = getint('IEQUA',0_ip,'#Index of equation')

                  !
                  ! Name
                  !
                  if (exists('NAME ')) then
                      mixedEq_eqs_chm(iequa) % name = getcha('NAME ','NONE ','#Equation name')
                  endif

                  !
                  ! Type
                  !
                  if (exists('TYPE ')) then
                      typ = getcha('TYPE ','NONE ','#Equation type')
                      call chm_mixedEq_setEqType(mixedEq_eqs_chm(iequa),typ)
                  endif

                  !
                  ! Source type
                  !
                  if (exists('TYPES')) then
                      typ = getcha('TYPES','NONE ','#Equation source type')
                      call chm_mixedEq_setSrcType(mixedEq_eqs_chm(iequa),typ)
                  endif

                  !
                  ! Source framework
                  !
                  if (exists('SRCFR')) then
                      mixedEq_eqs_chm(iequa) % kfl_source_fw = getint('SRCFR',0_ip,'#Equation source framework')
                  endif

                  !
                  ! Source Diffusion framework
                  !
                  if (exists('DIFFS')) then
                     mixedEq_eqs_chm(iequa) % kfl_diffsource_fw = getint('DIFFS',0_ip,'#Equation source framework diffusion')
                  endif

                  !
                  ! Source Premixed framework
                  !
                  if (exists('PREMS')) then
                     mixedEq_eqs_chm(iequa) % kfl_premsource_fw  = getint('PREMS',0_ip,'#Equation source framework premixed')
                  endif

                  !
                  ! Source ANN framework
                  !
                  if (exists('SRCAN')) then
                      mixedEq_eqs_chm(iequa) % kfl_source_ann = getint('SRCAN',0_ip,'#Equation source framework')
                  endif

                  !
                  ! Source column
                  !
                  if (exists('SRCCO')) then
                      mixedEq_eqs_chm(iequa) % kfl_source_col = getint('SRCCO',0_ip,'#Equation source column number in framework')
                  endif

                  !
                  ! Source Diffusion coloumn
                  !
                  if (exists('DIFFC')) then
                     mixedEq_eqs_chm(iequa) % kfl_diffsource_col = getint('DIFFC',0_ip,'#Equation source column &
                     & number in framework diffusion')
                  endif

                  !
                  ! Source Premixed coloumn
                  !
                  if (exists('PREMC')) then
                     mixedEq_eqs_chm(iequa) % kfl_premsource_col = getint('PREMC',0_ip,'#Equation source column &
                     & number in framework premixed')
                  endif

                  !
                  ! Source consumption column
                  !
                  if (exists('CSRCO')) then
                      mixedEq_eqs_chm(iequa) % kfl_consum_col = getint('CSRCO',0_ip,'#Equation source consumption column number in&
                          & framework')
                      !
                      ! If Consumption term exists, activate splitting source term
                      !
                      mixedEq_eqs_chm(iequa) % kfl_source_split = 1
                  endif

                  !
                  ! Fix diffusion flag
                  !
                  if (exists('FIXDI')) then
                      !
                      ! Constant diffusion activated
                      !
                      mixedEq_eqs_chm(iequa) % kfl_fix_diffusion = 1
                      mixedEq_eqs_chm(iequa) % diffusivity = getrea('FIXDI',0.0_rp,'#Diffusivity')
                  endif

                  !
                  ! Initial condition type
                  !
                  if (exists('INITY')) then
                      typ = getcha('INITY','NONE ','#Equation initialization type')
                      call chm_mixedEq_setIniType(mixedEq_eqs_chm(iequa),typ)
                  endif

                  !
                  ! Initial condition field: if present, type is automatically set to FIELD
                  !
                  if (exists('INIFI')) then
                      mixedEq_eqs_chm(iequa) % kfl_ini_field = getint('INIFI',0_ip,'#Field for initial condition')
                      call chm_mixedEq_setIniType(mixedEq_eqs_chm(iequa),'FIELD')
                  endif

                  !
                  ! Initial condition value: if present, type is automatically set to CONST
                  !
                  if (exists('INIVA')) then
                      mixedEq_eqs_chm(iequa) % ini_value = getrea('INIVA',0.0_rp,'#Initial value')
                      call chm_mixedEq_setIniType(mixedEq_eqs_chm(iequa),'CONST')
                  endif

                  !
                  ! Index of mean variable
                  !
                  if (exists('IEQME')) then
                      mixedEq_eqs_chm(iequa) % kfl_ieq_mean = getint('IEQME',0_ip,'#Index of mean quantitiy')
                  endif

                  !
                  ! Whether or not to postprocess as named unknown
                  !
                  if (exists('DOPOS')) then
                      mixedEq_eqs_chm(iequa) % kfl_do_post = 1_ip
                  endif

                  !
                  ! Lewis number
                  !
                  if (exists('LEWIS')) then
                      mixedEq_eqs_chm(iequa) % Lewis = getrea('LEWIS',1.0_rp,'#Lewis number')
                  endif


              end if
              call ecoute('chm_reaphy')
           end do equations

           do igrou = 1,ngrou_chm
              !
              ! Find control group
              !
              if (mixedEq_groups_chm(igrou) % kfl_grtype == CHM_GR_CONTROL) then
                  kfl_control_gr_chm = igrou
              endif
           enddo



        end if

     end do physical_problem


     !
     ! Definition of unknowns for chemic with soot model
     !
     !   Y*_k = Y_k + Y_s
     !
     if (kfl_soot_chm /= 0) then
         if (kfl_model_chm == 3 .or. &
            (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1)) then
            !
            ! Finite rate and CMC models
            !
            ngrou_chm = 2
            nspec_chm = nclas_chm
            nclas_chm = nspec_chm + nsect_ssm
         elseif(kfl_model_chm == 2) then
            !
            ! Mixed equations model
            !
            nspec_chm = nclas_chm - nsect_ssm
         endif

         !
         ! Copy variables to SSM copies
         !
         nspec_ssm = nspec_chm
         nclas_ssm = nclas_chm
     else
         nspec_chm = nclas_chm
     endif

     !
     ! Allocate equation structure, even if it is not explicitly used in the input file
     !
     if (kfl_model_chm /= 2) then
         call chm_memphy(9_ip)
         !
         ! Set start and end indeces
         !
         mixedEq_groups_chm(1) % nequa   = nspec_chm
         mixedEq_groups_chm(1) % i_start = 1
         mixedEq_groups_chm(1) % i_end   = nspec_chm
         !
         ! Set types
         !
         call chm_mixedEq_setGroupType(mixedEq_groups_chm(1),'REACT')

         !
         ! Set equations
         !
         do iequa = mixedEq_groups_chm(1) % i_start, mixedEq_groups_chm(1) % i_end
            call chm_mixedEq_setEqType(mixedEq_eqs_chm(iequa),'REACT')
            if ( kfl_field_chm(1) > 0 ) then
               mixedEq_eqs_chm(iequa) % kfl_ini_type = CHM_INI_OLDWAY
            endif
         enddo

         !
         ! SET SECTIONAL GROUP
         !
         if (kfl_soot_chm /= 0) then
            mixedEq_groups_chm(2) % nequa   = nsect_ssm
            mixedEq_groups_chm(2) % i_start = nspec_chm + 1
            mixedEq_groups_chm(2) % i_end   = nclas_chm

            call chm_mixedEq_setGroupType(mixedEq_groups_chm(2),'SECTI')
            mixedEq_groups_chm(2) % kfl_therm_phor = 1

            do iequa = mixedEq_groups_chm(2) % i_start, mixedEq_groups_chm(2) % i_end
               call chm_mixedEq_setEqType(mixedEq_eqs_chm(iequa),'SECTI')
               call chm_mixedEq_setSrcType(mixedEq_eqs_chm(iequa),'DSM  ')
               if (kfl_field_chm(4) > 0_ip) then
                  mixedEq_eqs_chm(iequa) % kfl_ini_type = CHM_INI_OLDWAY
               end if
            enddo
         endif
     endif

     !
     ! Information and actions for CMC model
     !
     if( kfl_model_chm == 4 ) then
        call chm_initial_actions_reaphy_CMC
     end if

  end if

end subroutine chm_reaphy
