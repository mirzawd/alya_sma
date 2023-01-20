!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_inivar(itask)
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_inivar
  ! NAME
  !    chm_inivar
  ! DESCRIPTION
  !    This routine initializes some variables
  ! USES
  ! USED BY
  !    chm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_master,                    only : momod, conce, ISLAVE, kfl_timco, kfl_timei, nspec, prthe, modul, postp, solve, speci
  use def_chemic,                    only : ADR_typ, typ_lookup_framework, table_fw, alpha_val_spec_CMC_chm,&
                                            av_Z_flux_chm, avd32_chm, avden_chm, avL2_chm, avL_chm, avmsk_chm, avposttab_chm,&
                                            avS0_chm, avS_chm, avtim_chm, avxYr_chm, avxYs_chm, avxZr_chm, avxZs_chm, avY2_chm,&
                                            avY_chm, avYv_chm, avZ2_chm, avZ_chm, avZv_chm, centroid_drop_chm, bemol_chm,&
                                            coeff_cp_k, compactness2_drop_chm, Corr_chm, cputi_chm, d32_chm, diameter_drop_chm,&
                                            diffu_chm, dt_chm, dt_rho_chm, dtinv_chm, dtmat_chm, elem_c, elem_h, elem_n, elem_o,&
                                            enthalp_CMC_chm, enthalpy_transport_nodes, grad_phi_levSet_chm, grad_phic_levSet_chm,&
                                            grad_T, grad_Yk, hrr_avg_chm, hrr_chm, iclaf_chm, iclai_chm, ittot_chm,&
                                            kfl_annfw_has_sources, kfl_control_gr_chm, kfl_cpCoefHT_end_tab_index_chm,&
                                            kfl_cpCoefHT_tab_index_chm, kfl_cpCoefLT_end_tab_index_chm, kfl_cpCoefLT_tab_index_chm,&
                                            kfl_DTRho_tab_index_chm, kfl_ellen_chm, kfl_fixbo_chm, kfl_fixno_chm, kfl_funno_chm,&
                                            kfl_funtn_chm, kfl_fw_has_sources, bvess_chm, kfl_grdif_chm, kfl_hrr_fw_chm,&
                                            kfl_icmean_chm, kfl_icvar_chm, kfl_izmean_chm, kfl_izvar_chm, kfl_k_tab_index_chm,&
                                            kfl_max_nvar_ann_in_chm, kfl_max_nvar_ann_out_chm, kfl_max_nvar_lookup_chm,&
                                            kfl_max_src_annfw_chm, kfl_max_srcfw_chm, kfl_min_src_annfw_chm, kfl_min_srcfw_chm,&
                                            kfl_model_chm, kfl_mu_tab_index_chm, kfl_post_fw_chm, kfl_shock_chm,&
                                            kfl_solve_cond_CMC_chm, kfl_stabi_chm, kfl_k_tab_index_chm, kfl_tab_fw_chm,&
                                            kfl_taust_chm, kfl_tiacc_chm, kfl_tiaor_chm, kfl_tibub_chm, kfl_timei_chm,&
                                            kfl_tisch_chm, kfl_ufpv_chm, kfl_varYc_chm, kfl_varZ_chm, kfl_W_tab_index_chm,&
                                            lap_phi_levSet_chm, Le_k, max_ann_fw, max_lookup_fw, mixedEq_eqs_chm,&
                                            mixedEq_groups_chm, mixfr_chm, nclas_chm, ncomp_chm, nspec_chm, nvar_CMC_chm, phi_chm,&
                                            kfl_T_tab_index_chm, prog_var_chm, prthe_chm, React_ind, shock_chm, Sigm0_chm,&
                                            Sigma_chm, src_chm, sum_reac_chm, temp_CMC_chm, volume_cluster_chm, volume_drop_chm,&
                                            W_k, xYr_chm, xYs_chm, xZr_chm, xZs_chm, Y_k_n, Yk_CMC_chm, zgradmax_chm,&
                                            kfl_field_chm, kfl_fw_src_equa_list, kfl_annfw_src_equa_list, ADR_chm, staco_chm,&
                                            kfl_tab_fw_chm_diff, kfl_tab_fw_chm_prem, kfl_fw_has_hybrid_sources,kfl_multimod_chm
  use def_kintyp,                    only : ip, rp
  use mod_ADR,                       only : ADR_initialize_type
  use mod_ADR,                       only : ADR_check_and_compute_data
  use mod_arrays,                    only : arrays_register
  use def_kermod,                    only : lookup_fw
  use def_kermod,                    only : ann_fw
  use mod_chm_rk_explicit,           only : chm_rk_explicit_initialization
  use mod_chm_mixedEq,               only : CHM_EQ_Z
  use mod_chm_mixedEq,               only : CHM_EQ_ZVAR
  use mod_chm_mixedEq,               only : CHM_EQ_ZZ
  use mod_chm_mixedEq,               only : CHM_EQ_YC
  use mod_chm_mixedEq,               only : CHM_EQ_YCVAR
  use mod_chm_mixedEq,               only : CHM_EQ_YCYC
  use mod_chm_mixedEq,               only : CHM_SRC_TABLE
  use mod_chm_mixedEq,               only : CHM_SRC_ANN
  use mod_chm_mixedEq,               only : chm_mixedEq_EqTyp2kfl
  use mod_chm_mixedEq,               only : chm_mixedEq_setEqType
  use mod_chm_mixedEq,               only : chm_mixedEq_setGroupType
  use mod_chm_mixedEq,               only : chm_mixedEq_setSrcType
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ispec,iclas
  integer(ip)             :: nclas_max
  integer(ip), save       :: nbou_set_vars
  type(ADR_typ)           :: ADR_read  ! ADR type for reading

  integer(ip)                         :: itab,idimt,loc_kfl_typ,iequa,ifw,ivari,ii,jj
  integer(ip), pointer                :: kfl_fw_loc
  integer(ip), pointer                :: kfl_fw_control_loc(:)
  type(typ_lookup_framework), pointer :: fw_loc

  external                            :: soldef

  select case(itask)

  case(0_ip)
     !
     ! Initialize modules
     !
     call chm_rk_explicit_initialization()
     !
     ! Postprocess
     !
     call arrays_register((/'CONCE','SCALA','NPOIN','PRIMA'/),conce,ENTITY_POSITION=1_ip,TIME_POSITION=3_ip,&
         EXCLUDE_RST=(/2_ip,3_ip/))
     call arrays_register((/'SOURC','SCALA','NPOIN','SECON'/))    ! Source terms combustion
     call arrays_register((/'VISCO','SCALA','NPOIN','SECON'/))    ! Viscosity
     call arrays_register((/'SPHEA','SCALA','NPOIN','SECON'/))    ! Specific heat
     call arrays_register((/'CONDU','SCALA','NPOIN','SECON'/))    ! Heat conductivity
     call arrays_register((/'ENTHA','SCALA','NPOIN','SECON'/))    ! Enthalpy
     call arrays_register((/'DIVEN','SCALA','NPOIN','SECON'/))    ! Enthalpy transport source term

     call arrays_register((/'SUMCO','SCALA','NPOIN','SECON'/))    ! Sum of concentration
     call arrays_register((/'MOLEC','SCALA','NPOIN','SECON'/))    ! Molecular weight
     call arrays_register((/'TEMPE','SCALA','NPOIN','SECON'/))    ! Temperature for low-mach Flamelet combustion model
     call arrays_register((/'AVY  ','SCALA','NPOIN','SECON'/))    ! Averaged reaction progress variable Yc or C
     call arrays_register((/'AVYV ','SCALA','NPOIN','SECON'/))    ! Averaged variance of reaction progress variable Yc or C
     call arrays_register((/'AVZV ','SCALA','NPOIN','SECON'/))    ! Averaged variance of mixture fraction Z

     call arrays_register((/'AVZ  ','SCALA','NPOIN','SECON'/))    ! Averaged mixture fraction Z
     call arrays_register((/'AVZ2 ','SCALA','NPOIN','SECON'/))    ! Averaged squared of mixture fraction
     call arrays_register((/'SPEC1','SCALA','NPOIN','SECON'/))    ! Species post-processing Flamelet model, radiation
     call arrays_register((/'SPEC2','SCALA','NPOIN','SECON'/))    ! Species post-processing Flamelet model, radiation
     call arrays_register((/'RADIA','SCALA','NPOIN','SECON'/))    ! Radiation source term Flamelet model
     call arrays_register((/'ZGRMA','SCALA','NPOIN','SECON'/))    ! Maximum mixture fraction variance in the direction
                                                                  !  perpendicular to the flame
     call arrays_register((/'PHI  ','SCALA','NPOIN','SECON'/))    ! Weighting parameter for the hybrid model
     call arrays_register((/'XYR  ','SCALA','NPOIN','SECON'/))    ! Scalar dissipation rate of Yc (resolved part)
     call arrays_register((/'XZR  ','SCALA','NPOIN','SECON'/))    ! Scalar dissipation rate of Z  (resolved part)
     call arrays_register((/'XYS  ','SCALA','NPOIN','SECON'/))    ! Scalar dissipation rate of Yc (subgrid part)
     call arrays_register((/'XZS  ','SCALA','NPOIN','SECON'/))    ! Scalar dissipation rate of Z  (subgrid part)
     call arrays_register((/'AVXYR','SCALA','NPOIN','SECON'/))    ! Average scalar dissipation rate of Yc (resolved part)
     call arrays_register((/'AVXZR','SCALA','NPOIN','SECON'/))    ! Average scalar dissipation rate of Z  (resolved part)
     call arrays_register((/'AVXYS','SCALA','NPOIN','SECON'/))    ! Average scalar dissipation rate of Yc (subgrid part)
     call arrays_register((/'AVXZS','SCALA','NPOIN','SECON'/))    ! Average scalar dissipation rate of Z  (subgrid part)
     call arrays_register((/'AVY2 ','SCALA','NPOIN','SECON'/))    ! Average progress variable squared Yc*Yc or C*C
     call arrays_register((/'AVL  ','SCALA','NPOIN','SECON'/))    ! Average liquid volume fraction phi_L
     call arrays_register((/'AVL2 ','SCALA','NPOIN','SECON'/))    ! Average liquid volume fraction squared phi_L*phi_L
     call arrays_register((/'AVS  ','SCALA','NPOIN','SECON'/))    ! Average interface surface density Sigma
     call arrays_register((/'AVS0 ','SCALA','NPOIN','SECON'/))    ! Average interface surface density Sigma_0
     call arrays_register((/'AVD32','SCALA','NPOIN','SECON'/))    ! Average Sauter mean diameter
     call arrays_register((/'AVDEN','SCALA','NPOIN','SECON'/))    ! Average density
     call arrays_register((/'SIGMA','SCALA','NPOIN','SECON'/))    ! Interface surface density Sigma
     call arrays_register((/'SIGM0','SCALA','NPOIN','SECON'/))    ! Interface surface density Sigma_0
     call arrays_register((/'D32  ','SCALA','NPOIN','SECON'/))    ! Sauter mean diameter
     call arrays_register((/'GRADY','VECTO','NPOIN','SECON'/))    ! Gradient species mass fractions
     call arrays_register((/'HTRAN','VECTO','NPOIN','SECON'/))    ! Enthalpy transport by diffusion
     call arrays_register((/'ELEMH','SCALA','NPOIN','SECON'/))    ! Elemental fraction H
     call arrays_register((/'ELEMO','SCALA','NPOIN','SECON'/))    ! Elemental fraction O
     call arrays_register((/'ELEMC','SCALA','NPOIN','SECON'/))    ! Elemental fraction C
     call arrays_register((/'MASSK','SCALA','NPOIN','SECON'/))    ! Mass source from partis
     call arrays_register((/'ENVIC','SCALA','NPOIN','SECON'/))    ! Entropy viscosity maximum among species
     call arrays_register((/'SCONC','SCALA','NPOIN','SECON'/))    ! Scaled control variables
     call arrays_register((/'NREAC','SCALA','NPOIN','SECON'/))    ! Sum of Reactions (finite Rate)
     call arrays_register((/'MIXFR','SCALA','NPOIN','SECON'/))    ! Mixture Fraction (finite rate)
     call arrays_register((/'HRR  ','SCALA','NPOIN','SECON'/))    ! Instantaneous heat release
     call arrays_register((/'AVHRR','SCALA','NPOIN','SECON'/))    ! Averaged heat release
     call arrays_register((/'AVMSK','SCALA','NPOIN','SECON'/))    ! Average mass source term from spray
     call arrays_register((/'POSTT','SCALA','NPOIN','SECON'/))    ! Posttab properties
     call arrays_register((/'AVPOT','SCALA','NPOIN','SECON'/))    ! Average posttab properties
     call arrays_register((/'MASCN','MATRI','NPOIN','PRIMA'/),Yk_CMC_chm,ENTITY_POSITION=2_ip)  ! Conditional mass fractions
                                                                                                !  for CMC model
     call arrays_register((/'MASUN','SCALA','NPOIN','SECON'/))    ! Unconditional mass fractions for CMC model
     call arrays_register((/'ENTCN','MATRI','NPOIN','PRIMA'/),enthalp_CMC_chm,ENTITY_POSITION=2_ip)  ! Conditional enthalpy
                                                                                                     !  for CMC model
     call arrays_register((/'ENTUN','SCALA','NPOIN','SECON'/))    ! Unconditional enthalpy for CMC model
     call arrays_register((/'TEMCN','MATRI','NPOIN','PRIMA'/),temp_CMC_chm,ENTITY_POSITION=2_ip)  ! Conditional temperature
                                                                                                  !  for CMC model
     call arrays_register((/'TEMUN','SCALA','NPOIN','SECON'/))    ! Unconditional temperature for CMC model
     call arrays_register((/'SRCCN','SCALA','NPOIN','SECON'/))    ! Conditional chemical source terms for CMC model
     call arrays_register((/'SRCUN','SCALA','NPOIN','SECON'/))    ! Unconditional chemical source terms for CMC model
     call arrays_register((/'VTHER','VECTO','NPOIN','SECON'/))    ! Thermophoretic velocity soot model
     call arrays_register((/'TTABL','SCALA','NPOIN','SECON'/))    ! Tabulated temperature
     call arrays_register((/'NAMED','SCALA','NPOIN','SECON'/))    ! Named variable
     call arrays_register((/'PDF  ','SCALA','NPOIN','SECON'/))    ! Unconditional chemical source terms for CMC model
     call arrays_register((/'ALPHA','MATRI','NPOIN','PRIMA'/),alpha_val_spec_CMC_chm,ENTITY_POSITION=2_ip)  ! Alpha values
                                                                                                            !  for HR for CMC
                                                                                                            !  models
     call arrays_register((/'YKSSM','SCALA','NPOIN','SECON'/))    ! Tabulated species for sectional soot model
     call arrays_register((/'QNUCT','SCALA','NPOIN','SECON'/))    ! Discrete sectional method: total soot source due to
                                                                  !  nucleation
     call arrays_register((/'QCOAT','SCALA','NPOIN','SECON'/))    ! Discrete sectional method: total soot source due to
                                                                  !  coagulation
     call arrays_register((/'QCONT','SCALA','NPOIN','SECON'/))    ! Discrete sectional method: total soot source due to
                                                                  !  condensation
     call arrays_register((/'QSURT','SCALA','NPOIN','SECON'/))    ! Discrete sectional method: total soot source due to
                                                                  !  surface growth
     call arrays_register((/'QTOTT','SCALA','NPOIN','SECON'/))    ! Discrete sectional method: total soot source due to
                                                                  !  all effects
     call arrays_register((/'DTRHO','SCALA','NPOIN','SECON'/))    ! Dt*rho
     call arrays_register((/'AVZFL','VECTO','NPOIN','PRIMA'/),av_Z_flux_chm,ENTITY_POSITION=2_ip)  ! Average Z mass flux
     call arrays_register((/'AVYCN','SCALA','NPOIN','SECON'/))    ! Temp. aveg. conditional mass fractions for CMC model
     call arrays_register((/'AVYUN','SCALA','NPOIN','SECON'/))    ! Temp. aveg. unconditional mass fractions for CMC model
     call arrays_register((/'AVHCN','SCALA','NPOIN','SECON'/))    ! Temp. aveg. conditional enthalpy for CMC model
     call arrays_register((/'AVHUN','SCALA','NPOIN','SECON'/))    ! Temp. aveg. unconditional enthalpy for CMC model
     call arrays_register((/'AVTCN','SCALA','NPOIN','SECON'/))    ! Temp. aveg. conditional temperature for CMC model
     call arrays_register((/'AVTUN','SCALA','NPOIN','SECON'/))    ! Temp. aveg. unconditional temperature for CMC model
     call arrays_register((/'AVVIS','SCALA','NPOIN','SECON'/))    ! Temp. aveg. unconditional viscosity for CMC model
     call arrays_register((/'SVF  ','SCALA','NPOIN','SECON'/))    ! Discrete sectional method: soot volume fraction
     call arrays_register((/'NDEN ','SCALA','NPOIN','SECON'/))    ! Discrete sectional method: soot volume fraction
     call arrays_register((/'AVNAM','MATRI','NPOIN','SECON'/))    ! Averaged Named variable
     call arrays_register((/'TCHMI','SCALA','NPOIN','SECON'/))    ! Time for chemical integration
     call arrays_register((/'PRODU','SCALA','NPOIN','SECON'/))    ! Production part of reactive source term
     call arrays_register((/'CONSU','SCALA','NPOIN','SECON'/))    ! Consumption part of reactive source term
     call arrays_register((/'FLAME','SCALA','NPOIN','SECON'/))    ! Flame index 
     call arrays_register((/'ZGRAD','SCALA','NPOIN','SECON'/))    ! Z gradient looked up from diffusion table

     !
     ! Nodal set variables
     !
     postp(1) % wonse (1)    = 'CONCE'

     !
     ! Elemental set variables
     !
     postp(1) % woese (1)    = 'MASS ' !  int_V rho dV
     postp(1) % woese (2)    = 'HEATR' !  int_V omega_h dV, assuming omega_h is in J/(m3 s)
     postp(1) % woese (3:6)  = 'CONCE' !  int_V rho Y_k dV

     !
     ! Boundary set variables
     !
     postp(1) % wobse (1)    = 'MASS ' ! int_S rho*u.n ds
     postp(1) % wobse (2:9)  = 'CONCE' ! <Yk> = int_S rho*u* Y_k dS / int_S rho*u dS, for k = 1,...,8
     nbou_set_vars = 1                 ! Need to set this for later reentry when we know nspec

     !
     ! Witness variables
     !
     postp(1) % wowit (1:8)  = 'CONCE' ! Species 1 to 8 (update the maximum number of species here)

     postp(1) % wowit (8+1)  = 'XZR  ' ! Resolved scalar dissipation of mixture fraction
     postp(1) % wowit (8+2)  = 'XZS  ' ! Subgrid scalar dissipation of mixture fraction
     postp(1) % wowit (8+3)  = 'XYR  ' ! Resolved scalar dissipation of progress variable
     postp(1) % wowit (8+4)  = 'XYS  ' ! Subgrid scalar dissipation of progress variable
     postp(1) % wowit (8+5)  = 'D32  ' ! Sauter mean diameter Eulerian atomization
     postp(1) % wowit (8+6)  = 'SIGMA' ! Liquid surface density S = S_0 + S'
     postp(1) % wowit (8+7)  = 'SIGM0' ! Mean liquid surface density S_0

     postp(1) % wowit ((8+8):(8+16))  = 'SCONC' ! Scaled control variable
     postp(1) % wowit (8+17)          = 'HRR'   ! Heat release rate (Finite rate)
     postp(1) % wowit (8+18)          = 'NREAC' ! Number of reactions (Finite rate)
     postp(1) % wowit (8+19)          = 'DENSI' ! Density
     !
     ! Solver
     !
     call soldef(-2_ip)
     solve(1) % wprob       = 'CONCENTRATION'
     solve(1) % kfl_solve   = 1

     solve(2) % wprob       = 'CONSISTENT'           ! Consistent matrix
     solve(2) % kfl_solve   = 1
     solve(2) % kfl_iffix   = 2

     !
     ! Others
     !
     cputi_chm = 0.0_rp  ! CPU times
     dtmat_chm = 0.0_rp  ! Matrix-based time step

     !
     ! Lookup frameworks
     !
     kfl_fw_has_sources = 0_ip
     kfl_annfw_has_sources = 0_ip

     !
     ! Nullify pointers
     !
     nullify(kfl_fixno_chm)
     nullify(kfl_fixbo_chm)
     nullify(kfl_funno_chm)
     nullify(kfl_funtn_chm)
     nullify(bvess_chm)

     !
     ! Lamianr and turbulent diffusion
     !
     nullify(Le_k)
     nullify(diffu_chm)

     !
     ! Flamelet combustion model
     !
     nullify(avY_chm)
     nullify(avYv_chm)
     nullify(avZ_chm)
     nullify(avZv_chm)
     nullify(avZ2_chm)
     nullify(avY2_chm)

     nullify(xYr_chm)
     nullify(xYs_chm)
     nullify(xZr_chm)
     nullify(xZs_chm)

     nullify(avxYr_chm)
     nullify(avxYs_chm)
     nullify(avxZr_chm)
     nullify(avxZs_chm)
     nullify(avposttab_chm)

     nullify(zgradmax_chm)
     nullify(phi_chm)

     !
     ! MixedEq model
     !
     nullify(mixedEq_groups_chm)
     nullify(mixedEq_eqs_chm)

     !
     ! ELSA model
     !
     nullify(avL_chm)
     nullify(avL2_chm)
     nullify(avS_chm)
     nullify(avS0_chm)
     nullify(avd32_chm)
     nullify(Sigma_chm)
     nullify(Sigm0_chm)
     nullify(d32_chm)
     nullify(lap_phi_levSet_chm)
     nullify(grad_phi_levSet_chm)
     nullify(grad_phic_levSet_chm)

     !
     ! Droplet Identification
     !
     nullify(volume_drop_chm)
     nullify(diameter_drop_chm)
     nullify(centroid_drop_chm)
     nullify(compactness2_drop_chm)
     nullify(volume_cluster_chm)

     !
     ! Others
     !
     nullify(avden_chm)
     nullify(av_Z_flux_chm)
     nullify(avmsk_chm)
     nullify(dt_rho_chm)
     nullify(dt_chm)

     nullify(grad_Yk)
     nullify(grad_T)
     nullify(enthalpy_transport_nodes)
     nullify(coeff_cp_k)
     nullify(W_k)
     nullify(Y_k_n)
     nullify(React_ind)
     nullify(elem_c)
     nullify(elem_h)
     nullify(elem_o)
     nullify(elem_n)
     nullify(src_chm)
     nullify(Corr_chm)
     nullify(mixfr_chm)
     nullify(prog_var_chm)
     nullify(sum_reac_chm)
     nullify(hrr_chm)
     nullify(hrr_avg_chm)

  case(2_ip)
     !
     ! AFTER ALLOCATING ARRAYS IN parall(2)
     !
     if ( kfl_model_chm == 1 ) then
        !
        ! Set equations
        !
        kfl_control_gr_chm = 1_ip
        call chm_mixedEq_setGroupType(mixedEq_groups_chm(kfl_control_gr_chm),'CONTR')
        mixedEq_groups_chm(kfl_control_gr_chm) % nequa   = nclas_chm
        mixedEq_groups_chm(kfl_control_gr_chm) % i_start = 1_ip
        mixedEq_groups_chm(kfl_control_gr_chm) % i_end   = nclas_chm
        !
        ! Control variable names for Flamlet model
        !
        speci(1_ip)%name = 'CMEAN'
        mixedEq_eqs_chm(1_ip) % name             = 'CMEAN'
        if (kfl_tab_fw_chm > 0) then
           call chm_mixedEq_setSrcType(mixedEq_eqs_chm(1_ip),'TABLE')
           mixedEq_eqs_chm (1_ip) % kfl_source_fw  = kfl_tab_fw_chm
           mixedEq_eqs_chm (1_ip) % kfl_source_col = 1_ip
        else
           call chm_mixedEq_setSrcType(mixedEq_eqs_chm(1_ip),'OFF  ')
        endif

        speci(2_ip)%name = 'CVAR'
        mixedEq_eqs_chm(2_ip) % name       = 'PASSI'
        call chm_mixedEq_setSrcType(mixedEq_eqs_chm(2_ip),'OFF  ')
        if (kfl_varYc_chm /= 0_ip) then
           if (kfl_varYc_chm == 1_ip) then
              mixedEq_eqs_chm(2_ip) % name = 'YCVAR'
           elseif (kfl_varYc_chm == 2_ip) then
              mixedEq_eqs_chm(2_ip) % name = 'YCYC'
           endif
           call chm_mixedEq_setSrcType(mixedEq_eqs_chm(2_ip),'TABLE')
           mixedEq_eqs_chm (2_ip) % kfl_source_fw  = kfl_tab_fw_chm
           mixedEq_eqs_chm (2_ip) % kfl_source_col = 17_ip
           mixedEq_eqs_chm (2_ip) % kfl_ieq_mean   = 1_ip
        endif

        if (nclas_chm > 2_ip) then
           speci(3_ip)%name = 'ZMEAN'
           mixedEq_eqs_chm(3_ip) % name = 'ZMEAN'
           call chm_mixedEq_setSrcType(mixedEq_eqs_chm(3_ip),'OFF  ')

           speci(4_ip)%name = 'ZVAR'
           mixedEq_eqs_chm(4_ip) % name = 'PASSI'
           if (kfl_varZ_chm /= 0_ip) then
              if (kfl_varZ_chm == 1_ip) then
                 mixedEq_eqs_chm(4_ip) % name = 'ZVAR'
              elseif (kfl_varZ_chm == 2_ip) then
                 mixedEq_eqs_chm(4_ip) % name = 'ZZ'
              endif
              mixedEq_eqs_chm(4_ip) % kfl_ieq_mean = 3_ip
           endif
           call chm_mixedEq_setSrcType(mixedEq_eqs_chm(4_ip),'OFF  ')
        endif

        !
        ! Control variable names for UFPV model
        !
        if (kfl_ufpv_chm > 0) then
           speci(2_ip)%name = 'CHIST'
           mixedEq_eqs_chm(2_ip) % name = 'NONE'
           call chm_mixedEq_setSrcType(mixedEq_eqs_chm(1_ip),'OFF  ')
           mixedEq_eqs_chm (2_ip) % kfl_source_fw = 0_ip

           !
           ! NO VARIANCE FOR Yc
           !
           kfl_varYc_chm = 0_ip
        endif

        do iequa = mixedEq_groups_chm(kfl_control_gr_chm) % i_start, mixedEq_groups_chm(kfl_control_gr_chm) % i_end
           !
           ! Set equation types
           !
           call chm_mixedEq_setEqType(mixedEq_eqs_chm(iequa), mixedEq_eqs_chm(iequa) % name)
           !
           ! Set initial fields
           !
           if (kfl_field_chm(1) > 0) &
              mixedEq_eqs_chm(iequa) % kfl_ini_field = kfl_field_chm(1) + (iequa-1_ip)
        enddo

     endif


     !
     ! Get lookup framework indexes for main table and postprocessing table
     !
     if (kfl_control_gr_chm > 0) then
        do itab = 1,5
           nullify(kfl_fw_loc)
           nullify(kfl_fw_control_loc)
           nullify(fw_loc)
           select case(itab)
           case(1)
               kfl_fw_loc              => kfl_tab_fw_chm
           case(2)
               kfl_fw_loc              => kfl_post_fw_chm
           case(3)
               kfl_fw_loc              => kfl_hrr_fw_chm
           case(4)
               kfl_fw_loc              => kfl_tab_fw_chm_diff 
           case(5)
               kfl_fw_loc              => kfl_tab_fw_chm_prem
           end select

           if (kfl_fw_loc > 0) then
              fw_loc              => lookup_fw( kfl_fw_loc )
              kfl_fw_control_loc  => fw_loc % kfl_chm_control
           endif

           !
           ! Check which unknown is to be passed to the lookup frameworks
           !
           if (kfl_fw_loc > 0) then
              do idimt=1,fw_loc % main_table % ndim
                 !
                 ! Go through table dimensions and see if
                 ! the name of the dimension corresponds to
                 ! some of the control variables
                 !
                 loc_kfl_typ = chm_mixedEq_EqTyp2kfl(fw_loc % main_table % coords(idimt) % name)
                 !
                 ! Deal with SQUARE of variable instead of variance
                 !
                 if (fw_loc % kfl_scale(idimt) <= -100) then
                     if (loc_kfl_typ == CHM_EQ_ZVAR)  loc_kfl_typ = CHM_EQ_ZZ
                     if (loc_kfl_typ == CHM_EQ_YCVAR) loc_kfl_typ = CHM_EQ_YCYC
                 endif

                 if (loc_kfl_typ /= 0_ip) then
                    do iequa = mixedEq_groups_chm(kfl_control_gr_chm) % i_start, mixedEq_groups_chm(kfl_control_gr_chm) % i_end
                       !
                       ! Find matching equation
                       !
                       if (loc_kfl_typ == mixedEq_eqs_chm(iequa) % kfl_eqtype) then
                           kfl_fw_control_loc(idimt) = iequa
                       endif
                    enddo
                 else
                    select case (fw_loc % main_table % coords(idimt) % name)
                    case ('CHIST')
                        kfl_fw_control_loc(idimt) = -2_ip
                    case ('IMEAN','I    ')
                        kfl_fw_control_loc(idimt) = -1_ip
                        if (fw_loc % kfl_scale(idimt) /= 3 .and. fw_loc % kfl_scale(idimt) /= -1 ) then
                            !
                            ! I or IMEAN coordinate of table is not constant or turned off
                            ! Thus enthalpy needs to be gathered for lookup
                            !
                            fw_loc % kfl_needs_enthalpy = 1_ip
                        endif
                    end select
                 endif
              enddo
           endif
        enddo
     endif



     !
     ! set index of some typical control variables
     !
     kfl_izmean_chm = 0_ip
     kfl_izvar_chm  = 0_ip
     kfl_icmean_chm = 0_ip
     kfl_icvar_chm  = 0_ip

     if (kfl_control_gr_chm > 0) then
        do iequa = mixedEq_groups_chm(kfl_control_gr_chm) % i_start, mixedEq_groups_chm(kfl_control_gr_chm) % i_end
           !
           ! Find matching equation
           !
           if     (CHM_EQ_Z     == mixedEq_eqs_chm(iequa) % kfl_eqtype) then
               kfl_izmean_chm = iequa
           elseif (CHM_EQ_ZVAR  == mixedEq_eqs_chm(iequa) % kfl_eqtype .or. &
                  &CHM_EQ_ZZ    == mixedEq_eqs_chm(iequa) % kfl_eqtype) then
               kfl_izvar_chm  = iequa
           elseif (CHM_EQ_YC    == mixedEq_eqs_chm(iequa) % kfl_eqtype) then
               kfl_icmean_chm = iequa
           elseif (CHM_EQ_YCVAR == mixedEq_eqs_chm(iequa) % kfl_eqtype .or. &
                  &CHM_EQ_YCYC  == mixedEq_eqs_chm(iequa) % kfl_eqtype) then
               kfl_icvar_chm  = iequa
           endif
        enddo
     endif

     !
     ! Maximum column index of lookups
     !
     kfl_max_nvar_lookup_chm = 1_ip
     kfl_max_nvar_ann_in_chm = 1_ip
     kfl_max_nvar_ann_out_chm = 1_ip

     !
     ! Determine which lookups to execute for getting source terms
     !
     kfl_min_srcfw_chm  = max_lookup_fw
     kfl_max_srcfw_chm  = 0_ip
     kfl_min_src_annfw_chm  = max_ann_fw
     kfl_max_src_annfw_chm  = 0_ip
     do iequa = 1,nclas_chm
        !
        ! Set initial fields
        !
        if ( mixedEq_eqs_chm(iequa) % kfl_source_type == CHM_SRC_TABLE ) then
            !
            ! Count number of equations related to each framework
            !
            do jj = 1,3
               ifw = 0_ip
               select case (jj)
                  case (1_ip)
                     if (kfl_multimod_chm == 0_ip ) then
                        ifw = mixedEq_eqs_chm(iequa) % kfl_source_fw
                        !
                        ! Count number of equations related to each framework
                        !
                        kfl_fw_has_sources(ifw) = 1_ip + kfl_fw_has_sources(ifw)
                        !
                        ! Save equation index 
                        !
                        kfl_fw_src_equa_list(ifw,kfl_fw_has_sources(ifw)) = iequa
                     end if 
                  case(2_ip)
                     if (kfl_multimod_chm == 1_ip ) then
                        ifw = mixedEq_eqs_chm(iequa) % kfl_diffsource_fw
                        !
                        ! Count number of equations related to each framework
                        !
                        kfl_fw_has_hybrid_sources(ifw) = 1_ip + kfl_fw_has_hybrid_sources(ifw)
                        !
                        ! Save equation index 
                        !
                        kfl_fw_src_equa_list(ifw,kfl_fw_has_hybrid_sources(ifw)) = iequa
                     end if 
                  case(3_ip)
                     if (kfl_multimod_chm == 1_ip ) then
                        ifw = mixedEq_eqs_chm(iequa) % kfl_premsource_fw
                        !
                        ! Count number of equations related to each framework
                        !
                        kfl_fw_has_hybrid_sources(ifw) = 1_ip + kfl_fw_has_hybrid_sources(ifw)
                        !
                        ! Save equation index 
                        !
                        kfl_fw_src_equa_list(ifw,kfl_fw_has_hybrid_sources(ifw)) = iequa
                     end if 
               end select
               if (ifw /= 0_ip) then 

                  !
                  ! Check lowest and highest used framework index
                  !
                  kfl_min_srcfw_chm = min(ifw,kfl_min_srcfw_chm)
                  kfl_max_srcfw_chm = max(ifw,kfl_max_srcfw_chm)

                  !
                  ! Maximum nvar for equation source terms
                  !
                  kfl_max_nvar_lookup_chm = max(kfl_max_nvar_lookup_chm, lookup_fw(ifw) % main_table % nvar )
               end if
            end do 

        elseif ( mixedEq_eqs_chm(iequa) % kfl_source_type == CHM_SRC_ANN ) then
            !
            ! Count number of equations related to each framework
            !
            ifw = mixedEq_eqs_chm(iequa) % kfl_source_ann
            kfl_annfw_has_sources(ifw) = 1_ip + kfl_annfw_has_sources(ifw)

            !
            ! Save equation index
            !
            kfl_annfw_src_equa_list(ifw,kfl_annfw_has_sources(ifw)) = iequa

            !
            ! Check lowest and highest used framework index
            !
            kfl_min_src_annfw_chm = min(ifw,kfl_min_src_annfw_chm)
            kfl_max_src_annfw_chm = max(ifw,kfl_max_src_annfw_chm)

            !
            ! Maximum nvar for equation source terms
            !
            kfl_max_nvar_ann_in_chm = max(kfl_max_nvar_ann_in_chm, ann_fw(ifw) % scaling_in % n_prod_shape )
            kfl_max_nvar_ann_out_chm = max(kfl_max_nvar_ann_out_chm, ann_fw(ifw) % scaling_out % n_prod_shape )

        endif
     enddo

     !
     ! Initialize association between chemic unknowns and ANNs from kernel for source term calculations
     !
     do ifw = kfl_min_src_annfw_chm, kfl_max_src_annfw_chm
        if (kfl_annfw_has_sources(ifw) > 0) then
            !
            ! This framework is used for source terms
            !
            do ii = 1,ann_fw(ifw) % scaling_in % n_prod_shape
                !
                ! Check for enthalpy or scalar dissipation rate
                !
                if (ann_fw(ifw) % scaling_in % names(ii) == "ENTHA" ) then
                    ann_fw(ifw) % scaling_in % iaux1(ii) = -1_ip
                elseif (ann_fw(ifw) % scaling_in % names(ii) == "CHIZ " ) then
                    ann_fw(ifw) % scaling_in % iaux1(ii) = -2_ip
                elseif (ann_fw(ifw) % scaling_in % names(ii) == "TEMPE" ) then
                    ann_fw(ifw) % scaling_in % iaux1(ii) = -3_ip
                endif


                !
                ! Check if any of the equation names match
                !
                do iequa = 1, nclas_chm
                    if (ann_fw(ifw) % scaling_in % names(ii) == mixedEq_eqs_chm(iequa) % name ) then
                        ann_fw(ifw) % scaling_in % iaux1(ii) = iequa
                    endif
                enddo
            enddo
        endif
     enddo


     !
     ! Determine which columns hold which looked up property
     !
     kfl_W_tab_index_chm            = 0_ip
     kfl_k_tab_index_chm            = 0_ip
     kfl_mu_tab_index_chm           = 0_ip
     kfl_cpCoefLT_tab_index_chm     = 0_ip
     kfl_cpCoefHT_tab_index_chm     = 0_ip
     kfl_cpCoefLT_end_tab_index_chm = 0_ip
     kfl_cpCoefHT_end_tab_index_chm = 0_ip
     kfl_T_tab_index_chm            = 0_ip
     kfl_DtRho_tab_index_chm        = 0_ip
     if (kfl_control_gr_chm > 0 .and. kfl_tab_fw_chm > 0) then
        if (kfl_multimod_chm == 1_ip) then
           do ivari=1,lookup_fw(kfl_tab_fw_chm_diff) % main_table % nvar
              select case(lookup_fw(kfl_tab_fw_chm_diff) % main_table % varname(ivari)) 
              case('W','WSTAR')
                  kfl_W_tab_index_chm       = ivari
              case('K')
                  kfl_k_tab_index_chm       = ivari
              case('MU')
                  kfl_mu_tab_index_chm      = ivari
              case('CP11')
                  kfl_cpCoefLT_tab_index_chm     = ivari
                  kfl_cpCoefLT_end_tab_index_chm = ivari+5
              case('CP21')
                  kfl_cpCoefHT_tab_index_chm     = ivari
                  kfl_cpCoefHT_end_tab_index_chm = ivari+5
              case('T','TSTAR')
                  kfl_T_tab_index_chm       = ivari
              case('DTRHO')
                  kfl_DTRho_tab_index_chm   = ivari
              end select
           enddo
        kfl_max_nvar_lookup_chm = max(kfl_max_nvar_lookup_chm, lookup_fw(kfl_tab_fw_chm_diff) % main_table % nvar )
        else
           do ivari=1,table_fw % main_table % nvar
              select case(table_fw % main_table % varname(ivari)) 
              case('W','WSTAR')
                  kfl_W_tab_index_chm       = ivari
              case('K')
                  kfl_k_tab_index_chm       = ivari
              case('MU')
                  kfl_mu_tab_index_chm      = ivari
              case('CP11')
                  kfl_cpCoefLT_tab_index_chm     = ivari
                  kfl_cpCoefLT_end_tab_index_chm = ivari+5
              case('CP21')
                  kfl_cpCoefHT_tab_index_chm     = ivari
                  kfl_cpCoefHT_end_tab_index_chm = ivari+5
              case('T','TSTAR')
                  kfl_T_tab_index_chm       = ivari
              case('DTRHO')
                  kfl_DTRho_tab_index_chm   = ivari
              end select
           enddo
        !
        ! Maximum nvar for properties
        !
        kfl_max_nvar_lookup_chm = max(kfl_max_nvar_lookup_chm, table_fw % main_table % nvar )
        end if 
     endif



     !
     ! Small workaround for CMC postprocessing
     !
     if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 0) then
         kfl_izmean_chm = 3_ip
         kfl_izvar_chm  = 4_ip
     endif


  case(3_ip)
     !
     ! AFTER READIN AND COMMUNICATING CONSTANTS (parall(1))
     ! BUT BEFORE ALLOCATION OF ARRAYS IN memall(2)
     !

     if(kfl_timei_chm==0) then                      ! Time integration
        dtinv_chm=1.0_rp
     else
        kfl_timei=1
     end if
     momod(modul) % kfl_stead = 0
     kfl_grdif_chm = 0


     !
     ! Dealing with nclas_chm / nspec_chm / etc
     !
     nspec     = nspec_chm
     iclai_chm = 1                                  ! Jacobi starting class
     iclaf_chm = nclas_chm                          ! Jacobi final class


     if( ISLAVE ) then
        do ispec = 1,nspec_chm
           speci(ispec)%name = ''
        enddo
     endif
     !
     ! Time integration definition
     !
     kfl_tiaor_chm = kfl_tiacc_chm                  ! Time accuracy: save original value
     if( kfl_timei_chm == 1 ) then

        if( kfl_tisch_chm == 1 ) then
           ncomp_chm = 3                            ! Trapezoidal rule
        else if( kfl_tisch_chm == 2 ) then
           ncomp_chm = 2+kfl_tiacc_chm              ! BDF scheme
        else if( kfl_tisch_chm == 3 ) then
           ncomp_chm = 2+2                          ! AB scheme
        else if( kfl_tisch_chm == 4 ) then
           ncomp_chm = 2+2                          ! RK scheme
        end if

     else
        ncomp_chm = 2
     end if

     ittot_chm = 0                                  ! Others

     !
     ! To solve equation of state when nastin is not activated
     !
     if ( (kfl_model_chm == 1  .or. kfl_model_chm == 2_ip .or. (kfl_model_chm==4 .and. kfl_solve_cond_CMC_chm == 0))&
         .and. prthe_chm /= 0 ) then
        prthe(1)      = prthe_chm                    ! Low-Mach: p0
        prthe(2)      = prthe_chm                    ! Low-Mach: p0^{n}
        prthe(3)      = prthe_chm                    ! Low-Mach: p0^{n-1}
        prthe(4)      = prthe_chm                    ! Low-Mach: p0^{0}
     end if

     avtim_chm = 0.0_rp                             ! Accumulated time for time-averaging variables
     !
     ! ADR type
     !
     call ADR_initialize_type(ADR_read)
     ADR_read % kfl_time_integration   =  kfl_timei_chm
     ADR_read % kfl_time_step_strategy =  kfl_timco
     ADR_read % kfl_stabilization      =  kfl_stabi_chm
     ADR_read % kfl_shock              =  kfl_shock_chm
     ADR_read % kfl_time_lumped        =  0
     ADR_read % kfl_tau_strategy       =  kfl_taust_chm
     ADR_read % kfl_laplacian          =  0
     ADR_read % kfl_nonlinear_sgs      =  0                ! Deprecated
     ADR_read % kfl_time_sgs           =  0                ! Deprecated
     ADR_read % kfl_time_bubble        =  kfl_tibub_chm
     ADR_read % kfl_time_scheme        =  kfl_tisch_chm
     ADR_read % kfl_time_order         =  kfl_tiacc_chm
     ADR_read % kfl_manufactured       =  0                ! Exact solutions not available in chemic
     ADR_read % kfl_length             =  kfl_ellen_chm
     ADR_read % kfl_first_order_sgs    =  1                ! Related to ADR_chm % kfl_time_order
     ADR_read % number_euler_steps     =  0                ! Deprecated
     ADR_read % lun_output4            =  int(momod(modul) % lun_outpu,4)
     ADR_read % bemol                  =  bemol_chm
     ADR_read % tau_parameters(1:3)    =  staco_chm(1:3)
     ADR_read % shock                  =  shock_chm

     if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1) then
        nclas_max = nvar_CMC_chm
     else
        nclas_max = nclas_chm
     end if

     !
     ! allocate
     !
     allocate(ADR_chm(nclas_max))
     !
     ! extension to ADR_chm(iclas)
     !
     do iclas=1,nclas_max
        call ADR_initialize_type(ADR_chm(iclas))
        ADR_chm(iclas) = ADR_read
        call ADR_check_and_compute_data(ADR_chm(iclas))
     end do

  end select

end subroutine chm_inivar
