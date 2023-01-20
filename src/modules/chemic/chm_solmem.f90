!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_solmem()
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_memall
  ! NAME
  !    chm_memall
  ! DESCRIPTION
  !    This routine allocates memory for the arrays needed
  ! USES
  ! USED BY
  !    chm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_master,                        only : INOTEMPTY, mem_modul, modul, radiative_heat, condu_gp, sphec_gp, visco_gp,&
                                                wmean_gp, tempe_gp, sphec_gp_ht, sphec_gp_lt, enthalpy_transport,&
                                                div_enthalpy_transport, condu_gp, condk, lescl, massk, visck, sphec, sphek, wmean,&
                                                solve
  use def_domain,                        only : ndime, nelem, npoin, ltype, ngaus
  use def_solver,                        only : solve_sol
  use def_kintyp,                        only : ip, rp
  use def_chemic,                        only : mass_gp, massConsumption_gp, DtRho_gp, Yk_ssm_gp, hk_gp, dummy_enthalpy, table_fw,&
                                                posttable_fw, sigma_gp_chm, sigma0_gp_chm, d32_gp_chm, aux_nodes, av_Named_Unk_chm,&
                                                avd32_chm, avden_chm, avL2_chm, avL_chm, avmsk_chm, avS0_chm, avS_chm, avxYr_chm,&
                                                avxYs_chm, avxZr_chm, avxZs_chm, avY2_chm, avY_chm, avYV_chm, avZ2_chm, bvess_chm,&
                                                coeff_cp_k, Corr_chm, d32_chm, dt_chm, dt_rho_chm, React_ind, avposttab_chm,&
                                                avZ_chm, avZV_chm, elem_c, elem_h, elem_n, elem_o, entha_chm,&
                                                enthalpy_transport_nodes, grad_H, grad_k_single, grad_phi_levSet_chm,&
                                                grad_phic_levSet_chm, grad_T, grad_Yk, grad_Yk_hk, hrr_avg_chm, hrr_chm,&
                                                kfl_cpCoefHT_tab_index_chm, kfl_cpCoefLT_tab_index_chm, kfl_DtRho_tab_index_chm,&
                                                kfl_fixno_chm, kfl_k_tab_index_chm, kfl_lookg_chm, kfl_model_chm,&
                                                kfl_mu_tab_index_chm, kfl_solve_cond_CMC_chm, kfl_soot_chm, kfl_spray_chm,&
                                                kfl_tab_fw_chm, kfl_table_ind_chm, kfl_W_tab_index_chm, kfl_yk_fw_ssm,&
                                                lap_phi_levSet_chm, mixfr_chm, nclas_chm, ncomp_chm, nreac_chm, nspec_chm,&
                                                nvar_CMC_chm, nZ_CMC_chm, phi_chm, prog_var_chm, ripts_chm, rspec_chm, Sigm0_chm,&
                                                Sigma_chm, src_chm, sum_reac_chm, sum_Yk_ext_cond_CMC_chm, Text_cond_CMC_chm, W_k,&
                                                xYr_chm, xYs_chm, xZr_chm, xZs_chm, Y_k_n, zgradmax_chm, zgrad_gp, flame_index_gp,&
                                                kfl_tab_fw_chm_diff, kfl_multimod_chm
  use mod_memory,                        only : memory_alloca
  use mod_chm_rk_explicit,               only : chm_rk_explicit_memory
  use mod_chm_operations_CMC,            only : chm_allocate_memory_CMC
  use mod_chm_sectional_soot_model_fast, only : chm_sectional_soot_memory_ssm
  use mod_chm_sectional_soot_model_fast, only : chm_initialization_ssm
  use mod_chm_sectional_soot_model_fast, only : chm_sectional_soot_memory_ssm_fast

  use mod_chm_thermophoretic,            only : chm_thermophoretic_memory
  use mod_output_postprocess,            only : output_postprocess_check_variable_postprocess
  use def_kermod,                        only : lookup_fw, kfl_soot_vect

  implicit none
  integer(ip) :: ielem,pelty,pgaus, ipoin,ireac

  external    :: soldef

  if( INOTEMPTY ) then

     if ((kfl_cpCoefLT_tab_index_chm > 0 .and. kfl_cpCoefHT_tab_index_chm > 0) .or. kfl_model_chm > 2 ) &
         call memory_alloca(mem_modul(1:2,modul),'SPHEC',    'chm_solmem',sphec,npoin,6_ip,2_ip)
     if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVDEN') ) &
         call memory_alloca(mem_modul(1:2,modul),'AVDEN_CHM','chm_solmem',avden_chm,npoin)
     if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVMSK') ) &
         call memory_alloca(mem_modul(1:2,modul),'AVMSK_CHM','chm_solmem',avmsk_chm,npoin)        ! Average mass source from spray
     if (kfl_lookg_chm == 0) &
         call memory_alloca(mem_modul(1:2,modul),'MASSK',    'chm_solmem',massk,npoin,nclas_chm)

     !
     ! Kernel shared variables
     !
     if ((kfl_W_tab_index_chm > 0 .and. kfl_lookg_chm == 0) .or. kfl_model_chm > 2 ) &
         call memory_alloca(mem_modul(1:2,modul),'WMEAN','chm_solmem',wmean,npoin,ncomp_chm)

     !
     ! Projection of dt/rho
     !
     call memory_alloca(mem_modul(1:2,modul),'DT_RHO_CHM','chm_solmem',dt_rho_chm,npoin)
     if (kfl_spray_chm /= 0_ip ) call memory_alloca(mem_modul(1:2,modul),'DT_CHM','chm_solmem',dt_chm,npoin)

     if (kfl_model_chm /= 4) then
        call memory_alloca(mem_modul(1:2,modul),'RADIATIVE_HEAT','chm_solmem',radiative_heat,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'RADIATIVE_HEAT % A','chm_solmem',radiative_heat(ielem)%a,pgaus,1_ip,1_ip)
        end do
     end if

     !
     ! Allocation transport properties
     !
     if (kfl_model_chm /= 4) then

        if (kfl_lookg_chm > 0) then

           if (kfl_k_tab_index_chm > 0 .or. kfl_model_chm > 2 ) then
              call memory_alloca(mem_modul(1:2,modul),'CONDU_GP','chm_solmem',condu_gp,nelem)
              do ielem = 1,nelem
                 pelty = ltype(ielem)
                 pgaus = ngaus(pelty)
                 call memory_alloca(mem_modul(1:2,modul),'CONDU_GP % A','chm_solmem',condu_gp(ielem)%a,pgaus,1_ip,1_ip)
              end do
           endif

           if (kfl_cpCoefLT_tab_index_chm > 0 .or. kfl_cpCoefHT_tab_index_chm > 0 .or. kfl_model_chm > 2 ) then
              call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP','chm_solmem',sphec_gp,nelem)
              do ielem = 1,nelem
                 pelty = ltype(ielem)
                 pgaus = ngaus(pelty)
                 call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP % A','chm_solmem',sphec_gp(ielem)%a,pgaus,1_ip,1_ip)
                 sphec_gp(ielem)%a=1.0_rp
              end do
           endif

           if (kfl_mu_tab_index_chm > 0 .or. kfl_model_chm > 2 ) then
              call memory_alloca(mem_modul(1:2,modul),'VISCO_GP','chm_solmem',visco_gp,nelem)
              do ielem = 1,nelem
                 pelty = ltype(ielem)
                 pgaus = ngaus(pelty)
                 call memory_alloca(mem_modul(1:2,modul),'VISCO_GP % A','chm_solmem',visco_gp(ielem)%a,pgaus,1_ip,1_ip)
              end do
           endif

           if (kfl_W_tab_index_chm > 0 .or. kfl_model_chm > 2 ) then
              call memory_alloca(mem_modul(1:2,modul),'WMEAN_GP','chm_solmem',wmean_gp,nelem)
              do ielem = 1,nelem
                 pelty = ltype(ielem)
                 pgaus = ngaus(pelty)
                 call memory_alloca(mem_modul(1:2,modul),'WMEAN_GP % A','chm_solmem',wmean_gp(ielem)%a,pgaus,4_ip,1_ip)
              end do
           endif

           call memory_alloca(mem_modul(1:2,modul),'TEMPE_GP','chm_solmem',tempe_gp,nelem)
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              call memory_alloca(mem_modul(1:2,modul),'TEMPE_GP % A','chm_solmem',tempe_gp(ielem)%a,pgaus,1_ip,1_ip)
              tempe_gp(ielem)%a=200.0_rp
           end do

           if (kfl_cpCoefHT_tab_index_chm > 0 .or. kfl_model_chm > 2 ) then
              call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP_HT','chm_solmem',sphec_gp_ht,nelem)
              do ielem = 1,nelem
                 pelty = ltype(ielem)
                 pgaus = ngaus(pelty)
                 call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP_HT % A','chm_solmem',sphec_gp_ht(ielem)%a,pgaus,6_ip,1_ip)
                 sphec_gp_ht(ielem)%a=1.0_rp
              end do
           endif

           if (kfl_cpCoefLT_tab_index_chm > 0 .or. kfl_model_chm > 2 ) then
              call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP_LT','chm_solmem',sphec_gp_lt,nelem)
              do ielem = 1,nelem
                 pelty = ltype(ielem)
                 pgaus = ngaus(pelty)
                 call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP_LT % A','chm_solmem',sphec_gp_lt(ielem)%a,pgaus,6_ip,1_ip)
                 sphec_gp_lt(ielem)%a=1.0_rp
              end do
           endif

           call memory_alloca(mem_modul(1:2,modul),'MASS_GP','chm_solmem',mass_gp,nelem)
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              call memory_alloca(mem_modul(1:2,modul),'MASS_GP % A','chm_solmem',mass_gp(ielem)%a,pgaus,nclas_chm,1_ip)
           end do

           call memory_alloca(mem_modul(1:2,modul),'ZGRAD_GP','chm_solmem',zgrad_gp,nelem)
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              call memory_alloca(mem_modul(1:2,modul),'ZGRAD_GP % A','chm_solmem',zgrad_gp(ielem)%a,pgaus,1_ip)
           end do

           call memory_alloca(mem_modul(1:2,modul),'FLAME_INDEX_GP','chm_solmem',flame_index_gp,nelem)
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              call memory_alloca(mem_modul(1:2,modul),'FLAME_INDEX_GP % A','chm_solmem',flame_index_gp(ielem)%a,pgaus,1_ip)
           end do

           call memory_alloca(mem_modul(1:2,modul),'MASSCONSUMPTION_GP','chm_solmem',massConsumption_gp,nelem)
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              call memory_alloca(mem_modul(1:2,modul),'MASSCONSUMPTION_GP % A','chm_solmem',massConsumption_gp(ielem)%a,pgaus,&
                  nclas_chm,1_ip)
           end do

           if (kfl_DtRho_tab_index_chm > 0 ) then
              call memory_alloca(mem_modul(1:2,modul),'DTRHO_GP','chm_solmem',DtRho_gp,nelem)
              do ielem = 1,nelem
                 pelty = ltype(ielem)
                 pgaus = ngaus(pelty)
                 call memory_alloca(mem_modul(1:2,modul),'DTRHO_GP % A','chm_solmem',DtRho_gp(ielem)%a,pgaus,1_ip,1_ip)
              end do
           endif

           if ( kfl_yk_fw_ssm > 0 ) then
              call memory_alloca(mem_modul(1:2,modul),'YK_SSM_GP','chm_solmem',Yk_ssm_gp,nelem)
              do ielem = 1,nelem
                 pelty = ltype(ielem)
                 pgaus = ngaus(pelty)
                 call memory_alloca(mem_modul(1:2,modul),'YK_SSM_GP','chm_solmem',Yk_ssm_gp(ielem)%a,pgaus,&
                     lookup_fw(kfl_yk_fw_ssm)%main_table%nvar,1_ip)
              end do
           end if


           if ( kfl_model_chm == 3 ) then
              call memory_alloca(mem_modul(1:2,modul),'HK_GP','chm_solmem',hk_gp,nelem)
              do ielem = 1,nelem
                 pelty = ltype(ielem)
                 pgaus = ngaus(pelty)
                 call memory_alloca(mem_modul(1:2,modul),'HK_GP % A','chm_solmem',hk_gp(ielem)%a,nclas_chm,pgaus,1_ip)
                 hk_gp(ielem)%a=0.0_rp
              end do
           end if

        else
           call memory_alloca(mem_modul(1:2,modul),'VISCK','chm_solmem',visck,npoin,nclas_chm)
           call memory_alloca(mem_modul(1:2,modul),'CONDK','chm_solmem',condk,npoin,nclas_chm)
           call memory_alloca(mem_modul(1:2,modul),'SPHEK','chm_solmem',sphek,npoin,nclas_chm)
        endif

     end if


     !
     ! Finite rate chemistry model variables
     !
     if ( kfl_model_chm == 3 ) then
        call memory_alloca(mem_modul(1:2,modul),'ENTHALPY_TRANSPORT','chm_solmem',enthalpy_transport,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'ENTHALPY_TRANSPORT % A','chm_solmem',enthalpy_transport(ielem)%a,pgaus,&
               ndime,1_ip)
           enthalpy_transport(ielem)%a=0.0_rp
        end do

        call memory_alloca(mem_modul(1:2,modul),'DIV_ENTHALPY_TRANSPORT','chm_solmem',div_enthalpy_transport,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'DIV_ENTHALPY_TRANSPORT % A','chm_solmem',div_enthalpy_transport(ielem)%a,pgaus,&
               1_ip,1_ip)
           div_enthalpy_transport(ielem)%a=0.0_rp
        end do

        call memory_alloca(mem_modul(1:2,modul),'dummy_enthalpy','chm_solmem',dummy_enthalpy,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'dummy_enthalpy % A','chm_solmem',dummy_enthalpy(ielem)%a,pgaus,1_ip,1_ip)
           dummy_enthalpy(ielem)%a=0.0_rp
        end do

        call memory_alloca(mem_modul(1:2,modul),'GRAD_YK',                 'chm_solmem',grad_Yk,                 nspec_chm,ndime,&
            npoin)
        call memory_alloca(mem_modul(1:2,modul),'GRAD_YK_HK',              'chm_solmem',grad_Yk_hk,              nspec_chm,ndime,&
            npoin)
        call memory_alloca(mem_modul(1:2,modul),'GRAD_T',                  'chm_solmem',grad_T,                  ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'AUX_NODES',               'chm_solmem',aux_nodes,               npoin)
        call memory_alloca(mem_modul(1:2,modul),'GRAD_H',                  'chm_solmem',grad_H,                  ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'GRAD_K_SINGLE',           'chm_solmem',grad_k_single,           ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'ENTHALPY_TRANSPORT_NODES','chm_solmem',enthalpy_transport_nodes,ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'Y_K_N',                   'chm_solmem',Y_k_n,                   npoin,nspec_chm,&
            ncomp_chm)
        call memory_alloca(mem_modul(1:2,modul),'REACT_IND',               'chm_solmem',React_ind,               npoin,nreac_chm)
        call memory_alloca(mem_modul(1:2,modul),'CORR_CHM',                'chm_solmem',Corr_chm,                npoin)
        call memory_alloca(mem_modul(1:2,modul),'SRC_CHM',                 'chm_solmem',src_chm,                 npoin,nclas_chm)
        call memory_alloca(mem_modul(1:2,modul),'ELEM_C',                  'chm_solmem',elem_c,                  npoin)
        call memory_alloca(mem_modul(1:2,modul),'ELEM_H',                  'chm_solmem',elem_h,                  npoin)
        call memory_alloca(mem_modul(1:2,modul),'ELEM_N',                  'chm_solmem',elem_n,                  npoin)
        call memory_alloca(mem_modul(1:2,modul),'ELEM_O',                  'chm_solmem',elem_o,                  npoin)
        call memory_alloca(mem_modul(1:2,modul),'MIXFR_CHM',               'chm_solmem',mixfr_chm,               npoin)
        call memory_alloca(mem_modul(1:2,modul),'PROG_VAR_CHM',            'chm_solmem',prog_var_chm,            npoin)
        call memory_alloca(mem_modul(1:2,modul),'SUM_REAC_CHM',            'chm_solmem',sum_reac_chm,            npoin)
        call memory_alloca(mem_modul(1:2,modul),'HRR_CHM',                 'chm_solmem',hrr_chm,                 npoin)
        call memory_alloca(mem_modul(1:2,modul),'HRR_AVG_CHM',             'chm_solmem',hrr_avg_chm,             npoin)
        call memory_alloca(mem_modul(1:2,modul),'ENTHA_CHM',               'chm_solmem',entha_chm,               nclas_chm,npoin)

        do ipoin = 1,npoin
           do ireac = 1,nreac_chm
              React_ind(ipoin,ireac)       = 1_ip
           enddo
           Corr_chm(ipoin)        = 1_ip
        enddo


     elseif( kfl_model_chm == 2 .or. kfl_model_chm == 1) then
        !
        ! Species fields for radiation for flamelet model
        !
        call memory_alloca(mem_modul(1:2,modul),'RSPEC_CHM'   ,'chm_solmem',rspec_chm,2_ip,npoin)
        call memory_alloca(mem_modul(1:2,modul),'ZGRADMAX_CHM','chm_solmem',zgradmax_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'PHI_CHM'     ,'chm_solmem',phi_chm,npoin)

        call memory_alloca(mem_modul(1:2,modul),'LESCL','chm_solmem',lescl,npoin)

        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVY  ') ) &
            call memory_alloca(mem_modul(1:2,modul),'AVY_CHM'  ,'chm_solmem',avY_chm, npoin)
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVZ  ') ) &
            call memory_alloca(mem_modul(1:2,modul),'AVZ_CHM'  ,'chm_solmem',avZ_chm, npoin)
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVYV ') ) &
            call memory_alloca(mem_modul(1:2,modul),'AVYV_CHM' ,'chm_solmem',avYV_chm,npoin)
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVZV ') ) &
            call memory_alloca(mem_modul(1:2,modul),'AVZV_CHM' ,'chm_solmem',avZV_chm,npoin)
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVZ2 ') ) &
            call memory_alloca(mem_modul(1:2,modul),'AVZ2_CHM' ,'chm_solmem',avZ2_chm,npoin)
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVY2 ') ) &
            call memory_alloca(mem_modul(1:2,modul),'AVY2_CHM' ,'chm_solmem',avY2_chm,npoin)
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVL  ') ) &
            call memory_alloca(mem_modul(1:2,modul),'AVL_CHM'  ,'chm_solmem',avL_chm, npoin)
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVL2 ') ) &
            call memory_alloca(mem_modul(1:2,modul),'AVL2_CHM' ,'chm_solmem',avL2_chm,npoin)

        !
        ! Allocation scalar dissipation rates
        !
        call memory_alloca(mem_modul(1:2,modul),'XYR_CHM','chm_solmem',xYr_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'XZR_CHM','chm_solmem',xZr_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'XYS_CHM','chm_solmem',xYs_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'XZS_CHM','chm_solmem',xZs_chm,npoin)

        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVXYR') ) &
            call memory_alloca(mem_modul(1:2,modul),'AVXYR_CHM',    'chm_solmem',avxYr_chm,     npoin)
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVXZR') ) &
            call memory_alloca(mem_modul(1:2,modul),'AVXZR_CHM',    'chm_solmem',avxZr_chm,     npoin)
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVXYS') ) &
            call memory_alloca(mem_modul(1:2,modul),'AVXYS_CHM',    'chm_solmem',avxYs_chm,     npoin)
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVXZS') ) &
            call memory_alloca(mem_modul(1:2,modul),'AVXZS_CHM',    'chm_solmem',avxZs_chm,     npoin)

        !
        ! Allocation averaged named unknowns mixed equation model
        !
        nullify(av_Named_Unk_chm)
        call memory_alloca(mem_modul(1:2,modul),'AV_NAMED_UNK_CHM','chm_solmem',av_Named_Unk_chm,npoin,nclas_chm)
        av_Named_Unk_chm = 0.0_rp
        !
        ! Allocate variable to save initial guess for lookup indeces
        !
        if (kfl_tab_fw_chm > 0) then
           if (kfl_multimod_chm == 1_ip) then
              call memory_alloca(mem_modul(1:2,modul),'KFL_TABLE_IND_CHM',    'chm_solmem',kfl_table_ind_chm,&
                   lookup_fw(kfl_tab_fw_chm_diff) % main_table % ndim, nelem)
           else
              call memory_alloca(mem_modul(1:2,modul),'KFL_TABLE_IND_CHM',    'chm_solmem',kfl_table_ind_chm,&
                  table_fw % main_table % ndim, nelem)
           end if 
           kfl_table_ind_chm = 1_ip
        endif


        if (associated(posttable_fw)) then
           call memory_alloca(mem_modul(1:2,modul),'AVPOSTTAB_CHM','chm_solmem',avposttab_chm, npoin,&
               posttable_fw % main_table % nvar)
        endif

     endif


     ! Allocate variables for CMC model
     if ( kfl_model_chm == 4 ) then
        call chm_allocate_memory_CMC
     end if


     !
     ! Allocate spray terms
     !
     if ( kfl_spray_chm /= 0 ) then

        call memory_alloca(mem_modul(1:2,modul),'AVS_CHM'  ,'chm_solmem',avS_chm  ,npoin)
        call memory_alloca(mem_modul(1:2,modul),'AVS0_CHM' ,'chm_solmem',avS0_chm ,npoin)
        call memory_alloca(mem_modul(1:2,modul),'AVD32_CHM','chm_solmem',avd32_chm,npoin)

        call memory_alloca(mem_modul(1:2,modul),'SIGMA_CHM','chm_solmem',Sigma_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'SIGM0_CHM','chm_solmem',Sigm0_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'D32_CHM'  ,'chm_solmem',d32_chm,  npoin)

        call memory_alloca(mem_modul(1:2,modul),'SIGMA_GP_CHM' ,'chm_solmem',sigma_gp_chm ,nelem)
        call memory_alloca(mem_modul(1:2,modul),'SIGMA0_GP_CHM','chm_solmem',sigma0_gp_chm,nelem)
        call memory_alloca(mem_modul(1:2,modul),'D32_GP_CHM'   ,'chm_solmem',d32_gp_chm   ,nelem)

        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'SIGMA_GP_CHM % A' ,'chm_solmem',sigma_gp_chm(ielem)%a ,pgaus,1_ip,1_ip)
           call memory_alloca(mem_modul(1:2,modul),'SIGMA0_GP_CHM % A','chm_solmem',sigma0_gp_chm(ielem)%a,pgaus,1_ip,1_ip)
           call memory_alloca(mem_modul(1:2,modul),'D32_GP_CHM % A'   ,'chm_solmem',d32_gp_chm(ielem)%a   ,pgaus,1_ip,1_ip)
        end do

     end if
     !
     ! Allocate variables for level set
     !
     if ( kfl_spray_chm == 2 ) then
        call memory_alloca(mem_modul(1:2,modul),'LAP_PHI_LEVSET_CHM' ,'chm_solmem',lap_phi_levSet_chm ,npoin)
        call memory_alloca(mem_modul(1:2,modul),'GRAD_PHI_LEVSET_CHM','chm_solmem',grad_phi_levSet_chm,ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'GRAD_PHIC_LEVSET_CHM','chm_solmem',grad_phic_levSet_chm,ndime,npoin)
     end if

  else
     !
     ! Projection of dt/rho
     !
     call memory_alloca(mem_modul(1:2,modul),'DT_RHO_CHM','chm_solmem',dt_rho_chm,1_ip)
     if (kfl_spray_chm /= 0_ip ) call memory_alloca(mem_modul(1:2,modul),'DT_CHM','chm_solmem',dt_chm,1_ip)

     !
     ! Allocation scalar dissipation rates
     !
     call memory_alloca(mem_modul(1:2,modul),'XYR_CHM','chm_solmem',xYr_chm,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'XZR_CHM','chm_solmem',xZr_chm,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'XYS_CHM','chm_solmem',xYs_chm,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'XZS_CHM','chm_solmem',xZs_chm,1_ip)

     call memory_alloca(mem_modul(1:2,modul),'AVXYR_CHM','chm_solmem',avxYr_chm,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'AVXZR_CHM','chm_solmem',avxZr_chm,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'AVXYS_CHM','chm_solmem',avxYs_chm,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'AVXZS_CHM','chm_solmem',avxZs_chm,1_ip)

     nullify(av_Named_Unk_chm)
     call memory_alloca(mem_modul(1:2,modul),'AV_NAMED_UNK_CHM','chm_solmem',av_Named_Unk_chm,1_ip,nclas_chm)
     av_Named_Unk_chm = 0.0_rp

     if (kfl_lookg_chm > 0) then
        call memory_alloca(mem_modul(1:2,modul),'CONDU_GP','chm_solmem',condu_gp,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP','chm_solmem',sphec_gp,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'VISCO_GP','chm_solmem',visco_gp,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'WMEAN_GP','chm_solmem',wmean_gp,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'TEMPE_GP','chm_solmem',tempe_gp,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP_HT','chm_solmem',sphec_gp_ht,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP_LT','chm_solmem',sphec_gp_lt,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'MASS_GP','chm_solmem',mass_gp,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'YK_SSM_GP','chm_solmem',Yk_ssm_gp,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'MASSCONSUMPTION_GP','chm_solmem',massConsumption_gp,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'ZGRAD_GP','chm_solmem',zgrad_gp,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'FLAME_INDEX_GP','chm_solmem',flame_index_gp,1_ip)

        call memory_alloca(mem_modul(1:2,modul),'ZGRAD_GP','chm_solmem',zgrad_gp(1)%a,1_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'FLAME_INDEX_GP','chm_solmem',flame_index_gp(1)%a,1_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'CONDU_GP','chm_solmem',condu_gp(1)%a,1_ip,1_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP','chm_solmem',sphec_gp(1)%a,1_ip,1_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'VISCO_GP','chm_solmem',visco_gp(1)%a,1_ip,1_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'WMEAN_GP','chm_solmem',wmean_gp(1)%a,1_ip,4_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'TEMPE_GP','chm_solmem',tempe_gp(1)%a,1_ip,4_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP_HT','chm_solmem',sphec_gp_ht(1)%a,1_ip,6_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP_LT','chm_solmem',sphec_gp_lt(1)%a,1_ip,6_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'MASS_GP','chm_solmem',mass_gp(1)%a,1_ip,nclas_chm,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'YK_SSM_GP','chm_solmem',Yk_ssm_gp(1)%a,1_ip,15_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'MASSCONSUMPTION_GP','chm_solmem',massConsumption_gp(1)%a,1_ip,nclas_chm,1_ip)
     endif

     !
     ! Allocate spray terms
     !
     if ( kfl_spray_chm /= 0 ) then

        call memory_alloca(mem_modul(1:2,modul),'AVS_CHM'  ,'chm_solmem',avS_chm  ,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'AVS0_CHM' ,'chm_solmem',avS0_chm ,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'AVD32_CHM','chm_solmem',avd32_chm,1_ip)

        call memory_alloca(mem_modul(1:2,modul),'SIGMA_GP_CHM' ,'chm_solmem',sigma_gp_chm ,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'SIGMA0_GP_CHM','chm_solmem',sigma0_gp_chm,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'D32_GP_CHM'   ,'chm_solmem',d32_gp_chm   ,1_ip)

     end if

     !
     ! Allocate variables for level set
     !
     if ( kfl_spray_chm == 2 ) then
        call memory_alloca(mem_modul(1:2,modul),'LAP_PHI_LEVSET_CHM' ,'chm_solmem',lap_phi_levSet_chm   ,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'GRAD_PHI_LEVSET_CHM','chm_solmem',grad_phi_levSet_chm  ,1_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'GRAD_PHIC_LEVSET_CHM','chm_solmem',grad_phic_levSet_chm,1_ip,1_ip)
     end if


     !
     ! Allocate variables for CMC
     !
     if ( kfl_model_chm == 4_ip .and. kfl_solve_cond_CMC_chm == 1_ip ) then
       call memory_alloca(mem_modul(1:2,modul),'EXTREME_COND_TEMP_CMC','chm_solmem',Text_cond_CMC_chm,2_ip*nZ_CMC_chm)
       call memory_alloca(mem_modul(1:2,modul),'EXTREME_COND_SUM_YK_CMC','chm_solmem',sum_Yk_ext_cond_CMC_chm,2_ip*nZ_CMC_chm)
       sum_Yk_ext_cond_CMC_chm = 1.0_rp
     end if
  end if

  if ( kfl_model_chm == 3 ) then
     call memory_alloca(mem_modul(1:2,modul),'COEFF_CP_K','chm_solmem',coeff_cp_k,nspec_chm,15_ip)
     call memory_alloca(mem_modul(1:2,modul),'W_K','chm_solmem',W_k,nspec_chm)
  endif

  !
  ! Allocate variables for soot model
  !
  if ( kfl_soot_chm /= 0 ) then
     if( kfl_soot_vect /= 0 ) then
        call chm_sectional_soot_memory_ssm_fast()
     else
        call chm_sectional_soot_memory_ssm()
     end if
     !!call chm_initialization_ssm()
  endif

  call chm_thermophoretic_memory(nclas_chm)

  !
  ! Memory (Solver)
  !
  if ( kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1 ) then
     solve(1) % ndofn = nvar_CMC_chm
     solve(2) % ndofn = nvar_CMC_chm
  else
     solve(1) % ndofn = nclas_chm
     solve(2) % ndofn = nclas_chm
  end if

  solve_sol => solve(1:)
  call soldef(4_ip)
  !
  ! Residual
  !
  if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1 ) then
     call memory_alloca(mem_modul(1:2,modul),'RIPTS_CHM','chm_solmem',ripts_chm,nvar_CMC_chm*nclas_chm)
  else
     call memory_alloca(mem_modul(1:2,modul),'RIPTS_CHM','chm_solmem',ripts_chm,nclas_chm*nclas_chm)
  end if

  !
  ! Boundary conditions
  !
  solve(1) % bvess     => bvess_chm
  solve(1) % kfl_fixno => kfl_fixno_chm
  solve(2) % bvess     => bvess_chm
  solve(2) % kfl_fixno => kfl_fixno_chm
  ! Set Dirichlet B.C. imposed by solver
  solve(1) % kfl_iffix =  1

  !
  ! Memory
  !
  call chm_rk_explicit_memory()

end subroutine chm_solmem

