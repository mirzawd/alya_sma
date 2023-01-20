!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_parall(order)
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_parall
  ! NAME
  !    chm_sendat
  ! DESCRIPTION
  !    This routine exchanges data
  ! USES
  ! USED BY
  !    chm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_master,         only : ID_CHEMIC, ID_NASTIN, INOTMASTER, ISEQUEN, ISLAVE, speci, kfl_htran,&
       prthe, kfl_coupl
  use def_chemic,         only : nsize_mech_name, nsize_red, Field_ind_chm, write_uncond_spec_CMC_chm, write_cond_spec_CMC_chm,&
       write_iZ_spec_CMC_chm, bemol_chm, bf_fuel_chm, bo_oxy_chm, chemical_time_factor, cotol_chm,&
       cutof_chm, dac_cor_chm, dac_crit_chm, dampi_chm, droplet_compactness_limit_chm,&
       droplet_h_factor_chm, droplet_max_diameter_chm, droplet_postprocess_frequency_chm, dtmax_chm,&
       dtmin_chm, epsht_chm, epstr_chm, index_N2, inv_exp_alpha_grid_CMC_chm, kfl_advec_chm,&
       kfl_allcl_chm, kfl_avg_cond_CMC_chm, kfl_bc_alpha_CMC_chm, kfl_bc_init_method_CMC_chm,&
       kfl_control_gr_chm, kfl_diffu_chm, kfl_disable_entpre_chm, kfl_droplet_id_chm,&
       kfl_dt_calc_CMC_chm, kfl_dtcri_chm, kfl_ellen_chm, kfl_entropy_chm, kfl_fields_scale_chm,&
       kfl_freq_chm, kfl_hrr_col_chm, kfl_hrr_fw_chm, kfl_incl_PDF_trans_CMC_chm, kfl_int_chm,&
       kfl_key_chm, kfl_limit_chm, kfl_lookg_chm, kfl_mesh_interp_CMC_chm, kfl_model_chm, kfl_negat_chm,&
       kfl_norma_chm, kfl_normc_chm, kfl_pfa_chm, kfl_posit_chm, kfl_post_fw_chm, kfl_post_gp_CMC_chm,&
       kfl_premix_chm, kfl_radia_chm, kfl_shock_chm, kfl_solve_cond_CMC_chm, kfl_solve_enth_CMC_chm,&
       kfl_soot_chm, kfl_spec_name_chm, kfl_split_CFD_CMC, kfl_split_chm, kfl_spray_chm, kfl_stabi_chm,&
       kfl_tab_fw_chm, kfl_taust_chm, kfl_tdac_write_chm, kfl_temli_chm, kfl_tiacc_chm, kfl_tibub_chm,&
       kfl_timei_chm, kfl_tisch_chm, kfl_trans_mxt_spc_CMC_chm, kfl_trans_phs_spc_CMC_chm,&
       kfl_transport_chm, kfl_ufpv_chm, kfl_varYc_chm, kfl_varZ_chm, kfl_wallc_chm, kfl_warni_chm,&
       kfl_weigh_in_eq_CMC_chm, kfl_yk_fw_ssm, kfl_z_chm, levelSet_threshold_chm, mechanism_path,&
       mixedEq_eqs_chm, mixedEq_groups_chm, nclas_chm, ngrou_chm, npara_chm, nreac_chm, nS_AMC_CMC_chm,&
       nsize_mech_name, nsize_red, nspec_chm, nspec_cond_write_CMC_chm, nspec_uncond_write_CMC_chm,&
       nvar_CMC_chm, nvar_therm_CMC_chm, nZ_AMC_CMC_chm, nZ_chm_int_CMC_chm, nZ_CMC_chm,&
       nZ_no_chm_int_CMC_chm, nZ_write_CMC_chm, order_RK_CMC_chm, pos_Zstq_CMC_chm, posttable_fw,&
       prthe_chm, radwt_chm, Red_spec, relax_chm, S_threshold, safet_chm, shock_chm, size_tncod_CMC_chm,&
       Smax_AMC_CMC_chm, sorad_chm, sstol_chm, strec_chm, surf_tension_chm, table_coords, table_fw,&
       table_tab, tbcod_chm, temli_chm, tncod_chm, write_iZ_spec_CMC_chm, Xtot_const_CMC_chm,&
       Zavg_const_CMC_chm, Zs_CMC_chm, Zstq_CMC_chm, Zvar_const_CMC_chm, kfl_field_chm, socen_chm,&
       stofu_chm, staco_chm, posZ_chem_integr_CMC_chm, extremes_Z_CMC_chm, ipara_chm,&
       rpara_chm, Z_CMC_chm, Xnormalized_prof_CMC_chm, PDF_min_CMC_chm, diff_Z_CMC_chm,&
       temp_inert_CMC_chm, temp_equil_CMC_chm, rscal_inert_CMC_chm, rscal_equil_CMC_chm, Z_AMC_CMC_chm,&
       Xintegrated_table_AMC_CMC_chm, S_AMC_CMC_chm, Zlim_clipZvar_CMC_chm, Slim_clipZvar_CMC_chm,&
       max_diff_Yk_inert_eq_CMC_chm, chem_int_iZ_CMC_chm, diffu_chm, Le_k, kfl_initi_chm, xinit_chm,&
       kfl_zg_fw_chm,kfl_zg_col_chm,kfl_tab_fw_chm_diff,kfl_tab_fw_chm_prem,chm_zmax,chm_zmin,&
       kfl_multimod_chm 
  use mod_memchk,         only : memchk
  use def_kintyp,         only : ip
  use mod_opebcs,         only : boundary_conditions_exchange
  use def_kermod,         only : gasco
  use def_kermod,         only : lookup_fw
  use mod_interp_tab,     only : tab_par_exchange
  use mod_chm_mixedEq,    only : chm_mixedEq_par_exchange
  use mod_exchange,           only : exchange_init
  use mod_exchange,           only : exchange_add
  use mod_exchange,           only : exchange_end
  use mod_solver,             only : solver_parall
  use mod_output_postprocess, only : output_postprocess_parall
#ifdef CANTERA
  use cantera,            only : importPhase
  use def_chemic,         only : gas_chm
#endif
  use mod_chm_sectional_soot_model_fast, only: nsect_ssm,nclas_ssm,nspec_ssm
  use mod_chm_sectional_soot_model_fast, only: RhoC_ssm,Vmax_ssm,RadSoot_ssm
  use mod_chm_sectional_soot_model_fast, only: gasCoupling_ssm, &
       ID_NUCL,ID_COND,ID_COAG,ID_SURF
  use mod_chm_operations_CMC,            only : chm_find_master_CMC, &
       chm_transfer_info_CMC_to_CFD, &
       chm_transfer_info_CFD_to_CMC, &
       chm_write_info_transfer_CMC
  use mod_output_postprocess,     only : output_postprocess_parall_old

  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: ki,ki2,iclas
  integer(ip)             :: r_master

  external                :: memerr
  external                :: chm_memphy
  external                :: chm_memnut
  external                :: chm_memous

  if( ISEQUEN ) return

  select case (order)

  case(1_ip)
     !
     ! Exchange data read in chm_reaphy, chm_reanut and chm_reaous
     !
     call exchange_init()

     !----------------------------------------------------------------
     !
     ! Exchange of chm_reaphy variables
     !
     !----------------------------------------------------------------

     call exchange_add(kfl_model_chm)
     call exchange_add(kfl_timei_chm)
     call exchange_add(kfl_advec_chm)
     call exchange_add(kfl_diffu_chm)
     call exchange_add(kfl_transport_chm)
     call exchange_add(kfl_norma_chm)
     call exchange_add(kfl_radia_chm)
     call exchange_add(kfl_pfa_chm)
     call exchange_add(kfl_z_chm)
     call exchange_add(kfl_key_chm)
     call exchange_add(kfl_freq_chm)
     call exchange_add(kfl_premix_chm)
     call exchange_add(kfl_varYc_chm)
     call exchange_add(kfl_varZ_chm)
     call exchange_add(kfl_control_gr_chm)
     call exchange_add(kfl_ufpv_chm)
     call exchange_add(kfl_tdac_write_chm)
     call exchange_add(kfl_lookg_chm)
     call exchange_add(kfl_tab_fw_chm)
     call exchange_add(kfl_post_fw_chm)
     call exchange_add(kfl_hrr_fw_chm)
     call exchange_add(kfl_hrr_col_chm)
     call exchange_add(kfl_zg_fw_chm)
     call exchange_add(kfl_zg_col_chm)
     call exchange_add(kfl_entropy_chm)
     call exchange_add(kfl_disable_entpre_chm)
     call exchange_add(kfl_spec_name_chm)
     call exchange_add(kfl_htran)
     call exchange_add(kfl_spray_chm)
     call exchange_add(kfl_field_chm(1))
     call exchange_add(kfl_field_chm(2))
     call exchange_add(kfl_field_chm(3))
     call exchange_add(kfl_field_chm(4))
     call exchange_add(kfl_tab_fw_chm_diff)
     call exchange_add(kfl_tab_fw_chm_prem)
     call exchange_add(kfl_multimod_chm)
     call exchange_add(chm_zmax)
     call exchange_add(chm_zmin)
     !
     ! Soot model
     !
     call exchange_add(kfl_soot_chm)

     call exchange_add(ID_NUCL)
     call exchange_add(ID_COND)
     call exchange_add(ID_COAG)
     call exchange_add(ID_SURF)
     call exchange_add(nsect_ssm)
     call exchange_add(gasCoupling_ssm)

     call exchange_add(RhoC_ssm)
     call exchange_add(Vmax_ssm)
     call exchange_add(RadSoot_ssm)

     call exchange_add(kfl_yk_fw_ssm)

     !
     ! End sectional soot model variables
     !

     call exchange_add(nclas_chm)
     call exchange_add(ngrou_chm)
     call exchange_add(nspec_chm)
     call exchange_add(nspec_ssm)
     call exchange_add(nclas_ssm)
     call exchange_add(nsect_ssm)
     call exchange_add(nsize_mech_name)
     call exchange_add(nsize_red)
     call exchange_add(prthe_chm)
     if (kfl_coupl(ID_CHEMIC,ID_NASTIN) == 0 ) then
        call exchange_add(prthe(1))
     end if
     call exchange_add(radwt_chm)
     call exchange_add(dac_crit_chm)
     call exchange_add(dac_cor_chm)
     call exchange_add(sorad_chm)
     call exchange_add(socen_chm(1))
     call exchange_add(socen_chm(2))
     call exchange_add(socen_chm(3))
     call exchange_add(surf_tension_chm)
     call exchange_add(bf_fuel_chm)
     call exchange_add(bo_oxy_chm)

     call exchange_add(nreac_chm)

     call exchange_add(kfl_droplet_id_chm)
     call exchange_add(droplet_postprocess_frequency_chm)
     call exchange_add(levelSet_threshold_chm)
     call exchange_add(droplet_compactness_limit_chm)
     call exchange_add(droplet_max_diameter_chm)
     call exchange_add(droplet_h_factor_chm)

     call exchange_add(stofu_chm(1))
     call exchange_add(stofu_chm(2))
     call exchange_add(stofu_chm(3))

     call exchange_add(gasco)

     call exchange_add(kfl_ellen_chm)
     call exchange_add(kfl_taust_chm)
     call exchange_add(kfl_shock_chm)
     call exchange_add(kfl_stabi_chm)
     call exchange_add(kfl_limit_chm)
     call exchange_add(kfl_tiacc_chm)
     call exchange_add(kfl_split_chm)
     call exchange_add(kfl_int_chm)
     call exchange_add(kfl_wallc_chm)
     call exchange_add(kfl_tibub_chm)
     call exchange_add(kfl_tisch_chm)
     call exchange_add(kfl_normc_chm)
     call exchange_add(kfl_dtcri_chm)
     call exchange_add(kfl_negat_chm)
     call exchange_add(kfl_posit_chm)
     call exchange_add(kfl_warni_chm)
     call exchange_add(kfl_temli_chm)

     call exchange_add(staco_chm(1))
     call exchange_add(staco_chm(2))
     call exchange_add(staco_chm(3))
     call exchange_add(shock_chm)
     call exchange_add(bemol_chm)
     call exchange_add(temli_chm)
     call exchange_add(cotol_chm)
     call exchange_add(safet_chm)
     call exchange_add(chemical_time_factor)
     call exchange_add(cutof_chm)
     call exchange_add(sstol_chm)
     call exchange_add(strec_chm)
     call exchange_add(dampi_chm)
     call exchange_add(epsht_chm)
     call exchange_add(epstr_chm)
     call exchange_add(dtmin_chm)
     call exchange_add(dtmax_chm)
     call exchange_add(relax_chm)

     !
     ! Variables for CMC combustion model
     !
     call exchange_add(kfl_weigh_in_eq_CMC_chm)
     call exchange_add(kfl_split_CFD_CMC)
     call exchange_add(kfl_solve_enth_CMC_chm)
     call exchange_add(kfl_solve_cond_CMC_chm)
     call exchange_add(kfl_bc_init_method_CMC_chm)
     call exchange_add(kfl_trans_phs_spc_CMC_chm)
     call exchange_add(kfl_trans_mxt_spc_CMC_chm)
     call exchange_add(kfl_mesh_interp_CMC_chm)
     call exchange_add(kfl_avg_cond_CMC_chm)
     call exchange_add(kfl_post_gp_CMC_chm)
     call exchange_add(kfl_incl_PDF_trans_CMC_chm)
     call exchange_add(kfl_dt_calc_CMC_chm)
     call exchange_add(kfl_bc_alpha_CMC_chm)
     call exchange_add(nZ_CMC_chm)
     call exchange_add(nvar_CMC_chm)
     call exchange_add(nvar_therm_CMC_chm)
     call exchange_add(nZ_AMC_CMC_chm)
     call exchange_add(nS_AMC_CMC_chm)
     call exchange_add(index_N2)
     call exchange_add(posZ_chem_integr_CMC_chm(1))
     call exchange_add(posZ_chem_integr_CMC_chm(2))
     call exchange_add(pos_Zstq_CMC_chm)
     call exchange_add(nspec_uncond_write_CMC_chm)
     call exchange_add(nspec_cond_write_CMC_chm)
     call exchange_add(nZ_write_CMC_chm)
     call exchange_add(nZ_chm_int_CMC_chm)
     call exchange_add(nZ_no_chm_int_CMC_chm)
     call exchange_add(size_tncod_CMC_chm)
     call exchange_add(Zs_CMC_chm)
     call exchange_add(Zstq_CMC_chm)
     call exchange_add(Smax_AMC_CMC_chm)
     call exchange_add(S_threshold)
     call exchange_add(extremes_Z_CMC_chm(1))
     call exchange_add(extremes_Z_CMC_chm(2))
     call exchange_add(order_RK_CMC_chm)

     !----------------------------------------------------------------
     !
     ! Exchange data read in chm_reaous
     !
     !----------------------------------------------------------------

     do ki=1,npara_chm
        call exchange_add(ipara_chm(ki))
     end do
     do ki=1,npara_chm
        call exchange_add(rpara_chm(ki))
     end do
     !
     ! Solver and postprocess
     !
     call solver_parall()
     call output_postprocess_parall()

     call exchange_end()
     
     !-------------------------------------------------------------------
     !
     ! Exchange tables
     !
     !-------------------------------------------------------------------
     
     call tab_par_exchange(table_coords,   table_tab   )

     if(kfl_tab_fw_chm > 0) then
        if (kfl_multimod_chm /= 1_ip) then
           table_fw => lookup_fw(kfl_tab_fw_chm)
        end if
     endif

     if(kfl_post_fw_chm > 0) then
        posttable_fw => lookup_fw(kfl_post_fw_chm)
     endif
     !
     ! Allocatable arrays
     !
     if( ISLAVE ) then
        call chm_memphy(1_ip) ! Allocate: Class dependent properties
        call chm_memphy(9_ip) ! Allocate equation structure
        call chm_memphy(10_ip) ! Allocate control variable association
        if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1) then
           ! Slaves allocate memory for the matrices whose size was obtained by the master
           call chm_memphy(3_ip)
           call chm_memphy(5_ip)
           call chm_memphy(6_ip)
           call chm_memphy(7_ip)
           call chm_memphy(8_ip)
           if (kfl_bc_init_method_CMC_chm == 2_ip) then
              call chm_memphy(11_ip)
           end if
           call chm_memnut(1_ip)
           call chm_memous(1_ip)
           call chm_memous(2_ip)
           call chm_memous(3_ip)
        end if
     end if

     !-------------------------------------------------------------------
     !
     ! Exchange equation data
     !
     !-------------------------------------------------------------------

     call chm_mixedEq_par_exchange(ngrou_chm,nclas_chm,mixedEq_groups_chm,mixedEq_eqs_chm)
     !
     ! Broadcast LREAC_CHM
     !
     if (kfl_model_chm == 3 .or. (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1)) then
        if( INOTMASTER ) &
             allocate(character(len=nsize_mech_name) :: mechanism_path)
     endif

     if (kfl_model_chm == 3) then
        if( INOTMASTER ) &
             allocate(character(len=nsize_red) :: Red_spec)
     endif

     !
     ! Allocate Field index in slaves
     !
     if (kfl_model_chm == 3 .and. kfl_spec_name_chm > 0_ip) then
        if( INOTMASTER ) &
             allocate(Field_ind_chm(kfl_spec_name_chm))
     endif

     !----------------------------------------------------------------
     !
     ! Exchange of chm_reaphy variables whose dimensions depend
     ! on what is read in chm_reaphy
     !
     !----------------------------------------------------------------

     !
     ! Finite-rate chemistry or CMC models
     !
     if (kfl_model_chm == 3 .or. (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1)) then
        call exchange_init()     
        call exchange_add(mechanism_path)
        call exchange_end()     
#ifdef CANTERA
        if( INOTMASTER ) &
             gas_chm = importPhase(mechanism_path)
#endif
     end if
     
     call exchange_init()    
     !
     ! Finite-rate chemistry model
     !
     if (kfl_model_chm == 3) then

        call exchange_add(Red_spec)

        do ki=1,kfl_spec_name_chm
           call exchange_add(Field_ind_chm(ki))
        end do
     end if
     !
     ! Variables from CMC combustion model
     !
     if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1) then
        if (kfl_bc_init_method_CMC_chm == 2_ip) then
           call exchange_add(inv_exp_alpha_grid_CMC_chm)
        end if
        do ki = 1,nZ_CMC_chm
           call exchange_add(Z_CMC_chm(ki))
           call exchange_add(Xnormalized_prof_CMC_chm(ki))
           call exchange_add(PDF_min_CMC_chm(ki))
        end do
        do ki = 1,nZ_CMC_chm-1
           call exchange_add(diff_Z_CMC_chm(ki))
        end do
        do ki = 1,nZ_CMC_chm
           call exchange_add(temp_inert_CMC_chm(ki))
           call exchange_add(temp_equil_CMC_chm(ki))
           do ki2 = 1,nspec_chm+1
              call exchange_add(rscal_inert_CMC_chm(ki,ki2))
              call exchange_add(rscal_equil_CMC_chm(ki,ki2))
           end do
        end do

        ! Variables for AMC model
        do ki = 1,nZ_AMC_CMC_chm
           call exchange_add(Z_AMC_CMC_chm(ki))
           do ki2 = 1,nS_AMC_CMC_chm
              call exchange_add(Xintegrated_table_AMC_CMC_chm(ki,ki2))
           end do
        end do
        do ki2 = 1,nS_AMC_CMC_chm
           call exchange_add(S_AMC_CMC_chm(ki2))
        end do

        call exchange_add(Zavg_const_CMC_chm)
        call exchange_add(Zvar_const_CMC_chm)
        call exchange_add(Xtot_const_CMC_chm)

        do ki = 1,2
           call exchange_add(Zlim_clipZvar_CMC_chm(ki))
           call exchange_add(Slim_clipZvar_CMC_chm(ki))
        end do

        do ki = 1,nspec_chm
           call exchange_add(max_diff_Yk_inert_eq_CMC_chm(ki))
        end do

        if (nspec_uncond_write_CMC_chm == nclas_chm) then
           if (ISLAVE)   write_uncond_spec_CMC_chm = (/( ki, ki=1_ip, nspec_uncond_write_CMC_chm )/)
        else
           do ki = 1, nspec_uncond_write_CMC_chm
              call exchange_add(write_uncond_spec_CMC_chm(ki))
           end do
        end if

        if (nspec_cond_write_CMC_chm == nclas_chm) then
           if (ISLAVE)   write_cond_spec_CMC_chm = (/( ki, ki=1_ip, nspec_cond_write_CMC_chm )/)
        else
           do ki = 1, nspec_cond_write_CMC_chm
              call exchange_add(write_cond_spec_CMC_chm(ki))
           end do
        end if

        if (nZ_write_CMC_chm == nZ_CMC_chm) then
           if (ISLAVE)   write_iZ_spec_CMC_chm = (/( ki, ki=1_ip, nZ_write_CMC_chm )/)
        else
           do ki = 1,nZ_write_CMC_chm
              call exchange_add(write_iZ_spec_CMC_chm(ki))
           end do
        end if

        do ki = 1, nZ_chm_int_CMC_chm
           call exchange_add(chem_int_iZ_CMC_chm(ki))
        end do

     end if

     do iclas=1,nspec_chm
        call exchange_add(diffu_chm(1,iclas))
        call exchange_add(diffu_chm(2,iclas))
        call exchange_add(Le_k(iclas))
     end do

     do iclas=1,nspec_chm ! Transfer of species type
        call exchange_add(speci(iclas)%visco)
        call exchange_add(speci(iclas)%lawvi)
        call exchange_add(speci(iclas)%weigh)
        call exchange_add(speci(iclas)%densi)
        call exchange_add(speci(iclas)%entha)
        call exchange_add(speci(iclas)%prand)
        call exchange_add(speci(iclas)%lewis)
        call exchange_add(speci(iclas)%trang)
        call exchange_add(speci(iclas)%cpcoe)
        call exchange_add(speci(iclas)%activ)
     end do
     call exchange_end()    

     if (kfl_model_chm == 4 .and. kfl_split_CFD_CMC == 1_ip) then
        call chm_find_master_CMC(1_ip,r_master)
        call chm_transfer_info_CMC_to_CFD(r_master)
        call chm_find_master_CMC(0_ip,r_master)
        call chm_transfer_info_CFD_to_CMC(r_master)
        call chm_write_info_transfer_CMC
     end if

  case(2_ip)

     !-------------------------------------------------------------------
     !
     ! Boundary conditions
     !
     !-------------------------------------------------------------------

     call boundary_conditions_exchange(tncod_chm,INCLUDE_CHARACTER=.false.)
     call boundary_conditions_exchange(tbcod_chm)

     call exchange_init()    
     call exchange_add(kfl_allcl_chm)
     call exchange_add(kfl_fields_scale_chm)
     call exchange_add(bo_oxy_chm)
     call exchange_add(bf_fuel_chm)

     do iclas = 1,nclas_chm
        call exchange_add(kfl_initi_chm(iclas))
     end do

     do iclas = 1,nclas_chm
        call exchange_add(xinit_chm(iclas,1))
     end do
     do iclas = 1,nclas_chm
        call exchange_add(xinit_chm(iclas,2))
     end do

     call exchange_end()

  end select

end subroutine chm_parall
