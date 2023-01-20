!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_sendat.f90
!> @author  houzeaux
!> @date    2020-11-27
!> @brief   Broadcast data
!> @details Broadcast data read from files by master
!> @} 
!-----------------------------------------------------------------------

subroutine nsi_sendat()

  use def_master
  use def_domain
  use def_nastin
  use mod_opebcs,             only : boundary_conditions_exchange
  use def_kermod,             only : gasco
  use mod_vortex,             only : vortex_sendat   
  use mod_exchange,           only : exchange_init
  use mod_exchange,           only : exchange_add
  use mod_exchange,           only : exchange_end
  use mod_solver,             only : solver_parall
  use mod_output_postprocess, only : output_postprocess_parall
  use mod_vortex,             only : vortex_parall

  implicit none
  
  !----------------------------------------------------------------
  !
  ! Exchange of nsi_reaphy variables 
  !
  !----------------------------------------------------------------

  call exchange_init()
  
  call exchange_add(kfl_timei_nsi)
  call exchange_add(kfl_advec_nsi)
  call exchange_add(kfl_convection_type_nsi)       
  call exchange_add(kfl_fvfua_nsi)
  call exchange_add(kfl_fvful_nsi)
  call exchange_add(kfl_cotem_nsi)
  call exchange_add(kfl_grtur_nsi)
  call exchange_add(kfl_visco_nsi)
  call exchange_add(kfl_colev_nsi)
  call exchange_add(kfl_regim_nsi)
  call exchange_add(kfl_dynco_nsi)
  call exchange_add(kfl_prthe_nsi)
  call exchange_add(kfl_surte_nsi)
  call exchange_add(kfl_force_nsi)
  call exchange_add(kfl_vegeo_time_nsi)
  call exchange_add(kfl_mfrco_nsi)
  call exchange_add(kfl_bnods_nsi)
  call exchange_add(kfl_hydro_gravity_nsi)
  call exchange_add(kfl_anipo_nsi)
  call exchange_add(kfl_fscon_nsi)
  call exchange_add(kfl_hydro_interface_nsi)
  call exchange_add(mfrse_nsi)
  call exchange_add(ntabf_nsi)
  call exchange_add(nbnod_nsi)
  call exchange_add(nbtim_nod_nsi)
  call exchange_add(nfiel_nsi)    

  call exchange_add(fcons_nsi)
  call exchange_add(fvins_nsi)
  call exchange_add(grnor_nsi)
  call exchange_add(fvnoa_nsi)
  call exchange_add(fanoa_nsi)
  call exchange_add(fvnol_nsi)
  call exchange_add(fanol_nsi)
  call exchange_add(bougr_nsi)
  call exchange_add(boube_nsi)
  call exchange_add(boutr_nsi)
  call exchange_add(lowtr_nsi)
  call exchange_add(turbu_nsi)
  call exchange_add(heihy_nsi)
  call exchange_add(surte_nsi)
  call exchange_add(tmass_nsi)
  call exchange_add(mfrub_nsi)
  call exchange_add(ubpre_nsi)
  call exchange_add(mfccf_nsi)

  call exchange_add(gravi_nsi)
  call exchange_add(gravb_nsi)
  call exchange_add(fvdia_nsi)
  call exchange_add(fvela_nsi)
  call exchange_add(fadia_nsi)
  call exchange_add(facca_nsi)
  call exchange_add(fvdil_nsi)
  call exchange_add(fvell_nsi)
  call exchange_add(fadil_nsi)
  call exchange_add(faccl_nsi)
  call exchange_add(frotc_nsi)
  call exchange_add(centr_nsi)
  call exchange_add(fvpaa_nsi)
  call exchange_add(fvpal_nsi)

  !----------------------------------------------------------------
  !
  ! Exchange of nsi_reanut variables 
  !
  !----------------------------------------------------------------

  call exchange_add(kfl_penal_nsi)
  call exchange_add(kfl_prepe_nsi)
  call exchange_add(kfl_dttyp_nsi)
  call exchange_add(kfl_ellen_nsi)
  call exchange_add(kfl_relax_nsi)
  call exchange_add(kfl_relap_nsi)
  call exchange_add(kfl_sgsco_nsi)
  call exchange_add(kfl_sgsti_nsi)
  call exchange_add(kfl_sgsac_nsi)
  call exchange_add(kfl_sgsli_nsi)
  call exchange_add(kfl_sgscp_nsi)
  call exchange_add(kfl_shock_nsi)
  call exchange_add(kfl_tiacc_nsi)
  call exchange_add(kfl_normc_nsi)
  call exchange_add(kfl_refer_nsi)
  call exchange_add(kfl_linea_nsi)
  call exchange_add(kfl_tisch_nsi)
  call exchange_add(kfl_algor_nsi)
  call exchange_add(kfl_predi_nsi)
  call exchange_add(kfl_taush_nsi)
  call exchange_add(kfl_ellsh_nsi)
  call exchange_add(kfl_updpr_nsi)
  call exchange_add(kfl_intpr_nsi)
  call exchange_add(kfl_assem_nsi)
  call exchange_add(kfl_asbou_nsi)
  call exchange_add(kfl_taust_nsi)
  call exchange_add(kfl_stabi_nsi)
  call exchange_add(kfl_limit_nsi)
  call exchange_add(kfl_trres_nsi)
  call exchange_add(kfl_prtre_nsi)
  call exchange_add(kfl_matdi_nsi)
  call exchange_add(kfl_intfo_nsi)
  call exchange_add(kfl_press_nsi)
  call exchange_add(kfl_bubbl_nsi)
  call exchange_add(kfl_fsgrb_nsi)
  call exchange_add(kfl_stab_corio_nsi)
  call exchange_add(momod(modul) % miinn)
  call exchange_add(kfl_stain_nsi)
  call exchange_add(kfl_immer_nsi)
  call exchange_add(kfl_velom_nsi)
  call exchange_add(kfl_grad_div_nsi)
  call exchange_add(misgs_nsi)
  call exchange_add(npica_nsi)
  call exchange_add(itinn(modul))
  call exchange_add(neule_nsi)
  call exchange_add(kfl_savco_nsi)
  call exchange_add(kfl_meshi_nsi)
  call exchange_add(kfl_corre_nsi)
  call exchange_add(kfl_sosch_nsi)
  call exchange_add(kfl_modfi_nsi)
  call exchange_add(kfl_expco_nsi)
  call exchange_add(kfl_addpr_nsi)
  call exchange_add(kfl_grvir_nsi)
  call exchange_add(kfl_hydro_nsi)
  call exchange_add(kfl_update_hydro_nsi)
  call exchange_add(kfl_hydro_interface_nsi)
  call exchange_add(mitri_nsi)
  call exchange_add(kfl_adres_nsi)
  call exchange_add(kfl_incre_nsi)
  call exchange_add(kfl_nota1_nsi)
  call exchange_add(kfl_ini_ts_guess_order_nsi)
  call exchange_add(kfl_vector_nsi)
  call exchange_add(kfl_press_stab_nsi)
  call exchange_add(kfl_stop_by_wit_nsi)
  call exchange_add(kfl_massm_nsi)
  call exchange_add(kfl_dampi_nsi)
  call exchange_add(k_dim_damp_nsi)

  call exchange_add(top_r_damp_nsi)
  call exchange_add(top_v_damp_nsi)
  call exchange_add(bot_r_damp_nsi)
  call exchange_add(bot_v_damp_nsi)
  call exchange_add(val_r_damp_nsi)
  call exchange_add(mul_v_damp_nsi)
  call exchange_add(v_geo_damp_nsi)

  call exchange_add(penal_nsi)
  call exchange_add(prepe_nsi)
  call exchange_add(dtcri_nsi)
  call exchange_add(staco_nsi)
  call exchange_add(shock_nsi)
  call exchange_add(safet_nsi)
  call exchange_add(bemol_nsi)
  call exchange_add(sstol_nsi)
  call exchange_add(cotol_nsi)
  call exchange_add(resid_nsi)
  call exchange_add(resip_nsi)
  call exchange_add(weigh_nsi)
  call exchange_add(relax_nsi)
  call exchange_add(relap_nsi)
  call exchange_add(relsg_nsi)
  call exchange_add(staco_corio_nsi)
  call exchange_add(tosgs_nsi)
  call exchange_add(strec_nsi)
  call exchange_add(dampi_nsi)
  call exchange_add(epsht_nsi)
  call exchange_add(epstr_nsi)
  call exchange_add(xfree_nsi)
  call exchange_add(safex_nsi)
  call exchange_add(adres_nsi)
  call exchange_add(toler_nsi)
  call exchange_add(safma_nsi)
  call exchange_add(safeo_nsi)
  call exchange_add(saflo_nsi)
  call exchange_add(gamma_nsi)
  call exchange_add(fsrot_nsi)

  call solver_parall()
  
  !----------------------------------------------------------------
  !
  ! Exchange data read in nsi_reabcs
  !
  !----------------------------------------------------------------

  call exchange_add(kfl_confi_nsi)
  call exchange_add(kfl_local_nsi)
  call exchange_add(kfl_conbc_nsi)
  call exchange_add(kfl_initi_nsi)
  call exchange_add(kfl_inico_nsi)
  call exchange_add(kfl_inipr_nsi)
  call exchange_add(kfl_nopen_nsi)
  call exchange_add(kfl_cadan_nsi)
  call exchange_add(kfl_syntu_nsi)
  call exchange_add(kfl_aiobo_nsi)
  call exchange_add(neddy_nsi)
  call exchange_add(itebc_nsi)
  call exchange_add(nodpr_global_nsi)
  call exchange_add(exfpr_nsi)
  call exchange_add(kfl_imppr_nsi)
  call exchange_add(kfl_divcorrec_nsi)
  call exchange_add(nflow_rates)
  call exchange_add(delta_nsi)
  call exchange_add(relbc_nsi)
  call exchange_add(valpr_nsi)
  call exchange_add(mulpr_nsi)
  call exchange_add(hydro_nsi)
  call exchange_add(velin_nsi)
  call exchange_add(velin_nsi)
  call exchange_add(poise_nsi)
  call exchange_add(divcorrec_nsi)
  call exchange_add(kfl_press_lev_nsi)
  call exchange_add(val_press_lev_nsi)
  call exchange_add(kfl_immersed_nsi)               
  call exchange_add(fact_immersed_nsi)  
  call exchange_add(inflow_cosine_nsi)
  call exchange_end()

  if( nflow_rates > 0 ) then
     call exchange_add(kfl_flow_rate_codes_nsi)
     call exchange_add(kfl_flow_rate_set_nsi)
     call exchange_add(kfl_flow_rate_normal_nsi)
     call exchange_add(kfl_flow_rate_stfn_nsi)
     call exchange_add(kfl_flow_rate_tfn_nsi)
     call exchange_add(kfl_flow_rate_pfn_nsi)
     call exchange_add(kfl_flow_rate_alg_nsi)
     call exchange_add(flow_rate_values_nsi)
     call exchange_add(flow_rate_press_nsi)
     call exchange_add(flow_rate_relax_nsi)
     call exchange_add(flow_rate_normals_nsi)
     call exchange_end()
  end if
  
  !----------------------------------------------------------------
  !
  ! Exchange data read in nsi_reaous
  !
  !----------------------------------------------------------------

  call exchange_add(kfl_exacs_nsi)
  call exchange_add(kfl_exfix_nsi)
  call exchange_add(kfl_inert_nsi)
  call exchange_add(kfl_psmat_nsi)
  call exchange_add(kfl_miniapp_nsi)
  call exchange_add(expar_nsi)
  call exchange_add(cloth_nsi)
  call exchange_add(metab_nsi)
  call exchange_add(wetme_nsi)
  call exchange_add(ambie_nsi)
  call exchange_add(radia_nsi)
  call exchange_add(relat_nsi)
  call exchange_add(avtim_nsi)
  call exchange_add(avste_nsi)
  call exchange_add(entim_nsi)  
  call output_postprocess_parall()
  call vortex_parall()
    
  call exchange_end()
  !
  ! Allocate properties memory for slaves
  !

  !
  ! Material force term
  !
  if( kfl_force_nsi == 1 .and. ISLAVE ) call nsi_memphy(5_ip) ! master was allocated before
  if( kfl_bnods_nsi == 1 .and. ISLAVE ) then
     call nsi_memphy(6_ip)
     call nsi_memphy(7_ip)
     call nsi_memphy(8_ip)
  end if

  !----------------------------------------------------------------
  !
  ! Exchange properties and material force read in nsi_reaphy
  !
  !----------------------------------------------------------------

  call exchange_init()
  call exchange_add(gasco)
  call exchange_add(sphea_nsi)
  call exchange_add(prthe_nsi)
  !
  ! Material force
  !
  if( kfl_force_nsi == 1 ) then
     call exchange_add(lforc_material_nsi)
     call exchange_add(ntabl_nsi)
     call exchange_add(ntabr_nsi)
     call exchange_add(xforc_material_nsi)
     call exchange_add(velta_nsi)
     call exchange_add(thrta_nsi)
     call exchange_add(powta_nsi)
     call exchange_add(veave_nsi)
     call exchange_add(radiu_nsi)
     call exchange_add(forcn_nsi)
     call exchange_add(forct_nsi)
  end if
  !
  ! exchange boundary nodes list
  !
  if( kfl_bnods_nsi == 1 ) then
     call exchange_add(nbnod_pos_nsi)
     call exchange_add(nbtim_pos_nsi)
     call exchange_add(nbtim_nod_pos_nsi)
     call exchange_add(nbtdt_nsi)
     call exchange_add(bntab_nsi)
     call exchange_add(bnval_nsi)
  end if

  call exchange_end()

  !----------------------------------------------------------------
  !
  ! Exchange boundary conditions
  !
  !----------------------------------------------------------------

  call boundary_conditions_exchange(tncod_nsi)
  call boundary_conditions_exchange(tscod_nsi)
  call boundary_conditions_exchange(tgcod_nsi)
  call boundary_conditions_exchange(tbcod_nsi)
  
end subroutine nsi_sendat
