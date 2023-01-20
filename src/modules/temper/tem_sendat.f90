!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_sendat()
  !-----------------------------------------------------------------------
  !****f* Temper/tem_sendat
  ! NAME
  !    tem_sendat
  ! DESCRIPTION
  !    This routine exchange data 
  ! USES
  ! USED BY
  !    tem_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_solver
  use def_temper
  use mod_opebcs
  use def_kermod,             only : gasco
  use mod_exchange,           only : exchange_init
  use mod_exchange,           only : exchange_add
  use mod_exchange,           only : exchange_end
  use mod_solver,             only : solver_parall
  use mod_output_postprocess, only : output_postprocess_parall

  implicit none

  call exchange_init()
  
  call exchange_add(kfl_timei_tem)
  call exchange_add(kfl_advec_tem)
  call exchange_add(kfl_joule_tem)
  call exchange_add(kfl_radia_tem)
  call exchange_add(kfl_sourc_tem)
  call exchange_add(kfl_sonum_tem)
  call exchange_add(kfl_adiab_tem)
  call exchange_add(kfl_condu_tem)
  call exchange_add(kfl_exint_tem)
  call exchange_add(kfl_inter_tem)
  call exchange_add(kfl_dynco_tem)
  call exchange_add(kfl_regim_tem)
  call exchange_add(kfl_lookg_tem)
  call exchange_add(kfl_diven_tem)
  call exchange_add(kfl_parti_tem)
  call exchange_add(idtem_tem)
  call exchange_add(kfl_skews_tem)
  call exchange_add(kfl_nudgi_tem)
  call exchange_add(kfl_forest_tem)
  call exchange_add(kfl_entropy_tem)
  call exchange_add(kfl_disable_entpre)
  call exchange_add(kfl_code_tem)

  call exchange_add(turbu_tem)
  call exchange_add(prthe_tem)
  call exchange_add(cfi_hmax_tem)
  call exchange_add(cfi_hmin_tem)
  call exchange_add(prtur_tem)
  call exchange_add(scond_tem)
  call exchange_add(react_tem)
  call exchange_add(sourc_tem)
  call exchange_add(nudgi_tem)
  call exchange_add(gasco)

  call exchange_add(lsour_material_tem)
  call exchange_add(xsour_material_tem)
  !
  ! Exchange of tem_reanut variables 
  !        
  call exchange_add(kfl_dttyp_tem)
  call exchange_add(kfl_ellen_tem)
  call exchange_add(kfl_sgsti_tem)
  call exchange_add(kfl_sgsno_tem)
  call exchange_add(kfl_taust_tem)
  call exchange_add(kfl_ortho_tem)
  call exchange_add(kfl_limit_tem)
  call exchange_add(kfl_shock_tem)
  call exchange_add(kfl_tiacc_tem)
  call exchange_add(kfl_tibub_tem)
  call exchange_add(kfl_assem_tem)
  call exchange_add(kfl_posit_tem)
  call exchange_add(kfl_negat_tem)
  call exchange_add(neule_tem)
  call exchange_add(kfl_tisch_tem)
  call exchange_add(kfl_normc_tem)
  call exchange_add(miinn_tem)
  call exchange_add(misgs_tem)
  call exchange_add(kfl_sgsac_tem)
  call exchange_add(kfl_meshi_tem)
  call exchange_add(kfl_discr_tem)
  call exchange_add(kfl_sgsli_tem)
  call exchange_add(kfl_plepp_tem)
  call exchange_add(kfl_explicit_tem)
  call exchange_add(kfl_rhs_scal_tem)

  call exchange_add(staco_tem)
  call exchange_add(shock_tem)
  call exchange_add(safet_tem)
  call exchange_add(source_safet_tem)
  call exchange_add(sstol_tem)
  call exchange_add(cotol_tem)
  call exchange_add(relax_tem)
  call exchange_add(bemol_tem)
  call exchange_add(relsg_tem)
  call exchange_add(tosgs_tem)
  !
  ! Exchange data read in tem_reaous
  !
  call exchange_add(kfl_splot_tem)
  call exchange_add(kfl_psmat_tem)
  call exchange_add(kfl_exacs_tem)
  call exchange_add(kfl_viewf_tem)
  call exchange_add(npp_bound_tem)
  call exchange_add(avtim_tem)
  call exchange_add(expar_tem)
  !
  ! Solver and postprocess
  !
  call solver_parall()
  call output_postprocess_parall()
  
  call exchange_end()

  !------------------------------------------------------------------- 
  !
  ! Variables read in reabcs
  !
  !------------------------------------------------------------------- 

  call boundary_conditions_exchange(tncod_tem)
  call boundary_conditions_exchange(tbcod_tem)

  call exchange_init()
  call exchange_add(kfl_conbc_tem) 
  call exchange_add(kfl_inidi_tem)
  call exchange_add(kfl_inico_tem)
  call exchange_add(kfl_intbc_tem)
  call exchange_add(npnat_tem)
  call exchange_add(delta_tem)   
  call exchange_add(initial_tem)
  call exchange_add(kfl_funty_tem)
  call exchange_add(funpa_tem)
  call exchange_end()

end subroutine tem_sendat
