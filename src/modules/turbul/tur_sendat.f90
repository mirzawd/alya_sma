!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_sendat(order)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_sendat
  ! NAME
  !    tur_sendat
  ! DESCRIPTION
  !    This routine exchange data 
  ! USES
  ! USED BY
  !    tur_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_solver
  use def_domain
  use def_turbul
  use def_inpout
  use mod_memchk
  use mod_opebcs
  use mod_exchange,           only : exchange_init
  use mod_exchange,           only : exchange_add
  use mod_exchange,           only : exchange_end
  use mod_solver,             only : solver_parall
  use mod_output_postprocess, only : output_postprocess_parall
  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: ji,jr,ki
  integer(4)              :: istat

  select case (order)

  case(1_ip)     
     !
     ! Exchange data read in tur_reaphy, tur_reanut and tur_reaous
     !
     call exchange_init()
     !
     ! Exchange of tur_reaphy variables 
     !
     call exchange_add(kfl_model_tur)
     call exchange_add(kfl_timei_tur)
     call exchange_add(kfl_cotem_tur)
     call exchange_add(kfl_advec_tur)
     call exchange_add(kfl_colev_tur)
     call exchange_add(kfl_ddesm_tur)
     call exchange_add(kfl_sasim_tur)
     call exchange_add(kfl_rotat_tur)
     call exchange_add(kfl_inifi_tur(1))           ! Initial fields 
     call exchange_add(kfl_inifi_tur(2))           ! Initial fields 
     call exchange_add(kfl_inifi_tur(3))           ! Initial fields
     do jr=1,3
        call exchange_add(nfiel_tur(jr))           ! Fields assignement
     end do

     call exchange_add(inits_tur)
     call exchange_add(nturb_tur)
     do ji=1,nipar_tur
        call exchange_add(ipara_tur(ji))
     end do
     call exchange_add(lawde_tur)
     call exchange_add(lawvi_tur)
     call exchange_add(kfl_kxmod_tur)
     call exchange_add(kfl_logva_tur)
     call exchange_add(kfl_discd_tur)
     call exchange_add(boube_tur)
     call exchange_add(grnor_tur)
     do jr=1,3
        call exchange_add(gravi_tur(jr))
     end do
     call exchange_add(prtur_tur)
     do jr=1,nrpar_tur
        call exchange_add(param_tur(jr))
     end do
     do jr=1,ncoef_tur
        call exchange_add(densi_tur(jr))
     end do
     do jr=1,ncoef_tur
        call exchange_add(visco_tur(jr))
     end do
     call exchange_add(densa_tur)
     call exchange_add(visca_tur)
     call exchange_add(cddes_tur)       
     call exchange_add(inv_l_max)       

     !
     ! Exchange of tur_reanut variables 
     !        
     call exchange_add(kfl_weigh_tur)
     call exchange_add(kfl_repro_tur)
     call exchange_add(kfl_taust_tur)
     call exchange_add(kfl_shock_tur)
     call exchange_add(kfl_algor_tur)
     call exchange_add(kfl_clipp_tur)
     call exchange_add(kfl_ellen_tur)
     call exchange_add(kfl_relax_tur)
     call exchange_add(kfl_tiacc_tur)
     call exchange_add(kfl_tisch_tur)
     call exchange_add(kfl_normc_tur)
     call exchange_add(kfl_walgo_tur)
     call exchange_add(kfl_assem_tur)
     call exchange_add(kfl_ortho_tur)
     call exchange_add(kfl_limit_tur)
     call exchange_add(kfl_produ_tur)
     call exchange_add(kfl_meshi_tur)

     call exchange_add(miinn_tur)
     call exchange_add(niter_tur)
     call exchange_add(neule_tur)

     call exchange_add(kfl_sgsti_tur)
     call exchange_add(kfl_sgsno_tur)
     call exchange_add(kfl_tibub_tur)
     call exchange_add(kfl_sgsac_tur)

     call exchange_add(staco_tur(1))
     call exchange_add(staco_tur(2))
     call exchange_add(staco_tur(3))
     call exchange_add(shock_tur)
     call exchange_add(sstol_tur)
     call exchange_add(safet_tur)
     call exchange_add(cotol_tur)
     call exchange_add(relax_tur)
     call exchange_add(safex_tur)
     call exchange_add(safma_tur)
     call exchange_add(safeo_tur) 
     call exchange_add(saflo_tur)
     call exchange_add(bemol_tur)

     call exchange_add(clipfac_tur)
     call exchange_add(kfl_lmaxi_tur)
     !
     ! Solver and postprocess
     !
     call solver_parall()
     call output_postprocess_parall()
     call exchange_end()     
     !
     ! Exchange of tur_reaphy variables 
     !
     if (kfl_discd_tur==1) then
        call exchange_init()
        do ji=1, nmate        
           call exchange_add(ldiss_material_tur(ji))
        end do
        call exchange_end()
     end if

  case(2_ip)     
     !
     ! Exchange data read in tur_reabcs
     !
     !
     ! Variables read in tur_reabcs
     !
     call exchange_init()
     !
     ! Exchange of tur_reabcs variables 
     !
     call exchange_add(kfl_inidi_tur)
     call exchange_add(kfl_wallw_tur)
     call exchange_add(kfl_infl1_tur)
     call exchange_add(kfl_infl2_tur)
     call exchange_add(kfl_usrbc_tur)
     call exchange_add(kfl_initi_tur)
     do ki = 1,4
        call exchange_add(kfl_valbc_tur(ki))
     end do
     call exchange_add(delta_tur)
     call exchange_add(turin_tur)
     call exchange_add(hdiam_tur)
     call exchange_add(turle_tur) 
     call exchange_add(nutnu_tur)
     call exchange_add(rebcs_tur)
     do ki = 1,4
        call exchange_add(xinit_tur(ki))
     end do
     call exchange_end()
     !
     ! Boundary codes
     !
     call boundary_conditions_exchange(tncod_tur)
     call boundary_conditions_exchange(tgcod_tur)
     call boundary_conditions_exchange(tbcod_tur)

  end select

end subroutine tur_sendat
