!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_sendat.f90
!> @author  Solidz Team
!> @date    November, 2017-Adds doxygen comments
!> @brief   Exchange data
!> @details It is used by sld_turnon
!> @}
!-----------------------------------------------------------------------

subroutine sld_sendat()

  use def_kintyp,                only : ip, rp
  use def_master,                only : IMASTER, ISLAVE
  use mod_memchk
  use mod_opebcs
  use mod_sld_cardiac_cycle,     only : sld_cardiac_cycle_send_data
  use mod_exchange,              only : exchange_init
  use mod_exchange,              only : exchange_add
  use mod_exchange,              only : exchange_end
  use mod_solver,                only : solver_parall
  use mod_output_postprocess,    only : output_postprocess_parall
  use mod_sld_cardiac_cycle,     only : sld_cardiac_cycle_send_data
  use mod_sld_post_reaction,     only : sld_post_reaction_parall
  use mod_sld_strain,            only : strain_basis
  use mod_sld_atm,               only : sld_atm_parall
  use def_solidz
#if defined COMMDOM && COMMDOM == 2
  use mod_sld_pdn_contact_plepp, only : sld_pdn_contact_parall
#endif
#ifndef PROPER_ELEM_PRIVATE_OFF
  use mod_ker_sysnet,            only : sysnet_sendat
#endif
  
  implicit none

  external :: sld_memphy
  external :: sld_membcs

  integer(ip) :: ifunc

  !----------------------------------------------------------------
  !
  ! Exchange variables read in sld_reaphy, sld_reanut and sld_reaous
  !
  !----------------------------------------------------------------

  call exchange_init()
  
  call exchange_add(kfl_timei_sld)    
  call exchange_add(kfl_rigid_sld)    
  call exchange_add(kfl_tange_sld)  
  call exchange_add(kfl_strai_sld)
  call exchange_add(kfl_sdvar_sld)    
  call exchange_add(kfl_fiber_sld)    
  call exchange_add(kfl_csysm_sld)
  call exchange_add(kfl_rmate_sld)
  call exchange_add(kfl_prdef_sld)
  call exchange_add(kfl_cohes_sld)    
  call exchange_add(kfl_nopassive_sld)    
  call exchange_add(ncalcio_sld)   
  call exchange_add(kfl_cshel_sld)
  call exchange_add(kfl_damag_sld)
  call exchange_add(kfl_savdt_sld)
  call exchange_add(kfl_dttyp_sld)
  call exchange_add(kfl_celen_sld)
  call exchange_add(kfl_vofor_sld)         
  call exchange_add(kfl_plane_sld)       
  call exchange_add(kfl_restr_sld)    
  call exchange_add(kfl_indis_sld)           
  call exchange_add(kfl_plast_sld)    
  call exchange_add(kfl_vecto_sld)
  call exchange_add(kfl_moduf_sld)
  call exchange_add(kfl_donna_sld)
  call exchange_add(nmate_sld)
  call exchange_add(oripa_sld)
  call exchange_add(gfibe_sld)
  call exchange_add(gravi_sld)
  call exchange_add(grnor_sld)
  call exchange_add(rbmas_sld)
  call exchange_add(thick_sld)
  call exchange_add(thiso_sld)
  call exchange_add(csysm_sld)
  call exchange_add(calmusico_sld)
  call sld_cardiac_cycle_send_data(1_ip) ! Cardiac cycle
  call exchange_add(kfl_xfeme_sld)
  call exchange_add(kfl_xfcra_sld)
  call exchange_add(kfl_stabi_sld)
  call exchange_add(kfl_resid_sld)
  call exchange_add(kfl_timet_sld)
  call exchange_add(kfl_ninex_sld)
  call exchange_add(kfl_tisch_sld)
  call exchange_add(kfl_serei_sld)
  call exchange_add(kfl_limit_sld)
  call exchange_add(kfl_safet_table_sld)
  call exchange_add(dtmin_sld)
  call exchange_add(    miinn_sld)
  call exchange_add(    minex_sld)
  call exchange_add(kfl_plane_sld)
  call exchange_add(kfl_gdepo)
  call exchange_add(kfl_ellen_sld)
  
  call exchange_add(tifac_sld(1))
  call exchange_add(tifac_sld(2))
  call exchange_add(tifac_sld(3))
  call exchange_add(tifac_sld(4))
  call exchange_add(tifac_sld(5))

  call exchange_add(cotol_sld)
  call exchange_add(resid_sld)
  call exchange_add(safet_sld)
  call exchange_add(safex_sld)              
  call exchange_add(safma_sld)  
  call exchange_add(epsex_sld)
  call exchange_add(safet_pseud_sld)
  call exchange_add(factor_penal_sld)
  call exchange_add(dafac_sld)
  call exchange_add(sstol_sld)
  call exchange_add(masss_sld)
  call exchange_add(nisaf_sld    )              ! Initial time step for variable cfl
  call exchange_add(kfl_pseud_sld)
  call exchange_add(kfl_penal_sld)
  call exchange_add(safet_table_sld)   ! Safety factor table: factor
  call exchange_add(rorig_sld)

#ifndef PROPER_ELEM_PRIVATE_OFF
  call sysnet_sendat(1_ip) ! Sysnet
#endif  
  
  call solver_parall()
  call strain_basis % send_data(end_exchange = .false.)

  call exchange_end()

  !----------------------------------------------------------------
  !
  ! Exchange material properties read in sld_reaphy
  !
  !----------------------------------------------------------------

  !
  ! Allocatable related to materials
  !
  if( ISLAVE ) call sld_memphy(1_ip)

  call exchange_init()
  
  call exchange_add(modfi_sld)
  call exchange_add(modor_sld)
  call exchange_add(kfl_dampi_sld)
  call exchange_add(lawco_sld)
  call exchange_add(lawst_sld)
  call exchange_add(lawho_sld)
  call exchange_add(lawch_sld)
  call exchange_add(lawpl_sld)
  call exchange_add(lawta_sld)

  call exchange_add(parsp_sld)
  call exchange_add(densi_sld)
  call exchange_add(parco_sld)
  call exchange_add(parch_sld)
  call exchange_add(parcf_sld)
  call exchange_add(dampi_sld)

  call exchange_add(stiff0_sld)
  
  call exchange_end()
  
  !-------------------------------------------------------------------
  !
  ! Cardiac cycle and sysnet: require variables exchanged previously
  !
  !-------------------------------------------------------------------

  call exchange_init()
  
#ifndef PROPER_ELEM_PRIVATE_OFF
  call sysnet_sendat(2_ip) ! Sysnet
#endif  
  
  call sld_cardiac_cycle_send_data(2_ip)
  call exchange_end()

  !-------------------------------------------------------------------
  !
  ! Arrays and dimensions read in sld_reaous
  !
  !-------------------------------------------------------------------
  
  call exchange_init()
  
  call exchange_add(kfl_exacs_sld)
  call exchange_add(kfl_foten_sld)
  call exchange_add(kfl_rotei_sld)
  call exchange_add(rorig_sld)
  call output_postprocess_parall()
  
  call exchange_end()
  
  call sld_post_reaction_parall()
  
  !-------------------------------------------------------------------
  !
  ! Dimensions read in sld_reabcs
  !
  !-------------------------------------------------------------------

  call exchange_init()
  
  call exchange_add(kfl_conbc_sld)
  call exchange_add(kfl_newbc_sld)
  call exchange_add(kfl_local_sld)
  call exchange_add(kfl_csysl_sld)
  call exchange_add(kfl_follo_sld)
  call exchange_add(kfl_bvesv_sld)
  call exchange_add(kfl_invel_sld)
  call exchange_add(kfl_insdv_sld)
  call exchange_add(kfl_funty_sld)
  call exchange_add(kfl_bodyf_sld)
  call exchange_add(kfl_conta_sld)
  call exchange_add(kfl_contf_sld)
  call exchange_add(kfl_mrele_sld)
  call exchange_add(nfunc_sld)
  call exchange_add(ncrak_sld)
  call exchange_add(contactbou_sld)
  call exchange_add(coupling_contact_its)
  call exchange_add(kfl_conta_stent)
  call exchange_add(conbo_sld)
  call exchange_add(r_fin_stent)
  call exchange_add(contact_friction_sld)
  call exchange_add(neumann_relax)
  call exchange_add(coupling_contact_tol)
  call exchange_add(vect_proje_sld)
  call exchange_add(petol_sld)
  call exchange_add(csysl_sld)
  call exchange_add(mtloa_sld)
  call exchange_add(rtico_sld)
  call exchange_add(fubcs_sld)
  call exchange_add(invel_sld)
#if defined COMMDOM && COMMDOM==2
  call sld_pdn_contact_parall()
#endif
  
  call exchange_end()

  !-------------------------------------------------------------------
  !
  ! Thermomechanical coupling
  !
  !-------------------------------------------------------------------
  
  call sld_atm_parall()
    
  !-------------------------------------------------------------------
  !
  ! Others
  !
  !-------------------------------------------------------------------

  if( ISLAVE .and. ncrak_sld > 0 ) call sld_membcs(31_ip) ! Cracks

  if( ISLAVE ) then
     call sld_membcs(3_ip)
     do ifunc=1,10
        if (mtloa_sld(ifunc) == 0) mtloa_sld(ifunc) = 1            ! done to allocate a default memory space
        call sld_membcs(10_ip + ifunc)  ! allocate the prescription time function vector for ifunc
     end do
  end if

  if( associated(tload_sld) ) then
     do ifunc = 1,size(tload_sld,kind=ip)
        call exchange_add(tload_sld(ifunc)%a)
     end do
  end if
  
  if( associated(crkco_sld) ) call exchange_add(crkco_sld)
  call exchange_end()

  if( IMASTER .and. ncrak_sld > 0 ) call sld_membcs(-31_ip) ! Cracks

  !----------------------------------------------------------------
  !
  ! Exchange boundary conditions
  !
  !----------------------------------------------------------------

  call boundary_conditions_exchange(tncod_sld)
  call boundary_conditions_exchange(tbcod_sld)
  
end subroutine sld_sendat
