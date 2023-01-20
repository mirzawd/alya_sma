!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine exm_sendat(itask)
   !-----------------------------------------------------------------------
   !****f* exmedi/exm_sendat
   ! NAME
   !    exm_sendat
   ! DESCRIPTION
   !    This routine exchange data 
   ! USES
   ! USED BY
   !    exm_turnon
   !***
   !-----------------------------------------------------------------------
   use      def_parame
   use      def_master
   use      def_domain
   use      def_inpout
   use      mod_memchk
   use      mod_opebcs
   use      def_exmedi
   use      mod_exchange,            only : exchange_init
   use      mod_exchange,            only : exchange_add
   use      mod_exchange,            only : exchange_end
   use      mod_solver,              only : solver_parall
   use      mod_output_postprocess,  only : output_postprocess_parall 
   use      mod_memory,              only : memory_alloca, memory_deallo, memory_size
   use      mod_exm_fitzhugh_nagumo, only : exm_fitzhugh_nagumo_SendData
   use      mod_exm_drugs,           only : exm_drugs_exchange
   use      mod_eccoupling 
   use      mod_exm_activation,      only : exm_stim_exchange, exm_stim_allocate
   use      mod_exm_ecg,             only : exm_ecg_exchange, exm_ecg_allocate
       
   implicit none
   integer(ip), intent(in) :: itask 


   select case(itask)
   case (0_ip)
      call exchange_add(kfl_appty_exm)
      call exchange_add(kfl_appva_exm)
      call exchange_add(kfl_gcoup_exm)
      call exchange_add(kfl_cemod_exm)
      call exchange_add(kfl_stree_exm)
      call exchange_add(kfl_nodif_exm)
      call exchange_add(kfl_paced_exm)
      call exchange_add(ndofn_exm)
      call exchange_add(ndof2_exm)
      call exchange_add(nmate_exm)
      call exchange_add(nevat_exm)
      call exchange_add(ncomp_exm)
      call exchange_add(kfl_user_specified_celltypes_exm)
      call exchange_add(kfl_ignore_steadystate_celltypes_exm)
      call exchange_add(kfl_steadystate_variable)
      call exchange_add(steady_tol_cellm)
      call exchange_add(timestep_cellm)
      
      call exchange_add(nrootecg_exm)
      call exchange_add(kcopeecg_exm)
      call exchange_add(nvint_exm)

      call exchange_add(kfl_fract_diffusion_exm)
      call exchange_add(kfl_ortho_diffusion_exm)
      call exchange_add(kfl_isac_exm)
      !call exchange_add(kfl_hfmod_exm)
      call exchange_add(kfl_hfmodmate_exm)
      call exchange_add(kfl_active_material)
      call exchange_add(kfl_reset_fisoc)
      call exchange_add(xmccmmate_exm)
      call exchange_add(moneclmate_exm)
      call exchange_add(ttparmate_exm)
           
      call exchange_add(nauxi_exm)
      call exchange_add(nconc_exm)
      call exchange_add(nicel_exm)
      call exchange_add(nstim_exm)
      call exchange_add(nstis_exm)
      call exchange_add(modab_exm)
      call exchange_add(modst_exm)
      call exchange_add(dtinv_exm)
      call exchange_add(fisoc_exm)               ! Isochrones trigger
      call exchange_add(gdiff_exm)

      call exchange_add(tcardiac_cycle)          ! Master variable, but exmedi controls it and solidz uses it

      call exchange_add(nstrb_exm)
      call exchange_add(strbo_exm)
      call exchange_add(fiaxe_exm)
      call exchange_add(stran_endo_exm)
      call exchange_add(stran_epi_exm)
      call exchange_add(vminimate_exm)
      call exchange_add(kfl_voini_exm)          
      call exchange_add(kfl_fract_diffusion_exm)
      call exchange_add(fract_diff_coef_exm)
      call exchange_add(fract_diff_nintp_exm)
      call exchange_end()

     
      call exm_fitzhugh_nagumo_SendData()
      call exm_drugs_exchange()

      !
      ! Exchange of exm_reanut variables 
      !
      call exchange_init()
      call exchange_add(kfl_genal_exm)              ! General alg. type
      call exchange_add(kfl_goite_exm)              ! Keep iterating
      call exchange_add(kfl_shock_exm)              ! 
      call exchange_add(kfl_comat_exm)              ! 
      call exchange_add(kfl_weigh_exm)              ! 
      call exchange_add(kfl_timet_exm)              ! 
      call exchange_add(kfl_tiacc_exm)              ! 
      call exchange_add(kfl_ticel_exm)              ! 
      call exchange_add(kfl_tisch_exm)              ! 
      call exchange_add(kfl_normc_exm)              ! 
      call exchange_add(kfl_adres_exm)              ! 
      call exchange_add(kfl_algso_exm)              ! 
      call exchange_add(kfl_repro_exm)              ! 
      call exchange_add(kfl_nolim_exm)              ! 
      call exchange_add(kfl_nolum_exm)              ! 
      call exchange_add(miinn_exm)              ! 
      call exchange_add(msste_exm)              ! 
      call exchange_add(mnoli_exm)              ! 
      call exchange_add(msoit_exm)              ! 
      call exchange_add(nkryd_exm)              ! 
      call exchange_add(itera_exm)              ! 
      call exchange_add(nunkn_exm)              ! 
      call exchange_add(dtcri_exm    )              ! 
      call exchange_add(shock_exm    )              ! 
      call exchange_add(sstol_exm    )              ! 
      call exchange_add(cotol_exm    )              ! 
      call exchange_add(corat_exm    )              ! 
      call exchange_add(dtext_exm    )              ! 
      call exchange_add(safet_exm    )              ! 
      call exchange_add(solco_exm    )              ! 
      call exchange_add(weigh_exm    )              ! 
      call exchange_add(tnoli_exm    )              ! 
      call exchange_add(err01_exm)           ! L1 error u
      call exchange_add(err02_exm)           ! L2 error u
      call exchange_add(err0i_exm)           ! Linf error u
      call exchange_add(err11_exm)           ! L1 error grad(u)
      call exchange_add(err12_exm)           ! L2 error grad(u)
      call exchange_add(err1i_exm)           ! Linf error grad(u)
      call exchange_add(cpu_exmed)           ! CPU for the exm problem
      call exchange_add(staco_exm)           ! 
      call exchange_add(resid_exm)           ! Residual for outer iterations (Mom,Cont,Ene,Glo)

      call exchange_add(kfl_exboc_exm)       ! Boundary conditions explicitly given
      
      !
      ! Exchange data read in exm_reaous
      !
      call exchange_add(kfl_save_convergence_cellmodel)               
      call exchange_add(kfl_save_init_cellmodel)    
      
      call output_postprocess_parall()        
      
      !
      ! Solvers
      !
      call solver_parall()
      

      call exchange_end()
      
      call eccou_send_cellmod()

      !
      ! Stumulus table, after exchanging nstim_exm, allocate arrays on slaves and exchange
      !
      if(ISLAVE) call exm_stim_allocate()
      call exm_stim_exchange()

      !
      ! ECG table
      !
      if(ISLAVE) call exm_ecg_allocate()
      call exm_ecg_exchange()


      if (kfl_exboc_exm == 1) then
          call boundary_conditions_exchange(tncod_exm)
          call boundary_conditions_exchange(tbcod_exm)
      end if

  case(1_ip) !after exm_inivar

    call exchange_init()
    call exchange_add(ncelltypes_per_model)
    call exchange_end()

  case default
    call runend("exm_sendat: unknown task id")
  end select

end subroutine exm_sendat
