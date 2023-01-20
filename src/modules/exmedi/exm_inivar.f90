!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_inivar.f90
!> @author  Mariano Vazquez
!> @brief   Initialize data
!> @date    16/11/1966
!> @details Initialize data
!> @} 
!-----------------------------------------------------------------------
subroutine exm_inivar(itask)
  !-----------------------------------------------------------------------
  !****f* Temper/exm_inivar
  ! NAME 
  !    exm_inicar
  ! DESCRIPTION
  !    This routine initializes some variables
  !    ITASK=1 ... When starting the run (from Turnon)
  !    ITASK=2 ... First time step. This is needed as some variables 
  !                are not initialized before
  ! USES
  ! USED BY
  !    exm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_exmedi 
  use def_solver
  use mod_iofile
  use mod_arrays,              only : arrays_register
  use mod_memory,              only : memory_alloca, memory_size
  use mod_exm_fitzhugh_nagumo, only : exm_fitzhugh_nagumo_GetNumberOfVariables, exm_fitzhugh_nagumo_InitVariables
  use mod_moduls_conf,         only : moduls_allocate_timing
  use mod_exm_cellmodel
  use mod_exm_drugs,           only : exm_drugs_allocate
  use mod_eccoupling
  use mod_moduls_conf,         only : moduls_allocate_timing
  use mod_exm_ohararudy,       only : exm_ohara_init

  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: iconc,iauxi
  integer(ip)             :: imate, icelltype, iionconcentr
  integer(ip)             :: nvars_activation, nvars_concentrations

  select case(itask)

  case(0_ip)
    !
    ! Register variables
    !
    call arrays_register((/'TAULO','SCALA','NPOIN','PRIMA'/),taulo         ,ENTITY_POSITION=1_ip)
    call arrays_register((/'FISOC','SCALA','NPOIN','PRIMA'/),fisoc         ,COMPONENT_POSITION=1_ip, ENTITY_POSITION=2_ip)
    call arrays_register((/'VMIMA','SCALA','NPOIN','PRIMA'/),elmag_minmax  ,COMPONENT_POSITION=1_ip, ENTITY_POSITION=2_ip)
    call arrays_register((/'ELMAG','SCALA','NPOIN','PRIMA'/),elmag         ,ENTITY_POSITION=1_ip,    TIME_POSITION=2_ip, EXCLUDE_RST=(/2_ip,3_ip/))
    call arrays_register((/'VCONC','SCALA','NPOIN','PRIMA'/),vconc         ,COMPONENT_POSITION=1_ip, ENTITY_POSITION=2_ip, TIME_POSITION=3_ip, EXCLUDE_RST=(/2_ip,3_ip/))
    ! 12_ip used to be an array of potentials for fitzhugh-nagumo, not used anymore
    call arrays_register((/'QNETT','SCALA','NPOIN','PRIMA'/),qneto_exm     ,ENTITY_POSITION=1_ip)
    call arrays_register((/'VAUXI','SCALA','NPOIN','PRIMA'/),vauxi_exm     ,COMPONENT_POSITION=1_ip, ENTITY_POSITION=2_ip, TIME_POSITION=3_ip, EXCLUDE_RST=(/2_ip,3_ip/))
    call arrays_register((/'VICEL','SCALA','NPOIN','PRIMA'/),vicel_exm     ,COMPONENT_POSITION=1_ip, ENTITY_POSITION=2_ip )
    call arrays_register((/'FIBER','VECTO','NPOIN','SECON'/) )
    call arrays_register((/'GRAFI','VECTO','NPOIN','SECON'/) )
    call arrays_register((/'IOCON','MULTI','NPOIN','SECON'/) )

    call arrays_register((/'CALCI','SCALA','NPOIN','SECON'/) )
    call arrays_register((/'POTAS','SCALA','NPOIN','SECON'/) )
    call arrays_register((/'EFLUX','VECTO','NPOIN','SECON'/) )
    call arrays_register((/'BIPOL','VECTO','NPOIN','SECON'/) )
    call arrays_register((/'INTRA','SCALA','NPOIN','SECON'/) )
    call arrays_register((/'EXTRA','SCALA','NPOIN','SECON'/) )
    call arrays_register((/'RECOV','SCALA','NPOIN','SECON'/) )
    call arrays_register((/'SHEET','SCALA','NPOIN','SECON'/) )
    call arrays_register((/'NORMA','SCALA','NPOIN','SECON'/) )
    call arrays_register((/'FREE4','SCALA','NPOIN','SECON'/) )
    call arrays_register((/'CECO1','VECTO','NPOIN','SECON'/) )
    call arrays_register((/'CECO2','VECTO','NPOIN','SECON'/) )
    call arrays_register((/'CECO3','VECTO','NPOIN','SECON'/) )
    call arrays_register((/'DISPL','VECTO','NPOIN','SECON'/) )
    call arrays_register((/'CURRE','MULTI','NPOIN','SECON'/) )
    call arrays_register((/'QNETO','SCALA','NPOIN','SECON'/) )
    call arrays_register((/'ISOCH','SCALA','NPOIN','SECON'/) )
    call arrays_register((/'REPOL','SCALA','NPOIN','SECON'/) )


    do iconc= 1,nconc_exm
      if (iconc<=9) then
         call arrays_register((/'VCO0'//trim(intost(iconc)),'SCALA','NPOIN','SECON'/) )
      else
         call arrays_register((/'VCO'//trim(intost(iconc)),'SCALA','NPOIN','SECON'/) )
      end if
    end do

    do iauxi= 1,nauxi_exm
      if (iconc<=9) then
         call arrays_register((/'VAU0'//trim(intost(iconc)),'SCALA','NPOIN','SECON'/) )
      else
         call arrays_register((/'VAU'//trim(intost(iconc)),'SCALA','NPOIN','SECON'/) )
      end if
    end do

    !
    ! Sets variables
    !

    postp(1) % woese (  1) = 'INTRA'
    postp(1) % woese (  2) = 'CALCI'         
    postp(1) % woese (  3) = 'SODIU'         
    postp(1) % woese (  4) = 'POTAS'         
    postp(1) % woese ( 17) = 'EFLUX'         
    postp(1) % woese ( 18) = 'BIPOL'         

    postp(1) % wonse (  1) = 'INTRA'
    postp(1) % wonse (  2) = 'CALCI'   
    postp(1) % wonse (  3) = 'SODIU'   
    postp(1) % wonse (  4) = 'POTAS'   
    postp(1) % wonse ( 17) = 'EFLUX'   
    postp(1) % wonse ( 18) = 'BIPOL'   

    !
    ! Witness variables
    !
    postp(1) % wowit  (1)     = 'INTRA'
    postp(1) % wowit  (2)     = 'COORX'
    postp(1) % wowit  (3)     = 'CAI  '  !VCONC(1:11) ! 2+1
    postp(1) % wowit  (4)     = 'NASS '
    postp(1) % wowit  (5)     = 'KI   '
    postp(1) % wowit  (6)     = 'KISS '
    postp(1) % wowit  (7)     = 'NAI  '
    postp(1) % wowit  (8)     = 'CASS '
    postp(1) % wowit  (9)     = 'CANSR'
    postp(1) % wowit (10)     = 'CAJSR'
    postp(1) % wowit (11)     = 'JRLNP' !Jrelnp
    postp(1) % wowit (12)     = 'JRLP ' !Jrelp
    postp(1) % wowit (13)     = 'CAMKT'
    postp(1) % wowit (14)     = 'INA  '  !VICEL(1:26) ! 13+1
    postp(1) % wowit (15)     = 'INAL '
    postp(1) % wowit (16)     = 'ITO  '
    postp(1) % wowit (17)     = 'ICAL '
    postp(1) % wowit (18)     = 'IKR  '
    postp(1) % wowit (19)     = 'IKS  '
    postp(1) % wowit (20)     = 'IK1  '
    postp(1) % wowit (21)     = 'INCI ' !INaCa_i
    postp(1) % wowit (22)     = 'INCSS' !INaCa_ss
    postp(1) % wowit (23)     = 'INAK '
    postp(1) % wowit (24)     = 'IKB  '
    postp(1) % wowit (25)     = 'INAB '
    postp(1) % wowit (26)     = 'ICAB '
    postp(1) % wowit (27)     = 'IPCA '
    postp(1) % wowit (28)     = 'JDI  ' !Jdiff
    postp(1) % wowit (29)     = 'JDINA' !JdiffNa
    postp(1) % wowit (30)     = 'JDIK '
    postp(1) % wowit (31)     = 'JUP '
    postp(1) % wowit (32)     = 'JLEAK'
    postp(1) % wowit (33)     = 'JTR  '
    postp(1) % wowit (34)     = 'JREL '
    postp(1) % wowit (35)     = 'CAMKA'
    postp(1) % wowit (36)     = 'STIM '
    postp(1) % wowit (37)     = 'ICANA'
    postp(1) % wowit (38)     = 'ICAK '
    postp(1) % wowit (39)     = 'CAMKB'
    postp(1) % wowit (40)     = 'ISOCH'
    postp(1) % wowit (41)     = 'REPOL'

!!!!!!     postp(1) % wowit (40)     = 'EFLUX'
!!!!!!     postp(1) % wowit (41)     = 'BIPOL'


    !
    ! Nullify arrays
    !
    nullify(idima_exm)  
    nullify(kgrfi_exm)        
    nullify(cedif_exm)
    nullify(grafi_exm)
    nullify(fiber_exm)
    nullify(sheet_exm)
    nullify(normal_exm)
    nullify(celty_exm)
    nullify(atbhe_exm)
    nullify(vdiag_exm)     
    nullify(fibe2_exm) 
    nullify(tncod_exm)      
    nullify(tbcod_exm)      
    nullify(amatr_auxi_exm) 
    nullify(appfi_exm)      
    nullify(vauxi_exm)
    nullify(ticel_exm)    
    nullify(jicel_exm)    
    nullify(vicel_exm)
    nullify(qneto_exm)

    nullify(kfl_voini_exm)
    !nullify(kfl_hfmod_exm)
    nullify(kfl_hfmodmate_exm)
    nullify(kfl_steadystate_variable)
    nullify(kfl_fract_diffusion_exm)
    nullify(kfl_isac_exm)
    nullify(xmccmmate_exm)
    nullify(fract_diff_coef_exm)
    nullify(steady_tol_cellm)
    nullify(timestep_cellm)
    nullify(fract_diff_nintp_exm)
    nullify(moneclmate_exm)
    nullify(vminimate_exm)
    nullify(gdiff_exm)
    nullify(ttparmate_exm)
    nullify(vauxi_exm_initial)
    nullify(vconc_initial)
    nullify(kfl_user_specified_celltypes_exm)
    nullify(kfl_ignore_steadystate_celltypes_exm)
    nullify(kfl_active_material)

    nullify(eflux_exm)

    nullify(bipol_exm)
    nullify(ncelltypes_per_model)

    ! Allocate arrays
    !These are used in reaphy to read data
    call memory_alloca(mem_modul(1:2,modul),'kfl_steadystate_variable',  'exm_inivar', kfl_steadystate_variable,  nmate)
    call memory_alloca(mem_modul(1:2,modul),'ttparmate_exm',             'exm_inivar', ttparmate_exm,             3_ip,  nvars_ttparmate_exm, nmate)
    call memory_alloca(mem_modul(1:2,modul),'steady_tol_cellm',          'exm_inivar', steady_tol_cellm,          nmate)
    call memory_alloca(mem_modul(1:2,modul),'timestep_cellm',            'exm_inivar', timestep_cellm,            nmate)
    call memory_alloca(mem_modul(1:2,modul),'kfl_voini_exm',             'exm_inivar', kfl_voini_exm,             nmate)
    call memory_alloca(mem_modul(1:2,modul),'kfl_fract_diffusion_exm',   'exm_inivar', kfl_fract_diffusion_exm,   nmate)
    
    !call memory_alloca(mem_modul(1:2,modul),'kfl_hfmod_exm',             'exm_inivar', kfl_hfmod_exm,             nmate)
    call memory_alloca(mem_modul(1:2,modul),'kfl_hfmodmate_exm',         'exm_inivar', kfl_hfmodmate_exm,         nmate)
    call memory_alloca(mem_modul(1:2,modul),'kfl_isac_exm',              'exm_inivar', kfl_isac_exm,              nmate)
    call memory_alloca(mem_modul(1:2,modul),'xmccmmate_exm',             'exm_inivar', xmccmmate_exm,             nmate)
    call memory_alloca(mem_modul(1:2,modul),'fract_diff_coef_exm',       'exm_inivar', fract_diff_coef_exm,       nmate)
    call memory_alloca(mem_modul(1:2,modul),'fract_diff_nintp_exm',      'exm_inivar', fract_diff_nintp_exm,      nmate)
    call memory_alloca(mem_modul(1:2,modul),'moneclmate_exm',            'exm_inivar', moneclmate_exm,            2_ip,         nmate)
    call memory_alloca(mem_modul(1:2,modul),'gdiff_exm',                 'exm_inivar', gdiff_exm,                 2_ip,  3_ip,  nmate)
    call memory_alloca(mem_modul(1:2,modul),'kfl_ortho_diffusion_exm',   'exm_inivar', kfl_ortho_diffusion_exm,   nmate)
    call memory_alloca(mem_modul(1:2,modul),'kfl_active_material',       'exm_inivar', kfl_active_material,       nmate)

    call memory_alloca(mem_modul(1:2,modul),'ncelltypes_per_model',      'exm_inivar', ncelltypes_per_model,      EXMSLD_CELL_MAXMODELID)
    call exm_init_ncelltypes(ncelltypes_per_model)

    call exm_drugs_allocate()


    call memory_alloca(mem_modul(1:2,modul),'vminimate_exm',                        'exm_inivar', vminimate_exm, &
                      ncelltypes_ecc, nmate)
    call memory_alloca(mem_modul(1:2,modul),'kfl_user_specified_celltypes_exm',     'exm_inivar', kfl_user_specified_celltypes_exm,&
                      nmate, ncelltypes_ecc) !TODO: flip dimensions for consistency
    call memory_alloca(mem_modul(1:2,modul),'kfl_ignore_steadystate_celltypes_exm', 'exm_inivar', kfl_ignore_steadystate_celltypes_exm, &
                      nmate, ncelltypes_ecc) !TODO: flip dimensions for consistency

    do icelltype=1_ip, size(vminimate_exm,1,kind=ip)
      do imate=1_ip, size(vminimate_exm,2,kind=ip)
          vminimate_exm(icelltype, imate) = 0.0_rp
      end do
    end do

 
 
 
    kfl_user_specified_celltypes_exm(:,:) = 0_ip


    !
    ! Pseudo-ecg
    !
    kcopeecg_exm = 0
    !
    ! Solvers
    !     
    call soldef(-1_ip)

    solve(1) % wprob     = 'ACTIV_POTENTIAL' ! length should be <26 otherwise timers fail
    solve(1) % kfl_solve = 1
    solve(1) % ndofn     = 1

    ! In the case of implicit, these values will be defined in nsa_reanut, when calling reasol
    ! These are default values, corresponding to explicit with global time step
    solve(1) % kfl_algso = SOL_SOLVER_RICHARDSON            
    solve(1) % kfl_preco = SOL_CLOSE_MASS_MATRIX

    ! Time variables
    kfl_timei = 1      ! exmedi is always transient

    !
    ! Materials
    !
    nmate_exm =  nmate
    !lcell_exm => lcell

    epres = 0.0_rp

    cpold_exm= 0.0_rp

    do imate=1,size(kfl_user_specified_celltypes_exm, 1, KIND=ip)
      do icelltype=1,size(kfl_user_specified_celltypes_exm, 2, KIND=ip)
          kfl_user_specified_celltypes_exm( imate, icelltype ) = 1_ip   !by deafult the ODE for all celltypes will be executed
          kfl_ignore_steadystate_celltypes_exm( imate, icelltype ) = 0_ip !do not ignore steady state by deafult
      end do
      kfl_steadystate_variable(imate) = EXM_CELL_STEADY_VOLTAGE !use voltage for steady state by default
      steady_tol_cellm(imate) = -1.0_rp                !let subroutine determine the optimal value
      timestep_cellm  (imate) = -1.0_rp                !let subroutine determine the optimal value
    end do


    !reset ttparmate_exm to 1, afterwards in reaphy it will be filled
    do icelltype=1,size(ttparmate_exm, 1, KIND=ip)
      do iionconcentr=1,size(ttparmate_exm, 2, KIND=ip)
          do imate=1,size(ttparmate_exm, 3, KIND=ip)
            ttparmate_exm(icelltype,iionconcentr,imate) = 1_rp
          end do
      end do
    end do
    ttparmate_exm(:, ohara_conductance_ikatp_row, :) = 0.0_rp   !IKatp current should only be non-zero if ischemia

    !
    ! Flags
    !
    kfl_save_convergence_cellmodel = 0_ip
    kfl_save_init_cellmodel = 0_ip
    kfl_reset_fisoc = .FALSE.

    !
    ! Exmedi needs to know if the electromechanical coupling is one or multi-level
    ! 
    if( kfl_exmsld_ecc )then
      call eccou_initialise_flags()
    endif

    !
    ! Initialize the variables of the cell models
    !
    call exm_fitzhugh_nagumo_InitVariables()
    !
    ! Timings
    !
    call moduls_allocate_timing(2_ip)
    call times(1) % register('stimuli'   ,'step beginning')
    call times(2) % register('time step' ,'step beginning')

  case ( 1_ip ) ! After reading .dat files
    !
    ! More initialisations after reading files
    !
    ncomp_exm = 2 + kfl_tiacc_exm
    ! 
    !
    ! FOR SUBCELLULAR IONIC CURRENTS MODELS (TT, LR, BR, ...)
    !


    nauxi_exm = 1
    nconc_exm = 1
    nicel_exm = 1

    do imate= 1,nmate

      select case(kfl_cellmod(imate))
      case (EXMSLD_CELL_NOMODEL) 
          nauxi_exm = max(nauxi_exm,1_ip)    ! Default value of number of variables used for activation/inactivation 
          nconc_exm = max(nconc_exm,1_ip)    ! Default value of number of variables used for concentrations
          nicel_exm = max(nicel_exm,1_ip)
      case (EXMSLD_CELL_FENTON) 
          call runend("EXM_INIVAR: FENTON KARMA MODELS NOT UPDATED")
          nconc_exm = 2
      case (EXMSLD_CELL_FITZHUGH) 
          call exm_fitzhugh_nagumo_GetNumberOfVariables(imate, nvars_activation, nvars_concentrations)
          nauxi_exm  = max( nauxi_exm, nvars_activation)  
          nconc_exm  = max( nconc_exm, nvars_concentrations)
          !else if (kfl_cellmod(imate) == CELL_TT2006_EXMEDI) then
          !   nauxi_exm = max(nauxi_exm,12_ip)  
          !   nconc_exm = max(nconc_exm,10_ip)  
          !   nicel_exm = max(nconc_exm,18_ip) 
      case (EXMSLD_CELL_OHARA, EXMSLD_CELL_OHARA_INAPA)
          nauxi_exm  = max( nauxi_exm,29_ip)  
          nconc_exm  = max( nconc_exm,11_ip)  
          nicel_exm  = max( nconc_exm,26_ip) 
      case (EXMSLD_CELL_TORORD) 
          nauxi_exm  = max( nauxi_exm,31_ip)
          nconc_exm  = max( nconc_exm,14_ip)
          nicel_exm  = max( nicel_exm,21_ip)
      case (EXMSLD_CELL_SCATRIA) 
          nauxi_exm  = max( nauxi_exm,14_ip)  
          nconc_exm  = max( nconc_exm,3_ip)  
          nicel_exm  = max( nconc_exm,15_ip)
      case (EXMSLD_CELL_SCVENTRI) 
          nauxi_exm = max(nauxi_exm,14_ip)  
          nconc_exm = max(nconc_exm,3_ip)  
          nicel_exm = max(nconc_exm,15_ip)
      case (EXMSLD_CELL_COURTE) 
          nauxi_exm = max(nauxi_exm,15_ip)                          
          nconc_exm = max(nconc_exm,5_ip)                           
          nicel_exm = max(nicel_exm,18_ip)     
          !     else if (kfl_cellmod(imate) == 6) then  !TTINA
          !        nauxi_exm = max(nauxi_exm,14_ip)  
          !        nconc_exm = max(nconc_exm,10_ip)  
          !        nicel_exm = max(nconc_exm,18_ip) 
      end select
    end do

    do imate=1,nmate
      if( kfl_cellmod(imate)>0 ) then
        if ( ncelltypes_ecc<ncelltypes_per_model(kfl_cellmod(imate)) ) then
          call runend("Cell model of material "//trim(intost(imate))//" requires "//trim(intost(ncelltypes_per_model(kfl_cellmod(imate))))//" celltypes instead of the specified in CELL_TYPE field ("//trim(intost(ncelltypes_ecc))// ")")
        end if
      end if
    end do

    call memory_alloca(mem_modul(1:2,modul),'vauxi_exm_initial', 'exm_inivar', vauxi_exm_initial, nauxi_exm, ncelltypes_ecc,  nmate)
    call memory_alloca(mem_modul(1:2,modul),'vconc_initial',     'exm_inivar', vconc_initial,     nconc_exm, ncelltypes_ecc,  nmate)

    !-----------------------------------------------------------------------
    !
    !   Initialize the models if needed
    !
    !-----------------------------------------------------------------------
    if( any( kfl_cellmod(1:nmate) == EXMSLD_CELL_OHARA ) .or. any( kfl_cellmod(1:nmate) == EXMSLD_CELL_OHARA_INAPA ) ) then 
        call exm_ohara_init()
    end if
    !-----------------------------------------------------------------------
    !
    !   End of Initialize the model if needed
    !
    !-----------------------------------------------------------------------

    ! Timers
    !call moduls_allocate_timing(2_ip)
    !call times(1) % register("precalc" ,"iteration other", "computation")
    !call times(2) % register("calc gates" ,"iteration other", "computation")


  end select

end subroutine exm_inivar
