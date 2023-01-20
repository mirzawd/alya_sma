!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_inivar(itask)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_inivar
  ! NAME 
  !    tem_inicar
  ! DESCRIPTION
  !    This routine initializes some variables
  !    ITASK=1 ... When starting the run (from Turnon)
  !    ITASK=2 ... First time step. This is needed as some variables 
  !                are not initialized before
  ! USES
  ! USED BY
  !    tem_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use def_solver
  use mod_arrays,          only : arrays_register
  use mod_ADR,             only : ADR_initialize_type
  use mod_tem_rk_explicit, only : tem_rk_explicit_initialization
  use mod_tem_explicit,    only : tem_explicit_initialization
  use def_kermod,          only : turmu_ker
  use mod_temper,          only : temper_initialization
  use mod_temper,          only : temper_memory_allocate
  implicit none
  integer(ip), intent(in) :: itask

  select case(itask)

  case(0_ip)
     !
     ! Intialization
     !
     call temper_initialization() 
     call ADR_initialize_type(ADR_tem)
     call tem_rk_explicit_initialization()
     call tem_explicit_initialization()
     !
     ! Postprocess
     !
     call arrays_register((/'TEMPE','SCALA','NPOIN','PRIMA'/),tempe           ,ENTITY_POSITION=1_ip,TIME_POSITION=2_ip,EXCLUDE_RST=(/2_ip,3_ip/))
     call arrays_register((/'THERM','SCALA','NPOIN','PRIMA'/),therm           ,ENTITY_POSITION=1_ip,TIME_POSITION=2_ip,EXCLUDE_RST=(/2_ip,3_ip/))
     call arrays_register((/'RESID','SCALA','NPOIN','PRIMA'/),teold_tem       ,ENTITY_POSITION=1_ip)
     call arrays_register((/'AVTEM','SCALA','NPOIN','PRIMA'/),avtem_tem       ,ENTITY_POSITION=1_ip)
     call arrays_register((/'FVTEM','SCALA','NPOIN','PRIMA'/),fvtem_tem       ,ENTITY_POSITION=1_ip)
     call arrays_register((/'FVTE2','SCALA','NPOIN','PRIMA'/),fvte2_tem       ,ENTITY_POSITION=1_ip)
     call arrays_register((/'AVTE2','SCALA','NPOIN','PRIMA'/),avte2_tem       ,ENTITY_POSITION=1_ip)                    ! average tempe*tempe
     call arrays_register((/'AVTEV','VECTO','NPOIN','PRIMA'/),avtev_tem       ,ENTITY_POSITION=2_ip)                    ! average veloc*tempe
     call arrays_register((/'AVDEN','SCALA','NPOIN','PRIMA'/),avden_tem       ,ENTITY_POSITION=1_ip)                    ! average density
     call arrays_register((/'FVVEL','VECTO','NPOIN','PRIMA'/),fvvel_tem       ,ENTITY_POSITION=2_ip)                    ! Favre average velocity
     call arrays_register((/'AVRES','SCALA','NPOIN','PRIMA'/),avres_tem       ,ENTITY_POSITION=1_ip)                    ! average Heat flux using residuals
     call arrays_register((/'PROJ1','SCALA','NPOIN','PRIMA'/),ADR_tem % proje1,ENTITY_POSITION=1_ip)
     call arrays_register((/'PROJ2','SCALA','NPOIN','PRIMA'/),ADR_tem % proje2,ENTITY_POSITION=1_ip)
     call arrays_register((/'TESG2','R3PVE','NELEM','PRIMA'/),ADR_tem % sgs   ,ENTITY_POSITION=1_ip,TIME_POSITION=3_ip)
     call arrays_register((/'BUBBT','SCALA','NELEM','PRIMA'/),ADR_tem % bubble,ENTITY_POSITION=1_ip,TIME_POSITION=2_ip)
     call arrays_register((/'AVHSK','SCALA','NPOIN','PRIMA'/),avhsk_tem       ,ENTITY_POSITION=1_ip)                    ! average spray heat source 

     call arrays_register((/'VELOC','VECTO','NPOIN','SECON'/))
     call arrays_register((/'TURVI','SCALA','NPOIN','SECON'/))
     call arrays_register((/'HEATF','SCALA','NPOIN','SECON'/))
     call arrays_register((/'TESGS','SCALA','NPOIN','SECON'/))
     call arrays_register((/'ERROR','SCALA','NPOIN','SECON'/))
     call arrays_register((/'GROUP','SCALA','NPOIN','SECON'/))
     call arrays_register((/'LIMIT','SCALA','NPOIN','SECON'/))
     call arrays_register((/'LINTE','SCALA','NPOIN','SECON'/))
     call arrays_register((/'TESTS','VECTO','NPOIN','SECON'/))
     call arrays_register((/'WATVA','SCALA','NPOIN','SECON'/))
     call arrays_register((/'FIXTE','SCALA','NPOIN','SECON'/))
     call arrays_register((/'TEMP2','SCALA','NPOIN','SECON'/))
     call arrays_register((/'HEATN','SCALA','NPOIN','SECON'/))
     call arrays_register((/'GRATE','VECTO','NPOIN','SECON'/))
     call arrays_register((/'CHARA','VECTO','NPOIN','SECON'/))
     call arrays_register((/'BVNAT','SCALA','NPOIN','SECON'/))
     call arrays_register((/'ENTHA','SCALA','NPOIN','SECON'/))
     call arrays_register((/'REACT','SCALA','NPOIN','SECON'/))
     call arrays_register((/'-----','SCALA','NPOIN','SECON'/))
     call arrays_register((/'-----','SCALA','NPOIN','SECON'/))
     call arrays_register((/'RESHE','SCALA','NPOIN','SECON'/))
     call arrays_register((/'NFIXN','SCALA','NPOIN','SECON'/))
     call arrays_register((/'-----','SCALA','NPOIN','SECON'/))
     call arrays_register((/'ENVIT','SCALA','NPOIN','SECON'/))
     call arrays_register((/'HEASK','SCALA','NPOIN','SECON'/))
     call arrays_register((/'HEATS','SCALA','NPOIN','SECON'/),heat_source)
     !
     ! Element sets
     !
     postp(1) % woese(1)     = 'MASS '  ! int_W rho dW:   mass in domain 
     postp(1) % woese(2)     = 'ENTHA'  ! int_W rho h dW: enthalpy in domain 
     postp(1) % woese(3)     = 'INTER'  ! int_W rho h - Pth*rho  dW: internal energy in domain 
     postp(1) % woese(4)     = 'HEATS'  ! int_W q: heat srouce

     !
     ! Boundary sets
     !
     postp(1) % wobse(1)     = 'MEANT'  ! temper  (boundary set)  
     postp(1) % wobse(2)     = 'MEANH'  ! heat flux (-k grad(T)·n) (boundary set)
     postp(1) % wobse(3)     = 'MHEAN'  ! variational heatflx (boundary set)
     postp(1) % wobse(4)     = 'DENSI'  ! density on boundary 

     !
     ! Nodal sets
     !
     postp(1) % wonse(1)     = 'TEMPE'  ! nodal set

     !
     ! Witness point variables
     !
     postp(1) % wowit(1)     = 'TEMPE'  
     postp(1) % wowit(2)     = 'FLUXX'
     postp(1) % wowit(3)     = 'FLUXY'
     postp(1) % wowit(4)     = 'FLUXZ'
     postp(1) % wowit(5)     = 'ENTHA'
     !
     ! Solver
     !     
     call soldef(-2_ip)
     solve(1) % wprob     = 'TEMPERATURE'          ! Equation name
     solve(1) % kfl_solve = 1                      ! Output flag

     solve(2) % wprob     = 'CONSISTENT'           ! Consistent matrix
     solve(2) % kfl_solve = 1
     solve(2) % kfl_iffix = 2
     !
     ! Allocate memory
     !
     call temper_memory_allocate('LSOUR_MATERIAL_TEM')
     call temper_memory_allocate('XSOUR_MATERIAL_TEM')

  case(1_ip)

     if(kfl_timei_tem==0) then                     ! Time integration
        dtinv_tem=1.0_rp
     else
        kfl_timei=1
        dtinv_tem=0.0_rp
     end if
     kfl_stead_tem=0

     if( turmu_ker % kfl_exist /= 0_ip ) then  ! Recommended not to use mu_t gradients 
        kfl_grdif_tem=0
     else 
        kfl_grdif_tem=1          ! It is recommnded for laminar flows
     end if

     kfl_tiaor_tem=kfl_tiacc_tem                    ! Time accuracy: save original value
     if(kfl_advec_tem==0) bemol_tem=0.0_rp          ! There is no advection: BEMOL_TEM=0
     if(kfl_timei_tem==1) then                      ! Number of velocity components
        if(kfl_tisch_tem==1) then
           ncomp_tem=3                              ! Trapezoidal rule
        else if(kfl_tisch_tem==2) then
           ncomp_tem=2+kfl_tiacc_tem                ! BDF scheme
        else if(kfl_tisch_tem==3) then
           ncomp_tem=2+2                              ! AB scheme
        else if(kfl_tisch_tem==4) then
           ncomp_tem=2+2                              ! RK scheme
        end if
     else
        ncomp_tem = 2     
     end if

     ittot_tem = 0                                  ! Others
     dpthe     = 0.0_rp                             ! Low-Mach: dp0/dt

     !
     ! To solve equation of state when nastin & chemic are not activated
     !
     if ( kfl_regim_tem >= 3 .and. prthe_tem /= 0 ) then 
        prthe(1)      = prthe_tem                       ! Low-Mach: p0
        prthe(2)      = prthe_tem                       ! Low-Mach: p0^{n}
        prthe(3)      = prthe_tem                       ! Low-Mach: p0^{n-1}
        prthe(4)      = prthe_tem                       ! Low-Mach: p0^{0}
     end if
     xmass_tem =0.0_rp                              ! Low-Mach: mass
     kfl_rstar_two= .false.                         ! restarting bdf file?     
     !
     ! Finite volume
     !
     if(      kfl_discr_tem == 0 ) then
        solve(1) % kfl_where = SOL_NODES
        nunkn_tem            = npoin
     else if( kfl_discr_tem == 1 ) then
        solve(1) % kfl_where = SOL_ELEMENTS
        nunkn_tem            = nelem
     end if
     !
     ! Clipping
     !
     posit_tem = -huge(1.0_rp)
     negat_tem =  huge(1.0_rp)     
     !
     ! FILL in ADR type
     !
     ADR_tem % kfl_time_integration   =  kfl_timei_tem
     ADR_tem % kfl_time_step_strategy =  kfl_timco
     ADR_tem % kfl_stabilization      =  kfl_ortho_tem
     ADR_tem % kfl_shock              =  kfl_shock_tem
     ADR_tem % kfl_time_lumped        =  0
     ADR_tem % kfl_tau_strategy       =  kfl_taust_tem
     ADR_tem % kfl_laplacian          =  0 
     ADR_tem % kfl_nonlinear_sgs      =  kfl_sgsno_tem
     ADR_tem % kfl_time_sgs           =  kfl_sgsti_tem
     ADR_tem % kfl_time_bubble        =  kfl_tibub_tem
     ADR_tem % kfl_time_scheme        =  kfl_tisch_tem
     ADR_tem % kfl_time_order         =  kfl_tiacc_tem
     ADR_tem % kfl_manufactured       =  kfl_exacs_tem
     ADR_tem % kfl_length             =  kfl_ellen_tem
     ADR_tem % kfl_discretization     =  kfl_discr_tem
     if( kfl_sgsac_tem /= 1 ) then
        ADR_tem % kfl_first_order_sgs = 0
     else
        ADR_tem % kfl_first_order_sgs = 1
     end if
     ADR_tem % number_euler_steps     =  neule_tem

     ADR_tem % lun_output4            =  int(momod(modul) % lun_outpu,4)
     ADR_tem % bemol                  =  bemol_tem
     ADR_tem % tau_parameters(1:3)    =  staco_tem(1:3)
     ADR_tem % shock                  =  shock_tem  
     ADR_tem % kfl_skewsymm           =  kfl_skews_tem

  end select

end subroutine tem_inivar
