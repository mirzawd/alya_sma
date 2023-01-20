!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Turbul
!> @{
!> @file    tur_reaphy.f90
!> @author  Guillaume
!> @brief   Read physical data
!> @details Read physical data
!> @} 
!-----------------------------------------------------------------------
subroutine tur_reaphy()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_reaphy
  ! NAME 
  !    tur_reaphy
  ! DESCRIPTION
  !    This routine reads the physical treatment 
  ! USED BY
  !    tur_turnon
  !***
  !-----------------------------------------------------------------------
  !.md<module>turbul
  !.md<input>case.tur.dat
  !.md<pos>0
  !.md<sec>
  use def_parame
  use def_inpout
  use def_master
  use def_turbul
  use def_domain
  use mod_memchk
  use def_kermod,   only : cmu_st, kfl_kxmod_ker, kfl_logva
  use mod_ecoute,   only : ecoute
  use mod_messages, only : messages_live
  implicit none
  integer(ip) :: imate
  real(rp)    :: dummr
  logical     :: rwind ! realizable using wind coeffficients
 
  if( INOTSLAVE ) then
     !
     ! Initializations (defaults)
     !
     kfl_model_tur = 0                                   ! No model
     kfl_timei_tur = 0                                   ! Stationary flow
     kfl_cotem_tur = 0                                   ! No temperature coupling
     kfl_advec_tur = 1                                   ! Advection
     kfl_colev_tur = 0                                   ! No coupling with LEVELS
     kfl_ddesm_tur = 0                                   ! No model DDES 
     kfl_sasim_tur = 0                                   ! No SAS model
     kfl_rotat_tur = 0                                   ! No rotating reference frame
     kfl_inifi_tur(1) = 0                                ! Initial fields
     kfl_inifi_tur(2) = 0                                ! Initial fields
     kfl_inifi_tur(3) = 0                                ! Initial fields
     inits_tur     = 0                                   ! Initial time step
     nturb_tur     = 1                                   ! # turbulence variables
     ipara_tur     = 0                                   ! Integer Parameters
     lawde_tur     = 0                                   ! Law for rho
     lawvi_tur     = 0                                   ! Law for mu 
     kfl_kxmod_tur     = 0                                   ! k eps modification
     kfl_logva_tur     = 0                                   ! works with logarithm of unknowns
     boube_tur     = 0.0_rp                              ! Beta
     grnor_tur     = 0.0_rp                              ! Gravity norm |g|
     gravi_tur     = 0.0_rp                              ! Gravity vector g
     prtur_tur     = 1.0_rp                              ! Turbulent Prandtl number
     param_tur     = 0.0_rp                              ! Real parameters
     densi_tur     = 0.0_rp                              ! Density
     visco_tur     = 0.0_rp                              ! Viscosity
     densa_tur     = 1.0_rp                              ! Air density
     visca_tur     = 1.0_rp                              ! Air viscosity
     cddes_tur     = 0.65_rp                             ! CDDES - smagorinsky-like constant for DES
     inv_l_max          = 0.0_rp                         ! Mixing length (k-eps, CFDWind2)
     kfl_discd_tur      = 0
     ldiss_material_tur = 0
     !
     ! Local variables
     !
   
     rwind         =.false.                              ! k eps realizable model using wind coefficients
     !
     ! Reach the section
     !
     rewind(lisda)
     call ecoute('tur_reaphy')
     do while(words(1)/='PHYSI')
        call ecoute('tur_reaphy')
     end do
     !
     ! Begin to read data
     ! 
     do while(words(1)/='ENDPH')
        call ecoute('tur_reaphy')

        if(words(1)=='PROBL') then
           !
           ! Problem definition data
           !
           !.md# Physical Properties Definition
           !.md<code>
           !.md<0>PHYSICAL_PROBLEM
           !
           call ecoute('tur_reaphy')
           do while(words(1)/='ENDPR')
              !
              !.md<1>PROBLEM_DEFINITION
              ! 
              if(words(1)=='MODEL') then
                 !.md<2>MODEL  SPALA |  STDKE | KOMEGA | SSTKO   RNG | REALI (WIND) | WIL88 | WIL98 | WIL06 (LowRe)| TKESGS  
                 !.md<field>MODEL
                 !.md<com>Selection of RANS model that will be used for calculation. The third word is the variation of model kind, for example KOMEGA WIL98 or STD_KE REALIZALE.
                 if(words(2)=='SPALA') then
                    kfl_model_tur =  1
                 else if (words(2)=='KXUCH') then
                    kfl_model_tur =  2
                 else if (words(2)=='STDKE') then
                    kfl_model_tur = 11
                    if (exists('RNG  ') ) kfl_kxmod_tur = 1
                    if (exists('REALI') ) then  ! realizable model
                       kfl_kxmod_tur = 2
                       if (exists('WIND ') )  rwind=.true.
                       param_tur(1)= 1.0_rp       !sigma_k
                       param_tur(2)= 1.2_rp       !sigma_e
                       param_tur(3)= 0.0_rp       !c1
                       param_tur(4)= 1.90_rp      !c2
                       param_tur(5)= 0.0_rp       !c3
                       param_tur(6)= 0.09_rp      !cmu
                       cmu_st      = param_tur(6) !copy only in master
                       if (rwind) then ! realizable for wind
                          param_tur(6)= 0.0333_rp      !cmu
                          cmu_st      = param_tur(6)   !copy (only in master, later in all subdomains (tur_sendat))
                          ! sigma_e to satisfy wall equilibrium
                          param_tur(2)= 0.41_rp*0.41_rp/( param_tur(4)*sqrt(cmu_st)-0.43_rp &
                               *sqrt(cmu_st/0.09_rp))    
                       end if
                       
                    end if
                    if (exists('KEFP ') ) kfl_kxmod_tur = 3 !  An improved k - ε model applied to a wind turbine wake in
                    !              atmospheric turbulence, Wind Energ. 2013; 00:1–19 DOI: 10.1002/we
                    if (exists('LOGVA') ) kfl_logva_tur = 1 ! works with logarithmic functions 
                    ! copy variable to kernel 
                    if (kfl_logva_tur==1) kfl_logva=1
                    if (kfl_logva==1)     kfl_logva_tur=1
                    kfl_kxmod_ker = kfl_kxmod_tur  ! copy variable to kernel 
                 else if (words(2)=='LAUND') then
                    kfl_model_tur = 12
                 else if (words(2)=='CHIEN') then
                    kfl_model_tur = 13
                 else if (words(2)=='LAMBR') then
                    kfl_model_tur = 14
                 else if (words(2)=='RODIT') then
                    kfl_model_tur = 15
                 else if (words(2)=='NATXU') then 
                    kfl_model_tur = 16
                 else if (words(2)=='KEPSV') then
                    kfl_model_tur = 17
                 else if (words(2)=='KEPSJ'.or.words(2)=='JAWHW') then
                    kfl_model_tur = 18
                 else if (words(2)=='KEPSN'.or.words(2)=='NAGAN') then
                    kfl_model_tur = 19
                 else if (words(2)=='KEPSP') then
                    kfl_model_tur = 20
                 else if (words(2)=='KOMEG') then
                    kfl_model_tur = 30
                    if( exists('ORIGI') ) ipara_tur(1) = 0 ! low Reynolds 98
                    if( exists('WIL88') ) then
                       ipara_tur(1) = 1_ip ! Wilcox 1988
                       param_tur(1) = 2.0_rp ! sigma_k
                       param_tur(2) = 2.0_rp ! sigma_w
                       call messages_live('RANS MODEL:  K-OMEGA WILCOX88')
                    else if( exists('WIL98') ) then
                       ipara_tur(1) = 2_ip ! Wilcox 1998
                       param_tur(1) = 2.0_rp ! sigma_k
                       param_tur(2) = 2.0_rp ! sigma_w
                       call messages_live('RANS MODEL: K-OMEGA WILCOX98')
                    else if( exists('WIL06') ) then 
                       ipara_tur(1) = 3_ip ! Wilcox 2006
                       param_tur(1) = 5.0_rp/3.0_rp ! sigma_k
                       param_tur(2) = 2.0_rp        ! sigma_w
                       call messages_live('RANS MODEL: K-OMEGA WILCOX06')
                       if (exists('LOWRE')) then
                          ipara_tur(1) = 4_ip ! low Reynolds number version of WIL06
                          call messages_live('RANS MODEL: Low Reynolds approximation')
                       end if
                    end if
                    ! copy flag variable (rans model kind) to kernel module
                    kfl_kxmod_tur = ipara_tur(1)
                    kfl_kxmod_ker = kfl_kxmod_tur  
                    
                 else if (words(2)=='BREDB') then
                    kfl_model_tur = 31
                 else if (words(2)=='SSTKO') then
                    kfl_model_tur = 32
                    if(exists('ROTAT'))  kfl_rotat_tur = 1! Rotating Reference frame
                 else if (words(2)=='TKE-S'.or.words(2)=='TKESG') then   ! LES model, based on ptognostic TKE equation
                    kfl_model_tur = 4   
                 end if
                 call tur_inivar(2_ip)
                 if(exists('START')) then
                    inits_tur=getint('START',0_ip,'#Start at time step')
                 end if

              else if(words(1)=='TEMPO') then                 ! Temporal evolution
                 !.md<2>TEMPORAL_DERIVATIVES: ON | OFF                                                !.md<field>TEMPORAL_DERIVATIVES
                 !.md<com>Decide if time derivative should be calculated in the RANS equations.
                 !.md<com>This term should be put to ON as for most of the problems, the convergence
                 !.md<com>is enhanced when passing through a stationary regime to reach a steady state.
                 if(exists('ON   ')) kfl_timei_tur = 1 
                 
              else if(words(1)=='DDES ') then                 ! DDES model flag
                 !.md<2>DDES:  IMPRO                                             
                 kfl_ddesm_tur = 1
                 if(exists('IMPRO ')) kfl_ddesm_tur = 2
                 cddes_tur     = getrea('CDDES ',0.65_rp,'#CDDES Smagorinski Constant')

              else if(words(1)=='SAS') then                 ! SAS model flag
                 if(exists('ON   ')) kfl_sasim_tur = 1

              else if(words(1)=='CONVE') then                 ! Advection
                 if(exists('NASTI')) kfl_advec_tur = 1 
                 if(exists('ON   ')) kfl_advec_tur = 1 
                 ! OFF option is not properly coded
!                 if(exists('OFF  ')) kfl_advec_tur = 0 

              else if(words(1)=='TEMPE'.and.words(2)=='BOUSS') then
                 !.md<2>TEMPERATURE: BOUSSINESQ | OFF  SOGACHEV  BETA=real G=real  GX=real GY=real GZ=real 
                 !.md<field>TEMPERATURE
                 !.md<com>Boussinesq coupling can be selected. BETA is the thermal expansion coefficients, which is 1/Tref in air or ideal gas.
                 !.md<com>G is the gravity field, which generates production of turbulent kinetic energy when density variations exist.
                 !.md<com>GX,GY,GZ denotes the direction of then gravity field.
                 !.md<com>SOGACHEV model should be selected, used for K-E in atmosphere
                 kfl_cotem_tur = 1
                 boube_tur    = getrea('BETA ',0.0_rp,'#Beta coefficient')
                 grnor_tur    = getrea('G    ',0.0_rp,'#Gravity norm')
                 gravi_tur(1) = getrea('GX   ',0.0_rp,'#x-component of G')
                 gravi_tur(2) = getrea('GY   ',0.0_rp,'#y-component of G')
                 gravi_tur(3) = getrea('GZ   ',0.0_rp,'#z-component of G')
                 call vecuni(3_ip,gravi_tur,dummr)
                 if (exists('SOGAC')) kfl_cotem_tur = 2 ! Sogachev's model

              else if(words(1)=='LEVEL') then
                 !.md<2>LEVEL: ON | OFF THICK 
                 if(words(2)=='ON   ') then
                    kfl_colev_tur=1
                    if(exists('THICK')) then
                       call runend('TUR_REAPHY: now thickness is read by level and stored in a thick that belongs to master')  
                    end if
                    if(exists('STAGG')) kfl_colev_tur=3
                 end if

              else if(words(1)=='GRAVI') then
                 !.md<2>GRAVITY: NORM=real  GX=real GY=real GZ=real    
                 grnor_tur     = getrea('NORM ',0.0_rp,'#Gravity norm')
                 gravi_tur(1)  = getrea('GX   ',0.0_rp,'#x-component of g')
                 gravi_tur(2)  = getrea('GY   ',0.0_rp,'#y-component of g')
                 gravi_tur(3)  = getrea('GZ   ',0.0_rp,'#z-component of g')
                 call vecuni(3_ip,gravi_tur,dummr)
              else if(words(1)=='DISCD') then ! extra dissipation term for actuator disc model
                 !.md<2>DISCDISS:
                 !.md<3>MATER = real
                 !.md<2>ENDDISS:
                 !.md<field>DISCDISS
                 !.md<com>Extra dissipation term for actuator disc model
                 kfl_discd_tur= 1
               
                 call ecoute('tur_reaphy')
                 do while( words(1) /= 'ENDDI' )
                    if( words(1) == 'MATER' ) then
                       imate = getint('MATER',1_ip,'#Material force number')
                       if( imate < 1 .or. imate > nmate ) call runend('TUR_REAPHY: WRONG MATERIAL NUMBER')
                       ldiss_material_tur(imate) = 1 ! materials for dissipation
                    else 
                       call runend('TUR_REAPHY: WRONG DISSIPATION FIELD') 
                    end if
                    if (max_mater_tur.lt.nmate) call runend('TUR_REAPHY: too many materials for ldiss_material')
                    call ecoute('tur_reaphy')
                 end do
              end if
              call ecoute('tur_reaphy')
              !md<1>END_PROBLEM_DEFINITION
           end do

        else if(words(1)=='PROPE') then
           !
           ! Properties and model constants
           !
           !.md<1>PROPERTIES:
           !.md<field>PROPERTIES
           !.md<com>Properties. Here, some properties could be required 
           !.md<com>module. The fluid properties (density, viscosity) should be defined in Kermod.
           do while(words(1)/='ENDPR')
              call ecoute('tur_reaphy')

              if(words(1)=='DENSI') then                      ! Density (rho)
                 densi_tur(1:10)=param(1:10)

              else if(words(1)=='VISCO') then                 ! Viscosity (mu)
                 visco_tur(1:10)=param(1:10)

              else if( words(1) == 'LMAXI' ) then
                 !.md<2>LMAXI
                 !.md<field>LMAXI
                 !.md<com>Maximum mixing length, used in K-E for atmospheric flows
                 !.md<com>Maximum mixing length is calclulated in terms of latitude and geostrophic pressure gradient
                 ! max. mixing length l_max
                 !
                 inv_l_max = getrea('LMAXI',0.0_rp,'#Max. mixing length (CFDWind model)')
                 ! The inverse of l_max  will be stored
                 if( inv_l_max < 1.0e8_rp.and.abs(inv_l_max)>1.0e-8_rp) then
                    inv_l_max= 1.0_rp/inv_l_max  
                 else if (inv_l_max >1.0e7_rp) then
                    inv_l_max= 0.0_rp
                 else
                    call runend ('tur_reaphy:invalid input data for l_max')
                 end if

              else if(words(1)=='LAWDE') then                 ! Law for rho
                 if(words(2)=='CONST') then
                    lawde_tur=0
                 else if(words(2)=='DENSI') then
                    lawde_tur=1
                 end if

              else if(words(1)=='LAWVI') then                 ! Law for mu
                 if(words(2)=='CONST') then
                    lawvi_tur=0
                 else if(words(2)=='VISCO') then
                    lawvi_tur=1
                 end if

              else if(words(1)=='TURBU') then                 ! turbulent Prandtl number
                 !.md<2>TURBULENT_PRANDTL  real
                 !.md<field>TURBULENT_PRANDTL
                 !.md<com> "real" is the prandtl number used in for turbulent kinetic energy production
                 prtur_tur = getrea('TURBU',0.0_rp,'#Turbulent Prandtl number')

              else if(words(1)=='REALP') then                 ! Model real parameters
                 !.md<2>REAL_PARAMETERS  real real real real real real | AUTOM
                 !.md<field>REAL_PARAMETERS 
                 !.md<com>Parameters used in the RANS model
                 !.md<com>This field is not mandatory and not used for all models.
                 !.md<com>AUTOM selects the parameters automatically
                 if(words(2)/='AUTOM') then
                    param_tur(1:11)=param(1:11)
                    ! Special k-e models with defined constants
                    if (kfl_model_tur == 11 ) then ! K-E model
                       if (kfl_kxmod_tur==1) then !overwrite for k-e RNG model 
                          param_tur(1)= 0.7179_rp    !sigma_k
                          param_tur(2)= 0.7179_rp    !sigma_e
                          param_tur(3)= 1.42_rp      !c1
                          param_tur(4)= 1.68_rp      !c2
                          param_tur(5)= 0.0_rp       !c3
                          param_tur(6)= 0.085_rp     !cmu
                          cmu_st      = param_tur(6) !copy only in master
                       else if(kfl_kxmod_tur==2) then !overwrite for k-e Realizable model
                          param_tur(1)= 1.0_rp       !sigma_k
                          param_tur(2)= 1.2_rp       !sigma_e
                          param_tur(3)= 0.0_rp       !c1
                          param_tur(4)= 1.90_rp      !c2
                          param_tur(5)= 0.0_rp       !c3
                          param_tur(6)= 0.09_rp      !cmu
                          cmu_st      = param_tur(6) !copy only in master
                          if (rwind) then ! realizable for wind
                             param_tur(6)= 0.0333_rp      !cmu
                             cmu_st      = param_tur(6)   !copy (only in master, later in all subdomains (tur_sendat))
                             ! sigma_e to satisfy wall equilibrium
                             param_tur(2)= 0.41_rp*0.41_rp/( param_tur(4)*sqrt(cmu_st)-0.43_rp &
                                  *sqrt(cmu_st/0.09_rp))    
                          end if
                       
                       end if
                    end if
                 end if

              else if(words(1)=='INTEG') then                 ! Model integer parameters
                 !.md<2>INTEGER_PARAMETERS  int int int int int | AUTOM
                 !.md<field>INTEGER_PARAMETERS 
                 !.md<com>Parameters used in the RANS model. Not mandatory.
                 !.md<com>AUTOM selects the parameters automatically.
                 if(words(2)/='AUTOM') then
                    ipara_tur(1:2)=int(param(1:2))
                 end if

              else if(words(1)=='AIRDE') then
                 densa_tur=getrea('AIRDE',1.0_rp,'#Air density')

              else if(words(1)=='AIRVI') then
                 visca_tur=getrea('AIRVI',1.0_rp,'#Air viscosity')

              else if(words(1)=='INITI') then
                 !
                 ! Initial Fields
                 !
                 if(words(2)=='CODES') then ! Initial Fields given by codes
                    call runend('TUR_REAPHY: "CODES" OPTION IN INITIAL_CONDITION IS NOT AVAILABLE')
                 else ! Initial Fields given on nodes
                    call ecoute('tur_reaphy')
                    do while(words(1)/='ENDIN')
                       if (words(1) == 'KEY  ') then
                          kfl_inifi_tur(1) = 1
                          nfiel_tur(1) = -getint('FIELD',1_ip,'#Field Number for kinetic energy')
                       else if (words(1) == 'NUFOR') then
                          kfl_inifi_tur(1) = 1
                          kfl_inifi_tur(2) = 1
                          nfiel_tur(1) = -getint('FIELD',1_ip,'#Field Number for kinetic eddy viscosity')
                       else if (words(1) == 'EPSIL') then
                          kfl_inifi_tur(2) = 1
                          nfiel_tur(2) = -getint('FIELD',1_ip,'#Field Number for second variable (epsilon)')
                       else if (words(1) == 'OMEGA') then
                          kfl_inifi_tur(2) = 2
                          nfiel_tur(2) = -getint('FIELD',1_ip,'#Field Number for second variable (omega)')
                       else if (words(1) == 'THIRD') then ! in case of third variable, not programmed yet
                          kfl_inifi_tur(3) = 1
                          nfiel_tur(3) = -getint('FIELD',1_ip,'#Field Number for third variable')
                       end if
                       call ecoute('tur_reaphy')
                    end do
                    if (kfl_inifi_tur(1) > 0 .and. kfl_inifi_tur(2) == 0) then
                       ! hay que ver aun...
                       call runend('TUR_REAPHY: KEY AND SECOND VARIABLE (EPSILON OR OMEGA) MUST BE GIVEN')
                    else if ((kfl_inifi_tur(1) + kfl_inifi_tur(2) + kfl_inifi_tur(3)) == 0) then
                       call runend('TUR_REAPHY: NO FIELDS WERE GIVEN')
                    end if
                 end if

              end if
              !.md<1>END_PROPERTIES
           end do
           !.md<0>END_PHYSICAL_PROPERTIES
        end if
     end do
     !
     ! Define number of turbulence variables
     !
     if(TUR_ONE_EQUATION_MODEL) then
        nturb_tur=1                                           ! One-equation models
     else if(TUR_K_EPS_V2_F) then
        nturb_tur=4                                           ! k-eps-v2-f
     else if(TUR_K_EPS_PHI_F) then
        nturb_tur=4                                           ! k-eps-phi-f
     else
        nturb_tur=2                                           ! Two-equation models
     end if
    
  end if

end subroutine tur_reaphy
