!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_reaphy.f90
!> @author  Guillaume Houzeaux
!> @brief   Read physical data
!> @details Read physical data
!> @} 
!-----------------------------------------------------------------------
subroutine nsi_reaphy()
!-----------------------------------------------------------------------
 !.md<module>nastin
 !.md<input>case.nsi.dat
 !.md<pos>0
 !.md<sec>

  use def_parame
  use def_inpout
  use def_master
  use def_nastin
  use def_domain
  use def_kermod, only : gasco, kfl_prope, anipo_ker
  use mod_ker_space_time_function, only : space_time_function_number
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip) :: imate, icoun, ii
  real(rp)    :: dummr

  if( INOTSLAVE ) then
     !
     ! Initializations (defaults)
     !
     kfl_timei_nsi = 0                                         ! Time integration off
     kfl_advec_nsi = 0                                         ! Convection is off
     kfl_convection_type_nsi = NSI_CONVECTION_NON_CONSERVATIVE ! Type of convection
     kfl_confi_nsi = 0                                         ! Flow is not confined -- This should not be initialized here it is initialized in nsi_reabcs
     kfl_visco_nsi = 1                                         ! Viscous term on
     kfl_colev_nsi = 0                                         ! No coupling with LEVELS
     kfl_regim_nsi = 0                                         ! Incompressible flow
     kfl_dynco_nsi = 0                                         ! Dynamical coupling
     kfl_prthe_nsi = 0                                         ! Thermodynamic pressure calculation
     kfl_surte_nsi = 0                                         ! Do not include surface tension
     kfl_force_nsi = 0                                         ! Force type    
     kfl_vegeo_time_nsi = 0                                    ! Time dependent force function
     kfl_mfrco_nsi = 0                                         ! Mass flow rate control activation flag (=0 OFF, =1 ON)
     kfl_bnods_nsi = 0                                         ! boundary nodes defined
     kfl_hydro_gravity_nsi = 0                                 ! No hydrostatic gravity
     kfl_anipo_nsi = 0                                         ! No anisotropic porosity
     
     fcons_nsi     = 0.0_rp                                    ! Non conservative form default
     fvins_nsi     = 1.0_rp                                    ! Divergence form default
     dtinv_nsi     = 1.0_rp                                    ! 1/dt=1.0
     nbnod_nsi     = 0                                         ! Number of nodes on boundary
     nbtim_nod_nsi = 0                                         ! total number of time instances and nodes for time-space boundary from file

     grnor_nsi     = 0.0_rp                                    ! Gravity norm
     gravi_nsi     = 0.0_rp                                    ! Gravity vector
     gravb_nsi     = 0.0_rp                                    ! Gravity vector for Boussinesq coupling
     fvnoa_nsi     = 0.0_rp                                    ! Frame angular velocity norm     
     fvdia_nsi     = 0.0_rp                                    ! Frame angular velocity direction 
     fvela_nsi     = 0.0_rp                                    ! Frame angular velocity vector    
     fvpaa_nsi     = 0.0_rp                                    ! Frame angular velocity parameters  
     kfl_fvfua_nsi = 0                                         ! Frame angular velocity function  
     fanoa_nsi     = 0.0_rp                                    ! Frame angular acceleration norm   
     fadia_nsi     = 0.0_rp                                    ! Frame angular acceleration direction 
     facca_nsi     = 0.0_rp                                    ! Frame angular acceleration vector       
     fvnol_nsi     = 0.0_rp                                    ! Frame linear velocity norm 
     fvdil_nsi     = 0.0_rp                                    ! Frame linear velocity direction     
     fvell_nsi     = 0.0_rp                                    ! Frame linear velocity vector        
     fvpal_nsi     = 0.0_rp                                    ! Frame linear velocity parameters  
     kfl_fvful_nsi = 0                                         ! Frame linear velocity function  
     fanol_nsi     = 0.0_rp                                    ! Frame linear acceleration norm      
     fadil_nsi     = 0.0_rp                                    ! Frame linear acceleration direction 
     faccl_nsi     = 0.0_rp                                    ! Frame linear acceleration vector 
     frotc_nsi     = 1.0e20_rp                                 ! Rotation center
     centr_nsi     = 1.0_rp                                    ! Centrifugal force

     kfl_cotem_nsi = 0                                         ! No coupling with temperature
     boube_nsi     = 0.0_rp                                    ! Beta
     bougr_nsi     = 0.0_rp                                    ! g
     boutr_nsi     = 0.0_rp                                    ! Tr
     lowtr_nsi    =  0.0_rp                                    ! Ref tempe for hydrodynamic pressure in Low Mach regime

     kfl_grtur_nsi = 0                                         ! Do not add grad(K)
     turbu_nsi     = 0.0_rp                                    ! Turbulence parameter
     heihy_nsi     = 0.0_rp                                    ! Height for hydrostatic pressure
     if (kfl_prope ==0 ) gasco     = 287.0_rp                  ! Gas constant R [J/K Kg], only if kernel does not carry it
     sphea_nsi     = 1006.0_rp                                 ! Specific heat Cp [J/K Kg]
     prthe_nsi     = 101325.0_rp                               ! Thermodynamic pressure
     imate         = 1

     surte_nsi = 0.0_rp                                        ! Surface tension coeficient (sigma)
     tmass_nsi = 0.0_rp                                        ! initial mean density for low Mach problem
     mfrub_nsi = 0.0_rp                                        ! Target bulk velocity when mass flow control activated
     ubpre_nsi = 0.0_rp                                        ! Bulk velocity from previous time-steo
     mfrse_nsi = 1                                             ! Set from which the mass flow rate is calculated (by default = 1)
     mfccf_nsi = 1.0_rp                                        ! Coefficient for the mass flow control formula (by default = 1)

     ntabf_nsi = 0_ip                                           ! Number of tables for spatial and temporal function tables


     nfiel_nsi(1:2) = 0                                   ! Fields assignement
     !
     ! Reach the section
     !
     rewind(lisda)
     call ecoute('nsi_reaphy')
     do while(words(1)/='PHYSI')
        call ecoute('nsi_reaphy')
     end do
     !--><group>
     !-->    <groupName>PHYSICAL_PROBLEM</groupName>
     !
     !.md# Physical Properties Definition
     !.md<code>
     !.md<0>PHYSICAL_PROBLEM
     !
     do while(words(1)/='ENDPH')

        call ecoute('nsi_reaphy')

        if(words(1) == 'PROBL' ) then

           call ecoute('nsi_reaphy')

           do while(words(1)/='ENDPR')
              !--><subGroup>
              !-->     <subGroupName>PROBLEM_DEFINITION</subGroupName>           
              !
              !.md<1>PROBLEM_DEFINITION
              ! 
              if(words(1)=='TEMPO') then          
                 !--><inputLine>
                 !-->      <inputLineName>TEMPORAL_DERIVATIVES</inputLineName>
                 !-->      <inputLineHelp>
                 !-->	   <![CDATA[TEMPORAL_DERIVATIVES:			     
                 !-->         Decide if time derivative should be calculated in the Navier-Stokes equations.
                 !-->         This term should be put to ON as for most of the problems, the convergence
                 !-->          is enhanced when passing through a stationary regime to reach a steady state.
                 !-->        ]]>
                 !-->		</inputLineHelp>
                 !--></inputLine> 
                 !
                 !.md<2>TEMPORAL_DERIVATIVES: ON | OFF                                                              $ Existence of temporal derivatives
                 !.md<field>TEMPORAL_DERIVATIVES
                 !.md<com>Decide if time derivative should be calculated in the Navier-Stokes equations.
                 !.md<com>This term should be put to ON as for most of the problems, the convergence
                 !.md<com>is enhanced when passing through a stationary regime to reach a steady state..$`x^2 \frac{1}{2}`$...
                 !
                 if(exists('ON   ')) kfl_timei_nsi = 1

              else if( words(1) == 'HYDRO' ) then 
                 !
                 ! Hydrostatic gravity term in the RHS of NS equations
                 !
                 if(exists('ON   ')) kfl_hydro_gravity_nsi = 1

              else if(words(1)=='ANISO') then
                 !
                 ! Anisotropic porosity on
                 !
                 if( option('ANISO') ) kfl_anipo_nsi = 1
                 
              else if(words(1)=='CONVE') then 
                 !--><inputLine>
                 !-->     <inputLineName>CONVECTIVE_TERM</inputLineName>
                 !-->     <inputLineHelp><![CDATA[Decide if the convective term should be calculated.]]></inputLineHelp>
                 !--> </inputLine>                         
                 !
                 !.md<2>CONVECTIVE_TERM:      ON/EMAC/EMAC2/NON_CONSERVATIVE | OFF        $ Existence of convective term    
                 !.md<field>CONVECTIVE_TERM
                 !.md<com>Decide if the convective term should be calculated.
                 !
                 if( words(2) /= 'OFF  ') kfl_advec_nsi = 1
                 if( words(2) == 'NONCO' ) then
                    fcons_nsi = 0.0_rp
                    kfl_convection_type_nsi = NSI_CONVECTION_NON_CONSERVATIVE
                 else if( words(2) == 'CONSE' ) then
                    fcons_nsi = 1.0_rp
                    kfl_convection_type_nsi = NSI_CONVECTION_CONSERVATIVE
                 else if(words(2) == 'SKEWS' ) then                 
                    fcons_nsi = 0.5_rp
                    kfl_convection_type_nsi = NSI_CONVECTION_SKEW
                 else if(words(2) == 'EMAC ' ) then                 
                    fcons_nsi = -1.0_rp
                    kfl_convection_type_nsi = NSI_CONVECTION_EMAC
                 else if(words(2) == 'EMAC2' ) then                 
                    fcons_nsi = -1.0_rp
                    kfl_convection_type_nsi = NSI_CONVECTION_EMAC2
                 end if
                 if( words(3) == 'VELOC' ) then
                    if( words(4) == 'NAVIE' ) then
                       kfl_advec_nsi = 1
                    else
                       kfl_advec_nsi = getint('VELOC',0_ip,'#Velocity function')
                    end if
                 end if

              else if(words(1)=='VISCO') then 
                 !--> <inputLine>
                 !-->      <inputLineName>VISCOUS_TERM</inputLineName>
                 !-->      <inputLineHelp><![CDATA[Form of the viscous term. Different forms are available dependending on the problem treated.
                 !-->                             
                 !-->                             -  LAPLACIAN: this form is the simplest one and it requires less operations than the others. It is however an approximation
                 !-->                              for low Mach regime or fluid with variable viscosity.
                 !-->                             -  DIVERGENCE: this form is the correct one for flows with variable viscosity 
                 !-->                             -  COMPLETE: this form should be used for low Mach regime as it includes the non-zero velocity divergence term. 
                 !-->                      ]]>
                 !-->      </inputLineHelp>
                 !-->      <inputElement>
                 !-->          <inputElementType>combo</inputElementType>
                 !-->          <item>
                 !-->              <itemName>OFF</itemName>
                 !-->          </item>
                 !-->          <item>
                 !-->              <itemName>LAPLACIAN</itemName>
                 !-->          </item>
                 !-->          <item>
                 !-->              <itemName>DIVERGENCE</itemName>
                 !-->          </item>
                 !-->          <item>
                 !-->              <itemName>COMPLETE</itemName>
                 !-->          </item>
                 !-->      </inputElement>
                 !-->  </inputLine>
                 !
                 !.md<2>VISCOUS_TERM:         OFF | LAPLACIAN | DIVERGENCE | COMPLETE                               $ Viscous term 
                 !.md<field>VISCOUS_TERM
                 !.md<com>Form of the viscous term. Different forms are available dependending on the problem treated.
                 !.md<com>
                 !.md<com>    -  LAPLACIAN: this form is the simplest one and it requires less operations than the others. It is however an approximation
                 !.md<com>     for low Mach regime or fluid with variable viscosity. 
                 !.md<com>    -  DIVERGENCE: this form is the correct one for flows with variable viscosity. 
                 !.md<com>    -  COMPLETE: this form should be used for low Mach regime as it includes the non-zero velocity divergence term. 
                 !.md<com> 
                 !
                 if(words(2)=='OFF  ') kfl_visco_nsi = 0
                 if(words(2)=='ON   ') fvins_nsi     = 0.0_rp                        ! Divergence form default
                 if(words(2)=='DIVER') fvins_nsi     = 1.0_rp                        ! Idem
                 if(words(2)=='LAPLA') fvins_nsi     = 0.0_rp                        ! Laplacian form
                 if(words(2)=='COMPL' .or. words(2)=='FULL ') fvins_nsi = 2.0_rp     ! Complete form

              else if( words(1) == 'REGIM' ) then
                 !--><inputLine>
                 !-->    <inputLineName>REGIME</inputLineName>
                 !-->    <inputLineHelp>Regime of the flow. If the flow is low mach, additional information is required. If system is open</inputLineHelp>
                 !-->    <inputElement>
                 !-->        <inputElementType>combo</inputElementType>
                 !-->        <item>
                 !-->            <itemName>INCOMPRESSIBLE</itemName>                            
                 !-->        </item>
                 !-->        <item>
                 !-->            <itemName>LOW_MACH</itemName>
                 !-->            <itemDependence>0</itemDependence>
                 !-->        </item>
                 !-->    </inputElement>
                 !-->    <inputElement>
                 !-->        <inputElementGroup>0</inputElementGroup>
                 !-->        <inputElementType>edit</inputElementType>
                 !-->        <inputElementValueType>REAL</inputElementValueType>                 
                 !-->        <inputLineEditName>Temperature</inputLineEditName>
                 !-->        <inputLineEditValue>5</inputLineEditValue>
                 !-->    </inputElement>
                 !-->    <inputElement>
                 !-->        <inputElementGroup>0</inputElementGroup>
                 !-->        <inputElementType>edit</inputElementType>
                 !-->        <inputElementValueType>REAL</inputElementValueType>                 
                 !-->        <inputLineEditName>Pressure</inputLineEditName>
                 !-->        <inputLineEditValue>5</inputLineEditValue>
                 !-->    </inputElement>
                 !--></inputLine>
                 !
                 !.md<2>REGIME:               INCOMPRESSIBLE | LOW_MACH [TEMPERATURE= real or PRESSURE= real]       $ Flow regime. If low mach, reference temp is needed
                 !.md<field>REGIME
                 !.md<com>Regime of the flow. If the flow is low mach, additional information is required.
                 !.md<com>If system is open: 
                 !
                 if(words(2)=='INCOM') then                 ! Incompressible
                    kfl_regim_nsi=0
                 else if(words(2)=='COMPR') then            ! Compressible
                    kfl_regim_nsi=1
                    if(exists('PRESS')) kfl_regim_nsi=1
                    if(exists('DENSI')) kfl_regim_nsi=2
                 else if(words(2)=='LOWMA') then            ! Low-Mach
                    kfl_regim_nsi=3
                    if(exists('TEMPE')) then
                       lowtr_nsi = getrea('TEMPE',0.0_rp,'#Reference temperature for hydrodynamic pressure')
                    endif
                 end if

              else if(words(1)=='DYNAM') then               ! Dynamical coupling
                 if(words(2)=='ON   ') kfl_dynco_nsi=1

              else if(words(1)=='GRAVI') then
                 !--><inputLine>
                 !-->   <inputLineName>GRAVITY</inputLineName>
                 !-->   <inputLineHelp> If this option is present, the gravity acceleration is added to the Navier-Stokes equations. Nastin does not use Kermod's gravity,
                 !-->                   as some modules could require gravity (Partis) but Nastin could run without.</inputLineHelp>
                 !-->  <inputElement>
                 !-->      <inputElementType>edit</inputElementType>
                 !-->          <inputElementValueType>REAL</inputElementValueType>                 
                 !-->      <inputLineEditName>Norm</inputLineEditName>
                 !-->       <inputLineEditValue>5</inputLineEditValue>
                 !-->   </inputElement>
                 !-->   <inputElement>
                 !-->       <inputElementType>edit</inputElementType>
                 !-->          <inputElementValueType>REAL</inputElementValueType>                 
                 !-->       <inputLineEditName>GX</inputLineEditName>
                 !-->       <inputLineEditValue>5</inputLineEditValue>
                 !-->   </inputElement>
                 !-->   <inputElement>
                 !-->       <inputElementType>edit</inputElementType>
                 !-->          <inputElementValueType>REAL</inputElementValueType>                 
                 !-->       <inputLineEditName>GY</inputLineEditName>
                 !-->       <inputLineEditValue>5</inputLineEditValue>
                 !-->   </inputElement>
                 !-->   <inputElement>
                 !-->       <inputElementType>edit</inputElementType>
                 !-->          <inputElementValueType>REAL</inputElementValueType>                 
                 !-->       <inputLineEditName>GZ</inputLineEditName>
                 !-->       <inputLineEditValue>5</inputLineEditValue>
                 !-->   </inputElement>
                 !--> </inputLine>
                 !
                 !.md<2>GRAVITY:              NORM= real, GX= real, GY= real, GZ= real                              $ Gravity accelaration
                 !.md<field>GRAVITY
                 !.md<com>If this option is present, the gravity acceleration is added to the 
                 !.md<com>Navier-Stokes equations. Nastin does not use Kermod's gravity, as some
                 !.md<com>modules could require gravity (Partis) but Nastin could run without.
                 !
                 grnor_nsi    = getrea('NORM ',0.0_rp,'#Gravity norm')
                 gravi_nsi(1) = getrea('GX   ',0.0_rp,'#x-component of g')
                 gravi_nsi(2) = getrea('GY   ',0.0_rp,'#y-component of g')
                 gravi_nsi(3) = getrea('GZ   ',0.0_rp,'#z-component of g')
                 call vecuni(3_ip,gravi_nsi,dummr)

              else if( words(1) == 'CENTR' ) then
                 !--><inputLine>
                 !-->    <inputLineName>CENTRIFUGAL_FORCE</inputLineName>
                 !-->    <inputLineHelp>GRAVITY: If this option is present, the gravity acceleration is added to the</inputLineHelp>
                 !--></inputLine>
                 !
                 !.md<2>CENTRIFUGAL_FORCE:    ON | OFF                                                              $ Centrifugal force
                 !.md<field>GRAVITY
                 !.md<com>If this option is present, the gravity acceleration is added to the 

                 if( words(2) == 'ON   ' ) then
                    centr_nsi = 1.0_rp
                 else if( words(2) == 'OFF  ' ) then
                    centr_nsi = 0.0_rp
                 end if

              else if(words(1)=='AXESR') then
                 !--><inputLine>
                 !-->    <inputLineName>Axes_rotation</inputLineName>
                 !-->    <inputLineHelp> Axes rotation vector definition.</inputLineHelp>
                 !-->    <inputElement>
                 !-->        <inputElementType>edit</inputElementType>
                 !-->          <inputElementValueType>REAL</inputElementValueType>                 
                 !-->        <inputLineEditName>Norm</inputLineEditName>
                 !-->        <inputLineEditValue>5</inputLineEditValue>
                 !-->    </inputElement>
                 !-->    <inputElement>
                 !-->        <inputElementType>edit</inputElementType>
                 !-->          <inputElementValueType>REAL</inputElementValueType>
                 !-->        <inputLineEditName>OX</inputLineEditName>
                 !-->        <inputLineEditValue>5</inputLineEditValue>
                 !-->    </inputElement>
                 !-->    <inputElement>
                 !-->        <inputElementType>edit</inputElementType>
                 !-->          <inputElementValueType>REAL</inputElementValueType>                 
                 !-->        <inputLineEditName>OY</inputLineEditName>
                 !-->        <inputLineEditValue>5</inputLineEditValue>
                 !-->    </inputElement>
                 !-->    <inputElement>
                 !-->        <inputElementType>edit</inputElementType>
                 !-->          <inputElementValueType>REAL</inputElementValueType>                 
                 !-->        <inputLineEditName>OZ</inputLineEditName>
                 !-->        <inputLineEditValue>5</inputLineEditValue>
                 !-->    </inputElement>
                 !--></inputLine>
                 !
                 !.md<2>AXES_ROTATION:        NORM= real, OX= real, OY= real, OZ= real                              $ Axes rotation vector definition
                 !.md<field>AXES_ROTATION
                 !.md<com>Axes rotation vector definition. 
                 !
                 fvnoa_nsi    = getrea('NORM ',0.0_rp,'#Axes rotation norm')
                 fvdia_nsi(1) = getrea('OX   ',0.0_rp,'#x-component of w')
                 fvdia_nsi(2) = getrea('OY   ',0.0_rp,'#y-component of w')
                 fvdia_nsi(3) = getrea('OZ   ',0.0_rp,'#z-component of w')
                 call vecuni(3_ip,fvdia_nsi,dummr)

              else if(words(1)=='AXESV') then
                 !--><inputLine>
                 !-->      <inputLineName>AXES_VELOCITY</inputLineName>
                 !-->      <inputLineHelp>Axes velocity vector definition.</inputLineHelp>
                 !-->      <inputElement>
                 !-->          <inputElementType>edit</inputElementType>
                 !-->          <inputElementValueType>REAL</inputElementValueType>
                 !-->          <inputLineEditName>Norm</inputLineEditName>
                 !-->          <inputLineEditValue>5</inputLineEditValue>
                 !-->      </inputElement>
                 !-->      <inputElement>
                 !-->          <inputElementType>edit</inputElementType>
                 !-->          <inputElementValueType>REAL</inputElementValueType>                 
                 !-->          <inputLineEditName>VX</inputLineEditName>
                 !-->          <inputLineEditValue>5</inputLineEditValue>
                 !-->      </inputElement>
                 !-->      <inputElement>
                 !-->          <inputElementType>edit</inputElementType>
                 !-->          <inputElementValueType>REAL</inputElementValueType>                 
                 !-->          <inputLineEditName>VY</inputLineEditName>
                 !-->          <inputLineEditValue>5</inputLineEditValue>
                 !-->      </inputElement>
                 !-->      <inputElement>
                 !-->          <inputElementType>edit</inputElementType>
                 !-->          <inputElementValueType>REAL</inputElementValueType>                 
                 !-->          <inputLineEditName>VZ</inputLineEditName>
                 !-->          <inputLineEditValue>5</inputLineEditValue>
                 !-->      </inputElement>
                 !-->  </inputLine>
                 !
                 !.md<2>AXES_VELOCITY:        NORM= real, VX= real, VY= real, VZ= real                              $ Axes velocity vector definition
                 !.md<field>AXES_VELOCITY
                 !.md<com>Axes velocity vector definition.
                 !
                 fvnol_nsi    = getrea('NORM ',0.0_rp,'#Axes velocity norm')
                 fvdil_nsi(1) = getrea('VX   ',0.0_rp,'#x-component of v')
                 fvdil_nsi(2) = getrea('VY   ',0.0_rp,'#y-component of v')
                 fvdil_nsi(3) = getrea('VZ   ',0.0_rp,'#z-component of v')
                 call vecuni(3_ip,fvdil_nsi,dummr)

              else if(words(1)=='ROTAT') then
                 !--><inputLine>
                 !-->    <inputLineName>ROTATION_FUNCTION</inputLineName>
                 !-->    <inputLineHelp><![CDATA[ROTATION_FUNCTION:
                 !-->     Rotation time function definition. It multiplies the AXES_ROTATION vector.
                 !-->     Parameters (up to 6 numbers after the type of function)
                 !-->     CONSTANT (no parameters, it implies a value of 1.0)
                 !-->     PARABOLIC (a t^2 + b t + c)    start_time end_time a b c 
                 !-->     SINUS (a cos( w t + phi) + c)  start_time end_time a w phi c]]></inputLineHelp>
                 !-->    <inputElement>
                 !-->        <inputElementType>combo</inputElementType>
                 !-->        <item>
                 !-->            <itemName>CONSTANT</itemName>
                 !-->        </item>
                 !-->        <item>
                 !-->            <itemName>PARABOLIC</itemName>
                 !-->        </item>
                 !-->        <item>
                 !-->            <itemName>SINUS</itemName>
                 !-->        </item>
                 !-->        <item>
                 !-->            <itemName>EXTERIOR</itemName>
                 !-->        </item>
                 !-->    </inputElement>
                 !--></inputLine>                 
                 !
                 !.md<2>ROTATION_FUNCTION:    CONSTANT | PARABOLIC | SINUS | EXTERIOR                               $ Rotation function definition
                 !.md<field>ROTATION_FUNCTION
                 !.md<com>Rotation time function definition. It multiplies the AXES_ROTATION vector.
                 !.md<com>Parameters (up to 6 numbers after the type of function)
                 !.md<com>CONSTANT (no parameters, it implies a value of 1.0)
                 !.md<com>PARABOLIC (a t^2 + b t + c)    start_time end_time a b c 
                 !.md<com>SINUS (a cos( w t + phi) + c)  start_time end_time a w phi c
                 !
                 if(exists('CONST')) kfl_fvfua_nsi = 0
                 if(exists('PARAB')) kfl_fvfua_nsi = 1 
                 if(exists('SINUS')) kfl_fvfua_nsi = 2 
                 if(exists('EXTER')) kfl_fvfua_nsi = 3 
                 fvpaa_nsi = param(1:6)

              else if(words(1)=='VELOC') then
                 !--><inputLine>
                 !-->      <inputLineName>Velocity_function</inputLineName>
                 !-->      <inputLineHelp>Velocity time function definition. It multiplies the AXES_VELOCITY vector.</inputLineHelp>
                 !-->      <inputElement>
                 !-->          <inputElementType>combo</inputElementType>
                 !-->          <item>
                 !-->              <itemName>CONSTANT</itemName>
                 !-->          </item>
                 !-->          <item>
                 !-->              <itemName>PARABOLIC</itemName>
                 !-->          </item>
                 !-->          <item>
                 !-->              <itemName>SINUS</itemName>
                 !-->          </item>
                 !-->          <item>
                 !-->              <itemName>EXTERIOR</itemName>
                 !-->          </item>
                 !-->      </inputElement>
                 !-->  </inputLine>           
                 !
                 !.md<2>VELOCITY_FUNCTION:    CONSTANT | PARABOLIC | SINUS | EXTERIOR                               $ Velocity function definition
                 !.md<field>VELOCITY_FUNCTION
                 !.md<com>Velocity time function definition. It multiplies the AXES_VELOCITY vector.
                 !
                 if(exists('CONST')) kfl_fvful_nsi = 0
                 if(exists('PARAB')) kfl_fvful_nsi = 1
                 if(exists('SINUS')) kfl_fvful_nsi = 2
                 if(exists('EXTER')) kfl_fvful_nsi = 3
                 fvpal_nsi = param(1:6)

              else if(words(1)=='CENTE') then
                 !--><inputLine>
                 !-->      <inputLineName>Center_rotation</inputLineName>
                 !-->      <inputLineHelp>Center of rotation definition. Used to compute the centrifugal force.</inputLineHelp>
                 !-->      <inputElement>
                 !-->          <inputElementType>edit2</inputElementType>
                 !-->          <inputElementValueType>REAL</inputElementValueType>                 
                 !-->          <inputLineEditValue>5</inputLineEditValue>
                 !-->      </inputElement>
                 !-->      <inputElement>
                 !-->          <inputElementType>edit2</inputElementType>
                 !-->          <inputElementValueType>REAL</inputElementValueType>                 
                 !-->          <inputLineEditValue>5</inputLineEditValue>
                 !-->      </inputElement>
                 !-->      <inputElement>
                 !-->          <inputElementType>edit2</inputElementType>
                 !-->          <inputElementValueType>REAL</inputElementValueType>                 
                 !-->          <inputLineEditValue>5</inputLineEditValue>
                 !-->      </inputElement>
                 !-->  </inputLine>
                 !
                 !.md<2>CENTER_ROTATION:      real1, real2, real3                                                   $ Center of rotation of coordinates
                 !.md<field>CENTER_ROTATION
                 !.md<com>Center of rotation definition. Used to compute the centrifugal force.
                 !
                 frotc_nsi(1:3) = param(1:3)

              else if((words(1)=='TEMPE'.and.words(2)=='BOUSS').or.&
                   &  (words(1)=='BOUSS'.and.exists('ON   '))) then
                 !--><inputLine>
                 !-->  <inputLineName>Temper_coupling</inputLineName>
                 !-->  <inputLineHelp><![CDATA[TEMPER_COUPLING:
                 !--> Temperature coupling definition. Only Boussinesq coupling is available here.
                 !--> Low Mach system should be chosen as a REGIME (see REGIME option).
                 !--> 
                 !-->     -  real1 = beta: thermal expansion coefficient [1/K].  
                 !-->     -  real2 = Tref: reference temperature [K].  
                 !-->     -  real3 = gb: gravity module [m/s2]. It can be useful for example if g_i /= 0 but |g| = 0.  
                 !--> 
                 !--> The force term added to the ith momentum equations is: f_i = rho * gb * g_i * beta * ( T - Tref ).]]> 
                 !--> </inputLineHelp>
                 !-->  <inputElement>
                 !-->      <inputElementType>combo</inputElementType>
                 !-->      <item>
                 !-->          <itemName>OFF</itemName>
                 !-->      </item>
                 !-->      <item>
                 !-->          <itemName>BOUSSINESQ</itemName>
                 !-->          <itemDependence>0</itemDependence>
                 !-->      </item>
                 !-->  </inputElement>
                 !-->  <inputElement>
                 !-->      <inputElementGroup>0</inputElementGroup>
                 !-->      <inputElementType>edit</inputElementType>
                 !-->          <inputElementValueType>REAL</inputElementValueType>                 
                 !-->      <inputLineEditName>BETA</inputLineEditName>
                 !-->      <inputLineEditValue>5</inputLineEditValue>
                 !-->  </inputElement>
                 !-->  <inputElement>
                 !-->      <inputElementGroup>0</inputElementGroup>
                 !-->      <inputElementType>edit</inputElementType>
                 !-->          <inputElementValueType>REAL</inputElementValueType>                 
                 !-->      <inputLineEditName>TR</inputLineEditName>
                 !-->      <inputLineEditValue>5</inputLineEditValue>
                 !-->  </inputElement>
                 !-->  <inputElement>
                 !-->      <inputElementGroup>0</inputElementGroup>
                 !-->      <inputElementType>edit</inputElementType>
                 !-->          <inputElementValueType>REAL</inputElementValueType>                 
                 !-->      <inputLineEditName>G</inputLineEditName>
                 !-->      <inputLineEditValue>5</inputLineEditValue>
                 !-->  </inputElement>
                 !--></inputLine>
                 !
                 !.md<2>TEMPER_COUPLING:      OFF | BOUSSINESQ [BETA=real1, TR=real2, G=real3]                      $ Coupling with temper module
                 !.md<field>TEMPER_COUPLING
                 !.md<com>Temperature coupling definition. Only Boussinesq coupling is available here.
                 !.md<com>Low Mach system should be chosen as a REGIME (see REGIME option).
                 !.md<com>
                 !.md<com>    -  real1 = beta: thermal expansion coefficient [1/K].  
                 !.md<com>    -  real2 = Tref: reference temperature [K].  
                 !.md<com>    -  real3 = gb: gravity module [m/s2]. It can be useful for example if g_i /= 0 but |g| = 0.  
                 !.md<com>
                 !.md<com>The force term added to the ith momentum equations is: f_i = rho * gb * g_i * beta * ( T - Tref ).
                 !
                 kfl_cotem_nsi = 1
                 boube_nsi     = getrea('BETA ',0.0_rp,'#Beta coefficient')                 
                 boutr_nsi     = getrea('TR   ',0.0_rp,'#Reference temperature')
                 bougr_nsi     = getrea('G    ',0.0_rp,'#Gravity acceleration')
                 gravb_nsi(1)  = getrea('GX   ',0.0_rp,'#x-component of G')
                 gravb_nsi(2)  = getrea('GY   ',0.0_rp,'#y-component of G')
                 gravb_nsi(3)  = getrea('GZ   ',0.0_rp,'#z-component of G')
                 call vecuni(3_ip,gravb_nsi,dummr) 
                 if( dummr < epsilon(1.0_rp) ) then
                    gravb_nsi(1) = gravi_nsi(1)
                    gravb_nsi(2) = gravi_nsi(2)
                    gravb_nsi(3) = gravi_nsi(3)
                 end if
                 if( abs(boube_nsi) < epsilon(1.0_rp) ) boube_nsi= 1.0_rp/boutr_nsi
              else if(words(1)=='TURBU') then

                 if(words(2)=='RANSD'.or.words(2)=='FROMT') then   ! We are using turmu
                    if(exists('TURPR')) kfl_grtur_nsi=1
                 end if

              else if(words(1)=='LEVEL') then
                 !--><inputLine>
                 !-->   <inputLineName>Levels_coupling</inputLineName>
                 !-->   <inputLineHelp>Coupling with Levels module. If this option is ON, then the air properties should also be provided.</inputLineHelp>
                 !--></inputLine>
                 !
                 !.md<2>LEVELS_COUPLING:      ON | OFF                                                              $ Coupling with Levels module
                 !.md<field>LEVELS_COUPLING
                 !.md<com>Coupling with Levels module. If this option is ON, then the air properties should
                 !.md<com>also be provided.
                 !
                 if(words(2)=='ON   ') then
                    kfl_colev_nsi=1
                    if(exists('THICK')) then
                       call runend('NSI_REAPHY: now thickness is read by level and stored in a thick that belongs to master')  
                    end if
                    if(exists('STAGG')) kfl_colev_nsi=3
                    !if(exists('NOHYD')) kfl_nohyd_nsi=1
                    !if(exists('ANALY')) then
                    !   kfl_nohyd_nsi=2
                    !   if( exists('HEIGH') ) &
                    !        heihy_nsi = getrea('HEIGH',0.0_rp,'#Height for hydrostatic pressure')
                    !end if
                 end if

              else if(words(1)=='SURFA') then
                 !--><inputLine>
                 !-->      <inputLineName>SURFACE_TENSION</inputLineName>
                 !-->      <inputLineHelp>If Levels is used, apply surface tension.</inputLineHelp>
                 !-->      <inputElement>
                 !-->          <inputElementType>combo</inputElementType>
                 !-->          <item>
                 !-->              <itemName>On</itemName>
                 !-->              <itemDependence>0</itemDependence>
                 !-->          </item>
                 !-->          <item>
                 !-->              <itemName>Off</itemName>
                 !-->          </item>
                 !-->      </inputElement>
                 !-->      <inputElement>
                 !-->          <inputElementGroup>0</inputElementGroup>
                 !-->          <inputElementType>edit</inputElementType>
                 !-->          <inputElementValueType>REAL</inputElementValueType>                 
                 !-->          <inputLineEditName>Coefficient</inputLineEditName>
                 !-->          <inputLineEditValue>5</inputLineEditValue>
                 !-->      </inputElement>
                 !--></inputLine>
                 !
                 !.md<2>SURFACE_TENSION:      ON, COEFFICIENT=real | OFF                                            $ If Levels is used, apply surface tension
                 !.md<field>SURFACE_TENSION
                 !.md<com>If Levels is used, apply surface tension.
                 !
                 if(words(2)=='ON   ') then
                    kfl_surte_nsi = 1
                    if( exists('COEFF') ) &
                         surte_nsi = getrea('COEFF',0.0_rp,'#Surface tension coeficient')
                 else if (words(2)=='ELSA ') then
                    kfl_surte_nsi = 2
                    if( exists('COEFF') ) then
                         surte_nsi = getrea('COEFF',0.0_rp,'#Surface tension coeficient')
                    endif

                 end if

              else if(words(1)=='FORCE') then  
                 !--><inputLine>
                 !-->     <inputLineName>Force_term</inputLineName>
                 !-->     <inputLineHelp><![CDATA[FORCE_TERM: 
                 !-->        int1 is the material number. char1 is the force model (only WAKE model is available).
                 !-->        The parameters are: Aref,dref,nx,ny,nz]]></inputLineHelp>
                 !-->     <inputElement>
                 !-->         <inputElementType>combo</inputElementType>
                 !-->         <item>
                 !-->             <itemName>MATERIAL</itemName>
                 !-->             <itemDependence>0</itemDependence>
                 !-->         </item>
                 !-->         <item>
                 !-->             <itemName>Off</itemName>
                 !-->         </item>
                 !-->     </inputElement>
                 !-->     <inputElement>
                 !-->         <inputElementType>list</inputElementType>
                 !-->         <inputElementGroup>0</inputElementGroup>
                 !-->         <item>
                 !-->             <itemName>MATERIALS</itemName>
                 !-->             <itemValueType>INT</itemValueType>
                 !-->         </item>
                 !-->         <item>
                 !-->             <itemName>MODEL</itemName>
                 !-->             <itemValueType>CUSTOM:WAKE</itemValueType>                 
                 !-->         </item>
                 !-->         <item>
                 !-->             <itemName>PARAMETERS</itemName>
                 !-->             <itemValueType>REALSEQ</itemValueType>
                 !-->         </item>
                 !-->     </inputElement>
                 !--> </inputLine> 
                 !
                 !.md<2>FORCE_TERM:           OFF | MATERIAL | SPACE=int1 | TIMEF=int2                                                   $ Force term
                 !.md<3>MATERIALS=int3, MODEL=WAKE, PARAMETERS=real1,real2,...                                  $ material force model and parameters 
                 !.md<3>...
                 !.md<2>END_FORCE_TERM
                 !.md<field>FORCE_TERM
                 !.md<com>int1 is the space time function number
                 !.md<com>int2 is the time dependent discrete function number
                 !.md<com>int3 is the material number. WAKE is the force model (only WAKE model is available).
                 !.md<com>The parameters are: Aref,dref,nx,ny,nz
                 !
                 if( words(2) == 'MATER' ) then
                    kfl_force_nsi = 1
                    call nsi_memphy(5_ip)
                    call ecoute('nsi_reaphy')
                    do while( words(1) /= 'ENDFO' )
                       if( words(1) == 'MATER' ) then
                          imate = getint('MATER',1_ip,'#Material force number')
                          if( imate < 1 .or. imate > nmate ) call runend('NSI_REAPHY: WRONG MATERIAL NUMBER')
                          if( exists('WAKE ') ) then
                             lforc_material_nsi(imate) = 1
                          else if (exists('CONST')) then
                             lforc_material_nsi(imate) = 2
                          else
                             call runend('NSI_REAPHY: WRONG MATERIAL FORCE MODEL')
                          end if
                          xforc_material_nsi(1:mforc_material_nsi,imate) = param(4:3+mforc_material_nsi)
                       else if (words(1) == 'TABLE') then  ! loads ct cp table
                          if( imate < 1 .or. imate > nmate ) call runend('NSI_REAPHY: WRONG MATERIAL NUMBER TO LOAD TABLE')
                          icoun = 0          
                          call ecoute('nsi_reaphy')
                          do while( words(1) /= 'ENDTA' )
                             icoun = icoun +1
                             velta_nsi(icoun, imate) = param(1) !velocity
                             thrta_nsi(icoun, imate) = param(2) !thrust coef.
                             powta_nsi(icoun, imate) = param(3) !power  coef.
                             veave_nsi(icoun, imate) = param(4) !simulated  averaged velocity at hub
                             call ecoute('nsi_reaphy')
                          end do
                          ntabl_nsi(imate)= icoun
                          if (icoun > mtabl_nsi) call runend('nsi_reaphy:number of data points in file ctcp distr exceeds maximum' )

                          ! velocity values in table is assumed to be in increasingorder
                          do icoun =1, ntabl_nsi(imate)-1
                             if (velta_nsi(icoun, imate) > velta_nsi(icoun+1, imate)) &
                                  call runend('NSI_REAPHY:WRONG CT TABLE,VELOCITY IS NOT IN INCREA')
                          end do
                       else if (words(1) == 'ROTAT') then  ! loads normal and tangential force distribution for blades
                          if( imate < 1 .or. imate > nmate ) call runend('NSI_REAPHY: WRONG MATERIAL NUMBER TO LOAD TABLE')
                          icoun = 0     
                          call ecoute('nsi_reaphy')
                          do while( words(1) /= 'ENDRO' )
                             if (abs(param(1)-0.0_rp) > 1.0e-7_rp) then
                                icoun = icoun +1
                                radiu_nsi(icoun, imate) = param(1) ! dimensionless radius
                                forcn_nsi(icoun, imate) = param(2) ! normal force
                                forct_nsi(icoun, imate) = param(3) ! tangential force
                             end if
                             call ecoute('nsi_reaphy')
                          end do
                          ntabr_nsi(imate)= icoun  ! when /=0 informs if there exists a force distribution  file
                          if (icoun.gt.mtabl_nsi) call runend('nsi_reaphy:number of data points in file force distr exceeds maximum' )                        
                          ! radial values in table are assumed to be in increasingorder
                          do icoun =1, ntabr_nsi(imate)-1
                             if (radiu_nsi(icoun, imate).gt.radiu_nsi(icoun+1, imate)) &
                                  call runend('NSI_REAPHY:WRONG ROTAT TABLE,RADIUS IS NOT IN INCREA')
                          end do
                       else
                          call runend('NSI_REAPHY: WRONG FORCE FIELD') 
                       end if
                       call ecoute('nsi_reaphy')                       
                    end do  !END FORCE
                 else if( words(2) == 'OFF  ' ) then
                    kfl_force_nsi = 0
                 else if( words(2) == 'SPACE' ) then
                    kfl_force_nsi = -space_time_function_number(getcha('SPACE','NONE ','#Space time function'))
                 else if( words(2) == 'TIMEF' ) then
                    kfl_vegeo_time_nsi =  getint('TIMEF',1_ip,'#Discrete function number')
                 end if
              else if(words(1)=='MASSF') then

                 !--><inputLine>
                 !-->     <inputLineName>MASS_FLOW_CONTROL</inputLineName>
                 !-->     <inputLineHelp><![CDATA[
                 !-->     Mass flow control formula using a PD controller.
                 !-->     The formula modifies the force term based on a target bulk velocity over a set. 
                 !-->     The three input parameters are:
                 !-->     
                 !-->         -   real1 is the target bulk velocity, positive values for inlet 
                 !-->         -   int1 is the set on which we want the target velocity 
                 !-->         -   real2 is the value of the coefficient b of the formula (optimal value is 0.5) 
                 !-->     
                 !-->     The formula is: Fnew=Fold+(b/dt)*(Utarget-2Ucurrent+Uprevious)]]>
                 !-->     </inputLineHelp>
                 !--> </inputLine> 
                 !
                 !.md<2>MASS_FLOW_CONTROL: real1, int1, real2
                 !.md<field>MASS_FLOW_CONTROL
                 !.md<com>real1 is the target bulk velocity, int1 is the set on which we want to impose the target velocity
                 !.md<com>real2 is the value of the coefficient b of the formula (optimal value is 0.5)
                 !
                 kfl_mfrco_nsi = 1
                 mfrub_nsi = param(1)
                 mfrse_nsi = int(param(2),ip)
                 mfccf_nsi = param(3)

              else if (words(1) == 'BOUND') then
                 call runend('NSI_REAPHY: OLD FORMAT TO PROVIDE A SINGLE INLET FUNCTION. DEPRECATED')

              else if (words(1) == 'VALUE') then
                 call runend('NSI_REAPHY: OLD FORMAT TO PROVIDE A SINGLE INLET FUNCTION. DEPRECATED')


              else if (words(1) == 'INLET') then

                 !.md<2>INLET_FUNCTIONS
                 !.md<3>SIZES, TABLES=N
                 !.md<4>TABLE=1, NODES=m_1, STEPS=n_1, TIME_STEP = dt_1
                 !.md<4>TABLE=2, NODES=m_2, STEPS=n_2, TIME_STEP = dt_2
                 !.md<4>...
                 !.md<4>TABLE=N, NODES=m_N, STEPS=n_N, TIME_STEP = dt_N
                 !.md<3>END_SIZES
                 !.md<3>
                 !.md<3>TABLE:1
                 !.md<4>BOUNDARY_NODES
                 !.md<5>1    b1_1
                 !.md<5>2    b2_1
                 !.md<5>...
                 !.md<5>m_1  bm1_1
                 !.md<4>END_BOUNDARY_NODES
                 !.md<4>VALUES_ON_BOUNDARY
                 !.md<5>1    1    vx1_t1     vy1_t1     vz1_t1     $$ First time step
                 !.md<5>...
                 !.md<5>m_1  1    vxm1_t1    vym1_t1    vzm1_t1    $$ First time step
                 !.md<5>1    2    vx1_t2     vy1_t2     vz1_t2     $$ Second time step
                 !.md<5>...
                 !.md<5>m_1  2    vxm1_t2    vym1_t2    vzm1_t2    $$ Second time step
                 !.md<5>...
                 !.md<5>m_1  n_1  vxm1_tn_1  vym1_tn_1  vzm1_tn_1  $$ Last time step
                 !.md<4>END_VALUES_ON_BOUNDARY
                 !.md<3>END_TABLE
                 !.md<3>
                 !.md<3>...
                 !.md<3>
                 !.md<3>TABLE:N
                 !.md<4>BOUNDARY_NODES
                 !.md<5>1    b1_N
                 !.md<5>2    b2_N
                 !.md<5>...
                 !.md<5>m_N  bm_N_N
                 !.md<4>END_BOUNDARY_NODES
                 !.md<4>VALUES_ON_BOUNDARY
                 !.md<5>1    1    vx1_t1     vy1_t1     vz1_t1
                 !.md<5>...
                 !.md<5>m_N  1    vxm_N_t1   vym_N_t1   vzm_N_t1
                 !.md<5>1    2    vx1_t2     vy1_t2     vz1_t2
                 !.md<5>...
                 !.md<5>m_N  2    vxm_N_t2   vym_N_t2   vzm_N_t2
                 !.md<5>...
                 !.md<5>m_N  n_N  vxm_N_tn_N vym_N_tn_N vzm_N_tn_N
                 !.md<4>END_VALUES_ON_BOUNDARY
                 !.md<3>END_TABLE
                 !.md<3>
                 !.md<2>END_INLET_FUNCTIONS
                 !.md<field>INLET_FUNCTIONS
                 !.md<com>Impose temporal and spatial varying velocities for several possible inlets from tabulated
                 !.md<com>fields.
                 !.md<com>   - SIZES: provide the number of tables (inlets with varying velocities) and for each table
                 !.md<com>            the number of nodes at the inlet, number of time steps and time step size.
                 !.md<com>   - TABLE: for each table list the nodes of the inlet (BOUNDARY_NODES) and the velocity
                 !.md<com>            components for each node and time step (VALUES_ON_BOUNDARY).

                 call ecoute('nsi_reaphy')
                 do while (words(1) /= 'ENDIN')

                    ! Read table sizes
                    if (words(1) == 'SIZES') then
                       if (words(2) == 'TABLE') then
                          !ntabf_nsi = int(param(2),ip)
                          ntabf_nsi = param(2)
                          call nsi_memphy(8_ip)
                          kfl_bnods_nsi = 1_ip
                       else
                          call runend('NSI_REAPHY: NUMBER OF TABLES NOT PROVIDED IN INLET FUNCTIONS')
                       end if

                       call ecoute('nsi_reaphy')
                       do while (words(1) /= 'ENDSI')
                          if (words(1) == 'TABLE') then
                             ii = int(param(1),ip)
                             if (ii > ntabf_nsi) call runend('NSI_REAPHY: NUMBER OF TABLE EXCEEDS THE MAXIMUM IN INLET FUNCTIONS')
                             if (words(2) == 'NODES') then
                                nbnod_pos_nsi(ii+1) = int(param(2),ip)
                             else
                                call runend('NSI_REAPHY: NUMBER OF NODES NOT PROVIDED IN INLET FUNCTIONS')
                             end if

                             if (words(3) == 'STEPS') then
                                nbtim_pos_nsi(ii+1) = int(param(3),ip)
                              else
                                call runend('NSI_REAPHY: NUMBER OF STEPS NOT PROVIDED IN INLET FUNCTIONS')
                             end if

                             if (words(4) == 'TIMES') then
                                nbtdt_nsi(ii) = param(4)
                             else
                                call runend('NSI_REAPHY: NUMBER OF STEPS NOT PROVIDED IN INLET FUNCTIONS')
                             end if

                          end if
                          call ecoute('nsi_reaphy')
                       end do

                       do ii = 2, ntabf_nsi+1_ip
                          nbtim_nod_pos_nsi(ii) = nbtim_nod_pos_nsi(ii-1_ip) + &
                             nbnod_pos_nsi(ii) * nbtim_pos_nsi(ii)
                          nbnod_pos_nsi(ii)     = nbnod_pos_nsi(ii-1_ip) + nbnod_pos_nsi(ii)
                          nbtim_pos_nsi(ii)     = nbtim_pos_nsi(ii-1_ip) + nbtim_pos_nsi(ii)
                       end do

                       nbnod_nsi = nbnod_pos_nsi(ntabf_nsi+1_ip)
                       nbtim_nod_nsi = nbtim_nod_pos_nsi(ntabf_nsi+1_ip) 

                       ! Allocate memory
                       call nsi_memphy(6_ip)
                       call nsi_memphy(7_ip)


                    ! Fill tables
                    else if (words(1) == 'TABLE') then
                       if (kfl_bnods_nsi /= 1_ip) call runend('NSI_REAPHY: SIZES FOR TABLES HAVE TO BE GIVEN FIRST IN INLET FUNCTIONS')
                       ii = param(1)  ! Table number
                       if (nbnod_pos_nsi(ii+1) == 0_ip) then
                          call runend('NSI_REAPHY: NO SIZES PROVIDED FOR TABLE IN INLET FUNCTIONS')
                       else if (ii > ntabf_nsi) then  
                          call runend('NSI_REAPHY: NUMBER OF TABLE EXCEEDS THE MAXIMUM IN INLET FUNCTIONS')
                       else
                          call ecoute('nsi_reaphy')
                          do while (words(1) /= 'ENDTA')
                             if( words(1) == 'BOUND' ) then
                                icoun = nbnod_pos_nsi(ii)
                                call ecoute('nsi_reaphy')
                                do while(words(1)/='ENDBO')
                                   icoun = icoun + 1_ip
                                   bntab_nsi(icoun,1_ip) = param(2_ip)
                                   bntab_nsi(icoun,2_ip) = param(2_ip)
                                   bntab_nsi(icoun,3_ip) = param(2_ip)
                                   call ecoute('nsi_reaphy')
                                end do
                                if (icoun /= nbnod_pos_nsi(ii+1)) call runend('NSI_REAPHY: MISMATCH BETWEEN THE NUMBER OF NODES IN FILE AND THE SPECIFICED ONE IN INLET FUNCTIONS')

                             else if( words(1) == 'VALUE' ) then
                                icoun = nbtim_nod_pos_nsi(ii)
                                call ecoute('nsi_reaphy')
                                do while(words(1)/='ENDVA')
                                   icoun = icoun + 1_ip
                                   bnval_nsi(icoun,1_ip) = param(3_ip)
                                   bnval_nsi(icoun,2_ip) = param(4_ip)
                                   bnval_nsi(icoun,3_ip) = param(5_ip)
                                   call ecoute('nsi_reaphy')
                                end do
                                if (icoun /= nbtim_nod_pos_nsi(ii+1)) call runend('NSI_REAPHY: MISMATCH BETWEEN THE NUMBER OF STEPS*NODES IN FILE AND THE SPECIFICED ONES IN INLET FUNCTIONS')
                             end if
                             call ecoute('nsi_reaphy')
                          end do
                       end if
                    end if

                    call ecoute('nsi_reaphy')

                 end do

              end if

              call ecoute('nsi_reaphy')
              !--> </subGroup>
              !
              !.md<1>END_PROBLEM_DEFINITION
              !
           end do

        else if(words(1)=='PROPE') then
           !--><subGroup>
           !-->     <subGroupName>PROPERTIES</subGroupName>
           !-->     <subGroupHelp>PROPERTIES:
           !-->           Properties. Here, some properties could be required if Nastin is coupled with another
           !-->           module. The fluid properties (density, viscosity) should be defined in Kermod.</subGroupHelp>
           !
           !.md<1>PROPERTIES:
           !.md<field>PROPERTIES
           !.md<com>Properties. Here, some properties could be required if Nastin is coupled with another
           !.md<com>module. The fluid properties (density, viscosity) should be defined in Kermod.
           !
           call ecoute('nsi_reaphy')
           imate=1
           do while(words(1)/='ENDPR')
              if( words(1) == 'MATER' ) then
                 imate = getint('MATER',1_ip,'#Current material')
                 if( imate > nmate ) then
                    print*,'nmate=',nmate
                    call runend('NSI_REAPHY: WRONG MATERIAL')
                 end if
              else if(words(1)=='GASCO') then               ! Gas constant (R)
                 gasco = getrea('GASCO',287.0_rp,'#Gas constant')
              else if(words(1)=='SPECI') then               ! Specific heat (Cp)
                 sphea_nsi = getrea('SPECI',1006.1_rp,'#Gas constantSpecific heat')

              else if(words(1)=='THERM') then               ! Thermodynamic pressure
                 !--><inputLine>
                 !-->      <inputLineName>THERMODYNAMIC_PRESSURE</inputLineName>
                 !-->      <inputLineHelp>Thermodynamic pressure "real" for Low Mach regime. - Closed system: choose MASS_CONSERVATION. Thermodynamics
                 !--->     pressure is computed from mass conservation equation. "real" is therefore the initial pressure. - Open system without.</inputLineHelp>
                 !-->      <inputElement>
                 !-->          <inputElementType>edit2</inputElementType>
                 !-->          <inputElementValueType>REAL</inputElementValueType>                 
                 !-->          <inputLineEditValue>5</inputLineEditValue>
                 !-->      </inputElement>
                 !-->      <inputElement>
                 !-->          <inputElementType>combo</inputElementType>
                 !-->          <item>
                 !-->              <itemName>MASS_CONSERVATION</itemName>
                 !-->          </item>
                 !-->          <item>
                 !-->              <itemName>TIME_DEPENDENT</itemName>
                 !-->          </item>
                 !-->      </inputElement>
                 !-->  </inputLine>
                 !
                 !.md<2>THERMODYNAMIC_PRESSURE: real, MASS_CONSERVATION | TIME_DEPENDENT
                 !.md<field>THERMODYNAMIC_PRESSURE
                 !.md<com>Thermodynamic pressure "real" for Low Mach regime.
                 !.md<com>    - Closed system: choose MASS_CONSERVATION. Thermodynamics pressure is computed from mass conservation equation.
                 !.md<com>      "real" is therefore the initial pressure.
                 !.md<com>    - Open system without.  
                 !
                 prthe_nsi=param(1) !getrea('THERM',101325.0_rp,'#Thermodynamic pressure')
                 if(exists('MASSC')) then
                    kfl_prthe_nsi=1 ! Computed from mass conservation, closed system
                    tmass_nsi = getrea('MASSC',0.0_rp,'#initial density')
                 else if(exists('TIMED')) then
                    kfl_prthe_nsi=2 ! Computed from evolution equation, closed system
                    tmass_nsi = getrea('TIMED',0.0_rp,'#initial density')
                 else
                    kfl_prthe_nsi=0
                 end if

              else if(words(1)=='INITI') then
                 !
                 ! Initial Fields
                 !
                 if(words(2)=='CODES') then ! Initial Fields given by codes
                    call runend('NSI_REAPHY: "CODES" OPTION IN INITIAL_CONDITION IS NOT AVAILABLE')
                 else ! Initial Fields given on nodes
                    call ecoute('nsi_reaphy')
                    do while(words(1)/='ENDIN')
                       if (words(1) == 'VEFOR') then
                          nfiel_nsi(2) = -getint('FIELD',1_ip,'#Field Number for forward velocities')
                       else if (words(1) == 'PRFOR') then
                          nfiel_nsi(1) = -getint('FIELD',1_ip,'#Field Number for forward pressure')
                       end if
                       call ecoute('nsi_reaphy')
                    end do
                 end if


              end if
              call ecoute('nsi_reaphy')
           end do
        end if
        !--></subGroup>
        !
        !.md<1>END_PROPERTIES
        !
     end do

     ! load anisotropic porosity on if it exists in the property in kernel
     if (anipo_ker % kfl_exist ==1 )  kfl_anipo_nsi = 1_ip
     !--></group>
     !
     !.md<0>END_PHYSICAL_PROBLEM
     !.md</code>
  end if
end subroutine nsi_reaphy
