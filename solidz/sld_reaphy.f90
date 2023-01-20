!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup SolidzInput
!> @ingroup    Solidz
!> @{
!> @file    sld_reaphy.f90
!> @author  Mariano Vazquez
!> @date    July, 2014
!>          - Subroutine written
!> @brief   Read physical data Solidz
!> @details Read physical data of the problem and assign defaults
!> @}
!-----------------------------------------------------------------------

subroutine sld_reaphy()
  !
  !.md<module>solidz
  !.md<input>case.sld.dat
  !.md<pos>0
  !.md<sec>
  !
  use def_master,                only: INOTSLAVE ,                  &
       intost, retost
  use def_master,                only: zeror
  use def_domain,                only: ndime
  use mod_ecoute,                only: ecoute
  use def_kintyp,                only: ip,rp
  use def_inpout,                only: nnpar, words,                &
       param, exists
  use def_inpout,                only: getint, getrea, getcha
  use mod_maths,                 only : maths_normalize_vector
  use mod_sld_stress_model_152,  only : sm152_precalculus
  use mod_sld_stress_model_154,  only : sm154_precalculus
  use mod_sld_interface_element, only : ELINT_LAW_TURON ,           &
       ELINT_LAW_TURON_2018
  use mod_sld_interface_element, only : ELINT_TANGENT_EXPLICIT,     &
       ELINT_TANGENT_ANALYTICAL,   &
       ELINT_TANGENT_NUMERICAL
  use mod_sld_interface_element, only : ELINT_TANGENT_SECANT
  use def_solidz,                only : kfl_timei_sld
  use def_solidz,                only : nmate_sld, ncoef_sld
  use def_solidz,                only : parco_sld, parch_sld,      &
       parcf_sld
  use def_solidz,                only : lawst_sld,                  &
       lawta_sld, lawco_sld,       &
       lawmo_sld, lawho_sld
  use def_solidz,                only : lawch_sld
  use def_solidz,                only : lawpl_sld
  use def_solidz,                only : densi_sld, rbmas_sld
  use def_solidz,                only : gravi_sld, grnor_sld, thick_sld
  use def_solidz,                only : modor_sld, modfi_sld
  use def_solidz,                only : kfl_rigid_sld
  use def_solidz,                only : kfl_restr_sld, kfl_indis_sld
  use def_solidz,                only : kfl_plane_sld, kfl_sdvar_sld, kfl_damag_sld
  use def_solidz,                only : kfl_cohes_sld
  use def_solidz,                only : kfl_cshel_sld
  use def_solidz,                only : SLD_CSHEL_EAS, SLD_CSHEL_ANS_SHEAR, SLD_CSHEL_ANS_TRAPEZOIDAL
  use def_solidz,                only : kfl_tange_sld, kfl_fiber_sld, kfl_prdef_sld
  use def_solidz,                only : kfl_vofor_sld, kfl_moduf_sld, kfl_isoch_sld
  use def_solidz,                only : kfl_strai_sld
  use def_solidz,                only : thiso_sld
  use def_solidz,                only : kfl_dampi_sld
  use def_solidz,                only : kfl_plast_sld
  use def_solidz,                only : kfl_donna_sld
  use def_solidz,                only : dampi_sld
  use def_solidz,                only : kfl_csysm_sld, csysm_sld, oripa_sld, kfl_rmate_sld
  use def_solidz,                only : SLD_CSYS_CARTESIAN, SLD_CSYS_CYLINDRICAL
  use def_solidz,                only : SLD_TANGENT_ANALYTICAL,    &
       SLD_TANGENT_CHECK,         &
       SLD_TANGENT_NUMERICAL,     &
       SLD_SECANT
  use def_solidz,                only : SLD_DYNAMIC_PROBLEM,       &
       SLD_STATIC_PROBLEM
  use def_solidz,                only : SLD_INFINITESIMAL,         &
       SLD_GREEN
  use mod_sld_stress_model_152,  only : sm152_precalculus
  use mod_sld_interface_element, only : ELINT_LAW_TURON ,          &
       ELINT_LAW_TURON_2018,      &
       ELINT_LAW_CONTACT
  use mod_sld_interface_element, only : ELINT_TANGENT_EXPLICIT,    &
       ELINT_TANGENT_ANALYTICAL,  &
       ELINT_TANGENT_NUMERICAL,   &
       ELINT_TANGENT_SECANT
  use mod_sld_fe2,               only : MAT_MICRO_NO_COUPLING,MAT_MICRO_ONE_WAY,MAT_MICRO_FULL
  use mod_sld_rbo,               only : sld_rbo_inipro
  use mod_messages,              only : messages_live,livinf
  use mod_sld_cardiac_cycle,     only : sld_cardiac_cycle_read_data
  use mod_sld_cardiac_cycle,     only : sld_cardiac_cycle_initialise_variables
  use mod_eccoupling,            only : eccou_show_warning
  use mod_biofibers,             only : biofib_display_warning
  use mod_sld_atm,               only : kfl_therm_sld
#ifndef PROPER_ELEM_PRIVATE_OFF
  use mod_ker_sysnet, only: sysnet_read_data
#endif

  implicit none

  external :: runend
  external :: sld_memphy
  external :: sm151_precalculus
  external :: sm153_precalculus

  integer(ip)     :: imate,nmatm,nposi
  real(rp)        :: K0,E,poiss
  character(300)  :: messa
  character(5)    :: wcsys, tcsys

  if( INOTSLAVE ) then

     kfl_timei_sld = SLD_DYNAMIC_PROBLEM
     kfl_rigid_sld = 0_ip                     ! Rigid body (0 = 0FF)
     kfl_tange_sld = SLD_TANGENT_ANALYTICAL   ! Analytical calculation of the tangent moduli (only for Implicit)
     kfl_strai_sld = SLD_GREEN                ! Large strains using Green-Lagrange
     kfl_sdvar_sld = 0_ip                     ! State dependent variables flag activation (0 = OFF; 1 = ON)
     kfl_damag_sld = 0_ip                     ! Damage model activation flag (0 = OFF; 1 = ON)
     kfl_cohes_sld = 0_ip                     ! Interface cohesive elements activation flag (0 = OFF; 1 = ON)
     kfl_fiber_sld = 0_ip                     ! Fiber field (0 if none)
     kfl_restr_sld = 0_ip                     ! Initial stress fields (ndime size)
     kfl_indis_sld = 0_ip                     ! Initial displacements
     kfl_prdef_sld = 0_ip                     ! Fibers field nodal- or element- wise
     kfl_vofor_sld = 0_ip                     ! External volume force
     kfl_moduf_sld = 0_ip                     ! No modulating fields
     kfl_plane_sld = 0_ip                     ! Default: plane strain assumption for 2D problems
     kfl_isoch_sld = 0_ip                     ! Isochrones
     kfl_cshel_sld = 0_ip                     ! Parameters continuum shell element
     kfl_rmate_sld = 0_ip                     ! Material axes rotation
     kfl_donna_sld = 0_ip                     ! Donnan Osmosis for porous hyperelastic swelling
     kfl_therm_sld = .false.
     thick_sld     = 1.0_rp                   ! Out-of-plane thickness (2-d problems)
     thiso_sld     = -10000000.0_rp
     grnor_sld     = 0.0_rp                   ! Gravity norm
     gravi_sld(:)  = 0.0_rp                   ! Gravity vector
     oripa_sld(:)  = 0_ip                     ! Orientation parameters
     csysm_sld(:)  = 0.0_rp                   ! Coordinate system point data
     rbmas_sld     = 1.0_rp                   ! Rigid body mass

     call sld_cardiac_cycle_initialise_variables()
     !
     ! Reach the section
     !
     call ecoute('sld_reaphy')
     do while(words(1)/='PHYSI')
        call ecoute('sld_reaphy')
     end do
     !
     ! Begin to read data
     !
     !
     !.md# Physical Properties Definition
     !.md<code>
     !.md<0>PHYSICAL_PROBLEM
     !
     messa = '        READING PHYSICS...'
     call livinf(0_ip,messa,1_ip)

     do while(words(1)/='ENDPH')
        call ecoute('sld_reaphy')

        if (      words(1) == 'PROBL' ) then
           !
           !.md<1>PROBLEM_DEFINITION
           !.md<com>
           !.md<com><strong>Problem definition</strong>
           !
           call sld_memphy(0_ip)
           call ecoute('sld_reaphy')
           do while(words(1)/='ENDPR')
              if(      words(1) == 'TEMPO' ) then
                 !
                 !.md<2>TEMPORAL_DERIVATIVES: DYNAMIC | STATIC                     $ Existence of temporal derivatives
                 !.md<field>TEMPORAL_DERIVATIVES
                 !.md<com>Set this parameter equal to DYNAMIC for solving transient problems or STATIC
                 !.md<com>for solving equilibrium problems. By default, Alya solves DYNAMIC (transient) problems.
                 !
                 if( words(2)=='DYNAM' ) then
                    kfl_timei_sld = SLD_DYNAMIC_PROBLEM
                    messa = '        SOLIDZ SOLVING DYNAMIC PROBLEM.'
                    call livinf(0_ip,messa,1_ip)
                 else
                    kfl_timei_sld = SLD_STATIC_PROBLEM
                    messa = '        SOLIDZ SOLVING STATIC PROBLEM.'
                    call livinf(0_ip,messa,1_ip)
                 end if

              else if( words(1) == 'RIGID' ) then
                 !
                 !.md<2>RIGID_BODY:           ON | OFF                             $ Rigid body
                 !.md<field>RIGID_BODY
                 !.md<com>Set this parameter equal to ON to consider your mesh as a rigid body.
                 !.md<com>By default this parameter is set to OFF.
                 !
                 if( words(2)=='ON   ' ) then
                    kfl_rigid_sld = 1_ip
                    call sld_rbo_inipro(0_ip)
                    messa = '        SOLIDZ SOLVING A RIGID BODY PROBLEM.'
                    call livinf(0_ip,messa,1_ip)
                 else if ( words(2)=='OFF  ' ) then
                    kfl_rigid_sld = 0_ip
                 end if

              else if( words(1) == 'GRAVI' ) then
                 !
                 !.md<2>GRAVITY:              NORM=real, GX=real, GY=real, GZ=real $ Gravity acceleration
                 !.md<field>GRAVITY
                 !.md<com>This option is used to specify gravity loads. Set NORM=<tt>real</tt> to define the vector's magnitude.
                 !.md<com>Set GX=<tt>real</tt>, GY=<tt>real</tt> and GZ=<tt>real</tt> to unity vectors which define
                 !.md<com>the orientation and the sense of the vector.
                 !
                 grnor_sld    = getrea('NORM ',0.0_rp,'#Gravity norm')
                 gravi_sld(1) = getrea('GX   ',0.0_rp,'#x-component of g')
                 gravi_sld(2) = getrea('GY   ',0.0_rp,'#y-component of g')
                 if( ndime == 3_ip ) then
                    gravi_sld(3) = getrea('GZ   ',0.0_rp,'#z-component of g')
                 end if
                 call maths_normalize_vector(ndime,gravi_sld)

              else if( words(1) == 'PLANE' ) then
                 !
                 !.md<2>PLANE:                STRESS | STRAIN                      $ Plane stress/strain assumptions (2-d problems)
                 !.md<field>PLANE
                 !.md<com>Plane stress/strain assumptions for 2-d problems. By default Alya assumes plane strain conditions.
                 !
                 if( words(2)=='STRES' )  then
                    kfl_plane_sld = 1_ip
                 else if( words(2)=='STRAI') then
                    kfl_plane_sld = 0_ip
                 end if

              else if( words(1) == 'THICK' ) then
                 !
                 !.md<2>THICKNESS_OUT_OF_PLANE= real                               $ Out-of-plane thickness (2-d problems)
                 !.md<field>THICKNESS_OUT_OF_PLANE
                 !.md<com>Out-of-plane thickness for 2-d problems. By default Alya assumes 1.0.
                 !
                 thick_sld = getrea('THICK',1.0_rp,'#Out-of-plane thickness')

              else if( words(1) == 'NLGEO' ) then
                 !
                 !.md<2>NLGEOM: ON | OFF                                           $ Nonlinear deflection effects
                 !.md<field>Activates large-deflection effects in a static or transient problems.
                 !.md<com>
                 !
                 if(      words(2) == 'ON   ' ) then
                    kfl_strai_sld = SLD_GREEN
                 else if( words(2) == 'OFF  ' ) then
                    kfl_strai_sld = SLD_INFINITESIMAL
                 end if
                 
              else if( words(1) == 'THERM' ) then
                 !
                 !.md<2>THERMAL_ANALYSIS: ON | OFF                                 $ Thermal analysis
                 !.md<field>Activates thermal effects for thermomechanical coupling.
                 !.md<com>
                 !              
                 if(      words(2) == 'ON   ' ) then
                    kfl_therm_sld = .true.
                 else if( words(2) == 'OFF  ' ) then
                    kfl_therm_sld = .false.
                 end if

              end if

              call ecoute('sld_reaphy')
           end do
           !
           !.md<1>END_PROBLEM_DEFINITION
           !
        else if ( words(1) == 'TANGE' ) then
           !
           ! ADO[1]> TANGENT_STIFFNESS_CALCULATION: CHECK | ANALYTICAL | NUMERICAL | SECANT
           ! ADO[d]> TANGENT:
           ! ADO[d]> This option is used to choose the calculation of the tangent stiffness matrix.
           ! ADO[d]> Set to ANALYTICAL to use the analytical calculation of the stiffness matrix.
           ! ADO[d]> Set to NUMERICAL to calculatet it numerically.
           ! ADO[d]> Set to SECANT to use the secant method.
           ! ADO[d]> This option can be used only for Implicit analysis. By default, the tangent
           ! ADO[d]> stiffness matrix is calculated analytically.
           if (exists('CHECK')) then
              kfl_tange_sld = SLD_TANGENT_CHECK
           else if (exists('ANALY')) then
              kfl_tange_sld = SLD_TANGENT_ANALYTICAL
           else if (exists('NUMER')) then
              kfl_tange_sld = SLD_TANGENT_NUMERICAL
           else if (exists('SECAN')) then
              kfl_tange_sld = SLD_SECANT
           else
              ! <GGU> Old way. To be deleted in the future
              kfl_tange_sld = getint('TANGE',1_ip,'Use numerical tangent')
           end if

        else if ( words(1) == 'PROPE' ) then
           !
           !.md<1>PROPERTIES
           !.md<com>
           !.md<com><strong>Properties</strong>
           !
           call sld_memphy(1_ip)
           call ecoute('sld_reaphy')

           imate = 1_ip
           messa = '        READING MATERIALS...'
           call livinf(0_ip,messa,1_ip)
           do while(words(1)/='ENDPR')

              if(       words(1) == 'MATER' ) then
                 !
                 !.md<2>MATERIAL= int                                  $ Material code
                 !.md<field>MATERIAL
                 !.md<com>Set this parameter equal to the material code number. By default, Alya assumes one material.
                 !
                 imate = getint('MATER=',1_ip,'#Current material')
                 messa = &
                      '        MATERIAL: '//trim(intost(imate))
                 call livinf(0_ip,messa,1_ip)

                 if( imate > nmate_sld ) then
                    call runend('SLD_REAPHY: THIS MODULE HAS MORE MATERIALS THAN THOSE DEFINED IN DOM.DAT')
                 end if

              else if ( words(1) == 'DENSI' ) then
                 !
                 !.md<2>DENSITY=  real                                 $ Density
                 !.md<field>DENSITY
                 !.md<com>Set this parameter equal to the density of the material. This option is required
                 !.md<com>for dynamic (transient) problems.
                 !
                 densi_sld(1:ncoef_sld,imate)=param(1:ncoef_sld)
                 messa = &
                      '           DENSITY= '//trim(retost(densi_sld(1,imate)))
                 call livinf(0_ip,messa,1_ip)

              else if ( words(1) == 'MASS ' ) then
                 !
                 !.md<2>MASS=     real                                 $ Mass
                 !.md<field>MASS
                 !.md<com>Set this parameter equal to the mass of the rigid body. Rigid body option
                 !.md<com>have to be activated.
                 !
                 if ( kfl_rigid_sld == 1_ip ) then
                    rbmas_sld = getrea('MASS ',1.0_rp,'Mass of the rigid body')
                    messa = '           MASS= '//trim(retost(rbmas_sld))
                    call livinf(0_ip,messa,1_ip)

                 end if

              else if ( words(1) == 'CONST' ) then
                 !
                 !.md<2>CONSTITUTIVE_MODEL:                            $ Material law
                 !.md<field>CONSTITUTIVE_MODEL
                 !.md<com>This option is to define the material model. The available material models
                 !.md<com>are the following:
                 !.md<com>

                 ! ADOC[d]> <li> NEOHOOK BELY (Neohookean hyperelastic): YOUNG, lambda, mu </li>
                 ! ADOC[d]> Corresponds to file 101.
                 !
                 ! ADOC[d]> <li> NEOHOOK POROU [DONNA] (Neohookean porous [with donnan osmotic swelling]): mu, ns0 [R, T, 
                 ! ADOC[d]> c_F0, c_b, t_ramp] </li>
                 ! ADOC[d]> Corresponds to file 400.
                 !
                 ! ADOC[d]> <li> NEOHOOK ANSYS (Neohookean hyperelastic): YOUNG, mu, K </li>
                 ! ADOC[d]> Corresponds to file 102.
                 !
                 ! ADOC[d]> <li> MOONEY 1/2/3 (Mooney-Rivlin): there exists 3 formulations </li>
                 ! ADOC[d]> MR1: YOUNG, C_1, C_2, D_1
                 ! ADOC[d]> MR2: YOUNG, C_1, C_2, D_1
                 ! ADOC[d]> MR3: YOUNG, C_1, C_, *, *, D_1
                 ! ADOC[d]> where C1,C2 and D1 are input parameters. If C2=0 we recover Neo-Hookean law.
                 ! ADOC[d]> This material model is useful for modeling certain types of isotropic biological tissues.
                 ! ADOC[d]> For example, the isotropic collagen matrix can be described quite well with this material model.
                 ! ADOC[d]> Corresponds to file 103.
                 !
                 ! ADOC[d]> <li> SPRIN (Constant model): lambda0 </li>
                 ! ADOC[d]> Corresponds to file 105.
                 !
                 ! ADOC[d]> <li> LINVI: (Yin & Lin with Peteron's EC coupling integrated) </li>
                 ! ADOC[d]> Corresponds to file 133.
                 !
                 ! ADOC[d]> <li> HOLZA (Transversaly isotropic hyperelastic): E_eq,k,a,b,af,bf,as,bs,afs,bfs,scalf.
                 ! ADOC[d]> E_eq is an equivalent rigidity used to calculate the time step. k is the compresibility factor.
                 ! ADOC[d]> scalf scales the stress tensor (defualt to 1) </li>
                 ! ADOC[d]> Corresponds to file 134
                 !
                 ! ADOC[d]> <li> GUCCI (Guccione): E, Kct, c1, c2, c3, c4 </li>
                 ! ADOC[d]> Corresponds to file 136.
                 !
                 !
                 ! ADOC[d]> <li> OLIVE (Isotropic Damage model of Oliver et al. (1990)): Exx, nuxx, Gxx, xT, xC, gF, idSur, idLaw. 
                 ! ADOC[d]> Corresponds to file 153.
                 !
                 ! ADOC[d]> <li> CSHEL (Continuum shell element):</li>
                 ! ADOC[d]> Constitutive models
                 !
                 ! ADOC[d]> <ul>
                 ! ADOC[d]> <li> ISOLI. Isotropic linear elastic sm100 for Continuum Shell. </li>
                 ! ADOC[d]> <li> ORTHO. Orthotropic linear elastic for Continuum Shell. </li>
                 ! ADOC[d]> <li> MAIMI. Maimi 3D damage model for Continuum Shell. </li>
                 ! ADOC[d]> </ul>
                 ! ADOC[d]> Corresponds to file 200.
                 !
                 ! ADOC[d]> <li> MICRO (microscopic fem model using micropp library):. </li>
                 ! ADOC[d]> <ul>
                 ! ADOC[d]> <li> INFIN. Uses Infinitesimal (Small) strain theory assumption. Strains and rotations are both small. 
                 ! ADOC[d]> </ul>
                 ! ADOC[d]> Corresponds to file MAT_MICRO.
                 ! ADOC[d]> </ul>
                 !
                 messa = '           MATERIAL MODEL: '
                 if(       words(2) == 'NEOHO' ) then
                    lawco_sld(imate)=1_ip
                    lawst_sld(imate)=101_ip   ! default, belytschko
                    parco_sld(1:ncoef_sld,imate)=param(2:ncoef_sld+1)

                    messa = adjustr(trim(messa))// ' NEOHOKEAN '

                    if(words(3)=='BELYT') then
                       messa = adjustr(trim(messa))// ' BELYTSCHKO '
                       lawst_sld(imate)=101_ip
                       parco_sld(1:ncoef_sld,imate)=param(3:ncoef_sld+2)
                    else if(words(3)=='ANSYS') then
                       messa = adjustr(trim(messa))// ' ANSYS '
                       lawst_sld(imate)=102_ip
                       parco_sld(1:ncoef_sld,imate)=param(3:ncoef_sld+2)
                    else if(words(3)=='CHAPE') then
                       messa = adjustr(trim(messa))// ' CHAPE '
                       lawst_sld(imate)=121_ip
                       parco_sld(1:ncoef_sld,imate)=param(3:ncoef_sld+2)
                    else if(words(3)=='POROU') then
                    !
                    !.md<com>    - NEOHOOKEAN POROUS [DONNAN FIBER]: Neohookean porous 
                    !.md<com>[with donnan osmotic swelling and/or fibers]. To use fibers 
                    !.md<com>without donnan, placeholders for DONNA and its parameters are required.
                    !.md<com>Corresponds to file 400.
                    !.md<com>The inputs are the following:
                    !.md<com><tt>G, ns_0</tt> 
                    !.md<com><tt>[[R, T, c_F0, c_b, t_ramp], [k_1, k_2]]</tt>
                    !.md<com>
                    !
                       messa = adjustr(trim(messa))// ' POROUS '
                       lawst_sld(imate)=400_ip
                       parco_sld(1:ncoef_sld,imate)=param(3:ncoef_sld+2)
                       if(words(4) == 'DONNA') then
                          messa = adjustr(trim(messa))// ' WITH SWELLING '
                          kfl_donna_sld = 1_ip
                          parco_sld(1:ncoef_sld,imate)=param(4:ncoef_sld+3)
                       end if
                       if(words(5) == 'FIBER') then
                          messa = adjustr(trim(messa))// ' WITH FIBER FIELDS '
                          parco_sld(1:ncoef_sld,imate)=param(5:ncoef_sld+4)
                       end if
                    end if

                 else if ( words(2) == 'ISOLI' .or. words(2) == 'ISOTR' ) then
                    !
                    !.md<com>    - ISOLIN: Isotropic linear elastic model.
                    !.md<com>Corresponds to file 100.
                    !.md<com>The inputs are the following:
                    !.md<com><tt>E, nu</tt>
                    !.md<com>
                    !
                    messa = adjustr(trim(messa))// ' ISOLINEAR '
                    lawco_sld(imate)=1_ip
                    lawst_sld(imate)=100_ip   ! isotropic linear elastic
                    parco_sld(1:ncoef_sld,imate)=param(2:ncoef_sld+1)
                    ! Deactivate thermal effects
                    if( .not. kfl_therm_sld ) parco_sld(3,imate) = 0.0_rp
                    
                    
                 else if ( words(2) == 'NITIN' ) then
                    !
                    !.md<com>    - NITINOL: Shape memory alloy
                    !.md<com>Corresponds to file 160.
                    !.md<com>The inputs are the following:
                    !.md<com><tt>E, nu</tt>
                    !.md<com>
                    !
                    messa = adjustr(trim(messa))// ' NITINOL '
                    lawco_sld(imate)=1_ip
                    kfl_sdvar_sld = 1_ip
                    lawst_sld(imate)=160_ip   ! isotropic linear elastic
                    parco_sld(1:ncoef_sld,imate)=param(2:ncoef_sld+1)

                    
                 else if ( words(2) == 'MOONE' ) then
                    messa = adjustr(trim(messa))// ' MOONEY-RIVLIN '
                    lawco_sld(imate)=1_ip
                    lawst_sld(imate)=103_ip   ! Mooney-Rivlin material (formulation Xiao and Belystcho by default)
                    parco_sld(1:ncoef_sld,imate)=param(2:ncoef_sld+1)
                    if(words(3)=='MR1') then
                       lawmo_sld(imate)=1_ip
                       parco_sld(1:ncoef_sld,imate)=param(3:ncoef_sld+2)
                    else if(words(3)=='MR2') then
                       lawmo_sld(imate)=2_ip
                       parco_sld(1:ncoef_sld,imate)=param(3:ncoef_sld+2)
                    else if(words(3)=='MR3') then
                       lawmo_sld(imate)=3_ip
                       parco_sld(1:ncoef_sld,imate)=param(3:ncoef_sld+2)
                    end if

                 else if ( words(2) == 'ORTLI' .or. words(2) == 'ORTHO' ) then
                    !
                    !.md<com>    - ORTHO: Orthotropic Saint-Venant Kirchoff model. Corresponds to file 151.
                    !.md<com>The inputs are the following:
                    !.md<com><tt>E11, E22, E33, nu12, nu13, nu23, G12, G13, G23</tt>
                    !.md<com><tt>alpha11, alpha22, alpha33</tt>
                    !.md<com>
                    !
                    messa = adjustr(trim(messa))// ' ORTHOTROPIC ELASTIC '
                    lawco_sld(imate)=1_ip
                    lawst_sld(imate)=151_ip   ! orthotropic linear elastic
                    nposi = 2_ip
                    !
                    ! Strain measure
                    if ( exists('INFIN') ) then
                       nposi = nposi + 1_ip
                       kfl_strai_sld = SLD_INFINITESIMAL
                    else if ( exists('GREEN') ) then
                       nposi = nposi + 1_ip
                       kfl_strai_sld = SLD_GREEN
                    end if
                    !
                    ! Read material properties
                    parco_sld(1:ncoef_sld,imate)=param(nposi:ncoef_sld+nposi-1)
                    !
                    ! Compute stiff0
                    call sm151_precalculus(imate)
                    !
                    ! Deactivate thermal effects
                    if( .not. kfl_therm_sld ) parco_sld(10:12,imate) = 0.0_rp

                 else if ( words(2) == 'MAIMI' ) then
                    !
                    !.md<com>    - MAIMI: 3D damage model for transversally isotropic materials. Corresponds to file 152.
                    !.md<com>The inputs are the following:
                    !.md<com><tt>E11, E22, G12, nu12, nu23,
                    !.md<com>xT,  fxT,  xC,  fxC,  muT, muLT,
                    !.md<com>yT,   yC,  sL,   nT,   nS,  nqT, nqS
                    !.md<com>gLT, fgT, gLC,  fgC,  g1T,  g2T, g2L
                    !.md<com>a11, a22, b11,  b22,   dM,  eta</tt>
                    !.md<com>
                    !
                    messa = adjustr(trim(messa))// ' ORTHOTROPIC DAMAGE MAIMI 2018'
                    lawco_sld(imate)=1_ip
                    lawst_sld(imate)=152_ip
                    kfl_damag_sld = 1_ip
                    kfl_sdvar_sld = 1_ip
                    nposi = 2_ip
                    !
                    ! Strain measure
                    if (words(3) == 'INFIN') then
                       nposi = nposi + 1_ip
                       kfl_strai_sld = SLD_INFINITESIMAL
                    else if (words(3) == 'GREEN') then
                       nposi = nposi + 1_ip
                       kfl_strai_sld = SLD_GREEN
                    else
                       call runend('SLD_REAPHY: STRAIN MEASURE FOR SM152 IS NOT DEFINED')
                    end if
                    !
                    ! Material tangent computation
                    if (exists('NUMER')) then
                       nposi = nposi + 1_ip
                       lawta_sld(imate)=SLD_TANGENT_NUMERICAL
                    else if (exists('SECAN')) then
                       nposi = nposi + 1_ip
                       lawta_sld(imate)=SLD_SECANT
                    else if (exists('ANALY')) then
                       nposi = nposi + 1_ip
                       lawta_sld(imate)=SLD_TANGENT_ANALYTICAL
                    else
                       lawta_sld(imate)=SLD_TANGENT_ANALYTICAL
                    end if
                    !
                    ! Read material properties
                    parco_sld(1:ncoef_sld,imate)=param(nposi:ncoef_sld+nposi-1)
                    !
                    ! Compute stiff0 without damage
                    call sm152_precalculus(imate)

                 else if ( words(2) == 'OLIVE' ) then
                    !
                    ! Isotropic damage model (Oliver 1990)
                    !
                    messa = adjustr(trim(messa))// ' DAMAGE OLIVER 1990'
                    lawco_sld(imate)=1_ip
                    lawst_sld(imate)=153_ip
                    kfl_damag_sld = 1_ip
                    kfl_sdvar_sld = 1_ip
                    parco_sld(1:ncoef_sld,imate)=param(2:ncoef_sld+1)
                    !
                    ! Compute stiff0 without damage
                    call sm153_precalculus(imate)

                 else if ( words(2) == 'BESSA' ) then
                    !
                    !.md<com>    - BESSA: Transversally isotropic damage model assuming plane stress.
                    !.md<com>Corresponds to file 154.
                    !.md<com>The inputs are the following:
                    !.md<com><tt>E11, E22, nu12, nu23, G12,
                    !.md<com>XT, fXT, XC, fXC, YT, YC, SL,
                    !.md<com>alpha0,
                    !.md<com>G1T, fGT, G1C, fGC, G2T, G2C, GSL,
                    !.md<com>s12p, Kp,
                    !.md<com>eta,
                    !.md<com>a11, a22, dT, b11, b22, dM
                    !.md<com>d1Max, d2Max, d3Max, d4Max, d5Max, d6Max</tt>
                    !.md<com>thick, muGSL, yBT, yBC, mufxc, E1c</tt>
                    !.md<com>
                    !
                    messa = adjustr(trim(messa))// ' ORTHOTROPIC DAMAGE BESSA 2019'
                    lawco_sld(imate)=1_ip
                    lawst_sld(imate)=154_ip
                    kfl_damag_sld = 1_ip
                    kfl_sdvar_sld = 1_ip
                    nposi = 2_ip
                    !
                    ! Strain measure
                    if ( exists('INFIN') ) then
                       nposi = nposi + 1_ip
                       kfl_strai_sld = SLD_INFINITESIMAL
                    else if ( exists('GREEN') ) then
                       nposi = nposi + 1_ip
                       kfl_strai_sld = SLD_GREEN
                    end if
                    !
                    ! Read material properties
                    parco_sld(1:ncoef_sld,imate)=param(nposi:ncoef_sld+nposi-1)
                    !
                    ! Compute stiff0 without damage
                    call sm154_precalculus(imate)
                    !
                    ! Deactivate thermal effects
                    if( .not. kfl_therm_sld ) parco_sld(24:29,imate) = 0.0_rp

                 else if ( words(2) == 'PLAST' ) then
                    !
                    ! Isotropic Hardening Plasticity
                    !
                    kfl_sdvar_sld = 1_ip
                    kfl_plast_sld = 1_ip
                    lawco_sld(imate) = 1_ip
                    lawst_sld(imate) = 0_ip
                    nposi = 2_ip
                    if ( words(3) == 'BISO ' ) then
                       lawpl_sld(imate) = 1_ip
                       nposi = nposi + 1_ip
                    end if
                    !
                    ! Strain measure
                    if ( exists('INFIN') ) then
                       nposi = nposi + 1_ip
                       kfl_strai_sld = SLD_INFINITESIMAL
                    else if ( exists('GREEN') ) then
                       nposi = nposi + 1_ip
                       kfl_strai_sld = SLD_GREEN
                    end if
                    !
                    ! Read material properties
                    parco_sld(1:ncoef_sld,imate)=param(nposi:ncoef_sld+nposi-1)

                 else if ( words(2) == 'CSHEL' ) then
                    !
                    ! Continuum shell elements
                    !
                    messa = adjustr(trim(messa))// ' CSHELL'
                    lawst_sld(imate)=200_ip
                    kfl_sdvar_sld = 1_ip
                    nposi = 2_ip
                    if ( words(3) == 'ISOLI' .or. words(2) == 'ISOTR' ) then
                       messa = adjustr(trim(messa))// ' ISOLINEAR '
                       nposi = nposi + 1_ip
                       lawco_sld(imate)=1_ip
                    else if ( words(3) == 'ORTLI' .or. words(3) == 'ORTHO' ) then
                       messa = adjustr(trim(messa))// ' ORTHOTROPIC '
                       nposi = nposi + 1_ip
                       lawco_sld(imate)=2_ip
                    else if ( words(3) == 'MAIMI' ) then
                       messa = adjustr(trim(messa))// ' ORTHOTROPIC DAMAGE MAIMI 2018'
                       nposi = nposi + 1_ip
                       lawco_sld(imate)=3_ip
                       kfl_damag_sld = 1_ip
                    else if ( words(3) == 'NEOHO' ) then
                       messa = adjustr(trim(messa))// ' NEO-HOOKEAN '
                       nposi = nposi + 1_ip
                       lawco_sld(imate)=4_ip
                    end if
                    !
                    ! Read material properties
                    parco_sld(1:ncoef_sld,imate)=param(nposi:ncoef_sld+nposi-1)
                    !
                    ! Compute stiff0 only for ORTHO
                    if ( words(3) == 'ORTHO' .or. words(3) == 'ORTLI' .or. words(3) == 'MAIMI' ) call sm151_precalculus(imate)

                 else if ( words(2) == 'SPRIN' ) then
                    lawco_sld(imate)=1_ip
                    lawst_sld(imate)=105_ip   ! constant spring
                    parco_sld(1:ncoef_sld,imate)=param(2:ncoef_sld+1)

                 else if ( words(2) == 'UCSAN' ) then
                    lawco_sld(imate)=1_ip
                    lawst_sld(imate)=135_ip
                    parco_sld(1:ncoef_sld,imate)=param(2:ncoef_sld+1)

                 else if ( words(2) == 'MICNC' ) then
                    lawco_sld(imate)=1_ip
                    lawst_sld(imate)=MAT_MICRO_NO_COUPLING
                    parco_sld(1:ncoef_sld,imate)=param(2:ncoef_sld+1)

                 else if ( words(2) == 'MICOW' ) then
                    lawco_sld(imate)=1_ip
                    lawst_sld(imate)=MAT_MICRO_ONE_WAY
                    parco_sld(1:ncoef_sld,imate)=param(2:ncoef_sld+1)

                 else if ( words(2) == 'MICFU' ) then
                    lawco_sld(imate)=1_ip
                    lawst_sld(imate)=MAT_MICRO_FULL
                    parco_sld(1:ncoef_sld,imate)=param(2:ncoef_sld+1)

                 else if ( words(2) == 'HOLZA' ) then
                    messa = adjustr(trim(messa))//' HOLZAPFEL'

                    lawco_sld(imate)=1_ip
                    lawho_sld(imate)=1_ip
                    lawst_sld(imate)=134_ip
                    parco_sld(1:ncoef_sld,imate)=param(2:ncoef_sld+1)
                    if (nnpar < 11_ip) parco_sld(11,imate)= 1.0_rp   ! required for retrocompatibility
                    ! Compute poisson

                    K0 = parco_sld(2,imate)
                    E = parco_sld(1,imate)
                    poiss = (3.0_rp*K0-E)/(6.0_rp*K0)

                    ! HGO for arteries
                    if(words(3)=='HGOAR') then
                       lawst_sld(imate)=137_ip
                       parco_sld(1:ncoef_sld,imate)=param(2:ncoef_sld+1)
                       ! compressible and isotropy
                       if (exists('COISO')) then
                          lawho_sld(imate)=4_ip
                          lawst_sld(imate)=137_ip
                          ! compressible and anisotropy
                       else if (exists('COANI')) then
                          lawho_sld(imate)=5_ip
                          lawst_sld(imate)=137_ip
                       end if
                    end if

                    messa = adjustr(trim(messa))// '--> POISSON='//trim(retost(poiss))

                 else if ( words(2) == 'GUCCI' ) then
                    lawco_sld(imate)=1_ip
                    lawst_sld(imate)=136_ip
                    parco_sld(1:ncoef_sld,imate)=param(2:ncoef_sld+1)

                 else if ( words(2) == 'ISOTH' ) then
                    messa = adjustr(trim(messa))// ' ISOLINEAR-THERMAL '
                    lawco_sld(imate)=1_ip
                    lawst_sld(imate)=170_ip
                    parco_sld(1:ncoef_sld,imate)=param(2:ncoef_sld+1)
                    if( .not. kfl_therm_sld ) parco_sld(3,imate) = 0.0_rp
                    
                 end if

                 call livinf(0_ip,messa,1_ip)

              else if ( words(1) == 'COHES' ) then
                 !
                 !.md<2>COHESIVE_MODEL:                                $ Cohesive model
                 !.md<field>COHESIVE_MODEL
                 !.md<com>This option is to define a cohesive zone model. The available laws
                 !.md<com>are the following:
                 !.md<com>
                 !
                 kfl_cohes_sld = 1_ip
                 kfl_sdvar_sld = 1_ip
                 nposi = 2_ip! Auxiliar variable to indicate at which position the parameters start
                 messa = '           COHESIVE MODEL: '
                 
                 if ( words(2) == 'DECRE' ) then
                    lawch_sld(imate)=900_ip

                 else if ( words(2) == 'SMITH' ) then
                    lawch_sld(imate)=901_ip

                 else if ( words(2) == 'EXPON' ) then
                    lawch_sld(imate)=902_ip

                 else if ( words(2) == 'POWER' ) then
                    lawch_sld(imate)=903_ip

                 else if ( words(2) == 'TURON' ) then

                    if ( exists('ORIGI') ) then
                       !
                       !.md<com>    - TURON, ORIGINAL: Interface cohesive law Turon et al. 2006-2010.
                       !.md<com>The inputs are the following:
                       !.md<com><tt>GIc, GIIc, tauI, Kp, eta, rho, dini, dmax</tt>
                       !.md<com>
                       !
                       nposi = nposi + 1_ip
                       lawch_sld(imate)=ELINT_LAW_TURON      ! idlaw = 904
                       messa = adjustr(trim(messa))// ' TURON 2006 '
                       
                    else if ( exists('CURRE') ) then
                       !
                       !.md<com>    - TURON, CURRENT: Interface cohesive law Turon et al. 2018.
                       !.md<com>The inputs are the following:
                       !.md<com><tt>GIc, GIIc, tauI, tauII, eta, Kp, rho, dini, dmax</tt>
                       !.md<com>
                       !
                       nposi = nposi + 1_ip
                       lawch_sld(imate)=ELINT_LAW_TURON_2018 ! idlaw = 905
                       messa = adjustr(trim(messa))// ' TURON 2018'
                       
                    else
                       ! Defaults
                       lawch_sld(imate)=ELINT_LAW_TURON
                       messa = adjustr(trim(messa))// ' TURON 2006 '
                       
                    end if
                    !
                    !.md<com><ul>
                    !.md<com>Calculation method of the tangent stiffness matrix. This option can be used only for Implicit analysis.
                    !.md<com>By default, the tangent stiffness matrix is calculated analytically.
                    !.md<com><ul>
                    !.md<com><li> Set to ANALYTICAL to use the analytical calculation of the stiffness matrix. </li>
                    !.md<com><li> Set to NUMERICAL  to calculate it numerically. </li>
                    !.md<com><li> Set to SECANT to use the secant method. </li>
                    !.md<com></ul>
                    !.md<com></ul>
                    !.md<com>
                    !
                    if ( exists('NUMER') ) then
                       nposi = nposi + 1_ip
                       lawta_sld(imate)=ELINT_TANGENT_NUMERICAL
                    else if ( exists('SECAN') ) then
                       nposi = nposi + 1_ip
                       lawta_sld(imate)=ELINT_TANGENT_SECANT
                    else if ( exists('ANALY') ) then
                       nposi = nposi + 1_ip
                       lawta_sld(imate)=ELINT_TANGENT_ANALYTICAL
                    else
                       ! Defaults
                       lawta_sld(imate)=ELINT_TANGENT_ANALYTICAL
                    end if

                 end if
                 !
                 ! Read material properties
                 parch_sld(1:ncoef_sld,imate)=param(nposi:ncoef_sld+nposi-1)

                 call livinf(0_ip,messa,1_ip)
                 
              else if ( words(1) == 'CONTA' ) then
                 !
                 !.md<2>CONTACT_MODEL:                                 $ Interface contact model
                 !.md<field>CONTACT_MODEL
                 !.md<com>This option is to define contact elements between instances. The available contact methods
                 !.md<com>are the following:
                 !.md<com>
                 !
                 kfl_cohes_sld = 1_ip
                 kfl_sdvar_sld = 1_ip
                 nposi = 2_ip ! Auxiliar variable to indicate at which position the parameters start
                 messa = '           CONTACT MODEL: '
                 
                 if ( words(2) == 'PENAL' ) then
                    !
                    !.md<com>    - PENALTY: Penalty method used by Turon et al. 2018
                    !.md<com>Corresponds to contact law 800. The inputs are the following:
                    !.md<com><tt>Kn, Kres</tt>
                    !.md<com>
                    !
                    lawch_sld(imate)=ELINT_LAW_CONTACT       ! idlaw = 800
                    messa = adjustr(trim(messa))// ' PENALTY'
                    !
                    !.md<com><ul>
                    !.md<com>Calculation method of the tangent stiffness matrix. This option can be used only for Implicit analysis.
                    !.md<com>By default, the tangent stiffness matrix is calculated analytically.
                    !.md<com><ul>
                    !.md<com><li> Set to ANALYTICAL to use the analytical calculation of the stiffness matrix. </li>
                    !.md<com><li> Set to NUMERICAL  to calculate it numerically. </li>
                    !.md<com><li> Set to SECANT to use the secant method. </li>
                    !.md<com></ul>
                    !.md<com></ul>
                    !.md<com>
                    !
                    if ( exists('NUMER') ) then
                       nposi = nposi + 1_ip
                       lawta_sld(imate)=ELINT_TANGENT_NUMERICAL
                    else if ( exists('SECAN') ) then
                       nposi = nposi + 1_ip
                       lawta_sld(imate)=ELINT_TANGENT_SECANT
                    else if ( exists('ANALY') ) then
                       nposi = nposi + 1_ip
                       lawta_sld(imate)=ELINT_TANGENT_ANALYTICAL
                    else
                       lawta_sld(imate)=ELINT_TANGENT_ANALYTICAL
                    end if

                 end if
                 !
                 ! Read properties
                 parcf_sld(1:ncoef_sld,imate)=param(nposi:ncoef_sld+nposi-1)

                 call messages_live(messa)
                 
              else if(words(1) == 'FIBER' .or. words(1) == 'ORTHO') then
                 !
                 ! Bio-Fiber model
                 !
                 call biofib_display_warning()

              else if( words(1) == 'CSYSM' ) then
                 !
                 !.md<2>CSYS_MATERIAL: FIELD=int, ORIENTATION | VECTOR $ Material coordinate system
                 !.md<field>CSYS_MATERIAL
                 !.md<com>This option is to define a material orientation using a coordinate system (element level).
                 !.md<com>The keyword FIELD=<tt>int</tt> is required in order to get the necessary data to build the material
                 !.md<com>coordinate system. If VECTOR (Default) is set the field must contain 2 or 3 vectors which define the 
                 !.md<com>the axes of the coordinate system. If ORIENTATION is set the field requires a list of orientations
                 !.md<com>and the coordinate system is created according to a user-defined method.
                 !.md<com>
                 !
                 !
                 ! Field code
                 oripa_sld(1) = getint('FIELD',1_ip,'#Field number for material orientations')

                 if ( exists('ORIEN') ) then
                    !
                    !.md<com><ul>
                    !.md<com><strong>Options (only for ORIENTATION):</strong>
                    !.md<com><ul>
                    !.md<com><li> Set TYPE=GLOBAL (Default) to use the global coordinate system as a material coordinate 
                    !.md<com>system.</li>
                    !.md<com><li> Set TYPE=USER to define a user coordinate system as a material coordinate system. In this case,
                    !.md<com>the coordinate system is defined in PARAMETERS section.</li>
                    !.md<com><li> Set TYPE=ELEMENT_NORMAL to define a material coordinate system using the element normal.
                    !.md<com>system.</li> In order to create this system a reference axis from the global coordinate system is
                    !.md<com>required using REFAXIS=<tt>int</tt>.
                    !.md<com>The default axis is 1 (x-direction).</li>
                    !.md<com><li> Use ROTAXIS=<tt>int1</tt> to select the rotation axis. By default, the rotation axis is 3 or the 
                    !.md<com>element normal.</li>
                    !.md<com></ul>
                    !.md<com></ul>
                    !
                    !
                    ! Method using ply angle orientations
                    !
                    if ( exists('TYPE ') ) then
                       tcsys = getcha('TYPE ','     ','#Type of coordinate system')
                       if ( tcsys == 'GLOBA' ) then         ! Using global CSYS
                          kfl_fiber_sld = 5_ip
                       else if ( tcsys == 'USER' ) then     ! Using user CSYS
                          kfl_fiber_sld = 6_ip
                       else if ( tcsys == 'ELEME' )  then   ! Element normal method
                          kfl_fiber_sld = 7_ip
                       end if
                    else
                       kfl_fiber_sld = 5_ip
                    end if
                    !
                    ! Optional parameters
                    if ( exists('ROTAX') ) then
                       oripa_sld(2) = getint('ROTAX',3_ip,'#Rotation axis')
                    end if
                    if ( exists('REFAX') ) then
                       oripa_sld(3) = getint('REFAX',1_ip,'#Reference axis for fiber orientation')
                    end if
                    
                 else if ( exists('VECTO') ) then
                    !
                    ! Method using vector fields
                    !
                    kfl_fiber_sld = 4_ip

                 else
                    !
                    ! Default
                    !
                    kfl_fiber_sld = 4_ip

                 end if
                 
              else if ( words(1) == 'CSHEL' ) then
                 !
                 ! Continuum shell elements (Reinoso et al. 2015)
                 !
                 ! ADOC[2]> CSHEL_PARAMETERS: EAS, ANS_SHEAR, ANS_TRAPEZOIDAL
                 ! ADOC[d]> CSHEL_PARAMETERS:
                 ! ADOC[d]> Parameters for continuum shell elements. These parameters are used
                 ! ADOC[d]> to remedey the typical locking patologies of continuum elements. EAS is not available in
                 ! ADOC[d]> Explicit analysis.
                 !
                 if ( exists('EAS  ') ) kfl_cshel_sld(1) = SLD_CSHEL_EAS
                 if ( exists('ANSSH') ) kfl_cshel_sld(2) = SLD_CSHEL_ANS_SHEAR
                 if ( exists('ANSTR') ) kfl_cshel_sld(3) = SLD_CSHEL_ANS_TRAPEZOIDAL

              else if ( words(1) == 'CONTA' ) then
                 ! Contact and friction model
                 parcf_sld(1:ncoef_sld,imate)=param(2:ncoef_sld+1)

              else if ( words(1) == 'STRES' ) then
                 ! STRESS_RESIDUAL_FIELD
                 !
                 ! ADO[2]> STRESS_RESIDUAL_FIELD: FIELD=int
                 ! ADO[d]> STRESS_RESIDUAL_FIELD:
                 ! ADO[d]> Residual initial stresses field in Voigt notation
                 if(words(2)=='FIELD') then
                    kfl_restr_sld = - getint('FIELD',-1_ip,'#Field number for the X row of the initial stress tensor')
                 end if

              else if ( words(1) == 'INDIS' ) then
                 !
                 ! ADO[2]> IN_DISPLACEMENT: FIELD=int | STRESSED
                 ! ADO[d]> IN_DISPLACEMENT:
                 ! ADO[d]> Initial displacements field.
                 ! ADO[d]> STRESSED means that the initial displacement field goes to the initial value of
                 ! ADO[d]> the displacement unknown, thus generating stress.
                 ! ADO[d]> Otherwise, when no additional tags, the displacement field goes to the
                 ! ADO[d]> adds to the coordinates (i.e. the reference system), with no stress generation.
                 !
                 if(words(2)=='FIELD') then
                    kfl_indis_sld(1) = - getint('FIELD',-1_ip,'#Field number for the initial displacements')
                    if (exists('STRES')) kfl_indis_sld(2)= 1_ip ! stressed displacement
                 end if

              else if ( words(1) == 'FORCE' ) then
                 !
                 ! ADOC[2]> FORCE: FIELD=int1                                       $ Nodal forces
                 ! ADOC[d]> FORCE:
                 ! ADOC[d]> This option is to define nodal forces. Nodal forces have to be defined according
                 ! ADOC[d]> to the global coordinate system. The parameter FIELD is set equal to the code
                 ! ADOC[d]> field <tt>int1</tt>, where nodal forces are defined.
                 ! This keyword should go to sld_reabcs because are BCs.
                 kfl_vofor_sld = getint('FIELD',1_ip,'#Field Number of the external volume force')

              else if ( words(1) == 'RAYLE' ) then
                 !
                 ! ADOC[2]> RAYLEIGH_DAMPING: RAYALPHA=real, RAYBETA=real           $ Rayleigh damping
                 ! ADOC[d]> RAYLEIGH_DAMPING:
                 ! ADOC[d]> This option is used to apply Rayleigh damping into the model. RAYALPHA=<tt>real</tt> is set equal to
                 ! ADOC[d]> the alpha value, while RAYBETA= <tt>real</tt> is set equal to a beta value.
                 !
                 kfl_dampi_sld(imate)= 1_ip
                 dampi_sld(1,imate)  = getrea('RAYAL',0.00_rp,'Rayleigh damping alpha factor')
                 dampi_sld(2,imate)  = getrea('RAYBE',0.00_rp,'Rayleigh damping beta factor')

              else if ( words(1) == 'MODUL' ) then
                 ! MODULATING FIELDS
                 !
                 ! ADOC[2]> MODULATING_FIELDS                                       $ Fields that modulate a given property
                 ! ADOC[3]> RIGIDITY                                            $ Modulating material rigidity
                 ! ADOC[3]> FIELD:  integer                                     $ Number of the modulating field
                 !
                 if (words(2) =='RIGID') then
                    kfl_moduf_sld(1) = getint('FIELD',1_ip,'#Field Number of the external volume force')
                 end if

              else if ( words(1) == 'COUPL' ) then
                 !
                 ! Coupling model (SOLIDZ parameters for EXMEDI coupling)
                 !
                 call eccou_show_warning()

              else if ( words(1) == 'ISOCH' ) then
                 kfl_isoch_sld=1_ip
                 if(words(2)=='DIISO') then !Displacement isochrones
                    thiso_sld(1)=param(2)
                 elseif(words(2)=='VMISO') then
                    thiso_sld(2)=param(2)
                 endif

              end if

              call ecoute('sld_reaphy')

           end do
           !
           !.md<1>END_PROPERTIES
           !
           !
           ! Final compatibility corrections
           !
           if (imate < nmate_sld) then
              ! When the total number of materials used by this module
              ! is smaller than the kernel materials:
              ! 1. if the nmatm == 1, this is the default one, so assign all properties to it
              ! 2. otherwise if the nmatm > 1 stop the problem and demand to reassign properties
              !    for this module
              nmatm = imate
              if (nmatm > 1) then
                 call runend('SLD_REAPHY: THIS MODULE HAS MORE THAN ONE MATERIAL, BUT LESS THAN THOSE AT DOM.DAT')
              end if

              do imate= nmatm, nmate_sld
                 densi_sld(1:ncoef_sld,imate) = densi_sld(1:ncoef_sld,1)
                 kfl_dampi_sld(imate)         = kfl_dampi_sld(1)
                 dampi_sld(1,imate)           = dampi_sld(1,1)
                 dampi_sld(2,imate)           = dampi_sld(2,1)
                 lawco_sld(imate)             = lawco_sld(1)
                 lawst_sld(imate)             = lawst_sld(1)
                 parco_sld(1:ncoef_sld,imate) = parco_sld(1:ncoef_sld,1)
                 parcf_sld(1:ncoef_sld,imate) = parcf_sld(1:ncoef_sld,1)

                 lawch_sld(imate)             = lawch_sld(1)
                 parch_sld(1:ncoef_sld,imate) = parch_sld(1:ncoef_sld,1)

                 lawpl_sld(imate)             = lawpl_sld(1)

                 modfi_sld(imate)             = modfi_sld(1)
                 modor_sld(1,imate)           = modor_sld(1,1)
                 modor_sld(2,imate)           = modor_sld(2,1)

              end do
           end if

        else if ( words(1) == 'PARAM' ) then

           !-------------------------------------------------------------
           !
           ! Parameters
           !
           !-------------------------------------------------------------
           !
           !.md<1>PARAMETERS
           !.md<com>
           !.md<com><strong>Parameters</strong>
           !
           call sld_memphy(2_ip)
           call ecoute('sld_reaphy')
           do while( words(1)/='ENDPA' )

              if( words(1) == 'COORD' ) then
                 !
                 !.md<2>COORDINATE_SYSTEM:                             $ User material coordinate system
                 !.md<field>COORDINATE_SYSTEM
                 !.md<com>User coordinate system for material orientation. This coordinate system is defined at
                 !.md<com>element level.
                 !.md<com>
                 !.md<com>    - Set BASIS=CARTESIAN (default) <tt>int1,int2,int3,int4,int5,int6,int7,int8,int9</tt>
                 !.md<com>a coordinate system by three points. The first point <tt>c(int1,int2,int3)</tt> corresponds to the center 
                 !.md<com>of the coordinate system;
                 !.md<com>The second point  <tt>a(int4,int5,int6)</tt> has to form the primary axis P; The last point
                 !.md<com>has to form the secondary axis S <tt>a(int7,int8,int9)</tt>. The normal axis N is calculated internally 
                 !.md<com>doing the vectorial product of N = P x S.
                 !.md<com>    - Set BASIS=CYLINDRICAL to define a cylindrical coordinate system
                 !.md<com>
                 !
                 ! Type of coordinate system
                 !
                 if( words(2) == 'BASIS' ) then
                    nposi = 2_ip
                    wcsys = getcha('BASIS','     ','#Type of material coordinate system')
                    if ( wcsys == 'CARTE' ) then
                       nposi = nposi + 1_ip
                       kfl_csysm_sld = SLD_CSYS_CARTESIAN
                    else if ( wcsys == 'CYLIN' ) then
                       nposi = nposi + 1_ip
                       kfl_csysm_sld = SLD_CSYS_CYLINDRICAL
                       call runend("SLD_REAPHY: CYLINDRICAL COORDINATE SYSTEM NOT IMPLEMENTED")
                    else
                       kfl_csysm_sld = SLD_CSYS_CARTESIAN
                    end if
                    !
                    ! Read parameters for local CSYS method
                    csysm_sld(1:ndime*3) = param(nposi:ndime*3+nposi-1)
                 end if

              else if( words(1) == 'CSYSM' ) then
                 !
                 !.md<2>CSYS_MATERIAL: FIELD=int, ORIENTATION | VECTOR $ Material coordinate system
                 !.md<field>CSYS_MATERIAL
                 !.md<com>This option is to define a material orientation using a coordinate system (element level).
                 !.md<com>The keyword FIELD=<tt>int</tt> is required in order to get the necessary data to build the material
                 !.md<com>coordinate system. If VECTOR (Default) is set the field must contain 2 or 3 vectors which define the 
                 !.md<com>the axes of the coordinate system. If ORIENTATION is set the field requires a list of orientations
                 !.md<com>and the coordinate system is created according to a user-defined method.
                 !.md<com>
                 !
                 !
                 ! Field code
                 oripa_sld(1) = getint('FIELD',1_ip,'#Field number for material orientations')
                 if( exists('ORIEN') ) then
                    !
                    !.md<com><ul>
                    !.md<com><strong>Options (only for ORIENTATION):</strong>
                    !.md<com><ul>
                    !.md<com><li> Set TYPE=GLOBAL (Default) to use the global coordinate system as a material coordinate 
                    !.md<com>system.</li>
                    !.md<com><li> Set TYPE=USER to define a user coordinate system as a material coordinate system. In this case,
                    !.md<com>the coordinate system is defined in PARAMETERS section.</li>
                    !.md<com><li> Set TYPE=ELEMENT_NORMAL to define a material coordinate system using the element normal.
                    !.md<com>system.</li> In order to create this system a reference axis from the global coordinate system is
                    !.md<com>required using REFAXIS=<tt>int</tt>.
                    !.md<com>The default axis is 1 (x-direction).</li>
                    !.md<com><li> Use ROTAXIS=<tt>int1</tt> to select the rotation axis. By default, the rotation axis is 3 or the 
                    !.md<com>element normal.</li>
                    !.md<com></ul>
                    !.md<com></ul>
                    !
                    ! Method using ply angle orientations
                    !
                    if ( exists('TYPE ') ) then
                       tcsys = getcha('TYPE ','     ','#Type of coordinate system')
                       if ( tcsys == 'GLOBA' ) then         ! Using global CSYS
                          kfl_fiber_sld = 5_ip
                       else if ( tcsys == 'USER' ) then     ! Using user CSYS
                          kfl_fiber_sld = 6_ip
                       else if ( tcsys == 'ELEME' )  then   ! Element normal method
                          kfl_fiber_sld = 7_ip
                       end if
                    else
                       kfl_fiber_sld = 5_ip
                    end if
                    !
                    ! Optional parameters
                    if ( exists('ROTAX') ) then
                       oripa_sld(2) = getint('ROTAX',3_ip,'#Rotation axis')
                    end if
                    if ( exists('REFAX') ) then
                       oripa_sld(3) = getint('REFAX',1_ip,'#Reference axis for fiber orientation')
                    end if

                 else if( exists('VECTO') ) then
                    !
                    ! Method using vector fields
                    !
                    kfl_fiber_sld = 4_ip

                 else
                    !
                    ! Default
                    !
                    kfl_fiber_sld = 4_ip

                 end if

              else if( words(1) == 'ROTAT' ) then
                 !
                 ! Material (Experimental)
                 !
                 if( words(2)=='ON   ' ) then
                    kfl_rmate_sld = 1_ip
                 else if ( words(2)=='OFF  ' ) then
                    kfl_rmate_sld = 0_ip
                 end if

              end if
              call ecoute('sld_reaphy')

           end do
           !
           !.md<1>END_PARAMETERS
           !
        else if ( words(1) == 'CARDI' ) then
           !.md<1>CARDIAC CYCLE
           !.md<com>
           !.md<com><strong>Parameters</strong>
           !
           call sld_cardiac_cycle_read_data()

        else if ( words(1) == 'SYSNE' ) then
           !.md<1>SYSNET (SYSTEM NETWORK) COUPLING
           !.md<com>
           !.md<com><strong>Parameters</strong>
           !
#ifndef PROPER_ELEM_PRIVATE_OFF
           call sysnet_read_data()
#else
           call runend("SLD_REAPHY: SYSNET IS AN ELEM PRIVATE MODULE. CONTACT YOUR BSC SUPERVISOR.")
#endif
        end if

     end do
     !
     !.md<0>END_PHYSICAL_PROBLEM
     !.md</code>
     !
  end if

end subroutine sld_reaphy

