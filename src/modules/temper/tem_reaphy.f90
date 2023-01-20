!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_reaphy()
  !------------------------------------------------------------------------
  !****f* Temper/tem_reaphy
  ! NAME 
  !    tem_reaphy
  ! DESCRIPTION
  !    This routine reads the physical problem definition for the
  !    temperature equation.
  ! USES
  ! USED BY
  !    tem_turnon
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_temper
  use def_domain
  use def_kermod
  use mod_ker_space_time_function
  use mod_ecoute, only : ecoute
  use def_kermod, only : gasco
  implicit none

  integer(ip) :: imate
  
  if( INOTSLAVE ) then
     !
     ! Initializations (defaults)
     !
     kfl_timei_tem = 0                                    ! Stationary flow
     kfl_advec_tem = 0                                    ! Convection is off
     kfl_joule_tem = 0                                    ! Joule effect is off
     kfl_radia_tem = 0                                    ! No radiation
     kfl_sourc_tem = 0                                    ! Sources are off
     kfl_adiab_tem = 0                                    ! Adiabatic mixing problem activation (= 1, adiabatic mixing) 
     kfl_condu_tem = 1                                    ! Conductivity is on
     kfl_exint_tem = 0                                    ! No properties interpolation
     kfl_inter_tem = 0                                    ! No interpolation of arrays
     kfl_dynco_tem = 0                                    ! Dynamical coupling
     kfl_regim_tem = 0                                    ! Regime
     kfl_parti_tem = 0                                    ! no particles in suspension
     kfl_skews_tem = 0                                    ! Skew symmetric convective term
     kfl_nudgi_tem = 0                                    ! nudging term to temper 
     kfl_forest_tem = 0                                   ! radiation in forest model
     turbu_tem     = 0.0_rp                               ! Turbulence parameters
     prthe_tem     = 0.0_rp                               ! Thermodynamic pressure (if NASTIN & CHEMIC not activated)
     prtur_tem     = 1.0_rp                               ! Turbulent Prandtl number Prt = 0
     scond_tem     = 0.0_rp                               ! S conductivity
     react_tem     = 0.0_rp                               ! Reaction term
     sourc_tem     = 1.0_rp                               ! Factor to scale the source term
     nudgi_tem     = 0.0_rp                               ! nudging coefficient in seconds. Term:(tempe -teref)/nudgi_tem
     !
     ! Reach the section
     !
     call ecoute('tem_reaphy')
     do while(words(1)/='PHYSI')
        call ecoute('tem_reaphy')
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDPH')
        call ecoute('tem_reaphy')
        if(words(1)=='PROBL') then
           !
           ! Problem definition data
           !
           call ecoute('tem_reaphy')
           do while(words(1)/='ENDPR')

              if(words(1)=='TEMPO') then                    ! Temporal evolution
                 if(words(2)=='ON   ') then
                    kfl_timei_tem = 1
                 else
                    kfl_timei_tem = 0
                 end if

              else if(words(1)=='REGIM') then              
                 if(words(2)=='INCOM') then
                    kfl_regim_tem=0                         ! Incompressible
                 else if(words(2)=='COMPR') then
                    kfl_regim_tem=1                         ! Compressible
                    if(exists('PRESS')) kfl_regim_tem=1
                    if(exists('DENSI')) kfl_regim_tem=2
                 else if(words(2)=='LOWMA') then
                    kfl_regim_tem = 3                       ! Low-Mach
                    if (words(3)=='ENTHA') then
                       kfl_regim_tem = 4                    ! Enthalpy equation
                       if(exists('ADIAB')) then             ! Input file --> REGIME:  LOWMACH, ENTHALPY, ADIABATIC: HMAX = 1000.0, HMIN = -300.0
                         cfi_hmax_tem = getrea('HMAX ',0.0_rp,'#Maximum enthalpy')
                         cfi_hmin_tem = getrea('HMIN ',0.0_rp,'#Minimum enthalpy')
                         kfl_adiab_tem = 1
                       end if
                    endif
                 end if

              else if(words(1)=='CONDU') then               ! Conductivity term
                 if(words(2)=='ON  ') then
                    kfl_condu_tem = 1 
                 else
                    kfl_condu_tem = 0
                 end if

              else if(words(1)=='DYNAM') then               ! Dynamical coupling
                 if(words(2)=='ON   ') kfl_dynco_tem=1

              else if(words(1)=='CONVE') then               ! Convective term
                 if(exists('ON   ')) then
                    kfl_advec_tem = 1
                    if(exists('FUNCT')) &
                         kfl_advec_tem = getint('FUNCT',0_ip,'#Velocity function')
                    if(words(3)=='VELOC') then
                       if(words(4)=='FUNCT') then
                          kfl_advec_tem = getint('FUNCT',0_ip,'#Velocity function')
                       else
                          kfl_advec_tem = 1                       
                       end if
                    else if(words(3)=='GRADI') then
                       kfl_advec_tem = -1                       
                    end if
                 end if
                 if(exists('SKEWS')) &   ! skew symmetric convective term
                      kfl_skews_tem = 1                

              else if( words(1) == 'SOURC' ) then               ! Source term

                 if( words(2) == 'MATER' ) then
                    !
                    ! Source by material
                    !
                    kfl_sourc_tem = SOURCE_TERM_MATERIAL
                    call ecoute('nsi_reaphy')
                    do while( words(1) /= 'ENDSO' )
                       if( words(1) == 'MATER' ) then
                          imate = getint('MATER',1_ip,'#Material force number')
                          if( imate < 1 .or. imate > nmate ) call runend('TEM_REAPHY: WRONG MATERIAL NUMBER')
                          if (exists('CONST')) then
                             lsour_material_tem(imate) = 2
                          else
                             call runend('TEM_REAPHY: WRONG MATERIAL FORCE MODEL')
                          end if
                          xsour_material_tem(1:msour_material_tem,imate) = param(4:3+msour_material_tem)
                       end if
                       call ecoute('tem_reaphy')                       
                    end do 
                    
                 else if( exists('SPACE') ) then
                    !
                    ! Space time source term
                    !
                    kfl_sourc_tem = SOURCE_TERM_SPACE_TIME
                    kfl_sonum_tem = space_time_function_number(getcha('SPACE','NONE ','#Space time function'))
                    
                 else if( exists('FIELD') ) then
                    !
                    ! Field source term
                    !
                    if( exists('NODE ') .or. exists('NODAL') ) then
                       kfl_sourc_tem = SOURCE_TERM_NODAL_FIELD
                       kfl_sonum_tem = getint('FIELD',1_ip,'#Element field number')
                    else
                       kfl_sourc_tem = SOURCE_TERM_FIELD
                       kfl_sonum_tem = getint('FIELD',1_ip,'#Element field number')
                    end if
                    
                 else if( exists('SPARE') ) then
                    !
                    ! Spare mesh source term
                    !
                    kfl_sourc_tem = SOURCE_TERM_SPARE_MESH
                    kfl_sonum_tem = getint('SPARE',1_ip,'#Spare mesh number')
                    
                 end if
                 
                 if( exists('JOULE') ) kfl_joule_tem = 1
                 if( exists('FACTO') ) sourc_tem = getrea('FACTO',1.0_rp,'#Source multiplicative factor')
                 
              else if(words(1)=='RADIA') then               ! Radiation
                 if(words(2)=='SURFA') kfl_radia_tem = 1

              else if(words(1)=='REACT') then               
                 !
                 ! Reaction term (s)
                 !
                 react_tem = getrea('REACT',0.0_rp,'#Reaction parameter')
              else if(words(1)=='NUDGI') then ! nudging with WRF, only if using tendencies
                 kfl_nudgi_tem = 1_ip 
                 nudgi_tem = getrea('NUDGI',0.0_rp,'#Nudging time (s.)')
                 nudgi_tem = 1.0_rp/ nudgi_tem  ! then coeff multiplies    
              else if(words(1)=='FORES') then ! nudging with WRF, only if using tendencies
                 if (exists('ON   ') ) kfl_forest_tem = 1_ip !#Radiation in forest on
                
              end if

              call ecoute('tem_reaphy')
           end do

        else if(words(1)=='PROPE') then
           !
           ! Allocate memory
           !  

           call ecoute('tem_reaphy')
           do while(words(1)/='ENDPR')
              if(words(1)=='REACT') then                    ! Reaction term (s)
                 react_tem = getrea('REACT',0.0_rp,'#Reaction parameter')

              else if(words(1)=='TURBU') then               ! Turbulent Prandtl number
                 prtur_tem = getrea('TURBU',1.0_rp,'#Turbulent Prandtl number')

              else if(words(1)=='THERM') then               ! Thermodynamic pressure
                 prthe_tem = getrea('THERM',1.0_rp,'#Thermodynamic pressure')

              else if(words(1)=='GASCO') then               ! Thermodynamic pressure
                 gasco = getrea('GASCO',1.0_rp,'#Gas constant')

              end if
              call ecoute('tem_reaphy')
           end do
        end if
     end do

  end if

end subroutine tem_reaphy
