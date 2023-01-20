!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    gus_reaphy.f90
!> @author  Guillaume Houzeaux
!> @brief   Read physical data
!> @details Read physical data
!> @} 
!-----------------------------------------------------------------------
subroutine gus_reaphy()
!-----------------------------------------------------------------------
 !.md<module>gusano
 !.md<input>case.gus.dat
 !.md<pos>0
 !.md<sec>

  use def_parame
  use def_inpout
  use def_master
  use def_gusano
  use def_domain
  use mod_ecoute, only : ecoute
  use mod_ecoute, only : ecoute_reach_section
  use mod_physics
  
  implicit none

  if( INOTSLAVE ) then
     !
     ! Initializations (defaults)
     !
     kfl_timei_gus = 0                 ! Time integration off
     kfl_regim_gus = LAMINAR_PIPE      ! Regime to determine velocity profile and friction
     kfl_bendi_gus = GUS_BEND_OFF      ! Bending model
     kfl_conve_gus = 1                 ! Convective term
     !
     ! Reach the section
     !
     rewind(lisda)
     call ecoute_reach_section('PHYSI')
     !--><group>
     !-->    <groupName>PHYSICAL_PROBLEM</groupName>
     !
     !.md# Physical Properties Definition
     !.md<code>
     !.md<0>PHYSICAL_PROBLEM
     !
     do while(words(1)/='ENDPH')

        call ecoute('gus_reaphy')

        if(words(1) == 'PROBL' ) then

           call ecoute('gus_reaphy')

           do while(words(1)/='ENDPR')
              !--><subGroup>
              !-->     <subGroupName>PROBLEM_DEFINITION</subGroupName>           
              !
              !.md<1>PROBLEM_DEFINITION
              ! 
              if( words(1) == 'CONVE' ) then

                 if( words(2) == 'OFF  ' ) then
                    kfl_conve_gus = 0
                 end if
              
              else if( words(1) == 'TEMPO' ) then          
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
                 !.md<com>is enhanced when passing through a stationary regime to reach a steady state..\f$ x^2 \frac{1}{2} \f$...
                 !
                 if(exists('ON   ')) kfl_timei_gus = 1

              else if(words(1)=='REGIM') then
                 !
                 ! Regime
                 !
                 if(      exists('LAMIN') .and. exists('PIPE ') ) then
                    kfl_regim_gus = LAMINAR_PIPE
                 else if( exists('TURBU') .and. exists('PIPE ') ) then
                    kfl_regim_gus = TURBULENT_PIPE
                 else if( exists('SMOOT') .and. exists('PIPE ') ) then
                    kfl_regim_gus = SMOOTH_PIPE
                 end if
                 
              else if(words(1)=='BENDI') then
                 !
                 ! Bending model
                 !
                 if(      exists('LAMIN') .and. exists('PIPE ') ) then
                    kfl_bendi_gus = LAMINAR_PIPE
                 !else if( exists('TURBU') .and. exists('PIPE ') ) then
                 !   kfl_bendi_gus = TURBULENT_PIPE
                 !else if( exists('SMOOT') .and. exists('PIPE ') ) then
                 !   kfl_bendi_gus = SMOOTH_PIPE
                 end if
                 
              end if
              call ecoute('gus_reaphy')
           end do
        end if
        !--></subGroup>
        !
        !.md<1>END_PROPERTIES
        !
     end do
     !--></group>
     !
     !.md<0>END_PHYSICAL_PROBLEM
     !.md</code>

  end if
end subroutine gus_reaphy
