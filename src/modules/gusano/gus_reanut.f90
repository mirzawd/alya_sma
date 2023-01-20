!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup GusanoInput
!> @{
!> @file    gus_reanut.f90
!> @author  Guillaume Houzeaux
!> @brief   Read numerical data
!> @details Read numerical data
!> @}
!-----------------------------------------------------------------------
subroutine gus_reanut()
  !.md<module>gusano
  !.md<input>case.gus.dat
  !.md<pos>1
  !.md<sec>
  use def_parame
  use def_inpout
  use def_master
  use def_gusano
  use def_domain
  use def_solver
  use mod_elmgeo, only : FLAT_BUBBLE
  use mod_elmgeo, only : QUADRATIC_BUBBLE
  use mod_elmgeo, only : FREE_SURFACE_BUBBLE
  use mod_ecoute, only : ecoute
  use mod_ecoute, only : ecoute_reach_section

  implicit none

  integer(ip) :: ivari
  
  if( INOTSLAVE ) then
     !
     !  Initializations (defaults)
     !
     kfl_stabi_gus           = GUS_ASGS
     cotol_gus               = 1.0e-2_rp
     safet_gus               = 1.0_rp
     sstol_gus               = 1.0e-5_rp                            ! Steady-satate tolerance
     kfl_algor_gus           = GUS_MONOLITHIC
     momod(modul) % miinn    = 1
     ivari                   = 1
     !
     ! Reach the section
     !
     call ecoute_reach_section('NUMER')
     !--><group>
     !-->       <groupName>Numerical_treatment</groupName>
     !
     !
     !.md# Numerical Treatment Definition
     !.md<code>
     !.md<0>NUMERICAL_TREATMENT
     !
     !-->          <subGroup>
     do while( words(1) /= 'ENDNU' )
        if( words(1) == 'ALGEB' ) then
           !
           ! Algebraic solvers
           !
           solve_sol => solve(ivari:)
           call reasol(1_ip)
           ivari = 1
           
        else if( words(1) == 'MAXIM' ) then
           !
           ! Inner iterations
           !
           momod(modul) % miinn = int(param(1))
           
        else if( words(1) == 'MOMEN' ) then
           !
           ! Momentum equation
           !
           ivari = 2
           
        else if( words(1) == 'CONTI' ) then
           !
           ! Continuity equation
           !
           ivari = 3
           
         else if( words(1) == 'CONVE' ) then
           !
           ! Convergence inner iteration tolerance
           !
           cotol_gus = param(1)
           
        else if( words(1) == 'ALGOR' ) then
           !
           ! Algorithm
           !
           if( words(2) == 'MONOL' ) then
              kfl_algor_gus = GUS_MONOLITHIC
           else if( words(2) == 'SCHUR' ) then
              kfl_algor_gus = GUS_SCHUR_COMPLEMENT
           else
              call runend('GUS_REANUT: UNKNOWN SOLUTION STRATEGY')
           end if
           
       else if( words(1) == 'STABI' ) then
           !--><inputLine>
           !-->            <inputLineName>STABILIZATION</inputLineName>
           !-->            <inputLineHelp><![CDATA[<b>STABILIZATION</b>: Determines the stabilization technique of the Navier-Stokes equations.
           !-->                    Both ASGS and FULL_OSS adds to the momentum and continuity equations a stabilization term.
           !-->                    The difference if the equation of the subgrid scale.
           !-->                    
           !-->                        -  ASGS:     Algebraic subgrid scale 
           !-->                        -  FULL_OSS: Orthogonal subgrid scale 
           !-->                        -  SPLIT:    Split orthogonal subgrid scale.
           !-->                    
           !-->                     ]]></inputLineHelp>
           !-->            <inputElement>
           !-->                <inputElementType>combo</inputElementType>
           !-->                <item><itemName>ASGS</itemName></item><item><itemName>FULL_OSS</itemName></item><item><itemDependence>0</itemDependence><itemName>SPLIT_OSS</itemName></item>
           !-->            </inputElement>
           !-->            <inputElement>
           !-->                <inputElementGroup>0</inputElementGroup>
           !-->                <inputElementType>combo</inputElementType>
           !-->                <item><itemName>EMPTY</itemName></item><item><itemName>FIRST_ORDER</itemName></item><item><itemName>SOTO</itemName></item><item><itemName>NO_LIMITER</itemName></item>
           !-->            </inputElement>
           !-->        </inputLine>
           !
           !.md<1>STABILIZATION:              ASGS | OSS | SPLIT_OSS
           !.md<field>STABILIZATION
           !.md<com> Determines the stabilization technique of the Navier-Stokes equations.
           !.md<com>Both ASGS and FULL_OSS adds to the momentum and continuity equations a stabilization term.
           !.md<com>The difference if the equation of the subgrid scale.
           !.md<com>
           !.md<com>    -  ASGS:     Algebraic subgrid scale 
           !.md<com>    -  FULL_OSS: Orthogonal subgrid scale 
           !.md<com>    -  SPLIT:    Split orthogonal subgrid scale.
           !.md<com>
           !
           if( words(2) == 'ASGS ' ) then
              kfl_stabi_gus = GUS_ASGS                                   ! ASGS

           else if( words(2) == 'FULLO' .or. words(2) == 'OSS  ' ) then
              kfl_stabi_gus = GUS_OSS                                    ! Full OSS

           else if( words(2) == 'SPLIT' ) then
              kfl_stabi_gus = GUS_SPLIT_OSS                              ! Split OSS

           end if
           
        else if( words(1) == 'SAFET' ) then
           !--><inputLine>
           !-->            <inputLineName>SAFETY_FACTOR</inputLineName>
           !-->            <inputLineHelp>When using an explicit scheme, the CFL condition determines a maximum time step to obtain a stable scheme. When using</inputLineHelp>
           !-->            <inputElement>
           !-->                <inputElementType>edit2</inputElementType>
           !-->                <inputElementValueType>REAL</inputElementValueType>
           !-->                <inputLineEditValue>5</inputLineEditValue>
           !-->            </inputElement>
           !--> </inputLine>
           !
           !.md<1>SAFETY_FACTOR=        real
           !.md<field>SAFETY_FACTOR
           !.md<com>When using an explicit scheme, the CFL condition determines a maximum time step
           !.md<com>to obtain a stable scheme. When using
           !
           safet_gus = param(1)

        else if( words(1) == 'STEAD' ) then
           !-->                <inputLine>
           !-->                    <inputLineName>STEADY_STATE_TOLERANCE</inputLineName>
           !-->                    <inputLineHelp>Tolerance to decide when the solution of the module should stop because the steady state have been reached.</inputLineHelp>
           !-->                    <inputElement>
           !-->                        <inputElementType>edit2</inputElementType>
           !-->                <inputElementValueType>REAL</inputElementValueType>
           !-->                        <inputLineEditValue>5</inputLineEditValue>
           !-->                    </inputElement>
           !-->                </inputLine>
           !
           !.md<1>STEADY_STATE_TOLERANCE=real                                                               $ Tolerance for steady state
           !.md<field>STEADY_STATE_TOLERANCE
           !.md<com>Tolerance to decide when the solution of the module should stop because the steady
           !.md<com>state have been reached.
           !
           sstol_gus = param(1)

        end if
        call ecoute('gus_reanut')
     end do
     !-->        </group>
     !
     !.md<0>END_NUMERICAL_TREATMENT
     !.md</code>
     !
  end if

end subroutine gus_reanut

