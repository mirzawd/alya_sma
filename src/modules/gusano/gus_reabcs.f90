!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup GusanoInput
!> @ingroup    Gusano
!> @{
!> @file    gus_reabcs.f90
!> @author  Guillaume Houzeaux
!> @brief   Read boundary conditions 
!> @details Read boundary conditions, initial conditions and parameters
!> @} 
!-----------------------------------------------------------------------
subroutine gus_reabcs()
  !-----------------------------------------------------------------------
  !.md<module>gusano
  !.md<input>case.gus.dat
  !.md<pos>3
  !.md<sec>
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_gusano 
  use mod_memchk
  use mod_opebcs
  use def_kermod
  use mod_ecoute, only : ecoute
  use mod_ecoute, only : ecoute_reach_section
  implicit none
  !
  ! Allocate memory
  !
  if( kfl_icodn > 0 ) then
     call opebcs_initialization_structure(1_ip,tncod_gus)
     call opebcs_initialization_variable (2_ip,tncod_gus(1))  ! 1D velocity and pressure
  end if
  if( kfl_icodb > 0 ) then
     call opebcs_initialization_structure(1_ip,tbcod_gus)     ! Momentum 
     call opebcs_initialization_variable (1_ip,tbcod_gus)
  end if
 
  if( INOTSLAVE ) then
     !
     ! Initialization global variables
     !
     !
     ! Reach the nodal-wise section
     !
     call ecoute_reach_section('BOUND')

     call ecoute('gus_reabcs')
     do while(words(1)/='ENDBO')

        if( words(1) == 'CODES' .and. exists('NODES') ) then

           !-------------------------------------------------------------
           !
           ! User-defined codes on nodes
           !
           !-------------------------------------------------------------

           tncod => tncod_gus(1:)
           call boundary_conditions_read_node_codes('VELOCITY')


        else if( words(1) == 'CODES' .and. exists('BOUND') ) then

           !-------------------------------------------------------------
           !
           ! User-defined codes on boundaries
           !          
           !-------------------------------------------------------------

           if( kfl_icodb <= 0 ) then
              call runend('GUS_REABCS: BOUNDARY CONDITIONS ON BOUNDARIES SHOULD BE DECLARED IN *.DOM.DAT FILE')
           else
              tbcod => tbcod_gus(1:)
              call boundary_conditions_read_boundary_codes('MOMENTUM')
           end if
           
        end if

        call ecoute('gus_reabcs')

     end do
     !
     !
  end if

end subroutine gus_reabcs
