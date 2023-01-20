!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_reabcs()

  !-----------------------------------------------------------------------
  ! Sources/modules/magnet/mag_reabcs.f90
  ! NAME
  !    	mag_reabcs
  ! DESCRIPTION
  !    	This routine reads the boundary conditions.
  ! USES
  ! USED BY
  !   	mag_turnon
  !-----------------------------------------------------------------------
	
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_magnet
  use mod_ecoute,              only : ecoute
  use mod_ecoute,              only : ecoute_reach_section
  use mod_boundary_conditions, only : boundary_conditions_read_edge_codes
  use mod_boundary_conditions, only : boundary_conditions_read_boundary_codes
  use mod_boundary_conditions, only : boundary_conditions_initialization_structure
  use mod_boundary_conditions, only : boundary_conditions_initialization_variable

  implicit none
  !
  ! Allocate memory
  !
  !if( kfl_icode > 0 ) then
  call boundary_conditions_initialization_structure(1_ip, tecod_mag)
  call boundary_conditions_initialization_variable(ndime, tecod_mag(1)) 
  call boundary_conditions_initialization_structure(1_ip, tbcod_mag)
  call boundary_conditions_initialization_variable(ndime, tbcod_mag(1)) 
  !end if
 
  if( INOTSLAVE ) then
     !
     ! Reach the nodal-wise section
     !
     call ecoute_reach_section('BOUND','mag_reabcs')
     !
     ! Loop over edges and/or boundaries
     !
     call ecoute('mag_reabcs')
     !
     do while(words(1)/='ENDBO')

        if( words(1) == 'CODES' .and. exists('EDGES') ) then
           !
           ! User-defined codes on edges
           !
           call boundary_conditions_read_edge_codes('MAGNET', tecod_mag)
           !
        else if( words(1) == 'CODES' .and. exists('BOUND') ) then
           !
           ! User-defined codes on boundaries
           !          
           call boundary_conditions_read_boundary_codes('MAGNET', tbcod_mag)
           !
        end if
        !
        call ecoute('mag_reabcs')
        !
     end do

  end if

end subroutine mag_reabcs
