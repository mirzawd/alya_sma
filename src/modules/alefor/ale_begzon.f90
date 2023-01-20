!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_begzon
  !-----------------------------------------------------------------------
  !*f Alefor/ale_begzon
  ! NAME 
  !    ale_begste
  ! DESCRIPTION
  !    This routine prepares for a new coupling iteration of the ALE formulation
  !    equation      
  ! USES
  !    ale_updunk
  ! USED BY
  !    Alefor
  !*
  !-----------------------------------------------------------------------
  use def_master,      only : ITASK_BEGZON,iblok,modul
  use def_coupli,      only : kfl_gozon
  use def_kintyp,      only : ip
  use mod_moduls_conf, only : moduls_in_block
  implicit none
  !
  ! At the beginning of a coupling iteration, the coordinates  
  ! must go back to the previous time step values in FSI
  !
  if ( kfl_gozon == 1_ip .and. moduls_in_block(iblok,modul) ) call ale_updunk( ITASK_BEGZON )

end subroutine ale_begzon
