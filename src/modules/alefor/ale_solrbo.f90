!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_solrbo()
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_solrbo
  ! NAME
  !    ale_solrbo
  ! DESCRIPTION
  !    This routines solves the rigid body motion ! similar to ibm_solite
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_alefor
  use mod_messages, only : livinf

  implicit none
  !
  ! Obtain the initial guess for inner iterations - This is here instead of nsi_begite 
  ! because nsi_begite is used for the basic ale - mesh deform
  !
  call ale_updunk(ITASK_BEGINN)
  !
  ! Add forces: gravity, magnetic: F^n+1 and T^n+1
  !     
  call ale_forces()
  !
  ! Obtain accel, veloc and displacements 
  !
  if ( kfl_crist_ale == 1 ) then
     call ale_sorbcr(dtime)

  else if( kfl_genco_ale == 1_ip )then
     call ale_sorbjc(dtime)   

  else
     call ale_sorbhh(dtime)

  end if
  !
  call ale_fixirb()

end subroutine ale_solrbo
