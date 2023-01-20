!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ker_iniunk()
  !-----------------------------------------------------------------------
  !****f* Kermod/ker_iniunk
  ! NAME 
  !    ker_iniunk
  ! DESCRIPTION
  !    This routine sets up the initial condition for the velocity.
  !    If this is a restart, initial condition are loaded somewhere else.
  !    Values are stored in position
  !    VELOC(:,:,NPREV_NSI) 
  !    PRESS(:,:,NPREV_NSI) 
  ! USED BY
  !    ker_begste
  !***
  !-----------------------------------------------------------------------
  use def_parame 
  use def_master
  use def_domain
  use def_kermod
  use mod_ker_proper
  use mod_eccoupling,    only : kfl_exmsld_ecc, eccou_set_state_variables
  use mod_biofibers,     only : kfl_biofibers, biofibers
  use mod_ker_updpro,    only : ker_updpro
  implicit none
  !
  ! Compute wall distance and wall normal
  !
  call ker_walgen(1_ip)
  call ker_walnor(1_ip)
  !
  ! Calculate roughness
  !
  call ker_roughn()
  !
  ! Calculate canopy height & height over terrain
  !
  call ker_canopy()
  !
  ! Velocity and temperature field
  !
  call ker_velfun(ITASK_INIUNK)
  call ker_temfun(ITASK_INIUNK)
  call ker_confun(ITASK_INIUNK)
  call ker_disfun(ITASK_INIUNK)
  call ker_arefun(ITASK_INIUNK)
  ! 
  ! Properties
  !  
  call ker_proper_allocate_properties()
  call ker_proper_check_properties()     
  call ker_updpro(ITASK_INIUNK)
  !
  ! Test properties
  !
  !call ker_tespro()
  !call ker_waexlo()
  !
  ! Elector-mechanical coupling
  !
  if( kfl_exmsld_ecc ) call eccou_set_state_variables()
  !
  ! Bio-fibers at the reference configuration
  !
  if( kfl_biofibers ) call biofibers % update_fibers_at_nodes(ITASK_INIUNK) 
  
end subroutine ker_iniunk
