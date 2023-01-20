!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Begste
!> @{
!> @file    Begste.f90
!> @author  Guillaume Houzeaux
!> @brief   Begin a time step
!> @details Begin a time step:
!>          - Modules update boundary conditions, unknowns, etc.
!> @} 
!-----------------------------------------------------------------------
subroutine Begste()
  use def_kintyp_basic,  only : ip
  use def_master,        only : iblok
  use def_master,        only : nblok
  use def_master,        only : itti2
  use def_master,        only : ITASK_BEGSTE
  use def_master,        only : ITASK_RECOVER_RESET
  use def_master,        only : ITASK_SAVE_RESET
  use mod_ker_proper
  use def_kermod,        only : kfl_reset
  use def_kermod,        only : kfl_reset_steps
  use def_domain,        only : kfl_domar_world
  use def_coupli,        only : mcoup
  use def_coupli,        only : coupling_driver_iteration
  use mod_couplings,     only : COU_TEMPORAL_PREDICTOR
  use mod_messages,      only : livinf
  use mod_moduls,        only : moduls 
  use mod_ker_updpro,    only : ker_updpro
  use mod_arrays,        only : arrays_recover_primary
  use mod_arrays,        only : arrays_save_primary
  implicit none 
  !
  ! Turn back reset flag to previous step
  !
  if( kfl_reset == 0 .and. mod(itti2,max(1_ip,kfl_reset_steps)) == 0 ) then
     call setgts(ITASK_SAVE_RESET)
     call arrays_save_primary()
  else if( kfl_reset == 1 .and. itti2 - kfl_reset_steps + 1 > 0 ) then
     call iniste(2_ip)
     call setgts(ITASK_RECOVER_RESET)
     call arrays_recover_primary()
     call ker_updpro()
     kfl_reset = 0
  end if
  !
  ! Update properties
  !
  !call ker_updpro(ITASK_BEGSTE)
  !
  ! Coupling
  !
  call cou_begste()
  !
  ! Begin a time step for each module
  !
  call Kermod(ITASK_BEGSTE)
  do iblok = 1,nblok
     call moduls(ITASK_BEGSTE)
  end do
  iblok = 1
  !
  ! Coupling: Put counters to zero
  ! 
  if( mcoup > 0 ) then
     coupling_driver_iteration(1:nblok) = 0
  end if
  !
  ! Moving meshes, recompute some things
  !
  if( kfl_domar_world == 1 ) call domarr(3_ip)

end subroutine Begste
