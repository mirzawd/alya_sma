!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_begste()
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_begste
  ! NAME 
  !    nsi_begste
  ! DESCRIPTION
  !    This routine prepares for a new time step of the incompressible NS
  !    equations.
  ! USES
  !    nsi_iniunk
  !    nsi_updtss
  !    nsi_updbcs
  !    nsi_updunk
  ! USED BY
  !    Nastin
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  use mod_nsi_low_Mach,      only : nsi_low_Mach_upd_mass_rho
  use mod_nsi_multi_step_fs, only : nsi_multi_step_fs_matrices
  use mod_nsi_multi_step_fs, only : nsi_multi_step_fs_solution
  use def_master,            only : times

  implicit none

  if( kfl_stead_nsi /= 1 )  then 
     !
     ! Update frame of reference
     !
     call nsi_updfor()
     !
     ! Update boundary conditions
     !
     call nsi_updbcs(ITASK_BEGSTE) 
     call nsi_updunk(ITASK_BEGSTE) ! (:,2) <= (:,1)
     !
     ! Initialize some variables before starting the time step
     !
     call nsi_inivar(3_ip)
     !
     ! Coupling with dynamic solver
     !
     call nsi_dyncou(1_ip)
     !
     ! Update mass_rho term for low-Mach, for drho/dt calculation
     !
     if( associated(mass_rho_nsi) ) &
          call nsi_low_Mach_upd_mass_rho()
     !
     ! Moving subdomains, recompute matrices
     !
     call times(8) % ini()
     if( associated(velom) .and. NSI_FRACTIONAL_STEP ) then
        call nsi_multi_step_fs_matrices()
     end if
     call times(8) % add()
     
  end if 
  !
  ! Surface tension
  !
  if( kfl_surte_nsi /= 0 ) call nsi_norcur

  call nsi_coupli(ITASK_BEGSTE)
  !
  ! Initialize coupling
  !
  call nsi_plugin_init()
  
end subroutine nsi_begste
