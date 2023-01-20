!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_timste
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_timste
  ! NAME 
  !    nsi_timste
  ! DESCRIPTION
  !    This routine computes the time step
  !    Initial solution is here because memory for matrix and RHS must 
  !    have been previously allocated by master
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
  implicit none 
  integer(ip), save :: miinn_old
  !
  ! Actualizes safety factor
  ! ittim is last time step number 
  ! Do not modify safet at first time step
  !  
  if( ittim > 0 ) safet_nsi = min(safeo_nsi*((safex_nsi)**ittim), safma_nsi)
  ! 
  ! Time step size
  !
  if( kfl_stead_nsi /= 1 ) call nsi_updtss()
  !
  ! Start inner iterations at time step kfl_stain_nsi
  !
  if( kfl_stain_nsi > 0 .and. ittim < kfl_stain_nsi ) then
     miinn_old = momod(modul) % miinn
     momod(modul) % miinn = 1
     kfl_stain_nsi        = -kfl_stain_nsi
  else if( kfl_stain_nsi < 0 .and. ittim >= abs(kfl_stain_nsi) ) then
     momod(modul) % miinn = miinn_old
     kfl_stain_nsi        = 0
  end if

end subroutine nsi_timste
