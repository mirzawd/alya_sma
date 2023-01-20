!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_updfor()
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_updfor
  ! NAME 
  !    nsi_updfor
  ! DESCRIPTION
  !    This routine updates the angular and linear velocity and
  !    acceleration of the frame of reference.
  ! USES
  ! USED BY
  !    nsi_begste
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_nastin
  implicit none
  real(rp), external :: funcre,funcrd 
  real(rp)           :: fvela_old(3),facca_old(3)
  real(rp)           :: fvell_old(3),faccl_old(3)
  !
  ! Second order: Old values
  !
  if( kfl_tiacc_nsi == 2 .and. kfl_tisch_nsi == 1 ) then
     fvela_old = fvela_nsi
     facca_old = facca_nsi
     fvell_old = fvell_nsi
     faccl_old = faccl_nsi
  end if
  !
  ! Values at n+1
  !
  fvela_nsi = fvdia_nsi * fvnoa_nsi * funcre(fvpaa_nsi,6_ip,kfl_fvfua_nsi,cutim) ! w
  facca_nsi = fadia_nsi * fanoa_nsi * funcrd(fvpaa_nsi,6_ip,kfl_fvfua_nsi,cutim) ! dw/dt
  fvell_nsi = fvdil_nsi * fvnol_nsi * funcre(fvpal_nsi,6_ip,kfl_fvful_nsi,cutim) ! v
  faccl_nsi = fadil_nsi * fanol_nsi * funcrd(fvpal_nsi,6_ip,kfl_fvful_nsi,cutim) ! dv/dt
  !
  ! Second order: values at n+theta
  !
  if( kfl_tiacc_nsi == 2 .and. kfl_tisch_nsi == 1 ) then
     fvela_nsi = 0.50_rp * ( fvela_old + fvela_nsi )
     facca_nsi = 0.50_rp * ( facca_old + facca_nsi )
     fvell_nsi = 0.50_rp * ( fvell_old + fvell_nsi )
     faccl_nsi = 0.50_rp * ( faccl_old + faccl_nsi )
  end if
  !
  ! Angular velocity norm |w|
  !
  call vecnor(fvela_nsi,3_ip,corio_nsi,2_ip)

end subroutine nsi_updfor
