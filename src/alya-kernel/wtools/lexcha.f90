!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lexcha( lg_vari )
  !-----------------------------------------------------------------------
  !****f* iexcha
  ! NAME
  !    iexcha
  ! DESCRIPTION
  !    This routine exchange integer data individually
  ! USES
  ! USED BY
  !    nsi_sendat
  !    nsa_sendat
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kintyp_basic,     only :  ip, lg
  implicit none

  logical(lg), intent(inout)    ::  lg_vari
  integer(ip)                   ::  ip_vari
  
  if( lg_vari )then
        ip_vari = 1_ip
  else
        ip_vari = 0_ip
  endif
  call iexcha( ip_vari )
  if( ip_vari == 1_ip )then
        lg_vari = .True.
  else
        lg_vari = .False.
  endif
end subroutine lexcha
