!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_turnof
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_turnof
  ! NAME 
  !    nsi_turnof
  ! DESCRIPTION
  !    This routine closes NASTIN module
  ! USES
  !    nsi_output
  !    nsi_openfi
  ! USED BY
  !    Nastin
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_solver
  use def_domain
  use def_nastin
  implicit none

  !-------------------------------------------------------------------
  !
  ! Interpolation from coarse to fine mesh
  !
  !-------------------------------------------------------------------
  !
  if(kfl_meshi_nsi /= 0_ip) call nsi_coarfine(2_ip)
  !
  ! Output miniapp
  !
  call nsi_output_miniapp()
  !
  ! llamada a Alya ADAN para cierre de archivos
  !
  call nsi_cadan(7_ip) 
  
end subroutine nsi_turnof

