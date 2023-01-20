!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Nastin 
!> @{
!> @file    nsi_plugin.f90
!> @date    14/10/2014
!> @author  Guillaume Houzeaux
!> @brief   Initialize coupling
!> @details Initialize some arrays for the coupling:
!>          1. Modify fixity array in order to impose UNKNOWN type
!>             interfaces automatically.
!>          2. Check all dofs are not prescribed on RESIDUAL type
!>             interfaces.
!>          3. Use keyword FORCE to force the imposition.
!> @}
!------------------------------------------------------------------------
subroutine nsi_plugin_init()
  !
  ! Mandatory variables => 
  ! 
  use def_kintyp,    only :  ip,rp
  use def_master,    only :  modul
  use mod_couplings, only :  COU_SET_FIXITY_ON_TARGET
  use mod_couplings, only :  COU_PUT_VALUE_ON_TARGET
  !
  ! <= Mandatory variables
  ! 
  use def_nastin,    only :  kfl_fixno_nsi
  use def_nastin,    only :  kfl_fixpp_nsi
  use def_nastin,    only :  kfl_fixpr_nsi
  implicit none

  call COU_SET_FIXITY_ON_TARGET('VELOC',modul,kfl_fixno_nsi)
  call COU_SET_FIXITY_ON_TARGET('PRESS',modul,kfl_fixpp_nsi)

  call COU_SET_FIXITY_ON_TARGET('UNKNO',modul,kfl_fixno_nsi)
  call COU_SET_FIXITY_ON_TARGET('UNKNO',modul,kfl_fixpp_nsi)

  call COU_SET_FIXITY_ON_TARGET('RESID',modul,kfl_fixno_nsi)
  call COU_SET_FIXITY_ON_TARGET('MOMEN',modul,kfl_fixno_nsi)

  call COU_SET_FIXITY_ON_TARGET('RESID',modul,kfl_fixpp_nsi)
  call COU_SET_FIXITY_ON_TARGET('CONTI',modul,kfl_fixpp_nsi)

  call COU_SET_FIXITY_ON_TARGET('RESID',modul,kfl_fixpr_nsi,'FREE BETWEEN SUBDOMAINS')
  call COU_SET_FIXITY_ON_TARGET('CONTI',modul,kfl_fixpr_nsi)
  !
  ! Free pressure on wet nodes
  !
  !call COU_PUT_VALUE_ON_TARGET(0_ip,kfl_fixpr_nsi)
  !call COU_PUT_VALUE_ON_TARGET(1_ip,kfl_fixpr_nsi,3_ip)

end subroutine nsi_plugin_init 
