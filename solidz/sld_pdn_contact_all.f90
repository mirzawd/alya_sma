!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_pdn_contact_all.f90
!> @author  gguillam
!> @date    2022-06-02
!> @brief   PDN contact algorithm version manager
!> @details PDN contact algorithm version manager
!>          
!> @} 
!-----------------------------------------------------------------------

subroutine sld_pdn_contact_all

#if defined COMMDOM && COMMDOM == 2
  
  use def_kintyp,                 only : ip,rp
  use def_solidz,                 only : kfl_rigid_sld
  use def_solidz,                 only : kfl_conta_sld
  use def_solidz,                 only : SLD_PDN_RBO_DEFORMABLE, SLD_PDN_UNILATERAL, SLD_PDN_BILATERAL
  use mod_sld_pdn_contact,        only : commdom_sld_plugin_pdn
  use mod_sld_pdn_contact,        only : sld_rbo_contact_plugin_rbo
  use mod_sld_pdn_contact_plepp,  only : SLD_PDN_RIGID_TO_DEFOR
  use mod_sld_pdn_contact_plepp,  only : sld_pdn_contact_plugin

  implicit none
  
  if( (kfl_conta_sld == SLD_PDN_UNILATERAL .or. kfl_conta_sld == SLD_PDN_BILATERAL) .and. kfl_rigid_sld == 0_ip ) then
     !
     ! Unilateral and Bilateral (Rivero Phd 2018)
     !
     call commdom_sld_plugin_pdn()

  else if( kfl_conta_sld == SLD_PDN_RBO_DEFORMABLE ) then
     !
     ! Rigid body - Deformable
     !
     call sld_rbo_contact_plugin_rbo()

  else if( kfl_conta_sld == SLD_PDN_RIGID_TO_DEFOR ) then
     !
     ! Rigid body - Deformable (v1.1)
     !
     call sld_pdn_contact_plugin()

  end if

#endif
  
end subroutine sld_pdn_contact_all


