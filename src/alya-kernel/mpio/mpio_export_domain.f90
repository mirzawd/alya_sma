!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup mpio
!> @{
!> @file    mpio_export_domain.f90
!> @author  houzeaux
!> @date    2019-05-03
!> @brief   Export domain
!> @details Export domain with MPIO
!> @} 
!-----------------------------------------------------------------------

subroutine mpio_export_domain()

  use def_kintyp,      only : ip
  use def_kermod,      only : kfl_oumes
  use def_kermod,      only : kfl_posdi
  use mod_mesh_type,   only : mesh_type_allocate_initialize
  use mod_mesh_type,   only : mesh_type_update_last_mesh
  use mod_mesh_type,   only : mesh_type_allocate_minimum
  use mod_messages,    only : messages_live
  use mod_output,      only : output_domain
  use def_domain,      only : npoin_own
  use def_master,      only : npoi3
  use mod_mpio_config, only : mpio_config
  implicit none
  
  if( mpio_config%output%post_process%export_only ) then

     call messages_live('MESH EXPORT IN MPIO FORMAT','START SECTION')
     
     kfl_oumes = 1                                           ! All mesh
     kfl_posdi = 0                                           ! Required by MPIO
     npoin_own = npoi3

     call mesh_type_allocate_initialize(0_ip)                ! MESHE
     call mesh_type_update_last_mesh(0_ip)                   ! MESHE(0)

     call output_domain(CURRENT_MESH=0_ip,ONLY_MESH=.true.)  ! Output domain

     call messages_live('MESH EXPORT IN MPIO FORMAT','END SECTION')
     call runend('O.K.!')

  end if

end subroutine mpio_export_domain
