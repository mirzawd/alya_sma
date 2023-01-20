!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_interp.f90
!> @author  houzeaux
!> @date    2019-06-17
!> @brief   Redistribution of arrays
!> @details Redistribution of arrays. mainly allocated nsi_memall
!> @} 
!-----------------------------------------------------------------------

subroutine nsi_interp()

  use def_master
  use def_domain
  use def_kermod
  use mod_parall
  use def_parall
  use mod_memory
  use mod_mesh_type
  use def_nastin
  use mod_AMR_interpolate,    only : AMR_interpolate
  use mod_output_postprocess, only : output_postprocess_check_variable_postprocess
  use mod_nsi_arrays,         only : nsi_arrays
  implicit none
  !
  ! Recompute boundary conditions
  !
  call nsi_inibcs()

!!$  call AMR_interpolate(kfl_fixno_nsi,    'NPOIN',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FIXNO_NSI')
!!$  call AMR_interpolate(kfl_fixbo_nsi,    'NBOUN',           MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FIXBO_NSI')
!!$  call AMR_interpolate(kfl_fixrs_nsi,    'NPOIN',           MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FIXRS_NSI')
!!$  call AMR_interpolate(kfl_fixpr_nsi,    'NPOIN',           MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FIXPR_NSI')
!!$  call AMR_interpolate(kfl_fixpp_nsi,    'NPOIN',           MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FIXPP_NSI')
!!$  call AMR_interpolate(bvess_nsi,        'NPOIN',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='BVESS_NSI')
!!$  call AMR_interpolate(bvnat_nsi,        'NBOUN',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='BVNAT_NSI')
!!$  call AMR_interpolate(bpess_nsi,        'NPOIN',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='BPESS_NSI')
!!$  call AMR_interpolate(kfl_wlawf_nsi,    'NPOIN',           MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_WLAWF_NSI')
!!$
!!$  call AMR_interpolate(kfl_funno_nsi,    'NPOIN',           MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FUNNO_NSI')
!!$  call AMR_interpolate(kfl_funbo_nsi,    'NBOUN',           MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FUNBO_NSI')
!!$  call AMR_interpolate(kfl_funtn_nsi,    'NPOIN',           MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FUNTN_NSI')
!!$  call AMR_interpolate(kfl_funtb_nsi,    'NBOUN',           MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FUNTB_NSI')
!!$
!!$  call AMR_interpolate(kfl_fixno_div_nsi,'NPOIN',           MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FIXNO_DIV_NSI')
  !
  ! Variables in memall
  !
  call nsi_arrays('INTERPOLATE')

  !----------------------------------------------------------------------
  !
  ! Some necessary extra-things
  !
  !----------------------------------------------------------------------

  if( kfl_regim_nsi == 2 ) then
     unk2n_nsi => densi
  else 
     unk2n_nsi => press
  end if
  !
  ! Update boundary conditions (for example, node with pressure prescription)
  !
  call nsi_updbcs(ITASK_TURNON)
  !
  ! Allocate minimum memory to be safe
  !
  call memory_alloca_min(mem_modul(1:2,modul),'VELOC','nsi_memall',veloc)
  call memory_alloca_min(mem_modul(1:2,modul),'PRESS','nsi_memall',press)

end subroutine nsi_interp
