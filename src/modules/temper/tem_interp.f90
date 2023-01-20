!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    tem_redist.f90
!> @author  houzeaux
!> @date    2019-06-17
!> @brief   Interpribution of arrays
!> @details Interpribution of arrays. mainly allocated tem_memall
!> @} 
!-----------------------------------------------------------------------

subroutine tem_interp()

  use def_master
  use def_domain
  use def_kermod
  use mod_parall
  use def_parall
  use mod_memory
  use mod_mesh_type
  use def_temper
  use mod_output_postprocess, only : output_postprocess_check_variable_postprocess
  use mod_tem_arrays,         only : tem_arrays
  use mod_AMR_interpolate,    only : AMR_interpolate
  use def_AMR
  implicit none
  !  
  ! Variables in tem_membcs: boundary conditions
  !  
  call AMR_interpolate(kfl_fixno_tem,'NPOIN',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FIXNO_TEM')
  call AMR_interpolate(kfl_fixbo_tem,'NBOUN',POSIT=1_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FIXBO_TEM')
  call AMR_interpolate(bvess_tem,    'NPOIN',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='BVESS_TEM')
  call AMR_interpolate(bvnat_tem,    'NBOUN',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='BVNAT_TEM')
  call AMR_interpolate(kfl_funno_tem,'NPOIN',POSIT=1_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FUNNO_TEM')
  call AMR_interpolate(kfl_funbo_tem,'NBOUN',POSIT=1_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FUNBO_TEM')
  call AMR_interpolate(kfl_funtn_tem,'NPOIN',POSIT=1_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FUNTN_TEM')
  call AMR_interpolate(kfl_funtb_tem,'NBOUN',POSIT=1_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FUNTB_TEM')
  if( kfl_regim_tem == 4 ) then
     call AMR_interpolate(bvtem_tem, 'NPOIN',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='BVTEM_TEM')
  end if
  !
  ! Variables in memall
  !
  call tem_arrays('INTERPOLATE')
  
end subroutine tem_interp 
