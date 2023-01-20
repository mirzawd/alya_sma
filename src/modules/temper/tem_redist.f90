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
!> @brief   Redistribution of arrays
!> @details Redistribution of arrays. mainly allocated tem_memall
!> @} 
!-----------------------------------------------------------------------

subroutine tem_redist()

  use def_master
  use def_domain
  use def_kermod
  use mod_redistribution
  use mod_parall
  use def_parall
  use mod_memory
  use mod_mesh_type
  use def_temper
  use mod_output_postprocess, only : output_postprocess_check_variable_postprocess
  use mod_tem_arrays,         only : tem_arrays
  implicit none
  !
  ! Variables in tem_membcs: boundary conditions
  !
  call redistribution_array(kfl_fixno_tem,'NPOIN',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FIXNO_TEM')
  call redistribution_array(kfl_fixbo_tem,'NBOUN',POSIT=1_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FIXBO_TEM')
  call redistribution_array(bvess_tem,    'NPOIN',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='BVESS_TEM')
  call redistribution_array(bvnat_tem,    'NBOUN',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='BVNAT_TEM')
  call redistribution_array(kfl_funno_tem,'NPOIN',POSIT=1_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FUNNO_TEM')
  call redistribution_array(kfl_funbo_tem,'NBOUN',POSIT=1_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FUNBO_TEM')
  call redistribution_array(kfl_funtn_tem,'NPOIN',POSIT=1_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FUNTN_TEM')
  call redistribution_array(kfl_funtb_tem,'NBOUN',POSIT=1_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FUNTB_TEM')
  if( kfl_regim_tem == 4 ) then
     call redistribution_array(bvtem_tem, 'NPOIN',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='BVTEM_TEM')
  end if
  !
  ! Variables in memall
  !  
  call tem_arrays('REDISTRIBUTE')

end subroutine tem_redist
