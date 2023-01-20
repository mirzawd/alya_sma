!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_redist.f90
!> @author  houzeaux
!> @date    2019-06-17
!> @brief   Redistribution of arrays
!> @details Redistribution of arrays. mainly allocated sld_memall
!> @}
!-----------------------------------------------------------------------

subroutine sld_redist()

  use def_kintyp,             only : ip
  use def_master,             only : mem_modul, modul
  use mod_redistribution,     only : redistribution_array
  use mod_memory,             only : memory_alloca_min
  use mod_memory,             only : memory_alloca
  use mod_output_postprocess, only : output_postprocess_check_variable_postprocess
  use mod_sld_arrays,         only : sld_arrays
  use def_solidz,             only : kfl_fixno_sld, kfl_fixbo_sld
  use def_solidz,             only : kfl_fixrs_sld
  use def_solidz,             only : kfl_funno_sld, kfl_funbo_sld
  use def_solidz,             only : bvess_sld, bvnat_sld
  use def_solidz,             only : jacrot_du_dq_sld, jacrot_dq_du_sld
  use def_solidz,             only : axis1_sld,axis2_sld,axis3_sld
  use def_solidz,             only : kfl_funtn_sld,kfl_funtb_sld
  use def_solidz,             only : kfl_contn_stent
  use mod_sld_cardiac_cycle,  only : cardiac_cycle_redist

  implicit none

  !----------------------------------------------------------------------
  !
  ! Variables in sld_memphy
  !
  !----------------------------------------------------------------------

  call redistribution_array(axis1_sld,'NELEM',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='AXIS1_SLD')
  call redistribution_array(axis2_sld,'NELEM',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='AXIS2_SLD')
  call redistribution_array(axis3_sld,'NELEM',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='AXIS3_SLD')

  !----------------------------------------------------------------------
  !
  ! Variables in sld_membcs
  !
  !----------------------------------------------------------------------

  call redistribution_array(kfl_fixno_sld,'NPOIN',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FIXNO_SLD')
  call redistribution_array(kfl_fixbo_sld,'NBOUN',           MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FIXBO_SLD')
  call redistribution_array(kfl_fixrs_sld,'NPOIN',           MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FIXRS_SLD')

  call redistribution_array(bvess_sld,    'NPOIN',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='BVESS_SLD')
  call redistribution_array(bvnat_sld,    'NBOUN',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='BVNAT_SLD')

  call redistribution_array(kfl_funno_sld,'NPOIN',           MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FUNNO_SLD')
  call redistribution_array(kfl_funbo_sld,'NBOUN',           MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FUNBO_SLD')
  call redistribution_array(kfl_funtn_sld,'NPOIN',           MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FUNTN_SLD')
  call redistribution_array(kfl_funtb_sld,'NBOUN',           MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FUNTB_SLD')

  call redistribution_array(kfl_contn_stent,'NPOIN',           MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_CONTN_STENT')

  call redistribution_array(jacrot_du_dq_sld,'NPOIN',POSIT=3_ip, MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='JACROT_DU_DQ_SLD')
  call redistribution_array(jacrot_dq_du_sld,'NPOIN',POSIT=3_ip, MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='JACROT_DQ_DU_SLD')

  !----------------------------------------------------------------------
  !
  ! Variables in sld_memall
  !
  !----------------------------------------------------------------------

  call sld_arrays('REDISTRIBUTE')

  !----------------------------------------------------------------------
  !
  ! Some necessary extra-things
  !
  !----------------------------------------------------------------------
  call cardiac_cycle_redist()

end subroutine sld_redist
