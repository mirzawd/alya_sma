!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Chemic
!> @{
!> @file    chm_redist.f90
!> @author  houzeaux
!> @date    2019-06-17
!> @brief   Redistribution of arrays
!> @details Redistribution of arrays. mainly allocated chm_memall
!> @}
!-----------------------------------------------------------------------

subroutine chm_redist()

  use def_master,             only : mem_modul, modul
  use mod_redistribution,     only : redistribution_array
  use def_kintyp,             only : ip
  use def_chemic,             only : bvess_chm, kfl_fixno_chm, kfl_fixbo_chm, kfl_funno_chm, kfl_funtn_chm
  use mod_chm_arrays,         only : chm_arrays
  implicit none
  !
  ! Variables in chm_membcs: boundary conditions
  !
  call redistribution_array(kfl_fixno_chm,'NPOIN',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FIXNO_CHM')
  call redistribution_array(kfl_fixbo_chm,'NBOUN',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FIXBO_CHM')
  call redistribution_array(bvess_chm,    'NPOIN',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='BVESS_CHM')
  call redistribution_array(kfl_funno_chm,'NPOIN',POSIT=1_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FUNNO_CHM')
  call redistribution_array(kfl_funtn_chm,'NPOIN',POSIT=1_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FUNTN_CHM')
  !
  ! Variables in memall
  !
  call chm_arrays('REDISTRIBUTE')

end subroutine chm_redist
