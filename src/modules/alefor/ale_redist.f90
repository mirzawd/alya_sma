!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Alefor
!> @{
!> @file    ale_redist.f90
!> @author  houzeaux
!> @date    2019-06-17
!> @brief   Redistribution of arrays
!> @details Redistribution of arrays. mainly allocated ale_memall
!> @} 
!-----------------------------------------------------------------------

subroutine ale_redist()

  use def_kintyp,         only : ip
  use def_master,         only : modul
  use def_master,         only : mem_modul
  use def_master,         only : kfl_fixno_ale
  use def_master,         only : bvess_ale
  use mod_redistribution, only : redistribution_array
  use mod_ale_arrays,     only : ale_arrays
  use def_alefor,         only : kfl_fixbo_ale
  use def_alefor,         only : kfl_fixrs_ale
  use def_alefor,         only : kfl_funno_ale
  use def_alefor,         only : kfl_funtn_ale
  use def_alefor,         only : kfl_funbo_ale

  implicit none
  !
  ! Variables in ale_membcs: boundary conditions
  !
  call redistribution_array(kfl_fixno_ale,'NPOIN',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FIXNO_ALE')
  call redistribution_array(kfl_fixbo_ale,'NBOUN',POSIT=1_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FIXBO_ALE')
  call redistribution_array(kfl_fixrs_ale,'NPOIN',POSIT=1_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FIXRS_ALE')
  call redistribution_array(bvess_ale,    'NPOIN',POSIT=2_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='BVESS_ALE')
  
  call redistribution_array(kfl_funno_ale,'NPOIN',POSIT=1_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FUNNO_ALE')
  call redistribution_array(kfl_funtn_ale,'NPOIN',POSIT=1_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FUNTN_ALE')
  call redistribution_array(kfl_funbo_ale,'NBOUN',POSIT=1_ip,MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KFL_FUNBO_ALE')
  !
  ! Variables in memall
  !  
  call ale_arrays('REDISTRIBUTE')
  
end subroutine ale_redist
