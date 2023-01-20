!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_redist.f90
!> @author  houzeaux
!> @date    2019-06-17
!> @brief   Redistribution of arrays
!> @details Redistribution of arrays. mainly allocated exm_memall
!> @} 
!-----------------------------------------------------------------------

subroutine exm_redist()

  use def_master, only : ip, mem_modul, modul
  use mod_redistribution
  use def_exmedi
  use mod_exm_arrays,         only : exm_arrays
  implicit none
  !
  ! Variables in exm_membcs: boundary conditions
  !
  !
  ! Variables in memall
  !
  call exm_arrays('REDISTRIBUTE')
  !
  ! Variables in exm_solmem: ???
  !
  call redistribution_array(vdiag_exm, 'NPOIN', POSIT=1_ip, MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='VDIAG_EXM') 
  call redistribution_array(idima_exm, 'NPOIN', POSIT=1_ip, MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='IDIMA_EXM') 
  call redistribution_array(ticel_exm, 'NPOIN', POSIT=1_ip, MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='TICEL_EXM') 
  call redistribution_array(jicel_exm, 'NPOIN', POSIT=1_ip, MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='JICEL_EXM')
  call redistribution_array(appfi_exm, 'NPOIN', POSIT=1_ip, MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='APPFI_EXM')
  call redistribution_array(cedif_exm, 'NPOIN', POSIT=3_ip, MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='CEDIF_EXM') 
  call redistribution_array(grafi_exm, 'NPOIN', POSIT=3_ip, MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='GRAFI_EXM')
  call redistribution_array(kgrfi_exm, 'NELEM', POSIT=1_ip, MEMOR=mem_modul(1:2,modul),VARIABLE_NAME='KGRFI_EXM')


end subroutine exm_redist
