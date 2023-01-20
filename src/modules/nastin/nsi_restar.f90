!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_restar.f90
!> @author  houzeaux
!> @date    2020-05-20
!> @brief   Restart
!> @details Nastin restart variables
!> @} 
!-----------------------------------------------------------------------

subroutine nsi_restar(itask)

  use def_kintyp_basic,   only : ip
  use def_master,         only : ITASK_READ_RESTART
  use def_master,         only : ITASK_WRITE_RESTART
  use def_nastin,         only : grnor_nsi
  use def_nastin,         only : ubpre_nsi
  use def_nastin,         only : avtim_nsi
  use def_nastin,         only : avste_nsi
  use def_nastin,         only : corio_nsi
  use mod_nsi_arrays,     only : nsi_arrays
  use mod_restart,        only : restart_ini
  use mod_restart,        only : restart_end
  use mod_restart,        only : restart_add
  implicit none

  integer(ip), intent(in) :: itask
  
  !----------------------------------------------------------------------
  !
  ! Primary arrays
  !
  !----------------------------------------------------------------------
  
  if( itask == ITASK_READ_RESTART ) then
     call nsi_arrays('READ RESTART')
  else if( itask == ITASK_WRITE_RESTART ) then
     call nsi_arrays('WRITE RESTART')
  end if
  
  !----------------------------------------------------------------------
  !
  ! Variables
  !
  !----------------------------------------------------------------------

  call restart_ini(itask)
  call restart_add(grnor_nsi,'grnor_nsi')
  call restart_add(ubpre_nsi,'ubpre_nsi')
  call restart_add(avtim_nsi,'avtim_nsi')
  call restart_add(avste_nsi,'avste_nsi')
  call restart_add(corio_nsi,'corio_nsi')
  call restart_end(itask)

end subroutine nsi_restar
 
