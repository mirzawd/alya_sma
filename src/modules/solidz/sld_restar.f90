!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_restar.f90
!> @author  Solidz Team
!> @date    September, 2017
!>          Adds state dependent variables for any model
!> @brief   This routine reads/writes values required for a restart
!> @details Displacement is  always read and written. Velocity and
!>          and acceleration only for dynamic problems.
!>
!>          \verbatim
!>          ITASK = 1 ... Reads the initial values from the restart file
!>                  2 ... Writes restart file
!>          \endverbatim
!> @}
!------------------------------------------------------------------------

subroutine sld_restar(itask)

  use def_kintyp,             only : ip,rp
  use def_master,             only : ITASK_READ_RESTART,ITASK_WRITE_RESTART
  use def_master,             only : ITER_K,displ
  use mod_memory,             only : memory_size
  use def_domain,             only : ndime,npoin
  use mod_restart,            only : restart_ini
  use mod_restart,            only : restart_end
  use mod_restart,            only : restart_add
  use mod_strings,            only : integer_to_string
  use def_solidz,             only : allie_sld,allwk_sld,allke_sld
  use def_solidz,             only : kfl_newbc_sld,nprev_sld
  use def_solidz,             only : veloc_sld,accel_sld
  use def_solidz,             only : kfl_rigid_sld
  use def_solidz,             only : grnor_sld
  use mod_sld_arrays,         only : sld_arrays
  use mod_sld_cardiac_cycle,  only : kfl_cardiac_cycle,sld_cardiac_cycle_manage_restart
  use mod_sld_rbo,            only : sld_rbo_restar
#ifndef PROPER_ELEM_PRIVATE_OFF
  use mod_ker_sysnet, only: kfl_sysnet, sysnet_manage_restart
#endif

  implicit none

  integer(ip), intent(in)    :: itask                    !< What to do
  integer(ip)                :: itime,ipoin

  !----------------------------------------------------------------------
  !
  ! Primary arrays
  !
  !----------------------------------------------------------------------

  if(      itask == ITASK_READ_RESTART ) then
     call sld_arrays('READ RESTART') 
  else if( itask == ITASK_WRITE_RESTART ) then
     call sld_arrays('WRITE RESTART')
  end if

  !----------------------------------------------------------------------
  !
  ! Variables
  !
  !----------------------------------------------------------------------

  call restart_ini(itask)

  call restart_add(grnor_sld,'grnor_sld')
  if( kfl_rigid_sld == 0 ) then
     call restart_add(allie_sld,'allie_sld')
     call restart_add(allwk_sld,'allwk_sld')
     call restart_add(allke_sld,'allke_sld')
  else
     call sld_rbo_restar()
  end if

  if( kfl_cardiac_cycle ) then
     call sld_cardiac_cycle_manage_restart()
  endif

#ifndef PROPER_ELEM_PRIVATE_OFF
  if( kfl_sysnet ) then
     call sysnet_manage_restart()
  endif
#endif

  call restart_end(itask)

  !----------------------------------------------------------------------
  !
  ! New boundary conditions
  !
  !----------------------------------------------------------------------
  
  if( kfl_newbc_sld == 1 .and. itask == ITASK_READ_RESTART ) then
     do itime = nprev_sld,1,-1
        do ipoin = 1,npoin
           displ(    1:ndime,ipoin,itime) = 0.0_rp
           veloc_sld(1:ndime,ipoin,itime) = 0.0_rp
           accel_sld(1:ndime,ipoin,itime) = 0.0_rp
        end do
     end do
  end if

end subroutine sld_restar
