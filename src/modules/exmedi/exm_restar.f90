!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_restar.f90
!> @date    04(03/2020
!> @author  Guillaume Houzeaux
!> @brief   Restart subroutine
!> @details Restart subroutine
!> @}
!------------------------------------------------------------------------
subroutine exm_restar(itask)
  
   use def_kintyp,             only : ip
   use def_master,             only : ITASK_READ_RESTART
   use def_master,             only : ITASK_WRITE_RESTART
   use mod_messages,           only : livinf
   use mod_exm_arrays,         only : exm_arrays
   use mod_restart,            only : restart_ini, restart_add, restart_end
   use mod_memory,             only : memory_size

   implicit none
 
   integer(ip), intent(in) :: itask 
 
   !----------------------------------------------------------------------
   !
   ! Primary arrays
   !
   !----------------------------------------------------------------------
 
   if( itask == ITASK_READ_RESTART ) then
      call exm_arrays('READ RESTART')
   else if( itask == ITASK_WRITE_RESTART ) then
      call exm_arrays('WRITE RESTART')
   end if

   !----------------------------------------------------------------------
   !
   ! Variables
   !
   !----------------------------------------------------------------------   

   !
   ! Do not save aptim to restart!
   ! It is very unlikely that in normal circumstances the restart will be used in [aplap_exm(istim), aplap_exm(istim)+ aptim(istim)]
   ! However not putting it in restart gives a flexibility to change activation points after the restart
   ! The latter is a more common scenario, since many times we will need to run a problem to steady state and then change activation
   !
   !nstim_exm_tmp = nstim_exm
   !
   !call restart_ini(itask)   
   !call restart_add(nstim_exm_tmp,'nstim_exm')
   !call restart_end(itask)
   !
   !if ( nstim_exm_tmp .ne. nstim_exm ) then
   !   call runend("EXMEDI: NUMBER OF STIMULI IS DIFFERENT IN .DAT AND IN RESTART")
   !end if
   !
   !if ( nstim_exm>0_ip ) then 
   !   call restart_ini(itask)
   !   call restart_add(aptim,'aptim')
   !   call restart_end(itask)
   !end if

end subroutine exm_restar
 
