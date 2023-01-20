!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ker_restar(itask)
  !------------------------------------------------------------------------
  !****f* Nastin/ker_restar
  ! NAME 
  !    ker_restar
  ! DESCRIPTION
  !    This routine:
  !    ITASK = 1 ... Reads the initial values from the restart file
  !            2 ... Writes restart file
  ! USES
  ! USED BY
  !    ker_turnon
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_kermod
  use mod_ker_arrays,         only : ker_arrays
  use mod_messages,           only : messages_live
  use mod_eccoupling,         only : eccou_manage_restart
  use mod_restart,            only : restart_ini
  use mod_restart,            only : restart_add
  use mod_restart,            only : restart_end
  use mod_strings,            only : integer_to_string
  
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ifunc

  !----------------------------------------------------------------------
  !
  ! Primary arrays
  !
  !----------------------------------------------------------------------

  if(       itask == ITASK_READ_RESTART ) then
     call ker_arrays('READ RESTART')
  else if( itask == ITASK_WRITE_RESTART ) then
     call ker_arrays('WRITE RESTART')
  end if
  
  !----------------------------------------------------------------------
  !
  ! Variables
  !
  !----------------------------------------------------------------------
 
  call restart_ini(itask)
  !
  ! Wind-Kessel
  !
  do ifunc = 1,max_windk_systems
     if( windk_systems(ifunc) % ndxs  > 0 ) then
        call restart_add(windk_systems(ifunc) % ndxs            ,'ndx_'             //integer_to_string(ifunc))
        call restart_add(windk_systems(ifunc) % stored_time_step,'stored_time_step_'//integer_to_string(ifunc))
        call restart_add(windk_systems(ifunc) % yprev           ,'yprev_'           //integer_to_string(ifunc)) 
        call restart_add(windk_systems(ifunc) % xprev           ,'xprev_'           //integer_to_string(ifunc))
        call restart_add(windk_systems(ifunc) % y_out           ,'y_out_'           //integer_to_string(ifunc))
        call restart_add(windk_systems(ifunc) % x_in            ,'x_in_'            //integer_to_string(ifunc))
     end if
  end do
  !
  ! Electro-mechanical coupling (exm-sld)
  !
  call eccou_manage_restart(itask)
  !
  ! End restart
  !
  call restart_end(itask)
 
end subroutine ker_restar
 
