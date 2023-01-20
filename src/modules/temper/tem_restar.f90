!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_restar(itask)
  !------------------------------------------------------------------------
  !****f* Temper/tem_restar
  ! NAME 
  !    tem_restar
  ! DESCRIPTION
  !    This routine:
  !    ITASK = 1 ... Reads the initial values from the restart file
  !            2 ... Writes restart file
  ! USES
  ! USED BY
  !    tem_turnon
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use def_kermod,         only : kfl_adj_prob
  use mod_tem_arrays,     only : tem_arrays
  use mod_communications, only : PAR_BROADCAST
  use mod_restart,        only : restart_ini
  use mod_restart,        only : restart_end
  use mod_restart,        only : restart_add
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,icomp
  
  !----------------------------------------------------------------------
  !
  ! Primary arrays
  !
  !----------------------------------------------------------------------

  if( itask == ITASK_READ_RESTART ) then
     call tem_arrays('READ RESTART')
  else if( itask == ITASK_WRITE_RESTART ) then
     call tem_arrays('WRITE RESTART')
  end if

  !----------------------------------------------------------------------
  !
  ! Restart file tem.rst
  !
  !----------------------------------------------------------------------

  call restart_ini(itask)
  call restart_add(vinvt_tem,'vinvt_tem')
  call restart_add(prthe,    'prthe')
  call restart_add(avtim_tem,'avtim_tem')
  call restart_end(itask)
  
  !----------------------------------------------------------------------
  !
  ! Assign constant tempe forward values for adjoint
  !
  !----------------------------------------------------------------------  
   
  if( itask == ITASK_READ_RESTART ) then
     icomp = min(3_ip,ncomp_tem)  
  else
     icomp = 1 
  end if
  if( itask == ITASK_READ_RESTART .and. kfl_adj_prob == 1 ) then
    do ipoin = 1, npoin
        tempe_forw(ipoin,1) = tempe(ipoin,icomp)    
        tempe(ipoin,icomp)  = 0.0_rp
    end do    
 endif
 
end subroutine tem_restar
 
