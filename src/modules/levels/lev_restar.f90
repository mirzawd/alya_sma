!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_restar(itask)
  !------------------------------------------------------------------------
  !****f* levels/lev_restar
  ! NAME 
  !    lev_restar
  ! DESCRIPTION
  !    This routine:
  !    ITASK = 1 ... Reads the initial values from the restart file
  !            2 ... Writes restart file
  ! USES
  ! USED BY
  !    lev_turnon
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_levels
  use mod_postpr
  use mod_memchk
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: icomp,iwopo,kfl_gores,ipoin
  !
  ! Check if restrt file should be read or written
  !
  call respre(itask,kfl_gores)
  if( kfl_gores == 0 ) return

  if( itask == READ_RESTART_FILE ) then
     icomp = 3
  else
     icomp = 1
  end if

  !----------------------------------------------------------------------
  !
  ! Level set
  !
  !----------------------------------------------------------------------

  iwopo =   1
  if( INOTMASTER ) gesca => fleve(:,icomp)
  call postpr(gesca,postp(1)%wopos(1:2,iwopo),ittim,cutim)
  if( INOTMASTER ) then
     do ipoin = 1,npoin
        bvess_lev(1,ipoin,1) = fleve(ipoin,icomp) 
     end do
  end if

  !----------------------------------------------------------------------
  !
  ! Finish
  !
  !----------------------------------------------------------------------

  call respre(3_ip,kfl_gores)

end subroutine lev_restar
 
