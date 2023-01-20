!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_chkpoi(ichec)
  !------------------------------------------------------------------------
  !****f* Parall/par_chkpoi
  ! NAME
  !    par_chkpoi
  ! DESCRIPTION
  !    Checkpoint. Retunrs ICHEC:
  !    ICHEC = 0 ... MPI REDUCE has performed well
  !
  !    
  ! OUTPUT
  !    ICHEC
  ! USED BY
  !    par_partit
  !***
  !------------------------------------------------------------------------
  use def_master
  use def_parall
  use mod_communications, only : PAR_SUM
  implicit none  
  integer(ip), intent(out) :: ichec
  integer(ip)              :: dummi
  !
  ! Checkpoint
  !
  if( kfl_paral>=0) then
     dummi = 1
     call PAR_SUM(dummi)
     ichec=nproc_par-dummi-1
  end if

end subroutine par_chkpoi

