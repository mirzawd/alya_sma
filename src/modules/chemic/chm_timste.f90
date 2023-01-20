!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_timste()
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_timste
  ! NAME
  !    chm_begste
  ! DESCRIPTION
  !    This routine computes the new time step
  ! USES
  !    chm_iniunk
  !    chm_updtss
  !    chm_updbcs
  !    chm_updunk
  !    chm_radvuf
  ! USED BY
  !    Chemic
  !***
  !-----------------------------------------------------------------------
  use def_master,  only : momod, modul
  implicit none

  external chm_updtss

  !
  ! Time step size
  !
  if(momod(modul) % kfl_stead/=1) call chm_updtss()

end subroutine chm_timste
