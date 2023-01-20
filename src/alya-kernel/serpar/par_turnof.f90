!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_turnof
  !-----------------------------------------------------------------------
  !****f* parall/par_turnof
  ! NAME
  !    par_turnof
  ! DESCRIPTION
  !    This subroutine turns off service
  ! USES
  ! USED BY
  !    Turnof
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_parall
  use mod_parall
  use mod_outfor, only : outfor
  use mod_run_config, only : run_config
  implicit none
  !
  ! Write CPU time heading and master's CPU time
  !
  call par_outcpu()
  !
  ! Write tail for formatted files
  !
  if(IMASTER) then
     call outfor( 26_ip,lun_outpu_par,' ')
  end if
  !
  ! Close files
  !
  call par_openfi(four)

end subroutine par_turnof
