!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_turnof()
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_turnof
  ! NAME
  !    chm_turnof
  ! DESCRIPTION
  !    This routine closes the run for the chemic module
  ! USES
  !    chm_outcpu
  ! USED BY
  !    Chemic
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only : ip
  implicit none

  external :: chm_openfi
  external :: chm_outcpu

  !
  ! Close files
  !
  call chm_openfi(4_ip)

  !
  ! Output CPU times
  !
  call chm_outcpu()

end subroutine chm_turnof

