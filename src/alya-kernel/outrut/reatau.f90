!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine reatau(kfl_taust)
  !-----------------------------------------------------------------------
  !****f* outrut/reatau
  ! NAME 
  !    reatau
  ! DESCRIPTION
  !    This routine reads the tau strategy used to compute tau in subroutines ADR_tau(new)  or/and in tauadr (old)
  ! INPUT
  ! OUTPUT
  ! USES
  ! USED BY
  !    ***_reanut
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip
  use def_inpout, only      :  words
  implicit none
  integer(ip),  intent(out) :: kfl_taust

  if(words(2)=='OFF  ') then ! No stabilization
     kfl_taust=0 
  else if(words(2)=='CODIN') then ! Codina
     kfl_taust=1
  else if(words(2)=='EXACT') then
     kfl_taust=2
  else if(words(2)=='SHAKI') then
     kfl_taust=3
  else if(words(2)=='INCLU') then ! Codina including DT in TAU
     kfl_taust=5
  else if(words(2)=='TIMES') then ! TAU =DT
     kfl_taust=6
  end if

end subroutine reatau
