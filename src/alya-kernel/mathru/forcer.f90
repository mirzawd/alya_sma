!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine forcer()

  !------------------------------------------------------------------------
  !
  ! This function results in an error so that alya stops and the traceback tells where you are
  !
  !------------------------------------------------------------------------

  use def_kintyp, only     :  ip,rp
  implicit none

  real(rp)                 :: hhhhh

  call cputim(hhhhh)

  hhhhh = sin(hhhhh)
  hhhhh = sqrt(hhhhh-5.0)
  
end subroutine forcer
 
