!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mixmax(slope,f,g,value)
  !-----------------------------------------------------------------------
  !****f* mathru/mixing
  ! NAME
  !   mixing
  ! DESCRIPTION
  !   Mix two functions f and g
  ! OUTPUT 
  !   WEIGH
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp 
  implicit none
  real(rp),    intent(in)  :: slope,f,g
  real(rp),    intent(out) :: value
  real(rp)                 :: H

  if(g==0.0_rp) then
     value  = f
  else if(f==0.0_rp) then
     value  = g
  else
     value  = f/g 
     H      = 0.5_rp*(tanh(slope*(value-1.0_rp))+1.0_rp)
     value  = f*H+(1.0_rp-H)*g
  end if

end subroutine mixmax
