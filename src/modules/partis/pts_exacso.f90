!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine pts_exacso(itask,t,a,u,x)

  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime
  use def_partis, only     :  kfl_exacs_pts
  implicit none
  integer(ip), intent(in)  :: itask
  real(rp),    intent(in)  :: t
  real(rp),    intent(out) :: a(ndime)
  real(rp),    intent(out) :: u(ndime)
  real(rp),    intent(out) :: x(ndime)
  real(rp)                 :: x_exa(3),u_exa(3),a_exa(3)

  a_exa = 0.0_rp
  u_exa = 0.0_rp
  x_exa = 0.0_rp

  select case ( kfl_exacs_pts ) 

  case ( 0_ip ) 

     return

  case ( 1_ip ) 

     a_exa(1) = -0.01_rp/(4.0_rp) * sin(t/2.0_rp)
     u_exa(1) =  0.01_rp/(2.0_rp) * cos(t/2.0_rp)
     x_exa(1) =  0.01_rp          * sin(t/2.0_rp) + 0.5_rp

  case ( 2_ip ) 

     a_exa(1) =  0.02_rp
     u_exa(1) =  0.02_rp*t
     x_exa(1) =  0.01_rp*t*t + 0.5_rp

  end select

  if( itask == 1 ) then
     !
     ! Initial coordinate
     !
     x(1:ndime) = x_exa(1:ndime)

  else if( itask == 2 ) then
     !
     ! Initial solution
     !
     x(1:ndime) = x_exa(1:ndime)
     u(1:ndime) = u_exa(1:ndime)
     a(1:ndime) = a_exa(1:ndime)

  else if( itask == 3 ) then
     !
     ! Acceleration 
     !
     a(1:ndime) = a_exa(1:ndime)

  end if

end subroutine pts_exacso
