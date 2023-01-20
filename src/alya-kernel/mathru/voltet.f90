!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine voltet(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,vol)

!-----------------------------------------------------------------------
!
! Computes the volume of a tetrahedra
!
!-----------------------------------------------------------------------
  use      def_kintyp
  use      def_master, only : zeror
  implicit none
  real(rp),    intent(in)  :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
  real(rp),    intent(out) :: vol

  real(rp)                 :: norm,pscal
  real(rp)                 :: xp,yp,zp,x1p,y1p,z1p

  xp = (y2-y1)*(z3-z1)-(z2-z1)*(y3-y1)
  yp = (z2-z1)*(x3-x1)-(x2-x1)*(z3-z1)
  zp = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)

  norm = sqrt(xp*xp+yp*yp+zp*zp)

  if (abs(norm) > zeror) then

     xp = xp/norm
     yp = yp/norm
     zp = zp/norm

     x1p = x4-x1
     y1p = y4-y1
     z1p = z4-z1

     pscal = xp*x1p+yp*y1p+zp*z1p

     vol = norm*abs(pscal)/6.0_rp

  else

     vol = 0.0_rp

  end if

end subroutine voltet
