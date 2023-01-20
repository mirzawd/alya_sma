!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine dipopl(px,py,pz,x1,y1,z1,x2,y2,z2,x3,y3,z3,d)

!-----------------------------------------------------------------------
!
! Computes the distance between a point and a plane
!
!-----------------------------------------------------------------------
  use      def_kintyp
  implicit none
  real(rp),    intent(in)  :: px,py,pz,x1,y1,z1,x2,y2,z2,x3,y3,z3
  real(rp),    intent(out) :: d

  integer(ip)              :: sames
  real(rp)                 :: x12,y12,z12,x13,y13,z13,norm,pscal,d1
  real(rp)                 :: cp1x,cp1y,cp1z,cp2x,cp2y,cp2z
  real(rp)                 :: xp,yp,zp,x1p,y1p,z1p
  real(rp)                 :: d23,d13,d12


  ! Compute the projection of the point on the triangle plan
  x12=x2-x1
  y12=y2-y1
  z12=z2-z1

  x13=x3-x1
  y13=y3-y1
  z13=z3-z1

  xp=y12*z13-z12*y13
  yp=z12*x13-x12*z13
  zp=x12*y13-y12*x13

  norm=sqrt(xp*xp+yp*yp+zp*zp)

  if ( norm /= 0.0_rp) then


     xp=xp/norm
     yp=yp/norm
     zp=zp/norm

     x1p=px-x1
     y1p=py-y1
     z1p=pz-z1

     pscal=xp*x1p+yp*y1p+zp*z1p

     xp=pscal*xp
     yp=pscal*yp
     zp=pscal*zp

     x1p=x1p-xp
     y1p=y1p-yp
     z1p=z1p-zp


     xp=x1p+x1
     yp=y1p+y1
     zp=z1p+z1

     d1=abs(pscal)

     ! Check the projection point is inside the triangle

     sames=0_ip

     cp1x=(y2-y1)*(zp-z1)-(z2-z1)*(yp-y1)
     cp1y=(z2-z1)*(xp-x1)-(x2-x1)*(zp-z1)
     cp1z=(x2-x1)*(yp-y1)-(y2-y1)*(xp-x1)

     cp2x=(y2-y1)*(z3-z1)-(z2-z1)*(y3-y1)
     cp2y=(z2-z1)*(x3-x1)-(x2-x1)*(z3-z1)
     cp2z=(x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)

     pscal=cp1x*cp2x+cp1y*cp2y+cp1z*cp2z      

     if(pscal>=0.0_rp) then
        sames=sames+1_ip
     endif

     cp1x=(y3-y2)*(zp-z2)-(z3-z2)*(yp-y2)
     cp1y=(z3-z2)*(xp-x2)-(x3-x2)*(zp-z2)
     cp1z=(x3-x2)*(yp-y2)-(y3-y2)*(xp-x2)

     cp2x=(y3-y2)*(z1-z2)-(z3-z2)*(y1-y2)
     cp2y=(z3-z2)*(x1-x2)-(x3-x2)*(z1-z2)
     cp2z=(x3-x2)*(y1-y2)-(y3-y2)*(x1-x2)

     pscal=cp1x*cp2x+cp1y*cp2y+cp1z*cp2z      

     if(pscal>=0.0_rp) then
        sames=sames+1_ip
     endif

     cp1x=(y1-y3)*(zp-z3)-(z1-z3)*(yp-y3)
     cp1y=(z1-z3)*(xp-x3)-(x1-x3)*(zp-z3)
     cp1z=(x1-x3)*(yp-y3)-(y1-y3)*(xp-x3)

     cp2x=(y1-y3)*(z2-z3)-(z1-z3)*(y2-y3)
     cp2y=(z1-z3)*(x2-x3)-(x1-x3)*(z2-z3)
     cp2z=(x1-x3)*(y2-y3)-(y1-y3)*(x2-x3)

     pscal=cp1x*cp2x+cp1y*cp2y+cp1z*cp2z      

     if(pscal>=0.0_rp) then
        sames = sames+1_ip
     endif

     if(sames==3_ip) then
        d = d1
     else 

        ! First minize with distance to each point
        d = sqrt((px-x1)*(px-x1)+(py-y1)*(py-y1)+(pz-z1)*(pz-z1))

        d1 = sqrt((px-x2)*(px-x2)+(py-y2)*(py-y2)+(pz-z2)*(pz-z2))

        if(d1<d) then
           d = d1
        endif

        d1=sqrt((px-x3)*(px-x3)+(py-y3)*(py-y3)+(pz-z3)*(pz-z3))

        if(d1<d) then
           d = d1
        endif

     endif

  else  ! the case where the 3 points do not form a plane, only a line or even they are all three the same 

     d23 = sqrt( (x2-x3)*(x2-x3) + (y2-y3)*(y2-y3) +  (z2-z3)*(z2-z3) )
     d13 = sqrt( (x13)*(x13) + (y13)*(y13) +  (z13)*(z13) )
     d12 = sqrt( (x12)*(x12) + (y12)*(y12) +  (z12)*(z12) )

     if ( (d23 >= d13) .and. (d23 >= d12) ) then  ! d23 is the maximum distance between 1,2&3  or they are all the same
           call dipoli(px,py,pz,x2,y2,z2,x3,y3,z3,d)
     else if ( (d13 >= d23) .and. (d13 >= d12) )  then  ! d13 is the maximum distance between 1,2&3  or they are all the same
           call dipoli(px,py,pz,x1,y1,z1,x3,y3,z3,d)
     else if ( (d12 >= d23) .and. (d12 >= d13) )  then  ! d12 is the maximum distance between 1,2&3  or they are all the same
           call dipoli(px,py,pz,x1,y1,z1,x2,y2,z2,d)
     end if

  end if

end subroutine dipopl




subroutine dipoli(px,py,pz,x1,y1,z1,x2,y2,z2,d)

!-----------------------------------------------------------------------
!
! Computes the distance between a point and a line can handle the degenerate case or p1=p2
!
!-----------------------------------------------------------------------
  use      def_kintyp
  implicit none
  real(rp),    intent(in)  :: px,py,pz,x1,y1,z1,x2,y2,z2
  real(rp),    intent(out) :: d

  real(rp)                 :: pscal,t

  if (sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) +  (z2-z1)*(z2-z1) ) == 0.0_rp) then  ! p1=p2
     d =  sqrt( (px-x1)*(px-x1) + (py-y1)*(py-y1) +  (pz-z1)*(pz-z1) )
  else
     pscal = (x1-px)*(x2-x1) + (y1-py)*(y2-y1) +  (z1-pz)*(z2-z1)
     t = - pscal / ( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) +  (z2-z1)*(z2-z1) )

     if (t<=0.0_rp) then
        d =  sqrt( (px-x1)*(px-x1) + (py-y1)*(py-y1) +  (pz-z1)*(pz-z1) )
     else if (t>=1.0_rp) then
        d =  sqrt( (px-x2)*(px-x2) + (py-y2)*(py-y2) +  (pz-z2)*(pz-z2) )
     else
        d =  sqrt( ((x1-px)+((x2-x1)*t))*((x1-px)+((x2-x1)*t))    &
             &   + ((y1-py)+((y2-y1)*t))*((y1-py)+((y2-y1)*t))    & 
             &   + ((z1-pz)+((z2-z1)*t))*((z1-pz)+((z2-z1)*t))  )
     end if
  end if

end subroutine dipoli
