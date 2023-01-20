!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine elmdel_der(pnode,ndime,elcod,gpdet_der,gpcar_der)
  !-----------------------------------------------------------------------
  !****f* Domain/elmdel_der
  ! NAME
  !    elmdel_der
  ! DESCRIPTION
  !    This routine calculates the spacial gradients of the Cartesian derivatives for 2D and 3D P1 element.
  !    
  ! USES
  ! USED BY
  ! SOURCE
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  implicit none
  integer(ip), intent(in)  :: pnode,ndime
  real(rp),    intent(in)  :: elcod(ndime,pnode)
  real(rp),    intent(out) :: gpdet_der(ndime,pnode),gpcar_der(ndime,pnode,ndime,pnode)
  real(rp)                 :: denom,gpdet,denom_der(ndime,pnode)
!  real(rp)                 :: b(3,3),t1,t2,t3
  integer(ip)              :: idime,inode
  real(rp)                 :: x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,aux1

  if( ndime == 2 .and. pnode == 3 ) then
     !
     ! 2D P1 element
     !
     gpdet      =  (-elcod(1,1)+elcod(1,2))*(-elcod(2,1)+elcod(2,3))&
          &       -(-elcod(2,1)+elcod(2,2))*(-elcod(1,1)+elcod(1,3))

     gpdet_der(1,1)      = -(-elcod(2,1)+elcod(2,3)) + (-elcod(2,1)+elcod(2,2))
     gpdet_der(2,1)      = -(-elcod(1,1)+elcod(1,2)) + (-elcod(1,1)+elcod(1,3))
     gpdet_der(1,2)      =  (-elcod(2,1)+elcod(2,3))
     gpdet_der(2,2)      = -(-elcod(1,1)+elcod(1,3))
     gpdet_der(1,3)      = -(-elcod(2,1)+elcod(2,2))
     gpdet_der(2,3)      =  (-elcod(1,1)+elcod(1,2))
          
          
          
     denom      =  1.0_rp/gpdet
     do idime = 1,ndime
       do inode = 1, pnode
         denom_der(idime,inode)      =  -gpdet_der(idime,inode)/(gpdet*gpdet)
       enddo
     enddo
     
!      gpcar(1,1) =  ( -elcod(2,3) + elcod(2,2) ) * denom
!      gpcar(1,2) =  ( -elcod(2,1) + elcod(2,3) ) * denom
!      gpcar(1,3) =  (  elcod(2,1) - elcod(2,2) ) * denom
!      gpcar(2,1) =  (  elcod(1,3) - elcod(1,2) ) * denom
!      gpcar(2,2) =  (  elcod(1,1) - elcod(1,3) ) * denom
!      gpcar(2,3) =  ( -elcod(1,1) + elcod(1,2) ) * denom
     
     do idime = 1,ndime
       do inode = 1, pnode
         gpcar_der(1,1,idime,inode) =  ( -elcod(2,3) + elcod(2,2) ) * denom_der(idime,inode)
         gpcar_der(1,2,idime,inode) =  ( -elcod(2,1) + elcod(2,3) ) * denom_der(idime,inode)
         gpcar_der(1,3,idime,inode) =  (  elcod(2,1) - elcod(2,2) ) * denom_der(idime,inode)
         gpcar_der(2,1,idime,inode) =  (  elcod(1,3) - elcod(1,2) ) * denom_der(idime,inode)
         gpcar_der(2,2,idime,inode) =  (  elcod(1,1) - elcod(1,3) ) * denom_der(idime,inode)
         gpcar_der(2,3,idime,inode) =  ( -elcod(1,1) + elcod(1,2) ) * denom_der(idime,inode)
       enddo
     enddo
     
     gpcar_der(1,1,2,3) = gpcar_der(1,1,2,3) - denom
     gpcar_der(1,1,2,2) = gpcar_der(1,1,2,2) + denom
     
     gpcar_der(1,2,2,1) = gpcar_der(1,2,2,1) - denom
     gpcar_der(1,2,2,3) = gpcar_der(1,2,2,3) + denom
     
     gpcar_der(1,3,2,1) = gpcar_der(1,3,2,1) + denom
     gpcar_der(1,3,2,2) = gpcar_der(1,3,2,2) - denom
     
     gpcar_der(2,1,1,3) = gpcar_der(2,1,1,3) + denom
     gpcar_der(2,1,1,2) = gpcar_der(2,1,1,2) - denom
     
     gpcar_der(2,2,1,1) = gpcar_der(2,2,1,1) + denom
     gpcar_der(2,2,1,3) = gpcar_der(2,2,1,3) - denom
     
     gpcar_der(2,3,1,1) = gpcar_der(2,3,1,1) - denom
     gpcar_der(2,3,1,2) = gpcar_der(2,3,1,2) + denom
   
     
  else if( ndime == 3 .and. pnode == 4 ) then
     !
     ! 3D P1 element
     !
!      gpcar(1,1) = elcod(1,2) - elcod(1,1)
!      gpcar(1,2) = elcod(1,3) - elcod(1,1)
!      gpcar(1,3) = elcod(1,4) - elcod(1,1)
!      gpcar(2,1) = elcod(2,2) - elcod(2,1)
!      gpcar(2,2) = elcod(2,3) - elcod(2,1)
!      gpcar(2,3) = elcod(2,4) - elcod(2,1)
!      gpcar(3,1) = elcod(3,2) - elcod(3,1)
!      gpcar(3,2) = elcod(3,3) - elcod(3,1)
!      gpcar(3,3) = elcod(3,4) - elcod(3,1)     
          
!      t1         = gpcar(2,2)*gpcar(3,3) - gpcar(3,2)*gpcar(2,3)
!      t2         =-gpcar(2,1)*gpcar(3,3) + gpcar(3,1)*gpcar(2,3)
!      t3         = gpcar(2,1)*gpcar(3,2) - gpcar(3,1)*gpcar(2,2)
     
!      gpdet      = (gpcar(1,1)*t1 + gpcar(1,2)*t2 + gpcar(1,3)*t3)
     
     gpdet_der(1,1)      = elcod(2,2)*(elcod(3,4)-elcod(3,3)) + elcod(2,3)*(elcod(3,2)-elcod(3,4)) &
                         + elcod(2,4)*(elcod(3,3)-elcod(3,2))
     
     gpdet_der(1,2)      = elcod(2,1)*(elcod(3,3)-elcod(3,4)) + elcod(2,3)*(elcod(3,4)-elcod(3,1)) &
                         + elcod(2,4)*(elcod(3,1)-elcod(3,3))
     
     gpdet_der(1,3)      = elcod(2,1)*(elcod(3,4)-elcod(3,2)) + elcod(2,2)*(elcod(3,1)-elcod(3,4)) &
                         + elcod(2,4)*(elcod(3,2)-elcod(3,1))
                         
     gpdet_der(1,4)      = elcod(2,1)*(elcod(3,2)-elcod(3,3)) + elcod(2,2)*(elcod(3,3)-elcod(3,1)) &
                         + elcod(2,3)*(elcod(3,1)-elcod(3,2))

     gpdet_der(2,1)      = elcod(1,2)*(elcod(3,3)-elcod(3,4)) + elcod(1,3)*(elcod(3,4)-elcod(3,2)) &
                         + elcod(1,4)*(elcod(3,2)-elcod(3,3))

     gpdet_der(2,2)      = elcod(1,1)*(elcod(3,4)-elcod(3,3)) + elcod(1,3)*(elcod(3,1)-elcod(3,4)) &
                         + elcod(1,4)*(elcod(3,3)-elcod(3,1))
                         
     gpdet_der(2,3)      = elcod(1,1)*(elcod(3,2)-elcod(3,4)) + elcod(1,2)*(elcod(3,4)-elcod(3,1)) &
                         + elcod(1,4)*(elcod(3,1)-elcod(3,2))
                         
     gpdet_der(2,4)      = elcod(1,1)*(elcod(3,3)-elcod(3,2)) + elcod(1,2)*(elcod(3,1)-elcod(3,3)) &
                         + elcod(1,3)*(elcod(3,2)-elcod(3,1))
                         
     gpdet_der(3,1)      = elcod(1,2)*(elcod(2,4)-elcod(2,3)) + elcod(1,3)*(elcod(2,2)-elcod(2,4)) &
                         + elcod(1,4)*(elcod(2,3)-elcod(2,2))

     gpdet_der(3,2)      = elcod(1,1)*(elcod(2,3)-elcod(2,4)) + elcod(1,3)*(elcod(2,4)-elcod(2,1)) &
                         + elcod(1,4)*(elcod(2,1)-elcod(2,3))
                         
     gpdet_der(3,3)      = elcod(1,1)*(elcod(2,4)-elcod(2,2)) + elcod(1,2)*(elcod(2,1)-elcod(2,4)) &
                         + elcod(1,4)*(elcod(2,2)-elcod(2,1))
                         
     gpdet_der(3,4)      = elcod(1,1)*(elcod(2,2)-elcod(2,3)) + elcod(1,2)*(elcod(2,3)-elcod(2,1)) &
                         + elcod(1,3)*(elcod(2,1)-elcod(2,2))
                         
!      denom      = 1.0_rp/gpdet
! 
!      b(1,1)     = t1*denom
!      b(2,1)     = t2*denom
!      b(3,1)     = t3*denom
!      b(2,2)     = ( gpcar(1,1) * gpcar(3,3) - gpcar(3,1) * gpcar(1,3)) * denom
!      b(3,2)     = (-gpcar(1,1) * gpcar(3,2) + gpcar(1,2) * gpcar(3,1)) * denom
!      b(3,3)     = ( gpcar(1,1) * gpcar(2,2) - gpcar(2,1) * gpcar(1,2)) * denom
!      b(1,2)     = (-gpcar(1,2) * gpcar(3,3) + gpcar(3,2) * gpcar(1,3)) * denom
!      b(1,3)     = ( gpcar(1,2) * gpcar(2,3) - gpcar(2,2) * gpcar(1,3)) * denom
!      b(2,3)     = (-gpcar(1,1) * gpcar(2,3) + gpcar(2,1) * gpcar(1,3)) * denom
! 
!      gpcar(1,1) = -b(1,1)-b(2,1)-b(3,1)
!      gpcar(1,2) =  b(1,1)
!      gpcar(1,3) =  b(2,1)
!      gpcar(1,4) =  b(3,1)
!      gpcar(2,1) = -b(1,2)-b(2,2)-b(3,2)
!      gpcar(2,2) =  b(1,2)
!      gpcar(2,3) =  b(2,2)
!      gpcar(2,4) =  b(3,2)
!      gpcar(3,1) = -b(1,3)-b(2,3)-b(3,3)
!      gpcar(3,2) =  b(1,3)
!      gpcar(3,3) =  b(2,3)
!      gpcar(3,4) =  b(3,3)
       
       x1 = elcod(1,1)
       y1 = elcod(2,1)
       z1 = elcod(3,1)
       x2 = elcod(1,2)
       y2 = elcod(2,2)
       z2 = elcod(3,2)
       x3 = elcod(1,3)
       y3 = elcod(2,3)
       z3 = elcod(3,3)
       x4 = elcod(1,4)
       y4 = elcod(2,4)
       z4 = elcod(3,4)
       
       aux1 = (-x2*y3*z1+x2*y4 *z1+x1 *y3* z2-x1 *y4 *z2+x2* y1* z3-x1* y2* z3+x1 *y4* z3-x2* y4* z3+x4* &
              (-y2* z1+y3* z1+y1* z2-y3* z2-y1* z3+y2* z3)-x2* y1* z4+x1* y2* z4-x1* y3* z4+x2* y3* z4+x3* & 
              (-y4* z1-y1* z2+y4* z2+y2* (z1-z4)+y1* z4))
              
       !!! gpcar(1,1)
       gpcar_der(1,1,1,1) = -(y4* (-z2+z3)+y3* (z2-z4)+y2* (-z3+z4))**2.0_rp/aux1**2.0_rp
       gpcar_der(1,1,1,2) =  (y4* (-z1+z3)+y3* (z1-z4)+y1* (-z3+z4))* (y4* (-z2+z3)+y3* (z2-z4)+y2* (-z3+z4))/aux1**2.0_rp
       gpcar_der(1,1,1,3) =  (y4* (-z1+z2)+y2* (z1-z4)+y1* (-z2+z4))* (y4* (z2-z3)+y2* (z3-z4)+y3* (-z2+z4))/aux1**2.0_rp
       gpcar_der(1,1,1,4) = -(y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))* (y4* (-z2+z3)+y3* (z2-z4)+y2* (-z3+z4))/aux1**2.0_rp

       gpcar_der(1,1,2,1) =  (x4* (z2-z3)+x2* (z3-z4)+x3* (-z2+z4))* (y4* (z2-z3)+y2* (z3-z4)+y3* (-z2+z4))/aux1**2.0_rp
       gpcar_der(1,1,2,2) =  (x4* (z2-z3)+x2* (z3-z4)+x3* (-z2+z4))* (y4* (-z1+z3)+y3* (z1-z4)+y1* (-z3+z4))/aux1**2.0_rp
       gpcar_der(1,1,2,3) =  (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))* (x4* (z2-z3)+x2* (z3-z4)+x3* (-z2+z4))/aux1**2.0_rp
       gpcar_der(1,1,2,4) = -(y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))* (x4* (z2-z3)+x2* (z3-z4)+x3* (-z2+z4))/aux1**2.0_rp
       
       gpcar_der(1,1,3,1) =  (x4* (-y2+y3)+x3* (y2-y4)+x2* (-y3+y4))* (y4* (z2-z3)+y2* (z3-z4)+y3* (-z2+z4))/aux1**2.0_rp
       gpcar_der(1,1,3,2) =  (x4* (-y2+y3)+x3* (y2-y4)+x2* (-y3+y4))* (y4* (-z1+z3)+y3* (z1-z4)+y1* (-z3+z4))/aux1**2.0_rp
       gpcar_der(1,1,3,3) =  (x4* (y2-y3)+x2* (y3-y4)+x3* (-y2+y4))* (y4* (-z1+z2)+y2* (z1-z4)+y1* (-z2+z4))/aux1**2.0_rp
       gpcar_der(1,1,3,4) =  (x4* (-y2+y3)+x3* (y2-y4)+x2* (-y3+y4))* (y3* (-z1+z2)+y2* (z1-z3)+y1* (-z2+z3))/aux1**2.0_rp

       !!! gpcar(1,2)
       gpcar_der(1,2,1,1) = -(y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))* (y4* (-z2+z3)+y3* (z2-z4)+y2* (-z3+z4))/aux1**2.0_rp
       gpcar_der(1,2,1,2) = -(y4* (-z1+z3)+y3* (z1-z4)+y1* (-z3+z4))**2/aux1**2.0_rp
       gpcar_der(1,2,1,3) =  (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux1**2.0_rp
       gpcar_der(1,2,1,4) =  (y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))* (y4* (-z1+z3)+y3* (z1-z4)+y1* (-z3+z4))/aux1**2.0_rp
       
       gpcar_der(1,2,2,1) =  (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))* (y4* (-z2+z3)+y3* (z2-z4)+y2* (-z3+z4))/aux1**2.0_rp
       gpcar_der(1,2,2,2) =  (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux1**2.0_rp
       gpcar_der(1,2,2,3) = -(x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))* (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))/aux1**2.0_rp
       gpcar_der(1,2,2,4) =  (y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))* (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))/aux1**2.0_rp
       
       gpcar_der(1,2,3,1) =  (x4* (-y1+y3)+x3* (y1-y4)+x1* (-y3+y4))* (y4* (-z2+z3)+y3* (z2-z4)+y2* (-z3+z4))/aux1**2.0_rp
       gpcar_der(1,2,3,2) =  (x4* (-y1+y3)+x3* (y1-y4)+x1* (-y3+y4))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux1**2.0_rp
       gpcar_der(1,2,3,3) =  (x4* (y1-y3)+x1* (y3-y4)+x3* (-y1+y4))* (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))/aux1**2.0_rp
       gpcar_der(1,2,3,4) =  (x4* (-y1+y3)+x3* (y1-y4)+x1* (-y3+y4))* (y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))/aux1**2.0_rp

       !!! gpcar(1,3)
       gpcar_der(1,3,1,1) = -(y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))* (y4* (z2-z3)+y2* (z3-z4)+y3* (-z2+z4))/aux1**2.0_rp
       gpcar_der(1,3,1,2) =  (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux1**2.0_rp
       gpcar_der(1,3,1,3) = -(y4* (-z1+z2)+y2* (z1-z4)+y1* (-z2+z4))**2/aux1**2.0_rp
       gpcar_der(1,3,1,4) =  (y3* (-z1+z2)+y2* (z1-z3)+y1* (-z2+z3))* (y4* (-z1+z2)+y2* (z1-z4)+y1* (-z2+z4))/aux1**2.0_rp

       gpcar_der(1,3,2,1) =  (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))* (y4* (z2-z3)+y2* (z3-z4)+y3* (-z2+z4))/aux1**2.0_rp
       gpcar_der(1,3,2,2) = -(x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux1**2.0_rp
       gpcar_der(1,3,2,3) =  (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))* (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))/aux1**2.0_rp
       gpcar_der(1,3,2,4) = -(y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))* (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))/aux1**2.0_rp
       
       gpcar_der(1,3,3,1) =  (x4* (-y1+y2)+x2* (y1-y4)+x1* (-y2+y4))* (y4* (z2-z3)+y2* (z3-z4)+y3* (-z2+z4))/aux1**2.0_rp
       gpcar_der(1,3,3,2) =  (x4* (y1-y2)+x1* (y2-y4)+x2* (-y1+y4))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux1**2.0_rp
       gpcar_der(1,3,3,3) =  (x4* (-y1+y2)+x2* (y1-y4)+x1* (-y2+y4))* (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))/aux1**2.0_rp
       gpcar_der(1,3,3,4) =  (x4* (y1-y2)+x1* (y2-y4)+x2* (-y1+y4))* (y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))/aux1**2.0_rp
       
       !!! gpcar(1,4)
       gpcar_der(1,4,1,1) = -(y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))* (y4* (-z2+z3)+y3* (z2-z4)+y2* (-z3+z4))/aux1**2.0_rp
       gpcar_der(1,4,1,2) =  (y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))* (y4* (-z1+z3)+y3* (z1-z4)+y1* (-z3+z4))/aux1**2.0_rp
       gpcar_der(1,4,1,3) =  (y3* (-z1+z2)+y2* (z1-z3)+y1* (-z2+z3))* (y4* (-z1+z2)+y2* (z1-z4)+y1* (-z2+z4))/aux1**2.0_rp
       gpcar_der(1,4,1,4) = -(y3* (-z1+z2)+y2* (z1-z3)+y1* (-z2+z3))**2/aux1**2.0_rp

       gpcar_der(1,4,2,1) = -(x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))* (y4* (z2-z3)+y2* (z3-z4)+y3* (-z2+z4))/aux1**2.0_rp
       gpcar_der(1,4,2,2) =  (x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux1**2.0_rp
       gpcar_der(1,4,2,3) = -(x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))* (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))/aux1**2.0_rp
       gpcar_der(1,4,2,4) =  (x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))* (y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))/aux1**2.0_rp
       
       gpcar_der(1,4,3,1) =  (x3* (y1-y2)+x1* (y2-y3)+x2* (-y1+y3))* (y4* (z2-z3)+y2* (z3-z4)+y3* (-z2+z4))/aux1**2.0_rp
       gpcar_der(1,4,3,2) =  (x3* (-y1+y2)+x2* (y1-y3)+x1* (-y2+y3))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux1**2.0_rp
       gpcar_der(1,4,3,3) =  (x3* (y1-y2)+x1* (y2-y3)+x2* (-y1+y3))* (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))/aux1**2.0_rp
       gpcar_der(1,4,3,4) =  (x3* (-y1+y2)+x2* (y1-y3)+x1* (-y2+y3))* (y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))/aux1**2.0_rp
       
       !!! gpcar(2,1)
       gpcar_der(2,1,1,1) = (x4* (z2-z3)+x2* (z3-z4)+x3* (-z2+z4))* (y4* (z2-z3)+y2* (z3-z4)+y3* (-z2+z4))/aux1**2.0_rp
       gpcar_der(2,1,1,2) = (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))* (y4* (-z2+z3)+y3* (z2-z4)+y2* (-z3+z4))/aux1**2.0_rp
       gpcar_der(2,1,1,3) = (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))* (y4* (z2-z3)+y2* (z3-z4)+y3* (-z2+z4))/aux1**2.0_rp
       gpcar_der(2,1,1,4) = -(x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))* (y4* (z2-z3)+y2* (z3-z4)+y3* (-z2+z4))/aux1**2.0_rp

       gpcar_der(2,1,2,1) = -(x4* (-z2+z3)+x3* (z2-z4)+x2* (-z3+z4))**2/aux1**2.0_rp
       gpcar_der(2,1,2,2) = (x4* (-z1+z3)+x3* (z1-z4)+x1* (-z3+z4))* (x4* (-z2+z3)+x3* (z2-z4)+x2* (-z3+z4))/aux1**2.0_rp
       gpcar_der(2,1,2,3) = (x4* (-z1+z2)+x2* (z1-z4)+x1* (-z2+z4))* (x4* (z2-z3)+x2* (z3-z4)+x3* (-z2+z4))/aux1**2.0_rp
       gpcar_der(2,1,2,4) = -(x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))* (x4* (-z2+z3)+x3* (z2-z4)+x2* (-z3+z4))/aux1**2.0_rp

       gpcar_der(2,1,3,1) = (x4* (y2-y3)+x2* (y3-y4)+x3* (-y2+y4))* (x4* (z2-z3)+x2* (z3-z4)+x3* (-z2+z4))/aux1**2.0_rp
       gpcar_der(2,1,3,2) = -(x4* (y2-y3)+x2* (y3-y4)+x3* (-y2+y4))* (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))/aux1**2.0_rp
       gpcar_der(2,1,3,3) = (x4* (y2-y3)+x2* (y3-y4)+x3* (-y2+y4))* (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))/aux1**2.0_rp
       gpcar_der(2,1,3,4) = (x4* (-y2+y3)+x3* (y2-y4)+x2* (-y3+y4))* (x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))/aux1**2.0_rp
       
       !!! gpcar(2,2)
       gpcar_der(2,2,1,1) = (x4* (z2-z3)+x2* (z3-z4)+x3* (-z2+z4))* (y4* (-z1+z3)+y3* (z1-z4)+y1* (-z3+z4))/aux1**2.0_rp
       gpcar_der(2,2,1,2) = (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux1**2.0_rp
       gpcar_der(2,2,1,3) = -(x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux1**2.0_rp
       gpcar_der(2,2,1,4) = (x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux1**2.0_rp

       gpcar_der(2,2,2,1) = -(x4* (z2-z3)+x2* (z3-z4)+x3* (-z2+z4))* (x4* (-z1+z3)+x3* (z1-z4)+x1* (-z3+z4))/aux1**2.0_rp
       gpcar_der(2,2,2,2) = -(x4* (-z1+z3)+x3* (z1-z4)+x1* (-z3+z4))**2/aux1**2.0_rp
       gpcar_der(2,2,2,3) = (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))* (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))/aux1**2.0_rp
       gpcar_der(2,2,2,4) = (x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))* (x4* (-z1+z3)+x3* (z1-z4)+x1* (-z3+z4))/aux1**2.0_rp

       gpcar_der(2,2,3,1) = -(x4* (y1-y3)+x1* (y3-y4)+x3* (-y1+y4))* (x4* (z2-z3)+x2* (z3-z4)+x3* (-z2+z4))/aux1**2.0_rp
       gpcar_der(2,2,3,2) = (x4* (y1-y3)+x1* (y3-y4)+x3* (-y1+y4))* (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))/aux1**2.0_rp
       gpcar_der(2,2,3,3) = (x4* (-y1+y3)+x3* (y1-y4)+x1* (-y3+y4))* (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))/aux1**2.0_rp
       gpcar_der(2,2,3,4) = (x4* (y1-y3)+x1* (y3-y4)+x3* (-y1+y4))* (x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))/aux1**2.0_rp

       !!! gpcar(2,3)
       gpcar_der(2,3,1,1) = (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))* (x4* (z2-z3)+x2* (z3-z4)+x3* (-z2+z4))/aux1**2.0_rp
       gpcar_der(2,3,1,2) = -(x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))* (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))/aux1**2.0_rp
       gpcar_der(2,3,1,3) = (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))* (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))/aux1**2.0_rp
       gpcar_der(2,3,1,4) = -(x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))* (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))/aux1**2.0_rp

       gpcar_der(2,3,2,1) = -(x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))* (x4* (z2-z3)+x2* (z3-z4)+x3* (-z2+z4))/aux1**2.0_rp
       gpcar_der(2,3,2,2) = (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))* (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))/aux1**2.0_rp
       gpcar_der(2,3,2,3) = -(x4* (-z1+z2)+x2* (z1-z4)+x1* (-z2+z4))**2/aux1**2.0_rp
       gpcar_der(2,3,2,4) = (x3* (-z1+z2)+x2* (z1-z3)+x1* (-z2+z3))* (x4* (-z1+z2)+x2* (z1-z4)+x1* (-z2+z4))/aux1**2.0_rp

       gpcar_der(2,3,3,1) = (x4* (y1-y2)+x1* (y2-y4)+x2* (-y1+y4))* (x4* (z2-z3)+x2* (z3-z4)+x3* (-z2+z4))/aux1**2.0_rp
       gpcar_der(2,3,3,2) = (x4* (-y1+y2)+x2* (y1-y4)+x1* (-y2+y4))* (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))/aux1**2.0_rp
       gpcar_der(2,3,3,3) = (x4* (y1-y2)+x1* (y2-y4)+x2* (-y1+y4))* (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))/aux1**2.0_rp
       gpcar_der(2,3,3,4) = (x4* (y1-y2)+x1* (y2-y4)+x2* (-y1+y4))* (x3* (-z1+z2)+x2* (z1-z3)+x1* (-z2+z3))/aux1**2.0_rp

       !!! gpcar(2,4)
       gpcar_der(2,4,1,1) = -(y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))* (x4* (z2-z3)+x2* (z3-z4)+x3* (-z2+z4))/aux1**2.0_rp
       gpcar_der(2,4,1,2) = (y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))* (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))/aux1**2.0_rp
       gpcar_der(2,4,1,3) = -(y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))* (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))/aux1**2.0_rp
       gpcar_der(2,4,1,4) = (x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))* (y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))/aux1**2.0_rp

       gpcar_der(2,4,2,1) = -(x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))* (x4* (-z2+z3)+x3* (z2-z4)+x2* (-z3+z4))/aux1**2.0_rp
       gpcar_der(2,4,2,2) = (x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))* (x4* (-z1+z3)+x3* (z1-z4)+x1* (-z3+z4))/aux1**2.0_rp
       gpcar_der(2,4,2,3) = (x3* (-z1+z2)+x2* (z1-z3)+x1* (-z2+z3))* (x4* (-z1+z2)+x2* (z1-z4)+x1* (-z2+z4))/aux1**2.0_rp
       gpcar_der(2,4,2,4) = -(x3* (-z1+z2)+x2* (z1-z3)+x1* (-z2+z3))**2/aux1**2.0_rp

       gpcar_der(2,4,3,1) = (x3* (y1-y2)+x1* (y2-y3)+x2* (-y1+y3))* (x4* (-z2+z3)+x3* (z2-z4)+x2* (-z3+z4))/aux1**2.0_rp
       gpcar_der(2,4,3,2) = (x3* (-y1+y2)+x2* (y1-y3)+x1* (-y2+y3))* (x4* (-z1+z3)+x3* (z1-z4)+x1* (-z3+z4))/aux1**2.0_rp
       gpcar_der(2,4,3,3) = (x3* (y1-y2)+x1* (y2-y3)+x2* (-y1+y3))* (x4* (-z1+z2)+x2* (z1-z4)+x1* (-z2+z4))/aux1**2.0_rp
       gpcar_der(2,4,3,4) = (x3* (y1-y2)+x1* (y2-y3)+x2* (-y1+y3))* (x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))/aux1**2.0_rp
       
       !!! gpcar(3,1)
       gpcar_der(3,1,1,1) = (x4* (-y2+y3)+x3* (y2-y4)+x2* (-y3+y4))* (y4* (z2-z3)+y2* (z3-z4)+y3* (-z2+z4))/aux1**2.0_rp
       gpcar_der(3,1,1,2) = (x4* (-y1+y3)+x3* (y1-y4)+x1* (-y3+y4))* (y4* (-z2+z3)+y3* (z2-z4)+y2* (-z3+z4))/aux1**2.0_rp
       gpcar_der(3,1,1,3) = (x4* (-y1+y2)+x2* (y1-y4)+x1* (-y2+y4))* (y4* (z2-z3)+y2* (z3-z4)+y3* (-z2+z4))/aux1**2.0_rp
       gpcar_der(3,1,1,4) = (x3* (y1-y2)+x1* (y2-y3)+x2* (-y1+y3))* (y4* (z2-z3)+y2* (z3-z4)+y3* (-z2+z4))/aux1**2.0_rp

       gpcar_der(3,1,2,1) = (x4* (y2-y3)+x2* (y3-y4)+x3* (-y2+y4))* (x4* (z2-z3)+x2* (z3-z4)+x3* (-z2+z4))/aux1**2.0_rp
       gpcar_der(3,1,2,2) = -(x4* (y1-y3)+x1* (y3-y4)+x3* (-y1+y4))* (x4* (z2-z3)+x2* (z3-z4)+x3* (-z2+z4))/aux1**2.0_rp
       gpcar_der(3,1,2,3) = (x4* (y1-y2)+x1* (y2-y4)+x2* (-y1+y4))* (x4* (z2-z3)+x2* (z3-z4)+x3* (-z2+z4))/aux1**2.0_rp
       gpcar_der(3,1,2,4) = (x3* (y1-y2)+x1* (y2-y3)+x2* (-y1+y3))* (x4* (-z2+z3)+x3* (z2-z4)+x2* (-z3+z4))/aux1**2.0_rp

       gpcar_der(3,1,3,1) = -(x4* (-y2+y3)+x3* (y2-y4)+x2* (-y3+y4))**2/aux1**2.0_rp
       gpcar_der(3,1,3,2) = (x4* (-y1+y3)+x3* (y1-y4)+x1* (-y3+y4))* (x4* (-y2+y3)+x3* (y2-y4)+x2* (-y3+y4))/aux1**2.0_rp
       gpcar_der(3,1,3,3) = (x4* (-y1+y2)+x2* (y1-y4)+x1* (-y2+y4))* (x4* (y2-y3)+x2* (y3-y4)+x3* (-y2+y4))/aux1**2.0_rp
       gpcar_der(3,1,3,4) = -(x3* (y1-y2)+x1* (y2-y3)+x2* (-y1+y3))* (x4* (-y2+y3)+x3* (y2-y4)+x2* (-y3+y4))/aux1**2.0_rp
       
       !!! gpcar(3,2)
       gpcar_der(3,2,1,1) = (x4* (-y2+y3)+x3* (y2-y4)+x2* (-y3+y4))* (y4* (-z1+z3)+y3* (z1-z4)+y1* (-z3+z4))/aux1**2.0_rp
       gpcar_der(3,2,1,2) = (x4* (-y1+y3)+x3* (y1-y4)+x1* (-y3+y4))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux1**2.0_rp
       gpcar_der(3,2,1,3) = (x4* (y1-y2)+x1* (y2-y4)+x2* (-y1+y4))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux1**2.0_rp
       gpcar_der(3,2,1,4) = (x3* (-y1+y2)+x2* (y1-y3)+x1* (-y2+y3))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux1**2.0_rp

       gpcar_der(3,2,2,1) = -(x4* (y2-y3)+x2* (y3-y4)+x3* (-y2+y4))* (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))/aux1**2.0_rp
       gpcar_der(3,2,2,2) = (x4* (y1-y3)+x1* (y3-y4)+x3* (-y1+y4))* (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))/aux1**2.0_rp
       gpcar_der(3,2,2,3) = (x4* (-y1+y2)+x2* (y1-y4)+x1* (-y2+y4))* (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))/aux1**2.0_rp
       gpcar_der(3,2,2,4) = (x3* (-y1+y2)+x2* (y1-y3)+x1* (-y2+y3))* (x4* (-z1+z3)+x3* (z1-z4)+x1* (-z3+z4))/aux1**2.0_rp

       gpcar_der(3,2,3,1) = -(x4* (y1-y3)+x1* (y3-y4)+x3* (-y1+y4))* (x4* (-y2+y3)+x3* (y2-y4)+x2* (-y3+y4))/aux1**2.0_rp
       gpcar_der(3,2,3,2) = -(x4* (-y1+y3)+x3* (y1-y4)+x1* (-y3+y4))**2/aux1**2.0_rp
       gpcar_der(3,2,3,3) = (x4* (y1-y2)+x1* (y2-y4)+x2* (-y1+y4))* (x4* (y1-y3)+x1* (y3-y4)+x3* (-y1+y4))/aux1**2.0_rp
       gpcar_der(3,2,3,4) = (x3* (y1-y2)+x1* (y2-y3)+x2* (-y1+y3))* (x4* (-y1+y3)+x3* (y1-y4)+x1* (-y3+y4))/aux1**2.0_rp

       !!! gpcar(3,3)
       gpcar_der(3,3,1,1) = (x4* (y2-y3)+x2* (y3-y4)+x3* (-y2+y4))* (y4* (-z1+z2)+y2* (z1-z4)+y1* (-z2+z4))/aux1**2.0_rp
       gpcar_der(3,3,1,2) = (x4* (y1-y3)+x1* (y3-y4)+x3* (-y1+y4))* (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))/aux1**2.0_rp
       gpcar_der(3,3,1,3) = (x4* (-y1+y2)+x2* (y1-y4)+x1* (-y2+y4))* (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))/aux1**2.0_rp
       gpcar_der(3,3,1,4) = (x3* (y1-y2)+x1* (y2-y3)+x2* (-y1+y3))* (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))/aux1**2.0_rp

       gpcar_der(3,3,2,1) = (x4* (y2-y3)+x2* (y3-y4)+x3* (-y2+y4))* (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))/aux1**2.0_rp
       gpcar_der(3,3,2,2) = (x4* (-y1+y3)+x3* (y1-y4)+x1* (-y3+y4))* (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))/aux1**2.0_rp
       gpcar_der(3,3,2,3) = (x4* (y1-y2)+x1* (y2-y4)+x2* (-y1+y4))* (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))/aux1**2.0_rp
       gpcar_der(3,3,2,4) = (x3* (y1-y2)+x1* (y2-y3)+x2* (-y1+y3))* (x4* (-z1+z2)+x2* (z1-z4)+x1* (-z2+z4))/aux1**2.0_rp

       gpcar_der(3,3,3,1) = -(x4* (y1-y2)+x1* (y2-y4)+x2* (-y1+y4))* (x4* (y2-y3)+x2* (y3-y4)+x3* (-y2+y4))/aux1**2.0_rp
       gpcar_der(3,3,3,2) = (x4* (y1-y2)+x1* (y2-y4)+x2* (-y1+y4))* (x4* (y1-y3)+x1* (y3-y4)+x3* (-y1+y4))/aux1**2.0_rp
       gpcar_der(3,3,3,3) = -(x4* (-y1+y2)+x2* (y1-y4)+x1* (-y2+y4))**2/aux1**2.0_rp
       gpcar_der(3,3,3,4) = (x3* (-y1+y2)+x2* (y1-y3)+x1* (-y2+y3))* (x4* (-y1+y2)+x2* (y1-y4)+x1* (-y2+y4))/aux1**2.0_rp
       
       !!! gpcar(3,4)
       gpcar_der(3,4,1,1) = (x4* (-y2+y3)+x3* (y2-y4)+x2* (-y3+y4))* (y3* (-z1+z2)+y2* (z1-z3)+y1* (-z2+z3))/aux1**2.0_rp
       gpcar_der(3,4,1,2) = (x4* (-y1+y3)+x3* (y1-y4)+x1* (-y3+y4))* (y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))/aux1**2.0_rp
       gpcar_der(3,4,1,3) = (x4* (y1-y2)+x1* (y2-y4)+x2* (-y1+y4))* (y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))/aux1**2.0_rp
       gpcar_der(3,4,1,4) = (x3* (-y1+y2)+x2* (y1-y3)+x1* (-y2+y3))* (y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))/aux1**2.0_rp

       gpcar_der(3,4,2,1) = (x4* (-y2+y3)+x3* (y2-y4)+x2* (-y3+y4))* (x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))/aux1**2.0_rp
       gpcar_der(3,4,2,2) = (x4* (y1-y3)+x1* (y3-y4)+x3* (-y1+y4))* (x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))/aux1**2.0_rp
       gpcar_der(3,4,2,3) = (x4* (y1-y2)+x1* (y2-y4)+x2* (-y1+y4))* (x3* (-z1+z2)+x2* (z1-z3)+x1* (-z2+z3))/aux1**2.0_rp
       gpcar_der(3,4,2,4) = (x3* (y1-y2)+x1* (y2-y3)+x2* (-y1+y3))* (x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))/aux1**2.0_rp

       gpcar_der(3,4,3,1) = -(x3* (y1-y2)+x1* (y2-y3)+x2* (-y1+y3))* (x4* (-y2+y3)+x3* (y2-y4)+x2* (-y3+y4))/aux1**2.0_rp
       gpcar_der(3,4,3,2) = (x3* (y1-y2)+x1* (y2-y3)+x2* (-y1+y3))* (x4* (-y1+y3)+x3* (y1-y4)+x1* (-y3+y4))/aux1**2.0_rp
       gpcar_der(3,4,3,3) = (x3* (-y1+y2)+x2* (y1-y3)+x1* (-y2+y3))* (x4* (-y1+y2)+x2* (y1-y4)+x1* (-y2+y4))/aux1**2.0_rp
       gpcar_der(3,4,3,4) = -(x3* (-y1+y2)+x2* (y1-y3)+x1* (-y2+y3))**2/aux1**2.0_rp
       
  end if

end subroutine elmdel_der 
!-----------------------------------------------------------------------
! NOTES
!
! P1 Element in 2D
! ----------------
!
! gpdet=(-x1+x2)*(-y1+y3)-(-y1+y2)*(-x1+x3)
!                _                       _    
!               | -y3+y2  -y1+y3   y1-y2  |   
! gpcar=1/gpdet |                         |   
!               |_ x3-x2   x1-x3  -x1+x2 _|
!
! P1 Element in 3D
! ----------------
!        _                     _                        
!       |  x2-x1  x3-x1  x4-x1  |
! xjacm=|  y2-y1  y3-y1  y4-y1  |
!       |_ z2-z1  z3-z1  z4-z1 _|
!
!            -1          
! xjaci=xjacm  
!        _                             _     _          _
!       |  dN1/ds dN2/ds dN3/ds dN4/ds  |   |  -1 1 0 0  |
! deriv=|  dN1/dt dN2/dt dN3/dt_dN4/dt  | = |  -1 0 1 0  |
!       |_ dN1/dz dN2/dz dN3/dz_dN4/dz _|   |_ -1 0 0 1 _|
!
!            t           
! cartd=xjaci *deriv
!
!***
!-----------------------------------------------------------------------
