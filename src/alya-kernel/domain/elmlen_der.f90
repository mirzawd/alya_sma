!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine elmlen_der(&
     ndime,pnode,deriv,elcod,hnatu,hleng_der)
  !-----------------------------------------------------------------------
  !****f* Domain/elmlen_der
  ! NAME
  !    elmlen_der
  ! DESCRIPTION
  !    Compute element lengths
  ! OUTPUT
  !    HLENG ... Element length with: 
  !              HLENG(1)     = Max length
  !              HLENG(NDIME) = Min length
  ! USED BY
  !    
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     : ip,rp
  use def_kermod, only       :  kfl_adj_prob
  implicit none 
  integer(ip), intent(in)  :: ndime,pnode
  real(rp),    intent(in)  :: hnatu,elcod(ndime,pnode)
  real(rp),    intent(in)  :: deriv(ndime,pnode)
  real(rp),    intent(out) :: hleng_der(ndime,pnode,ndime)
  integer(ip)              :: k,idime,inode
  real(rp)                 :: enor0,gpdet,denom,hleng(ndime),tragl(ndime,ndime)
!  real(rp)                 :: h_tem
  real(rp)                 :: xjacm(ndime,ndime),t1,t2,t3
  real(rp)                 :: xjacm_der(ndime,ndime,pnode),gpdet_der(ndime,pnode)
  real(rp)                 :: tragl_der(ndime,ndime,ndime,pnode),denom_der(ndime,pnode)
  real(rp)                 :: enor0_der(ndime,pnode),h_tem_der(ndime,pnode)
  real(rp)                 :: x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,aux
  
  
  if( kfl_adj_prob == 0 ) return
  
  tragl_der = 0.0_rp
  
  if( ndime == 1 ) then

     call runend('No esta implementado') 
!      xjacm(1,1) = 0.0_rp
!      do k = 1,pnode
!         xjacm(1,1) = xjacm(1,1) + elcod(1,k) * deriv(1,k)
!      end do
!      tragl(1,1) = 1.0_rp/xjacm(1,1)
!      enor0      = tragl(1,1) * tragl(1,1)
!      hleng(1)   = hnatu/sqrt(enor0)

  else if( ndime == 2 ) then

     xjacm(1,1) = 0.0_rp
     xjacm(1,2) = 0.0_rp
     xjacm(2,1) = 0.0_rp
     xjacm(2,2) = 0.0_rp
     do k = 1,pnode
        xjacm(1,1) = xjacm(1,1) + elcod(1,k) * deriv(1,k)
        xjacm(1,2) = xjacm(1,2) + elcod(1,k) * deriv(2,k)
        xjacm(2,1) = xjacm(2,1) + elcod(2,k) * deriv(1,k)
        xjacm(2,2) = xjacm(2,2) + elcod(2,k) * deriv(2,k)
     end do

     gpdet      =  xjacm(1,1) * xjacm(2,2) - xjacm(2,1) * xjacm(1,2)
     
     do k = 1,pnode
        xjacm_der(1,1,k) = deriv(1,k)
        xjacm_der(1,2,k) = deriv(2,k)
        xjacm_der(2,1,k) = deriv(1,k)
        xjacm_der(2,2,k) = deriv(2,k)
     end do

     do k = 1,pnode
        gpdet_der(1,k) = xjacm_der(1,1,k) * xjacm(2,2) - xjacm(2,1) * xjacm_der(1,2,k)
        gpdet_der(2,k) = xjacm(1,1) * xjacm_der(2,2,k) - xjacm_der(2,1,k) * xjacm(1,2)
     enddo
     
     denom      =  1.0_rp/gpdet
     do idime = 1,ndime
       do inode = 1, pnode
         denom_der(idime,inode)      =  -gpdet_der(idime,inode)/(gpdet*gpdet)
       enddo
     enddo
     
     tragl(1,1) =  xjacm(2,2) * denom
     tragl(2,2) =  xjacm(1,1) * denom
     tragl(2,1) = -xjacm(2,1) * denom
     tragl(1,2) = -xjacm(1,2) * denom  
     
     do idime = 1,ndime
       do inode = 1, pnode
         tragl_der(1,1,idime,inode) =  xjacm(2,2) * denom_der(idime,inode)
         tragl_der(2,2,idime,inode) =  xjacm(1,1) * denom_der(idime,inode)
         tragl_der(2,1,idime,inode) = -xjacm(2,1) * denom_der(idime,inode)
         tragl_der(1,2,idime,inode) = -xjacm(1,2) * denom_der(idime,inode)
       enddo
     enddo
     
     tragl_der(1,1,2,1) = tragl_der(1,1,2,1) + xjacm_der(2,2,1) * denom
     tragl_der(1,1,2,2) = tragl_der(1,1,2,2) + xjacm_der(2,2,2) * denom
     tragl_der(1,1,2,3) = tragl_der(1,1,2,3) + xjacm_der(2,2,3) * denom

     tragl_der(2,2,1,1) = tragl_der(2,2,1,1) + xjacm_der(1,1,1) * denom
     tragl_der(2,2,1,2) = tragl_der(2,2,1,2) + xjacm_der(1,1,2) * denom
     tragl_der(2,2,1,3) = tragl_der(2,2,1,3) + xjacm_der(1,1,3) * denom
     
     tragl_der(2,1,2,1) = tragl_der(2,1,2,1) - xjacm_der(2,1,1) * denom
     tragl_der(2,1,2,2) = tragl_der(2,1,2,2) - xjacm_der(2,1,2) * denom
     tragl_der(2,1,2,3) = tragl_der(2,1,2,3) - xjacm_der(2,1,3) * denom

     tragl_der(1,2,1,1) = tragl_der(1,2,1,1) - xjacm_der(1,2,1) * denom  
     tragl_der(1,2,1,2) = tragl_der(1,2,1,2) - xjacm_der(1,2,2) * denom  
     tragl_der(1,2,1,3) = tragl_der(1,2,1,3) - xjacm_der(1,2,3) * denom  
     
     
     enor0      =  tragl(1,1) * tragl(1,1) + tragl(1,2) * tragl(1,2) 
     hleng(1)   =  hnatu/sqrt(enor0)
     enor0      =  tragl(2,1) * tragl(2,1) + tragl(2,2) * tragl(2,2) 
     hleng(2)   =  hnatu/sqrt(enor0)
     
!      if( hleng(2) > hleng(1) ) then
!         h_tem    = hleng(2)
!         hleng(2) = hleng(1)
!         hleng(1) = h_tem
!      end if
     
     do idime = 1,ndime
       do inode = 1, pnode
         
         enor0_der(idime,inode)      =  2*tragl_der(1,1,idime,inode)*tragl(1,1) + 2*tragl_der(1,2,idime,inode)*tragl(1,2)
         enor0      =  tragl(1,1) * tragl(1,1) + tragl(1,2) * tragl(1,2)
         hleng_der(idime,inode,1)   =  -0.5_rp*hnatu*enor0_der(idime,inode)/(sqrt(enor0)*enor0)
         enor0_der(idime,inode)      =  2*tragl_der(2,1,idime,inode)*tragl(2,1) + 2*tragl_der(2,2,idime,inode)*tragl(2,2) 
         enor0      =  tragl(2,1) * tragl(2,1) + tragl(2,2) * tragl(2,2) 
         hleng_der(idime,inode,2)   =  -0.5_rp*hnatu*enor0_der(idime,inode)/(sqrt(enor0)*enor0)
         
         if( hleng(2) > hleng(1) ) then
           h_tem_der(idime,inode)   = hleng_der(idime,inode,2)
           hleng_der(idime,inode,2) = hleng_der(idime,inode,1)
           hleng_der(idime,inode,1) = h_tem_der(idime,inode)
         end if
         
       enddo
     enddo
     

  else if( ndime == 3 ) then

     
     xjacm(1,1) = 0.0_rp ! xjacm = elcod * deriv^t
     xjacm(1,2) = 0.0_rp ! tragl = xjacm^-1
     xjacm(1,3) = 0.0_rp 
     xjacm(2,1) = 0.0_rp
     xjacm(2,2) = 0.0_rp
     xjacm(2,3) = 0.0_rp
     xjacm(3,1) = 0.0_rp
     xjacm(3,2) = 0.0_rp
     xjacm(3,3) = 0.0_rp
     do k = 1,pnode
        xjacm(1,1) = xjacm(1,1) + elcod(1,k) * deriv(1,k)
        xjacm(1,2) = xjacm(1,2) + elcod(1,k) * deriv(2,k)
        xjacm(1,3) = xjacm(1,3) + elcod(1,k) * deriv(3,k)
        xjacm(2,1) = xjacm(2,1) + elcod(2,k) * deriv(1,k)
        xjacm(2,2) = xjacm(2,2) + elcod(2,k) * deriv(2,k)
        xjacm(2,3) = xjacm(2,3) + elcod(2,k) * deriv(3,k)
        xjacm(3,1) = xjacm(3,1) + elcod(3,k) * deriv(1,k)
        xjacm(3,2) = xjacm(3,2) + elcod(3,k) * deriv(2,k)
        xjacm(3,3) = xjacm(3,3) + elcod(3,k) * deriv(3,k)
     end do

     t1         =  xjacm(2,2) * xjacm(3,3) - xjacm(3,2) * xjacm(2,3)
     t2         = -xjacm(2,1) * xjacm(3,3) + xjacm(3,1) * xjacm(2,3)
     t3         =  xjacm(2,1) * xjacm(3,2) - xjacm(3,1) * xjacm(2,2)
     gpdet      =  xjacm(1,1) * t1 + xjacm(1,2) * t2 + xjacm(1,3) * t3
     denom      =  1.0_rp / gpdet
     tragl(1,1) =  t1 * denom
     tragl(2,1) =  t2 * denom
     tragl(3,1) =  t3 * denom
     tragl(2,2) =  ( xjacm(1,1) * xjacm(3,3) - xjacm(3,1) * xjacm(1,3)) * denom
     tragl(3,2) =  (-xjacm(1,1) * xjacm(3,2) + xjacm(1,2) * xjacm(3,1)) * denom
     tragl(3,3) =  ( xjacm(1,1) * xjacm(2,2) - xjacm(2,1) * xjacm(1,2)) * denom
     tragl(1,2) =  (-xjacm(1,2) * xjacm(3,3) + xjacm(3,2) * xjacm(1,3)) * denom
     tragl(1,3) =  ( xjacm(1,2) * xjacm(2,3) - xjacm(2,2) * xjacm(1,3)) * denom
     tragl(2,3) =  (-xjacm(1,1) * xjacm(2,3) + xjacm(2,1) * xjacm(1,3)) * denom
     
     if (pnode == 4) then
     
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
     
       aux = (-x2* y3* z1+x2* y4* z1+x1* y3* z2-x1* y4* z2+x2* y1* z3-x1* y2* z3+x1* y4* z3-&
              x2* y4* z3+x4* (-y2* z1+y3* z1+y1* z2-y3* z2-y1* z3+y2* z3)-x2* y1* z4+x1* y2* z4-x1* y3* z4+&
              x2* y3* z4+x3* (-y4* z1-y1* z2+y4* z2+y2* (z1-z4)+y1* z4))
       
       ! dtragl(1,1)
       tragl_der(1,1,1,1) = -(y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))* (y4* (-z2+z3)+y3* (z2-z4)+y2* (-z3+z4))/aux**2
       tragl_der(1,1,1,2) = -(y4* (-z1+z3)+y3* (z1-z4)+y1* (-z3+z4))**2/aux**2
       tragl_der(1,1,1,3) = (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux**2
       tragl_der(1,1,1,4) = (y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))* (y4* (-z1+z3)+y3* (z1-z4)+y1* (-z3+z4))/aux**2
     
       tragl_der(1,1,2,1) = (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))* (y4* (-z2+z3)+y3* (z2-z4)+y2* (-z3+z4))/aux**2
       tragl_der(1,1,2,2) = (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux**2
       tragl_der(1,1,2,3) = -(x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))* (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))/aux**2
       tragl_der(1,1,2,4) = (y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))* (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))/aux**2
     
       tragl_der(1,1,3,1) = (x4* (-y1+y3)+x3* (y1-y4)+x1* (-y3+y4))* (y4* (-z2+z3)+y3* (z2-z4)+y2* (-z3+z4))/aux**2
       tragl_der(1,1,3,2) = (x4* (-y1+y3)+x3* (y1-y4)+x1* (-y3+y4))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux**2
       tragl_der(1,1,3,3) = (x4* (y1-y3)+x1* (y3-y4)+x3* (-y1+y4))* (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))/aux**2
       tragl_der(1,1,3,4) = (x4* (-y1+y3)+x3* (y1-y4)+x1* (-y3+y4))* (y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))/aux**2
     
       ! dtragl(1,2)
       tragl_der(1,2,1,1) = (x4* (z2-z3)+x2* (z3-z4)+x3* (-z2+z4))* (y4* (-z1+z3)+y3* (z1-z4)+y1* (-z3+z4))/aux**2
       tragl_der(1,2,1,2) = (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux**2
       tragl_der(1,2,1,3) = -(x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux**2
       tragl_der(1,2,1,4) = (x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux**2
       
       tragl_der(1,2,2,1) = -(x4* (z2-z3)+x2* (z3-z4)+x3* (-z2+z4))* (x4* (-z1+z3)+x3* (z1-z4)+x1* (-z3+z4))/aux**2
       tragl_der(1,2,2,2) = -(x4* (-z1+z3)+x3* (z1-z4)+x1* (-z3+z4))**2/aux**2
       tragl_der(1,2,2,3) = (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))* (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))/aux**2
       tragl_der(1,2,2,4) = (x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))* (x4* (-z1+z3)+x3* (z1-z4)+x1* (-z3+z4))/aux**2
       
       tragl_der(1,2,3,1) = -(x4* (y1-y3)+x1* (y3-y4)+x3* (-y1+y4))* (x4* (z2-z3)+x2* (z3-z4)+x3* (-z2+z4))/aux**2
       tragl_der(1,2,3,2) = (x4* (y1-y3)+x1* (y3-y4)+x3* (-y1+y4))* (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))/aux**2
       tragl_der(1,2,3,3) = (x4* (-y1+y3)+x3* (y1-y4)+x1* (-y3+y4))* (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))/aux**2
       tragl_der(1,2,3,4) = (x4* (y1-y3)+x1* (y3-y4)+x3* (-y1+y4))* (x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))/aux**2
       
       ! dtragl(1,3)       
       tragl_der(1,3,1,1) = (x4* (-y2+y3)+x3* (y2-y4)+x2* (-y3+y4))* (y4* (-z1+z3)+y3* (z1-z4)+y1* (-z3+z4))/aux**2
       tragl_der(1,3,1,2) = (x4* (-y1+y3)+x3* (y1-y4)+x1* (-y3+y4))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux**2
       tragl_der(1,3,1,3) = (x4* (y1-y2)+x1* (y2-y4)+x2* (-y1+y4))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux**2
       tragl_der(1,3,1,4) = (x3* (-y1+y2)+x2* (y1-y3)+x1* (-y2+y3))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux**2
       
       tragl_der(1,3,2,1) = -(x4* (y2-y3)+x2* (y3-y4)+x3* (-y2+y4))* (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))/aux**2
       tragl_der(1,3,2,2) = (x4* (y1-y3)+x1* (y3-y4)+x3* (-y1+y4))* (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))/aux**2
       tragl_der(1,3,2,3) = (x4* (-y1+y2)+x2* (y1-y4)+x1* (-y2+y4))* (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))/aux**2
       tragl_der(1,3,2,4) = (x3* (-y1+y2)+x2* (y1-y3)+x1* (-y2+y3))* (x4* (-z1+z3)+x3* (z1-z4)+x1* (-z3+z4))/aux**2

       tragl_der(1,3,3,1) = -(x4* (y1-y3)+x1* (y3-y4)+x3* (-y1+y4))* (x4* (-y2+y3)+x3* (y2-y4)+x2* (-y3+y4))/aux**2
       tragl_der(1,3,3,2) = -(x4* (-y1+y3)+x3* (y1-y4)+x1* (-y3+y4))**2/aux**2
       tragl_der(1,3,3,3) = (x4* (y1-y2)+x1* (y2-y4)+x2* (-y1+y4))* (x4* (y1-y3)+x1* (y3-y4)+x3* (-y1+y4))/aux**2
       tragl_der(1,3,3,4) = (x3* (y1-y2)+x1* (y2-y3)+x2* (-y1+y3))* (x4* (-y1+y3)+x3* (y1-y4)+x1* (-y3+y4))/aux**2
       
       ! dtragl(2,1)       
       tragl_der(2,1,1,1) = -(y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))* (y4* (z2-z3)+y2* (z3-z4)+y3* (-z2+z4))/aux**2
       tragl_der(2,1,1,2) = (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux**2
       tragl_der(2,1,1,3) = -(y4* (-z1+z2)+y2* (z1-z4)+y1* (-z2+z4))/aux**2
       tragl_der(2,1,1,4) = (y3* (-z1+z2)+y2* (z1-z3)+y1* (-z2+z3))* (y4* (-z1+z2)+y2* (z1-z4)+y1* (-z2+z4))/aux**2
       
       tragl_der(2,1,2,1) = (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))* (y4* (z2-z3)+y2* (z3-z4)+y3* (-z2+z4))/aux**2
       tragl_der(2,1,2,2) = -(x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux**2
       tragl_der(2,1,2,3) = (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))* (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))/aux**2
       tragl_der(2,1,2,4) = -(y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))* (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))/aux**2
       
       tragl_der(2,1,3,1) = (x4* (-y1+y2)+x2* (y1-y4)+x1* (-y2+y4))* (y4* (z2-z3)+y2* (z3-z4)+y3* (-z2+z4))/aux**2
       tragl_der(2,1,3,2) = (x4* (y1-y2)+x1* (y2-y4)+x2* (-y1+y4))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux**2
       tragl_der(2,1,3,3) = (x4* (-y1+y2)+x2* (y1-y4)+x1* (-y2+y4))* (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))/aux**2
       tragl_der(2,1,3,4) = (x4* (y1-y2)+x1* (y2-y4)+x2* (-y1+y4))* (y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))/aux**2
       
       ! dtragl(2,2)       
       tragl_der(2,2,1,1) = (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))* (x4* (z2-z3)+x2* (z3-z4)+x3* (-z2+z4))/aux**2
       tragl_der(2,2,1,2) = -(x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))* (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))/aux**2
       tragl_der(2,2,1,3) = (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))* (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))/aux**2
       tragl_der(2,2,1,4) = -(x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))* (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))/aux**2
       
       tragl_der(2,2,2,1) = -(x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))* (x4* (z2-z3)+x2* (z3-z4)+x3* (-z2+z4))/aux**2
       tragl_der(2,2,2,2) = (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))* (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))/aux**2
       tragl_der(2,2,2,3) = -(x4* (-z1+z2)+x2* (z1-z4)+x1* (-z2+z4))**2/aux**2
       tragl_der(2,2,2,4) = (x3* (-z1+z2)+x2* (z1-z3)+x1* (-z2+z3))* (x4* (-z1+z2)+x2* (z1-z4)+x1* (-z2+z4))/aux**2

       tragl_der(2,2,3,1) = (x4* (y1-y2)+x1* (y2-y4)+x2* (-y1+y4))* (x4* (z2-z3)+x2* (z3-z4)+x3* (-z2+z4))/aux**2
       tragl_der(2,2,3,2) = (x4* (-y1+y2)+x2* (y1-y4)+x1* (-y2+y4))* (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))/aux**2
       tragl_der(2,2,3,3) = (x4* (y1-y2)+x1* (y2-y4)+x2* (-y1+y4))* (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))/aux**2
       tragl_der(2,2,3,4) = (x4* (y1-y2)+x1* (y2-y4)+x2* (-y1+y4))* (x3* (-z1+z2)+x2* (z1-z3)+x1* (-z2+z3))/aux**2
       
       ! dtragl(2,3)       
       tragl_der(2,3,1,1) = (x4* (y2-y3)+x2* (y3-y4)+x3* (-y2+y4))* (y4* (-z1+z2)+y2* (z1-z4)+y1* (-z2+z4))/aux**2
       tragl_der(2,3,1,2) = (x4* (y1-y3)+x1* (y3-y4)+x3* (-y1+y4))* (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))/aux**2
       tragl_der(2,3,1,3) = (x4* (-y1+y2)+x2* (y1-y4)+x1* (-y2+y4))* (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))/aux**2
       tragl_der(2,3,1,4) = (x3* (y1-y2)+x1* (y2-y3)+x2* (-y1+y3))* (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))/aux**2
       
       tragl_der(2,3,2,1) = (x4* (y2-y3)+x2* (y3-y4)+x3* (-y2+y4))* (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))/aux**2
       tragl_der(2,3,2,2) = (x4* (-y1+y3)+x3* (y1-y4)+x1* (-y3+y4))* (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))/aux**2
       tragl_der(2,3,2,3) = (x4* (y1-y2)+x1* (y2-y4)+x2* (-y1+y4))* (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))/aux**2
       tragl_der(2,3,2,4) = (x3* (y1-y2)+x1* (y2-y3)+x2* (-y1+y3))* (x4* (-z1+z2)+x2* (z1-z4)+x1* (-z2+z4))/aux**2

       tragl_der(2,3,3,1) = -(x4* (y1-y2)+x1* (y2-y4)+x2* (-y1+y4))* (x4* (y2-y3)+x2* (y3-y4)+x3* (-y2+y4))/aux**2
       tragl_der(2,3,3,2) = (x4* (y1-y2)+x1* (y2-y4)+x2* (-y1+y4))* (x4* (y1-y3)+x1* (y3-y4)+x3* (-y1+y4))/aux**2
       tragl_der(2,3,3,3) = -(x4* (-y1+y2)+x2* (y1-y4)+x1* (-y2+y4))/aux**2
       tragl_der(2,3,3,4) = (x3* (-y1+y2)+x2* (y1-y3)+x1* (-y2+y3))* (x4* (-y1+y2)+x2* (y1-y4)+x1* (-y2+y4))/aux**2
       
       ! dtragl(3,1)       
       tragl_der(3,1,1,1) = -(y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))* (y4* (-z2+z3)+y3* (z2-z4)+y2* (-z3+z4))/aux**2
       tragl_der(3,1,1,2) = (y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))* (y4* (-z1+z3)+y3* (z1-z4)+y1* (-z3+z4))/aux**2
       tragl_der(3,1,1,3) = (y3* (-z1+z2)+y2* (z1-z3)+y1* (-z2+z3))* (y4* (-z1+z2)+y2* (z1-z4)+y1* (-z2+z4))/aux**2
       tragl_der(3,1,1,4) = -(y3* (-z1+z2)+y2* (z1-z3)+y1* (-z2+z3))**2/aux**2
       
       tragl_der(3,1,2,1) = -(x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))* (y4* (z2-z3)+y2* (z3-z4)+y3* (-z2+z4))/aux**2
       tragl_der(3,1,2,2) = (x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux**2
       tragl_der(3,1,2,3) = -(x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))* (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))/aux**2
       tragl_der(3,1,2,4) = (x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))* (y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))/aux**2

       tragl_der(3,1,3,1) = (x3* (y1-y2)+x1* (y2-y3)+x2* (-y1+y3))* (y4* (z2-z3)+y2* (z3-z4)+y3* (-z2+z4))/aux**2
       tragl_der(3,1,3,2) = (x3* (-y1+y2)+x2* (y1-y3)+x1* (-y2+y3))* (y4* (z1-z3)+y1* (z3-z4)+y3* (-z1+z4))/aux**2
       tragl_der(3,1,3,3) = (x3* (y1-y2)+x1* (y2-y3)+x2* (-y1+y3))* (y4* (z1-z2)+y1* (z2-z4)+y2* (-z1+z4))/aux**2
       tragl_der(3,1,3,4) = (x3* (-y1+y2)+x2* (y1-y3)+x1* (-y2+y3))* (y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))/aux**2
       
       ! dtragl(3,2)       
       tragl_der(3,2,1,1) = -(y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))* (x4* (z2-z3)+x2* (z3-z4)+x3* (-z2+z4))/aux**2
       tragl_der(3,2,1,2) = (y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))* (x4* (z1-z3)+x1* (z3-z4)+x3* (-z1+z4))/aux**2
       tragl_der(3,2,1,3) = -(y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))* (x4* (z1-z2)+x1* (z2-z4)+x2* (-z1+z4))/aux**2
       tragl_der(3,2,1,4) = (x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))* (y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))/aux**2
       
       tragl_der(3,2,2,1) = -(x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))* (x4* (-z2+z3)+x3* (z2-z4)+x2* (-z3+z4))/aux**2
       tragl_der(3,2,2,2) = (x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))* (x4* (-z1+z3)+x3* (z1-z4)+x1* (-z3+z4))/aux**2
       tragl_der(3,2,2,3) = (x3* (-z1+z2)+x2* (z1-z3)+x1* (-z2+z3))* (x4* (-z1+z2)+x2* (z1-z4)+x1* (-z2+z4))/aux**2
       tragl_der(3,2,2,4) = -(x3* (-z1+z2)+x2* (z1-z3)+x1* (-z2+z3))**2/aux**2

       tragl_der(3,2,3,1) = (x3* (y1-y2)+x1* (y2-y3)+x2* (-y1+y3))* (x4* (-z2+z3)+x3* (z2-z4)+x2* (-z3+z4))/aux**2
       tragl_der(3,2,3,2) = (x3* (-y1+y2)+x2* (y1-y3)+x1* (-y2+y3))* (x4* (-z1+z3)+x3* (z1-z4)+x1* (-z3+z4))/aux**2
       tragl_der(3,2,3,3) = (x3* (y1-y2)+x1* (y2-y3)+x2* (-y1+y3))* (x4* (-z1+z2)+x2* (z1-z4)+x1* (-z2+z4))/aux**2
       tragl_der(3,2,3,4) = (x3* (y1-y2)+x1* (y2-y3)+x2* (-y1+y3))* (x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))/aux**2
       
       ! dtragl(3,3)       
       tragl_der(3,3,1,1) = (x4* (-y2+y3)+x3* (y2-y4)+x2* (-y3+y4))* (y3* (-z1+z2)+y2* (z1-z3)+y1* (-z2+z3))/aux**2
       tragl_der(3,3,1,2) = (x4* (-y1+y3)+x3* (y1-y4)+x1* (-y3+y4))* (y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))/aux**2
       tragl_der(3,3,1,3) = (x4* (y1-y2)+x1* (y2-y4)+x2* (-y1+y4))* (y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))/aux**2
       tragl_der(3,3,1,4) = (x3* (-y1+y2)+x2* (y1-y3)+x1* (-y2+y3))* (y3* (z1-z2)+y1* (z2-z3)+y2* (-z1+z3))/aux**2
       
       tragl_der(3,3,2,1) = (x4* (-y2+y3)+x3* (y2-y4)+x2* (-y3+y4))* (x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))/aux**2
       tragl_der(3,3,2,2) = (x4* (y1-y3)+x1* (y3-y4)+x3* (-y1+y4))* (x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))/aux**2
       tragl_der(3,3,2,3) = (x4* (y1-y2)+x1* (y2-y4)+x2* (-y1+y4))* (x3* (-z1+z2)+x2* (z1-z3)+x1* (-z2+z3))/aux**2
       tragl_der(3,3,2,4) = (x3* (y1-y2)+x1* (y2-y3)+x2* (-y1+y3))* (x3* (z1-z2)+x1* (z2-z3)+x2* (-z1+z3))/aux**2

       tragl_der(3,3,3,1) = -(x3* (y1-y2)+x1* (y2-y3)+x2* (-y1+y3))* (x4* (-y2+y3)+x3* (y2-y4)+x2* (-y3+y4))/aux**2
       tragl_der(3,3,3,2) = (x3* (y1-y2)+x1* (y2-y3)+x2* (-y1+y3))* (x4* (-y1+y3)+x3* (y1-y4)+x1* (-y3+y4))/aux**2
       tragl_der(3,3,3,3) = (x3* (-y1+y2)+x2* (y1-y3)+x1* (-y2+y3))* (x4* (-y1+y2)+x2* (y1-y4)+x1* (-y2+y4))/aux**2
       tragl_der(3,3,3,4) = -(x3* (-y1+y2)+x2* (y1-y3)+x1* (-y2+y3))**2/aux**2
     
     elseif (pnode == 8) then
       call runend('No esta implementado')
     endif
     
     !
     ! Element length HLENG
     !     
     do idime = 1,ndime
       do inode = 1, pnode

         enor0    =   tragl(1,1) * tragl(1,1) &
            &     + tragl(1,2) * tragl(1,2) &
            &     + tragl(1,3) * tragl(1,3)
         hleng(1) = hnatu/sqrt(enor0)
         enor0_der(idime,inode)    =   2.0_rp*tragl(1,1) * tragl_der(1,1,idime,inode) &
            &     + 2.0_rp*tragl(1,2) * tragl_der(1,2,idime,inode) &
            &     + 2.0_rp*tragl(1,3) * tragl_der(1,3,idime,inode)
         hleng_der(idime,inode,1)   =  -0.5_rp*hnatu*enor0_der(idime,inode)/(sqrt(enor0)*enor0)
       
       
         enor0    =   tragl(2,1) * tragl(2,1) &
              &     + tragl(2,2) * tragl(2,2) &
              &     + tragl(2,3) * tragl(2,3)
         hleng(2) = hnatu/sqrt(enor0)
         enor0_der(idime,inode)    =   2.0_rp*tragl(2,1) * tragl_der(2,1,idime,inode) &
              &     + 2.0_rp*tragl(2,2) * tragl_der(2,2,idime,inode) &
              &     + 2.0_rp*tragl(2,3) * tragl_der(2,3,idime,inode)
         hleng_der(idime,inode,2)   =  -0.5_rp*hnatu*enor0_der(idime,inode)/(sqrt(enor0)*enor0)
         
         
         enor0    =   tragl(3,1) * tragl(3,1) &
              &     + tragl(3,2) * tragl(3,2) &
              &     + tragl(3,3) * tragl(3,3) 
         hleng(3) = hnatu/sqrt(enor0) 
         enor0_der(idime,inode)    =   2.0_rp*tragl(3,1) * tragl_der(3,1,idime,inode) &
              &     + 2.0_rp*tragl(3,2) * tragl_der(3,2,idime,inode) &
              &     + 2.0_rp*tragl(3,3) * tragl_der(3,3,idime,inode) 
         hleng_der(idime,inode,3)   =  -0.5_rp*hnatu*enor0_der(idime,inode)/(sqrt(enor0)*enor0)
         
         !
         ! Sort hleng: hleng(1)=max; hleng(ndime)=min
         !     
         if( hleng(2) > hleng(1) ) then
           h_tem_der(idime,inode)   = hleng_der(idime,inode,2)
           hleng_der(idime,inode,2) = hleng_der(idime,inode,1)
           hleng_der(idime,inode,1) = h_tem_der(idime,inode)
         end if
         if( hleng(3) > hleng(1) ) then
           h_tem_der(idime,inode)   = hleng_der(idime,inode,3)
           hleng_der(idime,inode,3) = hleng_der(idime,inode,1)
           hleng_der(idime,inode,1) = h_tem_der(idime,inode)
         end if
         if( hleng(3) > hleng(2) ) then
           h_tem_der(idime,inode)   = hleng_der(idime,inode,3)
           hleng_der(idime,inode,3) = hleng_der(idime,inode,2)
           hleng_der(idime,inode,2) = h_tem_der(idime,inode)
         end if
              
       enddo
     enddo

  end if

end subroutine elmlen_der
