!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    jacobi.f90
!> @author  Guillaume Houzeaux
!> @date    27/01/2007
!> @brief   Computes the Jacobian, Jacobian determinant and
!           Cartesian derivatives of shape function of an element
!> @details The jacobian is
!                              _           _                    
!                             | dx/ds dx/dt |                t
!           Jacobian: XJACM = |             | = ELCOD * DERIV 
!                             |_dy/ds dy/dt_|
!           with
!                   _        _             _                    _
!                  | x1 x2 x3 |           | dN1/ds dN2/ds dN3/ds |
!          ELCOD = |          |,  DERIV = |                      |
!                  |_y1 y2 y3_|           |_dN1/dt dN2/dt dN3/dt_|
!
!           => Jacobian determinant: GPDET = det(XJACM)  
!!
!> @} 
!-----------------------------------------------------------------------

subroutine jacobi(ndime,pnode,elcod,deriv,xjacm,xjaci,gpcar,gpdet)

  use def_kintyp, only : ip,rp
  implicit none
  integer(ip), intent(in)  :: ndime,pnode
  real(rp),    intent(in)  :: elcod(ndime,pnode),deriv(ndime,pnode)
  real(rp),    intent(out) :: xjacm(ndime,ndime),xjaci(ndime,ndime)
  real(rp),    intent(out) :: gpcar(ndime,pnode),gpdet
  integer(ip)              :: j,k
  real(rp)                 :: t1,t2,t3,denom

  if( ndime == 2 .and. pnode == 3 ) then
     !
     ! 2D P1 element
     !
     gpdet = (-elcod(1,1)+elcod(1,2))*(-elcod(2,1)+elcod(2,3)) &
          & -(-elcod(2,1)+elcod(2,2))*(-elcod(1,1)+elcod(1,3))
     if( gpdet == 0.0_rp ) return     
     denom = 1.0_rp/gpdet

     gpcar(1,1) = ( -elcod(2,3) + elcod(2,2) ) * denom
     gpcar(1,2) = ( -elcod(2,1) + elcod(2,3) ) * denom
     gpcar(1,3) = (  elcod(2,1) - elcod(2,2) ) * denom
     gpcar(2,1) = (  elcod(1,3) - elcod(1,2) ) * denom
     gpcar(2,2) = (  elcod(1,1) - elcod(1,3) ) * denom
     gpcar(2,3) = ( -elcod(1,1) + elcod(1,2) ) * denom


!!$     xjacm(1,1) = 0.0_rp
!!$     xjacm(1,2) = 0.0_rp
!!$     xjacm(2,1) = 0.0_rp
!!$     xjacm(2,2) = 0.0_rp
!!$     do k = 1,pnode
!!$        xjacm(1,1) = xjacm(1,1) + elcod(1,k) * deriv(1,k)
!!$        xjacm(1,2) = xjacm(1,2) + elcod(1,k) * deriv(2,k)
!!$        xjacm(2,1) = xjacm(2,1) + elcod(2,k) * deriv(1,k)
!!$        xjacm(2,2) = xjacm(2,2) + elcod(2,k) * deriv(2,k)
!!$     end do
!!$
!!$     gpdet = xjacm(1,1) * xjacm(2,2) - xjacm(2,1) * xjacm(1,2)
!!$     if( gpdet == 0.0_rp ) return
!!$     denom = 1.0_rp/gpdet
!!$     xjaci(1,1) =  xjacm(2,2) * denom
!!$     xjaci(2,2) =  xjacm(1,1) * denom
!!$     xjaci(2,1) = -xjacm(2,1) * denom
!!$     xjaci(1,2) = -xjacm(1,2) * denom  




  else if( ndime == 3 .and. pnode == 4 ) then
     !
     ! 3D P1 element
     !
     gpcar(1,1) =  elcod(1,2) - elcod(1,1)
     gpcar(1,2) =  elcod(1,3) - elcod(1,1)
     gpcar(1,3) =  elcod(1,4) - elcod(1,1)
     gpcar(2,1) =  elcod(2,2) - elcod(2,1)
     gpcar(2,2) =  elcod(2,3) - elcod(2,1)
     gpcar(2,3) =  elcod(2,4) - elcod(2,1)
     gpcar(3,1) =  elcod(3,2) - elcod(3,1)
     gpcar(3,2) =  elcod(3,3) - elcod(3,1)
     gpcar(3,3) =  elcod(3,4) - elcod(3,1)
     t1         =  gpcar(2,2) * gpcar(3,3) - gpcar(3,2) * gpcar(2,3)
     t2         = -gpcar(2,1) * gpcar(3,3) + gpcar(3,1) * gpcar(2,3)
     t3         =  gpcar(2,1) * gpcar(3,2) - gpcar(3,1) * gpcar(2,2)
     gpdet      =  gpcar(1,1) * t1 + gpcar(1,2) * t2 + gpcar(1,3) * t3
     if( gpdet == 0.0_rp ) return
     denom      =  1.0_rp/gpdet

     xjaci(1,1) =  t1 * denom
     xjaci(2,1) =  t2 * denom
     xjaci(3,1) =  t3 * denom
     xjaci(2,2) = ( gpcar(1,1) * gpcar(3,3) - gpcar(3,1) * gpcar(1,3)) * denom
     xjaci(3,2) = (-gpcar(1,1) * gpcar(3,2) + gpcar(1,2) * gpcar(3,1)) * denom
     xjaci(3,3) = ( gpcar(1,1) * gpcar(2,2) - gpcar(2,1) * gpcar(1,2)) * denom
     xjaci(1,2) = (-gpcar(1,2) * gpcar(3,3) + gpcar(3,2) * gpcar(1,3)) * denom
     xjaci(1,3) = ( gpcar(1,2) * gpcar(2,3) - gpcar(2,2) * gpcar(1,3)) * denom
     xjaci(2,3) = (-gpcar(1,1) * gpcar(2,3) + gpcar(2,1) * gpcar(1,3)) * denom
     gpcar(1,1) = -xjaci(1,1)-xjaci(2,1)-xjaci(3,1)
     gpcar(1,2) =  xjaci(1,1)
     gpcar(1,3) =  xjaci(2,1)
     gpcar(1,4) =  xjaci(3,1)
     gpcar(2,1) = -xjaci(1,2)-xjaci(2,2)-xjaci(3,2)
     gpcar(2,2) =  xjaci(1,2)
     gpcar(2,3) =  xjaci(2,2)
     gpcar(2,4) =  xjaci(3,2)
     gpcar(3,1) = -xjaci(1,3)-xjaci(2,3)-xjaci(3,3)
     gpcar(3,2) =  xjaci(1,3)
     gpcar(3,3) =  xjaci(2,3)
     gpcar(3,4) =  xjaci(3,3)

  else if ( ndime == 1 ) then
     !
     ! 1D
     !
     xjacm(1,1) = 0.0_rp
     do k = 1,pnode
        xjacm(1,1) = xjacm(1,1) + elcod(1,k) * deriv(1,k)
     end do
     gpdet = xjacm(1,1)
     if( gpdet == 0.0_rp ) return
     xjaci(1,1) = 1.0_rp/xjacm(1,1)
     do j = 1,pnode
        gpcar(1,j) = xjaci(1,1) * deriv(1,j)
     end do

  else if ( ndime == 2 ) then
     !
     ! 2D
     !
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

     gpdet = xjacm(1,1) * xjacm(2,2) - xjacm(2,1) * xjacm(1,2)
     if( gpdet == 0.0_rp ) return
     denom = 1.0_rp/gpdet
     xjaci(1,1) =  xjacm(2,2) * denom
     xjaci(2,2) =  xjacm(1,1) * denom
     xjaci(2,1) = -xjacm(2,1) * denom
     xjaci(1,2) = -xjacm(1,2) * denom  

     do j = 1, pnode
        gpcar(1,j) =   xjaci(1,1) * deriv(1,j) &
             &       + xjaci(2,1) * deriv(2,j)

        gpcar(2,j) =   xjaci(1,2) * deriv(1,j) &
             &       + xjaci(2,2) * deriv(2,j)
     end do

  else if ( ndime == 3 ) then
     !
     ! 3D
     !
     ! xjacm = elcod * deriv^t 
     ! xjaci = xjacm^-1        
     ! gpcar = xjaci^t * deriv 
     !
     xjacm(1:3,1:3) = 0.0_rp
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

     t1    =  xjacm(2,2) * xjacm(3,3) - xjacm(3,2) * xjacm(2,3)
     t2    = -xjacm(2,1) * xjacm(3,3) + xjacm(3,1) * xjacm(2,3)
     t3    =  xjacm(2,1) * xjacm(3,2) - xjacm(3,1) * xjacm(2,2)
     gpdet =  xjacm(1,1) * t1 + xjacm(1,2) * t2 + xjacm(1,3) * t3
     if(gpdet == 0.0_rp ) return
     denom = 1.0_rp / gpdet
     xjaci(1,1) = t1*denom
     xjaci(2,1) = t2*denom
     xjaci(3,1) = t3*denom
     xjaci(2,2) = ( xjacm(1,1) * xjacm(3,3) - xjacm(3,1) * xjacm(1,3)) * denom
     xjaci(3,2) = (-xjacm(1,1) * xjacm(3,2) + xjacm(1,2) * xjacm(3,1)) * denom
     xjaci(3,3) = ( xjacm(1,1) * xjacm(2,2) - xjacm(2,1) * xjacm(1,2)) * denom
     xjaci(1,2) = (-xjacm(1,2) * xjacm(3,3) + xjacm(3,2) * xjacm(1,3)) * denom
     xjaci(1,3) = ( xjacm(1,2) * xjacm(2,3) - xjacm(2,2) * xjacm(1,3)) * denom
     xjaci(2,3) = (-xjacm(1,1) * xjacm(2,3) + xjacm(2,1) * xjacm(1,3)) * denom

     do j = 1, pnode
        gpcar(1,j) =   xjaci(1,1) * deriv(1,j) &
             &       + xjaci(2,1) * deriv(2,j) &
             &       + xjaci(3,1) * deriv(3,j)

        gpcar(2,j) =   xjaci(1,2) * deriv(1,j) &
             &       + xjaci(2,2) * deriv(2,j) &
             &       + xjaci(3,2) * deriv(3,j)

        gpcar(3,j) =   xjaci(1,3) * deriv(1,j) &
             &       + xjaci(2,3) * deriv(2,j) &
             &       + xjaci(3,3) * deriv(3,j)
     end do

  end if

end subroutine jacobi
