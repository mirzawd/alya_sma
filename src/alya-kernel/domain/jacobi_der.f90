!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    jacobi_der.f90
!> @author  Mohammad
!> @date    27/05/2016
!> @brief   Computes the spatial derivatives of jacobi_deran, jacobi_deran determinant and
!           Cartesian derivatives of shape function of an element
!> @details The jacobi_deran is
!                              _           _                    
!                             | dx/ds dx/dt |                t
!           jacobi_deran: XJACM = |             | = ELCOD * DERIV 
!                             |_dy/ds dy/dt_|
!           with
!                   _        _             _                    _
!                  | x1 x2 x3 |           | dN1/ds dN2/ds dN3/ds |
!          ELCOD = |          |,  DERIV = |                      |
!                  |_y1 y2 y3_|           |_dN1/dt dN2/dt dN3/dt_|
!
!           => jacobi_deran determinant: GPDET = det(XJACM)  
!!
!> @} 
!-----------------------------------------------------------------------

subroutine jacobi_der(ndime,pnode,elcod,deriv,gpdet_der,gpcar_der)

  use def_kintyp, only : ip,rp
  implicit none
  integer(ip), intent(in)  :: ndime,pnode
  real(rp),    intent(in)  :: elcod(ndime,pnode),deriv(ndime,pnode)
  real(rp)                 :: xjacm_der(ndime,ndime,pnode)
  real(rp),    intent(out) :: gpdet_der(ndime,pnode),gpcar_der(ndime,pnode,ndime,pnode)
  integer(ip)              :: j,k,i
  real(rp)                 :: denom,xjacm(ndime,ndime),denom_der(ndime,pnode),gpdet,xjaci_der(ndime,ndime,ndime,pnode)
!  real(rp)                 :: t1,t2,t3

  if( ndime == 2 .and. pnode == 3 ) then

  else if( ndime == 3 .and. pnode == 4 ) then

  else if ( ndime == 1 ) then
      call runend('No esta implementado') 
!      !
!      ! 1D
!      !
!      xjacm(1,1) = 0.0_rp
!      do k = 1,pnode
!         xjacm(1,1) = xjacm(1,1) + elcod(1,k) * deriv(1,k)
!      end do
!      gpdet = xjacm(1,1)
!      if( gpdet == 0.0_rp ) return
!      xjaci(1,1) = 1.0_rp/xjacm(1,1)
!      do j = 1,pnode
!         gpcar(1,j) = xjaci(1,1) * deriv(1,j)
!      end do

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
     
     do k = 1,pnode
        xjacm_der(1,1,k) = deriv(1,k)
        xjacm_der(1,2,k) = deriv(2,k)
        xjacm_der(2,1,k) = deriv(1,k)
        xjacm_der(2,2,k) = deriv(2,k)
     end do

     gpdet = xjacm(1,1) * xjacm(2,2) - xjacm(2,1) * xjacm(1,2)
     
     do k = 1,pnode
        gpdet_der(1,k) = xjacm_der(1,1,k) * xjacm(2,2) - xjacm(2,1) * xjacm_der(1,2,k)
        gpdet_der(2,k) = xjacm(1,1) * xjacm_der(2,2,k) - xjacm_der(2,1,k) * xjacm(1,2)
     enddo
     
     denom = 1.0_rp/gpdet
     do k = 1,pnode
       denom_der(1,k) = -gpdet_der(1,k)/(gpdet*gpdet)
       denom_der(2,k) = -gpdet_der(2,k)/(gpdet*gpdet)
     enddo
     
!      xjaci(1,1) =  xjacm(2,2) * denom
!      xjaci(2,2) =  xjacm(1,1) * denom
!      xjaci(2,1) = -xjacm(2,1) * denom
!      xjaci(1,2) = -xjacm(1,2) * denom

     do k = 1,pnode
       xjaci_der(1,1,1,k) =  denom_der(1,k)*xjacm(2,2)
       xjaci_der(1,1,2,k) =  denom_der(2,k)*xjacm(2,2) + denom*xjacm_der(2,2,k)
       
       xjaci_der(2,2,1,k) =  denom_der(1,k)*xjacm(1,1) + denom*xjacm_der(1,1,k)
       xjaci_der(2,2,2,k) =  denom_der(2,k)*xjacm(1,1)
       
       xjaci_der(1,2,1,k) =  -denom_der(1,k)*xjacm(1,2) - denom*xjacm_der(1,2,k)
       xjaci_der(1,2,2,k) =  -denom_der(2,k)*xjacm(1,2)

       xjaci_der(2,1,1,k) =  -denom_der(1,k)*xjacm(2,1)
       xjaci_der(2,1,2,k) =  -denom_der(2,k)*xjacm(2,1) - denom*xjacm_der(2,1,k)
       
     end do
! 
!      do j = 1, pnode
!         gpcar(1,j) =   xjaci(1,1) * deriv(1,j) &
!              &       + xjaci(2,1) * deriv(2,j)
! 
!         gpcar(2,j) =   xjaci(1,2) * deriv(1,j) &
!              &       + xjaci(2,2) * deriv(2,j)
!      end do


     do j = 1, pnode
       do i = 1,ndime
         do k = 1, pnode
           gpcar_der(1,j,i,k) = xjaci_der(1,1,i,k) * deriv(1,j) + xjaci_der(2,1,i,k) * deriv(2,j)
           gpcar_der(2,j,i,k) = xjaci_der(1,2,i,k) * deriv(1,j) + xjaci_der(2,2,i,k) * deriv(2,j)       
         enddo
       enddo
     end do



  else if ( ndime == 3 ) then
     !
     ! 3D
     !
     xjacm(1,1) = 0.0_rp ! xjacm = elcod * deriv^t
     xjacm(1,2) = 0.0_rp ! xjaci = xjacm^-1
     xjacm(1,3) = 0.0_rp ! gpcar = xjaci^t * deriv 
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

     do k = 1,pnode
        xjacm_der(1,1,k) = deriv(1,k)
        xjacm_der(1,2,k) = deriv(2,k)
        xjacm_der(1,3,k) = deriv(3,k)
        xjacm_der(2,1,k) = deriv(1,k)
        xjacm_der(2,2,k) = deriv(2,k)
        xjacm_der(2,3,k) = deriv(3,k)
        xjacm_der(3,1,k) = deriv(1,k)
        xjacm_der(3,2,k) = deriv(2,k)
        xjacm_der(3,3,k) = deriv(3,k)
     end do
     

     do k = 1,pnode
        gpdet_der(1,k) = xjacm_der(1,1,k) * ( xjacm(2,2) * xjacm(3,3) - xjacm(2,3) * xjacm(3,2) ) + &
                         xjacm_der(1,2,k) * ( xjacm(2,3) * xjacm(3,1) - xjacm(2,1) * xjacm(3,3) ) + &
                         xjacm_der(1,3,k) * ( xjacm(2,1) * xjacm(3,2) - xjacm(2,2) * xjacm(3,1) )

        gpdet_der(2,k) = xjacm_der(2,1,k) * ( xjacm(1,3) * xjacm(3,2) - xjacm(1,2) * xjacm(3,3) ) + &
                         xjacm_der(2,2,k) * ( xjacm(1,1) * xjacm(3,3) - xjacm(1,3) * xjacm(3,1) ) + &
                         xjacm_der(2,3,k) * ( xjacm(1,2) * xjacm(3,1) - xjacm(1,1) * xjacm(3,2) )
                         
        gpdet_der(3,k) = xjacm_der(3,1,k) * ( xjacm(1,2) * xjacm(2,3) - xjacm(1,3) * xjacm(2,2) ) + &
                         xjacm_der(3,2,k) * ( xjacm(1,3) * xjacm(2,1) - xjacm(1,1) * xjacm(2,3) ) + &
                         xjacm_der(3,3,k) * ( xjacm(1,1) * xjacm(2,2) - xjacm(1,2) * xjacm(2,1) )
     enddo
     
     
!      t1    =  xjacm(2,2) * xjacm(3,3) - xjacm(3,2) * xjacm(2,3)
!      t2    = -xjacm(2,1) * xjacm(3,3) + xjacm(3,1) * xjacm(2,3)
!      t3    =  xjacm(2,1) * xjacm(3,2) - xjacm(3,1) * xjacm(2,2)
!      
!      
!      
!      gpdet =  xjacm(1,1) * t1 + xjacm(1,2) * t2 + xjacm(1,3) * t3
!      if(gpdet == 0.0_rp ) return
!      denom = 1.0_rp / gpdet
!      xjaci(1,1) = t1*denom
!      xjaci(2,1) = t2*denom
!      xjaci(3,1) = t3*denom
!      xjaci(2,2) = ( xjacm(1,1) * xjacm(3,3) - xjacm(3,1) * xjacm(1,3)) * denom
!      xjaci(3,2) = (-xjacm(1,1) * xjacm(3,2) + xjacm(1,2) * xjacm(3,1)) * denom
!      xjaci(3,3) = ( xjacm(1,1) * xjacm(2,2) - xjacm(2,1) * xjacm(1,2)) * denom
!      xjaci(1,2) = (-xjacm(1,2) * xjacm(3,3) + xjacm(3,2) * xjacm(1,3)) * denom
!      xjaci(1,3) = ( xjacm(1,2) * xjacm(2,3) - xjacm(2,2) * xjacm(1,3)) * denom
!      xjaci(2,3) = (-xjacm(1,1) * xjacm(2,3) + xjacm(2,1) * xjacm(1,3)) * denom
! 
!      do j = 1, pnode
!         gpcar(1,j) =   xjaci(1,1) * deriv(1,j) &
!              &       + xjaci(2,1) * deriv(2,j) &
!              &       + xjaci(3,1) * deriv(3,j)
! 
!         gpcar(2,j) =   xjaci(1,2) * deriv(1,j) &
!              &       + xjaci(2,2) * deriv(2,j) &
!              &       + xjaci(3,2) * deriv(3,j)
! 
!         gpcar(3,j) =   xjaci(1,3) * deriv(1,j) &
!              &       + xjaci(2,3) * deriv(2,j) &
!              &       + xjaci(3,3) * deriv(3,j)
!      end do

  end if

end subroutine jacobi_der
