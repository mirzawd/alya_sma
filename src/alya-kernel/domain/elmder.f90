!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine elmder(pnode,ndime,deriv,elcod,gpcar,gpdet,xjacm,xjaci)
  !-----------------------------------------------------------------------
  !****f* domain/elmder
  ! NAME
  !    elmder
  ! DESCRIPTION
  !    This routine calculates the Cartesian derivatives GPCAR and 
  !  the determinant of the Jacobian GPDET.
  ! USES
  !    invmtx
  ! USED BY
  !    ***_elmope
  !    extnor
  ! SOURCE
  !-----------------------------------------------------------------------
  use def_parame, only     :  twopi
  use def_kintyp, only     :  ip,rp
  implicit none
  integer(ip), intent(in)  :: pnode,ndime
  real(rp),    intent(in)  :: deriv(ndime,pnode),elcod(ndime,pnode)
  real(rp),    intent(out) :: gpdet,gpcar(ndime,pnode)
  real(rp),    intent(out) :: xjacm(ndime,ndime)
  real(rp),    intent(out) :: xjaci(ndime,ndime)
  integer(ip)              :: j
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
     do j = 1,pnode
        xjacm(1,1) = xjacm(1,1) + elcod(1,j) * deriv(1,j)
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
     do j = 1,pnode
        xjacm(1,1) = xjacm(1,1) + elcod(1,j) * deriv(1,j)
        xjacm(1,2) = xjacm(1,2) + elcod(1,j) * deriv(2,j)
        xjacm(2,1) = xjacm(2,1) + elcod(2,j) * deriv(1,j)
        xjacm(2,2) = xjacm(2,2) + elcod(2,j) * deriv(2,j)
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
     xjacm(1,1) = 0.0_rp ! xjacm = elcod * deriv^t
     xjacm(1,2) = 0.0_rp ! xjaci = xjacm^-1
     xjacm(1,3) = 0.0_rp ! gpcar = xjaci^t * deriv 
     xjacm(2,1) = 0.0_rp
     xjacm(2,2) = 0.0_rp
     xjacm(2,3) = 0.0_rp
     xjacm(3,1) = 0.0_rp
     xjacm(3,2) = 0.0_rp
     xjacm(3,3) = 0.0_rp
     do j = 1,pnode
        xjacm(1,1) = xjacm(1,1) + elcod(1,j) * deriv(1,j)
        xjacm(1,2) = xjacm(1,2) + elcod(1,j) * deriv(2,j)
        xjacm(1,3) = xjacm(1,3) + elcod(1,j) * deriv(3,j)
        xjacm(2,1) = xjacm(2,1) + elcod(2,j) * deriv(1,j)
        xjacm(2,2) = xjacm(2,2) + elcod(2,j) * deriv(2,j)
        xjacm(2,3) = xjacm(2,3) + elcod(2,j) * deriv(3,j)
        xjacm(3,1) = xjacm(3,1) + elcod(3,j) * deriv(1,j)
        xjacm(3,2) = xjacm(3,2) + elcod(3,j) * deriv(2,j)
        xjacm(3,3) = xjacm(3,3) + elcod(3,j) * deriv(3,j)
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

end subroutine elmder
!-----------------------------------------------------------------------
! NOTES
! Try to remember... when life was so tender...
!        _        _           _                    _
!       | x1 x2 x3 |         | dN1/ds dN2/ds dN3/ds |
! elcod=|          |,  deriv=|                      |
!       |_y1 y2 y3_|         |_dN1/dt dN2/dt dN3/dt_|
!                  _           _                        
!                 | dx/ds dx/dt |            t
! Jacobian: xjacm=|             |=elcod*deriv 
!                 |_dy/ds dy/dt_|
!
! => Jacobian determinant = det(xjacm)
!
! o xjacm^t = deriv * elcod^t
!        _           _
!       | ds/dx ds/dy |      -1
! xjaci=|             |=xjacm  .This is because:
!       |_dt/dx dt/dy_|
!        _     _   _           _   _     _ 
!       | du/dx | | ds/dx dt/dx | | du/ds |    
!       |       |=|             | |       | and
!       |_du/dy_| |_ds/dy dt/dy_| |_du/dt_|     
!        _     _   _           _   _     _      _     _   _           _ -1  _     _ 
!       | du/ds | | dx/ds dy/ds | | du/dx |    | du/dx | | dx/ds dy/ds |   | du/ds |
!       |       |=|             | |       | => |       |=|             |   |       |
!       |_du/dt_| |_dx/dt dy/dt_| |_du/dy_|    |_du/dy_| |_dx/dt dy/dt_|   |_du/dt_| 
!                  _           _   _           _ -1
!                 | ds/dx dt/dx | | dx/ds dy/ds |      -1
!       Therefore |             |=|             |=xjacm  
!                 |_ds/dy dt/dy_| |_dx/dt dy/dt_|  
!        _                                                                               _    
!       | dN1/ds*ds/dx+dN1/dt*dt/dx  dN2/ds*ds/dx+dN2/dt*dt/dx dN3/ds*ds/dx+dN3/dt*dt/dx  |   
! gpcar=|                                                                                 |   
!       |_dN1/ds*ds/dy+dN1/dt*dt/dy  dN2/ds*ds/dy+dN2/dt*dt/dy dN3/ds*ds/dy+dN3/dt*dt/dy _|
!            t
!      =xjaci *deriv
!           _     _
!          | du/dx |                 t
! Example: |       |=gpcar*[U1 U2 U3]
!          |_du/dy_|        
!     
!***
!-----------------------------------------------------------------------
