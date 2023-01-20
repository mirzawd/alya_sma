!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine elmdel(pnode,ndime,elcod,gpcar,gpdet)
  !-----------------------------------------------------------------------
  !****f* Domain/elmdel
  ! NAME
  !    elmdel
  ! DESCRIPTION
  !    This routine calculates the Cartesian derivatives for 2D and
  !    3D P1 element.
  ! USES
  ! USED BY
  ! SOURCE
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  implicit none
  integer(ip), intent(in)  :: pnode,ndime
  real(rp),    intent(in)  :: elcod(ndime,pnode)
  real(rp),    intent(out) :: gpdet,gpcar(ndime,pnode)
  real(rp)                 :: b(3,3),t1,t2,t3,denom

  if( ndime == 2 .and. pnode == 3 ) then
     !
     ! 2D P1 element
     !
     gpdet      =  (-elcod(1,1)+elcod(1,2))*(-elcod(2,1)+elcod(2,3))&
          &       -(-elcod(2,1)+elcod(2,2))*(-elcod(1,1)+elcod(1,3))

     denom      =  1.0_rp/gpdet
     gpcar(1,1) =  ( -elcod(2,3) + elcod(2,2) ) * denom
     gpcar(1,2) =  ( -elcod(2,1) + elcod(2,3) ) * denom
     gpcar(1,3) =  (  elcod(2,1) - elcod(2,2) ) * denom
     gpcar(2,1) =  (  elcod(1,3) - elcod(1,2) ) * denom
     gpcar(2,2) =  (  elcod(1,1) - elcod(1,3) ) * denom
     gpcar(2,3) =  ( -elcod(1,1) + elcod(1,2) ) * denom

  else if( ndime == 3 .and. pnode == 4 ) then
     !
     ! 3D P1 element
     !
     gpcar(1,1) = elcod(1,2) - elcod(1,1)
     gpcar(1,2) = elcod(1,3) - elcod(1,1)
     gpcar(1,3) = elcod(1,4) - elcod(1,1)
     gpcar(2,1) = elcod(2,2) - elcod(2,1)
     gpcar(2,2) = elcod(2,3) - elcod(2,1)
     gpcar(2,3) = elcod(2,4) - elcod(2,1)
     gpcar(3,1) = elcod(3,2) - elcod(3,1)
     gpcar(3,2) = elcod(3,3) - elcod(3,1)
     gpcar(3,3) = elcod(3,4) - elcod(3,1)

     t1         = gpcar(2,2)*gpcar(3,3) - gpcar(3,2)*gpcar(2,3)
     t2         =-gpcar(2,1)*gpcar(3,3) + gpcar(3,1)*gpcar(2,3)
     t3         = gpcar(2,1)*gpcar(3,2) - gpcar(3,1)*gpcar(2,2)
     gpdet      = (gpcar(1,1)*t1 + gpcar(1,2)*t2 + gpcar(1,3)*t3)
     denom      = 1.0_rp/gpdet

     b(1,1)     = t1*denom
     b(2,1)     = t2*denom
     b(3,1)     = t3*denom
     b(2,2)     = ( gpcar(1,1) * gpcar(3,3) - gpcar(3,1) * gpcar(1,3)) * denom
     b(3,2)     = (-gpcar(1,1) * gpcar(3,2) + gpcar(1,2) * gpcar(3,1)) * denom
     b(3,3)     = ( gpcar(1,1) * gpcar(2,2) - gpcar(2,1) * gpcar(1,2)) * denom
     b(1,2)     = (-gpcar(1,2) * gpcar(3,3) + gpcar(3,2) * gpcar(1,3)) * denom
     b(1,3)     = ( gpcar(1,2) * gpcar(2,3) - gpcar(2,2) * gpcar(1,3)) * denom
     b(2,3)     = (-gpcar(1,1) * gpcar(2,3) + gpcar(2,1) * gpcar(1,3)) * denom

     gpcar(1,1) = -b(1,1)-b(2,1)-b(3,1)
     gpcar(1,2) =  b(1,1)
     gpcar(1,3) =  b(2,1)
     gpcar(1,4) =  b(3,1)
     gpcar(2,1) = -b(1,2)-b(2,2)-b(3,2)
     gpcar(2,2) =  b(1,2)
     gpcar(2,3) =  b(2,2)
     gpcar(2,4) =  b(3,2)
     gpcar(3,1) = -b(1,3)-b(2,3)-b(3,3)
     gpcar(3,2) =  b(1,3)
     gpcar(3,3) =  b(2,3)
     gpcar(3,4) =  b(3,3)

  end if

end subroutine elmdel 
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
