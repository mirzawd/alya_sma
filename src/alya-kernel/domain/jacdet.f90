!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    jacdet.f90
!> @author  Guillaume Houzeaux
!> @date    27/01/2007
!> @brief   Compute the Jacobian determinant of an element
!> @details The jacobian is
!!          \verbatim
!!                             _           _                    
!!                            | dx/ds dx/dt |                t
!!          Jacobian: XJACM = |             | = ELCOD * DERIV 
!!                            |_dy/ds dy/dt_|
!!          with
!!                   _        _             _                    _
!!                  | x1 x2 x3 |           | dN1/ds dN2/ds dN3/ds |
!!          ELCOD = |          |,  DERIV = |                      |
!!                  |_y1 y2 y3_|           |_dN1/dt dN2/dt dN3/dt_|
!!
!!          => Jacobian determinant: GPDET = det(XJACM)  
!!
!!          \endverbatim
!!
!> @} 
!-----------------------------------------------------------------------

subroutine jacdet(ndime,pnode,elcod,deriv,xjacm,gpdet)

  use def_kintyp, only : ip,rp
  implicit none
  integer(ip), intent(in)  :: ndime                          !< Space dimension
  integer(ip), intent(in)  :: pnode                          !< Element number of nodes
  real(rp),    intent(in)  :: elcod(ndime,pnode)             !< Element coordinates
  real(rp),    intent(in)  :: deriv(ndime,pnode)             !< Shape function derivatives
  real(rp),    intent(out) :: xjacm(ndime,ndime)             !< Jacobian
  real(rp),    intent(out) :: gpdet                          !> Jacobian determinant
  integer(ip)              :: inode
  real(rp)                 :: t1,t2,t3

  if( ndime == 2 .and. pnode == 3 ) then
     !
     ! 2D P1 element
     !
     gpdet =  ( -elcod(1,1) + elcod(1,2) ) * ( -elcod(2,1) + elcod(2,3) ) &
          & - ( -elcod(2,1) + elcod(2,2) ) * ( -elcod(1,1) + elcod(1,3) )

  else if( ndime == 3 .and. pnode == 4 ) then
     !
     ! 3D P1 element
     !
     xjacm(1,1) =  elcod(1,2) - elcod(1,1)
     xjacm(1,2) =  elcod(1,3) - elcod(1,1)
     xjacm(1,3) =  elcod(1,4) - elcod(1,1)
     xjacm(2,1) =  elcod(2,2) - elcod(2,1)
     xjacm(2,2) =  elcod(2,3) - elcod(2,1)
     xjacm(2,3) =  elcod(2,4) - elcod(2,1)
     xjacm(3,1) =  elcod(3,2) - elcod(3,1)
     xjacm(3,2) =  elcod(3,3) - elcod(3,1)
     xjacm(3,3) =  elcod(3,4) - elcod(3,1)
     t1         =  xjacm(2,2) * xjacm(3,3) - xjacm(3,2) * xjacm(2,3)
     t2         = -xjacm(2,1) * xjacm(3,3) + xjacm(3,1) * xjacm(2,3)
     t3         =  xjacm(2,1) * xjacm(3,2) - xjacm(3,1) * xjacm(2,2)
     gpdet      =  xjacm(1,1) * t1 + xjacm(1,2) * t2 + xjacm(1,3) * t3

  else if ( ndime == 3 ) then
     !
     ! 3D
     !
     xjacm(1,1) = 0.0_rp
     xjacm(1,2) = 0.0_rp
     xjacm(1,3) = 0.0_rp
     xjacm(2,1) = 0.0_rp
     xjacm(2,2) = 0.0_rp
     xjacm(2,3) = 0.0_rp
     xjacm(3,1) = 0.0_rp 
     xjacm(3,2) = 0.0_rp
     xjacm(3,3) = 0.0_rp
     do inode = 1,pnode
        xjacm(1,1) = xjacm(1,1) + elcod(1,inode) * deriv(1,inode)
        xjacm(1,2) = xjacm(1,2) + elcod(1,inode) * deriv(2,inode)
        xjacm(1,3) = xjacm(1,3) + elcod(1,inode) * deriv(3,inode)
        xjacm(2,1) = xjacm(2,1) + elcod(2,inode) * deriv(1,inode)
        xjacm(2,2) = xjacm(2,2) + elcod(2,inode) * deriv(2,inode)
        xjacm(2,3) = xjacm(2,3) + elcod(2,inode) * deriv(3,inode)
        xjacm(3,1) = xjacm(3,1) + elcod(3,inode) * deriv(1,inode)
        xjacm(3,2) = xjacm(3,2) + elcod(3,inode) * deriv(2,inode)
        xjacm(3,3) = xjacm(3,3) + elcod(3,inode) * deriv(3,inode)
     end do

     t1    =  xjacm(2,2) * xjacm(3,3) - xjacm(3,2) * xjacm(2,3)
     t2    = -xjacm(2,1) * xjacm(3,3) + xjacm(3,1) * xjacm(2,3)
     t3    =  xjacm(2,1) * xjacm(3,2) - xjacm(3,1) * xjacm(2,2)
     gpdet =  xjacm(1,1) * t1 + xjacm(1,2) * t2 + xjacm(1,3) * t3
     
  else if ( ndime == 2 ) then
     !
     ! 2D
     !
     xjacm(1,1) = 0.0_rp
     xjacm(1,2) = 0.0_rp
     xjacm(2,1) = 0.0_rp
     xjacm(2,2) = 0.0_rp
     do inode = 1,pnode
        xjacm(1,1) = xjacm(1,1) + elcod(1,inode) * deriv(1,inode)
        xjacm(1,2) = xjacm(1,2) + elcod(1,inode) * deriv(2,inode)
        xjacm(2,1) = xjacm(2,1) + elcod(2,inode) * deriv(1,inode)
        xjacm(2,2) = xjacm(2,2) + elcod(2,inode) * deriv(2,inode)
     end do
 
     gpdet = xjacm(1,1) * xjacm(2,2) - xjacm(2,1) * xjacm(1,2) 

  else if ( ndime == 1 ) then
     !
     ! 1D
     !
     gpdet = 0.0_rp
     do inode = 1,pnode
        gpdet = gpdet + elcod(1,inode) * deriv(1,inode)
     end do

  end if

end subroutine jacdet
