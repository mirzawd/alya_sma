!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Mathru
!> @{
!> @file    rotmat.f90
!> @author  Guillaume Houzeaux
!> @date    20/02/2013
!> @brief   Compute a rotation matrix
!> @details Rotation matrix from
!>          http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
!> @} 
!-----------------------------------------------------------------------

pure subroutine rotmat(ndime,rot_angle,rot_axis,rot_matrix)
  use def_kintyp, only : ip,rp
  implicit none
  integer(ip), intent(in)  :: ndime                       !< Problem dimension
  real(rp),    intent(in)  :: rot_angle                   !< Rot_angle of rotation
  real(rp),    intent(in)  :: rot_axis(3)                 !< Rot_axis of rotation
  real(rp),    intent(out) :: rot_matrix(ndime,ndime)     !< Rotation matrix
  real(rp)                 :: cosQ,sinQ
  real(rp)                 :: oneMinusCosQ
  !
  ! Cos and Sin of rot_angle
  !
  cosQ         = cos( rot_angle )
  sinQ         = sin( rot_angle )
  oneMinusCosQ = 1.0_rp - cosQ
  !
  ! Rotation matrix 
  !
  rot_matrix(1,1) = cosQ + rot_axis(1)**2          * oneMinusCosQ
  rot_matrix(1,2) =        rot_axis(1)*rot_axis(2) * oneMinusCosQ
  rot_matrix(2,1) =        rot_axis(1)*rot_axis(2) * oneMinusCosQ
  rot_matrix(2,2) = cosQ + rot_axis(2)**2          * oneMinusCosQ
  rot_matrix(1,2) =        rot_matrix (1,2)                       - rot_axis(3)*sinQ
  rot_matrix(2,1) =        rot_matrix (2,1)                       + rot_axis(3)*sinQ
  if( ndime == 3 ) then
     rot_matrix(3,1) =        rot_axis(1)*rot_axis(3) * oneMinusCosQ - rot_axis(2)*sinQ
     rot_matrix(2,3) =        rot_axis(2)*rot_axis(3) * oneMinusCosQ - rot_axis(1)*sinQ
     rot_matrix(1,3) =        rot_axis(1)*rot_axis(3) * oneMinusCosQ + rot_axis(2)*sinQ
     rot_matrix(3,2) =        rot_axis(2)*rot_axis(3) * oneMinusCosQ + rot_axis(1)*sinQ
     rot_matrix(3,3) = cosQ + rot_axis(3)**2          * oneMinusCosQ
  end if

end subroutine rotmat
