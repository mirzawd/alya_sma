!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Adaptivity
!> @{
!> @file    mod_objectiveFun.f90
!> @author  abel.gargallo
!> @date    2021-06-21
!> @brief   mod_objectiveFun
!> @details mod_objectiveFun
!>
!>          To add further details
!>
!>
!-----------------------------------------------------------------------

MODULE mod_objectiveFun
  
  
use def_kintyp_basic, only : ip,rp,lg

use mod_memory,             only: memory_alloca, memory_deallo
use def_adapt,              only: memor_adapt

use def_kintyp_mesh_basic,  only: mesh_type_basic
use mod_metric,             only: mesh_metric_type

implicit none

private

public :: set_fun_meshData, delete_fun_meshData
public :: set_fun_nodeData, delete_fun_nodeData
public :: f_node_distortion!, fgrad_node_distortion, fhess_node_distortion

type(mesh_type_basic),  pointer  :: theMesh
type(mesh_metric_type), pointer  :: theMetric
!
integer(ip)                      :: theNode
integer(ip),            pointer  :: adjElems(:)
!
real(rp), parameter :: delta_regularization = 0.0001_rp ! to compute safe numerical derivatives of distortion objective function
!
save
!  
!
!
CONTAINS
!
!
!
subroutine set_fun_meshData(mesh,metric)
  implicit none
  type(mesh_type_basic), target, intent(inout)  :: mesh
  type(mesh_metric_type),target, intent(inout)  :: metric
  !
  theMesh   => mesh
  theMetric => metric
  !
end subroutine set_fun_meshData
!
!
!
subroutine delete_fun_meshData()
  implicit none
  nullify(theMesh  )
  nullify(theMetric)
end subroutine delete_fun_meshData
!
!
!
subroutine set_fun_nodeData(inode,list_adjElems)
  implicit none
  integer(ip), intent(in)    :: inode
  integer(ip), intent(inout), target :: list_adjElems(:)
  theNode   =  inode
  adjElems  => list_adjElems
end subroutine set_fun_nodeData
!
!
!
subroutine delete_fun_nodeData()
  implicit none
  theNode   = 0_ip
  nullify(adjElems)
end subroutine delete_fun_nodeData
!
!
!
function f_node_distortion(dim,x) result(f)
  use mod_quality, only: compute_distortion_sizeshape_subset
  implicit none
  integer(ip), intent(in) :: dim
  real(rp),    intent(in) :: x(dim)
  !
  real(rp) :: f
  real(rp) :: x_safe(dim)
  
  x_safe = theMesh%coord(:,theNode)
  theMesh%coord(:,theNode) = x
  
  f = compute_distortion_sizeshape_subset(theMesh,theMetric,adjElems,delta_regularization)
  
  theMesh%coord(:,theNode) = x_safe
  
  return
end function f_node_distortion
!
!
!
END MODULE mod_objectiveFun

!> @}





!
!
!
! function fgrad_node_distortion(x) result(df)
!   use mod_quality, only: compute_distortion_sizeshape_subset
!   real(rp), intent(in) :: x(:)
!   real(rp) :: df(size(x))
!
!   integer(ip) :: idim
!   real(rp)    :: eps3, typ, ei, temp
!   real (rp)   :: x_min(size(x)),x_plus(size(x))
!   real(rp)    :: f_plus, f_min
!
!   eps3 = 4.64e-6_rp ! 1.0e-4 ! 4.64e-6
!   typ  = 0.5_rp
!
!   do idim=1,size(x)
!
!     ei   = eps3*max( abs(x(idim)),  typ )!*sign(1.0,x(1))  ! sign is not necessary (centered)
!     temp = x(idim) +  ei
!     ei   = temp -  x(idim)
!
!     x_min=x
!     x_min(idim)=x_min(idim)-ei
!     f_min= f_node_distortion(x_min)
!
!     x_plus=x
!     x_plus(idim)=x_plus(idim)+ei
!     f_plus= f_node_distortion(x_plus)
!
!     df(idim) = (f_plus-f_min)/(2.0_rp*ei)
!
!   end do
!
! end function fgrad_node_distortion
! !
! !
! !
! function fhess_node_distortion(x) result(H)
!   use mod_quality, only: compute_distortion_sizeshape_subset
!   real(rp), intent(in) :: x(:)
!
!   real(rp) :: H(size(x),size(x))
!
!   integer(ip) :: idim, jdim
!   real(rp)    :: eps3, typ, ei, ej, temp
!   real (rp)   :: x_aux(size(x))
!   real(rp)    :: fPiPj,fPiMj,fMiPj,fMiMj
!
!   eps3 = 4.64e-6_rp ! 1.0e-4 ! 4.64e-6
!   typ  = 0.5_rp
!
!   do idim=1,size(x)
!
!     ei   = eps3*max( abs(x(idim)),  typ )!*sign(1.0,x(1))  ! sign is not necessary (centered)
!     temp = x(idim) +  ei
!     ei   = temp -  x(idim)
!
!     do jdim=idim,size(x)
!
!       ej   = eps3*max( abs(x(idim)),  typ )!*sign(1.0,x(1))  ! sign is not necessary (centered)
!       temp = x(idim) +  ej
!       ej   = temp -  x(idim)
!
!       x_aux       = x
!       x_aux(idim) = x(idim)+ei
!       x_aux(jdim) = x(jdim)+ej
!       fPiPj       = f_node_distortion(x_aux)
!       x_aux(idim) = x(idim)+ei
!       x_aux(jdim) = x(jdim)-ej
!       fPiMj       = f_node_distortion(x_aux)
!       x_aux(idim) = x(idim)-ei
!       x_aux(jdim) = x(jdim)+ej
!       fMiPj       = f_node_distortion(x_aux)
!       x_aux(idim) = x(idim)-ei
!       x_aux(jdim) = x(jdim)-ej
!       fMiMj       = f_node_distortion(x_aux)
!
!       H(idim,jdim) = (fPiPj-fPiMj-fMiPj+fMiMj)/(4_rp*ei*ej)
!       H(jdim,idim) = H(idim,jdim)
!     end do
!
!   end do
!
! end function fhess_node_distortion