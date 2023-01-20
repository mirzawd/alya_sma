!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Adaptivity
!> @{
!> @file    mod_cavity.f90
!> @author  abel.gargallo
!> @date    2021-03-31
!> @brief   Cavity operators
!> @details Cavity operators
!>
!>          To add further details
!>
!>
!-----------------------------------------------------------------------

MODULE mod_cavity
!***************************************************************
!*
!*  Module for performing cavity operations
!*
!***************************************************************
use def_kintyp_basic, only: ip,rp,lg
use mod_meshTopology, only: get_boundary_faces, delete_boundary_faces

use mod_memory,             only: memory_alloca, memory_deallo
use def_adapt,              only: memor_adapt

implicit none

private

public :: cavity_insert_point
public :: cavity_reinsert
public :: cavity_reinsert_fromCavBoun
public :: delete_Tcav_remeshed

!save
!  
!
!
CONTAINS
!
!
!
! subroutine cavity_insert_point(newPoint,elems_cav,mesh,M, Tcav_remeshed, elems_cav_updated, is_remeshed)
!   implicit none
!
!   integer(ip),            intent(in) :: newPoint
!   integer(ip), pointer    intent(in) :: elems_cav
!   type(mesh_type_basic),  intent(in) :: mesh
!   type(mesh_metric_type), intent(in) :: M
!   !
!   integer(ip), pointer, intent(inout) :: Tcav_remeshed(:,:)
!   integer(ip), pointer, intent(inout) :: elems_cav_updated(:)
!   logical(lg),          intent(out) :: is_remeshed
subroutine cavity_insert_point(newPoint,Tcav,Tcav_remeshed)
  implicit none
  
  integer(ip), intent(in) :: newPoint
  integer(ip), intent(in) :: Tcav(:,:)
  !
  integer(ip), pointer, intent(inout) :: Tcav_remeshed(:,:)
  !
  integer(ip), pointer :: boundary_faces(:,:)
  integer(ip) :: iface
  !
  call get_boundary_faces(Tcav,boundary_faces)
  
  nullify(Tcav_remeshed)
  call memory_alloca(memor_adapt,'Tcav_remeshed',',mod_cavity',Tcav_remeshed,int(size(Tcav,1),ip),int(size(boundary_faces,2),ip))

  Tcav_remeshed(size(Tcav,1),:) = newPoint
  do iface =1,size(boundary_faces,2)
    Tcav_remeshed(1:size(boundary_faces,1),iface) = boundary_faces(:,iface)
  end do
  
  !call memory_deallo(memor_adapt,'boundary_faces','cavity_insert_point',boundary_faces)
  call delete_boundary_faces(boundary_faces)
  
  return 
end subroutine cavity_insert_point
!
!
!
subroutine delete_Tcav_remeshed(Tcav_remeshed)
  implicit none 
  
  integer(ip), pointer, intent(inout) :: Tcav_remeshed(:,:)
  
  call memory_deallo(memor_adapt,'Tcav_remeshed'  ,'cavity_insert_point',Tcav_remeshed)
  
end subroutine delete_Tcav_remeshed
!
!
!
subroutine cavity_reinsert(pointToReinsert,Tcav,Tcav_remeshed)
  implicit none
  
  integer(ip), intent(in) :: pointToReinsert
  integer(ip), intent(in) :: Tcav(:,:)
!   integer(ip), pointer, optional, intent(in):: boundary_faces_inp(:,:)
  !
  integer(ip), pointer, intent(inout) :: Tcav_remeshed(:,:)
  !
  integer(ip), pointer :: boundary_faces(:,:)
  integer(ip) :: newCavElems
  integer(ip), pointer :: Tcav_remeshed_big(:,:)
!   integer(ip) :: iface, inode, newCavElems
!   logical(lg), pointer :: isConsideredFace(:)

  call get_boundary_faces(Tcav,boundary_faces)
  
  nullify(Tcav_remeshed_big)
  call memory_alloca(memor_adapt,'Tcav_remeshed_big','cavity_insert_point',Tcav_remeshed_big,int(size(Tcav,1),ip),int(size(boundary_faces,2),ip))
  
  call cavity_reinsert_fromCavBoun(pointToReinsert,boundary_faces,Tcav_remeshed_big,newCavElems)
  
  nullify(Tcav_remeshed)
  call memory_alloca(memor_adapt,'Tcav_remeshed','mod_cavity',Tcav_remeshed,int(size(Tcav,1),ip),newCavElems)
  Tcav_remeshed = Tcav_remeshed_big(:,1:newCavElems)
  
  call memory_deallo(memor_adapt,'Tcav_remeshed_big'  ,'cavity_insert_point',Tcav_remeshed_big)
  call delete_boundary_faces(boundary_faces)

!   nullify(isConsideredFace)
!   call memory_alloca(memor_adapt,'isConsideredFace','cavity_insert_point',isConsideredFace,int(size(boundary_faces,2),ip))
!
!   newCavElems = 0
!   do iface =1,size(boundary_faces,2)
!     isConsideredFace(iface) = .true.
!     loop_face_nodes: do inode = 1,size(boundary_faces,1)
!       if( boundary_faces(inode,iface) == pointToReinsert ) then
!         isConsideredFace(iface) = .false.
!         exit loop_face_nodes
!       end if
!     end do loop_face_nodes
!     if( isConsideredFace(iface) ) then
!       newCavElems = newCavElems + 1_ip
!     end if
!   end do
!
!   nullify(Tcav_remeshed)
!   call memory_alloca(memor_adapt,'Tcav_remeshed','cavity_insert_point',Tcav_remeshed,int(size(Tcav,1),ip),newCavElems)
!
!   newCavElems = 0
!   Tcav_remeshed(size(Tcav,1),:) = pointToReinsert
!   do iface =1,size(boundary_faces,2)
!     if( isConsideredFace(iface) ) then
!       newCavElems = newCavElems +1
!       Tcav_remeshed(1:size(boundary_faces,1),newCavElems) = boundary_faces(:,iface)
!     end if
!   end do
! !   Tcav_remeshed(size(Tcav,1),:) = Tcav_remeshed(size(Tcav,1)-1,:)
! !   Tcav_remeshed(size(Tcav,1),:) = pointToReinsert
!
!   call memory_deallo(memor_adapt,'isConsideredFace','cavity_insert_point',isConsideredFace)
!   call memory_deallo(memor_adapt,'boundary_faces'  ,'cavity_insert_point',boundary_faces)
  
  return 
end subroutine cavity_reinsert
!
!
!
subroutine cavity_reinsert_fromCavBoun(pointToReinsert,boundary_faces,Tcav_remeshed,newCavElems)
  implicit none
  
  integer(ip), intent(in) :: pointToReinsert
  integer(ip), intent(in) :: boundary_faces(:,:)
  !
  integer(ip), pointer, intent(inout) :: Tcav_remeshed(:,:)
  integer(ip), intent(out) :: newCavElems
  !
  integer(ip) :: iface, inode
  logical(lg), pointer :: isConsideredFace(:)
  integer(ip) :: numNodElem
  
  integer(ip) :: posNodeFaces(size(boundary_faces,2))
  !
! !   nullify(isConsideredFace)
! !   call memory_alloca(memor_adapt,'isConsideredFace','cavity_insert_point',isConsideredFace,int(size(boundary_faces,2),ip))
!
! !   newCavElems = 0
! !   do iface =1,size(boundary_faces,2)
! !     isConsideredFace(iface) = .true.
! !     loop_face_nodes: do inode = 1,size(boundary_faces,1)
! !       if( boundary_faces(inode,iface) == pointToReinsert ) then
! !         isConsideredFace(iface) = .false.
! !         exit loop_face_nodes
! !       end if
! !     end do loop_face_nodes
! !     if( isConsideredFace(iface) ) then
! !       newCavElems = newCavElems + 1_ip
! !     end if
! !   end do
!   posNodeFaces=FINDLOC(boundary_faces, VALUE=pointToReinsert, DIM=1)
!   newCavElems = COUNT(posNodeFaces==0_ip)
!
!   numNodElem = size(boundary_faces,1)+1
!   nullify(Tcav_remeshed)
!   call memory_alloca(memor_adapt,'Tcav_remeshed','cavity_insert_point',Tcav_remeshed,numNodElem,newCavElems)
!
!   newCavElems = 0
! !   Tcav_remeshed(numNodElem,:) = pointToReinsert
!   do iface =1,size(boundary_faces,2)
! !     if( isConsideredFace(iface) ) then
! !       newCavElems = newCavElems +1
! !       Tcav_remeshed(1:size(boundary_faces,1),newCavElems) = boundary_faces(:,iface)
! !     end if
!     if(posNodeFaces(iface)==0_ip) then
!       newCavElems = newCavElems +1
!       Tcav_remeshed(1:size(boundary_faces,1),newCavElems) = boundary_faces(:,iface)
!       Tcav_remeshed(numNodElem              ,newCavElems) = pointToReinsert
!     end if
!   end do
!
! !   call memory_deallo(memor_adapt,'isConsideredFace','cavity_insert_point',isConsideredFace)

  newCavElems = 0
  numNodElem = size(boundary_faces,1)+1
  do iface =1,size(boundary_faces,2)
    if( .not.ANY(boundary_faces(:,iface)==pointToReinsert) ) then
      newCavElems = newCavElems +1 
      Tcav_remeshed(1:size(boundary_faces,1),newCavElems) = boundary_faces(:,iface)
      Tcav_remeshed(numNodElem              ,newCavElems) = pointToReinsert
    end if
  end do
  
  return 
end subroutine cavity_reinsert_fromCavBoun
!
!
!
! subroutine computeBounFacesNotConatiningNodes(listNodes,boundary_faces,numFaces_perNode,boundary_faces_perNode)
!
!   implicit none
!
!   integer(ip), intent(in)  :: listNodes(:)
!   integer(ip), intent(in)  :: boundary_faces(:,:)
!   !
!   integer(ip), pointer, intent(out) :: numFaces_perNode(:)
!   integer(ip), pointer, intent(out) :: boundary_faces_perNode(:,:,:) !numNodFace,numFacesNotContNod,numNod
!   !
!   integer(ip) :: iface, inode, newCavElems
!   logical(lg), pointer :: isConsideredFace(:)
!   integer(ip) :: numNodElem
!   !
!   nullify(isConsideredFace)
!   call memory_alloca(memor_adapt,'isConsideredFace','cavity_insert_point',isConsideredFace,int(size(boundary_faces,2),ip))
!
!   newCavElems = 0
!   do iface =1,size(boundary_faces,2)
!     isConsideredFace(iface) = .true.
!     loop_face_nodes: do inode = 1,size(boundary_faces,1)
!       if( boundary_faces(inode,iface) == pointToReinsert ) then
!         isConsideredFace(iface) = .false.
!         exit loop_face_nodes
!       end if
!     end do loop_face_nodes
!     if( isConsideredFace(iface) ) then
!       newCavElems = newCavElems + 1_ip
!     end if
!   end do
!
!   call memory_deallo(memor_adapt,'isConsideredFace','cavity_insert_point',isConsideredFace)
!
! end subroutine computeBounFacesNotConatiningNodes
!
!
!
END MODULE mod_cavity
!> @}








