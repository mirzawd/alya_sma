!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Adaptivity
!> @{
!> @file    mod_meshEdges.f90
!> @author  abel.gargallo
!> @date    2021-04-16
!> @brief   mod_meshEdges --> simplicies
!> @details mod_meshEdges 
!>
!>          To add further details
!>
!>
!-----------------------------------------------------------------------

MODULE mod_meshEdges
!***************************************************************
!*
!*  Module for metric computation and operations
!*
!***************************************************************
use def_kintyp_basic,      only: ip,rp,lg
use mod_metric,            only: mesh_metric_type
use def_kintyp_mesh_basic, only: mesh_type_basic

use def_adapt,              only : memor_adapt ! ask and replace by my own memory counter?
use mod_memory,             only : memory_alloca, memory_deallo

implicit none


private

public :: compute_mesh_edgeLengths, compute_mesh_edgeQuality
!  
!
!
CONTAINS
  !
  !
  !
  !-----------------------------------------------------------------------
  !
  !> @author  Abel Gargallo
  !> @date    18/04/2021
  !> @brief   compute_mesh_edgeLengths (mesh_type_basic, simplices)
  !> @details compute_mesh_edgeLengths
  !
  !-----------------------------------------------------------------------
subroutine compute_mesh_edgeLengths(mesh,metric,edge_to_node,edge_lengths)
  use mod_meshTopology, only: compute_mesh_edges
  use mod_debugTools,   only: out_performance, deb_meshEdgeLen,  deb_meshEdges, deb_allEdgesLen
  implicit none
  !
  type(mesh_type_basic),    intent(in) :: mesh
  type(mesh_metric_type),   intent(in) :: metric
  real(rp),    pointer,  intent(inout) :: edge_lengths(:)
  integer(ip), pointer,  intent(inout) :: edge_to_node(:,:)
  !
  integer(ip) :: iedge,n1,n2
  real(rp) :: t0,t1
  
  real(rp) :: l1(1,1), l2(1,1), x1x2(mesh%ndime,1), l, hmetric
  real(rp) :: len1,len2
  !
  if(mesh%npoin==0) return
  
  if(out_performance) call cpu_time(t0)
  
  call compute_mesh_edges(mesh,edge_to_node)
  
  if(out_performance) call cpu_time(t1)
  if(out_performance) deb_meshEdges = deb_meshEdges + (t1-t0)
  if(out_performance) call cpu_time(t0)
  
  nullify(edge_lengths)
  call memory_alloca(memor_adapt,'EDGE_LENGTHS','compute_mesh_edgeLengths',edge_lengths, int(size(edge_to_node,2),ip) )
  
  if(metric%isSizeField) then
    do iedge=1,size(edge_to_node,2)
      n1 = edge_to_node(1,iedge)
      n2 = edge_to_node(2,iedge)
    
      l = my_norm2( mesh%coord(:,n1)-mesh%coord(:,n2) )
      hmetric = (1.0_rp/metric%get_size_node(n1) + 1.0_rp/metric%get_size_node(n2))/2.0
      edge_lengths(iedge)  = l*hmetric
    end do
  else
    do iedge=1,size(edge_to_node,2)
      n1 = edge_to_node(1,iedge)
      n2 = edge_to_node(2,iedge)
    
      !x1x2(:,1) = mesh%coord(:,n1)-mesh%coord(:,n2)
      !l1 = sqrt( matmul( transpose(x1x2) , matmul( metric%get_metric_node(n1) , x1x2 )  )  )
      !l2 = sqrt( matmul( transpose(x1x2) , matmul( metric%get_metric_node(n2) , x1x2 )  )  )
      !edge_lengths(iedge)  = ( l1(1,1) + l2(1,1) )/2.0_rp
      edge_lengths(iedge) = compute_edge_length(&
       mesh%coord(:,n1),            mesh%coord(:,n2),&
       metric%get_metric_node(n1),  metric%get_metric_node(n2) )
    end do
  end if
  
  if(out_performance) call cpu_time(t1)
  if(out_performance) deb_allEdgesLen = deb_allEdgesLen + (t1-t0)
  if(out_performance) deb_meshEdgeLen = deb_allEdgesLen + deb_meshEdges
  
end subroutine compute_mesh_edgeLengths  
  !-----------------------------------------------------------------------
  !
  !> @author  Abel Gargallo
  !> @date    18/04/2021
  !> @brief   compute_mesh_edgeLengths (mesh_type_basic, simplices)
  !> @details compute_mesh_edgeLengths
  !
  !-----------------------------------------------------------------------
subroutine compute_mesh_edgeQuality(mesh,node_to_elems,q,edge_to_node,edge_q)
  use mod_meshTopology, only: compute_mesh_edges
  use mod_meshTopology, only: node_to_elems_type
  implicit none
  !
  type(mesh_type_basic),    intent(in) :: mesh
  type(node_to_elems_type), intent(in) :: node_to_elems
  real(rp), pointer    ,    intent(in) :: q(:)
  !
  real(rp),    pointer, intent(inout)  :: edge_q(:)    
  integer(ip), pointer, intent(inout)  :: edge_to_node(:,:)
  !
  integer(ip) :: iedge,n1,n2
  integer(ip),pointer :: edgeElems(:)
  !
  call compute_mesh_edges(mesh,edge_to_node)

  nullify(edge_q)
  call memory_alloca(memor_adapt,'edge_q','compute_mesh_edgeQuality',edge_q, int(size(edge_to_node,2),ip) )
  do iedge=1,size(edge_to_node,2)
    n1 = edge_to_node(1,iedge)
    n2 = edge_to_node(2,iedge)
    
    call node_to_elems%get_edgeContainers(n1,n2,edgeElems)
    
    edge_q(iedge) = compute_edgeQ_fromElemQ( q(edgeElems) )
    
    call node_to_elems%delete_elems(edgeElems)
  end do
  
end subroutine compute_mesh_edgeQuality
!
!
!
function compute_edgeQ_fromElemQ(q_edgeElems) result(q_edge)
  implicit none
  
  real(rp), intent(in) :: q_edgeElems(:)
  real(rp) :: q_edge
  real(rp), parameter :: p_norm = 4.0
  !
  q_edge = minval( q_edgeElems )!compute_minQuality_sizeshape_subset(mesh,metric,list_elems)
  !q_edge = 1.0/( sum(1.0/(q_edgeElems)**p_norm)/size(q_edgeElems) )**(1.0/p_norm) ! NORM p_norm, normalized by the num of elems
  !
end function compute_edgeQ_fromElemQ

!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    18/04/2021
!> @brief   compute_mesh_edgeLengths (mesh_type_basic, simplices)
!> @details compute_mesh_edgeLengths
!
!-----------------------------------------------------------------------
function compute_edge_length(x1,x2,M1,M2) result(l)
  use mod_metric, only: interpolation_law, law_geometric_h
  implicit none
  
  real(rp),    intent(in)  :: x1(:),x2(:),M1(:,:),M2(:,:)
  real(rp) :: l
  
  real(rp) :: l1(1,1), l2(1,1), x1x2(size(x1,1),1)
  real(rp) :: a,len1,len2

!   x1x2(:,1) = x1-x2
!
!   l1 = sqrt( matmul( transpose(x1x2) , matmul( M1 , x1x2 )  )  )
!   l2 = sqrt( matmul( transpose(x1x2) , matmul( M2 , x1x2 )  )  )
!
!   l = ( l1(1,1) + l2(1,1) )/2.0_rp
  
  !!!! The following is the mathematically correct.. but very expensive...
  ! From Alauzet, Size gradation
  l1 = sqrt(  matmul( transpose(x1x2) , matmul( M1 , x1x2 )  )  )
  l2 = sqrt(  matmul( transpose(x1x2) , matmul( M2 , x1x2 )  )  )

  if( abs(l1(1,1)-l2(1,1))/l1(1,1) < 1.0e-12_rp) then
    l=l1(1,1)
    return
  elseif( l1(1,1)>l2(1,1) ) then
    len1 = l1(1,1)
    len2 = l2(1,1)
  else
    len2 = l1(1,1)
    len1 = l2(1,1)
  end if
  
  a = len1/len2
  
  !!! Linear on h:
  !l = len1*log(a)/(a-1.0_rp)
  !!! Geometric on h
  if(interpolation_law.eq.law_geometric_h) then
    l = len1*(a-1.0_rp)/(a*log(a))
  else
    call runend('compute_edge_length: defaulted geometric law on h')
  end if
  !!! Linear interp on h^-1
  !l = 0.5_rp * len1*(a+1.0_rp)/a
  !!! Linear interp on eigen
  !l = len1*(2.0_rp/3.0_rp)*(a*a + a + 1.0_rp)/(a*(a+1.0_rp))
  
end function compute_edge_length
!
!
!
function my_norm2(v) result(v_norm2)
  implicit none
  !
  real(rp), intent(in) :: v(:)
  real(rp) :: v_norm2
  
  v_norm2 = sqrt(sum(v**2_ip))
  !
  return 
end function my_norm2
!
!
!
END MODULE mod_meshEdges

!> @}




