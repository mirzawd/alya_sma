!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Adaptivity
!> @{
!> @file    _adapt_swap.f90
!> @author  abel.gargallo
!> @date    2021-03-31
!> @brief   _adapt_swap
!> @details _adapt_swap: functions related with edge swap repair operations
!>
!>          The only required to use from outside is swap_mesh_edges
!>
!>
!-----------------------------------------------------------------------
!
!
!
! subroutine swap_mesh_edges_deprecated(mesh,metric,tol_qual)
!   use mod_meshTopology,       only: getIsBoundaryNode
!   implicit none
!   type(mesh_type_basic),  intent(inout)  :: mesh
!   type(mesh_metric_type), intent(in)  :: metric
!   real(rp),               intent(in)  :: tol_qual
!
!   logical(lg) :: isModifiedMesh_frozen
!   integer(ip) :: iter
!   integer(ip), parameter :: max_iter_frozen = 1000
!   integer(ip), parameter :: max_iter_global = 100
!
!   logical(lg), pointer   :: isBoundaryNode(:)
!
!   call getIsBoundaryNode(mesh,isBoundaryNode)
!
!   isModifiedMesh_frozen = .true.
!   iter = 0
!   do while ( isModifiedMesh_frozen.and.(iter<max_iter_frozen) )
!
!     iter = iter+1
!     call swap_mesh_edges_frozen(mesh,metric,isBoundaryNode,tol_qual,isModifiedMesh_frozen)
!
!   end do
!
!   call memory_deallo(memor_adapt,'isBoundaryNode',     'repair_mesh_edgeLengths',isBoundaryNode     )
!
!   return
! end subroutine swap_mesh_edges_deprecated
! !
! !
! !
! subroutine swap_mesh_edges_frozen_deprecated(mesh,metric,isBoundaryNode,tol_qual,isModifiedMesh)
!   use mod_meshTopology, only: node_to_elems_type, getNodesCav
!   use mod_quality,      only: compute_mesh_quality_sizeShape, compute_minQuality_sizeShape_cavity
!   use mod_cavity,       only: cavity_reinsert
!   implicit none
!   type(mesh_type_basic),  intent(inout) :: mesh
!   type(mesh_metric_type), intent(in)    :: metric
!   real(rp),               intent(in)    :: tol_qual
!   logical(lg), pointer,   intent(in)    :: isBoundaryNode(:)
!   logical(lg),            intent(out)   :: isModifiedMesh
!
!   type(node_to_elems_type) :: node_to_elems
!   integer(ip), pointer  :: perm_elems_byQ(:)
!   integer(ip) :: ielem, theElem, inode, jnode, n1, n2
!   logical(lg) :: isFrozenElem(mesh%nelem)
! !   logical(lg) :: isPerformedNode(mesh%npoin)
!   integer(ip), pointer :: elems_cavity(:)
!   real(rp), pointer :: q(:)
!   integer(ip) :: numNodCav
!   integer(ip),pointer :: nodesCav(:)
! !   real(rp) :: q_prev
!
!   real(rp) :: q_bestSwap, q_current
!   integer(ip), pointer :: T_bestSwap(:,:)
!   integer(ip), pointer :: cav_bestSwap(:)
!   integer(ip), pointer :: Tcav_remeshed(:,:)
!   integer(ip) :: knode, pointToReinsert
!   logical(lg) :: isImproved
!
!   logical(lg) :: isBoundEdge
!
!   integer(ip), pointer:: lnods_new(:,:)
!   integer(ip)         :: max_newElems, count_newElems, count_newPoints
!   real(rp), pointer   :: coords_new_fake(:,:)
!   !
!   call node_to_elems%set(mesh)
!
!   ! compute low quality elems
!   ! for each elem
!   !   for each edge
!   !     compute cavity of the edge
!   !     for each boundary node of the cavity
!   !       reinsert the node
!   !       compute quality of the new cavity
!   !       save something like (_q,_edge,_node)
!   !   Among all the prior combinations, perform the highest quality
!   !   that is perform the reinsertion of the _edge and cavity _node with highest _q  (_q,_edge,_node)
!
!   call compute_mesh_quality_sizeShape(mesh,metric,q)
!
!   call sort_perm_reals_increasing(q,perm_elems_byQ)
!
!   max_newElems = ceiling(0.2_rp*mesh%nelem)
!   nullify(lnods_new)
!   call memory_alloca(memor_adapt,'lnods_new','repair_mesh_edgeLengths_frozen',lnods_new,size(mesh%lnods,1),max_newElems)
!
!   isModifiedMesh = .false.
!   isFrozenElem   = .false.
!   count_newElems = 0_ip
!   loopElems: do ielem=1,mesh%nelem
!     theElem = perm_elems_byQ(ielem)
!     print*,ielem," ",q(theElem)
!     if( q(theElem)<tol_qual ) then
!       if( .not.isFrozenElem(theElem) ) then
!         !print*,"   q(theElem)",q(theElem)
!         !q_bestSwap = 0.0_rp
!         q_bestSwap = q(theElem)
!         isImproved = .false.
!         !!! ----  PERFORM SWAP ---- !!!!
!         do inode=1,mesh%mnode
!           do jnode=(inode+1_ip),mesh%mnode
!
!             n1 = mesh%lnods(inode,theElem)
!             n2 = mesh%lnods(jnode,theElem)
!
!             isBoundEdge = isBoundaryNode(n1).and.isBoundaryNode(n2)
!
!             if(.not.isBoundEdge) then
!               call node_to_elems%get_edgeContainers(n1,n2,elems_cavity)
!
!               if( .not.(any(isFrozenElem(elems_cavity))) ) then
!
!                 !!! ---- reinsert all nodes of cavity ---- !!!!
!                 call getNodesCav(mesh%lnods(:,elems_cavity),numNodCav,nodesCav)
!
!                 do knode =1,numNodCav
!                   pointToReinsert = nodesCav(knode)
!                   call cavity_reinsert(pointToReinsert,mesh%lnods(:,elems_cavity),Tcav_remeshed)
!                   q_current = compute_minQuality_sizeShape_cavity(mesh%coord,Tcav_remeshed,metric)
!                   !print*,"      q_current",q_current, "      ---> q_bestSwap: ",q_bestSwap
!
!                   if( q_current > q_bestSwap) then
!                     if( isImproved ) then
!                       call memory_deallo(memor_adapt,'T_bestSwap'  ,'swap_mesh_edges_frozen',T_bestSwap)
!                       call memory_deallo(memor_adapt,'cav_bestSwap','swap_mesh_edges_frozen',cav_bestSwap)
!                     end if
!
!                     isImproved = .true.
!                     q_bestSwap = q_current
!
!                     nullify(T_bestSwap)
!                     call memory_alloca(memor_adapt,'T_bestSwap','swap_mesh_edges_frozen',T_bestSwap,size(Tcav_remeshed,1),size(Tcav_remeshed,2))
!                     T_bestSwap = Tcav_remeshed(:,:)
!
!                     nullify(cav_bestSwap)
!                     call memory_alloca(memor_adapt,'cav_bestSwap','swap_mesh_edges_frozen',cav_bestSwap,size(elems_cavity))
!                     cav_bestSwap = elems_cavity(:)
!
!                   end if
!
!                   call memory_deallo(memor_adapt,'Tcav_remeshed','swap_mesh_edges_frozen',Tcav_remeshed)
!                 end do
!
!                 call memory_deallo(memor_adapt,'nodesCav','swap_mesh_edges_frozen',nodesCav)
!               end if !isFrozen
!               call memory_deallo(memor_adapt,'elems_cavity','swap_mesh_edges_frozen',elems_cavity)
!             end if !isBoundEdge
!           end do
!         end do
!
!         if( isImproved ) then
!           print*,'IS IMPROVEDDDDDDDDDDDD'
!           isModifiedMesh = .true.
!           isFrozenElem(cav_bestSwap) = .true.
!
!           if( (count_newElems+size(T_bestSwap,2)) > max_newElems ) then
!             go to 6666
!           end if
!
!           lnods_new(:,(count_newElems+1):(count_newElems+size(T_bestSwap,2))) = T_bestSwap
!           count_newElems = count_newElems+size(T_bestSwap,2)
!
!           call memory_deallo(memor_adapt,'T_bestSwap'  ,'swap_mesh_edges_frozen',  T_bestSwap)
!           call memory_deallo(memor_adapt,'cav_bestSwap','swap_mesh_edges_frozen',cav_bestSwap)
!         end if
!
!       end if
!     else
!       exit loopElems
!     end if
!
!   end do loopElems
!
!   6666 continue
!
!   if( isModifiedMesh ) then
!     count_newPoints = 0_ip
!     nullify(coords_new_fake)
!     call memory_alloca(memor_adapt,'coords_new_fake','swap_mesh_edges_frozen',coords_new_fake,mesh%ndime,count_newPoints)
!     !mesh = set_mesh_new(mesh,coords_new_fake,lnods_new,count_newPoints,count_newElems,isFrozenElem)
!     call set_mesh_new_modify(mesh,coords_new_fake,lnods_new,count_newPoints,count_newElems,isFrozenElem)
!   else
!     call memory_deallo(memor_adapt,'lnods_new'      ,'set_tagged_edges',lnods_new)
!   end if
!
!   call memory_deallo(memor_adapt,'q','swap_mesh_edges_frozen',q)
!
!   return
! end subroutine swap_mesh_edges_frozen_deprecated
! !
! !
! !
