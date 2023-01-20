!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Adaptivity
!> @{
!> @file    mod_meshTopology.f90
!> @author  abel.gargallo
!> @date    2021-04-16
!> @brief   mod_meshTopology --> simplicies
!> @details mod_meshTopology 
!>
!>          To add further details
!>
!>
!-----------------------------------------------------------------------

MODULE mod_meshTopology
!***************************************************************
!*
!*  Module for metric computation and operations
!*
!***************************************************************
use def_kintyp_basic,      only: ip,rp,lg
use def_kintyp_mesh_basic, only: mesh_type_basic
use mod_elmgeo,             only : element_type

use def_adapt,              only : memor_adapt ! ask and replace by my own memory counter?
use mod_memory,             only : memory_alloca, memory_deallo

use mod_edict, only: edict_type

implicit none

private

public :: compute_mesh_edges
public :: node_to_elems_type
public :: get_boundary_faces, delete_boundary_faces
public :: getIsBoundaryNode
public :: setIsBoundaryNodeFromIniMesh
public :: deleteIsBouNode
public :: getNodesCav, deleteNodesCav
public :: node_to_nodes_type
public :: delete_arrayElems

public :: edgeData_type
public :: setBoundaryEdgeData

type node_to_elems_type
  integer(ip), pointer :: n_to_e(:)
  integer(ip), pointer :: n_to_fisrt_e(:) ! size n+1
  !integer(ip), pointer :: n_to_last_e(:)
 contains
   procedure,   pass    :: set                => compute_node_to_elems
   procedure,   pass    :: get_elems          => get_elems_adj_to_node   ! get an element from the list
   procedure,   pass    :: get_edgeContainers => get_elems_containingEdge! get intersection of two
   procedure,   pass    :: get_edgeAdjacency  => get_elems_adjacentToEdge! get union of two
   procedure,   pass    :: deallo             => delete_node_to_elems
   procedure,   pass    :: delete_elems       => delete_elems_forType ! to delete elemes tracking memory
end type node_to_elems_type 

type node_to_nodes_type
  integer(ip), pointer :: n_to_n(:)
  integer(ip), pointer :: n_to_fisrt_n(:) ! size n+1
 contains
   procedure,   pass    :: set                => compute_node_to_nodes
   procedure,   pass    :: get_adj_nodes      => get_nodes_adj_to_node
   procedure,   pass    :: deallo             => delete_node_to_nodes
end type node_to_nodes_type 

type edgeData_type
    type(edict_type) :: edge_dict
  contains
    procedure, pass :: set      => setEdgeData
    procedure, pass :: setEdge
    procedure, pass :: isEdge
    procedure, pass :: alloca   => allocate_edgeData
    procedure, pass :: deallo   => deallo_edgeData
    procedure, pass :: isAlloca => isAllocated_edgeData
end type edgeData_type 
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
!> @date    28/04/2021
!> @brief   setBoundaryEdgeData
!> @details setBoundaryEdgeData: set boundary mesh edge
!
!-----------------------------------------------------------------------
function setBoundaryEdgeData(lnods) result(edgeData)
  implicit none
  integer(ip), intent(in) :: lnods(:,:)
  type(edgeData_type)     :: edgeData
  !
  integer(ip), pointer :: lnods_bound(:,:)
  !
  call get_boundary_faces(lnods,lnods_bound) 
  call edgeData%set(lnods_bound)
  call delete_boundary_faces(lnods_bound)
  !
end function setBoundaryEdgeData
!
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   setEdgeData
!> @details setEdgeData: set the edges of a mesh (using lnods, Topology)
!
!-----------------------------------------------------------------------
subroutine setEdgeData(edgeData,T)
  implicit none
  class(edgeData_type), intent(out) :: edgeData
  integer(ip),          intent(in)  :: T(:,:)
  !
  integer(ip) :: ielem, i, j
  integer(ip) :: sizeDict
  !
  sizeDict = size(T,2)*size(T,1)
  call edgeData%edge_dict%init(sizeDict)
  
  do ielem=1,size(T,2)
    do i=1,size(T,1)
      do j=i+1,size(T,1)
        call edgeData%edge_dict%setByNodes(T(i,ielem),T(j,ielem),.true.)
      end do
    end do
  end do
  !
end subroutine setEdgeData
!
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   setEdge
!> @details setEdge: include an edge in the dictionary
!
!-----------------------------------------------------------------------
subroutine setEdge(edgeData,n1,n2)
  implicit none
  class(edgeData_type), intent(inout) :: edgeData
  integer(ip),          intent(in)    :: n1,n2
  !
  call edgeData%edge_dict%setByNodes(n1,n2,.true.)
  !
end subroutine setEdge
!
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   alloca
!> @details allocate size and initialize edgeData
!
!-----------------------------------------------------------------------
subroutine allocate_edgeData(edgeData,sizeDict)
  implicit none
  class(edgeData_type), intent(out) :: edgeData
  integer(ip),          intent(in)  :: sizeDict
  !
  call edgeData%edge_dict%init(sizeDict)
  !
end subroutine allocate_edgeData
!
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   alloca
!> @details allocate size and initialize edgeData
!
!-----------------------------------------------------------------------
subroutine deallo_edgeData(edgeData)
  implicit none
  class(edgeData_type), intent(inout) ::edgeData
  !
  call edgeData%edge_dict%deallo()
  !
end subroutine deallo_edgeData
!
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   alloca
!> @details allocate size and initialize edgeData
!
!-----------------------------------------------------------------------
function isAllocated_edgeData(edgeData) result(isAllocated)
  implicit none
  class(edgeData_type), intent(in) :: edgeData
  logical(lg)                      :: isAllocated
  !
  isAllocated = edgeData%edge_dict%isAlloca()
  !
end function isAllocated_edgeData
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   isEdge
!> @details isEdge
!
!-----------------------------------------------------------------------
function isEdge(edgeData,n1,n2) result(isTrueEdge)
  implicit none
  class(edgeData_type), intent(in) :: edgeData
  integer(ip),          intent(in) :: n1,n2
  logical(lg) :: isTrueEdge
  !
  isTrueEdge = edgeData%edge_dict%areEdge(n1,n2)
  !
end function isEdge
!
!
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   get Nodes Cavity
!> @details get Nodes Cavity
!
!-----------------------------------------------------------------------
subroutine getNodesCav(Tcav,numNodCav,nodesCav)
  use mod_maths_sort,   only: maths_heap_sort
  implicit none
  integer(ip),          intent(in)    :: Tcav(:,:)
  integer(ip),          intent(out)   :: numNodCav
  integer(ip), pointer, intent(inout) :: nodesCav(:)
  
  numNodCav = size(Tcav,1)*size(Tcav,2)!elems_cavity)
  nullify(nodesCav)
  call memory_alloca(memor_adapt,'nodesCav','mod_meshTopology',nodesCav,numNodCav)
  nodesCav = reshape( Tcav ,(/numNodCav/) )
  
  call maths_heap_sort(itask=2_ip, nrows=numNodCav,ivin=nodesCav,ELIMINATE_REPLICATES=.true.) !ivo1=sort_perm,

!   integer(ip), intent(in) ::Tcav(:,:)
!   integer(ip), intent(out) :: numNodCav
!   integer(ip), pointer, intent(inout) :: nodesCav(:)
!
!   integer(ip) :: numNodCav_big, inode
!   integer(ip),pointer :: nodesCav_big(:)
!
!   numNodCav_big = size(Tcav,1)*size(Tcav,2)!elems_cavity)
!   nullify(nodesCav_big)
!   call memory_alloca(memor_adapt,'nodesCav_big','mod_meshTopology',nodesCav_big,numNodCav_big)
!   nodesCav_big = reshape( Tcav ,(/numNodCav_big/) )
!
!   call maths_heap_sort(itask=2_ip, nrows=numNodCav_big,ivin=nodesCav_big,ELIMINATE_REPLICATES=.true.) !ivo1=sort_perm,
!
!   numNodCav =numNodCav_big
!   !inode = 1_ip
!   !do while((nodesCav_big(inode)>0_ip).and.(inode<=numNodCav_big))
!   !  numNodCav = inode!numNodCav+1_ip
!   !  inode = inode+1_ip
!   !end do
!
!   nullify(nodesCav)
!   call memory_alloca(memor_adapt,'nodesCav','mod_meshTopology',nodesCav,numNodCav)
!   nodesCav(1:numNodCav) = nodesCav_big(1:numNodCav)
!
!   call memory_deallo(memor_adapt,'nodesCav_big','mod_meshTopology',nodesCav_big)
  
end subroutine getNodesCav
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   get Nodes Cavity
!> @details get Nodes Cavity
!
!-----------------------------------------------------------------------
subroutine deleteNodesCav(nodesCav)
  implicit none
  integer(ip), pointer, intent(inout) :: nodesCav(:)
  !
  call memory_deallo(memor_adapt,'nodesCav','mod_meshTopology',nodesCav)
  !
end subroutine deleteNodesCav
!
!
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   delete_node_to_elems data structure
!> @details delete_node_to_elems (deallocate and nullify pointers)
!
!-----------------------------------------------------------------------
subroutine delete_node_to_elems(node_to_elems)!node_to_elems,node_to_firstElem,node_to_lastElem)
  implicit none

  class(node_to_elems_type),       intent(inout)   :: node_to_elems

  call memory_deallo(memor_adapt,'n_to_e','delete_node_to_elems',node_to_elems%n_to_e)
  call memory_deallo(memor_adapt,'n_to_fisrt_e','delete_node_to_elems',node_to_elems%n_to_fisrt_e)
  
  nullify(node_to_elems%n_to_e)
  nullify(node_to_elems%n_to_fisrt_e)
  !
end subroutine delete_node_to_elems
!
!
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   compute_node_to_elems
!> @details compute_node_to_elems
!
!-----------------------------------------------------------------------
subroutine compute_node_to_elems(node_to_elems,meshe)!node_to_elems,node_to_firstElem,node_to_lastElem)
  use mod_debugTools, only: out_performance,deb_nodeToElem
  implicit none

  class(node_to_elems_type),      intent(inout) :: node_to_elems
  type(mesh_type_basic),          intent(in)    :: meshe
  !
  integer(ip)                          :: mpop2,ipoin
  integer(ip)                          :: ielem
  integer(ip)                          :: inode
  integer(ip)                          :: nlelp

  integer(ip),                 pointer :: nepoi(:)
  integer(ip),                 pointer :: pelpo(:)
  integer(ip),                 pointer :: lelpo(:)
  integer(ip),                 pointer :: ledgp(:)
  integer(ip),                 pointer :: pedgp(:)

  integer(ip),                 pointer :: lnods(:,:)
  integer(ip),                 pointer :: ltype(:)
  integer(ip)                          :: nelem
  integer(ip)                          :: npoin              
  
  integer(ip)                          :: lnnod_simplex
  
  real(rp) :: t0,t1
  !
  if(out_performance) call cpu_time(t0)
  !
  nullify(nepoi)
  nullify(pelpo)
  nullify(lelpo)
  nullify(ledgp)
  nullify(pedgp)

  nullify(lnods)
  nullify(ltype)

  lnods => meshe % lnods
  ltype => meshe % ltype
  nelem =  meshe % nelem
  npoin =  meshe % npoin
  
  lnnod_simplex = element_type(meshe%ltype(1))%number_nodes
  !
  ! Allocate memory for NEPOI and compute it
  ! 
  call memory_alloca(memor_adapt,'NEPOI','compute_node_to_elems',nepoi,npoin)
  mpop2 = 0
  do ielem = 1,nelem
     mpop2 = mpop2 + lnnod_simplex**2!lnnod(ielem)*lnnod(ielem)
     do inode = 1,lnnod_simplex!lnnod(ielem)
        ipoin = lnods(inode,ielem)
        nepoi(ipoin) = nepoi(ipoin) + 1
     end do
  end do
  !
  ! Allocate memory for PELPO and compute it
  !
  call memory_alloca(memor_adapt,'PELPO','compute_node_to_elems',pelpo,npoin+1_ip)
  pelpo(1) = 1
  do ipoin = 1,npoin
     pelpo(ipoin+1) = pelpo(ipoin) + nepoi(ipoin)
  end do
  !
  ! Allocate memory for LELPO and construct the list
  !
  nlelp = pelpo(npoin+1)!-1 ! this -1 is mine => I'll leave as it is
  call memory_alloca(memor_adapt,'LELPO','compute_node_to_elems',lelpo,nlelp)
  do ielem = 1,nelem
     do inode = 1,lnnod_simplex!lnnod(ielem)
        ipoin = lnods(inode,ielem)
        lelpo(pelpo(ipoin)) = ielem
        pelpo(ipoin) = pelpo(ipoin) + 1
     end do
  end do
  !
  ! Recompute PELPO and maximum number of element neighbors MEPOI
  !
  pelpo(1) = 1
  do ipoin = 1,npoin
     pelpo(ipoin+1) = pelpo(ipoin) + nepoi(ipoin)
  end do
  !
  ! Construct data structure
  !
  node_to_elems%n_to_e  => lelpo
  node_to_elems%n_to_fisrt_e  => pelpo
  
  call memory_deallo(memor_adapt,'NEPOI','compute_node_to_elems',nepoi)
  !
  
  if(out_performance) then
    call cpu_time(t1)
    deb_nodeToElem = deb_nodeToElem + (t1-t0)
  end if
  !
end subroutine compute_node_to_elems
!
!
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   get_elems_adj_to_node
!> @details get_elems_adj_to_node
!
!-----------------------------------------------------------------------
subroutine get_elems_adj_to_node(node_to_elems,inode,elems)
  implicit none
  
  class(node_to_elems_type), intent(in) :: node_to_elems
  integer(ip), intent(in) :: inode
  
  !integer(ip) :: elems( node_to_elems%n_to_last_e(inode) - node_to_elems%n_to_fisrt_e(inode) + 1  )
  integer(ip), pointer ,intent(inout) :: elems(:)
  
  integer(ip) :: numElems
  
  numElems =  node_to_elems%n_to_fisrt_e(inode+1) - node_to_elems%n_to_fisrt_e(inode)

  nullify(elems)
  call memory_alloca(memor_adapt,'elems','mod_meshTopology',elems,numElems)
  if(numElems>0_ip) then
    elems(:) = node_to_elems%n_to_e(&
      node_to_elems%n_to_fisrt_e(inode) : (node_to_elems%n_to_fisrt_e(inode+1)-1) )  
  end if
    
  !elems_pointer => elems
 
!   elems = node_to_elems%n_to_e(&
!           node_to_elems%n_to_fisrt_e(inode) : node_to_elems%n_to_last_e(inode) &
!           )

!   if( inode < node_to_elems%nn ) then
!      elems = node_to_elems%n_to_e(&
!        node_to_elems%n_to_fisrt_e(inode) : (node_to_elems%n_to_fisrt_e(inode+1)-1) )
!   else if ( inode == node_to_elems%nn ) then
!      elems = node_to_elems%n_to_e(&
!        node_to_elems%n_to_fisrt_e(inode):size(node_to_elems%n_to_e) )
!   else
!     call runend("inode larger than nodes in data structure")
!   end if
end subroutine get_elems_adj_to_node
!
!
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   get_elems_containingEdge
!> @details get_elems_containingEdge
!
!-----------------------------------------------------------------------
subroutine get_elems_containingEdge(node_to_elems,n1,n2,elems)
  implicit none
  
  class(node_to_elems_type),  intent(in)  :: node_to_elems
  integer(ip),                intent(in)  :: n1,n2
  integer(ip), pointer,       intent(inout) :: elems(:)
  
  integer(ip) :: numElems, ielem, jelem, theElem1, theElem2
  integer(ip) :: elems1( node_to_elems%n_to_fisrt_e(n1+1) - node_to_elems%n_to_fisrt_e(n1) )
  integer(ip) :: elems2( node_to_elems%n_to_fisrt_e(n2+1) - node_to_elems%n_to_fisrt_e(n2) )
  integer(ip) :: elems_max(&  
    node_to_elems%n_to_fisrt_e(n1+1) - node_to_elems%n_to_fisrt_e(n1) + 1 + &
    node_to_elems%n_to_fisrt_e(n2+1) - node_to_elems%n_to_fisrt_e(n2) + 1 )

  elems1(:) = node_to_elems%n_to_e(node_to_elems%n_to_fisrt_e(n1) : (node_to_elems%n_to_fisrt_e(n1+1)-1) )  
  elems2(:) = node_to_elems%n_to_e(node_to_elems%n_to_fisrt_e(n2) : (node_to_elems%n_to_fisrt_e(n2+1)-1) )  

  numElems = 0
  loopi: do ielem=1,size(elems1)
    theElem1 = elems1(ielem)
    loopj: do jelem=1,size(elems2)
      theElem2 = elems2(jelem)
      if(theElem1==theElem2) then
        numElems = numElems+1
        elems_max(numElems) = theElem1
        exit loopj
      end if
    end do loopj
  end do loopi
  
  nullify(elems)
  call memory_alloca(memor_adapt,'elems','mod_meshTopology',elems,numElems)
  elems(:) = elems_max(1:numElems)
  
end subroutine get_elems_containingEdge
!
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   get_elems_adjacentToEdge
!> @details get_elems_adjacentToEdge
!
!-----------------------------------------------------------------------
subroutine get_elems_adjacentToEdge(node_to_elems,n1,n2,elems)
  implicit none
  
  class(node_to_elems_type),  intent(in)  :: node_to_elems
  integer(ip),                intent(in)  :: n1,n2
  integer(ip), pointer,       intent(inout) :: elems(:)
  
  integer(ip) :: numElems, ielem, jelem, theElem1, theElem2
  integer(ip) :: elems1( node_to_elems%n_to_fisrt_e(n1+1) - node_to_elems%n_to_fisrt_e(n1) )
  integer(ip) :: elems2( node_to_elems%n_to_fisrt_e(n2+1) - node_to_elems%n_to_fisrt_e(n2) )
  integer(ip) :: elems_max(&  
    node_to_elems%n_to_fisrt_e(n1+1) - node_to_elems%n_to_fisrt_e(n1) + 1 + &
    node_to_elems%n_to_fisrt_e(n2+1) - node_to_elems%n_to_fisrt_e(n2) + 1 )
  logical(lg) :: isFound

  elems1(:) = node_to_elems%n_to_e(node_to_elems%n_to_fisrt_e(n1) : (node_to_elems%n_to_fisrt_e(n1+1)-1) )  
  elems2(:) = node_to_elems%n_to_e(node_to_elems%n_to_fisrt_e(n2) : (node_to_elems%n_to_fisrt_e(n2+1)-1) )  

  elems_max(1:size(elems2)) = elems2
  numElems = size(elems2)
  loopi: do ielem=1,size(elems1)
    theElem1 = elems1(ielem)
    isFound = .false.
    loopj: do jelem=1,size(elems2)
      theElem2 = elems2(jelem)
      if(theElem1==theElem2) then
        isFound = .true.
        exit loopj
      end if
    end do loopj
    if(.not.isFound) then
      numElems = numElems+1
      elems_max(numElems) = theElem1
    end if
  end do loopi
  
  nullify(elems)
  call memory_alloca(memor_adapt,'elems','mod_meshTopology',elems,numElems)

  elems(:) = elems_max(1:numElems)
  
end subroutine get_elems_adjacentToEdge
!
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   delete_arrayElems
!> @details delete_arrayElems
!
!-----------------------------------------------------------------------
subroutine delete_arrayElems(elems)
  implicit none
  integer(ip), pointer, intent(inout) :: elems(:)
  !
  call memory_deallo(memor_adapt,'elems','mod_meshTopology',elems)
  !
end subroutine delete_arrayElems
!
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   delete_elems_forType
!> @details delete_arrayElems
!
!-----------------------------------------------------------------------
subroutine delete_elems_forType(node_to_elems,elems)
  implicit none
  class(node_to_elems_type), intent(in) :: node_to_elems
  integer(ip), pointer ,intent(inout) :: elems(:)
  call memory_deallo(memor_adapt,'elems','mod_meshTopology',elems)
end subroutine delete_elems_forType
!-----------------------------------------------------------------------
!
!> @author  Guillaume Houzeaux (adapted by Abel Gargallo for simplices and mesh_type_basic)
!> @date    28/04/2021
!> @brief   Compute the list of edges of mesh type basic
!> @details Compute the list of edges
!
!-----------------------------------------------------------------------
subroutine compute_mesh_edges(meshe,edge_to_node)
  implicit none
  
  type(mesh_type_basic),      intent(in)    :: meshe
  integer(ip),      pointer , intent(inout) :: edge_to_node(:,:)         ! NEDGE

  integer(ip)  :: nedge

  integer(ip)                          :: mpop2,ipoin
  integer(ip)                          :: lsize,iedgg
  integer(ip)                          :: ielem
  integer(ip)                          :: inode,ilisn
  integer(ip)                          :: nlelp,jpoin
  integer(ip)                          :: ipoin_1,ipoin_2
  integer(ip)                          :: jj,iedge,pelty
  integer(ip)                          :: ielpo
  logical(lg)                          :: notfound

  integer(ip),                 pointer :: nepoi(:)
  integer(ip),                 pointer :: pelpo(:)
  integer(ip),                 pointer :: lelpo(:)
  integer(ip),                 pointer :: ledgp(:)
  integer(ip),                 pointer :: pedgp(:)

!   integer(ip),                 pointer :: lnnod(:)
  integer(ip),                 pointer :: lnods(:,:)
!   integer(ip),                 pointer :: lnnob(:)
!   integer(ip),                 pointer :: lnodb(:,:)
  integer(ip),                 pointer :: ltype(:)
  integer(ip)                          :: nelem
!  integer(ip)                          :: nboun
  integer(ip)                          :: npoin              
  
  integer(ip)                          :: lnnod_simplex

  if(meshe % npoin == 0) then
    return
  end if

  nullify(nepoi)
  nullify(pelpo)
  nullify(lelpo)
  nullify(ledgp)
  nullify(pedgp)

!   nullify(lnnod)
  nullify(lnods)
!   nullify(lnnob)
!   nullify(lnodb)
  nullify(ltype)

  !lnnod => meshe % lnnod
  lnods => meshe % lnods
  !lnnob => meshe % lnnob
  !lnodb => meshe % lnodb
  ltype => meshe % ltype
  nelem =  meshe % nelem
  !nboun =  meshe % nboun
  npoin =  meshe % npoin
  
  lnnod_simplex = element_type(meshe%ltype(1))%number_nodes
  !
  ! Allocate memory for NEPOI and compute it
  ! 
  call memory_alloca(memor_adapt,'NEPOI','compute_mesh_edges',nepoi,npoin)
  mpop2 = 0
  do ielem = 1,nelem
     mpop2 = mpop2 + lnnod_simplex**2!lnnod(ielem)*lnnod(ielem)
     do inode = 1,lnnod_simplex!lnnod(ielem)
        ipoin = lnods(inode,ielem)
!         if(ipoin<1_ip) then
!         !  print*,lnods
!         !  print*,lnods(inode,ielem)
!         end if
        nepoi(ipoin) = nepoi(ipoin) + 1
     end do
  end do
  !
  ! Allocate memory for PELPO and compute it
  !
  call memory_alloca(memor_adapt,'PELPO','compute_mesh_edges',pelpo,npoin+1_ip)
  pelpo(1) = 1
  do ipoin = 1,npoin
     pelpo(ipoin+1) = pelpo(ipoin) + nepoi(ipoin)
  end do
  !
  ! Allocate memory for LELPO and construct the list
  !
  nlelp = pelpo(npoin+1)
  call memory_alloca(memor_adapt,'LELPO','compute_mesh_edges',lelpo,nlelp)
  do ielem = 1,nelem
     do inode = 1,lnnod_simplex!lnnod(ielem)
        ipoin = lnods(inode,ielem)
        lelpo(pelpo(ipoin)) = ielem
        pelpo(ipoin) = pelpo(ipoin) + 1
     end do
  end do
  !
  ! Recompute PELPO and maximum number of element neighbors MEPOI
  !
  pelpo(1) = 1
  do ipoin = 1,npoin
     pelpo(ipoin+1) = pelpo(ipoin) + nepoi(ipoin)
  end do
  !
  ! Allocate memory
  !
  call memory_alloca(memor_adapt,'LEDGP','compute_mesh_edges',ledgp,mpop2)
  call memory_alloca(memor_adapt,'PEDGP','compute_mesh_edges',pedgp,npoin+1_ip)
  !
  ! Construct the array of indexes
  !     
  pedgp(1) = 1
  do ipoin = 1,npoin
     lsize = 0
     do ielpo = pelpo(ipoin),pelpo(ipoin+1)-1
        ielem = lelpo(ielpo)
        pelty = ltype(ielem)
        do iedge = 1,element_type(pelty) % number_edges
           ipoin_1 = lnods(element_type(pelty) % list_edges(1,iedge),ielem)
           ipoin_2 = lnods(element_type(pelty) % list_edges(2,iedge),ielem)
           if(      ipoin == ipoin_1 .and. ipoin_2 > ipoin ) then
              jpoin = ipoin_2
           else if( ipoin == ipoin_2 .and. ipoin_1 > ipoin ) then
              jpoin = ipoin_1
           else
              jpoin = 0
           end if
           if( jpoin /= 0 ) then
              jj = 0
              notfound = .true.
              do while( jj < lsize .and. notfound )
                 if( jpoin == ledgp(pedgp(ipoin)+jj) ) then
                    notfound = .false.
                 end if
                 jj = jj + 1
              end do
              if( notfound ) then
                 ledgp(pedgp(ipoin)+lsize) = jpoin
                 lsize = lsize + 1
              end if
           end if
        end do
     end do
     pedgp(ipoin+1) = pedgp(ipoin) + lsize
  end do
  nedge = pedgp(npoin+1) - 1       
  !
  ! EDGE_TO_NODE(1:2,:) ... Fill in edge table
  !
  nullify(edge_to_node)
  call memory_alloca(memor_adapt,' % EDGE_TO_NODE','compute_mesh_edges',edge_to_node,2_ip,nedge)
  iedgg = 0
  do ipoin = 1,npoin
     do ilisn = 1,pedgp(ipoin+1)-pedgp(ipoin)
        iedgg                         = iedgg + 1
        jpoin                         = ledgp(iedgg)
        if( ipoin < jpoin ) then
           edge_to_node(1,iedgg) = ipoin
           edge_to_node(2,iedgg) = jpoin
        else
           edge_to_node(1,iedgg) = jpoin
           edge_to_node(2,iedgg) = ipoin
        end if
     end do
  end do
  !
  ! Deallocate memory
  !
  call memory_deallo(memor_adapt,'LEDGP','compute_mesh_edges',ledgp)
  call memory_deallo(memor_adapt,'PEDGP','compute_mesh_edges',pedgp)
  call memory_deallo(memor_adapt,'LELPO','compute_mesh_edges',lelpo)
  call memory_deallo(memor_adapt,'PELPO','compute_mesh_edges',pelpo)
  call memory_deallo(memor_adapt,'NEPOI','compute_mesh_edges',nepoi)
  
  
end subroutine compute_mesh_edges
!
!
!
!-----------------------------------------------------------------------
!
!> @author  Guillaume Houzeaux (adapted by Abel Gargallo to lnods of simplices)
!> @date    28/04/2021
!> @brief   get_boundary_faces from element connectivity
!> @details get_boundary_faces: the idea is not duplicating code and using node_to_elems
!
!-----------------------------------------------------------------------
subroutine get_boundary_faces(lnods,lnods_bound)  
  use def_elmtyp,         only : BAR02,TRI03,TET04
  use mod_maths_sort,     only : maths_heap_sort

  implicit none

  integer(ip), intent(in)     :: lnods(:,:)
  integer(ip), pointer ,intent(inout)    :: lnods_bound(:,:)

  integer(ip)                           :: ielty,ielem,iface,inodf
  integer(ip)                           :: inode,jelem,jface,jelty,ipoin,pnodf
  integer(ip)                           :: ielpo,mface,mnodf,nboun
  integer(ip)                           :: kpoin
  !integer(8)                            :: memor_adapt(2)
  logical(lg)                           :: equal_faces  

  integer(ip)                           :: mepoi,nlelp
  integer(ip), pointer                  :: lelpo(:)
  integer(ip), pointer                  :: pelpo(:)
  integer(ip), pointer                  :: nepoi(:)

  integer(ip)                           :: nelem
  integer(ip)                           :: npoin
  integer(ip)                           :: mnode
  integer(ip)                           :: ltype
  integer(ip)                           :: ltype_face
  integer(ip)                           :: lnnod_simplex

  integer(ip)                           :: nfacg
  integer(ip), pointer                  :: facel(:,:,:)
  integer(ip), pointer                  :: permr(:)
  integer(ip), pointer                  :: invpr(:)

  !memor_adapt = 0_8

  ! remapejar lnods a 1:nombre de nodes per estalviar memoria..


  nelem = size(lnods,2)
  mnode = size(lnods,1)
  lnnod_simplex = size(lnods,1)
  npoin =  maxval(lnods)
  
  if (lnnod_simplex==3) then
    ltype = TRI03
    ltype_face = BAR02
  elseif(lnnod_simplex==4) then
    ltype = TET04
    ltype_face = TRI03
  else
    call runend("Not impelmented for this element type")
  end if
  

  nullify(lelpo)
  nullify(pelpo)
  nullify(nepoi)

  nullify(facel)
  nullify(permr)
  nullify(invpr)

  !----------------------------------------------------------------------
  !
  ! LELPO,PELPO: node-to-element graph
  !
  !----------------------------------------------------------------------
  
  mface = mnode    ! Simplices
  mnodf = mnode-1  ! Simplices

  call memory_alloca(memor_adapt,'NEPOI','get_boundary_faces',nepoi,npoin)
  do ielem = 1,nelem
     ielty = ltype
     do inode = 1,element_type(ielty) % number_nodes
        ipoin = lnods(inode,ielem)
        nepoi(ipoin) = nepoi(ipoin) + 1
     end do
  end do
  call memory_alloca(memor_adapt,'PELPO','get_boundary_faces',pelpo,npoin+1_ip)
  pelpo(1) = 1
  do ipoin = 1,npoin
     pelpo(ipoin+1) = pelpo(ipoin) + nepoi(ipoin)
  end do
  nlelp = pelpo(npoin+1)
  call memory_alloca(memor_adapt,'LELPO','get_boundary_faces',lelpo,nlelp)
  do ielem = 1,nelem
     ielty = ltype
     do inode = 1,element_type(ielty) % number_nodes
        ipoin = lnods(inode,ielem)
        lelpo(pelpo(ipoin)) = ielem
        pelpo(ipoin) = pelpo(ipoin)+1
     end do
  end do
  pelpo(1) =  1
  mepoi    = -1
  do ipoin = 1,npoin
     pelpo(ipoin+1) = pelpo(ipoin) + nepoi(ipoin)
     mepoi = max(mepoi,nepoi(ipoin))
  end do
  call memory_deallo(memor_adapt,'NEPOI','get_boundary_faces',nepoi)

  !----------------------------------------------------------------------
  !
  ! List of global faces
  !
  !----------------------------------------------------------------------
  !
  ! Allocate memory for lelfa, FACES, CFAEL AND NNODG
  !
  call memory_alloca(memor_adapt,'FACEL','get_boundary_faces',facel,mnodf+1_ip,mface,nelem)
  !
  ! Construct and sort FACES
  !
  do ielem = 1,nelem                                          
     ielty = ltype!abs(ltype(ielem))
     do iface = 1,element_type(ielty) % number_faces
        pnodf = element_type(ielty) % node_faces(iface) 
        do inodf = 1,pnodf 
           inode = element_type(ielty) % list_faces(inodf,iface) 
           facel(inodf,iface,ielem) = lnods(inode,ielem)
        end do
        facel(mnodf+1,iface,ielem) = 1
        call maths_heap_sort(2_ip,pnodf,facel(:,iface,ielem))
     end do
  end do
  !
  ! Compute FACES
  !
  nfacg=0_ip
  do ielem = 1,nelem                                            ! Compare the faces and 
     ielty = ltype!abs(ltype(ielem))                                  ! eliminate the repited faces
     do iface = 1,element_type(ielty) % number_faces
        if( facel(mnodf+1,iface,ielem) > 0 ) then
           nfacg = nfacg + 1
           ipoin = facel(1,iface,ielem)
           ielpo = pelpo(ipoin)-1
           do while( ielpo < pelpo(ipoin+1)-1 )
              ielpo = ielpo + 1
              jelem = lelpo(ielpo)
              if( jelem /= ielem ) then
                 jelty = ltype!abs(ltype(jelem))                      ! eliminate the repited faces
                 jface = 0
                 do while( jface < element_type(jelty) % number_faces )
                    jface = jface + 1
                    if( facel(mnodf+1,jface,jelem) > 0 ) then
                       equal_faces = .true.
                       inodf = 0
                       do while( equal_faces .and. inodf < element_type(jelty) % node_faces(jface) )
                          inodf = inodf + 1 
                          if( facel(inodf,iface,ielem) /= facel(inodf,jface,jelem) ) equal_faces = .false.
                       end do
                       if( equal_faces ) then
                          facel(mnodf+1,iface,ielem) =  jelem                              ! Keep IELEM face
                          facel(mnodf+1,jface,jelem) = -ielem                              ! Elminate JELEM face
                          facel(      1,iface,ielem) = -jface                              ! Remember IFACE face
                          jface                      =  element_type(jelty) % number_faces ! Exit JFACE do
                          ielpo                      =  pelpo(ipoin+1)                     ! Exit JELEM do  
                       end if
                    end if
                 end do
              end if
           end do
        end if
     end do
  end do
  call memory_deallo(memor_adapt,'LELPO','get_boundary_faces',lelpo)
  call memory_deallo(memor_adapt,'PELPO','get_boundary_faces',pelpo)
  
  !----------------------------------------------------------------------
  !
  ! Count boundaries and nodes that are involved
  !
  !----------------------------------------------------------------------
  call memory_alloca(memor_adapt,'PERMR','get_boundary_faces',permr,npoin)
  call memory_alloca(memor_adapt,'INVPR','get_boundary_faces',invpr,npoin)
  kpoin = 0
  nboun = 0
  mnode = 0
  do ielem = 1,nelem                      
     ielty = ltype!abs(ltype(ielem))                    
     do iface = 1,element_type(ielty) % number_faces
        if( facel(mnodf+1,iface,ielem) > 0 .and. facel(1,iface,ielem) > 0 ) then
           pnodf = element_type(ielty) % node_faces(iface)
           mnode = max(mnode,pnodf)
           nboun = nboun + 1
           do inodf = 1,pnodf 
              inode = element_type(ielty) % list_faces(inodf,iface) 
              ipoin = lnods(inode,ielem)
              if( permr(ipoin) == 0 ) then
                 kpoin        = kpoin + 1
                 permr(ipoin) = kpoin
                 invpr(kpoin) = ipoin
              end if
           end do
        end if
     end do
  end do
  
  !----------------------------------------------------------------------
  !
  ! Fill in faces_bound
  !
  !----------------------------------------------------------------------
  nullify(lnods_bound)
  call memory_alloca(memor_adapt,'lnods_bound','mod_meshTopology',lnods_bound,lnnod_simplex-1,nboun)

  nboun = 0
  do ielem = 1,nelem                      
     ielty = ltype!abs(ltype(ielem))
     do iface = 1,element_type(ielty) % number_faces
        if( facel(mnodf+1,iface,ielem) > 0 .and. facel(1,iface,ielem) > 0 ) then
           nboun = nboun + 1
           pnodf = element_type(ielty) % node_faces(iface)
           !mesh_boun % ltype(nboun) = element_type(ielty) % type_faces(iface)
           !mesh_boun % perme(nboun) = ielem
           do inodf = 1,pnodf 
              inode = element_type(ielty) % list_faces(inodf,iface) 
              ipoin = lnods(inode,ielem)
              lnods_bound(inodf,nboun) = ipoin!permr(ipoin)
           end do
        end if
     end do
  end do
  if(lnnod_simplex==4) then 
    ! for tets, a tet is face_inverted + node (look at mod_elmgeo)
    ! in contrast to tris: edge+node=tri
    lnods_bound((/2,1/),:) = lnods_bound((/1,2/),:)
  end if
  !
  ! Deallocate memory
  !
  call memory_deallo(memor_adapt,'FACEL','get_boundary_faces',facel)
  call memory_deallo(memor_adapt,'PERMR','get_boundary_faces',permr)
  call memory_deallo(memor_adapt,'INVPR','get_boundary_faces',invpr)

end subroutine get_boundary_faces
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   delete_boundary_faces
!> @details delete_boundary_faces: tracking memory correctly
!
!-----------------------------------------------------------------------
subroutine delete_boundary_faces(lnods_bound)
  implicit none
  integer(ip), pointer ,intent(inout)    :: lnods_bound(:,:)
  !
  call memory_deallo(memor_adapt,'lnods_bound','mod_meshTopology',lnods_bound)
  !
end subroutine delete_boundary_faces
!
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   getIsBoundaryNode
!> @details getIsBoundaryNode
!
!-----------------------------------------------------------------------
subroutine getIsBoundaryNode(mesh,isBoundaryNode)
  implicit none
  
  type(mesh_type_basic),     intent(in)  :: mesh
  logical(lg), pointer,      intent(inout) :: isBoundaryNode(:)

  type(mesh_type_basic)  :: mesh_bou
    
  nullify(isBoundaryNode)
  call mesh_bou % init('BOUNDARY')
  call mesh_bou % boundary_mesh(mesh)
  call memory_alloca(memor_adapt,'isBoundaryNode','mod_meshTopology',isBoundaryNode,mesh%npoin)
  isBoundaryNode(:) = .false.
  isBoundaryNode(mesh_bou%permn) = .true.
  call mesh_bou%deallo()
  
end subroutine getIsBoundaryNode
!
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   getIsBoundaryNode
!> @details getIsBoundaryNode
!
!-----------------------------------------------------------------------
subroutine setIsBoundaryNodeFromIniMesh(isBoundaryNode,isBouNode_iniMesh,mapNodes_new_to_old)
  use mod_debugTools, only: out_performance,deb_isBouFromIni
  implicit none
  !
  logical(lg), pointer,      intent(inout) :: isBoundaryNode(:)
  logical(lg), pointer,      intent(in   ) :: isBouNode_iniMesh(:)
  integer(ip), pointer,      intent(in   ) :: mapNodes_new_to_old(:)
  !
  integer(ip) :: inode, numNod_new, inode_iniMesh
  !
  real(rp) :: t0,t1
  !
  if(out_performance) call cpu_time(t0)
  !
  numNod_new = size(mapNodes_new_to_old)
  
  nullify(isBoundaryNode)
  call memory_alloca(memor_adapt,'isBoundaryNode','mod_meshTopology',isBoundaryNode,numNod_new)
  isBoundaryNode(:) = .false.
  
  do inode=1,numNod_new
    inode_iniMesh = mapNodes_new_to_old(inode)
    if(inode_iniMesh>0) isBoundaryNode(inode) = isBouNode_iniMesh(inode_iniMesh)
  end do
  
  if(out_performance) then
    call cpu_time(t1)
    deb_isBouFromIni = deb_isBouFromIni + (t1-t0)
  end if
end subroutine setIsBoundaryNodeFromIniMesh
!
!
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   deleteIsBouNode array
!> @details deleteIsBouNode (deallocate tracking memor_adapt)
!
!-----------------------------------------------------------------------
subroutine deleteIsBouNode(isBoundaryNode)
  implicit none
  logical(lg), pointer,      intent(inout) :: isBoundaryNode(:)
  call memory_deallo(memor_adapt,'isBoundaryNode','mod_meshTopology',isBoundaryNode)
end subroutine deleteIsBouNode
!
!
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   delete_node_to_elems data structure
!> @details delete_node_to_elems (deallocate and nullify pointers)
!
!-----------------------------------------------------------------------
subroutine delete_node_to_nodes(node_to_nodes)!node_to_elems,node_to_firstElem,node_to_lastElem)
  implicit none

  class(node_to_nodes_type),       intent(inout)   :: node_to_nodes

  call memory_deallo(memor_adapt,'n_to_n','delete_node_to_nodes',node_to_nodes%n_to_n)
  call memory_deallo(memor_adapt,'n_to_fisrt_n','delete_node_to_nodes',node_to_nodes%n_to_fisrt_n)
  
  nullify(node_to_nodes%n_to_n)
  nullify(node_to_nodes%n_to_fisrt_n)
  !
end subroutine delete_node_to_nodes
!
!
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   compute_node_to_nodes
!> @details compute_node_to_nodes
!
!-----------------------------------------------------------------------
subroutine compute_node_to_nodes(node_to_nodes,meshe,node_to_elems)!node_to_elems,node_to_firstElem,node_to_lastElem)
  implicit none

  class(node_to_nodes_type),      intent(inout) :: node_to_nodes
  type(mesh_type_basic),          intent(in)    :: meshe
  type(node_to_elems_type),       intent(in)    :: node_to_elems

  !
  print*,'compute_node_to_nodes: todo'
  print*,'copy idea form compute_node_to_elems, but using node_to_elems data'
  call runend("program this function to use it..")
  !
end subroutine compute_node_to_nodes
!
!
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   get_elems_adj_to_node
!> @details get_elems_adj_to_node
!
!-----------------------------------------------------------------------
subroutine get_nodes_adj_to_node(node_to_nodes,inode,adjNodes)
  implicit none
  
  class(node_to_nodes_type), intent(in) :: node_to_nodes
  integer(ip), intent(in) :: inode
  
  integer(ip), pointer ,intent(inout) :: adjNodes(:)
  
  integer(ip) :: numAdjNodes
  
  numAdjNodes =  node_to_nodes%n_to_fisrt_n(inode+1) - node_to_nodes%n_to_fisrt_n(inode)

  nullify(adjNodes)
  call memory_alloca(memor_adapt,'elems','get_nodes_adj_to_node',adjNodes,numAdjNodes)

  adjNodes(:) = node_to_nodes%n_to_n(&
    node_to_nodes%n_to_fisrt_n(inode) : (node_to_nodes%n_to_fisrt_n(inode+1)-1) )  
    
end subroutine get_nodes_adj_to_node
!
!
!
END MODULE mod_meshTopology

!> @}




