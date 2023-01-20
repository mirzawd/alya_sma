!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Meshing
!> @{
!> @file    mod_surface_mesh.f90
!> @author  Abel Gargallo and houzeaux
!> @date    2020-09-18
!> @brief   Surface mesh
!> @details Tollbox for surface meshing
!-----------------------------------------------------------------------

module mod_surface_mesh

  use def_kintyp_basic,      only : ip,rp,lg
  use def_kintyp_mesh_basic, only : mesh_type_basic
  use mod_memory_basic,      only : memory_alloca
  use mod_memory_basic,      only : memory_deallo
  use def_elmtyp,            only : TRI03,QUA04,TRI06,QUA09,QUA08
  use mod_elmgeo,            only : element_type

  implicit none

  private
  
  integer(ip)               :: nVor                 ! number of voronoi cells
  integer(ip)               :: num3DNodes
  integer(ip)               :: numCoords
  integer(ip)               :: numElems

  integer(ip)               :: numSurfaceNodes
  integer(ip)               :: numEdges
  integer(ip)               :: numCoarseElems
  integer(ip)               :: numCoarseNodes

  integer(ip), allocatable  :: numCellsVor(:)       !numVor
  integer(ip), allocatable  :: elemVor(:)           !numElem
  logical(lg), allocatable  :: edgeMask(:)
  
  integer(ip), allocatable  :: edgeInfo(:,:)        !(2,nboun*mnodb) ! contains the adjacent elements
  integer(ip), allocatable  :: edgeVor(:,:)         !3,nedges (isBoundVor,Vi,Vj)
  integer(ip), allocatable  :: elementEdges(:,:)    !numNodesElem, numElem
  integer(ip), allocatable  :: nodeEdges(:,:)       !maxNumEdgesNode,numSurfaceNodes
  integer(ip), allocatable  :: edgeIsBound(:)       !numEdges

  integer(ip), allocatable  :: globToLoc(:)         !numSurfaceNodes
  integer(ip), allocatable  :: locToGlob(:)         !numSurfNodes
  integer(ip), allocatable  :: nodeVor(:,:)         !maxNumVoronoiNode,numSurfaceNodes
!  integer(ip), allocatable  :: nodeVorCon(:,:)      !nVor,numSurfNodes
  integer(ip), allocatable  :: numNodeVor(:)        !numSurfNodes

  real(rp),    allocatable  :: SRho(:)              !nvor
  real(rp),    allocatable  :: SGamma(:,:)          !3,nvor

  real(rp),    allocatable  :: areaElem(:)          !nelem
  real(rp),    allocatable  :: cmElem(:,:)          !3,nelem
  real(rp),    allocatable  :: centroid(:,:)        !3,nVor
  real(rp),    allocatable  :: normalToVor(:,:)     !3,nVor
  integer(ip), allocatable  :: countNumEdgesNode(:) !numSurfNodes

  integer(ip), allocatable  :: isEliminated(:)

  integer(ip), allocatable  :: mapValidVor(:)
  
  real(rp),    allocatable  :: coord_coarse_ref(:,:)   !3,nVor
  real(rp),    allocatable  :: coord_coarse_closest(:,:)   !3,nVor
  real(rp),    allocatable  :: coord_coarse_furthest(:,:)   !3,nVor
  
  integer(ip), allocatable  :: coarseNodesToFine(:)
    
  logical(lg) :: doOutputMeshes = .false.
 
  public :: surface_mesh_coarsening

contains 

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-18
  !> @brief   Surface mesh coarsening
  !> @details Coarsening of a surface mesh
  !>          Coarsening type: either give the number of nodes or compute them as a npoin/factor
  !>          "divide": npoin_coarse = npoin/optInetger, "npoin_coarse" numPoints:n_poinCoarse approx optInteger
  !> 
  !-----------------------------------------------------------------------

  subroutine surface_mesh_coarsening(&
       mesh,mesh_coarse,coarseningType,optInteger,lmask,permutation_CoarseToFine)

    class(mesh_type_basic),                    intent(inout)    :: mesh                  !< Input surface mesh
    class(mesh_type_basic),                    intent(inout) :: mesh_coarse           !< Input surface mesh
    logical(lg),            optional, pointer, intent(in)    :: lmask(:)              !> Boundary mask
    character(*),           optional,          intent(in)    :: coarseningType
    integer(ip),            optional,          intent(in)    :: optInteger
    integer(ip) ,allocatable,optional,         intent(inout) :: permutation_CoarseToFine(:)
    !integer(ip) ,optional,          intent(out)   :: permutation_fineToCoarse(:)
    
    logical(lg)                                              :: isNotFinished
    integer(ip)                                              :: countEmptyCells,numBoundVorEdges
    integer(ip)                                              :: iter,ndime,mnodb,ipoin       
    integer(ip)                                              :: nboun,npoin,ielem       
    integer(ip),                      pointer                :: lnodb(:,:)    
    integer(ip)                                              :: nFactor != 2
    integer(ip),                      pointer                :: ltypb(:)      
    real(rp),                         pointer                :: coord(:,:)    
    
    type(mesh_type_basic)                                    :: mesh_bou   ! boundary of input mesh
    
!    real :: t0, t1

    ndime =  mesh % ndime
    mnodb =  mesh % mnode
    nboun =  mesh % nelem
    npoin =  mesh % npoin
    lnodb => mesh % lnods
    ltypb => mesh % ltype
    coord => mesh % coord

    call mesh_coarse % init(trim(mesh % name)//'_COARSE')
    mesh_coarse % ndime = ndime
    mesh_coarse % mnode = mnodb

    if( present(coarseningType) ) then
       if( coarseningType == 'divide' ) then
          nFactor = optInteger
          nVor    = nboun/nFactor  
       else if( coarseningType == 'npoin_coarse' ) then
          nFactor = 0
          nVor    = optInteger
       end if
    else
       nFactor = 10
       nVor    = nboun/nFactor      
    end if

    numElems   = nboun
    numCoords  = ndime
    num3DNodes = npoin

    if( nboun > 0 .and. nVor > 0 ) then
       !call cpu_time(t0)
       !
       ! set connectivity of the nodes, edges,...
       !
       call setSurfaceMeshInfo(ndime,mnodb,nboun,npoin,lnodb,ltypb,coord,lmask)
       !
       ! initialize each voronoi with a random mesh element
       !
       call setVoronoiSeeds(nboun,ltypb,numBoundVorEdges,countEmptyCells)

       isNotFinished = .true.
       iter = 1
       do while( isNotFinished .and. iter < 10 )
          !
          ! Swap elements between voronoi cells to minimize energy
          !
          call minimizeEnergy(numBoundVorEdges,countEmptyCells,ltypb)
          isNotFinished = (countEmptyCells>0)
          iter = iter+1
       end do
       !
       ! find to which voronois does each surface node belong
       !
       !print*,'findVoronoiNodes'
       call findVoronoiNodes()
       !
       ! centroids of CVD can be off the surface: find closest vertex and assign to it the centroid
       !
       call generateCoarseNodes(ndime,mnodb,nboun,npoin,lnodb,ltypb,coord,mesh_coarse % npoin,mesh_coarse % coord)
       !
       ! generate the topology to connect the centroids of the CVD
       !
       call generateCoarseTopology(coord,mesh_coarse % nelem,mesh_coarse % lnods)
       !
       ! Improve the location of the nodes of the coarse mesh
       !
       !call cpu_time(t1)
       !print*,"TIME up to  improve: ", (t1-t0)
       
       if (doOutputMeshes) then 
         print*,'Output CVD'
         call out_inp2D('./VoronoiDiagram.inp',coord(1,:),coord(2,:),coord(3,:),lnodb(1:3,:),npoin,int(3,ip),nboun,real(mapValidVor(elemVor),rp))

         print*,'Output Coarse Mesh'
         call out_inp2D('./VorCoarseMesh.inp',mesh_coarse % coord(1,:),mesh_coarse % coord(2,:),mesh_coarse % coord(3,:),&
         mesh_coarse % lnods,mesh_coarse % npoin,int(3,ip),mesh_coarse % nelem,real(mapValidVor(elemVor),rp))
       end if

       ! ---> the following global mesh improvement has been deprecated: TOO EXPENSIVE FOR VISUALIZATION
!        print*,'improveCoarseNodesLocation'
!        call cpu_time(t0)
!        call improveCoarseNodesLocation(ndime,mnodb,nboun,npoin,lnodb,ltypb,coord,&
!          mesh_coarse % npoin,mesh_coarse % coord,mesh_coarse % nelem,mesh_coarse % lnods)
!        call cpu_time(t1)
!        print*,"TIME up to  improve: ", (t1-t0)/60
       
       ! ---> Instead we will only improve the bounraies of the coarse mesh
       !print*,'improveCoarseNodesLocation_boundary'
       !call cpu_time(t0)
       call mesh_bou % init('BOUNDARY_MESH_EDGES')
       call mesh_bou % boundary_mesh(mesh)
       
       call improveCoarseNodesLocation_boundary(ndime,mnodb,nboun,npoin,lnodb,ltypb,coord,&
           mesh_coarse % npoin,mesh_coarse % coord,mesh_coarse % nelem,mesh_coarse % lnods, mesh_bou)
           
       !print*,"TIME up to  improve: ", (t1-t0)
       
       if (doOutputMeshes) then 
         print*,'Output Coarse Mesh'
         call out_inp2D('./VorCoarseMesh_imp.inp',mesh_coarse % coord(1,:),mesh_coarse % coord(2,:),mesh_coarse % coord(3,:),&
         mesh_coarse % lnods,mesh_coarse % npoin,int(3,ip),mesh_coarse % nelem,real(mapValidVor(elemVor),rp))
       end if
       
       !if( allocated(mesh_bou % lnods)) call mesh_bou % deallo()
       !if( mesh_bou%nelem > 0) call mesh_bou % deallo()
       call mesh_bou % deallo()
       !
       call mesh_coarse % alloca(LTYPE=.true.,LNINV_LOC=.true.,LEINV_LOC=.true.,PERME=.true.,PERMN=.true.)
       do ielem = 1,mesh_coarse % nelem
          mesh_coarse % ltype(ielem)     = TRI03
          mesh_coarse % leinv_loc(ielem) = ielem
          mesh_coarse % perme(ielem)     = ielem
       end do
       do ipoin = 1,mesh_coarse % npoin
          mesh_coarse % lninv_loc(ipoin) = ipoin
          mesh_coarse % permn(ipoin)     = ipoin
       end do
       
       if (present(permutation_CoarseToFine)) then
         allocate(permutation_CoarseToFine(mesh_coarse % npoin))
         permutation_CoarseToFine(:)=coarseNodesToFine
       end if
       !
       ! Deallocate
       !
       call deallocateModuleVariables()
       
    end if


  end subroutine surface_mesh_coarsening
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  Abel Gargallo
  !> @date    2021-01-15
  !> @brief   Output paraview triangle mesh
  !> @details Output paraview for debugging (not called if doOutputMeshes=false)
  !> 
  !-----------------------------------------------------------------------
  
  subroutine out_inp2D(fname,x,y,z,lnods,npoin,nnode,nelem,elemfield)
  !*************************************************
  !*
  !*     Output for Paraview of the surf mesh
  !*
  !*************************************************
       implicit none
  !
  !*** input variables
  !
       character(len=*),intent(in):: fname
       integer(ip)    , intent(in):: npoin,nnode,nelem
       integer(ip)    , intent(in):: lnods(:,:)
       real(rp)       , intent(in):: x(npoin),y(npoin),z(npoin)
       real(rp), intent(in)       :: elemfield(nelem)
  !
  !*** local variables
  !
       integer(ip) :: ielem,inode
       real(rp) :: qelem
       character(len=4) :: elemType
  !
  !*** file to write
  !
       open(90,file=TRIM(fname),status='unknown')
  !
  !*** header
  !    
       if((nnode.eq.4).or.(nnode.eq.3)) then
  			 write(90,'(2(i8,1x),3(1x,i1))') npoin,nelem,0,1,0
       else
          call runend(' out_inp: incorrect number of nodes')
       end if
  !
  !*** coordinates
  !
      do inode = 1,npoin
         write(90,'(i8,3(1x,E15.7))') inode,x(inode),y(inode),z(inode)
      end do
  !
  !*** elements
  !
      do ielem = 1,nelem
        if(nnode==3) then
          elemType = 'tri '
          write(90,'(i8,1x,i1,1x,A,3(1x,i8))') ielem,0,elemType,(lnods(inode,ielem),inode=1,nnode)
        else if(nnode==4) then
          elemType = 'quad'
          write(90,'(i8,1x,i1,1x,A,4(1x,i8))') ielem,0,elemType,(lnods(inode,ielem),inode=1,nnode)
        end if
      end do       
  !
  !*** write element quality
  !
      write(90,'(2(1x,i1))') 1,1
      write(90,'(A)') 'q, None'
      do ielem = 1,nelem
        qelem = elemfield(ielem)
        write(90,'(i8,1x,F9.3)') ielem,qelem
      end do
      
      close(90)

  return  
  end subroutine out_inp2D
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  Abel Gargallo
  !> @date    2020-09-18
  !> @brief   Mesh info
  !> @details Set globToLoc: gives you the local SuRFace id of each node of the 3D mesh
  !>          Computes and stores areas and CM of elements
  !> 
  !-----------------------------------------------------------------------

  subroutine setSurfaceMeshInfo(ndime,mnodb,nboun,npoin,lnodb,ltypb,coord,lmask)
  
    implicit none
    integer(ip),                     intent(in) :: ndime                 !< Dimension (1,2,3)
    integer(ip),                     intent(in) :: mnodb                 !< Max nb of nodes per boundary element
    integer(ip),                     intent(in) :: nboun                 !< Number of boundaries
    integer(ip),                     intent(in) :: npoin                 !< Number of nodes
    integer(ip),                     intent(in) :: lnodb(mnodb,nboun)    !< Boundary connectivity
    integer(ip),                     intent(in) :: ltypb(nboun)          !< Boundary type
    real(rp),                        intent(in) :: coord(ndime,npoin)    !< Boundary connectivity
    logical(lg),  pointer, optional, intent(in) :: lmask(:)              !< Boundary mask
    
    integer(ip),  allocatable                   :: edgeInfoBig(:,:)      ! delete afterwards
    integer(ip),  allocatable                   :: nodeEdgesAux(:,:)
    integer(ip)                                 :: aux_i, gnode, ielem, inode, lnode, newEdge, newNode, nextNode
    integer(ip)                                 :: maxNumEdgesNode, numElemNodes, createdEdge, iedge, jelem

    !***allocate related global variables
    allocate(elemVor(nboun))
    elemVor = 0
    allocate(globToLoc(npoin))
    globToLoc = 0
    allocate(elementEdges(mnodb,nboun))
    elementEdges = 0

    ! CODE
    allocate(areaElem(nboun))
    allocate(cmElem(numCoords,nboun))
    globToLoc = 0
    newNode   = 0
    do ielem = 1,nboun
       numElemNodes = element_type(ltypb(ielem)) % number_nodes

       do inode = 1,numElemNodes
          gnode = lnodb(inode,ielem)
          if( globToLoc(gnode) == 0 ) then
             newNode          = newNode+1
             globToLoc(gnode) = newNode
          end if
       end do

       ! precompute the area and centroid of the elements
       areaElem(ielem)   = computeArea(numElemNodes, coord(:,lnodb(1:numElemNodes,ielem)) )
       cmElem(:,ielem)   = sum( coord(:,lnodb(1:numElemNodes,ielem)) ,2) / real(numElemNodes,rp)
    end do
    numSurfaceNodes = newNode

    ! define the edges of the SRF mesh and save the elements that contain them
    ! save for each node to which edges does it belong
    allocate(locToGlob(numSurfaceNodes))
    locToGlob = 0 
    allocate(edgeInfoBig(4,nboun*mnodb))
    edgeInfoBig = 0
    allocate(countNumEdgesNode(numSurfaceNodes))
    countNumEdgesNode = 0
    maxNumEdgesNode = 10 ! we consider this is the maximum. If not, we rebuild the array
    allocate(nodeEdges(maxNumEdgesNode,numSurfaceNodes))
    nodeEdges = 0
    !allocate(nodeCon(numSurfaceNodes,numSurfaceNodes))
    !nodeCon =0

    newEdge = 0
    do ielem=1,nboun
       numElemNodes = element_type(ltypb(ielem)) % number_nodes 
       do inode = 1,numElemNodes
          gnode = lnodb(inode,ielem)
          if( (globToLoc(gnode)>0).and.(locToGlob(globToLoc(gnode))==0) ) then
             newNode = globToLoc(gnode)
             locToGlob(newNode) = gnode
          end if
          lnode = globToLoc(gnode)
          if(inode<numElemNodes) then
             nextNode = globToLoc(lnodb(inode+1_ip,ielem))
          elseif(inode==numElemNodes) then
             nextNode = globToLoc(lnodb(1,ielem))
          end if

          createdEdge = 0
          aux_i = 1
          findCreatedEdge : do while( (aux_i<=maxNumEdgesNode).and.(nodeEdges(aux_i,lnode)/=0) )
             if(edgeInfoBig( 3,nodeEdges(aux_i,lnode) ) ==nextNode) then
                createdEdge = nodeEdges(aux_i,lnode) 
                exit findCreatedEdge
             end if
             aux_i  = aux_i + 1
          end do findCreatedEdge

          !if(nodeCon(lnode,nextNode)==0) then
          if(createdEdge==0) then
             newEdge                     = newEdge+1
             createdEdge                 = newEdge
             !nodeCon(lnode,nextNode) = newEdge
             !nodeCon(nextNode,lnode) = newEdge
             edgeInfoBig(1,newEdge )     = ielem
             edgeInfoBig(3,newEdge )     = lnode
             edgeInfoBig(4,newEdge )     = nextNode
             countNumEdgesNode(lnode)    = countNumEdgesNode(lnode)+1
             countNumEdgesNode(nextNode) = countNumEdgesNode(nextNode)+1

             aux_i = max(countNumEdgesNode(lnode),countNumEdgesNode(nextNode))
             if(aux_i>maxNumEdgesNode) then
                !call runend('no en aquesta malla')
                allocate(nodeEdgesAux(aux_i,numSurfaceNodes))
                nodeEdgesAux = 0
                nodeEdgesAux(1:maxNumEdgesNode,:) = nodeEdges

                deallocate(nodeEdges)
                maxNumEdgesNode = aux_i
                allocate(nodeEdges(maxNumEdgesNode,numSurfaceNodes))
                nodeEdges = 0
                nodeEdges(1:maxNumEdgesNode,:) = nodeEdgesAux
                deallocate(nodeEdgesAux)
             end if
             nodeEdges(countNumEdgesNode(lnode),lnode) = newEdge
             nodeEdges(countNumEdgesNode(nextNode),nextNode) = newEdge
          else
             edgeInfoBig(2,createdEdge ) = ielem
          end if
          elementEdges(inode,ielem) =  createdEdge !nodeCon(lnode,nextNode)
       end do
    end do
    numEdges = newEdge

    allocate(edgeInfo(2,numEdges))
    edgeInfo = edgeInfoBig(1:2,1:numEdges)

    deallocate(edgeInfoBig)
    !deallocate(nodeCon)
    !
    ! Mask
    !
    if( present(lmask) ) then
       allocate(edgeMask(numEdges))
       do iedge = 1,numEdges
          ielem           = edgeInfo(1,iedge)
          jelem           = edgeInfo(2,iedge)
          if( jelem > 0 ) then
             edgeMask(iedge) = lmask(ielem) .and. lmask(jelem)
          else
             edgeMask(iedge) = lmask(ielem)
          end if
       end do
    end if
    !
    return
  end subroutine setSurfaceMeshInfo

  !-----------------------------------------------------------------------
  !> 
  !> @author  Abel Gargallo
  !> @date    2020-09-18
  !> @brief   Voronoi seed
  !> @details Randomly set each Vi as a certain element Cj
  !>          Set the initial information of each Vi (SGamma, SRho)
  !> 
  !-----------------------------------------------------------------------

  subroutine setVoronoiSeeds(nboun,ltypb,numBoundVorEdges, countEmptyCells)
 
    integer(ip),      intent(in)  :: nboun                 !< Number of boundaries
    integer(ip),      intent(in)  :: ltypb(nboun)          !< Boundary type
    integer(ip),      intent(out) :: numBoundVorEdges
    integer(ip),      intent(out) :: countEmptyCells

    integer(ip)                   :: aux_i, ielem, ivor
    integer(ip)                   :: edgeId, iboundEdge, iedge, numElemNodes

    !*** allocate related global variables  
    !nVor = nboun/nFactor  
    !   allocate( voronoi(nboun , nVor) )
    !   voronoi = 0
    allocate(numCellsVor(nVor))
    numCellsVor = 0
    allocate(SGamma(numCoords,nVor))
    SGamma = 0.0_rp
    allocate(SRho(nVor))
    SRho = 0.0_rp

    !CODE
    allocate(edgeVor(3,numEdges))
    edgeVor = 0 
    allocate(edgeIsBound(numEdges))
    edgeIsBound = 0

    countEmptyCells = nboun
    ielem = 1 
    iboundEdge = 0
    aux_i = nboun/(nVor+1)
    do ivor = 1,nVor
       ielem = 1+(ivor-1)*aux_i
       elemVor(ielem) = ivor
       numCellsVor(ivor) = numCellsVor(ivor)+1
       !     voronoi(ielem,ivor) = 1
       countEmptyCells = countEmptyCells-1

       numElemNodes = element_type(ltypb(ielem)) % number_nodes 
       do iedge = 1,numElemNodes
          edgeId =  elementEdges(iedge,ielem)
          if(edgeInfo(2,edgeId).ne.0) then ! if the edge does not belong to a boundary
             if(edgeVor(1,edgeId )==0) then
                edgeVor( 1, edgeId ) =1
                edgeVor( 2, edgeId ) =ivor
                edgeVor( 3, edgeId ) =0  ! neighboring no voronoi (is going to expand)

                edgeIsBound(edgeId) = 1
             else
                edgeVor( 3, edgeId ) =ivor
             end if
          end if
       end do
       call updateVariablesVoronoi(+1_ip,ivor,ielem)
    end do

    numBoundVorEdges = sum(edgeIsBound)!iboundEdge

    return
  end subroutine setVoronoiSeeds


  !-----------------------------------------------------------------------
  !> 
  !> @author  Abel Gargallo
  !> @date    2021-02-12
  !> @brief   Minimize energy to define optimal voronoi cells
  !> @details Loop on edges and until convergence expand and reduce each Vi
  !> 
  !-----------------------------------------------------------------------

  subroutine minimizeEnergy(numBoundVorEdges,countEmptyCells,ltypb)

    implicit none
    integer(ip),      intent(inout) :: numBoundVorEdges
    integer(ip),      intent(inout) :: countEmptyCells
    integer(ip),      intent(in)    :: ltypb(numElems)          !< Boundary type

    !  internal variables
    integer(ip) :: numElemNodes, aux_i, aux_j
!    integer(ip) :: newEdge
    integer(ip) :: iedge, jedge, edgeId,  edgeIdChange, vor1, vor2, vorWin, vorLost, elemChanged, iter, countChanges
    integer(ip) :: elem1,elem2
    logical(lg) :: isNotFinished, hasEvolved, if_edge
    real(rp)    :: sgamma1a(numCoords),sgamma2a(numCoords),srho1a,srho2a
    real(rp)    :: sgamma1b(numCoords),sgamma2b(numCoords),srho1b,srho2b
    real(rp)    :: F0,F1,F2

    integer(ip),  allocatable  :: edgeIsBoundNew(:)

    allocate(edgeIsBoundNew(numEdges))
    iter = 0 
    isNotFinished = .true.
    do while( isNotFinished )!.and.(iter<200))
       iter = iter+1
       !print*,'Voronoi iter:',iter

       edgeIsBoundNew = 0
       !newEdge = 0  
       countChanges = 0  

       do iedge = 1,numEdges!numBoundVorEdges

          if( allocated(edgeMask) ) then
             if_edge = edgeMask(iedge)
          else
             if_edge = .true.
          end if

          edgeId = iedge
          if(edgeIsBound(edgeId)==1) then ! it can happen that this edge has been already revised by another edge
             !         if(edgeVor(1,edgeId)==0) call runend('something wrong 1') ! this should not happen (check prev if)
             vor1 = edgeVor(2,edgeId )
             vor2 = edgeVor(3,edgeId )

             if(vor1==0) then
                vor1 = vor2
                vor2 = 0 
             end if

             !         if(vor1==vor2) call runend('somethiing else skipped')

             if(vor2==0) then ! directly get the neighboring cell for vor1
                ! get adjacent element and set that edges as new boundaries(if not in the current boundary)
                vorWin = vor1
                vorLost = 0
                if(elemVor(edgeInfo(1,edgeId))==0) then
                   elemChanged = edgeInfo(1,edgeId)
                elseif(elemVor(edgeInfo(2,edgeId))==0) then
                   elemChanged = edgeInfo(2,edgeId)
                   !           else
                   !             call runend('something wrong 2')
                end if
                elemVor(elemChanged) =  vorWin
                numCellsVor(vorWin)  = numCellsVor(vorWin)+1
                !numCellsVor(vorLost) = numCellsVor(vorLost)-1
                hasEvolved           = .true.
                countChanges         = countChanges+1

                call updateVariablesVoronoi(+1_ip,vorWin,elemChanged)

                countEmptyCells = countEmptyCells-1

                if(countEmptyCells<0) call runend('something wrong negative')
                !---!!! THIS SHOULD BE UNCOMMENT TO SEPARE EATING VOIDS FROM ENERGY MINIMIZATION: FutureWork,Decide
                !           if(countEmptyCells==0) print*,'Void elements removed'
                !
                ! !           if(sum(voronoi(elemChanged,:) )>1 ) call runend('error in void eating')
                !
             elseif((vor2.ne.0).and.(countEmptyCells.ne.0)) then
                ! until we have not distributed the empty cells, don't minimize F
                ! check solution without this elseif ! futureWork, assumption
                hasEvolved = .false.
!!!!!
             else
!!! decide how to move the boundary (we want minimum energy)
                ! elem1 is in v1, elem2 is in v2 and we have to compute the changes
                elem1 = edgeInfo(1,edgeId)
                elem2 = edgeInfo(2,edgeId)

                if(elemVor(elem1)==vor2) then
                   aux_i = elem1
                   elem1 = elem2
                   elem2 = aux_i
                end if

                !Framework 1: function without doing anything
                F0       = goalF( SGamma(:,vor1),SRho(vor1),SGamma(:,vor2),SRho(vor2) )

                !Framework 2: vor1 eats elem2, vor2 loses elem2
                srho1a   = SRho(vor1) + areaElem(elem2)
                sgamma1a = SGamma(:,vor1) + areaElem(elem2)*(cmElem(:,elem2))
                srho2a   = SRho(vor2) - areaElem(elem2)
                sgamma2a = SGamma(:,vor2) - areaElem(elem2)*(cmElem(:,elem2))
                F1       = goalF( sgamma1a,srho1a,sgamma2a,srho2a )

                !Framework 3: vor1 loses elem1, vor2 eats elem1
                srho1b   = SRho(vor1) - areaElem(elem1)
                sgamma1b = SGamma(:,vor1) - areaElem(elem1)*(cmElem(:,elem1))
                srho2b   = SRho(vor2) + areaElem(elem1)
                sgamma2b = SGamma(:,vor2) + areaElem(elem1)*(cmElem(:,elem1))
                F2       = goalF( sgamma1b,srho1b,sgamma2b,srho2b )

                aux_i = maxloc( (/F0,F1,F2/), 1 )
                if(aux_i==1) then
                   hasEvolved = .false.
                else
                   hasEvolved   = .true.
                   countChanges = countChanges+1
                   if(aux_i==2) then
                      vorWin         = vor1
                      vorLost        = vor2
                      elemChanged    = elem2
                      SRho(vor1)     = srho1a
                      SRho(vor2)     = srho2a
                      SGamma(:,vor1) = sgamma1a
                      SGamma(:,vor2) = sgamma2a
                   elseif(aux_i==3) then
                      vorWin         = vor2
                      vorLost        = vor1
                      elemChanged    = elem1
                      SRho(vor1)     = srho1b
                      SRho(vor2)     = srho2b
                      SGamma(:,vor1) = sgamma1b
                      SGamma(:,vor2) = sgamma2b
                   end if
                   elemVor(elemChanged) = vorWin
                   numCellsVor(vorWin)  = numCellsVor(vorWin)+1
                   numCellsVor(vorLost) = numCellsVor(vorLost)-1
                end if
             end if
             if(hasEvolved) then
                numElemNodes = element_type(ltypb(elemChanged)) % number_nodes
                do jedge = 1,numElemNodes !loop on the edges of the changed elements
                   edgeIdChange = elementEdges(jedge,elemChanged)
                   if(edgeInfo(2,edgeIdChange).ne.0) then ! if the edge does not belong to a boundary
                      ! casos
                      if(edgeIdChange==edgeId) then ! the same edge we are cheching
                         edgeVor(1,edgeId )     = 0
                         edgeVor(2,edgeId )     = 0
                         edgeVor(3,edgeId )     = 0
                         edgeIsBound(edgeId)    = 0 ! we are dealing with it -> not any more boundary
                         edgeIsBoundNew(edgeId) = 0 ! we are dealing with it -> not any more boundary 
                      else
                         if(edgeVor(1,edgeIdChange)==0) then ! if it was not voronoi boundary, now it is 
                            edgeIsBoundNew(edgeIdChange) = 1
                            edgeVor(1,edgeIdChange ) = 1
                            edgeVor(2,edgeIdChange ) = vorWin
                            edgeVor(3,edgeIdChange ) = vorLost
                         else ! it was already a voronoi boundary between vorLost and others
                            if(  edgeVor(2,edgeIdChange)==vorLost) then
                               aux_i = 2
                               aux_j = 3
                            elseif(edgeVor(3,edgeIdChange)==vorLost) then
                               aux_i = 3
                               aux_j = 2
                               !                   else
                               !                     call runend('something wrong 3')
                            end if
                            if(edgeVor(aux_j,edgeIdChange)==vorWin) then ! if it is neighboring with itself
                               edgeVor(1,edgeIdChange ) = 0
                               edgeVor(2,edgeIdChange ) = 0
                               edgeVor(3,edgeIdChange ) = 0  
                               edgeIsBound(edgeIdChange) = 0
                               edgeIsBoundNew(edgeIdChange) = 0
                            else
                               edgeVor(aux_i,edgeIdChange) = vorWin !where we had vorLos we have vorWin
                               if(edgeIsBoundNew(edgeIdChange)==0) then
                                  edgeIsBoundNew(edgeIdChange) = 1
                               end if
                               !                     if(edgeVor(1,edgeIdChange )==0) call runend('problema aqui 1')
                            end if
                         end if
                      end if

                   end if
                end do

             else ! we always keep it to check in te next iteration
                ! it does not evolve but it is still a voronoi boundary edge
                !newEdge = newEdge +1
                !vorBoundEdgesNew(newEdge) = edgeId
                edgeIsBoundNew(edgeId) = 1
                !           if(edgeVor(1,edgeId )==0) call runend('problema aqui 2')
             end if
          end if

       end do

       isNotFinished    = (countChanges >0) !newEdges>0
       !vorBoundEdges = vorBoundEdgesNew
       edgeIsBound      = edgeIsBoundNew
       numBoundVorEdges = sum(edgeIsBound)!newEdge


       if(iter>1000) then
          print*,'Warning: could not find a proper coarsening of the mesh'
          isNotFinished = .false.
       end if
    end do

    !deallocate(vorBoundEdgesNew)
    deallocate(edgeIsBoundNew)

  end subroutine minimizeEnergy
  !
  !-----------------------------------------------------------------------
  !> 
  !> @author  Abel Gargallo
  !> @date    2021-02-12
  !> @brief   Find a voronoi cell for each node
  !> @details Assign voronoi to mesh vertices
  !> 
  !-----------------------------------------------------------------------
  !
  subroutine findVoronoiNodes()

    integer(ip)  :: aux_i, aux_j, edgeId, iedge, inode, ivor
    integer(ip)  :: maxNumVoronoiNode
    integer(ip)  :: belongVor(10), countVor, vorFound

    !*** CODE
    maxNumVoronoiNode = 8
    allocate(nodeVor(maxNumVoronoiNode,numSurfaceNodes))
    nodeVor = 0 

    allocate(numNodeVor(numSurfaceNodes))
    numNodeVor = 0

    !   allocate(nodeVorCon(nVor,numSurfaceNodes))
    !   nodeVorCon = 0

    do inode =1,numSurfaceNodes
       belongVor = 0
       countVor = 0
       do iedge = 1,countNumEdgesNode(inode)
          edgeId = nodeEdges(iedge,inode)
          if(edgeIsBound(edgeId)==1) then
             do aux_i=1,2
                ivor = edgeVor(aux_i+1_ip,edgeId) !locations 2 and 3 save the voronoi diagrams
                !if(nodeVorCon(ivor,inode)==0) then
                vorFound = 0
                findVor : do aux_j=1,countVor
                   if(belongVor(aux_j)==ivor) then
                      vorFound = 1
                      exit findVor
                   end if
                end do findVor
                if(vorFound==0) then
                   !nodeVorCon(ivor,inode) = 1
                   countVor = countVor+1
                   belongVor(countVor) = ivor
                   numNodeVor(inode) = numNodeVor(inode)+1
                   nodeVor(numNodeVor(inode),inode) =ivor
                end if
             end do
          end if
       end do
    end do

    !   deallocate(nodeVorCon)

    !   if(maxval(numNodeVor)>5) call runend(&
    !   'Not implemented when a node colides with more than 5 voronoi cells. Implement it: n vornonoi -> n-2 new triangles')
    return
  end subroutine findVoronoiNodes
  !
  !-----------------------------------------------------------------------
  !> 
  !> @author  Abel Gargallo
  !> @date    2021-02-12
  !> @brief   Set the geometrical location of the nodes of the coarse mesh
  !> @details The natural location would be the center of mass of the
  !>          voronoi cell. To preserve the geometry, we have chosen 
  !>          instead the closest node from the initial mesh. It could be
  !>          modified to choose the best node to minimize distance,
  !>          preserve the volume, or any other desired heuristic.
  !>          Main tasks performed by the function:
  !>          - Project the voronoi centroids to the mesh
  !>          - Compute a normal for each voronoi cell
  !> 
  !-----------------------------------------------------------------------
  !
  subroutine generateCoarseNodes(ndime,mnodb,nboun,npoin,lnodb,ltypb,coord,npoin_coarse,coord_coarse)
    
    implicit none
    !*** input-output
    integer(ip),          intent(in)    :: ndime                 !< Dimension (1,2,3)
    integer(ip),          intent(in)    :: mnodb                 !< Max nb of nodes per boundary element
    integer(ip),          intent(in)    :: nboun                 !< Number of boundaries
    integer(ip),          intent(in)    :: npoin                 !< Number of nodes
    integer(ip),          intent(in)    :: lnodb(mnodb,nboun)    !< Boundary connectivity
    integer(ip),          intent(in)    :: ltypb(nboun)          !< Boundary type
    real(rp),             intent(in)    :: coord(ndime,npoin)    !< Boundary connectivity
    integer(ip),          intent(out)   :: npoin_coarse          !< Number of nodes
    real(rp),    pointer, intent(inout) :: coord_coarse(:,:)     !< Boundary connectivity

    !***internal variables
    integer(ip)  :: ivor, ielem, inode, numElemNodes, selectedNode(nVor), selectedNodeFurthest(nVor)
    real(rp)     :: minDist(nVor), nodeDist,maxDist(nVor)
    integer(ip)  :: countValidVor
    real(rp)     :: x1(3),x2(3)
    real(rp)     :: normalToVorAux(3,nVor)

    normalToVorAux = 0

    minDist = huge(1.0_rp)
    maxDist = 0.0
    do ielem=1,nboun
       ivor = elemVor(ielem)
       if( ivor > 0 ) then
          if(numCellsVor(ivor)>0) then
             numElemNodes = element_type(ltypb(ielem)) % number_nodes
             do inode = 1,numElemNodes
                nodeDist = computeVectorNorm( SGamma(:,ivor)/SRho(ivor) -coord(:,lnodb(inode,ielem)))
                if( minDist(ivor) > nodeDist ) then
                   minDist(ivor) = nodeDist
                   selectedNode(ivor) = lnodb(inode,ielem)
                end if
                if( maxDist(ivor) < nodeDist ) then
                   maxDist(ivor) = nodeDist
                   selectedNodeFurthest(ivor) = lnodb(inode,ielem)
                end if
             end do

             x1                     = coord(:,lnodb(2,ielem))-coord(:,lnodb(1,ielem))
             x2                     = coord(:,lnodb(3,ielem))-coord(:,lnodb(1,ielem))
             normalToVorAux(:,ivor) = normalToVorAux(:,ivor)+crossProduct(x1,x2)
          end if
       end if
    end do

    allocate(mapValidVor(nVor))
    countValidVor = 0
    do ivor =1,nVor
       if(numCellsVor(ivor)>0) then
          countValidVor = countValidVor+1
          mapValidVor(ivor) = countValidVor
       end if
    end do

    numCoarseNodes = countValidVor
    npoin_coarse = numCoarseNodes
    allocate(coord_coarse(numCoords,numCoarseNodes))
    allocate(centroid(numCoords,countValidVor))
    allocate(normalToVor(3,countValidVor))
    allocate(coord_coarse_ref(numCoords,numCoarseNodes))
    allocate(coord_coarse_closest(numCoords,numCoarseNodes))
    allocate(coord_coarse_furthest(numCoords,numCoarseNodes))
    allocate(coarseNodesToFine(npoin_coarse))
    do ivor =1,nVor
       if(numCellsVor(ivor)>0) then      
          !*** compute centroid
          centroid(:,mapValidVor(ivor)) = SGamma(:,ivor)/SRho(ivor)

          !*** generate coordinates corresponding node: futureWork, Assumption, Decide: on the initial mesh or not!!!
          !coord_coarse(:,mapValidVor(ivor)) = centroid(:,mapValidVor(ivor)) !can be off the surf but generates better quality mesh
          coord_coarse(:,mapValidVor(ivor)) = coord(:,selectedNode(ivor)) ! new coordinate on the initial mesh

          coord_coarse_ref(:,mapValidVor(ivor)) = centroid(:,mapValidVor(ivor))
          coord_coarse_closest(:,mapValidVor(ivor)) = coord(:,selectedNode(ivor))
          coord_coarse_furthest(:,mapValidVor(ivor)) = coord(:,selectedNodeFurthest(ivor))

          normalToVor(:,mapValidVor(ivor)) = normalToVorAux(:,ivor)/computeVectorNorm(normalToVorAux(:,ivor))
          
          coarseNodesToFine(mapValidVor(ivor)) = selectedNode(ivor)
       end if
    end do

    return
  end subroutine generateCoarseNodes
  !
  !-----------------------------------------------------------------------
  !> 
  !> @author  Abel Gargallo
  !> @date    2021-02-12
  !> @brief   Set the connectivity of the coarse mesh
  !> @details Build the coarse topology connecting the centroids of the CVD,
  !>          dual of the voronoi cells 
  !>          Abel: I have further improvements in my notes...
  !> 
  !-----------------------------------------------------------------------
  !
  subroutine generateCoarseTopology(coord,nboun_coarse,lnodb_coarse)
    implicit none
    ! input-output variables
    real(rp),        intent(in)  :: coord(numCoords,num3DNodes)    !< Boundary connectivity
    integer(ip),      intent(out) :: nboun_coarse          !< Number of boundaries
    integer(ip),   pointer,    intent(inout) :: lnodb_coarse(:,:)     !< Boundary connectivity
    ! internal variables
    integer(ip)      :: inode, aux_i, ivor, jnode, aux_k, minLoci!, countPos
    integer(ip)      :: newElem,auxPerm(5,6)
    real(rp)         :: x1(3),x2(3),aux_vec(3),angleMax,minDist,dist
    real(rp),  allocatable  :: angles(:)
    integer(ip),allocatable  :: reorderNodes(:)
    integer(ip),allocatable  :: marked(:)
!    logical        :: aux_log
    real(rp) :: tangent1(3), tangent2(3),normTangent1,normTangent2,normalToPatch(3)

    ! CODE
    numCoarseElems = 0
    do inode = 1,numSurfaceNodes
       if(numNodeVor(inode)>2) numCoarseElems = numCoarseElems + (numNodeVor(inode)-2)
    end do
    allocate(lnodb_coarse(3,numCoarseElems))
    newElem = 0
    do inode =1,numSurfaceNodes
       !     if(numNodeVor(inode)>2) then
       !       newElem = newElem+1
       !       lnodb_coarse(:,newElem) = nodeVor(1:3,inode)
       !end if  
       if(numNodeVor(inode)==3) then
          newElem = newElem+1
          lnodb_coarse(:,newElem) = mapValidVor( nodeVor(1:numNodeVor(inode),inode) )
          ! properly orientate
          aux_vec = 0.0_rp ! mean of the normals on the node
          do ivor=1,numNodeVor(inode)
             aux_vec=aux_vec+normalToVor(:, mapValidVor( nodeVor(ivor,inode) ) )
          end do
          aux_vec = aux_vec/computeVectorNorm(aux_vec)
          x1 = centroid(:,lnodb_coarse(2,newElem))-centroid(:,lnodb_coarse(1,newElem))
          x2 = centroid(:,lnodb_coarse(3,newElem))-centroid(:,lnodb_coarse(1,newElem))
          if(sum(crossProduct(x1,x2)*aux_vec)<0) lnodb_coarse(2:3,newElem) = lnodb_coarse((/3,2/),newElem)
       else if(numNodeVor(inode)>3) then 
          ! order the nodes
          ! compute the edges and compute the angles, and take as next the one with lowest angle

          ! compute the normal of the node to be able to orientate the elements with the proper normal
          aux_vec = 0.0_rp ! mean of the normals on the node
          do ivor=1,numNodeVor(inode)
             aux_vec=aux_vec+normalToVor(:,mapValidVor( nodeVor(ivor,inode) ))
          end do
          aux_vec = aux_vec/computeVectorNorm(aux_vec)
          normalToPatch = aux_vec

          ! Compute the best projection plane for the patch
          tangent1(1)   =   normalToPatch(3)
          tangent1(2)   =  0.0_rp
          tangent1(3)   = - normalToPatch(1)
          normTangent1=  computeVectorNorm(tangent1)
          if(normTangent1.le.(1.0e-6_rp)) then
             tangent1(1)   =   0.0_rp
             tangent1(2)   =  normalToPatch(3)
             tangent1(3)   = - normalToPatch(2)
             normTangent1=  computeVectorNorm(tangent1)
             if(normTangent1.le.(1.0e-6_rp)) then
                tangent1(1)   =   normalToPatch(2)
                tangent1(2)   = - normalToPatch(1)
                tangent1(3)   = 0.0_rp
                normTangent1=  computeVectorNorm(tangent1)
             end if
          end if
          tangent1 = tangent1/normTangent1   
          !call cross(normalToPatch,tangent1,tangent2)
          tangent2 = crossProduct(normalToPatch,tangent1)
          normTangent2=  computeVectorNorm(tangent2)
          tangent2 = tangent2/normTangent2

          !         allocate( marked( numNodeVor(inode)) )
          !         marked = 0
          !         allocate( reorderNodes( numNodeVor(inode)) )
          !         allocate(angles(numNodeVor(inode)))
          !         ! find optimal node to start with
          !         aux_log = .true.
          !         aux_i = 1
          !         do while(aux_log)
          !           x1 = centroid(:,mapValidVor( nodeVor(aux_i,inode) )) - coord(:,locToGlob(inode))
          !           x1 = x1/computeVectorNorm(x1)
          !           angles = -10
          !           countPos = 0
          !           do aux_k= 1 ,numNodeVor(inode)
          !             if(aux_k.ne.aux_i) then
          !               x2 = centroid(:,mapValidVor( nodeVor(aux_k,inode) )) - coord(:,locToGlob(inode))
          !               x2 = x2/computeVectorNorm(x2)
          !               if(sum(x1*x2)>0) countPos = countPos+1
          !             end if
          !           end do
          !           if(countPos.ge.2) then
          !             aux_log = .false.
          !           else
          !             aux_i = aux_i+1
          !           end if
          !         end do
          !         marked(aux_i) = 1
          !         reorderNodes(1) = mapValidVor( nodeVor(aux_i,inode) )    

          !get the ordered nodes of the patch
          if(numNodeVor(inode) ==4) then
             auxPerm(:,1) = (/1,2,3,4,1/)
             auxPerm(:,2) = (/1,2,4,3,1/)
             auxPerm(:,3) = (/1,3,2,4,1/)
             auxPerm(:,4) = (/1,3,4,2,1/)
             auxPerm(:,5) = (/1,4,2,3,1/)
             auxPerm(:,6) = (/1,4,3,2,1/)
             minDist = HUGE(real(rp))
             do aux_k= 1 ,6
                dist = 0.0_rp
                do aux_i= 2 ,5!(numNodeVor(inode)+1)
                   !             dist =dist+computeVectorNorm(&
                   !             centroid(:,mapValidVor( nodeVor(auxPerm(aux_i-1,aux_k),inode) )) -&
                   !             centroid(:,mapValidVor( nodeVor(auxPerm(aux_i   ,aux_k),inode) )) )
                   dist =dist+computeVectorNorm(&
                        projectToTangent(&
                        centroid(:,mapValidVor( nodeVor(auxPerm(aux_i-1,aux_k),inode) )) -&
                        centroid(:,mapValidVor( nodeVor(auxPerm(aux_i   ,aux_k),inode) )),&
                        tangent1,tangent2))
                end do
                if(minDist>dist) then
                   minLoci = aux_k
                   minDist = dist
                end if
             end do
             allocate( reorderNodes( 4 ) ) ! numNodeVor(inode)
             reorderNodes(:) = mapValidVor( nodeVor(auxPerm(1:4,minloci),inode) )
             !         minDist = HUGE(real(rp))
             !         do aux_k= 1 ,6
             !           dist = 0.0_rp
             !           do aux_i= 2 ,numNodeVor(inode)
             !             dist =dist+computeVectorNorm(&
             !             centroid(:,mapValidVor( nodeVor(auxPerm(aux_i-1,aux_k),inode) )) -&
             !             centroid(:,mapValidVor( nodeVor(auxPerm(aux_i   ,aux_k),inode) )) )
             !           end do
             !           if(minDist>dist) then
             !             minloci = aux_k
             !             minDist = dist
             !           end if
             !         end do
             !         allocate( reorderNodes( 4 ) ) ! numNodeVor(inode)
             !         reorderNodes(:) = mapValidVor( nodeVor(auxPerm(:,minloci),inode) )

          else
             allocate( marked( numNodeVor(inode)) )
             marked = 0
             marked(1) = 1
             allocate( reorderNodes( numNodeVor(inode)) )
             reorderNodes(1) = mapValidVor( nodeVor(1,inode) )
             allocate(angles(numNodeVor(inode)))
             do jnode = 2,numNodeVor(inode)
                x1 = centroid(:,reorderNodes(jnode-1)) - coord(:,locToGlob(inode))
                x1 = x1/computeVectorNorm(x1)
                x1 = projectToTangent(x1,tangent1,tangent2)
                angles = -10
                angleMax = -10
                do aux_k= 1 ,numNodeVor(inode)
                   if(marked(aux_k)==0) then
                      x2 = centroid(:,mapValidVor( nodeVor(aux_k,inode) )) - coord(:,locToGlob(inode))
                      x2 = x2/computeVectorNorm(x2)
                      x2 = projectToTangent(x2,tangent1,tangent2)
                      angles(aux_k) = sum(x1*x2)
                      if(angleMax<angles(aux_k)) then ! added to avoid the use of maxloc for JuanCarlos (problem in IBM)
                         angleMax=angles(aux_k)
                         aux_i=aux_k
                      end if
                   end if
                end do
                !aux_i = maxloc(angles,numNodeVor(inode))
                reorderNodes(jnode) = mapValidVor( nodeVor(aux_i,inode) )
                marked(aux_i) = 1
             end do
             deallocate(angles)
             deallocate(marked)
          end if

          ! generate new elements
          !       if((numNodeVor(inode)>4) ) then
          !       print*,centroid(:,reorderNodes(1))
          !       print*,centroid(:,reorderNodes(2))
          !       print*,centroid(:,reorderNodes(3))
          !       print*,centroid(:,reorderNodes(4))
          !       print*,centroid(:,reorderNodes(5))
          !       print*,coord(:,locToGlob(inode))
          !       end if
          do jnode = 2,(numNodeVor(inode)-1)
             newElem = newElem+1
             lnodb_coarse(:,newElem) =  (/ reorderNodes(1),reorderNodes(jnode),reorderNodes(jnode+1) /)
             ! reorder the element nodes so that the normal is well oriented
             x1 = centroid(:,lnodb_coarse(2,newElem))-centroid(:,lnodb_coarse(1,newElem))
             x2 = centroid(:,lnodb_coarse(3,newElem))-centroid(:,lnodb_coarse(1,newElem))
             if(sum(crossProduct(x1,x2)*aux_vec)<0.0_rp) lnodb_coarse(2:3,newElem) = lnodb_coarse((/3,2/),newElem)
          end do
          deallocate(reorderNodes)  

       end if

    end do

    nboun_coarse= numCoarseElems
    return  
  end subroutine generateCoarseTopology
  !
  !
  !-----------------------------------------------------------------------
  !> 
  !> @author  Abel Gargallo
  !> @date    2021-02-12
  !> @brief   Improve geometrical configuration of the nodes
  !> @details Loop on original nodes to find the best one to match
  !> 
  !-----------------------------------------------------------------------
  !
  subroutine improveCoarseNodesLocation_boundary(ndime,mnodb,nboun,npoin,lnodb,ltypb,coord,&
    npoin_coarse,coord_coarse,nelem_coarse,lnods_coarse,mesh_bou)
    
    implicit none
    !*** input-output
    integer(ip),          intent(in)    :: ndime                 !< Dimension (1,2,3)
    integer(ip),          intent(in)    :: mnodb                 !< Max nb of nodes per boundary element
    integer(ip),          intent(in)    :: nboun                 !< Number of boundaries
    integer(ip),          intent(in)    :: npoin                 !< Number of nodes
    integer(ip),          intent(in)    :: lnodb(mnodb,nboun)    !< Boundary connectivity
    integer(ip),          intent(in)    :: ltypb(nboun)          !< Boundary type
    real(rp),             intent(in)    :: coord(ndime,npoin)    !< Boundary connectivity
    integer(ip),          intent(in)    :: npoin_coarse          !< Number of nodes
    real(rp),    pointer, intent(inout) :: coord_coarse(:,:)     !< Coarse boundary coordinates
    integer(ip),          intent(in)    :: nelem_coarse          !< Number of nodes
    integer(ip), pointer, intent(in)    :: lnods_coarse(:,:)     !< Coarse boundary connectivity
    class(mesh_type_basic),intent(in)   :: mesh_bou              !< Boundary of the mesh lnodb
    !
    !***internal variables
    integer(ip) :: ielem_bou, ielem, ivor, numElemNodes!, n1, n2
    real(rp)    ::  function_to_minimize(nVor), my_fun, target_x(ndime)
    
    integer(ip) :: nodeToElems(100_ip,npoin), nodeToNumElems(npoin)
    integer(ip) :: inode, ipoin
    real(rp) :: prev_x(ndime)

    integer(ip) :: numImprovementLoops, iloop
    numImprovementLoops = 2
    
    if (ndime.ne.3) then
      return
    end if
    
    nodeToNumElems(:) = 0_ip
    do ielem=1,nelem_coarse
      do inode=1,3 !tri mesh
        ipoin = lnods_coarse(inode,ielem)
        nodeToNumElems(ipoin) = nodeToNumElems(ipoin) + 1
        if(nodeToNumElems(ipoin) > 100_ip) then
          call runend("This should not happen.. node in triangle mesh with over 100 elements connected...")
        end if
        nodeToElems(nodeToNumElems(ipoin),ipoin) = ielem
      end do      
    end do
        
    function_to_minimize(:) = huge(0.0)
    
    do iloop = 1,numImprovementLoops
      
    do ielem_bou=1,mesh_bou % nelem
      ielem = mesh_bou % perme(ielem_bou)
      ivor = mapValidVor( elemVor(ielem) )

!       n1 = mesh_bou % lnods(1,ielem_bou)
!       n2 = mesh_bou % lnods(2,ielem_bou)
!
!       !target_x =  (mesh_bou % coord(:, n1 ) + mesh_bou % coord(:, n2 )) / 2.0
!       target_x =  mesh_bou % coord(:, n1 )
      
      do inode = 1,2
        prev_x = coord_coarse(:,ivor)
        
        ipoin = mesh_bou % lnods(inode,ielem_bou)
        target_x = mesh_bou % coord(:, ipoin )
      
        !----- By distance
        !my_fun = norm2(coord_coarse_ref(:,ivor) - target_x)
        !! to maximize distance, use -d

        !---- By some criteria to help different subdomains get the same node
        !my_fun = sum((target_x))!abs
      
        !---- By increment of area<
        coord_coarse(:,ivor)=target_x
        my_fun = 0.0 !implemented in computeAreaMesh_noReference=> NOT WORKING, area of elements adjacent to node
        do ielem=1,nodeToNumElems(ivor)
          numElemNodes = size(lnods_coarse,1)
          my_fun = my_fun - computeArea( numElemNodes , coord_coarse(:,  lnods_coarse(:,  nodeToElems(ielem,ivor)  )  )  )
        end do

        !---- Discrete minimization of the desired function
        if( my_fun < function_to_minimize(ivor)) then
          function_to_minimize(ivor)=my_fun
          coord_coarse(:,ivor)=target_x
          coarseNodesToFine(ivor) =  mesh_bou % permn(ipoin)
        else
          coord_coarse(:,ivor) = prev_x
        end if
      end do 
      
    end do
    
    end do
    
    return
  end subroutine improveCoarseNodesLocation_boundary
  !-----------------------------------------------------------------------
  !> 
  !> @author  Abel Gargallo
  !> @date    2021-02-12
  !> @brief   Improve geometrical configuration of the nodes
  !> @details Loop on original nodes to find the best one to match
  !> 
  !-----------------------------------------------------------------------
  !
  subroutine improveCoarseNodesLocation(ndime,mnodb,nboun,npoin,lnodb,ltypb,coord,&
    npoin_coarse,coord_coarse,nelem_coarse,lnods_coarse)
    
    implicit none
    !*** input-output
    integer(ip),          intent(in)    :: ndime                 !< Dimension (1,2,3)
    integer(ip),          intent(in)    :: mnodb                 !< Max nb of nodes per boundary element
    integer(ip),          intent(in)    :: nboun                 !< Number of boundaries
    integer(ip),          intent(in)    :: npoin                 !< Number of nodes
    integer(ip),          intent(in)    :: lnodb(mnodb,nboun)    !< Boundary connectivity
    integer(ip),          intent(in)    :: ltypb(nboun)          !< Boundary type
    real(rp),             intent(in)    :: coord(ndime,npoin)    !< Boundary connectivity
    integer(ip),          intent(in)    :: npoin_coarse          !< Number of nodes
    real(rp),    pointer, intent(inout) :: coord_coarse(:,:)     !< Coarse boundary coordinates
    integer(ip),          intent(in)    :: nelem_coarse          !< Number of nodes
    integer(ip), pointer, intent(in)    :: lnods_coarse(:,:)     !< Coarse boundary connectivity
    !
    !***internal variables
    logical(lg) :: isModified
    integer(ip) :: numElemNodes, countAux
    integer(ip) :: iter, ipoin_coarse, ipoin, theVoronoi,ielem,inode
    real(rp)    :: prev_x(ndime), x_candidate(ndime)
    real(rp)    :: area_coarse, area_coarse_new, area_target
    logical(lg) :: isInverted
        
    !print*,"area_target"
    area_target = computeAreaMesh_noPoin(coord,lnodb,ltypb)
    !print*,area_target
    !print*,"area_coarse"
    area_coarse = computeAreaMesh(coord_coarse,lnods_coarse,isInverted)
    !print*,area_coarse
    
    
    isModified = .true.
    iter = 0
    do while( isModified .and. iter < 1 ) !should finalize in a few iters -> 1 for efficicency
      isModified = .false.
      
      countAux=0
      do ielem =1,nboun          
        theVoronoi = elemVor(ielem)
        ipoin_coarse = mapValidVor(theVoronoi)
          
        numElemNodes =  element_type(ltypb(ielem)) % number_nodes
        !x_candidate(:) =  sum(coord(:,lnodb(1:numElemNodes,ielem)),2)/numElemNodes
        
        do inode = 1,numElemNodes

          prev_x = coord_coarse(:,ipoin_coarse)
          ipoin = lnodb(inode,ielem)
          
          if (numNodeVor(ipoin).le.1) then
          
            x_candidate = coord(:,ipoin)
            coord_coarse(:,ipoin_coarse) = x_candidate
            area_coarse_new = computeAreaMesh(coord_coarse,lnods_coarse,isInverted)! TODO: improve with just elements adjacent to the node
        
            isModified = isImproved(area_coarse_new,area_coarse,area_target,isInverted) 
            if (isModified) then
              area_coarse = area_coarse_new
            else
              coord_coarse(:,ipoin_coarse) = prev_x
            end if
          end if
        end do
      end do
      iter = iter+1
    end do
    
    ! try to put it on the closest point on the surface (for planar improvement)
!     do ipoin_coarse=1,npoin_coarse
!       prev_x = coord_coarse(:,ipoin_coarse)
!       coord_coarse(:,ipoin_coarse)=coord_coarse_furthest(:,ipoin_coarse)
!       area_coarse_new = computeAreaMesh(coord_coarse,lnods_coarse,isInverted)! TODO: improve with just elements adjacent to the node
!       isModified = isImproved(area_coarse_new,area_coarse,area_target,isInverted)
!       if (isModified) then
!         area_coarse = area_coarse_new
!       else
!         coord_coarse(:,ipoin_coarse) = prev_x
!       end if
!     end do
    ! try to put it on the furthest point on the surface (for planar improvement)
!     do ipoin_coarse=1,npoin_coarse
!       prev_x = coord_coarse(:,ipoin_coarse)
!       coord_coarse(:,ipoin_coarse)=coord_coarse_closest(:,ipoin_coarse)
!       area_coarse_new = computeAreaMesh(coord_coarse,lnods_coarse,isInverted)! TODO: improve with just elements adjacent to the node
!       isModified = isImproved(area_coarse_new,area_coarse,area_target,isInverted)
!       if (isModified) then
!         area_coarse = area_coarse_new
!       else
!         coord_coarse(:,ipoin_coarse) = prev_x
!       end if
!     end do
    
    
    
    ! try center of mass of the elements
!     isModified = .true.
!     iter = 0
!     do while( isModified .and. iter < 100 ) !should finalize in a few iters
!       isModified = .false.
!
!       countAux=0
!       do ielem =1,nboun
!         theVoronoi = elemVor(ielem)
!         ipoin_coarse = mapValidVor(theVoronoi)
!
!         numElemNodes =  element_type(ltypb(ielem)) % number_nodes
!         x_candidate(:) =  sum(coord(:,lnodb(1:numElemNodes,ielem)),2)/numElemNodes
!
!         prev_x = coord_coarse(:,ipoin_coarse)
!         area_coarse_new = computeAreaMesh(coord_coarse,lnods_coarse,isInverted)! TODO: improve with just elements adjacent to the node
!
!         !isModified = isImproved(area_coarse_new,area_coarse,area_target,isInverted)
!         !if (isModified) then
!         if((area_coarse_new>area_coarse).and.(.not.isInverted)) then
!           area_coarse = area_coarse_new
!         else
!           coord_coarse(:,ipoin_coarse) = prev_x
!         end if
!
!       end do
!       iter = iter+1
!     end do
    
    !print*,"Area: ",area_coarse
    
    return
  end subroutine improveCoarseNodesLocation
  !
  !
  !
  function isImproved(area_coarse_new,area_coarse,area_target,isInverted) result(positiveImprovement)
    implicit none
    real(rp), intent(in)::area_coarse_new,area_coarse,area_target
    logical(lg), intent(in)::isInverted
    logical(lg) :: positiveImprovement
    
    !if ((area_coarse_new > area_coarse)) then
    !if ((area_coarse_new > area_coarse).and.(area_coarse_new.le.area_target)) then
    !if ( abs(area_target-area_coarse_new) < abs(area_target-area_coarse)) then
    
    positiveImprovement = ((abs(area_target-area_coarse_new) < abs(area_target-area_coarse)).and.(.not.isInverted))
    
    return
  end function isImproved
  !
  !
  !
  function computeAreaMesh_noPoin(coord,lnods,ltypb) result(area)
    
    implicit none
    !*** input-output
    real(rp),    intent(in) :: coord(:,:)     !< Coarse boundary coordinates
    integer(ip), intent(in) :: lnods(:,:)     !< Coarse boundary connectivity
    integer(ip),          intent(in)    :: ltypb(:)          !< Boundary type
    
    
    real(rp) :: area
    
    real(rp) :: coordsElem(size(coord,1),3)
    integer(ip) :: numElemNodes, ielem
    
    area = 0.0
    do ielem=1,size(lnods,2)
      numElemNodes = element_type(ltypb(ielem)) % number_nodes
      coordsElem = coord(:, lnods(1:numElemNodes,ielem) )
      area = area + computeArea( numElemNodes , coordsElem )
      !print*,ielem,"   -> area: ",computeArea( size(lnods,1) , coordsElem )
    end do
    
!     print*,"computeAreaMesh_noPoin: ",area
    
    return
  end function computeAreaMesh_noPoin
  
  function computeAreaMesh(coord,lnods,isInverted) result(area)
    
    implicit none
    !*** input-output
    real(rp),    pointer, intent(in) :: coord(:,:)     !< Coarse boundary coordinates
    integer(ip), pointer, intent(in) :: lnods(:,:)     !< Coarse boundary connectivity
    logical(lg), intent(out)::isInverted ! some element is inverted
    
    real(rp) :: area
    
    real(rp) :: coordsElem(size(coord,1),size(lnods,1)), areaElem, normalRef(3)
    integer(ip) :: numElemNodes, ielem
    
    area = 0.0
    isInverted = .false.
    do ielem=1,size(lnods,2)
      numElemNodes = size(lnods,1)
      
      coordsElem = coord_coarse_ref(:, lnods(:,ielem) )
      normalRef = computeNormal(numElemNodes,coordsElem) ! normal from the nice voronoi mesh
!       normalRef = sum(normalToVor(:, lnods(:,ielem)),2)  ! average normal from the neighbor cells
!       normalRef = normalRef/computeVectorNorm(normalRef)
      
      coordsElem = coord(:, lnods(:,ielem) )
!       areaElem = computeArea( numElemNodes , coordsElem )
      areaElem = computeArea_withSign( numElemNodes , coordsElem, normalRef )
      area = area + areaElem
      !print*,ielem,"   -> area: ",computeArea( size(lnods,1) , coordsElem )
      if (areaElem.le.0) then
        isInverted = .true.
      end if
    end do
    
!     print*,area
!     call runend('hiha un error aqui o algo')
    
    return
  end function computeAreaMesh
  !
  function computeAreaMesh_noReference(coord,lnods) result(area)
    
    implicit none
    !*** input-output
    real(rp),    pointer, intent(in) :: coord(:,:)     !< Coarse boundary coordinates
    integer(ip), pointer, intent(in) :: lnods(:,:)     !< Coarse boundary connectivity
    
    real(rp) :: area
    
    real(rp) :: coordsElem(size(coord,1),size(lnods,1)), areaElem
    integer(ip) :: numElemNodes, ielem
    
    print*,"size(lnods,2): ===>>> ", size(lnods,2)
    
    area = 0.0
    do ielem=1,size(lnods,2)
      numElemNodes = size(lnods,1)
      print*,"coordsElem"
      print*,"lnods(:,ielem)  ",lnods(:,ielem)
      print*,"sizecoord: ",size(coord,1)," ", size(coord,2)
      coordsElem = coord(:, lnods(:,ielem) )
      print*,"areaElem"
      areaElem = computeArea( numElemNodes , coordsElem )
      print*,"area"
      area = area + areaElem
    end do
    
    return
  end function computeAreaMesh_noReference
  !
  !
  !
  function projectToTangent(pk,tangent1,tangent2) result(newPk)
    !***
    !*  Projects the vector pk to the tangent vectors tangent1, tangent2 of the pseudoNormal plane
    !*** 
    implicit none

    real(rp), intent(in):: pk(3), tangent1(3), tangent2(3)
    real(rp)      :: newPk(3)

    newPk(1) = dot_product(pk,tangent1)
    newPk(2) = dot_product(pk,tangent2)
    newPk(3) = 0.0_rp
    
    !   newPk = dot_product(pk,tangent1)*tangent1 + dot_product(pk,tangent2)*tangent2
    return
  end function projectToTangent
  !
  !
  !
  subroutine checkClustering(countEmptyCells,mnodb,lnodb,ltypb)
    ! check if there is any voronoi cluster that has two or more disconnected parts
    ! if so, assign to one of the parts a zero voronoi cluster
    implicit none
    integer(ip),  intent(out):: countEmptyCells
    integer(ip),  intent(in)   :: mnodb                 !< Max nb of nodes per boundary element
    integer(ip),  intent(in)   :: lnodb(mnodb,numElems)    !< Boundary connectivity
    integer(ip),  intent(in)  :: ltypb(numElems)          !< Boundary type


    !***internal variables
    integer(ip)  :: ielem, ivor, iedge, numElemNodes, aux_i,edgeId
    integer(ip)  :: visited!, countVor
    integer(ip)  :: visitedElems(numElems),visitedList(numElems)
    integer(ip),allocatable  :: elemVorConn(:)  !maxNumElemVor,nVor


    countEmptyCells = 0_ip

    if(.not.allocated(isEliminated)) then
       allocate(isEliminated(nVor))
       isEliminated = 0
    end if
    allocate( elemVorConn(nVor) )
    do ielem = 1,numElems
       ivor = elemVor(ielem)
       !     if(elemVorConn(ivor)==0) then
       elemVorConn(ivor) = ielem
       !       countVor = countVor+1
       !     end if
       !     if(countVor==nVor) exit elemLoop
    end do

    do ivor = 1,nVor
       if( numCellsVor(ivor)>0 ) then
          !isEliminated(ivor) = 0

          visited = 0
          visitedElems = 0
          visitedList  = 0 

          ielem = elemVorConn(ivor)!1,ivor)

          call elementCheck(visited, ielem,ivor,visitedElems, visitedList)

          if(visited < numCellsVor(ivor)) then
             countEmptyCells = countEmptyCells + visited
             numCellsVor(ivor) = numCellsVor(ivor) - visited
             do aux_i = 1,visited
                ielem = visitedList(aux_i)
                elemVor(ielem) = 0
                !voronoi(ielem,ivor) = 0

                numElemNodes = element_type(ltypb(ielem)) % number_nodes
                do iedge = 1,numElemNodes
                   edgeId = elementEdges(iedge,ielem)
                   if(edgeVor(1,edgeId)==1) then!edgeIsBound(edgeId)
                      if(edgeVor(2,edgeId)==ivor) then
                         edgeVor(2,edgeId)=0
                      else
                         edgeVor(3,edgeId)=0
                      end if
                   end if
                end do
             end do
          end if
       else
          isEliminated(ivor) = 1
       end if
    end do

    deallocate(elemVorConn)

    return
  end subroutine checkClustering
  !
  !
  !
  recursive subroutine elementCheck(value,ielem,ivor,visitedElems, visitedList)

    implicit none
    integer(ip), intent(in)  :: ivor, ielem
    integer(ip) ,intent(inout)   :: visitedElems(numElems)
    integer(ip) ,intent(inout)   :: visitedList(numElems)
    integer(ip) ,intent(inout)   :: value

    integer(ip):: maxElemEdges,iedge, edgeId
    maxElemEdges = 4

    !   if( (visitedElems(ielem)==0).and.(voronoi(ielem,ivor)==1)) then
    !if( visitedElems(ielem)==0 ) then
    value = value+1
    visitedElems(ielem) = 1
    visitedList(value) = ielem
    do iedge =1 ,maxElemEdges
       edgeId = elementEdges(iedge,ielem)
       if((edgeId>0).and.(edgeInfo(2,edgeId)>0)) then
          if((edgeInfo(1,edgeId).ne.ielem).and.(visitedElems(edgeInfo(1,edgeId))==0).and.(elemVor(edgeInfo(1,edgeId))==ivor) ) then
             call elementCheck( value, edgeInfo(1,edgeId),ivor,visitedElems,visitedList)
          else if(  (visitedElems(edgeInfo(2,edgeId))==0).and.(elemVor(edgeInfo(2,edgeId))==ivor) ) then
             call elementCheck( value, edgeInfo(2,edgeId),ivor,visitedElems,visitedList)
          end if
       end if
    end do
    !end if  

    return
  end subroutine elementCheck
  !
  !
  !
  subroutine updateVariablesVoronoi(signi,ivor,ielem)
    implicit none

    integer(ip),  intent(in):: signi ! +-1 (depending if we add or remove the cell)
    integer(ip),  intent(in):: ivor 
    integer(ip),  intent(in):: ielem

    SRho(ivor) = SRho(ivor) + real(signi,rp)*areaElem(ielem)
    SGamma(:,ivor) = SGamma(:,ivor) + real(signi,rp)*areaElem(ielem)*cmElem(:,ielem)

    return
  end subroutine updateVariablesVoronoi
  !
  !
  !
  function goalF(sgamma1,srho1,sgamma2,srho2) result(F)

    implicit none
    real(rp),  intent(in):: sgamma1(numCoords),sgamma2(numCoords)
    real(rp),  intent(in):: srho1,srho2

    real(rp) :: F

    if((srho1==0).or.(srho2==0)) then
       F = HUGE(1.0_rp)
       return 
    end if

    F = ( computeVectorNorm(sgamma1)**2 ) / srho1 +&
         ( computeVectorNorm(sgamma2)**2 ) / srho2

    return
  end function goalF
  !
  !
  !
  function goalFGlobal(ivor,sgammaVor,srhoVor) result(F)

    implicit none

    integer(ip)  , intent(in)  :: ivor
    real(rp)    , intent(in)  :: sgammaVor(numCoords),srhoVor
    real(rp)           :: F

    real(rp)           ::gammaRho(numCoords)
    integer(ip)        ::ielem

    gammaRho  = sgammaVor/srhoVor

    F = 0.0_rp
    do ielem=1,numElems
       if(elemVor(ielem)==ivor) then
          !F = F + areaElem(ielem)*computeVectorNorm(  cmElem(:,ielem) - SGamma(:,ivor)/SRho(ivor)  )**2
          F = F + areaElem(ielem)*computeVectorNorm(  cmElem(:,ielem) - gammaRho  )**2
       end if
    end do

    return
  end function goalFGlobal
  !
  !
  ! 
  recursive function computeArea(numElemNodes,coord) result(area)
    implicit none

    integer(ip),  intent(in):: numElemNodes
    real(rp),    intent(in):: coord(numCoords,numElemNodes)

    real(rp)  :: area
    real(rp)  :: v1(3),v2(3)

    v1 = 0.0_rp
    v2 = 0.0_rp

    if(numElemNodes==3) then ! compute area of a triangle
       v1(1:numCoords) = coord(:,2) - coord(:,1)
       v2(1:numCoords) = coord(:,3) - coord(:,1)
       area = computeVectorNorm( crossProduct(v1,v2) ) / 2.0_rp
    elseif(numElemNodes==4) then ! compute area as the sum of two triangles
       area = 0.0_rp
       area = area + computeArea(3_ip,coord(:,(/1,2,3/)) )
       area = area + computeArea(3_ip,coord(:,(/1,3,4/)) )
    else 
       call runend('Only implemented for tris and quad boundary meshes')
    end if
    return
  end function computeArea
  !
  !
  ! 
  recursive function computeArea_withSign(numElemNodes,coord,n) result(area)
    implicit none

    integer(ip),  intent(in):: numElemNodes
    real(rp),    intent(in):: coord(numCoords,numElemNodes),n(3)

    real(rp)  :: area
    real(rp)  :: v1(3),v2(3),M(3,3),v3(3)

    v1 = 0.0_rp
    v2 = 0.0_rp

    if(numElemNodes==3) then ! compute area of a triangle
       v1(1:numCoords) = coord(:,2) - coord(:,1)
       v2(1:numCoords) = coord(:,3) - coord(:,1)
       v3(:) = crossProduct(v1,v2)
       area = computeVectorNorm( v3 ) / 2.0_rp
       
       M(:,1)=v1
       M(:,2)=v2
       M(:,3)=n!v3
       area=sign(area,determinant3D(M))
       
       !print*,area,"   ",determinant3D(M)
       
    else
       call runend('Only implemented for tri ')
    end if
    
!     if( area<0.0) then
!       call runend("epa, area negativa")
!     end if
    
    return
  end function computeArea_withSign
  !
  !
  !
  recursive function computeNormal(numElemNodes,coord) result(n)
    implicit none

    integer(ip),  intent(in):: numElemNodes
    real(rp),    intent(in):: coord(numCoords,numElemNodes)

    real(rp)  :: n(3)
    real(rp)  :: v1(3),v2(3)

    v1 = 0.0_rp
    v2 = 0.0_rp

    if(numElemNodes==3) then ! compute area of a triangle
       v1(1:numCoords) = coord(:,2) - coord(:,1)
       v2(1:numCoords) = coord(:,3) - coord(:,1)
       n = crossProduct(v1,v2)
       n = n/computeVectorNorm(n)
    else
       call runend('Only implemented for tri ')
    end if
    
    return
  end function computeNormal
  !
  !
  !
  function determinant3D(S) result(detS)
  	implicit none 

  	real(rp), intent(in) :: S(3,3)
  	real(rp) :: detS

  	detS = S(1,1)*S(2,2)*S(3,3)  &
      - S(1,1)*S(2,3)*S(3,2)  &
      - S(1,2)*S(2,1)*S(3,3)  &
      + S(1,2)*S(2,3)*S(3,1)  &
      + S(1,3)*S(2,1)*S(3,2)  &
      - S(1,3)*S(2,2)*S(3,1)
        return
  end function determinant3D
  !
  !
  !
  function crossProduct(a, b) result(v)
    implicit none
    real(rp), intent(in) :: a(3)
    real(rp), intent(in) :: b(3)
    real(rp)             :: v(3)

    v(1) = a(2) * b(3) - a(3) * b(2)
    v(2) = a(3) * b(1) - a(1) * b(3)
    v(3) = a(1) * b(2) - a(2) * b(1)
    return
  end function crossProduct
  !
  !
  !
  function computeVectorNorm(v) result(normV)
    implicit none
    real(rp), intent(in) :: v(numCoords)
    real(rp)             :: normV
    if(numCoords==3) then
       normV       = sqrt(v(1)*v(1) + v(2)*v(2)+ v(3)*v(3))
    else
       normV       = sqrt(v(1)*v(1) + v(2)*v(2))
    end if
    return
  end function computeVectorNorm
  !
  !
  !
  subroutine deallocateModuleVariables
    if( allocated(elemVor)           ) deallocate(elemVor)
    if( allocated(globToLoc)         ) deallocate(globToLoc)
    if( allocated(elementEdges)      ) deallocate(elementEdges)
    if( allocated(locToGlob)         ) deallocate(locToGlob)
    if( allocated(countNumEdgesNode) ) deallocate(countNumEdgesNode)
    if( allocated(nodeEdges)         ) deallocate(nodeEdges)
    if( allocated(areaElem)          ) deallocate(areaElem)
    if( allocated(cmElem)            ) deallocate(cmElem)
    if( allocated(edgeInfo)          ) deallocate(edgeInfo)
    if( allocated(numCellsVor)       ) deallocate(numCellsVor)
    if( allocated(SGamma)            ) deallocate(SGamma)
    if( allocated(SRho)              ) deallocate(SRho)
    if( allocated(edgeVor)           ) deallocate(edgeVor)
    if( allocated(edgeIsBound)       ) deallocate(edgeIsBound)
    if( allocated(nodeVor)           ) deallocate(nodeVor)
    if( allocated(numNodeVor)        ) deallocate(numNodeVor)
    if( allocated(centroid)          ) deallocate(centroid)
    if( allocated(normalToVor)       ) deallocate(normalToVor) 
    if( allocated(isEliminated)      ) deallocate(isEliminated)
    if( allocated(mapValidVor)       ) deallocate(mapValidVor)
    if( allocated(coord_coarse_ref)  ) deallocate(coord_coarse_ref)
    if( allocated(coord_coarse_closest)  ) deallocate(coord_coarse_closest)
    if( allocated(coord_coarse_furthest)  ) deallocate(coord_coarse_furthest)
    if( allocated(coarseNodesToFine)) deallocate(coarseNodesToFine)
  end subroutine deallocateModuleVariables

end module mod_surface_mesh
!> @}


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
