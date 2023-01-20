!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @defgroup Meshes_Toolbox
!> Toolbox for meshes manipulations, like creating submeshes,
!> surface mesh, etc.
!> @{
!> @name    ToolBox for meshes operations
!> @file    mod_meshes.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for meshes operations
!> @details ToolBox for meshes operations
!
!-----------------------------------------------------------------------

module mod_meshes

  use def_kintyp,            only : ip,rp,lg,i1p,i1pp,r3p
  use def_elmtyp,            only : ELFEM
  use def_elmtyp,            only : BAR02,TRI03
  use def_master,            only : INOTMASTER,kfl_paral
  use def_master,            only : intost
  use def_master,            only : NELEM_TYPE 
  use def_master,            only : NPOIN_TYPE 
  use def_domain,            only : memor_dom
  use def_domain,            only : mesh_type
  use def_domain,            only : mnode
  use mod_graphs,            only : graphs_elepoi
  use mod_graphs,            only : graphs_dealep
  use mod_graphs,            only : graphs_elepoi_deallocate
  use mod_maths,             only : maths_heap_sort
  use mod_elmgeo,            only : element_type
  use mod_memory,            only : memory_alloca 
  use mod_memory,            only : memory_deallo
  use mod_memory,            only : memory_copy
  use mod_memory,            only : memory_size
  use mod_communications,    only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications,    only: PAR_MAX
  use mod_optional_argument, only : optional_argument
  use mod_strings,           only : integer_to_string
  use mod_htable
  use mod_std
  implicit none
  private
  real(rp)                  :: epsil = epsilon(1.0_rp)
  character(200), parameter :: vacal = 'mod_meshes'
  
  interface meshes_submesh
     module procedure meshes_submesh_all,&
          &           meshes_submesh_s
  end interface meshes_submesh

  public :: meshes_submesh
  public :: meshes_surface_from_nodal_array
  public :: meshes_surface_from_nodal_array_deallocate
  public :: meshes_list_boundary_nodes
  public :: meshes_list_boundary_elements
  public :: meshes_check_mesh
  public :: meshes_gather_submesh
  public :: meshes_glue_two_meshes
  
contains 

  !-----------------------------------------------------------------------
  !
  !> @brief   Compute a submesh
  !> @details Compute a submesh given a permutation array
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

    subroutine meshes_submesh_s(&
       & mesh_in,mesh_sub,lperm,lnper,MESH_NAME)

    type(mesh_type),                      intent(in)    :: mesh_in
    type(mesh_type),                      intent(inout) :: mesh_sub
    integer(ip),       pointer,           intent(inout) :: lperm(:)
    integer(ip),       pointer, optional, intent(inout) :: lnper(:)
    character(len=*),           optional, intent(in)    :: MESH_NAME
    integer(ip)                                         :: ielem,inode,ipoin
    integer(ip)                                         :: pelty,pnode,mnode_loc
    integer(ip)                                         :: kelem,ipoin_sub
    integer(ip), pointer                                :: lpoin(:)
    character(20)                                       :: my_mesh_name
    
    if( present(MESH_NAME) ) then
       my_mesh_name = trim(MESH_NAME)
    else
       my_mesh_name = 'MESH'
    end if

    mnode_loc = size(mesh_in % lnods,1,KIND=ip)
    mesh_sub % npoin = 0
    if( .not. associated(lperm) ) then
       mesh_sub % nelem = 0
    else
       mesh_sub % nelem = size(lperm,KIND=ip)
    end if
    mesh_sub % ndime = mesh_in % ndime
    mesh_sub % mnode = mnode_loc
    
    if( mesh_sub % nelem /= 0 ) then

       call memory_alloca(memor_dom,trim(my_mesh_name)//' % LNODS','meshes_submesh',mesh_sub % lnods,mnode_loc,mesh_sub % nelem)
       call memory_alloca(memor_dom,trim(my_mesh_name)//' % LTYPE','meshes_submesh',mesh_sub % ltype,          mesh_sub % nelem)
       call memory_alloca(memor_dom,trim(my_mesh_name)//' % LNNOD','meshes_submesh',mesh_sub % lnnod,          mesh_sub % nelem)
       call memory_alloca(memor_dom,trim(my_mesh_name)//' % LEINV','meshes_submesh',mesh_sub % leinv_loc,      mesh_sub % nelem)

       allocate(lpoin(mesh_in % npoin))
       do ipoin = 1,mesh_in % npoin
          lpoin(ipoin) = 0
       end do
       !
       ! Renumber nodes
       !
       do kelem = 1,mesh_sub % nelem
          ielem                       = lperm(kelem)
          pelty                       = abs(mesh_in % ltype(ielem))
          pnode                       = mesh_in % lnnod(ielem)
          mesh_sub % ltype(kelem)     = pelty
          mesh_sub % lnnod(kelem)     = pnode
          mesh_sub % leinv_loc(kelem) = mesh_in % leinv_loc(ielem)
          do inode = 1,pnode
             ipoin = mesh_in % lnods(inode,ielem)
             if( lpoin(ipoin) == 0 ) then
                mesh_sub % npoin = mesh_sub % npoin + 1
                lpoin(ipoin)     = mesh_sub % npoin
             end if
          end do
       end do
       !
       ! Copy connectivity and coordinates
       !
       call memory_alloca(memor_dom,trim(my_mesh_name)//' % COORD','meshes_submesh',mesh_sub % coord,mesh_sub % ndime,mesh_sub % npoin)
       call memory_alloca(memor_dom,trim(my_mesh_name)//' % LNINV','meshes_submesh',mesh_sub % lninv_loc,mesh_sub % npoin)

       do kelem = 1,mesh_sub % nelem
          ielem = lperm(kelem)
          pelty = abs(mesh_in % ltype(ielem))
          pnode = mesh_in % lnnod(ielem)
          do inode = 1,pnode
             ipoin = mesh_in % lnods(inode,ielem)
             mesh_sub % lnods(inode,kelem) = lpoin(ipoin)
          end do
       end do
       do ipoin = 1,mesh_in % npoin
          ipoin_sub = lpoin(ipoin)
          if( ipoin_sub /= 0 ) then
             mesh_sub % lninv_loc(ipoin_sub) = mesh_in % lninv_loc(ipoin) 
             mesh_sub % coord(1:mesh_in % ndime,ipoin_sub) = mesh_in % coord(1:mesh_in % ndime,ipoin)
             if( present(lnper) ) lnper(ipoin_sub) = ipoin
          end if
       end do

       deallocate(lpoin)

    else

       mesh_sub % npoin = 0
       mesh_sub % nelem = 0

    end if

  end subroutine meshes_submesh_s

  subroutine meshes_submesh_all(                                                    &
       & ndime,    npoin    ,nelem    ,coord    ,lninv,lnods    ,ltype    ,lnnod,     leinv,&
       &           npoin_sub,nelem_sub,coord_sub,lninv_sub,lnods_sub,ltype_sub,lnnod_sub, leinv_sub,&
       &           lperm   ,lnper)

    integer(ip),          intent(in)              :: ndime
    integer(ip),          intent(in)              :: nelem
    integer(ip),          intent(in)              :: npoin
    real(rp),    pointer, intent(in)              :: coord(:,:)
    integer(ip), pointer, intent(in)              :: lninv(:)
    integer(ip), pointer, intent(in)              :: lnods(:,:)
    integer(ip), pointer, intent(in)              :: ltype(:)
    integer(ip), pointer, intent(in)              :: lnnod(:)
    integer(ip), pointer, intent(in)              :: leinv(:)
    integer(ip),          intent(out)             :: nelem_sub
    integer(ip),          intent(out)             :: npoin_sub
    real(rp),    pointer, intent(inout)           :: coord_sub(:,:)
    integer(ip), pointer, intent(inout)           :: lninv_sub(:)
    integer(ip), pointer, intent(inout)           :: lnods_sub(:,:)
    integer(ip), pointer, intent(inout)           :: ltype_sub(:)
    integer(ip), pointer, intent(inout)           :: lnnod_sub(:)
    integer(ip), pointer, intent(inout)           :: leinv_sub(:)
    integer(ip), pointer, intent(inout)           :: lperm(:)
    integer(ip), pointer, intent(inout), optional :: lnper(:)
    integer(ip)                                   :: ielem,inode,ipoin
    integer(ip)                                   :: pelty,pnode,mnode_loc
    integer(ip)                                   :: kelem,ipoin_sub
    integer(ip), pointer                          :: lpoin(:)

    mnode_loc = size(lnods,1,KIND=ip)
    npoin_sub = 0
    if( .not. associated(lperm) ) then
       nelem_sub = 0
    else
       nelem_sub = size(lperm,KIND=ip)
    end if

    if( nelem_sub /= 0 ) then

       call memory_alloca(memor_dom,'LNODS_SUB','meshes_submesh',lnods_sub,mnode_loc,nelem_sub)
       call memory_alloca(memor_dom,'LTYPE_SUB','meshes_submesh',ltype_sub,          nelem_sub)
       call memory_alloca(memor_dom,'LNNOD_SUB','meshes_submesh',lnnod_sub,          nelem_sub)
       call memory_alloca(memor_dom,'LEINV_SUB','meshes_submesh',leinv_sub,          nelem_sub)

       allocate(lpoin(npoin))
       do ipoin = 1,npoin
          lpoin(ipoin) = 0
       end do
       !
       ! Renumber nodes
       !
       do kelem = 1,nelem_sub
          ielem = lperm(kelem)
          pelty = abs(ltype(ielem))
          pnode = lnnod(ielem)
          ltype_sub(kelem) = pelty
          lnnod_sub(kelem) = pnode
          leinv_sub(kelem) = leinv(ielem)
          do inode = 1,pnode
             ipoin = lnods(inode,ielem)
             if( lpoin(ipoin) == 0 ) then
                npoin_sub = npoin_sub + 1
                lpoin(ipoin) = npoin_sub
             end if
          end do
       end do
       !
       ! Copy connectivity and coordinates
       !
       call memory_alloca(memor_dom,'COORD_SUB','meshes_submesh',coord_sub,ndime,npoin_sub)
       call memory_alloca(memor_dom,'LNINV_SUB','meshes_submesh',lninv_sub,npoin_sub)

       do kelem = 1,nelem_sub
          ielem = lperm(kelem)
          pelty = abs(ltype(ielem))
          pnode = lnnod(ielem)
          do inode = 1,pnode
             ipoin = lnods(inode,ielem)
             lnods_sub(inode,kelem) = lpoin(ipoin)
          end do
       end do
       do ipoin = 1,npoin
          ipoin_sub = lpoin(ipoin)
          if( ipoin_sub /= 0 ) then
             lninv_sub(ipoin_sub) = lninv(ipoin) 
             coord_sub(1:ndime,ipoin_sub) = coord(1:ndime,ipoin)
             if( present(lnper) ) lnper(ipoin_sub) = ipoin
          end if
       end do

       deallocate(lpoin)

    else

       npoin_sub = 0
       nelem_sub = 0

    end if

  end subroutine meshes_submesh_all

  !-----------------------------------------------------------------------
  !
  !> @brief   Construct a surface mesh
  !> @details Construct a surface mesh from a nodal array xarra
  !> @date    29/09/2015
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine meshes_surface_from_nodal_array_deallocate(&
       lnodb_sur,coord_sur,ltypb_sur,lelem_sur)
    implicit none
    integer(ip),     pointer,           intent(inout) :: lnodb_sur(:,:)
    real(rp),        pointer,           intent(inout) :: coord_sur(:,:)
    integer(ip),     pointer, optional, intent(inout) :: ltypb_sur(:)
    integer(ip),     pointer, optional, intent(inout) :: lelem_sur(:)

    if( INOTMASTER ) then
       call memory_deallo(memor_dom,'LNODB_SUR','meshes_surface_from_nodal_array',lnodb_sur)
       call memory_deallo(memor_dom,'COORD_SUR','meshes_surface_from_nodal_array',coord_sur)
       if( present(ltypb_sur) ) call memory_deallo(memor_dom,'LTYPB_SUR','meshes_surface_from_nodal_array',ltypb_sur)
       if( present(lelem_sur) ) call memory_deallo(memor_dom,'LELEM_SUR','meshes_surface_from_nodal_array',lelem_sur)
    end if

  end subroutine meshes_surface_from_nodal_array_deallocate

  !-----------------------------------------------------------------------
  !
  !> @brief   Construct a surface mesh
  !> @details Construct a surface mesh from a nodal array xarra
  !> @date    29/09/2015
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine meshes_surface_from_nodal_array(&
       xarra,meshe,npoin_sur,nboun_sur,lnodb_sur,coord_sur,ltypb_sur,lelem_sur)
    implicit none
    real(rp),        pointer,           intent(in)  :: xarra(:,:)                 !< 
    type(mesh_type),                    intent(in)  :: meshe                      !< Mesh type
    integer(ip),                        intent(out) :: npoin_sur
    integer(ip),                        intent(out) :: nboun_sur
    integer(ip),     pointer,           intent(inout) :: lnodb_sur(:,:)
    real(rp),        pointer,           intent(inout) :: coord_sur(:,:)
    integer(ip),     pointer, optional, intent(inout) :: ltypb_sur(:)
    integer(ip),     pointer, optional, intent(inout) :: lelem_sur(:)
    integer(ip)                                     :: ielem                      ! Indices and dimensions
    integer(ip)                                     :: pnode
    integer(ip)                                     :: inode,jnode,ipoin
    integer(ip)                                     :: compt,compli,compg
    integer(ip)                                     :: pnodb
    real(rp)                                        :: elcod(3,mnode)
    real(rp)                                        :: elarr(mnode)
    real(rp)                                        :: l1,lp,x1,y1,x2,y2
    real(rp)                                        :: norma(3),vecto(3)
    !
    ! Mesh arrays required
    !
    integer(ip)                                     :: nelem
    integer(ip)                                     :: npoin
    integer(ip)                                     :: ndime
    integer(ip),      pointer                       :: lnods(:,:)
    integer(ip),      pointer                       :: lnnod(:)
    integer(ip),      pointer                       :: lelch(:)
    real(rp),         pointer                       :: coord(:,:)

    if( INOTMASTER ) then
       !
       ! Mesh arrays
       !
       nelem =  meshe % nelem 
       npoin =  meshe % npoin 
       ndime =  meshe % ndime 
       lnods => meshe % lnods
       lnnod => meshe % lnnod
       lelch => meshe % lelch
       coord => meshe % coord
       compt =  0_ip

       if( ndime == 2 ) then
          !
          ! Loop over the elements to construct discrete interface 
          ! collection of segments (2d case)
          !
          pnodb = 2
          do ielem = 1,nelem
             if( lelch(ielem) == ELFEM ) then
                pnode = lnnod(ielem)
                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   elarr(inode) = xarra(ipoin,1)
                   elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                end do

                compli = 0_ip
                do inode = 1,pnode
                   if( inode == pnode ) then
                      jnode = 1 
                   else
                      jnode = inode + 1
                   endif
                   !
                   ! Count interface points 
                   !
                   if( elarr(inode)*elarr(jnode) < 0.0_rp ) then
                      compli = compli + 1
                      if( compli == 3 ) compli = 1         ! Count interface points 
                      if( compli == 1 ) compt = compt + 1 ! Count interface segments
                   endif
                end do
             end if
          end do
          !
          ! Registering of discrete interface vectors size
          !
          nboun_sur = compt
          npoin_sur = ndime * nboun_sur
          compt     = 0_ip
          compg     = 0_ip
          !
          ! Memory allocation for discrete interface vectors
          ! connectivity (lnodb_sur) and points coordinates (coord_sur)
          !
          call memory_deallo(memor_dom,'LNODB_SUR','meshes_surface_from_nodal_array',lnodb_sur)
          call memory_deallo(memor_dom,'COORD_SUR','meshes_surface_from_nodal_array',coord_sur)
          if( present(ltypb_sur) ) call memory_deallo(memor_dom,'LTYPB_SUR','meshes_surface_from_nodal_array',ltypb_sur)
          if( present(lelem_sur) ) call memory_deallo(memor_dom,'LELEM_SUR','meshes_surface_from_nodal_array',lelem_sur)

          call memory_alloca(memor_dom,'LNODB_SUR','meshes_surface_from_nodal_array',lnodb_sur,pnodb,nboun_sur)
          call memory_alloca(memor_dom,'COORD_SUR','meshes_surface_from_nodal_array',coord_sur,ndime,npoin_sur)
          if( present(ltypb_sur) ) call memory_alloca(memor_dom,'LTYPB_SUR','meshes_surface_from_nodal_array',ltypb_sur,nboun_sur)
          if( present(lelem_sur) ) call memory_alloca(memor_dom,'LELEM_SUR','meshes_surface_from_nodal_array',lelem_sur,nboun_sur)
          !
          ! filling of discrete interface vectors 
          ! connectivity (lnodb_sur) and points coordinates (coord_sur)
          !
          do ielem=1,nelem
             if( lelch(ielem) == ELFEM ) then
                pnode = lnnod(ielem)
                do inode = 1,pnode
                   ipoin                = lnods(inode,ielem)
                   elarr(inode)         = xarra(ipoin,1)
                   elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                end do

                compli = 0_ip

                do inode = 1,pnode
                   if( inode == pnode ) then
                      jnode = 1 
                   else
                      jnode = inode + 1
                   endif
                   if( elarr(inode)*elarr(jnode) < 0.0_rp ) then
                      compg = compg + 1
                      compli = compli + 1
                      if( compli == 3 ) compli = 1
                      if( compli == 1 ) compt = compt + 1
                      !
                      ! Compute the intersection of the elements with the surface 
                      !
                      l1 = abs(elarr(inode))
                      lp = abs(elarr(inode)-elarr(jnode))
                      x1 = elcod(1,inode)
                      x2 = elcod(1,jnode)
                      y1 = elcod(2,inode)
                      y2 = elcod(2,jnode)
                      lnodb_sur(compli,compt) = compg
                      coord_sur(1,compg)     =  x1 * (1.0_rp-l1/lp) + x2 * l1 / lp  
                      coord_sur(2,compg)     =  y1 * (1.0_rp-l1/lp) + y2 * l1 / lp
                      if( present(ltypb_sur) ) ltypb_sur(compt) = BAR02
                      if( present(lelem_sur) ) lelem_sur(compt) = ielem
                   end if
                end do
                !
                ! Orientate surface
                !
                if( compli == 2 ) then
                   x1       = coord_sur(1,compg-1)
                   y1       = coord_sur(2,compg-1)
                   x2       = coord_sur(1,compg)
                   y2       = coord_sur(2,compg)
                   norma(1) = y1-y2
                   norma(2) = x2-x1            
                   inode    = 1
                   vecto(1) = (elcod(1,inode)-0.5_rp*(x1+x2))
                   vecto(2) = (elcod(2,inode)-0.5_rp*(y1+y2))
                   l1       = norma(1)*vecto(1)+norma(2)*vecto(2)
                   if( l1 > 0.0_rp ) then
                      if( elarr(inode) >= 0.0_rp ) then
                         lnodb_sur(1,compt) = compg
                         lnodb_sur(2,compt) = compg-1
                      end if
                   else                      
                      if( elarr(inode) < 0.0_rp ) then
                         lnodb_sur(1,compt) = compg
                         lnodb_sur(2,compt) = compg-1
                      end if
                   end if
                end if

             end if
          end do

       else

          call runend('NOT CODED')

       end if

    end if

  end subroutine meshes_surface_from_nodal_array
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-31
  !> @brief   List of boundary nodes
  !> @details List of boundary noed of a mesh. No parallel exchange!
  !> 
  !-----------------------------------------------------------------------

  subroutine meshes_list_boundary_nodes(&
       nelem,npoin,lnods,ltype,&
       number_boundary_nodes,&
       list_boundary_nodes,&
       LIST_BOUNDARY_NODE_NAME)

    integer(ip),          intent(in)             :: nelem
    integer(ip),          intent(in)             :: npoin
    integer(ip), pointer, intent(in)             :: lnods(:,:)
    integer(ip), pointer, intent(in)             :: ltype(:)
    integer(ip),          intent(out)            :: number_boundary_nodes
    integer(ip), pointer, intent(inout)          :: list_boundary_nodes(:)
    character(len=*),     intent(in),   optional :: LIST_BOUNDARY_NODE_NAME
    integer(ip)                                  :: ielty,ielem,iface,inodf
    integer(ip)                                  :: inode,jelem,jface,jelty,ipoin,pnodf
    integer(ip)                                  :: ielpo,mface,mnodf,ii,mnode_loc
    logical(lg)                                  :: equal_faces  

    integer(ip)                                  :: mepoi_loc
    integer(ip), pointer                         :: lelpo_loc(:)
    integer(ip), pointer                         :: pelpo_loc(:)
    integer(ip), pointer                         :: lnnod_loc(:)

    integer(ip)                                  :: nfacg
    integer(ip), pointer                         :: facel(:,:,:)
    type(i1p),   pointer                         :: lelfa(:)

    logical(lg), pointer                         :: boundary_nodes(:)

    nullify(lelpo_loc)
    nullify(pelpo_loc)
    nullify(lnnod_loc)

    nullify(facel)
    nullify(lelfa)

    nullify(boundary_nodes)

    !----------------------------------------------------------------------
    !
    ! LELPO_LOC,PELPO_LOC: node-to-element graph
    !
    !----------------------------------------------------------------------

    if( INOTMASTER .and. nelem > 0 ) then

       call memory_alloca(memor_dom,'LNNOD_LOC','meshes_list_boundary_nodes',lnnod_loc,nelem)
       mnode_loc = size(lnods,1,KIND=ip)
       mface = 0
       mnodf = 0
       do ielem = 1,nelem
          ielty            = abs(ltype(ielem))
          mface            = max(mface,element_type(ielty) % number_faces)
          mnodf            = max(mnodf,element_type(ielty) % max_face_nodes)
          lnnod_loc(ielem) = element_type(ielty) % number_nodes       
       end do
       call graphs_elepoi(&
            npoin,nelem,mnode_loc,lnods,lnnod_loc,mepoi_loc,pelpo_loc,lelpo_loc,&
            PELPO_NAME='PELPO_LOC',LELPO_NAME='LELPO_LOC',memor=memor_dom)
       call memory_deallo(memor_dom,'LNNOD_LOC','meshes_list_boundary_nodes',lnnod_loc)

       !----------------------------------------------------------------------
       !
       ! LELFA: List of global faces
       !
       !----------------------------------------------------------------------
       !
       ! Allocate memory for lelfa, FACES, CFAEL AND NNODG
       !
       call memory_alloca(memor_dom,'LELFA','meshes_list_boundary_nodes',lelfa,nelem)
       call memory_alloca(memor_dom,'FACEL','meshes_list_boundary_nodes',facel,mnodf+1_ip,mface,nelem)
       !
       ! Construct and sort FACES
       !
       do ielem = 1,nelem                                          
          ielty = abs(ltype(ielem))
          call memory_alloca(memor_dom,'LELFA % L','meshes_list_boundary_nodes',lelfa(ielem)%l,element_type(ielty) % number_faces)
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
          ielty = abs(ltype(ielem))                                  ! eliminate the repited faces
          do iface = 1,element_type(ielty) % number_faces
             if( facel(mnodf+1,iface,ielem) > 0 ) then
                nfacg = nfacg + 1
                ipoin = facel(1,iface,ielem)
                ielpo = pelpo_loc(ipoin)-1
                do while( ielpo < pelpo_loc(ipoin+1)-1 )
                   ielpo = ielpo + 1
                   jelem = lelpo_loc(ielpo)
                   if( jelem /= ielem ) then
                      jelty = abs(ltype(jelem))                      ! eliminate the repited faces
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
                               ielpo                      =  pelpo_loc(ipoin+1)                 ! Exit JELEM do  
                            end if
                         end if
                      end do
                   end if
                end do
             end if
          end do
       end do

       call memory_alloca(memor_dom,'BOUNDARY_NODES','meshes_list_boundary_nodes',boundary_nodes,npoin)

       do ielem = 1,nelem                      
          ielty = abs(ltype(ielem))                    
          do iface = 1,element_type(ielty) % number_faces
             if( facel(mnodf+1,iface,ielem) > 0 .and. facel(1,iface,ielem) > 0 ) then
                pnodf = element_type(ielty) % node_faces(iface)
                do inodf = 1,pnodf
                   inode = element_type(ielty) % list_faces(inodf,iface) 
                   ipoin = lnods(inode,ielem)
                   boundary_nodes(ipoin) = .true.
                end do
             end if
          end do
       end do
       number_boundary_nodes = count(boundary_nodes,KIND=ip)
       if( associated(list_boundary_nodes) ) call runend('BOUNDARY NODE LIST ALREADY ASSOCIATED')
       if( present(LIST_BOUNDARY_NODE_NAME) ) then
          call memory_alloca(memor_dom,trim(LIST_BOUNDARY_NODE_NAME),'meshes_list_boundary_nodes',list_boundary_nodes,number_boundary_nodes)
       else          
          call memory_alloca(memor_dom,'LIST_BOUNDARY_NODES','meshes_list_boundary_nodes',list_boundary_nodes,number_boundary_nodes)
       end if
       if( number_boundary_nodes > 0 ) then
          ii = 0
          do ipoin = 1,npoin
             if( boundary_nodes(ipoin) ) then
                ii = ii + 1
                list_boundary_nodes(ii) = ipoin
             end if
          end do
       end if
       !
       ! Deallocate memory
       !
       call memory_deallo(memor_dom,'FACEL'         ,'meshes_list_boundary_nodes',facel)
       call memory_deallo(memor_dom,'LELFA'         ,'meshes_list_boundary_nodes',lelfa)
       call memory_deallo(memor_dom,'BOUNDARY_NODES','meshes_list_boundary_nodes',boundary_nodes)
       call graphs_elepoi_deallocate(&
            pelpo_loc,lelpo_loc,&
            PELPO_NAME='PELPO_LOC',LELPO_NAME='LELPO_LOC',memor=memor_dom)

    end if

  end subroutine meshes_list_boundary_nodes

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-31
  !> @brief   List of boundary nodes
  !> @details List of boundary noed of a mesh. No parallel exchange!
  !> 
  !-----------------------------------------------------------------------

  subroutine meshes_list_boundary_elements(&
       nelem,npoin,lnods,ltype,&
       nboun,lelbo,lnodb,ltypb,&
       MEMORY_COUNTER,&
       LELBO_NAME,&
       LNODB_NAME,&
       LTYPB_NAME)

    integer(ip),                         intent(in)    :: nelem
    integer(ip),                         intent(in)    :: npoin
    integer(ip),                pointer, intent(in)    :: lnods(:,:)
    integer(ip),                pointer, intent(in)    :: ltype(:)
    integer(ip),                         intent(out)   :: nboun
    integer(ip),                pointer, intent(inout) :: lelbo(:)
    integer(ip),      optional, pointer, intent(inout) :: lnodb(:,:)
    integer(ip),      optional, pointer, intent(inout) :: ltypb(:)
    integer(8),       optional,          intent(inout) :: MEMORY_COUNTER(2)
    character(len=*), optional,          intent(in)    :: LELBO_NAME
    character(len=*), optional,          intent(in)    :: LNODB_NAME
    character(len=*), optional,          intent(in)    :: LTYPB_NAME
    integer(ip)                                        :: ielty,ielem,iface,inodf
    integer(ip)                                        :: inode,jelem,jface,jelty,ipoin,pnodf
    integer(ip)                                        :: ielpo,mface,mnodf
    integer(ip)                                        :: iboun,mnode_loc
    integer(8)                                         :: memor_loc(2)
    logical(lg)                                        :: equal_faces  
    integer(ip)                                        :: mepoi_loc
    integer(ip),      pointer                          :: lelpo_loc(:)
    integer(ip),      pointer                          :: pelpo_loc(:)
    integer(ip),      pointer                          :: lnnod_loc(:)
    integer(ip)                                        :: nfacg
    integer(ip),      pointer                          :: facel(:,:,:)
    type(i1p),        pointer                          :: lelfa(:)
    character(50)                                      :: my_lelbo_name
    character(50)                                      :: my_lnodb_name
    character(50)                                      :: my_ltypb_name

    memor_loc     = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
    my_lelbo_name = optional_argument('BOUNDARY_ELEMENT'     ,LELBO_NAME)
    my_lnodb_name = optional_argument('BOUNDARY_CONNECTIVITY',LNODB_NAME)
    my_ltypb_name = optional_argument('BOUNDARY_TYPE',        LTYPB_NAME)

    nullify(lelpo_loc)
    nullify(pelpo_loc)
    nullify(lnnod_loc)

    nullify(facel)
    nullify(lelfa)

    !----------------------------------------------------------------------
    !
    ! LELPO_LOC,PELPO_LOC: node-to-element graph
    !
    !----------------------------------------------------------------------

    if( INOTMASTER ) then

       call memory_alloca(memor_loc,'LNNOD_LOC','meshes_list_boundary_elements',lnnod_loc,nelem)
       mnode_loc = size(lnods,1)
       mface = 0
       mnodf = 0
       do ielem = 1,nelem
          ielty            = abs(ltype(ielem))
          mface            = max(mface,element_type(ielty) % number_faces)
          mnodf            = max(mnodf,element_type(ielty) % max_face_nodes)
          lnnod_loc(ielem) = element_type(ielty) % number_nodes       
       end do

       call graphs_elepoi(&
            npoin,nelem,mnode_loc,lnods,lnnod_loc,mepoi_loc,pelpo_loc,lelpo_loc,&
            PELPO_NAME='PELPO_LOC',LELPO_NAME='LELPO_LOC',memor=memor_dom)

       call memory_deallo(memor_loc,'LNNOD_LOC','meshes_list_boundary_elements',lnnod_loc)

       !----------------------------------------------------------------------
       !
       ! LELFA: List of global faces
       !
       !----------------------------------------------------------------------
       !
       ! Allocate memory for lelfa, FACES, CFAEL AND NNODG
       !
       call memory_alloca(memor_loc,'LELFA','meshes_list_boundary_elements',lelfa,nelem)
       call memory_alloca(memor_loc,'FACEL','meshes_list_boundary_elements',facel,mnodf+1_ip,mface,nelem)
       !
       ! Construct and sort FACES
       !
       do ielem = 1,nelem                                          
          ielty = abs(ltype(ielem))
          call memory_alloca(memor_loc,'LELFA % L','meshes_list_boundary_elements',lelfa(ielem)%l,element_type(ielty) % number_faces)
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
          ielty = abs(ltype(ielem))                                  ! eliminate the repited faces
          do iface = 1,element_type(ielty) % number_faces
             if( facel(mnodf+1,iface,ielem) > 0 ) then
                nfacg = nfacg + 1
                ipoin = facel(1,iface,ielem)
                ielpo = pelpo_loc(ipoin)-1
                do while( ielpo < pelpo_loc(ipoin+1)-1 )
                   ielpo = ielpo + 1
                   jelem = lelpo_loc(ielpo)
                   if( jelem /= ielem ) then
                      jelty = abs(ltype(jelem))                      ! eliminate the repited faces
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
                               ielpo                      =  pelpo_loc(ipoin+1)                 ! Exit JELEM do  
                            end if
                         end if
                      end do
                   end if
                end do
             end if
          end do
       end do

       nboun = 0
       do ielem = 1,nelem                      
          ielty = abs(ltype(ielem))                    
          do iface = 1,element_type(ielty) % number_faces
             if( facel(mnodf+1,iface,ielem) > 0 .and. facel(1,iface,ielem) > 0 ) then
                nboun = nboun + 1
             end if
          end do
       end do
       !
       ! Lis of boundary elements (LELBO)
       !
       if( associated(lelbo) ) call runend('BOUNDARY NODE LIST ALREADY ASSOCIATED')
       call memory_alloca(memor_loc,trim(my_lelbo_name),'meshes_list_boundary_elements',lelbo,nboun)
       if( nboun > 0 ) then          
          iboun = 0
          do ielem = 1,nelem
             ielty = abs(ltype(ielem))                    
             do iface = 1,element_type(ielty) % number_faces
                if( facel(mnodf+1,iface,ielem) > 0 .and. facel(1,iface,ielem) > 0 ) then
                   iboun = iboun + 1
                   lelbo(iboun) = ielem
                end if
             end do
          end do
       end if
       !
       ! Boundary connectivity (LNODB)
       !
       if( present(lnodb) ) then
          call memory_alloca(memor_loc,trim(my_lnodb_name),'meshes_list_boundary_elements',lnodb,mnodf,nboun)
       end if
       if( present(ltypb) ) then
          call memory_alloca(memor_loc,trim(my_ltypb_name),'meshes_list_boundary_elements',ltypb,nboun)
       end if
       if( present(lnodb) .or. present(ltypb) ) then
          iboun = 0
          do ielem = 1,nelem
             ielty = abs(ltype(ielem))                    
             do iface = 1,element_type(ielty) % number_faces
                if( facel(mnodf+1,iface,ielem) > 0 .and. facel(1,iface,ielem) > 0 ) then
                   iboun = iboun + 1
                   pnodf = element_type(ielty) % node_faces(iface) 
                   do inodf = 1,pnodf 
                      inode = element_type(ielty) % list_faces(inodf,iface)
                      if( present(lnodb) ) & 
                           lnodb(inodf,iboun) = lnods(inode,ielem)
                   end do
                   if( present(ltypb) ) then
                      ltypb(iboun) = element_type(ielty) % type_faces(iface) 
                   end if
                end if
             end do
          end do
       end if
       !
       ! Deallocate memory
       !
       call memory_deallo(memor_loc,'FACEL'            ,'meshes_list_boundary_elements',facel)
       call memory_deallo(memor_loc,'LELFA'            ,'meshes_list_boundary_elements',lelfa)
       call graphs_elepoi_deallocate(pelpo_loc,lelpo_loc,&
            PELPO_NAME='PELPO_LOC',LELPO_NAME='LELPO_LOC',memor=memor_dom)

    end if

    if( present(MEMORY_COUNTER) ) then
       MEMORY_COUNTER = memor_loc
    else
       memor_dom = memor_loc 
    end if

  end subroutine meshes_list_boundary_elements

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-10-02
  !> @brief   Check element connectivity
  !> @details Check element connectvity and detect floating nodes
  !>          (nodes not belonging to any element)
  !> 
  !-----------------------------------------------------------------------
  
  subroutine meshes_check_mesh(mesh,message)

    type(mesh_type),               intent(in)    :: mesh
    character(len=:), allocatable, intent(inout) :: message
    integer(ip)                                  :: ierro_element_wrong_node
    integer(ip)                                  :: ierro_floating_node
    integer(ip)                                  :: ierro_boundary_wrong_node
    integer(ip)                                  :: ierro_local_node
    integer(ip)                                  :: ierro_lelbo
    integer(ip)                                  :: npoin
    integer(ip)                                  :: nelem
    integer(ip)                                  :: nboun
    integer(ip), pointer                         :: lnods(:,:)
    integer(ip), pointer                         :: lnodb(:,:)
    integer(ip), pointer                         :: ltype(:)
    integer(ip), pointer                         :: ltypb(:)
    integer(ip), pointer                         :: lelbo(:)
    integer(ip), pointer                         :: lbinv(:)
    integer(ip), pointer                         :: lninv(:)
    integer(ip), pointer                         :: leinv(:)
    integer(ip)                                  :: ielem,ipoin,inode,pnode
    integer(ip)                                  :: iboun,inodb,isize
    integer(ip), pointer                         :: touch_nodes(:)
    integer(ip)                                  :: ierro_all(5)
    logical(lg)                                  :: ifoun

    npoin =  mesh % npoin
    nelem =  mesh % nelem
    nboun =  mesh % nboun
    lnods => mesh % lnods
    lnodb => mesh % lnodb
    ltype => mesh % ltype
    ltypb => mesh % ltypb
    lelbo => mesh % lelbo
    lbinv => mesh % lbinv_loc
    lninv => mesh % lninv_loc
    leinv => mesh % lninv_loc

    nullify(touch_nodes)

    call memory_alloca(memor_dom,'TOUCH_NODES','meshes_check_element_connectivity',touch_nodes,npoin)

    ierro_element_wrong_node  = 0
    ierro_floating_node       = 0
    ierro_boundary_wrong_node = 0
    ierro_local_node          = 0
    ierro_lelbo               = 0
    !
    ! Types
    !
    !
    ! Element connectivity 
    !
    do ielem = 1,nelem
       do inode = 1,element_type(abs(ltype(ielem))) % number_nodes 
          ipoin = lnods(inode,ielem)
          if( ipoin < 1 .or. ipoin > npoin ) then
             ierro_element_wrong_node = leinv(ielem)
          else
             touch_nodes(ipoin) = 1 
          end if
       end do
    end do
    !
    ! Boundary connectivity
    !
    loop_iboun: do iboun = 1,nboun
       ielem = lelbo(iboun)
       if( ielem < 1 .or. ielem > nelem ) then
          ierro_lelbo = lbinv(iboun)
          exit loop_iboun
       end if
       pnode = element_type(abs(ltype(ielem))) % number_nodes
       do inodb = 1,element_type(abs(ltypb(iboun))) % number_nodes
          ipoin = lnodb(inodb,iboun)
          if( ipoin == 0 ) then
             ierro_boundary_wrong_node = 1
          else
             ifoun = .false.
             loop_inode: do inode = 1,pnode
                if( lnods(inode,ielem) == ipoin ) then
                   ifoun = .true.
                   exit loop_inode
                end if
             end do loop_inode
             if( .not. ifoun ) then
                !write(*,*) 'LOCAL NUMBERING NOT FOUND FOR BOUNDARY '//trim(intost(iboun))//' IN ELEMENT '//trim(intost(ielem))
                ierro_local_node = lbinv(iboun)
             end if
          end if
       end do
    end do loop_iboun

    call PAR_INTERFACE_NODE_EXCHANGE(touch_nodes,'SUM')

    do ipoin = 1,npoin
       if( touch_nodes(ipoin) == 0 ) then
          ierro_floating_node = lninv(ipoin)
       end if
    end do

    call memory_deallo(memor_dom,'TOUCH_NODES','meshes_check_element_connectivity',touch_nodes)

    ierro_all(1) = ierro_element_wrong_node
    ierro_all(2) = ierro_floating_node
    ierro_all(3) = ierro_boundary_wrong_node 
    ierro_all(4) = ierro_local_node   
    ierro_all(5) = ierro_lelbo
    isize        = size(ierro_all,KIND=ip)

    call PAR_MAX(isize,ierro_all)
    
    ierro_element_wrong_node  = ierro_all(1)
    ierro_floating_node       = ierro_all(2)
    ierro_boundary_wrong_node = ierro_all(3) 
    ierro_local_node          = ierro_all(4) 
    ierro_lelbo               = ierro_all(5) 

    if( ierro_floating_node       /= 0 ) call meshes_add_message(message,'THERE ARE NODES WITHOUT ELEMENTS, CHECK NODE '//integer_to_string(ierro_floating_node))
    if( ierro_local_node          /= 0 ) call meshes_add_message(message,'WRONG BOUNDARY CONNECTIVITY')
    if( ierro_element_wrong_node  /= 0 ) call meshes_add_message(message,'SOME ELEMENTS HAVE NODES OUT OF RANGE, CHECK ELEMENT '//integer_to_string(ierro_element_wrong_node))
    if( ierro_boundary_wrong_node /= 0 ) call meshes_add_message(message,'SOME BOUNDARY ELEMENTS HAVE NULL NODES')
    if( ierro_lelbo               /= 0 ) call meshes_add_message(message,'WRONG BOUNDARY/ELEMENT CONNECTIVITY FOR BOUNDARY '//integer_to_string(ierro_lelbo))

  end subroutine meshes_check_mesh

  subroutine meshes_add_message(message,message_in)
    
    character(len=:), allocatable, intent(inout) :: message
    character(len=*),              intent(in)    :: message_in

    if( allocated(message) ) then      
       message = message // '; '// message_in
    else
       message = message_in
    end if
    
  end subroutine meshes_add_message
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-07
  !> @brief   Gather a submesh
  !> @details Construcut a mesh in parallel arounbd a given point
  !> 
  !-----------------------------------------------------------------------

  subroutine meshes_gather_submesh(xx_in,diameter,NODE_NUMBER)

    use def_domain
    use def_master
    use mod_parall
    use mod_communications
    use mod_memory
    use mod_output
    use mod_iofile
    use mod_mesh_type
    use mod_messages, only : messages_live
    
    real(rp),       intent(in)           :: xx_in(3)
    real(rp),       intent(in)           :: diameter
    integer(ip),    intent(in), optional :: NODE_NUMBER
    integer(ip)                          :: ipoin,idime,knode,ielem,kpoin,kboun
    integer(ip)                          :: inode,ineig,iboun,kelem,nneig_gat
    integer(ip)                          :: jpoin,ndim1,ndim2,ifiel,nfiel_sav
    integer(ip)                          :: nelem_sub,npoin_sub,nboun_sub
    integer(ip)                          :: kneig
    real(rp)                             :: dummr,xx(3)
    type(r3p)                            :: xfiel_sub(mfiel)
    
    integer(ip),    pointer              :: mark_node(:)
    logical(lg),    pointer              :: mark_elem(:)
    logical(lg),    pointer              :: mark_boun(:)
    integer(ip),    pointer              :: permu_node(:)
    integer(ip),    pointer              :: permu_elem(:)
    integer(ip),    pointer              :: permu_boun(:)
    integer(ip),    pointer              :: invpe_node(:)
    integer(ip),    pointer              :: invpe_elem(:)

    integer(ip),    pointer              :: nelem_gat(:)
    integer(ip),    pointer              :: npoin_gat(:)
    integer(ip),    pointer              :: nboun_gat(:)
    integer(4),     pointer              :: lsiz4_gat(:)

    integer(ip),    pointer              :: lnods_sub(:,:)
    integer(ip),    pointer              :: ltype_sub(:)
    integer(ip),    pointer              :: lesub_sub(:)
    integer(ip),    pointer              :: lmate_sub(:)
    integer(ip),    pointer              :: lnnod_sub(:)
    integer(ip),    pointer              :: leset_sub(:)

    integer(ip),    pointer              :: lnodb_sub(:,:)
    integer(ip),    pointer              :: ltypb_sub(:)
    integer(ip),    pointer              :: lnnob_sub(:)
    integer(ip),    pointer              :: lelbo_sub(:)
    integer(ip),    pointer              :: codbo_sub(:)
    integer(ip),    pointer              :: lbset_sub(:)

    real(rp),       pointer              :: coord_sub(:,:)
    integer(ip),    pointer              :: lninv_sub(:)

    real(rp),       pointer              :: xarray(:)
    
    integer(ip)                          :: nunit_msh
    integer(ip)                          :: nunit_res
    integer(ip)                          :: nunit_dat
    type(mesh_type)                      :: meshe_sub

    call mesh_type_initialize(meshe_sub)
    nfiel_sav = nfiel
    
    nullify(mark_node)
    nullify(mark_elem)
    nullify(mark_boun)
    nullify(permu_node)
    nullify(permu_elem)
    nullify(permu_boun)
    nullify(invpe_node)
    nullify(invpe_elem)

    nullify(nelem_gat)
    nullify(npoin_gat)
    nullify(nboun_gat)
    nullify(lsiz4_gat)

    nullify(lnods_sub)
    nullify(ltype_sub)
    nullify(lesub_sub)
    nullify(lmate_sub)
    nullify(lnnod_sub)
    nullify(leset_sub)

    nullify(lnodb_sub)
    nullify(ltypb_sub)
    nullify(lnnob_sub)
    nullify(lelbo_sub)
    nullify(codbo_sub)
    nullify(lbset_sub)

    nullify(coord_sub)
    nullify(lninv_sub)

    nullify(xarray)
    
    npoin_sub = 0
    nelem_sub = 0
    nboun_sub = 0
    allocate(mark_node(npoin))
    allocate(mark_elem(nelem))
    allocate(mark_boun(nboun))

    if( present(NODE_NUMBER) ) then
       ipoin_loop: do ipoin = 1,npoin_own
          if( lninv_loc(ipoin) == NODE_NUMBER ) then
             xx(1:ndime) = coord(1:ndime,ipoin)
             exit ipoin_loop
          end if
       end do ipoin_loop
       call PAR_MAX(ndime,xx)
    else
       xx(1:ndime) = xx_in(1:ndime)
    end if
    
    !--------------------------------------------------------------------
    !
    ! Create submeshes from marked nodes
    !
    !--------------------------------------------------------------------
    
    if( INOTMASTER ) then
       !
       ! Mark nodes
       !
       do ipoin = 1,npoin
          dummr = 0.0_rp
          do idime = 1,ndime
             dummr = dummr + (coord(idime,ipoin)-xx(idime))*(coord(idime,ipoin)-xx(idime))
          end do
          if( sqrt(dummr) <= diameter ) then
             mark_node(ipoin) = 1
             npoin_sub     = 1
          else
             mark_node(ipoin) = 0
          end if
       end do
       if( npoin_sub == 1 ) then
          !
          ! Mark elements
          !
          nelem_sub = 0
          do ielem = 1,nelem
             mark_elem(ielem) = .false.
             loop_inode: do inode = 1,lnnod(ielem)
                ipoin = lnods(inode,ielem)
                if( mark_node(ipoin) == 1 ) then
                   do knode = 1,lnnod(ielem)
                      kpoin = lnods(knode,ielem)
                      if( mark_node(kpoin) == 0 ) mark_node(kpoin) = 2
                   end do
                   nelem_sub = nelem_sub + 1
                   mark_elem(ielem) = .true.
                   exit loop_inode
                end if
             end do loop_inode
          end do
          !
          ! Mark boundaries
          !
          nboun_sub = 0
          do iboun = 1,nboun
             ielem = lelbo(iboun)
             mark_boun(iboun) = mark_elem(ielem)
             if( mark_boun(iboun) ) nboun_sub = nboun_sub + 1
          end do
          !
          ! Renumber nodes
          !
          npoin_sub = 0
          do ipoin = 1,npoin
             if( mark_node(ipoin) /= 0 ) then
                npoin_sub = npoin_sub + 1
                mark_node(ipoin) = npoin_sub
             end if
          end do
          !
          ! Permutation
          !
          allocate(permu_node(npoin_sub))
          allocate(permu_elem(nelem_sub))
          allocate(permu_boun(nelem_sub))
          allocate(invpe_node(npoin))
          allocate(invpe_elem(nelem))
          npoin_sub = 0
          do ipoin = 1,npoin
             if( mark_node(ipoin) /= 0 ) then
                npoin_sub             = npoin_sub + 1
                permu_node(npoin_sub) = ipoin
                invpe_node(ipoin)     = npoin_sub
             end if
          end do
          nelem_sub = 0
          do ielem = 1,nelem
             if( mark_elem(ielem) ) then
                nelem_sub             = nelem_sub + 1
                permu_elem(nelem_sub) = ielem
                invpe_elem(ielem)     = nelem_sub
             end if
          end do
          nboun_sub = 0
          do iboun = 1,nboun
             if( mark_boun(iboun) ) then
                nboun_sub             = nboun_sub + 1
                permu_boun(nboun_sub) = iboun
             end if
          end do
          !
          ! Create submesh
          !
          call memory_alloca(memor_dom,'LNODS','mod_meshes',lnods_sub,mnode,nelem_sub)
          call memory_alloca(memor_dom,'LTYPE','mod_meshes',ltype_sub,nelem_sub)
          call memory_alloca(memor_dom,'LESUB','mod_meshes',lesub_sub,nelem_sub)
          call memory_alloca(memor_dom,'LMATE','mod_meshes',lmate_sub,nelem_sub)
          call memory_alloca(memor_dom,'LNNOD','mod_meshes',lnnod_sub,nelem_sub)
          if( neset > 0 ) call memory_alloca(memor_dom,'LESET','mod_meshes',leset_sub,nelem_sub)
         
          call memory_alloca(memor_dom,'LNODS','mod_meshes',lnodb_sub,mnodb,nboun_sub)
          call memory_alloca(memor_dom,'LTYPE','mod_meshes',ltypb_sub,nboun_sub)
          call memory_alloca(memor_dom,'LNNOB','mod_meshes',lnnob_sub,nboun_sub)
          call memory_alloca(memor_dom,'LELBO','mod_meshes',lelbo_sub,nboun_sub)
          call memory_alloca(memor_dom,'CODBO','mod_meshes',codbo_sub,nboun_sub)
          if( nbset > 0 ) call memory_alloca(memor_dom,'LBSET','mod_meshes',lbset_sub,nboun_sub)

          call memory_alloca(memor_dom,'COORD','mod_meshes',coord_sub,ndime,npoin_sub)
          call memory_alloca(memor_dom,'LNINV','mod_meshes',lninv_sub,npoin_sub)

          do ifiel = 1,nfiel 
             ndim1 = kfl_field(1,ifiel)
             ndim2 = kfl_field(4,ifiel)
             if(      kfl_field(2,ifiel) == NELEM_TYPE ) then
                call memory_alloca(memor_dom,'XFIEL','mod_meshes',xfiel_sub(ifiel) % a,ndim1,nelem_sub,ndim2)
             else if( kfl_field(2,ifiel) == NPOIN_TYPE ) then                                                                        
                call memory_alloca(memor_dom,'XFIEL','mod_meshes',xfiel_sub(ifiel) % a,ndim1,npoin_sub,ndim2)
             else if( kfl_field(2,ifiel) == NBOUN_TYPE ) then                                                                        
                call memory_alloca(memor_dom,'XFIEL','mod_meshes',xfiel_sub(ifiel) % a,ndim1,nboun_sub,ndim2)
             end if                       
          end do
          
          do ielem = 1,nelem_sub
             kelem              = permu_elem(ielem)
             lnods_sub(:,ielem) = invpe_node(lnods(:,kelem))
             ltype_sub(ielem)   = ltype(kelem)
             lesub_sub(ielem)   = lesub(kelem)
             lmate_sub(ielem)   = lmate(kelem)
             lnnod_sub(ielem)   = lnnod(kelem)
             if( neset > 0 ) leset_sub(ielem)   = leset(kelem)
          end do
          do iboun = 1,nboun_sub
             kboun              = permu_boun(iboun)
             lnodb_sub(:,iboun) = invpe_node(lnodb(:,kboun))
             ltypb_sub(iboun)   = ltypb(kboun)
             lnnob_sub(iboun)   = lnnob(kboun)
             lelbo_sub(iboun)   = invpe_elem(lelbo(kboun))
             codbo_sub(iboun)   = kfl_codbo(kboun)
             if( nbset > 0 ) lbset_sub(iboun)   = lbset(kboun)
          end do
          do ipoin = 1,npoin_sub
             kpoin              = permu_node(ipoin)
             coord_sub(:,ipoin) = coord(:,kpoin)
             lninv_sub(ipoin)   = lninv_loc(kpoin)
          end do
          do ifiel = 1,nfiel
             if(      kfl_field(2,ifiel) == NELEM_TYPE ) then
                do ielem = 1,nelem_sub
                   kelem = permu_elem(ielem)
                   xfiel_sub(ifiel) % a(:,ielem,:) = xfiel(ifiel) % a(:,kelem,:) 
                end do
             else if( kfl_field(2,ifiel) == NPOIN_TYPE ) then                                                                        
                do ipoin = 1,npoin_sub
                   kpoin = permu_node(ipoin)
                   xfiel_sub(ifiel) % a(:,ipoin,:) = xfiel(ifiel) % a(:,kpoin,:) 
                end do
             else if( kfl_field(2,ifiel) == NBOUN_TYPE ) then                                                                        
                do iboun = 1,nboun_sub
                   kboun = permu_boun(iboun)
                   xfiel_sub(ifiel) % a(:,iboun,:) = xfiel(ifiel) % a(:,kboun,:) 
                end do
             end if                       
          end do
          if( kfl_ngrou /= 0 ) then
             nfiel              = nfiel + 1
             ifiel              = nfiel
             kfl_field(1,ifiel) = 1
             kfl_field(2,ifiel) = NPOIN_TYPE
             kfl_field(4,ifiel) = 1
             ndim1              = kfl_field(1,ifiel)
             ndim2              = kfl_field(4,ifiel)
             call memory_alloca(memor_dom,'XFIEL','mod_meshes',xfiel_sub(ifiel) % a,ndim1,npoin_sub,ndim2)
             do ipoin = 1,npoin_sub
                kpoin = permu_node(ipoin)
                xfiel_sub(ifiel) % a(1,ipoin,1) = real(lgrou_dom(kpoin),rp) 
             end do
          end if
       end if

    end if

    !--------------------------------------------------------------------
    !
    ! Gather geometry
    !
    !--------------------------------------------------------------------
    
    if( INOTSLAVE ) then
       nneig_gat = PAR_WORLD_SIZE-1
       allocate( nelem_gat(0:PAR_WORLD_SIZE-1) )
       allocate( nboun_gat(0:PAR_WORLD_SIZE-1) )
       allocate( npoin_gat(0:PAR_WORLD_SIZE-1) )
       allocate( lsiz4_gat(0:PAR_WORLD_SIZE-1) )
    else
       nneig_gat = -1
    end if
    
    call PAR_GATHER(nelem_sub,nelem_gat,'IN THE WORLD') 
    call PAR_GATHER(nboun_sub,nboun_gat,'IN THE WORLD') 
    call PAR_GATHER(npoin_sub,npoin_gat,'IN THE WORLD')

    if( INOTSLAVE ) then
       nelem_sub = sum(nelem_gat)
       nboun_sub = sum(nboun_gat)
       npoin_sub = sum(npoin_gat)
       call memory_alloca(memor_dom,'LNODS','mod_meshes',lnods_sub,mnode,nelem_sub)
       call memory_alloca(memor_dom,'LTYPE','mod_meshes',ltype_sub,nelem_sub)
       call memory_alloca(memor_dom,'LESUB','mod_meshes',lesub_sub,nelem_sub)
       call memory_alloca(memor_dom,'LMATE','mod_meshes',lmate_sub,nelem_sub)
       call memory_alloca(memor_dom,'LNNOD','mod_meshes',lnnod_sub,nelem_sub)
       if( neset > 0 ) call memory_alloca(memor_dom,'LESET','mod_meshes',leset_sub,nelem_sub)

       call memory_alloca(memor_dom,'LNODB','mod_meshes',lnodb_sub,mnodb,nboun_sub)
       call memory_alloca(memor_dom,'LTYPB','mod_meshes',ltypb_sub,nboun_sub)
       call memory_alloca(memor_dom,'LNNOB','mod_meshes',lnnob_sub,nboun_sub)
       call memory_alloca(memor_dom,'LELBO','mod_meshes',lelbo_sub,nboun_sub)
       call memory_alloca(memor_dom,'CODBO','mod_meshes',codbo_sub,nboun_sub)
       if( nbset > 0 ) call memory_alloca(memor_dom,'LBSET','mod_meshes',lbset_sub,nboun_sub)

       call memory_alloca(memor_dom,'COORD','mod_meshes',coord_sub,ndime,npoin_sub)
       call memory_alloca(memor_dom,'LNINV','mod_meshes',lninv_sub,npoin_sub)
       
       do ifiel = 1,nfiel
          ndim1 = kfl_field(1,ifiel)
          ndim2 = kfl_field(4,ifiel)
          if(      kfl_field(2,ifiel) == NELEM_TYPE ) then
             call memory_alloca(memor_dom,'XFIEL','mod_meshes',xfiel_sub(ifiel) % a,ndim1,nelem_sub,ndim2)
          else if( kfl_field(2,ifiel) == NPOIN_TYPE ) then                                                                        
             call memory_alloca(memor_dom,'XFIEL','mod_meshes',xfiel_sub(ifiel) % a,ndim1,npoin_sub,ndim2)
          else if( kfl_field(2,ifiel) == NBOUN_TYPE ) then                                                                        
             call memory_alloca(memor_dom,'XFIEL','mod_meshes',xfiel_sub(ifiel) % a,ndim1,nboun_sub,ndim2)
          end if
       end do
       !
       ! Add field for groups
       !
       if( kfl_ngrou /= 0 ) then 
          nfiel = nfiel + 1
          ifiel = nfiel
          if( ifiel > mfiel ) call runend('MOD_MESHES: TOO MANY FIELDS!')
          kfl_field(1,ifiel) = 1
          kfl_field(2,ifiel) = NPOIN_TYPE
          kfl_field(4,ifiel) = 1
          ndim1              = kfl_field(1,ifiel)
          ndim2              = kfl_field(4,ifiel)
          call memory_alloca(memor_dom,'XFIEL','mod_meshes',xfiel_sub(ifiel) % a,ndim1,npoin_sub,ndim2)
       end if
    end if
    !
    ! ELement arrays
    !
    do ineig = 0,nneig_gat 
       lsiz4_gat(ineig) = int(mnode * nelem_gat(ineig),4)
    end do
    call PAR_GATHERV(lnods_sub,lnods_sub,lsiz4_gat,'IN THE WORLD')

    do ineig = 0,nneig_gat 
       lsiz4_gat(ineig) = int(nelem_gat(ineig),4)
    end do
    call PAR_GATHERV(ltype_sub,ltype_sub,lsiz4_gat,'IN THE WORLD')
    call PAR_GATHERV(lesub_sub,lesub_sub,lsiz4_gat,'IN THE WORLD')
    call PAR_GATHERV(lmate_sub,lmate_sub,lsiz4_gat,'IN THE WORLD')
    call PAR_GATHERV(lnnod_sub,lnnod_sub,lsiz4_gat,'IN THE WORLD')
    if( neset > 0 ) call PAR_GATHERV(leset_sub,leset_sub,lsiz4_gat,'IN THE WORLD')
    !
    ! Boundary arrays
    !
    do ineig = 0,nneig_gat 
       lsiz4_gat(ineig) = int(mnodb * nboun_gat(ineig),4)
    end do
    call PAR_GATHERV(lnodb_sub,lnodb_sub,lsiz4_gat,'IN THE WORLD')
    do ineig = 0,nneig_gat 
       lsiz4_gat(ineig) = int(nboun_gat(ineig),4)
    end do
    call PAR_GATHERV(ltypb_sub,ltypb_sub,lsiz4_gat,'IN THE WORLD')
    call PAR_GATHERV(lnnob_sub,lnnob_sub,lsiz4_gat,'IN THE WORLD')
    call PAR_GATHERV(lelbo_sub,lelbo_sub,lsiz4_gat,'IN THE WORLD')
    call PAR_GATHERV(codbo_sub,codbo_sub,lsiz4_gat,'IN THE WORLD')
    if( nbset > 0 ) call PAR_GATHERV(lbset_sub,lbset_sub,lsiz4_gat,'IN THE WORLD')
    !
    ! Node arrays
    !
    do ineig = 0,nneig_gat 
       lsiz4_gat(ineig) = int(ndime * npoin_gat(ineig),4)
    end do
    call PAR_GATHERV(coord_sub,coord_sub,lsiz4_gat,'IN THE WORLD')
    do ineig = 0,nneig_gat 
       lsiz4_gat(ineig) = int(npoin_gat(ineig),4)
    end do
    call PAR_GATHERV(lninv_sub,lninv_sub,lsiz4_gat,'IN THE WORLD')
    !
    ! Fields
    !
    do ifiel = 1,nfiel
       ndim1 = kfl_field(1,ifiel)
       ndim2 = kfl_field(4,ifiel)
       if(      kfl_field(2,ifiel) == NELEM_TYPE ) then
          do ineig = 0,nneig_gat 
             lsiz4_gat(ineig) = int(ndim1*ndim2*nelem_gat(ineig),4)
          end do
          call PAR_GATHERV(xfiel_sub(ifiel) % a,xfiel_sub(ifiel) % a,lsiz4_gat,'IN THE WORLD')
       else if( kfl_field(2,ifiel) == NPOIN_TYPE ) then                                                                        
          do ineig = 0,nneig_gat 
             lsiz4_gat(ineig) = int(ndim1*ndim2*npoin_gat(ineig),4)
          end do
          call PAR_GATHERV(xfiel_sub(ifiel) % a,xfiel_sub(ifiel) % a,lsiz4_gat,'IN THE WORLD')
       else if( kfl_field(2,ifiel) == NBOUN_TYPE ) then                                                                        
          do ineig = 0,nneig_gat 
             lsiz4_gat(ineig) = int(ndim1*ndim2*nboun_gat(ineig),4)
          end do
          call PAR_GATHERV(xfiel_sub(ifiel) % a,xfiel_sub(ifiel) % a,lsiz4_gat,'IN THE WORLD')
       end if
    end do
    !
    ! Renumber nodes
    !
    if( INOTSLAVE ) then
       ielem = 0
       iboun = 0
       kpoin = 0
       do ineig = 0,nneig_gat
          do kelem = 1,nelem_gat(ineig)
             ielem = ielem + 1
             lnods_sub(:,ielem) = lnods_sub(:,ielem) + kpoin
          end do
          do kboun = 1,nboun_gat(ineig)
             iboun = iboun + 1
             lnodb_sub(:,iboun) = lnodb_sub(:,iboun) + kpoin
          end do
          kpoin = kpoin + npoin_gat(ineig)
       end do
    end if
    !
    ! Renumber boundaries
    !
    if( INOTSLAVE ) then
       iboun = 0
       kelem = 0
       do ineig = 0,nneig_gat
          do kboun = 1,nboun_gat(ineig)
             iboun = iboun + 1
             lelbo_sub(iboun) = lelbo_sub(iboun) + kelem
          end do
          kelem = kelem + nelem_gat(ineig)
       end do
    end if
    !
    ! Eliminate duplicated nodes
    !
    if( INOTSLAVE ) then
       
       allocate(permu_node(npoin_sub))
       do ipoin = 1,npoin_sub
          permu_node(ipoin) = 0
       end do
       kpoin = 0
       do ipoin = 1,npoin_sub
          if( permu_node(ipoin) == 0 ) then
             kpoin = kpoin + 1
             permu_node(ipoin) = kpoin
             do jpoin = ipoin+1,npoin_sub
                if( lninv_sub(jpoin) == lninv_sub(ipoin) ) then
                   permu_node(jpoin) = -permu_node(ipoin)
                end if
             end do
          end if
       end do
       
       jpoin = 0
       do ipoin = 1,npoin_sub
          if( permu_node(ipoin) > 0 ) then
             jpoin = jpoin + 1
             coord_sub(:,jpoin) = coord_sub(:,ipoin)
          end if
       end do
       if( kfl_ngrou /= 0 ) then
          jpoin = 0
          do ipoin = 1,npoin_sub
             if( permu_node(ipoin) > 0 ) then
                jpoin = jpoin + 1
                xfiel_sub(nfiel) % a(1,jpoin,1) = xfiel_sub(nfiel) % a(1,ipoin,1)
             end if
          end do
       end if
       
       npoin_sub = kpoin
       do ielem = 1,nelem_sub
          lnods_sub(:,ielem) = abs(permu_node(lnods_sub(:,ielem)))
       end do
       do iboun = 1,nboun_sub
          lnodb_sub(:,iboun) = abs(permu_node(lnodb_sub(:,iboun)))
       end do
       
    end if
    
    !--------------------------------------------------------------------
    !
    ! Add a field for parallelization
    !
    !--------------------------------------------------------------------
    
     if( INOTSLAVE ) then
       nfiel = nfiel + 1
       ifiel = nfiel
       if( ifiel > mfiel ) call runend('MOD_MESHES: TOO MANY FIELDS!')
       kfl_field(1,ifiel) = 1
       kfl_field(2,ifiel) = NELEM_TYPE
       kfl_field(4,ifiel) = 1
       ndim1              = kfl_field(1,ifiel)
       ndim2              = kfl_field(4,ifiel)
       call memory_alloca(memor_dom,'XFIEL','mod_meshes',xfiel_sub(ifiel) % a,ndim1,nelem_sub,ndim2)
       ielem = 0
       kneig = 0
       do ineig = 0,nneig_gat
          if( nelem_gat(ineig) > 0 ) then
             kneig = kneig + 1
             do kelem = 1,nelem_gat(ineig)
                ielem = ielem + 1
                xfiel_sub(ifiel) % a(1,ielem,1) = real(kneig,rp)
             end do
          end if
       end do       
    end if
    if( IMASTER ) then
       call messages_live('NUMBER OF PARTITIONS INVOLVED IN THIS SUBMESH= ',INT_NUMBER=kneig)       
    end if
    
    !--------------------------------------------------------------------
    !
    ! Output mesh
    !
    !--------------------------------------------------------------------
    
    nunit_msh = 90
    nunit_res = 91
    nunit_dat = 92
    call iofile_open_unit(nunit_msh,trim(title)//'-submesh.post.msh','POSTPROCESS SUBMESH')
    call iofile_open_unit(nunit_res,trim(title)//'-submesh.post.res','POSTPROCESS SUBMESH')
    call iofile_open_unit(nunit_dat,trim(title)//'-submesh.dom.dat', 'POSTPROCESS SUBMESH')

    if( INOTSLAVE ) then
       meshe_sub % ndime     =  ndime
       meshe_sub % mnode     =  mnode
       meshe_sub % mnodb     =  mnodb
       meshe_sub % nfiel     =  nfiel
       meshe_sub % kfl_field =  kfl_field
       meshe_sub % nelem     =  nelem_sub
       meshe_sub % nboun     =  nboun_sub
       meshe_sub % npoin     =  npoin_sub
       if( kfl_ngrou /= 0 ) meshe_sub % kfl_ngrou =  nfiel-1
       meshe_sub % lnods     => lnods_sub
       meshe_sub % ltype     => ltype_sub
       meshe_sub % lesub     => lesub_sub
       meshe_sub % lmate     => lmate_sub
       meshe_sub % lnnod     => lnnod_sub
       meshe_sub % leset     => leset_sub
       meshe_sub % lnodb     => lnodb_sub
       meshe_sub % ltypb     => ltypb_sub
       meshe_sub % lnnob     => lnnob_sub
       meshe_sub % lelbo     => lelbo_sub
       meshe_sub % kfl_codbo => codbo_sub
       meshe_sub % lbset     => lbset_sub
       meshe_sub % coord     => coord_sub
       do ifiel = 1,nfiel
          meshe_sub % xfiel(ifiel) % a => xfiel_sub(ifiel) % a
       end do
       call output_mesh_gid_format   (meshe_sub,'SUBMESH',nunit_msh)
       call output_domain_alya_format(meshe_sub,nunit_dat)

       allocate(xarray(nelem_sub))
       do ielem = 1,nelem_sub
          xarray(ielem) = meshe_sub % xfiel(nfiel)%a(1,ielem,1)
       end do
       call output_result_gid_format (nunit_res,nelem_sub,xarray,NAME='PARTITION',wherein='ON ELEMENTS')
       deallocate(xarray)
       if( meshe_sub % kfl_ngrou > 0 ) then
          allocate(xarray(npoin_sub))
          do ipoin = 1,npoin_sub
             xarray(ipoin) = meshe_sub % xfiel(nfiel-1)%a(1,ipoin,1)
          end do
          call output_result_gid_format (nunit_res,npoin_sub,xarray,NAME='GROUPS')
          deallocate(xarray)
       end if
    end if
    
    !--------------------------------------------------------------------
    !
    ! Deallocate 
    !
    !--------------------------------------------------------------------
    
    call memory_deallo(memor_dom,'LNODS','mod_meshes',lnods_sub)
    call memory_deallo(memor_dom,'LTYPE','mod_meshes',ltype_sub)
    call memory_deallo(memor_dom,'LESUB','mod_meshes',lesub_sub)
    call memory_deallo(memor_dom,'LMATE','mod_meshes',lmate_sub)
    call memory_deallo(memor_dom,'LNNOD','mod_meshes',lnnod_sub)
    call memory_deallo(memor_dom,'LESET','mod_meshes',leset_sub)
    
    call memory_deallo(memor_dom,'LNODB','mod_meshes',lnodb_sub)
    call memory_deallo(memor_dom,'LTYPB','mod_meshes',ltypb_sub)
    call memory_deallo(memor_dom,'LNNOB','mod_meshes',lnnob_sub)
    call memory_deallo(memor_dom,'LELBO','mod_meshes',lelbo_sub)
    call memory_deallo(memor_dom,'CODBO','mod_meshes',codbo_sub)
    call memory_deallo(memor_dom,'LBSET','mod_meshes',lbset_sub)
    
    call memory_deallo(memor_dom,'COORD','mod_meshes',coord_sub)
    call memory_deallo(memor_dom,'LNINV','mod_meshes',lninv_sub)
    do ifiel = 1,mfiel
       call memory_deallo(memor_dom,'XFIEL','mod_meshes',xfiel_sub(ifiel) % a)
    end do
    
    if( associated(mark_node ) ) deallocate(mark_node )
    if( associated(mark_elem ) ) deallocate(mark_elem )
    if( associated(mark_boun ) ) deallocate(mark_boun )
    if( associated(permu_node) ) deallocate(permu_node)
    if( associated(permu_elem) ) deallocate(permu_elem)
    if( associated(permu_boun) ) deallocate(permu_boun)
    if( associated(invpe_node) ) deallocate(invpe_node)
    if( associated(invpe_elem) ) deallocate(invpe_elem)

    call iofile_close_unit(nunit_msh)
    call iofile_close_unit(nunit_dat)
    call iofile_close_unit(nunit_res)
    nfiel = nfiel_sav 

    call runend('O.K.!')
    
  end subroutine meshes_gather_submesh

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux and borrell
  !> @date    2020-02-07
  !> @brief   Glue two meshes
  !> @details Glue mesh 2 to mesh 1
  !> 
  !-----------------------------------------------------------------------

  subroutine meshes_glue_two_meshes(meshe_i,meshe_o,COLLAPSE,array_i,array_o)

    type(mesh_type),                    intent(inout) :: meshe_i         !< Input mesh (eventually deallocated)
    type(mesh_type),                    intent(inout) :: meshe_o         !< Output mesh
    logical(lg),     optional,          intent(in)    :: COLLAPSE        !< Nodes should be collapsed
    integer(ip),     optional, pointer, intent(inout) :: array_i(:,:)    !< Extra input array
    integer(ip),     optional, pointer, intent(inout) :: array_o(:,:)    !< Extra output array
    integer(ip)                                       :: nelem_3,npoin_3
    integer(ip),               pointer                :: lnods_3(:,:)
    integer(ip),               pointer                :: ltype_3(:)
    integer(ip),               pointer                :: leinv_3(:)
    integer(ip),               pointer                :: lninv_3(:)
    real(rp),                  pointer                :: coord_3(:,:)
    integer(ip),               pointer                :: lnnod_3(:)
    integer(ip),               pointer                :: array_3(:,:)
    type(r3p)                                         :: xfiel_3(max(1_ip,meshe_i % nfiel))
    type(hash_t)                                      :: htable_3 
    integer(ip)                                       :: lid,ipoin,ielem,kelem
    integer(ip)                                       :: inode
    integer(ip)                                       :: ndima,ndim1,ndim3,ifiel
    logical(lg)                                       :: isin
    
    nullify(lnods_3,ltype_3,leinv_3,lnnod_3,lninv_3,coord_3,array_3)

    if( meshe_i % nelem == 0 ) then

       return
       
    else if( meshe_o % nelem == 0 ) then
       !
       ! Copy first mesh
       !
       meshe_o % nelem = meshe_i % nelem
       meshe_o % npoin = meshe_i % npoin
       meshe_o % mnode = meshe_i % mnode
       meshe_o % ndime = meshe_i % ndime
       meshe_o % mnodb = meshe_i % mnodb
       
       call memory_alloca(memor_dom,'LNODS_O',trim(vacal),meshe_o % lnods,    meshe_o % mnode,meshe_o % nelem)
       call memory_alloca(memor_dom,'LTYPE_O',trim(vacal),meshe_o % ltype,    meshe_o % nelem)
       call memory_alloca(memor_dom,'LNNOD_O',trim(vacal),meshe_o % lnnod,    meshe_o % nelem)
       call memory_alloca(memor_dom,'LEINV_O',trim(vacal),meshe_o % leinv_loc,meshe_o % nelem)

       call memory_alloca(memor_dom,'COORD_O',trim(vacal),meshe_o % coord,    meshe_o % ndime,meshe_o % npoin)
       call memory_alloca(memor_dom,'LNINV_O',trim(vacal),meshe_o % lninv_loc,meshe_o % npoin)

       do ifiel = 1,meshe_i % nfiel
          if( associated(meshe_i % xfiel(ifiel) % a) ) then
             ndim1 = size(meshe_i % xfiel(ifiel) % a,1)
             ndim3 = size(meshe_i % xfiel(ifiel) % a,2)
             if(      meshe_i % kfl_field(2,ifiel) == NELEM_TYPE ) then
                call memory_alloca(memor_dom,'XFIEL_o % A',trim(vacal),meshe_o % xfiel(ifiel) % a,ndim1,meshe_o % nelem,ndim3)
             else if( meshe_i % kfl_field(2,ifiel) == NPOIN_TYPE ) then
                call memory_alloca(memor_dom,'XFIEL_O % A',trim(vacal),meshe_o % xfiel(ifiel) % a,ndim1,meshe_o % npoin,ndim3)                
             end if
          end if
       end do
    
       if( present(array_i) .and. present(array_o) ) then
          ndima = memory_size(array_i,1_ip)
          call memory_alloca(memor_dom,'ARRAY_O',trim(vacal),array_o,ndima,meshe_o % npoin)
       end if

       do ielem = 1,meshe_o % nelem
          meshe_o % lnods(:,ielem)   = meshe_i % lnods(:,ielem) 
          meshe_o % ltype(ielem)     = meshe_i % ltype(ielem) 
          meshe_o % lnnod(ielem)     = meshe_i % lnnod(ielem) 
          meshe_o % leinv_loc(ielem) = meshe_i % leinv_loc(ielem) 
       end do
       do ifiel = 1,meshe_i % nfiel
          if( associated(meshe_i % xfiel(ifiel) % a) .and. meshe_i % kfl_field(2,ifiel) == NELEM_TYPE ) then
             do ielem = 1,meshe_o % nelem
                meshe_o % xfiel(ifiel) % a(:,ielem,:) = meshe_i % xfiel(ifiel) % a(:,ielem,:) 
             end do
          end if
       end do       
       do ipoin = 1,meshe_o % npoin
          meshe_o % coord(:,ipoin)   = meshe_i % coord(:,ipoin) 
          meshe_o % lninv_loc(ipoin) = meshe_i % lninv_loc(ipoin) 
       end do
       do ifiel = 1,meshe_i % nfiel
          if( associated(meshe_i % xfiel(ifiel) % a) .and. meshe_i % kfl_field(2,ifiel) == NPOIN_TYPE ) then
             do ipoin = 1,meshe_o % npoin
                meshe_o % xfiel(ifiel) % a(:,ipoin,:) = meshe_i % xfiel(ifiel) % a(:,ipoin,:) 
             end do
          end if
       end do
       
       if( present(array_i) .and. present(array_o) ) then
          do ipoin = 1,meshe_o % npoin
             array_o(:,ipoin) = array_i(:,ipoin) 
          end do
          call memory_deallo(memor_dom,'ARRAY_I',trim(vacal),array_i)
       end if
       
       call memory_deallo(memor_dom,'LNODS_I',trim(vacal),meshe_i % lnods)
       call memory_deallo(memor_dom,'LTYPE_I',trim(vacal),meshe_i % ltype)
       call memory_deallo(memor_dom,'LNNOD_I',trim(vacal),meshe_i % lnnod)
       call memory_deallo(memor_dom,'LEINV_I',trim(vacal),meshe_i % leinv_loc)

       call memory_deallo(memor_dom,'COORD_I',trim(vacal),meshe_i % coord)
       call memory_deallo(memor_dom,'LNINV_I',trim(vacal),meshe_i % lninv_loc)

       do ifiel = 1,meshe_i % nfiel
          call memory_deallo(memor_dom,'XFIEL_I % A',trim(vacal),meshe_i % xfiel(ifiel) % a)
       end do
       
    else
       !
       ! Append mesh MESH_I to MESH_O
       !
       nelem_3 = meshe_o % nelem + meshe_i % nelem
       npoin_3 = meshe_o % npoin + meshe_i % npoin
       
       call htable_initialization(htable_3)
       call htades( htable_3, memor_opt=memor_dom  )
       call htaini( htable_3, npoin_3, lidson=.true., AUTOMATIC_SIZE=.true.,memor_opt=memor_dom)
       call htaadd( htable_3, meshe_o % lninv_loc, memor_opt=memor_dom)
       do ipoin = 1,meshe_i % npoin
          call htaadd(htable_3,meshe_i % lninv_loc(ipoin),lid,isin)
       end do
       npoin_3 = htable_3 % nelem
       !
       ! Allocate arrays
       !
       call memory_alloca(memor_dom,'LNODS_3',trim(vacal),lnods_3,meshe_i % mnode,nelem_3)
       call memory_alloca(memor_dom,'LTYPE_3',trim(vacal),ltype_3,nelem_3)
       call memory_alloca(memor_dom,'LEINV_3',trim(vacal),leinv_3,nelem_3)
       call memory_alloca(memor_dom,'LNNOD_3',trim(vacal),lnnod_3,nelem_3)
       
       call memory_alloca(memor_dom,'LNINV_3',trim(vacal),lninv_3,npoin_3)
       call memory_alloca(memor_dom,'COORD_3',trim(vacal),coord_3,meshe_i % ndime,npoin_3)
       if( present(array_i) .and. present(array_o) ) then
          ndima = memory_size(array_i,1_ip)
          call memory_alloca(memor_dom,'ARRAY_3',trim(vacal),array_3,ndima,npoin_3)
       end if
       
       do ifiel = 1,meshe_i % nfiel
          if( associated(meshe_i % xfiel(ifiel) % a) ) then
             ndim1 = size(meshe_i % xfiel(ifiel) % a,1)
             ndim3 = size(meshe_i % xfiel(ifiel) % a,2)
             if(      meshe_i % kfl_field(2,ifiel) == NELEM_TYPE ) then
                call memory_alloca(memor_dom,'XFIEL_3',trim(vacal),xfiel_3(ifiel) % a,ndim1,meshe_o % nelem,ndim3)
             else if( meshe_i % kfl_field(2,ifiel) == NPOIN_TYPE ) then
                call memory_alloca(memor_dom,'XFIEL_3',trim(vacal),xfiel_3(ifiel) % a,ndim1,meshe_o % npoin,ndim3)                
             end if
          end if
       end do
       !
       ! Element arrays
       !
       do ielem = 1,meshe_o % nelem
          lnods_3(:,ielem) = meshe_o % lnods(:,ielem)
          ltype_3(ielem)   = meshe_o % ltype(ielem)
          leinv_3(ielem)   = meshe_o % leinv_loc(ielem)
          lnnod_3(ielem)   = meshe_o % lnnod(ielem)
       end do
       kelem =  meshe_o % nelem
       do ielem = 1, meshe_i % nelem
          kelem = kelem + 1
          ltype_3(kelem) = meshe_i % ltype(ielem)
          leinv_3(kelem) = meshe_i % leinv_loc(ielem)
          lnnod_3(kelem) = meshe_i % lnnod(ielem)
          do inode = 1,meshe_i % lnnod(ielem)
             ipoin = meshe_i % lnods(inode,ielem)
             lid   = htalid(htable_3,meshe_i % lninv_loc(ipoin))          
             lnods_3(inode,kelem) = lid
          end do
       end do

       do ifiel = 1,meshe_i % nfiel
          if( associated(meshe_i % xfiel(ifiel) % a) .and. meshe_i % kfl_field(2,ifiel) == NELEM_TYPE ) then
             do ielem = 1,meshe_o % nelem
                xfiel_3(ifiel) % a(:,ielem,:) = meshe_o % xfiel(ifiel) % a(:,ielem,:)
             end do
             kelem =  meshe_o % nelem
             do ielem = 1, meshe_i % nelem
                kelem = kelem + 1
                xfiel_3(ifiel) % a(:,kelem,:) = meshe_o % xfiel(ifiel) % a(:,ielem,:)
             end do
          end if
       end do
       !
       ! Node arrays
       !
       do ipoin = 1,meshe_o % npoin
          lninv_3(ipoin)   = meshe_o % lninv_loc(ipoin)
          coord_3(:,ipoin) = meshe_o % coord(:,ipoin)
       end do
       do ipoin = 1,meshe_i % npoin
          lid = htalid(htable_3,meshe_i % lninv_loc(ipoin))
          coord_3(:,lid) = meshe_i % coord(:,ipoin)
          lninv_3(lid)   = meshe_i % lninv_loc(ipoin)
       end do
       if( present(array_i) .and. present(array_o) ) then
          do ipoin = 1,meshe_o % npoin
             array_3(:,ipoin) = array_o(:,ipoin)
          end do
          do ipoin = 1,meshe_i % npoin
             lid = htalid(htable_3,meshe_i % lninv_loc(ipoin))
             array_3(:,lid) = array_i(:,ipoin)
          end do
       end if
       do ifiel = 1,meshe_i % nfiel
          if( associated(meshe_i % xfiel(ifiel) % a) .and. meshe_i % kfl_field(2,ifiel) == NPOIN_TYPE ) then
             do ipoin = 1,meshe_o % npoin                
                xfiel_3(ifiel) % a(:,ipoin,:) = meshe_o % xfiel(ifiel) % a(:,ipoin,:)
             end do
             do ipoin = 1,meshe_i % npoin
                lid = htalid(htable_3,meshe_i % lninv_loc(ipoin))
                xfiel_3(ifiel) % a(:,lid,:) = meshe_i % xfiel(ifiel) % a(:,ipoin,:)
             end do
          end if
       end do
       
       call memory_deallo(memor_dom,'LNODS_O',trim(vacal),meshe_o % lnods)
       call memory_deallo(memor_dom,'LTYPE_O',trim(vacal),meshe_o % ltype)
       call memory_deallo(memor_dom,'LEINV_O',trim(vacal),meshe_o % leinv_loc)
       call memory_deallo(memor_dom,'LNNOD_O',trim(vacal),meshe_o % lnnod)

       call memory_deallo(memor_dom,'LNINV_O',trim(vacal),meshe_o % lninv_loc)
       call memory_deallo(memor_dom,'COORD_O',trim(vacal),meshe_o % coord)
       if( present(array_i) .and. present(array_o) ) then
          call memory_deallo(memor_dom,'ARRAY_O',trim(vacal),array_o)
       end if
       
       meshe_o % nelem = nelem_3
       meshe_o % npoin = npoin_3
       
       call memory_copy(memor_dom,'LNODS_O',trim(vacal),lnods_3,meshe_o % lnods)
       call memory_copy(memor_dom,'LTYPE_O',trim(vacal),ltype_3,meshe_o % ltype)
       call memory_copy(memor_dom,'LTYPE_O',trim(vacal),lnnod_3,meshe_o % lnnod)
       call memory_copy(memor_dom,'LEINV_O',trim(vacal),leinv_3,meshe_o % leinv_loc)
                                                                          
       call memory_copy(memor_dom,'COORD_O',trim(vacal),coord_3,meshe_o % coord)
       call memory_copy(memor_dom,'LNINV_O',trim(vacal),lninv_3,meshe_o % lninv_loc)

       if( present(array_i) .and. present(array_o) ) then
          call memory_copy(memor_dom,'ARRAY_O',trim(vacal),array_3,array_o)
       end if
       do ifiel = 1,meshe_i % nfiel
          if( associated(meshe_i % xfiel(ifiel) % a) ) then
             call memory_copy(memor_dom,'XFIEL % A',trim(vacal),xfiel_3(ifiel) % a,meshe_o % xfiel(ifiel) % a)
          end if
       end do
       
       call htades( htable_3 )

    end if

  end subroutine meshes_glue_two_meshes
    
end module mod_meshes
!> @}
