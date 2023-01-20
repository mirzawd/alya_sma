!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @defgroup Graphs_Toolbox
!> Toolbox for matrix operations
!> @{
!> @name    ToolBox for graphs and renumbering
!> @file    mod_graphs.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for graphs and renumbering.
!> @details ToolBox for graphs and renumbering. Uses METIS_NodeND,
!>          (Node dissection) for renumbering
!
!-----------------------------------------------------------------------

module mod_graphs

  use def_kintyp_basic,       only : ip,rp,lg,i1p,i2p
  use def_kintyp_mesh,        only : mesh_type
  use def_kintyp_mesh,        only : mesh_type_basic
  use def_elmtyp,             only : NOHOL
  use mod_memory,             only : memory_alloca
  use mod_memory,             only : memory_deallo
  use mod_memory,             only : memory_resize
  use mod_memory,             only : memory_size
  use mod_elmgeo,             only : element_type
  use mod_maths,              only : maths_heap_sort
  use mod_htable,             only : htaini
  use mod_htable,             only : htaadd
  use mod_htable,             only : htares
  use mod_htable,             only : htades
  use mod_htable,             only : hash_t
  use mod_alya2metis,         only : alya2metis_METIS_NodeND
  use mod_strings,            only : integer_to_string
  use mod_optional_argument,  only : optional_argument
  use mod_graphs_basic,       only : graphs_permut_metis_postordering
  use mod_graphs_basic,       only : graphs_permut_metis_postordering_deallocate
  use mod_graphs_basic,       only : graphs_rengra
  use mod_graphs_basic,       only : graphs_comper
  use mod_graphs_basic,       only : graphs_iniper
  use mod_graphs_basic,       only : graphs_permut
  use mod_graphs_basic,       only : graphs_postorder
  use mod_graphs_basic,       only : graphs_subgra
  use mod_graphs_basic,       only : graphs_number_to_linked_list
  use mod_memory_config,      only : memory_config
  
  implicit none
  
  private 

  interface graphs_deallocate
     module procedure graphs_deallocate_ip14,&
          &           graphs_deallocate_ip18,&
          &           graphs_deallocate_rp1
  end interface graphs_deallocate

  interface graphs_poipoi
     module procedure graphs_poipoi_0,&
          &           graphs_poipoi_1,&
          &           graphs_poipoi_2
  end interface graphs_poipoi

  interface graphs_elepoi
     module procedure graphs_elepoi_0,&
          &           graphs_elepoi_1
  end interface graphs_elepoi

  interface graphs_eleele_faces  
     module procedure graphs_eleele_faces_0,&
          &           graphs_eleele_faces_1
  end interface graphs_eleele_faces

  interface graphs_edges
     module procedure graphs_edges_arrays,&
          &           graphs_edges_type
  end interface graphs_edges

  public :: graphs_list_faces
  public :: graphs_deallocate_list_faces
  public :: graphs_comper
  public :: graphs_dealep
  public :: graphs_deaper
  public :: graphs_eleele
  public :: graphs_eleele_deallocate
  public :: graphs_elepoi
  public :: graphs_elepoi_deallocate
  public :: graphs_iniper
  public :: graphs_permut
  public :: graphs_postorder
  public :: graphs_rengra
  public :: graphs_subgra
  public :: graphs_coloring
  public :: graphs_coloring_greedy
  public :: graphs_coloring_deallocate
  public :: graphs_deallocate
  public :: graphs_copyij
  public :: graphs_copyij_deallocate
  public :: graphs_permut_metis_postordering
  public :: graphs_permut_metis_postordering_deallocate
  public :: graphs_number_along_vector
  public :: graphs_edges
  public :: graphs_number_to_linked_list
  public :: graphs_csr_to_coo                    ! Convert CSR to COO format
  public :: graphs_csr_to_ell                    ! Convert CSR to ELL format
  public :: graphs_poipoi
  public :: graphs_poipoi_deallocate
  public :: graphs_element_element_graph
  public :: graphs_output                        ! Output a graph in GiD format
  public :: graphs_domlis
  public :: graphs_compress                      ! Remove nodes from a graph
  public :: graphs_find_edge                     ! Find an edge in the CSR format
  public :: graphs_dual_graph                    ! Construct the dual graph
  public :: graphs_dual_graph_deallocate         ! Construct the dual graph
  public :: graphs_color_graph                   ! Color a graph
  public :: graphs_color_graph_deallocate
  public :: graphs_eleele_faces                  ! Element graph using faces
  
contains

  !-----------------------------------------------------------------------
  !
  !> @brief   Find the edge position in a graph
  !> @details Find the edge position IZ in a graph (IA,JA)
  !> @date    17/10/2017     
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine graphs_find_edge(ii,jj,ia,ja,iz,what)

    integer(ip),          intent(in)            :: ii       !< Vertex 1
    integer(ip),          intent(in)            :: jj       !< Vertex 2
    integer(ip), pointer, intent(in)            :: ia(:)    !< Input graph
    integer(ip), pointer, intent(in)            :: ja(:)    !< Input graph
    integer(ip),          intent(out)           :: iz       !< Edge number
    character(*),         intent(in), optional  :: what     !< Option
    integer(ip)                                 :: kk,i1,i2
    logical(lg)                                 :: kfl_directed

    if( ii == 0 .or. jj == 0 ) then
       iz = 0
    else
       kfl_directed = .true.
       if( present(what) ) then
          if(      trim(what) == 'DIRECTED' ) then
             kfl_directed = .true.
          else if( trim(what) == 'INDIRECTED' ) then
             kfl_directed = .false.
          else
             call runend('GRAPHS_FIND_EDGE: UNKNOWN OPTION')
          end if
       end if

       if( kfl_directed ) then
          i1 = ii
          i2 = jj
       else
          i1 = min(ii,jj)
          i2 = max(ii,jj)
       end if

       iz = 0
       kk = ia(i1)
       do while( kk <= ia(i1+1)-1 )
          if( ja(kk) == i2 ) then
             iz = kk
             return
          end if
          kk = kk + 1
       end do
    end if

  end subroutine graphs_find_edge

  !-----------------------------------------------------------------------
  !
  !> @brief   Remove nodes from a graph
  !> @details Remove nodes from a graph using a mask
  !>      
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine graphs_compress(nn_in,ia_in,ja_in,lmask,nn_out,ia_out,ja_out,invpr,memor)

    integer(ip),          intent(in)            :: nn_in     !< Number of nodes of the graph
    integer(ip), pointer, intent(in)            :: ia_in(:)  !< Input graph
    integer(ip), pointer, intent(in)            :: ja_in(:)  !< Input graph
    integer(ip), pointer, intent(in)            :: lmask(:)  !< Mask
    integer(ip),          intent(out)           :: nn_out    !< Number of output nodes
    integer(ip), pointer, intent(inout)           :: ia_out(:) !< Output graph
    integer(ip), pointer, intent(inout)           :: ja_out(:) !< Output graph
    integer(ip), pointer, intent(inout), optional :: invpr(:)  !< Inverse permutation OLD = INVPR(NEW)
    integer(8),           intent(inout)           :: memor(2)
    integer(ip), pointer                        :: permr(:)
    integer(ip)                                 :: nz_out,jz
    integer(ip)                                 :: ii,jj,iz

    nullify(permr)
    call memory_alloca(memor,'PERMR','graphs_compress',permr,nn_in)

    nn_out = 0
    do ii = 1,nn_in
       if( lmask(ii) >= 0 ) then
          nn_out = nn_out + 1
          permr(ii) = nn_out
       end if
    end do
    if( present(invpr) ) then
       call memory_alloca(memor,'INVPR','graphs_compress',invpr,nn_out)
       nn_out = 0
       do ii = 1,nn_in
          if( lmask(ii) >= 0 ) then
             nn_out = nn_out + 1
             invpr(nn_out) = ii
          end if
       end do
    end if

    nz_out = 0
    do ii = 1,nn_in
       if( permr(ii) > 0 ) then 
          do iz = ia_in(ii),ia_in(ii+1)-1
             if( permr(ja_in(iz)) > 0 ) nz_out = nz_out + 1
          end do
       end if
    end do
    call memory_alloca(memor,'IA_OUT','graphs_compress',ia_out,nn_out+1_ip)
    call memory_alloca(memor,'JA_OUT','graphs_compress',ja_out,nz_out)

    nz_out = 0
    if( nn_out > 0 ) ia_out(1) = 1
    nn_out = 0
    do ii = 1,nn_in
       if( permr(ii) > 0 ) then 
          jz = 0
          nn_out = nn_out + 1
          do iz = ia_in(ii),ia_in(ii+1)-1
             jj = ja_in(iz)
             if( permr(jj) > 0 ) then
                jz     = jz + 1
                nz_out = nz_out + 1
                ja_out(nz_out) = permr(jj) 
             end if
          end do
          ia_out(nn_out+1) = ia_out(nn_out) + jz
       end if
    end do

    call memory_deallo(memor,'PERMR','graphs_compress',permr)

  end subroutine graphs_compress

  !-----------------------------------------------------------------------
  !
  !> @brief   Compute a node-node graph
  !> @details Compute a node-node graph
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine graphs_poipoi_deallocate(ia,ja,memor)
    
    integer(ip), pointer, intent(inout), optional :: ia(:)                    !< Linked list of (node-node) pointer 
    integer(ip), pointer, intent(inout), optional :: ja(:)                    !< Linked list of (node-node) elements
    integer(8),           intent(inout)           :: memor(2)
    
    if( present(ia) ) then
       call memory_deallo(memor,'IA','graphs_poipoi_deallocate',ia)
    end if
    if( present(ja) ) then
       call memory_deallo(memor,'JA','graphs_poipoi_deallocate',ja)
    end if
    
  end subroutine graphs_poipoi_deallocate
  
  subroutine graphs_poipoi_0(&
       meshe,ia,ja,message,IA_NAME,JA_NAME,INCLUDE_HALOS,memor)
    
    type(mesh_type),      intent(inout)         :: meshe                    !< Mesh type
    integer(ip), pointer, intent(inout)         :: ia(:)                    !< Linked list of (node-node) pointer 
    integer(ip), pointer, intent(inout)         :: ja(:)                    !< Linked list of (node-node) elements
    character(len=*),     intent(in),  optional :: message                  !< Message for options
    character(len=*),     intent(in),  optional :: IA_NAME                  !< Graph name
    character(len=*),     intent(in),  optional :: JA_NAME                  !< Graph name
    logical(lg),          intent(in),  optional :: INCLUDE_HALOS            !< If halo should be included
    integer(8),           intent(inout)         :: memor(2)
    integer(ip)                                 :: bandw                    
    real(rp)                                    :: profi                    
    integer(ip), pointer                        :: pelpo(:)                 
    integer(ip), pointer                        :: lelpo(:)                 
    integer(ip)                                 :: mepoi
    integer(ip)                                 :: pelty,ielem
    integer(ip), pointer                        :: lnnod_loc(:)
    logical(lg)                                 :: if_halos
    integer(ip)                                 :: nelem_loc
    
    nullify(lnnod_loc)
    nullify(pelpo)
    nullify(lelpo)

    if( present(INCLUDE_HALOS) ) then
       if_halos = INCLUDE_HALOS      
    else
       if_halos = .false.
    end if
    if( if_halos ) then        
       nelem_loc = meshe % nelem_2
    else
       nelem_loc = meshe % nelem
    end if
    
    if( .not. associated(meshe % lnnod) ) then
       call memory_alloca(memor,'LNNOD_LOC','graphs_poipoi',lnnod_loc,nelem_loc)
       do ielem = 1,nelem_loc
          pelty = abs(meshe % ltype(ielem))
          lnnod_loc(ielem) = element_type(pelty) % number_nodes
       end do
    else
       lnnod_loc => meshe % lnnod
    end if
    
    if( if_halos ) then        
       call graphs_poipoi_go(&
            meshe % npoin_2,meshe % nelem_2,meshe % mnode,meshe % lnods,lnnod_loc,meshe % ltype,ia,ja,&
            bandw,profi,pelpo,lelpo,mepoi,message,&
            IA_NAME,JA_NAME,memor=memor)
    else
       call graphs_poipoi_go(&
            meshe % npoin,meshe % nelem,meshe % mnode,meshe % lnods,lnnod_loc,meshe % ltype,ia,ja,&
            bandw,profi,pelpo,lelpo,mepoi,message,&
            IA_NAME,JA_NAME,memor=memor)
    end if
    
    call memory_deallo(memor,'LELPO','graphs_poipoi',lelpo)
    call memory_deallo(memor,'PELPO','graphs_poipoi',pelpo)
    
    if( .not. associated(meshe % lnnod) ) then
       call memory_deallo(memor,'LNNOD_LOC','graphs_poipoi',lnnod_loc)
    end if
    
  end subroutine graphs_poipoi_0

  subroutine graphs_poipoi_1(&
       npoin,nelem,mnode,lnods,lnnod,ltype,ia,ja,message,&
       IA_NAME,JA_NAME,memor)
    implicit none
    integer(ip),          intent(in)            :: npoin                    !< Number of nodes
    integer(ip),          intent(in)            :: nelem                    !< Number of elements
    integer(ip),          intent(in)            :: mnode                    !< Max. number of nodes per element
    integer(ip), pointer, intent(in)            :: lnods(:,:)               !< Connectivity array
    integer(ip), pointer, intent(in)            :: lnnod(:)                 !< Array of number of element nodes
    integer(ip), pointer, intent(in)            :: ltype(:)                 !< Array of element types
    integer(ip), pointer, intent(inout)         :: ia(:)                    !< Linked list of (node-node) pointer 
    integer(ip), pointer, intent(inout)         :: ja(:)                    !< Linked list of (node-node) elements
    character(len=*),     intent(in),  optional :: message                  !< Message for options
    character(len=*),     intent(in),  optional :: IA_NAME                  !< Graph name
    character(len=*),     intent(in),  optional :: JA_NAME                  !< Graph name
    integer(8),           intent(inout)         :: memor(2)
    integer(ip)                                 :: bandw                    
    real(rp)                                    :: profi                    
    integer(ip), pointer                        :: pelpo(:)                 
    integer(ip), pointer                        :: lelpo(:)                 
    integer(ip)                                 :: mepoi
    integer(ip)                                 :: pelty,ielem
    integer(ip), pointer                        :: lnnod_loc(:)

    nullify(lnnod_loc)
    nullify(pelpo)
    nullify(lelpo)

    if( .not. associated(lnnod) ) then
       call memory_alloca(memor,'LNNOD_LOC','graphs_poipoi',lnnod_loc,nelem)
       do ielem = 1,nelem
          pelty = abs(ltype(ielem))
          lnnod_loc(ielem) = element_type(pelty) % number_nodes
       end do
    else
       lnnod_loc => lnnod
    end if
    call graphs_poipoi_go(&
         npoin,nelem,mnode,lnods,lnnod_loc,ltype,ia,ja,&
         bandw,profi,pelpo,lelpo,mepoi,message,&
         IA_NAME,JA_NAME,memor=memor)

    call memory_deallo(memor,'LELPO','graphs_poipoi',lelpo)
    call memory_deallo(memor,'PELPO','graphs_poipoi',pelpo)
    
    if( .not. associated(lnnod) ) then
       call memory_deallo(memor,'LNNOD_LOC','graphs_poipoi',lnnod_loc)
    end if
    
  end subroutine graphs_poipoi_1

  subroutine graphs_poipoi_2(&
       npoin,nelem,mnode,lnods,lnnod,ltype,ia,ja,&
       bandw,profi,pelpo,lelpo,mepoi,message,&
       IA_NAME,JA_NAME,memor)
    implicit none
    integer(ip),          intent(in)              :: npoin                    !< Number of nodes
    integer(ip),          intent(in)              :: nelem                    !< Number of elements
    integer(ip),          intent(in)              :: mnode                    !< Max. number of nodes per element
    integer(ip), pointer, intent(in)              :: lnods(:,:)               !< Connectivity array
    integer(ip), pointer, intent(in)              :: lnnod(:)                 !< Array of number of element nodes
    integer(ip), pointer, intent(in)              :: ltype(:)                 !< Array of element types
    integer(ip), pointer, intent(inout)           :: ia(:)                    !< Linked list of (node-node) pointer 
    integer(ip), pointer, intent(inout)           :: ja(:)                    !< Linked list of (node-node) elements
    integer(ip),          intent(out)             :: bandw                    !< Bandwidt
    real(rp),             intent(out)             :: profi                    !< Profile
    integer(ip), pointer, intent(inout), optional :: pelpo(:)                 !< Linked list of (element-node) pointer 
    integer(ip), pointer, intent(inout), optional :: lelpo(:)                 !< Linked list of (element-node) elements
    integer(ip),          intent(out),   optional :: mepoi                    !< Max number of element per node
    character(len=*),     intent(in),    optional :: message                  !< Message for options
    character(len=*),     intent(in),    optional :: IA_NAME                  !< Graph name
    character(len=*),     intent(in),    optional :: JA_NAME                  !< Graph name
    integer(8),           intent(inout)           :: memor(2)
    integer(ip)                                   :: pelty,ielem
    integer(ip), pointer                          :: lnnod_loc(:)

    nullify(lnnod_loc)

     if( .not. associated(lnnod) ) then
       call memory_alloca(memor,'LNNOD_LOC','graphs_poipoi',lnnod_loc,nelem)
       do ielem = 1,nelem
          pelty = abs(ltype(ielem))
          lnnod_loc(ielem) = element_type(pelty) % number_nodes
       end do
    else
       lnnod_loc => lnnod
    end if
    
    call graphs_poipoi_go(&
         npoin,nelem,mnode,lnods,lnnod_loc,ltype,ia,ja,&
         bandw,profi,pelpo,lelpo,mepoi,message,&
         IA_NAME,JA_NAME,memor=memor)
    
    if( .not. associated(lnnod) ) then
       call memory_deallo(memor,'LNNOD_LOC','graphs_poipoi',lnnod_loc)
    end if

  end subroutine graphs_poipoi_2

  subroutine graphs_poipoi_go(&
       npoin,nelem,mnode,lnods,lnnod,ltype,ia,ja,&
       bandw,profi,pelpo,lelpo,mepoi,message,&
       IA_NAME,JA_NAME,memor)
    implicit none
    integer(ip),          intent(in)              :: npoin                    !< Number of nodes
    integer(ip),          intent(in)              :: nelem                    !< Number of elements
    integer(ip),          intent(in)              :: mnode                    !< Max. number of nodes per element
    integer(ip),          intent(in)              :: lnods(mnode,nelem)       !< Connectivity array
    integer(ip),          intent(in)              :: lnnod(nelem)             !< Array of number of element nodes
    integer(ip),          intent(in)              :: ltype(nelem)             !< Array of element types
    integer(ip), pointer, intent(inout)           :: ia(:)                    !< Linked list of (node-node) pointer 
    integer(ip), pointer, intent(inout)           :: ja(:)                    !< Linked list of (node-node) elements
    integer(ip),          intent(out),   optional :: bandw                    !< Bandwidt
    real(rp),             intent(out),   optional :: profi                    !< Profile
    integer(ip), pointer, intent(inout), optional :: pelpo(:)                 !< Linked list of (element-node) pointer 
    integer(ip), pointer, intent(inout), optional :: lelpo(:)                 !< Linked list of (element-node) elements
    integer(ip),          intent(out),   optional :: mepoi                    !< Max number of element per node
    character(len=*),     intent(in),    optional :: message                  !< Message for options
    character(len=*),     intent(in),    optional :: IA_NAME                  !< Graph name
    character(len=*),     intent(in),    optional :: JA_NAME                  !< Graph name
    integer(8),           intent(inout)           :: memor(2)
    integer(ip)                                   :: ipoin,ielem,jelem,icoef
    integer(ip)                                   :: izdom,ncoef,nlelp,mtouc
    integer(ip)                                   :: lsize,mpopo,nz
    integer(ip)                                   :: which_graph,mepo2
    logical(lg)                                   :: only_edges
    integer(ip), pointer                          :: pelp2(:)                 
    integer(ip), pointer                          :: lelp2(:)                 
    integer(ip), pointer                          :: lista(:)
    logical(lg), pointer                          :: touch(:)
    character(5)                                  :: my_ia_name
    character(5)                                  :: my_ja_name

    if( present(IA_NAME) ) then
       my_ia_name = trim(IA_NAME)
    else
       my_ia_name = 'IA'
    end if
    if( present(JA_NAME) ) then
       my_ja_name = trim(JA_NAME)
    else
       my_ja_name = 'JA'
    end if
    !
    ! Initialization
    !
    nullify(lista)
    nullify(touch)
    nullify(pelp2)
    nullify(lelp2)
    which_graph = 0                                    ! 0:all, 1:all without diagonal, -1: inferior part, -2: edges
    only_edges  = .false.
    if( present(message) ) then
       if(      message == 'REMOVE DIAGONAL' ) then
          which_graph =  1                             !  1: Remvole diagonal
       else if( message == 'LOWER PART' ) then  
          which_graph = -1                             ! -1: only lower part
       else if( message == 'EDGES' ) then
          call runend('MOD_GRAPHS: NO LONGER SUPPORTED')
       else if( message == 'SQUARE REMOVE DIAGONAL' ) then  
          which_graph =  2                             !  2: only up to npoin
       end if
    end if
    !
    ! Check size
    !
    !if( size(lnods,1,KIND=ip) /= mnode ) call runend('WRONG LNODS SIZE')
    !
    ! Compute node-element connectivity if not already computed: PELP2, LELP2
    !
    if( present(pelpo) .and. present(lelpo) ) then
       if( .not. associated(pelpo) ) then
#ifdef __PGI
          call graphs_elepoi(npoin,nelem,mnode,lnods,lnnod,mepo2,pelpo,lelpo)
#else
          call graphs_elepoi(npoin,nelem,mnode,lnods,lnnod,mepo2,pelpo,lelpo,memor=memor)
#endif
       end if
       pelp2 => pelpo
       lelp2 => lelpo
    else
#ifdef __PGI
       call graphs_elepoi(npoin,nelem,mnode,lnods,lnnod,mepo2,pelp2,lelp2)
#else
       call graphs_elepoi(npoin,nelem,mnode,lnods,lnnod,mepo2,pelp2,lelp2,memor=memor)
#endif
    end if
    if( present(mepoi) ) mepoi = mepo2
    ! 
    ! Compute node-node graph: JA and IA
    !
    call memory_alloca(memor,trim(my_ia_name),'graphs_poipoi',ia,npoin+1)

    if( memory_config%high_memory ) then

       !-------------------------------------------------------------------
       !
       ! Strategy 1: slow but does not require lots of memory
       !
       !-------------------------------------------------------------------

       mtouc = 0
       do ipoin = 1,npoin
          mtouc = max(mtouc,pelp2(ipoin+1)-pelp2(ipoin))
       end do
       mtouc = mtouc * mnode
       nz    = 0

       call memory_alloca(memor,'TOUCH','graphs_poipoi',touch,mtouc)

       if( which_graph == 0 ) then
          !
          !*OMP  PARALLEL DO                                    &
          !*OMP  SCHEDULE     ( STATIC )                        & 
          !*OMP  DEFAULT      ( NONE )                          &
          !*OMP  PRIVATE      ( ipoin,nlelp,ncoef,icoef,touch ) &
          !*OMP  SHARED       ( npoin,pelp2,mnode,lnods,lnnod,  &
          !*OMP                 ltype,lelp2,only_edges )        &
          !*OMP  REDUCTION    ( +:nz )                    
          !
          do ipoin = 1,npoin
             nlelp = pelp2(ipoin+1) - pelp2(ipoin)
             ncoef = nlelp * mnode
             do icoef = 1,ncoef          
                touch(icoef) = .false.
             end do
             call graphs_nzecof(mnode,lnods,lnnod,ltype,nlelp,ncoef,nz,&
                  lelp2(pelp2(ipoin):pelp2(ipoin+1)-1),&
                  touch,&
                  0_ip,&
                  only_edges)
          end do
          !*OMP END PARALLEL DO
          !
          ! Construct the array of indexes
          ! 
          call memory_alloca(memor,trim(my_ja_name),'graphs_poipoi',ja,nz)
          izdom = 1
          do ipoin = 1,npoin
             nlelp = pelp2(ipoin+1) - pelp2(ipoin)
             ncoef = nlelp * mnode
             do icoef = 1,ncoef          
                touch(icoef) = .false.
             end do
             call graphs_arrind(mnode,lnods,lnnod,ltype,nlelp,ncoef,&
                  lelp2(pelp2(ipoin):pelp2(ipoin+1)-1),touch,izdom,ipoin,ia,ja,0_ip,&
                  only_edges)
          end do

       else if( which_graph == 2 ) then
          !
          ! Only include nodes up to npoin and remove diagonal
          !
          do ipoin = 1,npoin
             nlelp = pelp2(ipoin+1) - pelp2(ipoin)
             ncoef = nlelp * mnode
             do icoef = 1,ncoef          
                touch(icoef) = .false.
             end do
             call graphs_nzecof(mnode,lnods,lnnod,ltype,nlelp,ncoef,nz,&
                  lelp2(pelp2(ipoin):pelp2(ipoin+1)-1),&
                  touch,&
                  ipoin,&
                  only_edges,npoin)
          end do
          !
          ! Construct the array of indexes
          ! 
          call memory_alloca(memor,trim(my_ja_name),'graphs_poipoi',ja,nz)
          izdom = 1
          do ipoin = 1,npoin
             nlelp = pelp2(ipoin+1) - pelp2(ipoin)
             ncoef = nlelp * mnode
             do icoef = 1,ncoef          
                touch(icoef) = .false.
             end do
             call graphs_arrind(mnode,lnods,lnnod,ltype,nlelp,ncoef,&
                  lelp2(pelp2(ipoin):pelp2(ipoin+1)-1),touch,izdom,ipoin,ia,ja,ipoin,&
                  only_edges,npoin)
          end do

       else if( which_graph == 1 ) then

          do ipoin = 1,npoin
             nlelp = pelp2(ipoin+1) - pelp2(ipoin)
             ncoef = nlelp * mnode
             do icoef = 1,ncoef          
                touch(icoef) = .false.
             end do
             call graphs_nzecof(mnode,lnods,lnnod,ltype,nlelp,ncoef,nz,&
                  lelp2(pelp2(ipoin):pelp2(ipoin+1)-1),&
                  touch,&
                  ipoin,&
                  only_edges)
          end do
          !
          ! Construct the array of indexes
          ! 
          call memory_alloca(memor,trim(my_ja_name),'graphs_poipoi',ja,nz)
          izdom = 1
          do ipoin = 1,npoin
             nlelp = pelp2(ipoin+1) - pelp2(ipoin)
             ncoef = nlelp * mnode
             do icoef = 1,ncoef          
                touch(icoef) = .false.
             end do
             call graphs_arrind(mnode,lnods,lnnod,ltype,nlelp,ncoef,&
                  lelp2(pelp2(ipoin):pelp2(ipoin+1)-1),touch,izdom,ipoin,ia,ja,ipoin,&
                  only_edges)
          end do

       else if( which_graph == -1 ) then

          do ipoin = 1,npoin
             nlelp = pelp2(ipoin+1) - pelp2(ipoin)
             ncoef = nlelp * mnode
             do icoef = 1,ncoef          
                touch(icoef) = .false.
             end do
             call graphs_nzecof(mnode,lnods,lnnod,ltype,nlelp,ncoef,nz,&
                  lelp2(pelp2(ipoin):pelp2(ipoin+1)-1),touch,-ipoin,only_edges)
          end do
          !
          ! Construct the array of indexes
          ! 
          call memory_alloca(memor,trim(my_ja_name),'graphs_poipoi',ja,nz)
          izdom = 1
          do ipoin = 1,npoin
             nlelp = pelp2(ipoin+1) - pelp2(ipoin)
             ncoef = nlelp * mnode
             do icoef = 1,ncoef          
                touch(icoef) = .false.
             end do
             call graphs_arrind(mnode,lnods,lnnod,ltype,nlelp,ncoef,&
                  lelp2(pelp2(ipoin):pelp2(ipoin+1)-1),touch,izdom,ipoin,ia,ja,-ipoin,&
                  only_edges)
          end do

       end if

       ia(npoin+1) = nz + 1
       call memory_deallo(memor,'TOUCH','graphs_poipoi',touch)

    else

       !-------------------------------------------------------------------
       !
       ! Strategy 2: quick but requires lots of memory
       ! 
       !-------------------------------------------------------------------
       mpopo = 0_ip
       do ielem = 1,nelem
          mpopo = mpopo + lnnod(ielem)*lnnod(ielem)
       end do

       call memory_alloca(memor,'LISTA','graphs_poipoi',lista,mpopo)
       !
       ! Construct the array of indexes
       !     
       ia(1) = 1
       if( which_graph == 0 ) then
          do ipoin = 1, npoin
             lsize = 0
             do ielem = pelp2(ipoin), pelp2(ipoin+1)-1
                jelem = lelp2(ielem)
                call graphs_mergli( lista(ia(ipoin):), lsize, lnnod(jelem), &
                     lnods(:,jelem), ltype(jelem), 0_ip , only_edges)
             enddo
             ia(ipoin+1) = ia(ipoin) + lsize
          end do

       else if( which_graph == 1 ) then
          do ipoin = 1, npoin
             lsize = 0
             do ielem = pelp2(ipoin), pelp2(ipoin+1)-1
                jelem = lelp2(ielem)
                call graphs_mergli( lista(ia(ipoin):), lsize, lnnod(jelem), &
                     lnods(:,jelem), ltype(jelem), ipoin , only_edges)
             enddo
             ia(ipoin+1) = ia(ipoin) + lsize
          end do

       else if( which_graph == 2 ) then
          do ipoin = 1, npoin
             lsize = 0
             do ielem = pelp2(ipoin), pelp2(ipoin+1)-1
                jelem = lelp2(ielem)
                call graphs_mergli( lista(ia(ipoin):), lsize, lnnod(jelem), &
                     lnods(:,jelem), ltype(jelem), ipoin , only_edges, npoin)
             enddo
             ia(ipoin+1) = ia(ipoin) + lsize
          end do

       else if( which_graph == -1 ) then
          do ipoin = 1, npoin
             lsize = 0
             do ielem = pelp2(ipoin), pelp2(ipoin+1)-1
                jelem = lelp2(ielem)
                call graphs_mergli( lista(ia(ipoin):), lsize, lnnod(jelem), &
                     lnods(:,jelem), ltype(jelem), -ipoin , only_edges)
             enddo
             ia(ipoin+1) = ia(ipoin) + lsize
          end do
       end if

       nz = ia(npoin+1)-1
       call memory_alloca(memor,trim(my_ja_name),'graphs_poipoi',ja,nz)  
       do ipoin=1, nz
          ja(ipoin) = lista(ipoin)
       enddo
       call memory_deallo(memor,'LISTA','graphs_poipoi',lista)
    end if
    !
    ! Deallocate memory
    !
    if( present(pelpo) .and. present(lelpo) ) then
       continue
    else
       call memory_deallo(memor,'LELP2','graphs_poipoi',lelp2)
       call memory_deallo(memor,'PELP2','graphs_poipoi',pelp2)
    end if

    !-------------------------------------------------------------------
    !
    ! Order graph: the tst LSIZE > 0 is necessary to treat nodes
    ! without graphs. Example: master nodes comming from a neighboring
    ! subdomain.
    ! 
    !-------------------------------------------------------------------

    do ipoin = 1,npoin
       lsize = ia(ipoin+1) - ia(ipoin)
       if( lsize > 0 ) call heapsorti1(2_ip,lsize,ja(ia(ipoin)))
    end do

    !-------------------------------------------------------------------
    !
    ! Compute profile and bandwidth
    ! 
    !-------------------------------------------------------------------

    if( present(bandw) .and. present(profi) ) then
       call graphs_gtband(npoin,ja,ia,bandw,profi)
    end if

  end subroutine graphs_poipoi_go

  !-----------------------------------------------------------------------
  !
  !> @brief   Compute the node/element connectivity arrays.
  !> @details Compute the node/element connectivity arrays.
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine graphs_elepoi_deallocate(pelpo,lelpo,PELPO_NAME,LELPO_NAME,memor)

    integer(ip), pointer, intent(inout), optional :: pelpo(:)                 !< Linked list of (element-node) pointer
    integer(ip), pointer, intent(inout), optional :: lelpo(:)                 !< Linked list of (element-node) elements
    character(len=*),     intent(in),    optional :: PELPO_NAME
    character(len=*),     intent(in),    optional :: LELPO_NAME
    integer(8),           intent(inout), optional :: memor(2)
    integer(8)                                    :: memor_loc(2)

    memor_loc = optional_argument((/0_8,0_8/),memor)
    if( present(PELPO) ) then
       if( present(PELPO_NAME) ) then
          call memory_deallo(memor_loc,trim(PELPO_NAME),'graphs_elepoi_deallocate',pelpo)
       else
          call memory_deallo(memor_loc,'PELPO','graphs_elepoi_deallocate',pelpo)
       end if
    end if
    if( present(LELPO) ) then
       if( present(LELPO_NAME) ) then
          call memory_deallo(memor_loc,trim(LELPO_NAME),'graphs_elepoi_deallocate',lelpo)
       else
          call memory_deallo(memor_loc,'LELPO','graphs_elepoi_deallocate',lelpo)
       end if
    end if
    if( present(memor) ) memor = memor_loc

  end subroutine graphs_elepoi_deallocate
  
  subroutine graphs_elepoi_0(meshe,mepoi,pelpo,lelpo,PELPO_NAME,LELPO_NAME,memor)
    
    type(mesh_type_basic),          intent(in)    :: meshe          !< Mesh type
    integer(ip),          pointer,  intent(inout) :: pelpo(:)       !< Linked list of (element-node) pointer
    integer(ip),          pointer,  intent(inout) :: lelpo(:)       !< Linked list of (element-node) elements
    integer(ip),                    intent(out)   :: mepoi
    character(len=*),     optional, intent(in)    :: PELPO_NAME
    character(len=*),     optional, intent(in)    :: LELPO_NAME
    integer(8),           optional, intent(inout) :: memor(2)
    integer(ip)                                   :: ielem,pelty
    integer(ip), pointer                          :: lnnod_loc(:)
    integer(8)                                    :: memor_loc(2)

    if( meshe % nelem > 0 ) then
       nullify(lnnod_loc)
       call memory_alloca(memor_loc,'LNNOD_LOC','graphs_poipoi',lnnod_loc,meshe % nelem)
       do ielem = 1,meshe % nelem
          pelty = abs(meshe % ltype(ielem))
          lnnod_loc(ielem) = element_type(pelty) % number_nodes
       end do
       call graphs_elepoi_1(&
            meshe % npoin,meshe % nelem,meshe % mnode,meshe % lnods,lnnod_loc,mepoi,&
            pelpo,lelpo,PELPO_NAME=PELPO_NAME,LELPO_NAME=LELPO_NAME,memor=memor)       
       call memory_deallo(memor_loc,'LNNOD_LOC','graphs_poipoi',lnnod_loc)
    else
       mepoi = 0
    end if

  end subroutine graphs_elepoi_0
  
  subroutine graphs_elepoi_1(&
       npoin,nelem,mnode,lnods,lnnod,mepoi,pelpo,lelpo,leper,lninv,&
       PELPO_NAME,LELPO_NAME,memor)

    integer(ip),          intent(in)           :: npoin                    !< Number of nodes
    integer(ip),          intent(in)           :: nelem                    !< Number of elements
    integer(ip),          intent(in)           :: mnode                    !< Max. number of nodes per element
    integer(ip),          intent(in)           :: lnods(mnode,*)           !< Connectivity array
    integer(ip),          intent(in)           :: lnnod(*)                 !< Array of number of element nodes
    integer(ip),          intent(out)          :: mepoi                    !< Max. number of element per node
    integer(ip), pointer, intent(inout)        :: pelpo(:)                 !< Linked list of (element-node) pointer
    integer(ip), pointer, intent(inout)        :: lelpo(:)                 !< Linked list of (element-node) elements
    integer(ip), pointer, intent(in), optional :: leper(:)                 !< Element permutation array
    integer(ip), pointer, intent(in), optional :: lninv(:)                 !< Nodal permutation array
    character(len=*),     intent(in), optional :: PELPO_NAME
    character(len=*),     intent(in), optional :: LELPO_NAME
    integer(8),           intent(inout), optional :: memor(2)
    integer(ip), pointer                       :: nepoi(:)
    integer(ip)                                :: inode,ipoin,ielem
    integer(ip)                                :: jelem,jpoin,nlelp
    logical(lg)                                :: lperm,check_npoin_bound
    integer(8)                                 :: memor_loc(2)

    memor_loc = optional_argument((/0_8,0_8/),memor)

    if( nelem == 0 ) then
       mepoi = 0
    else
       !
       ! Initialization
       !
       nullify(nepoi)
       !
       ! Permutation?
       !
       lperm = present(leper) .and. present(lninv)
       !
       ! Check bounds. Useful for example when we want graph only up to
       ! a number of nodes npoin lower that the total number of nodes of the mesh
       !
       check_npoin_bound = .false.
       if( maxval(lnods(1:mnode,1:nelem)) > npoin ) check_npoin_bound = .true.
       !
       !
       ! Allocate memory for NEPOI and compute it
       !
       call memory_alloca(memor_loc,'NEPOI','elepoi',nepoi,npoin)
       !
       ! NEPOI: Number of element per node
       !
       if( check_npoin_bound ) then
          if( lperm ) then
             do ielem = 1,nelem
                jelem = leper(ielem)
                do inode = 1,lnnod(jelem)
                   ipoin = lnods(inode,jelem)
                   jpoin = lninv(ipoin)
                   if( jpoin <= npoin ) nepoi(jpoin) = nepoi(jpoin) + 1
                end do
             end do
          else
             do ielem = 1,nelem
                do inode = 1,lnnod(ielem)
                   ipoin = lnods(inode,ielem)
                   if( ipoin <= npoin ) nepoi(ipoin) = nepoi(ipoin) + 1
                end do
             end do
          end if
       else
          if( lperm ) then
             do ielem = 1,nelem
                jelem = leper(ielem)
                do inode = 1,lnnod(jelem)
                   ipoin = lnods(inode,jelem)
                   jpoin = lninv(ipoin)
                   nepoi(jpoin) = nepoi(jpoin) + 1
                end do
             end do
          else
             do ielem = 1,nelem
                do inode = 1,lnnod(ielem)
                   ipoin = lnods(inode,ielem)
                   nepoi(ipoin) = nepoi(ipoin) + 1
                end do
             end do
          end if
       end if
       !
       ! Allocate memory for PELPO and compute it
       !
       if( present(PELPO_NAME) ) then
          call memory_alloca(memor_loc,trim(PELPO_NAME),'elepoi',pelpo,npoin+1_ip)
       else
          call memory_alloca(memor_loc,'PELPO','elepoi',pelpo,npoin+1_ip)
       end if
       pelpo(1) = 1
       do ipoin = 1,npoin
          pelpo(ipoin+1) = pelpo(ipoin) + nepoi(ipoin)
       end do
       !
       ! Allocate memory for LELPO and construct the list
       !
       nlelp = pelpo(npoin+1)
       if( present(LELPO_NAME) ) then
          call memory_alloca(memor_loc,trim(LELPO_NAME),'elepoi',lelpo,nlelp)
       else
          call memory_alloca(memor_loc,'LELPO','elepoi',lelpo,nlelp)
       end if

       if( check_npoin_bound ) then
          if( lperm ) then
             do ielem = 1,nelem
                jelem = leper(ielem)
                do inode = 1,lnnod(jelem)
                   ipoin = lnods(inode,jelem)
                   jpoin = lninv(ipoin)
                   if( jpoin <= npoin ) then
                      lelpo(pelpo(jpoin)) = ielem
                      pelpo(jpoin) = pelpo(jpoin)+1
                   end if
                end do
             end do
          else
             do ielem = 1,nelem
                do inode = 1,lnnod(ielem)
                   ipoin = lnods(inode,ielem)
                   if( ipoin <= npoin ) then
                      lelpo(pelpo(ipoin)) = ielem
                      pelpo(ipoin) = pelpo(ipoin)+1
                   end if
                end do
             end do
          end if
       else
          if( lperm ) then
             do ielem = 1,nelem
                jelem = leper(ielem)
                do inode = 1,lnnod(jelem)
                   ipoin = lnods(inode,jelem)
                   jpoin = lninv(ipoin)
                   lelpo(pelpo(jpoin)) = ielem
                   pelpo(jpoin) = pelpo(jpoin)+1
                end do
             end do
          else
             do ielem = 1,nelem
                do inode = 1,lnnod(ielem)
                   ipoin = lnods(inode,ielem)
                   lelpo(pelpo(ipoin)) = ielem
                   pelpo(ipoin) = pelpo(ipoin)+1
                end do
             end do
          end if
       end if
       !
       ! Recompute PELPO and maximum number of element neighbors MEPOI
       !
       pelpo(1) =  1
       mepoi    = -1
       do ipoin = 1,npoin
          pelpo(ipoin+1) = pelpo(ipoin) + nepoi(ipoin)
          mepoi = max(mepoi,nepoi(ipoin))
       end do
       !
       ! Deallocate memory for temporary node/element connectivity
       !
       call memory_deallo(memor_loc,'NEPOI','elepoi',nepoi)

       if( present(memor) ) memor = memor_loc

    end if

  end subroutine graphs_elepoi_1

  !-----------------------------------------------------------------------
  !
  !> @brief   Compute a graph
  !> @details Compute the element-element graph using nodal connectivity
  !>    For example, given the following mesh, element 5 will have the
  !>    following neighbors: \n
  !>    @verbatim
  !>    +---+---+---+
  !>    | 1 | 2 | 3 |  1,2,3,4,6,7,8,9
  !>    +---+---+---+
  !>    | 4 | 5 | 6 |
  !>    +---+---+---+
  !>    | 7 | 8 | 9 |
  !>    +---+---+---+
  !>    @endverbatim
  !>    \n
  !>    Working arrays:\n
  !>    NEPOI:      Number of occurencies of a node in connectivity\n
  !>                Node ipoin belongs to nepoi(ipoin) elements\n
  !>    PELPO:      Pointer to node connectivities arrays lelpo\n
  !>    LELPO:      List of node to element connectivities\n
  !>    MELEL:      Maximum number of neighbors=max(nepoi)*mnode\n
  !>                                                  \n
  !>    Example:                                      \n
  !>    @verbatim
  !>                 1---2---3---4---5  
  !>                 |   |   |   |   | 
  !>                 6---7---8---9--10 
  !>                 |   |   |   |   |
  !>                11--12--13--14--15  
  !>    @endverbatim
  !>                                                  \n
  !>    @verbatim
  !>    nedge:  22                               
  !>    element #:  1  2  3  4  5 ...   
  !>    pelel:  1  3  8  9 13 ...                
  !>                |  |  |                       
  !>                |  |  +--------------+    
  !>                |  |                 |        
  !>                |  +--+              |          
  !>                |     |              |            
  !>    lelel:  2  6  1  3  6  7  8  2  4  7  8  9 ...
  !>    @endverbatim
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------
  
  subroutine graphs_eleele_deallocate(pelel,lelel,pelpo,lelpo,memor)

    integer(ip), pointer, intent(inout), optional :: pelel(:)                        !< Linked list of (element-node) pointer
    integer(ip), pointer, intent(inout), optional :: lelel(:)                        !< Linked list of (element-node) elements
    integer(ip), pointer, intent(inout), optional :: pelpo(:)                        !< Linked list of (element-node) pointer
    integer(ip), pointer, intent(inout), optional :: lelpo(:)                        !< Linked list of (element-node) elements
    integer(8),           intent(inout)           :: memor(2)
    
    if( present(PELEL) ) then
       call memory_deallo(memor,'PELEL','graphs_eleele_deallocate',pelel)
    end if
    if( present(LELEL) ) then
       call memory_deallo(memor,'LELEL','graphs_eleele_deallocate',lelel)
    end if
    if( present(PELPO) ) then
       call memory_deallo(memor,'PELPO','graphs_eleele_deallocate',pelpo)
    end if
    if( present(LELPO) ) then
       call memory_deallo(memor,'LELPO','graphs_eleele_deallocate',lelpo)
    end if
    
  end subroutine graphs_eleele_deallocate
  
  subroutine graphs_eleele(&
       nelem,npoin,mnode,mepoi,lnods,lnnod,&
       pelpo,lelpo,nedge,medge,pelel,lelel,&
       leper,lninv,memor)

    implicit none  
    integer(ip),          intent(in)              :: nelem                           !< Number of elements 
    integer(ip),          intent(in)              :: npoin                           !< Number of nodes 
    integer(ip),          intent(in)              :: mnode                           !< Max. number of nodes per element
    integer(ip),          intent(inout)           :: mepoi                           !< Max. number of element per node
    integer(ip),          intent(in)              :: lnods(mnode,*)                  !< Connectivity array
    integer(ip),          intent(in)              :: lnnod(*)                        !< Array of number of element nodes
    integer(ip), pointer, intent(inout), optional :: pelpo(:)                        !< Linked list of (element-node) pointer
    integer(ip), pointer, intent(inout), optional :: lelpo(:)                        !< Linked list of (element-node) elements
    integer(ip),          intent(out)             :: nedge                           !< Number of edges
    integer(ip),          intent(out)             :: medge                           !< Max. number of edges per element
    integer(ip), pointer, intent(inout)           :: pelel(:)                        !< Linked list of (element-element) pointer
    integer(ip), pointer, intent(inout)           :: lelel(:)                        !< Linked list of (element-element) elements
    integer(ip), pointer, intent(in),    optional :: leper(:)                        !< Element permutation array
    integer(ip), pointer, intent(in),    optional :: lninv(:)                        !< Nodal permutation array
    integer(8),           intent(inout), optional :: memor(2)
    integer(ip)                                   :: ielem,kelem,inode,jelem
    integer(ip)                                   :: melel,jnode,kpoin,pelem
    integer(ip)                                   :: ipoin
    integer(8)                                    :: msize
    type(hash_t)                                  :: ht
    integer(ip), pointer                          :: lista(:)
    integer(ip), pointer                          :: pelp2(:)
    integer(ip), pointer                          :: lelp2(:)
    logical(lg)                                   :: lperm,lelpo_local
    integer(8)                                    :: memor_loc(2)

    memor_loc = optional_argument((/0_8,0_8/),memor)
    !
    ! Initialization
    !
    nullify(lista)
    !
    ! Pelpo and lelpo available?
    !
    lelpo_local = .false.
    if( present(pelpo) .and. present(lelpo) ) then
       if( .not. associated(pelpo) ) then
#ifdef __PGI
          call graphs_elepoi(npoin,nelem,mnode,lnods,lnnod,mepoi,pelpo,lelpo)
#else
          call graphs_elepoi(npoin,nelem,mnode,lnods,lnnod,mepoi,pelpo,lelpo,memor=memor)
#endif
       end if
       pelp2 => pelpo
       lelp2 => lelpo
    else
       lelpo_local = .true.
#ifdef __PGI
       call graphs_elepoi(npoin,nelem,mnode,lnods,lnnod,mepoi,pelp2,lelp2)
#else
       call graphs_elepoi(npoin,nelem,mnode,lnods,lnnod,mepoi,pelp2,lelp2,memor=memor)
#endif
    end if
    !
    ! Permutation?
    !
    lperm = present(leper) .and. present(lninv)
    !
    ! MEPOI is ok?
    !
    if( mepoi <= 0 ) then
       mepoi = -1
       do ipoin = 1,npoin
          mepoi = max(mepoi,pelp2(ipoin+1)-pelp2(ipoin))
       end do
    end if
    !
    ! Allocate memory for PELEL
    !
    call memory_alloca(memor_loc,'PELEL','eleele',pelel,nelem+1_ip)
    melel = mepoi*mnode ! Optimize: it is overdimensionalized

    if( memory_config%high_memory ) then

       !-----------------------------------------------------------------
       !
       ! First option: we have a lot of memory
       !
       !-----------------------------------------------------------------

       msize = int(melel,8)*int(nelem,8)
       !
       ! Compute Hash table (initialize, reset, add and destroy)
       !
       call memory_alloca(memor_loc,'LISTA','eleele',lista,msize)
       call htaini( ht, melel )
       pelel(1) = 1

       if( lperm ) then

          do ielem= 1,nelem
             kelem = leper(ielem)
             call htares( ht, lista(pelel(kelem):) )
             do inode = 1,lnnod(kelem)
                jnode = lnods(inode,kelem)
                kpoin = lninv(jnode)
                do jelem = pelp2(kpoin), pelp2(kpoin+1)-1
                   kelem = lelp2(jelem)
                   if( kelem /= ielem ) then
                      call htaadd( ht, kelem )
                   end if
                end do
             end do
             pelel(ielem+1) = pelel(ielem) + ht % nelem
          end do

       else

          do ielem= 1, nelem
             call htares( ht, lista(pelel(ielem):) )
             do inode = 1, lnnod(ielem)
                jnode = lnods(inode,ielem)
                do jelem = pelp2(jnode), pelp2(jnode+1)-1
                   kelem = lelp2(jelem)
                   if( kelem /= ielem ) then
                      call htaadd( ht, kelem )
                   end if
                end do
             end do
             pelel(ielem+1) = pelel(ielem) + ht % nelem
          end do

       end if

       call htades( ht )
       nedge = pelel(nelem+1)-1
       !
       ! Allocate memory and compute list of adjacancies LELEL
       ! 
       call memory_alloca(memor_loc,'LELEL','eleele',lelel,nedge)
       do ielem = 1, nedge
          lelel(ielem) = lista(ielem)
       end do

    else

       !-----------------------------------------------------------------
       !
       ! Second option: we DO NOT have a lot of memory
       !
       !-----------------------------------------------------------------

       msize = int(melel,8) 
       call memory_alloca(memor_loc,'LISTA','eleele',lista,msize)

       call htaini( ht, melel )
       pelel(1) = 1

       if( lperm ) then

          do ielem = 1,nelem
             pelem = leper(ielem)
             call htares( ht, lista )
             do inode = 1,lnnod(pelem)
                jnode = lnods(inode,pelem)
                kpoin = lninv(jnode)
                do jelem = pelp2(kpoin), pelp2(kpoin+1)-1
                   kelem = lelp2(jelem)
                   if( kelem /= ielem ) then
                      call htaadd( ht, kelem )
                   end if
                end do
             end do
             pelel(ielem+1) = pelel(ielem) + ht%nelem
          end do
          nedge = pelel(nelem+1)-1    

       else

          do ielem = 1,nelem
             call htares( ht, lista )
             do inode = 1,lnnod(ielem)
                jnode = lnods(inode,ielem)
                do jelem = pelp2(jnode), pelp2(jnode+1)-1
                   kelem = lelp2(jelem)
                   if (kelem/=ielem) then
                      call htaadd( ht, kelem )
                   end if
                end do
             end do
             pelel(ielem+1) = pelel(ielem) + ht%nelem
          end do
          nedge = pelel(nelem+1)-1    

       end if

       call memory_alloca(memor_loc,'LELEL','eleele',lelel,nedge)

       if( nedge > 0 ) then

          if( lperm ) then

             do ielem = 1,nelem
                pelem = leper(ielem)
                call htares( ht, lelel(pelel(ielem):) )
                do inode = 1,lnnod(pelem)
                   jnode = lnods(inode,pelem)
                   kpoin = lninv(jnode)
                   do jelem = pelp2(kpoin), pelp2(kpoin+1)-1
                      kelem = lelp2(jelem)
                      if( kelem /= ielem ) then
                         call htaadd( ht, kelem )
                      end if
                   end do
                end do
             end do

          else

             do ielem = 1,nelem
                call htares( ht, lelel(pelel(ielem):) )
                do inode= 1, lnnod(ielem)
                   jnode = lnods(inode,ielem)
                   do jelem = pelp2(jnode), pelp2(jnode+1)-1
                      kelem = lelp2(jelem)
                      if (kelem/=ielem) then
                         call htaadd( ht, kelem )
                      end if
                   end do
                end do
             end do

          end if

          call htades( ht )

       end if
    end if
    !
    ! Deallocate LISTA
    !
    call memory_deallo(memor_loc,'LISTA','eleele',lista)
    !
    ! Maximum number of edges in the mesh
    !
    medge = 0
    do ielem = 1,nelem
       if( pelel(ielem+1)-pelel(ielem) > medge ) then
          medge = pelel(ielem+1)-pelel(ielem)
       end if
    end do
    !
    ! deallocate memory if necessary
    !
    if( lelpo_local ) then
       call memory_deallo(memor_loc,'LELP2','graphs_poipoi',lelp2)
       call memory_deallo(memor_loc,'PELP2','graphs_poipoi',pelp2)
    end if

    if( present(memor) ) memor = memor_loc

  end subroutine graphs_eleele

  !------------------------------------------------------------------------
  !
  !> @brief   Deallocate arrays
  !> @details Deallocate arrays
  !> @author  Guillaume Houzeaux
  !
  !------------------------------------------------------------------------

  subroutine graphs_deallocate_ip14(xarray,ARRAY_NAME,memor)
    implicit none  
    integer(4),  pointer, intent(inout)          :: xarray(:)
    character(len=*),     intent(in),   optional :: ARRAY_NAME
    integer(8),           intent(inout)          :: memor(2)
    if( present(ARRAY_NAME) ) then
       call memory_deallo(memor,trim(ARRAY_NAME),'graphs_deallocate',xarray)
    else
       call memory_deallo(memor,'XARRAY','graphs_deallocate',xarray)
    end if
  end subroutine graphs_deallocate_ip14
  subroutine graphs_deallocate_ip18(xarray,ARRAY_NAME,memor)
    implicit none  
    integer(8),  pointer, intent(inout) :: xarray(:)
    character(len=*),     intent(in),   optional :: ARRAY_NAME
    integer(8),           intent(inout)          :: memor(2)
    if( present(ARRAY_NAME) ) then
       call memory_deallo(memor,trim(ARRAY_NAME),'graphs_deallocate',xarray)
    else
       call memory_deallo(memor,'XARRAY','graphs_deallocate',xarray)
    end if
  end subroutine graphs_deallocate_ip18
  subroutine graphs_deallocate_rp1(xarray,ARRAY_NAME,memor)
    implicit none  
    real(rp),    pointer, intent(inout) :: xarray(:)
    character(len=*),     intent(in),   optional :: ARRAY_NAME
    integer(8),           intent(inout)          :: memor(2)
    if( present(ARRAY_NAME) ) then
       call memory_deallo(memor,trim(ARRAY_NAME),'graphs_deallocate',xarray)
    else
       call memory_deallo(memor,'XARRAY','graphs_deallocate',xarray)
    end if
  end subroutine graphs_deallocate_rp1

  !------------------------------------------------------------------------
  !
  !> @brief   Deallocate some variables needed to compute graph
  !> @details Deallocate some variables needed to compute graph
  !> @author  Guillaume Houzeaux
  !
  !------------------------------------------------------------------------

  subroutine graphs_dealep(ia,ja,IA_NAME,JA_NAME,memor)

    implicit none  
    integer(ip),      pointer, intent(inout)        :: ia(:)
    integer(ip),      pointer, intent(inout)        :: ja(:)
    character(len=*),          intent(in), optional :: IA_NAME
    character(len=*),          intent(in), optional :: JA_NAME
    integer(8),           intent(inout)          :: memor(2)

    if( present(IA_NAME) .and. present(JA_NAME) ) then
       call memory_deallo(memor,trim(IA_NAME),'graphs_dealep',ia)
       call memory_deallo(memor,trim(jA_NAME),'graphs_dealep',ja)
    else
       call memory_deallo(memor,'IA','graphs_dealep',ia)
       call memory_deallo(memor,'JA','graphs_dealep',ja)
    end if
    
  end subroutine graphs_dealep

  !------------------------------------------------------------------------
  !
  !> @brief   Deallocate a graph
  !> @details Deallocate a graph
  !> @author  Guillaume Houzeaux
  !
  !------------------------------------------------------------------------

  subroutine graphs_deagra(ia,ja,memor)

    implicit none  
    integer(ip), pointer, intent(inout) :: ia(:)
    integer(ip), pointer, intent(inout) :: ja(:)
    integer(8),           intent(inout) :: memor(2)

    call memory_deallo(memor,'IA','graphs_deagra',ia)
    call memory_deallo(memor,'JA','graphs_deagra',ja)

  end subroutine graphs_deagra

  !------------------------------------------------------------------------
  !
  !> @brief   Compute a subgraph from a graph using a condition
  !> @details This subroutine create a graph IA2, JA2 
  !> @author  Guillaume Houzeaux
  !
  !------------------------------------------------------------------------

  subroutine graphs_congra(nn,ia,ja,ifnode,ia2,ja2,memor)

    implicit none  
    integer(ip),          intent(in)    :: nn
    integer(ip),          intent(in)    :: ia(*)
    integer(ip),          intent(in)    :: ja(*)
    logical(lg), pointer, intent(in)    :: ifnode(:)
    integer(ip), pointer, intent(inout)   :: ia2(:)
    integer(ip), pointer, intent(inout)   :: ja2(:)
    integer(8),           intent(inout) :: memor(2)
    integer(ip)                         :: ii,jj,kk,jjz
    integer(ip)                         :: nnz2
    !
    ! Initialization
    !
    nullify(ia2)
    nullify(ja2)
    !
    ! Copy graph in temporary graphs IACPY, JACPY
    !
    call memory_alloca(memor,'IA2','graphs_congra',ia2,nn+1_ip)
    nnz2   = 0
    ia2(1) = 1
    do ii = 1,nn
       kk = 0
       if( ifnode(ii) ) then
          do jjz = ia(ii),ia(ii+1)-1
             jj = ja(jjz)
             if( ifnode(jj) ) kk = kk + 1
          end do
       end if
       ia2(ii+1) = ia2(ii) + kk 
       nnz2 = nnz2 + kk
    end do
    call memory_alloca(memor,'JA2','graphs_congra',ja2,nnz2)
    nnz2 = 0
    do ii = 1,nn
       if( ifnode(ii) ) then
          do jjz = ia(ii),ia(ii+1)-1
             jj = ja(jjz)
             if( ifnode(jj) ) then
                nnz2 = nnz2 + 1
                ja2(nnz2) = jj
             end if
          end do
       end if
    end do

  end subroutine graphs_congra


  subroutine graphs_deaper(permr,invpr,memor)
    !------------------------------------------------------------------------
    !****f* domain/graphs_deaper
    ! NAME
    !    graphs_deaper
    ! DESCRIPTION
    !    Deallocate permutation arrays
    ! USED BY
    !    par_partit
    !***
    !------------------------------------------------------------------------
    implicit none  
    integer(ip), pointer, intent(inout) :: permr(:)
    integer(ip), pointer, intent(inout) :: invpr(:)
    integer(8),           intent(inout)          :: memor(2)

    call memory_deallo(memor,'INVPR','graphs_permut',invpr)   
    call memory_deallo(memor,'PERMR','graphs_permut',permr)  

  end subroutine graphs_deaper

  subroutine graphs_renmat(ndof,nn,invpr,ia,ja,ianew,janew,an,memor)
    !------------------------------------------------------------------------
    !****f* domain/graphs_rrenmat
    ! NAME
    !    graphs_rrenmat
    ! DESCRIPTION
    !    This subroutine permutes a matrix
    ! USED BY
    !***
    !------------------------------------------------------------------------
    implicit none  
    integer(ip),          intent(in)    :: ndof
    integer(ip),          intent(in)    :: nn
    integer(ip), pointer, intent(in)    :: invpr(:)
    integer(ip), pointer, intent(in)    :: ia(:)
    integer(ip), pointer, intent(in)    :: ja(:)
    integer(ip), pointer, intent(in)    :: ianew(:)
    integer(ip), pointer, intent(in)    :: janew(:)
    real(rp),    pointer, intent(inout) :: an(:,:,:)
    integer(8),           intent(inout) :: memor(2)
    real(rp),    pointer                :: ancpy(:,:,:)
    integer(ip)                         :: ii,jj,kk,iiz
    integer(ip)                         :: idof,jdof,nnz
    integer(ip)                         :: iiznew,iinew,jjnew

    nullify(ancpy)
    nnz = size(an,3_ip,KIND=ip)
    if( nnz /= ia(nn+1) - 1 ) call runend('GRAPHS_RENMAT: WRONG MATRIX DIMENSIONS: '&
         //integer_to_string(nnz)//', '//integer_to_string(ia(nn+1)))
    nullify(ancpy)

    call memory_alloca(memor,'ANCPY','graphs_renmat',ancpy,ndof,ndof,nnz)

    do iiz = 1,nnz
       do idof = 1,ndof
          do jdof = 1,ndof
             ancpy(jdof,idof,iiz) = an(jdof,idof,iiz)
          end do
       end do
    end do

    do iinew = 1,nn
       ii = invpr(iinew)
       do iiznew = ianew(iinew),ianew(iinew+1)-1
          jjnew = janew(iiznew)
          jj    = invpr(jjnew)
          iiz   = ia(ii)
          kk    = ja(iiz)
          do while( jj /= kk )
             iiz = iiz + 1
             kk  = ja(iiz)
          end do
          if( iiz > ia(ii+1)-1 ) call runend('OUPS')
          do idof = 1,ndof
             do jdof = 1,ndof
                an(jdof,idof,iiznew) = ancpy(jdof,idof,iiz)
             end do
          end do
       end do
    end do

    call memory_deallo(memor,'ANCPY','graphs_renmat',ancpy)

  end subroutine graphs_renmat

  subroutine graphs_copyij_deallocate(iacpy,jacpy,memor_opt,IACPY_NAME,JACPY_NAME)
    !------------------------------------------------------------------------
    !****f* domain/graphs_copyij
    ! NAME
    !    graphs_cpoyij
    ! DESCRIPTION
    !    This subroutine copy a graph
    ! INPUT
    !    IA, JA ......... Graph
    ! OUTPUT
    !    IACPY = IA, JACPY = JA
    ! USED BY
    !***
    !------------------------------------------------------------------------

    integer(ip), pointer,  intent(inout)           :: iacpy(:)  
    integer(ip), pointer,  intent(inout)           :: jacpy(:)
    integer(8),            intent(inout)           :: memor_opt(2)
    character(len=*),      intent(in),    optional :: IACPY_NAME
    character(len=*),      intent(in),    optional :: JACPY_NAME
    !
    ! Copy graph in IACPY, JACPY
    !

    if( present(IACPY_NAME) ) then
       call memory_deallo(memor_opt,trim(IACPY_NAME),'graphs_copyij',iacpy)
    else
       call memory_deallo(memor_opt,'IACPY','graphs_copyij',iacpy)
    end if
    if( present(JACPY_NAME) ) then
       call memory_deallo(memor_opt,trim(JACPY_NAME),'graphs_copyij',jacpy)
    else
       call memory_deallo(memor_opt,'JACPY','graphs_copyij',jacpy)
    end if
    
  end subroutine graphs_copyij_deallocate

  subroutine graphs_copyij(nn,ia,ja,iacpy,jacpy,memor_opt,IACPY_NAME,JACPY_NAME)
    !------------------------------------------------------------------------
    !****f* domain/graphs_copyij
    ! NAME
    !    graphs_cpoyij
    ! DESCRIPTION
    !    This subroutine copy a graph
    ! INPUT
    !    IA, JA ......... Graph
    ! OUTPUT
    !    IACPY = IA, JACPY = JA
    ! USED BY
    !***
    !------------------------------------------------------------------------
    implicit none  
    integer(ip),           intent(in)              :: nn
    integer(ip), pointer,  intent(in)              :: ia(:)
    integer(ip), pointer,  intent(in)              :: ja(:)
    integer(ip), pointer,  intent(inout)           :: iacpy(:)  
    integer(ip), pointer,  intent(inout)           :: jacpy(:)
    integer(8),            intent(inout)           :: memor_opt(2)
    character(len=*),      intent(in),    optional :: IACPY_NAME
    character(len=*),      intent(in),    optional :: JACPY_NAME
    integer(ip)                                    :: ii,iiz,nnz
    !
    ! Copy graph in IACPY, JACPY
    !
    if( nn > 0 ) then
       nnz = memory_size(ja)
       if( nnz > 0 ) then
          if( nnz /= ia(nn+1) - 1 ) call runend('GRAPHS_COPYIJ: WRONG MATRIX DIMENSIONS'&
               //integer_to_string(nnz)//', '//integer_to_string(ia(nn+1)-1))

          if( present(IACPY_NAME) ) then
             call memory_alloca(memor_opt,trim(IACPY_NAME),'graphs_copyij',iacpy,nn+1_ip)
          else
             call memory_alloca(memor_opt,'IACPY','graphs_copyij',iacpy,nn+1_ip)
          end if
          if( present(JACPY_NAME) ) then
             call memory_alloca(memor_opt,trim(JACPY_NAME),'graphs_copyij',jacpy,nnz)
          else
             call memory_alloca(memor_opt,'JACPY','graphs_copyij',jacpy,nnz)
          end if
          
          do ii = 1,nn+1
             iacpy(ii) = ia(ii)
          end do
          do iiz = 1,nnz
             jacpy(iiz) = ja(iiz)
          end do
       end if
    end if

  end subroutine graphs_copyij

  !-----------------------------------------------------------------------
  !> @file    graphs_gtband.f90
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Compute the bandwth and profile of the mesh
  !> @details Compute the bandwth and profile of the mesh
  !-----------------------------------------------------------------------

  subroutine graphs_gtband(npoin,ja,ia,bandw,profi)
    implicit none
    integer(ip),          intent(in)  :: npoin
    integer(ip), pointer, intent(in)  :: ia(:)
    integer(ip), pointer, intent(in)  :: ja(:)
    integer(ip),          intent(out) :: bandw
    real(rp),             intent(out) :: profi
    integer(ip)                       :: ipoin,izdomin,band,izdom,bandloc,ipmax,jpmax  

    bandw =  0_ip
    profi =  0.0_rp
    !naved =  0      ! Average number of edges
    !nmied =  1e6    ! Min number of edges
    !nmaed = -1e6    ! Max number of edges

    do ipoin = 1,npoin
       !
       ! Initialize local bandwth
       ! 
       bandloc = 0_ip
       !
       ! Loop on neighbors
       !
       do izdom = ia(ipoin),ia(ipoin+1)-1
          izdomin = ja(izdom)
          if( ipoin /= izdomin ) then
             band = abs(izdomin-ipoin)
             !
             ! Test bandwth
             !
             if( band > bandw ) then
                bandw = band
                ipmax = ipoin
                jpmax = izdomin
             endif
             !
             ! Test profile
             !
             if( izdomin < ipoin ) then
                if( band > bandloc ) bandloc = band
             end if
          end if
       end do
       !
       ! Accumulate profile
       !
       profi = profi + real(bandloc,rp)

    end do

  end subroutine graphs_gtband

  !-----------------------------------------------------------------------
  !> @file    graphs_nzecof.f90
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   computes the number of non-zero coefficients graph
  !> @details computes the number of non-zero coefficients of a
  !>           mesh graph stored in compressed sparse row (CSR) format 
  !-----------------------------------------------------------------------

  subroutine graphs_nzecof(&
       mnode,lnods,lnnod,ltype,nlist,ncoef,nz,liste,&
       touch,me,only_edges,npoin_opt)
    implicit none
    integer(ip), intent(in)             :: mnode
    integer(ip), intent(in)             :: lnods(mnode,*)
    integer(ip), intent(in)             :: lnnod(*)
    integer(ip), intent(in)             :: ltype(*)
    integer(ip), intent(in)             :: nlist
    integer(ip), intent(in)             :: ncoef
    integer(ip), intent(in)             :: liste(nlist)
    integer(ip), intent(inout)          :: nz
    logical(lg), intent(inout)          :: touch(ncoef)
    integer(ip), intent(in)             :: me
    logical(lg), intent(in)             :: only_edges
    integer(ip), intent(in),   optional :: npoin_opt
    integer(ip)                         :: jelem,jnode,jpoin,nnodj,jposi,jlist
    integer(ip)                         :: kelem,knode,kpoin,nnodk,kposi,klist
!    integer(ip)                         :: jedge,jtype,ipoin,jnod1,jnod2

    if( only_edges ) then
       !
       ! Use only edges lower part to construct graph
       ! Consider only edges by imposing TOUCH = .FALSE.
       !
       call runend('GRAPHS_MNZECOF: NOT LONGER SUPPORTED')
!!$       ipoin = -me
!!$       do jlist = 1,nlist                              
!!$          jelem = liste(jlist)
!!$          nnodj = lnnod(jelem)
!!$          do jnode = 1,nnodj
!!$             jposi        = (jlist-1)*mnode+jnode
!!$             touch(jposi) = .true.
!!$          end do
!!$       end do
!!$       do jlist = 1,nlist   
!!$          jelem = liste(jlist)
!!$          jtype = abs(ltype(jelem))
!!$          
!!$          do jedge = 1,element_type(jtype) % number_edges
!!$             jnod1 = element_type(jtype) % list_edges(1,jedge)
!!$             jnod2 = element_type(jtype) % list_edges(2,jedge)
!!$             if( lnods(jnod1,jelem) == ipoin ) then
!!$                jposi        = (jlist-1)*mnode+jnod2
!!$                touch(jposi) = .false.                
!!$             else if( lnods(jnod2,jelem) == ipoin ) then
!!$                jposi        = (jlist-1)*mnode+jnod1
!!$                touch(jposi) = .false.
!!$             end if
!!$          end do
!!$          
!!$       end do
    end if

    if( me == 0 ) then
       !
       ! All graph
       !
       do jlist = 1,nlist                                      ! Loop over those elements 
          jelem = liste(jlist)                                 ! where the point is
          nnodj = lnnod(jelem)
          do jnode = 1,nnodj
             jpoin = lnods(jnode,jelem)
             jposi = (jlist-1)*mnode+jnode
             if( .not. touch(jposi) ) then                     ! Position not touched           
                do klist = 1,nlist                             ! Search other elements 
                   kelem = liste(klist)                        ! where JPOIN is and 
                   nnodk = lnnod(kelem)
                   do knode = 1,nnodk                          ! touch their position
                      kpoin = lnods(knode,kelem)
                      if( kpoin == jpoin ) then
                         kposi = (klist-1)*mnode+knode
                         touch(kposi) = .true.
                      end if
                   end do
                end do
                if( present(npoin_opt) ) then
                   if( jpoin <= npoin_opt ) nz = nz + 1
                else
                   nz = nz + 1
                end if
             end if
          end do
       end do

    else if( me < 0 ) then
       !
       ! Lower graph
       !
       do jlist = 1,nlist                                      ! Loop over those elements 
          jelem = liste(jlist)                                 ! where the point is
          nnodj = lnnod(jelem)
          do jnode = 1,nnodj
             jpoin = lnods(jnode,jelem)
             if( jpoin < -me ) then
                jposi = (jlist-1)*mnode+jnode
                if( .not. touch(jposi) ) then                     ! Position not touched           
                   do klist = 1,nlist                             ! Search other elements 
                      kelem = liste(klist)                        ! where JPOIN is and 
                      nnodk = lnnod(kelem)
                      do knode = 1,nnodk                          ! touch their position
                         kpoin = lnods(knode,kelem)
                         if( kpoin == jpoin ) then
                            kposi = (klist-1)*mnode+knode
                            touch(kposi) = .true.
                         end if
                      end do
                   end do
                   if( present(npoin_opt) ) then
                      if( jpoin <= npoin_opt ) nz = nz + 1
                   else                      
                      nz = nz + 1
                   end if
                end if
             end if
          end do
       end do

    else
       !
       ! All graph without diagonal
       !
       do jlist = 1,nlist                                      ! Loop over those elements 
          jelem = liste(jlist)                                 ! where the point is
          nnodj = lnnod(jelem)
          do jnode = 1,nnodj
             jpoin = lnods(jnode,jelem)
             if( jpoin /= me ) then
                jposi = (jlist-1)*mnode+jnode
                if( .not. touch(jposi) ) then                     ! Position not touched           
                   do klist = 1,nlist                             ! Search other elements 
                      kelem = liste(klist)                        ! where JPOIN is and 
                      nnodk = lnnod(kelem)
                      do knode = 1,nnodk                          ! touch their position
                         kpoin = lnods(knode,kelem)
                         if( kpoin == jpoin ) then
                            kposi = (klist-1)*mnode+knode
                            touch(kposi) = .true.
                         end if
                      end do
                   end do
                   if( present(npoin_opt) ) then
                      if( jpoin <= npoin_opt ) nz = nz + 1
                   else                      
                      nz = nz + 1
                   end if
                end if
             end if
          end do
       end do

    end if

  end subroutine graphs_nzecof

  !-----------------------------------------------------------------------
  !> @file    graphs_arrind.f90
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Constructs the arrays of indexes for a mesh graph
  !> @details Constructs the arrays of indexes for a mesh graph
  !>          These are organized as follows (CSR format):
  !>          IA(II) = coefficient of the graph where IPOIN starts,
  !>          JA(NZ) = column of the NZ coefficient of the graph
  !>          mesh graph stored in compressed sparse row (CSR) format 
  !-----------------------------------------------------------------------

  subroutine graphs_arrind(&
       mnode,lnods,lnnod,ltype,nlist,ncoef,liste,touch,&
       nz,ipoin,ia,ja,me,only_edges,npoin_opt)
    implicit none
    integer(ip),          intent(in)    :: mnode
    integer(ip),          intent(in)    :: lnods(mnode,*)
    integer(ip),          intent(in)    :: lnnod(*)
    integer(ip),          intent(in)    :: ltype(*)
    integer(ip),          intent(in)    :: nlist
    integer(ip),          intent(in)    :: ncoef
    integer(ip),          intent(in)    :: ipoin
    integer(ip),          intent(in)    :: liste(nlist)
    integer(ip),          intent(inout) :: nz
    logical(lg),          intent(inout) :: touch(ncoef)
    integer(ip), pointer, intent(inout) :: ia(:)
    integer(ip), pointer, intent(inout) :: ja(:)
    integer(ip),          intent(in)    :: me
    logical(lg),          intent(in)    :: only_edges
    integer(ip), optional,intent(in)    :: npoin_opt
    integer(ip)                         :: jelem,jnode,jpoin,nnodj,jposi,jlist
    integer(ip)                         :: kelem,knode,kpoin,nnodk,kposi,klist
    integer(ip)                         :: jtype,jedge,jnod1,jnod2

    ia(ipoin) = nz

    if( only_edges ) then
       !
       ! Use only edges lower part to construct graph
       !
       do jlist = 1,nlist                              
          jelem = liste(jlist)
          nnodj = lnnod(jelem)
          do jnode = 1,nnodj
             jposi = (jlist-1)*mnode+jnode
             touch(jposi) = .true.
          end do
       end do
       do jlist = 1,nlist   
          jelem = liste(jlist)          
          jtype = abs(ltype(jelem))

          do jedge = 1,element_type(jtype) % number_edges
             jnod1 = element_type(jtype) % list_edges(1,jedge)
             jnod2 = element_type(jtype) % list_edges(2,jedge)
             if( lnods(jnod1,jelem) == ipoin ) then
                jposi        = (jlist-1)*mnode+jnod2
                touch(jposi) = .false.                
             else if( lnods(jnod2,jelem) == ipoin ) then
                jposi        = (jlist-1)*mnode+jnod1
                touch(jposi) = .false.
             end if
          end do

          
       end do
    end if

    if( me == 0 ) then

       do jlist = 1,nlist
          jelem = liste(jlist)
          nnodj = lnnod(jelem)
          do jnode = 1,nnodj
             jpoin = lnods(jnode,jelem)
             jposi = (jlist-1) * mnode + jnode
             if( .not. touch(jposi) ) then
                do klist = 1,nlist
                   kelem = liste(klist)
                   nnodk = lnnod(kelem)
                   do knode = 1,nnodk
                      kpoin = lnods(knode,kelem)
                      if(kpoin==jpoin) then
                         kposi = (klist-1)*mnode+knode
                         touch(kposi) = .true.
                      end if
                   end do
                end do
                if( present(npoin_opt) ) then
                   if( jpoin <= npoin_opt ) then
                      ja(nz) = jpoin
                      nz = nz+1
                   end if
                else
                   ja(nz) = jpoin
                   nz = nz+1
                end if
             end if
          end do
       end do

    else if( me < 0 ) then

       do jlist = 1,nlist
          jelem = liste(jlist)
          nnodj = lnnod(jelem)
          do jnode = 1,nnodj
             jpoin = lnods(jnode,jelem)
             if( jpoin < -me ) then
                jposi = (jlist-1) * mnode + jnode
                if( .not. touch(jposi) ) then
                   do klist = 1,nlist
                      kelem = liste(klist)
                      nnodk = lnnod(kelem)
                      do knode = 1,nnodk
                         kpoin = lnods(knode,kelem)
                         if(kpoin==jpoin) then
                            kposi = (klist-1)*mnode+knode
                            touch(kposi) = .true.
                         end if
                      end do
                   end do
                   if( present(npoin_opt) ) then
                      if( jpoin <= npoin_opt ) then
                         ja(nz) = jpoin
                         nz = nz+1
                      end if
                   else
                      ja(nz) = jpoin
                      nz = nz+1
                   end if
                end if
             end if
          end do
       end do

    else

       do jlist = 1,nlist
          jelem = liste(jlist)
          nnodj = lnnod(jelem)
          do jnode = 1,nnodj
             jpoin = lnods(jnode,jelem)
             if( jpoin /= me ) then
                jposi = (jlist-1) * mnode + jnode
                if( .not. touch(jposi) ) then
                   do klist = 1,nlist
                      kelem = liste(klist)
                      nnodk = lnnod(kelem)
                      do knode = 1,nnodk
                         kpoin = lnods(knode,kelem)
                         if(kpoin==jpoin) then
                            kposi = (klist-1)*mnode+knode
                            touch(kposi) = .true.
                         end if
                      end do
                   end do
                   if( present(npoin_opt) ) then
                      if( jpoin <= npoin_opt ) then
                         ja(nz) = jpoin
                         nz = nz+1
                      end if
                   else
                      ja(nz) = jpoin
                      nz = nz+1
                   end if
                end if
             end if
          end do
       end do

    end if

  end subroutine graphs_arrind

  !-----------------------------------------------------------------------
  !> @file    graphs_mergli.f90
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Merges to list of nodes
  !> @details Merges to list of nodes
  !-----------------------------------------------------------------------

  subroutine graphs_mergli(lista,lsize,pnode,lnods,ltype,me,only_edges,npoin_opt)
    implicit none
    integer(ip), intent(inout)          :: lsize
    integer(ip), intent(inout)          :: lista(*)
    integer(ip), intent(in)             :: pnode
    integer(ip), intent(in)             :: lnods(pnode)
    integer(ip), intent(in)             :: ltype
    integer(ip), intent(in)             :: me
    logical(lg), intent(in)             :: only_edges
    integer(ip), intent(in),   optional :: npoin_opt
    integer(ip)                         :: ii, jj, n1, n2
    logical(lg)                         :: noEncontrado
    integer(ip)                         :: lnods_2(pnode)
    integer(ip)                         :: itype,i1,i2,iedge,meabs
    integer(ip)                         :: ipoin1,ipoin2

    if( only_edges ) then
       !
       ! Use only edges lower part to construct graph
       !
       meabs = abs(me)
       do ii = 1,pnode
          lnods_2(ii) = 0
       end do
       itype = abs(ltype)       
       
       do iedge = 1,element_type(itype) % number_edges
          i1     = element_type(itype) % list_edges(1,iedge) 
          i2     = element_type(itype) % list_edges(2,iedge) 
          ipoin1 = lnods(i1)
          ipoin2 = lnods(i2)
          if( ipoin1 == meabs .or. ipoin2 == meabs ) then
             lnods_2(i1) = lnods(i1)
             lnods_2(i2) = lnods(i2)
          end if
       end do
       
       do ii = 1, pnode
          n1 = lnods_2(ii)
          if( n1 < -me .and. n1 /= 0 ) then
             jj = 1
             noEncontrado = .true.
             do while( jj <= lsize .and. noEncontrado )
                n2 = lista(jj)
                if( n1 == n2 ) then
                   noEncontrado = .false.
                end if
                jj = jj + 1
             end do
             if( noEncontrado ) then
                if( present(npoin_opt) ) then
                   if( n1 <= npoin_opt ) then
                      lsize = lsize + 1
                      lista(lsize) = n1                      
                   end if
                else
                   lsize = lsize + 1
                   lista(lsize) = n1
                end if
             end if
          end if
       end do

    else 

       if( me == 0 ) then

          do ii = 1,pnode
             n1 = lnods(ii)
             jj = 1
             noEncontrado = .true.
             do while( jj <= lsize .and. noEncontrado )
                n2 = lista(jj)
                if ( n1 == n2 ) then
                   noEncontrado = .false.
                end if
                jj = jj + 1
             end do
             if ( noEncontrado ) then
                if( present(npoin_opt) ) then
                   if( n1 <= npoin_opt ) then
                      lsize = lsize + 1
                      lista(lsize) = n1                      
                   end if
                else
                   lsize = lsize + 1
                   lista(lsize) = n1
                end if
             end if
          end do

       else if( me < 0 ) then

          do ii = 1, pnode
             n1 = lnods(ii)
             if( n1 < -me ) then
                jj = 1
                noEncontrado = .true.
                do while( jj <= lsize .and. noEncontrado )
                   n2 = lista(jj)
                   if( n1 == n2 ) then
                      noEncontrado = .false.
                   end if
                   jj = jj + 1
                end do
                if( noEncontrado ) then
                   if( present(npoin_opt) ) then
                      if( n1 <= npoin_opt ) then
                         lsize = lsize + 1
                         lista(lsize) = n1                      
                      end if
                   else
                      lsize = lsize + 1
                      lista(lsize) = n1
                   end if
                end if
             end if
          end do

       else

          do ii = 1, pnode
             n1 = lnods(ii)
             if( n1 /= me ) then
                jj = 1
                noEncontrado = .true.
                do while( jj <= lsize .and. noEncontrado )
                   n2 = lista(jj)
                   if( n1 == n2 ) then
                      noEncontrado = .false.
                   end if
                   jj = jj + 1
                end do
                if( noEncontrado ) then
                   if( present(npoin_opt) ) then
                      if( n1 <= npoin_opt ) then
                         lsize = lsize + 1
                         lista(lsize) = n1                      
                      end if
                   else
                      lsize = lsize + 1
                      lista(lsize) = n1
                   end if
                end if
             end if
          end do

       end if

    end if

  end subroutine graphs_mergli

  !-----------------------------------------------------------------------
  !
  !> @brief   Compute a graph
  !> @details Compute the element-element graph using face connectivity
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine graphs_element_element_graph(&
       meshe,whow,whalo,wdiag,ia,ja,memor)
    
    type(mesh_type),                    intent(inout)         :: meshe          !< Mesh type
    character(*),                       intent(in)            :: whow
    character(*),                       intent(in)            :: whalo
    character(*),    optional,          intent(in)            :: wdiag
    integer(ip),     optional, pointer, intent(inout)         :: ia(:)
    integer(ip),     optional, pointer, intent(inout)         :: ja(:)
    integer(8),                         intent(inout)         :: memor(2)
    integer(ip)                                               :: mepoi_2
    integer(ip)                                               :: medge,nelem_2
    integer(ip),     pointer                                  :: lelpo_2(:)
    integer(ip),     pointer                                  :: pelpo_2(:)
    logical(lg)                                               :: use_r_dom_c_dom

    nullify(lelpo_2)
    nullify(pelpo_2)
    if( trim(whalo) == 'INCLUDING HALO' ) then        
       nelem_2 = meshe % nelem_2
    else
       nelem_2 = meshe % nelem       
    end if
    use_r_dom_c_dom = .false.

    if( trim(whalo) == 'INCLUDING HALO' ) then       

#ifdef __PGI
       call graphs_elepoi(&
            meshe % npoin_2,meshe % nelem_2,meshe % mnode,meshe % lnods,&
            meshe % lnnod,mepoi_2,pelpo_2,lelpo_2) 
#else
       call graphs_elepoi(&
            meshe % npoin_2,meshe % nelem_2,meshe % mnode,meshe % lnods,&
            meshe % lnnod,mepoi_2,pelpo_2,lelpo_2,memor=memor) 
#endif

    else

       if( associated(meshe % r_dom) .and. associated(meshe % c_dom) ) then
          use_r_dom_c_dom = .true.
          pelpo_2 => meshe % r_dom
          lelpo_2 => meshe % c_dom
       else
#ifdef __PGI
          call graphs_elepoi(&
               meshe % npoin,meshe % nelem,meshe % mnode,meshe % lnods,&
               meshe % lnnod,mepoi_2,pelpo_2,lelpo_2) 
#else
          call graphs_elepoi(&
               meshe % npoin,meshe % nelem,meshe % mnode,meshe % lnods,&
               meshe % lnnod,mepoi_2,pelpo_2,lelpo_2,memor=memor) 
#endif
       end if

    end if

    if(      trim(whow) == 'SHARING NODES' ) then

       if( present(ia) .and. present(ja) ) then
#ifdef __PGI
          call graphs_eleele(&
               meshe % nelem_2,meshe % npoin_2,meshe % mnode,mepoi_2,meshe % lnods,meshe % lnnod,&
               pelpo_2,lelpo_2,meshe % nzelm_2,medge,ia,ja)
#else
          call graphs_eleele(&
               meshe % nelem_2,meshe % npoin_2,meshe % mnode,mepoi_2,meshe % lnods,meshe % lnnod,&
               pelpo_2,lelpo_2,meshe % nzelm_2,medge,ia,ja,memor=memor)
#endif
       else
#ifdef __PGI
          call graphs_eleele(&
               meshe % nelem_2,meshe % npoin_2,meshe % mnode,mepoi_2,meshe % lnods,meshe % lnnod,&
               pelpo_2,lelpo_2,meshe % nzelm_2,medge,meshe % r_elm_2,meshe % c_elm_2)
#else
          call graphs_eleele(&
               meshe % nelem_2,meshe % npoin_2,meshe % mnode,mepoi_2,meshe % lnods,meshe % lnnod,&
               pelpo_2,lelpo_2,meshe % nzelm_2,medge,meshe % r_elm_2,meshe % c_elm_2,memor=memor)
#endif
       end if
       
    else if( trim(whow) == 'SHARING FACES' ) then

       if( present(ia) .and. present(ja) ) then
          call graphs_eleele_faces_1(&
               nelem_2        ,nelem_2      ,meshe % mnode,meshe % lnods  ,&
               meshe % ltype  ,pelpo_2,lelpo_2,&
               meshe % nzelm_2,medge        ,ia,&
               ja,wdiag,memor=memor)
       else
          call graphs_eleele_faces_1(&
               nelem_2        ,nelem_2      ,meshe % mnode,meshe % lnods  ,&
               meshe % ltype  ,pelpo_2,lelpo_2,&
               meshe % nzelm_2,medge        ,meshe % r_elm_2,&
               meshe % c_elm_2,wdiag,memor=memor)
       end if
       
    else
       call runend('GRAPHS_ELEMENT_ELEMENT_GRAPH: DO NOT KNOW WHAT TO DO 1')
    end if

    if( .not. use_r_dom_c_dom ) then
       call memory_deallo(memor,'LELPO_2','graphs_element_element_graph',lelpo_2)
       call memory_deallo(memor,'PELPO_2','graphs_element_element_graph',pelpo_2)
    end if

  end subroutine graphs_element_element_graph

  !-----------------------------------------------------------------------
  !
  !> @brief   Compute a graph
  !> @details Compute the element-element graph using face connectivity
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine graphs_eleele_faces_0(&    
       meshe,pelpo,lelpo,pelel,lelel,wdiag,memor)
    
    type(mesh_type_basic),                   intent(inout) :: meshe          !< Mesh type
    integer(ip),          optional, pointer, intent(in)    :: pelpo(:)       !< Linked list of (element-node) pointer
    integer(ip),          optional, pointer, intent(in)    :: lelpo(:)       !< Linked list of (element-node) elements
    integer(ip),                    pointer, intent(inout) :: pelel(:)       !< Linked list of (element-node) pointer
    integer(ip),                    pointer, intent(inout) :: lelel(:)       !< Linked list of (element-node) elements
    character(*),         optional,          intent(in)    :: wdiag
    integer(8),                              intent(inout) :: memor(2)
    integer(ip)                                            :: nedge,medge
    integer(ip)                                            :: mepoi
    integer(ip),                    pointer                :: pelpo_loc(:)  
    integer(ip),                    pointer                :: lelpo_loc(:)  

    nullify(lelpo_loc,pelpo_loc)
    
    if( present(lelpo) .and. present(pelpo) ) then
       lelpo_loc => lelpo
       pelpo_loc => pelpo       
    else       
       call graphs_elepoi_0(meshe,mepoi,pelpo_loc,lelpo_loc,memor=memor)
    end if
    
    call graphs_eleele_faces_1(& 
         meshe % nelem,meshe % nelem,meshe % mnode,meshe % lnods,meshe % ltype,&
         pelpo_loc,lelpo_loc,nedge,medge,pelel,lelel,wdiag,memor=memor)
    
    if( present(lelpo) .and. present(pelpo) ) then
       continue
    else
       call graphs_elepoi_deallocate(pelpo_loc,lelpo_loc,memor=memor)       
    end if
    
  end subroutine graphs_eleele_faces_0
   
  subroutine graphs_eleele_faces_1(& 
       nelem,nelem_2,mnode,lnods,ltype,pelpo,&
       lelpo,nedge,medge,pelel,lelel,wdiag,memor)
  
    integer(ip),          intent(in)           :: nelem                           !< Number of elements
    integer(ip),          intent(in)           :: nelem_2                         !< Number of elements including possibly halo
    integer(ip),          intent(in)           :: mnode                           !< Max. number of nodes per element
    integer(ip), pointer, intent(in)           :: lnods(:,:)                      !< Connectivity array
    integer(ip), pointer, intent(in)           :: ltype(:)                        !< Array of element types
    integer(ip), pointer, intent(in)           :: pelpo(:)                        !< Linked list of (element-node) pointer
    integer(ip), pointer, intent(in)           :: lelpo(:)                        !< Linked list of (element-node) elements
    integer(ip),          intent(out)          :: nedge                           !< Number of edges
    integer(ip),          intent(out)          :: medge                           !< Max. number of edges per element
    integer(ip), pointer, intent(inout)        :: pelel(:)                        !< Linked list of (element-element) pointer
    integer(ip), pointer, intent(inout)        :: lelel(:)                        !< Linked list of (element-element) elements
    character(*),         intent(in), optional :: wdiag
    integer(8),           intent(inout)        :: memor(2)
    integer(ip)                                :: ielem,inode
    integer(ip)                                :: iface,ielty
    integer(ip)                                :: dummi,inodb
    integer(ip)                                :: mface,pnodb
    integer(ip)                                :: pnodi,mnodb
    integer(ip), pointer                       :: faces(:,:,:)
    logical(lg)                                :: include_diag

    if( nelem <= 0 ) return
    !
    ! Check if diagonal should be included
    !
    include_diag = .false.
    if( present(wdiag) ) then
       if( trim(wdiag) == 'INCLUDING DIAGONAL' ) include_diag = .true.
    end if
    !
    ! Maximum number of faces and nodes per face
    !
    mface = 0
    mnodb = 0
    do ielem = 1,nelem_2                                         
       ielty = abs(ltype(ielem))
       mface = max(mface,element_type(ielty) % number_faces)
       do iface = 1,element_type(ielty) % number_faces
          mnodb = max(mnodb,element_type(ielty) % node_faces(iface))
       end do
    end do
    !
    ! Allocate memory
    !
    nullify(faces)       
    call memory_alloca(memor,'FACES','graphs_eleele_faces',faces,mnodb+1_ip,mface,nelem_2)
    call memory_alloca(memor,'PELEL','graphs_eleele_faces',pelel,nelem+1_ip )
    !
    ! Construct and sort FACES 
    !
    do ielem = 1,nelem_2                                        
       ielty = abs(ltype(ielem))
       do iface = 1,element_type(ielty) % number_faces
          pnodb = element_type(ielty) % node_faces(iface)
          do inodb = 1,pnodb
             inode = element_type(ielty) % list_faces(inodb,iface)
             faces(inodb,iface,ielem) = lnods(inode,ielem)
          end do
          call sortin(pnodb,faces(1,iface,ielem))
       end do
    end do
    !
    ! Compute FACES
    !
    call graphs_face_type(nelem,ltype,pelpo,lelpo,nedge,faces)
    !
    ! Allocate memory for adjacancies LELEL
    !
    if( include_diag ) nedge = nedge + nelem
    call memory_alloca(memor,'LELEL','graphs_eleele_faces',lelel,nedge)
    !
    ! Compute PELEL and LELEL
    !
    dummi    = 0
    pelel(1) = 1
    !
    !*OMP  PARALLEL DO SCHEDULE (GUIDED)      & 
    !*OMP  DEFAULT (NONE)                     &
    !*OMP  PRIVATE (dummi,ielem,ielty,iface)  &
    !*OMP  SHARED  (pelel,ltype,faces,lelel,nelem) 
    !
    if( include_diag ) then
       do ielem = 1,nelem
          dummi          = dummi + 1
          pelel(ielem+1) = pelel(ielem+1) + 1
          lelel(dummi)   = ielem
          ielty = abs(ltype(ielem))
          do iface = 1,element_type(ielty) % number_faces
             if( faces(1,iface,ielem) == 0 ) then
                pnodi          = element_type(ielty) % node_faces(iface)
                dummi          = dummi + 1
                pelel(ielem+1) = pelel(ielem+1) + 1
                !lelel(dummi)   = faces(pnodi+1,iface,ielem)
                lelel(dummi)   = faces(mnodb+1,iface,ielem)
             end if
          end do
          pelel(ielem+1) = pelel(ielem) + pelel(ielem+1)
       end do
    else
       do ielem = 1,nelem
          ielty = abs(ltype(ielem))
          do iface = 1,element_type(ielty) % number_faces
             if( faces(1,iface,ielem) == 0 ) then
                pnodi          = element_type(ielty) % node_faces(iface)
                dummi          = dummi + 1
                pelel(ielem+1) = pelel(ielem+1) + 1
                !lelel(dummi)   = faces(pnodi+1,iface,ielem)
                lelel(dummi)   = faces(mnodb+1,iface,ielem)
             end if
          end do
          pelel(ielem+1) = pelel(ielem) + pelel(ielem+1)
       end do
    end if
    !
    ! Deallocate memory of FACES
    !
    call memory_deallo(memor,'FACES','par_elmgra',faces)
    !
    !
    ! Order graph and compute maximum number of edges in the mesh
    !
    medge = 0
    do ielem = 1,nelem
       dummi = pelel(ielem+1) - pelel(ielem) 
       medge = max(medge,pelel(ielem+1)-pelel(ielem))
       if( dummi > 0 ) call maths_heap_sort(2_ip,dummi,lelel(pelel(ielem):))
    end do

  end subroutine graphs_eleele_faces_1

  !-----------------------------------------------------------------------
  !
  !> @brief   Compute a graph
  !> @details Compute the list of faces
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine graphs_list_faces(&
       nelem,mnode,mnodb,nelty,mface,lnods,ltype,&
       pelpo,lelpo,nfacg,lface,lelfa,&
       lchek,memor)
    
    integer(ip),          intent(in)           :: nelem                           !< Number of elements
    integer(ip),          intent(in)           :: mnode                           !< Max. number of nodes per element
    integer(ip),          intent(in)           :: mnodb                           !< Max. number of nodes per face
    integer(ip),          intent(in)           :: nelty                           !< Number of element types
    integer(ip),          intent(in)           :: mface                           !< Number of faces for all elements
    integer(ip),          intent(in)           :: lnods(mnode,*)                  !< Connectivity array
    integer(ip),          intent(in)           :: ltype(*)                        !< Array of element types
    integer(ip), pointer, intent(in)           :: pelpo(:)                        !< Linked list of (element-node) pointer
    integer(ip), pointer, intent(in)           :: lelpo(:)                        !< Linked list of (element-node) elements
    integer(ip),          intent(out)          :: nfacg                           !< Number of faces
    integer(ip), pointer, intent(inout)          :: lface(:,:)                      !< List of faces
    type(i1p),   pointer, intent(inout)          :: lelfa(:)                        !< List of alement faces
    logical(lg), pointer, intent(in), optional :: lchek(:)                        !< Linked list of (element-node) elements
    integer(8),           intent(inout), optional        :: memor(2)
    integer(ip)                                :: ielem,inode,jelem
    integer(ip)                                :: iface,jface,ielty,ipoin
    integer(ip)                                :: ilist,ielpo,inodb
    integer(ip)                                :: pepoi,jelty,pnodf
    integer(ip), pointer                       :: faces(:,:,:)
    logical(lg)                                :: equal_faces,ichek,jchek
    integer(8)                                    :: memor_loc(2)

    memor_loc = optional_argument((/0_8,0_8/),memor)
    !
    ! Allocate memory 
    !
    nullify(faces)
    call memory_alloca(memor_loc,'FACES','graphs_list_faces',faces,mnodb+1_ip,mface,nelem)
    call memory_alloca(memor_loc,'LELFA','graphs_list_faces',lelfa,nelem)
    !
    ! Construct and sort FACES 
    !
    do ielem = 1,nelem  
       ichek = .true.
       if( present(lchek) ) then
          ichek = lchek(ielem)
       end if
       if( ichek ) then
          ielty = abs(ltype(ielem))       
          call memory_alloca(memor_loc,'LELFA % L','lgface',lelfa(ielem)%l,element_type(ielty) % number_faces)
          do iface = 1,element_type(ielty) % number_faces
             do inodb = 1,element_type(ielty) % node_faces(iface) 
                inode = element_type(ielty) % list_faces(inodb,iface)
                faces(inodb,iface,ielem) = lnods(inode,ielem)
             end do
             faces(mnodb+1,iface,ielem) = 1
             call sortin(element_type(ielty) % node_faces(iface),faces(1,iface,ielem))
          end do
       end if
    end do
    !
    ! Compute FACES
    !
    nfacg = 0
    do ielem = 1,nelem                                                          ! Compare the faces and 
       ichek = .true.
       if( present(lchek) ) then
          ichek = lchek(ielem)
       end if
       if( ichek ) then
          ielty = abs(ltype(ielem))                                                ! eliminate the repited faces
          do iface = 1,element_type(ielty) % number_faces
             if( faces(mnodb+1,iface,ielem) > 0 ) then
                nfacg = nfacg + 1
                ipoin = faces(1,iface,ielem)
                if( ipoin /= 0 ) then
                   ilist = 1
                   pepoi = pelpo(ipoin+1)-pelpo(ipoin)
                   ielpo = pelpo(ipoin)-1
                   do while( ilist <= pepoi )
                      ielpo = ielpo+1
                      jelem = lelpo(ielpo)
                      if( jelem /= ielem .and. jelem <= nelem ) then
                         jchek = .true.
                         if( present(lchek) ) then
                            jchek = lchek(jelem)
                         end if
                         if( jchek ) then
                            jelty = abs(ltype(jelem))                                 ! eliminate the repited faces
                            jface = 0
                            do while( jface < element_type(jelty) % number_faces )
                               jface = jface + 1
                               if( faces(1,jface,jelem) /= 0 ) then
                                  equal_faces = .true.
                                  inodb       = 0
                                  do while( equal_faces .and. inodb /= element_type(jelty) % node_faces(jface) )  
                                     inodb = inodb + 1
                                     if( faces(inodb,iface,ielem) /= faces(inodb,jface,jelem) ) &
                                          equal_faces = .false.
                                  end do
                                  if( equal_faces ) then
                                     faces(mnodb+1,iface,ielem) =  jelem                              ! Keep IELEM face
                                     faces(mnodb+1,jface,jelem) = -ielem                              ! Elminate JELEM face
                                     faces(      1,iface,ielem) = -jface                              ! Remember IFACE face
                                     jface                      =  element_type(jelty) % number_faces ! Exit JFACE do
                                     ilist                      =  pepoi                              ! Exit JELEM do  
                                  end if
                               end if
                            end do
                         end if
                      end if
                      ilist = ilist + 1
                   end do
                end if
             end if
          end do
       end if
    end do
    !
    ! List of faces if required
    !
    call memory_alloca(memor_loc,'LFACE','graphs_list_faces',lface,4_ip,nfacg)     
    nfacg = 0
    do ielem = 1,nelem                                            ! Compare the faces and 
       ichek = .true.
       if( present(lchek) ) then
          ichek = lchek(ielem)
       end if
       if( ichek ) then
          ielty = abs(ltype(ielem))                                  ! eliminate the repeated faces
          do iface = 1,element_type(ielty) % number_faces
             if( faces(mnodb+1,iface,ielem) > 0 ) then
                nfacg = nfacg + 1
                pnodf = element_type(ielty) % node_faces(iface) 
                if( faces(1,iface,ielem) < 0 ) then
                   jelem                   =  faces(mnodb+1,iface,ielem)
                   jface                   = -faces(      1,iface,ielem)
                   lelfa(ielem) % l(iface) =  nfacg
                   lelfa(jelem) % l(jface) =  nfacg
                   lface(1,nfacg)          =  ielem
                   lface(2,nfacg)          =  jelem
                   lface(3,nfacg)          =  iface
                   lface(4,nfacg)          =  jface
                else
                   lelfa(ielem) % l(iface) =  nfacg
                   lface(1,nfacg)          =  ielem
                   lface(2,nfacg)          =  0
                   lface(3,nfacg)          =  iface
                   lface(4,nfacg)          =  0
                end if
             end if
          end do
       end if
    end do
    !
    ! Deallocate memory of FACES
    !
    call memory_deallo(memor_loc,'FACES','graphs_list_faces',faces)

    if( present(memor) ) memor = memor_loc

  end subroutine graphs_list_faces

  !------------------------------------------------------------------------
  !
  !> @brief   Color a graph
  !> @details Color a graph
  !> @author  Guillaume Houzeaux
  !
  !------------------------------------------------------------------------
  
  subroutine graphs_coloring_deallocate(lcolo,ia_color,ja_color,memor)
    
    integer(ip), intent(inout), pointer, optional :: lcolo(:)
    integer(ip), intent(inout), pointer, optional :: ia_color(:)
    integer(ip), intent(inout), pointer, optional :: ja_color(:)
    integer(8),  intent(inout)                    :: memor(2)
    
    if( present(lcolo) ) then
       call memory_deallo(memor,'LCOLO','graphs_coloring_deallocate',lcolo)
    end if
    if( present(ia_color) ) then
       call memory_deallo(memor,'IA_COLOR','graphs_coloring_deallocate',ia_color)
    end if
    if( present(ja_color) ) then
       call memory_deallo(memor,'jA_COLOR','graphs_coloring_deallocate',ja_color)
    end if
    
  end subroutine graphs_coloring_deallocate


  subroutine graphs_coloring(nvert,ia,ja,lcolo,ncolo,ia_color,ja_color,COLOR_NAME,IA_NAME,JA_NAME,memor)
    integer(ip), intent(in)                       :: nvert
    integer(ip), intent(in),    pointer           :: ia(:)
    integer(ip), intent(in),    pointer           :: ja(:)
    integer(ip), intent(inout), pointer           :: lcolo(:)
    integer(ip), intent(out),            optional :: ncolo
    integer(ip), intent(inout), pointer, optional :: ia_color(:)
    integer(ip), intent(inout), pointer, optional :: ja_color(:)
    character(len=*), intent(in),        optional :: COLOR_NAME
    character(len=*), intent(in),        optional :: IA_NAME
    character(len=*), intent(in),        optional :: JA_NAME
    integer(8),       intent(inout)               :: memor(2)
    integer(ip)                                   :: ivert,jvert,kvert,mvert
    integer(ip)                                   :: icolo,jcolo,izver,jzver
    integer(ip)                                   :: istack,ihuge,nvert0
    integer(ip)                                   :: nstack,ivert_last,color
    integer(ip),              pointer             :: lstack(:)

    nullify(lstack)
    if( .not. associated(lcolo) ) then
       if( present(COLOR_NAME) ) then
          call memory_alloca(memor,trim(COLOR_NAME),'graphs_coloring',lcolo,nvert)
       else
          call memory_alloca(memor,'LCOLO','graphs_coloring',lcolo,nvert)
       end if
    end if
    call memory_alloca(memor,'LSTACK','graphs_coloring',lstack,nvert)

    nvert0 = 0
    color  = 0
    !kvert=0

    do while( nvert0 /= nvert ) 
       !
       ! Goto next color
       !
       color      = color + 1
       istack     = 0
       nstack     = 0
       ivert_last = 0

       do 
          !
          ! Next seed point
          !
          if( istack == nstack ) then
             ivert = ivert_last
             do while( ivert < nvert )
                ivert = ivert + 1
                if( lcolo(ivert) == 0 ) then
                   ihuge         =  huge(1_ip)
                   lcolo(ivert)  =  color
                   lstack(1)     =  ivert     
                   nstack        =  1
                   istack        =  0
                   ivert_last    =  ivert
                   ivert         =  nvert  + 1
                end if
             end do
             if( ivert /= nvert + 1 ) exit
          end if
          !
          ! Next point in stack and choose color
          !
          istack = istack + 1   
          ivert  = lstack(istack)    
          if( lcolo(ivert) == color ) then
             icolo = ihuge
          else
             icolo = color
          end if
          !
          ! Color neighbors JVERT of IVERT with color ICOLO
          !
          do izver = ia(ivert),ia(ivert+1)-1
             jvert = ja(izver)
             if( lcolo(jvert) == 0 ) then
                jcolo = icolo
                if( icolo == color ) then
                   jzver = ia(jvert)-1
                   do while( jzver < ia(jvert+1)-1 )
                      jzver = jzver + 1
                      mvert = ja(jzver)
                      if( lcolo(mvert) == color ) then
                         jcolo = ihuge
                         jzver = ia(jvert+1)
                      end if
                   end do
                end if
                !      kvert          = kvert + 1
                lcolo(jvert)   = jcolo
                nstack         = nstack + 1
                lstack(nstack) = jvert                                
             end if
          end do

       end do
       !
       ! Paint IHUGE null color
       ! 
       do ivert = 1,nvert
          if( lcolo(ivert) == ihuge ) then
             lcolo(ivert) = 0
          else if( lcolo(ivert) == color ) then
             nvert0 = nvert0 + 1
          end if
       end do
       !
       ! If no color is prescribe, just pain with color 1
       !
       if( .not. present(ncolo) ) nvert0 = nvert

    end do
    !
    ! Number of colors
    !
    if( present(ncolo) ) ncolo = color
    call memory_deallo(memor,'LSTACK','graphs_coloring',lstack)
    !
    ! Color linked list
    !
    if( present(ncolo) .and. present(ia_color) .and. present(ja_color) ) then
       if( present(IA_NAME) ) then
          call memory_alloca(memor,trim(IA_NAME),'graphs_coloring',ia_color,ncolo+1)
       else
          call memory_alloca(memor,'IA_COLOR','graphs_coloring',ia_color,ncolo+1)
       end if
       if( present(JA_NAME) ) then
          call memory_alloca(memor,trim(JA_NAME),'graphs_coloring',ja_color,nvert)
       else
          call memory_alloca(memor,'JA_COLOR','graphs_coloring',ja_color,nvert)
       end if
       ia_color(1) = 1
       do icolo = 1,ncolo
          jvert = 0
          do ivert = 1,nvert
             if( lcolo(ivert) == icolo ) then
                kvert           = jvert + ia_color(icolo) 
                jvert           = jvert + 1
                ja_color(kvert) = ivert
             end if
          end do
          ia_color(icolo+1) = ia_color(icolo) + jvert
       end do
    end if

  end subroutine graphs_coloring

  !------------------------------------------------------------------------
  !
  !> @brief   Color a graph Greedy Descendent
  !> @details Color a graph
  !> @author  Antoni Artigues
  !
  !------------------------------------------------------------------------
  subroutine graphs_coloring_greedy(nvert,ia,ja,lcolo,ncolo,ia_color,ja_color,COLOR_NAME,IA_NAME,JA_NAME,memor)
    integer(ip), intent(in)                     :: nvert
    integer(ip), intent(in),  pointer           :: ia(:)
    integer(ip), intent(in),  pointer           :: ja(:)
    integer(ip), intent(inout), pointer           :: lcolo(:)
    integer(ip), intent(out),          optional :: ncolo
    integer(ip), intent(inout), pointer, optional :: ia_color(:)
    integer(ip), intent(inout), pointer, optional :: ja_color(:)
    character(len=*), intent(in),        optional :: COLOR_NAME
    character(len=*), intent(in),        optional :: IA_NAME
    character(len=*), intent(in),        optional :: JA_NAME
    integer(8),       intent(inout)             :: memor(2) 
    integer(ip)                                 :: ivert,jvert,kvert,mvert
    integer(ip)                                 :: icolo,izver
    integer(ip)                                 :: color
    !Nuevos
    integer(ip),              pointer           :: adjacentColors(:)
    integer(ip),              pointer           :: clrOrder(:)
    integer(ip),              pointer           :: vertexOrder(:)
    integer(ip),              pointer           :: vertexRank(:)
    integer(ip),              pointer           :: clrCnt(:)
    integer(ip)                                 :: clr, clrindex, aux
    integer(ip)                                 :: ii
    integer(ip)                                 :: maxAdjacent

    !--------------------------------
    !--INITIALIZATIONS 1--
    !--------------------------------
    nullify(clrOrder)
    nullify(vertexOrder)
    nullify(clrCnt)
    nullify(adjacentColors)
    nullify(vertexRank)
    if( .not. associated(lcolo) ) then
       if( present(COLOR_NAME) ) then
          call memory_alloca(memor,trim(COLOR_NAME),'graphs_coloring',lcolo,nvert)
       else
          call memory_alloca(memor,'LCOLO','graphs_coloring',lcolo,nvert)
       end if
    end if
    call memory_alloca(memor,'VERTEXORDER','graphs_coloring',vertexOrder,nvert)
    call memory_alloca(memor,'VERTEXRANK','graphs_coloring',vertexRank,nvert)


    !--------------------------------
    !--ORDER VERTEX IN ASCENDENT ORDER
    !--------------------------------
    do ivert = 1,nvert
       aux = (ia(ivert+1)-1) - ia(ivert)
       vertexRank(ivert) = aux
       vertexOrder(ivert) = ivert
    end do

    call orderVertex(vertexRank,vertexOrder)

    !--------------------------------
    !--INITIALIZATIONS 2--
    !--------------------------------
    !--Set maximum adjacency
    maxAdjacent = vertexRank(1)  + 1    
    call memory_deallo(memor,'LSTACK','graphs_coloring',vertexRank)
    call memory_alloca(memor,'CLRORDER','graphs_coloring',clrOrder,nvert)
    call memory_alloca(memor,'CLRCNT','graphs_coloring',clrCnt,nvert)
    call memory_alloca(memor,'ADJACENTCOLORS','graphs_coloring',adjacentColors,maxAdjacent)

    do ivert = 1, nvert        
       lcolo(ivert) = 0
       clrCnt(ivert) = 0
       clrOrder(ivert) = ivert
    end do

    !Write(*,*)'Graph Coloring Greedy: Vertex coloring'
    !--------------------------------    
    !--VERTEX COLORING IN ORDER
    !--------------------------------    
    clrindex = 1
    do ivert = 1,nvert                
       mvert = vertexOrder(ivert)
       !--get the adyacent colors of the actual vertex
       do jvert = 1, maxAdjacent
          adjacentColors(jvert) = 0
       end do
       ii = 1
       do izver = ia(mvert),ia(mvert+1)-1
          jvert = ja(izver)
          adjacentColors(ii) = lcolo(jvert) 
          ii = ii + 1            
       end do

       !--obtain vertex color
       call obtainVertexColor(adjacentColors, clrOrder, nvert,  maxAdjacent, clr, clrindex)

       !--Assign color
       lcolo(mvert) = clr
       clrCnt(clr) = clrCnt(clr) + 1
       !--reorder color list order
       call reorderColorOrder(clrCnt, clrOrder, clr, clrindex, nvert)       
    end do

    !Obtain the number of colors
    color = 0
    do ivert = 1, nvert
       if (clrCnt(ivert) > 0) then
          color = color + 1
       end if
    end do

    call memory_deallo(memor,'clrOrder','graphs_coloring',clrOrder)
    call memory_deallo(memor,'clrCnt','graphs_coloring',clrCnt)
    call memory_deallo(memor,'vertexOrder','graphs_coloring',vertexOrder)
    call memory_deallo(memor,'adjacentColors','graphs_coloring',adjacentColors)


    !Write(*,*)'Graph Coloring Greedy: Linked list coloring'
    !--------------------------------    
    !---COLOREADO DE LA LINKED LIST
    !--------------------------------    
    !
    ! Number of colors
    !
    if( present(ncolo) ) ncolo = color
    !call memory_deallo(memor_dom,'LSTACK','graphs_coloring',lstack)
    !
    ! Color linked list
    !
    if( present(ncolo) .and. present(ia_color) .and. present(ja_color) ) then
       if( present(IA_NAME) ) then
          call memory_alloca(memor,trim(IA_NAME),'graphs_coloring',ia_color,ncolo+1)
       else
          call memory_alloca(memor,'IA_COLOR','graphs_coloring',ia_color,ncolo+1)
       end if
       if( present(JA_NAME) ) then
          call memory_alloca(memor,trim(JA_NAME),'graphs_coloring',ja_color,nvert)
       else
          call memory_alloca(memor,'JA_COLOR','graphs_coloring',ja_color,nvert)
       end if
       ia_color(1) = 1
       do icolo = 1,ncolo
          jvert = 0
          do ivert = 1,nvert
             if( lcolo(ivert) == icolo ) then
                kvert           = jvert + ia_color(icolo) 
                jvert           = jvert + 1 
                ja_color(kvert) = ivert
             end if
          end do
          ia_color(icolo+1) = ia_color(icolo) + jvert
       end do
    end if

  end subroutine graphs_coloring_greedy


  !------------------------------------------------------------------------
  !
  !> @brief   QuickSort ordering the vertex in descendent order
  !> @details Auxiliary funcion of graphs_coloring_greedy
  !> @author  Antoni Artigues
  !
  !------------------------------------------------------------------------
  recursive subroutine orderVertex(A,B)
    integer(ip), intent(in out), dimension(:) :: A,B
    integer(ip) :: iq

    if(size(A,KIND=ip) > 1) then
       call Partition(A,B, iq)
       call orderVertex(A(:iq-1),B(:iq-1))
       call orderVertex(A(iq:),B(:iq-1))
    endif
  end subroutine orderVertex

  !------------------------------------------------------------------------
  !
  !> @brief   Partition of the QuickSort ordering algorithm
  !> @details Auxiliary funcion of orderVertex from graphs_coloring_greedy
  !> @author  Antoni Artigues
  !
  !------------------------------------------------------------------------
  subroutine Partition(A,B, marker)
    integer(ip), intent(in out), dimension(:) :: A,B
    integer(ip), intent(out) :: marker
    integer(ip) :: i, j
    integer(ip) :: temp
    integer(ip) :: x      !pivot point
    x = A(1)
    i= 0
    j= size(A,KIND=ip) + 1

    do
       j = j-1
       do
          if (A(j) >= x) exit
          j = j-1
       end do
       i = i+1
       do
          if (A(i) <= x) exit
          i = i+1
       end do
       if (i < j) then
          ! exchange A(i) and A(j)
          temp = A(i)
          A(i) = A(j)
          A(j) = temp

          temp = B(i)
          !print*,'par1'
          B(i) = B(j)
          !print*,'par2'
          B(j) = temp
       elseif (i == j) then
          marker = i+1
          return
       else
          marker = i
          return
       endif
    end do
  end subroutine Partition

  !------------------------------------------------------------------------
  !
  !> @brief   Obtain next color for a given adjacentColors array
  !> @details Auxiliary funcion of graphs_coloring_greedy
  !> @author  Antoni Artigues
  !
  !------------------------------------------------------------------------
  subroutine obtainVertexColor(adjacentColors, clrOrder, nvert, maxAdjacent, clr, idx )    
    integer(ip), intent(in),  pointer           :: adjacentColors(:)
    integer(ip), intent(in),  pointer           :: clrOrder(:)
    integer(ip), intent(in)                      :: nvert
    integer(ip), intent(in)                      :: maxAdjacent
    integer(ip), intent(out)                      :: clr
    integer(ip), intent(out)                      :: idx
    integer(ip)                                   :: found, i, j


    main_loop:do i = 1,nvert
       clr = clrOrder(i)

       !search the color in the adjacent colors list
       found = 0
       do j = 1, maxAdjacent
          if (adjacentColors(j) == clr) then
             found = 1
          end if
       end do

       if (found == 0) then
          idx = i
          exit main_loop
       end if
    end do main_loop

  end subroutine obtainVertexColor

  !------------------------------------------------------------------------
  !
  !> @brief   For a given idx and color reorders de list in ascendent order
  !> @details Auxiliary funcion of graphs_coloring_greedy
  !> @author  Antoni Artigues
  !
  !------------------------------------------------------------------------
  subroutine reorderColorOrder(clrCnt, clrOrder, clr, idx, nvert)
    integer(ip), intent(in),  pointer           :: clrCnt(:)
    integer(ip), intent(inout),  pointer        :: clrOrder(:)
    integer(ip), intent(in)                     :: clr
    integer(ip), intent(in)                     :: idx
    integer(ip), intent(in)                     :: nvert
    integer(ip)                                 :: countVer2, countVer3, ordAux, idxAux, i


    countVer2 = clrCnt(clrOrder(idx))
    idxAux = -1
    !Lo movemos hacia atrás
    main_loop:do i = idx -1, 1, -1
       countVer3 = clrCnt(clrOrder(i))
       !swap the two values if is greater than its predecesor
       if (countVer3 > countVer2) then
          ordAux = clrOrder(i+1)
          clrOrder(i+1) = clrOrder(i)
          clrOrder(i) = ordAux
       else
          exit main_loop             
       end if
    end do main_loop

    !Lo movemos hacia adelante
    main_loop2:do i = idx + 1, nvert
       countVer3 = clrCnt(clrOrder(i))
       !swap the two values if is less than its following
       if (countVer3 < countVer2 .and. countVer3 > 0) then
          ordAux = clrOrder(i-1)
          clrOrder(i-1) = clrOrder(i)
          clrOrder(i) = ordAux
       else
          exit main_loop2             
       end if
    end do main_loop2

  end subroutine reorderColorOrder

  !------------------------------------------------------------------------
  !
  !> @brief   Deallocate list of faces
  !> @details Deallocate list of faces
  !> @author  Guillaume Houzeaux
  !
  !------------------------------------------------------------------------

  subroutine graphs_deallocate_list_faces(lface,lelfa,lchek,memor)
    integer(ip), pointer, intent(inout)        :: lface(:,:)                      !< List of faces
    type(i1p),   pointer, intent(inout)        :: lelfa(:)                        !< List of alement faces
    logical(lg), pointer, intent(in), optional :: lchek(:)                        !< Linked list of (element-node) elements
    integer(8),           intent(inout), optional :: memor(2)
    integer(ip)                                   :: ielem
    logical(lg)                                   :: ichek
    integer(8)                                    :: memor_loc(2)

    memor_loc = optional_argument((/0_8,0_8/),memor)

    if( present(lchek) ) then
       do ielem = 1,size(lelfa,KIND=ip)
          ichek = lchek(ielem) 
          if( ichek ) then
             call memory_deallo(memor_loc,'LELFA % L','lgface',lelfa(ielem) % l)
          end if
       end do
    else
       do ielem = 1,size(lelfa,KIND=ip)
          call memory_deallo(memor_loc,'LELFA % L','lgface',lelfa(ielem) % l)
       end do
    end if
    call memory_deallo(memor_loc,'LELFA','graphs_list_faces',lelfa)
    call memory_deallo(memor_loc,'LFACE','graphs_list_faces',lface)  
   
    if( present(memor) ) memor = memor_loc

  end subroutine graphs_deallocate_list_faces

  !------------------------------------------------------------------------
  !
  !> @author  Paula Cordoba and Guillaume Houzeaux
  !> @brief   Deallocate list of faces
  !> @details Deallocate list of faces
  !>          \verbatim
  !>
  !>                            |vecto|
  !>          <--------------------------------------------
  !>
  !>          +--------------+-------------+--------------+-------------+
  !>          |              |             |              |             |
  !>          +--------------+-------------+--------------+-------------+
  !>          Imposed inflow  Free inflow   Interior       Outflow
  !>          (boundaries)    (boundaries)                 (Boundaries)
  !>
  !>          \endverbatim
  !
  !------------------------------------------------------------------------

  subroutine graphs_number_along_vector(meshe,vecto,angle_stream,permr,invpr,lgrou,kfl_fixno,ngrou_gs,message,memor)

    type(mesh_type), intent(in)                       :: meshe          !< Mesh type
    real(rp),        intent(in),    pointer           :: vecto(:,:)     !< Direction vector
    real(rp),        intent(in)                       :: angle_stream   !< Stream angle
    integer(ip),     intent(inout), pointer           :: permr(:)       !< Permutation   NEW = permr(OLD)
    integer(ip),     intent(inout), pointer           :: invpr(:)       !< Inverse perm. OLD = invpr = invpr(NEW)
    integer(ip),     intent(inout), pointer           :: lgrou(:)       !< Groups
    integer(ip),     intent(in),    pointer, optional :: kfl_fixno(:,:) !< Fixity type
    integer(ip),     intent(out),            optional :: ngrou_gs       !< Number of groups
    character(*),    intent(in),             optional :: message
    integer(8),      intent(inout)                    :: memor(2)
    integer(ip)                                       :: ipoin
!    integer(ip)                                       :: igrou
    integer(ip)                                       :: izdom,jpoin
    integer(ip)                                       :: kpoin
    integer(ip)                                       :: mpoin,ibopo
    integer(ip)                                       :: npinl
    integer(ip)                                       :: kbopo
    integer(ip)                                       :: nmarkt,nmark
    integer(ip)                                       :: ngrou
!    integer(ip)                                       :: tpoin,jelem
!    integer(ip)                                       :: unit1,unit2
    integer(ip)                                       :: npoin_ngrou
    real(rp)                                          :: cangl,cangm
    real(rp)                                          :: epsil,denom
    real(rp)                                          :: vpoin(3),alpha
    real(rp)                                          :: vrefe(3)
    !
    ! Local pointers
    !
    integer(ip),                  pointer             :: lorde(:) 
    logical(lg),                  pointer             :: lmark(:) 
    logical(lg),                  pointer             :: linle(:)
    real(rp),                     pointer             :: vmodu(:)
    !
    ! Mesh arrays required
    !
    integer(ip)                                       :: npoin
    integer(ip)                                       :: npoi1
    integer(ip)                                       :: pdime
    integer(ip),                  pointer             :: ia(:)
    integer(ip),                  pointer             :: ja(:)
    integer(ip),                  pointer             :: lnoch(:)
    integer(ip),                  pointer             :: lpoty(:)
    real(rp),                     pointer             :: coord(:,:)
    real(rp),                     pointer             :: exnor(:,:,:)

    if( present(message) ) then
       if( trim(message) == 'DEALLOCATE' ) then
          call memory_deallo(memor,'PERMR','graphs_numbering_along_vector',permr)
          call memory_deallo(memor,'INVPR','graphs_numbering_along_vector',invpr)
          call memory_deallo(memor,'LGROU','graphs_numbering_along_vector',lgrou)
          return
       end if
    end if

    !----------------------------------------------------------------------
    !
    ! Point to mesh arrays
    !      
    !----------------------------------------------------------------------

    npoin =  meshe % npoin 
    npoi1 =  meshe % npoi1 
    pdime =  meshe % ndime 
    ia    => meshe % r_dom
    ja    => meshe % c_dom
    lnoch => meshe % lnoch
    lpoty => meshe % lpoty
    coord => meshe % coord
    exnor => meshe % exnor

    !----------------------------------------------------------------------
    !
    ! Initialization
    !      
    !----------------------------------------------------------------------

    nullify(lorde)
    nullify(lmark)
    nullify(linle)
    nullify(vmodu)

    nmarkt = npoin           ! Number of marked points
    npinl  = 0               ! Number of inflow nodes
    epsil  = epsilon(1.0_rp) ! Very snall number
    vpoin  = 0.0_rp          ! Node velocity
    vrefe  = 0.0_rp          ! Reference velocity

    !----------------------------------------------------------------------
    !
    ! Allocate memory 
    !
    !----------------------------------------------------------------------
    !
    ! Allocate permutation arrays
    !
    if( .not. associated(permr) ) call memory_alloca(memor,'PERMR','graphs_numbering_along_vector',permr,npoin)
    if( .not. associated(invpr) ) call memory_alloca(memor,'INVPR','graphs_numbering_along_vector',invpr,npoin)
    if( .not. associated(lgrou) ) call memory_alloca(memor,'LGROU','graphs_numbering_along_vector',lgrou,npoin)
    permr(1:npoin) = 0
    invpr(1:npoin) = 0
    lgrou(1:npoin) = 0
    !do ipoin = 1,npoin
    !   invpr(ipoin) = ipoin
    !   permr(ipoin) = ipoin
    !end do
    !return
    !
    ! Local arrays
    ! 
    call memory_alloca(memor,'LORDE','graphs_numbering_along_vector',lorde,npoin)
    call memory_alloca(memor,'LMARK','graphs_numbering_along_vector',lmark,npoin)
    call memory_alloca(memor,'LINLE','graphs_numbering_along_vector',linle,npoin)
    call memory_alloca(memor,'VMODU','graphs_numbering_along_vector',vmodu,npoin)
    !
    ! Mark hole nodes (they are like imposed Dirichlet nodes) 
    !
    if( associated(lnoch) ) then
       do ipoin = 1,npoin
          if( lnoch(ipoin) == NOHOL ) then
             call runend('GRAPHS_RENUMBER_ALONG_VECTOR: NOT CODED')
             nmarkt       = nmarkt - 1
             lmark(ipoin) = .true.  
             lgrou(ipoin) = -huge(1_ip)        
             !                lgrou(ipoin) = -1
          end if
       end do
    end if
    !
    ! Interface node
    !
    do ipoin = npoi1+1,npoin
       nmarkt       = nmarkt - 1
       lmark(ipoin) = .true.
       lgrou(ipoin) = -huge(1_ip)
    end do
    !
    ! Velocity module, LINLE and INVPR
    !
    do ipoin = 1,npoi1
       vmodu(ipoin) = dot_product(vecto(1:pdime,ipoin),vecto(1:pdime,ipoin))
       invpr(ipoin) = ipoin
    end do
    !
    ! Put outflow node at the end
    !
    do ipoin = 1,npoi1
       ibopo = lpoty(ipoin)
       if( ibopo > 0 ) then
          cangl = dot_product(exnor(1:pdime,1,ibopo),vecto(1:pdime,ipoin))
          if( cangl > 0.0_rp ) vmodu(ipoin) = 0.0_rp
       end if
    end do
    !
    ! Order nodes by increasing module 
    !
    call heapsortri(1_ip,npoi1,vmodu,invpr)
    !
    ! Number first inflow prescribed fixity
    !
    if( present(kfl_fixno) ) then
       if( associated(kfl_fixno) ) then
          do jpoin = 1,npoi1
             ipoin = invpr(jpoin) ! OLD = INVPR(NEW)
             ibopo = lpoty(ipoin)
             if( maxval(kfl_fixno(:,ipoin)) > 0 .and. ibopo > 0 ) then
                cangl = dot_product(exnor(1:pdime,1,ibopo),vecto(1:pdime,ipoin))
                if( cangl < 0.0_rp ) then
                   npinl        = npinl + 1
                   lorde(npinl) = ipoin
                   linle(ipoin) = .true.
                end if
             end if
          end do
       end if
    end if
    !
    ! Then put non-prescribed inflow nodes 
    !
    do jpoin = 1,npoi1
       ipoin = invpr(jpoin) ! OLD = INVPR(NEW)
       if( .not. linle(ipoin) ) then
          ibopo = lpoty(ipoin)
          if( ibopo > 0 ) then
             cangl = dot_product(exnor(1:pdime,1,ibopo),vecto(1:pdime,ipoin))
             if( cangl < 0.0_rp ) then
                npinl        = npinl + 1
                lorde(npinl) = ipoin
                linle(ipoin) = .true.
             end if
          end if
       end if
    end do
    !
    ! Number remaining nodes
    !
    kpoin = npinl
    do jpoin = 1,npoi1
       ipoin = invpr(jpoin)
       if( .not. linle(ipoin) ) then
          kpoin = kpoin + 1
          lorde(kpoin) = ipoin
       end if
    end do
    !
    ! Interface node
    !
    do ipoin = npoi1+1,npoin
       kpoin = kpoin + 1
       invpr(kpoin) = ipoin
    end do

    !----------------------------------------------------------------------
    !
    ! Construct groups
    !
    !----------------------------------------------------------------------

    nmark   = 0
    ngrou   = 0     ! Change to avoid the ngrou = ngrou - 1 done below
    mpoin   = 1

    do while ( nmark < nmarkt )
       !
       ! Look for seed for group IGROU
       !
       ipoin = lorde(mpoin)
       do while( lmark(ipoin) )
          mpoin = mpoin + 1
          ipoin = lorde(mpoin)
       end do
       !
       ! Mark seed
       !
       ngrou        = ngrou + 1    ! Change to avoid the ngrou = ngrou - 1 below
       nmark        = nmark + 1
       lmark(ipoin) = .true.
       lgrou(ipoin) = ngrou  
       invpr(nmark) = ipoin
       npoin_ngrou  = 1

       loop_igrou: do
          cangl = 0.0_rp
          cangm = 0.0_rp
          kpoin = 0
          !
          ! Look for best neighbor KPOIN (cannot be an inflow node)
          !
          loop_izdom: do izdom = ia(ipoin),ia(ipoin+1)-1
             jpoin = ja(izdom)
             !if( ( .not. linle(jpoin) ) .and. ( .not. lmark(jpoin) ) .and. ipoin /= jpoin ) then
             !if( ( .not. lmark(jpoin) ) .and. ipoin /= jpoin ) then
             !
             ! NOT CLEAR WHAT IS BEST
             !
             if( ( .not. lmark(jpoin) ) .and. ipoin /= jpoin .and. jpoin <= npoi1 ) then
                vpoin(1:pdime) = coord(1:pdime,jpoin) - coord(1:pdime,ipoin)
                vrefe(1:pdime) = 0.5_rp * ( vecto(1:pdime,jpoin) + vecto(1:pdime,ipoin) )
                denom          = sqrt(dot_product(vpoin,vpoin)+epsil)*sqrt(dot_product(vrefe,vrefe)+epsil)
                cangl          = dot_product(vpoin,vrefe)/denom 
                if( cangl > cangm .and. cangl > angle_stream ) then                      
                   cangm = cangl
                   kpoin = jpoin 
                end if
             end if
          end do loop_izdom
          !
          ! Check if KPOIN is OK
          !
          if( kpoin /= 0 ) then
             ipoin        = kpoin
             nmark        = nmark + 1 
             lmark(ipoin) = .true.
             lgrou(ipoin) = ngrou
             invpr(nmark) = ipoin
             npoin_ngrou  = npoin_ngrou + 1
             !
             ! If KPOIN is outflow, go to next group
             !
             kbopo = lpoty(kpoin)
             if( kbopo > 0 ) then
                alpha = dot_product(exnor(1:pdime,1,kbopo),vecto(1:pdime,kpoin))
                if( alpha > 0.0_rp ) exit loop_igrou
             end if
          else
             exit loop_igrou
          end if
       end do loop_igrou
       !
       ! Goto next group
       !
       !          ngrou = ngrou + 1 !Guillaume
       !!          print*,'ngrou_mod',ngrou

       if( npoin_ngrou == 1 ) lgrou(ipoin) = -lgrou(ipoin) 

    end do
    !
    ! Save number of groups
    !
    !         ngrou = ngrou - 1 !Para los nodos internos, las boundaries van a parte y sino seria un grupo más cuando sale de grupo
    if( present(ngrou_gs) ) ngrou_gs = ngrou 
    !       print*,'ngrou',ngrou
    !
    ! PERMR: compute permutation array
    !
    do ipoin = 1,npoin
       jpoin = invpr(ipoin)
       permr(jpoin) = ipoin
    end do

    !-----------------------------------------------------------------
    !
    ! Comment following line (goto 10) to postprocess groups and stop
    !
    !-----------------------------------------------------------------

!!$       !if(ittim==10) then
!!$       !if( current_direction_neu /= 19 ) goto 10
!!$       goto 10
!!$       unit1 = (kfl_paral-1)*2+600
!!$       unit2 = (kfl_paral-1)*2+601
!!$       !      call postpr_right_now('PAULA','SCALA','NPOIN',lgrou)
!!$       open(unit=unit1,file='renumber-'//trim(intost(max(0_ip,kfl_paral)))//'.post.msh',status='unknown')
!!$       open(unit=unit2,file='renumber-'//trim(intost(max(0_ip,kfl_paral)))//'.post.res',status='unknown')
!!$       write(unit1,'(a,i1,a)') 'MESH RENUMBERING_GS dimension ',pdime,' Elemtype Linear Nnode 2'
!!$       jpoin = 0
!!$       write(unit1,*) 'coordinates'
!!$       do igrou = 1,ngrou
!!$          do kpoin = 1,npoin
!!$             ipoin = invpr(kpoin) 
!!$             if( abs(lgrou(ipoin)) == igrou ) then 
!!$                jpoin = jpoin + 1
!!$                write(unit1,'(i6,3(1x,e12.6))') jpoin,coord(1:pdime,ipoin)
!!$             end if
!!$          end do
!!$       end do
!!$       write(unit1,*) 'end coordinates'
!!$       jpoin = 0
!!$       jelem = 0
!!$       write(unit1,*) 'elements'
!!$       do igrou = 1,ngrou
!!$          tpoin = 0 
!!$          do kpoin = 1,npoin
!!$             ipoin = invpr(kpoin) 
!!$             if( abs(lgrou(ipoin)) == igrou ) then 
!!$                jpoin = jpoin + 1
!!$                tpoin = tpoin + 1
!!$                if( tpoin /= 1 ) then
!!$                   jelem = jelem + 1
!!$                   write(unit1,*) jelem,jpoin-1,jpoin
!!$                end if
!!$             end if
!!$          end do
!!$          if( tpoin == 1 ) then
!!$             jelem = jelem + 1
!!$             write(unit1,*) jelem,jpoin,jpoin             
!!$          end if
!!$       end do
!!$       write(unit1,*) 'end elements'
!!$       write(unit2,*) 'GiD Post Results File 1.0'
!!$       write(unit2,*) 'Result GROUPS_GS ALYA  0.00000000E+00 Scalar OnNodes'
!!$       write(unit2,*) 'ComponentNames GROUPS'
!!$       write(unit2,*) 'Values'
!!$       jpoin = 0
!!$       do igrou = 1,ngrou
!!$          do kpoin = 1,npoin
!!$             ipoin = invpr(kpoin) 
!!$             if( abs(lgrou(ipoin)) == igrou ) then 
!!$                jpoin = jpoin + 1
!!$                write(unit2,*) jpoin,lgrou(ipoin)
!!$             end if
!!$          end do
!!$       end do
!!$       write(unit2,*) 'end values'
!!$       write(unit2,*) 'Result VECTO_GS ALYA  0.00000000E+00 Vector OnNodes'
!!$       write(unit2,*) 'ComponentNames VECTO_GS_X,VECTO_GS_Y,VECTO_GS_Z'
!!$       write(unit2,*) 'Values'
!!$       jpoin = 0
!!$       do igrou = 1,ngrou
!!$          do kpoin = 1,npoin
!!$             ipoin = invpr(kpoin) 
!!$             if( abs(lgrou(ipoin)) == igrou ) then 
!!$                jpoin = jpoin + 1
!!$                write(unit2,*) jpoin,vecto(1:pdime,ipoin)
!!$             end if
!!$          end do
!!$       end do
!!$       write(unit2,*) 'end values'
!!$       close(unit=unit1)
!!$       close(unit=unit2)
!!$       !call runend('O.K.!')
!!$       !end if
!!$

10  continue

    !-----------------------------------------------------------------
    !
    ! Deallocate memory
    !
    !-----------------------------------------------------------------

    call memory_deallo(memor,'LORDE','graphs_numbering_along_vector',lorde)
    call memory_deallo(memor,'LMARK','graphs_numbering_along_vector',lmark)
    call memory_deallo(memor,'LINLE','graphs_numbering_along_vector',linle)
    call memory_deallo(memor,'VMODU','graphs_numbering_along_vector',vmodu)

    !end if

  end subroutine graphs_number_along_vector

  !------------------------------------------------------------------------
  !> @author  Guillaume Houzeaux
  !> @brief   List of subdomain connected to a node
  !> @details gives the list of subdomains DOMLI to which 
  !!          node ipoin belongs. The number of subdomains is NDOMI.
  !------------------------------------------------------------------------
  
  subroutine graphs_domlis(pelpo,lelpo,ipoin,domain,ndomi,domli)
  
    integer(ip), intent(in)  :: ipoin                   !< Node to inquire
    integer(ip), intent(in)  :: pelpo(*)                !< Node to element graph: pointer
    integer(ip), intent(in)  :: lelpo(*)                !< Node to element graph: list
    integer(ip), intent(in)  :: domain(*)               !< List of element subdomains
    integer(ip), intent(out) :: ndomi                   !< Number of subdomains node IPOIN belongs to
    integer(ip), intent(out) :: domli(*)                !< List of subdomain IPOIN belongs to
    integer(ip)              :: ii,jj,ielem,domin
    logical(lg)              :: liste

    ndomi = 0
    do ii = pelpo(ipoin),pelpo(ipoin+1)-1
       ielem = lelpo(ii)
       domin = domain(ielem)
       jj    = 1
       liste = .true.
       do while( jj <= ndomi .and. liste )
          if( domli(jj) == domin ) then
             liste = .false.
          end if
          jj = jj + 1
       end do
       if( liste ) then
          ndomi        = ndomi + 1
          domli(ndomi) = domin
       end if
    end do
  end subroutine graphs_domlis

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/01/1016
  !> @brief   Compute the list of edges
  !> @details Compute the list of edges
  !
  !-----------------------------------------------------------------------

  subroutine graphs_edges_type(meshe,memor)
    type(mesh_type), intent(inout) :: meshe          !< Mesh type
    integer(8),      intent(inout) :: memor(2)
    call graphs_edges_go(&
       meshe,meshe % medge,meshe % nedge,meshe % edge_to_node,&
       meshe % ledgs,meshe % lnned,meshe % ledgb,meshe % lnneb,&
       meshe % r_edg,meshe % c_edg,memor=memor)
    
  end subroutine graphs_edges_type

  subroutine graphs_edges_arrays(&
       meshe,medge,nedge,edge_to_node,&
       ledgs,lnned,ledgb,lnneb,r_edg,&
       c_edg,POINT_TO_TYPE,memor)
    type(mesh_type),                intent(inout) :: meshe          !< Mesh type
    integer(ip),                    intent(out)   :: medge
    integer(ip),                    intent(out)   :: nedge
    integer(ip),           pointer, intent(inout) :: edge_to_node(:,:)
    integer(ip),           pointer, intent(inout) :: ledgs(:,:)
    integer(ip),           pointer, intent(inout) :: lnned(:)
    integer(ip),           pointer, intent(inout) :: ledgb(:,:)
    integer(ip),           pointer, intent(inout) :: lnneb(:)
    integer(ip),           pointer, intent(inout) :: r_edg(:)
    integer(ip),           pointer, intent(inout) :: c_edg(:)
    logical(lg), optional,          intent(in)    :: POINT_TO_TYPE
    integer(8),                     intent(inout) :: memor(2)
    call graphs_edges_go(&
       meshe,medge,nedge,edge_to_node,&
       ledgs,lnned,ledgb,lnneb,r_edg,c_edg,memor=memor)
    if( present(POINT_TO_TYPE) ) then
       if( POINT_TO_TYPE ) then
          meshe % medge        =  medge
          meshe % nedge        =  nedge
          meshe % edge_to_node => edge_to_node
          meshe % ledgs        => ledgs
          meshe % lnned        => lnned
          meshe % ledgb        => ledgb
          meshe % lnneb        => lnneb
          meshe % r_edg        => r_edg
          meshe % c_edg        => c_edg 
       end if
    end if
    
  end subroutine graphs_edges_arrays
  
  subroutine graphs_edges_go(&
       meshe,medge,nedge,edge_to_node,&
       ledgs,lnned,ledgb,lnneb,r_edg,c_edg,&
       memor)
    
    type(mesh_type),           intent(inout) :: meshe          !< Mesh type
    integer(ip),               intent(out)   :: medge
    integer(ip),               intent(out)   :: nedge
    integer(ip),      pointer, intent(inout) :: edge_to_node(:,:)
    integer(ip),      pointer, intent(inout) :: ledgs(:,:)
    integer(ip),      pointer, intent(inout) :: lnned(:)
    integer(ip),      pointer, intent(inout) :: ledgb(:,:)
    integer(ip),      pointer, intent(inout) :: lnneb(:)
    integer(ip),      pointer, intent(inout) :: r_edg(:)
    integer(ip),      pointer, intent(inout) :: c_edg(:)
    integer(8),                intent(inout) :: memor(2)
    integer(ip)                              :: mpop2,ipoin
    integer(ip)                              :: lsize,iedgg
    integer(ip)                              :: ielem
    integer(ip)                              :: inode,ilisn
    integer(ip)                              :: nlelp,jpoin
    integer(ip)                              :: ipoin_1,ipoin_2
    integer(ip)                              :: jj,iedge,pelty
    integer(ip)                              :: pedge,ielpo
    integer(ip)                              :: iboun
    logical(lg)                              :: notfound

    integer(ip),      pointer                :: nepoi(:)
    integer(ip),      pointer                :: pelpo(:)
    integer(ip),      pointer                :: lelpo(:)
    integer(ip),      pointer                :: ledgp(:)
    integer(ip),      pointer                :: pedgp(:)

    integer(ip),                 pointer :: lnnod(:)
    integer(ip),                 pointer :: lnods(:,:)
    integer(ip),                 pointer :: lnnob(:)
    integer(ip),                 pointer :: lnodb(:,:)
    integer(ip),                 pointer :: ltype(:)
    integer(ip)                          :: nelem
    integer(ip)                          :: nboun
    integer(ip)                          :: npoin              

    nullify(nepoi)
    nullify(pelpo)
    nullify(lelpo)
    nullify(ledgp)
    nullify(pedgp)

    nullify(lnnod)
    nullify(lnods)
    nullify(lnnob)
    nullify(lnodb)
    nullify(ltype)

    lnnod => meshe % lnnod
    lnods => meshe % lnods
    lnnob => meshe % lnnob
    lnodb => meshe % lnodb
    ltype => meshe % ltype
    nelem =  meshe % nelem
    nboun =  meshe % nboun
    npoin =  meshe % npoin
    !
    ! Allocate memory for NEPOI and compute it
    ! 
    call memory_alloca(memor,'NEPOI','graphs_edges',nepoi,npoin)
    mpop2 = 0
    do ielem = 1,nelem
       mpop2 = mpop2 + lnnod(ielem)*lnnod(ielem)
       do inode = 1,lnnod(ielem)
          ipoin = lnods(inode,ielem)
          nepoi(ipoin) = nepoi(ipoin) + 1
       end do
    end do
    !
    ! Allocate memory for PELPO and compute it
    !
    call memory_alloca(memor,'PELPO','graphs_edges',pelpo,npoin+1_ip)
    pelpo(1) = 1
    do ipoin = 1,npoin
       pelpo(ipoin+1) = pelpo(ipoin) + nepoi(ipoin)
    end do
    !
    ! Allocate memory for LELPO and construct the list
    !
    nlelp = pelpo(npoin+1)
    call memory_alloca(memor,'LELPO','graphs_edges',lelpo,nlelp)
    do ielem = 1,nelem
       do inode = 1,lnnod(ielem)
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
    call memory_alloca(memor,'LEDGP','graphs_edges',ledgp,mpop2)
    call memory_alloca(memor,'PEDGP','graphs_edges',pedgp,npoin+1_ip)
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
    call memory_alloca(memor,'EDGE_TO_NODE','graphs_edges',edge_to_node,2_ip,nedge)
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
    ! Maximum number of edges
    !
    medge = 0
    do ielem = 1,nelem
       pelty = ltype(ielem)
       medge = max(medge,element_type(pelty) % number_edges)
    end do
    !
    ! LEDGS(1:PEDGE,1:NELEM) ... Fill in element to edge connectivity
    ! LLNED(IELEM) ............. Number of edges
    !
    call memory_alloca(memor,'LEDGS','graphs_edges',ledgs,medge,nelem)
    call memory_alloca(memor,'LNNED','graphs_edges',lnned,nelem)
    do ielem = 1,nelem
       pelty = ltype(ielem)
       pedge = element_type(pelty) % number_edges
       lnned(ielem) = pedge
       do iedge = 1,pedge
          ipoin_1 = lnods(element_type(pelty) % list_edges(1,iedge),ielem)
          ipoin_2 = lnods(element_type(pelty) % list_edges(2,iedge),ielem)
          if( ipoin_1 > ipoin_2 ) then
             ipoin   = ipoin_1
             ipoin_1 = ipoin_2
             ipoin_2 = ipoin
          end if
          iedgg_loop: do iedgg = pedgp(ipoin_1),pedgp(ipoin_1+1)-1                
             if( ledgp(iedgg) == ipoin_2 ) then
                ledgs(iedge,ielem) = iedgg
                exit iedgg_loop
             end if
          end do iedgg_loop
       end do
    end do
    !
    ! LEDGS(1:PEDGE,1:NELEM) ... Fill in element to edge connectivity
    ! LLNED(IELEM) ............. Number of edges
    !
    call memory_alloca(memor,'LEDGB','graphs_edges',ledgb,medge,nboun)
    call memory_alloca(memor,'LNNEB','graphs_edges',lnneb,nboun)
    do iboun = 1,nboun
       pelty = meshe % ltypb(iboun)
       pedge = element_type(pelty) % number_edges
       lnneb(iboun) = pedge
       do iedge = 1,pedge
          ipoin_1 = lnodb(element_type(pelty) % list_edges(1,iedge),iboun)
          ipoin_2 = lnodb(element_type(pelty) % list_edges(2,iedge),iboun)
          if( ipoin_1 > ipoin_2 ) then
             ipoin   = ipoin_1
             ipoin_1 = ipoin_2
             ipoin_2 = ipoin
          end if
          iedgb_loop: do iedgg = pedgp(ipoin_1),pedgp(ipoin_1+1)-1         
             if( ledgp(iedgg) == ipoin_2 ) then
                ledgb(iedge,iboun) = iedgg
                exit iedgb_loop
             end if
          end do iedgb_loop
       end do
    end do
    !
    ! Construct graph
    !
    call graphs_poipoi(&
         nedge,nelem,medge,ledgs,lnned,ltype,&
         r_edg,c_edg,memor=memor)
    !
    ! Deallocate memory
    !
    call memory_deallo(memor,'LEDGP','graphs_edges',ledgp)
    call memory_deallo(memor,'PEDGP','graphs_edges',pedgp)
    call memory_deallo(memor,'LELPO','graphs_edges',lelpo)
    call memory_deallo(memor,'PELPO','graphs_edges',pelpo)
    call memory_deallo(memor,'NEPOI','graphs_edges',nepoi)

  end subroutine graphs_edges_go

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    14/03/1016
  !> @brief   CSR to COO format
  !> @details Transform a CSR graph to an ELL graph
  !
  !-----------------------------------------------------------------------

  subroutine graphs_csr_to_coo(nn,ndofn,ia,ja,nz,row,column,message,memor)
    integer(ip),  intent(in)              :: nn
    integer(ip),  intent(in)              :: ndofn
    integer(ip),  intent(in)              :: ia(nn+1)
    integer(ip),  intent(in)              :: ja(*)
    integer(ip),  intent(inout), pointer  :: row(:)
    integer(ip),  intent(inout), pointer  :: column(:)
    integer(ip),  intent(out)             :: nz
    character(*), intent(in),    optional :: message
    integer(8),   intent(inout)           :: memor(2)
    integer(ip)                           :: ii,jj,iz,kz,idofn,jdofn
    integer(ip)                           :: ii_idofn,jj_jdofn,nnz
    logical(lg)                           :: symmetry
    !
    ! Symmetry
    !
    symmetry = .false.
    if( present(message) ) then
       if( trim(message) == 'SYMMETRIC' ) symmetry = .true.
    end if
    !
    ! Size of graph
    !
    nnz = ia(nn+1)-1
    if( symmetry ) then
       nz = ( nnz * ndofn * ndofn + nn * ndofn) / 2
    else
       nz = nnz * ndofn * ndofn
    end if
    !
    ! Allocate memory
    !
    if( .not. associated(row) ) then
       call memory_alloca(memor,'ROW'   ,'graphs_csr_to_coo',row,nz)
    else if( nz /= size(row,KIND=ip) ) then
       print*,'nz,size(row),kfl_paral',nz,size(row,KIND=ip) !,kfl_paral
       call runend('GRAPHS_CSR_TO_COO: TROUBLE IN SIZE')
    end if
    if( .not. associated(column) ) then
       call memory_alloca(memor,'COLUMN','graphs_csr_to_coo',column,nz)
    else if( nz /= size(column,KIND=ip) ) then
       print*,'nz,size(column),kfl_paral',nz,size(column,KIND=ip) !,kfl_paral
       call runend('GRAPHS_CSR_TO_COO: TROUBLE IN SIZE2')
    end if
    !
    ! Compute COO format including symmetric graph
    !
    if( ndofn == 1 ) then
       !
       ! One degree of freedom per node
       !
       if( symmetry ) then
          kz = 0
          do ii = 1,nn
             do iz = ia(ii),ia(ii+1)-1
                jj = ja(iz)
                if( jj >= ii ) then
                   kz         = kz + 1
                   row(kz)    = ii
                   column(kz) = jj
                end if
             end do
          end do
       else
          do ii = 1,nn
             do iz = ia(ii),ia(ii+1)-1
                jj         = ja(iz)
                row(iz)    = ii
                column(iz) = jj
             end do
          end do
       end if

    else if( ndofn > 1 ) then
       !
       ! Explode graph
       !
       if( symmetry ) then
          kz = 0
          do ii = 1,nn
             do idofn = 1,ndofn
                ii_idofn = (ii-1)*ndofn + idofn
                do iz = ia(ii),ia(ii+1)-1
                   jj = ja(iz)
                   do jdofn = 1,ndofn 
                      jj_jdofn  = (jj-1)*ndofn + jdofn
                      if( jj_jdofn >= ii_idofn ) then
                         kz        = kz + 1
                         row(kz)   = ii_idofn
                         column(kz)= jj_jdofn
                      end if
                   end do
                end do
             end do
          end do
       else
          kz = 0
          do ii = 1,nn
             do idofn = 1,ndofn
                ii_idofn = (ii-1)*ndofn + idofn
                do iz = ia(ii),ia(ii+1)-1
                   jj = ja(iz)
                   do jdofn = 1,ndofn 
                      jj_jdofn  = (jj-1)*ndofn + jdofn
                      kz        = kz + 1
                      row(kz)   = ii_idofn
                      column(kz)= jj_jdofn
                   end do
                end do
             end do
          end do
       end if
       if( kz /= nz ) then
          print*,'MATRIX SIZE=',kz,nz
          call runend('GRAPHS_CSR_TO_COO: PROBLEM WITH MATRIX SIZE')
       end if
    end if

  end subroutine graphs_csr_to_coo

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    14/03/1016
  !> @brief   XSR to COO format
  !> @details Transform a sum to a linked list
  !
  !-----------------------------------------------------------------------

  subroutine graphs_output(nn,ia,ja,coord,unit4,filename)
    
    integer(ip),  intent(in)            :: nn
    integer(ip),  intent(in)            :: ia(nn+1)
    integer(ip),  intent(in)            :: ja(*)
    real(rp),     intent(in), pointer   :: coord(:,:)
    integer(4),   intent(in)            :: unit4
    character(*), intent(in), optional  :: filename
    integer(ip)                         :: ii,jj,iz,pdime,ielem
    integer(4)                          :: unit1,unit2

    unit1 = unit4
    unit2 = unit4 + 1_4
    pdime = size(coord,1,KIND=ip)
    ielem = 0

    if( present(filename) ) then
       open(unit=unit1,file=trim(filename)//'.post.msh',status='unknown')
       open(unit=unit2,file=trim(filename)//'.post.res',status='unknown')       
    else
       open(unit=unit1,file='graph.post.msh',status='unknown')
       open(unit=unit2,file='graph.post.res',status='unknown')
    end if

    write(unit1,'(a)') 'MESH GRAPH dimension ',pdime,' Elemtype Linear Nnode 2'
    write(unit1,*) 'coordinates'
    do ii = 1,nn       
       write(unit1,'(i7,3(1x,e12.6))') ii,coord(1:pdime,ii)
    end do
    write(unit1,*) 'end coordinates'

    ielem = 0
    write(unit1,*) 'elements'
    do ii = 1,nn
       do iz = ia(ii),ia(ii+1)-1
          jj = ja(iz)
          if( ii /= jj ) then
             ielem = ielem + 1
             write(unit1,*) ielem,ii,jj      
          end if
       end do
    end do
    write(unit1,*) 'end elements'
    !
    ! Do not count myself
    !
    write(unit2,*) 'GiD Post Results File 1.0'
    write(unit2,*) 'Result NEIGHBORS ALYA  0.00000000E+00 Scalar OnNodes'
    write(unit2,*) 'ComponentNames NEIGHBORS'
    write(unit2,*) 'Values'
    do ii = 1,nn
       write(unit2,*) ii,ia(ii+1)-ia(ii)-1
    end do
    write(unit2,*) 'end values'

  end subroutine graphs_output


  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    14/03/1016
  !> @brief   CSR to COO format
  !> @details Transform a CSR graph to an ELL graph. Put NODFN > 1
  !>           to obtain an exploded graph
  !
  !-----------------------------------------------------------------------

  subroutine graphs_csr_to_ell(nn,ndofn,ia,ja,nz_ell,column,message,memor)
    integer(ip),  intent(in)              :: nn
    integer(ip),  intent(in)              :: ndofn
    integer(ip),  intent(in)              :: ia(nn+1)
    integer(ip),  intent(in)              :: ja(*)
    integer(ip),  intent(out)             :: nz_ell
    integer(ip),  intent(inout), pointer  :: column(:,:)
    character(*), intent(in),    optional :: message
    integer(8),   intent(inout)           :: memor(2)
    integer(ip)                           :: ii,jj,iz,icol
    !
    ! Size
    !
    nz_ell = 0
    do ii = 1,nn
       iz     = ia(ii+1)-ia(ii)
       nz_ell = max(nz_ell,iz)
    end do
    !
    ! Allocate
    !
    if( .not. associated(column) ) then
       call memory_alloca(memor,'COLUMN','graphs_csr_to_ell',column,nz_ell,nn)
    else
       call runend('GRAPHS_CSR_TO_ELL: GRAPH ALREADY ASSOCIATED')
    end if
    !
    ! Format
    !
    ! +-             -+    +-             -+
    ! |  1 3 6 7 0 0  |    |  1 3 6 7 7 7  |
    ! |  2 4 6 7 8 0  |    |  2 4 6 7 8 8  |
    ! |  3 4 5 6 0 0  | => |  3 4 5 6 6 6  |
    ! |  ...          |    |  ...          |
    ! +-             -+    +-             -+
    ! 
    do ii = 1,nn
       icol = 0
       do iz = ia(ii),ia(ii+1)-1
          jj = ja(iz)
          icol = icol + 1
          column(icol,ii) = jj
       end do
       !
       ! Fill in remaning part with last node instead of zeros
       !
       do iz = icol+1,nz_ell
          column(iz,ii) = jj
       end do
    end do
    !
    ! Total size
    !
    nz_ell = nz_ell * nn

  end subroutine graphs_csr_to_ell

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-06-01
  !> @brief   Compute a Dual graph
  !> @details Compute the dual of a graph
  !>
  !>          \verbatim
  !>
  !>            A simple example to illustrate arrays:
  !>
  !>            +----+------------+
  !>            |    |      2     |
  !>            | 1  +------+-----+
  !>            |    |      |     |
  !>            +----+      |  3  |
  !>            |           |     |
  !>            |     4     |     |
  !>            |           |     |
  !>            +-----------+-----+
  !>
  !>            - IADUAL and JADUAL: Dual graph
  !>            - NBDUAL: number of nodes of the dual graph
  !>            - TRANSDUAL: from adjacancies to edges of the dual graph
  !>
  !>                   1
  !>            (1)---------(2)
  !>              \         /|
  !>                \   3 /  |
  !>               2  \ /    | 4            
  !>                  / \    |              
  !>                /     \  |
  !>              /         \|
  !>             (3)--------(4)         
  !>                    5
  !>
  !>            NBDUAL = 5
  !>            I = 1 ... JADUAL(IADUAL(I):IADUAL(ii+1)-1) = 3,4,2
  !>            I = 2 ... JADUAL(IADUAL(I):IADUAL(ii+1)-1) = 4,5,1
  !>            I = 3 ... JADUAL(IADUAL(I):IADUAL(ii+1)-1) = 5,1,4
  !>            I = 4 ... JADUAL(IADUAL(I):IADUAL(ii+1)-1) = 2,5,1,3
  !>            I = 5 ... JADUAL(IADUAL(I):IADUAL(ii+1)-1) = 2,4,3
  !>
  !>            I = 1 ... JA(IA(I+0)) = 2   =>   TRANSDUAL(IA(I+0)) = 1 : subd1 and subd2: edge 1
  !>            I = 1 ... JA(IA(I+1)) = 4   =>   TRANSDUAL(IA(I+1)) = 2 : subd1 and subd2: edge 2
  !>            I = 2 ... JA(IA(I+0)) = 1   =>   TRANSDUAL(IA(I+0)) = 1 : subd2 and subd1: edge 1
  !>            I = 2 ... JA(IA(I+1)) = 3   =>   TRANSDUAL(IA(I+1)) = 3 : subd2 and subd3: edge 3
  !>            I = 2 ... JA(IA(I+2)) = 4   =>   TRANSDUAL(IA(I+2)) = 4 : subd2 and subd4: edge 4
  !>            I = 3 ... JA(IA(I+0)) = 2   =>   TRANSDUAL(IA(I+0)) = 3 : subd3 and subd2: edge 3
  !>            I = 3 ... JA(IA(I+1)) = 4   =>   TRANSDUAL(IA(I+1)) = 5 : subd3 and subd4: edge 5
  !>            I = 4 ... JA(IA(I+0)) = 1   =>   TRANSDUAL(IA(I+0)) = 2 : subd4 and subd2: edge 2
  !>            I = 4 ... JA(IA(I+1)) = 2   =>   TRANSDUAL(IA(I+1)) = 4 : subd4 and subd2: edge 4
  !>            I = 4 ... JA(IA(I+2)) = 3   =>   TRANSDUAL(IA(I+2)) = 5 : subd4 and subd3: edge 5
  !>
  !>          \endverbatim 
  !>
  !-----------------------------------------------------------------------

  subroutine graphs_dual_graph_deallocate(iaDual,jaDual,translDual,memor)

    integer(ip), intent(inout), pointer, optional  :: iaDual(:)
    integer(ip), intent(inout), pointer, optional  :: jaDual(:)
    integer(ip), intent(inout), pointer, optional  :: translDual(:)
    integer(8),  intent(inout)                     :: memor(2)

    if( present(iadual) ) then
       call memory_deallo(memor,'IADUAL'    ,'graphs_dual_graph',iadual)
    end if
    if( present(jadual) ) then
       call memory_deallo(memor,'jADUAL'    ,'graphs_dual_graph',jadual)
    end if
    if( present(translDual) ) then
       call memory_deallo(memor,'TRANSLDUAL','graphs_dual_graph',translDual)
    end if
    
  end subroutine graphs_dual_graph_deallocate
  
  subroutine graphs_dual_graph(nn,ia,ja,nnDual,iaDual,jaDual,translDual,memor)

    integer(ip), intent(in)              :: nn
    integer(ip), intent(in),    pointer  :: ia(:)
    integer(ip), intent(in),    pointer  :: ja(:)
    integer(ip), intent(out)             :: nnDual
    integer(ip), intent(inout), pointer  :: iaDual(:)
    integer(ip), intent(inout), pointer  :: jaDual(:)
    integer(ip), intent(inout), pointer  :: translDual(:)
    integer(8),  intent(inout)           :: memor(2)
    integer(ip)                          :: ii, jj, vv, ww, t1, t2
    logical(lg)                          :: fin
    integer(ip), pointer                 :: iwa(:)

    nullify(iwa)

    !
    ! Construir vector para asignar un id a cada arista y asi traducir las 
    ! aristas como nodos del nuevo grafo 
    !
    call memory_alloca(memor,'translDual','graphs_dual_graph',translDual,ia(nn+1)-1)

    nnDual = 0
    do vv = 1, nn
       do ii = ia(vv), ia(vv+1)-1
          ww = ja(ii)
          if (vv < ww) then
             nnDual = nnDual + 1
             translDual(ii) = nnDual
          else
             fin = .false.
             jj  = ia(ww)
             do while ((jj < ia(ww+1)) .and. (.not. fin))
                if (ja(jj) == vv) then
                   translDual(ii) = translDual(jj)
                   fin = .true.             
                endif
                jj = jj + 1
             end do
          end if
       end do
    end do
    !
    ! Compute size of jaDual (nodes)
    !
    call memory_alloca(memor,'iwa','graphs_dual_graph',iwa,nnDual)

    do ii = 1, nnDual
       iwa(ii) = 0
    enddo

    do ii = 1, ia(nn+1)-1
       vv = ja(ii)
       t1 = translDual(ii)
       do jj = ia(vv), ia(vv+1)-1
          if (translDual(jj) /= t1) then
             iwa(t1) = iwa(t1) + 1
          endif
       enddo
    enddo
    !
    ! Construct iaDual (nodes)
    !
    call memory_alloca(memor,'iaDual','graphs_dual_graph',iaDual,nnDual+1_ip)

    iaDual(1) = 1
    do ii = 1, nnDual
       iaDual(ii+1)   = iaDual(ii) + iwa(ii)
       iwa(ii) = iaDual(ii)
    enddo
    !
    ! Construct jaDual (not ordered)
    !
    call memory_alloca(memor,'jaDual','graphs_dual_graph',jaDual,iaDual(nnDual+1)-1)

    do ii = 1,ia(nn+1)-1
       vv = ja(ii)
       t1 = translDual(ii)
       do jj = ia(vv),ia(vv+1)-1
          t2 = translDual(jj)
          if (t1 /= t2) then
             jaDual(iwa(t1)) = t2
             iwa(t1)         = iwa(t1) + 1
          end if
       end do
    end do

    call memory_deallo(memor,'iwa'       ,'graphs_dual_graph',iwa)

  end subroutine graphs_dual_graph

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-06-01
  !> @brief   Color a graph
  !>          using the minimum number of colors so that each node of the
  !>          dual graph does not have the same color as its neighbors.
  !>          Remember that a node of the dual graph is an edge connecting
  !>          two adjacent subdomains. 
  !>
  !>          \verbatim
  !> 
  !>            A simple example to illustrate arrays:
  !>
  !>            +----+------------+
  !>            |    |      2     |
  !>            | 1  +------+-----+
  !>            |    |      |     |
  !>            +----+      |  3  |
  !>            |           |     |
  !>            |     4     |     |
  !>            |           |     |
  !>            +-----------+-----+
  !>
  !>            NBCOLORS: Number of colors
  !>            COLOURS: Colors of the dual graph
  !>            Colors are in brackets near the node of the dual graph:
  !>
  !>                   1[1]
  !>            (1)---------(2)
  !>              \    3[2] /|
  !>                \     /  |
  !>            2[2]  \ /    | 4[3]            
  !>                  / \    |              
  !>                /     \  |
  !>              /         \|
  !>            (3)---------(4)        
  !>                  5[1]
  !>
  !>            NBCOLORS = 3
  !>            I = 1 ... COLORS(I) = 1
  !>            I = 2 ... COLORS(I) = 2
  !>            I = 3 ... COLORS(I) = 2
  !>            I = 4 ... COLORS(I) = 3
  !>            I = 5 ... COLORS(I) = 1
  !>
  !>          \endverbatim
  !>
  !------------------------------------------------------------------------
  
  subroutine graphs_color_graph_deallocate(colors,memor)

    integer(ip), intent(inout), pointer, optional :: colors(:)
    integer(8),  intent(inout)                    :: memor(2)

    call memory_deallo(memor,'colors','graphs_color_graph_deallocate',colors)

  end subroutine graphs_color_graph_deallocate
  
  subroutine graphs_color_graph(nn,ia,ja,nncolors,colors,memor)

    integer(ip), intent(in)                 :: nn
    integer(ip), intent(in),    pointer     :: ia(:)
    integer(ip), intent(in),    pointer     :: ja(:)
    integer(ip), intent(out)                :: nncolors
    integer(ip), intent(inout), pointer     :: colors(:)
    integer(8),  intent(out)                :: memor(2)
    integer(ip),                parameter   :: maxColors = 1000_ip
    integer(ip)                             :: colNode, ii, vv, node,nn0,ifact
    logical(lg)                             :: colourFound
    integer(ip),                pointer     :: SortVec_node(:)
    integer(ip),                pointer     :: SortVec_nbAdj(:)
    integer(ip),                pointer     :: mask(:)
    integer(ip)                             :: umask

    if( nn == 0 ) then

       nncolors = 0

    else

       nullify(mask)
       nullify(SortVec_node)
       nullify(SortVec_nbAdj)

       ifact = 1_ip

       call memory_alloca(memor,'colors'       ,'graphs_color_graph',colors,nn)
       call memory_alloca(memor,'SortVec_node' ,'graphs_color_graph',SortVec_node, nn)
       call memory_alloca(memor,'SortVec_nbAdj','graphs_color_graph',SortVec_nbAdj,nn)

1      continue
       call memory_alloca(memor,'mask'         ,'graphs_color_graph',mask,ifact*maxColors+1_ip,lboun=0_ip)
       umask = ifact*maxColors

       do vv = 1,nn
          sortVec_node(vv)  = vv
          sortVec_nbAdj(vv) = ia(vv+1) - ia(vv)
       end do
       !
       ! Order sortVec, de momento burbuja
       !
       nn0 = nn
       call maths_heap_sort(2_ip,nn0,sortVec_nbAdj,ivo1=sortVec_node)
       !
       ! Initialize color for each dual graph node
       !
       colors   = 0
       nncolors = 0
       !
       ! Use the minimum number of colors so that each node of the graph does not have the same
       ! color as its neighbors
       !
       mask     = 0
       !
       ! Subdomain 1 has color 1... just to start with something!
       !
       nncolors                = 1
       colors(sortVec_node(1)) = 1

       do vv = 2,nn

          node = sortVec_node(vv)

          do ii = ia(node),ia(node+1)-1
             ! Should check that colors(ja(ii)) < maxColors: very unlikely to happen but we never know!...
             if( colors(ja(ii)) > umask ) then
                !write(*,*) "colors(ja(ii)) ",colors(ja(ii))
                call memory_deallo(memor,'mask'         ,'graphs_color_graph',mask)
                ifact = 2
                goto 1
             end if
             mask(colors(ja(ii))) = vv
          end do

          colourFound = .false.
          ii          = 1
          do while( ii <= nncolors .and. .not. colourFound )
             if( mask(ii) /= vv ) then
                colNode     = ii
                colourFound = .true.
             end if
             ii = ii + 1
          end do

          if ( .not. colourFound ) then
             nncolors = nncolors + 1
             colNode  = nncolors
          end if

          colors(sortVec_node(vv)) = colNode

       end do

       call memory_deallo(memor,'mask'         ,'graphs_color_graph',mask)
       call memory_deallo(memor,'SortVec_node' ,'graphs_color_graph',SortVec_node)
       call memory_deallo(memor,'SortVec_nbadj','graphs_color_graph',SortVec_nbadj)


    end if

  end subroutine graphs_color_graph

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-28
  !> @brief   Faces
  !> @details Compute element faces
  !> 
  !-----------------------------------------------------------------------
  
  subroutine graphs_face_type(nelem,ltype,pelpo,lelpo,nedge,faces)

    integer(ip),          intent(in)    :: nelem
    integer(ip), pointer, intent(in)    :: ltype(:)
    integer(ip), pointer, intent(in)    :: pelpo(:)
    integer(ip), pointer, intent(in)    :: lelpo(:)
    integer(ip),          intent(out)   :: nedge
    integer(ip), pointer, intent(inout) :: faces(:,:,:)
    integer(ip)                         :: ielem,ielty,ipoin
    integer(ip)                         :: pnodi,ilist,iface
    integer(ip)                         :: ielpo,jelem,jelty,jface
    integer(ip)                         :: inodb,pnodj,pepoi,ilast
    logical(lg)                         :: equal_faces
    !
    ! Compute FACES
    !
    ilast = size(faces,1)
    nedge = 0
    do ielem = 1,nelem                                          ! Compare the faces and 
       ielty = abs(ltype(ielem))                                ! eliminate the repited faces
       do iface=1,element_type(ielty) % number_faces
          ipoin = faces(1,iface,ielem)
          if( ipoin /= 0 ) then
             pnodi = element_type(ielty) % node_faces(iface)
             ilist = 1
             pepoi = pelpo(ipoin+1)-pelpo(ipoin)
             ielpo = pelpo(ipoin)-1
             do while( ilist <= pepoi )
                ielpo = ielpo+1
                jelem = lelpo(ielpo)
                if( jelem /= ielem ) then
                   jelty = abs(ltype(jelem))                              ! eliminate the repited faces
                   jface = 0
                   do while( jface /= element_type(jelty) % number_faces )
                      jface = jface + 1
                      if( faces(1,jface,jelem) /= 0 ) then
                         equal_faces = .true.
                         inodb       = 0
                         pnodj       = element_type(jelty) % node_faces(jface)
                         do while( equal_faces .and. inodb /= pnodj )
                            inodb = inodb + 1
                            if( faces(inodb,iface,ielem) /= faces(inodb,jface,jelem) ) &
                                 equal_faces = .false.
                         end do
                         if( equal_faces ) then
                            faces(    1,iface,ielem) = 0                                  ! IFACE and JFACE
                            faces(    1,jface,jelem) = 0                                  ! are eliminated
                            faces(ilast,iface,ielem) = jelem                              ! IFACE and JFACE
                            faces(ilast,jface,jelem) = ielem                              ! are eliminated
                            nedge                    = nedge + 2
                            jface                    = element_type(jelty) % number_faces ! Exit JFACE do
                            ilist                    = pepoi                              ! Exit JELEM do
                         end if
                      end if
                   end do
                end if
                ilist = ilist + 1
             end do
          end if
       end do
    end do

  end subroutine graphs_face_type

end module mod_graphs
!> @}
