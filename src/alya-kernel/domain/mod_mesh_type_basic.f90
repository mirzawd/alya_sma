!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @defgroup Mesh_Type_Toolbox
!> @{
!> @name    ToolBox for mesh type
!> @file    mod_mesh_type.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for basic mesh type
!> @details ToolBox for basic mesh type
!
!-----------------------------------------------------------------------

module mod_mesh_type_basic

  use def_kintyp_basic,      only : ip,rp,lg,i1p,r2p
  use def_kintyp_mesh,       only : mesh_type_basic
  use def_kintyp_mesh,       only : mesh_type
  use def_kintyp_comm,       only : comm_data_par_basic
  use def_maths_tree,        only : maths_octree
  use mod_maths,             only : maths_in_box
  use mod_maths,             only : maths_heap_sort
  use def_master,            only : cutim  
  use def_master,            only : ittim  
  use def_master,            only : namda  
  use def_master,            only : kfl_paral  
  use def_master,            only : npart
  use def_master,            only : INOTMASTER
  use def_master,            only : IPARALL
  use def_elmtyp,            only : element_max
  use mod_elmgeo,            only : element_type
  use mod_memory,            only : memory_alloca
  use mod_memory,            only : memory_deallo
  use mod_memory,            only : memory_copy
  use mod_parall,            only : PAR_COMM_MY_CODE
  use mod_parall,            only : PAR_COMM_WORLD
  use mod_parall,            only : PAR_MY_CODE_RANK
  use mod_parall,            only : commd
  use def_kintyp_basic,      only : ip,rp,lg
  use def_kintyp_mesh_basic, only : mesh_type_basic
  use def_kintyp_mesh_basic, only : NODE_TAG
  use def_domain,            only : memor_dom
  use mod_memory,            only : memory_alloca
  use mod_memory,            only : memory_deallo
  use mod_postpr,            only : postpr
  use mod_iofile,            only : iofile_open_unit
  use mod_iofile,            only : iofile_close_unit
  use mod_iofile,            only : iofile_available_unit
  use mod_iofile,            only : iofile_create_directory
  use mod_communications,    only : PAR_COMM_SPLIT 
  use mod_communications,    only : PAR_COMM_RANK_AND_SIZE
  use mod_communications,    only : PAR_ALLGATHER
  use mod_communications,    only : PAR_GATHER
  use mod_communications,    only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications,    only : PAR_MAX
  use mod_communications,    only : PAR_SEND_RECEIVE
  use mod_communications,    only : PAR_ALLTOALL
  use mod_communications,    only : PAR_BROADCAST
  use mod_arrays,            only : arrays_tag
  use mod_graphs,            only : graphs_eleele_faces
  use mod_graphs,            only : graphs_eleele_deallocate
  use mod_graphs,            only : graphs_number_to_linked_list
  use mod_clusters,          only : clusters_basic
  use mod_maths_basic,       only : maths_mapping_coord_to_3d
  use def_maths_bin,         only : maths_bin
  use def_maths_octbin,      only : maths_octbin
  use mod_messages,          only : messages_live
  use mod_optional_argument, only : optional_argument
  use def_mpi
#include "def_mpi.inc"
  implicit none

  private
  character(19), parameter :: vacal='mod_mesh_type_basic'

  
  public :: mesh_type_basic_copy_boundary_mesh        ! Copy a boundary mesh from the extended mesh format to basic format
  public :: mesh_type_basic_output                    ! Output of the domain of a basic mesh
  public :: mesh_type_basic_parallel                  ! Setup the parallel communication arrays
  public :: mesh_type_basic_valid_mesh                ! Convert a mesh to a valid mesh (with elements conncted through faces)
  public :: mesh_type_basic_global_numbering          ! Global numbering
  public :: mesh_type_basic_parallel_interface        ! Parallel interface
  public :: mesh_type_basic_send_recv                 ! Send and receive a mesh
  public :: mesh_type_basic_broadcast                 ! Broadcast a mesh
  public :: mesh_type_basic_octree_mesh               ! Create an octree mesh
  public :: mesh_type_basic_bin_mesh                  ! Create a bin mesh
  public :: mesh_type_basic_octbin_mesh               ! Create an octbin mesh
 
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-14
  !> @brief   Extract a boundary mesh
  !> @details Copy a boundary mesh from the extended mesh format to
  !>          basic format
  !> 
  !-----------------------------------------------------------------------
  
  subroutine mesh_type_basic_copy_boundary_mesh(mesh_bou,mesh,MEMORY_COUNTER,lmask)

    type(mesh_type_basic),                    intent(inout) :: mesh_bou
    type(mesh_type),                          intent(in)    :: mesh
    integer(8),            optional,          intent(inout) :: MEMORY_COUNTER(2)   !< Memory counter
    logical(lg),           optional, pointer, intent(inout) :: lmask(:)
    integer(ip)                                             :: kpoin,iboun,ipoin
    integer(ip)                                             :: inode,pelty,pnode
    integer(ip)                                             :: kelem,mnode_loc
    integer(ip),                     pointer                :: permr_nodes(:)
    integer(8)                                              :: memor_loc(2)
    logical(lg)                                             :: if_boundary

    nullify(permr_nodes)

    if( present(MEMORY_COUNTER) ) then
       memor_loc = MEMORY_COUNTER
    else
       memor_loc = 0_8
    end if
    !  
    ! Permutation
    !
    call memory_alloca(memor_loc,'PERMR_NODES','extract',permr_nodes,mesh % npoin)
    kpoin     = 0
    kelem     = 0
    mnode_loc = 0
    if( mesh % nboun > 0 ) mnode_loc = mesh % mnodb
    do iboun = 1,mesh % nboun
       if( present(lmask) ) then
          if_boundary = lmask(iboun)
       else
          if_boundary = .true.
       end if
       if( if_boundary ) then
          kelem     = kelem + 1
          pelty     = mesh % ltypb(iboun)
          pnode     = element_type(pelty) % number_nodes
          !mnode_loc = max(pnode,mnode_loc)
          do inode = 1,pnode
             ipoin = mesh % lnodb(inode,iboun)
             if( permr_nodes(ipoin) == 0 ) then
                kpoin              = kpoin + 1
                permr_nodes(ipoin) = kpoin
             end if
          end do
       end if
    end do
    !
    ! New mesh
    !
    if( trim(mesh_bou % name) == '' ) &
         mesh_bou % name  = trim(mesh % name) // '_extract'
    mesh_bou % ndime = mesh % ndime
    mesh_bou % nelem = kelem
    mesh_bou % npoin = kpoin
    mesh_bou % mnode = mnode_loc
    call mesh_bou % alloca(MEMORY_COUNTER=MEMORY_COUNTER)

    kelem = 0
    do iboun = 1,mesh % nboun 
       if( present(lmask) ) then
          if_boundary = lmask(iboun)
       else
          if_boundary = .true.
       end if
       if( if_boundary ) then
          pelty                           = mesh % ltypb(iboun)
          pnode                           = element_type(pelty) % number_nodes
          kelem                           = kelem + 1
          mesh_bou % ltype(kelem)         = mesh % ltypb    (iboun)
          mesh_bou % leinv_loc(kelem)     = mesh % lbinv_loc(iboun)             
          mesh_bou % lnods(1:pnode,kelem) = permr_nodes(mesh % lnodb(1:pnode,iboun))
       end if
    end do

    do ipoin = 1,mesh % npoin
       kpoin = permr_nodes(ipoin)
       if( kpoin /= 0 ) mesh_bou % permn(kpoin) = ipoin
    end do

    do kpoin = 1,mesh_bou % npoin
       ipoin                       = mesh_bou % permn(kpoin)
       mesh_bou % coord(:,kpoin)   = mesh % coord(:,ipoin)  
       mesh_bou % lninv_loc(kpoin) = mesh % lninv_loc(ipoin)
    end do

    call memory_deallo(memor_loc,'PERMR_NODES','extract',permr_nodes)

    if( present(MEMORY_COUNTER) )  MEMORY_COUNTER = memor_loc

  end subroutine mesh_type_basic_copy_boundary_mesh
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-14
  !> @brief   Output mesh
  !> @details Output mesh in mpio format
  !> 
  !-----------------------------------------------------------------------
  
  subroutine mesh_type_basic_output(mesh,MESH_OUTPUT)
    
    use def_master,          only : INOTSLAVE
    use mod_communications,  only : PAR_BARRIER

    class(mesh_type_basic),           intent(inout) :: mesh
    logical(lg),            optional, intent(in)    :: MESH_OUTPUT
    character(5)                                    :: wopos(5)
    integer(ip)                                     :: nsubd
    integer(ip),            pointer                 :: nelem_gat(:)
    integer(ip),            pointer                 :: npoin_gat(:) 
    integer(ip)                                     :: nunit,isubd

    nullify(nelem_gat)
    nullify(npoin_gat)

    if (INOTSLAVE) call iofile_create_directory(mesh % name)      
    call PAR_BARRIER()

    if( mesh % comm % rank4 >= 0_4 ) then
       !
       ! Parallelization log file
       !
       nsubd = int(mesh % comm % size4,ip)       
       if( nsubd > 0 ) then
          call memory_alloca(memor_dom,'NELEM_GAT','mod_domain_memory',nelem_gat,nsubd,lboun=0_ip)
          call memory_alloca(memor_dom,'NPOIN_GAT','mod_domain_memory',npoin_gat,nsubd,lboun=0_ip)
          call PAR_GATHER(mesh % nelem,nelem_gat,PAR_COMM_IN=mesh % comm % PAR_COMM_WORLD)
          call PAR_GATHER(mesh % npoin,npoin_gat,PAR_COMM_IN=mesh % comm % PAR_COMM_WORLD)
       end if

       if( mesh % comm % rank4 == 0_4 ) then

          nunit = iofile_available_unit()

          call iofile_open_unit(nunit,trim(mesh % name)//'/'//trim(namda)//'-'//trim(mesh % name)//'.post.alyapar','POSTPROCESS LOG FILE')          
          write(nunit,1) max(1_ip,nsubd)
          if( nsubd == 0 ) then
             write(nunit,2) &
                  &         1,&
                  &         mesh % nelem,&
                  &         mesh % npoin,&
                  &         0_ip
          else
             do isubd = 0,max(1_ip,nsubd)-1
                write(nunit,2) &
                     &         isubd+1,&
                     &         nelem_gat(isubd),&
                     &         npoin_gat(isubd),&
                     &         0_ip
             end do
          end if
          call iofile_close_unit(nunit)

       end if
       call memory_deallo(memor_dom,'NELEM_GAT','mod_domain_memory',nelem_gat)
       call memory_deallo(memor_dom,'NPOIN_GAT','mod_domain_memory',npoin_gat)

       if( optional_argument(.true.,MESH_OUTPUT) ) then
          !
          ! LTYPE  
          !
          wopos = arrays_tag('LTYPE',MODULE_NUMBER=0_ip)
          call postpr(mesh % ltype,wopos,ittim,cutim,MESH=mesh)
          !
          ! LNODS
          !
          wopos = arrays_tag('LNODS',MODULE_NUMBER=0_ip)
          call postpr(mesh % lnods,wopos,ittim,cutim,mesh % mnode,MESH=mesh)
          !
          ! LEINV_LOC
          !
          wopos = arrays_tag('LEINV',MODULE_NUMBER=0_ip)
          call postpr(mesh % leinv_loc,wopos,ittim,cutim,MESH=mesh) 
          !
          ! LNINV_LOC
          !
          wopos = arrays_tag('LNINV',MODULE_NUMBER=0_ip)
          call postpr(mesh % lninv_loc,wopos,ittim,cutim,MESH=mesh) 
          !
          ! COORD
          !
          wopos = arrays_tag('COORD',MODULE_NUMBER=0_ip)
          call postpr(mesh % coord,wopos,ittim,cutim,mesh % ndime,MESH=mesh) 

       end if

    end if

1   format(i9)
2   format(10(1x,i9))

  end subroutine mesh_type_basic_output
    
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-14
  !> @brief   Set parallel communication
  !> @details Set the parallel communcation for submeshes in the
  !>          original communicator
  !> 
  !-----------------------------------------------------------------------

  subroutine mesh_type_basic_parallel(mesh,PAR_COMM_INITIAL,PAR_RANK_INITIAL,COMM_INITIAL)

    class(mesh_type_basic),                      intent(inout) :: mesh
    MY_MPI_COMM   ,            optional,         intent(in)    :: PAR_COMM_INITIAL
    integer(ip),               optional,         intent(in)    :: PAR_RANK_INITIAL
    type(comm_data_par_basic), optional, target, intent(in)    :: COMM_INITIAL
    integer(ip)                                                :: comm_size
    integer(ip)                                                :: my_rank
    MY_MPI_COMM                                                :: PAR_COMM
    integer(ip)                                                :: icolor
    integer(ip)                                                :: ipoin,ielem,kpoin
    integer(ip)                                                :: offset_npoin
    integer(ip)                                                :: iparent
    integer(ip), pointer                                       :: npoin_gat(:)
    integer(ip), pointer                                       :: nelem_gat(:)
    integer(ip), pointer                                       :: owner(:)
    integer(ip)                                                :: my_npoin
    integer(ip)                                                :: parent_npoin
    MY_MPI_COMM                                                :: PAR_COMM_PARENT
    integer(ip)                                                :: PAR_RANK_PARENT
    type(comm_data_par_basic), pointer                         :: COMM_PARENT

    nullify(npoin_gat,nelem_gat,owner)
    nullify(COMM_PARENT)

    !--------------------------------------------------------------------
    !
    ! Parent communicator, rank and communication arrays
    !
    !--------------------------------------------------------------------
    !
    ! Check if Parent exists
    !
    if( associated(mesh % parent) ) then
       iparent      = 1
    else
       iparent      = 0
    end if
    !
    ! Guess initial/parent communiactor
    !
    if( present(PAR_COMM_INITIAL) ) then
       PAR_COMM_PARENT = PAR_COMM_INITIAL
    else if( iparent == 1 ) then
       PAR_COMM_PARENT = mesh % parent % comm % PAR_COMM_WORLD
    else
       PAR_COMM_PARENT = PAR_COMM_MY_CODE
    end if

    if( present(PAR_RANK_INITIAL) ) then
       PAR_RANK_PARENT = PAR_RANK_INITIAL
    else if( iparent == 1 ) then
       PAR_RANK_PARENT = mesh % parent % comm % RANK4
    else
       PAR_RANK_PARENT = PAR_MY_CODE_RANK
    end if

    if( present(COMM_INITIAL) ) then
       COMM_PARENT => COMM_INITIAL
    else if( iparent == 1 ) then
       COMM_PARENT => mesh % parent % comm
    end if

    if( PAR_RANK_PARENT >= 0 ) then
       
       !--------------------------------------------------------------------
       !
       ! Split communicator by consider only meshes with elements and get:
       !
       ! MESH % COMM % PAR_COMM_WORLD
       ! MESH % COMM % RANK4
       ! MESH % COMM % SIZE4
       !
       !--------------------------------------------------------------------
       
       icolor = min(mesh % nelem,1_ip)
       call PAR_COMM_SPLIT(icolor,PAR_COMM,my_rank,PAR_COMM_PARENT,PAR_RANK_PARENT)

       mesh % comm % PAR_COMM_WORLD = PAR_COMM
       mesh % comm % RANK4          = int(my_rank,4)

       if( my_rank >= 0 ) then
          call PAR_COMM_RANK_AND_SIZE(mesh % comm % PAR_COMM_WORLD,my_rank,comm_size)
          mesh % comm % SIZE4 = int(comm_size,4)            
       end if

       !--------------------------------------------------------------------
       !
       ! Node global numbering and offset
       !
       ! MESH % LNINV_LOC
       ! MESH % COMM % OFFSET_NPOIN
       !
       !--------------------------------------------------------------------
       
       if( my_rank >= 0 ) then
          call memory_alloca(memor_dom,'NPOIN_GAT','par_submsh',npoin_gat,comm_size,lboun=0_ip)
          call PAR_ALLGATHER(mesh % npoin,npoin_gat,PAR_COMM_IN=mesh % comm % PAR_COMM_WORLD)
          mesh % comm % offset_npoin = sum(npoin_gat(0:mesh % comm % rank4-1))
          call memory_deallo(memor_dom,'NPOIN_GAT','par_submsh',npoin_gat)
       end if

     if( iparent == 1 ) then
          !
          ! Define own interface nodes
          !
          parent_npoin = mesh % parent % npoin
          call memory_alloca(memor_dom,'OWNER','par_submsh',owner,parent_npoin)
          do kpoin = 1,parent_npoin
             owner(kpoin) = -1
          end do
         do ipoin = 1,mesh % npoin
             kpoin = mesh % permn(ipoin)
             owner(kpoin) = mesh % comm % RANK4
          end do
          call PAR_INTERFACE_NODE_EXCHANGE(owner,'MAX',COMM_PARENT)

          my_npoin = 0
          do ipoin = 1,mesh % npoin
             kpoin = mesh % permn(ipoin)
             if( owner(kpoin) == mesh % comm % RANK4 ) then
                my_npoin     = my_npoin + 1
             else
                owner(kpoin) = -1
             end if
          end do

          if( my_rank >= 0 ) then        
             call memory_alloca(memor_dom,'NPOIN_GAT','par_submsh',npoin_gat,comm_size,lboun=0_ip)
             call PAR_ALLGATHER(my_npoin,npoin_gat,PAR_COMM_IN=mesh % comm % PAR_COMM_WORLD)
             offset_npoin = sum(npoin_gat(0:mesh % comm % rank4-1))
             call memory_deallo(memor_dom,'NPOIN_GAT','par_submsh',npoin_gat)
             my_npoin = 0
             do ipoin = 1,mesh % npoin
                kpoin = mesh % permn(ipoin)
                if( owner(kpoin) == mesh % comm % RANK4 ) then
                   my_npoin     = my_npoin + 1
                   owner(kpoin) = my_npoin + offset_npoin
                end if
             end do
          end if
         
          call PAR_INTERFACE_NODE_EXCHANGE(owner,'MAX',COMM_PARENT)
          !
          ! Global numbering by eliminating duplicated nodes
          !
          do ipoin = 1,mesh % npoin
             kpoin = mesh % permn(ipoin)
             mesh % lninv_loc(ipoin) = owner(kpoin)
          end do
          call memory_deallo(memor_dom,'OWNER','par_submsh',owner)

       else if( my_rank >= 0 ) then
          !
          ! Do not eliminate interface nodes
          !       
          do ipoin = 1,mesh % npoin
             mesh % lninv_loc(ipoin) = ipoin + mesh % comm % offset_npoin
          end do

       end if

       !--------------------------------------------------------------------
       !
       ! Element global numbering and offset
       !
       ! MESH % LEINV_LOC
       ! MESH % COMM % OFFSET_NELEM
       !
       !--------------------------------------------------------------------

       if( my_rank >= 0 ) then
          !
          ! Compute offsets
          !
          call memory_alloca(memor_dom,'NELEM_GAT','par_submsh',nelem_gat,comm_size,lboun=0_ip)
          call PAR_ALLGATHER(mesh % nelem,nelem_gat,PAR_COMM_IN=mesh % comm % PAR_COMM_WORLD)
          mesh % comm % offset_nelem = sum(nelem_gat(0:mesh % comm % rank4-1))
          call memory_deallo(memor_dom,'NELEM_GAT','par_submsh',nelem_gat)
          !
          ! Change global numbering
          !    
          do ielem = 1,mesh % nelem
             mesh % leinv_loc(ielem) = ielem + mesh % comm % offset_nelem
          end do

       end if

    end if

  end subroutine mesh_type_basic_parallel

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-14
  !> @brief   Set parallel communication
  !> @details Set the parallel communcation for submeshes in the
  !>          original communicator
  !> 
  !-----------------------------------------------------------------------

  subroutine mesh_type_basic_global_numbering(mesh)

    class(mesh_type_basic),  intent(inout) :: mesh
    integer(ip)                            :: comm_size
    integer(ip)                            :: my_rank
    integer(ip)                            :: ipoin,ielem
    integer(ip), pointer                   :: npoin_gat(:)
    integer(ip), pointer                   :: nelem_gat(:)
    integer(ip), pointer                   :: owner(:)
    integer(ip)                            :: my_npoin

    nullify(npoin_gat,nelem_gat,owner)

    my_rank   = int(mesh % comm % RANK4,ip)
    comm_size = int(mesh % comm % SIZE4,ip)

    if( my_rank >= 0 ) then

       !--------------------------------------------------------------------
       !
       ! Node global numbering and offset
       !
       ! MESH % LNINV_LOC
       ! MESH % COMM % OFFSET_NPOIN
       !
       !--------------------------------------------------------------------
       !
       ! Define own interface nodes
       !
       call memory_alloca(memor_dom,'OWNER','par_submsh',owner,mesh % npoin)
       do ipoin = 1,mesh % npoin 
          owner(ipoin) = my_rank 
       end do

       call PAR_INTERFACE_NODE_EXCHANGE(owner,'MAX',mesh % comm)       

       my_npoin = 0
       do ipoin = 1,mesh % npoin
          if( owner(ipoin) == my_rank ) then
             my_npoin     = my_npoin + 1
          else
             owner(ipoin) = -1
          end if
       end do

       call memory_alloca(memor_dom,'NPOIN_GAT','par_submsh',npoin_gat,comm_size,lboun=0_ip)
       call PAR_ALLGATHER(my_npoin,npoin_gat,PAR_COMM_IN=mesh % comm % PAR_COMM_WORLD)
       mesh % comm % offset_npoin = sum(npoin_gat(0:my_rank-1))
       call memory_deallo(memor_dom,'NPOIN_GAT','par_submsh',npoin_gat)

       my_npoin = 0
       do ipoin = 1,mesh % npoin
          if( owner(ipoin) == my_rank ) then
             my_npoin     = my_npoin + 1
             owner(ipoin) = my_npoin + mesh % comm % offset_npoin
          end if
       end do

       call PAR_INTERFACE_NODE_EXCHANGE(owner,'MAX',mesh % comm)
       !
       ! Global numbering by eliminating duplicated nodes
       !
       do ipoin = 1,mesh % npoin
          mesh % lninv_loc(ipoin) = owner(ipoin)
       end do
       call memory_deallo(memor_dom,'OWNER','par_submsh',owner)

       !--------------------------------------------------------------------
       !
       ! Element global numbering and offset
       !
       ! MESH % LEINV_LOC
       ! MESH % COMM % OFFSET_NELEM
       !
       !--------------------------------------------------------------------

       call memory_alloca(memor_dom,'NELEM_GAT','par_submsh',nelem_gat,comm_size,lboun=0_ip)
       call PAR_ALLGATHER(mesh % nelem,nelem_gat,PAR_COMM_IN=mesh % comm % PAR_COMM_WORLD)
       mesh % comm % offset_nelem = sum(nelem_gat(0:mesh % comm % rank4-1))
       call memory_deallo(memor_dom,'NELEM_GAT','par_submsh',nelem_gat)

       do ielem = 1,mesh % nelem
          mesh % leinv_loc(ielem) = ielem + mesh % comm % offset_nelem
       end do

    end if

  end subroutine mesh_type_basic_global_numbering

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-28
  !> @brief   Valid mesh
  !> @details Extract a valid mesh by eliminating elements not connected
  !           through faces
  !> 
  !-----------------------------------------------------------------------

  subroutine mesh_type_basic_valid_mesh(mesh,mesh_not)

    type(mesh_type_basic),           intent(inout)   :: mesh
    type(mesh_type_basic), optional, intent(inout)   :: mesh_not
    integer(ip)                                      :: nclus,nn,ii
    integer(ip)                                      :: iclus,iclus_max
    integer(ip)                                      :: nn_max,nn_clus
    integer(ip),           pointer                   :: lclus(:)
    integer(ip),           pointer                   :: pelel(:)
    integer(ip),           pointer                   :: lelel(:)
    logical(lg),           pointer                   :: lmask(:)
    type(mesh_type_basic)                            :: mesh_cpy

    nullify(pelel,lelel,lclus,lmask)

    nn = mesh % nelem
    if( nn > 0 ) then
       !
       ! Eelemtn graph using common faces
       !
       call graphs_eleele_faces(mesh,PELEL=pelel,LELEL=lelel,memor=memor_dom)   
       !
       ! Look for clusters
       !
       call memory_alloca(memor_dom,'LCLUS','mesh_type_basic_valid_mesh',lclus,nn)
       call clusters_basic(nn,pelel,lelel,nclus,lclus)    
       call graphs_eleele_deallocate(pelel,lelel,memor=memor_dom)
       
       if( nclus > 1 ) then
          !
          ! Mark elements belonging to the larger cluster
          !
          iclus_max = 0
          nn_max    = 0
          do iclus = 1,nclus
             nn_clus = count(lclus==iclus)
             if( nn_clus > nn_max ) then
                iclus_max = iclus
                nn_max = nn_clus
             end if
         end do
          allocate(lmask(nn))
          do ii = 1,nn
             if( lclus(ii) == iclus_max ) then
                lmask(ii) = .true.
             else
                lmask(ii) = .false.
             end if
          end do
          call memory_deallo(memor_dom,'LCLUS','mesh_type_basic_valid_mesh',lclus)
          !
          ! Remove other clusters in mesh
          !
          call mesh_cpy % init() 
          call mesh_cpy % copy(mesh)
          call mesh     % deallo()
          call mesh     % extract(mesh_cpy,lmask)
          if( present(mesh_not) ) then
             lmask = .not. lmask
             call mesh_not % extract(mesh_cpy,lmask)
          end if
          call mesh_cpy % deallo()
          deallocate(lmask)

       end if
       
       call memory_deallo(memor_dom,'LCLUS','mesh_type_basic_valid_mesh',lclus)

    end if

  end subroutine mesh_type_basic_valid_mesh

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-08
  !> @brief   Create communication
  !> @details Create a basic communicator using global numbering
  !> 
  !-----------------------------------------------------------------------

  subroutine mesh_type_basic_parallel_interface(mesh)


    class(mesh_type_basic),  intent(inout) :: mesh
    integer(ip)                           :: COMM_SIZE
    MY_MPI_COMM                           :: PAR_COMM4
    integer(ip)                           :: npoin,ipart,ipoin
    integer(ip)                           :: kpoin,ndime,my_rank
    integer(ip)                           :: ii,kk,ineig
    real(rp)                              :: comin(3),comax(3)
    integer(ip)                           :: lnmin,lnmax
    real(rp)                              :: delta(3)
    real(rp),               pointer       :: coord(:,:)
    real(rp),               pointer       :: par_comin(:,:)
    real(rp),               pointer       :: par_comax(:,:)
    integer(ip),            pointer       :: par_lnmin(:)
    integer(ip),            pointer       :: par_lnmax(:)
    integer(ip),            pointer       :: npoin_send(:)
    integer(ip),            pointer       :: npoin_recv(:)
    integer(ip),            pointer       :: num_inte(:)
    integer(ip),            pointer       :: neights(:)
    integer(ip),            pointer       :: bound_size(:)
    type(i1p),              pointer       :: bound_perm(:)
    type(i1p),              pointer       :: lninv_send(:)
    type(i1p),              pointer       :: lninv_recv(:)
    logical(lg),            pointer       :: lmask(:)
    logical(lg),            pointer       :: linte(:)

    call messages_live('COMMUNICATION ARRAYS')
    
    nullify(par_comin,par_comax,npoin_send,npoin_recv,lninv_send,lninv_recv,coord)
    nullify(par_lnmin,par_lnmax)

    nullify(lmask)
    nullify(num_inte,linte,neights,bound_perm,bound_size)

    call memory_alloca(memor_dom,'LMASK',vacal,lmask,mesh % npoin)
    !
    ! Get possible boundary nodes
    !
    call mesh % boundary_nodes(lmask)

    coord               => mesh % coord    
    ndime               =  mesh % ndime
    npoin               =  mesh % npoin
    COMM_SIZE           =  int(mesh % comm % SIZE4,ip)
    PAR_COMM4           =  mesh % comm % PAR_COMM_WORLD
    my_rank             =  int(mesh % comm % RANK4,ip) !kfl_paral
    mesh % comm % nneig =  0
    !
    ! Bounding boxes
    !
    comin =  huge(1.0_rp)*0.1_rp
    comax = -huge(1.0_rp)*0.1_rp
    lnmin =  huge(1_ip)
    lnmax = -huge(1_ip)
    do ipoin = 1,npoin
       if( lmask(ipoin) ) then
          comin(1:ndime) = min(comin(1:ndime),coord(1:ndime,ipoin))
          comax(1:ndime) = max(comax(1:ndime),coord(1:ndime,ipoin))
          lnmin          = min(lnmin,mesh % lninv_loc(ipoin))
          lnmax          = max(lnmax,mesh % lninv_loc(ipoin))
       end if
    end do
    delta = (comax - comin) * 1.0e-3_rp    
    comin = comin - delta - epsilon(1.0_rp)
    comax = comax + delta + epsilon(1.0_rp)

    call memory_alloca(memor_dom,'PAR_COMIN',vacal,par_comin,3_ip,COMM_SIZE,'INITIALIZE',1_ip,0_ip)   
    call memory_alloca(memor_dom,'PAR_COMAX',vacal,par_comax,3_ip,COMM_SIZE,'INITIALIZE',1_ip,0_ip)    
    
    call memory_alloca(memor_dom,'PAR_LNMIN',vacal,par_lnmin,COMM_SIZE,'INITIALIZE',0_ip)   
    call memory_alloca(memor_dom,'PAR_LNMAX',vacal,par_lnmax,COMM_SIZE,'INITIALIZE',0_ip)    
    
    call PAR_ALLGATHER(comin,par_comin,3_4,PAR_COMM_IN=PAR_COMM4)
    call PAR_ALLGATHER(comax,par_comax,3_4,PAR_COMM_IN=PAR_COMM4)
    call PAR_ALLGATHER(lnmin,par_lnmin,    PAR_COMM_IN=PAR_COMM4)
    call PAR_ALLGATHER(lnmax,par_lnmax,    PAR_COMM_IN=PAR_COMM4)
    !
    ! Number of points to send
    !
    call memory_alloca(memor_dom,'NPOIN_SEND',vacal,npoin_send,COMM_SIZE,'INITIALIZE',0_ip)   
    call memory_alloca(memor_dom,'NPOIN_RECV',vacal,npoin_recv,COMM_SIZE,'INITIALIZE',0_ip)
    do ipoin = 1,npoin
       if( lmask(ipoin) ) then
          do ipart = 1,COMM_SIZE-1
             if( mesh % lninv_loc(ipoin) >= par_lnmin(ipart) .and. mesh % lninv_loc(ipoin) <= par_lnmax(ipart) ) then
                if( maths_in_box(mesh % ndime,coord(:,ipoin),par_comin(:,ipart),par_comax(:,ipart)) ) then
                   npoin_send(ipart) = npoin_send(ipart) + 1
                end if
             end if
          end do
       end if
    end do
    call PAR_ALLTOALL(1_ip,1_ip,npoin_send,npoin_recv,PAR_COMM_IN=PAR_COMM4)
    !
    ! Fill in numbering to send
    !
    call memory_alloca(memor_dom,'LNINV_SEND',vacal,lninv_send,COMM_SIZE,'INITIALIZE',lboun=0_ip)   
    call memory_alloca(memor_dom,'LNINV_RECV',vacal,lninv_recv,COMM_SIZE,'INITIALIZE',lboun=0_ip)

    do ipart = 1,COMM_SIZE-1
       call memory_alloca(memor_dom,'LNINV_SEND % L',vacal,lninv_send(ipart) % l,npoin_send(ipart))      
       call memory_alloca(memor_dom,'LNINV_RECV % L',vacal,lninv_recv(ipart) % l,npoin_recv(ipart))      
    end do

    do ipart = 1,COMM_SIZE-1
       kpoin = 0
       do ipoin = 1,npoin
          if( lmask(ipoin) ) then
             if( mesh % lninv_loc(ipoin) >= par_lnmin(ipart) .and. mesh % lninv_loc(ipoin) <= par_lnmax(ipart) ) then
                if( maths_in_box(mesh % ndime,coord(:,ipoin),par_comin(:,ipart),par_comax(:,ipart)) ) then
                   kpoin = kpoin + 1
                   lninv_send(ipart) % l(kpoin) = mesh % lninv_loc(ipoin)
                end if
             end if
          end if
       end do
    end do
    !
    ! Send receive coordinates
    !
    do ipart = 1,COMM_SIZE-1
       call PAR_SEND_RECEIVE(lninv_send(ipart) % l,lninv_recv(ipart) % l,DOM_I=ipart,PAR_COMM_IN=PAR_COMM4)       
    end do
    !
    ! Identify interface nodes
    !
    call memory_alloca(memor_dom,'NUM_INTE'  ,vacal,num_inte  ,COMM_SIZE  ,'INITIALIZE',lboun=0_ip)   
    call memory_alloca(memor_dom,'BOUND_PERM',vacal,bound_perm,COMM_SIZE  ,'INITIALIZE',lboun=0_ip)   
    call memory_alloca(memor_dom,'BOUND_SIZE',vacal,bound_size,COMM_SIZE  ,'INITIALIZE',lboun=0_ip)   
    call memory_alloca(memor_dom,'NEIGHTS'   ,vacal,neights   ,COMM_SIZE+1,'INITIALIZE') 
    call memory_alloca(memor_dom,'LINTE'     ,vacal,linte,npoin)

    call mesh % comm % deallo(memor_dom)
    call mesh % comm % init  (COMM_NAME=trim(mesh % name)//' % COMM')
    
    !call PAR_DEALLOCATE_COMMUNICATION_ARRAY(mesh % comm,memor_dom)
    !call PAR_INITIALIZE_COMMUNICATION_ARRAY(mesh % comm,COMM_NAME=trim(mesh % name)//' % COMM')

    do ipart = 1,COMM_SIZE-1
       if( ipart /= my_rank ) then
          do ipoin = 1,npoin
             linte(ipoin) = .false.
          end do
          do kk = 1,npoin_recv(ipart)
             kpoin = lninv_recv(ipart) % l(kk)
             do ipoin = 1,npoin
                if( lmask(ipoin) ) then
                   if( mesh % lninv_loc(ipoin) == kpoin ) then
                      num_inte(ipart) = num_inte(ipart) + 1
                      linte(ipoin) = .true.
                   end if
                end if
             end do
          end do
          if( num_inte(ipart) > 0 ) then
             mesh % comm % nneig          = mesh % comm % nneig + 1
             neights(mesh % comm % nneig) = ipart
             bound_size(ipart)            = count(linte) 
             call memory_alloca(memor_dom,'BOUND_PERM % L',vacal,bound_perm(ipart) % l,bound_size(ipart))
             kpoin = 0
             do ipoin = 1,npoin
                if( linte(ipoin) ) then
                   kpoin = kpoin + 1
                   bound_perm(ipart) % l(kpoin) = ipoin 
                end if
             end do
             if( kpoin > 0 ) call maths_heap_sort(2_ip,kpoin,bound_perm(ipart) % l,PERMUTATION=mesh % lninv_loc)         
          end if
       end if
    end do
    !
    ! Construct communication 
    !
    mesh % comm % PAR_COMM_WORLD = PAR_COMM_WORLD
    mesh % comm % RANK4          = int(kfl_paral,4_ip)
    mesh % comm % SIZE4          = npart+1
    mesh % comm % bound_dim      = sum(bound_size)
    mesh % comm % nneig          = count(neights/=0)
    call mesh % comm % alloca(memor_dom)
    !call PAR_ALLOCATE_COMMUNICATION_ARRAY(mesh % comm,memor_dom)!,COMM_NAME=trim(mesh % name)//' % COMM')
    kk = 0
    do ineig = 1,mesh % comm % nneig
       ipart                           = neights(ineig)
       mesh % comm % neights(ineig)    = ipart
       mesh % comm % bound_size(ineig) = bound_size(ipart)
       do ii = 1,bound_size(ipart)
          kk = kk + 1
          mesh % comm % bound_perm(kk) = bound_perm(ipart) % l(ii)
       end do
    end do
    if( mesh % comm % bound_dim > 0 ) call graphs_number_to_linked_list(mesh % comm % nneig,mesh % comm % bound_size) 
    !
    ! Deallocate
    !
    call memory_deallo(memor_dom,'LMASK'     ,vacal,lmask)
    call memory_deallo(memor_dom,'NUM_INTE'  ,vacal,num_inte)   
    call memory_deallo(memor_dom,'NEIGHTS'   ,vacal,neights) 
    call memory_deallo(memor_dom,'BOUND_PERM',vacal,bound_perm)   
    call memory_deallo(memor_dom,'BOUND_SIZE',vacal,bound_size)   
    call memory_deallo(memor_dom,'LINTE'     ,vacal,linte)
    call memory_deallo(memor_dom,'PAR_COMIN' ,vacal,par_comin)   
    call memory_deallo(memor_dom,'PAR_COMAX' ,vacal,par_comax)    
    call memory_deallo(memor_dom,'PAR_LNMIN' ,vacal,par_lnmin)   
    call memory_deallo(memor_dom,'PAR_LNMAX' ,vacal,par_lnmax)    
    call memory_deallo(memor_dom,'NPOIN_SEND',vacal,npoin_send)   
    call memory_deallo(memor_dom,'NPOIN_RECV',vacal,npoin_recv)
    call memory_deallo(memor_dom,'LNINV_SEND',vacal,lninv_send)   
    call memory_deallo(memor_dom,'LNINV_RECV',vacal,lninv_recv)

  end subroutine mesh_type_basic_parallel_interface
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-08
  !> @brief   Send receive
  !> @details Send and receive a mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine mesh_type_basic_send_recv(mesh_send,mesh_recv,ipart,PAR_COMM4)
    
    type(mesh_type_basic),  intent(inout) :: mesh_send
    type(mesh_type_basic),  intent(inout) :: mesh_recv
    integer(ip),            intent(in)    :: ipart
    MY_MPI_COMM   ,         intent(in)    :: PAR_COMM4
    integer(ip)                           :: itag

    call PAR_SEND_RECEIVE(mesh_send % coord    ,mesh_recv % coord,    'IN MY CODE',ipart,PAR_COMM_IN=PAR_COMM4)                   
    call PAR_SEND_RECEIVE(mesh_send % lninv_loc,mesh_recv % lninv_loc,'IN MY CODE',ipart,PAR_COMM_IN=PAR_COMM4)                            
    call PAR_SEND_RECEIVE(mesh_send % lnods    ,mesh_recv % lnods,    'IN MY CODE',ipart,PAR_COMM_IN=PAR_COMM4)                   
    call PAR_SEND_RECEIVE(mesh_send % ltype    ,mesh_recv % ltype,    'IN MY CODE',ipart,PAR_COMM_IN=PAR_COMM4)                   
    call PAR_SEND_RECEIVE(mesh_send % leinv_loc,mesh_recv % leinv_loc,'IN MY CODE',ipart,PAR_COMM_IN=PAR_COMM4)

    do itag = 1,mesh_send % ntags
       call PAR_SEND_RECEIVE(mesh_send % tags(itag) % values,mesh_recv % tags(itag) % values,'IN MY CODE',ipart,PAR_COMM_IN=PAR_COMM4)
    end do

  end subroutine mesh_type_basic_send_recv

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-08
  !> @brief   Send receive
  !> @details Send and receive a mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine mesh_type_basic_broadcast(mesh)
    
    type(mesh_type_basic),  intent(inout) :: mesh
    integer(ip)                           :: itag,ielty

    call PAR_BROADCAST(mesh % npoin    )                   
    call PAR_BROADCAST(mesh % nelem    )                   
    call PAR_BROADCAST(mesh % mnode    )                   
    call PAR_BROADCAST(mesh % ndime    )                   
    call PAR_BROADCAST(mesh % ntags    )                   

    if( INOTMASTER ) then
       call mesh % alloca()       
    end if

    call PAR_BROADCAST(mesh % coord    )                   
    call PAR_BROADCAST(mesh % lninv_loc)                            
    call PAR_BROADCAST(mesh % lnods    )                   
    call PAR_BROADCAST(mesh % ltype    )                   
    call PAR_BROADCAST(mesh % leinv_loc)

    do itag = 1,mesh % ntags
       call PAR_BROADCAST(mesh % tags(itag) % type  )
       if( INOTMASTER ) call mesh % alloca_tag(TAG_ID=itag)
       call PAR_BROADCAST(mesh % tags(itag) % values)
    end do
    !
    ! Quadrature and integration rule
    !
    do ielty = 1,size(mesh % quad)
       call PAR_BROADCAST(mesh % quad(ielty) % type)
       call PAR_BROADCAST(mesh % quad(ielty) % ngaus)
    end do

    if( INOTMASTER ) then
       do ielty = 1,size(element_type)
          if( mesh % quad(ielty) % ngaus > 0 ) then
             mesh % quad(ielty) % topo  = element_type(ielty) % topology
             mesh % quad(ielty) % ndime = element_type(ielty) % dimensions
             call mesh % quad(ielty) % set()
             call mesh % iso(ielty)  % set(&
                  element_type(ielty) % number_nodes,&
                  mesh % quad(ielty))
          end if
       end do
    end if
    
    
  end subroutine mesh_type_basic_broadcast

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-08
  !> @brief   Create an octree mesh
  !> @details Create an octree mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine mesh_type_basic_octree_mesh(mesh_oct,mesh,limit,OFFSET,CENTROID)
    
    type(mesh_type_basic),                   intent(inout) :: mesh_oct
    type(mesh_type_basic),                   intent(in)    :: mesh
    integer(ip),                             intent(in)    :: limit
    real(rp),             optional,          intent(in)    :: OFFSET
    real(rp),             optional, pointer, intent(inout) :: CENTROID(:,:)
    type(maths_octree)                                     :: octree 

    call octree   % init      ()
    call octree   % input     (LIMIT=limit)
    call octree   % fill      (mesh % coord)
    call octree   % mesh_dim  (mesh_oct % ndime,mesh_oct % mnode,mesh_oct % nelem,mesh_oct % npoin)
    call mesh_oct % alloca    ()
    call octree   % mesh      (mesh_oct % ndime,mesh_oct % mnode,mesh_oct % nelem,mesh_oct % npoin,&         
         &                     mesh_oct % lnods,mesh_oct % ltype,mesh_oct % coord,CENTROID=CENTROID)
    call octree % deallo      ()
    
  end subroutine mesh_type_basic_octree_mesh
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-08
  !> @brief   Create an bin mesh
  !> @details Create an bin mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine mesh_type_basic_bin_mesh(mesh_bin,mesh,boxes,OFFSET,CENTROID)
    
    type(mesh_type_basic),                   intent(inout) :: mesh_bin
    type(mesh_type_basic),                   intent(in)    :: mesh
    integer(ip),                             intent(in)    :: boxes(:)
    real(rp),             optional,          intent(in)    :: OFFSET
    real(rp),             optional, pointer, intent(inout) :: CENTROID(:,:)
    type(maths_bin)                                        :: bin

    call bin      % init      ()
    call bin      % input     (BOXES=boxes,RELATIVE_TOLERANCE=OFFSET)
    call bin      % fill      (mesh % coord)
    call bin      % mesh_dim  (mesh_bin % ndime,mesh_bin % mnode,mesh_bin % nelem,mesh_bin % npoin)
    call mesh_bin % alloca    ()
    call bin      % mesh      (mesh_bin % ndime,mesh_bin % mnode,mesh_bin % nelem,mesh_bin % npoin,&
         &                     mesh_bin % lnods,mesh_bin % ltype,mesh_bin % coord,CENTROID=CENTROID)
    call bin      % deallo    ()
    
  end subroutine mesh_type_basic_bin_mesh
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-08
  !> @brief   Create an octbin mesh
  !> @details Create an octbin mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine mesh_type_basic_octbin_mesh(mesh_octbin,mesh,boxes,limit,OFFSET,CENTROID)
    
    type(mesh_type_basic),                   intent(inout) :: mesh_octbin
    type(mesh_type_basic),                   intent(in)    :: mesh
    integer(ip),                             intent(in)    :: boxes(:)
    integer(ip),                             intent(in)    :: limit
    real(rp),             optional,          intent(in)    :: OFFSET
    real(rp),             optional, pointer, intent(inout) :: CENTROID(:,:)
    type(maths_octbin)                                     :: octbin

    call octbin      % init      ()
    call octbin      % input     (BOXES=boxes,LIMIT=limit,RELATIVE_TOLERANCE=OFFSET,ENABLE_GAP=.false.)
    call octbin      % fill      (mesh % coord)
    call octbin      % mesh_dim  (mesh_octbin % ndime,mesh_octbin % mnode,mesh_octbin % nelem,mesh_octbin % npoin)
    call mesh_octbin % alloca    ()
    call octbin      % mesh      (mesh_octbin % ndime,mesh_octbin % mnode,mesh_octbin % nelem,mesh_octbin % npoin,&
         &                        mesh_octbin % lnods,mesh_octbin % ltype,mesh_octbin % coord,CENTROID=CENTROID)
    call octbin      % deallo    ()

    
  end subroutine mesh_type_basic_octbin_mesh

 end module mod_mesh_type_basic
!> @}
