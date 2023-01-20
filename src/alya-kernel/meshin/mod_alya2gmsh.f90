!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




!-----------------------------------------------------------------------
!> @addtogroup Meshing
!> @{
!> @file    mod_alya2gmsh.f90
!> @author  houzeaux
!> @date    2020-05-26
!> @brief   Module to gmsh
!> @details Functions to link Alya to gmsh
!-----------------------------------------------------------------------

module mod_alya2gmsh

  use  def_kintyp_basic, only : ip,rp,lg
  use  def_kintyp_mesh,  only : mesh_type_basic
  use  mod_elmgeo,       only : element_type
  use  def_domain,       only : ndime
  use  def_elmtyp
  use  def_domain, only : ndime  
  use  mod_memory_basic
  use  mod_std
  use, intrinsic :: iso_c_binding
  
  implicit none
  private  

  type :: gmsh_type
     integer(8)                        :: nnodes8                ! Number of nodes
     integer(8),            pointer    :: nodes_tags8(:)         ! Node global numbering
     real(8),               pointer    :: nodes_coordinates8(:)  ! Node coordinates
     integer(8)                        :: ntypes8                ! Number of types
     integer(8),            pointer    :: types_tags8(:)         ! List of types
     integer(8),            pointer    :: nelements_bytype8(:)   ! Number of elements by type
     integer(8),            pointer    :: elements_tags8(:)      ! Element global numbering
     integer(8),            pointer    :: elements_nodes8(:)     ! Element connectivity
     integer(8)                        :: nelem8                 ! Number of elements: Not needed by gmsh
     integer(8)                        :: nelem_nodes8           ! Size of element connectivity: Not needed by gmsh
   contains
     procedure,             pass       :: init                   ! Initialize 
     procedure,             pass       :: alloca                 ! Allocate
     procedure,             pass       :: deallo                 ! Deallocate
  end type gmsh_type

  type :: gmsh_type_c
     integer(C_SIZE_T)                 :: nnodes8                ! Number of nodes
     type(C_PTR)                       :: nodes_tags8            ! Node global numbering
     type(C_PTR)                       :: nodes_coordinates8     ! Node coordinates
     integer(C_SIZE_T)                 :: ntypes8                ! Number of types
     type(C_PTR)                       :: types_tags8            ! List of types
     type(C_PTR)                       :: nelements_bytype8      ! Number of elements by type
     type(C_PTR)                       :: elements_tags8         ! Element global numbering
     type(C_PTR)                       :: elements_nodes8        ! Element connectivity
     integer(C_SIZE_T)                 :: nelem8                 ! Number of elements: Not needed by gmsh
     integer(C_SIZE_T)                 :: nelem_nodes8           ! Size of element connectivity: Not needed by gmsh
   contains
     procedure,             pass       :: deallo_c                ! Deallocate
  end type gmsh_type_c

  integer(ip)                         :: permr_a2g(100)         ! Alya to gmsh element type
  integer(ip)                         :: permr_g2a(100)         ! Gmsh to Alya element type
  integer(8)                          :: memor(2)               ! Local memory counter
  character(13), parameter            :: vacal='mod_alya2gmsh'

  INTERFACE

     !-----------------------------------------------------------------------
     !> 
     !> @author  houzeaux
     !> @date    2020-05-04
     !> @brief   Interface with c wrapper
     !> @details Interface with c wrapper to gmsh
     !> 
     !-----------------------------------------------------------------------

     SUBROUTINE reMeshing(&
          &    ndime,                                             &
          &    nnodes_vol, nodes_tags_vol, nodes_coordinates_vol, &
          &    ntypes_vol,types_tags_vol, nelements_bytype_vol,   & 
          &    elements_tags_vol,elements_nodes_vol,              &
          &    nnodes_sur, nodes_tags_sur, nodes_coordinates_sur, &
          &    ntypes_sur,types_tags_sur, nelements_bytype_sur,   & 
          &    elements_tags_sur,elements_nodes_sur,              &
          &    algorithm_3D, max_anisotropy,                      &
          &    min_size, max_size, mesh_size,                     &                    
          &    new_nnodes, new_nodes_tags,new_nodes_coordinates,  &
          &    new_ntypes, new_types_tags, new_nelements_bytype,  &
          &    new_elements_tags,new_elements_nodes, ipass)       BIND (c,NAME='reMeshing')

       use, intrinsic ::  ISO_C_BINDING
       implicit none       
       integer(C_INT),        value             :: ndime
       
       integer(C_SIZE_T),     value             :: nnodes_vol
       integer(C_SIZE_T),     dimension(*)      :: nodes_tags_vol
       real   (C_DOUBLE),     dimension(*)      :: nodes_coordinates_vol
       integer(C_SIZE_T),     value             :: ntypes_vol
       integer(C_SIZE_T),     dimension(*)      :: types_tags_vol
       integer(C_SIZE_T),     dimension(*)      :: nelements_bytype_vol
       integer(C_SIZE_T),     dimension(*)      :: elements_tags_vol
       integer(C_SIZE_T),     dimension(*)      :: elements_nodes_vol

       integer(C_SIZE_T),     value             :: nnodes_sur
       integer(C_SIZE_T),     dimension(*)      :: nodes_tags_sur
       real   (C_DOUBLE),     dimension(*)      :: nodes_coordinates_sur
       integer(C_SIZE_T),     value             :: ntypes_sur
       integer(C_SIZE_T),     dimension(*)      :: types_tags_sur
       integer(C_SIZE_T),     dimension(*)      :: nelements_bytype_sur
       integer(C_SIZE_T),     dimension(*)      :: elements_tags_sur
       integer(C_SIZE_T),     dimension(*)      :: elements_nodes_sur

       integer(C_INT),        value             :: algorithm_3D
       real   (C_DOUBLE),     value             :: max_anisotropy
       real   (C_DOUBLE),     value             :: min_size
       real   (C_DOUBLE),     value             :: max_size
       real   (C_DOUBLE),     dimension(*)      :: mesh_size       

       integer(C_SIZE_T),     intent(out)       :: new_nnodes
       type(C_PTR),           intent(out)       :: new_nodes_tags
       type(C_PTR),           intent(out)       :: new_nodes_coordinates
       integer(C_SIZE_T),     intent(out)       :: new_ntypes
       type(C_PTR),           intent(out)       :: new_types_tags
       type(C_PTR),           intent(out)       :: new_nelements_bytype
       type(C_PTR),           intent(out)       :: new_elements_tags
       type(C_PTR),           intent(out)       :: new_elements_nodes       
       integer(C_INT),        value             :: ipass
     END SUBROUTINE reMeshing

     SUBROUTINE freeArrays(&
          &    nodes_tags, nodes_coordinates, types_tags, &
          &    nelements_bytype, elements_tags, elements_nodes ) BIND (c,NAME='freeArrays')

       use, intrinsic ::  ISO_C_BINDING
       implicit none       
       type(C_PTR),           intent(out)       :: nodes_tags
       type(C_PTR),           intent(out)       :: nodes_coordinates
       type(C_PTR),           intent(out)       :: types_tags
       type(C_PTR),           intent(out)       :: nelements_bytype
       type(C_PTR),           intent(out)       :: elements_tags
       type(C_PTR),           intent(out)       :: elements_nodes

     END SUBROUTINE freeArrays


  END INTERFACE

  public :: alya2gmsh_initialization
  public :: alya2gmsh_remeshing
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-30
  !> @brief   Allocate a gmsh mesh
  !> @details DAllocate a gmsh mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine alloca(mesh_gmsh)

    class(gmsh_type), intent(inout) :: mesh_gmsh

    allocate(mesh_gmsh % nodes_tags8        (mesh_gmsh % nnodes8)      )
    allocate(mesh_gmsh % nodes_coordinates8 (mesh_gmsh % nnodes8*3)    )
    allocate(mesh_gmsh % types_tags8        (mesh_gmsh % ntypes8)      )
    allocate(mesh_gmsh % nelements_bytype8  (mesh_gmsh % ntypes8)      )
    allocate(mesh_gmsh % elements_tags8     (mesh_gmsh % nelem8)       )
    allocate(mesh_gmsh % elements_nodes8    (mesh_gmsh % nelem_nodes8) )

  end subroutine alloca

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-30
  !> @brief   Deallocate a gmsh mesh
  !> @details Deallocate a gmsh mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo(mesh_gmsh)

    class(gmsh_type), intent(inout) :: mesh_gmsh

    if( associated(mesh_gmsh % nodes_tags8)         ) deallocate(mesh_gmsh % nodes_tags8)     
    if( associated(mesh_gmsh % nodes_coordinates8)  ) deallocate(mesh_gmsh % nodes_coordinates8)
    !if( mesh_gmsh % ntypes8 > 1) then
       if( associated(mesh_gmsh % types_tags8)      ) deallocate(mesh_gmsh % types_tags8)
       if( associated(mesh_gmsh % nelements_bytype8)) deallocate(mesh_gmsh % nelements_bytype8)
       if( associated(mesh_gmsh % elements_tags8)      ) deallocate(mesh_gmsh % elements_tags8)
    !end if
    if( associated(mesh_gmsh % elements_nodes8)     ) deallocate(mesh_gmsh % elements_nodes8)  

  end subroutine deallo

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-30
  !> @brief   Deallocate a gmsh mesh in c format
  !> @details Deallocate a gmsh mesh in c format
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo_c(mesh_gmsh_c)

    class(gmsh_type_c), intent(inout) :: mesh_gmsh_c

#ifdef ALYA_GMSH
    call freeArrays(& 
         mesh_gmsh_c % nodes_tags8, &
         mesh_gmsh_c % nodes_coordinates8, &
         mesh_gmsh_c % types_tags8,        &           
         mesh_gmsh_c % nelements_bytype8,  &
         mesh_gmsh_c % elements_tags8,     &
         mesh_gmsh_c % elements_nodes8)
#endif

  end subroutine deallo_c

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-30
  !> @brief   Initialize a gmsh mesh
  !> @details Initialize a gmsh mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine init(mesh_gmsh)

    class(gmsh_type), intent(inout) :: mesh_gmsh

    mesh_gmsh % nnodes8      = 0_8
    mesh_gmsh % ntypes8      = 0_8
    mesh_gmsh % nelem8       = 0_8   
    mesh_gmsh % nelem_nodes8 = 0_8   

    nullify(mesh_gmsh % nodes_tags8)        
    nullify(mesh_gmsh % nodes_coordinates8)
    nullify(mesh_gmsh % types_tags8)        
    nullify(mesh_gmsh % nelements_bytype8)  
    nullify(mesh_gmsh % elements_tags8)     
    nullify(mesh_gmsh % elements_nodes8)    

  end subroutine init

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-30
  !> @brief   Initialization
  !> @details Initialization
  !>
  !>          1 is actually "Line 2" 
  !>          2 is actually "Triangle 3" 
  !>          3 is actually "Quadrilateral 4" 
  !>          4 is actually "Tetrahedron 4" 
  !>          5 is actually "Hexahedron 8" 
  !>          6 is actually "Prism 6" 
  !>          7 is actually "Pyramid 5" 
  !>          8 is actually "Line 3" 
  !>          9 is actually "Triangle 6" 
  !>         10 is actually "Quadrilateral 9" 
  !> 
  !-----------------------------------------------------------------------

  subroutine alya2gmsh_initialization()

    integer(ip) :: ii,jj

    memor            =  0_8

    permr_a2g        =  0_ip
    permr_a2g(BAR02) =  1_ip 
    permr_a2g(BAR03) =  8_ip
    permr_a2g(TRI03) =  2_ip 
    permr_a2g(TRI06) =  9_ip 
    permr_a2g(QUA04) =  3_ip  
    permr_a2g(QUA09) = 10_ip 
    permr_a2g(TET04) =  4_ip 
    permr_a2g(HEX08) =  5_ip 
    permr_a2g(PYR05) =  7_ip
    permr_a2g(PEN06) =  6_ip

    do ii = 1,size(permr_a2g)
       jj = permr_a2g(ii)
       if( jj /= 0 ) permr_g2a(jj) = ii
    end do

  end subroutine alya2gmsh_initialization

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-30
  !> @brief   Convert Alya type to gmsh type
  !> @details Convert Alya type to gmsh type
  !> 
  !-----------------------------------------------------------------------

  subroutine alya2gmsh_alya_to_gmsh(mesh_alya,mesh_gmsh)

    type(mesh_type_basic), intent(in)    :: mesh_alya
    type(gmsh_type),       intent(inout) :: mesh_gmsh
    integer(ip)                          :: ntypes
    integer(ip),           pointer       :: types_tags(:)
    integer(ip),           pointer       :: nelements_bytype(:)
    integer(ip)                          :: idime,ipoin,ielem,itype
    integer(ip)                          :: kelem,kpoin,inode,pnode
    integer(ip)                          :: pelty
    !
    ! List of element types
    !
    nullify(types_tags)
    nullify(nelements_bytype)
    call mesh_alya % types(ntypes,types_tags,nelements_bytype)
    !
    ! Dimensions
    !
    mesh_gmsh % nnodes8      = int  (mesh_alya % npoin,      KIND=8)
    mesh_gmsh % nelem8       = int  (mesh_alya % nelem,      KIND=8)
    mesh_gmsh % ntypes8      = int  (ntypes,                 KIND=8)
    mesh_gmsh % nelem_nodes8 = count(mesh_alya % lnods/=0)
    !
    ! Allocate
    !
    call mesh_gmsh % alloca()    
    !
    ! Convert node arrays
    !
    kpoin = 0
    do ipoin = 1,mesh_alya % npoin
       mesh_gmsh % nodes_tags8(ipoin) = int(ipoin,KIND=8) !int(mesh_alya % lninv_loc(ipoin),KIND=8)
       do idime = 1,mesh_alya % ndime
          kpoin = kpoin + 1
          mesh_gmsh % nodes_coordinates8(kpoin) = mesh_alya % coord(idime,ipoin)
       end do
       if( mesh_alya % ndime == 2 ) then
          kpoin = kpoin + 1
          mesh_gmsh % nodes_coordinates8(kpoin) = 0.0_rp
       end if
    end do
    !
    ! Convert element arrays
    !
    kelem = 0
    do ielem = 1,mesh_alya % nelem
       mesh_gmsh % elements_tags8(ielem) = ielem !int(mesh_alya % leinv_loc(ielem),KIND=8)
       pelty = mesh_alya % ltype(ielem)
       pnode = element_type(pelty) % number_nodes
       do inode = 1,pnode
          kelem = kelem + 1
          mesh_gmsh % elements_nodes8(kelem) = int(mesh_alya % lnods(inode,ielem),KIND=8)
       end do
    end do
    !
    ! Types
    !    
    do itype = 1,ntypes
       mesh_gmsh % nelements_bytype8(itype) = int(nelements_bytype(itype),KIND=8)
       mesh_gmsh % types_tags8(itype)       = int( permr_a2g(types_tags(itype)),KIND=8)
    end do
    deallocate(types_tags)
    deallocate(nelements_bytype)

  end subroutine alya2gmsh_alya_to_gmsh

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-30
  !> @brief   Convert Alya type to gmsh type
  !> @details Convert Alya type to gmsh type
  !> 
  !-----------------------------------------------------------------------

  subroutine alya2gmsh_gmsh_to_alya(mesh_gmsh,mesh_alya)

    type(gmsh_type),       intent(in)    :: mesh_gmsh
    type(mesh_type_basic), intent(inout) :: mesh_alya
    integer(ip)                          :: idime,ipoin,ielem,itype
    integer(ip)                          :: kelem,kpoin,inode,pnode
    integer(ip)                          :: nelem,pelty
    !
    ! Initialize alya mesh
    !
    
    do ipoin = 1,mesh_alya % npoin
       mesh_alya % lninv_loc(ipoin) = ipoin !int(mesh_gmsh % nodes_tags8(ipoin),KIND = ip)
       do idime = 1,mesh_alya % ndime
          kpoin = kpoin + 1
          mesh_alya % coord(idime,ipoin) = mesh_gmsh % nodes_coordinates8(kpoin)
       end do
    end do
    call mesh_alya % init()
    !
    ! Dimensions
    !
    mesh_alya % ndime = ndime
    mesh_alya % npoin = int(mesh_gmsh % nnodes8,KIND=ip)
    mesh_alya % nelem = int(mesh_gmsh % nelem8, KIND=ip)
    !mesh_alya % mnode = maxval(nnode(permr_g2a( int(mesh_gmsh % types_tags8(:), KIND =ip) )))
    mesh_alya % mnode = 0
    do itype = 1, int(mesh_gmsh % ntypes8,KIND =ip)
       pelty             = permr_g2a(int(mesh_gmsh % types_tags8(itype), KIND = ip))
       mesh_alya % mnode = max(mesh_alya % mnode, element_type(pelty) % number_nodes)       
    end do    
    !
    ! Allocate
    !
    call mesh_alya % alloca()
    !
    ! Convert node arrays
    !
    kpoin = 0
    do ipoin = 1,mesh_alya % npoin
       mesh_alya % lninv_loc(ipoin) = ipoin !int(mesh_gmsh % nodes_tags8(ipoin),KIND=ip)
       do idime = 1,mesh_alya % ndime
          kpoin = kpoin + 1
          mesh_alya % coord(idime,ipoin) = real(mesh_gmsh % nodes_coordinates8(kpoin),KIND = rp)
       end do
       if( mesh_alya % ndime == 2 ) kpoin = kpoin + 1
    end do

    
    !
    ! Convert element arrays
    !
    kelem = 1    
    do itype = 1,int(mesh_gmsh % ntypes8, KIND =ip)
       nelem = int(mesh_gmsh % nelements_bytype8(itype), KIND =ip)
       mesh_alya % ltype(kelem:kelem+nelem-1) = permr_g2a(int(mesh_gmsh % types_tags8(itype), KIND = ip))
       kelem = kelem + int(mesh_gmsh % nelements_bytype8(itype), KIND =ip)
    end do
    kelem = 0
    do ielem = 1,mesh_alya % nelem
       pelty = mesh_alya % ltype(ielem)
       pnode = element_type(pelty) % number_nodes
       mesh_alya % leinv_loc(ielem) = ielem !int(mesh_gmsh % elements_tags8(ielem),KIND = ip)
       do inode = 1,pnode
          kelem = kelem + 1
          mesh_alya % lnods(inode,ielem) = int(mesh_gmsh % elements_nodes8(kelem),ip)
       end do
    end do
  end subroutine alya2gmsh_gmsh_to_alya

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-30
  !> @brief   Remesh
  !> @details Remesh
  !> 
  !-----------------------------------------------------------------------

  subroutine alya2gmsh_remeshing(mesh_alya_new, mesh_alya_vol, mesh_alya_sur, mesh_size, min_size, max_size)

    class(mesh_type_basic), intent(inout) :: mesh_alya_new
    class(mesh_type_basic), intent(in)    :: mesh_alya_vol
    class(mesh_type_basic), intent(in)    :: mesh_alya_sur    
    real(rp),      pointer, intent(in)    :: mesh_size(:)
    real(rp),               intent(in)    :: max_size, min_size
    !integer(ip),            intent(inout) :: ipass
    integer(4)                            :: algorithm_3D
    integer(4)                            :: ndime4
    real(rp)                              :: max_anisotropy      
    type(gmsh_type)                       :: mesh_gmsh_vol    
    type(gmsh_type)                       :: mesh_gmsh_sur
    type(gmsh_type_c)                     :: mesh_gmsh_c
    integer(4), save                      :: ipass = 0_4
    
    ipass = ipass + 1_4

    if(mesh_alya_vol % nelem == 0 .or. mesh_alya_sur % nelem == 0) then
       !
       ! Empty mesh
       !
       call mesh_alya_new % init()

    else
       !
       ! Convert Alya mesh to gmsh mesh
       !
       call mesh_gmsh_vol % init()
       call alya2gmsh_alya_to_gmsh(mesh_alya_vol,mesh_gmsh_vol)
       call mesh_gmsh_sur % init()       
       call alya2gmsh_alya_to_gmsh(mesh_alya_sur,mesh_gmsh_sur)       
       !
       ! Remesh
       !
       algorithm_3D   = 1_4
       max_anisotropy = 1.0e+10_rp
       ndime4         = int(ndime,4)

       !do ipass = 1, 100000
       !   print *,""
       !   print *,"loop: ", ipass
          call reMeshing(&
               ndime4,                           &
               
               mesh_gmsh_vol   % nnodes8,            &
               mesh_gmsh_vol   % nodes_tags8,        &
               mesh_gmsh_vol   % nodes_coordinates8, &
               mesh_gmsh_vol   % ntypes8,            &
               mesh_gmsh_vol   % types_tags8,        &
               mesh_gmsh_vol   % nelements_bytype8,  &
               mesh_gmsh_vol   % elements_tags8,     &
               mesh_gmsh_vol   % elements_nodes8,    &
               
               mesh_gmsh_sur   % nnodes8,            &
               mesh_gmsh_sur   % nodes_tags8,        &
               mesh_gmsh_sur   % nodes_coordinates8, &
               mesh_gmsh_sur   % ntypes8,            &
               mesh_gmsh_sur   % types_tags8,        &
               mesh_gmsh_sur   % nelements_bytype8,  &
               mesh_gmsh_sur   % elements_tags8,     &
               mesh_gmsh_sur   % elements_nodes8,    &
               
               
               algorithm_3D,                     &
               real(max_anisotropy,  KIND=8),    &                        
               real(min_size,  KIND=8),          &
               real(max_size,  KIND=8),          &
               real(mesh_size, KIND=8),          &
               
               mesh_gmsh_c % nnodes8,            &
               mesh_gmsh_c % nodes_tags8,        &
               mesh_gmsh_c % nodes_coordinates8, &
               mesh_gmsh_c % ntypes8,            &
               mesh_gmsh_c % types_tags8,        &
               mesh_gmsh_c % nelements_bytype8,  &
               mesh_gmsh_c % elements_tags8,     &
               mesh_gmsh_c % elements_nodes8,    &
               ipass)
          
          !call mesh_gmsh_c % deallo_c()
       !end do
                 
       call mesh_gmsh_sur % deallo()
       call mesh_gmsh_vol % deallo()
       !call mesh_gmsh_vol % init()       
       !
       ! Convert gmhs_c to gmsh
       !   
       call alya2gmsh_gmsh_c_to_gmsh(mesh_gmsh_c,mesh_gmsh_vol)
       !
       ! Convert gmsh to Alya 
       !
       call mesh_alya_new % init()

       call alya2gmsh_gmsh_to_alya(mesh_gmsh_vol,mesh_alya_new)       
       !
       ! Deallocate gmsh_c
       !
       call mesh_gmsh_c % deallo_c()
       !call mesh_gmsh_vol % deallo()
    end if
    
  end subroutine alya2gmsh_remeshing

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-30
  !> @brief   Remesh
  !> @details Remesh
  !> 
  !-----------------------------------------------------------------------

  subroutine alya2gmsh_gmsh_c_to_gmsh(mesh_gmsh_c,mesh_gmsh)

    type(gmsh_type_c), intent(inout) :: mesh_gmsh_c
    type(gmsh_type),   intent(inout) :: mesh_gmsh
    integer(ip)                      :: itype,itype_gmsh,itype_alya
    integer(8)                       :: ncoords8,pnode8
    !
    ! Pass pointers from C to F90
    !
    mesh_gmsh % nnodes8 = mesh_gmsh_c % nnodes8 
    mesh_gmsh % ntypes8 = mesh_gmsh_c % ntypes8
    ncoords8            = mesh_gmsh   % nnodes8*3_8

    call C_F_POINTER(mesh_gmsh_c % nodes_tags8,       mesh_gmsh % nodes_tags8,        [mesh_gmsh % nnodes8])
    call C_F_POINTER(mesh_gmsh_c % nodes_coordinates8,mesh_gmsh % nodes_coordinates8, [ncoords8]           )
    call C_F_POINTER(mesh_gmsh_c % types_tags8,       mesh_gmsh % types_tags8,        [mesh_gmsh % ntypes8])
    call C_F_POINTER(mesh_gmsh_c % nelements_bytype8, mesh_gmsh % nelements_bytype8,  [mesh_gmsh % ntypes8])

    mesh_gmsh % nelem8 = sum(mesh_gmsh % nelements_bytype8)

    call C_F_POINTER(mesh_gmsh_c % elements_tags8,    mesh_gmsh % elements_tags8,     [mesh_gmsh % nelem8])

    mesh_gmsh % nelem_nodes8 = 0
    do itype = 1,int(mesh_gmsh % ntypes8,KIND=ip)
       itype_gmsh               = int(mesh_gmsh % types_tags8(itype),ip)
       itype_alya               = int(permr_g2a(itype_gmsh),KIND=ip)
       pnode8                   = int(element_type(itype_alya) % number_nodes,8)
       mesh_gmsh % nelem_nodes8 = mesh_gmsh % nelem_nodes8 + pnode8 * mesh_gmsh % nelements_bytype8(itype)
    end do
    
    call C_F_POINTER(mesh_gmsh_c % elements_nodes8,mesh_gmsh % elements_nodes8, [mesh_gmsh % nelem_nodes8])

  end subroutine alya2gmsh_gmsh_c_to_gmsh

end module mod_alya2gmsh
!> @}
