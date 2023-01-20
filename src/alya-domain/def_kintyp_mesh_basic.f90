!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    def_mesh_type.f90
!> @author  houzeaux
!> @date    2020-04-24
!> @brief   Mesh type
!> @details Mesh type definitions
!>
!>          Basic type
!>          ----------
!>          This type is intended to contain all element
!>          types. Thus boundaries are not represented specifically.
!>          This type is similar to gmsh structure.
!>
!>          |-> Basic procedures are defined here. More complex ones
!>              are in the associated module: mod_mesh_type_basic.f90
!>          !-> Postprocess is in mod_postpr_mesh.f90
!>
!>          Extended type
!>          -------------
!>          Lot more info about the mesh. For example, boundaries are
!>          explicitly declared.
!>          This is the native mesh type of Alya.
!>
!>
!-----------------------------------------------------------------------

module def_kintyp_mesh_basic

  use def_kintyp_basic,      only : ip,rp,i1p,i1pp,r3p,lg
  use def_kintyp_comm,       only : comm_data_par_basic
  use def_isoparametric,     only : isoparametric
  use def_parame,            only : pi
  use def_maths_bin,         only : maths_bin
  use def_maths_bin,         only : LINKED_LIST_BIN
  use mod_memory_basic,      only : memory_alloca
  use mod_memory_basic,      only : memory_deallo
  use mod_memory_basic,      only : memory_copy
  use mod_memory_basic,      only : memory_size
  use mod_memory_tools,      only : memory_counter_ini
  use mod_memory_tools,      only : memory_counter_end
  use mod_optional_argument, only : optional_argument
  use mod_elmgeo,            only : element_type
  use mod_elmgeo,            only : elmgeo_element_volume
  use mod_elmgeo,            only : elmgeo_element_distance
  use mod_elmgeo,            only : elmgeo_element_name_to_type
  use mod_elmgeo,            only : elmgeo_element_type  
  use def_elmgeo,            only : BAR02_TO_BAR02
  use def_elmgeo,            only : BAR03_TO_BAR02
  use def_elmgeo,            only : BAR03_TO_BAR03 
  use def_elmgeo,            only : BAR04_TO_BAR02
  use def_elmgeo,            only : BAR04_TO_BAR04 
  use def_elmgeo,            only : HEX08_TO_HEX08
  use def_elmgeo,            only : HEX27_TO_HEX08
  use def_elmgeo,            only : TRI03_TO_TRI03
  use def_elmgeo,            only : TRI06_TO_TRI03
  use def_elmgeo,            only : QUA04_TO_QUA04
  use def_elmgeo,            only : QUA09_TO_QUA04
  use def_elmgeo,            only : QUA16_TO_QUA04
  use def_elmgeo,            only : HEX08_TO_TET04
  use def_elmgeo,            only : PEN06_TO_PEN06
  use def_elmgeo,            only : PEN06_TO_TET04
  use def_elmgeo,            only : PYR05_TO_TET04
  use def_elmgeo,            only : PYR05_TO_PYR05
  use def_elmgeo,            only : TET04_TO_TET04
  use def_elmgeo,            only : TET10_TO_TET04
  use def_elmtyp,            only : POINT
  use def_elmtyp,            only : BAR02,BAR03,BAR04
  use def_elmtyp,            only : TRI03,TRI06
  use def_elmtyp,            only : QUA04,QUA09,QUA16
  use def_elmtyp,            only : HEX08,HEX27
  use def_elmtyp,            only : PEN06
  use def_elmtyp,            only : PYR05
  use def_elmtyp,            only : TET04,TET10
  use def_elmtyp,            only : BAR3D
  use def_elmtyp,            only : POI3D 
  use def_elmtyp,            only : element_num_ini
  use def_elmtyp,            only : element_num_end
  use def_elmtyp,            only : element_max 
  use mod_maths_sort,        only : maths_heap_sort
  use mod_maths_basic,       only : maths_normalize_vector
  use mod_maths_basic,       only : maths_rotation_matrix_3D
  use mod_maths_geometry,    only : maths_local_orthonormal_basis
  use mod_htable,            only : hash_t
  use mod_htable,            only : htaini
  use mod_htable,            only : htalid
  use mod_htable,            only : htaadd
  use mod_htable,            only : htades
  use mod_htable,            only : htable_initialization  
  use mod_ecoute,            only : ecoute
  use mod_ecoute,            only : ecoute_set_read_unit
  use mod_ecoute,            only : ecoute_reach_section
  use mod_strings,           only : integer_to_string
  use def_quadrature,        only : quadrature
  use def_quadrature,        only : GAUSS_LEGENDRE_RULE  
  use def_quadrature,        only : CLOSED_RULE 
  use def_quadrature,        only : TRAPEZOIDAL_RULE
  use def_quadrature,        only : CHEBYSHEV_RULE
  use def_inpout,            only : getint 
  use def_inpout,            only : param
  use def_inpout,            only : nnpar
  use def_inpout,            only : words 
  use def_inpout,            only : exists
  use def_inpout,            only : getcha
  use def_isoparametric,     only : LAGRANGE_INTERPOLATION 
  use def_isoparametric,     only : CHEBYSHEV_INTERPOLATION
  use mod_std
  implicit none

  private

  !----------------------------------------------------------------------
  !
  ! Basic mesh type
  !
  !----------------------------------------------------------------------

  character(20), parameter :: mesh_name_default = 'MESH'
  real(rp),      parameter :: epsil = epsilon(1.0_rp)
  integer(ip),   parameter :: NODE_TAG     = 1
  integer(ip),   parameter :: ELEMENT_TAG  = 2
  integer(ip),   parameter :: BOUNDARY_TAG = 3
  
  type tag
     character(len=:),      allocatable :: name
     integer(ip)                        :: type
     real(rp),                  pointer :: values(:)
     class(mesh_type_basic),    pointer :: mesh                   ! Parent mesh 
   contains
     procedure,                 pass    :: init   => init_tag
     procedure,                 pass    :: copy   => copy_tag
  end type tag
  
  type mesh_type_basic
     character(LEN=20)                  :: name                   ! Mesh name
     integer(ip)                        :: id                     ! An identifier
     integer(ip)                        :: ndime                  ! Space dimension
     integer(ip)                        :: mnode                  ! Max number of nodes per element
     integer(ip)                        :: nelem                  ! # elements
     integer(ip)                        :: npoin                  ! # nodes
     integer(ip)                        :: ntags                  ! # tags
     integer(8)                         :: memor(2)               ! Memory counter
     integer(ip),               pointer :: lnods(:,:)             ! Element connectivity
     integer(ip),               pointer :: ltype(:)               ! Element type
     integer(ip),               pointer :: leinv_loc(:)           ! Element type
     integer(ip),               pointer :: lninv_loc(:)           ! Node numbering
     real(rp),                  pointer :: coord(:,:)             ! Node coordinates
     type(tag),                 pointer :: tags(:)                ! Tags with types and values
  
     integer(ip),               pointer :: permn(:)               ! Permutation from a parent mesh
     integer(ip),               pointer :: perme(:)               ! Permutation from a parent mesh
     integer(ip),               pointer :: permn_grandpa(:)       ! Permutation from a grandparent mesh
     integer(ip),               pointer :: perme_grandpa(:)       ! Permutation from a grandparent mesh
     class(bmsh_type_basic),    pointer :: boundary               ! Associated boundary mesh
     class(mesh_type_basic),    pointer :: parent                 ! Parent mesh 
     type(comm_data_par_basic)          :: comm                   ! Communication arrays
     type(isoparametric)                :: iso(element_max)       ! Iso-parametric arrays
     type(quadrature)                   :: quad(element_max)      ! Quadrature rule parameters
   contains
     procedure,                 pass    :: init => init_basic     ! Initialize all
     procedure,                 pass    :: alloca                 ! Allocate 
     procedure,                 pass    :: alloca_tag             ! Allocate one single tag
     procedure,                 pass    :: alloca_tags            ! Allocate tags
     procedure,                 pass    :: alloca_com             ! Allocate communication
     procedure,                 pass    :: alloca_quad            ! Allocate quadrature
     procedure,                 pass    :: alloca_iso             ! Allocate iso-parametric
     procedure,                 pass    :: deallo                 ! Deallocate
     procedure,                 pass    :: deallo_tags            ! Deallocate tags
     procedure,                 pass    :: deallo_quad            ! Deallocate quadrature
     procedure,                 pass    :: deallo_iso             ! Deallocate iso-parametric
     procedure,                 pass    :: merge                  ! Merge two meshes
     procedure,                 pass    :: append                 ! Append two meshes
     procedure,                 pass    :: collapse               ! Collapse nodes of a mesh using coordinates
     procedure,                 pass    :: identical_coord        ! Find a list of nodes with identical coordinates
     procedure,                 pass    :: dims                   ! Get dimensions from pointer sizes
     procedure,                 pass    :: types                  ! Number of types and list of types
     procedure,                 pass    :: boundary_mesh          ! Create a boundary mesh
     procedure,                 pass    :: boundary_nodes         ! Identify boundary nodes
     procedure,                 pass    :: output                 ! Output a mesh
     procedure,                 pass    :: extract                ! Extract a submesh using a mask
     procedure,                 pass    :: extract_boundary       ! Extract a boundary mesh
     procedure,                 pass    :: cartesian_mesh         ! Generate a Cartesian mesh
     procedure,                 pass    :: identity               ! Put identity numbering
     procedure,                 pass    :: point_to               ! Point to a mesh
     procedure,                 pass    :: print_info             ! Print info about a mesh
     procedure,                 pass    :: copy                   ! Copy a mesh (same as equal but with optional args)
     procedure,                 pass    :: copy_comm              ! Copy communication arrays
     procedure,                 pass    :: cut_level              ! Cut a mesh given at a level field equal to zero
     procedure,                 pass    :: element_bb             ! Element bounding boxes
     procedure,                 pass    :: boundary_bb            ! Boundary bounding boxes
     procedure,                 pass    :: mesh_bb                ! Create a mesh from bounding boxes
     procedure,                 pass    :: mesh_from_BAR3D        ! Create a 3D mesh from a BAR3D mesh
     procedure,                 pass    :: assoc                  ! Associate a boundary mesh
     procedure,                 pass    :: disassoc               ! Disassociate a boundary mesh
     procedure,                 pass    :: renumber               ! Reneumber a mesh
     procedure,                 pass    :: null_arrays            ! Nullify mesh arrays
     procedure,                 pass    :: volume                 ! Compute the volume of a mesh
     procedure,                 pass    :: mesh_3Dcircle          ! Create a 3D circle mesh with QUA04 elements
     procedure,                 pass    :: mesh_3Dring            ! Create a 3D ring mesh with POI3D elements
     procedure,                 pass    :: mesh_3Dbox             ! Create a 3D ring mesh with POI3D elements
     procedure,                 pass    :: check_lnods            ! Check if lnods is correct
     procedure,                 pass    :: check                  ! Check mesh
     procedure,                 pass    :: centroid               ! Return the centroid coordinates
     procedure,                 pass    :: equal                  ! Copy a mesh
     procedure,                 pass    :: grandpa_permutation    ! Get grandparent permutation indices 
     procedure,                 pass    :: read_from_file         ! Read a mesh
     procedure,                 pass    :: linearize              ! Linearize a mesh
     procedure,                 pass    :: tet04_mesh             ! Transform a generic mesh into a TET04 mesh
     generic,   public                  :: assignment(=) => equal
     procedure,                 pass    :: results_r1             ! Output single results 
     procedure,                 pass    :: results_r2             ! Output single results 
     procedure,                 pass    :: results_i1             ! Output multiple results
     generic                            :: results =>  &
          &                                results_r1, &
          &                                results_r2, &
          &                                results_i1
   end type mesh_type_basic
  
  !----------------------------------------------------------------------
  !
  ! Basic boundary mesh type
  !
  !----------------------------------------------------------------------
  
  type, extends(mesh_type_basic) :: bmsh_type_basic
     integer(ip),               pointer :: lelbo(:)                  ! Element connect to boundary
     integer(ip),               pointer :: lboel(:,:)                ! Local node correspondance
     class(mesh_type_basic),    pointer :: mesh                      ! Corresponding volume mesh
   contains
     procedure,                 pass    :: init     => init_bmsh     ! Initialize all
     procedure,                 pass    :: centroid => centroid_bmsh ! Return the centroid coordinates
  end type bmsh_type_basic
    
  public :: mesh_type_basic
  public :: bmsh_type_basic
  public :: mesh_name_default
  public :: NODE_TAG
  public :: ELEMENT_TAG
  public :: BOUNDARY_TAG

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Initialization
  !> @details Initialization of a tag
  !> 
  !-----------------------------------------------------------------------

  subroutine init_tag(self,wname)
    class(tag)                             :: self
    character(len=*), optional, intent(in) :: wname

    self % name = optional_argument('TAG',wname)
    self % type = 0_ip
    nullify(self % values)
    
  end subroutine init_tag
  
  subroutine copy_tag(new_self,cpy_self,MEMORY_COUNTER,DO_ALLOCATE)
    
    class(tag),                      intent(inout) :: new_self
    class(tag),                      intent(in)    :: cpy_self
    integer(8),   optional,          intent(inout) :: MEMORY_COUNTER(2)
    logical(lg),  optional,          intent(in)    :: DO_ALLOCATE
    integer(ip)                                    :: n
    integer(8)                                     :: memor_loc(2)
    logical(lg)                                    :: if_allocate
    
    memor_loc   = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
    if_allocate = optional_argument(.true.,DO_ALLOCATE)
    
    n = memory_size(cpy_self % values,1_ip)
    if( n > 0 ) then
       if( if_allocate ) then
          call memory_alloca(memor_loc,trim(new_self % name)//' % FIELD % VALUES','copy_tag',new_self % values,n)
       end if
       n = min(n,memory_size(new_self % values,1_ip))
       new_self % values(1:n) = cpy_self % values(1:n)
    end if
    
    new_self % name = cpy_self % name
    new_self % type = cpy_self % type

    if( present(MEMORY_COUNTER) )  MEMORY_COUNTER = memor_loc
   
  end subroutine copy_tag
    
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Allocate tags
  !> @details Allocate tags
  !> 
  !-----------------------------------------------------------------------

  subroutine alloca_tags(mesh)
    class(mesh_type_basic)  :: mesh
    integer(ip)             :: itag
    
    allocate(mesh % tags(mesh % ntags))
    do itag = 1,mesh % ntags
       call mesh % tags(itag) % init()
    end do
    
  end subroutine alloca_tags
  
  subroutine alloca_tag(mesh,TAG_ID,MEMORY_COUNTER)
    
    class(mesh_type_basic)                 :: mesh
    integer(ip),  optional,  intent(in)    :: TAG_ID
    integer(8),   optional,  intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                             :: memor_loc(2)
    integer(ip)                            :: tag_ini,tag_end,itag
    
    call memory_counter_ini(memor_loc,mesh % memor,MEMORY_COUNTER)

    if( present(tag_id) ) then
       tag_ini = TAG_ID
       tag_end = TAG_ID
    else
       tag_ini = 1
       tag_end = mesh % ntags
    end if
    
    if( associated(mesh % tags) ) then
       do itag = tag_ini,tag_end
          if( itag > 0 .and. itag <= size(mesh % tags) ) then
             if(      mesh % tags(itag) % type == ELEMENT_TAG  ) then
                call memory_alloca(memor_loc,trim(mesh % name)//' % VALUES','alloca_tag',mesh % tags(itag) % values,mesh % nelem)
             else if( mesh % tags(itag) % type == NODE_TAG     ) then
                call memory_alloca(memor_loc,trim(mesh % name)//' % VALUES','alloca_tag',mesh % tags(itag) % values,mesh % npoin)
             else if( mesh % tags(itag) % type == BOUNDARY_TAG ) then
                call memory_alloca(memor_loc,trim(mesh % name)//' % VALUES','alloca_tag',mesh % tags(itag) % values,mesh % boundary % nelem)
             else
                call runend('ALLOCA_TAG: UNKNOWN TYPE OF TAG')
             end if
          end if
       end do
    end if
    
    call memory_counter_end(memor_loc,mesh % memor,MEMORY_COUNTER)

  end subroutine alloca_tag
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Allocate tags
  !> @details Allocate tags
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo_tags(mesh,MEMORY_COUNTER)
    
    class(mesh_type_basic)                 :: mesh
    integer(8),   optional,  intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                             :: memor_loc(2)
    integer(ip)                            :: itag
     
    call memory_counter_ini(memor_loc,mesh % memor,MEMORY_COUNTER)

    if( associated(mesh % tags) ) then
       do itag = 1,size(mesh % tags,KIND=ip)
          call memory_deallo(memor_loc,trim(mesh % name)//' % VALUES','alloca_tag',mesh % tags(itag) % values)
          if( allocated(mesh % tags(itag) % name) ) deallocate(mesh % tags(itag) % name)
       end do
       deallocate(mesh % tags)
    end if
   
    call memory_counter_end(memor_loc,mesh % memor,MEMORY_COUNTER)

  end subroutine deallo_tags

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Allocate quadrature
  !> @details Allocate quadrature
  !> 
  !-----------------------------------------------------------------------

  subroutine alloca_quad(mesh)
    class(mesh_type_basic)  :: mesh
    
    !allocate(mesh % quad(element_max))
    !do ii = 1,element_max
    !   call mesh % quad(ii) % alloca()
    !end do
    
  end subroutine alloca_quad

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Deallocate quadrature
  !> @details Deallocate quadrature
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo_quad(mesh)
    class(mesh_type_basic)  :: mesh
    integer(ip)             :: ii
    
    do ii = 1,element_max
       call mesh % quad(ii) % deallo()
    end do
    
  end subroutine deallo_quad

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Allocate iso
  !> @details Allocate iso
  !> 
  !-----------------------------------------------------------------------

  subroutine alloca_iso(mesh)
    class(mesh_type_basic)  :: mesh
    
    !do ii = 1,element_max
    !   call mesh % iso(ii) % alloca()
    !end do
    
  end subroutine alloca_iso

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Deallocate iso
  !> @details Deallocate iso
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo_iso(mesh)
    class(mesh_type_basic)  :: mesh
    integer(ip)             :: ii
    
    do ii = 1,element_max
       call mesh % iso(ii) % deallo()
    end do
    
  end subroutine deallo_iso

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Initialization
  !> @details Initialization. Should be called prior to any operations
  !>          on a mesh.
  !> 
  !-----------------------------------------------------------------------

  subroutine init_basic(mesh,wname)

    class(mesh_type_basic),          intent(inout) :: mesh
    character(len=*),      optional, intent(in)    :: wname
    integer(ip)                                    :: ii

    if( present(wname) ) then
       mesh % name = trim(wname)
    else
       mesh % name = mesh_name_default
    end if

    mesh % id     = 0_ip
    mesh % ndime  = 0_ip
    mesh % mnode  = 0_ip
    mesh % nelem  = 0_ip
    mesh % npoin  = 0_ip
    mesh % ntags  = 0_ip
    mesh % memor  = 0_8

    call mesh % null_arrays()

    do ii = 1,size(mesh % quad)
       call mesh % quad(ii) % init()
    end do
    do ii = 1,size(mesh % iso)
       call mesh % iso(ii) % init()
    end do   

    call mesh % comm   % init()

  end subroutine init_basic

  subroutine init_bmsh(mesh,wname)

    class(bmsh_type_basic),          intent(inout) :: mesh
    character(len=*),      optional, intent(in)    :: wname
    integer(ip)                                    :: ii

    if( present(wname) ) then
       mesh % name = trim(wname)
    else
       mesh % name = mesh_name_default
    end if

    mesh % id     = 0_ip
    mesh % ndime  = 0_ip
    mesh % mnode  = 0_ip
    mesh % nelem  = 0_ip
    mesh % npoin  = 0_ip
    mesh % ntags  = 0_ip
    mesh % memor  = 0_8

    call mesh % null_arrays()

    do ii = 1,size(mesh % quad)
       call mesh % quad(ii) % init()
    end do
    do ii = 1,size(mesh % iso)
       call mesh % iso(ii) % init()
    end do   

    call mesh % comm   % init()

  end subroutine init_bmsh

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Allocate 
  !> @details Allocate a mesh. Dimensions can be given as arguments.
  !>          If not, it is assumed they were previously defined.
  !>          If one of the array is presnet (LNODS, LTYPE, ETC.) only
  !>          those variables declared are allocated.
  !> 
  !-----------------------------------------------------------------------

  subroutine alloca(&
       mesh,ndime,mnode,nelem,npoin,ntags,&
       MEMORY_COUNTER,LNODS,LTYPE,PERME,LEINV_LOC,LNINV_LOC,COORD,&
       LELBO,LBOEL,PERMN,TAGS,ALL,MESH_NAME)

    class(mesh_type_basic),           intent(inout) :: mesh
    integer(ip),           optional,  intent(in)    :: ndime
    integer(ip),           optional,  intent(in)    :: mnode
    integer(ip),           optional,  intent(in)    :: nelem
    integer(ip),           optional,  intent(in)    :: npoin
    integer(ip),           optional,  intent(in)    :: ntags
    integer(8),            optional,  intent(inout) :: MEMORY_COUNTER(2)
    logical(lg),           optional,  intent(in)    :: LNODS
    logical(lg),           optional,  intent(in)    :: LTYPE
    logical(lg),           optional,  intent(in)    :: PERME
    logical(lg),           optional,  intent(in)    :: LEINV_LOC 
    logical(lg),           optional,  intent(in)    :: LNINV_LOC
    logical(lg),           optional,  intent(in)    :: COORD
    logical(lg),           optional,  intent(in)    :: LELBO
    logical(lg),           optional,  intent(in)    :: LBOEL
    logical(lg),           optional,  intent(in)    :: PERMN    
    logical(lg),           optional,  intent(in)    :: TAGS    
    logical(lg),           optional,  intent(in)    :: ALL    
    character(len=*),      optional,  intent(in)    :: MESH_NAME
    integer(ip)                                     :: ipoin,ielty
    integer(ip)                                     :: ielem
    integer(8)                                      :: memor_loc(2)
    logical(lg)                                     :: IF_LNODS
    logical(lg)                                     :: IF_LTYPE
    logical(lg)                                     :: IF_PERME
    logical(lg)                                     :: IF_LEINV_LOC
    logical(lg)                                     :: IF_LNINV_LOC
    logical(lg)                                     :: IF_COORD
    logical(lg)                                     :: IF_LELBO
    logical(lg)                                     :: IF_LBOEL
    logical(lg)                                     :: IF_PERMN
    logical(lg)                                     :: IF_TAGS
    logical(lg)                                     :: IF_DEFAULT
    character(LEN=:), allocatable                   :: my_mesh_name
    logical(lg)                                     :: if_all
    
    call memory_counter_ini(memor_loc,mesh % memor,MEMORY_COUNTER)

    my_mesh_name = optional_argument(mesh % name,MESH_NAME)

    if(    present(LNODS    ) .or. &
         & present(LTYPE    ) .or. &
         & present(PERME    ) .or. &
         & present(LEINV_LOC) .or. &
         & present(LNINV_LOC) .or. &
         & present(COORD    ) .or. &
         & present(LELBO    ) .or. &
         & present(LBOEL    ) .or. &
         & present(PERMN    ) .or. &
         & present(TAGS     ) ) then
       IF_DEFAULT = optional_argument(.false.,ALL)
       if_all     = .false.
    else
       IF_DEFAULT = .true.
       if_all     = .true.
    end if


    if( .not. if_all ) then
       !
       ! Only allocate specific arrays
       !
       IF_LNODS      = optional_argument(IF_DEFAULT,LNODS    )
       IF_LTYPE      = optional_argument(IF_DEFAULT,LTYPE    )
       IF_PERME      = optional_argument(IF_DEFAULT,PERME    )
       IF_LEINV_LOC  = optional_argument(IF_DEFAULT,LEINV_LOC)
       IF_LNINV_LOC  = optional_argument(IF_DEFAULT,LNINV_LOC)
       IF_COORD      = optional_argument(IF_DEFAULT,COORD    )
       IF_LELBO      = optional_argument(IF_DEFAULT,LELBO    )
       IF_LBOEL      = optional_argument(IF_DEFAULT,LBOEL    )
       IF_PERMN      = optional_argument(IF_DEFAULT,PERMN    )
       IF_TAGS       = optional_argument(IF_DEFAULT,TAGS     )

       if(IF_LNODS    ) call memory_alloca(memor_loc,trim(my_mesh_name)//' % LNODS'    ,'alloca',mesh % lnods    ,mesh % mnode,mesh % nelem)
       if(IF_LTYPE    ) call memory_alloca(memor_loc,trim(my_mesh_name)//' % LTYPE'    ,'alloca',mesh % ltype    ,mesh % nelem)             
       if(IF_PERME    ) call memory_alloca(memor_loc,trim(my_mesh_name)//' % PERME'    ,'alloca',mesh % perme    ,mesh % nelem)
       if(IF_LEINV_LOC) call memory_alloca(memor_loc,trim(my_mesh_name)//' % LEINV_LOC','alloca',mesh % leinv_loc,mesh % nelem)             
       if(IF_LNINV_LOC) call memory_alloca(memor_loc,trim(my_mesh_name)//' % LNINV_LOC','alloca',mesh % lninv_loc,mesh % npoin)             
       if(IF_COORD    ) call memory_alloca(memor_loc,trim(my_mesh_name)//' % COORD'    ,'alloca',mesh % coord    ,mesh % ndime,mesh % npoin)
       if(IF_PERMN    ) call memory_alloca(memor_loc,trim(my_mesh_name)//' % PERMN'    ,'alloca',mesh % permn    ,mesh % npoin)
       if(IF_TAGS     ) call mesh % alloca_tags()

       select type ( mesh )       
       class is ( bmsh_type_basic )
          if(IF_LELBO    ) call memory_alloca(memor_loc,trim(my_mesh_name)//' % LELBO'    ,'alloca',mesh % lelbo    ,mesh % nelem)
          if(IF_LBOEL    ) call memory_alloca(memor_loc,trim(my_mesh_name)//' % LBOEL'    ,'alloca',mesh % lboel    ,mesh % mnode,mesh % nelem)          
       end select

    else
       !
       ! Deallocate and reallocate all the mesh arrays
       !
       if( present(ndime) ) mesh % ndime = ndime
       if( present(mnode) ) mesh % mnode = mnode
       if( present(nelem) ) mesh % nelem = nelem
       if( present(npoin) ) mesh % npoin = npoin
       if( present(ntags) ) mesh % ntags = ntags

       call mesh % deallo(MEMORY_COUNTER=memor_loc)

       call memory_alloca(memor_loc,trim(my_mesh_name)//' % LNODS'    ,'alloca',mesh % lnods    ,mesh % mnode,mesh % nelem)
       call memory_alloca(memor_loc,trim(my_mesh_name)//' % LTYPE'    ,'alloca',mesh % ltype    ,mesh % nelem)             
       call memory_alloca(memor_loc,trim(my_mesh_name)//' % PERME'    ,'alloca',mesh % perme    ,mesh % nelem)
       call memory_alloca(memor_loc,trim(my_mesh_name)//' % LEINV_LOC','alloca',mesh % leinv_loc,mesh % nelem)             
       call memory_alloca(memor_loc,trim(my_mesh_name)//' % LNINV_LOC','alloca',mesh % lninv_loc,mesh % npoin)             
       call memory_alloca(memor_loc,trim(my_mesh_name)//' % COORD'    ,'alloca',mesh % coord    ,mesh % ndime,mesh % npoin)
       call memory_alloca(memor_loc,trim(my_mesh_name)//' % PERMN'    ,'alloca',mesh % permn    ,mesh % npoin)
       call mesh % alloca_tags()
       
       select type ( mesh )       
       class is ( bmsh_type_basic )
          call memory_alloca(memor_loc,trim(my_mesh_name)//' % LELBO'    ,'alloca',mesh % lelbo    ,mesh % nelem)
          call memory_alloca(memor_loc,trim(my_mesh_name)//' % LBOEL'    ,'alloca',mesh % lboel    ,mesh % mnode,mesh % nelem)                    
       end select

       do ielem = 1,mesh % nelem
          mesh % leinv_loc(ielem) = ielem
       end do
       do ipoin = 1,mesh % npoin
          mesh % lninv_loc(ipoin) = ipoin
       end do

       call mesh % comm % init()

    end if
    
    if( allocated(my_mesh_name) ) deallocate(my_mesh_name)
    call memory_counter_end(memor_loc,mesh % memor,MEMORY_COUNTER)

  end subroutine alloca
  
  subroutine alloca_com(mesh,MEMORY_COUNTER,MESH_NAME)

    class(mesh_type_basic),           intent(inout) :: mesh
    integer(8),            optional,  intent(inout) :: MEMORY_COUNTER(2)
    character(len=*),      optional,  intent(in)    :: MESH_NAME
    integer(8)                                      :: memor_loc(2)
    character(LEN=:), allocatable                   :: my_mesh_name
    
    call memory_counter_ini(memor_loc,mesh % memor,MEMORY_COUNTER)
    my_mesh_name = optional_argument(mesh % name,MESH_NAME)

    call memory_alloca(memor_loc,trim(my_mesh_name)//' % COMM % NEIGHTS'   ,'alloca_com',mesh % comm % neights   ,mesh % comm % nneig)
    call memory_alloca(memor_loc,trim(my_mesh_name)//' % COMM % BOUND_SIZE','alloca_com',mesh % comm % bound_size,mesh % comm % nneig+1)
    call memory_alloca(memor_loc,trim(my_mesh_name)//' % COMM % BOUND_PERM','alloca_com',mesh % comm % bound_perm,mesh % comm % bound_dim)

    if( present(MEMORY_COUNTER) )  MEMORY_COUNTER = memor_loc
    if( allocated(my_mesh_name) ) deallocate(my_mesh_name)
    
    call memory_counter_end(memor_loc,mesh % memor,MEMORY_COUNTER)

  end subroutine alloca_com

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Deallocate 
  !> @details Deallocate a mesh. Dimensions are kept just in case they
  !>          would be needed.
  !> 
  !-----------------------------------------------------------------------

  subroutine null_arrays(mesh)

    class(mesh_type_basic), intent(inout) :: mesh
    
    nullify(mesh % lnods     )
    nullify(mesh % ltype     )             
    nullify(mesh % perme     )
    nullify(mesh % leinv_loc )             
    nullify(mesh % lninv_loc )             
    nullify(mesh % coord     )
    nullify(mesh % permn     )
    nullify(mesh % tags      )
    nullify(mesh % permn_grandpa )
    nullify(mesh % perme_grandpa )

    nullify(mesh % boundary )
    nullify(mesh % parent   )

    select type ( mesh )
    class is (bmsh_type_basic)
       nullify(mesh % lelbo)
       nullify(mesh % lboel)
       nullify(mesh % mesh )
    end select

    call mesh % comm % null_arrays() 
    
  end subroutine null_arrays
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Deallocate 
  !> @details Deallocate a mesh. Dimensions are kept just in case they
  !>          would be needed.
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo(mesh,MEMORY_COUNTER,MESH_NAME)

    class(mesh_type_basic)                          :: mesh
    integer(8),            optional,  intent(inout) :: MEMORY_COUNTER(2)
    character(len=*),      optional,  intent(in)    :: MESH_NAME
    character(len=:),      allocatable              :: my_mesh_name
    integer(8)                                      :: memor_loc(2)
    integer(ip)                                     :: ii
    
    call memory_counter_ini(memor_loc,mesh % memor,MEMORY_COUNTER)
    my_mesh_name = optional_argument(mesh % name,MESH_NAME)
    
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LNODS'    ,'deallo',mesh % lnods    )
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LTYPE'    ,'deallo',mesh % ltype    )             
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % PERME'    ,'deallo',mesh % perme    )
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LEINV_LOC','deallo',mesh % leinv_loc)             
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LNINV_LOC','deallo',mesh % lninv_loc)             
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % COORD'    ,'deallo',mesh % coord    )
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % PERMN'    ,'deallo',mesh % permn    )

    nullify(  mesh % boundary )
    nullify(  mesh % parent   )

    call mesh % comm % deallo(memor_loc)
    call mesh % deallo_tags  (memor_loc)
    
    do ii = 1,size(mesh % quad)
       call mesh % quad(ii) % deallo()
    end do
    do ii = 1,size(mesh % iso)
       call mesh % iso(ii) % deallo()
    end do
    
    !if( allocated (mesh % name  ) ) deallocate(mesh % name)
    if( allocated (my_mesh_name ) ) deallocate(my_mesh_name)
 
    call memory_counter_end(memor_loc,mesh % memor,MEMORY_COUNTER)
    
  end subroutine deallo

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux and samaniego
  !> @date    2020-04-24
  !> @brief   Types
  !> @details Return the different element types involved in a mesh.
  !> 
  !-----------------------------------------------------------------------

  subroutine types(mesh,ntype,ltype,num_elements_per_type)

    class(mesh_type_basic), intent(in)    :: mesh
    integer(ip),            intent(out)   :: ntype
    integer(ip), pointer,   intent(inout) :: ltype(:)
    integer(ip), pointer,   intent(inout) :: num_elements_per_type(:)
    integer(ip)                           :: ktype,itype,ielem,iorde
    integer(ip),            allocatable   :: ltype_loc(:)
    integer(ip),            allocatable   :: ltype_sor(:)

    ntype = maxval(mesh % ltype)
    allocate(ltype_loc(ntype))
    allocate(ltype_sor(ntype))
    do itype = 1,ntype
       ltype_loc(itype) = 0
    end do
    
    iorde = 0_ip
    do ielem = 1,mesh % nelem
       itype = mesh % ltype(ielem)
       ! Count the number of elements for a type       
       ltype_loc(itype) = ltype_loc(itype) + 1
       ! Save the order if it is the first time that the type appears            
       if (ltype_loc(itype) == 1) then
          iorde = iorde + 1
          ltype_sor(itype) = iorde                 
       end if
    end do
   
    ntype = count(ltype_loc/=0,KIND=ip)
    allocate(num_elements_per_type(ntype))
    allocate(ltype(ntype))

    ! Add types according the order they appear
    do itype = 1,size(ltype_loc,KIND=ip)
       if( ltype_loc(itype) /= 0 ) then
          ktype = ltype_sor(itype)
          ltype(ktype) = itype
          num_elements_per_type(ktype) = ltype_loc(itype)
       end if
    end do
    deallocate(ltype_loc)
    deallocate(ltype_sor)

  end subroutine types

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Dimensions of a mesh
  !> @details Get the dimensions of a mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine dims(mesh,ndime,mnode,nelem,npoin)

    class(mesh_type_basic), intent(in)  :: mesh
    integer(ip),            intent(out) :: ndime
    integer(ip),            intent(out) :: mnode
    integer(ip),            intent(out) :: nelem
    integer(ip),            intent(out) :: npoin

    ndime = 0
    mnode = 0
    nelem = 0
    npoin = 0

    if(    associated(mesh % ltype)     .and. &
         & associated(mesh % lnods)     .and. &
         & associated(mesh % leinv_loc) ) then
       mnode = size(mesh % lnods,DIM=1_ip,KIND=ip)
       nelem = size(mesh % lnods,DIM=2_ip,KIND=ip)
    end if
    if(    associated(mesh % coord)     .and. &
         & associated(mesh % lninv_loc) ) then
       ndime = size(mesh % coord,DIM=1_ip,KIND=ip)
       npoin = size(mesh % coord,DIM=2_ip,KIND=ip)
    end if

  end subroutine dims

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Copy
  !> @details Copy a mesh NEW_MESH = CPY_MESH
  !> 
  !-----------------------------------------------------------------------
  
  subroutine equal(new_mesh,cpy_mesh)

    class(mesh_type_basic),           intent(in)    :: cpy_mesh
    class(mesh_type_basic),           intent(inout) :: new_mesh

    call new_mesh % copy(cpy_mesh)
    
  end subroutine equal
  
  subroutine copy(new_mesh,cpy_mesh,MEMORY_COUNTER,DO_ALLOCATE,nelem_in,npoin_in)

    class(mesh_type_basic),           intent(in)    :: cpy_mesh
    class(mesh_type_basic),           intent(inout) :: new_mesh
    integer(8),            optional,  intent(inout) :: MEMORY_COUNTER(2)
    logical(lg),           optional,  intent(in)    :: DO_ALLOCATE
    integer(ip),           optional,  intent(in)    :: nelem_in
    integer(ip),           optional,  intent(in)    :: npoin_in
    integer(ip)                                     :: ielem,ipoin
    integer(ip)                                     :: itag
    integer(ip)                                     :: pelem,ppoin
    integer(ip)                                     :: pdime,pnode
    integer(ip)                                     :: nelem_use
    integer(ip)                                     :: npoin_use
    logical(lg)                                     :: if_allocate
    integer(8)                                      :: memor_loc(2)
    
    call memory_counter_ini(memor_loc,new_mesh % memor,MEMORY_COUNTER)
    
    call new_mesh % dims(pdime,pnode,pelem,ppoin)

    if( present(DO_ALLOCATE) ) then
       if_allocate = DO_ALLOCATE
    else
       if_allocate = .true.
       if( ppoin >= cpy_mesh % npoin .and. pelem >= cpy_mesh % nelem ) then
          if( pdime == cpy_mesh % ndime .and. pnode >= cpy_mesh % mnode ) then 
             if_allocate = .false.
          end if
       end if
    end if
    !
    ! Nodes and elements to copy
    !
    if( present(npoin_in) ) then
       npoin_use  = npoin_in
    else
       npoin_use  = cpy_mesh % npoin
    end if
    if( present(nelem_in) ) then          
       nelem_use  = nelem_in
    else
       nelem_use  = cpy_mesh % nelem          
    end if
    !
    ! Allocate if required
    !
    if( if_allocate ) then

       call new_mesh % deallo(MEMORY_COUNTER=memor_loc)
       if(trim(new_mesh % name)==mesh_name_default) new_mesh % name   = cpy_mesh % name
       new_mesh % id     = cpy_mesh % id
       new_mesh % ndime  = cpy_mesh % ndime
       new_mesh % mnode  = cpy_mesh % mnode
       new_mesh % ntags  = cpy_mesh % ntags
       if( present(npoin_in) ) then
          new_mesh % npoin  = npoin_in 
       else
          new_mesh % npoin  = cpy_mesh % npoin
       end if
       if( present(nelem_in) ) then          
          new_mesh % nelem  = nelem_in
       else
          new_mesh % nelem  = cpy_mesh % nelem          
       end if
       
       call new_mesh % alloca(MEMORY_COUNTER=memor_loc)
       
    end if

    if( cpy_mesh % nelem > 0 ) pnode = min(size(cpy_mesh % lnods,DIM=1_ip,KIND=ip),size(new_mesh % lnods,DIM=1_ip,KIND=ip))

    do ielem = 1 , nelem_use 
       new_mesh % lnods(1:pnode,ielem) = cpy_mesh % lnods(1:pnode,ielem) 
       new_mesh % ltype(ielem)         = cpy_mesh % ltype(ielem)
       new_mesh % leinv_loc(ielem)     = cpy_mesh % leinv_loc(ielem)
    end do
    do ipoin = 1 , npoin_use 
       new_mesh % coord(:,ipoin)   = cpy_mesh % coord(:,ipoin) 
       new_mesh % lninv_loc(ipoin) = cpy_mesh % lninv_loc(ipoin)          
    end do
    !
    ! Tags
    !
    if( associated(new_mesh % tags) ) then
       do itag = 1,cpy_mesh % ntags
          call new_mesh % tags(itag) % copy(cpy_mesh % tags(itag),MEMORY_COUNTER=memor_loc,DO_ALLOCATE=if_allocate)
       end do
    end if
    !
    ! Permutations
    !
    if( associated(cpy_mesh % perme) ) then
       do ielem = 1 , nelem_use
          new_mesh % perme(ielem) = cpy_mesh % perme(ielem)           
       end do
    end if
    if( associated(cpy_mesh % permn) ) then
       do ipoin = 1 , npoin_use 
          new_mesh % permn(ipoin) = cpy_mesh % permn(ipoin) 
       end do
    end if

    call memory_counter_end(memor_loc,new_mesh % memor,MEMORY_COUNTER)

  end subroutine copy
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Copy communication
  !> @details Copy the communicaiton arrays 
  !> 
  !-----------------------------------------------------------------------

  subroutine copy_comm(new_mesh,cpy_mesh,MEMORY_COUNTER)
    
    class(mesh_type_basic),           intent(inout) :: cpy_mesh
    class(mesh_type_basic),           intent(inout) :: new_mesh
    integer(8),            optional,  intent(inout) :: MEMORY_COUNTER(2)
    !
    ! Communication
    !
    call new_mesh % comm % copy(cpy_mesh % comm,MEMORY_COUNTER)
    
  end subroutine copy_comm

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Merge
  !> @details Merge a mesh to another. Its like an append operation
  !>          but merging duplicated nodes:
  !>          MESH = MESH  U MESH_GLU
  !>
  !>          MESH_GLU is eventually deallocated 
  !> 
  !-----------------------------------------------------------------------

  subroutine merge(mesh,mesh_glu,MEMORY_COUNTER)

    class(mesh_type_basic),            intent(inout) :: mesh         !< Output mesh
    class(mesh_type_basic),            intent(inout) :: mesh_glu     !< Input mesh (eventually deallocated)
    integer(8),             optional,  intent(inout) :: MEMORY_COUNTER(2)
    type(hash_t)                                     :: htable_3 
    integer(ip)                                      :: lid,ipoin
    integer(ip)                                      :: ielem,kelem
    integer(ip)                                      :: inode,ielty
    integer(ip)                                      :: itag
    logical(lg)                                      :: isin
    type(mesh_type_basic)                            :: mesh_tmp
    integer(8)                                       :: memor_loc(2)

    call memory_counter_ini(memor_loc,mesh % memor,MEMORY_COUNTER)

    call mesh_tmp % init('MESH_TMP')

    if( mesh_glu % nelem == 0 .and. mesh_glu % npoin == 0 ) then

       return

    else if( mesh % nelem == 0 .and. mesh % npoin == 0 ) then
       !
       ! Copy first mesh
       !
       call mesh     % copy  (mesh_glu,MEMORY_COUNTER)
       call mesh_glu % deallo(MEMORY_COUNTER)

    else
       !
       ! Append mesh MESH_GLU to MESH_TMP
       !
       mesh_tmp % ndime = mesh % ndime
       mesh_tmp % mnode = mesh % mnode
       mesh_tmp % nelem = mesh % nelem + mesh_glu % nelem
       mesh_tmp % npoin = mesh % npoin + mesh_glu % npoin
       mesh_tmp % ntags = max(mesh % ntags,mesh_glu % ntags)

       call htable_initialization(htable_3)
       call htades( htable_3)
       call htaini( htable_3, mesh_tmp % npoin, lidson=.true., AUTOMATIC_SIZE=.true.)
       call htaadd( htable_3, mesh % lninv_loc)
       do ipoin = 1,mesh_glu % npoin
          call htaadd(htable_3,mesh_glu % lninv_loc(ipoin),lid,isin)
       end do
       mesh_tmp % npoin = htable_3 % nelem
       call mesh_tmp % alloca(MEMORY_COUNTER=memor_loc)       
       !
       ! Element arrays
       !
       do ielem = 1 , mesh % nelem 
          mesh_tmp % lnods(:,ielem)   = mesh % lnods(:,ielem)
          mesh_tmp % ltype(ielem)     = mesh % ltype(ielem)
          mesh_tmp % perme(ielem)     = mesh % perme(ielem)
          mesh_tmp % leinv_loc(ielem) = mesh % leinv_loc(ielem)
       end do
       kelem = mesh % nelem
       do ielem = 1, mesh_glu % nelem
          kelem                       = kelem + 1
          ielty                       = mesh_glu % ltype(ielem)
          mesh_tmp % ltype(kelem)     = mesh_glu % ltype(ielem)
          mesh_tmp % perme(kelem)     = mesh_glu % perme(ielem)
          mesh_tmp % leinv_loc(kelem) = mesh_glu % leinv_loc(ielem)
          do inode = 1,element_type(ielty) % number_nodes
             ipoin                       = mesh_glu % lnods(inode,ielem)
             if( ipoin > 0 ) then
                lid                         = htalid(htable_3,mesh_glu % lninv_loc(ipoin))          
                mesh_tmp % lnods(inode,kelem) = lid
             end if
          end do
       end do
       !
       ! Node arrays
       !
       do ipoin = 1 , mesh % npoin
          mesh_tmp % lninv_loc(ipoin) = mesh % lninv_loc(ipoin)
          mesh_tmp % coord(:,ipoin)   = mesh % coord(:,ipoin)
          mesh_tmp % permn(ipoin)     = mesh % permn(ipoin)
       end do
       do ipoin = 1,mesh_glu % npoin
          lid                       = htalid(htable_3,mesh_glu % lninv_loc(ipoin))
          mesh_tmp % coord(:,lid)   = mesh_glu % coord(:,ipoin)
          mesh_tmp % lninv_loc(lid) = mesh_glu % lninv_loc(ipoin)
          mesh_tmp % permn(lid)     = mesh_glu % permn(ipoin)
       end do
       !
       ! Tags
       !
       do itag = 1,mesh % ntags
          mesh_tmp % tags(itag) % type = mesh % tags(itag) % type
          call mesh_tmp % alloca_tag(itag,MEMORY_COUNTER=memor_loc)
          if( mesh % tags(itag) % type == ELEMENT_TAG ) then
             do ielem = 1,mesh % nelem
                mesh_tmp % tags(itag) % values(ielem) = mesh % tags(itag) % values(ielem)
             end do
             kelem = mesh % nelem
             do ielem = 1, mesh_glu % nelem
                kelem                       = kelem + 1
                mesh_tmp % tags(itag) % values(kelem) = mesh_glu % tags(itag) % values(ielem)
             end do
          else if( mesh % tags(itag) % type == NODE_TAG ) then
             do ipoin = 1,mesh % npoin
                mesh_tmp % lninv_loc(ipoin) = mesh % lninv_loc(ipoin)
                mesh_tmp % tags(itag) % values(ipoin) = mesh % tags(itag) % values(ipoin)
             end do
             do ipoin = 1,mesh_glu % npoin
                lid                                 = htalid(htable_3,mesh_glu % lninv_loc(ipoin))
                mesh_tmp % tags(itag) % values(lid) = mesh_glu % tags(itag) % values(ipoin)
             end do
          else if( mesh % tags(itag) % type == BOUNDARY_TAG ) then
             call runend('EXTRACT BOUNDARY FOR TAGS NOT CODED')
          end if
       end do
       !
       ! MESH = MESH_TMP
       !
       call mesh % copy(mesh_tmp,MEMORY_COUNTER=memor_loc)  
       !
       ! Deallocate
       !
       call mesh_glu % deallo(MEMORY_COUNTER)       
       call mesh_tmp % deallo(MEMORY_COUNTER=memor_loc)       
       call htades( htable_3 )

    end if

    call memory_counter_end(memor_loc,mesh % memor,MEMORY_COUNTER)

  end subroutine merge

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Extract a boundary mesh
  !> @details Extract a boundary mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine extract_boundary(mesh_boundary,mesh,MEMORY_COUNTER)

    class(mesh_type_basic),           intent(inout)  :: mesh_boundary       !< Output mesh (eventually deallocated)
    class(mesh_type_basic),           intent(in)     :: mesh                !< Input mesh
    integer(8),             optional, intent(inout)  :: MEMORY_COUNTER(2)   !< Memory counter
    integer(ip),            pointer                  :: list_types(:)
    integer(ip)                                      :: ntypes,itype,ii

    ntypes = element_num_end(mesh % ndime-1)-element_num_ini(mesh % ndime-1)+1
    allocate(list_types(ntypes))
    ii = 0
    do itype = element_num_ini(mesh % ndime-1),element_num_end(mesh % ndime-1)
       ii = ii + 1
       list_types(ii) = itype
    end do

    call mesh_boundary % deallo()
    call mesh_boundary % init()
    call extract(mesh_boundary,mesh,list_types=list_types,MEMORY_COUNTER=MEMORY_COUNTER)
    
    deallocate(list_types)
    
  end subroutine extract_boundary
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Extract a submesh
  !> @details Extract a submesh from MESH using a mask
  !>          Examples of use:
  !>          call mesh_out % extract(mesh,LIST_TYPES=(/TRI03/))
  !>          call mesh_out % extract(mesh,mask)
  !> 
  !-----------------------------------------------------------------------

  subroutine extract(mesh_out,mesh,mask,list_types,MEMORY_COUNTER,mesh_cmp)

    class(mesh_type_basic),                    intent(inout)  :: mesh_out            !< Output mesh
    class(mesh_type_basic),           target,  intent(in)     :: mesh                !< Input mesh
    logical(lg),            optional, pointer, intent(in)     :: mask(:)             !< Mask
    integer(ip),            optional,          intent(in)     :: list_types(:)       !< List of types to extract
    integer(8),             optional,          intent(inout)  :: MEMORY_COUNTER(2)   !< Memory counter
    class(mesh_type_basic), optional,          intent(inout)  :: mesh_cmp            !< Complementary mesh mesh
    integer(ip)                                               :: ipoin,ielem,kelem
    integer(ip)                                               :: inode,pelty,kpoin
    integer(ip)                                               :: pnode,ii,itype
    integer(ip)                                               :: kelem_cmp,kpoin_cmp
    integer(ip)                                               :: mnode_loc,itag
    logical(lg),            pointer                           :: lmask(:)
    integer(ip),            pointer                           :: permr_nodes(:)
    integer(ip),            pointer                           :: permr_cmp(:)
    integer(8)                                                :: memor_loc(2)

    nullify(permr_nodes,permr_cmp)

    call memory_counter_ini(memor_loc,mesh_out % memor,MEMORY_COUNTER)

    if( mesh % nelem == 0 ) then
       !
       ! Just initialize structure and allocate minimum
       !
       call mesh_out % alloca(MEMORY_COUNTER=memor_loc)       
       if( present(mesh_cmp) ) call mesh_cmp % alloca()

    else 
       !
       ! Define mask
       !
       if( present(mask) ) then
          lmask => mask
       else if( present(list_types) ) then
          nullify(lmask)
          call memory_alloca(memor_loc,'LMASK','extract',lmask,mesh % nelem)
          do ii = 1,size(list_types,KIND=ip)
             itype = list_types(ii)
             where( mesh % ltype == itype ) lmask = .true.
          end do
       else
          mesh_out = mesh
          if( present(mesh_cmp) ) call mesh_cmp % init()
          return
       end if
       !  
       ! Permutation
       !
       call memory_alloca(memor_loc,'PERMR_NODES','extract',permr_nodes,mesh % npoin)
       if( present(mesh_cmp) ) call memory_alloca(memor_loc,'PERMR_CMP'  ,'extract',permr_cmp  ,mesh % npoin)
       kpoin     = 0
       kelem     = 0
       kpoin_cmp = 0
       kelem_cmp = 0
       mnode_loc = 0

       do ielem = 1,mesh % nelem

          if( lmask(ielem) ) then

             kelem     = kelem + 1
             pelty     = mesh % ltype(ielem)
             pnode     = element_type(pelty) % number_nodes
             mnode_loc = max(pnode,mnode_loc)
             do inode = 1,pnode
                ipoin = mesh % lnods(inode,ielem)
                if( permr_nodes(ipoin) == 0 ) then
                   kpoin = kpoin + 1
                   permr_nodes(ipoin) = kpoin
                end if
             end do

          else if( present(mesh_cmp) ) then

             kelem_cmp = kelem_cmp + 1
             pelty     = mesh % ltype(ielem)
             pnode     = element_type(pelty) % number_nodes
             mnode_loc = max(pnode,mnode_loc)
             do inode = 1,pnode
                ipoin = mesh % lnods(inode,ielem)
                if( permr_cmp(ipoin) == 0 ) then
                   kpoin_cmp = kpoin_cmp + 1
                   permr_cmp(ipoin) = kpoin_cmp
                end if
             end do

          end if
       end do
       !
       ! New mesh
       !
       if( trim(mesh_out % name) == '' ) mesh_out % name  = trim(mesh % name) // '_extract'
       mesh_out % ndime = mesh % ndime
       mesh_out % nelem = kelem
       mesh_out % npoin = kpoin
       if( mesh_out % mnode == 0 ) mesh_out % mnode = mnode_loc
       mesh_out % ntags = mesh % ntags
       call mesh_out % alloca(MEMORY_COUNTER=memor_loc)

       kelem = 0
       do ielem = 1,mesh % nelem
          if( lmask(ielem) ) then
             pelty                           = mesh % ltype(ielem)
             pnode                           = element_type(pelty) % number_nodes
             kelem                           = kelem + 1
             mesh_out % perme(kelem)         = ielem
             mesh_out % ltype(kelem)         = mesh % ltype    (ielem)
             mesh_out % leinv_loc(kelem)     = mesh % leinv_loc(ielem)             
             mesh_out % lnods(1:pnode,kelem) = permr_nodes(mesh % lnods(1:pnode,ielem))
          end if
       end do
       !
       ! If this both meshes are boundary meshes
       !
       select type ( mesh_out )              
       class is ( bmsh_type_basic )
          select type ( mesh )              
          class is ( bmsh_type_basic )
             kelem = 0
             do ielem = 1,mesh % nelem
                if( lmask(ielem) ) then
                   pelty                     = mesh % ltype(ielem)
                   pnode                     = element_type(pelty) % number_nodes
                   kelem                     = kelem + 1
                   mesh_out % lelbo(kelem)   = mesh % lelbo(ielem)
                   mesh_out % lboel(:,kelem) = mesh % lboel(:,ielem)
                end if
             end do
             mesh_out % mesh => mesh % mesh 
          end select
       end select
       
       do ipoin = 1,mesh % npoin
          kpoin = permr_nodes(ipoin)
          if( kpoin /= 0 ) mesh_out % permn(kpoin) = ipoin
       end do

       do kpoin = 1,mesh_out % npoin
          ipoin                       = mesh_out % permn(kpoin)
          mesh_out % coord(:,kpoin)   = mesh     % coord(:,ipoin)  
          mesh_out % lninv_loc(kpoin) = mesh     % lninv_loc(ipoin)
       end do
       !
       ! Tags
       !
       if( mesh % ntags > 0 ) then
          do itag = 1,mesh % ntags
             mesh_out % tags(itag) % type = mesh % tags(itag) % type
             call mesh_out % alloca_tag(itag,MEMORY_COUNTER=memor_loc)
             if( mesh % tags(itag) % type == ELEMENT_TAG ) then
                kelem = 0
                do ielem = 1,mesh % nelem
                   if( lmask(ielem) ) then
                      kelem = kelem + 1
                      mesh_out % tags(itag) % values(kelem) = mesh % tags(itag) % values(ielem)
                   end if
                end do
             else if( mesh % tags(itag) % type == NODE_TAG ) then
                do kpoin = 1,mesh_out % npoin
                   ipoin = mesh_out % permn(kpoin)
                   mesh_out % tags(itag) % values(kpoin) = mesh % tags(itag) % values(ipoin)
                end do
             else if( mesh % tags(itag) % type == BOUNDARY_TAG ) then
                call runend('EXTRACT BOUNDARY FOR TAGS NOT CODED')
             end if
          end do
       end if       
       
       call memory_deallo(memor_loc,'PERMR_NODES','extract',permr_nodes)      
       !
       ! Complementary mesh
       !
       if( present(mesh_cmp) ) then
          if( trim(mesh_cmp % name) == '' ) mesh_cmp % name  = trim(mesh % name) // '_complement'
          mesh_cmp % ndime = mesh % ndime
          mesh_cmp % nelem = kelem_cmp
          mesh_cmp % npoin = kpoin_cmp
          mesh_cmp % mnode = mnode_loc
          mesh_cmp % ntags = mesh % ntags
          call mesh_cmp % alloca()

          kelem = 0
          do ielem = 1,mesh % nelem
             if( .not. lmask(ielem) ) then
                pelty                           = mesh % ltype(ielem)
                pnode                           = element_type(pelty) % number_nodes
                kelem                           = kelem + 1
                mesh_cmp % perme(kelem)         = ielem
                mesh_cmp % ltype(kelem)         = mesh % ltype    (ielem)
                mesh_cmp % leinv_loc(kelem)     = mesh % leinv_loc(ielem)             
                mesh_cmp % lnods(1:pnode,kelem) = permr_cmp(mesh % lnods(1:pnode,ielem))
             end if
          end do

          do ipoin = 1,mesh % npoin
             kpoin = permr_cmp(ipoin)
             if( kpoin /= 0 ) mesh_cmp % permn(kpoin) = ipoin
          end do

          do kpoin = 1,mesh_cmp % npoin
             ipoin                       = mesh_cmp % permn(kpoin)
             mesh_cmp % coord(:,kpoin)   = mesh % coord(:,ipoin)  
             mesh_cmp % lninv_loc(kpoin) = mesh % lninv_loc(ipoin)
          end do
          !
          ! Tags
          !
          if( mesh % ntags > 0 ) then
             do itag = 1,mesh % ntags
                mesh_cmp % tags(itag) % type = mesh % tags(itag) % type
                call mesh_cmp % alloca_tag(itag)
                if( mesh % tags(itag) % type == ELEMENT_TAG ) then
                   kelem = 0
                   do ielem = 1,mesh % nelem
                      if( .not. lmask(ielem) ) then
                         kelem = kelem + 1
                         mesh_cmp % tags(itag) % values(kelem) = mesh % tags(itag) % values(ielem)
                      end if
                   end do
                else if( mesh % tags(itag) % type == NODE_TAG ) then
                   do kpoin = 1,mesh_cmp % npoin
                      ipoin = mesh_cmp % permn(kpoin)
                      mesh_cmp % tags(itag) % values(kpoin) = mesh % tags(itag) % values(ipoin)
                   end do
                else if( mesh % tags(itag) % type == BOUNDARY_TAG ) then
                   call runend('EXTRACT BOUNDARY FOR TAGS NOT CODED')
                end if
             end do
          end if

          call memory_deallo(memor_loc,'PERMR_CMP','extract',permr_cmp)
          
       end if

    end if
    !
    ! Parent mesh
    !
    mesh_out % parent => mesh
    if( present(mesh_cmp) ) mesh_cmp % parent => mesh
    
    if( present(list_types) ) &
         call memory_deallo(memor_loc,'LMASK','extract',lmask)

    call memory_counter_end(memor_loc,mesh_out % memor,MEMORY_COUNTER)

  end subroutine extract
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Append
  !> @details Append a mesh to another. Duplicated nodes are not
  !>          eliminated
  !>          MESH = MESH // MESH_GLU
  !> 
  !-----------------------------------------------------------------------

  subroutine append(mesh,mesh_glu,MEMORY_COUNTER)

    class(mesh_type_basic),  intent(inout) :: mesh                !< Output mesh
    class(mesh_type_basic),  intent(inout) :: mesh_glu            !< Input mesh (eventually deallocated)
    integer(8),   optional,  intent(inout) :: MEMORY_COUNTER(2)   !< Memory counter
    integer(ip)                            :: ipoin
    integer(ip)                            :: ielem,kelem
    integer(ip)                            :: kpoin,inode
    integer(ip)                            :: npoin_lninv
    integer(ip)                            :: nelem_leinv
    integer(ip)                            :: ielty,itag
    type(mesh_type_basic)                  :: mesh_tmp
    integer(8)                             :: memor_loc(2)

    call memory_counter_ini(memor_loc,mesh % memor,MEMORY_COUNTER)
    
    if( mesh_glu % nelem == 0 ) then

       call mesh_glu % deallo(MEMORY_COUNTER)
       return

    else if( mesh % nelem == 0 ) then
       !
       ! Copy first mesh
       !
       call mesh     % copy(mesh_glu,MEMORY_COUNTER=memor_loc)
       call mesh_glu % deallo(MEMORY_COUNTER)

    else 
       !
       ! Append mesh MESH_GLU to MESH_TMP
       !
       call mesh_tmp % init()
       call mesh_tmp % copy(mesh,MEMORY_COUNTER=memor_loc)  
       mesh % ndime = mesh_tmp % ndime
       mesh % mnode = max(mesh_tmp % mnode,mesh_glu % mnode)
       mesh % nelem = mesh_tmp % nelem + mesh_glu % nelem
       mesh % npoin = mesh_tmp % npoin + mesh_glu % npoin
       mesh % ntags = max(mesh_tmp % ntags,mesh_glu % ntags)
       npoin_lninv  = maxval(mesh % lninv_loc)
       nelem_leinv  = maxval(mesh % leinv_loc)

       call mesh     % deallo(MEMORY_COUNTER=memor_loc)
       call mesh     % alloca(MEMORY_COUNTER=memor_loc)
       do itag = 1,mesh % ntags
          mesh % tags(itag) % type = mesh_glu % tags(itag) % type 
          mesh % tags(itag) % name = mesh_glu % tags(itag) % name 
          call mesh % alloca_tag(itag,MEMORY_COUNTER=memor_loc)
       end do
       call mesh     % copy(mesh_tmp,DO_ALLOCATE=.false.,MEMORY_COUNTER=memor_loc)      
       call mesh_tmp % deallo(MEMORY_COUNTER=memor_loc)
       !
       ! Element arrays
       !
       kelem = mesh_tmp % nelem
       kpoin = mesh_tmp % npoin
       do ielem = 1,mesh_glu % nelem
          kelem = kelem + 1
          ielty = mesh_glu % ltype(ielem)
          do inode = 1,element_type(ielty) % number_nodes
             ipoin = mesh_glu % lnods(inode,ielem)
             mesh % lnods(inode,kelem)  = mesh_glu % lnods(inode,ielem) + kpoin
          end do
          mesh % ltype(kelem)     = mesh_glu % ltype(ielem)
          mesh % leinv_loc(kelem) = mesh_glu % leinv_loc(ielem) !+ nelem_leinv
          mesh % perme(kelem)     = mesh_glu % perme(ielem)
       end do
       !
       ! Node arrays
       !
       do ipoin = 1,mesh_glu % npoin
          kpoin                   = kpoin + 1
          mesh % lninv_loc(kpoin) = mesh_glu % lninv_loc(ipoin) 
          mesh % coord(:,kpoin)   = mesh_glu % coord(:,ipoin)
          mesh % permn(kpoin)     = mesh_glu % permn(ipoin)
       end do 
       !
       ! Tags
       !
       if( mesh % ntags > 0 ) then
          do itag = 1,mesh % ntags
             if( mesh % tags(itag) % type == ELEMENT_TAG ) then
                kelem = mesh_tmp % nelem
                do ielem = 1,mesh_glu % nelem
                   kelem                             = kelem + 1
                   mesh % tags(itag) % values(kelem) = mesh_glu % tags(itag) % values(ielem)
                end do
             else if( mesh % tags(itag) % type == NODE_TAG ) then
                kpoin = mesh_tmp % npoin
                do ipoin = 1,mesh_glu % npoin
                   kpoin                             = kpoin + 1
                   mesh % tags(itag) % values(kpoin) = mesh_glu % tags(itag) % values(ipoin)
                end do
             else if( mesh % tags(itag) % type == BOUNDARY_TAG ) then
                call runend('APPEND BOUNDARY FOR TAGS NOT CODED')
             end if
          end do
       end if              
       !
       ! Deallocate
       !
       call mesh_glu % deallo(MEMORY_COUNTER)       
       call memory_counter_end(memor_loc,mesh % memor,MEMORY_COUNTER)

    end if

  end subroutine append

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-25
  !> @brief   Output
  !> @details Output mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine output(mesh,lunit,filename,material)

    class(mesh_type_basic),                    intent(inout) :: mesh                !< Mesh type
    integer(ip),            optional,          intent(in)    :: lunit               !< Output unit
    character(LEN=*),       optional,          intent(in)    :: filename            !< Filename
    integer(ip),            optional, pointer, intent(in)    :: material(:)
    integer(ip)                                              :: ifirs,inode
    integer(ip)                                              :: ielty,ielem
    integer(ip)                                              :: ipoin,ioerr
    integer(ip)                                              :: ipass,imate
    character(20)                                            :: intost
    character(200)                                           :: dumml
    character(200)                                           :: title
    logical(lg)                                              :: opened
    logical(lg)                                              :: if_output
    integer(4)                                               :: iunit4

    if( mesh % nelem <= 0 ) return

    if( present(filename) ) then
       title = trim(filename)       
    else if( trim(mesh % name) /= '' ) then
       title = adjustl(trim((mesh % name)))
    else
       title = 'MESH'
    end if

    if( present(lunit) ) then
       iunit4 = int(lunit,KIND=4_ip)
    else
       do iunit4 = 90_4,1000_4
          inquire(unit=iunit4,opened=opened,iostat=ioerr)
          if( ioerr /= 0 )  cycle
          if( .not. opened ) exit
       end do
       open(iunit4,file=trim(title)//'.post.msh',status='unknown')
    end if

    if( present(material) ) then
       if( .not. associated(material) ) then
          call runend('DEF_KINTYP_MESH_TYPE: MATERIAL NOT ASSOCIATED')
       else if( size(material,KIND=ip) < mesh % nelem ) then
          call runend('DEF_KINTYP_MESH_TYPE: WRONG MATERIAL SIZE')
       end if
    end if
    
    ifirs = 0
    
    do ipass = 1,2

       do ielty = 1,element_num_end(mesh % ndime)

          if_output = .false.
          if( ipass == 1 ) then
             if( any(abs(mesh % ltype)==ielty) ) if_output = .true.
          else if( ipass == 2 .and. associated(mesh % boundary) ) then
             if( any(abs(mesh % boundary % ltype)==ielty) ) if_output = .true.          
          end if

          if( if_output ) then
             !
             ! Header
             !
             write(intost,*) int(ielty,4)
             dumml = trim(title)//'_'//trim(adjustl(intost))
             write(iunit4,1)&
                  adjustl(trim(dumml)),mesh % ndime,&
                  adjustl(trim(element_type(ielty) % nametopo)),element_type(ielty) % number_nodes
                  !adjustl(trim(cetop(ielty))),element_type(ielty) % number_nodes
             !
             ! Coordinates
             !
             if( ifirs == 0 ) then
                ifirs = 1
                write(iunit4,2) 'coordinates'
                do ipoin = 1,mesh % npoin
                   write(iunit4,3) ipoin,mesh % coord(1:mesh % ndime,ipoin)
                end do
                write(iunit4,2) 'end coordinates'
             end if
             !
             ! Element connectivity
             !
             write(iunit4,2) 'elements'
             if( ipass == 1 ) then
                do ielem = 1,mesh % nelem
                   if( mesh % ltype(ielem) == ielty ) then
                      imate = 1
                      if( present(material) ) then
                         if( associated(material) ) imate = material(ielem)
                      end if
                      write(iunit4,4) ielem,(mesh % lnods(inode,ielem),inode=1,element_type(ielty) % number_nodes),imate
                   end if
                end do
             else
                do ielem = 1,mesh % boundary % nelem
                   if( mesh % boundary % ltype(ielem) == ielty ) then
                      imate = 1
                      if( present(material) ) then
                         if( associated(material) ) imate = material(ielem)
                      end if
                      write(iunit4,4) ielem,(mesh % boundary % lnods(inode,ielem),inode=1,element_type(ielty) % number_nodes)
                   end if
                end do
             end if
             write(iunit4,2) 'end elements'
             write(iunit4,2) ''

          end if

       end do
    end do
    !
    ! Close file
    !
    if( .not. present(lunit) ) close(iunit4)

1   format('MESH ',a,' dimension ',i1,' Elemtype ',a,' Nnode ',i2)
2   format(a)
3   format(i9, 3(1x,e16.8e3))
4   format(i9,50(1x,i9))

  end subroutine output

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-20
  !> @brief   Output results
  !> @details Output results 
  !> 
  !-----------------------------------------------------------------------
  
  subroutine results_r1(mesh,xx,names,lunit,filename,where)

    class(mesh_type_basic),           intent(inout) :: mesh                !< Mesh type
    real(rp),                pointer, intent(in)    :: xx(:)
    character(len=*),                 intent(in)    :: names    
    integer(ip),            optional, intent(in)    :: lunit               !< Output unit
    character(LEN=*),       optional, intent(in)    :: filename            !< Filename
    character(LEN=*),       optional, intent(in)    :: where               !< On nodes or on elements
    integer(4)                                      :: iunit4
    logical(lg)                                     :: on_elements

    if( associated(xx) ) then

      call results_header(mesh,iunit4,names,lunit,filename)
       
       on_elements = .true.
       if( present(where) ) then
          if( trim(where) == 'ON NODES' ) on_elements = .false.
       end if

       if( on_elements ) then
          write(iunit4,2) trim(names),'ALYA',0.0_rp,'Scalar','OnGaussPoints GP'
       else
          write(iunit4,2) trim(names),'ALYA',0.0_rp,'Scalar','OnNodes'
       end if
       write(iunit4,3) trim(names)
       write(iunit4,4)

       call results_single(iunit4,xx)
       write(iunit4,5)  
       
       close(iunit4)
    end if
    
2   format('Result ',a,' ',a,' ',e15.8,' ',a,' ',a)
3   format('ComponentNames ',a)
4   format('Values')
5   format('End values')

  end subroutine results_r1

  subroutine results_r2(mesh,xx,names,lunit,filename,where)

    class(mesh_type_basic),           intent(inout) :: mesh                !< Mesh type
    real(rp),                pointer, intent(in)    :: xx(:,:)
    character(len=5),                 intent(in)    :: names(:)    
    integer(ip),            optional, intent(in)    :: lunit               !< Output unit
    character(LEN=*),       optional, intent(in)    :: filename            !< Filename
    character(LEN=*),       optional, intent(in)    :: where               !< On nodes or on elements
    integer(ip)                                     :: ii
    integer(4)                                      :: iunit4
    logical(lg)                                     :: on_elements

    if( associated(xx) ) then

       call results_header(mesh,iunit4,names,lunit,filename)

       on_elements = .true.
       if( present(where) ) then
          if( trim(where) == 'ON NODES' ) on_elements = .false.
       end if
       
       do ii = 1,size(names,KIND=ip)
          if( on_elements ) then
             write(iunit4,2) trim(names(ii)),'ALYA',0.0_rp,'Scalar','OnGaussPoints GP'
          else
             write(iunit4,2) trim(names(ii)),'ALYA',0.0_rp,'Scalar','OnNodes'
          end if
          write(iunit4,3) trim(names(ii))
          write(iunit4,4)
          call results_single(iunit4,xx(:,ii))
          write(iunit4,5)        
       end do

       close(iunit4)
    end if
    
2   format('Result ',a,' ',a,' ',e15.8,' ',a,' ',a)
3   format('ComponentNames ',a)
4   format('Values')
5   format('End values')

  end subroutine results_r2

  subroutine results_i1(mesh,xx,names,lunit,filename,where)

    class(mesh_type_basic),           intent(inout) :: mesh                !< Mesh type
    integer(ip),             pointer, intent(in)    :: xx(:)
    character(len=5),                 intent(in)    :: names    
    integer(ip),            optional, intent(in)    :: lunit               !< Output unit
    character(LEN=*),       optional, intent(in)    :: filename            !< Filename
    character(LEN=*),       optional, intent(in)    :: where               !< On nodes or on elements
    integer(4)                                      :: iunit4
    logical(lg)                                     :: on_elements
    
    if( associated(xx) ) then

       call results_header(mesh,iunit4,names,lunit,filename)

       on_elements = .true.
       if( present(where) ) then
          if( trim(where) == 'ON NODES' ) on_elements = .false.
       end if
       if( on_elements ) then
          write(iunit4,2) trim(names),'ALYA',0.0_rp,'Scalar','OnGaussPoints GP'
       else
          write(iunit4,2) trim(names),'ALYA',0.0_rp,'Scalar','OnNodes'
       end if
       
       write(iunit4,3) trim(names)
       write(iunit4,4)
       call results_single(iunit4,xx)
       write(iunit4,5)        
       close(iunit4)
       
    end if
    
2   format('Result ',a,' ',a,' ',e15.8,' ',a,' ',a)
3   format('ComponentNames ',a)
4   format('Values')
5   format('End values')

  end subroutine results_i1

  subroutine results_single(iunit4,xx)
    class(*),      intent(in) :: xx(:)
    integer(4),    intent(in) :: iunit4
    integer(ip)               :: ii
    
    select type ( xx )
    type is ( integer(kind=ip) ) ; do ii = 1,size(xx,DIM=1_ip,KIND=ip) ; write(iunit4,1) ii,xx(ii) ; end do
    type is ( real   (kind=rp) ) ; do ii = 1,size(xx,DIM=1_ip,KIND=ip) ; write(iunit4,2) ii,xx(ii) ; end do
    end select
    
1   format(i9, 3(1x,i8))
2   format(i9, 3(1x,e13.6))

  end subroutine results_single

  subroutine results_header(mesh,iunit4,names,lunit,filename)

    class(mesh_type_basic),           intent(inout) :: mesh                !< Mesh type
    integer(4),                       intent(out)   :: iunit4              !< Unit
    character(len=5),                 intent(in)    :: names(*)
    integer(ip),            optional, intent(in)    :: lunit               !< Output unit
    character(LEN=*),       optional, intent(in)    :: filename            !< Filename
    integer(ip)                                     :: ioerr
    character(200)                                  :: title
    logical(lg)                                     :: opened
    logical(lg),            pointer                 :: lexis(:)
    integer(ip)                                     :: nelty,ielem,ielty
    
    if( present(filename) ) then
       title = trim(filename)       
    else if( trim(mesh % name) /= '' ) then
       title = adjustl(trim((mesh % name)))
    else
       title = 'MESH'
    end if
    
    if( present(lunit) ) then
       iunit4 = int(lunit,KIND=4_ip)
    else
       do iunit4 = 90_4,1000_4
          inquire(unit=iunit4,opened=opened,iostat=ioerr)
          if( ioerr /= 0 )  cycle
          if( .not. opened ) exit
       end do
       open(iunit4,file=trim(title)//'.post.res',status='unknown')
    end if

    write(iunit4,1) 
    write(iunit4,*)

    if( mesh % nelem > 0 ) then
       nelty = maxval(mesh % ltype)
       allocate(lexis(nelty))
       do ielty = 1,nelty
          lexis(ielty) = .false.
       end do
       do ielem = 1,mesh % nelem
          lexis(mesh % ltype(ielem)) = .true.
       end do       
       do ielty = 1,nelty
          if(lexis(ielty)) write(iunit4,7) element_type(ielty) % nametopo
       end do       
       deallocate(lexis)
    end if
    
1   format('GiD Post Results File 1.0')
7   format('GaussPoints GP Elemtype ',a,/,&
         & 'Number of Gauss Points: 1',/,&
         & 'Natural Coordinates: Internal',/,&
         & 'End GaussPoints')
   
  end subroutine results_header
    
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-31
  !> @brief   List of boundary nodes
  !> @details List of boundary nodes of a mesh. No parallel exchange!
  !> 
  !-----------------------------------------------------------------------

  subroutine boundary_nodes(mesh,lmask)
    
    class(mesh_type_basic),          intent(in)    :: mesh       
    logical(lg),            pointer, intent(inout) :: lmask(:)
    type(mesh_type_basic)                          :: mesh_boun
    integer(ip)                                    :: ipoin,iboun,inodb,pblty,jpoin
    integer(8)                                     :: memor_loc(2)

    memor_loc = 0_8

    call mesh_boun % init()
    call mesh_boun % boundary_mesh(mesh)
    
    if(.not. associated(lmask)) call memory_alloca(memor_loc,'LMASK','boundary_nodes',lmask,mesh % npoin)

    do ipoin = 1,mesh % npoin
       lmask(ipoin) = .false.
    end do
    
    do iboun = 1,mesh_boun % nelem
       pblty = mesh_boun % ltype(iboun)
       do inodb = 1,element_type(pblty) % number_nodes
          ipoin = mesh_boun % lnods(inodb,iboun)
          jpoin = mesh_boun % permn(ipoin)
          lmask(jpoin) = .true.
       end do
    end do

    call mesh_boun % deallo()
    
  end subroutine boundary_nodes
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-31
  !> @brief   List of boundary nodes
  !> @details List of boundary nodes of a mesh. No parallel exchange!
  !> 
  !-----------------------------------------------------------------------

  subroutine boundary_mesh(mesh_boun,mesh,mask)

    class(mesh_type_basic),           intent(inout) :: mesh_boun
    class(mesh_type_basic),           intent(in)    :: mesh       
    logical(lg),            optional, intent(in)    :: mask(:)
    integer(ip)                                     :: ielty,ielem,iface,inodf
    integer(ip)                                     :: inode,jelem,jface,jelty,ipoin,pnodf
    integer(ip)                                     :: ielpo,mface,mnodf,nboun
    integer(ip)                                     :: kpoin
    integer(8)                                      :: memor_loc(2)
    logical(lg)                                     :: equal_faces  

    integer(ip)                                     :: mepoi,nlelp
    integer(ip), pointer                            :: lelpo(:)
    integer(ip), pointer                            :: pelpo(:)
    integer(ip), pointer                            :: nepoi(:)

    integer(ip)                                     :: nelem
    integer(ip)                                     :: npoin
    integer(ip)                                     :: ndime
    integer(ip)                                     :: mnode
    integer(ip), pointer                            :: lnods(:,:)
    integer(ip), pointer                            :: ltype(:)

    integer(ip)                                     :: nfacg
    integer(ip), pointer                            :: facel(:,:,:)
    integer(ip), pointer                            :: permr(:)
    integer(ip), pointer                            :: invpr(:)
    logical(lg)                                     :: if_mask

    memor_loc = 0_8

    ndime =  mesh % ndime
    npoin =  mesh % npoin
    nelem =  mesh % nelem
    lnods => mesh % lnods
    ltype => mesh % ltype

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

    mface = 0
    mnodf = 0
    do ielem = 1,nelem
       ielty = abs(ltype(ielem))
       mface = max(mface,element_type(ielty) % number_faces)
       mnodf = max(mnodf,element_type(ielty) % max_face_nodes)
    end do

    if_mask = .true.

    call memory_alloca(memor_loc,'NEPOI','boundary_mesh',nepoi,npoin)
    do ielem = 1,nelem
       ielty = ltype(ielem)
       if( present(mask) ) if_mask = mask(ielem)
       if( if_mask ) then
          do inode = 1,element_type(ielty) % number_nodes
             ipoin = lnods(inode,ielem)
             nepoi(ipoin) = nepoi(ipoin) + 1
          end do
       end if
    end do
    call memory_alloca(memor_loc,'PELPO','boundary_mesh',pelpo,npoin+1_ip)
    pelpo(1) = 1
    do ipoin = 1,npoin
       pelpo(ipoin+1) = pelpo(ipoin) + nepoi(ipoin)
    end do
    nlelp = pelpo(npoin+1)
    call memory_alloca(memor_loc,'LELPO','boundary_mesh',lelpo,nlelp)
    do ielem = 1,nelem
       ielty = ltype(ielem)
       if( present(mask) ) if_mask = mask(ielem)
       if( if_mask ) then
          do inode = 1,element_type(ielty) % number_nodes
             ipoin = lnods(inode,ielem)
             lelpo(pelpo(ipoin)) = ielem
             pelpo(ipoin) = pelpo(ipoin)+1
          end do
       end if
    end do
    pelpo(1) =  1
    mepoi    = -1
    do ipoin = 1,npoin
       pelpo(ipoin+1) = pelpo(ipoin) + nepoi(ipoin)
       mepoi = max(mepoi,nepoi(ipoin))
    end do
    call memory_deallo(memor_loc,'NEPOI','boundary_mesh',nepoi)

    !----------------------------------------------------------------------
    !
    ! List of global faces
    !
    !----------------------------------------------------------------------
    !
    ! Allocate memory for lelfa, FACES, CFAEL AND NNODG
    !
    call memory_alloca(memor_loc,'FACEL','boundary_mesh',facel,mnodf+1_ip,mface,nelem)
    !
    ! Construct and sort FACES
    !
    do ielem = 1,nelem                                          
       ielty = abs(ltype(ielem))
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
       if( present(mask) ) if_mask = mask(ielem)
       if( if_mask ) then
          ielty = abs(ltype(ielem))                                  ! eliminate the repited faces
          do iface = 1,element_type(ielty) % number_faces
             if( facel(mnodf+1,iface,ielem) > 0 ) then
                nfacg = nfacg + 1
                ipoin = facel(1,iface,ielem)
                ielpo = pelpo(ipoin)-1
                do while( ielpo < pelpo(ipoin+1)-1 )
                   ielpo = ielpo + 1
                   jelem = lelpo(ielpo)
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
                               ielpo                      =  pelpo(ipoin+1)                     ! Exit JELEM do  
                            end if
                         end if
                      end do
                   end if
                end do
             end if
          end do
       end if
    end do
    call memory_deallo(memor_loc,'LELPO','boundary_mesh',lelpo)
    call memory_deallo(memor_loc,'PELPO','boundary_mesh',pelpo)

    !----------------------------------------------------------------------
    !
    ! Count boundaries and nodes that are involved
    !
    !----------------------------------------------------------------------

    call memory_alloca(memor_loc,'PERMR','boundary_mesh',permr,npoin)
    call memory_alloca(memor_loc,'INVPR','boundary_mesh',invpr,npoin)
    kpoin = 0
    nboun = 0
    mnode = 0
    do ielem = 1,nelem                      
       if( present(mask) ) if_mask = mask(ielem)
       if( if_mask ) then
          ielty = abs(ltype(ielem))                    
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
       end if
    end do

    !----------------------------------------------------------------------
    !
    ! Fill in MESH_BOUN
    !
    !----------------------------------------------------------------------

    if( trim(mesh_boun % name) == '' ) then
       mesh_boun % name  = trim(mesh % name) // '_BOUN'
    end if
    mesh_boun % id    = mesh % id
    mesh_boun % ndime = ndime
    mesh_boun % mnode = mnode
    mesh_boun % npoin = kpoin
    mesh_boun % nelem = nboun
    call mesh_boun % alloca()

    do kpoin = 1,mesh_boun % npoin
       ipoin = invpr(kpoin)
       mesh_boun % coord(1:ndime,kpoin) = mesh % coord(1:ndime,ipoin)
       mesh_boun % permn(kpoin)         = ipoin
    end do

    nboun = 0
    do ielem = 1,nelem                      
       if( present(mask) ) if_mask = mask(ielem)
       if( if_mask ) then
          ielty = abs(ltype(ielem))
          do iface = 1,element_type(ielty) % number_faces
             if( facel(mnodf+1,iface,ielem) > 0 .and. facel(1,iface,ielem) > 0 ) then
                nboun = nboun + 1
                pnodf = element_type(ielty) % node_faces(iface)
                mesh_boun % ltype(nboun) = element_type(ielty) % type_faces(iface)
                mesh_boun % perme(nboun) = ielem
                do inodf = 1,pnodf 
                   inode = element_type(ielty) % list_faces(inodf,iface) 
                   ipoin = lnods(inode,ielem)
                   mesh_boun % lnods(inodf,nboun) = permr(ipoin)
                end do
             end if
          end do
       end if
    end do
    !
    ! Deallocate memory
    !
    call memory_deallo(memor_loc,'FACEL','boundary_mesh',facel)
    call memory_deallo(memor_loc,'PERMR','boundary_mesh',permr)
    call memory_deallo(memor_loc,'INVPR','boundary_mesh',invpr)

  end subroutine boundary_mesh

  !-----------------------------------------------------------------------
  !>
  !> @author  David Oks
  !> @date    2022-03-15
  !> @brief   Permutation array with respect to grandparent mesh
  !> @details Permutation array with respect to grandparent mesh
  !>
  !>              M0: grandparent mesh
  !>              M1: parent mesh
  !>              M2: current (grandchild) mesh
  !>
  !>          Permutation arrays are given by:
  !>
  !>              P1:  M1 -> M0
  !>              P2:  M2 -> M1
  !>              P12: M2 -> M0 ( P12 := P1 o P2 )
  !>
  !>          Which are used as:
  !>
  !>              ipoin0 = P12n( ipoin2 ) = P1n( P2n( ipoin2 ) ) = P1n( ipoin1 )
  !>
  !>          And analogously for elements:
  !>
  !>              ielem0 = P12e( ielem2 ) = P1e( P2e( ielem2 ) ) = P1e( ielem1 )
  !>
  !>          This is useful to get the nodal or elemental indices of a
  !>          grandparent mesh (mesh0) that correpond to each grandchild (mesh2)
  !>          node or element.
  !>
  !-----------------------------------------------------------------------

  subroutine grandpa_permutation(mesh_2, mesh_1, MEMORY_COUNTER)

    class(mesh_type_basic), intent(inout) :: mesh_2
    class(mesh_type_basic), intent(in)    :: mesh_1
    integer(8),   optional, intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                            :: memor_loc(2)
    integer(ip)                           :: ipoin, kpoin, ielem, kelem

    memor_loc = optional_argument((/0_8,0_8/), MEMORY_COUNTER)

    ! Allocate arrays
    if(.not. associated(mesh_2 % permn_grandpa)) then
       call memory_alloca(memor_loc, trim(mesh_2 % name)//' % PERMN_GRANDPA', 'grandpa_permutation', mesh_2 % permn_grandpa, mesh_2 % npoin)
    end if
    if(.not. associated(mesh_2 % perme_grandpa)) then
       call memory_alloca(memor_loc,trim(mesh_2 % name)//' % PERME_GRANDPA','grandpa_permutation',mesh_2 % perme_grandpa,mesh_2 % nelem)
    end if

    ! Get grandparent node permutation array
    do kpoin = 1, mesh_2 % npoin
       ipoin = mesh_1 % permn( mesh_2 % permn(kpoin) )
       mesh_2 % permn_grandpa(kpoin) = ipoin
    end do

    ! Get grandparent element permutation array
    do kelem = 1, mesh_2 % nelem
       ielem = mesh_1 % perme(mesh_2 % perme(kelem))
       mesh_2 % perme_grandpa(kelem) = ielem
    end do

  end subroutine grandpa_permutation

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Cartesian Mesh
  !> @details Extract a submesh from MESH using a mask
  !>          Examples of use:
  !>          call mesh_out % extract(mesh,LIST_TYPES=(/TRI03/))
  !>          call mesh_out % extract(mesh,mask)
  !> 
  !-----------------------------------------------------------------------

  subroutine cartesian_mesh(mesh,ndime,boxes,comin,comax)

    class(mesh_type_basic), intent(inout)  :: mesh     !< Output mesh
    integer(ip),            intent(in)     :: ndime    !< Dimension
    integer(ip),            intent(in)     :: boxes(:) !< # elements
    real(rp),               intent(in)     :: comin(:) !< Min coordinates
    real(rp),               intent(in)     :: comax(:) !< Max coordinates
    integer(ip)                            :: ii,jj,kk
    integer(ip)                            :: ipoin
    integer(ip)                            :: kpoin
    integer(ip)                            :: ielem
    real(rp)                               :: boxer(3)
    real(rp)                               :: delta(3)
    real(rp)                               :: yy,zz

    mesh % ndime = ndime
    if( mesh % ndime == 2 ) then
       mesh % mnode = 4
       mesh % nelem = boxes(1)*boxes(2)
       mesh % npoin = (boxes(1)+1)*(boxes(2)+1)
    else if( mesh % ndime == 3 ) then
       mesh % mnode = 8
       mesh % nelem = boxes(1)*boxes(2)*boxes(3)        
       mesh % npoin = (boxes(1)+1)*(boxes(2)+1)*(boxes(3)+1)
    else
       return
    end if
    delta(1:ndime) = comax(1:ndime)-comin(1:ndime)
    boxer(1:ndime) = real(boxes(1:ndime),rp)
    
    call mesh % alloca()

    ipoin = 0
    ielem = 0
    
    if( ndime == 2 ) then
       do jj = 1,boxes(2)+1
          yy = comin(2) + real(jj-1,rp)/boxer(2) * delta(2)
          do ii = 1,boxes(1)+1
             ipoin = ipoin + 1
             mesh % coord(1,ipoin) = comin(1) + real(ii-1,rp)/boxer(1) * delta(1)
             mesh % coord(2,ipoin) = yy
          end do
       end do
       ipoin = 1
       do jj = 1,boxes(2)
          do ii = 1,boxes(1)
             ielem                 = ielem + 1
             mesh % lnods(1,ielem) = ipoin
             mesh % lnods(2,ielem) = ipoin+1
             mesh % lnods(3,ielem) = ipoin+boxes(1)+2
             mesh % lnods(4,ielem) = ipoin+boxes(1)+1
             mesh % ltype(ielem)   = QUA04
             ipoin                 = ipoin + 1
          end do
          ipoin = ipoin + 1
       end do
    else
       do kk = 1,boxes(3)+1
          zz = comin(3) + real(kk-1,rp)/boxer(3) * delta(3)
          do jj = 1,boxes(2)+1
             yy = comin(2) + real(jj-1,rp)/boxer(2) * delta(2)
             do ii = 1,boxes(1)+1
                ipoin = ipoin + 1
                mesh % coord(1,ipoin) = comin(1) + real(ii-1,rp)/boxer(1) * delta(1)
                mesh % coord(2,ipoin) = yy
                mesh % coord(3,ipoin) = zz
             end do
          end do
       end do
       ipoin = 1
       kpoin = (boxes(1)+1)*(boxes(2)+1)
       do kk = 1,boxes(3)
          do jj = 1,boxes(2)
             do ii = 1,boxes(1)
                ielem                 = ielem + 1
                
                mesh % lnods(1,ielem) = ipoin
                mesh % lnods(2,ielem) = ipoin+1
                mesh % lnods(3,ielem) = ipoin+boxes(1)+2
                mesh % lnods(4,ielem) = ipoin+boxes(1)+1
                mesh % lnods(5,ielem) = kpoin + ipoin
                mesh % lnods(6,ielem) = kpoin + ipoin+1
                mesh % lnods(7,ielem) = kpoin + ipoin+boxes(1)+2
                mesh % lnods(8,ielem) = kpoin + ipoin+boxes(1)+1
                
                mesh % ltype(ielem)   = HEX08
                ipoin                 = ipoin + 1
             end do
             ipoin = ipoin + 1
          end do
          ipoin = kpoin * kk + 1
       end do
    end if
    
  end subroutine cartesian_mesh

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Iag 
  !> @details Assign identity tag to mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine identity(mesh)

    class(mesh_type_basic), intent(inout) :: mesh
    integer(ip)                           :: ielem,ipoin
    
    do ielem = 1,mesh % nelem
       mesh % leinv_loc(ielem) = ielem
    end do
    do ipoin = 1,mesh % npoin
       mesh % lninv_loc(ipoin) = ipoin
    end do
    
  end subroutine identity

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Iag 
  !> @details Assign identity tag to mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine print_info(mesh)

    class(mesh_type_basic), intent(inout) :: mesh
    
    write(*,*) 'NAME=  ',trim(mesh % name)
    write(*,*) 'ID=    ',mesh % id
    write(*,*) 'NDIME= ',mesh % ndime
    write(*,*) 'MNODE= ',mesh % mnode
    write(*,*) 'NELEM= ',mesh % nelem
    write(*,*) 'NPOIN= ',mesh % npoin
    write(*,*) 'LNODS= ',memory_size(mesh % lnods,1_ip),memory_size(mesh % lnods,2_ip)
    write(*,*) 'LTYPE= ',memory_size(mesh % ltype,1_ip)
    write(*,*) 'COORD= ',memory_size(mesh % coord,1_ip),memory_size(mesh % coord,2_ip)
    
  end subroutine print_info

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-18
  !> @brief   Point to a mesh
  !> @details Point to a mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine point_to(mesh,mesh_target)
    
    class(mesh_type_basic), intent(inout) :: mesh
    class(mesh_type_basic), intent(inout) :: mesh_target

    mesh % ndime                 =   mesh_target % ndime
    mesh % mnode                 =   mesh_target % mnode
    mesh % npoin                 =   mesh_target % npoin
    mesh % nelem                 =   mesh_target % nelem
    mesh % lnods                 =>  mesh_target % lnods
    mesh % ltype                 =>  mesh_target % ltype
    mesh % leinv_loc             =>  mesh_target % leinv_loc
    mesh % lninv_loc             =>  mesh_target % lninv_loc
    mesh % coord                 =>  mesh_target % coord

    mesh % comm % RANK4          =   mesh_target % comm % RANK4
    mesh % comm % SIZE4          =   mesh_target % comm % SIZE4
    mesh % comm % PAR_COMM_WORLD =   mesh_target % comm % PAR_COMM_WORLD
    mesh % comm % bound_dim      =   mesh_target % comm % bound_dim
    mesh % comm % nneig          =   mesh_target % comm % nneig
    mesh % comm % neights        =>  mesh_target % comm % neights
    mesh % comm % bound_size     =>  mesh_target % comm % bound_size
    mesh % comm % bound_perm     =>  mesh_target % comm % bound_perm

  end subroutine point_to

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-03
  !> @brief   Compute the cut mesh given a level field
  !> @details Cut a mesh when the level changes sign
  !> 
  !-----------------------------------------------------------------------
  
  subroutine cut_level(mesh_out,mesh_in,level,MEMORY_COUNTER)

    class(mesh_type_basic),          intent(inout) :: mesh_out
    class(mesh_type_basic),          intent(in)    :: mesh_in
    real(rp),               pointer, intent(in)    :: level(:)
    integer(8),   optional,          intent(inout) :: MEMORY_COUNTER(2)
    integer(ip)                                    :: ielem,pelty,pnode
    integer(ip)                                    :: inode,jnode,ipoin
    integer(ip)                                    :: ii,ij,inod1,inod2        
    integer(ip)                                    :: signn,signp,sigtn
    integer(ip)                                    :: pdime,iedge,sigtp   
    integer(ip)                                    :: num_tet,compt,compl
    integer(ip)                                    :: ledgi(10),compg
    real(rp)                                       :: elcod(mesh_in % ndime,mesh_in % mnode)
    real(rp)                                       :: ellev(mesh_in % mnode),inter(3,4)
    real(rp)                                       :: l1,lp,x1,y1,x2,y2,z1,z2
    real(rp)                                       :: ledgr(10)
    logical(lg), pointer                           :: ifcut(:)
    integer(ip), pointer                           :: XXXXX_TO_TET04(:,:)
    integer(8)                                     :: memor_loc(2)
    character(9), parameter                        :: vacal='cut_level'

    if( present(MEMORY_COUNTER) ) then
       memor_loc = MEMORY_COUNTER
    else
       memor_loc = 0_8
    end if

    nullify(ifcut)

    call memory_alloca(memor_loc,'IFCUT',vacal,ifcut,mesh_in % nelem)
    compt = 0_ip
    pdime = mesh_in % ndime
    !
    ! Count number of cut elements
    !
    do ielem = 1,mesh_in % nelem

       pelty = mesh_in % ltype(ielem)
       pnode = element_type(pelty) % number_nodes
       do inode = 1,pnode
          ipoin        = mesh_in % lnods(inode,ielem)
          ellev(inode) = level(ipoin)
       end do

       if( pdime == 2 ) then
          !
          ! 2D elements
          !
          compl = 0_ip
          do iedge = 1,element_type(pelty) % number_edges
             inode = element_type(pelty) % list_edges(1,iedge)
             jnode = element_type(pelty) % list_edges(2,iedge)
             if( ellev(inode)*ellev(jnode) <= 0.0_rp ) compl = compl + 1
          end do

          if( compl >= 2 ) then
             compt = compt + 1
             ifcut(ielem) = .true. 
          end if

       else 
          !
          ! 3D elements
          !
          signn = 0_ip
          signp = 0_ip        
          do inode = 1,pnode
             if( ellev(inode) >= 0.0_rp) then
                signp = signp+1
             else 
                signn = signn+1
             end if
          end do
          if( signp /= pnode .and. signn /= pnode ) then
             if(      pelty == HEX08 ) then
                XXXXX_TO_TET04 => HEX08_TO_TET04
             else if( pelty == PEN06 ) then
                XXXXX_TO_TET04 => PEN06_TO_TET04
             else if( pelty == PYR05 ) then
                XXXXX_TO_TET04 => PYR05_TO_TET04
             else if( pelty == TET04 ) then
                XXXXX_TO_TET04 => TET04_TO_TET04
             else
                call runend('ELEMENT NOT CODED')
             end if
             num_tet = size(XXXXX_TO_TET04,DIM=2_ip,KIND=ip)
             !
             ! COMPT = Number of interface triangle
             !
             do ii = 1,num_tet
                sigtn = 0_ip
                sigtp = 0_ip                 
                do ij = 1,4
                   if(ellev(XXXXX_TO_TET04(ij,ii))>=0.0_rp) then
                      sigtp=sigtp+1
                   else 
                      sigtn=sigtn+1
                   endif
                end do
                if(      sigtp == 1 .or. sigtn == 1 ) then
                   compt        = compt + 1
                   ifcut(ielem) = .true.
                else if( sigtp == 2 .or. sigtn == 2 ) then
                   compt        = compt + 2
                   ifcut(ielem) = .true.
                endif
             end do
          end if
       end if
    end do
    !
    ! Allocate mesh
    !
    call mesh_out % init()
    mesh_out % nelem = compt
    mesh_out % npoin = pdime * compt
    mesh_out % ndime = pdime
    mesh_out % mnode = pdime
    call mesh_out % alloca() 

    compt = 0_ip
    compg = 0_ip

    do ielem = 1,mesh_in % nelem

       if( ifcut(ielem) ) then

          pelty = mesh_in % ltype(ielem)
          pnode = element_type(pelty) % number_nodes
          do inode = 1,pnode
             ipoin                = mesh_in % lnods(inode,ielem)
             ellev(inode)         = level(ipoin)
             elcod(1:pdime,inode) = mesh_in % coord(1:pdime,ipoin)
          end do

          if( pdime == 2 ) then

             compl = 0_ip           
             do iedge = 1,element_type(pelty) % number_edges
                inode        = element_type(pelty) % list_edges(1,iedge)
                jnode        = element_type(pelty) % list_edges(2,iedge)
                if( ellev(inode) * ellev(jnode) <= 0.0_rp ) then
                   compl        = compl + 1
                   ledgi(compl) = iedge
                   ledgr(compl) = abs(ellev(inode) * ellev(jnode))
                end if
             end do

             call maths_heap_sort(1_ip,compl,ledgr,ledgi)
             !
             ! Compute the intersection of the elements with the surface 
             !
             compt = compt + 1
             do compl = 1,2
                compg                         = compg + 1
                iedge                         = ledgi(compl)
                inode                         = element_type(pelty) % list_edges(1,iedge)
                jnode                         = element_type(pelty) % list_edges(2,iedge)              
                l1                            = abs(ellev(inode))
                lp                            = abs(ellev(inode)-ellev(jnode))
                x1                            = elcod(1,inode)
                x2                            = elcod(1,jnode)
                y1                            = elcod(2,inode)
                y2                            = elcod(2,jnode)
                mesh_out % lnods(compl,compt) = compg
                mesh_out % coord(1,compg)     = x1 * (1.0_rp-l1/lp) + x2 * l1 / lp  
                mesh_out % coord(2,compg)     = y1 * (1.0_rp-l1/lp) + y2 * l1 / lp
             end do

          else 

             signn = 0_ip
             signp = 0_ip
             compl = 0_ip
             inod1 = 0_ip

             do inode = 1,pnode
                if( ellev(inode) >= 0.0_rp ) then
                   signp = signp+1
                else 
                   signn = signn+1
                end if
             end do
             if(      pelty == HEX08 ) then
                XXXXX_TO_TET04 => HEX08_TO_TET04
             else if( pelty == PEN06 ) then
                XXXXX_TO_TET04 => PEN06_TO_TET04
             else if( pelty == PYR05 ) then
                XXXXX_TO_TET04 => PYR05_TO_TET04
             else if( pelty == TET04 ) then
                XXXXX_TO_TET04 => TET04_TO_TET04
             else
                call runend('ELEMENT NOT CODED')
             end if

             if( pelty == TET04 ) then
                !
                ! TET04
                !
                if(signn==1) then

                   !  research of one interface triangle
                   do ij=1,pnode
                      if(ellev(ij)<0.0_rp) then
                         inod1=ij
                      endif
                   end do

                   do ij=1,pnode
                      if(ij/=inod1) then
                         compg=compg+1
                         compl=compl+1
                         if(compl==1) then
                            compt=compt+1
                         endif
                         !
                         ! Compute the intersection of the elements with the surface 
                         !
                         l1=abs(ellev(inod1))
                         lp=abs(ellev(inod1)-ellev(ij))
                         x1=elcod(1,inod1)
                         x2=elcod(1,ij)
                         y1=elcod(2,inod1)
                         y2=elcod(2,ij)
                         z1=elcod(3,inod1)
                         z2=elcod(3,ij)

                         mesh_out % lnods(compl,compt) = compg
                         mesh_out % coord(1,compg)     = x1*(1-l1/lp)+x2*l1/lp  
                         mesh_out % coord(2,compg)     = y1*(1-l1/lp)+y2*l1/lp 
                         mesh_out % coord(3,compg)     = z1*(1-l1/lp)+z2*l1/lp 

                      endif

                   end do

                else if(signp==1) then
                   !  research of one interface triangle
                   do ij=1,pnode
                      if(ellev(ij)>=0.0_rp) then
                         inod1=ij
                      endif
                   end do

                   do ij=1,pnode
                      if(ij/=inod1) then
                         compg=compg+1
                         compl=compl+1
                         if(compl==1) then
                            compt=compt+1
                         endif
                         !
                         ! Compute the intersection of the elements with the surface 
                         !
                         l1=abs(ellev(inod1))
                         lp=abs(ellev(inod1)-ellev(ij))
                         x1=elcod(1,inod1)
                         x2=elcod(1,ij)
                         y1=elcod(2,inod1)
                         y2=elcod(2,ij)
                         z1=elcod(3,inod1)
                         z2=elcod(3,ij)

                         mesh_out % lnods(compl,compt) = compg
                         mesh_out % coord(1,compg)     = x1*(1-l1/lp)+x2*l1/lp  
                         mesh_out % coord(2,compg)     = y1*(1-l1/lp)+y2*l1/lp 
                         mesh_out % coord(3,compg)     = z1*(1-l1/lp)+z2*l1/lp 

                      endif
                   end do

                else if(signp==2) then

                   !  research of two interface triangles
                   do ij=1,4
                      if(ellev(ij)<0.0_rp) then
                         if(inod1==0) then
                            inod1=ij
                         else 
                            inod2=ij
                         endif
                      endif
                   end do

                   do ij=1,4

                      if(ij/=inod1.and.ij/=inod2) then
                         !
                         ! Compute the intersection of the elements with the surface 
                         !

                         compl=compl+1
                         l1=abs(ellev(inod1))
                         lp=abs(ellev(inod1)-ellev(ij))
                         x1=elcod(1,inod1)
                         x2=elcod(1,ij)
                         y1=elcod(2,inod1)
                         y2=elcod(2,ij)
                         z1=elcod(3,inod1)
                         z2=elcod(3,ij)

                         inter(1,compl)=  x1*(1-l1/lp)+x2*l1/lp  
                         inter(2,compl)=  y1*(1-l1/lp)+y2*l1/lp 
                         inter(3,compl)=  z1*(1-l1/lp)+z2*l1/lp 


                         compl=compl+1
                         l1=abs(ellev(inod2))
                         lp=abs(ellev(inod2)-ellev(ij))
                         x1=elcod(1,inod2)
                         x2=elcod(1,ij)
                         y1=elcod(2,inod2)
                         y2=elcod(2,ij)
                         z1=elcod(3,inod2)
                         z2=elcod(3,ij)

                         inter(1,compl)=  x1*(1-l1/lp)+x2*l1/lp  
                         inter(2,compl)=  y1*(1-l1/lp)+y2*l1/lp 
                         inter(3,compl)=  z1*(1-l1/lp)+z2*l1/lp 

                      endif

                   end do

                   compg                     = compg+1
                   compt                     = compt+1
                   mesh_out % lnods(1,compt) = compg
                   mesh_out % coord(1,compg) = inter(1,1)
                   mesh_out % coord(2,compg) = inter(2,1)
                   mesh_out % coord(3,compg) = inter(3,1) 

                   compg                     = compg+1
                   mesh_out % lnods(2,compt) = compg
                   mesh_out % coord(1,compg) = inter(1,2)
                   mesh_out % coord(2,compg) = inter(2,2)
                   mesh_out % coord(3,compg) = inter(3,2) 

                   compg                     = compg+1
                   mesh_out % lnods(3,compt) = compg
                   mesh_out % coord(1,compg) = inter(1,3)
                   mesh_out % coord(2,compg) = inter(2,3)
                   mesh_out % coord(3,compg) = inter(3,3)

                   compg                     = compg+1
                   compt                     = compt+1
                   mesh_out % lnods(1,compt) = compg
                   mesh_out % coord(1,compg) = inter(1,2)
                   mesh_out % coord(2,compg) = inter(2,2)
                   mesh_out % coord(3,compg) = inter(3,2) 

                   compg                     = compg+1
                   mesh_out % lnods(2,compt) = compg
                   mesh_out % coord(1,compg) = inter(1,3)
                   mesh_out % coord(2,compg) = inter(2,3)
                   mesh_out % coord(3,compg) = inter(3,3) 

                   compg                     = compg+1
                   mesh_out % lnods(3,compt) = compg
                   mesh_out % coord(1,compg) = inter(1,4)
                   mesh_out % coord(2,compg) = inter(2,4)
                   mesh_out % coord(3,compg) = inter(3,4) 

                end if

             else if( pelty == HEX08 ) then
                !
                ! HEX08
                !
                if( signp /= pnode .and. signn /= pnode ) then

                   do ii= 1,6
                      sigtn = 0_ip
                      sigtp = 0_ip
                      compl = 0_ip
                      inod1 = 0_ip
                      
                      do ij = 1,4
                         if( ellev(XXXXX_TO_TET04(ij,ii)) >= 0.0_rp ) then
                            sigtp = sigtp + 1
                         else 
                            sigtn = sigtn + 1
                         endif
                      end do

                      if(sigtn==1) then
                         !  research of one interface triangle
                         do ij=1,4
                            if(ellev(XXXXX_TO_TET04(ij,ii))<0.0_rp) then
                               inod1=ij
                            endif
                         end do

                         do ij = 1,4
                            if( ij /= inod1 ) then
                               compg = compg+1
                               compl = compl+1
                               if( compl == 1 ) then
                                  compt = compt + 1
                               endif                               
                               !
                               ! Compute the intersection of the elements with the surface 
                               !
                               l1 = abs(ellev(XXXXX_TO_TET04(inod1,ii)))
                               lp = abs(ellev(XXXXX_TO_TET04(inod1,ii))-ellev(XXXXX_TO_TET04(ij,ii)))
                               x1 = elcod(1,XXXXX_TO_TET04(inod1,ii))
                               x2 = elcod(1,XXXXX_TO_TET04(ij,ii))
                               y1 = elcod(2,XXXXX_TO_TET04(inod1,ii))
                               y2 = elcod(2,XXXXX_TO_TET04(ij,ii))
                               z1 = elcod(3,XXXXX_TO_TET04(inod1,ii))
                               z2 = elcod(3,XXXXX_TO_TET04(ij,ii))

                               mesh_out % lnods(compl,compt) = compg
                               mesh_out % coord(1,compg)     = x1 * (1.0_rp-l1/lp) + x2 * l1 / lp  
                               mesh_out % coord(2,compg)     = y1 * (1.0_rp-l1/lp) + y2 * l1 / lp 
                               mesh_out % coord(3,compg)     = z1 * (1.0_rp-l1/lp) + z2 * l1 / lp 
                            endif

                         end do

                      else if(sigtp==1) then
                         !  research of one interface triangle
                         do ij=1,4
                            if(ellev(XXXXX_TO_TET04(ij,ii))>=0.0_rp) then
                               inod1=ij
                            endif
                         end do

                         do ij=1,4
                            if(ij/=inod1) then
                               compg=compg+1
                               compl=compl+1
                               if(compl==1) then
                                  compt=compt+1
                               endif
                               !
                               ! Compute the intersection of the elements with the surface 
                               !
                               l1=abs(ellev(XXXXX_TO_TET04(inod1,ii)))
                               lp=abs(ellev(XXXXX_TO_TET04(inod1,ii))-ellev(XXXXX_TO_TET04(ij,ii)))
                               x1=elcod(1,XXXXX_TO_TET04(inod1,ii))
                               x2=elcod(1,XXXXX_TO_TET04(ij,ii))
                               y1=elcod(2,XXXXX_TO_TET04(inod1,ii))
                               y2=elcod(2,XXXXX_TO_TET04(ij,ii))
                               z1=elcod(3,XXXXX_TO_TET04(inod1,ii))
                               z2=elcod(3,XXXXX_TO_TET04(ij,ii))

                               mesh_out % lnods(compl,compt)= compg
                               mesh_out % coord(1,compg)=  x1 * (1.0_rp-l1/lp) + x2 * l1 / lp  
                               mesh_out % coord(2,compg)=  y1 * (1.0_rp-l1/lp) + y2 * l1 / lp 
                               mesh_out % coord(3,compg)=  z1 * (1.0_rp-l1/lp) + z2 * l1 / lp 

                            endif
                         end do

                      else if(sigtn==2) then

                         !  research of two interface triangles
                         do ij=1,4
                            if(ellev(XXXXX_TO_TET04(ij,ii))<0.0_rp) then
                               if(inod1==0) then
                                  inod1=ij
                               else 
                                  inod2=ij
                               endif
                            endif
                         end do

                         do ij=1,4
                            if(ij/=inod1.and.ij/=inod2) then
                               !
                               ! Compute the intersection of the elements with the surface 
                               !

                               compl=compl+1
                               l1=abs(ellev(XXXXX_TO_TET04(inod1,ii)))
                               lp=abs(ellev(XXXXX_TO_TET04(inod1,ii))-ellev(XXXXX_TO_TET04(ij,ii)))
                               x1=elcod(1,XXXXX_TO_TET04(inod1,ii))
                               x2=elcod(1,XXXXX_TO_TET04(ij,ii))
                               y1=elcod(2,XXXXX_TO_TET04(inod1,ii))
                               y2=elcod(2,XXXXX_TO_TET04(ij,ii))
                               z1=elcod(3,XXXXX_TO_TET04(inod1,ii))
                               z2=elcod(3,XXXXX_TO_TET04(ij,ii))

                               !mesh_out % lnods(compl,compt)= compg
                               !mesh_out % coord(1,compg)=  x1*(1-l1/lp)+x2*l1/lp  
                               !mesh_out % coord(2,compg)=  y1*(1-l1/lp)+y2*l1/lp 
                               !mesh_out % coord(3,compg)=  z1*(1-l1/lp)+z2*l1/lp
                               
                               inter(1,compl)=  x1*(1-l1/lp)+x2*l1/lp  
                               inter(2,compl)=  y1*(1-l1/lp)+y2*l1/lp 
                               inter(3,compl)=  z1*(1-l1/lp)+z2*l1/lp 

                               compl=compl+1
                               l1=abs(ellev(XXXXX_TO_TET04(inod2,ii)))
                               lp=abs(ellev(XXXXX_TO_TET04(inod2,ii))-ellev(XXXXX_TO_TET04(ij,ii)))
                               x1=elcod(1,XXXXX_TO_TET04(inod2,ii))
                               x2=elcod(1,XXXXX_TO_TET04(ij,ii))
                               y1=elcod(2,XXXXX_TO_TET04(inod2,ii))
                               y2=elcod(2,XXXXX_TO_TET04(ij,ii))
                               z1=elcod(3,XXXXX_TO_TET04(inod2,ii))
                               z2=elcod(3,XXXXX_TO_TET04(ij,ii))

                               !mesh_out % lnods(compl,compt)= compg
                               !mesh_out % coord(1,compg)= x1*(1-l1/lp)+x2*l1/lp     
                               !mesh_out % coord(2,compg)= y1*(1-l1/lp)+y2*l1/lp   
                               !mesh_out % coord(3,compg)= z1*(1-l1/lp)+z2*l1/lp  
                               
                               inter(1,compl)=  x1*(1-l1/lp)+x2*l1/lp  
                               inter(2,compl)=  y1*(1-l1/lp)+y2*l1/lp 
                               inter(3,compl)=  z1*(1-l1/lp)+z2*l1/lp 

                            endif

                         end do

                         compg                     = compg+1
                         compt                     = compt+1
                         
                         mesh_out % lnods(1,compt) = compg
                         mesh_out % coord(1,compg) = inter(1,1)
                         mesh_out % coord(2,compg) = inter(2,1)
                         mesh_out % coord(3,compg) = inter(3,1) 

                         compg                     = compg+1
                         mesh_out % lnods(2,compt) = compg
                         mesh_out % coord(1,compg) = inter(1,2)
                         mesh_out % coord(2,compg) = inter(2,2)
                         mesh_out % coord(3,compg) = inter(3,2) 

                         compg                     = compg+1
                         mesh_out % lnods(3,compt) = compg
                         mesh_out % coord(1,compg) = inter(1,3)
                         mesh_out % coord(2,compg) = inter(2,3)
                         mesh_out % coord(3,compg) = inter(3,3)

                         compg                     = compg+1
                         compt                     = compt+1
                         
                         mesh_out % lnods(1,compt) = compg
                         mesh_out % coord(1,compg) = inter(1,2)
                         mesh_out % coord(2,compg) = inter(2,2)
                         mesh_out % coord(3,compg) = inter(3,2) 

                         compg                     = compg+1
                         mesh_out % lnods(2,compt) = compg
                         mesh_out % coord(1,compg) = inter(1,3)
                         mesh_out % coord(2,compg) = inter(2,3)
                         mesh_out % coord(3,compg) = inter(3,3) 

                         compg                     = compg+1
                         mesh_out % lnods(3,compt) = compg
                         mesh_out % coord(1,compg) = inter(1,4)
                         mesh_out % coord(2,compg) = inter(2,4)
                         mesh_out % coord(3,compg) = inter(3,4) 

                      end if

                   end do

                end if

             else if( pelty == PEN06 ) then
                !
                ! PEN06
                !
                if( signp /= pnode .and. signn /= pnode ) then

                   do ii=1,3
                      sigtn=0_ip
                      sigtp=0_ip
                      compl=0_ip
                      inod1=0_ip

                      do ij=1,4
                         if(ellev(XXXXX_TO_TET04(ij,ii))>=0.0_rp) then
                            sigtp=sigtp+1
                         else 
                            sigtn=sigtn+1
                         endif
                      end do

                      if(sigtn==1) then
                         !  research of one interface triangle
                         do ij=1,4
                            if(ellev(XXXXX_TO_TET04(ij,ii))<0.0_rp) then
                               inod1=ij
                            endif
                         end do

                         do ij=1,4
                            if(ij/=inod1) then
                               compg=compg+1
                               compl=compl+1
                               if(compl==1) then
                                  compt=compt+1
                               endif
                               !
                               ! Compute the intersection of the elements with the surface 
                               !
                               l1=abs(ellev(XXXXX_TO_TET04(inod1,ii)))
                               lp=abs(ellev(XXXXX_TO_TET04(inod1,ii))-ellev(XXXXX_TO_TET04(ij,ii)))
                               x1=elcod(1,XXXXX_TO_TET04(inod1,ii))
                               x2=elcod(1,XXXXX_TO_TET04(ij,ii))
                               y1=elcod(2,XXXXX_TO_TET04(inod1,ii))
                               y2=elcod(2,XXXXX_TO_TET04(ij,ii))
                               z1=elcod(3,XXXXX_TO_TET04(inod1,ii))
                               z2=elcod(3,XXXXX_TO_TET04(ij,ii))

                               mesh_out % lnods(compl,compt) = compg
                               mesh_out % coord(1,compg)     =  x1 * (1.0_rp-l1/lp) + x2 * l1 / lp  
                               mesh_out % coord(2,compg)     =  y1 * (1.0_rp-l1/lp) + y2 * l1 / lp 
                               mesh_out % coord(3,compg)     =  z1 * (1.0_rp-l1/lp) + z2 * l1 / lp 

                            endif

                         end do

                      else if(sigtp==1) then
                         !  research of one interface triangle
                         do ij=1,4
                            if(ellev(XXXXX_TO_TET04(ij,ii))>=0.0_rp) then
                               inod1=ij
                            endif
                         end do

                         do ij=1,4
                            if(ij/=inod1) then
                               compg=compg+1
                               compl=compl+1
                               if(compl==1) then
                                  compt=compt+1
                               endif
                               !
                               ! Compute the intersection of the elements with the surface 
                               !
                               l1=abs(ellev(XXXXX_TO_TET04(inod1,ii)))
                               lp=abs(ellev(XXXXX_TO_TET04(inod1,ii))-ellev(XXXXX_TO_TET04(ij,ii)))
                               x1=elcod(1,XXXXX_TO_TET04(inod1,ii))
                               x2=elcod(1,XXXXX_TO_TET04(ij,ii))
                               y1=elcod(2,XXXXX_TO_TET04(inod1,ii))
                               y2=elcod(2,XXXXX_TO_TET04(ij,ii))
                               z1=elcod(3,XXXXX_TO_TET04(inod1,ii))
                               z2=elcod(3,XXXXX_TO_TET04(ij,ii))

                               mesh_out % lnods(compl,compt) = compg
                               mesh_out % coord(1,compg)     = x1 * (1.0_rp-l1/lp) + x2 * l1 / lp  
                               mesh_out % coord(2,compg)     = y1 * (1.0_rp-l1/lp) + y2 * l1 / lp 
                               mesh_out % coord(3,compg)     = z1 * (1.0_rp-l1/lp) + z2 * l1 / lp 

                            endif
                         end do

                      else if(sigtn==2) then

                         !  research of two interface triangles
                         do ij=1,4
                            if(ellev(XXXXX_TO_TET04(ij,ii))<0.0_rp) then
                               if(inod1==0) then
                                  inod1=ij
                               else 
                                  inod2=ij
                               endif
                            endif
                         end do

                         do ij=1,4
                            if(ij/=inod1.and.ij/=inod2) then
                               !
                               ! Compute the intersection of the elements with the surface 
                               !

                               compl=compl+1
                               l1=abs(ellev(XXXXX_TO_TET04(inod1,ii)))
                               lp=abs(ellev(XXXXX_TO_TET04(inod1,ii))-ellev(XXXXX_TO_TET04(ij,ii)))
                               x1=elcod(1,XXXXX_TO_TET04(inod1,ii))
                               x2=elcod(1,XXXXX_TO_TET04(ij,ii))
                               y1=elcod(2,XXXXX_TO_TET04(inod1,ii))
                               y2=elcod(2,XXXXX_TO_TET04(ij,ii))
                               z1=elcod(3,XXXXX_TO_TET04(inod1,ii))
                               z2=elcod(3,XXXXX_TO_TET04(ij,ii))

                               inter(1,compl)=  x1*(1-l1/lp)+x2*l1/lp  
                               inter(2,compl)=  y1*(1-l1/lp)+y2*l1/lp 
                               inter(3,compl)=  z1*(1-l1/lp)+z2*l1/lp 

                               compl=compl+1
                               l1=abs(ellev(XXXXX_TO_TET04(inod2,ii)))
                               lp=abs(ellev(XXXXX_TO_TET04(inod2,ii))-ellev(XXXXX_TO_TET04(ij,ii)))
                               x1=elcod(1,XXXXX_TO_TET04(inod2,ii))
                               x2=elcod(1,XXXXX_TO_TET04(ij,ii))
                               y1=elcod(2,XXXXX_TO_TET04(inod2,ii))
                               y2=elcod(2,XXXXX_TO_TET04(ij,ii))
                               z1=elcod(3,XXXXX_TO_TET04(inod2,ii))
                               z2=elcod(3,XXXXX_TO_TET04(ij,ii))

                               inter(1,compl)=  x1*(1-l1/lp)+x2*l1/lp  
                               inter(2,compl)=  y1*(1-l1/lp)+y2*l1/lp 
                               inter(3,compl)=  z1*(1-l1/lp)+z2*l1/lp 

                            endif

                         end do

                         compg                     = compg+1
                         compt                     = compt+1

                         mesh_out % lnods(1,compt) = compg
                         mesh_out % coord(1,compg) = inter(1,1)
                         mesh_out % coord(2,compg) = inter(2,1)
                         mesh_out % coord(3,compg) = inter(3,1) 

                         compg                     = compg+1
                         mesh_out % lnods(2,compt) = compg
                         mesh_out % coord(1,compg) = inter(1,2)
                         mesh_out % coord(2,compg) = inter(2,2)
                         mesh_out % coord(3,compg) = inter(3,2) 

                         compg                     = compg+1
                         mesh_out % lnods(3,compt) = compg
                         mesh_out % coord(1,compg) = inter(1,3)
                         mesh_out % coord(2,compg) = inter(2,3)
                         mesh_out % coord(3,compg) = inter(3,3) 

                         compg                     = compg+1
                         compt                     = compt+1

                         mesh_out % lnods(1,compt) = compg
                         mesh_out % coord(1,compg) = inter(1,2)
                         mesh_out % coord(2,compg) = inter(2,2)
                         mesh_out % coord(3,compg) = inter(3,2) 

                         compg                     = compg+1
                         mesh_out % lnods(2,compt) = compg
                         mesh_out % coord(1,compg) = inter(1,3)
                         mesh_out % coord(2,compg) = inter(2,3)
                         mesh_out % coord(3,compg) = inter(3,3) 

                         compg                     = compg+1
                         mesh_out % lnods(3,compt) = compg
                         mesh_out % coord(1,compg) = inter(1,4)
                         mesh_out % coord(2,compg) = inter(2,4)
                         mesh_out % coord(3,compg) = inter(3,4) 

                      endif


                   end do

                endif

             else if( pelty == PYR05 ) then
                !
                ! PYR05
                !
                 if( signp /= pnode .and. signn /= pnode ) then

                    do ii=1,2
                       sigtn=0_ip
                       sigtp=0_ip
                       compl=0_ip
                       inod1=0_ip

                       do ij=1,4
                          if(ellev(XXXXX_TO_TET04(ij,ii))>=0.0_rp) then
                             sigtp=sigtp+1
                          else 
                             sigtn=sigtn+1
                          endif
                       end do

                       if(sigtn==1) then
                          !  research of one interface triangle
                          do ij=1,4
                             if(ellev(XXXXX_TO_TET04(ij,ii))<0.0_rp) then
                                inod1=ij
                             endif
                          end do

                          do ij=1,4
                             if(ij/=inod1) then
                                compg=compg+1
                                compl=compl+1
                                if(compl==1) then
                                   compt=compt+1
                                endif
                                !
                                ! Compute the intersection of the elements with the surface 
                                !
                                l1=abs(ellev(XXXXX_TO_TET04(inod1,ii)))
                                lp=abs(ellev(XXXXX_TO_TET04(inod1,ii))-ellev(XXXXX_TO_TET04(ij,ii)))
                                x1=elcod(1,XXXXX_TO_TET04(inod1,ii))
                                x2=elcod(1,XXXXX_TO_TET04(ij,ii))
                                y1=elcod(2,XXXXX_TO_TET04(inod1,ii))
                                y2=elcod(2,XXXXX_TO_TET04(ij,ii))
                                z1=elcod(3,XXXXX_TO_TET04(inod1,ii))
                                z2=elcod(3,XXXXX_TO_TET04(ij,ii))

                                mesh_out % lnods(compl,compt) = compg
                                mesh_out % coord(1,compg)     =  x1 * (1.0_rp-l1/lp) + x2 * l1 / lp  
                                mesh_out % coord(2,compg)     =  y1 * (1.0_rp-l1/lp) + y2 * l1 / lp 
                                mesh_out % coord(3,compg)     =  z1 * (1.0_rp-l1/lp) + z2 * l1 / lp 

                             endif

                          end do

                       else if(sigtp==1) then
                          !  research of one interface triangle
                          do ij=1,4
                             if(ellev(XXXXX_TO_TET04(ij,ii))>=0.0_rp) then
                                inod1=ij
                             endif
                          end do

                          do ij=1,4
                             if(ij/=inod1) then
                                compg=compg+1
                                compl=compl+1
                                if(compl==1) then
                                   compt=compt+1
                                endif
                                !
                                ! Compute the intersection of the elements with the surface 
                                !
                                l1=abs(ellev(XXXXX_TO_TET04(inod1,ii)))
                                lp=abs(ellev(XXXXX_TO_TET04(inod1,ii))-ellev(XXXXX_TO_TET04(ij,ii)))
                                x1=elcod(1,XXXXX_TO_TET04(inod1,ii))
                                x2=elcod(1,XXXXX_TO_TET04(ij,ii))
                                y1=elcod(2,XXXXX_TO_TET04(inod1,ii))
                                y2=elcod(2,XXXXX_TO_TET04(ij,ii))
                                z1=elcod(3,XXXXX_TO_TET04(inod1,ii))
                                z2=elcod(3,XXXXX_TO_TET04(ij,ii))

                                mesh_out % lnods(compl,compt) = compg
                                mesh_out % coord(1,compg)     =  x1 * (1.0_rp-l1/lp) + x2 * l1 / lp  
                                mesh_out % coord(2,compg)     =  y1 * (1.0_rp-l1/lp) + y2 * l1 / lp 
                                mesh_out % coord(3,compg)     =  z1 * (1.0_rp-l1/lp) + z2 * l1 / lp 

                             endif
                          end do

                       else if(sigtn==2) then

                          !  research of two interface triangles
                          do ij=1,4
                             if(ellev(XXXXX_TO_TET04(ij,ii))<0.0_rp) then
                                if(inod1==0) then
                                   inod1=ij
                                else 
                                   inod2=ij
                                endif
                             endif
                          end do

                          do ij=1,4
                             if(ij/=inod1.and.ij/=inod2) then
                                !
                                ! Compute the intersection of the elements with the surface 
                                !

                                compl=compl+1
                                l1=abs(ellev(XXXXX_TO_TET04(inod1,ii)))
                                lp=abs(ellev(XXXXX_TO_TET04(inod1,ii))-ellev(XXXXX_TO_TET04(ij,ii)))
                                x1=elcod(1,XXXXX_TO_TET04(inod1,ii))
                                x2=elcod(1,XXXXX_TO_TET04(ij,ii))
                                y1=elcod(2,XXXXX_TO_TET04(inod1,ii))
                                y2=elcod(2,XXXXX_TO_TET04(ij,ii))
                                z1=elcod(3,XXXXX_TO_TET04(inod1,ii))
                                z2=elcod(3,XXXXX_TO_TET04(ij,ii))

                                inter(1,compl)=  x1*(1-l1/lp)+x2*l1/lp  
                                inter(2,compl)=  y1*(1-l1/lp)+y2*l1/lp 
                                inter(3,compl)=  z1*(1-l1/lp)+z2*l1/lp 


                                compl=compl+1
                                l1=abs(ellev(XXXXX_TO_TET04(inod2,ii)))
                                lp=abs(ellev(XXXXX_TO_TET04(inod2,ii))-ellev(XXXXX_TO_TET04(ij,ii)))
                                x1=elcod(1,XXXXX_TO_TET04(inod2,ii))
                                x2=elcod(1,XXXXX_TO_TET04(ij,ii))
                                y1=elcod(2,XXXXX_TO_TET04(inod2,ii))
                                y2=elcod(2,XXXXX_TO_TET04(ij,ii))
                                z1=elcod(3,XXXXX_TO_TET04(inod2,ii))
                                z2=elcod(3,XXXXX_TO_TET04(ij,ii))

                                inter(1,compl)=  x1*(1-l1/lp)+x2*l1/lp  
                                inter(2,compl)=  y1*(1-l1/lp)+y2*l1/lp 
                                inter(3,compl)=  z1*(1-l1/lp)+z2*l1/lp 

                             endif

                          end do

                          compg                     = compg+1
                          compt                     = compt+1

                          mesh_out % lnods(1,compt) = compg
                          mesh_out % coord(1,compg) = inter(1,1)
                          mesh_out % coord(2,compg) = inter(2,1)
                          mesh_out % coord(3,compg) = inter(3,1) 

                          compg                     = compg+1
                          mesh_out % lnods(2,compt) = compg
                          mesh_out % coord(1,compg) = inter(1,2)
                          mesh_out % coord(2,compg) = inter(2,2)
                          mesh_out % coord(3,compg) = inter(3,2) 

                          compg                     = compg+1
                          mesh_out % lnods(3,compt) = compg
                          mesh_out % coord(1,compg) = inter(1,3)
                          mesh_out % coord(2,compg) = inter(2,3)
                          mesh_out % coord(3,compg) = inter(3,3) 


                          compg                     = compg+1
                          compt                     = compt+1
                          mesh_out % lnods(1,compt) = compg
                          mesh_out % coord(1,compg) = inter(1,2)
                          mesh_out % coord(2,compg) = inter(2,2)
                          mesh_out % coord(3,compg) = inter(3,2) 

                          compg                     = compg+1
                          mesh_out % lnods(2,compt) = compg
                          mesh_out % coord(1,compg) = inter(1,3)
                          mesh_out % coord(2,compg) = inter(2,3)
                          mesh_out % coord(3,compg) = inter(3,3) 

                          compg                     = compg+1
                          mesh_out % lnods(3,compt) = compg
                          mesh_out % coord(1,compg) = inter(1,4)
                          mesh_out % coord(2,compg) = inter(2,4)
                          mesh_out % coord(3,compg) = inter(3,4) 


                       endif


                    end do

                 endif

              else

                 call runend('CUT_LEVEL: ELEMENT YPE NOT CODED')
                 
              end if
          end if
       end if
    end do
    !
    ! Give element types
    !
    if( pdime == 2 ) then
       do ielem = 1,mesh_out % nelem
          mesh_out % ltype(ielem) = BAR02 !BAR3D
       end do
    else
       do ielem = 1,mesh_out % nelem
          mesh_out % ltype(ielem) = TRI03 !SHELL
       end do       
    end if
    !
    ! Identity permutation
    !
    do ipoin = 1,mesh_out % npoin
       mesh_out % permn(ipoin) = ipoin
    end do
    do ielem = 1,mesh_out % nelem
       mesh_out % perme(ielem) = ielem
    end do
    call memory_deallo(memor_loc,'IFCUT',vacal,ifcut)
    
    if( present(MEMORY_COUNTER) )  MEMORY_COUNTER = memor_loc

  end subroutine cut_level

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-03
  !> @brief   Compute elements bounding boxes
  !> @details Compute elements bounding boxes
  !> 
  !-----------------------------------------------------------------------
  
  subroutine element_bb(mesh,bobox,MEMORY_COUNTER,ONLY_DEALLOCATE)

    class(mesh_type_basic),          intent(in)    :: mesh
    real(rp),               pointer, intent(inout) :: bobox(:,:,:)   
    integer(8),   optional,          intent(inout) :: MEMORY_COUNTER(2)
    logical(lg),  optional,          intent(in)    :: ONLY_DEALLOCATE
    integer(ip)                                    :: pelty,ielem,inode,ipoin
    real(rp)                                       :: comin(3),comax(3)
    real(rp)                                       :: element_distance
    integer(8)                                     :: memor_loc(2)
    real(rp)                                       :: elcod(mesh % ndime,mesh % mnode)

    memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)

    if( optional_argument(.false.,ONLY_DEALLOCATE) ) then
       call memory_deallo(memor_loc,'BOBOX','element_bb',bobox)
    else
       if( .not. associated(bobox) ) & 
            call memory_alloca(memor_loc,'BOBOX','element_bb',bobox,2_ip,mesh % ndime,mesh % nelem)

       select type ( mesh )

       class is ( mesh_type_basic )

          do ielem = 1,mesh % nelem
             pelty =  mesh % ltype(ielem)
             comin =  huge(1.0_rp)*0.1_rp
             comax = -huge(1.0_rp)*0.1_rp
             do inode = 1,element_type(pelty) % number_nodes
                ipoin                 = mesh % lnods(inode,ielem)
                comin(1:mesh % ndime) = min(comin(1:mesh % ndime),mesh % coord(1:mesh % ndime,ipoin))
                comax(1:mesh % ndime) = max(comax(1:mesh % ndime),mesh % coord(1:mesh % ndime,ipoin))
             end do
             bobox(1,1:mesh % ndime,ielem) = comin(1:mesh % ndime)
             bobox(2,1:mesh % ndime,ielem) = comax(1:mesh % ndime)
          end do

       class is ( bmsh_type_basic )

          do ielem = 1,mesh % nelem
             pelty =  mesh % ltype(ielem)
             comin =  huge(1.0_rp)*0.1_rp
             comax = -huge(1.0_rp)*0.1_rp
             do inode = 1,element_type(pelty) % number_nodes
                ipoin                       = mesh % lnods(inode,ielem)
                elcod(1:mesh % ndime,inode) = mesh % coord(1:mesh % ndime,ipoin)
                comin(1:mesh % ndime)       = min(comin(1:mesh % ndime),mesh % coord(1:mesh % ndime,ipoin))
                comax(1:mesh % ndime)       = max(comax(1:mesh % ndime),mesh % coord(1:mesh % ndime,ipoin))
             end do
             call elmgeo_element_distance(mesh % ndime,mesh % ltype(ielem),elcod,element_distance)
             bobox(1,1:mesh % ndime,ielem) = comin(1:mesh % ndime) - 1.0e-16_rp * element_distance - epsil
             bobox(2,1:mesh % ndime,ielem) = comax(1:mesh % ndime) + 1.0e-16_rp * element_distance + epsil          
          end do

       end select

    end if

    if( present(MEMORY_COUNTER) )  MEMORY_COUNTER = memor_loc

  end subroutine element_bb
    
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-03
  !> @brief   Compute boundary bounding boxes
  !> @details Compute boundary bounding boxes. Check for possible
  !>          zero bounding boxes
  !> 
  !-----------------------------------------------------------------------
  
  subroutine boundary_bb(mesh,bobox,MEMORY_COUNTER,ONLY_DEALLOCATE)

    class(mesh_type_basic),          intent(in)    :: mesh
    real(rp),               pointer, intent(inout) :: bobox(:,:,:)   
    integer(8),   optional,          intent(inout) :: MEMORY_COUNTER(2)
    logical(lg),  optional,          intent(in)    :: ONLY_DEALLOCATE
    integer(ip)                                    :: pelty,ielem,inode,ipoin
    real(rp)                                       :: comin(3),comax(3),rdime
    real(rp)                                       :: element_distance
    real(rp)                                       :: elcod(mesh % ndime,mesh % mnode)
    integer(8)                                     :: memor_loc(2)

    memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)

    if( optional_argument(.false.,ONLY_DEALLOCATE) ) then
       call memory_deallo(memor_loc,'BOBOX','element_bb',bobox)
    else
       rdime     = 1.0_rp / real(mesh % ndime,rp)

       if( associated(mesh % boundary) ) then
          if( .not. associated(bobox) ) & 
               call memory_alloca(memor_loc,'BOBOX','boundary_bb',bobox,2_ip,mesh % boundary % ndime,mesh % boundary % nelem)

          do ielem = 1,mesh % boundary % nelem
             pelty =  mesh % boundary % ltype(ielem)
             comin =  huge(1.0_rp)*0.1_rp
             comax = -huge(1.0_rp)*0.1_rp
             do inode = 1,element_type(pelty) % number_nodes
                ipoin                       = mesh % boundary % lnods(inode,ielem)
                elcod(1:mesh % ndime,inode) = mesh % boundary % coord(1:mesh % ndime,ipoin)
                comin(1:mesh % ndime)       = min(comin(1:mesh % ndime),mesh % boundary % coord(1:mesh % ndime,ipoin))
                comax(1:mesh % ndime)       = max(comax(1:mesh % ndime),mesh % boundary % coord(1:mesh % ndime,ipoin))
             end do
             call elmgeo_element_distance(mesh % ndime,mesh % boundary % ltype(ielem),elcod,element_distance)
             bobox(1,1:mesh % ndime,ielem) = comin(1:mesh % ndime) - 1.0e-16_rp * element_distance - epsil
             bobox(2,1:mesh % ndime,ielem) = comax(1:mesh % ndime) + 1.0e-16_rp * element_distance + epsil
          end do
       else
          call runend('DEF_KINTYP_MESH_BASIC: NO ASSOCIATED NOUDANRY TO MESH')       
       end if
    end if

    if( present(MEMORY_COUNTER) )  MEMORY_COUNTER = memor_loc

  end subroutine boundary_bb
    
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-08
  !> @brief   Create a circular mesh
  !> @details Create a circular mesh in 3D with QUA04 elements
  !> 
  !-----------------------------------------------------------------------

  subroutine mesh_3Dcircle(mesh,center,normal,radius,nelem_rad,nelem_cir,MEMORY_COUNTER)

    class(mesh_type_basic),         intent(inout) :: mesh
    real(rp),                       intent(in)    :: center(:)
    real(rp),                       intent(in)    :: normal(:)
    real(rp),                       intent(in)    :: radius
    integer(ip),                    intent(in)    :: nelem_rad
    integer(ip),                    intent(in)    :: nelem_cir
    integer(8),  optional,          intent(inout) :: MEMORY_COUNTER(2)
    integer(ip)                                   :: ielem_cir,ielem_rad
    integer(ip)                                   :: kpoin,kelem
    real(rp)                                      :: baloc(3,3)
    real(rp)                                      :: theta,v(3)
    real(rp)                                      :: rotma(3,3)
    
    mesh % ndime = 3
    mesh % nelem = nelem_cir*nelem_rad
    mesh % npoin = mesh % nelem + 1
    mesh % mnode = 4
    
    call mesh % alloca(MEMORY_COUNTER=MEMORY_COUNTER)
    
    baloc(1:3,1) = normal(1:3)
    call maths_normalize_vector       (3_ip,baloc(:,1))
    call maths_local_orthonormal_basis(3_ip,baloc(:,1))

    kpoin = 1
    mesh % coord(1:3,kpoin) = center(1:3)
    do ielem_cir = 1,nelem_cir
       theta = real(ielem_cir-1,rp)/real(nelem_cir,rp)*2.0_rp*pi
       do ielem_rad = 1,nelem_rad
          kpoin                   = kpoin + 1
          rotma                   = maths_rotation_matrix_3D(theta,V=baloc(1:3,1))
          v                       = matmul(rotma,baloc(1:3,2))
          mesh % coord(1:3,kpoin) = center(1:3) + radius*v(1:3)*real(ielem_rad,rp)/real(nelem_rad,rp)
       end do
    end do

    kelem = 0
    !
    ! TRI03 elements
    !
    do ielem_cir = 1,nelem_cir-1
       kelem                    = kelem + 1
       mesh % ltype(kelem)      = TRI03
       mesh % lnods(1,kelem)    = 1
       mesh % lnods(2,kelem)    = 1 + (ielem_cir-1)*nelem_rad+1
       mesh % lnods(3,kelem)    = 1 + (ielem_cir  )*nelem_rad+1
    end do
    kelem                       = kelem + 1
    mesh % ltype(kelem)         = TRI03
    ielem_cir                   = nelem_cir-1
    mesh % lnods(1,kelem)       = 1
    mesh % lnods(2,kelem)       = 1 + ielem_cir*nelem_rad+1
    ielem_cir                   = 1
    mesh % lnods(3,kelem)       = 1 + (ielem_cir-1)*nelem_rad+1
    !
    ! QUA04 elements
    !
    do ielem_cir = 1,nelem_cir-1
       do ielem_rad = 1,nelem_rad-1
          kelem                 = kelem + 1
          mesh % ltype(kelem)   = QUA04
          mesh % lnods(1,kelem) = 1 + (ielem_cir-1)*nelem_rad+ielem_rad
          mesh % lnods(2,kelem) = 1 + (ielem_cir-1)*nelem_rad+ielem_rad+1
          mesh % lnods(3,kelem) = 1 + (ielem_cir-1)*nelem_rad+ielem_rad+nelem_rad+1
          mesh % lnods(4,kelem) = 1 + (ielem_cir-1)*nelem_rad+ielem_rad+nelem_rad
       end do
    end do    
    do ielem_rad = 1,nelem_rad-1
       kelem                    = kelem + 1       
       ielem_cir                = nelem_cir-1
       mesh % ltype(kelem)      = QUA04
       mesh % lnods(2,kelem)    = 1 + (ielem_cir-1)*nelem_rad+ielem_rad+nelem_rad+1 
       mesh % lnods(1,kelem)    = 1 + (ielem_cir-1)*nelem_rad+ielem_rad+nelem_rad       
       ielem_cir                = 1
       mesh % lnods(4,kelem)    = 1 + (ielem_cir-1)*nelem_rad+ielem_rad
       mesh % lnods(3,kelem)    = 1 + (ielem_cir-1)*nelem_rad+ielem_rad+1
    end do
   
  end subroutine mesh_3Dcircle
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-08
  !> @brief   Create a ring mesh
  !> @details Create a ring mesh in 3D with POI3D elements
  !>          The final mesh will count with npoin_cir*npoin_rad + 1
  !>          node as a center node is added
  !> 
  !-----------------------------------------------------------------------

  subroutine mesh_3Dring(mesh,center,normal,radius,npoin_rad,npoin_cir,MEMORY_COUNTER)

    class(mesh_type_basic),         intent(inout) :: mesh
    real(rp),                       intent(in)    :: center(:)
    real(rp),                       intent(in)    :: normal(:)
    real(rp),                       intent(in)    :: radius
    integer(ip),                    intent(in)    :: npoin_rad
    integer(ip),                    intent(in)    :: npoin_cir
    integer(8),  optional,          intent(inout) :: MEMORY_COUNTER(2)
    integer(ip)                                   :: ipoin_cir,ipoin_rad
    integer(ip)                                   :: kpoin,ielem
    real(rp)                                      :: baloc(3,3)
    real(rp)                                      :: theta,v(3)
    real(rp)                                      :: rotma(3,3)
    
    mesh % ndime = 3
    mesh % npoin = npoin_cir*npoin_rad + 1 
    mesh % nelem = npoin_cir*npoin_rad + 1
    mesh % mnode = 1
    
    call mesh % alloca(MEMORY_COUNTER=MEMORY_COUNTER)
    
    baloc(1:3,1) = normal(1:3)
    call maths_normalize_vector       (3_ip,baloc(:,1))
    call maths_local_orthonormal_basis(3_ip,baloc(:,1))

    kpoin = 1
    mesh % coord(1:3,kpoin) = center(1:3)
    do ipoin_cir = 1,npoin_cir
       theta = real(ipoin_cir-1,rp)/real(npoin_cir,rp)*2.0_rp*pi
       do ipoin_rad = 1,npoin_rad
          kpoin                   = kpoin + 1
          rotma                   = maths_rotation_matrix_3D(theta,V=baloc(1:3,1))
          v                       = matmul(rotma,baloc(1:3,2))
          mesh % coord(1:3,kpoin) = center(1:3) + radius*v(1:3)*real(ipoin_rad,rp)/real(npoin_rad,rp)
       end do
    end do
    !
    ! POI3D elements
    !
    do ielem = 1,mesh % nelem
       mesh % ltype(ielem)      = POI3D
       mesh % lnods(1,ielem)    = ielem
    end do
    
  end subroutine mesh_3Dring
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-08
  !> @brief   Create a box mesh
  !> @details Create a box mesh in 3D with POI3D elements
  !>          The final mesh will count with npoin_cir*npoin_rad + 1
  !>          node as a center node is added
  !> 
  !-----------------------------------------------------------------------

  subroutine mesh_3Dbox(mesh,ndime,boxes,comin,comax,MEMORY_COUNTER)

    class(mesh_type_basic), intent(inout)  :: mesh              !< Output mesh
    integer(ip),            intent(in)     :: ndime             !< Dimension
    integer(ip),            intent(in)     :: boxes(:)          !< # elements
    real(rp),               intent(in)     :: comin(:)          !< Min coordinates
    real(rp),               intent(in)     :: comax(:)          !< Max coordinates
    integer(8),  optional,  intent(inout)  :: MEMORY_COUNTER(2) !< Memory counter
    integer(ip)                            :: ii,jj,kk
    integer(ip)                            :: ipoin
    integer(ip)                            :: ielem
    real(rp)                               :: boxer(3)
    real(rp)                               :: delta(3)
    real(rp)                               :: yy,zz

    mesh % ndime = ndime
    mesh % mnode = 1
    if( mesh % ndime == 2 ) then
       mesh % nelem = (boxes(1)+1)*(boxes(2)+1)
       mesh % npoin = (boxes(1)+1)*(boxes(2)+1)
    else if( mesh % ndime == 3 ) then
       mesh % nelem = (boxes(1)+1)*(boxes(2)+1)*(boxes(3)+1)
       mesh % npoin = (boxes(1)+1)*(boxes(2)+1)*(boxes(3)+1)
    else
       return
    end if
    delta(1:ndime) = comax(1:ndime)-comin(1:ndime)
    boxer(1:ndime) = real(boxes(1:ndime),rp)
    
    call mesh % alloca(MEMORY_COUNTER=MEMORY_COUNTER)
    ipoin = 0
    
    if( ndime == 2 ) then
       do jj = 1,boxes(2)+1
          yy = comin(2) + real(jj-1,rp)/boxer(2) * delta(2)
          do ii = 1,boxes(1)+1
             ipoin = ipoin + 1
             mesh % coord(1,ipoin) = comin(1) + real(ii-1,rp)/boxer(1) * delta(1)
             mesh % coord(2,ipoin) = yy
          end do
       end do
    else
       do kk = 1,boxes(3)+1
          zz = comin(3) + real(kk-1,rp)/boxer(3) * delta(3)
          do jj = 1,boxes(2)+1
             yy = comin(2) + real(jj-1,rp)/boxer(2) * delta(2)
             do ii = 1,boxes(1)+1
                ipoin = ipoin + 1
                mesh % coord(1,ipoin) = comin(1) + real(ii-1,rp)/boxer(1) * delta(1)
                mesh % coord(2,ipoin) = yy
                mesh % coord(3,ipoin) = zz
             end do
          end do
       end do
    end if    
    !
    ! POI3D elements
    !
    do ielem = 1,mesh % nelem
       mesh % ltype(ielem)      = POI3D
       mesh % lnods(1,ielem)    = ielem
    end do
    
  end subroutine mesh_3Dbox
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-08
  !> @brief   Create a 3D mesh from a BAR3D mesh
  !> @details Create a 3D mesh from a BAR3D mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine mesh_from_BAR3D(mesh_out,mesh_in,nelem_rad,nelem_cir,areas,perm,MEMORY_COUNTER)
    
    class(mesh_type_basic),         intent(inout) :: mesh_out
    class(mesh_type_basic),         intent(in)    :: mesh_in
    integer(ip),                    intent(in)    :: nelem_rad
    integer(ip),                    intent(in)    :: nelem_cir
    real(rp),              pointer, intent(in)    :: areas(:)
    integer(ip), optional, pointer, intent(inout) :: perm(:)
    integer(8),  optional,          intent(inout) :: MEMORY_COUNTER(2)
    integer(ip)                                   :: ielem,ipoin,kdime,kpoin,jelem
    integer(ip)                                   :: ipoin1,ipoin2,inode,num_bar3d
    integer(ip)                                   :: jpoin(2)
    integer(ip)                                   :: ielem_rad,ielem_cir,kelem
    integer(ip)                                   :: kpoin_ini(2)
    integer(8)                                    :: memor_loc(2)
    real(rp)                                      :: r1,r2,baloc(3,3),c(3,2)
    real(rp)                                      :: theta,rotma(3,3),v(3)

    if( present(MEMORY_COUNTER) ) then
       memor_loc = MEMORY_COUNTER
    else
       memor_loc = 0_8
    end if

    kdime     = mesh_in % ndime
    num_bar3d = 0
    do ielem = 1,mesh_in % nelem
       if( mesh_in % ltype(ielem) == BAR3D ) num_bar3d = num_bar3d + 1
    end do
    mesh_out % ndime = 3
    mesh_out % nelem = num_bar3d * nelem_rad * nelem_cir
    mesh_out % npoin = 2*num_bar3d*(nelem_rad * nelem_cir + 1)
    mesh_out % mnode = 8
    call mesh_out % alloca(MEMORY_COUNTER=MEMORY_COUNTER)
    if( present(perm) ) then
       call memory_alloca(memor_loc,'PERM','mesh_type_basic_valid_mesh',perm,mesh_out % npoin)       
    end if
    
    kpoin = 0
    kelem = 0
    do ielem = 1,mesh_in % nelem
       if( mesh_in % ltype(ielem) == BAR3D ) then
          ipoin1       = mesh_in % lnods(1,ielem)
          ipoin2       = mesh_in % lnods(2,ielem)
          r1           = sqrt(areas(ipoin1)/pi)
          r2           = sqrt(areas(ipoin2)/pi)
          c(1:3,1)     = mesh_in % coord(1:kdime,ipoin1)
          c(1:3,2)     = mesh_in % coord(1:kdime,ipoin2)
          baloc(1:3,1) = c(1:3,2)-c(1:3,1)
          call maths_normalize_vector(3_ip,baloc(:,1))
          call maths_local_orthonormal_basis(3_ip,baloc(:,1))
          do inode = 1,2
             kpoin = kpoin + 1
             kpoin_ini(inode) = kpoin
             mesh_out % coord(1:3,kpoin) = c(1:3,inode)
             if( present(perm) ) perm(kpoin) = mesh_in % lnods(inode,ielem)
             do ielem_cir = 1,nelem_cir
                theta = real(ielem_cir-1,rp)/real(nelem_cir,rp)*2.0_rp*pi
                do ielem_rad = 1,nelem_rad
                   kpoin = kpoin + 1
                   rotma = maths_rotation_matrix_3D(theta,V=baloc(1:3,1))
                   v      = matmul(rotma,baloc(1:3,2))
                   mesh_out % coord(1:3,kpoin) = c(1:3,inode) + r1*v(1:3)*real(ielem_rad,rp)/real(nelem_rad,rp)
                   if( present(perm) ) perm(kpoin) = mesh_in % lnods(inode,ielem)
                end do
             end do
          end do
          !
          ! PEN06 around centerline
          !
          jpoin(1) = kpoin_ini(1)
          jpoin(2) = kpoin_ini(2)
          do jelem = 1,nelem_cir-1
             kelem                     = kelem + 1
             mesh_out % ltype(kelem)   = PEN06
             mesh_out % lnods(1,kelem) = kpoin_ini(1)
             mesh_out % lnods(2,kelem) = kpoin_ini(1) + (jelem-1)*nelem_rad+1
             mesh_out % lnods(3,kelem) = kpoin_ini(1) + (jelem  )*nelem_rad+1
             mesh_out % lnods(4,kelem) = kpoin_ini(2)
             mesh_out % lnods(5,kelem) = kpoin_ini(2) + (jelem-1)*nelem_rad+1
             mesh_out % lnods(6,kelem) = kpoin_ini(2) + (jelem)  *nelem_rad+1
          end do
          kelem                     = kelem + 1
          mesh_out % ltype(kelem)   = PEN06
          jelem                     = nelem_cir-1
          mesh_out % lnods(1,kelem) = kpoin_ini(1)
          mesh_out % lnods(2,kelem) = kpoin_ini(1) + (jelem)*nelem_rad+1
          jelem                     = 1
          mesh_out % lnods(3,kelem) = kpoin_ini(1) + (jelem-1)*nelem_rad+1
          jelem                     = nelem_cir-1
          mesh_out % lnods(4,kelem) = kpoin_ini(2)
          mesh_out % lnods(5,kelem) = kpoin_ini(2) + (jelem)*nelem_rad+1
          jelem                     = 1
          mesh_out % lnods(6,kelem) = kpoin_ini(2) + (jelem-1)*nelem_rad+1
          !
          ! HEX08
          !
          do ielem_cir = 1,nelem_cir-1
             do ielem_rad = 1,nelem_rad-1
                kelem                     = kelem + 1
                mesh_out % ltype(kelem)   = HEX08
                mesh_out % lnods(1,kelem) = kpoin_ini(1) + (ielem_cir-1)*nelem_rad+ielem_rad
                mesh_out % lnods(2,kelem) = kpoin_ini(1) + (ielem_cir-1)*nelem_rad+ielem_rad+1
                mesh_out % lnods(3,kelem) = kpoin_ini(1) + (ielem_cir-1)*nelem_rad+ielem_rad+nelem_rad+1
                mesh_out % lnods(4,kelem) = kpoin_ini(1) + (ielem_cir-1)*nelem_rad+ielem_rad+nelem_rad
                mesh_out % lnods(5,kelem) = kpoin_ini(2) + (ielem_cir-1)*nelem_rad+ielem_rad
                mesh_out % lnods(6,kelem) = kpoin_ini(2) + (ielem_cir-1)*nelem_rad+ielem_rad+1
                mesh_out % lnods(7,kelem) = kpoin_ini(2) + (ielem_cir-1)*nelem_rad+ielem_rad+nelem_rad+1
                mesh_out % lnods(8,kelem) = kpoin_ini(2) + (ielem_cir-1)*nelem_rad+ielem_rad+nelem_rad
             end do
          end do
          do ielem_rad = 1,nelem_rad-1
             kelem                     = kelem + 1

             ielem_cir                 = nelem_cir-1
             mesh_out % ltype(kelem)   = HEX08
             mesh_out % lnods(2,kelem) = kpoin_ini(1) + (ielem_cir-1)*nelem_rad+ielem_rad+nelem_rad+1 
             mesh_out % lnods(1,kelem) = kpoin_ini(1) + (ielem_cir-1)*nelem_rad+ielem_rad+nelem_rad
             
             ielem_cir                 = 1
             mesh_out % lnods(4,kelem) = kpoin_ini(1) + (ielem_cir-1)*nelem_rad+ielem_rad
             mesh_out % lnods(3,kelem) = kpoin_ini(1) + (ielem_cir-1)*nelem_rad+ielem_rad+1
             
             ielem_cir                 = nelem_cir-1
             mesh_out % lnods(6,kelem) = kpoin_ini(2) + (ielem_cir-1)*nelem_rad+ielem_rad+nelem_rad+1
             mesh_out % lnods(5,kelem) = kpoin_ini(2) + (ielem_cir-1)*nelem_rad+ielem_rad+nelem_rad
             
             ielem_cir                 = 1
             mesh_out % lnods(8,kelem) = kpoin_ini(2) + (ielem_cir-1)*nelem_rad+ielem_rad
             mesh_out % lnods(7,kelem) = kpoin_ini(2) + (ielem_cir-1)*nelem_rad+ielem_rad+1
          end do          
       end if

    end do
    do ielem = 1,mesh_out % nelem
       mesh_out % perme(ielem) = ielem
    end do
    do ipoin = 1,mesh_out % npoin
       mesh_out % permn(ipoin) = ipoin
    end do
    
    if( present(MEMORY_COUNTER) )  MEMORY_COUNTER = memor_loc

  end subroutine mesh_from_BAR3D

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-03
  !> @brief   Associate a mesh
  !> @details Associate a mesh
  !> 
  !-----------------------------------------------------------------------
  
  subroutine assoc(mesh,ndime,mnode,nelem,npoin,lnods,ltype,coord,leinv_loc,lninv_loc,lelbo,lboel,BOUNDARY_MESH)

    class(mesh_type_basic),           target,  intent(inout) :: mesh
    integer(ip),                               intent(in)    :: ndime
    integer(ip),                               intent(in)    :: mnode
    integer(ip),                               intent(in)    :: nelem
    integer(ip),                               intent(in)    :: npoin
    integer(ip),                      pointer, intent(in)    :: lnods(:,:)
    integer(ip),                      pointer, intent(in)    :: ltype(:)
    real(rp),               optional, pointer, intent(in)    :: coord(:,:)
    integer(ip),            optional, pointer, intent(in)    :: leinv_loc(:)
    integer(ip),            optional, pointer, intent(in)    :: lninv_loc(:)
    integer(ip),            optional, pointer, intent(in)    :: lelbo(:)
    integer(ip),            optional, pointer, intent(in)    :: lboel(:,:)
    logical(lg),            optional,          intent(in)    :: BOUNDARY_MESH
    logical(lg)                                              :: if_boundary

    if_boundary = optional_argument(.false.,BOUNDARY_MESH)
    
    select type ( mesh )

    class is ( bmsh_type_basic )

       mesh % ndime =  ndime
       mesh % mnode =  mnode
       mesh % nelem =  nelem
       mesh % npoin =  npoin
       mesh % ltype => ltype
       mesh % lnods => lnods
       if( present(coord)     ) then
          mesh % coord     => coord
       else if( associated(mesh % mesh) ) then
          mesh % coord     => mesh % mesh % coord
       end if
       if( present(leinv_loc) ) mesh % leinv_loc => leinv_loc
       if( present(lninv_loc) ) mesh % lninv_loc => lninv_loc
       if( present(lelbo)     ) mesh % lelbo     => lelbo
       if( present(lboel)     ) mesh % lboel     => lboel 

    class is ( mesh_type_basic )

       if( if_boundary ) then
          if( .not. associated(mesh % boundary) ) allocate(mesh % boundary)
          
          mesh % boundary % ndime =  ndime
          mesh % boundary % mnode =  mnode
          mesh % boundary % nelem =  nelem
          mesh % boundary % npoin =  npoin
          mesh % boundary % ltype => ltype
          mesh % boundary % lnods => lnods
          mesh % boundary % mesh  => mesh
          if( present(coord)     ) then
             mesh % boundary % coord => coord
          else
             mesh % boundary % coord => mesh % coord
          end if
          if( present(leinv_loc) ) mesh % boundary % leinv_loc => leinv_loc
          if( present(lninv_loc) ) mesh % boundary % lninv_loc => lninv_loc
          if( present(lelbo)     ) mesh % boundary % lelbo     => lelbo
          if( present(lboel)     ) mesh % boundary % lboel     => lboel
          
       else

          mesh % ndime =  ndime
          mesh % mnode =  mnode
          mesh % nelem =  nelem
          mesh % npoin =  npoin
          mesh % ltype => ltype
          mesh % lnods => lnods
          if( present(coord)     ) mesh % coord     => coord
          if( present(leinv_loc) ) mesh % leinv_loc => leinv_loc
          if( present(lninv_loc) ) mesh % lninv_loc => lninv_loc

       end if
       
    end select

  end subroutine assoc

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-03
  !> @brief   Disassociate a mesh
  !> @details Disassociate a mesh
  !> 
  !-----------------------------------------------------------------------
  
  subroutine disassoc(mesh,BOUNDARY_MESH)
    
    class(mesh_type_basic),                    intent(inout) :: mesh
    logical(lg),            optional,          intent(in)    :: BOUNDARY_MESH
    logical(lg)                                              :: if_boundary

    if_boundary = optional_argument(.false.,BOUNDARY_MESH)

    select type ( mesh )
       
    class is ( bmsh_type_basic )
       
       nullify(mesh % lnods)
       nullify(mesh % ltype)
       nullify(mesh % leinv_loc)
       nullify(mesh % lninv_loc)
       nullify(mesh % coord)
       nullify(mesh % permn)
       nullify(mesh % perme)      
       nullify(mesh % tags )      
       nullify(mesh % lelbo)
       nullify(mesh % lboel)
       nullify(mesh % boundary)
       
    class is ( mesh_type_basic )

       if( if_boundary ) then
          
          nullify(mesh % boundary % lnods)
          nullify(mesh % boundary % ltype)
          nullify(mesh % boundary % leinv_loc)
          nullify(mesh % boundary % lninv_loc)
          nullify(mesh % boundary % coord)
          nullify(mesh % boundary % permn)
          nullify(mesh % boundary % perme)      
          nullify(mesh % boundary % tags )      
          nullify(mesh % boundary % lelbo)
          nullify(mesh % boundary % lboel)
          nullify(mesh % boundary % mesh )

       else

          nullify(mesh % lnods)
          nullify(mesh % ltype)
          nullify(mesh % leinv_loc)
          nullify(mesh % lninv_loc)
          nullify(mesh % coord)
          nullify(mesh % permn)
          nullify(mesh % perme) 
          nullify(mesh % tags ) 
          nullify(mesh % boundary)
          
       end if
       
    end select
  
  end subroutine disassoc

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-03
  !> @brief   Mesh volume
  !> @details Compute the volume of a mesh
  !> 
  !-----------------------------------------------------------------------

  function volume(mesh) result(vol)

    class(mesh_type_basic), intent(in)  :: mesh              !< Mesh
    real(rp)                            :: vol               !< Mesh volume
    integer(ip)                         :: pelty,ielem,pnode
    real(rp),               allocatable :: elcod(:,:)
    real(rp)                            :: element_volume

    vol = 0.0_rp
    allocate(elcod(mesh % ndime,mesh % mnode))
    do ielem = 1,mesh % nelem
       pelty = mesh % ltype(ielem)
       pnode = element_type(pelty) % number_nodes
       elcod(1:mesh % ndime,1:pnode) = mesh % coord(1:mesh % ndime,mesh % lnods(1:pnode,ielem))
       call elmgeo_element_volume(mesh % ndime,pelty,elcod,element_volume)
       vol = vol + element_volume
    end do
    deallocate(elcod)
    
  end function volume
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-03
  !> @brief   Bounding box mesh
  !> @details Create a mesh from bounding boxes
  !> 
  !-----------------------------------------------------------------------

  subroutine mesh_bb(mesh,bobox,MEMORY_COUNTER)

    class(mesh_type_basic),                   intent(inout) :: mesh              !< Mesh
    real(rp),                        pointer, intent(in)    :: bobox(:,:,:)      !< Bounding boxes
    integer(8),            optional,          intent(inout) :: MEMORY_COUNTER(2) !< Memory counter
    integer(ip)                                             :: ielem,ipoin
    integer(ip)                                             :: inode,pelty

    mesh % nelem = memory_size(bobox,3_ip)
    mesh % ndime = memory_size(bobox,2_ip)

    if( mesh % ndime == 2 ) then
       mesh % mnode = 4
       pelty        = QUA04
    else
       mesh % mnode = 8
       pelty        = HEX08
    end if
    mesh % npoin = mesh % mnode * mesh % nelem
    call mesh % alloca(MEMORY_COUNTER=MEMORY_COUNTER)

    ipoin = 0

    if( mesh % ndime == 2 ) then
       do ielem = 1,mesh % nelem
          mesh % ltype(ielem) = pelty
          do inode = 1,mesh % mnode
             mesh % lnods(inode,ielem) = ipoin + inode
          end do
          ipoin                 = ipoin + 1
          mesh % coord(1,ipoin) = bobox(1,1,ielem)
          mesh % coord(2,ipoin) = bobox(1,2,ielem)
          ipoin                 = ipoin + 1
          mesh % coord(1,ipoin) = bobox(2,1,ielem)
          mesh % coord(2,ipoin) = bobox(1,2,ielem)
          ipoin                 = ipoin + 1
          mesh % coord(1,ipoin) = bobox(2,1,ielem)
          mesh % coord(2,ipoin) = bobox(2,2,ielem)
          ipoin                 = ipoin + 1
          mesh % coord(1,ipoin) = bobox(1,1,ielem)
          mesh % coord(2,ipoin) = bobox(2,2,ielem)
       end do
    else
       do ielem = 1,mesh % nelem
          mesh % ltype(ielem) = pelty
          do inode = 1,mesh % mnode
             mesh % lnods(inode,ielem) = ipoin + inode
          end do
          ipoin                 = ipoin + 1
          mesh % coord(1,ipoin) = bobox(1,1,ielem)
          mesh % coord(2,ipoin) = bobox(1,2,ielem)
          mesh % coord(3,ipoin) = bobox(1,3,ielem)
          ipoin                 = ipoin + 1
          mesh % coord(1,ipoin) = bobox(2,1,ielem)
          mesh % coord(2,ipoin) = bobox(1,2,ielem)
          mesh % coord(3,ipoin) = bobox(1,3,ielem)
          ipoin                 = ipoin + 1
          mesh % coord(1,ipoin) = bobox(2,1,ielem)
          mesh % coord(2,ipoin) = bobox(2,2,ielem)
          mesh % coord(3,ipoin) = bobox(1,3,ielem)
          ipoin                 = ipoin + 1
          mesh % coord(1,ipoin) = bobox(1,1,ielem)
          mesh % coord(2,ipoin) = bobox(2,2,ielem)
          mesh % coord(3,ipoin) = bobox(1,3,ielem)

          ipoin                 = ipoin + 1
          mesh % coord(1,ipoin) = bobox(1,1,ielem)
          mesh % coord(2,ipoin) = bobox(1,2,ielem)
          mesh % coord(3,ipoin) = bobox(2,3,ielem)
          ipoin                 = ipoin + 1
          mesh % coord(1,ipoin) = bobox(2,1,ielem)
          mesh % coord(2,ipoin) = bobox(1,2,ielem)
          mesh % coord(3,ipoin) = bobox(2,3,ielem)
          ipoin                 = ipoin + 1
          mesh % coord(1,ipoin) = bobox(2,1,ielem)
          mesh % coord(2,ipoin) = bobox(2,2,ielem)
          mesh % coord(3,ipoin) = bobox(2,3,ielem)
          ipoin                 = ipoin + 1
          mesh % coord(1,ipoin) = bobox(1,1,ielem)
          mesh % coord(2,ipoin) = bobox(2,2,ielem)
          mesh % coord(3,ipoin) = bobox(2,3,ielem)
       end do
    end if

  end subroutine mesh_bb

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Renumber
  !> @details Renumber a mesh using element and node permutaiton arrays
  !>          This permutations can be given by argument or using specific
  !>          tags
  !> 
  !-----------------------------------------------------------------------
  
  subroutine renumber(mesh,ELEMENT_PERM,NODE_PERM,ELEMENT_USING_TAG,NODE_USING_TAG)
    
    class(mesh_type_basic),                    intent(inout) :: mesh
    integer(ip),            optional, pointer, intent(in)    :: ELEMENT_PERM(:)
    integer(ip),            optional, pointer, intent(in)    :: NODE_PERM(:)
    integer(ip),            optional,          intent(in)    :: ELEMENT_USING_TAG
    integer(ip),            optional,          intent(in)    :: NODE_USING_TAG
    integer(ip),                      pointer                :: perme_loc(:)
    integer(ip),                      pointer                :: permn_loc(:)
    integer(8)                                               :: memor_loc(2)
    integer(ip)                                              :: ielem,ipoin
    integer(ip)                                              :: jelem,jpoin,inode
    real(rp),                         pointer                :: coord(:,:)
    integer(ip),                      pointer                :: lninv_loc(:)
    integer(ip),                      pointer                :: permn(:)
    integer(ip),                      pointer                :: leinv_loc(:)
    integer(ip),                      pointer                :: ltype(:)
    integer(ip),                      pointer                :: lnods(:,:)
    integer(ip),                      pointer                :: perme(:)
    
    nullify(perme_loc)
    nullify(permn_loc)
    
    if( present(ELEMENT_PERM) ) then
       perme_loc => ELEMENT_PERM
    else if( present(ELEMENT_USING_TAG) ) then
       call memory_alloca(memor_loc,'PERME_LOC','renumber',perme_loc,mesh % nelem)
       do ielem = 1,mesh % nelem
          perme_loc(ielem) = int(mesh % tags(ELEMENT_USING_TAG) % values(ielem),ip)
       end do
    end if
       
    if( present(NODE_PERM) ) then
       permn_loc => NODE_PERM
    else if( present(NODE_USING_TAG) ) then
       call memory_alloca(memor_loc,'PERMN_LOC','renumber',permn_loc,mesh % npoin)
       do ipoin = 1,mesh % npoin
          permn_loc(ipoin) = int(mesh % tags(NODE_USING_TAG) % values(ipoin),ip)
       end do
    end if
    
    if( associated(perme_loc) .and. mesh % nelem > 0 ) then
       allocate(lnods    (mesh % mnode,mesh % nelem))
       allocate(leinv_loc(mesh % nelem))
       allocate(ltype    (mesh % nelem))
       lnods(1:mesh % mnode,1:mesh % nelem) = mesh % lnods    (1:mesh % mnode,1:mesh % nelem)
       ltype    (1:mesh % nelem)            = mesh % ltype    (1:mesh % nelem)
       leinv_loc(1:mesh % nelem)            = mesh % leinv_loc(1:mesh % nelem)
       do ielem = 1,mesh % nelem
          jelem                              = perme_loc(ielem)
          mesh % lnods(1:mesh % mnode,jelem) = lnods(1:mesh % mnode,ielem)
          mesh % leinv_loc(jelem)            = leinv_loc(ielem)
       end do
       deallocate(lnods,ltype,leinv_loc)
        if( associated(mesh % perme) ) then
          allocate(perme    (mesh % nelem))
          do ielem = 1,mesh % nelem
             jelem = perme_loc(ielem)
             mesh % perme(jelem) = perme(ielem)
          end do
          deallocate(perme)
       end if
       if( associated(permn_loc) .and. mesh % npoin > 0 ) then
          do ielem = 1,mesh % nelem
             do inode = 1,element_type(mesh % ltype(ielem)) % number_nodes
                ipoin = mesh % lnods(inode,ielem)
                mesh % lnods(inode,ielem) = permn_loc(ipoin)
             end do
          end do
       end if
    end if
    
    if( associated(permn_loc) .and. mesh % npoin > 0 ) then
       allocate(coord    (mesh % ndime,mesh % npoin))
       allocate(lninv_loc(mesh % npoin))
       coord(1:mesh % ndime,1:mesh % npoin) = mesh % coord    (1:mesh % ndime,1:mesh % npoin)
       lninv_loc(1:mesh % npoin)            = mesh % lninv_loc(1:mesh % npoin)
       do ipoin = 1,mesh % npoin
          jpoin                              = permn_loc(ipoin)
          mesh % coord(1:mesh % ndime,jpoin) = coord(1:mesh % ndime,ipoin)
          mesh % lninv_loc(jpoin)            = lninv_loc(ipoin)
       end do
       deallocate(coord,lninv_loc)
       if( associated(mesh % permn) ) then
          allocate(permn    (mesh % npoin))
          do ipoin = 1,mesh % npoin
             jpoin = permn_loc(ipoin)
             mesh % permn(jpoin) = permn(ipoin)
          end do
          deallocate(permn)
       end if
    end if
    
    if( present(ELEMENT_USING_TAG) ) call memory_deallo(memor_loc,'PERME_LOC','renumber',perme_loc)
    if( present(NODE_USING_TAG)    ) call memory_deallo(memor_loc,'PERMN_LOC','renumber',permn_loc)
        
  end subroutine renumber
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Collapse nodes
  !> @details Collpase nodes of a mesh using coordinates
  !> 
  !-----------------------------------------------------------------------
  
  subroutine collapse(mesh,MEMORY_COUNTER)

    class(mesh_type_basic),           intent(inout) :: mesh
    integer(8),             optional, intent(inout) :: MEMORY_COUNTER(2)
    type(maths_bin)                                 :: bin
    integer(ip)                                     :: ipoin
    integer(ip),            pointer                 :: list_nodes(:)
    integer(ip)                                     :: ii,jpoin,lninv_max
    integer(ip)                                     :: imast,npoin_new,kpoin
    integer(ip)                                     :: pelty,inode,ielem,itag
    integer(ip)                                     :: npoin_old
    integer(ip),            pointer                 :: lmast(:)
    real(rp),               pointer                 :: coord_tmp(:,:)
    integer(ip),            pointer                 :: permn_tmp(:)
    integer(ip),            pointer                 :: lninv_tmp(:)
    real(rp),               pointer                 :: tag_tmp(:)
    integer(8)                                      :: memor_loc(2)
   
    memor_loc   = optional_argument((/0_8,0_8/),MEMORY_COUNTER)

    nullify(lmast)
    nullify(coord_tmp)
    nullify(permn_tmp)
    nullify(lninv_tmp)
    nullify(tag_tmp)

    if( mesh % npoin > 0 ) then
 
       call bin % init  ()
       call bin % input (LIMIT=20_ip,MAX_BINS=1000000_ip)
       call bin % fill  (mesh % coord)
       !
       ! Keep nodes with maximum LNINV_LOC
       !
       call memory_alloca(memor_loc,'LMAST','mesh_type_basic_valid_mesh',lmast,mesh % npoin)
       do ipoin = 1,mesh % npoin
          if( lmast(ipoin) == 0 ) then
             list_nodes => mesh % identical_coord(mesh % coord(:,ipoin),bin)
             if( size(list_nodes) > 1 ) then
                lninv_max = maxval(mesh % lninv_loc(list_nodes(:)))
                !
                ! Look for master
                !
                imast = 0
                ii    = 0
                do while( imast == 0 .and. ii < size(list_nodes,KIND=ip) )
                   ii    = ii + 1
                   jpoin = list_nodes(ii)
                   if( mesh % lninv_loc(jpoin) == lninv_max ) then
                      imast = jpoin
                   end if
                end do
                !
                ! Set slave-master
                !
                do ii = 1,size(list_nodes,kind=ip)
                   jpoin = list_nodes(ii)
                   if( jpoin == imast ) then
                      lmast(jpoin) =  imast
                   else
                      lmast(jpoin) = -imast
                   end if
                end do
             end if
             deallocate(list_nodes)
          end if
       end do
       npoin_new = count(lmast>=0,KIND=ip)
       npoin_old = mesh % npoin

       if( npoin_new < npoin_old ) then
          !
          ! Check communication
          !
          if( associated(mesh % comm % bound_perm) ) then
             do ii = 1,size(mesh % comm % bound_perm,KIND=ip)
                ipoin = mesh % comm % bound_perm(ii)
                if( lmast(ipoin) < 0 ) call runend('WE HAVE ELIMINATED AN INTERFACE NODE!')
             end do
          end if
          !
          ! Renumber
          !
          kpoin = 0
          do ipoin = 1,mesh % npoin
             if( lmast(ipoin) >= 0 ) then
                kpoin        =  kpoin + 1
                lmast(ipoin) =  kpoin
             end if
          end do
          do ipoin = 1,mesh % npoin
             if( lmast(ipoin) < 0 ) then
                kpoin        = -lmast(ipoin)
                lmast(ipoin) =  lmast(kpoin)
             end if
          end do
          !
          ! Eliminate slave nodes
          !
          mesh % npoin = npoin_new
          call memory_copy  (memor_loc,trim(mesh % name)//' % COORD'    ,'mesh_type_basic_valid_mesh',mesh % coord    ,coord_tmp,COPY_NAME='COORD_TMP')
          call memory_copy  (memor_loc,trim(mesh % name)//' % PERMN'    ,'mesh_type_basic_valid_mesh',mesh % permn    ,permn_tmp,COPY_NAME='PERMN_TMP')
          call memory_copy  (memor_loc,trim(mesh % name)//' % LNINV_LOC','mesh_type_basic_valid_mesh',mesh % lninv_loc,lninv_tmp,COPY_NAME='LNINV_TMP')          
          call memory_alloca(memor_loc,trim(mesh % name)//' % COORD'    ,'mesh_type_basic_valid_mesh',mesh % coord    ,mesh % ndime,mesh % npoin)
          call memory_alloca(memor_loc,trim(mesh % name)//' % PERMN'    ,'mesh_type_basic_valid_mesh',mesh % permn    ,mesh % npoin)
          call memory_alloca(memor_loc,trim(mesh % name)//' % LNINV_LOC','mesh_type_basic_valid_mesh',mesh % lninv_loc,mesh % npoin)
          do ipoin = 1,npoin_old
             jpoin = lmast(ipoin)
             if( jpoin >= 0 ) then
                mesh % coord(1:mesh % ndime,jpoin) = coord_tmp(1:mesh % ndime,ipoin)
                mesh % permn(jpoin)                = permn_tmp(ipoin)
                mesh % lninv_loc(jpoin)            = lninv_tmp(ipoin)
             end if
          end do
          do ielem = 1,mesh % nelem
             pelty = mesh % ltype(ielem)
             do inode = 1,element_type(pelty) % number_nodes
                ipoin = mesh % lnods(inode,ielem)
                jpoin = lmast(ipoin)
                mesh % lnods(inode,ielem) = jpoin
             end do
          end do
          do itag = 1,mesh % ntags
             if( mesh % tags(itag) % type == NODE_TAG ) then
                call memory_copy  (memor_loc,trim(mesh % name)//' % TAGS % VALUES','mesh_type_basic_valid_mesh',mesh % tags(itag) % values,tag_tmp,COPY_NAME='TAG_TMP')
                call memory_alloca(memor_loc,trim(mesh % name)//' % TAGS % VALUES','mesh_type_basic_valid_mesh',mesh % tags(itag) % values,mesh % npoin)
                do ipoin = 1,npoin_old
                   jpoin = lmast(ipoin)
                   mesh % tags(itag) % values(jpoin) = tag_tmp(ipoin)
                end do
                call memory_deallo(memor_loc,'TAG_TMP','mesh_type_basic_valid_mesh',tag_tmp)          
             end if
          end do

          call memory_deallo(memor_loc,'COORD_TMP','mesh_type_basic_valid_mesh',coord_tmp)
          call memory_deallo(memor_loc,'PERMN_TMP','mesh_type_basic_valid_mesh',permn_tmp)
          call memory_deallo(memor_loc,'LNINV_TMP','mesh_type_basic_valid_mesh',lninv_tmp)          

       end if
       call memory_deallo(memor_loc,'LMAST','mesh_type_basic_valid_mesh',lmast)
       call bin % deallo()

    end if

    if( present(MEMORY_COUNTER) )  MEMORY_COUNTER = memor_loc
    
  end subroutine collapse

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Return a list of nodes with same coordinates as another
  !> @details Return a list of nodes with same coordinates as another
  !> 
  !-----------------------------------------------------------------------

  pure function identical_coord(mesh,xx,bin,TOLERANCE) result(res)

    class(mesh_type_basic),                   intent(in) :: mesh
    real(rp),                                 intent(in) :: xx(:)
    class(maths_bin),                         intent(in) :: bin
    real(rp),              optional,          intent(in) :: TOLERANCE
    integer(ip),           pointer                       :: res(:)
    integer(ip)                                          :: ll(3)
    integer(ip)                                          :: isize,ii,ib,num_nodes,idime
    integer(ip)                                          :: jpoin
    integer(ip),           allocatable                   :: list_tmp(:)
    real(rp)                                             :: epsil
    logical(lg)                                          :: ifoun
    
    epsil = optional_argument(epsilon(1.0_rp),TOLERANCE)
    nullify(res)
    
    call bin % find(xx,ll)

    isize = bin % num_points(ll)
    
    if( isize > 0 ) then
       num_nodes = 0
       allocate(list_tmp(isize))
       if( bin % fill_type == LINKED_LIST_BIN ) then
          ib = bin % box_num(ll)
          do ii = bin % ia(ib),bin % ia(ib+1)-1
             jpoin = bin % ja(ii)
             ifoun = .true.
             idime_loop1: do idime = 1,mesh % ndime
                if( abs(mesh % coord(idime,jpoin)-xx(idime)) > epsil ) then
                   ifoun = .false.
                   exit idime_loop1
                end if
             end do idime_loop1
             if( ifoun ) then
                num_nodes = num_nodes + 1
                list_tmp(num_nodes) = jpoin
             end if
          end do
       else
          do ii = 1,isize
             jpoin = bin % list(ll(1),ll(2),ll(3)) % l(ii)
             ifoun = .true.
             idime_loop2: do idime = 1,mesh % ndime
                if( abs(mesh % coord(idime,jpoin)-xx(idime)) > epsil ) then
                   ifoun = .false.
                   exit idime_loop2
                end if
             end do idime_loop2
             if( ifoun ) then
                num_nodes = num_nodes + 1
                list_tmp(num_nodes) = jpoin
             end if
          end do
       end if
       
       if( num_nodes > 0 ) then
          allocate( res(num_nodes) )
          res(1:num_nodes) = list_tmp(1:num_nodes)
       end if
       if( allocated(list_tmp) ) deallocate(list_tmp)
    end if

  end function identical_coord

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2021-10-22
  !> @brief   Check if mesh is correct
  !> @details Check if mesh is correct
  !>          Error=1 ... Wrong LTYPE
  !>          Error=2 ... Wrong LNODS
  !>          Error=3 ... Wrong COORD DIMENSIONS
  !> 
  !-----------------------------------------------------------------------

  pure function check(mesh) result(ierr)

    class(mesh_type_basic), intent(in) :: mesh
    integer(ip)                        :: ierr
    logical(lg)                        :: res
    integer(ip)                        :: ielem
    integer(ip)                        :: num_ini,num_end

    num_ini = element_num_ini(mesh % ndime)
    num_end = element_num_end(mesh % ndime)    
    ierr    = 0
    !
    ! Check LTYPE
    !
    do ielem = 1,mesh % nelem
       if( mesh % ltype(ielem) < num_ini .or. mesh % ltype(ielem) > num_end ) &
            ierr = 1
    end do
    !
    ! Check LNODS
    !
    if( ierr == 0 ) then
       res = mesh % check_lnods()
       if( .not. res ) ierr = 2
    end if
    !
    ! COORD
    !
    if( mesh % npoin > 0 .and. int(size(mesh % coord,1),ip) /= mesh % ndime ) ierr = 3
    
  end function check
     
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2021-10-22
  !> @brief   Check if lnods is correct
  !> @details Check if lnods is correct. Returns true of this is the case
  !> 
  !-----------------------------------------------------------------------

  pure function check_lnods(mesh) result(res)

    class(mesh_type_basic), intent(in) :: mesh
    logical(lg)                        :: res
    integer(ip)                        :: ielem,inode
    
    res = .true.
    do ielem = 1,mesh % nelem
       do inode = 1,element_type(abs(mesh % ltype(ielem))) % number_nodes
          if( mesh % lnods(inode,ielem) <= 0 .or. mesh % lnods(inode,ielem) > mesh % npoin ) then
             res = .false.
             return
          end if
       end do
    end do

  end function check_lnods

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Centroid
  !> @details Compute the element centroids of a mesh
  !> 
  !-----------------------------------------------------------------------
  
  subroutine centroid(mesh,xx,MEMORY_COUNTER,ONLY_DEALLOCATE)

    class(mesh_type_basic),                    intent(inout) :: mesh
    real(rp),                         pointer, intent(inout) :: xx(:,:)
    integer(8),             optional,          intent(inout) :: MEMORY_COUNTER(2)
    logical(lg),            optional,          intent(in)    :: ONLY_DEALLOCATE
    integer(8)                                               :: memor_loc(2)
    integer(ip)                                              :: pdime,ielem,pelty,pnode
    integer(ip)                                              :: ipoin,inode
    real(rp)                                                 :: yy(3)

    call memory_counter_ini(memor_loc,mesh % memor,MEMORY_COUNTER)

    if( optional_argument(.false.,ONLY_DEALLOCATE) ) then
       call memory_deallo(memor_loc,'XX','centroid',xx)
    else
       pdime = mesh % ndime
       call memory_alloca(memor_loc,'XX','centroid',xx,pdime,mesh % nelem)
       do ielem = 1,mesh % nelem
          pelty = mesh % ltype(ielem)
          pnode = element_type(pelty) % number_nodes       
          yy    = 0.0_rp
          do inode = 1,pnode
             ipoin       = mesh % lnods(inode,ielem)
             yy(1:pdime) = yy(1:pdime) + mesh % coord(1:pdime,ipoin)
          end do
          xx(1:pdime,ielem) = yy(1:pdime) / real(pnode,rp)
       end do
    end if
 
    call memory_counter_end(memor_loc,mesh % memor,MEMORY_COUNTER)
    
  end subroutine centroid

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Centroid
  !> @details Compute the element centroids of a mesh
  !> 
  !-----------------------------------------------------------------------
  
  subroutine centroid_bmsh(mesh,xx,MEMORY_COUNTER,ONLY_DEALLOCATE)

    class(bmsh_type_basic),                    intent(inout) :: mesh
    real(rp),                         pointer, intent(inout) :: xx(:,:)
    integer(8),             optional,          intent(inout) :: MEMORY_COUNTER(2)
    logical(lg),            optional,          intent(in)    :: ONLY_DEALLOCATE
    integer(8)                                               :: memor_loc(2)
    integer(ip)                                              :: pdime,ielem,pelty,pnode
    integer(ip)                                              :: ipoin,inode
    real(rp)                                                 :: yy(3)

    call memory_counter_ini(memor_loc,mesh % memor,MEMORY_COUNTER)

    if( optional_argument(.false.,ONLY_DEALLOCATE) ) then
       call memory_deallo(memor_loc,'XX','centroid',xx)
    else
       pdime = mesh % ndime
       call memory_alloca(memor_loc,'XX','centroid',xx,pdime,mesh % nelem)
       do ielem = 1,mesh % nelem
          pelty = mesh % ltype(ielem)
          pnode = element_type(pelty) % number_nodes       
          yy    = 0.0_rp
          do inode = 1,pnode
             ipoin       = mesh % lnods(inode,ielem)
             yy(1:pdime) = yy(1:pdime) + mesh % coord(1:pdime,ipoin)
          end do
          xx(1:pdime,ielem) = yy(1:pdime) / real(pnode,rp)
       end do
    end if

    call memory_counter_end(memor_loc,mesh % memor,MEMORY_COUNTER)
    
  end subroutine centroid_bmsh

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-24
  !> @brief   Read mesh
  !> @details Read mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine read_from_file(mesh,lunit,filename,CLOSE_UNIT)

    class(mesh_type_basic),           intent(inout) :: mesh
    integer(ip),            optional, intent(in)    :: lunit               !< Output unit
    character(LEN=*),       optional, intent(in)    :: filename            !< Filename
    logical(lg),            optional, intent(in)    :: CLOSE_UNIT
    integer(ip)                                     :: kunit,jelty
    integer(ip)                                     :: ipara,itypb
    integer(ip)                                     :: ielty,nelty
    integer(ip)                                     :: ielem,nboun
    integer(ip)                                     :: ipoin,inode
    integer(ip)                                     :: ioerr
    character(5)                                    :: welem
    integer(4)                                      :: iunit4
    logical(lg)                                     :: opened
    logical(lg)                                     :: lexis(element_max)
    logical(lg)                                     :: if_strategy
    !
    ! Set unit
    !
    if( present(lunit) .and. present(filename) ) then
       kunit  = lunit
       iunit4 = int(kunit,4)
       open(iunit4,file=filename,status='unknown')
    else if( present(lunit) ) then
       kunit  = lunit
       iunit4 = int(kunit,4)
    else if( present(filename) ) then
       do iunit4 = 90_4,1000_4
          inquire(unit=iunit4,opened=opened,iostat=ioerr)
          if( ioerr /= 0 )   cycle
          if( .not. opened ) exit
       end do
       kunit = int(iunit4,ip)
       open(iunit4,file=filename,status='unknown')
    else
       continue ! Use current unit
    end if
    !
    ! Reach beginning and set name
    !
    if( present(lunit) .or. present(filename) ) call ecoute_set_read_unit(kunit)
    call ecoute_reach_section('MESH ')
    if(      words(2) == 'NAME ' ) then
       mesh % name = getcha('NAME ','MESH ','#Mesh name')
    else if( words(2) /= '     ' ) then
       mesh % name = words(2)
    end if
    lexis       = .false.
    if_strategy = .false.
    nboun       = 0
    itypb       = 0
    
    do while( words(1) /= 'ENDME' )

       select case ( words(1) )

       case ( 'DIMEN' )

          !-----------------------------------------------------------------------
          !
          ! Read dimensions
          !
          !-----------------------------------------------------------------------         

          do while( words(1) /= 'ENDDI' )
             
             select case ( words(1) )

             case ( 'NODAL' ) ; mesh % npoin = getint('NODAL',0_ip,'*NUMBER OF NODAL POINTS')         
             case ( 'ELEME' ) ; mesh % nelem = getint('ELEME',0_ip,'*NUMBER OF ELEMENTS')         
             case ( 'MAXIM' ) ; mesh % mnode = getint('MAXIM',0_ip,'*MAXIMUM NUMBER OF NODES PER ELEMENT')         
             case ( 'SPACE' ) ; mesh % ndime = getint('SPACE',2_ip,'*SPACE DIMENSION')
             case ( 'BOUND' ) ; nboun        = getint('BOUND',0_ip,'*NUMBER OF BOUNDARIES')

             case ( 'TAGS ' )

                mesh % ntags = getint('TAGS ',1_ip,'*NUMBER OF TAGS')
                call mesh % alloca_tags()
                call ecoute('read_from_file')
                do while( words(1) /= 'ENDTA' )
                   if( words(1) == 'TAG  ' ) then
                      ipara = getint('TAG  ',1_ip,'*TAG NUMBER')
                      if(      exists('ELEME') ) then
                         mesh % tags(ipara) % type = ELEMENT_TAG
                      else if( exists('NODES') ) then
                         mesh % tags(ipara) % type = NODE_TAG
                      else if( exists('BOUND') ) then
                         mesh % tags(ipara) % type = BOUNDARY_TAG
                      else
                         call runend('READ_FROM_FILE: UNKNOWN TAG TYPE')
                      end if
                   end if
                   call ecoute('read_from_file')
                end do

             case ( 'TYPES' )

                nelty = size(element_type)
                if( nnpar == 0 ) then
                   ipara = 2
                   do while( trim(words(ipara)) /= '' )
                      call elmgeo_element_name_to_type(words(ipara),ielty)
                      if( ielty >= 2 .and. ielty <= nelty ) then
                         mesh % mnode = max(mesh % mnode,element_type(ielty) % number_nodes)
                         lexis(ielty) = .true.
                      end if
                      ipara = ipara + 1
                   end do
                else
                   do ipara = 1,nnpar
                      ielty = int(param(ipara),ip)
                      lexis(ielty) = .true.
                      if( ielty >= 2 .and. ielty <= nelty ) &
                           mesh % mnode = max(mesh % mnode,element_type(ielty) % number_nodes)
                   end do
                end if

             end select
             call ecoute('read_from_file')

          end do

       case ( 'STRAT' )

          !-----------------------------------------------------------------------
          !
          ! Read strategy
          !
          !-----------------------------------------------------------------------         

          if_strategy = .true.

          do while( words(1) /= 'ENDST' )

             select case ( words(1) )

             case( 'INTEG' )

                select case ( words(2) ) 
                case( 'OPEN ' , 'GAUSS' ) ; mesh % quad(:) % type = GAUSS_LEGENDRE_RULE  
                case( 'CLOSE'           ) ; mesh % quad(:) % type = CLOSED_RULE 
                case( 'TRAPE'           ) ; mesh % quad(:) % type = TRAPEZOIDAL_RULE
                case( 'CHEBY'           ) ; mesh % quad(:) % type = CHEBYSHEV_RULE
                case default              ; call runend('UNKNOWN QUADRATURE RULE')
                end select

             case( 'INTER' )

                select case ( words(2) ) 
                case( 'LAGRA' ) ; mesh % iso(:) % inter = LAGRANGE_INTERPOLATION 
                case( 'CHEBY' ) ; mesh % iso(:) % inter = CHEBYSHEV_INTERPOLATION 
                case default              ; call runend('UNKNOWN INTERPOLATION')
                end select

             case( 'DOMAI' )

                jelty = 0
                do ielty = 1,nelty
                   if( lexis(ielty) ) then
                      jelty = jelty + 1
                      mesh % quad(ielty) % ngaus = int(param(jelty),ip)
                   end if
                end do

             case( 'GAUSS' )  

                call ecoute('read_from_file')
                do while( words(1) /= 'ENDGA' )
                   call elmgeo_element_name_to_type(words(1),ielty)
                   if( ielty < 1 .or. ielty > nelty ) call runend('REASTR: WRONG ELEMENT TYPE')
                   mesh % quad(ielty) % ngaus = int(param(1),ip)
                   call ecoute('REAST')
                end do

             end select

             call ecoute('read_from_file')

          end do

       case ( 'GEOME' )

          !-----------------------------------------------------------------------
          !
          ! Read geometry
          !
          !-----------------------------------------------------------------------

          if( nboun > 0 ) then
             allocate(mesh % boundary)
             call mesh % boundary % init()
             mesh % boundary % nelem = nboun
             mesh % boundary % npoin = mesh % npoin
             mesh % boundary % ndime = mesh % ndime
             do ielty = 1,size(lexis)
                if( lexis(ielty) ) &
                     mesh % boundary % mnode = max(mesh % boundary % mnode,element_type(ielty) % max_face_nodes)
             end do
             call mesh % boundary % alloca()
          end if

          call mesh % alloca(ALL=.true.,TAGS=.false.)
          call mesh % alloca_tag()
         
          do while( words(1) /= 'ENDGE' )

             select case ( words(1) )

             case ( 'NODES' )

                call ecoute('read_from_file')
                do while( words(1) /= 'ENDNO' )
                   ielem               = int(param(1),ip)
                   inode               = int(param(2),ip)
                   ielty               = elmgeo_element_type(mesh % ndime,inode)
                   mesh % ltype(ielem) = ielty
                   call ecoute('read_from_file')
                end do

             case ( 'TYPES' )

                if( words(2) == 'ALL  ' ) then
                   !
                   ! All elements are of the same type
                   !
                   ielty = int(param(2),ip)
                   if( ielty == 0 ) then
                      welem = getcha('ALL  ','NULL ','#Element type')
                      call elmgeo_element_name_to_type(welem,ielty)
                   end if
                   do ielem = 1,mesh % nelem
                      mesh % ltype(ielem) = ielty
                   end do
                   call ecoute_reach_section('ENDTY')

                else
                   !
                   ! Alya type as defined in def_elmtyp
                   !
                   call ecoute('read_from_file')
                   do while( words(1) /= 'ENDTY' )
                      ielem               = int(param(1),ip)
                      ielty               = int(param(2),ip)
                      mesh % ltype(ielem) = ielty
                      call ecoute('read_from_file')
                   end do
                end if

             case ( 'COORD' )

                call ecoute('read_from_file')
                do while( words(1) /= 'ENDCO' )
                   ipoin                              = int(param(1),ip)
                   mesh % coord(1:mesh % ndime,ipoin) = param(2:1+mesh % ndime)
                   call ecoute('read_from_file')
                end do

             case ( 'ELEME' )

                call ecoute('read_from_file')
                do while( words(1) /= 'ENDEL' )
                   ielem                       = int(param(1),ip)
                   inode                       = element_type(mesh % ltype(ielem)) % number_nodes
                   mesh % lnods(1:inode,ielem) = int(param(2:1+inode),ip)             
                   call ecoute('read_from_file')
                end do

             case ( 'BOUND' )

                call ecoute('read_from_file')
                if( itypb == 0 ) then
                   do while( words(1) /= 'ENDBO' )
                      ielem                                  = int(param(1),ip)
                      inode                                  = nnpar-1
                      mesh % boundary % ltype(ielem)         = elmgeo_element_type(mesh % ndime-1,inode)
                      inode                                  = element_type(mesh % boundary % ltype(ielem)) % number_nodes
                      mesh % boundary % lnods(1:inode,ielem) = int(param(2:1+inode),ip)             
                      call ecoute('read_from_file')
                   end do
                else
                   do while( words(1) /= 'ENDBO' )
                      ielem                                  = int(param(1),ip)
                      inode                                  = element_type(mesh % boundary % ltype(ielem)) % number_nodes
                      mesh % boundary % lnods(1:inode,ielem) = int(param(2:1+inode),ip)             
                      call ecoute('read_from_file')
                   end do                   
                end if
                
             end select

             call ecoute('read_from_file')

          end do

       case ( 'TAGS ' )

          !-----------------------------------------------------------------------
          !
          ! Read tags
          !
          !-----------------------------------------------------------------------

          call ecoute('read_from_file')
          do while( words(1) /= 'ENDTA' )
             if( words(1) == 'TAG  ' ) then
                ipara = getint('TAG  ',1_ip,'*TAG NUMBER')
                do while( words(1) /= 'ENDTA' )
                   ielem = int(param(1),ip)
                   mesh % tags(ipara) % values(ielem) = param(2)
                   call ecoute('read_from_file')
                end do
             end if
             call ecoute('read_from_file')             
          end do

       case ( 'EOF  ' )

          !-----------------------------------------------------------------------
          !
          ! End of file
          !
          !-----------------------------------------------------------------------

          call runend('END OF FILE')

       end select

       call ecoute('read_from_file')

    end do
    !
    ! Integration rules and iso-parametric arrays
    !
    if( if_strategy ) then
       do ielty = 1,size(element_type)
          if( lexis(ielty) ) then
             if( mesh % quad(ielty) % ngaus == 0 ) &
                  mesh % quad(ielty) % ngaus = element_type(ielty) % number_nodes
             mesh % quad(ielty) % topo  = element_type(ielty) % topology
             mesh % quad(ielty) % ndime = element_type(ielty) % dimensions
             call mesh % quad(ielty) % set()
             call mesh % iso(ielty)  % set(element_type(ielty) % number_nodes,mesh % quad(ielty))
          end if
       end do
    end if
    !
    ! Close file
    !
    if( optional_argument(.false.,CLOSE_UNIT) ) close(iunit4)

  end subroutine read_from_file

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-08-01
  !> @brief   Transform a mesh
  !> @details Transform a generic mesh into a linear mesh
  !> 
  !-----------------------------------------------------------------------
  
  subroutine linearize(mesh,mesh_in)

    class(mesh_type_basic), intent(inout) :: mesh
    class(mesh_type_basic), intent(in)    :: mesh_in
    integer(ip)                           :: ielem,pelty,kelem
    integer(ip)                           :: ipoin,ii
    integer(ip),            allocatable   :: lnuty(:)
    integer(ip), pointer                  :: SOURCE_TO_TARGET(:,:)
    integer(ip)                           :: kelty,pelem,knode
    
    if( mesh_in % nelem > 0 ) then

       allocate(lnuty(size(element_type)))
       lnuty = 0
       do ielem = 1,mesh_in % nelem
          pelty        = mesh_in % ltype(ielem)
          lnuty(pelty) = lnuty(pelty) + 1
       end do
    
       mesh % mnode = 0
       mesh % ndime = mesh_in % ndime
       mesh % npoin = mesh_in % npoin
       mesh % nelem = lnuty(BAR02) + lnuty(BAR03) * 2 + lnuty(BAR04) * 3                                   + &
            &         lnuty(TRI03) + lnuty(TRI06) * 4 + lnuty(QUA04) + lnuty(QUA09) * 4 + lnuty(QUA16) * 9 + &
            &         lnuty(TET04) + lnuty(TET10) * 8                                                      + &
            &         lnuty(PYR05) + lnuty(PEN06)                                                          + &
            &         lnuty(HEX08) + lnuty(HEX27) * 8
       do pelty = 1,size(element_type)
          if( lnuty(pelty) > 0 ) mesh % mnode = max(mesh % mnode,element_type(pelty) % number_nodes)
       end do
       
       call mesh % alloca()       
       do ipoin = 1,mesh_in % npoin
          mesh % coord    (:,ipoin) = mesh_in % coord    (:,ipoin) 
          mesh % lninv_loc(  ipoin) = mesh_in % lninv_loc(  ipoin) 
       end do
       do kelem = 1,mesh % nelem
          mesh % leinv_loc(kelem) = kelem
       end do
       
       kelem = 0
       
       do ielem = 1,mesh_in % nelem
          pelty = mesh_in % ltype(ielem)
          select case ( pelty )
          case ( BAR02 ) ; SOURCE_TO_TARGET => BAR02_TO_BAR02 ; kelty = BAR02
          case ( BAR03 ) ; SOURCE_TO_TARGET => BAR03_TO_BAR02 ; kelty = BAR02
          case ( BAR04 ) ; SOURCE_TO_TARGET => BAR04_TO_BAR02 ; kelty = BAR02
          case ( TRI03 ) ; SOURCE_TO_TARGET => TRI03_TO_TRI03 ; kelty = TRI03
          case ( TRI06 ) ; SOURCE_TO_TARGET => TRI06_TO_TRI03 ; kelty = TRI03
          case ( QUA04 ) ; SOURCE_TO_TARGET => QUA04_TO_QUA04 ; kelty = QUA04
          case ( QUA09 ) ; SOURCE_TO_TARGET => QUA09_TO_QUA04 ; kelty = QUA04
          case ( QUA16 ) ; SOURCE_TO_TARGET => QUA16_TO_QUA04 ; kelty = QUA04
          case ( TET04 ) ; SOURCE_TO_TARGET => TET04_TO_TET04 ; kelty = TET04
          case ( TET10 ) ; SOURCE_TO_TARGET => TET10_TO_TET04 ; kelty = TET04
          case ( PYR05 ) ; SOURCE_TO_TARGET => PYR05_TO_PYR05 ; kelty = PYR05
          case ( PEN06 ) ; SOURCE_TO_TARGET => PEN06_TO_PEN06 ; kelty = PEN06
          case ( HEX08 ) ; SOURCE_TO_TARGET => HEX08_TO_HEX08 ; kelty = HEX08
          case ( HEX27 ) ; SOURCE_TO_TARGET => HEX27_TO_HEX08 ; kelty = HEX08
          end select
          knode = size(SOURCE_TO_TARGET,1)
          pelem = size(SOURCE_TO_TARGET,2)
          do ii = 1,pelem
             kelem                       = kelem + 1
             mesh % ltype(kelem)         = kelty
             mesh % lnods(1:knode,kelem) = mesh_in % lnods(SOURCE_TO_TARGET(:,ii),ielem)
          end do
       end do

       deallocate(lnuty)
       
    end if

  end subroutine linearize
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-08-01
  !> @brief   Transform a mesh
  !> @details Transform a generic mesh into a linear mesh
  !> 
  !-----------------------------------------------------------------------
  
  subroutine tet04_mesh(mesh,mesh_in)

    class(mesh_type_basic), intent(inout) :: mesh
    class(mesh_type_basic), intent(in)    :: mesh_in
    integer(ip)                           :: ielem,pelty,kelem
    integer(ip)                           :: ipoin,ii,jj
    integer(ip),            allocatable   :: lnuty(:)
    integer(ip), pointer                  :: SOURCE_TO_TARGET(:,:)
    integer(ip)                           :: pelem
    integer(ip)                           :: lnode(8)

    if( mesh_in % nelem > 0 ) then

       allocate(lnuty(size(element_type)))
       lnuty = 0
       do ielem = 1,mesh_in % nelem
          pelty        = mesh_in % ltype(ielem)
          lnuty(pelty) = lnuty(pelty) + 1
       end do
    
       mesh % mnode = 0
       mesh % ndime = mesh_in % ndime
       mesh % npoin = mesh_in % npoin
       mesh % nelem =   lnuty(TET04)                                                   &
            &         + lnuty(TET10) * size(TET10_TO_TET04,2)                          &
            &         + lnuty(PEN06) * size(PEN06_TO_TET04,2)                          &
            &         + lnuty(PYR05) * size(PYR05_TO_TET04,2)                          &
            &         + lnuty(HEX08) * size(HEX08_TO_TET04,2)                          &
            &         + lnuty(HEX27) * size(HEX27_TO_HEX08,2) * size(HEX08_TO_TET04,2)
       
       do pelty = 1,size(element_type)
          if( lnuty(pelty) > 0 ) mesh % mnode = max(mesh % mnode,element_type(pelty) % number_nodes)
       end do
       
       call mesh % alloca()       
       do ipoin = 1,mesh_in % npoin
          mesh % coord    (:,ipoin) = mesh_in % coord    (:,ipoin) 
          mesh % lninv_loc(  ipoin) = mesh_in % lninv_loc(  ipoin) 
       end do
       do kelem = 1,mesh % nelem
          mesh % leinv_loc(kelem) = kelem
       end do
       
       kelem = 0
       
       do ielem = 1,mesh_in % nelem
          pelty = mesh_in % ltype(ielem)
          select case ( pelty )
          case ( TET04 ) ; SOURCE_TO_TARGET => TET04_TO_TET04 
          case ( TET10 ) ; SOURCE_TO_TARGET => TET10_TO_TET04 
          case ( PYR05 ) ; SOURCE_TO_TARGET => PYR05_TO_TET04 
          case ( PEN06 ) ; SOURCE_TO_TARGET => PEN06_TO_TET04 
          case ( HEX08 ) ; SOURCE_TO_TARGET => HEX08_TO_TET04 
          end select

          if( mesh_in % ltype(ielem) == HEX27 ) then
             do ii = 1,size(HEX27_TO_HEX08,2)
                lnode(1:8) = mesh_in % lnods(HEX27_TO_HEX08(1:8,ii),ielem)
                do jj = 1,size(HEX08_TO_TET04,2)
                   kelem                   = kelem + 1
                   mesh % ltype(kelem)     = TET04                   
                   mesh % lnods(1:4,kelem) = lnode(HEX08_TO_TET04(1:4,jj))                   
                end do
             end do
          else             
             pelem = size(SOURCE_TO_TARGET,2)
             do ii = 1,pelem
                kelem                   = kelem + 1
                mesh % ltype(kelem)     = TET04
                mesh % lnods(1:4,kelem) = mesh_in % lnods(SOURCE_TO_TARGET(1:4,ii),ielem)
             end do
          end if
       end do
       
       deallocate(lnuty)
       
    end if

  end subroutine tet04_mesh
  
end module def_kintyp_mesh_basic
!> @}
