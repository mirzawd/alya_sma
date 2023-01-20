!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!>
!> @defgroup Elemental_Geometric_Toolbox
!> ToolBox for elemental and general geometrical operations
!> @{
!> @file    mod_elmgeo.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for elements
!> @details Different functions useful in finite element implementations
!>
!------------------------------------------------------------------------

module mod_elmgeo
 
#include "def_vector_size.inc"
  use def_kintyp_basic,      only : ip,rp,lg,i1p,i2p
#ifndef I_AM_NOT_ALYA
  use def_master,            only : kfl_paral
#endif
  use def_elmtyp,            only : POINT
  use def_elmtyp,            only : BAR02,BAR03,BAR04,TRI03,TRI06,TRI10,QUA04,QUA08,QUA09,QUA16
  use def_elmtyp,            only : TET04,TET10,TET20,PYR05,PYR14,PEN06,PEN15,PEN18,HEX08
  use def_elmtyp,            only : HEX27,HEX64,SHELL,BAR3D,POI3D
  use def_elmtyp,            only : element_num_ini 
  use def_elmtyp,            only : element_num_end 
  use def_elmtyp,            only : element_end 
  use mod_maths_solver,      only : maths_invert_matrix
  use mod_maths_geometry,    only : maths_vectorial_product
  use mod_maths_geometry,    only : maths_circumradius_tetrahedron     
  use mod_maths_geometry,    only : maths_circumradius_triangle        
  use mod_maths_geometry,    only : maths_circumradius_quadrilateral
  use mod_maths_geometry,    only : maths_distance_point_segment
  use mod_maths_arrays,      only : maths_maxloc_nonzero
  use mod_optional_argument, only : optional_argument
  use mod_bouder,            only : bouder
  use def_isoparametric,     only : LAGRANGE_INTERPOLATION
  use def_isoparametric,     only : set_isoparametric
  use def_elmgeo
  implicit none
  private

#ifdef I_AM_NOT_ALYA
#undef VECTOR_SIZE
  integer(ip), parameter :: VECTOR_SIZE = 1
  integer(ip), parameter :: kfl_paral   = 1
#endif
  
  real(rp),    parameter :: one3         = 1.0_rp/3.0_rp
  integer(ip), parameter :: perm1(3)     = (/2,1,1/)
  integer(ip), parameter :: perm2(3)     = (/3,3,2/)

  real(rp),    parameter :: epsilgeo_div = epsilon(1.0_rp) !epsilgeo used to avoid divisions by zero
  real(rp),    parameter :: epsilgeo_jac = 1.0e-20_rp      !epsilgeo used to avoid null jacobian in newrap
  real(rp),    parameter :: epsilgeo_alg = epsilon(1.0_rp) !epsilgeo used for algorithm branching
  real(rp),    parameter :: epsilgeo_def = epsilon(1.0_rp) !epsilgeo used as default tolerance when toler is an optional argument
  real(rp),    parameter :: epsilgeo_tol = epsilon(1.0_rp) !epsilgeo used as tolerance when it is not an argument
  !
  ! Parameters
  !
  integer(ip), parameter  :: FLAT_BUBBLE         = 1
  integer(ip), parameter  :: QUADRATIC_BUBBLE    = 2
  integer(ip), parameter  :: FREE_SURFACE_BUBBLE = 3
  !
  ! Scalar and vectorized versions
  !
  interface elmgeo_cartesian_derivatives_jacobian
     module procedure elmgeo_cartesian_derivatives_jacobian_scalar,&
          &           elmgeo_cartesian_derivatives_jacobian_vector
  end interface elmgeo_cartesian_derivatives_jacobian

  interface elmgeo_hessian
     module procedure elmgeo_hessian_scalar,&
          &           elmgeo_hessian_vector
  end interface elmgeo_hessian

  interface elmgeo_element_characteristic_length
     module procedure elmgeo_element_characteristic_length_scalar,&
          &           elmgeo_element_characteristic_length_vector
  end interface elmgeo_element_characteristic_length

  interface elmgeo_element_length
     module procedure elmgeo_element_length_scalar,&
          &           elmgeo_element_length_vector
  end interface elmgeo_element_length

  interface elmgeo_inside_element_bounding_box
     module procedure elmgeo_inside_element_bounding_box_1,&
          &           elmgeo_inside_element_bounding_box_2
  end interface elmgeo_inside_element_bounding_box

  interface elmgeo_gauss_to_element
     module procedure elmgeo_gauss_to_element_s,&
          &           elmgeo_gauss_to_element_1
  end interface elmgeo_gauss_to_element

  interface elmgeo_shapf_deriv_heslo
     module procedure elmgeo_shapf_deriv_heslo_1,&
          &           elmgeo_shapf_deriv_heslo_ngaus
  end interface elmgeo_shapf_deriv_heslo
  
  interface elmgeo_shapf_deriv_heslo_bubble
     module procedure elmgeo_shapf_deriv_heslo_bubble_1,&
          &           elmgeo_shapf_deriv_heslo_bubble_ngaus
  end interface elmgeo_shapf_deriv_heslo_bubble
  
  abstract interface
     subroutine elmgeo_jacobian(ndime,pnode,elcod,deriv,gpdet,xjaci,xjacm)
       import                             :: ip
       import                             :: rp
       integer(ip),           intent(in)  :: ndime               !< Dimension
       integer(ip),           intent(in)  :: pnode               !< Number of nodes
       real(rp),              intent(in)  :: elcod(ndime,pnode)  !< Element coordinates
       real(rp),              intent(in)  :: deriv(ndime,pnode)  !< Shape function derivatives
       real(rp),              intent(out) :: gpdet               !< Determinant
       real(rp),    optional, intent(out) :: xjaci(ndime,ndime)  !< Inverse Jacobian
       real(rp),    optional, intent(out) :: xjacm(ndime,ndime)  !< Jacobian
     end subroutine elmgeo_jacobian
  end interface
  
  public :: elmgeo_element_characteristic_length
  public :: elmgeo_element_length
  public :: elmgeo_segpla
  public :: elmgeo_segfac
  public :: elmgeo_bounor
  public :: elmgeo_natural_coordinates
  public :: elmgeo_shapf_deriv_heslo
  public :: elmgeo_inside_TET04
  public :: elmgeo_inside_TRI03_QUA04
  public :: elmgeo_newrap_norm
  public :: elmgeo_newrap
  public :: elmgeo_where_is
  public :: elmgeo_shapf_deriv_heslo_bubble
  public :: elmgeo_jacobian
  public :: elmgeo_jacobian_matrix
  public :: elmgeo_jacobian_boundary
  public :: elmgeo_natural_coordinates_on_boundaries
  public :: elmgeo_inside_element_bounding_box
  public :: elmgeo_inside_element_using_faces
  public :: elmgeo_element_type_initialization
  public :: elmgeo_nearest_point_on_element_faces
  public :: elmgeo_cartesian_derivatives
  public :: elmgeo_projection_on_a_face
  public :: elmgeo_intersection_segment_face
  public :: elmgeo_intersection_segment_QUA04
  public :: elmgeo_face_area
  public :: elmgeo_TET04_volume
  public :: elmgeo_PYR05_volume
  public :: elmgeo_element_volume
  public :: elmgeo_element_distance
  public :: elmgeo_bubble
  public :: elmgeo_nearest_point_on_QUA04
  public :: elmgeo_natural_coordinates_on_QUA04
  public :: elmgeo_nearest_intersection_point_on_element_faces
  public :: elmgeo_element_node_length
  public :: elem_typ
  public :: elmgeo_intersection_segment_TRI06
  public :: elmgeo_circumradius
  public :: elmgeo_gauss_to_element     ! Compute a value in an element from Gauss point values
  public :: elmgeo_number_nodes         ! Number of nodes
  public :: elmgeo_mm_elements          ! Number of subelements when using mesh multiplication
  public :: elmgeo_mm_interior_nodes    ! Number of interior nodes when using mesh multiplication
  public :: elmgeo_mm_face_nodes        ! Number of face nodes when using mesh multiplication
  public :: elmgeo_order_boundary_nodes ! Order face nodes
  public :: elmgeo_boundary_face        ! Face corresponding to a boundary
  public :: elmgeo_edge                 ! Return the edge number
  public :: elmgeo_face                 ! Return the face number
  public :: elmgeo_segment              ! Return a segment connectivity

  public :: list_faces_BAR02,type_faces_BAR02,list_edges_BAR02
  public :: list_faces_BAR03,type_faces_BAR03,list_edges_BAR03
  public :: list_faces_BAR04,type_faces_BAR04,list_edges_BAR04
  public :: list_faces_TRI03,type_faces_TRI03,list_edges_TRI03
  public :: list_faces_TRI06,type_faces_TRI06,list_edges_TRI06
  public :: list_faces_TRI10,type_faces_TRI10,list_edges_TRI10
  public :: list_faces_QUA04,type_faces_QUA04,list_edges_QUA04
  public :: list_faces_QUA08,type_faces_QUA08,list_edges_QUA08
  public :: list_faces_QUA09,type_faces_QUA09,list_edges_QUA09
  public :: list_faces_QUA16,type_faces_QUA16,list_edges_QUA16
  public :: list_faces_TET04,type_faces_TET04,list_edges_TET04
  public :: list_faces_TET10,type_faces_TET10,list_edges_TET10
  public :: list_faces_TET20,type_faces_TET20,list_edges_TET20
  public :: list_faces_PYR05,type_faces_PYR05,list_edges_PYR05
  public :: list_faces_PEN06,type_faces_PEN06,list_edges_PEN06
  public :: list_faces_PEN15,type_faces_PEN15,list_edges_PEN15
  public :: list_faces_PEN18,type_faces_PEN18,list_edges_PEN18
  public :: list_faces_HEX08,type_faces_HEX08,list_edges_HEX08
  public :: list_faces_HEX27,type_faces_HEX27,list_edges_HEX27
  public :: list_faces_HEX64,type_faces_HEX64,list_edges_HEX64
  public :: element_type
  public :: FLAT_BUBBLE
  public :: QUADRATIC_BUBBLE
  public :: FREE_SURFACE_BUBBLE
  public :: elmgeo_element_type                                  ! Guess element type
  public :: elmgeo_iso_parametric_element_coordinates            ! Coordinate of node of iso-parametric element
  public :: elmgeo_nearest_element_node                          ! Nearest element node
  public :: HEX08_TO_TET04                                       ! Subdivision of HEX08 into TET04
  public :: PEN06_TO_TET04                                       ! Subdivision of PEN06 into TET04
  public :: PYR05_TO_TET04                                       ! Subdivision of PYR05 into TET04
  public :: TET04_TO_TET04                                       ! Subdivision of TET04 into TET04

  
  public :: elmgeo_cartesian_derivatives_jacobian                ! Shape function derivatives, etc. (scalar and vectorized versions)
  public :: elmgeo_hessian                                       ! Shape function Hessian (scalar and vectorized versions)

  public :: elmgeo_output_jacobian                               ! Ouput the Jacobian of an element in GiD format
  public :: elmgeo_element_name_to_type

contains

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   Fill in element data base
  !> @details Element data base in mod_elmgeo. This redundant information
  !>          is required so that mod_elmgeo can be a stand-alone
  !>          module. Optional arguments can be used to test the
  !>          coherence between this database and Alya's one.
  !>
  !-----------------------------------------------------------------------

  subroutine elmgeo_element_type_initialization()

    integer(ip)            :: pnodf,porde
    integer(ip)            :: ielty,pface,pedge,pdime,pnode
    real(rp),    parameter :: x1ove3 = 1.0_rp/3.0_rp
    
    element_type(1:nelty) = element_type_init
    !
    ! top(col. 3) means topology. All qua and hex are 0, tri and tet 1 ...
    ! max(col. 6) means maximum number of nodes per face
    !
    !                                                            dim   ord  top  nnode    hess   npf nfac  nedg length                centroid
    element_type(POINT) = elem_typ(POINT,'POINT','Point'        ,0_ip,1_ip,-2_ip, 1_ip,.false., 0_ip,0_ip, 0_ip,1.0_rp,(/0.0_rp,0.0_rp,0.0_rp/),null(),null(),null(),null(),null(),null())
    element_type(BAR02) = elem_typ(BAR02,'BAR02','Linear'       ,1_ip,1_ip,-1_ip, 2_ip,.false., 1_ip,2_ip, 0_ip,2.0_rp,(/0.0_rp,0.0_rp,0.0_rp/),null(),null(),null(),null(),null(),null())
    element_type(BAR03) = elem_typ(BAR03,'BAR03','Linear'       ,1_ip,2_ip,-1_ip, 3_ip,.true. , 1_ip,2_ip, 0_ip,2.0_rp,(/0.0_rp,0.0_rp,0.0_rp/),null(),null(),null(),null(),null(),null())
    element_type(BAR04) = elem_typ(BAR04,'BAR04','Linear'       ,1_ip,3_ip,-1_ip, 4_ip,.true. , 1_ip,2_ip, 0_ip,2.0_rp,(/0.0_rp,0.0_rp,0.0_rp/),null(),null(),null(),null(),null(),null())
    element_type(TRI03) = elem_typ(TRI03,'TRI03','Triangle'     ,2_ip,1_ip, 1_ip, 3_ip,.false., 2_ip,3_ip, 3_ip,1.0_rp,(/x1ove3,x1ove3,0.0_rp/),null(),null(),null(),null(),null(),null())
    element_type(TRI06) = elem_typ(TRI06,'TRI06','Triangle'     ,2_ip,2_ip, 1_ip, 6_ip,.true. , 3_ip,3_ip, 3_ip,1.0_rp,(/x1ove3,x1ove3,0.0_rp/),null(),null(),null(),null(),null(),null())
    element_type(TRI10) = elem_typ(TRI10,'TRI10','Triangle'     ,2_ip,3_ip, 1_ip,10_ip,.false., 4_ip,3_ip, 3_ip,1.0_rp,(/x1ove3,x1ove3,0.0_rp/),null(),null(),null(),null(),null(),null())
    element_type(QUA04) = elem_typ(QUA04,'QUA04','Quadrilateral',2_ip,1_ip, 0_ip, 4_ip,.true. , 2_ip,4_ip, 4_ip,2.0_rp,(/0.0_rp,0.0_rp,0.0_rp/),null(),null(),null(),null(),null(),null())
    element_type(QUA08) = elem_typ(QUA08,'QUA08','Quadrilateral',2_ip,2_ip, 0_ip, 8_ip,.true. , 3_ip,4_ip, 4_ip,2.0_rp,(/0.0_rp,0.0_rp,0.0_rp/),null(),null(),null(),null(),null(),null())
    element_type(QUA09) = elem_typ(QUA09,'QUA09','Quadrilateral',2_ip,2_ip, 0_ip, 9_ip,.true. , 3_ip,4_ip, 4_ip,2.0_rp,(/0.0_rp,0.0_rp,0.0_rp/),null(),null(),null(),null(),null(),null())
    element_type(QUA16) = elem_typ(QUA16,'QUA16','Quadrilateral',2_ip,3_ip, 0_ip,16_ip,.true. , 4_ip,4_ip, 4_ip,2.0_rp,(/0.0_rp,0.0_rp,0.0_rp/),null(),null(),null(),null(),null(),null())
    element_type(TET04) = elem_typ(TET04,'TET04','Tetrahedra'   ,3_ip,1_ip, 1_ip, 4_ip,.false., 3_ip,4_ip, 6_ip,1.0_rp,(/x1ove3,x1ove3,x1ove3/),null(),null(),null(),null(),null(),null())
    element_type(TET10) = elem_typ(TET10,'TET10','Tetrahedra'   ,3_ip,2_ip, 1_ip,10_ip,.true. , 6_ip,4_ip, 6_ip,1.0_rp,(/x1ove3,x1ove3,x1ove3/),null(),null(),null(),null(),null(),null())
    element_type(TET20) = elem_typ(TET20,'TET20','Tetrahedra'   ,3_ip,3_ip, 1_ip,20_ip,.false.,10_ip,4_ip, 6_ip,1.0_rp,(/x1ove3,x1ove3,x1ove3/),null(),null(),null(),null(),null(),null())
    element_type(PYR05) = elem_typ(PYR05,'PYR05','Pyramid'      ,3_ip,1_ip, 3_ip, 5_ip,.true. , 4_ip,5_ip, 8_ip,1.0_rp,(/0.0_rp,0.0_rp,0.0_rp/),null(),null(),null(),null(),null(),null())
    element_type(PYR14) = elem_typ(PYR14,'PYR14','Pyramid'      ,3_ip,2_ip, 3_ip,14_ip,.true. , 0_ip,0_ip, 0_ip,1.0_rp,(/0.0_rp,0.0_rp,0.0_rp/),null(),null(),null(),null(),null(),null())
    element_type(PEN06) = elem_typ(PEN06,'PEN06','Prism'        ,3_ip,1_ip, 2_ip, 6_ip,.false., 4_ip,5_ip, 9_ip,1.0_rp,(/x1ove3,x1ove3,0.5_rp/),null(),null(),null(),null(),null(),null())
    element_type(PEN15) = elem_typ(PEN15,'PEN15','Prism'        ,3_ip,2_ip, 2_ip,15_ip,.true. , 8_ip,5_ip, 9_ip,1.0_rp,(/x1ove3,x1ove3,0.5_rp/),null(),null(),null(),null(),null(),null())
    element_type(PEN18) = elem_typ(PEN18,'PEN18','Prism'        ,3_ip,2_ip, 2_ip,18_ip,.true. , 9_ip,5_ip, 9_ip,1.0_rp,(/x1ove3,x1ove3,0.5_rp/),null(),null(),null(),null(),null(),null())
    element_type(HEX08) = elem_typ(HEX08,'HEX08','Hexahedra'    ,3_ip,1_ip, 0_ip, 8_ip,.true. , 4_ip,6_ip,12_ip,2.0_rp,(/0.0_rp,0.0_rp,0.0_rp/),null(),null(),null(),null(),null(),null())
    element_type(HEX27) = elem_typ(HEX27,'HEX27','Hexahedra'    ,3_ip,2_ip, 0_ip,27_ip,.true. , 9_ip,6_ip,12_ip,2.0_rp,(/0.0_rp,0.0_rp,0.0_rp/),null(),null(),null(),null(),null(),null())
    element_type(HEX64) = elem_typ(HEX64,'HEX64','Hexahedra'    ,3_ip,3_ip, 0_ip,64_ip,.true. ,16_ip,6_ip,12_ip,2.0_rp,(/0.0_rp,0.0_rp,0.0_rp/),null(),null(),null(),null(),null(),null())
    element_type(SHELL) = elem_typ(SHELL,'SHELL','Triangle'     ,2_ip,1_ip, 1_ip, 3_ip,.false., 2_ip,3_ip, 3_ip,1.0_rp,(/0.0_rp,0.0_rp,0.0_rp/),null(),null(),null(),null(),null(),null())
    element_type(BAR3D) = elem_typ(BAR3D,'BAR3D','Linear'       ,1_ip,1_ip,-1_ip, 2_ip,.false., 1_ip,2_ip, 1_ip,2.0_rp,(/0.0_rp,0.0_rp,0.0_rp/),null(),null(),null(),null(),null(),null())
    element_type(POI3D) = elem_typ(POI3D,'POI3D','Point'        ,0_ip,1_ip,-2_ip, 1_ip,.false., 1_ip,0_ip, 0_ip,1.0_rp,(/0.0_rp,0.0_rp,0.0_rp/),null(),null(),null(),null(),null(),null())

    do ielty = 1,nelty
       
       pface = element_type(ielty) % number_faces
       pnodf = element_type(ielty) % max_face_nodes
       pedge = element_type(ielty) % number_edges
       pdime = element_type(ielty) % dimensions
       pnode = element_type(ielty) % number_nodes
       porde = element_type(ielty) % order

       if( (.not. associated(element_type(ielty) % type_faces ).and. pface /= 0) ) allocate( element_type(ielty) % type_faces(max(1_ip,pface)) )
       if( (.not. associated(element_type(ielty) % node_faces ).and. pface /= 0) ) allocate( element_type(ielty) % node_faces(max(1_ip,pface)) )
       if( (.not. associated(element_type(ielty) % list_faces ).and. pface /= 0) ) allocate( element_type(ielty) % list_faces(max(1_ip,pnodf),max(1_ip,pface)) )
       if( (.not. associated(element_type(ielty) % type_edges ).and. pedge /= 0) ) allocate( element_type(ielty) % type_edges(pedge) )
       if( (.not. associated(element_type(ielty) % list_edges ).and. pedge /= 0) ) allocate( element_type(ielty) % list_edges(porde+1,pedge) )
       if(  .not. associated(element_type(ielty) % coord )    )                    allocate( element_type(ielty) % coord(pdime,pnode)  )
       
       if( pedge /= 0 ) element_type(ielty) % type_edges = 0

       if( ielty == BAR02 ) then
          !
          ! BAR02
          !
          element_type(ielty) % type_faces = type_faces_BAR02
          element_type(ielty) % node_faces = node_faces_BAR02
          element_type(ielty) % list_faces = list_faces_BAR02
          !element_type(ielty) % list_edges = list_edges_BAR02

       else if( ielty == BAR03 ) then
          !
          ! BAR03
          !
          element_type(ielty) % type_faces = type_faces_BAR03
          element_type(ielty) % node_faces = node_faces_BAR03
          element_type(ielty) % list_faces = list_faces_BAR03
          !element_type(ielty) % list_edges = list_edges_BAR03

       else if( ielty == BAR04 ) then
          !
          ! BAR04
          !
          element_type(ielty) % type_faces = type_faces_BAR04
          element_type(ielty) % node_faces = node_faces_BAR04
          element_type(ielty) % list_faces = list_faces_BAR04
          !element_type(ielty) % list_edges = list_edges_BAR04

       else if( ielty == TRI03 ) then
          !
          ! TRI03
          !
          element_type(ielty) % type_faces = type_faces_TRI03
          element_type(ielty) % node_faces = node_faces_TRI03
          element_type(ielty) % list_faces = list_faces_TRI03
          element_type(ielty) % type_edges = type_edges_TRI03
          element_type(ielty) % list_edges = list_edges_TRI03

       else if( ielty == QUA04 ) then
          !
          ! QUA04
          !
          element_type(ielty) % type_faces = type_faces_QUA04
          element_type(ielty) % node_faces = node_faces_QUA04
          element_type(ielty) % list_faces = list_faces_QUA04
          element_type(ielty) % type_edges = type_edges_QUA04
          element_type(ielty) % list_edges = list_edges_QUA04

       else if( ielty == TRI06 ) then
          !
          ! TRI06
          !
          element_type(ielty) % type_faces = type_faces_TRI06
          element_type(ielty) % node_faces = node_faces_TRI06
          element_type(ielty) % list_faces = list_faces_TRI06
          element_type(ielty) % type_edges = type_edges_TRI06
          element_type(ielty) % list_edges = list_edges_TRI06

       else if( ielty == QUA08 ) then
          !
          ! QUA08
          !
          element_type(ielty) % type_faces = type_faces_QUA08
          element_type(ielty) % node_faces = node_faces_QUA08
          element_type(ielty) % list_faces = list_faces_QUA08
          element_type(ielty) % type_edges = type_edges_QUA08
          element_type(ielty) % list_edges = list_edges_QUA08

       else if( ielty == QUA09 ) then
          !
          ! QUA09
          !
          element_type(ielty) % type_faces = type_faces_QUA09
          element_type(ielty) % node_faces = node_faces_QUA09
          element_type(ielty) % list_faces = list_faces_QUA09
          element_type(ielty) % type_edges = type_edges_QUA09
          element_type(ielty) % list_edges = list_edges_QUA09

       else if( ielty == TRI10 ) then
          !
          ! TRI10
          !
          element_type(ielty) % type_faces = type_faces_TRI10
          element_type(ielty) % node_faces = node_faces_TRI10
          element_type(ielty) % list_faces = list_faces_TRI10
          element_type(ielty) % list_edges = list_edges_TRI10

       else if ( ielty == QUA16 ) then
          !
          ! QUA16
          !
          element_type(ielty) % type_faces = type_faces_QUA16
          element_type(ielty) % node_faces = node_faces_QUA16
          element_type(ielty) % list_faces = list_faces_QUA16
          element_type(ielty) % type_edges = type_edges_QUA16
          element_type(ielty) % list_edges = list_edges_QUA16

       else if( ielty == TET04 ) then
          !
          ! TET04
          !
          element_type(ielty) % type_faces = type_faces_TET04
          element_type(ielty) % node_faces = node_faces_TET04
          element_type(ielty) % list_faces = list_faces_TET04
          element_type(ielty) % type_edges = type_edges_TET04
          element_type(ielty) % list_edges = list_edges_TET04

       else if( ielty == PYR05 ) then
          !
          ! PYR05
          !
          element_type(ielty) % type_faces = type_faces_PYR05
          element_type(ielty) % node_faces = node_faces_PYR05
          element_type(ielty) % list_faces = list_faces_PYR05
          element_type(ielty) % type_edges = type_edges_PYR05
          element_type(ielty) % list_edges = list_edges_PYR05

       else if( ielty == PEN06 ) then
          !
          ! PEN06
          !
          element_type(ielty) % type_faces = type_faces_PEN06
          element_type(ielty) % node_faces = node_faces_PEN06
          element_type(ielty) % list_faces = list_faces_PEN06
          element_type(ielty) % type_edges = type_edges_PEN06
          element_type(ielty) % list_edges = list_edges_PEN06

       else if( ielty == PEN15 ) then
          !
          ! PEN15
          !
          element_type(ielty) % type_faces = type_faces_PEN15
          element_type(ielty) % node_faces = node_faces_PEN15
          element_type(ielty) % list_faces = list_faces_PEN15
          !element_type(ielty) % list_edges = list_edges_PEN15

       else if( ielty == PEN18 ) then
          !
          ! PEN18
          !
          element_type(ielty) % type_faces = type_faces_PEN18
          element_type(ielty) % node_faces = node_faces_PEN18
          element_type(ielty) % list_faces = list_faces_PEN18
          !element_type(ielty) % list_edges = list_edges_PEN18

       else if( ielty == HEX08 ) then
          !
          ! HEX08
          !
          element_type(ielty) % type_faces = type_faces_HEX08
          element_type(ielty) % node_faces = node_faces_HEX08
          element_type(ielty) % list_faces = list_faces_HEX08
          element_type(ielty) % type_edges = type_edges_HEX08
          element_type(ielty) % list_edges = list_edges_HEX08

       else if( ielty == TET10 ) then
          !
          ! TET10
          !
          element_type(ielty) % type_faces = type_faces_TET10
          element_type(ielty) % node_faces = node_faces_TET10
          element_type(ielty) % list_faces = list_faces_TET10
          element_type(ielty) % type_edges = type_edges_TET10
          element_type(ielty) % list_edges = list_edges_TET10

       else if( ielty == HEX27 ) then
          !
          ! HEX27
          !
          element_type(ielty) % type_faces = type_faces_HEX27
          element_type(ielty) % node_faces = node_faces_HEX27
          element_type(ielty) % list_faces = list_faces_HEX27
          element_type(ielty) % type_edges = type_edges_HEX27
          element_type(ielty) % list_edges = list_edges_HEX27

       else if( ielty == HEX64 ) then
          !
          ! HEX64
          !
          element_type(ielty) % type_faces = type_faces_HEX64
          element_type(ielty) % node_faces = node_faces_HEX64
          element_type(ielty) % list_faces = list_faces_HEX64
          !element_type(ielty) % list_edges = list_edges_HEX64
          
       else if( ielty == TET20 ) then
          !
          ! TET20
          !
          element_type(ielty) % type_faces = type_faces_TET20
          element_type(ielty) % node_faces = node_faces_TET20
          element_type(ielty) % list_faces = list_faces_TET20
          !element_type(ielty) % list_edges = list_edges_TET20

       else if ( ielty == SHELL ) then
          !
          ! SHELL
          !
          element_type(ielty) % type_faces = type_faces_TRI03
          element_type(ielty) % node_faces = node_faces_TRI03
          element_type(ielty) % list_faces = list_faces_TRI03
          !element_type(ielty) % list_edges = list_edges_TRI03

       else if( ielty == BAR3D ) then
          !
          ! BAR3D
          !
          element_type(ielty) % type_faces = type_faces_BAR02
          element_type(ielty) % node_faces = node_faces_BAR02
          element_type(ielty) % list_faces = list_faces_BAR02
          !element_type(ielty) % list_edges = list_edges_BAR02

       else

          !call runend('FACES NOT PROGRAMMED')

       end if
       !
       ! Iso-parametric coordinates
       !
       call elmgeo_iso_parametric_element_coordinates(pdime,ielty,element_type(ielty) % coord)

    end do

  end subroutine elmgeo_element_type_initialization

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   Compute intersection between a segment and a plane
  !> @details Compute intersection between a segment and a plane
  !
  !-----------------------------------------------------------------------

  subroutine elmgeo_segpla(ndime,xplan,xcoo1,xcoo2,ifoun,plapo,toler_opt)
    implicit none
    integer(ip),           intent(in)  :: ndime                                      !< Dimension
    real(rp),              intent(in)  :: xplan(*)                                   !< Plane equation
    real(rp),              intent(in)  :: xcoo1(ndime)                               !< Segment first point
    real(rp),              intent(in)  :: xcoo2(ndime)                               !< Segment second point
    integer(ip),           intent(out) :: ifoun                                      !< If inside or not
    real(rp),    optional, intent(out) :: plapo(ndime)                               !< Intersection coordinates
    real(rp),    optional, intent(in)  :: toler_opt                                  !< Tolerance
    real(rp)                           :: t,toler

    ifoun = 0
    if( present(toler_opt) ) then
       toler = toler_opt
    else
       toler = epsilgeo_tol
    end if
    
    if( ndime == 2 ) then
       !
       ! 2D
       !
       t = xplan(1)*(xcoo2(1)-xcoo1(1)) + xplan(2)*(xcoo2(2)-xcoo1(2))

       if( t /= 0.0_rp ) then

          t = (- xplan(3) - xplan(1)*xcoo1(1) - xplan(2)*xcoo1(2) ) / t
          if( t >= toler .and. t <= 1.0_rp+toler ) then
             ifoun = 1
             if( present(plapo) ) then
                plapo(1) = t * (xcoo2(1)-xcoo1(1)) + xcoo1(1)
                plapo(2) = t * (xcoo2(2)-xcoo1(2)) + xcoo1(2)
             end if
          end if

       else if( abs( xplan(1) * xcoo1(1) + xplan(2) * xcoo1(2) + xplan(3) ) <= toler ) then

          ifoun = 1
          if( present(plapo) ) then
             plapo(1) = 0.5_rp * ( xcoo1(1) + xcoo2(1) )
             plapo(2) = 0.5_rp * ( xcoo1(2) + xcoo2(2) )
          end if

       end if

    else
       !
       ! Scalar product t = n.(P1,P2)
       !
       t = xplan(1)*(xcoo2(1)-xcoo1(1)) + xplan(2)*(xcoo2(2)-xcoo1(2)) + xplan(3)*(xcoo2(3)-xcoo1(3))

       if( t /= 0.0_rp ) then
          !
          ! Compute parametric coordinate t on P1-P2
          !
          t = (- xplan(4) - xplan(1)*xcoo1(1) - xplan(2)*xcoo1(2) - xplan(3)*xcoo1(3) ) / t
          if( t >= toler .and. t <= 1.0_rp+toler ) then
             ifoun = 1
             if( present(plapo) ) then
                plapo(1) = t * (xcoo2(1)-xcoo1(1)) + xcoo1(1)
                plapo(2) = t * (xcoo2(2)-xcoo1(2)) + xcoo1(2)
                plapo(3) = t * (xcoo2(3)-xcoo1(3)) + xcoo1(3)
             end if
          end if

       else if( abs( xplan(1) * xcoo1(1) + xplan(2) * xcoo1(2) + xplan(3) * xcoo1(3) + xplan(4) ) <= toler ) then
          !
          ! (P1,P2) is parallel to plane: check if P1 on plane
          !
          ifoun = 1
          if( present(plapo) ) then
             plapo(1) = 0.5_rp * ( xcoo1(1) + xcoo2(1) )
             plapo(2) = 0.5_rp * ( xcoo1(2) + xcoo2(2) )
             plapo(3) = 0.5_rp * ( xcoo1(3) + xcoo2(3) )
          end if
       end if
    end if

  end subroutine elmgeo_segpla

  !-----------------------------------------------------------------------
  !
  !> @brief   Compute intersection between a segment and a face
  !> @details Compute intersection between a segment and a face
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine elmgeo_segfac(ndime,pnodb,bocod,xcoo1,xcoo2,ifoun,plapo,toler_opt)
    implicit none
    integer(ip),           intent(in)  :: ndime                                      !< Dimension
    integer(ip),           intent(in)  :: pnodb                                      !< Number of boundary nodes
    real(rp),              intent(in)  :: bocod(ndime,pnodb)                         !< Boundary coordinates
    real(rp),              intent(in)  :: xcoo1(ndime)                               !< Segment first point
    real(rp),              intent(in)  :: xcoo2(ndime)                               !< Segment second point
    integer(ip),           intent(out) :: ifoun                                      !< If inside or not
    real(rp),              intent(out) :: plapo(ndime)                               !< Bandwidt
    real(rp),    optional, intent(in)  :: toler_opt                                  !< Tolerance
    integer(ip)                        :: idime,iboun,nboun,inodb
    logical(lg)                        :: kfl_qua04
    real(rp)                           :: boco2(ndime,3)
    real(rp)                           :: bouno(3),eps1,eps2
    real(rp)                           :: numer,denom,temdi,toler
    !
    ! Tolerance and initialization
    !
    if( present(toler_opt) ) then
       toler = toler_opt
    else
       toler = 1.0e-3_rp
    end if

    eps1  = -toler
    eps2  =  1.0_rp+toler
    ifoun =  0
    iboun =  0
    !
    ! Check if boundary is a QUA04 element
    !
    if( ndime == 3 .and. pnodb == 4 ) then
       kfl_qua04 = .true.
       nboun     = 2
    else
       kfl_qua04 = .false.
       nboun     = 1
       do inodb = 1,pnodb
          do idime = 1,ndime
             boco2(idime,inodb) = bocod(idime,inodb)
          end do
       end do
    end if

    do while( ifoun == 0 .and. iboun < nboun )
       iboun = iboun + 1
       if( kfl_qua04 ) then
          if( iboun == 1 ) then
             do idime = 1,3
                boco2(idime,1) = bocod(idime,1)
                boco2(idime,2) = bocod(idime,2)
                boco2(idime,3) = bocod(idime,4)
             end do
          else
             do idime = 1,3
                boco2(idime,1) = bocod(idime,2)
                boco2(idime,2) = bocod(idime,3)
                boco2(idime,3) = bocod(idime,4)
             end do
          end if
       end if
       !
       ! bouno: Exterior normal
       !

       call elmgeo_exttri(1_ip,ndime,pnodb,boco2,bouno)
       !
       ! Get characteristic dimensions
       !
       numer =  0.0_rp
       denom =  0.0_rp
       do idime = 1,ndime
          numer = numer + ( bouno(idime) * ( boco2(idime,1) - xcoo1(idime) )  )
          denom = denom + ( bouno(idime) * ( xcoo2(idime)   - xcoo1(idime) )  )
       end do

       if( denom /= 0.0_rp ) then
          !
          ! temdi: Normalized Distance measure from the first node
          !
          temdi = numer / denom
          !
          ! The segment intersects the plane that contains the iboun
          !
          if ( temdi >= eps1 .and. temdi <= eps2 ) then
             do idime = 1,ndime
                plapo(idime) = xcoo1(idime) + temdi * ( xcoo2(idime) - xcoo1(idime) )
             end do
             !
             ! Determine if the projection point on plane is inside the iboun
             !
             if( ndime == 3 ) then
                call elmgeo_instri(plapo,boco2,toler,ifoun)
             else if( ndime == 2 ) then
                call elmgeo_insbar(plapo,boco2,toler,ifoun)
             end if
          end if
       end if
    end do

  end subroutine elmgeo_segfac

  !-----------------------------------------------------------------------
  !
  !> @brief   Compute external normal to a boundary
  !> @details Compute external normal to a boundary
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine elmgeo_exttri(itask,ndime,pnodb,bocod,bouno)
    implicit none
    integer(ip), intent(in)  :: itask
    integer(ip), intent(in)  :: ndime
    integer(ip), intent(in)  :: pnodb
    real(rp),    intent(in)  :: bocod(ndime,pnodb)
    real(rp),    intent(out) :: bouno(ndime)
    integer(ip)              :: p1,p2,p3,idime
    real(rp)                 :: xfact,vec(3,3)

    if( ndime == 2 ) then

       if(itask == 1 ) then
          p1 = 1
          p2 = 2
       else
          p1 = 2
          p2 = 1
       end if

       vec(1,1) =  bocod(1,p2) - bocod(1,p1)
       vec(2,1) =  bocod(2,p2) - bocod(2,p1)
       bouno(1) =  vec(2,1)
       bouno(2) = -vec(1,1)

    else
       p1 = 1
       if( itask == 1 ) then
          p2 = 2
          p3 = 3
       else
          p2 = 3
          p3 = 2
       end if
       vec(1,1) = bocod(1,p2) - bocod(1,p1)
       vec(2,1) = bocod(2,p2) - bocod(2,p1)
       vec(3,1) = bocod(3,p2) - bocod(3,p1)
       vec(1,2) = bocod(1,p3) - bocod(1,p1)
       vec(2,2) = bocod(2,p3) - bocod(2,p1)
       vec(3,2) = bocod(3,p3) - bocod(3,p1)

       bouno(1) = vec(2,1) * vec(3,2) - vec(3,1) * vec(2,2)
       bouno(2) = vec(3,1) * vec(1,2) - vec(1,1) * vec(3,2)
       bouno(3) = vec(1,1) * vec(2,2) - vec(2,1) * vec(1,2)

    end if

    xfact = 0.0_rp
    do idime = 1,ndime
       xfact = xfact + bouno(idime) * bouno(idime)
    end do
    xfact = sqrt(xfact)
    if( xfact > 1.0e-12_rp ) then
       xfact = 1.0_rp / xfact
       do idime = 1,ndime
          bouno(idime) = xfact * bouno(idime)
       end do
    end if

  end subroutine elmgeo_exttri

  !-----------------------------------------------------------------------
  !
  !> @brief   Determine if a point is inside a triangle
  !> @details Determine if a point is inside a triangle using the same side technique
  !!          point: point facoo(dimension,vertices): triangle coordinates
  !!          ifoun: If is equal to 1, the point is inside the triangle, 0 otherwise
  !> @author  Cristobal Samaniego
  !
  !-----------------------------------------------------------------------

  subroutine elmgeo_instri(plapo,bocod,toler,ifoun)
    implicit none
    real(rp),    intent(in)  :: plapo(3)
    real(rp),    intent(in)  :: bocod(3,3)
    real(rp),    intent(in)  :: toler
    integer(ip), intent(out) :: ifoun
    real(rp)                 :: v0(3),v1(3),v2(3)
    real(rp)                 :: dot00,dot01,dot02,dot11,dot12
    real(rp)                 :: invDenom,bari1,bari2,xnorm
    !
    ! 3D
    !
    v0(1) = bocod(1,3) - bocod(1,1)
    v0(2) = bocod(2,3) - bocod(2,1)
    v0(3) = bocod(3,3) - bocod(3,1)
    xnorm = 1.0_rp/sqrt(dot_product(v0,v0))

    v1(1) = bocod(1,2) - bocod(1,1)
    v1(2) = bocod(2,2) - bocod(2,1)
    v1(3) = bocod(3,2) - bocod(3,1)

    v2(1) = plapo(1)   - bocod(1,1)
    v2(2) = plapo(2)   - bocod(2,1)
    v2(3) = plapo(3)   - bocod(3,1)
 
    v0    = v0 * xnorm
    v1    = v1 * xnorm
    v2    = v2 * xnorm
    
    dot00 = v0(1) * v0(1) + v0(2) * v0(2) + v0(3) * v0(3)
    dot01 = v0(1) * v1(1) + v0(2) * v1(2) + v0(3) * v1(3)
    dot02 = v0(1) * v2(1) + v0(2) * v2(2) + v0(3) * v2(3)
    dot11 = v1(1) * v1(1) + v1(2) * v1(2) + v1(3) * v1(3)
    dot12 = v1(1) * v2(1) + v1(2) * v2(2) + v1(3) * v2(3)
    !
    ! Compute barycentric coordinates
    !
    ifoun = 0
    if( abs(dot00 * dot11 - dot01 * dot01) > epsilgeo_alg) then
       invDenom = 1.0_rp / (dot00 * dot11 - dot01 * dot01)
       bari1    = (dot11 * dot02 - dot01 * dot12) * invDenom
       bari2    = (dot00 * dot12 - dot01 * dot02) * invDenom
       !
       ! Check
       !
       if( bari1 >= -toler .and. bari2 >= -toler .and. bari1 + bari2 <= 1.0_rp+toler ) then
          ifoun = 1
       end if
    end if

  end subroutine elmgeo_instri

  subroutine elmgeo_insbar(plapo,bocod,toler,ifoun)
    !-----------------------------------------------------------------------
    ! NAME
    !    segdis
    ! DESCRIPTION
    !    Minimun distance between a point and a segment
    !    plapo : point plapoinates
    !    coor1,coor2 : defines the segment
    !    ndime: dimension
    !    dista: distance
    !    proje: projection of the point on the segment
    !    ifoun = 1 the projection point is inside the segment
    !    ifoun = 0 the projection point is outside the segment
    ! USED BY
    !    pofadi
    !***
    !-----------------------------------------------------------------------
    implicit none
    real(rp),    intent(in)  :: plapo(2)
    real(rp),    intent(in)  :: bocod(2,2)
    real(rp),    intent(in)  :: toler
    integer(ip), intent(out) :: ifoun
    integer(ip)              :: idime
    real(rp)                 :: numer,denom,dsegm,eps1,eps2

    eps1  = -toler
    eps2  = 1.0_rp + toler
    numer = 0.0_rp
    denom = 0.0_rp
    do idime = 1,2
       numer = numer + (bocod(idime,2) - bocod(idime,1)) * (plapo(idime)   - bocod(idime,1))
       denom = denom + (bocod(idime,2) - bocod(idime,1)) * (bocod(idime,2) - bocod(idime,1))
    end do

    dsegm = numer / denom

    if( dsegm < eps1 ) then

       ifoun = 0_ip

    else if ( dsegm >= eps1 .and. dsegm <= eps2 ) then

       ifoun = 1_ip

    else

       ifoun = 0_ip

    end if

  end subroutine elmgeo_insbar

  !-----------------------------------------------------------------------
  !
  !> @brief   Compute the quality of a TET04 element
  !> @details Determine the quality based on gamma-team (INRIA) criterion
  !>          P.L. George, Improvement on Delaunay based 3D automatic mesh generator,
  !>          Finite Elements in Analysis and Design, Vol. 25, pp. 297--317, 1997
  !> @date    26/02/2013
  !> @author  Beatriz Eguzkitza
  !
  !-----------------------------------------------------------------------

  subroutine elmgeo_quatet_gamma(elcod,Q,signi)
    implicit none
    real(rp),    intent(in)  :: elcod(3,4)
    real(rp),    intent(out) :: Q
    integer(ip), intent(out) :: signi
    integer(ip)              :: idime
    real(rp)                 :: hmax,alpha,volum,sur1,sur2,sur3,sur4,surf
    real(rp)                 :: side1,side2,side3,side4,side5,side6,hbase,aux
    real(rp)                 :: cora1(3),cora2(3),cora4(3)
    real(rp)                 :: corb1(3),corb2(3),corb3(3)
    !
    ! Edges vectors
    !
    do idime = 1,3
       cora1(idime) = elcod(idime,1) - elcod(idime,2)  ! 2 -> 1
       corb1(idime) = elcod(idime,4) - elcod(idime,2)  ! 2 -> 4
       cora2(idime) = elcod(idime,4) - elcod(idime,1)  ! 1 -> 4
       corb2(idime) = elcod(idime,3) - elcod(idime,1)  ! 1 -> 3
       corb3(idime) = elcod(idime,3) - elcod(idime,2)  ! 2 -> 3
       cora4(idime) = elcod(idime,4) - elcod(idime,3)  ! 3 -> 4
    end do
    !
    ! hmax
    !
    side1 = 0.0_rp
    side2 = 0.0_rp
    side3 = 0.0_rp
    side4 = 0.0_rp
    side5 = 0.0_rp
    side6 = 0.0_rp

    do idime = 1,3
       side1 = side1 + cora1(idime) * cora1(idime)
       side2 = side2 + corb1(idime) * corb1(idime)
       side3 = side3 + cora2(idime) * cora2(idime)
       side4 = side4 + corb2(idime) * corb2(idime)
       side5 = side5 + corb3(idime) * corb3(idime)
       side6 = side6 + cora4(idime) * cora4(idime)
    end do

    hbase = min(side1,side3,side4)

    side1 = sqrt(side1)
    side2 = sqrt(side2)
    side3 = sqrt(side3)
    side4 = sqrt(side4)
    side5 = sqrt(side5)
    side6 = sqrt(side6)
    hmax  = max(side1,side2,side3,side4,side5,side6)
    !
    ! SURF= Total surface
    !
    sur1 =    ( cora1(1) * corb1(2) - cora1(2) * corb1(1)) * (cora1(1) * corb1(2) - cora1(2) * corb1(1) ) &
         &  + ( cora1(3) * corb1(1) - cora1(1) * corb1(3)) * (cora1(3) * corb1(1) - cora1(1) * corb1(3) ) &
         &  + ( cora1(2) * corb1(3) - cora1(3) * corb1(2)) * (cora1(2) * corb1(3) - cora1(3) * corb1(2) )

    sur2 =    ( cora2(1) * corb2(2) - corb2(1) * cora2(2)) * (cora2(1) * corb2(2) - corb2(1) * cora2(2) ) &
         &  + ( cora2(3) * corb2(1) - cora2(1) * corb2(3)) * (cora2(3) * corb2(1) - cora2(1) * corb2(3) ) &
         &  + ( cora2(2) * corb2(3) - cora2(3) * corb2(2)) * (cora2(2) * corb2(3) - cora2(3) * corb2(2) )

    sur3 =    ( cora1(1) * corb3(2) - corb3(1) * cora1(2)) * (cora1(1) * corb3(2) - corb3(1) * cora1(2) ) &
         &  + ( cora1(3) * corb3(1) - cora1(1) * corb3(3)) * (cora1(3) * corb3(1) - cora1(1) * corb3(3) ) &
         &  + ( cora1(2) * corb3(3) - cora1(3) * corb3(2)) * (cora1(2) * corb3(3) - cora1(3) * corb3(2) )

    sur4 =    ( cora4(1) * corb1(2) - corb1(1) * cora4(2)) * (cora4(1) * corb1(2) - corb1(1) * cora4(2) ) &
         &  + ( cora4(3) * corb1(1) - cora4(1) * corb1(3)) * (cora4(3) * corb1(1) - cora4(1) * corb1(3) ) &
         &  + ( cora4(2) * corb1(3) - cora4(3) * corb1(2)) * (cora4(2) * corb1(3) - cora4(3) * corb1(2) )

    surf = 0.5_rp * ( sqrt(sur1) + sqrt(sur2) + sqrt(sur3) + sqrt(sur4) )
    !
    ! VOLUM= Volume and SIGNI= signi
    !
    volum = (  cora1(1) * cora2(2) * corb2(3) + cora1(2) * cora2(3) * corb2(1) + cora1(3) * cora2(1) * corb2(2) &
         &    -cora1(3) * cora2(2) * corb2(1) - cora1(2) * cora2(1) * corb2(3) - cora1(1) * cora2(3) * corb2(2) ) / 6.0_rp
    if( volum <= 0.0_rp ) then
       signi = -1_ip
    else
       signi  = 1_ip
    end if
    !
    ! Q= Quality
    !
    alpha = sqrt(6.0_rp)/12.0_rp
    aux = alpha * hmax * surf
    Q     = aux / ( 3.0_rp * abs(volum) + epsilgeo_div * (aux+epsilgeo_div))

  end subroutine elmgeo_quatet_gamma

  !-----------------------------------------------------------------------
  !
  !> @brief   Compute the quality and aspect ratio of a TET04 element
  !> @details Determine the quality based on coondition number kappaS and
  !>          aspect ratio asrad. See Knupp criterion:
  !>          P.Knupp, "Algebraic Mesh Quality Metrics," SIAM J. Sci.
  !>          Comput., Vol. 23, No. 1, pp193-218, 2001.
  !> @date    26/02/2013
  !> @author  Beatriz Eguzkitza
  !
  !-----------------------------------------------------------------------

  subroutine elmgeo_quatet_kappar(ndime,elcod,kappaS,asrad)
    implicit none
    integer(ip), intent(in)           :: ndime
    real(rp),    intent(in)           :: elcod(ndime,*)
    real(rp),    intent(out)          :: kappaS
    real(rp),    intent(out)          :: asrad
    integer(ip)                       :: idime,inode,ii,jj,pnode
    real(rp)                          :: invW(ndime,3),tinvW(3,ndime),invWt(3,ndime)
    real(rp)                          :: A(ndime,3),S(ndime,3),invS(ndime,3)
    real(rp)                          :: normS,norminvS,StS(ndime,3)
    real(rp)                          :: aux(ndime,3),aux2(ndime,3)
    real(rp)                          :: detinvWt,detS,volum
    real(rp)                          :: detjm
    real(rp)                          :: gpcar(ndime,4),sumsi,t1,t2,t3
    real(rp)                          :: side1(ndime),side2(ndime),side3(ndime)
    real(rp)                          :: side4(ndime),side5(ndime),side6(ndime)
    !
    ! Condition number kappaS
    !
    detinvWt  =  0.0_rp
    detS      =  0.0_rp
    if( ndime == 3 ) then
       pnode     =  4
       invW(1,1) =  1.0_rp
       invW(1,2) = -1.0_rp/3.0_rp * sqrt(3.0_rp)         ! -5.77350269189626e-01
       invW(1,3) = -1.0_rp/6.0_rp * sqrt(6.0_rp)         ! -4.08248290463863e-01
       invW(2,1) =  0.0_rp
       invW(2,2) =  2.0_rp/3.0_rp * sqrt(3.0_rp)         !  1.15470053837925e+00
       invW(2,3) = -1.0_rp/6.0_rp * sqrt(6.0_rp)         ! -4.08248290463863e-01
       invW(3,1) =  0.0_rp
       invW(3,2) =  0.0_rp
       invW(3,3) =  1.0_rp/2.0_rp * sqrt(6.0_rp)         !  1.22474487139159e+00
    else
       pnode     =  3
       invW(1,1) =  1.0_rp
       invW(1,2) =  1.0_rp/2.0_rp
       invW(2,1) =  0.0_rp
       invW(2,2) =  sqrt(3.0_rp)/2.0_rp                  !  1.15470053837925e+00
    end if

    normS    =  0.0_rp
    norminvS =  0.0_rp

    do inode = 1,pnode-1
       do idime = 1,ndime
          A(idime,inode)    = 0.0_rp
          S(idime,inode)    = 0.0_rp
          StS (idime,inode) = 0.0_rp
          invS(idime,inode) = 0.0_rp
       end do
    end do
    do inode = 1,pnode-1
       do idime = 1,ndime
          tinvW (idime,inode) = invW(inode,idime)
       end do
    end do

    call maths_invert_matrix(ndime,tinvW,detinvWt,invWt)
    !call invmtx(tinvW,invWt,detinvWt,ndime)

    ii = 0
    do inode = 2,pnode
       ii = ii + 1
       do idime = 1,ndime
          A(idime,ii) = elcod(idime,inode) - elcod(idime,1)
       end do
    end do

    do inode = 1,pnode-1
       do idime = 1,ndime
          do jj = 1,pnode-1
             S(idime,inode) = S(idime,inode) + A(idime,jj) * invW(jj,inode)
          end do
       end do
    end do

    do inode = 1,pnode-1
       do idime = 1,ndime
          normS = normS + S(idime,inode) * S(idime,inode)
          do jj = 1,pnode-1
             StS (idime,inode) = StS(idime,inode) + S(jj,idime) * S(jj,inode)
          end do
       end do
    end do
    normS = sqrt(normS)

    do inode = 1,pnode-1
       do idime = 1,ndime
          aux(idime,inode) = S(idime,inode)
       end do
    end do
    call maths_invert_matrix(ndime,aux,detS,aux2)
    !call invmtx(aux,aux2,detS,ndime)
    do inode = 1,pnode-1
       do idime=1,ndime
          invS(idime,inode) = aux2(idime,inode)
       end do
    end do

    do inode = 1,pnode-1
       do idime=1,ndime
          norminvS = norminvS + invS(idime,inode) * invS(idime,inode)
       end do
    end do

    norminvS = sqrt(norminvS)
    kappaS   = normS *  norminvS
    kappaS   = kappaS /3.0_rp     ! So that ideal element has kappaS=1
    !
    ! VOLUM of the element
    !
    if( ndime == 2 ) then
       detjm      =  (-elcod(1,1)+elcod(1,2))*(-elcod(2,1)+elcod(2,3)) &
            &       -(-elcod(2,1)+elcod(2,2))*(-elcod(1,1)+elcod(1,3))
       volum      =  0.5_rp / detjm
    else
       gpcar(1,1) =  elcod(1,2) - elcod(1,1)
       gpcar(1,2) =  elcod(1,3) - elcod(1,1)
       gpcar(1,3) =  elcod(1,4) - elcod(1,1)
       gpcar(2,1) =  elcod(2,2) - elcod(2,1)
       gpcar(2,2) =  elcod(2,3) - elcod(2,1)
       gpcar(2,3) =  elcod(2,4) - elcod(2,1)
       gpcar(3,1) =  elcod(3,2) - elcod(3,1)
       gpcar(3,2) =  elcod(3,3) - elcod(3,1)
       gpcar(3,3) =  elcod(3,4) - elcod(3,1)
       t1         =  gpcar(2,2) * gpcar(3,3) - gpcar(3,2) * gpcar(2,3)
       t2         = -gpcar(2,1) * gpcar(3,3) + gpcar(3,1) * gpcar(2,3)
       t3         =  gpcar(2,1) * gpcar(3,2) - gpcar(3,1) * gpcar(2,2)
       detjm      =  gpcar(1,1) * t1 + gpcar(1,2) * t2 + gpcar(1,3) * t3
       volum      =  1.0_rp / ( 6.0_rp * detjm )
    end if

    if( ndime == 3 ) then
       do idime = 1,ndime
          side1(idime) = elcod(idime,4) - elcod(idime,1)
          side2(idime) = elcod(idime,3) - elcod(idime,1)
          side3(idime) = elcod(idime,2) - elcod(idime,1)
          side4(idime) = elcod(idime,4) - elcod(idime,2)
          side5(idime) = elcod(idime,3) - elcod(idime,2)
          side6(idime) = elcod(idime,4) - elcod(idime,3)
       end do
    else
       do idime = 1,ndime
          side1(idime) = elcod(idime,3)-elcod(idime,1)
          side2(idime) = elcod(idime,2)-elcod(idime,1)
          side3(idime) = elcod(idime,3)-elcod(idime,2)
          side4(idime) = 0.0_rp
          side5(idime) = 0.0_rp
          side6(idime) = 0.0_rp
       end do
    end if

    sumsi = 0.0_rp
    do idime = 1,ndime
       sumsi = sumsi + &
            side1(idime) * side1(idime) + side2(idime) * side2(idime) + &
            side3(idime) * side3(idime) + side4(idime) * side4(idime) + &
            side5(idime) * side5(idime) + side6(idime) * side6(idime)
    end do
    !
    ! Aspect ratio
    !
    asrad = sumsi / ( volum + epsilgeo_div * (sumsi + epsilgeo_div) )

  end subroutine elmgeo_quatet_kappar

  !-----------------------------------------------------------------------
  !> @author  Guillaume Houzeaux
  !> @date    18/10/2012
  !> @brief   Computes the boundary normals
  !> @details Comppute the boundary normals. Only valid for
  !>          TRI03, TRI06, QUA04, QUA08, QUA09 boundaries.
  !>          For QUA04 and QUA09, thye exterior normal is averaged
  !>          over the two triangles forming the boundary element.
  !>          According to the mesher used, we have two options:\n
  !>          - Clock wise numbering: normal points inwards
  !>          - Counterclock wise numbering: normal points outwards
  !>            - GiD: normal points inwards
  !>            - Fensap: normal points outwards
  !>            - Windmesh: normal points outwards
  !>          \verbatim
  !>          - BOUNO(1:NDIME,IBOUN) .... Exterior normal to IBOUN
  !>          - NINVE ................... Number of inverted normals
  !>                                      that were pointing inwards
  !>          \endverbatim
  !-----------------------------------------------------------------------

  subroutine elmgeo_bounor(&
       ichek,kboun,ndime,mnodb,mnode,lnodb,ltypb,&
       lelbo,ltype,lnods,nnode,coord,ninve,bouno)
    implicit none
    integer(ip), intent(in)            :: ichek                          !< Check or don't check
    integer(ip), intent(in)            :: kboun                          !< Number of boundaries
    integer(ip), intent(in)            :: ndime                          !< Space dimension
    integer(ip), intent(in)            :: mnodb                          !< Max number of nodes per boundary
    integer(ip), intent(in)            :: mnode                          !< Max number of nodes per element
    integer(ip), intent(in)            :: lnodb(mnodb,kboun)             !< Boundary connectivity
    integer(ip), intent(in)            :: ltypb(kboun)                   !< Boundary type
    integer(ip), intent(in), optional  :: lelbo(kboun)                   !< Boundary/element connectvity
    integer(ip), intent(in)            :: ltype(*)                       !< Element type
    integer(ip), intent(in)            :: nnode(*)                       !< Element type number of node
    integer(ip), intent(in)            :: lnods(mnode,*)                 !< Element connectivity
    real(rp),    intent(in)            :: coord(ndime,*)                 !< Node coordiantes
    integer(ip), intent(inout)         :: ninve                          !< Number of inverted boundaries
    real(rp),    intent(out)           :: bouno(ndime,kboun)             !< Boundary exterior normals
    integer(ip)                        :: iboun,p1,p2,p3,idime
    integer(ip)                        :: ipoin,inode,inodb,ielem
    integer(ip)                        :: pnodb,pblty,pnode
    real(rp)                           :: elcod(ndime,mnode)
    real(rp)                           :: bocod(ndime,mnodb)
    real(rp)                           :: vec(3,3),rmod

    if( ndime == 2 ) then

       !-------------------------------------------------------------------
       !
       ! 2D case
       !
       !-------------------------------------------------------------------

       do iboun = 1,kboun
          p1             =  lnodb(1,iboun)
          p2             =  lnodb(2,iboun)
          vec(1,1)       =  coord(1,p2) - coord(1,p1)
          vec(2,1)       =  coord(2,p2) - coord(2,p1)
          bouno(1,iboun) =  vec(2,1)
          bouno(2,iboun) = -vec(1,1)
       end do

    else if( ndime == 3 ) then

       !-------------------------------------------------------------------
       !
       ! 3D case
       !
       !-------------------------------------------------------------------

       do iboun = 1,kboun
          pblty = abs(ltypb(iboun))

          if( pblty == TRI03 .or.  pblty == TRI06 ) then

             p1                 = lnodb(1,iboun)
             p2                 = lnodb(2,iboun)
             p3                 = lnodb(3,iboun)
             call elmgeo_triangle_normal(p1,p2,p3,coord,vec,ndime)
             bouno(    1,iboun) = 0.5_rp*vec(1,3)
             bouno(    2,iboun) = 0.5_rp*vec(2,3)
             bouno(ndime,iboun) = 0.5_rp*vec(3,3)

          else if( pblty == QUA04 .or. pblty == QUA09 ) then

             p1             = lnodb(1,iboun)
             p2             = lnodb(2,iboun)
             p3             = lnodb(3,iboun)
             call elmgeo_triangle_normal(p1,p2,p3,coord,vec,ndime)
             bouno(1,iboun) = vec(1,3)
             bouno(2,iboun) = vec(2,3)
             bouno(3,iboun) = vec(3,3)

             p1             = lnodb(1,iboun)
             p2             = lnodb(3,iboun)
             p3             = lnodb(4,iboun)
             call elmgeo_triangle_normal(p1,p2,p3,coord,vec,ndime)
             bouno(1,iboun) = bouno(1,iboun) + vec(1,3)
             bouno(2,iboun) = bouno(2,iboun) + vec(2,3)
             bouno(3,iboun) = bouno(3,iboun) + vec(3,3)
             !
             ! HEX: We assume that ordering points towards inside
             !
             rmod = sqrt(   bouno(1,iboun) * bouno(1,iboun) &
                  &      + bouno(2,iboun) * bouno(2,iboun) &
                  &      + bouno(3,iboun) * bouno(3,iboun) )
             bouno(1,iboun) = bouno(1,iboun) / rmod
             bouno(2,iboun) = bouno(2,iboun) / rmod
             bouno(3,iboun) = bouno(3,iboun) / rmod

          end if

       end do
    end if
    !
    ! Ensure normal points outwards and normalize it
    !
    if( ichek == 1 ) then
       do iboun = 1,kboun
          pblty = abs(ltypb(iboun))
          pnodb = nnode(pblty)
          ielem = lelbo(iboun)
          pnode = nnode(ltype(ielem))
          do inodb = 1,pnodb
             ipoin = lnodb(inodb,iboun)
             do idime = 1,ndime
                bocod(idime,inodb) = coord(idime,ipoin)
             end do
          end do
          do inode = 1,pnode
             ipoin = lnods(inode,ielem)
             do idime = 1,ndime
                elcod(idime,inode) = coord(idime,ipoin)
             end do
          end do
          call elmgeo_check_boundary_orientation(ndime,pnode,pnodb,bouno(1,iboun),bocod,elcod,ninve)
          rmod = dot_product(bouno(1:ndime,iboun),bouno(1:ndime,iboun))
          if( rmod /= 0.0_rp ) bouno(1:ndime,iboun) = bouno(1:ndime,iboun) / sqrt(rmod)
       end do
    end if

  end subroutine elmgeo_bounor

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   Natural coordinates
  !> @details Evaluate the natural coordinate of a point in an element
  !>          or in a boundary element.
  !>          - Point COLOC(3) is in an element
  !>          - Point COLOC(3) is in a boundary element. In this case it
  !>            it should be on the curve
  !>          INPUT
  !>             NDIME ... Dimension
  !>             PTOPO ... Element topology
  !>                       -1 Bar
  !>                        0 Quadrilateral  Hexahedra
  !>                        1 Triangle       Tetrahedra
  !>                        2    -           Pentahedra (wedge-shaped)
  !>                        3    -           Pyramid
  !>             PNODE ... Number of element nodes
  !>             ELCOD ... Element node coordinates
  !>             COGLO ... Global coordinates of test point
  !>             LMINI ... Minimum local coordinate (0.0)
  !>             LMAXI ... Maximum local coordinate (1.0)
  !>          OUTPUT
  !>             IFOUN ... 1 if point is in element
  !>                       0 otherwise
  !>             COLOC ... Local coordinates of test point
  !>             SHAPF ... Shape function of test point in element
  !>             DERIV ... Shape function derivatives of test point in element
  !
  !-----------------------------------------------------------------------

  subroutine elmgeo_natural_coordinates(  &
       ndime,pelty,pnode,elcod,shapf,     &
       deriv,coglo,coloc,ifoun,toler_opt, &
       BOUNDING_BOX                       )
    
    integer(ip), intent(in)           :: ndime              !< Problem dimension
    integer(ip), intent(in)           :: pelty              !< Element type
    integer(ip), intent(in)           :: pnode              !< Element number of nodes
    real(rp),    intent(in)           :: coglo(ndime)       !< Test point coordinates
    real(rp),    intent(in)           :: elcod(ndime,pnode) !< Element coordinates
    integer(ip), intent(out)          :: ifoun              !< Test point is in/out element
    real(rp),    intent(out)          :: coloc(3)           !< Local coordinates of test point in element
    real(rp),    intent(out)          :: deriv(ndime,pnode) !< Shape funcitons derivatives of test point
    real(rp),    intent(out)          :: shapf(pnode)       !< Shape functions of test point
    real(rp),    intent(in), optional :: toler_opt          !< Tolerance rto decide it's in/ it's out
    logical(lg), intent(in), optional :: BOUNDING_BOX       !> If bounding box tests should be carried out
    real(rp)                          :: toler,lmini,lmaxi,tole2
    real(rp)                          :: time1,time2
    logical(lg)                       :: do_bounding_box_test
    !
    ! Define tolerance
    !
    if( present(toler_opt) ) then
       toler = toler_opt
    else
       toler = epsilgeo_def
    end if    
    lmini = - toler
    lmaxi = 1.0_rp + toler
    ifoun = 0
    !
    ! Bounding box test
    !
    if( present(BOUNDING_BOX) ) then
       do_bounding_box_test = BOUNDING_BOX
    else
       do_bounding_box_test = .true.
    end if
    
    if( do_bounding_box_test ) then
       do_bounding_box_test = elmgeo_inside_element_bounding_box(ndime,pnode,elcod,coglo,toler)
    else
       do_bounding_box_test = .true.
    end if
       
    if( do_bounding_box_test ) then
       !
       ! Point COGLO is in element bounding box
       !
       if( pelty < 0 ) then
          !
          ! Hole element
          !
          call runend('ELMGEO_NATURAL_COORDIANTES: DO NOT KNOW WHAT TO DO')

       else if( pelty == BAR02 ) then
          !
          ! BAR02
          !
          call cputim(time1) 
          call elmgeo_inside_BAR02(&
               pnode,lmaxi,elcod,coglo,coloc,shapf,deriv,ifoun)
          call cputim(time2)
          

       else if( ( pelty == TRI03 .or. pelty == QUA04 ) .and. ndime == 2 ) then
          !
          ! TRI03 and QUA04
          !
          call cputim(time1) 
          call elmgeo_inside_TRI03_QUA04(&
               pnode,lmini,lmaxi,elcod,coglo,coloc,shapf,deriv,ifoun)
          call cputim(time2)

       else if( pelty == QUA09 .and. ndime == 2 ) then
          !
          ! QUA09 - same restriction as in hex27
          !
          !tole2 = 10.0_rp * toler
          !if ( elmgeo_norm2(ndime,elcod(:, 5) - 0.50_rp * ( elcod(:,1) + elcod(:,2) ) ) > tole2 .or. &
          !     elmgeo_norm2(ndime,elcod(:, 6) - 0.50_rp * ( elcod(:,2) + elcod(:,3) ) ) > tole2 .or. &
          !     elmgeo_norm2(ndime,elcod(:, 7) - 0.50_rp * ( elcod(:,3) + elcod(:,4) ) ) > tole2 .or. &
          !     elmgeo_norm2(ndime,elcod(:, 8) - 0.50_rp * ( elcod(:,4) + elcod(:,1) ) ) > tole2 .or. &
          !     elmgeo_norm2(ndime,elcod(:, 9) - 0.25_rp * ( elcod(:,1) + elcod(:,2) + elcod(:,3) + elcod(:,4) ) ) > tole2  ) then
          !   write(*,*) 'elcod',elcod
          !   call runend('elmgeo_natural_coordinates: QUA09 not ready for elements without straight faces')
          !end if
          write(*,*) 'elmgeo_natural_coordinates: QUA09 not ready for elements without straight faces'
          !ifoun = 1
          coloc = 0.0_rp
          !call element_grid_test(&
          !     ndime,pnode,pelty,elcod,coglo,coloc)
          call elmgeo_newrap(&
               coglo,coloc,ndime,pnode,elcod,shapf,deriv,ifoun,&
               pelty,USER_INITIAL_GUESS=.true.)

          !print*,'coloc b=',coloc(1:2)
          !call elmgeo_inside_TRI03_QUA04(&
          !    4_ip,lmini,lmaxi,elcod,coglo,coloc,shapf,deriv,ifoun)

       else if( pelty == TET04 ) then
          !
          ! TET04
          !
          call cputim(time1)
          call elmgeo_inside_TET04(&
               lmini,lmaxi,elcod,coglo,coloc,&
               shapf,deriv,ifoun)
          call cputim(time2)

       else if( pelty == HEX27 ) then
          !
          ! HEX27 - Guillaume said using elmgeo_inside_element_using_faces
          ! He suggested usin newrap directly. I belive it will not work for the
          ! general case. I will do it as if it were an HEX08 if faces are 'flat'
          !
          ! Check that the nodes corresponding to quadratic are actually in the midlle
          ! Else runend
          !
          tole2 = 10.0_rp * toler
          !
          !call runend("ERROR EN BG, 'norm2' esto es estandar hasta Fortran 2008")
          !
          if ( elmgeo_norm2 (ndime,elcod(:, 9) - 0.5_rp * ( elcod(:,1) + elcod(:,2) ) ) > tole2 .or. &
               elmgeo_norm2 (ndime,elcod(:,10) - 0.5_rp * ( elcod(:,2) + elcod(:,3) ) ) > tole2 .or. &
               elmgeo_norm2 (ndime,elcod(:,11) - 0.5_rp * ( elcod(:,3) + elcod(:,4) ) ) > tole2 .or. &
               elmgeo_norm2 (ndime,elcod(:,12) - 0.5_rp * ( elcod(:,4) + elcod(:,1) ) ) > tole2 .or. &
               elmgeo_norm2 (ndime,elcod(:,13) - 0.5_rp * ( elcod(:,1) + elcod(:,5) ) ) > tole2 .or. &
               elmgeo_norm2 (ndime,elcod(:,14) - 0.5_rp * ( elcod(:,2) + elcod(:,6) ) ) > tole2 .or. &
               elmgeo_norm2 (ndime,elcod(:,15) - 0.5_rp * ( elcod(:,3) + elcod(:,7) ) ) > tole2 .or. &
               elmgeo_norm2 (ndime,elcod(:,16) - 0.5_rp * ( elcod(:,4) + elcod(:,8) ) ) > tole2 .or. &
               elmgeo_norm2 (ndime,elcod(:,17) - 0.5_rp * ( elcod(:,5) + elcod(:,6) ) ) > tole2 .or. &
               elmgeo_norm2 (ndime,elcod(:,18) - 0.5_rp * ( elcod(:,6) + elcod(:,7) ) ) > tole2 .or. &
               elmgeo_norm2 (ndime,elcod(:,19) - 0.5_rp * ( elcod(:,7) + elcod(:,8) ) ) > tole2 .or. &
               elmgeo_norm2 (ndime,elcod(:,20) - 0.5_rp * ( elcod(:,8) + elcod(:,5) ) ) > tole2 .or. &
               elmgeo_norm2 (ndime,elcod(:,21) - 0.25_rp *( elcod(:,1) + elcod(:,2) + elcod(:,3) + elcod(:,4) ) ) > tole2 .or. &
               elmgeo_norm2 (ndime,elcod(:,22) - 0.25_rp *( elcod(:,1) + elcod(:,2) + elcod(:,5) + elcod(:,6) ) ) > tole2 .or. &
               elmgeo_norm2 (ndime,elcod(:,23) - 0.25_rp *( elcod(:,2) + elcod(:,3) + elcod(:,6) + elcod(:,7) ) ) > tole2 .or. &
               elmgeo_norm2 (ndime,elcod(:,24) - 0.25_rp *( elcod(:,3) + elcod(:,4) + elcod(:,7) + elcod(:,8) ) ) > tole2 .or. &
               elmgeo_norm2 (ndime,elcod(:,25) - 0.25_rp *( elcod(:,1) + elcod(:,4) + elcod(:,5) + elcod(:,8) ) ) > tole2 .or. &
               elmgeo_norm2 (ndime,elcod(:,26) - 0.25_rp *( elcod(:,5) + elcod(:,6) + elcod(:,7) + elcod(:,8) ) ) > tole2 .or. &
               elmgeo_norm2 (ndime,elcod(:,27) - 0.125_rp*( elcod(:,1) + elcod(:,2) + elcod(:,3) + elcod(:,4) +                &
               elcod(:,5) + elcod(:,6) + elcod(:,7) + elcod(:,8)) ) > tole2) then
             write(*,*) 'elcod',elcod
             call runend('elmgeo_natural_coordinates: HEX27 not ready for elements without straight faces')
          end if
          !
          ! Compute local from global coordinates of an HEX08 - Change pelty to HEX08 & pnode to 8
          !
          if( elmgeo_inside_element_using_faces(HEX08,elcod,coglo) ) then
             call cputim(time1)
             call elmgeo_newrap(coglo,coloc,ndime,8_ip,elcod,shapf,deriv,ifoun,pelty)
             ifoun = 1
             call cputim(time2)

          else

             ifoun = 0
             return

          end if

       else
          !
          ! Compute local from global coordinates
          !
          if( elmgeo_inside_element_using_faces(pelty,elcod,coglo,toler) ) then
             
             call cputim(time1)
             call elmgeo_newrap(coglo,coloc,ndime,pnode,elcod,shapf,deriv,ifoun,pelty)
             call cputim(time2)

          else

             ifoun = 0
             return

          end if

       end if

    end if

  end subroutine elmgeo_natural_coordinates

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   Natural coordinates
  !> @details Check if point with global coordinates (x,y,z)=COGLO is inside
  !>          a TET04 element
  !>
  !>          +-                    -+ +-  -+   +-    -+
  !>          | -x1+x2 -x1+x3 -x1+x4 | | s1 |   | x-x1 |
  !>          | -y1+y2 -y1+y3 -y1+y3 | | s2 | = | y-y1 |
  !>          | -z1+z2 -z1+z3 -z1+z4 | | s3 |   | z-z1 |
  !>          +-                    -+ +-  -+   +-    -+
  !
  !-----------------------------------------------------------------------

  pure subroutine elmgeo_inside_TET04(&
       lmini,lmaxi,elcod,coglo,coloc,&
       shapf,deriv,ifoun)
    real(rp),    intent(in)  :: lmini
    real(rp),    intent(in)  :: lmaxi
    real(rp),    intent(in)  :: elcod(3,4)
    real(rp),    intent(in)  :: coglo(3)
    integer(ip), intent(out) :: ifoun
    real(rp),    intent(out) :: coloc(3)
    real(rp),    intent(out) :: shapf(4)
    real(rp),    intent(out) :: deriv(3,4)
    real(rp)                 :: deter,xjaci(3,3),xjacm(3,3)
    real(rp)                 :: rhsid(3),ezzzt,denom,t1,t2,t3
    
    real(rp) :: L
!    real(rp) :: t(3)
    
    ! t to center at origin, L to reescale first edge to unit length
    !t(:) = elcod(:,1)
    L    = sqrt( (elcod(1,2)-elcod(1,1))**2 + (elcod(2,2)-elcod(2,1))**2 + (elcod(3,2)-elcod(3,1))**2 )
    !L = 1_rp
    
    ifoun      = 0
    xjacm(1,1) = -elcod(1,1) + elcod(1,2)
    xjacm(2,1) = -elcod(2,1) + elcod(2,2)
    xjacm(3,1) = -elcod(3,1) + elcod(3,2)
    xjacm(1,2) = -elcod(1,1) + elcod(1,3)
    xjacm(2,2) = -elcod(2,1) + elcod(2,3)
    xjacm(3,2) = -elcod(3,1) + elcod(3,3)
    xjacm(1,3) = -elcod(1,1) + elcod(1,4)
    xjacm(2,3) = -elcod(2,1) + elcod(2,4)
    xjacm(3,3) = -elcod(3,1) + elcod(3,4)
    xjacm = xjacm/L

    t1         =  xjacm(2,2) * xjacm(3,3) - xjacm(3,2) * xjacm(2,3)
    t2         = -xjacm(2,1) * xjacm(3,3) + xjacm(3,1) * xjacm(2,3)
    t3         =  xjacm(2,1) * xjacm(3,2) - xjacm(3,1) * xjacm(2,2)
    deter      =  xjacm(1,1) * t1 + xjacm(1,2) * t2 + xjacm(1,3) * t3

    if( abs(deter) > epsilgeo_tol ) then
       denom       =  1.0_rp / deter
       xjaci(1,1)  =  t1 * denom
       xjaci(2,1)  =  t2 * denom
       xjaci(3,1)  =  t3 * denom
       xjaci(2,2)  =  ( xjacm(1,1) * xjacm(3,3) - xjacm(3,1) * xjacm(1,3)) * denom
       xjaci(3,2)  =  (-xjacm(1,1) * xjacm(3,2) + xjacm(1,2) * xjacm(3,1)) * denom
       xjaci(3,3)  =  ( xjacm(1,1) * xjacm(2,2) - xjacm(2,1) * xjacm(1,2)) * denom
       xjaci(1,2)  =  (-xjacm(1,2) * xjacm(3,3) + xjacm(3,2) * xjacm(1,3)) * denom
       xjaci(1,3)  =  ( xjacm(1,2) * xjacm(2,3) - xjacm(2,2) * xjacm(1,3)) * denom
       xjaci(2,3)  =  (-xjacm(1,1) * xjacm(2,3) + xjacm(2,1) * xjacm(1,3)) * denom

       rhsid(1)    =  (coglo(1) - elcod(1,1))/L
       rhsid(2)    =  (coglo(2) - elcod(2,1))/L
       rhsid(3)    =  (coglo(3) - elcod(3,1))/L
       coloc(1)    =  xjaci(1,1) * rhsid(1) + xjaci(1,2) * rhsid(2) + xjaci(1,3) * rhsid(3)
       coloc(2)    =  xjaci(2,1) * rhsid(1) + xjaci(2,2) * rhsid(2) + xjaci(2,3) * rhsid(3)
       coloc(3)    =  xjaci(3,1) * rhsid(1) + xjaci(3,2) * rhsid(2) + xjaci(3,3) * rhsid(3)

       if( ( coloc(1) >= lmini ) .and. ( coloc(1) <= lmaxi ) ) then
          if( ( coloc(2) >= lmini ) .and. ( coloc(2) <= lmaxi ) ) then
             if( ( coloc(3) >= lmini ) .and. ( coloc(3) <= lmaxi ) ) then
                ezzzt = 1.0_rp - coloc(1) - coloc(2) - coloc(3)
                if( ( ezzzt >= lmini ) .and. ( ezzzt <= lmaxi ) ) then
                   shapf(  1) =  1.0_rp-coloc(1)-coloc(2)-coloc(3)
                   shapf(  2) =  coloc(1)
                   shapf(  3) =  coloc(2)
                   shapf(  4) =  coloc(3)
                   deriv(1,1) = -1.0_rp
                   deriv(2,1) = -1.0_rp
                   deriv(3,1) = -1.0_rp
                   deriv(1,2) =  1.0_rp
                   deriv(2,3) =  1.0_rp
                   deriv(3,4) =  1.0_rp                   
                   ifoun      =  1
                end if
             end if
          end if
       end if

    end if

  end subroutine elmgeo_inside_TET04

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   Natural coordinates in BAR02 lement
  !> @details Check if a point x is inside a BAR02 element
  !>
  !>          x1      x2
  !>          o-------o
  !>                       (x-x1)
  !>          s = -1 + 2 * ------
  !>                       (x2-x1)
  !>
  !-----------------------------------------------------------------------

  subroutine elmgeo_inside_BAR02(&
       pnode,lmaxi,elcod,coglo,coloc,shapf,deriv,ifoun)

    integer(ip), intent(in)  :: pnode
    real(rp),    intent(in)  :: lmaxi
    real(rp),    intent(in)  :: elcod(1,pnode)
    real(rp),    intent(in)  :: coglo(1)
    integer(ip), intent(out) :: ifoun
    real(rp),    intent(out) :: coloc(*)
    real(rp),    intent(out) :: shapf(pnode)
    real(rp),    intent(out) :: deriv(1,pnode)

    coloc(1) = -1.0_rp + 2.0_rp * (coglo(1)-elcod(1,1))/(elcod(1,2)-elcod(1,1))

    if( coloc(1) >= -lmaxi .and. coloc(1) <= lmaxi ) then
       ifoun      =  1
       shapf(1)   =  0.5_rp*(1.0_rp-coloc(1))
       shapf(2)   =  0.5_rp*(1.0_rp+coloc(1))
       deriv(1,1) = -0.5_rp
       deriv(1,2) =  0.5_rp
    else
       ifoun      =  0
    end if

  end subroutine elmgeo_inside_BAR02

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   Natural coordinates
  !> @details Check if point with global coordinates (x,y)=COGLO is inside
  !>          a triangle P1 or a quadrilateral Q1. The Q1 element is
  !>          divided into two P1 elements. Returns the local coordinates
  !>          (s,t)=COLOC
  !>
  !>          For P1 triangles we have:
  !>          x = (1-s-t)*x1 + s*x2 + t*x3
  !>          y = (1-s-t)*y1 + s*y2 + t*y3
  !>
  !>          This linear problem is solved for (s,t):
  !>               (x3-x1)(y-y1) -(y3-y1)(x-x1)
  !>          s =  ----------------------------
  !>               (x3-x1)(y2-y1)-(y3-y1)(x2-x1)
  !>
  !>               (x-x1)(y2-y1) -(y-y1)(x2-x1)
  !>          t =  ----------------------------
  !>               (x3-x1)(y2-y1)-(y3-y1)(x2-x1)
  !
  !-----------------------------------------------------------------------

  subroutine elmgeo_inside_TRI03_QUA04(&
       pnode,lmini,lmaxi,elcod,coglo,coloc,&
       shapf,deriv,ifoun)
    integer(ip), intent(in)  :: pnode
    real(rp),    intent(in)  :: lmini
    real(rp),    intent(in)  :: lmaxi
    real(rp),    intent(in)  :: elcod(2,pnode)
    real(rp),    intent(in)  :: coglo(3)
    integer(ip), intent(out) :: ifoun
    real(rp),    intent(out) :: coloc(*)
    real(rp),    intent(out) :: shapf(pnode)
    real(rp),    intent(out) :: deriv(2,pnode)
    real(rp)                 :: deter,x2x1,y2y1,x3x1,y3y1,xx1,yy1

    ifoun = 0
    !
    ! P1 and Q1: Check if point is in first triangle: nodes 1-2-3
    !
    x2x1     = elcod(1,2) - elcod(1,1)
    y2y1     = elcod(2,2) - elcod(2,1)
    x3x1     = elcod(1,3) - elcod(1,1)
    y3y1     = elcod(2,3) - elcod(2,1)
    xx1      = coglo(1)   - elcod(1,1)
    yy1      = coglo(2)   - elcod(2,1)
    deter    = 1.0_rp / (x3x1*y2y1-y3y1*x2x1)
    coloc(1) = deter  * (x3x1*yy1 -y3y1*xx1)
    if( coloc(1) >= lmini .and. coloc(1) <= lmaxi ) then
       coloc(2) = deter * (y2y1*xx1-x2x1*yy1)
       if( coloc(2) >= lmini .and. coloc(1) + coloc(2) <= lmaxi ) then
          ifoun = 1
       end if
    end if

    if( pnode == 4 ) then
       !
       ! Q1: Check if point is in second triangle: nodes 1-3-4
       !
       if( ifoun == 0 ) then
          x2x1     = elcod(1,3)-elcod(1,1)
          y2y1     = elcod(2,3)-elcod(2,1)
          x3x1     = elcod(1,4)-elcod(1,1)
          y3y1     = elcod(2,4)-elcod(2,1)
          xx1      = coglo(1)  -elcod(1,1)
          yy1      = coglo(2)  -elcod(2,1)
          deter    = 1.0_rp / (x3x1*y2y1-y3y1*x2x1)
          coloc(1) = deter  * (x3x1*yy1 -y3y1*xx1)
          if( coloc(1) >= lmini .and. coloc(1) <= lmaxi ) then
             coloc(2) = deter * (y2y1*xx1-x2x1*yy1)
             if( coloc(2) >= lmini .and. coloc(1) + coloc(2) <= lmaxi ) then
                ifoun = 1
             end if
          end if
       end if

       if( ifoun == 1 ) then
          call elmgeo_newrap(&
               coglo,coloc,2_ip,4_ip,elcod,shapf,deriv)
       end if

    else if( pnode == 3 .and. ifoun == 1 ) then
       !
       ! P1: Compute shape function and derivatives
       !
       shapf(1)   =  1.0_rp - coloc(1) - coloc(2)
       shapf(2)   =  coloc(1)
       shapf(3)   =  coloc(2)
       deriv(1,1) = -1.0_rp
       deriv(1,2) =  1.0_rp
       deriv(1,3) =  0.0_rp
       deriv(2,1) = -1.0_rp
       deriv(2,2) =  0.0_rp
       deriv(2,3) =  1.0_rp

    end if

  end subroutine elmgeo_inside_TRI03_QUA04

  !-----------------------------------------------------------------------
  !
  !> @author  Edgar Olivares
  !> @brief   Natural coordinates
  !> @details Find root of the norm:
  !>          g(s,t,l) = g(x) = || sum_k Nk(s,t) - ( l*r + x0 ) || = 0
  !>          Let define,
  !>          f(s,t,l) = f(x) = sum_k Nk(s,t) - ( l*r + x0 )
  !>          Applying
  !>          Newton - Rapshon --> x^{i+1} = x^i - g(x)/g'(x)
  !>          g'(x) = J(x)*f(x)/g(x)
  !>          => x^{i+1} = xî - g(x)^2 / J(x)f(x)
  !>
  !-----------------------------------------------------------------------

  subroutine elmgeo_newrap_norm(&
       r,xx1,coloc,ndime,pnode,elcod,shapf,deriv)

    real(rp),     intent(in)            :: r(3)
    real(rp),     intent(in)            :: xx1(3)
    real(rp),     intent(out)           :: coloc(*)
    integer(ip),  intent(in)            :: ndime
    integer(ip),  intent(in)            :: pnode
    real(rp),     intent(in)            :: elcod(ndime,pnode)
    real(rp),     intent(out)           :: shapf(pnode)
    real(rp),     intent(out)           :: deriv(ndime,pnode)
    integer(ip)                         :: inode,iiter,jdime,maxit,idime
    real(rp)                            :: xjacm(ndime,ndime)
    real(rp)                            :: xjaci(ndime,ndime)
    real(rp)                            :: deltx(3),delts(3),detja,rhuge,toler
    real(rp)                            :: norm,normm

    !
    ! First guess: Center of the triangle (s,t,u) = (1/3,1/3,1/3)
    coloc(1:ndime) = 0.33_rp
    toler = 1.0e-08_rp
    rhuge = huge(1.0_rp)
    iiter = 0
    maxit = 200
    normm = 1.0_rp
    delts = 1.0_rp

    nr_iterations: do while( maxval(abs(delts(1:ndime))) > 0.03_rp .and. iiter < maxit )
       !
       ! Iteration counter
       !
       iiter = iiter + 1
       !
       ! Compute f^i = sum_k Nk(s^i) xk - l·r - xx1 = DELTX
       !
       call elmgeo_shapf_deriv_heslo(2_ip,pnode,coloc(1:2),shapf,deriv)
       deltx = 0.0_rp
       norm = 0.0_rp
       do inode = 1,pnode
          deltx(1:ndime) = deltx(1:ndime) + shapf(inode) * elcod(1:ndime,inode)
       end do
       deltx(1:ndime) = deltx(1:ndime) - coloc(3)*r(1:ndime) - xx1(1:ndime)
       norm = elmgeo_norm2(ndime,abs(deltx))

       xjacm = 0.0_rp
       do inode = 1,pnode
          do jdime = 1,2
             do idime = 1,ndime
                xjacm(idime,jdime) = xjacm(idime,jdime) + deriv(jdime,inode) * elcod(idime,inode)
             end do
          end do
       end do
       do idime = 1,ndime
          xjacm(idime,3) = - r(idime)
       end do

       call maths_invert_matrix(ndime,xjacm,detja,xjaci)
       !call invmtx(xjacm,xjaci,detja,ndime)
       if( abs(detja) == 0.0_rp ) then
          iiter = maxit + 1
          exit nr_iterations
          xjaci = 0.0_rp
          do idime = 1,ndime
             xjaci(idime,idime) = 0.1_rp
          end do
       end if

       delts = 0.0_rp
       do jdime = 1,ndime
          do idime = 1,ndime
             delts(idime) = delts(idime) - xjaci(idime,jdime) /( deltx(jdime)+0.0001_rp )
          end do
       end do
       delts(1:ndime) = delts(1:ndime)*norm**2
       coloc(1:ndime) = coloc(1:ndime) + delts(1:ndime)

       if ( coloc(1) > rhuge .or. coloc(2) > rhuge .or. norm > 100.0_rp) then
          iiter = maxit + 1
          exit nr_iterations
       end if

    end do nr_iterations

    if( iiter <= maxit ) then
       call elmgeo_shapf_deriv_heslo(ndime,pnode,coloc,shapf,deriv)
    end if

  end subroutine elmgeo_newrap_norm

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   Natural coordinates
  !> @details Calculate the inverse transformation (x,y,z)-->(s,t,r)
  !>
  !>          \verbatim
  !>
  !>          Solve the following system:
  !>
  !>          sum_k Nk(s,t) xk - x = 0 => f1(s,t) = 0
  !>          sum_k Nk(s,t) yk - y = 0 => f2(s,t) = 0
  !>
  !>          Apply Newton-Raphson:
  !>
  !>          +-               -+   +-           -+   +-             -+ +-  -+
  !>          | f1(s^i+1,t^i+1) |   | f1(s^i,t^i) |   | df1/ds df1/dt | | ds |
  !>          |                 | = |             | + |               | |    |
  !>          | f2(s^i+1,t^i+1) |   | f2(s^i,t^i) |   | df2/ds df2/dt | | dt |
  !>          +-               -+   +-           -+   +-             -+ +-  -+
  !>                                     DELTX              XJACM       DELTS
  !>
  !>          => DELTS = - XJACM^-1 . DELTX, with
  !>
  !>          df1/ds      = sum_k dNk(s,t)/ds xk
  !>          df1/dt      = sum_k dNk(s,t)/dt xk
  !>          df2/ds      = sum_k dNk(s,t)/ds yk
  !>          df2/dt      = sum_k dNk(s,t)/dt yk
  !>
  !>          s^i+1       = s^i + DELTS(1)
  !>          t^i+1       = t^i + DELTS(2)
  !>
  !>          f1(s^i,t^i) = sum_k Nk(s,t) xk - x = DELTX(1)
  !>          f2(s^i,t^i) = sum_k Nk(s,t) yk - y = DELTX(2)
  !>
  !>          \endverbatim
  !
  !------------------------------------------------------------------------

  subroutine elmgeo_newrap(&
       coglo_in,coloc,ndime,pnode,elcod_in,shapf,deriv,ifoun,pelty,USER_INITIAL_GUESS)

    integer(ip),  intent(in)            :: ndime
    real(rp),     intent(in)            :: coglo_in(ndime)!(*)
    real(rp),     intent(out)           :: coloc(*)
    integer(ip),  intent(in)            :: pnode
    real(rp),     intent(in)            :: elcod_in(ndime,pnode)
    real(rp),     intent(out)           :: shapf(pnode)
    real(rp),     intent(out)           :: deriv(ndime,pnode)
    integer(ip),  intent(out), optional :: ifoun
    integer(ip),  intent(in),  optional :: pelty
    logical(lg),  intent(in),  optional :: USER_INITIAL_GUESS
    integer(ip)                         :: inode,iiter,jdime,maxit,idime
    integer(ip)                         :: kdime,ielty
    real(rp)                            :: xjacm(ndime,ndime)
    real(rp)                            :: xjaci(ndime,ndime)
    real(rp)                            :: xjaci_old(ndime,ndime)
    real(rp)                            :: deltx(3),delts(3),detja,rhuge,toler

    integer(ip)                         :: lunit,ielem
    real(rp),   allocatable             :: points(:,:)
    real(rp),   allocatable             :: iso_element(:,:)
    
    real(rp)                            :: L
    real(rp)                            :: t(ndime)
    real(rp)                            :: coglo(ndime)
    real(rp)                            :: elcod(ndime,pnode)

    logical(lg), parameter              :: debug_mode          = .false.
    logical(lg), parameter              :: do_translateRescale = .false.
    !
    ! Variables to translate&rescale element and point
    !   t to center at origin, L to reescale first edge to unit length
    !   added  this if true and coglo_in & elcod_in for not being invasive with the code
    !   the cost of the copy of coglo_in & elcod_in is zero compared to NewtonRaphson cost
    !
    if(do_translateRescale) then ! TRANSLATE & RESCALE
      call runend('should be checked, makes unitt_second_try fail')
      t(:)  = elcod(:,1)
      L     = sqrt( sum( (elcod(:,2)-elcod(:,1))**2 ) )
      coglo = (coglo_in-t)/L
      do inode = 1,pnode
        elcod(:,inode) = (elcod_in(:,inode)-t)/L
      end do
    else 
      coglo = coglo_in
      elcod = elcod_in
    end if 
    !
    ! Initial condition s^0 = (0,0,0) or given by user
    !
    if( .not. present(USER_INITIAL_GUESS) ) then
       coloc(1:ndime) = 0.0_rp
       if( present(pelty) ) then
          coloc(1:ndime) = element_type(pelty) % cog(1:ndime)
       end if
    end if
    !
    ! Iterate for ds^i+1 = - (df/ds)^-1 . f(s^i)
    !
    toler = 1.0e-08_rp
    rhuge = huge(1.0_rp)
    iiter = 0
    maxit = 200
    delts = 1.0_rp
    ielty = elmgeo_element_type(ndime,pnode)
    !
    ! Debug mode: save points
    !
    if( debug_mode ) then
       lunit = 90
       allocate(points(ndime,0:maxit))
       allocate(iso_element(ndime,pnode))
       points(1:ndime,iiter) = coloc(1:ndime)
       ielty = elmgeo_element_type(ndime,pnode)
       call elmgeo_iso_parametric_element_coordinates(ndime,ielty,iso_element)
    end if

    nr_iterations: do while( maxval(abs(delts(1:ndime))) > toler .and. iiter < maxit )
       !
       ! Iteration counter
       !
       iiter = iiter + 1
       !write(90,*) iiter,maxval(abs(delts(1:ndime)))
       !
       ! Compute f^i = sum_k Nk(s^i) xk - x = DELTX
       !
       call elmgeo_shapf_deriv_heslo(ndime,pnode,coloc,shapf,deriv)
       deltx(1:ndime) = 0.0_rp
       do inode = 1,pnode
          deltx(1:ndime) = deltx(1:ndime) + shapf(inode) * elcod(1:ndime,inode)
       end do
       deltx(1:ndime) = deltx(1:ndime) - coglo(1:ndime)
       !
       !             +-             -+   +-                                   -+
       !             | df1/ds df1/dt |   | (sum_k DNk/ds xk) (sum_k DNk/dt xk) |
       ! Compute J = |               | = |                                     |
       !             | df2/ds df2/dt |   | (sum_k DNk/ds yk) (sum_k DNk/dt yk) |
       !             +-             -+   +-                                   -+
       !
       xjacm = 0.0_rp
       do inode = 1,pnode
          do jdime = 1,ndime
             do idime = 1,ndime
                xjacm(idime,jdime) = xjacm(idime,jdime) + deriv(jdime,inode) * elcod(idime,inode)
             end do
          end do
       end do
       !
       ! Compute J^{-1}
       !
       if( iiter  > 1 ) xjaci_old = xjaci
       call maths_invert_matrix(ndime,xjacm,detja,xjaci)       
       if( iiter == 1 ) xjaci_old = xjaci
       !
       ! Jacobian treatment.
       ! for PYR05, |J| goes to zero near apex 
       !
       ! Non-singular:     ds^i+1 = -J^-1.f(s^i)
       ! Almost singular:  PYR05: do a very last iteration (very very small)
       ! Singular:         PYR05: assume we one the apex
       !                   OTHER: exit NR
       !
       if( abs(detja) <= epsilgeo_jac .and. ielty == PYR05 ) then

          delts = 0.0_rp
          do jdime = 1,ndime
             delts(1:ndime) = delts(1:ndime) - xjaci(1:ndime,jdime) * deltx(jdime)
          end do
          coloc(1:ndime) = coloc(1:ndime) + delts(1:ndime)
          exit nr_iterations

       else if( abs(detja) == 0.0_rp ) then

          if( ielty == PYR05 ) then
             call elmgeo_iso_parametric_element_coordinates(ndime,ielty,coloc,5_ip)
          else
             iiter = maxit + 1
          end if
          
          exit nr_iterations
          
       end if
       !
       ! Compute ds = -J^{-1}.dx
       !
       delts = 0.0_rp
       do jdime = 1,ndime
          delts(1:ndime) = delts(1:ndime) - xjaci(1:ndime,jdime) * deltx(jdime)
       end do
       coloc(1:ndime) = coloc(1:ndime) + delts(1:ndime)
       
       if( debug_mode ) points(1:ndime,iiter) = coloc(1:ndime)

       if( maxval(abs(coloc(1:ndime))) > rhuge ) then
          iiter = maxit + 1
          exit nr_iterations
       end if

    end do nr_iterations
    !
    ! Debug mode
    !
    if( debug_mode ) then

       ielem = 1
       write(lunit,1)&
            'NEWTON_RAPHSON_'//element_type(ielty) % name,ndime,&
            adjustl(trim(element_type(ielty) % nametopo)),&
            element_type(ielty) % number_nodes
       write(lunit,2) 'coordinates'
       do inode = 1,pnode
          write(lunit,3) inode,iso_element(1:ndime,inode)
       end do
       do inode = 0,iiter
          write(lunit,3) inode+pnode+1,points(1:ndime,inode)
       end do

       write(lunit,2) 'end coordinates'
       write(lunit,2) 'elements'
       write(lunit,4) ielem,(kdime,kdime=1,pnode)
       write(lunit,2) 'end elements'

       ielty = BAR02
       write(lunit,1)&
            'NEWTON_RAPHSON_'//element_type(ielty) % name,ndime,&
            adjustl(trim(element_type(ielty) % nametopo)),&
            element_type(ielty) % number_nodes
       write(lunit,2) 'elements'
       inode = pnode
       do ielem = 0,iiter-1
          inode = inode + 1
          write(lunit,4) ielem+2,inode,inode+1
       end do
       write(lunit,2) 'end elements'

       deallocate(points,iso_element)
    end if
    
    if( iiter >= maxit ) then
       !
       ! NR has not converged
       !
       if( present(ifoun) ) ifoun = 0
       coloc(1) = 2.0_rp
       
       write(3000+kfl_paral,1)&
            'NEWTON_RAPHSON_'//element_type(ielty) % name,ndime,&
            adjustl(trim(element_type(ielty) % nametopo)),&
            element_type(ielty) % number_nodes
       write(3000+kfl_paral,2) 'coordinates'
       do inode = 1, pnode
          write(3000+kfl_paral,3) inode,elcod(1:ndime,inode)
       end do
       write(3000+kfl_paral,3) pnode+1,coglo(1:ndime)
       write(3000+kfl_paral,2) 'end coordinates'
       write(3000+kfl_paral,2) 'elements'
       !do ielem = 1,pnode
       ielem = 1
       write(3000+kfl_paral,'(10(1x,i3))') ielem,(inode,inode=1,pnode)
       !end do
       write(3000+kfl_paral,2) 'end elements'
       write(3000+kfl_paral,1)&
            'NEWTON_RAPHSON_PARTICLE',ndime,'Point',1_ip
       write(3000+kfl_paral,2) 'elements'
       write(3000+kfl_paral,'(10(1x,i3))') ielem,pnode+1
       write(3000+kfl_paral,2) 'end elements'
       flush (3000_ip+kfl_paral)
       call runend('ELMGEO_NEWRAP: NEWTON-RAPHSON HAS NOT CONVERGED')
       
    else
       !
       ! Final shape function and derivative
       !       
       if( present(ifoun) ) ifoun = 1       
       call elmgeo_shapf_deriv_heslo(ndime,pnode,coloc,shapf,deriv)
       
    end if

1   format('MESH ',a,' dimension ',i1,' Elemtype ',a,' Nnode ',i2)
2   format(a)
3   format(i9, 3(1x,e16.8e3))
4   format(i9,50(1x,i9))

  end subroutine elmgeo_newrap

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   Natural coordinates
  !> @details Evaluates shape functions, derivatives and Hessian
  !
  !-----------------------------------------------------------------------

  pure subroutine elmgeo_shapf_deriv_heslo_1(&
       ndime,nnode,posgp,shapf,deriv,heslo,ERROR,INTERPOLATION_FUNCTION)
    
    integer(ip), intent(in)            :: ndime
    integer(ip), intent(in)            :: nnode
    real(rp),    intent(in)            :: posgp(max(1_ip,ndime))
    real(rp),    intent(out)           :: shapf(nnode)
    real(rp),    intent(out), optional :: deriv(max(1_ip,ndime),nnode)
    real(rp),    intent(out), optional :: heslo(max(1_ip,3_ip*ndime-3_ip),nnode)
    integer(ip), intent(in),  optional :: INTERPOLATION_FUNCTION
    integer(ip), intent(out), optional :: ERROR
    integer(ip)                        :: pinte,ierro

    pinte = optional_argument(LAGRANGE_INTERPOLATION,INTERPOLATION_FUNCTION)
    ierro = 0

    call set_isoparametric(ndime,nnode,pinte,&
         posgp,shapf,deriv,heslo,ierro)

    if( present(ERROR) ) ERROR = ierro
    
  end subroutine elmgeo_shapf_deriv_heslo_1

  pure subroutine elmgeo_shapf_deriv_heslo_ngaus(&
       ndime,nnode,ngaus,posgp,shapf,deriv,heslo,ERROR,INTERPOLATION_FUNCTION)
    
    integer(ip), intent(in)            :: ndime
    integer(ip), intent(in)            :: nnode
    integer(ip), intent(in)            :: ngaus
    real(rp),    intent(in)            :: posgp(max(1_ip,ndime),ngaus)
    real(rp),    intent(out)           :: shapf(nnode,ngaus)
    real(rp),    intent(out), optional :: deriv(max(1_ip,ndime),nnode,ngaus)
    real(rp),    intent(out), optional :: heslo(max(1_ip,3_ip*ndime-3_ip),nnode,ngaus)
    integer(ip), intent(in),  optional :: INTERPOLATION_FUNCTION
    integer(ip), intent(out), optional :: ERROR
    integer(ip)                        :: pinte,ierro,igaus

    pinte = optional_argument(LAGRANGE_INTERPOLATION,INTERPOLATION_FUNCTION)
    ierro = 0

    do igaus = 1,ngaus
       call set_isoparametric(ndime,nnode,pinte,&
            posgp(:,igaus),shapf(:,igaus),deriv(:,:,igaus),heslo(:,:,igaus),ierro)
    end do
    
    if( present(ERROR) ) ERROR = ierro
    
  end subroutine elmgeo_shapf_deriv_heslo_ngaus

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   Natural coordinates
  !> @details Evaluates shape functions, derivatives and Hessian
  !>          for th bubble
  !
  !-----------------------------------------------------------------------

  subroutine elmgeo_shapf_deriv_heslo_bubble_1(&
       ndime,nnode,posgp,shapf,deriv,heslo)
    
    integer(ip), intent(in)            :: ndime
    integer(ip), intent(in)            :: nnode
    real(rp),    intent(in)            :: posgp(max(1_ip,ndime))
    real(rp),    intent(out)           :: shapf(nnode)
    real(rp),    intent(out), optional :: deriv(ndime,nnode)
    real(rp),    intent(out), optional :: heslo(max(1_ip,3_ip*ndime-3_ip),nnode)

    if( ndime == 0 ) then
       return
    else if( ndime == 1 ) then
       return
    else if( ndime == 2 ) then
       if( nnode == 3 ) then
          call elmgeo_shape_bubble_triangle(posgp(1),posgp(2),shapf,deriv,heslo)
       else if( nnode == 4 ) then
          call elmgeo_shape_bubble_quadrilateral(posgp(1),posgp(2),shapf,deriv,heslo)
       end if
    else if( ndime == 3 ) then
       return
       call runend('BUBBLE NOT CODED')
    end if

  end subroutine elmgeo_shapf_deriv_heslo_bubble_1

  subroutine elmgeo_shapf_deriv_heslo_bubble_ngaus(&
       ndime,nnode,ngaus,posgp,shapf,deriv,heslo)
    
    integer(ip), intent(in)            :: ndime
    integer(ip), intent(in)            :: nnode
    integer(ip), intent(in)            :: ngaus
    real(rp),    intent(in)            :: posgp(max(1_ip,ndime),ngaus)
    real(rp),    intent(out)           :: shapf(nnode,ngaus)
    real(rp),    intent(out), optional :: deriv(ndime,nnode,ngaus)
    real(rp),    intent(out), optional :: heslo(max(1_ip,3_ip*ndime-3_ip),nnode,ngaus)
    integer(ip)                        :: igaus
    
    if( ndime == 0 ) then
       return
    else if( ndime == 1 ) then
       return
    else if( ndime == 2 ) then
       if( nnode == 3 ) then
          do igaus = 1,ngaus
             call elmgeo_shape_bubble_triangle(posgp(1,igaus),posgp(2,igaus),shapf(:,igaus),deriv(:,:,igaus),heslo(:,:,igaus))
          end do
       else if( nnode == 4 ) then
          do igaus = 1,ngaus
             call elmgeo_shape_bubble_quadrilateral(posgp(1,igaus),posgp(2,igaus),shapf(:,igaus),deriv(:,:,igaus),heslo(:,:,igaus))
          end do
       end if
    else if( ndime == 3 ) then
       return
       call runend('BUBBLE NOT CODED')
    end if

  end subroutine elmgeo_shapf_deriv_heslo_bubble_ngaus

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux and beatriz Eguzkitza
  !> @date    01/10/2015
  !> @brief   Bubble function for triangle
  !> @details Bubble function for triangle
  !>
  !----------------------------------------------------------------------

  subroutine elmgeo_shape_bubble_triangle(s,t,shapf,deriv,heslo)
    
    real(rp),    intent(in)            :: s,t
    real(rp),    intent(out)           :: shapf(1)
    real(rp),    intent(out), optional :: deriv(2,1)
    real(rp),    intent(out), optional :: heslo(3,1)
    integer(ip)                        :: ii,jj
    real(rp)                           :: z

    if( present(heslo) ) then
       do ii = 1,1
          do jj = 1,3
             heslo(jj,ii) = 0.0_rp
          end do
       end do
    end if
    if( present(deriv) ) then
       do ii = 1,1
          do jj = 1,2
             deriv(jj,ii) = 0.0_rp
          end do
       end do
    end if
    z        =  1.0_rp-s-t
    shapf(1) =  27.0_rp*s*t*z
    if( present(deriv) ) then
       deriv(1,1) = 27.0_rp*(t*z-s*t)
       deriv(2,1) = 27.0_rp*(s*z-s*t)
    end if
    if( present(heslo) ) then
       !deriv(1,1) = 27.0_rp*(t*z-s*t)
       !deriv(2,1) = 27.0_rp*(s*z-s*t)
    end if

  end subroutine elmgeo_shape_bubble_triangle

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux and beatriz Eguzkitza
  !> @date    01/10/2015
  !> @brief   Bubble function for TET04
  !> @details Bubble function for TET04
  !>
  !----------------------------------------------------------------------

  subroutine elmgeo_shape_bubble_tet04(s,t,shapf,deriv,heslo)
    implicit none
    real(rp),    intent(in)            :: s,t
    real(rp),    intent(out)           :: shapf(1)
    real(rp),    intent(out), optional :: deriv(2,1)
    real(rp),    intent(out), optional :: heslo(3,1)
    integer(ip)                        :: ii,jj
    real(rp)                           :: z

    if( present(heslo) ) then
       do ii = 1,1
          do jj = 1,3
             heslo(jj,ii) = 0.0_rp
          end do
       end do
    end if
    if( present(deriv) ) then
       do ii = 1,1
          do jj = 1,2
             deriv(jj,ii) = 0.0_rp
          end do
       end do
    end if
    z        =  1.0_rp-s-t
    shapf(1) =  27.0_rp*s*t*z
    if( present(deriv) ) then
       deriv(1,1) = 27.0_rp*(t*z-s*t)
       deriv(2,1) = 27.0_rp*(s*z-s*t)
    end if
    if( present(heslo) ) then
       !deriv(1,1) = 27.0_rp*(t*z-s*t)
       !deriv(2,1) = 27.0_rp*(s*z-s*t)
    end if

  end subroutine elmgeo_shape_bubble_tet04

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux and beatriz Eguzkitza
  !> @date    01/10/2015
  !> @brief   Bubble function for triangle
  !> @details Bubble function for triangle
  !>
  !----------------------------------------------------------------------

  subroutine elmgeo_shape_bubble_quadrilateral(s,t,shapf,deriv,heslo)
    implicit none
    real(rp),    intent(in)            :: s,t
    real(rp),    intent(out)           :: shapf(1)
    real(rp),    intent(out), optional :: deriv(2,1)
    real(rp),    intent(out), optional :: heslo(3,1)
    integer(ip)                        :: ii,jj
!    real(rp)                           :: z

    if( present(heslo) ) then
       do ii = 1,1
          do jj = 1,3
             heslo(jj,ii) = 0.0_rp
          end do
       end do
    end if
    if( present(deriv) ) then
       do ii = 1,1
          do jj = 1,2
             deriv(jj,ii) = 0.0_rp
          end do
       end do
    end if

    shapf(1) = (1.0_rp-s)*(1.0_rp+s)*(1.0_rp-t)*(1.0_rp+t)

    if( present(deriv) ) then
       deriv(1,1) = (1.0_rp-t)*(1.0_rp+t)*2.0_rp*s
       deriv(2,1) = (1.0_rp-s)*(1.0_rp+s)*2.0_rp*t
    end if
    if( present(heslo) ) then
       !deriv(1,1) = 27.0_rp*(t*z-s*t)
       !deriv(2,1) = 27.0_rp*(s*z-s*t)
    end if

  end subroutine elmgeo_shape_bubble_quadrilateral

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Tells where a point is inside an element
  !> @details Tells where a point is inside an element
  !>          PTOPO        2D            3D
  !>          -----------------------------------
  !>          0        Quadrilateral  Hexahedra
  !>          1        Triangle       Tetrahedra
  !>          2               -       Pentahedra
  !>          3               -       Pyramid
  !>
  !----------------------------------------------------------------------

  subroutine elmgeo_where_is(&
       ndime,pnode,ptopo,pelty,elcod,coglo,winou,wwher,llist)
    implicit none
    integer(ip),   intent(in)            :: ndime
    integer(ip),   intent(in)            :: pnode
    integer(ip),   intent(in)            :: ptopo
    integer(ip),   intent(in)            :: pelty
    real(rp),      intent(in)            :: elcod(ndime,pnode)
    real(rp),      intent(in)            :: coglo(ndime)
    character(20), intent(out)           :: winou
    character(20), intent(out)           :: wwher
    integer(ip),   intent(out), optional :: llist(*)
    integer(ip)                          :: ifoun
    integer(ip)                          :: inode
    integer(ip)                          :: jnode
    real(rp)                             :: coloc(3)
    real(rp)                             :: deriv(3*64)
    real(rp)                             :: shapf(64)
    real(rp)                             :: lmaxi
    real(rp)                             :: lmini
    real(rp)                             :: zeror

    zeror =  1.0e-12_rp
    lmini = -zeror
    lmaxi =  1.0_rp+zeror
    call elmgeo_natural_coordinates(    &
         ndime,pelty,pnode,elcod,shapf, &
         deriv,coglo,coloc,ifoun,zeror)

    !call elsest_chkelm(&
    !     ndime,ptopo,pnode,elcod,shapf,deriv,&
    !     coglo,coloc,ifoun,lmini,lmaxi)

    if( ifoun == 0 ) then
       winou = 'OUTSIDE'
       wwher = 'NO_WHERE'
       return
    else
       winou = 'INSIDE'
    end if

    if( ifoun == 0 ) return

    if( ptopo == 0 ) then
       if( ndime == 2 ) then
          if( abs(abs(coloc(1)) - 1.0_rp) < zeror .and. abs(abs(coloc(2)) - 1.0_rp) < zeror ) then
             wwher = 'VERTEX'
          else if( abs(abs(coloc(1)) - 1.0_rp) < zeror .and. abs(coloc(2)) < zeror ) then
             wwher = 'MID_EDGE'
          else if( abs(abs(coloc(2)) - 1.0_rp) < zeror .and. abs(coloc(1)) < zeror ) then
             wwher = 'MID_EDGE'
          else if( abs(coloc(1)) < zeror .and. abs(coloc(2)) < zeror ) then
             wwher = 'COG'
          else
             wwher = 'ANYWHERE'
          end if
       else
          call runend('ELEMNT_WHERIS: NOT CODED')
       end if
    else
       call runend('ELEMNT_WHERIS: NOT CODED')
    end if

    if(present(llist)) then
       if( trim(wwher) == 'VERTEX' ) then
          inode = 0
          jnode = 0
          do while( inode < pnode )
             inode = inode + 1
             if(        abs(coglo(    1) - elcod(    1,inode)) <= zeror &
                  .and. abs(coglo(    2) - elcod(    2,inode)) <= zeror &
                  .and. abs(coglo(ndime) - elcod(ndime,inode)) <= zeror ) then
                jnode = inode
                inode = pnode
             end if
          end do
          llist(1) = jnode
       else if( trim(wwher) == 'MID_EDGE' ) then
          if(      abs(coloc(1)) < zeror        .and. abs(coloc(2)+1.0_rp) < zeror ) then
             llist(1) = 1
             llist(2) = 2
          else if( abs(coloc(1)-1.0_rp) < zeror .and. abs(coloc(2)) < zeror ) then
             llist(1) = 2
             llist(2) = 3
          else if( abs(coloc(1)) < zeror        .and. abs(coloc(2)-1.0_rp) < zeror ) then
             llist(1) = 3
             llist(2) = 4
          else
             llist(1) = 4
             llist(2) = 1
          end if
       else
          llist(1) = 0
       end if
    end if

  end subroutine elmgeo_where_is

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   Element length
  !> @details \verbatim
  !>          HLENG ... Element length with:
  !>                    HLENG(1)     = Max length
  !>                    HLENG(NDIME) = Min length
  !>          \endverbatim
  !>
  !-----------------------------------------------------------------------

  subroutine elmgeo_element_node_length(&
       ndime,pnode,elcod,hleng)
    implicit none
    integer(ip), intent(in)           :: ndime
    integer(ip), intent(in)           :: pnode
    real(rp),    intent(in)           :: elcod(ndime,pnode)
    real(rp),    intent(out)          :: hleng(ndime)
    integer(ip)                       :: inode,jnode,idime
    real(rp)                          :: dista

    hleng(1)     = -huge(1.0_rp)
    hleng(ndime) =  huge(1.0_rp)
    do inode = 1,pnode
       do jnode = inode+1,pnode
          dista = 0.0_rp
          do idime = 1,ndime
             dista = dista &
                  & + (elcod(idime,inode)-elcod(idime,jnode)) &
                  & * (elcod(idime,inode)-elcod(idime,jnode))
          end do
          dista = sqrt(dista)
          hleng(1)     = max(hleng(1),    dista)
          hleng(ndime) = min(hleng(ndime),dista)
       end do
    end do
    if( ndime == 3 ) hleng(2) = 0.0_rp

  end subroutine elmgeo_element_node_length

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   Element characteristic length
  !> @details \verbatim
  !>          HLENG ... Element length with:
  !>                    HLENG(1)     = Max length
  !>                    HLENG(NDIME) = Min length
  !>          \endverbatim
  !>
  !-----------------------------------------------------------------------

  subroutine elmgeo_element_characteristic_length_scalar(&
       ndime,pnode,deriv,elcod,hleng,hnatu_opt,tragl_opt)
    implicit none
    integer(ip), intent(in)           :: ndime,pnode
    real(rp),    intent(in)           :: elcod(ndime,pnode)
    real(rp),    intent(in)           :: deriv(ndime,pnode)
    real(rp),    intent(out)          :: hleng(ndime)
    real(rp),    intent(in), optional :: hnatu_opt
    real(rp),    intent(out),optional :: tragl_opt(ndime,ndime)
    integer(ip)                       :: inode
    real(rp)                          :: enor0,h_tem,gpdet,denom
    real(rp)                          :: xjacm(3,3),t1,t2,t3
    real(rp)                          :: tragl(ndime,ndime),hnatu

    if( present(hnatu_opt) ) then
       hnatu = hnatu_opt
    else
       hnatu = 1.0_rp
    end if

    if( ndime == 1 ) then

       xjacm(1,1) = 0.0_rp
       do inode = 1,pnode
          xjacm(1,1) = xjacm(1,1) + elcod(1,inode) * deriv(1,inode)
       end do
       tragl(1,1) = 1.0_rp/xjacm(1,1)
       enor0      = tragl(1,1) * tragl(1,1)
       hleng(1)   = hnatu/sqrt(enor0)

    else if( ndime == 2 ) then

       xjacm(1,1) = 0.0_rp
       xjacm(1,2) = 0.0_rp
       xjacm(2,1) = 0.0_rp
       xjacm(2,2) = 0.0_rp
       do inode = 1,pnode
          xjacm(1,1) = xjacm(1,1) + elcod(1,inode) * deriv(1,inode)
          xjacm(1,2) = xjacm(1,2) + elcod(1,inode) * deriv(2,inode)
          xjacm(2,1) = xjacm(2,1) + elcod(2,inode) * deriv(1,inode)
          xjacm(2,2) = xjacm(2,2) + elcod(2,inode) * deriv(2,inode)
       end do

       gpdet      =  xjacm(1,1) * xjacm(2,2) - xjacm(2,1) * xjacm(1,2)
       denom      =  1.0_rp/gpdet
       tragl(1,1) =  xjacm(2,2) * denom
       tragl(2,2) =  xjacm(1,1) * denom
       tragl(2,1) = -xjacm(2,1) * denom
       tragl(1,2) = -xjacm(1,2) * denom

       enor0      =  tragl(1,1) * tragl(1,1) + tragl(1,2) * tragl(1,2)
       hleng(1)   =  hnatu/sqrt(enor0)
       enor0      =  tragl(2,1) * tragl(2,1) + tragl(2,2) * tragl(2,2)
       hleng(2)   =  hnatu/sqrt(enor0)

       if( hleng(2) > hleng(1) ) then
          h_tem    = hleng(2)
          hleng(2) = hleng(1)
          hleng(1) = h_tem
       end if

    else if( ndime == 3 ) then

       xjacm(1,1) = 0.0_rp ! xjacm = elcod * deriv^t
       xjacm(1,2) = 0.0_rp ! tragl = xjacm^-1
       xjacm(1,3) = 0.0_rp
       xjacm(2,1) = 0.0_rp
       xjacm(2,2) = 0.0_rp
       xjacm(2,3) = 0.0_rp
       xjacm(3,1) = 0.0_rp
       xjacm(3,2) = 0.0_rp
       xjacm(3,3) = 0.0_rp
       do inode = 1,pnode
          xjacm(1,1) = xjacm(1,1) + elcod(1,inode) * deriv(1,inode)
          xjacm(1,2) = xjacm(1,2) + elcod(1,inode) * deriv(2,inode)
          xjacm(1,3) = xjacm(1,3) + elcod(1,inode) * deriv(3,inode)
          xjacm(2,1) = xjacm(2,1) + elcod(2,inode) * deriv(1,inode)
          xjacm(2,2) = xjacm(2,2) + elcod(2,inode) * deriv(2,inode)
          xjacm(2,3) = xjacm(2,3) + elcod(2,inode) * deriv(3,inode)
          xjacm(3,1) = xjacm(3,1) + elcod(3,inode) * deriv(1,inode)
          xjacm(3,2) = xjacm(3,2) + elcod(3,inode) * deriv(2,inode)
          xjacm(3,3) = xjacm(3,3) + elcod(3,inode) * deriv(3,inode)
       end do

       t1         =  xjacm(2,2) * xjacm(3,3) - xjacm(3,2) * xjacm(2,3)
       t2         = -xjacm(2,1) * xjacm(3,3) + xjacm(3,1) * xjacm(2,3)
       t3         =  xjacm(2,1) * xjacm(3,2) - xjacm(3,1) * xjacm(2,2)
       gpdet      =  xjacm(1,1) * t1 + xjacm(1,2) * t2 + xjacm(1,3) * t3
       denom      =  1.0_rp / gpdet
       tragl(1,1) =  t1 * denom
       tragl(2,1) =  t2 * denom
       tragl(3,1) =  t3 * denom
       tragl(2,2) =  ( xjacm(1,1) * xjacm(3,3) - xjacm(3,1) * xjacm(1,3)) * denom
       tragl(3,2) =  (-xjacm(1,1) * xjacm(3,2) + xjacm(1,2) * xjacm(3,1)) * denom
       tragl(3,3) =  ( xjacm(1,1) * xjacm(2,2) - xjacm(2,1) * xjacm(1,2)) * denom
       tragl(1,2) =  (-xjacm(1,2) * xjacm(3,3) + xjacm(3,2) * xjacm(1,3)) * denom
       tragl(1,3) =  ( xjacm(1,2) * xjacm(2,3) - xjacm(2,2) * xjacm(1,3)) * denom
       tragl(2,3) =  (-xjacm(1,1) * xjacm(2,3) + xjacm(2,1) * xjacm(1,3)) * denom
       !
       ! Element length HLENG
       !
       enor0    =   tragl(1,1) * tragl(1,1) &
            &     + tragl(1,2) * tragl(1,2) &
            &     + tragl(1,3) * tragl(1,3)
       hleng(1) = hnatu/sqrt(enor0)
       enor0    =   tragl(2,1) * tragl(2,1) &
            &     + tragl(2,2) * tragl(2,2) &
            &     + tragl(2,3) * tragl(2,3)
       hleng(2) = hnatu/sqrt(enor0)
       enor0    =   tragl(3,1) * tragl(3,1) &
            &     + tragl(3,2) * tragl(3,2) &
            &     + tragl(3,3) * tragl(3,3)
       hleng(3) = hnatu/sqrt(enor0)
       !
       ! Sort hleng: hleng(1)=max; hleng(ndime)=min
       !
       if( hleng(2) > hleng(1) ) then
          h_tem    = hleng(2)
          hleng(2) = hleng(1)
          hleng(1) = h_tem
       end if
       if( hleng(3) > hleng(1) ) then
          h_tem    = hleng(3)
          hleng(3) = hleng(1)
          hleng(1) = h_tem
       end if
       if( hleng(3) > hleng(2) ) then
          h_tem    = hleng(3)
          hleng(3) = hleng(2)
          hleng(2) = h_tem
       end if

    end if

    if( present(tragl_opt) ) then
       tragl_opt = tragl
    end if

  end subroutine elmgeo_element_characteristic_length_scalar

  subroutine elmgeo_element_characteristic_length_vector(&
       ndime,pnode,deriv,elcod,hleng,hnatu_opt,tragl_opt)
    implicit none
    integer(ip), intent(in)           :: ndime
    integer(ip), intent(in)           :: pnode
    real(rp),    intent(in)           :: deriv(ndime,pnode)
    real(rp),    intent(in)           :: elcod(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(out)          :: hleng(VECTOR_SIZE,ndime)
    real(rp),    intent(in), optional :: hnatu_opt
    real(rp),    intent(out),optional :: tragl_opt(VECTOR_SIZE,ndime,ndime)
    integer(ip)                       :: inode,ivect
    real(rp)                          :: enor0(VECTOR_SIZE)
    real(rp)                          :: tragl(VECTOR_SIZE,ndime,ndime)
    real(rp)                          :: xjacm(VECTOR_SIZE,3,3)
    real(rp)                          :: gpdet(VECTOR_SIZE)
    real(rp)                          :: denom(VECTOR_SIZE)
    real(rp)                          :: t1(VECTOR_SIZE)
    real(rp)                          :: t2(VECTOR_SIZE)
    real(rp)                          :: t3(VECTOR_SIZE)
    real(rp)                          :: h_tem,hnatu

    if( present(hnatu_opt) ) then
       hnatu = hnatu_opt
    else
       hnatu = 1.0_rp
    end if

    if( ndime == 1 ) then

       xjacm = 0.0_rp
       do inode = 1,pnode
          xjacm(1:VECTOR_SIZE,1,1) = xjacm(1:VECTOR_SIZE,1,1) + elcod(1:VECTOR_SIZE,1,inode) * deriv(1,inode)
       end do
       !tragl(1:VECTOR_SIZE,1,1) = 1.0_rp/xjacm(1:VECTOR_SIZE,1,1)
       tragl(1:VECTOR_SIZE,1,1) = 1.0_rp / (sign(1.0_rp,xjacm(1:VECTOR_SIZE,1,1))*max(abs(xjacm(1:VECTOR_SIZE,1,1)),epsilgeo_div))
       enor0(1:VECTOR_SIZE)     = tragl(1:VECTOR_SIZE,1,1) * tragl(1:VECTOR_SIZE,1,1)
       hleng(1:VECTOR_SIZE,1)   = hnatu/sqrt(enor0)

    else if( ndime == 2 ) then

       xjacm = 0.0_rp
       do inode = 1,pnode
          xjacm(1:VECTOR_SIZE,1,1) = xjacm(1:VECTOR_SIZE,1,1) + elcod(1:VECTOR_SIZE,1,inode) * deriv(1,inode)
          xjacm(1:VECTOR_SIZE,1,2) = xjacm(1:VECTOR_SIZE,1,2) + elcod(1:VECTOR_SIZE,1,inode) * deriv(2,inode)
          xjacm(1:VECTOR_SIZE,2,1) = xjacm(1:VECTOR_SIZE,2,1) + elcod(1:VECTOR_SIZE,2,inode) * deriv(1,inode)
          xjacm(1:VECTOR_SIZE,2,2) = xjacm(1:VECTOR_SIZE,2,2) + elcod(1:VECTOR_SIZE,2,inode) * deriv(2,inode)
       end do

       gpdet(1:VECTOR_SIZE)     =  xjacm(1:VECTOR_SIZE,1,1) * xjacm(1:VECTOR_SIZE,2,2) - xjacm(1:VECTOR_SIZE,2,1) * xjacm(1:VECTOR_SIZE,1,2)
       denom(1:VECTOR_SIZE)     =  1.0_rp / (sign(1.0_rp,gpdet(1:VECTOR_SIZE))*max(abs(gpdet(1:VECTOR_SIZE)),epsilgeo_div))
       !denom(1:VECTOR_SIZE)     =  1.0_rp / gpdet(1:VECTOR_SIZE)

       tragl(1:VECTOR_SIZE,1,1) =  xjacm(1:VECTOR_SIZE,2,2) * denom(1:VECTOR_SIZE)
       tragl(1:VECTOR_SIZE,2,2) =  xjacm(1:VECTOR_SIZE,1,1) * denom(1:VECTOR_SIZE)
       tragl(1:VECTOR_SIZE,2,1) = -xjacm(1:VECTOR_SIZE,2,1) * denom(1:VECTOR_SIZE)
       tragl(1:VECTOR_SIZE,1,2) = -xjacm(1:VECTOR_SIZE,1,2) * denom (1:VECTOR_SIZE)

       enor0(1:VECTOR_SIZE)     =  tragl(1:VECTOR_SIZE,1,1) * tragl(1:VECTOR_SIZE,1,1) + tragl(1:VECTOR_SIZE,1,2) * tragl(1:VECTOR_SIZE,1,2)
       denom(1:VECTOR_SIZE)    =   1.0_rp / max(sqrt(enor0(1:VECTOR_SIZE)),epsilgeo_div)
       hleng(1:VECTOR_SIZE,1)  =   hnatu * denom(1:VECTOR_SIZE)

       enor0(1:VECTOR_SIZE)     =  tragl(1:VECTOR_SIZE,2,1) * tragl(1:VECTOR_SIZE,2,1) + tragl(1:VECTOR_SIZE,2,2) * tragl(1:VECTOR_SIZE,2,2)
       denom(1:VECTOR_SIZE)    =   1.0_rp / max(sqrt(enor0(1:VECTOR_SIZE)),epsilgeo_div)
       hleng(1:VECTOR_SIZE,2)  =   hnatu * denom(1:VECTOR_SIZE)

       do ivect = 1,VECTOR_SIZE
          if( hleng(ivect,2) > hleng(ivect,1) ) then
             h_tem          = hleng(ivect,2)
             hleng(ivect,2) = hleng(ivect,1)
             hleng(ivect,1) = h_tem
          end if
       end do

    else if( ndime == 3 ) then
       !
       ! tragl = xjacm^-1
       ! xjacm = elcod * deriv^t
       !
       xjacm = 0.0_rp
       do inode = 1,pnode
          xjacm(1:VECTOR_SIZE,1,1) = xjacm(1:VECTOR_SIZE,1,1) + elcod(1:VECTOR_SIZE,1,inode) * deriv(1,inode)
          xjacm(1:VECTOR_SIZE,1,2) = xjacm(1:VECTOR_SIZE,1,2) + elcod(1:VECTOR_SIZE,1,inode) * deriv(2,inode)
          xjacm(1:VECTOR_SIZE,1,3) = xjacm(1:VECTOR_SIZE,1,3) + elcod(1:VECTOR_SIZE,1,inode) * deriv(3,inode)
          xjacm(1:VECTOR_SIZE,2,1) = xjacm(1:VECTOR_SIZE,2,1) + elcod(1:VECTOR_SIZE,2,inode) * deriv(1,inode)
          xjacm(1:VECTOR_SIZE,2,2) = xjacm(1:VECTOR_SIZE,2,2) + elcod(1:VECTOR_SIZE,2,inode) * deriv(2,inode)
          xjacm(1:VECTOR_SIZE,2,3) = xjacm(1:VECTOR_SIZE,2,3) + elcod(1:VECTOR_SIZE,2,inode) * deriv(3,inode)
          xjacm(1:VECTOR_SIZE,3,1) = xjacm(1:VECTOR_SIZE,3,1) + elcod(1:VECTOR_SIZE,3,inode) * deriv(1,inode)
          xjacm(1:VECTOR_SIZE,3,2) = xjacm(1:VECTOR_SIZE,3,2) + elcod(1:VECTOR_SIZE,3,inode) * deriv(2,inode)
          xjacm(1:VECTOR_SIZE,3,3) = xjacm(1:VECTOR_SIZE,3,3) + elcod(1:VECTOR_SIZE,3,inode) * deriv(3,inode)
       end do

       t1(1:VECTOR_SIZE)        =  xjacm(1:VECTOR_SIZE,2,2) * xjacm(1:VECTOR_SIZE,3,3) - xjacm(1:VECTOR_SIZE,3,2) * xjacm(1:VECTOR_SIZE,2,3)
       t2(1:VECTOR_SIZE)        = -xjacm(1:VECTOR_SIZE,2,1) * xjacm(1:VECTOR_SIZE,3,3) + xjacm(1:VECTOR_SIZE,3,1) * xjacm(1:VECTOR_SIZE,2,3)
       t3(1:VECTOR_SIZE)        =  xjacm(1:VECTOR_SIZE,2,1) * xjacm(1:VECTOR_SIZE,3,2) - xjacm(1:VECTOR_SIZE,3,1) * xjacm(1:VECTOR_SIZE,2,2)

       gpdet(1:VECTOR_SIZE)     =  xjacm(1:VECTOR_SIZE,1,1) * t1(1:VECTOR_SIZE) + xjacm(1:VECTOR_SIZE,1,2) * t2 (1:VECTOR_SIZE)+ xjacm(1:VECTOR_SIZE,1,3) * t3(1:VECTOR_SIZE)
       denom(1:VECTOR_SIZE)     =  1.0_rp / (sign(1.0_rp,gpdet(1:VECTOR_SIZE))*max(abs(gpdet(1:VECTOR_SIZE)),epsilgeo_div))
       !denom(1:VECTOR_SIZE)     =  1.0_rp / gpdet(1:VECTOR_SIZE)

       tragl(1:VECTOR_SIZE,1,1) =  t1(1:VECTOR_SIZE) * denom(1:VECTOR_SIZE)
       tragl(1:VECTOR_SIZE,2,1) =  t2(1:VECTOR_SIZE) * denom(1:VECTOR_SIZE)
       tragl(1:VECTOR_SIZE,3,1) =  t3(1:VECTOR_SIZE) * denom(1:VECTOR_SIZE)
       tragl(1:VECTOR_SIZE,2,2) =  ( xjacm(1:VECTOR_SIZE,1,1) * xjacm(1:VECTOR_SIZE,3,3) - xjacm(1:VECTOR_SIZE,3,1) * xjacm(1:VECTOR_SIZE,1,3)) * denom(1:VECTOR_SIZE)
       tragl(1:VECTOR_SIZE,3,2) =  (-xjacm(1:VECTOR_SIZE,1,1) * xjacm(1:VECTOR_SIZE,3,2) + xjacm(1:VECTOR_SIZE,1,2) * xjacm(1:VECTOR_SIZE,3,1)) * denom(1:VECTOR_SIZE)
       tragl(1:VECTOR_SIZE,3,3) =  ( xjacm(1:VECTOR_SIZE,1,1) * xjacm(1:VECTOR_SIZE,2,2) - xjacm(1:VECTOR_SIZE,2,1) * xjacm(1:VECTOR_SIZE,1,2)) * denom(1:VECTOR_SIZE)
       tragl(1:VECTOR_SIZE,1,2) =  (-xjacm(1:VECTOR_SIZE,1,2) * xjacm(1:VECTOR_SIZE,3,3) + xjacm(1:VECTOR_SIZE,3,2) * xjacm(1:VECTOR_SIZE,1,3)) * denom(1:VECTOR_SIZE)
       tragl(1:VECTOR_SIZE,1,3) =  ( xjacm(1:VECTOR_SIZE,1,2) * xjacm(1:VECTOR_SIZE,2,3) - xjacm(1:VECTOR_SIZE,2,2) * xjacm(1:VECTOR_SIZE,1,3)) * denom(1:VECTOR_SIZE)
       tragl(1:VECTOR_SIZE,2,3) =  (-xjacm(1:VECTOR_SIZE,1,1) * xjacm(1:VECTOR_SIZE,2,3) + xjacm(1:VECTOR_SIZE,2,1) * xjacm(1:VECTOR_SIZE,1,3)) * denom(1:VECTOR_SIZE)
       !
       ! Element length HLENG
       !
       enor0(1:VECTOR_SIZE)    =   tragl(1:VECTOR_SIZE,1,1) * tragl(1:VECTOR_SIZE,1,1) &
            &                    + tragl(1:VECTOR_SIZE,1,2) * tragl(1:VECTOR_SIZE,1,2) &
            &                    + tragl(1:VECTOR_SIZE,1,3) * tragl(1:VECTOR_SIZE,1,3)

       denom(1:VECTOR_SIZE)    =   1.0_rp / max(sqrt(enor0(1:VECTOR_SIZE)),epsilgeo_div)
       hleng(1:VECTOR_SIZE,1)  =   hnatu * denom(1:VECTOR_SIZE)

       enor0(1:VECTOR_SIZE)    =   tragl(1:VECTOR_SIZE,2,1) * tragl(1:VECTOR_SIZE,2,1) &
            &                    + tragl(1:VECTOR_SIZE,2,2) * tragl(1:VECTOR_SIZE,2,2) &
            &                    + tragl(1:VECTOR_SIZE,2,3) * tragl(1:VECTOR_SIZE,2,3)
       denom(1:VECTOR_SIZE)    =   1.0_rp / max(sqrt(enor0(1:VECTOR_SIZE)),epsilgeo_div)
       hleng(1:VECTOR_SIZE,2)  =   hnatu * denom(1:VECTOR_SIZE)

       enor0(1:VECTOR_SIZE)    =   tragl(1:VECTOR_SIZE,3,1) * tragl(1:VECTOR_SIZE,3,1) &
            &                    + tragl(1:VECTOR_SIZE,3,2) * tragl(1:VECTOR_SIZE,3,2) &
            &                    + tragl(1:VECTOR_SIZE,3,3) * tragl(1:VECTOR_SIZE,3,3)
       denom(1:VECTOR_SIZE)    =   1.0_rp / max(sqrt(enor0(1:VECTOR_SIZE)),epsilgeo_div)
       hleng(1:VECTOR_SIZE,3)  =   hnatu * denom(1:VECTOR_SIZE)
       !
       ! Sort hleng: hleng(1)=max; hleng(ndime)=min
       !
       do ivect = 1,VECTOR_SIZE
          if( hleng(ivect,2) > hleng(ivect,1) ) then
             h_tem          = hleng(ivect,2)
             hleng(ivect,2) = hleng(ivect,1)
             hleng(ivect,1) = h_tem
          end if
          if( hleng(ivect,3) > hleng(ivect,1) ) then
             h_tem          = hleng(ivect,3)
             hleng(ivect,3) = hleng(ivect,1)
             hleng(ivect,1) = h_tem
          end if
          if( hleng(ivect,3) > hleng(ivect,2) ) then
             h_tem          = hleng(ivect,3)
             hleng(ivect,3) = hleng(ivect,2)
             hleng(ivect,2) = h_tem
          end if
       end do

    end if

    if( present(tragl_opt) ) then
       tragl_opt = tragl
    end if

  end subroutine elmgeo_element_characteristic_length_vector

  !-----------------------------------------------------------------------
  !> @author  Guillaume Houzeaux
  !> @date    27/01/2007
  !> @brief   Computes the Jacobian, Jacobian determinant and
  !>          Cartesian derivatives of shape function of an element
  !> @details The jacobian is
  !>
  !-----------------------------------------------------------------------

  pure subroutine elmgeo_jacobian_boundary(ndime,pnode,elcod,deriv,gpdet,xjaci,xjacm)

    integer(ip),           intent(in)  :: ndime               !< Dimension
    integer(ip),           intent(in)  :: pnode               !< Number of nodes
    real(rp),              intent(in)  :: elcod(ndime,pnode)  !< Element coordinates
    real(rp),              intent(in)  :: deriv(ndime,pnode)  !< Shape function derivatives
    real(rp),              intent(out) :: gpdet               !< Determinant
    real(rp),    optional, intent(out) :: xjaci(ndime,ndime)  !< Inverse Jacobian
    real(rp),    optional, intent(out) :: xjacm(ndime,ndime)  !< Jacobian
    real(rp)                           :: xjac2(3,3)
    
    if( present(xjacm) ) then
       call bouder(pnode,ndime,ndime-1_ip,deriv,elcod,xjacm,gpdet)
    else
       call bouder(pnode,ndime,ndime-1_ip,deriv,elcod,xjac2,gpdet)       
    end if
    
  end subroutine elmgeo_jacobian_boundary
  
  !-----------------------------------------------------------------------
  !> @author  Guillaume Houzeaux
  !> @date    27/01/2007
  !> @brief   Computes the Jacobian, Jacobian determinant and
  !>          Cartesian derivatives of shape function of an element
  !> @details The jacobian is
  !>                             _           _
  !>                            | dx/ds dx/dt |                t
  !>          Jacobian: XJACM = |             | = ELCOD * DERIV
  !>                            |_dy/ds dy/dt_|
  !>          with
  !>                  _        _             _                    _
  !>                 | x1 x2 x3 |           | dN1/ds dN2/ds dN3/ds |
  !>         ELCOD = |          |,  DERIV = |                      |
  !>                 |_y1 y2 y3_|           |_dN1/dt dN2/dt dN3/dt_|
  !>
  !>           => Jacobian determinant: GPDET = det(XJACM)
  !>
  !-----------------------------------------------------------------------

  pure subroutine elmgeo_jacobian_matrix(ndime,pnode,elcod,deriv,gpdet,xjaci,xjacm)

    integer(ip),           intent(in)  :: ndime               !< Dimension
    integer(ip),           intent(in)  :: pnode               !< Number of nodes
    real(rp),              intent(in)  :: elcod(ndime,pnode)  !< Element coordinates
    real(rp),              intent(in)  :: deriv(ndime,pnode)  !< Shape function derivatives
    real(rp),              intent(out) :: gpdet               !< Determinant
    real(rp),    optional, intent(out) :: xjaci(ndime,ndime)  !< Inverse Jacobian
    real(rp),    optional, intent(out) :: xjacm(ndime,ndime)  !< Jacobian
    integer(ip)                        :: k
    real(rp)                           :: t1,t2,t3,denom
    real(rp)                           :: xjac2(3,3)          !< Jacobian

    if( ndime == 2 .and. pnode == 3 ) then
       !
       ! 2D P1 element
       !
       !gpdet = (-elcod(1,1)+elcod(1,2))*(-elcod(2,1)+elcod(2,3)) &
       !     & -(-elcod(2,1)+elcod(2,2))*(-elcod(1,1)+elcod(1,3))
       !if( gpdet == 0.0_rp ) return
       !denom = 1.0_rp/gpdet

       xjac2(1,1) = 0.0_rp
       xjac2(1,2) = 0.0_rp
       xjac2(2,1) = 0.0_rp
       xjac2(2,2) = 0.0_rp
       do k = 1,3
          xjac2(1,1) = xjac2(1,1) + elcod(1,k) * deriv(1,k)
          xjac2(1,2) = xjac2(1,2) + elcod(1,k) * deriv(2,k)
          xjac2(2,1) = xjac2(2,1) + elcod(2,k) * deriv(1,k)
          xjac2(2,2) = xjac2(2,2) + elcod(2,k) * deriv(2,k)
       end do

       gpdet = xjac2(1,1) * xjac2(2,2) - xjac2(2,1) * xjac2(1,2)
       if( gpdet == 0.0_rp ) return
       denom = 1.0_rp/gpdet

       if( present(xjaci) ) then
          xjaci(1,1) =  xjac2(2,2) * denom
          xjaci(2,2) =  xjac2(1,1) * denom
          xjaci(2,1) = -xjac2(2,1) * denom
          xjaci(1,2) = -xjac2(1,2) * denom
       end if
       
    else if( ndime == 3 .and. pnode == 4 ) then
       !
       ! 3D P1 element
       !
       xjac2(1,1) =  elcod(1,2) - elcod(1,1)
       xjac2(1,2) =  elcod(1,3) - elcod(1,1)
       xjac2(1,3) =  elcod(1,4) - elcod(1,1)
       xjac2(2,1) =  elcod(2,2) - elcod(2,1)
       xjac2(2,2) =  elcod(2,3) - elcod(2,1)
       xjac2(2,3) =  elcod(2,4) - elcod(2,1)
       xjac2(3,1) =  elcod(3,2) - elcod(3,1)
       xjac2(3,2) =  elcod(3,3) - elcod(3,1)
       xjac2(3,3) =  elcod(3,4) - elcod(3,1)
       t1         =  xjac2(2,2) * xjac2(3,3) - xjac2(3,2) * xjac2(2,3)
       t2         = -xjac2(2,1) * xjac2(3,3) + xjac2(3,1) * xjac2(2,3)
       t3         =  xjac2(2,1) * xjac2(3,2) - xjac2(3,1) * xjac2(2,2)
       gpdet      =  xjac2(1,1) * t1 + xjac2(1,2) * t2 + xjac2(1,3) * t3
       if( gpdet == 0.0_rp ) return
       denom      =  1.0_rp/gpdet

       if( present(xjaci) ) then
          xjaci(1,1) =  t1 * denom
          xjaci(2,1) =  t2 * denom
          xjaci(3,1) =  t3 * denom
          xjaci(2,2) = ( xjac2(1,1) * xjac2(3,3) - xjac2(3,1) * xjac2(1,3)) * denom
          xjaci(3,2) = (-xjac2(1,1) * xjac2(3,2) + xjac2(1,2) * xjac2(3,1)) * denom
          xjaci(3,3) = ( xjac2(1,1) * xjac2(2,2) - xjac2(2,1) * xjac2(1,2)) * denom
          xjaci(1,2) = (-xjac2(1,2) * xjac2(3,3) + xjac2(3,2) * xjac2(1,3)) * denom
          xjaci(1,3) = ( xjac2(1,2) * xjac2(2,3) - xjac2(2,2) * xjac2(1,3)) * denom
          xjaci(2,3) = (-xjac2(1,1) * xjac2(2,3) + xjac2(2,1) * xjac2(1,3)) * denom
       end if
       
    else if ( ndime == 1 ) then
       !
       ! 1D
       !
       xjac2(1,1) = 0.0_rp
       do k = 1,pnode
          xjac2(1,1) = xjac2(1,1) + elcod(1,k) * deriv(1,k)
       end do
       gpdet = xjac2(1,1)
       if( gpdet == 0.0_rp ) return
       if( present(xjaci) ) then
          xjaci(1,1) = 1.0_rp/xjac2(1,1)
       end if

    else if ( ndime == 2 ) then
       !
       ! 2D
       !
       xjac2(1,1) = 0.0_rp
       xjac2(1,2) = 0.0_rp
       xjac2(2,1) = 0.0_rp
       xjac2(2,2) = 0.0_rp
       do k = 1,pnode
          xjac2(1,1) = xjac2(1,1) + elcod(1,k) * deriv(1,k)
          xjac2(1,2) = xjac2(1,2) + elcod(1,k) * deriv(2,k)
          xjac2(2,1) = xjac2(2,1) + elcod(2,k) * deriv(1,k)
          xjac2(2,2) = xjac2(2,2) + elcod(2,k) * deriv(2,k)
       end do

       gpdet = xjac2(1,1) * xjac2(2,2) - xjac2(2,1) * xjac2(1,2)
       if( gpdet == 0.0_rp ) return
       denom = 1.0_rp/gpdet

       if( present(xjaci) ) then       
          xjaci(1,1) =  xjac2(2,2) * denom
          xjaci(2,2) =  xjac2(1,1) * denom
          xjaci(2,1) = -xjac2(2,1) * denom
          xjaci(1,2) = -xjac2(1,2) * denom
       end if
       
    else if ( ndime == 3 ) then
       !
       ! 3D
       !
       xjac2(1,1) = 0.0_rp ! xjac2 = elcod * deriv^t
       xjac2(1,2) = 0.0_rp ! xjaci = xjac2^-1
       xjac2(1,3) = 0.0_rp ! gpcar = xjaci^t * deriv
       xjac2(2,1) = 0.0_rp
       xjac2(2,2) = 0.0_rp
       xjac2(2,3) = 0.0_rp
       xjac2(3,1) = 0.0_rp
       xjac2(3,2) = 0.0_rp
       xjac2(3,3) = 0.0_rp
       do k = 1,pnode
          xjac2(1,1) = xjac2(1,1) + elcod(1,k) * deriv(1,k)
          xjac2(1,2) = xjac2(1,2) + elcod(1,k) * deriv(2,k)
          xjac2(1,3) = xjac2(1,3) + elcod(1,k) * deriv(3,k)
          xjac2(2,1) = xjac2(2,1) + elcod(2,k) * deriv(1,k)
          xjac2(2,2) = xjac2(2,2) + elcod(2,k) * deriv(2,k)
          xjac2(2,3) = xjac2(2,3) + elcod(2,k) * deriv(3,k)
          xjac2(3,1) = xjac2(3,1) + elcod(3,k) * deriv(1,k)
          xjac2(3,2) = xjac2(3,2) + elcod(3,k) * deriv(2,k)
          xjac2(3,3) = xjac2(3,3) + elcod(3,k) * deriv(3,k)
       end do

       t1    =  xjac2(2,2) * xjac2(3,3) - xjac2(3,2) * xjac2(2,3)
       t2    = -xjac2(2,1) * xjac2(3,3) + xjac2(3,1) * xjac2(2,3)
       t3    =  xjac2(2,1) * xjac2(3,2) - xjac2(3,1) * xjac2(2,2)
       gpdet =  xjac2(1,1) * t1 + xjac2(1,2) * t2 + xjac2(1,3) * t3
       if(gpdet == 0.0_rp ) return

       if( present(xjaci) ) then       
          denom = 1.0_rp / gpdet
          xjaci(1,1) = t1*denom
          xjaci(2,1) = t2*denom
          xjaci(3,1) = t3*denom
          xjaci(2,2) = ( xjac2(1,1) * xjac2(3,3) - xjac2(3,1) * xjac2(1,3)) * denom
          xjaci(3,2) = (-xjac2(1,1) * xjac2(3,2) + xjac2(1,2) * xjac2(3,1)) * denom
          xjaci(3,3) = ( xjac2(1,1) * xjac2(2,2) - xjac2(2,1) * xjac2(1,2)) * denom
          xjaci(1,2) = (-xjac2(1,2) * xjac2(3,3) + xjac2(3,2) * xjac2(1,3)) * denom
          xjaci(1,3) = ( xjac2(1,2) * xjac2(2,3) - xjac2(2,2) * xjac2(1,3)) * denom
          xjaci(2,3) = (-xjac2(1,1) * xjac2(2,3) + xjac2(2,1) * xjac2(1,3)) * denom
       end if
       
    end if

    if( present(xjacm) ) xjacm(1:ndime,1:ndime) = xjac2(1:ndime,1:ndime)
    
  end subroutine elmgeo_jacobian_matrix

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux and J.C. Cajas
  !> @date    09/11/2015
  !> @brief   Local coordinates on a boundary
  !> @details Determine the local coordinates on a 3D boundary
  !>
  !-----------------------------------------------------------------------

  recursive subroutine elmgeo_natural_coordinates_on_boundaries(&
       ndime,pblty,pnodb,bocod,shapb,derib,coglo,coloc,ifoun,toler_opt,NEAREST_POINT)

    integer(ip),           intent(in)    :: ndime
    integer(ip),           intent(in)    :: pblty
    integer(ip),           intent(in)    :: pnodb
    real(rp),              intent(in)    :: bocod(ndime,pnodb)
    real(rp),              intent(out)   :: shapb(pnodb)
    real(rp),              intent(out)   :: derib(ndime-1,pnodb)
    real(rp),              intent(in)    :: coglo(ndime)
    real(rp),              intent(out)   :: coloc(ndime-1)
    integer(ip),           intent(out)   :: ifoun
    real(rp),    optional, intent(in)    :: toler_opt
    real(rp),    optional, intent(inout) :: NEAREST_POINT(:)
    integer(ip)                          :: ndimb
    real(rp)                             :: toler
    real(rp)                             :: lmini,lmaxi,ezzzt
    integer(ip)                          :: ii,iedge,p1,p2
    real(rp)                             :: pp(3),d,dimin,pp_min(3)
    
    if( present(toler_opt) ) then
       toler = toler_opt
    else
       toler = 0.01_rp
    end if
    lmini = - toler
    lmaxi = 1.0_rp + toler
    ndimb = ndime - 1
    ifoun = 0

    if( pblty == BAR02 ) then
       !
       ! BAR02 element
       !
       call elmgeo_natural_coordinates_on_BAR02(bocod,coglo,coloc)
       call elmgeo_shapf_deriv_heslo(1_ip,2_ip,coloc,shapb,derib)
       !call elmgeo_shape1(coloc(1),2_ip,shapb,derib)
       if( coloc(1) >= -lmaxi .and. coloc(1) <= lmaxi ) ifoun = 1

    else if( pblty == TRI03 ) then
       !
       ! TRI03 element
       !
       call elmgeo_natural_coordinates_on_TRI03(bocod,coglo,coloc)
       call elmgeo_shapf_deriv_heslo(2_ip,3_ip,coloc(1:2),shapb,derib)

       !call elmgeo_shape2(coloc(1),coloc(2),3_ip,shapb,derib)

       if(    coloc(1) >= lmini .and. coloc(1) <= lmaxi ) then
          if( coloc(2) >= lmini .and. coloc(2) <= lmaxi ) then
             ezzzt = 1.0_rp-coloc(1)-coloc(2)
             if( ezzzt >= lmini .and. ezzzt <= lmaxi ) then
                ifoun = 1
             end if
          end if
       end if

    else if( pblty == QUA04 ) then
       !
       ! QUA04 element
       !
       call elmgeo_natural_coordinates_on_QUA04(bocod,coglo,coloc)
       call elmgeo_shapf_deriv_heslo(2_ip,4_ip,coloc(1:2),shapb,derib)
       !call elmgeo_shape2(coloc(1),coloc(2),4_ip,shapb,derib)

       if(   coloc(1) >= -lmaxi .and. coloc(1) <= lmaxi ) then
          if( coloc(2)>= -lmaxi .and. coloc(2) <= lmaxi ) then
             ifoun = 1
          end if
       end if

    else

       call runend('ELMGEO_NATURAL_COORDINATES_ON_BOUNDARIES: PRORGRAMALO SI TE ATREVES!')

    end if
    !
    ! Look for nearest point
    !
    if( ifoun == 0 .and. present(NEAREST_POINT) ) then
       dimin = huge(1.0_rp)
       if( element_type(pblty) % dimensions == 1 ) then
          do ii = 1,pnodb
             d = dot_product(bocod(1:ndime,ii)-coglo(1:ndime),bocod(1:ndime,ii)-coglo(1:ndime))
             if( d < dimin ) then
                dimin = d
                pp_min(1:ndime) = bocod(1:ndime,ii)
             end if
          end do
       else
          do iedge = 1,element_type(pblty) % number_faces
             p1 = element_type(pblty) % list_faces(1,iedge)
             p2 = element_type(pblty) % list_faces(2,iedge)
             call maths_distance_point_segment(bocod(1:ndime,p1),bocod(1:ndime,p2),coglo,pp,d)
             if( d < dimin ) then
                pp_min(1:ndime) = pp(1:ndime)
                dimin = d
             end if
          end do
       end if
       call elmgeo_natural_coordinates_on_boundaries(&
            ndime,pblty,pnodb,bocod,shapb,derib,pp_min,coloc,ifoun,toler_opt)
       NEAREST_POINT(1:ndime) = pp_min(1:ndime)
       
    else if( present(NEAREST_POINT) ) then
       
       NEAREST_POINT(1:ndime) = coglo(1:ndime)
       
    end if

  end subroutine elmgeo_natural_coordinates_on_boundaries

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   Check if an element is in element bounding box
  !> @details Check if an element is in element bounding box
  !
  !------------------------------------------------------------------------

  pure function elmgeo_inside_element_bounding_box_1(ndime,pnode,elcod,coglo,toler)
    
    integer(ip), intent(in)           :: ndime
    integer(ip), intent(in)           :: pnode
    real(rp),    intent(in)           :: elcod(ndime,pnode)
    real(rp),    intent(in)           :: coglo(ndime)
    real(rp),    intent(in)           :: toler
    real(rp)                          :: xdist(3),xmini(3),xmaxi(3)
    logical(lg)                       :: elmgeo_inside_element_bounding_box_1
    
    if( ndime == 2 ) then

       xmini(1)   = minval(elcod(1,1:pnode))
       xmaxi(1)   = maxval(elcod(1,1:pnode))
       xmini(2)   = minval(elcod(2,1:pnode))
       xmaxi(2)   = maxval(elcod(2,1:pnode))
       xdist(1:2) = toler * (xmaxi(1:2)-xmini(1:2))
       xmini(1:2) = xmini(1:2) - xdist(1:2)
       xmaxi(1:2) = xmaxi(1:2) + xdist(1:2)

       if(      coglo(1) > xmaxi(1)  ) then
          elmgeo_inside_element_bounding_box_1 = .false.
       else if( coglo(2) > xmaxi(2)  ) then
          elmgeo_inside_element_bounding_box_1 = .false.
       else if( coglo(1) < xmini(1)  ) then
          elmgeo_inside_element_bounding_box_1 = .false.
       else if( coglo(2) < xmini(2)  ) then
          elmgeo_inside_element_bounding_box_1 = .false.
       else
          elmgeo_inside_element_bounding_box_1 = .true.
       end if

    else if( pnode == 4 ) then

       xmini(1)   = minval(elcod(1,1:4))
       xmaxi(1)   = maxval(elcod(1,1:4))
       xmini(2)   = minval(elcod(2,1:4))
       xmaxi(2)   = maxval(elcod(2,1:4))
       xmini(3)   = minval(elcod(3,1:4))
       xmaxi(3)   = maxval(elcod(3,1:4))

       xdist(1:3) = toler * (xmaxi(1:3)-xmini(1:3))
       xmini(1:3) = xmini(1:3) - xdist(1:3)
       xmaxi(1:3) = xmaxi(1:3) + xdist(1:3)
       
       if(      coglo(1) > xmaxi(1)  ) then
          elmgeo_inside_element_bounding_box_1 = .false.
       else if( coglo(2) > xmaxi(2)  ) then
          elmgeo_inside_element_bounding_box_1 = .false.
       else if( coglo(3) > xmaxi(3)  ) then
          elmgeo_inside_element_bounding_box_1 = .false.
       else if( coglo(1) < xmini(1)  ) then
          elmgeo_inside_element_bounding_box_1 = .false.
       else if( coglo(2) < xmini(2)  ) then
          elmgeo_inside_element_bounding_box_1 = .false.
       else if( coglo(3) < xmini(3)  ) then
          elmgeo_inside_element_bounding_box_1 = .false.
       else
          elmgeo_inside_element_bounding_box_1 = .true.
       end if

    else

       xmini(1)   = minval(elcod(1,1:pnode))
       xmaxi(1)   = maxval(elcod(1,1:pnode))
       xmini(2)   = minval(elcod(2,1:pnode))
       xmaxi(2)   = maxval(elcod(2,1:pnode))
       xmini(3)   = minval(elcod(3,1:pnode))
       xmaxi(3)   = maxval(elcod(3,1:pnode))

       if( toler > 0.0_rp ) then
          xdist(1:3) = toler * (xmaxi(1:3)-xmini(1:3))
          xmini(1:3) = xmini(1:3) - xdist(1:3)
          xmaxi(1:3) = xmaxi(1:3) + xdist(1:3)
       end if
       
       if(      coglo(1) > xmaxi(1)  ) then
          elmgeo_inside_element_bounding_box_1 = .false.
       else if( coglo(2) > xmaxi(2)  ) then
          elmgeo_inside_element_bounding_box_1 = .false.
       else if( coglo(3) > xmaxi(3)  ) then
          elmgeo_inside_element_bounding_box_1 = .false.
       else if( coglo(1) < xmini(1)  ) then
          elmgeo_inside_element_bounding_box_1 = .false.
       else if( coglo(2) < xmini(2)  ) then
          elmgeo_inside_element_bounding_box_1 = .false.
       else if( coglo(3) < xmini(3)  ) then
          elmgeo_inside_element_bounding_box_1 = .false.
       else
          elmgeo_inside_element_bounding_box_1 = .true.
       end if

    end if
    
  end function elmgeo_inside_element_bounding_box_1

  pure function elmgeo_inside_element_bounding_box_2(ndime,xmini,xmaxi,coglo,pelty)
    
    integer(ip), intent(in)           :: ndime
    real(rp),    intent(in)           :: xmini(ndime)
    real(rp),    intent(in)           :: xmaxi(ndime)
    real(rp),    intent(in)           :: coglo(ndime)
    integer(ip), intent(in), optional :: pelty
    logical(lg)                       :: elmgeo_inside_element_bounding_box_2
!    real(rp)                          :: time1

    !if( present(pelty) ) call cputim(time1)
    
    if( ndime == 2 ) then

       if(      coglo(1) > xmaxi(1)  ) then
          elmgeo_inside_element_bounding_box_2 = .false.
       else if( coglo(2) > xmaxi(2)  ) then
          elmgeo_inside_element_bounding_box_2 = .false.
       else if( coglo(1) < xmini(1)  ) then
          elmgeo_inside_element_bounding_box_2 = .false.
       else if( coglo(2) < xmini(2)  ) then
          elmgeo_inside_element_bounding_box_2 = .false.
       else
          elmgeo_inside_element_bounding_box_2 = .true.
       end if

    else
       
       if(      coglo(1) > xmaxi(1)  ) then
          elmgeo_inside_element_bounding_box_2 = .false.
       else if( coglo(2) > xmaxi(2)  ) then
          elmgeo_inside_element_bounding_box_2 = .false.
       else if( coglo(3) > xmaxi(3)  ) then
          elmgeo_inside_element_bounding_box_2 = .false.
       else if( coglo(1) < xmini(1)  ) then
          elmgeo_inside_element_bounding_box_2 = .false.
       else if( coglo(2) < xmini(2)  ) then
          elmgeo_inside_element_bounding_box_2 = .false.
       else if( coglo(3) < xmini(3)  ) then
          elmgeo_inside_element_bounding_box_2 = .false.
       else
          elmgeo_inside_element_bounding_box_2 = .true.
       end if

    end if
    
  end function elmgeo_inside_element_bounding_box_2

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   Transform coordinates from one tetra to another
  !> @details Transform coordinates in a tetrahedron to the coordinates
  !>          in the tetrahedra coming from the decomposition of
  !>          pyramid, prisms and hexahedra
  !
  !------------------------------------------------------------------------

  subroutine elmgeo_tetrahedra_to_tetrahedra(pelty,itetr,coloc)

    integer(ip), intent(in)    :: pelty         !< Element type
    integer(ip), intent(in)    :: itetr         !< Tetrahedra number
    real(rp),    intent(inout) :: coloc(3)      !< Final coordinate
    real(rp)                   :: elcod(3,4)
    real(rp)                   :: xjacm(3,3)
    real(rp)                   :: coloc_cpy(3)

    if( pelty == PYR05 .and. itetr == 1 ) then
       !
       ! First tetrahedra of PYR05 (1,2,4,5)
       !
       elcod(1,1) = -1.0_rp
       elcod(2,1) = -1.0_rp
       elcod(3,1) = -1.0_rp

       elcod(1,2) =  1.0_rp
       elcod(2,2) = -1.0_rp
       elcod(3,2) = -1.0_rp

       elcod(1,3) = -1.0_rp
       elcod(2,3) =  1.0_rp
       elcod(3,3) = -1.0_rp

       elcod(1,4) =  0.0_rp
       elcod(2,4) =  0.0_rp
       elcod(3,4) =  1.0_rp

    else if( pelty == PYR05 .and. itetr == 2 ) then
       !
       ! Second tetrahedra of PYR05 (2,3,4,5)
       !
       elcod(1,1) =  1.0_rp
       elcod(2,1) = -1.0_rp
       elcod(3,1) = -1.0_rp

       elcod(1,2) =  1.0_rp
       elcod(2,2) =  1.0_rp
       elcod(3,2) = -1.0_rp

       elcod(1,3) = -1.0_rp
       elcod(2,3) =  1.0_rp
       elcod(3,3) = -1.0_rp

       elcod(1,4) =  0.0_rp
       elcod(2,4) =  0.0_rp
       elcod(3,4) =  1.0_rp

    else if( pelty == PEN06 .and. itetr == 1 ) then
       !
       ! Second tetrahedra of PEN06 (1,2,3,5)
       !
       elcod(1,1) =  0.0_rp
       elcod(2,1) =  0.0_rp
       elcod(3,1) =  0.0_rp

       elcod(1,2) =  1.0_rp
       elcod(2,2) =  0.0_rp
       elcod(3,2) =  0.0_rp

       elcod(1,3) =  0.0_rp
       elcod(2,3) =  1.0_rp
       elcod(3,3) =  0.0_rp

       elcod(1,4) =  1.0_rp
       elcod(2,4) =  0.0_rp
       elcod(3,4) =  1.0_rp

    else if( pelty == PEN06 .and. itetr == 2 ) then
       !
       ! Second tetrahedra of PEN06 (4,1,6,5)
       !
       elcod(1,1) =  0.0_rp
       elcod(2,1) =  0.0_rp
       elcod(3,1) =  1.0_rp

       elcod(1,2) =  0.0_rp
       elcod(2,2) =  0.0_rp
       elcod(3,2) =  0.0_rp

       elcod(1,3) =  0.0_rp
       elcod(2,3) =  1.0_rp
       elcod(3,3) =  1.0_rp

       elcod(1,4) =  1.0_rp
       elcod(2,4) =  0.0_rp
       elcod(3,4) =  1.0_rp

    else if( pelty == PEN06 .and. itetr == 3 ) then
       !
       ! Second tetrahedra of PEN06 (3,6,1,5)
       !
       elcod(1,1) =  0.0_rp
       elcod(2,1) =  1.0_rp
       elcod(3,1) =  0.0_rp

       elcod(1,2) =  0.0_rp
       elcod(2,2) =  1.0_rp
       elcod(3,2) =  1.0_rp

       elcod(1,3) =  0.0_rp
       elcod(2,3) =  0.0_rp
       elcod(3,3) =  0.0_rp

       elcod(1,4) =  1.0_rp
       elcod(2,4) =  0.0_rp
       elcod(3,4) =  1.0_rp

    end if

    coloc_cpy    = coloc
    xjacm(1:3,1) = elcod(1:3,2) - elcod(1:3,1)
    xjacm(1:3,2) = elcod(1:3,3) - elcod(1:3,1)
    xjacm(1:3,3) = elcod(1:3,4) - elcod(1:3,1)
    coloc(1)     = elcod(1,1)   + dot_product(xjacm(1,1:3),coloc_cpy(1:3))
    coloc(2)     = elcod(2,1)   + dot_product(xjacm(2,1:3),coloc_cpy(1:3))
    coloc(3)     = elcod(3,1)   + dot_product(xjacm(3,1:3),coloc_cpy(1:3))

  end subroutine elmgeo_tetrahedra_to_tetrahedra

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   Nearest point on element faces
  !> @details Look for the nearest point XX element faces
  !
  !-----------------------------------------------------------------------

  subroutine elmgeo_nearest_point_on_element_faces(&
       ndime,pelty,elcod,coglo,coloc,shapf,deriv,dista,toler_opt)
    integer(ip), intent(in)            :: ndime
    integer(ip), intent(in)            :: pelty
    real(rp),    intent(in)            :: elcod(ndime,*)
    real(rp),    intent(in)            :: coglo(ndime)
    real(rp),    intent(out)           :: coloc(*)             !< Parametrical coordinates
    real(rp),    intent(out)           :: shapf(*)             !< Shape functions of test point
    real(rp),    intent(out)           :: deriv(ndime,*)       !< Shape funcitons derivatives of test point
    real(rp),    intent(out)           :: dista                !< Original distance to element
    real(rp),    intent(in),  optional :: toler_opt            !< Tolerance
    integer(ip)                        :: iface,pnodf,inodf
    integer(ip)                        :: inode,pflty
    integer(ip)                        :: pnode,ifoun,jnode
    integer(ip)                        :: jfoun
    real(rp)                           :: bocod(ndime,16)
    real(rp)                           :: nn(3),xx(3),toler
    real(rp)                           :: ezzzt,lmini,lmaxi
    real(rp)                           :: xdist,xx_min(3)


    dista = huge(1.0_rp)

    if( present(toler_opt) ) then
       toler = toler_opt
    else
       toler = epsilgeo_def
    end if
    lmini = -toler
    lmaxi = 1.0_rp + toler
    pnode = element_type(pelty) % number_nodes
    ifoun = 0
    xdist = huge(1.0_rp)

    loop_faces: do iface = 1,element_type(pelty) % number_faces

       jfoun = 0
       pnodf = element_type(pelty) % node_faces(iface)
       pflty = element_type(pelty) % type_faces(iface)
       do inodf = 1,pnodf
          inode = element_type(pelty) % list_faces(inodf,iface)
          bocod(1:ndime,inodf) = elcod(1:ndime,inode)
       end do
       !
       ! Projection point XX on faces
       !
       if( pflty == BAR02 ) then
          !
          ! BAR02
          !
          call elmgeo_nearest_point_on_BAR02(bocod,coglo,xx,nn)
          call elmgeo_natural_coordinates_on_BAR02(bocod,xx,coloc)
          if( coloc(1) >= -lmaxi .and. coloc(1) <= lmaxi ) then
             jfoun = 1
          end if

       else if( pflty == TRI03 ) then
          !
          ! TRI03
          !
          call elmgeo_nearest_point_on_TRI03(bocod,coglo,xx,nn)
          call elmgeo_natural_coordinates_on_TRI03(bocod,xx,coloc)

          if(    coloc(1) >= lmini .and. coloc(1) <= lmaxi ) then
             if( coloc(2) >= lmini .and. coloc(2) <= lmaxi ) then
                ezzzt = 1.0_rp-coloc(1)-coloc(2)
                if( ezzzt >= lmini .and. ezzzt <= lmaxi ) then
                   jfoun = 1
                end if
             end if
          end if

       else if( pflty == QUA04 ) then
          !
          ! QUA04
          !
          call elmgeo_nearest_point_on_QUA04(bocod,coglo,xx,nn)
          call elmgeo_natural_coordinates_on_QUA04(bocod,xx,coloc)

          if(   coloc(1) >= -lmaxi .and. coloc(1) <= lmaxi ) then
             if( coloc(2)>= -lmaxi .and. coloc(2) <= lmaxi ) then
                jfoun = 1
             end if
          end if

       else

          call runend('ELMGEO_INSIDE_ELEMENT_USING_FACES: UNKNOWN TYPE OF FACE')

       end if

       if( jfoun == 1 ) then
          ifoun = 1
          xdist = sqrt(dot_product(xx(1:ndime)-coglo(1:ndime),xx(1:ndime)-coglo(1:ndime)))
          if( xdist < dista ) then
             dista = xdist
             xx_min(1:ndime) = xx(1:ndime)
          end if
       end if

    end do loop_faces
    !
    ! Find natural coordinate of the projection XX in the element
    !
    if( ifoun == 1 ) then
       call elmgeo_natural_coordinates(    &
            ndime,pelty,pnode,elcod,shapf, &
            deriv,xx_min,coloc,ifoun,toler)
    end if
    !
    ! If we are not the element look for nearest node
    !
    if( ifoun == 0 ) then
       dista = huge(1.0_rp)
       do inode = 1,pnode
          xdist = dot_product( coglo(1:ndime)-elcod(1:ndime,inode),coglo(1:ndime)-elcod(1:ndime,inode) )
          if( xdist <= dista ) then
             jnode = inode
             dista = xdist
          end if
       end do
       xx_min(1:ndime) = elcod(1:ndime,jnode)
       dista = sqrt(dista)
       call elmgeo_natural_coordinates(    &
            ndime,pelty,pnode,elcod,shapf, &
            deriv,xx_min,coloc,ifoun,toler)
    end if

  end subroutine elmgeo_nearest_point_on_element_faces

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   Nearest point on a face
  !> @details Porject point COGLO on a face with coordinates BOCOD
  !
  !-----------------------------------------------------------------------

  pure subroutine elmgeo_projection_on_a_face(&
       ndime,pflty,bocod,coglo,projection)
    integer(ip), intent(in)            :: ndime
    integer(ip), intent(in)            :: pflty
    real(rp),    intent(in)            :: bocod(ndime,*)
    real(rp),    intent(in)            :: coglo(ndime)
    real(rp),    intent(out)           :: projection(ndime)     !< Coordinates of the projection point
    real(rp)                           :: nn(3)

    if( pflty == BAR02 ) then
       !
       ! BAR02
       !
       call elmgeo_nearest_point_on_BAR02(bocod,coglo,projection,nn)

    else if( pflty == TRI03 ) then
       !
       ! TRI03
       !
       call elmgeo_nearest_point_on_TRI03(bocod,coglo,projection,nn)

    else if( pflty == QUA04 ) then
       !
       ! QUA04
       !
       call elmgeo_nearest_point_on_QUA04(bocod,coglo,projection,nn)

    else

       !call runend('ELMGEO_PROJECTION_ON_A_FACE: UNKNOWN TYPE OF FACE')

    end if

  end subroutine elmgeo_projection_on_a_face

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   Projection on a line
  !> @details Look for the nearest point XX to COGLO on a line
  !
  !-----------------------------------------------------------------------

  pure subroutine elmgeo_nearest_point_on_BAR02(bocod,coglo,xx,nn)
    real(rp),    intent(in)  :: bocod(2,2)
    real(rp),    intent(in)  :: coglo(2)
    real(rp),    intent(out) :: xx(2)
    real(rp),    intent(out) :: nn(2)
    real(rp)                 :: t

    nn(1)   = - ( bocod(2,2) - bocod(2,1) )
    nn(2)   =   ( bocod(1,2) - bocod(1,1) )
    nn(1:2) =   nn(1:2) / sqrt( nn(1)*nn(1) + nn(2)*nn(2) )
    t       =   dot_product(nn(1:2),bocod(1:2,1)-coglo(1:2))
    xx(1:2) =   coglo(1:2) + t * nn(1:2)

  end subroutine elmgeo_nearest_point_on_BAR02

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   Projection on a linear surface
  !> @details Look for the nearest point XX to COGLO on a linear
  !>          surface
  !
  !-----------------------------------------------------------------------

  pure function elmgeo_side_of_TRI03(bocod,xx1,xx2)
    real(rp),    intent(in) :: bocod(3,3)
    real(rp),    intent(in) :: xx1(3)
    real(rp),    intent(in) :: xx2(3)
    logical(lg)             :: elmgeo_side_of_TRI03
    real(rp)                :: nn(3)
    real(rp)                :: a(3),b(3),dot1,dot2

    a(1:3)  =   bocod(1:3,2) - bocod(1:3,1)
    b(1:3)  =   bocod(1:3,3) - bocod(1:3,1)
    nn(1)   =   a(2) * b(3) - a(3) * b(2)
    nn(2)   = - a(1) * b(3) + a(3) * b(1)
    nn(3)   =   a(1) * b(2) - a(2) * b(1)
    nn(1:3) =   nn(1:3) / sqrt( dot_product(nn,nn) )

    dot1    =   dot_product(bocod(1:3,1) - xx1(1:3),nn)
    dot2    =   dot_product(bocod(1:3,1) - xx2(1:3),nn)

    if( dot1 * dot2 >= 0.0_rp ) then
       elmgeo_side_of_TRI03 = .true.
    else
       elmgeo_side_of_TRI03 = .false.
    end if

  end function elmgeo_side_of_TRI03

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   Projection on a linear surface
  !> @details Look for the nearest point XX to COGLO on a linear
  !>          surface
  !
  !-----------------------------------------------------------------------

  pure subroutine elmgeo_nearest_point_on_TRI03(bocod,coglo,xx,nn)
    real(rp),    intent(in)  :: bocod(3,4)
    real(rp),    intent(in)  :: coglo(3)
    real(rp),    intent(out) :: xx(3)
    real(rp),    intent(out) :: nn(3)
    real(rp)                 :: a(3),b(3),t

    a(1:3)  =   bocod(1:3,2) - bocod(1:3,1)
    b(1:3)  =   bocod(1:3,3) - bocod(1:3,1)
    nn(1)   =   a(2) * b(3) - a(3) * b(2)
    nn(2)   = - a(1) * b(3) + a(3) * b(1)
    nn(3)   =   a(1) * b(2) - a(2) * b(1)
    t       =   sqrt( nn(1)*nn(1) + nn(2)*nn(2) + nn(3)*nn(3) )
    nn(1)   =   nn(1)/t
    nn(2)   =   nn(2)/t
    nn(3)   =   nn(3)/t
    t       =   nn(1)*(bocod(1,1)-coglo(1)) + nn(2)*(bocod(2,1)-coglo(2)) + nn(3)*(bocod(3,1)-coglo(3))  
    xx(1)   =   coglo(1) + t * nn(1)
    xx(2)   =   coglo(2) + t * nn(2)
    xx(3)   =   coglo(3) + t * nn(3)
    
    !nn(1:3) =   nn(1:3) / sqrt( nn(1)*nn(1) + nn(2)*nn(2) + nn(3)*nn(3) )
    !t       =   dot_product(nn(1:3),bocod(1:3,1)-coglo(1:3))
    !xx(1:3) =   coglo(1:3) + t * nn(1:3)

  end subroutine elmgeo_nearest_point_on_TRI03

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   Projection on a bilinear surface
  !> @details Look for the nearest point XX to COGLO on a bilinear
  !>          surface
  !>          Minimize the distance from x to the plane:
  !>          D = [x-P(u,v)]^2 with
  !>          P(u,v) = a*u+b*v+c*u*v+d
  !>          by solving
  !>              +-     -+   +-                           -+
  !>              | dD/du |   | (a+c*v).(a*u+b*v+c*u*v+d-x) |
  !>          F:= |       | = |                             | = 0
  !>              | dD/dv |   | (b+c*u).(a*u+b*v+c*u*v+d-x) |
  !>              +-     -+   +-                           -+
  !>          using a Newton-Raphson
  !>
  !>          U^{i+1} = U^{i} -J^-1.F,  U=(u,v)^t
  !>
  !>              +-             -+
  !>              | dF1/du dF1/dv |
  !>          J = |               |
  !>              | dF2/dv dF2/dv |
  !>              +-             -+
  !>
  !>          dF1/du = (a+c*v).(a+c*v)
  !>          dF1/dv = c.(a*u+b*v+c*u*v+d-x) + (a+c*v).(b+c*u)
  !>          dF2/du = dF1/dv
  !>          dF2/dv = (b+c*u).(b+c*u)
  !>
  !>
  !-----------------------------------------------------------------------

  pure subroutine elmgeo_nearest_point_on_QUA04(bocod,coglo,xx,nn)

    real(rp),    intent(in)  :: bocod(3,4)
    real(rp),    intent(in)  :: coglo(3)
    real(rp),    intent(out) :: xx(3)
    real(rp),    intent(out) :: nn(3)
    integer(ip)              :: iiter,maxit
    real(rp)                 :: a(3),b(3),c(3),d(3),u,v,toler,denom
    real(rp)                 :: a2,ab,ac,ad,ax,bc,bd,bx,b2,c2,cd,cx,f(2)
    real(rp)                 :: dfdu(2,2),dfduinv(2,2),vec1(3),vec2(3)
    real(rp)                 :: abcdcx,fnorm
    real(rp)                 :: n1(3),n2(3),t

    dfdu(:,:) = 0.0_rp

    vec1(1:3) =   bocod(1:3,2) - bocod(1:3,1)
    vec2(1:3) =   bocod(1:3,3) - bocod(1:3,1)
    n1(1)     =   vec1(2) * vec2(3) - vec1(3) * vec2(2)
    n1(2)     = - vec1(1) * vec2(3) + vec1(3) * vec2(1)
    n1(3)     =   vec1(1) * vec2(2) - vec1(2) * vec2(1)
    n1(1:3)   =   n1(1:3) / sqrt( dot_product(n1,n1) )

    vec1(1:3) =   vec2(1:3)
    vec2(1:3) =   bocod(1:3,4) - bocod(1:3,1)
    n2(1)     =   vec1(2) * vec2(3) - vec1(3) * vec2(2)
    n2(2)     = - vec1(1) * vec2(3) + vec1(3) * vec2(1)
    n2(3)     =   vec1(1) * vec2(2) - vec1(2) * vec2(1)
    n2(1:3)   =   n2(1:3) / sqrt( dot_product(n2,n2) )

    nn        =   0.5_rp * ( n1 + n2 )

    if( maxval(abs(n1-n2)) < 1.0e-8_rp ) then
       !
       ! Point are co-planar
       !     nx*(d-x) + ny*(e-y) + nz*(f-z)
       ! t = ------------------------------
       !               n.n
       !
       t       = ( n1(1)*(bocod(1,1)-coglo(1)) + n1(2)*(bocod(2,1)-coglo(2)) + n1(3)*(bocod(3,1)-coglo(3)) )
       xx(1:3) = coglo(1:3) + t*n1(1:3)

    else
       !
       ! P(u,v) = (1-u)*[(1-v)*P00+v*P01] + u*[(1-v)*P10+v*P11]
       !
       a(1:3)    = - bocod(1:3,1) + bocod(1:3,2)
       b(1:3)    = - bocod(1:3,1) + bocod(1:3,4)
       c(1:3)    =   bocod(1:3,1) - bocod(1:3,4) - bocod(1:3,2) + bocod(1:3,3)
       d(1:3)    =   bocod(1:3,1)

       fnorm     = dot_product(vec2(1:3),vec2(1:3))

       a2        = dot_product(a,a)
       b2        = dot_product(b,b)
       ab        = dot_product(a,b)
       ac        = dot_product(a,c)
       ad        = dot_product(a,d)
       bd        = dot_product(b,d)
       bc        = dot_product(b,c)
       c2        = dot_product(c,c)
       cd        = dot_product(c,d)
       ax        = dot_product(a,coglo)
       bx        = dot_product(b,coglo)
       cx        = dot_product(c,coglo)
       abcdcx    = ab + cd - cx

       toler     = 1.0e-8_rp * fnorm
       iiter     = 0
       maxit     = 100
       denom     = 1.0_rp
       u         = 0.5_rp
       v         = 0.5_rp
       f         = 1.0_rp

       do while( maxval(abs(f(1:2))) > toler .and. iiter < maxit )
          !TODO: dfdu may be incorrectly initialized
          dfdu(2,2)    = b2 + 2.0_rp*bc*u + c2*u*u

          denom        =  dfdu(1,1) * dfdu(2,2) - dfdu(2,1) * dfdu(1,2)

          denom        =  1.0_rp / denom
          dfduinv(1,1) =  dfdu(2,2) * denom
          dfduinv(2,2) =  dfdu(1,1) * denom
          dfduinv(2,1) = -dfdu(2,1) * denom
          dfduinv(1,2) = -dfdu(1,2) * denom

          u            = u - ( dfduinv(1,1) * f(1) + dfduinv(1,2) * f(2) )
          v            = v - ( dfduinv(2,1) * f(1) + dfduinv(2,2) * f(2) )

       end do

       !print*,'caca=',iiter,abs(f(1)),abs(f(2))
       !if( iiter >= maxit ) then
       !   write(2000+kfl_paral,*)  bocod(1:3,1)
       !   write(2000+kfl_paral,*)  bocod(1:3,2)
       !   write(2000+kfl_paral,*)  bocod(1:3,3)
       !   write(2000+kfl_paral,*)  bocod(1:3,4)
       !   write(2000+kfl_paral,*)  abs(f(1)),abs(f(2))
       !   write(2000+kfl_paral,*)  coglo(1:3)
       !   flush(2000+kfl_paral)
       !   stop
       !end if
       !
       ! x = u*(-P00+P10) + v*(-P00+P01) + uv*(P00-P01-P10) + P00
       !
       xx(1:3) = a(1:3)*u + b(1:3)*v + c(1:3)*u*v + d(1:3)

    end if

  end subroutine elmgeo_nearest_point_on_QUA04

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   Check if a point is inside an element
  !> @details A ray casting strategy is being used, following the even-odd
  !>          rule for the number of intersections ray-polyhedra:
  !>          even (outside), odd (inside)
  !>
  !-----------------------------------------------------------------------

  function elmgeo_inside_element_using_faces(pelty,elcod,coglo,toler_opt)

    integer(ip), intent(in)           :: pelty
    real(rp),    intent(in)           :: elcod(3,*)
    real(rp),    intent(in)           :: coglo(3)
    real(rp),    intent(in), optional :: toler_opt
    logical(lg)                       :: elmgeo_inside_element_using_faces
    integer(ip)                       :: iface,pnodf,inodf,inode,pflty
    integer(ip)                       :: ifoun,num_intersect,pnode,idime
    integer(ip)                       :: inod1,inod2,iedge,istat,idire
    real(rp)                          :: bocod(3,16)
    real(rp)                          :: toler
    real(rp)                          :: xx(3)
    real(rp)                          :: xx2(3),segment(3),ynorm
    real(rp)                          :: xx_center(3),n1(3),maxn1,norma_segment
    real(rp)                          :: direction(3,10)
    real(rp)                          :: time1,time2
    !print *," NOU ELEMENT: " !DBG

    call cputim(time1)
    !
    ! Tolerance
    !
    if( present(toler_opt) ) then
       toler = toler_opt
    else
       toler = 1.0e-10_rp
    end if
    !
    ! Optional directions: to be used if the previos direction produces
    ! a ray intersecting an edge (first directions defined with the mass
    ! center)
    !
    direction(1,1) = 0.49755959009261719_rp
    direction(2,1) = 6.6824707611273348E-002_rp
    direction(3,1) = 0.46591537549612494_rp

    direction(1,2) = 0.24792768547143218_rp
    direction(2,2) = -0.13260910262524428_rp
    direction(3,2) = -1.9363101245268521E-002_rp

    direction(1,3) = -0.42624573636601548_rp
    direction(2,3) = -0.49464477072227275_rp
    direction(3,3) = -0.15291871148198455_rp

    direction(1,4) = -0.15775618392716495_rp
    direction(2,4) = -0.28204827366152740_rp
    direction(3,4) = -0.36683958998634070_rp

    direction(1,5) = 0.40052451442189241_rp
    direction(2,5) = -0.11323398954259645_rp
    direction(3,5) = -5.4517710612151937E-002_rp

    direction(1,6) = 0.16193218089584283_rp
    direction(2,6) = -0.48389169956944067_rp
    direction(3,6) = 0.15085483610391681_rp

    direction(1,7) = 0.14640882548382539_rp
    direction(2,7) =  -0.17701270905944422_rp
    direction(3,7) = 0.35569240288533133_rp

    direction(1,8) = -9.8713080636186112E-002_rp
    direction(2,8) = -0.29312567078124308_rp
    direction(3,8) = 0.46853946421659987_rp

    direction(1,9) = 9.8399534318134640E-002_rp
    direction(2,9) = 0.17298073277626325_rp
    direction(3,9) = -4.3117689327037967E-002_rp

    direction(1,10) = -0.16998487642662530_rp
    direction(2,10) = -0.39961707349782272_rp
    direction(3,10) = 0.25545330475972683_rp
    !
    ! Check if we are on an edge to remove pathological cases
    !
    maxn1 = 2*toler
    do iedge = 1, element_type(pelty) % number_edges
       inod1 = element_type(pelty) % list_edges(1,iedge)
       inod2 = element_type(pelty) % list_edges(2,iedge)
       n1(1) = sqrt(dot_product(coglo(1:3)-elcod(1:3,inod1),coglo(1:3)-elcod(1:3,inod1)))
       n1(2) = sqrt(dot_product(coglo(1:3)-elcod(1:3,inod2),coglo(1:3)-elcod(1:3,inod2)))
       n1(3) = sqrt(dot_product(elcod(1:3,inod1)-elcod(1:3,inod2),elcod(1:3,inod1)-elcod(1:3,inod2)))
       maxn1 = max(maxn1,n1(3))
       if( abs(n1(1)+n1(2)-n1(3))/n1(3) <= toler ) then
          elmgeo_inside_element_using_faces = .true.
          return
       end if
    end do
    !
    ! Center of gravity... check if I don't collapse, we never know!
    !
    pnode = element_type(pelty) % number_nodes
    do idime = 1,3
       xx_center(idime) = sum(elcod(idime,1:pnode)) / real(pnode,rp)
    end do
    ynorm = sqrt(dot_product(coglo-xx_center,coglo-xx_center))
    if( ynorm <= toler * n1(3) ) then
       elmgeo_inside_element_using_faces = .true.
       return
    end if
    segment = (xx_center-coglo) / ynorm
    istat=4
    idire = 0

    do while (istat==4 .and. idire < 10)
       idire = idire + 1
       !print *,"idire: ",idire !DGB
       !
       ! Choose second segment point far
       !
       ! +----------+
       ! |          |
       ! | x--c =========> xx2
       ! |          |
       ! +----------+
       !
       xx2 = coglo + segment  * maxn1 * 10.0_rp
       !
       ! Loop over faces
       !
       num_intersect = 0
       FACES: do iface = 1,element_type(pelty) % number_faces
          pnodf = element_type(pelty) % node_faces(iface)
          pflty = element_type(pelty) % type_faces(iface)

          do inodf = 1,pnodf
             inode = element_type(pelty) % list_faces(inodf,iface)
             bocod(1:3,inodf) = elcod(1:3,inode)
          end do
          !if( element_type(pelty) % number_faces == 4 .and. pnodf == 6) then
          !   call elmgeo_intersection_segment_TRI06(&
          !        bocod,coglo,xx2,toler,xx,ifoun,istat)
          !else
          call elmgeo_intersection_segment_face(&
                  3_ip,pflty,bocod,coglo,xx2,xx,ifoun,toler,istat)
          !end if

          if( istat == 4) then
             !
             ! Ray intersects edge (try with next direction)
             !
             norma_segment = 0.0_rp
             segment(:) = direction(:,idire)
             norma_segment = sqrt(dot_product(segment,segment))
             segment = segment * 1.0_rp / norma_segment
             exit FACES
          end if
          if( ifoun /= 0 .and. ( istat == 1 .or. istat == 3 ) ) then
             !
             ! COGLO is on face!
             !
             elmgeo_inside_element_using_faces = .true.
             return

          else if( ifoun > 0 ) then
             !
             ! Add number of intersections
             !
             num_intersect = num_intersect + ifoun
          end if
       end do FACES
    end do
    if( mod(num_intersect,2_ip) == 0 ) then
       elmgeo_inside_element_using_faces = .false.
    else
       elmgeo_inside_element_using_faces = .true.
    end if

    call cputim(time2)
    
  end function elmgeo_inside_element_using_faces

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   Check if a point is inside an element
  !> @details Compute the scalar product with the face normals to check
  !>          if a point is inside and element
  !>
  !-----------------------------------------------------------------------

  subroutine elmgeo_nearest_intersection_point_on_element_faces(ndime,pelty,elcod,xx1,xx2,xx_intersection,ifoun,toler_opt)
    integer(ip), intent(in)           :: ndime
    integer(ip), intent(in)           :: pelty
    real(rp),    intent(in)           :: elcod(ndime,*)
    real(rp),    intent(in)           :: xx1(ndime)
    real(rp),    intent(in)           :: xx2(ndime)
    real(rp),    intent(out)          :: xx_intersection(ndime)
    integer(ip), intent(out)          :: ifoun
    real(rp),    intent(in), optional :: toler_opt
    integer(ip)                       :: iface,pnodf,inodf,inode,pflty,kfoun
    real(rp)                          :: bocod(ndime,16),toler
    real(rp)                          :: xdist,xnorm,xx(3)

    if( present(toler_opt) ) then
       toler = toler_opt
    else
       toler = 1.0e-08_rp
    end if
    !
    ! Check
    !
    xdist = huge(1.0_rp)
    ifoun = 0
    !
    ! Loop over faces
    !
    do iface = 1,element_type(pelty) % number_faces

       pnodf = element_type(pelty) % node_faces(iface)
       pflty = element_type(pelty) % type_faces(iface)
       do inodf = 1,pnodf
          inode = element_type(pelty) % list_faces(inodf,iface)
          bocod(1:ndime,inodf) = elcod(1:ndime,inode)
       end do

       call elmgeo_intersection_segment_face(&
            ndime,pflty,bocod,xx1,xx2,xx,kfoun,toler)

       if( kfoun > 0 ) then
          ifoun = ifoun + kfoun
          xnorm = dot_product(xx1(:ndime)-xx(1:ndime),xx1(:ndime)-xx(1:ndime))
          if( xnorm < xdist ) then
             xdist = xnorm
             xx_intersection(1:ndime) = xx(1:ndime)
          end if
       end if

    end do

  end subroutine elmgeo_nearest_intersection_point_on_element_faces

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   Natural coordinates on a bilinear surface
  !> @details
  !>        4(-1,1)    3 (1,1)
  !>          o----------o
  !>          |          |
  !>          |          |
  !>          |          |
  !>          o----------o
  !>        1(-1,-1)   2 (1,-1)
  !>
  !>          +-           -+ +- -+   +-  -+
  !>          | -1 -1  1  1 | | a |   | x1 |
  !>          |  1 -1 -1  1 | | b |   | x2 |
  !>          |  1  1  1  1 | | c | = | x3 |
  !>          | -1  1 -1  1 | | d |   | x4 |
  !>          +-           -+ +- -+   +-  -+
  !>
  !>          +- -+     +-           -+ +-  -+
  !>          | a |     | -1  1  1 -1 | | x1 |
  !>          | b |   1 | -1 -1  1  1 | | x2 |
  !>          | c | = - |  1 -1  1 -1 | | x3 |
  !>          | d |   4 | -1  1  1  1 | | x4 |
  !>          +- -+     +-           -+ +-  -+
  !>
  !-----------------------------------------------------------------------

  pure subroutine elmgeo_natural_coordinates_on_QUA04(bocod,coglo,coloc)

    real(rp),    intent(in)  :: bocod(3,4)
    real(rp),    intent(in)  :: coglo(3)
    real(rp),    intent(out) :: coloc(2)
    integer(ip)              :: iiter,maxit,ii,jj,iimax(1),jjmax(1)
    real(rp)                 :: a(3),b(3),c(3),d(3),toler,denom,s,t,fnorm
    real(rp)                 :: f(2),df(2,2),dfinv(2,2),bcpy(1:3)
    !
    ! x = a*u + b*v + c*u*v + d
    !
    a(1:3) = 0.25_rp * ( - bocod(1:3,1) + bocod(1:3,2) + bocod(1:3,3) - bocod(1:3,4) )
    b(1:3) = 0.25_rp * ( - bocod(1:3,1) - bocod(1:3,2) + bocod(1:3,3) + bocod(1:3,4) )
    c(1:3) = 0.25_rp * (   bocod(1:3,1) - bocod(1:3,2) + bocod(1:3,3) - bocod(1:3,4) )
    d(1:3) = 0.25_rp * (   bocod(1:3,1) + bocod(1:3,2) + bocod(1:3,3) + bocod(1:3,4) )
    bcpy   = b
    !
    ! Choose coordinates with highest change to invert system
    !
    iimax    =  maxloc(abs(a(1:3)))
    ii       =  iimax(1)
    bcpy(ii) =  0.0_rp
    jjmax    =  maxloc(abs(bcpy(1:3)))
    jj       =  jjmax(1)
    !
    ! f(u,v) = a*u + b*v + c*u*v + d - x
    !
    fnorm  = sqrt( dot_product(bocod(1:3,3)-bocod(1:3,1),bocod(1:3,3)-bocod(1:3,1)) )
    toler  = 1.0e-08_rp * fnorm
    iiter  = 0
    maxit  = 100
    s      = 0.0_rp
    t      = 0.0_rp
    f      = 1.0_rp
    !
    ! Newton-Raphson: s^i+1 = s - [df/ds]^-1.f(s)
    !
    do while( maxval(abs(f(1:2))) > toler .and. iiter <= maxit )
       iiter      =  iiter + 1
       f(1)       =  a(ii)*s + b(ii)*t + c(ii)*s*t + d(ii) - coglo(ii)
       f(2)       =  a(jj)*s + b(jj)*t + c(jj)*s*t + d(jj) - coglo(jj)

       df(1,1)    =  a(ii) + c(ii)*t
       df(1,2)    =  b(ii) + c(ii)*s
       df(2,1)    =  a(jj) + c(jj)*t
       df(2,2)    =  b(jj) + c(jj)*s

       denom      =  df(1,1) * df(2,2) - df(2,1) * df(1,2)
       denom      =  1.0_rp  / denom
       dfinv(1,1) =  df(2,2) * denom
       dfinv(2,2) =  df(1,1) * denom
       dfinv(2,1) = -df(2,1) * denom
       dfinv(1,2) = -df(1,2) * denom

       s          =  s - ( dfinv(1,1) * f(1) + dfinv(1,2) * f(2) )
       t          =  t - ( dfinv(2,1) * f(1) + dfinv(2,2) * f(2) )

       !write(90,*) maxval(abs(f(1:2))),toler
    end do

    if( iiter > maxit ) then
       !call runend('elmgeo_natural_coordinates_on_QUA04: NEWTON-RAPHSON NOT CONVERGED')
       !stop 1
       coloc = huge(1.0_rp)
    else
       coloc(1) = s
       coloc(2) = t
    end if

  end subroutine elmgeo_natural_coordinates_on_QUA04

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   Natural coordinates on a segment
  !> @details Natural coordinates on a segment
  !>
  !-----------------------------------------------------------------------

  pure subroutine elmgeo_natural_coordinates_on_BAR02(bocod,coglo,coloc)

    real(rp),    intent(in)  :: bocod(2,2)
    real(rp),    intent(in)  :: coglo(2)
    real(rp),    intent(out) :: coloc(1)
    integer(ip)              :: ii,iimax(1)
    real(rp)                 :: uu(3)

    uu(1:2)  =  bocod(1:2,2)-bocod(1:2,1)
    iimax    =  maxloc(abs(uu(1:2)))
    ii       =  iimax(1)
    coloc(1) = -1.0_rp + 2.0_rp * ( coglo(ii)-bocod(ii,1) ) / uu(ii)

  end subroutine elmgeo_natural_coordinates_on_BAR02

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   Natural coordinates on a linear surface
  !> @details
  !>          3
  !>          o
  !>          |  \
  !>          |   \
  !>          o----o
  !>          1    2
  !>
  !>         Equation of the triangle: x = (x2-x1) * s + (x3-x1) * t + x1
  !>         One equation is redundant so take the tywo one corresponding
  !>         to the lowest normal components:
  !>         +-            -+ +- -+    +-    -+
  !>         | x2-x1  x3-x1 | | s |    | x-x1 |
  !>         |              | |   | =  |      |
  !>         | y2-y1  y3-y1 | | t |    | y-y1 |
  !>         +-            -+ +- -+    +-    -+
  !>
  !-----------------------------------------------------------------------

  pure subroutine elmgeo_natural_coordinates_on_TRI03(bocod,coglo,coloc)

    real(rp),    intent(in)  :: bocod(3,3)
    real(rp),    intent(in)  :: coglo(3)
    real(rp),    intent(out) :: coloc(2)
    integer(ip)              :: kkmax(1),ii,jj
    real(rp)                 :: mat(2,2),invmat(2,2)
    real(rp)                 :: a(3),b(3),rhs(3),deter,n(3)
    !
    ! Do not consider direction of maximum normal component
    !
    a(1:3)   =  bocod(1:3,2) - bocod(1:3,1)
    b(1:3)   =  bocod(1:3,3) - bocod(1:3,1)
    n(1)     =  a(2) * b(3) - a(3) * b(2)
    n(2)     = -a(1) * b(3) + a(3) * b(1)
    n(3)     =  a(1) * b(2) - a(2) * b(1)
    kkmax    =  maxloc(abs(n))
    ii       =  perm1(kkmax(1))
    jj       =  perm2(kkmax(1))

    mat(1,1) =  a(ii)
    mat(1,2) =  b(ii)
    rhs(1)   =  coglo(ii) - bocod(ii,1)

    mat(2,1) =  a(jj)
    mat(2,2) =  b(jj)
    rhs(2)   =  coglo(jj) - bocod(jj,1)

    call maths_invert_matrix(2_ip,mat,deter,invmat)
    !call invmtx(mat,invmat,deter,2_ip)

    coloc(1) = invmat(1,1) * rhs(1) + invmat(1,2) * rhs(2)
    coloc(2) = invmat(2,1) * rhs(1) + invmat(2,2) * rhs(2)

  end subroutine elmgeo_natural_coordinates_on_TRI03

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   Test subroutine
  !> @details This subroutine tests elmgeo_natural_coordinates_on_boundaries
  !>
  !-----------------------------------------------------------------------

  subroutine elmgeo_test_natural_coordinates_on_boundaries(&
       ndime,mnode,nelem,ltype,lnods,coord)

    integer(ip), intent(in) :: ndime
    integer(ip), intent(in) :: mnode
    integer(ip), intent(in) :: nelem
    integer(ip), intent(in) :: ltype(*)
    integer(ip), intent(in) :: lnods(mnode,*)
    real(rp),    intent(in) :: coord(ndime,*)
    integer(ip)             :: ifoun,pelty,ielem
    integer(ip)             :: iface,inodf,inode,pflty,pnodf,ii
    real(rp)                :: bocod(ndime,64),coglo2(3)
    real(rp)                :: shapb(64),derib(ndime,64)
    real(rp)                :: coglo(3),coloc(3),xdiff

    do ielem = 1,nelem

       pelty = ltype(ielem)

       do iface = 1,element_type(pelty) % number_faces
          pnodf = element_type(pelty) % node_faces(iface)
          pflty = element_type(pelty) % type_faces(iface)
          do inodf = 1,pnodf
             inode = element_type(pelty) % list_faces(inodf,iface)
             bocod(1:ndime,inodf) = coord(1:ndime,lnods(inode,ielem))
          end do

          do inodf = 1,pnodf
             inode = element_type(pelty) % list_faces(inodf,iface)
             coglo(1:ndime) = bocod(1:ndime,inodf)

             call elmgeo_natural_coordinates_on_boundaries(&
                  ndime,pflty,pnodf,bocod,shapb,derib,coglo,coloc,ifoun)

             coglo2 = 0.0_rp
             do ii = 1,pnodf
                coglo2(1:ndime) = coglo2(1:ndime) + shapb(ii)*bocod(1:ndime,ii)
             end do
             xdiff = sqrt( dot_product(coglo2(1:ndime)-coglo(1:ndime),coglo2(1:ndime)-coglo(1:ndime)) )
             if( abs(xdiff) > epsilgeo_tol ) call runend('PROBLEM WITH INTERPOLATION')
          end do

       end do
    end do

  end subroutine elmgeo_test_natural_coordinates_on_boundaries

  !-----------------------------------------------------------------------
  !> @author  Guillaume Houzeaux
  !> @date    27/01/2007
  !> @brief   Computes the Jacobian, Jacobian determinant and
  !           Cartesian derivatives of shape function of an element
  !> @details The jacobian is
  !                              _           _
  !                             | dx/ds dx/dt |                t
  !           Jacobian: XJACM = |             | = ELCOD * DERIV
  !                             |_dy/ds dy/dt_|
  !           with
  !                   _        _             _                    _
  !                  | x1 x2 x3 |           | dN1/ds dN2/ds dN3/ds |
  !          ELCOD = |          |,  DERIV = |                      |
  !                  |_y1 y2 y3_|           |_dN1/dt dN2/dt dN3/dt_|
  !
  !           => Jacobian determinant: GPDET = det(XJACM)
  !-----------------------------------------------------------------------

  pure subroutine elmgeo_cartesian_derivatives(ndime,pnode,elcod,deriv,gpcar,gpdet_opt)

    integer(ip), intent(in)            :: ndime,pnode
    real(rp),    intent(in)            :: elcod(ndime,pnode),deriv(ndime,pnode)
    real(rp),    intent(out)           :: gpcar(ndime,pnode)
    real(rp),    intent(out), optional :: gpdet_opt
    integer(ip)                        :: j,k
    real(rp)                           :: t1,t2,t3,denom,gpdet
    real(rp)                           :: xjacm(ndime,ndime),xjaci(ndime,ndime)

    if( ndime == 2 .and. pnode == 3 ) then
       !
       ! 2D P1 element
       !
       gpdet = (-elcod(1,1)+elcod(1,2))*(-elcod(2,1)+elcod(2,3)) &
            & -(-elcod(2,1)+elcod(2,2))*(-elcod(1,1)+elcod(1,3))
       if( gpdet == 0.0_rp ) goto 10
       denom = 1.0_rp/gpdet

       gpcar(1,1) = ( -elcod(2,3) + elcod(2,2) ) * denom
       gpcar(1,2) = ( -elcod(2,1) + elcod(2,3) ) * denom
       gpcar(1,3) = (  elcod(2,1) - elcod(2,2) ) * denom
       gpcar(2,1) = (  elcod(1,3) - elcod(1,2) ) * denom
       gpcar(2,2) = (  elcod(1,1) - elcod(1,3) ) * denom
       gpcar(2,3) = ( -elcod(1,1) + elcod(1,2) ) * denom

    else if( ndime == 3 .and. pnode == 4 ) then
       !
       ! 3D P1 element
       !
       gpcar(1,1) =  elcod(1,2) - elcod(1,1)
       gpcar(1,2) =  elcod(1,3) - elcod(1,1)
       gpcar(1,3) =  elcod(1,4) - elcod(1,1)
       gpcar(2,1) =  elcod(2,2) - elcod(2,1)
       gpcar(2,2) =  elcod(2,3) - elcod(2,1)
       gpcar(2,3) =  elcod(2,4) - elcod(2,1)
       gpcar(3,1) =  elcod(3,2) - elcod(3,1)
       gpcar(3,2) =  elcod(3,3) - elcod(3,1)
       gpcar(3,3) =  elcod(3,4) - elcod(3,1)
       t1         =  gpcar(2,2) * gpcar(3,3) - gpcar(3,2) * gpcar(2,3)
       t2         = -gpcar(2,1) * gpcar(3,3) + gpcar(3,1) * gpcar(2,3)
       t3         =  gpcar(2,1) * gpcar(3,2) - gpcar(3,1) * gpcar(2,2)
       gpdet      =  gpcar(1,1) * t1 + gpcar(1,2) * t2 + gpcar(1,3) * t3
       if( gpdet == 0.0_rp ) goto 10
       denom      =  1.0_rp/gpdet

       xjaci(1,1) =  t1 * denom
       xjaci(2,1) =  t2 * denom
       xjaci(3,1) =  t3 * denom
       xjaci(2,2) = ( gpcar(1,1) * gpcar(3,3) - gpcar(3,1) * gpcar(1,3)) * denom
       xjaci(3,2) = (-gpcar(1,1) * gpcar(3,2) + gpcar(1,2) * gpcar(3,1)) * denom
       xjaci(3,3) = ( gpcar(1,1) * gpcar(2,2) - gpcar(2,1) * gpcar(1,2)) * denom
       xjaci(1,2) = (-gpcar(1,2) * gpcar(3,3) + gpcar(3,2) * gpcar(1,3)) * denom
       xjaci(1,3) = ( gpcar(1,2) * gpcar(2,3) - gpcar(2,2) * gpcar(1,3)) * denom
       xjaci(2,3) = (-gpcar(1,1) * gpcar(2,3) + gpcar(2,1) * gpcar(1,3)) * denom
       gpcar(1,1) = -xjaci(1,1)-xjaci(2,1)-xjaci(3,1)
       gpcar(1,2) =  xjaci(1,1)
       gpcar(1,3) =  xjaci(2,1)
       gpcar(1,4) =  xjaci(3,1)
       gpcar(2,1) = -xjaci(1,2)-xjaci(2,2)-xjaci(3,2)
       gpcar(2,2) =  xjaci(1,2)
       gpcar(2,3) =  xjaci(2,2)
       gpcar(2,4) =  xjaci(3,2)
       gpcar(3,1) = -xjaci(1,3)-xjaci(2,3)-xjaci(3,3)
       gpcar(3,2) =  xjaci(1,3)
       gpcar(3,3) =  xjaci(2,3)
       gpcar(3,4) =  xjaci(3,3)

    else if ( ndime == 1 ) then
       !
       ! 1D
       !
       xjacm(1,1) = 0.0_rp
       do k = 1,pnode
          xjacm(1,1) = xjacm(1,1) + elcod(1,k) * deriv(1,k)
       end do
       gpdet = xjacm(1,1)
       if( gpdet == 0.0_rp ) goto 10
       xjaci(1,1) = 1.0_rp/xjacm(1,1)
       do j = 1,pnode
          gpcar(1,j) = xjaci(1,1) * deriv(1,j)
       end do

    else if ( ndime == 2 ) then
       !
       ! 2D
       !
       xjacm(1,1) = 0.0_rp
       xjacm(1,2) = 0.0_rp
       xjacm(2,1) = 0.0_rp
       xjacm(2,2) = 0.0_rp
       do k = 1,pnode
          xjacm(1,1) = xjacm(1,1) + elcod(1,k) * deriv(1,k)
          xjacm(1,2) = xjacm(1,2) + elcod(1,k) * deriv(2,k)
          xjacm(2,1) = xjacm(2,1) + elcod(2,k) * deriv(1,k)
          xjacm(2,2) = xjacm(2,2) + elcod(2,k) * deriv(2,k)
       end do

       gpdet = xjacm(1,1) * xjacm(2,2) - xjacm(2,1) * xjacm(1,2)
       if( gpdet == 0.0_rp ) goto 10
       denom = 1.0_rp/gpdet
       xjaci(1,1) =  xjacm(2,2) * denom
       xjaci(2,2) =  xjacm(1,1) * denom
       xjaci(2,1) = -xjacm(2,1) * denom
       xjaci(1,2) = -xjacm(1,2) * denom

       do j = 1, pnode
          gpcar(1,j) =   xjaci(1,1) * deriv(1,j) &
               &       + xjaci(2,1) * deriv(2,j)

          gpcar(2,j) =   xjaci(1,2) * deriv(1,j) &
               &       + xjaci(2,2) * deriv(2,j)
       end do

    else if ( ndime == 3 ) then
       !
       ! 3D
       !
       xjacm(1,1) = 0.0_rp ! xjacm = elcod * deriv^t
       xjacm(1,2) = 0.0_rp ! xjaci = xjacm^-1
       xjacm(1,3) = 0.0_rp ! gpcar = xjaci^t * deriv
       xjacm(2,1) = 0.0_rp
       xjacm(2,2) = 0.0_rp
       xjacm(2,3) = 0.0_rp
       xjacm(3,1) = 0.0_rp
       xjacm(3,2) = 0.0_rp
       xjacm(3,3) = 0.0_rp
       do k = 1,pnode
          xjacm(1,1) = xjacm(1,1) + elcod(1,k) * deriv(1,k)
          xjacm(1,2) = xjacm(1,2) + elcod(1,k) * deriv(2,k)
          xjacm(1,3) = xjacm(1,3) + elcod(1,k) * deriv(3,k)
          xjacm(2,1) = xjacm(2,1) + elcod(2,k) * deriv(1,k)
          xjacm(2,2) = xjacm(2,2) + elcod(2,k) * deriv(2,k)
          xjacm(2,3) = xjacm(2,3) + elcod(2,k) * deriv(3,k)
          xjacm(3,1) = xjacm(3,1) + elcod(3,k) * deriv(1,k)
          xjacm(3,2) = xjacm(3,2) + elcod(3,k) * deriv(2,k)
          xjacm(3,3) = xjacm(3,3) + elcod(3,k) * deriv(3,k)
       end do

       t1    =  xjacm(2,2) * xjacm(3,3) - xjacm(3,2) * xjacm(2,3)
       t2    = -xjacm(2,1) * xjacm(3,3) + xjacm(3,1) * xjacm(2,3)
       t3    =  xjacm(2,1) * xjacm(3,2) - xjacm(3,1) * xjacm(2,2)
       gpdet =  xjacm(1,1) * t1 + xjacm(1,2) * t2 + xjacm(1,3) * t3
       if(gpdet == 0.0_rp ) goto 10
       denom = 1.0_rp / gpdet
       xjaci(1,1) = t1*denom
       xjaci(2,1) = t2*denom
       xjaci(3,1) = t3*denom
       xjaci(2,2) = ( xjacm(1,1) * xjacm(3,3) - xjacm(3,1) * xjacm(1,3)) * denom
       xjaci(3,2) = (-xjacm(1,1) * xjacm(3,2) + xjacm(1,2) * xjacm(3,1)) * denom
       xjaci(3,3) = ( xjacm(1,1) * xjacm(2,2) - xjacm(2,1) * xjacm(1,2)) * denom
       xjaci(1,2) = (-xjacm(1,2) * xjacm(3,3) + xjacm(3,2) * xjacm(1,3)) * denom
       xjaci(1,3) = ( xjacm(1,2) * xjacm(2,3) - xjacm(2,2) * xjacm(1,3)) * denom
       xjaci(2,3) = (-xjacm(1,1) * xjacm(2,3) + xjacm(2,1) * xjacm(1,3)) * denom

       do j = 1, pnode
          gpcar(1,j) =   xjaci(1,1) * deriv(1,j) &
               &       + xjaci(2,1) * deriv(2,j) &
               &       + xjaci(3,1) * deriv(3,j)

          gpcar(2,j) =   xjaci(1,2) * deriv(1,j) &
               &       + xjaci(2,2) * deriv(2,j) &
               &       + xjaci(3,2) * deriv(3,j)

          gpcar(3,j) =   xjaci(1,3) * deriv(1,j) &
               &       + xjaci(2,3) * deriv(2,j) &
               &       + xjaci(3,3) * deriv(3,j)
       end do

    end if
    
10 continue
    if( present(gpdet_opt) ) gpdet_opt = gpdet
    
  end subroutine elmgeo_cartesian_derivatives

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   Intersection segement-face
  !> @details Compute the intersection point between a segment (xx1,xx2)
  !>          and a face
  !>          IFOUN =  0 ... Segment does not intersect with element
  !>                =  1 ... Segment intersects element once
  !>                =  2 ... Segment intersects element twice
  !>          ISTAT =  0 ... None of intersection points is on interface
  !>                =  1 ... First point xx1 on surface
  !>                =  2 ... Second point xx2 on surface
  !>                =  3 ... Both points
  !>                =  4 ... intersection point inside an edge and differnt
  !>                         than first and second point
  !
  !-----------------------------------------------------------------------

  subroutine elmgeo_intersection_segment_face(&
       ndime,pflty,bocod,xx1,xx2,intersection,ifoun,toler,istat)
    integer(ip), intent(in)            :: ndime                   !< Problem dimension
    integer(ip), intent(in)            :: pflty                   !< Face type
    real(rp),    intent(in)            :: bocod(ndime,*)          !< Face node coordinates
    real(rp),    intent(in)            :: xx1(ndime)              !< Segment 1st point
    real(rp),    intent(in)            :: xx2(ndime)              !< Segment 2nd point
    real(rp),    intent(out)           :: intersection(ndime)     !< Coordinates of the intersection point
    integer(ip), intent(out)           :: ifoun                   !< Intersects or not
    real(rp),    intent(in),  optional :: toler                   !< Tolerance
    integer(ip), intent(out), optional :: istat                   !< Stat on intersection points
    real(rp)                           :: eps

    if( present(toler) ) then
       eps = toler
    else
       eps = 0.0_rp
    end if

    ifoun = 0
    if( present(istat) ) istat = 0

    if( pflty == BAR02 ) then
       !
       ! BAR02
       !
       call elmgeo_intersection_segment_BAR02(bocod,xx1,xx2,eps,intersection,ifoun)
       !print *, " Interesect BAR02 ",ifoun !DBG

    else if( pflty == TRI03 ) then
       !
       ! TRI03
       !
       call elmgeo_intersection_segment_TRI03(bocod,xx1,xx2,eps,intersection,ifoun,istat)
       !print *, " Interesect TRI03 ",ifoun !DBG

    else if( pflty == QUA04 ) then
       !
       ! QUA04
       !
       call elmgeo_intersection_segment_QUA04(bocod,xx1,xx2,eps,intersection,ifoun,istat)
       !print *, " Interesect QUA04 ",ifoun !DBG

    else if( pflty == BAR03 ) then
      !
      ! BAR03
      !
      call elmgeo_intersection_segment_BAR03(bocod,xx1,xx2,eps,intersection,ifoun)
    
    else if( pflty == TRI06 ) then
       !
       ! TRI06
       !
       call elmgeo_intersection_segment_TRI06(bocod,xx1,xx2,eps,intersection,ifoun,istat)

    else

       call runend('ELMGEO_INTERSECTION_ON_A_FACE: UNKNOWN TYPE OF FACE')

    end if

  end subroutine elmgeo_intersection_segment_face

  !-----------------------------------------------------------------------
  ! 
  !> @author  Edgar Olivares
  !> @brief   Intersection quadratic segment: BAR03
  !> @details Compute the intersection point between a ray and a
  !>          quadratic segment
  !
  !-----------------------------------------------------------------------
  subroutine elmgeo_intersection_segment_BAR03(bocod,xx1,xx2,toler,intersection1,ifoun)
    real(rp),    intent(in)            :: bocod(2,3)          !> FACE NODE COORDINATES
    real(rp),    intent(in)            :: xx1(2)
    real(rp),    intent(in)            :: xx2(2)
    real(rp),    intent(in)            :: toler
    real(rp),    intent(out)           :: intersection1(2)
    integer(ip), intent(out)           :: ifoun
    real(rp)                           :: r(2),L(2),M(2),N(2),intersection2(2) 
    real(rp)                           :: lmini,lmaxi,a,b,c,delta,s1,s2,t1,t2
    real(rp)                           :: shapf(3),deriv(1,3),rr(2),tangent

    ifoun = 0
    lmini = -toler
    lmaxi = 1.0_rp + toler

    r = xx2 - xx1
    L = xx1 - bocod(:,3)
    M = (bocod(:,1) + bocod(:,2))/2.0_rp
    N = (bocod(:,1) + bocod(:,2))/2.0_rp - bocod(:,3)
    ! 2nd order equation values (a,b,c) such as ax^2+bx+c = 0
    a =  N(2) - r(2)/(r(1)+toler)*N(1)
    b =  M(2) + r(2)/(r(1)+toler)*M(1)
    c = -L(2) + r(2)/(r(1)+toler)*L(1)
    delta = b**2-4.0_rp*a*c
    if( a /= 0.0_rp .and. delta >= 0.0_rp )then
      s1 = -b+sqrt(delta)
      s1 = s1/(2.0_rp*a)
      s2 = -b-sqrt(delta)
      s2 = s2/(2.0_rp*a)
    end if
    if( a==0.0_rp ) then
      s1 = -c/b
      s2 = -c/b
    end if
    if( delta < 0.0_rp ) then
      ! No solution => no intersection
      s1 = 1000.0_rp
      s2 = 1000.0_rp
    end if

    ! Calculate intersections if 's' is defined inside the element
    if( s1 >= lmini .and. s1 <= lmaxi ) then
      t1 = -( L(1) + s1*M(1)-s1**2*N(1) )
      intersection1 = xx1 + r*t1
      ! Let's check if the intersection is tangencial.
      ! In such case, we should not count this as a valid interseciton.
      call elmgeo_shapf_deriv_heslo(1_ip,3_ip,[s1],shapf,deriv)

      !call elmgeo_shape1(s1,3_ip,shapf,deriv)
      rr = deriv(1,1)*bocod(:,1) + deriv(1,2)*bocod(:,2) + deriv(1,3)*bocod(:,3)
      tangent = r(2)*rr(1) - r(1)*rr(2)
      if( tangent > toler .or. tangent < -toler ) ifoun = ifoun + 1
    end if
    if( s2 >= lmini .and. s2 <= lmaxi ) then
      t2 = -( L(1) + s2*M(1)-s2**2*N(1) )
      intersection2 = xx1 + r*t1
      ! Let's check if the intersection is tangencial.
      ! In such case, we should not count this as a valid interseciton.
      call elmgeo_shapf_deriv_heslo(1_ip,3_ip,[s2],shapf,deriv)
      !call elmgeo_shape1(s2,3_ip,shapf,deriv)
      rr(:) = deriv(1,1)*bocod(:,1) + deriv(1,2)*bocod(:,2) + deriv(1,3)*bocod(:,3)
      tangent = r(2)*rr(1) - r(1)*rr(2)
      if( tangent > toler .or. tangent < -toler ) ifoun = ifoun + 1
    end if 
 
  end subroutine elmgeo_intersection_segment_BAR03

  ! >> UNDER CONSTRUCTION <<
  !-----------------------------------------------------------------------
  !
  !> @author  Edgar Olivares
  !> @brief   Intersection quadratic face: TRI06
  !> @details Compute the intersection point between a segment (xx1,xx2)
  !
  !-----------------------------------------------------------------------
  subroutine elmgeo_intersection_segment_TRI06(bocod,xx1,xx2,toler,intersection,ifoun,istat)
    real(rp),    intent(in)            :: bocod(3,6)          !> FACE NODE COORDINATES
    !real(rp),    intent(in)            :: elcod(3,10)
    real(rp),    intent(in)            :: xx1(3)
    real(rp),    intent(in)            :: xx2(3)
    real(rp),    intent(in)            :: toler
    real(rp),    intent(out)           :: intersection(3)
    integer(ip), intent(out)           :: ifoun
    integer(iP), intent(out), optional :: istat
    real(rp)                           :: u(3,8),v(3,8),n(3),deriv(3,6)
    real(rp)                           :: r(3),w(3),coloc(3),shapf(10)
    real(rp)                           :: lmini,lmaxi
    real(rp)                           :: normn,nn
    real(rp)                           :: u_dot_u,v_dot_v
    real(rp)                           :: u_dot_w,v_dot_w
    real(rp)                           :: u_dot_v,s,t,z,a,b
    integer(ip)                        :: jstat,i

    ifoun  = 0
    jstat  = 0
    lmini  = -toler
    lmaxi  = 1.0_rp + toler

    r      = xx2 - xx1

    ! NR approach
    if( 0 == 1 )then

       call elmgeo_newrap_norm(r,xx1,coloc,3_ip,6_ip,bocod,shapf,deriv)

       intersection = xx1 + coloc(3)*r
       s = coloc(1)
       z = coloc(2)
       t = coloc(3)

       if( s >= lmini .and. s <= lmaxi ) then
          if( z >= lmini .and. s+z <= lmaxi ) then
             ifoun = 1
             if( abs(t) <= toler ) then
                jstat = 1
             else if( abs(t-1.0_rp) <= toler ) then
                jstat = 2
             else if( abs(s) <= toler .or. abs(z) <= toler .or. abs(s+t-1.0_rp)<= toler) then
                jstat = 4
             end if
          end if
       end if
    ! Multiple-triangles apporach
    else

       u(:,1) = bocod(:,4) - bocod(:,1)
       v(:,1) = bocod(:,6) - bocod(:,1)

       u(:,2) = bocod(:,4) - bocod(:,2)
       v(:,2) = bocod(:,5) - bocod(:,2)

       u(:,3) = bocod(:,5) - bocod(:,3)
       v(:,3) = bocod(:,6) - bocod(:,3)

       u(:,4) = bocod(:,5) - bocod(:,4)
       v(:,4) = bocod(:,6) - bocod(:,4)

       !r      = xx2 - xx1

       !
       !> Loop over the 4 sub-triangles building the main one
       !
       TRIANGLES: do i = 1,4 !8
          n(1)   = u(2,i)*v(3,i) - u(3,i)*v(2,i)
          n(2)   = u(3,i)*v(1,i) - u(1,i)*v(3,i)
          n(3)   = u(1,i)*v(2,i) - u(2,i)*v(1,i)
          normn    = sqrt(dot_product(n(:),n(:)))
          nn       = 1.0_rp/normn
          if( normn < epsilgeo_alg) then
             !
             ! Quadratic triangle degenerates
             !
             call runend("elmgeo_intersection_segment_TRI06: quadratic triangle degenerates")

          else
             n(1:3) = n(1:3) * nn
             b        = dot_product(n(:),r)              * sqrt(nn)
             a        = dot_product(n(:),bocod(:,i)-xx1) * sqrt(nn)
             !if( i > 6 ) then
             !  a = dot_product(n(:),bocod(:,3)-xx1) * sqrt(nn)
             !else
             !  a = dot_product(n(:),bocod(:,4)-xx1) * sqrt(nn)
             !end if

             if( abs(b) < epsilgeo_alg .and. abs(a) < epsilgeo_alg ) then
                call runend('elmgeo_intersection_segment_TRI06: CODE THIS ESPECIE DE VAGO!')

             else
                t = a/b
                if( t >= lmini .and. t <= lmaxi ) then
                   !
                   ! Segment intersects with plane
                   !
                   intersection = xx1 + t * r
                   w            = intersection - bocod(:,i)
                   !if( i > 6 ) then
                   !   w = intersection - bocod(:,3)
                   !else
                   !   w = intersection - bocod(:,4)
                   !end if

                   u_dot_u      = dot_product(u(:,i),u(:,i))
                   v_dot_v      = dot_product(v(:,i),v(:,i))
                   u_dot_v      = dot_product(u(:,i),v(:,i))
                   u_dot_w      = dot_product(u(:,i),w(:))
                   v_dot_w      = dot_product(v(:,i),w(:))
                   a            = 1.0_rp / ( u_dot_v*u_dot_v - u_dot_u*v_dot_v )
                   s            = ( u_dot_v*v_dot_w - v_dot_v*u_dot_w ) * a

                   if( s >= lmini .and. s <= lmaxi ) then
                      z = ( u_dot_v*u_dot_w - u_dot_u*v_dot_w ) * a

                      if( z >= lmini .and. s+z <= lmaxi ) then

                         ifoun = 1
                         if( abs(t) <= toler ) then
                            jstat = 1
                         else if( abs(t-1.0_rp) <= toler ) then
                            jstat = 2
                         else if( abs(s) <= toler .or. abs(z) <= toler .or. abs(s+t-1.0_rp)<= toler) then
                            jstat = 4
                         end if
                         exit TRIANGLES
                      end if
                   end if
                end if
             end if
          end if

       end do TRIANGLES
    end if
    if( present(istat) ) istat = jstat

  end subroutine elmgeo_intersection_segment_TRI06

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   Intersection segement-BAR02
  !> @details Compute the intersection point between a segment (xx1,xx2)
  !>          and a BAR02 element
  !>          IFOUN = 0 ... Segment does not intersect with element
  !>                = 1 ... Segment intersects element
  !>                = 2 ... Segment overlaps element
  !>
  !>          \verbatim
  !>
  !>          q+s
  !>          _
  !>         | \     p+r
  !>            \    _
  !>             \   /|
  !>              \ /
  !>               o  p + t*r = q + u*s
  !>              / |
  !>             /   |
  !>            /     |
  !>          p=xx1  q=bocod(:,1)
  !>
  !>          \endverbatim
  !
  !-----------------------------------------------------------------------

  subroutine elmgeo_intersection_segment_BAR02(&
       bocod,xx1,xx2,toler,intersection,ifoun)
    real(rp),    intent(in)   :: bocod(2,2)      !< Face node coordinates
    real(rp),    intent(in)   :: xx1(2)          !< Segment 1st point
    real(rp),    intent(in)   :: xx2(2)          !< Segment 2nd point
    real(rp),    intent(in)   :: toler           !< Tolerance
    real(rp),    intent(out)  :: intersection(2) !< Coordinates of the intersection point
    integer(ip), intent(out)  :: ifoun           !< Intersects or not
    real(rp)                  :: r(2)
    real(rp)                  :: s(2)
    real(rp)                  :: qmp(2)
    real(rp)                  :: r_cross_s
    real(rp)                  :: qmp_cross_r
    real(rp)                  :: qmp_cross_s
    real(rp)                  :: t0,t1,t,u
    real(rp)                  :: tmin,tmax
    real(rp)                  :: lmini,lmaxi
    real(rp)                  :: r_dot_r,s_dot_s
    real(rp)                  :: rzero

    ifoun       = 0
    lmini       = -toler
    lmaxi       = 1.0_rp + toler
    r(1:2)      = xx2(1:2)     - xx1(1:2)
    s(1:2)      = bocod(1:2,2) - bocod(1:2,1)
    qmp(1:2)    = bocod(1:2,1) - xx1(1:2)

    r_dot_r     = dot_product(r,r)
    s_dot_s     = dot_product(s,s)
    r_cross_s   = r(1)*s(2)   - r(2)*s(1)
    qmp_cross_r = qmp(1)*r(2) - qmp(2)*r(1)

    rzero = min(epsilgeo_alg, epsilgeo_alg*r_dot_r*s_dot_s)

    if(r_dot_r == 0_rp) then
       call runend("elmgeo_intersection_segment_BAR02: segment with norm 0")
    endif
    if(s_dot_s == 0_rp) then
       call runend("BAR02 with norm 0")
    endif

    if( abs(r_cross_s) < rzero .and. abs(qmp_cross_r) < rzero ) then
       !
       ! Two lines are colinear
       !
       r_dot_r = dot_product(r,r)
       t0      = dot_product(qmp,r)   / r_dot_r
       t1      = dot_product(qmp+s,r) / r_dot_r
       tmin    = min(t0,t1)
       tmax    = max(t0,t1)
       if( tmin <= lmaxi .and. tmax >= lmini ) then
          ifoun = 2
          t     = 0.5_rp * (max(0.0_rp,tmin)+min(1.0_rp,tmax))
       end if

    else if( abs(r_cross_s) < rzero .and. abs(qmp_cross_r) >= rzero ) then
       !
       ! Lines are parallel but non-intersecting
       !
       continue

    else
       !
       ! Lines intersect
       !
       qmp_cross_s = qmp(1)*s(2) - qmp(2)*s(1)
       t           = qmp_cross_s / r_cross_s
       u           = qmp_cross_r / r_cross_s
       if(     t >= -toler .and. t <= 1.0_rp + toler .and. &
            &  u >= -toler .and. u <= 1.0_rp + toler ) then
          ifoun = 1
       end if

    end if

    if( ifoun > 0 ) then
       intersection = xx1 + t * r
    else
       intersection = 0.0_rp
    end if

  end subroutine elmgeo_intersection_segment_BAR02

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Hou/zeaux
  !> @brief   Intersection segement-BAR02
  !> @details Compute the intersection point between a segment (xx1,xx2)
  !>          and a TRI03 element
  !>          IFOUN =  0 ... Segment does not intersect with element
  !>                =  1 ... Segment intersects element
  !>                =  2 ... Segment overlaps element
  !>          ISTAT =  0 ... None of intersection points is on interface
  !>                =  1 ... First point xx1 on surface
  !>                =  2 ... Second point xx2 on surface
  !>                =  3 ... Both points: not coded yet!
  !>                =  4 ... intersection point inside an edge and differnt
  !>                         than first and second point
  !
  !-----------------------------------------------------------------------
  subroutine elmgeo_intersection_segment_TRI03(&
       bocod,xx1,xx2,toler,intersection,ifoun,istat)
    real(rp),    intent(in)            :: bocod(3,3)      !< Face node coordinates
    real(rp),    intent(in)            :: xx1(3)          !< Segment 1st point
    real(rp),    intent(in)            :: xx2(3)          !< Segment 2nd point
    real(rp),    intent(in)            :: toler           !< Tolerance
    real(rp),    intent(out)           :: intersection(3) !< Coordinates of the intersection point
    integer(ip), intent(out)           :: ifoun           !< Intersects or not
    integer(ip), intent(out), optional :: istat           !< Intersects or not
    real(rp)                           :: u(3)
    real(rp)                           :: v(3)
    real(rp)                           :: w(3)
    real(rp)                           :: n(3)
    real(rp)                           :: r(3)
    real(rp)                           :: lmini,lmaxi
    real(rp)                           :: u_dot_u,v_dot_v
    real(rp)                           :: u_dot_w,v_dot_w,nn,normn
    real(rp)                           :: u_dot_v,s,t,z,a,b
    integer(ip)                        :: jstat

    ifoun =  0
    lmini =  -toler
    lmaxi =  1.0_rp + toler
    u     =  bocod(:,2) - bocod(:,1)
    v     =  bocod(:,3) - bocod(:,1)
    r     =  xx2        - xx1
    n(1)  =  u(2)*v(3)  - u(3)*v(2)
    n(2)  =  u(3)*v(1)  - u(1)*v(3)
    n(3)  =  u(1)*v(2)  - u(2)*v(1)          ! n=[m^2]

    normn =  sqrt(dot_product(n,n))
    jstat = 0

    nn    =  1.0_rp / normn
    if( normn < epsilgeo_alg) then
       !
       ! Triangle degenerates
       !
       call runend("elmgeo_intersection_segment_TRI03: triangle degenerates")

    else
       n     =  n * nn
       b     =  dot_product(n,r) * sqrt(nn)
       a     =  dot_product(n,bocod(:,1)-xx1) * sqrt(nn)

       if( abs(b) < epsilgeo_alg ) then
          !
          ! Segment parallel to triangle plane
          !
          if( abs(a) < epsilgeo_alg ) then
             ifoun = 2 ! Segment lies in triangle plane
             jstat = 4
             call runend('elmgeo_intersection_segment_TRI03: CODE THIS ESPECIE DE VAGO!')
          else
             ifoun = 0 ! Segment disjoint from plane
          end if
       else
          t = a / b
          if( t >= lmini .and. t <= lmaxi ) then
             !
             ! Segment intersects with plane
             !
             intersection = xx1 + t * r
             w            = intersection - bocod(:,1)
             u_dot_u      = dot_product(u,u)
             v_dot_v      = dot_product(v,v)
             u_dot_v      = dot_product(u,v)
             u_dot_w      = dot_product(u,w)
             v_dot_w      = dot_product(v,w)
             a            = 1.0_rp / ( u_dot_v*u_dot_v - u_dot_u*v_dot_v )
             s            = ( u_dot_v*v_dot_w - v_dot_v*u_dot_w ) * a
             if( s >= lmini .and. s <= lmaxi ) then
                z = ( u_dot_v*u_dot_w - u_dot_u*v_dot_w ) * a
                if( z >= lmini .and. s+z <= lmaxi ) then
                   ifoun = 1
                   if( abs(t) <= toler ) then
                      jstat = 1
                   else if( abs(t-1.0_rp) <= toler ) then
                      jstat = 2
                   else if( abs(z) <= toler .or. abs(s) <= toler .or. abs(s+z-1.0_rp)<= toler) then
                      jstat = 4
                   end if
                end if
             end if
          end if
       end if
    end if

    if( present(istat) ) istat = jstat

  end subroutine elmgeo_intersection_segment_TRI03

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   Intersection segement-BAR02
  !> @details Compute the intersection point between a segment (xx1,xx2)
  !>          and a QUA04 element. See
  !>          S.D. Ramsey, K. Potter and C. Hansen.
  !>          "Ray Bilinear Patch Intersections"
  !>          IFOUN =  0 ... Segment does not intersect with element
  !>                =  1 ... Segment intersects element once
  !>                =  2 ... Segment intersects element twice
  !>          ISTAT =  0 ... None of intersection points is on interface
  !>                =  1 ... First point xx1 on surface
  !>                =  2 ... Second point xx2 on surface
  !>                =  3 ... Both points
  !>                =  4 ... intersection point inside an edge and differnt
  !>                         than first and second point
  !>          \verbatim
  !>
  !>           /| v
  !>            |
  !>           p01           p11
  !>            +------------+
  !>            |            |
  !>            |            |
  !>            |            |
  !>            +------------+-----> u
  !>           p00           p10
  !>
  !>
  !>            Bilinear patch equation: u*v*a + u*b + v*c + d = 0   !
  !>            Ray: p(t) = r + t*q
  !>
  !>
  !>            We have 3 equations with 3 unknowns
  !>
  !>            qx * t = u*v * ax + u * bx + v * cx + dx - rx
  !>            qy * t = u*v * ay + u * by + v * cy + dy - ry
  !>            qz * t = u*v * az + u * bz + v * cz + dz - rz
  !>
  !>            Eliminate t from the equation where q is max
  !>
  !>
  !>          \endverbatim
  !>
  !-----------------------------------------------------------------------

  subroutine elmgeo_intersection_segment_QUA04(&
       bocod,xx1,xx2,toler,intersection,ifoun,istat)
    real(rp),    intent(in)            :: bocod(3,4)      !< Face node coordinates
    real(rp),    intent(in)            :: xx1(3)          !< Segment 1st point
    real(rp),    intent(in)            :: xx2(3)          !< Segment 2nd point
    real(rp),    intent(in)            :: toler           !< Tolerance
    real(rp),    intent(out)           :: intersection(3) !< Coordinates of the intersection point
    integer(ip), intent(out)           :: ifoun           !< Intersects or not
    integer(ip), intent(out), optional :: istat           !< Stat about intersection points
    real(rp)                           :: a(3)
    real(rp)                           :: b(3)
    real(rp)                           :: c(3)
    real(rp)                           :: d(3)
    real(rp)                           :: q(3)
    real(rp)                           :: r(3)
    real(rp)                           :: A1,B1,C1,D1
    real(rp)                           :: A2,B2,C2,D2
    real(rp)                           :: aa,bb,cc,delta,e1(3),e2(3)
    real(rp)                           :: v1,v2,p1(3),p2(3)
    real(rp)                           :: t1,t2,u1,u2,xnorm
    integer(ip)                        :: jstat,imax(1)
    logical(lg)                        :: patch1,patch2, inedge
    real(rp)                           :: toler2          ! Tolerance to determine that a ray is not convenient

    toler2 = 0.001_rp
    !
    ! /| v
    !  |
    ! p01           p11
    !  +------------+
    !  |            |
    !  |            |
    !  |            |
    !  +------------+-----> u
    ! p00           p10
    !
    ! Compute max diagonal norm to adimesionalize the equations
    !
    e1(1:3) = bocod(1:3,1) - bocod(1:3,3)                               ! e1 = p00-p11
    e2(1:3) = bocod(1:3,2) - bocod(1:3,4)                               ! e2 = p10-p01
    xnorm   = 1.0_rp / max(dot_product(e1,e1),dot_product(e2,e2))       !    = 1 / max(e1.e1,e2.e2) max diagonal
    !
    ! Bilinear patch equation: u*v*a + u*b + v*c + d = 0
    !
    a(1:3)  = bocod(1:3,3) - bocod(1:3,2) - bocod(1:3,4) + bocod(1:3,1) ! a  = p11-p10-p01+p00
    b(1:3)  = bocod(1:3,2) - bocod(1:3,1)                               ! b  = p10-p00
    c(1:3)  = bocod(1:3,4) - bocod(1:3,1)                               ! c  = p01-p00
    d(1:3)  = bocod(1:3,1)                                              ! d  = p00
    !
    ! Ray: p(t) = r + t*q
    !
    q(1:3) = xx2 - xx1
    r(1:3) = xx1
    imax   = maxloc(abs(q))
    !
    ! We have 3 equations with 3 unknowns
    !
    ! qx * t = u*v * ax + u * bx + v * cx + dx - rx
    ! qy * t = u*v * ay + u * by + v * cy + dy - ry
    ! qz * t = u*v * az + u * bz + v * cz + dz - rz
    !
    ! Eliminate t from the equation where q is max
    !
    if( imax(1) == 3 ) then
       !
       ! Eliminate z-equation
       !
       A1 = a(1)*q(3) - a(3)*q(1)
       B1 = b(1)*q(3) - b(3)*q(1)
       C1 = c(1)*q(3) - c(3)*q(1)
       D1 = ( d(1)-r(1) )*q(3) - ( d(3)-r(3) )*q(1)

       A2 = a(2)*q(3) - a(3)*q(2)
       B2 = b(2)*q(3) - b(3)*q(2)
       C2 = c(2)*q(3) - c(3)*q(2)
       D2 = ( d(2)-r(2) )*q(3) - ( d(3)-r(3) )*q(2)

    else if( imax(1) == 2 ) then

       A1 = a(1)*q(2) - a(2)*q(1)
       B1 = b(1)*q(2) - b(2)*q(1)
       C1 = c(1)*q(2) - c(2)*q(1)
       D1 = ( d(1)-r(1) )*q(2) - ( d(2)-r(2) )*q(1)

       A2 = a(3)*q(2) - a(2)*q(3)
       B2 = b(3)*q(2) - b(2)*q(3)
       C2 = c(3)*q(2) - c(2)*q(3)
       D2 = ( d(3)-r(3) )*q(2) - ( d(2)-r(2) )*q(3)

    else if( imax(1) == 1 ) then

       A1 = a(2)*q(1) - a(1)*q(2)
       B1 = b(2)*q(1) - b(1)*q(2)
       C1 = c(2)*q(1) - c(1)*q(2)
       D1 = ( d(2)-r(2) )*q(1) - ( d(1)-r(1) )*q(2)

       A2 = a(3)*q(1) - a(1)*q(3)
       B2 = b(3)*q(1) - b(1)*q(3)
       C2 = c(3)*q(1) - c(1)*q(3)
       D2 = ( d(3)-r(3) )*q(1) - ( d(1)-r(1) )*q(3)

    else

       write(3000+kfl_paral,*) '###', xx1, '###'
       write(3000+kfl_paral,*) '###', xx2, '###'
       write(3000+kfl_paral,*) '<<<', toler, '>>>'
       flush (3000+kfl_paral)
       call runend(" elmgeo_intersection_segment_QUA04(): norm(segment) < toler fort.(3000+kfl_paral)")

    end if
    !
    ! Normalize coefficients
    !
    A1 = A1 * xnorm
    B1 = B1 * xnorm
    C1 = C1 * xnorm
    D1 = D1 * xnorm

    A2 = A2 * xnorm
    B2 = B2 * xnorm
    C2 = C2 * xnorm
    D2 = D2 * xnorm
    !
    ! Quadratic equation for v: aa*v^2 + bb*v + cc = 0
    !
    aa     = A2*C1 - A1*C2
    bb     = A2*D1 - A1*D2 + B2*C1 - B1*C2
    cc     = B2*D1 - B1*D2
    delta  = bb*bb - 4.0_rp*aa*cc

    t1     = -10.0_rp
    t2     = -10.0_rp
    patch1 = .false.
    patch2 = .false.
    ifoun  = 0
    jstat  = 0
    inedge = .false.

    if( abs(aa) <= epsilgeo_alg ) then
       !
       ! Linear equation
       !
       if( abs(bb) <= epsilgeo_alg) then
          patch1 = .false.
       else
          !
          ! One real solution:
          !
          !      -c
          ! v1 = ---
          !       b
          !
          v1 = -cc/bb
          call elmgeo_solve_patch(A1,B1,C1,D1,A2,B2,C2,D2,a,b,c,d,q,r,v1,u1,t1,p1,patch1,toler)
          if (patch1 .and. (abs(v1)<toler2 .or. abs(v1-1.0_rp)<toler2 .or. abs(u1)<toler2 .or. abs(u1-1.0_rp)<toler2)) inedge = .true.
       end if


    else if( delta >= 0.0_rp ) then

       if( delta <= epsilgeo_alg ) then
          !
          ! One real solution:
          !
          !      -b
          ! v1 = --
          !      2a
          !
          v1 = -bb/(2.0_rp*aa)
          call elmgeo_solve_patch(A1,B1,C1,D1,A2,B2,C2,D2,a,b,c,d,q,r,v1,u1,t1,p1,patch1,toler)
          if (patch1 .and. (abs(v1)<toler2 .or. abs(v1-1.0_rp)<toler2 .or. abs(u1)<toler2 .or. abs(u1-1.0_rp)<toler2 )) inedge = .true.
       else
          !
          ! Two real solutions
          !
          !if( 4.0_rp*aa*cc < 1.0e-6_rp * bb**2 ) then
          !   !
          !   ! To avoid errors when 4*a*c << b^2, Taylor Series: (1-x)^{1/2} = 1 - 1/2*x - 1/8*x^2
          !   !
          !   v1 = 0.5_rp*(-bb+abs(bb))/aa - cc/bb * ( 1.0_rp + aa*cc/(bb*bb) )
          !   v2 = 0.5_rp*(-bb-abs(bb))/aa - cc/bb * ( 1.0_rp + aa*cc/(bb*bb) )
          !   call elmgeo_solve_patch(A1,B1,C1,D1,A2,B2,C2,D2,a,b,c,d,q,r,v1,u1,t1,p1,patch1,toler)
          !   if (patch1 .and. (abs(v1)<toler2 .or. abs(v1-1.0_rp)<toler2 .or. abs(u1)<toler2 .or. abs(u1-1.0_rp)<toler2)) inedge = .true.
          !   call elmgeo_solve_patch(A1,B1,C1,D1,A2,B2,C2,D2,a,b,c,d,q,r,v2,u2,t2,p2,patch2,toler)
          !   if ( patch2 .and. (abs(v2)<toler2 .or. abs(v2-1.0_rp)<toler2 .or. abs(u2)<toler2 .or. abs(u2-1.0_rp)<toler2)) inedge = .true.
          !else
          !
          !      -b+sqrt(delta)       -b-sqrt(delta)
          ! v1 = --------------, v2 = --------------. For robustness, we rather compute
          !            2a                    2a
          !
          !      -b-sign(b)*sqrt(delta)        c
          ! v1 = ----------------------, v2 = ----
          !               2a                  a*v1
          !
          delta = sqrt(delta)
          if( bb >= 0.0_rp ) then
             v1 = (-bb-delta) / (2.0_rp*aa)
             if( v1 /= 0.0_rp ) then
                v2 = cc / (aa*v1)
             else
                v2 = (-bb+delta) / (2.0_rp*aa)
             end if
          else
             v1 = (-bb+delta) / (2.0_rp*aa)
             if( v1 /= 0.0_rp ) then
                v2 = cc / (aa*v1)
             else
                v2 = (-bb-delta) / (2.0_rp*aa)
             end if
          end if
          !v1    = (-bb+delta) / (2.0_rp*aa)
          !v2    = (-bb-delta) / (2.0_rp*aa)
          call elmgeo_solve_patch(A1,B1,C1,D1,A2,B2,C2,D2,a,b,c,d,q,r,v1,u1,t1,p1,patch1,toler)
          if ( patch1 .and. (abs(v1)<toler2 .or. abs(v1-1.0_rp)<toler2 .or. abs(u1)<toler2 .or. abs(u1-1.0_rp)<toler2)) inedge = .true.
          call elmgeo_solve_patch(A1,B1,C1,D1,A2,B2,C2,D2,a,b,c,d,q,r,v2,u2,t2,p2,patch2,toler)
          if ( patch2 .and. (abs(v2)<toler2 .or. abs(v2-1.0_rp)<toler2 .or. abs(u2)<toler2 .or. abs(u2-1.0_rp)<toler2)) inedge = .true.
          !end if

       end if
    end if

    if( .not. patch1 .and. patch2 ) then
       !
       ! Only one intersection point
       !
       ifoun        = 1
       intersection = p2
       if( abs(t2) <= toler ) then
          jstat = 1
       else if( abs(t2-1.0_rp) <= toler ) then
          jstat = 2
       end if

    else if( .not. patch2 .and. patch1 ) then
       !
       ! Only one intersection point
       !
       ifoun        = 1
       intersection = p1
       if( abs(t1) <= toler ) then
          jstat = 1
       else if( abs(t1-1.0_rp) <= toler ) then
          jstat = 2
       end if

    else if( patch1 .and. patch2 ) then
       !
       ! Two intersection points
       !
       if( t1 < t2 ) then
          ifoun        = 2
          intersection = p1
       else
          ifoun        = 2
          intersection = p2
       end if

       if(    ( abs(t2) <= toler .and. abs(t1-1.0_rp) <= toler ) .or. &
            & ( abs(t1) <= toler .and. abs(t2-1.0_rp) <= toler ) ) then
          jstat = 3
       end if

    end if
    if(inedge .and. jstat/=1 .and. jstat/= 2 .and. jstat /=3) jstat = 4
    if( present(istat) ) istat = jstat

  end subroutine elmgeo_intersection_segment_QUA04

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   Patch test
  !> @details Patch test
  !>
  !-----------------------------------------------------------------------

  subroutine elmgeo_solve_patch(A1,B1,C1,D1,A2,B2,C2,D2,a,b,c,d,q,r,v,u,t,p,patch,toler)
    real(rp),    intent(in)  :: A1,B1,C1,D1
    real(rp),    intent(in)  :: A2,B2,C2,D2
    real(rp),    intent(in)  :: a(3),b(3),c(3),d(3)
    real(rp),    intent(in)  :: q(3)
    real(rp),    intent(in)  :: r(3)
    real(rp),    intent(in)  :: v
    real(rp),    intent(in)  :: toler
    real(rp),    intent(out) :: u
    real(rp),    intent(out) :: t
    real(rp),    intent(out) :: p(3)
    logical(lg), intent(out) :: patch
    real(rp)                 :: aa,bb
    real(rp)                 :: toler2

    toler2 = toler
    patch  = .false.
    u      = 1000.0_rp

    if( v >= -toler2 .and. v <= 1.0_rp+toler2 ) then
       aa = v*A2 + B2
       bb = v*(A2 - A1) + B2 - B1
       if( abs(bb) >= abs(aa) ) then
          u = (v*(C1-C2)+D1-D2) / bb
       else
          u = (-v*C2-D2) / aa
       end if
       if( u >= -toler2 .and. u <= 1.0_rp+toler2 ) then
          p = u*v*a + u*b + v*c + d
          if( abs(q(1)) >= abs(q(2)) .and. abs(q(1)) >= abs(q(3)) ) then
             t = (p(1)-r(1)) / q(1)
          else if( abs(q(2)) > abs(q(3)) ) then
             t = (p(2)-r(2)) / q(2)
          else
             t = (p(3)-r(3)) / q(3)
          end if
          if( t >= -toler2 .and. t <= 1.0_rp+toler2 ) patch = .true.
       end if
    end if

  end subroutine elmgeo_solve_patch

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   Area of a face
  !> @details Area of a face
  !>
  !-----------------------------------------------------------------------

  subroutine elmgeo_face_area(ndime,pflty,bocod,face_area,normal)

    integer(ip), intent(in)            :: ndime                   !< Problem dimension
    integer(ip), intent(in)            :: pflty                   !< Face type
    real(rp),    intent(in)            :: bocod(ndime,*)          !< Face node coordinates
    real(rp),    intent(out)           :: face_area
    real(rp),    intent(out), optional :: normal(ndime)
    integer(ip)                        :: pelty
    real(rp)                           :: vect1(3),vect2(3),vect3(3)

    pelty = abs(pflty)

    if(       pelty == BAR02 ) then
       !
       ! BAR02
       !
       vect1(1:ndime) = bocod(1:ndime,2) - bocod(1:ndime,1)
       face_area      = sqrt(dot_product(vect1(1:ndime),vect1(1:ndime)))

       if( present(normal) ) then
          if(face_area /= 0_rp) then
             normal(1) = ( bocod(2,1)-bocod(2,2) )  / face_area
             normal(2) = ( bocod(1,2)-bocod(1,1) )  / face_area
          else
             call runend("elmgeo_face_area: face area is zero, normal not defined ")
          endif
       end if

    else if ( pelty == TRI03 ) then
       !
       ! TRI03
       !
       vect1(1:3) = bocod(1:3,2) - bocod(1:3,1)
       vect2(1:3) = bocod(1:3,3) - bocod(1:3,1)
       vect3(1)   = vect1(2)*vect2(3) - vect1(3)*vect2(2)
       vect3(2)   = vect1(3)*vect2(1) - vect1(1)*vect2(3)
       vect3(3)   = vect1(1)*vect2(2) - vect1(2)*vect2(1)
       face_area  = 0.5_rp * sqrt(dot_product(vect3(1:3),vect3(1:3)))
       if( present(normal) ) then
          if(face_area /= 0_rp) then
             normal(1:3) = vect3(1:3) / face_area
          else
             call runend("elmgeo_face_area: face area is zero, normal not defined ")
          endif
       end if

    else

       call runend('ELMGEO_FACE_AREA: TU TAMBIEN LO PUEDES PROGRAMAR!')

    end if

  end subroutine elmgeo_face_area
  
  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   Length of an element
  !> @details Length of an element
  !>
  !-----------------------------------------------------------------------
  
  subroutine elmgeo_element_distance(ndime,pflty,elcod,element_length)

    integer(ip), intent(in)            :: ndime                   !< Problem dimension
    integer(ip), intent(in)            :: pflty                   !< Face type
    real(rp),    intent(in)            :: elcod(ndime,*)          !< Face node coordinates
    real(rp),    intent(out)           :: element_length
    integer(ip)                        :: pelty
    real(rp)                           :: element_volume 

    pelty = abs(pflty)
    call elmgeo_element_volume(ndime,pflty,elcod,element_volume)

    if(      pelty >= element_num_ini(1) .and. pelty <= element_num_end(1) ) then
       element_length = element_volume
    else if( pelty >= element_num_ini(2) .and. pelty <= element_num_end(2) ) then
       element_length = sqrt(element_volume)
    else 
       element_length = element_volume ** (1.0_rp/3.0_rp)
    end if
 
  end subroutine elmgeo_element_distance
  
  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   Volume of an element
  !> @details Volume of an element
  !>
  !-----------------------------------------------------------------------
  
  pure subroutine elmgeo_BAR02_volume(ndime,v1,v2,element_volume)
    
    integer(ip), intent(in)            :: ndime              !< Problem dimension
    real(rp),    intent(in)            :: v1(ndime)          !< Face node coordinates
    real(rp),    intent(in)            :: v2(ndime)          !< Face node coordinates
    real(rp),    intent(out)           :: element_volume

    element_volume = sqrt(dot_product(v1(1:ndime)-v2(1:ndime),v1(1:ndime)-v2(1:ndime)))
    
  end subroutine elmgeo_BAR02_volume
  
  pure subroutine elmgeo_TRI03_volume(ndime,v1,v2,v3,element_volume)
    
    integer(ip), intent(in)            :: ndime              !< Problem dimension
    real(rp),    intent(in)            :: v1(ndime)          !< Face node coordinates
    real(rp),    intent(in)            :: v2(ndime)          !< Face node coordinates
    real(rp),    intent(in)            :: v3(ndime)          !< Face node coordinates
    real(rp),    intent(out)           :: element_volume
    real(rp)                           :: vect1(3),vect2(3)
    real(rp)                           :: vect3(3)

    vect1(1) = v2(1) - v1(1)
    vect1(2) = v2(2) - v1(2)
    vect1(3) = v2(3) - v1(3)

    vect2(1) = v3(1) - v1(1)
    vect2(2) = v3(2) - v1(2)
    vect2(3) = v3(3) - v1(3)
    
    vect3(1) =   vect1(2) * vect2(3) - vect1(3) * vect2(2) ! vect1 x vect2
    vect3(2) = - vect1(1) * vect2(3) + vect1(3) * vect2(1)
    vect3(3) =   vect1(1) * vect2(2) - vect1(2) * vect2(1)
    
    element_volume = sqrt(dot_product(v3(1:ndime),v3(1:ndime)))
    
  end subroutine elmgeo_TRI03_volume
  
  pure subroutine elmgeo_TET04_volume(ndime,v1,v2,v3,v4,element_volume)

    integer(ip), intent(in)            :: ndime              !< Problem dimension
    real(rp),    intent(in)            :: v1(ndime)          !< Face node coordinates
    real(rp),    intent(in)            :: v2(ndime)          !< Face node coordinates
    real(rp),    intent(in)            :: v3(ndime)          !< Face node coordinates
    real(rp),    intent(in)            :: v4(ndime)          !< Face node coordinates
    real(rp),    intent(out)           :: element_volume
    real(rp)                           :: vect1(3),vect2(3)
    real(rp)                           :: vect3(3),vect4(3)

    ! vol = ([(v4-v2) x (v4-v3)] · [v4-v1] )/6.0 where v_i are the vertices

    vect1(1) = v4(1) - v2(1)
    vect1(2) = v4(2) - v2(2)
    vect1(3) = v4(3) - v2(3)

    vect2(1) = v4(1) - v3(1)
    vect2(2) = v4(2) - v3(2)
    vect2(3) = v4(3) - v3(3)

    vect3(1) = v4(1) - v1(1)
    vect3(2) = v4(2) - v1(2)
    vect3(3) = v4(3) - v1(3)

    vect4(1) =   vect1(2) * vect2(3) - vect1(3) * vect2(2) ! vect1 x vect2
    vect4(2) = - vect1(1) * vect2(3) + vect1(3) * vect2(1)
    vect4(3) =   vect1(1) * vect2(2) - vect1(2) * vect2(1)

    element_volume = (1.0_rp/6.0_rp)*abs(dot_product(vect4(1:3), vect3(1:3)))

  end subroutine elmgeo_TET04_volume

  pure subroutine elmgeo_PYR05_volume(ndime,v1,v2,v3,v4,v5,element_volume)

    integer(ip), intent(in)            :: ndime              !< Problem dimension
    real(rp),    intent(in)            :: v1(ndime)          !< Face node coordinates
    real(rp),    intent(in)            :: v2(ndime)          !< Face node coordinates
    real(rp),    intent(in)            :: v3(ndime)          !< Face node coordinates
    real(rp),    intent(in)            :: v4(ndime)          !< Face node coordinates
    real(rp),    intent(in)            :: v5(ndime)          !< Face node coordinates
    real(rp),    intent(out)           :: element_volume
    real(rp)                           :: volume1, volume2

    ! A PYR05 is two merged TET04...mmm... what about bilinear face

    call elmgeo_TET04_volume(ndime,v1,v2,v3,v5,volume1)
    call elmgeo_TET04_volume(ndime,v1,v4,v3,v5,volume2)

    element_volume = volume1 + volume2

  end subroutine elmgeo_PYR05_volume

  subroutine elmgeo_element_volume(ndime,pflty,elcod,element_volume,element_center)

    integer(ip), intent(in)            :: ndime                   !< Problem dimension
    integer(ip), intent(in)            :: pflty                   !< Face type
    real(rp),    intent(in)            :: elcod(ndime,*)          !< Face node coordinates
    real(rp),    intent(out)           :: element_volume
    real(rp),    intent(out), optional :: element_center(ndime)
    integer(ip)                        :: pelty,idime
    real(rp)                           :: xfact

    pelty = abs(pflty)

    if( pelty == BAR02 ) then
       !
       ! BAR02
       !
       call elmgeo_BAR02_volume(ndime,elcod(1:ndime,1),elcod(1:ndime,2),element_volume)          
       
    else if( pelty == TRI03 ) then
       !
       ! TRI03
       !
       if( ndime == 2 ) then
          element_volume =  0.5_rp*abs(  (elcod(1,2)-elcod(1,1))*(elcod(2,3)-elcod(2,2))    &
               &                        -(elcod(2,2)-elcod(2,1))*(elcod(1,3)-elcod(1,2)) )
       else
          call elmgeo_TRI03_volume(ndime,elcod(1:ndime,1),elcod(1:ndime,2),elcod(1:ndime,3),element_volume)          
       end if
       
    else if( pelty == QUA04 ) then
       !
       ! QUA04
       !
       element_volume =  0.5_rp*abs(  (elcod(1,2)-elcod(1,1))*(elcod(2,3)-elcod(2,2))    &
            &                        -(elcod(2,2)-elcod(2,1))*(elcod(1,3)-elcod(1,2)) )  &
            &          + 0.5_rp*abs(  (elcod(1,1)-elcod(1,3))*(elcod(2,4)-elcod(2,1))    &
            &                        -(elcod(2,1)-elcod(2,3))*(elcod(1,4)-elcod(1,1)) )

    else if ( pelty == TET04 ) then

       call elmgeo_TET04_volume(ndime,elcod(1:ndime,1),elcod(1:ndime,2),elcod(1:ndime,3),elcod(1:ndime,4),element_volume)

    else if ( pelty == PYR05 ) then

       call elmgeo_PYR05_volume(ndime,elcod(1:ndime,1),elcod(1:ndime,2),elcod(1:ndime,3),elcod(1:ndime,4),elcod(1:ndime,5),element_volume)

    else

       call runend('ELMGEO_ELEMENT_VOLUME: VAGO, TU TAMBIEN LO PUEDES PROGRAMAR!')

    end if

    if( present(element_center) ) then ! El centroide no tendria que ser el centro de massa? Es equivalente??
       xfact = 1.0_rp / real(element_type(pelty) % number_nodes,rp)
       do idime = 1,ndime
          element_center(idime) = sum(elcod(idime,1:element_type(pelty) % number_nodes)) * xfact
       end do
    end if

  end subroutine elmgeo_element_volume

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   Natural coordinates
  !> @details Evaluates shape functions, derivatives and Hessian
  !>          of a bubble
  !
  !-----------------------------------------------------------------------

  subroutine elmgeo_bubble(&
       itask,ndime,mnode,pnode,pgaus,elcod,gpsha,deriv,gpcar,&
       gpsha_bub,gpcar_bub,elfun)

    integer(ip), intent(in)           :: itask
    integer(ip), intent(in)           :: ndime
    integer(ip), intent(in)           :: mnode
    integer(ip), intent(in)           :: pnode
    integer(ip), intent(in)           :: pgaus
    real(rp),    intent(in)           :: elcod(ndime,pnode)
    real(rp),    intent(in)           :: gpsha(pnode,pgaus)
    real(rp),    intent(in)           :: deriv(ndime,pnode,pgaus)
    real(rp),    intent(in)           :: gpcar(ndime,mnode,pgaus)
    real(rp),    intent(out)          :: gpsha_bub(1,pgaus)
    real(rp),    intent(out)          :: gpcar_bub(ndime,pgaus)
    real(rp),    intent(in), optional :: elfun(pnode)
    real(rp)                          :: deriv_bub(ndime,1,pgaus)
    real(rp)                          :: heslo_bub(3*ndime-3,1,pgaus)
    integer(ip)                       :: igaus,idime,inode
    real(rp)                          :: xjaci(ndime,ndime)
    real(rp)                          :: posgp(3)
    real(rp)                          :: gpdet,xplus,xminu,fmini,fmaxi
    real(rp)                          :: f(pnode),xsign

    if(      itask == FLAT_BUBBLE ) then
       !
       ! Constant bubble
       !
       gpsha_bub = 1.0_rp
       gpcar_bub = 0.0_rp

    else if( itask == QUADRATIC_BUBBLE ) then
       !
       ! Quadratic bubble
       !
       do igaus = 1,pgaus
          do idime = 1,ndime
             posgp(idime) = dot_product(gpsha(1:pnode,igaus),elcod(idime,1:pnode))
          end do
          call elmgeo_shapf_deriv_heslo_bubble(&
               ndime,pnode,posgp,gpsha_bub(:,igaus),deriv_bub(:,:,igaus),heslo_bub(:,:,igaus))
       end do

       do igaus = 1,pgaus
          call elmgeo_jacobian_matrix(&
               ndime,pnode,elcod,deriv(1,1,igaus),gpdet,xjaci)
          do idime = 1,ndime
             gpcar_bub(idime,igaus) = dot_product(xjaci(1:ndime,idime),deriv_bub(1:ndime,1,igaus))
          end do
       end do

    else if( itask == FREE_SURFACE_BUBBLE ) then
       !
       ! Free surface bubble
       ! Ne = sum_i Ni | phi_i | - | sum_i N_i * phi_i |
       !    = xplus - xminu
       ! if xplus > 0
       !   dNe/dx =   sum_i dN_i/dx [ phi_i - abs(phi_i) ]
       ! else
       !   dNe/dx = - sum_i dN_i/dx [ phi_i + abs(phi_i) ]
       ! endif
       !
       !#
       !# EXAMPLE FOR TRI03
       !#
       !set xrange[0:1] ; set yrange[0:1] ; set pm3d ; set hidden3d
       !set samples    50 ; set isosamples 50
       !phi1       = -0.3
       !phi2       =  0.4
       !phi3       =  0.5
       !n1(x,y)    =  1.0-x-y
       !n2(x,y)    =  x
       !n3(x,y)    =  y
       !dn1(x,y)   = -1.0
       !dn2(x,y)   =  1.0
       !dn3(x,y)   =  0.0
       !xplus(x,y) =  n1(x,y)*phi1 + n2(x,y)*phi2 + n3(x,y)*phi3
       !g(x,y)     =  n1(x,y)*abs(phi1) + n2(x,y)*abs(phi2) + n3(x,y)*abs(phi3) - abs(xplus(x,y))
       !dgp(x,y)   =  dn1(x,y)*(abs(phi1)-phi1) + dn2(x,y)*(abs(phi2)-phi2) + dn3(x,y)*(abs(phi3)-phi3)
       !dgm(x,y)   =  dn1(x,y)*(abs(phi1)+phi1) + dn2(x,y)*(abs(phi2)+phi2) + dn3(x,y)*(abs(phi3)+phi3)
       !dg(x,y)    =  (xplus(x,y)>0) ? dgp(x,y) : dgm(x,y)
       !df(x,y)    =  (x+y<1) ? dg(x,y) : 0.0
       !f(x,y)     =  (x+y<1) ?  g(x,y) : 0.0
       !splot f(x,y)
       !
       gpcar_bub = 0.0_rp
       fmini     = minval(elfun)
       fmaxi     = maxval(elfun)
!!!!!!!!!!f         = 2.0_rp*( elfun - fmini ) / ( fmaxi - fmini ) -1.0_rp  !!!OJOOO
       f         = elfun
       do igaus = 1,pgaus
          xplus              = dot_product(gpsha(1:pnode,igaus),f(1:pnode))
          xminu              = dot_product(gpsha(1:pnode,igaus),abs(f(1:pnode)))
          gpsha_bub(1,igaus) = xminu - abs(xplus)
          xsign              = sign(1.0_rp,xplus)
          do inode = 1,pnode
             gpcar_bub(1:ndime,igaus) = gpcar_bub(1:ndime,igaus) &
                  + gpcar(1:ndime,inode,igaus) * ( abs(f(inode)) - xsign * f(inode) )
          end do
          !if(gpsha_bub(igaus)>100.0_rp)then
          !   print*,'ojoooooooo'
          !   print*,'gpsha',gpsha(1:pnode,igaus)
          !   print*,'level',f(1:pnode)
          !   stop
          !end if
       end do
    else

       call runend('ELMGEO_BUBBLE: BUBBLE NOT CODE')

    end if
  end subroutine elmgeo_bubble

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   L2 norm of an array
  !> @details L2 norm of an array
  !
  !------------------------------------------------------------------------

  function elmgeo_norm2(ndime,array)
    integer(ip), intent(in) :: ndime
    real(rp),    intent(in) :: array(ndime)
    real(rp)                :: elmgeo_norm2

    elmgeo_norm2 = sqrt(dot_product(array,array))

  end function elmgeo_norm2

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    15/09/2016
  !> @brief   Computes the Jacobian, Jacobian determinant and
  !>          Cartesian derivatives of shape function of an element
  !> @details The jacobian is
  !>                             _           _
  !>                            | dx/ds dx/dt |                t
  !>          Jacobian: XJACM = |             | = ELCOD * DERIV
  !>                            |_dy/ds dy/dt_|
  !>          with
  !>                   _        _             _                    _
  !>                  | x1 x2 x3 |           | dN1/ds dN2/ds dN3/ds |
  !>          ELCOD = |          |,  DERIV = |                      |
  !>                  |_y1 y2 y3_|           |_dN1/dt dN2/dt dN3/dt_|
  !>
  !>          => Jacobian determinant: GPDET = det(XJACM)
  !>
  !>          P1 Element in 2D
  !>          ----------------
  !>
  !>          gpdet=(-x1+x2)*(-y1+y3)-(-y1+y2)*(-x1+x3)
  !>                         _                       _
  !>                        | -y3+y2  -y1+y3   y1-y2  |
  !>          gpcar=1/gpdet |                         |
  !>                        |_ x3-x2   x1-x3  -x1+x2 _|
  !>
  !>          P1 Element in 3D
  !>          ----------------
  !>                 _                     _
  !>                |  x2-x1  x3-x1  x4-x1  |
  !>          xjacm=|  y2-y1  y3-y1  y4-y1  |
  !>                |_ z2-z1  z3-z1  z4-z1 _|
  !>
  !>                     -1
  !>          xjaci=xjacm
  !>                 _                             _     _          _
  !>                |  dN1/ds dN2/ds dN3/ds dN4/ds  |   |  -1 1 0 0  |
  !>          deriv=|  dN1/dt dN2/dt dN3/dt_dN4/dt  | = |  -1 0 1 0  |
  !>                |_ dN1/dz dN2/dz dN3/dz_dN4/dz _|   |_ -1 0 0 1 _|
  !>
  !>                     t
  !>          cartd=xjaci *deriv
  !>
  !-----------------------------------------------------------------------

  subroutine elmgeo_cartesian_derivatives_jacobian_vector(ndime,mnode,pnode,pgaus,elcod,deriv,xjaci,gpcar,gpdet)

    integer(ip),           intent(in)    :: ndime
    integer(ip),           intent(in)    :: mnode
    integer(ip),           intent(in)    :: pnode
    integer(ip),           intent(in)    :: pgaus
    real(rp),              intent(in)    :: elcod(VECTOR_SIZE,ndime,pnode)
    real(rp),              intent(in)    :: deriv(ndime,pnode,pgaus)
    real(rp),              intent(inout) :: xjaci(VECTOR_SIZE,ndime,ndime,pgaus)
    real(rp),              intent(out)   :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)
    real(rp),              intent(out)   :: gpdet(VECTOR_SIZE,pgaus)
    integer(ip)                          :: j,k,idime,inode,igaus,jdime,ivect
    real(rp)                             :: xjacm(VECTOR_SIZE,ndime,ndime)
    real(rp)                             :: t1(VECTOR_SIZE)
    real(rp)                             :: t2(VECTOR_SIZE)
    real(rp)                             :: t3(VECTOR_SIZE)
    real(rp)                             :: denom(VECTOR_SIZE)
    
    if( ndime == 2 .and. pnode == 3 ) then
       !
       ! 2D P1 element
       !
       gpdet(1:VECTOR_SIZE,1)     = (    -elcod(1:VECTOR_SIZE,1,1) + elcod(1:VECTOR_SIZE,1,2)) * (-elcod(1:VECTOR_SIZE,2,1) + elcod(1:VECTOR_SIZE,2,3) ) &
            &                         -( -elcod(1:VECTOR_SIZE,2,1) + elcod(1:VECTOR_SIZE,2,2)) * (-elcod(1:VECTOR_SIZE,1,1) + elcod(1:VECTOR_SIZE,1,3) )

       denom(1:VECTOR_SIZE)       = 1.0_rp / (sign(1.0_rp,gpdet(1:VECTOR_SIZE,1))*max(abs(gpdet(1:VECTOR_SIZE,1)),epsilgeo_div))

       gpcar(1:VECTOR_SIZE,1,1,1) = ( -elcod(1:VECTOR_SIZE,2,3) + elcod(1:VECTOR_SIZE,2,2) ) * denom(1:VECTOR_SIZE)
       gpcar(1:VECTOR_SIZE,1,2,1) = ( -elcod(1:VECTOR_SIZE,2,1) + elcod(1:VECTOR_SIZE,2,3) ) * denom(1:VECTOR_SIZE)
       gpcar(1:VECTOR_SIZE,1,3,1) = (  elcod(1:VECTOR_SIZE,2,1) - elcod(1:VECTOR_SIZE,2,2) ) * denom(1:VECTOR_SIZE)
       gpcar(1:VECTOR_SIZE,2,1,1) = (  elcod(1:VECTOR_SIZE,1,3) - elcod(1:VECTOR_SIZE,1,2) ) * denom(1:VECTOR_SIZE)
       gpcar(1:VECTOR_SIZE,2,2,1) = (  elcod(1:VECTOR_SIZE,1,1) - elcod(1:VECTOR_SIZE,1,3) ) * denom(1:VECTOR_SIZE)
       gpcar(1:VECTOR_SIZE,2,3,1) = ( -elcod(1:VECTOR_SIZE,1,1) + elcod(1:VECTOR_SIZE,1,2) ) * denom(1:VECTOR_SIZE)

       !TODO: xjaci may be uninitialized
       do igaus = 2,pgaus
          gpdet(1:VECTOR_SIZE,igaus) = gpdet(1:VECTOR_SIZE,1)
          do jdime = 1,2
             do idime = 1,2
                xjaci(1:VECTOR_SIZE,idime,jdime,igaus) = xjaci(1:VECTOR_SIZE,idime,jdime,1)
             end do
          end do
          do inode = 1,3
             do idime = 1,2
                gpcar(1:VECTOR_SIZE,idime,inode,igaus) = gpcar(1:VECTOR_SIZE,idime,inode,1)
             end do
          end do
       end do
      
    else if( ndime == 3 .and. pnode == 4 ) then
       !
       ! 3D P1 element
       !
       do ivect = 1,VECTOR_SIZE
          gpcar(ivect,1,1,1) =  elcod(ivect,1,2)   - elcod(ivect,1,1)
          gpcar(ivect,1,2,1) =  elcod(ivect,1,3)   - elcod(ivect,1,1)
          gpcar(ivect,1,3,1) =  elcod(ivect,1,4)   - elcod(ivect,1,1)
          gpcar(ivect,2,1,1) =  elcod(ivect,2,2)   - elcod(ivect,2,1)
          gpcar(ivect,2,2,1) =  elcod(ivect,2,3)   - elcod(ivect,2,1)
          gpcar(ivect,2,3,1) =  elcod(ivect,2,4)   - elcod(ivect,2,1)
          gpcar(ivect,3,1,1) =  elcod(ivect,3,2)   - elcod(ivect,3,1)
          gpcar(ivect,3,2,1) =  elcod(ivect,3,3)   - elcod(ivect,3,1)
          gpcar(ivect,3,3,1) =  elcod(ivect,3,4)   - elcod(ivect,3,1)
          t1(ivect)          =  gpcar(ivect,2,2,1) * gpcar(ivect,3,3,1) - gpcar(ivect,3,2,1) * gpcar(ivect,2,3,1)
          t2(ivect)          = -gpcar(ivect,2,1,1) * gpcar(ivect,3,3,1) + gpcar(ivect,3,1,1) * gpcar(ivect,2,3,1)
          t3(ivect)          =  gpcar(ivect,2,1,1) * gpcar(ivect,3,2,1) - gpcar(ivect,3,1,1) * gpcar(ivect,2,2,1)
          gpdet(ivect,1)     =  gpcar(ivect,1,1,1) * t1(ivect) + gpcar(ivect,1,2,1) * t2(ivect) + gpcar(ivect,1,3,1) * t3(ivect)

          denom(ivect)       =  1.0_rp / (sign(1.0_rp,gpdet(ivect,1))*max(abs(gpdet(ivect,1)),epsilgeo_div))

          xjaci(ivect,1,1,1) =  t1(ivect) * denom(ivect)
          xjaci(ivect,2,1,1) =  t2(ivect) * denom(ivect)
          xjaci(ivect,3,1,1) =  t3(ivect) * denom(ivect)
          xjaci(ivect,2,2,1) = ( gpcar(ivect,1,1,1) * gpcar(ivect,3,3,1) - gpcar(ivect,3,1,1) * gpcar(ivect,1,3,1)) * denom(ivect)
          xjaci(ivect,3,2,1) = (-gpcar(ivect,1,1,1) * gpcar(ivect,3,2,1) + gpcar(ivect,1,2,1) * gpcar(ivect,3,1,1)) * denom(ivect)
          xjaci(ivect,3,3,1) = ( gpcar(ivect,1,1,1) * gpcar(ivect,2,2,1) - gpcar(ivect,2,1,1) * gpcar(ivect,1,2,1)) * denom(ivect)
          xjaci(ivect,1,2,1) = (-gpcar(ivect,1,2,1) * gpcar(ivect,3,3,1) + gpcar(ivect,3,2,1) * gpcar(ivect,1,3,1)) * denom(ivect)
          xjaci(ivect,1,3,1) = ( gpcar(ivect,1,2,1) * gpcar(ivect,2,3,1) - gpcar(ivect,2,2,1) * gpcar(ivect,1,3,1)) * denom(ivect)
          xjaci(ivect,2,3,1) = (-gpcar(ivect,1,1,1) * gpcar(ivect,2,3,1) + gpcar(ivect,2,1,1) * gpcar(ivect,1,3,1)) * denom(ivect)

          gpcar(ivect,1,1,1) = -xjaci(ivect,1,1,1) - xjaci(ivect,2,1,1) - xjaci(ivect,3,1,1)
          gpcar(ivect,1,2,1) =  xjaci(ivect,1,1,1)
          gpcar(ivect,1,3,1) =  xjaci(ivect,2,1,1)
          gpcar(ivect,1,4,1) =  xjaci(ivect,3,1,1)
          gpcar(ivect,2,1,1) = -xjaci(ivect,1,2,1) - xjaci(ivect,2,2,1) - xjaci(ivect,3,2,1)
          gpcar(ivect,2,2,1) =  xjaci(ivect,1,2,1)
          gpcar(ivect,2,3,1) =  xjaci(ivect,2,2,1)
          gpcar(ivect,2,4,1) =  xjaci(ivect,3,2,1)
          gpcar(ivect,3,1,1) = -xjaci(ivect,1,3,1) - xjaci(ivect,2,3,1) - xjaci(ivect,3,3,1)
          gpcar(ivect,3,2,1) =  xjaci(ivect,1,3,1)
          gpcar(ivect,3,3,1) =  xjaci(ivect,2,3,1)
          gpcar(ivect,3,4,1) =  xjaci(ivect,3,3,1)

          do igaus = 2,pgaus
             gpdet(ivect,igaus) = gpdet(ivect,1)
             do jdime = 1,3
                do idime = 1,3
                   xjaci(ivect,idime,jdime,igaus) = xjaci(ivect,idime,jdime,1)
                end do
             end do
             do inode = 1,4
                do idime = 1,3
                   gpcar(ivect,idime,inode,igaus) = gpcar(ivect,idime,inode,1)
                end do
             end do
          end do
       end do

    else if ( ndime == 1 ) then
       !
       ! 1D
       !
       do igaus = 1,pgaus
          xjacm(1:VECTOR_SIZE,1,1) = 0.0_rp
          do k = 1,pnode
             xjacm(1:VECTOR_SIZE,1,1) = xjacm(1:VECTOR_SIZE,1,1) + elcod(1:VECTOR_SIZE,1,k) * deriv(1,k,igaus)
          end do
          gpdet(1:VECTOR_SIZE,igaus) = xjacm(1:VECTOR_SIZE,1,1)
          xjaci(1:VECTOR_SIZE,1,1,igaus) = 1.0_rp / xjacm(1:VECTOR_SIZE,1,1)
          do j = 1,pnode
             gpcar(1:VECTOR_SIZE,1,j,igaus) = xjaci(1:VECTOR_SIZE,1,1,igaus) * deriv(1,j,igaus)
          end do
       end do
       
    else if ( ndime == 2 ) then
       !
       ! 2D
       !
       do igaus = 1,pgaus

          xjacm(1:VECTOR_SIZE,1,1) = 0.0_rp
          xjacm(1:VECTOR_SIZE,1,2) = 0.0_rp
          xjacm(1:VECTOR_SIZE,2,1) = 0.0_rp
          xjacm(1:VECTOR_SIZE,2,2) = 0.0_rp
          do k = 1,pnode
             xjacm(1:VECTOR_SIZE,1,1) = xjacm(1:VECTOR_SIZE,1,1) + elcod(1:VECTOR_SIZE,1,k) * deriv(1,k,igaus)
             xjacm(1:VECTOR_SIZE,1,2) = xjacm(1:VECTOR_SIZE,1,2) + elcod(1:VECTOR_SIZE,1,k) * deriv(2,k,igaus)
             xjacm(1:VECTOR_SIZE,2,1) = xjacm(1:VECTOR_SIZE,2,1) + elcod(1:VECTOR_SIZE,2,k) * deriv(1,k,igaus)
             xjacm(1:VECTOR_SIZE,2,2) = xjacm(1:VECTOR_SIZE,2,2) + elcod(1:VECTOR_SIZE,2,k) * deriv(2,k,igaus)
          end do

          gpdet(1:VECTOR_SIZE,igaus) = xjacm(1:VECTOR_SIZE,1,1) * xjacm(1:VECTOR_SIZE,2,2) - xjacm(1:VECTOR_SIZE,2,1) * xjacm(1:VECTOR_SIZE,1,2)

          denom(1:VECTOR_SIZE)       = 1.0_rp / (sign(1.0_rp,gpdet(1:VECTOR_SIZE,igaus))*max(abs(gpdet(1:VECTOR_SIZE,igaus)),epsilgeo_div))

          xjaci(1:VECTOR_SIZE,1,1,igaus) =  xjacm(1:VECTOR_SIZE,2,2) * denom(1:VECTOR_SIZE)
          xjaci(1:VECTOR_SIZE,2,2,igaus) =  xjacm(1:VECTOR_SIZE,1,1) * denom(1:VECTOR_SIZE)
          xjaci(1:VECTOR_SIZE,2,1,igaus) = -xjacm(1:VECTOR_SIZE,2,1) * denom(1:VECTOR_SIZE)
          xjaci(1:VECTOR_SIZE,1,2,igaus) = -xjacm(1:VECTOR_SIZE,1,2) * denom(1:VECTOR_SIZE)

          do j = 1, pnode
             gpcar(1:VECTOR_SIZE,1,j,igaus) =   xjaci(1:VECTOR_SIZE,1,1,igaus) * deriv(1,j,igaus) &
                  &                           + xjaci(1:VECTOR_SIZE,2,1,igaus) * deriv(2,j,igaus)

             gpcar(1:VECTOR_SIZE,2,j,igaus) =   xjaci(1:VECTOR_SIZE,1,2,igaus) * deriv(1,j,igaus) &
                  &                           + xjaci(1:VECTOR_SIZE,2,2,igaus) * deriv(2,j,igaus)
          end do
          
       end do

    else if ( ndime == 3 ) then
       !
       ! 3D
       !
       ! xjacm = elcod * deriv^t
       ! xjaci = xjacm^-1
       ! gpcar = xjaci^t * deriv
       !
       do igaus = 1,pgaus

          do ivect = 1,VECTOR_SIZE

             xjacm(ivect,1:3,1:3) = 0.0_rp
             do k = 1,pnode
                xjacm(ivect,1,1) = xjacm(ivect,1,1) + elcod(ivect,1,k) * deriv(1,k,igaus)
                xjacm(ivect,1,2) = xjacm(ivect,1,2) + elcod(ivect,1,k) * deriv(2,k,igaus)
                xjacm(ivect,1,3) = xjacm(ivect,1,3) + elcod(ivect,1,k) * deriv(3,k,igaus)
                xjacm(ivect,2,1) = xjacm(ivect,2,1) + elcod(ivect,2,k) * deriv(1,k,igaus)
                xjacm(ivect,2,2) = xjacm(ivect,2,2) + elcod(ivect,2,k) * deriv(2,k,igaus)
                xjacm(ivect,2,3) = xjacm(ivect,2,3) + elcod(ivect,2,k) * deriv(3,k,igaus)
                xjacm(ivect,3,1) = xjacm(ivect,3,1) + elcod(ivect,3,k) * deriv(1,k,igaus)
                xjacm(ivect,3,2) = xjacm(ivect,3,2) + elcod(ivect,3,k) * deriv(2,k,igaus)
                xjacm(ivect,3,3) = xjacm(ivect,3,3) + elcod(ivect,3,k) * deriv(3,k,igaus)
             end do

             t1(ivect)              =  xjacm(ivect,2,2) * xjacm(ivect,3,3) - xjacm(ivect,3,2) * xjacm(ivect,2,3)
             t2(ivect)              = -xjacm(ivect,2,1) * xjacm(ivect,3,3) + xjacm(ivect,3,1) * xjacm(ivect,2,3)
             t3(ivect)              =  xjacm(ivect,2,1) * xjacm(ivect,3,2) - xjacm(ivect,3,1) * xjacm(ivect,2,2)
             gpdet(ivect,igaus)     =  xjacm(ivect,1,1) * t1(ivect) + xjacm(ivect,1,2) * t2(ivect) + xjacm(ivect,1,3) * t3(ivect)

             denom(ivect)           = 1.0_rp / (sign(1.0_rp,gpdet(ivect,igaus))*max(abs(gpdet(ivect,igaus)),epsilgeo_div))
             xjaci(ivect,1,1,igaus) = t1(ivect) * denom(ivect)
             xjaci(ivect,2,1,igaus) = t2(ivect) * denom(ivect)
             xjaci(ivect,3,1,igaus) = t3(ivect) * denom(ivect)
             xjaci(ivect,2,2,igaus) = ( xjacm(ivect,1,1) * xjacm(ivect,3,3) - xjacm(ivect,3,1) * xjacm(ivect,1,3)) * denom(ivect)
             xjaci(ivect,3,2,igaus) = (-xjacm(ivect,1,1) * xjacm(ivect,3,2) + xjacm(ivect,1,2) * xjacm(ivect,3,1)) * denom(ivect)
             xjaci(ivect,3,3,igaus) = ( xjacm(ivect,1,1) * xjacm(ivect,2,2) - xjacm(ivect,2,1) * xjacm(ivect,1,2)) * denom(ivect)
             xjaci(ivect,1,2,igaus) = (-xjacm(ivect,1,2) * xjacm(ivect,3,3) + xjacm(ivect,3,2) * xjacm(ivect,1,3)) * denom(ivect)
             xjaci(ivect,1,3,igaus) = ( xjacm(ivect,1,2) * xjacm(ivect,2,3) - xjacm(ivect,2,2) * xjacm(ivect,1,3)) * denom(ivect)
             xjaci(ivect,2,3,igaus) = (-xjacm(ivect,1,1) * xjacm(ivect,2,3) + xjacm(ivect,2,1) * xjacm(ivect,1,3)) * denom(ivect)

             do j = 1, pnode
                gpcar(ivect,1,j,igaus) =   xjaci(ivect,1,1,igaus) * deriv(1,j,igaus) &
                     &                           + xjaci(ivect,2,1,igaus) * deriv(2,j,igaus) &
                     &                           + xjaci(ivect,3,1,igaus) * deriv(3,j,igaus)

                gpcar(ivect,2,j,igaus) =   xjaci(ivect,1,2,igaus) * deriv(1,j,igaus) &
                     &                           + xjaci(ivect,2,2,igaus) * deriv(2,j,igaus) &
                     &                           + xjaci(ivect,3,2,igaus) * deriv(3,j,igaus)

                gpcar(ivect,3,j,igaus) =   xjaci(ivect,1,3,igaus) * deriv(1,j,igaus) &
                     &                           + xjaci(ivect,2,3,igaus) * deriv(2,j,igaus) &
                     &                           + xjaci(ivect,3,3,igaus) * deriv(3,j,igaus)
             end do

          end do

       end do
    end if
 
  end subroutine elmgeo_cartesian_derivatives_jacobian_vector

  subroutine elmgeo_cartesian_derivatives_jacobian_scalar(ndime,mnode,pnode,pgaus,elcod,deriv,xjaci,gpcar,gpdet)

    integer(ip), intent(in)  :: ndime
    integer(ip), intent(in)  :: mnode
    integer(ip), intent(in)  :: pnode
    integer(ip), intent(in)  :: pgaus
    real(rp),    intent(in)  :: elcod(ndime,pnode)
    real(rp),    intent(in)  :: deriv(ndime,pnode,pgaus)
    real(rp),    intent(out) :: xjaci(ndime,ndime,pgaus)
    real(rp),    intent(out) :: gpcar(ndime,mnode,pgaus)
    real(rp),    intent(out) :: gpdet(pgaus)
    integer(ip)              :: j,k,idime,inode,igaus
    real(rp)                 :: xjacm(ndime,ndime)
    real(rp)                 :: t1
    real(rp)                 :: t2
    real(rp)                 :: t3
    real(rp)                 :: denom

    if( ndime == 2 .and. pnode == 3 ) then
       !
       ! 2D P1 element
       !
       gpdet(1)     = (    -elcod(1,1) + elcod(1,2)) * (-elcod(2,1) + elcod(2,3) ) &
            &                         -( -elcod(2,1) + elcod(2,2)) * (-elcod(1,1) + elcod(1,3) )
       denom       = 1.0_rp / (sign(1.0_rp,gpdet(1))*max(abs(gpdet(1)),epsilgeo_div))

       gpcar(1,1,1) = ( -elcod(2,3) + elcod(2,2) ) * denom
       gpcar(1,2,1) = ( -elcod(2,1) + elcod(2,3) ) * denom
       gpcar(1,3,1) = (  elcod(2,1) - elcod(2,2) ) * denom
       gpcar(2,1,1) = (  elcod(1,3) - elcod(1,2) ) * denom
       gpcar(2,2,1) = (  elcod(1,1) - elcod(1,3) ) * denom
       gpcar(2,3,1) = ( -elcod(1,1) + elcod(1,2) ) * denom

       xjaci = 0.0_rp

       do igaus = 2,pgaus
          gpdet(igaus) = gpdet(1)
          do inode = 1,3
             do idime = 1,2
                gpcar(idime,inode,igaus) = gpcar(idime,inode,1)
             end do
          end do
       end do

    else if( ndime == 3 .and. pnode == 4 ) then
       !
       ! 3D P1 element
       !
       gpcar(1,1,1) =  elcod(1,2)   - elcod(1,1)
       gpcar(1,2,1) =  elcod(1,3)   - elcod(1,1)
       gpcar(1,3,1) =  elcod(1,4)   - elcod(1,1)
       gpcar(2,1,1) =  elcod(2,2)   - elcod(2,1)
       gpcar(2,2,1) =  elcod(2,3)   - elcod(2,1)
       gpcar(2,3,1) =  elcod(2,4)   - elcod(2,1)
       gpcar(3,1,1) =  elcod(3,2)   - elcod(3,1)
       gpcar(3,2,1) =  elcod(3,3)   - elcod(3,1)
       gpcar(3,3,1) =  elcod(3,4)   - elcod(3,1)
       t1          =  gpcar(2,2,1) * gpcar(3,3,1) - gpcar(3,2,1) * gpcar(2,3,1)
       t2          = -gpcar(2,1,1) * gpcar(3,3,1) + gpcar(3,1,1) * gpcar(2,3,1)
       t3          =  gpcar(2,1,1) * gpcar(3,2,1) - gpcar(3,1,1) * gpcar(2,2,1)
       gpdet(1)     =  gpcar(1,1,1) * t1 + gpcar(1,2,1) * t2 + gpcar(1,3,1) * t3

       denom       =  1.0_rp / (sign(1.0_rp,gpdet(1))*max(abs(gpdet(1)),epsilgeo_div))

       xjaci(1,1,1) =  t1 * denom
       xjaci(2,1,1) =  t2 * denom
       xjaci(3,1,1) =  t3 * denom
       xjaci(2,2,1) = ( gpcar(1,1,1) * gpcar(3,3,1) - gpcar(3,1,1) * gpcar(1,3,1)) * denom
       xjaci(3,2,1) = (-gpcar(1,1,1) * gpcar(3,2,1) + gpcar(1,2,1) * gpcar(3,1,1)) * denom
       xjaci(3,3,1) = ( gpcar(1,1,1) * gpcar(2,2,1) - gpcar(2,1,1) * gpcar(1,2,1)) * denom
       xjaci(1,2,1) = (-gpcar(1,2,1) * gpcar(3,3,1) + gpcar(3,2,1) * gpcar(1,3,1)) * denom
       xjaci(1,3,1) = ( gpcar(1,2,1) * gpcar(2,3,1) - gpcar(2,2,1) * gpcar(1,3,1)) * denom
       xjaci(2,3,1) = (-gpcar(1,1,1) * gpcar(2,3,1) + gpcar(2,1,1) * gpcar(1,3,1)) * denom

       gpcar(1,1,1) = -xjaci(1,1,1) - xjaci(2,1,1) - xjaci(3,1,1)
       gpcar(1,2,1) =  xjaci(1,1,1)
       gpcar(1,3,1) =  xjaci(2,1,1)
       gpcar(1,4,1) =  xjaci(3,1,1)
       gpcar(2,1,1) = -xjaci(1,2,1) - xjaci(2,2,1) - xjaci(3,2,1)
       gpcar(2,2,1) =  xjaci(1,2,1)
       gpcar(2,3,1) =  xjaci(2,2,1)
       gpcar(2,4,1) =  xjaci(3,2,1)
       gpcar(3,1,1) = -xjaci(1,3,1) - xjaci(2,3,1) - xjaci(3,3,1)
       gpcar(3,2,1) =  xjaci(1,3,1)
       gpcar(3,3,1) =  xjaci(2,3,1)
       gpcar(3,4,1) =  xjaci(3,3,1)

       xjaci = 0.0_rp
       do igaus = 2,pgaus
          gpdet(igaus) = gpdet(1)
          do inode = 1,4
             do idime = 1,3
                gpcar(idime,inode,igaus) = gpcar(idime,inode,1)
             end do
          end do
       end do

    else if ( ndime == 1 ) then
       !
       ! 1D
       !
       do igaus = 1,pgaus
          xjacm(1,1) = 0.0_rp
          do k = 1,pnode
             xjacm(1,1) = xjacm(1,1) + elcod(1,k) * deriv(1,k,igaus)
          end do
          gpdet(igaus) = xjacm(1,1)
          xjaci(1,1,igaus) = 1.0_rp / (sign(1.0_rp,xjacm(1,1))*max(abs(xjacm(1,1)),epsilgeo_div))
          do j = 1,pnode
             gpcar(1,j,igaus) = xjaci(1,1,igaus) * deriv(1,j,igaus)
          end do
       end do

    else if ( ndime == 2 ) then
       !
       ! 2D
       !
       do igaus = 1,pgaus

          xjacm(1,1) = 0.0_rp
          xjacm(1,2) = 0.0_rp
          xjacm(2,1) = 0.0_rp
          xjacm(2,2) = 0.0_rp
          do k = 1,pnode
             xjacm(1,1) = xjacm(1,1) + elcod(1,k) * deriv(1,k,igaus)
             xjacm(1,2) = xjacm(1,2) + elcod(1,k) * deriv(2,k,igaus)
             xjacm(2,1) = xjacm(2,1) + elcod(2,k) * deriv(1,k,igaus)
             xjacm(2,2) = xjacm(2,2) + elcod(2,k) * deriv(2,k,igaus)
          end do

          gpdet(igaus) = xjacm(1,1) * xjacm(2,2) - xjacm(2,1) * xjacm(1,2)

          denom = 1.0_rp / (sign(1.0_rp,gpdet(igaus))*max(abs(gpdet(igaus)),epsilgeo_div))

          xjaci(1,1,igaus) =  xjacm(2,2) * denom
          xjaci(2,2,igaus) =  xjacm(1,1) * denom
          xjaci(2,1,igaus) = -xjacm(2,1) * denom
          xjaci(1,2,igaus) = -xjacm(1,2) * denom

          do j = 1, pnode
             gpcar(1,j,igaus) =   xjaci(1,1,igaus) * deriv(1,j,igaus) &
                  &                           + xjaci(2,1,igaus) * deriv(2,j,igaus)

             gpcar(2,j,igaus) =   xjaci(1,2,igaus) * deriv(1,j,igaus) &
                  &                           + xjaci(2,2,igaus) * deriv(2,j,igaus)
          end do
       end do

    else if ( ndime == 3 ) then
       !
       ! 3D
       !
       ! xjacm = elcod * deriv^t
       ! xjaci = xjacm^-1
       ! gpcar = xjaci^t * deriv
       !
       do igaus = 1,pgaus

          xjacm(1:3,1:3) = 0.0_rp
          do k = 1,pnode
             xjacm(1,1) = xjacm(1,1) + elcod(1,k) * deriv(1,k,igaus)
             xjacm(1,2) = xjacm(1,2) + elcod(1,k) * deriv(2,k,igaus)
             xjacm(1,3) = xjacm(1,3) + elcod(1,k) * deriv(3,k,igaus)
             xjacm(2,1) = xjacm(2,1) + elcod(2,k) * deriv(1,k,igaus)
             xjacm(2,2) = xjacm(2,2) + elcod(2,k) * deriv(2,k,igaus)
             xjacm(2,3) = xjacm(2,3) + elcod(2,k) * deriv(3,k,igaus)
             xjacm(3,1) = xjacm(3,1) + elcod(3,k) * deriv(1,k,igaus)
             xjacm(3,2) = xjacm(3,2) + elcod(3,k) * deriv(2,k,igaus)
             xjacm(3,3) = xjacm(3,3) + elcod(3,k) * deriv(3,k,igaus)
          end do

          t1              =  xjacm(2,2) * xjacm(3,3) - xjacm(3,2) * xjacm(2,3)
          t2              = -xjacm(2,1) * xjacm(3,3) + xjacm(3,1) * xjacm(2,3)
          t3              =  xjacm(2,1) * xjacm(3,2) - xjacm(3,1) * xjacm(2,2)
          gpdet(igaus)     =  xjacm(1,1) * t1 + xjacm(1,2) * t2 + xjacm(1,3) * t3

          denom           =  1.0_rp / (sign(1.0_rp,gpdet(igaus))*max(abs(gpdet(igaus)),epsilgeo_div))
          xjaci(1,1,igaus) = t1 * denom
          xjaci(2,1,igaus) = t2 * denom
          xjaci(3,1,igaus) = t3 * denom
          xjaci(2,2,igaus) = ( xjacm(1,1) * xjacm(3,3) - xjacm(3,1) * xjacm(1,3)) * denom
          xjaci(3,2,igaus) = (-xjacm(1,1) * xjacm(3,2) + xjacm(1,2) * xjacm(3,1)) * denom
          xjaci(3,3,igaus) = ( xjacm(1,1) * xjacm(2,2) - xjacm(2,1) * xjacm(1,2)) * denom
          xjaci(1,2,igaus) = (-xjacm(1,2) * xjacm(3,3) + xjacm(3,2) * xjacm(1,3)) * denom
          xjaci(1,3,igaus) = ( xjacm(1,2) * xjacm(2,3) - xjacm(2,2) * xjacm(1,3)) * denom
          xjaci(2,3,igaus) = (-xjacm(1,1) * xjacm(2,3) + xjacm(2,1) * xjacm(1,3)) * denom

          do j = 1, pnode
             gpcar(1,j,igaus) =   xjaci(1,1,igaus) * deriv(1,j,igaus) &
                  &                           + xjaci(2,1,igaus) * deriv(2,j,igaus) &
                  &                           + xjaci(3,1,igaus) * deriv(3,j,igaus)

             gpcar(2,j,igaus) =   xjaci(1,2,igaus) * deriv(1,j,igaus) &
                  &                           + xjaci(2,2,igaus) * deriv(2,j,igaus) &
                  &                           + xjaci(3,2,igaus) * deriv(3,j,igaus)

             gpcar(3,j,igaus) =   xjaci(1,3,igaus) * deriv(1,j,igaus) &
                  &                           + xjaci(2,3,igaus) * deriv(2,j,igaus) &
                  &                           + xjaci(3,3,igaus) * deriv(3,j,igaus)
          end do

       end do

    end if

  end subroutine elmgeo_cartesian_derivatives_jacobian_scalar

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    15/09/2016
  !> @brief   Hessian
  !> @details This routine computes the Gphesan matrix in the Cartesian system
  !>          of coordinates according to the rule
  !>          d^2 N / d x_i d x_j
  !>          = (d^2 N / d s_k d s_l) (d s_k/ d x_i)(d s_l/ d x_j)
  !>          + (d N / d s_k) (d^2 s_k / d x_i d x_j)
  !>
  !-----------------------------------------------------------------------

  subroutine elmgeo_hessian_vector(&
       ndime,mnode,ntens,pnode,pgaus,heslo,gphes,xjaci,deriv,elcod)

    integer(ip), intent(in)  :: ndime
    integer(ip), intent(in)  :: mnode
    integer(ip), intent(in)  :: ntens
    integer(ip), intent(in)  :: pnode
    integer(ip), intent(in)  :: pgaus
    real(rp),    intent(in)  :: xjaci(VECTOR_SIZE,ndime,ndime,pgaus)
    real(rp),    intent(in)  :: deriv(ndime,pnode,pgaus)
    real(rp),    intent(in)  :: elcod(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(in)  :: heslo(ntens,pnode,pgaus)
    real(rp),    intent(out) :: gphes(VECTOR_SIZE,ntens,mnode,pgaus)
    real(rp)                 :: d2sdx(VECTOR_SIZE,ndime,ndime,ndime)
    integer(ip)              :: inode,igaus
    real(rp)                 :: heslo1, heslo2, heslo3, heslo4, heslo5, heslo6
    real(rp)                 :: xjaci11(VECTOR_SIZE), xjaci12(VECTOR_SIZE), xjaci13(VECTOR_SIZE)
    real(rp)                 :: xjaci21(VECTOR_SIZE), xjaci22(VECTOR_SIZE), xjaci23(VECTOR_SIZE)
    real(rp)                 :: xjaci31(VECTOR_SIZE), xjaci32(VECTOR_SIZE), xjaci33(VECTOR_SIZE)

    if( ndime == 1 ) then
       !
       ! 1D
       !
       do igaus = 1,pgaus
          do inode = 1,pnode
             gphes(1:VECTOR_SIZE,1,inode,igaus) = xjaci(1:VECTOR_SIZE,1,1,igaus) * heslo(1,inode,igaus) * xjaci(1:VECTOR_SIZE,1,1,igaus)&
                  - deriv(1,inode,igaus) * xjaci(1:VECTOR_SIZE,1,1,igaus) * heslo(1,inode,igaus) * elcod(1:VECTOR_SIZE,1,inode)
          end do
       end do

    else if( ndime == 2 ) then
       !
       ! 2D
       !
       do igaus = 1,pgaus

          xjaci11(1:VECTOR_SIZE) = xjaci(1:VECTOR_SIZE,1,1,igaus)
          xjaci12(1:VECTOR_SIZE) = xjaci(1:VECTOR_SIZE,1,2,igaus)
          xjaci21(1:VECTOR_SIZE) = xjaci(1:VECTOR_SIZE,2,1,igaus)
          xjaci22(1:VECTOR_SIZE) = xjaci(1:VECTOR_SIZE,2,2,igaus)

          do inode = 1,pnode

             gphes(1:VECTOR_SIZE,1,inode,igaus) = 0.0_rp
             gphes(1:VECTOR_SIZE,2,inode,igaus) = 0.0_rp
             gphes(1:VECTOR_SIZE,3,inode,igaus) = 0.0_rp

             heslo1                             = heslo(1,inode,igaus)
             heslo2                             = heslo(2,inode,igaus)
             heslo3                             = heslo(3,inode,igaus)

             gphes(1:VECTOR_SIZE,1,inode,igaus) = gphes(1:VECTOR_SIZE,1,inode,igaus) + xjaci11(1:VECTOR_SIZE) * heslo1 * xjaci11(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,1,inode,igaus) = gphes(1:VECTOR_SIZE,1,inode,igaus) + xjaci11(1:VECTOR_SIZE) * heslo3 * xjaci21(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,1,inode,igaus) = gphes(1:VECTOR_SIZE,1,inode,igaus) + xjaci21(1:VECTOR_SIZE) * heslo3 * xjaci11(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,1,inode,igaus) = gphes(1:VECTOR_SIZE,1,inode,igaus) + xjaci21(1:VECTOR_SIZE) * heslo2 * xjaci21(1:VECTOR_SIZE)

             gphes(1:VECTOR_SIZE,3,inode,igaus) = gphes(1:VECTOR_SIZE,3,inode,igaus) + xjaci11(1:VECTOR_SIZE) * heslo1 * xjaci12(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,3,inode,igaus) = gphes(1:VECTOR_SIZE,3,inode,igaus) + xjaci11(1:VECTOR_SIZE) * heslo3 * xjaci22(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,3,inode,igaus) = gphes(1:VECTOR_SIZE,3,inode,igaus) + xjaci21(1:VECTOR_SIZE) * heslo3 * xjaci12(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,3,inode,igaus) = gphes(1:VECTOR_SIZE,3,inode,igaus) + xjaci21(1:VECTOR_SIZE) * heslo2 * xjaci22(1:VECTOR_SIZE)

             gphes(1:VECTOR_SIZE,2,inode,igaus) = gphes(1:VECTOR_SIZE,2,inode,igaus) + xjaci12(1:VECTOR_SIZE) * heslo1 * xjaci12(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,2,inode,igaus) = gphes(1:VECTOR_SIZE,2,inode,igaus) + xjaci12(1:VECTOR_SIZE) * heslo3 * xjaci22(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,2,inode,igaus) = gphes(1:VECTOR_SIZE,2,inode,igaus) + xjaci22(1:VECTOR_SIZE) * heslo3 * xjaci12(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,2,inode,igaus) = gphes(1:VECTOR_SIZE,2,inode,igaus) + xjaci22(1:VECTOR_SIZE) * heslo2 * xjaci22(1:VECTOR_SIZE)

          end do

          d2sdx(1:VECTOR_SIZE,1,1,1) = 0.0_rp
          d2sdx(1:VECTOR_SIZE,2,1,1) = 0.0_rp
          d2sdx(1:VECTOR_SIZE,1,1,2) = 0.0_rp
          d2sdx(1:VECTOR_SIZE,2,1,2) = 0.0_rp
          d2sdx(1:VECTOR_SIZE,1,2,2) = 0.0_rp
          d2sdx(1:VECTOR_SIZE,2,2,2) = 0.0_rp

          do inode = 1,pnode
             d2sdx(1:VECTOR_SIZE,1,1,1) = d2sdx(1:VECTOR_SIZE,1,1,1) - xjaci11(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,1,inode,igaus) * elcod(1:VECTOR_SIZE,1,inode)
             d2sdx(1:VECTOR_SIZE,1,1,1) = d2sdx(1:VECTOR_SIZE,1,1,1) - xjaci12(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,1,inode,igaus) * elcod(1:VECTOR_SIZE,2,inode)
             d2sdx(1:VECTOR_SIZE,2,1,1) = d2sdx(1:VECTOR_SIZE,2,1,1) - xjaci21(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,1,inode,igaus) * elcod(1:VECTOR_SIZE,1,inode)
             d2sdx(1:VECTOR_SIZE,2,1,1) = d2sdx(1:VECTOR_SIZE,2,1,1) - xjaci22(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,1,inode,igaus) * elcod(1:VECTOR_SIZE,2,inode)
             d2sdx(1:VECTOR_SIZE,1,1,2) = d2sdx(1:VECTOR_SIZE,1,1,2) - xjaci11(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,3,inode,igaus) * elcod(1:VECTOR_SIZE,1,inode)
             d2sdx(1:VECTOR_SIZE,1,1,2) = d2sdx(1:VECTOR_SIZE,1,1,2) - xjaci12(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,3,inode,igaus) * elcod(1:VECTOR_SIZE,2,inode)
             d2sdx(1:VECTOR_SIZE,2,1,2) = d2sdx(1:VECTOR_SIZE,2,1,2) - xjaci21(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,3,inode,igaus) * elcod(1:VECTOR_SIZE,1,inode)
             d2sdx(1:VECTOR_SIZE,2,1,2) = d2sdx(1:VECTOR_SIZE,2,1,2) - xjaci22(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,3,inode,igaus) * elcod(1:VECTOR_SIZE,2,inode)
             d2sdx(1:VECTOR_SIZE,1,2,2) = d2sdx(1:VECTOR_SIZE,1,2,2) - xjaci11(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,2,inode,igaus) * elcod(1:VECTOR_SIZE,1,inode)
             d2sdx(1:VECTOR_SIZE,1,2,2) = d2sdx(1:VECTOR_SIZE,1,2,2) - xjaci12(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,2,inode,igaus) * elcod(1:VECTOR_SIZE,2,inode)
             d2sdx(1:VECTOR_SIZE,2,2,2) = d2sdx(1:VECTOR_SIZE,2,2,2) - xjaci21(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,2,inode,igaus) * elcod(1:VECTOR_SIZE,1,inode)
             d2sdx(1:VECTOR_SIZE,2,2,2) = d2sdx(1:VECTOR_SIZE,2,2,2) - xjaci22(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,2,inode,igaus) * elcod(1:VECTOR_SIZE,2,inode)
          end do
          !
          ! Computes the second Cartesian derivatives of the shape functions
          !
          do inode = 1,pnode

             gphes(1:VECTOR_SIZE,1,inode,igaus) = gphes(1:VECTOR_SIZE,1,inode,igaus) + deriv(1,inode,igaus) * d2sdx(1:VECTOR_SIZE,1,1,1)
             gphes(1:VECTOR_SIZE,1,inode,igaus) = gphes(1:VECTOR_SIZE,1,inode,igaus) + deriv(2,inode,igaus) * d2sdx(1:VECTOR_SIZE,2,1,1)

             gphes(1:VECTOR_SIZE,3,inode,igaus) = gphes(1:VECTOR_SIZE,3,inode,igaus) + deriv(1,inode,igaus) * d2sdx(1:VECTOR_SIZE,1,1,2)
             gphes(1:VECTOR_SIZE,3,inode,igaus) = gphes(1:VECTOR_SIZE,3,inode,igaus) + deriv(2,inode,igaus) * d2sdx(1:VECTOR_SIZE,2,1,2)

             gphes(1:VECTOR_SIZE,2,inode,igaus) = gphes(1:VECTOR_SIZE,2,inode,igaus) + deriv(1,inode,igaus) * d2sdx(1:VECTOR_SIZE,1,2,2)
             gphes(1:VECTOR_SIZE,2,inode,igaus) = gphes(1:VECTOR_SIZE,2,inode,igaus) + deriv(2,inode,igaus) * d2sdx(1:VECTOR_SIZE,2,2,2)

          end do

       end do

    else
       !
       ! 3D
       !
       do igaus  =  1,pgaus

          xjaci11(1:VECTOR_SIZE)  =  xjaci(1:VECTOR_SIZE,1,1,igaus)
          xjaci12(1:VECTOR_SIZE)  =  xjaci(1:VECTOR_SIZE,1,2,igaus)
          xjaci13(1:VECTOR_SIZE)  =  xjaci(1:VECTOR_SIZE,1,3,igaus)
          xjaci21(1:VECTOR_SIZE)  =  xjaci(1:VECTOR_SIZE,2,1,igaus)
          xjaci22(1:VECTOR_SIZE)  =  xjaci(1:VECTOR_SIZE,2,2,igaus)
          xjaci23(1:VECTOR_SIZE)  =  xjaci(1:VECTOR_SIZE,2,3,igaus)
          xjaci31(1:VECTOR_SIZE)  =  xjaci(1:VECTOR_SIZE,3,1,igaus)
          xjaci32(1:VECTOR_SIZE)  =  xjaci(1:VECTOR_SIZE,3,2,igaus)
          xjaci33(1:VECTOR_SIZE)  =  xjaci(1:VECTOR_SIZE,3,3,igaus)

          do inode=1,pnode

             gphes(1:VECTOR_SIZE,1,inode,igaus) = 0.0_rp
             gphes(1:VECTOR_SIZE,4,inode,igaus) = 0.0_rp
             gphes(1:VECTOR_SIZE,5,inode,igaus) = 0.0_rp
             gphes(1:VECTOR_SIZE,2,inode,igaus) = 0.0_rp
             gphes(1:VECTOR_SIZE,6,inode,igaus) = 0.0_rp
             gphes(1:VECTOR_SIZE,3,inode,igaus) = 0.0_rp

             heslo1                             = heslo(1,inode,igaus)
             heslo2                             = heslo(2,inode,igaus)
             heslo3                             = heslo(3,inode,igaus)
             heslo4                             = heslo(4,inode,igaus)
             heslo5                             = heslo(5,inode,igaus)
             heslo6                             = heslo(6,inode,igaus)

             gphes(1:VECTOR_SIZE,1,inode,igaus) = gphes(1:VECTOR_SIZE,1,inode,igaus) + xjaci11(1:VECTOR_SIZE) * heslo1 * xjaci11(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,1,inode,igaus) = gphes(1:VECTOR_SIZE,1,inode,igaus) + xjaci11(1:VECTOR_SIZE) * heslo4 * xjaci21(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,1,inode,igaus) = gphes(1:VECTOR_SIZE,1,inode,igaus) + xjaci11(1:VECTOR_SIZE) * heslo5 * xjaci31(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,1,inode,igaus) = gphes(1:VECTOR_SIZE,1,inode,igaus) + xjaci21(1:VECTOR_SIZE) * heslo4 * xjaci11(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,1,inode,igaus) = gphes(1:VECTOR_SIZE,1,inode,igaus) + xjaci21(1:VECTOR_SIZE) * heslo2 * xjaci21(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,1,inode,igaus) = gphes(1:VECTOR_SIZE,1,inode,igaus) + xjaci21(1:VECTOR_SIZE) * heslo6 * xjaci31(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,1,inode,igaus) = gphes(1:VECTOR_SIZE,1,inode,igaus) + xjaci31(1:VECTOR_SIZE) * heslo5 * xjaci11(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,1,inode,igaus) = gphes(1:VECTOR_SIZE,1,inode,igaus) + xjaci31(1:VECTOR_SIZE) * heslo6 * xjaci21(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,1,inode,igaus) = gphes(1:VECTOR_SIZE,1,inode,igaus) + xjaci31(1:VECTOR_SIZE) * heslo3 * xjaci31(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,4,inode,igaus) = gphes(1:VECTOR_SIZE,4,inode,igaus) + xjaci11(1:VECTOR_SIZE) * heslo1 * xjaci12(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,4,inode,igaus) = gphes(1:VECTOR_SIZE,4,inode,igaus) + xjaci11(1:VECTOR_SIZE) * heslo4 * xjaci22(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,4,inode,igaus) = gphes(1:VECTOR_SIZE,4,inode,igaus) + xjaci11(1:VECTOR_SIZE) * heslo5 * xjaci32(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,4,inode,igaus) = gphes(1:VECTOR_SIZE,4,inode,igaus) + xjaci21(1:VECTOR_SIZE) * heslo4 * xjaci12(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,4,inode,igaus) = gphes(1:VECTOR_SIZE,4,inode,igaus) + xjaci21(1:VECTOR_SIZE) * heslo2 * xjaci22(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,4,inode,igaus) = gphes(1:VECTOR_SIZE,4,inode,igaus) + xjaci21(1:VECTOR_SIZE) * heslo6 * xjaci32(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,4,inode,igaus) = gphes(1:VECTOR_SIZE,4,inode,igaus) + xjaci31(1:VECTOR_SIZE) * heslo5 * xjaci12(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,4,inode,igaus) = gphes(1:VECTOR_SIZE,4,inode,igaus) + xjaci31(1:VECTOR_SIZE) * heslo6 * xjaci22(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,4,inode,igaus) = gphes(1:VECTOR_SIZE,4,inode,igaus) + xjaci31(1:VECTOR_SIZE) * heslo3 * xjaci32(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,5,inode,igaus) = gphes(1:VECTOR_SIZE,5,inode,igaus) + xjaci11(1:VECTOR_SIZE) * heslo1 * xjaci13(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,5,inode,igaus) = gphes(1:VECTOR_SIZE,5,inode,igaus) + xjaci11(1:VECTOR_SIZE) * heslo4 * xjaci23(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,5,inode,igaus) = gphes(1:VECTOR_SIZE,5,inode,igaus) + xjaci11(1:VECTOR_SIZE) * heslo5 * xjaci33(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,5,inode,igaus) = gphes(1:VECTOR_SIZE,5,inode,igaus) + xjaci21(1:VECTOR_SIZE) * heslo4 * xjaci13(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,5,inode,igaus) = gphes(1:VECTOR_SIZE,5,inode,igaus) + xjaci21(1:VECTOR_SIZE) * heslo2 * xjaci23(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,5,inode,igaus) = gphes(1:VECTOR_SIZE,5,inode,igaus) + xjaci21(1:VECTOR_SIZE) * heslo6 * xjaci33(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,5,inode,igaus) = gphes(1:VECTOR_SIZE,5,inode,igaus) + xjaci31(1:VECTOR_SIZE) * heslo5 * xjaci13(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,5,inode,igaus) = gphes(1:VECTOR_SIZE,5,inode,igaus) + xjaci31(1:VECTOR_SIZE) * heslo6 * xjaci23(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,5,inode,igaus) = gphes(1:VECTOR_SIZE,5,inode,igaus) + xjaci31(1:VECTOR_SIZE) * heslo3 * xjaci33(1:VECTOR_SIZE)

             gphes(1:VECTOR_SIZE,2,inode,igaus) = gphes(1:VECTOR_SIZE,2,inode,igaus) + xjaci12(1:VECTOR_SIZE) * heslo1 * xjaci12(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,2,inode,igaus) = gphes(1:VECTOR_SIZE,2,inode,igaus) + xjaci12(1:VECTOR_SIZE) * heslo4 * xjaci22(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,2,inode,igaus) = gphes(1:VECTOR_SIZE,2,inode,igaus) + xjaci12(1:VECTOR_SIZE) * heslo5 * xjaci32(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,2,inode,igaus) = gphes(1:VECTOR_SIZE,2,inode,igaus) + xjaci22(1:VECTOR_SIZE) * heslo4 * xjaci12(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,2,inode,igaus) = gphes(1:VECTOR_SIZE,2,inode,igaus) + xjaci22(1:VECTOR_SIZE) * heslo2 * xjaci22(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,2,inode,igaus) = gphes(1:VECTOR_SIZE,2,inode,igaus) + xjaci22(1:VECTOR_SIZE) * heslo6 * xjaci32(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,2,inode,igaus) = gphes(1:VECTOR_SIZE,2,inode,igaus) + xjaci32(1:VECTOR_SIZE) * heslo5 * xjaci12(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,2,inode,igaus) = gphes(1:VECTOR_SIZE,2,inode,igaus) + xjaci32(1:VECTOR_SIZE) * heslo6 * xjaci22(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,2,inode,igaus) = gphes(1:VECTOR_SIZE,2,inode,igaus) + xjaci32(1:VECTOR_SIZE) * heslo3 * xjaci32(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,6,inode,igaus) = gphes(1:VECTOR_SIZE,6,inode,igaus) + xjaci12(1:VECTOR_SIZE) * heslo1 * xjaci13(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,6,inode,igaus) = gphes(1:VECTOR_SIZE,6,inode,igaus) + xjaci12(1:VECTOR_SIZE) * heslo4 * xjaci23(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,6,inode,igaus) = gphes(1:VECTOR_SIZE,6,inode,igaus) + xjaci12(1:VECTOR_SIZE) * heslo5 * xjaci33(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,6,inode,igaus) = gphes(1:VECTOR_SIZE,6,inode,igaus) + xjaci22(1:VECTOR_SIZE) * heslo4 * xjaci13(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,6,inode,igaus) = gphes(1:VECTOR_SIZE,6,inode,igaus) + xjaci22(1:VECTOR_SIZE) * heslo2 * xjaci23(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,6,inode,igaus) = gphes(1:VECTOR_SIZE,6,inode,igaus) + xjaci22(1:VECTOR_SIZE) * heslo6 * xjaci33(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,6,inode,igaus) = gphes(1:VECTOR_SIZE,6,inode,igaus) + xjaci32(1:VECTOR_SIZE) * heslo5 * xjaci13(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,6,inode,igaus) = gphes(1:VECTOR_SIZE,6,inode,igaus) + xjaci32(1:VECTOR_SIZE) * heslo6 * xjaci23(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,6,inode,igaus) = gphes(1:VECTOR_SIZE,6,inode,igaus) + xjaci32(1:VECTOR_SIZE) * heslo3 * xjaci33(1:VECTOR_SIZE)

             gphes(1:VECTOR_SIZE,3,inode,igaus) = gphes(1:VECTOR_SIZE,3,inode,igaus) + xjaci13(1:VECTOR_SIZE) * heslo1 * xjaci13(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,3,inode,igaus) = gphes(1:VECTOR_SIZE,3,inode,igaus) + xjaci13(1:VECTOR_SIZE) * heslo4 * xjaci23(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,3,inode,igaus) = gphes(1:VECTOR_SIZE,3,inode,igaus) + xjaci13(1:VECTOR_SIZE) * heslo5 * xjaci33(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,3,inode,igaus) = gphes(1:VECTOR_SIZE,3,inode,igaus) + xjaci23(1:VECTOR_SIZE) * heslo4 * xjaci13(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,3,inode,igaus) = gphes(1:VECTOR_SIZE,3,inode,igaus) + xjaci23(1:VECTOR_SIZE) * heslo2 * xjaci23(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,3,inode,igaus) = gphes(1:VECTOR_SIZE,3,inode,igaus) + xjaci23(1:VECTOR_SIZE) * heslo6 * xjaci33(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,3,inode,igaus) = gphes(1:VECTOR_SIZE,3,inode,igaus) + xjaci33(1:VECTOR_SIZE) * heslo5 * xjaci13(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,3,inode,igaus) = gphes(1:VECTOR_SIZE,3,inode,igaus) + xjaci33(1:VECTOR_SIZE) * heslo6 * xjaci23(1:VECTOR_SIZE)
             gphes(1:VECTOR_SIZE,3,inode,igaus) = gphes(1:VECTOR_SIZE,3,inode,igaus) + xjaci33(1:VECTOR_SIZE) * heslo3 * xjaci33(1:VECTOR_SIZE)

          end do

          d2sdx(1:VECTOR_SIZE,1,1,1) = 0.0_rp
          d2sdx(1:VECTOR_SIZE,2,1,1) = 0.0_rp
          d2sdx(1:VECTOR_SIZE,3,1,1) = 0.0_rp
          d2sdx(1:VECTOR_SIZE,1,1,2) = 0.0_rp
          d2sdx(1:VECTOR_SIZE,2,1,2) = 0.0_rp
          d2sdx(1:VECTOR_SIZE,3,1,2) = 0.0_rp
          d2sdx(1:VECTOR_SIZE,1,2,2) = 0.0_rp
          d2sdx(1:VECTOR_SIZE,2,2,2) = 0.0_rp
          d2sdx(1:VECTOR_SIZE,3,2,2) = 0.0_rp
          d2sdx(1:VECTOR_SIZE,1,1,3) = 0.0_rp
          d2sdx(1:VECTOR_SIZE,2,1,3) = 0.0_rp
          d2sdx(1:VECTOR_SIZE,3,1,3) = 0.0_rp
          d2sdx(1:VECTOR_SIZE,1,2,3) = 0.0_rp
          d2sdx(1:VECTOR_SIZE,2,2,3) = 0.0_rp
          d2sdx(1:VECTOR_SIZE,3,2,3) = 0.0_rp
          d2sdx(1:VECTOR_SIZE,1,3,3) = 0.0_rp
          d2sdx(1:VECTOR_SIZE,2,3,3) = 0.0_rp
          d2sdx(1:VECTOR_SIZE,3,3,3) = 0.0_rp

          do inode = 1,pnode
             d2sdx(1:VECTOR_SIZE,1,1,1) = d2sdx(1:VECTOR_SIZE,1,1,1) - xjaci11(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,1,inode,igaus) * elcod(1:VECTOR_SIZE,1,inode)
             d2sdx(1:VECTOR_SIZE,1,1,1) = d2sdx(1:VECTOR_SIZE,1,1,1) - xjaci12(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,1,inode,igaus) * elcod(1:VECTOR_SIZE,2,inode)
             d2sdx(1:VECTOR_SIZE,1,1,1) = d2sdx(1:VECTOR_SIZE,1,1,1) - xjaci13(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,1,inode,igaus) * elcod(1:VECTOR_SIZE,3,inode)
             d2sdx(1:VECTOR_SIZE,2,1,1) = d2sdx(1:VECTOR_SIZE,2,1,1) - xjaci21(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,1,inode,igaus) * elcod(1:VECTOR_SIZE,1,inode)
             d2sdx(1:VECTOR_SIZE,2,1,1) = d2sdx(1:VECTOR_SIZE,2,1,1) - xjaci22(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,1,inode,igaus) * elcod(1:VECTOR_SIZE,2,inode)
             d2sdx(1:VECTOR_SIZE,2,1,1) = d2sdx(1:VECTOR_SIZE,2,1,1) - xjaci23(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,1,inode,igaus) * elcod(1:VECTOR_SIZE,3,inode)
             d2sdx(1:VECTOR_SIZE,3,1,1) = d2sdx(1:VECTOR_SIZE,3,1,1) - xjaci31(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,1,inode,igaus) * elcod(1:VECTOR_SIZE,1,inode)
             d2sdx(1:VECTOR_SIZE,3,1,1) = d2sdx(1:VECTOR_SIZE,3,1,1) - xjaci32(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,1,inode,igaus) * elcod(1:VECTOR_SIZE,2,inode)
             d2sdx(1:VECTOR_SIZE,3,1,1) = d2sdx(1:VECTOR_SIZE,3,1,1) - xjaci33(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,1,inode,igaus) * elcod(1:VECTOR_SIZE,3,inode)
             d2sdx(1:VECTOR_SIZE,1,1,2) = d2sdx(1:VECTOR_SIZE,1,1,2) - xjaci11(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,4,inode,igaus) * elcod(1:VECTOR_SIZE,1,inode)
             d2sdx(1:VECTOR_SIZE,1,1,2) = d2sdx(1:VECTOR_SIZE,1,1,2) - xjaci12(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,4,inode,igaus) * elcod(1:VECTOR_SIZE,2,inode)
             d2sdx(1:VECTOR_SIZE,1,1,2) = d2sdx(1:VECTOR_SIZE,1,1,2) - xjaci13(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,4,inode,igaus) * elcod(1:VECTOR_SIZE,3,inode)
             d2sdx(1:VECTOR_SIZE,2,1,2) = d2sdx(1:VECTOR_SIZE,2,1,2) - xjaci21(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,4,inode,igaus) * elcod(1:VECTOR_SIZE,1,inode)
             d2sdx(1:VECTOR_SIZE,2,1,2) = d2sdx(1:VECTOR_SIZE,2,1,2) - xjaci22(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,4,inode,igaus) * elcod(1:VECTOR_SIZE,2,inode)
             d2sdx(1:VECTOR_SIZE,2,1,2) = d2sdx(1:VECTOR_SIZE,2,1,2) - xjaci23(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,4,inode,igaus) * elcod(1:VECTOR_SIZE,3,inode)
             d2sdx(1:VECTOR_SIZE,3,1,2) = d2sdx(1:VECTOR_SIZE,3,1,2) - xjaci31(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,4,inode,igaus) * elcod(1:VECTOR_SIZE,1,inode)
             d2sdx(1:VECTOR_SIZE,3,1,2) = d2sdx(1:VECTOR_SIZE,3,1,2) - xjaci32(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,4,inode,igaus) * elcod(1:VECTOR_SIZE,2,inode)
             d2sdx(1:VECTOR_SIZE,3,1,2) = d2sdx(1:VECTOR_SIZE,3,1,2) - xjaci33(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,4,inode,igaus) * elcod(1:VECTOR_SIZE,3,inode)
             d2sdx(1:VECTOR_SIZE,1,2,2) = d2sdx(1:VECTOR_SIZE,1,2,2) - xjaci11(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,2,inode,igaus) * elcod(1:VECTOR_SIZE,1,inode)
             d2sdx(1:VECTOR_SIZE,1,2,2) = d2sdx(1:VECTOR_SIZE,1,2,2) - xjaci12(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,2,inode,igaus) * elcod(1:VECTOR_SIZE,2,inode)
             d2sdx(1:VECTOR_SIZE,1,2,2) = d2sdx(1:VECTOR_SIZE,1,2,2) - xjaci13(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,2,inode,igaus) * elcod(1:VECTOR_SIZE,3,inode)
             d2sdx(1:VECTOR_SIZE,2,2,2) = d2sdx(1:VECTOR_SIZE,2,2,2) - xjaci21(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,2,inode,igaus) * elcod(1:VECTOR_SIZE,1,inode)
             d2sdx(1:VECTOR_SIZE,2,2,2) = d2sdx(1:VECTOR_SIZE,2,2,2) - xjaci22(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,2,inode,igaus) * elcod(1:VECTOR_SIZE,2,inode)
             d2sdx(1:VECTOR_SIZE,2,2,2) = d2sdx(1:VECTOR_SIZE,2,2,2) - xjaci23(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,2,inode,igaus) * elcod(1:VECTOR_SIZE,3,inode)
             d2sdx(1:VECTOR_SIZE,3,2,2) = d2sdx(1:VECTOR_SIZE,3,2,2) - xjaci31(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,2,inode,igaus) * elcod(1:VECTOR_SIZE,1,inode)
             d2sdx(1:VECTOR_SIZE,3,2,2) = d2sdx(1:VECTOR_SIZE,3,2,2) - xjaci32(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,2,inode,igaus) * elcod(1:VECTOR_SIZE,2,inode)
             d2sdx(1:VECTOR_SIZE,3,2,2) = d2sdx(1:VECTOR_SIZE,3,2,2) - xjaci33(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,2,inode,igaus) * elcod(1:VECTOR_SIZE,3,inode)
             d2sdx(1:VECTOR_SIZE,1,1,3) = d2sdx(1:VECTOR_SIZE,1,1,3) - xjaci11(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,5,inode,igaus) * elcod(1:VECTOR_SIZE,1,inode)
             d2sdx(1:VECTOR_SIZE,1,1,3) = d2sdx(1:VECTOR_SIZE,1,1,3) - xjaci12(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,5,inode,igaus) * elcod(1:VECTOR_SIZE,2,inode)
             d2sdx(1:VECTOR_SIZE,1,1,3) = d2sdx(1:VECTOR_SIZE,1,1,3) - xjaci13(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,5,inode,igaus) * elcod(1:VECTOR_SIZE,3,inode)
             d2sdx(1:VECTOR_SIZE,2,1,3) = d2sdx(1:VECTOR_SIZE,2,1,3) - xjaci21(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,5,inode,igaus) * elcod(1:VECTOR_SIZE,1,inode)
             d2sdx(1:VECTOR_SIZE,2,1,3) = d2sdx(1:VECTOR_SIZE,2,1,3) - xjaci22(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,5,inode,igaus) * elcod(1:VECTOR_SIZE,2,inode)
             d2sdx(1:VECTOR_SIZE,2,1,3) = d2sdx(1:VECTOR_SIZE,2,1,3) - xjaci23(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,5,inode,igaus) * elcod(1:VECTOR_SIZE,3,inode)
             d2sdx(1:VECTOR_SIZE,3,1,3) = d2sdx(1:VECTOR_SIZE,3,1,3) - xjaci31(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,5,inode,igaus) * elcod(1:VECTOR_SIZE,1,inode)
             d2sdx(1:VECTOR_SIZE,3,1,3) = d2sdx(1:VECTOR_SIZE,3,1,3) - xjaci32(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,5,inode,igaus) * elcod(1:VECTOR_SIZE,2,inode)
             d2sdx(1:VECTOR_SIZE,3,1,3) = d2sdx(1:VECTOR_SIZE,3,1,3) - xjaci33(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,5,inode,igaus) * elcod(1:VECTOR_SIZE,3,inode)
             d2sdx(1:VECTOR_SIZE,1,2,3) = d2sdx(1:VECTOR_SIZE,1,2,3) - xjaci11(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,6,inode,igaus) * elcod(1:VECTOR_SIZE,1,inode)
             d2sdx(1:VECTOR_SIZE,1,2,3) = d2sdx(1:VECTOR_SIZE,1,2,3) - xjaci12(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,6,inode,igaus) * elcod(1:VECTOR_SIZE,2,inode)
             d2sdx(1:VECTOR_SIZE,1,2,3) = d2sdx(1:VECTOR_SIZE,1,2,3) - xjaci13(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,6,inode,igaus) * elcod(1:VECTOR_SIZE,3,inode)
             d2sdx(1:VECTOR_SIZE,2,2,3) = d2sdx(1:VECTOR_SIZE,2,2,3) - xjaci21(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,6,inode,igaus) * elcod(1:VECTOR_SIZE,1,inode)
             d2sdx(1:VECTOR_SIZE,2,2,3) = d2sdx(1:VECTOR_SIZE,2,2,3) - xjaci22(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,6,inode,igaus) * elcod(1:VECTOR_SIZE,2,inode)
             d2sdx(1:VECTOR_SIZE,2,2,3) = d2sdx(1:VECTOR_SIZE,2,2,3) - xjaci23(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,6,inode,igaus) * elcod(1:VECTOR_SIZE,3,inode)
             d2sdx(1:VECTOR_SIZE,3,2,3) = d2sdx(1:VECTOR_SIZE,3,2,3) - xjaci31(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,6,inode,igaus) * elcod(1:VECTOR_SIZE,1,inode)
             d2sdx(1:VECTOR_SIZE,3,2,3) = d2sdx(1:VECTOR_SIZE,3,2,3) - xjaci32(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,6,inode,igaus) * elcod(1:VECTOR_SIZE,2,inode)
             d2sdx(1:VECTOR_SIZE,3,2,3) = d2sdx(1:VECTOR_SIZE,3,2,3) - xjaci33(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,6,inode,igaus) * elcod(1:VECTOR_SIZE,3,inode)
             d2sdx(1:VECTOR_SIZE,1,3,3) = d2sdx(1:VECTOR_SIZE,1,3,3) - xjaci11(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,3,inode,igaus) * elcod(1:VECTOR_SIZE,1,inode)
             d2sdx(1:VECTOR_SIZE,1,3,3) = d2sdx(1:VECTOR_SIZE,1,3,3) - xjaci12(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,3,inode,igaus) * elcod(1:VECTOR_SIZE,2,inode)
             d2sdx(1:VECTOR_SIZE,1,3,3) = d2sdx(1:VECTOR_SIZE,1,3,3) - xjaci13(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,3,inode,igaus) * elcod(1:VECTOR_SIZE,3,inode)
             d2sdx(1:VECTOR_SIZE,2,3,3) = d2sdx(1:VECTOR_SIZE,2,3,3) - xjaci21(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,3,inode,igaus) * elcod(1:VECTOR_SIZE,1,inode)
             d2sdx(1:VECTOR_SIZE,2,3,3) = d2sdx(1:VECTOR_SIZE,2,3,3) - xjaci22(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,3,inode,igaus) * elcod(1:VECTOR_SIZE,2,inode)
             d2sdx(1:VECTOR_SIZE,2,3,3) = d2sdx(1:VECTOR_SIZE,2,3,3) - xjaci23(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,3,inode,igaus) * elcod(1:VECTOR_SIZE,3,inode)
             d2sdx(1:VECTOR_SIZE,3,3,3) = d2sdx(1:VECTOR_SIZE,3,3,3) - xjaci31(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,3,inode,igaus) * elcod(1:VECTOR_SIZE,1,inode)
             d2sdx(1:VECTOR_SIZE,3,3,3) = d2sdx(1:VECTOR_SIZE,3,3,3) - xjaci32(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,3,inode,igaus) * elcod(1:VECTOR_SIZE,2,inode)
             d2sdx(1:VECTOR_SIZE,3,3,3) = d2sdx(1:VECTOR_SIZE,3,3,3) - xjaci33(1:VECTOR_SIZE) * gphes(1:VECTOR_SIZE,3,inode,igaus) * elcod(1:VECTOR_SIZE,3,inode)
          end do
          !
          ! Computes the second Cartesian derivatives of the shape functions
          !
          do inode  =  1,pnode

             gphes(1:VECTOR_SIZE,1,inode,igaus) = gphes(1:VECTOR_SIZE,1,inode,igaus) + deriv(1,inode,igaus) * d2sdx(1:VECTOR_SIZE,1,1,1)
             gphes(1:VECTOR_SIZE,1,inode,igaus) = gphes(1:VECTOR_SIZE,1,inode,igaus) + deriv(2,inode,igaus) * d2sdx(1:VECTOR_SIZE,2,1,1)
             gphes(1:VECTOR_SIZE,1,inode,igaus) = gphes(1:VECTOR_SIZE,1,inode,igaus) + deriv(3,inode,igaus) * d2sdx(1:VECTOR_SIZE,3,1,1)
             gphes(1:VECTOR_SIZE,4,inode,igaus) = gphes(1:VECTOR_SIZE,4,inode,igaus) + deriv(1,inode,igaus) * d2sdx(1:VECTOR_SIZE,1,1,2)
             gphes(1:VECTOR_SIZE,4,inode,igaus) = gphes(1:VECTOR_SIZE,4,inode,igaus) + deriv(2,inode,igaus) * d2sdx(1:VECTOR_SIZE,2,1,2)
             gphes(1:VECTOR_SIZE,4,inode,igaus) = gphes(1:VECTOR_SIZE,4,inode,igaus) + deriv(3,inode,igaus) * d2sdx(1:VECTOR_SIZE,3,1,2)
             gphes(1:VECTOR_SIZE,5,inode,igaus) = gphes(1:VECTOR_SIZE,5,inode,igaus) + deriv(1,inode,igaus) * d2sdx(1:VECTOR_SIZE,1,1,3)
             gphes(1:VECTOR_SIZE,5,inode,igaus) = gphes(1:VECTOR_SIZE,5,inode,igaus) + deriv(2,inode,igaus) * d2sdx(1:VECTOR_SIZE,2,1,3)
             gphes(1:VECTOR_SIZE,5,inode,igaus) = gphes(1:VECTOR_SIZE,5,inode,igaus) + deriv(3,inode,igaus) * d2sdx(1:VECTOR_SIZE,3,1,3)

             gphes(1:VECTOR_SIZE,2,inode,igaus) = gphes(1:VECTOR_SIZE,2,inode,igaus) + deriv(1,inode,igaus) * d2sdx(1:VECTOR_SIZE,1,2,2)
             gphes(1:VECTOR_SIZE,2,inode,igaus) = gphes(1:VECTOR_SIZE,2,inode,igaus) + deriv(2,inode,igaus) * d2sdx(1:VECTOR_SIZE,2,2,2)
             gphes(1:VECTOR_SIZE,2,inode,igaus) = gphes(1:VECTOR_SIZE,2,inode,igaus) + deriv(3,inode,igaus) * d2sdx(1:VECTOR_SIZE,3,2,2)
             gphes(1:VECTOR_SIZE,6,inode,igaus) = gphes(1:VECTOR_SIZE,6,inode,igaus) + deriv(1,inode,igaus) * d2sdx(1:VECTOR_SIZE,1,2,3)
             gphes(1:VECTOR_SIZE,6,inode,igaus) = gphes(1:VECTOR_SIZE,6,inode,igaus) + deriv(2,inode,igaus) * d2sdx(1:VECTOR_SIZE,2,2,3)
             gphes(1:VECTOR_SIZE,6,inode,igaus) = gphes(1:VECTOR_SIZE,6,inode,igaus) + deriv(3,inode,igaus) * d2sdx(1:VECTOR_SIZE,3,2,3)

             gphes(1:VECTOR_SIZE,3,inode,igaus) = gphes(1:VECTOR_SIZE,3,inode,igaus) + deriv(1,inode,igaus) * d2sdx(1:VECTOR_SIZE,1,3,3)
             gphes(1:VECTOR_SIZE,3,inode,igaus) = gphes(1:VECTOR_SIZE,3,inode,igaus) + deriv(2,inode,igaus) * d2sdx(1:VECTOR_SIZE,2,3,3)
             gphes(1:VECTOR_SIZE,3,inode,igaus) = gphes(1:VECTOR_SIZE,3,inode,igaus) + deriv(3,inode,igaus) * d2sdx(1:VECTOR_SIZE,3,3,3)

          end do

       end do

    end if

  end subroutine elmgeo_hessian_vector

  subroutine elmgeo_hessian_scalar(&
       ndime,mnode,ntens,pnode,pgaus,heslo,gphes,xjaci,deriv,elcod)

    integer(ip), intent(in)  :: ndime
    integer(ip), intent(in)  :: mnode
    integer(ip), intent(in)  :: ntens
    integer(ip), intent(in)  :: pnode
    integer(ip), intent(in)  :: pgaus
    real(rp),    intent(in)  :: xjaci(ndime,ndime,pgaus)
    real(rp),    intent(in)  :: deriv(ndime,pnode,pgaus)
    real(rp),    intent(in)  :: elcod(ndime,pnode)
    real(rp),    intent(in)  :: heslo(ntens,pnode,pgaus)
    real(rp),    intent(out) :: gphes(ntens,mnode,pgaus)
    real(rp)                 :: d2sdx(ndime,ndime,ndime)
    integer(ip)              :: inode,igaus
    real(rp)                 :: heslo1, heslo2, heslo3, heslo4, heslo5, heslo6
    real(rp)                 :: xjaci11, xjaci12, xjaci13
    real(rp)                 :: xjaci21, xjaci22, xjaci23
    real(rp)                 :: xjaci31, xjaci32, xjaci33

    if( ndime == 1 ) then
       !
       ! 1D
       !
       do igaus = 1,pgaus
          do inode = 1,pnode
             gphes(1,inode,igaus) = xjaci(1,1,igaus)  *  heslo(1,inode,igaus)  *  xjaci(1,1,igaus)&
                  - deriv(1,inode,igaus)  *  xjaci(1,1,igaus)  *  heslo(1,inode,igaus)  *  elcod(1,inode)
          end do
       end do

    else if( ndime == 2 ) then
       !
       ! 2D
       !
       do igaus = 1,pgaus

          xjaci11 = xjaci(1,1,igaus)
          xjaci12 = xjaci(1,2,igaus)
          xjaci21 = xjaci(2,1,igaus)
          xjaci22 = xjaci(2,2,igaus)

          do inode = 1,pnode

             gphes(1,inode,igaus) = 0.0_rp
             gphes(2,inode,igaus) = 0.0_rp
             gphes(3,inode,igaus) = 0.0_rp

             heslo1               = heslo(1,inode,igaus)
             heslo2               = heslo(2,inode,igaus)
             heslo3               = heslo(3,inode,igaus)

             gphes(1,inode,igaus) = gphes(1,inode,igaus) + xjaci11 * heslo1 * xjaci11
             gphes(1,inode,igaus) = gphes(1,inode,igaus) + xjaci11 * heslo3 * xjaci21
             gphes(1,inode,igaus) = gphes(1,inode,igaus) + xjaci21 * heslo3 * xjaci11
             gphes(1,inode,igaus) = gphes(1,inode,igaus) + xjaci21 * heslo2 * xjaci21

             gphes(3,inode,igaus) = gphes(3,inode,igaus) + xjaci11 * heslo1 * xjaci12
             gphes(3,inode,igaus) = gphes(3,inode,igaus) + xjaci11 * heslo3 * xjaci22
             gphes(3,inode,igaus) = gphes(3,inode,igaus) + xjaci21 * heslo3 * xjaci12
             gphes(3,inode,igaus) = gphes(3,inode,igaus) + xjaci21 * heslo2 * xjaci22

             gphes(2,inode,igaus) = gphes(2,inode,igaus) + xjaci12 * heslo1 * xjaci12
             gphes(2,inode,igaus) = gphes(2,inode,igaus) + xjaci12 * heslo3 * xjaci22
             gphes(2,inode,igaus) = gphes(2,inode,igaus) + xjaci22 * heslo3 * xjaci12
             gphes(2,inode,igaus) = gphes(2,inode,igaus) + xjaci22 * heslo2 * xjaci22

          end do

          d2sdx(1,1,1) = 0.0_rp
          d2sdx(2,1,1) = 0.0_rp
          d2sdx(1,1,2) = 0.0_rp
          d2sdx(2,1,2) = 0.0_rp
          d2sdx(1,2,2) = 0.0_rp
          d2sdx(2,2,2) = 0.0_rp

          do inode = 1,pnode
             d2sdx(1,1,1) = d2sdx(1,1,1) - xjaci11 * gphes(1,inode,igaus) * elcod(1,inode)
             d2sdx(1,1,1) = d2sdx(1,1,1) - xjaci12 * gphes(1,inode,igaus) * elcod(2,inode)
             d2sdx(2,1,1) = d2sdx(2,1,1) - xjaci21 * gphes(1,inode,igaus) * elcod(1,inode)
             d2sdx(2,1,1) = d2sdx(2,1,1) - xjaci22 * gphes(1,inode,igaus) * elcod(2,inode)
             d2sdx(1,1,2) = d2sdx(1,1,2) - xjaci11 * gphes(3,inode,igaus) * elcod(1,inode)
             d2sdx(1,1,2) = d2sdx(1,1,2) - xjaci12 * gphes(3,inode,igaus) * elcod(2,inode)
             d2sdx(2,1,2) = d2sdx(2,1,2) - xjaci21 * gphes(3,inode,igaus) * elcod(1,inode)
             d2sdx(2,1,2) = d2sdx(2,1,2) - xjaci22 * gphes(3,inode,igaus) * elcod(2,inode)
             d2sdx(1,2,2) = d2sdx(1,2,2) - xjaci11 * gphes(2,inode,igaus) * elcod(1,inode)
             d2sdx(1,2,2) = d2sdx(1,2,2) - xjaci12 * gphes(2,inode,igaus) * elcod(2,inode)
             d2sdx(2,2,2) = d2sdx(2,2,2) - xjaci21 * gphes(2,inode,igaus) * elcod(1,inode)
             d2sdx(2,2,2) = d2sdx(2,2,2) - xjaci22 * gphes(2,inode,igaus) * elcod(2,inode)
          end do
          !
          ! Computes the second Cartesian derivatives of the shape functions
          !
          do inode = 1,pnode

             gphes(1,inode,igaus) = gphes(1,inode,igaus) + deriv(1,inode,igaus) * d2sdx(1,1,1)
             gphes(1,inode,igaus) = gphes(1,inode,igaus) + deriv(2,inode,igaus) * d2sdx(2,1,1)

             gphes(3,inode,igaus) = gphes(3,inode,igaus) + deriv(1,inode,igaus) * d2sdx(1,1,2)
             gphes(3,inode,igaus) = gphes(3,inode,igaus) + deriv(2,inode,igaus) * d2sdx(2,1,2)

             gphes(2,inode,igaus) = gphes(2,inode,igaus) + deriv(1,inode,igaus) * d2sdx(1,2,2)
             gphes(2,inode,igaus) = gphes(2,inode,igaus) + deriv(2,inode,igaus) * d2sdx(2,2,2)

          end do

       end do

    else
       !
       ! 3D
       !
       do igaus  =  1,pgaus

          xjaci11  =  xjaci(1,1,igaus)
          xjaci12  =  xjaci(1,2,igaus)
          xjaci13  =  xjaci(1,3,igaus)
          xjaci21  =  xjaci(2,1,igaus)
          xjaci22  =  xjaci(2,2,igaus)
          xjaci23  =  xjaci(2,3,igaus)
          xjaci31  =  xjaci(3,1,igaus)
          xjaci32  =  xjaci(3,2,igaus)
          xjaci33  =  xjaci(3,3,igaus)

          do inode = 1,pnode

             gphes(1,inode,igaus) =  0.0_rp
             gphes(4,inode,igaus) =  0.0_rp
             gphes(5,inode,igaus) =  0.0_rp
             gphes(2,inode,igaus) =  0.0_rp
             gphes(6,inode,igaus) =  0.0_rp
             gphes(3,inode,igaus) =  0.0_rp

             heslo1               =  heslo(1,inode,igaus)
             heslo2               =  heslo(2,inode,igaus)
             heslo3               =  heslo(3,inode,igaus)
             heslo4               =  heslo(4,inode,igaus)
             heslo5               =  heslo(5,inode,igaus)
             heslo6               =  heslo(6,inode,igaus)

             gphes(1,inode,igaus) = gphes(1,inode,igaus) + xjaci11 * heslo1 * xjaci11
             gphes(1,inode,igaus) = gphes(1,inode,igaus) + xjaci11 * heslo4 * xjaci21
             gphes(1,inode,igaus) = gphes(1,inode,igaus) + xjaci11 * heslo5 * xjaci31
             gphes(1,inode,igaus) = gphes(1,inode,igaus) + xjaci21 * heslo4 * xjaci11
             gphes(1,inode,igaus) = gphes(1,inode,igaus) + xjaci21 * heslo2 * xjaci21
             gphes(1,inode,igaus) = gphes(1,inode,igaus) + xjaci21 * heslo6 * xjaci31
             gphes(1,inode,igaus) = gphes(1,inode,igaus) + xjaci31 * heslo5 * xjaci11
             gphes(1,inode,igaus) = gphes(1,inode,igaus) + xjaci31 * heslo6 * xjaci21
             gphes(1,inode,igaus) = gphes(1,inode,igaus) + xjaci31 * heslo3 * xjaci31
             gphes(4,inode,igaus) = gphes(4,inode,igaus) + xjaci11 * heslo1 * xjaci12
             gphes(4,inode,igaus) = gphes(4,inode,igaus) + xjaci11 * heslo4 * xjaci22
             gphes(4,inode,igaus) = gphes(4,inode,igaus) + xjaci11 * heslo5 * xjaci32
             gphes(4,inode,igaus) = gphes(4,inode,igaus) + xjaci21 * heslo4 * xjaci12
             gphes(4,inode,igaus) = gphes(4,inode,igaus) + xjaci21 * heslo2 * xjaci22
             gphes(4,inode,igaus) = gphes(4,inode,igaus) + xjaci21 * heslo6 * xjaci32
             gphes(4,inode,igaus) = gphes(4,inode,igaus) + xjaci31 * heslo5 * xjaci12
             gphes(4,inode,igaus) = gphes(4,inode,igaus) + xjaci31 * heslo6 * xjaci22
             gphes(4,inode,igaus) = gphes(4,inode,igaus) + xjaci31 * heslo3 * xjaci32
             gphes(5,inode,igaus) = gphes(5,inode,igaus) + xjaci11 * heslo1 * xjaci13
             gphes(5,inode,igaus) = gphes(5,inode,igaus) + xjaci11 * heslo4 * xjaci23
             gphes(5,inode,igaus) = gphes(5,inode,igaus) + xjaci11 * heslo5 * xjaci33
             gphes(5,inode,igaus) = gphes(5,inode,igaus) + xjaci21 * heslo4 * xjaci13
             gphes(5,inode,igaus) = gphes(5,inode,igaus) + xjaci21 * heslo2 * xjaci23
             gphes(5,inode,igaus) = gphes(5,inode,igaus) + xjaci21 * heslo6 * xjaci33
             gphes(5,inode,igaus) = gphes(5,inode,igaus) + xjaci31 * heslo5 * xjaci13
             gphes(5,inode,igaus) = gphes(5,inode,igaus) + xjaci31 * heslo6 * xjaci23
             gphes(5,inode,igaus) = gphes(5,inode,igaus) + xjaci31 * heslo3 * xjaci33

             gphes(2,inode,igaus) = gphes(2,inode,igaus) + xjaci12 * heslo1 * xjaci12
             gphes(2,inode,igaus) = gphes(2,inode,igaus) + xjaci12 * heslo4 * xjaci22
             gphes(2,inode,igaus) = gphes(2,inode,igaus) + xjaci12 * heslo5 * xjaci32
             gphes(2,inode,igaus) = gphes(2,inode,igaus) + xjaci22 * heslo4 * xjaci12
             gphes(2,inode,igaus) = gphes(2,inode,igaus) + xjaci22 * heslo2 * xjaci22
             gphes(2,inode,igaus) = gphes(2,inode,igaus) + xjaci22 * heslo6 * xjaci32
             gphes(2,inode,igaus) = gphes(2,inode,igaus) + xjaci32 * heslo5 * xjaci12
             gphes(2,inode,igaus) = gphes(2,inode,igaus) + xjaci32 * heslo6 * xjaci22
             gphes(2,inode,igaus) = gphes(2,inode,igaus) + xjaci32 * heslo3 * xjaci32
             gphes(6,inode,igaus) = gphes(6,inode,igaus) + xjaci12 * heslo1 * xjaci13
             gphes(6,inode,igaus) = gphes(6,inode,igaus) + xjaci12 * heslo4 * xjaci23
             gphes(6,inode,igaus) = gphes(6,inode,igaus) + xjaci12 * heslo5 * xjaci33
             gphes(6,inode,igaus) = gphes(6,inode,igaus) + xjaci22 * heslo4 * xjaci13
             gphes(6,inode,igaus) = gphes(6,inode,igaus) + xjaci22 * heslo2 * xjaci23
             gphes(6,inode,igaus) = gphes(6,inode,igaus) + xjaci22 * heslo6 * xjaci33
             gphes(6,inode,igaus) = gphes(6,inode,igaus) + xjaci32 * heslo5 * xjaci13
             gphes(6,inode,igaus) = gphes(6,inode,igaus) + xjaci32 * heslo6 * xjaci23
             gphes(6,inode,igaus) = gphes(6,inode,igaus) + xjaci32 * heslo3 * xjaci33

             gphes(3,inode,igaus) = gphes(3,inode,igaus) + xjaci13 * heslo1 * xjaci13
             gphes(3,inode,igaus) = gphes(3,inode,igaus) + xjaci13 * heslo4 * xjaci23
             gphes(3,inode,igaus) = gphes(3,inode,igaus) + xjaci13 * heslo5 * xjaci33
             gphes(3,inode,igaus) = gphes(3,inode,igaus) + xjaci23 * heslo4 * xjaci13
             gphes(3,inode,igaus) = gphes(3,inode,igaus) + xjaci23 * heslo2 * xjaci23
             gphes(3,inode,igaus) = gphes(3,inode,igaus) + xjaci23 * heslo6 * xjaci33
             gphes(3,inode,igaus) = gphes(3,inode,igaus) + xjaci33 * heslo5 * xjaci13
             gphes(3,inode,igaus) = gphes(3,inode,igaus) + xjaci33 * heslo6 * xjaci23
             gphes(3,inode,igaus) = gphes(3,inode,igaus) + xjaci33 * heslo3 * xjaci33

          end do

          d2sdx(1,1,1) = 0.0_rp
          d2sdx(2,1,1) = 0.0_rp
          d2sdx(3,1,1) = 0.0_rp
          d2sdx(1,1,2) = 0.0_rp
          d2sdx(2,1,2) = 0.0_rp
          d2sdx(3,1,2) = 0.0_rp
          d2sdx(1,2,2) = 0.0_rp
          d2sdx(2,2,2) = 0.0_rp
          d2sdx(3,2,2) = 0.0_rp
          d2sdx(1,1,3) = 0.0_rp
          d2sdx(2,1,3) = 0.0_rp
          d2sdx(3,1,3) = 0.0_rp
          d2sdx(1,2,3) = 0.0_rp
          d2sdx(2,2,3) = 0.0_rp
          d2sdx(3,2,3) = 0.0_rp
          d2sdx(1,3,3) = 0.0_rp
          d2sdx(2,3,3) = 0.0_rp
          d2sdx(3,3,3) = 0.0_rp

          do inode = 1,pnode
             d2sdx(1,1,1) = d2sdx(1,1,1) - xjaci11 * gphes(1,inode,igaus) * elcod(1,inode)
             d2sdx(1,1,1) = d2sdx(1,1,1) - xjaci12 * gphes(1,inode,igaus) * elcod(2,inode)
             d2sdx(1,1,1) = d2sdx(1,1,1) - xjaci13 * gphes(1,inode,igaus) * elcod(3,inode)
             d2sdx(2,1,1) = d2sdx(2,1,1) - xjaci21 * gphes(1,inode,igaus) * elcod(1,inode)
             d2sdx(2,1,1) = d2sdx(2,1,1) - xjaci22 * gphes(1,inode,igaus) * elcod(2,inode)
             d2sdx(2,1,1) = d2sdx(2,1,1) - xjaci23 * gphes(1,inode,igaus) * elcod(3,inode)
             d2sdx(3,1,1) = d2sdx(3,1,1) - xjaci31 * gphes(1,inode,igaus) * elcod(1,inode)
             d2sdx(3,1,1) = d2sdx(3,1,1) - xjaci32 * gphes(1,inode,igaus) * elcod(2,inode)
             d2sdx(3,1,1) = d2sdx(3,1,1) - xjaci33 * gphes(1,inode,igaus) * elcod(3,inode)
             d2sdx(1,1,2) = d2sdx(1,1,2) - xjaci11 * gphes(4,inode,igaus) * elcod(1,inode)
             d2sdx(1,1,2) = d2sdx(1,1,2) - xjaci12 * gphes(4,inode,igaus) * elcod(2,inode)
             d2sdx(1,1,2) = d2sdx(1,1,2) - xjaci13 * gphes(4,inode,igaus) * elcod(3,inode)
             d2sdx(2,1,2) = d2sdx(2,1,2) - xjaci21 * gphes(4,inode,igaus) * elcod(1,inode)
             d2sdx(2,1,2) = d2sdx(2,1,2) - xjaci22 * gphes(4,inode,igaus) * elcod(2,inode)
             d2sdx(2,1,2) = d2sdx(2,1,2) - xjaci23 * gphes(4,inode,igaus) * elcod(3,inode)
             d2sdx(3,1,2) = d2sdx(3,1,2) - xjaci31 * gphes(4,inode,igaus) * elcod(1,inode)
             d2sdx(3,1,2) = d2sdx(3,1,2) - xjaci32 * gphes(4,inode,igaus) * elcod(2,inode)
             d2sdx(3,1,2) = d2sdx(3,1,2) - xjaci33 * gphes(4,inode,igaus) * elcod(3,inode)
             d2sdx(1,2,2) = d2sdx(1,2,2) - xjaci11 * gphes(2,inode,igaus) * elcod(1,inode)
             d2sdx(1,2,2) = d2sdx(1,2,2) - xjaci12 * gphes(2,inode,igaus) * elcod(2,inode)
             d2sdx(1,2,2) = d2sdx(1,2,2) - xjaci13 * gphes(2,inode,igaus) * elcod(3,inode)
             d2sdx(2,2,2) = d2sdx(2,2,2) - xjaci21 * gphes(2,inode,igaus) * elcod(1,inode)
             d2sdx(2,2,2) = d2sdx(2,2,2) - xjaci22 * gphes(2,inode,igaus) * elcod(2,inode)
             d2sdx(2,2,2) = d2sdx(2,2,2) - xjaci23 * gphes(2,inode,igaus) * elcod(3,inode)
             d2sdx(3,2,2) = d2sdx(3,2,2) - xjaci31 * gphes(2,inode,igaus) * elcod(1,inode)
             d2sdx(3,2,2) = d2sdx(3,2,2) - xjaci32 * gphes(2,inode,igaus) * elcod(2,inode)
             d2sdx(3,2,2) = d2sdx(3,2,2) - xjaci33 * gphes(2,inode,igaus) * elcod(3,inode)
             d2sdx(1,1,3) = d2sdx(1,1,3) - xjaci11 * gphes(5,inode,igaus) * elcod(1,inode)
             d2sdx(1,1,3) = d2sdx(1,1,3) - xjaci12 * gphes(5,inode,igaus) * elcod(2,inode)
             d2sdx(1,1,3) = d2sdx(1,1,3) - xjaci13 * gphes(5,inode,igaus) * elcod(3,inode)
             d2sdx(2,1,3) = d2sdx(2,1,3) - xjaci21 * gphes(5,inode,igaus) * elcod(1,inode)
             d2sdx(2,1,3) = d2sdx(2,1,3) - xjaci22 * gphes(5,inode,igaus) * elcod(2,inode)
             d2sdx(2,1,3) = d2sdx(2,1,3) - xjaci23 * gphes(5,inode,igaus) * elcod(3,inode)
             d2sdx(3,1,3) = d2sdx(3,1,3) - xjaci31 * gphes(5,inode,igaus) * elcod(1,inode)
             d2sdx(3,1,3) = d2sdx(3,1,3) - xjaci32 * gphes(5,inode,igaus) * elcod(2,inode)
             d2sdx(3,1,3) = d2sdx(3,1,3) - xjaci33 * gphes(5,inode,igaus) * elcod(3,inode)
             d2sdx(1,2,3) = d2sdx(1,2,3) - xjaci11 * gphes(6,inode,igaus) * elcod(1,inode)
             d2sdx(1,2,3) = d2sdx(1,2,3) - xjaci12 * gphes(6,inode,igaus) * elcod(2,inode)
             d2sdx(1,2,3) = d2sdx(1,2,3) - xjaci13 * gphes(6,inode,igaus) * elcod(3,inode)
             d2sdx(2,2,3) = d2sdx(2,2,3) - xjaci21 * gphes(6,inode,igaus) * elcod(1,inode)
             d2sdx(2,2,3) = d2sdx(2,2,3) - xjaci22 * gphes(6,inode,igaus) * elcod(2,inode)
             d2sdx(2,2,3) = d2sdx(2,2,3) - xjaci23 * gphes(6,inode,igaus) * elcod(3,inode)
             d2sdx(3,2,3) = d2sdx(3,2,3) - xjaci31 * gphes(6,inode,igaus) * elcod(1,inode)
             d2sdx(3,2,3) = d2sdx(3,2,3) - xjaci32 * gphes(6,inode,igaus) * elcod(2,inode)
             d2sdx(3,2,3) = d2sdx(3,2,3) - xjaci33 * gphes(6,inode,igaus) * elcod(3,inode)
             d2sdx(1,3,3) = d2sdx(1,3,3) - xjaci11 * gphes(3,inode,igaus) * elcod(1,inode)
             d2sdx(1,3,3) = d2sdx(1,3,3) - xjaci12 * gphes(3,inode,igaus) * elcod(2,inode)
             d2sdx(1,3,3) = d2sdx(1,3,3) - xjaci13 * gphes(3,inode,igaus) * elcod(3,inode)
             d2sdx(2,3,3) = d2sdx(2,3,3) - xjaci21 * gphes(3,inode,igaus) * elcod(1,inode)
             d2sdx(2,3,3) = d2sdx(2,3,3) - xjaci22 * gphes(3,inode,igaus) * elcod(2,inode)
             d2sdx(2,3,3) = d2sdx(2,3,3) - xjaci23 * gphes(3,inode,igaus) * elcod(3,inode)
             d2sdx(3,3,3) = d2sdx(3,3,3) - xjaci31 * gphes(3,inode,igaus) * elcod(1,inode)
             d2sdx(3,3,3) = d2sdx(3,3,3) - xjaci32 * gphes(3,inode,igaus) * elcod(2,inode)
             d2sdx(3,3,3) = d2sdx(3,3,3) - xjaci33 * gphes(3,inode,igaus) * elcod(3,inode)
          end do
          !
          ! Computes the second Cartesian derivatives of the shape functions
          !
          do inode  =  1,pnode

             gphes(1,inode,igaus) = gphes(1,inode,igaus) + deriv(1,inode,igaus) * d2sdx(1,1,1)
             gphes(1,inode,igaus) = gphes(1,inode,igaus) + deriv(2,inode,igaus) * d2sdx(2,1,1)
             gphes(1,inode,igaus) = gphes(1,inode,igaus) + deriv(3,inode,igaus) * d2sdx(3,1,1)
             gphes(4,inode,igaus) = gphes(4,inode,igaus) + deriv(1,inode,igaus) * d2sdx(1,1,2)
             gphes(4,inode,igaus) = gphes(4,inode,igaus) + deriv(2,inode,igaus) * d2sdx(2,1,2)
             gphes(4,inode,igaus) = gphes(4,inode,igaus) + deriv(3,inode,igaus) * d2sdx(3,1,2)
             gphes(5,inode,igaus) = gphes(5,inode,igaus) + deriv(1,inode,igaus) * d2sdx(1,1,3)
             gphes(5,inode,igaus) = gphes(5,inode,igaus) + deriv(2,inode,igaus) * d2sdx(2,1,3)
             gphes(5,inode,igaus) = gphes(5,inode,igaus) + deriv(3,inode,igaus) * d2sdx(3,1,3)

             gphes(2,inode,igaus) = gphes(2,inode,igaus) + deriv(1,inode,igaus) * d2sdx(1,2,2)
             gphes(2,inode,igaus) = gphes(2,inode,igaus) + deriv(2,inode,igaus) * d2sdx(2,2,2)
             gphes(2,inode,igaus) = gphes(2,inode,igaus) + deriv(3,inode,igaus) * d2sdx(3,2,2)
             gphes(6,inode,igaus) = gphes(6,inode,igaus) + deriv(1,inode,igaus) * d2sdx(1,2,3)
             gphes(6,inode,igaus) = gphes(6,inode,igaus) + deriv(2,inode,igaus) * d2sdx(2,2,3)
             gphes(6,inode,igaus) = gphes(6,inode,igaus) + deriv(3,inode,igaus) * d2sdx(3,2,3)

             gphes(3,inode,igaus) = gphes(3,inode,igaus) + deriv(1,inode,igaus) * d2sdx(1,3,3)
             gphes(3,inode,igaus) = gphes(3,inode,igaus) + deriv(2,inode,igaus) * d2sdx(2,3,3)
             gphes(3,inode,igaus) = gphes(3,inode,igaus) + deriv(3,inode,igaus) * d2sdx(3,3,3)

          end do

       end do

    end if

  end subroutine elmgeo_hessian_scalar

  pure subroutine elmgeo_element_length_vector(&
       ndime,pnode,porde,tragl,hleng,elcod,elvel,chave,chale,&
       hnatu,kfl_advec,kfl_ellen)
    !-----------------------------------------------------------------------
    !****f* Domain/elmchl
    ! NAME
    !   elmchl
    ! DESCRIPTION
    !   This routine computes the characteristic element lengths CHALE
    !   according to a given strategy. CHALE is divided by two for
    !   quadratic elements:
    !   KFL_ELLEN = 0 ... CHALE(1) = Minimum element length
    !                 ... CHALE(2) = Minimum element length
    !   KFL_ELLEN = 1 ... CHALE(1) = Maximum element length
    !                 ... CHALE(2) = Maximum element length
    !   KFL_ELLEN = 2 ... CHALE(1) = Average element length
    !                 ... CHALE(2) = Average element length
    !   KFL_ELLEN = 3 ... IF KFL_ADVEC = 1:
    !                     CHALE(1) = Flow direction
    !                     CHALE(2) = Flow direction
    !                     ELSE IF KFL_ADVEC =0:
    !                     CHALE(1) = Minimum element length
    !                     CHALE(2) = Minimum element length
    !   KFL_ELLEN = 4 ... CHALE(1) = Approx. diameter=sqrt(hmin*hmax)
    !                 ... CHALE(2) = Approx. diameter=sqrt(hmin*hmax)
    !   KFL_ELLEN = 5 ... CHALE(1) = Length in flow direction
    !                 ... CHALE(2) = Minimum element kength
    ! OUTPUT
    !   CHALE
    ! USES
    ! USED BY
    !***
    !-----------------------------------------------------------------------

    integer(ip), intent(in)  :: ndime
    integer(ip), intent(in)  :: pnode
    integer(ip), intent(in)  :: porde
    integer(ip), intent(in)  :: kfl_advec
    integer(ip), intent(in)  :: kfl_ellen
    real(rp),    intent(in)  :: hnatu
    real(rp),    intent(out) :: chale(VECTOR_SIZE,2)
    real(rp),    intent(in)  :: tragl(VECTOR_SIZE,ndime,ndime)
    real(rp),    intent(in)  :: hleng(VECTOR_SIZE,ndime)
    real(rp),    intent(in)  :: elcod(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(in)  :: elvel(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(out) :: chave(VECTOR_SIZE,ndime,2)
    integer(ip)              :: idime,inode,ivect,kdime
    real(rp)                 :: elno1(VECTOR_SIZE)
    real(rp)                 :: elno2(VECTOR_SIZE)

    if( kfl_ellen == 0 ) then
       !
       ! Minimum element length
       !
       chale(1:VECTOR_SIZE,1) = hleng(1:VECTOR_SIZE,ndime)
       chale(1:VECTOR_SIZE,2) = chale(1:VECTOR_SIZE,1)

    else if( kfl_ellen == 1 ) then
       !
       ! Maximum element length
       !
       chale(1:VECTOR_SIZE,1) = hleng(1:VECTOR_SIZE,1)
       chale(1:VECTOR_SIZE,2) = chale(1:VECTOR_SIZE,1)

    else if( kfl_ellen == 2 ) then
       !
       ! Average length
       !
       chale(1:VECTOR_SIZE,1) = 0.0_rp
       do idime = 1,ndime
          chale(1:VECTOR_SIZE,1) = chale(1:VECTOR_SIZE,1) + hleng(1:VECTOR_SIZE,idime)
       end do
       chale(1:VECTOR_SIZE,1) = chale(1:VECTOR_SIZE,1) / real(ndime,rp)
       chale(1:VECTOR_SIZE,2) = chale(1:VECTOR_SIZE,1)

    else if( kfl_ellen == 3 ) then
       !
       ! Length in flow direction
       !
       return
       !call runend('ELMGEO_ELEMENT_LENGTH: NOT CODED')
       if( kfl_advec /= 0 ) then
          !
          ! Characteristic element velocity (average)
          !
          chave = 0.0_rp
          do idime = 1,ndime
             do inode = 1,pnode
                chave(1:VECTOR_SIZE,idime,1) = chave(1:VECTOR_SIZE,idime,1) + elvel(1:VECTOR_SIZE,idime,inode)
             end do
          end do
          do idime = 1,ndime
             chave(1:VECTOR_SIZE,idime,1) = chave(1:VECTOR_SIZE,idime,1) / real(pnode,rp)
          end do
          !
          ! Characteristic element length u^l = J^(-t) u^g
          !
          do idime = 1,ndime
             chave(1:VECTOR_SIZE,idime,2) = 0.0_rp
             do kdime = 1,ndime
                chave(1:VECTOR_SIZE,idime,2) = chave(1:VECTOR_SIZE,idime,2) + tragl(1:VECTOR_SIZE,idime,kdime) * chave(1:VECTOR_SIZE,kdime,1)
             end do
          end do
          elno1 = 0.0_rp
          elno2 = 0.0_rp
          do idime = 1,ndime
             elno1(1:VECTOR_SIZE) = elno1(1:VECTOR_SIZE) + chave(1:VECTOR_SIZE,idime,1) * chave(1:VECTOR_SIZE,idime,1)
             elno2(1:VECTOR_SIZE) = elno2(1:VECTOR_SIZE) + chave(1:VECTOR_SIZE,idime,2) * chave(1:VECTOR_SIZE,idime,2)
          end do
          elno1(1:VECTOR_SIZE)  = sqrt(elno1(1:VECTOR_SIZE))
          elno2(1:VECTOR_SIZE)  = sqrt(elno2(1:VECTOR_SIZE))

          do ivect = 1,VECTOR_SIZE
             if( elno2(ivect) > 1.0e-16_rp.and. elno1(ivect) > 1.0e-16_rp) then
                chale(ivect,1) = hnatu * elno1(ivect) / elno2(ivect)
             else
                chale(ivect,1) = hleng(ivect,ndime)
             end if
          end do

          if( ndime == 3 ) then
             chale(1:VECTOR_SIZE,2) = (hleng(1:VECTOR_SIZE,ndime)*hleng(1:VECTOR_SIZE,2)*hleng(1:VECTOR_SIZE,1))**(1.0_rp/3.0_rp)
          else if (ndime==2) then
             chale(1:VECTOR_SIZE,2) = sqrt(hleng(1:VECTOR_SIZE,2)*hleng(1:VECTOR_SIZE,1))
          end if

       else
          chale(1:VECTOR_SIZE,1) = hleng(1:VECTOR_SIZE,ndime)
          chale(1:VECTOR_SIZE,2) = hleng(1:VECTOR_SIZE,ndime)
       end if

    else if( kfl_ellen == 4 ) then
       !
       ! sqrt(hmin*hmax)
       !
       chale(1:VECTOR_SIZE,1) = sqrt(hleng(1:VECTOR_SIZE,1)*hleng(1:VECTOR_SIZE,ndime))
       chale(1:VECTOR_SIZE,2) = chale(1:VECTOR_SIZE,1)

    else if( kfl_ellen == 5 ) then
       !
       ! Along velocity direction
       !
!!$       call velchl(pnode,elcod,elvel,chale,hleng)

    else if( kfl_ellen == 6 ) then
       !
       ! Mixed element length - hmin for tau1, hmax for tau2 - here we only obtain the values for tau1 - tau2 directly in nsi_elmsgs
       !
       chale(1:VECTOR_SIZE,1) = hleng(1:VECTOR_SIZE,ndime)
       chale(1:VECTOR_SIZE,2) = chale(1:VECTOR_SIZE,1)

    end if
    !
    ! Divide h by 2 for quadratic elements and 3 for cubic elements
    !
    chale(1:VECTOR_SIZE,1) = chale(1:VECTOR_SIZE,1) / real(porde,rp)
    chale(1:VECTOR_SIZE,2) = chale(1:VECTOR_SIZE,2) / real(porde,rp)

  end subroutine elmgeo_element_length_vector

  subroutine elmgeo_element_length_scalar(&
       ndime,pnode,porde,tragl,hleng,elcod,elvel,chave,chale,&
       hnatu,kfl_advec,kfl_ellen)
    !-----------------------------------------------------------------------
    !****f* Domain/elmchl
    ! NAME
    !   elmchl
    ! DESCRIPTION
    !   This routine computes the characteristic element lengths CHALE
    !   according to a given strategy. CHALE is divided by two for
    !   quadratic elements:
    !   KFL_ELLEN = 0 ... CHALE(1) = Minimum element length
    !                 ... CHALE(2) = Minimum element length
    !   KFL_ELLEN = 1 ... CHALE(1) = Maximum element length
    !                 ... CHALE(2) = Maximum element length
    !   KFL_ELLEN = 2 ... CHALE(1) = Average element length
    !                 ... CHALE(2) = Average element length
    !   KFL_ELLEN = 3 ... IF KFL_ADVEC = 1:
    !                     CHALE(1) = Flow direction
    !                     CHALE(2) = Flow direction
    !                     ELSE IF KFL_ADVEC =0:
    !                     CHALE(1) = Minimum element length
    !                     CHALE(2) = Minimum element length
    !   KFL_ELLEN = 4 ... CHALE(1) = Approx. diameter=sqrt(hmin*hmax)
    !                 ... CHALE(2) = Approx. diameter=sqrt(hmin*hmax)
    !   KFL_ELLEN = 5 ... CHALE(1) = Length in flow direction
    !                 ... CHALE(2) = Minimum element kength
    ! OUTPUT
    !   CHALE
    ! USES
    ! USED BY
    !***
    !-----------------------------------------------------------------------

    integer(ip), intent(in)  :: ndime
    integer(ip), intent(in)  :: pnode
    integer(ip), intent(in)  :: porde
    integer(ip), intent(in)  :: kfl_advec
    integer(ip), intent(in)  :: kfl_ellen
    real(rp),    intent(in)  :: hnatu
    real(rp),    intent(out) :: chale(2)
    real(rp),    intent(in)  :: tragl(ndime,ndime)
    real(rp),    intent(in)  :: hleng(ndime)
    real(rp),    intent(in)  :: elcod(ndime,pnode)
    real(rp),    intent(in)  :: elvel(ndime,pnode)
    real(rp),    intent(out) :: chave(ndime,2)
    integer(ip)              :: idime,inode,kdime
    real(rp)                 :: elno1
    real(rp)                 :: elno2

    if( kfl_ellen == 0 ) then
       !
       ! Minimum element length
       !
       chale(1) = hleng(ndime)
       chale(2) = chale(1)

    else if( kfl_ellen == 1 ) then
       !
       ! Maximum element length
       !
       chale(1) = hleng(1)
       chale(2) = chale(1)

    else if( kfl_ellen == 2 ) then
       !
       ! Average length
       !
       chale(1) = 0.0_rp
       do idime = 1,ndime
          chale(1) = chale(1) + hleng(idime)
       end do
       chale(1) = chale(1) / real(ndime,rp)
       chale(2) = chale(1)

    else if( kfl_ellen == 3 ) then
       !
       ! Length in flow direction
       !
       call runend('ELMGEO_ELEMENT_LENGTH: NOT CODED')
       if( kfl_advec /= 0 ) then
          !
          ! Characteristic element velocity (average)
          !
          chave = 0.0_rp
          do idime = 1,ndime
             do inode = 1,pnode
                chave(idime,1) = chave(idime,1) + elvel(idime,inode)
             end do
          end do
          do idime = 1,ndime
             chave(idime,1) = chave(idime,1) / real(pnode,rp)
          end do
          !
          ! Characteristic element length u^l = J^(-t) u^g
          !
          do idime = 1,ndime
             chave(idime,2) = 0.0_rp
             do kdime = 1,ndime
                chave(idime,2) = chave(idime,2) + tragl(idime,kdime) * chave(kdime,1)
             end do
          end do
          elno1 = 0.0_rp
          elno2 = 0.0_rp
          do idime = 1,ndime
             elno1 = elno1 + chave(idime,1) * chave(idime,1)
             elno2 = elno2 + chave(idime,2) * chave(idime,2)
          end do
          elno1  = sqrt(elno1)
          elno2  = sqrt(elno2)

          if( elno2 > 1.0e-16_rp.and. elno1 > 1.0e-16_rp) then
             chale(1) = hnatu * elno1 / elno2
          else
             chale(1) = hleng(ndime)
          end if

          if( ndime == 3 ) then
             chale(2) = (hleng(ndime)*hleng(2)*hleng(1))**(1.0_rp/3.0_rp)
          else if (ndime==2) then
             chale(2) = sqrt(hleng(2)*hleng(1))
          end if

       else
          chale(1) = hleng(ndime)
          chale(2) = hleng(ndime)
       end if

    else if( kfl_ellen == 4 ) then
       !
       ! sqrt(hmin*hmax)
       !
       chale(1) = sqrt(hleng(1)*hleng(ndime))
       chale(2) = chale(1)

    else if( kfl_ellen == 5 ) then
       !
       ! Along velocity direction
       !
!!$       call velchl(pnode,elcod,elvel,chale,hleng)

    else if( kfl_ellen == 6 ) then
       !
       ! Mixed element length - hmin for tau1, hmax for tau2 - here we only obtain the values for tau1 - tau2 directly in nsi_elmsgs
       !
       chale(1) = hleng(ndime)
       chale(2) = chale(1)

    end if
    !
    ! Divide h by 2 for quadratic elements and 3 for cubic elements
    !
    chale(1) = chale(1) / real(porde,rp)
    chale(2) = chale(2) / real(porde,rp)

  end subroutine elmgeo_element_length_scalar

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    12/01/2018
  !> @brief   Guess element type
  !> @details Find the element type as a function of number of nodes
  !>          and dimension
  !>
  !-----------------------------------------------------------------------

  function elmgeo_element_type(ndime,pnode)

    integer(ip), intent(in) :: ndime
    integer(ip), intent(in) :: pnode
    integer(ip)             :: elmgeo_element_type

    elmgeo_element_type = 0

    if( ndime == 3 ) then
       if(pnode == 3) then
          elmgeo_element_type = SHELL
       else if(pnode == 4) then
          elmgeo_element_type = TET04
       else if(pnode == 5) then
          elmgeo_element_type = PYR05
       else if(pnode == 6) then
          elmgeo_element_type = PEN06
       else if(pnode == 8) then
          elmgeo_element_type = HEX08
       else if(pnode == 10) then
          elmgeo_element_type = TET10
       else if(pnode == 14) then
          elmgeo_element_type = PYR14
       else if(pnode == 15) then
          elmgeo_element_type = PEN15
       else if(pnode == 18) then
          elmgeo_element_type = PEN18
       else if(pnode == 20) then
          !elmgeo_element_type = HEX20
          elmgeo_element_type = TET20
       else if(pnode == 27) then
          elmgeo_element_type = HEX27
       else if(pnode == 64) then
          elmgeo_element_type = HEX64
       end if
    else if(ndime == 2) then
       if(pnode == 3) then
          elmgeo_element_type = TRI03
       else if(pnode == 4) then
          elmgeo_element_type = QUA04
       else if(pnode == 6) then
          elmgeo_element_type = TRI06
       else if(pnode == 8) then
          elmgeo_element_type = QUA08
       else if(pnode == 9) then
          elmgeo_element_type = QUA09
       else if(pnode == 10) then
          elmgeo_element_type = TRI10
       else if(pnode == 16) then
          elmgeo_element_type = QUA16
       end if
    else if( ndime == 1 ) then
       if(pnode == 2) then
          elmgeo_element_type = BAR02
       else if(pnode == 3) then
          elmgeo_element_type = BAR03
       else if(pnode == 4) then
          elmgeo_element_type = BAR04
       end if
    else if( ndime == 0 ) then
       elmgeo_element_type=POINT
    end if

  end function elmgeo_element_type


  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    12/01/2018
  !> @brief   Coordinates of iso-parametric element
  !> @details Coordinates of iso-parametric element
  !>
  !-----------------------------------------------------------------------

    subroutine elmgeo_iso_parametric_element_coordinates(ndime,ielty,iso_element_coord,NODE)

    integer(ip), intent(in)           :: ndime
    integer(ip), intent(in)           :: ielty
    real(rp),    intent(out)          :: iso_element_coord(ndime,*)
    integer(ip), intent(in), optional :: NODE
    integer(ip)                       :: pnode
    real(rp),    allocatable          :: iso_element_coord_loc(:,:)
    
    pnode = element_type(ielty) % number_nodes
    allocate( iso_element_coord_loc(ndime,pnode) )

    select case ( ielty )
       
    case ( BAR02 )

       iso_element_coord_loc(1:1,1) = (/ -1.0_rp /)
       iso_element_coord_loc(1:1,2) = (/  1.0_rp /)

    case ( BAR03 )

       iso_element_coord_loc(1:1,1) = (/ -1.0_rp /)
       iso_element_coord_loc(1:1,2) = (/  1.0_rp /)
       iso_element_coord_loc(1:1,3) = (/  0.0_rp /)

    case ( BAR04 )

       iso_element_coord_loc(1:1,1) = (/ -1.0_rp        /)
       iso_element_coord_loc(1:1,2) = (/  1.0_rp        /)
       iso_element_coord_loc(1:1,3) = (/ -1.0_rp/3.0_rp /)
       iso_element_coord_loc(1:1,4) = (/  1.0_rp/3.0_rp /)

    case ( TRI03 )

       iso_element_coord_loc(1:2,1) = (/  0.0_rp,  0.0_rp /)
       iso_element_coord_loc(1:2,2) = (/  1.0_rp,  0.0_rp /)
       iso_element_coord_loc(1:2,3) = (/  0.0_rp,  1.0_rp /)

    case ( TRI06 )

       iso_element_coord_loc(1:2,1) = (/  0.0_rp,  0.0_rp /)
       iso_element_coord_loc(1:2,2) = (/  1.0_rp,  0.0_rp /)
       iso_element_coord_loc(1:2,3) = (/  0.0_rp,  1.0_rp /)
       iso_element_coord_loc(1:2,4) = (/  0.5_rp,  0.0_rp /)
       iso_element_coord_loc(1:2,5) = (/  0.5_rp,  0.5_rp /)
       iso_element_coord_loc(1:2,6) = (/  0.0_rp,  0.5_rp /)

    case ( TRI10 )

       iso_element_coord_loc(1:2, 1) = (/  0.0_rp,  0.0_rp /)
       iso_element_coord_loc(1:2, 2) = (/  1.0_rp,  0.0_rp /)
       iso_element_coord_loc(1:2, 3) = (/  0.0_rp,  1.0_rp /)
       iso_element_coord_loc(1:2, 4) = (/  1.0_rp/3.0_rp,  0.0_rp /)
       iso_element_coord_loc(1:2, 5) = (/  2.0_rp/3.0_rp,  0.0_rp /)
       iso_element_coord_loc(1:2, 6) = (/  2.0_rp/3.0_rp,  1.0_rp/3.0_rp /)
       iso_element_coord_loc(1:2, 7) = (/  1.0_rp/3.0_rp,  2.0_rp/3.0_rp /)
       iso_element_coord_loc(1:2, 8) = (/  0.0_rp,  2.0_rp/3.0_rp /)
       iso_element_coord_loc(1:2, 9) = (/  0.0_rp,  1.0_rp/3.0_rp /)
       iso_element_coord_loc(1:2,10) = (/ 1.0_rp/3.0_rp,  1.0_rp/3.0_rp /)

    case ( QUA04 )

       iso_element_coord_loc(1:2,1) = (/ -1.0_rp, -1.0_rp /)
       iso_element_coord_loc(1:2,2) = (/  1.0_rp, -1.0_rp /)
       iso_element_coord_loc(1:2,3) = (/  1.0_rp,  1.0_rp /)
       iso_element_coord_loc(1:2,4) = (/ -1.0_rp,  1.0_rp /)

    case ( QUA09 )

       iso_element_coord_loc(1:2,1) = (/ -1.0_rp, -1.0_rp /)
       iso_element_coord_loc(1:2,2) = (/  1.0_rp, -1.0_rp /)
       iso_element_coord_loc(1:2,3) = (/  1.0_rp,  1.0_rp /)
       iso_element_coord_loc(1:2,4) = (/ -1.0_rp,  1.0_rp /)
       iso_element_coord_loc(1:2,5) = (/  0.0_rp, -1.0_rp /)
       iso_element_coord_loc(1:2,6) = (/  1.0_rp,  0.0_rp /)
       iso_element_coord_loc(1:2,7) = (/  0.0_rp,  1.0_rp /)
       iso_element_coord_loc(1:2,8) = (/ -1.0_rp,  0.0_rp /)
       iso_element_coord_loc(1:2,9) = (/  0.0_rp,  0.0_rp /)

    case ( QUA16 )
       
       iso_element_coord_loc(1:2, 1) = (/ -1.0_rp, -1.0_rp /)
       iso_element_coord_loc(1:2, 2) = (/  1.0_rp, -1.0_rp /)
       iso_element_coord_loc(1:2, 3) = (/  1.0_rp,  1.0_rp /)
       iso_element_coord_loc(1:2, 4) = (/ -1.0_rp,  1.0_rp /)
       iso_element_coord_loc(1:2, 5) = (/   -one3, -1.0_rp /)
       iso_element_coord_loc(1:2, 6) = (/    one3, -1.0_rp /)
       iso_element_coord_loc(1:2, 7) = (/  1.0_rp,   -one3 /)
       iso_element_coord_loc(1:2, 8) = (/  1.0_rp,    one3 /)
       iso_element_coord_loc(1:2, 9) = (/    one3,  1.0_rp /)
       iso_element_coord_loc(1:2,10) = (/   -one3,  1.0_rp /)
       iso_element_coord_loc(1:2,11) = (/ -1.0_rp,    one3 /)
       iso_element_coord_loc(1:2,12) = (/ -1.0_rp,   -one3 /)
       iso_element_coord_loc(1:2,13) = (/   -one3,   -one3 /)
       iso_element_coord_loc(1:2,14) = (/    one3,   -one3 /)
       iso_element_coord_loc(1:2,15) = (/    one3,    one3 /)
       iso_element_coord_loc(1:2,16) = (/   -one3,    one3 /)

    case ( TET04 )

       iso_element_coord_loc(1:3,1) = (/  0.0_rp,  0.0_rp,  0.0_rp /)
       iso_element_coord_loc(1:3,2) = (/  1.0_rp,  0.0_rp,  0.0_rp /)
       iso_element_coord_loc(1:3,3) = (/  0.0_rp,  1.0_rp,  0.0_rp /)
       iso_element_coord_loc(1:3,4) = (/  0.0_rp,  0.0_rp,  1.0_rp /)

    case ( TET10 )

       iso_element_coord_loc(1:3, 1) = (/  0.0_rp,  0.0_rp,  0.0_rp /)
       iso_element_coord_loc(1:3, 2) = (/  1.0_rp,  0.0_rp,  0.0_rp /)
       iso_element_coord_loc(1:3, 3) = (/  0.0_rp,  1.0_rp,  0.0_rp /)
       iso_element_coord_loc(1:3, 4) = (/  0.0_rp,  0.0_rp,  1.0_rp /)
       iso_element_coord_loc(1:3, 5) = (/  0.5_rp,  0.0_rp,  0.0_rp /)
       iso_element_coord_loc(1:3, 6) = (/  0.5_rp,  0.5_rp,  0.0_rp /)
       iso_element_coord_loc(1:3, 7) = (/  0.0_rp,  0.5_rp,  0.0_rp /)
       iso_element_coord_loc(1:3, 8) = (/  0.0_rp,  0.0_rp,  0.5_rp /)
       iso_element_coord_loc(1:3, 9) = (/  0.5_rp,  0.0_rp,  0.5_rp /)
       iso_element_coord_loc(1:3,10) = (/  0.0_rp,  0.5_rp,  0.5_rp /)

    case ( TET20 )

       iso_element_coord_loc(1:3, 1) = (/  0.0_rp,         0.0_rp,         0.0_rp        /)
       iso_element_coord_loc(1:3, 2) = (/  1.0_rp,         0.0_rp,         0.0_rp        /)
       iso_element_coord_loc(1:3, 3) = (/  0.0_rp,         1.0_rp,         0.0_rp        /)
       iso_element_coord_loc(1:3, 4) = (/  0.0_rp,         0.0_rp,         1.0_rp        /)
       iso_element_coord_loc(1:3, 5) = (/  1.0_rp/3.0_rp,  0.0_rp,         0.0_rp        /)
       iso_element_coord_loc(1:3, 6) = (/  2.0_rp/3.0_rp,  0.0_rp,         0.0_rp        /)
       iso_element_coord_loc(1:3, 7) = (/  2.0_rp/3.0_rp,  1.0_rp/3.0_rp,  0.0_rp        /)
       iso_element_coord_loc(1:3, 8) = (/  1.0_rp/3.0_rp,  2.0_rp/3.0_rp,  0.0_rp        /)
       iso_element_coord_loc(1:3, 9) = (/  0.0_rp,         2.0_rp/3.0_rp,  0.0_rp        /)
       iso_element_coord_loc(1:3,10) = (/  0.0_rp,         1.0_rp/3.0_rp,  0.0_rp        /)
       iso_element_coord_loc(1:3,11) = (/  0.0_rp,         0.0_rp,         2.0_rp/3.0_rp /)
       iso_element_coord_loc(1:3,12) = (/  0.0_rp,         0.0_rp,         1.0_rp/3.0_rp /)
       iso_element_coord_loc(1:3,13) = (/  0.0_rp,         1.0_rp/3.0_rp,  2.0_rp/3.0_rp /)
       iso_element_coord_loc(1:3,14) = (/  0.0_rp,         2.0_rp/3.0_rp,  1.0_rp/3.0_rp /)
       iso_element_coord_loc(1:3,15) = (/  1.0_rp/3.0_rp,  0.0_rp,         2.0_rp/3.0_rp /)
       iso_element_coord_loc(1:3,16) = (/  2.0_rp/3.0_rp,  0.0_rp,         1.0_rp/3.0_rp /)
       iso_element_coord_loc(1:3,17) = (/  1.0_rp/3.0_rp,  1.0_rp/3.0_rp,  0.0_rp        /)
       iso_element_coord_loc(1:3,18) = (/  1.0_rp/3.0_rp,  0.0_rp,         1.0_rp/3.0_rp /)
       iso_element_coord_loc(1:3,19) = (/  0.0_rp,         1.0_rp/3.0_rp,  1.0_rp/3.0_rp /)
       iso_element_coord_loc(1:3,20) = (/  1.0_rp/3.0_rp,  1.0_rp/3.0_rp,  1.0_rp/3.0_rp /)

    case ( PYR05 )

       iso_element_coord_loc(1:3,1) = (/ -1.0_rp, -1.0_rp, -1.0_rp /)
       iso_element_coord_loc(1:3,2) = (/  1.0_rp, -1.0_rp, -1.0_rp /)
       iso_element_coord_loc(1:3,3) = (/  1.0_rp,  1.0_rp, -1.0_rp /)
       iso_element_coord_loc(1:3,4) = (/ -1.0_rp,  1.0_rp, -1.0_rp /)
       iso_element_coord_loc(1:3,5) = (/  0.0_rp,  0.0_rp,  1.0_rp /)

    case ( PEN06 )

       iso_element_coord_loc(1:3,1) = (/  0.0_rp,  0.0_rp,  0.0_rp /)
       iso_element_coord_loc(1:3,2) = (/  1.0_rp,  0.0_rp,  0.0_rp /)
       iso_element_coord_loc(1:3,3) = (/  0.0_rp,  1.0_rp,  0.0_rp /)
       iso_element_coord_loc(1:3,4) = (/  0.0_rp,  0.0_rp,  1.0_rp /)
       iso_element_coord_loc(1:3,5) = (/  1.0_rp,  0.0_rp,  1.0_rp /)
       iso_element_coord_loc(1:3,6) = (/  0.0_rp,  1.0_rp,  1.0_rp /)

    case ( PEN18 )

       iso_element_coord_loc(1:3, 1) = (/  0.0_rp,  0.0_rp,  0.0_rp /)
       iso_element_coord_loc(1:3, 2) = (/  1.0_rp,  0.0_rp,  0.0_rp /)
       iso_element_coord_loc(1:3, 3) = (/  0.0_rp,  1.0_rp,  0.0_rp /)
       iso_element_coord_loc(1:3, 4) = (/  0.0_rp,  0.0_rp,  1.0_rp /)
       iso_element_coord_loc(1:3, 5) = (/  1.0_rp,  0.0_rp,  1.0_rp /)
       iso_element_coord_loc(1:3, 6) = (/  0.0_rp,  1.0_rp,  1.0_rp /)
       iso_element_coord_loc(1:3, 7) = (/  0.5_rp,  0.0_rp,  0.0_rp /)
       iso_element_coord_loc(1:3, 8) = (/  0.5_rp,  0.5_rp,  0.0_rp /)
       iso_element_coord_loc(1:3, 9) = (/  0.0_rp,  0.5_rp,  0.0_rp /)
       iso_element_coord_loc(1:3,10) = (/  0.0_rp,  0.0_rp,  0.5_rp /)
       iso_element_coord_loc(1:3,11) = (/  1.0_rp,  0.0_rp,  0.5_rp /)
       iso_element_coord_loc(1:3,12) = (/  0.0_rp,  1.0_rp,  0.5_rp /)
       iso_element_coord_loc(1:3,13) = (/  0.5_rp,  0.0_rp,  1.0_rp /)
       iso_element_coord_loc(1:3,14) = (/  0.5_rp,  0.5_rp,  1.0_rp /)
       iso_element_coord_loc(1:3,15) = (/  0.0_rp,  0.5_rp,  1.0_rp /)
       iso_element_coord_loc(1:3,16) = (/  0.5_rp,  0.0_rp,  0.5_rp /)
       iso_element_coord_loc(1:3,17) = (/  0.5_rp,  0.5_rp,  0.5_rp /)
       iso_element_coord_loc(1:3,18) = (/  0.0_rp,  0.5_rp,  0.5_rp /)

    case ( HEX08 )

       iso_element_coord_loc(1:3,1) = (/ -1.0_rp, -1.0_rp, -1.0_rp /)
       iso_element_coord_loc(1:3,2) = (/  1.0_rp, -1.0_rp, -1.0_rp /)
       iso_element_coord_loc(1:3,3) = (/  1.0_rp,  1.0_rp, -1.0_rp /)
       iso_element_coord_loc(1:3,4) = (/ -1.0_rp,  1.0_rp, -1.0_rp /)
       iso_element_coord_loc(1:3,5) = (/ -1.0_rp, -1.0_rp,  1.0_rp /)
       iso_element_coord_loc(1:3,6) = (/  1.0_rp, -1.0_rp,  1.0_rp /)
       iso_element_coord_loc(1:3,7) = (/  1.0_rp,  1.0_rp,  1.0_rp /)
       iso_element_coord_loc(1:3,8) = (/ -1.0_rp,  1.0_rp,  1.0_rp /)

    case ( HEX27 )

       iso_element_coord_loc(1:3, 1) = (/ -1.0_rp , -1.0_rp , -1.0_rp /)     
       iso_element_coord_loc(1:3, 2) = (/  1.0_rp , -1.0_rp , -1.0_rp /)     
       iso_element_coord_loc(1:3, 3) = (/  1.0_rp ,  1.0_rp , -1.0_rp /)     
       iso_element_coord_loc(1:3, 4) = (/ -1.0_rp ,  1.0_rp , -1.0_rp /)     
       iso_element_coord_loc(1:3, 5) = (/ -1.0_rp , -1.0_rp ,  1.0_rp /)     
       iso_element_coord_loc(1:3, 6) = (/  1.0_rp , -1.0_rp ,  1.0_rp /)     
       iso_element_coord_loc(1:3, 7) = (/  1.0_rp ,  1.0_rp ,  1.0_rp /)     
       iso_element_coord_loc(1:3, 8) = (/ -1.0_rp ,  1.0_rp ,  1.0_rp /)     
       iso_element_coord_loc(1:3, 9) = (/  0.0_rp , -1.0_rp , -1.0_rp /)     
       iso_element_coord_loc(1:3,10) = (/  1.0_rp ,  0.0_rp , -1.0_rp /)     
       iso_element_coord_loc(1:3,11) = (/  0.0_rp ,  1.0_rp , -1.0_rp /)     
       iso_element_coord_loc(1:3,12) = (/ -1.0_rp ,  0.0_rp , -1.0_rp /)     
       iso_element_coord_loc(1:3,13) = (/ -1.0_rp , -1.0_rp ,  0.0_rp /)     
       iso_element_coord_loc(1:3,14) = (/  1.0_rp , -1.0_rp ,  0.0_rp /)     
       iso_element_coord_loc(1:3,15) = (/  1.0_rp ,  1.0_rp ,  0.0_rp /)     
       iso_element_coord_loc(1:3,16) = (/ -1.0_rp ,  1.0_rp ,  0.0_rp /)     
       iso_element_coord_loc(1:3,17) = (/  0.0_rp , -1.0_rp ,  1.0_rp /)     
       iso_element_coord_loc(1:3,18) = (/  1.0_rp ,  0.0_rp ,  1.0_rp /)     
       iso_element_coord_loc(1:3,19) = (/  0.0_rp ,  1.0_rp ,  1.0_rp /)     
       iso_element_coord_loc(1:3,20) = (/ -1.0_rp ,  0.0_rp ,  1.0_rp /)     
       iso_element_coord_loc(1:3,21) = (/  0.0_rp ,  0.0_rp , -1.0_rp /)     
       iso_element_coord_loc(1:3,22) = (/  0.0_rp , -1.0_rp ,  0.0_rp /)     
       iso_element_coord_loc(1:3,23) = (/  1.0_rp ,  0.0_rp ,  0.0_rp /)     
       iso_element_coord_loc(1:3,24) = (/  0.0_rp ,  1.0_rp ,  0.0_rp /)     
       iso_element_coord_loc(1:3,25) = (/ -1.0_rp ,  0.0_rp ,  0.0_rp /)    
       iso_element_coord_loc(1:3,26) = (/  0.0_rp ,  0.0_rp ,  1.0_rp /)    
       iso_element_coord_loc(1:3,27) = (/  0.0_rp ,  0.0_rp ,  0.0_rp /)    

    case ( HEX64 ) 

       iso_element_coord_loc(1:3, 1) = (/ -1.0_rp, -1.0_rp, -1.0_rp /)
       iso_element_coord_loc(1:3, 2) = (/  1.0_rp, -1.0_rp, -1.0_rp /)    
       iso_element_coord_loc(1:3, 3) = (/  1.0_rp,  1.0_rp, -1.0_rp /)    
       iso_element_coord_loc(1:3, 4) = (/ -1.0_rp,  1.0_rp, -1.0_rp /)    
       iso_element_coord_loc(1:3, 5) = (/ -1.0_rp, -1.0_rp,  1.0_rp /)    
       iso_element_coord_loc(1:3, 6) = (/  1.0_rp, -1.0_rp,  1.0_rp /)    
       iso_element_coord_loc(1:3, 7) = (/  1.0_rp,  1.0_rp,  1.0_rp /)    
       iso_element_coord_loc(1:3, 8) = (/ -1.0_rp,  1.0_rp,  1.0_rp /)    
       iso_element_coord_loc(1:3, 9) = (/   -one3, -1.0_rp, -1.0_rp /)    
       iso_element_coord_loc(1:3,10) = (/    one3, -1.0_rp, -1.0_rp /)    
       iso_element_coord_loc(1:3,11) = (/  1.0_rp,   -one3, -1.0_rp /)    
       iso_element_coord_loc(1:3,12) = (/  1.0_rp,    one3, -1.0_rp /)    
       iso_element_coord_loc(1:3,13) = (/    one3,  1.0_rp, -1.0_rp /)
       iso_element_coord_loc(1:3,14) = (/   -one3,  1.0_rp, -1.0_rp /)    
       iso_element_coord_loc(1:3,15) = (/ -1.0_rp,    one3, -1.0_rp /)    
       iso_element_coord_loc(1:3,16) = (/ -1.0_rp,   -one3, -1.0_rp /)    
       iso_element_coord_loc(1:3,17) = (/ -1.0_rp, -1.0_rp,   -one3 /)      
       iso_element_coord_loc(1:3,18) = (/  1.0_rp, -1.0_rp,   -one3 /)    
       iso_element_coord_loc(1:3,19) = (/  1.0_rp,  1.0_rp,   -one3 /)    
       iso_element_coord_loc(1:3,20) = (/ -1.0_rp,  1.0_rp,   -one3 /)    
       iso_element_coord_loc(1:3,21) = (/ -1.0_rp, -1.0_rp,    one3 /)
       iso_element_coord_loc(1:3,22) = (/  1.0_rp, -1.0_rp,    one3 /)    
       iso_element_coord_loc(1:3,23) = (/  1.0_rp,  1.0_rp,    one3 /) 
       iso_element_coord_loc(1:3,24) = (/ -1.0_rp,  1.0_rp,    one3 /)    
       iso_element_coord_loc(1:3,25) = (/   -one3, -1.0_rp,  1.0_rp /)    
       iso_element_coord_loc(1:3,26) = (/    one3, -1.0_rp,  1.0_rp /)    
       iso_element_coord_loc(1:3,27) = (/  1.0_rp,   -one3,  1.0_rp /)    
       iso_element_coord_loc(1:3,28) = (/  1.0_rp,    one3,  1.0_rp /)    
       iso_element_coord_loc(1:3,29) = (/    one3,  1.0_rp,  1.0_rp /)    
       iso_element_coord_loc(1:3,30) = (/   -one3,  1.0_rp,  1.0_rp /)    
       iso_element_coord_loc(1:3,31) = (/ -1.0_rp,    one3,  1.0_rp /)    
       iso_element_coord_loc(1:3,32) = (/ -1.0_rp,   -one3,  1.0_rp /)    
       iso_element_coord_loc(1:3,33) = (/   -one3,   -one3, -1.0_rp /)   
       iso_element_coord_loc(1:3,34) = (/    one3,   -one3, -1.0_rp /)    
       iso_element_coord_loc(1:3,35) = (/    one3,    one3, -1.0_rp /)    
       iso_element_coord_loc(1:3,36) = (/   -one3,    one3, -1.0_rp /)    
       iso_element_coord_loc(1:3,37) = (/   -one3, -1.0_rp,   -one3 /)    
       iso_element_coord_loc(1:3,38) = (/    one3, -1.0_rp,   -one3 /)    
       iso_element_coord_loc(1:3,39) = (/  1.0_rp,   -one3,   -one3 /)    
       iso_element_coord_loc(1:3,40) = (/  1.0_rp,    one3,   -one3 /)    
       iso_element_coord_loc(1:3,41) = (/    one3,  1.0_rp,   -one3 /)    
       iso_element_coord_loc(1:3,42) = (/   -one3,  1.0_rp,   -one3 /)    
       iso_element_coord_loc(1:3,43) = (/ -1.0_rp,    one3,   -one3 /)      
       iso_element_coord_loc(1:3,44) = (/ -1.0_rp,   -one3,   -one3 /)
       iso_element_coord_loc(1:3,45) = (/   -one3, -1.0_rp,    one3 /)    
       iso_element_coord_loc(1:3,46) = (/    one3, -1.0_rp,    one3 /)    
       iso_element_coord_loc(1:3,47) = (/  1.0_rp,   -one3,    one3 /)    
       iso_element_coord_loc(1:3,48) = (/  1.0_rp,    one3,    one3 /)    
       iso_element_coord_loc(1:3,49) = (/    one3,  1.0_rp,    one3 /)    
       iso_element_coord_loc(1:3,50) = (/   -one3,  1.0_rp,    one3 /)           
       iso_element_coord_loc(1:3,51) = (/ -1.0_rp,    one3,    one3 /)    
       iso_element_coord_loc(1:3,52) = (/ -1.0_rp,   -one3,    one3 /)
       iso_element_coord_loc(1:3,53) = (/   -one3,   -one3,  1.0_rp /)    
       iso_element_coord_loc(1:3,54) = (/    one3,   -one3,  1.0_rp /)    
       iso_element_coord_loc(1:3,55) = (/    one3,    one3,  1.0_rp /)    
       iso_element_coord_loc(1:3,56) = (/   -one3,    one3,  1.0_rp /)    
       iso_element_coord_loc(1:3,57) = (/   -one3,   -one3,   -one3 /)
       iso_element_coord_loc(1:3,58) = (/    one3,   -one3,   -one3 /)    
       iso_element_coord_loc(1:3,59) = (/    one3,    one3,   -one3 /)    
       iso_element_coord_loc(1:3,60) = (/   -one3,    one3,   -one3 /)    
       iso_element_coord_loc(1:3,61) = (/   -one3,   -one3,    one3 /)    
       iso_element_coord_loc(1:3,62) = (/    one3,   -one3,    one3 /)    
       iso_element_coord_loc(1:3,63) = (/    one3,    one3,    one3 /)    
       iso_element_coord_loc(1:3,64) = (/   -one3,    one3,    one3 /)    

    case default
       
       !call runend('UNKNON ELEMENT: CANNOT COMPUTE ISO-PARAMETRIC COORDINATES')

    end select

    if( present(NODE) ) then
       if( NODE > pnode .or. NODE < 1 ) call runend('WONR NUMBER OF NODE FOR ISO PARAMETRIC ELEMENT')
       iso_element_coord(1:ndime,1)       = iso_element_coord_loc(1:ndime,NODE)
    else
       iso_element_coord(1:ndime,1:pnode) = iso_element_coord_loc(1:ndime,1:pnode)
    end if

    deallocate(iso_element_coord_loc)
       
  end subroutine elmgeo_iso_parametric_element_coordinates

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    18/10/2018
  !> @brief   Check if normal is outwards 
  !> @details Check if normal is outwards
  !>
  !>          \verbatim
  !>
  !>          o----------o
  !>          |          |    cge = element center of gravity
  !>          |          |    cgb = boundary center of gravity
  !>          |    cge   |    v   = vector(cgb,cge)
  !>          |    /\    |
  !>          |  v ||    |
  !>          o----cgb---o <- iboun
  !>               ||
  !>               \/ n
  !>
  !>          \endverbatim
  !>
  !>          The procedure is the following:
  !>          - Compute element center of gravity COCOG
  !>         - Compute boundary center of gravity COCOB
  !>          - v = COCOG - COCOB
  !>          - Compute n.v
  !>          - Invert BALOC if necessary
  !>
  !-----------------------------------------------------------------------

  subroutine elmgeo_check_boundary_orientation(ndime,pnode,pnodb,baloc,bocod,elcod,ninve)
    
    integer(ip), intent(in)    :: ndime
    integer(ip), intent(in)    :: pnode,pnodb
    real(rp),    intent(in)    :: bocod(ndime,pnodb),elcod(ndime,pnode)
    real(rp),    intent(inout) :: baloc(ndime)
    integer(ip), intent(inout) :: ninve
    integer(ip)                :: inode,inodb
    real(rp)                   :: produ,cocog(3),cocob(3),dummr,dummb
    !
    ! Coordinates center of gravity
    !
    dummr = 1.0_rp / real(pnode,rp)
    dummb = 1.0_rp / real(pnodb,rp)

    if( ndime == 1 ) then

       !-------------------------------------------------------------------
       !
       ! 1D
       !
       !-------------------------------------------------------------------

       cocog(1)=0.0_rp 
       do inode=1,pnode
          cocog(1)=cocog(1)+elcod(1,inode)
       end do
       cocog(1)=cocog(1)*dummr

       produ = baloc(1)*baloc(1)
       if( produ /= 0.0_rp ) baloc(1)=baloc(1)/sqrt(produ)

       produ=(cocog(1)-bocod(1,1))*baloc(1)

       if( produ > 0.0_rp ) then
          baloc(1) = -baloc(1)         ! n=-n                      
       end if

    else if( ndime == 2 ) then

       !-------------------------------------------------------------------
       !
       ! 2D
       !
       !-------------------------------------------------------------------

       cocog(1)=0.0_rp
       cocog(2)=0.0_rp
       do inode=1,pnode
          cocog(1)=cocog(1)+elcod(1,inode)
          cocog(2)=cocog(2)+elcod(2,inode)
       end do
       cocog(1)=cocog(1)*dummr
       cocog(2)=cocog(2)*dummr

       cocob(1)=0.0_rp
       cocob(2)=0.0_rp
       do inodb=1,pnodb
          cocob(1)=cocob(1)+bocod(1,inodb)
          cocob(2)=cocob(2)+bocod(2,inodb)
       end do
       cocob(1)=cocob(1)*dummb
       cocob(2)=cocob(2)*dummb

       produ= baloc(1)*baloc(1) &
            + baloc(2)*baloc(2)
       if( produ /= 0.0_rp ) then
          produ    = 1.0_rp/sqrt(produ)
          baloc(1) = produ*baloc(1)
          baloc(2) = produ*baloc(2)
       end if

       produ=(cocog(1)-cocob(1))*baloc(1)&
            +(cocog(2)-cocob(2))*baloc(2)

       if(produ>0.0_rp) then
          baloc(1) = -baloc(1)     ! n =-n                    
          baloc(2) = -baloc(2)     ! n =-n
       end if

    else

       !-------------------------------------------------------------------
       !
       ! 3D
       !
       !-------------------------------------------------------------------

       cocog(1) = 0.0_rp
       cocog(2) = 0.0_rp
       cocog(3) = 0.0_rp
       do inode = 1,pnode
          cocog(1) = cocog(1) + elcod(1,inode)
          cocog(2) = cocog(2) + elcod(2,inode)
          cocog(3) = cocog(3) + elcod(3,inode)
       end do
       cocog(1) = cocog(1) * dummr
       cocog(2) = cocog(2) * dummr
       cocog(3) = cocog(3) * dummr

       cocob(1) = 0.0_rp
       cocob(2) = 0.0_rp
       cocob(3) = 0.0_rp
       do inodb = 1,pnodb
          cocob(1) = cocob(1) + bocod(1,inodb)
          cocob(2) = cocob(2) + bocod(2,inodb)
          cocob(3) = cocob(3) + bocod(3,inodb)
       end do
       cocob(1) = cocob(1) * dummb
       cocob(2) = cocob(2) * dummb
       cocob(3) = cocob(3) * dummb

       produ = baloc(1)*baloc(1) &
            +  baloc(2)*baloc(2) &
            +  baloc(3)*baloc(3)

       if( produ /= 0.0_rp ) then
          produ    = 1.0_rp/sqrt(produ)
          baloc(1) = produ*baloc(1)
          baloc(2) = produ*baloc(2)
          baloc(3) = produ*baloc(3)
       end if

       produ =   ( cocog(1) - cocob(1) ) * baloc(1) &
            &  + ( cocog(2) - cocob(2) ) * baloc(2) &
            &  + ( cocog(3) - cocob(3) ) * baloc(3)

       if( produ > 0.0_rp ) then                   
          ninve    =  ninve + 1
          baloc(1) = -baloc(1)     ! n =-n                  
          baloc(2) = -baloc(2)     ! n =-n                  
          baloc(3) = -baloc(3)     ! n =-n
       end if

    end if

  end subroutine elmgeo_check_boundary_orientation

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    18/10/2018
  !> @brief   This routine computes the boundary normals
  !> @details This routine computes the boundary normals
  !>
  !-----------------------------------------------------------------------

  subroutine elmgeo_triangle_normal(p1,p2,p3,coord_in,vec,ndime)

    integer(ip), intent(in)  :: p1,p2,p3,ndime
    real(rp),    intent(in)  :: coord_in(ndime,*)
    real(rp),    intent(out) :: vec(3,3)

    vec(1,1) = coord_in(1,p2) - coord_in(1,p1)
    vec(2,1) = coord_in(2,p2) - coord_in(2,p1)
    vec(3,1) = coord_in(3,p2) - coord_in(3,p1)
    vec(1,2) = coord_in(1,p3) - coord_in(1,p1)
    vec(2,2) = coord_in(2,p3) - coord_in(2,p1)
    vec(3,2) = coord_in(3,p3) - coord_in(3,p1)
    
    call maths_vectorial_product(vec(1,1),vec(1,2),vec(1,3),ndime)

  end subroutine elmgeo_triangle_normal

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    18/10/2018
  !> @brief   Test element jacobian
  !> @details This subroutine creates a mesh in GiD format to visualize
  !>          the values of the Jacobian in an alement. Output can be
  !>          given in local or global coordinates using variable
  !>          local_coord
  !>
  !-----------------------------------------------------------------------

  subroutine elmgeo_output_jacobian(&
       ndime,pnode,elcod)

    integer(ip),  intent(in)            :: ndime
    integer(ip),  intent(in)            :: pnode
    real(rp),     intent(in)            :: elcod(ndime,pnode)
    real(rp)                            :: shapf(pnode)
    real(rp)                            :: deriv(ndime,pnode)
    integer(ip)                         :: inode,jdime,maxit,idime
    integer(ip)                         :: jnode,ielty,ii,jj,kk,ll
    real(rp)                            :: coloc(3),coglo(3)
    real(rp)                            :: xjacm(ndime,ndime)
    real(rp)                            :: xjaci(ndime,ndime)
    real(rp)                            :: detja
    real(rp)                            :: xmin,delta
    real(rp),   allocatable             :: iso_element(:,:)
    integer(ip)                         :: lumsh,lures,ielem
    logical(lg)                         :: local_coord

    xmin        = -2.0_rp
    delta       = -xmin*2.0_rp
    local_coord = .false.
    !
    ! Iterate for ds^i+1 = - (df/ds)^-1 . f(s^i)
    !
    lumsh = 90
    lures = 91
    open(unit=lumsh,file='mesh.post.msh',status='unknown')
    open(unit=lures,file='mesh.post.res',status='unknown')
    
    ielty = elmgeo_element_type(ndime,pnode)
    allocate(iso_element(ndime,pnode))
    call elmgeo_iso_parametric_element_coordinates(ndime,ielty,iso_element)

    maxit = 100
    
    ll    = 0

    write(lumsh,1)&
         'JACOBIAN',ndime,'Hexahedra',8_ip
    
    write(lumsh,2) 'coordinates'
    do kk = 1,maxit
       do jj = 1,maxit
          do ii = 1,maxit
             ll       = ll + 1
             coloc(1) = xmin + delta*real(ii-1,rp)/real(maxit-1,rp)
             coloc(2) = xmin + delta*real(jj-1,rp)/real(maxit-1,rp)
             coloc(3) = xmin + delta*real(kk-1,rp)/real(maxit-1,rp)
             if( local_coord ) then
                write(lumsh,3) ll,coloc(1:ndime)
             else
                call elmgeo_shapf_deriv_heslo(ndime,pnode,coloc,shapf,deriv) 
                coglo(1) = dot_product(shapf(1:pnode),elcod(1,1:pnode))
                coglo(2) = dot_product(shapf(1:pnode),elcod(2,1:pnode))
                coglo(3) = dot_product(shapf(1:pnode),elcod(3,1:pnode))
                write(lumsh,3) ll,coglo(1:ndime)
             end if
          end do
       end do
    end do
    
    jnode = ll
    if( local_coord ) then
       do ii = 1,pnode
          ll = ll + 1
          write(lumsh,3) ll,iso_element(1:ndime,ii)
       end do
    else
       do ii = 1,pnode
          ll = ll + 1
          write(lumsh,3) ll,elcod(1:ndime,ii)
       end do       
    end if
    write(lumsh,2) 'end coordinates'

    ll    = 0
    ielem = 0
    write(lumsh,2) 'elements'
    do kk = 1,maxit-1
       do jj = 1,maxit-1
          do ii = 1,maxit-1
             ll = ll + 1
             ielem = ielem + 1
             write(lumsh,4) ielem,&
                  &         ll,            ll+1,            ll+maxit+1,ll+maxit,&
                  &         ll+maxit*maxit,ll+maxit*maxit+1,ll+maxit*maxit+maxit+1,ll+maxit*maxit+maxit
          end do
          ll = ll + 1
       end do
       ll = ll + maxit
    end do
    write(lumsh,2) 'end elements'

    write(lumsh,1)&
         'ELEMENT_'//element_type(ielty) % name,ndime,&
         adjustl(trim(element_type(ielty) % nametopo)),&
         element_type(ielty) % number_nodes
    write(lumsh,2) 'elements'
    write(lumsh,4) ielem+1,(jnode+inode,inode=1,pnode)
    write(lumsh,2) 'end elements'
    
    write(lures,'(a)') 'GiD Post Results File 1.0'
    write(lures,'(a)') ' '

    write(lures,'(a)') 'Result JACOBIAN ALYA 0.0 Scalar OnNodes'
    write(lures,'(a)') 'ComponentNames JACOBIAN'
    write(lures,'(a)') 'values'
    ll    = 0
    do kk = 1,maxit ; do jj = 1,maxit ; do ii = 1,maxit
       ll       = ll + 1
       coloc(1) = xmin + delta*real(ii-1,rp)/real(maxit-1,rp)
       coloc(2) = xmin + delta*real(jj-1,rp)/real(maxit-1,rp)
       coloc(3) = xmin + delta*real(kk-1,rp)/real(maxit-1,rp)             
       call elmgeo_shapf_deriv_heslo(ndime,pnode,coloc,shapf,deriv)
       xjacm = 0.0_rp
       do inode = 1,pnode
          do jdime = 1,ndime
             do idime = 1,ndime
                xjacm(idime,jdime) = xjacm(idime,jdime) + deriv(jdime,inode) * elcod(idime,inode)
             end do
          end do
       end do
       call maths_invert_matrix(ndime,xjacm,detja,xjaci)             
       write(lures,*) ll,detja              
    end do; end do; end do    
    write(lures,'(a)') 'end values'

!!$    write(lures,'(a)') 'Result COORD ALYA 0.0 Vector OnNodes'
!!$    write(lures,'(a)') 'ComponentNames COORD_X,COORD_Y,COORD_Z'
!!$    write(lures,'(a)') 'values'
!!$    ll    = 0
!!$    do kk = 1,maxit ; do jj = 1,maxit ; do ii = 1,maxit
!!$       ll       = ll + 1
!!$       coloc(1) = xmin + delta*real(ii-1,rp)/real(maxit-1,rp)
!!$       coloc(2) = xmin + delta*real(jj-1,rp)/real(maxit-1,rp)
!!$       coloc(3) = xmin + delta*real(kk-1,rp)/real(maxit-1,rp)             
!!$       call elmgeo_shapf_deriv_heslo(ndime,pnode,coloc,shapf,deriv)
!!$       xjacm = 0.0_rp
!!$       do inode = 1,pnode
!!$          do jdime = 1,ndime
!!$             do idime = 1,ndime
!!$                xjacm(idime,jdime) = xjacm(idime,jdime) + deriv(jdime,inode) * elcod(idime,inode)
!!$             end do
!!$          end do
!!$       end do
!!$       call maths_invert_matrix(ndime,xjacm,detja,xjaci)             
!!$       write(lures,*) ll,coloc(1:ndime)              
!!$    end do; end do; end do    
!!$    write(lures,'(a)') 'end values'

1   format('MESH ',a,' dimension ',i1,' Elemtype ',a,' Nnode ',i2)
2   format(a)
3   format(i9, 3(1x,e16.8e3))
4   format(i9,50(1x,i9))
    
  end subroutine elmgeo_output_jacobian

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-26
  !> @brief   Nearest element node to a point
  !> @details Nearest element node to a point
  !> 
  !-----------------------------------------------------------------------

  subroutine elmgeo_nearest_element_node( &
       ndime,pelty,pnode,elcod,coglo,     &
       inode_nearest,coloc,shapf,deriv)
    
    integer(ip), intent(in)            :: ndime              !< Problem dimension
    integer(ip), intent(in)            :: pelty              !< Element type
    integer(ip), intent(in)            :: pnode              !< Element number of nodes
    real(rp),    intent(in)            :: elcod(ndime,pnode) !< Element coordinates
    real(rp),    intent(in)            :: coglo(ndime)       !< Test point coordinates
    integer(ip), intent(out)           :: inode_nearest      !< Nearest node
    real(rp),    intent(out), optional :: coloc(3)           !< Local coordinates of test point in element
    real(rp),    intent(out), optional :: shapf(pnode)       !< Shape functions of test point
    real(rp),    intent(out), optional :: deriv(ndime,pnode) !< Shape funcitons derivatives of test point

    integer(ip)                        :: inode
    real(rp)                           :: dista_nearest,dista

    dista_nearest = huge(1.0_rp)
    do inode = 1,pnode
       dista = dot_product(coglo(1:ndime)-elcod(1:ndime,inode),coglo(1:ndime)-elcod(1:ndime,inode))
       if( dista <= dista_nearest ) then
          inode_nearest = inode
          dista_nearest = dista
       end if
    end do
    dista_nearest = sqrt(dista_nearest)

    if( present(coloc) ) then
       call elmgeo_iso_parametric_element_coordinates(ndime,pelty,coloc,inode_nearest)
       if( present(shapf) ) call elmgeo_shapf_deriv_heslo(ndime,pnode,coloc,shapf,deriv)
    end if
    
  end subroutine elmgeo_nearest_element_node

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-25
  !> @brief   Name to type
  !> @details Give the element type given its name
  !> 
  !-----------------------------------------------------------------------

  subroutine elmgeo_element_name_to_type(welty,ielty)

    character(5), intent(in)  :: welty
    integer(ip),  intent(out) :: ielty

    if(      welty == 'BAR02' ) then
       ielty = BAR02
    else if( welty == 'BAR03' )  then
       ielty = BAR03   
    else if( welty == 'BAR04' )  then
       ielty = BAR04   
    !else if( welty == 'SB105' )  then
    !   ielty = SB105
    else if( welty == 'TRI03' )  then
       ielty = TRI03   
    else if( welty == 'TRI06' )  then
       ielty = TRI06   
    else if( welty == 'TRI10' )  then
       ielty = TRI10
    else if( welty == 'QUA04' )  then
       ielty = QUA04   
    else if( welty == 'QUA08' )  then
       ielty = QUA08   
    else if( welty == 'QUA09' )  then
       ielty = QUA09   
    else if( welty == 'QUA16' ) then
       ielty = QUA16  
    !else if( welty == 'SQ203' ) then
    !   ielty = SQ203  
    !else if( welty == 'SQ205' ) then
    !   ielty = SQ205  
    else if( welty == 'TET04' )  then
       ielty = TET04  
    else if( welty == 'TET10' ) then
       ielty = TET10 
    else if( welty == 'TET20' ) then
       ielty = TET20
    else if( welty == 'PYR05' )  then
       ielty = PYR05   
    else if( welty == 'PYR14' ) then
       ielty = PYR14 
    else if( welty == 'PEN06' )  then
       ielty = PEN06   
    else if( welty == 'PEN15' ) then
       ielty = PEN15  
    else if( welty == 'PEN18' ) then
       ielty = PEN18  
    else if( welty == 'HEX08' )  then
       ielty = HEX08   

    !else if( welty == 'HEX20' ) then
    !   ielty = HEX20  

    else if( welty == 'HEX27' ) then
       ielty = HEX27    
    else if( welty == 'HEX64' ) then
       ielty = HEX64
    !else if( welty == 'SH303' ) then
    !   ielty = SH303 
    !else if( welty == 'SH305' ) then
    !   ielty = SH305 
    else if( welty == 'SHELL' ) then
       ielty = SHELL
    else if( welty == 'BAR3D' ) then
       ielty = BAR3D
    else if( welty == 'POI3D' ) then
       ielty = POI3D
    else
       ielty = 0
    end if

  end subroutine elmgeo_element_name_to_type

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-05
  !> @brief   Circumradius
  !> @details Circumradius of an element
  !>  
  !-----------------------------------------------------------------------
 
  real(rp) function elmgeo_circumradius(ielty,xx) result(R)

    integer(ip), intent(in) :: ielty
    real(rp),    intent(in) :: xx(:,:)
    integer(ip)             :: inode,jnode,iedge
    
    if( ielty == TRI03 ) then
       R = maths_circumradius_triangle(xx(1:2,1),xx(1:2,2),xx(1:2,3))
    else if( ielty == TET04 ) then
       R = maths_circumradius_tetrahedron(xx(1:3,1),xx(1:3,2),xx(1:3,3),xx(1:3,4))
    else if( ielty == QUA04 ) then
       R = maths_circumradius_quadrilateral(xx(1:2,1),xx(1:2,2),xx(1:2,3),xx(1:2,4))
    else
       R  = 0.0_rp
       do iedge = 1,element_type(ielty) % number_edges
          inode = element_type(ielty) % list_edges(1,iedge)
          jnode = element_type(ielty) % list_edges(2,iedge)
          R     = R + sqrt(dot_product( xx(:,inode)-xx(:,jnode),xx(:,inode)-xx(:,jnode) ))
       end do
       R = R / real(element_type(ielty) % number_edges,rp)
    end if
       
  end function elmgeo_circumradius

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-19
  !> @brief   Compute a value in an element given by a Gauss point value
  !> @details Uses extrapolation function
  !> 
  !-----------------------------------------------------------------------

  pure subroutine elmgeo_gauss_to_element_s(pnode,pgaus,shapf,shaga,xx,yy,INITIALIZE)

    integer(ip),           intent(in)    :: pnode
    integer(ip),           intent(in)    :: pgaus
    real(rp),              intent(in)    :: shaga(pgaus,pnode)
    real(rp),              intent(in)    :: shapf(pnode)
    real(rp),              intent(in)    :: xx(pgaus)
    real(rp),              intent(inout) :: yy
    logical(lg), optional, intent(in)    :: INITIALIZE
    integer(ip)                          :: igaus,inode
    real(rp)                             :: zz(pnode)
    logical(lg)                          :: if_initialize

    zz = 0.0_rp
    do inode = 1,pnode
       do igaus = 1,pgaus
          zz(inode) = zz(inode) + shaga(igaus,inode) * xx(igaus)
       end do
    end do

    if( present(INITIALIZE) ) then
       if_initialize = INITIALIZE
    else       
       if_initialize = .true.
    end if
    if( if_initialize ) yy = 0.0_rp
    
    do inode = 1,pnode
       yy = yy + zz(inode) * shapf(inode)
    end do
    
  end subroutine elmgeo_gauss_to_element_s
  
  pure subroutine elmgeo_gauss_to_element_1(pnode,pgaus,shapf,shaga,ndofn,xx,yy,INITIALIZE)

    integer(ip),           intent(in)    :: pnode
    integer(ip),           intent(in)    :: pgaus
    real(rp),              intent(in)    :: shaga(pgaus,pnode)
    real(rp),              intent(in)    :: shapf(pnode)
    integer(ip),           intent(in)    :: ndofn
    real(rp),              intent(in)    :: xx(ndofn,pgaus)
    real(rp),              intent(inout) :: yy(ndofn)
    logical(lg), optional, intent(in)    :: INITIALIZE
    integer(ip)                          :: igaus,inode
    real(rp)                             :: zz(ndofn,pnode)
    logical(lg)                          :: if_initialize

    zz = 0.0_rp
    do inode = 1,pnode
       do igaus = 1,pgaus
          zz(:,inode) = zz(:,inode) + shaga(igaus,inode) * xx(:,igaus)
       end do
    end do

    if( present(INITIALIZE) ) then
       if_initialize = INITIALIZE
    else       
       if_initialize = .true.
    end if
    if( if_initialize ) yy = 0.0_rp
    
    do inode = 1,pnode
       yy(:) = yy(:) + zz(:,inode) * shapf(inode)
    end do
    
  end subroutine elmgeo_gauss_to_element_1

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-19
  !> @brief   Number of nodes
  !> @details Return number of nodes
  !> 
  !-----------------------------------------------------------------------
  
  integer(ip) pure function elmgeo_number_nodes(pelty,lnods)
    
    integer(ip),           intent(in) :: pelty
    integer(ip), optional, intent(in) :: lnods(:)
    integer(ip)                       :: qelty

    qelty = abs(pelty) 
    if(  qelty <= element_end ) then
       elmgeo_number_nodes = element_type(qelty) % number_nodes       
    else if( present(lnods) ) then
       elmgeo_number_nodes = maths_maxloc_nonzero(lnods(:))
    else
       elmgeo_number_nodes = 0
    end if
    
  end function elmgeo_number_nodes

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-19
  !> @brief   Number of subelements
  !> @details Number of subelements generated by mesh multiplication
  !> 
  !-----------------------------------------------------------------------
  
  function elmgeo_mm_elements(pelty) result(howdi)

    integer(ip), intent(in) :: pelty
    integer(ip)             :: howdi

    if(      pelty < 0 ) then
       howdi = 1   ! Negative type
    else if( pelty == POINT  ) then
       howdi = 1   ! POINT
    else if( pelty == PYR05  ) then
       howdi = 10  ! PYR05
    else if( pelty == BAR3D  ) then
       howdi = 2   ! BAR3D
    else if( pelty == SHELL  ) then
       howdi = 4   ! SHELL
    else if( pelty >= element_num_ini(1) .and. pelty <= element_num_end(1) ) then
       howdi = 2   ! 1D elements
    else if( pelty >= element_num_ini(2) .and. pelty <= element_num_end(2) ) then
       howdi = 4   ! 2D elements
    else if( pelty >= element_num_ini(3) .and. pelty <= element_num_end(3) ) then
       howdi = 8   ! 3D elements
    else
       howdi = 0
    end if

  end function elmgeo_mm_elements
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-19
  !> @brief   Number of interior nodes
  !> @details Number of interior nodes generated by mesh multiplication
  !> 
  !-----------------------------------------------------------------------
  
  function elmgeo_mm_interior_nodes(pelty) result(newno)

    integer(ip), intent(in) :: pelty
    integer(ip)             :: newno

    select case ( pelty )
       
    case ( BAR02 ) ; newno =  1
    case ( BAR03 ) ; newno =  2
    case ( BAR04 ) ; newno =  3
    case ( TRI03 ) ; newno =  0
    case ( TRI06 ) ; newno =  3
    case ( QUA04 ) ; newno =  1
    case ( QUA09 ) ; newno =  8
    case ( QUA16 ) ; newno = 21
    case ( TET10 ) ; newno =  1
    case ( HEX08 ) ; newno =  1
    case ( HEX27 ) ; newno = 26
    case default   ; newno =  0
       
    end select

  end function elmgeo_mm_interior_nodes
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-19
  !> @brief   Number of interior nodes
  !> @details Number of interior nodes generated by mesh multiplication
  !> 
  !-----------------------------------------------------------------------
  
  function elmgeo_mm_face_nodes(pelty) result(newfa)

    integer(ip), intent(in) :: pelty
    integer(ip)             :: newfa

    select case ( pelty )
       
    case ( TRI03 ) ; newfa =  1
    case ( TRI06 ) ; newfa =  3
    case ( QUA04 ) ; newfa =  1
    case ( QUA09 ) ; newfa =  8
    case ( QUA16 ) ; newfa = 21
    case default   ; newfa =  0
       
    end select

  end function elmgeo_mm_face_nodes
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-19
  !> @brief   Number of subelements
  !> @details Number of subelements generated by mesh multiplication
  !> 
  !-----------------------------------------------------------------------

  subroutine elmgeo_order_boundary_nodes(pelty,lnods,lnodb,lboel,istat)

    integer(ip),           intent(in)    :: pelty
    integer(ip),           intent(in)    :: lnods(:)
    integer(ip),           intent(inout) :: lnodb(:)
    integer(ip),           intent(out)   :: lboel(:)
    integer(ip), optional, intent(out)   :: istat
    integer(ip)                          :: pos,num,iface,inodb,pnodb,ipoin

    loop_faces: do iface = 1,element_type(pelty) % number_faces
       pnodb = element_type(pelty) % node_faces(iface)
       if( pnodb <= size(lnodb) ) then
          num = 0
          do inodb = 1,pnodb
             ipoin = lnods(element_type(pelty) % list_faces(inodb,iface))
             pos   = elmgeo_findloc(lnodb,ipoin)
             num   = num + min(1_ip,pos)
          end do
          if( num == pnodb ) then
             do inodb = 1,pnodb
                ipoin        = lnods(element_type(pelty) % list_faces(inodb,iface))
                lnodb(inodb) = ipoin
                lboel(inodb) = elmgeo_findloc(lnods,ipoin)
             end do
             if( present(istat) ) istat = 0
             return
          end if
       end if
    end do loop_faces

    if( present(istat) ) istat = 1
    
  end subroutine elmgeo_order_boundary_nodes
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-19
  !> @brief   Find a face
  !> @details Find the face number corresponding to a boundary
  !>          connectivity
  !> 
  !-----------------------------------------------------------------------

  function elmgeo_boundary_face(pelty,lnods,lnodb,istat) result(iface)

    integer(ip),           intent(in)    :: pelty
    integer(ip),           intent(in)    :: lnods(:)
    integer(ip),           intent(inout) :: lnodb(:)
    integer(ip), optional, intent(out)   :: istat
    integer(ip)                          :: pos,num,iface,inodb,pnodb,ipoin

    iface = 0
    loop_faces: do iface = 1,element_type(pelty) % number_faces
       pnodb = element_type(pelty) % node_faces(iface)
       if( pnodb <= size(lnodb) ) then
          num = 0
          do inodb = 1,pnodb
             ipoin = lnods(element_type(pelty) % list_faces(inodb,iface))
             pos   = elmgeo_findloc(lnodb,ipoin)
             num   = num + min(1_ip,pos)
          end do
          if( num == pnodb ) then
             return
          end if
       end if
    end do loop_faces

    if( present(istat) ) istat = 1
    
  end function elmgeo_boundary_face
  
  function elmgeo_findloc(xx,ii) result(pos)
    integer(ip), intent(in) :: xx(:)
    integer(ip), intent(in) :: ii
    integer(ip)             :: pos,kk

    kk = 0
    pos = 0
    do while( kk < size(xx) )
       kk = kk + 1
       if( xx(kk) == ii ) then
          pos = kk
          return
       end if       
    end do

  end function elmgeo_findloc

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-19
  !> @brief   Find an edge
  !> @details Find the edge number corresponding to an edge
  !>          connectivity
  !> 
  !-----------------------------------------------------------------------

  function elmgeo_edge(edges,inod1,inod2,istat) result(iedge)

    integer(ip),           intent(in)    :: edges(:,:)
    integer(ip),           intent(in)    :: inod1
    integer(ip),           intent(in)    :: inod2
    integer(ip), optional, intent(out)   :: istat
    integer(ip)                          :: iedge,ii

    if( present(istat) ) istat = 0
    do ii = 1,size(edges,2)
       if(    ( inod1 == edges(1,ii) .and. inod2 == edges(2,ii) ) .or. &
            & ( inod1 == edges(2,ii) .and. inod2 == edges(1,ii) ) ) then
          iedge = ii
          return
       end if
    end do
    
    if( present(istat) ) istat = 1
    
  end function elmgeo_edge  
 
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-19
  !> @brief   Find a face
  !> @details Find the face number corresponding to a face
  !>          connectivity
  !> 
  !-----------------------------------------------------------------------

  function elmgeo_face(pelty,nodes,istat) result(iface)

    integer(ip),           intent(in)    :: pelty
    integer(ip),           intent(in)    :: nodes(:)
    integer(ip), optional, intent(out)   :: istat
    integer(ip)                          :: iface,inodf,inode,knode

    if( present(istat) ) istat = 0    
    
    do iface = 1,element_type(pelty) % number_faces
       knode = 0
       loop_inodf: do inodf = 1,element_type(pelty) % node_faces(iface)
          inode = element_type(pelty) % list_faces(inodf,iface)
          if( any(nodes==inode) ) then
             knode = knode + 1
             exit loop_inodf
          end if
       end do loop_inodf
       if( knode == size(nodes) ) return
    end do
    
    if( present(istat) ) istat = 1
    
  end function elmgeo_face
 
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-19
  !> @brief   Find a segment in local numbering
  !> @details Find a segment in local numbering
  !> 
  !-----------------------------------------------------------------------

  function elmgeo_segment(ielty,iedge,isegm) result(lnode)

    integer(ip), intent(in) :: ielty
    integer(ip), intent(in) :: iedge
    integer(ip), intent(in) :: isegm
    integer(ip)             :: lnode(2)
    integer(ip)             :: pnodb,pblty

    pblty = element_type(ielty) % type_edges(iedge)
    pnodb = element_type(pblty) % number_nodes

    if( isegm == 1 ) then
       lnode = [element_type(ielty) % list_edges(1_ip,iedge),element_type(ielty) % list_edges(min(pnodb,3_ip),iedge)]
    else if( isegm == pnodb-1 ) then
       lnode = [element_type(ielty) % list_edges(pnodb,iedge),element_type(ielty) % list_edges(2_ip,iedge)]
    else
       lnode = [element_type(ielty) % list_edges(isegm+1,iedge),element_type(ielty) % list_edges(isegm+2,iedge)]
    end if
       
  end function elmgeo_segment

end module mod_elmgeo
!> @}
