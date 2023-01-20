!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @defgroup Maths_Toolbox
!> Toolbox for mathematical solvers
!> @{
!> @name    ToolBox for mathematics operations
!> @file    mod_maths.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for basic geometry
!> @details ToolBox for basic geometry
!
!-----------------------------------------------------------------------

module mod_maths_geometry

  use def_kintyp_basic, only : ip,rp,lg
  use def_parame,       only : pi
  use mod_maths_basic,  only : maths_cross_product
  use def_maths,        only : memor
  use mod_memory_basic, only : memory_alloca
  use mod_memory_basic, only : memory_deallo
  implicit none
  private

  real(rp), parameter :: zeror = epsilon(1.0_rp)

  public :: maths_point_face_distance        ! Point to face distance
  public :: maths_point_plane_distance       ! Point to plane distance
  public :: maths_point_line_distance        ! Point to line distance
  public :: maths_point_segment_distance     ! Point to segment distance
  public :: maths_point_in_triangle          ! Check if a point is inside a triangle
  public :: maths_point_in_box               ! Point in a box
  public :: maths_vectorial_product          ! Vectorial product
  public :: maths_normal_to_triangle         ! Normal to a triangle
  public :: maths_circumradius_tetrahedron   ! Circumradius of a tetra
  public :: maths_circumradius_quadrilateral ! Circumradius of a quadrilateral
  public :: maths_circumradius_triangle      ! Circumradius of a triangle
  public :: maths_area_triangle              ! Area of a triangle
  public :: maths_volume_tetrahedron         ! Volume of a treta
  public :: maths_local_orthonormal_basis    ! Construct an orthonormal basis
  public :: maths_min_max_box_vertices       ! Compute the min and max distance to a box vertices
  public :: maths_min_box_vertices           ! Compute the min distance to a box vertices
  public :: maths_max_box_vertices           ! Compute the max distance to a box vertices
  public :: maths_distance_point_segment     ! Distance and projection of a point to a segment
  public :: maths_box_measure                ! Compute a box measure (length, area, vol)
  public :: maths_boxes_union_measure        ! Compute the measure of a union of boxes

  public :: maths_normalize_basis            ! Basis nornalization
  public :: maths_angle_of_a_vector          ! Angle of a vector

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-03
  !> @brief   Distance to a box
  !> @details Compute the min and max distance to a box vertices
  !> 
  !-----------------------------------------------------------------------

  pure function maths_min_max_box_vertices(x,bobox) result(dista)

    real(rp), intent(in) :: x(:)
    real(rp), intent(in) :: bobox(:,:)
    real(rp)             :: dista(2)
    real(rp)             :: temp
    integer(ip)          :: ii

    dista = 0.0_rp
    if( size(x) == 2 ) then
       do ii = 1,2
          temp     = max(0.0_rp,bobox(1,ii)-x(ii))  
          temp     = max(temp,x(ii)-bobox(2,ii))  
          dista(1) = dista(1) + temp * temp
          temp     = max(0.0_rp,bobox(2,ii)-x(ii))  
          temp     = max(temp,x(ii)-bobox(1,ii))  
          dista(2) = dista(2) + temp * temp
       end do
    else
       do ii = 1,3
          temp     = max(0.0_rp,bobox(1,ii)-x(ii))  
          temp     = max(temp,x(ii)-bobox(2,ii))  
          dista(1) = dista(1) + temp * temp
          temp     = max(0.0_rp,bobox(2,ii)-x(ii))  
          temp     = max(temp,x(ii)-bobox(1,ii))  
          dista(2) = dista(2) + temp * temp
       end do
    end if
    dista = sqrt(dista)
    
  end function maths_min_max_box_vertices

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-03
  !> @brief   Distance to a box
  !> @details Compute the min distance to a box vertices
  !> 
  !-----------------------------------------------------------------------

  pure function maths_min_box_vertices(x,bobox) result(dista)

    real(rp), intent(in) :: x(:)
    real(rp), intent(in) :: bobox(:,:)
    real(rp)             :: dista
    real(rp)             :: temp
    integer(ip)          :: ii

    dista = 0.0_rp
    if( size(x) == 2 ) then
       do ii = 1,2
          temp  = max(0.0_rp,bobox(1,ii)-x(ii))  
          temp  = max(temp,x(ii)-bobox(2,ii))  
          dista = dista + temp * temp
       end do
    else
       do ii = 1,3
          temp  = max(0.0_rp,bobox(1,ii)-x(ii))  
          temp  = max(temp,x(ii)-bobox(2,ii))  
          dista = dista + temp * temp
       end do
    end if
    dista = sqrt(dista)
    
  end function maths_min_box_vertices

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-03
  !> @brief   Distance to a box
  !> @details Compute the min distance to a box vertices
  !> 
  !-----------------------------------------------------------------------

  pure function maths_max_box_vertices(x,bobox) result(dista)

    real(rp), intent(in) :: x(:)
    real(rp), intent(in) :: bobox(:,:)
    real(rp)             :: dista
    real(rp)             :: temp
    integer(ip)          :: ii

    dista = 0.0_rp
    if( size(x) == 2 ) then
       do ii = 1,2
          temp  = max(0.0_rp,bobox(2,ii)-x(ii))  
          temp  = max(temp,x(ii)-bobox(1,ii))  
          dista = dista + temp * temp
       end do
    else
       do ii = 1,3
          temp  = max(0.0_rp,bobox(2,ii)-x(ii))  
          temp  = max(temp,x(ii)-bobox(1,ii))  
          dista = dista + temp * temp
       end do
    end if
    dista = sqrt(dista)

  end function maths_max_box_vertices

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-03
  !> @brief   Distance point to plane
  !> @details Compute the distance of a point to a plane
  !>          Handles point to line distance as well
  !> 
  !-----------------------------------------------------------------------

  pure function maths_point_plane_distance(x,a,b,c,d) result(dista)

    real(rp), intent(in) :: a
    real(rp), intent(in) :: b
    real(rp), intent(in) :: c
    real(rp), intent(in) :: d
    real(rp), intent(in) :: x(3)
    real(rp)             :: dista

    dista = (a*x(1)+b*x(2)+c*x(3)+d)/sqrt(a*a+b*b+c*c)

  end function maths_point_plane_distance

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-03
  !> @brief   Distance point to plane
  !> @details Compute the distance of a point to a plane
  !> 
  !-----------------------------------------------------------------------

  pure function maths_point_line_distance(x,a,b,c) result(dista)

    real(rp), intent(in) :: a
    real(rp), intent(in) :: b
    real(rp), intent(in) :: c
    real(rp), intent(in) :: x(2)
    real(rp)             :: dista

    dista = (a*x(1)+b*x(2)+c)/sqrt(a*a+b*b)

  end function maths_point_line_distance

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-21
  !> @brief   ???
  !> @details Determine if a point is inside a triangle using the same side technique
  !>          point: point
  !>          facoo(dimension,vertices): triangle coordinates
  !>          ifoun: If is equal to 1, the point is inside the triangle, 0 otherwise
  !> 
  !-----------------------------------------------------------------------

  subroutine maths_point_in_triangle(pnodb,point,facoo,ifoun,bari1,bari2,ntria)

    integer(ip),   intent(in)     :: pnodb
    integer(ip),   intent(out)    :: ifoun
    real(rp),      intent(in)     :: point(3)
    real(rp),      intent(in)     :: facoo(3,*)
    real(rp),      intent(out)    :: bari1,bari2
    integer(ip),   intent(out)    :: ntria
    real(rp)                      :: v0(3),v1(3),v2(3)
    real(rp)                      :: dot00,dot01,dot02,dot11,dot12
    real(rp)                      :: invDenom
    !
    ! 3D
    !
    v0(1) = facoo(1,3) - facoo(1,1)
    v0(2) = facoo(2,3) - facoo(2,1)
    v0(3) = facoo(3,3) - facoo(3,1)

    v1(1) = facoo(1,2) - facoo(1,1)
    v1(2) = facoo(2,2) - facoo(2,1)
    v1(3) = facoo(3,2) - facoo(3,1)

    v2(1) = point(1)   - facoo(1,1)
    v2(2) = point(2)   - facoo(2,1)
    v2(3) = point(3)   - facoo(3,1)

    dot00 = v0(1)*v0(1) + v0(2)*v0(2) + v0(3)*v0(3)
    dot01 = v0(1)*v1(1) + v0(2)*v1(2) + v0(3)*v1(3)
    dot02 = v0(1)*v2(1) + v0(2)*v2(2) + v0(3)*v2(3)
    dot11 = v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3)
    dot12 = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)  
    !
    ! Compute barycentric coordinates
    !
    ntria = 0
    if (abs(dot00 * dot11 - dot01 * dot01) > zeror) then
       invDenom = 1.0_rp / (dot00 * dot11 - dot01 * dot01)
       bari1 = (dot11 * dot02 - dot01 * dot12) * invDenom
       bari2 = (dot00 * dot12 - dot01 * dot02) * invDenom
       !
       ! Check
       !
       if( bari1 >= 0.0_rp .and. bari2 >= 0.0_rp .and. bari1 + bari2 <= 1.0_rp ) then
          ntria = 1
          ifoun = 1
       else
          ifoun = 0
       end if
    else
       ifoun = 0
    end if
    if( pnodb == 4 .and. ifoun == 0 ) then
       v0(1) = facoo(1,4) - facoo(1,1)
       v0(2) = facoo(2,4) - facoo(2,1)
       v0(3) = facoo(3,4) - facoo(3,1)

       v1(1) = facoo(1,3) - facoo(1,1)
       v1(2) = facoo(2,3) - facoo(2,1)
       v1(3) = facoo(3,3) - facoo(3,1)

       v2(1) = point(1)   - facoo(1,1)
       v2(2) = point(2)   - facoo(2,1)
       v2(3) = point(3)   - facoo(3,1)


       dot00 = v0(1)*v0(1) + v0(2)*v0(2) + v0(3)*v0(3)
       dot01 = v0(1)*v1(1) + v0(2)*v1(2) + v0(3)*v1(3)
       dot02 = v0(1)*v2(1) + v0(2)*v2(2) + v0(3)*v2(3)
       dot11 = v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3)
       dot12 = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)  

       if (abs(dot00 * dot11 - dot01 * dot01) > zeror) then
          invDenom = 1.0_rp / (dot00 * dot11 - dot01 * dot01)
          bari1 = (dot11 * dot02 - dot01 * dot12) * invDenom
          bari2 = (dot00 * dot12 - dot01 * dot02) * invDenom
          !
          ! Check
          !
          if( bari1 >= 0.0_rp .and. bari2 >= 0.0_rp .and. bari1 + bari2 <= 1.0_rp ) then
             ifoun = 1
             ntria = 2
          else
             ifoun = 0
          end if
       else
          ifoun = 0
       end if
    end if

  end subroutine maths_point_in_triangle

  !-----------------------------------------------------------------------
  !>
  !> @date    08/01/2018
  !> @author  Guillaume Houzeaux
  !> @brief   Different elements in an array
  !> @details Two and three-dimensional vectorial product of two vectors  v3 = v1 x v2.
  !>          The same pointer as for v1 or v2 may be used for v3. If N = 2, it is
  !>          assumed that v1 = (0,0,v1_3) and v2 = (v2_1,v2_2,0).
  !>
  !-----------------------------------------------------------------------

  pure subroutine maths_vectorial_product(v1,v2,v3,n)

    integer(ip), intent(in)  :: n
    real(rp),    intent(in)  :: v1(n)
    real(rp),    intent(in)  :: v2(n)
    real(rp),    intent(out) :: v3(3)

    if( n == 2 ) then
       v3(1:2) = 0.0_rp
       v3(3)   = v1(1)*v2(2)-v1(2)*v2(1)
    else if( n == 3 ) then
       v3(1)   = v1(2)*v2(3)-v1(3)*v2(2)
       v3(2)   = v1(3)*v2(1)-v1(1)*v2(3)
       v3(3)   = v1(1)*v2(2)-v1(2)*v2(1)
    end if

  end subroutine maths_vectorial_product

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-21
  !> @brief   Normal to triangle
  !> @details Compute the normal to a triangle
  !> 
  !-----------------------------------------------------------------------

  function maths_normal_to_triangle(p1,p2,p3,coord) result( normal )

    integer(ip), intent(in)  :: p1,p2,p3
    real(rp),    intent(in)  :: coord(3,*)
    real(rp)                 :: normal(3)
    real(rp)                 :: v1(3),v2(3)

    v1(1:3) = coord(1:3,p2) - coord(1:3,p1)
    v2(1:3) = coord(1:3,p3) - coord(1:3,p1)

    call maths_vectorial_product(v1,v2,normal,3_ip)

  end function maths_normal_to_triangle

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-05
  !> @brief   Circumradius
  !> @details Circumradius of a triangle
  !> 
  !-----------------------------------------------------------------------

  real(rp) pure function maths_circumradius_triangle(x1,x2,x3) result(R)

    real(rp), intent(in) :: x1(2)
    real(rp), intent(in) :: x2(2)
    real(rp), intent(in) :: x3(2)
    real(rp)             :: v(2,3),a,b,c,S,cross(3)

    v(1:2,1)  = x1(1:2)-x2(1:2) ! a: 1-2
    v(1:2,2)  = x1(1:2)-x3(1:2) ! b: 1-3
    v(1:2,3)  = x2(1:2)-x3(1:2) ! c: 2-3
    a         = sqrt(dot_product(v(1:2,1),v(1:2,1)))
    b         = sqrt(dot_product(v(1:2,2),v(1:2,2)))
    c         = sqrt(dot_product(v(1:2,3),v(1:2,3)))
    cross     = maths_cross_product(v(:,1),v(:,2),2_ip) 
    S         = 0.5_rp * sqrt(dot_product(cross,cross)) 
    R         = a*b*c/(4.0_rp*S)

  end function maths_circumradius_triangle

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-05
  !> @brief   Circumradius
  !> @details Circumradius of a quadrilateral
  !> 
  !-----------------------------------------------------------------------

  real(rp) pure function maths_circumradius_quadrilateral(x1,x2,x3,x4) result(R)

    real(rp), intent(in) :: x1(2)
    real(rp), intent(in) :: x2(2)
    real(rp), intent(in) :: x3(2)
    real(rp), intent(in) :: x4(2)
    real(rp)             :: v(2,4),a,b,c,d,s

    v(1:2,1)  = x1(1:2)-x2(1:2) ! a: 1-2
    v(1:2,2)  = x2(1:2)-x3(1:2) ! b: 2-3
    v(1:2,3)  = x3(1:2)-x4(1:2) ! c: 3-4
    v(1:2,4)  = x4(1:2)-x1(1:2) ! c: 1-4
    a         = sqrt(dot_product(v(1:2,1),v(1:2,1)))
    b         = sqrt(dot_product(v(1:2,2),v(1:2,2)))
    c         = sqrt(dot_product(v(1:2,3),v(1:2,3)))
    d         = sqrt(dot_product(v(1:2,4),v(1:2,4)))
    s         = 0.5_rp*(a+b+c+d)
    R         = 0.25_rp*sqrt(( (a*c+b*d)*(a*d+b*c)*(a*b+c*d) ) / ( (s-a)*(s-b)*(s-c)*(s-d) ) )

  end function maths_circumradius_quadrilateral

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-05
  !> @brief   Circumradius
  !> @details Circumradius of a tet
  !> 
  !-----------------------------------------------------------------------

  real(rp) pure function maths_circumradius_tetrahedron(x1,x2,x3,x4) result(R)

    real(rp), intent(in) :: x1(3)
    real(rp), intent(in) :: x2(3)
    real(rp), intent(in) :: x3(3)
    real(rp), intent(in) :: x4(3)
    real(rp)             :: v(3,6),aa,bb,cc,Vol,edgel(6)
    real(rp)             :: f1,f2,f3,f4

    Vol        = maths_volume_tetrahedron(x1,x2,x3,x4)

    v(1:3,1)   = x1(1:3)-x2(1:3) ! a: 1-2
    v(1:3,2)   = x1(1:3)-x3(1:3) ! b: 1-3
    v(1:3,3)   = x1(1:3)-x4(1:3) ! c: 1-4

    v(1:3,4)   = x3(1:3)-x4(1:3) ! A: 3-4
    v(1:3,5)   = x2(1:3)-x4(1:3) ! B: 2-4
    v(1:3,6)   = x2(1:3)-x3(1:3) ! C: 2-3

    edgel(1)   = sqrt(dot_product(v(1:3,1),v(1:3,1)))
    edgel(2)   = sqrt(dot_product(v(1:3,2),v(1:3,2)))
    edgel(3)   = sqrt(dot_product(v(1:3,3),v(1:3,3)))
    edgel(4)   = sqrt(dot_product(v(1:3,4),v(1:3,4)))
    edgel(5)   = sqrt(dot_product(v(1:3,5),v(1:3,5)))
    edgel(6)   = sqrt(dot_product(v(1:3,6),v(1:3,6)))

    aa         = edgel(1) * edgel(4)
    bb         = edgel(2) * edgel(5)
    cc         = edgel(3) * edgel(6)          
    f1         =   aa + bb + cc
    f2         =   aa + bb - cc
    f3         =   aa - bb + cc 
    f4         = - aa + bb + cc

    R          = sqrt( f1*f2*f3*f4 ) / (24.0_rp*Vol)

  end function maths_circumradius_tetrahedron

  !-----------------------------------------------------------------------
  !>
  !> @date    08/01/2018
  !> @author  Guillaume Houzeaux
  !> @brief   Triangle area
  !> @details Compute the area of a triangle
  !>
  !-----------------------------------------------------------------------

  pure function maths_area_triangle(x1,x2,x3,n) result(S)

    integer(ip), intent(in) :: n
    real(rp),    intent(in) :: x1(n)
    real(rp),    intent(in) :: x2(n)
    real(rp),    intent(in) :: x3(n)
    real(rp)                :: S
    real(rp)                :: v1(3)
    real(rp)                :: v2(3)
    real(rp)                :: v3(3)

    if( n == 2 ) then
       v1(1:2) = x2(1:2)-x1(1:2)
       v2(1:2) = x3(1:2)-x1(1:2)
       v3(1:2) = 0.0_rp
       v3(3)   = v1(1)*v2(2)-v1(2)*v2(1)
    else if( n == 3 ) then
       v1(1:3) = x2(1:3)-x1(1:3)
       v2(1:3) = x3(1:3)-x1(1:3)
       v3(1)   = v1(2)*v2(3)-v1(3)*v2(2)
       v3(2)   = v1(3)*v2(1)-v1(1)*v2(3)
       v3(3)   = v1(1)*v2(2)-v1(2)*v2(1)
    end if

    S = 0.5_rp * sqrt(dot_product(v3(1:3),v3(1:3)))

  end function maths_area_triangle

  !-----------------------------------------------------------------------
  !>
  !> @date    08/01/2018
  !> @author  Guillaume Houzeaux
  !> @brief   Triangle area
  !> @details Compute the volume of a tetrahedron
  !>          https://en.wikipedia.org/wiki/Tetrahedron#Volume
  !>
  !-----------------------------------------------------------------------

  pure function maths_volume_tetrahedron(x1,x2,x3,x4) result(V)

    real(rp),    intent(in) :: x1(3)
    real(rp),    intent(in) :: x2(3)
    real(rp),    intent(in) :: x3(3)
    real(rp),    intent(in) :: x4(3)
    real(rp)                :: V
    real(rp)                :: v1(3)
    real(rp)                :: v2(3)
    real(rp)                :: v3(3)
    real(rp)                :: v4(3)

    v1(1:3) = x1(1:3)-x4(1:3)
    v2(1:3) = x2(1:3)-x4(1:3)
    v3(1:3) = x3(1:3)-x4(1:3)
    v4      = maths_cross_product(v2,v3,3_ip)
    V       = abs(dot_product(v1(1:3),v4(1:3)))/6.0_rp

  end function maths_volume_tetrahedron

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-21
  !> @brief   ???
  !> @details Minimun distance between a point and a face
  !>          point:  point
  !>          bocod:  triangle coordinates
  !>          norma:  nornmalized exterior normal to the triangle
  !>          dista:  minimum distance between the point and the triangle
  !>          factor: inside or outside
  !>          proje:  nearest triangle point projection  
  !> 
  !-----------------------------------------------------------------------

  subroutine maths_point_face_distance(&
       pdime,pnodb,xcoor,bocod,fabox,norma,dista,factor,proje,TOLERANCE)

    integer(ip),             intent(in)  :: pdime
    integer(ip),             intent(in)  :: pnodb
    real(rp),                intent(in)  :: xcoor(pdime)
    real(rp),                intent(in)  :: norma(pdime)
    real(rp),                intent(in)  :: bocod(pdime,pnodb)
    real(rp),                intent(in)  :: fabox(2,pdime)
    real(rp),                intent(out) :: dista
    real(rp),                intent(out) :: factor
    real(rp),                intent(out) :: proje(pdime)
    real(rp),      optional, intent(in)  :: TOLERANCE
    integer(ip)                          :: idime,ifoun,ntria
    real(rp)                             :: pladi,plapo(pdime)
    real(rp)                             :: dist1,dist2,dist3,dist4
    real(rp)                             :: proj1(3),proj2(3),proj3(3),proj4(3)
    real(rp)                             :: bari1,bari2
    !
    ! Distance between the point and the plane formed by the face
    ! Dot product between the exterior normal and the vector build from the first triangle vertex
    ! to the point  
    !
    pladi = 0.0_rp
    do idime = 1,pdime
       pladi = pladi + norma(idime) * ( xcoor(idime) - bocod(idime,1) )
    end do
    !
    ! Point projection on the plane
    !  
    plapo(1)     = xcoor(1)     - pladi*norma(1)
    plapo(2)     = xcoor(2)     - pladi*norma(2)
    plapo(pdime) = xcoor(pdime) - pladi*norma(pdime)

    if ( pladi < 0.0_rp ) then
       factor = -1.0_rp
    else
       factor =  1.0_rp
    end if
    ifoun = 0

    if( pdime == 3 ) then
       !
       ! Check if we are in bounding box
       !
       if( plapo(1) >= fabox(1,1) .and. plapo(1) <= fabox(2,1) ) then
          if( plapo(2) >= fabox(1,2) .and. plapo(2) <= fabox(2,2) ) then
             if( plapo(3) >= fabox(1,3) .and. plapo(3) <= fabox(2,3) ) then
                !
                ! Determine if the projection point on plane is inside the triangle
                !
               call maths_point_in_triangle(pnodb,plapo,bocod,ifoun,bari1,bari2,ntria)
             end if
          end if
       end if
       !
       ! The projection point is inside the triangle
       !
       if( ifoun == 1 ) then
          
          dista    = pladi        
          proje(1) = plapo(1)
          proje(2) = plapo(2)
          proje(3) = plapo(3)
          
       else
          
          if( pnodb == 4 ) then
                    
             call maths_point_segment_distance(pdime,xcoor,bocod(1,1),bocod(1,2),dist1,proj1,ifoun,TOLERANCE)
             call maths_point_segment_distance(pdime,xcoor,bocod(1,2),bocod(1,4),dist2,proj2,ifoun,TOLERANCE)
             call maths_point_segment_distance(pdime,xcoor,bocod(1,3),bocod(1,4),dist3,proj3,ifoun,TOLERANCE)         
             call maths_point_segment_distance(pdime,xcoor,bocod(1,1),bocod(1,4),dist4,proj4,ifoun,TOLERANCE)         

             if( dist1 < dist2 .and. dist1 < dist3 .and. dist1 < dist4 ) then              
                dista    = dist1 
                proje(1) = proj1(1)
                proje(2) = proj1(2)
                proje(3) = proj1(3)
             else if( dist2 < dist1 .and. dist2 < dist3 .and. dist2 < dist4 ) then        
                dista    = dist2
                proje(1) = proj2(1)
                proje(2) = proj2(2)
                proje(3) = proj2(3)
             else if( dist3 < dist1 .and. dist3 < dist2 .and. dist3 < dist4 ) then        
                dista    = dist3 
                proje(1) = proj3(1)
                proje(2) = proj3(2)
                proje(3) = proj3(3)              
             else
                dista    = dist4
                proje(1) = proj4(1)
                proje(2) = proj4(2)
                proje(3) = proj4(3)
             end if

          else
             
             call maths_point_segment_distance(pdime,xcoor,bocod(1,1),bocod(1,2),dist1,proj1,ifoun,TOLERANCE)
             call maths_point_segment_distance(pdime,xcoor,bocod(1,2),bocod(1,3),dist2,proj2,ifoun,TOLERANCE)
             call maths_point_segment_distance(pdime,xcoor,bocod(1,1),bocod(1,3),dist3,proj3,ifoun,TOLERANCE)
             
             if( dist1 < dist2 .and. dist1 < dist3 ) then              
                dista    = dist1
                proje(1) = proj1(1)
                proje(2) = proj1(2)
                proje(3) = proj1(3)
             else if( dist2 < dist1 .and. dist2 < dist3 ) then
                dista    = dist2 
                proje(1) = proj2(1)
                proje(2) = proj2(2)
                proje(3) = proj2(3)
             else
                dista    = dist3 
                proje(1) = proj3(1)
                proje(2) = proj3(2)
                proje(3) = proj3(3)
             end if
          end if
       end if

    else if( pdime == 2 ) then

       call maths_point_segment_distance(pdime,xcoor,bocod(1,1),bocod(1,2),dista,proje,ifoun,TOLERANCE)           

    end if

  end subroutine maths_point_face_distance

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-21
  !> @brief   Point to segment distance
  !> @details Minimun distance between a point and a segment
  !>          xcoor : point xcoorinates 
  !>          coor1,coor2 : defines the segment
  !>          ndime: dimension
  !>          dista: distance
  !>          proje: projection of the point on the segment
  !>          ifoun = 1 the projection point is inside the segment
  !>          ifoun = 0 the projection point is outside the segment
  !> 
  !-----------------------------------------------------------------------

  pure subroutine maths_point_segment_distance(pdime,xcoor,coor1,coor2,dista,proje,ifoun,TOLERANCE)

    integer(ip),             intent(in)  :: pdime
    real(rp),                intent(in)  :: xcoor(pdime)
    real(rp),                intent(in)  :: coor1(pdime)
    real(rp),                intent(in)  :: coor2(pdime)
    real(rp),                intent(out) :: dista
    real(rp),                intent(out) :: proje(pdime)
    integer(ip),             intent(out) :: ifoun
    real(rp),      optional, intent(in)  :: TOLERANCE
    real(rp)                             :: numer,denom,dsegm,toler

    if( present(TOLERANCE) ) then
       toler = TOLERANCE
    else
       toler = 0.01_rp
    end if
    
    numer = dot_product(coor2(1:pdime) - coor1(1:pdime),xcoor(1:pdime) - coor1(1:pdime))
    denom = dot_product(coor2(1:pdime) - coor1(1:pdime),coor2(1:pdime) - coor1(1:pdime)) 
    dsegm = numer / (denom+zeror)

    if( dsegm < -toler ) then
       
       ifoun        = 0_ip
       proje(1)     = coor1(1)
       proje(2)     = coor1(2)
       proje(pdime) = coor1(pdime)

    else if ( dsegm >= -toler .and. dsegm <= 1.0_rp+toler ) then

       ifoun        = 1_ip
       proje(1)     = coor1(1)     + dsegm * (coor2(1)     - coor1(1))
       proje(2)     = coor1(2)     + dsegm * (coor2(2)     - coor1(2))
       proje(pdime) = coor1(pdime) + dsegm * (coor2(pdime) - coor1(pdime))

    else

       ifoun        = 0_ip
       proje(1)     = coor2(1)
       proje(2)     = coor2(2)
       proje(pdime) = coor2(pdime)

    end if

    dista = sqrt(dot_product(xcoor(1:pdime)-proje(1:pdime),xcoor(1:pdime)-proje(1:pdime)))

  end subroutine maths_point_segment_distance

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Location in a bin
  !> @details Find the location of a point in a bin
  !> 
  !-----------------------------------------------------------------------

  pure logical(lg) function maths_point_in_box(comin,comax,coord) 

    real(rp), intent(in)  :: comin(:)
    real(rp), intent(in)  :: comax(:)
    real(rp), intent(in)  :: coord(:) !< Coordinate of the test point

    maths_point_in_box = .true.

    select case ( size(coord) )

    case ( 1_ip )
       
       if( coord(1) > comax(1) .or. coord(1) < comin(1) ) maths_point_in_box = .false.
       
    case ( 2_ip )
       
       if(     coord(1) > comax(1) ) then
          maths_point_in_box = .false. ; return
       else if( coord(1) < comin(1) ) then
          maths_point_in_box = .false. ; return
       else if( coord(2) > comax(2) ) then
          maths_point_in_box = .false. ; return
       else if( coord(2) < comin(2) ) then
          maths_point_in_box = .false. ; return
       end if
       
     case ( 3_ip )
      
       if(     coord(1) > comax(1) ) then
          maths_point_in_box = .false. ; return
       else if( coord(1) < comin(1) ) then
          maths_point_in_box = .false. ; return
       else if( coord(2) > comax(2) ) then
          maths_point_in_box = .false. ; return
       else if( coord(2) < comin(2) ) then
          maths_point_in_box = .false. ; return
       else if( coord(3) > comax(3) ) then
          maths_point_in_box = .false. ; return
       else if( coord(3) < comin(3) ) then
          maths_point_in_box = .false. ; return
       end if
       
    end select

  end function maths_point_in_box
    
  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    04/12/2015
  !> @brief   Computes tangent vectors
  !> @details Computes tangent vectors given a normal vector which form
  !>          an orthogonal basis
  !>
  !----------------------------------------------------------------------

  pure subroutine maths_normalize_basis(basis)
    
    real(rp), intent(inout) :: basis(:,:)  !< Local basis: normal is BASIS(1:NDIME,1)
    integer(ip)             :: idime

    do idime = 1,int(size(basis,2),ip)
       basis(:,idime) = basis(:,idime) / sqrt(dot_product(basis(:,idime),basis(:,idime)))       
    end do

  end subroutine maths_normalize_basis
  
  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    04/12/2015
  !> @brief   Computes tangent vectors
  !> @details Computes tangent vectors given a normal vector which form
  !>          an orthogonal basis
  !>
  !----------------------------------------------------------------------

  subroutine maths_local_orthonormal_basis(ndime,basis,ierro)

    integer(ip), intent(in)              :: ndime               !< Dimension
    real(rp),    intent(inout)           :: basis(ndime,ndime)  !< Local basis: normal is BASIS(1:NDIME,1)
    integer(ip), intent(inout), optional :: ierro               !< Error counter
    integer(ip)                          :: idime,kdime
    real(rp)                             :: xnor1,xnor2
    real(rp)                             :: xnor3,xnoma
    real(rp)                             :: exwor(3)

    if( ndime == 2 ) then

       basis(1,2) = -basis(2,1)
       basis(2,2) =  basis(1,1)

    else if( ndime == 3 ) then
       !
       ! Look for e_k such that n x e_k is maximum
       !
       xnor1    = basis(1,1) * basis(1,1)
       xnor2    = basis(2,1) * basis(2,1)
       xnor3    = basis(3,1) * basis(3,1)
       exwor(1) = xnor2 + xnor3
       exwor(2) = xnor1 + xnor3
       exwor(3) = xnor1 + xnor2
       xnoma    = 0.0_rp
       kdime    = 0
       do idime = 1,3
          if( exwor(idime) > xnoma ) then
             xnoma = exwor(idime)
             kdime = idime
          end if
       end do
       xnoma = 1.0_rp / sqrt(xnoma)
       !
       ! Set t_1 = e_k x n, first tangent vector
       !
       if( kdime == 1 ) then
          basis(1,2) =  0.0_rp
          basis(2,2) = -basis(3,1) * xnoma
          basis(3,2) =  basis(2,1) * xnoma
       else if( kdime == 2 ) then
          basis(1,2) =  basis(3,1) * xnoma
          basis(2,2) =  0.0_rp
          basis(3,2) = -basis(1,1) * xnoma
       else if( kdime == 3 ) then
          basis(1,2) = -basis(2,1) * xnoma
          basis(2,2) =  basis(1,1) * xnoma
          basis(3,2) =  0.0_rp
       else
          if( present(ierro) ) then
             ierro = ierro + 1
          else
             call runend('MATHS_LOCAL_BASIS: WRONG NORMAL')
          end if
       end if
       !
       ! Set t_2 = n x t_1, second tangent vector
       !
       basis(1,3) = basis(2,1) * basis(3,2) - basis(3,1) * basis(2,2)
       basis(2,3) = basis(3,1) * basis(1,2) - basis(1,1) * basis(3,2)
       basis(3,3) = basis(1,1) * basis(2,2) - basis(2,1) * basis(1,2)

    end if

  end subroutine maths_local_orthonormal_basis

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    04/12/2015
  !> @brief   Shortest distance point segment
  !> @details Shortest distance of a point to a segment
  !>
  !----------------------------------------------------------------------

  pure subroutine maths_distance_point_segment(x1,x2,xx,pp,d) 

    real(rp), intent(in)  :: x1(:)
    real(rp), intent(in)  :: x2(:)
    real(rp), intent(in)  :: xx(:)
    real(rp), intent(out) :: pp(:)
    real(rp), intent(out) :: d
    integer(ip)           :: nn
    real(rp)              :: t,xn

    nn       = int(size(xx),ip)
    pp(1:nn) = x2(1:nn)-x1(1:nn)
    xn       = dot_product(pp(1:nn),pp(1:nn))
    if( nn == 2 ) then
       t     =  ( (xx(1) - x1(1)) * pp(1) + (xx(2) - x1(2)) * pp(2) ) / xn
    else
       t     =  ( (xx(1) - x1(1)) * pp(1) + (xx(2) - x1(2)) * pp(2) +  (xx(3) - x1(3)) * pp(3) ) / xn       
    end if
    t        = min(1.0_rp,t)
    t        = max(0.0_rp,t)

    pp(1:nn) = x1(1:nn) + t * pp(1:nn)
    d        = sqrt(dot_product(pp(1:nn) - xx(1:nn),pp(1:nn) - xx(1:nn)))
    
  end subroutine maths_distance_point_segment
  
  !----------------------------------------------------------------------
  !>
  !> @author  Ricard Borrell
  !> @date    04/02/2021
  !> @brief   Computes box measure
  !> @details Computes the length, area or volume of a box
  !>
  !----------------------------------------------------------------------
  
  real(rp) pure function maths_box_measure(ndime,box)

     implicit none 
     integer(ip), intent(in)  :: ndime
     real(rp),    intent(in)  :: box(2*ndime)
     integer(ip)              :: ii 

     maths_box_measure = 1.0_rp
     do ii = 1, ndime
        maths_box_measure = maths_box_measure * (box(ndime+ii)-box(ii)) 
     end do

   end function maths_box_measure
   
  !----------------------------------------------------------------------
  !>
  !> @author  Ricard Borrell
  !> @date    04/02/2021
  !> @brief   Computes measure of boxes union
  !> @details Computes measure (length, area or volume) of a union of boxes.
  !>          Of course overlapping reagions are counted only once
  !>          ndime: space dimention
  !>          nbox:  number of boxes
  !>          lbox:  lists of boxes 
  !>
  !----------------------------------------------------------------------
  real(rp) function maths_boxes_union_measure(ndime,nbox,lbox)
     
     implicit none 
     integer(ip)       :: ndime
     integer(ip)       :: nbox
     real(rp), pointer     :: lbox(:,:)

     real(rp), pointer     :: buffin(:,:)
     real(rp), pointer     :: buffout(:,:)
     real(rp), pointer     :: lsubbox(:,:)
     integer(ip)           :: nbin,nbout
     integer(ip)           :: iin,iout,kk
     integer(ip)           :: nsubb,ndimb
     integer(ip)           :: bdim, sdim
     real(rp)              :: res

     bdim = nbox*(3**ndime)
     sdim = 2*(3**(ndime-1))
     nullify(buffin,buffout,lsubbox)
     call memory_alloca(memor,'buffin' ,"maths_boxes_union_measure",buffin, 2*ndime,bdim)
     call memory_alloca(memor,'buffout',"maths_boxes_union_measure",buffout,2*ndime,bdim)
     call memory_alloca(memor,'lsubbox',"maths_boxes_union_measure",lsubbox,2*ndime,sdim)

     nbin   = nbox
     nbout  = 0
     iin    = 1
     ndimb  = 2*ndime
     buffin(1:ndimb,1:nbox)=lbox(1:ndimb,1:nbox) 
     res    = 0.0_rp

     do while (iin <= nbin)
        nsubb = 0
        iout    = 1
        do while (nsubb == 0 .and. iout<=nbout)
           call maths_conform_boxes(ndime,buffin(1:ndimb,iin),buffout(1:ndimb,iout),nsubb,lsubbox)
           !
           ! Add to buffin the new boxes
           !
           do kk=1,nsubb
              nbin = nbin + 1
              if(nbout > bdim ) then
                 call runend('MATHS_BOXES_UNION_MEASURE: buffer too small, make bdim larger')
              endif 
              buffin(1:ndimb,nbin) = lsubbox(1:ndimb,kk)
           enddo
           if(nsubb == -1) then
              !
              ! Box in buffin contains box in buffout
              ! In this case, box in buffout is eliminated
              !
              res = res - maths_box_measure(ndime,buffout(1:ndimb,iout))
              do kk = iout,nbout-1
                 buffout(1:ndimb,kk)=buffout(1:ndimb,kk+1)
              enddo
              nbout  = nbout -1
              nsubb = 0
           else
              iout = iout + 1
           endif
        enddo
        !
        ! Box without intersections can be added to buffout
        !
        if(nsubb == 0) then
           nbout = nbout + 1
           if(nbout > bdim ) then
                 call runend('MATHS_BOXES_UNION_MEASURE: buffer too small, make bdim larger')
           endif 
           buffout(1:ndimb,nbout) = buffin(1:ndimb,iin)  
           res = res + maths_box_measure(ndime,buffout(1:ndimb,nbout))
        endif
        iin = iin+1
     enddo
     call memory_deallo(memor,'buffin' ,"maths_boxes_union_measure",buffin)
     call memory_deallo(memor,'buffout',"maths_boxes_union_measure",buffout)
     call memory_deallo(memor,'lsubbox',"maths_boxes_union_measure",lsubbox)

     maths_boxes_union_measure = res
     return 

  end function
  !----------------------------------------------------------------------
  !>
  !> @author  Ricard Borrell
  !> @date    04/02/2021
  !> @brief   Utility for boxes_union_measure function
  !> @details Given two boxes evaluates the intersection and 
  !>          decomposes box1 into conforming sub-boxes sush as 
  !>          one of the sub-boxes is the intersection
  !>          nsubb  > 0:  number of subboxes generated
  !>          nsubb = -1:  box2 contained within box1
  !>          nsubb = -2:  box1 contained within box2 
  !>
  !----------------------------------------------------------------------
  
  pure subroutine maths_conform_boxes(ndime,box1,box2,nsubb,lsubb)

     integer(ip),       intent(in)    :: ndime
     real(rp),          intent(in)    :: box1(2*ndime)
     real(rp),          intent(in)    :: box2(2*ndime)
     integer(ip),       intent(out)   :: nsubb
     real(rp), pointer, intent(inout) :: lsubb(:,:)

     real(rp)                         :: boxp(2*ndime)
     real(rp)                         :: M1,M2,Mp,tol
     real(rp)                         :: lim(4,3)
     integer(ip)                      :: nlim(3)
     integer(ip)                      :: ii,jj,kk

     tol = 1.e-06_rp
     !
     ! Evaluates intesecction stored in boxp
     !
     do ii=1,ndime
        do jj=1,2
           boxp(ii+(jj-1)*ndime) = max(box2(ii+(jj-1)*ndime),box1(ii))
           boxp(ii+(jj-1)*ndime) = min(boxp(ii+(jj-1)*ndime),box1(ii+ndime))
        enddo
     enddo
     !
     ! Evaluate meaures
     !
     M1 = maths_box_measure(ndime,box1)
     M2 = maths_box_measure(ndime,box2)
     Mp = maths_box_measure(ndime,boxp)
     !
     ! Analyse cases
     !
     if(Mp < tol) then
        !
        ! No intersection
        !
        nsubb = 0

     else if(abs(M1-Mp)/M1 < tol) then
        !
        ! box2 contains box1
        !
        nsubb = -2

     else if(abs(M2-Mp)/M2 < tol) then 
        !
        ! box2 included within box1
        !
        nsubb = -1

     else 
        !
        ! Generate new subboxes limits
        !
        nlim(:)=0
        lim(:,:)=0
        do ii =1,ndime
           nlim(ii) = nlim(ii) + 1
           lim(nlim(ii),ii)=box1(ii)
           do kk = 0, 1
              if(  boxp(ii+kk*ndime)-lim(nlim(ii),ii) > tol .and. &
                 & box1(ii+ndime)- boxp(ii+kk*ndime)  > tol) then
                 nlim(ii) = nlim(ii) + 1
                 lim(nlim(ii),ii)=boxp(ii+kk*ndime)
              endif
           enddo
           nlim(ii) = nlim(ii) + 1
           lim(nlim(ii),ii)=box1(ii+ndime)
        enddo
        !
        ! Generate new subboxes
        !
        nsubb = 0
        do ii = 1,max(nlim(1)-1_ip,1_ip)
           do jj = 1,max(nlim(2)-1_ip,1_ip)
              do kk = 1,max(nlim(3)-1_ip,1_ip)
                 nsubb = nsubb+1
                 lsubb(1,nsubb)       = lim(ii,1)
                 lsubb(1+ndime,nsubb) = lim(ii+1,1)
                 if(ndime>1) then
                    lsubb(2,nsubb)       = lim(jj,2)
                    lsubb(2+ndime,nsubb) = lim(jj+1,2)
                    if(ndime>2) then
                       lsubb(3,nsubb)       = lim(kk,3)
                       lsubb(3+ndime,nsubb) = lim(kk+1,3)
                    endif
                 endif
              enddo
           enddo
        enddo
     endif
  end subroutine


   !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    03/03/2016
  !> @brief   Angle
  !> @details Compute the angle formed by a vecfir
  !>
  !-----------------------------------------------------------------------

  subroutine maths_angle_of_a_vector(ndime,vecto,angl1,angl2)

    integer(ip), intent(in)  :: ndime        !< Dimension
    real(rp),    intent(in)  :: vecto(ndime) !< Vector
    real(rp),    intent(out) :: angl1        !< Angle 1
    real(rp),    intent(out) :: angl2        !< Angle 2
    real(rp)                 :: vect2(3)
    real(rp)                 :: cosa,sina

    if( ndime == 2 ) then
       vect2(1:ndime) = vecto(1:ndime) / sqrt(dot_product(vecto,vecto)+zeror)
       cosa  = vect2(1)
       sina  = vect2(2)
       angl1 = acos(cosa)*180.0_rp/pi
       if( sina < 0.0_rp ) angl1 = -angl1
    else
       angl1 = 0.0_rp
       angl2 = 0.0_rp
    end if

  end subroutine maths_angle_of_a_vector

  
end module mod_maths_geometry
!> @}
