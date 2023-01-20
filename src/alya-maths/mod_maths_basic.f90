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
!> @brief   ToolBox for basic maths
!> @details ToolBox for basic maths
!
!-----------------------------------------------------------------------

module mod_maths_basic

  use def_maths
  implicit none
  private
  
  integer(ip), parameter,dimension(16)  :: postb2 = &
       (/ 0,2,3,1,0,1,3,2,3,2,0,1,3,1,0,2/)
  integer(ip), parameter,dimension(16)  :: postc2 = &
       (/ 0,3,1,2,0,1,3,2,2,3,1,0,2,1,3,0/)
  integer(2),  parameter,dimension(16)  :: typtb2 = &
       (/ 1_2,0_2,0_2,2_2,0_2,1_2,1_2,3_2,3_2,2_2,2_2,0_2,2_2,3_2,3_2,1_2/)
  integer(ip), parameter,dimension(192) :: postb3 = (/ &
       0_ip,2_ip,3_ip,1_ip,5_ip,7_ip,6_ip,4_ip,0_ip,4_ip,5_ip,1_ip,3_ip,7_ip,6_ip,2_ip,&
       0_ip,4_ip,6_ip,2_ip,3_ip,7_ip,5_ip,1_ip,5_ip,4_ip,0_ip,1_ip,3_ip,2_ip,6_ip,7_ip,&
       3_ip,7_ip,5_ip,1_ip,0_ip,4_ip,6_ip,2_ip,6_ip,4_ip,5_ip,7_ip,3_ip,1_ip,0_ip,2_ip,&
       0_ip,2_ip,6_ip,4_ip,5_ip,7_ip,3_ip,1_ip,0_ip,1_ip,5_ip,4_ip,6_ip,7_ip,3_ip,2_ip,&
       6_ip,2_ip,0_ip,4_ip,5_ip,1_ip,3_ip,7_ip,5_ip,7_ip,6_ip,4_ip,0_ip,2_ip,3_ip,1_ip,&
       3_ip,2_ip,6_ip,7_ip,5_ip,4_ip,0_ip,1_ip,3_ip,2_ip,0_ip,1_ip,5_ip,4_ip,6_ip,7_ip,&
       5_ip,7_ip,3_ip,1_ip,0_ip,2_ip,6_ip,4_ip,6_ip,7_ip,3_ip,2_ip,0_ip,1_ip,5_ip,4_ip,&
       6_ip,2_ip,3_ip,7_ip,5_ip,1_ip,0_ip,4_ip,5_ip,1_ip,3_ip,7_ip,6_ip,2_ip,0_ip,4_ip,&
       0_ip,1_ip,3_ip,2_ip,6_ip,7_ip,5_ip,4_ip,3_ip,1_ip,0_ip,2_ip,6_ip,4_ip,5_ip,7_ip,&
       6_ip,4_ip,0_ip,2_ip,3_ip,1_ip,5_ip,7_ip,3_ip,7_ip,6_ip,2_ip,0_ip,4_ip,5_ip,1_ip,&
       5_ip,4_ip,6_ip,7_ip,3_ip,2_ip,0_ip,1_ip,6_ip,7_ip,5_ip,4_ip,0_ip,1_ip,3_ip,2_ip,&
       3_ip,1_ip,5_ip,7_ip,6_ip,4_ip,0_ip,2_ip,5_ip,1_ip,0_ip,4_ip,6_ip,2_ip,3_ip,7_ip/)
  integer(ip), parameter,dimension(192) :: postc3 = (/ &
       0_ip,3_ip,1_ip,2_ip,7_ip,4_ip,6_ip,5_ip,0_ip,3_ip,7_ip,4_ip,1_ip,2_ip,6_ip,5_ip,&
       0_ip,7_ip,3_ip,4_ip,1_ip,6_ip,2_ip,5_ip,2_ip,3_ip,5_ip,4_ip,1_ip,0_ip,6_ip,7_ip,&
       4_ip,3_ip,7_ip,0_ip,5_ip,2_ip,6_ip,1_ip,6_ip,5_ip,7_ip,4_ip,1_ip,2_ip,0_ip,3_ip,&
       0_ip,7_ip,1_ip,6_ip,3_ip,4_ip,2_ip,5_ip,0_ip,1_ip,7_ip,6_ip,3_ip,2_ip,4_ip,5_ip,&
       2_ip,5_ip,1_ip,6_ip,3_ip,4_ip,0_ip,7_ip,4_ip,7_ip,5_ip,6_ip,3_ip,0_ip,2_ip,1_ip,&
       6_ip,7_ip,1_ip,0_ip,5_ip,4_ip,2_ip,3_ip,2_ip,3_ip,1_ip,0_ip,5_ip,4_ip,6_ip,7_ip,&
       4_ip,3_ip,5_ip,2_ip,7_ip,0_ip,6_ip,1_ip,4_ip,5_ip,3_ip,2_ip,7_ip,6_ip,0_ip,1_ip,&
       6_ip,5_ip,1_ip,2_ip,7_ip,4_ip,0_ip,3_ip,6_ip,1_ip,5_ip,2_ip,7_ip,0_ip,4_ip,3_ip,&
       0_ip,1_ip,3_ip,2_ip,7_ip,6_ip,4_ip,5_ip,2_ip,1_ip,3_ip,0_ip,5_ip,6_ip,4_ip,7_ip,&
       2_ip,5_ip,3_ip,4_ip,1_ip,6_ip,0_ip,7_ip,4_ip,7_ip,3_ip,0_ip,5_ip,6_ip,2_ip,1_ip,&
       6_ip,7_ip,5_ip,4_ip,1_ip,0_ip,2_ip,3_ip,4_ip,5_ip,7_ip,6_ip,3_ip,2_ip,0_ip,1_ip,&
       6_ip,1_ip,7_ip,0_ip,5_ip,2_ip,4_ip,3_ip,2_ip,1_ip,5_ip,6_ip,3_ip,0_ip,4_ip,7_ip/)
  integer(2), parameter,dimension(192) :: typtb3 = (/ &
       1_2,  6_2, 6_2,11_2,11_2,12_2,12_2,14_2, 0_2, 2_2, 2_2, 3_2, 3_2, 4_2, 4_2, 5_2,&
       16_2, 1_2, 1_2,18_2,18_2,19_2,19_2,20_2,12_2,20_2,20_2, 1_2, 1_2,11_2,11_2,18_2,&
       11_2,19_2,19_2,12_2,12_2, 1_2, 1_2,21_2,14_2,18_2,18_2,20_2,20_2,22_2,22_2, 1_2,&
       7_2,  0_2, 0_2, 8_2, 8_2, 9_2, 9_2,10_2, 6_2,16_2,16_2,23_2,23_2,21_2,21_2,22_2,&
       21_2,14_2,14_2, 6_2, 6_2,23_2,23_2,11_2,23_2,12_2,12_2,21_2,21_2, 6_2, 6_2,19_2,&
       22_2,11_2,11_2,14_2,14_2,20_2,20_2, 6_2, 4_2,10_2,10_2, 0_2, 0_2, 3_2, 3_2, 8_2,&
       3_2,  9_2, 9_2, 4_2, 4_2, 0_2, 0_2,13_2,18_2,21_2,21_2,19_2,19_2,16_2,16_2,12_2,&
       5_2,  8_2, 8_2,10_2,10_2,15_2,15_2, 0_2,20_2,23_2,23_2,22_2,22_2,14_2,14_2,16_2,&
       2_2,  7_2, 7_2,17_2,17_2,13_2,13_2,15_2,19_2,22_2,22_2,16_2,16_2,18_2,18_2,23_2,&
       13_2, 5_2, 5_2, 2_2, 2_2,17_2,17_2, 3_2,17_2, 4_2, 4_2,13_2,13_2, 2_2, 2_2, 9_2,&
       15_2, 3_2, 3_2, 5_2, 5_2,10_2,10_2, 2_2, 8_2,13_2,13_2, 9_2, 9_2, 7_2, 7_2, 4_2,&
       10_2,17_2,17_2,15_2,15_2, 5_2, 5_2, 7_2, 9_2,15_2,15_2, 7_2, 7_2, 8_2, 8_2,17_2/)

  interface maths_equalize_arrays
     module procedure maths_equalize_arrays_RP_11,&
          &           maths_equalize_arrays_RP_12,&
          &           maths_equalize_arrays_RP_22,&
          &           maths_equalize_arrays_RP_ndim1,&
          &           maths_equalize_arrays_RP_ndim1_ndim2_22,&
          &           maths_equalize_arrays_RP_ndim1_ndim2_32,&
          &           maths_equalize_arrays_RP_ndim1_ndim2_23,&
          &           maths_equalize_arrays_RP_ndim1_ndim2_33
  end interface maths_equalize_arrays

  interface maths_mapping_coord_to_3d
     module procedure maths_mapping_coord_to_3d_1,&
          &           maths_mapping_coord_to_3d_2
  end interface maths_mapping_coord_to_3d

  interface maths_sfc_1d_to2d3d_tab
     module procedure maths_sfc_d2xy_tab,&
          &           maths_sfc_d2xyz_tab
  end interface maths_sfc_1d_to2d3d_tab

  interface maths_rotation_matrix_2D
     module procedure maths_rotation_matrix_2D_scalar,&
          &           maths_rotation_matrix_2D_vector
  end interface maths_rotation_matrix_2D

  interface maths_angle_vector_2D
     module procedure maths_angle_vector_2D_scalar,&
          &           maths_angle_vector_2D_vector
  end interface maths_angle_vector_2D
  
  public :: maths_MULT_MxV
  public :: maths_MULT_VxMxV
  public :: maths_linear_regression
  public :: maths_weighted_linear_regression
  public :: maths_equalize_arrays
  public :: maths_mapping_coord_to_3d
  public :: maths_in_box
  public :: maths_mapping_3d_to_1d
  public :: maths_mapping_1d_to_3d
  public :: maths_mapping_1d_to_3d_x
  public :: maths_mapping_1d_to_3d_y
  public :: maths_mapping_1d_to_3d_z
  public :: maths_rank_table
  public :: maths_sfc_1d_to2d3d_tab
  public :: maths_sfc_d2xy_tab  
  public :: maths_sfc_d2xyz_tab  
  public :: maths_sfc_ExaIndex
  public :: maths_swap_rp                       ! Swap two real numbers
  public :: maths_rotation_matrix_2D            ! 2D rotation matrix
  public :: maths_rotation_matrix_3D            ! 3D rotation matrix
  public :: maths_angle_vector_2D               ! Angle of a vector with x axis
  public :: maths_cross_product
  public :: maths_normalize_vector
  public :: maths_heaviside                     ! Heaviside
  public :: maths_inverse_permutation           ! Compute an inverse permutation
  public :: maths_inverse                       ! Defined inverse
  !
  ! We cannot have an interface for this,
  ! because with PGI we have qp=rp
  !
  public :: maths_rotation_matrix_2D_qp         ! 2D rotation matrix
  public :: maths_rotation_matrix_3D_qp         ! 3D rotation matrix
  public :: maths_angle_vector_2D_qp            ! Angle of a vector with x axis
  
contains


   !----------------------------------------------------------------------
   !
   !> @author  ???
   !> @date    ???
   !> @brief   Calculate M*V'
   !
   !----------------------------------------------------------------------
   pure function maths_MULT_MxV(M,v,ai) result(r)
      ! -------------------------------
      implicit none
      ! -------------------------------
      integer(ip),        intent(in) :: ai
      real(rp),           intent(in) :: M(ai,ai)
      real(rp),           intent(in) :: v(ai)
      real(rp)                       :: r(ai)
      integer(ip)                    :: ii, jj
      ! -------------------------------
      r(:) = 0.0_rp
      do jj = 1, ai
         do ii = 1, ai
            r(ii) =  r(ii) + M(ii,jj)*v(jj)
         enddo
      enddo
      ! -------------------------------
   end function maths_MULT_MxV

  !----------------------------------------------------------------------
  !
  !> @author  ???
  !> @date    ???
  !> @brief   Calculate V*M*V'
  !> @details Project tensor M on vector V. Moved from sld_stress_model_134.f90
  !
  !----------------------------------------------------------------------
   pure function maths_MULT_VxMxV(a,ai,M,b,bj) result(r)
      ! -------------------------------
      implicit none
      ! -------------------------------
      integer(ip), intent(in)                :: ai, bj
      real(rp), dimension(ai),    intent(in) :: a
      real(rp), dimension(ai,bj), intent(in) :: M
      real(rp), dimension(bj),    intent(in) :: b
      real(rp)                               :: r
      integer(ip)                            :: i, j
      ! -------------------------------
      r = 0.0_rp
      do j = 1, ai
         do i = 1, bj
               r = r + a(i)*M(i,j)*b(j)
         enddo
      enddo
      ! -------------------------------
   end function maths_MULT_VxMxV

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/10/2014
  !> @brief   Equalize some arrays
  !> @details Equalize some arrays
  !
  !----------------------------------------------------------------------

  pure subroutine maths_equalize_arrays_RP_ndim1(ndim1,x_in,x_out)
    integer(ip),       intent(in)  :: ndim1
    real(rp),          intent(in)  :: x_in(*)
    real(rp),          intent(out) :: x_out(*)
    integer(ip)                    :: idim1
    do idim1 = 1,ndim1
       x_out(idim1) = x_in(idim1)
    end do
  end subroutine maths_equalize_arrays_RP_ndim1
  
  pure subroutine maths_equalize_arrays_RP_ndim1_ndim2_22(ndim1,ndim2,x_in,x_out)
    integer(ip),       intent(in)  :: ndim1
    integer(ip),       intent(in)  :: ndim2
    real(rp),          intent(in)  :: x_in(ndim1,ndim2)
    real(rp),          intent(out) :: x_out(ndim1,ndim2)
    integer(ip)                    :: idim1,idim2
    do idim2 = 1,ndim2
       do idim1 = 1,ndim1
          x_out(idim1,idim2) = x_in(idim1,idim2)
       end do
    end do
  end subroutine maths_equalize_arrays_RP_ndim1_ndim2_22
  
  pure subroutine maths_equalize_arrays_RP_ndim1_ndim2_32(ndim1,ndim2,x_in,x_out)
    integer(ip),       intent(in)  :: ndim1
    integer(ip),       intent(in)  :: ndim2
    real(rp),          intent(in)  :: x_in(ndim1,ndim2,*)
    real(rp),          intent(out) :: x_out(ndim1,ndim2)
    integer(ip)                    :: idim1,idim2
    do idim2 = 1,ndim2
       do idim1 = 1,ndim1
          x_out(idim1,idim2) = x_in(idim1,idim2,1)
       end do
    end do
  end subroutine maths_equalize_arrays_RP_ndim1_ndim2_32
  
  pure subroutine maths_equalize_arrays_RP_ndim1_ndim2_23(ndim1,ndim2,x_in,x_out)
    integer(ip),       intent(in)  :: ndim1
    integer(ip),       intent(in)  :: ndim2
    real(rp),          intent(in)  :: x_in(ndim1,ndim2)
    real(rp),          intent(out) :: x_out(ndim1,ndim2,*)
    integer(ip)                    :: idim1,idim2
    do idim2 = 1,ndim2
       do idim1 = 1,ndim1
          x_out(idim1,idim2,1) = x_in(idim1,idim2)
       end do
    end do
  end subroutine maths_equalize_arrays_RP_ndim1_ndim2_23
  
  pure subroutine maths_equalize_arrays_RP_ndim1_ndim2_33(ndim1,ndim2,x_in,x_out)
    integer(ip),       intent(in)  :: ndim1
    integer(ip),       intent(in)  :: ndim2
    real(rp),          intent(in)  :: x_in(ndim1,ndim2,*)
    real(rp),          intent(out) :: x_out(ndim1,ndim2,*)
    integer(ip)                    :: idim1,idim2
    do idim2 = 1,ndim2
       do idim1 = 1,ndim1
          x_out(idim1,idim2,1) = x_in(idim1,idim2,1)
       end do
    end do
  end subroutine maths_equalize_arrays_RP_ndim1_ndim2_33

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/10/2014
  !> @brief   Equalize some arrays
  !> @details Equalize some arrays
  !
  !----------------------------------------------------------------------

  pure subroutine maths_equalize_arrays_RP_11(x_in,x_out)
    real(rp), pointer, intent(in)    :: x_in(:)
    real(rp), pointer, intent(inout) :: x_out(:)
    integer(ip)                      :: idim1
    integer(ip)                      :: ndim1

    ndim1 = min(size(x_in,KIND=ip),size(x_out,KIND=ip))

    do idim1 = 1,ndim1
       x_out(idim1) = x_in(idim1)
    end do

  end subroutine maths_equalize_arrays_RP_11

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/10/2014
  !> @brief   Equalize some arrays
  !> @details Equalize some arrays
  !
  !----------------------------------------------------------------------

  subroutine maths_equalize_arrays_RP_12(x_in,x_out)
    
    real(rp), pointer, intent(in)    :: x_in(:)
    real(rp), pointer, intent(inout) :: x_out(:,:)
    integer(ip)                      :: idim1,idim2,idime
    integer(ip)                      :: ndim1,ndim2

    ndim1 = size(x_out,1,KIND=ip)
    ndim2 = size(x_out,2,KIND=ip)

    if( size(x_in,KIND=ip) /= ndim1*ndim2 ) then
       call runend('maths_equalize_arrays: ARRAYS ARE NOT OF SAME DIMENSIONS')
    else
       idime = 0
       do idim2 = 1,ndim2
          do idim1 = 1,ndim1
             idime = idime + 1
             x_out(idim1,idim2) = x_in(idime)
          end do
       end do
    end if

  end subroutine maths_equalize_arrays_RP_12

  subroutine maths_equalize_arrays_RP_22(x_in,x_out)
    real(rp), pointer, intent(in)    :: x_in(:,:)
    real(rp), pointer, intent(inout) :: x_out(:,:)
    integer(ip)                      :: idim1,idim2
    integer(ip)                      :: ndim1,ndim2

    ndim1 = size(x_out,1,KIND=ip)
    ndim2 = size(x_out,2,KIND=ip)

    if( size(x_in,1,KIND=ip) /= ndim1 .and. size(x_in,2,KIND=ip) /= ndim2 ) then
       call runend('maths_equalize_arrays: ARRAYS ARE NOT OF SAME DIMENSIONS')
    else
       do idim2 = 1,ndim2
          do idim1 = 1,ndim1
             x_out(idim1,idim2) = x_in(idim1,idim2)
          end do
       end do
    end if

  end subroutine maths_equalize_arrays_RP_22

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Mapping (x,y,z) => s
  !> @details Mapping 3D array <=> 1D array
  !>          x: 1 to mx
  !>          y: 1 to my
  !>          z: 1 to mz
  !>          s: 1 to mx*my*mx
  !>
  !>          3D to 1D: s = x + mx*[ (y-1) + (z-1)*my ]
  !>          1D to 3D: x = mod(s-1,mx) + 1
  !>                    y = mod((s-1)/mx,my) + 1
  !>                    z = (s-1)/(mx*my) + 1
  !
  !----------------------------------------------------------------------

  elemental function maths_mapping_3d_to_1d(mx,my,mz,xx,yy,zz)
    integer(ip), intent(in) :: mz,my,mx
    integer(ip), intent(in) :: zz,yy,xx
    integer(ip)             :: maths_mapping_3d_to_1d
    maths_mapping_3d_to_1d = xx+mx*( yy-1 + (zz-1)*my)
  end function maths_mapping_3d_to_1d

  elemental function maths_mapping_1d_to_3d_x(mx,my,mz,ss)
    integer(ip), intent(in) :: mx,my,mz,ss
    integer(ip)             :: maths_mapping_1d_to_3d_x
    maths_mapping_1d_to_3d_x = modulo(ss-1,mx)+1
  end function maths_mapping_1d_to_3d_x

  elemental function maths_mapping_1d_to_3d_y(mx,my,mz,ss)
    integer(ip), intent(in) :: mx,my,mz,ss
    integer(ip)             :: maths_mapping_1d_to_3d_y
    maths_mapping_1d_to_3d_y = modulo((ss-1)/mx,my)+1
  end function maths_mapping_1d_to_3d_y

  elemental function maths_mapping_1d_to_3d_z(mx,my,mz,ss)
    integer(ip), intent(in) :: mx,my,mz,ss
    integer(ip)             :: maths_mapping_1d_to_3d_z
    maths_mapping_1d_to_3d_z = (ss-1) / (mx*my)+1
  end function maths_mapping_1d_to_3d_z

  elemental subroutine maths_mapping_1d_to_3d(ss,mx,my,mz,xx,yy,zz)
    integer(ip), intent(in)  :: ss,mz,my,mx
    integer(ip), intent(out) :: zz,yy,xx

    xx = modulo(ss-1,mx)+1
    yy = modulo((ss-1)/mx,my)+1
    zz = (ss-1) / (mx*my)+1

  end subroutine maths_mapping_1d_to_3d

  !-----------------------------------------------------------------------
  !
  !> @date    27/02/2014
  !> @author  Guillaume Houzeaux
  !> @brief   Find the box in a bin structure
  !> @details Find the box (II,JJ,KK) the point COORD is located in.
  !>          To detect if the point is outside, put the message
  !>          'DETECT OUTSIDE' as an argument. Then,
  !>          II = 0 if points is out of the bin
  !>
  !
  !-----------------------------------------------------------------------

  pure subroutine maths_mapping_coord_to_3d_1(ndime,boxes,comin,comax,coord,xx) !,message)
    
    integer(ip), intent(in)  :: ndime          !< Dimension of problem
    integer(ip), intent(in)  :: boxes(ndime)   !< # boxes in each dimension
    real(rp),    intent(in)  :: comin(ndime)   !< Minimum bin coordinates
    real(rp),    intent(in)  :: comax(ndime)   !< maximum bin coordinates
    real(rp),    intent(in)  :: coord(ndime)   !< Coordinate of the test point
    integer(ip), intent(out) :: xx(3)          !< Box in x direction

    select case ( ndime )

    case ( 1_ip )

       xx(1) = int( ( (coord(1)-comin(1)-epsil) / (comax(1)-comin(1)) )*real(boxes(1),rp), ip ) + 1
       xx(1) = min(max(1_ip,xx(1)),boxes(1))
       xx(2) = 1
       xx(3) = 1

    case ( 2_ip )

       xx(1:2) = int( ( (coord(1:2)-comin(1:2)-epsil) / (comax(1:2)-comin(1:2)) )*real(boxes(1:2),rp), ip ) + 1
       xx(1:2) = min(max(1_ip,xx(1:2)),boxes(1:2))
       xx(3)   = 1

    case ( 3_ip )

       xx(1:3) = int( ( (coord(1:3)-comin(1:3)-epsil) / (comax(1:3)-comin(1:3)) )*real(boxes(1:3),rp), ip ) + 1
       xx(1:3) = min(max(1_ip,xx(1:3)),boxes(1:3))

    end select
    !
    ! Out of bin?
    !
    !detect_outside = .false.
    !if( present(message) ) then
    !   if( trim(message) == 'DETECT OUTSIDE') detect_outside = .true.
    !end if
    !if( detect_outside ) then
    !   if( ii < 1 .or. ii > boxes(1) ) ii = 0
    !   if( ndime >= 2 .and. ( jj < 1 .or. jj > boxes(2) ) ) ii = 0
    !   if( ndime == 3 .and. ( kk < 1 .or. kk > boxes(3) ) ) ii = 0
    !else
    !   ii = min(max(1_ip,ii),boxes(1))
    !   if( ndime >= 2 ) jj = min(max(1_ip,jj),boxes(2))
    !   if( ndime == 3 ) kk = min(max(1_ip,kk),boxes(3))
    !end if

  end subroutine maths_mapping_coord_to_3d_1

  pure subroutine maths_mapping_coord_to_3d_2(ndime,boxes,comin,comax,coord,ii,jj,kk) !,message)
    
    integer(ip), intent(in)  :: ndime          !< Dimension of problem
    integer(ip), intent(in)  :: boxes(ndime)   !< # boxes in each dimension
    real(rp),    intent(in)  :: comin(ndime)   !< Minimum bin coordinates
    real(rp),    intent(in)  :: comax(ndime)   !< maximum bin coordinates
    real(rp),    intent(in)  :: coord(ndime)   !< Coordinate of the test point
    integer(ip), intent(out) :: ii             !< Box in x direction
    integer(ip), intent(out) :: jj             !< Box in y direction
    integer(ip), intent(out) :: kk             !< Box in z direction
    !character(*),intent(in), optional :: message        !< Box in z direction
    !logical(lg)                       :: detect_outside
    integer(ip)                       :: xx(3)

    select case ( ndime )

    case ( 1_ip )

       xx(1) = int( ( (coord(1)-comin(1)-epsil) / (comax(1)-comin(1)) )*real(boxes(1),rp), ip ) + 1
       ii    = min(max(1_ip,xx(1)),boxes(1))
       jj    = 1
       kk    = 1

    case ( 2_ip )

       xx(1:2) = int( ( (coord(1:2)-comin(1:2)-epsil) / (comax(1:2)-comin(1:2)) )*real(boxes(1:2),rp), ip ) + 1
       ii      = min(max(1_ip,xx(1)),boxes(1))
       jj      = min(max(1_ip,xx(2)),boxes(2))
       kk      = 1

    case ( 3_ip )

       xx(1:3) = int( ( (coord(1:3)-comin(1:3)-epsil) / (comax(1:3)-comin(1:3)) )*real(boxes(1:3),rp), ip ) + 1
       ii      = min(max(1_ip,xx(1)),boxes(1))
       jj      = min(max(1_ip,xx(2)),boxes(2))
       kk      = min(max(1_ip,xx(3)),boxes(3))

    end select
    !
    ! Out of bin?
    !
    !detect_outside = .false.
    !if( present(message) ) then
    !   if( trim(message) == 'DETECT OUTSIDE') detect_outside = .true.
    !end if
    !if( detect_outside ) then
    !   if( ii < 1 .or. ii > boxes(1) ) ii = 0
    !   if( ndime >= 2 .and. ( jj < 1 .or. jj > boxes(2) ) ) ii = 0
    !   if( ndime == 3 .and. ( kk < 1 .or. kk > boxes(3) ) ) ii = 0
    !else
    !   ii = min(max(1_ip,ii),boxes(1))
    !   if( ndime >= 2 ) jj = min(max(1_ip,jj),boxes(2))
    !   if( ndime == 3 ) kk = min(max(1_ip,kk),boxes(3))
    !end if

  end subroutine maths_mapping_coord_to_3d_2

  !-----------------------------------------------------------------------
  !
  !> @date    23/05/2014
  !> @author  Guillaume Houzeaux
  !> @brief   Heaviside
  !> @details Heaviside H(x+x_0)
  !
  !-----------------------------------------------------------------------

  pure elemental function maths_heaviside(x,x_0) result(H)

    real(rp), intent(in) :: x
    real(rp), intent(in) :: x_0
    real(rp)             :: H
    
    if( x > x_0 ) then
       H = 1.0_rp
    else
       H = 0.0_rp
    end if
    
  end function maths_heaviside
  
  !-----------------------------------------------------------------------
  !
  !> @date    23/05/2014
  !> @author  Guillaume Houzeaux
  !> @brief   If in a box
  !> @details If in a box
  !
  !-----------------------------------------------------------------------

  pure function maths_in_box(ndime,xx,box_comin,box_comax)
    
    integer(ip), intent(in) :: ndime
    real(rp),    intent(in) :: xx(ndime)
    real(rp),    intent(in) :: box_comin(ndime)
    real(rp),    intent(in) :: box_comax(ndime)
    logical(lg)             :: maths_in_box
    integer(ip)             :: idime

    maths_in_box = .true.
    do idime = 1,ndime
       if( xx(idime) < box_comin(idime) .or. xx(idime) > box_comax(idime) ) then
          maths_in_box = .false.
          return
       end if
    end do

  end function maths_in_box

  !----------------------------------------------------------------------
  !>
  !> @author  Herbert Owen
  !> @date    14/06/2016
  !> @brief   obtain rank table from index table
  !> @details See Figure 8.4.1 from numerical recipes (sorting)
  !>          adapted from NR - Given indx(1:n) as output from the routine indexx, this routine returns an array irank(1:n),
  !>          the corresponding table of ranks.
  !>          In order to obtain indx you have to give maths_heap_sort   on input the vector 1,2,3,....n to ivo1
  !>
  !----------------------------------------------------------------------

  pure subroutine maths_rank_table(n,indx,irank)
    
    integer(ip), intent(in)            :: n               !< Dimension
    integer(ip), intent(in)            :: indx(n)         !< index table
    integer(ip), intent(out)           :: irank(n)        !< rank table
    integer(ip)                        :: j

    do  j=1,n
       irank(indx(j))=j
    end do
    
  end subroutine maths_rank_table

  !-----------------------------------------------------------------------
  !
  !> @date    28/10/2014
  !> @author  Ricard Borrell
  !> @brief   Converts 1D Hilbert coordinate to 2D cartesian coordinates
  !> @details Kernel based on bit operations
  !           n: number of boxes in each direction
  !           d: the hilbert coordinate (1<=d<=n*n)
  !           x,y: 1<=x,y<=n
  !
  !-----------------------------------------------------------------------
  
  elemental subroutine maths_sfc_d2xy_bit(n,d,x,y)

    integer(ip), intent(in)  :: n,d
    integer(ip), intent(out) :: x,y
    integer(ip)              :: rx,ry,s,t

    x = 0_ip
    y = 0_ip
    t = d-1_ip
    s = 1_ip

    do while( s<n )
       rx = mod(t/2,2_ip)
       if( rx==0 ) then
          ry = mod(t,2_ip)
       else
          ry = mod(ieor(t,rx),2_ip)
       endif
       call maths_sfc_rot_bit(s,x,y,rx,ry)
       x = x+s*rx
       y = y+s*ry
       t = t/4
       s=s*2
    end do

    x=x+1_ip
    y=y+1_ip

  end subroutine maths_sfc_d2xy_bit
  
  !-----------------------------------------------------------------------
  !
  !> @date    28/10/2014
  !> @author  Ricard Borrell
  !> @brief   Rotates and flips a quadrant appropritely
  !
  !---------------------------------------------------------------------
  
  elemental subroutine maths_sfc_rot_bit(n,x,y,rx,ry)

    integer(ip), intent(in)    :: n,rx,ry
    integer(ip), intent(inout) :: x,y
    integer(ip)                :: t

    if(ry==0) then

       !Reflect
       if( rx == 1_ip ) then
          x = n-1_ip-x
          y = n-1_ip-y
       endif

       !Flip
       t = x
       x = y
       y = t
    endif

  end subroutine maths_sfc_rot_bit
  
  !-----------------------------------------------------------------------
  !
  !> @date    28/10/2014
  !> @author  Ricard Borrell
  !> @brief   Converts 1D Hilbert coordinate to 2D cartesian coordinates
  !> @details Kernel based on the table of recursive divisions
  !           n: number of boxes in each direction
  !           d: the hilbert coordinate (1<=d<=n*n)
  !           x,y: 1<=x,y<=n
  !
  !-----------------------------------------------------------------------
  
  subroutine maths_sfc_d2xy_tab(n,d,x,y,typ)

    integer(ip), intent(in)             :: n,d
    integer(ip), intent(out)            :: x,y
    integer(2),  optional,intent(inout) :: typ
    integer(ip)                         :: ntot,naxi
    integer(ip)                         :: pos1,pos2,t
    integer(2)                          :: typ1

    x    = 0_ip
    y    = 0_ip
    t    = d-1_ip
    typ1 = 0_2
    if(present(typ)) then
       typ1 = typ
    endif

    ntot = n*n
    naxi = n
    do while( ntot> 1_ip )
       naxi = naxi/2_ip
       pos1 = (t*4_ip)/ntot
       pos2  = int(postb2(int(typ1,ip)*4_ip+pos1+1_ip),KIND=ip)
       typ1  = typtb2(typ1*4+pos1+1)
       if(pos2==1_ip .or. pos2==3_ip) then
          x = x + naxi
       endif
       if(pos2>1_ip) then
          y = y + naxi
       endif
       t = t-(pos1*ntot)/4_ip
       ntot = ntot/4_ip
    end do

    x=x+1_ip
    y=y+1_ip
    if(present(typ)) typ = typ1

  end subroutine maths_sfc_d2xy_tab

  !-----------------------------------------------------------------------
  !
  !> @date    3/11/2014
  !> @author  Ricard Borrell
  !> @brief   Converts 1D Hilbert coordinate to 3D cartesian coordinates
  !> @details Kernel based on the table of recursive divisions
  !           n: number of boxes in each direction
  !           d: the hilbert coordinate (1<=d<=n*n*n)
  !           x,y,z: 1<=x,y,z<=n
  !
  !-----------------------------------------------------------------------
  
  subroutine maths_sfc_d2xyz_tab(n,d,x,y,z,typ)

    integer(ip), intent(in)             :: n,d
    integer(ip), intent(out)            :: x,y,z
    integer(2), optional,intent(inout)  :: typ
    integer(ip)                         :: ntot,naxi,t
    integer(ip)                         :: pos1,pos2,auxp
    integer(2)                          :: typ1

    x   = 0_ip
    y   = 0_ip
    z   = 0_ip
    t   = d-1_ip
    typ1 = 0_2
    if(present(typ)) typ1 = typ

    ntot = n*n*n
    naxi = n
    do while( ntot> 1_ip )
       naxi = naxi/2_ip
       pos1 = (t*8_ip)/ntot
       pos2  = postb3(int(typ1,ip)*8_ip+pos1+1_ip)
       typ1  = typtb3(int(typ1,ip)*8_ip+pos1+1_ip)
       auxp=mod(pos2,4_ip)
       if(auxp == 1_ip .or. auxp == 3_ip) x = x + naxi
       if(auxp > 1_ip) y = y + naxi
       if(pos2 > 3_ip) z = z + naxi
       t = t-(pos1*ntot)/8_ip
       ntot = ntot/8_ip
    end do

    x=x+1_ip
    y=y+1_ip
    z=z+1_ip
    if(present(typ)) typ = typ1

  end subroutine maths_sfc_d2xyz_tab
  
  !-----------------------------------------------------------------------
  !
  !> @date    7/10/2020
  !> @author  Ricard Borrell
  !> @brief   Calculate Exa SFC-index
  !> @details Given the bounding box (min_coord-max_coord)
  !>          the SFC orientation (typ) evaluate the SFC index (ind)
  !>          corresponding to the bin containing coord, for an SFC of 
  !>          level 20 (with 20**sdim bins)
  !
  !-----------------------------------------------------------------------
  subroutine maths_sfc_ExaIndex(coord,min_coord,max_coord,sdim,ind,maxIdx,typ)

     real(rp),                intent(in)    :: coord(3)
     real(rp),                intent(in)    :: min_coord(3)
     real(rp),                intent(in)    :: max_coord(3)
     integer(ip),             intent(in)    :: sdim
     integer(8),              intent(out)   :: ind      
     integer(8),              intent(out)   :: maxIdx     
     integer(2),  optional,   intent(inout) :: typ

     integer(ip)          :: ilev,idcar,idsfc,ii
     real(rp)             :: rcoord(3)
     integer(ip)          :: cacoor(3)
     integer(2)           :: typ1
     integer(ip)          :: lev
     integer(8)           :: auxi
     integer(8)           :: base

     typ1 = 0_2
     if(present(typ)) typ1 = typ
     do ii=1,sdim
        rcoord(ii)=(coord(ii)-min_coord(ii))/(max_coord(ii)-min_coord(ii))
     enddo
     ind  = 0_8
     lev  = 20
     base = 2_8**lev

     if(sdim==2) then

        maxIdx=base*base
        auxi = maxIdx
        do ilev=1, lev
           auxi = auxi/4_8
           do ii=1,sdim
              cacoor(ii)=1_ip
              if(rcoord(ii) > 0.5_rp) then
                 rcoord(ii) = rcoord(ii)-0.5_rp
                 cacoor(ii)=2_ip
              endif
              rcoord(ii)=rcoord(ii)/0.5_rp
           enddo
           idcar = (cacoor(2)-1)*2+cacoor(1)
           idsfc = postc2(int(typ1,ip)*4_ip+idcar)+1
           typ1  = typtb2(int(typ1,ip)*4_ip+idsfc)
           ind   = ind + (idsfc-1) * auxi
        enddo

     else if (sdim ==3) then

        maxIdx=base*base*base
        auxi = maxIdx
        do ilev=1, lev
           auxi = auxi/8_8
           do ii=1,sdim
              cacoor(ii)=1_ip
              if(rcoord(ii) > 0.5_rp) then
                 rcoord(ii) = rcoord(ii)-0.5_rp
                 cacoor(ii)=2_ip
              endif
              rcoord(ii)=rcoord(ii)/0.5_rp
           enddo
           idcar = (cacoor(3)-1)*4+(cacoor(2)-1)*2+cacoor(1)
           idsfc = postc3(int(typ1,ip)*8_ip+idcar)+1
           typ1  = typtb3(int(typ1,ip)*8_ip+idsfc)
           ind   = ind + (idsfc-1) * auxi
        enddo
     else
        call runend("aths_sfc_index_relative: sdim must be 1 or 2")
     endif
     ind = ind + 1_8
     if(present(typ)) typ = typ1

  end subroutine maths_sfc_ExaIndex

  !-----------------------------------------------------------------------
  !
  !> @date    3/11/2014
  !> @author  Ricard Borrell
  !> @brief   Converts 1D lexicografical order to 3D cartesian coordinates
  !> @details n: number of boxes in each direction
  !>          d: position in lexicografical order (1<=d<=n*n*n)
  !>          x,y: 1<=x,y,z<=n
  !
  !-----------------------------------------------------------------------
  
  subroutine maths_sfc_d2xyz_lex(n,d,x,y,z)
    
    integer(ip), intent(in)  :: n,d
    integer(ip), intent(out) :: x,y,z
    x=n
    y=n
    z=(d-1_ip)/(n*n)
    if(d-z*n*n >0_ip) then
       y=((d-z*n*n)-1_ip)/n
       x=d-z*n*n-y*n
    endif
    z=z+1_ip
    y=y+1_ip

  end subroutine maths_sfc_d2xyz_lex
  

  !-----------------------------------------------------------------------
  !> 
  !> @author  Adria Quintanas-Corominas
  !> @date    2020-05-20
  !> @brief   Swap to real values
  !> @details Swap to real values
  !> 
  !-----------------------------------------------------------------------

  elemental subroutine maths_swap_rp(lhs,rhs)

     real(rp), intent(inout) :: lhs,rhs
     real(rp)                :: temp

     temp = lhs
     lhs  = rhs
     rhs  = temp

  end subroutine maths_swap_rp

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-19
  !> @brief   Angle of a vector with respect to x
  !> @details Angle of a vector with respect to x using the cos
  !> 
  !-----------------------------------------------------------------------

  pure function maths_angle_vector_2D_scalar(vect) result(phi)

    real(rp), intent(in) :: vect(2)
    real(rp)             :: phi
    real(rp)             :: vectn
    real(rp)             :: vectu(2)
    real(rp)             :: cos_phi
    real(rp)             :: sin_phi

    vectn = sqrt(dot_product(vect(1:2),vect(1:2)))
    if( vectn /= 0.0_rp ) then
       vectu(1:2) = vect(1:2) / vectn
       cos_phi    = max(-1.0_rp,min(1.0_rp,vectu(1)))
       sin_phi    = max(-1.0_rp,min(1.0_rp,vectu(2)))
       phi        = acos(cos_phi)
       if( sin_phi < 0.0_rp ) phi = 2.0_rp*pi-phi
    else
       phi     = 0.0_rp       
    end if
       
  end function maths_angle_vector_2D_scalar
  
  pure function maths_angle_vector_2D_vector(vect) result(phi)

    real(rp), intent(in) :: vect(:,:)
    real(rp)             :: phi(size(vect,1))
    real(rp)             :: vectn
    real(rp)             :: vectu(2)
    real(rp)             :: cos_phi
    real(rp)             :: sin_phi
    integer(ip)          :: ii

    
    do ii = 1,size(vect,1)
       vectn = sqrt(dot_product(vect(ii,1:2),vect(ii,1:2)))
       if( vectn /= 0.0_rp ) then
          vectu(1:2) = vect(ii,1:2) / vectn
          cos_phi    = max(-1.0_rp,min(1.0_rp,vectu(1)))
          sin_phi    = max(-1.0_rp,min(1.0_rp,vectu(2)))
          phi(ii)    = acos(cos_phi)
          if( sin_phi < 0.0_rp ) phi(ii) = 2.0_rp*pi-phi(ii)
       else
          phi(ii) = 0.0_rp       
       end if
    end do
    
  end function maths_angle_vector_2D_vector
  
  pure function maths_angle_vector_2D_qp(vect) result(phi)

    real(qp), intent(in) :: vect(2)
    real(qp)             :: phi
    real(qp)             :: vectn
    real(qp)             :: vectu(2)
    real(qp)             :: cos_phi
    real(qp)             :: sin_phi

    vectn = sqrt(dot_product(vect(1:2),vect(1:2)))
    if( vectn /= 0.0_qp ) then
       vectu(1:2) = vect(1:2) / vectn
       cos_phi    = max(-1.0_qp,min(1.0_qp,vectu(1)))
       sin_phi    = max(-1.0_qp,min(1.0_qp,vectu(2)))
       phi        = acos(cos_phi)
       if( sin_phi < 0.0_qp ) phi = 2.0_qp*pi-phi
    else
       phi     = 0.0_qp       
    end if
       
  end function maths_angle_vector_2D_qp
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-19
  !> @brief   Rotation matrix in 2D
  !> @details Rotation matrix in 2D
  !> 
  !-----------------------------------------------------------------------

  pure function maths_rotation_matrix_2D_scalar(phi) result(R)

    real(rp), intent(in) :: phi
    real(rp)             :: R(2,2)

    R(1,1) =  cos(phi)
    R(1,2) = -sin(phi)
    R(2,1) =  sin(phi)
    R(2,2) =  cos(phi)
    
  end function maths_rotation_matrix_2D_scalar
  
  pure function maths_rotation_matrix_2D_vector(phi) result(R)

    real(rp), intent(in) :: phi(:)
    real(rp)             :: R(size(phi),2,2)
    integer(ip)          :: ii

    do ii = 1,size(phi)
       R(ii,1,1) =  cos(phi(ii))
       R(ii,1,2) = -sin(phi(ii))
       R(ii,2,1) =  sin(phi(ii))
       R(ii,2,2) =  cos(phi(ii))
    end do
    
  end function maths_rotation_matrix_2D_vector
  
  pure function maths_rotation_matrix_2D_qp(phi) result(R)

    real(qp), intent(in) :: phi
    real(qp)             :: R(2,2)

    R(1,1) =  cos(phi)
    R(1,2) = -sin(phi)
    R(2,1) =  sin(phi)
    R(2,2) =  cos(phi)
    
  end function maths_rotation_matrix_2D_qp
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-19
  !> @brief   Rotation matrix in 3D
  !> @details Rotation matrix in 3D
  !> 
  !-----------------------------------------------------------------------

  pure function maths_rotation_matrix_3D(phi,axis,v) result(R)

    real(rp),              intent(in) :: phi
    integer(ip), optional, intent(in) :: axis
    real(rp),    optional, intent(in) :: v(3)
    real(rp)                          :: R(3,3)
    real(rp)                          :: cos1
    real(rp)                          :: cosp
    real(rp)                          :: sinp
    
    R = 0.0_rp

    if( present(v) ) then

       cosp   = cos(phi)
       sinp   = sin(phi)
       cos1   = 1.0_rp - cos(phi)
       R(1,1) = v(1)*v(1)*cos1 + cosp 
       R(1,2) = v(1)*v(2)*cos1 - v(3)*sinp
       R(1,3) = v(1)*v(3)*cos1 + v(2)*sinp
       R(2,1) = v(1)*v(2)*cos1 + v(3)*sinp
       R(2,2) = v(2)*v(2)*cos1 + cosp
       R(2,3) = v(2)*v(3)*cos1 - v(1)*sinp
       R(3,1) = v(1)*v(3)*cos1 - v(2)*sinp
       R(3,2) = v(2)*v(3)*cos1 + v(1)*sinp
       R(3,3) = v(3)*v(3)*cos1 + cosp
       
    else if( present(axis) ) then
       if(      axis == 1 ) then
          R(1,1) =  1.0_rp
          R(2,2) =  cos(phi)
          R(2,3) = -sin(phi)
          R(3,2) =  sin(phi)
          R(3,3) =  cos(phi)
       else if( axis == 2 ) then
          R(2,2) =  1.0_rp
          R(1,1) =  cos(phi)
          R(1,3) =  sin(phi)
          R(3,1) = -sin(phi)
          R(3,3) =  cos(phi)
       else if( axis == 3 ) then
          R(3,3) =  1.0_rp
          R(1,1) =  cos(phi)
          R(1,2) = -sin(phi)
          R(2,1) =  sin(phi)
          R(2,2) =  cos(phi)
       end if
    end if
    
  end function maths_rotation_matrix_3D

  pure function maths_rotation_matrix_3D_qp(phi,axis) result(R)

    real(qp),    intent(in) :: phi
    integer(ip), intent(in) :: axis
    real(qp)                :: R(3,3)

    R = 0.0_qp
    
    if(      axis == 1 ) then
       R(1,1) =  1.0_qp
       R(2,2) =  cos(phi)
       R(2,3) = -sin(phi)
       R(3,2) =  sin(phi)
       R(3,3) =  cos(phi)
    else if( axis == 2 ) then
       R(2,2) =  1.0_qp
       R(1,1) =  cos(phi)
       R(1,3) =  sin(phi)
       R(3,1) = -sin(phi)
       R(3,3) =  cos(phi)
    else if( axis == 3 ) then
       R(3,3) =  1.0_qp
       R(1,1) =  cos(phi)
       R(1,2) = -sin(phi)
       R(2,1) =  sin(phi)
       R(2,2) =  cos(phi)
    end if

 end function maths_rotation_matrix_3D_qp

 !-----------------------------------------------------------------------
 !> 
 !> @author  borrell
 !> @date    29-06-2020
 !> @brief   Weighted Linear regression
 !> @details Perform a weighted linear regression
 !> 
 !-----------------------------------------------------------------------
 subroutine maths_weighted_linear_regression(xx,yy,N,a,b,ww_,windo_)

    real(rp),      pointer,             intent(in)     :: xx(:)
    real(rp),      pointer,             intent(in)     :: yy(:)
    integer(ip),                        intent(in)     :: N
    real(rp),                           intent(out)    :: a
    real(rp),                           intent(out)    :: b
    real(rp),      pointer, optional,   intent(in)     :: ww_(:)
    integer(ip),            optional,   intent(in)     :: windo_

    real(rp),      pointer :: ww(:)
    integer(ip)            :: windo
    real(rp)               :: sx, sy, sxx, sxy, rNW
    real(rp)               :: numer,denom
    integer(ip)            :: ii
    integer(ip)            :: fi,ind1

    !
    ! Process optional arguments
    !
    if(present(windo_)) then
       windo = windo_
    else
       windo = N
    endif

    nullify(ww)
    if(present(ww_)) then
       ww => ww_
    else
       call memory_alloca(memor,'ww',"maths_weighted_linear_regression",ww,windo)
    endif

    !
    ! Evaluate regression coeficients
    !
    sx  = 0.0_rp
    sy  = 0.0_rp
    sxx = 0.0_rp
    sxy = 0.0_rp

    fi = min(N,windo)
    ind1 = mod(N,windo)+1

    rNW  = 0.0_rp
    do ii = 1,fi
       sx  = sx  + ww(ii)*xx(ii)
       sy  = sy  + ww(ii)*yy(ii)
       rNW = rNW + ww(ii)
    enddo
    sx = sx/rNW
    sy = sy/rNW

    numer = 0.0_rp
    denom = 0.0_rp
    do ii = 1,fi
       numer = numer + ww(ii)*(yy(ii)-sy)*(xx(ii)-sx)
       denom = denom + ww(ii)*(xx(ii)-sx)*(xx(ii)-sx)
    enddo

    if(abs(denom)<1.e-8) then
       if(denom < 0) then
          denom = denom - 1.e-8
       else  
          denom = denom + 1.e-8
       endif
    endif

    a = numer / denom
    b = (sy     -  a*sx)

    !
    ! Memory
    !
    if(.not. present(ww_)) then
       call memory_deallo(memor,'ww',"maths_weighted_linear_regression",ww)
    endif

 end subroutine maths_weighted_linear_regression

 !-----------------------------------------------------------------------
 !> 
 !> @author  borrell
 !> @date    29-06-2020
 !> @brief   Linear regression
 !> @details Perform a linear regression
 !> 
 !-----------------------------------------------------------------------
 subroutine maths_linear_regression(xx,yy,N,a,b,windo)

    real(rp), pointer,    intent(in)     :: xx(:)
    real(rp), pointer,    intent(in)     :: yy(:)
    integer(ip),          intent(in)     :: N
    real(rp),             intent(out)    :: a
    real(rp),             intent(out)    :: b
    integer(ip), optional,intent(in)     :: windo

    real(rp)  :: sx, sy, sxx, sxy, rNW,denom
    integer(ip)   :: ii
    integer(ip)   :: fi,ind1

    sx  = 0.0_rp
    sy  = 0.0_rp
    sxx = 0.0_rp
    sxy = 0.0_rp

    fi = min(N,windo)
    ind1 = mod(N,windo)+1
    rNW  = min(windo-1,fi)

    do ii = 1,fi

       if(ii/=ind1) then
          sx  = sx  + xx(ii)
          sy  = sy  + yy(ii)
          sxx = sxx + xx(ii) * xx(ii)
          sxy = sxy + xx(ii) * yy(ii)  
       endif
    enddo
    denom = rNW*sxx-sx*sx
    if(abs(denom)<1.e-8) then
       if(denom < 0) then
          denom = denom - 1.e-8
       else  
          denom = denom + 1.e-8
       endif
    endif

    a = (rNW*sxy - sx*sy) / denom
    b = (sy     -  a*sx) / rNW

 end subroutine maths_linear_regression
    
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

  pure function maths_cross_product(v1,v2,n) result(v3)

    integer(ip), intent(in)  :: n
    real(rp),    intent(in)  :: v1(n)
    real(rp),    intent(in)  :: v2(n)
    real(rp)                 :: v3(3)

    if( n == 2 ) then
       v3(1:2) = 0.0_rp
       v3(3)   = v1(1)*v2(2)-v1(2)*v2(1)
    else if( n == 3 ) then
       v3(1)   = v1(2)*v2(3)-v1(3)*v2(2)
       v3(2)   = v1(3)*v2(1)-v1(1)*v2(3)
       v3(3)   = v1(1)*v2(2)-v1(2)*v2(1)
    end if

  end function maths_cross_product

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-04-09
  !> @brief   ???
  !> @details This routine computes the length of vector V and converts it to
  !> a unit one
  !>
  !-----------------------------------------------------------------------

  pure subroutine maths_normalize_vector(n,v,norm)

    integer(ip), intent(in)            :: n
    real(rp),    intent(inout)         :: v(n)
    real(rp),    intent(out), optional :: norm
    real(rp)                           :: rmod

    rmod = sqrt(dot_product(v(1:n),v(1:n)))
    if( rmod > epsilon(1.0_rp) ) v = v / rmod
    if( present(norm) ) norm = rmod

  end subroutine maths_normalize_vector

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-04-09
  !> @brief   Permutation
  !> @details Compute an inverse permutation
  !>
  !-----------------------------------------------------------------------

  subroutine maths_inverse_permutation(permu,invpe,MEMORY_COUNTER) 

    integer(ip),           pointer, intent(in)    :: permu(:)
    integer(ip),           pointer, intent(inout) :: invpe(:)
    integer(8),  optional,          intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                                    :: memor_loc(2)
    integer(ip)                                   :: nn,ii,kk

    if( associated(permu) ) then
       memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
       nn = maxval(permu)
       call memory_alloca(memor_loc,'MASK','mod_maths_basic',invpe,nn)
       do ii = 1,size(permu)
          kk        = permu(ii)
          if( kk /= 0 ) invpe(kk) = ii
       end do
       if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc
    end if

  end subroutine maths_inverse_permutation
  
  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-04-09
  !> @brief   Defined inverse
  !> @details Defined inverse
  !>
  !-----------------------------------------------------------------------

  elemental function maths_inverse(x,eps) result(x_inv)
    real(rp), intent(in) :: x
    real(rp), intent(in) :: eps
    real(rp)             :: x_inv

    if( abs(x) <= eps ) then
       x_inv = 0.0_rp
    else
       x_inv = 1.0_rp/x
    end if
    
  end function maths_inverse
  
end module mod_maths_basic
!> @}
