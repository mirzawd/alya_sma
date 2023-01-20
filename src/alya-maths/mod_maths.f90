!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @defgroup Maths_Toolbox
!> Toolbox for mathematical functions and subroutines
!> @{
!> @name    ToolBox for mathematics operations
!> @file    mod_maths.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for maths
!> @details ToolBox for maths
!
!-----------------------------------------------------------------------

module mod_maths
  
  use def_maths
  use mod_maths_basic
  use mod_maths_sort
  use mod_maths_solver
  use mod_maths_geometry
  use mod_maths_arrays
  implicit none

  private

  interface maths_matrix_vector_multiplication
     module procedure maths_matrix_vector_multiplication_0,&
          &           maths_matrix_vector_multiplication_1
  end interface maths_matrix_vector_multiplication

  interface maths_array_permutation
     module procedure maths_array_permutation_IP_1,&
          &           maths_array_permutation_IP_2
  end interface maths_array_permutation
  
  interface maths_sfc_par_part_2
     module procedure maths_sfc_par_part_2_ip,&
          &           maths_sfc_par_part_2_rp
  end interface maths_sfc_par_part_2
     
  interface maths_sfc_part_2
     module procedure maths_sfc_part_2_ip,&
          &           maths_sfc_part_2_rp
  end interface maths_sfc_part_2

  public :: maths_local_orthonormal_basis
  public :: maths_vector_to_new_basis
  public :: maths_vector_from_new_basis
  public :: maths_array_permutation
  public :: maths_time_unit
  public :: maths_outer_product
  public :: maths_norm2                         ! L2 norm of an array
  public :: maths_backsu
  public :: maths_matrix_multiplication
  public :: maths_matrix_vector_multiplication
  public :: maths_matrix_transpose_vector_multiplication
  public :: maths_economy_qrhousehold
  public :: maths_sfc_part
  public :: maths_sfc_par_part
  public :: maths_sfc_par_part_2
  public :: maths_list_different_elements       ! List the different elements of an array
  public :: maths_unit                          ! Scaling and basic unit of a value
  public :: maths_vectorial_product             ! Vectorial product
  public :: maths_normalize_vector              ! Normalize a vector
  public :: maths_merge_ordered_lists           ! Merge ordered list
  public :: maths_quadratic_equation            ! Solve a quadratic equation
  public :: maths_Szudzik_pairing_function
  public :: maths_distribute_iterations         ! Dsitribute iteration among a list of candidates
  public :: maths_cross_product                 ! Cross products
  !
  ! From maths_geometry
  !
  public :: maths_circumradius_tetrahedron      ! Circumradius of a tet
  public :: maths_circumradius_triangle         ! Circumradius of a triangle
  public :: maths_circumradius_quadrilateral    ! Circumradius of a quadrilateral
  public :: maths_area_triangle                 ! Area of a triangle
  public :: maths_volume_tetrahedron            ! Volume of a tet
  public :: maths_normal_to_triangle            ! Compute the normal to a triangle
  !
  ! From maths_basic
  !
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
  public :: maths_swap_rp                       ! Swap two real numbers
  public :: maths_rotation_matrix_2D            ! 2D rotation matrix
  public :: maths_rotation_matrix_3D            ! 3D rotation matrix
  public :: maths_angle_vector_2D               ! Angle of a vector with x axis
  !
  ! From maths_solver
  !
  public :: maths_invert_matrix                 ! Matrix inversion
  public :: maths_conjugate_gradient            ! Maths CG solver
  public :: maths_direct                        ! Maths direct solver
  public :: maths_schur_complement              ! Schur complement
  public :: maths_eigen_symmetric_matrix        ! Bridgr to 2D and 3D
  public :: maths_eigen_3x3_symmetric_matrix    ! Eigenvalues and vectors for real symmetric 3x3 system
  public :: maths_eigen_2x2_symmetric_matrix    ! Eigenvalues and vectors for real symmetric 2x2 system
  public :: maths_solve_overdetermined_system   ! Overdetermined system
  !
  ! From maths_sort
  !
  public :: maths_heap_sort
  public :: maths_quick_sort
  public :: maths_geometrical_sort_using_coordinates
  public :: maths_geometrical_sort_using_sfc
  !
  ! From maths_arrays
  !
  public :: maths_maxloc_nonzero                ! Last non-zero position
  public :: maths_findloc                       ! First position of value in array

contains

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    04/12/2015
  !> @brief   Rotate a vector
  !> @details Rotate a vector to a new basis, global -> local
  !>
  !>          Let ei' be the cartesian basis
  !>          e1' = [ 1 0 0 ]
  !>          e2' = [ 0 1 0 ]
  !>          e3' = [ 0 0 1 ]
  !>          and ei the new basis where
  !>          e1  = [ basis(1,1) basis(2,1) basis(3,1) ]
  !>          e2  = [ basis(1,2) basis(2,2) basis(3,2) ]
  !>          e3  = [ basis(1,3) basis(2,3) basis(3,3) ]
  !>
  !>          The base change is the following:
  !>
  !>          +-  -+   +-                                -+ +-   -+
  !>          | x1 |   | basis(1,1) basis(2,1) basis(3,1) | | x1' |
  !>          | x2 | = | basis(1,2) basis(2,2) basis(3,2) | | x2' |
  !>          | x3 |   | basis(1,3) basis(2,3) basis(3,3) | | x3' |
  !>          +-  -+   +-                                -+ +-   -+
  !>
  !----------------------------------------------------------------------

  subroutine maths_vector_to_new_basis(ndime,basis,vv,ww)

    integer(ip), intent(in)              :: ndime               !< Dimension
    real(rp),    intent(in)              :: basis(ndime,ndime)  !< Local basis: normal is BASIS(1:NDIME,1)
    real(rp),    intent(inout)           :: vv(ndime)           !< Vector 1
    real(rp),    intent(inout), optional :: ww(ndime)           !< Vector 2
    real(rp)                             :: vv_tmp(ndime)
    real(rp)                             :: ww_tmp(ndime)
    integer(ip)                          :: ii,jj

    if( present(ww) ) then
       vv_tmp(1:ndime) = vv(1:ndime)
       ww_tmp(1:ndime) = ww(1:ndime)
       do ii = 1,ndime
          vv(ii) = 0.0_rp
          ww(ii) = 0.0_rp
          do jj = 1,ndime
             vv(ii) = vv(ii) + basis(jj,ii) * vv_tmp(jj)
             ww(ii) = ww(ii) + basis(jj,ii) * ww_tmp(jj)
          end do
       end do
    else
       vv_tmp(1:ndime) = vv(1:ndime)
       do ii = 1,ndime
          vv(ii) = 0.0_rp
          do jj = 1,ndime
             vv(ii) = vv(ii) + basis(jj,ii) * vv_tmp(jj)
          end do
       end do
    end if

  end subroutine maths_vector_to_new_basis

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    04/12/2015
  !> @brief   Rotate a vector
  !> @details Rotate a vector from a new basis: local -> global
  !>
  !>                         n1         t1         t2
  !>          +-   -+   +-                                -+ +-  -+
  !>          | x1' |   | basis(1,1) basis(1,2) basis(1,3) | | x1 |
  !>          | x2' | = | basis(2,1) basis(2,2) basis(2,3) | | x2 |
  !>          | x3' |   | basis(3,1) basis(3,2) basis(3,3) | | x3 |
  !>          +-   -+   +-                                -+ +-  -+
  !>
  !----------------------------------------------------------------------

  subroutine maths_vector_from_new_basis(ndime,basis,vv,ww)

    integer(ip), intent(in)              :: ndime               !< Dimension
    real(rp),    intent(in)              :: basis(ndime,ndime)  !< Local basis: normal is BASIS(1:NDIME,1)
    real(rp),    intent(inout)           :: vv(ndime)           !< Vector 1
    real(rp),    intent(inout), optional :: ww(ndime)           !< Vector 2
    real(rp)                             :: vv_tmp(ndime)
    real(rp)                             :: ww_tmp(ndime)
    integer(ip)                          :: ii,jj

    if( present(ww) ) then
       vv_tmp(1:ndime) = vv(1:ndime)
       ww_tmp(1:ndime) = ww(1:ndime)
       do ii = 1,ndime
          vv(ii) = 0.0_rp
          ww(ii) = 0.0_rp
          do jj = 1,ndime
             vv(ii) = vv(ii) + basis(ii,jj) * vv_tmp(jj)
             ww(ii) = ww(ii) + basis(ii,jj) * ww_tmp(jj)
          end do
       end do
    else
       vv_tmp(1:ndime) = vv(1:ndime)
       do ii = 1,ndime
          vv(ii) = 0.0_rp
          do jj = 1,ndime
             vv(ii) = vv(ii) + basis(ii,jj) * vv_tmp(jj)
          end do
       end do
    end if

  end subroutine maths_vector_from_new_basis

  !-----------------------------------------------------------------------
  !> @author  Guillaume Houzeaux
  !> @date    15/02/2016
  !> @brief   Permute arrays
  !> @details Permute arrays
  !-----------------------------------------------------------------------

  subroutine maths_array_permutation_IP_1(invpr,array)
    
    integer(ip), intent(in) ,   pointer :: invpr(:)   !< NEW  = INVPR(OLD)
    integer(ip), intent(inout), pointer :: array(:)   !< Array
    integer(ip)                         :: ii,jj,nn
    integer(ip), allocatable            :: array_tmp(:)

    if( .not. associated(invpr) ) return
    nn = size(invpr,KIND=ip)
    if( nn <= 0 ) return
    if( nn /= size(array,KIND=ip) ) call runend('MATHS_ARRAY_PERMUTATION: WROND DIMENSIONS')

    allocate(array_tmp(nn))
    array_tmp = array
    do ii = 1,nn
       jj        = invpr(ii)
       array(jj) = array_tmp(ii)
    end do
    deallocate(array_tmp)

  end subroutine maths_array_permutation_IP_1

  subroutine maths_array_permutation_IP_2(invpr,array)
    integer(ip), intent(in) ,   pointer :: invpr(:)   !< NEW  = INVPR(OLD)
    integer(ip), intent(inout), pointer :: array(:,:) !< Array
    integer(ip)                         :: ndofn
    integer(ip)                         :: ii,jj,nn
    integer(ip), allocatable            :: array_tmp(:,:)

    if( .not. associated(invpr) ) return
    nn    = size(invpr,KIND=ip)
    ndofn = size(array,1,KIND=ip)
    if( nn /= size(array,2,KIND=ip) ) call runend('MATHS_ARRAY_PERMUTATION: WROND DIMENSIONS')

    allocate(array_tmp(ndofn,nn))
    array_tmp = array
    do ii = 1,nn
       jj          = invpr(ii)
       array(:,jj) = array_tmp(:,ii)
    end do
    deallocate(array_tmp)

  end subroutine maths_array_permutation_IP_2

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    03/03/2016
  !> @brief   Time units
  !> @details Returns time scaling info according to a given value
  !>
  !-----------------------------------------------------------------------

  subroutine maths_time_unit(time_value,time_char,time_factor)
    real(rp),      intent(in)  :: time_value  !< Time in seconds
    character(*),  intent(out) :: time_char   !< Time unit character
    real(rp),      intent(out) :: time_factor !< Time scaling

    time_char = ' '

    if(      time_value < 1.0e-12_rp ) then  ! femto seconds
       time_factor = 1.0e15_rp
       time_char   = 'fs'
    else if( time_value < 1.0e-9_rp  ) then  ! pico seconds
       time_factor = 1.0e12_rp
       time_char   = 'ps'
    else if( time_value < 1.0e-6_rp  ) then  ! nano seconds
       time_factor = 1.0e9_rp
       time_char   = 'ns'
    else if( time_value < 1.0e-3_rp  ) then  ! micro seconds
       time_factor = 1.0e6_rp
       time_char   = 'mus'
    else if( time_value < 1.0_rp     ) then  ! milli seconds
       time_factor = 1.0e3_rp
       time_char   = 'ms '
    else if( time_value < 3600.0_rp  ) then  ! seconds
       time_factor = 1.0_rp
       time_char   = 's  '
    else                                     ! hours
       time_factor = 1.0_rp/3600.0_rp
       time_char   = 'h  '
    end if

  end subroutine maths_time_unit

  !----------------------------------------------------
  !>
  !> @author  Alfonso Santiago
  !> @date    15/06/2016
  !> @brief   economy size qrdecomposition with household reflexions for rectangular matrices
  !> @details compute a QR=A decomposition were Q is orthogonal and R upper triangular
  !>
  ! This subroutine decomposes the matrix A(m,n)
  ! in two matrices Q and R, where Q is orthogonal
  ! and R is upper triangular.
  !
  ! as Q is orthogonal Q^T*Q=I so  the inversion of
  ! this matrix is trivial.
  !
  ! Also:
  !       Q*R=A
  !
  ! In the original version of the algorithm
  ! Q(m,m) as it is orthogonal and R(m,n) for
  ! consistency in the size of the matrix.
  !
  ! But in the economy-size version, as R is upper
  ! triangular, you can build an R(n,n) matrix
  ! (because you know that the rest of the elements
  ! are zero and an Q(m,n), because the columns
  ! n->m are going to be multiplied by zero. In this
  ! algorithm A=Q*R holds.
  !                 n
  !            +--+--+--+
  !            |  |  |  |
  !            +--+--+--+
  !            |  |  |  |
  !  A(m,n)=   +--+--+--+ m
  !            |  |  |  |
  !            +--+--+--+
  !            |  |  |  |
  !            +--+--+--+
  !
  !
  !                n
  !            +--+--+--+
  !            |  |  |  |                   n
  !            +--+--+--+               +--+--+--+
  !            |  |  |  |               |  |  |  |
  !  Q(m,n)=   +--+--+--+ m             +--+--+--+
  !            |  |  |  |      R(m,n)=  |  |  |  | n
  !            +--+--+--+               +--+--+--+
  !            |  |  |  |               |  |  |  |
  !            +--+--+--+               +--+--+--+
  !
  !                       Q(m,n)*R(n,n) = A(m,n)
  !
  !
  !
  !
  !----------------------------------------------------
  
  subroutine maths_economy_qrhousehold(A,Q,R,mmax,nmax)

    implicit none

    real(rp), intent (in)               :: A(:,:)
    real(rp), intent (out)              :: Q(:,:)
    real(rp), intent (out)              :: R(:,:)
    integer(ip), intent(in)             :: mmax, nmax !! maximum values of the matrix to decompose
    !  integer(ip), intent(in)             :: m, n       !! actuall full values of the matrix
    integer(ip)                         :: i_col, i, j, k
    real(rp)                            :: norm_column, re_aux
    real(rp)                            :: v(mmax,nmax)
    real(rp)                            :: A_aux(mmax,nmax)

    !  m=size(A,1,KIND=ip)
    !  n=size(A,2,KIND=ip)
    !  if( ( m .lt. 0)             .or. &       ! A is pointing to memory
    !      ( n .lt. 0)             .or. &       !
    !      ( size(Q,1,KIND=ip) .lt. mmax)  .or. &       ! Q should be as big as mmax*mmax
    !      ( size(Q,2,KIND=ip) .lt. mmax)  .or. &       !
    !      ( size(R,1,KIND=ip) .lt. mmax)  .or. &       ! R shold be as big as mmax*nmax
    !      ( size(R,2,KIND=ip) .lt. nmax)  .or. &       !
    !      ( mmax .gt. size(A,1,KIND=ip) ) .or. &       ! mmax is always smaller than m
    !      ( nmax .gt. size(A,2,KIND=ip) ) .or. ) then  ! nmax is always smaller than n
    !      call runend('MATRIX DECOMPOSITION WRONG DIMENSIONS')
    !  endif

    Q=0.0_rp
    A_aux=A

    columns: do i_col=1, nmax
       ! a_i=A(i_col,i_col:mmax)
       v(:,i_col)=0.0_rp
       v(i_col:mmax,i_col) = A_aux(i_col:mmax,i_col)

       ! ||alpha|| = sqrt(sum(A(i)^2))
       !
       norm_column = 0.0_rp
       do i=i_col,mmax
          norm_column =norm_column + v(i,i_col) * v(i,i_col)
       enddo
       norm_column = sqrt(norm_column)

       ! u=a_i-alpha*e_1
       !
       !! IMPORTANT THING HERE. Some algoritmhs (i.e. matlab)
       !! uses a_i + alpha* e_i. This algorithm also converges
       !! but slightly different. The minus sign has been
       !! left here trivially.
       v(i_col,i_col) = v(i_col,i_col) - norm_column

       ! ||u|| = sqrt(sum(u(i)^2))
       !
       norm_column = 0.0_rp
       do i=i_col,mmax
          norm_column =norm_column + v(i,i_col) * v(i,i_col)
       enddo
       norm_column = sqrt(norm_column)

       ! v=u/||u||
       !
       do i=i_col,mmax
          v(i,i_col)=v(i,i_col)/norm_column
       enddo

       ! Q_i=I-2*v*v^T
       !
       ! obtained in function Q_mat(v,i_ini,i,j)

       ! A_i+1=Q_i * A
       !
       if(.not.(i_col.eq.nmax))then

          ! Tener la matrix A_i es indispensable no se puede tener solo
          ! el vector a_i ya que para el paso siguiente necesitaremos
          ! el a_ip1 y para el siguiente el a_ip2 y asi sucesivamente
          ! Por lo que se ha de tener toda la matrix A

          call vQ_times_Qaux(v(:,i_col),i_col,A_aux)

       endif

    enddo columns

    ! Q=Q_1*Q_2*Q_3*....Q_n
    !

    Q=0.0_rp
    Q=obtain_Q_mat(v(:,nmax),nmax,mmax)

    do i_col=nmax-1,1,-1
       call vQ_times_Qaux(v(:,i_col),i_col,Q)
    enddo

    ! A=Q*R  -> R=Q^T*A
    !
    !

    do i=1,nmax !swipe in rows
       do j=1,nmax
          re_aux=0.0_rp
          do k=1,mmax !swipe in columns
             re_aux = re_aux + Q(k,i) * A(k,j) !contract columns. Q has inversed indexes because it is the transpost
          enddo
          R(i,j)=re_aux
       enddo
    enddo

  contains
    !--------------------------------------
    ! Function that gives me the result of
    ! the matrix without having the matrix
    !--------------------------------------
    pure function Q_mat(v,i_ini,i,j) result(q_ij)
      implicit none

      real(rp), intent(in)        :: v(:)
      integer(ip), intent(in)     :: i_ini,i,j
      real(rp)                    :: q_ij

      if (i .lt. i_ini) then
         if(i .eq. j) then
            q_ij = 1.0_rp
         else
            q_ij = 0.0_rp
         endif
      else
         if(i .eq. j) then
            q_ij = 1.0_rp - 2.0_rp*v(i)*v(j)
         else
            q_ij = -2.0_rp*v(i)*v(j)
         endif
      endif

    end function Q_mat
    !--------------------------------------
    ! Function that gives me the result of
    ! the matrix without having the matrix
    !--------------------------------------
    function obtain_Q_mat(v,i_ini,max_col) result(Q)
      implicit none

      real(rp), intent(in)        :: v(:)
      integer(ip), intent(in)     :: i_ini, max_col
      integer(ip)                 :: i, j, m
      real(rp)                    :: Q(size(v,1,KIND=ip),size(v,1,KIND=ip))

      m=size(v,KIND=ip)

      if ( (max_col .gt. m)     .or. &
           (max_col .lt. i_ini) ) then
         ! call runend('mod_maths: problem with obtain_Q_mat. Wrong max column number')
         write(6,*) 'OBTAIN_Q_MAT ERROR: MAXIMUM COLUMN NUMBER PROBLEM '
      endif
      do i=1,m
         do j=1,max_col
            Q(i,j) = Q_mat(v,i_ini,i,j)
         enddo
      enddo

    end function obtain_Q_mat

    !--------------------------------------
    ! Function that gives me the result of
    ! the matrix without having the matrix
    !--------------------------------------
    subroutine vQ_times_Qaux(v,i_ini,Qaux)
      implicit none

      real(rp), intent(in)        :: v(:)
      integer(ip), intent(in)     :: i_ini
      real(rp), intent(inout)     :: Qaux(:,:)
      real(rp), allocatable       :: mat_aux(:,:)
      real(rp)                    :: re_aux
      integer(ip)                 :: i, j, mv,mQ,nQ

      mv=size(v,KIND=ip)
      mQ=size(Qaux,1,KIND=ip)
      nQ=size(Qaux,2,KIND=ip)

      if ( (mv .ne. mQ ) ) write(6,*) 'vQ_times_Qaux: wrong dimensions' !call runend(mod_maths: vQ_times_Qaux: wrong dimensions')

      allocate(mat_aux(mQ,nQ))

      do i=1,mv ! NWET swipe in rows
         do j=1,nQ
            re_aux=0.0_rp
            do k=1,mv !swipe in columns
               re_aux = re_aux + Q_mat(v, i_ini, i, k) * Qaux(k,j) !contract columns
            enddo
            mat_aux(i, j) = re_aux
         enddo
      enddo

      Qaux=mat_aux

      deallocate(mat_aux)

    end subroutine vQ_times_Qaux

    !--------------------------------------
    ! Function that multiplies two Q
    ! special matrices
    !--------------------------------------
    pure function vQ_times_vQ(v1, i_ini_1, v2, i_ini_2) result(Q)
      implicit none

      real(rp), intent(in)        :: v1(:), v2(:)
      integer(ip), intent(in)     :: i_ini_1, i_ini_2
      real(rp)                    :: Q(size(v1,1,KIND=ip),size(v1,1,KIND=ip))
      real(rp)                    :: rea_aux
      integer(ip)                 :: m,i,j,k

      m=size(v1,1,KIND=ip)
      Q=0.0_rp

      do i=1,m !swipe in rows
         do j=1,m
            rea_aux=0.0_rp
            do k=1,m !swipe in columns
               rea_aux = rea_aux + Q_mat(v1, i_ini_1, i, k) * Q_mat(v2, i_ini_2, k, j)!contract columns
            enddo
            Q(i,j)=rea_aux
         enddo
      enddo

    end function vQ_times_vQ
  end subroutine maths_economy_qrhousehold

  !----------------------------------------------------
  !>
  !> @author  Alfonso Santiago
  !> @date    08/06/2016
  !> @brief   matrix multiplication
  !> @details Multiplies two rectangular matrices
  !>
  !----------------------------------------------------
  
  subroutine maths_matrix_multiplication(A,B,C,nn)
    !
    ! A(m,n)*B(n,p)=C(m,p)
    !
    real(rp),              intent(in)  :: A(:,:)
    real(rp),              intent(in)  :: B(:,:)
    real(rp),              intent(out) :: C(:,:)
    integer(ip), optional, intent(in)  :: nn
    integer(ip)                        :: i,j,k
    integer(ip)                        :: m,n,r

    C = 0.0_rp

    if( present(nn) ) then
       do i = 1,nn 
          do j = 1,nn
             C(i,j) = 0.0_rp
             do k = 1,nn 
                C(i,j) = C(i,j) + A(i,k)*B(k,j) 
             end do
          end do
       end do
    else
       m = size(A,1,KIND=ip) ! Rows
       n = size(B,2,KIND=ip) ! Columns
       r = size(A,2,KIND=ip) ! reduce
       if( r /= size(B,1,KIND=ip) ) then
          CALL RUNEND('MATRIX_MULTIPLIUCATION: WRONG DIMENSIONS IN MATRIX')
       end if
       do i = 1,m !swipe in rows
          do j = 1,n
             C(i,j) = 0.0_rp
             do k = 1,r !swipe in columns
                C(i,j) = C(i,j) + A(i,k)*B(k,j) !contract columns
             end do
          end do
       end do
    end if

  end subroutine maths_matrix_multiplication

  !----------------------------------------------------
  !>
  !> @author  Alfonso Santiago
  !> @date    08/06/2016
  !> @brief   matrix multiplication
  !> @details Multiplies two rectangular matrices
  !>
  !
  !
  !----------------------------------------------------

  subroutine maths_matrix_vector_multiplication_0(A,B,C,maxn)
    !
    ! A(m,n)*B(n,1)=C(m,1)
    !
    real(rp),    intent(in)  :: A(:,:)
    real(rp),    intent(in)  :: B(:)
    real(rp),    intent(out) :: C(:)
    integer(ip), intent(in)  :: maxn
    integer(ip)              :: m,n
    integer(ip)              :: i,j

    m = size(A,1,KIND=ip) ! Rows
    n = size(A,2,KIND=ip) ! Columns

    if(size(C,1,KIND=ip) == 0_ip) return

    if( (n /= size(B,1,KIND=ip)) .or. (maxn > n) ) then
       CALL RUNEND('MATRIX_VECTOR_MULTIPLIUCATION: WRONG DIMENSIONS IN MATRIX')
    endif

    C = 0.0_rp
    do i = 1,m !swipe in rows
       C(i) = 0.0_rp
        do j = 1,maxn
          C(i) = C(i) + A(i,j) * B(j) !contract columns
       end do
    end do

  end subroutine maths_matrix_vector_multiplication_0

  pure subroutine maths_matrix_vector_multiplication_1(m,n,A,B,C)
    !
    ! A(m,n)*B(n,1)=C(m,1)
    !
    integer(ip), intent(in)  :: m
    integer(ip), intent(in)  :: n
    real(rp),    intent(in)  :: A(m,n)
    real(rp),    intent(in)  :: B(n)
    real(rp),    intent(out) :: C(m)
    integer(ip)              :: i,j

    C = 0.0_rp
    do i = 1,m 
       C(i) = 0.0_rp
        do j = 1,n
          C(i) = C(i) + A(i,j) * B(j)
       end do
    end do

  end subroutine maths_matrix_vector_multiplication_1
  
  pure subroutine maths_matrix_transpose_vector_multiplication(m,n,A,B,C)
    !
    ! A(m,n)*B(n,1)=C(m,1)
    !
    integer(ip), intent(in)  :: m
    integer(ip), intent(in)  :: n
    real(rp),    intent(in)  :: A(m,n)
    real(rp),    intent(in)  :: B(n)
    real(rp),    intent(out) :: C(m)
    integer(ip)              :: i,j
    
    C = 0.0_rp
    do i = 1,m 
       C(i) = 0.0_rp
        do j = 1,n
          C(i) = C(i) + A(j,i) * B(j)
       end do
    end do

  end subroutine maths_matrix_transpose_vector_multiplication

  !----------------------------------------------------
  !>
  !> @author  Alfonso Santiago
  !> @date    08/06/2016
  !> @brief   Back substitution
  !> @details backsubtitutes the sistem Ax=B were A is upper triangular
  !>
  !
  !
  !----------------------------------------------------
  
  subroutine maths_backsu(A,x,b, nmax)

    real(rp), intent(in)        :: A(:,:)
    real(rp), intent(in)        :: b(:)
    integer(ip), intent(in)     :: nmax !bounds of the operation
    real(rp), intent(out)       :: x(:)
    integer(ip)                 :: i,j
    real(rp)                    :: aux

    x=0.0_rp

    if (size(b,1,KIND=ip).eq.0_ip) return

    do i=nmax,1,-1

       aux=0.0_rp

       do j=nmax,i+1,-1
          aux=aux+A(i,j)*x(j)
       enddo

       if(A(i,i).eq. 0.0_rp) call runend('mod_maths backsubstitution: diagonal equal to zero in backsubstitution')

       x(i)=(b(i)-aux)/A(i,i)

    enddo

  end subroutine maths_backsu

  !-----------------------------------------------------------
  !>
  !> @author  J.C. Cajas
  !> @date    22/06/2016
  !> @brief   Outer product
  !> @details Outer product between two vectors A(m) and B(n)
  !> @        result is C(mxn)
  !>
  !>       _    _
  !>      |  a1  |
  !> A =  |  a2  |
  !>      |   .  |
  !>      |   .  |
  !>      |   .  |
  !>      |_ am _|
  !>
  !> B = [b1,b2,...,bn]
  !>
  !>           --                      --
  !>          | a1 b1   a1 b2  ... a1 bn |
  !>          | a2 b1   a2 b2      a2 bn |
  !>          |  .                       |
  !> A xo B = |  .                       |
  !>          |  .                       |
  !>          | am b1     ...      am bn |
  !>           --                      --
  !------------------------------------------------------------
  
  subroutine maths_outer_product(A,B,C)
    real(rp), intent(in) ,   pointer :: A(:),B(:)   !< Input vectors
    real(rp), intent(inout), pointer :: C(:,:)      !< Output matrix
    integer(ip)                      :: mm,nn       !< Dimensions of the vectors
    integer(ip)                      :: ii,jj

    mm = size(A,KIND=ip)
    nn = size(B,KIND=ip)

    if( mm /= size(C,1_ip,KIND=ip) .or. nn /= size(C,2_ip,KIND=ip) )call runend('mod_maths, outer product: Wrong dimensions')

    do jj = 1_ip, nn
       do ii = 1_ip, mm
          C(ii,jj) = A(ii) * B(jj)
       end do
    end do

  end subroutine maths_outer_product

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   L2 norm of an array
  !> @details L2 norm of an array
  !
  !------------------------------------------------------------------------

  pure function maths_norm2(ndime,array)
    integer(ip), intent(in) :: ndime
    real(rp),    intent(in) :: array(ndime)
    real(rp)                :: maths_norm2

    maths_norm2 = sqrt(dot_product(array,array))

  end function maths_norm2

  !-----------------------------------------------------------------------
  !
  !> @date    27/10/2014
  !> @author  Ricard Borrell
  !> @brief   Partition of bin boxes by means of space filling curves
  !> @details npart can be a real number
  !
  !-----------------------------------------------------------------------
  
  subroutine maths_sfc_part(ndime,boxes,weigh,npart,parts,fronw,irank,nrank,ityp)

    integer(ip),          intent(in)     :: ndime         !< Dimension of problem
    integer(ip),          intent(in)     :: boxes(ndime)  !< # boxes in each dimension
    real(rp),    pointer, intent(in)     :: weigh(:)      !< boxes weights
    real(rp),             intent(in)     :: npart         !< # partitions
    integer(ip), pointer, intent(inout)  :: parts(:)      !< Partition to which each box is assigned
    real(rp),    optional,intent(in)     :: fronw         !< weigh of first partition
    integer(ip), optional,intent(in)     :: irank         !> first output partitions rank
    integer(ip), optional,intent(in)     :: nrank         !< # ranks used in the partition
    integer(2), optional,intent(inout)   :: ityp          ! orientation type

    integer(ip)                       :: irank_,nrank_
    integer(ip)                       :: aux,icont,ipart
    integer(ip)                       :: nbox,ibox,frank
    integer(ip)                       :: hcoor(ndime),coo1d
    integer(2)                        :: auxtyp,auxtyp2
    real(rp)                          :: fronw_
    real(rp)                          :: sump,iweig,wsum
    real(rp)                          :: wprev,waux,wobj

    !
    ! Check arguments
    !

    if(ndime/= 2_ip .and. ndime/= 3_ip)&
         call runend("maths_sfc_part: wrong ndime value")
    if(boxes(1) <= 0_ip) call runend("maths_sfc_part: boxes(1)<=0")
    nbox = boxes(1)
    do icont =2_ip,ndime
       if(boxes(icont) /= boxes(icont-1_ip)) call runend("maths_sfc_part: different boxes sizes")
       nbox = nbox * boxes(1)
    enddo
    if(iand(nbox,nbox-1_ip)/=0)&
         call runend("math_sfc_part: nbox is not power of two")
    if( size(weigh,KIND=ip) /= nbox )&
         call runend("math_sfc_part: Incongruent number of boxes and weigh dim")

    !
    ! Check or optional inputs
    !

    if(present(fronw))then
       if(fronw <= 0_rp .or. fronw> 1_rp) then
          call runend("maths_sfc_part: fronw <= 0 or fronw > 1.0!")
       endif
       fronw_=fronw
    else
       fronw_ = 1_rp
    endif

    if(present(irank)) then
       if(irank<=0_ip) call runend("maths_sfc_part: irank<=1")
       irank_=irank
    else
       irank_ = 1_ip
    endif

    aux = int(npart,ip)-int(fronw_,KIND=ip)
    if( npart-int(fronw_,KIND=ip)-aux > 0_ip) aux = aux + 1_ip
    aux = aux + 1
    if(present(nrank)) then
       if(nrank > aux .or. nrank < aux-1) then
          call runend("maths_sfc_part: nrank out of rank")
       endif
       nrank_=nrank
    else
       nrank_ = aux
    endif

    if(present(ityp)) then
       auxtyp = ityp
    else
       auxtyp = 0_2
    endif

    !
    ! Allocate parts
    !

    nullify(parts)
    call memory_alloca(memor,'parts',"maths_sfc_part",parts,nbox)

    !
    ! calculate sum of weights and initial obj. weight
    !
    frank = nrank_ + irank_ - 1_ip
    wsum = 0_rp
    do ibox=1,nbox
       wsum=wsum+weigh(ibox)
    enddo
    wobj = (wsum/npart)*fronw_
    sump=fronw_
    !
    ! assign a partition to each box
    !
    iweig = 0_rp
    wprev = 0_rp
    ipart = irank_

    do ibox=1_ip,nbox

       if(ndime==2_ip) then
          auxtyp2=auxtyp
          call maths_sfc_d2xy_tab(boxes(1),ibox,hcoor(1),hcoor(2),auxtyp2)
          coo1d = (hcoor(2)-1)*boxes(1)+hcoor(1)
       else
          auxtyp2=auxtyp
          call maths_sfc_d2xyz_tab(boxes(1),ibox,hcoor(1),hcoor(2),hcoor(3),auxtyp2)
          coo1d = (hcoor(3)-1)*boxes(2)*boxes(1)+(hcoor(2)-1)*boxes(1)+hcoor(1)
       endif
       waux  = weigh(coo1d)

       if( abs(iweig+waux-wobj) <= abs(iweig-wobj) .or. frank == ipart ) then
          parts(coo1d)=ipart
          iweig = iweig+waux
          wprev = wprev+waux
       else
          wobj  = (wsum-wprev)/real(npart-sump,rp)
          ipart = ipart + 1_ip
          parts(coo1d) = ipart
          iweig = waux
          wprev = wprev + waux
          sump = sump + 1_ip
       endif

    enddo

  end subroutine maths_sfc_part
  
  !-----------------------------------------------------------------------
  !
  !> @date    27/10/2014
  !> @author  Ricard Borrell
  !> @brief   Partition of bin boxes by means of space filling curves
  !> @details npart can be a real number
  !
  !-----------------------------------------------------------------------

  subroutine maths_sfc_part_2_ip(ndime,boxes,weigh,parts,ityp,irank,nrank,wrank,PART_NAME)

    integer(ip),                   intent(in)     :: ndime         !< Dimension of problem
    integer(ip),                   intent(in)     :: boxes(ndime)  !< # boxes in each dimension
    integer(ip), pointer,          intent(in)     :: weigh(:)      !< boxes weights
    integer(ip), pointer,          intent(inout)  :: parts(:)      !< Partition to which each box is assigned
    integer(2),           optional,intent(inout)  :: ityp          ! orientation type
    integer(ip),          optional,intent(in)     :: irank         !> first output partitions rank
    integer(ip),          optional,intent(in)     :: nrank         !< # ranks used in the partition
    real(rp),    pointer, optional,intent(in)     :: wrank(:)      !< relative weight per rank
    character(*),         optional,intent(in)     :: PART_NAME
    integer(ip)                                   :: irank_,nrank_
    integer(ip)                                   :: icont,ipart,ilpart
    integer(ip)                                   :: nbox,ibox,frank
    integer(ip)                                   :: hcoor(ndime),coo1d
    integer(2)                                    :: auxtyp,auxtyp2
    real(rp)                                      :: sump,iweig,wsum
    real(rp)                                      :: wprev,waux,wobj
    real(rp),     pointer                         :: wrank_(:)
    real(rp)                                      :: npart        

    character(100), PARAMETER :: vacal = "maths_sfc_part_2"
    !
    ! TODO: I supose by now that all optional arguments are present...
    !

    !
    ! Check arguments
    !

    if(ndime/= 2_ip .and. ndime/= 3_ip)&
         call runend("maths_sfc_part: wrong ndime value")
    if(boxes(1) <= 0_ip) call runend("maths_sfc_part: boxes(1)<=0")
    nbox = boxes(1)
    do icont =2_ip,ndime
       if(boxes(icont) /= boxes(icont-1_ip)) call runend("maths_sfc_part: different boxes sizes")
       nbox = nbox * boxes(1)
    enddo
    if(iand(nbox,nbox-1_ip)/=0)&
         call runend("math_sfc_part: nbox is not power of two")
    if( size(weigh,KIND=ip) /= nbox )&
         call runend("math_sfc_part: Incongruent number of boxes and weigh dim")

    ! rick: check combinations of elements wich have to be present...

    !
    ! Check optional inputs
    !

    if(present(ityp)) then
       auxtyp = ityp
    else
       auxtyp = 0_2
    endif

    if(present(irank)) then
       if(irank<=0_ip) call runend("maths_sfc_part: irank<=1")
       irank_=irank
    else
       irank_ = 1_ip
    endif

    if(present(nrank)) then
       nrank_=nrank
       nullify(wrank_)
       allocate(wrank_(nrank_))
       wrank_(1:nrank)=wrank(1:nrank)
    endif
    !
    ! Allocate parts
    !
    if( present(PART_NAME) ) then
       call memory_alloca(memor,trim(PART_NAME),"maths_sfc_part",parts,nbox)
    else
       call memory_alloca(memor,'parts',"maths_sfc_part",parts,nbox)
    end if
    !
    ! calculate sum of weights and initial obj. weight
    !
    npart = sum(wrank_)
    frank = nrank_ + irank_ - 1_ip  !! final rank?
    wsum = 0.0_rp
    do ibox=1,nbox
       wsum=wsum+real(weigh(ibox),rp)
    enddo
    
    wobj = (wsum/npart)*wrank_(1)
    sump=wrank_(1)
    !
    ! assign a partition to each box
    !
    iweig = 0.0_rp
    wprev = 0.0_rp
    ipart = irank_
    ilpart = 1

    do ibox=1_ip,nbox

       if(ndime==2_ip) then
          auxtyp2=auxtyp
          call maths_sfc_d2xy_tab(boxes(1),ibox,hcoor(1),hcoor(2),auxtyp2)
          coo1d = (hcoor(2)-1)*boxes(1)+hcoor(1)
       else
          auxtyp2=auxtyp
          call maths_sfc_d2xyz_tab(boxes(1),ibox,hcoor(1),hcoor(2),hcoor(3),auxtyp2)
          coo1d = (hcoor(3)-1)*boxes(2)*boxes(1)+(hcoor(2)-1)*boxes(1)+hcoor(1)
       endif
       waux  = real(weigh(coo1d),rp)

       if( abs(iweig+waux-wobj) <= abs(iweig-wobj) .or. frank == ipart ) then
          parts(coo1d)=ipart
          iweig = iweig+waux
          wprev = wprev+waux
       else
          ilpart = ilpart + 1_ip
          wobj  = ((wsum-wprev)/real(npart-sump,rp))*wrank_(ilpart)
          ipart = ipart + 1_ip
          parts(coo1d) = ipart
          iweig = waux
          wprev = wprev + waux
          sump = sump + wrank_(ilpart)
       endif

    enddo

    !rick: see what happends when no present
    deallocate(wrank_)

  end subroutine maths_sfc_part_2_ip
  
  subroutine maths_sfc_part_2_rp(ndime,boxes,weigh,parts,ityp,irank,nrank,wrank,PART_NAME)

    integer(ip),          intent(in)     :: ndime         !< Dimension of problem
    integer(ip),          intent(in)     :: boxes(ndime)  !< # boxes in each dimension
    real(rp),    pointer, intent(in)     :: weigh(:)      !< boxes weights
    integer(ip), pointer, intent(inout)  :: parts(:)      !< Partition to which each box is assigned

    integer(2),  optional,intent(inout)  :: ityp          ! orientation type
    integer(ip), optional,intent(in)     :: irank         !> first output partitions rank
    integer(ip), optional,intent(in)     :: nrank         !< # ranks used in the partition
    real(rp),    pointer,optional,intent(in)  :: wrank(:) !< relative weight per rank
    character(*),         optional,intent(in) :: PART_NAME
    integer(ip)                       :: irank_,nrank_
    integer(ip)                       :: icont,ipart,ilpart
    integer(ip)                       :: nbox,ibox,frank
    integer(ip)                       :: hcoor(ndime),coo1d
    integer(2)                        :: auxtyp,auxtyp2
    real(rp)                          :: sump,iweig,wsum
    real(rp)                          :: wprev,waux,wobj
    real(rp), pointer                 :: wrank_(:)
    real(rp)                          :: npart         !< # partitions

    character(100), PARAMETER :: vacal = "maths_sfc_part_2"
    !
    ! TODO: I supose by now that all optional arguments are present...
    !

    !
    ! Check arguments
    !

    if(ndime/= 2_ip .and. ndime/= 3_ip)&
         call runend("maths_sfc_part: wrong ndime value")
    if(boxes(1) <= 0_ip) call runend("maths_sfc_part: boxes(1)<=0")
    nbox = boxes(1)
    do icont =2_ip,ndime
       if(boxes(icont) /= boxes(icont-1_ip)) call runend("maths_sfc_part: different boxes sizes")
       nbox = nbox * boxes(1)
    enddo
    if(iand(nbox,nbox-1_ip)/=0)&
         call runend("math_sfc_part: nbox is not power of two")
    if( size(weigh,KIND=ip) /= nbox )&
         call runend("math_sfc_part: Incongruent number of boxes and weigh dim")

    ! rick: check combinations of elements wich have to be present...

    !
    ! Check optional inputs
    !

    if(present(ityp)) then
       auxtyp = ityp
    else
       auxtyp = 0_2
    endif

    if(present(irank)) then
       if(irank<=0_ip) call runend("maths_sfc_part: irank<=1")
       irank_=irank
    else
       irank_ = 1_ip
    endif

    if(present(nrank)) then
       nrank_=nrank
       nullify(wrank_)
       call memory_alloca(memor,'wrank_',vacal,wrank_,nrank_)
       wrank_(1:nrank)=wrank(1:nrank)
    endif

    !
    ! Allocate parts
    !

    !nullify(parts)
    if( present(PART_NAME) ) then
       call memory_alloca(memor,trim(PART_NAME),"maths_sfc_part",parts,nbox)!,REALLOCATE=.true.)
    else
       call memory_alloca(memor,'parts',"maths_sfc_part",parts,nbox)!,REALLOCATE=.true.)
    end if
    !
    ! calculate sum of weights and initial obj. weight
    !
    npart = sum(wrank_)
    frank = nrank_ + irank_ - 1_ip  !! final rank?
    wsum = 0_rp
    do ibox=1,nbox
       wsum=wsum+weigh(ibox)
    enddo

    wobj = (wsum/npart)*wrank_(1)
    sump=wrank_(1)

    ! assign a partition to each box
    !
    iweig = 0_rp
    wprev = 0_rp
    ipart = irank_
    ilpart = 1

    do ibox=1_ip,nbox

       if(ndime==2_ip) then
          auxtyp2=auxtyp
          call maths_sfc_d2xy_tab(boxes(1),ibox,hcoor(1),hcoor(2),auxtyp2)
          coo1d = (hcoor(2)-1)*boxes(1)+hcoor(1)
       else
          auxtyp2=auxtyp
          call maths_sfc_d2xyz_tab(boxes(1),ibox,hcoor(1),hcoor(2),hcoor(3),auxtyp2)
          coo1d = (hcoor(3)-1)*boxes(2)*boxes(1)+(hcoor(2)-1)*boxes(1)+hcoor(1)
       endif
       waux  = weigh(coo1d)

       if( abs(iweig+waux-wobj) <= abs(iweig-wobj) .or. frank == ipart ) then
          parts(coo1d)=ipart
          iweig = iweig+waux
          wprev = wprev+waux
       else
          ilpart = ilpart + 1_ip
          wobj  = ((wsum-wprev)/real(npart-sump,rp))*wrank_(ilpart)
          ipart = ipart + 1_ip
          parts(coo1d) = ipart
          iweig = waux
          wprev = wprev + waux
          sump = sump + wrank_(ilpart)
       endif

    enddo

    !rick: see what happends when no present
    call memory_deallo(memor,'wrank_',vacal,wrank_)

  end subroutine maths_sfc_part_2_rp
  
  !-----------------------------------------------------------------------
  !
  !> @date    27/10/2014
  !> @author  Ricard Borrell
  !> @brief   Parallel Partition of bin boxes by means of space filling curves
  !> @details ...
  !
  !-----------------------------------------------------------------------
  
  subroutine maths_sfc_par_part(ndime,loc_boxes,loc_weigh,npar,proc_boxes,mypro,proc_wdist,loc_part)

    integer(ip),          intent(in)     :: ndime             !< Dimension of problem
    integer(ip),          intent(in)     :: loc_boxes(ndime)  !< # boxes in each dimension (local)
    real(rp),    pointer, intent(in)     :: loc_weigh(:)      !< boxes weights (local)
    integer(ip),          intent(in)     :: npar              !< # partitions (global)
    integer(ip),          intent(in)     :: proc_boxes(ndime) !< # processes performing partition
    integer(ip),          intent(in)     :: mypro             !< # my process rank
    real(rp), pointer,    intent(in)     :: proc_wdist(:)     !< weig accumulated by each process
    integer(ip), pointer, intent(inout)  :: loc_part(:)       !< Partition to which each box is assigned (local)

    real(rp),    pointer                 :: loc_npar(:)
    real(rp),    pointer                 :: fronw(:)
    integer(ip), pointer                 :: irank(:)
    integer(ip), pointer                 :: reord(:)
    integer(ip), pointer                 :: invor(:)
    integer(2),  pointer                 :: types(:)
    integer(ip)                          :: coo1d
    integer(2)                           :: auxty
    real(rp)                             :: wsum,npsum
    integer(ip)                          :: iproc,x,y,z
    integer(ip)                          :: nproc,idime
    integer(ip)                          :: invmy

    character(100), PARAMETER :: vacal = "maths_sfc_par_part"

    !
    ! Check proc_boxes input
    !
    do idime=2_ip,ndime
       if(proc_boxes(idime)/=proc_boxes(idime-1_ip))&
            call runend("maths_sfc_par_part: wrong proc_boxes")
    enddo
    nproc = proc_boxes(1_ip)**(ndime)
    if(iand(nproc,nproc-1)/=0)&
         call runend("math_sfc_part: boxes(1) is not a power of 2")
    !
    ! Eval the new order of proc. boxes and its orientation type
    !
    nullify(loc_npar,fronw,irank,reord,invor,types)
    call memory_alloca(memor,'loc_npar',vacal,loc_npar,nproc)
    call memory_alloca(memor,'fronw',vacal,fronw,nproc)
    call memory_alloca(memor,'irank',vacal,irank,nproc)
    call memory_alloca(memor,'reord',vacal,reord,nproc)
    call memory_alloca(memor,'invor',vacal,invor,nproc)
    allocate(types(nproc))
    wsum  = 0_rp
    do iproc=1,nproc
       wsum = wsum + proc_wdist(iproc)
       auxty = 0_ip
       if(ndime==2_ip) then
          call maths_sfc_d2xy_tab(proc_boxes(1),iproc,x,y,auxty)
          coo1d=(y-1)*proc_boxes(1)+x
       else
          call maths_sfc_d2xyz_tab(proc_boxes(1),iproc,x,y,z,auxty)
          coo1d=(z-1)*proc_boxes(1)*proc_boxes(1)+(y-1)*proc_boxes(1)+x
       endif
       reord(iproc) = coo1d
       types(iproc) = auxty
       invor(coo1d) = iproc
    enddo

    if(wsum == 0_rp) return
    !
    ! Perform partition
    !
    npsum = 0_rp
    do iproc=1,nproc
       loc_npar(iproc) = (proc_wdist(reord(iproc))/wsum) * real(npar,rp)
       irank(iproc) = int(npsum) + 1_ip
       fronw(iproc) = real(irank(iproc),rp) - npsum
       npsum = npsum + loc_npar(iproc)
    enddo

    invmy = invor(mypro)
    if(loc_npar(invmy) /= 0_rp) then
       call maths_sfc_part(ndime,loc_boxes,loc_weigh,loc_npar(invmy),loc_part,&
            fronw=fronw(invmy),irank=irank(invmy),ityp=types(invmy))
    endif
    call memory_deallo(memor,'loc_npar',vacal,loc_npar)
    call memory_deallo(memor,'fronw',vacal,fronw)
    call memory_deallo(memor,'irank',vacal,irank)
    call memory_deallo(memor,'reord',vacal,reord)
    call memory_deallo(memor,'invor',vacal,invor)
    deallocate(types)

  end subroutine maths_sfc_par_part
  
  !-----------------------------------------------------------------------
  !
  !> @date    27/10/2014
  !> @author  Ricard Borrell
  !> @brief   Parallel Partition of bin boxes by means of space filling curves
  !> @details ...
  !
  !-----------------------------------------------------------------------
  
  subroutine maths_sfc_par_part_2_ip(ndime,loc_boxes,loc_weigh,npar,proc_boxes,&
       mypro,proc_wdist,loc_part,corr_dist,totalw,PART_NAME)

    integer(ip),                    intent(in)    :: ndime                 !< Dimension of problem
    integer(ip),                    intent(in)    :: loc_boxes(ndime)      !< # boxes in each dimension (local)
    integer(ip), pointer,           intent(in)    :: loc_weigh(:)          !< boxes weights (local)
    integer(ip),                    intent(in)    :: npar                  !< # partitions (global)
    integer(ip),                    intent(in)    :: proc_boxes(ndime)     !< # processes performing partition
    integer(ip),                    intent(in)    :: mypro                 !< # my process rank
    integer(ip), pointer,           intent(in)    :: proc_wdist(:)         !< weig accumulated by each process
    integer(ip), pointer,           intent(inout) :: loc_part(:)           !< Partition to which each box is assigned (local)
    real(rp),    pointer, optional, intent(inout) :: corr_dist(:) !< correction factors for the distribution
    real(rp),             optional, intent(in)    :: totalw
    character(*),         optional, intent(in)    :: PART_NAME
    real(rp),    pointer                          :: loc_npar(:)
    real(rp),    pointer                          :: dist_corr(:)
    integer(ip), pointer                          :: irank(:)
    integer(ip), pointer                          :: reord(:)
    integer(ip), pointer                          :: invor(:)
    integer(2),  pointer                          :: types(:)
    integer(ip)                                   :: coo1d
    integer(2)                                    :: auxty
    real(rp)                                      :: wsum
    integer(ip)                                   :: iproc,x,y,z
    integer(ip)                                   :: nproc,idime
    integer(ip)                                   :: nrank
    real(rp),    pointer                          :: wrank(:)
    logical                                       :: auxlg
    integer(ip)                                   :: iaux

    !
    ! Check proc_boxes input
    !
    do idime=2_ip,ndime
       if(proc_boxes(idime)/=proc_boxes(idime-1_ip))&
            call runend("maths_sfc_par_part: wrong proc_boxes")
    enddo
    nproc = proc_boxes(1_ip)**(ndime)
    if(iand(nproc,nproc-1)/=0)&
         call runend("math_sfc_part: boxes(1) is not a power of 2")
    !
    ! Eval the new order of proc. boxes and its orientation type
    !
    nullify(loc_npar,irank,reord,invor,types,dist_corr)
    allocate(loc_npar(nproc))
    allocate(dist_corr(npar))
    allocate(irank(nproc))
    allocate(reord(nproc))
    allocate(invor(nproc))
    allocate(types(nproc))
    if(.NOT. present(corr_dist)) then
       dist_corr(:) = 1_rp
    else
       iaux = min(size(dist_corr,KIND=ip),size(corr_dist,KIND=ip))
       dist_corr(1:iaux) = corr_dist(1:iaux)
    endif
    !
    ! Ordering of partition processess acording to coarse SFC
    !
    wsum  = 0.0_rp
    if( ndime == 2 ) then
       do iproc=1,nproc
          wsum         = wsum + real(proc_wdist(iproc),rp)
          auxty        = 0_ip
          call maths_sfc_d2xy_tab(proc_boxes(1),iproc,x,y,auxty)
          coo1d        = (y-1)*proc_boxes(1)+x
          reord(iproc) = coo1d
          types(iproc) = auxty
          invor(coo1d) = iproc
       enddo
    else
       do iproc=1,nproc
          wsum         = wsum + real(proc_wdist(iproc),rp)
          auxty        = 0_ip
          call maths_sfc_d2xyz_tab(proc_boxes(1),iproc,x,y,z,auxty)
          coo1d        = (z-1)*proc_boxes(1)*proc_boxes(1)+(y-1)*proc_boxes(1)+x
          reord(iproc) = coo1d
          types(iproc) = auxty
          invor(coo1d) = iproc
       enddo
    end if
    if(present(totalw)) wsum=totalw;
    if(wsum == 0.0_rp) return
    !
    ! Perform partition
    !
    nullify(wrank)
    allocate(wrank(npar))

    iaux  = 1_ip
    wrank = 0.0_rp

    do iproc=1,nproc

       nrank = 0
       loc_npar(iproc) = ( real(proc_wdist(iproc),rp)/wsum) * real(npar,rp)
       if(loc_npar(iproc) > 0.0_rp) then

          irank(iproc) = iaux
          auxlg = .true.

          do while(auxlg .and. nrank < npar .and. iaux <= npar)
             loc_npar(iproc) = loc_npar(iproc) - dist_corr(iaux)
             if(loc_npar(iproc) < 0) then
                auxlg = .false.
                wrank(nrank+1)= loc_npar(iproc) + dist_corr(iaux)
                dist_corr(iaux) = -loc_npar(iproc)
             else
                wrank(nrank+1)=dist_corr(iaux)
                iaux = iaux + 1
             endif
             nrank = nrank + 1    ! si queda 0 el resido petara - repassar aquest cas
             if(loc_npar(iproc) == 0) auxlg = .false.
          end do
          if(iproc == nproc) then
             nrank = npar-irank(iproc)+1
          endif

          !loc_npar(iproc) = (proc_wdist(reord(iproc))/wsum) * real(npar,rp)
          loc_npar(iproc) = (real(proc_wdist(iproc),rp)/wsum) * real(npar,rp)

          !if(invor(mypro) == iproc) then
          if(mypro == iproc) then
             call maths_sfc_part_2(ndime,loc_boxes,loc_weigh,loc_part,types(iproc),&
                  irank(iproc),nrank,wrank,PART_NAME)
             exit
          endif
       endif
    enddo
    !
    ! deallocate pointers
    !
    deallocate(wrank,loc_npar,dist_corr,irank,reord,invor,types)

  end subroutine maths_sfc_par_part_2_ip
  
  subroutine maths_sfc_par_part_2_rp(ndime,loc_boxes,loc_weigh,npar,proc_boxes,&
       mypro,proc_wdist,loc_part,corr_dist,totalw,PART_NAME)

    integer(ip),          intent(in)              :: ndime                 !< Dimension of problem
    integer(ip),          intent(in)              :: loc_boxes(ndime)      !< # boxes in each dimension (local)
    real(rp),    pointer, intent(in)              :: loc_weigh(:)          !< boxes weights (local)
    integer(ip),          intent(in)              :: npar                  !< # partitions (global)
    integer(ip),          intent(in)              :: proc_boxes(ndime)     !< # processes performing partition
    integer(ip),          intent(in)              :: mypro                 !< # my process rank
    real(rp),    pointer, intent(in)              :: proc_wdist(:)         !< weig accumulated by each process
    integer(ip), pointer, intent(inout)           :: loc_part(:)           !< Partition to which each box is assigned (local)
    real(rp),    pointer, intent(inout), optional :: corr_dist(:) !< correction factors for the distribution
    real(rp),             intent(in),    optional :: totalw
    character(*),         intent(in),    optional :: PART_NAME

    real(rp),    pointer                          :: loc_npar(:)
    real(rp),    pointer                          :: dist_corr(:)
    integer(ip), pointer                          :: irank(:)
    integer(ip), pointer                          :: reord(:)
    integer(ip), pointer                          :: invor(:)
    integer(2),  pointer                          :: types(:)
    integer(ip)                                   :: coo1d
    integer(2)                                    :: auxty
    real(rp)                                      :: wsum
    integer(ip)                                   :: iproc,x,y,z
    integer(ip)                                   :: nproc,idime
    integer(ip)                                   :: nrank
    real(rp),    pointer                          :: wrank(:)
    logical                                       :: auxlg
    integer(ip)                                   :: iaux

    character(100), PARAMETER :: vacal = "maths_sfc_par_part_2"

    !
    ! Check proc_boxes input
    !
    do idime=2_ip,ndime
       if(proc_boxes(idime)/=proc_boxes(idime-1_ip))&
            call runend("maths_sfc_par_part: wrong proc_boxes")
    enddo
    nproc = proc_boxes(1_ip)**(ndime)
    if(iand(nproc,nproc-1)/=0)&
         call runend("math_sfc_part: boxes(1) is not a power of 2")
    !
    ! Eval the new order of proc. boxes and its orientation type
    !
    nullify(loc_npar,irank,reord,invor,types,dist_corr)
    call memory_alloca(memor,'loc_npar' ,vacal,loc_npar,nproc)
    call memory_alloca(memor,'dist_corr',vacal,dist_corr,npar)
    call memory_alloca(memor,'irank'    ,vacal,irank,nproc)
    call memory_alloca(memor,'reord'    ,vacal,reord,nproc)
    call memory_alloca(memor,'invor'    ,vacal,invor,nproc)
    allocate(types(nproc))
    if(.NOT. present(corr_dist)) then
       dist_corr(:) = 1_rp
    else
       iaux = min(size(dist_corr,KIND=ip),size(corr_dist,KIND=ip))
       dist_corr(1:iaux) = corr_dist(1:iaux)
    endif
    !
    ! Ordering of partition processess acording to coarse SFC
    !
    wsum  = 0.0_rp
    if( ndime == 2 ) then
       do iproc=1,nproc
          wsum         = wsum + proc_wdist(iproc)
          auxty        = 0_ip
          call maths_sfc_d2xy_tab(proc_boxes(1),iproc,x,y,auxty)
          coo1d        = (y-1)*proc_boxes(1)+x
          reord(iproc) = coo1d
          types(iproc) = auxty
          invor(coo1d) = iproc
       enddo
    else
       do iproc=1,nproc
          wsum         = wsum + proc_wdist(iproc)
          auxty        = 0_ip
          call maths_sfc_d2xyz_tab(proc_boxes(1),iproc,x,y,z,auxty)
          coo1d        = (z-1)*proc_boxes(1)*proc_boxes(1)+(y-1)*proc_boxes(1)+x
          reord(iproc) = coo1d
          types(iproc) = auxty
          invor(coo1d) = iproc
       enddo
    end if
    if(present(totalw)) wsum=totalw;
    if(wsum == 0.0_rp) return
    !
    ! Perform partition
    !
    iaux = 1_ip
    nullify(wrank)
    call memory_alloca(memor,'wrank',vacal,wrank,npar)

    do iproc=1,nproc

       nrank = 0
       !loc_npar(iproc) = (proc_wdist(reord(iproc))/wsum) * real(npar,rp)
       loc_npar(iproc) = (proc_wdist(iproc)/wsum) * real(npar,rp)
       if(loc_npar(iproc) > 0.0_rp) then

          irank(iproc) = iaux
          auxlg = .true.

          do while(auxlg .and. nrank < npar .and. iaux <= npar)
             loc_npar(iproc) = loc_npar(iproc) - dist_corr(iaux)
             if(loc_npar(iproc) < 0) then
                auxlg = .false.
                wrank(nrank+1)= loc_npar(iproc) + dist_corr(iaux)
                dist_corr(iaux) = -loc_npar(iproc)
             else
                wrank(nrank+1)=dist_corr(iaux)
                iaux = iaux + 1
             endif
             nrank = nrank + 1    ! si queda 0 el resido petara - repassar aquest cas
             if(loc_npar(iproc) == 0) auxlg = .false.
          end do
          if(iproc == nproc) then
             nrank = npar-irank(iproc)+1
          endif

          !loc_npar(iproc) = (proc_wdist(reord(iproc))/wsum) * real(npar,rp)
          loc_npar(iproc) = (proc_wdist(iproc)/wsum) * real(npar,rp)

          !if(invor(mypro) == iproc) then
          if(mypro == iproc) then
             call maths_sfc_part_2(ndime,loc_boxes,loc_weigh,loc_part,types(iproc),&
                  irank(iproc),nrank,wrank,PART_NAME)
             exit
          endif
       endif
    enddo
    !
    ! deallocate pointers
    !
    call memory_deallo(memor,'wrank',vacal,wrank)
    call memory_deallo(memor,'loc_npar',vacal,loc_npar)
    call memory_deallo(memor,'dist_corr',vacal,dist_corr)
    call memory_deallo(memor,'irank',vacal,irank)
    call memory_deallo(memor,'reord',vacal,reord)
    call memory_deallo(memor,'invor',vacal,invor)
    deallocate(types)

  end subroutine maths_sfc_par_part_2_rp  

  !-----------------------------------------------------------------------
  !>
  !> @date    08/01/2019
  !> @author  Guillaume Houzeaux
  !> @brief   Different elements in an array
  !> @details Extract the different elements in an array
  !>
  !-----------------------------------------------------------------------

  subroutine maths_list_different_elements(list_in,num_diff,list_diff,OPTION)

    integer(ip), intent(in),    pointer           :: list_in(:)
    integer(ip), intent(out)                      :: num_diff
    integer(ip), intent(inout), pointer, optional :: list_diff(:)
    character(*),                        optional :: OPTION
    logical(lg),                pointer           :: mask(:)
    integer(ip),                pointer           :: index_vector(:)
    integer(ip)                                   :: i,nn

    nn = size(list_in,KIND=ip)
    if( present(list_diff) ) nullify(list_diff)

    if( nn <= 0 .or. .not. associated(list_in) ) then

       num_diff = 0

    else

       allocate(mask(nn))
       mask = .true.

       do i = nn,2,-1
          mask(i) = .not.(any(list_in(:i-1)==list_in(i)))
       end do

       if( present(OPTION) ) then
          if( trim(OPTION) == 'STRICTLY POSITIVE' ) then
             where( list_in <= 0 ) mask = .false.
          end if
       end if

       ! Count

       num_diff = size(pack([(i,i=1,nn)],mask),KIND=ip)

       if( present(list_diff) .and. num_diff > 0 ) then
          !
          ! Make an index vector
          !
          allocate(index_vector(num_diff))
          index_vector=pack([(i,i=1,nn)],mask)
          !
          ! Now copy the unique elements of list_in to list_diff
          !
          allocate(list_diff(num_diff))

          list_diff=list_in(index_vector)
          deallocate(index_vector)

       end if

       deallocate(mask)

    end if

  end subroutine maths_list_different_elements

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    24/01/2018
  !> @brief   Units
  !> @details Compute unit to scale a variable
  !>
  !-----------------------------------------------------------------------

  subroutine maths_unit(unit_value,unit_char,unit_factor)
    implicit none
    real(rp),      intent(in)  :: unit_value  !< Memory in bytes
    character(1),  intent(out) :: unit_char   !< Memory unit character
    real(rp),      intent(out) :: unit_factor !< Memory scaling

    if(      unit_value >= 1.0e12_rp ) then
       unit_factor = 1.0e-12_rp
       unit_char   = 'P'
    else if( unit_value >= 1.0e9_rp  ) then
       unit_factor = 1.0e-9_rp
       unit_char   = 'G'
    else if( unit_value >= 1.0e6_rp  ) then
       unit_factor = 1.0e-6_rp
       unit_char   = 'M'
    else if( unit_value >= 1.0e3_rp  ) then
       unit_factor = 1.0e3_rp
       unit_char   = 'k'
    else if( unit_value >= 1.0_rp ) then
       unit_factor = 1.0_rp
       unit_char   = 'm'
    else if( unit_value >= 1.0e-3_rp ) then
       unit_factor = 1.0e3_rp
       unit_char   = 'm'
    else if( unit_value >= 1.0e-6_rp ) then
       unit_factor = 1.0e6_rp
       unit_char   = 'u'
    else if( unit_value >= 1.0e-9_rp ) then
       unit_factor = 1.0e9_rp
       unit_char   = 'n'
    else if( unit_value >= 1.0e-12_rp ) then
       unit_factor = 1.0e12_rp
       unit_char   = 'p'
    else
       unit_factor = 1.0_rp
       unit_char   = ' '
    end if

  end subroutine maths_unit

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-11-02
  !> @brief   Merge arrays
  !> @details Merge ordered arrays:
  !>          JA_IN  = [ 2 6 9 ]
  !>          JA_OUT = [ 3 6 8 10 ]
  !>          Result:
  !>          JA_OUT = [ 2 3 6 8 9 10 ]
  !>
  !-----------------------------------------------------------------------

  subroutine maths_merge_ordered_lists(nz_in,ja_in,nz_out,ja_out,MEMORY_COUNTER,RESIZE_GRAPH,LIST_NAME)

    integer(ip),      intent(in)              :: nz_in
    integer(ip),      intent(in)              :: ja_in(*)
    integer(ip),      intent(out)             :: nz_out
    integer(ip),      intent(inout), pointer  :: ja_out(:)
    integer(8),       intent(inout), optional :: MEMORY_COUNTER(2)
    logical(lg),      intent(in),    optional :: RESIZE_GRAPH
    character(len=*), intent(in),    optional :: LIST_NAME
    integer(ip)                               :: ipoin,iz_out,last_position,iz_in
    integer(8)                                :: memor(2)
    character(20)                             :: my_list_name

    if( present(LIST_NAME) ) then
       my_list_name = trim(LIST_NAME)
    else
       my_list_name = 'JA_OUT'
    end if
    
    if( present(MEMORY_COUNTER) ) then
       memor = MEMORY_COUNTER
    else
       memor = 0_8
    end if
    !
    ! Look for last first available position
    !
    nz_out = size(ja_out)
    last_position = nz_out
    if( ja_out(last_position) /= 0 ) then
       nz_out = int(1.5_rp*real(size(ja_out,KIND=ip),rp))
       call memory_resize(memor,trim(my_list_name),'maths_merge_ordered_lists',ja_out,nz_out)
    end if
    do while( ja_out(last_position) == 0 )
       last_position = last_position - 1
    end do
    last_position = last_position + 1
    iz_out = 1
    !
    ! Scan entries of JA_IN and possibly add it to JA_OUT
    !
    do iz_in = 1,nz_in
       ipoin = ja_in(iz_in)
       if( ipoin > 0 ) then
          !
          ! Cjeck if entry already exists
          !
          loop_ipoin: do while( ipoin > ja_out(iz_out) .and. ja_out(iz_out) /= 0 )
             iz_out = iz_out + 1
             if( iz_out > nz_out ) exit loop_ipoin
          end do loop_ipoin
          !
          ! Reallocate if necessay
          !
          if( iz_out > nz_out .or. last_position > nz_out ) then
             nz_out = int(1.5_rp*real(size(ja_out,KIND=ip),rp))
             call memory_resize(memor,trim(my_list_name),'maths_merge_ordered_lists',ja_out,nz_out)
          end if
          !
          ! Add new entry at position LAST_POSITION
          !
          if( ipoin /= ja_out(iz_out) ) then
             ja_out(last_position) = ipoin
             last_position = last_position + 1
          end if

       end if
    end do
    !
    ! Order JA_OUT
    !
    last_position = last_position - 1
    nz_out = last_position
    call maths_heap_sort(2_ip,last_position,ja_out)
    !
    ! Resize graph
    !
    if( present(RESIZE_GRAPH) ) then
       if( RESIZE_GRAPH ) &
            call memory_resize(memor,trim(my_list_name),'maths_merge_ordered_lists',ja_out,nz_out)
    end if

    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor

  end subroutine maths_merge_ordered_lists

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2019-03-11
  !> @brief   Solve a quedratic equation
  !> @details Find roots of a*x^2+b*x+c=0
  !>          All the computations can be done in 16-byte precision:
  !>          Solutions are nevertheless converted into 8-byte precision
  !>          at the end.
  !>
  !-----------------------------------------------------------------------

  subroutine maths_quadratic_equation(aa,bb,cc,x1,x2,num_solutions,QUAD_REAL)

    real(rp),              intent(in)  :: aa             !< Polynomial coef
    real(rp),              intent(in)  :: bb             !< Polynomial coef
    real(rp),              intent(in)  :: cc             !< Polynomial coef
    real(rp),              intent(out) :: x1             !< First root
    real(rp),              intent(out) :: x2             !< Second root
    integer(ip),           intent(out) :: num_solutions  !< Number of solutions
    logical(lg), optional, intent(in)  :: QUAD_REAL      !< If 16-bytes precision is required
    real(rp)                           :: delta

    real(qp)                           :: aa_qp
    real(qp)                           :: bb_qp
    real(qp)                           :: cc_qp
    real(qp)                           :: x1_qp
    real(qp)                           :: delta_qp
    logical(lg)                        :: if_quad_real

    if( present(QUAD_REAL) ) then
       if_quad_real = QUAD_REAL
    else
       if_quad_real = .false.
    end if

    if( if_quad_real ) then

       !----------------------------------------------------------------
       !
       ! Quad-real precision
       !
       !----------------------------------------------------------------

       aa_qp = real(aa,qp)
       bb_qp = real(bb,qp)
       cc_qp = real(cc,qp)
       if( abs(aa_qp) > 0.0_qp ) then
          delta_qp  = ((bb_qp/aa_qp)*bb_qp - 4.0_qp*cc_qp)/aa_qp
          if( delta_qp == 0.0_qp ) then
             !
             ! One real solution:
             !
             !      -b
             ! x1 = --
             !      2a
             !
             x1 = real(-bb_qp/(2.0_qp*aa_qp),rp)
             x2 = x1                             ! Ensures no garbage values are used for unit test
             num_solutions = 1_ip
          else if( delta_qp < 0.0_qp ) then
             !
             ! No real solution
             !
             num_solutions = 0_ip
          else
             !
             ! Two real solutions
             !
             delta_qp = sqrt(delta_qp)
             if( bb_qp >= 0.0_qp ) then
                x1_qp = 0.5_qp * (-bb_qp/aa_qp-delta_qp)
                x1    = real(x1_qp,rp)
                if( x1_qp /= 0.0_qp ) then
                   x2 = real(cc_qp / (aa_qp*x1_qp),rp)
                else
                   x2 = real(0.5_qp * (-bb_qp/aa_qp+delta_qp),rp)
                end if
             else
                x1_qp = 0.5_qp * (-bb_qp/aa_qp+delta_qp)
                x1    = real(x1_qp,rp)
                if( x1_qp /= 0.0_qp ) then
                   x2 = real(cc_qp / (aa_qp*x1_qp),rp)
                else
                   x2 = real(0.5_qp * (-bb_qp/aa_qp-delta_qp),rp)
                end if
             end if
             num_solutions = 2_ip
          end if

       else
          !
          ! a=0
          !
          if( abs(bb_qp) > 0.0_qp ) then
             x1 = real(-cc_qp/bb_qp,rp)
             x2 = x1
             num_solutions = 1_ip
          else
             num_solutions = 0_ip
          end if

       end if

    else

       !----------------------------------------------------------------
       !
       ! RP real precision
       !
       !----------------------------------------------------------------

       if( abs(aa) > 0.0_rp ) then
          delta  = ((bb/aa)*bb - 4.0_rp*cc)/aa
          if( delta == 0.0_rp ) then
             !
             ! One real solution:
             !
             !      -b
             ! x1 = --
             !      2a
             !
             x1 = -bb/(2.0_rp*aa)
             x2 = x1
             num_solutions = 1_ip
          else if( delta < 0.0_rp ) then
             !
             ! No real solution
             !
             num_solutions = 0_ip
          else
             !
             ! Two real solutions
             !
             delta = sqrt(delta)
             if( bb >= 0.0_rp ) then
                x1 = 0.5_rp * (-bb/aa-delta)
                if( x1 /= 0.0_rp ) then
                   x2 = cc / (aa*x1)
                else
                   x2 = 0.5_rp * (-bb/aa+delta)
                end if
             else
                x1 = 0.5_rp * (-bb/aa+delta)
                if( x1 /= 0.0_rp ) then
                   x2 = cc / (aa*x1)
                else
                   x2 = 0.5_rp * (-bb/aa-delta)
                end if
             end if
             num_solutions = 2_ip
          end if

       else
          !
          ! a=0
          !
          if( abs(bb) > 0.0_rp ) then
             x1 = -cc/bb
             x2 = x1
             num_solutions = 1_ip
          else
             num_solutions = 0_ip
          end if

       end if
    end if

  end subroutine maths_quadratic_equation

  subroutine maths_quadratic_equation_unity_test()

    real(rp)    :: aa,bb,cc,x1,x2,x3,x4,xtmp,delta
    real(qp)    :: xe1,xe2
    integer(ip) :: num1,num2,solution

    solution = 4

    select case ( solution )

    case ( 1 )

       aa  =  2.0_rp
       bb  = -3.0_rp
       cc  =  1.0_rp
       xe1 =  0.5_qp
       xe2 =  1.0_qp

    case ( 2 )

       aa  =  1.0_rp
       bb  =  200.0_rp
       cc  = -0.000015_rp
       xe1 =-200.000000075_qp
       xe2 =   0.000000075_qp

    case ( 3 )

       xe1 = 1.000000000000000_qp
       xe2 = 1.000000028975958_qp
       aa  = 94906265.625_rp
       bb  = -189812534.0_rp
       cc  = 94906268.375_rp

    case ( 4 )

       aa  = 1.0_rp
       bb  = 1.0e8_rp
       cc  = 1.0_rp
       xe1 = -99999999.9999999899999999999999990024_qp
       xe2 = -1.0000000000000001000000000000000e-0008_qp

    end select

    x1 = 0.0_rp
    x2 = 0.0_rp
    x3 = 0.0_rp
    x4 = 0.0_rp
    !
    ! Our algorithm
    !
    call maths_quadratic_equation(aa,bb,cc,x1,x2,num1) !,QUAD_REAL=.true.)
    !
    ! Original algorithm
    !
    delta = bb**2-4.0_rp*aa*cc
    if( delta > 0.0_rp ) then
       x3   = (-bb-sqrt(delta))/(2.0_rp*aa)
       x4   = (-bb+sqrt(delta))/(2.0_rp*aa)
       num2 = 2
    else if( delta == 0.0_rp ) then
       x3   = -bb/(2.0_rp*aa)
       x4   = 0.0_rp
       num2 = 1
    end if
    !
    ! Order solution
    !
    if( num2 == 2 ) then
       if( x4 < x3 ) then
          xtmp = x3
          x3   = x4
          x4   = xtmp
       end if
    end if
    if( num1 == 2 ) then
       if( x2 < x1 ) then
          xtmp = x1
          x1   = x2
          x2   = xtmp
       end if
    end if

    print*,' '
    print*,'Modified algorithm:'
    print*,'  - Number solution= ',num1
    print*,'  - Solution 1=      ',x1, abs(real(x1,qp)-xe1)/abs(xe1)*100.0_qp
    print*,'  - Solution 2=      ',x2, abs(real(x2,qp)-xe2)/abs(xe2)*100.0_qp
    print*,' '
    print*,'Original algorithm:'
    print*,'  - Number solution= ',num2
    print*,'  - Solution 1=      ',x3, abs(real(x3,qp)-xe1)/abs(xe1)*100.0_qp
    print*,'  - Solution 2=      ',x4, abs(real(x4,qp)-xe2)/abs(xe2)*100.0_qp

  end subroutine maths_quadratic_equation_unity_test

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2019-04-04
  !> @brief   Szudzik pairing function
  !> @details Provides a bijective 2 to 1 mapping
  !>
  !-----------------------------------------------------------------------

  integer(ip) elemental function maths_Szudzik_pairing_function(ii,jj)

    integer(ip), intent(in) :: ii
    integer(ip), intent(in) :: jj

    if( ii /= max(ii,jj) ) then
       maths_Szudzik_pairing_function = jj*jj+ii
    else
       maths_Szudzik_pairing_function = ii*(ii+1)+jj
    end if

  end function maths_Szudzik_pairing_function

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-07
  !> @brief   Distribute iterations
  !> @details Distribute n iterations in npart parts
  !> 
  !-----------------------------------------------------------------------

  subroutine maths_distribute_iterations(npart,niter,liter,DEFAULT_VALUE)

    integer(ip),            intent(in)    :: npart
    integer(ip),            intent(in)    :: niter
    integer(ip), pointer,   intent(inout) :: liter(:)
    integer(ip), optional,  intent(in)    :: DEFAULT_VALUE
    integer(ip)                           :: ipart,iiter

    if( niter > 0 ) then
       
       if( .not. associated(liter) ) call memory_alloca(memor,'liter',"maths_distribute_iterations",liter,niter)
       
       if( present(DEFAULT_VALUE) .and. npart == 0 ) then
          liter = DEFAULT_VALUE
       else
          iiter = 0
          do while( iiter < niter ) 
             ipart = 0
             do while( iiter < niter .and. ipart < npart )
                iiter        = iiter + 1
                ipart        = ipart + 1
                liter(iiter) = ipart
             end do
          end do
       end if
    end if
    
  end subroutine maths_distribute_iterations
      
end module mod_maths
!> @}
