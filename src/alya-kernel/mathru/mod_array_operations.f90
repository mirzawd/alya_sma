!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @addtogroup Array_Operations_Toolbox
!> Toolbox for array operations, like axpy, dot product, imitating BLAS.
!> @{
!> @name    ToolBox for array operations
!> @file    mod_array_operations.f90
!> @date    22/05/2015
!> @author  Guillaume Houzeaux
!> @brief   Array operations ("a la" BLAS)
!> @details The following subroutines are available:
!>          \verbatim
!>          AXPY ........ y = y + alpha * x
!>          AXPBY ....... y = beta * y + alpha * x
!>          COPY ........ y = x
!>          AXOZ ........ y = alpha * x (*,/,+,-) z
!>          CONST ....... y = alpha
!>          \endverbatim
!>          The following functions are available:
!>          \verbatim
!>          NORM2 ....... alpha = sqrt(x.x)
!>          DOT ......... alpha = x.y
!>          \endverbatim
!>          Different operations are also coded, like the residual norm
!>          or the min and max of an array
!
!-----------------------------------------------------------------------

module mod_array_operations

  use def_kintyp,         only : ip,rp,lg
  use def_master,         only : zeror
  use def_domain,         only : npoin_own
  use def_master,         only : INOTMASTER
  use def_master,         only : IPARALL
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_MAX
  implicit none
  private 

  interface array_operations_residual_norm
     module procedure array_operations_residual_norm_1_1,&
          &           array_operations_residual_norm_1_2,&
          &           array_operations_residual_norm_2_2,&
          &           array_operations_residual_norm_1_3,&
          &           array_operations_residual_norm_3_3,&
          &           array_operations_residual_norm_3_1
  end interface array_operations_residual_norm

  interface array_operations_min_max
     module procedure array_operations_min_max_rp_1,&
          &           array_operations_min_max_rp_2,&
          &           array_operations_min_max_rp_3
  end interface array_operations_min_max

  interface array_operations_axpy
     module procedure array_operations_axpy_1     , &
          &           array_operations_axpy_2     , &
          &           array_operations_axpy_1n
  end interface array_operations_axpy
  interface array_operations_axpby
     module procedure array_operations_axpby_1    , &
          &           array_operations_axpby_2    , &
          &           array_operations_axpby_3    , &
          &           array_operations_axpby_1n   , &
          &           array_operations_axpby_2n   , &
          &           array_operations_axpby_12n  , &
          &           array_operations_axpby_21n 
  end interface array_operations_axpby

  interface array_operations_copy
     module procedure array_operations_copy_11  , &
          &           array_operations_copy_22  , &
          &           array_operations_copy_33  , &
          &           array_operations_copy_12  , &
          &           array_operations_copy_21  , &
          &           array_operations_copy_23  , &
          &           array_operations_copy_p1  , &
          &           array_operations_copy_p2    
  end interface array_operations_copy

  interface array_operations_axoz
     module procedure array_operations_axoz_11  , &
          &           array_operations_axoz_p1
  end interface array_operations_axoz

  interface array_operations_const
     module procedure array_operations_const_p1  , &
          &           array_operations_const_p2
  end interface array_operations_const

  interface array_operations_norm2
     module procedure array_operations_norm2_p1  , &
          &           array_operations_norm2_p2
  end interface array_operations_norm2

  interface array_operations_dot
     module procedure array_operations_dot_p1    , &
          &           array_operations_dot_p2
  end interface array_operations_dot

  interface array_operations_initialization
     module procedure array_operations_initialization_1, &
          &           array_operations_initialization_2, &
          &           array_operations_initialization_3
  end interface array_operations_initialization

  public :: array_operations_axpy
  public :: array_operations_axpby
  public :: array_operations_copy
  public :: array_operations_axoz
  public :: array_operations_const
  public :: array_operations_norm2
  public :: array_operations_dot
  public :: array_operations_initialization
  public :: array_operations_residual_norm
  public :: array_operations_min_max

contains

  !-----------------------------------------------------------------------
  !
  !> @brief   x = 0
  !> @details Array initialization
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine array_operations_initialization_1(xx)
    real(rp),   pointer :: xx(:)
    integer(ip)         :: ii

    if( associated(xx) ) then
       do ii = 1,size(xx,KIND=ip)
          xx(ii) = 0.0_rp
       end do
    end if
    
  end subroutine array_operations_initialization_1

  subroutine array_operations_initialization_2(xx)
    real(rp),   pointer :: xx(:,:)
    integer(ip)         :: ii,jj

    if( associated(xx) ) then
       do jj = 1,size(xx,2,KIND=ip)
          do ii = 1,size(xx,1,KIND=ip)
             xx(ii,jj) = 0.0_rp
          end do
       end do
    end if
    
  end subroutine array_operations_initialization_2

  subroutine array_operations_initialization_3(xx)
    real(rp),   pointer :: xx(:,:,:)
    integer(ip)         :: ii,jj,kk

    if( associated(xx) ) then
       do kk = 1,size(xx,3,KIND=ip)
          do jj = 1,size(xx,2,KIND=ip)
             do ii = 1,size(xx,1,KIND=ip)
                xx(ii,jj,kk) = 0.0_rp
             end do
          end do
       end do
    end if
    
  end subroutine array_operations_initialization_3

  !-----------------------------------------------------------------------
  !
  !> @brief   y = y + alpha*x
  !> @details Array operation using or not using OMP: y = y + alpha*x
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine array_operations_axpy_generic(nn,alpha,xx,yy,use_omp)
    integer(ip), intent(in)    :: nn
    real(rp),    intent(in)    :: alpha
    real(rp),    intent(in)    :: xx(*)
    real(rp),    intent(inout) :: yy(*)
    logical(lg), intent(in)    :: use_omp
    integer(ip)                :: ii

    if( INOTMASTER ) then
       if( use_omp ) then
          !$OMP PARALLEL DO SCHEDULE (STATIC)  &
          !$OMP DEFAULT  ( NONE )              &
          !$OMP SHARED   ( xx, yy, alpha, nn ) &
          !$OMP PRIVATE  ( ii )            
          do ii = 1,nn
             yy(ii) = yy(ii) + alpha * xx(ii)
          end do
          !$OMP END PARALLEL DO
       else
          do ii = 1,nn
             yy(ii) = yy(ii) + alpha * xx(ii)
          end do
       end if
    end if

  end subroutine array_operations_axpy_generic

  subroutine array_operations_axpy_1(alpha,xx,yy,worder)
    real(rp),               intent(in)    :: alpha
    real(rp),     pointer,  intent(in)    :: xx(:)
    real(rp),     pointer,  intent(inout) :: yy(:)
    character(*), optional, intent(in)    :: worder
    integer(ip)                           :: nn
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( size(xx,1,KIND=ip) /= size(yy,1,KIND=ip) ) &
            call runend('ARRAY_OPERATIONS_AXPY: WRONG DIMENSIONS')
       nn = size(xx,KIND=ip)       
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_axpy_generic(nn,alpha,xx,yy,use_omp)
    end if

  end subroutine array_operations_axpy_1

  subroutine array_operations_axpy_1n(nn,alpha,xx,yy,worder)
    integer(ip),            intent(in)    :: nn
    real(rp),               intent(in)    :: alpha
    real(rp),               intent(in)    :: xx(*)
    real(rp),               intent(inout) :: yy(*)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_axpy_generic(nn,alpha,xx,yy,use_omp)
    end if

  end subroutine array_operations_axpy_1n

  subroutine array_operations_axpy_2(alpha,xx,yy,worder)
    real(rp),     pointer,  intent(in)    :: xx(:,:)
    real(rp),     pointer,  intent(inout) :: yy(:,:)
    real(rp),               intent(in)    :: alpha
    character(*), optional, intent(in)    :: worder
    integer(ip)                           :: nn
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       nn = size(xx,1,KIND=ip)*size(xx,2,KIND=ip)       
       if( size(xx,1,KIND=ip) /= size(yy,1,KIND=ip) .or. size(xx,2,KIND=ip) /= size(yy,2,KIND=ip) ) &
            call runend('ARRAY_OPERATIONS_AXPY: WRONG DIMENSIONS')
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_axpy_generic(nn,alpha,xx,yy,use_omp)
    end if

  end subroutine array_operations_axpy_2

  !-----------------------------------------------------------------------
  !
  !> @brief   y = alpha*x + beta*y
  !> @details Array operation using or not using OMP: y = alpha*x + beta*y
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine array_operations_axpby_generic(nn,alpha,beta,xx,yy,use_omp)
    integer(ip), intent(in)    :: nn
    real(rp),    intent(in)    :: alpha
    real(rp),    intent(in)    :: beta
    real(rp),    intent(in)    :: xx(*)
    real(rp),    intent(inout) :: yy(*)
    logical(lg), intent(in)    :: use_omp
    integer(ip)                :: ii

    if( use_omp ) then
       !$OMP PARALLEL DO SCHEDULE (STATIC)        &
       !$OMP DEFAULT  ( NONE )                    &
       !$OMP SHARED   ( xx, yy, alpha, beta, nn ) &
       !$OMP PRIVATE  ( ii )            
       do ii = 1,nn
          yy(ii) = alpha * xx(ii) + beta * yy(ii)
       end do
       !$OMP END PARALLEL DO
    else
       do ii = 1,nn
          yy(ii) = alpha * xx(ii) + beta * yy(ii)
       end do
    end if

  end subroutine array_operations_axpby_generic

  subroutine array_operations_axpby_1(alpha,beta,xx,yy,worder)
    real(rp),               intent(in)    :: alpha
    real(rp),               intent(in)    :: beta
    real(rp),     pointer,  intent(in)    :: xx(:)
    real(rp),     pointer,  intent(inout) :: yy(:)
    character(*), optional, intent(in)    :: worder
    integer(ip)                           :: nn
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( size(xx,1,KIND=ip) /= size(yy,1,KIND=ip) ) &
            call runend('ARRAY_OPERATIONS_AXPY: WRONG DIMENSIONS')
       nn = size(xx,KIND=ip)       
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_axpby_generic(nn,alpha,beta,xx,yy,use_omp)
    end if

  end subroutine array_operations_axpby_1

  subroutine array_operations_axpby_2(alpha,beta,xx,yy,worder)
    real(rp),               intent(in)    :: alpha
    real(rp),               intent(in)    :: beta
    real(rp),     pointer,  intent(in)    :: xx(:,:)
    real(rp),     pointer,  intent(inout) :: yy(:,:)
    character(*), optional, intent(in)    :: worder
    integer(ip)                           :: nn
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       nn = size(xx,1,KIND=ip)*size(xx,2,KIND=ip)       
       if( size(xx,1,KIND=ip) /= size(yy,1,KIND=ip) .or. size(xx,2,KIND=ip) /= size(yy,2,KIND=ip) ) &
            call runend('ARRAY_OPERATIONS_AXPY: WRONG DIMENSIONS')
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_axpby_generic(nn,alpha,beta,xx,yy,use_omp)
    end if

  end subroutine array_operations_axpby_2

  subroutine array_operations_axpby_3(alpha,beta,xx,yy,worder)
    real(rp),               intent(in)    :: alpha
    real(rp),               intent(in)    :: beta
    real(rp),     pointer,  intent(in)    :: xx(:,:,:)
    real(rp),     pointer,  intent(inout) :: yy(:,:,:)
    character(*), optional, intent(in)    :: worder
    integer(ip)                           :: nn
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       nn = size(xx,1,KIND=ip)*size(xx,2,KIND=ip)       
       if( size(xx,1,KIND=ip) /= size(yy,1,KIND=ip) .or. size(xx,2,KIND=ip) /= size(yy,2,KIND=ip) ) &
            call runend('ARRAY_OPERATIONS_AXPY: WRONG DIMENSIONS')
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_axpby_generic(nn,alpha,beta,xx,yy,use_omp)
    end if

  end subroutine array_operations_axpby_3

  subroutine array_operations_axpby_1n(n1,alpha,beta,xx,yy,worder)
    integer(ip),            intent(in)    :: n1
    real(rp),               intent(in)    :: alpha
    real(rp),               intent(in)    :: beta
    real(rp),               intent(in)    :: xx(*)
    real(rp),               intent(inout) :: yy(*)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_axpby_generic(n1,alpha,beta,xx,yy,use_omp)
    end if

  end subroutine array_operations_axpby_1n

  subroutine array_operations_axpby_2n(n1,n2,alpha,beta,xx,yy,worder)
    integer(ip),            intent(in)    :: n1
    integer(ip),            intent(in)    :: n2
    real(rp),               intent(in)    :: alpha
    real(rp),               intent(in)    :: beta
    real(rp),               intent(in)    :: xx(n1,*)
    real(rp),               intent(inout) :: yy(n1,*)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_axpby_generic(n1*n2,alpha,beta,xx,yy,use_omp)
    end if

  end subroutine array_operations_axpby_2n

  subroutine array_operations_axpby_12n(n1,n2,alpha,beta,xx,yy,worder)
    integer(ip),            intent(in)    :: n1
    integer(ip),            intent(in)    :: n2
    real(rp),               intent(in)    :: xx(*)
    real(rp),               intent(inout) :: yy(n2,*)
    real(rp),               intent(in)    :: alpha
    real(rp),               intent(in)    :: beta
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_axpby_generic(n1*n2,alpha,beta,xx,yy,use_omp)
    end if

  end subroutine array_operations_axpby_12n

  subroutine array_operations_axpby_21n(n1,n2,alpha,beta,xx,yy,worder)
    integer(ip),            intent(in)    :: n1
    integer(ip),            intent(in)    :: n2
    real(rp),               intent(in)    :: xx(n1,*)
    real(rp),               intent(inout) :: yy(*)
    real(rp),               intent(in)    :: alpha
    real(rp),               intent(in)    :: beta
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_axpby_generic(n1*n2,alpha,beta,xx,yy,use_omp)
    end if

  end subroutine array_operations_axpby_21n

  !-----------------------------------------------------------------------
  !
  !> @brief   y = x
  !> @details Array operation using or not using OMP: y = x
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine array_operations_copy_generic(nn,xx,yy,use_omp)
    integer(ip), intent(in)    :: nn
    real(rp),    intent(in)    :: xx(*)
    real(rp),    intent(inout) :: yy(*)
    logical(lg), intent(in)    :: use_omp
    integer(ip)                :: ii

    if( use_omp ) then
       !$OMP PARALLEL DO SCHEDULE (STATIC)  &
       !$OMP DEFAULT  ( NONE )              &
       !$OMP SHARED   ( xx, yy, nn )        &
       !$OMP PRIVATE  ( ii )            
       do ii = 1,nn
          yy(ii) = xx(ii)
       end do
       !$OMP END PARALLEL DO
    else
       do ii = 1,nn
          yy(ii) = xx(ii)
       end do
    end if

  end subroutine array_operations_copy_generic

  subroutine array_operations_copy_11(n1,xx,yy,worder)
    integer(ip),            intent(in)    :: n1
    real(rp),               intent(in)    :: xx(*)
    real(rp),               intent(inout) :: yy(*)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_copy_generic(n1,xx,yy,use_omp)
    end if

  end subroutine array_operations_copy_11

  subroutine array_operations_copy_22(n1,n2,xx,yy,worder)
    integer(ip),            intent(in)    :: n1
    integer(ip),            intent(in)    :: n2
    real(rp),               intent(in)    :: xx(n1,*)
    real(rp),               intent(out)   :: yy(n1,*)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_copy_generic(n1*n2,xx,yy,use_omp)
    end if

  end subroutine array_operations_copy_22

  subroutine array_operations_copy_33(n1,n2,n3,xx,yy,worder)
    integer(ip),            intent(in)    :: n1
    integer(ip),            intent(in)    :: n2
    integer(ip),            intent(in)    :: n3
    real(rp),               intent(in)    :: xx(n1,n2,*)
    real(rp),               intent(inout) :: yy(n1,n2,*)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_copy_generic(n1*n2*n3,xx,yy,use_omp)
    end if

  end subroutine array_operations_copy_33

  subroutine array_operations_copy_12(n1,n2,xx,yy,worder)
    integer(ip),            intent(in)    :: n1
    integer(ip),            intent(in)    :: n2
    real(rp),               intent(in)    :: xx(*)
    real(rp),               intent(inout) :: yy(n1,*)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_copy_generic(n1*n2,xx,yy,use_omp)
    end if

  end subroutine array_operations_copy_12

  subroutine array_operations_copy_21(n1,n2,xx,yy,worder)
    integer(ip),            intent(in)    :: n1
    integer(ip),            intent(in)    :: n2
    real(rp),               intent(in)    :: xx(n1,*)
    real(rp),               intent(inout) :: yy(*)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_copy_generic(n1*n2,xx,yy,use_omp)
    end if

  end subroutine array_operations_copy_21

  subroutine array_operations_copy_23(n1,n2,xx,yy,worder)
    integer(ip),            intent(in)    :: n1
    integer(ip),            intent(in)    :: n2
    real(rp),               intent(in)    :: xx(n1,n2,*)
    real(rp),               intent(inout) :: yy(n1,n2,*)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_copy_generic(n1*n2,xx,yy,use_omp)
    end if

  end subroutine array_operations_copy_23

  subroutine array_operations_copy_p1(xx,yy,worder)
    real(rp),     pointer,  intent(in)    :: xx(:)
    real(rp),     pointer,  intent(inout) :: yy(:)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp
    integer(ip)                           :: nn

    if( INOTMASTER ) then
       nn = size(xx,1,KIND=ip)
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_copy_generic(nn,xx,yy,use_omp)
    end if

  end subroutine array_operations_copy_p1

  subroutine array_operations_copy_p2(xx,yy,worder)
    real(rp),     pointer,   intent(in)    :: xx(:,:)
    real(rp),     pointer,   intent(inout) :: yy(:,:)
    character(*), optional,  intent(in)    :: worder
    logical(lg)                            :: use_omp
    integer(ip)                            :: n1,n2

    if( INOTMASTER ) then
       n1 = size(xx,1,KIND=ip)
       n2 = size(xx,2,KIND=ip)
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_copy_generic(n1*n2,xx,yy,use_omp)
    end if

  end subroutine array_operations_copy_p2

  !-----------------------------------------------------------------------
  !
  !> @brief   y = alpha*x (*,/,+,-) z
  !> @details Array operation using or not using OMP: y = alpha*x (*,/,+,-) z
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine array_operations_axoz_generic(nn,alpha,woperation,xx,zz,yy,use_omp)
    integer(ip),  intent(in)    :: nn
    real(rp),     intent(in)    :: alpha
    character(1), intent(in)    :: woperation
    real(rp),     intent(in)    :: xx(*)
    real(rp),     intent(in)    :: zz(*)
    real(rp),     intent(inout) :: yy(*)
    logical(lg),  intent(in)    :: use_omp
    integer(ip)                 :: ii

    if( INOTMASTER ) then

       if( woperation == '*' ) then
          !
          ! y = alpha * x * z
          !
          if( use_omp ) then
             !$OMP PARALLEL DO SCHEDULE (STATIC)      &
             !$OMP DEFAULT  ( NONE )                  &
             !$OMP SHARED   ( xx, zz, yy, alpha, nn ) &
             !$OMP PRIVATE  ( ii )            
             do ii = 1,nn
                yy(ii) = alpha * xx(ii) * zz(ii)
             end do
             !$OMP END PARALLEL DO
          else
             do ii = 1,nn
                yy(ii) = alpha * xx(ii) * zz(ii)
             end do
          end if

       else if( woperation == '/' ) then
          !
          ! y = alpha * x / z
          !
          if( use_omp ) then
             !$OMP PARALLEL DO SCHEDULE (STATIC)      &
             !$OMP DEFAULT  ( NONE )                  &
             !$OMP SHARED   ( xx, zz, yy, alpha, nn ) &
             !$OMP PRIVATE  ( ii )            
             do ii = 1,nn
                yy(ii) = alpha * xx(ii) / zz(ii)
             end do
             !$OMP END PARALLEL DO
          else
             do ii = 1,nn
                yy(ii) = alpha * xx(ii) / zz(ii)
             end do
          end if

       else if( woperation == '+' ) then
          !
          ! y = alpha * x + z
          !
          if( use_omp ) then
             !$OMP PARALLEL DO SCHEDULE (STATIC)      &
             !$OMP DEFAULT  ( NONE )                  &
             !$OMP SHARED   ( xx, zz, yy, alpha, nn ) &
             !$OMP PRIVATE  ( ii )            
             do ii = 1,nn
                yy(ii) = alpha * xx(ii) + zz(ii)
             end do
             !$OMP END PARALLEL DO
          else
             do ii = 1,nn
                yy(ii) = alpha * xx(ii) + zz(ii)
             end do
          end if

       else if( woperation == '-' ) then
          !
          ! y = alpha * x - z
          !
          if( use_omp ) then
             !$OMP PARALLEL DO SCHEDULE (STATIC)      &
             !$OMP DEFAULT  ( NONE )                  &
             !$OMP SHARED   ( xx, zz, yy, alpha, nn ) &
             !$OMP PRIVATE  ( ii )            
             do ii = 1,nn
                yy(ii) = alpha * xx(ii) - zz(ii)
             end do
             !$OMP END PARALLEL DO
          else
             do ii = 1,nn
                yy(ii) = alpha * xx(ii) - zz(ii)
             end do
          end if

       end if

    end if

  end subroutine array_operations_axoz_generic

  subroutine array_operations_axoz_11(nn,alpha,woperation,xx,zz,yy,worder)
    integer(ip),            intent(in)    :: nn
    real(rp),               intent(in)    :: alpha
    character(1),           intent(in)    :: woperation
    real(rp),               intent(in)    :: xx(*)
    real(rp),               intent(in)    :: zz(*)
    real(rp),               intent(inout) :: yy(*)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_axoz_generic(nn,alpha,woperation,xx,zz,yy,use_omp)
    end if

  end subroutine array_operations_axoz_11

  subroutine array_operations_axoz_p1(alpha,woperation,xx,zz,yy,worder)
    real(rp),               intent(in)    :: alpha
    character(1),           intent(in)    :: woperation
    real(rp),     pointer,  intent(in)    :: xx(:)
    real(rp),     pointer,  intent(in)    :: zz(:)
    real(rp),     pointer,  intent(inout) :: yy(:)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp
    integer(ip)                           :: nn

    if( INOTMASTER ) then
       nn = size(xx,KIND=ip)
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       call array_operations_axoz_generic(nn,alpha,woperation,xx,zz,yy,use_omp)
    end if

  end subroutine array_operations_axoz_p1

  subroutine array_operations_const_p1(alpha,yy,worder)
    real(rp),               intent(in)    :: alpha
    real(rp),     pointer,  intent(inout) :: yy(:)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp
    integer(ip)                           :: ii

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       if( use_omp ) then
          !$OMP PARALLEL DO SCHEDULE (STATIC) &
          !$OMP DEFAULT  ( NONE )             &
          !$OMP SHARED   ( yy, alpha )        &
          !$OMP PRIVATE  ( ii )            
          do ii = 1,size(yy,KIND=ip)
             yy(ii) = alpha 
          end do
          !$OMP END PARALLEL DO
       else
          do ii = 1,size(yy,KIND=ip)
             yy(ii) = alpha 
          end do
       end if
    end if

  end subroutine array_operations_const_p1

  subroutine array_operations_const_p2(alpha,yy,worder)
    real(rp),               intent(in)    :: alpha
    real(rp),     pointer,  intent(inout) :: yy(:,:)
    character(*), optional, intent(in)    :: worder
    logical(lg)                           :: use_omp
    integer(ip)                           :: ii,jj,n1,n2

    if( INOTMASTER ) then
       if( present(worder) ) then
          if( worder == 'USE OPENMP' ) then
             use_omp = .true.
          else if( worder == 'DO NOT USE OPENMP' ) then
             use_omp = .false.
          end if
       else
          use_omp = .true.
       end if
       n1 = size(yy,1,KIND=ip)
       n2 = size(yy,2,KIND=ip)
       if( use_omp ) then
          !$OMP PARALLEL DO SCHEDULE (STATIC)  &
          !$OMP DEFAULT  ( NONE )              &
          !$OMP SHARED   ( yy, alpha, n1, n2 ) &
          !$OMP PRIVATE  ( ii, jj )            
          do jj = 1,n2
             do ii = 1,n1
                yy(ii,jj) = alpha
             end do
          end do
          !$OMP END PARALLEL DO
       else
          do jj = 1,n2
             do ii = 1,n1
                yy(ii,jj) = alpha
             end do
          end do
       end if
    end if

  end subroutine array_operations_const_p2

  subroutine array_operations_norm2_generic(nn,yy,yynorm2,use_omp)
    integer(ip), intent(in)  :: nn
    real(rp),    intent(in)  :: yy(*)
    real(rp),    intent(out) :: yynorm2
    logical(lg), intent(in)  :: use_omp
    integer(ip)              :: ii

    yynorm2 = 0.0_rp

    if( INOTMASTER ) then    
       if( use_omp ) then
          !$OMP PARALLEL DO SCHEDULE (STATIC) &
          !$OMP DEFAULT  ( NONE )             &
          !$OMP PRIVATE  ( ii )               &
          !$OMP SHARED   ( yy, nn )           &
          !$OMP REDUCTION (+:yynorm2)  
          do ii = 1,nn
             yynorm2 = yynorm2 + yy(ii) * yy(ii)
          end do
          !$OMP END PARALLEL DO
       else
          do ii = 1,nn
             yynorm2 = yynorm2 + yy(ii) * yy(ii)
          end do
       end if
    end if
    call PAR_SUM(yynorm2,'IN MY CODE')
    yynorm2 = sqrt(yynorm2)

  end subroutine array_operations_norm2_generic

  function array_operations_norm2_p1(yy,worder,DOFS)
    real(rp),     pointer,  intent(in) :: yy(:) 
    character(*), optional, intent(in) :: worder
    integer(ip),  optional, intent(in) :: DOFS
    real(rp)                           :: array_operations_norm2_p1
    logical(lg)                        :: use_omp
    integer(ip)                        :: nn
    real(rp), dimension(1)             :: yy_tmp

    if( present(worder) ) then
       if( worder == 'USE OPENMP' ) then
          use_omp = .true.
       else if( worder == 'DO NOT USE OPENMP' ) then
          use_omp = .false.
       end if
    else
       use_omp = .true.
    end if

    if( associated(yy) )then
       if( present(DOFS) ) then
          nn = npoin_own*DOFS
       else
          nn = npoin_own
       end if
       call array_operations_norm2_generic(nn,yy,array_operations_norm2_p1,use_omp)
    else
       nn = 0_ip
       yy_tmp = 0.0_rp
       call array_operations_norm2_generic(nn,yy_tmp,array_operations_norm2_p1,use_omp)
    endif
    
  end function array_operations_norm2_p1

  function array_operations_norm2_p2(yy,worder,DOFS)
    real(rp),     pointer,  intent(in) :: yy(:,:) 
    character(*), optional, intent(in) :: worder
    integer(ip),  optional, intent(in) :: DOFS
    real(rp)                           :: array_operations_norm2_p2
    logical(lg)                        :: use_omp
    integer(ip)                        :: nn
    real(rp), dimension(1)             :: yy_tmp

    if( present(worder) ) then
       if( worder == 'USE OPENMP' ) then
          use_omp = .true.
       else if( worder == 'DO NOT USE OPENMP' ) then
          use_omp = .false.
       end if
    else
       use_omp = .true.
    end if

    if( associated(yy) )then
       if( present(DOFS) ) then
          nn = npoin_own*DOFS
       else
          nn = npoin_own*size(yy,1,KIND=ip)
       end if
       call array_operations_norm2_generic(nn,yy,array_operations_norm2_p2,use_omp)
    else
        nn = 0_ip
        yy_tmp = 0.0_rp
        call array_operations_norm2_generic(nn,yy_tmp,array_operations_norm2_p2,use_omp)
    endif

  end function array_operations_norm2_p2

  subroutine array_operations_dot_generic(nn,xx,yy,xdoty,use_omp)
    integer(ip), intent(in)  :: nn
    real(rp),    intent(in)  :: xx(*)
    real(rp),    intent(in)  :: yy(*)
    real(rp),    intent(out) :: xdoty
    logical(lg), intent(in)  :: use_omp
    integer(ip)              :: ii

    xdoty = 0.0_rp

    if( INOTMASTER ) then    
       if( use_omp ) then
          !$OMP PARALLEL DO SCHEDULE (STATIC) &
          !$OMP DEFAULT  ( NONE )             &
          !$OMP PRIVATE  ( ii )               &
          !$OMP SHARED   ( yy, nn )           &
          !$OMP REDUCTION (+:xdoty)  
          do ii = 1,nn
             xdoty = xdoty + xx(ii) * yy(ii)
          end do
          !$OMP END PARALLEL DO
       else
          do ii = 1,nn
             xdoty = xdoty + xx(ii) * yy(ii)
          end do
       end if
    end if
    call PAR_SUM(xdoty,'IN MY CODE')

  end subroutine array_operations_dot_generic

  function array_operations_dot_p1(xx,yy,worder)
    real(rp),     pointer,  intent(in) :: xx(:) 
    real(rp),     pointer,  intent(in) :: yy(:) 
    character(*), optional, intent(in) :: worder
    real(rp)                           :: array_operations_dot_p1
    logical(lg)                        :: use_omp
    integer(ip)                        :: nn

    if( present(worder) ) then
       if( worder == 'USE OPENMP' ) then
          use_omp = .true.
       else if( worder == 'DO NOT USE OPENMP' ) then
          use_omp = .false.
       end if
    else
       use_omp = .true.
    end if
    nn = size(yy,KIND=ip)
    call array_operations_dot_generic(nn,xx,yy,array_operations_dot_p1,use_omp)

  end function array_operations_dot_p1

  function array_operations_dot_p2(xx,yy,worder)
    real(rp),     pointer,  intent(in) :: xx(:,:) 
    real(rp),     pointer,  intent(in) :: yy(:,:) 
    character(*), optional, intent(in) :: worder
    real(rp)                           :: array_operations_dot_p2
    logical(lg)                        :: use_omp
    integer(ip)                        :: nn

    if( present(worder) ) then
       if( worder == 'USE OPENMP' ) then
          use_omp = .true.
       else if( worder == 'DO NOT USE OPENMP' ) then
          use_omp = .false.
       end if
    else
       use_omp = .true.
    end if
    nn = size(yy,KIND=ip)
    call array_operations_dot_generic(nn,xx,yy,array_operations_dot_p2,use_omp)

  end function array_operations_dot_p2


  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-02
  !> @brief   Norms of residual
  !> @details Compute the L2, L1, Linf difference (with relaxation) between
  !>          two vectors: 
  !>          redif = || [r*v1+(1-r)*v2] - v2|| / ||r*v1+(1-r)*v2||  
  !>
  !>          Arrays are of the type vi(nni,npoin,*)
  !>          kcomi is the offset to the array i to start looping. The loop
  !>          over nni starts at (ii-1)*nni + kcomi
  !>
  !-----------------------------------------------------------------------  

  real(rp) function array_operations_residual_norm_1_1(norm,nn1,nn2,v1,v2,kcom1,kcom2,kdime,RELAXATION)

    integer(ip), intent(in)           :: norm       !< Norm type (Linf, L1, L2)
    integer(ip), intent(in)           :: nn1        !< Dimension of array v1 (1D, 2D, 3D)
    integer(ip), intent(in)           :: nn2        !< Dimension of array v2 (1D, 2D, 3D)
    real(rp),    intent(in), pointer  :: v1(:)      !< First array
    real(rp),    intent(in), pointer  :: v2(:)      !< Second array
    integer(ip), intent(in)           :: kcom1      !< Component (iteration or time step)
    integer(ip), intent(in)           :: kcom2      !< Component (iteration or time step)
    integer(ip), intent(in)           :: kdime      !< Dimensions to compute the norm over
    real(rp),    intent(in), optional :: RELAXATION !< Relaxation
    real(rp)                          :: redif      !< Norm
    real(rp)                          :: v1_tmp(2)
    real(rp)                          :: v2_tmp(2)
    real(rp)                          :: relax

    if( present(RELAXATION) ) then
       relax = RELAXATION
    else
       relax = 1.0_rp
    end if
    
    if( ( .not. associated(v1) ) .and. ( .not. associated(v2) ) ) then
       call array_operations_residual_norm_all(norm,nn1,nn2,v1_tmp,v2_tmp,kcom1,kcom2,kdime,relax,redif)  
    else if( ( .not. associated(v1) ) ) then
       call array_operations_residual_norm_all(norm,nn1,nn2,v1_tmp,v2,kcom1,kcom2,kdime,relax,redif)  
    else if( ( .not. associated(v2) ) ) then
       call array_operations_residual_norm_all(norm,nn1,nn2,v1,v2_tmp,kcom1,kcom2,kdime,relax,redif)  
    else
       call array_operations_residual_norm_all(norm,nn1,nn2,v1,v2,kcom1,kcom2,kdime,relax,redif)  
    end if

    array_operations_residual_norm_1_1 = redif
    
  end function array_operations_residual_norm_1_1

  real(rp) function array_operations_residual_norm_1_2(norm,nn1,nn2,v1,v2,kcom1,kcom2,kdime,RELAXATION)

    integer(ip), intent(in)           :: norm       !< Norm type (Linf, L1, L2)
    integer(ip), intent(in)           :: nn1        !< Dimension of array v1 (1D, 2D, 3D)
    integer(ip), intent(in)           :: nn2        !< Dimension of array v2 (1D, 2D, 3D)
    real(rp),    intent(in), pointer  :: v1(:)      !< First array
    real(rp),    intent(in), pointer  :: v2(:,:)    !< Second array
    integer(ip), intent(in)           :: kcom1      !< Component (iteration or time step)
    integer(ip), intent(in)           :: kcom2      !< Component (iteration or time step)
    integer(ip), intent(in)           :: kdime      !< Dimensions to compute the norm over
    real(rp),    intent(in), optional :: RELAXATION !< Relaxation
    real(rp)                          :: redif      !< Norm
    real(rp)                          :: v1_tmp(2)
    real(rp)                          :: v2_tmp(2)
    real(rp)                          :: relax

    if( present(RELAXATION) ) then
       relax = RELAXATION
    else
       relax = 1.0_rp
    end if
    
    if( ( .not. associated(v1) ) .and. ( .not. associated(v2) ) ) then
       call array_operations_residual_norm_all(norm,nn1,nn2,v1_tmp,v2_tmp,kcom1,kcom2,kdime,relax,redif)  
    else if( ( .not. associated(v1) ) ) then
       call array_operations_residual_norm_all(norm,nn1,nn2,v1_tmp,v2,kcom1,kcom2,kdime,relax,redif)  
    else if( ( .not. associated(v2) ) ) then
       call array_operations_residual_norm_all(norm,nn1,nn2,v1,v2_tmp,kcom1,kcom2,kdime,relax,redif)  
    else
       call array_operations_residual_norm_all(norm,nn1,nn2,v1,v2,kcom1,kcom2,kdime,relax,redif)  
    end if

    array_operations_residual_norm_1_2 = redif
    
  end function array_operations_residual_norm_1_2

   real(rp) function array_operations_residual_norm_2_2(norm,nn1,nn2,v1,v2,kcom1,kcom2,kdime,RELAXATION)

    integer(ip), intent(in)           :: norm       !< Norm type (Linf, L1, L2)
    integer(ip), intent(in)           :: nn1        !< Dimension of array v1 (1D, 2D, 3D)
    integer(ip), intent(in)           :: nn2        !< Dimension of array v2 (1D, 2D, 3D)
    real(rp),    intent(in), pointer  :: v1(:,:)    !< First array
    real(rp),    intent(in), pointer  :: v2(:,:)    !< Second array
    integer(ip), intent(in)           :: kcom1      !< Component (iteration or time step)
    integer(ip), intent(in)           :: kcom2      !< Component (iteration or time step)
    integer(ip), intent(in)           :: kdime      !< Dimensions to compute the norm over
    real(rp),    intent(in), optional :: RELAXATION !< Relaxation
    real(rp)                          :: redif      !< Norm
    real(rp)                          :: v1_tmp(2)
    real(rp)                          :: v2_tmp(2)
    real(rp)                          :: relax

    if( present(RELAXATION) ) then
       relax = RELAXATION
    else
       relax = 1.0_rp
    end if
    
    if( ( .not. associated(v1) ) .and. ( .not. associated(v2) ) ) then
       call array_operations_residual_norm_all(norm,nn1,nn2,v1_tmp,v2_tmp,kcom1,kcom2,kdime,relax,redif)  
    else if( ( .not. associated(v1) ) ) then
       call array_operations_residual_norm_all(norm,nn1,nn2,v1_tmp,v2,kcom1,kcom2,kdime,relax,redif)  
    else if( ( .not. associated(v2) ) ) then
       call array_operations_residual_norm_all(norm,nn1,nn2,v1,v2_tmp,kcom1,kcom2,kdime,relax,redif)  
    else
       call array_operations_residual_norm_all(norm,nn1,nn2,v1,v2,kcom1,kcom2,kdime,relax,redif)  
    end if

    array_operations_residual_norm_2_2 = redif
    
  end function array_operations_residual_norm_2_2

  real(rp) function array_operations_residual_norm_1_3(norm,nn1,nn2,v1,v2,kcom1,kcom2,kdime,RELAXATION)

    integer(ip), intent(in)           :: norm       !< Norm type (Linf, L1, L2)
    integer(ip), intent(in)           :: nn1        !< Dimension of array v1 (1D, 2D, 3D)
    integer(ip), intent(in)           :: nn2        !< Dimension of array v2 (1D, 2D, 3D)
    real(rp),    intent(in), pointer  :: v1(:)      !< First array
    real(rp),    intent(in), pointer  :: v2(:,:,:)  !< Second array
    integer(ip), intent(in)           :: kcom1      !< Component (iteration or time step)
    integer(ip), intent(in)           :: kcom2      !< Component (iteration or time step)
    integer(ip), intent(in)           :: kdime      !< Dimensions to compute the norm over
    real(rp),    intent(in), optional :: RELAXATION !< Relaxation
    real(rp)                          :: redif      !< Norm
    real(rp)                          :: v1_tmp(2)
    real(rp)                          :: v2_tmp(2)
    real(rp)                          :: relax

    if( present(RELAXATION) ) then
       relax = RELAXATION
    else
       relax = 1.0_rp
    end if
    
    if( ( .not. associated(v1) ) .and. ( .not. associated(v2) ) ) then
       call array_operations_residual_norm_all(norm,nn1,nn2,v1_tmp,v2_tmp,kcom1,kcom2,kdime,relax,redif)  
    else if( ( .not. associated(v1) ) ) then
       call array_operations_residual_norm_all(norm,nn1,nn2,v1_tmp,v2,kcom1,kcom2,kdime,relax,redif)  
    else if( ( .not. associated(v2) ) ) then
       call array_operations_residual_norm_all(norm,nn1,nn2,v1,v2_tmp,kcom1,kcom2,kdime,relax,redif)  
    else
       call array_operations_residual_norm_all(norm,nn1,nn2,v1,v2,kcom1,kcom2,kdime,relax,redif)  
    end if

    array_operations_residual_norm_1_3 = redif
    
  end function array_operations_residual_norm_1_3

  real(rp) function array_operations_residual_norm_3_1(norm,nn1,nn2,v1,v2,kcom1,kcom2,kdime,RELAXATION)

    integer(ip), intent(in)           :: norm       !< Norm type (Linf, L1, L2)
    integer(ip), intent(in)           :: nn1        !< Dimension of array v1 (1D, 2D, 3D)
    integer(ip), intent(in)           :: nn2        !< Dimension of array v2 (1D, 2D, 3D)
    real(rp),    intent(in), pointer  :: v1(:,:,:)  !< First array
    real(rp),    intent(in), pointer  :: v2(:)      !< Second array
    integer(ip), intent(in)           :: kcom1      !< Component (iteration or time step)
    integer(ip), intent(in)           :: kcom2      !< Component (iteration or time step)
    integer(ip), intent(in)           :: kdime      !< Dimensions to compute the norm over
    real(rp),    intent(in), optional :: RELAXATION !< Relaxation
    real(rp)                          :: redif      !< Norm
    real(rp)                          :: v1_tmp(2)
    real(rp)                          :: v2_tmp(2)
    real(rp)                          :: relax

    if( present(RELAXATION) ) then
       relax = RELAXATION
    else
       relax = 1.0_rp
    end if
    
    if( ( .not. associated(v1) ) .and. ( .not. associated(v2) ) ) then
       call array_operations_residual_norm_all(norm,nn1,nn2,v1_tmp,v2_tmp,kcom1,kcom2,kdime,relax,redif)  
    else if( ( .not. associated(v1) ) ) then
       call array_operations_residual_norm_all(norm,nn1,nn2,v1_tmp,v2,kcom1,kcom2,kdime,relax,redif)  
    else if( ( .not. associated(v2) ) ) then
       call array_operations_residual_norm_all(norm,nn1,nn2,v1,v2_tmp,kcom1,kcom2,kdime,relax,redif)  
    else
       call array_operations_residual_norm_all(norm,nn1,nn2,v1,v2,kcom1,kcom2,kdime,relax,redif)  
    end if

    array_operations_residual_norm_3_1 = redif
    
  end function array_operations_residual_norm_3_1

  real(rp) function array_operations_residual_norm_3_3(norm,nn1,nn2,v1,v2,kcom1,kcom2,kdime,RELAXATION)

    integer(ip), intent(in)           :: norm       !< Norm type (Linf, L1, L2)
    integer(ip), intent(in)           :: nn1        !< Dimension of array v1 (1D, 2D, 3D)
    integer(ip), intent(in)           :: nn2        !< Dimension of array v2 (1D, 2D, 3D)
    real(rp),    intent(in), pointer  :: v1(:,:,:)  !< First array
    real(rp),    intent(in), pointer  :: v2(:,:,:)  !< Second array
    integer(ip), intent(in)           :: kcom1      !< Component (iteration or time step)
    integer(ip), intent(in)           :: kcom2      !< Component (iteration or time step)
    integer(ip), intent(in)           :: kdime      !< Dimensions to compute the norm over
    real(rp),    intent(in), optional :: RELAXATION !< Relaxation
    real(rp)                          :: redif      !< Norm
    real(rp)                          :: v1_tmp(2)
    real(rp)                          :: v2_tmp(2)
    real(rp)                          :: relax

    if( present(RELAXATION) ) then
       relax = RELAXATION
    else
       relax = 1.0_rp
    end if
    
    if( ( .not. associated(v1) ) .and. ( .not. associated(v2) ) ) then
       call array_operations_residual_norm_all(norm,nn1,nn2,v1_tmp,v2_tmp,kcom1,kcom2,kdime,relax,redif)  
    else if( ( .not. associated(v1) ) ) then
       call array_operations_residual_norm_all(norm,nn1,nn2,v1_tmp,v2,kcom1,kcom2,kdime,relax,redif)  
    else if( ( .not. associated(v2) ) ) then
       call array_operations_residual_norm_all(norm,nn1,nn2,v1,v2_tmp,kcom1,kcom2,kdime,relax,redif)  
    else
       call array_operations_residual_norm_all(norm,nn1,nn2,v1,v2,kcom1,kcom2,kdime,relax,redif)  
    end if

    array_operations_residual_norm_3_3 = redif
    
  end function array_operations_residual_norm_3_3

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-17
  !> @brief   Residual norm
  !> @details Compute the L2, L1, Linf difference (with relaxation) between
  !>          two vectors: 
  !>          redif = || [r*v1+(1-r)*v2] - v2|| / ||r*v1+(1-r)*v2||  
  !> 
  !>          V1(NN1,*)
  !>          V2(NN2,*)
  !>          residual over the KDIME dimensions, starting from KCOM1 and KCOM2.
  !> 
  !-----------------------------------------------------------------------

  subroutine array_operations_residual_norm_all(norm,nn1,nn2,v1,v2,kcom1,kcom2,kdime,relax,redif)

    integer(ip), intent(in)  :: nn1
    integer(ip), intent(in)  :: nn2
    integer(ip), intent(in)  :: norm
    integer(ip), intent(in)  :: kcom1
    integer(ip), intent(in)  :: kcom2
    integer(ip), intent(in)  :: kdime
    real(rp),    intent(in)  :: v1(*)
    real(rp),    intent(in)  :: v2(*)
    real(rp),    intent(in)  :: relax
    real(rp),    intent(out) :: redif
    integer(ip)              :: ii,i1,i2,idime
    real(rp)                 :: va,vo,vd,ra,ro,resin(2)
    real(rp)                 :: resid(2)
    integer(ip)              :: istar1
    integer(ip)              :: istar2

    istar1 = kcom1 + 1_ip
    istar2 = kcom2 + 1_ip

    ra    = relax
    ro    = 1.0_rp-ra
    resid = 0.0_rp
    
    select case( norm )

    case( 0_ip )
       !
       ! Linf norm
       !
       do ii = 1,npoin_own
          i1    = (ii-1)*nn1 + istar1
          i2    = (ii-1)*nn2 + istar2
          resin = 0.0_rp
          do idime = 0,kdime-1
             vo       = v2(i2+idime)
             va       = ra*v1(i1+idime)+ro*vo
             vd       = va-vo
             resin(1) = resin(1) + vd * vd
             resin(2) = resin(2) + va * va                 
          end do
          resin    = sqrt(resin)
          resid(1) = max(resid(1),resin(1))
          resid(2) = max(resid(2),resin(2))              
       end do

       call PAR_MAX(2_ip,resid)

       if( resid(2) > zeror ) then
          redif = resid(1)/resid(2)
       else
          redif = resid(1)
       end if

    case ( 1_ip )
       !
       ! L1 norm
       !
       do ii = 1,npoin_own
          i1 = (ii-1)*nn1 + istar1
          i2 = (ii-1)*nn2 + istar2
          do idime = 0,kdime-1                 
             vo       = v2(i2+idime)
             va       = ra*v1(i1+idime)+ro*vo                              
             resid(1) = resid(1) + abs(va-vo)
             resid(2) = resid(2) + abs(va)
          end do
       end do
       
       call PAR_SUM(2_ip,resid) 

       if( resid(2) > zeror ) then
          redif = resid(1)/resid(2)
       else
          redif = resid(1)
       end if

    case ( 2_ip )
       !
       ! L2 norm
       !
       if( ra == 1.0_rp ) then
          do ii = 1,npoin_own
             i1 = (ii-1)*nn1 + istar1
             i2 = (ii-1)*nn2 + istar2             
             do idime = 0,kdime-1 
                va       = v1(i1+idime)
                vd       = va-v2(i2+idime)
                resid(1) = resid(1) + vd*vd
                resid(2) = resid(2) + va*va
             end do
          end do
       else
          do ii = 1,npoin_own
             i1 = (ii-1)*nn1 + istar1
             i2 = (ii-1)*nn2 + istar2
             do idime = 0,kdime-1                 
                vo       = v2(i2+idime)
                va       = ra*v1(i1+idime)+ro*vo    
                vd       = va-vo
                resid(1) = resid(1) + vd*vd
                resid(2) = resid(2) + va*va
             end do
          end do
       end if

       call PAR_SUM(2_ip,resid) 

       if( resid(2) > zeror ) then
          redif = sqrt(resid(1)/resid(2))
       else
          redif = sqrt(resid(1))
       end if

    case default

       redif = 0.0_rp
       
    end select

  end subroutine array_operations_residual_norm_all

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-20
  !> @brief   Compute min and max of an array
  !> @details Compute min and max of an array
  !> 
  !-----------------------------------------------------------------------

  subroutine array_operations_min_max_rp_1(ndim1,ndim2,ndim3,vecto,vemin,vemax)
    
    integer(ip),          intent(in)    :: ndim1,ndim2,ndim3
    real(rp),    pointer, intent(inout) :: vecto(:)
    real(rp),             intent(out)   :: vemin,vemax
    real(rp)                            :: vecto_tmp(1)
    
    if( associated(vecto) ) then
       call array_operations_min_max_all(ndim1,ndim2,ndim3,vecto,vemin,vemax)
    else
       call array_operations_min_max_all(ndim1,ndim2,ndim3,vecto_tmp,vemin,vemax)
    end if
    
  end subroutine array_operations_min_max_rp_1
  
  subroutine array_operations_min_max_rp_2(ndim1,ndim2,ndim3,vecto,vemin,vemax)
    
    integer(ip),          intent(in)    :: ndim1,ndim2,ndim3
    real(rp),    pointer, intent(inout) :: vecto(:,:)
    real(rp),             intent(out)   :: vemin,vemax
    real(rp)                            :: vecto_tmp(1,1)
    
    if( associated(vecto) ) then
       call array_operations_min_max_all(ndim1,ndim2,ndim3,vecto,vemin,vemax)
    else
       call array_operations_min_max_all(ndim1,ndim2,ndim3,vecto_tmp,vemin,vemax)
    end if
    
  end subroutine array_operations_min_max_rp_2
  
  subroutine array_operations_min_max_rp_3(ndim1,ndim2,ndim3,vecto,vemin,vemax)
    
    integer(ip),          intent(in)    :: ndim1,ndim2,ndim3
    real(rp),    pointer, intent(inout) :: vecto(:,:,:)
    real(rp),             intent(out)   :: vemin,vemax
    real(rp)                            :: vecto_tmp(1,1,1)
    
    if( associated(vecto) ) then
       call array_operations_min_max_all(ndim1,ndim2,ndim3,vecto,vemin,vemax)
    else
       call array_operations_min_max_all(ndim1,ndim2,ndim3,vecto_tmp,vemin,vemax)
    end if
    
  end subroutine array_operations_min_max_rp_3
  
  subroutine array_operations_min_max_all(ndim1,ndim2,ndim3,vecto,vemin,vemax)

    integer(ip), intent(in)  :: ndim1,ndim2,ndim3
    real(rp),    intent(in)  :: vecto(ndim1,ndim2)
    real(rp),    intent(out) :: vemin,vemax
    integer(ip)              :: idim1,idim2,ndim4
    real(rp)                 :: uvalu
    real(rp)                 :: vminmax(2)

    vemin= huge(1.0_rp)
    vemax=-huge(1.0_rp)

    if(ndim1==1) then
       do idim2=1,ndim2
          if(vecto(1,idim2)>vemax) vemax=vecto(1,idim2)
          if(vecto(1,idim2)<vemin) vemin=vecto(1,idim2)
       end do
    else
       ndim4=min(ndim1,ndim3)
       do idim2=1,ndim2
          uvalu=0.0_rp
          do idim1=1,ndim4
             uvalu=uvalu+vecto(idim1,idim2)*vecto(idim1,idim2)
          end do
          uvalu=sqrt(uvalu)
          if(uvalu>vemax) vemax=uvalu
          if(uvalu<vemin) vemin=uvalu
       end do
    end if

    if( IPARALL ) then
       vminmax(1) = -vemin
       vminmax(2) =  vemax
       call PAR_MAX(2_ip,vminmax)
       vemin = -vminmax(1) 
       vemax =  vminmax(2)
    end if

  end subroutine array_operations_min_max_all
  
end module mod_array_operations
!> @}
