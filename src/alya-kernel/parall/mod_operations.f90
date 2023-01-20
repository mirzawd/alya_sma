!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @name    Parallelization toolbox
!> @file    mod_operations.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!> @brief   ToolBox for parallel operations
!> @details ToolBox for parallel operations
!>
!------------------------------------------------------------------------

module mod_operations
  use def_kintyp_comm,    only : comm_data_par
  use def_kintyp,         only : ip,rp,lg
  use def_domain,         only : npoin
  use def_master,         only : ISEQUEN,INOTMASTER
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_DEFINE_COMMUNICATOR
  use mod_memory,         only : memory_size
  use def_mpi
#include "def_mpi.inc"
  
  implicit none
  private

  interface operations_parallel_vector_L2norm
     module procedure operations_parallel_vector_L2norm_0,&
          &           operations_parallel_vector_L2norm_1,&
          &           operations_parallel_vector_L2norm_2
  end interface operations_parallel_vector_L2norm

  interface operations_parallel_dot_product
     module procedure operations_parallel_dot_product_0,&
          &           operations_parallel_dot_product_1,&
          &           operations_parallel_dot_product_2
  end interface operations_parallel_dot_product

  public :: operations_parallel_vector_L2norm
  public :: operations_parallel_dot_product

contains

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    06/03/2014
  !> @brief   Computes the 2-norm (Eucledian norm) of a vector XX
  !> @details Computes the 2-norm (Eucledian norm) of a vector XX
  !
  !----------------------------------------------------------------------

  subroutine operations_parallel_vector_L2norm_0(ndofn,xx,sumxx,wherein)
    implicit none
    integer(ip),                   intent(in)  :: ndofn
    real(rp),                      intent(in)  :: xx(ndofn,*) !< Vector
    real(rp),                      intent(out) :: sumxx       !< Eucledian norm
    character(*),        optional, intent(in)  :: wherein
    type(comm_data_par), pointer               :: commu
    MY_MPI_COMM                                :: PAR_COMM_TO_USE
    integer(ip)                                :: ii

    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
    end if

    sumxx = 0.0_rp

    if( ISEQUEN ) then
       if( ndofn == 1 ) then
          sumxx = dot_product(xx(1,1:npoin),xx(1,1:npoin))
       else
          do ii = 1,ndofn
             sumxx = sumxx + dot_product(xx(ii,1:npoin),xx(ii,1:npoin))
          end do
       end if
    else if( INOTMASTER ) then
       if( ndofn == 1 ) then
          sumxx = dot_product(xx(1,1:commu % npoi3),xx(1,1:commu % npoi3))
       else
          do ii = 1,ndofn
             sumxx = sumxx + dot_product(xx(ii,1:commu % npoi3),xx(ii,1:commu % npoi3))
          end do
       end if
    end if

    if( present(wherein) ) then
       call PAR_SUM(sumxx,wherein)
    else
       call PAR_SUM(sumxx)
    end if

    sumxx = sqrt(sumxx)

  end subroutine operations_parallel_vector_L2norm_0

  subroutine operations_parallel_vector_L2norm_1(xx,sumxx,wherein)
    implicit none
    real(rp),            pointer,  intent(in)  :: xx(:)       !< Vector
    real(rp),                      intent(out) :: sumxx       !< Eucledian norm
    character(*),        optional, intent(in)  :: wherein
    type(comm_data_par), pointer               :: commu
    MY_MPI_COMM                                :: PAR_COMM_TO_USE

    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
    end if

    sumxx = 0.0_rp

    if( ISEQUEN ) then
       sumxx = dot_product(xx(1:npoin),xx(1:npoin))
    else if( INOTMASTER ) then
       sumxx = dot_product(xx(1:commu % npoi3),xx(1:commu % npoi3))
    end if

    if( present(wherein) ) then
       call PAR_SUM(sumxx,wherein)
    else
       call PAR_SUM(sumxx)
    end if

    sumxx = sqrt(sumxx)

  end subroutine operations_parallel_vector_L2norm_1

  subroutine operations_parallel_vector_L2norm_2(xx,sumxx,wherein)
    implicit none
    real(rp),            pointer,  intent(in)  :: xx(:,:)     !< Vector
    real(rp),                      intent(out) :: sumxx       !< Eucledian norm
    character(*),        optional, intent(in)  :: wherein
    type(comm_data_par), pointer               :: commu
    MY_MPI_COMM                                :: PAR_COMM_TO_USE
    integer(ip)                                :: ii

    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
    end if

    sumxx = 0.0_rp

    if( ISEQUEN ) then
       if( memory_size(xx,1_ip) == 1 ) then
          sumxx = dot_product(xx(1,1:npoin),xx(1,1:npoin))
       else
          do ii = 1,memory_size(xx,1_ip)
             sumxx = sumxx + dot_product(xx(ii,1:npoin),xx(ii,1:npoin))
          end do
       end if
    else if( INOTMASTER ) then
       if( memory_size(xx,1_ip) == 1 ) then
          sumxx = dot_product(xx(1,1:commu % npoi3),xx(1,1:commu % npoi3))
       else
          do ii = 1,memory_size(xx,1_ip)
             sumxx = sumxx + dot_product(xx(ii,1:commu % npoi3),xx(ii,1:commu % npoi3))
          end do
       end if
    end if

    if( present(wherein) ) then
       call PAR_SUM(sumxx,wherein)
    else
       call PAR_SUM(sumxx)
    end if

    sumxx = sqrt(sumxx)

  end subroutine operations_parallel_vector_L2norm_2

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    06/03/2014
  !> @brief   Dot product
  !> @details Dot product of two vectors x and y
  !
  !----------------------------------------------------------------------

  subroutine operations_parallel_dot_product_0(ndofn,xx,yy,xdoty,wherein)
    implicit none
    integer(ip),                   intent(in)  :: ndofn
    real(rp),                      intent(in)  :: xx(ndofn,*) !< Vector
    real(rp),                      intent(in)  :: yy(ndofn,*) !< Vector
    real(rp),                      intent(out) :: xdoty       !< Eucledian norm
    character(*),        optional, intent(in)  :: wherein
    type(comm_data_par), pointer               :: commu
    MY_MPI_COMM                                :: PAR_COMM_TO_USE
    integer(ip)                                :: ii

    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
    end if

    xdoty = 0.0_rp

    if( ISEQUEN ) then
       if( ndofn == 1 ) then
          xdoty = dot_product(xx(1,1:npoin),yy(1,1:npoin))
       else
          do ii = 1,ndofn
             xdoty = xdoty + dot_product(xx(ii,1:npoin),yy(ii,1:npoin))
          end do
       end if
    else if( INOTMASTER ) then
       if( ndofn == 1 ) then
          xdoty = dot_product(xx(1,1:commu % npoi3),yy(1,1:commu % npoi3))
       else
          do ii = 1,ndofn
             xdoty = xdoty + dot_product(xx(ii,1:commu % npoi3),yy(ii,1:commu % npoi3))
          end do
       end if
    end if

    if( present(wherein) ) then
       call PAR_SUM(xdoty,wherein)
    else
       call PAR_SUM(xdoty)
    end if

  end subroutine operations_parallel_dot_product_0

  subroutine operations_parallel_dot_product_1(xx,yy,xdoty,wherein)
    implicit none
    real(rp),            pointer,  intent(in)  :: xx(:)       !< Vector
    real(rp),            pointer,  intent(in)  :: yy(:)       !< Vector
    real(rp),                      intent(out) :: xdoty       !< Eucledian norm
    character(*),        optional, intent(in)  :: wherein
    type(comm_data_par), pointer               :: commu
    MY_MPI_COMM                                :: PAR_COMM_TO_USE

    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
    end if

    xdoty = 0.0_rp

    if( associated(xx) .and. associated(yy) ) then
       if( ISEQUEN ) then
          xdoty = dot_product(xx(1:npoin),yy(1:npoin))
       else if( INOTMASTER ) then
          xdoty = dot_product(xx(1:commu % npoi3),yy(1:commu % npoi3))
       end if
    end if

    if( present(wherein) ) then
       call PAR_SUM(xdoty,wherein)
    else
       call PAR_SUM(xdoty)
    end if

  end subroutine operations_parallel_dot_product_1

  subroutine operations_parallel_dot_product_2(xx,yy,xdoty,wherein)
    implicit none
    real(rp),            pointer,  intent(in)  :: xx(:,:)     !< Vector
    real(rp),            pointer,  intent(in)  :: yy(:,:)     !< Vector
    real(rp),                      intent(out) :: xdoty       !< Eucledian norm
    character(*),        optional, intent(in)  :: wherein
    type(comm_data_par), pointer               :: commu
    MY_MPI_COMM                                :: PAR_COMM_TO_USE
    integer(ip)                                :: ii

    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
    end if

    xdoty = 0.0_rp

    if( associated(xx) .and. associated(yy) ) then
       if( ISEQUEN ) then
          if( memory_size(xx,1_ip) == 1 ) then
             xdoty = dot_product(xx(1,1:npoin),yy(1,1:npoin))
          else
             do ii = 1,memory_size(xx,1_ip)
                xdoty = xdoty + dot_product(xx(ii,1:npoin),yy(ii,1:npoin))
             end do
          end if
       else if( INOTMASTER ) then
          if( memory_size(xx,1_ip) == 1 ) then
             xdoty = dot_product(xx(1,1:commu % npoi3),yy(1,1:commu % npoi3))
          else
             do ii = 1,memory_size(xx,1_ip)
                xdoty = xdoty + dot_product(xx(ii,1:commu % npoi3),yy(ii,1:commu % npoi3))
             end do
          end if
       end if
    end if

    if( present(wherein) ) then
       call PAR_SUM(xdoty,wherein)
    else
       call PAR_SUM(xdoty)
    end if

  end subroutine operations_parallel_dot_product_2

end module mod_operations
!> @}
