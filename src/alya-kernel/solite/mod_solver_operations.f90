!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



  
  !-----------------------------------------------------------------------
  !> @addtogroup Solver
  !> @{
  !> @file    mod_solver_operations.f90
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Solver operations
  !> @details Solver operations
  !-----------------------------------------------------------------------
  
module mod_solver_operations

  use def_kintyp_basic,                   only : ip,rp,lg
  use def_mat,                            only : mat
  use def_mat_csr,                        only : mat_csr
  use def_kintyp_solvers,                 only : soltyp
  use mod_communications_point_to_point,  only : PAR_INTERFACE_EXCHANGE
  use mod_parall,                         only : PAR_COMM_MY_CODE_ARRAY
  use mod_optional_argument,              only : optional_argument
  use mod_solver,                         only : solver_parallel_double_scalar_product
  use mod_solver,                         only : solver_parallel_scalar_product
  use mod_solver,                         only : solver_parallel_vector_L2norm
  implicit none

  private

  public :: parallel_mv
  public :: parallel_L2norm
  public :: parallel_exchange
  public :: parallel_double_dot_product
  
contains
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Parallel MV
  !> @details Parallel MV
  !> 
  !-----------------------------------------------------------------------
  
  subroutine parallel_mv(&
       solve,aa,xx,yy,&
       MPI,OPENMP,INITIALIZATION,TIMING,&
       MY_TIMING,ASYNCHRONISM)

    type(soltyp),                    intent(inout) :: solve            !< Solver type
    class(mat),                      intent(in)    :: aa               !< Matrix type
    real(rp),               pointer, intent(inout) :: xx(:)            !< Multiplicand
    real(rp),               pointer, intent(inout) :: yy(:)            !< Result vector
    logical(lg),  optional,          intent(in)    :: MPI              !< If MPI should be used or not
    logical(lg),  optional,          intent(in)    :: OPENMP           !< If OpenMP should be used or not
    logical(lg),  optional,          intent(in)    :: INITIALIZATION   !< If result vector should be initialized
    logical(lg),  optional,          intent(in)    :: TIMING           !< If timing should be activated 
    real(rp),     optional,          intent(out)   :: MY_TIMING(2)     !< Register timing in this array
    logical(lg),  optional,          intent(in)    :: ASYNCHRONISM     !< Blocking or non-blocking communications
    integer(ip)                                    :: n1,n2,n3,n4
    logical(lg)                                    :: if_mpi

    if( associated(xx) .and. associated(yy) ) then

       if_mpi = optional_argument(.true.,MPI)
       n1     = 1
       n2     = solve % nequ1
       n3     = solve % nequ1+1 
       n4     = solve % nequa

       call aa % mv(xx,yy,n3,n4,INITIALIZATION,OPENMP,solve % omp_chunk_size,solve % omp_schedule)
       if( if_mpi ) call parallel_exchange(solve,yy)
       call aa % mv(xx,yy,n1,n2,INITIALIZATION,OPENMP,solve % omp_chunk_size,solve % omp_schedule)

    end if

  end subroutine parallel_mv

  subroutine parallel_exchange(solve,yy)
    
    type(soltyp),          intent(inout) :: solve            !< Solver type
    real(rp),     pointer, intent(inout) :: yy(:)            !< Result vector

    if( associated(yy) ) then
       call pararr('SOL',3_ip,solve % nunkn,yy)
       !call PAR_INTERFACE_EXCHANGE(solve % ndofn,yy,solve % comm)
       !call PAR_INTERFACE_EXCHANGE(solve % ndofn,yy,PAR_COMM_MY_CODE_ARRAY(1))
    end if
    
  end subroutine parallel_exchange
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Parallel MV
  !> @details Parallel MV
  !> 
  !-----------------------------------------------------------------------
  
  subroutine parallel_L2norm(solve,xx,xnorm)
  
    type(soltyp),                    intent(inout) :: solve            !< Solver type
    real(rp),               pointer, intent(inout) :: xx(:)            !< Multiplicand
    real(rp),                        intent(out)   :: xnorm
    real(rp)                                       :: xx_tmp(2)
    
    if( associated(xx) ) then
       call solver_parallel_vector_L2norm(solve,xx,xnorm)
    else
       call solver_parallel_vector_L2norm(solve,xx_tmp,xnorm)
    end if
  
  end subroutine parallel_L2norm

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Double dot product
  !> @details Double dot product
  !> 
  !-----------------------------------------------------------------------

  subroutine parallel_double_dot_product(solve,xx,yy,zz,xxdotyy,xxdotzz,MY_TIMING,OPENMP)
    type(soltyp),          intent(inout)         :: solve
    real(rp),     pointer, intent(in)            :: xx(:)
    real(rp),     pointer, intent(in)            :: yy(:)
    real(rp),     pointer, intent(in)            :: zz(:)
    real(rp),              intent(out)           :: xxdotyy
    real(rp),              intent(out)           :: xxdotzz
    real(rp),              intent(out), optional :: MY_TIMING(2)
    logical(lg),           intent(in),  optional :: OPENMP
    real(rp)                                     :: xx_tmp(2)
    
    if( associated(xx) ) then
       call solver_parallel_double_scalar_product(solve,xx,yy,zz,xxdotyy,xxdotzz,MY_TIMING,OPENMP)
    else
       call solver_parallel_double_scalar_product(solve,xx_tmp,xx_tmp,xx_tmp,xxdotyy,xxdotzz,MY_TIMING,OPENMP)
    end if

  end subroutine parallel_double_dot_product
  
end module mod_solver_operations
!> @}
