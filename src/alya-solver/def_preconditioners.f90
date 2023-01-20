!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




!-----------------------------------------------------------------------
!> @addtogroup Solvers
!> @{
!> @file    def_iterative_solvers.f90
!> @author  houzeaux
!> @date    2022-03-01
!> @brief   Preconditioners
!> @details Preconditioners
!>         
!-----------------------------------------------------------------------

module def_preconditioners 
  
  use def_kintyp_basic,                  only : ip,rp,lg
  use def_mat,                           only : mat
  use def_mat_dia,                       only : mat_dia
  use def_solvers,                       only : solver
  use def_iterative_solvers,             only : iterative_solver
  use def_iterative_solvers,             only : preconditioner
  use def_direct_solvers,                only : direct_solver
  use mod_memory_basic,                  only : memory_alloca
  use mod_memory_basic,                  only : memory_deallo
  use mod_memory_basic,                  only : memory_size
  use mod_maths_basic,                   only : maths_inverse
  use mod_optional_argument,             only : optional_argument
  use mod_communications_point_to_point, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications_global,         only : PAR_SUM
  use def_all_solvers
  use def_solver
  implicit none

  real(rp),      parameter :: epsil=epsilon(1.0_rp)
  character(19), parameter :: vacal = 'def_preconditioners'

  private

  public :: preconditioner_setup
  public :: preconditioner_alloca

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-07
  !> @brief   Set preconditioner
  !> @details Set preconditioner
  !> 
  !-----------------------------------------------------------------------

  subroutine preconditioner_alloca(solve)

    class(solver),                   intent(inout) :: solve
    class(preconditioner),   pointer               :: preco

    select type ( solve )
    class is ( iterative_solver )

       if( solve % input % kfl_preco /= SOL_NO_PRECOND ) then
          if( .not. associated(solve % preco) ) then
             allocate(solve % preco(1))
          end if
          preco => solve % preco(1)
          
          select case ( solve % input % kfl_preco )
             
          case ( SOL_SOLVER_DIAGONAL       ) ; allocate( gauss_elimination :: preco % p)
          case ( SOL_SOLVER_SOR            ) ; allocate( sor               :: preco % p)
          case ( SOL_SOLVER_SSOR           ) ; allocate( ssor              :: preco % p)
          case ( SOL_SOLVER_APPROX_INVERSE ) ; allocate( approx_inv        :: preco % p)
          case ( SOL_SOLVER_LINELET        ) ; allocate( linelet           :: preco % p)
          case ( SOL_SOLVER_RAS            ) ; allocate( ras               :: preco % p)
          case default                       ; call runend('OK: UNKNOWN PRECONDITIONER')
          end select
          
          call preco % p % init()
          preco % p % comm              => solve % comm
          preco % p % comm_is_allocated = .false.

       end if
       
    end select

  end subroutine preconditioner_alloca
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-07
  !> @brief   Set preconditioner
  !> @details Set preconditioner
  !> 
  !-----------------------------------------------------------------------

  recursive subroutine preconditioner_setup(solve,a)

    class(solver),                             intent(inout) :: solve
    class(mat),                        target, intent(inout) :: a
    class(preconditioner),             pointer               :: preco

    select type ( solve )
    class is ( iterative_solver )
       if( associated(solve % preco) ) then

          preco => solve % preco(1)
          preco % a_is_allocated            = .false.
          preco % p % input % kfl_symm      = solve % input % kfl_symm
          preco % p % input % verbose       = solve % input % verbose
          preco % p % input % kfl_full_rows = solve % input % kfl_full_rows
          preco % p % iwrite                = solve % iwrite

          select case ( solve % input % kfl_preco )

          case ( SOL_SOLVER_DIAGONAL )
             !
             ! Diagonal preconditioning
             !
             preco % a_is_allocated = .true.
             allocate(mat_dia :: preco % a)
             call preco % p % dim(a)
             call preco % p % setup(a)
             call preco % a % init          ()
             call a     % diag              (preco % a)
             call solve % parallel_exchange (preco % a)
             select type ( v => preco % a )
             class is ( mat_dia ) ; call v % inverse()
             end select

          case ( SOL_SOLVER_SOR, SOL_SOLVER_SSOR ) 
             !
             ! Gauss-Seidel (SOR) and SSOR preconditionings
             !
             preco % a => a
             preco % p % input % miter = -1
             preco % p % input % relax =  1.0_rp
             call preco % p % dim(a)
             call preco % p % setup(a)

          case ( SOL_SOLVER_APPROX_INVERSE ) 
             !
             ! Approximate inverse preconditioning
             !
             preco % a => a
             preco % p % input % levels_sparsity = solve % input % levels_sparsity 
             call preco % p % dim(a)
             call preco % p % setup(a)


          case ( SOL_SOLVER_RAS ) 
             !
             ! RAS preconditioning
             !
             preco % a => a
             preco % p % input % kfl_block_ras   = solve % input % kfl_block_ras
             preco % p % input % num_subd_ras    = solve % input % num_subd_ras
             call preco % p % dim(a)
             call preco % p % setup(a)

          case default

             preco % a => a
             call preco % p % dim(a)
             call preco % p % setup(a)

          end select

       end if

    end select
    
  end subroutine preconditioner_setup
  
end module def_preconditioners
!> @}
