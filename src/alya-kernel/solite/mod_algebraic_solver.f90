!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solver
!> @{
!> @file    mod_algebraic_solver.f90
!> @author  guillaume
!> @date    2021-03-29
!> @brief   Solvers
!> @details Driver for algebraic solvers
!-----------------------------------------------------------------------

module mod_algebraic_solver

  use def_kintyp_basic, only : ip,rp
  use mod_matrix,       only : matrix_CSR_SpMV
  use mod_memory_basic, only : memory_alloca
  use mod_memory_basic, only : memory_deallo
  use mod_block_solver, only : block_solver
  use mod_solver,       only : solver_preprocess_0
  use def_solver
  implicit none

  private

  public :: algebraic_solver
  
contains
  
  subroutine algebraic_solver(solve,A,b,u,M)

    type(soltyp),                                 intent(inout) :: solve
    real(rp),               contiguous, pointer,  intent(inout) :: A(:)
    real(rp),               contiguous, pointer,  intent(inout) :: b(:)
    real(rp),               contiguous, pointer,  intent(inout) :: u(:)
    real(rp),     optional, contiguous, pointer,  intent(in)    :: M(:)
    real(rp)                                                    :: dummr

    call solver_preprocess_0(solve,A,b,u)
    
    select case ( solve % kfl_block )

    case ( 0_ip )

       !-----------------------------------------------------------------
       !
       ! Normal solvers
       !
       !-----------------------------------------------------------------

       !call solver_old_to_new(solve,A,As)
       if( present(M) ) then
         if( associated(M) ) then
           call solver(b,u,A,M)
         else
           call solver(b,u,A,dummr)
         end if
       else
         call solver(b,u,A,dummr)
       end if
       
    case ( 1_ip )

       !-----------------------------------------------------------------
       !
       ! Block solvers for systems
       !
       !-----------------------------------------------------------------

       call block_solver(solve,A,b,u,M)
       
    end select

  end subroutine algebraic_solver

end module mod_algebraic_solver
!> @}
  
