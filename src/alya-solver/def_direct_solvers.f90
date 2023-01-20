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
!> @brief   Krylov solvers
!> @details Krylov solvers
!>         
!-----------------------------------------------------------------------

module def_direct_solvers
  
  use def_kintyp_basic,   only : ip,rp
  use def_solvers,        only : solver
  implicit none

  private

  !----------------------------------------------------------------------
  !
  ! Direct solver class
  !
  !----------------------------------------------------------------------

  type, extends(solver), abstract :: direct_solver
   contains
  end type direct_solver

  !----------------------------------------------------------------------
  !
  ! Factorization-based solvers
  !
  ! Synmbolical factorization in done in setup
  ! Numerical factorization should be explicitly done before solving
  !
  !----------------------------------------------------------------------

  type, extends(direct_solver), abstract :: direct_factorization
   contains
     procedure (solver_symbolical),        pass, deferred :: symbolical     
     procedure (solver_numerical),         pass, deferred :: numerical     
     procedure (solver_cleaning),          pass, deferred :: cleaning     
     procedure (solver_partial_cleaning),  pass, deferred :: partial_cleaning 
     procedure (solver_set),               pass, deferred :: set
  end type direct_factorization
  
  abstract interface
          
     subroutine solver_symbolical(self,a)
       import :: direct_factorization
       !import :: mat
       class(direct_factorization), intent(inout) :: self
       class(*),                  intent(in)    :: a    !< Matrix A
     end subroutine solver_symbolical
  
     subroutine solver_numerical(self,a)
       import :: direct_factorization
       !import :: mat
       class(direct_factorization), intent(inout) :: self
       class(*),                    intent(in)    :: a    !< Matrix A
     end subroutine solver_numerical
  
     subroutine solver_set(self,a)
       import :: direct_factorization
       !import :: mat
       class(direct_factorization), intent(inout) :: self
       class(*),                    intent(in)    :: a    !< Matrix A
     end subroutine solver_set
  
     subroutine solver_cleaning(self)
       import :: direct_factorization
       class(direct_factorization), intent(inout) :: self
     end subroutine solver_cleaning
  
     subroutine solver_partial_cleaning(self)
       import :: direct_factorization
       class(direct_factorization), intent(inout) :: self
     end subroutine solver_partial_cleaning
  
  end interface
  
  public :: direct_solver 
  public :: direct_factorization
  
end module def_direct_solvers
