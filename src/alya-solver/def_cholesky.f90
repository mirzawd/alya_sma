!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




!-----------------------------------------------------------------------
!> @addtogroup Solvers
!> @{
!> @file    def_lu.f90
!> @author  houzeaux
!> @date    2022-03-01
!> @brief   LU solvers
!> @details LU solvers
!>          Initilization => setup (reordering) => symbolical => numerical
!>          => partial cleaning or cleaning
!>         
!-----------------------------------------------------------------------

module def_cholesky

  use def_kintyp,             only : ip,rp
  use def_mat,                only : mat
  use def_mat_sky,            only : mat_sky
  use def_mat_csr,            only : mat_csr
  use mod_memory,             only : memory_alloca
  use mod_memory,             only : memory_deallo
  use def_direct_solvers,     only : direct_factorization
  use mod_alya_direct_solver, only : alya_cholesky_solver
  implicit none
  private
  
  type, extends(direct_factorization) :: cholesky
     type(alya_cholesky_solver) :: dir
   contains
     procedure,           pass :: init     
     procedure,           pass :: alloca   
     procedure,           pass :: deallo   
     procedure,           pass :: solve    
     procedure,           pass :: setup
     procedure,           pass :: set
     procedure,           pass :: numerical 
     procedure,           pass :: symbolical 
     procedure,           pass :: cleaning     
     procedure,           pass :: partial_cleaning     
  end type cholesky

  character(8), parameter :: vacal = 'cholesky'

  public :: cholesky
    
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Initialization
  !> @details Initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine init(self)
    
    class(cholesky), intent(inout) :: self

    call self % dir % init()

  end subroutine init
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Partial cleaning
  !> @details Partial cleaning
  !> 
  !-----------------------------------------------------------------------

  subroutine partial_cleaning(self)
    
    class(cholesky), intent(inout) :: self

    call self % dir % deallo()

  end subroutine partial_cleaning
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Partial cleaning
  !> @details Partial cleaning
  !> 
  !-----------------------------------------------------------------------

  subroutine cleaning(self)
    
    class(cholesky), intent(inout) :: self

    call self % dir % deallo()
    call self % deallo()
 
  end subroutine cleaning
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Deallocate
  !> @details Deallocate
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo(self)
    
    class(cholesky), intent(inout) :: self

  end subroutine deallo

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Symbolical factorization
  !> @details Symbolical factorization. Nothing to do!
  !> 
  !-----------------------------------------------------------------------

  subroutine symbolical(self,a)
    
    class(cholesky), intent(inout) :: self
    class(*),        intent(in)    :: a    !< Matrix A

  end subroutine symbolical
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Numerical factorization
  !> @details Numerical factorization
  !> 
  !-----------------------------------------------------------------------

  subroutine numerical(self,a)
    
    class(cholesky), intent(inout) :: self
    class(*),        intent(in)    :: a    !< Matrix A

    select type ( a )
    class is ( mat_sky ) ; call self % dir % factor_numerical(a)
    class default        ; call runend('DEF_CHOLESKY: WRONG MATRIX FORMAT')
    end select
    
  end subroutine numerical
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Allocate
  !> @details Allocate
  !> 
  !-----------------------------------------------------------------------

  subroutine alloca(self)
    
    class(cholesky), intent(inout) :: self
        
  end subroutine alloca

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Set
  !> @details Set
  !> 
  !-----------------------------------------------------------------------

  subroutine set(self,a)

    class(cholesky),                    intent(inout) :: self          !< Solver
    class(*),                           intent(in)    :: a             !< Matrix A
    
  end subroutine set
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Setup
  !> @details Setup and reorder to minimize fillin
  !> 
  !-----------------------------------------------------------------------

  subroutine setup(self,a)

    class(cholesky),                    intent(inout) :: self          !< Solver
    class(mat),                         intent(inout) :: a             !< Matrix A

    call self % numerical(a)
    
  end subroutine setup
  
  !----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Solve
  !> @details Factorize, solve and partial clean
  !> 
  !-----------------------------------------------------------------------

  subroutine solve(self,b,x,a)

    class(cholesky),          intent(inout) :: self          !< Solver
    real(rp),        pointer, intent(in)    :: b(:)          !< RHS b
    real(rp),        pointer, intent(inout) :: x(:)          !< Solve x
    class(mat),               intent(in)    :: a             !< Matrix A

    call self % dir % solve(b,x,a)

  end subroutine solve

end module def_cholesky
!> @}
