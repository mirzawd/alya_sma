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

module def_lu_factorization

  use def_kintyp,             only : ip,rp
  use def_mat,                only : mat
  use def_mat_csr,            only : mat_csr
  use mod_memory,             only : memory_alloca
  use mod_memory,             only : memory_deallo
  use def_direct_solvers,     only : direct_factorization
  use mod_alya_direct_solver, only : alya_direct_solver
  use mod_graphs_basic,       only : graphs_permut_metis_postordering
  implicit none
  private
  
  type, extends(direct_factorization) :: lu_factorization
     type(alya_direct_solver)  :: dir
     integer(ip),  pointer     :: ia(:)   ! Renumbered graph
     integer(ip),  pointer     :: ja(:)   ! Renumbered graph
     integer(ip),  pointer     :: perm(:) ! Permutation
     integer(ip),  pointer     :: invp(:) ! Inverse permutation
   contains
     procedure,           pass :: init     
     procedure,           pass :: alloca   
     procedure,           pass :: deallo   
     procedure,           pass :: solve    
     procedure,           pass :: set    
     procedure,           pass :: setup    
     procedure,           pass :: symbolical     
     procedure,           pass :: numerical     
     procedure,           pass :: reordering     
     procedure,           pass :: cleaning     
     procedure,           pass :: partial_cleaning     
  end type lu_factorization

  character(16), parameter :: vacal = 'lu_factorization'

  public :: lu_factorization
    
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
    
    class(lu_factorization), intent(inout) :: self

    call self % dir % init()
    nullify(self % ia)
    nullify(self % ja)
    nullify(self % perm)
    nullify(self % invp)

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
    
    class(lu_factorization), intent(inout) :: self

    call self % dir % deallo_numerical()

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
    
    class(lu_factorization), intent(inout) :: self

    call self % dir % deallo_numerical ()
    call self % dir % deallo_symbolical()
    call self % deallo                 ()
 
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
    
    class(lu_factorization), intent(inout) :: self
    
    call memory_deallo(self % memor,'SELF % IA'  ,vacal,self % ia)
    call memory_deallo(self % memor,'SELF % JA'  ,vacal,self % ja)
    call memory_deallo(self % memor,'SELF % PERM',vacal,self % perm)
    call memory_deallo(self % memor,'SELF % INVP',vacal,self % invp)

  end subroutine deallo

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Symbolical factorization
  !> @details Symbolical factorization
  !> 
  !-----------------------------------------------------------------------

  subroutine symbolical(self,a)
    
    class(lu_factorization), intent(inout) :: self
    class(*),                intent(in)    :: a    !< Matrix A

    select type ( a )
    class is ( mat ) 
       call self % dir % factor_symbolical(a,self % ia,self % ja)
    end select
    
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
    
    class(lu_factorization), intent(inout) :: self
    class(*),                intent(in)    :: a    !< Matrix A

    select type ( a )
    class is ( mat_csr )
       call self % dir % factor_numerical(a,self % invp,self % perm,self % ia,self % ja)
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
    
    class(lu_factorization), intent(inout) :: self
        
  end subroutine alloca
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Reordering
  !> @details Reordering
  !> 
  !-----------------------------------------------------------------------

  subroutine reordering(self,a)
    
    class(lu_factorization),            intent(inout) :: self          !< Solver
    class(mat),                         intent(in)    :: a             !< Matrix A
    integer(ip)                                       :: n,nz,i
    
    select type ( a )
    class is ( mat_csr )
       
       n  = a % nrows 
       nz = a % nz
       call memory_alloca(self % memor,'SELF % IA'  ,vacal,self % ia  ,n+1_ip)
       call memory_alloca(self % memor,'SELF % JA'  ,vacal,self % ja  ,nz)
       call memory_alloca(self % memor,'SELF % PERM',vacal,self % perm,n)
       call memory_alloca(self % memor,'SELF % INVP',vacal,self % invp,n)
       do i = 1,n+1
          self % ia(i) = a % ia(i)
       end do
       do i = 1,nz
          self % ja(i) = a % ja(i)
       end do
       
       call graphs_permut_metis_postordering(&
            a % nrows   , a % nz,      &
            self % ia   , self % ja,   &
            self % perm , self % invp, &
            memor=self % memor)

    end select

  end subroutine reordering
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Set
  !> @details Reorder to minimize fillin
  !> 
  !-----------------------------------------------------------------------

  subroutine set(self,a)

    class(lu_factorization),            intent(inout) :: self          !< Solver
    class(*),                           intent(in)    :: a             !< Matrix A

    select type ( a )
    class is ( mat ) 
       call self % reordering  (a)
       call self % symbolical  (a)
    end select
    
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

    class(lu_factorization),            intent(inout) :: self          !< Solver
    class(mat),                         intent(inout) :: a             !< Matrix A

    call self % numerical(a)
 
  end subroutine setup
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Solve
  !> @details Factorize, solve and partial clean
  !> 
  !-----------------------------------------------------------------------

  subroutine solve(self,b,x,a)

    class(lu_factorization),  intent(inout) :: self          !< Solver
    real(rp),        pointer, intent(in)    :: b(:)          !< RHS b
    real(rp),        pointer, intent(inout) :: x(:)          !< Solve x
    class(mat),               intent(in)    :: a             !< Matrix A

    call self % dir % solve(b,x,a,self % invp)
    
  end subroutine solve

end module def_lu_factorization
!> @}
