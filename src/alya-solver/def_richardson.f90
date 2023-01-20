!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




!-----------------------------------------------------------------------
!> @addtogroup Solvers
!> @{
!> @file    def_krylov_solvers.f90
!> @author  houzeaux
!> @date    2022-03-01
!> @brief   Krylov solvers
!> @details Krylov solvers
!>         
!-----------------------------------------------------------------------

module def_richardson

  use def_kintyp,            only : ip,rp
  use def_mat,               only : mat
  use def_mat_dia,           only : mat_dia
  use def_mat_csr,           only : mat_csr
  use def_iterative_solvers, only : stationary
  use mod_memory_basic,      only : memory_alloca
  use mod_memory_basic,      only : memory_deallo
  implicit none
  
  type, extends(stationary)       :: richardson
   contains
     procedure,               pass :: init     
     procedure,               pass :: solve    
     procedure,               pass :: setup    
     procedure,               pass :: deallo   
     procedure,               pass :: alloca     
     procedure,               pass :: mv
  end type richardson
  
  public :: richardson
  
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
    
    class(richardson), intent(inout) :: self
    
  end subroutine init  

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Deallocate
  !> @details Deallocate
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo(self)
    
    class(richardson), intent(inout) :: self

  end subroutine deallo

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Allocate
  !> @details Allocate
  !> 
  !-----------------------------------------------------------------------

  subroutine alloca(self)
    
    class(richardson),            intent(inout) :: self
    
  end subroutine alloca

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Setup
  !> @details Setup
  !> 
  !-----------------------------------------------------------------------

  subroutine setup(self,a)

    class(richardson),                  intent(inout) :: self          !< Solver
    class(mat),                         intent(inout) :: a             !< Matrix A

  end subroutine setup
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Solve
  !> @details Solve system
  !> 
  !-----------------------------------------------------------------------

  subroutine solve(self,b,x,a)

    class(richardson),        intent(inout) :: self          !< Solver
    real(rp),        pointer, intent(in)    :: b(:)          !< RHS b
    real(rp),        pointer, intent(inout) :: x(:)          !< Solve x
    class(mat),               intent(in)    :: a             !< Matrix A

    integer(ip)                             :: i,nn,ierro
    real(rp)                                :: relax
    real(rp)                                :: resid
    real(rp),        pointer                :: r(:)
    real(rp),        pointer                :: z(:)

    nullify(r,z)
    
    nn = self % nn

    call memory_alloca(self % memor,'R','solve',r,nn)
    call memory_alloca(self % memor,'Z','solve',z,nn)

    call self % start (b,x,a,r,z,ierro)    

    if( ierro /= 0 ) goto 1   
    
    resid = sqrt(self % parallel_dot(r,z))
    relax = self % input % relax
    call self % init_iterations(resid)                           ! Residuals and tolerance

    do while ( self % output % iters < self % input % miter .and. resid > self % output % toler )

       do i = 1,nn
          x(i) = x(i) + relax * z(i) 
       end do
       call self % parallel_residual(a,x,b,r)
       call self % preconditioning  (r,z)
      
       resid                       = sqrt(self % parallel_dot(r,z))
       self % output % resip_old   = self % output % resip_final
       self % output % resip_final = resid * self % output % bnorp_inv
       self % output % iters       = self % output % iters + 1
       call self % outcvg(b,x,a,r)

    end do 

1   continue

    call self % end(b,x,a,r)

    call memory_deallo(self % memor,'R','solve',r)
    call memory_deallo(self % memor,'Z','solve',z)

  end subroutine solve

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-08
  !> @brief   Multiply by the preconditioner
  !> @details Multiply y = D x
  !> 
  !-----------------------------------------------------------------------
  
  subroutine mv(self,x,y,a)

    class(richardson),        intent(in)    :: self          !< Solver
    real(rp),        pointer, intent(in)    :: x(:)          !< RHS b
    real(rp),        pointer, intent(inout) :: y(:)          !< Solve x
    class(mat),               intent(in)    :: a             !< Matrix A

  end subroutine mv
  
end module def_richardson
!> @}
