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

module def_gauss_elimination

  use def_kintyp,         only : ip,rp
  use def_mat,            only : mat
  use def_mat_dia,        only : mat_dia
  use def_mat_den,        only : mat_den
  use def_direct_solvers, only : direct_solver
  implicit none
  private
  
  type, extends(direct_solver) :: gauss_elimination
   contains
     procedure, pass :: init     
     procedure, pass :: alloca   
     procedure, pass :: deallo   
     procedure, pass :: solve    
     procedure, pass :: setup    
  end type gauss_elimination

  public :: gauss_elimination
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Diag solver
  !> @details Diag solver Ax=b
  !> 
  !-----------------------------------------------------------------------

  subroutine init(self)
    
    class(gauss_elimination), intent(inout) :: self
    
  end subroutine init
  
  subroutine deallo(self)
    
    class(gauss_elimination), intent(inout) :: self
    
  end subroutine deallo
  
  subroutine alloca(self)
    
    class(gauss_elimination), intent(inout) :: self
        
  end subroutine alloca
  
  subroutine setup(self,a)

    class(gauss_elimination),           intent(inout) :: self          !< Solver
    class(mat),                         intent(inout) :: a             !< Matrix A

    self % n    = a % nrows
    self % ndof = a % ndof1
    self % nn   = self % n * self % ndof

  end subroutine setup
  
  subroutine solve(self,b,x,a)

    class(gauss_elimination), intent(inout) :: self          !< Solver
    real(rp),        pointer, intent(in)    :: b(:)          !< RHS b
    real(rp),        pointer, intent(inout) :: x(:)          !< Solve x
    class(mat),               intent(in)    :: a             !< Matrix A

    if( associated(b) .and. associated(x) ) then
       
       select type ( a )
          
       class is ( mat_dia ) ; call a % solve(b,x)
       class is ( mat_den ) ; call a % solve(b,x)
          
       end select

    end if
    
  end subroutine solve

end module def_gauss_elimination
!> @}
