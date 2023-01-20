!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




!-----------------------------------------------------------------------
!> @addtogroup Solvers
!> @{
!> @file    def_sor.f90
!> @author  houzeaux
!> @date    2022-03-01
!> @brief   General SOR
!> @details SSOR solvers with relaxation parameter w. This is a
!>          non-symmetric preconditioner.
!>          Solve x^k+1 = x^k + M^-1 (b-A x^k)
!>          M^-1 = (D+wL)^-1
!>          If w=1, this is the Gauss-Sediel method:
!>          M^-1 = (D+L)^-1
!>         
!-----------------------------------------------------------------------

module def_sor

  use def_kintyp,            only : ip,rp
  use def_mat,               only : mat
  use def_mat_dia,           only : mat_dia
  use def_mat_csr,           only : mat_csr
  use def_iterative_solvers, only : stationary
  use mod_memory_basic,      only : memory_alloca
  use mod_memory_basic,      only : memory_deallo
  use mod_maths_basic,       only : maths_inverse
  implicit none
  
  type, extends(stationary)       :: sor
     type(mat_dia)                :: dia
  contains
    procedure,               pass :: init     
    procedure,               pass :: solve    
    procedure,               pass :: setup    
    procedure,               pass :: deallo   
    procedure,               pass :: alloca   
    procedure,               pass :: iteration   
    procedure,               pass :: mv   
  end type sor                
 
  public :: sor
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-02
  !> @brief   Sor solver
  !> @details Sor solver
  !> 
  !-----------------------------------------------------------------------

  subroutine init(self)
    
    class(sor), intent(inout) :: self

    call self % dia % init()

  end subroutine init  

  subroutine deallo(self)
    
    class(sor), intent(inout) :: self

    call self % dia % deallo()
    
  end subroutine deallo

  subroutine alloca(self)
    
    class(sor), intent(inout) :: self
    
  end subroutine alloca

  subroutine setup(self,a)

    class(sor),                         intent(inout) :: self          !< Solver
    class(mat),                         intent(inout) :: a             !< Matrix A
    !
    ! Diagonal matrix
    !
    call a    % diag             (self % dia)
    call self % parallel_exchange(self % dia)
    call self % dia % inverse    ()

  end subroutine setup

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-09
  !> @brief   Gauss-Seidel method
  !> @details A = L + D + U
  !>          (D+wL) x^k+1 = w b - [wU +(w-1)D] x^k
  !>                 x^k+1 = x^k + (D+wL) (b - A x^k)
  !> 
  !-----------------------------------------------------------------------
  
  subroutine solve(self,b,x,a)

    class(sor),               intent(inout) :: self          !< Solver
    real(rp),        pointer, intent(in)    :: b(:)          !< RHS b
    real(rp),        pointer, intent(inout) :: x(:)          !< Solve x
    class(mat),               intent(in)    :: a             !< Matrix A

    integer(ip)                             :: i,nn
    integer(ip)                             :: ierro
    real(rp)                                :: relax
    real(rp)                                :: resid
    real(rp),        pointer                :: r(:)
    real(rp),        pointer                :: z(:)

    nullify(r,z)
    nn    = a % nrows * a % ndof1
    relax = self % input % relax
    call memory_alloca(self % memor,'R','solve',r,nn)

    if( self % input % miter == -1 ) then
       !
       ! Gauss-Seidel as a preconditioner
       !
       do i = 1,nn
          x(i) = 0.0_rp
       end do
       call self % iteration(b,x,a,r)
          
    else
       !
       ! Gauss-Seidel as a solver
       !
       call memory_alloca(self % memor,'Z','solve',z,nn)
       call self % start (b,x,a,r,z,ierro)    

       if( ierro /= 0 ) goto 1   

       resid                       = self % parallel_L2norm(z)
       self % output % resip_init  = resid * self % output % bnorp_inv
       self % output % resip_old   = self % output % resip_init 
       self % output % resip_final = self % output % resip_init 
       self % output % toler       = self % tolerance()
 
       do while ( self % output % iters < self % input % miter .and. resid > self % output % toler )

          call self % iteration(b,x,a,r)

          resid                       = self % parallel_L2norm(r)
          self % output % resip_old   = self % output % resip_final
          self % output % resip_final = resid * self % output % bnorp_inv
          self % output % iters       = self % output % iters + 1
          call self % outcvg(b,x,a,r)

       end do

1      continue

       call memory_deallo(self % memor,'Z','solve',z)
       call self % end(b,x,a)

    end if
    
    call memory_deallo(self % memor,'R','solve',r)
    
  end subroutine solve

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-08
  !> @brief   One iteration
  !> @details One single iteration of Gauss-Seidel
  !>                                1  +-                                                  -+
  !>          xi^k+1 = (1-w)xi^k + --- | bi - sum_j=1^i-1 aij xj^k+1 - sum_j=i+1^n aij xj^k |
  !>                               aii +-                                                  -+
  !>          Add: aii xi^k - aii xi^k in parenthesis
  !>
  !>                           w  +-                                                -+
  !>          xi^k+1 = xi^k + --- | bi - sum_j=1^i-1 aij xj^k+1 - sum_j=i^n aij xj^k |
  !>                          aii +-                                                -+
  !>
  !-----------------------------------------------------------------------
  
  subroutine iteration(self,b,x,a,r)

    class(sor),               intent(inout) :: self          !< Solver
    real(rp),        pointer, intent(in)    :: b(:)          !< RHS b
    real(rp),        pointer, intent(inout) :: x(:)          !< Solve x
    class(mat),               intent(in)    :: a             !< Matrix A
    real(rp),        pointer                :: r(:)
    integer(ip)                             :: i,idof,k
    !
    ! Jacobi on interface
    !
    call self % parallel_residual(a,x,b,r,ON_INTERFACE=.true.)
    i = self % nint * self % ndof
    do k = self % nint+1,self % n
       do idof = 1,self % ndof
          i = i + 1
          call self % dia % scale(r,i,idof)
          x(i) = x(i) + r(i)
       end do
    end do
    !
    ! Gauss-Seidel on other nodes
    !
    i = 0
    do k = 1,self % nint
       do idof = 1,self % ndof
          i    = i + 1
          r(i) = a % mv_row(x,i,idof)
          r(i) = b(i) - r(i)
          call self % dia % scale(r,i,idof)
          x(i) = x(i) + self % input % relax * r(i)
       end do
    end do
    
  end subroutine iteration
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-08
  !> @brief   Multiply by the preconditioner
  !> @details Multiply y = (L+D) x
  !> 
  !-----------------------------------------------------------------------
  
  subroutine mv(self,x,y,a)

    class(sor),               intent(in)    :: self          !< Solver
    real(rp),        pointer, intent(in)    :: x(:)          !< RHS b
    real(rp),        pointer, intent(inout) :: y(:)          !< Solve x
    class(mat),               intent(in)    :: a             !< Matrix A
    integer(ip)                             :: i,idof,k    
    !
    ! Jacobi on interface
    !
    i = 0
    do k = self % nint+1,self % n
       do idof = 1,self % ndof
          i = i + 1
          y(i) = x(i) * maths_inverse(self % dia % vA(idof,k),0.0_rp)
       end do
    end do
    !
    ! Gauss-Seidel on other nodes
    !
    call a % mv_lower(x,y,1_ip,self % nint)
    
  end subroutine mv
  
end module def_sor
!> @}
