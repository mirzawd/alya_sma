!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




!-----------------------------------------------------------------------
!> @addtogroup Solvers
!> @{
!> @file    def_gmres.f90
!> @author  houzeaux
!> @date    2022-03-01
!> @brief   Conjugate gradient
!> @details Conjugate gradient
!>         
!-----------------------------------------------------------------------

module def_cg

  use def_kintyp,            only : ip,rp
  use def_mat,               only : mat
  use def_mat_dia,           only : mat_dia
  use def_mat_csr,           only : mat_csr
  use def_iterative_solvers, only : krylov_symmetric
  use mod_memory_basic,      only : memory_alloca
  use mod_memory_basic,      only : memory_deallo
  implicit none
  private
  
  type, extends(krylov_symmetric) :: cg
   contains
     procedure,         pass :: init 
     procedure,         pass :: alloca
     procedure,         pass :: deallo
     procedure,         pass :: solve
     procedure,         pass :: setup
  end type cg

  character(2), parameter :: vacal = 'cg'

  public :: cg
  
contains

  subroutine init(self)
    
    class(cg), intent(inout) :: self
    
  end subroutine init  
  
  subroutine deallo(self)
    
    class(cg), intent(inout) :: self
    
  end subroutine deallo
  
  subroutine alloca(self)

    class(cg),                intent(inout) :: self          !< Solver
    
  end subroutine alloca
 
  subroutine setup(self,a)

    class(cg),                          intent(inout) :: self          !< Solver
    class(mat),                         intent(inout) :: a             !< Matrix A

  end subroutine setup
 
  subroutine solve(self,b,x,a)
    class(cg),                intent(inout) :: self              !< Solver
    real(rp),        pointer, intent(in)    :: b(:)              !< RHS b
    real(rp),        pointer, intent(inout) :: x(:)              !< Solve x
    class(mat),               intent(in)    :: a                 !< Matrix A
    integer(ip)                             :: n,nn
    integer(ip)                             :: ndof
    integer(ip)                             :: ierro,i
    real(rp)                                :: alpha,beta,rho
    real(rp)                                :: newrho,resid
    real(rp)                                :: denom

    real(rp),        pointer                :: r(:)
    real(rp),        pointer                :: p(:)
    real(rp),        pointer                :: q(:)

    nullify(r,p,q)
    ndof =  a % ndof1
    n    =  a % nrows
    nn   =  ndof * n
    
    call memory_alloca(self % memor,'R',vacal,r,nn)
    call memory_alloca(self % memor,'P',vacal,p,nn)
    call memory_alloca(self % memor,'Q',vacal,q,nn)

    call self % start (b,x,a,r,q,ierro)    
    if( ierro /= 0 ) goto 1
    
    newrho = self % parallel_dot(r,q)
    if( newrho <= 0.0_rp ) goto 1

    resid = sqrt(newrho)    
    call self % init_iterations(resid)                           ! Residuals and tolerance
    do i = 1,nn
       p(i) = q(i)
    end do
    
    do while( self % output % iters < self % input % miter .and. resid > self % output % toler )
       !
       ! q^{k+1} =  A p^{k+1}
       !
       call self % parallel_mv(a,p,q)
       !
       ! alpha = rho^k / <p^{k+1},q^{k+1}>
       !       
       denom = self % parallel_dot(p,q)
       if( denom <= 0.0_rp ) then
          ierro = 2
          goto 1
       end if
       rho   = newrho
       alpha = newrho / denom
       !
       ! x^{k+1} = x^k + alpha * p^{k+1}
       ! r^{k+1} = r^k - alpha * q^{k+1}
       !
       do i = 1,nn
          x(i) = x(i) + alpha * p(i)
          r(i) = r(i) - alpha * q(i)
       end do
       !
       !  L z^{k+1} = r^{k+1}
       !
       call self % preconditioning(r,q)
       newrho = self % parallel_dot(r,q)
       !
       ! beta  = rho^k / rho^{k-1}
       !
       beta = newrho / rho
       !
       ! p^{k+1} = z^k + beta*p^k
       !
       do i = 1,nn
          p(i) = q(i) + beta * p(i)
       end do
       !
       ! Update and output residuals
       !
       resid                       = sqrt(newrho)
       self % output % resip_old   = self % output % resip_final
       self % output % resip_final = resid * self % output % bnorp_inv
       self % output % iters       = self % output % iters + 1
       call self % outcvg(b,x,a,r)
       
    end do

1   continue

    call self % end(b,x,a,r)
    call memory_deallo(self % memor,'R',vacal,r)
    call memory_deallo(self % memor,'P',vacal,p)
    call memory_deallo(self % memor,'Q',vacal,q)

    if( self % input % lun_cvgso /= 0 .and. self % iwrite .and. ierro > 0 ) then
       if( ierro == 2 ) write(self % input % lun_cvgso,2) self % output % iters
    end if
    
2   format(&
       & '# Error at iteration ',i6,&
       & 'Dividing by zero: alpha = rho^k / <p^{k+1},q^{k+1}>')
    
  end subroutine solve

end module def_cg
!> @}
