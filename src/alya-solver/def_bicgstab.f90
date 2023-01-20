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

module def_bicgstab

  use def_kintyp,            only : ip,rp
  use def_mat,               only : mat
  use def_mat_dia,           only : mat_dia
  use def_mat_csr,           only : mat_csr
  use def_iterative_solvers, only : krylov_unsymmetric
  use mod_memory_basic,      only : memory_alloca
  use mod_memory_basic,      only : memory_deallo
  implicit none
  private
  
  type, extends(krylov_unsymmetric) :: bicgstab
   contains
     procedure,         pass :: init 
     procedure,         pass :: alloca
     procedure,         pass :: deallo
     procedure,         pass :: solve
     procedure,         pass :: setup
  end type bicgstab

  character(8), parameter :: vacal = 'bicgstab'

  public :: bicgstab
  
contains

  subroutine init(self)
    
    class(bicgstab), intent(inout) :: self

  end subroutine init  
  
  subroutine deallo(self)
    
    class(bicgstab), intent(inout) :: self
        
  end subroutine deallo
  
  subroutine alloca(self)

    class(bicgstab),          intent(inout) :: self          !< Solver
    
  end subroutine alloca
 
  subroutine setup(self,a)

    class(bicgstab),                    intent(inout) :: self          !< Solver
    class(mat),                         intent(inout) :: a             !< Matrix A

  end subroutine setup
 
  subroutine solve(self,b,x,a)
    class(bicgstab),          intent(inout) :: self              !< Solver
    real(rp),        pointer, intent(in)    :: b(:)              !< RHS b
    real(rp),        pointer, intent(inout) :: x(:)              !< Solve x
    class(mat),               intent(in)    :: a                 !< Matrix A
    integer(ip)                             :: n,nn
    integer(ip)                             :: ndof
    integer(ip)                             :: ierro,i
    real(rp)                                :: alpha,beta,rho
    real(rp)                                :: newrho,resid
    real(rp)                                :: denom,omega,dot(2)

    real(rp),        pointer                :: r0(:)
    real(rp),        pointer                :: r(:)
    real(rp),        pointer                :: p(:)
    real(rp),        pointer                :: q(:)
    real(rp),        pointer                :: s(:)
    real(rp),        pointer                :: t(:)
    real(rp),        pointer                :: w(:)

    nullify(r0,r,p,q,s,t,w)
    ndof =  a % ndof1
    n    =  a % nrows
    nn   =  ndof * n
    
    call memory_alloca(self % memor,'R0',vacal,r0,nn)
    call memory_alloca(self % memor,'R' ,vacal,r ,nn)
    call memory_alloca(self % memor,'P' ,vacal,p ,nn)
    call memory_alloca(self % memor,'Q' ,vacal,q ,nn)
    call memory_alloca(self % memor,'S' ,vacal,s ,nn)
    call memory_alloca(self % memor,'T' ,vacal,t ,nn)
    call memory_alloca(self % memor,'W' ,vacal,w ,nn)

    call self % start (b,x,a,r,r0,ierro)    
    if( ierro /= 0 ) goto 1
    
    newrho = self % parallel_dot(r0,r0)
    if( newrho <= 0.0_rp ) goto 1
    
    resid = sqrt(newrho)    
    alpha = 0.0_rp
    omega = 1.0_rp
    rho   = 1.0_rp
    do i = 1,nn
       r(i) = r0(i)
    end do
    call self % init_iterations(resid)                           ! Residuals and tolerance
    
    do while( self % output % iters < self % input % miter .and. resid > self % output % toler )
       !
       ! beta = (rho^k/rho^{k-1})*(alpha/omega)
       !
       if( rho /= 0.0_rp .and. omega /= 0.0_rp ) then
          beta = (newrho/rho) * (alpha/omega)
       else
          ierro = 1
          goto 1
       end if
       !
       ! p^{k+1} = r^{k} + beta*(p^k - omega*q^k)
       !
       do i = 1,nn
          p(i) = r(i) + beta * (p(i) - omega * q(i))
       end do
       !
       ! L q^{k+1} = A ( R^{-1} p^{k+1} )
       !
       call self % preconditioning(p,q,w,a)        
       !
       ! alpha = rho^k / <r0,q^{k+1}>
       !
       alpha = self % parallel_dot(r0,q)
       if( alpha == 0.0_rp ) then
          ierro = 2
          goto 1
       end if
       alpha = newrho / alpha
       !
       ! s = r^k - alpha*q^{k+1}
       !
       do i = 1,nn
          s(i) = r(i) - alpha * q(i)
       end do
       !
       ! L t = A ( R^{-1} s )
       !
       call self % preconditioning(s,t,w,a)
       !
       ! omega = <t,s> / <t,t>
       !
       dot   = self % parallel_dots(t,s,t)
       omega = dot(1)
       denom = dot(2)
       if( denom == 0.0_rp ) then
          ierro = 2
          goto 1
       end if
       omega = omega / denom
       !
       ! x^{k+1} = x^k + alpha*p^{k+1} + omega*s
       !
       do i = 1,nn
          x(i) = x(i) + alpha * p(i) + omega * s(i)
       end do
       !
       ! r^{k+1} = s - omega*t
       !
       do i = 1,nn
          r(i) = s(i) - omega * t(i)
       end do

       !
       ! rho^k = <r0,r^k> and || r^{k+1} ||
       !
       rho    = newrho
       dot    = self % parallel_dots(r,r,r0)
       resid  = dot(1)
       newrho = dot(2)
       !
       ! Update and output residuals
       !
       resid                       = sqrt(resid)
       self % output % resip_old   = self % output % resip_final
       self % output % resip_final = resid * self % output % bnorp_inv
       self % output % iters       = self % output % iters + 1
       call self % outcvg(b,x,a)
       
    end do

1   continue

    call self % end(b,x,a,r)

    if( self % input % lun_cvgso /= 0 .and. self % iwrite .and. ierro > 0 ) then
       if( ierro == 1 ) write(self % input % lun_cvgso,201) self % output % iters
       if( ierro == 2 ) write(self % input % lun_cvgso,202) self % output % iters
    end if

    call memory_deallo(self % memor,'W' ,vacal,w)
    call memory_deallo(self % memor,'T' ,vacal,t)
    call memory_deallo(self % memor,'S' ,vacal,s)
    call memory_deallo(self % memor,'R' ,vacal,r)
    call memory_deallo(self % memor,'P' ,vacal,p)
    call memory_deallo(self % memor,'Q' ,vacal,q)
    call memory_deallo(self % memor,'R0',vacal,r0)

100 format(i7,1x,e12.6)
110 format(i5,18(2x,e12.6))
201 format(&
       & '# Error at iteration ',i6,&
       & 'Dividing by zero: alpha = beta = (rho^k/rho^{k-1})*(alpha/omega)')
202 format(&
       & '# Error at iteration ',i6,&
       & 'Dividing by zero: alpha = beta = alpha = rho^k / <r0,q^{k+1}>')
    
  end subroutine solve

end module def_bicgstab
!> @}
