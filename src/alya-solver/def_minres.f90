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

module def_minres

  use def_kintyp,            only : ip,rp
  use def_mat,               only : mat
  use def_mat_dia,           only : mat_dia
  use def_mat_csr,           only : mat_csr
  use def_iterative_solvers, only : krylov_symmetric
  use mod_memory_basic,      only : memory_alloca
  use mod_memory_basic,      only : memory_deallo
  implicit none
  private

  type, extends(krylov_symmetric) :: minres
   contains
     procedure,         pass :: init 
     procedure,         pass :: alloca
     procedure,         pass :: deallo
     procedure,         pass :: solve
     procedure,         pass :: setup
  end type minres

  character(6), parameter :: vacal = 'minres'

  public :: minres

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-14
  !> @brief   Init
  !> @details Init
  !> 
  !-----------------------------------------------------------------------

  subroutine init(self)

    class(minres), intent(inout) :: self

  end subroutine init

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-14
  !> @brief   Deallocate
  !> @details Deallocate
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo(self)

    class(minres), intent(inout) :: self

  end subroutine deallo

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-14
  !> @brief   Allocate
  !> @details Allocate
  !> 
  !-----------------------------------------------------------------------

  subroutine alloca(self)

    class(minres),                intent(inout) :: self          !< Solver

  end subroutine alloca

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-14
  !> @brief   Setup
  !> @details Setup
  !> 
  !-----------------------------------------------------------------------

  subroutine setup(self,a)

    class(minres),                      intent(inout) :: self          !< Solver
    class(mat),                         intent(inout) :: a             !< Matrix A

  end subroutine setup

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-14
  !> @brief   Solve
  !> @details Solve
  !>          Creadits to https://web.stanford.edu/group/SOL/software/minres/
  !> 
  !-----------------------------------------------------------------------

  subroutine solve(self,b,x,a)
    class(minres),            intent(inout) :: self              !< Solver
    real(rp),        pointer, intent(in)    :: b(:)              !< RHS b
    real(rp),        pointer, intent(inout) :: x(:)              !< Solve x
    class(mat),               intent(in)    :: a                 !< Matrix A
    integer(ip)                             :: n,nn
    integer(ip)                             :: ndof
    integer(ip)                             :: ierro,i
    real(rp)                                :: resid

    real(rp)                                :: Anorm,gmin,gmax
    real(rp)                                :: tnorm2,ynorm2,denom,delta
    real(rp)                                :: ynorm,beta,oldb,qrnorm
    real(rp)                                :: rnorm,gamma,cs,sn,phi,phibar
    real(rp)                                :: gbar,dbar,rhs1,rhs2
    real(rp)                                :: oldeps,Acond,alfa,beta1,root,epsln
    real(rp)                                :: s,z,Arnorm,rnorml
    real(rp),        pointer                :: r1(:)
    real(rp),        pointer                :: r2(:)
    real(rp),        pointer                :: v(:)
    real(rp),        pointer                :: w(:)
    real(rp),        pointer                :: w1(:)
    real(rp),        pointer                :: w2(:)
    real(rp),        pointer                ::  y(:)

    nullify(r1,r2,v,w,w1,w2,y)
    ndof =  a % ndof1
    n    =  a % nrows
    nn   =  ndof * n

    call memory_alloca(self % memor,'R1',vacal,r1,nn)
    call memory_alloca(self % memor,'R2',vacal,r2,nn)
    call memory_alloca(self % memor,'V' ,vacal, v,nn)
    call memory_alloca(self % memor,'W2',vacal,w2,nn)
    call memory_alloca(self % memor,'W1',vacal,w1,nn)
    call memory_alloca(self % memor,'W' ,vacal, w,nn)
    call memory_alloca(self % memor,'Y ',vacal, y,nn)

    call self % start (b,x,a,r1,y,ierro)    
    if( ierro /= 0 ) goto 1

    beta1 = self % parallel_dot(r1,y)
    if( beta1 <= 0.0_rp ) goto 1
    beta1 = sqrt(beta1)
    resid = beta1
    call self % init_iterations(resid) ! Residuals and tolerance

    oldb   =  0.0_rp
    beta   =  beta1
    dbar   =  0.0_rp
    epsln  =  0.0_rp
    qrnorm =  beta1
    phibar =  beta1
    rhs1   =  beta1
    rhs2   =  0.0_rp
    tnorm2 =  0.0_rp
    ynorm2 =  0.0_rp
    cs     = -1.0_rp
    sn     =  0.0_rp

    do i = 1,nn
       r2(i) = r1(i)
    end do
    
    do while( self % output % iters < self % input % miter .and. resid > self % output % toler )
       !
       ! Obtain first Lanczos
       !
       s = 1.0_rp / beta       
       do i = 1,nn
          v(i) = y(i) * s
       end do
       !
       ! Compute v = A z and delta = (Az,z)
       !
       call self % parallel_mv(a,v,y)

       if( self % output % iters >= 1 ) then
          do i = 1,nn
             y(i) = y(i) - (beta/oldb) * r1(i)
          end do
       end if
       alfa = self % parallel_dot(v,y)
       !
       ! y^{k+1} = y^k - (alfa/beta) * r1^k - (gama^i/gama^{i-1}) * v^{i-1}
       !
       do i = 1,nn
          y(i)  = y(i) - (alfa/beta) * r2(i)
          r1(i) = r2(i)
          r2(i) = y(i)
       end do
       !
       !  L y^{k+1} = r2^{k+1} 
       !
       call self % preconditioning(r2,y)
       !
       ! gama^{i+1} = (z^{i+1},v^{i+1})^{1/2}  
       !
       oldb = beta
       beta = self % parallel_dot(r2,y)
       if( beta < 0.0_rp ) then
          stop
          ierro = 6
          go to 1
       end if
       beta = sqrt(beta)
       !
       ! Update constants
       !
       tnorm2 = alfa**2 + oldb**2 + beta**2
       if( self % output % iters == 0 ) then
          gmax = abs(alfa)
          gmin = gmax
       end if

       oldeps = epsln
       delta  = cs * dbar + sn * alfa       ! delta1 = 0         deltak
       gbar   = sn * dbar - cs * alfa       ! gbar 1 = alfa1     gbar k
       epsln  = sn * beta                   ! epsln2 = 0         epslnk+1
       dbar   =           - cs * beta       ! dbar 2 = beta2     dbar k+1

       gamma  = sqrt( gbar**2 + beta**2 )   ! gammak
       cs     = gbar / gamma                ! ck
       sn     = beta / gamma                ! sk
       phi    = cs * phibar                 ! phik
       phibar = sn * phibar                 ! phibark+1
       denom  = 1.0_rp/gamma

       do i = 1, n
          w1(i) = w2(i)
          w2(i) = w(i)
          w(i)  = ( v(i) - oldeps*w1(i) - delta*w2(i) ) * denom
          x(i)  =   x(i) + phi * w(i)
       end do
       !
       ! Go round again.
       !
       gmax   = max( gmax, gamma )
       gmin   = min( gmin, gamma )
       z      = rhs1 / gamma
       ynorm2 = z**2  +  ynorm2
       rhs1   = rhs2  -  delta * z
       rhs2   =       -  epsln * z
       !
       ! Estimate various norms and test for convergence.
       !
       Anorm  = sqrt( tnorm2 )
       ynorm  = sqrt( ynorm2 )
       qrnorm = phibar
       rnorm  = qrnorm                     ! ||L^-1 r||
       rnorml = rnorm
       root   = sqrt( gbar**2 + dbar**2)
       Arnorm = rnorml * root
       !
       ! Estimate  cond(A).
       ! In this version we look at the diagonals of  R  in the
       ! factorization of the lower Hessenberg matrix,  Q * H = R,
       ! where H is the tridiagonal matrix from Lanczos with one
       ! extra row, beta(k+1) e_k^T.
       !
       Acond  = gmax / gmin
       !
       ! Update and output residuals
       !
       resid                       = rnorm
       self % output % resip_old   = self % output % resip_final
       self % output % resip_final = resid * self % output % bnorp_inv
       self % output % iters       = self % output % iters + 1
       call self % outcvg(b,x,a)

    end do

1   continue

    call self % end(b,x,a)
    call memory_deallo(self % memor,'R1',vacal,r1)
    call memory_deallo(self % memor,'R2',vacal,r2)
    call memory_deallo(self % memor,'V' ,vacal, v)
    call memory_deallo(self % memor,'W2',vacal,w2)
    call memory_deallo(self % memor,'W1',vacal,w1)
    call memory_deallo(self % memor,'W' ,vacal, w)
    call memory_deallo(self % memor,'Y ',vacal, y)

    if( self % input % lun_cvgso /= 0 .and. self % iwrite .and. ierro > 0 ) then
       if( ierro == 6 ) write(self % input % lun_cvgso,2) 
    end if

2   format('# Preconditioner is not positive definite')

  end subroutine solve

end module def_minres
!> @}
