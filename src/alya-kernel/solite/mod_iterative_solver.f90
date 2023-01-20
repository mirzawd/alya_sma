!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @defgroup   Sequential_Iterative_Solvers
!> Toolbox for sequential iterative solvers
!> @ingroup    Algebraic_Solver
!> @{
!> @file    mod_iterative_solver.f90
!> @author  houzeaux
!> @date    2018-08-10
!> @brief   Iterative solvers
!> @details Krylov solvers
!-----------------------------------------------------------------------

module mod_iterative_solver
  
  use def_kintyp, only : ip,rp,lg
  use def_solver, only : memit,soltyp
  use mod_solver, only : solver_diagonal
  use mod_matrix, only : matrix_CSR_SpMV
  use mod_matrix, only : matrix_diagonal_CSR
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_deallo
  implicit none

  private

  integer(8) :: memor(2)

  public :: iterative_solver_conjugate_gradient
  public :: iterative_solver_initialization
  public :: iterative_solver_polynomial
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-08-10
  !> @brief   Initialization
  !> @details Initialization of module variables
  !> 
  !-----------------------------------------------------------------------

  subroutine iterative_solver_initialization()

    memor = 0_8
    
  end subroutine iterative_solver_initialization
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-08-10
  !> @brief   Conjugate gradient
  !> @details Conjugate gradient
  !> 
  !-----------------------------------------------------------------------

  subroutine iterative_solver_conjugate_gradient(ndofn,nn,ia,ja,aa,bb,xx,TOLERANCE,ITERATIONS,RELATIVE_TOLERANCE)

    integer(ip),          intent(in)             :: ndofn
    integer(ip),          intent(in)             :: nn
    integer(ip), pointer, intent(in)             :: ia(:)
    integer(ip), pointer, intent(in)             :: ja(:)
    real(rp),    pointer, intent(in)             :: aa(:)
    real(rp),    pointer, intent(in)             :: bb(:)
    real(rp),    pointer, intent(inout)          :: xx(:)
    real(rp),             intent(in),   optional :: TOLERANCE
    integer(ip),          intent(in),   optional :: ITERATIONS
    logical(lg),                        optional :: RELATIVE_TOLERANCE
    integer(ip)                                  :: nsize,maxiter,iiter
    real(rp),    pointer                         :: rr(:),pp(:),zz(:),qq(:),invdiag(:)
    real(rp)                                     :: alpha,beta,rho,newrho,toler,resid,bnorm

    nsize = ndofn * nn
    nullify(rr,pp,zz,qq,invdiag)
    
    call memory_alloca(memit,'RR',     'krylov_solver_conjugate_gradient',rr,nsize)
    call memory_alloca(memit,'PP',     'krylov_solver_conjugate_gradient',pp,nsize)
    call memory_alloca(memit,'ZZ',     'krylov_solver_conjugate_gradient',zz,nsize)
    call memory_alloca(memit,'QQ',     'krylov_solver_conjugate_gradient',qq,nsize)
    call memory_alloca(memit,'INVDIAG','krylov_solver_conjugate_gradient',invdiag,nsize)

    call matrix_diagonal_CSR(nn,ndofn,0_ip,ia,ja,aa,invdiag)
    invdiag(1:nsize) = 1.0_rp/invdiag(1:nsize)
    
    call matrix_CSR_SpMV(1_ip,nn,ndofn,ndofn,ndofn,ndofn,ia,ja,aa,xx,rr,OPENMP=.false.,INITIALIZATION=.true.)
    rr(1:nsize) = bb(1:nsize) - rr(1:nsize) 
    zz(1:nsize) = invdiag(1:nsize) * rr(1:nsize)
    pp(1:nsize) = zz(1:nsize)
    qq(1:nsize) = invdiag(1:nsize) * bb(1:nsize)
    newrho      = dot_product(rr,zz)       
    bnorm       = sqrt(dot_product(bb(1:nsize),qq(1:nsize)))    
    
    toler       = 1.0e-6_rp
    maxiter     = 1000
    resid       = sqrt(newrho)/bnorm
    iiter       = 0
    if( present(TOLERANCE)  ) toler   = TOLERANCE 
    if( present(ITERATIONS) ) maxiter = ITERATIONS
    if( present(RELATIVE_TOLERANCE) ) then
       if( RELATIVE_TOLERANCE ) toler = TOLERANCE * resid
    end if
    
    do while( iiter < maxiter .and. resid > toler )
       iiter = iiter + 1
       call matrix_CSR_SpMV(1_ip,nn,ndofn,ndofn,ndofn,ndofn,ia,ja,aa,pp,qq,OPENMP=.false.,INITIALIZATION=.true.)
       alpha       = dot_product(pp,qq)
       if( alpha <= 0.0_rp ) goto 10
       alpha       = newrho / alpha
       xx(1:nsize) = xx(1:nsize) + alpha * pp(1:nsize)
       rr(1:nsize) = rr(1:nsize) - alpha * qq(1:nsize)
       zz(1:nsize) = invdiag(1:nsize) * rr(1:nsize)
       rho         = newrho
       newrho      = dot_product(rr,zz)
       if( newrho <= 0.0_rp ) goto 10
       beta        = newrho / rho
       resid       = sqrt(newrho)/bnorm
       !write(90,*) iiter,resid 
       pp(1:nsize) = zz(1:nsize) + beta * pp(1:nsize)
    end do

10  continue
    
    call memory_deallo(memit,'RR',     'krylov_solver_conjugate_gradient',rr)
    call memory_deallo(memit,'PP',     'krylov_solver_conjugate_gradient',pp)
    call memory_deallo(memit,'ZZ',     'krylov_solver_conjugate_gradient',zz)
    call memory_deallo(memit,'QQ',     'krylov_solver_conjugate_gradient',qq)
    call memory_deallo(memit,'INVDIAG','krylov_solver_conjugate_gradient',invdiag)

  end subroutine iterative_solver_conjugate_gradient

  subroutine iterative_solver_polynomial(solve,aa,bb,uu)

    use mod_solver, only : solver_parallel_SpMV
    type(soltyp), intent(inout)  :: solve(1)
    real(rp),     intent(in)     :: aa(*)
    real(rp),     intent(in)     :: bb(*)
    real(rp),     intent(inout)  :: uu(*)
    real(rp),     pointer        :: vv(:)
    real(rp),     pointer        :: yy(:)
    real(rp),     pointer        :: invdiag(:)
    integer(ip)                  :: ii,kk

    nullify(vv)
    nullify(yy)
    nullify(invdiag)

    call memory_alloca(memit,'INVDIAG','solver_approximate_inverse',invdiag,solve(1) % nequa)

    call solver_diagonal(solve,aa,invdiag)
    do ii = 1,solve(1) % nequa
       invdiag(ii) = invdiag(ii) 
    end do
    do ii = 1,solve(1) % nequa
       invdiag(ii) = 1.0_rp / invdiag(ii) 
    end do

    call memory_alloca(memit,'VV','solver_approximate_inverse',vv,solve(1) % nequa)
    call memory_alloca(memit,'YY','solver_approximate_inverse',yy,solve(1) % nequa)

    do ii = 1,solve(1) % nequa
       if( solve(1) % kfl_fixno(1,ii) <= 0 ) then
          uu(ii) = bb(ii) * invdiag(ii)
          vv(ii) = uu(ii)
       else
          vv(ii) = 0.0_rp
       end if
    end do

    do kk = 1,10
       call solver_parallel_SpMV(solve(1),aa,vv,yy,OPENMP=.true.,TIMING=.false.)
       do ii = 1,solve(1) % nequa
          if( solve(1) % kfl_fixno(1,ii) <= 0 ) then
             uu(ii) = uu(ii) + yy(ii) 
             vv(ii) = yy(ii)
          end if
       end do
    end do

    call memory_deallo(memit,'VV','solver_approximate_inverse',vv)
    call memory_deallo(memit,'YY','solver_approximate_inverse',yy)
    call memory_deallo(memit,'INVDIAG','solver_approximate_inverse',invdiag)

  end subroutine iterative_solver_polynomial
  
end module mod_iterative_solver
!> @}
