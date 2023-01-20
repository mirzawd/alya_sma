!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




    !-------------------------------------------------------------------
    !
    ! MINRES  is designed to solve the system of linear equations
    !
    !    Ax = b
    !
    ! or the least-squares problem
    !
    !    min ||Ax - b||_2,
    !
    ! where A is an n by n symmetric matrix and b is a given vector.
    ! The matrix A may be indefinite and/or singular.
    !
    ! 1. If A is known to be positive definite, the Conjugate Gradient
    ! Method might be preferred, since it requires the same number
    ! of iterations as MINRES but less work per iteration.
    !
    ! 2. If A is indefinite but Ax = b is known to have a solution
    ! (e.g. if A is nonsingular), SYMMLQ might be preferred,
    ! since it requires the same number of iterations as MINRES
    ! but slightly less work per iteration.
    !
    ! The matrix A is intended to be large and sparse.  It is accessed
    ! by means of a subroutine call of the form
    ! SYMMLQ development:
    !
    !    call Aprod ( n, x, y )
    !
    ! which must return the product y = Ax for any given vector x.
    !
    !
    ! More generally, MINRES is designed to solve the system
    !
    !    (A - shift*I) x = b
    ! or
    !    min ||(A - shift*I) x - b||_2,
    !
    ! where  shift  is a specified scalar value.  Again, the matrix
    ! (A - shift*I) may be indefinite and/or singular.
    ! The work per iteration is very slightly less if  shift = 0.
    !
    ! Note: If  shift  is an approximate eigenvalue of  A
    ! and  b  is an approximate eigenvector,  x  might prove to be
    ! a better approximate eigenvector, as in the methods of
    ! inverse iteration and/or Rayleigh-quotient iteration.
    ! However, we're not yet sure on that -- it may be better to use SYMMLQ.
    !
    ! A further option is that of preconditioning, which may reduce
    ! the number of iterations required.  If M = C C' is a positive
    ! definite matrix that is known to approximate  (A - shift*I)
    ! in some sense, and if systems of the form  My = x  can be
    ! solved efficiently, the parameters precon and Msolve may be
    ! used (see below).  When  precon = .true., MINRES will
    ! implicitly solve the system of equations
    !
    !    P (A - shift*I) P' xbar  =  P b,
    !
    ! i.e.             Abar xbar  =  bbar
    ! where                    P  =  C**(-1),
    !                       Abar  =  P (A - shift*I) P',
    !                       bbar  =  P b,
    !
    ! and return the solution       x  =  P' xbar.
    ! The associated residual is rbar  =  bbar - Abar xbar
    !                                  =  P (b - (A - shift*I)x)
    !                                  =  P r.
    !
    ! In the discussion below, eps refers to the machine precision.
    !
    ! Parameters
    ! ----------
    !
    ! n       input      The dimension of the matrix A.
    ! b(n)    input      The rhs vector b.
    ! x(n)    output     Returns the computed solution x.
    !
    ! Aprod   external   A subroutine defining the matrix A.
    !                       call Aprod ( n, x, y )
    !                    must return the product y = Ax
    !                    without altering the vector x.
    !
    ! Msolve  external   An optional subroutine defining a
    !                    preconditioning matrix M, which should
    !                    approximate (A - shift*I) in some sense.
    !                    M must be positive definite.
    !
    !                       call Msolve( n, x, y )
    !
    !                    must solve the linear system My = x
    !                    without altering the vector x.
    !
    !                    In general, M should be chosen so that Abar has
    !                    clustered eigenvalues.  For example,
    !                    if A is positive definite, Abar would ideally
    !                    be close to a multiple of I.
    !                    If A or A - shift*I is indefinite, Abar might
    !                    be close to a multiple of diag( I  -I ).
    !
    ! checkA  input      If checkA = .true., an extra call of Aprod will
    !                    be used to check if A is symmetric.  Also,
    !                    if precon = .true., an extra call of Msolve
    !                    will be used to check if M is symmetric.
    !
    ! precon  input      If precon = .true., preconditioning will
    !                    be invoked.  Otherwise, subroutine Msolve
    !                    will not be referenced; in this case the
    !                    actual parameter corresponding to Msolve may
    !                    be the same as that corresponding to Aprod.
    !
    ! shift   input      Should be zero if the system Ax = b is to be
    !                    solved.  Otherwise, it could be an
    !                    approximation to an eigenvalue of A, such as
    !                    the Rayleigh quotient b'Ab / (b'b)
    !                    corresponding to the vector b.
    !                    If b is sufficiently like an eigenvector
    !                    corresponding to an eigenvalue near shift,
    !                    then the computed x may have very large
    !                    components.  When normalized, x may be
    !                    closer to an eigenvector than b.
    !
    ! nout    input      A file number.
    !                    If nout > 0, a summary of the iterations
    !                    will be printed on unit nout.
    !
    ! itnlim  input      An upper limit on the number of iterations.
    !
    ! rtol    input      A user-specified tolerance.  MINRES terminates
    !                    if it appears that norm(rbar) is smaller than
    !                       rtol * norm(Abar) * norm(xbar),
    !                    where rbar is the transformed residual vector,
    !                       rbar = bbar - Abar xbar.
    !
    !                    If shift = 0 and precon = .false., MINRES
    !                    terminates if norm(b - A*x) is smaller than
    !                       rtol * norm(A) * norm(x).
    !
    ! istop   output     An integer giving the reason for termination...
    !
    !          -1        beta2 = 0 in the Lanczos iteration; i.e. the
    !                    second Lanczos vector is zero.  This means the
    !                    rhs is very special.
    !                    If there is no preconditioner, b is an
    !                    eigenvector of A.
    !                    Otherwise (if precon is true), let My = b.
    !                    If shift is zero, y is a solution of the
    !                    generalized eigenvalue problem Ay = lambda My,
    !                    with lambda = alpha1 from the Lanczos vectors.
    !
    !                    In general, (A - shift*I)x = b
    !                    has the solution         x = (1/alpha1) y
    !                    where My = b.
    !
    !           0        b = 0, so the exact solution is x = 0.
    !                    No iterations were performed.
    !
    !           1        Norm(rbar) appears to be less than
    !                    the value  rtol * norm(Abar) * norm(xbar).
    !                    The solution in  x  should be acceptable.
    !
    !           2        Norm(rbar) appears to be less than
    !                    the value  eps * norm(Abar) * norm(xbar).
    !                    This means that the residual is as small as
    !                    seems reasonable on this machine.
    !
    !           3        Norm(Abar) * norm(xbar) exceeds norm(b)/eps,
    !                    which should indicate that x has essentially
    !                    converged to an eigenvector of A
    !                    corresponding to the eigenvalue shift.
    !
    !           4        Acond (see below) has exceeded 0.1/eps, so
    !                    the matrix Abar must be very ill-conditioned.
    !                    x may not contain an acceptable solution.
    !
    !           5        The iteration limit was reached before any of
    !                    the previous criteria were satisfied.
    !
    !           6        The matrix defined by Aprod does not appear
    !                    to be symmetric.
    !                    For certain vectors y = Av and r = Ay, the
    !                    products y'y and r'v differ significantly.
    !
    !           7        The matrix defined by Msolve does not appear
    !                    to be symmetric.
    !                    For vectors satisfying My = v and Mr = y, the
    !                    products y'y and r'v differ significantly.
    !
    !           8        An inner product of the form  x' M**(-1) x
    !                    was not positive, so the preconditioning matrix
    !                    M does not appear to be positive definite.
    !
    !                    If istop >= 5, the final x may not be an
    !                    acceptable solution.
    !
    ! itn     output     The number of iterations performed.
    !
    ! Anorm   output     An estimate of the norm of the matrix operator
    !                    Abar = P (A - shift*I) P',   where P = C**(-1).
    !
    ! Acond   output     An estimate of the condition of Abar above.
    !                    This will usually be a substantial
    !                    under-estimate of the true condition.
    !
    ! rnorm   output     An estimate of the norm of the final
    !                    transformed residual vector,
    !                       P (b  -  (A - shift*I) x).
    !
    ! ynorm   output     An estimate of the norm of xbar.
    !                    This is sqrt( x'Mx ).  If precon is false,
    !                    ynorm is an estimate of norm(x).
    !-------------------------------------------------------------------
    ! MINRES is an implementation of the algorithm described in
    ! the following reference:
    !
    ! C. C. Paige and M. A. Saunders (1975),
    ! Solution of sparse indefinite systems of linear equations,
    ! SIAM J. Numer. Anal. 12(4), pp. 617-629.
    !-------------------------------------------------------------------
    !
    !
    ! MINRES development:
    !    1972: First version, similar to original SYMMLQ.
    !          Later lost @#%*!
    !    Oct 1995: Tried to reconstruct MINRES from
    !              1995 version of SYMMLQ.
    ! 30 May 1999: Need to make it more like LSQR.
    !              In middle of major overhaul.
    ! 19 Jul 2003: Next attempt to reconstruct MINRES.
    !              Seems to need two vectors more than SYMMLQ.  (w1, w2)
    !              Lanczos is now at the top of the loop,
    !              so the operator Aprod is called in just one place
    !              (not counting the initial check for symmetry).
    ! 22 Jul 2003: Success at last.  Preconditioning also works.
    !              minres.f added to http://www.stanford.edu/group/SOL/.
    !
    ! 16 Oct 2007: Added a stopping rule for singular systems,
    !              as derived in Sou-Cheng Choi's PhD thesis.
    !              Note that ||Ar|| small => r is a null vector for A.
    !              Subroutine minrestest2 in minresTestModule.f90
    !              tests this option.  (NB: Not yet working.)
    !-------------------------------------------------------------------

subroutine minres(  &
     nbnodes, nbvar, idprecon, maxiter, &
     eps, an, pn, kfl_cvgso, lun_cvgso, kfl_solve, &
     lun_outso, ja, ia, bb, xx )

  use def_kintyp, only         : ip,rp,lg
  use def_master, only         : IMASTER
  use def_solver, only         : memit
  use def_solver, only         : resi1,resi2,solve_sol
  use mod_memory, only         : memory_alloca
  use mod_memory, only         : memory_deallo
  use mod_solver, only         : solver_SpMV
  use mod_solver, only         : solver_parallel_scalar_product
  implicit none

  integer(ip), intent(in)    :: nbnodes,nbvar,idprecon,maxiter
  integer(ip), intent(in)    :: kfl_cvgso,lun_cvgso
  integer(ip), intent(in)    :: kfl_solve,lun_outso
  real(rp),    intent(in)    :: eps
  real(rp),    intent(in)    :: an(nbvar,nbvar,*),pn(*)
  integer(ip), intent(in)    :: ja(*),ia(*)
  real(rp),    intent(in)    :: bb(*)
  real(rp),    intent(inout) :: xx(nbvar*nbnodes)
  integer(ip)                :: ii,nrows,ierr,npoin
  real(rp)                   :: stopcri,resid
  real(rp)                   :: invnb,dummr

  real(rp),    pointer       :: v0(:),v1(:),w0(:),w1(:),vv(:)
  real(rp),    pointer       :: invdiag(:),zz(:),ww(:),p0(:),z1(:)

  real(rp)                   :: gama1,gama0,neta,s0,s1,c0,c1,delta
  real(rp)                   :: alpha0,alpha1,alpha2,alpha3,gama

  if( IMASTER ) then
     nrows = 1                ! Minimum memory for working arrays
     npoin = 0                ! master does not perform any loop
  else
     nrows = nbnodes * nbvar
     npoin = nbnodes
  end if
  !
  ! Sequential and slaves: Working arrays
  ! 
  nullify(v0)
  nullify(v1)
  nullify(vv)
  nullify(w0)
  nullify(w1)
  nullify(z1)
  nullify(zz)
  nullify(ww)
  nullify(p0)
  nullify(invdiag)
  call memory_alloca(memit,'V0'     ,'deflcg',v0     ,nrows)
  call memory_alloca(memit,'V1'     ,'deflcg',v1     ,nrows)
  call memory_alloca(memit,'VV'     ,'deflcg',vv     ,nrows)
  call memory_alloca(memit,'W0'     ,'deflcg',w0     ,nrows)
  call memory_alloca(memit,'W1'     ,'deflcg',w1     ,nrows)
  call memory_alloca(memit,'Z1'     ,'deflcg',z1     ,nrows)
  call memory_alloca(memit,'ZZ'     ,'deflcg',zz     ,nrows)
  call memory_alloca(memit,'WW'     ,'deflcg',ww     ,nrows)
  call memory_alloca(memit,'P0'     ,'deflcg',p0     ,nrows)
  call memory_alloca(memit,'INVDIAG','deflcg',invdiag,nrows)

  if( IMASTER ) nrows = 0 ! Master does not perform any loop

  !----------------------------------------------------------------------
  !
  ! Initial computations
  ! v^1    = b - Ax^0
  ! M z^1  = v^1
  ! gama1 = (v^1,z^1)^1/2
  !
  !----------------------------------------------------------------------

  call solope(&
       1_ip, nbvar, idprecon, eps, an, pn, ja, ia, bb, xx ,    &
       ierr, stopcri, gama1, resid, invnb, v1, z1, dummr, ww, &
       invdiag, p0 )

  if( ierr /= 0 ) goto 10
 
  gama1 = sqrt(gama1)
  gama0 = 1.0_rp
  neta   = gama1
  s0     = 0.0_rp
  s1     = 0.0_rp
  c0     = 1.0_rp
  c1     = 1.0_rp

  do ii = 1,nrows
     v0(ii) = 0.0_rp
     w0(ii) = 0.0_rp
     w1(ii) = 0.0_rp
  end do

  !-----------------------------------------------------------------
  !
  !  MAIN LOOP
  !
  !-----------------------------------------------------------------

  do while( solve_sol(1) % iters < maxiter .and. resid > stopcri )
     !
     ! Obtain first Lanczos
     !
     do ii = 1,nrows
        z1(ii) = z1(ii) / gama1
     end do
     !
     ! Compute v=Az and delta = (Az,z)
     !
     call solver_SpMV(solve_sol(1),an,z1,vv)
     call solver_parallel_scalar_product(solve_sol(1),vv,z1,delta)
     !
     ! v^{i+1} = A.z^i - (delta/gama^i) * v^i - (gama^i/gama^{i-1}) * v^{i-1}
     !
     do ii = 1,nrows
        vv(ii) = vv(ii) - (delta/gama1) * v1(ii) - (gama1/gama0) * v0(ii)
     end do 
     !
     ! M z^{i+1}  = v^{i+1} 
     !
     call precon(&
          3_ip,nbvar,nbnodes,nrows,solve_sol(1) % kfl_symme,idprecon,ia,ja,an,&
          pn,invdiag,ww,vv,zz)
     !
     ! gama^{i+1} = (z^{i+1},v^{i+1})^{1/2}  
     !
     call solver_parallel_scalar_product(solve_sol(1),zz,vv,gama)
     if( gama < 0.0_rp ) then
        stop
        ierr = 6
        go to 10
     end if
     gama = sqrt(gama)
     !
     ! Update constants
     !
     alpha0 = c1 * delta - c0 * s1 * gama1
     alpha1 = sqrt( alpha0**2 + gama**2 )
     alpha2 = s1 * delta + c0 * c1 * gama1
     alpha3 = s0 * gama1
     !
     ! Save old values
     !
     gama0 = gama1
     gama1 = gama
     c0     = c1
     s0     = s1

     c1     = alpha0 / alpha1
     s1     = gama  / alpha1
     !
     ! Update solution
     !
     do ii = 1,nrows        
        ww(ii) = ( z1(ii) - alpha3 * w0(ii) - alpha2 * w1(ii) ) / alpha1
     end do
     do ii = 1,nrows        
        xx(ii) = xx(ii) + c1 * neta * ww(ii)
     end do
     do ii = 1,nrows        
        w0(ii) = w1(ii)
        w1(ii) = ww(ii)
        v0(ii) = v1(ii)
        v1(ii) = vv(ii)
        z1(ii) = zz(ii) 
     end do
     neta = -s1 * neta  

     call solver_SpMV(solve_sol(1),an,xx,ww) 
     do ii = 1,nrows
        ww(ii) = bb(ii) - ww(ii)
     end do
     call precon(&
          3_ip,nbvar,nbnodes,nrows,solve_sol(1) % kfl_symme,idprecon,ia,ja,an,&
          pn,invdiag,p0,ww,vv)
     !call norm2x(nbvar,ww,resid)
     !do ii = 1,nrows
     !   ww(ii) = xx(ii)-s1*neta*w1(ii)/c1
     !end do
     call solver_parallel_scalar_product(solve_sol(1),ww,vv,resid)
     resid = sqrt(resid)
     !resid = resid * invnb
     !
     ! Residual
     !
     !resid  = gama1
     resi2  = resi1
     resi1  = resid * invnb
     solve_sol(1) % iters = solve_sol(1) % iters + 1

     !
     ! Convergence post process:
     ! kk    = iteration number
     ! resi1 = preconditioned residual
     !
     if( kfl_cvgso == 1 ) call outcso(an,bb,xx)

     !if(imaster) write(90,*) 'h=',solve_sol(1) % iters,sqrt(resid),abs(c1),gama1
     !if( solve_sol(1) % iters == 20 ) stop

  end do

  !-----------------------------------------------------------------
  !
  !  END MAIN LOOP
  !
  !-----------------------------------------------------------------

10 continue
  call solope(&
       2_ip, nbvar, idprecon, dummr, an, dummr, ja, ia, bb, xx , &
       ierr, dummr, dummr, resi1, dummr,zz, p0, ww, dummr, &
       invdiag, dummr )

  if( kfl_solve == 1 ) then
     if( ierr == 8 ) write(lun_outso,120) solve_sol(1) % iters
  end if
  !
  ! Deallocate memory
  !
  call memory_deallo(memit,'V0'     ,'deflcg',v0     )
  call memory_deallo(memit,'V1'     ,'deflcg',v1     )
  call memory_deallo(memit,'VV'     ,'deflcg',vv     )
  call memory_deallo(memit,'W0'     ,'deflcg',w0     )
  call memory_deallo(memit,'W1'     ,'deflcg',w1     )
  call memory_deallo(memit,'Z1'     ,'deflcg',z1     )
  call memory_deallo(memit,'ZZ'     ,'deflcg',zz     )
  call memory_deallo(memit,'WW'     ,'deflcg',ww     )
  call memory_deallo(memit,'P0'     ,'deflcg',p0     )
  call memory_deallo(memit,'INVDIAG','deflcg',invdiag)
  !
  ! Display final status
  !
120 format(&
       & '# Preconditioner is not positive definite')

end subroutine minres
