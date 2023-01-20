!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Algebraic_Solver
!> @addtogroup Krylov_Solver
!> @{
!> @file    cgrpls.f90
!> @author  Guillaume Houzeaux
!> @brief   Conjugate gradient solver
!> @details Conjugate gradient solver:
!>          \verbatim
!>          Compute r0 := b - Ax0, z0 = M-1r0, and p0 := z0
!>          For j = 0, 1,..., until convergence Do:
!>            vj   = Apj
!>            aj   = (rj, zj)/(v, pj)
!>            xj+1 = xj + aj pj
!>            rj+1 = rj - aj vj
!>            zj+1 = M-1 rj+1
!>            bj   = (rj+1, zj+1)/(rj, zj)
!>            pj+1 = zj+1 + bjpj
!>          EndDo
!>          \endverbatim
!> @}
!-----------------------------------------------------------------------
subroutine cgrpls( &
     nbnodes, nbvar, idprecon, maxiter, &
     eps, an, pn, kfl_cvgso, lun_cvgso, kfl_solve, &
     lun_outso, ja, ia, bb, xx )

  use mod_couplings

  use def_kintyp, only       :  ip,rp,lg
  use def_master, only       :  IMASTER
  use def_master, only       :  NPOIN_TYPE
  use def_solver, only       :  memit,SOL_NO_PRECOND,SOL_LINELET
  use def_solver, only       :  SOL_SQUARE,SOL_DIAGONAL,SOL_MATRIX
  use def_solver, only       :  resi1,resi2,solve_sol
  use mod_solver, only       :  solver_krylov_subspace_save
  use mod_solver, only       :  solver_parallel_scalar_product
!  use def_domain, only       :  npoin                      ! CHANGE
  use mod_memory, only       :  memory_alloca
  use mod_memory, only       :  memory_deallo
  use mod_solver, only       :  solver_invdiag
  use mod_solver, only       :  solver_SpMV
  implicit none
  integer(ip), intent(in)    :: nbnodes, nbvar, idprecon, maxiter
  integer(ip), intent(in)    :: kfl_cvgso, lun_cvgso
  integer(ip), intent(in)    :: kfl_solve, lun_outso
  real(rp),    intent(in)    :: eps
  real(rp),    intent(in)    :: an(nbvar,nbvar,*), pn(*)
  integer(ip), intent(in)    :: ja(*), ia(*)
  real(rp),    intent(in)    :: bb(*)
  real(rp),    intent(inout) :: xx(nbvar*nbnodes)
  integer(ip)                :: ii, nrows, ierr, ncols
!  integer(ip)                :: kcoup, kk
  real(rp)                   :: alpha, beta, rho, stopcri, resid
  real(rp)                   :: invnb, newrho, dummr, denom
  real(rp),    pointer       :: rr(:), pp(:), zz(:), invdiag(:)
  real(rp),    pointer       :: ww(:), p0(:)

#ifdef EVENT
  call mpitrace_user_function(1)
#endif

  nullify(rr,pp,p0,zz,ww,invdiag)

  if( IMASTER ) then
     nrows = 0
     ncols = 0
  else
     nrows = nbnodes * nbvar
     ncols = solve_sol(1) % ncols * nbvar
  end if
  !
  ! Sequential and slaves: Working arrays
  !
  call memory_alloca(memit,'RR'     ,'cgrpls',rr,     max(1_ip,ncols,nrows))
  call memory_alloca(memit,'PP'     ,'cgrpls',pp,     max(1_ip,ncols,nrows))
  call memory_alloca(memit,'P0'     ,'cgrpls',p0,     max(1_ip,ncols,nrows))
  call memory_alloca(memit,'ZZ'     ,'cgrpls',zz,     max(1_ip,ncols,nrows))
  call memory_alloca(memit,'WW'     ,'cgrpls',ww,     max(1_ip,ncols,nrows))
  call memory_alloca(memit,'INVDIAG','cgrpls',invdiag,max(1_ip,nrows))

  solver_invdiag => invdiag

  !----------------------------------------------------------------------
  !
  ! Initial computations
  !
  !----------------------------------------------------------------------

  call solope(&
       1_ip, nbvar, idprecon, eps, an, pn, ja, ia, bb, xx , &
       ierr, stopcri, newrho, resid, invnb, rr, zz, pp, ww, &
       invdiag, p0 )

  if( ierr /= 0 ) goto 10

  !-----------------------------------------------------------------
  !
  !  MAIN LOOP
  !
  !-----------------------------------------------------------------
100 continue

  do while( solve_sol(1) % iters < maxiter .and. resid > stopcri )
     !
     ! q^{k+1} =  A p^{k+1}
     !
     call precon(&
          4_ip,nbvar,nbnodes,nrows,solve_sol(1)%kfl_symme,idprecon,ia,ja,an,&
          pn,invdiag,ww,pp,zz)
     !
     ! alpha = rho^k / <p^{k+1},q^{k+1}>
     !
     call solver_parallel_scalar_product(solve_sol(1),pp,zz,denom)

     if( denom == 0.0_rp ) then
        ierr = 2
        goto 10
     end if
     rho   = newrho
     alpha = newrho / denom
     !
     ! x^{k+1} = x^k + alpha * p^{k+1}
     ! r^{k+1} = r^k - alpha * q^{k+1}
     !
#ifdef BLAS
     call DAXPY(nrows, alpha,pp,1_ip,xx,1_ip)
     call DAXPY(nrows,-alpha,zz,1_ip,rr,1_ip)
#else
     do ii = 1,nrows
        xx(ii) = xx(ii) + alpha * pp(ii)
        rr(ii) = rr(ii) - alpha * zz(ii)
     end do
#endif     
     !if( mod(solve_sol(1) % iters,50)== 0 ) then
     !   call bcsrax( 1_ip, npoin, nbvar, an, ja, ia, xx, rr )
     !   do ii = 1,nrows
     !      rr(ii) = bb(ii) - rr(ii)
     !   end do       
     !end if
     !
     ! Save Krylov subspace
     !
     if( solve_sol(1) % kfl_save_krylov > 0 ) then
        call solver_krylov_subspace_save(solve_sol(1),pp,denom)
     end if
     !
     !  L z^{k+1} = r^{k+1}
     !
     call precon(&
          3_ip,nbvar,nbnodes,nrows,solve_sol(1)%kfl_symme,idprecon,ia,ja,an,&
          pn,invdiag,ww,rr,zz)
     call solver_parallel_scalar_product(solve_sol(1),rr,zz,newrho)
     !
     ! beta  = rho^k / rho^{k-1}
     !
     beta = newrho / rho
     !
     ! p^{k+1} = z^k + beta*p^k
     !
#ifdef BLAS
     call DAXPY(nrows,beta,pp,1_ip,zz,1_ip)  ! z =  z + beta*p
     call DCOPY(nrows,zz,1_ip,pp,1_ip)       ! p <= z
#else
     do ii = 1,nrows
        pp(ii) = zz(ii) + beta * pp(ii)
     end do
#endif

     resid  = sqrt(newrho)
     resi2  = resi1
     resi1  = resid * invnb
     solve_sol(1) % iters = solve_sol(1) % iters + 1
!!$     !
!!$     ! Orthogonality
!!$     !
!!$     if( solve_sol(1) % kfl_schum == 1 ) then
!!$        call bcsrax_schur( 1_ip, nbnodes, nbvar, solve_sol(1) % ndofn_A3 , &
!!$             solve_sol(1) % A1, solve_sol(1) % A2, solve_sol(1) % invA3, solve_sol(1) % A4,&
!!$             invdiag, ja, ia, pp, ww )
!!$     else
!!$        call solver_SpMV(solve_sol(1),an,pp,ww) 
!!$     end if
!!$     call cosixy(nbvar,nbnodes,ww,p0,solve_sol(1) % xorth)
!!$     if( abs(solve_sol(1) % xorth) > 1.0e-10_rp ) then
!!$        ! r0 := b - Ax
!!$        call bcsrax( 1_ip, nbnodes, nbvar, an, ja, ia, xx, rr )
!!$        do ii = 1,nrows
!!$           rr(ii) = bb(ii) - rr(ii)
!!$        end do
!!$        ! p0 = u0 := M^-1 r0 
!!$        call precon(&
!!$             3_ip,nbvar,nbnodes,nrows,solve_sol(1) % kfl_symme,idprecon,ia,ja,an,&
!!$             pn,invdiag,ww,rr,pp)
!!$        do ii = 1,nrows
!!$           p0(ii) = pp(ii)
!!$        end do
!!$         ! rho = (r0,p0)
!!$        call prodxy(nbvar,nbnodes,pp,rr,newrho)
!!$        solve_sol(1) % iters = 0
!!$        goto 100        
!!$     end if
     !
     ! Convergence post process:
     ! kk    = iteration number
     ! resi1 = preconditioned residual
     !
     if( kfl_cvgso == 1 ) call outcso(an,bb,xx)
     !if(modul==7.and.ittim==2.and.solve_sol(1) % iters==112) then
     !   print*,'a1=',solve_sol(1) % iters ,maxiter ,resid , stopcri
     !   print*,'a2=',solve_sol(1) % iters < maxiter .and. resid > stopcri
     !end if

  end do

  !-----------------------------------------------------------------
  !
  !  END MAIN LOOP
  !
  !-----------------------------------------------------------------

10 continue

  call solope(&
       2_ip, nbvar, idprecon, dummr, an, dummr, ja, ia, bb, xx , &
       ierr, dummr, dummr, resi1, dummr, rr, p0, pp, dummr, &
       invdiag, dummr )

  if( kfl_solve == 1 ) then
     if( ierr > 0 ) then
        write(lun_outso,120) solve_sol(1) % iters
        call runend('CGRPLS: SOLVER HAS NOT CONVERGED! YOUR MATRIX MAY BE ILL CONDITIONED')
     end if
  end if

  !mcoup = kcoup

  !-----------------------------------------------------------------
  !
  ! Deallocate memory
  !
  !-----------------------------------------------------------------

  call memory_deallo(memit,'RR'     ,'cgrpls',rr)
  call memory_deallo(memit,'PP'     ,'cgrpls',pp)
  call memory_deallo(memit,'P0'     ,'cgrpls',p0)
  call memory_deallo(memit,'ZZ'     ,'cgrpls',zz)
  call memory_deallo(memit,'WW'     ,'cgrpls',ww)
  call memory_deallo(memit,'INVDIAG','cgrpls',invdiag)

120 format(&
       & '# Error at iteration ',i6,&
       & 'Dividing by zero: alpha = rho^k / <p^{k+1},q^{k+1}>')

#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine cgrpls
