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
!> @file    deflcg.f90
!> @author  Guillaume Houzeaux
!> @brief   Deflated conjugate gradient solver
!> @details Deflated conjugate gradient solver:
!>          @verbatim
!>          1. Compute r0 := b - Ax0, z0 = M-1r0, and p0 := z0
!>          2. For j = 0, 1,..., until convergence Do:
!>          3. aj   := (rj, zj)/(Apj, pj)
!>          4. xj+1 := xj + ajpj
!>          5. rj+1 := rj - ajApj
!>          6. zj+1 := M-1 rj+1
!>          7. bj   := (rj+1, zj+1)/(rj, zj)
!>          8. pj+1 := zj+1 + bjpj
!>          9. EndDo
!>          @endverbatim
!> @}
!-----------------------------------------------------------------------
subroutine deflcg( &
     nbnodes, nbvar, idprecon, maxiter, &
     eps, an, pn, kfl_cvgso, lun_cvgso, kfl_solve, &
     lun_outso, ja, ia, bb, xx )

  use def_kintyp, only       :  ip,rp,lg
  use def_master, only       :  IMASTER,INOTMASTER
  use def_solver, only       :  memit,SOL_NO_PRECOND,SOL_LINELET,&
       &                        SOL_SQUARE,SOL_DIAGONAL,SOL_MATRIX,&
       &                        resi1,resi2,solve_sol,SOL_MATRIX_HAS_CHANGED,&
       &                        SOL_SOLVER_A_DEF2
  use mod_memory
  use mod_csrdir 
  use mod_matrix
  use mod_solver,           only : solver_krylov_subspace_save
  use mod_solver,           only : solver_define_groups_from_field
  use mod_solver,           only : solver_SpMV
  use mod_solver,           only : solver_parallel_scalar_product
  use mod_direct_solver,    only : direct_solver_partialcleaning
  use mod_direct_solver,    only : direct_solver_factorization
  use mod_direct_solver,    only : direct_solver_solution
  use mod_direct_solver,    only : direct_solver_matrix_size
  use mod_direct_solver,    only : direct_allocate_temporary_matrix
  use mod_iterative_solver, only : iterative_solver_conjugate_gradient
  use mod_deflated_cg,      only : matgro
  use mod_deflated_cg,      only : matgr2
  use mod_deflated_cg,      only : wtvect
  use mod_deflated_cg,      only : wvect
  implicit none
  
  integer(ip), intent(in)    :: nbnodes,nbvar,idprecon,maxiter
  integer(ip), intent(in)    :: kfl_cvgso,lun_cvgso
  integer(ip), intent(in)    :: kfl_solve,lun_outso
  real(rp),    intent(in)    :: eps
  real(rp),    intent(in)    :: an(nbvar,nbvar,*),pn(*)
  integer(ip), intent(in)    :: ja(*),ia(*)
  real(rp),    intent(in)    :: bb(*)
  real(rp),    intent(inout) :: xx(nbvar*nbnodes)
  integer(ip)                :: ii,nrows,ierr,ngrou,ncols,nsize
  integer(ip)                :: acoarse_size
  real(rp)                   :: alpha,beta,rho,stopcri,resid
  real(rp)                   :: invnb,newrho,dummr,time1,time2,denom
  real(rp),    pointer       :: rr(:),pp(:),zz(:),invdiag(:)
  real(rp),    pointer       :: ww(:),mu(:),p0(:),wk(:)
  real(rp),    pointer       :: murhs(:),acoarse(:)

#ifdef EVENT
  call mpitrace_user_function(1)
#endif

  ngrou = solve_sol(1) % ngrou
  !
  ! Skyline, CSR
  !
  if( ngrou > 0 ) then
     acoarse_size  = direct_solver_matrix_size(solve_sol(1) % direct_solver_Deflation)
  end if

  if( IMASTER ) then
     nrows = 0
     ncols = 0
  else
     nrows = nbnodes * nbvar
     ncols = solve_sol(1) % ncols * nbvar
  end if
  nsize = max(1_ip,ncols,nrows)
  !
  ! Sequential and slaves: Working arrays
  !
  nullify(rr)
  nullify(pp)
  nullify(p0)
  nullify(zz)
  nullify(ww)
  nullify(invdiag)
  nullify(mu)
  nullify(murhs)
  nullify(acoarse)
  nullify(wk)
  call memory_alloca(memit,'RR'     ,'deflcg',rr     ,nsize)
  call memory_alloca(memit,'PP'     ,'deflcg',pp     ,nsize)
  call memory_alloca(memit,'P0'     ,'deflcg',p0     ,nsize)
  call memory_alloca(memit,'ZZ'     ,'deflcg',zz     ,nsize)
  call memory_alloca(memit,'WW'     ,'deflcg',ww     ,nsize)
  call memory_alloca(memit,'INVDIAG','deflcg',invdiag,nsize)
  if( solve_sol(1) % kfl_algso == SOL_SOLVER_A_DEF2 ) call memory_alloca(memit,'WK'     ,'deflcg',wk     ,nsize)

  if( ngrou > 0 ) then

     !----------------------------------------------------------------------
     !
     ! Allocate memory for groups
     !
     !----------------------------------------------------------------------

     call memory_alloca(memit,'MU'     ,'deflcg',mu      , nbvar*ngrou)
     call memory_alloca(memit,'MURHS'  ,'deflcg',murhs   , nbvar*ngrou+1_ip)
     call memory_alloca(memit,'ACOARSE','deflcg',acoarse , acoarse_size)

     !----------------------------------------------------------------------
     !
     ! Compute A'= W^T.A.W and factorize A'
     !
     !----------------------------------------------------------------------

     if( solve_sol(1) % kfl_update_precond == 1 ) then
        call cputim(time1)
        if( solve_sol(1) % kfl_symme == 1 ) then
           call matgro(ngrou,nbnodes,acoarse_size,nbvar,ia,ja,an,acoarse)
        else
           call matgr2(ngrou,nbnodes,acoarse_size,nbvar,ia,ja,an,acoarse)
        end if
        if( INOTMASTER .and. solve_sol(1) % kfl_defso == 0 ) then
           call direct_solver_factorization(solve_sol(1) % direct_solver_Deflation,acoarse)
        end if
        call cputim(time2)
        solve_sol(1) % cputi(2) = solve_sol(1) % cputi(2) + time2 - time1
     end if

     !----------------------------------------------------------------------
     !
     ! X_0: Compute pre-initial guess
     !
     !----------------------------------------------------------------------

     if( solve_sol(1) % kfl_schum == 1 ) then
        call bcsrax_schur( 1_ip,nbnodes, nbvar, solve_sol(1) % ndofn_A3 , &
             solve_sol(1) % A1, solve_sol(1) % A2, solve_sol(1) % invA3, &
             solve_sol(1) % A4, ja, ia, xx, rr )
     else
        call solver_SpMV(solve_sol(1),an,xx,rr)                                      ! A.x_{-1}
     end if
     do ii= 1, nrows
        rr(ii) = bb(ii) - rr(ii)                                                     ! r_{-1} = b-A.x_{-1}
     end do
     
     call wtvect(nbnodes,ngrou,nbvar,murhs,rr)                                       ! W^T.r_{-1}

     if( INOTMASTER ) then
       if( solve_sol(1) % kfl_defso == 0 ) then
           call direct_solver_solution(solve_sol(1) % direct_solver_Deflation,murhs,mu)
        else
           call iterative_solver_conjugate_gradient(&
                nbvar,ngrou,solve_sol(1) % direct_solver_Deflation % ia,&
                solve_sol(1) % direct_solver_Deflation % ja,acoarse,murhs,mu,&
                TOLERANCE=1.0e-5_rp,ITERATIONS=1000_ip,RELATIVE_TOLERANCE=.true.)
        end if
     end if

     call wvect(nbnodes,nbvar,mu,rr)                                                 ! W.mu

#ifdef BLAS
     if( INOTMASTER ) call DAXPY(nrows,1.0_rp,rr,1_ip,xx,1_ip)
#else
     do ii = 1,nrows
        xx(ii) = xx(ii) + rr(ii)                                                        ! x0 = x_{-1} + W.mu
     end do
#endif

  end if
  
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

  !----------------------------------------------------------------------
  !
  ! Deflated: Modify initial pp
  !
  !----------------------------------------------------------------------

  if( ngrou > 0 ) then
     !
     ! Solve A'.mu = W^T.A.r^0
     !
     if( solve_sol(1) % kfl_schum == 1 ) then
        call bcsrax_schur( 1_ip, nbnodes, nbvar, solve_sol(1) % ndofn_A3 , &
             solve_sol(1) % A1, solve_sol(1) % A2, solve_sol(1) % invA3, &
             solve_sol(1) % A4, ja, ia, zz, ww )
     else
        call solver_SpMV(solve_sol(1),an,zz,ww)                   ! A.z^0
     end if
    if( solve_sol(1) % kfl_algso == SOL_SOLVER_A_DEF2 ) then
        do ii = 1,nrows
           ww(ii) = ww(ii) - rr(ii)                               ! A.z^0-r^0
        end do
     end if

     call wtvect(nbnodes,ngrou,nbvar,murhs,ww)

     if( INOTMASTER ) then
        if( solve_sol(1) % kfl_defso == 0 ) then
           call direct_solver_solution(solve_sol(1) % direct_solver_Deflation,murhs,mu)
        else
           call iterative_solver_conjugate_gradient(&
                nbvar,ngrou,solve_sol(1) % direct_solver_Deflation % ia,&
                solve_sol(1) % direct_solver_Deflation % ja,acoarse,murhs,mu,&
                TOLERANCE=1.0e-5_rp,ITERATIONS=1000_ip,RELATIVE_TOLERANCE=.true.)
        end if
     end if

     call wvect(nbnodes,nbvar,mu,ww)
     !
     ! Initial p^{k+1} = z^k - W.mu^k
     !
#ifdef BLAS
     if( INOTMASTER ) call DAXPY(nrows,-1.0_rp,ww,1_ip,pp,1_ip)
     if( INOTMASTER ) call DCOPY(nrows,pp,1_ip,p0,1_ip)
#else
     do ii = 1,nrows
        pp(ii) = pp(ii) - ww(ii)
        p0(ii) = pp(ii)
     end do
     if( solve_sol(1) % kfl_algso == SOL_SOLVER_A_DEF2 ) then
        call solver_parallel_scalar_product(solve_sol(1),rr,ww,newrho)
     end if
#endif
  end if

  !-----------------------------------------------------------------
  !
  !  MAIN LOOP
  !
  !-----------------------------------------------------------------

  do while( solve_sol(1) % iters < maxiter .and. resid > stopcri )
     !
     ! q^{k+1} =  A R^-1 p^{k+1}
     !
     call precon(&
          4_ip,nbvar,nbnodes,nrows,solve_sol(1) % kfl_symme,idprecon,ia,ja,an,&
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
     ! x^{k+1} = x^k + alpha*p^{k+1}
     !
#ifdef BLAS
     if( INOTMASTER ) call DAXPY(nrows,alpha,pp,1_ip,xx,1_ip)
#else
     do ii= 1, nrows
        xx(ii) = xx(ii) + alpha * pp(ii)
     end do
#endif
     !
     ! r^{k+1} = r^k - alpha*q^{k+1}
     !
#ifdef BLAS
     if( INOTMASTER ) call DAXPY(nrows,-alpha,zz,1_ip,rr,1_ip)
#else
     do ii = 1,nrows
        rr(ii) = rr(ii) - alpha * zz(ii)
     end do
#endif
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
          3_ip,nbvar,nbnodes,nrows,solve_sol(1) % kfl_symme,idprecon,ia,ja,an,&
          pn,invdiag,ww,rr,zz)
     !
     ! Solve A'.mu = W^T.A.z^k
     !
     if( ngrou > 0 ) then
        if( solve_sol(1) % kfl_algso == SOL_SOLVER_A_DEF2 ) then
           call solver_SpMV(solve_sol(1),an,zz,wk) 
           do ii = 1,nrows
              wk(ii) = wk(ii) - rr(ii)
           end do
           call wtvect(nbnodes,ngrou,nbvar,murhs,wk)
           if( INOTMASTER ) then
              call direct_solver_solution(solve_sol(1) % direct_solver_Deflation,murhs,mu)
           end if
           call wvect(nbnodes,nbvar,mu,ww)
           do ii = 1,nrows
              wk(ii) = zz(ii) - ww(ii)
           end do
           call solver_parallel_scalar_product(solve_sol(1),rr,wk,newrho)

        else
           ! A.z^k, rho=r.z, W^T.A.z^k
           call bsyjes(nbnodes, nbvar, an, zz, ww ,rr, ngrou, newrho, murhs)
           
           call cputim(time1)
           
           if( INOTMASTER ) then
              if( solve_sol(1) % kfl_defso == 0 ) then
                 call direct_solver_solution(solve_sol(1) % direct_solver_Deflation,murhs,mu)
              else
                 call iterative_solver_conjugate_gradient(&
                      nbvar,ngrou,solve_sol(1) % direct_solver_Deflation % ia,&
                      solve_sol(1) % direct_solver_Deflation % ja,acoarse,murhs,mu,&
                      TOLERANCE=1.0e-5_rp,ITERATIONS=1000_ip,RELATIVE_TOLERANCE=.true.)
              end if              
           end if
           call wvect(nbnodes,nbvar,mu,ww)
           
           call cputim(time2)
           solve_sol(1) % cputi(4) = solve_sol(1) % cputi(4) + time2 - time1
        end if
        
     else
        call solver_parallel_scalar_product(solve_sol(1),rr,zz,newrho)
     end if
     !
     ! beta  = rho^k / rho^{k-1}
     !
     beta = newrho / rho
     !
     ! p^{k+1} = z^k + beta*p^k - W.mu^k
     !
     if( ngrou > 0 ) then
#ifdef BLAS
        if( INOTMASTER ) call DAXPY(nrows,   beta,pp,1_ip,zz,1_ip)
        if( INOTMASTER ) call DAXPY(nrows,-1.0_rp,ww,1_ip,zz,1_ip)
        if( INOTMASTER ) call DCOPY(nrows,zz,1_ip,pp,1_ip)
#else
        if( solve_sol(1) % kfl_algso == SOL_SOLVER_A_DEF2 ) then
           do ii= 1, nrows
              pp(ii) = wk(ii) + beta * pp(ii) 
           end do
        else
           do ii= 1, nrows
              pp(ii) = zz(ii) + beta * pp(ii) - ww(ii)
           end do
        end if
#endif
     else
#ifdef BLAS
        if( INOTMASTER ) call DAXPY(nrows,   beta,pp,1_ip,zz,1_ip)
        if( INOTMASTER ) call DCOPY(nrows,zz,1_ip,pp,1_ip)
#else
        do ii= 1, nrows
           pp(ii) = zz(ii) + beta * pp(ii)
        end do
#endif
     end if

     resid  = sqrt(newrho)
     resi2  = resi1
     resi1  = resid * invnb
     solve_sol(1) % iters  = solve_sol(1) % iters + 1
     !
     ! Convergence post process:
     ! kk    = iteration number
     ! resi1 = preconditioned residual
     !
     if( kfl_cvgso == 1 ) call outcso(an,bb,xx)

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
     if( ierr > 0 ) write(lun_outso,120) solve_sol(1) % iters
  end if

  !-----------------------------------------------------------------
  !
  ! Deallocate memory
  !
  !-----------------------------------------------------------------

  if( ngrou > 0 ) then

     call memory_deallo(memit,'ACOARSE','deflcg',acoarse)
     call memory_deallo(memit,'MU'     ,'deflcg',mu)
     call memory_deallo(memit,'MURHS'  ,'deflcg',murhs)
     if( solve_sol(1) % kfl_clean_precond == 1 ) &
          call direct_solver_partialcleaning(solve_sol(1) % direct_solver_Deflation)

  end if

  call memory_deallo(memit,'INVDIAG','deflcg',invdiag)
  call memory_deallo(memit,'WW'     ,'deflcg',ww)
  call memory_deallo(memit,'ZZ'     ,'deflcg',zz)
  call memory_deallo(memit,'P0'     ,'deflcg',p0)
  call memory_deallo(memit,'PP'     ,'deflcg',pp)
  call memory_deallo(memit,'RR'     ,'deflcg',rr)
  if( solve_sol(1) % kfl_algso == SOL_SOLVER_A_DEF2 ) call memory_deallo(memit,'WK'     ,'deflcg',wk)

120 format(&
       & '# Error at iteration ',i6,&
       & 'Dividing by zero: alpha = rho^k / <p^{k+1},q^{k+1}>')

#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine deflcg
