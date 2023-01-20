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
!> @file    bcgpls.f90
!> @author  Guillaume Houzeaux
!> @brief   BICGSTAB preconditioned solver
!> @details BICGSTAB preconditioned solver
!> @}
!-----------------------------------------------------------------------
subroutine bcgpls(&
     nbnodes, nbvar, idprecon, maxiter, eps, an, pn, &
     kfl_cvgso, lun_cvgso, kfl_solve, lun_outso,&
     ja, ia, bb, xx )

  use def_kintyp, only       :  ip,rp,lg
  use def_master, only       :  IMASTER,INOTSLAVE
  use def_solver, only       :  memit
  use def_solver, only       :  resi1,resi2,solve_sol
  use mod_solver, only       :  solver_parallel_scalar_product
  use mod_memory, only       :  memory_alloca
  use mod_memory, only       :  memory_deallo
  implicit none
  integer(ip), intent(in)    :: nbnodes,nbvar,idprecon,maxiter
  integer(ip), intent(in)    :: kfl_cvgso,lun_cvgso
  integer(ip), intent(in)    :: kfl_solve,lun_outso
  real(rp),    intent(inout) :: eps
  real(rp),    intent(in)    :: an(nbvar,nbvar,*),pn(*)
  integer(ip), intent(in)    :: ja(*),ia(*)
  real(rp),    intent(in)    :: bb(*)
  real(rp),    intent(inout) :: xx(*)
  real(rp),    pointer       :: rr(:),r0(:),pp(:),tt(:),qq(:),ss(:)
  real(rp),    pointer       :: ww(:),invdiag(:)
  integer(ip)                :: ii,nrows,ierr,ncols
  real(rp)                   :: alpha,rho,newrho,beta,omega
  real(rp)                   :: raux,stopcri,invnb,resid,dummr
  
  nullify(rr,r0,pp,tt,qq,ss,ww,invdiag)

  if( IMASTER ) then
     nrows = 0
     ncols = 0
  else
     nrows = nbnodes * nbvar
     ncols = solve_sol(1) % ncols * nbvar
  end if
  !
  ! Allocate memory for working arrays
  !
  call memory_alloca(memit,'RR'     ,'bcgpls',rr,     max(1_ip,ncols,nrows))
  call memory_alloca(memit,'R0'     ,'bcgpls',r0,     max(1_ip,ncols,nrows))
  call memory_alloca(memit,'PP'     ,'bcgpls',pp,     max(1_ip,ncols,nrows))
  call memory_alloca(memit,'TT'     ,'bcgpls',tt,     max(1_ip,ncols,nrows))
  call memory_alloca(memit,'QQ'     ,'bcgpls',qq,     max(1_ip,ncols,nrows))
  call memory_alloca(memit,'SS'     ,'bcgpls',ss,     max(1_ip,ncols,nrows))
  call memory_alloca(memit,'WW'     ,'bcgpls',ww,     max(1_ip,ncols,nrows))
  call memory_alloca(memit,'INVDIAG','bcgpls',invdiag,max(1_ip,nrows))

  !----------------------------------------------------------------------
  !
  ! Initial computations
  !
  !----------------------------------------------------------------------

  call solope(&
       1_ip, nbvar, idprecon, eps, an, pn, ja, ia, bb, xx , &
       ierr, stopcri, newrho, resid, invnb, rr, r0, pp, ww, &
       invdiag, dummr )

  if( ierr /= 0 ) goto 10

  alpha = 0.0_rp
  omega = 1.0_rp
  rho   = 1.0_rp
  do ii = 1,nrows
     pp(ii) = 0.0_rp
     qq(ii) = 0.0_rp
  end do

  !----------------------------------------------------------------------
  !
  ! MAIN LOOP
  !
  !----------------------------------------------------------------------

  do while( solve_sol(1) % iters < maxiter .and. resid > stopcri )
     !
     ! beta = (rho^k/rho^{k-1})*(alpha/omega)
     !
     if( rho /= 0.0_rp .and. omega /= 0.0_rp ) then
        beta = (newrho/rho) * (alpha/omega)
     else
        ierr = 1
        goto 10
     end if
     !
     ! p^{k+1} = r^{k} + beta*(p^k - omega*q^k)
     !
     do ii = 1,nrows
        pp(ii) = rr(ii) + beta * (pp(ii) - omega * qq(ii))
     end do
     !do ii = 1,nrows
     !   if( solve_sol(1) % kfl_fixno(1,ii) > 0 ) then
     !      rr(ii) = 0.0_rp
     !      pp(ii) = 0.0_rp
     !      qq(ii) = 0.0_rp
     !      ss(ii) = 0.0_rp
     !      tt(ii) = 0.0_rp
     !   end if
     !end do
     !
     ! L q^{k+1} = A ( R^{-1} p^{k+1} )
     !
     call precon(&
          1_ip,nbvar,nbnodes,nrows,solve_sol(1)%kfl_symme,idprecon,ia,ja,an,&
          pn,invdiag,ww,pp,qq)
     !
     ! alpha = rho^k / <r0,q^{k+1}>
     !
     call solver_parallel_scalar_product(solve_sol(1),r0,qq,alpha)

     if( alpha == 0.0_rp ) then
        ierr = 2
        goto 10
     end if
     alpha = newrho / alpha
     !
     ! s = r^k - alpha*q^{k+1}
     !
     do ii = 1,nrows
        ss(ii) = rr(ii) - alpha * qq(ii)
     end do
     !
     ! L t = A ( R^{-1} s )
     !
     call precon(&
          1_ip,nbvar,nbnodes,nrows,solve_sol(1)%kfl_symme,idprecon,ia,ja,an,& ! OJO HOY
          pn,invdiag,ww,ss,tt)
     !
     ! omega = <t,s> / <t,t>
     !
     call solver_parallel_scalar_product(solve_sol(1),tt,ss,tt,omega,raux)
     
     if( raux == 0.0_rp ) then
        ierr = 2
        goto 10
     end if
     omega = omega / raux
     !
     ! x^{k+1} = x^k + alpha*p^{k+1} + omega*s
     !
     do ii =1,nrows
        xx(ii) = xx(ii) + alpha * pp(ii) + omega * ss(ii)
     end do
     !
     ! r^{k+1} = s - omega*t
     !
     do ii = 1,nrows
        rr(ii) = ss(ii) - omega * tt(ii)
     end do
     !
     ! rho^k = <r0,r^k> and || r^{k+1} ||
     !
     rho = newrho
     call solver_parallel_scalar_product(solve_sol(1),rr,rr,r0,resid,newrho)

     resid = sqrt(resid)
     resi2 = resi1
     resi1 = resid * invnb
     solve_sol(1) % iters = solve_sol(1) % iters + 1
     !
     ! Convergence output
     !
     if( kfl_cvgso == 1 ) call outcso(an,bb,xx)

  end do 

  !----------------------------------------------------------------------
  !
  ! END MAIN LOOP
  !
  !----------------------------------------------------------------------

10 continue
  call solope(&
       2_ip, nbvar, idprecon, dummr, an, dummr, ja, ia, &
       bb,xx ,ierr, dummr, dummr, dummr, invnb, rr, dummr, &
       dummr, dummr, invdiag, dummr )

  if( kfl_solve == 1 .and. INOTSLAVE ) then
     if( ierr == 1 ) write(lun_outso,201) solve_sol(1) % iters
     if( ierr == 2 ) write(lun_outso,202) solve_sol(1) % iters
  end if

  call memory_deallo(memit,'RR'     ,'bcgpls',rr)
  call memory_deallo(memit,'R0'     ,'bcgpls',r0)
  call memory_deallo(memit,'PP'     ,'bcgpls',pp)
  call memory_deallo(memit,'TT'     ,'bcgpls',tt)
  call memory_deallo(memit,'QQ'     ,'bcgpls',qq)
  call memory_deallo(memit,'SS'     ,'bcgpls',ss)
  call memory_deallo(memit,'WW'     ,'bcgpls',ww)
  call memory_deallo(memit,'INVDIAG','bcgpls',invdiag)

100 format(i7,1x,e12.6)
110 format(i5,18(2x,e12.6))
201 format(&
       & '# Error at iteration ',i6,&
       & 'Dividing by zero: alpha = beta = (rho^k/rho^{k-1})*(alpha/omega)')
202 format(&
       & '# Error at iteration ',i6,&
       & 'Dividing by zero: alpha = beta = alpha = rho^k / <r0,q^{k+1}>')

end subroutine bcgpls
