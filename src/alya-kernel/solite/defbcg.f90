!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine defbcg( &
     nbnodes, nbvar, idprecon, maxiter, eps, an, pn, &
     kfl_cvgso_sol, lun_cvgso, kfl_solve_sol, lun_outso,&
     ja, ia, bb, xx )

   !-----------------------------------------------------------------------
  ! Objective: Solve    [A] [M]^-1  [M] x = b
  !                        [A']        x' = b'   by the BiCGstab method.
  !
  !            Three preconditioners are possible:
  !            idprecon = 0 => [M]^-1 = [M] = I
  !            idprecon = 2 => [M]^-1 = diag([A])
  !            idprecon = 3 => [M]^-1 = SPAI

  !            The diagonal terms of [A] must be the first in each row.
  !
  !            If Diag. Scaling is selected the preconditioned system is:
  !                    D^-1/2 [A] D^-1/2   D^1/2 x  =  D^-1/2 b
  !                          [A']              x'   =      b'
  !
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp,lg
  use def_master, only       :  kfl_paral,&
       &                        NPOIN_REAL_12DI
  use def_solver, only       :  memit,SOL_NO_PRECOND,SOL_LINELET,&
       &                        SOL_SQUARE,SOL_DIAGONAL,SOL_MATRIX,&
       &                        resin,resfi,resi1,resi2,solve_sol
  use mod_solver, only       :  solver_SpMV
  use mod_solver, only       :  solver_parallel_vector_L2norm
  use mod_solver, only       :  solver_parallel_scalar_product
  use mod_deflated_cg,  only :  wvect
  use mod_deflated_cg,  only :  wtvect
  use mod_skyline,      only :  lusolv
  use mod_memchk
  !use mod_postpr
  implicit none
  integer(ip), intent(in)    :: nbnodes, nbvar, idprecon, maxiter
  integer(ip), intent(in)    :: kfl_cvgso_sol, lun_cvgso
  integer(ip), intent(in)    :: kfl_solve_sol, lun_outso
  real(rp),    intent(in)    :: eps
  real(rp),    intent(in)    :: an(nbvar, nbvar, *),pn(*)
  integer(ip), intent(in)    :: ja(*), ia(*)
  real(rp),    intent(inout) :: bb(*)
  real(rp),    intent(inout) :: xx(nbvar*nbnodes)
  integer(ip)                :: ii, jj, kk, ll, nrows, ierr, npopo, ngrou
  integer(ip)                :: nskyl , info
  integer(4)                 :: istat
  real(rp)                   :: alpha, beta, rho, stopcri, resid
  real(rp)                   :: invnb, newrho ,raux , omega 
  real(rp),    pointer       :: rr(:), r0(:), pp(:), tt(:), qq(:), ss(:)
  real(rp),    pointer       :: wa1(:), wa2(:), mu(:), askyl(:)
!  real(rp)                   :: kaka=1.0_rp
!  character(5) :: wopos(2)

  ngrou =  solve_sol(1)%ngrou
  nskyl =  solve_sol(1)%nskyl
  ierr  =  0
  resi1 =  1.0_rp
  if(kfl_paral==0) then
     nrows   = 1                ! Minimum memory for working arrays
     npopo   = 0                ! master does not perform any loop
  else
     npopo   = nbnodes
     nrows   = nbnodes * nbvar
  end if
  !
  ! Sequential and slaves: Working arrays
  !
  allocate(rr(nrows),stat=istat) 
  call memchk(0_ip,istat,memit,'RR','defbcg',rr)
  allocate(r0(nrows),stat=istat)
  call memchk(0_ip,istat,memit,'R0','defbcg',r0)
  allocate(pp(nrows),stat=istat)  
  call memchk(0_ip,istat,memit,'PP','defbcg',pp)
  allocate(tt(nrows),stat=istat)
  call memchk(0_ip,istat,memit,'TT','defbcg',tt)
  allocate(qq(nrows),stat=istat) 
  call memchk(0_ip,istat,memit,'QQ','defbcg',qq)
  allocate(ss(nrows),stat=istat)
  call memchk(0_ip,istat,memit,'SS','defbcg',ss)
  allocate(wa1(nrows),stat=istat)
  call memchk(0_ip,istat,memit,'WA1','defbcg',wa1)
  allocate(wa2(nrows),stat=istat) 
  call memchk(0_ip,istat,memit,'WA2','defbcg',wa2)
  

  if(ngrou/=0) then
     allocate(mu(ngrou*nbvar),stat=istat)
     call memchk(0_ip,istat,memit,'MU','defbcg',mu)
     allocate(askyl(nskyl),stat=istat)
     call memchk(0_ip,istat,memit,'ASKYL','defbcg',askyl)
  end if
  if(kfl_paral==0) nrows = 0 ! Master does not perform any loop

  if(idprecon==SOL_SQUARE.or.idprecon==SOL_DIAGONAL) then
     !
     ! wa1 = D^1/2
     !
     do ii= 1, npopo
        jj = ia(ii)
        ll = -1
        do while (jj< ia(ii+1) .and. ll ==-1)
           if(ja(jj)==ii) then
              ll = jj
           end if
           jj = jj+1
        end do
        if(ll/=-1) then
           jj = (ii-1) * nbvar
           do kk= 1, nbvar
              wa1(jj+kk)=an(kk,kk,ll)
           end do
        end if
     end do
     if(idprecon==SOL_SQUARE) then
        do ii= 1, nrows
           wa1(ii) = sqrt(abs(wa1(ii)))
        end do
     end if
     !
     ! Periodicity and Parallelization
     !
     call rhsmod(nbvar,wa1)
     !
     ! WA2 = M^{-1} = 1/diag(A)
     !
     do ii= 1, nrows
        wa2(ii) = 1.0_rp / wa1(ii)
     end do
   endif 
  !
  ! Compute A'= W^T.A.W
  !
  call magru2(ngrou,npopo,nskyl,nbvar,ia,ja,an,askyl,idprecon,wa2)
  !
  ! X_0: Compute pre-initial guess 
  !
  if(ngrou/=0) then

     call solver_SpMV(solve_sol(1),an,xx,rr)
          
     do ii= 1, nrows
        rr(ii) = bb(ii) - rr(ii)                           ! r_{-1} = b-A.x_{-1}
     end do 
    
     call wtvect(npopo,ngrou,nbvar,mu,rr)                        ! W^T.r_{-1}
     if(kfl_paral/=0) then
        call LUsolv(&                                      ! A'.mu = W^T.r_{-1} 
             solve_sol(1)%ngrou*nbvar,solve_sol(1)%nskyl,&
             solve_sol(1)%iskyl,1_ip,askyl,&
             mu,solve_sol(1)%ngrou*nbvar,solve_sol(1)%idiag,info)
        if(info/=0) call runend('MAGRU2: COULD NOT SOLVE INITIAL SYSTEM')
     end if
     call wvect(npopo,nbvar,mu,rr)                               ! W.mu
     do ii=1,nrows
        xx(ii) = xx(ii) + rr(ii)                           ! x0 = x_{-1} + W.mu
     end do
  end if
  !wopos(1)='SOLUT'
  !wopos(2)='SCALA'
  !call postpr(xx,wopos,1_ip,kaka)


  !
  ! raux = <b',b'> and some computations needed for Diag. Scal.
  !
  if(idprecon==SOL_SQUARE.or.idprecon==SOL_DIAGONAL) then
       !
     ! PP = M^{-1} b
     !
     do ii= 1, nrows
        pp(ii) = wa2(ii) * bb(ii)
     end do

     call solver_parallel_vector_L2norm(solve_sol(1),pp,raux)

  else if(idprecon==SOL_MATRIX) then

     call solver_SpMV(solve_sol(1),pn,bb,pp)

  else

     call solver_parallel_vector_L2norm(solve_sol(1),bb,raux)

  end if
  !
  ! ||M^{-1} b||=0: Trivial solution XX=0
  !
  if(raux==0.0_rp) then
     do ii= 1, nrows
        xx(ii) = 0.0_rp
     end do
     goto 10
  end if

!
  ! Stop criterion STOPCRI = ||L^{-1}b||*eps
  !
  stopcri = eps * raux
  invnb   = 1.0_rp / raux
  kk      = 0
  !
  ! Initial residual: r0 = L^{-1}( b - A x^0 )
  !
  call solver_SpMV(solve_sol(1),an,xx,rr)
  do ii= 1, nrows
     rr(ii) = bb(ii) - rr(ii)
  end do
  if(idprecon==SOL_SQUARE.or.idprecon==SOL_DIAGONAL) then
     do ii= 1, nrows
        rr(ii) = rr(ii) * wa2(ii)
     end do
  end if
  !
  ! rho^0 = <r0,r0> = ||r0||^2
  !
  call solver_parallel_vector_L2norm(solve_sol(1),rr,resid)
  newrho = resid*resid
  resin  = resid*invnb

print*,'resid=',resid,stopcri
  if(resid<=stopcri) goto 10
!
  ! Initial x0 = R x
  !
  if(idprecon==SOL_SQUARE) then
     do ii= 1, nrows
        xx(ii) = wa1(ii) * xx(ii)
     end do
  end if

  !
  ! Initial residual r0 = L^{-1}( b - A x^0 )
  !
  do ii= 1, nrows
     r0(ii) = rr(ii)
  end do

  !-----------------------------------------------------------------
  !
  !  MAIN LOOP
  !
  !-----------------------------------------------------------------

  do while( kk<maxiter .and. resid>stopcri )
    if(kk==0) then       
        !
        ! Initial p^{k+1} = z^k - W.mu^k
        !
          do ii= 1, nrows
              pp(ii) = rr(ii)
           end do

     else 
        !
        ! beta = (rho^k/rho^{k-1})*(alpha/omega)
        !
        if((rho/=0.0_rp).and.(omega/=0.0_rp)) then
           beta = (newrho/rho) * (alpha/omega)
        else
           ierr = 1
           goto 10
        end if
        !
        ! p^{k+1} = r^{k} + beta*(p^k - omega*q^k)
        !
        do ii= 1, nrows
           pp(ii) = rr(ii) + beta * (pp(ii) - omega * qq(ii))
        end do
     end if
     !
     ! Solve A'.mu = W^T.A.p
     !
     if(ngrou/=0) then

      if(idprecon==SOL_NO_PRECOND) then

         call solver_SpMV(solve_sol(1),an,pp,wa1)

     else if(idprecon==SOL_SQUARE) then

        do ii= 1, nrows
           wa1(ii) = wa2(ii) * pp(ii)                           ! R^{-1} p^{k+1})
        end do
        call solver_SpMV(solve_sol(1),an,wa1,qq)                  ! A (...)
        do ii= 1, nrows
           wa1(ii) = wa2(ii) * qq(ii)                           ! L^{-1} [ ... ]
        end do

     else if(idprecon==SOL_DIAGONAL) then

        call solver_SpMV(solve_sol(1),an,pp,wa1)                   ! A p^{k+1}
        do ii= 1, nrows
           wa1(ii) = wa2(ii) * wa1(ii)                           ! L^{-1} [ ... ]
        end do

     endif 

        call wtvect(npopo,ngrou,nbvar,mu,wa1)
        if(kfl_paral/=0) then
           call LUsolv(&                                       ! A'.mu = W^T.A.p
                solve_sol(1)%ngrou*nbvar,solve_sol(1)%nskyl,solve_sol(1)%iskyl,&
                1_ip,askyl,mu,solve_sol(1)%ngrou*nbvar,solve_sol(1)%idiag,info)
           if(info/=0) call runend('MATGRO: COULD NOT SOLVE SYSTEM')
        end if
        call wvect(npopo,nbvar,mu,wa1)
        !
        !    p=p-Z.(Z^T.A.Z)^{-1}.Z^T.A.p
        !
        do ii= 1, nrows
           pp(ii) = pp(ii) -wa1(ii)
        end do
     end if
     !
     ! q^{k+1} = L^{-1} [ A (R^{-1} p^{k+1}) ]
     !
     if(idprecon==SOL_NO_PRECOND) then

        call solver_SpMV(solve_sol(1),an,pp,qq)                  

     else if(idprecon==SOL_SQUARE) then

        do ii= 1, nrows
           wa1(ii) = wa2(ii) * pp(ii)                           ! R^{-1} p^{k+1})
        end do
        call solver_SpMV(solve_sol(1),an,wa1,qq)                  ! A (...)
        do ii= 1, nrows
           qq(ii) = wa2(ii) * qq(ii)                            ! L^{-1} [ ... ]
        end do

     else if(idprecon==SOL_DIAGONAL) then

        call solver_SpMV(solve_sol(1),an,pp,qq)                   ! A p^{k+1}
        do ii= 1, nrows
           qq(ii) = wa2(ii) * qq(ii)                            ! L^{-1} [ ... ]
        end do

     else if(idprecon==SOL_MATRIX) then

        call solver_SpMV(solve_sol(1),an,pp,wa1)                  ! A p^{k+1}
        call solver_SpMV(solve_sol(1),pn,wa1,qq)                  ! P^{-1} [ ... ]

     end if
     !
     ! alpha = rho^k / <r0,q^{k+1}>
     !
     call solver_parallel_scalar_product(solve_sol(1),r0,qq,alpha)

     if (alpha==0.0_rp) then
        ierr = 2
        goto 10
     end if
     alpha = newrho / alpha
     !
     ! s = r^k - alpha*q^{k+1}
     !
     do ii= 1, nrows
        ss(ii) = rr(ii) - alpha * qq(ii)
     end do
     !
     ! Solve A'.mu = W^T.A.s
     !
     if(ngrou/=0) then 
       if(idprecon==SOL_NO_PRECOND) then

          call solver_SpMV(solve_sol(1),an,ss,wa1) 

     else if(idprecon==SOL_SQUARE) then

        do ii= 1, nrows
           wa1(ii) = wa2(ii) * ss(ii)                           ! R^{-1} p^{k+1})
        end do
        call solver_SpMV(solve_sol(1),an,wa1,qq)                  ! A (...)
        do ii= 1, nrows
           wa1(ii) = wa2(ii) * qq(ii)                            ! L^{-1} [ ... ]
        end do

     else if(idprecon==SOL_DIAGONAL) then

        call solver_SpMV(solve_sol(1),an,ss,wa1)                   ! A p^{k+1}
        do ii= 1, nrows
           wa1(ii) = wa2(ii) * wa1(ii)                            ! L^{-1} [ ... ]
        end do

     endif

        call wtvect(npopo,ngrou,nbvar,mu,wa1)
        if(kfl_paral/=0) then
           call LUsolv(&                                       ! A'.mu = W^T.A.s
                solve_sol(1)%ngrou*nbvar,solve_sol(1)%nskyl,solve_sol(1)%iskyl,&
                1_ip,askyl,mu,solve_sol(1)%ngrou*nbvar,solve_sol(1)%idiag,info)
           if(info/=0) call runend('MATGRO: COULD NOT SOLVE SYSTEM')
        end if
        call wvect(npopo,nbvar,mu,wa1)
!
        do ii= 1, nrows
           ss(ii) = ss(ii) - wa1(ii)
        end do
! 
     end if
!
     ! t = L^{-1} [ A (R^{-1} s) ]
     !
     if(idprecon==SOL_NO_PRECOND) then

        call solver_SpMV(solve_sol(1),an,ss,tt) 

     else if(idprecon==SOL_SQUARE) then

        do ii= 1, nrows
           wa1(ii) = wa2(ii) * ss(ii)
        end do
        call solver_SpMV(solve_sol(1),an,wa1,tt) 
        do ii= 1, nrows
           tt(ii) = wa2(ii) * tt(ii)
        end do

     else if(idprecon==SOL_DIAGONAL) then

        call solver_SpMV(solve_sol(1),an,ss,tt) 
        do ii= 1, nrows
           tt(ii) = wa2(ii) * tt(ii)
        end do

     else if(idprecon==SOL_MATRIX) then

        call solver_SpMV(solve_sol(1),an,ss,wa1) 
        call solver_SpMV(solve_sol(1),pn,wa1,tt) 

     end if
          !
     ! omega = <t,s> / <t,t>
     !
     call solver_parallel_scalar_product(solve_sol(1),tt,tt,ss,raux,omega)

     if (raux==0.0_rp) then
        ierr = 2
        goto 10
     end if
     omega = omega / raux
     !
     ! x^{k+1} = x^k + alpha*p^{k+1} + omega*s
     !
     do ii= 1, nrows
        xx(ii)  = xx(ii) + alpha * pp(ii) + omega * ss(ii)
     end do
     !
     ! r^{k+1} = s - omega*t
     !
     do ii= 1, nrows
        rr(ii) = rr(ii)-alpha *qq(ii) - omega * tt(ii)
     end do
     !
     ! rho^k = <r0,r^k> and || r^{k+1} ||
     !
     rho   = newrho
     call solver_parallel_scalar_product(solve_sol(1),rr,rr,r0,resid,newrho)
     resid = sqrt(resid)
  
     resi2 = resi1
     resi1 = resid*invnb
     kk    = kk+1
     !
     ! Convergence post process
     !
     if(kfl_cvgso_sol==1) write(lun_cvgso,100) kk,resi1

  end do

  !-----------------------------------------------------------------
  !
  !  END MAIN LOOP
  !
  !-----------------------------------------------------------------

  !
  ! x = [M]^-1 x'
  !
  if(idprecon==SOL_SQUARE) then
     do ii= 1, nrows
        xx(ii) = wa2(ii) * xx(ii)
     end do
  end if

10 resfi = resi1 
  solve_sol(1) % iters = kk

  if(kfl_solve_sol==1) then
     if(ierr/=0) write(lun_outso,120) kk
  end if

  if(ngrou/=0) then
     call memchk(2_ip,istat,memit,'ASKYL','defbcg',askyl)
     deallocate(askyl,stat=istat)
     if(istat/=0) call memerr(2_ip,'ASKYL','defbcg',0_ip)
     call memchk(2_ip,istat,memit,'MU','defbcg',mu)
     deallocate(mu,stat=istat)
     if(istat/=0) call memerr(2_ip,'MU','defbcg',0_ip)
  end if

  call memchk(2_ip,istat,memit,'WA2','defbcg',wa2)
  deallocate(wa2,stat=istat)
  if(istat/=0) call memerr(2_ip,'WA2','defbcg',0_ip)

  call memchk(2_ip,istat,memit,'WA1','defbcg',wa1)
  deallocate(wa1,stat=istat)
  if(istat/=0) call memerr(2_ip,'WA1','defbcg',0_ip)

  call memchk(2_ip,istat,memit,'SS','defbcg',ss)
  deallocate(ss,stat=istat)
  if(istat/=0) call memerr(2_ip,'SS','defbcg',0_ip)

  call memchk(2_ip,istat,memit,'QQ','defbcg',qq)
  deallocate(qq,stat=istat)
  if(istat/=0) call memerr(2_ip,'QQ','defbcg',0_ip)
  
  call memchk(2_ip,istat,memit,'TT','defbcg',tt)
  deallocate(tt,stat=istat)
  if(istat/=0) call memerr(2_ip,'TT','defbcg',0_ip)

  call memchk(2_ip,istat,memit,'PP','defbcg',pp)
  deallocate(pp,stat=istat)
  if(istat/=0) call memerr(2_ip,'PP','defbcg',0_ip)
 
  call memchk(2_ip,istat,memit,'R0','defbcg',r0)
  deallocate(r0,stat=istat)
  if(istat/=0) call memerr(2_ip,'R0','defbcg',0_ip)
 
  call memchk(2_ip,istat,memit,'RR','defbcg',rr)
  deallocate(rr,stat=istat)
  if(istat/=0) call memerr(2_ip,'RR','defbcg',0_ip)


100 format(i7,1x,e12.6)
110 format(i5,18(2x,e12.6))
120 format(&
       & '# Error at iteration ',i6,&
       & 'Dividing by zero: alpha = rho^k / <p^{k+1},q^{k+1}>')


end subroutine defbcg
