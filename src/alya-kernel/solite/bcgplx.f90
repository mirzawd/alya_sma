!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine bcgplx(npoin,nbvar,maxiter,precond,eps,an,ja,ia,bb,xx)

  !-----------------------------------------------------------------------------------------------------------
  ! Sources/kernel/solite/bcgplx.f90
  ! NAME
  !    bcgplx
  ! DESCRIPTION
  !    This routine solves a right preconditioned complex sparse linear system of algebraic equations:
  !                                     [A] [M]^-1  [M] x = b
  !                                       [A']         x' = b'
  !    using the Complex Biconjugate Gradient Stabilized algorithm with a preconditioner [M].
  ! INPUT ARGUMENTS
  !    NPOIN ..... Number of nodes in the interior of a mesh
  !    NBVAR ..... Number of unknowns in each node
  !    MAXITER ... Maximal number of iterations
  !    PRECOND ... Preconditioner 
  !    EPS ....... Tolerance
  !    AN ........ Sparse complex matrix of the original system in BCSR (Blocked Compressed Sparse Row) format
  !    JA ........ Compressed Sparse format: index vector for column numbers
  !    IA ........ Compressed Sparse format: index vector for beginning of a row block
  !    BB ........ RHS complex vector of the original system
  !    XX ........ Vector of complex unknowns of the original system - 
  !                input argument is an initial guess passed to the routine from outside
  ! OUTPUT ARGUMENTS 
  !    XX ........ Vector of complex unknowns of the original system 
  !-----------------------------------------------------------------------------------------------------------

  use def_kintyp, only       :  ip,rp
  use def_solver, only       :  memit,SOL_NO_PRECOND,SOL_DIAGONAL,SOL_DEFLATED,solve_sol
  use def_master, only       :  IMASTER,INOTMASTER,INOTSLAVE
!  use def_master, only       :  kfl_paral
  use mod_memchk                               
  use mod_csrdirx
  use mod_precon

  !Declaration statements
  implicit none

  real(rp),    parameter      :: MAXVALUE = 1.0e30_rp       !Value used for divergence check

  !Dummy arguments
  integer(ip), intent(in)     :: npoin,nbvar,maxiter,precond
  real(rp),    intent(in)     :: eps
  complex(rp), intent(in)     :: an(*)
  integer(ip), intent(in)     :: ja(*),ia(*)
  complex(rp), intent(in)     :: bb(*)
  complex(rp), intent(in out) :: xx(*)

  !Local variables
  character(len=20)        :: filename
  logical                  :: FINISH 
  integer(ip)              :: ii,ierr,ierror,info
  integer(ip)              :: nbnodes,totnode
  integer(4)               :: istat
  integer(ip)              :: ngrou,ncors
  real(rp)                 :: raub,raux,stopcri,resid,omegar
!  real(rp)                 :: cpu_syma1,cpu_syma2,cpu_syma
!  real(rp)                 :: cpu_smal1,cpu_smal2,cpu_smal
  complex(rp)              :: alpha,beta,omega,omega1,rho,newrho
  complex(rp), pointer     :: rr(:),rrT(:),pp(:),ppH(:),qq(:),ss(:),ssH(:),tt(:)
  complex(rp), pointer     :: vv(:),diag(:),invdiag(:)

  integer(ip), pointer     :: iL(:),jL(:)
  complex(rp), pointer     :: Ln(:)
  integer(ip), pointer     :: iU(:),jU(:)
  complex(rp), pointer     :: Un(:)
  complex(rp), pointer     :: pc(:),qc(:),ancors(:)
  integer(ip), pointer     :: invpR(:),invpC(:)
  integer(ip), pointer     :: iagro(:),jagro(:)

  !Information needed for deflation
  ngrou =  solve_sol(1) % ngrou                        !Number of groups 
  !write(*,*)'ngrou',solve_sol(1) % ngrou
  !write(*,*)'precond',precond
  !write(*,*)'csr',solve_sol(1) % kfl_defas
  ncors = (solve_sol(1) % nzgro) * nbvar * nbvar       !Number of elements of the deflated matrix     
  iagro => solve_sol(1) % iagro                        !Index vector for beginning of a row block for the deflated matrix
  jagro => solve_sol(1) % jagro                        !Index vector for column numbers for the deflated matrix

  if (precond == SOL_DEFLATED) then
     if (ngrou == 0) call runend('Error! Number of groups for deflated preconditioning must be larger than zero.')
  endif

  !Information needed for memory allocation for working arrays
  if (IMASTER) then
     nbnodes = 0_ip       !Master will allocate minimum memory for working arrays since it does not perform any calculations
     totnode = 1_ip
  else
     nbnodes = npoin
     totnode = npoin * nbvar
  endif

  !Compute Euclidean norm of b: raub = ||b||  
  call no2plx(nbnodes,nbvar,bb,raub)
  IF (raub == 0.0_rp) THEN
     !If ||b|| = 0 then x = 0: the trivial solution
     if (INOTMASTER) then
        do ii = 1,totnode
           xx(ii) = (0.0_rp,0.0_rp)
        enddo
     endif
     solve_sol(1) % iters =  0_ip
     ierr  = -1_ip
  ELSE
     !Compute the stop criterion: stopcri = eps * ||b|| 
     stopcri = eps * raub

     if (INOTSLAVE) then
        filename = 'resid-fwd.res'
        open(unit=11,file=filename,iostat=ierror)
        write(11,*) raub, eps
     endif

     !Allocate memory for working arrays
     allocate(rr(totnode), stat=istat)
     call memchk(0_ip,istat,memit,'RR','bcgplx',rr) 
     allocate(rrT(totnode),stat=istat)
     call memchk(0_ip,istat,memit,'RRT','bcgplx',rrT)
     allocate(pp(totnode),stat=istat)
     call memchk(0_ip,istat,memit,'PP','bcgplx',pp) 
     allocate(ppH(totnode),stat=istat)
     call memchk(0_ip,istat,memit,'PPH','bcgplx',ppH)   
     allocate(qq(totnode), stat=istat) 
     call memchk(0_ip,istat,memit,'QQ','bcgplx',qq) 
     allocate(ss(totnode), stat=istat)
     call memchk(0_ip,istat,memit,'SS','bcgplx',ss) 
     allocate(ssH(totnode),stat=istat)
     call memchk(0_ip,istat,memit,'SSH','bcgplx',ssH) 
     allocate(tt(totnode), stat=istat)
     call memchk(0_ip,istat,memit,'TT','bcgplx',tt) 
     allocate(vv(totnode), stat=istat)
     call memchk(0_ip,istat,memit,'VV','bcgplx',vv)
     allocate(diag(totnode),stat=istat) 
     call memchk(0_ip,istat,memit,'DIAG','bcgplx',diag)
     allocate(invdiag(totnode),stat=istat)
     call memchk(0_ip,istat,memit,'INVDIAG','bcgplx',invdiag)  

     if (IMASTER) totnode = 0_ip       !Master does not perform any loop

     !If Jacobi (diagonal) preconditioner is used
     if (precond == SOL_DIAGONAL .or. precond == SOL_DEFLATED) then

        call diagpr(nbnodes,nbvar,an,ja,ia,diag,invdiag)
     end if

     !If deflation is performed 

     IF (precond == SOL_DEFLATED) THEN  
        !Allocate memory for groups
        allocate(pc((ngrou+1)*nbvar),stat=istat)       !+1 because last dimension used in dot product
        call memchk(0_ip,istat,memit,'PC','bcgplx',pc)

        allocate(qc((ngrou+1)*nbvar),stat=istat)         
        call memchk(0_ip,istat,memit,'QC','bcgplx',qc)

        allocate(ancors(ncors),stat=istat)
        call memchk(0_ip,istat,memit,'ANCORS','bcgplx',ancors)

        !!call cputim(cpu_syma1)
        !Compute coarse matrix: A'= W^(T) * A * W 
        call matgrx2(ngrou,nbnodes,ncors,nbvar,ia,ja,an,ancors)
        !call PAR_BARRIER() 
        !!call cputim(cpu_syma2)
        !!if (INOTSLAVE) then
        !!cpu_syma = cpu_syma2 - cpu_syma1
        !!write(*,*)'Time to build A:',cpu_syma, kfl_paral
        !!endif  
        !Factorize coarse matrix A' using LU factorization: A' = L * U
        !!call cputim(cpu_syma1)
        if (INOTMASTER) then
           nullify(iL,jL,Ln,iU,jU,Un)                   
           nullify(invpR,invpC)       !Permutation arrays
           call CSR_LU_Factorizationx(ngrou,nbvar,iagro,jagro,ancors,iL,jL,Ln,iU,jU,Un,info)
           if (info /= 0) call runend('DEFLATION: SINGULAR MATRIX')
        elseif (IMASTER) then
           allocate(iL(1),stat=istat)
           allocate(jL(1),stat=istat)
           allocate(Ln(1),stat=istat)
           allocate(iU(1),stat=istat)
           allocate(jU(1),stat=istat)
           allocate(Un(1),stat=istat)
        endif
        !call PAR_BARRIER() 
        !!call cputim(cpu_syma2)
        !!if (INOTSLAVE) then
        !!cpu_syma = cpu_syma2 - cpu_syma1
        !!write(*,*)'Time to factorize A:',cpu_syma,kfl_paral
        !!endif  
     ENDIF

     !Calculate initial residual r0   

     !Calculate Euclidean norm of x0: raux = ||x0||
     call no2plx(nbnodes,nbvar,xx,raux)
     if (raux == 0.0_rp) then
        !If x0 = 0 then r0 = b 
        do ii = 1,totnode
           rr(ii) = bb(ii)                               !r0 = b
        enddo
     else
        call runend('x0 is not zero')
        !If x0 /= 0 then r0 = b - A * x0 
        call bcsplx(nbnodes,nbvar,an,ja,ia,xx,rr)       !A * x0   
        do ii = 1,totnode
           rr(ii) = bb(ii) - rr(ii)                      !r0 = b - A * x0
        enddo
        !If x0 /= 0 then x0' = M * x0
        if (precond == SOL_DIAGONAL) then
!!$omp parallel do                 &
!!$omp            default(shared)  &
!!$omp            private(ii)      &
!!$omp            schedule(static) 
           do ii = 1,totnode
              xx(ii) = diag(ii) * xx(ii)                  !x0' = D * x0
           enddo
        elseif (precond == SOL_DEFLATED) then
           write(*,*)'x = M * x'                         !x0' = M * x0
        endif
     endif
     !Calculate Euclidean norm of r0: resid = ||r0|| (for right preconditioning: r'= r)
     call no2plx(nbnodes,nbvar,rr,resid) 

     !Test if the initial guess is the solution: if ||r0|| <= eps * ||b|| then x = x0 else calculate x1...
     if (resid <= stopcri) then
        FINISH = .true.
     else
        FINISH = .false.

        !r0tilde = r0
        ! !$omp parallel do                 &
        ! !$omp            default(shared)  &
        ! !$omp            private(ii)      &
        ! !$omp            schedule(static) 
        do ii = 1,totnode
           rrT(ii) = rr(ii)
        enddo
        ! !$omp end parallel do
        !Calculate inner product: newrho = <r0, r0tilde> = <r0, r0> = ||r0||^2
        newrho = cmplx(resid * resid,0.0_rp,kind=rp)
        !Initialize rho = alpha = omega = 1
        rho   = (1.0_rp,0.0_rp)
        alpha = (1.0_rp,0.0_rp)
        omega = (1.0_rp,0.0_rp)
     endif

     solve_sol(1) % iters = 0_ip
     ierr  = 0_ip

     !-------------------------------------------------
     !
     !       MAIN LOOP
     !
     !------------------------ -------------------------

     DO WHILE (.not. FINISH)
        !!call cputim(cpu_syma1)
        IF ((rho /= (0.0_rp,0.0_rp)).and.(omega /= (0.0_rp,0.0_rp))) THEN       !if rho /= 0 and omega /= 0
           if (solve_sol(1) % iters == 0_ip) then 
              !Initial p0 = r0
              ! !$omp parallel do                 &
              ! !$omp            default(shared)  &
              ! !$omp            private(ii)      &
              ! !$omp            schedule(static) 
              do ii = 1,totnode
                 pp(ii) = rr(ii)
              enddo
              ! !$omp end parallel do
           else
              !rho = <r(k-1), r0tilde>  
              rho = newrho
              !Calculate inner product: newrho = <r(k), r0tilde> 
              call proplx(nbnodes,nbvar,rr,rrT,newrho)
              !beta(k) = (<r(k), r0tilde> / <r(k-1), r0tilde>) * (alpha(k-1) / omega(k-1))
              beta = (newrho / rho) * (alpha / omega)       
              !p(k) = r(k) + beta(k) * (p(k-1) - omega(k-1) * q(k-1))
              ! !$omp parallel do                 &
              ! !$omp            default(shared)  &
              ! !$omp            private(ii)      &
              ! !$omp            schedule(static) 
              do ii = 1,totnode
                 pp(ii) = rr(ii) + beta * (pp(ii) - omega * qq(ii))
              enddo
              ! !$omp end parallel do
           endif

           !Calculate alpha: alpha(k) = <r(k), r0tilde> / <A * M^(-1) * p(k), r0tilde>

           if (precond == SOL_NO_PRECOND) then
              !q(k) = A * p(k) 
              call bcsplx(nbnodes,nbvar,an,ja,ia,pp,qq)

           else if (precond == SOL_DIAGONAL) then
              !phat(k) = D^(-1) * p(k)
              ! !$omp parallel do                 &
              ! !$omp            default(shared)  &
              ! !$omp            private(ii)      &
              ! !$omp            schedule(static) 
              do ii = 1,totnode
                 ppH(ii) = invdiag(ii) * pp(ii)
              enddo
              ! !$omp end parallel do
              !q(k) = A * phat(k)
              call bcsplx(nbnodes,nbvar,an,ja,ia,ppH,qq)

           else if (precond == SOL_DEFLATED) then 
              !phat(k) = M^(-1) * p(k)
              !call precon1(nbnodes,ngrou,nbvar,invpR,invpC,iL,jL,Ln,iU,jU,Un,pc,qc,pp,ppH)
              !!call cputim(cpu_smal1)
              call precon2(nbnodes,ngrou,nbvar,an,ja,ia,diag,invdiag,invpR,invpC,iL,jL,Ln,iU,jU,Un,pc,qc,pp,ppH)
              !call PAR_BARRIER() 
              !!call cputim(cpu_smal2)
              !!if (INOTSLAVE) then
              !!cpu_smal = cpu_smal2 - cpu_smal1
              !!write(*,*)'Time to deflate 1:',cpu_smal,kfl_paral
              !!endif
              !q(k) = A * phat(k)
              call bcsplx(nbnodes,nbvar,an,ja,ia,ppH,qq)

           else
              call runend('Error! Unknown preconditioner.')
           endif

           !alpha = <q(k), r0tilde>
           call proplx(nbnodes,nbvar,qq,rrT,alpha)       
           IF (alpha /= (0.0_rp,0.0_rp)) THEN            !if <q(k), r0tilde> /= 0
              alpha = newrho / alpha

              !Calculate s
              !s(k) = r(k) - alpha(k) * q(k)
              ! !$omp parallel do                 &
              ! !$omp            default(shared)  &
              ! !$omp            private(ii)      &
              ! !$omp            schedule(static) 
              do ii = 1,totnode
                 ss(ii) = rr(ii) - alpha * qq(ii)
              enddo
              ! !$omp end parallel do
              !Calculate omega:omega(k) = <s(k), A * M^(-1) * s(k)> / <A * M^(-1) * s(k), A * M^(-1) * s(k)>

              if (precond == SOL_NO_PRECOND) then
                 !t(k) = A * s(k) 
                 call bcsplx(nbnodes,nbvar,an,ja,ia,ss,tt) 

              elseif (precond == SOL_DIAGONAL) then
                 !shat(k) = D^(-1) * s(k)
                 ! !$omp parallel do                 &
                 ! !$omp            default(shared)  &
                 ! !$omp            private(ii)      &
                 ! !$omp            schedule(static) 
                 do ii = 1,totnode
                    ssH(ii) = invdiag(ii) * ss(ii)
                 enddo
                 ! !$omp end parallel do
                 !t(k) = A * shat(k)
                 call bcsplx(nbnodes,nbvar,an,ja,ia,ssH,tt)

              else if (precond == SOL_DEFLATED) then 
                 !shat(k) = M^(-1) * s(k)
                 !call precon1(nbnodes,ngrou,nbvar,invpR,invpC,iL,jL,Ln,iU,jU,Un,pc,qc,ss,ssH)
                 !!call cputim(cpu_smal1)
                 call precon2(nbnodes,ngrou,nbvar,an,ja,ia,diag,invdiag,invpR,invpC,iL,jL,Ln,iU,jU,Un,pc,qc,ss,ssH) 
                 !call PAR_BARRIER() 
                 !!call cputim(cpu_smal2)
                 !!if (INOTSLAVE) then
                 !!cpu_smal = cpu_smal2 - cpu_smal1
                 !!write(*,*)'Time to deflate 2:',cpu_smal,kfl_paral
                 !!endif       
                 !t(k) = A * shat(k)
                 call bcsplx(nbnodes,nbvar,an,ja,ia,ssH,tt)

              else
                 call runend('Error! Unknown preconditioner.')
              endif

              !omega = <t(k), t(k)>
              call no2plx(nbnodes,nbvar,tt,omegar)       
              omega = cmplx(omegar * omegar,0.0_rp,kind=rp)
              IF (omega /= (0.0_rp,0.0_rp)) THEN         !if <t(k), t(k)> /= 0
                 call proplx(nbnodes,nbvar,ss,tt,omega1)
                 omega = omega1 / omega

                 !Calculate new approximation of the solution x         
                 !x(k+1) = x(k) + alpha(k) * p(k) + omega(k) * s(k)
                 ! !$omp parallel do                 &
                 ! !$omp            default(shared)  &
                 ! !$omp            private(ii)      &
                 ! !$omp            schedule(static) 
                 do ii = 1,totnode
                    xx(ii) = xx(ii) + alpha * pp(ii) + omega * ss(ii)
                 enddo
                 ! !$omp end parallel do
                 !Calculate residual r
                 !r(k+1) = s(k) - omega(k) * t(k)
                 ! !$omp parallel do                 &
                 ! !$omp            default(shared)  &
                 ! !$omp            private(ii)      &
                 ! !$omp            schedule(static) 
                 do ii = 1,totnode
                    rr(ii) = ss(ii) - omega * tt(ii)
                 enddo
                 ! !$omp end parallel do
                 !Calculate Euclidean norm of r(k+1): resid = ||r(k+1)||
                 call no2plx(nbnodes,nbvar,rr,resid) 
                 if (INOTSLAVE) write(11,*) solve_sol(1) % iters+1,resid

                 if (INOTSLAVE) write(solve_sol(1) % lun_cvgso,100) solve_sol(1) % iters,resid

                 if (resid <= stopcri) then           !Iterative method converges
                    FINISH = .true.
                 elseif (resid > MAXVALUE) then       !Iterative method diverges
                    FINISH = .true.
                    ierr = 3_ip
                 else
                    if (solve_sol(1) % iters >= maxiter) then         !Iterative method does not convege in maxiter iterations
                       FINISH = .true.
                       ierr = 1_ip
                    endif
                 endif
                 solve_sol(1) % iters = solve_sol(1) % iters + 1_ip
              ELSE                                !if <t(k), t(k)> == 0
                 FINISH = .true.
                 ierr = 2_ip
              ENDIF
           ELSE                                 !if <q(k), r0tilde> == 0 
              FINISH = .true.
              ierr = 2_ip
           ENDIF
        ELSE                                  !if rho == 0 or omega == 0
           FINISH = .true.
           ierr = 2_ip
        ENDIF
        !call PAR_BARRIER() 
        !!call cputim(cpu_syma2)
        !!if (INOTSLAVE) then
        !!cpu_syma = cpu_syma2 - cpu_syma1
        !!write(*,*)'Time of one iteration:',cpu_syma,kfl_paral
        !!endif

     ENDDO
     !-------------------------------------------------
     !
     !                   END MAIN LOOP
     !
     !-------------------------------------------------

     !x = D^(-1) * x'
     if (precond == SOL_DIAGONAL) then
        ! !$omp parallel do                 &
        ! !$omp            default(shared)  &
        ! !$omp            private(ii)      &
        ! !$omp            schedule(static) 
        do ii = 1,totnode
           xx(ii) = invdiag(ii) * xx(ii)        
        enddo
        ! !$omp end parallel do
        !x = M^(-1) * x'
     elseif (precond == SOL_DEFLATED) then 
        !call precon1(nbnodes,ngrou,nbvar,invpR,invpC,iL,jL,Ln,iU,jU,Un,pc,qc,xx,vv)
        !!call cputim(cpu_smal1)
        call precon2(nbnodes,ngrou,nbvar,an,ja,ia,diag,invdiag,invpR,invpC,iL,jL,Ln,iU,jU,Un,pc,qc,xx,vv)
        !call PAR_BARRIER() 
        !!call cputim(cpu_smal2)
        !!if (INOTSLAVE) then
        !!cpu_smal = cpu_smal2 - cpu_smal1
        !!write(*,*)'Time to deflate 3:',cpu_smal,kfl_paral
        !!endif
        ! !$omp parallel do                 &
        ! !$omp            default(shared)  &
        ! !$omp            private(ii)      &
        ! !$omp            schedule(static) 
        do ii = 1,totnode
           xx(ii) = vv(ii)       
        enddo
        ! !$omp end parallel do
     endif
     !write(*,*)'reziduo 1', resid
     call no2plx(nbnodes,nbvar,xx,resid) 
     !write(*,*)'resenje', resid
     call bcsplx(nbnodes,nbvar,an,ja,ia,xx,vv)       !A * x0   
     ! !$omp parallel do                 &
     ! !$omp            default(shared)  &
     ! !$omp            private(ii)      &
     ! !$omp            schedule(static) 
     do ii = 1,totnode
        vv(ii) = bb(ii) - vv(ii)                      !r0 = b - A * x0
     enddo
     ! !$omp end parallel do
     call no2plx(nbnodes,nbvar,vv,resid) 
     !write(*,*)'reziduo 2', resid
     if (INOTSLAVE) close (unit=11)

     !Deallocate memory for groups   
     IF (precond == SOL_DEFLATED) THEN
        call memchk(2_ip,istat,memit,'ANCORS','bcgplx',ancors)
        deallocate(ancors,stat=istat)
        if(istat/=0) call memerr(2_ip,'ANCORS','bcgplx',0_ip)

        call memchk(2_ip,istat,memit,'QC','bcgplx',qc)
        deallocate(qc,stat=istat)
        if (istat /= 0) call memerr(2_ip,'QC','bcgplx',0_ip)

        call memchk(2_ip,istat,memit,'PC','bcgplx',pc)
        deallocate(pc,stat=istat)
        if(istat/=0) call memerr(2_ip,'PC','bcgplx',0_ip)

        if (INOTMASTER) then         
           call CSR_LUfinx(IL,JL,LN,IU,JU,UN)      
        endif
     ENDIF
     !Deallocate memory of working arrays
     call memchk(2_ip,istat,memit,'INVDIAG','bcgplx',invdiag)
     deallocate(invdiag,stat=istat)
     if ( istat /= 0_ip ) call memerr(2_ip,'INVDIAG','bcgplx',0_ip)

     call memchk(2_ip,istat,memit,'DIAG','bcgplx',diag)
     deallocate(diag,stat=istat)
     if ( istat /= 0_ip ) call memerr(2_ip,'DIAG','bcgplx',0_ip)

     call memchk(2_ip,istat,memit,'VV','bcgplx',vv)
     deallocate(vv,stat=istat)
     if ( istat /= 0_ip ) call memerr(2_ip,'VV','bcgplx',0_ip)

     call memchk(2_ip,istat,memit,'TT','bcgplx',tt)
     deallocate(tt,stat=istat)
     if ( istat /= 0_ip ) call memerr(2_ip,'TT','bcgplx',0_ip)

     call memchk(2_ip,istat,memit,'SSH','bcgplx',ssH)
     deallocate(ssH,stat=istat)
     if ( istat /= 0_ip ) call memerr(2_ip,'SSH','bcgplx',0_ip)

     call memchk(2_ip,istat,memit,'SS','bcgplx',ss)
     deallocate(ss,stat=istat)
     if ( istat /= 0_ip ) call memerr(2_ip,'SS','bcgplx',0_ip)

     call memchk(2_ip,istat,memit,'QQ','bcgplx',qq)
     deallocate(qq,stat=istat)
     if ( istat /= 0_ip ) call memerr(2_ip,'QQ','bcgplx',0_ip)

     call memchk(2_ip,istat,memit,'PPH','bcgplx',ppH)
     deallocate(ppH,stat=istat)
     if ( istat /= 0_ip ) call memerr(2_ip,'PPH','bcgplx',0_ip)

     call memchk(2_ip,istat,memit,'PP','bcgplx',pp)
     deallocate(pp,stat=istat)
     if ( istat /= 0_ip ) call memerr(2_ip,'PP','bcgplx',0_ip)

     call memchk(2_ip,istat,memit,'RRT','bcgplx',rrT)
     deallocate(rrT,stat=istat)
     if ( istat /= 0_ip ) call memerr(2_ip,'RRT','bcgplx',0_ip)

     call memchk(2_ip,istat,memit,'RR','bcgplx',rr)
     deallocate(rr,stat=istat)
     if ( istat /= 0_ip ) call memerr(2_ip,'RR','bcgplx',0_ip)
  ENDIF

  !if (INOTSLAVE) then
  !  write (*,*) ' ', ierr, ' '
  !  write (*,*) ' ', solve_sol(1) % iters, ' '
  !endif

100 format(i7,2(1x,e12.6))

end subroutine bcgplx

