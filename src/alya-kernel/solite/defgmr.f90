!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine defgmr( &
     nbnodes, nbvar,&
     idprecon, maxiter, kryldim, &
     kfl_cvgso_sol, lun_cvgso,&
     kfl_solve_sol, lun_outso,&
     eps, an, pn, ja, ia, bb, xx )
  !-----------------------------------------------------------------------
  !
  ! Objective: Solve    [L]^-1 [A] [R]^-1  [R] xx  = [L]^-1 bb
  !                            [A']            xx' =        bb'
  ! by the GMRES(m) method.
  !
  !            Three preconditioners are possible:
  !            idprecon = 0 => [L] = I, [R] = I
  !            idprecon = 1 => [L] = [R] = sqrt(abs(diag([A])))
  !            idprecon = 2 => [L] = diag([A]), [R] = I
  !
  !
  !            If Diag. Scaling is selected the preconditioned system is:
  !                    D^-1/2 [A] D^-1/2   D^1/2 x  =  D^-1/2 b
  !                          [A']              x'   =      b'
  !
  !
  !  kryldim : Dimension of the Krylov basis
  !  kryl    : Working matrix with the Krylov basis.
  !            It must be allocated with:  nbnodes*nbvar*(kryldim+1)
  !  hh      : [H]. It is stored by columns. Its dimension is:
  !            (kryldim+1) * kryldim. It must be allocated with:
  !            SUM (i=1, kryldim) {ii + 1} =
  !                                    kryldim + (kryldim * kryldim+1) / 2
  !  cc, ss  : Givens rotations. They must be allocated with:  kryldim
  !  rs      : Right Hand side of the LSQ problem.
  !            It must be allocated with: kryldim+1
  !
  !  solve_sol(1) % iters : present iteration
  !  jj    : present iteration relative to [1 ... kryldim]
  !
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp,lg
  use def_master, only     :  kfl_paral,intost
!  use def_master, only     :  INOTMASTER
  use def_solver, only     :  SOL_NO_PRECOND,SOL_SQUARE,SOL_DIAGONAL,SOL_MATRIX,resi1,solve_sol
!  use def_solver, only     :  memit,SOL_NO_PRECOND
!  use def_solver, only     :  resin,resfi,resi1,resi2,solve_sol
  use mod_memchk
  use mod_postpr
  implicit none
  integer(ip), intent(in)  :: nbnodes, nbvar, idprecon, maxiter, kryldim
  integer(ip), intent(in)  :: kfl_cvgso_sol, lun_cvgso, kfl_solve_sol, lun_outso
  real(rp),    intent(in)  :: eps
  real(rp),    intent(in)  :: an(nbvar, nbvar, *)
  real(rp),    intent(in)  :: pn(nbvar, nbvar, *)
  integer(ip), intent(in)  :: ja(*), ia(*)
  real(rp),    intent(in)  :: bb(*)
  real(rp),    intent(out) :: xx(nbvar*nbnodes)
  logical(lg)              :: convergence, fin2
  integer(ip)              :: nrows, npopo,ngrou,nskyl
!  integer(ip)              :: ii, jj, kk, ll, jj1, idx, nrows, npopo,ngrou,nskyl,info,ibvar
!  integer(ip)              :: igrou,ipoin
!  integer(4)               :: istat
!  real(rp)                 :: raux, gamma, stopcri, invnb
  real(rp),    parameter   :: epsmac = 1.0d-16

!  real(rp),    allocatable :: wa1(:), wa3(:), wa4(:), kryl(:,:)
!  real(rp),    pointer     :: wa2(:),mu(:),askyl(:)
!  real(rp),    allocatable :: cc(:), ss(:), rs(:), hh(:)
!  character(20)            :: mess1,mess2
!  character(150)           :: messa
  integer(ip),pointer          :: lgrou(:)
!  character(5)              :: wopos(2) ! OJO
!  real(rp) :: kaka=0.0_rp ! OJO
#ifdef EVENT
  call mpitrace_user_function(1)
#endif
print*,'popo'
  ngrou =  solve_sol(1)%ngrou
  lgrou => solve_sol(1)%lgrou
  nskyl =  solve_sol(1)%nskyl
  solve_sol(1) % iters       = 0
  convergence = .false.
  fin2        = .false.
  resi1       = 1.0_rp
  if(kfl_paral==0) then
     nrows   = 1                ! Minimum memory for working arrays
     npopo   = 0                ! master does not perform any loop
  else
     npopo   = nbnodes
     nrows   = nbnodes * nbvar
  end if
!!$  !
!!$  ! Allocate memory of arrays with Krylov dimension dependent
!!$  !
!!$  allocate(cc(kryldim),stat=istat)
!!$  call memchk(0_ip,istat,memit,'CC','gmrpls',cc)
!!$  allocate(ss(kryldim),stat=istat)
!!$  call memchk(0_ip,istat,memit,'SS','gmrpls',ss)
!!$  allocate(rs(kryldim+1),stat=istat)
!!$  call memchk(0_ip,istat,memit,'RS','gmrpls',rs)
!!$  allocate(hh(kryldim+((kryldim*(kryldim+1))/2)),stat=istat)
!!$  call memchk(0_ip,istat,memit,'HH','gmrpls',hh)
!!$  !
!!$  ! Allocate memory for large arrays
!!$  !
!!$  allocate(wa1(nrows),stat=istat)
!!$  call memchk(0_ip,istat,memit,'WA1','gmrpls',wa1)
!!$  allocate(wa2(nrows),stat=istat)
!!$  call memchk(0_ip,istat,memit,'WA2','gmrpls',wa2)
!!$  allocate(wa3(nrows),stat=istat)
!!$  call memchk(0_ip,istat,memit,'WA3','gmrpls',wa3)
!!$  allocate(kryl(nrows,kryldim+1),stat=istat)
!!$  call memchk(0_ip,istat,memit,'KRYL','gmrpls',kryl)
!!$  if(idprecon>=SOL_MATRIX) then
!!$     allocate(wa4(nrows),stat=istat)
!!$     call memchk(0_ip,istat,memit,'WA4','gmrpls',wa4)
!!$  end if
!!$
!!$  if(ngrou/=0) then
!!$     allocate(mu(ngrou*nbvar),stat=istat)
!!$     call memchk(0_ip,istat,memit,'MU','defgmr',mu)
!!$     allocate(askyl(nskyl),stat=istat)
!!$     call memchk(0_ip,istat,memit,'ASKYL','defgmr',askyl)
!!$  end if
!!$
!!$  if(kfl_paral==0) nrows=0 ! Master does not perform any loop
!!$  !
!!$  ! Solve A'mu= W^T.b => x = W.mu
!!$  !
!!$print*,'a'
!!$  call matgru(ngrou,npopo,nskyl,nbvar,ia,ja,an,askyl,wa1)
!!$print*,'b'
!!$  if(ngrou/=0) then
!!$     call wtvect(npopo,ngrou,nbvar,mu,bb)                   ! W^T.b
!!$print*,'c'
!!$     if( INOTMASTER ) then
!!$        call LUsolv(&                                       ! A'.mu = W^T.r_{-1}
!!$             solve_sol(1)%ngrou*nbvar,solve_sol(1)%nskyl,&
!!$             solve_sol(1)%iskyl,1_ip,askyl,&
!!$             mu,solve_sol(1)%ngrou*nbvar,solve_sol(1)%idiag,info)
!!$        if( info /= 0 ) call runend('MATGRU: COULD NOT SOLVE INITIAL SYSTEM')
!!$     end if
!!$print*,'d'
!!$     call wvect2(npopo,nbvar,mu,xx)   
!!$print*,'e'
!!$print*,mu
!!$     !do ipoin = 1,npopo
!!$     !   igrou = solve_sol(1)%lgrou(ipoin)
!!$     !   jj    = (igrou-1)*nbvar
!!$     !   ii    = (ipoin-1)*nbvar
!!$     !   do ibvar = 1,nbvar
!!$     !      jj = jj + 1
!!$     !      ii = ii + 1
!!$     !      if( igrou > 0 ) then
!!$     !         xx(ii) = mu(jj)
!!$     !      end if
!!$     !   enddo
!!$     !end do
!!$
!!$    call memchk(2_ip,istat,memit,'ASKYL','gmrpls',askyl)
!!$    deallocate(askyl,stat=istat)
!!$    if(istat/=0) call memerr(2_ip,'ASKYL','gmrpls',0_ip)
!!$
!!$    call memchk(2_ip,istat,memit,'MU','gmrpls',mu)
!!$    deallocate(mu,stat=istat)
!!$    if(istat/=0) call memerr(2_ip,'MU','gmrpls',0_ip)
!!$
!!$  end if
!!$  return
!!$
!!$  !
!!$  ! Compute diagonal DIAG (and exchange in parallel)
!!$  !
!!$  if(idprecon==SOL_SQUARE.or.idprecon==SOL_DIAGONAL) then
!!$
!!$     if(nbvar==1) then
!!$        do ii= 1, npopo 
!!$           jj = ia(ii)
!!$           ll = -1
!!$           do while (jj< ia(ii+1) .and. ll ==-1)
!!$              if(ja(jj)==ii) then
!!$                 ll = jj
!!$              end if
!!$              jj = jj+1
!!$           end do
!!$           if(ll/=-1) wa2(ii)=an(1,1,ll)
!!$        end do
!!$     else
!!$        do ii= 1, npopo 
!!$           jj = ia(ii)
!!$           ll = -1
!!$           do while (jj< ia(ii+1) .and. ll ==-1)
!!$              if(ja(jj)==ii) then
!!$                 ll = jj
!!$              end if
!!$              jj = jj+1
!!$           end do
!!$           if(ll/=-1) then
!!$              jj = (ii-1) * nbvar
!!$              do kk= 1, nbvar
!!$                 wa2(jj+kk)=an(kk,kk,ll)
!!$              end do
!!$           end if
!!$        end do
!!$     end if
!!$     !
!!$     ! Periodicity: RHS Master-slave strategy 
!!$     !
!!$     call rhsmod(nbvar,wa2) 
!!$
!!$  end if
!!$  !
!!$  ! Initial residual
!!$  !
!!$  call resvec(nbvar,npopo,bb,wa1,resin)
!!$  !
!!$  ! raux = <b',b'> and some computations needed for Diag. Scal.
!!$  !
!!$  if(idprecon==SOL_SQUARE) then
!!$     !
!!$     ! wa1 = D^1/2    wa2 = D^-1/2
!!$     !
!!$     do ii= 1, npopo
!!$        jj = ia(ii)
!!$        ll = -1
!!$        do while (jj< ia(ii+1) .and. ll ==-1)
!!$           if(ja(jj)==ii) then
!!$              ll = jj
!!$           end if
!!$           jj = jj+1
!!$        end do
!!$
!!$        if(ll/=-1) then
!!$           jj = (ii-1) * nbvar
!!$
!!$           do kk= 1, nbvar
!!$              if(wa2(jj+kk)/=0.0_rp) then
!!$                 wa1(jj+kk) = sqrt( abs(wa2(jj+kk)) )
!!$                 wa2(jj+kk) = 1.0_rp / wa1(jj+kk)
!!$              else
!!$                 mess1=intost(ii)
!!$                 mess2=intost(kk)
!!$                 messa='GMRES: ZERO DIAGONAL FOUND AT (NODE,DOF)=('&
!!$                      //trim(mess1)//','//trim(mess2)//')'                    
!!$                 if(kfl_solve_sol/=0)&
!!$                      write(lun_outso,*) trim(messa)
!!$                 call runend(messa)
!!$              end if
!!$           end do
!!$        else
!!$           mess1=intost(ii)
!!$           messa='GMRES: NO DIAGONAL FOUND FOR NODE='//trim(mess1)              
!!$           if(kfl_solve_sol/=0)&
!!$                write(lun_outso,*) trim(messa)
!!$           call runend(messa)                   
!!$        end if
!!$     end do
!!$     !
!!$     ! b' = D^-1/2 b
!!$     !
!!$     do ii= 1, nrows
!!$        wa3(ii) = wa2(ii) * bb(ii)
!!$     end do
!!$
!!$  else  if(idprecon==SOL_DIAGONAL) then
!!$     !
!!$     ! wa2 = D^-1
!!$     !
!!$     do ii= 1, npopo
!!$        jj = ia(ii)
!!$        ll = -1
!!$        do while (jj< ia(ii+1) .and. ll ==-1)
!!$           if(ja(jj)==ii) then
!!$              ll = jj
!!$           end if
!!$           jj = jj+1
!!$        end do
!!$
!!$        if(ll/=-1) then
!!$           jj = (ii-1) * nbvar
!!$
!!$           do kk= 1, nbvar
!!$              if(wa2(jj+kk)/=0.0_rp) then
!!$                 wa2(jj+kk) = 1.0_rp /  wa2(jj+kk)
!!$              else
!!$                 mess1=intost(ii)
!!$                 mess2=intost(kk)
!!$                 messa='GMRES: ZERO DIAGONAL FOUND AT (NODE,DOF)=('&
!!$                      //trim(mess1)//','//trim(mess2)//')'                    
!!$                 if(kfl_solve_sol/=0)&
!!$                      write(lun_outso,*) trim(messa)
!!$                 call runend(messa)
!!$              end if
!!$           end do
!!$        else
!!$           mess1=intost(ii)
!!$           messa='GMRES: NO DIAGONAL FOUND FOR NODE= '//trim(mess1)              
!!$           if(kfl_solve_sol/=0)&
!!$                write(lun_outso,*) trim(messa)
!!$           call runend(messa)                                 
!!$        end if
!!$     end do
!!$     !
!!$     ! b' = D^-1 b
!!$     !
!!$     do ii= 1, nrows
!!$        wa3(ii) = wa2(ii) * bb(ii)
!!$     end do
!!$
!!$  else if(idprecon>=SOL_MATRIX) then
!!$     !
!!$     !  b' = b
!!$     !
!!$     call bcsrax( 1_ip, npopo, nbvar, pn, ja, ia, bb, wa3 )
!!$  else
!!$     !
!!$     ! b' = b
!!$     !
!!$     do ii= 1, nrows
!!$        wa3(ii) = bb(ii)
!!$     end do
!!$
!!$  end if
!!$  !
!!$  ! ||b'|| and stop criterion
!!$  !
!!$  call norm2x(nbvar,wa3,raux)
!!$
!!$  invnb   = 1.0_rp / raux
!!$  stopcri = eps * raux
!!$  !
!!$  ! Initial x' = D^1/2 x
!!$  !
!!$  if(idprecon==SOL_SQUARE) then
!!$     do ii= 1, nrows
!!$        xx(ii) = wa1(ii) * xx(ii)
!!$     end do
!!$  end if
!!$
!!$  !----------------------------------------------------------------------
!!$  ! 
!!$  ! MAIN LOOP
!!$  !
!!$  !----------------------------------------------------------------------
!!$
!!$  do while(.not.convergence)
!!$     !
!!$     ! Initial residual: kryl(*,1) = b' - [A]' x'
!!$     !
!!$#ifdef EVENT_SOLVER
!!$     call mpitrace_eventandcounters(300,1)
!!$#endif
!!$     !
!!$     ! X_0: Compute pre-initial guess
!!$     !
!!$     if(ngrou/=0) then
!!$        call bcsrax( 1_ip, npopo, nbvar, an, ja, ia, xx, wa1 ) ! A.x_{-1}
!!$        do ii= 1, nrows
!!$           wa1(ii) = bb(ii) - wa1(ii)                           ! r_{-1} = b-A.x_{-1}
!!$        end do
!!$        call wtvect(npopo,ngrou,nbvar,mu,wa1) ! W^T.r_{-1}
!!$        if(kfl_paral/=0) then
!!$           call LUsolv(&                                      ! A'.mu = W^T.r_{-1}
!!$                solve_sol(1)%ngrou*nbvar,solve_sol(1)%nskyl,&
!!$                solve_sol(1)%iskyl,1_ip,askyl,&
!!$                mu,solve_sol(1)%ngrou*nbvar,solve_sol(1)%idiag,info)
!!$           if(info/=0) call runend('MATGRU: COULD NOT SOLVE INITIAL SYSTEM')
!!$        end if
!!$        call wvect(npopo,nbvar,mu,wa1)                               ! W.mu
!!$        do ii=1,nrows
!!$           xx(ii) = xx(ii) + wa1(ii)                           ! x0 = x_{-1} + W.mu
!!$        end do
!!$        kaka=kaka+1.0_rp
!!$        !wopos(1)='SOLUT'
!!$        !wopos(2)='SCALA'
!!$        !call postpr(xx,wopos,1_ip,kaka)
!!$
!!$     end if
!!$
!!$     if(idprecon==SOL_NO_PRECOND) then
!!$        call bcsrax( 1_ip, npopo, nbvar, an, ja, ia, xx, kryl(1,1))
!!$
!!$     else if(idprecon==SOL_SQUARE) then 
!!$        do kk= 1, nrows
!!$           wa1(kk) = wa2(kk) * xx(kk)
!!$        end do
!!$
!!$        call bcsrax( 1_ip, npopo, nbvar, an, ja, ia, wa1,kryl(1,1))
!!$
!!$        do kk= 1, nrows
!!$           kryl(kk,1) = wa2(kk) * kryl(kk,1)
!!$        end do
!!$     else if(idprecon==SOL_DIAGONAL) then 
!!$        call bcsrax( 1_ip, npopo, nbvar, an, ja, ia, xx,kryl(1,1))
!!$
!!$        do kk= 1, nrows
!!$           kryl(kk,1) = wa2(kk) * kryl(kk,1)
!!$        end do
!!$     else if(idprecon>=SOL_MATRIX) then 
!!$        call bcsrax( 1_ip, npopo, nbvar, an, ja, ia, xx,kryl(1,1))
!!$        call bcsrax( 1_ip, npopo, nbvar, pn, ja, ia, kryl(1,1) ,wa4 )
!!$
!!$        do kk= 1, nrows
!!$           kryl(kk,1) = wa4(kk)
!!$        end do
!!$
!!$     end if
!!$
!!$     do kk= 1, nrows
!!$        kryl(kk,1) = wa3(kk) - kryl(kk,1)
!!$     end do
!!$     !
!!$     ! raux = ||kryl(*,1)||
!!$     !
!!$     call norm2x(nbvar,kryl,raux)
!!$
!!$     resi2=resi1
!!$     resi1=raux*invnb
!!$     if(kfl_cvgso_sol==1) write(lun_cvgso,100) solve_sol(1) % iters,raux*invnb
!!$
!!$     if(raux<=stopcri) then
!!$        !
!!$        ! The initial guess is the solution
!!$        !
!!$        convergence = .true.
!!$     else
!!$        !
!!$        ! Initialize 1-st term of the rhs of hessenberg system
!!$        !
!!$        rs(1) = raux
!!$        !
!!$        ! Ortonormalize kryl(*,1)
!!$        !
!!$        raux = 1.0_rp / raux
!!$        do kk= 1, nrows
!!$           kryl(kk,1) = kryl(kk,1) * raux
!!$        end do
!!$     endif
!!$
!!$     jj   = 0
!!$     idx  = -1
!!$     fin2 = convergence
!!$
!!$#ifdef EVENT_SOLVER
!!$     call mpitrace_eventandcounters(300,0)
!!$#endif
!!$     !
!!$     ! Inner loop. Restarted each kryldim iterations
!!$     !
!!$     do while(.not.fin2)
!!$
!!$        solve_sol(1) % iters = solve_sol(1) % iters + 1
!!$        jj    = jj + 1
!!$        jj1   = jj + 1
!!$        idx   = idx + jj
!!$        !
!!$        !  kryl(*,jj1) = [A]' kryl(*,jj)
!!$        !
!!$        !
!!$        ! Solve A'.mu = W^T.A.p
!!$        !
!!$        if(ngrou/=0) then
!!$
!!$           call bcsrax( 1_ip, npopo, nbvar, an, ja, ia, kryl(1,jj), wa1 ) ! A.p
!!$           call wtvect(npopo,ngrou,nbvar,mu,wa1)
!!$           if(kfl_paral/=0) then
!!$              call LUsolv(&                                       ! A'.mu = W^T.A.p
!!$                   solve_sol(1)%ngrou*nbvar,solve_sol(1)%nskyl,solve_sol(1)%iskyl,&
!!$                   1_ip,askyl,mu,solve_sol(1)%ngrou*nbvar,solve_sol(1)%idiag,info)
!!$              if(info/=0) call runend('MATGRO: COULD NOT SOLVE SYSTEM')
!!$           end if
!!$           call wvect(npopo,nbvar,mu,wa1)
!!$           !
!!$           !    p=p-Z.(Z^T.A.Z)^{-1}.Z^T.A.p
!!$           !
!!$           do ii= 1, nrows
!!$              wa1(ii) = kryl(ii,jj)-wa1(ii)
!!$           end do
!!$        else
!!$
!!$           do ii=1,nrows
!!$              wa1(ii) = kryl(ii,jj)                           ! x0 = x_{-1} + W.mu
!!$           end do
!!$
!!$        endif
!!$
!!$#ifdef EVENT_SOLVER
!!$        call mpitrace_eventandcounters(300,2)
!!$#endif
!!$        if(idprecon==SOL_NO_PRECOND) then
!!$           call bcsrax( 1_ip, npopo, nbvar, an, ja, ia,&
!!$                wa1, kryl(1,jj1) )
!!$
!!$        else if(idprecon==SOL_SQUARE) then
!!$           do kk= 1, nrows
!!$              wa1(kk) = wa2(kk) * kryl(kk,jj)
!!$           end do
!!$
!!$           call bcsrax( 1_ip, npopo, nbvar, an, ja, ia, &
!!$                wa1,  kryl(1,jj1) )
!!$           do kk= 1, nrows
!!$              kryl(kk,jj1) = wa2(kk) * kryl(kk,jj1)
!!$           end do
!!$        else if(idprecon==SOL_DIAGONAL) then
!!$           call bcsrax( 1_ip, npopo, nbvar, an, ja, ia,& 
!!$                kryl(1,jj),  kryl(1,jj1) )
!!$           do kk= 1, nrows
!!$              kryl(kk,jj1) = wa2(kk) * kryl(kk,jj1)
!!$           end do
!!$        else if(idprecon>=SOL_MATRIX) then
!!$           call bcsrax( 1_ip, npopo, nbvar, an, ja, ia,kryl(1,jj), kryl(1,jj1) )
!!$           call bcsrax( 1_ip, npopo, nbvar, pn, ja, ia,kryl(1,jj1),wa4 )
!!$           do kk= 1, nrows
!!$              kryl(kk,jj1) = wa4(kk)
!!$           end do
!!$        end if
!!$
!!$#ifdef EVENT_SOLVER
!!$        call mpitrace_eventandcounters(300,0)
!!$#endif
!!$        !
!!$        ! Modified Gram-Schmidt
!!$        ! For i= 1, j
!!$        !     H(i,j) = <v_i, v_j1>
!!$        !       v_j1 = v_j1 - H(i,j) * v_i
!!$        !
!!$#ifdef EVENT_SOLVER
!!$        call mpitrace_eventandcounters(300,3)
!!$#endif
!!$        do ii=1,jj
!!$
!!$           call prodxy(nbvar,npopo,kryl(1,ii),kryl(1,jj1),raux)
!!$
!!$           hh(idx+ii) = raux
!!$
!!$           do kk= 1, nrows
!!$              kryl(kk,jj1) = kryl(kk,jj1) - raux * kryl(kk,ii)
!!$           end do
!!$
!!$        end do
!!$#ifdef EVENT_SOLVER
!!$        call mpitrace_eventandcounters(300,0)
!!$#endif
!!$#ifdef EVENT_SOLVER
!!$        call mpitrace_eventandcounters(300,4)
!!$#endif
!!$        !
!!$        ! H(jj1,jj) = ||kryl(*,jj1)||
!!$        !
!!$        call norm2x(nbvar,kryl(1,jj1),raux)
!!$
!!$        hh(idx+jj1) = raux
!!$
!!$        if(raux==0.0_rp) then
!!$           fin2        = .true.
!!$           convergence = .true.
!!$           idx = idx - jj
!!$           jj = jj - 1
!!$           solve_sol(1) % iters = solve_sol(1) % iters - 1
!!$           if(kfl_solve_sol/=0)&
!!$                write(lun_outso,*) '||kryl(*,jj1)|| = 0  ', jj1, solve_sol(1) % iters
!!$           goto 10
!!$        else
!!$           !
!!$           ! Ortonormalize kryl(*,jj1)
!!$           !
!!$           raux = 1.0_rp / raux
!!$           do kk= 1, nrows
!!$              kryl(kk,jj1) = kryl(kk,jj1) * raux
!!$           end do
!!$        end if
!!$#ifdef EVENT_SOLVER
!!$        call mpitrace_eventandcounters(300,0)
!!$#endif
!!$        !
!!$        ! Update factorization of H. Perform previous 
!!$        ! transformations on jj-th column of H
!!$        !
!!$#ifdef EVENT_SOLVER
!!$        call mpitrace_eventandcounters(300,5)
!!$#endif
!!$        do ii= 1, jj-1
!!$           kk         =   ii + 1
!!$           raux       =   hh(idx+ii)
!!$           hh(idx+ii) =   cc(ii) * raux + ss(ii) * hh(idx+kk)
!!$           hh(idx+kk) = - ss(ii) * raux + cc(ii) * hh(idx+kk)
!!$        end do
!!$        gamma = hh(idx+jj)*hh(idx+jj) + hh(idx+jj1)*hh(idx+jj1)
!!$        gamma = sqrt(gamma)
!!$        !
!!$        ! if gamma is zero then take any small
!!$        ! value will affect only residual estimate
!!$        !
!!$        if(gamma==0.0_rp) then
!!$           gamma = epsmac
!!$           if(kfl_solve_sol/=0) write(lun_outso,*) 'gamma==0.0  ', solve_sol(1) % iters
!!$        end if
!!$#ifdef EVENT_SOLVER
!!$        call mpitrace_eventandcounters(300,0)
!!$#endif
!!$#ifdef EVENT_SOLVER
!!$        call mpitrace_eventandcounters(300,6)
!!$#endif
!!$        !
!!$        ! Get next plane rotation
!!$        !
!!$        gamma      =   1.0_rp / gamma
!!$        cc(jj)     =   hh(idx+jj)  * gamma
!!$        ss(jj)     =   hh(idx+jj1) * gamma
!!$        hh(idx+jj) =   cc(jj) * hh(idx+jj) + ss(jj) * hh(idx+jj1)
!!$        !
!!$        ! Update the rhs of the LS problem
!!$        !
!!$        rs(jj1)    = - ss(jj) * rs(jj)
!!$        rs(jj)     =   cc(jj) * rs(jj)
!!$        !
!!$        ! Convergence Test
!!$        !
!!$        raux       =   abs( rs(jj1) )
!!$
!!$        if(kfl_cvgso_sol==1) write(lun_cvgso,100) solve_sol(1) % iters,raux*invnb
!!$
!!$        if(raux<=stopcri) then
!!$           convergence = .true.
!!$           fin2        = .true.
!!$        else
!!$           if(solve_sol(1) % iters>=maxiter) then
!!$              convergence = .true.
!!$              fin2        = .true.
!!$           else
!!$              if(jj>=kryldim) then
!!$                 fin2 = .true.
!!$              end if
!!$           end if
!!$        end if
!!$
!!$#ifdef EVENT_SOLVER
!!$        call mpitrace_eventandcounters(300,0)
!!$#endif
!!$
!!$10      continue
!!$     end do
!!$
!!$
!!$     !----------------------------------------------------------------------
!!$     ! 
!!$     ! END INNER LOOP 
!!$     !
!!$     !----------------------------------------------------------------------
!!$
!!$#ifdef EVENT_SOLVER
!!$     call mpitrace_eventandcounters(300,7)
!!$#endif
!!$
!!$     !
!!$     ! Compute y => Solve upper triangular system
!!$     !
!!$     do ii= jj, 2, -1
!!$        rs(ii) = rs(ii) / hh(idx+ii)
!!$        raux = rs(ii)
!!$
!!$        do kk= 1, ii-1
!!$           rs(kk) = rs(kk) - hh(idx+kk) * raux
!!$        end do
!!$
!!$        idx = idx - ii
!!$     end do
!!$
!!$     rs(1) = rs(1) / hh(1)
!!$     !
!!$     ! Linear combination of kryl(*,jj)'s to get the solution.
!!$     !
!!$     do kk= 1, nrows
!!$        wa1(kk)=0_rp
!!$     enddo
!!$
!!$     do ii= 1, jj
!!$        raux = rs(ii)
!!$        do kk= 1, nrows
!!$           wa1(kk) = wa1(kk) + raux * kryl(kk,ii) 
!!$        end do
!!$     end do
!!$     !
!!$     ! Solve A'.mu = W^T.A.p
!!$     !
!!$     if(ngrou/=0) then
!!$
!!$        call bcsrax( 1_ip, npopo, nbvar, an, ja, ia, wa1, wa2 ) ! A.p
!!$        call wtvect(npopo,ngrou,nbvar,mu,wa2)
!!$        if(kfl_paral/=0) then
!!$           call LUsolv(&                                       ! A'.mu = W^T.A.p
!!$                solve_sol(1)%ngrou*nbvar,solve_sol(1)%nskyl,solve_sol(1)%iskyl,&
!!$                1_ip,askyl,mu,solve_sol(1)%ngrou*nbvar,solve_sol(1)%idiag,info)
!!$           if(info/=0) call runend('MATGRO: COULD NOT SOLVE SYSTEM')
!!$        end if
!!$        call wvect(npopo,nbvar,mu,wa2)
!!$        !
!!$        !    p=p-Z.(Z^T.A.Z)^{-1}.Z^T.A.p
!!$        !
!!$        do ii= 1, nrows
!!$           xx(ii) = xx(ii) + wa1(ii) - wa2(ii)
!!$        end do
!!$
!!$     else
!!$
!!$        do ii= 1, nrows
!!$           xx(ii) = xx(ii) + wa1(ii)
!!$        end do
!!$
!!$     end if
!!$#ifdef EVENT_SOLVER
!!$     call mpitrace_eventandcounters(300,0)
!!$#endif
!!$
!!$  end do
!!$
!!$  !-----------------------------------------------------------------
!!$  !
!!$  !  END MAIN LOOP
!!$  !
!!$  !-----------------------------------------------------------------
!!$  !
!!$  !  Final x = [D]^-1/2 x'
!!$  !
!!$  if(idprecon==SOL_SQUARE) then
!!$     do kk= 1, nrows
!!$        xx(kk) = wa2(kk) * xx(kk)
!!$     end do
!!$  end if
!!$
!!$20 continue
!!$
!!$  call bcsrax( 1_ip, npopo, nbvar, an, ja, ia, xx, wa1 )
!!$  do kk= 1, nrows
!!$     wa1(kk) = bb(kk) - wa1(kk)
!!$  end do
!!$
!!$  call norm2x(nbvar,wa1,raux)
!!$  call norm2x(nbvar,bb,gamma)
!!$  resfi = raux/gamma
!!$
!!$  if(idprecon>=SOL_MATRIX) then
!!$     call memchk(2_ip,istat,memit,'WA4','gmrpls',wa4)
!!$     deallocate(wa4,stat=istat)
!!$     if(istat/=0) call memerr(2_ip,'WA4','gmrpls',0_ip)
!!$  end if
!!$  call memchk(2_ip,istat,memit,'HH','gmrpls',hh)
!!$  deallocate(hh,stat=istat)
!!$  if(istat/=0) call memerr(2_ip,'HH','gmrpls',0_ip)
!!$  call memchk(2_ip,istat,memit,'RS','gmrpls',rs)
!!$  deallocate(rs,stat=istat)
!!$  if(istat/=0) call memerr(2_ip,'RS','gmrpls',0_ip)
!!$  call memchk(2_ip,istat,memit,'SS','gmrpls',ss)
!!$  deallocate(ss,stat=istat)
!!$  if(istat/=0) call memerr(2_ip,'SS','gmrpls',0_ip)
!!$  call memchk(2_ip,istat,memit,'CC','gmrpls',cc)
!!$  deallocate(cc,stat=istat)
!!$  if(istat/=0) call memerr(2_ip,'CC','gmrpls',0_ip)
!!$  call memchk(2_ip,istat,memit,'KRYL','gmrpls',kryl)
!!$  deallocate(kryl,stat=istat)
!!$  if(istat/=0) call memerr(2_ip,'KRYL','gmrpls',0_ip)
!!$  call memchk(2_ip,istat,memit,'WA3','gmrpls',wa3)
!!$  deallocate(wa3,stat=istat)
!!$  if(istat/=0) call memerr(2_ip,'WA3','gmrpls',0_ip)
!!$  call memchk(2_ip,istat,memit,'WA2','gmrpls',wa2)
!!$  deallocate(wa2,stat=istat)
!!$  if(istat/=0) call memerr(2_ip,'WA2','gmrpls',0_ip)
!!$  call memchk(2_ip,istat,memit,'WA1','gmrpls',wa1)
!!$  deallocate(wa1,stat=istat)
!!$  if(istat/=0) call memerr(2_ip,'WA1','gmrpls',0_ip)
!!$
!!$100 format(i7,1x,e12.6)
!!$110 format(i5,18(2x,e12.6))
!!$
!!$#ifdef EVENT
!!$  call mpitrace_user_function(0)
!!$#endif

end subroutine defgmr


