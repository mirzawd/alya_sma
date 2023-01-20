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
!> @file    pipelined_CG_rr.f90
!> @author  Guillaume Houzeaux
!> @brief   Conjugate gradient solver
!> @details Deflated pipelined CG  
!>          \verbatim
!>
!>          Initial solution x_-1
!>          if( deflation ) then
!>            Factorize A'
!>            r_-1 = b - A x_-1
!>            x0   = x_-1 + ( W A'^-1 W^T ) w_0 ) r_-1
!>            u0'  = ( W A'^-1 W^T ) w_0
!>          else
!>            x0   = x_-1
!>            u0'  = 0
!>          end if
!>          r_0   = b - A x_0
!>          M u_0 = r_0
!>          w_0   = A u_0
!>          u_0   = u_0 - u0'
!>
!>          for i = 0,...
!>
!>            gama_i = (r_i,u_i)
!>            delta   = (w_i,u_i)
!>            w'_i    = W^T w_i
!>            M m_i   = w_i
!>            n_i     = A m_i
!>            if( i > 0 ) then
!>              beta_i  = gama_i / gama_i-1
!>              alpha_i = gama_i / ( delta - beta_i * gama_i / alpha_i-1 )
!>            else
!>              beta_i  = 0
!>              alpha_i = gama_i / delta
!>            end if
!>            if( deflation ) then
!>               A' d'_i = w'_i
!>               d_i     = W d'_i
!>            else
!>               d_i     = 0
!>            end if
!>            z_i = n_i + beta_i  * z_i-1 - A M^-1 A d_i
!>            q_i = m_i + beta_i  * q_i-1 - M^-1 A d_i
!>            s_i = w_i + beta_i  * s_i-1 - A d_i
!>            p_i = u_i + beta_i  * p_i-1 - d_i
!>            x_i = x_i + alpha_i * p_i
!>            r_i = r_i - alpha_i * s_i
!>            u_i = u_i - alpha_i * q_i
!>            w_i = w_i - alpha_i * z_i
!>
!>          end for
!>
!>          \endverbatim
!> @} 
!-----------------------------------------------------------------------
subroutine pipelined_CG_rr( &
     nbnodes, nbvar, idprecon, maxiter, &
     eps, an, pn, kfl_cvgso, lun_cvgso, kfl_solve, &
     lun_outso, ja, ia, bb, xx )
  use def_kintyp,         only :  ip,rp,lg
  use def_master,         only :  IMASTER,INOTMASTER
  use def_domain,         only :  npoin_own
  use def_solver,         only :  memit,resi1,resi2,solve_sol
  use mod_solver,         only :  solver_SpMV
  use mod_solver,         only :  solver_parallel_vector_L2norm
  use mod_memory,         only :  memory_alloca
  use mod_memory,         only :  memory_deallo 
  use mod_csrdir,         only :  CSR_LU_Factorization
  use mod_csrdir,         only :  CSR_LUsol
  use mod_communications_global, only :  PAR_WAITALL_REDUCE
  use mod_communications, only :  PAR_SUM
  use mod_communications, only :  PAR_GATHER, PAR_MAX, PAR_SUM
  use mod_communications, only :  PAR_ALLGATHER
  use mod_direct_solver,  only :  direct_solver_partialcleaning
  use mod_direct_solver,  only :  direct_solver_factorization
  use mod_direct_solver,  only :  direct_solver_solution
  use mod_direct_solver,  only :  direct_allocate_temporary_matrix
  use mod_direct_solver,  only :  direct_solver_matrix_size
  use mod_communications, only :  PAR_BROADCAST
  use mod_deflated_cg,    only :  matgro
  use mod_deflated_cg,    only :  matgr2
  use mod_deflated_cg,    only :  wtvect
  use mod_deflated_cg,    only :  wvect
  use mod_deflated_cg,    only :  wtvect_without_all_reduce
  implicit none
  integer(ip), intent(in)      :: nbnodes
  integer(ip), intent(in)      :: nbvar
  integer(ip), intent(in)      :: idprecon
  integer(ip), intent(in)      :: maxiter
  integer(ip), intent(in)      :: kfl_cvgso
  integer(ip), intent(in)      :: lun_cvgso
  integer(ip), intent(in)      :: kfl_solve
  integer(ip), intent(in)      :: lun_outso
  real(rp),    intent(in)      :: eps
  real(rp),    intent(in)      :: an(nbvar,nbvar,*), pn(*)
  integer(ip), intent(in)      :: ja(*), ia(*)
  real(rp),    intent(in)      :: bb(*)
  real(rp),    intent(inout)   :: xx(nbvar*nbnodes)
  integer(ip)                  :: ii,nrows,ierr,ngrou,ntotn
  integer(ip)                  :: igama,idelta,acoarse_size
  integer(ip)                  :: ncols,nsize
  real(rp)                     :: alpha, beta, stopcri, resid
  real(rp)                     :: invnb, newrho, dummr, gama, delta
  real(rp)                     :: gama_old, alpha_old
  real(rp),    pointer         :: rr(:), pp(:), p0(:), zz(:)
  real(rp),    pointer         :: ww(:), ss(:), uu(:), qq(:), dd(:)
  real(rp),    pointer         :: aa(:), maa(:), aamaa(:), w0(:), u0(:)
  real(rp),    pointer         :: mm(:), nn(:), invdiag(:), acoarse(:)
  real(rp),    pointer         :: murhs(:),mu_gama_delta(:),mu(:)
  real(rp),    pointer         :: xo(:),uo(:),wo(:)
  logical(lg)                  :: non_blocking
  !
  ! Rdound off error corrections
  !
  logical(lg)                  :: replace
  logical(lg)                  :: roe_correction
  real(rp)                     :: zeta, tau, theta, mu_ben, Ainf, tmp
  real(rp)                     :: rho, chi, pi, sigma, xi, omega, phi, psi, nu
  real(rp)                     :: e_f, e_h, e_g, e_j, f, g, h, j, sqrt_real_n, eps_rr
  real(rp)                     :: n, rho_old, sigma_old, pi_old, phi_old, psi_old, beta_old, f_old
  integer(ip)                  :: jj, hh, ll,mu_gama_delta_size,kac_kere
  integer(ip)                  :: ichi,ipi,isigma,ixi,iomega,iphi,ipsi,inu,irho
  real(rp),     pointer        :: sums(:)
  !
  !
  ! r00  = b - A x00
  ! A'd0 = W^t r00
  ! x0   = x00 + W d0
  ! r0   = b - A x0
  ! M u0 = r0
  ! 
  ! w0   = A u0
  ! A'd  = W^t w0
  ! p0   = - W d + u0
  !
  ! s     = Ap
  ! gama = (r,u)
  ! caca  = (p,s)
  ! alpha = gama / caca
  !
  ! x   = x + alpha*p 
  ! r   = r - alpha*s
  ! M u = r
  !
  ! gama = (r,u) 
  ! beta  = gama / gama_old
  ! A'd   = W^t A u
  ! 
  ! p     = u + beta*p - W d
  !
#ifdef EVENT
  call mpitrace_user_function(1)
#endif 

  non_blocking   = .true.
  if( solve_sol(1) % kfl_roe_correction /= 0 ) then
     roe_correction = .true.
  else
     roe_correction = .false.
  end if
  ngrou          = max(0_ip,solve_sol(1) % ngrou )
  !
  ! Skyline, CSR or dense format when using deflation
  !   
  if( ngrou > 0 ) then
     acoarse_size = direct_solver_matrix_size(solve_sol(1) % direct_solver_Deflation)
  end if
  igama = nbvar * ngrou +  1
  idelta = nbvar * ngrou +  2
  ichi   = nbvar * ngrou +  3
  ipi    = nbvar * ngrou +  4
  isigma = nbvar * ngrou +  5
  ixi    = nbvar * ngrou +  6
  iomega = nbvar * ngrou +  7
  iphi   = nbvar * ngrou +  8
  ipsi   = nbvar * ngrou +  9
  inu    = nbvar * ngrou + 10
  irho   = nbvar * ngrou + 11
  if( roe_correction ) then
     mu_gama_delta_size = 11
  else
     mu_gama_delta_size = 2
  end if
  if( IMASTER ) then
     nrows   = 0 
     ncols   = 0
  else
     nrows   = nbnodes * nbvar
     ncols   = solve_sol(1) % ncols * nbvar
  end if
  ntotn = ngrou * nbvar
  nsize = max(1_ip,ncols,nrows)  
  !
  ! Initialize variables
  !
  nullify(rr)
  nullify(pp)
  nullify(p0)
  nullify(zz)
  nullify(dd)
  nullify(ww)
  nullify(ss)
  nullify(uu)
  nullify(qq)
  nullify(mm)
  nullify(nn)
  nullify(w0)
  nullify(u0)
  nullify(invdiag)
  nullify(mu_gama_delta)

  nullify(aa)
  nullify(maa)
  nullify(aamaa)
  nullify(murhs)
  nullify(mu)
  nullify(acoarse)

  nullify(xo)
  nullify(uo)
  nullify(wo)
  nullify(sums)
  !
  ! Sequential and slaves: Working arrays
  !
  call memory_alloca(memit,'RR'            ,'pipelined_CG_rr',rr            ,nsize)
  call memory_alloca(memit,'PP'            ,'pipelined_CG_rr',pp            ,nsize)
  call memory_alloca(memit,'P0'            ,'pipelined_CG_rr',p0            ,nsize)
  call memory_alloca(memit,'ZZ'            ,'pipelined_CG_rr',zz            ,nsize)
  call memory_alloca(memit,'DD'            ,'pipelined_CG_rr',dd            ,nsize)
  call memory_alloca(memit,'WW'            ,'pipelined_CG_rr',ww            ,nsize)
  call memory_alloca(memit,'SS'            ,'pipelined_CG_rr',ss            ,nsize)
  call memory_alloca(memit,'UU'            ,'pipelined_CG_rr',uu            ,nsize)
  call memory_alloca(memit,'QQ'            ,'pipelined_CG_rr',qq            ,nsize)
  call memory_alloca(memit,'MM'            ,'pipelined_CG_rr',mm            ,nsize)
  call memory_alloca(memit,'NN'            ,'pipelined_CG_rr',nn            ,nsize)
  call memory_alloca(memit,'W0'            ,'pipelined_CG_rr',w0            ,nsize)
  call memory_alloca(memit,'U0'            ,'pipelined_CG_rr',u0            ,nsize)  
  call memory_alloca(memit,'INVDIAG'       ,'pipelined_CG_rr',invdiag       ,nsize)
  call memory_alloca(memit,'MU_GAMA_DELTA','pipelined_CG_rr',mu_gama_delta,ntotn+mu_gama_delta_size)

  call memory_alloca(memit,'AA'            ,'pipelined_CG_rr',aa            ,nsize)
  call memory_alloca(memit,'MAA'           ,'pipelined_CG_rr',maa           ,nsize)
  call memory_alloca(memit,'AAMAA'         ,'pipelined_CG_rr',aamaa         ,nsize)
  call memory_alloca(memit,'MU'            ,'pipelined_CG_rr',mu            ,ntotn)

  if( roe_correction ) then    
     call memory_alloca(memit,'XO'            ,'pipelined_CG_rr',xo         ,nsize)   ! x^{i-1}
     call memory_alloca(memit,'UO'            ,'pipelined_CG_rr',uo         ,nsize)   ! u^{i-1} 
     call memory_alloca(memit,'WO'            ,'pipelined_CG_rr',wo         ,nsize)   ! w^{i-1}
     call memory_alloca(memit,'SUMS'          ,'pipelined_CG_rr',sums       ,nbvar)
  end if

  if( ngrou > 0 ) then

     !----------------------------------------------------------------------
     !
     ! Allocate memory for groups
     !
     !----------------------------------------------------------------------

     call memory_alloca(memit , 'MURHS'   , 'pipelined_CG_rr' , murhs   , ntotn)
     call memory_alloca(memit , 'ACOARSE' , 'pipelined_CG_rr' , acoarse , acoarse_size)

     !----------------------------------------------------------------------
     !
     ! Compute A'= W^T.A.W and factorize A'
     !
     !----------------------------------------------------------------------

     if( solve_sol(1) % kfl_symme == 1 ) then
        call matgro(ngrou,nbnodes,acoarse_size,nbvar,ia,ja,an,acoarse)
     else
        call matgr2(ngrou,nbnodes,acoarse_size,nbvar,ia,ja,an,acoarse)
     end if
     if( INOTMASTER ) call direct_solver_factorization(solve_sol(1) % direct_solver_Deflation,acoarse) 

     !----------------------------------------------------------------------
     !
     ! X_0: Compute pre-initial guess 
     !
     !----------------------------------------------------------------------

     if( solve_sol(1) % kfl_schum == 1 ) then
        call bcsrax_schur( 1_ip, nbnodes, nbvar, solve_sol(1) % ndofn_A3 , &
             solve_sol(1) % A1, solve_sol(1) % A2, solve_sol(1) % invA3, &
             solve_sol(1) % A4, ja, ia, xx, rr )   
     else
        call solver_SpMV(solve_sol(1),an,xx,rr)                         ! A.x_{-1}
     end if
     do ii= 1, nrows
        rr(ii) = bb(ii) - rr(ii)                                      ! r_{-1} = b-A.x_{-1}
     end do

     call wtvect(nbnodes,ngrou,nbvar,mu_gama_delta,rr)                 ! W^T.r_{-1}
     if( INOTMASTER ) then   
        murhs(1:ntotn) = mu_gama_delta(1:ntotn)
        call direct_solver_solution(solve_sol(1) % direct_solver_Deflation,murhs,mu_gama_delta) 
     end if
     call wvect(nbnodes,nbvar,mu_gama_delta,rr)                        ! W.mu
     do ii = 1,nrows
        xx(ii) = xx(ii) + rr(ii)                                      ! x0 = x_{-1} + W.mu
     end do
  end if

  !----------------------------------------------------------------------
  !
  ! Initial computations
  ! r0    = b - A x0
  ! M u_0 = r0
  ! w0    = A u0
  !
  !----------------------------------------------------------------------
  
  call solope(&
       1_ip, nbvar, idprecon, eps, an, pn, ja, ia, bb, xx , &
       ierr, stopcri, newrho, resid, invnb, rr, uu, pp, ww, &
       invdiag, p0 )
  call solver_SpMV(solve_sol(1),an,uu,ww)
  
  if( ierr /= 0 ) goto 10
  !
  ! Deflation
  ! 
  if( ngrou > 0 ) then
     !
     ! Needed to compute alpha_0 = (r0,r0) / (Ap0,p0) (p0 includes deflation)
     !
     call wtvect(nbnodes,ngrou,nbvar,mu_gama_delta,ww)
     if( INOTMASTER ) then
        murhs(1:ngrou*nbvar) = mu_gama_delta(1:ngrou*nbvar)
        call direct_solver_solution(solve_sol(1) % direct_solver_Deflation,murhs,mu_gama_delta) 
     end if
     call wvect(nbnodes,nbvar,mu_gama_delta,dd)
     do ii = 1,nrows
        u0(ii) = uu(ii) - dd(ii)
     end do
  else
     !
     ! No deflation
     !
     do ii = 1,nrows
        u0(ii) = uu(ii) 
     end do
  end if
  !
  ! w0 = A u0 => delta0 = (Ap0,p0)
  !
  call solver_SpMV(solve_sol(1),an,u0,w0)

  !----------------------------------------------------------------------
  !
  ! ROUND-OFF ERROR CONTROL
  !
  ! zeta    = norm(b,2);
  ! tau     = sqrt(epsilon);
  ! n       = length(b);
  ! theta   = sqrt(n) * L_inf(A)
  ! mu      = L_inf(A) = maximum number of nonzeros in any row of A
  ! replace = false
  !
  !----------------------------------------------------------------------

  if( roe_correction ) then
     call solver_parallel_vector_L2norm(solve_sol(1),bb,zeta)
     eps_rr      = 1.0e-32_rp
     tau         = sqrt(eps_rr)
     n           = real(npoin_own,rp)
     call PAR_SUM(n,'IN MY CODE')
     sqrt_real_n = sqrt(n)
     !
     ! Linf norm of the matrix
     !
     do ii = 1,nbvar
        sums(ii) = 0.0_rp
     end do

     Ainf    = 0.0_rp
     tmp     = 0.0_rp
     mu_ben  = 0.0_rp
     replace = .false.

     do ii = 1,nrows
        mu_ben = max(mu_ben,real(ia(ii+1)-ia(ii),rp))     
        do jj = ia(ii),ia(ii+1)-1
           do hh = 1,nbvar
              do ll = 1,nbvar
                 sums(ll) = sums(ll) + abs(an(hh,ll,jj))
              end do
           end do
        end do
        tmp = maxval(sums)
        if( tmp > Ainf ) then
           Ainf = tmp
        end if
        do ll = 1,nbvar
           sums(ll) = 0.0_rp
        end do
     end do
     mu_ben = mu_ben * real(nbvar,rp)
     call PAR_MAX(mu_ben, 'IN MY CODE')
     call PAR_MAX(Ainf,   'IN MY CODE')

     theta = sqrt_real_n * Ainf
  end if

100 continue
  
  if( roe_correction ) then
     do ii = 1,nrows
        xo(ii) = xx(ii)
        uo(ii) = uu(ii)
        wo(ii) = ww(ii)
     end do
     f          = 0.0_rp
     g          = 0.0_rp
     h          = 0.0_rp
     j          = 0.0_rp
     e_g        = 0.0_rp
     e_j        = 0.0_rp
     e_f        = 0.0_rp
     e_h        = 0.0_rp
     chi        = 0.0_rp
     pi         = 0.0_rp
     sigma      = 0.0_rp
     xi         = 0.0_rp
     omega      = 0.0_rp
     phi        = 0.0_rp
     psi        = 0.0_rp
     nu         = 0.0_rp
     rho        = 0.0_rp
     sigma_old  = 0.0_rp
     phi_old    = 0.0_rp
     pi_old     = 0.0_rp
     psi_old    = 0.0_rp
     kac_kere   = 0
  end if
  alpha = 0.0_rp
  gama = 0.0_rp
  
  !-----------------------------------------------------------------
  !
  !  MAIN LOOP
  !
  !-----------------------------------------------------------------

  do while( solve_sol(1) % iters < maxiter .and. resid > stopcri )
     !
     ! Deflation
     ! gama = (r,u)
     ! delta = (w,u)
     ! rho   = norm(r,2)
     !
     if( ngrou > 0 ) call wtvect_without_all_reduce(nbnodes,ngrou,nbvar,mu_gama_delta,ww)

     if( roe_correction ) then
        call pipelined_CG_rr_norms_rr(nbvar,uu,rr,ww,xo,pp,ss,uo,wo,qq,zz,mm,mu_gama_delta(ntotn+1:))
     else
        call pipelined_CG_rr_norms(nbvar,uu,rr,ww,mu_gama_delta(ntotn+1:))
     end if
     if( non_blocking ) then
        call PAR_SUM(ntotn+mu_gama_delta_size,mu_gama_delta,'IN MY ZONE','NON BLOCKING')
     else
        call PAR_SUM(ntotn+mu_gama_delta_size,mu_gama_delta,'IN MY ZONE')
     end if
     !
     ! L m^i = w^i
     ! n^i   = A m^i
     !
     call precon(&
          3_ip,nbvar,nbnodes,nrows,solve_sol(1) % kfl_symme,idprecon,ia,ja,an,&
          pn,invdiag,w0,ww,mm)
     call solver_SpMV(solve_sol(1),an,mm,nn)
     ! 
     ! Wait all for the all reduce
     !
     if( non_blocking ) call PAR_WAITALL_REDUCE()

     gama_old             = gama
     alpha_old            = alpha
     gama                 = mu_gama_delta(igama)    ! gama
     delta                = mu_gama_delta(idelta)   ! delta
     solve_sol(1) % xorth = 0.0_rp                  ! Orthogonality check
     !
     ! Round off error correction variable
     !
     if( roe_correction .and. solve_sol(1) % iters > 0 ) then
        rho_old   = rho     
        sigma_old = sigma
        pi_old    = pi
        phi_old   = phi
        psi_old   = psi
        chi       = sqrt(mu_gama_delta(ichi))
        pi        = sqrt(mu_gama_delta(ipi))
        sigma     = sqrt(mu_gama_delta(isigma))    
        xi        = sqrt(mu_gama_delta(ixi))
        omega     = sqrt(mu_gama_delta(iomega))   
        phi       = sqrt(mu_gama_delta(iphi))    
        psi       = sqrt(mu_gama_delta(ipsi))      
        nu        = sqrt(mu_gama_delta(inu))
        rho       = sqrt(mu_gama_delta(irho))       
     end if
     !
     ! Beta and Alpha
     !
     if( solve_sol(1) % iters > 0 ) then
        beta_old = beta
        beta     = gama / gama_old
        alpha    = 1.0_rp / ( delta / gama - beta / alpha_old )
        !
        ! To obtain same iterate as original DCG, uncomment following lines
        !
        !if( ngrou /= 0 ) then 
        !   call bcsrax( 1_ip, nbnodes, nbvar, an, ja, ia, uu, maa)
        !   call wtvect(nbnodes,ngrou,nbvar,mu,maa)
        !   if( INOTMASTER ) then
        !      if( solve_sol(1) % kfl_defas == 0 ) then
        !         call chosol(&
        !              solve_sol(1) % ngrou*nbvar,solve_sol(1) % acoarse_size,&
        !              solve_sol(1) % iskyl,1_ip,solve_sol(1) % askyldef,mu,&
        !              solve_sol(1) % ngrou*nbvar,info) 
        !      else
        !         do igrou = 1,ngrou*nbvar
        !            murhs(igrou) = mu(igrou)
        !         end do
        !         call CSR_LUsol(                                          &
        !              ngrou                   , nbvar                   , &
        !              solve_sol(1) % invpRdef , solve_sol(1) % invpCdef , &
        !              solve_sol(1) % ILdef    , solve_sol(1) % JLdef    , &
        !              solve_sol(1) % LNdef    , solve_sol(1) % IUdef    , &
        !              solve_sol(1) % JUdef    , solve_sol(1) % UNdef    , &
        !              murhs                   , mu                      )   
        !      end if
        !   end if
        !   call wvect(nbnodes,nbvar,mu,dd)
        !end if
        !do ii = 1,nrows
        !   u0(ii) = uu(ii) - dd(ii)             ! u^{i+1}
        !end do
        !call bcsrax( 1_ip, nbnodes, nbvar, an, ja, ia, u0,w0)
        !call prodxy(nbvar,nbnodes,u0,w0,delta)
        !alpha = 1.0_rp / ( delta / gama - beta / alpha_old )
     else
        beta  = 0.0_rp
        alpha = gama / delta
     end if  
     !
     ! Solve A'.mu = W^T.A.u^k, d = W mu
     !
     if( ngrou > 0 ) then 
        if( INOTMASTER ) then          
           murhs(1:ngrou*nbvar) = mu_gama_delta(1:ngrou*nbvar)           
           call direct_solver_solution(solve_sol(1) % direct_solver_Deflation,murhs,mu_gama_delta) 
        end if
        call wvect(nbnodes,nbvar,mu_gama_delta,dd)
        !
        ! Arrays needed for deflation
        !
        call solver_SpMV(solve_sol(1),an,dd,aa)                       ! A d
        call precon(&
             3_ip,nbvar,nbnodes,nrows,solve_sol(1) % kfl_symme,&    ! M^{-1} A d
             idprecon,ia,ja,an,pn,invdiag,w0,aa,maa) 
        call solver_SpMV(solve_sol(1),an,maa, aamaa )                 ! A M^{-1} A d
     end if
     !
     ! Save old solution
     !
     if( roe_correction ) then
        do ii = 1,nrows
           xo(ii) = xx(ii)
           uo(ii) = uu(ii)
           wo(ii) = ww(ii)
        end do
     end if
     !
     ! Same iterate as original DCG
     !
     !do ii = 1,nrows
     !   zz(ii) = nn(ii) + beta  * zz(ii) - aamaa(ii) ! z^i
     !   qq(ii) = mm(ii) + beta  * qq(ii) - maa(ii)   ! q^i
     !   ss(ii) = ww(ii) + beta  * ss(ii) - aa(ii)    ! s^i
     !   pp(ii) = uu(ii) + beta  * pp(ii) - dd(ii)    ! p^i
     !end do
     !call prodxy(nbvar,nbnodes,pp,ss,raux1)
     !alpha = gama / raux1
     !if(imaster) print*,abs(gama / raux1-alpha)     
     !
     ! Updates
     !
#ifdef BLAS
     if( INOTMASTER ) then
        call DAXPY(nrows,beta,zz,1_ip,nn,1_ip)             ! n =  n + beta*z
        call DCOPY(nrows,nn,1_ip,zz,1_ip)                  ! z <= n
        call DAXPY(nrows,beta,qq,1_ip,mm,1_ip)             ! m =  m + beta*q
        call DCOPY(nrows,mm,1_ip,qq,1_ip)                  ! q <= m

        call DCOPY(nrows,ww,1_ip,w0,1_ip)
        call DAXPY(nrows,beta,ss,1_ip,w0,1_ip)             ! w =  w + beta*s
        call DCOPY(nrows,w0,1_ip,ss,1_ip)                  ! s <= w

        call DCOPY(nrows,uu,1_ip,w0,1_ip)
        call DAXPY(nrows,beta,pp,1_ip,w0,1_ip)             ! u =  u + beta*p
        call DCOPY(nrows,w0,1_ip,pp,1_ip)                  ! p <= u

        call DAXPY(nrows, alpha,pp,1_ip,xx,1_ip)           ! x =  x + alpha*p
        call DAXPY(nrows,-alpha,ss,1_ip,rr,1_ip)           ! r =  r - alpha*s
        call DAXPY(nrows,-alpha,qq,1_ip,uu,1_ip)           ! u =  u - alpha*q
        call DAXPY(nrows,-alpha,zz,1_ip,ww,1_ip)           ! w =  w -alpha*z
        if( ngrou > 0 ) then
           call DAXPY(nrows,-1.0_rp,aamaa,1_ip,zz,1_ip) 
           call DAXPY(nrows,-1.0_rp,maa  ,1_ip,qq,1_ip)
           call DAXPY(nrows,-1.0_rp,aa   ,1_ip,ss,1_ip)
           call DAXPY(nrows,-1.0_rp,dd   ,1_ip,pp,1_ip)
        end if
     end if
#else
     if( ngrou > 0 ) then
        do ii = 1,nrows 
           zz(ii) = nn(ii) + beta  * zz(ii) - aamaa(ii)    ! z^i
           qq(ii) = mm(ii) + beta  * qq(ii) - maa(ii)      ! q^i
           ss(ii) = ww(ii) + beta  * ss(ii) - aa(ii)       ! s^i
           pp(ii) = uu(ii) + beta  * pp(ii) - dd(ii)       ! p^i
           xx(ii) = xx(ii) + alpha * pp(ii)                ! x^{i+1}
           rr(ii) = rr(ii) - alpha * ss(ii)                ! r^{i+1}
           uu(ii) = uu(ii) - alpha * qq(ii)                ! u^{i+1}
           ww(ii) = ww(ii) - alpha * zz(ii)                ! w^{i+1}
        end do
     else
        do ii = 1,nrows 
           zz(ii) = nn(ii) + beta  * zz(ii)                ! z^i
           qq(ii) = mm(ii) + beta  * qq(ii)                ! q^i
           ss(ii) = ww(ii) + beta  * ss(ii)                ! s^i
           pp(ii) = uu(ii) + beta  * pp(ii)                ! p^i
           xx(ii) = xx(ii) + alpha * pp(ii)                ! x^{i+1}
           rr(ii) = rr(ii) - alpha * ss(ii)                ! r^{i+1}
           uu(ii) = uu(ii) - alpha * qq(ii)                ! u^{i+1}
           ww(ii) = ww(ii) - alpha * zz(ii)                ! w^{i+1}
        end do
     end if
#endif
     !
     ! ROUND-OFF ERROR CORRECTION
     !
     if( roe_correction .and. solve_sol(1) % iters > 0 ) then

        e_f = theta * chi + 2.0_rp * abs(alpha_old) * theta * pi  + rho   + 2.0_rp * abs(alpha_old) * sigma
        e_h = theta * xi  + 2.0_rp * abs(alpha_old) * theta * phi + omega + 2.0_rp * abs(alpha_old) * psi

        if (solve_sol(1) % iters > 1) then
           e_g = theta * xi + 2.0_rp * abs(beta_old) * theta * pi_old + omega + 2.0_rp * abs(beta_old) * sigma_old
           e_j = (mu_ben * sqrt_real_n + 2.0_rp) * theta * nu + 2.0_rp * abs(beta_old) * theta * phi_old + 2.0_rp * abs(beta_old) * psi_old
        end if

        if( solve_sol(1) % iters == 1 .or. replace ) then
           f_old   = f
           f       = eps_rr * sqrt((mu_ben * sqrt_real_n + 1.0_rp) * theta * chi + zeta) + eps_rr * sqrt(abs(alpha_old) * mu_ben * sqrt_real_n * theta * pi) + sqrt(e_f) * eps_rr
           g       = eps_rr * sqrt(mu_ben * sqrt_real_n * theta * pi)
           h       = eps_rr * sqrt(mu_ben * sqrt_real_n * theta * xi) + eps_rr * sqrt(abs(alpha_old) * mu_ben * sqrt_real_n * theta * phi) + sqrt(e_h) * eps_rr
           j       = eps_rr * sqrt(mu_ben * sqrt_real_n * theta * phi)
           replace = .false.
        else
           f_old = f
           f     = f + abs(alpha_old) * abs(beta_old) * g + abs(alpha_old) * h + sqrt(e_f) * eps_rr + abs(alpha_old) * sqrt(e_g) * eps_rr
           g     = abs(beta_old) * g + h + sqrt(e_g) * eps_rr
           h     = h + abs(alpha_old) * abs(beta_old) * j + sqrt(e_h) * eps_rr + abs(alpha_old) * sqrt(e_j) * eps_rr
           j     = abs(beta_old) * j + sqrt(e_j) * eps_rr
        end if

        if( f_old < tau*rho_old .and. f > tau*rho ) then
           !
           ! Orthogonality check with p0= <Ap,p0> / ( ||Ap|| ||p0|| )
           !
           !if( solve_sol(1) % kfl_symme == 1 ) then
           !   call bsymax( 1_ip, nbnodes, nbvar, an, ja, ia, pp, w0 )
           !else if( solve_sol(1) % kfl_schum == 1 ) then
           !   call bcsrax_schur( 1_ip, nbnodes, nbvar, solve_sol(1) % ndofn_A3 , &
           !        solve_sol(1) % A1, solve_sol(1) % A2, solve_sol(1) % invA3, solve_sol(1) % A4,&
           !        invdiag, ja, ia, pp, w0 )
           !else
           !   call bcsrax( 1_ip, nbnodes, nbvar, an, ja, ia, pp, w0 )
           !end if
           !call cosixy(nbvar,nbnodes,w0,p0,solve_sol(1) % xorth)
           !print*,'popo=',solve_sol(1) % xorth 
           !
           ! Recompute r = b - Ax, Mu = r, w = Au
           !
           kac_kere = kac_kere + 1
           if( IMASTER ) print *, "iter_no : ", solve_sol(1) % iters

           call solver_SpMV(solve_sol(1),an,xx,rr)    
           do ii = 1,nrows
              rr(ii) = bb(ii) - rr(ii)
           end do
           call precon(&
                3_ip,nbvar,nbnodes,nrows,solve_sol(1) % kfl_symme,idprecon,ia,ja,an,&
                pn,invdiag,w0,rr,uu)
           call solver_SpMV(solve_sol(1),an,uu,ww)  

           if( solve_sol(1) % kfl_roe_correction == 1 ) then
              !
              ! Restart
              !
              do ii = 1,nrows
                 p0(ii) = uu(ii)
              end do
              solve_sol(1) % iters = 0
              replace = .false.
              goto 100
              
           else if( solve_sol(1) % kfl_roe_correction == 2 ) then
              !
              ! Incomplete restart, recompute: s = Ap, Mq = s, z = Aq
              !
              replace = .true.
              call solver_SpMV(solve_sol(1),an,pp,ss)                
              call precon(&
                   3_ip,nbvar,nbnodes,nrows,solve_sol(1) % kfl_symme,idprecon,ia,ja,an,&
                   pn,invdiag,w0,ss,qq)
              call solver_SpMV(solve_sol(1),an,qq,zz)  
              
           end if

        end if

     end if
     !
     ! Orthogonality check with p0= <Ap,p0> / ( ||Ap|| ||p0|| )
     !
     !if( solve_sol(1) % kfl_schum == 1 ) then
     !   call bcsrax_schur( 1_ip, nbnodes, nbvar, solve_sol(1) % ndofn_A3 , &
     !        solve_sol(1) % A1, solve_sol(1) % A2, solve_sol(1) % invA3, solve_sol(1) % A4,&
     !        invdiag, ja, ia, pp, w0 )
     !else
     !   call solver_SpMV(solve_sol(1),an,pp,w0)  
     !end if
     !call cosixy(nbvar,nbnodes,w0,p0,solve_sol(1) % xorth)

     resid  = sqrt(gama)
     resi2  = resi1
     resi1  = resid * invnb
     solve_sol(1) % iters  = solve_sol(1) % iters + 1
     !
     ! Convergence post process:
     ! kk    = iteration number
     ! resi1 = preconditioned residual
     !
     if( kfl_cvgso == 1 ) &
          call outcso(nbvar,nrows,nbnodes,invdiag,an,ja,ia,bb,xx)
  end do

  !-----------------------------------------------------------------
  !
  !  END MAIN LOOP
  !
  !-----------------------------------------------------------------

  
  !if (IMASTER ) then
  !   print *, "entered : " ,kac_kere, " times!!"
  !endif
  
10 continue
  call solope(&
       2_ip, nbvar, idprecon, dummr, an, dummr, ja, ia, bb, xx , &
       ierr, dummr, dummr, resi1, dummr, rr, w0, pp, dummr, &
       invdiag, dummr ) 

  if( kfl_solve == 1 ) then
     if( ierr > 0 ) write(lun_outso,120) solve_sol(1) % iters
  end if

  !-----------------------------------------------------------------
  !
  ! Deallocate memory
  !
  !-----------------------------------------------------------------

  call direct_solver_partialcleaning(solve_sol(1) % direct_solver_Deflation) 

  call memory_deallo(memit,'RR'            ,'pipelined_CG_rr',rr)
  call memory_deallo(memit,'PP'            ,'pipelined_CG_rr',pp)
  call memory_deallo(memit,'P0'            ,'pipelined_CG_rr',p0)
  call memory_deallo(memit,'ZZ'            ,'pipelined_CG_rr',zz)
  call memory_deallo(memit,'DD'            ,'pipelined_CG_rr',dd)
  call memory_deallo(memit,'WW'            ,'pipelined_CG_rr',ww)
  call memory_deallo(memit,'SS'            ,'pipelined_CG_rr',ss)
  call memory_deallo(memit,'UU'            ,'pipelined_CG_rr',uu)
  call memory_deallo(memit,'QQ'            ,'pipelined_CG_rr',qq)
  call memory_deallo(memit,'MM'            ,'pipelined_CG_rr',mm)
  call memory_deallo(memit,'NN'            ,'pipelined_CG_rr',nn)
  call memory_deallo(memit,'W0'            ,'pipelined_CG_rr',w0)
  call memory_deallo(memit,'U0'            ,'pipelined_CG_rr',u0)
  call memory_deallo(memit,'XO'            ,'pipelined_CG_rr',xo)
  call memory_deallo(memit,'UO'            ,'pipelined_CG_rr',uo)
  call memory_deallo(memit,'WO'            ,'pipelined_CG_rr',wo)
  call memory_deallo(memit,'INVDIAG'       ,'pipelined_CG_rr',invdiag)
  call memory_deallo(memit,'AA'            ,'pipelined_CG_rr',aa)
  call memory_deallo(memit,'MAA'           ,'pipelined_CG_rr',maa)
  call memory_deallo(memit,'AAMAA'         ,'pipelined_CG_rr',aamaa)
  call memory_deallo(memit,'MU_GAMA_DELTA' ,'pipelined_CG_rr',mu_gama_delta)
  call memory_deallo(memit,'MURHS'         ,'pipelined_CG_rr',murhs)
  call memory_deallo(memit,'ACOARSE'       ,'pipelined_CG_rr',acoarse)
  call memory_deallo(memit,'SUMS'          ,'pipelined_CG_rr',sums)

120 format(&
       & '# Error at iteration ',i6,&
       & 'Dividing by zero: alpha = gama / ( delta - beta * gama / alpha_old )')

#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine pipelined_CG_rr

subroutine pipelined_CG_rr_norms_rr(nbvar,uu,rr,ww,xo,pp,ss,uo,wo,qq,zz,mm,dotprod)
  use def_kintyp, only : ip,rp
  use def_domain, only : npoin_own
  use def_master, only : INOTMASTER
  implicit none
  
  integer(ip),  intent(in)            :: nbvar
  real(rp),     intent(in)            :: uu(*)
  real(rp),     intent(in)            :: rr(*)
  real(rp),     intent(in)            :: ww(*)
  real(rp),     intent(in)            :: xo(*)
  real(rp),     intent(in)            :: pp(*)
  real(rp),     intent(in)            :: ss(*)
  real(rp),     intent(in)            :: uo(*)
  real(rp),     intent(in)            :: wo(*)
  real(rp),     intent(in)            :: qq(*)
  real(rp),     intent(in)            :: zz(*)
  real(rp),     intent(in)            :: mm(*)
  real(rp),     intent(out)           :: dotprod(11)
  integer(ip)                         :: nn
#ifdef BLAS
  real(rp) :: DDOT
  external DDOT
#endif
  
  dotprod = 0.0_rp
  nn      = nbvar * npoin_own
  
  if( INOTMASTER ) then
#ifdef BLAS
     dotprod(11) = DDOT( nn,rr,1_ip,rr,1_ip )
     dotprod( 1) = DDOT( nn,uu,1_ip,rr,1_ip )
     dotprod( 2) = DDOT( nn,uu,1_ip,ww,1_ip )
     dotprod( 3) = DDOT( nn,xo,1_ip,xo,1_ip )
     dotprod( 4) = DDOT( nn,pp,1_ip,pp,1_ip )
     dotprod( 5) = DDOT( nn,ss,1_ip,ss,1_ip )
     dotprod( 6) = DDOT( nn,uo,1_ip,uo,1_ip )
     dotprod( 7) = DDOT( nn,wo,1_ip,wo,1_ip )
     dotprod( 8) = DDOT( nn,qq,1_ip,qq,1_ip )
     dotprod( 9) = DDOT( nn,zz,1_ip,zz,1_ip )
     dotprod(10) = DDOT( nn,mm,1_ip,mm,1_ip )
#else
     dotprod(11) = dot_product( rr(1:nn) , rr(1:nn) )
     dotprod( 1) = dot_product( uu(1:nn) , rr(1:nn) )
     dotprod( 2) = dot_product( uu(1:nn) , ww(1:nn) )
     dotprod( 3) = dot_product( xo(1:nn) , xo(1:nn) )
     dotprod( 4) = dot_product( pp(1:nn) , pp(1:nn) )
     dotprod( 5) = dot_product( ss(1:nn) , ss(1:nn) )
     dotprod( 6) = dot_product( uo(1:nn) , uo(1:nn) )
     dotprod( 7) = dot_product( wo(1:nn) , wo(1:nn) )
     dotprod( 8) = dot_product( qq(1:nn) , qq(1:nn) )
     dotprod( 9) = dot_product( zz(1:nn) , zz(1:nn) )
     dotprod(10) = dot_product( mm(1:nn) , mm(1:nn) )
#endif
  end if

end subroutine pipelined_CG_rr_norms_rr

subroutine pipelined_CG_rr_norms(nbvar,uu,rr,ww,dotprod)
  use def_kintyp, only : ip,rp
  use def_domain, only : npoin_own
  use def_master, only : INOTMASTER
  implicit none
  
  integer(ip),  intent(in)            :: nbvar
  real(rp),     intent(in)            :: uu(*)
  real(rp),     intent(in)            :: rr(*)
  real(rp),     intent(in)            :: ww(*)
  real(rp),     intent(out)           :: dotprod(2)
  integer(ip)                         :: nn
#ifdef BLAS
  real(rp) :: DDOT
  external DDOT
#endif
  
  dotprod = 0.0_rp
  nn      = nbvar * npoin_own
  

  if( INOTMASTER ) then
#ifdef BLAS
     dotprod( 1) = DDOT( nn,uu,1_ip,rr,1_ip )
     dotprod( 2) = DDOT( nn,uu,1_ip,ww,1_ip )
#else
     dotprod( 1) = dot_product( uu(1:nn) , rr(1:nn) )
     dotprod( 2) = dot_product( uu(1:nn) , ww(1:nn) )
#endif
  end if

end subroutine pipelined_CG_rr_norms
