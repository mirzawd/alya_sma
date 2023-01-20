!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @defgroup Krylov_Solver
!> @ingroup  Algebraic_Solver
!> @{
!> @file    gmrpls.f90
!> @author  Guillaume Houzeaux
!> @brief   GMRES right preconditioned solver
!> @details GMRES right preconditioned solver:
!!          @verbatim
!!          1. Start: Choose x0 and a dimension m
!!          2. Arnoldi process:
!!             - Compute r0 = b - Ax0, b = ||r0||_2 and v1 = r0/b.
!!             - For j = 1,...,m do
!!               - Compute zj := M -1vj
!!               - Compute w := Azj
!!               - for i = 1,...,j, do :
!!                   hi,j := (w, vi)
!!                   w    := w - hi,jvi
!!                 enddo
!!             - hj+1,1 = ||w||2
!!             - vj+1 = w/hj+1,1
!!             - Define Vm := [v1,....,vm] and Hm = {hi,j}.
!!          3. Form the approximate solution: xm = x0 + M -1Vmym
!!             where ym = argmin_y||b*e1-Hmy||2 and e1 = [1, 0, . . . , 0]T .
!!          4. Restart: If satisfied stop, else set x0 <= xm and goto 2.
!!          @endverbatim
!> @}
!-----------------------------------------------------------------------
subroutine gmrpls( &
     nbnodes, nbvar,&
     idprecon, maxiter, kryldim, kfl_ortho, &
     kfl_cvgso, lun_cvgso,&
     kfl_solve, lun_outso,&
     eps, an, pn, ja, ia, bb, xx )

  use def_kintyp,    only  :  ip,rp,lg 
  use def_master,    only  :  IMASTER
  use def_solver,    only  :  memit,resi1,resi2,solve_sol
  use mod_memory,    only  :  memory_alloca
  use mod_memory,    only  :  memory_deallo
  use mod_couplings, only  :  COU_INTERPOLATE_NODAL_VALUES_go
  use mod_couplings, only  :  I_AM_IN_COUPLING
  use mod_couplings, only  :  I_AM_INVOLVED_IN_A_COUPLING_TYPE
  use mod_parall,    only  :  PAR_NODE_NUMBER_OF_PARTITIONS
  use def_coupli,    only  :  mcoup
  use def_coupli,    only  :  coupling_type
  use def_coupli,    only  :  DIRICHLET_IMPLICIT
  use def_coupli,    only  :  BETWEEN_SUBDOMAINS
  use def_master,    only  :  INOTMASTER
  use def_coupli,    only  :  coupling_type
  use mod_solver,    only  :  solver_krylov_subspace_save
  use mod_solver,    only  :  solver_parallel_vector_L2norm
  use mod_solver,    only  :  solver_parallel_scalar_product
  implicit none
  integer(ip), intent(in)    :: nbnodes,nbvar,idprecon,maxiter,kryldim
  integer(ip), intent(in)    :: kfl_ortho
  integer(ip), intent(in)    :: kfl_cvgso,lun_cvgso,kfl_solve,lun_outso
  real(rp),    intent(in)    :: eps
  real(rp),    intent(inout) :: an(nbvar, nbvar, *)
  real(rp),    intent(inout) :: pn(nbvar, nbvar, *)
  integer(ip), intent(in)    :: ja(*),ia(*)
  real(rp),    intent(in)    :: bb(*)
  real(rp),    intent(inout) :: xx(*)
  logical(lg)                :: convergence,fin2
  integer(ip)                :: ii,jj,kk,jj1,idx,nrows,ierr
  real(rp)                   :: resid,gama,stopcri,invnb,dummr
  real(rp)                   :: time1,time2
  real(rp),    parameter     :: epsmac=1.0e-16_rp

  real(rp),    pointer       :: wa1(:), wa3(:), kryl(:,:)
  real(rp),    pointer       :: invdiag(:)
  real(rp),    pointer       :: cc(:), ss(:), rs(:), hh(:)

  integer(ip), pointer       :: lneig(:)
  real(rp),    pointer       :: xx_tmp(:)
  integer(ip)                :: icoup,ncols

  nullify(wa1,wa3,kryl,invdiag,cc,ss,rs,hh,lneig,xx_tmp)
  convergence = .false.
  fin2        = .false.

  !---------------------------------------------------------------------
  !
  ! Memory
  !
  !----------------------------------------------------------------------

  if( IMASTER ) then
     nrows = 0
     ncols = 0
  else
     nrows = nbnodes * nbvar
     ncols = solve_sol(1) % ncols * nbvar
  end if

  call memory_alloca(memit,'CC'     ,'gmrpls',cc     ,kryldim)
  call memory_alloca(memit,'SS'     ,'gmrpls',ss     ,kryldim)
  call memory_alloca(memit,'RS'     ,'gmrpls',rs     ,kryldim+1_ip)
  call memory_alloca(memit,'HH'     ,'gmrpls',hh     ,kryldim+((kryldim*(kryldim+1_ip))/2_ip),'DO_NOT_INITIALIZE')
  call memory_alloca(memit,'WA1'    ,'gmrpls',wa1    ,max(1_ip,ncols,nrows))
  call memory_alloca(memit,'WA3'    ,'gmrpls',wa3    ,max(1_ip,ncols,nrows))
  call memory_alloca(memit,'KRYL'   ,'gmrpls',kryl   ,max(1_ip,nrows),kryldim+1_ip,'DO_NOT_INITIALIZE')
  call memory_alloca(memit,'INVDIAG','gmrpls',invdiag,max(1_ip,ncols,nrows))

  if( INOTMASTER .and. I_AM_INVOLVED_IN_A_COUPLING_TYPE(BETWEEN_SUBDOMAINS,DIRICHLET_IMPLICIT) ) then
     nullify(lneig)
     allocate(lneig(nbnodes))
     call PAR_NODE_NUMBER_OF_PARTITIONS(nbnodes,lneig)
  end if

  !---------------------------------------------------------------------
  !
  ! Initial computations
  !
  !----------------------------------------------------------------------

  call solope(&
       1_ip, nbvar, idprecon, eps, an, pn, ja, ia, bb, xx , &
       ierr, stopcri, dummr, resid, invnb, kryl, dummr, dummr, wa1, &
       invdiag , wa3 )

  if( ierr /= 0 ) goto 20

  !----------------------------------------------------------------------
  !
  ! MAIN LOOP
  !
  !----------------------------------------------------------------------

  do while( .not. convergence )
     !
     ! Initial residual: kryl(1) = L^-1 ( b - A x )
     !
     !
     ! L kryl(1) = A R^-1 x
     !
     call precon(&
          1_ip,nbvar,nbnodes,nrows,solve_sol(1)%kfl_symme,idprecon,ia,ja,an,&
          pn,invdiag,wa1,xx,kryl(1,1))

     do kk = 1,nrows
        kryl(kk,1) = wa3(kk) - kryl(kk,1)
     end do
     !
     !call couplings_check_dirichlet(1_ip,kryl(:,1),solve_sol(1))
     !
     ! resid = ||kryl(1)||
     !
     call solver_parallel_vector_L2norm(solve_sol(1),kryl(:,1),resid)

     resi2 = resi1
     resi1 = resid * invnb

     if( kfl_cvgso == 1 ) call outcso(an,bb,xx)

     if( resid <= stopcri ) then
        !
        ! The initial guess is the solution
        !
        convergence = .true.
     else
        !
        ! Initialize 1-st term of the rhs of hessenberg system
        !
        rs(1) = resid
        !
        ! Ortonormalize kryl(*,1)
        !
        resid = 1.0_rp / resid
#ifdef BLAS
        if( INOTMASTER ) call DSCAL(nrows,resid,kryl,1_ip)
#else
        do kk = 1,nrows
           kryl(kk,1) = kryl(kk,1) * resid
        end do
#endif
     end if

     jj   = 0
     idx  = -1
     fin2 = convergence
     !
     ! Inner loop. Restarted each kryldim iterations
     !
     do while( .not. fin2 ) 

        solve_sol(1) % iters = solve_sol(1) % iters + 1
        jj    = jj + 1
        jj1   = jj + 1
        idx   = idx + jj
        !
        !  L kryl(jj1) = A R^-1 kryl(jj)
        !
        call precon(&
             1_ip,nbvar,nbnodes,nrows,solve_sol(1)%kfl_symme,idprecon,ia,ja,an,&
             pn,invdiag,wa1,kryl(1,jj),kryl(1,jj1))
        !
        ! Orthogonalization
        ! For i= 1, j
        !     H(i,j) = <v_i, v_j1>
        !       v_j1 = v_j1 - H(i,j) * v_i
        !
        call cputim(time1)
        if( kfl_ortho == 0 ) then
           !
           ! Gram-Schmidt
           !
           call gmrort(nbvar,nbnodes,jj,jj1,kryl,hh(idx+1))
           do ii = 1,jj

              resid = hh(idx+ii) ! hh(i,j)
#ifdef BLAS
              if( INOTMASTER ) call DAXPY(nrows,-resid,kryl(1,ii),1_ip,kryl(1,jj1),1_ip)
#else
              do kk = 1,nrows
                 kryl(kk,jj1) = kryl(kk,jj1) - resid * kryl(kk,ii)
              end do
#endif
           end do
        else
           !
           ! Modified Gram-Schmidt
           !
           do ii = 1,jj
              call solver_parallel_scalar_product(solve_sol(1),kryl(:,ii),kryl(:,jj1),resid)

              hh(idx+ii) = resid
#ifdef BLAS
              if( INOTMASTER ) call DAXPY(nrows,-resid,kryl(1,ii),1_ip,kryl(1,jj1),1_ip)
#else
              do kk= 1,nrows
                 kryl(kk,jj1) = kryl(kk,jj1) - resid * kryl(kk,ii)
              end do
#endif
           end do
        end if

        call cputim(time2)
        solve_sol(1) % cputi(5) = solve_sol(1) % cputi(5) + time2 - time1
        !
        ! H(jj1,jj) = ||kryl(*,jj1)||
        !
        call solver_parallel_vector_L2norm(solve_sol(1),kryl(:,jj1),resid)
        !call norm2x(nbvar,kryl(1,jj1),resid)

        hh(idx+jj1) = resid

        if( resid == 0.0_rp ) then
           fin2        = .true.
           convergence = .true.
           idx         = idx - jj
           jj          = jj - 1
           solve_sol(1) % iters       = solve_sol(1) % iters - 1
           if( kfl_solve /= 0 ) write(lun_outso,*) '||kryl(*,jj1)|| = 0 ',jj1,solve_sol(1) % iters
           goto 10
        else
           !
           ! Ortonormalize kryl(*,jj1)
           !
           resid = 1.0_rp / resid
#ifdef BLAS
           if( INOTMASTER ) call DSCAL(nrows,resid,kryl(1,jj1),1_ip)
#else
           do kk = 1,nrows
              kryl(kk,jj1) = kryl(kk,jj1) * resid
           end do
#endif
        end if
        !
        ! Update factorization of H. Perform previous
        ! transformations on jj-th column of H
        !
        do ii = 1,jj-1
           kk         =   ii + 1
           resid      =   hh(idx+ii)
           hh(idx+ii) =   cc(ii) * resid + ss(ii) * hh(idx+kk)
           hh(idx+kk) = - ss(ii) * resid + cc(ii) * hh(idx+kk)
        end do
        gama = hh(idx+jj)*hh(idx+jj) + hh(idx+jj1)*hh(idx+jj1)
        gama = sqrt(gama)
        !
        ! if gama is zero then take any small
        ! value will affect only residual estimate
        !
        if( gama == 0.0_rp ) then
           gama = epsmac
           if( kfl_solve /= 0 ) write(lun_outso,*) 'gama==0.0 ',solve_sol(1) % iters
        end if
        !
        ! Get next plane rotation
        !
        gama      =   1.0_rp / gama
        cc(jj)     =   hh(idx+jj)  * gama
        ss(jj)     =   hh(idx+jj1) * gama
        hh(idx+jj) =   cc(jj) * hh(idx+jj) + ss(jj) * hh(idx+jj1)
        !
        ! Update the rhs of the LS problem
        !
        rs(jj1)    = - ss(jj) * rs(jj)
        rs(jj)     =   cc(jj) * rs(jj)
        !
        ! Convergence Test
        !
        resid      =   abs( rs(jj1) )

        resi2 = resi1
        resi1 = resid * invnb
        if( kfl_cvgso == 1 ) call outcso(an,bb,xx)

        if( resid <= stopcri ) then
           convergence = .true.
           fin2        = .true.
        else
           if( solve_sol(1) % iters >= maxiter ) then
              convergence = .true.
              fin2        = .true.
           else
              if( jj >= kryldim ) then
                 fin2 = .true.
              end if
           end if
        end if

10      continue
     end do
     
     !----------------------------------------------------------------------
     !
     ! END INNER LOOP
     !
     !----------------------------------------------------------------------
     !
     ! Compute y => Solve upper triangular system
     !
     do ii= jj,2,-1
        rs(ii) = rs(ii) / hh(idx+ii)
        resid = rs(ii)

        do kk = 1,ii-1
           rs(kk) = rs(kk) - hh(idx+kk) * resid
        end do

        idx = idx - ii
     end do

     if( hh(1) /= 0.0_rp ) rs(1) = rs(1) / hh(1)
     !
     ! Linear combination of kryl(*,jj)'s to get the solution.
     !
     do ii = 1,jj
        resid = rs(ii)
#ifdef BLAS
        if( INOTMASTER ) call DAXPY(nrows,resid,kryl(1,ii),1_ip,xx,1_ip)
#else
        do kk = 1,nrows
           xx(kk) = xx(kk) + resid * kryl(kk,ii)
        end do
        !print*,xx(204)
#endif
     end do

     !call couplings_check_dirichlet(1_ip,xx,solve_sol(1))

     !do ii = 1,nrows
     !   if( lninv_loc(ii) == 19  ) print*,'a=',lninv_loc(ii),xx(ii)
     !   if( lninv_loc(ii) == 141 ) print*,'a=',lninv_loc(ii),xx(ii)
     !   if( lninv_loc(ii) == 126 ) print*,'a=',lninv_loc(ii),xx(ii)
     !end do

     !
     ! Implicit Dirichlet/Dirichlet
     !
     if( mcoup > 0 ) then
        if( INOTMASTER .and. I_AM_INVOLVED_IN_A_COUPLING_TYPE(BETWEEN_SUBDOMAINS,DIRICHLET_IMPLICIT) ) then
           nullify(xx_tmp)
           allocate(xx_tmp(nrows))
           do ii = 1,nrows
              xx_tmp(ii) = xx(ii)
           end do
           !
           ! Interpolate new value. RHS should be modified as well
           ! and multiplied by its number of partitions (as this is the case
           ! for the matrix)
           !
           do icoup = 1,mcoup
              if(    I_AM_IN_COUPLING(icoup)                            .and. &
                   & coupling_type(icoup) % kind == BETWEEN_SUBDOMAINS  .and. &
                   & coupling_type(icoup) % what == DIRICHLET_IMPLICIT  ) then
                 call COU_INTERPOLATE_NODAL_VALUES_go(icoup,1_ip,xx,xx_tmp)
                 do kk = 1,coupling_type(icoup) % wet % npoin_wet
                    ii     = coupling_type(icoup) % wet % lpoin_wet(kk)
                    call runend('GMRPLS: NOT CODED')
                    !bb(ii) = xx(ii) * real(lneig(ii),rp)
                 end do
              end if
           end do
           !
           ! wa3 = L^-1 b should be recomputed for proper restart
           !
           call precon(&
                3_ip,nbvar,nbnodes,nrows,solve_sol(1) % kfl_symme,idprecon,ia,ja,an,&
                pn,invdiag,wa1,bb,wa3)
           deallocate(xx_tmp)
        end if
     end if

  end do

  !-----------------------------------------------------------------
  !
  !  END MAIN LOOP
  !
  !-----------------------------------------------------------------

20 continue
  !
  ! Save Krylov subspace
  !
  if( solve_sol(1) % kfl_save_krylov > 0 ) then
     do ii = 1,kryldim
        !jj = (ii-1)*kryldim+1
        !call solver_krylov_subspace_save(solve_sol(1),kryl(:,ii),hh(ii))
     end do
  end if
  !
  ! Last operations
  !
  call solope(&
       2_ip, nbvar, idprecon, dummr, an, dummr, ja, ia, &
       bb, xx , ierr, dummr, dummr, resi1, dummr, kryl, dummr, &
       dummr, dummr, invdiag, dummr )

  if( INOTMASTER .and. I_AM_INVOLVED_IN_A_COUPLING_TYPE(BETWEEN_SUBDOMAINS,DIRICHLET_IMPLICIT) ) then
     deallocate(lneig)
  end if

  call memory_deallo(memit,'INVDIAG','gmrpls',invdiag)
  call memory_deallo(memit,'KRYL'   ,'gmrpls',kryl   )
  call memory_deallo(memit,'WA3'    ,'gmrpls',wa3    )
  call memory_deallo(memit,'WA1'    ,'gmrpls',wa1    )
  call memory_deallo(memit,'HH'     ,'gmrpls',hh     )
  call memory_deallo(memit,'RS'     ,'gmrpls',rs     )
  call memory_deallo(memit,'SS'     ,'gmrpls',ss     )
  call memory_deallo(memit,'CC'     ,'gmrpls',cc     )

end subroutine gmrpls
 
