!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Algebraic_Solver
!> @{                                                                   
!> @file    precon.f90
!> @author  Guillaume Houzeaux
!> @date    15/07/2015
!> @brief   Apply a preconditioning to basic operations
!> @details Perform the following operations:
!>
!>          ITASK = 1 ... Solve system             L q = A R^-1 p 
!>          ITASK = 2 ... Solve system             L q = q 
!>          ITASK = 3 ... Solve system             L q = p 
!>          ITASK = 4 ... Solve system               q = A R^-1 p
!>          ITASK = 5 ... Recover final solution   R x = x'
!>          ITASK = 6 ... Compute initial solution   x = R x'
!>
!>          L = left preconditioner
!>          R = right preconditioner
!>
!>          INVDIAG ... D^-1
!>          WW ........ Working vector
!>          Q ......... Output vector
!>          P ......... Input vector
!>          AN ........ Matrix
!>          PN ........ Preconditioner Matrix
!>
!>          If there is a coarse grid preconditioner C
!>          Let P = L or P = R
!>          The preconditioner with coarse grid correction is:
!>          M^-1 = P^-1 + C^-1
!>
!>          At some point, we have to solve:
!>          M q = p so that
!>          q    = (P^-1 + C^-1) p = P^-1 p + C^-1 p
!>
!>          1. Compute coarse grid correction: p' = C^-1 p 
!>             then q = P^-1 p + p'
!>          2. Solve: P q1 = p => q1 = P^-1 p
!>          3. Update:  q  = q1 + p'
!>
!> @}                                                                   
!-----------------------------------------------------------------------

recursive subroutine precon(&
     itask,nbvar,npoin,nrows,kfl_symme,idprecon,ia,ja,an,&
     pn,invdiag,ww,pp,qq)
  
  use def_kintyp,         only  : ip,rp
  use def_master,         only  : INOTMASTER,npoi1,NPOIN_TYPE,zeror
  use def_domain,         only  : invpr_aii
  
  use def_solver,         only  : block_diagonal
  use def_solver,         only  : SOL_NO_PRECOND,SOL_SQUARE,SOL_DIAGONAL
  use def_solver,         only  : SOL_MATRIX,SOL_GAUSS_SEIDEL,SOL_BIDIAGONAL,SOL_LINELET
  use def_solver,         only  : SOL_ORTHOMIN,solve_sol,SOL_AII,SOL_RAS
  use def_solver,         only  : SOL_MULTIGRID,SOL_NEUMANN
  use def_solver,         only  : SOL_LEFT_PRECONDITIONING
  use def_solver,         only  : SOL_RIGHT_PRECONDITIONING
  use def_solver,         only  : SOL_BLOCK_DIAGONAL
  use mod_solver,         only  : solver_SpMV
  use mod_solver,         only  : solver_ras_right_preconditioning
  use mod_preconditioner, only  : preconditioner_gauss_seidel
  use mod_preconditioner, only  : preconditioner_bidiagonal_initialize
  
  use mod_matrix,         only  : matrix_CSR_SpMV
  use mod_matrix

  use def_coupli,         only  : RESIDUAL,UNKNOWN
  use def_coupli,         only  : BETWEEN_SUBDOMAINS 

  use mod_couplings,      only  : COU_PRESCRIBE_DIRICHLET_IN_MATRIX
  use mod_couplings,      only  : COU_INTERPOLATE_NODAL_VALUES
  use mod_couplings,      only  : couplings_impose_dirichlet

  implicit none
  integer(ip), intent(in)     :: itask             !< What to do
  integer(ip), intent(in)     :: nbvar             !< Number of dof per node                            
  integer(ip), intent(in)     :: npoin             !< Numnber of nodes    
  integer(ip), intent(in)     :: nrows             !< Totoal number of dof=nbvar*npoin                  
  integer(ip), intent(in)     :: kfl_symme         !< If matrix is assembled with symmetric structure   
  integer(ip), intent(in)     :: idprecon          !< Preconditioner ID                                 
  integer(ip), intent(in)     :: ia(*)             !< Matrix graph                                      
  integer(ip), intent(in)     :: ja(*)             !< Matrix graph                                      
  real(rp),    intent(in)     :: an(nbvar,nbvar,*) !< Matrix of system                                  
  real(rp),    intent(in)     :: pn(*)             !< Preconditioner matrix                             
  real(rp),    intent(in)     :: invdiag(*)        !< Inverse diagonal                                  
  real(rp),    intent(inout)  :: pp(*)             !< RHS vector                                        
  real(rp),    intent(inout)  :: ww(*)             !< Working array                                     
  real(rp),    intent(inout)  :: qq(*)             !< Preconditioned vector                             
  integer(ip)                 :: ii,jj,kk 
  integer(ip)                 :: iijj,iikk
  integer(ip)                 :: iin,jjn
  integer(ip)                 :: ki,kj
  real(rp)                    :: time1,time2
  real(rp)                    :: time3,time4
  real(rp),    pointer        :: w2(:),Ap(:)

  call cputim(time1)
  time3 = 0.0_rp
  time4 = 0.0_rp
  
  select case ( itask )

  case ( 1_ip ) 

     !------------------------------------------------------------------- 
     !
     ! Solve system L q = A R^-1 p
     !
     !------------------------------------------------------------------- 
     !
     ! Left: w = A p
     ! Right w = p
     !
     allocate(Ap(max(1_ip,nrows)))
     if( solve_sol(1) % kfl_leftr == SOL_LEFT_PRECONDITIONING ) then
        call cputim(time3)        
        call solver_SpMV(solve_sol(1),an,pp,Ap)
        call cputim(time4)
     else if( solve_sol(1) % kfl_leftr == SOL_RIGHT_PRECONDITIONING ) then
#ifdef BLAS
        if( INOTMASTER ) call DCOPY(nrows,pp,1_ip,Ap,1_ip) 
#else
        do ii = 1,nrows
           Ap(ii) = pp(ii)
        end do
#endif
     end if 
     !
     ! Preconditioning
     !
     call all_precon(&
          nbvar,npoin,nrows,kfl_symme,idprecon,ia,ja,an,&
          pn,invdiag,ww,Ap,qq)
     !
     ! Coupling
     !
     call couplings_impose_dirichlet(solve_sol(1),qq)
     !
     ! At this point:
     ! Left:  L^{-1} q = A p
     ! Right:        q = R^{-1} p => then to q = A (R^{-1} p)
     !
     ! Coarse grid correction
     !
     if( solve_sol(1) % kfl_coarse == 1 ) then
        call coafin(nbvar,1.0_rp,Ap,qq)
     end if
     !
     ! Right preconditioning
     !
     deallocate(Ap)
     if( solve_sol(1) % kfl_leftr == SOL_RIGHT_PRECONDITIONING ) then
#ifdef BLAS
        if( INOTMASTER ) call DCOPY(nrows,qq,1_ip,ww,1_ip) 
#else
        do ii = 1,nrows 
           ww(ii) = qq(ii)
        end do
#endif
        call cputim(time3)
        call solver_SpMV(solve_sol(1),an,ww,qq)
        call cputim(time4)
        
     end if

  case ( 2_ip )

     !------------------------------------------------------------------- 
     !
     ! Solve L q = q
     !
     !------------------------------------------------------------------- 

     if( solve_sol(1) % kfl_leftr == SOL_LEFT_PRECONDITIONING ) then
        !
        ! Left preconditioning
        !
        allocate(w2(max(1_ip,nrows)))
#ifdef BLAS
        if( INOTMASTER ) call DCOPY(nrows,qq,1_ip,w2,1_ip) 
#else
        do ii = 1,nrows
           w2(ii) = qq(ii)
        end do
#endif
        call all_precon(&
             nbvar,npoin,nrows,kfl_symme,idprecon,ia,ja,an,&
             pn,invdiag,ww,w2,qq)
        !
        ! Coupling
        !
        call couplings_impose_dirichlet(solve_sol(1),qq)
        !
        ! Coarse grid correction
        !
        if( solve_sol(1) % kfl_coarse == 1 ) then
           call coafin(nbvar,1.0_rp,w2,qq)
        end if
        deallocate(w2)
     end if

  case ( 3_ip )

     !------------------------------------------------------------------- 
     !
     ! Solve L q = p (same as 2 but input is p)
     !
     !------------------------------------------------------------------- 

     if(    idprecon == SOL_NO_PRECOND .or. &
          & solve_sol(1) % kfl_leftr == SOL_RIGHT_PRECONDITIONING ) then
        !
        ! No preconditioner or Right preconditioning: do not do anything
        !
#ifdef BLAS
        if( INOTMASTER ) call DCOPY(nrows,pp,1_ip,qq,1_ip) 
#else
        do ii = 1,nrows
           qq(ii) = pp(ii) 
        end do
#endif

     else if( solve_sol(1) % kfl_leftr == SOL_LEFT_PRECONDITIONING ) then
        !
        ! Left preconditioning
        !
        call all_precon(&
             nbvar,npoin,nrows,kfl_symme,idprecon,ia,ja,an,&
             pn,invdiag,ww,pp,qq)
        !
        ! Coupling
        !
        call couplings_impose_dirichlet(solve_sol(1),qq)
        !
        ! Coarse grid correction
        ! 
        if( solve_sol(1) % kfl_coarse == 1 ) then
           call coafin(nbvar,1.0_rp,pp,qq)
        end if

     end if

  case( 4_ip ) 

     !------------------------------------------------------------------- 
     !
     ! Compute q = A R^-1 p
     !
     !------------------------------------------------------------------- 

     if( solve_sol(1) % kfl_leftr == SOL_LEFT_PRECONDITIONING ) then
        !
        ! Left preconditioning: q = A p
        !
        call cputim(time3)
        call solver_SpMV(solve_sol(1),an,pp,qq)
        call cputim(time4)
     else
        !
        ! Right preconditioning
        ! R w2 = p, q = A w
        !
        call all_precon(&
             nbvar,npoin,nrows,kfl_symme,idprecon,ia,ja,an,&
             pn,invdiag,ww,pp,qq)
        !
        ! Coupling
        !
        call couplings_impose_dirichlet(solve_sol(1),qq)
        !
        ! Coarse grid correction
        !
        if( solve_sol(1) % kfl_coarse /= 0 ) then
           call coafin(nbvar,1.0_rp,pp,qq)
        end if
        !
        ! q = A w
        !
        if( INOTMASTER ) then
           allocate(w2(nrows))
#ifdef BLAS
           if( INOTMASTER ) call DCOPY(nrows,qq,1_ip,w2,1_ip) 
#else
           do ii = 1,nrows
              w2(ii) = qq(ii)
           end do
#endif
           call cputim(time3)
           call solver_SpMV(solve_sol(1),an,w2,qq)
           call cputim(time4)
           deallocate(w2)
        end if

     end if

  case ( 5_ip )

     !------------------------------------------------------------------- 
     !
     ! Recover final solution: q <= R^-1 q 
     !
     !------------------------------------------------------------------- 
     !
     ! Coarse grid correction: save q
     !
     if ( solve_sol(1) % kfl_leftr == SOL_RIGHT_PRECONDITIONING ) then
        allocate(w2(max(1_ip,nrows)))
#ifdef BLAS
        if( INOTMASTER ) call DCOPY(nrows,qq,1_ip,w2,1_ip) 
#else
        do ii = 1,nrows
           w2(ii) = qq(ii)
        end do
#endif
        call all_precon(&
             nbvar,npoin,nrows,kfl_symme,idprecon,ia,ja,an,&
             pn,invdiag,ww,w2,qq)
        !
        ! Coupling
        !
        call couplings_impose_dirichlet(solve_sol(1),qq)
        !
        ! Coarse grid correction
        !
        if( solve_sol(1) % kfl_coarse /= 0 ) then
           call coafin(nbvar,1.0_rp,w2,qq)
        end if
        deallocate(w2)
     end if

  case ( 6_ip )

     !------------------------------------------------------------------- 
     !
     ! Compute initial solution: q <= R q
     !
     !------------------------------------------------------------------- 

     if ( idprecon == SOL_SQUARE ) then 

        do ii = 1,nrows
           qq(ii) = qq(ii) / invdiag(ii)
        end do

     else if( solve_sol(1) % kfl_leftr == SOL_RIGHT_PRECONDITIONING ) then

        if ( idprecon == SOL_DIAGONAL ) then
           !
           ! Diagonal 
           !
           do ii = 1,nrows
              if( abs(invdiag(ii)) < zeror ) then
                 qq(ii) = qq(ii)!0.0_rp
              else
                 qq(ii) = qq(ii) / invdiag(ii)
              end if
           end do

        else if ( idprecon == SOL_ORTHOMIN ) then
           !
           ! Orthomin(1)
           !
           do ii = 1,nrows
              qq(ii) = qq(ii) / invdiag(ii) 
           end do

        else if ( idprecon == SOL_GAUSS_SEIDEL .and. INOTMASTER ) then
           !
           ! Gauss-Seidel
           !
           do ii = 1,nrows
              ww(ii) = qq(ii)
           end do
           call preconditioner_gauss_seidel(&
                2_ip,npoin,nbvar,an,ja,ia,invdiag,qq,ww) 

        else if ( idprecon == SOL_BIDIAGONAL .and. INOTMASTER ) then
           !
           ! Bidiagonal
           !
           do ii = 1,nrows
              ww(ii) = qq(ii)
           end do
           call preconditioner_bidiagonal_initialize(&
                2_ip,npoin,nbvar,invdiag,qq,ww)

        else if( idprecon == SOL_BLOCK_DIAGONAL .and. INOTMASTER ) then
           !
           ! Block diagonal
           !
           allocate( w2(nrows) ) 
#ifdef BLAS
           if( INOTMASTER ) call DCOPY(nrows,qq,1_ip,w2,1_ip) 
#else
           do ii = 1,nrows
              w2(ii) = qq(ii)
           end do
#endif
           do ii = 1,npoin        
              do jj = 1,nbvar
                 iijj = (ii-1) * nbvar + jj
                 qq(iijj) = 0.0_rp
                 do kk = 1,nbvar
                    iikk = (ii-1) * nbvar + kk
                    qq(iijj) = qq(iijj) + block_diagonal(jj,kk,ii) * w2(iikk)
                 end do
              end do
           end do
           deallocate( w2 )

        else if( idprecon == SOL_AII .and. INOTMASTER ) then
           !
           ! Block LU
           !  
           ! +-  -+   +-        -+ +-  -+    
           ! | qi |   | Aii  Aib | | qi |    
           ! |    | = |          | |    | 
           ! | qb |   |  0    D  | | qb |     
           ! +-  -+   +-        -+ +-  -+    
           !
           allocate( w2(nrows) )
           call matrix_CSR_SpMV(&
                1_ip,solve_sol(1) % direct_solver_block_LU % nn, nbvar, nbvar , nbvar, nbvar , &
                solve_sol(1) % direct_solver_block_LU % iA, solve_sol(1) % direct_solver_block_LU % jA   , &
                solve_sol(1) % direct_solver_block_LU % aa, qq , w2, invpr=invpr_aii)

           do ii = 1,npoi1
              do kk = ia(ii),ia(ii+1)-1
                 jj = ja(kk)
                 if( jj > npoi1 ) then  
                    do ki = 1,nbvar
                       iin = (ii-1) * nbvar + ki
                       do kj = 1,nbvar
                          jjn = (jj-1) * nbvar + kj
                          w2(iin) = w2(iin) + an(kj,ki,kk) * qq(jjn) 
                       end do
                    end do
                 end if
              end do
           end do
#ifdef BLAS
           if( INOTMASTER ) call DCOPY(nrows,w2,1_ip,qq,1_ip) 
#else
           do ii = 1,npoi1*nbvar
              qq(ii) = w2(ii)
           end do
#endif
           deallocate( w2 )

           do ii = npoi1*nbvar+1,npoin*nbvar
              qq(ii) = qq(ii) / invdiag(ii) 
           end do

        else if( idprecon == SOL_RAS ) then
           !
           ! RAS
           !
           call solver_ras_right_preconditioning(solve_sol(1),qq)

        end if

     end if
     !
     ! Coupling
     !
     call couplings_impose_dirichlet(solve_sol(1),qq)

  end select

  call cputim(time2)
  solve_sol(1) % cputi(3) = solve_sol(1) % cputi(3) + time2 + time3 - ( time1 + time4 )

end subroutine precon
