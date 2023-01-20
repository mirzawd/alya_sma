!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Algebraic_Solver
!> @addtogroup Preconditioner
!> @{
!> @file    all_precon.f90
!> @author  Guillaume Houzeaux
!> @date    15/07/2015
!> @brief   Solve a preconditioned system
!> @details Solve M q = p where M is the preconditioner
!> @}
!-----------------------------------------------------------------------

subroutine all_precon(&
     nbvar,npoin,nrows,kfl_symme,idprecon,ia,ja,an,&
     pn,invdiag,ww,pp,qq)
  use def_kintyp,          only :  ip,rp
  use def_master,          only :  INOTMASTER,npoi1
  use def_master,          only :  NPOIN_TYPE
  use def_solver,          only :  block_invdiagonal
  use def_solver,          only :  SOL_NO_PRECOND,SOL_SQUARE,&
       &                           SOL_DIAGONAL,SOL_MATRIX,&
       &                           SOL_GAUSS_SEIDEL,SOL_BIDIAGONAL,SOL_LINELET,&
       &                           SOL_ORTHOMIN,solve_sol,&
       &                           SOL_AII,SOL_RAS,&
       &                           SOL_MULTIGRID,SOL_NEUMANN,&
       &                           SOL_LEFT_PRECONDITIONING,&
       &                           SOL_RIGHT_PRECONDITIONING,&
       &                           SOL_BLOCK_DIAGONAL
  use mod_csrdir,         only  :  CSR_LUsol
  use mod_preconditioner, only  :  preconditioner_gauss_seidel
  use mod_preconditioner, only  :  preconditioner_bidiagonal_initialize
  use mod_direct_solver,  only  :  direct_solver_solution
  use mod_parall,         only  :  PAR_COMM_MY_CODE_ARRAY
  use mod_linelet,        only  :  linelet_solution
  use mod_solver,         only  :  solver_SpMV
  use mod_solver,         only  :  solver_parallel_scalar_product

  implicit none
  integer(ip), intent(in)     :: nbvar               !< Number of dof per node
  integer(ip), intent(in)     :: npoin               !< Numnber of nodes
  integer(ip), intent(in)     :: nrows               !< Totoal number of dof=nbvar*npoin
  integer(ip), intent(in)     :: kfl_symme           !< If matrix is assembled with symmetric structure
  integer(ip), intent(in)     :: idprecon            !< Preconditioner ID
  integer(ip), intent(in)     :: ia(*)               !< Matrix graph
  integer(ip), intent(in)     :: ja(*)               !< Matrix graph
  real(rp),    intent(in)     :: an(nbvar,nbvar,*)   !< Matrix of system
  real(rp),    intent(in)     :: pn(*)               !< Preconditioner matrix
  real(rp),    intent(in)     :: invdiag(*)          !< Inverse diagonal
  real(rp),    intent(inout)  :: pp(*)               !< RHS vector
  real(rp),    intent(inout)  :: ww(*)               !< Working array
  real(rp),    intent(inout)  :: qq(*)               !< Preconditioned vector
  integer(ip)                 :: ii,jj,kk
  integer(ip)                 :: iijj,iikk,ki,kj
  integer(ip)                 :: iin,jjn,ntotn,nn,iras
  integer(ip)                 :: ll,igrou
  real(rp)                    :: numer,denom
  real(rp)                    :: alpha
  real(rp),    pointer        :: r0(:)
  real(rp),    pointer        :: d1r0(:)
  real(rp),    pointer        :: Ad1r0(:)
  real(rp),    pointer        :: pp_1(:)
  real(rp),    pointer        :: qq_1(:)

  if( idprecon == SOL_NO_PRECOND ) then

     !------------------------------------------------------------------
     !
     ! No preconditioning
     !
     !------------------------------------------------------------------

#ifdef BLAS
     if( INOTMASTER ) call DCOPY(nrows,pp,1_ip,qq,1_ip)
#else
     do ii = 1,nrows
        qq(ii) = pp(ii)
     end do
#endif

  else if( idprecon == SOL_DIAGONAL ) then

     !------------------------------------------------------------------
     !
     ! Diagonal preconditioning
     !
     !------------------------------------------------------------------

#ifdef BLAS
     if( INOTMASTER ) call DSBMV('L',nrows,0_ip,1.0_rp,invdiag,1_ip,pp,1_ip,0.0_rp,qq,1_ip)
#else
     do ii = 1,nrows
        qq(ii) = invdiag(ii) * pp(ii)
     end do
#endif

  else if( idprecon == SOL_BLOCK_DIAGONAL ) then

     !------------------------------------------------------------------
     !
     ! Block diagonal
     !
     !------------------------------------------------------------------

     do ii = 1,npoin
        do jj = 1,nbvar
           iijj = (ii-1) * nbvar + jj
           qq(iijj) = 0.0_rp
           do kk = 1,nbvar
              iikk = (ii-1) * nbvar + kk
              qq(iijj) = qq(iijj) + block_invdiagonal(jj,kk,ii) * pp(iikk)
           end do
        end do
     end do

  else if( idprecon == SOL_MATRIX ) then

     !------------------------------------------------------------------
     !
     ! Matrix PN preconditioning
     !
     !------------------------------------------------------------------

     call solver_SpMV(solve_sol(1),pn,pp,qq)

  else if( idprecon == SOL_LINELET ) then

     !------------------------------------------------------------------
     !
     ! Linelet
     !
     !------------------------------------------------------------------

#ifdef BLAS
     if( INOTMASTER ) call DCOPY(nrows,pp,1_ip,qq,1_ip)
#else
     do ii = 1,nrows
        qq(ii) = pp(ii)
     end do
#endif
     call linelet_solution(nbvar,qq,ww,invdiag)

  else if( idprecon == SOL_GAUSS_SEIDEL .and. INOTMASTER ) then

     !------------------------------------------------------------------
     !
     ! Gauss-Seidel
     !
     !------------------------------------------------------------------

     if( kfl_symme == 1 ) then
        call runend('PRECON: NOT CODED')
     else
        call preconditioner_gauss_seidel(&
             1_ip,npoin,nbvar,an,ja,ia,invdiag,qq,pp)
     end if

  else if( idprecon == SOL_BIDIAGONAL ) then

     !-------------------------------------------------------------------
     !
     ! Bidiagonal
     !
     !-------------------------------------------------------------------

     if( kfl_symme == 1 ) then
        call runend('PRECON: NOT CODED')
     else
        call preconditioner_bidiagonal_initialize(&
             1_ip,npoin,nbvar,invdiag,qq,pp)
     end if

  else if( idprecon == SOL_ORTHOMIN ) then

     !------------------------------------------------------------------
     !
     ! Orthomin(1)
     !
     !------------------------------------------------------------------
     !
     ! q1    = q0 + alpha * D^-1 ( p - A q0 )
     ! q0    = D-1 p
     ! r0    = p - A q0
     !
     !            ( r0 , A D^-1 r0 )
     ! alpha = -------------------------
     !         ( A D^-1 r0 , A D^-1 r0 )
     !
     allocate( r0(   max(1_ip,nrows)) )
     allocate( d1r0( max(1_ip,nrows)) )
     allocate( Ad1r0(max(1_ip,nrows)) )

     do ii = 1,nrows
        ww(ii) = invdiag(ii) * pp(ii)                              ! q0 = D-1 p
     end do

     call solver_SpMV(solve_sol(1),an,ww,r0)
     do ii = 1,nrows
        r0(ii)   = pp(ii) - r0(ii)                                 ! r0 = p - A q0
        d1r0(ii) = invdiag(ii) * r0(ii)                            ! D^-1 r0
     end do

     call solver_SpMV(solve_sol(1),an,d1r0,Ad1r0)                    ! A D^-1 r0 

     call solver_parallel_scalar_product(solve_sol(1),r0,   Ad1r0,numer)
     call solver_parallel_scalar_product(solve_sol(1),Ad1r0,Ad1r0,denom)
     
     if( denom == 0.0_rp ) then
        alpha = 1.0_rp
     else
        alpha = numer / denom
     end if

     do ii = 1,nrows
        qq(ii) = ww(ii) + alpha * d1r0(ii)
     end do

     deallocate(r0,d1r0,Ad1r0)

  else if( idprecon == SOL_AII .and. INOTMASTER ) then

     !------------------------------------------------------------------
     !
     ! Block LU
     !
     !------------------------------------------------------------------
     !
     ! Solve:
     !
     ! +-        -+ +-  -+    +-  -+
     ! | Aii  Aib | | qi |    | pi |
     ! |          | |    | =  |    | => Aii qi = pi - Aib D^{-1} pb
     ! |  0    D  | | qb |    | pb |
     ! +-        -+ +-  -+    +-  -+
     !
     ! NB: for symmetric problems, one should put block Aib to zero
     !
     ! qb = D^{-1} pb
     !
     do ii = npoi1*nbvar+1,npoin*nbvar
        qq(ii) = invdiag(ii) * pp(ii)
     end do

     if( solve_sol(1) % kfl_symeq == 1 ) then
        call direct_solver_solution(solve_sol(1) % direct_solver_block_LU,pp,qq)
     else
        !call direct_solver_solution(solve_sol(1) % direct_solver_block_LU,pp,qq)
        !return
        !
        ! ww <= qb = Aib qb
        ! NB: for symmetric problems, comment this loop
        !
        ntotn = npoi1*nbvar
        ww(1:ntotn) = 0.0_rp
        do ii = 1,npoi1
           do kk = ia(ii),ia(ii+1)-1
              jj = ja(kk)
              if( jj > npoi1 ) then
                 do ki = 1,nbvar
                    iin = (ii-1) * nbvar + ki
                    do kj = 1,nbvar
                       jjn = (jj-1) * nbvar + kj
                       ww(iin) = ww(iin) + an(kj,ki,kk) * qq(jjn)
                    end do
                 end do
              end if
           end do
        end do
        call pararr('SLX',NPOIN_TYPE,nrows,ww)
        !
        ! w = pi - Aib ( D^{-1} pb )
        !
        ww(1:ntotn) = pp(1:ntotn) - ww(1:ntotn)
        !
        ! Solve Aii qi = wi
        !
        call direct_solver_solution(solve_sol(1) % direct_solver_block_LU,ww,qq)
     end if

  else if( idprecon == SOL_RAS .and. INOTMASTER ) then

     !------------------------------------------------------------------
     !
     ! RAS
     !
     !------------------------------------------------------------------

     if( solve_sol(1) % kfl_block_ras == 0 ) then
        if( solve_sol(1) % num_subd_ras == 1 ) then
           call direct_solver_solution(solve_sol(1) % direct_solver_RAS(1),pp,qq)
        else
           do iras = 1,solve_sol(1) % num_subd_ras
              nn = solve_sol(1) % direct_solver_RAS(iras) % nn
              allocate(qq_1(nn*nbvar)) 
              allocate(pp_1(nn*nbvar))
              do ii = 1,nn
                 kk       = solve_sol(1) % permr_ras(iras) % l(ii)
                 do ll = 1,nbvar
                    pp_1((ii-1)*nbvar+ll) = pp((kk-1)*nbvar+ll)
                 end do
              end do
              call direct_solver_solution(solve_sol(1) % direct_solver_RAS(iras),pp_1,qq_1)
              do ii = 1,nn
                 kk       = solve_sol(1) % permr_ras(iras) % l(ii)
                 do ll = 1,nbvar
                    qq((kk-1)*nbvar+ll)   = qq_1((ii-1)*nbvar+ll)
                 end do
              end do
              deallocate(qq_1) 
              deallocate(pp_1)              
           end do
        end if
     else
        if( solve_sol(1) % num_subd_ras == 1 ) then
           allocate( qq_1(npoin),pp_1(npoin) )
           do ii = 1,nbvar
              do jj = 1,npoin
                 pp_1(jj) = pp( (jj-1)*nbvar+ii )
              end do
              call direct_solver_solution(solve_sol(1) % direct_solver_RAS(ii),pp_1,qq_1)
              do jj = 1,npoin
                 qq( (jj-1)*nbvar+ii ) = qq_1(jj)
              end do
           end do
           deallocate( qq_1,pp_1 )
        else
           iras = 0
           do igrou = 1,solve_sol(1) % num_subd_ras
              do ll = 1,nbvar
                 iras = iras + 1
                 nn   = solve_sol(1) % direct_solver_RAS(iras) % nn
                 allocate(qq_1(nn)) 
                 allocate(pp_1(nn))
                 do ii = 1,nn
                    kk       = solve_sol(1) % permr_ras(iras) % l(ii)
                    pp_1(ii) = pp((kk-1)*nbvar+ll)
                 end do
                 call direct_solver_solution(solve_sol(1) % direct_solver_RAS(iras),pp_1,qq_1)
                 do ii = 1,nn
                    kk                  = solve_sol(1) % permr_ras(iras) % l(ii)
                    qq((kk-1)*nbvar+ll) = qq_1(ii)
                 end do
                 deallocate(qq_1) 
                 deallocate(pp_1)              
              end do
           end do
        end if
     end if
     !
     ! Only one side is responsible
     !
     !do ii = solve_sol(1) % nequa_own*nbvar+1,solve_sol(1) % nequa*nbvar
     !  qq(ii) = 0.0_rp
     !end do
     !call pararr('SLX',NPOIN_TYPE,nrows,qq)
     !
     ! Perform average on interface nodes
     !
     call pararr('SLX',NPOIN_TYPE,nrows,qq)
     !
     ! RAS: perform average, do not do anything for AS to maintain symmetry
     !
     if( solve_sol(1) % kfl_add_schwarz == 0 ) then
        do ii = npoi1+1,npoin
           jj = (ii-1) * nbvar
           do kk = 1,nbvar
              jj     = jj + 1
              qq(jj) = qq(jj) / real(PAR_COMM_MY_CODE_ARRAY(1) % bound_multiplicity(ii-npoi1),rp)
           end do
        end do
     end if
     
  else if( idprecon == SOL_MULTIGRID ) then

     !------------------------------------------------------------------
     !
     ! Multigrid
     !
     !------------------------------------------------------------------

     call coajac(&
          nbvar,npoin,nrows,ia,ja,an,invdiag,pp,qq)

  else if( idprecon == SOL_NEUMANN ) then

     !------------------------------------------------------------------
     !
     ! Neumann preconditioner
     !
     !------------------------------------------------------------------
     !call matrix_neumann_preconditioner(&
     !     itask,nbvar,npoin,kfl_symme,ia,ja,an,invdiag,pp,s,w,memit,qq)

  end if

end subroutine all_precon
