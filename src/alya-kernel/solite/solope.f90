!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Algebraic_Solver
!> @{
!> @file    solope.f90
!> @author  Guillaume Houzeaux
!> @date    15/07/2015
!> @brief   This routine computes some preliminary arrays for iterative solvers
!>          which solve L^-1 A x = L^-1 b or A R^-1 R x = b
!> @details
!>          \verbatim
!>
!>          ITASK == 1: Initial computations
!>          --------------------------------
!>          BB ........... RHS b
!>          NBNODES ...... Number of nodes
!>          XX ........... Initial solution
!>          PN ........... Preconditionning matrix
!>          IDPRECON ..... Preconditionner
!>          EPS .......... Solver tolerance
!>          AN ........... Matrix A
!>          IA, JA ....... Matrix graph in BCSR format
!>          INVDIAG ...... D^-1/2 or D^-1
!>          RR ........... b-Ax or L^-1 (b-Ax)
!>          ZZ ........... L^-1 (b-Ax)
!>          WW ........... Temporary array
!>          BP ........... Save first search vector to compute loss of
!>                         orthogonality
!>          INVNB ........ 1 / || L^-1 b ||
!>          NEWRHO ....... < L^-1 r , L^-1 r > or < L^-1 r , r >
!>          RESID ........ sqrt(NEWRHO)
!>          RESIN ........ RESID * INVNB
!>          IERRO = -1 ... Trivial solution x = 0
!>                = -2 ... Initial guess is solution
!>
!>          In the case of right preconditioning, solve:
!>          A R^-1 R x  = b
!>          Solve for x': A R^-1 x' = b
!>          so x'0 = R.x0
!>          and at the end recover x by solving R x = x'
!>
!>          ITASK == 2: End of solver
!>          --------------------------------
!>
!>          SOLVE_SOL(1) % RESF2 ... Final non-preconditionned residual
!>          Deallocate memory
!>
!>          \endverbatim
!>
!> @}
!-----------------------------------------------------------------------
subroutine solope(&
     itask, nbvar, idprecon, eps, an, pn, ja, ia, &
     bb, xx , ierro, stopcri, newrho, resid, invnb, rr, zz,&
     pp, ww, invdiag , bp )

  use def_kintyp,         only : ip,rp
  use def_elmtyp,         only : NOHOL
  use def_master,         only : IMASTER,zeror,INOTMASTER,npoi1
!  use def_master,         only : INOTSLAVE
  use def_solver,         only : memit,SOL_AII
  use def_solver,         only : SOL_SQUARE
  use def_solver,         only : SOL_GAUSS_SEIDEL,SOL_BIDIAGONAL,SOL_LINELET
  use def_solver,         only : SOL_MULTIGRID,SOL_RAS,SOL_SOLVER_CG
  use def_solver,         only : SOL_SOLVER_DEFLATED_CG
  use def_solver,         only : SOL_SOLVER_A_DEF2
  use def_solver,         only : SOL_SOLVER_PIPELINED_CG
  use def_solver,         only : SOL_NODES, SOL_EDGES
  use def_solver,         only : resin,resfi,resi1,resi2
  use def_solver,         only : solve_sol
  use def_solver,         only : SOL_SOLVER_DEFLATED_GMRES
  use def_solver,         only : SOL_SOLVER_GMRES
  use def_solver,         only : SOL_LEFT_PRECONDITIONING
  use def_solver,         only : SOL_BLOCK_DIAGONAL,block_diagonal,block_invdiagonal
  use def_solver,         only : SOL_SOLVER_BICGSTAB,SOL_SOLVER_PIPELINED_DEFLATED_CG
  use def_solver,         only : SOL_SOLVER_DEFLATED_BICGSTAB,SOL_MATRIX_HAS_CHANGED
  use def_solver,         only : SOL_SOLVER_MINRES
  use mod_solver,         only : solver_krylov_subspace_initial_guess
  use mod_solver,         only : solver_diagonal
  use mod_solver,         only : solver_condition_number
  use mod_solver,         only : solver_parallel_scalar_product
  use mod_solver,         only : solver_SpMV
  use mod_solver,         only : solver_parallel_vector_L2norm
  use mod_solver,         only : solver_impose_dirichlet_condition
  use def_domain,         only : lnoch,meshe,npoin,memor_dom
  use def_kermod,         only : ndivi
  use mod_memory,         only : memory_alloca         
  use mod_memory,         only : memory_deallo    
  use mod_matrix,         only : matrix_sparsify
  use mod_matrix,         only : matrix_copy_matrix
  use mod_matrix,         only : matrix_invdia
  use mod_matrix,         only : matrix_copy_matrix_block
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications, only : PAR_INTERFACE_MATRIX_EXCHANGE
  use mod_communications, only : PAR_INTERFACE_OWN_NODE_EXCHANGE
  use mod_communications, only : PAR_MAX
  use mod_maths,          only : maths_invert_matrix
  use def_kermod,         only : kfl_detection
  use mod_ker_detection,  only : ker_detection_maximum_residual
  use mod_graphs,         only : graphs_permut_metis_postordering
  use mod_graphs,         only : graphs_permut
  use mod_graphs,         only : graphs_number_along_vector
  use mod_preconditioner, only : preconditioner_gauss_seidel
  use mod_preconditioner, only : preconditioner_bidiagonal_initialize
  use mod_preconditioner, only : preconditioner_bidiagonal_arrays
  use mod_direct_solver,  only : direct_solver_factorization
  use mod_direct_solver,  only : direct_solver_partialcleaning
  use mod_direct_solver,  only : direct_allocate_temporary_matrix
  use mod_direct_solver,  only : direct_solver_cleaning
  use mod_direct_solver,  only : direct_solver_initialization
  use mod_direct_solver,  only : direct_solver_solution
  use mod_couplings,      only : COU_PRESCRIBE_DIRICHLET_IN_MATRIX
  use mod_couplings,      only : COU_INTERPOLATE_NODAL_VALUES
  use mod_linelet,        only : linelet_factorization
  use mod_postpr,         only : postpr_right_now
  use mod_deflated_cg,    only : matgro
  use mod_deflated_cg,    only : matgr2
  use mod_matrix,         only : matrix_chksym
  use mod_solver,         only : solver_ras_factorization

  implicit none
  integer(ip), intent(in)    :: itask
  integer(ip), intent(in)    :: nbvar
  integer(ip), intent(in)    :: idprecon
  real(rp),    intent(inout) :: eps
  real(rp),    intent(inout) :: an(nbvar,nbvar,*)
  real(rp),    intent(in)    :: pn(*)
  integer(ip), intent(in)    :: ja(*)
  integer(ip), intent(in)    :: ia(*)
  real(rp),    intent(in)    :: bb(*)
  real(rp),    intent(inout) :: xx(*)
  real(rp),    intent(inout) :: bp(*)
  integer(ip), intent(inout) :: ierro
  real(rp),    intent(out)   :: stopcri
  real(rp),    intent(out)   :: newrho
  real(rp),    intent(out)   :: resid
  real(rp),    intent(out)   :: invnb
  real(rp),    intent(out)   :: rr(*)
  real(rp),    intent(out)   :: zz(*)
  real(rp),    intent(out)   :: pp(*)
  real(rp),    intent(out)   :: ww(*)
  real(rp),    intent(out)   :: invdiag(*)
  integer(ip)                :: ii,jj,kk,izdom,ipoin,jpoin
  integer(ip)                :: nskyl
  integer(ip)                :: max_residual_node
  integer(ip)                :: nbnodes,nrows,ncols
  real(rp)                   :: bnorm,xnorm,rnorm,dummr
  real(rp)                   :: time1,time2
!  real(rp)                   :: xsymm
  real(rp)                   :: max_residual
  real(rp),    save          :: invb2
  
  nbnodes = solve_sol(1) % nequa
  nrows   = solve_sol(1) % nunkn
  ncols   = solve_sol(1) % ncols * solve_sol(1) % ndofn
  if( IMASTER ) then
     nrows = 0        ! Master does not perform any loop
     ncols = 0
  end if

  select case ( itask )

  case ( 1_ip )

     !----------------------------------------------------------------------
     !
     ! Initial computations
     !
     !----------------------------------------------------------------------

     call cputim(time1)

     ierro                = 0
     resi1                = 1.0_rp
     resi2                = 1.0_rp
     resin                = 1.0_rp
     solve_sol(1) % iters = 0
     solve_sol(1) % resin = resin
     solve_sol(1) % xorth = 0.0_rp
     !
     call penali(nbvar,nbnodes,ia,ja,an,bb,xx)

     !----------------------------------------------------------------------
     !
     ! Check symmetry of the matrix
     !
     !----------------------------------------------------------------------

!!$     xsymm = 0.0_rp
!!$     if( INOTMASTER ) then
!!$        call matrix_chksym(nbnodes,nbvar,ia,ja,an,xsymm,memit)
!!$     end if
!!$     xsymm=abs(xsymm)
!!$     call PAR_MAX(xsymm)
!!$     if( xsymm > 1.0e-12_rp ) then
!!$        if( INOTSLAVE ) then
!!$           print*,' '
!!$           print*,' '
!!$           print*,'XSYMM=',xsymm
!!$           print*,' '
!!$        end if
!!$        call runend('O.K.!')
!!$     end if

     !----------------------------------------------------------------------
     !
     ! Matrix comes from Schur complement: Compute A3^-1 = diag(A3)^-1
     !
     !----------------------------------------------------------------------

     if( solve_sol(1) % kfl_schum == 1 .and. INOTMASTER .and. nbnodes > 0 ) then
        call matrix_invdia(&
             nbnodes,solve_sol(1) % ndofn_A3,solve_sol(1) % kfl_symme,ia,ja,&
             solve_sol(1) % A3,solve_sol(1) % invA3,memit)
     end if

     !----------------------------------------------------------------------
     !
     ! INVDIAG = D^-1/2 or D^-1
     !
     !----------------------------------------------------------------------

     call solver_diagonal(solve_sol,an,invdiag)

     if( idprecon == SOL_SQUARE ) then
        if( nbvar == 1 ) then
           do ii = 1,nbnodes
              invdiag(ii) = sqrt( abs(invdiag(ii)) )
           end do
        else
           do ii = 1,nbnodes
              jj = (ii-1) * nbvar
              do kk = 1,nbvar
                 invdiag(jj+kk) = sqrt( abs(invdiag(jj+kk)) )
              end do
           end do
        end if
     end if
     !
     ! Put zero inverse diagonal on hole nodes
     !
     if( solve_sol(1) % kfl_leftr == SOL_LEFT_PRECONDITIONING .and. solve_sol(1) % kfl_where == SOL_NODES ) then
        jj = 0
        do ii= 1,nbnodes
           if( lnoch(ii) == NOHOL ) then
              do kk = 1,nbvar
                 jj = jj + 1
                 invdiag(jj) = 0.0_rp
              end do
           else
              do kk = 1,nbvar
                 jj = jj + 1
                 if( abs(invdiag(jj)) < zeror ) then
                    invdiag(jj) = 0.0_rp ! OJO
                 else
                    invdiag(jj) = 1.0_rp / invdiag(jj)
                 end if
              end do
           end if
        end do
     else if ( solve_sol(1) % kfl_leftr == SOL_LEFT_PRECONDITIONING .and. solve_sol(1) % kfl_where == SOL_EDGES ) then
        jj = 0
        do ii= 1,nbnodes
           do kk = 1,nbvar
              jj = jj + 1
              if(invdiag(jj)/=0.0_rp) then
                 invdiag(jj) = 1.0_rp / invdiag(jj)
              else
                 invdiag(jj) = 1.0_rp
              end if
           end do
        end do
     else
        jj = 0
        do ii= 1,nbnodes
           do kk = 1,nbvar
              jj = jj + 1
              if( abs(invdiag(jj)) < zeror ) then
                 invdiag(jj) = 0.0_rp
              else
                 invdiag(jj) = 1.0_rp / invdiag(jj)
              end if
           end do
        end do
     end if

     !----------------------------------------------------------------------
     !
     ! Coarse solver: compute and factorize coarse matrix
     !
     !----------------------------------------------------------------------

     if( solve_sol(1) % kfl_coarse == 1 ) then

        nskyl = solve_sol(1) % nzgro * nbvar * nbvar
        call direct_allocate_temporary_matrix(solve_sol(1) % direct_solver_coarse)
        if( solve_sol(1) % kfl_symme == 1 ) then
           call matgro(solve_sol(1) % ngrou,nbnodes,nskyl,nbvar,ia,ja,an,solve_sol(1) % direct_solver_coarse % aa)
        else
           call matgr2(solve_sol(1) % ngrou,nbnodes,nskyl,nbvar,ia,ja,an,solve_sol(1) % direct_solver_coarse % aa)
        end if
        if( INOTMASTER ) then
           call direct_solver_factorization(solve_sol(1) % direct_solver_coarse)
        end if

     end if

     !----------------------------------------------------------------------
     !
     ! Factorize preconditioners
     !
     !----------------------------------------------------------------------

     if( solve_sol(1) % kfl_update_precond == 1 ) then

        if( idprecon == SOL_LINELET .or. ( idprecon == SOL_MULTIGRID .and. solve_sol(1) % kfl_defpr == SOL_LINELET ) ) then
           !
           ! Linelet
           !
           call linelet_factorization(solve_sol(1) % ndofn,an,invdiag)

        else if( idprecon == SOL_GAUSS_SEIDEL .and. solve_sol(1) % kfl_renumbered_gs == 1 .and. INOTMASTER ) then
           !
           ! Streamwise Gauss-Seidel
           !
           call graphs_number_along_vector(&
                meshe(ndivi),solve_sol(1) % vecto_gs,solve_sol(1) % angle_stream,solve_sol(1) % permr_gs,&
                solve_sol(1) % invpr_gs,solve_sol(1) % lgrou_gs,solve_sol(1) % kfl_fixno,solve_sol(1)%ngrou_gs,&
                memor=memor_dom)

        else if( idprecon == SOL_BIDIAGONAL .and. solve_sol(1) % kfl_renumbered_gs == 1 .and. INOTMASTER ) then
           !
           ! Streamwise Bidiagonal
           !
           !
           ! Memory allocation for diagonal column (idiag) and matrix sub-diagonal elements (adiag)
           !
           !call memory_alloca(memit,'adiag1','solope',solve_sol(1) % adiag1,nbvar,nbvar,npoin)
           !call memory_alloca(memit,'idiag1','solope',solve_sol(1) % idiag1,nbvar,npoi1)

           call graphs_number_along_vector(&
                meshe(ndivi),solve_sol(1) % vecto_gs,solve_sol(1) % angle_stream,solve_sol(1) % permr_gs,&
                solve_sol(1) % invpr_gs,solve_sol(1) % lgrou_gs,solve_sol(1) % kfl_fixno,solve_sol(1)%ngrou_gs,&
                memor=memor_dom)
           !
           ! Memory allocation for diagonal column (idiag) and matrix sub-diagonal elements (adiag)
           !
           call memory_alloca(memit,'SOLVE % ADIAG1','solope',solve_sol(1) % adiag1,nbvar,nbvar,npoin)
           call memory_alloca(memit,'SOLVE % IDIAG1','solope',solve_sol(1) % idiag1,nbvar,npoi1)
           !
           ! Initialize Preconditioner
           !
           call preconditioner_bidiagonal_arrays(npoin,nbvar,an,ja,ia)

        else if( idprecon == SOL_AII .and. INOTMASTER .and. &
             ( solve_sol(1) % kfl_assem == SOL_MATRIX_HAS_CHANGED .or. solve_sol(1) % kfl_clean_precond == 1 ) ) then
           !
           ! Block LU: complete LU for interior nodes
           !
           call direct_allocate_temporary_matrix(solve_sol(1) % direct_solver_Block_LU)
           !
           ! Copy interior matrix
           !
           call matrix_copy_matrix(&
                npoi1,nbvar,solve_sol(1) % ia,solve_sol(1) % ja,an,&
                & solve_sol(1) % direct_solver_Block_LU % aa,solve_sol(1) % direct_solver_Block_LU % ia,&
                & solve_sol(1) % direct_solver_Block_LU % ja)
           !
           ! Sparsify matrix before doing the factorization
           !
           if( solve_sol(1) % threshold >= 0.0_rp ) then
              !dummi = solve_sol(1) % direct_solver_Block_LU % nz
              !dummr = real(solve_sol(1) % direct_solver_Block_LU % fillin,rp)
              call direct_solver_cleaning(solve_sol(1) % direct_solver_Block_LU)
              call matrix_sparsify(&
                   solve_sol(1) % threshold,solve_sol(1) % direct_solver_Block_LU % nn,nbvar,&
                   solve_sol(1) % direct_solver_Block_LU % nz,solve_sol(1) % direct_solver_Block_LU % ia,&
                   solve_sol(1) % direct_solver_Block_LU % ja,solve_sol(1) % direct_solver_Block_LU % aa,&
                   memor_opt = solve_sol(1) % direct_solver_Block_LU % memor )
              call direct_solver_initialization(solve_sol(1) % direct_solver_Block_LU)
           end if
           call direct_solver_factorization(solve_sol(1) % direct_solver_Block_LU)

        else if( idprecon == SOL_RAS .and. INOTMASTER .and. &
             ( solve_sol(1) % kfl_assem == SOL_MATRIX_HAS_CHANGED .or. solve_sol(1) % kfl_clean_precond == 1 ) ) then
           !
           ! RAS
           !
           call solver_ras_factorization(solve_sol(1),an)

        else if( idprecon == SOL_MULTIGRID ) then
           !
           ! Deflated multigrid
           !
           call direct_allocate_temporary_matrix(solve_sol(1) % direct_solver_AMG)
           solve_sol(1) % kfl_defas = 1
           nskyl = nbvar * nbvar * solve_sol(1) % nzgro
           if( solve_sol(1) % kfl_symme == 1 ) then
              call matgro(solve_sol(1) % ngrou,nbnodes,nskyl,nbvar,ia,ja,an,solve_sol(1) % direct_solver_AMG % aa)
           else
              call matgr2(solve_sol(1) % ngrou,nbnodes,nskyl,nbvar,ia,ja,an,solve_sol(1) % direct_solver_AMG % aa)
           end if
           if( INOTMASTER ) then
              call direct_solver_factorization(solve_sol(1) % direct_solver_AMG)
           end if

        else if( idprecon == SOL_BLOCK_DIAGONAL .and. INOTMASTER ) then
           !
           ! Block diagonal
           !
           nullify( block_diagonal )
           nullify( block_invdiagonal )
           call memory_alloca(memit,'BLOCK_DIAGONAL'   ,'solope',block_diagonal,   nbvar,nbvar,nbnodes)
           call memory_alloca(memit,'BLOCK_IBVDIAGONAL','solope',block_invdiagonal,nbvar,nbvar,nbnodes)
           do ipoin = 1,nbnodes
              izdom = ia(ipoin)
              do while( izdom <= ia(ipoin+1)-1 )
                 jpoin = ja(izdom)
                 if( ipoin == jpoin ) then
                    do jj = 1,nbvar
                       do ii = 1,nbvar
                          block_diagonal(ii,jj,ipoin) = an(jj,ii,izdom)
                       end do
                    end do
                    izdom = ia(ipoin+1)
                 end if
                 izdom = izdom + 1
              end do
           end do
           call PAR_INTERFACE_NODE_EXCHANGE(block_diagonal,'SUM','IN MY ZONE','SYNCHRONOUS')
           do ipoin = 1,nbnodes
              block_invdiagonal(1:nbvar,1:nbvar,ipoin) = block_diagonal(1:nbvar,1:nbvar,ipoin)
              call maths_invert_matrix(nbvar,block_invdiagonal(1:nbvar,1:nbvar,ipoin))
           end do

        end if

     end if

     !----------------------------------------------------------------------
     !
     ! INVNB = 1 / || L^-1 b ||
     !
     !----------------------------------------------------------------------
     !
     ! w = L^-1 b, BNORM = || L^-1 w ||
     !
     call precon(&
          3_ip,nbvar,nbnodes,nrows,solve_sol(1) % kfl_symme,idprecon,ia,ja,an,&
          pn,invdiag,rr,bb,ww)

     if(    solve_sol(1) % kfl_algso == SOL_SOLVER_CG                    .or. &
          & solve_sol(1) % kfl_algso == SOL_SOLVER_DEFLATED_CG           .or. &
          & solve_sol(1) % kfl_algso == SOL_SOLVER_PIPELINED_CG          .or. &
          & solve_sol(1) % kfl_algso == SOL_SOLVER_PIPELINED_DEFLATED_CG .or. &
          & solve_sol(1) % kfl_algso == SOL_SOLVER_MINRES                .or. &
          & solve_sol(1) % kfl_algso == SOL_SOLVER_A_DEF2                ) then
        !
        ! CG and DCG: normalize residual by || b || := (b , L^-1 b )^1/2
        !
        call solver_parallel_scalar_product(solve_sol(1),bb,ww,bnorm)
        if( bnorm < 0.0_rp ) call runend('SOLOPE: PRECONDITIONER M IS NOT SPD AS (b,M^-1 b) < 0 !')
        bnorm = sqrt(bnorm)

     else
        !
        ! Other solvers: normalize residual by || b || := (L^-1 b , L^-1 b )^1/2
        !
        call solver_parallel_vector_L2norm(solve_sol(1),ww,bnorm)

     end if
     !
     ! || L^-1 b || = 0: trivial solution x = 0
     !
     if( bnorm <= solve_sol(1) % bnorm_min .and. solve_sol(1) % solco >= 0.0_rp ) then
       do ii = 1,nrows
           xx(ii) = 0.0_rp
        end do
        resid                =  0.0_rp
        resin                =  0.0_rp
        resi1                =  0.0_rp
        invnb                =  0.0_rp
        invb2                =  0.0_rp
        solve_sol(1) % resin =  0.0_rp
        solve_sol(1) % resi2 =  0.0_rp
        ierro                = -1
        if( solve_sol(1) % kfl_cvgso == 1 ) call outcso(an,bb,xx)
        goto 10
     end if
     !
     ! User normalization BNORM=n_user. For left preconditioning, BNORM = ||n_user|| x ||L^-1||
     ! To compute ||L^-1||, we solve L x = 1 => x = L^-1.1 and ||L^-1||=||x|| which consists
     ! of computing the L^2 norm of the lumped inverse preconditioner matrix
     !
     if( solve_sol(1) % kfl_normalization /= 0 ) then
        bnorm = solve_sol(1) % normalization + zeror
        if( solve_sol(1) % kfl_leftr == SOL_LEFT_PRECONDITIONING ) then
           do ii = 1,nrows
              rr(ii) = 1.0_rp
           end do
           call precon(&
                2_ip,nbvar,nbnodes,nrows,solve_sol(1) % kfl_symme,idprecon,ia,ja,an,&
                pn,invdiag,ww,dummr,rr)
           call solver_parallel_vector_L2norm(solve_sol(1),rr,dummr)
           bnorm = max(dummr,zeror) * bnorm
        end if
     end if

     invnb = 1.0_rp / bnorm
     !
     ! GMRES needs L^-1 b
     !
     if(    solve_sol(1) % kfl_algso == SOL_SOLVER_GMRES .or. &
          & solve_sol(1) % kfl_algso == SOL_SOLVER_DEFLATED_GMRES ) then
        do ii= 1, nrows
           bp(ii) = ww(ii)
        end do
     end if
     !
     ! Initial x <= R x
     !
     !call precon(&
     !     6_ip,nbvar,nbnodes,nrows,solve_sol(1) % kfl_symme,idprecon,ia,ja,an,&
     !     pn,invdiag,ww,dummr,xx)

     !----------------------------------------------------------------------
     !
     ! Initial guess using previous Krylov supspace
     !
     !----------------------------------------------------------------------

     if(  solve_sol(1) % kfl_save_krylov > 0 .and. &
          solve_sol(1) % krylov_size     > 0 ) then
        !
        ! R = residual
        ! 
        if( solve_sol(1) % kfl_schum == 1 ) then
           call bcsrax_schur( 1_ip, nbnodes, nbvar, solve_sol(1) % ndofn_A3 , &
                solve_sol(1) % A1, solve_sol(1) % A2, solve_sol(1) % invA3, solve_sol(1) % A4,&
                invdiag, ja, ia, xx, rr )
        else
           call solver_SpMV(solve_sol(1),an,xx,rr)
        end if
        do ii = 1,nrows
           rr(ii) = bb(ii) - rr(ii)
        end do
        !
        ! X = Initial solution
        !
        call solver_krylov_subspace_initial_guess(solve_sol(1),xx,rr)

     end if

     !----------------------------------------------------------------------
     !
     ! Initial residual: r = b - Ax
     !
     !----------------------------------------------------------------------
     !
     ! XNORM = ||x||
     !
     call norm1x(nbvar,xx,xnorm)

     if( xnorm == 0.0_rp ) then
        !
        ! Initial x is zero: r = b
        !
        do ii = 1,nrows
           rr(ii) = bb(ii)
        end do

     else
        !
        ! r = b - A x
        !
        if( solve_sol(1) % kfl_schum == 1 ) then
           call bcsrax_schur( 1_ip, nbnodes, nbvar, solve_sol(1) % ndofn_A3 , &
                solve_sol(1) % A1, solve_sol(1) % A2, solve_sol(1) % invA3, solve_sol(1) % A4,&
                invdiag, ja, ia, xx, rr )
        else
           call solver_SpMV(solve_sol(1),an,xx,rr)
        end if

        call solver_parallel_vector_L2norm(solve_sol(1),rr,rnorm)

        do ii = 1,nrows
           rr(ii) = bb(ii) - rr(ii)
        end do
     end if

     !print*,'a=',rr(209)
     !do ii = 1,nrows
     !   if( solve_sol(1) % kfl_fixno(1,ii) > 0 ) rr(ii) = 0.0_rp
     !end do
     !print*,xx(204),rr(204),invdiag(204)
     !do ii = 1,nrows 
     !   if( lninv_loc(ii) == 19 ) print*,'RR 0=',rr(ii)
     !end do
     !call runend('O.K.!')
     !
     ! Residual = 0 on Dirichlet nodes
     !
     !if( solve_sol(1) % kfl_iffix > 0 ) then
     !   ii = 0
     !   do jj = 1,nbnodes
     !      do kk = 1,nbvar
     !         ii = ii + 1
     !         if( solve_sol(1) % kfl_fixno(kk,jj) > 0 ) rr(ii) = 0.0_rp
     !      end do
     !   end do
     !end if

     call solver_parallel_vector_L2norm(solve_sol(1),rr,rnorm)
     call solver_parallel_vector_L2norm(solve_sol(1),bb,invb2)

     solve_sol(1) % bnorm = invb2

     if( invb2 /= 0.0_rp ) then
        invb2 = 1.0_rp / invb2
     end if
     solve_sol(1) % resi2 = rnorm * invb2           ! Normalized residual
     
     !----------------------------------------------------------------------
     !
     ! Right preconditioning: initial x <= R x
     !
     !----------------------------------------------------------------------

     call precon(&
          6_ip,nbvar,nbnodes,nrows,solve_sol(1) % kfl_symme,idprecon,ia,ja,an,&
          pn,invdiag,ww,dummr,xx)

     !----------------------------------------------------------------------
     !
     ! Preconditioned residual: z = L^{-1} r
     !
     !----------------------------------------------------------------------

     if(  solve_sol(1) % kfl_algso == SOL_SOLVER_GMRES .or. &
          solve_sol(1) % kfl_algso == SOL_SOLVER_DEFLATED_GMRES ) then
        call precon(&
             2_ip,nbvar,nbnodes,nrows,solve_sol(1) % kfl_symme,idprecon,ia,ja,an,&
             pn,invdiag,ww,dummr,rr)
     else
        call precon(&
             3_ip,nbvar,nbnodes,nrows,solve_sol(1) % kfl_symme,idprecon,ia,ja,an,&
             pn,invdiag,ww,rr,zz)
     end if

     !----------------------------------------------------------------------
     !
     ! Initial rho: BCGSTAB: rho = < L^-1 r , L^-1 r >
     !              GMRES:   rho = < L^-1 r , L^-1 r >
     !              CG:      rho = < L^-1 r ,      r >
     !              MINRES:  rho = sqrt(< L^-1 r ,      r >) (gamma1)
     !
     !----------------------------------------------------------------------

     if(    solve_sol(1) % kfl_algso == SOL_SOLVER_BICGSTAB .or. &
          & solve_sol(1) % kfl_algso == SOL_SOLVER_DEFLATED_BICGSTAB ) then
        !
        ! BICGSTAB (r0'=r0)
        !
        do ii = 1,nrows
           rr(ii) = zz(ii)
        end do
        call solver_parallel_scalar_product(solve_sol(1),zz,zz,newrho)

     else if( solve_sol(1) % kfl_algso == SOL_SOLVER_GMRES .or. &
          &   solve_sol(1) % kfl_algso == SOL_SOLVER_DEFLATED_GMRES ) then
        !
        ! GMRES
        !
        call solver_parallel_scalar_product(solve_sol(1),rr,rr,newrho)
        
     else if( solve_sol(1) % kfl_algso == 11 ) then
        !
        ! RICHARDSON 
        !
        call solver_parallel_scalar_product(solve_sol(1),zz,zz,newrho)

     else
        ! 
        ! CG, DCG, pipelined CG, MINRES
        !
        call solver_parallel_scalar_product(solve_sol(1),rr,zz,newrho)
        if( newrho < 0.0_rp ) call runend('SOLOPE: PRECONDITIONER M IS NOT SPD AS (r,M^-1 r) < 0 !')

     end if

     resid                = sqrt(newrho)
     resin                = resid * invnb
     resi1                = resin
     solve_sol(1) % resin = resin
     !
     ! Convergence output. For Pipeline, initial residual included in algorithm
     !
     if( solve_sol(1) % kfl_cvgso == 1 ) then
        if(  solve_sol(1) % kfl_algso /= SOL_SOLVER_PIPELINED_CG          .and. &
             solve_sol(1) % kfl_algso /= SOL_SOLVER_PIPELINED_DEFLATED_CG ) then
           call outcso(an,bb,xx)
        end if
     end if
     !
     ! Adaptive residual
     !
     if( solve_sol(1) % kfl_adres == 1 ) then
        eps = max( solve_sol(1) % resin * solve_sol(1) % adres , solve_sol(1) % solmi )
     end if
     !
     ! stop criterion = ||b|| * eps
     !
     stopcri = eps * bnorm
     !
     ! Test if the initial guess is the solution
     !
     if ( resid <= stopcri ) then
        ierro = -2
        resfi =  resin
        goto 10
     end if
     !
     ! Initial direction: not needed by GMRES nor RICHARDSON nor MINRES
     !
     if(    solve_sol(1) % kfl_algso /= SOL_SOLVER_GMRES          .and. &
          & solve_sol(1) % kfl_algso /= SOL_SOLVER_DEFLATED_GMRES .and. &
          & solve_sol(1) % kfl_algso /= SOL_SOLVER_MINRES         ) then
        do ii = 1,nrows
           pp(ii) = zz(ii)
        end do
     end if
     !
     ! Save first direction to check orthogonality: CG type
     !
     if(  solve_sol(1) % kfl_algso == SOL_SOLVER_CG                    .or. &
          solve_sol(1) % kfl_algso == SOL_SOLVER_DEFLATED_CG           .or. &
          solve_sol(1) % kfl_algso == SOL_SOLVER_PIPELINED_CG          .or. &
          solve_sol(1) % kfl_algso == SOL_SOLVER_PIPELINED_DEFLATED_CG .or. &
          solve_sol(1) % kfl_algso == SOL_SOLVER_MINRES                .or. &
          solve_sol(1) % kfl_algso == SOL_SOLVER_A_DEF2                ) then
        do ii = 1,nrows
           bp(ii) = zz(ii)
        end do
     end if

     !----------------------------------------------------------------------
     !
     ! Compute condition number of the matrix
     !
     !----------------------------------------------------------------------

     if( solve_sol(1) % kappa > 0 ) then
        call solver_condition_number(&
             solve_sol(1),an,solve_sol(1) % kappa,&
             solve_sol(1) % lambda_min,solve_sol(1) % lambda_max)
     end if

10   continue

     call cputim(time2)
     solve_sol(1) % cputi(2) = solve_sol(1) % cputi(2) + time2 - time1

  case ( 2_ip )

     !----------------------------------------------------------------------
     !
     ! End of solver and non-preconditionned residual
     ! r = || b - Ax || / || b ||
     !
     !----------------------------------------------------------------------
     !
     ! Final R x = x'
     !
     if( ierro /= -1 ) then
        call precon(&
             5_ip,nbvar,nbnodes,nrows,solve_sol(1) % kfl_symme,idprecon,ia,ja,an,&
             pn,invdiag,rr,bb,xx)
     end if
     !
     ! Impose Dirichlet boundary conditions
     !
     call solver_impose_dirichlet_condition(solve_sol(1),xx)
     !
     ! Recover the solution on halo nodes
     !
     if( solve_sol(1) % kfl_full_rows == 1 .and. INOTMASTER ) then
        call PAR_INTERFACE_OWN_NODE_EXCHANGE(solve_sol(1) % ndofn,xx)
     end if
     !
     ! Final non-preconditioned residual
     !
     resfi = resi1
     solve_sol(1) % resfi = resi1

     if( solve_sol(1) % kfl_schur == 1 ) then
        !call par_schuax(&
        !     nbvar,zi,xb,aib,abi,abb,r_dom_aib,c_dom_aib,r_dom_abb,&
        !     c_dom_abb,r_dom_abi,c_dom_abi,askyl,iskyl,nskyl,rr)
     else
        if( solve_sol(1) % kfl_schum == 1 ) then
           call bcsrax_schur( 1_ip, nbnodes, nbvar, solve_sol(1) % ndofn_A3 , &
                solve_sol(1) % A1, solve_sol(1) % A2, solve_sol(1) % invA3, solve_sol(1) % A4,&
                invdiag, ja, ia, xx, rr )
        else
           call solver_SpMV(solve_sol(1),an,xx,rr)
        end if
     end if

     do ii = 1,nrows
        rr(ii) = bb(ii) - rr(ii)
     end do

     call solver_parallel_vector_L2norm(solve_sol(1),rr,rnorm)

     solve_sol(1) % resf2 = rnorm * invb2
     !
     ! Anomaly detection
     !
     if( solve_sol(1) % iters >= solve_sol(1) % miter .and. kfl_detection /= 0 ) then
        max_residual = -1.0_rp
        max_residual_node = 0
        do ii = 1,nbnodes
           rnorm = 0.0_rp
           jj = (ii-1)*nbvar
           do kk = 1,nbvar
              jj = jj + 1
              rnorm = rnorm + rr(jj) * rr(jj)
           end do
           if( rnorm > max_residual ) then
              max_residual = sqrt(rnorm)
              max_residual_node = ii
           end if
        end do
        call ker_detection_maximum_residual(&
             max_residual_node,max_residual,nbvar,xx,&
             'SOLVER_NOT_CONVERGED',trim(solve_sol(1) % wprob),'NO_COMMENT')
     end if

     if( solve_sol(1) % kfl_clean_precond == 1 ) then
        !
        ! Deallocate memory for Coarse solver
        !
        if( solve_sol(1) % kfl_coarse == 1 ) then
           call direct_solver_partialcleaning(solve_sol(1) % direct_solver_coarse)
        end if
        !
        ! Deallocate memory for block LU
        !
        if( idprecon == SOL_AII .and. INOTMASTER ) then
           call direct_solver_partialcleaning(solve_sol(1) % direct_solver_Block_LU)
        end if
        !
        if( idprecon == SOL_RAS .and. INOTMASTER .and. solve_sol(1) % kfl_clean_precond == 1 ) then
           !solve_sol(1) % num_solves == solve_sol(1) % num_multiple_solves ) then
           do ii = 1,size(solve_sol(1) % direct_solver_RAS,KIND=ip)
              call direct_solver_partialcleaning(solve_sol(1) % direct_solver_RAS(ii))
           end do
        end if
        !
        ! Deallocate memory for multigrid preconditioner
        !
        if( idprecon == SOL_MULTIGRID ) then
           call direct_solver_partialcleaning(solve_sol(1) % direct_solver_AMG)
        end if
        !
        ! Deallocate memory of streamwise bidiagonal preconditioner
        !
        if (idprecon == SOL_BIDIAGONAL .and. solve_sol(1) % kfl_renumbered_gs == 1 .and. INOTMASTER ) then
           !        call memory_deallo(memit, 'adiag1', 'solope', solve_sol(1) % adiag1)
           !        call memory_deallo(memit, 'idiag1', 'solope', solve_sol(1) % idiag1)
           call graphs_number_along_vector(&
                meshe(ndivi),solve_sol(1) % vecto_gs,solve_sol(1) % angle_stream,solve_sol(1) % permr_gs,&
                solve_sol(1) % invpr_gs,solve_sol(1) % lgrou_gs,solve_sol(1) % kfl_fixno,&
                solve_sol(1)%ngrou_gs,'DEALLOCATE',memor=memor_dom)

           call memory_deallo(memit, 'SOLVE % ADIAG1', 'solope', solve_sol(1) % adiag1)
           call memory_deallo(memit, 'SOLVE % IDIAG1', 'solope', solve_sol(1) % idiag1)

           !call preconditioner_bidiagonal_arrays(npoin,nbvar,an,ja,ia)

        end if
        !
        ! Deallocate memory of streamwise Gauss-Seidel preconditioner
        !
        if( idprecon == SOL_GAUSS_SEIDEL .and. solve_sol(1) % kfl_renumbered_gs == 1 .and. INOTMASTER ) then
           call graphs_number_along_vector(&
                meshe(ndivi),solve_sol(1) % vecto_gs,solve_sol(1) % angle_stream,solve_sol(1) % permr_gs,&
                solve_sol(1) % invpr_gs,solve_sol(1) % lgrou_gs,solve_sol(1) % kfl_fixno,&
                solve_sol(1)%ngrou_gs,'DEALLOCATE',memor=memor_dom)
        end if
     end if
     !
     ! Compute cos of angle between last p and p0 to check orthogonality
     ! <Ap,p0> / ( ||Ap|| ||p0|| ); zz is p0
     !
     if(  solve_sol(1) % kfl_algso == SOL_SOLVER_CG           .or. &
          solve_sol(1) % kfl_algso == SOL_SOLVER_DEFLATED_CG  .or. &
          solve_sol(1) % kfl_algso == SOL_SOLVER_PIPELINED_CG .or. &
          solve_sol(1) % kfl_algso == SOL_SOLVER_MINRES       .or. &
          solve_sol(1) % kfl_algso == SOL_SOLVER_A_DEF2       ) then
        if( solve_sol(1) % kfl_schum == 1 ) then
           call bcsrax_schur( 1_ip, nbnodes, nbvar, solve_sol(1) % ndofn_A3 , &
                solve_sol(1) % A1, solve_sol(1) % A2, solve_sol(1) % invA3, solve_sol(1) % A4,&
                invdiag, ja, ia, pp, rr )
        else
           call solver_SpMV(solve_sol(1),an,pp,rr)
        end if
        call cosixy(nbvar,nbnodes,rr,zz,solve_sol(1)%xorth)
     end if
     !
     ! Block diagonal
     !
     if( idprecon == SOL_BLOCK_DIAGONAL .and. INOTMASTER ) then
        call memory_deallo(memit,'BLOCK_DIAGONAL',   'solope',block_diagonal)
        call memory_deallo(memit,'BLOCK_INVDIAGONAL','solope',block_invdiagonal)
     end if
     !
     ! Matrix from a Schur
     !
     if( solve_sol(1) % kfl_schum == 1 .and. INOTMASTER ) then
        call memory_deallo(memit,'INVA3','solope',solve_sol(1) % invA3)
     end if

  case ( 3_ip )

     !----------------------------------------------------------------------
     !
     ! Initial guess using previous Krylov supspace
     !
     !----------------------------------------------------------------------

     if(  solve_sol(1) % kfl_save_krylov > 0 .and. &
          solve_sol(1) % krylov_size > 0     ) then
        !
        ! R = residual
        !
        if( solve_sol(1) % kfl_schum == 1 ) then
           call bcsrax_schur( 1_ip, nbnodes, nbvar, solve_sol(1) % ndofn_A3 , &
                solve_sol(1) % A1, solve_sol(1) % A2, solve_sol(1) % invA3, solve_sol(1) % A4,&
                invdiag, ja, ia, xx, rr )
        else
           call solver_SpMV(solve_sol(1),an,xx,rr)
        end if
        do ii = 1,nrows
           rr(ii) = bb(ii) - rr(ii)
        end do
        !
        ! X = Initial solution
        !
        call solver_krylov_subspace_initial_guess(solve_sol(1),xx,rr)

     end if

  end select

end subroutine solope
