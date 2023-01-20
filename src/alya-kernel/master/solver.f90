!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Algebraic_Solver
!> @{
!> @file    solver.f90
!> @author  Guillaume Houzeaux
!> @brief   Bridge to solver
!> @details This routine calls the solvers
!>
!>          For diagonal solve which uses vmass, amatr must NOT be modified
!>
!>          Residual RHSID
!>          --------------
!>
!>          solve_sol(1) % kfl_recov = 0 ... Do not recover local residual
!>                                   = 1 ... Recover local residual
!>                                   = 2 ... Residual is already global
!>
!>          Tolerance SOLCO
!>          ---------------
!>
!>          solve_sol(1) % kfl_adres = 0 ... Solver tolerance is given by user
!>                                   = 1 ... Solver tolerance is adaptive
!>
!>          CPU times computed
!>          ------------------
!>
!>          SOLVE_SOL(1) % CPUTI(1) ... All operations
!>          SOLVE_SOL(1) % CPUTI(2) ... Initial operations
!>          SOLVE_SOL(1) % CPUTI(3) ... Preconditioning
!>          SOLVE_SOL(1) % CPUTI(4) ... Coarse system solver for DCG
!>          SOLVE_SOL(1) % CPUTI(5) ... Orthonormalization for GMRES
!>          SOLVE_SOL(1) % CPUTI(6) ... Full row matrix copy and exchange
!>
!>          Preconditioner and coarse solver update
!>          ---------------------------------------
!>
!>          When a module has assembled a matrix, the option should be written:
!>          solve_sol(1) % kfl_assem = SOL_MATRIX_HAS_CHANGED
!>
!>          Then, after solving the problem, the solver initializes
!>          solve_sol(1) % kfl_assem to:
!>
!>          solve_sol(1) % kfl_assem = SOL_SYSTEM_HAS_BEEN_SOLVED
!>
!>          It is very important to tell the solver when the matrix has changed
!>          so that all the preconditioners or coarse solvers are recomputed and
!>          refactorized.
!>
!>          If you have a doubt, activate the following option:
!>
!>          PRECONDITIONER: ALWAYS_COMPUTE
!>
!>          Size of unknown UNKOW
!>          ---------------------
!>
!>          UNKNO(1:NZRHS) ... Solution is always required up to NZRHS
!>
!>          NZRHS ............ Size of unknown and RHS array
!>                             = NPOIN      * NDOFN
!>          NEQUA ............ Number of nodes to solve
!>                             = NPOIN      ........... partial rows
!>                             = NPOIN_OWN  ........... full rows
!>                             = NEQUA_HALO ........... number halo nodes
!>          NUNKN ............ Number of equation to solve
!>                             = NPOIN      * NDOFN ... partial rows
!>                             = NPOIN_OWN  * NDOFN ... full rows
!>          NCOLS ............ Number of columns required to perform SpMV
!>                             = NPOIN      * NDOFN ... partial rows
!>                             = NPOIN_HALO * NDOFN ... full rows
!>
!> @}
!-----------------------------------------------------------------------

subroutine solver(rhsid,unkno,amatr,pmatr)

  use def_kintyp,           only :  ip,rp,lg
  use def_master,           only :  INOTMASTER
  use def_master,           only :  NPOIN_TYPE,nul1r
  use def_master,	          only :  NEDGE_TYPE

  use def_domain,           only :  npoin,c_dom,r_dom,c_sym,r_sym
  use def_solver,           only :  SOL_EDGES, SOL_NODES
  use def_solver,           only :  solve_sol,memma,cpu_solve,memit
  use def_solver,           only :  SOL_SOLVER_RICHARDSON
  use def_solver,           only :  SOL_NO_SOLVER
  use def_solver,           only :  SOL_SOLVER_SPARSE_DIRECT
  use def_solver,           only :  SOL_SOLVER_MAPHYS_UNSYMMETRIC
  use def_solver,           only :  SOL_SOLVER_MAPHYS_SYMMETRIC
  use def_solver,           only :  SOL_SOLVER_MAPHYS_SPD
  use def_solver,           only :  SOL_SOLVER_PETSC
  use def_solver,           only :  SOL_SYSTEM_HAS_BEEN_SOLVED
  use def_solver,           only :  SOL_MATRIX_HAS_CHANGED

  use mod_couplings,        only :  COU_PRESCRIBE_DIRICHLET_IN_MATRIX
  use mod_couplings,        only :  COU_INTERPOLATE_NODAL_VALUES
  use mod_couplings,        only :  couplings_impose_dirichlet 

  use mod_parall,           only :  PAR_COMM_MY_CODE_ARRAY
  use mod_communications,   only :  PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications,   only :  PAR_BARRIER
  use mod_communications,   only :  PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE
  use mod_communications,   only :  PAR_INTERFACE_OWN_NODE_EXCHANGE
  use mod_communications,   only :  PAR_INTERFACE_EDGE_EXCHANGE

  use mod_matrix,           only :  matrix_blokgs
  use mod_matrix,           only :  matrix_partial_to_full_row_matrix
  use mod_memory,           only :  memory_alloca
  use mod_memory,           only :  memory_deallo
  use mod_memory,           only :  memory_size
  use mod_iofile,           only :  iofile
  use mod_direct_solver,    only :  direct_solver_factorization_solution_partialcleaning
  use mod_matrix,           only :  matrix_permute_and_copy_matrix
  use mod_ker_timeline,     only :  ker_timeline
  use mod_outfor,           only :  outfor

  use mod_graphs,           only :  graphs_csr_to_coo
  implicit none
  real(rp),    intent(inout)   :: unkno(*)
  real(rp),    intent(inout)   :: amatr(*)
  real(rp),    intent(in)      :: pmatr(*)
  real(rp),    intent(inout)   :: rhsid(*)
  real(rp)                     :: time1,time2
  integer(ip)                  :: icomp,ndofn
  integer(ip)                  :: kblok,nblok
  real(rp),    pointer         :: rhscp(:)
  real(rp),    pointer         :: aa_idofn(:)     ! Block Gauss-Seidel
  real(rp),    pointer         :: xx_idofn(:)     ! Block Gauss-Seidel
  real(rp),    pointer         :: bb_idofn(:)     ! Block Gauss-Seidel
  integer(ip)                  :: auxdof          ! Auxiliar variable to recober the amoung of original dof
  logical(lg)                  :: exchange_RHS

  real(rp),    pointer         :: amatr_ndof(:)
  real(rp),    pointer         :: faux(:)
  real(rp),    pointer         :: xaux(:)
  integer(ip), pointer         :: ia_ndof(:)
  integer(ip), pointer         :: ja_ndof(:)
  integer(ip)                  :: nz,ii
  real(rp),    pointer         :: unkno_cpy(:)
  real(rp),    pointer         :: amatr_full(:)
  real(rp),    pointer         :: unkno_full(:)
  integer(ip)                  :: nrows,ncols,nzrhs  

  nullify(rhscp)
  nullify(aa_idofn)
  nullify(xx_idofn)
  nullify(bb_idofn)

  nullify(amatr_ndof)
  nullify(faux)
  nullify(xaux)
  nullify(ia_ndof)
  nullify(ja_ndof)

  nullify(unkno_cpy)

  nullify(unkno_full)
  nullify(amatr_full)

  !----------------------------------------------------------------------
  !
  ! Define some local flags
  !
  !----------------------------------------------------------------------
  
  exchange_RHS = .true.
  if(  solve_sol(1) % kfl_algso == SOL_SOLVER_MAPHYS_UNSYMMETRIC .or. &
       solve_sol(1) % kfl_algso == SOL_SOLVER_MAPHYS_SYMMETRIC   .or. &
       solve_sol(1) % kfl_algso == SOL_SOLVER_MAPHYS_SPD         .or. &
       solve_sol(1) % kfl_algso == SOL_SOLVER_PETSC ) then
     exchange_RHS = .false.
  end if
  if( solve_sol(1) % kfl_exchange == 0 ) exchange_RHS = .false.

  if( solve_sol(1) % nunkn > 0 ) then

     !----------------------------------------------------------------------
     !
     ! Dirichlet bc.
     !
     !----------------------------------------------------------------------

     if( solve_sol(1) % kfl_version == 0 ) then
        if( solve_sol(1) % kfl_algso == SOL_SOLVER_RICHARDSON ) then
           call presol(4_ip,solve_sol(1) % ndofn,solve_sol(1) % ndofn,rhsid,amatr,unkno)
        else
           call presol(1_ip,solve_sol(1) % ndofn,solve_sol(1) % ndofn,rhsid,amatr,unkno)
        end if
     end if

     !----------------------------------------------------------------------
     !
     ! Coupling with dirichlet
     !
     !----------------------------------------------------------------------

     call COU_PRESCRIBE_DIRICHLET_IN_MATRIX(solve_sol(1) % ndofn,npoin,r_dom,c_dom,amatr,rhsid,unkno)

     !----------------------------------------------------------------------
     !
     ! Limiting step
     !
     !----------------------------------------------------------------------

     if( solve_sol(1) % kfl_symme == 0 ) then
        call limite(solve_sol(1) % ndofn,r_dom,c_dom,amatr)
     else
        call limite(solve_sol(1) % ndofn,r_sym,c_sym,amatr)
     end if

     !----------------------------------------------------------------------
     !
     ! Others
     !
     !----------------------------------------------------------------------
     !
     ! Copy of RHS: useful if local (to each slave) RHS is
     ! needed after the solver
     !
     if( INOTMASTER .and. solve_sol(1) % kfl_recov == 1 ) then
        call memory_alloca(memma,'RHSCP','solver',rhscp,solve_sol(1) % nunkn)
        do icomp = 1,solve_sol(1) % nunkn
           rhscp(icomp) = rhsid(icomp)
        end do
     end if

  end if

  !----------------------------------------------------------------------
  !
  ! Determine of preconditioner and coarse solvers should be recomputed
  !
  !----------------------------------------------------------------------

  solve_sol(1) % kfl_update_precond = 1
  if( solve_sol(1) % kfl_clean_precond == 0 .and. solve_sol(1) % nsolv >= 1 ) then
     if( solve_sol(1) % kfl_assem /= SOL_MATRIX_HAS_CHANGED ) then
        !
        ! Matrix is not constant
        !
        solve_sol(1) % kfl_update_precond = 0
     else
        !
        ! Preconditioner is maintaine constant
        !
        solve_sol(1) % kfl_update_precond = 0
     end if
  end if
  !
  ! Headers
  !
  if( solve_sol(1) % kfl_algso /= 9 ) then
     if( solve_sol(1) % heade == 0 ) then
        solve_sol(1) % heade = 1
        call outfor(-39_ip,0_ip,' ')
        if( solve_sol(1) % kfl_algso == 2 ) call outdef(1_ip)
     end if
     call outfor(5_ip,0_ip,' ')
  end if

  !----------------------------------------------------------------------
  !
  ! Coupling
  !
  !----------------------------------------------------------------------

  call couplings_impose_dirichlet(solve_sol(1),unkno)

  !----------------------------------------------------------------------
  !
  ! Zones: Put 1 on diagonal
  !
  !----------------------------------------------------------------------

  call soldod(solve_sol(1) % ndofn,amatr,rhsid,unkno)

  !----------------------------------------------------------------------
  !
  ! Full row strategy: exchange matrix for own nodes. Matrix AMATR will
  ! be overwritten and recovered after the solver.
  ! The unknown must also be overdimensioned: when doing matrix-vector
  ! product in order to take into account halo nodes and not overwrites
  ! the unknown
  !
  !----------------------------------------------------------------------

  if( solve_sol(1) % kfl_full_rows == 1 ) then
     !
     ! Copy matrix
     !
     if( INOTMASTER ) then
        call cputim(time1)
        ndofn = solve_sol(1) % ndofn
        nrows = solve_sol(1) % nunkn
        ncols = solve_sol(1) % ncols
        nzrhs = solve_sol(1) % nzrhs

        call memory_alloca(memit,'AMATR_FULL','solver',amatr_full,solve_sol(1) % nnz * ndofn * ndofn)
        call memory_alloca(memit,'UNKNO_FULL','solver',unkno_full,max(nrows,ncols))
        !
        ! Convert from partial to full row
        !
        call matrix_partial_to_full_row_matrix(&
             solve_sol(1) % nequa,ndofn,solve_sol(1) % ia,solve_sol(1) % ja,amatr,&
             solve_sol(1) % ia_full,solve_sol(1) % ja_full,amatr_full)
        !
        ! Unknown is copied up to own nodes
        !
        unkno_full(1:nrows) = unkno(1:nrows)
        !
        ! Exchange matrix
        !
        call PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE(ndofn,amatr,amatr_full,PAR_COMM_MY_CODE_ARRAY(1))
        call PAR_INTERFACE_OWN_NODE_EXCHANGE(ndofn,unkno_full)
        call cputim(time2)
        solve_sol(1) % cputi(6) = solve_sol(1) % cputi(6) + time2 - time1
     else
        call memory_alloca(memit,'AMATR_FULL','solver',amatr_full,1_ip)
        call memory_alloca(memit,'UNKNO_FULL','solver',unkno_full,1_ip)
     end if
  end if

  !----------------------------------------------------------------------
  !
  ! Modify RHS due to Parall service: RHSID
  ! MAPHYS requires non-assembled (local) RHS
  !
  !----------------------------------------------------------------------

  if( solve_sol(1) % kfl_recov /= 2 .and. exchange_RHS ) then
     if(      solve_sol(1) % kfl_where == SOL_NODES .and. npoin > 0 ) then
        call pararr('SOL',NPOIN_TYPE,solve_sol(1) % ndofn*npoin,rhsid)
     else if( solve_sol(1) % kfl_where == SOL_EDGES) then
        call PAR_INTERFACE_EDGE_EXCHANGE(solve_sol(1) % ndofn, rhsid, 'SUM', 'IN MY CODE')
     end if
  end if

  !----------------------------------------------------------------------
  !
  ! Initialize some variables
  !
  !----------------------------------------------------------------------

  solve_sol(1) % iters = 0

  !----------------------------------------------------------------------
  !
  ! Algebraic solver
  !
  !----------------------------------------------------------------------

  call PAR_BARRIER()
  call cputim(time1)
  call ker_timeline('INI_SOLVER_'//trim(solve_sol(1) % wprob))

  if( solve_sol(1) % kfl_blogs == 1 .and. solve_sol(1) % ndofn > 1 ) then
     !
     ! Block Gauss-Seidel
     !
     auxdof = solve_sol(1) % ndofn
     nblok  = abs(solve_sol(1) % nblok)
     if (nblok == solve_sol(1) % ndofn) then
        solve_sol(1) % nzmat = solve_sol(1) % nzmat / solve_sol(1) % ndof2
        solve_sol(1) % nzrhs = solve_sol(1) % nzrhs / solve_sol(1) % ndofn
        solve_sol(1) % ndofn = 1
        solve_sol(1) % ndof2 = solve_sol(1) % ndofn * solve_sol(1) % ndofn
     else
        solve_sol(1) % nzmat = solve_sol(1) % nzmat / (nblok**2)
        solve_sol(1) % nzrhs = solve_sol(1) % nzrhs / nblok
        solve_sol(1) % ndofn = solve_sol(1) % ndofn_per_block
        solve_sol(1) % ndof2 = solve_sol(1) % ndofn * solve_sol(1) % ndofn
     end if

     do kblok = 1,nblok

        if( INOTMASTER ) then

           if( solve_sol(1) % nblok < 0 ) then
              if( solve_sol(1) % kfl_full_rows == 0 ) then
                 call matrix_blokgs(&
                      1_ip,solve_sol(1) % nequa,solve_sol(1) % ndofn_per_block,nblok,kblok,r_dom,c_dom,amatr,rhsid,unkno,&
                      aa_idofn,bb_idofn,xx_idofn,memma, solve_sol(1) % lperm_block)
              else
                 call matrix_blokgs(&
                      1_ip,solve_sol(1) % nequa,solve_sol(1) % ndofn_per_block,nblok,kblok,r_dom,c_dom,amatr_full,rhsid,unkno,&
                      aa_idofn,bb_idofn,xx_idofn,memma, solve_sol(1) % lperm_block)
              end if
           else
              if( solve_sol(1) % kfl_full_rows == 0 ) then
                 call matrix_blokgs(&
                      1_ip,solve_sol(1) % nequa,solve_sol(1) % ndofn_per_block,nblok,kblok,r_dom,c_dom,amatr,rhsid,unkno,&
                      aa_idofn,bb_idofn,xx_idofn,memma)
              else
                 call matrix_blokgs(&
                      1_ip,solve_sol(1) % nequa,solve_sol(1) % ndofn_per_block,nblok,kblok,r_dom,c_dom,amatr_full,rhsid,unkno,&
                      aa_idofn,bb_idofn,xx_idofn,memma)
              end if
           end if

        else
           aa_idofn => nul1r
           bb_idofn => nul1r
           xx_idofn => nul1r
        end if

        if( solve_sol(1) % kfl_algso == -999 ) then
           !
           ! No solver defined
           !
           call runend('SOLVER: NO SOLVER HAS BEEN DEFINED '&
                //'FOR PROBLEM '//trim(solve_sol(1) % wprob))

        else if( solve_sol(1) % kfl_algso <= 0 .or. solve_sol(1) % kfl_algso == SOL_SOLVER_SPARSE_DIRECT ) then
           !
           ! Direct solver Alya
           !
           if( solve_sol(1) % kfl_full_rows == 0 ) then
              call direct_solver_factorization_solution_partialcleaning(solve_sol(1) % direct_solver,amatr,rhsid,unkno)
           else
              call direct_solver_factorization_solution_partialcleaning(solve_sol(1) % direct_solver,amatr_full,rhsid,unkno)
           end if
        else
           !
           ! Iterative solvers
           !
           call solite(bb_idofn,xx_idofn,aa_idofn,pmatr)

        end if


        if( INOTMASTER ) then

           if( solve_sol(1) % nblok < 0 ) then
              if( solve_sol(1) % kfl_full_rows == 0 ) then
                 call matrix_blokgs(&
                      2_ip,solve_sol(1) % nequa,solve_sol(1) % ndofn_per_block,nblok,kblok,r_dom,c_dom,amatr,rhsid,unkno,&
                      aa_idofn,bb_idofn,xx_idofn,memma, solve_sol(1) % lperm_block)
              else
                 call matrix_blokgs(&
                      2_ip,solve_sol(1) % nequa,solve_sol(1) % ndofn_per_block,nblok,kblok,r_dom,c_dom,amatr_full,rhsid,unkno,&
                      aa_idofn,bb_idofn,xx_idofn,memma, solve_sol(1) % lperm_block)
              end if
           else
              if( solve_sol(1) % kfl_full_rows == 0 ) then
                 call matrix_blokgs(&
                      2_ip,solve_sol(1) % nequa,solve_sol(1) % ndofn_per_block,nblok,kblok,r_dom,c_dom,amatr,rhsid,unkno,&
                      aa_idofn,bb_idofn,xx_idofn,memma)
              else
                 call matrix_blokgs(&
                      2_ip,solve_sol(1) % nequa,solve_sol(1) % ndofn_per_block,nblok,kblok,r_dom,c_dom,amatr_full,rhsid,unkno,&
                      aa_idofn,bb_idofn,xx_idofn,memma)
              end if
           end if
        end if

     end do

     if( INOTMASTER ) then
        if( solve_sol(1) % nblok < 0 ) then
           if( solve_sol(1) % kfl_full_rows == 0 ) then
              call matrix_blokgs(&
                   3_ip,solve_sol(1) % nequa,solve_sol(1) % ndofn_per_block,nblok,kblok,r_dom,c_dom,amatr,rhsid,unkno,&
                   aa_idofn,bb_idofn,xx_idofn,memma, solve_sol(1) % lperm_block)
           else
              call matrix_blokgs(&
                   3_ip,solve_sol(1) % nequa,solve_sol(1) % ndofn_per_block,nblok,kblok,r_dom,c_dom,amatr_full,rhsid,unkno,&
                   aa_idofn,bb_idofn,xx_idofn,memma, solve_sol(1) % lperm_block)
           end if
        else
           if( solve_sol(1) % kfl_full_rows == 0 ) then
              call matrix_blokgs(&
                   3_ip,solve_sol(1) % nequa,solve_sol(1) % ndofn_per_block,nblok,kblok,r_dom,c_dom,amatr,rhsid,unkno,&
                   aa_idofn,bb_idofn,xx_idofn,memma)
           else
              call matrix_blokgs(&
                   3_ip,solve_sol(1) % nequa,solve_sol(1) % ndofn_per_block,nblok,kblok,r_dom,c_dom,amatr_full,rhsid,unkno,&
                   aa_idofn,bb_idofn,xx_idofn,memma)
           end if
        end if
     end if

     if( nblok == solve_sol(1) % ndofn ) then
        solve_sol(1) % ndofn = nblok
        solve_sol(1) % ndof2 = solve_sol(1) % ndofn * solve_sol(1) % ndofn
        solve_sol(1) % nzmat = solve_sol(1) % nzmat * solve_sol(1) % ndof2
        solve_sol(1) % nzrhs = solve_sol(1) % nzrhs * solve_sol(1) % ndofn
     else
        solve_sol(1) % ndofn = auxdof
        solve_sol(1) % ndof2 = solve_sol(1) % ndofn * solve_sol(1) % ndofn
        solve_sol(1) % nzmat = solve_sol(1) % nzmat * (nblok**2)
        solve_sol(1) % nzrhs = solve_sol(1) % nzrhs * nblok
     end if

  else
     !
     ! Monolithic
     !
     if(      solve_sol(1) % kfl_algso == SOL_NO_SOLVER ) then
        !
        ! No solver defined
        !
        call runend('SOLVER: NO SOLVER HAS BEEN DEFINED FOR PROBLEM '//trim(solve_sol(1) % wprob))

     else if( solve_sol(1) % kfl_algso <= 0 .or. solve_sol(1) % kfl_algso == SOL_SOLVER_SPARSE_DIRECT ) then
        !
        ! Direct solver Alya
        !
        if( solve_sol(1) % kfl_full_rows == 0 ) then
           call direct_solver_factorization_solution_partialcleaning(solve_sol(1) % direct_solver,amatr,rhsid,unkno)
        else
           call direct_solver_factorization_solution_partialcleaning(solve_sol(1) % direct_solver,amatr_full,rhsid,unkno)
        end if

     else
        !
        ! Iterative solvers
        !
        if( solve_sol(1) % kfl_full_rows == 0 ) then
           call solite(rhsid,unkno,amatr,pmatr)
        else
           call solite(rhsid,unkno_full,amatr_full,pmatr)           
        end if
        
     end if


  end if

  !----------------------------------------------------------------------
  !
  ! Full row strategy: recover original matrix
  !
  !----------------------------------------------------------------------

  if( solve_sol(1) % kfl_full_rows == 1 ) then
     if( INOTMASTER ) then
        if( nzrhs > memory_size(unkno_full) ) call runend('SOLVER: WE ARE IN TROUBLE, WRONG SIZE')
        do ii = 1,nzrhs
           unkno(ii) = unkno_full(ii)
        end do
     end if
     call memory_deallo(memit,'UNKNO_FULL','solver',unkno_full)
     call memory_deallo(memit,'AMATR_FULL','solver',amatr_full)
  end if

  !----------------------------------------------------------------------
  !
  ! Finish solver: perdiocity, internal force
  !
  !----------------------------------------------------------------------

  if( solve_sol(1) % kfl_version == 0 ) then
     !
     ! Internal reaction
     !
     call presol(5_ip,solve_sol(1) % ndofn,solve_sol(1) % ndofn,rhsid,amatr,unkno)

  else if( INOTMASTER ) then

  end if
  !
  ! The system has been solved
  ! Option can be useful if we solve once more with the same matrix (factorization, etc.)
  !
  solve_sol(1) % kfl_assem = SOL_SYSTEM_HAS_BEEN_SOLVED
  !
  ! Cleaning was forced by someone outside
  !
  if( solve_sol(1) % kfl_clean_precond == -1 ) solve_sol(1) % kfl_clean_precond = 0
  !
  ! Recover local residual RHSID in parallel
  !
  if( INOTMASTER  .and. solve_sol(1) % kfl_recov == 1 ) then
     do icomp = 1,solve_sol(1) % nunkn
        rhsid(icomp) = rhscp(icomp)
     end do
     call memory_deallo(memma,'RHSCP','solver',rhscp)
  end if

  !----------------------------------------------------------------------
  !
  ! Output matrix in MatrixMarket format
  !
  !----------------------------------------------------------------------

  if( ( solve_sol(1) % kfl_marke == 1 .or. solve_sol(1) % kfl_marke == 2 ) .and. INOTMASTER ) then

     if( solve_sol(1) % ndofn == 1 ) then

        ia_ndof => r_dom
        ja_ndof => c_dom
        call output_matrix(rhsid,unkno,amatr,ia_ndof,ja_ndof)
        nullify(ia_ndof)
        nullify(ja_ndof)

     else

        nz  = r_dom(npoin+1) - 1

        allocate( amatr_ndof(solve_sol(1) % ndofn * solve_sol(1) % ndofn * nz) )
        allocate( ja_ndof   (solve_sol(1) % ndofn * solve_sol(1) % ndofn * nz) )
        allocate( ia_ndof   (solve_sol(1) % ndofn * npoin + 1) )
        allocate( faux      (solve_sol(1) % ndofn * npoin    ) )
        allocate( xaux      (solve_sol(1) % ndofn * npoin    ) )

        call matreorder(rhsid,unkno,amatr,r_dom,c_dom,solve_sol(1) % ndofn,npoin,1_ip,&  ! 1_ip means 1x,1y,2x,2y
             amatr_ndof,ja_ndof,ia_ndof,faux,xaux)

        call output_matrix(faux,xaux,amatr_ndof,ia_ndof,ja_ndof)

        deallocate( amatr_ndof )
        deallocate( ja_ndof )
        deallocate( ia_ndof )
        deallocate( faux )
        deallocate( xaux )

     end if

  end if

  call cputim(time2)
  !
  ! Update some variables
  !
  solve_sol(1) % num_solves = solve_sol(1) % num_solves + 1
  solve_sol(1) % nsolv      = solve_sol(1) % nsolv      + 1
  !
  ! Solver statistics
  !
  if( solve_sol(1) % iters < solve_sol(1) % itsol(1)) solve_sol(1) % itsol(1) = solve_sol(1) % iters
  if( solve_sol(1) % iters > solve_sol(1) % itsol(2)) solve_sol(1) % itsol(2) = solve_sol(1) % iters
  solve_sol(1) % itsol(3)   = solve_sol(1) % itsol(3) + solve_sol(1) % iters
  solve_sol(1) % cputi(1)   = solve_sol(1) % cputi(1) + (time2-time1)
  cpu_solve                 = cpu_solve               + (time2-time1)
  if( solve_sol(1) % num_solves >= solve_sol(1) % num_multiple_solves ) solve_sol(1) % num_solves = 0
  !
  ! Timeline
  !
  call ker_timeline('END_SOLVER_'//trim(solve_sol(1) % wprob),solve_sol(1) % iters)

end subroutine solver

subroutine output_matrix(rhsid,unkno,amatr,ia,ja)

  !----------------------------------------------------------------------
  !
  ! Output a matrix in different formats:
  ! matrix market
  ! ASCII format
  !
  !----------------------------------------------------------------------

  use def_kintyp,         only :  ip,rp
  use def_solver,         only :  solve_sol
  use def_master,         only :  INOTMASTER,ISLAVE,kfl_paral
  use def_master,         only :  lninv_loc,intost
  use mod_iofile,         only :  iofile
  use def_domain,         only :  npoin
  implicit none
  real(rp),    intent(in)      :: amatr(*)
  real(rp),    intent(in)      :: rhsid(*)
  real(rp),    intent(in)      :: unkno(*)
  integer(ip), intent(in)      :: ia(*),ja(*)

  integer(ip)                  :: izdom,jzdom,ipoin,jpoin
  character(150)               :: file_market_aa
  character(150)               :: file_market_bb
  character(150)               :: file_market_cc
  character(150)               :: file_market_permu

  integer(ip),save             :: ipass = 0_ip

  ipass = ipass+1_ip

  if( solve_sol(1) % kfl_marke == 1 .and. INOTMASTER ) then
     !
     ! Matrix market format: one file per slave
     !
     file_market_aa    = 'MATRIXMARKET_MATRIX'//trim(intost(kfl_paral))//'_'//trim(intost(ipass))//'.txt'
     file_market_bb    = 'MATRIXMARKET_RHS'   //trim(intost(kfl_paral))//'_'//trim(intost(ipass))//'.txt'
     file_market_cc    = 'MATRIXMARKET_UNKNO' //trim(intost(kfl_paral))//'_'//trim(intost(ipass))//'.txt'

     call iofile(0_ip,90_ip,file_market_aa,   'MATRIXMARKET_MATRIX')
     call iofile(0_ip,91_ip,file_market_bb,   'MATRIXMARKET_RHS')
     call iofile(0_ip,93_ip,file_market_cc,   'MATRIXMARKET_UNKNO')

     write(90,'(a)') '%%MatrixMarket matrix coordinate real general'
     write(91,'(a)') '%%MatrixMarket matrix array real general'
     write(93,'(a)') '%%MatrixMarket matrix array real general'

     write(90,'(3(1x,i12))') npoin * solve_sol(1) % ndofn, npoin * solve_sol(1) % ndofn, &
          solve_sol(1) % nzmat
     write(91,'(3(1x,i12))') npoin * solve_sol(1) % ndofn,1_ip
     write(93,'(3(1x,i12))') npoin * solve_sol(1) % ndofn,1_ip

     do ipoin = 1,npoin * solve_sol(1) % ndofn
        write(91,'(1(1x,i9),1x,e13.6)') ipoin,rhsid(ipoin)
        write(93,'(1(1x,i9),1x,e13.6)') ipoin,unkno(ipoin)
        do izdom = ia(ipoin),ia(ipoin+1)-1
           jpoin = ja(izdom)
           write(90,'(2(1x,i9),1x,e13.6)') ipoin,jpoin,amatr(izdom)
        end do
     end do

     call iofile(2_ip,90_ip,file_market_aa,   'MATRIXMARKET_MATRIX')
     call iofile(2_ip,91_ip,file_market_bb,   'MATRIXMARKET_RHS')
     call iofile(2_ip,93_ip,file_market_cc,   'MATRIXMARKET_UNKNO')

     if( ISLAVE ) then

        file_market_permu = 'MATRIXMARKET_PERMUTATION'//trim(intost(kfl_paral))//'_'//trim(intost(ipass))//'.txt'
        call iofile(0_ip,92_ip,file_market_permu,'MATRIXMARKET_PERMUTATION')
        write(92,'(a)') '%%MatrixMarket matrix array integer general'
        write(92,'(3(1x,i12))') npoin * solve_sol(1) % ndofn,1_ip

        if( solve_sol(1) % ndofn == 1 ) then
           do ipoin = 1,npoin * solve_sol(1) % ndofn
              write(92,'(2(1x,i9))') ipoin,lninv_loc(ipoin)
           end do
        else  ! ndofn /=1 is only ready in sequential else some equivalent to lninv_loc would be needed
           do ipoin = 1,npoin * solve_sol(1) % ndofn
              write(92,'(2(1x,i9))') ipoin,ipoin
           end do
        end if

        call iofile(2_ip,92_ip,file_market_permu,'MATRIXMARKET_PERMUTATION')

     end if

  else if( solve_sol(1) % kfl_marke == 2 .and. INOTMASTER ) then
     !
     ! Alya ASCII format
     !
     file_market_aa    = 'MATRIX_ASCII'//trim(intost(kfl_paral))//'.txt'
     call iofile(0_ip,90_ip,file_market_aa,'MATRIX_ASCII')

     do ipoin = 1,npoin * solve_sol(1) % ndofn
        do jpoin = 1,npoin * solve_sol(1) % ndofn
           izdom = ia(ipoin)
           jzdom = 0
           izdom_loop: do while( izdom <= ia(ipoin+1)-1 )
              if( ja(izdom) == jpoin ) then
                 jzdom = izdom
                 exit izdom_loop
              end if
              izdom = izdom + 1
           end do izdom_loop
           if( jzdom /= 0 ) then
              write(90,'(1x,e13.6)',advance='no') amatr(jzdom)
           else
              write(90,'(1x,e13.6)',advance='no') 0.0_rp
           end if
        end do
        write(90,*)
     end do

  end if

end subroutine output_matrix
