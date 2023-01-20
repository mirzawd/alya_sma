!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine soldef(itask)
  !-----------------------------------------------------------------------
  !****f* master/soldef
  ! NAME
  !    soldef
  ! DESCRIPTION
  !    This subroutine deos the following:
  !    ITASK = 0 ... Initialize the solver type
  !    ITASK = 1 ... Bridge between modules and parall service
  !    ITASK = 4 ... Define matrices, RHS and preconditioner sizes
  ! USES
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_domain
  use def_master
  use def_solver
  use def_inpout
  use def_kermod,        only : ndivi 
  use mod_solver,        only : solver_initialization
  use mod_solver,        only : solver_errors
  use mod_direct_solver, only : direct_solver_type_initialization  
  use mod_memory_basic,  only : memory_alloca
  use mod_memory,        only : memory_size
  use mod_parall,        only : commd
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ivari,nvari,cont,idime,i1,i2,nn
  integer(ip)             :: ipoin,pos,ndofn_per_block,ii,kblok,n1,n2

  if( itask <= 0_ip ) then

     !-------------------------------------------------------------------
     !
     ! Initialization
     !
     !-------------------------------------------------------------------

     if( itask < 0 ) then
        nvari = -itask
        allocate(momod(modul) % solve(nvari))                ! Allocate memory
        solve       => momod(modul) % solve
        solve_sol   => momod(modul) % solve
     end if

     do ivari = 1,size(solve_sol,KIND=ip)
        call solver_initialization(solve_sol(ivari))
        if( modul == ID_TEMPER ) &
             solve_sol(ivari) % kfl_version = 1
        if( modul == ID_SOLIDZ ) &
             solve_sol(ivari) % kfl_version = 1
        if( modul == ID_SOLFE2 ) &
             solve_sol(ivari) % kfl_version = 1
        if( modul == ID_NASTIN ) &
             solve_sol(ivari) % kfl_version = 1
        if( modul == ID_ALEFOR ) &
             solve_sol(ivari) % kfl_version = 1
        if( modul == ID_PARTIS ) &
             solve_sol(ivari) % kfl_version = 1
        if( modul == ID_NEUTRO ) &
             solve_sol(ivari) % kfl_version = 1
        if( modul == ID_GUSANO ) &
             solve_sol(ivari) % kfl_version = 1
    !    if( modul == ID_TURBUL ) &
    !         solve_sol(ivari) % kfl_version = 1
     end do

  else if( itask == 1_ip ) then

     !-------------------------------------------------------------------
     !
     ! Used for Parall service
     !
     !-------------------------------------------------------------------

     do ivari=1,size(solve_sol,KIND=ip)
        call iexcha(solve_sol(ivari) % nrhss)
        call iexcha(solve_sol(ivari) % kfl_algso)
        call iexcha(solve_sol(ivari) % kfl_symme)
        call iexcha(solve_sol(ivari) % kfl_symeq)
        call iexcha(solve_sol(ivari) % kfl_recov)
        call iexcha(solve_sol(ivari) % kfl_exchange) 
        call iexcha(solve_sol(ivari) % kfl_ortho)
        call iexcha(solve_sol(ivari) % kfl_limit)
        call iexcha(solve_sol(ivari) % kfl_adres)
        call iexcha(solve_sol(ivari) % kfl_cmplx)
        call iexcha(solve_sol(ivari) % kfl_assem)
        call iexcha(solve_sol(ivari) % kfl_schur)
        call iexcha(solve_sol(ivari) % kfl_schum)
        call iexcha(solve_sol(ivari) % miter)
        call iexcha(solve_sol(ivari) % nkryd)
        call iexcha(solve_sol(ivari) % kfl_coarse)
        call iexcha(solve_sol(ivari) % kfl_version)
        call iexcha(solve_sol(ivari) % kfl_format)
        call iexcha(solve_sol(ivari) % kfl_blogs)
        call iexcha(solve_sol(ivari) % nblok) 
        !
        ! Deflated CG
        !
        call iexcha(solve_sol(ivari) % ngrou)
        !call iexcha(solve_sol(ivari) % nskyl)
        !call iexcha(solve_sol(ivari) % icoml)
        call iexcha(solve_sol(ivari) % ifbop)
        call iexcha(solve_sol(ivari) % kfl_gathe)
        call iexcha(solve_sol(ivari) % kfl_defas)
        call iexcha(solve_sol(ivari) % kfl_defso)
        !call iexcha(solve_sol(ivari) % limpo)       ! List of imposed nodes (deflated CG)
        !call iexcha(solve_sol(ivari) % lgrou)       ! Group list for deflated CG
        !call iexcha(solve_sol(ivari) % iskyl)       !
        !call iexcha(solve_sol(ivari) % idiag)       !
        !
        ! Linelet
        !
        !call iexcha(solve_sol(ivari) % npntr)       ! Linelet: tolerance aspect ratio
        !call iexcha(solve_sol(ivari) % nlpntr       ! Linelet: # points in linelet
        !call iexcha(solve_sol(ivari) % nline        ! Linelet: # linelets
        call iexcha(solve_sol(ivari) % kfl_factl)    ! Linelet: Factorization flag
        !call iexcha(solve_sol(ivari) % lpntr(:)     ! Linelet: tridiag. to original sym. matrix
        !call iexcha(solve_sol(ivari) % lrenu(:)     ! Linelet: point position in tridiag.
        !call iexcha(solve_sol(ivari) % lrenup(:)    ! Linelet: inverse permutation
        !call iexcha(solve_sol(ivari) % limli(:)     ! Linelet: list of imposed node
        !call iexcha(solve_sol(ivari) % lline(:)     ! Linelet: Pointer for each linelet
        !call rexcha(solve_sol(ivari) % trima)       ! Linelet: Tridiagonal matrix
        call rexcha(solve_sol(ivari) % toler)        ! Linelet: tolerance aspect ratio
        call iexcha(solve_sol(ivari) % kfl_linty)    ! Linelet: linelet type
        !
        ! Renumbered Gauss-Seidel
        !
        call iexcha(solve_sol(ivari) % kfl_renumbered_gs) ! Renumbered Gauss-Seidel permutation
        call rexcha(solve_sol(ivari) % angle_stream)      ! Minimum angle for renumbering
        !
        ! Penalization
        !
        call iexcha(solve_sol(ivari) % kfl_penal)          ! Penalization: method
        call rexcha(solve_sol(ivari) % penal)              ! Penalization: parameter

        call rexcha(solve_sol(ivari) % relax)              ! Relaxation parameter
        
        call iexcha(solve_sol(ivari) % kfl_scpre)          ! Schur preconditioner solver: skyline
        call iexcha(solve_sol(ivari) % kfl_scaii)          ! Schur matrix solver: skyline

        call iexcha(solve_sol(ivari) % kfl_defpr)          ! Multigrid precond: smoother preconditioner

        call rexcha(solve_sol(ivari) % solco)              ! Solver tolerance
        call rexcha(solve_sol(ivari) % xdiag)              ! Coefficient for diagonal solver
        call rexcha(solve_sol(ivari) % resin)              !
        call rexcha(solve_sol(ivari) % resfi)              !
        call rexcha(solve_sol(ivari) % resi2)              !
        call rexcha(solve_sol(ivari) % resf2)              !
        call rexcha(solve_sol(ivari) % adres)              ! Adaptive residual
        call rexcha(solve_sol(ivari) % solmi)              ! Minimum solver tolerance
        call rexcha(solve_sol(ivari) % xorth)              ! Orthogonality
        call iexcha(solve_sol(ivari) % heade)              ! Header
        call iexcha(solve_sol(ivari) % itsol(1))           ! Max. # of solver iterations
        call iexcha(solve_sol(ivari) % itsol(2))           ! Min. # of solver iterations
        call iexcha(solve_sol(ivari) % itsol(3))           ! Ave. # of solver iterations
        call iexcha(solve_sol(ivari) % nsolv)              ! # of solves
        call iexcha(solve_sol(ivari) % kfl_cvgso)          ! Output Convergence flag
        call iexcha(solve_sol(ivari) % lun_cvgso)          ! Convergence unit
        call iexcha(solve_sol(ivari) % kfl_exres)          ! Exact residual
        call iexcha(solve_sol(ivari) % kfl_solve)          ! Output solve flag
        call iexcha(solve_sol(ivari) % lun_solve)          ! Output unit
        call iexcha(solve_sol(ivari) % lun_exsol)          ! External solver output unit
        call iexcha(solve_sol(ivari) % kfl_preco)          ! Preconditioner type
        call iexcha(solve_sol(ivari) % kfl_leftr)          ! Preconditioner left or right
        call iexcha(solve_sol(ivari) % itpre)              ! Preconditioner iteration
        call iexcha(solve_sol(ivari) % kfl_marke)          ! Output in market format
        call iexcha(solve_sol(ivari) % kfl_force)          ! Force continuity after solver (in Parallel)
        call iexcha(solve_sol(ivari) % kfl_where)          ! Where are the unknowns
        call iexcha(solve_sol(ivari) % kfl_clean_precond)  ! Preconditioner/coarse solver computation
        call iexcha(solve_sol(ivari) % kfl_normalization)  ! Residual normalization strategy
        call iexcha(solve_sol(ivari) % kfl_block_ras)      ! Full/block RAS
        call iexcha(solve_sol(ivari) % max_dof_ras)        ! Maximum number of dofs of RAS
        call iexcha(solve_sol(ivari) % num_subd_ras)       ! Number of RAS subdomains
        call rexcha(solve_sol(ivari) % normalization)      ! Residual normalization value
        call iexcha(solve_sol(ivari) % num_multiple_solves)! Multiple solves using same matrix
        call iexcha(solve_sol(ivari) % kfl_full_rows)      ! By default, matrix are partiall assembled
        call iexcha(solve_sol(ivari) % kfl_save_krylov)    ! Krylov subdpace saving
        call iexcha(solve_sol(ivari) % kfl_kappa)          ! Condition number
        call iexcha(solve_sol(ivari) % kfl_roe_correction) ! Round off error corrections
        call rexcha(solve_sol(ivari) % bnorm_min)          ! Minimum RHS norm
        call rexcha(solve_sol(ivari) % threshold)          ! Threshold for direct solvers
        call iexcha(solve_sol(ivari) % omp_schedule)       ! OMP schedule
        call iexcha(solve_sol(ivari) % omp_chunk_size)     ! OMP chunk size
        call iexcha(solve_sol(ivari) % omp_interface)      ! OMP interface chunk size (=0: do not use OMP)
        call iexcha(solve_sol(ivari) % kfl_iffix)          ! Fixity treated by solver
        call iexcha(solve_sol(ivari) % kfl_bvnat)          ! Natural b.c.
        call iexcha(solve_sol(ivari) % kfl_dirichlet)      ! OMP interface chunk size (=0: do not use OMP)
        call iexcha(solve_sol(ivari) % kfl_exp_method)     ! Explicit method
        call iexcha(solve_sol(ivari) % kfl_exp_order)      ! Explicit order of approx inverse
 
        !
        ! Block structure
        !
        call iexcha(solve_sol(ivari) % kfl_block)
        call iexcha(solve_sol(ivari) % num_blocks)
        call iexcha(solve_sol(ivari) % block_num)
        do ii = 1,size(solve_sol(ivari) % block_dimensions,KIND=ip)
           call iexcha(solve_sol(ivari) % block_dimensions(ii))
        end do
        call iexcha(solve_sol(ivari) % kfl_react)
        nullify(solve_sol(ivari) % reaction)
        nullify(solve_sol(ivari) % lpoin_block)
        nullify(solve_sol(ivari) % lpoin_reaction)
        !
        ! Dimension of A3
        !
        call iexcha(solve_sol(ivari) % ndofn_A3)     ! Dimension of A3 matrix
        !
        ! Characters
        !
        if( nparc > len(parch) ) call runend('SOLDEF: TOO MANY CHARACTERS')
        if( parii == 2 .and. IMASTER ) parch(nparc+1:nparc+50) = solve_sol(ivari) % conf_file(1:50)
        if( parii == 2 .and. ISLAVE  ) solve_sol(ivari) % conf_file(1:50) = parch(nparc+1:nparc+50)
        nparc = nparc + 50

     end do

  else if( itask == 4 ) then

     !-------------------------------------------------------------------
     !
     ! Solver dimensions
     !
     !-------------------------------------------------------------------

     do ivari = 1,size(solve_sol,KIND=ip)

        if( solve_sol(ivari) % kfl_algso == -999 ) then

           continue
           !call runend('SOLDEF: NO SOLVER HAS BEEN DEFINED '&
           !     //'FOR PROBLEM '//trim(solve_sol(ivari) % wprob))

        else
           !
           ! Detect possible errors
           !
           call solver_errors(solve_sol(ivari))
           !
           ! Initialize direct solvers
           !
           call direct_solver_type_initialization(solve_sol(ivari) % direct_solver_coarse    )
           call direct_solver_type_initialization(solve_sol(ivari) % direct_solver_block_LU  )
           call direct_solver_type_initialization(solve_sol(ivari) % direct_solver           )
           call direct_solver_type_initialization(solve_sol(ivari) % direct_solver_Deflation )
           !
           ! Full row does not make sense in sequential
           !
           if( ISEQUEN ) solve_sol(ivari) % kfl_full_rows = 0
           !
           ! Special case RAS: allocate memory before initializing
           !
           if( solve_sol(ivari) % kfl_block_ras == 0 ) then
              allocate(solve_sol(ivari) % direct_solver_RAS(1))
           else 
              allocate(solve_sol(ivari) % direct_solver_RAS(solve_sol(ivari) % ndofn))
           end if
           call direct_solver_type_initialization(solve_sol(ivari) % direct_solver_RAS)

           if( INOTMASTER ) then

              if(      solve_sol(ivari) % kfl_where == SOL_NODES ) then

                 if( solve_sol(ivari) % kfl_full_rows == 0 ) then
                    solve_sol(ivari) % nequa = meshe(ndivi) % npoin
                    solve_sol(ivari) % nunkn = solve_sol(ivari) % ndofn * solve_sol(ivari) % nequa
                    solve_sol(ivari) % ncols = solve_sol(ivari) % ndofn * solve_sol(ivari) % nequa
                    solve_sol(ivari) % nequ1 = meshe(ndivi) % npoi1
                    solve_sol(ivari) % nequ2 = meshe(ndivi) % npoi2
                    solve_sol(ivari) % nequ3 = meshe(ndivi) % npoi3
                 else
                    solve_sol(ivari) % nequa = meshe(ndivi) % npoin_own
                    solve_sol(ivari) % nunkn = solve_sol(ivari) % ndofn * solve_sol(ivari) % nequa
                    solve_sol(ivari) % ncols = solve_sol(ivari) % ndofn * meshe(ndivi) % npoin_halo
                    solve_sol(ivari) % nequ1 = meshe(ndivi) % npoi1
                    solve_sol(ivari) % nequ2 = meshe(ndivi) % npoi2
                    solve_sol(ivari) % nequ3 = meshe(ndivi) % npoi3
                 end if
                 solve_sol(ivari) % nequa_own  = meshe(ndivi) % npoin_own
                 solve_sol(ivari) % nequa_halo = meshe(ndivi) % npoin_halo

              else if( solve_sol(ivari) % kfl_where == SOL_EDGES ) then

                 solve_sol(ivari) % nequa      = meshe(ndivi) % nedge
                 solve_sol(ivari) % nequ1      = meshe(ndivi) % nedge
                 solve_sol(ivari) % nequ2      = meshe(ndivi) % nedg2
                 solve_sol(ivari) % nequ3      = meshe(ndivi) % nedg3
                 solve_sol(ivari) % nunkn      = solve_sol(ivari) % ndofn * solve_sol(ivari) % nequa
                 solve_sol(ivari) % ncols      = solve_sol(ivari) % ndofn * solve_sol(ivari) % nequa
                 solve_sol(ivari) % nequa_halo = meshe(ndivi) % nedge
                 !
                 ! This should be fixed
                 !
                 if( IPARALL ) then
                    !solve_sol(ivari) % nequ2     = meshe(ndivi) commd % nedg2
                    !solve_sol(ivari) % nequ3     = meshe(ndivi) commd % nedg3
                    !solve_sol(ivari) % nequa_own = meshe(ndivi) commd % nedg3
                    solve_sol(ivari) % nequ2     = commd % nedg2
                    solve_sol(ivari) % nequ3     = commd % nedg3
                    solve_sol(ivari) % nequa_own = commd % nedg3
                 else
                    solve_sol(ivari) % nequ2     = 0
                    solve_sol(ivari) % nequ3     = -1
                    solve_sol(ivari) % nequa_own = meshe(ndivi) % nedge
                 end if

              else if( solve_sol(ivari) % kfl_where == SOL_ELEMENTS ) then

                 solve_sol(ivari) % nequa = meshe(ndivi) % nelem
                 solve_sol(ivari) % nunkn = solve_sol(ivari) % ndofn * solve_sol(ivari) % nequa
                 solve_sol(ivari) % ncols = solve_sol(ivari) % ndofn * nelem_2
                 solve_sol(ivari) % nequ1 = meshe(ndivi) % nelem
                 solve_sol(ivari) % nequ2 = 0
                 solve_sol(ivari) % nequ3 = -1
                 solve_sol(ivari) % nequa_own = solve_sol(ivari) % nequa

              end if

              solve_sol(ivari) % ndof2 = solve_sol(ivari) % ndofn**2

              !----------------------------------------------------------
              !
              ! Matrix
              !
              !----------------------------------------------------------
              
              if( solve_sol(ivari) % kfl_algso == SOL_SOLVER_RICHARDSON ) then
                 !
                 ! Richardson solver
                 !
                 continue

              else if( solve_sol(ivari) % kfl_algso == SOL_SOLVER_EXPLICIT ) then
                 !
                 ! Explicit solver
                 !
                 continue

              else
                 !
                 ! Iterative solvers
                 !
                 if(    solve_sol(ivari) % kfl_format == SOL_CSR_FORMAT .or. &
                      & solve_sol(ivari) % kfl_format == SOL_COO_FORMAT ) then
                    !
                    ! CSR and COO formats
                    !
                    if( solve_sol(ivari) % kfl_symme == 1 ) then
                       if(      solve_sol(ivari) % kfl_where == SOL_NODES ) then
                          solve_sol(ivari) % nzmat = nzsym * solve_sol(ivari) % ndof2
                       else if( solve_sol(ivari) % kfl_where == SOL_EDGES ) then
                          call runend('NOT CODED')
                       else if( solve_sol(ivari) % kfl_where == SOL_ELEMENTS ) then
                          call runend('NOT CODED')
                       else
                          call runend('solve_sol%kfl_where?')
                       end if
                    else
                       if(      solve_sol(ivari) % kfl_where == SOL_NODES ) then
                          !if( solve_sol(ivari) % kfl_full_rows == 0 ) then
                          !   solve_sol(ivari) % nzmat = nzsol * solve_sol(ivari) % ndof2
                          !else
                          !   !
                          !   ! Partial matrix is assembled anyway
                          !   ! Just in case, take the max with the full row matrix dimensions
                          !   !
                          !   solve_sol(ivari) % nzmat     = nzsol     * solve_sol(ivari) % ndof2
                          !end if
                          solve_sol(ivari) % nzmat     = nzsol     * solve_sol(ivari) % ndof2
                          solve_sol(ivari) % nzmat_own = nzdom_own * solve_sol(ivari) % ndof2
                       else if( solve_sol(ivari) % kfl_where == SOL_EDGES ) then
                          solve_sol(ivari) % nzmat = memory_size(meshe(ndivi) % c_edg) * solve_sol(ivari) % ndof2
                       else if( solve_sol(ivari) % kfl_where == SOL_ELEMENTS ) then
                          solve_sol(ivari) % nzmat =  meshe(ndivi) % nzelm_2 * solve_sol(ivari) % ndof2
                       else
                          call runend('solve_sol%kfl_where')
                       end if
                    end if
                    
                 else if( solve_sol(ivari) % kfl_format == SOL_ELL_FORMAT ) then
                    !
                    ! ELL format
                    !
                    solve_sol(ivari) % nzmat = meshe(ndivi) % nzdom_ell * solve_sol(ivari) % ndof2
                    
                 end if
                 
              end if

              !----------------------------------------------------------
              !
              ! RHS
              !
              !----------------------------------------------------------

              if(      solve_sol(ivari) % kfl_where == SOL_NODES ) then
                 solve_sol(ivari) % nzrhs = solve_sol(ivari) % ndofn * meshe(ndivi) % npoin
              else if( solve_sol(ivari) % kfl_where == SOL_EDGES ) then
                 solve_sol(ivari) % nzrhs = solve_sol(ivari) % ndofn * meshe(ndivi) % nedge
              else if( solve_sol(ivari) % kfl_where == SOL_ELEMENTS ) then
                 solve_sol(ivari) % nzrhs = solve_sol(ivari) % ndofn * nelem_2
              end if

              !----------------------------------------------------------
              !
              ! Block matrix: dimensions of each block
              !
              !----------------------------------------------------------

              if( solve_sol(ivari) % kfl_block == 0 ) then 
                 do kblok = 1,solve(1) % num_blocks
                    solve_sol(ivari) % block_dimensions(kblok) = solve(1+kblok-1) % ndofn
                 end do
                 if( solve_sol(ivari) % num_blocks == 1 ) then
                    solve_sol(ivari) % block_array(1) % bvess     => solve_sol(ivari) % bvess
                    solve_sol(ivari) % block_array(1) % kfl_fixno => solve_sol(ivari) % kfl_fixno
                    solve_sol(ivari) % block_array(1) % bvnat     => solve_sol(ivari) % bvnat
                 end if
              end if
              
              !----------------------------------------------------------
              !
              ! Preconditioning
              !
              !----------------------------------------------------------

              if( solve_sol(ivari) % kfl_preco == SOL_MATRIX ) then
                 !
                 ! Matrix preconditioning
                 !
                 solve_sol(ivari) % nzpre = solve_sol(ivari) % ndof2*nzsol

              else if( solve_sol(ivari) % kfl_preco == SOL_MASS_MATRIX .or. solve_sol(ivari) % kfl_preco == SOL_CLOSE_MASS_MATRIX ) then
                 !
                 ! Diagonal mass
                 !
                 continue

              else if( solve_sol(ivari) % kfl_preco == SOL_LOCAL_DIAGONAL .or. solve_sol(ivari) % kfl_preco == SOL_AUGMENTED_DIAGONAL ) then
                 !
                 ! Diagonal preconditioning for Richardson iterations
                 !
                 solve_sol(ivari) % nzpre = solve_sol(ivari) % ndofn * npoin

              end if
              !
              ! Only left preconditioning for conjugate gradient family
              !
              if(    solve_sol(ivari) % kfl_algso  == SOL_SOLVER_CG                   .or. &
                   & solve_sol(ivari) % kfl_algso  == SOL_DEFLATED_CG                 .or. &
                   & solve_sol(ivari) % kfl_algso  == SOL_SOLVER_A_DEF2               .or. &
                   & solve_sol(ivari) % kfl_algso  == SOL_GDECG                       .or. &
                   & solve_sol(ivari) % kfl_algso  == SOL_GCG                         .or. &
                   & solve_sol(ivari) % kfl_algso  == SOL_GCGNOPREC                         .or. &
                   & solve_sol(ivari) % kfl_algso  == SOL_SOLVER_PIPELINED_DEFLATED_CG ) then
                 solve_sol(ivari) % kfl_leftr = SOL_LEFT_PRECONDITIONING
              end if

              !----------------------------------------------------------
              !
              ! Block Gauss-Seidel
              !
              !----------------------------------------------------------

              if( solve_sol(1) % kfl_blogs == 1 .and. solve_sol(1) % ndofn > 1 ) then
                 !solve_sol(1) % nblok = abs(solve_sol(1) % nblok)
                 ! Check that ndofn / nblok is an integer
                 if (solve_sol(1) % nblok /= 0) then
                    if (mod(solve_sol(1) % ndofn,solve_sol(1) % nblok)/=0) then
                       call runend('SOLDEF: NDOF/NBLOK NOT AN INTEGER: WRONG NUMBER OF DOF PER BLOCK IN BGS')
                    end if
                 end if

                 if( solve_sol(1) % nblok == 0 ) then
                    !
                    ! One block per dof
                    !
                    solve_sol(1) % nblok = solve_sol(1) % ndofn
                    solve_sol(ivari) % ndofn_per_block = solve_sol(1) % ndofn / solve_sol(1) % nblok
                 else if( solve_sol(1) % nblok > 0 ) then
                    !
                    ! Contigous permutation
                    !
                    solve_sol(ivari) % ndofn_per_block = solve_sol(1) % ndofn / solve_sol(1) % nblok
                 else if( solve_sol(1) % nblok < 0 ) then
                    !
                    ! Uncontigous permutation
                    !
                    solve_sol(1) % ndofn_per_block = abs(solve_sol(1) % ndofn / solve_sol(1) % nblok)
                    allocate( solve_sol(1) % lperm_block(2,ndime*npoin) )
                    ndofn_per_block = -solve_sol(1) % ndofn / solve_sol(1) % nblok
                    allocate( solve_sol(1) % linvp_block(ndofn_per_block,nblok) )


                    ! lperm_block
                    do kblok=1,abs(solve_sol(1) % nblok)
                       cont = 0
                       do ipoin=1,npoin
                          do idime=1,ndime
                             pos = (ipoin-1)*ndime + idime
                             solve_sol(1) % lperm_block(kblok,pos) = idime + 2*cont*ndime + 2*(kblok-1)
                          end do
                          cont = cont+1
                       end do
                    end do

                    ! linvp_block
                 end if
              end if
              !
              ! Size of matrices and RHS
              !
             if( solve_sol(ivari) % kfl_cmplx == 0 ) then
                 nzmat = max( nzmat , solve_sol(ivari) % nzmat , solve_sol(ivari) % nzmat_own )                 
                 nzrhs = max( nzrhs , solve_sol(ivari) % nzrhs * solve_sol(ivari) % nrhss )
                 nzpre = max( nzpre , solve_sol(ivari) % nzpre )
              else 
                 nzmax = max( nzmax , solve_sol(ivari) % nzmat , solve_sol(ivari) % nzmat_own )
                 nzrhx = max( nzrhx , solve_sol(ivari) % nzrhs * solve_sol(ivari) % nrhss )
                 nzprx = max( nzprx , solve_sol(ivari) % nzpre )
              end if
              !
              ! Matrix has a block structure
              !
!!$              if( solve_sol(ivari) % block_number == 1 ) then
!!$                 !solve_sol(ivari) % block_pointers(1,1) = 1
!!$              else
!!$              end if
!!$                 last_pointer = 1
!!$                 kblok = 0
!!$                 allocate( block_size(solve_sol(ivari) % block_number * solve_sol(ivari) % block_number) )
!!$                 do kblok = 1,solve_sol(ivari) % block_number
!!$                    do jblok = 1,solve_sol(ivari) % block_number
!!$                       kblok = kblok + 1
!!$                       block_size(kblok) = &
!!$                            solve_sol(ivari) % block_dimensions(kblok) * solve_sol(ivari) % block_dimensions(jblok)
!!$                    end do
!!$                 end do
!!$                 ! block_poinetrs(:,:) or block_poinetrs(:) ???
!!$                 do kblok = 1,solve_sol(ivari) % block_number
!!$                    do jblok = 1,solve_sol(ivari) % block_number
!!$                       if( kblok == 1 .and. jblok == 1 ) then
!!$
!!$                       else
!!$
!!$                       end if
!!$                    end do
!!$                 end do
!!$                 deallocate( block_size )
!!$              end if

           end if
           !
           ! Name of solver
           !          
           if( solve_sol(ivari) % kfl_algso ==  SOL_SOLVER_DIRECT ) solve_sol(ivari) % wsolv = 'DIRECT LDU'
           if( solve_sol(ivari) % kfl_algso ==  SOL_SOLVER_PASTIX ) solve_sol(ivari) % wsolv = 'PASTIX - SPARSE DIRECT'
           if( solve_sol(ivari) % kfl_algso ==  SOL_SOLVER_MUMPS  ) solve_sol(ivari) % wsolv = 'MUMPS - SPARSE DIRECT'
           if( solve_sol(ivari) % kfl_algso ==  SOL_SOLVER_WSMP   ) solve_sol(ivari) % wsolv = 'WSMP - SPARSE DIRECT'
           if( solve_sol(ivari) % kfl_algso ==  SOL_SOLVER_PWSMP  ) solve_sol(ivari) % wsolv = 'PWSMP - SPARSE DIRECT'
           if( solve_sol(ivari) % kfl_algso ==  SOL_SOLVER_GPUQR  ) solve_sol(ivari) % wsolv = 'GPU - QR DIRECT SOLVER'
           if( solve_sol(ivari) % kfl_algso ==  SOL_SOLVER_CG     ) then
              if( solve_sol(ivari) % kfl_schur == 0 ) then
                 solve_sol(ivari) % wsolv = 'CONJUGATE GRADIENT'
              else
                 solve_sol(ivari) % wsolv = 'CONJUGATE GRADIENT FOR SCHUR COMPLEMENT'
              end if
           end if
           if( solve_sol(ivari) % kfl_algso == SOL_SOLVER_DEFLATED_CG           ) solve_sol(ivari) % wsolv = 'DEFLATED CONJUGATE GRADIENT'
           if( solve_sol(ivari) % kfl_algso == SOL_GDECG                        ) solve_sol(ivari) % wsolv = 'GPU - DEFLATED CONJUGATE GRADIENT'
           if( solve_sol(ivari) % kfl_algso == SOL_GGMRS                        ) solve_sol(ivari) % wsolv = 'GPU - GMRES'
           if( solve_sol(ivari) % kfl_algso == SOL_GCG                          ) solve_sol(ivari) % wsolv = 'GPU - CONJUGATE GRADIENT'
           if( solve_sol(ivari) % kfl_algso == SOL_GPCG                         ) solve_sol(ivari) % wsolv = 'GPU - PIPELINED CONJUGATE GRADIENT'
           if( solve_sol(ivari) % kfl_algso == SOL_GAMGX                        ) solve_sol(ivari) % wsolv = 'GPU - NVIDIA MULTIGRID'
           if( solve_sol(ivari) % kfl_algso == SOL_GCGNOPREC                    ) solve_sol(ivari) % wsolv = 'GPU - CONJUGATE GRADIENT NO PREC'
           if( solve_sol(ivari) % kfl_algso == SOL_SOLVER_BICGSTAB              ) solve_sol(ivari) % wsolv = 'BiCGSTAB'
           if( solve_sol(ivari) % kfl_algso == SOL_SOLVER_GMRES                 ) solve_sol(ivari) % wsolv = 'GMRES'
           if( solve_sol(ivari) % kfl_algso == SOL_SOLVER_RICHARDSON            ) solve_sol(ivari) % wsolv = 'RICHARDSON'
           if( solve_sol(ivari) % kfl_algso == SOL_SOLVER_MATRIX_DIAGONAL       ) solve_sol(ivari) % wsolv = 'MATRIX_BASED_RICHARDSON'
           if( solve_sol(ivari) % kfl_algso == SOL_SOLVER_DEFLATED_BICGSTAB     ) solve_sol(ivari) % wsolv = 'DEFLATED BiCGSTAB'
           if( solve_sol(ivari) % kfl_algso == SOL_SOLVER_DEFLATED_GMRES        ) solve_sol(ivari) % wsolv = 'DEFLATED GMRES'
           if( solve_sol(ivari) % kfl_algso == SOL_SOLVER_SPARSE_DIRECT         ) solve_sol(ivari) % wsolv = 'SPARSE DIRECT SOLVER'
           if( solve_sol(ivari) % kfl_algso == SOL_SOLVER_STEEPEST_DESCENT      ) solve_sol(ivari) % wsolv = 'STEEPEST DESCENT'
           if( solve_sol(ivari) % kfl_algso == SOL_SOLVER_PIPELINED_CG          ) solve_sol(ivari) % wsolv = 'PIPELINED CONJUGATE GRADIENT'
           if( solve_sol(ivari) % kfl_algso == SOL_SOLVER_PIPELINED_DEFLATED_CG ) solve_sol(ivari) % wsolv = 'PIPELINED DEFLATED CONJUGATE GRADIENT'
           if( solve_sol(ivari) % kfl_algso == SOL_SOLVER_MAPHYS_UNSYMMETRIC    ) solve_sol(ivari) % wsolv = 'MAPHYS with Pack GMRES'
           if( solve_sol(ivari) % kfl_algso == SOL_SOLVER_MAPHYS_SYMMETRIC      ) solve_sol(ivari) % wsolv = 'MAPHYS with Pack CG'
           if( solve_sol(ivari) % kfl_algso == SOL_SOLVER_MAPHYS_SPD            ) solve_sol(ivari) % wsolv = 'MAPHYS with Pack CG for SPD matrix'
           if( solve_sol(ivari) % kfl_algso == SOL_SOLVER_AGMG                  ) solve_sol(ivari) % wsolv = 'AGMG - if krylov==1 Flexible CG is used else GCR'
           if( solve_sol(ivari) % kfl_algso == SOL_SOLVER_PSBLAS                ) solve_sol(ivari) % wsolv = 'PSBLAS'
           if( solve_sol(ivari) % kfl_algso == SOL_SOLVER_PETSC                 ) solve_sol(ivari) % wsolv = 'PETSC'
           if( solve_sol(ivari) % kfl_algso == SOL_SOLVER_MINRES                ) solve_sol(ivari) % wsolv = 'MINRES'
           if( solve_sol(ivari) % kfl_algso == SOL_SOLVER_EXPLICIT              ) solve_sol(ivari) % wsolv = 'EXPLICIT'
           if( solve_sol(ivari) % kfl_algso == SOL_SOLVER_POLYNOMIAL            ) solve_sol(ivari) % wsolv = 'POLYNOMIAL APPROX OF THE INVERSE'
           !
           ! Name of preconditioner
           !
           if( solve_sol(ivari) % kfl_preco == SOL_NO_PRECOND )         solve_sol(ivari) % wprec = 'NONE'
           if( solve_sol(ivari) % kfl_preco == SOL_SQUARE )             solve_sol(ivari) % wprec = 'SYMMETRIC: sqrt[|diag(A)|]'
           if( solve_sol(ivari) % kfl_preco == SOL_DIAGONAL )           solve_sol(ivari) % wprec = 'DIAGONAL diag(A)'
           if( solve_sol(ivari) % kfl_preco == SOL_MATRIX )             solve_sol(ivari) % wprec = 'MATRIX'
           if( solve_sol(ivari) % kfl_preco == SOL_LINELET )            solve_sol(ivari) % wprec = 'LINELET'
           if( solve_sol(ivari) % kfl_preco == SOL_MASS_MATRIX )        solve_sol(ivari) % wprec = 'MASS MATRIX (DIAGONAL)'
           if( solve_sol(ivari) % kfl_preco == SOL_GAUSS_SEIDEL ) then
              if(      solve_sol(ivari) % kfl_renumbered_gs == SYMMETRIC_GAUSS_SEIDEL  ) then
                 solve_sol(ivari) % wprec = 'SYMMETRIC GAUSS-SEIDEL'
              else if( solve_sol(ivari) % kfl_renumbered_gs == STREAMWISE_GAUSS_SEIDEL ) then
                 solve_sol(ivari) % wprec = 'STREAMWISE GAUSS-SEIDEL'
              else if( solve_sol(ivari) % kfl_renumbered_gs == CLASSICAL_GAUSS_SEIDEL  ) then
                 solve_sol(ivari) % wprec = 'GAUSS-SEIDEL'
              end if
           end if
           if( solve_sol(ivari) % kfl_preco == SOL_BIDIAGONAL ) then
               if(     solve_sol(ivari) % kfl_renumbered_gs == STREAMWISE_BIDIAGONAL   ) then
                  solve_sol(ivari) % wprec = 'STREAMWISE BIDIAGONAL'
               else if( solve_sol(ivari) % kfl_renumbered_gs == BIDIAGONAL             ) then
                  solve_sol(ivari) % wprec = 'BIDIAGONAL'
               end if
!              if( solve_sol(ivari) % kfl_renumbered_gs == 1 )           solve_sol(ivari) % wprec = 'STREAMWISE BIDIAGONAL'
!              if( solve_sol(ivari) % kfl_renumbered_gs == 0 )           solve_sol(ivari) % wprec = 'BIDIAGONAL'
           end if
           if( solve_sol(ivari) % kfl_preco == SOL_CLOSE_MASS_MATRIX )  solve_sol(ivari) % wprec = 'CLOSE MASS MATRIX'
           if( solve_sol(ivari) % kfl_preco == SOL_RICHARDSON )         solve_sol(ivari) % wprec = 'RICHARDSON'
           if( solve_sol(ivari) % kfl_preco == SOL_ORTHOMIN )           solve_sol(ivari) % wprec = 'ORTHOMIN(1)'
           if( solve_sol(ivari) % kfl_preco == SOL_APPROXIMATE_SCHUR)   solve_sol(ivari) % wprec = 'APPROXIMATE SCHUR=ABB-ABI diag(AII)^{-1} AIB'
           if( solve_sol(ivari) % kfl_preco == SOL_ABB )                solve_sol(ivari) % wprec = 'ABB'
           if( solve_sol(ivari) % kfl_preco == SOL_MOD_DIAGONAL )       solve_sol(ivari) % wprec = 'diag(ABB) - diag( ABI diag(AII)^{-1} AIB)'
           if( solve_sol(ivari) % kfl_preco == SOL_AII )                solve_sol(ivari) % wprec = 'AII^{-1}'
           if( solve_sol(ivari) % kfl_preco == SOL_MULTIGRID )          solve_sol(ivari) % wprec = 'MULTIGRID DEFLATED'
           if( solve_sol(ivari) % kfl_preco == SOL_LOCAL_DIAGONAL )     solve_sol(ivari) % wprec = 'LOCAL DIAGONAL'
           if( solve_sol(ivari) % kfl_preco == SOL_AUGMENTED_DIAGONAL ) solve_sol(ivari) % wprec = 'AUGMENTED DIAGONAL P: (M/dt+P) u = M/dt u+b'
           if( solve_sol(ivari) % kfl_preco == SOL_RAS ) then
              if( solve_sol(ivari) %kfl_block_ras == 0 ) then
                 if( solve_sol(ivari) % kfl_add_schwarz == 1 ) then
                    solve_sol(ivari) % wprec = 'ADDITIVE SCHWARZ (AS)'
                 else
                    solve_sol(ivari) % wprec = 'RESTRICTED ADDITIVE SCHWARZ (RAS)'
                 end if
              else
                 if( solve_sol(ivari) % kfl_add_schwarz == 1 ) then
                    solve_sol(ivari) % wprec = 'BLOCKWISE ADDITIVE SCHWARZ (AS)'
                 else
                    solve_sol(ivari) % wprec = 'BLOCKWISE RESTRICTED ADDITIVE SCHWARZ (RAS)'
                 end if
              end if
           end if
           !
           ! Complex solver
           !
           if( solve_sol(ivari) % kfl_cmplx == 1 ) &
                solve_sol(ivari) % wsolv = 'COMPLEX '//trim(solve_sol(ivari) % wsolv)

        end if
        
     end do
     
     !----------------------------------------------------------
     !
     ! Block matrix: define sizes
     !
     !----------------------------------------------------------
     
     do ivari = 1,size(solve_sol,KIND=ip)
        if( solve_sol(ivari) % num_blocks == 2 ) then
           i1                               = ivari
           i2                               = i1+1
           n1                               = solve_sol(i1) % nzmat / solve_sol(i1) % ndofn * solve_sol(i2) % ndofn
           n2                               = solve_sol(i2) % nzmat / solve_sol(i2) % ndofn * solve_sol(i1) % ndofn
           solve_sol(i1) % block_sizes(1,1) = solve_sol(i1) % nzmat
           solve_sol(i1) % block_sizes(1,2) = max(n1,n2)
           solve_sol(i1) % block_sizes(2,1) = max(n1,n2)
           solve_sol(i1) % block_sizes(2,2) = solve_sol(2) % nzmat
           nn                               = solve_sol(i1) % block_sizes(1,1)+solve_sol(i1) % block_sizes(1,2) + &
                &                             solve_sol(i1) % block_sizes(2,1)+solve_sol(i1) % block_sizes(2,2)
           !poa11_nsi                       =  1
           !poa12_nsi                       =  poa11_nsi + solve_sol(i1) % block_sizes(1,1)
           !poa21_nsi                       =  poa12_nsi + solve_sol(i1) % block_sizes(1,2)
           !poa22_nsi                       =  poa21_nsi + solve_sol(i1) % block_sizes(2,1)
           nzmat                            =  max(nzmat,nn)
           nzrhs                            =  max(nzrhs,solve_sol(i1) % nunkn+solve_sol(i2) % nunkn)
        end if
     end do
     
  end if

end subroutine soldef

subroutine solcpy(ipro1,ipro2)
  use def_master
  use def_solver
  implicit none
  integer(ip),  intent(in) :: ipro1,ipro2
  type(soltyp)             :: solve_tmp
  !
  ! Copy solver
  !
  solve_tmp        = solve_sol(ipro2)
  solve_sol(ipro2) = solve_sol(ipro1)
  !
  ! Recover selected orginal values
  !
  solve_sol(ipro2) % wprob = solve_tmp%wprob

end subroutine solcpy
