!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine inivar(itask)
  !-----------------------------------------------------------------------
  !****f* master/inivar
  ! NAME
  !    inivar
  ! DESCRIPTION
  !    This routine initialize some variables
  ! USES
  ! USED BY
  !    Reapro
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_elmtyp
  use def_parame
  use def_postpr
  use def_solver
  use mod_alya2maphys,     only : alya2maphys_initialization
  use mod_memory,          only : memory_alloca_min
  use mod_memory,          only : memory_alloca
  use mod_memory,          only : memory_copy
  use mod_parall,          only : PAR_COMM_MY_CODE_ARRAY
  use mod_parall,          only : par_omp_num_threads
  use def_coupli,          only : RESIDUAL
  use def_coupli,          only : BETWEEN_ZONES
  use mod_couplings,       only : COU_LIST_SOURCE_NODES
  use mod_communications,  only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_graphs,          only : graphs_permut_metis_postordering
  use mod_graphs,          only : graphs_copyij
  use mod_direct_solver,   only : direct_solver_initialization
  use mod_alya2agmg,       only : alya2agmg_initialization
  use mod_messages,        only : messages_live
  use mod_std
#ifdef INC_PSBLAS
  use mod_alya2psblas,     only : alya2psblas_initialization
#endif
  use mod_linelet,         only : linelet_initialization
  use mod_direct_solver,   only : direct_solver_type_initialization
  use mod_couplings,       only : couplings_initialize_solver
  use mod_solver,          only : solver_explicit_initialization
  use mod_solver,          only : solver_ras_initialization
  use mod_groups,          only : groups_check_prescription
  use mod_moduls_conf,     only : moduls_allocate_timing

  implicit none
  integer(ip), intent(in)      :: itask
  integer(ip)                  :: ipoin,idofn,iflin,ivari
  integer(ip)                  :: ndofn,iflim,pelty,num_blocks
  integer(ip)                  :: jzdom,ndofn_block
  integer(ip)                  :: jblok,ndofn_kblok,ndofn_jblok,kblok

  select case( itask )

  case ( 1_ip )
     !
     ! Postprocess
     !
     varnu_pos=0
     varna_pos='NULL'
     !
     ! Timings
     !
     call moduls_allocate_timing(11_ip,CURRENT_MODULE=0_ip)
     
     call times( 1) % register('amr'                   ,'kernel','cc')
     call times( 2) % register('adaptation'            ,'amr'   ,'computation')
     call times( 3) % register('interface motion'      ,'amr'   ,'cc')
     call times( 4) % register('repartition mesh'      ,'amr'   ,'cc')
     call times( 5) % register('parallel interface'    ,'amr'   ,'cc')
     call times( 6) % register('global numbering'      ,'amr'   ,'cc')
     call times( 7) % register('mesh interpolation'    ,'amr'   ,'cc')
     call times( 8) % register('size interpolation'    ,'amr'   ,'computation')

     call times( 9) % register('coupling'             ,'kernel')
     call times(10) % register('initialization'       ,'coupling')
     call times(11) % register('driver'               ,'coupling')

  case ( 2_ip )
     !
     ! Cancel Hessian calculation
     !
     if( kfl_lapla == 0 ) then
        do pelty = 1,nelty
           llapl(pelty) = 0
        end do
     end if
     !
     ! Initial time step
     !
     ittim_ini = ittim

     !----------------------------------------------------------
     !
     ! Coupling for solvers
     !
     !----------------------------------------------------------

     call messages_live('SOLVER INITIALIZATION','START SECTION')
     call couplings_initialize_solver()
     !
     ! Prepare module solvers
     !
     do modul = 1,mmodu
        if( kfl_modul(modul) == 1 ) then
           if( associated(momod(modul) % solve) ) then
              current_zone = lzone(modul)
              do ivari = 1,size(momod(modul) % solve,KIND=ip)
                 solve_sol => momod(modul) % solve(ivari:)

                 if( solve_sol(1) % kfl_algso /= SOL_NO_SOLVER ) then

                    call messages_live(namod(modul)//': '//trim(solve_sol(1) % wsolv)//' FOR '//trim(solve_sol(1) % wprob))
                    !
                    ! Krylov subspace saving
                    !
                    if( solve_sol(1) % kfl_save_krylov /= 0 ) then
                       if( solve_sol(1) % kfl_save_krylov < 0 ) then
                          if( solve_sol(1) % kfl_algso == SOL_SOLVER_GMRES ) then
                             solve_sol(1) % kfl_save_krylov = solve_sol(1) % nkryd
                          else
                             solve_sol(1) % kfl_save_krylov = solve_sol(1) % miter
                          end if
                       end if
                       if( INOTMASTER ) then
                          call memory_alloca(mem_modul(1:2,modul),'SOLVE % KRYLOV_DIR','inivar',&
                               solve_sol(1) % krylov_dir,solve_sol(1) % nunkn,solve_sol(1) % kfl_save_krylov)
                          call memory_alloca(mem_modul(1:2,modul),'SOLVE % KRYLOV_DOT_PRODUCTS','inivar',&
                               solve_sol(1) % krylov_dot_products,solve_sol(1) % kfl_save_krylov)
                       end if
                    end if
                    !
                    ! ELL format
                    !
                    solve_sol(1) % nnz_ell =   meshe(ndivi) % nzdom_ell
                    solve_sol(1) % cols_ell => meshe(ndivi) % ell_cols
                    !
                    ! Schur based solver
                    !
                    if( solve_sol(1) % kfl_schur == 1 ) then
                       solve_sol(1) % poaii     =  1
                       solve_sol(1) % poaib     =  solve_sol(1) % poaii + nzdom_aii * solve_sol(1) % ndof2
                       solve_sol(1) % poabi     =  solve_sol(1) % poaib + nzdom_aib * solve_sol(1) % ndof2
                       solve_sol(1) % poabb     =  solve_sol(1) % poabi + nzdom_abi * solve_sol(1) % ndof2
                       solve_sol(1) % r_dom_aii => r_dom_aii
                       solve_sol(1) % c_dom_aii => c_dom_aii 
                       solve_sol(1) % r_dom_aib => r_dom_aib
                       solve_sol(1) % c_dom_aib => c_dom_aib 
                       solve_sol(1) % r_dom_abi => r_dom_abi 
                       solve_sol(1) % c_dom_abi => c_dom_abi 
                       solve_sol(1) % r_dom_abb => r_dom_abb 
                       solve_sol(1) % c_dom_abb => c_dom_abb 
                    end if
                    !
                    ! Point to corresponding graph
                    !
                    if(      solve_sol(1) % kfl_where == SOL_NODES ) then

                       if( solve_sol(1) % kfl_full_rows == 0 ) then
                          solve_sol(1) % ia   => meshe(ndivi) % r_dom
                          solve_sol(1) % ja   => meshe(ndivi) % c_dom
                          solve_sol(1) % rows => meshe(ndivi) % coo_rows
                          solve_sol(1) % cols => meshe(ndivi) % coo_cols
                          solve_sol(1) % nnz  =  meshe(ndivi) % nzdom
                       else
                          solve_sol(1) % ia   => meshe(ndivi) % r_dom
                          solve_sol(1) % ja   => meshe(ndivi) % c_dom
                          solve_sol(1) % rows => meshe(ndivi) % coo_rows
                          solve_sol(1) % cols => meshe(ndivi) % coo_cols
                          solve_sol(1) % nnz  =  meshe(ndivi) % nzdom_own
                       end if
                       solve_sol(1) % ia_full     => meshe(ndivi) % r_dom_own
                       solve_sol(1) % ja_full     => meshe(ndivi) % c_dom_own
                       solve_sol(1) % ia_full_end => meshe(ndivi) % r_dom_end
                       solve_sol(1) % ia_full_ini => meshe(ndivi) % r_dom_ini

                    else if( solve_sol(1) % kfl_where == SOL_ELEMENTS ) then
                       solve_sol(1) % ia   => meshe(ndivi) % r_elm_2
                       solve_sol(1) % ja   => meshe(ndivi) % c_elm_2
                       solve_sol(1) % nnz  =  meshe(ndivi) % nzelm_2
                    else if( solve_sol(1) % kfl_where == SOL_EDGES    ) then
                       solve_sol(1) % ia   => meshe(ndivi) % r_edg
                       solve_sol(1) % ja   => meshe(ndivi) % c_edg
                       solve_sol(1) % nnz  =  meshe(ndivi) % nzedg
                    end if

                    if(  solve_sol(1) % kfl_algso == SOL_SOLVER_EXPLICIT ) then
                       !
                       ! Explicit solvers
                       !
                       call solver_explicit_initialization(solve_sol(1))

                    else if(  solve_sol(1) % kfl_algso == SOL_SOLVER_DIRECT        .or. &
                         solve_sol(1) % kfl_algso == SOL_SOLVER_PASTIX        .or. &
                         solve_sol(1) % kfl_algso == SOL_SOLVER_MUMPS         .or. &
                         solve_sol(1) % kfl_algso == SOL_SOLVER_WSMP          .or. &
                         solve_sol(1) % kfl_algso == SOL_SOLVER_PWSMP         .or. &
                         solve_sol(1) % kfl_algso == SOL_SOLVER_GPUQR         .or. &
                         solve_sol(1) % kfl_algso == SOL_SOLVER_SPARSE_DIRECT ) then
                       !
                       ! Direct solver
                       !
                       call direct_solver_type_initialization(solve_sol(1) % direct_solver)
                       if( solve_sol(1) % kfl_algso == SOL_SOLVER_PWSMP ) then
                          solve_sol(1) % direct_solver % ia          => solve_sol(1) % ia_full
                          solve_sol(1) % direct_solver % ja          => solve_sol(1) % ja_full
                          solve_sol(1) % direct_solver % ndof        =  solve_sol(1) % ndofn
                          solve_sol(1) % direct_solver % nn          =  npoin_own
                          if( INOTMASTER ) solve_sol(1) % direct_solver % nz =  solve_sol(1) % ia_full(npoin_own+1)-1
                       else
                          solve_sol(1) % direct_solver % ia          => solve_sol(1) % ia
                          solve_sol(1) % direct_solver % ja          => solve_sol(1) % ja
                          solve_sol(1) % direct_solver % ndof        =  solve_sol(1) % ndofn
                          solve_sol(1) % direct_solver % nn          =  solve_sol(1) % nequa
                          solve_sol(1) % direct_solver % nz          =  solve_sol(1) % nnz
                       end if
                       solve_sol(1) % direct_solver % kfl_paral   =  1_ip
                       solve_sol(1) % direct_solver % num_threads =  max(1_ip,par_omp_num_threads)
                       solve_sol(1) % direct_solver % nn_own      =  solve_sol(1) % nequa_own
                       solve_sol(1) % direct_solver % nz_own      =  solve_sol(1) % nzmat_own / (solve_sol(1) % ndofn ** 2)

                       if(       solve_sol(1) % kfl_algso == SOL_SOLVER_DIRECT ) then

                          solve_sol(1) % direct_solver % kfl_solver =  SOL_DIRECT_SOLVER_ALYA

                       else if( solve_sol(1) % kfl_algso == SOL_SOLVER_PASTIX ) then

                          solve_sol(1) % direct_solver % kfl_solver =  SOL_DIRECT_SOLVER_PASTIX

                       else if( solve_sol(1) % kfl_algso == SOL_SOLVER_MUMPS  ) then

                          solve_sol(1) % direct_solver % kfl_solver =  SOL_DIRECT_SOLVER_MUMPS

                       else if(solve_sol(1) % kfl_algso == SOL_SOLVER_WSMP ) then

                          solve_sol(1) % direct_solver % kfl_solver =  SOL_DIRECT_SOLVER_WSMP

                       else if(solve_sol(1) % kfl_algso == SOL_SOLVER_PWSMP ) then

                          solve_sol(1) % direct_solver % kfl_solver = SOL_DIRECT_SOLVER_PWSMP

                       else if(solve_sol(1) % kfl_algso == SOL_SOLVER_GPUQR ) then

                          solve_sol(1) % direct_solver % kfl_solver =  SOL_DIRECT_SOLVER_GPUQR

                       else

                          solve_sol(1) % direct_solver % kfl_solver =  kfl_direct_solver

                       end if

                       call direct_solver_initialization(solve_sol(1) % direct_solver)
                       solve_sol(1) % kfl_algso = SOL_SOLVER_SPARSE_DIRECT

                    else if( solve_sol(1) % kfl_algso == SOL_SOLVER_MAPHYS_UNSYMMETRIC .or. &
                         &   solve_sol(1) % kfl_algso == SOL_SOLVER_MAPHYS_SYMMETRIC   .or. &
                         &   solve_sol(1) % kfl_algso == SOL_SOLVER_MAPHYS_SPD         ) then
                       !
                       ! MAPHYS solver
                       !
                       if( ISEQUEN ) then
                          call runend('MAPHYS CAN ONLY RUN IN PARALLEL')
                       else
                          call alya2maphys_initialization(meshe(ndivi),PAR_COMM_MY_CODE_ARRAY(1),solve_sol(1))
                       end if

                    else if( solve_sol(1) % kfl_algso == SOL_SOLVER_AGMG ) then
                       !
                       ! AGMG solver
                       !
                       if( kfl_paral == -1 ) then
                          !call runend('AGMG CAN ONLY RUN IN PARALLEL')  !Now the serial version should also work
                       else
#ifdef PARAL_AGMG
                          call alya2agmg_initialization(meshe(ndivi),solve_sol(1))
#endif
                       end if

                    else if( solve_sol(1) % kfl_algso == SOL_SOLVER_PSBLAS ) then
                       !
                       ! PSBLAS solver
                       !
                       if( kfl_paral == -1 ) then
                          !call runend('PSBLAS CAN ONLY RUN IN PARALLEL')
                       else
#ifdef INC_PSBLAS
                          call alya2psblas_initialization(meshe(ndivi),solve_sol(1))
#endif
                       endif

                    else if( solve_sol(1) % kfl_algso == SOL_SOLVER_PETSC ) then
                       !
                       ! PETSc solver library
                       !
#ifdef PETSC
                     block
                     use mod_alya2petsc, only: alya2petsc_initialiseLinearSolver
                     use mod_alya2petsc, only: alya2petsc_createLinearSolver
                     call alya2petsc_initialiseLinearSolver(solve_sol(1) % petsc)
                     call alya2petsc_createLinearSolver(solve_sol(1) % petsc, solve_sol(1))
                     endblock
#endif
                    else

                       iflin     =  0
                       iflim     =  0
                       !
                       ! linelet
                       !
                       if(    solve_sol(1) % kfl_preco == SOL_LINELET   .or.  &
                            ( solve_sol(1) % kfl_preco == SOL_MULTIGRID .and. &
                            & solve_sol(1) % kfl_defpr == SOL_LINELET   ) ) then
                          iflin = 1
                          if( INOTMASTER ) then
                             ndofn = size(solve_sol(1) % kfl_fixno,1,KIND=ip)
                             !if( ndofn > 1 ) then
                             call memory_alloca(mem_modul(1:2,modul),'SOLVE % LIMLI','inivar',solve_sol(1) % limli,npoin)
                             do ipoin = 1,npoin
                                do idofn = 1,ndofn
                                   if( momod(modul) % solve(ivari) % kfl_fixno(idofn,ipoin) > 0 ) solve_sol(1) % limli(ipoin) = 1
                                end do
                             end do
                             !else
                             !   solve_sol(1) % limli => momod(modul) % solve(ivari) % kfl_fixno(1,1:)
                             !end if
                          end if
                          call linelet_initialization()
                       end if
                       !
                       ! Possible errors
                       !
                       if( solve_sol(1) % kfl_preco == SOL_RAS .and. kfl_graph == 0 ) then
                          call runend('INIVAR: TO USE RAS, ADD EXTENDED_GRAPH=ON TO THE MESH SECTION OF THE ker.dat FILE')
                       end if
                       !if( solve_sol(1) % kfl_preco == SOL_RAS .and. ISEQUEN ) then
                       !   call runend('INIVAR: RAS PRECONDITIONER DOES NOT MAKE SENSE IN SEQUENTIAL!')
                       !end if
                       !
                       ! Streamwise Gauss-Seidel preconditioner: point to advection field
                       !
                       if(  solve_sol(1) % kfl_preco         == SOL_GAUSS_SEIDEL .and. &
                            solve_sol(1) % kfl_renumbered_gs == STREAMWISE_GAUSS_SEIDEL ) then
                          solve_sol(1) % vecto_gs => advec(:,:,1)
                       end if
                       !
                       ! Streamwise Bidiagonal preconditioner: point to advection field
                       !
                       if( solve_sol(1) % kfl_preco         == SOL_BIDIAGONAL .and. &
                            solve_sol(1) % kfl_renumbered_gs == STREAMWISE_BIDIAGONAL ) then
                          solve_sol(1) % vecto_gs => advec(:,:,1)
                       end if
                       !
                       ! Deflated CG / multigrid
                       !
                       if(    solve_sol(1) % kfl_algso  == SOL_DEFLATED_CG                  .or. &
                            & solve_sol(1) % kfl_algso  == SOL_GDECG                        .or. &
                            & solve_sol(1) % kfl_preco  == SOL_MULTIGRID                    .or. &
                            & solve_sol(1) % kfl_algso  == SOL_SOLVER_PIPELINED_DEFLATED_CG .or. &
                            & solve_sol(1) % kfl_algso  == SOL_SOLVER_A_DEF2                .or. &
                            & solve_sol(1) % kfl_coarse == 1    ) then
                          if( INOTMASTER ) then

                             if( associated(solve_sol(1) % kfl_fixno) ) then
                                ndofn = size(solve_sol(1) % kfl_fixno,1,KIND=ip)
                                if( iflin == 1 .and. ndofn > 1 ) then
                                   call memory_copy(mem_modul(1:2,modul),'SOLVE % LIMPO','inivar',momod(modul) % solve(ivari) % limli,solve_sol(1) % limpo)
                                else
                                   call memory_alloca(mem_modul(1:2,modul),'SOLVE % LIMPO','inivar',solve_sol(1) % limpo,npoin)
                                   do ipoin = 1,npoin
                                      solve_sol(1) % limpo(ipoin) = 0
                                      do idofn = 1,ndofn
                                         if( solve_sol(1) % kfl_fixno(idofn,ipoin) > 0 ) then
                                            solve_sol(1) % limpo(ipoin) = 1
                                         end if
                                      end do
                                   end do
                                end if
                             else
                                nullify(solve_sol(1) % limpo)
                             end if
                          end if
                          call groups_check_prescription(solve_sol(1) % limpo)
                          if( kfl_ngrou == 0 ) then
                             call runend('INIVAR: NUMBER OF GROUPS AND METHOD SHOULD BE DECLARED IN STRATEGY FIELD')
                          else
                             call cregro()
                          end if
                       end if
                       !
                       ! Multigrid preconditioner
                       !
                       if( solve_sol(1) % kfl_preco  == SOL_MULTIGRID ) then
                          call direct_solver_type_initialization(solve_sol(1) % direct_solver_AMG)
                          solve_sol(1) % direct_solver_AMG % ia                  => solve_sol(1) % iagro
                          solve_sol(1) % direct_solver_AMG % ja                  => solve_sol(1) % jagro
                          solve_sol(1) % direct_solver_AMG % ndof                =  solve_sol(1) % ndofn
                          solve_sol(1) % direct_solver_AMG % kfl_solver          =  kfl_direct_solver
                          solve_sol(1) % direct_solver_AMG % nn                  =  solve_sol(1) % ngrou
                          solve_sol(1) % direct_solver_AMG % nz                  =  solve_sol(1) % nzgro
                          solve_sol(1) % direct_solver_AMG % num_threads         =  1
                          call direct_solver_initialization(solve_sol(1) % direct_solver_AMG)
                       end if
                       !
                       ! Deflation
                       !
                       if(  solve_sol(1) % kfl_algso == SOL_GDECG         .or. &
                            solve_sol(1) % kfl_algso == SOL_SOLVER_A_DEF2 .or. &
                            solve_sol(1) % kfl_algso == SOL_DEFLATED_CG   .or. &
                            solve_sol(1) % kfl_algso == SOL_SOLVER_PIPELINED_DEFLATED_CG ) then
                          call direct_solver_type_initialization(solve_sol(1) % direct_solver_Deflation)
                          solve_sol(1) % direct_solver_Deflation % ia            => solve_sol(1) % iagro
                          solve_sol(1) % direct_solver_Deflation % ja            => solve_sol(1) % jagro
                          solve_sol(1) % direct_solver_Deflation % idiag         => solve_sol(1) % idiag
                          solve_sol(1) % direct_solver_Deflation % iskyl         => solve_sol(1) % iskyl
                          solve_sol(1) % direct_solver_Deflation % ndof          =  solve_sol(1) % ndofn
                          if( solve_sol(1) % kfl_defas == 0 ) then
                             solve_sol(1) % direct_solver_Deflation % kfl_solver =  SOL_DIRECT_SOLVER_SKYLINE_ALYA
                          else
                             solve_sol(1) % direct_solver_Deflation % kfl_solver =  kfl_direct_solver
                          end if
                          solve_sol(1) % direct_solver_Deflation % nn            =  solve_sol(1) % ngrou
                          solve_sol(1) % direct_solver_Deflation % nn_own        =  solve_sol(1) % ngrou
                          solve_sol(1) % direct_solver_Deflation % nz            =  solve_sol(1) % nzgro
                          solve_sol(1) % direct_solver_Deflation % nz_own        =  solve_sol(1) % nzgro
                          solve_sol(1) % direct_solver_Deflation % nskyl         =  solve_sol(1) % nskyl
                          solve_sol(1) % direct_solver_Deflation % num_threads   =  max(1_ip,par_omp_num_threads)
                          call direct_solver_initialization(solve_sol(1) % direct_solver_Deflation)
                       end if
                       !
                       ! AII preconditioner
                       !
                       if( solve_sol(1) % kfl_preco == SOL_AII ) then
                          call direct_solver_type_initialization(solve_sol(1) % direct_solver_block_LU)
                          solve_sol(1) % direct_solver_block_LU % ia          => r_dom_aii
                          solve_sol(1) % direct_solver_block_LU % ja          => c_dom_aii
                          solve_sol(1) % direct_solver_block_LU % ndof        =  solve_sol(1) % ndofn
                          solve_sol(1) % direct_solver_block_LU % kfl_solver  =  kfl_direct_solver
                          solve_sol(1) % direct_solver_block_LU % nn          =  npoi1
                          solve_sol(1) % direct_solver_block_LU % nz          =  nzdom_aii
                          solve_sol(1) % direct_solver_block_LU % num_threads =  max(1_ip,par_omp_num_threads)
                          call direct_solver_initialization(solve_sol(1) % direct_solver_block_LU,'NOT MASTER')
                       end if
                       !
                       ! RAS preconditioner
                       !
                       if( solve_sol(1) % kfl_preco  == SOL_RAS ) then
                          call solver_ras_initialization(solve_sol(1),kfl_direct_solver)
                       end if
                       !
                       ! Coarse solver
                       !
                       if( solve_sol(1) % kfl_coarse == 1 ) then
                          call direct_solver_type_initialization(solve_sol(1) % direct_solver_coarse)
                          solve_sol(1) % direct_solver_coarse % ia          => solve_sol(1) % iagro
                          solve_sol(1) % direct_solver_coarse % ja          => solve_sol(1) % jagro
                          solve_sol(1) % direct_solver_coarse % ndof        =  solve_sol(1) % ndofn
                          solve_sol(1) % direct_solver_coarse % kfl_solver  =  kfl_direct_solver
                          solve_sol(1) % direct_solver_coarse % nn          =  solve_sol(1) % ngrou
                          solve_sol(1) % direct_solver_coarse % nz          =  solve_sol(1) % nzgro
                          solve_sol(1) % direct_solver_coarse % num_threads =  1
                          call direct_solver_initialization(solve_sol(1) % direct_solver_coarse)
                       end if
                       !
                       ! Symmetry of direct solvers
                       !
                       if( solve_sol(1) % kfl_symeq == 1 ) then
                          solve_sol(1) % direct_solver_AMG       % kfl_symmetric = 1
                          solve_sol(1) % direct_solver_Deflation % kfl_symmetric = 1
                          solve_sol(1) % direct_solver_block_LU  % kfl_symmetric = 1
                          solve_sol(1) % direct_solver_coarse    % kfl_symmetric = 1
                          if( associated(solve_sol(1) % direct_solver_RAS) ) solve_sol(1) % direct_solver_RAS % kfl_symmetric = 1
                          solve_sol(1) % direct_solver           % kfl_symmetric = 1
                       end if

                    end if

                 end if

              end do

              !----------------------------------------------------------
              !
              ! Determine if and where reaction forces should be computed
              !
              !----------------------------------------------------------

#ifndef COMMDOM
              
              do ivari = 1,size(momod(modul) % solve,KIND=ip)
                 solve_sol => momod(modul) % solve(ivari:)

                 if( solve_sol(1) % kfl_algso /= SOL_NO_SOLVER ) then

                    ndofn      = solve_sol(1) % ndofn
                    num_blocks = solve_sol(1) % num_blocks

                    if( solve_sol(1) % block_num == 1 .and. INOTMASTER ) then
                       !
                       ! Allocate reaction arrays
                       !
                       if( solve_sol(1) % kfl_react > 0 ) then
                          call memory_alloca(mem_modul(1:2,modul),'SOLVE % REACTION','inivar',solve_sol(1) % reaction,ndofn,npoin)
                          if( .not. associated(solve_sol(1) % lpoin_reaction) )&
                               call memory_alloca(mem_modul(1:2,modul),'SOLVE % LPOIN_REACTION','inivar',solve_sol(1) % lpoin_reaction,npoin)
                          do kblok = 2,num_blocks
                             ndofn_block = solve_sol(1) % block_dimensions(kblok)
                             call memory_alloca(mem_modul(1:2,modul),'SOLVE % REACTION','inivar',solve_sol(kblok) % reaction,ndofn_block,npoin)
                          end do
                       end if
                       !
                       ! Compute reaction on all Dirichlet nodes
                       !
                       if( solve_sol(1) % kfl_react == 2 .and. associated(solve_sol(1) % kfl_fixno) ) then
                          do ipoin = 1,npoin
                             if( maxval(solve_sol(1) % kfl_fixno(:,ipoin)) > 0 ) solve_sol(1) % lpoin_reaction(ipoin) = .true.
                          end do
                       end if
                       !
                       ! Compute reaction on all nodes
                       !
                       if( solve_sol(1) % kfl_react == 3 ) then
                          do ipoin = 1,npoin
                             solve_sol(1) % lpoin_reaction(ipoin) = .true.
                          end do
                       end if

                       if( associated(solve_sol(1) % lpoin_reaction) ) then
                          call PAR_INTERFACE_NODE_EXCHANGE(solve_sol(1) % lpoin_reaction,'OR','IN CURRENT ZONE')

                          if( num_blocks == 1 .and. solve_sol(1) % kfl_react > 0 ) then
                             !
                             ! Monolithic system
                             !
                             allocate( solve_sol(1) % lpoin_block(npoin) )
                             do ipoin = 1,npoin
                                nullify(solve_sol(1) % lpoin_block(ipoin) % block1_num)
                                nullify(solve_sol(1) % lpoin_block(ipoin) % block2_num)
                             end do
                             ndofn = solve_sol(1) % ndofn
                             do ipoin = 1,npoin
                                if( solve_sol(1) % lpoin_reaction(ipoin) ) then
                                   jzdom = r_dom(ipoin+1) - r_dom(ipoin)

                                   allocate( solve_sol(1) % lpoin_block(ipoin) % block2_num(1,1)                             )
                                   allocate( solve_sol(1) % lpoin_block(ipoin) % block1_num(1)                               )

                                   nullify(  solve_sol(1) % lpoin_block(ipoin) % block1_num(1)   % reaction                  )
                                   nullify(  solve_sol(1) % lpoin_block(ipoin) % block1_num(1)   % rhs                       )
                                   nullify(  solve_sol(1) % lpoin_block(ipoin) % block1_num(1)   % bvess                     )
                                   nullify(  solve_sol(1) % lpoin_block(ipoin) % block1_num(1)   % kfl_fixno                 )
                                   nullify(  solve_sol(1) % lpoin_block(ipoin) % block1_num(1)   % bvnat                     )
                                   nullify(  solve_sol(1) % lpoin_block(ipoin) % block2_num(1,1) % matrix                    )

                                   allocate( solve_sol(1) % lpoin_block(ipoin) % block1_num(1)   % rhs(ndofn)                )
                                   allocate( solve_sol(1) % lpoin_block(ipoin) % block2_num(1,1) % matrix(ndofn,ndofn,jzdom) )

                                   solve_sol(1) % lpoin_block(ipoin) % block1_num(1)   % rhs    = 0.0_rp
                                   solve_sol(1) % lpoin_block(ipoin) % block2_num(1,1) % matrix = 0.0_rp
                                end if
                             end do

                          else if( solve_sol(1) % kfl_react > 0 ) then
                             !
                             ! ( n x n ) block system
                             !
                             allocate( solve_sol(1) % lpoin_block(npoin) )
                             do ipoin = 1,npoin
                                nullify(solve_sol(1) % lpoin_block(ipoin) % block1_num)
                                nullify(solve_sol(1) % lpoin_block(ipoin) % block2_num)
                             end do
                             do ipoin = 1,npoin
                                if( solve_sol(1) % lpoin_reaction(ipoin) ) then
                                   allocate( solve_sol(1) % lpoin_block(ipoin) % block2_num(num_blocks,num_blocks) )
                                   allocate( solve_sol(1) % lpoin_block(ipoin) % block1_num(num_blocks) )
                                   do iblok = 1,num_blocks
                                      nullify(  solve_sol(1) % lpoin_block(ipoin) % block1_num(iblok)   % reaction            )
                                      nullify(  solve_sol(1) % lpoin_block(ipoin) % block1_num(iblok)   % rhs                 )
                                      nullify(  solve_sol(1) % lpoin_block(ipoin) % block1_num(iblok)   % bvess               )
                                      nullify(  solve_sol(1) % lpoin_block(ipoin) % block1_num(iblok)   % kfl_fixno           )
                                      nullify(  solve_sol(1) % lpoin_block(ipoin) % block1_num(iblok)   % bvnat               )
                                      do jblok = 1,num_blocks
                                         nullify(  solve_sol(1) % lpoin_block(ipoin) % block2_num(jblok,iblok)   % matrix     )
                                      end do
                                   end do
                                   jzdom = r_dom(ipoin+1) - r_dom(ipoin)
                                   do kblok = 1,num_blocks
                                      ndofn_kblok = solve_sol(1) % block_dimensions(kblok)
                                      allocate( solve_sol(1) % lpoin_block(ipoin) % block1_num(kblok) % rhs(ndofn_kblok) )
                                      solve_sol(1) % lpoin_block(ipoin) % block1_num(kblok) % rhs = 0.0_rp
                                      do jblok = 1,num_blocks
                                         ndofn_jblok = solve_sol(1) % block_dimensions(jblok)
                                         allocate( solve_sol(1) % lpoin_block(ipoin) % block2_num(kblok,jblok) % matrix(ndofn_jblok,ndofn_kblok,jzdom) )
                                         solve_sol(1) % lpoin_block(ipoin) % block2_num(kblok,jblok) % matrix = 0.0_rp
                                      end do
                                   end do
                                end if
                             end do
                          end if
                       end if
                    end if
                 end if

              end do

              !----------------------------------------------------------------------------
              !
              ! Allocate minimum memory to some arrays
              !
              !----------------------------------------------------------------------------

              do ivari = 1,size(momod(modul) % solve,KIND=ip)
                 solve_sol => momod(modul) % solve(ivari:)
                 if( solve_sol(1) % kfl_algso /= SOL_NO_SOLVER ) then
                    call memory_alloca_min(mem_modul(1:2,modul),'SOLVE % REACTION','inivar',solve_sol(1) % reaction)
                    call memory_alloca_min(mem_modul(1:2,modul),'SOLVE % BVNAT'   ,'inivar',solve_sol(1) % bvnat)
                    if( solve_sol(1) % block_num == 1 ) then
                       do iblok = 1,solve_sol(1) % num_blocks
                          call memory_alloca_min(mem_modul(1:2,modul),'SOLVE % REACTION','inivar',solve_sol(iblok) % reaction)
                          call memory_alloca_min(mem_modul(1:2,modul),'SOLVE % BLOCK_ARRAY % BVNAT','inivar',solve_sol(1) % block_array(iblok) % bvnat)
                       end do
                    end if
                 end if
              end do
#endif
           end if
        end if
     end do

     call messages_live('SOLVER INITIALIZATION','END SECTION')

     modul = 0
     current_zone = 0

  end select

end subroutine inivar
