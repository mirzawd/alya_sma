!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup ???
!> @{
!> @file    reasol.f90
!> @author  houzeaux
!> @date    2020-12-02
!> @brief   This routine reads the solver data
!> @details All solver options
!> @} 
!-----------------------------------------------------------------------

subroutine reasol(itask)
 !.md<module>kernel
 !.md<input>case.sol.dat
 !.md<pos>0
 !.md<sec>

  use def_kintyp, only    :  ip
  use def_solver, only    :  SOL_SOLVER_DEFLATED_CG,SOL_GDECG
  use def_solver, only    :  SOL_SOLVER_PIPELINED_DEFLATED_CG
  use def_solver, only    :  SOL_SOLVER_A_DEF2
  use mod_ecoute, only    :  ecoute
  use def_inpout
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: kfl_check

  kfl_check = 0

  if( itask == 1 ) then

     if( trim(words(2)) /= '' .and. words(2) /= 'VARIA' ) then

        call intsol(1_ip)
        call intsol(3_ip)
        call intsol(4_ip)
        call intsol(5_ip)

     else

        kfl_check = 1
        do while( words(1) /= 'ENDAL' )
           if( words(1) == 'SOLVE' ) then
              call intsol(1_ip)
           else if( words(1) == 'PRECO' ) then
              call intsol(2_ip)
           else if( words(1) == 'OUTPU' ) then
              call intsol(3_ip)
           else if( words(1) == 'CONVE' ) then
              call intsol(4_ip)
           else if( words(1) == 'OPTIO' ) then
              call intsol(5_ip)
           else if( words(1) == 'SMOOT' ) then
              call intsol(6_ip)
           else if( words(1) == 'COARS' ) then
              call intsol(7_ip)
           else if( words(1) == 'RESID' ) then
              call intsol(8_ip)
           else if( words(1) == 'CONFI' ) then
              call intsol(9_ip)
           else if( words(1) == 'RHSMI' ) then
              call intsol(10_ip)
           else if( words(1) == 'MATRI' ) then
              call intsol(11_ip)
           else if( words(1) == 'OPENM' ) then
              call intsol(12_ip)
           else if( words(1) == 'BLOCK' ) then
              call intsol(13_ip)
           else if( words(1) == 'RELAX' ) then
              call intsol(14_ip)
           end if
           call ecoute('REASOL')
        end do

     end if

  end if


  if( kfl_check == 0 .and. ( itask == 1 .or. itask == 2 ) ) then

     if( trim(words(1)) /= '' ) then
        call intsol(2_ip)
     end if
  end if

end subroutine reasol

subroutine intsol(itask)
  !-----------------------------------------------------------------------
  !****f* outrut/reasol
  ! NAME
  !    reasol
  ! DESCRIPTION
  !    Interprte solver
  ! INPUT
  ! OUTPUT
  ! USES
  ! USED BY
  !    ***_reanut
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only    :  ip
  use def_solver
  use def_master, only    :  kfl_symgr,kfl_schur,kfl_aiipr
  use def_inpout
  implicit none
  integer(ip), intent(in) :: itask
  !
  !.md# Algebraic solver definition
  !.md<code>
  !.md<0>ALGEBRAIC_SOLVER
  !
  if( itask == 1 ) then

     !-------------------------------------------------------------------
     !
     ! Solver
     !
     !-------------------------------------------------------------------
     !
     !.md<1>SOLVER: CG | BICGSTAB | GMRES, KRYLOV=x | DEFLATED_CG, COARSE=SKYLINE | SPARSE 
     !.md<field>SOLVER
     !.md<com>Solver type:
     !.md<com>    - CG:           Conjugate gradient
     !.md<com>    - BICGSTAB:     BiCGStab
     !.md<com>    - GMRES:        GMRES. Krylov dimension should be given
     !.md<com>    - DEFLATED_CG:  Deflated conjugate gradient. Coarse solver can be based skyline or sparse formats. Groups are given in case.dom.dat file.
     !.md<com>
     !.md<com>
   
     if(      words(2) == 'DIREC' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_DIRECT                 ! ALYA   - SPARSE DIRECT (now=14)
     else if( words(2) == 'PASTI' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_PASTIX                 ! PASTIX - SPARSE DIRECT
     else if( words(2) == 'WSMP ' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_WSMP                   ! WSMP   - SPARSE DIRECT
     else if( words(2) == 'PWSMP ' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_PWSMP                  ! PWSMP   - SPARSE DIRECT     
     else if( words(2) == 'MUMPS' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_MUMPS                  ! MUMPS  - SPARSE DIRECT
     else if( words(2) == 'GPUQR' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_GPUQR                  ! ALYA   - GPU DIRECT QR SOLVER
     else if( words(2) == 'CG   ' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_CG                     ! CONJUGATE GRADIENT
     else if( words(2) == 'CONJU' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_CG                     ! CONJUGATE GRADIENT
     else if( words(2) == 'DEFLA' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_DEFLATED_CG            ! DEFLATED CONJUGATE GRADIENT
     else if( words(2) == 'ADEF2' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_A_DEF2                 ! DEFLATED CONJUGATE GRADIENT: A-DEF2 ALGORITHM
     else if( words(2) == 'BICGS' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_BICGSTAB               ! BiCGSTAB
     else if( words(2) == 'GMRES' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_GMRES                  ! GMRES
     else if( words(2) == 'MINRE' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_MINRES                 ! MINRES
     else if( words(2) == 'DIAGO' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_RICHARDSON             ! RICHARDSON WITH MASS
     else if( words(2) == 'RICHA' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_RICHARDSON             ! RICHARDSON
     else if( words(2) == 'MATRI' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_MATRIX_RICHARDSON      ! RICHARDSON WITH MASS AND MATRIX
     else if( words(2) == 'DEFBC' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_DEFLATED_BICGSTAB      ! DEFLATED BiCGSTAB
     else if( words(2) == 'DEFGM' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_DEFLATED_GMRES         ! DEFLATED GMRES
     else if( words(2) == 'SPARS' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_SPARSE_DIRECT          ! SPARSE DIRECT SOLVER
     else if( words(2) == 'STEEP' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_STEEPEST_DESCENT       ! STEEPEST DESCENT
     else if( words(2) == 'JACOB' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_JACOBI                 ! JACOBI
     else if( words(2) == 'PCG  ' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_PIPELINED_CG           ! PIPELINED CONJUGATE GRADIENT
     else if( words(2) == 'PDCG ' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_PIPELINED_DEFLATED_CG  ! PIPELINED DEFLATED CONJUGATE GRADIENT
     else if( words(2) == 'MAPHY' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_MAPHYS_UNSYMMETRIC     ! MAPHYS UNSYMMETRIC
     else if( words(2) == 'OFF  ' ) then
        solve_sol(1) % kfl_algso = SOL_NO_SOLVER                     ! NO SOLVER
     else if( words(2) == 'GCG  ' ) then
        solve_sol(1) % kfl_algso = SOL_GCG                           ! GPU CONJUGATE GRADIENT
     else if( words(2) == 'GCGNP ' ) then
        solve_sol(1) % kfl_algso = SOL_GCGNOPREC                     ! GPU CONJUGATE GRADIENT
     else if( words(2) == 'GGMR ' ) then
        solve_sol(1) % kfl_algso = SOL_GGMRS                         ! GPU GMRES
     else if( words(2) == 'GDECG' ) then
        solve_sol(1) % kfl_algso = SOL_GDECG                         ! GPU DEFCG
     else if( words(2) == 'GPCG ' ) then
        solve_sol(1) % kfl_algso = SOL_GPCG                          ! GPU DEFCG
     else if( words(2) == 'GAMGX' ) then
        solve_sol(1) % kfl_algso = SOL_GAMGX                         ! GPU AMGX
     else if( words(2) == 'AGMG ' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_AGMG                   ! AGMG     
     else if( words(2) == 'PSBLA' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_PSBLAS                 ! PSBLAS
     else if( words(2) == 'EXPLI' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_EXPLICIT               ! Explicit solver 
     else if( words(2) == 'POLYN' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_POLYNOMIAL             ! Polynomial solver 
     else if( words(2) == 'ORTHO' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_ORTHOMIN               ! Orthomin solver
     else if( words(2) == 'GAUSS' ) then
        solve_sol(1) % kfl_algso = SOL_SOLVER_GAUSS_SEIDEL           ! Gauss-Seidel solver
      else if( words(2) == 'PETSC' ) then
#ifdef PETSC
         solve_sol(1) % kfl_algso = SOL_SOLVER_PETSC                 ! PETSC
#else
         call runend('REASOL: PETSC NOT LINKED')
#endif
     else
        call runend('REASOL: NON-EXISTING SOLVER')
     end if
     !
     ! For MAPHYS, determine if solver is unsymmetric, SPD or symmetric
     !
     if( solve_sol(1) % kfl_algso == SOL_SOLVER_MAPHYS_UNSYMMETRIC ) then
        if(      exists('UNSYM') ) then
           solve_sol(1) % kfl_algso = SOL_SOLVER_MAPHYS_UNSYMMETRIC
        else if( exists('SYMME') ) then
           solve_sol(1) % kfl_algso = SOL_SOLVER_MAPHYS_SYMMETRIC
        else if( exists('SPD  ') ) then
           solve_sol(1) % kfl_algso = SOL_SOLVER_MAPHYS_SPD
        else
           call runend('REASOL: CHOOSE AN OPTION FOR MAPHYS UNSYMMETRIC/SPD/SYMMETRIC')
        end if
     end if
     !
     ! Explicit solver options
     !
     if( solve_sol(1) % kfl_algso == SOL_SOLVER_EXPLICIT ) then
        if(      words(3) == 'LUMPE' ) then
           solve_sol(1) % kfl_exp_method = SOL_EXPLICIT_LUMPED
        else if( words(3) == 'CONSI' ) then
           solve_sol(1) % kfl_exp_method = SOL_EXPLICIT_CONSISTENT
        else if( words(3) == 'APPRO' ) then
           solve_sol(1) % kfl_exp_method = SOL_EXPLICIT_APPROX_INVERSE
           solve_sol(1) % kfl_exp_order  = getint('ORDER',1_ip,'#Order of approximate inverse')
        end if
     end if
     !
     ! Schur or not Schur
     !
     if( exists('SCHUR') ) then
        solve_sol(1) % kfl_schur = 1
        kfl_schur = max(1_ip,kfl_schur)
        if(exists('CSR  ')) solve_sol(1) % kfl_scaii = 1
        if(exists('SPARS')) solve_sol(1) % kfl_scaii = 1
        if(exists('SKYLI')) solve_sol(1) % kfl_scaii = 0
        if(exists('CHOLE')) solve_sol(1) % kfl_scaii = 0
     end if
     !
     ! Define symmetry of the graph. Force matrice to be always fully assembled
     !
     if( solve_sol(1) % kfl_algso == 1 .or. solve_sol(1) % kfl_algso == 2 ) then
        solve_sol(1) % kfl_symme = 0
     end if
     !
     ! Define symmetry of the equation
     !
     if(  solve_sol(1) % kfl_schur == 1                                .or. &
          solve_sol(1) % kfl_algso == SOL_GCG                          .or. &
          solve_sol(1) % kfl_algso == SOL_GPCG                         .or. &
          solve_sol(1) % kfl_algso == SOL_GDECG                        .or. &
          solve_sol(1) % kfl_algso == SOL_SOLVER_CG                    .or. &
          solve_sol(1) % kfl_algso == SOL_SOLVER_DEFLATED_CG           .or. &
          solve_sol(1) % kfl_algso == SOL_SOLVER_PIPELINED_DEFLATED_CG .or. &
          solve_sol(1) % kfl_algso == SOL_SOLVER_MAPHYS_SYMMETRIC      .or. &
          solve_sol(1) % kfl_algso == SOL_SOLVER_MAPHYS_SPD            .or. &
          solve_sol(1) % kfl_algso == SOL_SOLVER_MINRES                ) then
        solve_sol(1) % kfl_symeq = 1
     else
        solve_sol(1) % kfl_symeq = 0
     end if
     !
     ! Round off error corecctions
     !
     if( exists('ROUND') ) then
        solve_sol(1) % kfl_roe_correction = 2
        if( exists('RESTA') ) then
           solve_sol(1) % kfl_roe_correction = 1
        else if( exists('ORIGI') ) then
           solve_sol(1) % kfl_roe_correction = 2
        end if
     end if
     !
     ! Deflated: define options
     !
     if(  solve_sol(1) % kfl_algso == SOL_SOLVER_DEFLATED_CG           .or. &
          solve_sol(1) % kfl_algso == SOL_GDECG                        .or. &
          solve_sol(1) % kfl_algso == SOL_SOLVER_PIPELINED_DEFLATED_CG ) then
        if( exists('ALLRE') ) solve_sol(1) % kfl_gathe = 0   ! DCG: All reduce
        if( exists('ALLGA') ) solve_sol(1) % kfl_gathe = 1   ! DCG: All gather
        if( exists('SENDR') ) solve_sol(1) % kfl_gathe = 2   ! Send/receive
        if( exists('COARS') ) then
           if( exists('DIREC') ) then
              solve_sol(1) % kfl_defso = 0
           else if( exists('ITERA') ) then
              solve_sol(1) % kfl_defso = 1
           end if
           if( exists('CSR  ') .or. exists('SPARS') ) then
              solve_sol(1) % kfl_defas = 1
           else if( exists('SKYLI') .or. exists('CHOLE') ) then
              solve_sol(1) % kfl_defas = 0
           else if( exists('DENSE') ) then
              solve_sol(1) % kfl_defas = 2
           end if
           if( solve_sol(1) % kfl_defso == 1 .and. solve_sol(1) % kfl_defas /= 2 ) then
              solve_sol(1) % kfl_defas = 1
           end if
        end if
     end if
     !
     ! Complex solver
     !
     if( exists('COMPL') ) solve_sol(1) % kfl_cmplx = 1
     !
     ! Deflated CG
     !
     if(exists('GROUP')) then
        call runend('GROUPS OF DEFALTED SHOULD BE DECLARED IN MASTER')
     end if
     !
     ! GMRES
     !
     if(  solve_sol(1) % kfl_algso == SOL_SOLVER_GMRES              .or. &
          solve_sol(1) % kfl_algso == SOL_SOLVER_MAPHYS_UNSYMMETRIC .or. &
          solve_sol(1) % kfl_algso == SOL_SOLVER_AGMG               .or. &
          solve_sol(1) % kfl_algso == SOL_SOLVER_PSBLAS             .or. &
          solve_sol(1) % kfl_algso == SOL_GGMRS                     ) then
        if(exists('CLASS')) solve_sol(1) % kfl_ortho = 0                                                   ! GMRES: Classical Gram-Schmidt
        if(exists('MODIF')) solve_sol(1) % kfl_ortho = 1                                                   ! GMRES: Modified  Gram-Schmidt
        if(exists('ITERA')) solve_sol(1) % kfl_ortho = 3                                                   ! GMRES: Iterative classical Gram-Schmidt
        solve_sol(1) % nkryd = getint('KRYLO',10_ip,'#Solver Krylov dimension')
     end if
     !
     ! Symmetric graph
     !
     if( solve_sol(1) % kfl_symme == 1 ) kfl_symgr = 1

  else if( itask == 2 ) then

     !-------------------------------------------------------------------
     !
     ! Preconditioner
     !
     !-------------------------------------------------------------------
     !
     !.md<1>PRECONDITIONER: DIAGONAL | LINELET | RAS | AS
     !.md<field>PRECONDITIONER
     !.md<com>Preconditioner
     !.md<com>    - DIAGONAL:     simple diagonal preconditioner
     !.md<com>    - LINELET:      linelet preconditioner. Options: BLAYER/BOUNDARY, TOLERANCE=10
     !.md<com>      - BOUNDARY/BLAYER: from where linelets are generated
     !.md<com>        - BOUNDARY: means that linelet are created from boundaries.
     !.md<com>        - BLAYER:   means that linelet are created from boundaries and only expend in boundary layer elements (QUA04, HEX08, PEN06, etc.). 
     !.md<com>      - TOLERANCE:  tolerance to create the linelet. A high number means small linelets. A value of 10 is usually taken
     !.md<com>    - RAS and AS:   restricted and classical additive schwarz preconditioners. One RAS or AS subdomain per MPI.
     !.md<com>      - BLOCK:      for equations with mutliple dofs, BLOCK solves each dof independently. Convergence is slower but factorization faster!
     !.md<com>      - MAXIMUM=x:  if you want more RAS subdomains than MPI processes, indicate here the maximum numer of nodes you want per RAS subdomain
     !.md<com>In addition you can prescribe where to impose the preconditioner:
     !.md<com>    - LEFT/RIGHT: choose right or left preconditioning. For symmetric solvers (CG, Deflated CG), this option should not be given as the preconditioning is always implemented as a left)
     !.md<com>
     !.md<com>

     if( exists('LEFT ') ) then
        solve_sol(1) % kfl_leftr = 0
     else if( exists('RIGHT') ) then
        solve_sol(1) % kfl_leftr = 1
     end if

     if( exists('THRES') ) then
        solve_sol(1) % threshold = getrea('THRES',1.0e-3_rp,'#Threshold')
     end if
     
     if(exists('SQUAR')) then
        solve_sol(1) % kfl_preco = 1
        solve_sol(1) % wprec     = 'LEFT-RIGHT SYMMETRIC: sqrt[|diag(A)|]'
        !call runend('REASOL: NOT AVAILABLE')

     else if(exists('DIAGO')) then
        solve_sol(1) % kfl_preco = SOL_DIAGONAL
        solve_sol(1) % wprec     = 'LEFT: DIAGONAL diag(A)'

     else if(exists('MATRI')) then
        solve_sol(1) % kfl_preco = SOL_MATRIX
        solve_sol(1) % wprec     = 'LEFT: MATRIX'

     else if(exists('LINEL').and. (.not.exists('MULTI'))) then
        solve_sol(1) % kfl_preco = SOL_LINELET
        solve_sol(1) % wprec     = 'LINELET'
        if( exists('TOLER') ) then
           solve_sol(1) % toler = getrea('TOLER',10.0_rp,'#Linelet tolerance')
        end if
        if(      exists('BOUND') ) then
           solve_sol(1) % kfl_linty = 1
        else if( exists('BLAYE') ) then
           solve_sol(1) % kfl_linty = 3
        else if( exists('INTER') ) then
           solve_sol(1) % kfl_linty = 0
        else if( exists('PRESC') ) then
           solve_sol(1) % kfl_linty = 2
        else if( exists('FIELD') ) then
           solve_sol(1) % kfl_linty = -getint('FIELD',-1_ip,'#Field for linelet')
        else if( exists('NATUR') ) then
           solve_sol(1) % kfl_linty = -999
        end if

     else if(exists('MASSM')) then
        solve_sol(1) % kfl_preco = SOL_MASS_MATRIX
        solve_sol(1) % wprec     = 'MASS MATRIX (DIAGONAL)'

     else if(exists('GAUSS')) then
        solve_sol(1) % kfl_preco = SOL_GAUSS_SEIDEL
        if( exists('SYMME') ) then
           solve_sol(1) % kfl_renumbered_gs = -1
           solve_sol(1) % wprec             = 'SYMMETRIC GAUSS-SEIDEL'
        else if( exists('STREA') ) then
           solve_sol(1) % angle_stream      =  getrea('ANGLE',0.0_rp,'Angle of streamwise')
           solve_sol(1) % kfl_renumbered_gs =  1
           solve_sol(1) % wprec             = 'STREAMWISE GAUSS-SEIDEL'
        else
           solve_sol(1) % kfl_renumbered_gs =  0
           solve_sol(1) % wprec             = 'GAUSS-SEIDEL'
        end if

     else if(exists('BIDIAG')) then
        solve_sol(1) % kfl_preco = SOL_BIDIAGONAL
        if( exists('STREA') ) then
           if( exists( 'ANGLE' ) )solve_sol(1) % angle_stream = getrea('ANGLE',0.0_rp,'Angle of streamwise')
           solve_sol(1) % kfl_renumbered_gs = 1
           solve_sol(1) % wprec             = 'STREAMWISE BIDIAGONAL'
        else
           solve_sol(1) % kfl_renumbered_gs = 0
           solve_sol(1) % wprec             = 'BIDIAGONAL'
        end if

     else if(exists('SSOR ')) then
        solve_sol(1) % kfl_preco            = SOL_GAUSS_SEIDEL
           solve_sol(1) % kfl_renumbered_gs = -1
           solve_sol(1) % wprec             = 'LEFT: SYMMETRIC GAUSS-SEIDEL'

     else if(exists('CLOSE')) then
        solve_sol(1) % kfl_preco = SOL_CLOSE_MASS_MATRIX
        solve_sol(1) % wprec     = 'CLOSE MASS MATRIX (DIAGONAL)'

     else if(exists('RICHA')) then
        solve_sol(1) % kfl_preco = SOL_RICHARDSON
        solve_sol(1) % wprec     = 'LEFT: RICHARDSON'

     else if(exists('ORTHO')) then
        solve_sol(1) % kfl_preco = SOL_ORTHOMIN
        solve_sol(1) % wprec     = 'LEFT: ORTHOMIN(1)'

     else if(exists('OFF  ').or.exists('NONE ')) then
        solve_sol(1) % kfl_preco = 0
        solve_sol(1) % wprec     = 'NONE'

     else if(exists('APPRO')) then
        solve_sol(1) % kfl_preco = SOL_APPROXIMATE_SCHUR
        solve_sol(1) % wprec     = 'LEFT: APPROXIMATE SCHUR= ABB - ABI diag(AII)^{-1} AIB'
        kfl_schur                = max(3_ip,kfl_schur)

     else if(exists('ABB  ')) then
        solve_sol(1) % kfl_preco = SOL_ABB
        solve_sol(1) % wprec     = 'LEFT: ABB'
        kfl_schur                = max(2_ip,kfl_schur)
        if(exists('CSR  ')) solve_sol(1) % kfl_scpre = 1
        if(exists('SPARS')) solve_sol(1) % kfl_scpre = 1
        if(exists('SKYLI')) solve_sol(1) % kfl_scpre = 0
        if(exists('CHOLE')) solve_sol(1) % kfl_scpre = 0

     else if(exists('MODDI')) then
        solve_sol(1) % kfl_preco = SOL_MOD_DIAGONAL
        solve_sol(1) % wprec     = 'LEFT: diag(ABB) - diag( ABI diag(AII)^{-1} AIB)'
        kfl_schur                = max(2_ip,kfl_schur)

     else if(exists('AII  ')) then
        if( exists('COARS') ) then
           solve_sol(1) % kfl_preco = SOL_COARSE_AII
           solve_sol(1) % wprec     = 'LEFT: COARSE AII^{-1}'
        else
           solve_sol(1) % kfl_preco = SOL_AII
           solve_sol(1) % wprec     = 'LEFT: AII^{-1}'
           kfl_aiipr                = 1
        end if

     else if(exists('AUGME')) then
        solve_sol(1) % kfl_preco = SOL_AUGMENTED_DIAGONAL
        solve_sol(1) % wprec     = 'AUGMENTED DIAGONAL'

     else if(exists('MULTI')) then
        solve_sol(1) % kfl_preco = SOL_MULTIGRID
        solve_sol(1) % wprec     = 'MULTIGRID DEFLATED'
        if( exists('DIREC') ) then
           solve_sol(1) % kfl_defso = 0
        else
           solve_sol(1) % kfl_defso = 1
           solve_sol(1) % kfl_defas = 1
        end if
        if( exists('CSR  ') .or. exists('SPARS') ) then
           solve_sol(1) % kfl_defas = 1
        else if( exists('SKYLI') .or. exists('CHOLE') ) then
           solve_sol(1) % kfl_defas = 0
        end if
        if( solve_sol(1) % kfl_defso == 1 ) then
           solve_sol(1) % kfl_defas = 1
        end if
        if( exists('LINEL') ) then
           solve_sol(1) % kfl_defpr = 4
           solve_sol(1) % toler     = getrea('TOLER',10.0_rp,'#Linelet tolerance')
           if( exists('BOUND') ) then
              solve_sol(1) % kfl_linty = 1
           else if( exists('INTER') ) then
              solve_sol(1) % kfl_linty = 0
           else if( exists('PRESC') ) then
              solve_sol(1) % kfl_linty = 2
           else if( exists('FIELD') ) then
              solve_sol(1) % kfl_linty = -getint('FIELD',-1_ip,'#Field for linelet')
           end if
        else if( exists('DIAGO') .or.  exists('JACOB') ) then
           solve_sol(1) % kfl_defpr = 2
        else if( exists('GAUSS') ) then
           solve_sol(1) % kfl_defpr = 6
        else
           solve_sol(1) % kfl_defpr = 2
        end if

     else if(exists('NEUMA')) then
        solve_sol(1) % kfl_preco = SOL_NEUMANN
        solve_sol(1) % wprec     = 'NEUMANN'

     else if(exists('LOCDI')) then
        solve_sol(1) % kfl_preco = SOL_LOCAL_DIAGONAL
        solve_sol(1) % wprec     = 'LOCAL DIAGONAL: DIAGONAL WITH DIFFERENT ENTRIES'

     else if(exists('RAS  ') .or. exists('AS   ') .or. exists('ADDIT') ) then
        solve_sol(1) % kfl_preco = SOL_RAS
        solve_sol(1) % wprec     = 'RAS: RESTRICTED ADDITIVE SCHWARZ'
        if( exists('AS   ') .or. exists('ADDIT') ) solve_sol(1) % kfl_add_schwarz = 1
        if( exists('BLOCK') ) then
           solve_sol(1) % kfl_block_ras = 1
        else
           solve_sol(1) % kfl_block_ras = 0
        end if
        if( exists('MAXIM') ) then
           solve_sol(1) % max_dof_ras = getint('MAXIM',1_ip,'#Max number of dofs per subdomain')
        end if
        
     else if(exists('BLOCK')) then
        solve_sol(1) % kfl_preco = SOL_BLOCK_DIAGONAL
        solve_sol(1) % wprec     = 'BLOCK DIAGONAL'

     end if
     !
     ! Always clean preconditioner. By default, preconditioner is only computed/factorized
     ! if matrix has changed
     !
     if( exists('ALWAY') ) then
        solve_sol(1) % kfl_clean_precond = 1
     else if( exists('NEVER') ) then
        solve_sol(1) % kfl_clean_precond = 0
     end if
     !
     ! For CG solvers, do not permit right preconditioning
     !
     if(    solve_sol(1) % kfl_algso == SOL_SOLVER_CG                    .or. &
          & solve_sol(1) % kfl_algso == SOL_SOLVER_DEFLATED_CG           .or. &
          & solve_sol(1) % kfl_algso == SOL_SOLVER_PIPELINED_CG          .or. &
          & solve_sol(1) % kfl_algso == SOL_SOLVER_PIPELINED_DEFLATED_CG ) then
        if( solve_sol(1) % kfl_leftr == SOL_RIGHT_PRECONDITIONING ) then
           solve_sol(1) % kfl_leftr = SOL_LEFT_PRECONDITIONING
        end if
     end if

  else if( itask == 3 ) then

     !-------------------------------------------------------------------
     !
     ! Output
     !
     !-------------------------------------------------------------------
     !
     !.md<1>OUTPUT: CONVERGENCE | NON_PRECONDITIONED
     !.md<field>OUTPUT
     !.md<com>Outputs
     !.md<com>    - CONVERGENCE:  output convergence in case-*.cso file. This could be time consuming...
     !.md<com>    - FORCE:        NON_PRECONDITIONED`: output the non-preconditioned residual in case-*.cso file.
     !.md<com>                    This implies one extra SpMV unless you are using right preconditioning. Do not abuse!
     !.md<com>    - MARKET:       Output matrices in matrix market format
     !.md<com>

     if(exists('CONVE')) solve_sol(1) % kfl_cvgso = 1                                                   ! Convergence
     if(exists('MARKE')) solve_sol(1) % kfl_marke = 1                                                   ! Output matrix in market format
     if(exists('MATRI')) solve_sol(1) % kfl_marke = 2                                                   ! Output matrix in ASCII format
     if(exists('NONPR')) solve_sol(1) % kfl_exres = 1                                                   ! Output non-preconditioned residual
     if(exists('OUTPU')) solve_sol(1) % kfl_solve = 1                                                   ! Output solver information

  else if( itask == 4 ) then

     !-------------------------------------------------------------------
     !
     ! Convergence
     !
     !-------------------------------------------------------------------

     if(exists('ADAPT')) solve_sol(1) % kfl_adres = 1                                                   ! Adaptive tolerance
     if(exists('TOLER')) solve_sol(1) % solco = getrea('TOLER',0.01_rp,'#Solver tolerance')             ! Tolerance
     if(exists('ITERA')) solve_sol(1) % miter = getint('ITERA',1_ip, '#Solver number of iterations')    ! Max, # iterations
     if(exists('RATIO')) solve_sol(1) % adres = getrea('RATIO',0.1_rp, '#Ratio wr to initial residual') ! Adpative residual
     solve_sol(1) % solmi = solve_sol(1) % solco                                                        ! Minimum tolerance

  else if( itask == 5 ) then

     !-------------------------------------------------------------------
     !
     ! Options
     !
     !-------------------------------------------------------------------
     !
     !.md<1>OPTIONS: FIXITY | FORCE
     !.md<field>OPTIONS
     !.md<com>Options
     !.md<com>    - FIXITY:       let the solver impose your Dirichlet conditions (useful!). Your module shoudl handle this.
     !.md<com>    - ZERO_FIXITY:  let the solver impose your zero Dirichlet conditions (useful!). Your module shoudl handle this.
     !.md<com>    - FORCE:        force continuity of the solution after the solver (if you think round-off errors have been large)
     !.md<com>    - FULL_ROW:     use the full row format (rectangular matrices). Alya uses by default the partial row format (square matrices)
     !.md<com>    - REACTION:     the solver computes the algebraic reaction on Dirichlet nodes
     !.md<com>    - RECOVER_RHS:  by default, the RHS is exchanged (assembled on interface nodes) at the end aof the solver.
     !.md<com>                    To recover the original one, use this option
     !.md<com>    - SAVE_KRYLOV:  save Krylov subspace to project new initial solution at next iterations. To restrict the dimension
     !.md<com>                    of the subspace to x, use SIZE=x.
     !.md<com>
     !.md<com>

     if(exists('UNSYM')) solve_sol(1) % kfl_symme = 0                                                   ! Force non-symmetric assembly
     if(exists('NONSY')) solve_sol(1) % kfl_symme = 0                                                   ! Force non-symmetric assembly
     if(exists('SYMME')) solve_sol(1) % kfl_symme = 1                                                   ! Force symmetric assembly
     if(exists('PENAL')) then
        if( exists('DUALT') ) solve_sol(1) % kfl_penal = 1
        solve_sol(1) % penal = getrea('PARAM',1.0_rp, '#Penalization parameter')
     end if
     if(exists('FORCE')) solve_sol(1) % kfl_force = 1                                                   ! Force continuity after solver
     if(exists('LIMIT')) solve_sol(1) % kfl_limit = 1                                                   ! Algebraic limiter
     if(exists('CSR  ')) solve_sol(1) % kfl_defas = 1                                                   ! DCG: Sparse assembly
     if(exists('ALLRE')) solve_sol(1) % kfl_gathe = 0                                                   ! DCG: All reduce
     if(exists('ALLGA')) solve_sol(1) % kfl_gathe = 1                                                   ! DCG: All gather
     if(exists('SENDR')) solve_sol(1) % kfl_gathe = 2                                                   ! Send/receive
     if(exists('FIXIT')) solve_sol(1) % kfl_iffix = 1                                                   ! Fixity imposed by solver
     if(exists('ZEROF')) solve_sol(1) % kfl_iffix = 2                                                   ! Fixity imposed by solver with zero value
     if(exists('UNKNO')) solve_sol(1) % kfl_iffix = 3                                                   ! Use unknown to prescribe Dirichlet
     if(exists('ONLYE')) solve_sol(1) % kfl_dirichlet = 2                                               ! Only eliminate Dirichlet, do not enforce it
     if(exists('MODIF')) solve_sol(1) % kfl_ortho = 1                                                   ! GMRES: Modified  Gram-Schmidt
     if(exists('CLASS')) solve_sol(1) % kfl_ortho = 0                                                   ! GMRES: Classical Gram-Schmidt
     if(exists('SCHUR')) solve_sol(1) % kfl_schum = 1                                                   ! If matrix comes from Schur complement
     if(exists('REACT')) solve_sol(1) % kfl_react = 2                                                   ! Compute reaction at Dirichlet nodes
     if(exists('ALLRE')) solve_sol(1) % kfl_react = 3                                                   ! Compute reaction on all nodes
     if(exists('FULLR')) solve_sol(1) % kfl_full_rows = 1                                               ! Full-row assembly
     if(exists('RECOV')) solve_sol(1) % kfl_recov = 1                                                   ! Recover RHS
     if(exists('SAVEK')) then
        solve_sol(1) % kfl_save_krylov = -1                                                             ! Save Krylov subspace
        if( exists('SIZE ') ) then
           solve_sol(1) % kfl_save_krylov = getint('SIZE ',1_ip,'#Size of saved Krylov subspace')
        end if
     end if
     !
     ! Level the reaction
     !
     if( ( exists('REACT') .or. exists('ALLRE') ) .and. exists('LEVEL') ) then
        solve_sol(1) % reaction_level = getrea('LEVEL',0.0_rp,'#Reaction level')
     end if
     !
     ! Example:    ndofn = 6, block = 3, [ u1 u2 u3 u4 u5 u6 ] =>
     ! Contiguous:   [ u1 u2 ] [ u3 u4 ] [ u5 u6 ]
     ! Uncontiguous: [ u1 u3 ] [ u2 u5 ] [ u3 u6 ]
     !     
     if(exists('BLOCK')) then
        solve_sol(1) % kfl_blogs = 1                                                                    ! Block Gauss-Seidel treatment
        if( exists('NUMBL') ) then
           solve_sol(1) % nblok = getint('NUMBL',1_ip,'#Number of blocks Gauss-Seidel')
           if( exists('UNCON') ) then
              solve_sol(1) % nblok = -abs(solve_sol(1) % nblok)
           else if( exists('CONTI') ) then
              solve_sol(1) % nblok =  abs(solve_sol(1) % nblok)
           end if
        end if
     end if

  else if( itask == 6 ) then

     !-------------------------------------------------------------------
     !
     ! Smoother: only for multigrid solver
     !
     !-------------------------------------------------------------------

     if( exists('ITERA') ) then
        solve_sol(1) % itpre = getint('ITERA',1_ip,'#Jacobi iterations')
     end if
     if( exists('LINEL') ) then
        solve_sol(1) % kfl_defpr = 4
        solve_sol(1) % toler     = getrea('TOLER',10.0_rp,'#Linelet tolerance')
        if( exists('BOUND') ) then
           solve_sol(1) % kfl_linty = 1
        else if( exists('INTER') ) then
           solve_sol(1) % kfl_linty = 0
        else if( exists('PRESC') ) then
           solve_sol(1) % kfl_linty = 2
        else if( exists('FIELD') ) then
           solve_sol(1) % kfl_linty = -getint('FIELD',-1_ip,'#Field for linelet')
        end if
     else if( exists('DIAGO') ) then
        solve_sol(1) % kfl_defpr = 2
     else if( exists('GAUSS') ) then
        solve_sol(1) % kfl_defpr = 6
     else if( exists('BIDIAG') ) then
        solve_sol(1) % kfl_defpr = 24
     else
        solve_sol(1) % kfl_defpr = 2
     end if

  else if( itask == 7 ) then

     !-------------------------------------------------------------------
     !
     ! Coarse solver
     !
     !-------------------------------------------------------------------

     if( words(2) == 'ON   ' .or. words(2) == 'YES  ' ) then
        solve_sol(1) % kfl_coarse = 1
        solve_sol(1) % kfl_defas  = 1
     end if

  else if( itask == 8 ) then

     !-------------------------------------------------------------------
     !
     ! Residual normalization
     ! Decide on how to normalize the residual norm:
     !       || b-Ax ||
     ! eps = ----------
     !         || n ||
     !
     ! This can be useful if the user wants to impose the normalization
     ! criteria. By default, || n || = || b ||.
     ! But in the case of using Jacobian (Newton-Raphson), || b ||
     ! tends to zero and another criteria is better.
     !
     ! For example:  A Dx = b, x^{i+1} = Dx + x^{i}
     ! => then one can choose || n || = A  x^{i}
     !
     !-------------------------------------------------------------------

     if(      words(2) == 'RHS  ' ) then
        solve_sol(1) % kfl_normalization = 0
     else if( words(2) == 'USER ' ) then
        solve_sol(1) % kfl_normalization = 1
     else if( words(2) == 'CONST' ) then
        solve_sol(1) % kfl_normalization = 2
        solve_sol(1) % normalization = getrea('VALUE',1.0_rp,'#Residual normalization value')
     end if

  else if( itask == 9 ) then

     !-------------------------------------------------------------------
     !
     ! Configuration file
     !
     !-------------------------------------------------------------------

     call getbigcha(solve_sol(1) % conf_file,'CONFI','NULL ','#Configuration file')

  else if( itask == 10 ) then

     !-------------------------------------------------------------------
     !
     ! RHS minimum norm
     !
     !-------------------------------------------------------------------

     solve_sol(1) % bnorm_min = getrea('RHSMI',1.0e-12_rp,'#RHS minimum norm')

  else if( itask == 11 ) then

     !-------------------------------------------------------------------
     !
     ! Matrix format: CSR, COO, ELL
     !
     !-------------------------------------------------------------------

     if(      words(2) == 'CSR  ' ) then
        solve_sol(1) % kfl_format = 1
     else if( words(2) == 'COO  ' ) then
        solve_sol(1) % kfl_format = 2
     else if( words(2) == 'ELL  ' ) then
        solve_sol(1) % kfl_format = 3
     end if
     
  else if( itask == 12 ) then

     !-------------------------------------------------------------------
     !
     ! OpenMP options
     !
     !-------------------------------------------------------------------

     if(      words(2) == 'OFF  ' ) then
        !
        ! OMP off
        !
        solve_sol(1) % omp_schedule = SOL_OMP_OFF
        
     else if(  words(2) == 'STATI' ) then
        !
        ! OMP with static scheduling
        !
        solve_sol(1) % omp_schedule = SOL_OMP_STATIC
        
     else if( words(2) == 'GUIDE' ) then
        !
        ! OMP with guided scheduling
        !
        solve_sol(1) % omp_schedule = SOL_OMP_GUIDED
        
     else if( words(2) == 'DYNAM' ) then
        !
        ! OMP with dynamic scheduling
        !
        solve_sol(1) % omp_schedule = SOL_OMP_DYNAMIC
        if( exists('CHUNK') ) then
           !
           ! Chunk size: prescribed or auto-tuning
           !
           if( exists('AUTOT') .or. exists('AUTOM') ) then
              solve_sol(1) % omp_chunk_size = -1
           else
              solve_sol(1) % omp_chunk_size = getint('CHUNK',1000_ip,'#OMP chunk size')
           end if
        else
           call runend('REASOL: CHUNK SIZE IS MANDATORY WHEN USING OMP DYNAMIC SCHEDULE')
        end if
        
     else if( words(2) == 'INTER' ) then
        !
        ! Interface chunks. 0 to desactivate OMP
        !
        solve_sol(1) % omp_interface = getint('INTER',0_ip,'#OMP chunk size for interface')
        
     else if( words(2) == 'AUTOT' ) then
        !
        ! Auto-tuning: select automatically scheduling and chunk size
        !
        solve_sol(1) % omp_chunk_size = -2
     else
        call runend('REASOL: UNDEFINED OMP SCHEDULE')
     end if

  else if( itask == 13 ) then

     !-------------------------------------------------------------------
     !
     ! Block solver
     !
     !-------------------------------------------------------------------

     if( words(2) == 'ON   ' .or. words(2) == 'YES  ' ) solve_sol(1) % kfl_block = 1
     if( exists('SCHUR') ) then
        solve_sol(1) % kfl_block_schur = getint('SCHUR',2_ip,'#Solver based on shcur complement')
     end if
     
  else if( itask == 14 ) then

     !-------------------------------------------------------------------
     !
     ! Relaxation
     !
     !-------------------------------------------------------------------
     !
     !.md<1>RELAXATION= real
     !.md<field>RELAXATION
     !.md<com>RELAXATION
     !.md<com>    Relax the unknown in the solver. Only available for very few solvers
     !.md<com>

     if( exists('RELAX') ) then
        solve_sol(1) % relax = getrea('RELAX',1.0_rp,'#Relaxation parameter')
     end if
     
  end if
  !
  !.md<0>END_ALGEBRAICS_SOLVER
  !.md</code>
  !

end subroutine intsol

