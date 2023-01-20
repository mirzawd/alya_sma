!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_inisol()

  use def_kintyp, only: ip, rp
  use def_master, only: solve
  use def_solver

  implicit none
  !
  ! Allocate memory for solvers
  !
  call soldef(-1_ip)
  !
  ! Output flag
  !
  solve(1) % kfl_solve = 1
  !
  ! Consider unsymmetric matrix
  !
  solve(1) % kfl_symme = 0
  !
  ! Number of degrees of freedom per node/edge
  !
  solve(1) % ndofn = 1
  !
  ! Solver name
  !
  solve(1) % wprob = 'MAGNET'
  !
  ! Degrees of freedom in element edges
  !
  solve(1) % kfl_where = SOL_EDGES
  solve(1) % kfl_cvgso = 1
  !
  ! Direct solver: SOL_SOLVER_MUMPS, SOL_SOLVER_SPARSE_DIRECT (no funciona)
  !
  ! Algo no funcion con edge elements para direct solver
  ! Aparece un error GRAPHS_PERMUT_METIS_POSTORDERING
  ! Probablemente esta relacionado con solpre.f90
  !
  ! Conjugate gradient
  !
  solve(1) % kfl_algso = SOL_SOLVER_CG
  !
  ! Preconditioner: SOL_DIAGONAL, SOL_GAUSS_SEIDEL, SOL_NO_PRECOND, SOL_SQUARE, SOL_RASS, SOL_MULTIGRID
  !
  ! Precondicionador diagonal no funciona bien con rectangulos
  ! Se extrae la diagonal de la matrix correctamente?
  ! Check solite/precon.f90
  !
  solve(1) % kfl_preco = SOL_NO_PRECOND
  !    solve(1) % kfl_preco = SOL_DIAGONAL
  !     solve(1) % kfl_preco = SOL_GAUSS_SEIDEL
  !     solve(1) % kfl_renumbered_gs = SYMMETRIC_GAUSS_SEIDEL
  !
  ! Maximum number of iterations
  !
  solve(1) % miter = 2000
  !
  !
  ! Tolerance
  !
  solve_sol(1) % solco = 1.0e-6_rp
  !
end subroutine mag_inisol
