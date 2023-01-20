!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_solite_mod()
  !------------------------------------------------------------------------
  !****f* magnet/mag_solite
  ! NAME
  !    mag_solite
  ! DESCRIPTION
  !    This routine solves one iteration
  ! USES
  !    mag_matrix
  ! USED BY
  !    mag_doiter
  !***
  !------------------------------------------------------------------------

  use def_master
  use def_solver, only: SOL_MATRIX_HAS_CHANGED, SOL_SYSTEM_HAS_BEEN_SOLVED, resi1, resi2
  use def_magnet
  use mod_solver, only : solver_solve, solver_parallel_scalar_product, solver_postprocess
  use mod_communications, only : PAR_INTERFACE_EDGE_EXCHANGE
  use mod_mag_linalg, only: lusolver
  use mod_mag_inpdat

  implicit none

  integer(ip) :: iconstr, jconstr

  real(rp) :: fimod_mag, femod_mag

  call mag_updbcs(ITASK_DOITER)

  a0_mag = bdfode_mag % ai0
  Hpav_mag => bdfode_mag % ypi

  !---------------------------------------------------------
  ! Assemble matrix and RHS
  !---------------------------------------------------------
  call mag_matrix()
  !---------------------------------------------------------


  !---------------------------------------------------------
  ! Correct rhsid and compute residual rmult
  !---------------------------------------------------------
  !
  do iconstr = 1, constr_total
    !
    ! rH' = b - A * H - sum_j lambda_j * C_j = rH - sum_j lambda_j * C_j
    !
    rhsid = rhsid - constr_mag(iconstr) % multc * constr_mag(iconstr) % cnstr
    !
  end do
  !
  !---------------------------------------------------------

  !---------------------------------------------------------
  !
  ! Alya solver_parallel_scalar_product
  !
  ! nequ1 =
  ! nequ2 = start own boundary
  ! nequ3 = end own boundary
  ! Para comprobar que el producto escalar funciona en paralelo
  ! Haz rhsid = a
  ! Producto escalar = (nequ2 - 1) * a**2 + (nequ3 - nequ2 + 1) * (2 * a)**2
  ! Si no se produce el intercambio el resultado sera
  ! Producto escalar = nequ3 * a**2
  !
  ! Update some parameters for the solver
  !
  solve_sol(1) % bnorm_min = 0.0_rp
  solve_sol(1) % kfl_assem = SOL_MATRIX_HAS_CHANGED
  !
  ! Subroutine 'solver_solve' in 'solite/mod_solver.f90'
  ! 'solver_solve' calls 'solver_solve_system' in 'solite/mod_solver.f90'
  ! 'solver_solve_system' calls 'solver' in 'master/solver.f90'
  ! 'solver' calls 'solite' in 'solite/solite.f90'
  ! 'solite' calls 'cgrpls' in 'solite/cgrpls.f90'
  !
  ! Solve linear system!
  ! There is something wrong about the CG when running in parallel
  ! When running in parallel the master manages to leave the CG loop
  ! The max. number of iterations in the master is 0 actually...
  ! The slaves however keep iterating
  ! max. number of iterations = nedge!
  !
  ! Estamos en solite/precon.f90 case(4)
  ! solver_SpMV no esta funcionando correctamente
  ! No calcula correctamente q = A p
  !
  !---------------------------------------------------------
  ! Solve field
  ! dH = An \ rH'
  !---------------------------------------------------------
  !
  ! Initial guess
  !
  unkno = 0.0_rp
  !
  call solver_solve(momod(modul) % solve, amatr, rhsid, unkno)
  !
  lsiter0_mag = solve_sol(1) % iters
  lsresi1_mag = resi1
  lsresi2_mag = resi2
  !
  ! You may check if the system has been solved
  !
  if (solve_sol(1) % kfl_assem /= SOL_SYSTEM_HAS_BEEN_SOLVED) call runend('mag_solite: system could not be solved')
  !
  if (INOTMASTER) dH_mag = unkno
  !
  !--------------------------------------------------------- 
  !
  do iconstr = 1, constr_total
    !------------------------------------------------------------
    !
    ! Solve A * unkno = C_i
    !
    !------------------------------------------------------------
    !
    ! Initial guess
    !
    constr_mag(iconstr) % unkno = 0.0_rp
    !
    ! Linear Solver
    !
    call solver_solve(momod(modul) % solve, amatr, constr_mag(iconstr) % cnstr, constr_mag(iconstr) % unkno)
    !
    ! Number of iterations
    !
    constr_mag(iconstr) % lsite = solve_sol(1) % iters
    !
    if (solve_sol(1) % kfl_assem /= SOL_SYSTEM_HAS_BEEN_SOLVED) call runend('mag_solite: system could not be solved')
    !
    !------------------------------------------------------------
    !
    ! Compute Lagrange multiplier residuals: rmult_i = I_i - C_i * H
    ! C_i values have been exchanged already among different processors
    !
    !------------------------------------------------------------
    !
    call solver_parallel_scalar_product(solve_sol(1), constr_mag(iconstr) % cnstr, He_mag, constr_mag(iconstr) % valui)
    !
    constr_mag(iconstr) % rmult = constr_mag(iconstr) % valuc - constr_mag(iconstr) % valui
    !
    !------------------------------------------------------------
    !
    ! Compute constraint matrix components: C_j A \ C_i = C_j * unkno
    !
    !------------------------------------------------------------
    do jconstr = 1, constr_total
      !
      call solver_parallel_scalar_product(solve_sol(1), constr_mag(jconstr) % cnstr, constr_mag(iconstr) % unkno, cntmat_mag(jconstr, iconstr))
      !
    end do
    !------------------------------------------------------------
    !
    ! Compute constraint rhs components: C_i A \ rH - rmult_i
    !
    !------------------------------------------------------------
    !
    call solver_parallel_scalar_product(solve_sol(1), constr_mag(iconstr) % cnstr, dH_mag, cntrhs_mag(iconstr))
    !
    cntrhs_mag(iconstr) = cntrhs_mag(iconstr) - constr_mag(iconstr) % rmult
    !
  end do
  !--------------------------------------------------------------
  !
  ! Solve constraint system to obtain multd_i
  !
  !--------------------------------------------------------------
  if (constr_total > 0) call lusolver(cntunk_mag, cntmat_mag, cntrhs_mag)
  !--------------------------------------------------------------
  do iconstr = 1, constr_total
    !
    constr_mag(iconstr) % multd = cntunk_mag(iconstr)
    !
    ! Update Lagrange multiplier
    !
    constr_mag(iconstr) % multc = constr_mag(iconstr) % multc + constr_mag(iconstr) % multd
    !
    ! Update dH_mag
    !
    dH_mag = dH_mag - constr_mag(iconstr) % multd * constr_mag(iconstr) % unkno
    !
  end do
  !--------------------------------------------------------------
  if (INOTMASTER) He_mag = He_mag + dH_mag
  !--------------------------------------------------------------
  !
  !--------------------------------------------------------------
  ! Compute residual norms
  !--------------------------------------------------------------
  !
  ! RHS l2-norm
  !
  call PAR_INTERFACE_EDGE_EXCHANGE(solve_sol(1) % ndofn, bhsid, 'SUM', 'IN MY CODE')
  call solver_parallel_scalar_product(solve_sol(1), bhsid, bhsid, femod_mag)
  !
  call PAR_INTERFACE_EDGE_EXCHANGE(solve_sol(1) % ndofn, fhsid, 'SUM', 'IN MY CODE')
  call solver_parallel_scalar_product(solve_sol(1), fhsid, fhsid, fimod_mag)
  !
  refres_mag(1) = max(sqrt(femod_mag), sqrt(fimod_mag))
  !
  ! Residual l2-norm
  !
  call solver_parallel_scalar_product(solve_sol(1), rhsid, rhsid, residu_mag(1))
  residu_mag(1) = sqrt(residu_mag(1))
  !
  ! Lagrange multiplier residuals
  !
  do iconstr = 1_ip, constr_total
    refval_mag(iconstr + 1_ip) = abs(constr_mag(iconstr) % multc)
    delval_mag(iconstr + 1_ip) = abs(constr_mag(iconstr) % multd)
    refres_mag(iconstr + 1_ip) = max( abs(constr_mag(iconstr) % valuc), abs(constr_mag(iconstr) % valui) )
    residu_mag(iconstr + 1_ip) = abs(constr_mag(iconstr) % rmult)
  end do
  !
  ! Solution variation norm
  !
  call solver_parallel_scalar_product(solve_sol(1), dH_mag, dH_mag, delval_mag(1))
  delval_mag(1) = sqrt(delval_mag(1))
  !
  ! Solution norm
  !
  call solver_parallel_scalar_product(solve_sol(1), He_mag, He_mag, refval_mag(1))
  refval_mag(1) = sqrt(refval_mag(1))
  !
end subroutine mag_solite_mod
