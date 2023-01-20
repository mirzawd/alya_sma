!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




!-----------------------------------------------------------------------
!> @addtogroup Solvers
!> @{
!> @file    def_all_solvers.f90
!> @author  houzeaux
!> @date    2022-03-01
!> @brief   All solvers
!> @details All solvers
!>         
!-----------------------------------------------------------------------

module def_all_solvers

  use def_solvers,              only : solver
  use def_direct_solvers,       only : direct_solver
  use def_direct_solvers,       only : direct_factorization
  use def_iterative_solvers,    only : iterative_solver
  ! Krylov solvers
  use def_cg,                   only : cg
  use def_dcg,                  only : dcg
  use def_bicgstab,             only : bicgstab
  use def_gmres,                only : gmres
  use def_minres,               only : minres
  ! Stationary solvers
  use def_jacobi,               only : jacobi
  use def_richardson,           only : richardson
  use def_ssor,                 only : ssor
  use def_sor,                  only : sor
  ! DDM solvers
  use def_ras,                  only : ras
  ! Direct solvers
  use def_gauss_elimination,    only : gauss_elimination
  use def_approx_inv,           only : approx_inv
  use def_linelet,              only : linelet
  use def_lu_factorization,     only : lu_factorization
  use def_cholesky,             only : cholesky
  use def_solver
  implicit none

end module def_all_solvers
!> @}
