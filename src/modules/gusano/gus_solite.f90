!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    gus_solite.f90
!> @author  houzeaux
!> @date    2020-10-20
!> @brief   Solve
!> @details Solve Gusano
!> @} 
!-----------------------------------------------------------------------

subroutine gus_solite()

  use def_elmtyp
  use def_domain
  use def_master
  use def_gusano
  use mod_gus_element_operations,     only : gus_element_operations
  use mod_gus_boundary_operations,    only : gus_boundary_operations
  use mod_gus_projections,            only : gus_projections_solve
  use mod_gus_projections,            only : gus_projections_initialization
  use mod_solver,                     only : solver_initialize_matrix_and_rhs
  use mod_solver,                     only : solver_solve
  use mod_gus_transmssion_conditions, only : gus_transmssion_conditions
  use mod_gus_bending,                only : gus_bending
  use mod_memory,                     only : memory_size
  use mod_solver,                     only : solver_preprocess_dirichlet_1by1
  implicit none
  integer(ip) :: iz
  
  itinn(modul) = itinn(modul) + 1
  !
  ! Projections
  !
  call gus_projections_initialization()
  !
  ! Initialize solver matrix and RHS
  !
  call solver_initialize_matrix_and_rhs(solve,amatr,rhsid)
  do iz = 1,memory_size(schur_gus)
     schur_gus(iz) = 0.0_rp
  end do
  !
  ! Assemble equations
  !
  call gus_element_operations()
  !
  ! Impost Neumann condition
  !
  call gus_boundary_operations()
  !
  ! Assemble Dirichlet transmission conditions
  !
  call gus_transmssion_conditions(amatr,rhsid)
  !
  ! Bending pressure loss
  !
  call gus_bending(amatr)
  !
  ! Boundary conditions
  !
  if( kfl_algor_gus == GUS_SCHUR_COMPLEMENT .and. memory_size(schur_gus) > 0 ) then
     call solver_preprocess_dirichlet_1by1('MATRIX',solve(3),1_ip,schur_gus,FIXNO=kfl_fixsc_gus)
  end if
  !
  ! Solve system
  !
  !call solver_solve(momod(modul) % solve,amatr,rhsid,unkno,pmatr)
  block
    use mod_algebraic_solver
    call algebraic_solver(momod(modul) % solve(1),amatr,rhsid,unkno,schur_gus)
  end block
  !print*,'unkn=',unkno(12:20)
  !
  ! Solve projections
  !
  call gus_projections_solve()
 
end subroutine gus_solite
