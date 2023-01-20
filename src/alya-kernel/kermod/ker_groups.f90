!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ker_groups()

  use def_parame
  use def_elmtyp
  use def_master
  use def_kermod
  use def_domain
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_deallo
  use def_solver, only : solve_sol 
  use mod_solver, only : solver_define_groups_from_field
  use mod_solver, only : solver_initialize_matrix_and_rhs
  use mod_solver, only : solver_solve
  use mod_ADR,    only : ADR_assemble_laplacian
  use mod_opebcs, only : boundary_conditions_impose_node_codes
  implicit none

  integer(ip), pointer :: kfl_fixno_loc(:,:)
  real(rp),    pointer :: bvess_loc(:,:)
  
  if( ngrou_dom < -4 ) then

     nullify(kfl_fixno_loc)
     nullify(bvess_loc)
     
     solve_sol => momod(ID_KERMOD) % solve(6:)

     if( INOTMASTER ) then
        call memory_alloca(mem_modul(1:2,modul),'FIXNO','ker_groups',solve_sol(1) % kfl_fixno,1_ip,npoin)
        call memory_alloca(mem_modul(1:2,modul),'BVESS','ker_groups',solve_sol(1) % bvess,1_ip,npoin)
        call  boundary_conditions_impose_node_codes(tncod_ker(6),solve_sol(1) % kfl_fixno,solve_sol(1) % bvess)
     end if
     
     call solver_initialize_matrix_and_rhs(solve_sol,amatr,rhsid)
     call ADR_assemble_laplacian(meshe(ndivi),elmar,amatr)
     call solver_solve(solve_sol,amatr,rhsid,unkno)
          call runend('O.K.!')

     call solver_define_groups_from_field(ngrou_dom,npoin,unkno,solve_sol(1) % lgrou)

     call memory_deallo(mem_modul(1:2,modul),'FIXNO','ker_groups',solve_sol(1) % kfl_fixno)
     call memory_deallo(mem_modul(1:2,modul),'BVESS','ker_groups',solve_sol(1) % bvess)

  end if 

end subroutine ker_groups
