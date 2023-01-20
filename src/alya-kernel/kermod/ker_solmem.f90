!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ker_solmem()

  use def_master
  use def_domain
  use def_kermod
  use def_solver
  use mod_maths
  use mod_memory
  use mod_chktyp, only : check_type
  implicit none

  !----------------------------------------------------------------------
  !
  ! Memory for solver
  !
  !----------------------------------------------------------------------

  if( kfl_extro /= 0 ) then  
     solve_sol => solve(1:1)  ! Roughness extension
     call soldef(4_ip)
  end if
  if( kfl_walld /= 0 ) then
     solve_sol => solve(2:2)  ! Roughness extension
     call soldef(4_ip)
  end if
  if( kfl_suppo /= 0 ) then
     solve_sol => solve(3:3)  ! Support geometry for mesh multiplication 
     call soldef(4_ip)
  end if
  if( kfl_walln /= 0 ) then
     solve_sol => solve(5:5)  ! Wall normal
     call soldef(4_ip)
  end if
  ! if( kfl_defor /= 0 ) then
  !    solve_sol => solve(4:4)  ! Mesh deformation
  !    call soldef(4_ip)
  ! end if
  !
  ! Solver fixity
  !
  solve(1) % kfl_fixno => kfl_fixno_rough_ker
  solve(2) % kfl_fixno => kfl_fixno_walld_ker
  solve(3) % kfl_fixno => kfl_fixno_suppo_ker
  solve(5) % kfl_fixno => kfl_fixno_walln_ker

  !if( kfl_conma /= 0 ) then
  !   call memory_alloca(memor_dom,'CMASS','mod_mass_matrix',cmass,nzdom)
  !end if
  !if( kfl_conma_weighted /= 0 ) then
  !   call memory_alloca(memor_dom,'CMASS_WEIGHTED','mod_mass_matrix',cmass_weighted,nzdom)
  !end if


end subroutine ker_solmem
