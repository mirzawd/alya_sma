!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Temper
!> @{
!> @file    tem_solmem.f90
!> @date    14/06/2019
!> @author  Guillaume Houzeaux
!> @brief   Allocate memory
!> @details Allocate memory for solver and variables that should be
!>          reallocated and not redistributed in case of repartitioning
!> @}
!-----------------------------------------------------------------------

subroutine tem_solmem()

  use def_master
  use def_temper
  use def_domain,          only : ndime,npoin
  use mod_memory,          only : memory_alloca
  use mod_tem_explicit,    only : tem_explicit_memory
  use mod_tem_rk_explicit, only : tem_rk_explicit_memory

  implicit none
  !
  ! Explicit
  !
  if(kfl_explicit_tem == 1) then
     call memory_alloca(mem_modul(1:2,modul),'DT_RHO_CP_TEM','tem_memall',dt_rho_cp_tem,npoin)
  end if
  !
  ! GRTEM: Temperature gradients
  !
  if(kfl_ellen_tem==-1) then
     call memory_alloca(mem_modul(1:2,modul),'GRTEM_TEM','tem_memall',grtem_tem,ndime,nunkn_tem)
  end if
  !
  ! Water vapor concentration gradients
  !
  call memory_alloca(mem_modul(1:2,modul),'GRADC_TEM','tem_memall',gradc_tem,ndime,nunkn_tem)
  !
  ! Solver memory
  !
  solve_sol => solve(1:)
  call soldef(4_ip)
  !
  ! Subrgid scale
  !
  tesgs => ADR_tem % sgs
  !
  ! Solver fixity
  !
  solve(1) % bvess     => bvess_tem(:,:,1)
  solve(1) % kfl_fixno => kfl_fixno_tem
  !
  ! Consistent mass matrix fixity
  !
  solve(2) % kfl_fixno => kfl_fixno_tem
  !
  ! Solution strategies
  !
  if (kfl_explicit_tem == 1) then
     if(kfl_tisch_tem==3) then 
        call tem_explicit_memory()
     else
        call tem_rk_explicit_memory()
     end if
  else 
     continue
  end if
  
end subroutine tem_solmem
