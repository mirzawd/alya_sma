!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_solver.f90
!> @date    14/06/2019
!> @author  Guillaume Houzeaux
!> @brief   Allocate memory
!> @details Allocate memory for solver and variables that should be
!>          reallocated and not redistributed in case of repartitioning
!> @}
!-----------------------------------------------------------------------

subroutine nsi_solmem()
  
  use def_master,               only : INOTEMPTY
  use def_master,               only : solve
  use def_master,               only : modul
  use def_master,               only : mem_modul
  use def_domain,               only : nzrhs
  use def_domain,               only : nzmat
  use def_domain,               only : ndime
  use def_domain,               only : npoin
  use def_domain,               only : nelem 
  use def_domain,               only : ltype 
  use def_domain,               only : ngaus
  use def_kermod,               only : kfl_twola_ker
  use def_kermod,               only : kfl_noslw_ker
  use def_kermod,               only : kfl_waexl_ker
!  use def_solver,               only : solve_sol
  use def_coupli,               only : kfl_dimbou
  use mod_solver,               only : solver_allocate_system_memory
  use mod_memory,               only : memory_alloca
  use mod_nsi_multi_step_fs,    only : nsi_multi_step_fs_memory
  use mod_nsi_fractional_step,  only : nsi_fractional_step_memory
  use mod_nsi_semi_implicit,    only : nsi_semi_implicit_memory
  use mod_nsi_schur_complement, only : nsi_schur_complement_memory
  use mod_output_postprocess,   only : output_postprocess_check_variable_postprocess
  use mod_nsi_algebraic_forces, only : nsi_algebraic_forces_allocate
  use def_coupli,               only : kfl_immer
  use def_coupli,               only : kfl_efect
  use def_nastin
  
  implicit none
  integer(ip) :: ielem,pgaus,pelty,npoin_min,nelem_min

  npoin_min = max(1_ip,npoin)
  nelem_min = max(1_ip,nelem)

  !----------------------------------------------------------------------
  !
  ! Solver
  !
  !----------------------------------------------------------------------
  !
  ! Tell the solver that although we have defined
  ! 2 solvers (momentum and continuity), the
  ! matrix used for these two solvers has a block
  ! structure, owned by the first solver!
  !
  if( NSI_MONOLITHIC ) then
     solve(1) % num_blocks          = 1
     solve(1) % block_dimensions(1) = ndime+1
  else if( NSI_SCHUR_COMPLEMENT ) then
     solve(1) % num_blocks          = 2
     solve(1) % block_dimensions(1) = ndime
     solve(1) % block_dimensions(2) = 1
     solve(1) % block_num           = 1
     solve(2) % block_num           = 2
  else if( NSI_FRACTIONAL_STEP .or. NSI_SEMI_IMPLICIT ) then
     solve(1) % num_blocks          = 2
     solve(1) % block_dimensions(1) = ndime
     solve(1) % block_dimensions(2) = 1
     solve(1) % block_num           = 1
     solve(2) % block_num           = 2
  end if
  !
  ! Memory
  !  solve_sol => solve
  call soldef(4_ip)
  !
  ! Boundary conditions
  !
  solve(1) % bvess     => bvess_nsi(:,:,1)   ! Momentum
  solve(1) % kfl_fixno => kfl_fixno_nsi

  solve(9) % kfl_iffix =  1
  solve(9) % bvess     => bvess_nsi(:,:,1)   ! Viscous term
  solve(9) % kfl_fixno => kfl_fixno_nsi

  if( solve(1) % kfl_iffix /= 0 ) then       ! Velocity correction
     solve(4) % kfl_iffix =  2
     solve(4) % kfl_fixno => kfl_fixno_nsi
     nullify(solve(4) % bvess)
  end if

  solve(2) % kfl_fixno => kfl_fixpr_nsi      ! Pressure schur complement
  nullify(solve(2) % bvess)
  !
  ! Navier-Stokes is block-wise:
  ! First block:  momentum equations
  ! Second block: continuity equation
  !
  solve(1) % block_array(1) % bvess     => bvess_nsi(:,:,1)
  solve(1) % block_array(1) % kfl_fixno => kfl_fixno_nsi
  solve(1) % block_array(1) % bvnat     => solve(1) % bvnat

  solve(1) % block_array(2) % bvess     => bpess_nsi(:,:,1)
  solve(1) % block_array(2) % kfl_fixno => kfl_fixpp_nsi
  solve(1) % block_array(2) % bvnat     => solve(2) % bvnat
  !
  ! Divergence free correction
  !
  if( kfl_divcorrec_nsi /= 0 ) then
     solve(6) % kfl_iffix =  2
     solve(6) % kfl_fixno => kfl_fixno_div_nsi
     nullify(solve(6) % bvess)
  end if
  !
  ! Mass correction  - actually I guess that this lines would not be necesary if the dirichlet bcs are not imposed by the solver
  !
  if( kfl_corre_nsi == 3 ) then
     solve(4) % kfl_iffix =  0               ! not done by the solver
     nullify(solve(4) % bvess)
     solve(4) % num_blocks          = 1
     solve(4) % block_dimensions(1) = ndime
  end if
  if ( kfl_bnods_nsi == 1 ) then
     call memory_alloca(mem_modul(1:2,modul),'IBOUN_NSI','nsi_solmem',iboun_nsi,npoin,2_ip)
  end if
  !
  ! Consistent mass
  !
  solve(8) % bvess     => bvess_nsi(:,:,1)
  solve(8) % kfl_fixno => kfl_fixno_nsi
  !
  ! Schur complement system: modify matrix and RHIS size
  ! Velocity and pressure are consecutively in UNKNO even
  ! if split scheme is used
  !
  if( INOTEMPTY ) then

     if( NSI_SCHUR_COMPLEMENT ) then
        !
        ! Schur complement solver
        !
        nmauu_nsi =  solve(1) % nzmat
        nmaup_nsi =  solve(1) % nzmat/ndime
        nmapu_nsi =  solve(1) % nzmat/ndime
        nmapp_nsi =  solve(2) % nzmat
        poauu_nsi =  1
        poaup_nsi =  poauu_nsi + nmauu_nsi
        poapu_nsi =  poaup_nsi + nmaup_nsi
        poapp_nsi =  poapu_nsi + nmapu_nsi
        nzmat     =  max(nzmat,nmauu_nsi+nmaup_nsi+nmapu_nsi+nmapp_nsi)
        nzrhs     =  max(nzrhs,(ndime+1_ip)*npoin)

     else if( NSI_FRACTIONAL_STEP .or. NSI_SEMI_IMPLICIT ) then
        !
        ! Fractional step: Auu and App are not assembled
        !
        if( NSI_FRACTIONAL_STEP ) then
           nmauu_nsi =  solve(1) % nzmat
        else
           nmauu_nsi =  solve(9) % nzmat
        end if
        nmaup_nsi =  max(solve(1) % nzmat/ndime,solve(2) % nzmat * ndime)
        nmapu_nsi =  max(solve(1) % nzmat/ndime,solve(2) % nzmat * ndime)
        nmapp_nsi =  solve(2) % nzmat
        poauu_nsi =  1
        poaup_nsi =  poauu_nsi + nmauu_nsi
        poapu_nsi =  poaup_nsi + nmaup_nsi
        poapp_nsi =  poapu_nsi + nmapu_nsi
        nzmat     =  max(nzmat,nmauu_nsi+nmaup_nsi+nmapu_nsi+nmapp_nsi)
        nzrhs     =  max(nzrhs,(ndime+1_ip)*npoin)
        !
        ! 
        !
        solve(1) % kfl_iffix = 0
        solve(2) % kfl_iffix = 0
        solve(8) % kfl_iffix = 2
     !   kfl_matdi_nsi        = NSI_DIRICHLET_ALGORITHM 
     end if

  else

     nmauu_nsi =  1
     nmaup_nsi =  1
     nmapu_nsi =  1
     nmapp_nsi =  1
     poauu_nsi =  1
     poaup_nsi =  1
     poapu_nsi =  1
     poapp_nsi =  1

  end if
  !
  ! Dimensions
  !
  ndbgs_nsi = ndime * npoin
  !
  ! Laplacian matrix 
  !
  if( kfl_predi_nsi /= 0 ) then
     call solver_allocate_system_memory(solve(NSI_SOLVER_CONTINUITY), lapla_nsi,MEMORY_COUNTER=mem_modul(1:2,modul))
  end if
  !
  ! Viscous matrix 
  !
  if( NSI_SEMI_IMPLICIT ) then
     call solver_allocate_system_memory(solve(NSI_SOLVER_VISCOUS_TERM),visco_nsi,MEMORY_COUNTER=mem_modul(1:2,modul))
  end if
  !
  ! Consistent mass matrix
  !
  if( kfl_corre_nsi == 3 ) then
     call memory_alloca(mem_modul(1:2,modul),'CMAMA_NSI','nsi_solmem',cmama_nsi,solve(4) % nzmat)
     write(*,*) 'CONTROL LUEGO QUITAR: solve(2) % nzmat, solve(4) % nzmat',solve(2) % nzmat, solve(4) % nzmat
     if(solve(2) % nzmat*ndime*ndime /= solve(4) % nzmat) &
          call runend('NSI_SOLMEM:solve(2) % nzmat*ndime*ndime /= solve(4) % nzmat - PERHAPS you forgot MASSCORRECTION solver') !QUITARLO una vez probado
  end if
  !
  ! Save linear matrix
  !
  if( kfl_savco_nsi == 1 ) then
     call memory_alloca(mem_modul(1:2,modul),'AMATR_NSI','nsi_solmem',amatr_nsi,solve(1) % nzmat)
  end if
  !
  ! Projections of dt/rho and tau
  !
  if( kfl_predi_nsi == 7 .or. kfl_predi_nsi == 8 .or. kfl_predi_nsi == 9 ) then
     call memory_alloca(mem_modul(1:2,modul),'DT_RHO_NSI'  ,'nsi_solmem',dt_rho_nsi  ,npoin_min)
     call memory_alloca(mem_modul(1:2,modul),'TAU_NSI  '   ,'nsi_solmem',tau_nsi     ,npoin_min)
  end if

  !----------------------------------------------------------------------
  !
  ! Solution strategies
  !
  !----------------------------------------------------------------------

  if( NSI_FRACTIONAL_STEP ) then
     !
     ! Fraction step
     !
     if(kfl_tisch_nsi == 4) then 
        call nsi_multi_step_fs_memory()     ! Multi step Fractional Step
     else
        call nsi_fractional_step_memory()   ! Single step Fractional Step    
     end if

  else if( NSI_SEMI_IMPLICIT ) then
     !
     ! Semi-implicit
     !
     call nsi_semi_implicit_memory()
     
  else if( NSI_SCHUR_COMPLEMENT ) then
     !
     ! Schur complement
     !
     call nsi_schur_complement_memory()
     
  end if

  !----------------------------------------------------------------------
  !
  ! Others
  !
  !----------------------------------------------------------------------
  !
  ! Low Mach 
  !
  if( kfl_regim_nsi == 3 ) then
     call memory_alloca(mem_modul(1:2,modul),'DRHODT_NSI','nsi_solmem',drhodt_nsi,npoin)
  end if
  !
  ! Surface tension
  !
  if( kfl_surte_nsi /= 0 ) then
     call memory_alloca(mem_modul(1:2,modul),'NORLE_NSI','nsi_solmem',norle_nsi,ndime,npoin)
     call memory_alloca(mem_modul(1:2,modul),'CURLE_NSI','nsi_solmem',curle_nsi,npoin)
  end if
  !
  ! Hydrostatic level set: used to compute rho_{hyd} to add to NS' RHS
  !
  if( kfl_hydro_gravity_nsi /= 0 ) then
     call memory_alloca(mem_modul(1:2,modul),'HYDRO_DENSITY_NSI','nsi_solmem',hydro_density_nsi,nelem_min)
     do ielem=1,nelem
        pelty = abs(ltype(ielem))
        pgaus = ngaus(pelty)
        call memory_alloca(mem_modul(1:2,modul),'HYDRO_DENSITY_NSI % A','nsi_solmem',hydro_density_nsi(ielem)%a,pgaus)
     end do
  end if
  !
  ! Boundary mass
  !
  if( kfl_waexl_ker /= 0 .or. kfl_twola_ker /= 0 .or. kfl_noslw_ker /= 0 ) then
     call memory_alloca(mem_modul(1:2,modul),'MASSB_NSI','nsi_solmem',massb_nsi,npoin)
  end if
  !
  ! Immersed boundary coupling for deformable bodies 
  !
  if( kfl_dimbou .or. kfl_immer ) then
     call memory_alloca(mem_modul(1:2,modul),'FSIFO_NSI','nsi_solmem',fsifo_nsi,ndime,npoin)
  end if
  !
  ! Velocity boundary condition used in Embedded Finite Element Coupling Technique (EFECT)
  !
  if( kfl_efect ) then
     call memory_alloca(mem_modul(1:2,modul),'VEFIX',    'nsi_solmem',vefix,ndime,npoin)
  end if

  !----------------------------------------------------------------------
  !
  ! Postprocess variables
  !
  !----------------------------------------------------------------------
  !
  ! RESCH_NSI: Schur complement residual
  !
  if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='SCHUR') ) then
     call memory_alloca(mem_modul(1:2,modul),'RESCH_NSI','nsi_solmem',resch_nsi,npoin)
  end if
  !
  ! REMOM_NSI: Schur complement residual
  !
  if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='MOMEN') ) then
     call memory_alloca(mem_modul(1:2,modul),'REMOM_NSI','nsi_solmem',remom_nsi,ndime,npoin)
  end if
  !
  ! TURMU_NSI: LES turbulent viscosity (postprocess)
  !
  if(       output_postprocess_check_variable_postprocess(VARIABLE_NAME='TURMU') &
       .or. output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVMUT') &
       .or. output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVSTX') &
       .or. output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVSTY') &
       .or. output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVSTZ') ) then
     call memory_alloca(mem_modul(1:2,modul),'TURMU_NSI','nsi_solmem',turmu_nsi,nelem)
     do ielem = 1,nelem
        pelty = abs(ltype(ielem))
        pgaus = ngaus(pelty)
        call memory_alloca(mem_modul(1:2,modul),'TURMU_NSI % A','nsi_solmem',turmu_nsi(ielem)%a,1_ip,pgaus,1_ip)
     end do
  end if
  !
  ! Forces due to porous media
  !
  if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='PORFO') )  then
     call memory_alloca(mem_modul(1:2,modul),'BUPOR_NSI','nsi_solmem',bupor_nsi,ndime,npoin)
  end if
  !
  ! Matrix copy for internal force
  !
  call nsi_algebraic_forces_allocate()
  
end subroutine nsi_solmem
