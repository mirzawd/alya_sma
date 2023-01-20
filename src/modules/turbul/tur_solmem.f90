!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_solmem()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_memall
  ! NAME 
  !    tur_memall
  ! DESCRIPTION
  !    This routine allocates memory for the arrays needed to solve the
  !    turbulence equations      
  ! USES
  ! USED BY
  !    tur_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_solver
  use def_turbul
  use mod_memory
  use mod_tur_arrays, only : tur_arrays
  use mod_ADR, only : FULL_OSS
  use mod_ADR, only : A_OSS  
  use mod_ADR, only : AR_OSS 
  use mod_ADR, only : BUBBLE
  use mod_ADR, only : ADR_initialize_type
  use mod_ADR, only : ADR_check_and_compute_data
  use mod_ADR, only : ADR_allocate_projections_bubble_sgs 
  use mod_output_postprocess, only : output_postprocess_check_variable_postprocess
  implicit none
  integer(ip) :: iturb

  !
  ! Solver memory
  !
  solve_sol => solve(1:1)
  call soldef(4_ip)
  solve_sol => solve(3:3)
  call soldef(4_ip)
  if(kfl_algor_tur==1) then
     solve_sol => solve(2:nturb_tur)
     call soldef(4_ip)
  end if
  !
  ! ADR type
  ! 
  do iturb = 1,nturb_tur
     call ADR_initialize_type(ADR_tur(iturb))
     ADR_tur(iturb) % kfl_time_integration   =  kfl_timei_tur
     ADR_tur(iturb) % kfl_time_step_strategy =  kfl_timco
     ADR_tur(iturb) % kfl_stabilization      =  kfl_ortho_tur
     ADR_tur(iturb) % kfl_shock              =  kfl_shock_tur
     ADR_tur(iturb) % kfl_time_lumped        =  0
     ADR_tur(iturb) % kfl_tau_strategy       =  kfl_taust_tur
     ADR_tur(iturb) % kfl_laplacian          =  0 
     ADR_tur(iturb) % kfl_nonlinear_sgs      =  kfl_sgsno_tur 
     ADR_tur(iturb) % kfl_time_sgs           =  kfl_sgsti_tur
     ADR_tur(iturb) % kfl_time_bubble        =  kfl_tibub_tur
     ADR_tur(iturb) % kfl_time_scheme        =  kfl_tisch_tur
     ADR_tur(iturb) % kfl_time_order         =  kfl_sgsac_tur
     ADR_tur(iturb) % kfl_manufactured       =  kfl_exacs_tur
     ADR_tur(iturb) % kfl_length             =  kfl_ellen_tur
     ADR_tur(iturb) % number_euler_steps     =  neule_tur

     ADR_tur(iturb) % lun_output4            =  int(momod(modul) % lun_outpu,4)
     ADR_tur(iturb) % bemol                  =  bemol_tur
     ADR_tur(iturb) % tau_parameters(1:3)    =  staco_tur(1:3)
     ADR_tur(iturb) % shock                  =  shock_tur  
     call ADR_check_and_compute_data(ADR_tur(iturb))
     call ADR_allocate_projections_bubble_sgs(ADR_tur(iturb))
  end do

end subroutine tur_solmem
