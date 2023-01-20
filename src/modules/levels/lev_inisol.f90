!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_inisol(iprob)
!-----------------------------------------------------------------------
!   
! This routine loads the solver data for the temperature equation.
! In general, it may change from time step to time step or even
! from iteration to iteration.
!
!-----------------------------------------------------------------------
  use def_domain
  use def_master
  use def_levels
  use def_solver
  implicit none
  integer(ip), intent(in) :: iprob

  solve_sol => solve(iprob:)
  
  if( solve(iprob)%kfl_algso == 9 ) solve(iprob)%xdiag = 1.0_rp/dtred_lev

!  pense que iba a poder adaptar lo de tur_inivar   pero ahora creo que no  notar que ahí el Kfl_fixno_tur es el del amster   
!  en cambio mis icupt_lev son de cada subdominio , de ultima puedo lograr que ande  solo en serie pensar
!
  if(iprob==3_ip) then
     !
     ! Deflated Solver: groups composed by master  borrowed from tur_inivar
     !
     if(solve_sol(1)%kfl_algso==2) then    ! if Deflated Solver

        call runend('lev_inisol: DEFLATED for POISSON not ready')
        !
        ! Keep only wall condition as positive to enable the creation
        ! of the groups for the deflated solver
        !
        solve_sol(1)%limpo => icupt_lev
        call cregro()
     end if
  end if



end subroutine lev_inisol
