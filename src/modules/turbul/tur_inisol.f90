!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_inisol(itask)
  !-----------------------------------------------------------------------
  !   
  ! This routine loads the solver data for the turbulence equations.
  ! In general, it may change from time step to time step or even
  ! from iteration to iteration.
  !
  !-----------------------------------------------------------------------
  use def_domain
  use def_turbul
  use def_master
  use def_solver
  implicit none
  integer(ip), intent(in) :: itask

  if(itask==2) then
     !
     ! Wall distance
     !
     solve_sol => solve(5:)                 ! Solver type
  else
     !
     ! Turbulence variable
     !
     if(kfl_algor_tur==1) then
        solve_sol => solve(iunkn_tur:)      ! Solver type
     else
        solve_sol => solve(1:)
     end if
  end if

end subroutine tur_inisol
