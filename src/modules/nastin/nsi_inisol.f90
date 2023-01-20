!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_inisol(itask)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_inisol
  ! NAME 
  !    nsi_inisol
  ! DESCRIPTION
  !    This routine loads the solver data for the incomp. NS equations.
  !    In general, it may change from time step to time step or even
  !    from iteration to iteration.
  ! USED BY
  !    nsi_begite
  !    nsi_elmcor
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  use def_solver
  implicit none
  integer(ip), intent(in) :: itask

  if( dtinv /= 0.0_rp ) solve(1) % xdiag = 1.0_rp/dtinv
  ! 
  ! Which equation(s) is(are) to be solved
  !
  if( itask == 1 .or. itask == 2 .or. itask == 4 ) then
     ivari_nsi  = ivari_nsi_mom                           ! Momentum (+Continuity) equation  
  else if( itask == 3 ) then
     ivari_nsi  = ivari_nsi_cont                          ! Continuity
  end if

  if( itask == 4 ) then
     !
     ! Check boundary conditions
     !
     solve_sol => solve(3:)

  else if( itask == 5 ) then          
     !                     
     ! Diagonal solver for mass correction (nsi_elmcor)
     !
     solve_sol => solve(4:)
     
  else
     !
     ! Current variable
     !
     solve_sol => solve(ivari_nsi:)
    
  end if

end subroutine nsi_inisol
