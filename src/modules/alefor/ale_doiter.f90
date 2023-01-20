!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_doiter
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_doiter
  ! NAME 
  !    ale_doiter
  ! DESCRIPTION
  !    Check if ALEFOR should be solved.
  !    It is not solved if:
  !    - All d.o.f. are prescribed
  !    - All prescribed nodes are prescribed to zero
  ! USES
  !    ale_begite
  !    ale_solite
  !    ale_endite
  ! USED BY
  !    Aelfor
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_solver
  use def_alefor
  use mod_timings, only : timings_ini
  use mod_timings, only : timings_end
  implicit none
  integer(ip) :: kfl_fmale_ale,kfl_bvess_ale ! ,kfl_solve_ale
    
  kfl_fmale_ale = 0_ip ! kfl_fmale_ale initialization
  kfl_solve_ale = 1_ip ! kfl_solve_ale initialization, alefor will always solve when called
  kfl_bvess_ale = 1_ip ! kfl_bvess_ale_initialization, always use the boundary conditions in bvess_ale

  !------------------------------------------------------
  !
  ! Rigid body motion
  !
  !------------------------------------------------------
  
  if( kfl_rigid_ale /= 0_ip ) call ale_solrbo()
  
  !-------------------------------------------------------------------
  !
  ! Solve mesh deformation and smoothing
  !
  !-------------------------------------------------------------------

  call timings_ini()
  call ale_begite()
  call timings_end(ITASK_BEGITE)
  call ale_smodef()
  
  !-------------------------------------------------------------------
  !
  ! Recompute some domain variables and output domain mesh
  !
  !-------------------------------------------------------------------
  
  if( kfl_fmale_ale == 0_ip ) kfl_domar = 1_ip

  !-------------------------------------------------------------------
  !
  ! End of iterations
  !
  !-------------------------------------------------------------------

  call timings_ini()
  call ale_endite()
  call timings_end(ITASK_ENDITE)

end subroutine ale_doiter

 
