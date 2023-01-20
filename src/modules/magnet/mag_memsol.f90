!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_memsol()

  use def_kintyp, only: ip
  use def_master, only: solve
  use def_solver, only: solve_sol
  use def_magnet
  
  implicit none
  !
  ! The master needs to run the solver memory setup?
  !
  !-----------------------------------------------------------------------------
  ! Solver Memory Setup
  !-----------------------------------------------------------------------------
  !
  ! Enterate bien de que por que se hace esto...
  !
  solve_sol                => solve(1:1)
  solve_sol(1) % bvess     => bvess_mag(:,:,1)
  solve_sol(1) % kfl_fixno => kfl_fixno_mag
  !
  ! Module fails here...
  ! Something needs to be done regarding edge elements for matrices allocation
  ! Dimensions are saved in 'master/soldef.f90'
  ! Memory is allocated in 'master/memunk'
  !
  call soldef(4_ip)
  !
end subroutine mag_memsol
