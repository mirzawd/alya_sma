!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_updcon(itask)

  use def_master
  use def_magnet
  use def_solver, only: solve_sol
  use mod_solver, only: solver_parallel_scalar_product
  use mod_communications, only: PAR_INTERFACE_EDGE_EXCHANGE
  use mod_mag_inpdat

  implicit none

  integer(ip), intent(in) :: itask

  integer(ip)             :: iconstr

  select case (itask)

  case (ITASK_BEGSTE)

    if (INOTMASTER .and. kfl_lagr_mag) call mag_asscon()
    !
    do iconstr = 1, constr_total
      !
      ! Save previous Lagrange multiplier
      !
      constr_mag(iconstr) % multp = constr_mag(iconstr) % multc
      !
      ! Constraint value: current, previous, difference
      !
      constr_mag(iconstr) % valuc = mag_constr(myclocktime_mag + dt_mag, iconstr)
      constr_mag(iconstr) % valup = mag_constr(myclocktime_mag, iconstr)
      constr_mag(iconstr) % valud = constr_mag(iconstr) % valuc - constr_mag(iconstr) % valup
      !
      ! Exchange values at ghost edges
      !
      call PAR_INTERFACE_EDGE_EXCHANGE(solve_sol(1) % ndofn, constr_mag(iconstr) % cnstr, 'SUM', 'IN MY CODE')
      !
      ! Compute squared norm
      !
      call solver_parallel_scalar_product(solve_sol(1), constr_mag(iconstr) % cnstr, constr_mag(iconstr) % cnstr, constr_mag(iconstr) % normc)
      !
      ! Update magnetic field variable with initial guess
      !
      He_mag = He_mag + constr_mag(iconstr) % valud / constr_mag(iconstr) % normc * constr_mag(iconstr) % cnstr
      !
    end do

  end select

end subroutine mag_updcon
