!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_doiter()
  !-----------------------------------------------------------------------
  !****f* Magnet/mag_doiter
  ! NAME
  !    mag_doiter
  ! DESCRIPTION
  !    This routine controls the non-linear loop of the H-formulation of 
  !    Maxwell's equations
  ! USES
  !    mag_solite
  ! USED BY
  !    Magnet
  !-----------------------------------------------------------------------

  use def_master, only: INOTSLAVE, timef !, kfl_reset
  use def_magnet
  use def_kermod, only: kfl_reset

  implicit none

  integer(ip) :: iconstr
  !
  ! Start non-linear loop
  !
  kfl_goiter_mag = .true.
  !
  ! Initialize non-linear counter
  !
  nonlCount_mag = 0_ip
  !
  do while (kfl_goiter_mag .and. nonlCount_mag < nlite_mag)
    !
    ! Assemble and solve system of equations
    !
    call mag_solite()
    !
    ! Increase non-linear counter
    !
    nonlCount_mag = nonlCount_mag + 1_ip
    !
    ! Give some hints to the user about how is everything going...
    ! Originally using solve_sol(1)
    !
    if (INOTSLAVE) then

      !write(*, 38) nonlCount_mag, lsiter0_mag, rmult_mag, gloRes_mag, relRes_mag, solVar_mag, Hnorm_mag
      !38 format ('  ', i3, '. ls Iter. = ', i5, '  Curr. res. = ', e11.4, '  gloRes = ', e11.4, '  relRes = ', e11.4, '  solVar = ', e11.4, '  Hnorm = ', e11.4)

      ! WIP
      write(*, 39) nonlCount_mag, lsiter0_mag, lsresi1_mag, gloRes_mag, relRes_mag, solVar_mag, Hnorm_mag
      39 format ('  ', i3, '. ls Iter. = ', i5, '   ls Res. = ', e11.4, '  gloRes = ', e11.4, '  relRes = ', e11.4, '  solVar = ', e11.4, '  Hnorm = ', e11.4)
      if (constr_total > 0) write(*, 40)
      do iconstr = 1, constr_total
        write(*, 41) iconstr, constr_mag(iconstr) % rmult, constr_mag(iconstr) % lsite, constr_mag(iconstr) % valuc
      end do
      40 format ('          CONSTRAINT        ABS. RESIDUAL        L.S. ITERATIONS        CONSTR. VALUE')
      41 format ('            ', i3, '              ', e11.4, '             ', i5, '              ', e11.4)
      ! WIP

    end if
    !
    ! Check wether we need to continue iterating: Max. iterations reached? Residual under tolerance?
    !
    if ( absRes_mag <= (0.01_rp * nltol_mag * absRes0_mag) .and. dHnorm_mag <= (nltol_mag * Hnorm_mag) ) then
      kfl_goiter_mag = .false.

      do iconstr = 1, constr_total
        if (abs(constr_mag(iconstr) % rmult) > (1.0e-12_rp + 1.0e-10_rp * abs(constr_mag(iconstr) % valuc) )) kfl_goiter_mag = .true.
      end do

    end if
    !
  end do
  !
  ! Convergence was not reached: repeat time step 
  !
  if (kfl_goiter_mag) then
    !
    ! Time step did not converge:
    !
    ! 1. Reduce time step
    !
    dtp_mag = dt_mag
    dt_mag = min( max( dtmin_mag, min( dtmax_mag, dt_mag * 0.5_rp ) ), timef - myclocktime_mag )
    !
    ! 2. Restore previous field and Lagrange multipliers
    !
    He_mag = Hp_mag
    !
    do iconstr = 1_ip, constr_total
      constr_mag(iconstr) % multc = constr_mag(iconstr) % multp
    end do
    !
    ! 3. Tell Alya to repeat the time step
    !
    kfl_reset = 1_ip
    kfl_reset_mag = .true.
    !
    ! 4. Tell the user that the time step is going to be repeated
    !
    if (INOTSLAVE) then
      write(*, 38) dt_mag
      38 format ('  Convergence was not reached. Reduced time step will be used, dt = ', e11.4)
    end if
    !
    ! You need to call here mag_begste since Alya does not call the module begste when kfl_reset = 1
    !
    call mag_timste()
    call mag_begste()
    !
  else
    !
    ! Time step converged. Tell Alya to continue
    !
    kfl_reset = -1_ip
    kfl_reset_mag = .false.
    !
  end if
  !
end subroutine mag_doiter
