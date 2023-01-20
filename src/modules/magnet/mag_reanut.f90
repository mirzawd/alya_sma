!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_reanut()

  use def_inpout
  use def_master
  use def_magnet
  use mod_ecoute, only: ecoute

  implicit none

  if (INOTSLAVE) then
    ! 
    ! Initializations (defaults)
    !
    dtmax_mag = 1.0e-5_rp
    dtmin_mag = 5.0e-9_rp
    theta_mag = 1.0_rp
    nltol_mag = 1.0e-4_rp
    reltol_mag = nltol_mag
    abstol_mag = 1.0e-12_rp
    nlite_mag = 15_ip
    nlide_mag = 5_ip
    gslin_mag = 4_ip
    gstri_mag = 4_ip
    gsqua_mag = 4_ip
    gstet_mag = 4_ip
    gshex_mag = 3_ip
    !
    ! Reach the NUMERICAL_TREATMENT section
    !
    call ecoute('mag_reanut')
    do while (words(1) /= 'NUMER')
      call ecoute('mag_reanut')
    end do
    !
    ! Start reading data
    !
    do while (words(1) /= 'ENDNU')
      call ecoute('mag_reanut')

      if (words(1) == 'ALGEB') then
        !
        ! Algebraic solver
        !
        solve_sol => solve(1:)
        call reasol(1_ip)

      else if (words(1) == 'DTMIN') then

        dtmin_mag = param(1)

      else if (words(1) == 'DTMAX') then

        dtmax_mag = param(1)

      else if (words(1) == 'THETA') then

        theta_mag = param(1)

      else if (words(1) == 'NLTOL') then

        nltol_mag = param(1)
        reltol_mag = nltol_mag

      else if (words(1) == 'RETOL') then

        reltol_mag = param(1)

      else if (words(1) == 'ABTOL') then

        abstol_mag = param(1)

      else if (words(1) == 'NLITE') then

        nlite_mag = int(param(1), ip)

      else if (words(1) == 'NLIDE') then

        nlide_mag = int(param(1), ip)

      else if (words(1) == 'MADOF') then

        !maxdof_mag = int(param(1), ip)

      else if (words(1) == 'MAGAU') then

        !maxgau_mag = int(param(1), ip)

      else if (words(1) == 'GSLIN') then

        gslin_mag = int(param(1), ip)

      else if (words(1) == 'GSTRI') then

        gstri_mag = int(param(1), ip)

      else if (words(1) == 'GSQUA') then

        gsqua_mag = int(param(1), ip)

      else if (words(1) == 'GSTET') then

        gstet_mag = int(param(1), ip)

      else if (words(1) == 'GSHEX') then

        gshex_mag = int(param(1), ip)

      else if (words(1) == 'BDFOR') then

        bdfode_mag % s = int(param(1), ip)

      else if (words(1) == 'STRUC') then

        struct_mag = option('STRUC')

      else if (words(1) == 'GAUSS') then

        postev_mag = option('GAUSS')

      else if (words(1) == 'CENTR') then

        postce_mag = option('CENTR')

      end if

    end do

  end if

end subroutine mag_reanut
