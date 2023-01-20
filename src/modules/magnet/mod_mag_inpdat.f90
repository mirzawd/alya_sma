!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_mag_inpdat

  use def_kintyp, only: ip, rp
  use def_parame, only: pi
  use def_magnet, only: mu0_mag
  use def_master, only: namda
  
  implicit none

contains

  !##############################################################
  function mag_inifie(x) result(H)

    implicit none

    real(rp), intent(in) :: x(:)

    real(rp) :: H(size(x))

    H = 0.0_rp

    if (namda(1:3) == 'TS-' .AND. size(x) == 2_ip) then
      if (namda(1:8) == 'TS-TAPE-') then
        H(2) = 0.5_rp / mu0_mag
      end if
    elseif (namda(1:3) == 'TS-' .AND. size(x) == 3_ip) then
      if (namda(1:8) == 'TS-TAPE-') then
        H(2) = 0.5_rp / mu0_mag
      end if
    else
    end if

  end function mag_inifie
  !##############################################################


  !##############################################################
  function mag_dirfie(t, x) result(field)

    implicit none

    real(rp), intent(in) :: &
      t

    real(rp), intent(in) :: &
      x(:)

    real(rp) :: &
      field(size(x))

    real(rp) :: &
      r,    &
      cost,    &
      sint

    if (namda(1:3) == 'TS-' .AND. size(x) == 2_ip) then

      if (namda(1:8) == 'TS-WIRE-') then

        r = sqrt( dot_product(x, x) )

        cost = x(1) / r
        sint = x(2) / r

        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        ! WIRE
        ! Rw = 1e-3 
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        field = [ -sint, cost ] * mag_curren(t) / (2_rp * pi * 0.001_rp)
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      elseif (namda(1:8) == 'TS-SLAB-') then
      
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        ! SLAB
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if (x(1) <= 0.d0) then
          field = [ 0.d0, -1.d0 ]
        elseif (x(1) >= 1.0e-3_rp) then
          field = [ 0.d0, 1.d0 ]
        elseif (x(2) <= 0.d0) then
          field = [ 0.d0, 0.d0 ]
        elseif (x(2) >= 1.0e-4_rp) then
          field = [ 0.d0, 0.d0 ]
        else
          call runend("mag_dirfie: wrong boundary conditions")
        end if

        field = field * 5.0e4_rp * sin(100_rp * pi * t)
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      elseif (namda(1:8) == 'TS-AIRC-') then

        r = sqrt( dot_product(x, x) )

        cost = x(1) / r
        sint = x(2) / r

        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        ! WIRE_AIR
        ! Close Boundary: Rw = 1e-3, Ra = 1.2e-3
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        field = [ -sint, cost ] * mag_curren(t) / (2_rp * pi * 0.0012_rp)
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      elseif (namda(1:8) == 'TS-AIRF-') then

        r = sqrt( dot_product(x, x) )

        cost = x(1) / r
        sint = x(2) / r
      
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        ! WIRE_AIR
        ! Far Boundary: Rw = 1e-3, Ra = 5e-3
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        field = [ -sint, cost ] * mag_curren(t) / (2_rp * pi * 0.005_rp)
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      else if (namda(1:8) == 'TS-SELF-') then

        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        ! SELF FIELD
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        field = [ 0.0_rp, 0.0_rp ]
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      else if (namda(1:8) == 'TS-TAPE-') then

        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        ! TAPE2D EXTERNAL FIELD
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        field = 0.0_rp
        if (t <= 0.005_rp) then
          field = 0.5_rp / mu0_mag * cos(100.0_rp * pi * t) * [0.0_rp, 1.0_rp]
        end if
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      else if (namda(1:8) == 'TS-RING-') then

        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        ! CURRENT LOOP
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        field = 0.0_rp
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      else

        call runend("mag_dirfie: unknown case")

      end if

    elseif (namda(1:3) == 'TS-' .AND. size(x) == 3_ip) then

      if (namda(1:8) == 'TS-SLAB-') then

        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        ! BENCH3D
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        field = [-sqrt(3.0_rp) / 2.0_rp, 0.0_rp, -0.5_rp] * 0.2_rp / mu0_mag * sin(100.0_rp * pi * t)
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      elseif (namda(1:8) == 'TS-TAPE-') then

        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        ! TAPE3D
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        field = 0.0_rp
        if (t <= 0.005_rp) then
          field(2) = 0.5_rp / mu0_mag * cos(100.0_rp * pi * t)
        end if
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      end if

    else
 
      field = 0.0_rp

      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      ! REAL CASE
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      if (t <= 0.005_rp) then
        field = [0.0_rp, 0.5_rp] * 1.0_rp / mu0_mag * cos(100.0_rp * pi * t)
      end if
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    end if

  end function mag_dirfie
  !##############################################################


  !##############################################################
  function mag_curren(t) result(I)

    implicit none

    real(rp), intent(in) :: &
      t

    real(rp) :: &
      I,    &
      I0 = 400.d0,    &
      freq = 50.d0

    I = I0 * sin(2_rp * pi * freq * t)

  end function mag_curren
  !##############################################################


  !##############################################################
  function mag_constr(t, iconstr) result(constr)

    implicit none

    integer(ip), intent(in) :: &
      iconstr

    real(rp), intent(in) :: &
      t

    real(rp) :: &
      constr

    if (iconstr == 1) then
      constr = mag_curren(t)
    elseif (iconstr == 2) then
      constr = 0.0_rp !mag_curren(t)
    else
      stop "mag_constr: undefined constraint"
    end if

  end function mag_constr
  !##############################################################


  !##############################################################
  function mag_source(t, x) result(src)

    implicit none

    real(rp), intent(in) :: t, x(:)

    real(rp) :: src(size(x))

!    src = [0.0_rp, -1.0_rp / mu0_mag * 2.0_rp * pi * 50.0_rp * cos(2.0_rp * pi * 50.0_rp * t)] * 0.0_rp

    src = 0.0_rp

  end function mag_source
  !##############################################################

end module mod_mag_inpdat
