!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_mag_bdfode

  use def_kintyp, only: ip, rp
  use mod_mag_linalg

  implicit none

  type stcbdf
    integer(ip) :: s = 1_ip, sc
    real(rp) :: ai(7), ae(7)
    real(rp) :: erri, erre
    real(rp) :: ai0, ae0
    ! dt(1) = t^{n+1} - ^t{n}
    ! dt(2) = t^{n} - t^{n-1}
    ! dt(3) = t^{n-1} - t^{n-2}
    ! dt(4) = t^{n-2} - t^{n-3}
    ! dt(5) = t^{n-3} - t^{n-4}
    ! dt(6) = t^{n-4} - t^{n-5}
    real(rp) :: dt(6) 
    ! y1 = y^{n}
    ! y2 = y^{n-1}
    ! y3 = y^{n-2}
    ! y4 = y^{n-3}
    ! y5 = y^{n-4}
    ! y6 = y^{n-5}
    real(rp), pointer :: y1(:), y2(:), y3(:), y4(:), y5(:), y6(:)
    real(rp), pointer :: ypi(:), ype(:)
  end type stcbdf

contains

  !##############################################################
  subroutine bdfodeTestCoefficients()

    implicit none

    real(rp) :: ai(7), ae(7), dt(6)

    integer(ip) :: s

    dt = [0.1_rp, 0.2_rp, 0.16_rp, 0.03_rp, 0.15_rp, 0.22_rp]
    dt = 1.0_rp

    do s = 1_ip, 6_ip
      ai = bdfodeImplicitCoefficients(s, dt)
      write(*, '(A,I3,A,7F10.6,A)') "s=", s, "   a=[", ai(1), ai(2), ai(3), ai(4), ai(5), ai(6), ai(7), "]"
    end do

    do s = 1_ip, 6_ip
      ae = bdfodeExplicitCoefficients(s, dt)
      write(*, '(A,I3,A,7F10.6,A)') "s=", s, "   a=[", ae(1), ae(2), ae(3), ae(4), ae(5), ae(6), ae(7), "]"
    end do

    do s = 1_ip, 6_ip
      write(*, '(A,F10.6)') "s!=", factorial(s)
    end do

    do s = 1_ip, 6_ip
      ai = bdfodeImplicitCoefficients(s, dt)
      write(*, '(A,F10.6)') "C_{s+1}=", bdfodeImplicitExpansion(s + 1_ip, s, dt, ai)
    end do

    do s = 1_ip, 6_ip
      ae = bdfodeExplicitCoefficients(s, dt)
      write(*, '(A,F10.6)') "C_{s+1}=", bdfodeExplicitExpansion(s + 1_ip, s, dt, ae)
    end do

  end subroutine bdfodeTestCoefficients
  !##############################################################


  !##############################################################
  subroutine bdfodeUpdate(bdf, dt, y)

    implicit none

    type (stcbdf) :: bdf

    real(rp), intent(in) :: dt
    real(rp), intent(in) :: y(:)
    !
    ! Save previous time steps
    !
    if (bdf % sc >= 6_ip) bdf % dt(6) = bdf % dt(5)
    if (bdf % sc >= 5_ip) bdf % dt(5) = bdf % dt(4)
    if (bdf % sc >= 4_ip) bdf % dt(4) = bdf % dt(3)
    if (bdf % sc >= 3_ip) bdf % dt(3) = bdf % dt(2)
    if (bdf % sc >= 2_ip) bdf % dt(2) = bdf % dt(1)
    if (bdf % sc >= 1_ip) bdf % dt(1) = dt
    !
    ! BDF coefficients: implicit and explicit
    !
    bdf % ai = bdfodeImplicitCoefficients(bdf % sc, bdf % dt)
    bdf % ae = bdfodeExplicitCoefficients(bdf % sc, bdf % dt)
    !
    ! BDF error coefficients: implicit and explicit
    !
    bdf % erri = bdfodeImplicitExpansion(bdf % s + 1_ip, bdf % s, bdf % dt, bdf % ai)
    bdf % erre = bdfodeExplicitExpansion(bdf % s + 1_ip, bdf % s, bdf % dt, bdf % ae)
    !
    ! Contributions to ypi and ype
    !
    bdf % ypi = - bdf % ai(2) * y !(1:size(bdf % ypi))
    bdf % ype = - bdf % ae(2) * y !(1:size(bdf % ype))
    !
    if (bdf % sc >= 6_ip) then
      bdf % y6 = bdf % y5
      bdf % ypi = bdf % ypi - bdf % ai(7) * bdf % y6
      bdf % ype = bdf % ype - bdf % ae(7) * bdf % y6
    end if
    if (bdf % sc >= 5_ip) then
      bdf % y5 = bdf % y4
      bdf % ypi = bdf % ypi - bdf % ai(6) * bdf % y5
      bdf % ype = bdf % ype - bdf % ae(6) * bdf % y5
    end if
    if (bdf % sc >= 4_ip) then
      bdf % y4 = bdf % y3
      bdf % ypi = bdf % ypi - bdf % ai(5) * bdf % y4
      bdf % ype = bdf % ype - bdf % ae(5) * bdf % y4
    end if
    if (bdf % sc >= 3_ip) then
      bdf % y3 = bdf % y2
      bdf % ypi = bdf % ypi - bdf % ai(4) * bdf % y3
      bdf % ype = bdf % ype - bdf % ae(4) * bdf % y3
    end if
    if (bdf % sc >= 2_ip) then
      bdf % y2 = bdf % y1
      bdf % ypi = bdf % ypi - bdf % ai(3) * bdf % y2
      bdf % ype = bdf % ype - bdf % ae(3) * bdf % y2
    end if
    if (bdf % sc >= 1_ip) then
      bdf % y1 = y !(1:size(bdf % y1))
    end if
    !
    bdf % ai0 = bdf % ai(1)
    bdf % ae0 = bdf % ae(1)
    !
  end subroutine bdfodeUpdate
  !##############################################################


  !##############################################################
  function bdfodeExplicitCoefficients(s, dt) result(a)
    !
    ! Returns the coefficients for BDF explicit formulas
    ! Valid coefficients for variable time steps are implemented for s = 1, 2
    !
    implicit none

    integer(ip), intent(in) ::                  s

    real(rp), optional, intent(in) ::           dt(6)

    real(rp) ::    &
      a(7),    &
      As(s, s),    &
      bs(s),    &
      xs(s),    &
      dtaux

    integer(ip) ::    &
      ii,    &
      jj

    a = 0.0_rp

    if (s == 1_ip) then
      a(1) = 1.0_rp
      a(2) = -1.0_rp
    elseif (s > 6_ip) then
      call runend("bdfodeExplicitCoefficients: invalid input s")
    else
      dtaux = 0.0_rp
      As(:, 1) = 1.0_rp
      bs = 0.0_rp
      bs(1) = 1.0_rp
      do ii = 2_ip, s
        dtaux = dtaux + dt(ii) / dt(1)
        As(1, ii) = - dtaux
        do jj = 2_ip, s
          As(jj, ii) = - As(jj-1, ii) * dtaux
        end do
      end do

      call lusolver(xs, As, bs)

      a(1) = xs(1)
      a(3:s+1) = xs(2:s)
      a(2) = -sum(a(1:s+1))
      
    end if

  end function bdfodeExplicitCoefficients
  !##############################################################


  !##############################################################
  function bdfodeImplicitCoefficients(s, dt) result(a)
    !
    ! Returns the coefficients for BDF implicit formulas 1 <= s <= 6
    ! Valid coefficients for variable time steps are implemented for s = 1, 2, 3, 4, 5, 6
    !
    implicit none

    integer(ip), intent(in) ::                  s

    real(rp), optional, intent(in) ::           dt(6)

    real(rp) ::  &
         a(7),    &
         As(s, s),    &
         bs(s),    &
         xs(s),    &
         dtaux          

    integer(ip) :: ii, jj

    a = 0.0_rp

    if (s == 1_ip) then
       a(1) = 1.0_rp
       a(2) = -1.0_rp
    elseif (s > 6_ip) then
       call runend("bdfodeImplicitCoefficients: invalid input s")
    else
       if (present(dt)) then
          dtaux = 0.0_rp
          As(:, 1) = 1.0_rp
          bs = 0.0_rp
          bs(1) = 1.0_rp
          do ii = 2_ip, s
             dtaux = dtaux + dt(ii) / dt(1)
             bs(ii) = bs(ii-1) + 1.0_rp
             As(1, ii) = -dtaux
             do jj = 2_ip, s
                As(jj, ii) = - As(jj-1, ii) * dtaux 
             end do
          end do

          call lusolver(xs, As, bs)

          a(1) = xs(1)
          a(3:s+1) = xs(2:s)
          a(2) = -sum(a(1:s+1))
       else

          select case (s)

          case (2_ip)
             a(1) = 1.5_rp
             a(2) = -2.0_rp
             a(3) = 0.5_rp

          case (3_ip)
             a(1) = 11.0_rp / 6.0_rp
             a(2) = -3.0_rp
             a(3) = 1.5_rp
             a(4) = -1.0_rp / 3.0_rp

          case (4_ip)
             a(1) = 25.0_rp / 12.0_rp
             a(2) = -4.0_rp
             a(3) = 3.0_rp
             a(4) = -4.0_rp / 3.0_rp
             a(5) = 0.25_rp

          case (5_ip)
             a(1) = 137.0_rp / 60.0_rp
             a(2) = -5.0_rp
             a(3) = 5.0_rp
             a(4) = -10.0_rp / 3.0_rp
             a(5) = 1.25_rp
             a(6) = -0.25_rp

          case (6_ip)
             a(1) = 147.0_rp / 60.0_rp
             a(2) = -6.0_rp
             a(3) = 7.5_rp
             a(4) = -20.0_rp / 3.0_rp
             a(5) = 3.75_rp
             a(6) = -1.20_rp
             a(7) = 1.0_rp / 6.0_rp

          end select

       end if

    end if

  end function bdfodeImplicitCoefficients
  !##############################################################


  !##############################################################
  function bdfodeExplicitExpansion(k, s, dt, a) result(ck)
    !
    ! Returns the error constant of the implicit BDF formula of order s + 1
    ! error = O(h**(s + 1))
    !
    implicit none

    integer(ip), intent(in) :: k, s

    real(rp), intent(in) :: dt(6), a(7)

    real(rp) :: ck, dtaux

    integer(ip) :: ii

    ck = 0.0_rp
    dtaux = 0.0_rp

    do ii = 2_ip, s
      dtaux = dtaux + dt(ii) / dt(1)
      ck = ck + dtaux**k * a(ii+1)
    end do

    if (k == 1_ip) then
      ck = (ck * (-1.0_rp)**k + a(1) -1.0_rp) / factorial(k)
    else
      ck = (ck * (-1.0_rp)**k + a(1)) / factorial(k)
    end if

    ck = ck / a(1)

  end function bdfodeExplicitExpansion
  !##############################################################


  !##############################################################
  function bdfodeImplicitExpansion(k, s, dt, a) result(ck)
    !
    ! Returns the error constant of the implicit BDF formula of order s + 1
    ! error = O(h**(s + 1))
    !
    implicit none

    integer(ip), intent(in) :: k, s

    real(rp), intent(in) :: dt(6), a(7)

    real(rp) :: ck, dtaux

    integer(ip) :: ii

    ck = 0.0_rp
    dtaux = 0.0_rp

    do ii = 2_ip, s
      dtaux = dtaux + dt(ii) / dt(1)
      ck = ck + dtaux**k * a(ii+1)
    end do

    ck = (ck * (-1.0_rp)**k + a(1) - real(k,rp)) / factorial(k)

    ck = ck / a(1)

  end function bdfodeImplicitExpansion
  !##############################################################


  !##############################################################
  function factorial(n) result(f)

    implicit none

    integer(ip), intent(in) :: n

    real(rp) :: f

    integer(ip) :: ii

    f = 1.0_rp

    do ii = 1_ip, n
      f = f * real(ii,rp)
    end do

  end function factorial
  !##############################################################

end module mod_mag_bdfode
