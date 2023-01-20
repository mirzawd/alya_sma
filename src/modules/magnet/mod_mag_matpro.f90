!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_mag_matpro

  use def_kintyp, only: ip, rp, lg
  use def_magnet, only: rho_mag, Ec0_mag, nc0_mag, Jc0_mag, B0k_mag, mu0_mag, scalinOpt_mag, resistOpt_mag, &
                        ZER_MAG, resmat_mag, Jcrmat_mag, rhomat_mag, Jc0mat_mag, Ecrmat_mag, ncrmat_mag, &
                        Bc0mat_mag, Tc0mat_mag, intp1_mag, interpLinear1, murmat_mag

  implicit none

contains

  !##############################################################
  function mag_scalin(material, x, Hmag, Tmag, icomp) result(Jc)

    implicit none

    integer(ip), intent(in) ::    &
      material

    real(rp), intent(in) ::    &
      x(:),    &
      Hmag(:),    &
      Tmag

    integer(ip), intent(in), optional ::    &
      icomp

    real(rp) ::    &
      Jc,    &
      theta

    integer(ip) ::    &
      i

    if (present(icomp)) then
      i = icomp
    else
      i = 3_ip
    end if

    select case (Jcrmat_mag(i, material))
    case (0_ip)
      Jc = Jc0mat_mag(i, material)
    case (1_ip)
      Jc = Jc0mat_mag(i, material) * Bc0mat_mag(i, material) / (Bc0mat_mag(i, material) + mu0_mag * SQRT( dot_product(Hmag, Hmag) ))
    case (9_ip)
!      Jc = mag_fitted(x) !interpLinear1(intp1_mag(1), 0.2_rp)
      Jc = (189.53456_rp / 0.0000012_rp) * (1.0_rp + ((sqrt(dot_product(Hmag,Hmag)) * mu0_mag) / (0.05184_rp)))**(-0.50016_rp) * (1.0_rp - ((sqrt(dot_product(Hmag,Hmag)) * mu0_mag) / (7.22788_rp)))**2.0_rp

      if (Hmag(2) == 0.0_rp) then
        if (Hmag(1) < 0.0_rp)then
          theta = (3.0_rp * 3.141592_rp) / 2.0_rp
        else
          theta = 3.141592_rp / 2.0_rp
        end if
      else
        theta = atan(Hmag(1) / Hmag(2))
        if (Hmag(1) > 0.0_rp)then
          if (Hmag(2) > 0.0_rp)then
            theta = theta
          elseif (Hmag(2) < 0.0_rp)then
            theta = 3.141592_rp + theta
          end if
        elseif (Hmag(1) < 0.0_rp)then
          if (Hmag(2) < 0.0_rp)then
            theta = 3.141592_rp + theta
          elseif (Hmag(2) > 0.0_rp)then
            theta = 2.0_rp * 3.141592_rp + theta
          end if
        end if
      end if
     
      if (theta > 3.141592_rp)then
        theta = theta - 3.141592_rp
      end if

      Jc = Jc * interpLinear1(intp1_mag(1), theta * 180.0_rp / 3.141592_rp) / interpLinear1(intp1_mag(1), 30.00_rp)

    case (-1_ip) 
      Jc = 0.0_rp
    case default
      Jc = 0.0_rp
      call runend('mag_scalin: unknown scaling law')
    end select
      
  end function mag_scalin
  !##############################################################


  !##############################################################
  function mag_resist(material, x, Jmag, Hmag, Tmag, axsym, icomp) result(rho)

    implicit none

    integer(ip), intent(in) ::    &
         material

    real(rp), intent(in) ::    &
         x(:),    &
         Jmag(:),    &
         Hmag(:),    &
         Tmag

    logical(lg), intent(in) ::    &
         axsym

    integer(ip), intent(in), optional ::    &
         icomp

    real(rp) ::    &
         Jc1,    &
         rho,    &
         Jn

    integer(ip) ::    &
         i

    if (present(icomp)) then
      i = icomp
    else
      i = 3_ip
    end if

    select case (resmat_mag(i, material))
    case (0_ip)
      rho = rhomat_mag(i, material)
    case (1_ip)
      Jc1 = mag_scalin(material, x, Hmag, Tmag, i)
      Jn = SQRT( dot_product(Jmag, Jmag) ) / Jc1
      !
      ! The following if-else statement is needed to avoid overflow
      ! Otherwise the code will not work properly with debug flags
      !
      if (Jn <= 100.0_rp * ZER_MAG**(1.0_rp / ncrmat_mag(i, material))) then
        rho = 0.0_rp
      else
        rho = Ecrmat_mag(i, material) / Jc1 * Jn**(ncrmat_mag(i, material) - 1.0_rp)
      end if
    case (-1_ip)
      rho = 0.0_rp
    case default
      rho = 0.0_rp
      call runend('mag_resist: unknown resistivity function')
    end select
    
    if (axsym) then
      rho = rho * x(1)
    end if

  end function mag_resist
  !##############################################################


  !##############################################################
  function mag_resdif(material, x, Jmag, Hmag, Tmag, axsym, icomp) result(rhodif)

    implicit none

    integer(ip), intent(in) ::    &
      material

    integer(ip), intent(in), optional ::    &
      icomp

    real(rp), intent(in) ::    &
      x(:),    &
      Jmag(:),    &
      Hmag(:),    &
      Tmag

    logical(lg), intent(in) ::    &
      axsym

    real(rp) :: &
      rhodif,    &
      Jc1, &
      Jn

    integer(ip) ::    &
      i

    if (present(icomp)) then
      i = icomp
    else
      i = 3_ip
    end if

    select case (resmat_mag(i, material))
    case (0_ip)
      rhodif = 0.0_rp
    case (1_ip) 
      Jc1 = mag_scalin(material, x, Hmag, Tmag, i)
      Jn = SQRT( dot_product(Jmag, Jmag) ) / Jc1
      !
      ! The following if-else statement is needed to avoid overflow
      ! Otherwise the code will not work properly with debug flags
      !
      if (Jn <= 100.0_rp * ZER_MAG**(1.0_rp / ncrmat_mag(i, material))) then
        rhodif = 0.0_rp
      else
        rhodif = (ncrmat_mag(i, material) - 1.0_rp) * Ecrmat_mag(i, material) / &
                 (Jc1**2.0_rp) * Jn**(ncrmat_mag(i, material) - 2.0_rp)
      end if
    case (-1_ip)
      rhodif = 0.0_rp
    case default
      rhodif = 0.0_rp
      call runend('mag_resdif: unknown resistivity function')
    end select

    if (axsym) then
      rhodif = rhodif * x(1)
    end if

  end function mag_resdif
  !##############################################################


  !##############################################################
  function mag_permea(material, x, Hmag, axsym, icomp) result(mu)

    implicit none

    integer(ip), intent(in) ::    &
      material

    integer(ip), intent(in) ::    &
      icomp

    real(rp), intent(in) ::    &
      x(:),    &
      Hmag(:)

    logical(lg), intent(in) ::    &
      axsym

    real(rp) ::    &
      mu

    mu = murmat_mag(icomp, material)

    if (axsym) then
      mu = mu * x(1)
    end if

  end function mag_permea
  !##############################################################


  !##############################################################
  function mag_fitted(x) result(Jc)

    implicit none

    real(rp), intent(in) ::    &
      x(:)

    real(rp) ::    &
      Jc,    &
      a(8),    &
      b(8),    &
      c(8),    &
      xn

    integer(ip) ::    &
      i

    a = [0.1822_rp, 0.5522_rp, 0.9064_rp, 0.2437_rp, 0.02993_rp, 0.4885_rp, 0.3275_rp, 0.3233_rp]
    b = [0.7413_rp, -0.7203_rp, -0.2007_rp, -0.9981_rp, -0.03299_rp, 0.3284_rp, 0.5306_rp, 0.6692_rp]
    c = [0.07313_rp, 0.2969_rp, 0.5274_rp, 0.1505_rp, 0.8686_rp, 0.2943_rp, 0.1768_rp, 0.08369_rp]

    xn = (x(1) - 0.0011009_rp) / 0.006245_rp

    Jc = 0.0_rp

    do i = 1_ip, 8_ip
      Jc = Jc + a(i)*exp(-((xn - b(i))/c(i))**2)
    end do

    Jc = Jc * 7.5785d9 / 50.0_rp

  end function mag_fitted
  !##############################################################

end module mod_mag_matpro
