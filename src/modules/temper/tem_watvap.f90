!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Temper
!> @{
!> @file    tem_watvap.f90
!> @author  Guillaume Houzeaux
!> @brief   Water vapor model
!> @details Water vapor model
!>          @verbatim
!>          -----------------------------------
!>          => Air
!>          
!>          -----------------------------------
!>          ----------------------------------- <= Interface
!>          Membrane  (h_memb)
!>          -----------------------------------
!>          Base tissue
!>
!>          @endverbatim
!>          - qair          =  qL - k_memb/h_memb * (Tair-Tbase)  [W / m^2]
!>          - k_memb        = Membrane conductivity               [W / m K]
!>          - h_memb        = Membrane thickness                  [m]
!>          - qL            = Mass/Area/sec * Dh                  [W / m^2]
!>          - Mass/Area/sec = Mass flux of water vapor            [Kg / m^2s]
!>                          = rho_wv * kap_wv * gradc_wv.n
!>          - rho_wv        = water vapor density                 [Kg / m^3]
!>          - kap_wv        = water vapor diffusivity             [m^2 / s]
!>          - Dh            = Latent heat of evaporation          [J / Kg] 
!>          Units
!>          - J   = Joule                                         [W s] = [kg m^2 / s^2]
!>          - k   = thermal conductivity                          [W / m K]
!>          - W   = Watts                                         [kg m^2 / s]
!>          - kap = diffusivity                                   [m^2 / s]
!> @} 
!-----------------------------------------------------------------------
subroutine tem_watvap(&
     pnodb,ndime,baloc,bogrc,gbsha,T,qrobi,arobi,trobi)
  use def_kintyp, only : ip,rp
  implicit none
  integer(ip), intent(in)  :: pnodb
  integer(ip), intent(in)  :: ndime
  real(rp),    intent(in)  :: baloc(ndime) 
  real(rp),    intent(in)  :: bogrc(ndime,pnodb)
  real(rp),    intent(in)  :: gbsha(pnodb)
  real(rp),    intent(in)  :: T
  real(rp),    intent(out) :: qrobi
  real(rp),    intent(out) :: arobi
  real(rp),    intent(out) :: trobi
  integer(ip)              :: idime,inodb
  real(rp)                 :: rhowv,gpgrc(3)
  real(rp)                 :: hmemb,kmemb,kapwv,DH,Tbase
  real(rp)                 :: gradc_dot_n,qL
!  real(rp)                 :: a,b,c
  !
  ! Coefficients
  !
  !a     = 2.0e-5_rp
  !b     = 0.0003_rp
  !c     = 0.0025_rp
  !rhowv = a * T**2 + b * T + c                          ! Saturation fraction of water in air: [kg/m^3]
  rhowv = 0.804_rp                                      ! Water vapor density                  [kg/m^3]
  kapwv = 2.338e-5_rp                                   ! Water vapor diffusivity              [m^2/s]
  DH    = 2.411e6_rp                                    ! Latent heat evaporation              [J/Kg]
  hmemb = 0.001_rp                                      ! Membrane width                       [m]
  kmemb = 0.58_rp                                       ! Membrane conductivity (water)        [W/mK]
  Tbase = 34.0_rp                                       ! Base temperature                     [K]
  !
  ! Latent heat of evaporation qL
  ! qL = rho_wv * kap_wv * |grad(c).n| 
  !
  gpgrc = 0.0_rp
  do inodb = 1,pnodb
     do idime = 1,ndime
        gpgrc(idime) = gpgrc(idime) + gbsha(inodb) * bogrc(idime,inodb)
     end do
  end do
  gradc_dot_n = 0.0_rp
  do idime = 1,ndime
     gradc_dot_n = gradc_dot_n + gpgrc(idime) * baloc(idime)
  end do
  qL =  rhowv * kapwv * gradc_dot_n * Dh
  !
  ! Robin condition
  ! qair = k*grad(T).n = qrobi + arobi       * (T - trobi)
  !                   =  qL    - kmemb/hmemb * (T - Tbase)
  !
  qrobi =  qL
  arobi = -kmemb / hmemb
  trobi =  Tbase

end subroutine tem_watvap
