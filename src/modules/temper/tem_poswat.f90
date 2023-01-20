!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Temper
!> @{
!> @file    tem_poswat.f90
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
!>          - qair   = -qL - k_memb/h_memb * (Tair-Tbase)
!>          - qL     = w/Area * Dh 
!>          - w/Area = rho_wv k_wv * gradc_wv.n 
!>          - rho_wv = a*T^2 + b*T + c
!>            with a=2.0e-5,b=0.0003,c=0.0025
!>          - qair   = -(a*T^2 + b*T + c) k_wv * gradc_wv.n *Dh
!>                     - k_memb/h_memb * (T-Tbase)
!>                   = [ -(a*T+b)*k_wv * gradc_wv.n *Dh - k_memb/h_memb ] T
!>                     - c * k_wv * gradc_wv.n *Dh + k_memb/h_memb * Tbase
!>          Units
!>          J   = W s = [kg m^2 / s^2]
!>          W   = [kg m^2 / s]
!>          k   = thermal conductivity  [W/mK] = [kg m / K s^3]
!>          kap = thermal diffusivity   [m^2/s]
!> @} 
!-----------------------------------------------------------------------
subroutine tem_poswat()
  use def_kintyp
  use def_master
  use def_domain
  use def_temper
  implicit none
  real(rp)    :: T
  real(rp)    :: qrobi
  real(rp)    :: arobi
!  real(rp)    :: trobi
  integer(ip) :: idime,ibopo,ipoin
  real(rp)    :: hmemb,kmemb,kapwv,DH,Tbase
  real(rp)    :: a,b,c,d,gradc_dot_n,rhowv

  do ipoin = 1,npoin

     ibopo = lpoty(ipoin)

     if( ibopo > 0 ) then
        T     = tempe(ipoin,1) 
        a     = 2.0e-5_rp
        b     = 0.0003_rp
        c     = 0.0025_rp
        rhowv = a * T**2 + b * T + c                          ! Water vapor density:          [kg/m^3]
        kapwv = 2.338e-5_rp                                   ! Water vapor diffusivity       [m^2/s]
        DH    = 2.411e6_rp                                    ! Latent heat evaporation       [J/Kg]
        hmemb = 0.001_rp                                      ! Membrane width                [m]
        kmemb = 0.58_rp                                       ! Membrane conductivity (water) [W/mK]
        Tbase = 34.0_rp                                       ! Base temperature              [K]

        gradc_dot_n = 0.0_rp
        do idime = 1,ndime
           gradc_dot_n = gradc_dot_n + gradc_tem(idime,ipoin) * exnor(idime,1,ibopo)
        end do
        !
        ! Robin condition
        ! qair = k*grad(T).n = qrobi + arobi       * (T   -trobi)
        !                   = -qL    - kmemb/hmemb * (Tair-Tbase)
        !
        d     =  gradc_dot_n * kapwv * DH
        arobi = -(a*T+b)*d-kmemb/hmemb
        qrobi = -c*d + kmemb/hmemb*Tbase

        gesca(ipoin) = qrobi + arobi * T

     end if

  end do

end subroutine tem_poswat
