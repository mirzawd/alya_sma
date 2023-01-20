!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_nsi_low_Mach

  use def_kintyp_basic,  only: ip
  use def_domain,        only: npoin
  use def_nastin,        only: drhodt_nsi
  use def_nastin,        only: mass_rho_nsi
  use def_nastin,        only: dtinv_nsi
  implicit none
  private

  public :: nsi_low_Mach_drho_dt
  public :: nsi_low_Mach_upd_mass_rho

contains 
  !-----------------------------------------------------------------------
  !> 
  !> @author  aboth
  !> @date    2021-02-25
  !> @brief   Density
  !> @details Time derivative of density
  !> 
  !-----------------------------------------------------------------------
  
  subroutine nsi_low_Mach_drho_dt(i_ord) !, fact_substep)
    integer(ip), intent(in)     :: i_ord
    !real(rp),    intent(in)     :: fact_substep
    integer(ip)                 :: ipoin

    select case (i_ord)
       
    case (-1_ip)
       !
       ! DO NOT CALCULATE
       return
       
    case (0_ip)
       !
       ! Time derivative is computed from two old values
       !
       do ipoin = 1,npoin
          drhodt_nsi(ipoin) = ( mass_rho_nsi(ipoin,1) - mass_rho_nsi(ipoin,2) ) * dtinv_nsi ! time derivative at n-0.5
       end do
       
    case (1_ip)
       ! Time derivative linearly extrapolated to current timestep: 
       ! drho/dt(s) = drho/dt(n-0.5) + (t(s)-t(n-0.5)) * d^2rho/dt^2  
       ! on paper this is equivalent to Nicoud's formulation in the final timestep 
       ! drho_dt = (((dt_inv0+dt_inv00)**2-dt_inv0**2)*Nul(ipoin) - ((dt_inv0+dt_inv00)**2)*Nu0l(ipoin) + (dt_inv0**2)*Nu00l(ipoin))/(dt_inv0*dt_inv00*(dt_inv0+dt_inv00))
       ! (F. Nicoud, Conservative high-order finite-difference schemes for low-Mach number flows,
       ! Journal of Computational Physics 158 (2000) 71-97)
       !do ipoin = 1,npoin
       !   drho_dt0   = (Nul(ipoin) -Nu0l(ipoin) ) * dtinv_nsi    ! time derivative at n-0.5
       !   drho_dt00  = (Nu0l(ipoin)-Nu00l(ipoin)) / dt_inv00     ! time derivative at n-1.5
       !   drhodt_nsi(ipoin) = drho_dt0 + (fact_substep/dtinv_nsi - 0.5_rp/dtinv_nsi)  *  (drho_dt0 - drho_dt00) / ( 0.5_rp * (dt_inv00 + 1.0_rp/dtinv_nsi) )
       !enddo

    case (2_ip)
       ! Parabola through the 3 drho/dt with t_loc=0 at n-0.5
       ! drho/dt (t) = a t^2 + b t + c
       ! 
       !           dt_inv000        dt_inv00      1/dtinv_nsi
       !       |================|==============|===============|
       !      n-3              n-2            n-1              n
       !             n-2.5           n-1.5           n-0.5    
       !          drho_dt000       drho_dt00       drho_dt0
       !         t = -del_t_2     t = -del_t_1       t = 0
       ! 
       !do ipoin = 1,npoin
       !   drho_dt0   = (Nul(ipoin) -Nu0l(ipoin) ) * dtinv_nsi    ! time derivative at n-0.5
       !   drho_dt00  = (Nu0l(ipoin)-Nu00l(ipoin)) / dt_inv00     ! time derivative at n-1.5
       !   drho_dt000 = (Nu00l(ipoin)-Nu000l(ipoin)) / dt_inv000  ! time derivative at n-2.5

       !   aa = 0.0_rp
       !   bb = 0.0_rp
       !   cc = 0.0_rp
       !   del_t_1 = 0.5_rp*(dt_inv00 + 1.0_rp/dtinv_nsi )
       !   del_t_2 = 0.5_rp*(dt_inv000 + 1.0_rp/dtinv_nsi ) + dt_inv00
       !   cc = drho_dt0
       !   aa = ( (drho_dt00-cc) - del_t_1/del_t_2 * (drho_dt000-cc) ) / ( del_t_1**2 - del_t_1*del_t_2 )
       !   bb = ( (drho_dt000-cc) - del_t_2**2 * aa  ) / (-1.0_rp * del_t_2)
       !   del_t_step = fact_substep/dtinv_nsi - 0.5_rp/dtinv_nsi
       !   drhodt_nsi(ipoin) = aa * del_t_step**2 + bb * del_t_step + cc
       !enddo
    end select
  end subroutine nsi_low_Mach_drho_dt



  subroutine nsi_low_Mach_upd_mass_rho() 
    !
    ! Save old value of mass_rho_nsi
    !
    integer(ip)                 :: ipoin

    do ipoin = 1,npoin
       mass_rho_nsi(ipoin,2) = mass_rho_nsi(ipoin,1) 
    end do
  end subroutine nsi_low_Mach_upd_mass_rho


end module mod_nsi_low_Mach
