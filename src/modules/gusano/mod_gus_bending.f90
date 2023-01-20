!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    gus_bending.f90
!> @author  houzeaux
!> @date    2020-10-26
!> @brief   Impose bending
!> @details Impose bending pressure drop
!-----------------------------------------------------------------------

module mod_gus_bending

  use def_kintyp_basic, only : ip,rp
  use def_master
  use def_domain
  use def_parame
  use def_gusano
  
  implicit none
  private

  public :: gus_bending

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-26
  !> @brief   Bending model
  !> @details Bending model
  !>          Curved laminar pipe. Formula by Ito:
  !>
  !>                               21.5 De
  !>          K / K_straight = ----------------
  !>                           [ 1.56 + log10 De ]^5.73
  !>
  !>          [1] Ito, H. (1969). Laminar flow in curved pipes. Z.
  !>              Angew Math Mech. 49 653-663. 
  !>          [2] Ito, H. (1959). Friction factors for turbulent flow in
  !>              curved pipes. Trans. ASME J. Basic Eng. 81D 123-134
  !>          [3] https://oatao.univ-toulouse.fr/20040/1/Spedding_20040.pdf
  !>
  !-----------------------------------------------------------------------

  subroutine gus_bending_go(aa)

    real(rp),   intent(inout) :: aa(2,2,*)
    integer(ip)               :: iz,ipoin
    real(rp)                  :: r,De,Re,Rc,u,rho,mu,K,D
    real(rp)                  :: alpha
    !
    ! Compute properties
    !
    if( kfl_bendi_gus == GUS_BEND_LAMINAR_PIPE ) then
       
       do ipoin = 1,npoin
          if( abs(angle_gus(ipoin)) > 0.0_rp ) then
             loop_iz1: do iz = r_dom(ipoin),r_dom(ipoin+1)-1                
                if( c_dom(iz) == ipoin ) then
                   rho        = densi_gus(ipoin)
                   mu         = visco_gus(ipoin)
                   u          = abs(vel1d(ipoin,1))
                   r          = sqrt(areas(ipoin,1)/pi)
                   D          = 2.0_rp * r
                   Rc         = bendi_gus(1,ipoin)
                   Re         = rho*u*D/mu
                   De         = Re * sqrt(r/Rc)
                   if( De > 13.5_rp .and. De < 2000.0_rp ) then
                      alpha      = 21.5_rp * De / (1.56_rp+log10(De))**5.73_rp
                      K          = 32.0_rp * mu / D
                      aa(1,1,iz) = aa(1,1,iz) + 0.5_rp*rho*u*K
                   end if
                   exit loop_iz1
                end if
             end do loop_iz1
          end if
       end do
       
    else if( kfl_bendi_gus == GUS_BEND_TURBULENT_PIPE ) then
       
       do ipoin = 1,npoin
          if( abs(angle_gus(ipoin)) > 0.0_rp ) then
             loop_iz2: do iz = r_dom(ipoin),r_dom(ipoin+1)-1                
                if( c_dom(iz) == ipoin ) then
                   rho        = densi_gus(ipoin)
                   mu         = visco_gus(ipoin)
                   u          = abs(vel1d(ipoin,1))
                   r          = sqrt(areas(ipoin,1)/pi)
                   Rc         = bendi_gus(1,ipoin)
                   D          = 2.0_rp * r
                   Re         = rho*u*D/mu
                   De         = Re * sqrt(r/Rc)
                   alpha      = 0.003625_rp + 0.038_rp * ( Re*r*r/(R*R) )**(-0.25_rp)
                   K          = alpha * sqrt(r/R)
                   aa(1,1,iz) = aa(1,1,iz) + 0.5_rp*rho*u*K  
                   exit loop_iz2
                end if
             end do loop_iz2
          end if
       end do
    end if

  end subroutine gus_bending_go

  subroutine gus_bending(aa)
    
    real(rp), pointer, intent(inout) :: aa(:)

    if( associated(aa) ) call gus_bending_go(aa)
    
  end subroutine gus_bending
  
end module mod_gus_bending
!> @}
