!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @defgroup Integration_Rules_Toolbox
!> @{
!> @name    ToolBox for integration rules
!> @file    mod_elmgeo.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for integration rules 
!> @details Different integration rules to compute integrals 
!
!-----------------------------------------------------------------------

module mod_integration_rules

  use def_kintyp, only : ip,rp,lg
  implicit none
  private

  public :: integration_rules_trapezoidal

contains

  subroutine integration_rules_trapezoidal(ndime,ngaus,posgp,weigp,ierro)
    integer(ip), intent(in)            :: ndime
    integer(ip), intent(in)            :: ngaus
    real(rp),    intent(out)           :: posgp(ndime,ngaus)
    real(rp),    intent(out)           :: weigp(ngaus)
    integer(ip), intent(out), optional :: ierro
    integer(ip)                        :: igaus,jgaus,kgaus,ngaus1
    real(rp)                           :: xfact1

    if( present(ierro) ) ierro = 0

    if( ndime == 1 ) then
       if( ngaus == 1 ) then
          posgp(1,1) = 0.0_rp
          weigp(1)   = 2.0_rp
       else
          xfact1 = 2.0_rp / real(ngaus-1,rp)       
          do igaus = 1,ngaus
             posgp(1,igaus) = -1.0_rp + real(igaus-1,rp) * xfact1
             weigp(igaus)   =  xfact1
          end do
          weigp(1)     = 0.5_rp * xfact1 
          weigp(ngaus) = 0.5_rp * xfact1
       end if
    else if( ndime == 2 ) then
       if( ngaus == 1 ) then
          posgp(1,1) = 0.0_rp
          weigp(1)   = 2.0_rp
       else
          ngaus1 = int(sqrt(real(ngaus,rp)),ip)
          xfact1 = 2.0_rp / real(ngaus1-1,rp)    
          if( ngaus1**2 /= ngaus ) call runend('INTEGATION_RULES_TRAPEZOIDAL: WRONG NUMBER OF GAUSS POINTS') 
          kgaus  = 0
          do igaus = 1,ngaus1
             do jgaus = 1,ngaus1
                kgaus = kgaus + 1
                posgp(1,kgaus) = -1.0_rp + real(igaus-1,rp) * xfact1
                posgp(2,kgaus) = -1.0_rp + real(jgaus-1,rp) * xfact1
                weigp(kgaus)   =  xfact1
                if( jgaus == 1 .or. jgaus == ngaus1 )  weigp(kgaus) = 0.5_rp * xfact1 
             end do
             if( jgaus == 1 .or. jgaus == ngaus1 )  weigp(kgaus) = 0.5_rp * xfact1 
          end do
          if( igaus == 1 .or. igaus == ngaus1 )  weigp(kgaus) = 0.5_rp * xfact1 
       end if
    else
       call runend('INTEGATION_RULES_TRAPEZOIDAL: WRONG INTEGRATION RULE')
    end if

  end subroutine integration_rules_trapezoidal

end module mod_integration_rules
!> @}
