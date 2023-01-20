!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Doiter
!> @{
!> @file    frivel.f90
!> @author  Guillaume Houzeaux
!> @brief   Wall law: Compute friction velocity
!> @details Wall law: compute the friction velocity \f$ u_* \f$ according to the 
!>          selected wall law:\n
!>
!>          - Reichardt's law: \n
!>            \f$  u^+ = \frac{1}{\kappa} \ln{(1+0.4 {y}^+)}
!>               + 7.8 \left[ 1 -\exp{\left(-\frac{{y}^+}{11}\right)}
!>               - \frac{{y}^+}{11} \exp{(-0.33 {y}^+)} \right] \f$
!>          - ABL wall law: \n
!>            \f$ u^+ = 1/ 0.41 \ln{[ (y+k_0 ) / k_0 ]} \f$ \n
!>          \n
!>          where
!>          \n
!>          - \f$  y \f$ is the wall distance;
!>          - \f$  u^+ = u/u_* \f$ is the non-dimensional velocity;
!>          - \f$  y^+ = y u_* / \nu \f$ is the non-dimensional wall distance;
!>          - \f$  k_0 \f$ is the surface roughness.\n
!>          \n
!>          The Reichardt's law is solved using the Newton-Raphson method.
!-----------------------------------------------------------------------

module mod_frivel

  use def_kintyp, only  :  ip,rp 
  implicit none

contains
  
  pure subroutine frivel(kfl_ustar,y,k,u,nu,ustar)
    !$acc routine seq
    integer(ip), intent(in)  :: kfl_ustar
    real(rp),    intent(in)  :: y          !< Wall distance
    real(rp),    intent(in)  :: u          !< Velocity
    real(rp),    intent(in)  :: nu         !< Kinematic viscosity
    real(rp),    intent(in)  :: k          !< Roughness		
    real(rp),    intent(out) :: ustar      !< Friction velocity    
    integer(ip)              :: itera
    real(rp)                 :: xmuit,fdvfr,devfr
    real(rp)                 :: vkinv,diffd,parco,yplus,onovu,yplu2
    real(rp)                 :: ypele,expye,expyt,oneoe,firsl,ypel2
    real(rp)                 :: pplus,densi,gradp,grpr2,py,sq,inv,ln4,uplus

    if( kfl_ustar == 0 ) then
       !
       ! Wall law taking into account the buffer zone and log layer:
       ! u+ =1/k*ln(1+0.4*y+)+7.8*[1-exp(-y+/11)-y+/11*exp(-0.33*y+)]   
       !   
       if( y > 0.0_rp .and. u /= 0.0_rp ) then

          ustar = sqrt( u * nu / y )
          if( ustar * y / nu > 5.0_rp ) then
             vkinv = 1.0_rp / 0.41_rp
             onovu = 1.0_rp / u
             xmuit = y / nu                                                    ! D/nu
             itera = 0                                                         ! i=0
             !          ustar = u                                                         ! U_*^0    This initialization does not make sense
             ! ustar has already been initialized correctly before
             parco = 1.0_rp
             oneoe = 1.0_rp / 11.0_rp
             do while( parco >= 1.0e-6_rp .and. itera < 100 )
                itera = itera + 1                                              ! i=i+1
                ustar = max(ustar,0.0_rp)
                yplus = ustar * xmuit
                ypele = yplus * oneoe
                ypel2 = min(ypele,20.0_rp)
                expye = exp( -ypel2 )
                yplu2 = min(yplus,70.0_rp)
                expyt = exp( -yplu2 * 0.33_rp ) 
                firsl = vkinv * log( 1.0_rp + 0.4_rp * yplus )
                fdvfr = ustar*(firsl+7.8_rp*(1.0_rp-expye-ypele*expyt))-u      ! F(U_*^i)
                diffd = firsl + vkinv*0.4_rp*yplus/(1.0_rp+0.4_rp*yplus)&      ! DF(U_*^i)
                     &  + 7.8_rp*(1.0_rp-expye*(1.0_rp-ypele)&
                     &  - ypele*expyt*(2.0_rp-yplus*0.33_rp))               
                devfr = -fdvfr / diffd                                         ! y(U_*^0)
                parco = abs( devfr * onovu )                                   ! U_*^i=U_*^i-1
                ustar = ustar + devfr                                          ! +y(U_*)
             end do

          end if

       else
          !
          ! If y = 0 or u = 0, the friction velocity is set to zero.
          !
          ustar = 0.0_rp

       end if

    else if ( ( kfl_ustar == 1 ) .or.  ( kfl_ustar == 2 ) ) then
       !
       ! ABL model: u+ = 1/k * ln( (y+y0) / y0 ); y0 is the roughness
       !
       ustar = u * 0.41_rp / log( (y+k) / k )

    else if ( kfl_ustar == 3 ) then   ! Still to be revised, tested and used
       ! It would be good to test convergence with an external program 
       ! Moreover we could use that program to test the convergence 
       ! convergence criteria of kfl_ustar==0  
       ! is abs( devfr / ustar  )   or abs( devfr * onovu  )  better 
       !
       ! Wall law with pressure gradient
       ! Continuous Formulation of Wall Function with Adverse Pressure Gradient - Thomas Rober 2013
       ! For thr moment only the logaritmic region
       ! The pressure gradient for the moment I leave it fixed to -336  -- Gimenez Channel Ret = 2003
       ! Moreover I need the density that for the moment I leave it fixed not to modify input parameters
       !
       densi = 1.14_rp
       gradp = -336.0_rp
       grpr2 = gradp * nu /densi
       ln4 = log (4.0_rp)
       if( y > 0.0_rp .and. u /= 0.0_rp ) then

          ustar = sqrt( u * nu / y )
          onovu = 1.0_rp / u
          vkinv = 1.0_rp / 0.41_rp
          itera = 0                                                         ! i=0
          ustar = u                                                         ! U_*^0
          parco = 1.0_rp
          do while( parco >= 1.0e-6_rp .and. itera < 100 )
             itera = itera + 1                                              ! i=i+1
             ustar = max(ustar,0.0_rp)
             yplus = ustar * xmuit
             pplus = grpr2 / ( ustar ** 3_ip )
             py    = pplus * yplus 
             sq    = sqrt ( py + 1.0_rp ) 
             inv   = 1.0_rp / sq
             uplus = vkinv * ( log(yplus) + 2.0_rp * (sq - 1.0_rp) + ln4 - log( 2.0_rp + py + 2.0_rp*sq ) ) + 5.5_rp
             fdvfr = ustar * uplus - u                                      ! F(U_*^i)
             diffd = uplus + vkinv * ( 1.0_rp -( 2.0_rp * py * ( inv -  &   ! DF(U_*^i)
                  ( ( 1.0_rp + inv ) / ( 2.0_rp + py + 2.0_rp * sq ) ) ) ) ) 

             devfr = -fdvfr / diffd                                         ! delta_ustar
             parco = abs( devfr / ustar  )                                  ! The convergence crieria depends on delta_ustar/ustar
             ! for kfl_ustar==0 delta_ustar/u is used 
             ! Guillaume was not totally sure of the reason 
             ustar = ustar + devfr                                          ! +y(U_*)
          end do

       else
          !
          ! If y = 0 or u = 0, the friction velocity is set to zero.
          !
          ustar = 0.0_rp

       end if

    end if

  end subroutine frivel

end module mod_frivel
!> @} 
