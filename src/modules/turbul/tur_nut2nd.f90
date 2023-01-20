!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_nut2nd(itask,ipoin,nu,nut,k,u)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_nut2nd
  ! NAME
  !   tur_nut2nd
  ! DESCRIPTION
  !    ITASK= 1 ... Compute U (2nd turb. variable) from NUT and K
  !    ITASK= 2 ... Compute NUT
  !    Useful relations:
  !    y+  = y*(u*)/nu
  !    y*  = y*sqrt(k)/nu
  !    ReT = k^2/(nu*eps)
  !    Ry  = sqrt(k)*y/nu
  ! USES
  ! USED BY
  !    tur_outvar
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_turbul
  use mod_ker_regularization, only : regul_k, regul_e, kfl_regularization

  implicit none
  integer(ip), intent(in)    :: itask,ipoin
  real(rp),    intent(in)    :: nu
  real(rp),    intent(inout) :: u,nut,k
  integer(ip)                :: ii,jbopo,nmax
  real(rp)                   :: ReT,as,uold,resid,fmu,fog,onReT,y,F2,a1
  real(rp)                   :: phi,f,g,H,yplus,Cmu,T,Ry,eps,w,toler
  real(rp)                   :: ystar,vv,ystarv,lmu,mu_inn,mu_out,wemix
  
  Cmu  = param_tur(6)

  if( itask == 1 ) then
     
     !-------------------------------------------------------------------
     !
     ! Compute U (2nd turb. variable) from NUT and K
     !
     !-------------------------------------------------------------------

     u     = 0.0_rp            ! Default value
     resid = 1.0_rp            ! Residual iterative process
     toler = 1.e-4_rp          ! Tolerance iterative process
     nmax  = 100               ! Max. number of iterations
     ii    = 0                 ! Iteration counter

  else if( itask == 2 .or. itask == 3 ) then

     !-------------------------------------------------------------------
     !
     ! Compute NUT
     !
     !-------------------------------------------------------------------

     nut = 0.0_rp              ! Default value
     k   = untur(1,ipoin,1)

     if( nturb_tur > 1 ) then
        eps = untur(2,ipoin,1)
        w   = eps
     end if
     
  end if

  if( TUR_SPALART_ALLMARAS ) then  

     !-------------------------------------------------------------------
     !
     ! Spalart-Allmaras turbulence model
     !
     ! nut = fv1*rho*nu'
     ! fv1 = X^3/(X^3+cv1^3), X = nu'/nu
     !
     !-------------------------------------------------------------------

     if( itask == 1 ) then
        call tur_spanut(1_ip,nu,u,nut)
     else
        call tur_spanut(2_ip,nu,k,nut)  ! k is nu'
     end if

  else if( TUR_K_XU_CHIEN ) then

     !-------------------------------------------------------------------
     !
     ! Xu-Chen-Nieuwstadt k-model for natural convection
     !
     ! nut = sqrt(vv)*lmu
     !
     ! vv  = k*(7.19e-3*y*-4.33e-5*(y*)**2+8.8e-8*(y*)**3)
     ! y*v = y*sqrt(vv)/nu
     ! lmu = 0.544*y/(1.0+5.025e-04*(y*v**1.65))
     !
     !-------------------------------------------------------------------

     if( itask == 1 ) then
        call runend('TUR_NUT2ND: NOT PROGRAMMED')
     else
        y      = walld(ipoin)
        ystar  = y*sqrt(k)/nu
        vv     = k*(7.19e-3_rp*ystar-4.33e-5_rp*ystar**2+8.8e-8_rp*ystar**3_rp)
        ystarv = y*sqrt(vv)/nu
        lmu    = 0.544_rp*y/(1.0_rp+5.025e-04_rp*(ystarv**1.65_rp))
        nut    = sqrt(vv)*lmu
     end if

  else if(TUR_K_EPS_STD) then

     !-------------------------------------------------------------------
     !
     ! Standard k-eps model 
     !
     ! nut = Cmu*k^2/eps
     !
     !-------------------------------------------------------------------
     if (kfl_logva==0) then
        if( itask == 1 ) then
           if(nut/=0.0_rp)   u = Cmu*k*k/nut 
        else
           if( eps > 1.0e-9_rp ) nut = Cmu * k * k / eps
        end if
        if (kfl_regularization) &
             nut = Cmu * regul_k(k)*regul_k(k)/regul_e(eps)
     else if (kfl_logva==1) then
        nut = Cmu * exp(2.0_rp*k - eps)
     end if

  else if(TUR_K_EPS_LAUNDER_SHARMA) then

     !-------------------------------------------------------------------
     !
     ! Launder-Sharma k-eps model 
     !
     ! nut = fmu*Cmu*k^2/eps
     ! fmu = exp[-3.4/(1+ReT/50)^2]
     !
     !-------------------------------------------------------------------

     if( itask == 1 ) then
        uold  = Cmu*k*k/nut
        u     = uold
        do while(ii<nmax.and.resid>toler)
           ii=ii+1
           if(u/=0.0_rp) then
              ReT   = k*k/(nu*u)
              fmu   = (1.0_rp+0.02_rp*ReT)*(1.0_rp+0.02_rp*ReT)
              fmu   = exp(-3.4_rp/fmu)
              u     = Cmu*fmu*k*k/nut
              resid = abs(uold-u)/u
              uold  = u
           else
              return
           end if
        end do
     else
        if(eps/=0.0_rp) then
           ReT = k*k/(nu*eps) 
           fmu = (1.0_rp+0.02_rp*ReT)**2
           fmu = exp(-3.4_rp/fmu)
           nut = Cmu*fmu*k*k/eps
        end if
     end if

  else if(TUR_K_EPS_CHIEN) then 

     !-------------------------------------------------------------------
     !
     ! Chien k-eps model 
     ! 
     ! nut = fmu*Cmu*k^2/eps
     ! fmu = 1-exp(-0.0115*y+)
     !
     !-------------------------------------------------------------------

     jbopo = lwnei_tur(ipoin)
     yplus = walld(ipoin)*ustar_tur(jbopo)/nu
     fmu   = 1.0_rp-exp(-0.0115_rp*yplus)
     if( itask == 1 ) then
        if(nut/=0.0_rp)  u   = Cmu*fmu*k*k/nut
     else
        if(eps/=0.0_rp)  nut = Cmu*fmu*k*k/eps
     end if

  else if(TUR_K_EPS_NAGANO_TAGAWA) then 

     !-------------------------------------------------------------------
     !
     ! Nagano-Tagawa k-eps model 
     !
     ! nut = fmu*Cmu*k^2/eps
     ! fmu = [ 1-exp(-y+/26+) ]^2 * ( 1+4.1/ReT^{3/4} )
     !
     !-------------------------------------------------------------------

     if( itask == 1 ) then
        uold  = Cmu*k*k/nut
        u     = uold
        do while(ii<nmax.and.resid>toler)
           ii=ii+1
           if(u/=0.0_rp.and.k/=0.0_rp) then
              jbopo = lwnei_tur(ipoin)
              yplus = walld(ipoin)*ustar_tur(jbopo)/nu
              onReT = (nu*eps)/(k*k)**0.75_rp
              fmu   = (1.0_rp-exp(-yplus/26.0_rp))**2*(1.0_rp+4.1_rp*onReT)
              u     = Cmu*fmu*k*k/nut
              resid = abs(uold-u)/u
              uold  = u
           else
              return
           end if
        end do
     else       
        if(k/=0.0_rp.and.eps/=0.0_rp) then
           jbopo = lwnei_tur(ipoin)
           yplus = walld(ipoin)*ustar_tur(jbopo)/nu
           onReT = (nu*eps)/(k*k)**0.75_rp
           fmu   = (1.0_rp-exp(-yplus/26.0_rp))**2*(1.0_rp+4.1_rp*onReT)
           nut   = Cmu*fmu*k*k/eps
        end if
     end if

  else if(TUR_K_EPS_LAM_BREMHORST) then 

     !-------------------------------------------------------------------
     !
     ! Lam-Bremhorst k-eps model 
     !
     ! nut = fmu*Cmu*k^2/eps
     ! fmu = [ 1.0-exp(-0.0165*Ry) ]^2 * ( 1.0+20.5/ReT )
     ! nut = Cmu*[ 1.0-exp(-0.0165*Ry) ]^2 * ( k^2/eps+20.5*nu )
     !
     !-------------------------------------------------------------------

     if( itask == 1 ) then
        uold  = Cmu*k*k/nut
        u     = uold
        do while(ii<nmax.and.resid>toler)
           ii=ii+1
           if(k/=0.0_rp) then
              onReT = (nu*u)/(k*k)
              Ry    = sqrt(k)*walld(ipoin)/nu
              fmu   = (1.0_rp-exp(-0.0165_rp*Ry))**2.0_rp*(1.0_rp+20.5_rp*onReT)
              u     = Cmu*fmu*k*k/nut
              resid = abs(uold-u)/u
              uold  = u
           else
              return
           end if
        end do
     else
        if(eps/=0.0_rp) then
           Ry  = sqrt(k)*walld(ipoin)/nu
           nut = Cmu*(1.0_rp-exp(-0.0165_rp*Ry))**2.0_rp*(k*k/eps + 20.5_rp*nu)
        end if
     end if

  else if(TUR_K_EPS_JAW_HWANG) then 

     !-------------------------------------------------------------------
     !
     ! Jaw-Hwang k-eps model 
     !
     ! nut = fmu*Cmu*T*k
     ! fmu = [ 1-exp(-Ry/70) ]^1.75
     ! T   = max( k/eps , 6.0_rp*sqrt(nu/eps) ) using a blending function
     !
     !-------------------------------------------------------------------

     if( itask == 1 ) then
        uold  = Cmu*k*k/nut
        u     = uold
        do while(ii<nmax.and.resid>toler)
           ii=ii+1
           if(u/=0.0_rp) then
              f     = k
              g     = 6.0_rp*sqrt(nu*u)
              fog   = f/g
              call mixing(1_ip,H,4.0_rp,1.0_rp,fog)
              T     = f*H+(1.0_rp-H)*g
              Ry    = sqrt(k)*walld(ipoin)/nu
              fmu   = (1.0_rp-exp(-Ry/70.0_rp))**1.75_rp          
              u     = Cmu*fmu*k*k/nut
              resid = abs(uold-u)/u
              uold  = u
           else
              call runend('NULL EPSILON')
           end if
        end do
     else
        if(eps/=0.0_rp) then
           f    = k/eps
           g    = 6.0_rp*sqrt(nu/eps)
           fog  = f/g
           call mixing(1_ip,H,4.0_rp,1.0_rp,fog)
           T    = f*H+(1.0_rp-H)*g
           Ry   = sqrt(k)*walld(ipoin)/nu
           fmu  = (1.0_rp-exp(-Ry/70.0_rp))**1.75_rp
           nut  = Cmu*fmu*T*k 
        end if
     end if

  else if(TUR_K_EPS_PHI_F) then 

     !-------------------------------------------------------------------
     !
     ! Laurence, Uribe, Utyuzhnikov k-eps-phi-f model 
     !
     ! nut = Cmu*phi*k*T
     !
     !-------------------------------------------------------------------

     if( itask == 1 ) then
        uold  = Cmu*k*k/nut
        u     = uold
        do while(ii<nmax.and.resid>toler)
           ii=ii+1
           if(uold/=0.0_rp) then
              phi = untur(3,ipoin,1) 
              f     = k
              g     = 6.0_rp*sqrt(nu*u)
              call mixmax(4.0_rp,f,g,T)
              u     = Cmu*k*T*phi/nut
              resid = abs(uold-u)/u
              uold  = u
           else
              return
           end if
        end do
     else
        if(ittim>=inits_tur) then
           if(eps/=0.0_rp) then
              phi = untur(3,ipoin,1) 
              f   = k/eps
              g   = 6.0_rp*sqrt(nu/eps)
              call mixmax(4.0_rp,f,g,T)
              nut = Cmu*k*phi*T 
           end if
        else
           if(eps/=0.0_rp) then
              Ry  = sqrt(k)*walld(ipoin)/nu
              nut = Cmu*(1.0_rp-exp(-0.0165_rp*Ry))**2.0_rp&
                   &  *(k*k/eps + 20.5_rp*nu)
           end if
        end if
     end if

  else if(TUR_TWO_LAYER_RODI) then

     !-------------------------------------------------------------------
     !
     ! Two-layer Rodi k-eps model 
     !
     ! Outer layer: nut = Cmu*k^2/eps
     ! Inner layer: nut = sqrt(vv)*lmu
     !
     !-------------------------------------------------------------------

     if( itask == 1 ) then
        call runend('TUR_NUT2ND: NOT PROGRAMMED')
     else        
        y      = walld(ipoin)
        jbopo  = lwnei_tur(ipoin)
        yplus  = y*ustar_tur(jbopo)/nu
        ystar  = y*sqrt(k)/nu
        vv     = k*0.4_rp*(1.0_rp-exp(-ystar**2/4200.0_rp))
        ystarv = walld(ipoin)*sqrt(vv)/nu
        lmu    = 0.33_rp*y/(1.0_rp+0.0005025_rp*(ystarv**1.53_rp))
        call mixing(1_ip,wemix,param_tur(8),param_tur(9),ystar)
        mu_inn = sqrt(vv)*lmu
        mu_out = 0.0_rp
        if(eps/=0.0_rp) mu_out = Cmu*k*k/eps
        nut    = wemix*mu_out+(1.0_rp-wemix)*mu_inn
     end if

  else if(TUR_TWO_LAYER_XU_CHEN) then

     !-------------------------------------------------------------------
     !
     ! Two-layer Xu-Chen k-eps model for natural convection
     !
     ! Outer layer: nut = Cmu*k^2/eps
     ! Inner layer: nut = sqrt(vv)*lmu
     !
     !-------------------------------------------------------------------

     if( itask == 1 ) then
        call runend('TUR_NUT2ND: NOT PROGRAMMED')
     else        
        jbopo  = lwnei_tur(ipoin)
        y      = walld(ipoin)
        yplus  = y*ustar_tur(jbopo)/nu
        ystar  = y*sqrt(untur(1,ipoin,1))/nu
        vv     = k*(7.19e-03_rp*ystar-4.33e-05_rp*ystar**2+8.8e-08_rp*ystar**3)
        ystarv = y*sqrt(vv)/nu
        lmu    = 0.544_rp*y/(1.0_rp+0.0005025_rp*(ystarv**1.65_rp))
        call mixing(1_ip,wemix,param_tur(8),param_tur(9),ystar)
        mu_inn = sqrt(vv)*lmu
        mu_out = 0.0_rp
        if(eps/=0.0_rp) mu_out = Cmu*k*k/eps
        nut    = wemix*mu_out+(1.0_rp-wemix)*mu_inn
     end if

  else if(TUR_K_OMEGA) then         
     
     !-------------------------------------------------------------------
     !
     ! Wilcox k-w model 
     !
     ! nut = (alpha*)*k/w
     ! a*  = ( 1/40+ReT/6 ) / ( 1+ReT/6 )
     !
     !-------------------------------------------------------------------

     if( ipara_tur(1) == 1 ) then
        if( itask == 1 ) then
           u  = k/nut
        else
           nut = k/w
        end if
     else
        if( itask == 1 ) then
           uold  = k/nut
           u     = uold
           do while(ii<nmax.and.resid>toler)
              ii=ii+1
              if(u/=0.0_rp) then
                 ReT   = k/(nu*u)
                 as    = (param_tur(10)+ReT/param_tur(7))/(1.0_rp+ReT/param_tur(7))
                 u     = as*k/nut
                 resid = abs(uold-u)/u
                 uold  = u
              else
                 return
              end if
           end do
        else
           if(w/=0.0_rp) then
              ReT = k/(nu*w)
              as  = (param_tur(10)+ReT/param_tur(7))/(1.0_rp+ReT/param_tur(7))
              nut = as*k/w
           end if
        end if
     end if

  else if(TUR_K_OMEGA_BREDBERG) then

     !-------------------------------------------------------------------
     !
     ! Bredberg-Peng-Davidson k-w model 
     !
     ! nut = Cmu*fmu*k/w
     ! fmu = 0.09 + [ 0.91 + 1/(ReT^3) ]*[ 1 -exp( -(ReT/25)^2.75 ) ]
     !
     !-------------------------------------------------------------------

     if( itask == 1 ) then
        uold  = k/nut
        u     = uold
        do while(ii<nmax.and.resid>toler)
           ii=ii+1
           if(uold/=0.0_rp) then
              ReT   = k/(nu*u)
              if(ReT>0.001_rp) then
                 fmu = 0.09_rp+(0.91_rp+1.0_rp/(ReT*ReT*ReT))&
                      &  *(1.0_rp-exp(-(ReT/25.0_rp)**2.75_rp))
              else
                 fmu = 0.09_rp
              end if
              as    = (param_tur(10)+ReT/param_tur(7))/(1.0_rp+ReT/param_tur(7))
              u     = fmu*param_tur( 6)*k/nut
              resid = abs(uold-u)/u
              uold  = u
           else
              return
           end if
        end do
     else
        if(k/=0.0_rp.and.w/=0.0_rp) then 
           ReT = k/(nu*w)
           if(ReT>0.001_rp) then
              fmu = 0.09_rp+(0.91_rp+1.0_rp/(ReT*ReT*ReT))&
                   *(1.0_rp-exp(-(ReT/25.0_rp)**2.75_rp))
           else
              fmu = 0.09_rp
           end if
           nut = fmu*param_tur( 6)*k/w
        end if
     end if

  else if(TUR_SST_K_OMEGA) then         
     
     !-------------------------------------------------------------------
     !
     ! Menter SST k-w model 
     !
     ! nut = rho*a1*k/max(a1*w,Omega*F2)
     !
     !-------------------------------------------------------------------

     if( itask == 1 ) then 

        u = k/nut

     else if( itask == 2 ) then 

        y  = walld(ipoin) 
        a1 = 0.31_rp
        if( y > 0.0_rp ) then   
           F2 = max( 2.0_rp * sqrt(max(0.0_rp,k)) / ( 0.09_rp * w * y ) , 500.0_rp * nu / ( y * y * w ) )
           F2 = tanh( F2 * F2 )
           as = max( a1 * w , vorti_tur(ipoin) * F2 )
           if( ( as > 1.0e-10_rp ) .or. ( kfl_clipp_tur > 1 ) ) then
              nut = a1 * k / as
           end if
           sstf2_tur(ipoin) = F2
        else
           sstf2_tur(ipoin) = 0.0_rp
        end if

     else if( itask == 3 ) then 

        y  = walld(ipoin) 
        a1 = 0.31_rp
        if( y > 0.0_rp ) then   
           F2 = max( 2.0_rp * sqrt(max(0.0_rp,k)) / ( 0.09_rp * w * y ) , 500.0_rp * nu / ( y * y * w ) )
           F2 = tanh( F2 * F2 )
           as = max( a1 * w , vorti_tur(ipoin) * F2 )
           u  = as
        else
           u  = 0.0_rp
        end if

     end if

  end if

end subroutine tur_nut2nd
