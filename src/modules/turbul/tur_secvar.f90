!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_secvar(itask,works_tur)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_secvar
  ! NAME
  !   tur_secvar
  ! DESCRIPTION
  !    Compute turbulence secondary variables
  ! USES
  ! USED BY
  !    tur_outvar
  !***
  !-----------------------------------------------------------------------
  use  def_parame
  use  def_master
  use  def_domain
  use  def_turbul
  use mod_ker_proper
  implicit none
  integer(ip), intent(in)  :: itask
  real(rp),    intent(out) :: works_tur(npoin)
  integer(ip)              :: ipoin,ibopo,jpoin,dummi
  real(rp)                 :: k,w,bs,ReT,eps,leps,y,ystar,ystarv,vv,nu,rho(1),mu(1)

  select case (itask)

  case(1_ip)  !when TUR_FAMILY_K_OMEGA
     !
     ! eps from w
     !
     do ipoin=1,npoin
        k = untur(1,ipoin,1)
        w = untur(2,ipoin,1)        
        call ker_proper('DENSI','IPOIN',ipoin,dummi,rho(1:))
        call ker_proper('VISCO','IPOIN',ipoin,dummi,mu(1:))
      
        if(w.gt.epsilon(1.0_rp)) then
           ReT=rho(1)*k/(w*mu(1))
        else
           ReT=0.0_rp
        end if
        bs    = param_tur(6) !&
!             &  *(5.0_rp/18.0_rp+(ReT/param_tur(9))**4.0_rp)&
!             &  /(1.0_rp        +(ReT/param_tur(9))**4.0_rp)
        works_tur(ipoin)=bs*w*k
     end do

  case(2_ip)
     !
     ! w from eps
     !
     do ipoin=1,npoin
        k=untur(1,ipoin,1)
        eps=untur(2,ipoin,1)
        if(k>zetur) works_tur(ipoin)=eps/(param_tur(6)*k)
     end do

  case(3_ip)
     !
     ! l from eps
     !
     do ipoin=1,npoin
        k=untur(1,ipoin,1)
        eps=untur(2,ipoin,1)
        if(eps>zetur) works_tur(ipoin)=param_tur(6)*k**1.5_rp/eps
     end do

  case(4_ip)
     !
     ! l from w
     !
     do ipoin=1,npoin
        k=untur(1,ipoin,1)
        w=untur(2,ipoin,1)
        if(w>zetur) works_tur(ipoin)=param_tur(6)*sqrt(k)/w
     end do

  case(5_ip)
     !
     ! eps from leps: eps=(vv)^{1/2)*k/leps
     !
     call runend('TUR_SECAVR: NOT PROGRAMMED')
     do ipoin=1,npoin
        ibopo  = lpoty(ipoin)
        y      = walld(ipoin)
        if(y<1.0e-10_rp.and.ibopo/=0) then
           !jpoin = lfiel_tur(ibopo)
           y     = walld(jpoin)
        else
           jpoin = ipoin
        end if
        call ker_proper('DENSI','IPOIN',ipoin,dummi,rho(1:))
        call ker_proper('VISCO','IPOIN',ipoin,dummi,mu(1:))
       
        k      = untur(1,jpoin,1) 
        nu     = mu(1)/rho(1)
        ystar  = y*sqrt(k)/nu                                 ! y* = y*K^{1/2}/nu
        vv     = k*(7.19e-03_rp*ystar-4.33e-05_rp*ystar**2_rp&
             &          +8.8e-08_rp*ystar**3_rp)         
        ystarv = y*sqrt(vv)/nu                                ! (yv*) = y*vv^{1/2}/nu        
        if(ystarv<1.0e-10_rp) then
           works_tur(ipoin)=0.0_rp
        else
           leps  = 8.8_rp*y/(1.0_rp+10.0_rp/ystarv&           ! leps = 8.8*y/[1+10/(yv*)+0.0515*(yv*)]
                &  +0.0515_rp*ystarv)
           works_tur(ipoin)=sqrt(vv)*k/leps
        end if
     end do

  case(6_ip)
     !
     ! leps
     !
     do ipoin=1,npoin
        ibopo  = lpoty(ipoin)
        y      = walld(ipoin)
        if(y<1.0e-10.and.ibopo/=0) then
           !jpoin = lfiel_tur(ibopo)
           y     = walld(jpoin)
        else
           jpoin = ipoin
        end if
        k      = untur(1,jpoin,1) 
      
        call ker_proper('DENSI','IPOIN',ipoin,dummi,rho(1:))
        call ker_proper('VISCO','IPOIN',ipoin,dummi,mu(1:))
        
        nu     = mu(1)/rho(1)
        ystar  = y*sqrt(k)/nu                                 ! y* = y*K^{1/2}/nu
        vv     = k*(7.19e-03_rp*ystar-4.33e-05_rp*ystar**2_rp&
             &          +8.8e-08_rp*ystar**3_rp)         
        ystarv = y*sqrt(vv)/nu                                ! (yv*) = y*vv^{1/2}/nu        
        if(ystarv<1.0e-10_rp) then
           works_tur(ipoin)=0.0_rp
        else
           leps  = 8.8_rp*y/(1.0_rp+10.0_rp/ystarv&           ! leps = 8.8*y/[1+10/(yv*)+0.0515*(yv*)]
                &  +0.0515_rp*ystarv)
           works_tur(ipoin)=leps
        end if
     end do

  case(7_ip)
     !
     ! yv*
     !
     do ipoin=1,npoin
        y     = walld(ipoin)
        k     = untur(1,ipoin,1)
       
        call ker_proper('DENSI','IPOIN',ipoin,dummi,rho(1:))
        call ker_proper('VISCO','IPOIN',ipoin,dummi,mu(1:))
        nu    = mu(1)/rho(1)
        ystar = y*sqrt(k)/nu
        vv    = k*(7.19e-03_rp*ystar-4.33e-05_rp*ystar**2_rp+8.8e-08_rp*ystar**3_rp)         
        works_tur(ipoin)=y*sqrt(vv)/nu
     end do

  case(8_ip)
     !
     ! eps = eps0 + eps'
     !
     do ipoin=1,npoin
        if(walld(ipoin)/=0.0_rp) then
           call ker_proper('DENSI','IPOIN',ipoin,dummi,rho(1:))
           call ker_proper('VISCO','IPOIN',ipoin,dummi,mu(1:))
            works_tur(ipoin)=untur(2,ipoin,1)&
                +2.0_rp*mu(1)/rho(1)&
                *untur(1,ipoin,1)/(walld(ipoin)*walld(ipoin))
        else
           works_tur(ipoin)=untur(2,ipoin,1)
        end if
     end do

  case(9_ip)
     ! 
     ! y+: dimensionless distance to the wall
     !  
     do ipoin=1,npoin
        ibopo=lpoty(ipoin)
        if(ibopo/=0.and.kfl_fixno_tur(1,ipoin,1)==3) then
           call ker_proper('DENSI','IPOIN',ipoin,dummi,rho(1:))
           call ker_proper('VISCO','IPOIN',ipoin,dummi,mu(1:))
           works_tur(ipoin) = walld(ipoin)*rho(1)*ustar_tur(ibopo)/mu(1)
        else
           works_tur(ipoin)=0.0_rp
        end if
     end do
     
  end select

end subroutine tur_secvar
