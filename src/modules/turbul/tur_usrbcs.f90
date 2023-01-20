!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_usrbcs(iprob) 
  !------------------------------------------------------------------------
  !****f* Turbul/tur_usrbcs
  ! NAME 
  !    tur_usrbcs
  ! DESCRIPTION
  !    This routine sets user boundary contitions. There exists some
  !    predefined problems:
  !    IPROB ... 0: program yourself
  !          ... 1: flow over cylinder (0,0) of radius 0.5
  !                 x=[-6,12], y=[-12,12]
  ! USES
  ! USED BY
  !    tur_reabcs
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_turbul
  use mod_memchk
  use mod_ker_proper
  implicit none
  integer(ip), intent(in) :: iprob
  integer(ip)             :: ipoin,dummi
  real(rp)                :: ustar,nu,y,yplus,Res,kap,mu(1),rho(1),k,eps

  select case(iprob)

  case(0_ip)
     !
     ! Program yourself
     !

  case(1_ip)
     !
     ! Program yourself
     !


  case(2_ip)
     !
     ! Fully developed channel flow: [0,2]
     !
     do ipoin=1,npoin
        call ker_proper('DENSI','IPOIN',ipoin,dummi,rho(1:))
        call ker_proper('VISCO','IPOIN',ipoin,dummi,mu(1:))
       
        nu = mu(1)/rho(1)
        if(coord(2,ipoin)<1.0_rp) then
           y = coord(2,ipoin)
        else
           y = 2.0_rp - coord(2,ipoin)
        end if
        ustar = 1.0_rp
        yplus = y*ustar/nu

        Res   = rho(1)*ustar*hdiam_tur/mu(1)                           ! Re* = rho*(u*)*D/mu
        k     = ustar**2*(0.07_rp*yplus**2*exp(-yplus/8.0_rp)&   ! k   = u*^2*[0.07y+*exp(-y+/8)
             &  +4.5_rp*(1.0_rp-exp(-yplus/20.0_rp))&            !       +4.5*(1-exp(-y+/20))
             &  /(1.0_rp+4.0_rp*yplus/Res))                      !       /(1+4y+/Re*)
        eps   = ustar**4/(nu*kap&                                ! eps = u*^4/(nu*kap*(y+^4+15^4)^{1/4})
             &  *(yplus**4_rp+15.0_rp**4_rp)**0.25_rp)

     end do

  case(3_ip)
     !
     ! NACA0015_6
     !    
     do ipoin = 1,npoin
        if( lpoty(ipoin) /= 0 ) then
           if ( coord(1,ipoin) > -0.1_rp .and. coord(1,ipoin) < 0.5_rp .and. &
                &    coord(2,ipoin) > -0.1_rp .and. coord(2,ipoin) < 0.1_rp ) then             
              kfl_fixno_tur(1,ipoin,1) = 4      
           else if( coord(1,ipoin) <= 0.296001_rp ) then
              kfl_fixno_tur(1,ipoin,1) = 6           
           end if
        end if
     end do

  case(4_ip)
     !
     ! K1
     !
     do ipoin = 1,npoin
        if( lpoty(ipoin) /= 0 ) then

           if (   coord(1,ipoin) > -3.69457_rp .and. coord(2,ipoin) < -3.7869_rp ) then     
        
              kfl_fixno_tur(1,ipoin,:) = 0        

           else if (   coord(1,ipoin) > -3.69457_rp .and. coord(2,ipoin) > 3.2509_rp ) then     
        
              kfl_fixno_tur(1,ipoin,:) = 0             
              
           end if
        end if
     end do

  end select

end subroutine tur_usrbcs
