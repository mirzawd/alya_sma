!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> 
!> @author  houzeaux
!> @date    2018-09-25
!> @brief   Brownian motion
!> @details Compute some variable to impose Brownian motion as
!>          a displacement or as a force
!> 
!-----------------------------------------------------------------------

subroutine pts_brown(itint,ielem,itype,pnode,shapf,dt_k,r,mu,rho,Cc,D,epsi1,epsi2,epsi3)

  use def_master, only : therm
  use def_kintyp, only : ip,rp
  use def_parame, only : kb,pi
  use def_domain, only : lnods,mnode,ndime
  use mod_random, only : random_generate_number
  use def_partis, only : parttyp

  implicit none

  integer(ip), intent(in)    :: ielem,itype,itint
  integer(ip), intent(in)    :: pnode
  real(rp),    intent(in)    :: dt_k,r,Cc,mu,rho
  real(rp),    intent(in)    :: shapf(mnode)
  real(rp),    intent(inout) :: D
  real(rp),    intent(out)   :: epsi1,epsi2,epsi3
  real(rp)                   :: Temp,CBrown,denpa,diame,U1,U2,U3,U4,S_0
!  real(rp)                   :: mass
  real(rp),           save   :: epsi1_save
  real(rp),           save   :: epsi2_save
  real(rp),           save   :: epsi3_save

  Temp  = 293.15_rp   ! 20 degrees

  denpa = parttyp(itype) % denpa
  diame = parttyp(itype) % diame 
  !
  ! Random numbers must only be generated during the
  ! first iteration in order not to modify the statistcs
  !
  if( itint == 1 )then
     if( ndime == 2 )then
        U1    = random_generate_number(broadcast_seed=.false.)
        U2    = random_generate_number(broadcast_seed=.false.)
        if(U1==0.0_rp)then
           epsi1 = 0.0_rp
           epsi2 = 0.0_rp
           epsi3  = 0.0_rp
        else
           epsi1  = sqrt( -2.0_rp * log(U1) ) * cos(2.0_rp*pi*U2)
           epsi2  = sqrt( -2.0_rp * log(U1) ) * sin(2.0_rp*pi*U2)
           epsi3  = 0.0_rp
        end if
     else
        U1    = random_generate_number(broadcast_seed=.false.)
        U2    = random_generate_number(broadcast_seed=.false.)
        U3    = random_generate_number(broadcast_seed=.false.)
        U4    = random_generate_number(broadcast_seed=.false.)
        if( U1==0.0_rp .or. U3==0.0_rp)then
           epsi1  = 0.0_rp
           epsi2  = 0.0_rp
           epsi3  = 0.0_rp
        else
           epsi1  = sqrt( -2.0_rp * log(U1) ) * cos(2.0_rp*pi*U2)
           epsi2  = sqrt( -2.0_rp * log(U1) ) * sin(2.0_rp*pi*U2)
           epsi3  = sqrt( -2.0_rp * log(U3) ) * cos(2.0_rp*pi*U4)
        end if
     end if
     epsi1_save = epsi1 
     epsi2_save = epsi2 
     epsi3_save = epsi3 
  else
     epsi1 = epsi1_save 
     epsi2 = epsi2_save
     epsi3 = epsi3_save    
  end if
  !
  ! If necessary, calculate diffusivity
  ! kbolt = 1.38065e-23
  !
  if( D == 0.0_rp )then
     if( associated(therm) )then
        Temp = dot_product(shapf(1:pnode),therm(lnods(1:pnode,ielem),1))
     end if
     D  = kb * Temp * Cc / (6.0_rp*pi*mu*r)                                
  end if
  !
  ! Calculte the Brownian diffusion coefficient
  !
  if( parttyp(itype) % kfl_brown == 1 ) then
     !   
     ! Brownian Motion (from Einstein 1905; rms displacement)
     !
     CBrown = sqrt(2.0_rp*D*dt_k)
     ! modif from Longest and Xi for 2D 
     !CBrown = sqrt(4.0_rp*D*dt_k)
  else
     !
     ! Brownian Force (from Longest and Xi, Li and Ahmadi)
     !
     !!!!!! CBrown = 3.0_rp/(4.0_rp*pi*r**3*denpa) * sqrt(2.0_rp*kb**2*Temp**2/(D*dt_k))
     S_0 = 216.0_rp*mu/rho*kb*Temp/(pi*pi*rho*diame**5*(denpa/rho)**2*Cc)
     Cbrown = sqrt(pi*S_0/dt_k)

     !mass   = pi*diame**3/6.0_rp
     !S_0    = sqrt(2.0_rp*D)/(4.0*mass)
     !Cbrown = S_0 * dt_k**(-1.5_rp)
     
  end if
  if( ndime == 2 )then
     epsi1 = epsi1 * CBrown
     epsi2 = epsi2 * CBrown
     epsi3 = 0.0_rp
  else
     epsi1 = epsi1 * CBrown
     epsi2 = epsi2 * CBrown
     epsi3 = epsi3 * Cbrown
  end if

end subroutine pts_brown
