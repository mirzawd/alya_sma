!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine extefo(imode,coord,denpa,spher,denfl,visfl,t,accel)
  !-----------------------------------------------------------------------
  !****f* mathru/extefo
  ! NAME
  !    extefo
  ! DESCRIPTION
  !    Compute an external acceleration
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_parame, only     :  pi
  use def_domain, only     :  ndime
  use def_kermod, only     :  grnor,gravi
  implicit none
  integer(ip), intent(in)  :: imode
  real(rp),    intent(in)  :: coord(ndime)
  real(rp),    intent(in)  :: denpa
  real(rp),    intent(in)  :: spher
  real(rp),    intent(in)  :: denfl
  real(rp),    intent(in)  :: visfl
  real(rp),    intent(in)  :: t
  real(rp),    intent(out) :: accel(ndime)
  integer(ip)              :: idime
  real(rp)                 :: B0,M,w,z0,rhoap,xl
!  real(rp)                 :: xr
  real(rp)                 :: xfact

  accel(1:ndime) = 0.0_rp
 
  select case ( imode )

  case ( 0_ip )
     !
     ! No force
     !
     return

  case ( 1_ip )
     !
     ! Magnet for vasito experiment: rhoap is apparent density
     !
     B0    = 0.6_rp
     M     = 300.0_rp
     z0    = 0.04_rp
     w     = 0.24_rp
     rhoap = 2.0_rp*pi/w*B0*M/grnor*exp(-2.0_rp*pi*(coord(2)+z0)/w)
     do idime = 1,ndime
        accel(idime) = - grnor * gravi(idime) * rhoap / denpa
     end do
     accel    = 0.0_rp
     accel(2) = 1.0_rp / denpa

  case ( 2_ip )
     !
     ! Magnet: apparent water density goes down (Bin's experiment)
     ! Equilibrium hiehgt is at 0.076338
     !
     !denfl = 1007.0_rp
     !denpa = 905.0_rp
     B0    = 0.6_rp
     M     = 470.0_rp
     z0    = 0.2_rp
     w     = 0.24_rp
     rhoap  = - 2.0_rp*pi/w*B0*M/grnor*exp(2.0_rp*pi*(coord(2)-z0)/w)
     do idime = 1,ndime
        accel(idime) = - grnor * gravi(idime) * rhoap / denpa
     end do

  case ( 3_ip )
     !
     ! Magnet: prototype (old MDS)
     !
     !denfl = 1007.0_rp
     !denpa = 905.0_rp
     B0    =  0.6_rp
     M     =  470.0_rp
     z0    =  0.165_rp
     w     =  0.24_rp
     xl    = -0.09_rp
     !xr    =  0.7_rp
     !xr    =  2.7_rp
     !xfact =  ( 1.0_rp+exp( -(coord(1)-xl)/(z0-coord(2)) ) ) * ( 1.0_rp+exp( (coord(1)-xr)/(z0-coord(2)) ) )
     xfact =   1.0_rp
     rhoap = - 2.0_rp*pi/w*B0*M/grnor*exp(2.0_rp*pi*(coord(2)-z0)/w)/xfact
     do idime = 1,ndime
        accel(idime) = - grnor * gravi(idime) * rhoap / denpa
     end do

  case ( 4_ip )
     !
     ! Magnet for vasito experiment: rhoap is apparent density
     !
     B0    = 0.6_rp
     M     = 300.0_rp
     z0    = 0.04_rp
     w     = 0.24_rp
     rhoap = B0*M*exp(-2.0_rp*pi*(coord(2)+z0)/w)
     !accel(1) = rhoap
     accel(1) = coord(2)

  case ( 5_ip )
     !
     ! Magnet: new prototype (new MDS)
     !
     !denfl = 1007.0_rp
     !denpa = 905.0_rp
     B0    =  0.6_rp
     M     =  470.0_rp
     z0    =  0.06_rp
     w     =  0.24_rp
     xl    = -0.09_rp
     !xr    =  0.7_rp
     !xr    =  2.7_rp
     !xfact =  ( 1.0_rp+exp( -(coord(1)-xl)/(z0-coord(2)) ) ) * ( 1.0_rp+exp( (coord(1)-xr)/(z0-coord(2)) ) )
     xfact =   1.0_rp
     rhoap = - 2.0_rp*pi/w*B0*M/grnor*exp(2.0_rp*pi*(coord(2)-z0)/w)/xfact
     do idime = 1,ndime
        accel(idime) = - grnor * gravi(idime) * rhoap / denpa
     end do

  case ( 6_ip )
     !
     ! Magnet Damping test in 2D
     ! B0 [T] = [N/Am]
     ! M  [A/m]
     !
     B0       =   0.6_rp
     M        = 597.7_rp
     z0       =  -0.0105_rp
     w        =   0.24_rp
     accel(2) =   2.0_rp*pi/w*B0*M*exp(-2.0_rp*pi*(coord(2)-z0)/w) / denpa
     !accel(3) =   2.0_rp*pi/w*B0*M*exp(-2.0_rp*pi*(coord(3)-z0)/w) / denpa
  case default
     !
     ! Others
     !
     call runend('EXTEFO: NON-EXISTING MODEL')

  end select

end subroutine extefo
