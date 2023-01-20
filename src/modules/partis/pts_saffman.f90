!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Partis
!! @{
!> @name    Saffman force
!! @file    pts_saffma.f90
!> @author  Vincent Boyer
!> @date    14/07/2013
!! @brief   This routine calculates Saffman's force
!! @details Determine for each configuration which model should be used for Saffman force
!! 
!! Formula based on paper "On the role of lift force in turbulence simulation of particule deposition" by WANG, SQUIRES, CHEN & Mc LAUGHLIN
!!
!! First detemine if particle is close to wall : just check if it is in element belonging to boundary
!!
!------------------------------------------------------------------------

subroutine pts_saffma(v_pt,v_fl,d,dens_fl,visc_k,t_relax,epsi,saff_denom,in_wall,fsaff,wdist)
  use def_parame
  use def_master
  use def_kermod
  use def_partis
  use def_domain
  use mod_ker_proper
  use mod_random
  use mod_elmgeo
  use mod_memory
  implicit none
  
  integer(ip)              :: idime,jdime
!  integer(ip)              :: pnode
!  integer(ip)              :: sielem,sielel,sjelem,spelty,sptopo,sinode
  real(rp),    intent(in)  :: v_pt(3)
  real(rp),    intent(in)  :: v_fl(3)
  real(rp),    intent(in)  :: d
  real(rp),    intent(in)  :: dens_fl         ! fluid density
  real(rp),    intent(in)  :: visc_k
  real(rp),    intent(in)  :: t_relax
  real(rp),    intent(in)  :: saff_denom
  real(rp),    intent(in)  :: epsi(3,3) 
  integer(ip), intent(in)  :: in_wall
  real(rp),    intent(inout) :: fsaff(3)          ! Saffman's force return as result
  real(rp),    intent(in)  :: wdist
  real(rp)                 :: Re_s(3)           ! Reynolds of the particule
  real(rp)                 :: Re_g(3)           ! Reynolds based on shear field            
  real(rp)                 :: Stlen(3)          ! Stokes lengthscale
  real(rp)                 :: Salen(3)          ! Saffman lengthscale
  real(rp)                 :: ka              
  real(rp)                 :: lambda(3)
!  real(rp),    intent(in)  :: wdist             ! particle distance to wall
  real(rp)                 :: velor(3)          ! absolute value of relative velocity of the fluid seen by the particle
  real(rp)                 :: e(3)
  real(rp)                 :: G(3)              ! gradient of fluid velocity
!-------------------------------------------------------------------------------------
!
! Definition des grandeurs utiles
!
!-------------------------------------------------------------------------------------
  G = 0.0_rp
  velor = 0.0_rp
  do idime = 1,ndime
    velor(idime) = abs( v_fl(idime) - v_pt(idime) )
    do jdime = 1,ndime
      G(idime) = G(idime) + epsi(idime,jdime)
    end do
    if (saff_denom /= 0.0) then
      G(idime) = G(idime) / saff_denom
    else
      G(idime) = 0.0_rp
    end if 
  end do
!write(2000,*) epsi,saff_denom
  do idime = 1,ndime
    Re_g(idime) = 0.0_rp
    Re_s(idime) = 0.0_rp
    Re_g(idime) = velor(idime) * d / visc_k
    Re_s(idime) = abs(G(idime)) * (d * d) / visc_k
  end do
 
  do idime = 1,ndime
    e(idime) = 0.0_rp
    if ( velor(idime) /= 0.0 ) then
      e(idime) = sqrt( abs(G(idime)) * visc_k ) / velor(idime)
!     write(1000,*) e(idime),G,velor
    end if
  end do
!write(1000,*) sqrt(Re_s),Re_g
!write(2000,*) e,sqrt(Re_s)-Re_g


  do idime = 1,ndime
    Stlen(idime) = 0.0_rp
    Salen(idime) = 0.0_rp
    lambda(idime) = 0.0_rp
    if ( G(idime) /= 0.0_rp) then
      Salen(idime) = sqrt(visc_k / abs(G(idime)))
      lambda(idime) = (d / 2.0_rp) * velor(idime) / G(idime)
    end if
    if (velor(idime) /= 0.0 ) then
      Stlen = visc_k / abs( velor(idime))
    end if
  end do

!  wdist = 0.0_rp
!  do inode = 1,pnode
!    ipoin = lnods(inode,ielem)
!    wdist = wdist + shapf(inode) * walld(ipoin)  
!  end do

  ka = 0.0_rp
  !if ( wdist /= 0.0) then
    ka = (d / 2.0_rp) / wdist
  !end if

!--------------------------------------------------------------------------------------!
!										       !
!  Treat particle normaly with general Saffman formula for particle Reynolds up to 3   !
!                                                                                      !
!--------------------------------------------------------------------------------------
  fsaff = 0.0_rp
  do idime = 1,ndime
     !fsaff(idime) = ( 3.0_rp / 32.0_rp ) * pi * dens_fl * (( d / 2.0_rp ) ** 2) * velor(idime) * ( 6 * velor(idime) - 11 * G(idime) * wdist )
     fsaff(idime) = (d / 2.0_rp) * visc_k * (abs(v_fl(idime)) - abs(v_pt(idime))) * Re_g(idime) * ( ( 1.7631_rp + (0.3561_rp * ka) - (1.1837_rp * (ka**2)) +&
                                 & (0.845163_rp * (ka**3)) ) - ( (3.24139_rp / ka) + 2.6760_rp + (0.8248_rp * ka) - (0.4616_rp*(ka**2)) )*lambda(idime) +&
                                 & (1.8081_rp + (0.879585 * ka) - (1.9009_rp * (ka**2)) + (0.98149_rp * (ka**3)) )* (lambda(idime)**2) )
  end do  

!write(6666,*) fsaff
!flush(6666)

end subroutine pts_saffma
