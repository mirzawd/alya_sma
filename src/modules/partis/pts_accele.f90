!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    pts_accele.f90
!> @author  Guillaume Houzeaux
!> @brief   Read physical data
!> @details Read physical data
!> @} 
!-----------------------------------------------------------------------
subroutine pts_accele(                                           &
  ndime,kfl_drafo,kfl_extfo,grafo,buofo,gravi,veloc_fluid,       &
  visco_fluid,densi_fluid,veloc_parti,coord_parti,densi_parti,   &
  diame_parti,spher_parti,dista,t,accel)

  use def_kintyp,  only : ip,rp
  use mod_physics, only : physics_drag_force 
  implicit none
  integer(ip), intent(in)  :: ndime
  integer(ip), intent(in)  :: kfl_drafo
  integer(ip), intent(in)  :: kfl_extfo
  real(rp),    intent(in)  :: grafo
  real(rp),    intent(in)  :: buofo
  real(rp),    intent(in)  :: gravi(*)
  real(rp),    intent(in)  :: veloc_fluid(*)
  real(rp),    intent(in)  :: visco_fluid
  real(rp),    intent(in)  :: densi_fluid
  real(rp),    intent(in)  :: veloc_parti(*)
  real(rp),    intent(in)  :: coord_parti(*)
  real(rp),    intent(in)  :: densi_parti
  real(rp),    intent(in)  :: diame_parti
  real(rp),    intent(in)  :: spher_parti
  real(rp),    intent(in)  :: dista
  real(rp),    intent(in)  :: t
  real(rp),    intent(out) :: accel(*)
  integer(ip)              :: idime
  real(rp)                 :: CdRe,Cd,Re,dCddRe,alpha,Du
  real(rp)                 :: acext(3)

  do idime = 1,ndime
     accel(idime) = 0.0_rp
  end do
  !
  ! Drag force
  !
  if( kfl_drafo /= 0 ) then

     Du = 0.0_rp
     do idime = 1,ndime
        Du =  Du + ( veloc_parti(idime) - veloc_fluid(idime) ) ** 2
     end do
     Du = sqrt( Du )

     call physics_drag_force(kfl_drafo,Du,visco_fluid,densi_fluid,diame_parti,Cd,Re,CdRe,dCddRe,spher_parti,dista)                      

     alpha = - 0.75_rp * visco_fluid / ( densi_parti * diame_parti * diame_parti ) * CdRe 

     do idime = 1,ndime
        accel(idime) =  accel(idime) + alpha * ( veloc_parti(idime) - veloc_fluid(idime) ) 
     end do

  end if
  !
  ! External force 
  !
  if( kfl_extfo /= 0 ) then
     call extefo(& 
          kfl_extfo,coord_parti,densi_parti,spher_parti,densi_fluid,visco_fluid,t,acext)
     do idime = 1,ndime
        accel(idime) =  accel(idime) + acext(idime)
     end do
  end if
  !
  ! Gravity and buoyancy
  !
  alpha = ( grafo - buofo * densi_fluid / densi_parti )
  do idime = 1,ndime
     accel(idime) = accel(idime) + gravi(idime) * alpha
  end do
  
end subroutine pts_accele
