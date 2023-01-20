!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_elmexf.f90
!> @author  Guillaume Houzeaux
!> @brief   Add a force term to the momentum equations
!> @details 
!> @} 
!-----------------------------------------------------------------------
subroutine nsi_elmexf(                     &
     ndime,pgaus,pnode,lforc_material_nsi,gpden, &
     xforc_material_nsi,gprhs, gpvel,elcod,gpsha,pmate)

  use def_kintyp,            only : ip,rp
  use def_parame,            only : pi
  use def_nastin,            only : ntabr_nsi

#ifndef NASTIN_PRIVATE_OFF
  use mod_nsi_actuator_disc, only : nsi_forcent
#endif

  implicit none
  integer(ip), intent(in)  :: ndime
  integer(ip), intent(in)  :: pgaus
  integer(ip), intent(in)  :: pnode
  integer(ip), intent(in)  :: lforc_material_nsi
  integer(ip), intent(in)  :: pmate
  real(rp),    intent(in)  :: gpden(pgaus)
  real(rp),    intent(in)  :: xforc_material_nsi(*), gpvel(ndime, pgaus)
  real(rp),    intent(in)  :: elcod(ndime,pnode) 
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(out) :: gprhs(ndime,pgaus)
  integer(ip)              :: idime,igaus,inode
  real(rp)                 :: n_disk(3),Ct,Cp,uref2,uref3,V_disk
  real(rp)                 :: d_disk, veinf, Force_normal, Force_tange
  real(rp)                 :: norma, tange, radiu, radia(ndime), tangd(ndime), rcent(ndime)
  real(rp)                 :: gpcod(ndime), proje, omega, norm

  ! rotational velocity
  ! from RPM to rad/s
  omega = 9.0_rp *pi/30.0_rp
  
  if( lforc_material_nsi == 1 ) then
     !
     ! Wake model
     !
     n_disk(1) = xforc_material_nsi( 3)
     n_disk(2) = xforc_material_nsi( 4)
     n_disk(3) = xforc_material_nsi( 5)
     V_disk    = xforc_material_nsi(15)
     d_disk    = xforc_material_nsi( 1)
     Veinf     = xforc_material_nsi(18) ! infinite velocity
     ! traction coefficient
     Ct        =  xforc_material_nsi(16)
     Cp        =  xforc_material_nsi(17)
     uref2     =  veinf*veinf
     uref3     =  veinf*veinf*veinf
     !  modification of Ct to use local velocity
     !  assuming theoretical relationship Ct=4a*(1-a)
     !     Ct         = 4.0_rp*Ct/(2.0_rp+2.0_rp*sqrt(1.0_rp-Ct)-Ct)
     if (ntabr_nsi(pmate)==0) then !Uniform loaded disc, non normal and tangential forces distibution
   
        do igaus = 1,pgaus
           !      
           ! local velocity !
           !
           !       uref2 = gpvel(1, igaus) * n_disk(1) &
           !            + gpvel(2, igaus) * n_disk(2) & 
           !            + gpvel(3, igaus) * n_disk(3) 
           !       uref2 = uref2*uref2
           
           Force_normal = 0.5_rp * gpden(igaus) * Ct * uref2/d_disk
           ! force assembly to rhs
           gprhs(1:ndime,igaus) = gprhs(1:ndime,igaus) + n_disk(1:ndime) * Force_normal
        end do
     else   ! normal and tangential forces distibution
        ! center of disc       
        rcent(1:3) = xforc_material_nsi( 6:8) 
        do igaus = 1,pgaus
           ! gauss point coordinates
           gpcod(1:ndime) =0.0_rp
           do inode =1, pnode
              gpcod(1:ndime) = gpcod(1:ndime) + gpsha(inode, igaus)*elcod(1:ndime, inode)
           end do
           ! radial coordinate (respect to center of disc)
           radia(1:ndime) = gpcod(1:ndime) - rcent(1:ndime)
           ! projection on to normal
           proje = dot_product(radia(1:ndime),n_disk(1:ndime))
           ! radial projection
           radia(1:ndime) = radia(1:ndime)- proje*n_disk(1:ndime) 
           ! radial coordinate
           radiu = sqrt(radia(1)*radia(1) +  radia(2)*radia(2)+ radia(3)*radia(3))
           ! normal and tangential forces
#ifndef NASTIN_PRIVATE_OFF
           call nsi_forcent(radiu, norma, tange, pmate)
#endif
           !normal foce 
           if (radiu > 1.0e-7_rp) then
              Force_normal = 0.5_rp* gpden(igaus) * Ct * uref2 * norma/(2.0_rp*pi*radiu*d_disk)
              !tangential force
              Force_tange  = 0.5_rp* gpden(igaus) * Cp * uref3 * tange/(2.0_rp*pi*radiu*d_disk*omega)
              ! tangent direction
              call vecpro(radia,n_disk,tangd,ndime)
              ! tangential direction unit vector
              norm = sqrt(tangd(1)*tangd(1) + tangd(2)*tangd(2)+tangd(3)*tangd(3))
              tangd(1:ndime)=tangd(1:ndime)/norm
           else
              Force_normal = 0.0_rp
              Force_tange  = 0.0_rp
              tangd= 0.0_rp
           end if
           
           ! tangential direction radia x norma
           call vecpro(radia,n_disk,tangd,ndime)
           ! tangential direction unit vector
           norm = sqrt(tangd(1)*tangd(1) + tangd(2)*tangd(2)+tangd(3)*tangd(3))
           tangd(1:ndime)=tangd(1:ndime)/norm
           ! rhside total force
           gprhs(1:ndime,igaus) = gprhs(1:ndime,igaus) + n_disk(1:ndime) * Force_normal - tangd(1:ndime)*Force_tange
           
        end do
        
     end if
  else if ( lforc_material_nsi == 2) then
     do igaus = 1,pgaus
        do idime=1, ndime
           gprhs(idime, igaus) =  gprhs(idime, igaus) +  xforc_material_nsi(idime)
        end do
     end do
  end if

!!$  if (kfl_tendencies_ker) then
!!$     fcori = 2.0_rp*abs(fvela_nsi(3))
!!$      do igaus = 1,pgaus
!!$         gprhs(1, igaus) =  gprhs(1, igaus) + gpden(igaus)*fcori*(-gpten_geo(2, igaus) + gpten_adv(1,igaus))
!!$         gprhs(2, igaus) =  gprhs(2, igaus) + gpden(igaus)*fcori*( gpten_geo(1, igaus) + gpten_adv(2,igaus))
!!$      end do
!!$  end if
  
end subroutine nsi_elmexf

