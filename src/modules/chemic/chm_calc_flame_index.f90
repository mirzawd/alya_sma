!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_calc_flame_index(&
     pnode,pgaus,elconZ,elconYc,gpcar,&
     flame_index,gp_zgrad,gp_z,gp_c) 
  
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_calc_flame_index
  ! NAME 
  !    chm_calc_flame_index
  ! DESCRIPTION
  !    Compute Flame index at Gauss points for flamelet combustion model
  ! USES
  ! USED BY
  !    chm_post_scalar_dissipation_rate_fast
  !***
  !-----------------------------------------------------------------------
  use def_master, only      :  zeror
  use def_kintyp, only      :  ip,rp
  use def_domain, only      :  ndime,mnode
  use def_chemic, only      :  chm_zmax,chm_zmin

  implicit none
  integer(ip),  intent(in)    :: pnode,pgaus
  real(rp),     intent(in)    :: elconZ(pnode)
  real(rp),     intent(in)    :: elconYc(pnode)
  real(rp),     intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),     intent(in)    :: gp_zgrad(pgaus)
  real(rp),     intent(in)    :: gp_z(pgaus)
  real(rp),     intent(in)    :: gp_c(pgaus)
  real(rp),     intent(inout) :: flame_index(pgaus)
  real(rp)                    :: grad_Z(pgaus,ndime),grad_Yc_grad_Z(pgaus)
  real(rp)                    :: zeta(pgaus),maggradYC(pgaus)
  real(rp)                    :: grad_Yc(pgaus,ndime)
  real(rp)                    :: zmid
  integer(ip)                 :: igaus,inode,idime
  
  zeta = 0.0_rp
  grad_Yc_grad_Z = 0.0_rp
  grad_Yc     = 0.0_rp
  grad_Z     = 0.0_rp

  do igaus = 1,pgaus
     do inode = 1,pnode
        do idime = 1,ndime
           grad_Yc(igaus,idime) = grad_Yc(igaus,idime) + gpcar(idime,inode,igaus) * elconZ(inode)
           grad_Z(igaus,idime) = grad_Z(igaus,idime) + gpcar(idime,inode,igaus) * elconYc(inode)
        end do
     end do
  end do

  !
  ! Dani's Implementation
  !
  do igaus = 1,pgaus
     !
     ! zeta = dot(gradZ.gradYC)
     !
     if (ndime == 2_ip) then
        zeta(igaus) = dot_product(grad_Z(igaus,1:2),grad_Yc(igaus,1:2))
        zeta(igaus) = sqrt(zeta(igaus)* zeta(igaus))
        maggradYC(igaus) = sqrt( grad_Yc(igaus,1) **2 + grad_Yc(igaus,2) **2 )
     else 
        zeta(igaus) = dot_product(grad_Z(igaus,1:3),grad_Yc(igaus,1:3))
        zeta(igaus) = sqrt(zeta(igaus)* zeta(igaus))
        maggradYC(igaus) = sqrt( grad_Yc(igaus,1) **2 + grad_Yc(igaus,2) **2 + grad_Yc(igaus,3) **2)
     end if 

     if (gp_c(igaus) < 0.00001_rp ) then 
       zeta(igaus) = 0.0_rp
     else if ( maggradYC(igaus) < 0.00001_rp ) then
       zeta(igaus) = max(0.0_rp, ( 1.0_rp - zeta(igaus) / (gp_zgrad(igaus) + zeror) ) )
     else 
       zeta(igaus) = max(0.0_rp, (1.0_rp - (zeta(igaus) / (gp_zgrad(igaus) * maggradYC(igaus)))))
     end if  
     zmid = chm_zmax / ( 1.0_rp + chm_zmax - chm_zmin) 
     if (gp_z(igaus) < zmid) then
        zeta(igaus) = zeta(igaus) * min(1.0_rp,(gp_z(igaus)/chm_zmin))
     else 
        zeta(igaus) = zeta(igaus) * min(1.0_rp, ( 1.0_rp-gp_z(igaus)/(1.0_rp - chm_zmax)))
     end if 

     flame_index(igaus) = max(0.0_rp,zeta(igaus))
  end do
  !
  ! My Implementation
  !
  !!do igaus = 1,pgaus
  !!   !
  !!   ! Normal to the YC gradient
  !!   !
  !!   if (ndime == 2_ip) then
  !!     norm_grad_Yc(igaus,:) = grad_Yc(igaus,:) / sqrt( grad_Yc(igaus,1) **2 + grad_Yc(igaus,2) **2 )
  !!     grad_Yc_grad_Z(igaus) =  dot_product(norm_grad_Yc(igaus,1:2),grad_Z(igaus,1:2))
  !!   else 
  !!     norm_grad_Yc(igaus,:) = grad_Yc(igaus,:) / sqrt( grad_Yc(igaus,1) **2 + grad_Yc(igaus,2) **2  + grad_Yc(igaus,3) **2)
  !!     grad_Yc_grad_Z(igaus) =  dot_product(norm_grad_Yc(igaus,1:3),grad_Z(igaus,1:3))
  !!   end if
  !!   
  !!   if ( ( gp_z(igaus) > chm_zmax ) .or. (gp_z(igaus) < chm_zmin ) )then
  !!      if  (gp_z(igaus) > chm_zmax) then
  !!         !k = 1.0_rp - ( (gp_z(igaus) - chm_zmax ) / ( 1.0_rp - chm_zmax ) )
  !!         k = min(1.0_rp, (1.0_rp - gp_z(igaus)/( 1.0_rp - chm_zmax  )  ) ) 
  !!      elseif (gp_z(igaus) < chm_zmin ) then
  !!         k = min( 1.0_rp, (gp_z(igaus) / chm_zmin) )
  !!      end if 
  !!   else
  !!      k = 1.0_rp
  !!   endif 
  !!   flame_index(igaus) = k * ( 1.0_rp - ( grad_Yc_grad_Z(igaus) / abs(gp_zgrad(igaus) + zeror) ) )
  !!end do
end subroutine chm_calc_flame_index
