!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_elmcoe(&
     pnode,gpvis,gpden,gpsha,gpcar,elvel,eltur,eledd,&
     elwal,elust,eltem,elgr2,elsqk,elgrp,gpvol,gpdif,gpgrd,&
     gprea,gpvel,gprhs,gpgrv,gptur)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_elmcoe
  ! NAME
  !   tur_elmcoe
  ! DESCRIPTION
  !    Compute coefficient of the equation of Spalart-Almmaras model
  !    Useful relations:
  !    eps = Cmu*k^{3/2}/l
  !    w   = sqrt(k)/l = Cmu*eps/k
  !    nut = sqrt(k)*l = Cmu*k^2/eps = (alpha*)*k/w   
  ! OUTPUT 
  !    GPREA .......... r 
  !    GPDIF .......... k 
  !    GPRHS .......... f 
  !    GPGRD(NDIME) ... grad(k) coefficient
  !    GPVEL .......... Advection
  ! USES
  ! USED BY
  !    tur_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime
  use def_turbul, only     :  nturb_tur,iunkn_tur,kfl_timei_tur,&
       &                      kfl_cotem_tur,kfl_grve2_tur,&
       &                      param_tur,kfl_walld_tur,kfl_ustar_tur,&
       &                      grnor_tur,gravi_tur,boube_tur,prtur_tur,&
       &                      ipara_tur,kfl_inidi_tur,kfl_algor_tur,&
       &                      TUR_SPALART_ALLMARAS,TUR_K_XU_CHIEN,&            
       &                      TUR_K_EPS_LAUNDER_SHARMA,&   
       &                      TUR_K_EPS_CHIEN,TUR_K_EPS_LAM_BREMHORST,&   
       &                      TUR_TWO_LAYER_RODI,TUR_TWO_LAYER_XU_CHEN,&      
       &                      TUR_K_OMEGA,TUR_K_OMEGA_BREDBERG,&
       &                      TUR_K_EPS_V2_F,TUR_K_EPS_JAW_HWANG,&
       &                      TUR_K_EPS_NAGANO_TAGAWA,kfl_grsqk_tur,&
       &                      TUR_K_EPS_PHI_F, kfl_kxmod_tur, TUR_K_EPS_STD, kfl_kxmod_tur
  implicit none
  integer(ip), intent(in)  :: pnode
  real(rp),    intent(in)  :: gpvis,gpden
  real(rp),    intent(in)  :: gpsha(pnode),gpcar(ndime,pnode)
  real(rp),    intent(in)  :: elvel(ndime,pnode)
  real(rp),    intent(in)  :: eltur(nturb_tur,pnode,3)
  real(rp),    intent(in)  :: eledd(pnode),elwal(pnode)
  real(rp),    intent(in)  :: elust(pnode),eltem(pnode),gpvol
  real(rp),    intent(in)  :: elgr2(pnode),elsqk(pnode)
  real(rp),    intent(in)  :: elgrp(ndime,pnode)
  real(rp),    intent(out) :: gpdif(kfl_algor_tur)
  real(rp),    intent(out) :: gpgrd(kfl_algor_tur,ndime)
  real(rp),    intent(out) :: gprea(kfl_algor_tur,kfl_algor_tur)
  real(rp),    intent(out) :: gpvel(ndime)
  real(rp),    intent(out) :: gprhs(kfl_algor_tur)
  real(rp),    intent(out) :: gpgrv(ndime,ndime)
  real(rp),    intent(out) :: gptur(nturb_tur,3)
  integer(ip)              :: idime,inode,iturb,jturb,jdime
  real(rp)                 :: gpmut,gpgra,gppro,gpwal,gpgte(3)
  real(rp)                 :: gpust,gphev,gppr2,gpsqk,gpgrp(3),gplap, eta, cmu
  !
  ! Initialization
  !
  do iturb=1,kfl_algor_tur
     gprhs(iturb)=0.0_rp
     do jturb=1,kfl_algor_tur
        gprea(jturb,iturb)=0.0_rp
     end do
  end do
  !
  ! Turbulence variables GPTUR(NTURB_TUR,3)
  !
  do iturb=1,nturb_tur
     gptur(iturb,1)=0.0_rp
  end do
  do inode=1,pnode
     do iturb=1,nturb_tur
        gptur(iturb,1)=gptur(iturb,1)&
             +gpsha(inode)*eltur(iturb,inode,1)
     end do
  end do
  if(kfl_timei_tur/=0) then
     do iturb=1,nturb_tur
        gptur(iturb,3)=0.0_rp
     end do
     do inode=1,pnode
        do iturb=1,nturb_tur
           gptur(iturb,3)=gptur(iturb,3)&
                +gpsha(inode)*eltur(iturb,inode,3)
        end do
     end do
  end if
  !
  ! Eddy viscosity GPMUT = mu_t
  !
  gpmut=0.0_rp
  do inode=1,pnode
     gpmut=gpmut+gpsha(inode)*eledd(inode)
  end do
  !
  ! Velocity GPVEL(NDIME) = a
  !
  do idime=1,ndime
     gpvel(idime)=0.0_rp
  end do
  if(TUR_K_EPS_PHI_F.and.iunkn_tur==4) then ! do not need velocity
     continue                                   
  else
     if(kfl_inidi_tur==0) then
        do inode=1,pnode
           do idime=1,ndime
              gpvel(idime)=gpvel(idime)+gpsha(inode)*elvel(idime,inode)
           end do
        end do
     end if
  end if
  !
  ! Velocity gradient GPGRV(j,i) = grad(a) = da(i)/dx(j)
  !
  do idime=1,ndime
     do jdime=1,ndime
        gpgrv(jdime,idime)=0.0_rp
     end do
  end do
  do inode=1,pnode
     do idime=1,ndime
        do jdime=1,ndime
           gpgrv(jdime,idime)=gpgrv(jdime,idime)&
                +gpcar(jdime,inode)*elvel(idime,inode)
        end do
     end do
  end do
  !
  ! Wall distance GPWAL. Set to small value if all nodes of current
  ! element are located on the wall
  !
  if(kfl_walld_tur/=0) then
     gpwal=0.0_rp
     do inode=1,pnode
        gpwal=gpwal+gpsha(inode)*elwal(inode)
     end do
     if(gpwal==0.0_rp) then
        gpwal=1.0e-6_rp*gpvol**(1.0_rp/real(ndime,rp))
     end if
  end if
  !
  ! Wall shear stress GPUST = U*
  !
  if(kfl_ustar_tur/=0) then
     gpust=0.0_rp
     do inode=1,pnode
        gpust=gpust+gpsha(inode)*elust(inode)
     end do
  end if
  !
  ! Shear production GPPRO = mu_t*(dui/dxj+duj/dxi)*dui/dxj
  !
  if(.not.TUR_SPALART_ALLMARAS) then
     gppr2=0.0_rp
     do idime=1,ndime
        do jdime=1,ndime
           gppr2=gppr2+(gpgrv(jdime,idime)+gpgrv(idime,jdime))&
                *gpgrv(jdime,idime)
        end do
     end do
     gppro=gpmut*gppr2
  end if
  !
  ! Gravity production GPGRA=beta*mut/Prt*g.grad(T)
  !
  gpgra=0.0_rp
  if(kfl_cotem_tur/=0) then
     do idime=1,ndime
        gpgte(idime)=0.0_rp
     end do
     do inode=1,pnode
        do idime=1,ndime
           gpgte(idime)=gpgte(idime)+gpcar(idime,inode)*eltem(inode)
        end do
     end do
     do idime=1,ndime
        gpgra=gpgra+gpgte(idime)*grnor_tur*gravi_tur(idime)
     end do
     gpgra=boube_tur*gpgra/prtur_tur*gpmut
  end if
  !
  ! Velocity 2nd order derivative
  !
  if(kfl_grve2_tur/=0) then
     gphev=0.0_rp
     do inode=1,pnode
        gphev=gphev+elgr2(inode)*gpsha(inode)
     end do
  end if
  !
  ! Velocity 2nd order derivative
  !
  if(kfl_grsqk_tur/=0) then
     gpsqk=0.0_rp
     do inode=1,pnode
        gpsqk=gpsqk+elsqk(inode)*gpsha(inode)
     end do
  end if
  !
  ! Lapl(phi)
  !
  if(TUR_K_EPS_PHI_F.and.iunkn_tur==4) then
     gplap=0.0_rp
     gpgrp=0.0_rp
     do inode=1,pnode
        do idime=1,ndime
           gpgrp(idime)=gpgrp(idime)+gpcar(idime,inode)*elgrp(idime,inode)
        end do
     end do
     gplap=0.0_rp
     do idime=1,ndime
        gplap=gplap+gpgrp(idime)
     end do
  end if

  !------------------------------------------------------------------------
  !
  ! Compute coefficients for the different models
  !
  !------------------------------------------------------------------------

  if(TUR_SPALART_ALLMARAS) then
     ! 
     ! Spalart-Allmaras
     ! 
     call tur_spaalm(&
          ndime,pnode,nturb_tur,ipara_tur,param_tur,&
          eltur,gptur,gpvis,gpden,gpgrv,gpwal,gpcar,&
          gprea,gpdif,gprhs,gpgrd,gpvel)

  else if(TUR_K_XU_CHIEN) then
     !
     ! Xu-Chen-Nieuwstadt k model for natural convection
     !
     call tur_xuchen(& 
          ndime,pnode,param_tur,gptur,gpden,gpvis,&
          gpmut,gpwal,gpgra,eledd,gpcar,gppro,&
          gprea,gpdif,gprhs,gpgrd)

  else if(TUR_K_EPS_STD.or.TUR_TWO_LAYER_RODI.or.TUR_TWO_LAYER_XU_CHEN) then
     !
     ! Standard k-e model, two-layer Rodi, two-layer Xu-Chen
     ! for natural convection
     !
     call tur_stdkep(& 
          pnode,gptur,gpden,gpvis,gpmut,gpgra,eledd,& 
          gpcar,gppro,gppr2,gprea,gpdif,gprhs,gpgrd,kfl_kxmod_tur, eta, cmu)

  else if(TUR_K_EPS_LAUNDER_SHARMA) then
     !
     ! Launder-Sharma k-e model
     !
     call tur_lausha(& 
          pnode,gptur,gpden,gpvis,gpmut,gpgra,gphev,&
          eledd,gpcar,gppro,gpsqk,gprea,gpdif,gprhs,&
          gpgrd)

  else if(TUR_K_EPS_CHIEN) then
     !
     ! Chien k-e model
     !
     call tur_kchien(& 
          ndime,pnode,gptur,gpden,gpvis,gpmut,gpwal,&
          gpgra,eledd,gpust,gpcar,gppro,gppr2,gprea,&
          gpdif,gprhs,gpgrd)

  else if(TUR_K_EPS_LAM_BREMHORST) then
     !
     ! Lam-Bremhorst k-e model
     !
     call tur_lambre(& 
          ndime,pnode,gptur,gpden,gpvis,gpmut,gpwal,&
          gpgra,eledd,gpcar,gppro,gppr2,gprea,gpdif,&
          gprhs,gpgrd)

  else if(TUR_K_EPS_JAW_HWANG) then
     !
     ! Jaw-Hwang k-e model
     !
     call tur_kjawhw(&
          ndime,pnode,gptur,gpden,gpvis,gpmut,gpwal,&
          gpgra,eledd,gpust,gpcar,gppro,gppr2,gprea,&
          gpdif,gprhs,gpgrd)

  else if(TUR_K_EPS_NAGANO_TAGAWA) then
     !
     ! Nagano-Tagawa k-e model
     !
     call tur_nagano(&
          ndime,pnode,gptur,gpden,gpvis,gpmut,gpwal,&
          gpgra,eledd,gpust,gpcar,gppro,gppr2,gprea,&
          gpdif,gprhs,gpgrd)

  else if(TUR_K_EPS_V2_F) then
     !
     ! Lien-Durbin k-e-v2-f model
     !
     call tur_kepsv2(& 
          ndime,pnode,gptur,gpden,gpvis,gpmut,gpwal,&
          gpgra,eledd,eltur,gpcar,gppro,gppr2,gplap,&
          gprea,gpdif,gprhs,gpgrd)
     
  else if(TUR_K_OMEGA) then
     !
     ! Standard k-w
     !
     call tur_komega(&
          ndime,pnode,gptur,gpden,gpvis,gpmut,gpgra,eledd,&
          gpcar,gppro,gppr2,gprea,gpdif,gprhs,gpgrd)

  else if(TUR_K_OMEGA_BREDBERG) then
     !
     ! Bredberg-Peng-Davidson k-w
     !
     call tur_brepen(&
          ndime,pnode,gptur,gpden,gpvis,gpmut,gpgra,&
          eledd,eltur,gpcar,gppro,gppr2,gprea,gpdif,&
          gprhs,gpgrd)

  else if(TUR_K_EPS_PHI_F) then
     !
     ! Laurence-Uribe-Utyuzhnikov k-e-phi-f model
     !
     call tur_kepspf(& 
          ndime,pnode,gptur,gpden,gpvis,gpmut,gpwal,&
          gpgra,eledd,eltur,gpcar,gppro,gplap,gprea,&
          gpdif,gprhs,gpgrd)

  end if
  !
  ! Two layer models: compute inner layer coefficients
  !
  if(TUR_TWO_LAYER_RODI) then
     !
     ! Two-layer Rodi
     !
     call tur_tworod(& 
          ndime,pnode,nturb_tur,iunkn_tur,param_tur,&
          gptur,gpden,gpvis,gpmut,gpwal,gpgra,eledd,&
          gpust,gpsha,gpcar,gppro,gprea,gpdif,gprhs,&
          gpgrd,gpvel)

  else if(TUR_TWO_LAYER_XU_CHEN) then
     !
     ! Two-layer Xu-Chen for natural convection
     !
     call tur_twonat(& 
          ndime,pnode,nturb_tur,iunkn_tur,param_tur,&
          gptur,gpden,gpvis,gpmut,gpwal,gpgra,eledd,&
          gpcar,gppro,gprea,gpdif,gprhs,gpgrd,gpvel)
  end if

end subroutine tur_elmcoe
