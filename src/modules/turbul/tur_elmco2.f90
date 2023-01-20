!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_elmco2(&
     pnode,gpvis,gpden,gpsha,gpcar,elvel,elmsh,eltur,eledd,&
     elwal,elust,eltem,elgr2,elsqk,elgrp,elpro,elprd,gpvol,&
     chale,gpdif,gpgrd,gprea,gpvel,gprhs,gpgrv,gptur,&
     gpprj,hleng,fddes,gddes,sstf1,sasso, gpcan, sreac, gpprr, &
     elprr, gppgr, elpgr, pmate, turvi, ipass, conve, gptvi, ellmax)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_elmco2
  ! NAME
  !   tur_elmco2
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
  use def_master, only     :  kfl_coupl,ID_TURBUL,ID_ALEFOR
  use def_domain, only     :  ndime
  use def_kermod, only     :  cmu_st, kfl_logva
  use def_turbul, only     :  nturb_tur,iunkn_tur,kfl_timei_tur,&
       &                      kfl_cotem_tur,kfl_grve2_tur,&
       &                      param_tur,kfl_walld_tur,kfl_ustar_tur,&
       &                      grnor_tur,gravi_tur,boube_tur,prtur_tur,&
       &                      ipara_tur,kfl_inidi_tur,&
       &                      kfl_ortho_tur,kfl_produ_tur,&
       &                      TUR_SPALART_ALLMARAS,TUR_K_XU_CHIEN,&            
       &                      TUR_K_EPS_LAUNDER_SHARMA,&   
       &                      TUR_K_EPS_CHIEN,TUR_K_EPS_LAM_BREMHORST,&   
       &                      TUR_TWO_LAYER_RODI,TUR_TWO_LAYER_XU_CHEN,&      
       &                      TUR_K_OMEGA,TUR_K_OMEGA_BREDBERG,&
       &                      TUR_K_EPS_V2_F,TUR_K_EPS_JAW_HWANG,&
       &                      TUR_K_EPS_NAGANO_TAGAWA,kfl_grsqk_tur,&
       &                      TUR_K_EPS_PHI_F,TUR_SST_K_OMEGA, &
       &                      TUR_K_EPS_STD, TUR_TKE_SGS,  &
       &                      zetur, kfl_kxmod_tur, &
       &                      kfl_shock_tur, nbdfp_tur
  use mod_ker_regularization, only : kfl_regularization, regul_k, regul_e
  implicit none
  integer(ip), intent(in)  :: pnode
  real(rp),    intent(in)  :: gpvis,gpden
  real(rp),    intent(in)  :: gpsha(pnode),gpcar(ndime,pnode)
  real(rp),    intent(in)  :: elvel(ndime,pnode)
  real(rp),    intent(in)  :: elmsh(ndime,pnode)
  real(rp),    intent(in)  :: eltur(nturb_tur,pnode,nbdfp_tur+1)
  real(rp),    intent(in)  :: eledd(pnode),elwal(pnode)
  real(rp),    intent(in)  :: elust(pnode),eltem(pnode)
  real(rp),    intent(in)  :: gpvol,chale(2)
  real(rp),    intent(in)  :: elgr2(pnode),elsqk(pnode)
  real(rp),    intent(in)  :: elgrp(ndime,pnode)
  real(rp),    intent(in)  :: elpro(pnode)
  real(rp),    intent(in)  :: elprd(pnode)
  real(rp),    intent(in)  :: elprr(pnode), elpgr(ndime, pnode), ellmax(pnode)
  real(rp),    intent(out) :: fddes, gddes
  real(rp),    intent(out) :: sstf1, sasso
  real(rp),    intent(out) :: gpdif
  real(rp),    intent(inout) :: gpgrd(ndime)
  real(rp),    intent(out) :: gprea
  real(rp),    intent(out) :: sreac
  real(rp),    intent(out) :: gpvel(ndime), conve(ndime) ! velocity and convective terms
  real(rp),    intent(out) :: gprhs
  real(rp),    intent(out) :: gpgrv(ndime,ndime)
  real(rp),    intent(out) :: gptur(nturb_tur,nbdfp_tur+1)
  real(rp),    intent(out) :: gpprj, gpprr, gppgr(ndime)
  real(rp),    intent(in)  :: hleng(ndime)
  real(rp),    intent(in)  :: gpcan
  real(rp),    intent(in)  :: gptvi
  integer(ip), intent(in)  :: pmate, ipass  
  integer(ip)              :: idime,inode,iturb,jdime, kdime, itime
  real(rp)                 :: gpmut(2),gpgra,gppro,gpwal,gpgte(3), phi, As, cmu,uaste
  real(rp)                 :: gpust,gphev,gppr2,gpsqk,gpgrp(3),gplap, eta, seci4, A0, W, divve, simgr(3,3), grunk(ndime, 2)
  real(rp)                 :: a1, F2, Cr, f0, f02, Wsqr6, turvi(2)
  real(rp)                 :: reguk, regue, sigmr, fpfac, gplmax
  !
  ! Initialization
  !
  gprhs = 0.0_rp
  gprea = 0.0_rp
  !
  ! Turbulence variables GPTUR
  !
  do iturb = 1,nturb_tur
     gptur(iturb,1) = 0.0_rp
     gptur(iturb,2) = 0.0_rp
  end do 
  do inode = 1,pnode
     do iturb = 1,nturb_tur
        gptur(iturb,1) = gptur(iturb,1) &
             + gpsha(inode) * eltur(iturb,inode,1)
        gptur(iturb,2) = gptur(iturb,2) &
             + gpsha(inode) * eltur(iturb,inode,2)
     end do
  end do
  !
  ! Time integration
  !
  if( kfl_timei_tur /= 0_ip ) then
     do itime = 3_ip,nbdfp_tur+1_ip
        do iturb = 1,nturb_tur
           gptur(iturb,itime) = 0.0_rp
        end do
        do inode = 1,pnode
           do iturb = 1,nturb_tur
              gptur(iturb,itime) = gptur(iturb,itime) &
                   + gpsha(inode) * eltur(iturb,inode,itime)
           end do
        end do
     end do
  end if
  !
  ! Velocity gradient GPGRV(j,i) = grad(a) = da(i)/dx(j)
  !
  gpgrv = 0.0_rp
  do inode = 1,pnode
     do idime = 1,ndime
        do jdime = 1,ndime
           gpgrv(jdime,idime) = gpgrv(jdime,idime) &
                + gpcar(jdime,inode) * elvel(idime,inode)
        end do
     end do
  end do
  divve =0.0_rp

  do idime = 1,ndime
     divve = divve + gpgrv(idime, idime)
  end do

  if( TUR_K_EPS_STD) then
     if (kfl_kxmod_tur==2.or.kfl_kxmod_tur==3) then ! Realizable or kefp model
        A0 = 4.04_rp 
        A0 = (1.0_rp - 3.0_rp*sqrt(0.5_rp*cmu_st))/cmu_st
        do idime = 1,ndime
!           divve = divve + gpgrv(idime, idime)
           do jdime = 1,ndime
              simgr(idime,jdime) = 0.5_rp*(gpgrv(idime, jdime) &
                   +  gpgrv(jdime, idime))
           end do
        end do
        seci4 =0.0_rp
        W =0.0_rp
        uaste=0.0_rp
        do idime = 1,ndime
           simgr(idime, idime ) = simgr(idime, idime) - divve/3.0_rp
        end do
        do idime = 1,ndime
           do jdime = 1,ndime
              uaste = uaste +  gpgrv(jdime,idime) &   !  D_ij : D_ij
                   &          * gpgrv(jdime,idime)
              seci4 = seci4 +  simgr (jdime,idime) &   !  S_ij : S_ij
                   &          * (simgr(jdime,idime))
              do kdime =1, ndime
                 W = W +  simgr(jdime,idime)* &
                      simgr(jdime,kdime)* &
                      simgr(kdime,idime)
              end do
           end do
        end do
        uaste = sqrt(uaste -divve*divve /3.0_rp)
        if (seci4.gt.zetur) W = W/sqrt(seci4*seci4*seci4)
        Wsqr6 = sqrt(6.0_rp)*W        
        Wsqr6 = min( Wsqr6, 1.0_rp)
        Wsqr6 = max( Wsqr6, -1.0_rp)
        phi = 1.0_rp/3.0_rp*acos( Wsqr6)
        As= sqrt(6.0_rp)*cos(phi)
        if (kfl_kxmod_tur==2) then      ! realizable
           cmu =1.0_rp/(A0+As*gptur(1,2)/gptur(2,2)*uaste)
        else if (kfl_kxmod_tur==3) then ! kefp model
           Cr = 4.5_rp
           f0 = Cr/(Cr-1.0_rp)
           f02= 4.0_rp*f0*(f0-1.0_rp)
           if (kfl_regularization) then
              reguk  = regul_k(gptur(1,2))
              regue  = regul_e(gptur(2,2))
           else
              reguk  = gptur(1,2)
              regue  = gptur(2,2)
           end if
           sigmr = cmu_st*(uaste*reguk/regue)* (uaste*reguk/regue)
           fpfac = 2.0_rp*f0 / (1.0_rp + sqrt(1.0_rp+ f02*sigmr ))
           cmu= cmu_st*fpfac
        end if

        gpmut(2) = max(gptur(1,2)*gptur(1,2)*cmu*gpden/gptur(2,2), gpvis)
        gpmut(1) = max(gptur(1,1)*gptur(1,1)*cmu*gpden/gptur(2,1), gpvis)
        
     else if (kfl_logva==1) then
        gpmut(2) = gpden*param_tur(6)*exp(2.0_rp*gptur(1,2)- gptur(2,2))
        gpmut(1) = gpden*param_tur(6)*exp(2.0_rp*gptur(1,1)- gptur(2,1))
       
        gpmut(2) = gpmut(2) + gpvis ! max(gpmut(2), gpvis)
        gpmut(1) = gpmut(1) + gpvis ! max(gpmut(1), gpvis)
        ! when eps tends to zero, 
        gpmut(2) = min(gpmut(2), 1.0e10_rp)
        gpmut(1) = min(gpmut(1), 1.0e10_rp)
        if (ipass==0) then ! initialization of global variable (turvi(iiter, igaus, ielem))
           turvi(2) = gpmut(2)
           turvi(1) = turvi(2)
        end if
     else
        ! standard k eps
        gpmut(2) = gptur(1,2)*gptur(1,2)*param_tur(6)*gpden/gptur(2,2)
        gpmut(1) = gptur(1,1)*gptur(1,1)*param_tur(6)*gpden/gptur(2,1)
        gpmut(2) = max(gpmut(2), gpvis)
        gpmut(1) = max(gpmut(1), gpvis)
        ! when eps tends to zero, 
        gpmut(1) = min(gpmut(1), 1.0e20_rp)
        gpmut(2) = min(gpmut(2), 1.0e20_rp)
     end if

  else if (TUR_SST_K_OMEGA) then ! K omega with open integration
     gpmut(1) =0.0_rp
     gpmut(2) = gpmut(1)
     gpwal = 0.0_rp
     do inode = 1,pnode
        gpwal = gpwal + gpsha(inode) * elwal(inode)
     end do
     if (abs(gpwal).gt.1.0e-8_rp) then
        do idime = 1,ndime
           do jdime = 1,ndime
!               simgr(idime,jdime) = 0.5_rp*(gpgrv(idime, jdime) &
!                    +  gpgrv(jdime, idime))
              !!!! CORRECTION
              simgr(idime,jdime) = 0.5_rp*(gpgrv(idime, jdime) &
                   - gpgrv(jdime, idime))
           end do
        end do
        seci4 =0.0_rp
        do idime = 1,ndime
           do jdime = 1,ndime
              seci4 = seci4 + 2.0_rp * simgr (jdime,idime) &   ! 2 W_ij : W_ij
                   &          * (simgr(jdime,idime))
           end do
        end do
        seci4 =sqrt(seci4)
        a1=0.31_rp       
        F2 = max( 2.0_rp *sqrt(max(0.0_rp,gptur(1,1))) /(0.09_rp*gptur(2,1)*gpwal), 500.0_rp * gpvis / (gpden*gpwal*gpwal*gptur(2,1))) 
        F2 = tanh( F2 * F2 )
        as = max( a1 * gptur(2,1),  seci4 * F2 )         
        gpmut(1) = gpden*a1*gptur(1,1)/as
        gpmut(2) = gpmut(1)
        gpmut(2) = max(gpmut(2), 1.0e-1_rp*gpvis)
        gpmut(1) = max(gpmut(1), 1.0e-1_rp*gpvis)
        gpmut(2) = min(gpmut(2), 1.0e20_rp)
        gpmut(1) = min(gpmut(1), 1.0e20_rp)
     end if

     ! not SST komega not k eps  
  else   ! properties from ker_proper                                
     gpmut(1) = gptvi
     gpmut(2) = gptvi
     gpgrd(1:ndime)= gpgrd(1:ndime)/param_tur(iunkn_tur)
  end if

  
!!$     !
  ! Diffusion gradient GPGRD = grad(mu_t)/sig
  !  
!!$     do idime = 1,ndime
!!$        gpgrd(idime) = 0.0_rp
!!$     end do
!!$     do inode = 1,pnode
!!$        do idime = 1,ndime
!!$           gpgrd(idime) = gpgrd(idime) + gpcar(idime,inode) * eledd(inode)
!!$        end do
!!$     end do
!!$     fact = 1/param_tur(iunkn_tur)
!!$     do idime = 1,ndime 
!!$        gpgrd(idime) = gpgrd(idime) * fact
!!$     end do
  !  zero viscosity gradient: (because it deteriorates inner linearization)
  gpgrd = 0.0_rp   
     
  !
  ! Velocity GPVEL(NDIME) = a
  !
  do idime = 1,ndime
     gpvel(idime) = 0.0_rp
  end do
  if( TUR_K_EPS_PHI_F .and. iunkn_tur == 4 ) then ! do not need velocity
     continue                                   
  else
     if( kfl_inidi_tur == 0 ) then
        if( kfl_coupl(ID_TURBUL,ID_ALEFOR) /= 0 ) then    ! ale : substract mesh velocity
           do inode = 1,pnode
              do idime = 1,ndime
                 gpvel(idime) = gpvel(idime) + gpsha(inode) * ( elvel(idime,inode) - elmsh(idime,inode) )
              end do
           end do
        else
           do inode = 1,pnode
              do idime = 1,ndime
                 gpvel(idime) = gpvel(idime) + gpsha(inode) * elvel(idime,inode)
              end do
           end do
        end if
     end if
  end if
  ! convective term equals velocity
  conve(1:ndime) = gpvel(1:ndime) 

  !
  ! Wall distance GPWAL. Set to small value if all nodes of current
  ! element are located on the wall
  !
  if( kfl_walld_tur /= 0.or.TUR_TKE_SGS ) then
     gpwal = 0.0_rp
     do inode = 1,pnode
        gpwal = gpwal + gpsha(inode) * elwal(inode)
     end do
     if( abs(gpwal).lt.epsilon(1.0_rp) ) then 
        gpwal = 1.0e-6_rp *gpvol ** (1.0_rp/real(ndime,rp))
     end if
  end if
  !
  ! Wall shear stress GPUST = U*
  !
  if( kfl_ustar_tur /= 0 ) then
     gpust = 0.0_rp
     do inode = 1,pnode
        gpust = gpust + gpsha(inode) * elust(inode)
     end do
  end if
  !
  ! Shear production GPPRO = mu_t*(dui/dxj+duj/dxi)*dui/dxj
  !
  if( .not. TUR_SPALART_ALLMARAS ) then
     if( kfl_produ_tur == 0 ) then
        gppr2 = 0.0_rp
        do idime = 1,ndime
           do jdime = 1,ndime
              gppr2 = gppr2 + ( gpgrv(jdime,idime) + gpgrv(idime,jdime) ) &
                   &          * gpgrv(jdime,idime) ! 2 S_ij : S_ij
           end do
        end do
        gppr2 = gppr2 - 2.0_rp*divve*divve/3.0_rp
        gppro = gpmut(2) * gppr2
        if (TUR_K_EPS_STD) &
             eta = sqrt(gppr2)*gptur(1,2)/gptur(2,2) ! parameter for RNG model 
     else
        !
        ! GPPRO = 2*mut*Sij*Sij
        !
        gppr2 = 0.0_rp
        do inode = 1,pnode
           gppr2 = gppr2 + gpsha(inode) * elprd(inode)
        end do
        gppro = 2.0_rp * gpmut(2) * gppr2
     end if
  end if
  !
  ! Gravity production GPGRA=beta*mut/Prt*g.grad(T)
  !
  gpgra = 0.0_rp
  if( kfl_cotem_tur /= 0_ip ) then 
     do idime = 1,ndime
        gpgte(idime) = 0.0_rp
     end do
     do inode = 1,pnode
        do idime = 1,ndime
           gpgte(idime) = gpgte(idime) + gpcar(idime,inode) * eltem(inode)
        end do
     end do
     do idime = 1,ndime
        gpgra = gpgra + gpgte(idime) * grnor_tur * gravi_tur(idime)
     end do
     gpgra = boube_tur * gpgra / prtur_tur * gpmut(2)
  end if
  !
  ! Velocity 2nd order derivative
  !
  if( kfl_grve2_tur /= 0_ip ) then
     gphev = 0.0_rp
     do inode = 1,pnode
        gphev = gphev + elgr2(inode) * gpsha(inode)
     end do
  end if
  !
  ! Velocity 2nd order derivative
  !
  if( kfl_grsqk_tur /= 0_ip ) then
     gpsqk = 0.0_rp
     do inode = 1,pnode
        gpsqk = gpsqk + elsqk(inode) * gpsha(inode)
     end do
  end if
  !
  ! Lapl(phi)
  !
  if( TUR_K_EPS_PHI_F .and. iunkn_tur == 4_ip ) then
     gplap = 0.0_rp
     gpgrp = 0.0_rp
     do inode = 1,pnode
        do idime = 1,ndime
           gpgrp(idime) = gpgrp(idime) + gpcar(idime,inode) * elgrp(idime,inode)
        end do
     end do
     gplap = 0.0_rp
     do idime = 1,ndime
        gplap = gplap + gpgrp(idime)
     end do
  end if
  !
  ! Orthogonal SGS
  !
  if( kfl_ortho_tur == 1_ip ) then ! FULL OSS
     gpprj = 0.0_rp
     do inode = 1,pnode
        gpprj = gpprj + gpsha(inode) * elpro(inode)
     end do
  else if ( kfl_ortho_tur == 2_ip ) then ! SPLIT  OSS
     gpprj = 0.0_rp
     gpprr = 0.0_rp
     do inode = 1,pnode
        gpprj = gpprj + gpsha(inode) * elpro(inode)
        gpprr = gpprr + gpsha(inode) * elprr(inode)
     end do
  end if
  !
  ! Orthogonal Shock capturing
  !
  if (kfl_shock_tur.ne.0) then
     ! projection of turbulence unknown  gradient
     gppgr(1:ndime) = 0.0_rp
     do inode =1, pnode
        gppgr(1:ndime) = gppgr(1:ndime) + gpsha(inode)*elpgr(1:ndime, inode)
     end do
  end if
  !------------------------------------------------------------------------
  !
  ! Compute coefficients for the different models
  !
  !------------------------------------------------------------------------
 
  if( TUR_SPALART_ALLMARAS ) then
     ! 
     ! Spalart-Allmaras
     ! 
     call tur_spaalm(&
          ndime,pnode,nturb_tur,ipara_tur,param_tur,&
          eltur,gptur,gpvis,gpden,gpgrv,gpwal,gpcar,&
          gprea,gpdif,gprhs,gpvel)

  else if( TUR_K_XU_CHIEN ) then
     !
     ! Xu-Chen-Nieuwstadt k model for natural convection
     !
     call tur_xuchen(& 
          ndime,pnode,param_tur,gptur,gpden,gpvis,&
          gpmut,gpwal,gpgra,eledd,gpcar,gppro,&
          gprea,gpdif,gprhs,gpgrd)

  else if( TUR_K_EPS_STD .or. TUR_TWO_LAYER_RODI .or. TUR_TWO_LAYER_XU_CHEN ) then
     !
     ! Standard k-e model, two-layer Rodi, two-layer Xu-Chen
     ! for natural convection
     !
     if (kfl_logva==1.or.kfl_regularization ) then ! unknown is logarithmic of k and log eps
        do idime =1, ndime
           grunk(idime, 1:2) = 0.0_rp
           do inode =1, pnode
              grunk(idime,1)=grunk(idime,1) +  eltur(iunkn_tur,inode,1)*gpcar(idime,inode)
              grunk(idime,2)=grunk(idime,2) +  eltur(iunkn_tur,inode,2)*gpcar(idime,inode)
           end do
        end do
     end if
     !
     ! maximum mixing length 
     !
     gplmax=0.0_rp
     do inode = 1,pnode
        gplmax = gplmax &
             + gpsha(inode) * ellmax(inode)
     end do
     gpgra=gpgra/gpmut(2)
   
     call tur_stdkep(& 
          pnode,gptur,gpden,gpvis,gpmut,gpgra,&
          gpcar,gppr2,gprea,gpdif,gprhs,&
          kfl_kxmod_tur, eta, cmu, gpcan, sreac, &
          pmate, grunk, conve, gplmax)

     gpvel(1:ndime)= conve(1:ndime) ! to stabilize convective terms
  else if( TUR_K_EPS_LAUNDER_SHARMA ) then
     !
     ! Launder-Sharma k-e model
     !
     call tur_lausha(& 
          pnode,gptur,gpden,gpvis,gpmut,gpgra,gphev,&
          eledd,gpcar,gppro,gpsqk,gprea,gpdif,gprhs,&
          gpgrd)

  else if( TUR_K_EPS_CHIEN ) then
     !
     ! Chien k-e model
     !
     call tur_kchien(& 
          ndime,pnode,gptur,gpden,gpvis,gpmut,gpwal,&
          gpgra,eledd,gpust,gpcar,gppro,gppr2,gprea,&
          gpdif,gprhs,gpgrd)

  else if( TUR_K_EPS_LAM_BREMHORST ) then
     !
     ! Lam-Bremhorst k-e model
     !
     call tur_lambre(& 
          ndime,pnode,gptur,gpden,gpvis,gpmut,gpwal,&
          gpgra,eledd,gpcar,gppro,gppr2,gprea,gpdif,&
          gprhs,gpgrd)

  else if( TUR_K_EPS_JAW_HWANG ) then
     !
     ! Jaw-Hwang k-e model
     !
     call tur_kjawhw(&
          ndime,pnode,gptur,gpden,gpvis,gpmut,gpwal,&
          gpgra,eledd,gpust,gpcar,gppro,gppr2,gprea,&
          gpdif,gprhs,gpgrd)

  else if( TUR_K_EPS_NAGANO_TAGAWA ) then
     !
     ! Nagano-Tagawa k-e model
     !
     call tur_nagano(&
          ndime,pnode,gptur,gpden,gpvis,gpmut,gpwal,&
          gpgra,eledd,gpust,gpcar,gppro,gppr2,gprea,&
          gpdif,gprhs,gpgrd)

  else if( TUR_K_EPS_V2_F ) then
     !
     ! Lien-Durbin k-e-v2-f model
     !
     call tur_kepsv2(& 
          ndime,pnode,gptur,gpden,gpvis,gpmut,gpwal,&
          gpgra,eledd,eltur,gpcar,gppro,gppr2,gplap,&
          gprea,gpdif,gprhs,gpgrd)
     
  else if( TUR_K_OMEGA ) then
     !
     ! Standard k-w
     !
     if (ipara_tur(1)==0) then ! old version
        call tur_komega(&
             ndime,pnode,gptur,gpden,gpvis,gpmut,gpgra,eledd,&
             gpcar,gppro,gppr2,gpgrv,eltur,gprea,gpdif,gprhs,&
             gpgrd)
     else ! new version k omega Wilcox
        call tur_kowilc(&
             pnode,gptur,gpden,gpvis,gpgra,&
             gpcar,gppr2,gpgrv,eltur,gprea,&
             gpdif,gprhs,conve)
 
! if next line is commented, stabilizes only velocity, if not stabilizes modification of velocity. This last option is right but more nonlinear      
!        gpvel(1:ndime)= conve(1:ndime) ! to stabilize convective terms

     end if

  else if( TUR_K_OMEGA_BREDBERG ) then
     !
     ! Bredberg-Peng-Davidson k-w
     !
     call tur_brepen(&
          ndime,pnode,gptur,gpden,gpvis,gpmut,gpgra,&
          eledd,eltur,gpcar,gppro,gppr2,gprea,gpdif,&
          gprhs,gpgrd)

  else if( TUR_SST_K_OMEGA ) then
     !
     ! Menter SST k-w
     !
     call tur_sstkom(&
          ndime,pnode,gptur,gpden,gpvis,gpmut,gpgra,eledd,&
          eltur,gpcar,gppro,gpwal,gppr2,gprea,gpdif,gprhs,&
          gpgrd,gpgrv,hleng,fddes,gddes,sstf1,sasso)

  else if( TUR_K_EPS_PHI_F ) then
     !
     ! Laurence-Uribe-Utyuzhnikov k-e-phi-f model
     !
     call tur_kepspf(& 
          ndime,pnode,gptur,gpden,gpvis,gpmut,gpwal,&
          gpgra,eledd,eltur,gpcar,gppro,gplap,gprea,&
          gpdif,gprhs,gpgrd)

  else if(TUR_TKE_SGS) then
     !
     ! TKE-SGS model
     !
     call tur_tkesgs(&
          gptur,gpden,gpvis,gpmut,gpgra,&
          gppr2, grnor_tur,hleng, gpwal,gprea,gpdif,gprhs,&
          gpcan ,sreac)
              
  end if
  !
  ! Two layer models: compute inner layer coefficients
  !
  if( TUR_TWO_LAYER_RODI ) then
     !
     ! Two-layer Rodi
     !
     call tur_tworod(& 
          ndime,pnode,nturb_tur,iunkn_tur,param_tur,&
          gptur,gpden,gpvis,gpmut,gpwal,gpgra,eledd,&
          gpust,gpsha,gpcar,gppro,gprea,gpdif,gprhs,&
          gpgrd,gpvel)

  else if( TUR_TWO_LAYER_XU_CHEN ) then
     !
     ! Two-layer Xu-Chen for natural convection
     !
     call tur_twonat(& 
          ndime,pnode,nturb_tur,iunkn_tur,param_tur,&
          gptur,gpden,gpvis,gpmut,gpwal,gpgra,eledd,&
          gpcar,gppro,gprea,gpdif,gprhs,gpgrd,gpvel)
  end if

end subroutine tur_elmco2
