!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_der_co2_all(&
     pnode,gpvis,gpden,gpsha,gpcar,elvel,elmsh,eltur,eledd,&
     elwal,elust,eltem,elrg2,elsqk,elgrp,elpro,elprd,gpvol,&
     chale,gpdif,gpvel,gpgrv,gptur,&
     hleng,gpcan, elprr, pmate, &
     drea_dtur,dsou_dtur, &
     dsou_dvel,gpmut,dgpmut,dgpmut_dvel)

     
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_der_co2_all
  ! NAME
  !   tur_der_co2_all
  ! DESCRIPTION
  !    Compute derivative of coefficient of the equation of SST K-W
  ! OUTPUT 
  ! USES
  ! USED BY
  !    tur_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  ID_TURBUL,ID_ALEFOR
  use def_domain, only     :  ndime
  use def_turbul, only     :  nturb_tur,iunkn_tur,&
       &                      kfl_produ_tur,&
       &                      TUR_SPALART_ALLMARAS,&
       &                      TUR_SST_K_OMEGA,  &
       &                      kfl_clipp_tur
  implicit none
  integer(ip), intent(in)  :: pnode
  real(rp),    intent(in)  :: gpvis,gpden
  real(rp),    intent(in)  :: gpsha(pnode),gpcar(ndime,pnode)
  real(rp),    intent(in)  :: elvel(ndime,pnode)
  real(rp),    intent(in)  :: elmsh(ndime,pnode)
  real(rp),    intent(in)  :: eltur(nturb_tur,pnode)
  real(rp),    intent(in)  :: eledd(pnode),elwal(pnode)
  real(rp),    intent(in)  :: elust(pnode),eltem(pnode)
  real(rp),    intent(in)  :: gpvol,chale(2)
  real(rp),    intent(in)  :: elrg2(pnode),elsqk(pnode)
  real(rp),    intent(in)  :: elgrp(ndime,pnode)
  real(rp),    intent(in)  :: elpro(pnode)
  real(rp),    intent(in)  :: elprd(pnode)
  real(rp),    intent(in)  :: elprr(pnode)
  real(rp),    intent(in)  :: gptur(nturb_tur)
  real(rp),    intent(in)  :: gpgrv(ndime,ndime)
  real(rp),    intent(in)  :: gpdif
  real(rp),    intent(in)  :: gpvel(ndime)
  real(rp),    intent(in)  :: gpmut
  real(rp),    intent(in)  :: dgpmut(nturb_tur,pnode)
  real(rp),    intent(in)  :: dgpmut_dvel(ndime,pnode)
  
  real(rp),    intent(out) :: drea_dtur(nturb_tur,pnode)
  real(rp),    intent(out) :: dsou_dtur(nturb_tur,pnode)
  real(rp),    intent(out) :: dsou_dvel(ndime,pnode)
  
  real(rp),    intent(in)  :: hleng(ndime)
  real(rp),    intent(in)  :: gpcan
  integer(ip), intent(in)  :: pmate
  
  integer(ip)              :: idime,inode,jdime
  real(rp)                 :: gppro,gpwal
  real(rp)                 :: gppr2
  real(rp)                 :: gpkin, gpome, dgpkin(pnode), dgpome(pnode),dgppr2_dvel(pnode,ndime),dgppro_dvel(pnode,ndime)
  real(rp)                 :: F1,CD
  real(rp)                 :: gradk(ndime),gradw(ndime),grakw,c1,c2,c3,F11,divve
  real(rp)                 :: dgradk(ndime,pnode), dgradw(ndime,pnode),dgrakw_dk(pnode),dgrakw_dw(pnode)
  real(rp)                 :: sigk1,sigw1,beta1,sigk2,sigw2,beta,gama
  real(rp)                 :: beta2,betas,gama1,gama2,arg1,xsmal
!  real(rp)                 :: P
  real(rp)                 :: Crc
  
  real(rp)                 :: dCD_dk(pnode), dCD_dw(pnode),darg1_dk(pnode),darg1_dw(pnode)
  real(rp)                 :: dc1_dk(pnode),dc1_dw(pnode),dc2_dk(pnode),dc2_dw(pnode),dc3_dk(pnode),dc3_dw(pnode)
  real(rp)                 :: dF1_dk(pnode), dF1_dw(pnode), dF11_dk(pnode), dF11_dw(pnode)
  real(rp)                 :: dgama_dk(pnode), dgama_dw(pnode), dbeta_dk(pnode), dbeta_dw(pnode)
  real(rp)                 :: dgpmut_dkin(pnode),dgpmut_dome(pnode)

  if (.not. TUR_SST_K_OMEGA) return
  
  !
  ! Initialization
  !
  drea_dtur = 0.0_rp
  dsou_dtur = 0.0_rp
  dsou_dvel = 0.0_rp
  divve = 0.0_rp
  dgpkin = 0.0_rp
  dgpome = 0.0_rp
  dgppr2_dvel = 0.0_rp
  !
  ! Definitions
  !
  sigk1 = 0.85_rp          
  sigw1 = 0.5_rp               
  beta1 = 0.075_rp
  sigk2 = 1.0_rp
  sigw2 = 0.856_rp            
  beta2 = 0.0828_rp
  betas = 0.09_rp
  xsmal = 1.0e-10_rp
  gama1 = 5.0_rp/9.0_rp
  gama2 = 0.44_rp
  Crc = 1.4_rp  ! Magic number for rotating reference frame

  gradk = 0.0_rp
  gradw = 0.0_rp
  grakw = 0.0_rp
  dgrakw_dk = 0.0_rp
  dgrakw_dw = 0.0_rp
  dCD_dw = 0.0_rp
  dCD_dk = 0.0_rp
  darg1_dk = 0.0_rp
  darg1_dw = 0.0_rp
  dF1_dk = 0.0_rp
  dF1_dw = 0.0_rp
  dF11_dk = 0.0_rp
  dF11_dw = 0.0_rp
  dgama_dk = 0.0_rp
  dgama_dw = 0.0_rp
  dbeta_dk = 0.0_rp
  dbeta_dw = 0.0_rp

  !
  ! Derivatives of gauss point values (k,w) w.r.t nodal values
  !
  gpkin = max(0.0_rp, gptur(1))
  gpome = gptur(2)
  do inode=1,pnode
    if (gptur(1) > 0.0_rp) dgpkin(inode) = gpsha(inode)
    dgpome(inode) = gpsha(inode)
  enddo

  do inode = 1,pnode
    dgpmut_dkin(inode) = dgpmut(1,inode)
    dgpmut_dome(inode) = dgpmut(2,inode)
  enddo
  !
  ! Gradients of gauss point values (k,w) and its derivatives w.r.t nodal values
  !
  do inode = 1,pnode
    do idime = 1,ndime
	gradk(idime) = gradk(idime) + eltur(1,inode) * gpcar(idime,inode)
	gradw(idime) = gradw(idime) + eltur(2,inode) * gpcar(idime,inode)
	dgradk(idime,inode) = gpcar(idime,inode)
	dgradw(idime,inode) = gpcar(idime,inode)
    end do
  end do    

  do idime = 1,ndime
    grakw = grakw + gradk(idime) * gradw(idime)
    do inode = 1,pnode
      dgrakw_dk(inode) = dgrakw_dk(inode) + dgradk(idime,inode) * gradw(idime)
      dgrakw_dw(inode) = dgrakw_dw(inode) + gradk(idime) * dgradw(idime,inode)
    enddo
  end do
  !
  ! calculate gpwall
  !
  gpwal = 0.0_rp
  do inode = 1,pnode
    gpwal = gpwal + gpsha(inode) * elwal(inode)
  end do

  !
  ! Calculate CD and its derivatives w.r.t nodal values
  !
  if( ( gpome > xsmal ) .or. ( kfl_clipp_tur > 1 ) ) then
    CD = max( 1.0e-10_rp , 2.0_rp*gpden*sigw2*grakw/gpome )
    if (2.0_rp*gpden*sigw2*grakw/gpome>1.0e-10_rp ) then
      do inode = 1,pnode
	dCD_dk(inode) = 2.0_rp*gpden*sigw2*dgrakw_dk(inode)/gpome
	dCD_dw(inode) = 2.0_rp*gpden*sigw2*(dgrakw_dw(inode)/gpome-grakw*dgpome(inode)/gpome**2.0_rp)
      enddo
    else
      do inode = 1,pnode
	dCD_dk(inode) = 0.0_rp
	dCD_dw(inode) = 0.0_rp
      enddo
    endif
  else
    CD = max( 1.0e-10_rp , 2.0_rp*gpden*sigw2*grakw/xsmal )
    if (2.0_rp*gpden*sigw2*grakw/xsmal>1.0e-10_rp ) then
      do inode = 1,pnode
	dCD_dk(inode) = 2.0_rp*gpden*sigw2*dgrakw_dk(inode)/xsmal
	dCD_dw(inode) = 2.0_rp*gpden*sigw2*dgrakw_dw(inode)/xsmal
      enddo
    else
      do inode = 1,pnode
	dCD_dk(inode) = 0.0_rp
	dCD_dw(inode) = 0.0_rp
      enddo
    endif
  end if
  
  !
  ! Calculate arg1 and its derivatives w.r.t nodal values
  !
  if( ( gpome > xsmal ) .or. ( kfl_clipp_tur > 1 ) ) then
    c1   = sqrt( gpkin ) / ( betas * gpome * gpwal )
    c2   = 500.0_rp * gpvis / ( gpden * gpwal * gpwal * gpome )
    c3   = 4.0_rp * gpden * sigw2 * gpkin / ( CD * gpwal * gpwal )
    do inode = 1,pnode
      dc1_dk(inode) = -0.5_rp*dgpkin(inode) / ( sqrt( gpkin )*betas * gpome * gpwal )
      dc1_dw(inode) = - dgpome(inode)*sqrt( gpkin ) / ( betas * gpome**2.0_rp * gpwal )
      dc2_dk(inode) = 0.0_rp
      dc2_dw(inode) = -500.0_rp * dgpome(inode) * gpvis / ( gpden * gpwal * gpwal * gpome**2.0 )
      dc3_dk(inode) = 4.0_rp * gpden * sigw2 * (dgpkin(inode) / ( CD * gpwal * gpwal ) - gpkin * dCD_dk(inode)/ ( CD*CD * gpwal * gpwal ))
      dc3_dw(inode) = -4.0_rp * gpden * sigw2 * gpkin * dCD_dw(inode) / ( CD*CD * gpwal * gpwal )
    enddo
    arg1 = min ( max (  c1 , c2 ) , c3 )
    do inode = 1,pnode
      if (c1 >= c2) then
	if ( c1 >= c3) then
	  darg1_dk(inode) = dc3_dk(inode)
	  darg1_dw(inode) = dc3_dw(inode)
	else
	  darg1_dk(inode) = dc1_dk(inode)
	  darg1_dw(inode) = dc1_dw(inode)
	endif    
      else
	if ( c2 >= c3) then
	  darg1_dk(inode) = dc3_dk(inode)
	  darg1_dw(inode) = dc3_dw(inode)
	else
	  darg1_dk(inode) = dc2_dk(inode)
	  darg1_dw(inode) = dc2_dw(inode)
	endif
      endif
    enddo !inode
  else
    do inode = 1,pnode
      darg1_dk(inode) = dc3_dk(inode)
      darg1_dw(inode) = dc3_dw(inode)
    enddo
  end if

  !
  ! Calculate F1 and F11 and its derivatives w.r.t nodal values
  !
  F1 = tanh( arg1**4.0_rp )
  do inode = 1,pnode
    dF1_dk(inode) = 4.0_rp*arg1**3.0_rp*darg1_dk(inode)*(1.0_rp-F1*F1)
    dF1_dw(inode) = 4.0_rp*arg1**3.0_rp*darg1_dw(inode)*(1.0_rp-F1*F1)
  enddo  
  F11          = 1.0_rp - F1
  do inode = 1,pnode
    dF11_dk(inode) = -dF1_dk(inode)
    dF11_dw(inode) = -dF1_dw(inode)
  enddo
  
  !
  ! Calculate gama and beta and its derivatives w.r.t nodal values
  !
  gama         = F1 * gama1 + F11 * gama2 ! gama
  beta         = F1 * beta1 + F11 * beta2 ! beta
  do inode = 1,pnode
    dgama_dk(inode)         = dF1_dk(inode) * gama1 + dF11_dk(inode) * gama2
    dgama_dw(inode)         = dF1_dw(inode) * gama1 + dF11_dw(inode) * gama2
    dbeta_dk(inode)         = dF1_dk(inode) * beta1 + dF11_dk(inode) * beta2
    dbeta_dw(inode)         = dF1_dw(inode) * beta1 + dF11_dw(inode) * beta2
  enddo
 
  !
  ! Calculate gppr2 = (dui/dxj+duj/dxi)*dui/dxj for shear production GPPRO = mu_t*(dui/dxj+duj/dxi)*dui/dxj
  !
  if( .not. TUR_SPALART_ALLMARAS ) then
      if( kfl_produ_tur == 0 ) then 
        gppr2 = 0.0_rp
        do idime = 1,ndime
          do jdime = 1,ndime
             gppr2 = gppr2 + (gpgrv(jdime,idime) + gpgrv(idime,jdime))* gpgrv(jdime,idime) ! 2 S_ij : S_ij
          end do
        end do
        gppr2 = gppr2 - divve*divve/3.0_rp 
        gppro = gpmut * gppr2 ! just to check the limits
      else
          print *, "DERIVATION IS NOT IMPLEMENTED"
      end if
  end if
    
  !
  ! Calculate derivatives of gppr2 and gppro w.r.t nodal velocities
  !
  if( .not. TUR_SPALART_ALLMARAS ) then
      if( kfl_produ_tur == 0 ) then 
        ! calculate dgppr2_dvel
        if (ndime == 2) then
          do inode=1,pnode
	    dgppr2_dvel(inode,1) = 4.0_rp*gpcar(1,inode)*gpgrv(1,1) + 2.0_rp*gpcar(2,inode)*(gpgrv(1,2) + gpgrv(2,1))
	    dgppr2_dvel(inode,2) = 4.0_rp*gpcar(2,inode)*gpgrv(2,2) + 2.0_rp*gpcar(1,inode)*(gpgrv(1,2) + gpgrv(2,1))
          enddo
        else
          do inode=1,pnode
	    dgppr2_dvel(inode,1) = 4.0_rp*gpcar(1,inode)*gpgrv(1,1) + 2.0_rp*gpcar(2,inode)*(gpgrv(1,2) + gpgrv(2,1)) + &
								      2.0_rp*gpcar(3,inode)*(gpgrv(1,3) + gpgrv(3,1))
	    dgppr2_dvel(inode,2) = 4.0_rp*gpcar(2,inode)*gpgrv(2,2) + 2.0_rp*gpcar(1,inode)*(gpgrv(1,2) + gpgrv(2,1)) + &
								      2.0_rp*gpcar(3,inode)*(gpgrv(2,3) + gpgrv(3,2))
	    dgppr2_dvel(inode,3) = 4.0_rp*gpcar(3,inode)*gpgrv(3,3) + 2.0_rp*gpcar(1,inode)*(gpgrv(1,3) + gpgrv(3,1)) + &
								      2.0_rp*gpcar(2,inode)*(gpgrv(2,3) + gpgrv(3,2))
          enddo
        endif ! idime
        ! calculate dgppro_dvel
        do idime=1,ndime
          do inode=1,pnode
	    dgppro_dvel(inode,idime) = gpmut*dgppr2_dvel(inode,idime) + dgpmut_dvel(idime,inode)*gppr2
          enddo
        enddo
      else
          print *, "DERIVATION IS NOT IMPLEMENTED"
      end if
  end if

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Here we calculate dsou_dtur => the derivatives of turbulence source terms w.r.t turbulence unknowns (k and w) !!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(iunkn_tur == 1) then ! FOR K
  
    !
    !         P     = min( gppro , 10.0_rp*gpden*betas*gpome*gpkin )
    !         gprhs = P
    if (gppro < 10.0_rp*gpden*betas*gpome*gpkin) then ! check the limits
      do inode=1,pnode
	dsou_dtur(1,inode) = dgpmut_dkin(inode) * gppr2
	dsou_dtur(2,inode) = dgpmut_dome(inode) * gppr2
      end do !inode
    else
      do inode=1,pnode
	dsou_dtur(1,inode) = 10.0_rp*gpden*betas*gpome*dgpkin(inode)
	dsou_dtur(2,inode) = 10.0_rp*gpden*betas*dgpome(inode)*gpkin
      end do !inode
    endif
      
  else ! FOR W
  
    ! 
    !        P     = min( gppr2 , 10.0_rp*gpden*betas*gpome*gpkin/gpmut )
    !        gprhs = gama*gpden*P
    if ( gppr2 < 10.0_rp*gpden*betas*gpome*gpkin/gpmut) then ! check the limits
      do inode=1,pnode
	dsou_dtur(1,inode) = dgama_dk(inode)*gpden*gppr2
	dsou_dtur(2,inode) = dgama_dw(inode)*gpden*gppr2
      end do !inode
    else
      do inode=1,pnode
	dsou_dtur(1,inode) = 10.0_rp*gama*gpden*gpden*betas*gpome* &
			    (dgpkin(inode)/gpmut - gpkin*dgpmut_dkin(inode)/gpmut**2.0_rp) + &
			    10.0_rp*gpden*gpden*gpome*gpkin/gpmut*betas*dgama_dk(inode)
	dsou_dtur(2,inode) = 10.0_rp*gama*gpden*gpden*betas*gpkin* &
			    (dgpome(inode)/gpmut - gpome*dgpmut_dome(inode)/gpmut**2.0_rp) + &
			    10.0_rp*gpden*gpden*gpome*gpkin/gpmut*betas*dgama_dw(inode)
      end do !inode
    endif
    !
    !  gprhs  = gprhs + gprhs  = gprhs + 2.0_rp*F11*gpden*sigw2*grakw/max(gpome,xsmal)
    !
    if (gpome < xsmal) then ! check the limits
      do inode=1,pnode
	dsou_dtur(1,inode) = dsou_dtur(1,inode) + 2.0_rp*gpden*sigw2/xsmal*(dF11_dk(inode)*grakw + F11*dgrakw_dk(inode))
	dsou_dtur(2,inode) = dsou_dtur(2,inode) + 2.0_rp*gpden*sigw2/xsmal*(dF11_dw(inode)*grakw + F11*dgrakw_dw(inode))
      end do !inode
    else
      do inode=1,pnode
	dsou_dtur(1,inode) = dsou_dtur(1,inode) + 2.0_rp*gpden*sigw2/gpome*(dF11_dk(inode)*grakw + F11*dgrakw_dk(inode))
	dsou_dtur(2,inode) = dsou_dtur(2,inode) + 2.0_rp*gpden*sigw2/gpome*(dF11_dw(inode)*grakw + F11*dgrakw_dw(inode))  &
						    - 2.0_rp*gpden*sigw2*F11*grakw*dgpome(inode)/gpome**2.0_rp
      end do !inode	
    endif
  endif ! iunkn_tur

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Here we calculate drea_dtur => the derivatives of turbulence reaction terms w.r.t turbulence unknowns (k and w) !!
  ! by assuming kfl_rotat_tur = 0.0                                                                                         !!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  if(iunkn_tur == 1) then ! FOR K
    !         gprea = gpden * betas * gpome * gpkin
    do inode=1,pnode
      drea_dtur(1,inode) = gpden*betas*gpome*dgpkin(inode)
      drea_dtur(2,inode) = gpden*betas*gpkin*dgpome(inode)
    end do !inode
  else  ! FOR W
    !           gprea = gpden * beta * gpome * gpome
    do inode=1,pnode
      drea_dtur(1,inode) = gpden * dbeta_dk(inode) * gpome * gpome
      drea_dtur(2,inode) = gpden * dbeta_dw(inode) * gpome * gpome + 2.0_rp*gpden*beta*gpome*dgpome(inode)
    end do !inode
  endif ! iunkn_tur
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !!!!!!!! Here we calculate dsou_dvel => the derivatives of turbulence source terms w.r.t nodal velocities !!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(iunkn_tur == 1) then ! FOR K
  
    !
    !         P     = min( gppro , 10.0_rp*gpden*betas*gpome*gpkin )
    !         gprhs = P
    !
    if (gppro < 10.0_rp*gpden*betas*gpome*gpkin) then ! check the limits
      do inode=1,pnode
	do idime=1,ndime
	  dsou_dvel(idime, inode) = dgppro_dvel(inode,idime)
	enddo
      end do !inode
    endif
      
  else ! FOR W
  
    ! 
    !        P     = min( gppr2 , 10.0_rp*gpden*betas*gpome*gpkin/gpmut )
    !        gprhs = gama*gpden*P
    !
    if ( gppr2 < 10.0_rp*gpden*betas*gpome*gpkin/gpmut) then ! check the limits
      do inode=1,pnode
	do idime=1,ndime
	  dsou_dvel(idime, inode) = gama*gpden*dgppr2_dvel(inode,idime)
	enddo
      end do !inode
    else
      do inode=1,pnode
	do idime=1,ndime 
	  dsou_dvel(idime, inode) = dsou_dvel(idime, inode) -10.0_rp*gama*gpden*dgpmut_dvel(idime,inode)* &
	                                  gpden*betas*gpome*gpkin/gpmut**2.0_rp
	enddo
      end do !inode
    endif
    
  endif ! iunkn_tur    
    
    
   
end subroutine tur_der_co2_all
