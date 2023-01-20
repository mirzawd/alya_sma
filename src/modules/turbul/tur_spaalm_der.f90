!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_spaalm_der(&
     itask,pnode,gpvis,gpden,gpsha,gpcar,elvel,elmsh,eltur,eledd,&
     gpwal,elust,eltem,elrg2,elsqk,elgrp,elpro,elprd,gpvol,&
     chale,gpdif,gpvel,gpgrv,gptur,gpcar_der,gpwal_der,&
     drea_dtur,dsou_dtur,dsou_dvel,gprea_diff,gprhs_diff,gpreawal_diff,gprhswal_diff)

     
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_spaalm_der
  ! NAME
  !   tur_spaalm_der
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
  use def_kermod, only     :  kfl_dvar_type
  use def_turbul, only     :  nturb_tur,&
       &                      param_tur,&
       &                      ipara_tur
  implicit none
  integer(ip), intent(in)  :: pnode,itask
  real(rp),    intent(in)  :: gpvis,gpden
  real(rp),    intent(in)  :: gpsha(pnode),gpcar(ndime,pnode)
  real(rp),    intent(in)  :: elvel(ndime,pnode)
  real(rp),    intent(in)  :: elmsh(ndime,pnode)
  real(rp),    intent(in)  :: eltur(nturb_tur,pnode)
  real(rp),    intent(in)  :: eledd(pnode),gpwal
  real(rp),    intent(in)  :: elust(pnode),eltem(pnode)
  real(rp),    intent(in)  :: gpvol,chale(2)
  real(rp),    intent(in)  :: elrg2(pnode),elsqk(pnode)
  real(rp),    intent(in)  :: elgrp(ndime,pnode)
  real(rp),    intent(in)  :: elpro(pnode)
  real(rp),    intent(in)  :: elprd(pnode)
  real(rp),    intent(in)  :: gptur
  real(rp),    intent(in)  :: gpgrv(ndime,ndime)
  real(rp),    intent(in)  :: gpdif
  real(rp),    intent(in)  :: gpvel(ndime)
  real(rp),    intent(in)  :: gpcar_der(ndime,pnode,ndime,pnode)
  real(rp),    intent(in)  :: gpwal_der(ndime,pnode)
  
  real(rp),    intent(out) :: drea_dtur(nturb_tur,pnode)
  real(rp),    intent(out) :: dsou_dtur(nturb_tur,pnode)
  real(rp),    intent(out) :: dsou_dvel(ndime,pnode)
  real(rp),    intent(out) :: gprea_diff(ndime,pnode)
  real(rp),    intent(out) :: gprhs_diff(ndime,pnode)
  real(rp),    intent(out) :: gpreawal_diff(ndime,pnode)
  real(rp),    intent(out) :: gprhswal_diff(ndime,pnode)
  
  integer(ip)              :: iprod,idime,jdime,inode,kdime,knode
  real(rp)                 :: fv1,fv2,fw,X,g,r,S,Stild
  real(rp)                 :: cb1,cb2,cv1,sigma,cw1,cw2,cw3,vonka
  real(rp)                 :: siinv,ft2
  real(rp)                 :: k2,d2,Xto3,k2d2,gto6
  real(rp)                 :: ono6,cw3t6,fact1
  real(rp)                 :: gpgrt(3),gprot
  
  real(rp)                 :: sign1,gto5
  real(rp)                 :: dS_dvel(pnode,ndime),dStild_dvel(pnode,ndime)
  real(rp)                 :: dr_dvel(pnode,ndime),dg_dvel(pnode,ndime),dfw_dvel(pnode,ndime)
  
  real(rp)                 :: dX_dtur(pnode),dft2_dtur(pnode),dfv1_dtur(pnode),dfv2_dtur(pnode)
  real(rp)                 :: dr_dtur(pnode),dg_dtur(pnode),dfw_dtur(pnode),dStild_dtur(pnode),dgpgr2_dtur,gpnve
  
  real(rp)                 :: gpgrv_diff(ndime,ndime,ndime,pnode),S_diff(ndime,pnode),omega_diff(ndime,pnode)
  real(rp)                 :: Stild_diff(ndime,pnode),fact1_diff(ndime,pnode),r_diff(ndime,pnode),fw_diff(ndime,pnode)
  real(rp)                 :: g_diff(ndime,pnode),gpgrt_diff(ndime,ndime,pnode),gprot_diff,gpgr2_diff
  real(rp)                 :: fact1wal_diff(ndime,pnode),rwal_diff(ndime,pnode),gwal_diff(ndime,pnode)
  real(rp)                 :: Stildwal_diff(ndime,pnode),fwwal_diff(ndime,pnode),gpwal_ownder(ndime,pnode)
  
  !
  ! Variables
  !
  cb1   = param_tur(1)
  cb2   = param_tur(2)
  sigma = param_tur(3)
  cv1   = param_tur(4)
  cw1   = param_tur(8)
  cw2   = param_tur(5)
  cw3   = param_tur(6)
  vonka = param_tur(7)
  iprod = ipara_tur(1)
  !
  ! Initialization
  !
  drea_dtur = 0.0_rp
  dsou_dtur = 0.0_rp
  dsou_dvel = 0.0_rp
  gprea_diff = 0.0_rp
  gprhs_diff = 0.0_rp  
  !
  ! grad(nu')
  !
  gpgrt(1:ndime) = 0.0_rp
  do inode = 1,pnode
     gpgrt(1:ndime) = gpgrt(1:ndime) + gpcar(1:ndime,inode)*eltur(1,inode)
  end do
  !
  ! sqrt(2*O_ij*O_ij)
  !
  S = 0.0_rp
  do idime = 1,ndime
     do jdime = 1,ndime
        gprot = 0.50_rp*(gpgrv(idime,jdime)-gpgrv(jdime,idime))
        S     = S + gprot*gprot
     end do
  end do
  S = sqrt(2.0_rp*S)
  !
  ! Transition function FT2
  !
  X     = gptur/(gpvis/gpden)                              ! nu'/nu
  Xto3  = X*X*X                                            ! X^3
  siinv = 1.0_rp/sigma                                     ! 1/sigma
  k2    = vonka*vonka                                      ! k^2
  d2    = gpwal*gpwal                                      ! d^2
  k2d2  = k2*d2                                            ! k^2*d^2
  ono6  = 1.0_rp/6.0_rp                                    ! 1/6
  

  fv1   = Xto3 / ( Xto3 + cv1 * cv1 * cv1 )           
  fv2   = 1.0_rp - X / ( 1.0_rp + X * fv1 )
  Stild = S + gptur*fv2/k2d2
  fact1 = Stild*k2d2
  cw3t6 = cw3*cw3*cw3*cw3*cw3*cw3                          ! cw3^6
  r     = min( gptur/fact1 , 10.0_rp )                   ! nu'/(S' k^2 d^2)
  g     = r + cw2*(r*r*r*r*r*r-r)                       ! r + cw2(r^6-r)          
  gto6  = g*g*g*g*g*g
  fw    = g*( (1.0_rp+cw3t6)/(gto6  +cw3t6))**ono6
  ft2   = 1.2_rp*exp(-0.5_rp*X*X)    

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! EXACT LINEARIZATION W.R.T. VELOCITY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!! Here we calculate dsou_dvel => the derivatives of turbulence source terms w.r.t nodal velocities !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Calculate derivatives of S w.r.t nodal velocities
  !
  call vecnor(gpvel,ndime,gpnve,2_ip)
!   if (gpnve /= 0.0) then
  if (gpnve /= 0.0_rp .and. S /= 0.0_rp ) then
  
    if (ndime == 2) then
      ! S = |dux/dy - duy/dx| 
      sign1 = abs(gpgrv(2,1)-gpgrv(1,2))/(gpgrv(2,1)-gpgrv(1,2))
      do inode=1,pnode
        dS_dvel(inode,1) = sign1*(  gpcar(2,inode) )
        dS_dvel(inode,2) = sign1*( -gpcar(1,inode) )
      enddo
    else
      ! S = sqrt( (dux/dy - duy/dx)^2 + (duy/dz - duz/dy)^2 + (duz/dx - dux/dz)^2 )
      do inode=1,pnode
        dS_dvel(inode,1) = (  gpcar(2,inode)*(gpgrv(2,1)-gpgrv(1,2)) -gpcar(3,inode)*(gpgrv(1,3)-gpgrv(3,1)) )/S
        dS_dvel(inode,2) = ( -gpcar(1,inode)*(gpgrv(2,1)-gpgrv(1,2)) +gpcar(3,inode)*(gpgrv(3,2)-gpgrv(2,3)) )/S
        dS_dvel(inode,3) = ( -gpcar(2,inode)*(gpgrv(3,2)-gpgrv(2,3)) +gpcar(1,inode)*(gpgrv(1,3)-gpgrv(3,1)) )/S
      enddo
    endif
  else
    dS_dvel(:,:) = 0.0_rp
  endif
  
  
  !
  ! Calculate derivatives of Stild w.r.t nodal velocities
  !
  do idime = 1,ndime
    do inode=1,pnode
      dStild_dvel(inode,idime) = dS_dvel(inode,idime)
    enddo
  enddo
  !
  ! Calculate derivatives of r=min( gptur/fact1 , 10.0_rp )  w.r.t.  nodal velocities
  !
  if (gptur/fact1 < 10.0_rp) then
    do idime = 1,ndime
      do inode=1,pnode
        dr_dvel(inode,idime) = -dStild_dvel(inode,idime)*gptur/(Stild*Stild*k2d2)
      enddo
    enddo
  else
    dr_dvel(:,:) = 0.0_rp
  endif
  !
  ! Calculate derivatives of g = r + cw2*(r*r*r*r*r*r-r)  w.r.t.  nodal velocities
  !
  do idime = 1,ndime
    do inode=1,pnode
      dg_dvel(inode,idime) = dr_dvel(inode,idime) + cw2*(6*dr_dvel(inode,idime)*r**5 - dr_dvel(inode,idime))
    enddo
  enddo
  !
  ! Calculate derivatives of fw = g*( (1.0_rp+cw3t6)/(gto6  +cw3t6))**ono6  w.r.t.  nodal velocities
  !
  gto5  = g*g*g*g*g
  do idime = 1,ndime
    do inode=1,pnode
      dfw_dvel(inode,idime) = dg_dvel(inode,idime)*((1.0_rp+cw3t6)/(gto6  +cw3t6))**ono6 + &
                              g*ono6* (-6.0_rp*dg_dvel(inode,idime)*gto5*(1.0_rp+cw3t6)/(gto6  +cw3t6)**2.0_rp) * &
                              ((1.0_rp+cw3t6)/(gto6  +cw3t6))**(ono6-1.0_rp)
    enddo
  enddo
  !
  ! Calculate derivatives of source term   -cb1*(1-ft2)*S_tild*nut + cw1*fw*nut²/d²    w.r.t.  nodal velocities
  !
  do inode=1,pnode
	do idime=1,ndime
	  dsou_dvel(idime, inode) = -gpden*cb1*(1.0_rp-ft2)*dStild_dvel(inode,idime)*gptur + gpden*cw1*dfw_dvel(inode,idime)*gptur*gptur/d2
! 	  dsou_dvel(idime, inode) = -gpden*cb1*(1.0_rp-ft2)*dStild_dvel(inode,idime)*gptur + 0.0_rp*gpden*cw1*dfw_dvel(inode,idime)*gptur*gptur/d2
	enddo
  end do 
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! EXACT LINEARIZATION W.R.T. NUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!! Here we calculate the derivatives of turbulence terms w.r.t nodal nut !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !
  ! Calculate derivatives of X  w.r.t. nodal nut
  !
  do inode=1,pnode
    dX_dtur(inode)  = gpsha(inode)/(gpvis/gpden)
  enddo
  !
  ! Calculate derivatives of ft2 = 1.2_rp*exp(-0.5_rp*X*X) w.r.t. nodal nut
  !
  do inode=1,pnode
    dft2_dtur(inode)  = -1.2_rp*0.5_rp*2.0_rp*X*dX_dtur(inode)*exp(-0.5_rp*X*X)
  enddo
  !
  ! Calculate derivatives of fv1   = Xto3 / ( Xto3 + cv1 * cv1 * cv1 ) w.r.t. nodal nut
  !
  do inode=1,pnode
    dfv1_dtur(inode)  = ( 3.0_rp*X*X*dX_dtur(inode)*(Xto3+cv1*cv1*cv1) - 3.0_rp*X*X*dX_dtur(inode)*Xto3 )/ ( Xto3 + cv1 * cv1 * cv1 )**2.0_rp
  enddo
  !
  ! Calculate derivatives of fv2   = 1.0_rp - X / ( 1.0_rp + X * fv1 ) w.r.t. nodal nut
  !
  do inode=1,pnode
    dfv2_dtur(inode)  = - ( dX_dtur(inode)*(1.0_rp+X*fv1) - (dX_dtur(inode)*fv1 + X*dfv1_dtur(inode))*X )/( 1.0_rp + X * fv1 )**2.0_rp
  enddo
  !
  ! Calculate derivatives of Stild = S + gptur*fv2/k2d2 w.r.t. nodal nut
  !
  do inode=1,pnode
    dStild_dtur(inode)  = gpsha(inode)*fv2/k2d2 + gptur*dfv2_dtur(inode)/k2d2
  enddo
  
  !
  ! Calculate derivatives of r=min( gptur/fact1 , 10.0_rp )  w.r.t.  nodal nut
  !
  if (gptur/fact1 < 10.0_rp) then
    do inode=1,pnode
      dr_dtur(inode) = ( gpsha(inode)*Stild - dStild_dtur(inode)*gptur )/(Stild*Stild*k2d2)
!       dr_dtur(inode) = ( gpsha(inode)*Stild - 0.0_rp*dStild_dtur(inode)*gptur )/(Stild*Stild*k2d2)
    enddo
  else
    dr_dtur(:) = 0.0_rp
  endif
  !
  ! Calculate derivatives of g = r + cw2*(r*r*r*r*r*r-r)  w.r.t.  nodal nut
  !
  do inode=1,pnode
    dg_dtur(inode) = dr_dtur(inode) + cw2*(6*dr_dtur(inode)*r**5 - dr_dtur(inode))
  enddo
  !
  ! Calculate derivatives of fw = g*( (1.0_rp+cw3t6)/(gto6  +cw3t6))**ono6  w.r.t.  nodal nut
  !
  gto5  = g*g*g*g*g
  do inode=1,pnode
    dfw_dtur(inode) = dg_dtur(inode)*((1.0_rp+cw3t6)/(gto6+cw3t6))**ono6 + &
                            g*ono6* (-6.0_rp*dg_dtur(inode)*gto5*(1.0_rp+cw3t6)/(gto6+cw3t6)**2.0_rp) * &
                            ((1.0_rp+cw3t6)/(gto6+cw3t6))**(ono6-1.0_rp)
  enddo
  !
 ! Calculate derivatives of -rho*cb1*(1-ft2)*Stild*nut w.r.t. nodal nut (some terms are already calculated in gprea)
  !
  do inode=1,pnode
    drea_dtur(1,inode)  = -gpden*cb1*(1.0_rp-ft2)*dStild_dtur(inode)*gptur + gpden*cb1*Stild*dft2_dtur(inode)*gptur - &
                           gpden*cb1*Stild*gpsha(inode)    
  enddo
  !
  ! Calculate derivatives of rho*(cw1*fw - cb1/k²*ft2)*nut²/d² w.r.t. nodal nut (some terms are already calculated in gprea)
  !
  do inode=1,pnode
    drea_dtur(1,inode)  = drea_dtur(1,inode) &
                      - gpden*cb1/k2d2*dft2_dtur(inode)*gptur*gptur - 2.0_rp*gpden*cb1/k2d2*ft2*gptur*gpsha(inode) &
                      + gpden*cw1*dfw_dtur(inode)*gptur*gptur/d2 + gpden*cw1*fw*gptur*gpsha(inode)/d2
  enddo
  !
  ! Calculate derivatives of -1/sigma*rho*cb2*dnut/dxi*dnut/dxi w.r.t. nodal nut 
  !
  do inode=1,pnode
    dgpgr2_dtur = 0.0_rp
    do idime = 1,ndime
      dgpgr2_dtur = dgpgr2_dtur + 2.0_rp*gpgrt(idime)*gpcar(idime,inode)
    end do
    drea_dtur(1,inode)  = drea_dtur(1,inode) - gpden*siinv*cb2*dgpgr2_dtur
  enddo
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! EXACT LINEARIZATION W.R.T. COORDINATES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  if (kfl_dvar_type == 5 .and. itask == 30) then
  
    !
    ! derivatives of Velocity gradient w.r.t the coordinates d(GPGRV(j,i))/dX
    !
    gpgrv_diff = 0.0_rp
    do kdime = 1,ndime
      do knode=1,pnode

        do inode = 1,pnode
          do idime = 1,ndime
            do jdime = 1,ndime
              gpgrv_diff(jdime,idime,kdime,knode) = gpgrv_diff(jdime,idime,kdime,knode) &
                                                + gpcar_der(jdime,inode,kdime,knode) * elvel(idime,inode)
            end do
          end do
        end do
      
      enddo
    enddo
  
    S_diff = 0.0_rp
    omega_diff = 0.0_rp
    do kdime = 1,ndime
      do knode=1,pnode
      
        do idime = 1,ndime
          do jdime = 1,ndime
            gprot = 0.50_rp*(gpgrv(idime,jdime)-gpgrv(jdime,idime))
            gprot_diff = 0.50_rp*(gpgrv_diff(idime,jdime,kdime,knode)-gpgrv_diff(jdime,idime,kdime,knode))
            omega_diff(kdime,knode) = omega_diff(kdime,knode) + 2.0_rp*gprot*gprot_diff
          end do
        end do
        if (S /= 0.0_rp) S_diff(kdime,knode) = omega_diff(kdime,knode)/S
      
      enddo
    enddo

    do kdime = 1,ndime
      do knode=1,pnode
        Stild_diff(kdime,knode) = S_diff(kdime,knode)
        fact1_diff(kdime,knode) = Stild_diff(kdime,knode)*k2d2
      enddo
    enddo
  
    if (gptur/fact1 < 10.0_rp) then
      do kdime = 1,ndime
        do knode=1,pnode
          r_diff(kdime,knode) = -gptur*fact1_diff(kdime,knode)/(fact1*fact1)
        enddo
      enddo
    else
      r_diff(:,:) = 0.0_rp
    endif
  
    do kdime = 1,ndime
      do knode=1,pnode
        g_diff(kdime,knode) = r_diff(kdime,knode) + cw2*( 6.0_rp*r_diff(kdime,knode)*r*r*r*r*r - r_diff(kdime,knode) )
        fw_diff(kdime,knode) = g_diff(kdime,knode)*((1.0_rp+cw3t6)/(gto6  +cw3t6))**ono6 + &
                              g*ono6* (-6.0_rp*g_diff(kdime,knode)*gto5*(1.0_rp+cw3t6)/(gto6 + cw3t6)**2.0_rp) * ((1.0_rp+cw3t6)/(gto6  +cw3t6))**(ono6-1.0_rp)
      enddo
    enddo
  
    gpgrt_diff = 0.0_rp
    do kdime = 1,ndime
      do knode=1,pnode
        do inode = 1,pnode
          gpgrt_diff(1:ndime,kdime,knode) = gpgrt_diff(1:ndime,kdime,knode) + gpcar_der(1:ndime,inode,kdime,knode)*eltur(1,inode)
        end do
      enddo
    enddo

    do kdime = 1,ndime
      do knode=1,pnode
        gpgr2_diff = 0.0_rp
        do idime = 1,ndime
          gpgr2_diff = gpgr2_diff + 2.0_rp*gpgrt_diff(idime,kdime,knode)*gpgrt(idime)
        end do
        gprhs_diff(kdime,knode) = gpden*( siinv*cb2*gpgr2_diff + cb1*gptur*Stild_diff(kdime,knode) )
        gprea_diff(kdime,knode) = gpden*(cw1*fw_diff(kdime,knode)*gptur/d2+cb1*ft2*Stild_diff(kdime,knode))
      enddo
    enddo

    !
    ! Derivatives of the wall distance w.r.t. coordinates of its own nodes
    !
    do kdime = 1,ndime
      do knode=1,pnode
        gpwal_ownder(kdime,knode) = -gpwal_der(kdime,knode)
      enddo
    enddo
    
    do kdime = 1,ndime
      do knode=1,pnode
        Stildwal_diff(kdime,knode) = -2.0_rp*gptur*fv2*gpwal_ownder(kdime,knode)/(k2d2*gpwal)
        fact1wal_diff(kdime,knode) = Stildwal_diff(kdime,knode)*k2d2 + 2.0_rp*gpwal_ownder(kdime,knode)*Stild*k2*gpwal
      enddo
    enddo
  
    if (gptur/fact1 < 10.0_rp) then
      do kdime = 1,ndime
        do knode=1,pnode
          rwal_diff(kdime,knode) = -gptur*fact1wal_diff(kdime,knode)/(fact1*fact1)
        enddo
      enddo
    else
      rwal_diff(:,:) = 0.0_rp
    endif
  
    do kdime = 1,ndime
      do knode=1,pnode
        gwal_diff(kdime,knode) = rwal_diff(kdime,knode) + cw2*( 6.0_rp*rwal_diff(kdime,knode)*r*r*r*r*r - rwal_diff(kdime,knode) )
        fwwal_diff(kdime,knode) = gwal_diff(kdime,knode)*((1.0_rp+cw3t6)/(gto6  +cw3t6))**ono6 + &
                              g*ono6* (-6.0_rp*gwal_diff(kdime,knode)*gto5*(1.0_rp+cw3t6)/(gto6 + cw3t6)**2.0_rp) * ((1.0_rp+cw3t6)/(gto6  +cw3t6))**(ono6-1.0_rp)
      enddo
    enddo
  
    do kdime = 1,ndime
      do knode=1,pnode
        gprhs_diff(kdime,knode) = gprhs_diff(kdime,knode) + gpden*( cb1*gptur*Stildwal_diff(kdime,knode) - 2.0_rp*cb1*ft2*gptur*gptur*gpwal_ownder(kdime,knode)/(k2d2*gpwal) )
        gprea_diff(kdime,knode) = gprea_diff(kdime,knode) + gpden*( cw1*fwwal_diff(kdime,knode)*gptur/d2 - 2.0_rp*cw1*fw*gptur*gpwal_ownder(kdime,knode)/(d2*gpwal) & 
                                         + cb1*ft2*Stildwal_diff(kdime,knode) )
      enddo
    enddo    
    !
    ! Derivatives of the wall distance w.r.t. coordinates of the nearest wall nodes
    !
    do kdime = 1,ndime
      do knode=1,pnode
        Stildwal_diff(kdime,knode) = -2.0_rp*gptur*fv2*gpwal_der(kdime,knode)/(k2d2*gpwal)
        fact1wal_diff(kdime,knode) = Stildwal_diff(kdime,knode)*k2d2 + 2.0_rp*gpwal_der(kdime,knode)*Stild*k2*gpwal
      enddo
    enddo
  
    if (gptur/fact1 < 10.0_rp) then
      do kdime = 1,ndime
        do knode=1,pnode
          rwal_diff(kdime,knode) = -gptur*fact1wal_diff(kdime,knode)/(fact1*fact1)
        enddo
      enddo
    else
      rwal_diff(:,:) = 0.0_rp
    endif
  
    do kdime = 1,ndime
      do knode=1,pnode
        gwal_diff(kdime,knode) = rwal_diff(kdime,knode) + cw2*( 6.0_rp*rwal_diff(kdime,knode)*r*r*r*r*r - rwal_diff(kdime,knode) )
        fwwal_diff(kdime,knode) = gwal_diff(kdime,knode)*((1.0_rp+cw3t6)/(gto6  +cw3t6))**ono6 + &
                              g*ono6* (-6.0_rp*gwal_diff(kdime,knode)*gto5*(1.0_rp+cw3t6)/(gto6 + cw3t6)**2.0_rp) * ((1.0_rp+cw3t6)/(gto6  +cw3t6))**(ono6-1.0_rp)
      enddo
    enddo
  
    do kdime = 1,ndime
        do knode=1,pnode  
        gprhswal_diff(kdime,knode) = gpden*( cb1*gptur*Stildwal_diff(kdime,knode) - 2.0_rp*cb1*ft2*gptur*gptur*gpwal_der(kdime,knode)/(k2d2*gpwal) )
        gpreawal_diff(kdime,knode) = gpden*( cw1*fwwal_diff(kdime,knode)*gptur/d2 - 2.0_rp*cw1*fw*gptur*gpwal_der(kdime,knode)/(d2*gpwal) & 
                                         + cb1*ft2*Stildwal_diff(kdime,knode) )
      enddo
    enddo
  
  endif !itask
  

  if (kfl_dvar_type == 6 .and. itask == 30) then
  
!     gprhs_diff(1,1) = gpden*cb1*gptur*Stild
!     gprhs_diff(1,1) = gpden*cb1*gptur*fw*gptur/k2d2
!     gprea_diff(1,1) = gpden*(cb1*ft2*Stild)
!     gprea_diff(1,1) = 1.0_rp !gpden*(cb1*1.00000_rp/k2d2) 
    
  endif
    
    
    
    
    
    
    
    
    
    
end subroutine tur_spaalm_der
