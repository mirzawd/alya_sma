!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_spaalm(&
     ndime,pnode,nturb_tur,ipara_tur,param_tur,eltur,gptur,gpvis,&
     gpden,gpgrv,gpwal,gpcar,gprea,gpdif,gprhs,gpvel)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_spaalm
  ! NAME
  !   tur_spaalm
  ! DESCRIPTION
  !    Compute coefficient of the equation of Spalart-Almmaras model
  ! OUTPUT 
  !    GPREA .......... r 
  !    GPDIF .......... k 
  !    GPRHS .......... f 
  !    GPGRD(NDIME) ... grad(k) coefficient
  ! USES
  ! USED BY
  !    tur_elmcoe
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_kermod, only     :  kfl_adj_prob
  implicit none
  integer(ip), intent(in)  :: ndime,pnode,nturb_tur
  integer(ip), intent(in)  :: ipara_tur(*)
  real(rp),    intent(in)  :: param_tur(*)
  real(rp),    intent(in)  :: eltur(nturb_tur,pnode,3)
  real(rp),    intent(in)  :: gptur
  real(rp),    intent(in)  :: gpvis,gpden,gpwal
  real(rp),    intent(in)  :: gpcar(ndime,pnode)
  real(rp),    intent(in)  :: gpgrv(ndime,ndime)
  real(rp),    intent(out) :: gpdif,gprea
  real(rp),    intent(out) :: gprhs
  real(rp),    intent(inout) :: gpvel(ndime)
  integer(ip)              :: iprod,idime,jdime,inode,imeth
  real(rp)                 :: fv1,fv2,fw,X,g,r,S,Stild,fv3
  real(rp)                 :: cb1,cb2,cv1,sigma,cw1,cw2,cw3,vonka
  real(rp)                 :: siinv,ft2,ct3,ct4
  real(rp)                 :: k2,d2,Xto3,k2d2,gto6
!  real(rp)                 :: Xocv2
  real(rp)                 :: ono6,cw3t6,fact1
!  real(rp)                 :: cv2
  real(rp)                 :: gpgrt(3),gprot,gpgr2
!  real(rp)                 :: gpgrt(3),gprot,gpgr2,fact2
!  real(rp)                 :: Sbar
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
  ! grad(nu')
  !
  gpgrt(1:ndime) = 0.0_rp
  do inode = 1,pnode
     gpgrt(1:ndime) = gpgrt(1:ndime) + gpcar(1:ndime,inode)*eltur(1,inode,1)
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
  ct3   = 1.1_rp
  ct4   = 2.0_rp
  ft2   = ct3*exp(max(-1.0e2_rp,-ct4*X*X))
  !
  ! Some usefull constants     
  !
  Xto3  = X*X*X                                            ! X^3
  siinv = 1.0_rp/sigma                                     ! 1/sigma
  k2    = vonka*vonka                                      ! k^2
  d2    = gpwal*gpwal                                      ! d^2
  k2d2  = k2*d2                                            ! k^2*d^2
  ono6  = 1.0_rp/6.0_rp                                    ! 1/6
  
  
  
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !
  ! Viscous functions FV2 and FV3 according to the model 
  !
!   if( iprod == 0 ) then                                    ! Original SA + low Reynolds correction
!      fv1   = Xto3 / ( Xto3 + cv1 * cv1 * cv1 )           
!      fv2   = 1.0_rp - X / ( 1.0_rp + X * fv1 )
!      fv3   = 1.0_rp
!   else if( iprod == 1 ) then                               ! Geuzaine, Delanaye and Liu corrections
!      fv1   = Xto3 / ( Xto3 + cv1 * cv1 * cv1 )    
!      cv2   = 5.0_rp
!      Xocv2 = 1.0_rp + X/cv2
!      fv2   = 1.0_rp / ( Xocv2 * Xocv2 * Xocv2 )
!      if( X < 1.0e-3 ) then
!         fv3 = 3.0_rp / ( cv2 + 3.0_rp * X )
!      else
!         fv3 = ( 1.0_rp + X * fv1 ) * ( 1.0_rp - fv2 ) / X
!      end if
!   else                                                     ! Original SA
!      fv1 = 1.0_rp
!      fv2 = 0.0_rp
!      fv3 = 1.0_rp
!      ft2 = 0.0_rp
!   end if
!   !
!   ! Destruction (wall) function FW
!   !
!   Stild = S*fv3 + gptur*fv2/k2d2                           ! S*fv3 + nu' fv2/(k^2 d^2)
!   fact1 = Stild*k2d2
!   cw3t6 = cw3*cw3*cw3*cw3*cw3*cw3                          ! cw3^6
!   if( gptur >= fact1*3.0_rp ) then
!      fw    = (1.0_rp+cw3t6)**ono6
!   else
!      r     = min( gptur/fact1 , 3.0_rp )                   ! nu'/(S' k^2 d^2)
!      g     = r + cw2*(r*r*r*r*r*r-r)                       ! r + cw2(r^6-r)          
!      gto6  = g*g*g*g*g*g
!      fw    = g*( (1.0_rp+cw3t6)/(gto6+cw3t6))**ono6
!   end if

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  fv1   = Xto3 / ( Xto3 + cv1 * cv1 * cv1 )           
  fv2   = 1.0_rp - X / ( 1.0_rp + X * fv1 )
  fv3   = 1.0_rp
  Stild = S + gptur*fv2/k2d2
  fact1 = Stild*k2d2
  cw3t6 = cw3*cw3*cw3*cw3*cw3*cw3                          ! cw3^6
  r     = min( gptur/fact1 , 10.0_rp )                   ! nu'/(S' k^2 d^2)
  g     = r + cw2*(r*r*r*r*r*r-r)                       ! r + cw2(r^6-r)          
  gto6  = g*g*g*g*g*g
  fw    = g*( (1.0_rp+cw3t6)/(gto6  +cw3t6))**ono6
  ft2   = 1.2_rp*exp(-0.5_rp*X*X)    
  
!   fv1   = Xto3 / ( Xto3 + cv1 * cv1 * cv1 )           
!   fv2   = 1.0_rp - X / ( 1.0_rp + X * fv1 )
!   fv3   = 1.0_rp
!   c2 = 0.7_rp
!   c3 = 0.9_rp
!   Sbar = gptur*fv2/k2d2
!   if (Sbar >= -c2*S) then
!     Stild = S + Sbar
!   else
!     Stild = S + S*(c2**2*S + c3*Sbar)/((c3-2*c2)*S-Sbar)
!   endif
!   fact1 = Stild*k2d2
!   cw3t6 = cw3*cw3*cw3*cw3*cw3*cw3                          ! cw3^6
!   r     = min( gptur/fact1 , 10.0_rp )                   ! nu'/(S' k^2 d^2)
!   g     = r + cw2*(r*r*r*r*r*r-r)                       ! r + cw2(r^6-r)          
!   gto6  = g*g*g*g*g*g
!   fw    = g*( (1.0_rp+cw3t6)/(gto6  +cw3t6))**ono6
!   ft2   = 1.2_rp*exp(-0.5_rp*X*X)    
  
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Equation coefficients 
  !
  gprea = gpden*(cw1*fw*gptur/d2+cb1*ft2*Stild)            ! rho*(cw1*fw*nu'/d^2+cb1*ft2*S')  
  gpdif = siinv*(gpvis+gpden*gptur)                        ! rho*(nu+nu')/sigma
  
  if (kfl_adj_prob == 0) then
    imeth = 2
  else
    imeth = 1
  endif
  
  if( imeth == 1 ) then
     gpgr2 = 0.0_rp
     do idime = 1,ndime
        gpgr2 = gpgr2 + gpgrt(idime)*gpgrt(idime)
     end do
     gprhs = gpden*(siinv*cb2*gpgr2+cb1*gptur&             ! rho*(cb1*nu'*(S'+ft2*nu'/(k^2*d^2)
          &  *(Stild+ft2*gptur/k2d2))                      ! +cb2/sig*grad(nu')^2)
     
  else
     
     gpgr2 = 0.0_rp
     do idime = 1,ndime
        gpgr2 = gpgr2 + gpgrt(idime)*gpgrt(idime)
     end do
     gprhs = gpden*(siinv*cb2*gpgr2+cb1*gptur&             ! rho*(cb1*nu'*(S'+ft2*nu'/(k^2*d^2)
          &  *(Stild+ft2*gptur/k2d2))                      ! +cb2/sig*grad(nu')^2)
          
!      do idime = 1,ndime
!         gpvel(idime) = gpvel(idime)-siinv*cb2*gpgrt(idime) ! a-cb2/sig*grad(nu')
!      end do  
!      gprhs = gpden*(cb1*gptur*(Stild+ft2*gptur/k2d2))

     gprea = gprea + gpden*cw1*fw*gptur/d2
     gprhs = gprhs + gpden*cw1*fw*gptur/d2*gptur
     
  end if
  
  !
  ! Compressibility terms
  !
  !fact1 = 0.0_rp
  !fact2 = 0.0_rp
  !do idime = 1,ndime
  !   fact1 = fact1 + gpgrt(idime) * grden(idime)
  !   fact2 = fact2 + grden(idime) * grden(idime)
  !end do
  !do idime = 1,ndime
  !   gprhs = gprhs + siinv*cb2* ( gptur * fact1 + gptur*gptur/(4.0_rp*gpden)*fact2 )
  !end do
  
!   gprea = gprea + gpden*cw1*fw*gptur/d2
!   gprhs = gprhs + gpden*cw1*fw*gptur/d2*gptur
  
!   gprea = gprea - 2.0_rp*gpden*cb1*ft2*gptur/k2d2
!   gprhs = gprhs - gpden*cb1*ft2*gptur*gptur/k2d2
! 
!   gprea = gprea - gpden*cb1*Stild
!   gprhs = gprhs - gpden*cb1*gptur*Stild
  
end subroutine tur_spaalm
