!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_kowilc(&
     pnode,gptur,gpden,gpvis,gpgra,&
     gpcar,gppr2,gpgrv,eltur,gprea,&
     gpdif,gprhs, conve)
  !------------------------------------------------------------------------
  !****f* Turbul/tur_kowilc
  ! NAME
  !   tur_komega
  ! DESCRIPTION
  !    Compute coefficients of the equation of k-w model developed by Wilcox in 1988, 1998 and 2006
  !    The formulation is detailed in this link https://turbmodels.larc.nasa.gov/wilcox.html 
  !    Transport equations
  !    -------------------
  !   
  !       dK         dK                         d +-          dK  -+
  !    rho-- + rho*ui--- = P - beast*rho*K*w + ---|(mu+mut/sk)---  |
  !       dt         dxi                       dxi+-          dxi -+
  !
  !       dw                       dK   dw           w                     d +-          dw  -+
  !    rho-- + rho*(ui-rho*sigde/w --- )--- =  gamma - P - beta*rho*w^2 + ---|(mu+mut/sw)---  |
  !       dt                       dxi  dxi          K                    dxi+-          dxi -+
  !
  !    for k-w 88 and k-w 98 gamma w P/K = gamma*rho*gppr2   
  !
  !    ipara_tur(1)= 1  k-w 88
  !    ipara_tur(1)= 2  k-w 98
  !    ipara_tur(1)= 3  k-w 2006
  !    ipara_tur(1)= 4  k-w 2006 low Reynolds
  !    Auxiliary relations
  !    -------------------
  !   mixing length and dissipation
  !    L = k^{1/2}/w, eps = b* w k
  !    mu_t=rho*(alpha*)*k/w
  !
  ! OUTPUT 
  !    GPREA .......... r 
  !    GPDIF .......... k 
  !    GPRHS .......... f 
  ! USES
  ! USED BY
  !    tur_elmco2
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
  use def_turbul, only       :  nturb_tur,iunkn_tur,param_tur,ipara_tur, inv_l_max
  implicit none
  integer(ip), intent(in)    :: pnode 
  real(rp),    intent(in)    :: gptur(nturb_tur, 2)
  real(rp),    intent(in)    :: gpgra,gpden,gpvis,gppr2
  real(rp),    intent(in)    :: gpgrv(ndime,ndime) ! velocity gradient
  real(rp),    intent(in)    :: eltur(nturb_tur,pnode)
  real(rp),    intent(out)   :: gpdif
  real(rp),    intent(inout) :: gprea
  real(rp),    intent(in)    :: gpcar(ndime,pnode)
  real(rp),    intent(out)   :: gprhs
  real(rp),    intent(inout) :: conve(ndime)                       !< conective term
  real(rp)                   :: gppro   !< Gauss point Shear production GPPRO = mu_t*(dui/dxj+duj/dxi)*dui/dxj! TKE Production
  real(rp)                   :: gpmut(2),gpkin,gpome,gpkio,gpomo, divve
  real(rp)                   :: bast0, beast,beta0, clim,gamma
  real(rp)                   :: alphast,ReT
  real(rp)                   :: omehat, omehao, beta, sigde, gammap, lm
  real(rp)                   :: gradk(3),gradw(3),grakw,xk,fbs,fb,xw
  real(rp)                   :: S(ndime,ndime),W(ndime,ndime), diffu, linco
  integer(ip)                :: idime,jdime,kdime,inode
  logical                    :: WIL88 =.false., WIL98=.false., WIL06=.false., lowRe=.false.
  
  bast0 = 0.09_rp ! beta coefficient (Cmu)
  if (ipara_tur(1) == 1_ip) then ! k-w Wilcox '88
     WIL88=.true.
     gamma = 5.0_rp/9.0_rp 
     beta  = 3.0_rp/40.0_rp     
  elseif (ipara_tur(1) == 2_ip) then  ! k-w Wilcox '98
     WIL98=.true.
     gamma = 13.0_rp/25.0_rp  
     beta0 = 9.0_rp/125.0_rp
  elseif (ipara_tur(1) == 3_ip.or.ipara_tur(1)==4_ip) then   ! k-w Wilcox 2006
     WIL06=.true.
     gamma = 13.0_rp/25.0_rp  
     beta0 = 0.0708_rp    
  end if
  if (ipara_tur(1) == 4_ip) LowRe =.true.  ! Low Reynolds version of k-w 2006

  !gauss point unknowns evaluated at last and previous iterations
  gpkin = gptur(1,1) 
  gpome = gptur(2,1)
  !unknows values frozen along inner iterations
  gpkio = gptur(1,2)
  gpomo = gptur(2,2)
 

  if (lowre) then  ! Low Reynolds version of k-w 2006
     ReT = gpden*gpkin/(gpvis*gpome)
     alphast = (beta0/3.0_rp + ReT/6.0_rp)/(1.0_rp+ ReT/6.0_rp ) 
     bast0 = bast0*(beta0/0.27_rp + (ReT/8.0_rp)**4.0_rp )/(1.0_rp + (ReT/8.0_rp)**4.0_rp  )
     if (iunkn_tur==2) gamma = gamma* (1.0_rp/9.0_rp + ReT/2.61_rp ) /((1.0_rp + ReT/2.61_rp  )*alphast)
  else
     alphast =1.0_rp
  end if

  if (WIL06 ) then ! only for k-w 06 
     Clim = 7.0_rp/8.0_rp
  else 
     Clim= -999.0_rp ! no effect
  end if
    
  ! turbulent visco
  omehat = max(gpome, Clim*sqrt(alphast*gppr2/bast0))  ! = gpome for k-w 88 and 98
  omehao = max(gpomo, Clim*sqrt(alphast*gppr2/bast0)) 
  gpmut(1) = max(alphast*gpden*gpkin/omehat, 0.000001_rp*gpvis) ! varies along iterations (nonlinear term)
  gpmut(2) = max(alphast*gpden*gpkio/omehao, 0.000001_rp*gpvis) ! frozen along iterations
  ! TKE production (limited for WIL06)
  gppro = gpmut(2) * gppr2

  ! diffusive term = rho *K/w is not gpmut for WIL06
  if (.not.WIL06) then
     diffu = gpmut(iunkn_tur)
  else if (WIL06) then   ! diffusive term
     diffu = max(alphast*gpden*gptur(1, iunkn_tur)/gptur(2,iunkn_tur), 0.000001_rp*gpvis) 
  end if
  
  ! symmetic and antisymmetric strain rate tensor S and W
  divve = 0.0_rp
  do idime = 1,ndime
     divve = divve + gpgrv(idime, idime)
  end do
  do idime = 1,ndime
     do jdime = 1,ndime
        S(idime,jdime) = 0.50_rp * ( gpgrv(jdime,idime) + gpgrv(idime,jdime) )
        W(idime,jdime) = 0.50_rp * ( gpgrv(jdime,idime) - gpgrv(idime,jdime) )
     end do
     S(idime,idime) =  S(idime,idime)  - 0.5_rp * divve
  end do
  
  
  if(iunkn_tur==1) then ! TKE equation
    
     if( WIL88.or.WIL06) then
        beast = bast0  ! beta_ast
     else if ( WIL98 ) then !load beta_ast
        gradk = 0.0_rp
        gradw = 0.0_rp
        do inode = 1,pnode
           do idime = 1,ndime
              gradk(idime) = gradk(idime) + eltur(1,inode) * gpcar(idime,inode)
              gradw(idime) = gradw(idime) + eltur(2,inode) * gpcar(idime,inode)
           end do
        end do
        grakw = 0.0_rp ! = gradk*gradw
        do idime = 1,ndime
           grakw = grakw + gradk(idime) * gradw(idime)
        end do
        xk = grakw/(gpome*gpome*gpome)
        if( xk <= 0.0_rp ) then
           fbs = 1.0_rp
        else
           fbs = ( 1.0_rp + 680.0_rp * xk * xk ) / ( 1.0_rp + 400.0_rp * xk * xk )
        end if      
        beast  = bast0 * fbs  
     end if
     ! coefficients
     gprea = gpden*beast*gpome                      ! (b*)*rho*w
!     sreac = gprea                                   !  stabilized
     gprhs = gppro + gpgra                           !  P+G
!     print *, 'balance', gprhs - gprea*gpkin, gprhs
     
     
  else if (iunkn_tur==2) then  ! OMEGA EQUATION 
     !mixing length
     lm = sqrt(gpkin)/((bast0**0.25_rp)*gpomo)
     if (lm.lt.0.0_rp) lm = 1.0_rp/inv_l_max
     
     ! values for gamma, beta, 
     if ( WIL88) then
!        gamma = 5.0_rp/9.0_rp 
!        beta  = 3.0_rp/40.0_rp
        gammap = gamma + (beta/bast0 - gamma)*lm*inv_l_max
        gprhs = gammap*gpden*gppr2 
     else if ( WIL98.or.WIL06) then
!        gamma = 13.0_rp/25.0_rp  
        xw = 0.0_rp
        do idime = 1,ndime
           do jdime = 1,ndime
              do kdime = 1,ndime
                 xw = xw + W(idime,jdime)*W(jdime,kdime)*S(kdime,idime)
              end do
           end do
        end do
        xw  = abs(xw/( (bast0*gpome)**3.0_rp ))
        ! closure functions and coefficients
        if ( WIL98) then
!           beta0 = 9.0_rp/125.0_rp
           fb  = ( 1.0_rp + 70.0_rp * xw ) / ( 1.0_rp + 80.0_rp * xw )          
           beta = beta0*fb
           ! max. mixing lentgh method
           gammap = gamma + (beta/bast0 - gamma)*lm*inv_l_max
           gprhs  = gammap*gpden*gppr2    
        else if ( WIL06) then
!           beta0 = 0.0708_rp
           fb  = ( 1.0_rp + 85.0_rp * xw ) / ( 1.0_rp + 100.0_rp * xw )
           beta = beta0*fb  
           gradk = 0.0_rp
           gradw = 0.0_rp
           do inode = 1,pnode
              do idime = 1,ndime
                 gradk(idime) = gradk(idime) + eltur(1,inode) * gpcar(idime,inode)
                 gradw(idime) = gradw(idime) + eltur(2,inode) * gpcar(idime,inode)
              end do
           end do
           grakw = 0.0_rp ! = gradk*gradw
           do idime = 1,ndime
              grakw = grakw + gradk(idime) * gradw(idime)
           end do
           if (grakw.gt.0.0_rp) then
              sigde = 1.0_rp/8.0_rp
              ! convective term 
              conve(1:ndime) = conve(1:ndime)  - sigde/gpomo*gradk(1:ndime) 
           end if
           gammap = gamma + (beta/bast0 - gamma)*lm*inv_l_max
           gprhs  = gammap*gppro*gpomo/gpkio
        end if
!        print *, 'Dgamma, lm, lmax',gammap-gamma, lm, 1.0/inv_lmax 
     end if
     ! Newton Raphson integration  w^2 approx 2* wi*w(i+1) - wi
     linco  = gpden*beta*gpome ! linear coeff
     gprea  = 2.0_rp*linco                           ! b0*fb*w
!     sreac  = gprea                                 ! stabilization of reactive term
     gprhs  = gprhs +linco*gpome              ! linearization
  end if
  !
  ! GPDIF diffusion term 
  ! 
!  call tur_elmdif(pnode,gpvis,diffu, gpcar,eledd,gpdif,gpgrd)
  gpdif = gpvis + diffu/param_tur(iunkn_tur)  ! diffusion/sigma_x
end subroutine tur_kowilc
