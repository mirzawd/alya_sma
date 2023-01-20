!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_updbcs(itask)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_updbcs
  ! NAME 
  !    tur_updbcs
  ! DESCRIPTION
  !    This routine loads the turbulence boundary conditions
  !    1. Before a time step begins
  !    2. Before a global iteration begins
  !    3. Before an inner iteration begins
  ! USED BY
  !    tur_begste
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_parame, only       :  pi
  use def_kermod
  use mod_ker_proper
  use def_domain
  use def_turbul
  use mod_ker_space_time_function
  use mod_ker_regularization, only : regul_k, regul_e, inv_regul_k, inv_regul_e,  kfl_regularization
  use mod_frivel,             only : frivel

  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,ibopo,jbopo,iunkn,i1,i2,dummi,befit
  integer(ip)             :: ifunc,itype
  integer(ip)             :: kk
  real(rp)                :: xx
  real(rp)                :: nu,ustar,tveno,mu(1),rho(1),Sr,krp,bs,ustar2
  real(rp)                :: wemix,ystar,ystarv,leps,vv,eps_inn,l
  real(rp)                :: k,w,y,eps,wLRN,wHRN,f,Cmu,kap,cv1t3
  real(rp)                :: dupdy,nuT,lmix,Re,lambda,u,xfact
  real(rp)                :: newnu,yplus,wfirs,Res,rebc1_tur,xsmal,umbra
  real(rp)                :: dummy,bvaux, reguk
  

  if( IMASTER ) return
  umbra=1.0e-6_rp
  kap       = 0.41_rp
  Cmu       = 0.09_rp
  bs        = 0.09_rp
  xsmal     = 1.0e-10_rp
  rebc1_tur = 1.0_rp-rebcs_tur
  if(TUR_K_EPS_STD.or.TUR_FAMILY_K_OMEGA) Cmu = param_tur(6)

  if(kfl_algor_tur==1.and.itask/=TUR_BEFORE_GLOBAL_ITERATION.and.itask/=TUR_BEFORE_TIME_STEP) then
     i1=iunkn_tur
     i2=iunkn_tur        
  else
     i1=1
     i2=nturb_tur
  end if

  turbulent_variables: do iunkn=i1,i2
     
     !
     ! Transient fields
     !
     if( itask==TUR_BEFORE_TIME_STEP .and. kexist_tran_fiel > 0 ) then   ! needs checkin if done only before time step is enough
        do ipoin=1,npoin

           ifunc = kfl_funno_tur(ipoin,iunkn)
           itype = kfl_funtn_tur(ipoin,iunkn)

           if( itype == FUNCTION_SPACE_TIME ) then
              !
              ! Space time
              !
              kk = k_tran_fiel(ifunc)
              xx = x_tran_fiel(ifunc)
              bvaux = xfiel(ifunc) % a(1,ipoin,kk) * xx + xfiel(ifunc) % a(1,ipoin,kk+1) * (1.0_rp-xx)
              if (kfl_logva==1) bvaux = log(bvaux)
           end if

           if( itype /= 0 ) then
              if( kfl_timei_tur /= 0 .and. kfl_tiacc_tur == 2 .and. kfl_tisch_tur == 1 ) then
                 bvess_tur(1,ipoin,iunkn) = 0.50_rp * ( bvaux + untur(iunkn,ipoin,nprev_tur) )
              else
                 bvess_tur(1,ipoin,iunkn) = bvaux
              end if
           end if
           
        end do  
     end if     
     if( itask==TUR_BEFORE_TIME_STEP .and. kfl_exist_fixi7_tur==1 ) then   ! needs checkin if done only before time step is enough
!!$     To prescribe inflow outflow following velocity bcs
!!$          do ipoin = 1,npoin
!!$
!!$           ifunc = kfl_funno_nsi(ipoin)
!!$           itype = kfl_funtn_nsi(ipoin)
!!$
!!$           call ker_functions(ipoin,ifunc,itype,bvess_nsi(:,ipoin,2),vefun)
!!$           !
!!$           ! Velocity update according to time discretization
!!$           !
!!$           if( ifunc /= 0 ) then
!!$              if( kfl_timei_nsi /= 0 .and. kfl_tiacc_nsi == 2 .and. kfl_tisch_nsi == 1 ) then
!!$                 do idime = 1,ndime
!!$                    veold = veloc(idime,ipoin,ncomp_nsi)
!!$                    bvess_nsi(idime,ipoin,1) = 0.50_rp*(vefun(idime)+veold)
!!$                 end do
!!$              else
!!$                 do idime = 1,ndime
!!$                    bvess_nsi(idime,ipoin,1) = vefun(idime)
!!$                 end do
!!$              end if
!!$           end if
!!$        end do
!!$      
        !
        ! Open or Close preliminary steps - the rest is done later
        !
        call memgen(1_ip,npoin,0_ip)   !allocate memory for integer gisca
        call open_close(2_ip,kfl_fixno_tur(1,1,iunkn),dummy,1_ip)   ! The last parameter (1_ip) corresponds to the first dimension of kfl_fixno_tur
        call memgen(3_ip,npoin,0_ip)   !deallocate memory for integer gisca
     end if

     if(iunkn==2.and.itask==TUR_BEFORE_INNER_ITERATION) then
        call tur_projec()
     end if
     if (kfl_ustar==2) then
        befit = TUR_BEFORE_INNER_ITERATION
     else
        befit = TUR_BEFORE_GLOBAL_ITERATION
     end if
     nodes: do ipoin=1,npoin
 
        ibopo = lpoty(ipoin)
        call ker_proper('DENSI','IPOIN',ipoin,dummi,rho(1:))
        call ker_proper('VISCO','IPOIN',ipoin,dummi,mu(1:))
      
        nu = mu(1)/rho(1) 
        if( kfl_rough > 0 ) then  ! variable roughness
           rough_dom = rough(ipoin) 
!        else  ! unifform or not roughness use rough_dom value
        end if
        

        type_bc: if(kfl_fixno_tur(1,ipoin,iunkn) == 3.and.itask==BEFIT) then

           !----------------------------------------------------------------------
           !
           ! Wall law
           !
           !----------------------------------------------------------------------

           if (kfl_delta == 1 ) then
              y = ywalp(ibopo)
           else
              y = delta_tur
           end if
           call vecnor(advec(1,ipoin,1),ndime,tveno,2_ip)
           call frivel(kfl_ustar,y,rough_dom,tveno,nu,ustar)

           yplus = y*ustar/nu

           if(TUR_SPALART_ALLMARAS) then
              !
              ! Spalart-Allmaras
              !
              cv1t3 = param_tur(4)**3
              dupdy = &                                             ! |W+|=|dU+/dy+|
                   0.4_rp/(kap*(1.0_rp+0.4_rp*yplus))&              ! U+ being known
                   +7.8_rp/11.0_rp*( exp(-yplus/11.0_rp)&           ! from the law of the wall
                   +exp(-0.33_rp*yplus)*(0.33_rp*yplus-1.0_rp))               
              lmix  = kap*yplus*(1.0_rp-exp(-yplus/26.0_rp))        ! l = ky+(1-e(-y+/26))
              nuT   = nu*lmix*lmix*abs(dupdy)                       ! nuT=nu l^2 du/dy    
              call tur_spanut(1_ip,nu,newnu,nuT)                    ! Obtain nu' knowing nu_t  
              bvess_tur(1,ipoin,1) = newnu                          ! nu_t+ = k^2 y+^2 [1-exp(-y+/26)]^2 dU+/dy+  
              !bvess_tur(1,ipoin,1) = max(1.0e-7_rp,newnu)           ! nu_t+ = clipping
           
           else if(TUR_FAMILY_K_EPS.or.TUR_K_XU_CHIEN) then
              !
              ! k-eps
              !
              if (kfl_ustar==2.and.iunkn==2) then ! two scale velocites for wall law
                 ! ustar ^2 in terms of k

                 if (kfl_regularization) then
                    reguk = regul_k(untur(1, ipoin,1))
                 else
                    reguk = untur(1, ipoin,1)
                 end if

                 ustar2 = reguk*sqrt(Cmu)         ! u*=k^0.5*cmu^0.25      
                 y   = walld(ipoin)

                 !                    if (ustar2.gt.1.0e-2_rp) then   ! this is to introduce higher eps.
                 ustar = sqrt(ustar2)
                 !                    else
                 !                       ustar2 = ustar*ustar
                 !                    end if


                 !                  please, note that it is not good that epsilon depends on veloc, 
                 !                  when velocity increases
                 eps = ustar2*ustar/(kap*(y+rough_dom))     ! eps=U*^3/(kappa*y)=Cmu^{3/4}*k^{3/2}/(kap*y)

                 if (kfl_regularization)  eps= inv_regul_e(eps)

                 bvess_tur(1,ipoin,2) = eps ! max(eps,xsmal)
              else if (kfl_ustar.eq.1) then   ! one velocity scale for wal lar
                 y   = walld(ipoin)
                 k   = ustar*ustar/sqrt(Cmu)                           ! k=U*^2/sqrt(Cmu)
                 eps = ustar*ustar*ustar/(kap*(y+rough_dom))           ! eps=U*^3/(kappa*y)=Cmu^{3/4}*k^{3/2}/(kap*y)                 
                 if(iunkn==1) bvess_tur(1,ipoin,1) = k
                 if(iunkn==2) bvess_tur(1,ipoin,2) = max(eps,xsmal)
              end if

           else if(TUR_FAMILY_K_OMEGA) then
              !
              ! k-w
              !
               if (kfl_ustar==2.and.iunkn==2) then ! two scale velocities
                
                 ! ustar ^2 in terms of k
                    
                 reguk = untur(1, ipoin,1)
                 
                 ustar2 = reguk*sqrt(Cmu)         ! u*=k^0.5*cmu^0.25      
                 y   = walld(ipoin)
                 ustar = sqrt(ustar2)
                 w = ustar/(sqrt(cmu)*kap*(y+rough_dom))     ! eps=U*^3/(kappa*y)=Cmu^{3/4}*k^{3/2}/(kap*y)
                 bvess_tur(1,ipoin,2) = w ! max(eps,xsmal)
              else if (kfl_ustar.eq.1) then ! one velocity scale
                 
                 !
                 ! k-w
                 !
                 y   = walld(ipoin)
                 k = ustar*ustar/sqrt(cmu)                              ! k=U*^2/sqrt(B*)
                 w = ustar/(sqrt(cmu)*kap*(y+rough_dom))                ! w=U*/(sqrt(beta*)*kap*y)
                 
                 if(iunkn==1) bvess_tur(1,ipoin,1) = k
                 if(iunkn==2) bvess_tur(1,ipoin,2) = max(w,xsmal)
              end if
              
           end if

        else if( kfl_fixno_tur(1,ipoin,iunkn) == 1 .and. itask == TUR_BEFORE_TIME_STEP ) then

           !----------------------------------------------------------------------
           !
           ! Prescribed value
           !
           !----------------------------------------------------------------------

           if( kfl_funtn_tur(ipoin,iunkn) == FUNCTION_SPACE_TIME  ) then 
              ifunc = kfl_funno_tur(ipoin,iunkn)            
              call ker_space_time_function(&
                   ifunc,coord(1,ipoin),coord(2,ipoin),coord(ndime,ipoin),cutim,xfact)
              bvess_tur(1,ipoin,iunkn) = xfact
           end if 
           if( kfl_adj_prob == 1 ) then 
              bvess_tur(1,ipoin,iunkn) = 0.0_rp
           end if 

        else if(kfl_fixno_tur(1,ipoin,iunkn) == 4.and.itask==TUR_BEFORE_INNER_ITERATION) then

           !----------------------------------------------------------------------
           !
           ! No-slip wall
           !
           !----------------------------------------------------------------------

           if(TUR_SPALART_ALLMARAS) then
              !
              ! Spalart-Allmaras
              !
              bvess_tur(1,ipoin,1) = 0.0_rp

           else if(TUR_K_EPS_CHIEN.or.TUR_K_XU_CHIEN.or.TUR_K_EPS_LAUNDER_SHARMA) then
              !
              ! Chien k-eps
              !
              if(iunkn==1) bvess_tur(1,ipoin,1) = 0.0_rp              ! k   = 0
              if(iunkn==2) bvess_tur(1,ipoin,2) = 0.0_rp              ! eps = 0

           else if(TUR_FAMILY_K_EPS) then
              !
              ! k-eps
              !
              if(iunkn==1) then
                 bvess_tur(1,ipoin,1) = 0.0_rp                        ! k = 0
              else if(iunkn==2) then
                 eps = 2.0_rp*nu*grk12_tur(ibopo)*grk12_tur(ibopo)    ! eps = 2*nu*(d sqrt(K)/dy )^2
                 bvess_tur(1,ipoin,2)=rebc1_tur*bvess_tur(1,ipoin,2)& !     = 2*nu*[(sqrt(K)_{first_node}-0)/y]^2
                      +rebcs_tur*eps
              end if

           else if(TUR_FAMILY_K_OMEGA) then
              !
              ! k-w
              !
              if(iunkn==1) then
                 bvess_tur(1,ipoin,1)=0.0_rp
              else
                 if(rough_dom<=zetur.and.rough_dom>=-zetur) then    ! Smooth wall

                    if (kfl_model_tur == 32)  then ! Special treatment for sst komega - it does not use the general value for k-omega
                                                   ! For the moment we have no adaptive treament similar to Bredberg
                                                   ! The cases with roughness are treated by the general k omega

                       w = 800.0_rp*nu*grono_tur(ibopo)   !  wfirst = (60/0.075)*nu/(y^2)
                       w = max(1.0e-10_rp,w)                         ! Added this minimum limit as is done in the general k-omega  

                    else   ! rest of cases

                       wLRN=2.0_rp*nu/bs*grono_tur(ibopo)              ! wfirst = 2*nu/((b*)*y^2)

                       if(kfl_wallw_tur==0) then
                          wHRN  = 0.0_rp
                          f     = 1.0_rp
                       else
                          !
                          ! Wfirst = f*wLRN + (1-f)*wHRN
                          ! J. Bredberg and L. Davidson
                          ! Low-Reynolds Number Turbulence Models: An Approach
                          ! for Reducing Mesh Sensitivity
                          ! Transactions of the ASME, 126, 14-21 (2004)
                          !
                          wHRN  = grk12_tur(ibopo)/(Cmu**0.25_rp*kap)  ! (k^{1/2}/y)/(Cmu^{1/4}*kappa)
                          ystar = greps_tur(ibopo)/nu**0.75_rp         ! y*=u_eps*y/nu
                          f     = exp(-ystar/17.0_rp)                  ! f=exp(-y*/17)

                       end if

                       wfirs = f*wLRN+(1.0_rp-f)*wHRN
                       w     = max(1.0e-10_rp,10.0_rp*wfirs)
                    end if

                 else if(rough_dom<-zetur) then
                    ustar = ustar_tur(lpoty(ipoin))
                    krp   = abs(rough_dom)
                    if(krp<25) then
                       Sr=(50.0_rp/krp)**2
                    else
                       Sr=100.0_rp/krp
                    end if
                    w = ustar*ustar/nu*Sr

                 else
                    ustar = ustar_tur(lpoty(ipoin))             
                    krp   = ustar*rough_dom/nu
                    if(krp<25) then
                       w = (50.0_rp/rough_dom)**2*nu
                    else
                       w = 100.0_rp*ustar/rough_dom
                    end if

                 end if

                 bvess_tur(1,ipoin,2)=rebc1_tur*bvess_tur(1,ipoin,2)+rebcs_tur*w

              end if

           end if

        else if(( kfl_fixno_tur(1,ipoin,iunkn) == 6 .or. kfl_fixno_tur(1,ipoin,iunkn) == 8 ) &
             .and. itask == TUR_BEFORE_GLOBAL_ITERATION) then

           !----------------------------------------------------------------------
           !
           ! Inflow condition
           !
           ! Internal flows: I = 1% to 10%
           ! External flows: I down to 0.05%
           !
           ! For fully dvped pipe flow
           ! -------------------------
           ! Re     3,000  5,000  10,000  15,000  20,000  >>
           ! nut/nu  11.6   16.5    26.7    34.0    50.1  100.0
           !
           ! For freestream
           ! --------------
           ! nut/nu=0.1-0.2
           !
           !----------------------------------------------------------------------

           if( kfl_infl1_tur == 2 .or. kfl_infl1_tur == 3 .or. kfl_infl2_tur == 2 ) then  
              !
              ! U* is needed
              !
              Re  = avvel_tur*hdiam_tur/nu
              if(Re==0.0_rp) then
                 lambda = 0.0_rp
              else if(Re<=2000.0_rp) then
                 lambda = 64.0_rp/Re
              else if(Re<=4000.0_rp) then
                 lambda = 0.021377_rp+5.3155e-6_rp*Re
              else
                 lambda = 1.0_rp/(1.8_rp*log10(Re)-1.64_rp)**2
              end if
              ustar = avvel_tur*sqrt(lambda/8.0_rp)

           else if( kfl_infl1_tur == 4 .or. kfl_infl2_tur == 4 ) then
              y   = walld(ipoin)                                       ! y              
              k   = usref*usref/sqrt(Cmu)                              ! k = u*^2/sqrt(Cmu) 
              eps = usref*usref*usref/(kap*(y+k_ref))                  ! eps=U*^3/(kappa*y)
           end if
           !
           ! 1st variable
           !
           if(kfl_infl1_tur==1) then                                   ! 1. Intensity
              k     = 1.5_rp*(turin_tur*avvel_tur)**2                  ! k = 3/2*(I*u)**2

           else if(kfl_infl1_tur==2) then                              ! 2. Channel flow
              y     = walld(ipoin)                                     ! y
              yplus = y*ustar/nu                                       ! y+  = y*u*/nu
              Res   = ustar*hdiam_tur/nu                               ! Re* = rho*(u*)*D/mu
              k     = ustar**2*(0.07_rp*yplus**2*exp(-yplus/8.0_rp)&   ! k   = u*^2*[0.07y+*exp(-y+/8)
                   &  +4.5_rp*(1.0_rp-exp(-yplus/20.0_rp))&            !       +4.5*(1-exp(-y+/20))
                   &  /(1.0_rp+4.0_rp*yplus/Res))                      !       /(1+4y+/Re*)

           else if(kfl_infl1_tur==3) then                              ! 3. Smooth regime
              k     = ustar*ustar/sqrt(Cmu)                            ! k = u*^2/sqrt(Cmu) 

           end if
           !
           ! 2nd variable
           !
           if(kfl_infl2_tur==1) then                                   ! 1. Length scale
              if(turle_tur==0.0_rp) then
                 l = 0.1_rp*Cmu**0.25_rp*kap*hdiam_tur                 ! l   = Cmu^{1/4}*kap*Dh/10
              else if(turle_tur>0.0_rp) then
                 l = kap*turle_tur                                     ! l = kap*L
              else
                 l = kap*max(walld(ipoin),1.0e-6_rp)                   ! l = kap*y
              end if
              eps  = (Cmu*k*k)**0.75_rp/l                              ! eps = Cmu^{3/4} k^{3/2} / l
              w    = sqrt(k)/l                                         ! w   = sqrt(k)/l
              nut  = sqrt(k)*l                                         ! nut = sqrt(k)*l

           else if(kfl_infl2_tur==2) then                              ! 2. Channel                 
              y     = walld(ipoin)                                     ! y
              yplus = y*ustar/nu                                       ! y+  = y*u*/nu
              eps  = ustar**4/(nu*kap&                                 ! 
                   &  *(yplus**4_rp+15.0_rp**4_rp)**0.25_rp)           ! eps = u*^4/(nu*kap*(y+^4+15^4)^{1/4})
              w    = eps/(Cmu*(max(k,1.0e-6_rp)))                      ! w   = eps/(Cmu*k)
              nut  = Cmu*k*k/eps                                       ! nut = Cmu*k^2/eps

           else if(kfl_infl2_tur==3) then                              ! 3. Ratio viscosities
              xfact = nu*nutnu_tur
              call tur_nut2nd(1_ip,ipoin,nu,xfact,k,u)
              eps  = u                                                 ! eps = Cmu*k^2/(nu * nut/nu )
              w    = u                                                 ! w   = k/(nu * nut/nu )
              nut  = u                                                 ! nut = nu*(nut/nu)

           end if

           if(TUR_SPALART_ALLMARAS) then
              if( kfl_infl2_tur /= 0 ) bvess_tur(1,ipoin,1) = nut
           else if( TUR_FAMILY_K_EPS .or. TUR_K_XU_CHIEN ) then
              if( iunkn == 1 .and. kfl_infl1_tur /= 0 ) bvess_tur(1,ipoin,1) = k
              if( iunkn == 2 .and. kfl_infl2_tur /= 0 ) bvess_tur(1,ipoin,2) = eps
           else if(TUR_FAMILY_K_OMEGA) then
              if( iunkn == 1 .and. kfl_infl1_tur /= 0 ) bvess_tur(1,ipoin,1) = k
              if( iunkn == 2 .and. kfl_infl2_tur /= 0 ) bvess_tur(1,ipoin,2) = w
           end if

        else if(kfl_fixno_tur(1,ipoin,iunkn) == 5 .and. itask==TUR_BEFORE_INNER_ITERATION) then

           !----------------------------------------------------------------------
           !
           ! Inner-layer
           !
           !----------------------------------------------------------------------

           if(TUR_TWO_LAYER_RODI) then
              !
              ! Rodi two-layer k-eps
              !
              if(iunkn==2) then
                 jbopo  = lwnei_tur(ipoin)
                 yplus  = walld(ipoin)*ustar_tur(jbopo)/nu
                 ystar  = walld(ipoin)*sqrt(untur(1,ipoin,1))/nu            ! y*   = y*k^{1/2}/nu
                 call mixing(1_ip,wemix,param_tur(8),param_tur(9),ystar)
                 vv     = untur(1,ipoin,1)*ystar*(0.0004_rp+0.0000465_rp*ystar) ! vv   = k*(y*)*(4e-04+4.65e-5*y*)
                 ystarv = walld(ipoin)*sqrt(vv)/nu                          ! y*v  = y*vv^{1/2}/nu
                 if(ystarv<1.0e-10_rp) then
                    eps_inn = 0.0_rp
                 else
                    leps    = 1.3_rp*walld(ipoin)/(1.0_rp+2.12_rp/ystarv)   ! leps = 1.3*y/(1+2.12/y*v) 
                    eps_inn = sqrt(vv)*untur(1,ipoin,1)/leps
                 end if
                 !untur(2,ipoin,1)=wemix*bvess_tur(2,ipoin)+(1.0_rp-wemix)*eps_inn
                 !untur(2,ipoin,1)=eps_inn
                 !if(wemix>0.5_rp) then
                 !   kfl_fixno_tur(junkn_tur,ipoin)=5
                 !else
                 !   kfl_fixno_tur(junkn_tur,ipoin)=6
                 !end if
              end if

           end if

        else if(abs(kfl_fixno_tur(1,ipoin,iunkn)) == 7 .and. itask==TUR_BEFORE_TIME_STEP ) then   ! hhh beware I am nit totaly sure about itask           
           !
           ! Open or closes point depending on angle between normal and VELOC
           !
           if( gisca(ipoin) > 0 ) then
              kfl_fixno_tur(1,ipoin,iunkn) = 7
           else
              kfl_fixno_tur(1,ipoin,iunkn) = -7
           end if

        end if type_bc

     end do nodes ! ipoin

     if( itask==TUR_BEFORE_TIME_STEP .and. kfl_exist_fixi7_tur==1 ) call memgen(3_ip,npoin,0_ip)   !deallocate gisca


  end do turbulent_variables
  !
  ! When solving a manufactured solution, impose exact Dirichlet bc
  !
  !if( kfl_exacs_tur /= 0 ) then
  !   call tur_manufactured_Dirichlet_condition()
  !end if
  
end subroutine tur_updbcs
