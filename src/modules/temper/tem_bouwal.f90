!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_bouwal(&
     lboel,iboun,gbsha,pnodb,pnode,&
     acden,acvis,acsph,actco,prtur,&
     twall,rough,kinen,baloc,velex, veave,&
     arobi,trobi,qrobi, gbvel,ustar)
  
!------------------------------------------------------------------------
!
! Wall law for temperature
!
! This routine updates the surface parameters for the heat equation at
! a given integration point of a boundary IBOUN received by argument.
! The parameters will be saved in the array TPARA, as follows:
!
! q_presc=alfa*(T-Ta)+q
!
! Parameters alfa, Ta and q are computed using the law of the wall
! for the temperature
!
! Pr    = Prandtl number
! Prt   = Turbulent Prandtl number
! Rt    = Pr/Prt
! Dt    = exp(-y+/11 Rt^beta)      
!          +-      
!          | 1.1   for Rt<0.1
! beta  = -+
!          | 0.333 for Rt>0.1
!          +-
! Pt    = 9.24(Rt^{0.75}-1)(1+0.28 exp(-0.007Rt))      
! T+    = Pr u+ Dt + Prt(1-Dt)(u+ + Pt)
!
! The wall heat flux qw is given by:
! 
!         rho cp(Tw-T)u*
! qw    = --------------, therefore
!              T+
!
! alfa  = -rho cp u*/T+
! Ta    =  Tw
! q     =  0      
!
! NB: In the limit y+ -> 0, we have T+=Pr y+=Pr y u*/nu so that
!         rho cp(Tw-T)nu
! qw    = --------------
!             Pr y
! 
!------------------------------------------------------------------------
  use def_kintyp
  use def_kermod, only : kfl_delta, kfl_ustar, dexlo_ker, kfl_waexl_ker, &
       cmu_st, kfl_wlaav_ker
  use def_temper, only     : delta_tem
  use def_domain, only     : ndime,ywalb
  use mod_frivel, only     : frivel
  implicit  none
  integer(ip), intent(in)  :: pnodb,pnode,iboun
  integer(ip), intent(in)  :: lboel(pnodb)
  real(rp),    intent(in)  :: actco,prtur,acvis,acden,twall
  real(rp),    intent(in)  :: gbsha(pnodb) 
  real(rp),    intent(in)  :: baloc(ndime,ndime),rough, kinen
  real(rp),    intent(in)  :: velex(3), veave(3), gbvel(ndime)
  real(rp),    intent(out) :: arobi,trobi,qrobi, ustar
  real(rp)                 :: yplus,uplus,Pr,Rt
  real(rp)                 :: alpha,Dt,Pt,Tplus,vonka
  real(rp)                 :: tveno,vikin,hconv,acsph, velfr
  real(rp)                 :: vewal(ndime), delta_aux, tvelo(ndime)
  integer(ip)              :: idime, ldime

  arobi = 0.0_rp
  trobi = 0.0_rp
  
  !
  ! wall distance
  !
  if ( kfl_waexl_ker == 0_ip ) then  !classical behaviour
     if( kfl_delta == 1 ) then
        delta_aux = ywalb(iboun)                                  ! variable wall distance
     else
        delta_aux = delta_tem                                     ! fixed wall distance
     end if
  else                              ! Exchange location
     if( kfl_delta == 1 ) then
        delta_aux = ywalb(iboun)                                  ! variable wall distance
     else
        delta_aux = dexlo_ker                                     ! fixed wall distance
     end if
  end if
  !
  ! Boundary velocity
  !
  if (kfl_wlaav_ker==0) then
     if ( kfl_waexl_ker == 0_ip ) then 
        vewal(1:ndime) = gbvel(1:ndime)
     else ! using temp and veloc from exchange location
        vewal(1:ndime) = velex(1:ndime)
     end if
     !
     ! Tangent velocity
     !
     tveno = 0.0_rp
     do idime =1, ndime
        tvelo(idime) = vewal(idime)
        do ldime = 1,ndime
           tvelo(idime)= tvelo(idime)   &
                - baloc(idime,ndime) &
                * baloc(ldime,ndime) *vewal(ldime)
        end do
        tveno = tveno+tvelo(idime)*tvelo(idime)
     end do
  else ! average velocity for LES simulations
     call vecnor(veave,ndime,tveno,2_ip)        
  end if
     
  if(tveno>1.0e-12_rp) tveno=sqrt(tveno)

  vonka=0.41_rp 
  vikin=acvis/acden 

  call frivel(kfl_ustar,delta_aux,rough,tveno,vikin,velfr)             ! U*
 
  
  yplus = delta_aux*velfr/vikin 

  if (kfl_ustar >0  ) then ! ABL wall law
     
     if (kfl_ustar==2) & !redefine ustar in terms of tke
          velfr = kinen**0.5_rp*cmu_st**0.25_rp 
     ! convective factor
     
     hconv = acden*acsph*velfr*vonka/(prtur*log(1.0_rp + delta_aux/rough)) 
  else  if(yplus<1.0_rp) then
     Pr    =  acsph*acvis/actco
     ! convective factor
     hconv = vikin/(Pr*delta_aux) 
     
  else
     if(velfr>1.0e-12_rp) then
        uplus = tveno/velfr
     else
        uplus = 0.0_rp
     end if
     Pr    = acsph*acvis/actco
     Rt    = Pr/prtur
     if(Rt<0.1_rp) then
        alpha = 1.100_rp
     else
        alpha = 0.333_rp
     end if
     Dt    = exp(-yplus/11.0_rp*Rt**alpha)
     Pt    = 9.24_rp*(Rt**0.75_rp-1.0_rp)*(1.0_rp+0.28_rp*exp(-0.007_rp*Rt))
     Tplus = Pr*uplus*Dt+prtur*(1.0_rp-Dt)*(uplus+Pt)
     ! convective factor
     if(Tplus>1.0e-12_rp) then
        hconv  = acden*acsph*velfr/Tplus
     else
        hconv  = 0.0_rp
     end if
  end if


  arobi=-hconv 
  trobi= twall
  qrobi= 0.0_rp
  ustar= velfr
end subroutine tem_bouwal
