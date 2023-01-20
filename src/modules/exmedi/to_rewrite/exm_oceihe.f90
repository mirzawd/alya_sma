!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!subroutine exm_oceihe(n)
  !-----------------------------------------------------------------------
  !****f* Exmedi/exm_ceiHET  ***Updated by JAS FOR TT HETEROGENEOUS
  ! NAME 
  !    exm_ceicur
  ! DESCRIPTION
  !   This routine computes cell ionic currents of Ten Tusscher model (TT). 
  !
  !   [pA/pF] --> I_Na, I_CaL, I_to, I_Ks, I_Kr, I_K1, I_NaCa, I_NaK
  !   [pA/pF] --> I_pCa, I_pK, I_bNa, I_bCa, I_stim, I_ax
  !   [mM/ms] --> I_leak, I_up, I_rel, I_xfer
  !
  ! USES
  !   def_parame
  !   def_master
  !   def_domain
  !   def_exmedi
  ! USED BY
  !    exm_doiter
  !***
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  ! The variables are 17:
  !   - xina:   I_Na    -->  viclo_exm(1,1)  -->  fast current of sodium (Page H1575 TT paper)
  !   - xical:  I_CaL   -->  viclo_exm(2,1)  -->  current of L-Type calcium (Page H1576 TT paper)
  !   - xito:   I_to    -->  viclo_exm(3,1)  -->  transient outward current (Page H1577 TT paper)
  !   - xiks:   I_Ks    -->  viclo_exm(4,1)  -->  slow delayed rectifier current (Page H1578 TT paper)
  !   - xikr:   I_Kr    -->  viclo_exm(5,1)  -->  rapid delayed rectifier current (Page H1579 TT paper)
  !   - xik1:   I_K1    -->  viclo_exm(6,1)  -->  inward rectifier K current (Page H1580 TT paper)
  !   - xinaca: I_NaCa  -->  viclo_exm(7,1)  -->  sodium/calcium exchanger current (Page H1580 TT paper)
  !   - xinak:  I_NaK   -->  viclo_exm(8,1)  -->  sodium/potassium pump current (Page H1580 TT paper)
  !   - xipca:  I_pCa   -->  viclo_exm(9,1)  -->  calcium plateau current (Page H1580 TT paper)
  !   - xipk:   I_pK    -->  viclo_exm(10,1) -->  potassium plateau current (Page H1580 TT paper)
  !   - xibna:  I_bNa   -->  viclo_exm(11,1) -->  sodium background current (Page H1580 TT paper)
  !   - xibca:  I_bCa   -->  viclo_exm(12,1) -->  calcium background current (Page H1581 TT paper)
  !   - xileak: I_leak  -->  viclo_exm(13,1) -->  leak from SR to the cytoplasm (Page H1581 TT paper)
  !   - xiup:   I_up    -->  viclo_exm(14,1) -->  pump taking up Ca in the SR (Page H1581 TT paper)
  !   - xirel:  I_rel   -->  viclo_exm(15,1) -->  Ca-induced Ca release (CICR) (Page H1581 TT paper)
  !   - xistim: I_stim  -->  viclo_exm(16,1) -->  external stimulus current (Page H1581 TT paper)
  !   - xiax:   I_ax    -->  viclo_exm(17,1) -->  axial current flow (Page H1581 TT paper)
  !   - xixfer: I_xfer  -->  viclo_exm(18,1) -->  diffusive Ca current between diadic subspace and bulk cytoplasm (Page H1581 TT paper)
  ! The order is important because it is related with the order of nicel_exm variables (exm_memall.f90):
  !
  !-----------------------------------------------------------------------

!  use      def_parame
!  use      def_master
!  use      def_elmtyp
!  use      def_domain
!  use      def_exmedi

  ! DEFINITION OF VARIABLES
!  implicit none
!  integer(ip), intent(in) :: n
  !integer(ip), dimension(1), intent(in) :: nmate
  !integer(ip) :: nbeat, kfl_odets_exm, t, tt, i, j, itim
  !real(rp)    :: xina, xical, xito, xiks, xikr, xik1, xinaca, xinak, xipca, xipk
  !real(rp)    :: xibna, xibca, xileak, xiup, xirel, xiax
  !real(rp)    :: enapot, ecapot, ekpot, ekspot 
  !real(rp)    :: groupcon, xgroupcon
  !real(rp)    :: vaux1, vaux2, vaux3, vaux4, vaux5
  !real(rp)   ::  celltem, rgascon, faracon, xtrkcon, xtrnacon, xtrcacon
  !real(rp)   ::  nagmax, dosis1, dosis2, dosis3, dosis4, i50cal, i50kr, i50na, i50ito
  !real(rp)   ::  calgmax, gcaldrug, gkrdrug, gnadrug, gitodrug
  !real(rp)   ::  togmaxmce, togmaxepi, togmaxend
  !real(rp)   ::  ksgmaxend, ksgmaxepi, ksgmaxmce, ksnaperm
  !real(rp)   ::  krgmax
  !real(rp)   ::  k1gmax
  !real(rp)   ::  nacamax, gamanaca, kmcanaca, kmnanaca, satunaca, alfanaca
  !real(rp)   ::  nakmax, kmknak, kmnanak
  !real(rp)   ::  pcagmax, kpcapca
  !real(rp)   ::  pkgmax
  !real(rp)   ::  bnagmax
  !real(rp)   ::  bcagmax
  !real(rp)   ::  upvmax, upkcon
  !real(rp)   ::  leakvmax
  !real(rp)   ::  bufcyto, bufsr, bufss, kbufcyto, kbufsr, kbufss
  !real(rp)    :: caitot, casrto, naitot, kitot, cassto, rprime  
  !real(rp)   ::  volcyto, volsr, capacit, volss
  !real(rp)   ::  faradco, dtimon, k1, k2, k3, k4
  !real(rp)    :: thold, thnew, atim, batim, dur, aux1, val0, rhsx 
  !real(rp)    :: vinf, alphx, alph1, alph2, betax, beta1, beta2, gammax, taux, tau1, tau2, tau3, xitaux, vaux0


!  call runend("EXM_OCEIHE: exm_oceihe currently is not working correctly!")
!  return 

!   !n=nodemat(ipoin)
! !!!!  DRUG DEFINITION TO IKR, INA OR ICAL
!  if(kfl_drugsmate_exm(n)==1_ip)then
!     dosis1 = drugdmate_exm(1,n)
!     i50cal = drugdmate_exm(2,n)
!     dosis2 = drugdmate_exm(3,n)
!     i50kr = drugdmate_exm(4,n)
!     dosis3 = drugdmate_exm(5,n)
!     i50na = drugdmate_exm(6,n)
!     dosis4 = drugdmate_exm(7,n)
!     i50ito = drugdmate_exm(8,n)
!     gcaldrug = 1.0_rp / (1.0_rp + (dosis1/i50cal))
!     gkrdrug = 1.0_rp / (1.0_rp + (dosis2/i50kr))
!     gnadrug = 1.0_rp / (1.0_rp + (dosis3/i50na))
!     gitodrug = 1.0_rp / (1.0_rp + (dosis4/i50ito))
!  else
!     gcaldrug = 1.0_rp
!     gkrdrug = 1.0_rp
!     gnadrug = 1.0_rp
!     gitodrug = 1.0_rp
!  endif 

!   ! CONSTANTS DEFINITION  
!   rgascon = 8314.472_rp !ttparmate_exm(1)%value(1)      ! R [mJ/(K*mol)] --> ideal gas constant
!   celltem = 310.0_rp !ttparmate_exm(1)%value(2)      ! T [K] --> cell system temperature 
!   faracon = 96485.3415_rp !ttparmate_exm(1)%value(3)      ! F [C/mol] --> faraday constant
!   xtrkcon = 5.4_rp !ttparmate_exm(1)%value(4)      ! K_o [mM] --> extracellular concentration of potassium (K+)
!   xtrnacon = 140.0_rp !ttparmate_exm(1)%value(5)     ! Na_o [mM] --> extracellular concentration of sodium (Na+)
!   xtrcacon = 2.0_rp !ttparmate_exm(1)%value(6)     ! Ca_o [mM] --> extracellular concentration of calcium (Ca++)
!   volcyto = 0.016404_rp !ttparmate_exm(1)%value(7)      ! V_c [nL] --> cytoplasmic volume
!   volsr = 0.001094_rp * ttparmate_exm(1,5,n)        ! V_SR [nL] --> sarcoplasmic reticulum (SR) volume
!   capacit = 0.185_rp !ttparmate_exm(1)%value(9)      ! C [nF] --> Capacitance used in Ca_i_tot, Na_i_toy and K_i_tot
!   volss = 0.00005468_rp !ttparmate_exm(1)%value(10)        ! V_SS
!   faradco = 96485.3415_rp !ttparmate_exm(1)%value(3)      ! F [C/mol] --> faraday constant  

!   nagmax = 14.838_rp * gnadrug       ! G_Na [nS/pF] --> maximal conductance of I_Na current

!   calgmax = 0.0000398_rp * gcaldrug      ! G_CaL [L/(mF*s)]--> maximal conductance of I_CaL current

!   togmaxmce = 0.294_rp * ttparmate_exm(1,1,n)    ! G_to_M [nS/pF] --> maximal conductance of I_to current for M Cell
!   togmaxepi = 0.294_rp * ttparmate_exm(1,1,n)    ! G_to_epi [nS/pF] --> max conductance of I_to current for epicardium
!   togmaxend = 0.073_rp * ttparmate_exm(1,1,n)   ! G_to_end [nS/pF] --> max conductance of I_to current for endocardium

!   ksgmaxmce = 0.098_rp * ttparmate_exm(1,2,n)    ! G_Ks_M [nS/pF] --> max conductance of I_Ks current for M Cell
!   ksgmaxepi = 0.392_rp * ttparmate_exm(1,2,n)     ! G_Ks_epi [nS/pF] --> max conductance of I_Ks current for epicardium
!   ksgmaxend = 0.392_rp * ttparmate_exm(1,2,n)     ! G_Ks_end [nS/pF] --> max conductance of I_Ks current for endocardium
!   ksnaperm = 0.03_rp !ttparmate_exm(5)%value(4)     ! p_KNa [adim] --> relative I_Ks current permeability to Na+

!   krgmax = 0.153_rp * gkrdrug       ! G_Kr [nS/pF] --> max conductance of I_Kr current

!   k1gmax = 5.405_rp * ttparmate_exm(1,3,n)       ! G_K1 [nS/pF] --> max conductance of I_K1 current

!   nacamax = 1000.0_rp !ttparmate_exm(8)%value(1)      ! k_NaCa [pA/pF] --> max I_NaCa current
!   gamanaca = 0.35_rp !ttparmate_exm(8)%value(2)     ! gama [adim] --> voltage dependent parameter of I_NaCa current
!   kmcanaca = 1.38_rp !ttparmate_exm(8)%value(3)     ! K_mCa [mM] --> Cai half-saturation constant for I_NaCa current
!   kmnanaca = 87.5_rp !ttparmate_exm(8)%value(4)     ! K_mNai [mM] --> Nai half-saturation constant for I_NaCa current
!   satunaca = 0.1_rp ! ttparmate_exm(8)%value(5)    ! k_sat [adim] --> saturation factor I_NaCa current
!   alfanaca = 2.5_rp !ttparmate_exm(8)%value(6)     ! alfa [adim] --> factor enhancing outward nature of I_NaCa current

!   nakmax = 2.724_rp !ttparmate_exm(9)%value(1)       ! P_NaK [pA/pF] --> max I_NaK current
!   kmknak = 1.0_rp ! ttparmate_exm(9)%value(2)       ! K_mK [mM] --> K_o half-saturation constant of I_NaK current
!   kmnanak = 40.0_rp !ttparmate_exm(9)%value(3)      ! K_mNa [mM] --> Nai half-saturation constant of I_NaK current

!   pcagmax = 0.1238_rp !ttparmate_exm(10)%value(1)     ! G_pCa [pA/pF] --> max conductance of I_pCa current
!   kpcapca = 0.0005_rp !ttparmate_exm(10)%value(2)     ! K_pCa [mM] --> Cai half-saturation constant of I_pCa current

!   pkgmax = 0.0146_rp !ttparmate_exm(11)%value(1)      ! G_pK [nS/pF] --> max conductance of I_pK current

!   bnagmax = 0.00029_rp !ttparmate_exm(12)%value(1)     ! G_bNa [nS/pF] --> max conductance of I_bNa current

!   bcagmax = 0.000592_rp !ttparmate_exm(13)%value(1)     ! G_bCa [nS/pF] --> max conductance of I_bCa current

!   upvmax = 0.006375_rp !ttparmate_exm(14)%value(1)      ! V_maxup [mM/ms] --> max I_up current
!   upkcon = 0.00025_rp !ttparmate_exm(14)%value(2)      ! K_up [mM] --> half-saturation constant of I_up current

!   aux1 = ttparmate_exm(1,4,n)
!   leakvmax = 0.00036_rp * aux1    ! V_leak [1/ms] --> max I_leak current

!   bufcyto = 0.2_rp !ttparmate_exm(16)%value(1)     ! Buf_c [mM] --> total cytoplasmic buffer concentration
!   kbufcyto = 0.001_rp !ttparmate_exm(16)%value(2)    ! K_Bufc [mM] --> Cai half-saturation constant for cytoplasmic buffer
!   bufsr = 10.0_rp !ttparmate_exm(16)%value(3)       ! Buf_sr [mM] --> total sarcoplasmic buffer concentration
!   kbufsr = 0.3_rp !ttparmate_exm(16)%value(4)      ! K_Bufsr [mM] --> Ca_SR half-saturation constant for sarcoplasmic buffer
!   bufss = 0.4_rp !ttparmate_exm(16)%value(5)       ! Buf_sr [mM] --> total sarcoplasmic buffer concentration
!   kbufss = 0.00025_rp !ttparmate_exm(16)%value(6)      ! K_Bufsr [mM] --> Ca_SR half-saturation constant for sarcoplasmic buffer

!   groupcon = rgascon * celltem / faracon          ! Group of constants: R * T / F  [mV] 
!   xgroupcon = 1.0_rp / groupcon

!   dtimon = 0.001_rp

!   !if(INOTSLAVE) then

!   thnew = 0.5_rp            
!   thold = 1.0_rp - thnew       
!   kfl_odets_exm = 1
!   nbeat = 0
!   itim = 0
!   batim = 0
!   t = 1
!   tt = 0
!   j = 0

!   atim = 0.0_rp

!   do while (batim < (moneclmate_exm(1,n) * moneclmate_exm(2,n)))
!      aux1 = 0.0_rp
!      atim = dtimon * aux1
!      do while (atim < moneclmate_exm(2,n))
!         !if (atim < moneclmate_exm(2) )then            
!         aux1 = aux1 + 1
!         atim = dtimon * aux1
!         dur = 1.0_rp ! duration time (ms)
!         t = t + 1
!         tt = t - 1
!         if (atim < dur) then 
!            ! External stimulus current I_stim [pA/pF]
!            viclo_exm(16,t) = -52.0_rp
!            viclo_exm(16,tt) = -52.0_rp     ! Value of I_stim current [pA/pF]
!         else
!            viclo_exm(16,t) = 0.0_rp
!            viclo_exm(16,tt) = 0.0_rp
!         end if

!         ! Fast current of sodium I_Na [pA/pF]
!         ! nicel_exm = 1
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1585 (col 2)
!         enapot = groupcon * log (xtrnacon / vcolo_exm(3,tt))    ! Reversal potential of Na [mV]
!         vaux1 = vaulo_exm(1,tt) * vaulo_exm(1,tt) * vaulo_exm(1,tt)
!         vaux1 = nagmax * vaux1 * vaulo_exm(2,tt) * vaulo_exm(3,tt)
!         vaux2 = elmlo_exm(tt) - enapot
!         xina = vaux1 * vaux2  ! Variable xina

!         viclo_exm(1,t) = xina       ! Value of I_Na current [pA/pF]


!         ! Current of L-Type calcium I_CaL   [pA/pF]
!         ! nicel_exm = 2
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1586 (col 1)
!         !i_CaL = ((( g_CaL*d*f*f2*fCass*4*(V - 15)*pow(F, 2))/(R*T))*(0.25_rp*Ca_ss*(exp((( 2*(V - 15)*F)/(R*T)))) - Ca_o)) /
!         !        ((exp((( 2*(V - 15)*F)/(R*T)))) - 1)

!         vaux1 = calgmax * vaulo_exm(4,tt) * vaulo_exm (5,tt) * vaulo_exm(12,tt)
!         vaux1 = vaux1 * vaulo_exm(6,tt) * 4.0_rp * faracon * faracon * (elmlo_exm(tt) - 15.0_rp)
!         vaux1 = vaux1 / (rgascon * celltem)
!         vaux2 = exp((2.0_rp * (elmlo_exm(tt) - 15.0_rp) * faracon) / (rgascon * celltem))
!         vaux3 = (0.25_rp * vcolo_exm(5,tt) * vaux2) - xtrcacon
!         vaux3 = vaux1 * vaux3 
!         xical = vaux3 / (vaux2 - 1.0_rp) ! Variable xical

!         viclo_exm(2,t) = xical       ! Value of I_CaL current [pA/pF]

!         ! Transient outward current I_to [pA/pF]
!         ! nicel_exm = 3
!         ! AND Slow delayed rectifier current I_Ks [pA/pF]
!         ! nicel_exm = 4
!         !    i_to = g_to*r*s*(V - E_K)
!         !    i_Ks = g_Ks*(pow(Xs,2))*(V - E_Ks)
!         ekpot = groupcon * log(xtrkcon / vcolo_exm(4,tt))    ! Reversal potential of K [mV]

!         vaux1 = xtrkcon + (ksnaperm * xtrnacon)
!         vaux2 = vcolo_exm(4,tt) + (ksnaperm * vcolo_exm(3,tt))               
!         ekspot = groupcon * log( vaux1 / vaux2 )    ! Reversal potential of Ks [mV]
!         ! Variable xito
!         vaux1 =  vaulo_exm(7,tt) * vaulo_exm(8,tt)!
!         vaux1 = vaux1 * (elmlo_exm(tt) - ekpot)
!         ! Variable xiks                    
!         vaux2 = vaulo_exm(9,tt) * vaulo_exm(9,tt)
!         vaux2 = vaux2 * (elmlo_exm(tt) - ekspot)

!         if(ituss_exm == 3) then   !epicardial    
!            xito = vaux1 * togmaxepi
!            xiks = vaux2 * ksgmaxepi
!         else if(ituss_exm == 1) then !endocardial
!            xito = vaux1 * togmaxend
!            xiks = vaux2 * ksgmaxend
!         else if(ituss_exm == 2) then  !mid
!            xito = vaux1 * togmaxmce
!            xiks = vaux2 * ksgmaxmce
!         end if
!         viclo_exm(3,t) = xito  ! Variable xIto
!         viclo_exm(4,t) = xiks  ! Variable xIks

!         ! Rapid delayed rectifier current I_Kr [pA/pF]
!         ! nicel_exm = 5
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1586 (col 2)

!         vaux1 = krgmax * sqrt(xtrkcon / 5.4_rp) * vaulo_exm(10,tt) * vaulo_exm(11,tt)
!         xikr = vaux1 * (elmlo_exm(tt) - ekpot)
!         viclo_exm(5,t) = xikr  ! Value of I_Kr current [pA/pF]

!         ! Inward rectifier potassium current I_K1 [pA/pF]
!         ! nicel_exm = 6
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1)
!         !  beta_K1 = ( 3*exp( 0.0002_rp*(V - E_K+100))+(exp(( 0.1_rp*((V - E_K) - 10))))) /
!         !   (1+(exp((-0.5_rp*(V - E_K)))))
!         vaux1 = 3.0_rp * exp(0.0002_rp * (100.0_rp + elmlo_exm(tt) - ekpot))
!         vaux1 = vaux1 + exp(0.1_rp * (-10.0_rp + elmlo_exm(tt) - ekpot))
!         vaux2 = 1.0_rp + exp(-0.5_rp * (elmlo_exm(tt) - ekpot))
!         vaux2 = 1.0_rp / vaux2
!         vaux1 = vaux1 * vaux2                  ! Value of beta_K1
!         vaux2 = 1.0_rp + exp(0.06_rp * (elmlo_exm(tt) - ekpot - 200.0_rp))
!         vaux2 = 0.1_rp / vaux2                 ! Value of alpha_K1
!         vaux3 = vaux2 / (vaux1 + vaux2)        ! Value of K1_inf
!         !   i_K1 = g_K1*xK1_inf* pow((K_o/5.4_rp), (1/2))*(V - E_K)     
!         vaux1 = sqrt(xtrkcon / 5.4_rp) * vaux3 * (elmlo_exm(tt) - ekpot)
!         xik1 = k1gmax * vaux1
!         viclo_exm(6,t) = xik1       ! Value of I_K1 current [pA/pF]

!         ! Sodium/calcium exchanger current I_NaCa [pA/pF]
!         ! nicel_exm = 7
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1) 
!         !     xinaca = ( K_NaCa.*( (exp((( gamma_NaCa.*V.*F)./( R.*T)))).*(Nai .^ 3.0).*Cao -  (exp((( (gamma_NaCa - 1.0).*V.*F)./( R.*T)))).
!         !   *(Nao .^ 3.0).*Cai.*alpha_NaCa))./( ((Km_Nai .^ 3.0)+(Nao .^ 3.0)).*(Km_Ca+Cao).*(1.0+ K_sat.*(exp((( (gamma_NaCa - 1.0).*V.*F)./( R.*T))))));

!         vaux1 = exp((gamanaca-1.0_rp) * elmlo_exm(tt) * xgroupcon)
!         vaux2 = (kmnanaca * kmnanaca * kmnanaca) + (xtrnacon * xtrnacon * xtrnacon)
!         vaux2 = vaux2 * (kmcanaca + xtrcacon)
!         vaux2 = vaux2 * (1.0_rp + (satunaca * vaux1))     ! current denominator 
!         vaux2 = 1.0_rp / vaux2                          ! inverse of current denominator

!         vaux3 = -vaux1 * xtrnacon * xtrnacon * xtrnacon * vcolo_exm(1,tt) * alfanaca         
!         vaux1 = exp(gamanaca * elmlo_exm(tt) * xgroupcon)
!         vaux3 = vaux3 + (vaux1 * vcolo_exm(3,tt) * vcolo_exm(3,tt) * vcolo_exm(3,tt) * xtrcacon)          ! current numerator

!         xinaca = nacamax * vaux3 * vaux2 

!         viclo_exm(7,t) = xinaca       ! Value of I_NaCa current [pA/pF]

!         ! Sodium/potassium pump current I_NaK [pA/pF]
!         ! nicel_exm = 8
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1)
!         !i_NaK = (((( P_NaK*K_o)/(K_o+K_mk))*Na_i)/(Na_i+K_mNa))/
!         !        (1+ 0.1245*(exp(((-0.1*V*F)/( R*T))))+ 0.0353*(exp(((-V*F)/( R*T))))) 
!         vaux1 = ((nakmax * xtrkcon) / (xtrkcon + kmknak)) * vcolo_exm(3,tt) 
!         vaux1 = vaux1 / (vcolo_exm(3,tt) + kmnanak)             ! current denominator 
!         vaux2 = 1.0_rp + 0.1245_rp * exp(-0.1_rp * elmlo_exm(tt) * xgroupcon)
!         vaux2 = vaux2 + 0.0353_rp * exp(-elmlo_exm(tt) * xgroupcon)
!         xinak =  vaux1 / vaux2 

!         viclo_exm(8,t) = xinak       ! Value of I_NaK current [pA/pF]

!         ! Calcium plateau current I_pCa [pA/pF]
!         ! nicel_exm = 9
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1)

!         vaux1 = kpcapca + vcolo_exm(1,tt)
!         vaux1 = 1.0_rp / vaux1             ! inverse of current denominator
!         xipca = pcagmax * vaux1 * vcolo_exm(1,tt)

!         viclo_exm(9,t) = xipca       ! Value of I_pCa current [pA/pF]

!         ! Potassium plateau current I_pK [pA/pF]
!         ! nicel_exm = 10
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1)
!         !    i_p_K = (g_pK*(V - E_K))/(1+(exp(((25 - V)/5.98_rp))))
!         vaux1 = 1.0_rp + exp((25.0_rp - elmlo_exm(tt)) / 5.98_rp)
!         vaux1 = 1.0_rp / vaux1             ! inverse of current denominator
!         xipk = pkgmax * vaux1 * (elmlo_exm(tt) - ekpot)

!         viclo_exm(10,t) = xipk       ! Value of I_pK current [pA/pF]

!         ! Sodium background current I_bNa [pA/pF]
!         ! nicel_exm = 11
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1)

!         xibna = bnagmax * (elmlo_exm(tt) - enapot)

!         viclo_exm(11,t) = xibna      ! Value of I_bNa current [pA/pF]

!         ! Calcium background current I_bCa [pA/pF]
!         ! nicel_exm = 12
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1)

!         ecapot = 0.50_rp * groupcon  * log (xtrcacon / vcolo_exm(1,tt))    ! Reversal potential of Ca [mV]
!         xibca = bcagmax * (elmlo_exm(tt) - ecapot)

!         viclo_exm(12,t) = xibca      ! Value of I_bCa current [pA/pF]

!         ! Leak current from sarcoplasmic reticulum (SR) to the cytoplasm I_leak [mM/ms]
!         ! nicel_exm = 13
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1)

!         xileak = leakvmax * (vcolo_exm(2,tt) - vcolo_exm(1,tt))

!         viclo_exm(13,t) = xileak      ! Value of I_leak current [mM/ms]

!         ! Pump current taking up calcium in the sarcoplasmic reticulum (SR) I_up [mM/ms]
!         ! nicel_exm = 14
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1)
!         ! i_up = Vmax_up/(1+pow(K_up, 2)/pow(Ca_i, 2)) 
!         vaux1 = vcolo_exm(1,tt) * vcolo_exm(1,tt)
!         vaux1 = 1.0_rp / vaux1
!         vaux2 = 1.0_rp + (upkcon * upkcon * vaux1)
!         vaux2 = 1.0_rp / vaux2            ! inverse of current denominator
!         xiup = upvmax * vaux2

!         viclo_exm(14,t) = xiup      ! Value of I_up current [mM/ms]

!         ! Calcium-induced calcium release (CICR) current I_rel [mM/ms]
!         ! nicel_exm = 15

!         ! Variable of activation O (related with I_rel current)
!         ! nconc = 10
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1)
!         !  kcasr = max_sr - (max_sr - min_sr)/(1+(pow((EC/Ca_SR), 2)))

!         vaux1 = 1.0_rp + ((1.5_rp / vcolo_exm(2,tt))*(1.5_rp / vcolo_exm(2,tt)))  !conc2 = CaSRfree = CaSR
!         vaux1 = 1.0_rp / vaux1
!         vaux1 = 2.5_rp - (1.5_rp * vaux1)  !kCaSR
!         !     O = ( k1*pow(Ca_ss, 2)*R_prime)/(k3+ k1*pow(Ca_ss, 2))      
!         vaux3 = 0.15_rp / vaux1   !K1
!         vaux2 = 0.045_rp * vaux1  !K2


!         rhsx = 1.0_rp / (0.06_rp + (vaux3*vcolo_exm(5,tt)*vcolo_exm(5,tt)))      !O for calcium dynamics
!         rhsx = vaux3 * vcolo_exm(5,tt) * vcolo_exm(5,tt) * vcolo_exm(9,tt) * rhsx !O
!         vcolo_exm(10,t) = rhsx
!         ! Variable of activation O (related with I_rel current)
!         ! nconc= 9
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1)
!         !  dR_prime_dt = - k2*Ca_ss*R_prime+ k4*(1 - R_prime)

!         vaux3 = - vaux2 * vcolo_exm(5,tt) * vcolo_exm(9,tt) 
!         rhsx = vaux3 + (0.005_rp * (1.0_rp - vcolo_exm(9,tt)))

!         val0 = vcolo_exm(9,tt)    ! Value of Rprime in previous dtimon step 

!         k1 = rhsx 
!         k2 = rhsx + 0.5_rp * dtimon* k1
!         k3 = rhsx + 0.5_rp * dtimon * k2 
!         k4 = rhsx + dtimon * k3 
!         rprime = val0 + ((dtimon / (6.0_rp)) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4))
! !!!!rprime = val0 + dtimon/capacit * rhsx
!         vcolo_exm(9,t) = rprime

!         ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1)

!         xirel = 0.102_rp * vcolo_exm(10,t) * (vcolo_exm(2,tt) - vcolo_exm(5,tt))  

!         viclo_exm(15,t) = xirel      ! Value of I_rel current [mM/ms]

!         ! IMPORTANTE: ESTO (I_ax) EN REALIDAD NO SE DEBERIA DEFINIR ACA
!         !             HABRIA QUE INTRODUCIRLO DESDE AFUERA
!         !             SE USA PARA EL CALCULO DE LAS CONCENTRACIONES

!         ! Axial current flow I_ax [pA/pF]
!         ! nicel_exm = 17
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1581

!         xiax = 0.0_rp

!         viclo_exm(17,t) = xiax      ! Value of I_ax current [pA/pF]

!         !Variable xixfer
!         viclo_exm(18,t) = 0.0038_rp * (vcolo_exm(5,tt) - vcolo_exm(1,tt))

!         !!  VOLTAGE INTEGRATION    
!         elmlo_exm(t) = 0.0_rp
!         do i=1,12
!            elmlo_exm(t)= elmlo_exm(t) + viclo_exm(i,t)
!         end do
!         elmlo_exm(t) =  -(elmlo_exm(t) + viclo_exm(16,t))
!         !!  Runge Kutta
!         k1 = elmlo_exm(t) 
!         k2 = elmlo_exm(t) + 0.5_rp * dtimon * k1
!         k3 = elmlo_exm(t) + 0.5_rp * dtimon * k2 
!         k4 = elmlo_exm(t) + dtimon * k3 
!         elmlo_exm(t) = elmlo_exm(tt) + ((dtimon/6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4))

!         !elmlo_exm(tt)=elmlo_exm(t)   

! !!!!!  CONCENTRATIONS!!!!            
!         ! Concentration of total calcium in the Cytoplasm (Ca_i_tot)
!         ! nconc_exm = 1
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 2)

!         !   Ca_i_bufc = 1/(1+( Buf_c*K_buf_c)/(pow((Ca_i+K_buf_c), 2)))
!         vaux1 = vcolo_exm(1,tt) + kbufcyto
!         vaux1 = 1.0_rp / (vaux1 * vaux1)
!         vaux1 = vaux1 * kbufcyto * bufcyto    !Variable caibuf
!         vaux1 = 1.0_rp / (1.0_rp + vaux1)
!         vcolo_exm(6,t) = vaux1

!         ! d_caitot/d_t = - (I_CaL + I_bCa + I_pCa - 2 * I_NaCa) * C / (2 * V_c * F) + I_leak - I_up + I_rel
!         ! dCa_i_dt = Ca_i_bufc*((( (i_leak - i_up)*V_sr)/V_c+i_xfer) - (1*((i_b_Ca+i_p_Ca) - 2*i_NaCa)*Cm)/( 2*1*V_c*F))

!         vaux1 = - 2.0_rp * volcyto * faradco
!         vaux1 = 1.0_rp / vaux1  

!         vaux3 = (viclo_exm(12,t) + viclo_exm(9,t)) - (2.0_rp *  viclo_exm(7,t) )
!         vaux3 = vaux1 * vaux3 * capacit
!         vaux2 = (viclo_exm(13,t) - viclo_exm(14,t)) * (volsr / volcyto) + viclo_exm(18,t)
!         rhsx =  vcolo_exm(6,t) * (vaux2 + vaux3) !

!         val0 = vcolo_exm(1,tt)

!         k1 = rhsx 
!         k2 = rhsx + 0.5_rp * dtimon * k1
!         k3 = rhsx + 0.5_rp * dtimon * k2 
!         k4 = rhsx + dtimon* k3 
!         caitot = val0 + ((dtimon/(6.0_rp)) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4))

! !!!!!!caitot = val0 + dtimon/capacit * rhsx
!         !caitot =  val0 + (0.5_rp * dtimon * (k1 + k2))
!         vcolo_exm(1,t) = caitot       ! Value of Ca_i_tot concentration

!         ! Concentration of total calcium in the Sarcoplasmic Reticulum (Ca_sr_tot)
!         ! nconc_exm = 2
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 2)

!         vaux1 = vcolo_exm(2,tt)  + kbufsr
!         vaux1 = 1.0_rp / (vaux1 * vaux1)
!         vaux1 = vaux1 * kbufsr * bufsr
!         vaux1 = 1.0_rp / (1.0_rp + vaux1)
!         vcolo_exm(7,t) = vaux1

!         ! d_casrto/d_t = Ca_sr_bufsr*(i_up - (i_rel+i_leak))

!         vaux1 = viclo_exm(14,t) - viclo_exm(13,t) - viclo_exm(15,t)

!         rhsx = vaux1 * vcolo_exm(7,t)  !capacit *

!         val0 = vcolo_exm(2,tt)
!         k1 = rhsx 
!         k2 = rhsx + 0.5_rp * dtimon * k1
!         k3 = rhsx + 0.5_rp * dtimon * k2 
!         k4 = rhsx + dtimon * k3 
!         casrto = val0 + ((dtimon / ( 6.0_rp)) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4))

! !!!!casrto = val0 + dtimon/capacit * rhsx
!         vcolo_exm(2,t) = casrto       ! Value of Ca_sr_tot concentration

!         ! Concentration of total calcium in the diadic space (Ca_ss_tot)
!         ! nconc_exm = 5
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 2)
!         vaux1 = vcolo_exm(5,tt)  + kbufss
!         vaux1 = 1.0_rp / (vaux1 * vaux1)
!         vaux1 = vaux1 * kbufss * bufss
!         vaux1 = 1.0_rp / (1.0_rp + vaux1)
!         vcolo_exm(8,t) = vaux1                    !Variable cassbufss

!         ! d_cass/d_t = Ca_ss_bufss*(((-1*i_CaL)/( 2*1*V_ss*F)+( i_rel*V_sr)/V_ss) - ( i_xfer*V_c)/V_ss)

!         vaux1 = 1.0_rp / (2.0_rp * faradco * volss)
!         vaux1 = vaux1 * viclo_exm(2,t) * capacit
!         vaux2 = viclo_exm(15,t) * volsr / volss
!         vaux3 = viclo_exm(18,t) * volcyto / volss

!         rhsx =  vcolo_exm(8,t) * (- vaux1 + vaux2 - vaux3)

!         val0 = vcolo_exm(5,tt)

!         k1 = rhsx 
!         k2 = rhsx + 0.5_rp * dtimon * k1
!         k3 = rhsx + 0.5_rp * dtimon* k2 
!         k4 = rhsx + dtimon * k3 
!         cassto = val0 + ((dtimon / ( 6.0_rp)) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4))
! !!!!cassto = val0 + dtimon/capacit * rhsx
!         vcolo_exm(5,t) = cassto       ! Value of Ca_ss_tot concentration

!         ! Intracellular sodium concentration (Na_i_tot)
!         ! nconc_exm = 3
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 2)

!         ! d_naitot/d_t =  ((- 1*(i_Na+i_b_Na+ 3*i_NaK+ 3*i_NaCa))/( 1*V_c*F))*Cm

!         vaux1 = volcyto * faradco
!         vaux1 = 1.0_rp / vaux1 !* capacit

!         vaux3 = (3.0_rp * viclo_exm(8,t)) + (3.0_rp * viclo_exm(7,t)) 
!         vaux3 = vaux3 + viclo_exm(1,t) + viclo_exm(11,t)
!         rhsx = - capacit * vaux1 * vaux3

!         val0 = vcolo_exm(3,tt)

!         k1 = rhsx 
!         k2 = rhsx + 0.5_rp * dtimon * k1
!         k3 = rhsx + 0.5_rp * dtimon * k2 
!         k4 = rhsx + dtimon * k3 
!         naitot = val0 + ((dtimon / ( 6.0_rp)) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4))
! !!!!naitot = val0 + dtimon/xmccm_exm(nomat_exm(ipoin)) * rhsx
!         vcolo_exm(3,t) = naitot       ! Value of Na_i_tot concentration


!         ! Intracellular potassium concentration (K_i_tot)
!         ! nconc_exm = 4
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 2)

!         ! d_kitot/d_t = ((- 1*((i_K1+i_to+i_Kr+i_Ks+i_p_K+i_Stim) -  2*i_NaK))/( 1*V_c*F))*Cm

!         vaux1 = - volcyto * faradco
!         vaux1 = capacit / vaux1 !* capacit  

!         vaux3 = (viclo_exm(6,t) + viclo_exm(3,t) + viclo_exm(5,t) &
!              + viclo_exm(4,t) + viclo_exm(10,t) + viclo_exm(16,t)) - (2.0_rp * viclo_exm(8,t)) 

!         rhsx = vaux1 * vaux3

!         val0 = vcolo_exm(4,tt)

!         k1 = rhsx 
!         k2 = rhsx + 0.5_rp * dtimon * k1
!         k3 = rhsx + 0.5_rp * dtimon * k2 
!         k4 = rhsx + dtimon * k3 
!         kitot = val0 + ((dtimon / (6.0_rp)) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4))
! !!!!kitot = val0 + dtimon/capacit * rhsx
!         vcolo_exm(4,t) = kitot          ! Value of K_i_tot concentration

! !!!!! AUXILIARY:  GATES  VAUXI

!         vinf = 1.0_rp + exp((-56.86_rp - elmlo_exm(tt)) / 9.03_rp)
!         vinf = vinf * vinf
!         vinf = 1.0_rp / vinf              ! Value of m_infinite

!         alphx = 1.0_rp + exp((-60.00_rp - elmlo_exm(tt)) / 5.0_rp)
!         alphx = 1.0_rp / alphx            ! Value of alpha_m
!         !   beta_m = 0.1_rp/(1+(exp(((V+35)/5))))+0.1_rp/(1+(exp(((V - 50)/200))))
!         beta1 = 1.0_rp + exp(( 35.00_rp + elmlo_exm(tt)) / 5.0_rp)
!         beta1 = 0.1_rp / beta1
!         beta2 = 1.0_rp + exp((-50.00_rp + elmlo_exm(tt)) / 200.0_rp)
!         beta2 = 0.1_rp / beta2
!         betax = beta1 + beta2             ! Value of beta_m

!         taux  = betax * alphx             ! Value of tau_m
!         xitaux = 1.0_rp / taux            ! Value of 1/tau_m 

!         vaux0 = vaulo_exm(1,tt)      ! Value of m in previous time step 


!         ! Scheme for numerical integration
!         if (kfl_odets_exm == 1) then      ! Runge - Kutta 4

!            vaux1 = (vinf - vaux0) * xitaux;
!            vaux2 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux1)) * xitaux;
!            vaux3 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux2)) * xitaux;
!            vaux4 = (vinf - (vaux0 + dtimon * vaux3)) * xitaux;
!            vaux5 = vaux0 + (dtimon/ 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4); 
!            vaux3 = vaux5

!         else if (kfl_odets_exm == 2) then      ! Crank - Nicholson

!            vaux1 = dtimon * xitaux
!            vaux2 = 1.0_rp + thnew * vaux1
!            vaux2 = 1.0_rp / vaux2
!            vaux3 = (1.0_rp - thold * vaux1) * vaux0 + vaux1 * vinf
!            vaux3 = vaux3 * vaux2       

!         else if (kfl_odets_exm == 3) then      ! Rush - Larsen

!            vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/capacit * xitaux)    

!         end if

!         vaulo_exm(1,t) = vaux3      ! Value of variable m


!         ! Variable of activation h (related with I_Na current)
!         ! nauxi_exm = 2
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1586 (col 1)

!         vinf = 1.0_rp + exp((71.55_rp + elmlo_exm(tt)) / 7.43_rp)
!         vinf = vinf * vinf
!         vinf = 1.0_rp / vinf              ! Value of h_infinite

!         if (elmlo_exm(tt) >= -40.0_rp) then

!            alphx = 0.0_rp                     ! Value of alpha_h

!            betax = 0.13_rp * ( 1.0_rp + exp((-10.66_rp - elmlo_exm(tt)) / 11.1_rp) )
!            betax = 0.77_rp / betax         ! Value of beta_h

!         else

!            alphx = 0.057_rp * exp((-80.0_rp - elmlo_exm(tt)) / 6.8_rp)        ! Value of alpha_h
!            !   ( 2.7*(exp(( 0.079*V)))+ 310000*(exp((0.3485*V)))) 
!            beta1 = 2.7_rp * exp(0.079_rp * elmlo_exm(tt))
!            beta2 = 310000.0_rp * exp(0.3485_rp * elmlo_exm(tt))
!            betax = beta1 + beta2           ! Value of beta_h

!         end if

!         xitaux  = betax + alphx          ! Value of 1/tau_h

!         vaux0 = vaulo_exm(2,tt)      ! Value of h in previous time step 


!         ! Scheme for numerical integration
!         if (kfl_odets_exm == 1) then      ! Runge - Kutta 4

!            vaux1 = (vinf - vaux0) * xitaux;
!            vaux2 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux1)) * xitaux;
!            vaux3 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux2)) * xitaux;
!            vaux4 = (vinf - (vaux0 + dtimon * vaux3)) * xitaux;
!            vaux5 = vaux0 + (dtimon / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4); 
!            vaux3 = vaux5

!         else if (kfl_odets_exm == 2) then      ! Crank - Nicholson

!            vaux1 = dtimon * xitaux
!            vaux2 = 1.0_rp + thnew * vaux1
!            vaux2 = 1.0_rp / vaux2
!            vaux3 = (1.0_rp - thold * vaux1) * vaux0 + vaux1 * vinf
!            vaux3 = vaux3 * vaux2       

!         else if (kfl_odets_exm == 3) then      ! Rush - Larsen

!            vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/capacit * xitaux)    

!         end if
!         vaulo_exm(2,t) = vaux3      ! Value of variable h

!         ! Variable of activation j (related with I_Na current)
!         ! nauxi_exm = 3
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1586 (col 1)

!         vinf = 1.0_rp + exp((71.55_rp + elmlo_exm(tt)) / 7.43_rp)
!         vinf = vinf * vinf
!         vinf = 1.0_rp / vinf              ! Value of j_infinite

!         if (elmlo_exm(tt) >= -40.0_rp) then

!            alphx = 0                       ! Value of alpha_j

!            beta1 = 1.0_rp + exp(-0.1_rp * (32.0_rp + elmlo_exm(tt)))
!            beta1 = 1.0_rp / beta1
!            beta2 = 0.6_rp * exp(0.057_rp * elmlo_exm(tt))
!            betax = beta1 * beta2           ! Value of beta_j

!         else
!            !    (( ( (-25428*exp( 0.2444*V)) -  (6.948e-6*exp(-0.04391*V)) )*(V+37.78) ) / (1+exp( 0.311*(V+79.23)) ))
!            alph1 = 1.0_rp + exp(0.311_rp * (79.23_rp + elmlo_exm(tt)))
!            alph1 = 1.0_rp / alph1
!            alph2 = - 25428.0_rp * exp(0.2444_rp * elmlo_exm(tt))
!            alph2 = alph2 - (0.000006948_rp * exp(-0.04391_rp * elmlo_exm(tt)))
!            alph2 = alph2 * (elmlo_exm(tt) + 37.78_rp)
!            alphx = alph1 * alph2           ! Value of alpha_j

!            beta1 = 1.0_rp + exp(-0.1378_rp * (40.14_rp + elmlo_exm(tt)))
!            beta1 = 1.0_rp / beta1
!            beta2 = 0.02424_rp * exp(-0.01052_rp * elmlo_exm(tt))
!            betax = beta1 * beta2           ! Value of beta_j

!         end if

!         xitaux  = betax + alphx          ! Value of 1/tau_j

!         vaux0 = vaulo_exm(3,tt)      ! Value of j in previous time step 


!         ! Scheme for numerical integration
!         if (kfl_odets_exm == 1) then      ! Runge - Kutta 4

!            vaux1 = (vinf - vaux0) * xitaux;
!            vaux2 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux1)) * xitaux;
!            vaux3 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux2)) * xitaux;
!            vaux4 = (vinf - (vaux0 + dtimon * vaux3)) * xitaux;
!            vaux5 = vaux0 + (dtimon / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4); 
!            vaux3 = vaux5

!         else if (kfl_odets_exm == 2) then      ! Crank - Nicholson

!            vaux1 = dtimon * xitaux
!            vaux2 = 1.0_rp + thnew * vaux1
!            vaux2 = 1.0_rp / vaux2
!            vaux3 = (1.0_rp - thold * vaux1) * vaux0 + vaux1 * vinf
!            vaux3 = vaux3 * vaux2       

!         else if (kfl_odets_exm == 3) then      ! Rush - Larsen

!            vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/capacit * xitaux)    

!         end if


!         vaulo_exm(3,t) = vaux3      ! Value of variable j


!         ! Variable of activation d (related with I_CaL current)
!         ! nauxi_exm = 4
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1586 (col 1)

!         vinf = 1.0_rp + exp((-8.0_rp - elmlo_exm(tt)) / 7.5_rp)
!         vinf = 1.0_rp / vinf              ! Value of d_infinite

!         alphx = 1.0_rp + exp((-35.0_rp - elmlo_exm(tt)) / 13.0_rp)
!         alphx = (1.4_rp / alphx) + 0.25_rp  ! Value of alpha_d

!         betax = 1.0_rp + exp((5.0_rp + elmlo_exm(tt)) / 5.0_rp)
!         betax = 1.4_rp / betax            ! Value of beta_d

!         gammax = 1.0_rp + exp((50.0_rp - elmlo_exm(tt)) / 20.0_rp)
!         gammax = 1.0_rp / gammax          ! Value of gamma_d

!         taux = (alphx * betax) + gammax     ! Value of tau_d
!         xitaux = 1.0_rp / taux            ! Value of 1/tau_d

!         vaux0 = vaulo_exm(4,tt)      ! Value of d in previous time step 


!         ! Scheme for numerical integration
!         if (kfl_odets_exm == 1) then      ! Runge - Kutta 4

!            vaux1 = (vinf - vaux0) * xitaux;
!            vaux2 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux1)) * xitaux;
!            vaux3 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux2)) * xitaux;
!            vaux4 = (vinf - (vaux0 + dtimon * vaux3)) * xitaux;
!            vaux5 = vaux0 + (dtimon / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4); 
!            vaux3 = vaux5

!         else if (kfl_odets_exm == 2) then      ! Crank - Nicholson

!            vaux1 = dtimon * xitaux
!            vaux2 = 1.0_rp + thnew * vaux1
!            vaux2 = 1.0_rp / vaux2
!            vaux3 = (1.0_rp - thold * vaux1) * vaux0 + vaux1 * vinf
!            vaux3 = vaux3 * vaux2       

!         else if (kfl_odets_exm == 3) then      ! Rush - Larsen

!            vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/capacit * xitaux)    

!         end if

!         vaulo_exm(4,t) = vaux3      ! Value of variable d

!         ! Variable of activation f (related with I_CaL current)
!         ! nauxi_exm = 5
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1586 (col 1)

!         vinf = 1.0_rp + exp((20.0_rp + elmlo_exm(tt)) / 7.0_rp)
!         vinf = 1.0_rp / vinf              ! Value of f_infinite
!         ! tau_f = 1102.5*(exp(( - (pow((V+27),2))/225)))+200/(1+(exp(((13 - V)/10))))+180/(1+(exp(((V+30)/10))))+20
!         tau1 = 1102.5_rp * exp(-((27.0_rp + elmlo_exm(tt)) * (27.0_rp + elmlo_exm(tt))) / 225.0_rp)
!         tau2 = 1.0_rp + exp((13.0_rp - elmlo_exm(tt)) / 10.0_rp)
!         tau2 = 200.0_rp / tau2
!         tau3 = 1.0_rp + exp((30.0_rp + elmlo_exm(tt)) / 10.0_rp)
!         tau3 = 180.0_rp / tau3
!         taux = tau1 + tau2 + tau3 + 20.0_rp      ! Value of tau_f
!         xitaux = 1.0_rp / taux            ! Value of 1/tau_f

!         vaux0 = vaulo_exm(5,tt)      ! Value of f in previous time step 

!         ! Scheme for numerical integration
!         if (kfl_odets_exm == 1) then      ! Runge - Kutta 4

!            vaux1 = (vinf - vaux0) * xitaux;
!            vaux2 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux1)) * xitaux;
!            vaux3 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux2)) * xitaux;
!            vaux4 = (vinf - (vaux0 + dtimon * vaux3)) * xitaux;
!            vaux5 = vaux0 + (dtimon / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4); 
!            vaux3 = vaux5

!         else if (kfl_odets_exm == 2) then      ! Crank - Nicholson

!            vaux1 = dtimon * xitaux
!            vaux2 = 1.0_rp + thnew * vaux1
!            vaux2 = 1.0_rp / vaux2
!            vaux3 = (1.0_rp - thold * vaux1) * vaux0 + vaux1 * vinf
!            vaux3 = vaux3 * vaux2       

!         else if (kfl_odets_exm == 3) then      ! Rush - Larsen

!            vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/capacit * xitaux)    

!         end if

!         vaulo_exm(5,t) = vaux3      ! Value of variable f

!         ! Variable of activation f2 (related with I_CaL current)
!         ! nauxi_exm = 12
!         ! "A model for ventricular tissue", Ten Tusscher, 2006 paper
!         ! modified by JAS

!         vinf = 1.0_rp + exp((35.0_rp + elmlo_exm(tt)) / 7.0_rp)
!         vinf = (0.67_rp / vinf ) + 0.33_rp            ! Value of f2_infinite

!         tau1 = 600.0_rp * exp((-(25.0_rp + elmlo_exm(tt)) * (25.0_rp + elmlo_exm(tt))) / 170.0_rp)
!         tau2 = 1.0_rp + exp((25.0_rp - elmlo_exm(tt)) / 10.0_rp)
!         tau2 = 31.0_rp / tau2
!         tau3 = 1.0_rp + exp((30.0_rp + elmlo_exm(tt)) / 10.0_rp)
!         tau3 = 16.0_rp / tau3
!         taux = tau1 + tau2 + tau3      ! Value of tau_f2
!         xitaux = 1.0_rp / taux            ! Value of 1/tau_f2

!         vaux0 = vaulo_exm(12,tt)      ! Value of f2 in previous time step 


!         ! Scheme for numerical integration
!         if (kfl_odets_exm == 1) then      ! Runge - Kutta 4

!            vaux1 = (vinf - vaux0) * xitaux;
!            vaux2 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux1)) * xitaux;
!            vaux3 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux2)) * xitaux;
!            vaux4 = (vinf - (vaux0 + dtimon * vaux3)) * xitaux;
!            vaux5 = vaux0 + (dtimon / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4); 
!            vaux3 = vaux5

!         else if (kfl_odets_exm == 2) then      ! Crank - Nicholson

!            vaux1 = dtimon * xitaux
!            vaux2 = 1.0_rp + thnew * vaux1
!            vaux2 = 1.0_rp / vaux2
!            vaux3 = (1.0_rp - thold * vaux1) * vaux0 + vaux1 * vinf
!            vaux3 = vaux3 * vaux2       

!         else if (kfl_odets_exm == 3) then      ! Rush - Larsen

!            vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/capacit * xitaux)    

!         end if

!         vaulo_exm(12,t) = vaux3      ! Value of variable f2

!         ! Variable of activation fCass (related to I_CaL current)
!         ! nauxi_exm = 6
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1586 (col 2)

!         vinf = 1.0_rp + ((vcolo_exm(5,tt) / 0.05_rp)*(vcolo_exm(5,tt) / 0.05_rp))
!         vinf = (0.6_rp / vinf ) + 0.4_rp            ! Value of fCass_infinite

!         tau1 = 1.0_rp + ((vcolo_exm(5,tt) / 0.05_rp)*(vcolo_exm(5,tt) / 0.05_rp))
!         tau1 = (80.0_rp / tau1) + 2.0_rp           ! Value of tau_fCass    
!         xitaux = 1.0_rp / tau1            ! Value of 1/tau_f

!         vaux0 = vaulo_exm(6,tt)      ! Value of f in previous time step 

!         ! Scheme for numerical integration
!         if (kfl_odets_exm == 1) then      ! Runge - Kutta 4

!            vaux1 = (vinf - vaux0) * xitaux;
!            vaux2 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux1)) * xitaux;
!            vaux3 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux2)) * xitaux;
!            vaux4 = (vinf - (vaux0 + dtimon * vaux3)) * xitaux;
!            vaux5 = vaux0 + (dtimon / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4); 
!            vaux3 = vaux5

!         else if (kfl_odets_exm == 2) then      ! Crank - Nicholson

!            vaux1 = dtimon * xitaux
!            vaux2 = 1.0_rp + thnew * vaux1
!            vaux2 = 1.0_rp / vaux2
!            vaux3 = (1.0_rp - thold * vaux1) * vaux0 + vaux1 * vinf
!            vaux3 = vaux3 * vaux2       

!         else if (kfl_odets_exm == 3) then      ! Rush - Larsen

!            vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/capacit * xitaux)    

!         end if

!         vaulo_exm(6,t) = vaux3      ! Value of variable fCass

!         ! Variable of activation r (related with I_to current)
!         ! nauxi_exm = 7
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1586 (col 2)

!         vinf = 1.0_rp + exp((20.0_rp - elmlo_exm(tt)) / 6.0_rp)
!         vinf = 1.0_rp / vinf              ! Value of r_infinite

!         tau1 = 9.5_rp * exp(-((40.0_rp + elmlo_exm(tt)) * (40.0_rp + elmlo_exm(tt))) / 1800.0_rp)
!         taux = tau1 + 0.8_rp              ! Value of tau_r
!         xitaux = 1.0_rp / taux            ! Value of 1/tau_r

!         vaux0 = vaulo_exm(7,tt)      ! Value of r in previous time step 

!         ! Scheme for numerical integration
!         if (kfl_odets_exm == 1) then      ! Runge - Kutta 4

!            vaux1 = (vinf - vaux0) * xitaux;
!            vaux2 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux1)) * xitaux;
!            vaux3 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux2)) * xitaux;
!            vaux4 = (vinf - (vaux0 + dtimon * vaux3)) * xitaux;
!            vaux5 = vaux0 + (dtimon / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4); 
!            vaux3 = vaux5

!         else if (kfl_odets_exm == 2) then      ! Crank - Nicholson

!            vaux1 = dtimon * xitaux
!            vaux2 = 1.0_rp + thnew * vaux1
!            vaux2 = 1.0_rp / vaux2
!            vaux3 = (1.0_rp - thold * vaux1) * vaux0 + vaux1 * vinf
!            vaux3 = vaux3 * vaux2       

!         else if (kfl_odets_exm == 3) then      ! Rush - Larsen

!            vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/capacit * xitaux)    

!         end if

!         vaulo_exm(7,t) = vaux3      ! Value of variable r

!         ! Variable of activation s (related with I_to current)
!         ! nauxi_exm = 8
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1586 (col 2)

!         !           ituss_exm==  !(3=epicardial, 1=endocardial, 2=M cells)
!         if(ituss_exm == 3)  then

!            vinf = 1.0_rp + exp((20.0_rp + elmlo_exm(tt)) / 5.0_rp)
!            vinf = 1.0_rp / vinf              ! Value of s_infinite

!            tau1 = 85.0_rp * exp(-((45.0_rp + elmlo_exm(tt)) * (45.0_rp + elmlo_exm(tt))) / 320.0_rp)
!            tau2 = 1.0_rp + exp((-20.0_rp + elmlo_exm(tt)) / 5.0_rp)
!            tau2 = 5.0_rp / tau2
!            taux = tau1 + tau2 + 3.0_rp !)*1.6_rp       ! Value of tau_s
!            xitaux = 1.0_rp / taux            ! Value of 1/tau_s

!         else if (ituss_exm == 1) then

!            vinf = 1.0_rp + exp((28.0_rp + elmlo_exm(tt)) / 5.0_rp)
!            vinf = 1.0_rp / vinf              ! Value of s_infinite

!            tau1 = 1000.0_rp * exp(-((67.0_rp + elmlo_exm(tt)) * (67.0_rp + elmlo_exm(tt))) / 1000.0_rp)
!            taux = tau1 + 8.0_rp !)*1.6_rp       ! Value of tau_s
!            xitaux = 1.0_rp / taux            ! Value of 1/tau_s 

!         else if (ituss_exm == 2) then

!            vinf = 1.0_rp + exp((20.0_rp + elmlo_exm(tt)) / 5.0_rp)
!            vinf = 1.0_rp / vinf              ! Value of s_infinite

!            tau1 = 85.0_rp * exp(-((45.0_rp + elmlo_exm(tt)) * (45.0_rp + elmlo_exm(tt))) / 320.0_rp)
!            tau2 = 1.0_rp + exp((-20.0_rp + elmlo_exm(tt)) / 5.0_rp)
!            tau2 = 5.0_rp / tau2
!            taux = tau1 + tau2 + 3.0_rp !)*1.6_rp       ! Value of tau_s
!            xitaux = 1.0_rp / taux            ! Value of 1/tau_s 

!         end if

!         vaux0 = vaulo_exm(8,tt)      ! Value of s in previous time step 

!         ! Scheme for numerical integration
!         if (kfl_odets_exm == 1) then      ! Runge - Kutta 4

!            vaux1 = (vinf - vaux0) * xitaux;
!            vaux2 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux1)) * xitaux;
!            vaux3 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux2)) * xitaux;
!            vaux4 = (vinf - (vaux0 + dtimon * vaux3)) * xitaux;
!            vaux5 = vaux0 + (dtimon / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4); 
!            vaux3 = vaux5

!         else if (kfl_odets_exm == 2) then      ! Crank - Nicholson

!            vaux1 = dtimon * xitaux
!            vaux2 = 1.0_rp + thnew * vaux1
!            vaux2 = 1.0_rp / vaux2
!            vaux3 = (1.0_rp - thold * vaux1) * vaux0 + vaux1 * vinf
!            vaux3 = vaux3 * vaux2       

!         else if (kfl_odets_exm == 3) then      ! Rush - Larsen

!            vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/capacit * xitaux)    

!         end if

!         vaulo_exm(8,t) = vaux3      ! Value of variable s

!         ! Variable of activation xs (related with I_Ks current)
!         ! nauxi_exm = 9
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1586 (col 2)

!         vinf = 1.0_rp + exp((-5.0_rp - elmlo_exm(tt)) / 14.0_rp)
!         vinf = 1.0_rp / vinf              ! Value of xs_infinite

!         alphx = 1.0_rp + exp((5.0_rp - elmlo_exm(tt)) / 6.0_rp)
!         alphx = 1400.0_rp / sqrt(alphx)   ! Value of alpha_xs

!         betax = 1.0_rp + exp((-35.0_rp + elmlo_exm(tt)) / 15.0_rp)
!         betax = 1.0_rp / betax            ! Value of beta_xs

!         taux = (alphx * betax) + 80.0_rp             ! Value of tau_xs
!         xitaux = 1.0_rp / taux            ! Value of 1/tau_xs

!         vaux0 = vaulo_exm(9,tt)      ! Value of xs in previous time step 

!         ! Scheme for numerical integration
!         if (kfl_odets_exm == 1) then      ! Runge - Kutta 4

!            vaux1 = (vinf - vaux0) * xitaux;
!            vaux2 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux1)) * xitaux;
!            vaux3 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux2)) * xitaux;
!            vaux4 = (vinf - (vaux0 + dtimon * vaux3)) * xitaux;
!            vaux5 = vaux0 + (dtimon / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4); 
!            vaux3 = vaux5

!         else if (kfl_odets_exm == 2) then      ! Crank - Nicholson

!            vaux1 = dtimon * xitaux
!            vaux2 = 1.0_rp + thnew * vaux1
!            vaux2 = 1.0_rp / vaux2
!            vaux3 = (1.0_rp - thold * vaux1) * vaux0 + vaux1 * vinf
!            vaux3 = vaux3 * vaux2       

!         else if (kfl_odets_exm == 3) then      ! Rush - Larsen

!            vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/capacit * xitaux)    

!         end if

!         vaulo_exm(9,t) = vaux3      ! Value of variable xs

!         ! Variable of activation xr1 (related with I_Kr1 current)
!         ! nauxi_exm = 10
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1586 (col 2)

!         vinf = 1.0_rp + exp((-26.0_rp - elmlo_exm(tt)) / 7.0_rp)
!         vinf = 1.0_rp / vinf              ! Value of xr1_infinite

!         alphx = 1.0_rp + exp((-45.0_rp - elmlo_exm(tt)) / 10.0_rp)
!         alphx = 450.0_rp / alphx          ! Value of alpha_xr1

!         betax = 1.0_rp + exp((30.0_rp + elmlo_exm(tt)) / 11.5_rp)
!         betax = 6.0_rp / betax            ! Value of beta_xr1

!         taux = alphx * betax              ! Value of tau_xr1
!         xitaux = 1.0_rp / taux            ! Value of 1/tau_xr1

!         vaux0 = vaulo_exm(10,tt)     ! Value of xr1 in previous time step 

!         ! Scheme for numerical integration
!         if (kfl_odets_exm == 1) then      ! Runge - Kutta 4

!            vaux1 = (vinf - vaux0) * xitaux;
!            vaux2 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux1)) * xitaux;
!            vaux3 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux2)) * xitaux;
!            vaux4 = (vinf - (vaux0 + dtimon * vaux3)) * xitaux;
!            vaux5 = vaux0 + (dtimon / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4); 
!            vaux3 = vaux5

!         else if (kfl_odets_exm == 2) then      ! Crank - Nicholson

!            vaux1 = dtimon * xitaux
!            vaux2 = 1.0_rp + thnew * vaux1
!            vaux2 = 1.0_rp / vaux2
!            vaux3 = (1.0_rp - thold * vaux1) * vaux0 + vaux1 * vinf
!            vaux3 = vaux3 * vaux2       

!         else if (kfl_odets_exm == 3) then      ! Rush - Larsen

!            vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/capacit * xitaux)    

!         end if

!         vaulo_exm(10,t) = vaux3     ! Value of variable xr1

!         ! Variable of activation xr2 (related with I_Kr2 current)
!         ! nauxi_exm = 11
!         ! "A model for ventricular tissue", Ten Tusscher, Page H1586 (col 2)

!         vinf = 1.0_rp + exp((88.0_rp + elmlo_exm(tt)) / 24.0_rp)
!         vinf = 1.0_rp / vinf              ! Value of xr2_infinite

!         alphx = 1.0_rp + exp((-60.0_rp - elmlo_exm(tt)) / 20.0_rp)
!         alphx = 3.0_rp / alphx            ! Value of alpha_xr2

!         betax = 1.0_rp + exp((-60.0_rp + elmlo_exm(tt)) / 20.0_rp)
!         betax = 1.12_rp / betax           ! Value of beta_xr2

!         taux = alphx * betax              ! Value of tau_xr2
!         xitaux = 1.0_rp / taux            ! Value of 1/tau_xr2

!         vaux0 = vaulo_exm(11,tt)     ! Value of xr2 in previous time step 

!         ! Scheme for numerical integration
!         if (kfl_odets_exm == 1) then      ! Runge - Kutta 4

!            vaux1 = (vinf - vaux0) * xitaux;
!            vaux2 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux1)) * xitaux;
!            vaux3 = (vinf - (vaux0 + 0.5_rp * dtimon * vaux2)) * xitaux;
!            vaux4 = (vinf - (vaux0 + dtimon * vaux3)) * xitaux;
!            vaux5 = vaux0 + (dtimon / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4); 
!            vaux3 = vaux5

!         else if (kfl_odets_exm == 2) then      ! Crank - Nicholson

!            vaux1 = dtimon * xitaux
!            vaux2 = 1.0_rp + thnew * vaux1
!            vaux2 = 1.0_rp / vaux2
!            vaux3 = (1.0_rp - thold * vaux1) * vaux0 + vaux1 * vinf
!            vaux3 = vaux3 * vaux2       

!         else if (kfl_odets_exm == 3) then      ! Rush - Larsen

!            vaux3 = vinf - (vinf - vaux0) * exp(-dtimon/capacit * xitaux)    

!         end if

!         vaulo_exm(11,t) = vaux3     ! Value of variable xr2


!         elmlo_exm(tt)=elmlo_exm(t) 
!         itim = itim + 1
! !!        if (j>=10) then
! !!           write(996,*) vcolo_exm(1,2)
! !!           write(997,*) elmlo_exm(2)
! !!           j = 0
! !!        end if
!         j = j + 1
!         vcolo_exm(:,tt)=vcolo_exm(:,t)
!         vaulo_exm(:,tt)=vaulo_exm(:,t)
!         viclo_exm(:,tt)=viclo_exm(:,t)
!         t = 1
!         !write(993,*) vcolo_exm(:,2)  !JAS
!         !write(994,*) vaulo_exm(:,2)
!         !write(995,*) elmlo_exm(2)
!      end do
!      !call exm_analyz(atim)
!      batim = real(itim) * dtimon
!   end do
!   !write(993,*) vcolo_exm(:,2)  !JAS
!   !write(994,*) vaulo_exm(:,2)
!   !write(995,*) elmlo_exm(2)

!   if(ituss_exm==3) then !EPICARDIUM
!      !vauxi_exm_initial(:,1,1)=vaulo_exm(:,tt)
!      !vconc_initial(:,1,1) = vcolo_exm(:,tt)
!      !vminimate_exm(1,1) = elmlo_exm(tt)
!      vauxi_exm_initial(:,3,1)=vaulo_exm(:,t)
!      vconc_initial(:,3,1) = vcolo_exm(:,t)
!      vminimate_exm(3,1) = elmlo_exm(t)
!   else if (ituss_exm==1) then !ENDOCARDIUM
!      !vauxi_exm_initial(:,2,1)=vaulo_exm(:,tt)
!      !vconc_initial(:,2,1) = vcolo_exm(:,tt)
!      !vminimate_exm(1,1) = elmlo_exm(tt)
!      vauxi_exm_initial(:,1,1)=vaulo_exm(:,t)
!      vconc_initial(:,1,1) = vcolo_exm(:,t)
!      vminimate_exm(1,1) = elmlo_exm(t)
!   else if (ituss_exm==2) then !MIDMYOCARDIUM
!      !vauxi_exm_initial(:,3,1)=vaulo_exm(:,tt)
!      !vconc_initial(:,3,1) = vcolo_exm(:,tt)
!      !vminimate_exm(1,1) = elmlo_exm(tt)
!      vauxi_exm_initial(:,2,1)=vaulo_exm(:,t)
!      vconc_initial(:,2,1) = vcolo_exm(:,t)
!      vminimate_exm(2,1) = elmlo_exm(t)
!   end if


  !end if

!end subroutine exm_oceihe
