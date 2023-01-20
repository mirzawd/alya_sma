!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!!------------------------------------------------------------------------
!!> @addtogroup ExmediIonicurrents
!!> @{
!!> @file    exm_tt2006_ionicurrents.f90
!!> @author  Jazmin Aguado-Sierra
!!> @brief   Currents calculation for TenTuscher-Panfilov 2006 heterogeneous model
!!> @details Calculate every CURRENT at every time step\n
!!! The variables are 18:\n
!!!   - xina:   I_Na    -->  vicel_exm(1,ipoin,1)  -->  fast current of sodium \n
!!!   - xical:  I_CaL   -->  vicel_exm(2,ipoin,1)  -->  current of L-Type calcium \n
!!!   - xito:   I_to    -->  vicel_exm(3,ipoin,1)  -->  transient outward current \n
!!!   - xiks:   I_Ks    -->  vicel_exm(4,ipoin,1)  -->  slow delayed rectifier current \n
!!!   - xikr:   I_Kr    -->  vicel_exm(5,ipoin,1)  -->  rapid delayed rectifier current \n
!!!   - xik1:   I_K1    -->  vicel_exm(6,ipoin,1)  -->  inward rectifier K current \n
!!!   - xinaca: I_NaCa  -->  vicel_exm(7,ipoin,1)  -->  sodium/calcium exchanger current \n
!!!   - xinak:  I_NaK   -->  vicel_exm(8,ipoin,1)  -->  sodium/potassium pump current \n
!!!   - xipca:  I_pCa   -->  vicel_exm(9,ipoin,1)  -->  calcium plateau current \n
!!!   - xipk:   I_pK    -->  vicel_exm(10,ipoin,1) -->  potassium plateau current \n
!!!   - xibna:  I_bNa   -->  vicel_exm(11,ipoin,1) -->  sodium background current \n
!!!   - xibca:  I_bCa   -->  vicel_exm(12,ipoin,1) -->  calcium background current \n
!!!   - xileak: I_leak  -->  vicel_exm(13,ipoin,1) -->  leak from SR to the cytoplasm \n
!!!   - xiup:   I_up    -->  vicel_exm(14,ipoin,1) -->  pump taking up Ca in the SR \n
!!!   - xirel:  I_rel   -->  vicel_exm(15,ipoin,1) -->  Ca-induced Ca release (CICR) \n
!!!   - xistim: I_stim  -->  vicel_exm(16,ipoin,1) -->  external stimulus current \n
!!!   - xiax:   I_ax    -->  vicel_exm(17,ipoin,1) -->  UNUSED IN TTHET\n
!!!   - xixfer: I_xfer  -->  vicel_exm(18,ipoin,1) -->  diffusive Ca current between diadic subspace and bulk cytoplasm \n
!!! The order is important because it is related with the order of nicel_exm variables (exm_memall.f90):\n
!!> @} 
!!!-----------------------------------------------------------------------
!
!subroutine exm_tt2006_ionicurrents(ipoin,xioni,dioni)
!
!  use      def_parame
!  use      def_master
!  use      def_elmtyp
!  use      def_domain
!  use      def_exmedi
!  use      def_kermod,           only :  kfl_celltype_fun
!
!
!  ! DEFINITION OF VARIABLES
!  implicit none
!  integer(ip), intent(in) :: ipoin
!  real(rp), intent(out) :: xioni
!  real(rp), intent(out)   :: dioni
!  integer(ip) :: iicel, ituss_exm
!  real(rp)    :: xina, xical, xito, xiks, xikr, xik1, xinaca, xinak, xipca, xipk
!  real(rp)    :: xibna, xibca, xileak, xiup, xirel, xistim, xiax
!  real(rp)    :: enapot, ecapot, ekpot, ekspot 
!  real(rp)    :: groupcon, xgroupcon, rhsx
!  real(rp)   ::  nagmax, dosis1, dosis2, dosis3, i50cal, i50kr, i50na
!  real(rp)   ::  calgmax, gcaldrug, gkrdrug, gnadrug, gitodrug
!  real(rp)   ::  celltem, rgascon, faracon, xtrkcon, xtrnacon, xtrcacon
!  real(rp)   ::  togmaxmce, togmaxepi, togmaxend
!  real(rp)   ::  ksgmaxend, ksgmaxepi, ksgmaxmce, ksnaperm
!  real(rp)   ::  krgmax, k1, k2, k3, k4, dtimeEP
!  real(rp)   ::  k1gmax
!  real(rp)   ::  nacamax, gamanaca, kmcanaca, kmnanaca, satunaca, alfanaca
!  real(rp)   ::  nakmax, kmknak, kmnanak
!  real(rp)   ::  pcagmax, kpcapca
!  real(rp)   ::  pkgmax
!  real(rp)   ::  bnagmax
!  real(rp)   ::  bcagmax
!  real(rp)   ::  upvmax, upkcon
!  real(rp)   ::  leakvmax, aux1
!  real(rp)   ::  bufcyto, bufsr, bufss, kbufcyto, kbufsr, kbufss
!  integer(ip) :: kfl_odets_exm, clmat, n, imate
!  real(rp)    :: thold, thnew
!  real(rp)    :: vinf, alphx, alph1, alph2, betax, beta1, beta2, gammax, taux, tau1, tau2, tau3, xitaux
!  real(rp)    :: vaux0, vaux1, vaux2, vaux3, vaux4, vaux5
!  real(rp)    :: caitot, casrto, naitot, kitot, cassto, rprime
!  real(rp)    :: val0
!  real(rp)   ::  volcyto, volsr, capacit, volss
!  real(rp)   ::  faradco
!
!  xito = 0.0_rp
!  xiks = 0.0_rp
!
!print *,"================== SUCCESS =================="
!
!   if(kfl_celltype_fun == 0_ip) then !if there is no celltype field
!      call runend("exm_tt2006_ionicurrents: Cell type filed is required")
!   end if
!
!
!  if(INOTMASTER) then
!
!     clmat= 0
!     n=0
!     do imate= 1,nmate_exm
!        clmat=kfl_cellmod(imate)
!        if(clmat == 4) then
!           n=imate
!        end if
!        !Vector TT of Imate => 5-3 == 2
!     end do
!!!!!  DRUG DEFINITION TO IKR, INA OR ICAL
!     if(kfl_drugsmate_exm(n)==1_ip)then
!        dosis1 = drugdmate_exm(1,n)
!        i50cal = drugdmate_exm(2,n)
!        dosis2 = drugdmate_exm(3,n)
!        i50kr = drugdmate_exm(4,n)
!        dosis3 = drugdmate_exm(5,n)
!        i50na = drugdmate_exm(6,n)
!        gcaldrug = 1.0_rp / (1.0_rp + (dosis1/i50cal))
!        gkrdrug = 1.0_rp / (1.0_rp + (dosis2/i50kr))
!     else
!        gcaldrug = 1.0_rp
!        gkrdrug = 1.0_rp
!        gnadrug = 1.0_rp
!        gitodrug = 1.0_rp
!     endif
!
!     ! CONSTANTS DEFINITION  
!     rgascon = 8314.472_rp       ! R [mJ/(K*mol)] --> ideal gas constant
!     celltem = 310.0_rp          ! T     [K] --> cell system temperature 
!     faracon = 96485.3415_rp     ! F [C/mol] --> faraday constant
!     xtrkcon = 5.4_rp            ! K_o  [mM] --> extracellular concentration of potassium (K+)                !EPBench =5.4_rp [mM]
!     ! xtrkcon = 9.0             ! K_o  [mM] --> extracellular concentration of potassium (K+) for Ischemia    !EPBench
!     xtrnacon = 140.0_rp         ! Na_o [mM] --> extracellular concentration of sodium (Na+)
!     xtrcacon = 2.0_rp           ! Ca_o [mM] --> extracellular concentration of calcium (Ca++)
!
!     nagmax = 14.838_rp * gnadrug             ! G_Na [nS/pF]     --> maximal conductance of I_Na current
!
!     calgmax = 0.0000398_rp * gcaldrug        ! G_CaL [L/(mF*s)] --> maximal conductance of I_CaL current
!
!     togmaxmce = 0.294_rp * ttparmate_exm(1,1,n)    ! G_to_M   [nS/pF] --> maximal conductance of I_to current for M Cell
!     togmaxepi = 0.294_rp * ttparmate_exm(1,1,n)    ! G_to_epi [nS/pF] --> max conductance of I_to current for epicardium
!     togmaxend = 0.073_rp * ttparmate_exm(1,1,n)    ! G_to_end [nS/pF] --> max conductance of I_to current for endocardium
!
!     ksgmaxmce = 0.098_rp * ttparmate_exm(1,2,n)    ! G_Ks_M   [nS/pF] --> max conductance of I_Ks current for M Cell
!     ksgmaxepi = 0.392_rp * ttparmate_exm(1,2,n)    ! G_Ks_epi [nS/pF] --> max conductance of I_Ks current for epicardium
!     ksgmaxend = 0.392_rp * ttparmate_exm(1,2,n)    ! G_Ks_end [nS/pF] --> max conductance of I_Ks current for endocardium
!     ksnaperm = 0.03_rp                       ! p_KNa [adim] --> relative I_Ks current permeability to Na+
!
!     krgmax = 0.153_rp * gkrdrug              ! G_Kr [nS/pF] --> max conductance of I_Kr current
!
!     k1gmax = 5.405_rp * ttparmate_exm(1,3,n)       ! G_K1 [nS/pF] --> max conductance of I_K1 current
!
!     nacamax = 1000.0_rp     ! k_NaCa [pA/pF] --> max I_NaCa current
!     gamanaca = 0.35_rp      ! gama    [adim] --> voltage dependent parameter of I_NaCa current
!     kmcanaca = 1.38_rp      ! K_mCa     [mM] --> Cai half-saturation constant for I_NaCa current
!     kmnanaca = 87.5_rp      ! K_mNai    [mM] --> Nai half-saturation constant for I_NaCa current
!     satunaca = 0.1_rp       ! k_sat   [adim] --> saturation factor I_NaCa current
!     alfanaca = 2.5_rp       ! alfa    [adim] --> factor enhancing outward nature of I_NaCa current
!
!     nakmax = 2.724_rp       ! P_NaK [pA/pF] --> max I_NaK current
!     kmknak = 1.0_rp         ! K_mK  [mM] --> K_o half-saturation constant of I_NaK current
!     kmnanak = 40.0_rp       ! K_mNa [mM] --> Nai half-saturation constant of I_NaK current
!
!     pcagmax = 0.1238_rp     ! G_pCa [pA/pF] --> max conductance of I_pCa current
!     kpcapca = 0.0005_rp     ! K_pCa [mM] --> Cai half-saturation constant of I_pCa current
!
!     pkgmax = 0.0146_rp      ! G_pK  [nS/pF] --> max conductance of I_pK current
!
!     bnagmax = 0.00029_rp    ! G_bNa [nS/pF] --> max conductance of I_bNa current
!
!     bcagmax = 0.000592_rp   ! G_bCa [nS/pF] --> max conductance of I_bCa current
!
!     upvmax = 0.006375_rp    ! V_maxup [mM/ms] --> max I_up current
!     upkcon = 0.00025_rp     ! K_up [mM] --> half-saturation constant of I_up current
!
!     aux1 = ttparmate_exm(1,4,n)
!     leakvmax = 0.00036_rp * aux1     ! V_leak [1/ms] --> max I_leak current
!
!     bufcyto = 0.2_rp       ! Buf_c   [mM] --> total cytoplasmic buffer concentration
!     kbufcyto = 0.001_rp    ! K_Bufc  [mM] --> Cai half-saturation constant for cytoplasmic buffer
!     bufsr = 10.0_rp        ! Buf_sr  [mM] --> total sarcoplasmic buffer concentration
!     kbufsr = 0.3_rp        ! K_Bufsr [mM] --> Ca_SR half-saturation constant for sarcoplasmic buffer
!     bufss = 0.4_rp         ! Buf_sr  [mM] --> total sarcoplasmic buffer concentration
!     kbufss = 0.00025_rp    ! K_Bufsr [mM] --> Ca_SR half-saturation constant for sarcoplasmic buffer
!
!     groupcon = rgascon * celltem / faracon          ! Group of constants: R * T / F  [mV] 
!     xgroupcon = 1.0_rp / groupcon
!
!     volcyto = 0.016404_rp        ! V_c [nL] --> cytoplasmic volume
!     volsr = 0.001094_rp  * ttparmate_exm(1,5,n)       ! V_SR [nL] --> sarcoplasmic reticulum (SR) volume
!     capacit = 0.185_rp           ! C [microF] --> Capacitance used in Ca_i_tot, Na_i_toy and K_i_tot
!     volss = 0.00005468_rp        ! V_SS
!
!     bufcyto = 0.2_rp             ! Buf_c   [mM] --> total cytoplasmic buffer concentration
!     kbufcyto = 0.001_rp          ! K_Bufc  [mM] --> Cai half-saturation constant for cytoplasmic buffer
!     bufsr = 10.0_rp              ! Buf_sr  [mM] --> total sarcoplasmic buffer concentration
!     kbufsr = 0.3_rp              ! K_Bufsr [mM] --> Ca_SR half-saturation constant for sarcoplasmic buffer
!     bufss = 0.4_rp               ! Buf_sr  [mM] --> total sarcoplasmic buffer concentration
!     kbufss = 0.00025_rp          ! K_Bufsr [mM] --> Ca_SR half-saturation constant for sarcoplasmic buffer
!
!     faradco = 96485.3415_rp      ! F [C/mol] --> faraday constant
!
!!!! First calculate the auxiliaries, Currents, then concentrations and then calculate new voltage,      
!
!     dtimeEP=dtime * 1000.0_rp
!
!     if (kfl_cellmod(nodemat(ipoin)) == 4) then
!
!        ituss_exm = int(celty_exm(ipoin), kind=ip)
!
!        thnew = 0.5_rp            
!        thold = 1.0_rp - thnew       
!        kfl_odets_exm = 1
!        
!        ! Concentration of total calcium in the Cytoplasm (Ca_i_tot)
!        ! nconc_exm = 1
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 2)
!
!        !   Ca_i_bufc = 1/(1+( Buf_c*K_buf_c)/(pow((Ca_i+K_buf_c), 2)))
!        vaux1 = vconc(1,ipoin,2) + kbufcyto
!        vaux1 = 1.0_rp / (vaux1 * vaux1)
!        vaux1 = vaux1 * kbufcyto * bufcyto    !Variable caibuf
!        vaux1 = 1.0_rp / (1.0_rp + vaux1)
!        vconc(6,ipoin,1) = vaux1
!        vconc(6,ipoin,2) = vaux1
!        ! d_caitot/d_t = - (I_CaL + I_bCa + I_pCa - 2 * I_NaCa) * C / (2 * V_c * F) + I_leak - I_up + I_rel
!        ! dCa_i_dt = Ca_i_bufc*((( (i_leak - i_up)*V_sr)/V_c+i_xfer) - (1*((i_b_Ca+i_p_Ca) - 2*i_NaCa)*Cm)/( 2*1*V_c*F))
!
!        vaux1 = - 2.0_rp * volcyto * faradco
!        vaux1 = 1.0_rp / vaux1  
!
!        vaux3 = (vicel_exm(12,ipoin,1) + vicel_exm(9,ipoin,1)) - (2.0_rp *  vicel_exm(7,ipoin,1) )
!        vaux3 = vaux1 * vaux3 * capacit
!        vaux2 = (vicel_exm(13,ipoin,1) - vicel_exm(14,ipoin,1)) * (volsr / volcyto) + vicel_exm(18,ipoin,1)
!        rhsx =  vconc(6,ipoin,2) * (vaux2 + vaux3) !
!
!        val0 = vconc(1,ipoin,2)
!
!        k1 = rhsx 
!        k2 = rhsx + 0.5_rp * dtimeEP * k1
!        k3 = rhsx + 0.5_rp * dtimeEP * k2 
!        k4 = rhsx + dtimeEP * k3 
!        caitot = val0 + ((dtimeEP/6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4))
!
!        !caitot = val0 + dtimeEP/xmccmmate_exm(nodemat(ipoin)) * rhsx
!        !caitot =  val0 + (0.5_rp * dtimeEP * (k1 + k2))
!        vconc(1,ipoin,1) = caitot       ! Value of Ca_i_tot concentration
!
!        ! Concentration of total calcium in the Sarcoplasmic Reticulum (Ca_sr_tot)
!        ! nconc_exm = 2
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 2)
!
!        vaux1 = vconc(2,ipoin,2)  + kbufsr
!        vaux1 = 1.0_rp / (vaux1 * vaux1)
!        vaux1 = vaux1 * kbufsr * bufsr
!        vaux1 = 1.0_rp / (1.0_rp + vaux1)
!        vconc(7,ipoin,1) = vaux1
!        vconc(7,ipoin,2) = vaux1
!
!        ! d_casrto/d_t = Ca_sr_bufsr*(i_up - (i_rel+i_leak))
!
!        vaux1 = vicel_exm(14,ipoin,1) - vicel_exm(13,ipoin,1) - vicel_exm(15,ipoin,1)
!
!        rhsx = vaux1 * vconc(7,ipoin,2)  !capacit *
!
!        val0 = vconc(2,ipoin,2)
!        k1 = rhsx 
!        k2 = rhsx + 0.5_rp * dtimeEP * k1
!        k3 = rhsx + 0.5_rp * dtimeEP * k2 
!        k4 = rhsx + dtimeEP * k3 
!        casrto = val0 + ((dtimeEP / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4))
!
!        !casrto = val0 + dtimeEP/xmccmmate_exm(nodemat(ipoin))  * rhsx
!        vconc(2,ipoin,1) = casrto       ! Value of Ca_sr_tot concentration
!
!        ! Concentration of total calcium in the diadic space (Ca_ss_tot)
!        ! nconc_exm = 5
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 2)
!        vaux1 = vconc(5,ipoin,2)  + kbufss
!        vaux1 = 1.0_rp / (vaux1 * vaux1)
!        vaux1 = vaux1 * kbufss * bufss
!        vaux1 = 1.0_rp / (1.0_rp + vaux1)
!        vconc(8,ipoin,1) = vaux1                    !Variable cassbufss
!        vconc(8,ipoin,2) = vaux1
!
!        ! d_cass/d_t = Ca_ss_bufss*(((-1*i_CaL)/( 2*1*V_ss*F)+( i_rel*V_sr)/V_ss) - ( i_xfer*V_c)/V_ss)
!
!        vaux1 = 1.0_rp / (2.0_rp * faradco * volss)
!        vaux1 = vaux1 * vicel_exm(2,ipoin,1) * capacit
!        vaux2 = vicel_exm(15,ipoin,1) * volsr / volss
!        vaux3 = vicel_exm(18,ipoin,1) * volcyto / volss
!
!        rhsx =  vconc(8,ipoin,2) * (- vaux1 + vaux2 - vaux3)
!
!        val0 = vconc(5,ipoin,2)
!
!        k1 = rhsx 
!        k2 = rhsx + 0.5_rp * dtimeEP * k1
!        k3 = rhsx + 0.5_rp * dtimeEP * k2 
!        k4 = rhsx + dtimeEP * k3 
!        cassto = val0 + ((dtimeEP / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4))
!        !cassto = val0 + dtimeEP/xmccmmate_exm(nodemat(ipoin))  * rhsx
!        vconc(5,ipoin,1) = cassto       ! Value of Ca_ss_tot concentration
!
!        ! Intracellular sodium concentration (Na_i_tot)
!        ! nconc_exm = 3
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 2)
!
!        ! d_naitot/d_t =  ((- 1*(i_Na+i_b_Na+ 3*i_NaK+ 3*i_NaCa))/( 1*V_c*F))*Cm
!
!        vaux1 = volcyto * faradco
!        vaux1 = 1.0_rp / vaux1 !* capacit
!
!        vaux3 = (3.0_rp * vicel_exm(8,ipoin,1)) + (3.0_rp * vicel_exm(7,ipoin,1)) 
!        vaux3 = vaux3 + vicel_exm(1,ipoin,1) + vicel_exm(11,ipoin,1)
!        rhsx = - capacit * vaux1 * vaux3
!
!        val0 = vconc(3,ipoin,2)
!
!        k1 = rhsx 
!        k2 = rhsx + 0.5_rp * dtimeEP * k1
!        k3 = rhsx + 0.5_rp * dtimeEP * k2 
!        k4 = rhsx + dtimeEP * k3 
!        naitot = val0 + ((dtimeEP / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4))
!        !naitot = val0 + dtimeEP/xmccmmate_exm(nodemat(ipoin))  * rhsx
!        vconc(3,ipoin,1) = naitot       ! Value of Na_i_tot concentration
!
!
!        ! Intracellular potassium concentration (K_i_tot)
!        ! nconc_exm = 4
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 2)
!
!        ! d_kitot/d_t = ((- 1*((i_K1+i_to+i_Kr+i_Ks+i_p_K+i_Stim) -  2*i_NaK))/( 1*V_c*F))*Cm
!
!        vaux1 = - volcyto * faradco
!        vaux1 = capacit / vaux1 !* capacit  
!
!        vaux3 = (vicel_exm(6,ipoin,1) + vicel_exm(3,ipoin,1) + vicel_exm(5,ipoin,1) &
!             + vicel_exm(4,ipoin,1) + vicel_exm(10,ipoin,1) + vicel_exm(16,ipoin,1)) - (2.0_rp * vicel_exm(8,ipoin,1)) 
!
!        rhsx = vaux1 * vaux3
!
!        val0 = vconc(4,ipoin,2)
!
!        k1 = rhsx 
!        k2 = rhsx + 0.5_rp * dtimeEP * k1
!        k3 = rhsx + 0.5_rp * dtimeEP * k2 
!        k4 = rhsx + dtimeEP * k3 
!        kitot = val0 + ((dtimeEP / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4))
!        !kitot = val0 + dtimeEP/xmccmmate_exm(nodemat(ipoin))  * rhsx
!        vconc(4,ipoin,1) = kitot          ! Value of K_i_tot concentration
!
!        ! Variable of activation O (related with I_rel current)
!        ! nconc = 10
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1)
!        !  kcasr = max_sr - (max_sr - min_sr)/(1+(pow((EC/Ca_SR), 2)))
!        vaux1 = 1.0_rp + ((1.5_rp / vconc(2,ipoin,2))*(1.5_rp / vconc(2,ipoin,2)))  !conc2 = CaSRfree = CaSR
!        vaux1 = 1.0_rp / vaux1
!        vaux1 = 2.5_rp - (1.5_rp * vaux1)  !kCaSR
!        !     O = ( k1*pow(Ca_ss, 2)*R_prime)/(k3+ k1*pow(Ca_ss, 2))      
!        vaux3 = 0.15_rp / vaux1   !K1
!        vaux2 = 0.045_rp * vaux1  !K2
!
!
!        rhsx = 1.0_rp / (0.06_rp + (vaux3*vconc(5,ipoin,2)*vconc(5,ipoin,2)))      !O for calcium dynamics
!        rhsx = vaux3 * vconc(5,ipoin,2) * vconc(5,ipoin,2) * vconc(9,ipoin,2) * rhsx !O
!        vconc(10,ipoin,1) = rhsx 
!        vconc(10,ipoin,2) = rhsx
!
!        ! Variable of activation O (related with I_rel current)
!        ! nconc= 9
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1)
!        !  dR_prime_dt = - k2*Ca_ss*R_prime+ k4*(1 - R_prime)
!
!        vaux3 = - vaux2 * vconc(5,ipoin,2) * vconc(9,ipoin,2) 
!        rhsx = vaux3 + (0.005_rp * (1.0_rp - vconc(9,ipoin,2)))
!
!        val0 = vconc(9,ipoin,2)    ! Value of Rprime in previous dtimeEP step 
!
!        k1 = rhsx 
!        k2 = rhsx + 0.5_rp * dtimeEP * k1
!        k3 = rhsx + 0.5_rp * dtimeEP * k2 
!        k4 = rhsx + dtimeEP * k3 
!        rprime = val0 + ((dtimeEP / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4))
!        !rprime = val0 + dtimeEP/xmccmmate_exm(nodemat(ipoin))  * rhsx
!        vconc(9,ipoin,1) = rprime          
!
!        ! Variable of activation m (related with I_Na current)
!        ! nauxi_exm = 1
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1585 (col 2)
!
!        vinf = 1.0_rp + exp((-56.86_rp - elmag(ipoin,ITER_K)) / 9.03_rp)
!        vinf = vinf * vinf
!        vinf = 1.0_rp / vinf              ! Value of m_infinite
!
!        alphx = 1.0_rp + exp((-60.00_rp - elmag(ipoin,ITER_K)) / 5.0_rp)
!        alphx = 1.0_rp / alphx            ! Value of alpha_m
!        !   beta_m = 0.1/(1+(exp(((V+35)/5))))+0.1/(1+(exp(((V - 50)/200))))
!        beta1 = 1.0_rp + exp(( 35.00_rp + elmag(ipoin,ITER_K)) / 5.0_rp)
!        beta1 = 0.1_rp / beta1
!        beta2 = 1.0_rp + exp((-50.00_rp + elmag(ipoin,ITER_K)) / 200.0_rp)
!        beta2 = 0.1_rp / beta2
!        betax = beta1 + beta2             ! Value of beta_m
!
!        taux  = betax * alphx             ! Value of tau_m
!        xitaux = 1.0_rp / taux            ! Value of 1/tau_m 
!
!        vaux0 = vauxi_exm(1,ipoin,2)      ! Value of m in previous time step 
!
!
!        ! Scheme for numerical integration
!        if (kfl_odets_exm == 1) then      ! Runge - Kutta 4
!
!           vaux1 = (vinf - vaux0) * xitaux
!           vaux2 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux1)) * xitaux
!           vaux3 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux2)) * xitaux
!           vaux4 = (vinf - (vaux0 + dtimeEP * vaux3)) * xitaux
!           vaux5 = vaux0 + (dtimeEP / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4) 
!           vaux3 = vaux5
!
!        else if (kfl_odets_exm == 2) then      ! Crank - Nicholson
!
!           vaux1 = dtimeEP * xitaux
!           vaux2 = 1.0_rp + thnew * vaux1
!           vaux2 = 1.0_rp / vaux2
!           vaux3 = (1.0_rp - thold * vaux1) * vaux0 + vaux1 * vinf
!           vaux3 = vaux3 * vaux2       
!
!        else if (kfl_odets_exm == 3) then      ! Rush - Larsen
!
!           vaux3 = vinf - (vinf - vaux0) * exp(-dtimeEP/xmccmmate_exm(nodemat(ipoin)) * xitaux)    
!
!        end if
!
!        vauxi_exm(1,ipoin,1) = vaux3      ! Value of variable m
!
!
!        ! Variable of activation h (related with I_Na current)
!        ! nauxi_exm = 2
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1586 (col 1)
!
!        vinf = 1.0_rp + exp((71.55_rp + elmag(ipoin,ITER_K)) / 7.43_rp)
!        vinf = vinf * vinf
!        vinf = 1.0_rp / vinf              ! Value of h_infinite
!
!        if (elmag(ipoin,ITER_K) >= -40.0_rp) then
!
!           alphx = 0.0_rp                  ! Value of alpha_h
!
!           betax = 0.13_rp * ( 1.0_rp + exp((-10.66_rp - elmag(ipoin,ITER_K)) / 11.1_rp) )
!           betax = 0.77_rp / betax         ! Value of beta_h
!
!        else
!
!           alphx = 0.057_rp * exp((-80.0_rp - elmag(ipoin,ITER_K)) / 6.8_rp)        ! Value of alpha_h
!           !   ( 2.7*(exp(( 0.079*V)))+ 310000*(exp((0.3485*V)))) 
!           beta1 = 2.7_rp * exp(0.079_rp * elmag(ipoin,ITER_K))
!           beta2 = 310000.0_rp * exp(0.3485_rp * elmag(ipoin,ITER_K))
!           betax = beta1 + beta2           ! Value of beta_h
!
!        end if
!
!        xitaux  = betax + alphx          ! Value of 1/tau_h
!
!        vaux0 = vauxi_exm(2,ipoin,2)     ! Value of h in previous time step 
!
!
!        ! Scheme for numerical integration
!        if (kfl_odets_exm == 1) then     ! Runge - Kutta 4
!
!           vaux1 = (vinf - vaux0) * xitaux
!           vaux2 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux1)) * xitaux
!           vaux3 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux2)) * xitaux
!           vaux4 = (vinf - (vaux0 + dtimeEP * vaux3)) * xitaux
!           vaux5 = vaux0 + (dtimeEP / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4) 
!           vaux3 = vaux5
!
!        else if (kfl_odets_exm == 2) then      ! Crank - Nicholson
!
!           vaux1 = dtimeEP * xitaux
!           vaux2 = 1.0_rp + thnew * vaux1
!           vaux2 = 1.0_rp / vaux2
!           vaux3 = (1.0_rp - thold * vaux1) * vaux0 + vaux1 * vinf
!           vaux3 = vaux3 * vaux2       
!
!        else if (kfl_odets_exm == 3) then      ! Rush - Larsen
!
!           vaux3 = vinf - (vinf - vaux0) * exp(-dtimeEP/xmccmmate_exm(nodemat(ipoin))* xitaux)    
!
!        end if
!        vauxi_exm(2,ipoin,1) = vaux3      ! Value of variable h
!
!        ! Variable of activation j (related with I_Na current)
!        ! nauxi_exm = 3
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1586 (col 1)
!
!        vinf = 1.0_rp + exp((71.55_rp + elmag(ipoin,ITER_K)) / 7.43_rp)
!        vinf = vinf * vinf
!        vinf = 1.0_rp / vinf              ! Value of j_infinite
!
!        if (elmag(ipoin,ITER_K) >= -40.0_rp ) then
!
!           alphx = 0                       ! Value of alpha_j
!
!           beta1 = 1.0_rp + exp(-0.1_rp * (32.0_rp + elmag(ipoin,ITER_K)))
!           beta1 = 1.0_rp / beta1
!           beta2 = 0.6_rp * exp(0.057_rp * elmag(ipoin,ITER_K))
!           betax = beta1 * beta2           ! Value of beta_j
!
!        else
!           !    (( ( (-25428*exp( 0.2444*V)) -  (6.948e-6*exp(-0.04391*V)) )*(V+37.78) ) / (1+exp( 0.311*(V+79.23)) ))
!           alph1 = 1.0_rp + exp(0.311_rp * (79.23_rp + elmag(ipoin,ITER_K)))
!           alph1 = 1.0_rp / alph1
!           alph2 = - 25428.0_rp * exp(0.2444_rp * elmag(ipoin,ITER_K))
!           alph2 = alph2 - (0.000006948_rp * exp(-0.04391_rp * elmag(ipoin,ITER_K)))
!           alph2 = alph2 * (elmag(ipoin,ITER_K) + 37.78_rp)
!           alphx = alph1 * alph2           ! Value of alpha_j
!
!           beta1 = 1.0_rp + exp(-0.1378_rp * (40.14_rp + elmag(ipoin,ITER_K)))
!           beta1 = 1.0_rp / beta1
!           beta2 = 0.02424_rp * exp(-0.01052_rp * elmag(ipoin,ITER_K))
!           betax = beta1 * beta2           ! Value of beta_j
!
!        end if
!
!        xitaux  = betax + alphx           ! Value of 1/tau_j
!
!        vaux0 = vauxi_exm(3,ipoin,2)      ! Value of j in previous time step 
!
!
!        ! Scheme for numerical integration
!        if (kfl_odets_exm == 1) then      ! Runge - Kutta 4
!
!           vaux1 = (vinf - vaux0) * xitaux
!           vaux2 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux1)) * xitaux
!           vaux3 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux2)) * xitaux
!           vaux4 = (vinf - (vaux0 + dtimeEP * vaux3)) * xitaux
!           vaux5 = vaux0 + (dtimeEP / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4) 
!           vaux3 = vaux5
!
!        else if (kfl_odets_exm == 2) then      ! Crank - Nicholson
!
!           vaux1 = dtimeEP * xitaux
!           vaux2 = 1.0_rp + thnew * vaux1
!           vaux2 = 1.0_rp / vaux2
!           vaux3 = (1.0_rp - thold * vaux1) * vaux0 + vaux1 * vinf
!           vaux3 = vaux3 * vaux2       
!
!        else if (kfl_odets_exm == 3) then      ! Rush - Larsen
!
!           vaux3 = vinf - (vinf - vaux0) * exp(-dtimeEP/xmccmmate_exm(nodemat(ipoin))     * xitaux)    
!
!        end if
!
!
!        vauxi_exm(3,ipoin,1) = vaux3      ! Value of variable j
!
!
!        ! Variable of activation d (related with I_CaL current)
!        ! nauxi_exm = 4
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1586 (col 1)
!
!        vinf = 1.0_rp + exp((-8.0_rp - elmag(ipoin,ITER_K)) / 7.5_rp)
!        vinf = 1.0_rp / vinf              ! Value of d_infinite
!
!        alphx = 1.0_rp + exp((-35.0_rp - elmag(ipoin,ITER_K)) / 13.0_rp)
!        alphx = (1.4_rp / alphx) + 0.25_rp  ! Value of alpha_d
!
!        betax = 1.0_rp + exp((5.0_rp + elmag(ipoin,ITER_K)) / 5.0_rp)
!        betax = 1.4_rp / betax            ! Value of beta_d
!
!        gammax = 1.0_rp + exp((50.0_rp - elmag(ipoin,ITER_K)) / 20.0_rp)
!        gammax = 1.0_rp / gammax          ! Value of gamma_d
!
!        taux = (alphx * betax) + gammax   ! Value of tau_d
!        xitaux = 1.0_rp / taux            ! Value of 1/tau_d
!
!        vaux0 = vauxi_exm(4,ipoin,2)      ! Value of d in previous time step 
!
!
!        ! Scheme for numerical integration
!        if (kfl_odets_exm == 1) then      ! Runge - Kutta 4
!
!           vaux1 = (vinf - vaux0) * xitaux
!           vaux2 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux1)) * xitaux
!           vaux3 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux2)) * xitaux
!           vaux4 = (vinf - (vaux0 + dtimeEP * vaux3)) * xitaux
!           vaux5 = vaux0 + (dtimeEP / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4) 
!           vaux3 = vaux5
!
!        else if (kfl_odets_exm == 2) then      ! Crank - Nicholson
!
!           vaux1 = dtimeEP * xitaux
!           vaux2 = 1.0_rp + thnew * vaux1
!           vaux2 = 1.0_rp / vaux2
!           vaux3 = (1.0_rp - thold * vaux1) * vaux0 + vaux1 * vinf
!           vaux3 = vaux3 * vaux2       
!
!        else if (kfl_odets_exm == 3) then      ! Rush - Larsen
!
!           vaux3 = vinf - (vinf - vaux0) * exp(-dtimeEP/xmccmmate_exm(nodemat(ipoin))     * xitaux)    
!
!        end if
!
!        vauxi_exm(4,ipoin,1) = vaux3      ! Value of variable d
!
!        ! Variable of activation f (related with I_CaL current)
!        ! nauxi_exm = 5
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1586 (col 1)
!
!        vinf = 1.0_rp + exp((20.0_rp + elmag(ipoin,ITER_K)) / 7.0_rp)
!        vinf = 1.0_rp / vinf              ! Value of f_infinite
!        ! tau_f = 1102.5*(exp(( - (pow((V+27),2))/225)))+200/(1+(exp(((13 - V)/10))))+180/(1+(exp(((V+30)/10))))+20
!        tau1 = 1102.5_rp * exp(-((27.0_rp + elmag(ipoin,ITER_K)) * (27.0_rp + elmag(ipoin,ITER_K))) / 225.0_rp)
!        tau2 = 1.0_rp + exp((13.0_rp - elmag(ipoin,ITER_K)) / 10.0_rp)
!        tau2 = 200.0_rp / tau2
!        tau3 = 1.0_rp + exp((30.0_rp + elmag(ipoin,ITER_K)) / 10.0_rp)
!        tau3 = 180.0_rp / tau3
!        taux = tau1 + tau2 + tau3 + 20.0_rp      ! Value of tau_f
!        xitaux = 1.0_rp / taux            ! Value of 1/tau_f
!
!        vaux0 = vauxi_exm(5,ipoin,2)      ! Value of f in previous time step 
!
!        ! Scheme for numerical integration
!        if (kfl_odets_exm == 1) then      ! Runge - Kutta 4
!
!           vaux1 = (vinf - vaux0) * xitaux
!           vaux2 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux1)) * xitaux
!           vaux3 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux2)) * xitaux
!           vaux4 = (vinf - (vaux0 + dtimeEP * vaux3)) * xitaux
!           vaux5 = vaux0 + (dtimeEP / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4) 
!           vaux3 = vaux5
!
!        else if (kfl_odets_exm == 2) then      ! Crank - Nicholson
!
!           vaux1 = dtimeEP * xitaux
!           vaux2 = 1.0_rp + thnew * vaux1
!           vaux2 = 1.0_rp / vaux2
!           vaux3 = (1.0_rp - thold * vaux1) * vaux0 + vaux1 * vinf
!           vaux3 = vaux3 * vaux2       
!
!        else if (kfl_odets_exm == 3) then      ! Rush - Larsen
!
!           vaux3 = vinf - (vinf - vaux0) * exp(-dtimeEP/xmccmmate_exm(nodemat(ipoin))     * xitaux)    
!
!        end if
!
!        vauxi_exm(5,ipoin,1) = vaux3      ! Value of variable f
!
!        ! Variable of activation f2 (related with I_CaL current)
!        ! nauxi_exm = 12
!        ! "A model for ventricular tissue", Ten Tusscher, 2006 paper
!        ! modified by JAS
!
!        vinf = 1.0_rp + exp((35.0_rp + elmag(ipoin,ITER_K)) / 7.0_rp)
!        vinf = (0.67_rp / vinf ) + 0.33_rp            ! Value of f2_infinite
!
!        tau1 = 600.0_rp * exp((-(25.0_rp + elmag(ipoin,ITER_K)) * (25.0_rp + elmag(ipoin,ITER_K))) / 170.0_rp)
!        tau2 = 1.0_rp + exp((25.0_rp - elmag(ipoin,ITER_K)) / 10.0_rp)
!        tau2 = 31.0_rp / tau2
!        tau3 = 1.0_rp + exp((30.0_rp + elmag(ipoin,ITER_K)) / 10.0_rp)
!        tau3 = 16.0_rp / tau3
!        taux = tau1 + tau2 + tau3      ! Value of tau_f2
!        xitaux = 1.0_rp / taux         ! Value of 1/tau_f2
!
!        vaux0 = vauxi_exm(12,ipoin,2)  ! Value of f2 in previous time step 
!
!
!        ! Scheme for numerical integration
!        if (kfl_odets_exm == 1) then      ! Runge - Kutta 4
!
!           vaux1 = (vinf - vaux0) * xitaux
!           vaux2 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux1)) * xitaux
!           vaux3 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux2)) * xitaux
!           vaux4 = (vinf - (vaux0 + dtimeEP * vaux3)) * xitaux
!           vaux5 = vaux0 + (dtimeEP / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4) 
!           vaux3 = vaux5
!
!        else if (kfl_odets_exm == 2) then      ! Crank - Nicholson
!
!           vaux1 = dtimeEP * xitaux
!           vaux2 = 1.0_rp + thnew * vaux1
!           vaux2 = 1.0_rp / vaux2
!           vaux3 = (1.0_rp - thold * vaux1) * vaux0 + vaux1 * vinf
!           vaux3 = vaux3 * vaux2       
!
!        else if (kfl_odets_exm == 3) then      ! Rush - Larsen
!
!           vaux3 = vinf - (vinf - vaux0) * exp(-dtimeEP/xmccmmate_exm(nodemat(ipoin))     * xitaux)    
!
!        end if
!
!        vauxi_exm(12,ipoin,1) = vaux3      ! Value of variable f2
!
!        ! Variable of activation fCass (related to I_CaL current)
!        ! nauxi_exm = 6
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1586 (col 2)
!
!        vinf = 1.0_rp + ((vconc(5,ipoin,2) / 0.05_rp)*(vconc(5,ipoin,2) / 0.05_rp)) !changed from vconc(5,ipoin,2)
!        vinf = (0.6_rp / vinf ) + 0.4_rp            ! Value of fCass_infinite
!
!        tau1 = 1.0_rp + ((vconc(5,ipoin,2) / 0.05_rp)*(vconc(5,ipoin,2) / 0.05_rp)) !changed from vconc(5,ipoin,2)
!        tau1 = (80.0_rp / tau1) + 2.0_rp           ! Value of tau_fCass    
!        xitaux = 1.0_rp / tau1            ! Value of 1/tau_f
!
!        vaux0 = vauxi_exm(6,ipoin,2)      ! Value of f in previous time step 
!
!        ! Scheme for numerical integration
!        if (kfl_odets_exm == 1) then      ! Runge - Kutta 4
!
!           vaux1 = (vinf - vaux0) * xitaux
!           vaux2 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux1)) * xitaux
!           vaux3 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux2)) * xitaux
!           vaux4 = (vinf - (vaux0 + dtimeEP * vaux3)) * xitaux
!           vaux5 = vaux0 + (dtimeEP / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4) 
!           vaux3 = vaux5
!
!        else if (kfl_odets_exm == 2) then      ! Crank - Nicholson
!
!           vaux1 = dtimeEP * xitaux
!           vaux2 = 1.0_rp + thnew * vaux1
!           vaux2 = 1.0_rp / vaux2
!           vaux3 = (1.0_rp - thold * vaux1) * vaux0 + vaux1 * vinf
!           vaux3 = vaux3 * vaux2       
!
!        else if (kfl_odets_exm == 3) then      ! Rush - Larsen
!
!           vaux3 = vinf - (vinf - vaux0) * exp(-dtimeEP/xmccmmate_exm(nodemat(ipoin))     * xitaux)    
!
!        end if
!
!        vauxi_exm(6,ipoin,1) = vaux3      ! Value of variable fCass
!
!        ! Variable of activation r (related with I_to current)
!        ! nauxi_exm = 7
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1586 (col 2)
!
!        vinf = 1.0_rp + exp((20.0_rp - elmag(ipoin,ITER_K)) / 6.0_rp)
!        vinf = 1.0_rp / vinf              ! Value of r_infinite
!
!        tau1 = 9.5_rp * exp(-((40.0_rp + elmag(ipoin,ITER_K)) * (40.0_rp + elmag(ipoin,ITER_K))) / 1800.0_rp)
!        taux = tau1 + 0.8_rp              ! Value of tau_r
!        xitaux = 1.0_rp / taux            ! Value of 1/tau_r
!
!        vaux0 = vauxi_exm(7,ipoin,2)      ! Value of r in previous time step 
!
!        ! Scheme for numerical integration
!        if (kfl_odets_exm == 1) then      ! Runge - Kutta 4
!
!           vaux1 = (vinf - vaux0) * xitaux
!           vaux2 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux1)) * xitaux
!           vaux3 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux2)) * xitaux
!           vaux4 = (vinf - (vaux0 + dtimeEP * vaux3)) * xitaux
!           vaux5 = vaux0 + (dtimeEP / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4) 
!           vaux3 = vaux5
!
!        else if (kfl_odets_exm == 2) then      ! Crank - Nicholson
!
!           vaux1 = dtimeEP * xitaux
!           vaux2 = 1.0_rp + thnew * vaux1
!           vaux2 = 1.0_rp / vaux2
!           vaux3 = (1.0_rp - thold * vaux1) * vaux0 + vaux1 * vinf
!           vaux3 = vaux3 * vaux2       
!
!        else if (kfl_odets_exm == 3) then      ! Rush - Larsen
!
!           vaux3 = vinf - (vinf - vaux0) * exp(-dtimeEP/xmccmmate_exm(nodemat(ipoin))     * xitaux)    
!
!        end if
!
!        vauxi_exm(7,ipoin,1) = vaux3      ! Value of variable r
!
!        ! Variable of activation s (related with I_to current)
!        ! nauxi_exm = 8
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1586 (col 2)
!
!        !           ituss_exm==  !(3=epicardial, 1=endocardial, 2=M cells)
!        if(ituss_exm == EXM_CELLTYPE_EPI)  then
!
!           vinf = 1.0_rp + exp((20.0_rp + elmag(ipoin,ITER_K)) / 5.0_rp)
!           vinf = 1.0_rp / vinf              ! Value of s_infinite
!
!           tau1 = 85.0_rp * exp(-((45.0_rp + elmag(ipoin,ITER_K)) * (45.0_rp + elmag(ipoin,ITER_K))) / 320.0_rp)
!           tau2 = 1.0_rp + exp((-20.0_rp + elmag(ipoin,ITER_K)) / 5.0_rp)
!           tau2 = 5.0_rp / tau2
!           taux = (tau1 + tau2 + 3.0_rp)*1.6_rp       ! Value of tau_s
!           xitaux = 1.0_rp / taux            ! Value of 1/tau_s
!
!        else if (ituss_exm == EXM_CELLTYPE_ENDO) then
!
!           vinf = 1.0_rp + exp((28.0_rp + elmag(ipoin,ITER_K)) / 5.0_rp)
!           vinf = 1.0_rp / vinf              ! Value of s_infinite
!
!           tau1 = 1000.0_rp * exp(-((67.0_rp + elmag(ipoin,ITER_K)) * (67.0_rp + elmag(ipoin,ITER_K))) / 1000.0_rp)
!           taux = (tau1 + 8.0_rp)*1.6_rp     ! Value of tau_s
!           xitaux = 1.0_rp / taux            ! Value of 1/tau_s 
!
!        else if (ituss_exm == EXM_CELLTYPE_MID) then
!
!           vinf = 1.0_rp + exp((20.0_rp + elmag(ipoin,ITER_K)) / 5.0_rp)
!           vinf = 1.0_rp / vinf              ! Value of s_infinite
!
!           tau1 = 85.0_rp * exp(-((45.0_rp + elmag(ipoin,ITER_K)) * (45.0_rp + elmag(ipoin,ITER_K))) / 320.0_rp)
!           tau2 = 1.0_rp + exp((-20.0_rp + elmag(ipoin,ITER_K)) / 5.0_rp)
!           tau2 = 5.0_rp / tau2
!           taux = (tau1 + tau2 + 3.0_rp)*1.6_rp       ! Value of tau_s
!           xitaux = 1.0_rp / taux            ! Value of 1/tau_s 
!
!        end if
!
!        vaux0 = vauxi_exm(8,ipoin,2)      ! Value of s in previous time step 
!
!        ! Scheme for numerical integration
!        if (kfl_odets_exm == 1) then      ! Runge - Kutta 4
!
!           vaux1 = (vinf - vaux0) * xitaux
!           vaux2 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux1)) * xitaux
!           vaux3 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux2)) * xitaux
!           vaux4 = (vinf - (vaux0 + dtimeEP * vaux3)) * xitaux
!           vaux5 = vaux0 + (dtimeEP / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4) 
!           vaux3 = vaux5
!
!        else if (kfl_odets_exm == 2) then      ! Crank - Nicholson
!
!           vaux1 = dtimeEP * xitaux
!           vaux2 = 1.0_rp + thnew * vaux1
!           vaux2 = 1.0_rp / vaux2
!           vaux3 = (1.0_rp - thold * vaux1) * vaux0 + vaux1 * vinf
!           vaux3 = vaux3 * vaux2       
!
!        else if (kfl_odets_exm == 3) then      ! Rush - Larsen
!
!           vaux3 = vinf - (vinf - vaux0) * exp(-dtimeEP/xmccmmate_exm(nodemat(ipoin))     * xitaux)    
!
!        end if
!
!        vauxi_exm(8,ipoin,1) = vaux3      ! Value of variable s
!
!        ! Variable of activation xs (related with I_Ks current)
!        ! nauxi_exm = 9
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1586 (col 2)
!
!        vinf = 1.0_rp + exp((-5.0_rp - elmag(ipoin,ITER_K)) / 14.0_rp)
!        vinf = 1.0_rp / vinf              ! Value of xs_infinite
!
!        alphx = 1.0_rp + exp((5.0_rp - elmag(ipoin,ITER_K)) / 6.0_rp)
!        alphx = 1400.0_rp / sqrt(alphx)   ! Value of alpha_xs
!
!        betax = 1.0_rp + exp((-35.0_rp + elmag(ipoin,ITER_K)) / 15.0_rp)
!        betax = 1.0_rp / betax            ! Value of beta_xs
!
!        taux = (alphx * betax) + 80.0_rp             ! Value of tau_xs
!        xitaux = 1.0_rp / taux            ! Value of 1/tau_xs
!
!        vaux0 = vauxi_exm(9,ipoin,2)      ! Value of xs in previous time step 
!
!        ! Scheme for numerical integration
!        if (kfl_odets_exm == 1) then      ! Runge - Kutta 4
!
!           vaux1 = (vinf - vaux0) * xitaux
!           vaux2 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux1)) * xitaux
!           vaux3 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux2)) * xitaux
!           vaux4 = (vinf - (vaux0 + dtimeEP * vaux3)) * xitaux
!           vaux5 = vaux0 + (dtimeEP / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4) 
!           vaux3 = vaux5
!
!        else if (kfl_odets_exm == 2) then      ! Crank - Nicholson
!
!           vaux1 = dtimeEP * xitaux
!           vaux2 = 1.0_rp + thnew * vaux1
!           vaux2 = 1.0_rp / vaux2
!           vaux3 = (1.0_rp - thold * vaux1) * vaux0 + vaux1 * vinf
!           vaux3 = vaux3 * vaux2       
!
!        else if (kfl_odets_exm == 3) then      ! Rush - Larsen
!
!           vaux3 = vinf - (vinf - vaux0) * exp(-dtimeEP/xmccmmate_exm(nodemat(ipoin))     * xitaux)    
!
!        end if
!
!        vauxi_exm(9,ipoin,1) = vaux3      ! Value of variable xs
!
!        ! Variable of activation xr1 (related with I_Kr1 current)
!        ! nauxi_exm = 10
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1586 (col 2)
!
!        vinf = 1.0_rp + exp((-26.0_rp - elmag(ipoin,ITER_K)) / 7.0_rp)
!        vinf = 1.0_rp / vinf              ! Value of xr1_infinite
!
!        alphx = 1.0_rp + exp((-45.0_rp - elmag(ipoin,ITER_K)) / 10.0_rp)
!        alphx = 450.0_rp / alphx          ! Value of alpha_xr1
!
!        betax = 1.0_rp + exp((30.0_rp + elmag(ipoin,ITER_K)) / 11.5_rp)
!        betax = 6.0_rp / betax            ! Value of beta_xr1
!
!        taux = alphx * betax              ! Value of tau_xr1
!        xitaux = 1.0_rp / taux            ! Value of 1/tau_xr1
!
!        vaux0 = vauxi_exm(10,ipoin,2)     ! Value of xr1 in previous time step 
!
!        ! Scheme for numerical integration
!        if (kfl_odets_exm == 1) then      ! Runge - Kutta 4
!
!           vaux1 = (vinf - vaux0) * xitaux
!           vaux2 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux1)) * xitaux
!           vaux3 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux2)) * xitaux
!           vaux4 = (vinf - (vaux0 + dtimeEP * vaux3)) * xitaux
!           vaux5 = vaux0 + (dtimeEP / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4) 
!           vaux3 = vaux5
!
!        else if (kfl_odets_exm == 2) then      ! Crank - Nicholson
!
!           vaux1 = dtimeEP * xitaux
!           vaux2 = 1.0_rp + thnew * vaux1
!           vaux2 = 1.0_rp / vaux2
!           vaux3 = (1.0_rp - thold * vaux1) * vaux0 + vaux1 * vinf
!           vaux3 = vaux3 * vaux2       
!
!        else if (kfl_odets_exm == 3) then      ! Rush - Larsen
!
!           vaux3 = vinf - (vinf - vaux0) * exp(-dtimeEP/xmccmmate_exm(nodemat(ipoin))     * xitaux)    
!
!        end if
!
!        vauxi_exm(10,ipoin,1) = vaux3     ! Value of variable xr1
!
!        ! Variable of activation xr2 (related with I_Kr2 current)
!        ! nauxi_exm = 11
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1586 (col 2)
!
!        vinf = 1.0_rp + exp((88.0_rp + elmag(ipoin,ITER_K)) / 24.0_rp)
!        vinf = 1.0_rp / vinf              ! Value of xr2_infinite
!
!        alphx = 1.0_rp + exp((-60.0_rp - elmag(ipoin,ITER_K)) / 20.0_rp)
!        alphx = 3.0_rp / alphx            ! Value of alpha_xr2
!
!        betax = 1.0_rp + exp((-60.0_rp + elmag(ipoin,ITER_K)) / 20.0_rp)
!        betax = 1.12_rp / betax           ! Value of beta_xr2
!
!        taux = alphx * betax              ! Value of tau_xr2
!        xitaux = 1.0_rp / taux            ! Value of 1/tau_xr2
!
!        vaux0 = vauxi_exm(11,ipoin,2)     ! Value of xr2 in previous time step 
!
!        ! Scheme for numerical integration
!        if (kfl_odets_exm == 1) then      ! Runge - Kutta 4
!
!           vaux1 = (vinf - vaux0) * xitaux
!           vaux2 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux1)) * xitaux
!           vaux3 = (vinf - (vaux0 + 0.5_rp * dtimeEP * vaux2)) * xitaux
!           vaux4 = (vinf - (vaux0 + dtimeEP * vaux3)) * xitaux
!           vaux5 = vaux0 + (dtimeEP/ 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4) 
!           vaux3 = vaux5
!
!        else if (kfl_odets_exm == 2) then      ! Crank - Nicholson
!
!           vaux1 = dtimeEP * xitaux
!           vaux2 = 1.0_rp + thnew * vaux1
!           vaux2 = 1.0_rp / vaux2
!           vaux3 = (1.0_rp - thold * vaux1) * vaux0 + vaux1 * vinf
!           vaux3 = vaux3 * vaux2       
!
!        else if (kfl_odets_exm == 3) then      ! Rush - Larsen
!
!           vaux3 = vinf - (vinf - vaux0) * exp(-dtimeEP/xmccmmate_exm(nodemat(ipoin))     * xitaux)    
!
!        end if
!
!        vauxi_exm(11,ipoin,1) = vaux3     ! Value of variable xr2               
!
!        ! Fast current of sodium I_Na [pA/pF]
!        ! nicel_exm = 1
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1585 (col 2)
!        enapot = groupcon * log (xtrnacon / vconc(3,ipoin,1))    ! Reversal potential of Na [mV]
!        vaux1 = vauxi_exm(1,ipoin,1) * vauxi_exm(1,ipoin,1) * vauxi_exm(1,ipoin,1)
!        vaux1 = nagmax * vaux1 * vauxi_exm(2,ipoin,1) * vauxi_exm(3,ipoin,1)
!        vaux2 = elmag(ipoin,ITER_K) - enapot
!        xina = vaux1 * vaux2  ! Variable xina
!
!        vicel_exm(1,ipoin,1) = xina       ! Value of I_Na current [pA/pF]
!
!
!        ! Current of L-Type calcium I_CaL   [pA/pF]
!        ! nicel_exm = 2
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1586 (col 1)
!        !i_CaL = ((( g_CaL*d*f*f2*fCass*4*(V - 15)*pow(F, 2))/(R*T))*(0.25*Ca_ss*(exp((( 2*(V - 15)*F)/(R*T)))) - Ca_o)) /
!        !        ((exp((( 2*(V - 15)*F)/(R*T)))) - 1)
!
!        vaux1 = calgmax * vauxi_exm(4,ipoin,1) * vauxi_exm (5,ipoin,1) * vauxi_exm(12,ipoin,1)
!        vaux1 = vaux1 * vauxi_exm(6,ipoin,1) * 4.0_rp * faracon * faracon * (elmag(ipoin,ITER_K) - 15.0_rp)
!        vaux1 = vaux1 / (rgascon * celltem)
!        vaux2 = exp((2.0_rp * (elmag(ipoin,ITER_K) - 15.0_rp) * faracon) / (rgascon * celltem))
!        vaux3 = (0.25_rp * vconc(5,ipoin,1) * vaux2) - xtrcacon
!        vaux3 = vaux1 * vaux3 
!        xical = vaux3 / (vaux2 - 1.0_rp) ! Variable xical
!
!        vicel_exm(2,ipoin,1) = xical       ! Value of I_CaL current [pA/pF]
!
!        ! Transient outward current I_to [pA/pF]
!        ! nicel_exm = 3
!        ! AND Slow delayed rectifier current I_Ks [pA/pF]
!        ! nicel_exm = 4
!        !    i_to = g_to*r*s*(V - E_K)
!        !    i_Ks = g_Ks*(pow(Xs,2))*(V - E_Ks)
!        ekpot = groupcon * log(xtrkcon / vconc(4,ipoin,1))    ! Reversal potential of K [mV]
!
!        vaux1 = xtrkcon + (ksnaperm * xtrnacon)
!        vaux2 = vconc(4,ipoin,1) + (ksnaperm * vconc(3,ipoin,1))               
!        ekspot = groupcon * log( vaux1 / vaux2 )    ! Reversal potential of Ks [mV]
!        ! Variable xito
!        vaux1 =  vauxi_exm(7,ipoin,1) * vauxi_exm(8,ipoin,1)!
!        vaux1 = vaux1 * (elmag(ipoin,ITER_K) - ekpot)
!        ! Variable xiks                    
!        vaux2 = vauxi_exm(9,ipoin,1) * vauxi_exm(9,ipoin,1)
!        vaux2 = vaux2 * (elmag(ipoin,ITER_K) - ekspot)
!
!        if(ituss_exm == EXM_CELLTYPE_EPI) then   !epicardial    
!           xito = vaux1 * togmaxepi
!           xiks = vaux2 * ksgmaxepi
!        else if(ituss_exm == EXM_CELLTYPE_ENDO) then !endocardial
!           xito = vaux1 * togmaxend
!           xiks = vaux2 * ksgmaxend
!        else if(ituss_exm == EXM_CELLTYPE_MID) then  !mid
!           xito = vaux1 * togmaxmce
!           xiks = vaux2 * ksgmaxmce
!        end if
!        vicel_exm(3,ipoin,1) = xito  ! Variable xIto
!        vicel_exm(4,ipoin,1) = xiks  ! Variable xIks
!
!        ! Rapid delayed rectifier current I_Kr [pA/pF]
!        ! nicel_exm = 5
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1586 (col 2)
!
!        vaux1 = krgmax * sqrt(xtrkcon / 5.4_rp) * vauxi_exm(10,ipoin,1) * vauxi_exm(11,ipoin,1)
!        xikr = vaux1 * (elmag(ipoin,ITER_K) - ekpot)
!        vicel_exm(5,ipoin,1) = xikr  ! Value of I_Kr current [pA/pF]
!
!        ! Inward rectifier potassium current I_K1 [pA/pF]
!        ! nicel_exm = 6
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1)
!        !  beta_K1 = ( 3*exp( 0.0002*(V - E_K+100))+(exp(( 0.1*((V - E_K) - 10))))) /
!        !   (1+(exp((-0.5*(V - E_K)))))
!        vaux1 = 3.0_rp * exp(0.0002_rp * (100.0_rp + elmag(ipoin,ITER_K) - ekpot))
!        vaux1 = vaux1 + exp(0.1_rp * (-10.0_rp + elmag(ipoin,ITER_K) - ekpot))
!        vaux2 = 1.0_rp + exp(-0.5_rp * (elmag(ipoin,ITER_K) - ekpot))
!        vaux2 = 1.0_rp / vaux2
!        vaux1 = vaux1 * vaux2                  ! Value of beta_K1
!        vaux2 = 1.0_rp + exp(0.06_rp * (elmag(ipoin,ITER_K) - ekpot - 200.0_rp))
!        vaux2 = 0.1_rp / vaux2                 ! Value of alpha_K1
!        vaux3 = vaux2 / (vaux1 + vaux2)        ! Value of K1_inf
!        !   i_K1 = g_K1*xK1_inf* pow((K_o/5.4), (1/2))*(V - E_K)     
!        vaux1 = sqrt(xtrkcon / 5.4_rp) * vaux3 * (elmag(ipoin,ITER_K) - ekpot)
!        xik1 = k1gmax * vaux1
!        vicel_exm(6,ipoin,1) = xik1       ! Value of I_K1 current [pA/pF]
!
!        ! Sodium/calcium exchanger current I_NaCa [pA/pF]
!        ! nicel_exm = 7
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1) 
!        !     xinaca = ( K_NaCa.*( (exp((( gamma_NaCa.*V.*F)./( R.*T)))).*(Nai .^ 3.0).*Cao -  (exp((( (gamma_NaCa - 1.0).*V.*F)./( R.*T)))).
!        !   *(Nao .^ 3.0).*Cai.*alpha_NaCa))./( ((Km_Nai .^ 3.0)+(Nao .^ 3.0)).*(Km_Ca+Cao).*(1.0+ K_sat.*(exp((( (gamma_NaCa - 1.0).*V.*F)./( R.*T))))));
!
!        vaux1 = exp((gamanaca-1.0_rp) * elmag(ipoin,ITER_K) * xgroupcon)
!        vaux2 = (kmnanaca * kmnanaca * kmnanaca) + (xtrnacon * xtrnacon * xtrnacon)
!        vaux2 = vaux2 * (kmcanaca + xtrcacon)
!        vaux2 = vaux2 * (1.0_rp + (satunaca * vaux1))     ! current denominator 
!        vaux2 = 1.0_rp / vaux2                          ! inverse of current denominator
!
!        vaux3 = -vaux1 * xtrnacon * xtrnacon * xtrnacon * vconc(1,ipoin,1) * alfanaca         
!        vaux1 = exp(gamanaca * elmag(ipoin,ITER_K) * xgroupcon)
!        vaux3 = vaux3 + (vaux1 * vconc(3,ipoin,1) * vconc(3,ipoin,1) * vconc(3,ipoin,1) * xtrcacon)          ! current numerator
!
!        xinaca = nacamax * vaux3 * vaux2 
!
!        vicel_exm(7,ipoin,1) = xinaca       ! Value of I_NaCa current [pA/pF]
!
!        ! Sodium/potassium pump current I_NaK [pA/pF]
!        ! nicel_exm = 8
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1)
!        !i_NaK = (((( P_NaK*K_o)/(K_o+K_mk))*Na_i)/(Na_i+K_mNa))/
!        !        (1+ 0.1245*(exp(((-0.1*V*F)/( R*T))))+ 0.0353*(exp(((-V*F)/( R*T))))) 
!        vaux1 = ((nakmax * xtrkcon) / (xtrkcon + kmknak)) * vconc(3,ipoin,1) 
!        vaux1 = vaux1 / (vconc(3,ipoin,1) + kmnanak)             ! current denominator 
!        vaux2 = 1.0_rp + 0.1245_rp * exp(-0.1_rp * elmag(ipoin,ITER_K) * xgroupcon)
!        vaux2 = vaux2 + 0.0353_rp * exp(-elmag(ipoin,ITER_K) * xgroupcon)
!        xinak =  vaux1 / vaux2 
!
!        vicel_exm(8,ipoin,1) = xinak       ! Value of I_NaK current [pA/pF]
!
!        ! Calcium plateau current I_pCa [pA/pF]
!        ! nicel_exm = 9
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1)
!
!        vaux1 = kpcapca + vconc(1,ipoin,1)
!        vaux1 = 1.0_rp / vaux1             ! inverse of current denominator
!        xipca = pcagmax * vaux1 * vconc(1,ipoin,1)
!
!        vicel_exm(9,ipoin,1) = xipca       ! Value of I_pCa current [pA/pF]
!
!        ! Potassium plateau current I_pK [pA/pF]
!        ! nicel_exm = 10
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1)
!        !    i_p_K = (g_pK*(V - E_K))/(1+(exp(((25 - V)/5.98))))
!        vaux1 = 1.0_rp + exp((25.0_rp - elmag(ipoin,ITER_K)) / 5.98_rp)
!        vaux1 = 1.0_rp / vaux1             ! inverse of current denominator
!        xipk = pkgmax * vaux1 * (elmag(ipoin,ITER_K) - ekpot)
!
!        vicel_exm(10,ipoin,1) = xipk       ! Value of I_pK current [pA/pF]
!
!        ! Sodium background current I_bNa [pA/pF]
!        ! nicel_exm = 11
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1)
!
!        xibna = bnagmax * (elmag(ipoin,ITER_K) - enapot)
!
!        vicel_exm(11,ipoin,1) = xibna      ! Value of I_bNa current [pA/pF]
!
!        ! Calcium background current I_bCa [pA/pF]
!        ! nicel_exm = 12
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1)
!
!        ecapot = 0.50_rp * groupcon  * log (xtrcacon / vconc(1,ipoin,1))    ! Reversal potential of Ca [mV]
!
!!!!!   JAZMIN VER ACAAAA     write (1300+ipoin,*) xtrcacon,vconc(1,ipoin,1)
!
!        xibca = bcagmax * (elmag(ipoin,ITER_K) - ecapot)
!
!        vicel_exm(12,ipoin,1) = xibca      ! Value of I_bCa current [pA/pF]
!
!        ! Leak current from sarcoplasmic reticulum (SR) to the cytoplasm I_leak [mM/ms]
!        ! nicel_exm = 13
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1)
!
!        xileak = leakvmax * (vconc(2,ipoin,1) - vconc(1,ipoin,1))
!
!        vicel_exm(13,ipoin,1) = xileak      ! Value of I_leak current [mM/ms]
!
!        ! Pump current taking up calcium in the sarcoplasmic reticulum (SR) I_up [mM/ms]
!        ! nicel_exm = 14
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1)
!        ! i_up = Vmax_up/(1+pow(K_up, 2)/pow(Ca_i, 2)) 
!        vaux1 = vconc(1,ipoin,1) * vconc(1,ipoin,1)
!        vaux1 = 1.0_rp / vaux1
!        vaux2 = 1.0_rp + (upkcon * upkcon * vaux1)
!        vaux2 = 1.0_rp / vaux2            ! inverse of current denominator
!        xiup = upvmax * vaux2
!
!        vicel_exm(14,ipoin,1) = xiup      ! Value of I_up current [mM/ms]
!
!        ! Calcium-induced calcium release (CICR) current I_rel [mM/ms]
!        ! nicel_exm = 15
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1)
!
!        xirel = 0.102_rp * vconc(10,ipoin,1) * (vconc(2,ipoin,1) - vconc(5,ipoin,1))  
!
!        vicel_exm(15,ipoin,1) = xirel      ! Value of I_rel current [mM/ms]
!
!        ! External stimulus current I_stim [pA/pF]
!        ! nicel_exm = 16
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1581
!
!        xistim = 0.0_rp !appfi_exm(ipoin)
!        vicel_exm(16,ipoin,1) = xistim      ! Value of I_stim current [pA/pF]
!
!        ! IMPORTANTE: ESTO (I_ax) EN REALIDAD NO SE DEBERIA DEFINIR ACA
!        !             HABRIA QUE INTRODUCIRLO DESDE AFUERA
!        !             SE USA PARA EL CALCULO DE LAS CONCENTRACIONES
!
!        ! Axial current flow I_ax [pA/pF]
!        ! nicel_exm = 17
!        ! "A model for ventricular tissue", Ten Tusscher, Page H1581
!
!        xiax = 0.0_rp
!
!        vicel_exm(17,ipoin,1) = xiax      ! Value of I_ax current [pA/pF]
!
!        !Variable xixfer
!        vicel_exm(18,ipoin,1) = 0.0038_rp * (vconc(5,ipoin,1) - vconc(1,ipoin,1))
!
!
!        xioni= 0.0_rp
!        xioni= 0.0_rp
!        do iicel=1,12
!           xioni= xioni + vicel_exm(iicel,ipoin,1)
!        end do
!
!        !!  Runge Kutta, now done outside
!        
!        !vicel_exm(:,ipoin,2) = vicel_exm(:,ipoin,1)
!        vauxi_exm(:,ipoin,3) = vauxi_exm(:,ipoin,2)       
!        vauxi_exm(:,ipoin,2) = vauxi_exm(:,ipoin,1)
!
!        !Moved to exm_updunk
!        !vconc(:,ipoin,3)     = vconc(:,ipoin,2)
!        !vconc(:,ipoin,2)     = vconc(:,ipoin,1)
!!/xmccmmate_exm(n)
!        !xioni = xioni * 0.000003335640952_rp
!        dioni = 0.0_rp
!     end if
!  end if
!  !write(994,*) elmag(1,:,2)
!end subroutine exm_tt2006_ionicurrents
