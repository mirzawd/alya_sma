!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_exm_torord_model
   use def_kintyp, only : ip, rp
   use def_master, only : intost
   use mod_eccoupling

   implicit none

   contains
!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_torord_model.f90
!> @author  Zhinuo Jenny Wang
!> @brief   The ToR-ORd model (2019)
!> @date    14/FEB/2020
!> @details A single implementation of the model which takes the state variables as input
!> @}
!!-----------------------------------------------------------------------

subroutine exm_torord_model(cell_type,U,tconc_exm,tgates_exm,ti_exm,dt,a2bas,U_n,flag_land,flag_3D,flag_border,flag_drug,GKr_scaling,drugd,statvar,ipoin)

   use      def_exmedi
   use      def_domain
   use      mod_eccoupling, only : get_active_traction_LANDS

   ! Definition of variables
   implicit none
   integer(ip), intent(in) :: cell_type
   real(rp), intent(in) :: U, dt, a2bas, GKr_scaling
   real(rp), intent(out) :: U_n
   logical, intent(in) :: flag_land,flag_3D,flag_border,flag_drug
   real(rp), intent(inout) :: tgates_exm(31,2),ti_exm(21,1),tconc_exm(14,2)
   real(rp), intent(inout), optional :: statvar(7,2)
   real(rp), intent(in) :: drugd(12)
   integer(ip), intent(in), optional :: ipoin
   real(rp) :: ENa, EK, EKs, ECl, Jrel, Jup, Jtr, Jdiff, JdiffNa, JdiffK, Jleak, Farad
!   real(rp) :: JdiffCl
   real(rp) ::    CaMKa, CaMKb, vffrt, vfrt, mss, tm, hss, th, jss, tj, hssp, tjp
   real(rp) :: GNa, fINap, mLss, tmL, hLss, thL, hLssp, thLp, fINaLp, ass, ta
   real(rp) :: iss_Ito, delta_epiIto, tiF, tiS, AiF, AiS, i_to, assp, dti_develop
   real(rp) ::     dti_recover, tiFp, tiSp, i_p, fItop, dss, td, fss, tff, tfs, Aff
   real(rp) :: Afs, f, fcass, tfcaf, tfcas, Afcaf, Afcas, fca, tjca, tffp, fp
   real(rp) :: tfcafp, fcap, Kmn, k2n, km2n, ancaiss, PhiCaL_i, PhiCaNa_i
   real(rp) ::     PhiCaK_i, zca, PCap, PCaNa, PCaK, PCaNap, PCaKp, fICaLp, xs1ss
   real(rp) :: txs1, xs2ss, txs2, KsCa, kna1, kna2, kna3, kasymm, wna, wca,wnaca
   real(rp) ::     kcaon, kcaoff, qna, qca, hca, hna, h1, h2, h3, h4, h5, h6, h7, h8
   real(rp) ::     h9, h10, h11, h12, k1, k2, k3p, k3pp, k3, k4p, k4pp, k4, k5, k6
   real(rp) ::     k7, k8, x1, x2, x3, x4, E1, E2, E3, E4, KmCaAct, allo, zna
   real(rp) ::     JncxNa, JncxCa, k1p, k1m, k2p, k2m, k3m, k4m, Knai0, Knao0, delta
   real(rp) ::   Knai, Knao, Kki, Kko, MgADP, MgATP, Kmgatp, HH, eP, Khp, Knap
   real(rp) ::     Kxkur, P, a1, b1, a2, b2, a3, b3, a4, b4, zk, JnakNa, JnakK, xkb
   real(rp) ::     GKb, PNab, PCab, GpCa, bt, a_rel, Jrel_inf, tau_rel, btp, a_relp
   real(rp) ::     Jrel_infp, tau_relp, fJrelp, Jupnp, Jupp, fJupp, Bcai, Bcass
   real(rp) ::     Bcajsr, zcl, clo, cli, ah, bh, aj, bj, dielConstant, Io, Iss
!   real(rp) ::     clss
   real(rp) :: constA, gamma_cass, gamma_cao, gamma_nass, gamma_nao
   real(rp) :: gamma_kss, gamma_ko, ancass, PhiCaL_ss,PhiCaNa_ss, PhiCaK_ss
   real(rp) :: ICaL_ss, ICaNa_ss, ICaK_ss, Ii, gamma_cai, gamma_nai, gamma_ki
   real(rp) :: ICaL_i, ICaNa_i, ICaK_i, alpha_1, beta_1, alpha, beta, alpha_2
   real(rp) :: beta_2, alpha_i, beta_i, alpha_C2ToI, beta_ItoC2, dC1, dC2, dC3
   real(rp) :: ddO_IKr, dI_IKr, aK1, bK1, K1ss, INaCa_fractionSS, GClb, Fjunc
   real(rp) :: GClCa, KdClCa, IClCa_junc, IClCa_sl, A_atp, K_atp, K_o_n, fkatp
   real(rp) :: gkatp, akik, bkik, Jrel_b, jcass, ICaL_fractionSS, tiF_b, tiS_b
   real(rp) :: ICaL_temp1, ICaL_temp2, Iion_temp, nao, cao, ko, temp
   real(rp) :: nai, nass, ki, kss, cai, cass, cansr, cajsr, CaMKt, KmCaMK, KmCaM, Jrel_n
   real(rp) :: Jrel_p, INa, INaL, Ito, ICaL, ICaNa, ICaK, IKr, IKs, IK1, INaCa_i, INaCa_ss, INaK, INab, ICab
   real(rp) :: IKb, IpCa, IClCa, IClb, I_katp, Istim, m, h, j, hp, jp, mL, hL, hLp, a, i_F, iS, iFp, iSp
   real(rp) :: d, ff, fs, fcaf, fcas, jca, ffp, fcafp, nca_ss, nca_i, C1, C2, C3, O_IKr, I_IKr, xs1, xs2
   real(rp) :: GNaL, Gto, Jrel_inf0, Jrel_infp0, upScale, cmdnmax, PCa, GKr, GKs, GK1, Gncx, Pnak
   real(rp) :: vcell, trpnmax, vjsr, vmyo, vnsr, vss, Acap, Ageo, ap, aCaMK, bCaMK, BSLmax, BSRmax,CaMKo
   real(rp) :: csqnmax,dU, Iion, KmBSL, KmBSR, kmcmdn, kmcsqn, kmtrpn,L,R,T,rad, INaCa
   real(rp) :: lambda,dCaTRPN,Ca50,k_TRPN,n_TRPN,cbracket
   real(rp) :: dcatrpndt,sac,tau_relp_scaling,aCaMK_scaling
   real(rp) :: dosis1,dosis2,dosis3,dosis4,i50cal,i50kr,i50na,i50ito,gcaldrug,gkrdrug,gnadrug,gitodrug
   real(rp) :: Tact,dTact
   real(rp), target :: props(17)

   ! Initialise all variables from input -----------------------------------------------------------
   Istim = ti_exm(20,1)

   ! Intracellular concentrations [mM]
   nai = tconc_exm(1,1)
   nass = tconc_exm(2,1)
   ki = tconc_exm(3,1)
   kss = tconc_exm(4,1)
   cai = tconc_exm(5,1)
   cass = tconc_exm(6,1)
   cansr = tconc_exm(7,1)
   cajsr = tconc_exm(8,1)
   cli = tconc_exm(9,1)
   !clss = tconc_exm(10,1)
   CaMKt = tconc_exm(10,1)
   KmCaMK = tconc_exm(11,1)
   KmCaM = tconc_exm(12,1)
   Jrel_n = tconc_exm(13,1)
   Jrel_p = tconc_exm(14,1)

   ! Gates
   m      = tgates_exm(1,1)
   h      = tgates_exm(2,1)
   j      = tgates_exm(3,1)
   hp     = tgates_exm(4,1)
   jp     = tgates_exm(5,1)
   mL     = tgates_exm(6,1)
   hL     = tgates_exm(7,1)
   hLp    = tgates_exm(8,1)
   a      = tgates_exm(9,1)
   i_F    = tgates_exm(10,1)
   iS     = tgates_exm(11,1)
   ap     = tgates_exm(12,1)
   iFp    = tgates_exm(13,1)
   iSp    = tgates_exm(14,1)
   d      = tgates_exm(15,1)
   ff     = tgates_exm(16,1)
   fs     = tgates_exm(17,1)
   fcaf   = tgates_exm(18,1)
   fcas   = tgates_exm(19,1)
   jca    = tgates_exm(20,1)
   ffp    = tgates_exm(21,1)
   fcafp  = tgates_exm(22,1)
   nca_ss = tgates_exm(23,1)
   nca_i  = tgates_exm(24,1)
   C1     = tgates_exm(25,1)
   C2     = tgates_exm(26,1)
   C3     = tgates_exm(27,1)
   O_IKr  = tgates_exm(28,1)
   I_IKr  = tgates_exm(29,1)
   xs1    = tgates_exm(30,1)
   xs2    = tgates_exm(31,1)

  ! Drug effects
   if (flag_drug) then
     dosis1 = drugd(1)
     i50cal = drugd(2)
     dosis2 = drugd(3)
     i50kr = drugd(4)
     dosis3 = drugd(5)
     i50na = drugd(6)
     dosis4 = drugd(7)
     i50ito = drugd(8)
     gcaldrug = 1.0_rp / (1.0_rp + (dosis1/i50cal))
     gkrdrug = 1.0_rp / (1.0_rp + (dosis2/i50kr))
     gnadrug = 1.0_rp / (1.0_rp + (dosis3/i50na))
     gitodrug = 1.0_rp / (1.0_rp + (dosis4/i50ito))
   else
     gcaldrug = 1.0_rp
     gkrdrug = 1.0_rp
     gnadrug = 1.0_rp
     gitodrug = 1.0_rp
   end if

   ! Ionic conductances that vary with cell type [mS/muF]
   GNaL = 0.0279_rp
   GNa = 11.7802_rp * gnadrug
   Gto = 0.16_rp * gitodrug
   Jrel_inf0 = 1.0_rp
   Jrel_infp0 = 1.0_rp
   upScale = 1.0_rp
   cmdnmax = 0.05_rp
   PCa = 8.3757e-5_rp * gcaldrug
   GKr = 0.0321_rp * gkrdrug
   GKs = 0.0011_rp * a2bas
   GK1 = 0.6992_rp
   Gncx = 0.0034_rp
   Pnak = 15.4509_rp
   PCab = 5.9194e-8_rp
   GKb = 0.0189_rp


   if(cell_type == 3) then !!epi
      GNaL = GNaL * 0.6_rp
      Gto = Gto * 2.0_rp
      Jrel_inf0 = Jrel_inf0
      Jrel_infp0 = Jrel_infp0
      upScale = upScale * 1.3_rp
      cmdnmax = cmdnmax * 1.3_rp
      PCa = PCa * 1.2_rp
      GKr = GKr * 1.3_rp
      GKs = GKs * 1.4_rp
      GK1 = GK1 * 1.2_rp
      Gncx = Gncx * 1.1_rp
      Pnak = Pnak * 0.9_rp
      GKb = GKb * 0.6_rp
   else if(cell_type == 2) then !!MID
      Gto = Gto * 2.0_rp
      Jrel_inf0 = Jrel_inf0 * 1.7_rp
      Jrel_infp0 = Jrel_infp0 * 1.7_rp
      PCa = PCa * 2.0_rp
      GKr = GKr * 0.8_rp
      GK1 = GK1 * 1.3_rp
      Gncx = Gncx * 1.4_rp
      Pnak = Pnak * 0.7_rp
   else if(cell_type == 1) then  !ENDO
   end if

   if (flag_border) then ! Infarction borderzone remodelling
      GNa = 0.4_rp * GNa
      Gto = 0.0_rp
      PCa = 0.64_rp * PCa
      GK1 = 0.6_rp * GK1
      PCab = 1.33_rp * PCab
      aCaMK_scaling = 1.5_rp
      tau_relp_scaling = 6.0_rp
      GKr = GKr_scaling * GKr
      if (GKr_scaling < 0.5_rp) then
         print *, 'gkr scaling', GKr_scaling
      end if
   else
      aCaMK_scaling = 1.0_rp
      tau_relp_scaling = 1.0_rp
   end if

   ! -----------------------------------------------------------------------------------------------
   ! Extracellular concentrations
   nao = 140.0_rp ! [mM]
   cao = 1.8_rp   ! [mM]
   ko = 5.0_rp    ! [mM]
   clo = 150.0_rp    ! [mM]
   ! Physical constants (if no units indicated then it is dimensionless)
   R = 8314.0_rp ! [J/kmol*K]
   T = 310.0_rp  ! [K]
   Farad = 96485.0_rp ! [coulomb/mol]
   zna = 1.0_rp
   zca = 2.0_rp
   zk = 1.0_rp
   zcl = -1.0_rp
   ! Cell geometry
   L = 0.01_rp    ! [cm]
   rad = 0.0011_rp ! [cm]
   vcell = 1000.0_rp * 3.14_rp * rad * rad * L ![mL]
   Ageo = 2.0_rp * 3.14_rp * rad * rad + 2.0_rp * 3.14_rp * rad * L ![cm^2]
   Acap = 2.0_rp * Ageo     ! [cm^2]
   ! Sub-compartments
   vmyo = 0.68_rp * vcell   ! [mL]
   vnsr = 0.0552_rp * vcell ! [mL]
   vjsr = 0.0048_rp * vcell ! [mL]
   vss = 0.02_rp * vcell     ! [mL]
   !--Updated CaMK----------------------------------------------------------------------------------
   aCaMK = 0.05_rp    ! [mM^-1 ms^-1]
   aCaMK = aCaMK_scaling * aCaMK
   bCaMK = 0.00068_rp ! [ms^-1]
   CaMKo = 0.05_rp    ! dimensionless
   CaMKb=CaMKo*(1.0_rp-CaMKt)/(1.0_rp+KmCaM/cass)
   CaMKa=CaMKb+CaMKt
   CaMKt=CaMKt+dt*(aCaMK*CaMKb*(CaMKb+CaMKt)-bCaMK*CaMKt)
   !--Reversal Potential----------------------------------------------------------------------------
   zcl=-1.0_rp
   clo=150.0_rp
   !cli=24.0_rp
   ENa=(R*T/(zna*Farad))*log(nao/nai)
   EK=(R*T/(zk*Farad))*log(ko/ki)
   EKs=(R*T/(zk*Farad))*log((ko+0.01833_rp*nao)/(ki+0.01833_rp*nai))
   ECl=R*T/(zcl*Farad)*log(clo/cli)
   !EClss=R*T/(zcl*Farad)*log(clo/clss)
   !--Physical features-----------------------------------------------------------------------------
   vffrt=U*Farad*Farad/(R*T)
   vfrt=U*Farad/(R*T)
   !! Intracellular concentrations
   !--Fast Sodium Current INa-----------------------------------------------------------------------
   mss = 1.0_rp/((1.0_rp+safe_exp(-(U+56.86_rp)/9.03_rp))*(1.0_rp+safe_exp(-(U+56.86_rp)/9.03_rp)))
   tm =0.1292_rp*safe_exp(-((U+45.79_rp)/15.54_rp)**2.0_rp)+0.06487_rp*safe_exp(-((U-4.823_rp)/51.12_rp)**2.0_rp)
   m=mss-(mss-m)*safe_exp(-dt/tm)
   hss = 1.0_rp/((1.0_rp+safe_exp((U+71.55_rp)/7.43_rp))*(1.0_rp+safe_exp((U+71.55_rp)/7.43_rp)))
   if (U >= -40.0_rp) then
     ah = 0.0_rp
   else
     ah = 0.057_rp*safe_exp(-(U+80.0_rp)/6.8_rp)
   endif
   if (U >= -40.0_rp) then
     bh = 0.77_rp/(0.13_rp*(1.0_rp+safe_exp(-(U+10.66_rp)/11.1_rp)))
   else
     bh = 2.7_rp*safe_exp(0.079_rp*U)+3.1e5_rp*safe_exp(0.3485_rp*U)
   endif
   th = 1.0_rp/(ah+bh)
   h=hss-(hss-h)*safe_exp(-dt/th)
   if (U >= -40.0_rp) then
     aj = 0.0_rp
   else
     aj = (-2.5428e4_rp*safe_exp(0.2444_rp*U)-6.948e-6_rp*safe_exp(-0.04391_rp*U))*(U+37.78_rp)/(1.0_rp+safe_exp(0.311_rp*(U+79.23_rp)))
   endif

   if (U >= -40.0_rp) then
     bj = 0.6_rp*safe_exp(0.057_rp*U)/(1.0_rp+safe_exp(-0.1_rp*(U+32.0_rp)))
   else
     bj = 0.02424_rp*safe_exp(-0.01052_rp*U)/(1.0_rp+safe_exp(-0.1378_rp*(U+40.14_rp)))
   endif
   jss=hss
   tj = 1.0_rp/(aj+bj)
   j=jss-(jss-j)*safe_exp(-dt/tj)
   hssp=1.0_rp/((1.0_rp+safe_exp((U+77.55_rp)/7.43_rp))*(1.0_rp+safe_exp((U+77.55_rp)/7.43_rp)))
   hp=hssp-(hssp-hp)*safe_exp(-dt/th)
   tjp=1.46_rp*tj
   jp=jss-(jss-jp)*safe_exp(-dt/tjp)
   fINap=(1.0_rp/(1.0_rp+KmCaMK/CaMKa))
   INa=GNa*(U-ENa)*m*m*m*((1.0_rp-fINap)*h*j+fINap*hp*jp)
   !--Late Sodium Current INaL----------------------------------------------------------------------
   mLss=1.0_rp/(1.0_rp+safe_exp((-(U+42.85_rp))/5.264_rp))
   tmL=tm
   mL=mLss-(mLss-mL)*safe_exp(-dt/tmL)
   hLss=1.0_rp/(1.0_rp+safe_exp((U+87.61_rp)/7.488_rp))
   thL=200.0_rp
   hL=hLss-(hLss-hL)*safe_exp(-dt/thL)
   hLssp=1.0_rp/(1.0_rp+safe_exp((U+93.81_rp)/7.488_rp))
   thLp=3.0_rp*thL
   hLp=hLssp-(hLssp-hLp)*safe_exp(-dt/thLp)
   fINaLp=(1.0_rp/(1.0_rp+KmCaMK/CaMKa))
   INaL=GNaL*(U-ENa)*mL*((1.0_rp-fINaLp)*hL+fINaLp*hLp)
   !--Transient Outward Potassium Current Ito-------------------------------------------------------
   ass=1.0_rp/(1.0_rp+safe_exp((-(U-14.34_rp))/14.82_rp))
   ta=1.0515_rp/(1.0_rp/(1.2089_rp*(1.0_rp+safe_exp(-(U-18.4099_rp)/29.3814_rp)))+3.5_rp/(1.0_rp+safe_exp((U+100.0_rp)/29.3814_rp)))
   a=ass-(ass-a)*safe_exp(-dt/ta)
   iss_Ito=1.0_rp/(1.0_rp+safe_exp((U+43.94_rp)/5.711_rp))
   if (cell_type == 3) then
      delta_epiIto = 1.0_rp - 0.95_rp / (1.0_rp + safe_exp((U + 70.0_rp)/5.0_rp))
   else
      delta_epiIto = 1.0_rp
   endif
   tiF_b=4.562_rp+1.0_rp/(0.3933_rp*safe_exp((-(U+100.0_rp))/100.0_rp)+0.08004_rp*safe_exp((U+50.0_rp)/16.59_rp))
   tiS_b=23.62_rp+1.0_rp/(0.001416_rp*safe_exp((-(U+96.52_rp))/59.05_rp)+1.780e-8_rp*safe_exp((U+114.1_rp)/8.079_rp))
   tiF=tiF_b*delta_epiIto
   tiS=tiS_b*delta_epiIto
   AiF=1.0_rp/(1.0_rp+safe_exp((U-213.6_rp)/151.2_rp))
   AiS=1.0_rp-AiF
   i_F=iss_Ito-(iss_Ito-i_F)*safe_exp(-dt/tiF)
   iS=iss_Ito-(iss_Ito-iS)*safe_exp(-dt/tiS)
   i_to=AiF*i_F+AiS*iS
   assp=1.0_rp/(1.0_rp+safe_exp((-(U-24.34_rp))/14.82_rp))
   ap=assp-(assp-ap)*safe_exp(-dt/ta)
   dti_develop=1.354_rp+1.0e-4_rp/(safe_exp((U-167.4_rp)/15.89_rp)+safe_exp(-(U-12.23_rp)/0.2154_rp))
   dti_recover=1.0_rp-0.5_rp/(1.0_rp+safe_exp((U+70.0_rp)/20.0_rp))
   tiFp=dti_develop*dti_recover*tiF
   tiSp=dti_develop*dti_recover*tiS
   iFp=iss_Ito-(iss_Ito-iFp)*safe_exp(-dt/tiFp)
   iSp=iss_Ito-(iss_Ito-iSp)*safe_exp(-dt/tiSp)
   i_p=AiF*iFp+AiS*iSp
   fItop=(1.0_rp/(1.0_rp+KmCaMK/CaMKa))
   Ito=Gto*(U-EK)*((1.0_rp-fItop)*a*i_to+fItop*ap*i_p)
   !--L type Calcium Current ICaL-------------------------------------------------------------------
   Aff=0.6_rp
   ICaL_fractionSS = 0.8_rp
   Kmn=0.002_rp
   dielConstant = 74.0_rp
   k2n=500.0_rp
   tjca=75.0_rp
   if (U >= 31.4978_rp) then
     dss = 1.0_rp
   else
     temp = -1.007_rp * safe_exp(-0.0829_rp*U)
     if (temp < -100.0_rp) then
        dss = 0.0_rp
     else
        dss = 1.0763_rp*safe_exp(temp)
     end if
   endif
   td=0.6_rp+1.0_rp/(safe_exp(-0.05_rp*(U+6.0_rp))+safe_exp(0.09_rp*(U+14.0_rp)))
   d=dss-(dss-d)*safe_exp(-dt/td)
   d=safe_value(d)
   fss=1.0_rp/(1.0_rp+safe_exp((U+19.58_rp)/3.696_rp))
   tff=7.0_rp+1.0_rp/(0.0045_rp*safe_exp(-(U+20.0_rp)/10.0_rp)+0.0045_rp*safe_exp((U+20.0_rp)/10.0_rp))
   tfs=1000.0_rp+1.0_rp/(0.000035_rp*safe_exp(-(U+5.0_rp)/4.0_rp)+0.000035_rp*safe_exp((U+5.0_rp)/6.0_rp))
   Afs=1.0_rp-Aff
   ff=fss-(fss-ff)*safe_exp(-dt/tff)
   fs=fss-(fss-fs)*safe_exp(-dt/tfs)
   f=Aff*ff+Afs*fs
   fcass=fss
   tfcaf=7.0_rp+1.0_rp/(0.04_rp*safe_exp(-(U-4.0_rp)/7.0_rp)+0.04_rp*safe_exp((U-4.0_rp)/7.0_rp))
   tfcas=100.0_rp+1.0_rp/(0.00012_rp*safe_exp(-U/3.0_rp)+0.00012_rp*safe_exp(U/7.0_rp))
   Afcaf=0.3_rp+0.6_rp/(1.0_rp+safe_exp((U-10.0_rp)/10.0_rp))
   Afcas=1.0_rp-Afcaf
   fcaf=fcass-(fcass-fcaf)*safe_exp(-dt/tfcaf)
   fcas=fcass-(fcass-fcas)*safe_exp(-dt/tfcas)
   fca=Afcaf*fcaf+Afcas*fcas
   jcass=1.0_rp/(1.0_rp+safe_exp((U+18.08_rp)/2.7916_rp))
   jca=jcass-(jcass-jca)*safe_exp(-dt/tjca)
   tffp=2.5_rp*tff
   ffp=fss-(fss-ffp)*safe_exp(-dt/tffp)
   fp=Aff*ffp+Afs*fs
   tfcafp=2.5_rp*tfcaf
   fcafp=fcass-(fcass-fcafp)*safe_exp(-dt/tfcafp)
   fcap=Afcaf*fcafp+Afcas*fcas
   km2n=jca*1.0_rp
   ancass=1.0_rp/(k2n/km2n+(1.0_rp+Kmn/cass)**4.0_rp)
   nca_ss=ancass*k2n/km2n-(ancass*k2n/km2n-nca_ss)*safe_exp(-km2n*dt)
   Ii = 0.5_rp*(nai+ki+cli+4.0_rp*cai)/1000.0_rp
   Iss = 0.5_rp*(nass+kss+cli+4.0_rp*cass)/1000.0_rp
   Io = 0.5_rp*(nao+ko+clo+4.0_rp*cao)/1000.0_rp
   constA = 1.82e6_rp*(dielConstant*T)**(-1.5_rp)
   gamma_cai = safe_exp(-constA*4.0_rp*(sqrt(Ii)/(1.0_rp+sqrt(Ii))-0.3_rp*Ii))
   gamma_cass = safe_exp(-constA*4.0_rp*(sqrt(Iss)/(1.0_rp+sqrt(Iss))-0.3_rp*Iss))
   gamma_cao = safe_exp(-constA*4.0_rp*(sqrt(Io)/(1.0_rp+sqrt(Io))-0.3_rp*Io))
   gamma_nai = safe_exp(-constA*1.0_rp*(sqrt(Ii)/(1.0_rp+sqrt(Ii))-0.3_rp*Ii))
   gamma_nass = safe_exp(-constA*1.0_rp*(sqrt(Iss)/(1.0_rp+sqrt(Iss))-0.3_rp*Iss))
   gamma_nao = safe_exp(-constA*1.0_rp*(sqrt(Io)/(1.0_rp+sqrt(Io))-0.3_rp*Io))
   gamma_ki = safe_exp(-constA*1.0_rp*(sqrt(Ii)/(1.0_rp+sqrt(Ii))-0.3_rp*Ii))
   gamma_kss = safe_exp(-constA*1.0_rp*(sqrt(Iss)/(1.0_rp+sqrt(Iss))-0.3_rp*Iss))
   gamma_ko = safe_exp(-constA*1.0_rp*(sqrt(Io)/(1.0_rp+sqrt(Io))-0.3_rp*Io))
   zca=2.0_rp
   PCap=1.1_rp*PCa
   PCaNa=0.00125_rp*PCa
   PCaK=3.574e-4_rp*PCa
   PCaNap=0.00125_rp*PCap
   PCaKp=3.574e-4_rp*PCap
   fICaLp=1.0_rp/(1.0_rp+KmCaMK/CaMKa)
   PhiCaL_ss=4.0_rp*vffrt*(gamma_cass*cass*safe_exp(2.0_rp*vfrt)-gamma_cao*cao)/(safe_exp(2.0_rp*vfrt)-1.0_rp)
   PhiCaNa_ss=1.0_rp*vffrt*(gamma_nass*nass*safe_exp(1.0_rp*vfrt)-gamma_nao*nao)/(safe_exp(1.0_rp*vfrt)-1.0_rp)
   PhiCaK_ss=1.0_rp*vffrt*(gamma_kss*kss*safe_exp(1.0_rp*vfrt)-gamma_ko*ko)/(safe_exp(1.0_rp*vfrt)-1.0_rp)
   ICaL_temp1=(1.0_rp-fICaLp)*PCa*PhiCaL_ss*d*(f*(1.0_rp-nca_ss)+jca*fca*nca_ss)
   ICaL_temp2=fICaLp*PCap*PhiCaL_ss*d*(fp*(1.0_rp-nca_ss)+jca*fcap*nca_ss)
   ICaL_ss = ICaL_fractionSS*(ICaL_temp1+ICaL_temp2)
   ICaL_temp1=(1.0_rp-fICaLp)*PCaNa*PhiCaNa_ss*d*(f*(1.0_rp-nca_ss)+jca*fca*nca_ss)
   ICaL_temp2=fICaLp*PCaNap*PhiCaNa_ss*d*(fp*(1.0_rp-nca_ss)+jca*fcap*nca_ss)
   ICaNa_ss = ICaL_fractionSS*(ICaL_temp1+ICaL_temp2)
   ICaL_temp1=(1.0_rp-fICaLp)*PCaK*PhiCaK_ss*d*(f*(1.0_rp-nca_ss)+jca*fca*nca_ss)
   ICaL_temp2=fICaLp*PCaKp*PhiCaK_ss*d*(fp*(1.0_rp-nca_ss)+jca*fcap*nca_ss)
   ICaK_ss = ICaL_fractionSS*(ICaL_temp1+ICaL_temp2)
   ICaK_ss = safe_value(ICaK_ss)
   !--
   ancaiss = 1.0_rp/(k2n/km2n+(1.0_rp+Kmn/cai)**4.0_rp)
   nca_i=ancaiss*k2n/km2n-(ancaiss*k2n/km2n-nca_i)*safe_exp(-km2n*dt)
   PhiCaL_i = 4.0_rp*vffrt*(gamma_cai*cai*safe_exp(2.0_rp*vfrt)-gamma_cao*cao)/(safe_exp(2.0_rp*vfrt)-1.0_rp)
   PhiCaNa_i = 1.0_rp*vffrt*(gamma_nai*nai*safe_exp(1.0_rp*vfrt)-gamma_nao*nao)/(safe_exp(1.0_rp*vfrt)-1.0_rp)
   PhiCaK_i = 1.0_rp*vffrt*(gamma_ki*ki*safe_exp(1.0_rp*vfrt)-gamma_ko*ko)/(safe_exp(1.0_rp*vfrt)-1.0_rp)
   ICaL_temp1=(1.0_rp-fICaLp)*PCa*PhiCaL_i*d*(f*(1.0_rp-nca_i)+jca*fca*nca_i)
   ICaL_temp2=fICaLp*PCap*PhiCaL_i*d*(fp*(1.0_rp-nca_i)+jca*fcap*nca_i)
   ICaL_i = (1.0_rp-ICaL_fractionSS)*(ICaL_temp1+ICaL_temp2)
   ICaL_temp1=(1.0_rp-fICaLp)*PCaNa*PhiCaNa_i*d*(f*(1.0_rp-nca_i)+jca*fca*nca_i)
   ICaL_temp2=fICaLp*PCaNap*PhiCaNa_i*d*(fp*(1.0_rp-nca_i)+jca*fcap*nca_i)
   ICaNa_i = (1.0_rp-ICaL_fractionSS)*(ICaL_temp1+ICaL_temp2)
   ICaL_temp1=(1.0_rp-fICaLp)*PCaK*PhiCaK_i*d*(f*(1.0_rp-nca_i)+jca*fca*nca_i)
   ICaL_temp2=fICaLp*PCaKp*PhiCaK_i*d*(fp*(1.0_rp-nca_i)+jca*fcap*nca_i)
   ICaK_i = (1.0_rp-ICaL_fractionSS)*(ICaL_temp1+ICaL_temp2)
   ICaL = ICaL_ss+ICaL_i
   ICaNa = ICaNa_ss+ICaNa_i
   ICaK = ICaK_ss+ICaK_i
   !--Rapid delayed Potassium Rectifier IKr---------------------------------------------------------
   alpha_1 = 0.154375_rp
   beta_1 = 0.1911_rp
   alpha = 0.1161_rp*safe_exp(0.299_rp*vfrt)
   beta = 0.2442_rp*safe_exp(-1.604_rp*vfrt)
   alpha_2 = 0.0578_rp*safe_exp(0.971_rp*vfrt)
   beta_2 = 0.349e-3_rp*safe_exp(-1.062_rp*vfrt)
   alpha_i = 0.2533_rp*safe_exp(0.5953_rp*vfrt)
   beta_i = 0.06525_rp*safe_exp(-0.8209_rp*vfrt)
   alpha_C2ToI = 0.52e-4_rp*safe_exp(1.525_rp*vfrt)
   beta_ItoC2 = beta_2*beta_i*alpha_C2ToI/(alpha_2*alpha_i)
   dC3 = beta*C2-alpha*C3
   dC2 = alpha*C3+beta_1*C1-(beta+alpha_1)*C2
   dC1 = alpha_1*C2+beta_2*O_IKr+beta_ItoC2*I_IKr-(beta_1+alpha_2+alpha_C2ToI)*C1
   ddO_IKr = alpha_2*C1+beta_i*I_IKr-(beta_2+alpha_i)*O_IKr
   dI_IKr = alpha_C2ToI*C1+alpha_i*O_IKr-(beta_ItoC2+beta_i)*I_IKr
   C1=C1+dt*dC1
   C2=C2+dt*dC2
   C3=C3+dt*dC3

   if (C1 < 0.0_rp) then
     C1 = 0.0_rp
     print *, 'C1 < 0'
   end if
   if (C2 <0.0_rp) then
     C2 = 0.0_rp
     print *, 'C2 < 0'
   end if
   if (C3 <0.0_rp) then
     C3 = 0.0_rp
     print * ,'C3 < 0'
   end if
   if (C1 > 1.0_rp) then
     C1 = 1.0_rp
     print *, 'C1 > 1'
   end if
   if (C2 > 1.0_rp) then
     C2 = 1.0_rp
     print *, 'C2 > 1'
   end if
   if (C3 > 1.0_rp) then
     C3 = 1.0_rp
     print *, 'C3 > 1'
   end if

   O_IKr=O_IKr+dt*ddO_IKr
   I_IKr=I_IKr+dt*dI_IKr
   IKr=GKr*sqrt(ko/5.0_rp)*O_IKr*(U-EK)
   !--Slow delayed Potassium Rectifier IKs----------------------------------------------------------
   xs1ss=1.0_rp/(1.0_rp+safe_exp((-(U+11.60_rp))/8.932_rp))
   txs1=817.3_rp+1.0_rp/(2.326e-4_rp*safe_exp((U+48.28_rp)/17.80_rp)+0.001292_rp*safe_exp((-(U+210.0_rp))/230.0_rp))
   xs1=xs1ss-(xs1ss-xs1)*safe_exp(-dt/txs1)
   xs2ss=xs1ss
   txs2=1.0_rp/(0.01_rp*safe_exp((U-50.0_rp)/20.0_rp)+0.0193_rp*safe_exp((-(U+66.54_rp))/31.0_rp))
   xs2=xs2ss-(xs2ss-xs2)*safe_exp(-dt/txs2)
   KsCa=1.0_rp+0.6_rp/(1.0_rp+(3.8e-5_rp/cai)**1.4_rp)
   IKs=GKs*KsCa*xs1*xs2*(U-EKs)
   !--Inward Potassium Rectifier IK1----------------------------------------------------------------
   aK1 = 4.094_rp/(1.0_rp+safe_exp(0.1217_rp*(U-EK-49.934_rp)))
   bK1 = (15.72_rp*safe_exp(0.0674_rp*(U-EK-3.257_rp))+safe_exp(0.0618_rp*(U-EK-594.31_rp)))/(1.0_rp+safe_exp(-0.1629_rp*(U-EK+14.207_rp)))
   K1ss = aK1/(aK1+bK1)
   IK1=GK1*sqrt(ko/5.0_rp)*K1ss*(U-EK)
   !--Sodium-Calcium Exchange current INCX ---------------------------------------------------------
   INaCa_fractionSS = 0.35_rp
   kna1=15.0_rp
   kna2=5.0_rp
   kna3=88.12_rp
   kasymm=12.5_rp
   wna=6.0e4_rp
   wca=6.0e4_rp
   wnaca=5.0e3_rp
   kcaon=1.5e6_rp
   kcaoff=5.0e3_rp
   qna=0.5224_rp
   qca=0.1670_rp
   hca=safe_exp(qca*vfrt)
   hna=safe_exp(qna*vfrt)
   h1=1.0_rp+nai/kna3*(1.0_rp+hna)
   h2=(nai*hna)/(kna3*h1)
   h3=1.0_rp/h1
   h4=1.0_rp+nai/kna1*(1+nai/kna2)
   h5=nai*nai/(h4*kna1*kna2)
   h6=1.0_rp/h4
   h7=1.0_rp+nao/kna3*(1.0_rp+1.0_rp/hna)
   h8=nao/(kna3*hna*h7)
   h9=1.0_rp/h7
   h10=kasymm+1.0_rp+nao/kna1*(1.0_rp+nao/kna2)
   h11=nao*nao/(h10*kna1*kna2)
   h12=1.0_rp/h10
   k1=h12*cao*kcaon
   k2=kcaoff
   k3p=h9*wca
   k3pp=h8*wnaca
   k3=k3p+k3pp
   k4p=h3*wca/hca
   k4pp=h2*wnaca
   k4=k4p+k4pp
   k5=kcaoff
   k6=h6*cai*kcaon
   k7=h5*h2*wna
   k8=h8*h11*wna
   x1=k2*k4*(k7+k6)+k5*k7*(k2+k3)
   x2=k1*k7*(k4+k5)+k4*k6*(k1+k8)
   x3=k1*k3*(k7+k6)+k8*k6*(k2+k3)
   x4=k2*k8*(k4+k5)+k3*k5*(k1+k8)
   E1=x1/(x1+x2+x3+x4)
   E2=x2/(x1+x2+x3+x4)
   E3=x3/(x1+x2+x3+x4)
   E4=x4/(x1+x2+x3+x4)
   KmCaAct=150.0e-06_rp
   allo=1.0_rp/(1.0_rp+(KmCaAct/cai)*(KmCaAct/cai))
   zna=1.0_rp
   JncxNa=3.0_rp*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp
   JncxCa=E2*k2-E1*k1
   INaCa_i=(1.0_rp-INaCa_fractionSS)*Gncx*allo*(zna*JncxNa+zca*JncxCa)
   !-
   h1=1.0_rp+nass/kna3*(1.0_rp+hna)
   h2=(nass*hna)/(kna3*h1)
   h3=1.0_rp/h1
   h4=1.0_rp+nass/kna1*(1.0_rp+nass/kna2)
   h5=nass*nass/(h4*kna1*kna2)
   h6=1.0_rp/h4
   h7=1.0_rp+nao/kna3*(1.0_rp+1.0_rp/hna)
   h8=nao/(kna3*hna*h7)
   h9=1.0_rp/h7
   h10=kasymm+1.0_rp+nao/kna1*(1.0_rp+nao/kna2)
   h11=nao*nao/(h10*kna1*kna2)
   h12=1.0_rp/h10
   k1=h12*cao*kcaon
   k2=kcaoff
   k3p=h9*wca
   k3pp=h8*wnaca
   k3=k3p+k3pp
   k4p=h3*wca/hca
   k4pp=h2*wnaca
   k4=k4p+k4pp
   k5=kcaoff
   k6=h6*cass*kcaon
   k7=h5*h2*wna
   k8=h8*h11*wna
   x1=k2*k4*(k7+k6)+k5*k7*(k2+k3)
   x2=k1*k7*(k4+k5)+k4*k6*(k1+k8)
   x3=k1*k3*(k7+k6)+k8*k6*(k2+k3)
   x4=k2*k8*(k4+k5)+k3*k5*(k1+k8)
   E1=x1/(x1+x2+x3+x4)
   E2=x2/(x1+x2+x3+x4)
   E3=x3/(x1+x2+x3+x4)
   E4=x4/(x1+x2+x3+x4)
   allo=1.0_rp/(1.0_rp+(KmCaAct/cass)*(KmCaAct/cass))
   JncxNa=3.0_rp*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp
   JncxCa=E2*k2-E1*k1
   INaCa_ss=INaCa_fractionSS*Gncx*allo*(zna*JncxNa+zca*JncxCa)
   INaCa=INaCa_i+INaCa_ss
   !--Sodium-Potassium Pump current INaK -----------------------------------------------------------
   k1p=949.5_rp
   k1m=182.4_rp
   k2p=687.2_rp
   k2m=39.4_rp
   k3p=1899.0_rp
   k3m=79300.0_rp
   k4p=639.0_rp
   k4m=40.0_rp
   Knai0=9.073_rp
   Knao0=27.78_rp
   delta=-0.1550_rp
   Knai=Knai0*safe_exp(delta*vfrt/3.0_rp)
   Knao=Knao0*safe_exp((1.0_rp-delta)*vfrt/3.0_rp)
   Kki=0.5_rp
   Kko=0.3582_rp
   MgADP=0.05_rp
   MgATP=9.8_rp
   Kmgatp=1.698e-7_rp
   HH=1.0e-7_rp
   eP=4.2_rp
   Khp=1.698e-7_rp
   Knap=224.0_rp
   Kxkur=292.0_rp
   P=eP/(1.0_rp+HH/Khp+nai/Knap+ki/Kxkur)
   a1=safe_value((k1p*(nai/Knai)**3.0_rp)/((1.0_rp+nai/Knai)**3.0_rp+(1.0_rp+ki/Kki)**2.0_rp-1.0_rp))
   b1=safe_value(k1m*MgADP)
   a2=safe_value(k2p)
   b2=safe_value((k2m*(nao/Knao)**3.0_rp)/((1.0_rp+nao/Knao)**3.0_rp+(1.0_rp+ko/Kko)**2.0_rp-1.0_rp))
   a3=safe_value((k3p*(ko/Kko)**2.0_rp)/((1.0_rp+nao/Knao)**3.0_rp+(1.0_rp+ko/Kko)**2.0_rp-1.0_rp))
   b3=safe_value((k3m*P*HH)/(1.0_rp+MgATP/Kmgatp))
   a4=safe_value((k4p*MgATP/Kmgatp)/(1.0_rp+MgATP/Kmgatp))
   b4=safe_value((k4m*(ki/Kki)**2.0_rp)/((1.0_rp+nai/Knai)**3.0_rp+(1.0_rp+ki/Kki)**2.0_rp-1.0_rp))
   x1=a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2
   x2=b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4
   x3=a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1
   x4=b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1
   E1=x1/(x1+x2+x3+x4)
   E2=x2/(x1+x2+x3+x4)
   E3=x3/(x1+x2+x3+x4)
   E4=x4/(x1+x2+x3+x4)
   zk=1.0_rp
   JnakNa=3.0_rp*(E1*a3-E2*b3)
   JnakK=2.0_rp*(E4*b1-E3*a1)
   INaK=Pnak*(zna*JnakNa+zk*JnakK)
   !-- IClCa & IClb --------------------------------------------------------------------------------
   GClb = 1.98e-3_rp
   Fjunc = 1.0_rp
   GClCa = 0.2843_rp
   KdClCa = 0.1_rp
   IClCa_junc = Fjunc*GClCa/(1.0_rp+KdClCa/cass)*(U-ECl)
   IClCa_sl = (1.0_rp-Fjunc)*GClCa/(1.0_rp+KdClCa/cai)*(U-ECl)
   IClCa = IClCa_junc+IClCa_sl
   IClb = GClb*(U-ECl)
   !-- Ikatp ---------------------------------------------------------------------------------------
   A_atp = 2.0_rp
   K_atp = 0.25_rp
   K_o_n = 5.0_rp
   fkatp = 0.0_rp
   gkatp = 4.3195_rp
   akik = (ko/K_o_n)**0.24_rp
   bkik = 1.0_rp/(1.0_rp+(A_atp/K_atp)*(A_atp/K_atp))
   I_katp = fkatp*gkatp*akik*bkik*(U-EK)
   !--Background currents: IKb, INab, ICab ---------------------------------------------------------
   xkb=1.0_rp/(1.0_rp+safe_exp(-(U-10.8968_rp)/23.9871_rp))
   IKb=GKb*xkb*(U-EK)
   PNab=1.9239e-9_rp
   INab=PNab*vffrt*(nai*safe_exp(vfrt)-nao)/(safe_exp(vfrt)-1.0_rp)
   ICab=PCab*4.0_rp*vffrt*(gamma_cai*cai*safe_exp(2.0_rp*vfrt)-gamma_cao*cao)/(safe_exp(2.0_rp*vfrt)-1.0_rp)
   GpCa=0.0005_rp
   IpCa=GpCa*cai/(0.0005_rp+cai)
   !--Stretch activated ionic current---------------------------------------------------------------
   sac = 0.0_rp
   !-- Overall currents ----------------------------------------------------------------------------
   Iion_temp=INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa+INaK+INab+IKb+IpCa+ICab
   Iion=Iion_temp+IClCa+IClb+I_katp+Istim+sac
   !-- Diffusion Fluxes ----------------------------------------------------------------------------
   JdiffNa=(nass-nai)/2.0_rp
   JdiffK=(kss-ki)/2.0_rp
   Jdiff=(cass-cai)/0.2_rp
   !-- RyR3 Ca Release------------------------------------------------------------------------------
   Jrel_b=1.5378_rp
   bt=4.75_rp
   a_rel=0.5_rp*bt
   tau_rel=max(bt/(1.0_rp+0.0123_rp/cajsr),0.001_rp)
   Jrel_inf=Jrel_inf0*a_rel*(-ICaL_ss)/(1.0_rp+(1.7_rp/cajsr)**8.0_rp)
   Jrel_n=Jrel_inf-(Jrel_inf-Jrel_n)*safe_exp(-dt/tau_rel)
   !Jrel_n=Jrel_n+dt*(Jrel_inf-Jrel_n)/tau_rel
   btp=1.25_rp*bt
   a_relp=0.5_rp*btp
   tau_relp=max(tau_relp_scaling*btp/(1.0_rp+0.0123_rp/cajsr),0.001_rp)
   Jrel_infp=Jrel_infp0*a_relp*(-ICaL_ss)/(1.0_rp+(1.7_rp/cajsr)**8.0_rp)
   Jrel_p=Jrel_infp-(Jrel_infp-Jrel_p)*safe_exp(-dt/tau_relp)
   !Jrel_p=Jrel_p+dt*(Jrel_infp-Jrel_p)/tau_relp
   fJrelp=1.0_rp/(1.0_rp+KmCaMK/CaMKa)
   Jrel=Jrel_b*((1.0_rp-fJrelp)*Jrel_n+fJrelp*Jrel_p)
   !-- Ca uptake via SERCA--------------------------------------------------------------------------
   Jupnp=upScale*0.005425_rp*cai/(cai+0.00092_rp)
   Jupp=upScale*2.75_rp*0.005425_rp*cai/(cai+0.00092_rp-0.00017_rp)
   fJupp=(1.0_rp/(1.0_rp+KmCaMK/CaMKa))
   Jleak=0.0048825_rp*cansr/15.0_rp
   Jup=(1.0_rp-fJupp)*Jupnp+fJupp*Jupp-Jleak
   !-- Tranlocation SR Fluxes ----------------------------------------------------------------------
   Jtr=(cansr-cajsr)/60.0_rp
   !-- Ionic Concentrations ------------------------------------------------------------------------
   !Na
   nai=nai+dt*(-(INa+INaL+3.0_rp*INaCa_i+ICaNa_i+3.0_rp*INaK+INab)*Acap/(Farad*vmyo)+JdiffNa*vss/vmyo)
   nass=nass+dt*(-(ICaNa_ss+3.0_rp*INaCa_ss)*Acap/(Farad*vss)-JdiffNa)
   !Ki
   ki=ki+dt*(-(Ito+IKr+IKs+IK1+IKb+I_katp+Istim-2.0_rp*INaK+ICaK_i)*Acap/(Farad*vmyo)+JdiffK*vss/vmyo)
   !print *, 'ki', ki
   !print *, 'Istim', Istim
   !print *, 'IKr', IKr
   kss=kss+dt*(-ICaK_ss*Acap/(Farad*vss)-JdiffK)
   !Cli
   !cli=cli+dt*((IClb+IClCa_sl)*Acap/(Farad*vmyo)+JdiffCl*vss/vmyo)
   !clss=clss+dt*(IClCa_junc*Acap/(Farad*vss)-JdiffCl)
   !Cai
   kmcmdn=0.00238_rp
   trpnmax=0.07_rp
   kmtrpn=0.0005_rp
   BSRmax=0.047_rp
   KmBSR=0.00087_rp
   BSLmax=1.124_rp
   KmBSL=0.0087_rp
   csqnmax=10.0_rp
   kmcsqn=0.8_rp
   cbracket = -(ICaL_i+IpCa+ICab-2.0_rp*INaCa_i)*Acap/(2.0_rp*Farad*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo
   ! Coupling with Land
   sac = 0.0_rp
   if (flag_land) then
      Bcai=1.0_rp/(1.0_rp+cmdnmax*kmcmdn/((kmcmdn+cai)**2.0_rp))

      if (flag_3D) then
         ! 3D case
         call exm_torland_calcium(ipoin,sac,dCaTRPN)
         ! No need to call Land model since SOLIDZ will handle it
      else
         ! Single cell case
         lambda = 1.0_rp
         statvar(7,:) = 1.0_rp ! Previous lambda
         ! Call Land model in single cell model
         props = [100.0_rp, 2.0_rp, 0.805_rp, 1000.0_rp, 5.0_rp, 0.35_rp, &
           & 182.0_rp, 12.0_rp, 0.5_rp, 0.25_rp, 0.0085_rp, 0.615_rp, 2.23_rp, &
           & 25.0_rp, 2.3_rp, -2.4_rp, 120.0_rp] ! Land model parameters
         call get_active_traction_LANDS(props,lambda,cai,Tact,dTact,statvar)
         dCaTRPN=k_TRPN*((cai/Ca50)**n_TRPN*(1.0_rp-statvar(3,2))-statvar(3,2))
      end if
      dcatrpndt=trpnmax*dCaTRPN
      cbracket=cbracket-dcatrpndt
   else
      Bcai=1.0_rp/(1.0_rp+cmdnmax*kmcmdn/((kmcmdn+cai)**2.0_rp)+trpnmax*kmtrpn/((kmtrpn+cai)**2.0_rp))
   end if
   cai=cai+dt*Bcai*cbracket
   Bcass=1.0_rp/(1.0_rp+BSRmax*KmBSR/((KmBSR+cass)**2.0_rp)+BSLmax*KmBSL/((KmBSL+cass)**2.0_rp))
   cass=cass+dt*(Bcass*(-(ICaL_ss-2.0_rp*INaCa_ss)*Acap/(2.0_rp*Farad*vss)+Jrel*vjsr/vss-Jdiff))
   cansr=cansr+dt*(Jup-Jtr*vjsr/vnsr)
   Bcajsr=1.0_rp/(1.0_rp+csqnmax*kmcsqn/((kmcsqn+cajsr)**2.0_rp))
   cajsr=cajsr+dt*(Bcajsr*(Jtr-Jrel))

   ! Evaluate updated membrane potential
   dU = -Iion
   U_n = U + dU * dt

   ! Intracellular concentrations [mM]
   tconc_exm(1,2) = safe_value(nai)
   tconc_exm(2,2) = safe_value(nass)
   tconc_exm(3,2) = safe_value(ki)
   tconc_exm(4,2) = safe_value(kss)
   tconc_exm(5,2) = safe_value(cai)
   tconc_exm(6,2) = safe_value(cass)
   tconc_exm(7,2) = safe_value(cansr)
   tconc_exm(8,2) = safe_value(cajsr)
   tconc_exm(9,2) = safe_value(cli)
   !tconc_exm(10,2) = safe_value(clss)
   tconc_exm(10,2) = safe_value(CaMKt)
   tconc_exm(11,2) = safe_value(KmCaMK)
   tconc_exm(12,2) = safe_value(KmCaM)
   tconc_exm(13,2) = safe_value(Jrel_n)
   tconc_exm(14,2) = safe_value(Jrel_p)

   ! Ionic currents [muA/muF]
   ti_exm(1,1) = safe_value(INa)
   ti_exm(2,1) = safe_value(INaL)
   ti_exm(3,1) = safe_value(Ito)
   ti_exm(4,1) = safe_value(ICaL)
   ti_exm(5,1) = safe_value(ICaNa)
   ti_exm(6,1) = safe_value(ICaK)
   ti_exm(7,1) = safe_value(IKr)
   ti_exm(8,1) = safe_value(IKs)
   ti_exm(9,1) = safe_value(IK1)
   ti_exm(10,1) = safe_value(INaCa_i)
   ti_exm(11,1) = safe_value(INaCa_ss)
   ti_exm(12,1) = safe_value(INaK)
   ti_exm(13,1) = safe_value(INab)
   ti_exm(14,1) = safe_value(IKb)
   ti_exm(15,1) = safe_value(ICab)
   ti_exm(16,1) = safe_value(IpCa)
   ti_exm(17,1) = safe_value(IClCa)
   ti_exm(18,1) = safe_value(IClb)
   ti_exm(19,1) = safe_value(I_katp)
   ti_exm(20,1) = safe_value(Istim)
   ti_exm(21,1) = safe_value(sac)

   ! Gates
   tgates_exm(1,2)= safe_value(m)
   tgates_exm(2,2)= safe_value(h)
   tgates_exm(3,2)= safe_value(j)
   tgates_exm(4,2)= safe_value(hp)
   tgates_exm(5,2)= safe_value(jp)
   tgates_exm(6,2)= safe_value(mL)
   tgates_exm(7,2)= safe_value(hL)
   tgates_exm(8,2)= safe_value(hLp)
   tgates_exm(9,2)= safe_value(a)
   tgates_exm(10,2)= safe_value(i_F)
   tgates_exm(11,2)= safe_value(iS)
   tgates_exm(12,2)= safe_value(ap)
   tgates_exm(13,2)= safe_value(iFp)
   tgates_exm(14,2)= safe_value(iSp)
   tgates_exm(15,2)= safe_value(d )
   tgates_exm(16,2)= safe_value(ff)
   tgates_exm(17,2)= safe_value(fs)
   tgates_exm(18,2)= safe_value(fcaf)
   tgates_exm(19,2)= safe_value(fcas)
   tgates_exm(20,2)= safe_value(jca )
   tgates_exm(21,2)= safe_value(ffp )
   tgates_exm(22,2)= safe_value(fcafp)
   tgates_exm(23,2)= safe_value(nca_ss)
   tgates_exm(24,2)= safe_value(nca_i )
   tgates_exm(25,2)= safe_value(C1 )
   tgates_exm(26,2)= safe_value(C2 )
   tgates_exm(27,2)= safe_value(C3 )
   tgates_exm(28,2)= safe_value(O_IKr)
   tgates_exm(29,2)= safe_value(I_IKr)
   tgates_exm(30,2)= safe_value(xs1 )
   tgates_exm(31,2)= safe_value(xs2 )

contains
#include "exm_safe_exp.f90.inc"

pure function safe_value(x) result(r)
   use def_kintyp, only: rp
   implicit none
   real(rp), intent(in) :: x
   real(rp)             :: r
   real(rp), parameter  :: tol = 1.0e-20_rp
   r = sign(max(abs(x),tol),x)
end function safe_value

end subroutine exm_torord_model
end module mod_exm_torord_model
