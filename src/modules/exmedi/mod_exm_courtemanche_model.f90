!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Exmedi
!>
!> @file    exm_courtemanche_ionicurrents.f90
!> @author  Eva Casoni Rero
!> @brief   Multiple cell run for Initial condition setup for Courtemanche  heterogeneous model
!> @date    29/MAR/2020
!> @details 
!> @
!!-----------------------------------------------------------------------

!!#############################
    !voltage: elmag
    
    !Ca_i:vconc(1,ipoin,1:3)
    !Ca_rel:vconc(2,ipoin,1:3)  
    !Ca_up:vconc(3,ipoin,1:3)
    !K_i:vconc(4,ioin,1:3)
    !Na_i:vconc(5,ipoin,1:3)  
      
    !Gate variables
    !u:vauxi_exm(1,:)  ;  v:vauxi_exm(2,:)  ;     !w:vauxi_exm(3,:)
    !d:vauxi_exm(4,:)  ;  f_Ca:vauxi_exm(5,:) ;   !f:vauxi_exm(6,:)
    !h:vauxi_exm(7,:)  ;  j:vauxi_exm(8,:)    ;   !m:vauxi_exm(9,:)
    !xr:vauxi_exm(10,:)  ;  xs:vauxi_exm(11,:)  ;  !oa:vauxi_exm(12,:)
    !oi:vauxi_exm(13,:)  ;  ua:vauxi_exm(14,:)  ;  !ui: vauxi_exm(15,:)
    
    !I_Na:vicel_exm(1,:)  ;  I_K1:vicel_exm(2,:)  ;  !I_to:vicel_exm(3,:)
    !I_kur:vicel_exm(4,:)  ;  I_Kr:vicel_exm(5,:)  ;  !I_Ks:vicel_exm(6,:)
    !I_Ca_L:vicel_exm(7,:)  ;  I_NaK:vicel_exm(8,:)  ;  !I_NaCa:vicel_exm(9,:)
    !I_CaP:vicel_exm(10,:)  ;  I_B_Na:vicel_exm(11,:)  ;  !I_B_Ca:vicel_exm(12,:)
    !I_B_K:vicel_exm(13,:)  ;  I_rel:vicel_exm(14,:)  ;  I_tr:vicel_exm(15,:)
    !I_up_leak:vicel_exm(16,:)  ;  I_up:vicel_exm(17,:)  ;  I_stim:vicel_exm(18,:)
!!#############################
!! main code
!!#############################


module mod_exm_courtemanche_model
   use def_kintyp, only : ip, rp
   use def_master, only : intost
   use mod_eccoupling

   implicit none

   contains
   

!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_torord_model.f90
!> @author  Eva Casoni
!> @brief   The Courtemanche model (1998)
!> @date    29/MAR/2021
!> @details A single implementation of the model which takes the state variables as input
!> @}
!------------------------------------------------------------------------

!subroutine exm_courtemanche_model(cell_type,elmag(ipoin,ITER_K),vconc_exm(:,ipoin,:),vauxi_exm(:,ipoin,:),vicel_exm(:,ipoin),dtimeEP,qneto_exm(ipoin),flag_land,flag_3D,flag_border,statvar,ipoin)

subroutine exm_courtemanche_model(cell_type,U,tconc_exm,tgates_exm,ti_exm,dt,U_n,flag_land,flag_3D,flag_isac,statvar,ipoin)

  use      def_parame
  use      def_domain
  use      def_exmedi
  use      mod_eccoupling


  integer(ip), intent(in) :: cell_type
  real(rp), intent(in) :: dt
  real(rp), intent(out) :: U_n
  logical, intent(in) :: flag_land,flag_3D, flag_isac
  real(rp), intent(inout) :: tgates_exm(15,2),ti_exm(18),tconc_exm(5,2), U(2)
  real(rp), intent(inout), optional :: statvar(7,2)
  integer(ip), intent(in), optional :: ipoin

!  integer(ip) :: imate
  real(rp)    :: Iion, dU
  real(rp)    :: Ca_CMDN, Ca_CSQN, Ca_TRPN   ! Ca buffers (millimolar)
  real(rp)    :: i_up_leak                   ! in Ca_leak_current_by the NSR (millimolar per millisecond)
  real(rp)    :: tau_u, tau_v, tau_w         ! in Ca_relseare_current_from_JSR (millisecond)
  real(rp)    :: u_infty, v_infty, w_infty   ! (dimensionless)
  real(rp)    :: Fn                          ! in Ca_release_current_from_JSR (dimensionless)
  real(rp)    :: i_rel, i_up                 ! in Ca_release_current_JSR and Ca_uptake_current_NSR (millimolar per millisecond)
  real(rp)    :: tau_d, tau_f_Ca, tau_f      ! in L_type_Ca_channel_d_gate, L_type_Ca_channel_f_Ca_gate, L_type_channel_f_gate (millisecond)
  real(rp)    :: d_infty, f_Ca_infty, f_infty! in L_type_Ca_channel_d_gate, L_type_Ca_channel_f_Ca_gate, L_type_channer_f_gate (dimensionless)
  real(rp)    :: i_Ca_L, i_NaCa              ! in L_type_Ca_channel and L_type_Ca_exchanger_current (picoA)
  real(rp)    :: E_Ca                        ! in background_currents (millivolt)
  real(rp)    :: i_B_Ca, i_B_K, i_B_Na       ! in background_currents (picoA)
  real(rp)    :: alpha_h, beta_h             ! in fast_sodium_current_h_gate (per _millisecond)
  real(rp)    :: alpha_j, beta_j             ! in fast_sodium_current_j_gate (per _millisecond)
  real(rp)    :: alpha_m, beta_m             ! in fast_sodium_current_m_gate (per _millisecond)
  real(rp)    :: h_inf, j_inf, m_inf         ! in fast_sodium_current_(h,j,m)_gate (dimensionless)
  real(rp)    :: tau_h, tau_j, tau_m         ! in fast_sodium_current_(h,j,m)_gate (millisecond)
  real(rp)    :: E_Na                        ! in fast_sodium_current (millivolt)
  real(rp)    :: i_Na                        ! in fast_sodium_current (picoA)
  real(rp)    :: B1, B2                      ! in intracellular_ion_concentratrion (millimolar_per_millisecond, dimensionless)
  real(rp)    :: i_st                        ! in membrane (picoA)
  real(rp)    :: alpha_xr, beta_xr, tau_xr   ! in rapid_delayed_rectifier_K_current_xr_gate (per _millisecond, per_millisecond, milliseond)
  real(rp)    :: xr_infty                    ! in rapid_delayed_rectifier_K_current_xr_gate (dimensionless)
  real(rp)    :: i_Kr, i_CaP                 ! in rapid_delayed rectifier_K_current (picoA), in rapid_sarcolemnal_calcium_pump_currend (picoA)
  real(rp)    :: alpha_xs, beta_xs, tau_xs   ! in slow_delayed_rectifier_K_current_xs_gate (per_millisecond, per_millisecond, millisecond)
  real(rp)    :: xs_infty                    ! in slow_delayed_rectifier_K_current_xs_gate (dimensionless)
  real(rp)    :: i_Ks                        ! in slow_elayed_rectifier_K_current (picoA)
  real(rp)    :: f_NaK, i_NaK                ! in sodium_potassium_pump (dimensionless, picoA)
  real(rp)    :: sigma                       ! in sodium_potassium_pump (dimensionless)
  real(rp)    :: E_K, i_K1                   ! in time_intdependent_potassium_current (millivolt, picoA)
  real(rp)    :: i_tr                        ! in transfer_current_from_NSR_to_JSR (millimolar_per_millisecond)
  real(rp)    :: alpha_oa, beta_oa, tau_oa   ! in transient_outward_K_current_oa_gate (per_millisecond)
  real(rp)    :: oa_infty                    ! in transiend_outward_K_current_oa_gate (dimensionless)
  real(rp)    :: alpha_oi, beta_oi, tau_oi   ! in transient_outward_K_current_oi_gate (per_millisecond, per_millisecond, millisecond)
  real(rp)    :: oi_infty                    ! in transient_outward_K_current_oi_gate (dimensionless)
  real(rp)    :: i_to                        ! in transient_outward_K_current (picoA)
  real(rp)    :: alpha_ua, beta_ua, tau_ua   ! in ultrarapid_delayed_rectifier_K_current_ua_gate (per_millisec, per_millisec, millisec)
  real(rp)    :: ua_infty                    ! in ultrarapid_delayed_rectifier_K_current_ua_gate (dimensionless)
  real(rp)    :: alpha_ui, beta_ui, tau_ui   ! in ultrarapid_delayed_rectifier_K_current_ui_gate (per_millisec, per_millisec, millisec)
  real(rp)    :: ui_infty                    ! in ultrarapid_delayed_rectifier_K_current_ui_gate (dimensionless)
  real(rp)    :: g_Kur, i_Kur                ! in ultrarapid_delayed_rectifier_K_current (nanoS_per_picoF, picoA)
  real(rp)    :: g_sac                       ! conductance of channel for ISAC
  real(rp)    :: A_Na, B_Na, C_Na            ! in I_SAC_Na
  real(rp)    :: A_K, B_K, C_K               ! in I_SAC_K
  real(rp)    :: A_Ca, B_Ca, C_Ca            ! in I_SAC_Ca
  real(rp)    :: I_SAC_Na, I_SAC_Ca, I_SAC_K ! currents for ISAC
  real(rp)    :: I_SAC
  real(rp)    :: a, b

  real(rp), parameter :: Cm = 100.0_rp           ! in membrane (millimolar)
  real(rp), parameter :: CMDN_max = 0.05_rp      ! in Ca_buffers (millimolar)
  real(rp), parameter :: CSQN_max = 10.0_rp      ! in Ca_buffers (millimolar)
  real(rp), parameter :: Km_CMDN = 0.00238_rp    ! in Ca_buffers (millimolar)
  real(rp), parameter :: Km_CSQN = 0.8_rp        ! in Ca_buffers (millimolar)
  real(rp), parameter :: Km_TRPN = 0.0005_rp     ! in Ca_buffers (millimolar)
  real(rp), parameter :: TRPN_max = 0.07_rp      ! in Ca_buffers (millimolar)
  real(rp), parameter :: Ca_up_max = 15.0_rp     ! in Ca_leak_current_by_the_NSR (millimolar_per_millisecond)
  real(rp), parameter :: K_rel = 30.0_rp         ! in Ca_release_current_from_JSR (per_millisecond -1)
  real(rp), parameter :: I_up_max = 0.005_rp     ! in Ca_uptake_current_by_the_NSR (millimolar_per_millisecond)
  real(rp), parameter :: K_up = 0.00092_rp       ! in Ca_uptake_current_by_the_NSR (millimolar)
  real(rp), parameter :: I_NaCa_max = 1600.0_rp  ! in Na_Ca_exchanger_current (picoA_per_picoF)
  real(rp), parameter :: K_mCa = 1.38_rp         ! in Na_Ca_exchanger_current (millmolar)
  real(rp), parameter :: K_mNa = 87.5_rp         ! in Na_Ca_exchanger_current (millimolar)
  real(rp), parameter :: K_sat = 0.1_rp          ! in Na_Ca_exchanger_current (dimensionless)
  real(rp), parameter :: gama = 0.35_rp          ! in Na_Ca_exchanger_current (dimensionless)
  real(rp), parameter :: g_B_Ca = 0.001131_rp    ! in backbround_currents (nanoS_per_picoF)
  real(rp), parameter :: g_B_K = 0.0_rp          ! in background_currents (nanoS_per_picoF)
  real(rp), parameter :: g_B_Na = 0.0006744375_rp  ! in background_currents (nanoS_per_picoF)
  real(rp), parameter :: g_Na = 7.8_rp           ! in fast_sodium_current (nanoS_per_picoF)
  real(rp), parameter :: V_cell = 20100.0_rp     ! in intracellular_ion_concentrations (micrometre_3), 
  real(rp), parameter :: V_i = 13668.0_rp        ! in intracellular_ion_concentratrion (micrometre_3), V_cell*0.068
  real(rp), parameter :: V_up = 1109.52_rp       ! in intracellular_ion_concentratrion (micrometre_3), V_cell*0.0552
  real(rp), parameter :: V_rel = 96.48_rp        ! in intracellular_ion_concentratrion (micrometre_3), V_cell*0.0048
  real(rp), parameter :: F = 96.4867_rp          ! in membrane (coulomb_per_millimole)
  real(rp), parameter :: R = 8.3143_rp           ! in membrane (joule_per_mole_kelvin)
  real(rp), parameter :: T = 310.0_rp            ! in membrane (kelvin)
  real(rp), parameter :: g_Kr_0 = 0.0294_rp      ! in rapid_delayed_rectifier_K_current (nanoS_per_picoF)
  real(rp), parameter :: i_CaP_max = 0.275_rp    ! in sarcolemmal_calcium_pump_current (picoA_per_picoF)
  real(rp), parameter :: g_Ks_0 = 0.12941176_rp  ! in slow_delayed_rectifier_K_current (nanoS_per_picoF)
  real(rp), parameter :: Km_K_o = 1.5_rp         ! in sodium_potassium_pump (millimolar)
  real(rp), parameter :: Km_Na_i = 10.0_rp       ! in sodium_potassium_pump (millimolar)
  real(rp), parameter :: i_NaK_max = 0.59933874_rp ! in sodium_potassium_pump (picoA_per_picoF)
  real(rp), parameter :: Ca_o = 1.8_rp           ! in standard_ionic_concentrations (millimolar)
  real(rp), parameter :: K_o = 5.4_rp            ! in standard_ionic_concentrations (millimolar)
  real(rp), parameter :: Na_o = 140.0_rp         ! in standard_ionic_concentrations (millimolar)
  real(rp), parameter :: g_K1_0 = 0.09_rp        ! in time_independent_potassium_current (nanoS_per_picoF)
  real(rp), parameter :: tau_tr = 180.0_rp       ! in transfer_current_from_NSR_to_JSR (millisecond)
  real(rp), parameter :: K_Q10 = 3.0_rp          ! in transient_outward_K_current (dimensionless)
  real(rp), parameter :: g_to_0 = 0.1652_rp      ! in transient_outward_K_current (nanoS_per_picoF)
  real(rp), parameter :: g_Ca_L_0 = 0.12375_rp   ! in L_type_Ca_channel (nanoS_per_picoF)
  ! these are for ISAC
  real(rp), parameter :: K_sac = 100_rp          ! amount of current when cell is not stretch
  real(rp), parameter :: GG_sac = 0.05_rp         ! maximum conductance of membrane
  real(rp), parameter :: alpha_sac = 3.0_rp      ! non-lineal dependence of stretch
  real(rp), parameter :: lambda = 1.0_rp         ! ratio of stretch
  real(rp), parameter :: P_na = 1.0_rp
  real(rp), parameter :: P_k = 1.0_rp
  real(rp), parameter :: P_ca = 1.0_rp
  real(rp), parameter :: Z_na = 1.0_rp
  real(rp), parameter :: Z_k = 1.0_rp
  real(rp), parameter :: Z_ca = 2.0_rp

  
  real(rp) :: g_Kr, g_Ks, g_K1, g_to, g_Ca_L                    ! currents to be modified by the heterogenity of the region
  real(rp) :: vaux1, vaux3, vaux4, rhsx1, rhsx2, rhsx, dCaTRPN  ! auxiliar variables for Calcium-Land coupling
  real(rp) :: Ca50_ref, Ca50_fact, n_TRPN                       ! parameters for the Land model
  real(rp) :: Y(21), dY(21)     ! membrance currents: vauxi_exm, vconc, elmag
  integer(ip) :: i


   ! Activate regions
   select case(cell_type)
 
     case(1_ip)
         g_K1 = 1.0_rp*g_K1_0
         g_to = 1.0_rp*g_to_0
         g_Kr = 1.0_rp*g_Kr_0
         g_Ks = 1.0_rp*g_Ks_0
         g_Ca_L = 1.0_rp*g_Ca_L_0

     case(2_ip)      
         g_K1 = 1.0_rp*g_K1_0
         g_to = 1.0_rp*g_to_0
         g_Kr = 1.0_rp*g_Kr_0
         g_Ks = 1.0_rp*g_Ks_0
         g_Ca_L = 1.67_rp*g_Ca_L_0

     case(3_ip)
         g_K1 = 1.0_rp*g_K1_0
         g_to = 1.0_rp*g_to_0
         g_Kr = 1.6_rp*g_Kr_0
         g_Ks = 1.0_rp*g_Ks_0
         g_Ca_L = 1.67_rp*g_Ca_L_0
   
     case(4_ip)
         g_K1 = 1.0_rp*g_K1_0
         g_to = 1.0_rp*g_to_0
         g_Kr = 1.3_rp*g_Kr_0
         g_Ks = 1.0_rp*g_Ks_0
         g_Ca_L = 0.8_rp*g_Ca_L_0

     case(5_ip)
         g_K1 = 1.0_rp*g_K1_0
         g_to = 1.0_rp*g_to_0
         g_Kr = 2.6_rp*g_Kr_0
         g_Ks = 1.0_rp*g_Ks_0
         g_Ca_L = 0.8_rp*g_Ca_L_0

     case(6_ip)
         g_K1 = 1.0_rp*g_K1_0
         g_to = 0.68_rp*g_to_0
         g_Kr = 1.1_rp*g_Kr_0
         g_Ks = 1.0_rp*g_Ks_0
         g_Ca_L = 0.9_rp*g_Ca_L_0

     case(7_ip)
         g_K1 = 1.0_rp*g_K1_0
         g_to = 1.0_rp*g_to_0
         g_Kr = 2.2_rp*g_Kr_0
         g_Ks = 1.0_rp*g_Ks_0
         g_Ca_L = 0.9_rp*g_Ca_L_0

     case(8_ip)
         g_K1 = 1.0_rp*g_K1_0
         g_to = 1.0_rp*g_to_0
         g_Kr = 2.0_rp*g_Kr_0
         g_Ks = 1.0_rp*g_Ks_0
         g_Ca_L = 0.9_rp*g_Ca_L_0

     case(9_ip)
         g_K1 = 0.9_rp*g_K1_0
         g_to = 0.9_rp*g_to_0
         g_Kr = 2.5_rp*g_Kr_0
         g_Ks = 1.9_rp*g_Ks_0
         g_Ca_L = 0.8_rp*g_Ca_L_0

    end select


     ! Declare some parameters of the Land model
     Ca50_ref = 0.805_rp !props_ecc(3 , imate) !0.805_rp
     Ca50_fact = 1.0_rp  !props_ecc(20, imate)!1.0_rp by defaults
     n_TRPN = 2.0_rp     

     ! Get membrane currents
     Y(1) = tgates_exm(1,1)
     Y(2) = tgates_exm(2,1)
     Y(3) = tgates_exm(3,1)
     Y(4) = tgates_exm(4,1)
     Y(5) = tgates_exm(5,1)
     Y(6) = tgates_exm(6,1)
     Y(7) = tgates_exm(7,1)
     Y(8) = tgates_exm(8,1)
     Y(9) = tgates_exm(9,1)
     Y(10) = tconc_exm(1,1)
     Y(11) = tconc_exm(2,1)
     Y(12) = tconc_exm(3,1)
     Y(13) = tconc_exm(4,1)
     Y(14) = tconc_exm(5,1)
     Y(15) = U(1)
     Y(16) = tgates_exm(10,1)
     Y(17) = tgates_exm(11,1)
     Y(18) = tgates_exm(12,1)
     Y(19) = tgates_exm(13,1)
     Y(20) = tgates_exm(14,1)
     Y(21) = tgates_exm(15,1)

     i_st = ti_exm(18)

         !************************************************************************************************** 
         ! HERE STARTS THE CODE OF COURTEMANCHE MODEL (matlab based routine, based on the article cited above)
         !************************************************************************************************** 
        
         ! Equilibrium potentials
         E_Na = R*T/F*log(Na_o/Y(14))           !Eq.28
         E_Ca = R*T/(2.0_rp*F)*log(Ca_o/Y(10))  !Eq.28/2
         E_K = R*T/F*log(K_o/Y(13))             !Eq.28

         ! Fast Na+ current
         !print *,"======================================================"
         !print *, Cm
         !print *, g_Na
         !print *, Y(9)
         !print *, 'Y(7) =', Y(7)
         !print *, Y(8)
         !print *, Y(15)
         !print *, E_Na
         !print *,"======================================================"

         if ( abs(Y(7)) < epsilon(1.0_rp) ) then
            Y(7) = 0.0_rp
         end if
         i_Na = Cm*g_Na*Y(9)**(3.0_rp)*Y(7)*Y(8)*(Y(15)-E_Na) !Eq.29

         ! m
         if (Y(15) == -47.13_rp) then  !Eq.30
            alpha_m = 3.2_rp
         else
            alpha_m = 0.32_rp*(Y(15)+47.13_rp)/(1.0_rp-safe_exp(-0.1_rp*(Y(15)+47.13_rp)))
         end if
         beta_m = 0.08_rp*safe_exp(-Y(15)/11.0_rp)  !Eq.30
         m_inf = alpha_m/(alpha_m+beta_m)      !Eq.34
         tau_m = (alpha_m+beta_m)**(-1.0_rp)   !Eq.34

         ! h
         if (Y(15) .lt. -40.0_rp) then         !Eq,31
            alpha_h = 0.135_rp*safe_exp((Y(15)+80.0_rp)/(-6.8_rp))
         else
            alpha_h = 0.0_rp
         end if
         if (Y(15) .lt. -40_rp) then           !Eq.31    
            beta_h = 3.56_rp*safe_exp(0.079_rp*Y(15))+3.1e5_rp*safe_exp(0.35_rp*Y(15))
         else
            beta_h = 1.0_rp/(0.13_rp*(1.0_rp+safe_exp((Y(15)+10.66_rp)/(-11.1_rp))))
         end if          
         h_inf = alpha_h/(alpha_h+beta_h)      !Eq.34
         tau_h = 1.0_rp/(alpha_h+beta_h)       !Eq.34

         ! j 
         if (Y(15) .lt. -40.0_rp) then        !Eq.32
            alpha_j = (-1.2714e5_rp*safe_exp(0.2444_rp*Y(15))-3.47e-5_rp*safe_exp(-0.04391*Y(15))) & 
                      & *(Y(15)+37.78_rp)/(1.0_rp+safe_exp(0.311*(Y(15)+79.23_rp)))
         else
            alpha_j = 0.0_rp
         end if
         if (Y(15) .lt. -40.0_rp) then         !Eq.33
            beta_j = 0.1212_rp*safe_exp(-0.01052_rp*Y(15))/(1.0_rp+safe_exp(-0.1378_rp*(Y(15)+40.14_rp)))
         else
            beta_j = 0.3_rp*safe_exp(-2.535e-7_rp*Y(15))/(1.0_rp+safe_exp(-0.1_rp*(Y(15)+32.0_rp)))
         end if
         j_inf = alpha_j/(alpha_j+beta_j)      !Eq.34
         tau_j = 1.0_rp/(alpha_j+beta_j)       !Eq.34

         ! Time-independent K+ current
         i_K1 = Cm*g_K1*(Y(15)-E_K)/(1.0_rp+safe_exp(0.07_rp*(Y(15)+80.0_rp))) !Eq.35

         ! Transient Outward K+ current
         i_to = Cm*g_to*Y(18)**3.0_rp*Y(19)*(Y(15)-E_K)   !Eq.36

         ! oa
         alpha_oa = 0.65_rp*(safe_exp((Y(15)+10.0_rp)/(-8.5_rp))+safe_exp((Y(15)-30.0_rp)/(-59.0_rp)))**(-1.0_rp)  !Eq.37
         beta_oa = 0.65_rp*(2.5_rp+safe_exp((Y(15)+82.0_rp)/17.0_rp))**(-1.0_rp)                         !Eq.37
         tau_oa = ((alpha_oa+beta_oa)**(-1.0_rp))/K_Q10                                            !Eq.38
         oa_infty = (1.0_rp+safe_exp((Y(15)+20.47_rp)/(-17.54_rp)))**(-1.0_rp)                            !Eq.38 
  
         ! oi
         alpha_oi = (18.53_rp+1.0_rp*safe_exp((Y(15)+113.7_rp)/10.95_rp))**(-1.0_rp)                    !Eq.39
         beta_oi = (35.56_rp+1.0_rp*safe_exp((Y(15)+1.26_rp)/(-7.44)))**(-1.0_rp)                         !Eq.39
         tau_oi = ((alpha_oi+beta_oi)**(-1.0_rp))/K_Q10                                            !Eq.40
         oi_infty = (1.0_rp+safe_exp((Y(15)+43.1_rp)/5.3_rp))**(-1.0_rp)                                !Eq.40 

         ! Ultrarapid delayed rectifier K+ current
         g_Kur = (0.005_rp+0.05_rp/(1.0_rp+safe_exp((Y(15)-15.0_rp)/(-13.0_rp))))   !Eq.42
         i_Kur = Cm*g_Kur*Y(20)**3.0_rp*Y(21)*(Y(15)-E_K)                    !Eq.41
 
         ! ua
         alpha_ua = 0.65_rp*(safe_exp((Y(15)+10.0_rp)/(-8.5_rp))+safe_exp((Y(15)-30.0_rp)/(-59.0_rp)))**(-1.0_rp)  !Eq.43
         beta_ua = 0.65_rp*(2.5_rp+safe_exp((Y(15)+82.0_rp)/17.0_rp))**(-1.0_rp)                        !Eq.43
         tau_ua = ((alpha_ua+beta_ua)**(-1.0_rp))/K_Q10                                            !Eq.44
         ua_infty = (1.0_rp+safe_exp((Y(15)+30.03_rp)/(-9.6_rp)))**(-1.0_rp)                              !Eq.44 
      
         ! ui
         alpha_ui = (21.0_rp+1.0_rp*safe_exp((Y(15)-185.0_rp)/(-28.0_rp)))**(-1.0_rp)                     !Eq.45
         beta_ui = safe_exp((Y(15)-158.0_rp)/(16.0_rp))                                                 !Eq.45
         tau_ui = (alpha_ui+beta_ui)**(-1.0_rp)/K_Q10                                            !Eq.46
         ui_infty = (1.0_rp+safe_exp((Y(15)-99.45_rp)/27.48_rp))**(-1.0_rp)                             !Eq.46 

         ! Rapid delayed outward rectifier K+ current
         i_Kr = Cm*g_Kr*Y(16)*(Y(15)-E_K)/(1.0_rp+safe_exp((Y(15)+15.0_rp)/22.4))   !Eq.47
 
         ! xr
         if (abs(Y(15)+14.1_rp) .lt. 1.0e-10_rp) then   !Eq.48
            alpha_xr = 0.0015_rp
         else
            alpha_xr = 0.0003_rp*(Y(15)+14.1_rp)/(1.0_rp-safe_exp((Y(15)+14.1_rp)/(-5.0_rp)))
         end if
         if (abs(Y(15)-3.3328_rp) .lt. 1.0e-10_rp) then !Eq.48
            beta_xr = 3.7836118e-4_rp
         else 
            beta_xr = 0.000073898_rp*(Y(15)-3.3328_rp)/(safe_exp((Y(15)-3.3328_rp)/5.1237_rp)-1.0_rp)
         end if
         tau_xr = (alpha_xr+beta_xr)**(-1.0_rp)                       !Eq.49
         xr_infty = (1.0_rp+safe_exp((Y(15)+14.1_rp)/(-6.5_rp)))**(-1-0_rp)  !Eq.49


         ! Slow dealyed outward rectifier K+ current
         i_Ks = Cm*g_Ks*Y(17)**2.0_rp*(Y(15)-E_K)    !Eq.50
 
         ! xs
         if (abs(Y(15)-19.9_rp) .lt. 1.0e-10_rp) then   !Eq.51
            alpha_xs = 0.00068_rp
         else
            alpha_xs = 0.00004_rp*(Y(15)-19.9_rp)/(1.0_rp-safe_exp((Y(15)-19.9_rp)/(-17.0_rp)))
         end if
         if (abs(Y(15)-19.9_rp) .lt. 1.0e-10_rp) then   !Eq.51
            beta_xs = 0.000315_rp
         else
            beta_xs = 0.000035_rp*(Y(15)-19.9_rp)/(safe_exp((Y(15)-19.9_rp)/9.0_rp)-1.0_rp)
         end if
         tau_xs = 0.5_rp*(alpha_xs+beta_xs)**(-1.0_rp)                 !Eq.52
         xs_infty = (1.0_rp+safe_exp((Y(15)-19.9_rp)/(-12.7_rp)))**(-0.5_rp)  !Eq.52


         ! L-Type Ca2+ current
         i_Ca_L = Cm*g_Ca_L*Y(4)*Y(6)*Y(5)*(Y(15)-65.0_rp)   !Eq.53

         ! d
         if (abs(Y(15)+10.0_rp) .lt. 1.0e-10_rp) then  !Eq.54
            tau_d = 4.579_rp/(1.0_rp+safe_exp((Y(15)+10.0_rp)/(-6.24_rp)))
         else 
            tau_d = (1.0_rp-safe_exp((Y(15)+10.0_rp)/(-6.24_rp)))/(0.035_rp*(Y(15)+10.0_rp)*(1.0_rp+safe_exp((Y(15)+10.0_rp)/(-6.24_rp))))
         end if
         d_infty = (1.0_rp+safe_exp((Y(15)+10.0_rp)/(-8.0_rp)))**(-1.0_rp)  !Eq.54

         ! f
         tau_f = 9.0_rp*(0.0197_rp*safe_exp(-0.0337_rp**2.0_rp*(Y(15)+10.0_rp)**2.0_rp)+0.02_rp)**(-1.0_rp)  !Eq.55
         f_infty = (1.0_rp+safe_exp((Y(15)+28.0_rp)/6.9_rp))**(-1.0_rp)

         ! f_Ca
         tau_f_Ca = 2.0_rp  !Eq.56
         f_Ca_infty = (1.0_rp+Y(10)/0.00035_rp)**(-1.0_rp)  !Eq.56


         ! Na+-K- pump Current
         sigma = 1.0_rp/7.0_rp*(safe_exp(Na_o/67.3_rp)-1.0_rp)   !Eq.59
         f_NaK = (1.0_rp+0.1245_rp*safe_exp(-0.1_rp*F*Y(15)/(R*T)) + 0.0365_rp*sigma*safe_exp(-F*Y(15)/(R*T)))**(-1.0_rp)  !Eq.58
         i_NaK = Cm*i_NaK_max*f_NaK*(1.0_rp/(1.0_rp+(Km_Na_i/Y(14))**1.5_rp))*(K_o/(K_o+Km_K_o))   !Eq.57

         ! Na+/Ca2+ exchanger current
         i_NaCa = Cm*I_NaCa_max*(safe_exp(gama*F*Y(15)/(R*T))*Y(14)**3.0_rp*Ca_o-safe_exp((gama-1.0_rp)*F*Y(15)/(R*T))*Na_o**3.0_rp*Y(10)) &
                  & /((K_mNa**3.0_rp+Na_o**3.0_rp)*(K_mCa+Ca_o)*(1.0_rp+K_sat*safe_exp((gama-1.0_rp)*Y(15)*F/(R*T))))   !Eq.60 

         ! Background currents
         i_B_Na = Cm*g_B_Na*(Y(15)-E_Na)  !Eq.62
         i_B_Ca = Cm*g_B_Ca*(Y(15)-E_Ca)  !Eq.61
         i_B_K = Cm*g_B_K*(Y(15)-E_K)     !does not exist

         ! Ca2+ pump current
         i_CaP = Cm*i_CaP_max*Y(10)/(0.0005_rp+Y(10))

         ! Ca2+ release current from JSR
         i_rel = K_rel*Y(1)**2.0_rp*Y(2)*Y(3)*(Y(11)-Y(10))   !Eq.64
         Fn = 1e-12_rp*V_rel*i_rel - ((5e-13_rp)/F)*(0.5_rp*i_Ca_L-0.2_rp*i_NaCa)  !Eq.68 

         ! u 
         tau_u = 8.0_rp   !Eq.65
         if ( (1.0_rp/(1.0_rp+safe_exp(-(Fn-3.4175e-13_rp)/13.67e-16_rp))) .lt. 1.0e-25_rp) then  !Eq.65
            u_infty = 0.0_rp
         else 
            u_infty = 1.0_rp/(1.0_rp+safe_exp(-(Fn-3.4175e-13_rp)/13.67e-16_rp))
         end if

         ! v
         tau_v = 1.91_rp+2.09_rp*(1.0_rp+safe_exp(-(Fn-3.4175e-13_rp)/13.67e-16_rp))**(-1.0_rp)   !Eq.66
         v_infty = 1.0_rp-(1.0_rp+safe_exp(-(Fn-6.835e-14_rp)/13.67e-16_rp))**(-1.0_rp)            !Eq.66

         ! w
         if (abs(Y(15)-7.9_rp) .lt. 1.0e-10_rp) then  !Eq.67
            tau_w = 6.0_rp*0.2_rp/1.3_rp
         else 
            tau_w = 6.0_rp*(1.0_rp-safe_exp(-(Y(15)-7.9_rp)/5.0_rp)) & 
                    & /((1.0_rp+0.3_rp*safe_exp(-(Y(15)-7.9_rp)/5.0_rp))*1.0_rp*(Y(15)-7.9_rp))
         end if
         w_infty = 1.0_rp-(1.0_rp+safe_exp(-(Y(15)-40.0_rp)/17.0_rp))**(-1.0_rp)   !Eq.67

         ! Transfer current from NSR to JSR
         i_tr = (Y(12)-Y(11))/tau_tr   !Eq.69

         ! Ca2+ leak current by the NSR
         i_up_leak = I_up_max*Y(12)/Ca_up_max   !Eq.72

         ! Ca2- uptake current by the NSAR
         i_up = I_up_max/(1.0_rp+K_up/Y(10))   !Eq.71

         if (flag_isac) then

            a = GG_sac
            b=(1+K_sac*exp(-(alpha_sac*(lambda-1))))
            g_sac=a/b         !conductancia del canal     

            A_Na= P_na*g_sac*Z_na**2.0_rp*F**2.0_rp*Y(15)/(R*T)
            B_Na= (Y(14)-Na_o*exp(-Z_na*F*Y(15)/(R*T)))
            C_Na= (1-exp(-Z_na*F*Y(15)/(R*T)))
            I_SAC_Na= A_Na*B_Na/C_Na

            A_K= P_k*g_sac*Z_k**2.0_rp*F**2.0_rp*Y(15)/(R*T)
            B_K= (Y(13)-K_o*exp(-Z_k*F*Y(15)/(R*T)))
            C_K= (1-exp(-Z_k*F*Y(15)/(R*T)))
            I_SAC_K=A_K*B_K/C_K

            A_Ca= P_ca*g_sac*Z_ca**2.0_rp*F**2.0_rp*Y(15)/(R*T)
            B_Ca= (Y(10)-Ca_o*exp(-Z_ca*F*Y(15)/(R*T)))
            C_Ca= (1-exp(-Z_ca*F*Y(15)/(R*T)))
            I_SAC_Ca=A_Ca*B_Ca/C_Ca

            I_SAC =I_SAC_Na + I_SAC_K + I_SAC_Ca

         else 
            I_SAC_Na = 0.0_rp
            I_SAC_Ca = 0.0_rp
            I_SAC_K = 0.0_rp

            I_SAC =I_SAC_Na + I_SAC_K + I_SAC_Ca
     
         end if

            ! Total ion current: update ionic currents
            Iion=i_Na+i_K1+i_to+i_Kur+i_Kr+i_Ks+i_B_Na+i_B_Ca+i_NaK+i_CaP+i_NaCa+i_Ca_L+I_SAC;

         ! Ca2+ buffers
         Ca_CMDN = CMDN_max*Y(10)/(Y(10)+Km_CMDN)  !Eq.73
         Ca_TRPN = TRPN_max*Y(10)/(Y(10)+Km_TRPN)  !Eq.74
         Ca_CSQN = CSQN_max*Y(11)/(Y(11)+Km_CSQN)  !Eq.75


         !--------------------------------------------------------------------------------
         ! Sitmulus current
         !--------------------------------------------------------------------------------
         
         ! Update membrance currents 
         ti_exm(1) = i_Na
         ti_exm(2) = i_K1
         ti_exm(3) = i_to
         ti_exm(4) = i_Kur
         ti_exm(5) = i_Kr
         ti_exm(6) = i_Ks
         ti_exm(7) = i_Ca_L
         ti_exm(8) = i_NaK
         ti_exm(9) = i_NaCa
         ti_exm(10) = i_CaP
         ti_exm(11) = i_B_Na
         ti_exm(12) = i_B_Ca
         ti_exm(13) = i_B_K
         ti_exm(14) = i_rel
         ti_exm(15) = i_tr     
         ti_exm(16) = i_up_leak      
         ti_exm(17) = i_up  
         ti_exm(18) = i_st 


         !--------------------------------------------------------------------------------
         ! Differential equations
         !--------------------------------------------------------------------------------
        
         dY(1) = (u_infty-Y(1))/tau_u
         dY(2) = (v_infty-Y(2))/tau_v
         dY(3) = (w_infty-Y(3))/tau_w
         dY(4) = (d_infty-Y(4))/tau_d
         dY(5) = (f_Ca_infty-Y(5))/tau_f_Ca 
         dY(6) = (f_infty-Y(6))/tau_f
         dY(7) = (h_inf-Y(7))/tau_h
         dY(8) = (j_inf-Y(8))/tau_j
         dY(9) = (m_inf-Y(9))/tau_m


         if (flag_land) then
           vaux1 = (Km_CMDN+Y(10)) * (Km_CMDN+Y(10))
           vaux3 = 1.0_rp/(1.0_rp + (CMDN_max*Km_CMDN/vaux1))
           vaux4 = (Y(10)/(Ca50_ref*Ca50_fact))**n_TRPN * (1.0_rp-Ca_TRPN) 
           dCaTRPN = Km_TRPN *( vaux4 - Ca_TRPN)  
           rhsx1 = 2.0_rp*i_NaCa - i_CaP - i_Ca_L - i_B_Ca
           rhsx2 = V_up*(i_up_leak - i_up) + i_rel*V_rel
           rhsx = vaux3*(rhsx1/(2.0_rp*F*V_i) + rhsx2/V_i- dCaTRPN)
           dY(10) = rhsx
         else
           B1 = (2.0_rp*i_NaCa-(i_CaP+i_Ca_L+i_B_Ca+I_SAC_Ca))/(2.0_rp*V_i*F)+(V_up*(i_up_leak-i_up)+i_rel*V_rel)/V_i   !Eq.24
           B2 = 1.0_rp+TRPN_max*Km_TRPN/(Y(10)+Km_TRPN)**2.0_rp+CMDN_max*Km_CMDN/(Y(10)+Km_CMDN)**2.0_rp   !Eq.25       
           dY(10) = B1/B2    !Eq.23   !EVA: ECUACION DEL CA_i
         end if
         dY(11) = (i_tr-i_rel)*(1.0_rp+CSQN_max*Km_CSQN/(Y(11)+Km_CSQN)**2.0_rp)**(-1.0_rp)   !Eq.27
         dY(12) = i_up-(i_up_leak+i_tr*V_rel/V_up)   !Eq.26
         dY(13) = (2.0_rp*i_NaK-(i_st+i_K1+i_to+i_Kur+i_Kr+i_Ks+i_B_K+I_SAC_K))/(V_i*F)   !Eq.22
         dY(14) = (-3.0_rp*i_NaK-(3.0_rp*i_NaCa+i_B_Na+i_Na+I_SAC_Na))/(V_i*F)   !Eq.21
         dY(15) = -(Iion+i_st)/Cm    !Cm=xmccmmate_exm
         dY(16) = (xr_infty-Y(16))/tau_xr
         dY(17) = (xs_infty-Y(17))/tau_xs
         dY(18) = (oa_infty-Y(18))/tau_oa
         dY(19) = (oi_infty-Y(19))/tau_oi
         dY(20) = (ua_infty-Y(20))/tau_ua
         dY(21) = (ui_infty-Y(21))/tau_ui


         ! Solve equations: voltage integration
         ! y^n+1 = y^n + dY/dt

         ! Euler forward
         do i=1,21
            Y(i) = Y(i) + dt*dY(i)
         end do 

         !--------------------------------------------------------------------------------
         ! Update variables
         !--------------------------------------------------------------------------------

         ! Intracellular ion concentrations
         tconc_exm(1,2) = Y(10)  !Ca_i 
         tconc_exm(2,2) = Y(11)  !Ca_rel      
         tconc_exm(3,2) = Y(12)  !Ca_up      
         tconc_exm(4,2) = Y(13)  !K_i      
         tconc_exm(5,2) = Y(14)  !Na_i

         ! Gate variables     
         tgates_exm(1,2) = Y(1)  !u
         tgates_exm(2,2) = Y(2)  !v
         tgates_exm(3,2) = Y(3)  !w
         tgates_exm(4,2) = Y(4)  !d
         tgates_exm(5,2) = Y(5)  !f_Ca
         tgates_exm(6,2) = Y(6)  !f
         tgates_exm(7,2) = Y(7)  !h
         tgates_exm(8,2) = Y(8)  !j
         tgates_exm(9,2) = Y(9)  !m
         tgates_exm(10,2) = Y(16)  !xr
         tgates_exm(11,2) = Y(17)  !xs
         tgates_exm(12,2) = Y(18)  !oa
         tgates_exm(13,2) = Y(19)  !oi
         tgates_exm(14,2) = Y(20)  !ua
         tgates_exm(15,2) = Y(21)  !ui

         U(2) = Y(15) !V

        !----------------------------------------------------------------------------------
        ! Outputs
        !-----------------------------------------------------------------------------------
        Iion = 0.0_rp
        do i=1,(nicel_exm-1_ip)
          Iion = + ti_exm(i)
        end do

        ! Evaluate updated membrane potential
        dU = -Iion
        U_n = U(2) + dU*dt

contains
#include "exm_safe_exp.f90.inc"
end subroutine exm_courtemanche_model
end module mod_exm_courtemanche_model
