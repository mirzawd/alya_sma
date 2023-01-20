!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!########################### CONSTANTS #####################
module mod_exm_hipscatriaconst
!! Constants
  use def_kintyp, only : ip,rp

  implicit none

!ionic concentrations
real(rp), parameter :: &
Na_o = 151.0_rp,  K_o = 5.4_rp,  Ca_o = 1.8_rp,  K_i = 150.0_rp

!!Atrial
!Cell size
real(rp), parameter :: &
C_m = 78.6671e-12_rp,    V_c = 7012.0_rp,    V_sr = 465.20_rp, &
!Maximum conductances and currents
g_Na = 6.646185_rp*(10.0_rp**3.0_rp),    g_to = 59.8077_rp,    g_K1 = 19.1925_rp,    P_NaK = 1.4731392_rp,    K_NaCa = 2450.0_rp,    Vmax_up = 0.22_rp


!Maximum conductances and currents
real(rp), parameter :: &
g_CaL = 8.635702e-5_rp,  g_Kr = 29.8667_rp,  g_Ks = 2.041_rp,  g_f = 30.10312_rp,  a_rel = 16.464_rp, &
b_rel = 0.25_rp,  c_rel = 8.232_rp,  V_leak = 4.4444e-4_rp,  g_pCa = 0.4125_rp,  g_bNa = 0.9_rp,  g_bCa = 0.69264_rp

!Other constants
real(rp), parameter :: &
Buf_c = 0.25_rp,  Buf_sr = 10.0_rp,  K_Buf_c = 0.001_rp,  K_Buf_sr = 0.3_rp,  K_up = 0.00025_rp, &
K_pCa = 0.0005_rp,  F = 96485.3415_rp,  Rcons = 8.314472_rp,  Temp = 310.0_rp,  L_0 = 0.025_rp,  P_kna = 0.03_rp, &
K_sat = 0.1_rp,  Km_Ca = 1.38_rp,  Km_Na_i = 87.5_rp,  Alpha = 2.8571432_rp,  Gamma = 0.35_rp, &
K_mNa = 40.0_rp,  K_mk = 1.0_rp
end module mod_exm_hipscatriaconst

!#####################################################################################
