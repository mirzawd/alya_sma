!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!!#############################
    !voltage: elmag
    
    !Ca_i:vconc(1,:)
    !Ca_SR:vconc(2,:)
    !Na_i:vconc(3,:)
    
    !h:vauxi_exm(1,:)  ;  j:vauxi_exm(2,:)  ;  !m:vauxi_exm(3,:)
    !d:vauxi_exm(4,:)  ;  f_Ca:vauxi_exm(5,:)  ;  !f1:vauxi_exm(6,:)
    !f2:vauxi_exm(7,:)  ;  r:vauxi_exm(8,:)  ;  !q:vauxi_exm(9,:)
    !xr1:vauxi_exm(10,:)  ;  xr2:vauxi_exm(11,:)  ;  !Xs:vauxi_exm(12,:)
    !Xf:vauxi_exm(13,:)  ;  g:vauxi_exm(14,:)
    
    !I_Na:vicel_exm(1,:)  ;  I_Cal:vicel_exm(2,:)  ;  !I_to:vicel_exm(3,:)
    !I_kr:vicel_exm(4,:)  ;  I_ks:vicel_exm(5,:)  ;  !I_k1:vicel_exm(6,:)
    !I_f:vicel_exm(7,:)  ;  I_NaK:vicel_exm(8,:)  ;  !I_NaCa:vicel_exm(9,:)
    !I_rel:vicel_exm(10,:)  ;  I_up:vicel_exm(11,:)  ;  !I_leak:vicel_exm(12,:)
    !I_pCa:vicel_exm(13,:)  ;  I_bNa:vicel_exm(14,:)  ;  I_bCa:vicel_exm(15,:)
!!#############################
!! main code
module mod_exm_scvent_ionicurrents
    implicit none
    private

    public :: exm_scvent_ionicurrents

    contains
    
    

#include "exm_safe_exp.f90.inc"

subroutine exm_scvent_ionicurrents(ipoin,xioni,dioni)

  use      def_parame
  use      def_master
  use      def_elmtyp
  use      def_domain
  use      def_exmedi
  use      mod_exm_hipscventriconst


  ! definition of variables
  implicit none
  integer(ip), intent(in) :: ipoin !< node
  real(rp), intent(out) :: xioni   !< current
  real(rp), intent(out)   :: dioni !< current derivative

  real(rp) :: tt, h, j, m, d, f_Ca, f1, f2, r, q, xr1, xr2, xs, xf, g, Ca_i, Ca_SR, Na_i, V
    
  real(rp) :: E_Na, E_K, E_Ks, E_f, E_Ca, Alpha_K1, beta_K1, x_K1_inf
  !! RK constants declaration
  real(rp) :: K1h, K1j, K1m, K1d, K1f_Ca, K1f1, K1f2, &
  K1r, K1q, K1xr1, K1xr2, K1Xs, K1Xf, K1g, K1Ca_i, K1Ca_SR, K1Na_i, K1V

  real(rp) :: K2h, K2j, K2m, K2d, K2f_Ca, K2f1, K2f2, &
  K2r, K2q, K2xr1, K2xr2, K2Xs, K2Xf, K2g, K2Ca_i, K2Ca_SR, K2Na_i, K2V

  real(rp) :: K3h, K3j, K3m, K3d, K3f_Ca, K3f1, K3f2, &
  K3r, K3q, K3xr1, K3xr2, K3Xs, K3Xf, K3g, K3Ca_i, K3Ca_SR, K3Na_i, K3V

  real(rp) :: K4h, K4j, K4m, K4d, K4f_Ca, K4f1, K4f2, &
  K4r, K4q, K4xr1, K4xr2, K4Xs, K4Xf, K4g, K4Ca_i, K4Ca_SR, K4Na_i, K4V
  !! !! defining functions and Solving the system of ODE using Runge-kutta of 4th order

  real(rp) :: df1 = -1.0_rp

  V = elmag(ipoin,ITER_K)/1000.0_rp
  tt = cutim
  h = vauxi_exm(1,ipoin,2)
  j = vauxi_exm(2,ipoin,2)
  m = vauxi_exm(3,ipoin,2)
  d = vauxi_exm(4,ipoin,2)
  f_Ca = vauxi_exm(5,ipoin,2)
  f1 = vauxi_exm(6,ipoin,2)
  f2 = vauxi_exm(7,ipoin,2)
  r = vauxi_exm(8,ipoin,2)
  Q = vauxi_exm(9,ipoin,2)
  xr1 = vauxi_exm(10,ipoin,2)
  xr2 = vauxi_exm(11,ipoin,2)
  xs = vauxi_exm(12,ipoin,2)
  xf = vauxi_exm(13,ipoin,2)
  g = vauxi_exm(14,ipoin,2)
  Ca_i = vconc(1,ipoin,2)
  Ca_SR = vconc(2,ipoin,2)
  Na_i = vconc(3,ipoin,2)
  
  E_Na = (Rcons*Temp/F)*log(Na_o/vconc(3,ipoin,2))
  E_K = (Rcons*Temp/F)*log(K_o/K_i)
  E_Ks = (Rcons*Temp/F)*log((K_o+P_kna*Na_o)/(K_i+P_kna*vconc(3,ipoin,2)))
  E_f = -0.017_rp
  E_Ca = (0.5_rp*Rcons*Temp/F)*log(Ca_o/vconc(1,ipoin,2))
  Alpha_K1 = 3.91_rp/(1.0_rp+safe_exp(0.5942_rp*(V*1000.0_rp-E_K*1000.0_rp-200.0_rp)))
  beta_K1 = (-1.509_rp*safe_exp(0.0002_rp*(V*1000.0_rp-E_K*1000.0_rp+100.0_rp))+safe_exp(0.5886_rp*(V*1000.0_rp-E_K*1000.0_rp-10.0_rp)))/ &
            (1.0_rp+safe_exp(0.4547_rp*(V*1000.0_rp-E_K*1000.0_rp)))
  x_K1_inf = Alpha_K1/(Alpha_K1+beta_K1)
    
  vicel_exm(1,ipoin) = g_Na*(vauxi_exm(3,ipoin,2)**3.0_rp)*vauxi_exm(1,ipoin,2)*vauxi_exm(2,ipoin,2)*(V-E_Na)
  vicel_exm(2,ipoin) = (g_CaL*4.0_rp*V*F**2.0_rp)*(vconc(1,ipoin,2)*safe_exp(2.0_rp*V*F/(Rcons*Temp))-0.341_rp*Ca_o)*vauxi_exm(4,ipoin,2)*vauxi_exm(6,ipoin,2)*vauxi_exm(7,ipoin,2)*f_Ca/ & 
                   (Rcons*Temp*(safe_exp(2.0_rp*V*F/(Rcons*Temp))-1.0_rp))
  vicel_exm(3,ipoin) = g_to*vauxi_exm(8,ipoin,2)*vauxi_exm(9,ipoin,2)*(V-E_K)
  vicel_exm(4,ipoin) = g_Kr*sqrt(K_o/5.4_rp)*vauxi_exm(10,ipoin,2)*vauxi_exm(11,ipoin,2)*(V-E_K)
  vicel_exm(5,ipoin) = g_Ks*vauxi_exm(12,ipoin,2)**2.0_rp*(1.0_rp+0.6_rp/(1.0_rp+(3.8_rp*1e-5_rp/vconc(1,ipoin,2))**1.4_rp))*(V-E_Ks)
  vicel_exm(6,ipoin) = g_K1*x_K1_inf*sqrt(K_o/5.4_rp)*(V-E_K)
  vicel_exm(7,ipoin) = g_f*vauxi_exm(13,ipoin,2)*(V-E_f)
  vicel_exm(8,ipoin) = ((P_NaK*K_o*vconc(3,ipoin,2)/(K_o+K_mk))/(vconc(3,ipoin,2)+K_mNa))/ &
                   (1.0_rp+0.1245_rp*safe_exp(-0.1_rp*V*F/(Rcons*Temp))+0.353_rp*safe_exp(-V*F/(Rcons*Temp)))
  vicel_exm(9,ipoin) = K_NaCa*(safe_exp(Gamma*V*F/(Rcons*Temp))*vconc(3,ipoin,2)**3.0_rp*Ca_o-safe_exp((Gamma-1.0_rp)*V*F/(Rcons*Temp)) &
                   *Na_o**3.0_rp*vconc(1,ipoin,2)*Alpha)/((Km_Na_i**3.0_rp+Na_o**3.0_rp)*(Km_Ca+Ca_o)* &
                   (1.0_rp+K_sat*safe_exp((Gamma-1.0_rp)*V*F/(Rcons*Temp))))
  vicel_exm(10,ipoin) = ((a_rel*(vconc(2,ipoin,2)**2.0_rp)/(b_rel**2.0_rp+vconc(2,ipoin,2)**2.0_rp))+c_rel)*vauxi_exm(4,ipoin,2)*vauxi_exm(14,ipoin,2)*(0.0556_rp)
  vicel_exm(11,ipoin) = Vmax_up/(1.0_rp+K_up**2.0_rp/vconc(1,ipoin,2)**2.0_rp)
  vicel_exm(12,ipoin) = V_leak*(vconc(2,ipoin,2)-vconc(1,ipoin,2))
  vicel_exm(13,ipoin) = g_pCa*vconc(1,ipoin,2)/(vconc(1,ipoin,2)+K_pCa)
  vicel_exm(14,ipoin) = g_bNa*(V-E_Na)
  vicel_exm(15,ipoin) = g_bCa*(V-E_Ca)
  
  !! Applying Runge-Kutta
  !open(unit=1, file='voltage.dat')

    !! K1
    K1h = dtime*func_hV(tt,V,h)
    K1j = dtime*func_jV(tt,V,j)
    K1m = dtime*func_mV(tt,V,m)
    K1d = dtime*func_dV(tt,V,d)
    K1f_Ca = dtime*func_fCaV(tt,f_Ca,Ca_i,V)
    K1f1 = dtime*func_f1V(tt,V,f1,Ca_i,df1)
    K1f2 = dtime*func_f2V(tt,V,f2)
    K1r = dtime*func_rV(tt,V,r)
    K1q = dtime*func_qV(tt,V,q)
    K1xr1 = dtime*func_xr1V(tt,V,xr1)
    K1xr2 = dtime*func_xr2V(tt,V,xr2)
    K1Xs = dtime*func_XsV(tt,V,xs)
    K1Xf = dtime*func_XfV(tt,V,xf)
    K1g = dtime*func_gV(tt,Ca_i,g,V)
    K1Ca_i = dtime*func_CaiV(tt,Ca_i,Ca_SR,d,g,V,Na_i,f1,f2,f_Ca)
    K1Ca_SR = dtime*func_CaSRV(tt,Ca_i,Ca_SR,d,g)
    K1Na_i = dtime*func_NaiV(tt,m,h,j,V,Na_i,Ca_i)
    K1V = dtime*func_VV(tt,V,r,q,xr1,xr2,xs,Ca_i,d,f1,f2,f_Ca,Na_i,m,h,j,xf)

    !! K2
    K2h =  dtime*func_hV(tt+dtime/2.0_rp,  V+K1V/2.0_rp,  h+K1h/2.0_rp)
    K2j =  dtime*func_jV(tt+dtime/2.0_rp,  V+K1V/2.0_rp,  j+K1j/2.0_rp)
    K2m =     dtime*func_mV(tt+dtime/2.0_rp,  V+K1V/2.0_rp,  m+K1m/2.0_rp)
    K2d =     dtime*func_dV(tt+dtime/2.0_rp,  V+K1V/2.0_rp,  d+K1d/2.0_rp)
    K2f_Ca =  dtime*func_fCaV(tt+dtime/2.0_rp,  f_Ca+K1f_Ca/2.0_rp,  Ca_i+K1Ca_i/2.0_rp ,V+K1V/2.0_rp)
    K2f1 = dtime*func_f1V(tt+dtime/2.0_rp,  V+K1V/2.0_rp,  f1+K1f1/2.0_rp  ,Ca_i+K1f1/2.0_rp ,df1)
    K2f2 =    dtime*func_f2V(tt+dtime/2.0_rp,  V+K1V/2.0_rp,  f2+K1f2/2.0_rp)
    K2r =     dtime*func_rV(tt+dtime/2.0_rp,  V+K1V/2.0_rp,  r+K1r/2.0_rp)
    K2q =     dtime*func_qV(tt+dtime/2.0_rp,  V+K1V/2.0_rp,  q+K1q/2.0_rp)
    K2xr1 =   dtime*func_xr1V(tt+dtime/2.0_rp,V+K1V/2.0_rp,  xr1+K1xr1/2.0_rp)
    K2xr2 =   dtime*func_xr2V(tt+dtime/2.0_rp,  V+K1V/2.0_rp,  xr2+K1xr2/2.0_rp)
    K2Xs =    dtime*func_XsV(tt+dtime/2.0_rp,  V+K1V/2.0_rp,  xs+K1Xs/2.0_rp)
    K2Xf =    dtime*func_XfV(tt+dtime/2.0_rp,  V+K1V/2.0_rp,  xf+K1Xf/2.0_rp)
    K2g =     dtime*func_gV(tt+dtime/2.0_rp,  Ca_i+K1Ca_i/2.0_rp,  g+K1g/2.0_rp ,V+K1V/2.0_rp)
    K2Ca_i =  dtime*func_CaiV(tt+dtime/2.0_rp,  Ca_i+K1Ca_i/2.0_rp,  Ca_SR+K1Ca_SR/2.0_rp,  d+K1d/2.0_rp,  g+K1g/2.0_rp,  V+K1V/2.0_rp, Na_i+K1Na_i/2.0_rp,  f1+K1f1/2.0_rp  ,f2+K1f2/2.0_rp  ,f_Ca+K1f_Ca/2.0_rp)
    K2Ca_SR = dtime*func_CaSRV(tt+dtime/2.0_rp,  Ca_i+K1Ca_i/2.0_rp,  Ca_SR+K1Ca_SR/2.0_rp,  d+K1d/2.0_rp,  g+K1g/2.0_rp)
    K2Na_i =  dtime*func_NaiV(tt+dtime/2.0_rp,  m+K1m/2.0_rp,  h+K1h/2.0_rp,  j+K1j/2.0_rp,  V+K1V/2.0_rp,  Na_i+K1Na_i/2.0_rp,  Ca_i+K1Ca_i/2.0_rp)
    K2V =     dtime*func_VV(tt+dtime/2.0_rp,  V+K1V/2.0_rp,  r+K1r/2.0_rp,  q+K1q/2.0_rp,  &
            xr1+K1xr1/2.0_rp,  xr2+K1xr2/2.0_rp,  xs+K1Xs/2.0_rp, &
              Ca_i+K1Ca_i/2.0_rp,  d+K1d/2.0_rp,  f1+K1f1/2.0_rp,  f2+K1f2/2.0_rp,  &
              f_Ca+K1f_Ca/2.0_rp,  Na_i+K1Na_i/2.0_rp,  m+K1m/2.0_rp,  h+K1h/2.0_rp,  j+K1j/2.0_rp  ,xf+K1Xf/2.0_rp)

    !! K3
    K3h = dtime*func_hV(tt+dtime/2.0_rp,  V+K2V/2.0_rp,  h+K2h/2.0_rp)
    K3j = dtime*func_jV(tt+dtime/2.0_rp,  V+K2V/2.0_rp,  j+K2j/2.0_rp)
    K3m =     dtime*func_mV(tt+dtime/2.0_rp,  V+K2V/2.0_rp,  m+K2m/2.0_rp)
    K3d =     dtime*func_dV(tt+dtime/2.0_rp,  V+K2V/2.0_rp,  d+K2d/2.0_rp)
    K3f_Ca =  dtime*func_fCaV(tt+dtime/2.0_rp,  f_Ca+K2f_Ca/2.0_rp,  Ca_i+K2Ca_i/2.0_rp ,V+K2V/2.0_rp)
    K3f1 = dtime*func_f1V(tt+dtime/2.0_rp,  V+K2V/2.0_rp,  f1+K2f1/2.0_rp, Ca_i+K2f1/2.0_rp,df1)
    K3f2 =    dtime*func_f2V(tt+dtime/2.0_rp,  V+K2V/2.0_rp,  f2+K2f2/2.0_rp)
    K3r =     dtime*func_rV(tt+dtime/2.0_rp,  V+K2V/2.0_rp,  r+K2r/2.0_rp)
    K3q =     dtime*func_qV(tt+dtime/2.0_rp,  V+K2V/2.0_rp,  q+K2q/2.0_rp)
    K3xr1 =   dtime*func_xr1V(tt+dtime/2.0_rp,V+K2V/2.0_rp,  xr1+K2xr1/2.0_rp)
    K3xr2 =  dtime* func_xr2V(tt+dtime/2.0_rp,  V+K2V/2.0_rp,  xr2+K2xr2/2.0_rp)
    K3Xs =    dtime*func_XsV(tt+dtime/2.0_rp,  V+K2V/2.0_rp,  xs+K2Xs/2.0_rp)
    K3Xf =    dtime*func_XfV(tt+dtime/2.0_rp,  V+K2V/2.0_rp,  xf+K2Xf/2.0_rp)
    K3g =     dtime*func_gV(tt+dtime/2.0_rp,  Ca_i+K2Ca_i/2.0_rp,  g+K2g/2.0_rp ,V+K2V/2.0_rp)
    K3Ca_i =  dtime*func_CaiV(tt+dtime/2.0_rp,  Ca_i+K2Ca_i/2.0_rp,  Ca_SR+K2Ca_SR/2.0_rp,  d+K2d/2.0_rp,  g+K2g/2.0_rp,  V+K2V/2.0_rp,  Na_i+K2Na_i/2.0_rp,  f1+K2f1/2.0_rp  ,f2+K2f2/2.0_rp  ,f_Ca+K2f_Ca/2.0_rp)
    K3Ca_SR = dtime*func_CaSRV(tt+dtime/2.0_rp,  Ca_i+K2Ca_i/2.0_rp,  Ca_SR+K2Ca_SR/2.0_rp,  d+K2d/2.0_rp,  g+K2g/2.0_rp)
    K3Na_i =  dtime*func_NaiV(tt+dtime/2.0_rp,  m+K2m/2.0_rp,  h+K2h/2.0_rp,  j+K2j/2.0_rp,  V+K2V/2.0_rp,  Na_i+K2Na_i/2.0_rp,  Ca_i+K2Ca_i/2.0_rp)
    K3V =     dtime*func_VV(tt+dtime/2.0_rp,  V+K2V/2.0_rp,  r+K2r/2.0_rp,  q+K2q/2.0_rp,  xr1+K2xr1/2.0_rp,  &
            xr2+K2xr2/2.0_rp,  xs+K2Xs/2.0_rp,  Ca_i+K2Ca_i/2.0_rp,  d+K2d/2.0_rp,  &
            f1+K2f1/2.0_rp,  f2+K2f2/2.0_rp,  f_Ca+K2f_Ca/2.0_rp,  Na_i+K2Na_i/2.0_rp,  &
            m+K2m/2.0_rp,  h+K2h/2.0_rp,  j+K2j/2.0_rp  ,xf+K2Xf/2.0_rp)

    !! K4
    K4h = dtime*func_hV(tt+dtime,  V+K3V,  h+K3h)
    K4j = dtime*func_jV(tt+dtime,  V+K3V,  j+K3j)
    K4m = dtime*func_mV(tt+dtime,  V+K3V,  m+K3m)
    K4d = dtime*func_dV(tt+dtime,  V+K3V,  d+K3d)
    K4f_Ca =  dtime*func_fCaV(tt+dtime,  f_Ca+K3f_Ca,  Ca_i+K3Ca_i ,V+K3V)
    K4f1 = dtime*func_f1V(tt+dtime,  V+K3V,  f1+K3f1, Ca_i+K3f1 ,df1)
    K4f2 =    dtime*func_f2V(tt+dtime,  V+K3V,  f2+K3f2)
    K4r =     dtime*func_rV(tt+dtime,  V+K3V,  r+K3r)
    K4q =     dtime*func_qV(tt+dtime,  V+K3V,  q+K3q)
    K4xr1 =   dtime*func_xr1V(tt+dtime, V+K3V,  xr1+K3xr1)
    K4xr2 =   dtime*func_xr2V(tt+dtime,  V+K3V,  xr2+K3xr2)
    K4Xs =    dtime*func_XsV(tt+dtime,  V+K3V,  xs+K3Xs)
    K4Xf =    dtime*func_XfV(tt+dtime,  V+K3V,  xf+K3Xf)
    K4g =     dtime*func_gV(tt+dtime,  Ca_i+K3Ca_i,  g+K3g ,V+K3V)
    K4Ca_i =  dtime*func_CaiV(tt+dtime,  Ca_i+K3Ca_i,  Ca_SR+K3Ca_SR,  d+K3d,  g+K3g,  V+K3V,  Na_i+K3Na_i,  f1+K3f1  ,f2+K3f2  ,f_Ca+K3f_Ca)
    K4Ca_SR = dtime*func_CaSRV(tt+dtime,  Ca_i+K3Ca_i,  Ca_SR+K3Ca_SR,  d+K3d,  g+K3g)
    K4Na_i =  dtime*func_NaiV(tt+dtime,  m+K3m,  h+K3h,  j+K3j,  V+K3V,  Na_i+K3Na_i,  Ca_i+K3Ca_i)
    K4V =     dtime*func_VV(tt+dtime,  V+K3V,  r+K3r,  q+K3q,  xr1+K3xr1,  xr2+K3xr2,  xs+K3Xs,  Ca_i+K3Ca_i,  d+K3d,  f1+K3f1,  f2+K3f2,  f_Ca+K3f_Ca,  Na_i+K3Na_i,  m+K3m,  h+K3h,  j+K3j  ,xf+K3Xf)

!!Gates, Concentrations and Voltage
    vauxi_exm(1,ipoin,1) = h + (K1h + 2.0_rp*K2h + 2.0_rp*K3h + K4h)*(1.0_rp/6.0_rp)
    vauxi_exm(2,ipoin,1) = j + (K1j + 2.0_rp*K2j + 2.0_rp*K3j + K4j)*(1.0_rp/6.0_rp)
    vauxi_exm(3,ipoin,1) = m + (K1m + 2.0_rp*K2m + 2.0_rp*K3m + K4m)*(1.0_rp/6.0_rp)
    vauxi_exm(4,ipoin,1) = d + (K1d + 2.0_rp*K2d + 2.0_rp*K3d + K4d)*(1.0_rp/6.0_rp)
    vauxi_exm(5,ipoin,1) = f_Ca + (K1f_Ca + 2.0_rp*K2f_Ca + 2.0_rp*K3f_Ca + K4f_Ca)*(1.0_rp/6.0_rp)
    vauxi_exm(6,ipoin,1) = f1 + (K1f1 + 2.0_rp*K2f1 + 2.0_rp*K3f1 + K4f1)*(1.0_rp/6.0_rp)
    vauxi_exm(7,ipoin,1) = f2 + (K1f2 + 2.0_rp*K2f2 + 2.0_rp*K3f2 + K4f2)*(1.0_rp/6.0_rp)
    vauxi_exm(8,ipoin,1) = r + (K1r + 2.0_rp*K2r + 2.0_rp*K3r + K4r)*(1.0_rp/6.0_rp)
    vauxi_exm(9,ipoin,1) = q + (K1q + 2.0_rp*K2q + 2.0_rp*K3q + K4q)*(1.0_rp/6.0_rp)
    vauxi_exm(10,ipoin,1) = xr1 + (K1xr1 + 2.0_rp*K2xr1 + 2.0_rp*K3xr1 + K4xr1)*(1.0_rp/6.0_rp)
    vauxi_exm(11,ipoin,1) = xr2 + (K1xr2 + 2.0_rp*K2xr2 + 2.0_rp*K3xr2 + K4xr2)*(1.0_rp/6.0_rp)
    vauxi_exm(12,ipoin,1) = xs + (K1Xs + 2.0_rp*K2Xs + 2.0_rp*K3Xs + K4Xs)*(1.0_rp/6.0_rp)
    vauxi_exm(13,ipoin,1) = xf + (K1Xf + 2.0_rp*K2Xf + 2.0_rp*K3Xf + K4Xf)*(1.0_rp/6.0_rp)
    vauxi_exm(14,ipoin,1) = g + (K1g + 2.0_rp*K2g + 2.0_rp*K3g + K4g)*(1.0_rp/6.0_rp)
    vconc(1,ipoin,1) = Ca_i + (K1Ca_i + 2.0_rp*K2Ca_i + 2.0_rp*K3Ca_i + K4Ca_i)*(1.0_rp/6.0_rp)
    vconc(2,ipoin,1) = Ca_SR + (K1Ca_SR + 2.0_rp*K2Ca_SR + 2.0_rp*K3Ca_SR + K4Ca_SR)*(1.0_rp/6.0_rp)
    vconc(3,ipoin,1) = Na_i + (K1Na_i + 2.0_rp*K2Na_i + 2.0_rp*K3Na_i + K4Na_i)*(1.0_rp/6.0_rp)
    !xioni = V + (K1V + 2.0_rp*K2V + 2.0_rp*K3V + K4V)*(1.0_rp/6.0_rp)
    xioni = -func_VV(tt,V,r,q,xr1,xr2,xs,Ca_i,d,f1,f2,f_Ca,Na_i,m,h,j,xf)
    xioni = xioni * 1000.0_rp
    df1 = (vauxi_exm(6,ipoin,1)-vauxi_exm(6,ipoin,2))/dtime
    dioni = 0.0_rp
    !write(*,*) df1
  
  !open(unit=1, file='voltage_Ventricular.dat')
  !do z=1,counter
  !    write(1,*) tspan(z), elmag(z)
  !end do
  !close(unit=1)
    !Moved to exm_updunk
    !vconc(:,ipoin,3)=vconc(:,ipoin,2)         
    !vconc(:,ipoin,2)=vconc(:,ipoin,1) 
    vauxi_exm(:,ipoin,2)=vauxi_exm(:,ipoin,1)   

end subroutine exm_scvent_ionicurrents


!################################# functions are defined here#########################
!#####################################################################################
real(rp) function func_CaiV(tt,Ca_i,Ca_SR,d,g,V,Na_i,f1,f2,f_Ca)
  use mod_exm_hipscventriconst


  real(rp), intent(in) :: tt, Ca_i, Ca_SR, d, g, V, Na_i, f1, f2, f_Ca
  real(rp) :: Ca_ibufc, E_Ca, I_CaL, I_NaCa, I_bCa, I_pCa, I_rel, I_up, I_leak

  !! Reversal potentials
  E_Ca = (0.5_rp*Rcons*Temp/F)*log(Ca_o/Ca_i)

  !! L-type Ca2+ current, I_CaL
  I_CaL = (g_CaL*4.0_rp*V*F**2.0_rp)*(Ca_i*safe_exp(2.0_rp*V*F/(Rcons*Temp))-0.341_rp*Ca_o)*d*f1*f2*f_Ca/(Rcons*Temp*(safe_exp(2.0_rp*V*F/(Rcons*Temp))-1.0_rp))  
  !! Na+/Ca2+ exchanger current, I_NaCa
  I_NaCa = K_NaCa*(safe_exp(Gamma*V*F/(Rcons*Temp))*Na_i**3.0_rp*Ca_o-safe_exp((Gamma-1.0_rp)*V*F/(Rcons*Temp))*Na_o**3.0_rp*Ca_i*Alpha)/((Km_Na_i**3.0_rp+Na_o**3.0_rp)*(Km_Ca+Ca_o)*(1.0_rp+K_sat*safe_exp((Gamma-1.0_rp)*V*F/(Rcons*Temp))))
    
  !! background currents
  I_bCa = g_bCa*(V-E_Ca)

  !! Ca2+ pump current, I_pCa
  I_pCa = g_pCa*Ca_i/(Ca_i+K_pCa)
  !! Ca2+ dynamics
  !Ventricular
  I_rel = ((a_rel*Ca_SR**2.0_rp/(b_rel**2.0_rp+Ca_SR**2.0_rp))+c_rel)*d*g*(0.0411_rp)


  I_up = Vmax_up/(1.0_rp+K_up**2.0_rp/Ca_i**2.0_rp)
  I_leak = V_leak*(Ca_SR-Ca_i)

  Ca_ibufc = 1.0_rp/(1.0_rp+Buf_c*K_Buf_c/(Ca_i+K_Buf_c)**2.0_rp)
  func_CaiV = Ca_ibufc*(I_leak-I_up+I_rel-(I_CaL+I_bCa+I_pCa-2.0_rp*I_NaCa)*C_m/(2.0_rp*V_c*F*1.0e-18_rp));
end function

!#####################################################################################
real(rp) function func_CaSRV(tt,Ca_i,Ca_SR,d,g)
  use mod_exm_hipscventriconst

  real(rp), intent(in) :: tt, Ca_i, Ca_SR, d, g
  real(rp) :: Ca_srbufsr, I_rel, I_up, I_leak

  !! Ca2+ dynamics
  !Ventricular
  I_rel = ((a_rel*Ca_SR**2.0_rp/(b_rel**2.0_rp+Ca_SR**2.0_rp))+c_rel)*d*g*(0.0411_rp)

  I_up = Vmax_up/(1.0_rp+K_up**2.0_rp/Ca_i**2.0_rp)
  I_leak = V_leak*(Ca_SR-Ca_i)

  !!
  Ca_srbufsr = 1.0_rp/(1.0_rp+Buf_sr*K_Buf_sr/(Ca_SR+K_Buf_sr)**2.0_rp)
  func_CaSRV = Ca_srbufsr*V_c*(I_up-(I_rel+I_leak))/V_sr

end function

!#####################################################################################
real(rp) function func_dV(tt,V,d)
  use mod_exm_hipscventriconst

  real(rp), intent(in) :: tt, V, d

  real(rp) :: d_inf, Alpha_d, beta_d, Gamma_d, taw_d
  !I_CaL, d gate
  !Ventricular
  d_inf = 1.0_rp/(1.0_rp+safe_exp(-(V*1000.0_rp+9.1_rp)/7.0_rp))

  Alpha_d = (1.4_rp/(1.0_rp+safe_exp((-35.0_rp-V*1000.0_rp)/13.0_rp)))+0.25_rp
  beta_d = 1.4_rp/(1.0_rp+safe_exp((V*1000.0_rp+5.0_rp)/5.0_rp))
  Gamma_d = 1.0_rp/(1.0_rp+safe_exp((50.0_rp-V*1000.0_rp)/20.0_rp))
  taw_d = (Alpha_d*beta_d+Gamma_d)/1000.0_rp

  func_dV = (d_inf-d)/taw_d;
end function 

!#####################################################################################
real(rp) function  func_f1V(tt,V,f1,Ca_i,df1)
  use mod_exm_hipscventriconst

  real(rp), intent(in) :: tt, V, f1, Ca_i, df1
  
  real(rp) :: f1_inf, taw_f1

  !! I_CaL, f1 gate
  !Ventricular
  f1_inf = 1.0_rp/(1.0_rp+safe_exp((V*1000.0_rp+26.0_rp)/3.0_rp))

  if (df1>0.0_rp) then
    taw_f1 = ((1102.5_rp*(safe_exp(-(((V*1000.0_rp+27.0_rp)**2.0_rp)/15.0_rp)**2.0_rp))+200.0_rp/(1.0_rp+safe_exp((13.0_rp-V*1000.0_rp)/10.0_rp))+180.0_rp/(1.0_rp+safe_exp((30.0_rp+V*1000.0_rp)/10.0_rp))+20.0_rp)*(1.0_rp+1433.0_rp*(Ca_i-50.0_rp*1e-6_rp)))/1000.0_rp
  else
    taw_f1 = ((1102.58_rp*(safe_exp(-(((V*1000.0_rp+27.0_rp)**2.0_rp)/15.0_rp)**2.0_rp))+200.0_rp/(1.0_rp+safe_exp((13.0_rp-V*1000.0_rp)/10.0_rp))+180.0_rp/(1.0_rp+safe_exp((30.0_rp+V*1000.0_rp)/10.0_rp))+20.0_rp))/1000.0_rp
  end if

  func_f1V = (f1_inf-f1)/taw_f1;

end function


!#####################################################################################
real(rp) function func_f2V(tt,V,f2)
  use mod_exm_hipscventriconst

  real(rp), intent(in) :: tt, V, f2
  
  real(rp) :: f2_inf, taw_f2

  !! I_CaL, f2 gate
  !Ventricular
  f2_inf = (0.67_rp/(1.0_rp+safe_exp((V*1000.0_rp+35.0_rp)/4.0_rp)))+0.33_rp
  taw_f2 =  ((600.0_rp*(safe_exp((-(V*1000.0_rp+25.0_rp)**2.0_rp)/170.0_rp))+31.0_rp/(1.0_rp+safe_exp((25.0_rp-V*1000.0_rp)/10.0_rp))+16.0_rp/(1.0_rp+safe_exp((30.0_rp+V*1000.0_rp)/10.0_rp))))/1000.0_rp

  func_f2V = (f2_inf-f2)/taw_f2;
end function

!#####################################################################################
! From here on, functions V = A (HISPCM)
!#####################################################################################

real(rp) function func_fCaV(tt,f_Ca,Ca_i,V)
  use mod_exm_hipscventriconst

  real(rp), intent(in) :: tt, f_Ca, Ca_i, V
  
  real(rp) :: Alpha_fCa, beta_fCa, Gamma_fCa, f_Ca_inf, constf_Ca
  real(rp), parameter :: taw_fCa = 0.002_rp
  
  !! I_CaL, f_Ca gate
  Alpha_fCa = 1.0_rp/(1.0_rp+(Ca_i/0.0006_rp)**8.0_rp)
  beta_fCa = 0.1_rp/(1.0_rp+safe_exp((Ca_i-0.0009_rp)/0.0001_rp))
  Gamma_fCa = 0.3_rp/(1.0_rp+safe_exp((Ca_i-0.00075_rp)/0.0008_rp))
  f_Ca_inf = (Alpha_fCa+beta_fCa+Gamma_fCa)/1.3156_rp
  
  if ((f_Ca_inf>f_Ca) .and. (V>-0.06_rp)) then
    constf_Ca = 0.0_rp
  else
    constf_Ca = 1.0_rp
  end if

  func_fCaV = constf_Ca*(f_Ca_inf-f_Ca)/taw_fCa
end function

!#####################################################################################
real(rp) function  func_XsV(tt,V,Xs)
  use mod_exm_hipscventriconst

  real(rp), intent(in) :: tt, V, Xs
  
  real(rp) :: xs_inf, Alpha_xs, beta_xs, taw_xs

  !I_Ks, Xs gate
  xs_inf = 1.0_rp/(1.0_rp+safe_exp((-V*1000.0_rp-20.0_rp)/16.0_rp))
  Alpha_xs =1100.0_rp/sqrt(1.0_rp+safe_exp((-V*1000.0_rp-10.0_rp)/6.0_rp))
  beta_xs = 1.0_rp/(1.0_rp+safe_exp((V*1000.0_rp-60.0_rp)/20.0_rp))
  taw_xs = (Alpha_xs*beta_xs)/1000.0_rp

  func_XsV = (xs_inf-Xs)/taw_xs
end function

!#####################################################################################
real(rp) function func_gV(tt,Ca_i,g,V)
  use mod_exm_hipscventriconst

  real(rp), intent(in) :: tt, Ca_i, g, V
  real(rp) :: g_inf
  real(rp), parameter :: taw_g = 0.002_rp
  real(rp) :: constii

  if (Ca_i<=0.00035_rp) then
    g_inf = 1.0_rp/(1.0_rp+(Ca_i/0.00035_rp)**6.0_rp)
  else
    g_inf = 1.0_rp/(1.0_rp+(Ca_i/0.00035_rp)**16.0_rp)
  end if

  if ((g_inf>g) .and. (V>-0.06_rp)) then
    constii = 0.0_rp
  else
    constii = 1.0_rp
  end if

  func_gV = constii*(g_inf-g)/taw_g
end function

!#####################################################################################
real(rp) function func_hV(tt,V,h)
  use mod_exm_hipscventriconst

  real(rp), intent(in) :: tt, V, h
  
  real(rp) :: h_inf, Alpha_h, beta_h, taw_h

  !! I_Na, h gate
  h_inf = 1.0_rp/(sqrt(1.0_rp+safe_exp((V*1000.0_rp+72.1_rp)/5.7_rp)))

  if (V<-0.04_rp) then
    Alpha_h = 0.057_rp*safe_exp(-(V*1000.0_rp+80.0_rp)/6.8_rp)
    beta_h = 2.7_rp*safe_exp(0.079_rp*V*1000.0_rp)+3.1_rp*(10.0_rp**5.0_rp)*safe_exp(0.3485_rp*V*1000.0_rp)
  else
    Alpha_h = 0.0_rp
    beta_h = 0.77_rp/(0.13_rp*(1.0_rp+safe_exp((1000.0_rp*V+10.66_rp)/(-11.1_rp))))
  end if
  if (V<-0.04_rp) then
    taw_h = 1.5_rp/((Alpha_h+beta_h)*1000.0_rp)
  else
    taw_h = 2.542_rp/1000.0_rp
  end if
  func_hV = (h_inf-h)/taw_h
end function

!#####################################################################################
real(rp) function func_jV(tt,V,j)
  use mod_exm_hipscventriconst

  real(rp), intent(in) :: tt, V, j
  
  real(rp) :: j_inf, Alpha_j, beta_j, taw_j

  !! I_Na, j gate
  j_inf = 1.0_rp/(sqrt(1.0_rp+safe_exp((1000.0_rp*V+72.1_rp)/5.7_rp)))

  if (V<-0.04_rp) then
    Alpha_j = (-25428.0_rp*safe_exp(0.2444_rp*V*1000.0_rp)-6.948_rp*(1e-6_rp)*safe_exp(-0.04391_rp*V*1000.0_rp))*(V*1000.0_rp+37.78_rp)/(1.0_rp+safe_exp(0.311_rp*(V*1000.0_rp+79.23_rp)))
    beta_j =  0.02424_rp*safe_exp(-0.01052_rp*V*1000.0_rp)/(1.0_rp+safe_exp(-0.1378_rp*(V*1000.0_rp+40.14_rp)))
  else
    Alpha_j = 0.0_rp
    beta_j =  0.6_rp*safe_exp(0.057_rp*V*1000.0_rp)/(1.0_rp+safe_exp(-0.1_rp*(V*1000.0_rp+32.0_rp)))
  end if

  taw_j = 7.0_rp/((Alpha_j+beta_j)*1000.0_rp)
  func_jV = (j_inf-j)/taw_j
end function

!#####################################################################################
real(rp) function func_mV(tt,V,m)
  use mod_exm_hipscventriconst

  real(rp), intent(in) :: tt, V, m
  
  real(rp) :: m_inf, Alpha_m, beta_m, taw_m

  !! I_Na, m gate
  m_inf =  (1.0_rp+safe_exp((-34.1_rp-V*1000.0_rp)/5.9_rp))**(-1.0_rp/3.0_rp)
  Alpha_m =  1.0_rp/(1.0_rp+safe_exp((-60.0_rp-V*1000.0_rp)/5.0_rp))
  beta_m =  0.1_rp/(1.0_rp+safe_exp((V*1000.0_rp+35.0_rp)/5.0_rp))+0.1_rp/(1.0_rp+safe_exp((V*1000.0_rp-50.0_rp)/200.0_rp))
  taw_m = (Alpha_m*beta_m)/1000.0_rp

  func_mV = (m_inf-m)/taw_m
end function


!#####################################################################################
real(rp) function func_NaiV(tt,m,h,j,V,Na_i,Ca_i)
  use mod_exm_hipscventriconst

  real(rp), intent(in) :: tt, m, h, j, V, Na_i, Ca_i
  real(rp) :: E_Na, I_bNa, I_Na, I_NaK, I_NaCa
 
  !! Reversal potentials
  E_Na = (Rcons*Temp/F)*log(Na_o/Na_i)

  !! background currents
  I_bNa = g_bNa*(V-E_Na)
  I_Na = g_Na*(m**3.0_rp)*h*j*(V-E_Na)

  !! Na+/K+ pump current, I_NaK
  I_NaK = ((P_NaK*K_o*Na_i/(K_o+K_mk))/(Na_i+K_mNa))/(1.0_rp+0.1245_rp*safe_exp(-0.1_rp*V*F/(Rcons*Temp))+0.353_rp*safe_exp(-V*F/(Rcons*Temp)))

  !! Na+/Ca2+ exchanger current, I_NaCa
  I_NaCa = K_NaCa*(safe_exp(Gamma*V*F/(Rcons*Temp))*(Na_i**3.0_rp)*Ca_o-safe_exp((Gamma-1.0_rp)*V*F/(Rcons*Temp))*(Na_o**3.0_rp)*Ca_i*Alpha)/((Km_Na_i**3.0_rp+Na_o**3.0_rp)*(Km_Ca+Ca_o)*(1.0_rp+K_sat*safe_exp((Gamma-1.0_rp)*V*F/(Rcons*Temp))))

  !! sodium dynamics
  func_NaiV = -C_m*(I_Na+I_bNa+3.0_rp*I_NaK+3.0_rp*I_NaCa)/(F*V_c*1e-18_rp)
end function

!#####################################################################################
real(rp) function func_qV(tt,V,q)
  use mod_exm_hipscventriconst

  real(rp), intent(in) :: tt, V, q
  
  real(rp) :: q_inf, taw_q

  !! I_to, q gate
  q_inf = 1.0_rp/(1.0_rp+safe_exp((V*1000.0_rp+53.0_rp)/13.0_rp))
  taw_q = (39.102_rp/(0.57_rp*safe_exp(-0.08_rp*(V*1000.0_rp+44.0_rp))+0.065_rp*safe_exp(0.1_rp*(V*1000.0_rp+45.93_rp)))+6.06_rp)/1000.0_rp

  func_qV = (q_inf-q)/taw_q
end function

!#####################################################################################
real(rp) function func_rV(tt,V,r)
  use mod_exm_hipscventriconst

  real(rp), intent(in) :: tt, V, r
  
  real(rp) :: r_inf, taw_r

  !I_to, r gate
  r_inf = 1.0_rp/(1.0_rp+safe_exp((22.3_rp-V*1000.0_rp)/18.75_rp))
  taw_r = (14.405516_rp/(1.037_rp*safe_exp(0.09_rp*(V*1000.0_rp+30.61_rp))+0.369_rp*safe_exp(-0.12_rp*(V*1000.0_rp+23.84_rp)))+2.75352_rp)/1000.0_rp

  func_rV = (r_inf-r)/taw_r
end function

!#####################################################################################
real(rp) function func_VV(tt,V,r,q,xr1,xr2,Xs,Ca_i,d,f1,f2,f_Ca,Na_i,m,h,j,Xf)
  use mod_exm_hipscventriconst

  real(rp), intent(in) :: tt, V, r, q, xr1, xr2, &
  Xs, Ca_i, d, f1, f2, f_Ca, Na_i, m, h, j, Xf
  
  real(rp) :: Alpha_K1, beta_K1, x_K1_inf, E_Na, E_K, E_Ks, E_Ca, E_f, &
      I_CaL, I_to, I_Kr, I_Ks, I_K1, I_f, I_NaK, I_NaCa, I_bNa, I_bCa, I_Na, &
      I_pCa

  !! Reversal potentials
  E_Na = (Rcons*Temp/F)*log(Na_o/Na_i)
  E_K = (Rcons*Temp/F)*log(K_o/K_i)
  E_Ks = (Rcons*Temp/F)*log((K_o+P_kna*Na_o)/(K_i+P_kna*Na_i))
  E_Ca = (0.5_rp*Rcons*Temp/F)*log(Ca_o/Ca_i)
  E_f = -0.017_rp


  !! L-type Ca2+ current, I_CaL
  I_CaL = (g_CaL*4.0_rp*V*F**2.0_rp)*(Ca_i*safe_exp(2.0_rp*V*F/(Rcons*Temp))-0.341_rp*Ca_o)*d*f1*f2*f_Ca/(Rcons*Temp*(safe_exp(2.0_rp*V*F/(Rcons*Temp))-1.0_rp))
  !!
  !! Transient outward current, I_to
  I_to = g_to*r*q*(V-E_K)

  !! Rapid delayed rectifier K+ current, I_Kr
  I_Kr = g_Kr*sqrt(K_o/5.4_rp)*xr1*xr2*(V-E_K)
  !! Slow delayed rectifier K+ current, I_Ks
  I_Ks = g_Ks*Xs**2.0_rp*(1.0_rp+0.6_rp/(1.0_rp+(3.8_rp*1e-5_rp/Ca_i)**1.4_rp))*(V-E_Ks)

  !! Inward rectifier K+ current, I_K1
  Alpha_K1 = 3.91_rp/(1.0_rp+safe_exp(0.5942_rp*(V*1000.0_rp-E_K*1000.0_rp-200.0_rp)))
  beta_K1 = (-1.509_rp*safe_exp(0.0002_rp*(V*1000.0_rp-E_K*1000.0_rp+100.0_rp))+safe_exp(0.5886_rp*(V*1000.0_rp-E_K*1000.0_rp-10.0_rp)))/(1.0_rp+safe_exp(0.4547_rp*(V*1000.0_rp-E_K*1000.0_rp)))
  x_K1_inf = Alpha_K1/(Alpha_K1+beta_K1)
  I_K1 = g_K1*x_K1_inf*sqrt(K_o/5.4_rp)*(V-E_K)

  !! Hyperpolarization activated funny current, I_f
  I_f = g_f*Xf*(V-E_f)

  !! Na+/K+ pump current, I_NaK
  I_NaK = ((P_NaK*K_o*Na_i/(K_o+K_mk))/(Na_i+K_mNa))/(1.0_rp+0.1245_rp*safe_exp(-0.1_rp*V*F/(Rcons*Temp))+0.353_rp*safe_exp(-V*F/(Rcons*Temp)));

  !! Na+/Ca2+ exchanger current, I_NaCa
  I_NaCa = K_NaCa*(safe_exp(Gamma*V*F/(Rcons*Temp))*Na_i**3.0_rp*Ca_o-safe_exp((Gamma-1.0_rp)*V*F/(Rcons*Temp))*Na_o**3.0_rp*Ca_i*Alpha)/((Km_Na_i**3.0_rp+Na_o**3.0_rp)*(Km_Ca+Ca_o)*(1.0_rp+K_sat*safe_exp((Gamma-1.0_rp)*V*F/(Rcons*Temp))))

  !! background currents
  I_bNa = g_bNa*(V-E_Na)
  I_bCa = g_bCa*(V-E_Ca)
  I_Na = g_Na*m**3.0_rp*h*j*(V-E_Na)

  !! Ca2+ pump current, I_pCa
  I_pCa = g_pCa*Ca_i/(Ca_i+K_pCa)

  !! Membrane potential
  func_VV = -(I_K1+I_to+I_Kr+I_Ks+I_CaL+I_NaK+I_Na+I_NaCa+I_pCa+I_f+I_bNa+I_bCa)
end function

!#####################################################################################

real(rp) function func_XfV(tt,V,Xf)
  use mod_exm_hipscventriconst

  real(rp), intent(in) :: tt, V, Xf
  
  real(rp) :: xf_inf, taw_f

  !I_f, X_f gate
  xf_inf = 1.0_rp/(1.0_rp+safe_exp((V*1000.0_rp+77.85_rp)/5.0_rp))
  taw_f = 1900.0_rp/((1.0_rp+safe_exp((V*1000.0_rp+15.0_rp)/10.0_rp))*1000.0_rp)

  func_XfV = (xf_inf-Xf)/taw_f
end function

!#####################################################################################
real(rp) function func_xr1V(tt,V,xr1)
  use mod_exm_hipscventriconst

  real(rp), intent(in) :: tt, V, xr1

  real(rp), parameter :: Qcons = 2.3_rp
  real(rp) :: V_half, xr1_inf, Alpha_xr1, beta_xr1, taw_xr1

  V_half = 1000.0_rp*(-(Rcons*Temp/(F*Qcons))*log(((1.0_rp+Ca_o/2.6_rp)**4.0_rp)/(L_0*(1.0_rp+Ca_o/0.58_rp)**4.0_rp))-0.019_rp)
  xr1_inf = 1.0_rp/(1.0_rp+safe_exp((V_half - V*1000.0_rp)/4.9_rp))
  Alpha_xr1 = 450.0_rp/(1.0_rp+safe_exp((-45.0_rp-V*1000.0_rp)/10.0_rp))
  beta_xr1 = 6.0_rp/(1.0_rp+safe_exp((V*1000.0_rp+30.0_rp)/11.5_rp))
  taw_xr1 = (Alpha_xr1*beta_xr1)/1000.0_rp

  func_xr1V = (xr1_inf-xr1)/taw_xr1
end function

!#####################################################################################
real(rp) function func_xr2V(tt,V,xr2)
  use mod_exm_hipscventriconst

  real(rp), intent(in) :: tt, V, xr2
  
  real(rp) :: xr2_inf, Alpha_xr2, beta_xr2, taw_xr2

  !! I_Kr, xr2 gate
  xr2_inf = 1.0_rp/(1.0_rp+safe_exp((V*1000.0_rp+88.0_rp)/50.0_rp))
  Alpha_xr2 = 3.0_rp/(1.0_rp+safe_exp((-V*1000.0_rp-60.0_rp)/20.0_rp))
  beta_xr2 = 1.12_rp/(1.0_rp+safe_exp((V*1000.0_rp-60.0_rp)/20.0_rp))
  taw_xr2 = (Alpha_xr2*beta_xr2)/1000.0_rp

  func_xr2V = (xr2_inf-xr2)/taw_xr2
end function


end module mod_exm_scvent_ionicurrents
