!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_courteland_calcium.f90
!> @author  Eva Casoni
!> @brief   Cellular electro-mechanical coupling mechanisms with Land model 
!> @date   16/NOV/1966
!> @details Bidirectional coupling between the models of electrophysiology and active stress  \n
!!    replace the algebraic formulation of the troponin in Courtemanche model by the evoluion of calcium \n
!!    bound to troponin from the Land model \n
!!    Reference: Gerach, T. et al: "Electro-Mechanical Whole-Heart Digital Twins: A Fully
!!               Coupled Multi-Physics Approach" (2021)
!> @} 
!!-----------------------------------------------------------------------
subroutine exm_courteland_calcium(kmcmdn,cmdnmax,trpnmax,dtimeEP,ipoin,cai,sac,volt,ituss_exm,imate)

  use      def_parame
  use      def_master
  use      def_elmtyp
  use      def_domain
  use      def_exmedi
  use      mod_eccoupling, only: troponin_ecc, props_ecc

  implicit none

  real(rp), intent(in) :: kmcmdn,cmdnmax,trpnmax,dtimeEP,volt
  real(rp), intent(out) :: cai, sac
  real(rp) :: vaux1,vaux3,vaux4,rhsx1,rhsx2,rhsx,val0,k_TRPN,n_TRPN,Ca50,Ca50_ref,beta_1,lambda,  &
   C_tens(3,3),dCaTRPN,lambda0,gsac,esac,Ca50_fact
  real(rp) :: V_up, V_rel, V_i, F   !parameters from Courtemanche model

  integer(ip) , intent(in) :: ipoin, imate, ituss_exm
  integer(ip) :: idime,jdime,kdime

  ! Declare parameters needed from Courtemanche model
  V_up = 1109.52_rp
  V_rel = 96.48_rp
  V_i = 13668.0_rp
  F = 96.4867_rp

  ! Declare some parameters of the Land model
  k_TRPN =   props_ecc(1 , imate) * 0.001_rp !0.1_rp
  n_TRPN =   props_ecc(2 , imate) !2.0_rp
  Ca50_ref = props_ecc(3 , imate)/1000.0_rp !0.805_rp
  beta_1 =   props_ecc(16, imate) !-2.4_rp
  Ca50_fact = props_ecc(20, imate)!1.0_rp by defaults

  ! Recover the right Cauchy-Green strain tensor
  C_tens = 0.0_rp
  do idime = 1,ndime
    do jdime = 1,ndime
      do kdime = 1,ndime
        C_tens(idime,jdime) = C_tens(idime,jdime)+gdepo(kdime,idime,ipoin)*                   &
         gdepo(kdime,jdime,ipoin)
      end do
    end do
  end do

  ! Recover the stretch in the fibre direction
  lambda = 0.0_rp
  do idime = 1,ndime
    do jdime = 1,ndime
      lambda = lambda+fiber(idime,ipoin)*C_tens(idime,jdime)*fiber(jdime,ipoin)
    end do
  end do
  lambda = SQRT(lambda)

  ! Calculate calcium sensitivity
  lambda0 = lambda
  lambda = MIN(1.2_rp,lambda)
  Ca50 = (Ca50_ref+beta_1*(lambda-1.0_rp))*Ca50_fact


  ! Calcium calcium in the coupled model
  vaux1 = (kmcmdn+vconc(1,ipoin,2)) * (kmcmdn+vconc(1,ipoin,2))
  vaux3 = 1.0_rp/(1.0_rp + (cmdnmax*kmcmdn/vaux1))

  vaux4 = (vconc(1,ipoin,2)/(Ca50_ref*Ca50_fact))**n_TRPN * (1.0_rp-(troponin_ecc(ipoin,1)*0.001_rp))
  dCaTRPN = k_TRPN *( vaux4 - (troponin_ecc(ipoin,1)*0.001_rp))

  rhsx1 = 2.0_rp*vicel_exm(9,ipoin) - vicel_exm(10,ipoin) - vicel_exm(7,ipoin) - vicel_exm(12,ipoin)
  rhsx2 = V_up*(vicel_exm(16,ipoin) - vicel_exm(17,ipoin)) + vicel_exm(14,ipoin)*V_rel
  rhsx = vaux3*(rhsx1/(2.0_rp*F*V_i) + rhsx2/V_i- dCaTRPN)
  val0 = vconc(1,ipoin,2)
  cai = val0 + dtimeEP*rhsx

  ! Calculate SAC current
  if (lambda0 >= 1.0_rp) then
    gsac = 0.1_rp
    esac = 0.0_rp
    sac = gsac*((lambda0-1.0_rp)/(1.1_rp-1.0_rp))*(volt-esac)
  else
    sac = 0.0_rp  
  end if

end subroutine exm_courteland_calcium
