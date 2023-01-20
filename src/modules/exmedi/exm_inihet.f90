!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!!------------------------------------------------------------------------
!!> @addtogroup Exmedi
!!> @{
!!> @file    exm_inihet.f90
!!> @author  Jazmin Aguado-Sierra
!!> @brief   Initial condition setup for TenTuscher-Panfilov 2006 heterogeneous model
!!> @details Obtains initial conditions of Normal or Heart Failure cell at  70 bpm or 857 ms \n
!!!    Otherwise, it obtains the initial conditions from single cell runs in ONECEL.F90 AND OCEIHE.F90 \n
!!> @} 
!!!-----------------------------------------------------------------------
!subroutine exm_inihet(ipoin,mat)
!!subroutine exm_inihet(kmodel_ipoin,ipoin,mat)
!
!  use      def_master
!  use      def_elmtyp
!  use      def_domain
!  use      def_kermod,           only :  kfl_celltype_fun
!
!  use      def_exmedi
!
!  implicit none
!  integer(ip), intent(in) :: ipoin, mat
!  integer(ip)             :: ituss_exm
!  real(rp)    ::  vaux1, vaux2, vaux3, rhsx
!  !integer(ip) ::  mat
!
!  if(kfl_celltype_fun == 0_ip) then
!      call runend("exm_inihet: Cell type field is required")
!  end if 
!
!     !
!     ! Load initial conditions for the potentials
!     !
!     !mat=nodemat(ipoin)
!     ituss_exm = int(celty_exm(ipoin),KIND=ip)
!     !write(*,*) ituss_exm
!     if(ituss_exm == EXM_CELLTYPE_EPI) then   !epicardial
!        elmag(ipoin,1:3) = vminimate_exm(ituss_exm,mat)          ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
!        vauxi_exm(1,ipoin,1:3)  = vauxi_exm_initial(1, ituss_exm,mat)          !Variable m
!        vauxi_exm(2,ipoin,1:3)  = vauxi_exm_initial(2, ituss_exm,mat)          !Variable h 
!        vauxi_exm(3,ipoin,1:3)  = vauxi_exm_initial(3, ituss_exm,mat)          !Variable j
!        vauxi_exm(4,ipoin,1:3)  = vauxi_exm_initial(4, ituss_exm,mat)          !Variable d
!        vauxi_exm(5,ipoin,1:3)  = vauxi_exm_initial(5, ituss_exm,mat)          !Variable f
!        vauxi_exm(12,ipoin,1:3) = vauxi_exm_initial(12,ituss_exm,mat)        !Variable f2
!        vauxi_exm(6,ipoin,1:3)  = vauxi_exm_initial(6, ituss_exm,mat)          ! Variable fCa
!        vauxi_exm(7,ipoin,1:3)  = vauxi_exm_initial(7, ituss_exm,mat)          !Variable r
!        vauxi_exm(8,ipoin,1:3)  = vauxi_exm_initial(8, ituss_exm,mat)          !Variable s
!        vauxi_exm(9,ipoin,1:3)  = vauxi_exm_initial(9, ituss_exm,mat)          !Variable xs 
!        vauxi_exm(10,ipoin,1:3) = vauxi_exm_initial(10,ituss_exm,mat)        !Variable xr1
!        vauxi_exm(11,ipoin,1:3) = vauxi_exm_initial(11,ituss_exm,mat)        !Variable xr2 
!        vconc(1,ipoin,1:3) = vconc_initial(1,ituss_exm,mat)              !Variable Cai
!        vconc(2,ipoin,1:3) = vconc_initial(2,ituss_exm,mat)              !Variable CaSR
!        vconc(3,ipoin,1:3) = vconc_initial(3,ituss_exm,mat)              !Variable Nai
!        vconc(4,ipoin,1:3) = vconc_initial(4,ituss_exm,mat)              !Variable Ki
!        vconc(5,ipoin,1:3) = vconc_initial(5,ituss_exm,mat)              !Variable CaSS
!        vconc(9,ipoin,1:3) = vconc_initial(9,ituss_exm,mat)              !Variable Rprime
!     else if(ituss_exm == EXM_CELLTYPE_ENDO) then  !endocardial
!        elmag(ipoin,1:3) = vminimate_exm(ituss_exm,mat)          ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
!        vauxi_exm(1,ipoin,1:3) = vauxi_exm_initial(1,  ituss_exm,mat)          !Variable m
!        vauxi_exm(2,ipoin,1:3) = vauxi_exm_initial(2,  ituss_exm,mat)          !Variables h 
!        vauxi_exm(3,ipoin,1:3) = vauxi_exm_initial(3,  ituss_exm,mat)          !Variable j
!        vauxi_exm(4,ipoin,1:3) = vauxi_exm_initial(4,  ituss_exm,mat)          !Variable d
!        vauxi_exm(5,ipoin,1:3) = vauxi_exm_initial(5,  ituss_exm,mat)          !Variable f
!        vauxi_exm(12,ipoin,1:3) = vauxi_exm_initial(12,ituss_exm,mat)        !Variable f2
!        vauxi_exm(6,ipoin,1:3) = vauxi_exm_initial(6,  ituss_exm,mat)          !Variable fCa
!        vauxi_exm(7,ipoin,1:3) = vauxi_exm_initial(7,  ituss_exm,mat)          !Variable r
!        vauxi_exm(8,ipoin,1:3) = vauxi_exm_initial(8,  ituss_exm,mat)          !Variable s
!        vauxi_exm(9,ipoin,1:3) = vauxi_exm_initial(9,  ituss_exm,mat)          !Variable xs 
!        vauxi_exm(10,ipoin,1:3) = vauxi_exm_initial(10,ituss_exm,mat)        !Variable xr1
!        vauxi_exm(11,ipoin,1:3) = vauxi_exm_initial(11,ituss_exm,mat)        !Variable xr2 
!        vconc(1,ipoin,1:3) = vconc_initial(1,ituss_exm,mat)              !Variable Cai
!        vconc(2,ipoin,1:3) = vconc_initial(2,ituss_exm,mat)              !Variable CaSR
!        vconc(3,ipoin,1:3) = vconc_initial(3,ituss_exm,mat)              !Variable Nai
!        vconc(4,ipoin,1:3) = vconc_initial(4,ituss_exm,mat)              !Variable Ki
!        vconc(5,ipoin,1:3) = vconc_initial(5,ituss_exm,mat)              !Variable CaSS
!        vconc(9,ipoin,1:3) = vconc_initial(9,ituss_exm,mat)              !Variable Rprime
!     else if(ituss_exm == EXM_CELLTYPE_MID) then  !mid
!        elmag(ipoin,1:3) = vminimate_exm(ituss_exm,mat)          ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
!        vauxi_exm(1,ipoin,1:3) =  vauxi_exm_initial(1, ituss_exm,mat)          !Variable m
!        vauxi_exm(2,ipoin,1:3) =  vauxi_exm_initial(2, ituss_exm,mat)          !Variables h 
!        vauxi_exm(3,ipoin,1:3) =  vauxi_exm_initial(3, ituss_exm,mat)          !Variable j
!        vauxi_exm(4,ipoin,1:3) =  vauxi_exm_initial(4, ituss_exm,mat)          !Variable d
!        vauxi_exm(5,ipoin,1:3) =  vauxi_exm_initial(5, ituss_exm,mat)          !Variable f
!        vauxi_exm(12,ipoin,1:3) = vauxi_exm_initial(12,ituss_exm,mat)        !Variable f2
!        vauxi_exm(6,ipoin,1:3) =  vauxi_exm_initial(6, ituss_exm,mat)          ! Variable fCa
!        vauxi_exm(7,ipoin,1:3) =  vauxi_exm_initial(7, ituss_exm,mat)          !Variable r
!        vauxi_exm(8,ipoin,1:3) =  vauxi_exm_initial(8, ituss_exm,mat)          !Variable s
!        vauxi_exm(9,ipoin,1:3) =  vauxi_exm_initial(9, ituss_exm,mat)          !Variable xs 
!        vauxi_exm(10,ipoin,1:3) = vauxi_exm_initial(10,ituss_exm,mat)        !Variable xr1
!        vauxi_exm(11,ipoin,1:3) = vauxi_exm_initial(11,ituss_exm,mat)        !Variable xr2 
!        vconc(1,ipoin,1:3) = vconc_initial(1,ituss_exm,mat)              !Variable Cai
!        vconc(2,ipoin,1:3) = vconc_initial(2,ituss_exm,mat)              !Variable CaSR
!        vconc(3,ipoin,1:3) = vconc_initial(3,ituss_exm,mat)              !Variable Nai
!        vconc(4,ipoin,1:3) = vconc_initial(4,ituss_exm,mat)              !Variable Ki
!        vconc(5,ipoin,1:3) = vconc_initial(5,ituss_exm,mat)              !Variable CaSS
!        vconc(9,ipoin,1:3) = vconc_initial(9,ituss_exm,mat)              !Variable Rprime
!     end if
!     ! Variable of activation O (related with I_rel current)
!     ! nconc = 10
!     ! "A model for ventricular tissue", Ten Tusscher, Page H1587 (col 1)
!     !  kcasr = max_sr - (max_sr - min_sr)/(1+(pow((EC/Ca_SR), 2)))
!     vaux1 = 1.0_rp + ((1.5_rp / vconc(2,ipoin,2))*(1.5_rp / vconc(2,ipoin,2)))  !conc2 = CaSRfree = CaSR
!     vaux1 = 1.0_rp / vaux1
!     vaux1 = 2.5_rp - (1.5_rp * vaux1)  !kCaSR
!     !     O = ( k1*pow(Ca_ss, 2)*R_prime)/(k3+ k1*pow(Ca_ss, 2))      
!     vaux3 = 0.15_rp / vaux1   !K1
!     vaux2 = 0.045_rp * vaux1  !K2
!     rhsx = 1.0_rp / (0.06_rp + (vaux3*vconc(5,ipoin,2)*vconc(5,ipoin,2)))      !O for calcium dynamics
!     rhsx = vaux3 * vconc(5,ipoin,2) * vconc(5,ipoin,2) * vconc(9,ipoin,2) * rhsx !O
!     vconc(10,ipoin,1:3) = rhsx 
!
!  !write(991,*) vconc(1,:,1)  !JAS
!end subroutine exm_inihet
