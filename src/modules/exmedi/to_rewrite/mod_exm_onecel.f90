!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    mod_exm_onecel.f90
!> @author  mixed
!> @date    2019-11-16
!> @brief   mod_exm_onecel
!> @details mod_exm_onecel
!> @}
!-----------------------------------------------------------------------
!module mod_exm_onecel
!
!   use def_master
!   use def_elmtyp
!   use def_domain
!   use def_exmedi
!
!   implicit none
!
!   integer(ip), parameter  :: &
!      EXM_ONECEL_NORMAL_70BPM = 0_ip, &
!      EXM_ONECEL_HF_70BPM = 1_ip
!
!   public :: exm_onecel
!   private :: exm_init_voltages
!
!contains
!
!   subroutine exm_onecel(mat)
!      !------------------------------------------------------------------------
!      !> @addtogroup Exmedi
!      !> @{
!      !> @file    exm_onecel.f90
!      !> @date    03/10/2012
!      !> @author  Mariano Vazquez
!      !> @brief   One cell simulations to obtain initial conditions
!      !> @details Runs for changes of cycle length and drug administration
!      !> @}
!      !------------------------------------------------------------------------
!
!      implicit none
!      integer(ip), intent(in) :: mat
!      integer(ip) :: n
!
!      !call runend("exm_onecel needs to be recalculated before use")
!
!      n = mat
!      if (kfl_hfmodmate_exm(mat) == 0) then  !if it's Normal
!
!         if ((moneclmate_exm(2, n) == 1000) .and. (kfl_drugsmate_exm(n) == 0)) then   !Normal at 70 bpm
!            call exm_init_voltages(EXM_ONECEL_NORMAL_70BPM, mat)
!
!            !return
!         else if ((moneclmate_exm(2, n) /= 857) .or. (kfl_drugsmate_exm(n) == 1)) then      !Normal with different heart rate
!
!            call runend("EXM_ONECEL: exm_oceihe currently is not working correctly!")
!
!            !        elmlo_exm(1) = -85.23_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
!            !        elmlo_exm(2) = -85.23_rp
!            !     !do nbeat = 1, beat
!            !        ituss_exm = 3    !epicardial
!            !        vaulo_exm(1,1) = 0.00172_rp           !Variable m
!            !        vaulo_exm(2,1) = 0.7444_rp            !Variable h
!            !        vaulo_exm(3,1) = 0.7045_rp            !Variable j
!            !        vaulo_exm(4,1) = 0.00003373_rp        !Variable d
!            !        vaulo_exm(5,1) = 0.7888_rp            !Variable f
!            !        vaulo_exm(12,1) = 0.9755_rp           !Variable f2
!            !        vaulo_exm(6,1) = 0.9953_rp            !Variable fCa
!            !        vaulo_exm(7,1) = 0.0000000242_rp      !Variable r
!            !        vaulo_exm(8,1) = 0.999998_rp          !Variable s
!            !        vaulo_exm(9,1) = 0.0095_rp            !Variable xs
!            !        vaulo_exm(10,1) = 0.00621_rp          !Variable xr1
!            !        vaulo_exm(11,1) = 0.4712_rp           !Variable xr2
!            !        vcolo_exm(1,1) = 0.000126_rp          !Variable Cai
!            !        vcolo_exm(2,1) = 3.64_rp              !Variable CaSR
!            !        vcolo_exm(3,1) = 8.608_rp             !Variable Nai
!            !        vcolo_exm(4,1) = 136.89_rp            !Variable Ki
!            !        vcolo_exm(5,1) = 0.00036_rp           !Variable CaSS
!            !        vcolo_exm(9,1) = 0.9073_rp            !Variable Rprime
!            !        !call exm_oceihe
!
!            !        ituss_exm = 1    !endocardial
!            !        vaulo_exm(1,1) = 0.0016_rp           !Variable m
!            !        vaulo_exm(2,1) = 0.7573_rp           !Variables h
!            !        vaulo_exm(3,1) = 0.7225_rp           !Variable j
!            !        vaulo_exm(4,1) = 0.00003164_rp       !Variable d
!            !        vaulo_exm(5,1) = 0.8009_rp           !Variable f
!            !        vaulo_exm(12,1) = 0.9778_rp          ! Variable f2
!            !        vaulo_exm(6,1) = 0.99530_rp          ! Variable fCa
!            !        vaulo_exm(7,1) = 0.00000002235_rp    !Variable r
!            !        vaulo_exm(8,1) = 0.3212_rp           !Variable s
!            !        vaulo_exm(9,1) = 0.0087_rp           !Variable xs
!            !        vaulo_exm(10,1) = 0.0045_rp          !Variable xr1
!            !        vaulo_exm(11,1) = 0.4760_rp          !Variable xr2
!            !        vcolo_exm(1,1) = 0.00013_rp          !Variable Cai
!            !        vcolo_exm(2,1) = 3.7150_rp           !Variable CaSR
!            !        vcolo_exm(3,1) = 10.355_rp           !Variable Nai
!            !        vcolo_exm(4,1) = 138.4_rp            !Variable Ki
!            !        vcolo_exm(5,1) = 0.00036_rp          !Variable CaSS
!            !        vcolo_exm(9,1) = 0.9068_rp           !Variable Rprime
!            !        !call exm_oceihe
!
!            !        ituss_exm = 2    ! MIDMYOCARDIAL
!            !        vaulo_exm(1,1) = 0.0017_rp           !Variable m
!            !        vaulo_exm(2,1) = 0.749_rp            !Variables h
!            !        vaulo_exm(3,1) = 0.6788_rp           !Variable j
!            !        vaulo_exm(4,1) = 0.00003288_rp       !Variable d
!            !        vaulo_exm(5,1) = 0.7026_rp           !Variable f
!            !        vaulo_exm(12,1) = 0.9526_rp          !Variable f2
!            !        vaulo_exm(6,1) = 0.99420_rp          ! Variable fCa
!            !        vaulo_exm(7,1) = 0.00000002347_rp    !Variable r
!            !        vaulo_exm(8,1) = 1.0_rp              !Variable s
!            !        vaulo_exm(9,1) = 0.0174_rp           !Variable xs
!            !        vaulo_exm(10,1) = 0.0165_rp          !Variable xr1
!            !        vaulo_exm(11,1) = 0.473_rp           !Variable xr2
!            !        vcolo_exm(1,1) = 0.000153_rp         !Variable Cai
!            !        vcolo_exm(2,1) = 4.2720_rp           !Variable CaSR
!            !        vcolo_exm(3,1) = 10.132_rp           !Variable Nai
!            !        vcolo_exm(4,1) = 138.52_rp           !Variable Ki
!            !        vcolo_exm(5,1) = 0.00042_rp          !Variable CaSS
!            !        vcolo_exm(9,1) = 0.8978_rp           !Variable Rprime
!            !        !call exm_oceihe
!            !            !end do
!            !        return
!         end if
!
!      else if (kfl_hfmodmate_exm(mat) == 1) then
!         if ((moneclmate_exm(2, n) == 857) .and. (kfl_drugsmate_exm(n) == 0)) then   !Heart Failure at 70 bpm
!
!            call exm_init_voltages(EXM_ONECEL_HF_70BPM, mat)
!            !return
!
!         else if ((moneclmate_exm(2, n) /= 857) .or. (kfl_drugsmate_exm(n) == 1)) then
!
!            call runend("EXM_ONECEL: exm_oceihe currently is not working correctly!")
!
!            !elmlo_exm(1) = -83.4949_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
!            !elmlo_exm(2) = -83.4949_rp
!            !do nbeat = 1, beat
!            !ituss_exm = 3    !epicardial
!            !EPICARDIAL HEART FAILURE
!            !VOLTAGE -82.5277658119425
!            !vaulo_exm(1,1) = 0.00303276091526701_rp
!            !vaulo_exm(2,1) = 0.661659742738716_rp
!            !vaulo_exm(3,1) = 0.634628244886488_rp
!            !vaulo_exm(4,1) = 4.83560182302662e-05_rp
!            !vaulo_exm(5,1) = 0.926032127446219_rp
!            !vaulo_exm(12,1) = 0.999241154724418_rp
!            !vaulo_exm(6,1) = 0.999869780216330_rp
!            !vaulo_exm(7,1) = 3.79981651056992e-08_rp
!            !vaulo_exm(8,1) = 0.999996272388211_rp
!            !vaulo_exm(9,1) = 0.00455055453605926_rp
!            !vaulo_exm(10,1) = 0.00168040788756380_rp
!            !vaulo_exm(11,1) = 0.443217253596480_rp
!            !vcolo_exm(1,1) = 8.63303349957043e-05_rp
!            !vcolo_exm(2,1) = 2.97266869543131_rp
!            !vcolo_exm(3,1) = 9.72525639401306_rp
!            !vcolo_exm(4,1) = 135.848336580048_rp
!            !vcolo_exm(5,1) = 0.000257031579830260_rp
!            !vcolo_exm(9,1) = 0.980074819847184_rp
!            !call exm_oceihe
!
!            !do nbeat = 1, beat
!            !    ituss_exm = 1    !endocardial
!            !             !!!  ENDOCARDIAL  HEART FAILURE
!            !     !  VOLTAGE -84.2920820453902
!            !    vaulo_exm(1,1) = 0.00209247866395593_rp
!            !    vaulo_exm(2,1) = 0.717586845558180_rp
!            !    vaulo_exm(3,1) = 0.703610324031878_rp
!            !    vaulo_exm(4,1) = 3.82182239461387e-05_rp
!            !    vaulo_exm(5,1) = 0.927116903038717_rp
!            !    vaulo_exm(12,1) = 0.997941159192524_rpOTHERBPM
!            !    vaulo_exm(6,1) = 0.999857972863601_rp
!            !    vaulo_exm(7,1) = 2.82975518708235e-08_rp
!            !    vaulo_exm(8,1) = 0.499007377734646_rp
!            !    vaulo_exm(9,1) = 0.00404264224963596_rp
!            !    vaulo_exm(10,1) = 0.000568838457062049_rp
!            !    vaulo_exm(11,1) = 0.461434647116189_rp
!            !    vcolo_exm(1,1) = 0.000132911621580043_rp
!            !    vcolo_exm(2,1) = 3.97092446663378_rp
!            !    vcolo_exm(3,1) = 10.1377699148182_rp
!            !    vcolo_exm(4,1) = 138.592632537711_rp
!            !    vcolo_exm(5,1) = 0.000283334727151014_rp
!            !    vcolo_exm(9,1) = 0.979004485438106_rp
!            !    !call exm_oceihe
!
!            ! !do nbeat = 1, beat
!            !    ituss_exm = 2    !midmyocardial
!            !              !mid HEART FAILURE
!            !     ! Voltage:  -83.6649251312024
!            !    vaulo_exm(1,1) = 0.00238853756757240_rp           !Variable m
!            !    vaulo_exm(2,1) = 0.698129858756220_rp        !Variables h
!            !    vaulo_exm(3,1) = 0.670007606728020_rp        !Variable j
!            !    vaulo_exm(4,1) = 4.15544209009521e-05_rp           !Variable d
!            !    vaulo_exm(5,1) = 0.882815089054957_rp         !Variable f
!            !    vaulo_exm(12,1) = 0.995047915633200_rp          !Variable f2
!            !    vaulo_exm(6,1) = 0.999809858618127_rp         ! Variable fCa
!            !    vaulo_exm(7,1) = 3.14471296989901e-08_rp           !Variable r
!            !    vaulo_exm(8,1) = 0.999997027982031_rp           !Variable s
!            !    vaulo_exm(9,1) = 0.00539523790336693_rp        !Variable xs
!            !    vaulo_exm(10,1) = 0.00221286069613379_rp        !Variable xr1
!            !    vaulo_exm(11,1) = 0.454935810543141_rp       !Variable xr2
!            !    vcolo_exm(1,1) = 0.000148758208418421_rp                 !Variable Cai
!            !    vcolo_exm(2,1) = 4.50025313681586_rp                 !Variable CaSR
!            !    vcolo_exm(3,1) = 9.73203055605199_rp                 !Variable Nai
!            !    vcolo_exm(4,1) = 138.900841456135_rp                 !Variable Ki
!            !    vcolo_exm(5,1) = 0.000316819832530177_rp                 !Variable CaSS
!            !    vcolo_exm(9,1) = 0.976140456252427_rp        !Variable Rprime
!            !    !call exm_oceihe
!            !    return
!         end if
!      end if
!
!   end subroutine exm_onecel
!
!   subroutine exm_init_voltages(condition_type, mat)
!      !------------------------------------------------------------------------
!      !> @addtogroup Exmedi
!      !> @{
!      !> @file    exm_onecel.f90
!      !> @date    03/10/2012
!      !> @author  Mariano Vazquez
!      !> @brief   One cell simulations to obtain initial conditions
!      !> @details Runs for changes of cycle length and drug administration
!      !> @}
!      !------------------------------------------------------------------------
!
!      implicit none
!
!      integer(ip), intent(in) :: condition_type, mat
!
!      select case (condition_type)
!      case (EXM_ONECEL_NORMAL_70BPM)
!
!         vminimate_exm(EXM_CELLTYPE_EPI, mat) = -85.23_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL    !EPBench =-85.23  [mV]  <-
!         vminimate_exm(EXM_CELLTYPE_ENDO, mat) = -85.23_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL    !EPBench =-85.23  [mV]  <-
!         vminimate_exm(EXM_CELLTYPE_MID, mat) = -85.23_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL    !EPBench =-85.23  [mV]  <-
!         !epicardial
!         vauxi_exm_initial(1,  EXM_CELLTYPE_ENDO, mat) = 0.00172_rp       !Variable m
!         vauxi_exm_initial(2,  EXM_CELLTYPE_ENDO, mat) = 0.7444_rp        !Variable h                           !EPBench = 0.7444_rp
!         vauxi_exm_initial(3,  EXM_CELLTYPE_ENDO, mat) = 0.7045_rp        !Variable j                           !EPBench = 0.7045_rp
!         vauxi_exm_initial(4,  EXM_CELLTYPE_ENDO, mat) = 0.00003373_rp    !Variable d                           !EPBench = 0.00003373_rp
!         vauxi_exm_initial(5,  EXM_CELLTYPE_ENDO, mat) = 0.7888_rp        !Variable f                           !EPBench = 0.7888_rp
!         vauxi_exm_initial(12, EXM_CELLTYPE_ENDO, mat) = 0.9755_rp       !Variable f2                          !EPBench = 0.9755_rp
!         vauxi_exm_initial(6,  EXM_CELLTYPE_ENDO, mat) = 0.9953_rp        !Variable fCa                         !EPBench = 0.9953_rp
!         vauxi_exm_initial(7,  EXM_CELLTYPE_ENDO, mat) = 0.0000000242_rp  !Variable r                           !EPBench = 0.0000000242_rp
!         vauxi_exm_initial(8,  EXM_CELLTYPE_ENDO, mat) = 0.999998_rp      !Variable s                           !EPBench = 0.999998_rp
!         vauxi_exm_initial(9,  EXM_CELLTYPE_ENDO, mat) = 0.0095_rp        !Variable xs                          !EPBench = 0.0095_rp
!         vauxi_exm_initial(10, EXM_CELLTYPE_ENDO, mat) = 0.00621_rp      !Variable xr1                         !EPBench = 0.00621_rp
!         vauxi_exm_initial(11, EXM_CELLTYPE_ENDO, mat) = 0.4712_rp       !Variable xr2                         !EPBench = 0.4712_rp
!         vconc_initial(1,  EXM_CELLTYPE_ENDO, mat) = 0.000126_rp      !Variable Cai                         !EPBench = 0.000126_rp [mM]
!         vconc_initial(2,  EXM_CELLTYPE_ENDO, mat) = 3.64_rp          !Variable CaSR                        !EPBench = 3.64_rp     [mM]
!         vconc_initial(3,  EXM_CELLTYPE_ENDO, mat) = 8.608_rp         !Variable Nai                         !EPBench = 8.604_rp    [mM] <-
!         vconc_initial(4,  EXM_CELLTYPE_ENDO, mat) = 136.89_rp        !Variable Ki                          !EPBench = 136.89_rp   [mM]
!         vconc_initial(5,  EXM_CELLTYPE_ENDO, mat) = 0.00036_rp       !Variable CaSSi                       !EPBench = 0.0036_rpi  [mM]
!         vconc_initial(9,  EXM_CELLTYPE_ENDO, mat) = 0.9073_rp        !Variable Rprime                      !EPBench = 0.9073_rp
!         !endocardial
!         vauxi_exm_initial(1,  EXM_CELLTYPE_MID, mat) = 0.0016_rp        !Variable m
!         vauxi_exm_initial(2,  EXM_CELLTYPE_MID, mat) = 0.7573_rp        !Variables h
!         vauxi_exm_initial(3,  EXM_CELLTYPE_MID, mat) = 0.7225_rp        !Variable j
!         vauxi_exm_initial(4,  EXM_CELLTYPE_MID, mat) = 0.00003164_rp    !Variable d
!         vauxi_exm_initial(5,  EXM_CELLTYPE_MID, mat) = 0.8009_rp        !Variable f
!         vauxi_exm_initial(12, EXM_CELLTYPE_MID, mat) = 0.9778_rp       !Variable f2
!         vauxi_exm_initial(6,  EXM_CELLTYPE_MID, mat) = 0.99530_rp       !Variable fCa
!         vauxi_exm_initial(7,  EXM_CELLTYPE_MID, mat) = 0.00000002235_rp !Variable r
!         vauxi_exm_initial(8,  EXM_CELLTYPE_MID, mat) = 0.3212_rp        !Variable s
!         vauxi_exm_initial(9,  EXM_CELLTYPE_MID, mat) = 0.0087_rp        !Variable xs
!         vauxi_exm_initial(10, EXM_CELLTYPE_MID, mat) = 0.0045_rp       !Variable xr1
!         vauxi_exm_initial(11, EXM_CELLTYPE_MID, mat) = 0.4760_rp       !Variable xr2
!         vconc_initial(1,  EXM_CELLTYPE_MID, mat) = 0.00013_rp       !Variable Cai
!         vconc_initial(2,  EXM_CELLTYPE_MID, mat) = 3.7150_rp        !Variable CaSR
!         vconc_initial(3,  EXM_CELLTYPE_MID, mat) = 10.355_rp        !Variable Nai
!         vconc_initial(4,  EXM_CELLTYPE_MID, mat) = 138.4_rp         !Variable Ki
!         vconc_initial(5,  EXM_CELLTYPE_MID, mat) = 0.00036_rp       !Variable CaSS
!         vconc_initial(9,  EXM_CELLTYPE_MID, mat) = 0.9068_rp        !Variable Rprime
!         !midmyocardial
!         vauxi_exm_initial(1,  EXM_CELLTYPE_EPI, mat) = 0.0017_rp        !Variable m
!         vauxi_exm_initial(2,  EXM_CELLTYPE_EPI, mat) = 0.749_rp         !Variables h
!         vauxi_exm_initial(3,  EXM_CELLTYPE_EPI, mat) = 0.6788_rp        !Variable j
!         vauxi_exm_initial(4,  EXM_CELLTYPE_EPI, mat) = 0.00003288_rp    !Variable d
!         vauxi_exm_initial(5,  EXM_CELLTYPE_EPI, mat) = 0.7026_rp        !Variable f
!         vauxi_exm_initial(12, EXM_CELLTYPE_EPI, mat) = 0.9526_rp       !Variable f2
!         vauxi_exm_initial(6,  EXM_CELLTYPE_EPI, mat) = 0.99420_rp       ! Variable fCa
!         vauxi_exm_initial(7,  EXM_CELLTYPE_EPI, mat) = 0.00000002347_rp !Variable r
!         vauxi_exm_initial(8,  EXM_CELLTYPE_EPI, mat) = 1.0_rp           !Variable s
!         vauxi_exm_initial(9,  EXM_CELLTYPE_EPI, mat) = 0.0174_rp        !Variable xs
!         vauxi_exm_initial(10, EXM_CELLTYPE_EPI, mat) = 0.0165_rp       !Variable xr1
!         vauxi_exm_initial(11, EXM_CELLTYPE_EPI, mat) = 0.473_rp        !Variable xr2
!         vconc_initial(1,  EXM_CELLTYPE_EPI, mat) = 0.000153_rp      !Variable Cai
!         vconc_initial(2,  EXM_CELLTYPE_EPI, mat) = 4.2720_rp        !Variable CaSR
!         vconc_initial(3,  EXM_CELLTYPE_EPI, mat) = 10.132_rp        !Variable Nai
!         vconc_initial(4,  EXM_CELLTYPE_EPI, mat) = 138.52_rp        !Variable Ki
!         vconc_initial(5,  EXM_CELLTYPE_EPI, mat) = 0.00042_rp       !Variable CaSS
!         vconc_initial(9,  EXM_CELLTYPE_EPI, mat) = 0.8978_rp        !Variable Rprime
!
!      case (EXM_ONECEL_HF_70BPM)
!
!         vminimate_exm(EXM_CELLTYPE_ENDO, mat) = -83.4949_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
!         vminimate_exm(EXM_CELLTYPE_MID, mat) = -85.23_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL    !EPBench =-85.23  [mV]  <-
!         vminimate_exm(EXM_CELLTYPE_EPI, mat) = -85.23_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL    !EPBench =-85.23  [mV]  <-
!
!         !EPICARDIAL HEART FAILURE
!         !VOLTAGE -82.5277658119425
!         vauxi_exm_initial(1,  EXM_CELLTYPE_ENDO, mat) = 0.00303276091526701_rp
!         vauxi_exm_initial(2,  EXM_CELLTYPE_ENDO, mat) = 0.661659742738716_rp
!         vauxi_exm_initial(3,  EXM_CELLTYPE_ENDO, mat) = 0.634628244886488_rp
!         vauxi_exm_initial(4,  EXM_CELLTYPE_ENDO, mat) = 4.83560182302662e-05_rp
!         vauxi_exm_initial(5,  EXM_CELLTYPE_ENDO, mat) = 0.926032127446219_rp
!         vauxi_exm_initial(12, EXM_CELLTYPE_ENDO, mat) = 0.999241154724418_rp
!         vauxi_exm_initial(6,  EXM_CELLTYPE_ENDO, mat) = 0.999869780216330_rp
!         vauxi_exm_initial(7,  EXM_CELLTYPE_ENDO, mat) = 3.79981651056992e-08_rp
!         vauxi_exm_initial(8,  EXM_CELLTYPE_ENDO, mat) = 0.999996272388211_rp
!         vauxi_exm_initial(9,  EXM_CELLTYPE_ENDO, mat) = 0.00455055453605926_rp
!         vauxi_exm_initial(10, EXM_CELLTYPE_ENDO, mat) = 0.00168040788756380_rp
!         vauxi_exm_initial(11, EXM_CELLTYPE_ENDO, mat) = 0.443217253596480_rp
!         vconc_initial(1,  EXM_CELLTYPE_ENDO, mat) = 8.63303349957043e-05_rp
!         vconc_initial(2,  EXM_CELLTYPE_ENDO, mat) = 2.97266869543131_rp
!         vconc_initial(3,  EXM_CELLTYPE_ENDO, mat) = 9.72525639401306_rp
!         vconc_initial(4,  EXM_CELLTYPE_ENDO, mat) = 135.848336580048_rp
!         vconc_initial(5,  EXM_CELLTYPE_ENDO, mat) = 0.000257031579830260_rp
!         vconc_initial(9,  EXM_CELLTYPE_ENDO, mat) = 0.980074819847184_rp
!
!                !!!  ENDOCARDIAL  HEART FAILURE
!         !  VOLTAGE -84.2920820453902
!         vauxi_exm_initial(1,  EXM_CELLTYPE_MID, mat) = 0.00209247866395593_rp
!         vauxi_exm_initial(2,  EXM_CELLTYPE_MID, mat) = 0.717586845558180_rp
!         vauxi_exm_initial(3,  EXM_CELLTYPE_MID, mat) = 0.703610324031878_rp
!         vauxi_exm_initial(4,  EXM_CELLTYPE_MID, mat) = 3.82182239461387e-05_rp
!         vauxi_exm_initial(5,  EXM_CELLTYPE_MID, mat) = 0.927116903038717_rp
!         vauxi_exm_initial(12, EXM_CELLTYPE_MID, mat) = 0.997941159192524_rp
!         vauxi_exm_initial(6,  EXM_CELLTYPE_MID, mat) = 0.999857972863601_rp
!         vauxi_exm_initial(7,  EXM_CELLTYPE_MID, mat) = 2.82975518708235e-08_rp
!         vauxi_exm_initial(8,  EXM_CELLTYPE_MID, mat) = 0.499007377734646_rp
!         vauxi_exm_initial(9,  EXM_CELLTYPE_MID, mat) = 0.00404264224963596_rp
!         vauxi_exm_initial(10, EXM_CELLTYPE_MID, mat) = 0.000568838457062049_rp
!         vauxi_exm_initial(11, EXM_CELLTYPE_MID, mat) = 0.461434647116189_rp
!         vconc_initial(1,  EXM_CELLTYPE_MID, mat) = 0.000132911621580043_rp
!         vconc_initial(2,  EXM_CELLTYPE_MID, mat) = 3.97092446663378_rp
!         vconc_initial(3,  EXM_CELLTYPE_MID, mat) = 10.1377699148182_rp
!         vconc_initial(4,  EXM_CELLTYPE_MID, mat) = 138.592632537711_rp
!         vconc_initial(5,  EXM_CELLTYPE_MID, mat) = 0.000283334727151014_rp
!         vconc_initial(9,  EXM_CELLTYPE_MID, mat) = 0.979004485438106_rp
!
!         !mid HEART FAILUREEXM_ONECEL_HF_70BPM 0.00238853756757240_rp           !Variable m
!         vauxi_exm_initial(2,  EXM_CELLTYPE_EPI, mat) = 0.698129858756220_rp        !Variables h
!         vauxi_exm_initial(3,  EXM_CELLTYPE_EPI, mat) = 0.670007606728020_rp        !Variable j
!         vauxi_exm_initial(4,  EXM_CELLTYPE_EPI, mat) = 4.15544209009521e-05_rp           !Variable d
!         vauxi_exm_initial(5,  EXM_CELLTYPE_EPI, mat) = 0.882815089054957_rp         !Variable f
!         vauxi_exm_initial(12, EXM_CELLTYPE_EPI, mat) = 0.995047915633200_rp          !Variable f2
!         vauxi_exm_initial(6,  EXM_CELLTYPE_EPI, mat) = 0.999809858618127_rp         ! Variable fCa
!         vauxi_exm_initial(7,  EXM_CELLTYPE_EPI, mat) = 3.14471296989901e-08_rp           !Variable r
!         vauxi_exm_initial(8,  EXM_CELLTYPE_EPI, mat) = 0.999997027982031_rp           !Variable s
!         vauxi_exm_initial(9,  EXM_CELLTYPE_EPI, mat) = 0.00539523790336693_rp        !Variable xs
!         vauxi_exm_initial(10, EXM_CELLTYPE_EPI, mat) = 0.00221286069613379_rp        !Variable xr1
!         vauxi_exm_initial(11, EXM_CELLTYPE_EPI, mat) = 0.454935810543141_rp       !Variable xr2
!         vconc_initial(1,  EXM_CELLTYPE_EPI, mat) = 0.000148758208418421_rp                 !Variable Cai
!         vconc_initial(2,  EXM_CELLTYPE_EPI, mat) = 4.50025313681586_rp                 !Variable CaSR
!         vconc_initial(3,  EXM_CELLTYPE_EPI, mat) = 9.73203055605199_rp                 !Variable Nai
!         vconc_initial(4,  EXM_CELLTYPE_EPI, mat) = 138.900841456135_rp                 !Variable Ki
!         vconc_initial(5,  EXM_CELLTYPE_EPI, mat) = 0.000316819832530177_rp                 !Variable CaSS
!         vconc_initial(9,  EXM_CELLTYPE_EPI, mat) = 0.976140456252427_rp        !Variable Rprime
!
!      case default
!         call runend("EXM_INIT_VOLATGES: Undefined condition type.")
!      end select
!
!   end subroutine exm_init_voltages
!
!end module mod_exm_onecel
!