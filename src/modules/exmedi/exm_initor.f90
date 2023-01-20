!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_onetor.f90
!> @date    12/04/2013
!> @author  Zhinuo Jenny Wang
!> @brief   Initial condition setup for ToR-ORd 2019 heterogeneous model\n
!!          The init_toggle controls whether the initialisation is run as single cell steady-state\n
!!          or during a 2D/3D simulation\n
!> @}
!------------------------------------------------------------------------

subroutine exm_initor(ipoin, mat)

   use def_master
   use def_domain
   use def_elmtyp
   use def_exmedi 
   use def_kermod,           only :  kfl_celltype_fun
   use mod_exm_cellmodel

   
   implicit none 
   integer(ip), intent(in) :: ipoin, mat
   integer(ip) :: ituss_exm

   if(kfl_celltype_fun == 0_ip) then
      call runend("exm_initor: Cell type field is required")
   end if 
  
   
   ! Initialise single cell states
   vicel_exm(1:21,ipoin) = 0.0_rp ! [muA/muF]

   ! Initialise each node 
   ! Initialise to values after 30 minutes pacing at cycle length of 1000 ms  
   ituss_exm = int(celty_exm(ipoin), kind=ip)

   elmag(ipoin,1:3) = vminimate_exm(ituss_exm,mat)
   
   if(ituss_exm == EXM_CELLTYPE_EPI) then   !epicardial
      elmag(ipoin,1:3) = vminimate_exm(ituss_exm,mat)
      vconc(1:14, ipoin, 1) = vconc_initial(1:14,ituss_exm,mat)
      vconc(1:14, ipoin, 2) = vconc_initial(1:14,ituss_exm,mat)
      vconc(1:14, ipoin, 3) = vconc_initial(1:14,ituss_exm,mat)
      vauxi_exm(1:31, ipoin, 1) = vauxi_exm_initial(1:31,ituss_exm,mat) 
      vauxi_exm(1:31, ipoin, 2) = vauxi_exm_initial(1:31,ituss_exm,mat)
      vauxi_exm(1:31, ipoin, 3) = vauxi_exm_initial(1:31,ituss_exm,mat)
      
   else if (ituss_exm == EXM_CELLTYPE_MID) then ! mid-myocardial
      elmag(ipoin,1:3) = vminimate_exm(ituss_exm,mat)
      vconc(1:14, ipoin, 1) = vconc_initial(1:14,ituss_exm,mat)
      vconc(1:14, ipoin, 2) = vconc_initial(1:14,ituss_exm,mat)
      vconc(1:14, ipoin, 3) = vconc_initial(1:14,ituss_exm,mat)
      vauxi_exm(1:31, ipoin, 1) = vauxi_exm_initial(1:31,ituss_exm,mat)
      vauxi_exm(1:31, ipoin, 2) = vauxi_exm_initial(1:31,ituss_exm,mat)
      vauxi_exm(1:31, ipoin, 3) = vauxi_exm_initial(1:31,ituss_exm,mat)
   
   else if (ituss_exm == EXM_CELLTYPE_ENDO) then ! endocardial
      elmag(ipoin,1:3) = vminimate_exm(ituss_exm,mat)
      vconc(1:14, ipoin, 1) = vconc_initial(1:14,ituss_exm,mat)
      vconc(1:14, ipoin, 2) = vconc_initial(1:14,ituss_exm,mat)
      vconc(1:14, ipoin, 3) = vconc_initial(1:14,ituss_exm,mat)
      vauxi_exm(1:31, ipoin, 1) = vauxi_exm_initial(1:31,ituss_exm,mat)
      vauxi_exm(1:31, ipoin, 2) = vauxi_exm_initial(1:31,ituss_exm,mat)
      vauxi_exm(1:31, ipoin, 3) = vauxi_exm_initial(1:31,ituss_exm,mat)
      
   end if

end subroutine exm_initor

