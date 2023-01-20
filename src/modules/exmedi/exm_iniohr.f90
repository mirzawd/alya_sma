!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_iniohr.f90
!> @author  Jazmin Aguado-Sierra
!> @brief   Initial condition setup for Ohara-Rudy heterogeneous model
!> @details Obtains initial conditions of Normal or Heart Failure cell at  70 bpm or 857 ms \n
!!   C++ Implementation of the O'Hara-Rudy dynamic (ORd) model for the \n
!!   undiseased human ventricular action potential and calcium transient \n
!!  \n
!!   The ORd model is described in the article "Simulation of the Undiseased \n
!!   Human Cardiac Ventricular Action Potential: Model Formulation and \n
!!   Experimental Validation" \n
!!   by Thomas O'Hara, Laszlo Virag, Andras Varro, and Yoram Rudy \n
!!   \n
!!   The article and supplemental materails are freely available in the \n
!!   Open Access jounal PLoS Computational Biology \n
!!   Link to Article: \n
!!   http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1002061 \n
!!    \n
!!   Email: tom.ohara@gmail.com / rudy@wustl.edu \n
!!   Web: http://rudylab.wustl.edu \n
!> @} 
!!-----------------------------------------------------------------------
subroutine exm_iniohr(ipoin,mat)
!subroutine exm_iniohr(kmodel_ipoin,ipoin,mat)

   use      def_master
   use      def_domain
   use      def_elmtyp
   use      def_exmedi
   use      def_kermod,           only :  kfl_celltype_fun
   use      mod_exm_cellmodel

   implicit none
   integer(ip), intent(in) :: ipoin, mat
   integer(ip)   :: ituss_exm
   !ituss = 0_ip  !!%endo = 1, epi = 0, M = 2  

   if(kfl_celltype_fun == 0_ip) then
      call runend("exm_iniohr: Cell type field is required")
   end if 

   ituss_exm = int(celty_exm(ipoin), kind=ip)

   if(ituss_exm == EXM_CELLTYPE_EPI) then   !epicardial
      elmag(ipoin,1:3) = vminimate_exm(ituss_exm,mat)     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
      vicel_exm(1:26,ipoin) = 0.0_rp
      vconc(1:11,ipoin,1) = vconc_initial(1:11,ituss_exm,mat) !8.17202753912090_rp      ! nai=2
      vconc(1:11,ipoin,2) = vconc_initial(1:11,ituss_exm,mat) !8.17202753912090_rp      ! nai=2
      vconc(1:11,ipoin,3) = vconc_initial(1:11,ituss_exm,mat) !8.17202753912090_rp      ! nai=2
      vauxi_exm(1:29,ipoin,1) = vauxi_exm_initial(1:29,ituss_exm,mat) 
      vauxi_exm(1:29,ipoin,2) = vauxi_exm_initial(1:29,ituss_exm,mat) 
      vauxi_exm(1:29,ipoin,3) = vauxi_exm_initial(1:29,ituss_exm,mat) 

   else if(ituss_exm == EXM_CELLTYPE_ENDO) then   !endocardial
      elmag(ipoin,1:3) = vminimate_exm(ituss_exm,mat)     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
      vicel_exm(1:26,ipoin) = 0.0_rp
      vconc(1:11,ipoin,1) = vconc_initial(1:11,ituss_exm,mat)      ! nai=2
      vconc(1:11,ipoin,2) = vconc_initial(1:11,ituss_exm,mat)       ! nai=2
      vconc(1:11,ipoin,3) = vconc_initial(1:11,ituss_exm,mat)      ! nai=2
      vauxi_exm(1:29,ipoin,1) = vauxi_exm_initial(1:29,ituss_exm,mat) 
      vauxi_exm(1:29,ipoin,2) = vauxi_exm_initial(1:29,ituss_exm,mat) 
      vauxi_exm(1:29,ipoin,3) = vauxi_exm_initial(1:29,ituss_exm,mat) 

   else if(ituss_exm == EXM_CELLTYPE_MID) then   !MIDmyocardial
      elmag(ipoin,1:3) = vminimate_exm(ituss_exm,mat)     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
      vicel_exm(1:26,ipoin) = 0.0_rp
      vconc(1:11,ipoin,1) = vconc_initial(1:11,ituss_exm,mat)       ! nai=2
      vconc(1:11,ipoin,2) = vconc_initial(1:11,ituss_exm,mat)       ! nai=2
      vconc(1:11,ipoin,3) = vconc_initial(1:11,ituss_exm,mat)       ! nai=2
      vauxi_exm(1:29,ipoin,1) = vauxi_exm_initial(1:29,ituss_exm,mat) 
      vauxi_exm(1:29,ipoin,2) = vauxi_exm_initial(1:29,ituss_exm,mat) 
      vauxi_exm(1:29,ipoin,3) = vauxi_exm_initial(1:29,ituss_exm,mat) 

   end if
end subroutine exm_iniohr
