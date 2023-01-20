!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_torord_ionicurrents.f90
!> @author  Zhinuo Jenny Wang
!> @brief   Wrapper for ToR-ORd 2020 model in 3D
!> @date    14/FEB/2020
!> @details 
!> @}
!!-----------------------------------------------------------------------
subroutine exm_torord_ionicurrents(ipoin, xioni, dioni, cai)

   use      def_parame
   use      def_master
   use      def_elmtyp
   use      def_domain
   use      def_exmedi
   use      mod_exm_torord_model, only: exm_torord_model
   use      def_kermod,           only : kfl_celltype_fun
   use      mod_eccoupling
   use      mod_exm_drugs

   ! Definition of variables
   implicit none
   integer(ip), intent(in) :: ipoin !< node
   real(rp), intent(out) :: xioni   !< current
   real(rp), intent(out)   :: dioni !< current derivative
   real(rp), intent(out)  :: cai    !< Intracellular calcium concentration
   real(rp) :: a2bas                !< Apex-to-base scaling factor for Gks
   real(rp) :: dtimeEP              !< Time step for solving cell model
   real(rp) :: gkr_scaling          !< Possible GKr scaling for BZ
   real(rp) :: drugd(24_ip)             !< Drug effects scaling
   integer(ip) :: i, n,imate,ituss_exm
   logical :: flag_land,flag_3D,flag_border,flag_drug
   real(rp) :: statvar(7,2)

   if(kfl_celltype_fun == 0_ip) then !if there is no celltype field
      call runend("exm_torord_ionicurrents: Cell type filed is required")
   end if


   if (INOTMASTER) then

      ! Get apex to base gradient of gKs if it is defined
      if (modab_exm == 0_ip) then ! no gradient
         a2bas = 1.0_rp
      else                        ! gradient
         a2bas = atbhe_exm(1,ipoin)
      end if

      if (kfl_cellmod(nodemat(ipoin)) == EXMSLD_CELL_TORORD) then

         ituss_exm = int(celty_exm(ipoin), kind=ip)

         n = nodemat(ipoin)

         if (n==0) then
            call runend('Elements with material 0 are present')
         end if

         dtimeEP = dtime * 1000.0_rp ! Convert to [ms]

         ! Provide stimulus
         vicel_exm(20,ipoin) = appfi_exm(ipoin)
         flag_land = .FALSE.
         if (kfl_exmsld_ecc) then
            do imate = 1,nmate_exm
               if (kfl_eccty(imate) == 4_ip) flag_land=.TRUE.
            end do
         end if
         flag_3D = .TRUE.

         !if (kfl_borde_exm(n)==1) then
         !   flag_border = .TRUE.
         !   gkr_scaling = border_gkrsc_exm(n)
         !else
            flag_border = .FALSE.
            gkr_scaling = 1.0_rp
         !end if

         if ( exm_drugs_ondrugs(n) ) then
            flag_drug = .TRUE.
            !TODO: drugs need redoing here in this model, it's completely outdated
            call exm_drugs_get_drugdmate(n, drugd, size(drugd, 1_ip, kind=ip)/4_ip, 2_ip)
         else
            flag_drug = .FALSE.
            drugd = 1.0_rp
         end if

         statvar(:,:) = 0.0_rp

         call exm_torord_model(ituss_exm, elmag(ipoin,ITER_K),vconc(:,ipoin,:),vauxi_exm(:,ipoin,:),vicel_exm(:,ipoin),dtimeEP, a2bas, &
           & qneto_exm(ipoin),flag_land,flag_3D,flag_border,flag_drug,gkr_scaling,drugd,statvar,ipoin)

         dioni = 0.0_rp
         xioni = 0.0_rp
         cai = vconc(5,ipoin,2) ! Output intracellular calcium concentration
         do i=1,19
            xioni = xioni + vicel_exm(i,ipoin) ! Istim is added in mod_exm_ionicurrents
         end do
         xioni = xioni + vicel_exm(21,ipoin) ! Add stretch-activated current

         vconc(:,ipoin,1) = vconc(:,ipoin,2)
         vauxi_exm(:,ipoin,1) = vauxi_exm(:,ipoin,2)
      end if
   end if

end subroutine exm_torord_ionicurrents
