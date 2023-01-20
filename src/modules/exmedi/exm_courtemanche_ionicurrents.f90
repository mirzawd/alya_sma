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
!> @date    15/NOV/2020
!> @details 
!> @
!!-----------------------------------------------------------------------

!!#############################
    !voltage: elmag
    
    !Ca_i:vconc(1,:)
    !Ca_rel:vconc(2,:)  
    !Ca_up:vconc(3,:)
    !K_i:vconc(4,:)
    !Na_i:vconc(5,:)  
      
    !Gate variables
    !u:vauxi_exm(1,:)  ;  v:vauxi_exm(2,:)  ;     !w:vauxi_exm(3,:)
    !d:vauxi_exm(4,:)  ;  f_Ca:vauxi_exm(5,:) ;   !f:vauxi_exm(6,:)
    !h:vauxi_exm(7,:)  ;  j:vauxi_exm(8,:)    ;   !m:vauxi_exm(9,:)
    !xr:vauxi_exm(10,:)  ;  xs:vauxi_exm(11,:)  ;  !oa:vauxi_exm(12,:)
    !oi:vauxi_exm(13,:)  ;  ua:vauxi_exm(14,:)  ;  !ui: vauxi_exm(15,:)
    
    !I_Na:vicel_exm(1,:)  ;  I_Cal:vicel_exm(2,:)  ;  !I_to:vicel_exm(3,:)
    !I_kur:vicel_exm(4,:)  ;  I_kr:vicel_exm(5,:)  ;  !I_ks:vicel_exm(6,:)
    !I_k1:vicel_exm(7,:)  ;  I_NaK:vicel_exm(8,:)  ;  !I_NaCa:vicel_exm(9,:)
    !I_bCa:vicel_exm(10,:)  ;  I_bNa:vicel_exm(11,:)  ;  !I_rel:vicel_exm(12,:)
    !I_tr:vicel_exm(13,:)  ;  I_up:vicel_exm(14,:)  ;  I_stim:vicel_exm(15,:)
!!#############################
!! main code

subroutine exm_courtemanche_ionicurrents(ipoin,xioni,dioni,cai)

  use      def_parame
  use      def_master
  use      def_elmtyp
  use      def_domain
  use      def_exmedi
  use      mod_eccoupling
  use      def_kermod
  use      mod_exm_drugs
  use      mod_exm_courtemanche_model

  ! definition of variables
  implicit none
  integer(ip), intent(in) :: ipoin !< node
  real(rp), intent(out) :: xioni   !< current
  real(rp), intent(out)   :: dioni !< current derivative
  real(rp), intent(out)   :: cai   !< current calcium
  real(rp) :: dtimeEP              !< Time step fol solving cell model (ms)
  integer(ip) :: i, imate, ituss_exm
!  integer(ip) :: n
  logical :: flag_land, flag_3D, flag_isac
  real(rp) :: statvar(7,2)
!  real(rp) :: kmcmdn, cmdnmax, trpnmax, sac


  if(INOTMASTER) then

     if (kfl_cellmod(nodemat(ipoin)) == EXMSLD_CELL_COURTE) then 

         ituss_exm = int(celty_exm(ipoin), kind=ip)  

         imate = nodemat(ipoin)
 
         dtimeEP=dtime * 1000.0_rp  ! in (ms)

         ! Provide stimulus
         vicel_exm(18,ipoin) = appfi_exm(ipoin)

         flag_land = .FALSE.
         if (kfl_exmsld_ecc) then
            if (kfl_eccty(imate) == EXMSLD_EMEC_LAND .or. kfl_eccty(imate) == EXMSLD_EMEC_LAND_BIDIR) flag_land = .TRUE. 
         end if
         flag_3D = .TRUE.
       
         ! if ( exm_drugs_ondrugs(n) ) then
         !   drug_effects = exm_drugs_calculate_effect(n)
         !   ! define drug variables 
         !end if
   
         flag_isac = .FALSE.
         do imate = 1,nmate_exm
           if (kfl_isac_exm(imate)==1_ip) flag_isac = .TRUE.
         end do

         statvar(:,:) = 0.0_rp
 
         ! here call the funciton of the model
         call exm_courtemanche_model(ituss_exm,elmag(ipoin,:),vconc(:,ipoin,:),vauxi_exm(:,ipoin,:),vicel_exm(:,ipoin),dtimeEP, &
           & qneto_exm(ipoin),flag_land,flag_3D,flag_isac,statvar,ipoin)

        ! if (flag_land) then
        !   call exm_courteland_calcium(kmcmdn,cmdnmax,trpnmax,dtimeEP,ipoin,cai,sac,elmag(ipoin,ITER_K),ituss_exm,imate)
        !   vconc(1,ipoin,1) = cai
        ! end if

        !----------------------------------------------------------------------------------
        ! Outputs
        !-----------------------------------------------------------------------------------
        dioni = 0.0_rp
        xioni = 0.0_rp
        cai = vconc(1,ipoin,1)
        do i=1,17
          xioni = + vicel_exm(i,ipoin)
        end do

        ! Uptdate variables
        vconc(:,ipoin,1) = vconc(:,ipoin,2)
        vauxi_exm(:,ipoin,1) = vauxi_exm(:,ipoin,2)
        elmag(ipoin,1) = elmag(ipoin,2)

   end if
 
  end if

end subroutine exm_courtemanche_ionicurrents
