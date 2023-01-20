!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_vect_csm.f90
!> @author  Adria Quintanas-Corominas
!> @date    November 2021
!> @brief   Toolbox with solid mechanics operations for the vectorised assembly
!> @details Toolbox with solid mechanics operations for the vectorised assembly
!------------------------------------------------------------------

module mod_sld_vect_csm

#include "def_solidz_vect_dim.inc" 
   use def_kintyp, only  : ip, rp, lg
   use def_domain, only  : ndime
   use mod_sld_vect_maths

   implicit none

   real(rp),    parameter :: PI = 16.0_rp*ATAN(1.0_rp/5.0_rp) - 4.0_rp*ATAN(1.0_rp/239.0_rp)
   integer(ip), parameter :: Idd(3,3) = reshape( [ 1_ip, 0_ip, 0_ip  , &
                                                   0_ip, 1_ip, 0_ip  , &
                                                   0_ip, 0_ip, 1_ip ], &
                                                  [ 3, 3 ] )
   real(rp),    parameter :: Inn(3,3) = reshape( [ 1.0_rp, 0.0_rp, 0.0_rp  , &
                                                   0.0_rp, 1.0_rp, 0.0_rp  , &
                                                   0.0_rp, 0.0_rp, 1.0_rp ], & 
                                                  [ 3, 3 ] )
   contains


      subroutine sld_vcsm_get_det_tensor_at_GP(pgaus,gpgdi,gpdet)

         integer(ip), intent(in)    :: pgaus
         real(rp),    intent(in)    :: gpgdi(DVS,ZDIME,ZDIME,ZGAUS)
         real(rp),    intent(out)   :: gpdet(DVS,ZGAUS)
         integer(ip)                :: ig

         do ig = 1, ZGAUS
            gpdet(:,ig) = vmath_DET(ZDIME,gpgdi(:,:,:,ig))
         enddo

      end subroutine sld_vcsm_get_det_tensor_at_GP


      subroutine sld_vcsm_get_right_cauchy_strain_tensor_at_GP(pgaus,gpgdi,gpcau)

         integer(ip), intent(in)    :: pgaus
         real(rp),    intent(in)    :: gpgdi(DVS,ZDIME,ZDIME,ZGAUS)
         real(rp),    intent(out)   :: gpcau(DVS,ZDIME,ZDIME,ZGAUS)
         integer(ip)                :: ig

         do ig = 1, ZGAUS
            gpcau(:,:,:,ig) = TxM_ndime(gpgdi(:,:,:,ig),gpgdi(:,:,:,ig))
         enddo

      end subroutine sld_vcsm_get_right_cauchy_strain_tensor_at_GP


      subroutine sld_vcsm_get_lagrange_strain_tensors_at_GP(pgaus,gpgdi,gpcau,gpgre)

         integer(ip), intent(in)    :: pgaus
         real(rp),    intent(in)    :: gpgdi(DVS,ZDIME,ZDIME,ZGAUS)
         real(rp),    intent(out)   :: gpcau(DVS,ZDIME,ZDIME,ZGAUS)
         real(rp),    intent(out)   :: gpgre(DVS,ZDIME,ZDIME,ZGAUS)
         integer(ip)                :: ii, ig

         do ig = 1, ZGAUS 
            ! Right Cauchy strain tensor
            gpcau(:,:,:,ig) = TxM_ndime(gpgdi(:,:,:,ig),gpgdi(:,:,:,ig))
            ! Green-Lagrangen strain tensor
            gpgre(:,:,:,ig) = 0.5_rp * gpcau(:,:,:,ig)
            do ii = 1, ZDIME
               gpgre(1:DVS,ii,ii,ig) = gpgre(1:DVS,ii,ii,ig) - 0.5_rp
            enddo
         enddo

      end subroutine sld_vcsm_get_lagrange_strain_tensors_at_GP

      
      subroutine sld_vcsm_get_infinitesimal_strain_tensors_at_GP(pgaus,gpgdi,gpinf)

         integer(ip), intent(in)    :: pgaus
         real(rp),    intent(in)    :: gpgdi(DVS,ZDIME,ZDIME,ZGAUS)
         real(rp),    intent(out)   :: gpinf(DVS,ZDIME,ZDIME,ZGAUS)
         integer(ip)                :: ii, jj, ig
         
         do ig = 1, ZGAUS
            gpinf(1:DVS,:,:,ig) = 0.0_rp
            do ii = 1, ZDIME
               do jj = 1, ZDIME
                  gpinf(1:DVS,ii,jj,ig) = &
                       gpinf(1:DVS,ii,jj,ig) + 0.5_rp*(gpgdi(1:DVS,ii,jj,ig) + gpgdi(1:DVS,jj,ii,ig)) - Inn(ii,jj)
               enddo
            enddo
         enddo

      end subroutine sld_vcsm_get_infinitesimal_strain_tensors_at_GP

      !-----------------------------------------------------------------------
      !> 
      !> @author  aquintanas
      !> @date    2022-11-16
      !> @brief   From 2nd Piola-Kirchhoff to 1st Piola-Kirchhoff stress tensor
      !> @details 
      !>          PK1_{ij} = F_{ik} S_{kj}
      !>          P_{ij}   = S_{ik} F_{kj}^{T} (Nominal stress in Bely's book)
      !>
      !-----------------------------------------------------------------------

      subroutine sld_vcsm_transform_PK2_to_PK1_at_GP(pgaus,gpgdi,gpstr,gppio)

         integer(ip), intent(in)    :: pgaus
         real(rp),    intent(in)    :: gpgdi(DVS,ZDIME,ZDIME,ZGAUS) !< Deformation gradient tensor
         real(rp),    intent(in)    :: gpstr(DVS,ZDIME,ZDIME,ZGAUS) !< Second Piola-Kirchhoff tensor (symmetric)
         real(rp),    intent(out)   :: gppio(DVS,ZDIME,ZDIME,ZGAUS) !< First Piola-Kirchhoff tensor (non-symmetric)
         integer(ip)                :: ig

         do ig = 1, ZGAUS
            gppio(:,:,:,ig) = AxB_dime(gpgdi(:,:,:,ig),gpstr(:,:,:,ig))
         enddo

      end subroutine sld_vcsm_transform_PK2_to_PK1_at_GP

      !-----------------------------------------------------------------------
      !> 
      !> @author  aquintanas
      !> @date    2022-11-16
      !> @brief   From Cauchy to Nominal stress tensor
      !> @details 
      !>          PK1_{ij} = J Sigma_{ik} F_{kj}^{-T}
      !>          P_{ij}   = J F_{ik}^{-1} Sigma_{kj} (Nominal stress in Bely's book)
      !>
      !-----------------------------------------------------------------------
      
      subroutine sld_vcsm_transform_sigma_to_PK1_at_GP(pgaus,gpidg,gpdet,gpsig,gppio)

         integer(ip), intent(in)    :: pgaus
         real(rp),    intent(in)    :: gpidg(DVS,ZDIME,ZDIME,ZGAUS) !< Inverse of Deformation gradient tensor
         real(rp),    intent(in)    :: gpdet(DVS,ZGAUS)             !< Determinant of the deformantion gradient tensor
         real(rp),    intent(in)    :: gpsig(DVS,ZDIME,ZDIME,ZGAUS) !< Cauchy stress tensor
         real(rp),    intent(out)   :: gppio(DVS,ZDIME,ZDIME,ZGAUS) !< First Piola-Kirchhoff tensor
         integer(ip)                :: ig

         do ig = 1, ZGAUS
            gppio(:,:,:,ig) = SxAxB_dime(gpdet(:,ig),gpidg(:,:,:,ig),gpsig(:,:,:,ig))
         enddo
         
      end subroutine sld_vcsm_transform_sigma_to_PK1_at_GP
      
      subroutine sld_vcsm_transform_dPK2dE_to_dPK1dF_at_GP(pgaus,gpgdi,gpstr,gpdds,gptmo)

         integer(ip), intent(in)    :: pgaus
         real(rp),    intent(in)    :: gpgdi(DVS,ZDIME,ZDIME,ZGAUS)             !< Deformation gradient tensor
         real(rp),    intent(in)    :: gpstr(DVS,ZDIME,ZDIME,ZGAUS)             !< Second Piola-Kirchoff tensor
         real(rp),    intent(in)    :: gpdds(DVS,ZDIME,ZDIME,ZDIME,ZDIME,ZGAUS) !< Second elasticity tensor
         real(rp),    intent(out)   :: gptmo(DVS,ZDIME,ZDIME,ZDIME,ZDIME,ZGAUS) !< First elasticity tensor
         integer(ip)                :: ig, ii, jj, kk, ll, mm, nn
         real(rp),    parameter     :: one = 1.0_rp, zero = 0.0_rp
         real(rp),    parameter     :: tkron(3,3) = reshape( &
                                          (/ one, zero, zero, zero, one, zero, zero, zero, one /) , &
                                          (/ 3, 3 /) )

         do ig = 1, ZGAUS
            ! Geometric contribution
            do ii = 1, ZDIME
               do jj = 1, ZDIME
                  do kk = 1, ZDIME
                     do ll = 1, ZDIME
                        gptmo(1:DVS,ii,jj,kk,ll,ig) = tkron(ii,kk) * gpstr(1:DVS,jj,ll,ig)
                     enddo
                  enddo
               enddo
            enddo
            ! Material contribution
            do ii = 1, ZDIME
               do jj = 1, ZDIME
                  do kk = 1, ZDIME
                     do ll = 1, ZDIME
                        do mm = 1, ZDIME
                           do nn = 1, ZDIME
                              gptmo(1:DVS,ii,jj,kk,ll,ig) = gptmo(1:DVS,ii,jj,kk,ll,ig) + &
                                 gpdds(1:DVS,mm,jj,nn,ll,ig) * gpgdi(1:DVS,ii,mm,ig) * gpgdi(1:DVS,kk,nn,ig)
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo

      end subroutine sld_vcsm_transform_dPK2dE_to_dPK1dF_at_GP


      subroutine sld_vcsm_compact_stress_at_GP(pgaus,gpstr,gpsig)

         integer(ip), intent(in)    :: pgaus
         real(rp),    intent(in)    :: gpstr(DVS,ZDIME,ZDIME,ZGAUS) !< Stress tensor
         real(rp),    intent(out)   :: gpsig(DVS,ZVOIG,ZGAUS)       !< Stress tensor compacted
         integer(ip)                :: ig

         if( ZDIME == 2_ip )then
            do ig = 1, ZGAUS
               gpsig(:,:,ig) = vcsm_compact_stress_tensor_2D(gpstr(:,:,:,ig)) 
            enddo
         else
            do ig = 1, ZGAUS
               gpsig(:,:,ig) = vcsm_compact_stress_tensor_3D(gpstr(:,:,:,ig))
            enddo
         endif         

      end subroutine sld_vcsm_compact_stress_at_GP


      subroutine sld_vcsm_uncompact_stress_at_GP(pgaus,gpsig,gpstr)

         integer(ip), intent(in)    :: pgaus
         real(rp),    intent(in)    :: gpsig(DVS,ZVOIG,ZGAUS)       !< Stress tensor compacted
         real(rp),    intent(out)   :: gpstr(DVS,ZDIME,ZDIME,ZGAUS) !< Stress tensor
         integer(ip)                :: ig

         if( ZDIME == 2_ip )then
            do ig = 1, ZGAUS
               gpstr(:,:,:,ig) = vcsm_uncompact_stress_tensor_2D(gpsig(:,:,ig))
            enddo
         else
            do ig = 1, ZGAUS
               gpstr(:,:,:,ig) = vcsm_uncompact_stress_tensor_3D(gpsig(:,:,ig))
            enddo
         endif

      end subroutine sld_vcsm_uncompact_stress_at_GP


      subroutine sld_vcsm_compact_strain_at_GP(pgaus,gpgre,gpeps)

         integer(ip), intent(in)    :: pgaus
         real(rp),    intent(in)    :: gpgre(DVS,ZDIME,ZDIME,ZGAUS) !< Strain tensor
         real(rp),    intent(out)   :: gpeps(DVS,ZVOIG,ZGAUS)       !< Strain tensor compacted
         integer(ip)                :: ig

         if( ZDIME == 2_ip )then
            do ig = 1, ZGAUS
               gpeps(:,:,ig) = vcsm_compact_strain_tensor_2D(gpgre(:,:,:,ig))
            enddo
         else
            do ig = 1, ZGAUS
               gpeps(:,:,ig) = vcsm_compact_strain_tensor_3D(gpgre(:,:,:,ig))
            enddo
         endif

      end subroutine sld_vcsm_compact_strain_at_GP


      subroutine sld_vcsm_uncompact_strain_at_GP(pgaus,gpeps,gpgre)

         integer(ip), intent(in)    :: pgaus
         real(rp),    intent(in)    :: gpeps(DVS,ZVOIG,ZGAUS)       !< Strain tensor compacted
         real(rp),    intent(out)   :: gpgre(DVS,ZDIME,ZDIME,ZGAUS) !< Strain tensor
         integer(ip)                :: ig

         if( ZDIME == 2_ip )then
            do ig = 1, ZGAUS
               gpgre(:,:,:,ig) = vcsm_uncompact_strain_tensor_2D(gpeps(:,:,ig))
            enddo
         else
            do ig = 1, ZGAUS
               gpgre(:,:,:,ig) = vcsm_uncompact_strain_tensor_3D(gpeps(:,:,ig))
            enddo
         endif

      end subroutine sld_vcsm_uncompact_strain_at_GP


      subroutine sld_vcsm_compact_stiffness_tensor_at_GP(pgaus,ten,voi)

         integer(ip), intent(in)    :: pgaus
         real(rp),    intent(in)    :: ten(DVS,ZDIME,ZDIME,ZDIME,ZDIME,ZGAUS) !< Stiffness tensor
         real(rp),    intent(out)   :: voi(DVS,ZVOIG,ZVOIG,ZGAUS)             !< Stiffness tensor compacted
         integer(ip)                :: ig

         do ig = 1, ZGAUS
            voi(:,:,:,ig) = vcsm_compact_stiffness_tensor(ten(:,:,:,:,:,ig))
         enddo

      end subroutine sld_vcsm_compact_stiffness_tensor_at_GP


      subroutine sld_vcsm_uncompact_stiffness_tensor_at_GP(pgaus,voi,ten)

         integer(ip), intent(in)    :: pgaus
         real(rp),    intent(in)    :: voi(DVS,ZVOIG,ZVOIG,ZGAUS)             !< Stiffness tensor compacted
         real(rp),    intent(out)   :: ten(DVS,ZDIME,ZDIME,ZDIME,ZDIME,ZGAUS) !< Stiffness tensor
         integer(ip)                :: ig

         do ig = 1, ZGAUS
            ten(:,:,:,:,:,ig) = vcsm_uncompact_stiffness_tensor(voi(:,:,:,ig))
         enddo

      end subroutine sld_vcsm_uncompact_stiffness_tensor_at_GP


      subroutine sld_vcsm_compact_transformation_stress_at_GP(itask,pgaus,basis,gpORI,gpTRA)

         integer(ip), intent(in)    :: itask
         integer(ip), intent(in)    :: pgaus
         real(rp),    intent(in)    :: basis(DVS,ZDIME,ZDIME)
         real(rp),    intent(in)    :: gpORI(DVS,ZVOIG,ZGAUS) !< Original compacted strain tensor 
         real(rp),    intent(out)   :: gpTRA(DVS,ZVOIG,ZGAUS) !< Transfor compacted strain tensor
         integer(ip)                :: ig
         real(rp)                   :: voiba(DVS,ZVOIG,ZVOIG)
         integer(ip), parameter     :: GLOB_to_MAT = 0_ip, MAT_to_GLOB = 1_ip

         if( ZDIME == 2_ip )then
            !call runend('sld_vcsm_compact_rotation_matrix_at_GP: 2D not programmed')
         else
            if(     itask == MAT_to_GLOB )then
               voiba = vcsm_compact_transformation_matrix_Tg_3D(basis)
               do ig = 1, ZGAUS
                  gpTRA(:,:,ig) =  TxV_voigt(voiba(:,:,:),gpORI(:,:,ig))
               enddo
            elseif( itask == GLOB_to_MAT )then
               voiba = vcsm_compact_transformation_matrix_T_3D(basis)
               do ig = 1, ZGAUS
                  gpTRA(:,:,ig) =  MxV_voigt(voiba(:,:,:),gpORI(:,:,ig))
               enddo
            endif
         endif

      end subroutine sld_vcsm_compact_transformation_stress_at_GP


      subroutine sld_vcsm_compact_transformation_strain_at_GP(itask,pgaus,basis,gpORI,gpTRA)

         integer(ip), intent(in)    :: itask
         integer(ip), intent(in)    :: pgaus
         real(rp),    intent(in)    :: basis(DVS,ZDIME,ZDIME)
         real(rp),    intent(in)    :: gpORI(DVS,ZVOIG,ZGAUS) !< Original compacted strain tensor 
         real(rp),    intent(out)   :: gpTRA(DVS,ZVOIG,ZGAUS) !< Transfor compacted strain tensor
         integer(ip)                :: ig
         real(rp)                   :: voiba(DVS,ZVOIG,ZVOIG)
         integer(ip), parameter     :: GLOB_to_MAT = 0_ip, MAT_to_GLOB = 1_ip

         if( ZDIME == 2_ip )then
            !call runend('sld_vcsm_compact_rotation_matrix_at_GP: 2D not programmed')
         else
            if(     itask == MAT_to_GLOB )then
               voiba = vcsm_compact_transformation_matrix_T_3D(basis)
               do ig = 1, ZGAUS
                  gpTRA(:,:,ig) =  TxV_voigt(voiba(:,:,:),gpORI(:,:,ig)) 
               enddo
            elseif( itask == GLOB_to_MAT )then
               voiba = vcsm_compact_transformation_matrix_Tg_3D(basis)
               do ig = 1, ZGAUS
                  gpTRA(:,:,ig) =  MxV_voigt(voiba(:,:,:),gpORI(:,:,ig)) 
               enddo
            endif
         endif

      end subroutine sld_vcsm_compact_transformation_strain_at_GP


      subroutine sld_vcsm_compact_transformation_stiffness_at_GP(itask,pgaus,basis,gpORI,gpTRA)

         integer(ip), intent(in)    :: itask
         integer(ip), intent(in)    :: pgaus
         real(rp),    intent(in)    :: basis(DVS,ZDIME,ZDIME)
         real(rp),    intent(in)    :: gpORI(DVS,ZVOIG,ZVOIG,ZGAUS) !< Original compacted stiffness tensor 
         real(rp),    intent(out)   :: gpTRA(DVS,ZVOIG,ZVOIG,ZGAUS) !< Transfor compacted stiffness tensor
         integer(ip)                :: ig, ii, jj, kk, ll
         real(rp)                   :: voiba(DVS,ZVOIG,ZVOIG)
         integer(ip), parameter     :: GLOB_to_MAT = 0_ip, MAT_to_GLOB = 1_ip

         if( ZDIME == 2_ip )then
            !call runend('sld_vcsm_compact_rotation_matrix_at_GP: 2D not programmed')
         else
            if(     itask == GLOB_to_MAT )then
               voiba = vcsm_compact_transformation_matrix_T_3D(basis)
               do ig = 1, ZGAUS 
                  gpTRA(1:DVS,:,:,ig) = 0.0_rp 
                  do ii = 1, 6; do jj = 1, 6; do kk = 1, 6; do ll = 1, 6
                     gpTRA(1:DVS,ii,jj,ig) = gpTRA(1:DVS,ii,jj,ig) + voiba(1:DVS,ii,kk) * voiba(1:DVS,jj,ll) * &
                        gpORI(1:DVS,kk,ll,ig)
                 enddo; enddo; enddo; enddo
               enddo
            elseif( itask == MAT_to_GLOB )then
               voiba = vcsm_compact_transformation_matrix_Tg_3D(basis)
               do ig = 1, ZGAUS
                  gpTRA(1:DVS,:,:,ig) = 0.0_rp 
                  do ii = 1, 6; do jj = 1, 6; do kk = 1, 6; do ll = 1, 6
                     gpTRA(1:DVS,ii,jj,ig) = gpTRA(1:DVS,ii,jj,ig) + voiba(1:DVS,kk,ii) * voiba(1:DVS,ll,jj) * &
                        gpORI(1:DVS,kk,ll,ig)
                  enddo; enddo; enddo; enddo
               enddo
            endif
         endif

      end subroutine sld_vcsm_compact_transformation_stiffness_at_GP

      !
      ! VECTORISED CONTINUUM SOLID MECHANICS FUNCTIONS
      !
      pure function vcsm_compact_stress_tensor_2D(ten) result(voi)

         real(rp),    intent(in)    :: ten(DVS,ZDIME,ZDIME) !< Stress tensor compacted
         real(rp)                   :: voi(DVS,ZVOIG)       !< Stress tensor

         voi(1:DVS,1) = ten(1:DVS,1,1) 
         voi(1:DVS,2) = ten(1:DVS,2,2) 
         voi(1:DVS,3) = ten(1:DVS,1,2) 

      end function vcsm_compact_stress_tensor_2D


      pure function vcsm_compact_stress_tensor_3D(ten) result(voi)

         real(rp),    intent(in)    :: ten(DVS,ZDIME,ZDIME) !< Stress tensor compacted
         real(rp)                   :: voi(DVS,ZVOIG)       !< Stress tensor

         voi(1:DVS,1) = ten(1:DVS,1,1)  
         voi(1:DVS,2) = ten(1:DVS,2,2)  
         voi(1:DVS,3) = ten(1:DVS,3,3)  
         voi(1:DVS,4) = ten(1:DVS,2,3)  
         voi(1:DVS,5) = ten(1:DVS,1,3)  
         voi(1:DVS,6) = ten(1:DVS,1,2)  

      end function vcsm_compact_stress_tensor_3D


      pure function vcsm_compact_strain_tensor_2D(ten) result(voi)

         real(rp),    intent(in)    :: ten(DVS,ZDIME,ZDIME) !< Strain tensor compacted
         real(rp)                   :: voi(DVS,ZVOIG)       !< Strain tensor

         voi(1:DVS,1) = ten(1:DVS,1,1)  
         voi(1:DVS,2) = ten(1:DVS,2,2)  
         voi(1:DVS,3) = ten(1:DVS,1,2) + ten(1:DVS,2,1)  

      end function vcsm_compact_strain_tensor_2D


      pure function vcsm_compact_strain_tensor_3D(ten) result(voi)
         implicit none
         real(rp),    intent(in)    :: ten(DVS,ZDIME,ZDIME) !< Strain tensor compacted
         real(rp)                   :: voi(DVS,ZVOIG)       !< Strain tensor

         voi(1:DVS,1) = ten(1:DVS,1,1) 
         voi(1:DVS,2) = ten(1:DVS,2,2)
         voi(1:DVS,3) = ten(1:DVS,3,3)
         voi(1:DVS,4) = ten(1:DVS,2,3) + ten(1:DVS,3,2)
         voi(1:DVS,5) = ten(1:DVS,1,3) + ten(1:DVS,3,1)
         voi(1:DVS,6) = ten(1:DVS,1,2) + ten(1:DVS,2,1)

      end function vcsm_compact_strain_tensor_3D


      pure function vcsm_uncompact_stress_tensor_2D(voi) result(ten)
         implicit none 
         real(rp),    intent(in)    :: voi(DVS,ZVOIG)       !< Strain tensor compacted
         real(rp)                   :: ten(DVS,ZDIME,ZDIME) !< Strain tensor

         ten(1:DVS,1,1) = voi(1:DVS,1)
         ten(1:DVS,2,2) = voi(1:DVS,2)
         ten(1:DVS,1,2) = voi(1:DVS,3)
         ten(1:DVS,2,1) = voi(1:DVS,3)

      end function vcsm_uncompact_stress_tensor_2D

  
      pure function vcsm_uncompact_stress_tensor_3D(voi) result(ten)

         real(rp),    intent(in)    :: voi(DVS,ZVOIG)       !< Strain tensor compacted
         real(rp)                   :: ten(DVS,ZDIME,ZDIME) !< Strain tensor

         ten(1:DVS,1,1) = voi(1:DVS,1)
         ten(1:DVS,2,2) = voi(1:DVS,2)
         ten(1:DVS,3,3) = voi(1:DVS,3)
         ten(1:DVS,2,3) = voi(1:DVS,4)
         ten(1:DVS,3,2) = voi(1:DVS,4)
         ten(1:DVS,1,3) = voi(1:DVS,5)
         ten(1:DVS,3,1) = voi(1:DVS,5)
         ten(1:DVS,1,2) = voi(1:DVS,6)
         ten(1:DVS,2,1) = voi(1:DVS,6)

      end function vcsm_uncompact_stress_tensor_3D


      pure function vcsm_uncompact_strain_tensor_2D(voi) result(ten)

         real(rp),    intent(in)    :: voi(DVS,ZVOIG)       !< Strain tensor compacted
         real(rp)                   :: ten(DVS,ZDIME,ZDIME) !< Strain tensor

         ten(1:DVS,1,1) = voi(1:DVS,1)
         ten(1:DVS,2,2) = voi(1:DVS,2)
         ten(1:DVS,1,2) = voi(1:DVS,3) * 0.5_rp
         ten(1:DVS,2,1) = voi(1:DVS,3) * 0.5_rp

      end function vcsm_uncompact_strain_tensor_2D


      pure function vcsm_uncompact_strain_tensor_3D(voi) result(ten)

         real(rp),    intent(in)    :: voi(DVS,ZVOIG)       !< Strain tensor compacted
         real(rp)                   :: ten(DVS,ZDIME,ZDIME) !< Strain tensor

         ten(1:DVS,1,1) = voi(1:DVS,1)
         ten(1:DVS,2,2) = voi(1:DVS,2)
         ten(1:DVS,3,3) = voi(1:DVS,3)
         ten(1:DVS,2,3) = voi(1:DVS,4) * 0.5_rp
         ten(1:DVS,3,2) = voi(1:DVS,4) * 0.5_rp
         ten(1:DVS,1,3) = voi(1:DVS,5) * 0.5_rp
         ten(1:DVS,3,1) = voi(1:DVS,5) * 0.5_rp
         ten(1:DVS,1,2) = voi(1:DVS,6) * 0.5_rp
         ten(1:DVS,2,1) = voi(1:DVS,6) * 0.5_rp

      end function vcsm_uncompact_strain_tensor_3D


      pure function vcsm_compact_stiffness_tensor(ten) result(voi)

         real(rp),    intent(in)    :: ten(DVS,ZDIME,ZDIME,ZDIME,ZDIME) !< Stiffness tensor
         real(rp)                   :: voi(DVS,ZVOIG,ZVOIG)             !< Stiffness tensor compacted
         integer(ip)                :: i, j, k, l, p, q

         voi = 0.0_rp
         do l = 1, ZDIME
            do k = 1, ZDIME
               do j = 1, ZDIME
                  do i = 1, ZDIME
                     p = i * Idd(i,j) + (1_ip - Idd(i,j)) * (ZVOIG + 3_ip - i - j)
                     q = k * Idd(k,l) + (1_ip - Idd(k,l)) * (ZVOIG + 3_ip - k - l)
                     voi(1:DVS,p,q) = ten(1:DVS,i,j,k,l)
                  enddo
               enddo
            enddo
         enddo
       
      end function vcsm_compact_stiffness_tensor


      pure function vcsm_uncompact_stiffness_tensor(voi) result(ten)

         real(rp),    intent(in)    :: voi(DVS,ZVOIG,ZVOIG)             !< Stiffness tensor
         real(rp)                   :: ten(DVS,ZDIME,ZDIME,ZDIME,ZDIME) !< Stiffness tensor compacted
         integer(ip)                :: i, j, k, l, p, q

         ten = 0.0_rp
         do l = 1, ZDIME
            do k = 1, ZDIME
               do j = 1, ZDIME
                  do i = 1, ZDIME
                     p = i * Idd(i,j) + (1_ip - Idd(i,j)) * (ZVOIG + 3_ip - i - j)
                     q = k * Idd(k,l) + (1_ip - Idd(k,l)) * (ZVOIG + 3_ip - k - l)
                     ten(1:DVS,i,j,k,l) = voi(1:DVS,p,q)
                  enddo
               enddo
            enddo
         enddo

      end function vcsm_uncompact_stiffness_tensor


      pure function vcsm_compact_transformation_matrix_T_3D(ten) result(voi)

         real(rp),    intent(in)    :: ten(DVS,ZDIME,ZDIME) !< Orthonormal coordinates transformation matrix 
         real(rp)                   :: voi(DVS,ZVOIG,ZVOIG) !< COmpacted transformation matrix for stress tensor

         voi(1:DVS,1,1) = ten(1:DVS,1,1)*ten(1:DVS,1,1)
         voi(1:DVS,1,2) = ten(1:DVS,1,2)*ten(1:DVS,1,2)
         voi(1:DVS,1,3) = ten(1:DVS,1,3)*ten(1:DVS,1,3)
         voi(1:DVS,1,4) = ten(1:DVS,1,2)*ten(1:DVS,1,3) * 2.0_rp
         voi(1:DVS,1,5) = ten(1:DVS,1,1)*ten(1:DVS,1,3) * 2.0_rp
         voi(1:DVS,1,6) = ten(1:DVS,1,1)*ten(1:DVS,1,2) * 2.0_rp
         voi(1:DVS,2,1) = ten(1:DVS,2,1)*ten(1:DVS,2,1)
         voi(1:DVS,2,2) = ten(1:DVS,2,2)*ten(1:DVS,2,2)
         voi(1:DVS,2,3) = ten(1:DVS,2,3)*ten(1:DVS,2,3)
         voi(1:DVS,2,4) = ten(1:DVS,2,2)*ten(1:DVS,2,3) * 2.0_rp
         voi(1:DVS,2,5) = ten(1:DVS,2,1)*ten(1:DVS,2,3) * 2.0_rp
         voi(1:DVS,2,6) = ten(1:DVS,2,1)*ten(1:DVS,2,2) * 2.0_rp
         voi(1:DVS,3,1) = ten(1:DVS,3,1)*ten(1:DVS,3,1)
         voi(1:DVS,3,2) = ten(1:DVS,3,2)*ten(1:DVS,3,2)
         voi(1:DVS,3,3) = ten(1:DVS,3,3)*ten(1:DVS,3,3)
         voi(1:DVS,3,4) = ten(1:DVS,3,2)*ten(1:DVS,3,3) * 2.0_rp
         voi(1:DVS,3,5) = ten(1:DVS,3,1)*ten(1:DVS,3,3) * 2.0_rp
         voi(1:DVS,3,6) = ten(1:DVS,3,1)*ten(1:DVS,3,2) * 2.0_rp
         voi(1:DVS,4,1) = ten(1:DVS,2,1)*ten(1:DVS,3,1)
         voi(1:DVS,4,2) = ten(1:DVS,2,2)*ten(1:DVS,3,2)
         voi(1:DVS,4,3) = ten(1:DVS,2,3)*ten(1:DVS,3,3)
         voi(1:DVS,4,4) = ten(1:DVS,2,2)*ten(1:DVS,3,3) + ten(1:DVS,2,3)*ten(1:DVS,3,2)
         voi(1:DVS,4,5) = ten(1:DVS,2,1)*ten(1:DVS,3,3) + ten(1:DVS,2,3)*ten(1:DVS,3,1)
         voi(1:DVS,4,6) = ten(1:DVS,2,1)*ten(1:DVS,3,2) + ten(1:DVS,2,2)*ten(1:DVS,3,1)
         voi(1:DVS,5,1) = ten(1:DVS,1,1)*ten(1:DVS,3,1)
         voi(1:DVS,5,2) = ten(1:DVS,1,2)*ten(1:DVS,3,2)
         voi(1:DVS,5,3) = ten(1:DVS,1,3)*ten(1:DVS,3,3)
         voi(1:DVS,5,4) = ten(1:DVS,1,2)*ten(1:DVS,3,3) + ten(1:DVS,1,3)*ten(1:DVS,3,2)
         voi(1:DVS,5,5) = ten(1:DVS,1,1)*ten(1:DVS,3,3) + ten(1:DVS,1,3)*ten(1:DVS,3,1)
         voi(1:DVS,5,6) = ten(1:DVS,1,1)*ten(1:DVS,3,2) + ten(1:DVS,1,2)*ten(1:DVS,3,1)
         voi(1:DVS,6,1) = ten(1:DVS,1,1)*ten(1:DVS,2,1)
         voi(1:DVS,6,2) = ten(1:DVS,1,2)*ten(1:DVS,2,2)
         voi(1:DVS,6,3) = ten(1:DVS,1,3)*ten(1:DVS,2,3)
         voi(1:DVS,6,4) = ten(1:DVS,1,2)*ten(1:DVS,2,3) + ten(1:DVS,1,3)*ten(1:DVS,2,2)
         voi(1:DVS,6,5) = ten(1:DVS,1,1)*ten(1:DVS,2,3) + ten(1:DVS,1,3)*ten(1:DVS,2,1)
         voi(1:DVS,6,6) = ten(1:DVS,1,1)*ten(1:DVS,2,2) + ten(1:DVS,1,2)*ten(1:DVS,2,1)

     end function vcsm_compact_transformation_matrix_T_3D


      pure function vcsm_compact_transformation_matrix_Tg_3D(ten) result(voi)

         real(rp),    intent(in)    :: ten(DVS,ZDIME,ZDIME) !< Orthonormal coordinates transformation matrix
         real(rp)                   :: voi(DVS,ZVOIG,ZVOIG) !< Compacted tranformation matrix for strain tensor

         voi(1:DVS,1,1) = ten(1:DVS,1,1)*ten(1:DVS,1,1)
         voi(1:DVS,1,2) = ten(1:DVS,1,2)*ten(1:DVS,1,2)
         voi(1:DVS,1,3) = ten(1:DVS,1,3)*ten(1:DVS,1,3)
         voi(1:DVS,1,4) = ten(1:DVS,1,2)*ten(1:DVS,1,3)
         voi(1:DVS,1,5) = ten(1:DVS,1,1)*ten(1:DVS,1,3)
         voi(1:DVS,1,6) = ten(1:DVS,1,1)*ten(1:DVS,1,2)
         voi(1:DVS,2,1) = ten(1:DVS,2,1)*ten(1:DVS,2,1)
         voi(1:DVS,2,2) = ten(1:DVS,2,2)*ten(1:DVS,2,2)
         voi(1:DVS,2,3) = ten(1:DVS,2,3)*ten(1:DVS,2,3)
         voi(1:DVS,2,4) = ten(1:DVS,2,2)*ten(1:DVS,2,3)
         voi(1:DVS,2,5) = ten(1:DVS,2,1)*ten(1:DVS,2,3)
         voi(1:DVS,2,6) = ten(1:DVS,2,1)*ten(1:DVS,2,2)
         voi(1:DVS,3,1) = ten(1:DVS,3,1)*ten(1:DVS,3,1)
         voi(1:DVS,3,2) = ten(1:DVS,3,2)*ten(1:DVS,3,2)
         voi(1:DVS,3,3) = ten(1:DVS,3,3)*ten(1:DVS,3,3)
         voi(1:DVS,3,4) = ten(1:DVS,3,2)*ten(1:DVS,3,3)
         voi(1:DVS,3,5) = ten(1:DVS,3,1)*ten(1:DVS,3,3)
         voi(1:DVS,3,6) = ten(1:DVS,3,1)*ten(1:DVS,3,2)
         voi(1:DVS,4,1) = ten(1:DVS,2,1)*ten(1:DVS,3,1) * 2.0_rp
         voi(1:DVS,4,2) = ten(1:DVS,2,2)*ten(1:DVS,3,2) * 2.0_rp
         voi(1:DVS,4,3) = ten(1:DVS,2,3)*ten(1:DVS,3,3) * 2.0_rp
         voi(1:DVS,4,4) = ten(1:DVS,2,2)*ten(1:DVS,3,3) + ten(1:DVS,2,3)*ten(1:DVS,3,2)
         voi(1:DVS,4,5) = ten(1:DVS,2,1)*ten(1:DVS,3,3) + ten(1:DVS,2,3)*ten(1:DVS,3,1)
         voi(1:DVS,4,6) = ten(1:DVS,2,1)*ten(1:DVS,3,2) + ten(1:DVS,2,2)*ten(1:DVS,3,1)
         voi(1:DVS,5,1) = ten(1:DVS,1,1)*ten(1:DVS,3,1) * 2.0_rp
         voi(1:DVS,5,2) = ten(1:DVS,1,2)*ten(1:DVS,3,2) * 2.0_rp
         voi(1:DVS,5,3) = ten(1:DVS,1,3)*ten(1:DVS,3,3) * 2.0_rp
         voi(1:DVS,5,4) = ten(1:DVS,1,2)*ten(1:DVS,3,3) + ten(1:DVS,1,3)*ten(1:DVS,3,2)
         voi(1:DVS,5,5) = ten(1:DVS,1,1)*ten(1:DVS,3,3) + ten(1:DVS,1,3)*ten(1:DVS,3,1)
         voi(1:DVS,5,6) = ten(1:DVS,1,1)*ten(1:DVS,3,2) + ten(1:DVS,1,2)*ten(1:DVS,3,1)
         voi(1:DVS,6,1) = ten(1:DVS,1,1)*ten(1:DVS,2,1) * 2.0_rp
         voi(1:DVS,6,2) = ten(1:DVS,1,2)*ten(1:DVS,2,2) * 2.0_rp
         voi(1:DVS,6,3) = ten(1:DVS,1,3)*ten(1:DVS,2,3) * 2.0_rp
         voi(1:DVS,6,4) = ten(1:DVS,1,2)*ten(1:DVS,2,3) + ten(1:DVS,1,3)*ten(1:DVS,2,2)
         voi(1:DVS,6,5) = ten(1:DVS,1,1)*ten(1:DVS,2,3) + ten(1:DVS,1,3)*ten(1:DVS,2,1)
         voi(1:DVS,6,6) = ten(1:DVS,1,1)*ten(1:DVS,2,2) + ten(1:DVS,1,2)*ten(1:DVS,2,1)

     end function vcsm_compact_transformation_matrix_Tg_3D

end module mod_sld_vect_csm
!> @}
