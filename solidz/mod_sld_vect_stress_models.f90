!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_vect_stress_models.f90
!> @author  Adria Quintanas-Corominas
!> @date    November 2021
!> @brief   Toolbox with some stress models for the vectorised assembly
!> @details Toolbox with some stress models for the vectorised assembly
!------------------------------------------------------------------

module mod_sld_vect_stress_models

#include "def_solidz_vect_dim.inc" 
   use def_kintyp, only : ip, rp, lg
   use def_domain, only : ndime
   use def_solidz, only : parco_sld, kfl_plane_sld
   use mod_sld_vect_csm 
   use mod_sld_vect_maths
   
   implicit none

   integer(ip), parameter :: SM_ISOLIN = 100_ip
   integer(ip), parameter :: SM_ORTLIN = 151_ip
   integer(ip), parameter :: SM_BESSA  = 154_ip
   
   real(rp), dimension(:,:,:,:), pointer :: vdsde_sld

   contains

      subroutine sld_vect_stress_model_100(itask,pmate,pgaus,gpgre,gpstr,gpdds)
         ! Sij = Cijkl Ekl = lambda * trace(E) delta_ij + 2 mu * Eij
         ! Cijkl = lambda * dij * dkl  +  mu * ( dik * djl + dil * djk )

         integer(ip), intent(in)    :: itask
         integer(ip), intent(in)    :: pgaus
         integer(ip), intent(in)    :: pmate
         real(rp),    intent(in)    :: gpgre(DVS,ZDIME,ZDIME,ZGAUS)
         real(rp),    intent(out)   :: gpstr(DVS,ZDIME,ZDIME,ZGAUS)
         real(rp),    intent(out)   :: gpdds(DVS,ZDIME,ZDIME,ZDIME,ZDIME,ZGAUS)
         integer(ip)                :: iv, ig, ii, jj, kk, ll
         real(rp)                   :: yng, nu, lmb, mu
         real(rp)                   :: lmbN(DVS), muN(DVS)
         real(rp)                   :: trace(DVS,ZGAUS)
         real(rp),    parameter     :: one = 1.0_rp, zero = 0.0_rp
         real(rp),    parameter     :: tkron(3,3) = reshape( &
                                          (/ one, zero, zero, zero, one, zero, zero, zero, one /) , &
                                          (/ 3, 3 /) )

         ! Recorver properties
         yng = parco_sld(1,pmate)
         nu  = parco_sld(2,pmate)
         
         lmb = (yng*nu)/((1.0_rp+nu)*(1.0_rp-2.0_rp*nu)) 
         mu  = yng/(2.0_rp*(1.0_rp+nu))
         
         ! Modify according to kinematic assumption
         if( kfl_plane_sld == 1_ip )then 
            ! 2D plane stress assumption
            do iv = 1, DVS
               lmbN(iv) = lmb * (1.0_rp - 2.0_rp * nu) / (1.0_rp - nu)
               muN(iv) = mu
            enddo
         else
            ! 2D plane strain assumption 
            do iv = 1, DVS
               lmbN(iv) = lmb
               muN(iv) = mu
            enddo
         endif
         ! Trace
         do ig = 1, ZGAUS
            trace(:,ig) = vmath_TRACE(gpgre(:,:,:,ig))
         enddo
         ! Second Piola-Kirchooff Stress
         do ig = 1, ZGAUS
            do ii = 1, ZDIME
               do jj = 1, ZDIME
                  gpstr(1:DVS,ii,jj,ig) = lmbN(1:DVS) * trace(1:DVS,ig) * tkron(ii,jj) + &
                     2.0_rp * muN(1:DVS) * gpgre(1:DVS,ii,jj,ig)
               enddo
            enddo
         enddo
         ! Second Elasticity tesnor
         if( itask == 2_ip )then
            do ig = 1, ZGAUS
               do ii = 1, ZDIME
                  do jj = 1, ZDIME
                     do kk = 1, ZDIME
                        do ll = 1, ZDIME
                           gpdds(1:DVS,ii,jj,kk,ll,ig) = lmbN(1:DVS) * tkron(ii,jj) * tkron(kk,ll) + &
                              muN(1:DVS) * ( tkron(ii,kk) * tkron(jj,ll) + tkron(ii,ll) * tkron(jj,kk) )
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         else
            gpdds = 0.0_rp
         endif

      end subroutine sld_vect_stress_model_100


      subroutine sld_vect_stress_model_151(itask,pmate,pgaus,elsys,gpgre,gpstr,gpdds)

         integer(ip), intent(in)    :: itask
         integer(ip), intent(in)    :: pgaus
         integer(ip), intent(in)    :: pmate
         real(rp),    intent(in)    :: elsys(DVS,ZDIME,ZDIME)
         real(rp),    intent(in)    :: gpgre(DVS,ZDIME,ZDIME,ZGAUS)
         real(rp),    intent(out)   :: gpstr(DVS,ZDIME,ZDIME,ZGAUS)
         real(rp),    intent(out)   :: gpdds(DVS,ZDIME,ZDIME,ZDIME,ZDIME,ZGAUS)
         real(rp)                   :: gpepsXYZ(DVS,ZVOIG,ZGAUS) 
         real(rp)                   :: gpeps123(DVS,ZVOIG,ZGAUS) 
         real(rp)                   :: gpsig123(DVS,ZVOIG,ZGAUS) 
         real(rp)                   :: gpsigXYZ(DVS,ZVOIG,ZGAUS) 
         real(rp)                   :: gpdse123(DVS,ZVOIG,ZVOIG,ZGAUS) 
         real(rp)                   :: gpdseXYZ(DVS,ZVOIG,ZVOIG,ZGAUS)
         integer(ip)                :: ig

         ! Compact E -> eps
         call sld_vcsm_compact_strain_at_GP(ZGAUS,gpgre,gpepsXYZ)       
         ! Rotate eps -> eps_123
         call sld_vcsm_compact_transformation_strain_at_GP(0_ip,ZGAUS,elsys,gpepsXYZ,gpeps123)
         ! Compute sig_123 
         do ig = 1, ZGAUS
             ! Recover dsde_123 (material stiffness tensor)
             gpdse123(:,:,:,ig) = vdsde_sld(:,:,:,pmate)
             ! Compute sig_123
             gpsig123(:,:,ig) = vmath_MxV(gpdse123(:,:,:,ig),gpeps123(:,:,ig))
         enddo
         ! Rotate sig_123 -> sig_xyz
         call sld_vcsm_compact_transformation_stress_at_GP(1_ip,ZGAUS,elsys,gpsig123,gpsigXYZ)
         ! Uncompact sig_xyz -> S
         call sld_vcsm_uncompact_stress_at_GP(ZGAUS,gpsigXYZ,gpstr)
         ! Second elasticity tensor
         if( itask == 2_ip )then
            ! Rotate K_123 -> K_xyz
            call sld_vcsm_compact_transformation_stiffness_at_GP(1_ip,ZGAUS,elsys,gpdse123,gpdseXYZ)
            ! Uncompact K_xyz -> dSdE
            call sld_vcsm_uncompact_stiffness_tensor_at_GP(ZGAUS,gpdseXYZ,gpdds)
         else
            gpdds = 0.0_rp
         endif

      end subroutine sld_vect_stress_model_151


      subroutine sld_vect_stress_model_151b(itask,pmate,pgaus,elsys,gpgre,gpstr,gpdds)

         integer(ip), intent(in)    :: itask
         integer(ip), intent(in)    :: pgaus
         integer(ip), intent(in)    :: pmate
         real(rp),    intent(in)    :: elsys(DVS,ZDIME,ZDIME)
         real(rp),    intent(in)    :: gpgre(DVS,ZDIME,ZDIME,ZGAUS)
         real(rp),    intent(out)   :: gpstr(DVS,ZDIME,ZDIME,ZGAUS)
         real(rp),    intent(out)   :: gpdds(DVS,ZDIME,ZDIME,ZDIME,ZDIME,ZGAUS)
         real(rp)                   :: gpepsXYZ(DVS,ZVOIG) 
         real(rp)                   :: gpeps123(DVS,ZVOIG) 
         real(rp)                   :: gpsig123(DVS,ZVOIG) 
         real(rp)                   :: gpsigXYZ(DVS,ZVOIG) 
         real(rp)                   :: gpTg(DVS,ZVOIG,ZVOIG) 
         real(rp)                   :: gpdse123(DVS,ZVOIG,ZVOIG) 
         real(rp)                   :: gpdseXYZ(DVS,ZVOIG,ZVOIG) 
         integer(ip)                :: ig, ii, jj, kk, ll

         ! Get transformation matrix for compacted notation ratations
         gpTg = vcsm_compact_transformation_matrix_Tg_3D(elsys)
         ! Recover the elastic material stiffness at the material CSYS
         gpdse123(:,:,:) = vdsde_sld(:,:,:,pmate)
         ! Compute the PK2 stress tensor (S)
         do ig = 1, ZGAUS
             gpepsXYZ(:,:) = vcsm_compact_strain_tensor_3D(gpgre(:,:,:,ig))
             gpeps123(:,:) = vmath_MxV(gpTg,gpepsXYZ)
             gpsig123(:,:) = vmath_MxV(gpdse123,gpeps123)
             gpsigXYZ(:,:) = vmath_TXV(gpTg,gpsig123)
             gpstr(:,:,:,ig) = vcsm_uncompact_stress_tensor_3D(gpsigXYZ)
         enddo
         ! Second elasticity tensor (dS/dE)
         if( itask == 2_ip )then
            do ig = 1, ZGAUS
               ! vmath_TxMxM
               gpdseXYZ = 0.0_rp
               do ii = 1, 6; do jj = 1, 6; do kk = 1, 6; do ll = 1, 6
                  gpdseXYZ(1:DVS,ii,jj) = gpdseXYZ(1:DVS,ii,jj) + &
                        gpTg(1:DVS,kk,ii) * gpTg(1:DVS,ll,jj) * gpdse123(1:DVS,kk,ll)
               enddo; enddo; enddo; enddo
               gpdds(:,:,:,:,:,ig) = vcsm_uncompact_stiffness_tensor(gpdseXYZ) 
            enddo
         else
            gpdds = 0.0_rp
         endif

      end subroutine sld_vect_stress_model_151b


      subroutine sld_stress_model_151_compute_properties(pmate)
         use def_master, only : IEMPTY
         use def_solidz, only : nmate_sld

         integer(ip), intent(in) :: pmate
         real(rp)                :: E11, E22, E33, v12, v13, v23, G12, G13, G23, v21, v31, v32
         real(rp)                :: delta, auxS1

         if( IEMPTY ) return

         if( .not. associated(vdsde_sld) )then
            allocate(vdsde_sld(DVS,ZVOIG,ZVOIG,nmate_sld))
         endif
         
         vdsde_sld = 0.0_rp

         E11 = parco_sld(1,pmate)
         E22 = parco_sld(2,pmate)
         E33 = parco_sld(3,pmate)
         v12 = parco_sld(4,pmate)
         v13 = parco_sld(5,pmate)
         v23 = parco_sld(6,pmate)
         G12 = parco_sld(7,pmate)
         G13 = parco_sld(8,pmate)
         G23 = parco_sld(9,pmate)
         v21 = (v12*E22)/E11
         v31 = (v13*E33)/E11
         v32 = (v23*E33)/E22
         delta = (1.0_rp - v12*v21 - v23*v32 - v31*v13 - 2.0_rp*v12*v23*v31)/(E11*E22*E33)
         auxS1 = E22*E33*delta
         vdsde_sld(1:DVS,1,1,pmate) = (1.0_rp - v23*v32)/auxS1
         vdsde_sld(1:DVS,1,2,pmate) = (v21 + v31*v23)/auxS1
         vdsde_sld(1:DVS,1,3,pmate) = (v31 + v21*v32)/auxS1
         auxS1 = E33*E11*delta
         vdsde_sld(1:DVS,2,1,pmate) = (v12 + v13*v32)/auxS1
         vdsde_sld(1:DVS,2,2,pmate) = (1.0_rp - v31*v13)/auxS1
         vdsde_sld(1:DVS,2,3,pmate) = (v32 + v31*v12)/auxS1
         auxS1 = E11*E22*delta
         vdsde_sld(1:DVS,3,1,pmate) = (v13 + v12*v23)/auxS1
         vdsde_sld(1:DVS,3,2,pmate) = (v23 + v13*v21)/auxS1
         vdsde_sld(1:DVS,3,3,pmate) = (1.0_rp - v12*v21)/auxS1
         vdsde_sld(1:DVS,4,4,pmate) = G23
         vdsde_sld(1:DVS,5,5,pmate) = G13
         vdsde_sld(1:DVS,6,6,pmate) = G12

      end subroutine sld_stress_model_151_compute_properties

end module mod_sld_vect_stress_models
!> @}
