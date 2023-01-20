!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_builtin_materials.f90
!> @author  Solidz team
!> @date    2020-04-02
!> @brief   Interface for calling built-in subroutines
!> @details Interface for calling built-in subroutines
!   INPUT
!      PGAUS ... Number of Gauss points
!      PMATE ... Material number
!      TEMP0 ... Previous temperature ............................ temp(n)
!      TEMP1 ... Updated temperature ............................. temp(n+1)
!      GPGI0 ... Previous deformation gradient tensor ............ F(n)
!      GPGDI ... Updated deformation gradient tensor ............. F(n+1)
!      GPIGD ... Inverse of updated deformation gradient tensor .. F^{-1}(n+1)
!      GPCAU ... Updated right Cauchy-Green tensor ............... C(n+1)
!      GPDET ... Updated jacobian ................................ J = det(F(n+1))
!      GPRAT ... Rate of Deformation tensor ...................... Fdot = grad(phidot)
!      FLAGT ... Flag for tangent moduli calculations ............ 1 if yes; 0 if no
!      FLAGS ... Flag for strain gradient calculations ........... 1 if yes; 0 if no
!   INPUT/OUTPUT
!      GPPIO ... 1st Piola-Kirchhoff stress tensor at t(n)/t(n+1)  P(n)/P(n+1)
!      GPPIO_EPS ... Perturbed 1st Piola-Kirchhoff stress tensor at t(n)/t(n+1)  P(n)/P(n+1)
!      GPSTR ... 2nd Piola-Kirchhoff stress tensor at t(n)/t(n+1)  S(n)/S(n+1)
!      GPENE ... Stored energy function .......................... W
!      GPGRE ... Green-Lagrange         .......................... E
!   OUTPUT
!      GPTMO ... Tangent moduli at t(n+1) ........................ dP/dF(n+1)
!      GGSGO ... Deformation gradient gradient tensor ............ dF/dX(n+1)
!> @}
!-----------------------------------------------------------------------
subroutine sld_builtin_materials( &
   pgaus,pmate,gpvol,gptmp,gpgi0,gpgdi,gpigd,gpcau,gpdet,gprat,gppio,gpstr,&
   gpene,flagt,gptmo,flags,ggsgo,ielem,gpfib,elcod,pnode,lnods,gpsha,gpdds,&
   gpigd_eps,gpgdi_eps,gpdet_eps,gppio_eps,gpmof)
   ! -------------------------------
   use def_kintyp,               only: ip,rp
   use def_domain,               only: ndime,mnode
   use mod_sld_fe2,              only: MAT_MICRO_NO_COUPLING,MAT_MICRO_ONE_WAY,MAT_MICRO_FULL
   use def_solidz,               only: kfl_tange_sld, SLD_TANGENT_NUMERICAL, &
   &                                   lawst_sld, lawho_sld, lawpl_sld, svegm_sld
   use mod_sld_stress_model_152, only: sld_stress_model_152
   use mod_sld_stress_model_154, only: sld_stress_model_154
   use mod_sld_stress_model_400, only: sld_stress_model_400
   use sld_stress_model_160,     only: eval_sld_stress_model_160
   use mod_sld_plastic_model,    only: sld_plastic_model_biso
   use mod_sld_stress_model_170, only: sm170_get_stress
   use mod_sld_fe2,              only: fe2_get_stress_and_ctan
   use mod_sld_atm,              only: ntmp_atm
   use mod_biofibers,            only: biofibers
   use def_solidz,               only: parco_sld
   ! -------------------------------
   implicit none
   ! -------------------------------
   integer(ip), intent(in)        :: pgaus                                ! Number of gauss points 
   integer(ip), intent(in)        :: pmate                                ! Number of material
   integer(ip), intent(in)        :: ielem                                ! Number of element
   integer(ip), intent(in)        :: pnode                                ! Number of nodes
   integer(ip), intent(in)        :: lnods(pnode)                         ! List of nodes
   real(rp),    intent(in)        :: elcod(ndime,mnode)                   ! Element coordinates 
   real(rp),    intent(in)        :: gpsha(pnode,pgaus)                   ! Shape functions
   real(rp),    intent(in)        :: gptmp(ntmp_atm,pgaus)                ! Temperature
   real(rp),    intent(in)        :: gpgi0(ndime,ndime,pgaus)             ! Deformation gradient tensor (previous)
   real(rp),    intent(in)        :: gpgdi(ndime,ndime,pgaus)             ! Deformation gradient tensor (current)
   real(rp),    intent(in)        :: gpigd(ndime,ndime,pgaus)             ! Inverse of the deformation gradient tensor  (current)
   real(rp),    intent(in)        :: gpcau(ndime,ndime,pgaus)             ! Right Cauchy straint ensor (current)
   real(rp),    intent(in)        :: gpdet(pgaus),gpmof(pgaus)            ! Determinant of the deformation gradient tensor (current)
   real(rp),    intent(in)        :: gprat(ndime,ndime,pgaus)             ! Rate of the deformation tensor (current)
   real(rp),    intent(in)        :: gpvol(pgaus)                         ! Volume 
   integer(ip), intent(in)        :: flagt                                ! Material tangent module required? 
   integer(ip), intent(in)        :: flags                                ! Other flags
   real(rp),    intent(inout)     :: gpene(pgaus)                         ! Stored energy function
   real(rp),    intent(inout)     :: gpstr(ndime,ndime,pgaus,2)           ! Stress tensor computed by the material model
   real(rp),    intent(out)       :: gpdds(ndime,ndime,ndime,ndime,pgaus) ! Elasticity tensor computed by the material model 
   real(rp),    intent(inout)     :: gppio(ndime,ndime,pgaus)             ! Stress tensor needed by Alya (P)
   real(rp),    intent(out)       :: gptmo(ndime,ndime,ndime,ndime,pgaus) ! Elasticity tensor needed by Alya (dPdF)
   real(rp),    intent(out)       :: gpfib(ndime,pgaus)                   ! Fiber direction used by the material model
   real(rp),    intent(inout)     :: gppio_eps(ndime,ndime,pgaus)         ! Stress tensor P for the perturbated state
   real(rp),    intent(in)        :: gpgdi_eps(ndime,ndime,pgaus)         ! Displacement Gradient F for the perturbed state
   real(rp),    intent(in)        :: gpigd_eps(ndime,ndime,pgaus)         ! Displacement Gradient Inverse F^{-1} for the pert. state
   real(rp),    intent(in)        :: gpdet_eps(pgaus)                     ! Determinant Gradient F for the perturbed state
   real(rp),    intent(out)       :: ggsgo(ndime,ndime)                   ! Deformation gradient gradient tensor
   ! -------------------------------
   integer(ip)                    :: igaus, ii, jj, kk, ll, mm, nn
   integer(ip)                    :: kfl_STR                              ! Flag for stress tensor
   integer(ip)                    :: kfl_MAT                              ! Flag for tangent material tensor
   real(rp)                       :: tkron(ndime,ndime)
   real(rp)                       :: gpfsn(ndime,ndime,pgaus)             ! Fiber-Sheet-Normal basis for bio-materials
   real(rp)                       :: temp0(pgaus)                         ! Temperature  (previous)
   real(rp)                       :: temp1(pgaus)                         ! Temperatrue  (current)
   ! -------------------------------
   integer(ip), parameter         :: STR_FROM_PK2_TO_PK1   = 1_ip
   integer(ip), parameter         :: MAT_NOT_COMPUTED      = 0_ip
   integer(ip), parameter         :: MAT_FROM_DSDE_TO_DPDF = 1_ip
   ! -------------------------------
   !
   ! Initialisation
   !
   ! Flag fro debugging
   kfl_STR = -1_ip
   kfl_MAT = -1_ip 
   ! Kronecker delta
   tkron = 0.0_rp
   do ii = 1, ndime
       tkron(ii,ii) = 1.0_rp
   end do

   !
   ! Tangent moduli (only for implicit computations)
   !
   if( (kfl_tange_sld == SLD_TANGENT_NUMERICAL .or. kfl_tange_sld == 2_ip) .and. flagt==1_ip )then
      !
      ! Compute stress and elasticity tensors, but the tangent moduli numerical
      ! 
      ggsgo(1,1)       = 0.0_rp
      gptmo(1,1,1,1,1) = 0.0_rp
      temp0(:) = gptmp(3,:)
      temp1(:) = gptmp(1,:)
      call sld_dertan(pgaus,pmate,temp0,temp1,&
         gpgi0,gpgdi,gpigd,gpcau,gpdet,gprat,gppio,gpstr,&
         gpene,flagt,gptmo,flags,ggsgo,ielem,gpfib,elcod,pnode,lnods,gpsha,gpdds,gpmof)

   else
      !
      ! Compute stress and elasticity tensors calling the built-in material law
      !
      if(      lawst_sld(pmate) == 100_ip )then
         !
         ! ISOLIN : Linear isotropic model
         !   {Ref:  Belytschko 2n edition book}
         kfl_STR = STR_FROM_PK2_TO_PK1
         kfl_MAT = MAT_FROM_DSDE_TO_DPDF
         call sld_stress_model_100(pgaus,pmate,gptmp,gpgdi,gpstr,gpcau,gpdet, &
         &   gpigd,ielem,elcod,flagt,gpdds,gpigd_eps,gpgdi_eps,   &
         &   gpdet_eps,gpmof)

      else if( lawst_sld(pmate) == 101_ip )then
         !                                
         ! Neo-Hookean
         !   {Ref:  Belytschko 2n edition book}
         kfl_STR = STR_FROM_PK2_TO_PK1
         kfl_MAT = MAT_FROM_DSDE_TO_DPDF
         call sld_stress_model_101(pgaus,pmate,gpcau,gpgdi,gpene,gpstr, &
         &   gpdet,flagt,gpdds,gpigd_eps,gpgdi_eps,gpdet_eps,gpmof)

      else if( lawst_sld(pmate) == 102_ip )then
         !
         ! Neo-Hookean Ansys
         !   {Ref: Ansys documentation}
         kfl_STR = STR_FROM_PK2_TO_PK1
         kfl_MAT = MAT_FROM_DSDE_TO_DPDF
         call sld_stress_model_102(pgaus,pmate,gpgdi,gpstr,gpdet,flagt, &
         &   gpdds,gpigd_eps,gpgdi_eps,gpdet_eps,gpmof)

      else if( lawst_sld(pmate) == 103_ip )then
         !
         ! Mooney-Rivlin 
         !   {Ref: Belytschko notation)
         kfl_STR = STR_FROM_PK2_TO_PK1
         kfl_MAT = MAT_FROM_DSDE_TO_DPDF
         call sld_stress_model_103(pgaus,pmate,gpgdi,gpstr,gpdet,flagt, &
         &   gpdds,gpigd_eps,gpgdi_eps,gpdet_eps,gpmof)

      else if( lawst_sld(pmate) == 105_ip )then
         !
         ! Constant spring
         kfl_STR = STR_FROM_PK2_TO_PK1
         kfl_MAT = MAT_NOT_COMPUTED
         call sld_stress_model_105(pgaus,pmate,gpstr,flagt,gpdds,&
         &   gpigd_eps,gpgdi_eps,gpdet_eps,gpmof)

      else if( lawst_sld(pmate) == 121_ip  )then
         !
         ! Hyperelastic model from Chappelle and co.
         kfl_STR = STR_FROM_PK2_TO_PK1
         kfl_MAT = MAT_NOT_COMPUTED
         call sld_stress_model_121(pgaus,pmate,gpcau,gpgdi,gpene,&
         &  gpstr,gprat,gpdet) 

      else if( lawst_sld(pmate) == MAT_MICRO_NO_COUPLING .or. &
            &    lawst_sld(pmate) == MAT_MICRO_ONE_WAY     .or. &
            &    lawst_sld(pmate) == MAT_MICRO_FULL       )then
         !
         ! MicroPP microscale model
         kfl_STR = STR_FROM_PK2_TO_PK1
         kfl_MAT = MAT_FROM_DSDE_TO_DPDF
         call fe2_get_stress_and_ctan(ielem, pgaus, gpstr, gpdds)

      else if( lawst_sld(pmate) == 134_ip )then
            !
            ! Holza & Ogden
            kfl_STR = STR_FROM_PK2_TO_PK1
            kfl_MAT = MAT_FROM_DSDE_TO_DPDF
            call biofibers % get_reference_basis_at_gp(&
            &   ielem,pnode,pgaus,gpsha,gpfsn)
            call sld_stress_model_134(&
            &   pgaus,pmate,flagt,gpdet,gpcau,gpfsn,gpstr,gpdds)
      
      else if( lawst_sld(pmate) == 135_ip )then
         !
         ! Usyk 
         kfl_STR = STR_FROM_PK2_TO_PK1
         kfl_MAT = MAT_NOT_COMPUTED
         call sld_stress_model_135(pgaus,pmate,gpcau,gpgdi,gpigd,gpstr,  &
         &   gpdet,pnode,elcod,lnods,ielem,gpsha,gpfib,&
         &   gpigd_eps,gpgdi_eps,gpdet_eps,gpmof)

      else if( lawst_sld(pmate) == 136_ip )then
         !
         ! Guccione
         kfl_STR = STR_FROM_PK2_TO_PK1
         kfl_MAT = MAT_FROM_DSDE_TO_DPDF
         call sld_stress_model_136(pgaus,pmate,gpcau,gpgdi,gpigd,gpstr,  &
         &   gpdet,gpfib,ielem,elcod,pnode,lnods,gpsha,flagt,      &
         &   gpdds,&
         &   gpigd_eps,gpgdi_eps,gpdet_eps,gpmof)

      else if( lawst_sld(pmate) == 137_ip )then
         !
         ! Holza & Ogden for arteries
         kfl_STR = STR_FROM_PK2_TO_PK1
         kfl_MAT = MAT_FROM_DSDE_TO_DPDF
         if (lawho_sld(pmate)==4) then
            call sld_stress_model_137hgoc(pgaus,pmate,gpvol,gpcau,gpgdi,&
            &   gpigd,gpstr,gpdet,gpfib,ielem,elcod,pnode,lnods,  &
            &   gpsha,flagt,gpdds,&
            &   gpigd_eps,gpgdi_eps,gpdet_eps,gpmof)

         else if (lawho_sld(pmate)==5) then
            call sld_stress_model_137ma(pgaus,pmate,gpvol,gpcau,gpgdi,  &
            &   gpigd,gpstr,gpdet,gpfib,ielem,elcod,pnode,lnods,  &
            &   gpsha,flagt,gpdds,gpigd_eps,gpgdi_eps,gpdet_eps,gpmof)
         else
            call sld_stress_model_137hgoc(pgaus,pmate,gpvol,gpcau,gpgdi,&
            &   gpigd,gpstr,gpdet,gpfib,ielem,elcod,pnode,lnods,  &
            &   gpsha,flagt,gpdds,&
            &   gpigd_eps,gpgdi_eps,gpdet_eps,gpmof)
         end if

      else if( lawst_sld(pmate) == 140_ip )then
            !
            ! Holza & Ogden with compressibility modification (Nolan, 2014)
            kfl_STR = STR_FROM_PK2_TO_PK1
            kfl_MAT = MAT_FROM_DSDE_TO_DPDF
            call sld_stress_model_140(pgaus,pmate,gpvol,gpcau,gpgdi,gpigd,  &
            &   gpstr,gpdet,gpfib,ielem,elcod,pnode,lnods,gpsha,flagt,&
            &   gpdds, &
            &   gpigd_eps,gpgdi_eps,gpdet_eps,gpmof)

      else if( lawst_sld(pmate) == 151_ip )then
            !
            ! Orthotropic linear-elastic model
            !   {Ref : Belytschko book }
            kfl_STR = STR_FROM_PK2_TO_PK1
            kfl_MAT = MAT_FROM_DSDE_TO_DPDF
            call sld_stress_model_151(pgaus,pmate,gpgdi,gpigd,gpdet,gptmp,ielem,  &
            &   flagt,gpstr,gpdds)

      else if( lawst_sld(pmate) == 152_ip )then
            !
            ! Damage model: transversally-isotropic 3D damage model
            !   {Ref : Quintanas-Corominas et al. 2018}
            kfl_STR = STR_FROM_PK2_TO_PK1
            kfl_MAT = MAT_FROM_DSDE_TO_DPDF
            call sld_stress_model_152(pgaus,pmate,gpgdi,gpigd,gpdet,gpstr,  &
            &   ielem,flagt,gpdds)
            
      else if( lawst_sld(pmate) == 153_ip )then
            !
            ! Damage model: isotropic-linear
            !   {Ref : Oliver et al. 1990}
            kfl_STR = STR_FROM_PK2_TO_PK1
            kfl_MAT = MAT_FROM_DSDE_TO_DPDF
            call sld_stress_model_153(pgaus,pmate,gpgdi,gpstr,ielem,flagt,  &
            &   gpdds)

      else if( lawst_sld(pmate) == 154_ip )then
            ! 
            ! Damage model: transversally-isotropic 2D/3D damage model
            !   {Ref : Maimí et al. 2007a/b}
            kfl_STR = STR_FROM_PK2_TO_PK1
            kfl_MAT = MAT_NOT_COMPUTED
            call sld_stress_model_154(pgaus,pmate,gpgdi,gpigd,gpdet,ielem,  &
            &   flagt,gpstr,gpdds)

      else if( lawst_sld(pmate) == 170_ip )then
         ! 
         ! Thermo-elasticity model: isotropy
         !
         kfl_STR = STR_FROM_PK2_TO_PK1
         kfl_MAT = MAT_FROM_DSDE_TO_DPDF
         call sm170_get_stress(flagt, pgaus, parco_sld(:,pmate), gpgdi, gptmp, & 
              &    gpstr(:,:,:,1), gpdds)
         
      else if( lawst_sld(pmate) == 400_ip )then
         ! 
         ! Poroelasticity
         !    {Ref:  Schroeder et. al 2008}
         !                                
         kfl_STR = STR_FROM_PK2_TO_PK1
         kfl_MAT = MAT_FROM_DSDE_TO_DPDF
         call sld_stress_model_400(&
              & ielem, pgaus, pnode, pmate, gpgdi, gpstr, gpdet, flagt, gpdds, &
              & gpene, gpcau, gpsha)
              
      else if( lawst_sld(pmate) == 160_ip )then
            !
            ! Shape memory alloy model
            !   { Ref_A : 10.1016/S0045-7825(96)01147-4
            !     Ref_B : 10.1016/S0749-6419(00)00050-4
            !     Ref_D : 10.1007/s00466-017-1518-9}

            kfl_STR = STR_FROM_PK2_TO_PK1
            kfl_MAT = MAT_FROM_DSDE_TO_DPDF
            


            call eval_sld_stress_model_160( &
            pgaus, &
            parco_sld(1:18, pmate), &
            gpgdi(1:ndime, 1:ndime, 1:pgaus), &
            temp1(1:pgaus), &
            svegm_sld(ielem)%a(1:ndime, 1:pgaus, 1:2), &
            gpstr(1:ndime, 1:ndime, 1:pgaus, 1), &
            gpdds(1:ndime, 1:ndime, 1:ndime, 1:ndime, 1:pgaus), &
            flagt)

      else if( lawpl_sld(pmate) == 1_ip   .and. &
               lawst_sld(pmate) == 0_ip  ) then
            !
            ! Plasticity models
            kfl_STR = STR_FROM_PK2_TO_PK1
            kfl_MAT = MAT_FROM_DSDE_TO_DPDF
            call sld_plastic_model_biso(pgaus,gpgdi,gpigd,gpdet,gpstr,      &
            &   ielem,flagt,gpdds)
      else
         call runend('SLD_BUILTIN_MATERIALS : laws_sld(pmate) not defined') 

      end if

   endif
   !
   ! Stress and Elasticity tensors transformations
   !
   if( kfl_STR == STR_FROM_PK2_TO_PK1 )then
       !
       ! Calculate 1st Piola-Kirchhoff from 2nd Piola-Kirchhoff
       !   P_jK = F_jI * S_IK
       !   {Ref : Belytschko book}
       gppio(:,:,:) = 0.0_rp
       do igaus=1,pgaus
           do kk=1,ndime
               do jj=1,ndime
                   do ii=1,ndime
                       gppio(jj,kk,igaus)= gppio(jj,kk,igaus) + &
                       &   gpgdi(jj,ii,igaus)*gpstr(ii,kk,igaus,1)
                   enddo
               enddo
           enddo
       end do

   endif

   !
   ! Elasticity tensor transormations
   !
   if( flagt == 1_ip .and. kfl_MAT == MAT_NOT_COMPUTED )then
       !
       ! We need a tangent but the model didnt compute one
       call runend('SLD_BUILTIN_MATERIALS: STRESS TANGENT NOT PROGRAMMED')

   elseif( flagt == 1_ip .and. kfl_MAT == MAT_FROM_DSDE_TO_DPDF ) then
       !
       ! Calculate derivative of 1st Piola-Kirchhoff stress tensor wrt.
       ! deformation gradient (dPdF) from the derivative of the 2nd 
       ! Piola-Kirchhoff stress tensor wrt. the Green-Lagrange (dSdE)
       !   dPdF_iJkL = delta_ik * S_JL + F_iM * FkN * dSdE_MJNL
       !   {Ref : Belytschko book}
       gptmo(:,:,:,:,:) = 0.0_rp
       do igaus=1,pgaus
           do ii=1,ndime
               do jj=1,ndime
                   do kk=1,ndime 
                       do ll=1,ndime
                           do mm=1,ndime
                               do nn=1,ndime
                                   gptmo(ii,jj,kk,ll,igaus)= gptmo(ii,jj,kk,ll,igaus) + &
                                   &    gpdds(mm,jj,nn,ll,igaus)*gpgdi(ii,mm,igaus)*gpgdi(kk,nn,igaus)
                               enddo
                           enddo
                           gptmo(ii,jj,kk,ll,igaus) = gptmo(ii,jj,kk,ll,igaus) + &
                           &   tkron(ii,kk)*gpstr(jj,ll,igaus,1)
                       enddo
                   enddo
               enddo
           enddo
       enddo

   else
      !
      ! Set dPdS to 0.0 (just in case)
      !
      gptmo(:,:,:,:,:) = 0.0_rp

   endif
   ! -------------------------------
end subroutine sld_builtin_materials

