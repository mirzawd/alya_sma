!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!----------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_stress_model_151.f90
!> @author  Adria Quintanas (adria.quintanas@udg.edu)
!> @date    November, 2015
!>          - Subroutine written
!> @author  Gerard Guillamet
!> @date    December, 2017
!>          - Adds strain measures corrections
!> @author  Gerard Guillamet
!> @date    June, 2022
!>          - Adds rotation matrix for material orientations
!> @brief   Sant Venant - Kirchoff Orthotropic material model
!>
!> @details
!>
!>          References:\n
!>          T. Belytschko, W. K. Liu, B. Moran, K. I. Elkhodary
!>          Nonlinear Finite elements for Continua and Structures\n
!>          E. J. Barbero. Introduction to Composite Materials Design\n
!>
!> @}
!----------------------------------------------------------------------------

subroutine sld_stress_model_151(pgaus,pmate,gpgdi,gpidg,gpdet,gptmp,ielem,flagt,gpstr,gpdds)

  use def_kintyp,                  only : ip,rp
  use def_domain,                  only : ndime
  use def_solidz,                  only : nvoig_sld,stiff0_sld,rmate_sld
  use def_solidz,                  only : kfl_strai_sld, SLD_INFINITESIMAL, SLD_GREEN
  use mod_sld_stress_model_comput, only : sm_rotate_basis_creation
  use mod_sld_stress_model_comput, only : sm_rotate_voigt_second
  use mod_sld_stress_model_comput, only : sm_rotate_voigt_fourth
  use mod_sld_stress_model_comput, only : sm_strain_tensor
  use mod_sld_stress_model_comput, only : sm_stress_tensor
  use mod_sld_stress_model_comput, only : sm_tensor_to_voigt_fourth
  use mod_sld_stress_model_comput, only : sm_tensor_to_voigt_second, SM_voigt_to_tensor_second
  use mod_sld_stress_model_comput, only : SM_stress_transport
  use mod_sld_stress_model_comput, only : SM_stiffness_transport
  use mod_sld_stress_model_comput, only : SM_PULLBACK
  use mod_sld_atm,                 only : ntmp_atm

  implicit none

  ! ================================================================================
  ! INIT
  ! --------------------------------------------------------------------------------
  integer(ip), intent(in)  :: pgaus                                !< No. gauss points
  integer(ip), intent(in)  :: pmate                                !< Material code
  integer(ip), intent(in)  :: ielem                                !< Current element number
  integer(ip), intent(in)  :: flagt                                !< Integration scheme flag: 0: explicit / 1: implicit
  real(rp),    intent(in)  :: gpgdi(ndime,ndime,pgaus)             !< Displacement Gradient
  real(rp),    intent(in)  :: gpidg(ndime,ndime,pgaus)             !< Inverse of updated deformation gradient tensor
  real(rp),    intent(in)  :: gpdet(pgaus)                         !< Updated Jacobian
  real(rp),    intent(in)  :: gptmp(ntmp_atm,pgaus)                !< Temperature
  real(rp),    intent(out) :: gpstr(ndime,ndime,pgaus)             !< 2nd Piola-Kirchoff stresses in tensor form
  real(rp),    intent(out) :: gpdds(ndime,ndime,ndime,ndime,pgaus) !< 2nd elasticity tensor in tensor form

  ! --------------------------------------------------------------------------------
  integer(ip)                             :: &
       igaus
  real(rp)                                :: &
       auxMA33(ndime,ndime),                 & ! Auxiliar 2nd Piola-Kirchoff stresses in tensor form
       gpgre(ndime,ndime,pgaus),             & ! Green-Lagrange strains in tensor form
       gprot(ndime,ndime),                   & ! Rotation matrix for material orientations
       vogre(nvoig_sld),                     & ! Strains  in Voigt notation
       vostr(nvoig_sld),                     & ! Stresses in Voigt notation
       stiff(nvoig_sld,nvoig_sld),           & ! Stiffness in Voigt form
       eltan(nvoig_sld,nvoig_sld),           & ! Second elasticity tensor in Voig form
       auxV1(nvoig_sld),                     & ! Auxiliary stuff
       auxMA3333(ndime,ndime,ndime,ndime)
  real(rp)                                 :: E11, E22, E33              !< Young Moduli               
  real(rp)                                 :: v12, v13, v23              !< Poisson ratio               
  real(rp)                                 :: G12, G13, G23              !< Shear moduli        
  real(rp)                                 :: alpha11, alpha22, alpha33  !< Thermal coefficients
  real(rp)                                 :: temp1, temp0
  
  ! =============================================================|    INIT    |=====

  ! ================================================================================
  ! MAIN
  ! --------------------------------------------------------------------------------
  ! INITIALIZE VARIABLES
  !
  ! Get material properties
  !
  call sm151_getproperties(pmate,E11,E22,E33,v12,v13,v23,G12,G13,G23,alpha11,alpha22,alpha33)
  !
  ! Undamaged stiffness tensor
  !
  stiff(:,:) = stiff0_sld(:,:,pmate)  
  !
  !
  ! Initialise
  !
  gpstr = 0.0_rp
  gpdds = 0.0_rp
  ! --------------------------------------------------------------------------------
  ! LOOP OVER GAUSS POINTS
  !
  !...| Do gauss points |...........................................................
  do igaus = 1,pgaus
     !
     ! -----------------------------------------------------------------------------
     ! STRAIN TENSOR
     !
     ! Strain tensor
     !   0 - Infinitesimal tensor
     !   1 - Green-Lagrange tensor
     !
     if (kfl_strai_sld == SLD_INFINITESIMAL) then
        call SM_strain_tensor(0_ip, gpgdi(:,:,igaus), gpgre(:,:,igaus))
     else if (kfl_strai_sld == SLD_GREEN) then
        call SM_strain_tensor(1_ip, gpgdi(:,:,igaus), gpgre(:,:,igaus))
     end if
     !
     ! From tensor to Voigt notation
     !
     call SM_tensor_to_voigt_second(ndime, gpgre(:,:,igaus), vogre(:))
     !
     ! Rotate from the global to the material coordinate system
     !
     auxV1 = vogre
     call SM_rotate_basis_creation(ielem, gpgdi(:,:,igaus), gprot(:,:))
     call SM_rotate_voigt_second(1_ip, gprot(:,:), auxV1(:), vogre(:))
     !
     ! Thermal effects
     !
     temp0 = gptmp(4,igaus) ! Past Time-value or Initial value
     temp1 = gptmp(1,igaus) ! Current iteration value
     !
     vogre(1) = vogre(1) - alpha11*(temp1 - temp0)
     vogre(2) = vogre(2) - alpha22*(temp1 - temp0)
     vogre(3) = vogre(3) - alpha33*(temp1 - temp0)
     !
     ! -----------------------------------------------------------------------------
     ! EFFECTIVE STRESS TENSOR
     !
     call SM_stress_tensor(0_ip, vogre(:), stiff(:,:), vostr(:))
     !
     ! Rotate from the local to the global coordinate system
     !
     auxV1 = vostr
     call SM_rotate_voigt_second(2_ip, gprot(:,:), vostr(:), auxV1(:))
     !
     ! From Voigt notation to tensorial notation
     !
     call SM_voigt_to_tensor_second(ndime, auxMA33(:,:), vostr(:))
     !
     ! Stress tensor due to strain measure
     !
     if (kfl_strai_sld == SLD_INFINITESIMAL) then
        !
        ! Assuming only Infinitesimal
        !
        call SM_stress_transport(SM_PULLBACK, ndime, gpdet(igaus), gpgdi(:,:,igaus), gpidg(:,:,igaus), &
             auxMA33(:,:), gpstr(:,:,igaus))

     else if (kfl_strai_sld == SLD_GREEN) then
        !
        ! Assuming full Green Lagrange
        !
        gpstr(:,:,igaus) = auxMA33(:,:)

     end if

     !
     ! ------------------------------------------------------------------------------
     ! SECOND ELASTICITY TENSOR (dS/dE = C^T)
     !
     ! [dS/dE] = [H]^-1 - [M]
     !
     ! ...| if implicit scheme |......................................................
     if (flagt == 1_ip) then
        !
        ! Rotate from material to global coordinate system
        !
        call SM_rotate_voigt_fourth(2_ip, gprot(:,:), eltan(:,:), stiff(:,:))
        !
        ! From Voigt to tensor
        !
        call SM_tensor_to_voigt_fourth(ndime, 1_ip, auxMA3333(:,:,:,:), eltan(:,:))
        !
        ! Second elasticity tensor due to strain measure
        !
        if (kfl_strai_sld == SLD_INFINITESIMAL) then
           !
           ! Assuming only Infinitesimal
           !
           call SM_stiffness_transport(SM_PULLBACK, ndime, gpgdi(:,:,igaus),gpidg(:,:,igaus), &
                auxMA3333(:,:,:,:), gpdds(:,:,:,:,igaus))

        else if (kfl_strai_sld == SLD_GREEN) then
           !
           ! Assuming full Green Lagrange
           !
           gpdds(:,:,:,:,igaus) = auxMA3333(:,:,:,:)
        end if

     end if
     ! .....................................................| if implicit scheme |...
     
     ! ---------------------------------------------------
     ! STORING VARIABLES
     rmate_sld(ielem) % a(:,:,igaus) = gprot(:,:)
     
  end do

end subroutine sld_stress_model_151

!----------------------------------------------------------------------------
!> @brief   Pre-calculus
!>
!> @details Get material properties and built stiffness matrix for material
!>          model 151
!----------------------------------------------------------------------------
subroutine sm151_precalculus(imate)

  ! ================================================================================
  ! INIT
  ! --------------------------------------------------------------------------------
  use def_parame, only                    :  &
       ip, rp
  use def_solidz, only                    :  &
       stiff0_sld
  ! --------------------------------------------------------------------------------
  implicit none
  ! --------------------------------------------------------------------------------
  integer(ip), intent(in) :: imate             !< Material code
  ! --------------------------------------------------------------------------------
  real(rp)                                :: &
       E11, E22, E33,                        & ! Material properties
       v12, v13, v23,                        &
       G12, G13, G23,                        &
       v21, v31, v32,                        & ! Not read, calculated
       delta, auxS1                            ! Auxiliary variables
  real(rp)                                :: &
       alpha11,alpha22,alpha33
  !
  ! =============================================================|    INIT    |=====

  ! ================================================================================
  ! MAIN
  ! --------------------------------------------------------------------------------
  ! MATERIAL PROPERTIES
  ! 
  call sm151_getproperties(imate,E11,E22,E33,v12,v13,v23,G12,G13,G23,alpha11,alpha22,alpha33)
  !
  v21 = (v12*E22)/E11
  v31 = (v13*E33)/E11
  v32 = (v23*E33)/E22
  !
  ! --------------------------------------------------------------------------------
  ! STIFFNESS TENSOR
  stiff0_sld(:,:,imate) = 0.0_rp
  delta = (1.0_rp - v12*v21 - v23*v32 - v31*v13 - 2.0_rp*v12*v23*v31)/(E11*E22*E33)
  !
  ! Stiff(1,:)
  !
  auxS1 = E22*E33*delta
  stiff0_sld(1,1,imate) = (1.0_rp - v23*v32)/auxS1
  stiff0_sld(1,2,imate) = (v21 + v31*v23)/auxS1
  stiff0_sld(1,3,imate) = (v31 + v21*v32)/auxS1
  !
  ! Stiff(2,:)
  !
  auxS1 = E33*E11*delta
  stiff0_sld(2,1,imate) = (v12 + v13*v32)/auxS1    ! It must be equal to stiff(1,2)
  stiff0_sld(2,2,imate) = (1.0_rp - v31*v13)/auxS1
  stiff0_sld(2,3,imate) = (v32 + v31*v12)/auxS1
  !
  ! Stiff(3,:)
  !
  auxS1 = E11*E22*delta
  stiff0_sld(3,1,imate) = (v13 + v12*v23)/auxS1    ! It must be equal to stiff(1,3)
  stiff0_sld(3,2,imate) = (v23 + v13*v21)/auxS1    ! It must be equal to stiff(2,3)
  stiff0_sld(3,3,imate) = (1.0_rp - v12*v21)/auxS1
  !
  ! Stiff(4,:), stiff(5,:) and stiff(6,:)
  !
  stiff0_sld(4,4,imate) = G23
  stiff0_sld(5,5,imate) = G13
  stiff0_sld(6,6,imate) = G12
  !
  ! ============================================================|    MAIN     |=====

end subroutine sm151_precalculus

!-----------------------------------------------------------------------
!> 
!> @author  gguillam
!> @date    2022-10-05
!> @brief   Get properties for sm151
!> @details Get properties for sm151
!> 
!-----------------------------------------------------------------------

subroutine sm151_getproperties(pmate,E11,E22,E33,v12,v13,v23,G12,G13,G23,alpha11,alpha22,alpha33)

  use def_kintyp, only : ip,rp 
  use def_solidz, only : parco_sld
  
  implicit none
  
  integer(ip), intent(in)            :: pmate                    !< Material code                              
  real(rp),    intent(out)           :: E11, E22, E33            !< Young Moduli               
  real(rp),    intent(out)           :: v12, v13, v23            !< Poisson ratio               
  real(rp),    intent(out)           :: G12, G13, G23            !< Shear moduli        
  real(rp),    intent(out) :: alpha11,alpha22,alpha33  !< Thermal coefficients                       

  E11     = parco_sld(1,pmate) 
  E22     = parco_sld(2,pmate)
  E33     = parco_sld(3,pmate)
  v12     = parco_sld(4,pmate) 
  v13     = parco_sld(5,pmate)
  v23     = parco_sld(6,pmate)
  G12     = parco_sld(7,pmate) 
  G13     = parco_sld(8,pmate)
  G23     = parco_sld(9,pmate)
  alpha11 = parco_sld(10,pmate)
  alpha22 = parco_sld(11,pmate)
  alpha33 = parco_sld(12,pmate)

end subroutine sm151_GetProperties
