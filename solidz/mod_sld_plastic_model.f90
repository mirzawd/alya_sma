!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!---------------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_plastic_model.f90
!> @author  Eva Casoni
!> @date    February, 2018
!>          - Subroutine written
!> @date    19 February, 2018
!>          - Adds strain measures assumptions
!> @brief   Isotropic hardening plasticity (bilinear model)
!> @}
!---------------------------------------------------------------------------------

module mod_sld_plastic_model

  use def_kintyp,         only: ip, rp, lg
  use def_master,         only: dtime, ittim, itinn, modul
  use def_master,         only: ITER_K_STATE, TIME_N_STATE
  use def_domain,         only: ndime
  use def_solidz,         only: nvoig_sld
  use def_solidz,         only: SLD_IMPLICIT_SCHEME
  use def_solidz,         only: kfl_strai_sld, SLD_INFINITESIMAL, SLD_GREEN
  use def_solidz,         only: svegm_sld, parco_sld, lmate_sld
  use mod_sld_stress_model_comput

  implicit none
  !
  !===============================================================================
  !===============================================================================
  ! PUBLIC / PRIVATE
  ! -------------------------------------------------------------------------------
  public              :: sld_plastic_model_biso

  private

  !===============================================================================
  !===============================================================================
  ! CONTAINS
  ! -------------------------------------------------------------------------------
contains

  !-------------------------------------------------------------------
  !> @brief  Main subroutine for plasticity model isotropic hardening
  !> @details
  !-------------------------------------------------------------------
  subroutine sld_plastic_model_biso(pgaus,gpgdi,gpidg,gpdet,gpstr,ielem,flagImpl,gpdds)

    implicit none

    integer(ip), intent(in)                            :: &
         pgaus,                                           & !< Number of gauss points
         ielem,                                           & !< Current element number
         flagImpl                                           !< Integration scheme flag

    real(rp),    intent(in)                            :: &
         gpgdi(ndime,ndime,pgaus),                        & !< F       : Deformation gradient tensor
         gpidg(ndime,ndime,pgaus),                        & !< F^{1}   : Inverse of the deformation gradient tensor
         gpdet(pgaus)                                       !< J = |F| : Determinant of the deformation gradient tensor

    real(rp),    intent(out)                           :: &
         gpstr(ndime,ndime,pgaus),                        & !< 2nd Piola-Kirchoff stresses tensor
         gpdds(ndime,ndime,ndime,ndime,pgaus)               !< 2nd elasticity tensor

    integer(ip)                                        :: &
         igaus, idime, jdime,                             & ! Index
         i, j,                                            & !
         pmate                                              !  Current material number

    real(rp)                                           :: &
         E, nu,                                           & ! Material properties (- elastic)
         mu, mu2, lambda,bulk3,                           & !
         Ep,                                              & ! Plastic parameters (young modulus plastic, sigma yield)
         gpgre(ndime,ndime,pgaus),                        & ! Strains tensor (mechanical)
                                ! strain measures
         strain1(6),                                      &
         epsPl(6),                                        & ! Plastic strain in voigt notation (it is a state variable)
         epsEl(6),                                        &
         epsPlEq,                                         & ! Equivalent plastic strain (scalar)
         depsPlEq,                                        & ! Equivalent plastic strain (scalar)
                                ! stress measures
         gpcau(ndime,ndime,pgaus),                        & ! Cauchy stress tensor
         stress(6,pgaus),                                 & ! Sigma in voigt notaion (cauchy)
         sigTr(6),                                        & ! Trial stress
         sigDev(6),                                       & ! Deviatoric stress: sigDev = sigTr-sigHydro*deltaij
         sigVM,                                           & ! Equivalent Von Mises stress (scalar)
         sigHydro,                                        & ! Hydrostatic stress O 1/3*(sig11+sig22+sig33)
         sigYeld,                                         & ! yelding stress
         sigY0,                                           & ! Initial yelding stress
                                ! auxiliar variables
         ftr,                                             & ! Yelding function
                                ! tangent stiffness matrices
         vodds_iso(6,6),                                  & ! Tangent stiffness for the elastic part
         vodds_plas(6,6),                                 & ! Tangent stiffness for the plastic part
         vodds(6,6),                                      & ! Tangent stiffness in voigt notation
         gpdds_aux(ndime,ndime,ndime,ndime,pgaus),        & ! Tangent stiffness tensor 4th order
         flow(6),                                         &
         effg,effg2,effg3,efflam,effhrd                     ! Auxiliar scalar quantities for solution purpose

    !
    ! Reading material properties
    !
    pmate = lmate_sld(ielem)
    E = parco_sld(1,pmate)
    nu = parco_sld(2,pmate)
    sigY0 = parco_sld(3,pmate)
    Ep = parco_sld(4,pmate)
    !
    ! Compute elastic properties
    !
    bulk3 = E/(1.0_rp-2.0_rp*nu)
    mu2 = E/(1.0_rp+nu)
    mu = E/(2.0_rp*(1.0_rp+nu))
    lambda = (bulk3-mu2)/3.0_rp

    !
    ! Initialize variables
    !
    gpdds = 0.0_rp
    gpstr = 0.0_rp

    !----| Loop on gauss points |---------
    material_points: do igaus = 1,pgaus
       !
       ! Strain tensor
       !   0 - Infinitesimal tensor
       !   1 - Green-Lagrange tensor
       !
       if (      kfl_strai_sld == SLD_INFINITESIMAL ) then
          call SM_strain_tensor(0_ip, gpgdi(:,:,igaus), gpgre(:,:,igaus))
       else if ( kfl_strai_sld == SLD_GREEN ) then
          call SM_strain_tensor(1_ip, gpgdi(:,:,igaus), gpgre(:,:,igaus))
       end if
       call SM_tensor_to_voigt_second(ndime, gpgre(:,:,igaus), strain1(:))
       !
       ! Compute isotropic tangent stiffness in voigt notation
       !
       vodds_iso = 0.0_rp
       vodds_iso(1,1) = lambda+mu2
       vodds_iso(2,2) = lambda+mu2
       vodds_iso(3,3) = lambda+mu2
       vodds_iso(4,4) = mu
       vodds_iso(5,5) = mu
       vodds_iso(6,6) = mu
       vodds_iso(1,2) = lambda
       vodds_iso(1,3) = lambda
       vodds_iso(2,3) = lambda
       vodds_iso(2,1) = lambda
       vodds_iso(3,1) = lambda
       vodds_iso(3,2) = lambda
       !
       ! Recover Elastic and plastic state variables
       !
       do i=1,6
          epsPl(i) = svegm_sld(ielem)%a(i,  igaus,TIME_N_STATE)
          epsEl(i) = svegm_sld(ielem)%a(i+6,igaus,TIME_N_STATE)
       end do
       epsPlEq = svegm_sld(ielem)%a(13,igaus,TIME_N_STATE)

       !
       ! Calculate predictor stress and elastic strain (voigt notation)
       ! sigma_trial = C:epsilon_elas = C:(eps_n+1 - epsPl_n)
       sigTr = 0.0_rp
       do i=1,6
          do j=1,6
             sigTr(j) = sigTr(j) + vodds_iso(j,i)*(strain1(i)-epsPl(i))
          end do
          epsEl(i) = epsEl(i) + strain1(i) - epsPl(i)
       end do

       !
       ! Compute deviatoric stress and equivalent von misses stress
       !
       sigVM = (sigTr(1)-sigTr(2))**2 + (sigTr(2)-sigTr(3))**2 + (sigTr(3)-sigTr(1))**2
       do i=4,6
          sigVM = sigVM + 6.0_rp*sigTr(i)**2
       end do
       sigVM = sqrt(sigVM/2.0_rp)

       !
       ! Compute current yeld stress (scalar values) and yelding function
       !
       sigYeld = sigY0 + Ep*epsPlEq
       ftr  = sigVM - sigYeld

       !
       ! Yelding criteria
       !
       if (ftr < 0) then
          !
          ! elastic step
          !
          stress(:,igaus) = sigTr
          vodds = vodds_iso

       else
          !
          ! Plastic step for bilinear isotopric hardening
          !
          sigHydro = (sigTr(1)+sigTr(2)+sigTr(3))/3.0_rp
          sigDev = sigTr - sigHydro
          do i=1,ndime
             flow(i) = sigDev(i)/sigVM
          end do
          do i=4,6
             flow(i) = sigTr(i)/sigVM
          end do
          !
          ! Return mapping iteration: particular case of bilineal law
          !
          depsPlEq = ftr/(3.0_rp*mu+Ep)
          !
          ! Update stresses, plastic strains
          !
          do idime=1,ndime
             stress(idime,igaus) = flow(idime)*sigYeld + sigHydro
             epsPl(idime) = epsPl(idime) + (3.0_rp/2.0_rp)*depsPlEq*flow(idime)
             epsEl(idime) = epsEl(idime) - (3.0_rp/2.0_rp)*depsPlEq*flow(idime)
          end do
          do idime=4,6
             stress(idime,igaus) = flow(idime)*sigYeld
             epsPl(idime) = epsPl(idime) + 3.0_rp*depsPlEq*flow(idime)
             epsEl(idime) = epsEl(idime) - 3.0_rp*depsPlEq*flow(idime)
          end do
          epsPlEq = epsPlEq + depsPlEq
          !
          ! Compute plastic stiffness matrix
          !
          effg = mu*sigYeld/sigVM
          effg2 = 2.0_rp*effg
          effg3 = (3.0_rp/2.0_rp)*effg2
          efflam = (bulk3-effg2)/3.0_rp
          effhrd = 3.0_rp*mu*Ep/(3.0_rp*mu + Ep) - effg3

          vodds_plas = 0.0_rp
          if ( flagImpl == 1 ) then
             do idime=1,ndime
                do jdime=1,ndime
                   vodds_plas(jdime,idime) = efflam
                end do
                vodds_plas(idime,idime) = effg2+efflam
             end do
             do idime=4,6
                vodds_plas(idime,idime) = effg
             end do
             do idime=1,6
                do jdime=1,6
                   vodds_plas(jdime,idime) = vodds_plas(jdime,idime) + effhrd*flow(jdime)*flow(idime)
                end do
             end do
          end if
          vodds = vodds_plas
       end if

       call SM_voigt_to_tensor_second(ndime, gpdds_aux(:,:,:,:,igaus), vodds)
       call SM_voigt_to_tensor_second(ndime, gpcau(:,:,igaus),         stress(:,igaus))
       !
       ! Corrections stress and 2n elasticity tensors in accordance to the strain measure assumptions
       !
       if (      kfl_strai_sld == SLD_INFINITESIMAL ) then
          call SM_stress_transport(SM_PULLBACK, ndime, gpdet(igaus), gpgdi(:,:,igaus), gpidg(:,:,igaus), &
               gpcau(:,:,igaus), gpstr(:,:,igaus))
          call SM_stiffness_transport(SM_PULLBACK, ndime, gpgdi(:,:,igaus),gpidg(:,:,igaus), &
               gpdds_aux(:,:,:,:,igaus), gpdds(:,:,:,:,igaus))
       else if ( kfl_strai_sld == SLD_GREEN ) then
          gpstr(:,:,igaus)     = gpcau(:,:,igaus)
          gpdds(:,:,:,:,igaus) = gpdds_aux(:,:,:,:,igaus)
       end if

       ! ------------------------------------------------------
       ! STORING STATE VARIABLES
       !
       do i = 1,6
          svegm_sld(ielem)%a(i,  igaus,ITER_K_STATE) = epsPl(i)
          svegm_sld(ielem)%a(i+6,igaus,ITER_K_STATE) = epsEl(i)
       end do
       svegm_sld(ielem)%a(13,igaus,ITER_K_STATE) = epsPlEq

    end do material_points

  end subroutine sld_plastic_model_biso

end module mod_sld_plastic_model
