!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!---------------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_stress_model_153.f90
!> @author  Adria Quintanas (adria.quintanas@udg.edu)
!> @date    November, 2015
!>          - Subroutine written
!> @brief   Isotropic damage model (Oliver et. al 1990)
!> @details Isotropic damage model (Oliver et. al 1990)
!> @}
!---------------------------------------------------------------------------------

subroutine sld_stress_model_153(pgaus, pmate, gpgdi, gpstr, ielem, flagt, gpdds)

  ! ================================================================================
  ! INIT
  ! --------------------------------------------------------------------------------
  use def_kintyp,                  only       : &
       ip, rp, lg
  use def_master,                  only       : &
       itinn, ittim, modul,                     &
       ITER_K_STATE, TIME_N_STATE
  use def_domain,                  only       : &
       ndime
  use def_solidz,                  only       : &
       stiff0_sld,                              &
       svegm_sld, celen_sld
  use mod_sld_stress_model_comput, only       : &
       SM_strain_tensor,                        &
       SM_stress_tensor,                        &
       SM_tensor_to_voigt_second,               &
       SM_voigt_to_tensor_second,               &
       SM_tensor_to_voigt_fourth

  implicit none

  ! --------------------------------------------------------------------------------
  integer(ip), intent(in)                    :: &
       pgaus,                                   & ! Number of gauss points
       pmate,                                   & ! Current material number
       ielem,                                   & ! Current element number
       flagt                                      ! Integration scheme flag: 0=explicit; 1=implicit

  real(rp),    intent(in)                    :: &
       gpgdi(ndime,ndime,pgaus)                   ! Displacement gradient

  real(rp),    intent(out)                   :: &
       gpstr(ndime,ndime,pgaus,2),              & ! 2nd Piola-Kirchoff stresses tensor
       gpdds(ndime,ndime,ndime,ndime,pgaus)       ! 2nd elasticity tensor

  ! --------------------------------------------------------------------------------
  integer(ip)                                :: &
       igaus, i, j,                             & ! Index
       idSur, idLaw, flagr

  real(rp)                                   :: &
       gpstr_aux(ndime,ndime,pgaus),            & ! 2nd Piola-Kirchoff stresses tensor
       Exx, vxx, Gxx, xT, xC, gF,               & ! Material properties
       cel,                                     & ! Characterised Element length
       phi, r, d,                               & ! Damage variables
       thefo,                                   & ! Thermodynamical force
       disen,                                   & ! Dissipated energy
       gpgre(ndime,ndime,pgaus),                & ! Strains tensor (mechanical)
       vogre(6),                                & ! Strains tensor in Voigt notation
       vostr(6,2),                              & ! Stresses tensor in Voigt notation:
                                !   vostr(:,1) -> effective
                                !   vostr(:,2) -> nominal
       stiff(6,6,2),                            & ! Stiffness tensor in Voigt form:
                                !   stiff(:,1) -> undamaged (precal)
                                !   stiff(:,2) -> damaged
       dddr, drde(6), eltan(6,6),               & ! 2nd elasticity tensor stuff
       auxS1, auxS2,                            &  ! Auxiliary variables
       auxM2(3,3)

  ! =============================================================|    INIT    |=====

  ! ================================================================================
  ! MAIN
  ! --------------------------------------------------------------------------------
  ! INITIALIZE VARIABLES
  !
  ! Get the material properties
  !
  call sm153_get_properties(pmate, Exx, vxx, Gxx, xT, xC, gF, idSur, idLaw)
  !
  ! Characteristics element length (only computed once)
  !
  if (ittim == 1_ip .and. itinn(modul) == 1_ip ) then
     call sm153_get_length(ielem, cel)
     celen_sld(ielem) = cel
  endif
  cel = celen_sld(ielem)
  !
  ! Undamaged stiffness tensor
  !
  stiff(:, :, 1) = stiff0_sld(:, :, pmate)
  !
  ! Initialise
  !
  gpstr_aux = 0.0_rp
  gpdds = 0.0_rp
  ! --------------------------------------------------------------------------------
  ! LOOP OVER GAUSS POINTS
  !
  !...| Do gauss points |...........................................................
  do igaus = 1, pgaus
     !
     ! Initialise variables
     !
     auxS1 = 0.0_rp
     auxS2 = 0.0_rp
     auxM2 = 0.0_rp
     !
     ! ------------------------------------------------------------------------------
     ! STRAIN TENSOR
     !
     ! Strain tensor
     !   0 - Infinitesimal tensor
     !   1 - Green-Lagrange tensor
     !
     call SM_strain_tensor(1_ip, gpgdi(:,:,igaus), gpgre(:,:,igaus))
     !
     ! From tensor to Voigt notation
     !
     call SM_tensor_to_voigt_second(ndime, gpgre(:,:,igaus), vogre(:))
     !
     ! ------------------------------------------------------------------------------
     ! EFFECTIVE STRESS TENSOR
     !
     call SM_stress_tensor(0_ip, vogre(:), stiff(:,:,1), vostr(:,1))
     !
     ! ------------------------------------------------------------------------------
     ! LOADING FAILURE FUNCTION (phi)
     !
     ! Thermodinamical force (Y)
     !
     thefo = 0.0_rp
     do i = 1, 6
        thefo = thefo + vogre(i)*vostr(i, 1)
     enddo
     !
     ! Constant
     !
     select case(idSur)
        !
     case(1_ip)
        !
        ! Id Surface = 1 : Compression equal to tensile strength
        !
        auxS1 = 1.0_rp/xT
        !
     case(2_ip)
        !
        ! Id Surface = 2 : Compression greater than tensile strength
        !
        call runend('SM153: IdSurf 2 hast not been implemented yet')
        !
     case default
        call runend('SM153: No idSur')
     end select
     !
     ! Loading failure function
     !
     phi = auxS1*sqrt(Exx*thefo)
     !
     ! ------------------------------------------------------------------------------
     ! INTERNAL DAMAGE LAWS (r)
     flagr = 0_ip
     r = max(1.0_rp, svegm_sld(ielem)%a(1, igaus, TIME_N_STATE))
     if (phi .gt. r) then
        r = phi
        flagr = 1_ip
     endif
     !
     ! ------------------------------------------------------------------------------
     ! DAMAGE VARIABLES (d)
     d = 0.0_rp
     select case (idLaw)
        !
     case(1_ip)
        !
        ! IdLaw 1, here gF = rmax
        !
        d = (r - 1.0_rp)/(gF - 1.0_rp)
        !
     case(2_ip)
        !
        ! IdLaw 2
        !
        auxS1 = 1.0_rp/((Exx*gF)/(cel*(xT**2.0_rp)) - 0.5_rp)
        if (auxS1 <= 0.001_rp) then
           call runend('SM153: ELEMENT LENGTH TO LARGE TO AVOID SNAPBACK')
        end if
        d = 1.0_rp - (1.0_rp/r)*exp(auxS1*(1.0_rp - r))
        !
     case(3_ip)
        !
        ! IdLaw 3
        !
        auxS1 = (Exx*(xT**2.0_rp)*cel)/(((xT**2.0_rp)*cel) - (2.0_rp*Exx*gF))
        d = 1.0_rp - (auxS1/Exx)+(auxS1/Exx - 1.0_rp)*(1.0_rp/r)
        !
     case default
        call runend('SM153: No idLaw defined')
     end select
     !
     if ( d > 1.0_rp ) then
        d = 1.0_rp
     endif
     !
     ! ------------------------------------------------------------------------------
     ! DAMAGED STIFFNESS TENSOR
     call sm153_stiffness_tensor(stiff(:,:,2), Exx, vxx, Gxx, d)
     !
     ! ------------------------------------------------------------------------------
     ! STRESS TENSOR
     call SM_stress_tensor(0_ip, vogre(:), stiff(:,:,2), vostr(:,2))
     !
     ! From Voigt notation to tensorial notation
     !
     call SM_voigt_to_tensor_second(ndime, gpstr_aux(:,:,igaus), vostr(:,2))
     ! adria, esto te lo agregue yo, no tiene que cambiar nada, es solo para agregar una
     ! columna mas al gpstr que voy a usar luego
     gpstr(1:ndime,1:ndime,igaus,1)= gpstr_aux(1:ndime,1:ndime,igaus)
     !
     ! ------------------------------------------------------------------------------
     ! SECOND ELASTICITY TENSOR (dS/dE = C^T)
     !
     ! [dS/dE] = [H]^-1 - [M]
     !
     ! ...| if implicit scheme |......................................................
     if (flagt .eq. 1_ip) then
        !
        ! [H]^-1
        !
        do i = 1, 6
           do j = 1, 6
              eltan(i, j) = stiff(i, j, 2)
           enddo
        enddo
        !
        ! [H]^-1 - [M]
        !
        if (flagr == 1_ip .and. d < 1.0_rp) then
           !
           ! dd/dr
           !
           dddr = 0.0_rp

           select case (idLaw)
              !
           case(1_ip)
              !
              ! Id law = 1
              !
              dddr = 1.0_rp/(gF - 1.0_rp)
              !
           case(2_ip)
              !
              ! Id law = 2
              !
              auxS2 = Exx*gF/(cel*xT**2.0_rp) - 0.5_rp
              dddr = (auxS2*r + 1.0_rp)/(r**2.0_rp)*exp(auxS2*(1.0_rp - r))
              !
           case(3_ip)
              !
              ! Id law = 3
              !
              dddr = ((1.0_rp/(r**2.0_rp))*(1.0_rp - (auxS1/Exx)))
              !
           case default
              call runend('SM153: No idLaw defined')
           end select
           !
           ! dr/de
           !
           select case(idSur)
              !
           case(1_ip)
              !
              ! Id surface 1
              !
              do i = 1, 6
                 drde(i) = Exx*vostr(i,1)/((xT**2.0_rp)*phi)
              end do
              !
           case(2_ip)
              !
              ! Id surface 2
              !
              call runend('SM153: IdSurf 2 hast not been implemented yet')
              !
           case default
              call runend('SM153: No idLaw defined')
           end select
           !
           ! [H]^-1 - [M]
           !
           do i = 1, 6
              do j = 1, 6
                 eltan(i,j) = eltan(i,j) - dddr*vostr(i,1)*drde(j)
              enddo
           enddo
        end if
        !
        ! From Voigt to tensorial notation
        !
        call SM_tensor_to_voigt_fourth(ndime, 1_ip, gpdds(:,:,:,:,igaus), eltan(:,:))
        !
     endif
     ! .....................................................| if implicit scheme |...
     !
     ! --------------------------------
     ! DISSIPATED ENERGY
     disen = svegm_sld(ielem)%a(3, igaus, 2)
     auxS1 = 0.5_rp*(d - svegm_sld(ielem)%a(2, igaus, TIME_N_STATE))
     do i = 1, 6
        disen = disen + auxS1*vostr(i, 1)*vogre(i)
     end do
     !
     ! ------------------------------------------------------------------------------
     ! STORING STATE VARIABLES
     svegm_sld(ielem)%a(1, igaus, ITER_K_STATE) = r
     svegm_sld(ielem)%a(2, igaus, ITER_K_STATE) = d
     svegm_sld(ielem)%a(3, igaus, ITER_K_STATE) = disen

  enddo
  ! ..........................................................| Do gauss points |...

  ! ============================================================|    MAIN     |=====

end subroutine sld_stress_model_153
! ##################################################################################

! ##################################################################################
! SUBROUTINE SM153_PRECALCULUS(imate)
subroutine sm153_precalculus(imate)

  ! ================================================================================
  ! INIT
  ! --------------------------------------------------------------------------------
  use def_kintyp, only                     :  &
      ip, rp, lg
  use def_solidz, only                     :  &
      stiff0_sld
  !
  ! --------------------------------------------------------------------------------
  implicit none
  ! --------------------------------------------------------------------------------
  integer(ip), intent(in)                  :: &
      imate
  ! -------------------------------------------------------------------------------
  integer(ip)                              :: &
      idSur, idLaw                              ! Material properties
  real(rp)                                 :: & ! Material properties
      Exx, vxx, Gxx, xT, xC, gF
  !
  ! =============================================================|    INIT    |=====

  ! ================================================================================
  ! MAIN
  ! --------------------------------------------------------------------------------
  ! GET PROPERTIES
  !
  call sm153_get_properties(imate, Exx, vxx, Gxx, xT, xC, gF, idSur, idLaw)
  !
  ! --------------------------------------------------------------------------------
  ! UNDAMAGE STIFF TENSOR
  !
  call sm153_stiffness_tensor(stiff0_sld(:, :, imate), Exx, vxx, Gxx, 0.0_rp)
  !
  ! ============================================================|    MAIN     |=====

end subroutine sm153_precalculus
! ##################################################################################

! ##################################################################################
! SUBROUTINE SM153_STIFFNESS_TENSOR
subroutine sm153_stiffness_tensor(stiff, Exx, vxx, Gxx, d)

  ! ================================================================================
  ! INIT
  ! --------------------------------------------------------------------------------
   use def_kintyp, only       :  ip, rp
  ! usa def_domain, only       :  ndime
  ! --------------------------------------------------------------------------------
   implicit none
  ! --------------------------------------------------------------------------------
  real(rp),    intent(in)                 :: &
      Exx, vxx, Gxx, d
  real(rp),    intent(out)                :: &
      stiff(6,6)                               ! Elastic modulus tensor
  !
  ! =============================================================|    INIT    |=====

  ! ================================================================================
  ! MAIN
  !
  ! --------------------------------------------------------------------------------
  ! STIFFNESS TENSOR
  stiff = 0.0_rp
  stiff(1,1) = (1.0_rp - d)*Exx*(1.0_rp - vxx)/(1.0_rp + vxx)/(1.0_rp - 2.0_rp*vxx)
  stiff(2,2) = stiff(1,1)
  stiff(3,3) = stiff(1,1)
  stiff(1,2) = (1.0_rp - d)*Exx*vxx/(1.0_rp + vxx)/(1.0_rp - 2.0_rp*vxx)
  stiff(1,3) = stiff(1,2)
  stiff(2,1) = stiff(1,2)
  stiff(2,3) = stiff(1,2)
  stiff(3,1) = stiff(1,2)
  stiff(3,2) = stiff(1,2)
  stiff(4,4) = (1.0_rp - d)*Gxx
  stiff(5,5) = stiff(4,4)
  stiff(6,6) = stiff(4,4)
  !
  ! ============================================================|    MAIN     |=====

end subroutine sm153_stiffness_tensor
! ##################################################################################

! ##################################################################################
! SUBROUTINE SM153_GET_PROPERTIES
!
subroutine sm153_get_properties(imate, Exx, vxx, Gxx, xT, xC, gF, idSur, idLaw)

  ! ================================================================================
  ! INIT
  ! --------------------------------------------------------------------------------
  use def_kintyp, only                     :  &
      ip, rp, lg
  use def_solidz, only                     :  &
      parco_sld
  ! --------------------------------------------------------------------------------
  implicit none
  ! --------------------------------------------------------------------------------
  integer(ip),  intent(in)                 :: &
     imate
  integer(ip),  intent(out)                :: &
     idSur, idLaw
  real(rp),     intent(out)                :: &
     Exx, vxx, Gxx, xT, xC, gF
  !
  ! ============================================================|    INIT     |=====

  ! ================================================================================
  ! MAIN
  ! --------------------------------------------------------------------------------
  !
  ! Elastic properties
  !
  Exx = parco_sld( 1, imate)
  vxx = parco_sld( 2, imate)
  Gxx = Exx/2.0_rp/(1.0_rp + vxx)
  !
  ! Strength properties
  !
  xT  = parco_sld( 3, imate)
  xC  = parco_sld( 4, imate)
  gF  = parco_sld( 5, imate)
  !
  ! Options of loading surface and law
  !
  idSur = int(parco_sld( 6, imate), ip)
  idLaw = int(parco_sld( 7, imate), ip)
  !
  ! ============================================================|    MAIN     |=====

end subroutine sm153_get_properties

! ##################################################################################
! SUBROUTINE SM153_GET_LENGTH
!
subroutine sm153_get_length(ielem, celen)

  ! ================================================================================
  ! INIT
  ! --------------------------------------------------------------------------------
  use def_kintyp, only                     :  &
      ip, rp, lg
  use def_domain, only                     :  &
      ltype, lnnod, lnods, ndime,             &
      elmar, hnatu, mnode, coord
  ! --------------------------------------------------------------------------------
  implicit none
  ! --------------------------------------------------------------------------------
  integer(ip),  intent(in)                 :: &
      ielem
  real(rp),     intent(out)                :: &
      celen
  ! --------------------------------------------------------------------------------
  integer(ip)                              :: &
      idime,                                  & ! Index
      pelty, pnode, inode, ipoin                ! Element length required indices
  real(rp)                                 :: &
      tragl(9), hleng(3), elcod(ndime,mnode)    ! Element length
  ! ============================================================|    INIT     |=====

  ! ================================================================================
  ! MAIN

  !
  ! Element length
  !
  pelty = ltype(ielem)
  pnode = lnnod(ielem)
  do inode = 1, lnnod(ielem)
     ipoin = lnods(inode, ielem)
     do idime = 1, ndime
        elcod(idime, inode) = coord(idime, ipoin)
     enddo
  enddo
  call elmlen(ndime, pnode, elmar(pelty)%dercg, tragl, elcod, hnatu(pelty), hleng)
  celen = hleng(1)
  ! ============================================================|    MAIN     |=====

end subroutine sm153_get_length
