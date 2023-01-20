!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!---------------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_stress_model_152.f90
!> @author  Adria Quintanas (adria.quintanas@udg.edu)
!> @date    August, 2016
!>          - Subroutine written
!> @author  Gerard Guillamet
!> @date    December, 2017
!>          - Adds strain measures assumptions
!> @date    January, 2019
!>          - Activation viscosity
!>
!> @brief   Damage model for transversally orthotropic materials.
!>          Based on the model of Pere Maimi (pere.maimi@udg.edu)
!>
!> @details Reference:\n
!>
!>          A. Quintanas-Corominas, P. Maimi, E. Casoni, J.A. Mayugo, A. Turon,
!>          M. Vazquez. A 3D transversally isotropic constitutive model for
!>          advanced composites implemented in a high-performance computation code.
!>          European Journal of Mechanics A/solids, 2018\n
!>
!>          You are invited to use the subroutine for academic research       
!>          purposes only. If you are going to use the subroutine for         
!>          industrial purposes, please notice to the authors.                
!>                                                                      
!>          Please cite the previous papers in your work if you are using the subroutine.
!>          If you have any comment/suggestion you are welcome to send it to the authors.
!>          
!>          To share is to improve.
!>
!> @todo    To do list:\n
!>
!>          - Find another way to save the characteristic length
!>          - Add Thermal and hygrothermal effects
!>
!> @}
!---------------------------------------------------------------------------------

module mod_sld_stress_model_152
  !
  ! ============================================================================================
  ! INIT
  ! --------------------------------------------------------------------------------------------
  use def_kintyp, only                                 :  &
       ip, rp, lg
  use def_master
  !
  implicit none
  !
  real(rp), parameter                                 ::  &
       supOne =  1.0000001_rp,                            &
       infOne =  0.9999999_rp,                            &
       supZer =  0.0001000_rp,                            &
       infZer = -0.0000001_rp,                            &
       dlimi  =  1.0000000_rp,                            & ! Damage limit
       dtlimi =  0.9800000_rp                               ! Damage limit for tan. stiff.
  !
  !==============================================================================| init |=======
  !=============================================================================================
  ! PUBLIC / PRIVATE
  ! --------------------------------------------------------------------------------------------
  public                                              ::  &
       sld_stress_model_152,                              &
       sm152_precalculus
  !
  private                                             ::  &
       sm152_model_properties,                            &
       sm152_stiffness_tensor,                            &
       sm152_loading_functions,                           &
       sm152_damage_evolution,                            &
       sm152_damage_state,                                &
       sm152_tangent_analytical,                          &
       sm152_characteristics_length
  !
  !============================================================================| private |======
  !=============================================================================================
  ! CONTAINS
  ! --------------------------------------------------------------------------------------------
  contains
  !
  !----------------------------------------------------------------------------
  !> @brief   Subroutine model properties
  !> @details Get material properties and perform some checks
  !----------------------------------------------------------------------------
  subroutine sm152_model_properties( &
       properties, E11, E22, EKT, EGT, G12, v12, v21, v23,  xT, xC, yT, yC, sL, sT, aT,         &
       bT, nT, nS, nTQ, nSQ, gLT, gLC, g1T, g2T, g2L, hLT1, sLT1, hLT2, sLT2, rLTc, hLC1,       &
       sLC1, hLC2, sLC2, rLCc, hG1, sG1, h61, s61, a11, a22, b11, b22, length, A, B,            &
       flaVi, eta)

    ! -----------------------------------------------------------------------------------------
    ! INIT
    ! -----------------------------------------------------------------------------------------
    !
    implicit none
    !
    real(rp), intent(in)                               :: &
         properties(:),                                   & !< Properties vector
         length                                             !< Characteristic element length

    real(rp), intent(out)                              :: &
         E11, E22, EKT, EGT, G12, v12, v21, v23,          & ! Elastic
         xT, xC, yT, yC, sL, sT,                          & ! Strength
         aT, bT, nT, nS, nTQ, nSQ,                        & ! Loading functions
         gLT, gLC, g1T, g2T, g2L,                         & ! Fracture toughens
         hLT1, hLT2, sLT1, sLT2, rLTc,                    & ! Slopes, ordinates and inflection
         hLC1, hLC2, sLC1, sLC2, rLCc,                    &
         hG1, sG1, h61, s61,                              &
         a11, a22, b11, b22,                              & ! Hygrothermal
         A, B,                                            & ! Adjusting parameters
         eta                                                ! Viscosity

    logical(lg), intent(out)                           :: &
         flaVi

    real(rp)                                           :: &
         fxT, fxC, fgT, fgC, aux1, aux2
    !
    ! --------------------------------------------------------------------------|  INIT  |-----

    ! -----------------------------------------------------------------------------------------
    ! GET PROPERTIES
    ! -----------------------------------------------------------------------------------------
    !
    ! Elastic properties
    !
    E11 = properties(1)
    E22 = properties(2)
    G12 = properties(3)
    v12 = properties(4)
    v23 = properties(5)
    EKT = E22/(2.0_rp - 2.0_rp*v23)
    EGT = E22/(2.0_rp + 2.0_rp*v23)
    v21 = v12*(E22/E11)

    !
    ! Strength properties
    !
    xT  = properties(6)
    fxT = properties(7)
    xC  = properties(8)
    fxC = properties(9)
    aT  = properties(10)
    bT  = properties(11)
    yT  = properties(12)
    yC  = properties(13)
    sL  = properties(14)
    nT  = properties(15)
    nS  = properties(16)
    nTQ = properties(17)
    nSQ = properties(18)
    sT  = (yT*yC*sqrt(1.0_rp + aT)/(yC + yT))

    !
    ! Cohesive law
    !
    gLT  = properties(19)
    fgT  = properties(20)
    hLT1 = (xT**2)/(gLT*fgT*2.0_rp)
    sLT1 = xT
    hLT2 = ((fxT*sLT1)**2)/((1.0_rp + fgT*(fxT**2 - 1.0_rp))*2.0_rp*gLT)
    sLT2 = ((1.0_rp + fgT*(fxT - 1.0_rp))*fxT*xT)/(1.0_rp + fgT*(fxT**2 - 1.0_rp))
    rLTc = ((sLT1 - sLT2)*E11 + (hLT1*sLT2 - hLT2*sLT1)*length)/((hLT1 - hLT2)*length*sLT1)
    !
    gLC  = properties(21)
    fgC  = properties(22)
    hLC1 = (xC**2)/(gLC*fgC*2.0_rp)
    sLC1 = xC
    hLC2 = ((fxC*sLC1)**2)/((1.0_rp + fgC*(fxC**2 - 1.0_rp))*2.0_rp*gLC)
    sLC2 = ((1.0_rp + fgC*(fxC - 1.0_rp))*fxC*xC)/(1.0_rp + fgC*(fxC**2 - 1.0_rp))
    aux1 = 1.0_rp/(((hLC1*sLC2 - hLC2*sLC1)*length)/(E11*(sLC1 - sLC2))+ 1.0_rp)
    aux2 = 2.0_rp*v12*v21 + v23 - 1.0_rp
    rLCc = -((sLC1 - sLC2)*E11 + (hLC1*sLC2 - hLC2*sLC1)*length)/((hLC1 - hLC2)*length*sLC1)*(1.0_rp  &
         /aux2)*(sqrt(aux2**2 - 4.0_rp*v12*v21*aux2*aux1+ (nTQ + 4.0_rp*v12**2)*(v21         &
         *aux1)**2) - nT*v21*aux1)
    !
    g1T  = properties(23)
    g2T  = properties(24)
    sG1  = sT
    hG1  = (sG1**2)/(g2T*2.0_rp)
    !
    g2L  = properties(25)
    s61  = sL
    h61  = (s61**2)/(g2L*2.0_rp)

    !
    ! Thermal effects
    !
    a11  = properties(26)
    a22  = properties(27)
    b11  = properties(28)
    b22  = properties(29)

    !
    ! Adjusting parameters
    !
    A = 111.5_rp ! To assume G1T = 0.277 kJ/m2
    B = EGT/EKT

    !
    ! Viscostiy
    !
    eta = properties(30)
    if( eta > 0.0_rp ) then
        flaVi = .true.
    else
        flaVi = .false.
    end if

    !
    ! Strength properties reduction
    !
    ! - Fiber. Tensile
    !
    if ((2.0_rp*gLT*E11 - length*xT**2) < supZer) then
       print*, 'WARNING: SOLIDZ ELEMENT LENGTH FOR SM152 IS TOO LARGE TO AVOID SNAPBACK LT'
    end if
    !
    ! - Fiber. Compression
    !
    if ((2.0_rp*gLC*E11 - length*xC**2) < supZer) then
       print*, 'WARNING: SOLIDZ ELEMENT LENGTH FOR SM152 IS TOO LARGE TO AVOID SNAPBACK LC'
    end if
    !
    ! - Matrix. Transversal mode II. dG
    !
    if ((2.0_rp*g2T*E22 - length*sT**2) < supZer) then
       print*, 'WARNING: SOLIDZ ELEMENT LENGTH FOR SM152 IS TOO LARGE TO AVOID SNAPBACK TG'
    end if
    !
    ! - Matrix. d6
    !
    if ((2.0_rp*g2L*G12 - length*sL**2) < supZer) then
       print*, 'WARNING: SOLIDZ ELEMENT LENGTH FOR SM152 IS TOO LARGE TO AVOID SNAPBACK T6'
    end if
    !
    ! -----------------------------------------------------------------| GET PROPERTIES |-----
    !
  end subroutine sm152_model_properties

  !
  !----------------------------------------------------------------------------
  !> @brief   Subroutine stiffness tensor
  !> @details Calculation of the stiffness tensor
  !----------------------------------------------------------------------------
  subroutine sm152_stiffness_tensor( &
       stiff, E11, EKT, EGT, G12, v12, d1, dK, dG, d6)

    ! ----------------------------------------------------------------------------------------
    ! INIT
    ! ----------------------------------------------------------------------------------------
    use def_kintyp, only                  :  &
         ip, rp, lg
    ! ----------------------------------------------------------------------------------------
    implicit none
    !
    real(rp),    intent(in)              :: E11        !< Young Modulus fiber direction
    real(rp),    intent(in)              :: EGT        !< Young Modulus
    real(rp),    intent(in)              :: EKT        !< Youn Modulus
    real(rp),    intent(in)              :: G12        !< In-plane shear modulus
    real(rp),    intent(in)              :: v12        !< In-plane poisson ratio
    real(rp),    intent(in)              :: d1         !< Damage variable for fiber tension/compression
    real(rp),    intent(in)              :: dK         !< Damage variable for matrix mode I
    real(rp),    intent(in)              :: dG         !< Damage variable for matrix mode II
    real(rp),    intent(in)              :: d6         !< Damage variable for matrix shear

    real(rp),    intent(out)             :: stiff(6,6) !< Elastic modulus tensor

    real(rp)                             :: &
         auxS1, auxS2                                  ! Auxiliary variables
    !
    ! -------------------------------------------------------------------------|  INIT  |-----

    ! ----------------------------------------------------------------------------------------
    ! STIFNESS TENSOR
    ! ----------------------------------------------------------------------------------------
    !
    ! INITIALIZE VARIABLES
    stiff(:,:) = 0.0_rp

    !
    ! STIFFNESS TENSOR
    auxS1 = E11 + 4.0_rp*EKT*(v12**2)*(d1*(1.0_rp - dK) + dK - 1.0_rp)
    auxS2 = EKT*(1.0_rp - dK) - (4.0_rp*(EKT**2)*(v12**2)*(d1 - 1.0_rp) &
         *(dK - 1.0_rp)**2)/auxS1
    stiff(1,1) = (E11**2)*(1.0_rp - d1)/auxS1
    stiff(1,2) = 2.0_rp*E11*EKT*v12*(d1 - 1.0_rp)*(dK - 1.0_rp)/auxS1
    stiff(1,3) = stiff(1,2)
    stiff(2,1) = stiff(1,2)
    stiff(2,2) =  EGT*(1.0_rp - dG) + auxS2
    stiff(2,3) = -EGT*(1.0_rp - dG) + auxS2
    stiff(3,1) = stiff(1,2)
    stiff(3,2) = stiff(2,3)
    stiff(3,3) = stiff(2,2)
    stiff(4,4) = EGT*(1.0_rp - dG)
    stiff(5,5) = G12*(1.0_rp - d6)
    stiff(6,6) = stiff(5,5)
    !
    ! ----------------------------------------------------------------------------------------
  end subroutine sm152_stiffness_tensor

  !
  !----------------------------------------------------------------------------
  !> @brief   Subroutine loading functions
  !> @details Define loading failure functions
  !>
  !> @param[in] s1 Stress invariant
  !> @param[in] pT Stress invariant
  !> @param[in] tT Stress invariant
  !----------------------------------------------------------------------------
  subroutine sm152_loading_functions( &
       s1, pT, tT, tL, v12, xT, xC, yT, yC, sL, aT, bT, nT, nS, nTQ, nSQ, phiLT, phiLC, phiT)

    ! -----------------------------------------------------------------------------------------
    ! INIT
    ! -----------------------------------------------------------------------------------------
    use def_kintyp, only                                :  &
         ip, rp, lg
    ! -----------------------------------------------------------------------------------------
    implicit none
    !
    real(rp), intent(in)                                :: &
         s1, pT, tT, tL,                                   & !< Stress invariants
         v12,                                              & !< Elastic properties
         xT, xC, yT, yC, sL,                               & !< Strengths properties
         aT, bT, nT, nS, nTQ, nSQ

    real(rp), intent(out)                               :: &
         phiLT, phiLC, phiT                                  !< Loading failure functions

    !
    ! -----------------------------------------------------------------------------------------

    ! -----------------------------------------------------------------------------------------
    ! Loading failure functions
    !
    ! - Longitudinal
    if (s1 > supZer) then
       ! tensile
       phiLT = (s1 - 2.0_rp*v12*pT)/xT
       phiLC = 0.0_rp
    else if (s1 < -xC*0.1_rp) then
       ! compression
       phiLT = 0.0_rp
       phiLC = (sqrt(s1**2 + nTQ*pT**2 + nSQ*tL**2) + nT*pT  + nS*tL)/xC
    else
       !
       phiLT = 0.0_rp
       phiLC = 0.0_rp
    end if
    ! - Transversal
    phiT = sqrt(((yT + yC)/(yT*yC))**2*(tT**2 + aT*pT**2)/(1.0_rp + aT)+ (bT*tL/sL) &
         **2) + ((yC - yT)*pT)/(yT*yC) + (1.0_rp - bT)*tL/sL
    !
    ! -----------------------------------------------------------------------------------------
  end subroutine sm152_loading_functions

  !
  !----------------------------------------------------------------------------
  !> @brief   Subroutine damage evolution
  !> @details Definition of the threshold variables for longitudinal and transverse loadings
  !----------------------------------------------------------------------------
  subroutine sm152_damage_evolution( &
       flagVisco, eta, dtime, phiLT, phiLC, phiT, rLTold, rLCold, rTold, rLTnew, rLCnew, rTnew,   &
       flagLT, flagLC, flagT)

    ! -----------------------------------------------------------------------------------------
    ! INIT
    ! -----------------------------------------------------------------------------------------
    !
    implicit none
    !
    logical(lg), intent(in)                               :: &
         flagVisco

    real(rp),    intent(in)                               :: &
         eta, dtime,                                          & ! Viscosity parameters
         phiLT, phiLC, phiT,                                  & ! Loading state variables
         rLTold, rLCold, rTold

    real(rp),    intent(out)                              :: &
         rLTnew, rLCnew, rTnew                                  ! Damage threshold variables

    logical(lg), intent(out)                              :: &
         flagLT, flagLC, flagT
    !
    ! ----------------------------------------------------------------------|    INIT    |-----

    ! -----------------------------------------------------------------------------------------
    ! THRESHOLD VARIABLES
    ! -----------------------------------------------------------------------------------------

    !
    ! LONGITUDINAL LOADING
    !
    ! Tensile
    if (phiLT > rLTold) then
       flagLT = .true.
       if (flagVisco) then
          rLTnew = rLTold*eta/(eta + dtime) + phiLT*dtime/(eta + dtime)
       else
          rLTnew = phiLT
       end if
    else
       flagLT = .false.
       rLTnew = rLTold
    end if

    !
    ! Compression
    if (phiLC > rLCold) then
       flagLC = .true.
       if (flagVisco) then
          rLCnew = rLCold*eta/(eta + dtime) + phiLC*dtime/(eta + dtime)
       else
          rLCnew = phiLC
       end if
       ! Copling with the longitudinal tensile damage
       if (rLCnew > rLTnew) then
          rLTnew = rLCnew
       end if
    else
       flagLC = .false.
       rLCnew = rLCold
    end if

    !
    ! TRANSVERSE LOADING
    !
    ! Tensile
    if (phiT > rTold) then
       flagT = .true.
       rTnew = phiT
    else
       flagT = .false.
       rTnew = rTold
    end if
    ! -----------------------------------------------|    DAMAGE THRESHOLD VARIABLE      |-----
    return

  end subroutine sm152_damage_evolution

  !
  !----------------------------------------------------------------------------
  !> @brief   Subroutine damage state
  !> @details Definition of each damage state d1, dG, dK and d6
  !----------------------------------------------------------------------------
  subroutine sm152_damage_state( &
       s1, E11, EGT, G12, v12, v21, v23, nT, nTQ, hLT1, sLT1, hLT2, sLT2, rLTc, hLC1, sLC1, hLC2, &
       sLC2, rLCc, hG1, sG1, h61, s61, A, B, length, rLT, rLC, rT, d1, dG, dK, d6)

    ! -----------------------------------------------------------------------------------------
    ! INIT
    ! -----------------------------------------------------------------------------------------
    implicit none
    !
    real(rp), intent(in)                                :: &
         s1,                                               & !< Stress invariants
         E11, EGT, G12, v12, v21, v23,                     & !< Elastic properties
         nT, nTQ,                                          & !< Strength properties
         hLT1, sLT1, hLT2, sLT2, rLTc,                     & !< Softening law properties (slopes and ordinates)
         hLC1, sLC1, hLC2, sLC2, rLCc,                     &
         hG1, sG1, h61, s61,                               &
         A, B,                                             & !< dK adjusting parameters
         length,                                           & !< Characteristic element length
         rLT, rLC, rT                                        !< Loading failure functions and damage threshold variables

    real(rp), intent(out)                               :: &
         d1, dG, dK, d6                                      !< Scalar damage variables

    real(rp)                                            :: &
         hLT, sLT, hLC, sLC, hG, sG, h6, s6,               & !< Auxiliar variables
         aux1, aux2, aux3, aux4, aux5
    !
    ! ----------------------------------------------------------------------|    INIT    |-----

    ! -----------------------------------------------------------------------------------------
    ! DAMAGE STATE VARIABLES
    ! -----------------------------------------------------------------------------------------

    !
    ! LONGITUDINAL DIRECTIOON
    !
    d1 = 0.0_rp
    !
    ! Tensile
    if (s1 > supZer .and. rLT > supOne) then
       !
       ! 1s tram
       if (rLT < rLTc) then
          hLT = hLT1
          sLT = sLT1
          ! 2nd tram
       else
          hLT = hLT2
          sLT = sLT2
       end if
       !
       ! Damage computation
       d1 = min(dlimi, (E11/(E11 - hLT*length))*(1.0_rp - sLT/(rLT*sLT1)))

       !
       ! Compression
    elseif (s1 < infZer .and. rLC > supOne) then
       !
       ! 1st segment
       if (rLC < rLCc) then
          hLC = hLC1
          sLC = sLC1
          !
          ! 2nd segment
       else
          hLC = hLC2
          sLC = sLC2
       end if
       !
       ! Damage computation
       aux1 = 1.0_rp - length*hLC/E11
       aux2 = 2.0_rp*v12*v21 + v23 - 1.0_rp
       aux3 = sLC1/sLC
       aux4 = (nT*v21 + aux1*aux2*aux3*rLC)**2 - (4.0_rp*(v12**2) + nTQ)*(v21**2)
       aux5 = 2.0_rp*aux2*(2.0_rp*v12*v21 - (nT*v21 + aux1*aux2*aux3*rLC)*aux3*rLC)
       d1 = min(dlimi, (-aux5 - sqrt((aux5**2) - 4.0_rp*aux4*((aux2**2)*((aux3*rLC)    &
            **2 - 1.0_rp))))/(2.0_rp*aux4))
    end if

    !
    ! TRANSVERESAL DAMAGES
    !
    dK = 0.0_rp
    dG = 0.0_rp
    d6 = 0.0_rp
    if (rT > supOne) then
       !
       ! dG
       hG = hG1
       sG = sG1
       dG = min(dlimi, (EGT/(EGT - hG*length))*(1.0_rp - sG/(rT*sG1)))

       !
       ! dK
       dK = min(dlimi, 1.0_rp - (1.0_rp - dG)*(B*(rT - 1.0_rp) + A + 1.0_rp)/(rT + A))
       !
       ! d6
       h6 = h61
       s6 = s61
       d6 = min(dlimi, (G12/(G12 - h6*length))*(1.0_rp - s6/(rT*s61)))
    end if
    !
    ! --------------------------------------------------------------|    DAMAGE STATE    |-----

  end subroutine sm152_damage_state

  !
  !----------------------------------------------------------------------------
  !> @brief   Subroutine tangent analytical
  !> @details Calculation of the tangent stiffness matrix analytically
  !----------------------------------------------------------------------------
  subroutine sm152_tangent_analytical( &
       E11, EKT, EGT, G12, v12, v21, v23, xT, xC, yT, yC, sL,                                   &
       aT, bT, nT, nS, nTQ, nSQ, hLT1, sLT1, hLT2, sLT2,                                        &
       rLTc, hLC1, sLC1, hLC2, sLC2, rLCc, hG1, sG1, h61, s61, A, B, length, flaVi, eta, dTime, &
       streff, strnom, stieff, stinom, rLT, rLC, rT, flaLT, flaLC, flaT, d1, dG, dK, d6, tange)

    ! -----------------------------------------------------------------------------------------
    ! INIT
    ! -----------------------------------------------------------------------------------------
    !
    implicit none
    !
    logical(lg), intent(in)                             :: &
         flaLT, flaLC, flaT, flaVi                           !< Flags

    real(rp), intent(in)                                :: &
         E11, EKT, EGT, G12, v12, v21, v23,                & !< Material properties
         xT, xC, yT, yC, sL,                               &
         aT, bT, nT, nS, nTQ, nSQ,                         &
         hLT1, sLT1, hLT2, sLT2, rLTc,                     &
         hLC1, sLC1, hLC2, sLC2, rLCc,                     &
         hG1, sG1, h61, s61,                               &
         A, B,                                             &
         eta, dTime,                                       & !< Viscosity
         length,                                           & !< Element length
         rLT, rLC, rT,                                     & !< Damage threshold variables
         streff(6),                                        & !< Effective stress tensor
         strnom(6),                                        & !< Nominal stress tensor
         stieff(6,6),                                      & !< Effective stiffness tensor
         stinom(6,6)                                         !< Nominal stiffness tensor

    real(rp), intent(in)                                :: &
         d1, dG, dK, d6                                      !< Scalar damage variables

    real(rp), intent(out)                               :: &
         tange(6,6)                                          !< Tangent stiffness matrix

    integer(ip)                                         :: &
         idime

    real(rp)                                            :: &
         s1, pT, tT, tL,                                   & ! Stress invariants
         d1t, dGt, dKt, d6t,                               & ! Damage variables Tang. stiffness
         hLT, sLT, hLC, sLC, hG, sG, h6, s6,               &
         M(6,6), coef, dddr, drds(6), drde(6),             & ! Tangent stiffness tensor
         aux1, aux2, aux3, aux4, aux5, aux6, aux7, aux8,   & ! Auxiliar variables
         aux9
    !
    ! --------------------------------------------------------------------------|  INIT  |-----

    ! -----------------------------------------------------------------------------------------
    ! ANALYTICAL TANGENT STIFFNESS
    !------------------------------------------------------------------------------------------

    !
    ! Invariants
    !
    s1 = streff(1)
    pT = 0.5_rp*(streff(2) + streff(3))
    if (abs(pT) < supZer) then
       pT = 0.0_rp
    end if
    tT = 0.5_rp*sqrt((streff(2) - streff(3))**2 + 4.0_rp*streff(4)**2)
    if (abs(tT) < supZer) then
       tT = 0.0_rp
    end if
    tL = sqrt(streff(5)**2 + streff(6)**2)
    if (abs(tL) < supZer) then
       tL = 0.0_rp
    end if

    !
    ! Damage limitation
    !
    d1t = min(dtlimi, d1)
    dGt = min(dtlimi, dG)
    dKt = min(dtlimi, dK)
    d6t = min(dtlimi, d6)

    !
    ! M(:,:)
    !
    M  = 0.0_rp
    hG = 0.0_rp
    sG = 0.0_rp
    !  - M(1,:)
    if (flaLT) then
       if (d1t < 1.0_rp) then
          ! drds = drL/dseff
          drds(1) = 1.0_rp/xT
          drds(2) = -v12/xT
          drds(3) = -v12/xT
          drds(4) = 0.0_rp
          drds(5) = 0.0_rp
          drds(6) = 0.0_rp
          ! drde = (drL/dseff)*(dseff/de)
          drde(:) = matmul(drds(:),stieff(:,:))
          !  - visocisty
          if (flaVi) then
             drde(:) = drde(:)*(dTime/(eta + dTime))
          end if
          ! dddr = dd1/drL
          if (rLT < rLTc) then
             hLT = hLT1
             sLT = sLT1
          else
             hLT = hLT2
             sLT = sLT2
          end if
          dddr = (E11*sLT)/(xT*(E11 - hLT*length)*(rLT**2))
          ! coef = (s11*dH/dd1)
          coef = strnom(1)/(E11*((1.0_rp - d1t)**2))
          ! M(1,:) = (s11*dH/dd1)*(dd1/drL)*(drL/dseff)*(dseff/de)
          M(1,:) = coef*dddr*drde(:)
       end if
    else if (flaLC) then
       if (d1t < 1.0_rp) then
          ! drds = drL/dseff
          aux1 = sqrt((streff(1)**2) + nTQ*(pT**2) + nSQ*(tL**2))
          drds(:) = 0.0_rp
          if (aux1 > supZer) then
             if (abs(streff(1)) > supZer) then
                drds(1) = streff(1)/(xC*aux1)
             end if
             drds(2) = 0.5_rp*(nT + nTQ*pT/aux1)/xC
             drds(3) = 0.5_rp*(nT + nTQ*pT/aux1)/xC
             if (abs(streff(5)) > supZer) then
                drds(5) = (nS/tL + nSQ/(2.0_rp*aux1))*streff(5)/xC
             end if
             if (abs(streff(6)) > supZer) then
                drds(6) = (nS/tL + nSQ/(2.0_rp*aux1))*streff(6)/xC
             end if
          end if
          ! drde = (drL/dseff)*(dseff/der)
          drde(:) = matmul(drds(:),stieff(:,:))
          if (flaVi) then
             drde(:) = drde(:)*(dTime/(eta + dTime))
          end if
          ! dddr = dd1/drL
          if (rLC <= rLCc) then
             hLC = hLC1
             sLC = sLC1
          else
             hLC = hLC2
             sLC = sLC2
          end if
          aux1 = 1.0_rp - length*hLC/E11
          aux2 = v23 + 2.0_rp*v12*v21 - 1.0_rp
          aux3 = sLC1/sLC
          aux4 = (nT*v21 + aux1*aux2*aux3*rLC)**2 - (4.0_rp*v12**2 + nTQ)*(v21**2)
          aux5 = 2.0_rp*aux2*(2.0_rp*v12*v21 - (nT*v21 + aux1*aux2*aux3*rLC)*aux3*rLC)
          aux6 = (aux2**2)*((aux3*rLC)**2 - 1.0_rp)
          aux7 = 2.0_rp*aux1*aux2*aux3*(aux1*aux2*aux3*rLC + nT*v21)
          aux8 = -2.0_rp*aux2*aux3*(2.0_rp*aux1*aux2*aux3*rLC + nT*v21)
          aux9 = 2.0_rp*((aux2*aux3)**2)*rLC
          dddr = (aux4*(-aux8 + (2.0_rp*aux6*aux7 - aux5*aux8 + 2.0_rp*aux4*aux9)/sqrt(aux5**2 &
               - 4.0_rp*aux4*aux6)) + (aux5 + sqrt(aux5**2 - 4.0_rp*aux4*aux6))*aux7)/(2.0_rp   &
               *aux4**2)
          coef = strnom(1)/(E11*((1.0_rp - d1t)**2))
          ! M(1,:) = (s11*dH/dd1)*(dd1/drL)*(drL/dseff)*(dseff/de)
          M(1,:) = coef*dddr*drde(:)
       end if
    end if
    !  - M(2:6,:)
    if (flaT) then
       ! drds = drT/dsef
       aux1 = ((yT + yC)/(yT*yC))**2
       aux2 = (yC - yT)/(yT*yC)
       aux3 = sqrt(aux1*((tT**2) + aT*(pT**2))/(1.0_rp + aT) + (bT*tL &
            /sL)**2)
       drds(:) = 0.0_rp
       if (aux3 > supZer) then
          if (abs(streff(2)) > supZer .or. abs(streff(3)) > supZer) then
             drds(2) = 0.5_rp*(aux1*(aT*pT + (streff(2) - streff(3))*0.5_rp)/((aT + 1.0_rp)      &
                  *aux3) + aux2)
             drds(3) = 0.5_rp*(aux1*(aT*pT + (streff(3) - streff(2))*0.5_rp)/((aT + 1.0_rp)      &
                  *aux3) + aux2)
          end if
          if (abs(streff(4)) > supZer) then
             drds(4) = (aux1/((1.0_rp + aT)*aux3))*streff(6)
          end if
          if (abs(streff(5)) > supZer) then
             drds(5) = ((bT**2)/(sL*aux3) + ((1.0_rp - bT)/tL))*(streff(5)/sL)
          end if
          if (abs(streff(6)) > supZer) then
             drds(6) = ((bT**2)/(sL*aux3) + ((1.0_rp - bT)/tL))*(streff(6)/sL)
          end if
       end if
       ! drde = (drL/dseff)*(dseff/de)
       drde(:) = matmul(drds(:),stieff(:,:))
       ! - M(2,:)
       ! - M(3,:)
       if (pT > 0.0_rp .and. dGt < 1.0_rp) then
          ! dddr = ddG/drT
          hG = hG1
          sG = sG1
          dddr = (EGT*sG)/(sG1*(EGT - hG*length)*(rT**2))
          ! coef = (s22 - s33)*(dH/ddG)
          coef = (strnom(2) - strnom(3))/(4.0_rp*EGT*((1.0_rp - dGt)**2))
          ! aux1 = coef*ddr
          aux1 = coef*dddr
          !
          ! coef = (s33 + s22)*(dH/ddG)
          coef = (strnom(3) - strnom(2))/(4.0_rp*EGT*((1.0_rp - dGt)**2))
          ! aux1 = coef*dddr
          aux2 = coef*dddr
       else
          aux1 = 0.0_rp
          aux2 = 0.0_rp
       end if
       if (dKt < 1.0_rp) then
          ! dddr = ddk/drT
          dddr = ((A + rT)*(A + B*(rT - 1.0_rp) + 1.0_rp)*(EGT*sG)/(sG1*(EGT - hG*length)    &
               *(rT**2)) + (A + 1.0_rp)*(B - 1.0_rp)*dGt - (A + 1.0_rp)*(B - 1.0_rp))/((A   &
               + rT)**2)
          ! coef = (s22 + s33)*(dH/ddK)
          coef = (strnom(2) + strnom(3))/(4.0_rp*EKT*((1.0_rp - dKt)**2))
          ! aux 3 = coef*dddr
          aux3 = coef*dddr
          ! coef = (s22 + s33)*(dH/ddK)
          coef = (strnom(2) + strnom(3))/(4.0_rp*EKT*((1.0_rp - dKt)**2))
          ! aux 3 = coef*dddr
          aux4 = coef*dddr
       else
          aux3 = 0.0_rp
          aux4 = 0.0_rp
       endif
       ! M(2,:) = ((s22 - s33)*dH/ddG)*ddG/drT + (s22 + s33)*dH/ddK*ddK/drT)*drT/dseff*dseff/de
       M(2,:) = (aux1 + aux3)*drde(:)
       ! M(3,:) = ((s33 - s22)*dH/ddG)*ddG/drT + (s22 + s33)*dH/ddK*ddK/drT)*drT/dseff*dseff/de
       M(3,:) = (aux2 + aux3)*drde(:)
       ! - M(4,:)
       if (dGt < 1.0_rp) then
          ! dddr = ddG/drT
          hG = hG1
          sG = sG1
          dddr = (EGT*sG)/(sG1*(EGT - hG*length)*(rT**2))
          ! coef = (s23*dH/ddG)
          coef = strnom(4)/(EGT*((1.0_rp - dGt)**2))
          ! M(4,:) = (s23*dH/ddG)*(ddG/drT)*(drL/dseff)*(dseff/de)
          M(4,:) = coef*dddr*drde(:)
       end if
       ! - M(5,:)
       ! - M(6,:)
       if (d6t < 1.0_rp) then
          ! dddr = dd6/drT
          h6 = h61
          s6 = s61
          dddr = (G12*s6)/(s61*(G12 - h6*length)*(rT**2))
          !
          ! coef = (s12*dH/dd6)
          coef = strnom(5)/(G12*((1.0_rp - d6t)**2))
          ! M(5,:) = (s12*dH/dd6)*(dd6/drT)*(drL/dseff)*(dseff/de)
          M(5,:) = coef*dddr*drde(:)
          !
          ! coef = (s13*dH/dd6)
          coef = strnom(6)/(G12*((1.0_rp - d6t)**2))
          ! M(6,:) = (s12*dH/dd6)*(dd6/drT)*(drL/dseff)*(dseff/de)
          M(6,:) = coef*dddr*drde(:)
       end if
    end if

    !
    ! (I - M)
    !
    do idime = 1, 6
       M(idime,idime) = 1.0_rp - M(idime,idime)
    end do

    !
    ! DDSDDE = H^(-1) : (I - M)
    !
    tange(:,:) = 0.0_rp
    tange(:,:) = matmul(stinom(:,:),M(:,:))
    !
    ! ============================================================|    MAIN     |========

  end subroutine sm152_tangent_analytical

  !
  !----------------------------------------------------------------------------
  !> @brief   Subroutine characteristic element lenght
  !> @details Calculation element characteristic lenght (minimum)
  !----------------------------------------------------------------------------
  subroutine sm152_characteristics_length( &
       ielem, length)

    ! ================================================================================
    ! INIT
    ! --------------------------------------------------------------------------------
    use def_domain, only                     :  &
         ltype, lnnod, lnods, ndime,             &
         elmar, hnatu, mnode, coord
    !
    implicit none
    !
    integer(ip),  intent(in)                 :: &
         ielem
    real(rp),     intent(out)                :: &
         length
    !
    integer(ip)                              :: &
         idime,                                  & ! Index
         pelty, pnode, inode, ipoin                ! Element length required indices
    real(rp)                                 :: &
         tragl(9), hleng(3), elcod(ndime,mnode)    ! Element length
    !
    ! =============================================================|    INIT    |=====

    ! ================================================================================
    ! ELEMENT LENGTH
    !
    pelty = ltype(ielem)
    pnode = lnnod(ielem)
    !
    ! Element coordinates
    !
    do inode = 1, lnnod(ielem)
       ipoin = lnods(inode, ielem)
       do idime = 1, ndime
          elcod(idime, inode) = coord(idime, ipoin)
       end do
    end do
    ! Element lenght
    call elmlen(ndime, pnode, elmar(pelty)%dercg, tragl, elcod, hnatu(pelty), hleng)
    ! Minimum element lenght
    length = hleng(ndime)
    !
    ! ===================================================|  ELEMENT LENGTH  |=========

  end subroutine sm152_characteristics_length

  !
  !----------------------------------------------------------------------------
  !> @brief   Pre-calculus for stress model 152
  !> @details Get properties and define undamage stiffness tensor
  !----------------------------------------------------------------------------
  subroutine sm152_precalculus( &
       imate)

    ! ================================================================================
    ! INIT
    ! --------------------------------------------------------------------------------
    use def_solidz, only                     :  &
         stiff0_sld, parco_sld
    ! --------------------------------------------------------------------------------
    implicit none
    !
    integer(ip), intent(in)                  :: & !< Material code
         imate

    real(rp)                                 :: & ! Material properties
         E11, E22, EKT, EGT, G12, v12, v23!,    &
    !v21
    ! =============================================================|    INIT    |=====

    ! ================================================================================
    ! MAIN
    ! --------------------------------------------------------------------------------
    ! GET PROPERTIES
    !
    E11 = parco_sld(1,imate)
    E22 = parco_sld(2,imate)
    G12 = parco_sld(3,imate)
    v12 = parco_sld(4,imate)
    v23 = parco_sld(5,imate)
    EKT = E22/(2.0_rp - 2.0_rp*v23)
    EGT = E22/(2.0_rp + 2.0_rp*v23)

    ! --------------------------------------------------------------------------------
    ! UNDAMAGE STIFF TENSOR
    !
    call sm152_stiffness_tensor(stiff0_sld(:, :,imate), E11, EKT, EGT, G12, v12, 0.0_rp, 0.0_rp,   &
         0.0_rp, 0.0_rp)
    ! --------------------------------------------------------------------------------
    ! FIND A & B. INTERPOLATING POLYNOMIAL
    !
    !call sm152_find_A_B(imate, EKT, EGT, v12, v21, v23, aT, yT, yC, g1T, g2T, sT)                  ! <<<< DEBUG >>>>

    ! ============================================================|    MAIN     |=====

  end subroutine sm152_precalculus

  !
  !----------------------------------------------------------------------------
  !> @brief   Main subroutine for stress model 152
  !> @details Get properties and define undamage stiffness tensor
  !----------------------------------------------------------------------------
  subroutine sld_stress_model_152( &
       pgaus, pmate, gpgdi, gpidg, gpdet, gpstr, ielem, flagImpl, gpdds)
    ! -----------------------------------------------------------------------------------------
    ! INIT
    ! -----------------------------------------------------------------------------------------
    use def_kintyp, only                               :  &
         ip, rp, lg
    use def_master, only                               :  &
         dtime, itinn, ittim, modul,                      &
         ITER_K_STATE, TIME_N_STATE
    use def_domain, only                               :  &
         ndime
    use def_solidz, only                               :  &
         stiff0_sld, parco_sld,                           &
         svegm_sld, celen_sld,                            &
         kfl_strai_sld, SLD_GREEN, SLD_INFINITESIMAL,     &
         lawta_sld, SLD_TANGENT_ANALYTICAL,           &
         SLD_TANGENT_NUMERICAL, SLD_SECANT
    use mod_sld_stress_model_comput, only              :  &
         SM_strain_tensor,                                &
         SM_stress_tensor,                                &
         SM_stress_transport,                             &
         SM_stiffness_transport,                          &
         SM_tensor_to_voigt_second,                       &
         SM_voigt_to_tensor_second,                       &
         SM_tensor_to_voigt_fourth,                       &
         SM_rotate_basis_creation,                        &
         SM_rotate_voigt_second,                          &
         SM_rotate_voigt_fourth,                          &
         SM_PULLBACK
    !
    implicit none
    !
    integer(ip), intent(in)                            :: &
         pgaus,                                           & !< Number of gauss points
         pmate,                                           & !< Current material number
         ielem,                                           & !< Current element number
         flagImpl                                           !< Integration scheme flag

    real(rp),    intent(in)                            :: &
         gpgdi(ndime,ndime,pgaus),                        & !< F       : Deformation gradient tensor
         gpidg(ndime,ndime,pgaus),                        & !< F^{1}   : Inverse of the deformation gradient tensor
         gpdet(pgaus)                                       !< J = |F| : Determinant of the deformation gradient tensor

    real(rp),    intent(out)                           :: &
         gpstr(ndime,ndime,pgaus),                        & !< 2nd Piola-Kirchoff stresses tensor
         gpdds(ndime,ndime,ndime,ndime,pgaus)               !< 2nd elasticity tensor

    !
    logical(lg)                                        :: &
         flagLT, flagLC, flagT, flagVisc                     ! Flags (damLT, damLC, damT, visc)

    integer(ip)                                        :: &
         igaus, ivoig

    real(rp)                                           :: &
         E11, E22, EKT, EGT, G12, v12, v21, v23,          & ! Material properties (- elastic)
         xT, xC, yT, yC, sL, sT,                          & !   - Strength
         aT, bT, nT, nS, nTQ, nSQ,                        & !   - Friction and adjust
         gLT, gLC, g1T, g2T, g2L,                         & !   - Fracture toughens
         hLT1, hLT2, sLT1, sLT2, rLTc,                    & !   - Slopes, ordinates and inflection
         hLC1, hLC2, sLC1, sLC2, rLCc,                    &
         hG1, sG1, h61, s61,                              &
         a11, a22, b11, b22,                              & !   - Hygrothermal
         s1, pT, tT, tL, pTnom,                           & ! Stress invariants
         phiLT, phiLC, phiT,                              & ! Loading functions
         rLTold, rLCold, rTold,                           & ! Internal variables
         rLTnew, rLCnew, rTnew,                           &
         d1, dK, dG, d6,                                  & ! Damage variables
         d1aux, dKaux, dGaux, d6aux,                      &
         A, B,                                            & ! Auxiliar damage variables
         length,                                          & ! Element length
         eta,                                             & ! Viscosity
         traMa(3,3),                                      & ! Transformation matrix
         alpha(3,3), beta(3,3),                           & ! Hygrothermal tensor in Voigt form.
         gpgre(ndime,ndime,pgaus),                        & ! Strains tensor (mechanical)
         greGlo(6),                                       & ! Strains tensor in Voigt form (global CSYS)
         greMat(6),                                       & ! Strains tensor in Voigt form (material CSYS)
         strGloNom(6),                                    & ! Nominal stress tensor in Voigt form (global CSYS)
         strMatEff(6),                                    & ! Stress tensor in Voigt form (material CSYS)
         strMatNom(6),                                    & ! Nominal stress tensor in Voigt form (material CSYS)
         stiMatEff(6,6),                                  & ! Effective tiffness tensor in Voigt form (material CSYS)
         stiMatNom(6,6),                                  & ! Nominal stiffness tensor in Voigt form (material CSYS)
         tanGlo(6,6),                                     & ! Material tangent stiffness (global CSYS)
         tanMat(6,6),                                     & ! Material tangent stiffness (material CSYS)
         auxS1, auxS2, auxS3, auxMA33(3,3),               & ! Auxiliar stuff
         auxMA3333(3,3,3,3)
    !
    ! --------------------------------------------------------------------------|  INIT  |-----

    ! -----------------------------------------------------------------------------------------
    ! STRESS TENSOR AND SECOND ELASTICITY TENSOR
    ! -----------------------------------------------------------------------------------------
    ! GET PROPERTIES
    !
    ! Characteristics element length (computed only once)
    !
    if (ittim == 1_ip .and. itinn(modul) == 1_ip ) then
       call sm152_characteristics_length(ielem, length)
       celen_sld(ielem) = length
    endif
    length = celen_sld(ielem)
    !
    ! Material properties
    !
    call sm152_model_properties( &
         parco_sld(:,pmate), E11, E22, EKT, EGT, G12, v12, v21, v23,  xT, xC, yT, yC, sL, sT,   &
         aT, bT, nT, nS, nTQ, nSQ, gLT, gLC, g1T, g2T, g2L, hLT1, sLT1, hLT2, sLT2, rLTc, hLC1, &
         sLC1, hLC2, sLC2, rLCc, hG1, sG1, h61, s61, a11, a22, b11, b22, length, A, B,      &
         flagVisc, eta)

    !
    ! Undamaged stiffness tensor
    !
    stiMatEff(:,:) = stiff0_sld(:,:,pmate)
    !
    ! Thermal expansion tensor  <<<< NOT IMPLEMENTED >>>>
    !
    alpha = 0.0_rp
    alpha(1,1) = a11
    alpha(2,2) = a22
    alpha(3,3) = a22
    !
    ! Moisture tensor  <<<< NOT IMPLEMENTED >>>>
    !
    beta = 0.0_rp
    beta(1,1) = b11
    beta(2,2) = b22
    beta(3,3) = b22
    !
    ! Loop over material points
    !
    gpstr(:,:,:) = 0.0_rp
    gpdds(:,:,:,:,:) = 0.0_rp

    ! ....| material_points |................
    material_points: do igaus = 1, pgaus

       ! ------------------------------------------------------------
       ! STRAIN TENSOR
       !
       greMat(:) = 0.0_rp
       greGlo(:) = 0.0_rp
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
       call SM_tensor_to_voigt_second(ndime, gpgre(:,:,igaus), greGlo(:))
       !
       ! Rotate from the global to the material coordinate system
       !
       call SM_rotate_basis_creation(ielem, gpgdi(:,:,igaus), traMa(:,:))
       call SM_rotate_voigt_second(1_ip, traMa(:,:), greGlo(:), greMat(:))

       ! ------------------------------------------------------------
       ! EFFECTIVE STRESS TENSOR AND INVARIANTS
       !
       strMatEff(:) = 0.0_rp
       !
       ! Effective stress tensor
       !
       strMatEff(:) = matmul(stiMatEff(:,:),greMat(:))
       !
       ! Invariants
       !
       s1 = strMatEff(1)
       pT = 0.5_rp*(strMatEff(2) + strMatEff(3))
       if (abs(pT) < supZer) then
          pT = 0.0_rp
       end if
       tT = 0.5_rp*sqrt((strMatEff(2) - strMatEff(3))**2 + 4.0_rp*strMatEff(4)**2)
       if (abs(tT) < supZer) then
          tT = 0.0_rp
       end if
       tL = sqrt(strMatEff(5)**2 + strMatEff(6)**2)
       if (abs(tL) < supZer) then
          tL = 0.0_rp
       end if

       ! ------------------------------------------------------------
       ! DAMAGE STATE
       !
       ! Loading functions
       !
       call sm152_loading_functions( &
            s1, pT, tT, tL, v12, xT, xC, yT, yC, sL, aT, bT, nT, nS, nTQ, nSQ, phiLT, phiLC,  &
            phiT)
       !
       ! Damage threshold
       !
       rLTold = max(infOne, svegm_sld(ielem)%a(1,igaus,TIME_N_STATE))
       rLCold = max(infOne, svegm_sld(ielem)%a(2,igaus,TIME_N_STATE))
       rTold  = max(infOne, svegm_sld(ielem)%a(3,igaus,TIME_N_STATE))
       !
       call sm152_damage_evolution( &
            flagVisc, eta, dtime, phiLT, phiLC, phiT, rLTold, rLCold, rTold, rLTnew, rLCnew,   &
            rTnew, flagLT, flagLC, flagT)
       !
       ! Damage state
       !
       call sm152_damage_state( &
            s1, E11, EGT, G12, v12, v21, v23, nT, nTQ, hLT1, sLT1, hLT2, sLT2, rLTc, hLC1,     &
            sLC1, hLC2, sLC2, rLCc, hG1, sG1, h61, s61, A, B, length, rLTnew, rLCnew, rTnew,   &
            d1, dG, dK, d6)

       ! ------------------------------------------------------------
       ! NOMINAL STRESS TENSOR
       !
       ! Nominal transversal pressure (pT)
       !
       auxS1 = E11 + 4.0_rp*EKT*(v12**2)*(d1*(1.0_rp - dK) + dK - 1.0_rp)
       auxS2 = 2.0_rp*E11*EKT*v12*(d1 - 1.0_rp)*(dK - 1.0_rp)/auxS1
       auxS3 = 2.0_rp*(EKT*(1.0_rp - dK) - (4.0_rp*(EKT**2)*(v12**2)*(d1 - 1.0_rp)  &
            *(dK - 1.0_rp)**2)/auxS1)
       pTnom = 0.5_rp*(2.0_rp*auxS2*greMat(1) + auxs3*(greMat(2) + greMat(3)))
       if (pTnom < -supZer) then
          dK = 0.0_rp
       end if
       !
       ! Stiffness tensor
       !
       call sm152_stiffness_tensor(stiMatNom(:,:), E11, EKT, EGT, G12, v12, d1, dK, dG,      &
            d6)
       !
       ! Nominal stress tensor (material csys)
       !
       strMatNom(:) = matmul(stiMatNom(:,:),greMat(:))
       do ivoig = 1, 6
          if( abs(strMatEff(ivoig)) > 0.0_rp .and. strMatNom(ivoig) == 0.0_rp ) then
             strMatNom(ivoig) = -1.0E-5_rp
          end if
       end do
       !
       ! Nominal stress tensor (global csys)
       !
       call SM_rotate_voigt_second(2_ip, traMa(:,:), strGloNom(:), strMatNom(:))
       !
       ! From Voigt to tensor
       !
       auxMA33(:,:) = 0.0_rp
       call SM_voigt_to_tensor_second(ndime, auxMA33(:,:), strGloNom(:))
       !
       ! PK2 due to strain measure
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

       ! ---------------------------------------------------
       ! SECOND ELASTICITY TENSOR
       !
       ! ....| Implicit scheme |............
       if (flagImpl .eq. 1_ip) then
          !
          tanGlo(:,:) = 0.0_rp
          tanMat(:,:) = 0.0_rp
          !
          if ( lawta_sld(pmate) == SLD_TANGENT_ANALYTICAL) then
             !
             ! Analytical calculation
             !
             call sm152_tangent_analytical( &
                  E11, EKT, EGT, G12, v12, v21, v23, xT, xC, yT, yC, sL,                          &
                  aT, bT, nT, nS, nTQ, nSQ, hLT1, sLT1, hLT2, sLT2,                               &
                  rLTc, hLC1, sLC1, hLC2, sLC2, rLCc, hG1, sG1, h61, s61, A, B, length, flagVisc, &
                  eta, dtime, strMatEff(:), strMatNom(:), stiMatEff(:,:), stiMatNom(:,:),         &
                  rLTnew, rLCnew, rTnew, flagLT, flagLC, flagT, d1, dG, dK, d6, tanMat(:,:))

          else if ( lawta_sld(pmate) == SLD_TANGENT_NUMERICAL) then
             !
             ! Numerical calculation
             !
             call runend('SLS_SM152: NUMERICAL CALCULATION OF THE TANGENT NOT IMPLEMENTED')

          else if ( lawta_sld(pmate) == SLD_SECANT) then
             !
             ! By using the Secant
             !
             d1aux = min(0.9_rp, d1)
             dKaux = min(0.9_rp, dK)
             dGaux = min(0.9_rp, dG)
             d6aux = min(0.9_rp, d6)
             call sm152_stiffness_tensor(tanMat(:,:), E11, EKT, EGT, G12, v12, d1aux, dKaux, dGaux, d6aux)
          end if
          !
          ! Rotate from material to global coordinate system
          !
          call SM_rotate_voigt_fourth(2_ip, traMa(:,:), tanGlo(:,:), tanMat(:,:))
          !
          ! Voigt to tensor form
          !
          call SM_tensor_to_voigt_fourth(ndime, 1_ip, auxMA3333(:,:,:,:), tanGlo(:,:))
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
       !...........| Implicit scheme |.....

       ! ---------------------------------------------------
       ! STORING STATE VARIABLES
       svegm_sld(ielem)%a(1,igaus,ITER_K_STATE) = rLTnew
       svegm_sld(ielem)%a(2,igaus,ITER_K_STATE) = rLCnew
       svegm_sld(ielem)%a(3,igaus,ITER_K_STATE) = rTnew
       svegm_sld(ielem)%a(4,igaus,ITER_K_STATE) = d1
       svegm_sld(ielem)%a(5,igaus,ITER_K_STATE) = dG
       svegm_sld(ielem)%a(6,igaus,ITER_K_STATE) = dK
       svegm_sld(ielem)%a(7,igaus,ITER_K_STATE) = d6

    end do material_points
    ! .................|material_points |....
  end subroutine sld_stress_model_152

end module mod_sld_stress_model_152

