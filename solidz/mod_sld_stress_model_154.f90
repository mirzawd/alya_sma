!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_stress_model_154.f90
!> @author  Adria Quintanas (adria.quintanas@udg.edu)
!> @author  Gerard Guillamet (gerard.guillamet@bsc.es)
!> @date    May, 2019
!>          - Subroutine written
!> @date    September, 2021
!>          - Adding dMax as input for the user
!>          - Fix min/max limits d1,d2,d66
!>          - Postprocess at node point and strenght reduction properties
!>          - Strength reductions is applied for each time step
!> @brief   Damage model for transversally orthotropic materials.
!>
!> @details
!>
!>          State Dependent Variables:\n
!>
!>          SDV01 := rFTnew (Damage threshold variables)
!>          SDV02 := rFCnew
!>          SDV03 := rMTnew
!>          SDV04 := rMCnew
!>          SDV05 := d1 (Damage index fiber)  
!>          SDV06 := d2 (Damage index matrix)
!>          SDV07 := d3 (Damage index coupling matrix/fiber compression)
!>          SDV08 := d4 (Damage index matrix shear)
!>          SDV09 := d5 (Damage index fiber tension)
!>          SDV10 := d6 (Damage index matrix shear + coupling with fiber tension)
!>          SDV11 := grePlast (Plastic strain)
!>          SDV12 := greInela (Elastic strain)
!>         
!>          Strength reduction properties:\n
!>
!>          SRP01: lenchar (characteristic element length (max.))
!>          SRP02: xT
!>          SRP03: xT0
!>          SRP04: xC
!>          SRP05: xC0
!>          SRP06: yT
!>          SRP07: yC
!>          SRP08: sL
!>          
!> @details References:\n
!>          Soto, A. and Gonzalez, E.V. Maimi, P. Mayugo, J.A. Pasquali, P.R. Camanho, P.P.\n
!>          A methodology to simulate low velocity impact and compression after impact in \n
!>          large composite stiffened panels, 2018 \n
!>          Maimi,P., Camanho, P.P., Mayugo, J.A., Davila,C.G. A continuum damage model for \n
!>          composite laminates: Part II - Computational, \n
!>          implementation and validation, 2007 \n
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
!>          - Improve calculation of lenchar
!>          - Add Thermal and hygrothermal effects
!>
!> @}
!-----------------------------------------------------------------------

module mod_sld_stress_model_154
  !
  ! ============================================================================================
  ! INIT
  ! --------------------------------------------------------------------------------------------
  use def_kintyp, only : ip, rp, lg
  
  implicit none
  
  real(rp)               :: d1Max, d2Max, d3Max, d4Max, d5Max, d6Max
  real(rp),    parameter :: zero    = 0.0_rp
  real(rp),    parameter :: zeroapp = 0.001_rp
  real(rp),    parameter :: one     = 1.0_rp             
  real(rp),    parameter :: damMax  = 0.99999_rp    ! Bessa uses 0.999_rp. Warning when using dMax
  real(rp),    parameter :: minG_i  = 0.001_rp      ! Strength reduction limits (Bessa)
  real(rp),    parameter :: minG_j  = 0.0001_rp                                        
  !
  !==============================================================================| init |=======
  !=============================================================================================
  ! PUBLIC / PRIVATE
  ! --------------------------------------------------------------------------------------------
  public                                                  ::  &
       sld_stress_model_154,                                  &
       sm154_precalculus
  !
  private                                                 ::  &
       sm154_model_properties,                                &
       sm154_strain_tensor,                                   &
       sm154_stiffness_tensor,                                &
       sm154_failure_criteria,                                &
       sm154_damage_evolution,                                &
       sm154_damage_state,                                    &
       sm154_characteristics_length

  !
  !============================================================================| private |======
  !=============================================================================================
  ! CONTAINS
  ! --------------------------------------------------------------------------------------------
contains
  !
  ! Get properties
  subroutine sm154_model_properties( &
       props, E11, E22, v12, v21, v23, G12, G23, xT, xTO, xC, xCO, yT, yC, sL, &
       alpha0, gxT, gxTO, gxC, gxCO, gyT, gyC, gsL, gG, flagPlast, sP, kP, flagVisco, eta, &
       flagHygro, a11, a22, incrTemp, b11, b22, incrMois)

    use def_solidz, only : ncoef_sld

    ! -----------------------------------------------------------------------------------------
    ! INIT
    ! -----------------------------------------------------------------------------------------
    !
    implicit none
    !
    real(rp),    intent(in)                               :: &
         props(ncoef_sld)

    logical(lg), intent(out)                              :: &
         flagPlast, flagVisco, flagHygro

    real(rp),    intent(out)                              :: &
         E11, E22, v12, v21, v23, G12, G23,                  & ! Elastic
         xT, xTO, xC, xCO, yT, yC, sL, alpha0,               & ! Strength
         gxT, gxTO, gxC, gxCO, gyT, gyC, gsL, gG,            & ! Fracture thoughness
         sP, kP,                                             & ! Plasticity
         eta,                                                & ! Viscosity
         a11, a22, incrTemp, b11, b22, incrMois                ! Hygrothermal
    !
    ! --------------------------------------------------------------------------|  INIT  |-----

    ! -----------------------------------------------------------------------------------------
    ! GET PROPERTIES
    ! -----------------------------------------------------------------------------------------
    !
    ! Elastic material properties
    !
    E11 = props(1)
    E22 = props(2)
    v12 = props(3)
    v23 = props(4)
    v21 = E22*v12/E11
    G12 = props(5)
    G23 = E22/(1.0_rp + v23)/2.0_rp
    !
    ! Strengths
    !
    xT     = props(6)
    xTO    = props(7)*xT
    xC     = props(8)
    xCO    = props(9)*xC
    yT     = props(10)
    yC     = props(11)
    sL     = props(12)
    alpha0 = props(13)
    !
    ! Fracture toughness
    !
    gxT  = props(14)*props(15)             ! gxT  = fGT*G
    gxTO = (1.0_rp - props(15))*props(14)  ! gxTO =(1 - fG)*G
    gxC  = props(17)*props(16)             ! where fG = props(18) (or props(20) for C)
    gxCO = (1.0_rp - props(17))*props(16)  !        G = props(17) (or props(19) for C)
    gyT  = props(18)
    gyC  = props(19)
    gsL  = props(20)
    gG   = gyT/gsL
    !
    ! Plastic properties
    !
    sP = props(21)
    kP = props(22)
    if (kP > 0.0_rp) then
       flagPlast = .true.
    else
       flagPlast = .false.
    end if
    !
    ! Viscosity properties
    !
    eta = props(23)
    if (eta > 0.0_rp) then
       flagVisco = .true.
    else
       flagVisco = .false.
    end if
    !
    ! Hygrothermal elastic properties
    !
    a11 = props(24)
    a22 = props(25)
    incrTemp = props(26)
    b11 = props(27)
    b22 = props(28)
    incrMois = props(29)
    if (incrTemp > 0.0_rp .or. incrMois > 0.0_rp) then
       flagHygro = .true.
    else
       flagHygro = .false.
    end if
    !
    ! Max. damage variables
    !
    d1Max = abs(props(30))                                     
    d2Max = abs(props(31))                                 
    d3Max = abs(props(32))
    d4Max = abs(props(33))
    d5Max = abs(props(34))
    d6Max = abs(props(35))
    !
    ! -----------------------------------------------------------------| GET PROPERTIES |-----
    !
  end subroutine sm154_model_properties

  !-----------------------------------------------------------------------
  !>
  !> @author  Adria Quintanas
  !> @date    2019-05-15
  !> @brief   Strength reduction
  !> @details Strength reduction
  !>
  !-----------------------------------------------------------------------

  subroutine sm154_strength_reduction( &
       E11, E22, G12, gxT, gxTO, gxC, gXCO, gyT, gyC, gsL, lenchar, ielem, xT, xTO, xC, xCO, yT,   &
       yC, sL)

    use def_solidz, only : srpro_sld
    
    ! -----------------------------------------------------------------------------------------
    ! INIT
    ! -----------------------------------------------------------------------------------------
    !
    implicit none
    !
    integer(ip), intent(in)                               :: &
         ielem
    real(rp), intent(in)                                  :: &
         E11, E22, G12,                                      & ! Elastic
         gxT, gxTO, gxC, gXCO, gyT, gyC, gsL,                & ! Fracture thoughness
         lenchar
    real(rp), intent(inout)                               :: &
         xT, xTO, xC, xCO, yT, yC, sL                          ! Strengths
    real(rp)                                              :: &
         minGA, minGB, minGC, lim, d1ch

    !
    ! -----------------------------------------------------------------------------------------

    ! -----------------------------------------------------------------------------------------
    ! STRENGTH REUDCTION
    ! -----------------------------------------------------------------------------------------
    !
    ! LONGITUDINAL
    !
    ! Tension
    minGA = 2.0_rp*(gxT + gxTO)*E11 - lenchar*(xT**2)
    !
    if (minGA < minG_i) then
       lim = -10.0_rp*E11
       xT  = sqrt(-2.0_rp*(gxT + gxTO)*E11*lim/lenchar/(E11 - lim))
       xTO = xT*gxTO/(gxT + gxTO)
    else
       minGB = 2.0_rp*gXT*E11 - lenchar*xT*(xT - xTO)
       if (minGB < minG_j) then
          xTO = 1.01_rp*(lenchar*(xT**2) - 2.0_rp*gxT*E11)/(lenchar*xT)
       end if
       d1ch = 2.0_rp*gxT*E11/(xTO*xT*lenchar + 2.0_rp*gxT*E11)
       !
       minGC = 2.0_rp*gXTO*E11*(1.0_rp - d1ch) - lenchar*(xTO**2)
       if (minGC < minG_j) then
          xTO = 0.95_rp*sqrt(gxTO*E11*(1.0_rp - d1ch)/lenchar)
       end if
    end if
    !
    ! Compression
    minGA = 2.0_rp*(gxC + gxCO)*E11 - lenchar*(xC**2)
    !
    if (minGA < minG_i) then
       lim = -10.0_rp*E11
       xC  = sqrt(-2.0_rp*(gxC + gxCO)*E11*lim/lenchar/(E11 - lim))
       xCO = xC*gxCO/(gxC + gxCO)
    else
       minGB = 2.0_rp*gXC*E11 - lenchar*xC*(xC - xCO)
       if (minGB < minG_j) then
          xCO = 1.01_rp*(lenchar*(xC**2) - 2.0_rp*gxC*E11)/(lenchar*xC)
       end if
       d1ch = 2.0_rp*gxC*E11/(xCO*xC*lenchar + 2.0_rp*gxC*E11)
       !
       minGC = 2.0_rp*gXCO*E11*(1.0_rp - d1ch) - lenchar*(xCO**2)
       if (minGC < minG_j) then
          xCO = 0.95_rp*sqrt(gxCO*E11*(1.0_rp - d1ch)/lenchar)
       end if
    end if
    !
    ! TRANSVERSE
    !
    ! Tensile
    minGA =  2.0_rp*gyT*E22 - lenchar*(yT**2)
    !
    if (minGA < minG_j) then
       lim = -10.0_rp*E22
       yT = sqrt(-2.0_rp*gyT*E22*lim/lenchar/(E22 - lim))
    endif
    !
    ! Compressive
    minGA =  2.0_rp*gyC*E22 - lenchar*(yC**2)
    !
    if (minGA < minG_j) then
       lim = -10.0_rp*E22
       yC = sqrt(-2.0_rp*gyC*E22*lim/lenchar/(E22 - lim))
    endif
    !
    ! SHEAR
    !
    minGA = 2.0_rp*gsL*G12 - lenchar*(sL**2)
    !
    if (minGA < minG_j) then
       lim = -10.0_rp*G12
       sL = sqrt(-2.0_rp*gsL*G12*lim/lenchar/(G12 - lim))
    end if
    !
    ! Save strength reduction properties 
    !
    srpro_sld(1,ielem) = lenchar
    srpro_sld(2,ielem) = xT
    srpro_sld(3,ielem) = xTO
    srpro_sld(4,ielem) = xC
    srpro_sld(5,ielem) = xCO
    srpro_sld(6,ielem) = yT
    srpro_sld(7,ielem) = yC
    srpro_sld(8,ielem) = sL    
    ! -------------------------------------------------------------|  STRENGTH REDUCTION |-----

  end subroutine sm154_strength_reduction

  !-----------------------------------------------------------------------
  !>
  !> @author  Adria Quintanas
  !> @date    2019-05-15
  !> @brief   Strain tensor
  !> @details Strain tensor
  !>
  !-----------------------------------------------------------------------

  subroutine sm154_strain_tensor( &
       G12, flagHygro, alpha, incTemp, beta, incMois,     &
       sP, kP, flagDamMT, greTotal, greElast, grePlast, greInela, flagPlast)

    ! -----------------------------------------------------------------------------------------
    ! INIT
    ! -----------------------------------------------------------------------------------------
    !
    implicit none
    !
    logical(lg), intent(in)                               :: &
         flagHygro, flagDamMT                                  ! Flags

    real(rp), intent(in)                                  :: &
         G12,                                                & ! Elastic properties
         alpha(3), incTemp, beta(3), incMois,                & ! Hygrothermal properties
         sP, kP                                                ! Plasticity properties
    !
    logical(lg), intent(inout)                            :: &
         flagPlast
    !
    real(rp), intent(inout)                               :: &
         greTotal(6),                                        & ! Total strain tensor
         greElast(6),                                        & ! Elastic strain tensor
         grePlast,                                           & ! Plastic strain
         greInela                                              ! Inelastic strain
    !
    integer(ip)                                           :: &
         idof

    real(rp)                                              :: &
         gre12, ftrial, edgreInela, greTrial
    !
    ! -----------------------------------------------------------------------------------------

    ! -----------------------------------------------------------------------------------------
    ! STRAIN TENSORS
    ! -----------------------------------------------------------------------------------------

    !
    ! TOTAL STRAIN
    !
    ! Hygrothermal effects on longitudinal strains
    ! <AQC> Bessa said that should be only alpha*tempNew
    if (flagHygro) then
       do idof = 1, 3
          greTotal(idof) = greTotal(idof) - alpha(idof)*incTemp - beta(idof)*incMois 
       end do
    end if
    !
    ! PLATIC RESPONSE
    !
    edgreInela = 0.0_rp
    !
    if (flagPlast .and. .not. flagDamMT) then
       gre12 = greTotal(6)
       ftrial = abs(gre12 - grePlast) - sP/G12 - kP*greInela
       if (ftrial >= 0.0_rp) then
          edgreInela = ftrial/(1.0_rp + kP)
          greInela = greInela + edgreInela
          greTrial = gre12 - grePlast
          grePlast = grePlast + edgreInela*sign(1.0_rp, greTrial)
          flagPlast = .true.
       else
          flagPlast = .false.
       end if
    else
       flagPlast = .false.
    end if

    !
    ! ELATIC RESPONSE
    !
    greElast(:) = greTotal(:)
    greElast(6) = greElast(6) - grePlast
    ! -----------------------------------------------------------------------------------------
    !return

  end subroutine sm154_strain_tensor

  !-----------------------------------------------------------------------
  !>
  !> @author  Adria Quintanas
  !> @date    2019-05-15
  !> @brief   Stiffness tensor
  !> @details Stiffness tensor
  !>
  !-----------------------------------------------------------------------

  subroutine sm154_stiffness_tensor(&
       E11, E22, G12, G23, v12, v23, d1, d2, d3, d4, d5, d6, stiff)

    use def_master, only : isnain
    ! -----------------------------------------------------------------------------------------
    ! INIT
    ! -----------------------------------------------------------------------------------
    !
    implicit none
    !
    real(rp),    intent(in)                               :: &
         E11, E22, G12, G23, v12, v23,                       &  ! Material properties
         d1, d2, d3, d4, d5, d6                                 ! Damage variables

    real(rp),    intent(out)                              :: &
         stiff(6,6)                                             ! Elastic modulus tensor
    real(rp)                                              :: &
         H11, H22, H33, H12, H23, delta                         ! Auxiliary variables
    !
    ! --------------------------------------------------------------------------|  INIT  |-----

    ! -----------------------------------------------------------------------------------------
    ! STIFFNESS TENSOR
    ! -----------------------------------------------------------------------------------------
    stiff(:,:) = 0.0_rp
    !
    ! Inversion flexibility matrix (3x3)
    if (d1>=damMax .and. d2<damMax .and. d3<damMax) then
       H22 = 1.0_rp/E22/(1.0_rp - d2)
       H33 = 1.0_rp/E22/(1.0_rp - d3)
       H23 = -V23/E22
       stiff(2,2) =  H33/(H33*H22 - H23**2)
       stiff(3,3) =  H22/(H33*H22 - H23**2)
       stiff(2,3) = -H23/(H33*H22 - H23**2)
    else if (d1<damMax .and. d2>=damMax .and. d3<damMax) then
       H11 = 1.0_rp/E11/(1.0_rp - d1)
       H33 = 1.0_rp/E22/(1.0_rp - d3)
       H12 = -V12/E11
       stiff(1,1) =  H33/(H33*H11 - H12**2)
       stiff(3,3) =  H11/(H33*H11 - H12**2)
       stiff(1,3) = -H12/(H33*H11 - H12**2)
    else if (d1<damMax .and. d2<damMax .and. d3>=damMax) then
       H11 = 1.0_rp/E11/(1.0_rp-d1)
       H22 = 1.0_rp/E22/(1.0_rp-d2)
       H12 = -V12/E11
       stiff(1,1) =  H22/(H22*H11 - H12**2)
       stiff(2,2) =  H11/(H22*H11 - H12**2)
       stiff(1,2) = -H12/(H22*H11 - H12**2)
    else if (d1<damMax .and. d2>=damMax .and. d3>=damMax) then
       stiff(1,1) = E11*(1.0_rp - d1)
    else if (d1>=damMax .and. d2<damMax .and. d3>=damMax) then
       stiff(2,2) = E22*(1.0_rp - d2)
    else if (d1>=damMax .and. d2>=damMax .and. d3<damMax) then
       stiff(3,3) = E22*(1.0_rp - d3)
    else if (d1<damMax .and. d2<damMax .and. d3<damMax) then
       H11 = 1.0_rp/E11/(1.0_rp - d1)
       H22 = 1.0_rp/E22/(1.0_rp - d2)
       H33 = 1.0_rp/E22/(1.0_rp - d3)
       H12 = -V12/E11
       H23 = -V23/E22
       delta = H11*(H23**2 - H22*H33) + (H22 + H33-2.0_rp*H23)*(H12**2)
       stiff(1,1) = (H23**2 - H22*H33)/delta
       stiff(1,2) =  H12*(H33 - H23)/delta
       stiff(1,3) =  H12*(H22 - H23)/delta
       stiff(2,2) = (H12**2 - H11*H33)/delta
       stiff(2,3) = (H11*H23 - H12**2)/delta
       stiff(3,3) = (H12**2 - H11*H22)/delta
    end if
    stiff(2,1) = stiff(1,2)
    stiff(3,1) = stiff(1,3)
    stiff(3,2) = stiff(2,3)
    !
    ! Inversion of shear contributions
    stiff(4,4) = (1.0_rp - d4)*G23
    stiff(5,5) = (1.0_rp - d5)*G12
    stiff(6,6) = (1.0_rp - d6)*G12
    !
    ! --------------------------------------------------------------|  STIFFNESS TENSOR  |-----

  end subroutine sm154_stiffness_tensor

  !-----------------------------------------------------------------------
  !>
  !> @author  Adria Quintanas
  !> @date    2019-05-15
  !> @brief   Failure criteria
  !> @details Failure criteria
  !>
  !-----------------------------------------------------------------------

  subroutine sm154_failure_criteria( &
       eps11, eps22, sig11, sig22, sig12, E11, E22, xT, xC, yT, yC, nL, nT, gG, sL, sT, alpha0,   &
       phiFT, phiFC, phiMT, phiMC)

    ! -----------------------------------------------------------------------------------------
    ! INIT
    ! -----------------------------------------------------------------------------------------
    implicit none
    !
    real(rp),    intent(in)                               :: &
         eps11, eps22,                                       & ! Longitudinal strain
         sig11, sig22, sig12,                                & ! Effective stresses
         E11, E22,                                           & ! Elastic material properties
         xT, xC, yT, yC, nL, nT, sL, sT,                     & ! Strength properties
         gG,                                                 & ! Fracture toughness ratio
         alpha0                                                ! Fracture angle

    real(rp),    intent(out)                              :: &
         phiFT,                                              & ! Fiber tensile
         phiFC,                                              & ! Fiber compression
         phiMT,                                              & ! Matrix tensile
         phiMC                                                 ! Matrix compression

    real(rp)                                              :: &
         auxs1, str12abs,                                    & ! Auxiliary variables
         phiC, cosPhiC, sinPhiC, sinAlpha, cosAlpha,         &
         theta, s22T, s12T, teffT, teffL
    !
    ! --------------------------------------------------------------------------|  INIT  |-----

    ! -----------------------------------------------------------------------------------------
    ! LARC03 & LARC04 FAILURE CRITERIA
    ! -----------------------------------------------------------------------------------------
    phiFT = 0.0_rp
    phiFC = 0.0_rp
    phiMT = 0.0_rp
    phiMC = 0.0_rp
    str12abs = abs(sig12)
    !
    ! FIBER DAMAGE
    !
    !
    ! Compressive case
    if (sig11 < zero) then
       ! Larc04 criteria
       auxs1 = abs(sL/xC) + nL
       phiC = atan((1.0_rp - sqrt(max(0.0_rp, 1.0_rp - 4.0_rp*auxs1*sL/xC)))/(2.0_rp*auxs1))
       cosPhiC = cos(phiC)
       sinPhiC = sin(phiC)
       s22T =  sig11*(sinPhiC**2) + sig22*(cosPhiC**2) - 2.0_rp*str12abs*sinPhiC*cosPhiC
       s12T = (sig22 - sig11)*sinPhiC*cosPhic + str12abs*(cosPhiC**2 - sinPhiC**2)
       !
       phiFC = min(max(0.0_rp, abs(s12T) + nL*s22T)/sL, abs(eps11)*E11/yC)
       phiFT = 0.0_rp
       !
       ! Tensile case
    else
       ! Larc03 criteria: max strain
       phiFT = eps11*E11/xT
       phiFC = 0.0_rp

    end if

    !
    ! MATRIX DAMAGE
    !
    ! Compressive case
    if (sig22 < -zeroapp) then
       ! Larc03
       sinAlpha = sin(alpha0)
       cosAlpha = cos(alpha0)
       theta = atan(-str12abs/(sig22*sinAlpha))
       teffT = max(0.0_rp, sig22*cosAlpha*(nT*cosAlpha*cos(theta) - sinAlpha))
       teffL = max(0.0_rp, cosAlpha*(str12abs + nL*cosAlpha*sin(theta)*sig22))
       !
       phiMC = min(sqrt((teffT/sT)**2 + (teffL/sL)**2), E22*(-eps22)/yT)
       phiMT = max(0.0_rp, str12abs + nL*sig22)/sL
       !
       ! Tensile case
    else
       ! Larc03
       phiMT = (1.0_rp - gG)*(sig22/yT) + gG*(sig22/yT)**2 + (sig12/sL)**2
       if (phiMT > zero) then
          phiMT = sqrt(phiMT)
          phiMC = 0.0_rp
       else
          phiMT = 0.0_rp
          phiMC = 0.0_rp
       end if
    end if
    ! ---------------------------------------------------------------------|  LARC03-04  |-----

  end subroutine sm154_failure_criteria

  !-----------------------------------------------------------------------
  !>
  !> @author  Adria Quintanas
  !> @date    2019-05-15
  !> @brief   Damage evolution
  !> @details Damage evolution
  !>
  !-----------------------------------------------------------------------

  subroutine sm154_damage_evolution( &
       flagVisco, eta, dt, phiFT, phiFC, phiMT, phiMC, rFTold, rFCold, rMTold, rMCold, rFTnew, &
       rFCnew, rMTnew, rMCnew, flagFT, flagFC, flagMT, flagMC)

    ! -----------------------------------------------------------------------------------------
    ! INIT
    ! -----------------------------------------------------------------------------------------
    !
    implicit none
    !
    logical(lg), intent(in)                               :: &
         flagVisco

    real(rp),    intent(in)                               :: &
         eta, dt,                                            & ! Viscosity parameters
         phiFT, phiFC, phiMT, phiMC,                         & ! Loading state variables
         rFTold, rFCold, rMTold, rMCold

    real(rp),    intent(out)                              :: &
         rFTnew, rFCnew, rMTnew, rMCnew                         ! Damage threshold variables

    logical(lg), intent(out)                              :: &
         flagFT, flagFC, flagMT, flagMC
    !
    ! --------------------------------------------------------------------------|  INIT  |-----

    ! -----------------------------------------------------------------------------------------
    ! THRESHOLD VARIABLES
    ! -----------------------------------------------------------------------------------------

    !
    ! LONGITUDINAL LOADING
    !
    ! Tensile
    if (phiFT > rFTold) then
       flagFT = .true.
       if (flagVisco) then
          rFTnew = rFTold*eta/(eta + dt) + phiFT*dt/(eta + dt)
       else
          rFTnew = phiFT
       end if
    else
       flagFT = .false.
       rFTnew = rFTold
    end if

    !
    ! Compresssion
    if (phiFC > rFCold) then
       flagFC = .true.
       if (flagVisco) then
          rFCnew = rFCold*eta/(eta + dt) + phiFC*dt/(eta + dt)
       else
          rFCnew = phiFC
       end if
       ! Coupling with the longitudinal tensile damage
       if (rFCnew > rFTnew) then
          rFTnew = rFCnew
       end if
    else
       flagFC = .false.
       rFCnew = rFCold
    end if

    !
    ! TRANSVERSE LOADING
    !
    ! Tensile
    if (phiMT>rMTold) then
       flagMT = .true.
       rMTnew = phiMT
    else
       flagMT = .false.
       rMTnew = rMTold
    end if

    !
    ! Compressive (transversal damage at 53º)
    if (phiMC>rMCold) then
       flagMC = .true.
       rMCnew = phiMC
       ! Coupling with the transverse tensile damage
       if (rMCnew > rMTnew) then
          rMTnew = rMCnew
       end if
    else
       flagMC = .false.
       rMCnew = rMCold
    end if
    
    ! -----------------------------------------------------|  DAMAGE TRHESHOLD VARIABLE  |-----

  end subroutine sm154_damage_evolution

  !
  ! Damage state
  subroutine sm154_damage_state( &
       E11, E22, v12, v21, G12, xT, xTO, xC, xCO, yT, yC, sL, sT, nL, nT, gxT, gxTO, gxC, gxCO,   &
       gyT, gyC, gsL, gG, alpha0, kP, flagPlast, lenchar, eps11, eps22, rFT, rFC, rMT, rMC,       &
       d1, d2, d3, d4, d5, d6)

    ! -----------------------------------------------------------------------------------------
    ! INIT
    ! -----------------------------------------------------------------------------------
    !
    implicit none
    !
    logical(lg), intent(in)                                :: &
         flagPlast

    real(rp), intent(in)                                   :: &
         E11, E22, v12, v21, G12,                             & ! Elastic properties
         xT, xTO, xC, xCO, yT, yC, sL, sT, nL, nT,            & ! Strength properties
         gxT, gxTO, gxC, gxCO, gyT, gyC, gsL, gG,             & ! Energy properties
         alpha0,                                              & ! Fracture angle
         kP,                                                  & ! Plastic property
         lenchar,                                             & ! Characteristic element length
         eps11, eps22,                                        & ! Longitudinal and transversal strain
         rFT, rFC, rMT, rMC                                     ! Damage threshold variables

    real(rp), intent(out)                                  :: &
         d1, d2, d3, d4, d5, d6                                 ! Damage state variables

    real(rp)                                               :: &
         dFTch, rFTch, mxT, nxT,                              & ! Longitudinal direction tensile
         dFCch, rFCch, mxC, nxC,                              & ! Longtindunal direction compression
         myT, nyT,                                            & ! Transverse direction tensil
         myC, nyC,                                            & ! Tranvserse direction compression
         msL, nsL,                                            & ! Shear direction
         s1, s2, s3, a1, a2, a3, auxs1, invGreEff,            & ! Auxiliar parameters
         AYC, BYC,                                            &
         d1T, d1C, d2T, d2C, d66,                             & ! Temporal damage variables
         phiC, cosPhiC, sinPhiC                                 ! Auxilair material properties

    !
    ! --------------------------------------------------------------------------|  INIT  |-----

    ! -----------------------------------------------------------------------------------------
    ! DAMAGE STATE
    ! -----------------------------------------------------------------------------------------
    !
    ! LONGITUDINAL DIRECTION
    !
    ! Tensile
    d1T = 0.0_rp
    if (rFT > one) then
       dFTch = 2.0_rp*gXT*E11/(xT*xTO*lenchar + 2.0_rp*gXT*E11)
       rFTch = xTO/(xT*(1.0_rp - dFTch))
       !
       ! First tram
       if (rFT <= rFTch ) then
          mxT = lenchar*E11*xT*(xT - xTO)/(lenchar*xT*(xT - xTO) - 2.0_rp*gXT*E11)
          nxT = xT*(1.0_rp - mxT/E11)
          !
          ! Second tram
       else
          mxT = lenchar*E11*(1.0_rp - dFTch)*(xTO**2)/(lenchar*(xTO**2) - &
               2.0_rp*gXTO*E11*(1.0_rp - dFTch))
          nxT = xTO*(1.0_rp - mxT/(E11*(1.0_rp - dFTch)))
       end if
       !
       ! Damage variable
       d1T = 1.0_rp - mxT/E11 - nxT/(XT*rFT)
       d1T = max(0.0_rp,d1T)
       d1T = min(1.0_rp,d1T)
    end if
    !
    ! Compression
    d1C = 0.0_rp
    if (rFC > one) then
       dFCch  = 2.0_rp*gXC*E11/(xC*xCO*lenchar + 2.0_rp*gXC*E11)
       !
       mxC = lenchar*xC*E11*(xC - xCO)/(lenchar*xC*(xC - xCO) - 2.0_rp*gXC*E11)
       nxC = XC*(1.0_rp - mxC/E11)
       !
       auxs1 = abs(sL/xC) + nL
       phiC = atan((1.0_rp - sqrt(max(0.0_rp, 1.0_rp - 4.0_rp*auxs1*sL/xC)))/(2.0_rp*auxs1))
       cosPhiC = cos(phiC)
       sinPhiC = sin(phiC)
       !
       s1 = E11 - mxC*v12*v21
       s2 = (E11 - mxC)*v21
       a1 = (s2 - s1)*sinPhiC*cosPhiC + nL*(s1*(sinPhiC**2) + s2*(cosPhiC)**2)
       a2 = nxC*v21*((1.0_rp - v12)*sinPhiC*cosPhiC + nL*((cosPhiC)**2 + v12*(sinPhiC**2)))
       rFCch = ((-xCO/E11/(1.0_rp - dFCch))*a1 + a2)/((1.0_rp - v12*v21)*sL)
       !
       ! First tram
       if (rFC <= rFCch ) then
          invGreEff = a1/(E11*(sL*(1.0_rp - v12*v21)*rFC - a2))
          !
          ! Second tram
       else
          mxC = lenchar*E11*(1.0_rp - dFCch)*(xCO**2)/(lenchar*(xCO**2) - 2.0_rp*gXCO*E11*(1.0_rp - dFCch))
          nxC = xCO*(1.0_rp - mxC/(E11*(1.0_rp - dFCch)))
          !
          s1 = E11 - mxC*v12*v21
          s2 = (E11 - mxC)*v21
          a1 = (s2 - s1)*sinPhiC*cosPhiC + nL*(s1*(sinPhiC**2) + s2*(cosPhiC)**2)
          a2 = nxC*v21*((1.0_rp - v12)*sinPhiC*cosPhiC + nL*((cosPhiC)**2 + v12*(sinPhiC**2)))
          invGreEff = a1/(E11*(sL*(1.0_rp - v12*v21)*rFC - a2))
       end if
       !
       ! Damage variable
       d1C = 1.0_rp - mxC/E11 + nxC*invGreEff
       d1C = max(0.0_rp,d1C)
       d1C = min(1.0_rp,d1C)
    end if
    !
    ! TRANSVERSE DIRECTION
    !
    ! Tensile
    d2T = 0.0_rp
    if (rMT > one) then
       myT = lenchar*(yT**2)*E22/(lenchar*(yT**2) - 2.0_rp*gyT*E22)
       nyT = yT*(1.0_rp - myT/E22)
       !
       s1 = (E22 - myT*v12*v21)/yT/(1.0_rp - v12*v21)
       s2 = nyT*v12*v21/yT/(1.0_rp - v12*v21)
       s3 = s2*(s2*gG + 1.0_rp - gG)
       !
       a1 = gG*(s1**2)
       a2 = s1*(gG - 1.0_rp - 2.0_rp*gG*s2)
       a3 = sqrt((a2**2) - 4.0_rp*a1*(s3 - rMT**2))
       !
       ! Damage variable
       d2T = 1.0_rp - myT/E22 - (2.0_rp*a1*nyT)/(E22*(a3 + a2))
       d2T = max(0.0_rp,d2T)
       d2T = min(1.0_rp,d2T)
    end if
    !
    ! Compressive
    d2C = 0.0_rp
    if (rMC > one) then
       myC = lenchar*(yC**2)*E22/(lenchar*(yC**2) - 2.0_rp*gyC*E22)
       nyC = yC*(1.0_rp - myC/E22)
       !
       AYC = -nyC*v12*v21
       BYC = (v12*v21 - 1.0_rp)/(cos(alpha0)*(sin(alpha0) - nT*cos(alpha0)))
       !
       ! Damage variable
       d2C = 1.0_rp - myC/E22 + nyC*(E22 - myC*v12*v21)/(E22*(AYC + BYC*rMC*sT))
       d2C = max(0.0_rp,d2C)
       d2C = min(1.0_rp,d2C)
    end if
    !
    ! SHEAR DIRECTION
    !
    d66 = 0.0_rp
    if (rMT > one) then
       msL = lenchar*(sL**2)*G12/(lenchar*(sL**2) - 2.0_rp*gsL*G12)
       if (flagPlast) then
          msL = msL*(1.0_rp + kP)/kP
       end if
       nsL = sL*(1.0_rp - msL/G12)
       !
       ! Damage variable
       d66 = 1.0_rp - msL/G12 - nsL/(sL*rMT)
       d66 = max(0.0_rp,d66)
       d66 = min(1.0_rp,d66)
    end if
    !
    ! ACTIVE DAMAGE VARIABLES
    !
    ! Longitudinal direction
    if (eps11 > (d2C - 1.0_rp)*v21*eps22) then
       d1 = d1T
    else if (eps11 < (d2T - 1.0_rp)*v21*eps22) then
       d1 = d1C
    else
       d1 = 0.5_rp*(d1T + d1C)
    end if
    !
    ! Transverse direction
    if (eps22 > (d1 - 1.0_rp)*v12*eps11) then
       d2 = d2T
    else
       d2 = d2C
    end if
    !
    ! DEPENDEND DIRECTIONS
    !
    d3 = 1.0_rp - (1.0_rp - d2C)*(1.0_rp - d1C)
    d4 = d66
    d5 = d1T
    d6 = 1.0_rp - (1.0_rp - d66)*(1.0_rp - d1T)
    !
    ! Limit damage variables (Max. values)
    !
    d1 = min(d1Max,d1)
    d2 = min(d2Max,d2)
    d3 = min(d3Max,d3)
    d4 = min(d4Max,d4)
    d5 = min(d5Max,d5)
    d6 = min(d6Max,d6)
    ! -----------------------------------------------------------------------------------------

  end subroutine sm154_damage_state

  !-----------------------------------------------------------------------
  !>
  !> @author  Adria Quintanas
  !> @date    2019-05-15
  !> @brief   Characteristic length
  !> @details Characteristic length
  !> @note    We should have a lenght for each damage mode
  !>
  !-----------------------------------------------------------------------

  subroutine sm154_characteristics_length( &
       ielem, length)

    ! ================================================================================
    ! INIT
    ! --------------------------------------------------------------------------------
    use def_domain, only                     :  &
         ltype, lnnod, lnods,                   &
         elmar, mnode, coord, lorde,            &
         ngaus, mgaus
    use def_domain, only : hnatu
    !
    implicit none
    !
    integer(ip),  intent(in)                 :: &
         ielem
    real(rp),     intent(out)                :: &
         length
    !
    integer(ip)                              :: &
         pelty, pnode, porde, inode, ipoin,     &  ! Element length required indices
         igaus, pgaus
    real(rp)                                 :: &
         elcod(3,mnode),                        &  ! Element length
         gpcar(3,mnode,mgaus),                  &
         xjaci(3,3),xjacm(3,3)
    real(rp)                                 :: &
         tragl(9), hleng(3)
    real(rp)                                 :: &
         gpdet,evolu
    !
    ! =============================================================|    INIT    |=====

    ! ================================================================================
    ! ELEMENT LENGTH
    !
    pelty = ltype(ielem)
    pnode = lnnod(ielem)
    porde = lorde(pelty)
    pgaus = ngaus(pelty)
    !
    ! Element coordinates
    !
    do inode = 1,pnode
       ipoin = lnods(inode, ielem)
       elcod(1:3, inode) = coord(1:3, ipoin)
    end do
    !
    ! Element volum
    !
    evolu = 0.0_rp
    do igaus = 1,pgaus
       call elmder(pnode,3_ip,elmar(pelty)%deriv(1,1,igaus),&  ! Cartesian derivative
                   elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)    ! and Jacobian
       evolu = evolu + elmar(pelty)%weigp(igaus)*gpdet          ! dV:=|J|*wg
    end do
    !
    ! Element length
    !
    call elmlen(3_ip,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),hleng)
    !
    ! Characteristic element length
    !
    ! Option 1:
    !length = evolu**(1.0_rp/3.0_rp)
    ! Option 2:
    length = evolu/(hleng(1)*hleng(1)) ! It should be divided by the largest face area!
    ! Save lenchar
    length = length/real(porde,rp)
    ! ===================================================|  ELEMENT LENGTH  |=========

  end subroutine sm154_characteristics_length

  !-----------------------------------------------------------------------
  !>
  !> @author  Gerard Guillamet
  !> @date    2019-05-15
  !> @brief   Pre-calculus for stress model 154
  !> @details Get properties and define undamage stiffness tensor
  !>
  !-----------------------------------------------------------------------

  subroutine sm154_precalculus( &
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
         E11, E22, v12, v23, v21, G12, G23

    ! =============================================================|    INIT    |=====

    ! ================================================================================
    ! MAIN
    ! --------------------------------------------------------------------------------
    ! GET PROPERTIES
    !
    E11 = parco_sld(1,imate)
    E22 = parco_sld(2,imate)
    v12 = parco_sld(3,imate)
    v23 = parco_sld(4,imate)
    v21 = E22*v12/E11
    G12 = parco_sld(5,imate)
    G23 = E22/(1.0_rp + v23)/2.0_rp

    ! --------------------------------------------------------------------------------
    ! UNDAMAGE STIFF TENSOR
    !
    call sm154_stiffness_tensor(E11, E22, G12, G23, v12, v23, &
         0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, stiff0_sld(:,:,imate))

  end subroutine sm154_precalculus

  !-----------------------------------------------------------------------
  !>
  !> @author  Adria Quintanas
  !> @date    2019-05-15
  !> @brief   Alya stress model
  !> @details Alya stress model
  !>
  !-----------------------------------------------------------------------

  subroutine sld_stress_model_154( &
       pgaus, pmate, gpgdi, gpigd, gpdet, ielem, flagImpl, gpstr, gpdds)

    ! -----------------------------------------------------------------------------------------
    ! INIT
    ! -----------------------------------------------------------------------------------------
    use def_kintyp, only                                :  &
         ip, rp, lg
    use def_master, only                                :  &
         dtime, itinn, ittim, modul,                       &
         ITER_K_STATE, TIME_N_STATE
    use def_domain, only                                :  &
         xfiel
    use def_solidz, only                                :  &
         stiff0_sld, svegm_sld, parco_sld, celen_sld,      &
         kfl_strai_sld, SLD_GREEN, SLD_INFINITESIMAL,      &
         kfl_insdv_sld, rmate_sld
    use mod_sld_stress_model_comput, only               :  &
         SM_strain_tensor,                                 &
         SM_stress_tensor,                                 &
         SM_stress_transport,                              &
         SM_tensor_to_voigt_second,                        &
         SM_voigt_to_tensor_second,                        &
         SM_tensor_to_voigt_fourth,                        &
         SM_rotate_basis_creation,                         &
         SM_rotate_voigt_second,                           &
         SM_rotate_voigt_fourth,                           &
         SM_PULLBACK

    implicit none

    integer(ip), intent(in)                             :: &
         pgaus,                                            & ! Number of gauss points
         pmate,                                            & ! Current material number
         ielem,                                            & ! Current element number
         flagImpl                                            ! Integration scheme flag

    real(rp),    intent(in)                             :: &
         gpgdi(3,3,pgaus),                                 & ! Displacement gradient
         gpigd(3,3,pgaus),                                 & ! Inverse of deformation gradient tensor
         gpdet(pgaus)                                        ! Determinant of the deformation gradient tensor

    real(rp),    intent(out)                            :: &
         gpstr(3,3,pgaus),                                 & ! 2nd Piola-Kirchoff stresses tensor
         gpdds(3,3,3,3,pgaus)                                ! 2nd elasticity tensor

    logical(lg)                                         :: &
         flagFT, flagFC, flagMT, flagMC, flagVisco,        & ! Flags (damage, damge, damge, damage, viscosity, plasticity, damage)
         flagPlast, flagHygro, flagDamMT

    integer(ip)                                         :: &
         igaus                                               ! Index

    real(rp)                                            :: &
         E11, E22, G12, G23, v12, v21, v23,                & ! Material properties (- elastic)
         xT, xTO, xC, xCO, yT, yC, sL, sT, alpha0,         & !   - Strength
         nL, nT,                                           & !
         gxT, gxTO, gxC, gXCO, gYT, gYC, gSL, gG,          & !   - Fracture toughness parameters
         sP, kP,                                           & !   - Plasticity parameters
         eta, dt,                                          & !   - Viscosity parameters
         a11, a22, incrTemp, b11, b22, incrMois,           & !   - Hygrothermal parameters
         str11, str22, str12,                              & ! Stresses
         gre11, gre22, grePlast, greInela,                 & ! Strains
         phiFT, phiFC, phiMT, phiMC,                       & ! Loading functions
         rFTold, rFCold, rMTold, rMCold,                   & ! Internal variables
         rFTnew, rFCnew, rMTnew, rMCnew,                   &
         d1, d2, d3, d4, d5, d6,                           & ! Damage variables
         lenchar,                                          & ! Characteristic element legth
         alpha(3,3), beta(3,3),                            & ! Hygrothermal tensor in Voigt form.
         gpgre(3,3,pgaus),                                 & ! Strains tensor (mechanical)
         gprot(3,3),                                       & ! Rotation matrix for material orientations
         greGloTot(6),                                     & ! Total strain tensor in Voigt form (global CSYS)
         greMatTot(6),                                     & ! Total strain tensor in Voigt form (material CSYS)
         greMatEla(6),                                     & ! Elastic strain tensor in Voigt form (material CSYS)
         strGloNom(6),                                     & ! Nominal stress tensor in Voigt form (global CSYS)
         strMatEff(6),                                     & ! Stress tensor in Voigt form (material CSYS)
         strMatNom(6),                                     & ! Nominal stress tensor in Voigt form (material CSYS)
         stiMatEff(6,6),                                   & ! Effective tiffness tensor in Voigt form (material CSYS)
         stiMatNom(6,6),                                   & ! Nominal stiffness tensor in Voigt form (material CSYS)
         auxMA33(3,3)
    !
    ! --------------------------------------------------------------------------|  INIT  |-----

    ! -----------------------------------------------------------------------------------------
    ! STRESS TENSOR AND SECOND ELASTICITY TENSOR
    ! -----------------------------------------------------------------------------------------
    ! GET PROPERTIES
    !
    ! Material properties
    !
    call sm154_model_properties( &
         parco_sld(:,pmate),E11, E22, v12, v21, v23, G12, G23, xT, xTO, xC, xCO, yT, yC, sL, &
         alpha0, gxT, gxTO, gxC, gxCO, gyT, gyC, gsL, gG, flagPlast, sP, kP, &
         flagVisco, eta, flagHygro, a11, a22, incrTemp, b11, b22, incrMois)
    !
    ! Characteristics element length (computed at the beginning) 
    !
    if( ittim == 1_ip .and. itinn(modul) == 1_ip ) then
       call sm154_characteristics_length(ielem, lenchar)
       celen_sld(ielem) = lenchar
    end if
    lenchar = celen_sld(ielem)
    !
    ! Strength reduction properties
    !
    call sm154_strength_reduction(&
         E11, E22, G12, gxT, gxTO, gxC, gXCO, gyT, gyC, gsL, lenchar, ielem, &
         xT, xTO, xC, xCO, yT, yC, sL)
    !
    ! Re-computation of some material properties due to strength reduction
    !
    sT = yC*cos(alpha0)*(sin(alpha0) + cos(alpha0)/tan(2.0_rp*alpha0))
    nT = -1.0_rp/tan(2.0_rp*alpha0)
    nL = -sL*cos(2.0_rp*alpha0)/(yC*cos(alpha0)**2)
    !
    ! Undamaged stiffness tensor
    !
    stiMatEff(:,:) = stiff0_sld(:,:,pmate)
    !
    ! Thermal expansion tensor  <<<< NOT IMPLEMENTED >>>>
    !
    alpha(:,:) = 0.0_rp
    alpha(1,1) = a11
    alpha(2,2) = a22
    alpha(3,3) = a22
    !
    ! Moisture tensor  <<<< NOT IMPLEMENTED >>>>
    !
    beta(:,:) = 0.0_rp
    beta(1,1) = b11
    beta(2,2) = b22
    beta(3,3) = b22
    !
    ! ----------------------------------------------------------------
    ! LOOP OVER MATERIAL POINTS
    !
    gpstr(:,:,:) = 0.0_rp
    gpdds(:,:,:,:,:) = 0.0_rp

    ! ....| material_points |................
    material_points: do igaus = 1, pgaus

       ! ------------------------------------------------------------
       ! GET STATE VARIABLES
       !
       rFTold   = max(1.0_rp, svegm_sld(ielem)%a( 1, igaus, TIME_N_STATE))
       rFCold   = max(1.0_rp, svegm_sld(ielem)%a( 2, igaus, TIME_N_STATE))
       rMTold   = max(1.0_rp, svegm_sld(ielem)%a( 3, igaus, TIME_N_STATE))
       rMCold   = max(1.0_rp, svegm_sld(ielem)%a( 4, igaus, TIME_N_STATE))
       grePlast =             svegm_sld(ielem)%a(11, igaus, TIME_N_STATE)
       greInela =             svegm_sld(ielem)%a(12, igaus, TIME_N_STATE)
       !
       if (rMTold > 1.0_rp) then
          flagDamMT = .true.
       else
          flagDamMT = .false.
       end if
       
       ! ------------------------------------------------------------
       ! STRAIN TENSOR
       !
       greMatTot(:) = 0.0_rp
       greGloTot(:) = 0.0_rp
       greMatEla(:) = 0.0_rp
       !
       ! Total strain tensor
       !   0 - Infinitesimal tensor
       !   1 - Green-Lagrange tensor
       !
       if(      kfl_strai_sld == SLD_INFINITESIMAL ) then
          call SM_strain_tensor(0_ip, gpgdi(:,:,igaus), gpgre(:,:,igaus))
       else if( kfl_strai_sld == SLD_GREEN ) then
          call SM_strain_tensor(1_ip, gpgdi(:,:,igaus), gpgre(:,:,igaus))
       end if
       !
       ! From tensor to Voigt notation
       !
       call SM_tensor_to_voigt_second(3_ip, gpgre(:,:,igaus), greGloTot(:))
       !
       ! Rotate from the global to the material coordinate system
       !
       call SM_rotate_basis_creation(ielem, gpgdi(:,:,igaus), gprot(:,:))
       call SM_rotate_voigt_second(1_ip, gprot(:,:), greGloTot(:), greMatTot(:))
       !
       ! Elastic strain tensor
       !
       call sm154_strain_tensor(G12, flagHygro, alpha, incrTemp, beta, incrMois, sP, kP, flagDamMT, &
            greMatTot, greMatEla, grePlast, greInela, flagPlast)

       ! ------------------------------------------------------------
       ! EFFECTIVE STRESS TENSOR AND INVARIANTS
       !
       strMatEff(:) = 0.0_rp
       !
       ! Effective stress tensor
       !
       strMatEff(:) = matmul(stiMatEff(:,:),greMatEla(:))
       !
       ! Invariants
       !
       gre11 = greMatEla(1)
       gre22 = greMatEla(2)
       str11 = strMatEff(1)
       str22 = strMatEff(2)
       str12 = strMatEff(6)

       ! ------------------------------------------------------------
       ! DAMAGE STATE
       !
       ! Failure criteria / loading functions
       !
       call sm154_failure_criteria( &
            gre11, gre22, str11, str22, str12, E11, E22, xT, xC, yT, yC, nL, nT, gG, sL, sT, &
            alpha0, phiFT, phiFC, phiMT, phiMC)
       !
       ! Damage threshold
       !
       dt = dtime
       call sm154_damage_evolution( &
            flagVisco, eta, dt, phiFT, phiFC, phiMT, phiMC, rFTold, rFCold, rMTold, &
            rMCold, rFTnew, rFCnew, rMTnew, rMCnew, flagFT, flagFC, flagMT, flagMC)
       !
       ! Damage state
       !
       call sm154_damage_state( &
            E11, E22, v12, v21, G12, xT, xTO, xC, xCO, yT, yC, sL, sT, nL, nT, gxT, gxTO, gxC, &
            gxCO, gyT, gyC, gsL, gG, alpha0, kP, flagPlast, lenchar, gre11, gre22,             &
            rFTnew, rFCnew, rMTnew, rMCnew, d1, d2, d3, d4, d5, d6)
       !
       ! Initial condition for damage state
       !
       if( kfl_insdv_sld > 0_ip ) then
          d1 = min(d1Max,max(d1,xfiel(kfl_insdv_sld) % a(1,ielem,1)))
          d2 = min(d2Max,max(d2,xfiel(kfl_insdv_sld) % a(2,ielem,1)))
          d3 = min(d3Max,max(d3,xfiel(kfl_insdv_sld) % a(3,ielem,1)))
          d4 = min(d4Max,max(d4,xfiel(kfl_insdv_sld) % a(4,ielem,1)))
          d5 = min(d5Max,max(d5,xfiel(kfl_insdv_sld) % a(5,ielem,1)))
          d6 = min(d6Max,max(d6,xfiel(kfl_insdv_sld) % a(6,ielem,1)))
       end if
       !
       ! ------------------------------------------------------------
       ! NOMINAL STRESS TENSOR
       !
       ! Stiffness tensor
       !
       call sm154_stiffness_tensor(E11, E22, G12, G23, v12, v23, &
               d1, d2, d3, d4, d5, d6, stiMatNom(:,:))
       !
       ! Nominal stress tensor (material csys)
       !
       strMatNom(:) = matmul(stiMatNom(:,:),greMatEla(:))
       !
       ! Nominal stress tensor (global csys)
       !
       call SM_rotate_voigt_second(2_ip, gprot(:,:), strGloNom(:), strMatNom(:))
       !
       ! From Voigt to tensorial notation
       !
       auxMA33(:,:) = 0.0_rp
       call SM_voigt_to_tensor_second(3_ip, auxMA33(:,:), strGloNom(:))
       !
       ! PK2 due to strain measure
       !
       if(     kfl_strai_sld == SLD_INFINITESIMAL) then
          !
          ! Assuming only Infinitesimal
          !
          call SM_stress_transport(SM_PULLBACK, 3_ip, gpdet(igaus), gpgdi(:,:,igaus), gpigd(:,:,igaus), &
               auxMA33(:,:), gpstr(:,:,igaus))

       else if( kfl_strai_sld == SLD_GREEN) then
          !
          ! Assuming full Green Lagrange
          !
          gpstr(:,:,igaus) = auxMA33(:,:)

       end if
       !
       ! ------------------------------------------------------------
       ! ENERGIES
       !
       !
       ! call sm154_energies() ! <AQC> Not implemented yet
       ! ---------------------------------------------------
       ! SECOND ELASTICITY TENSOR
       !
       ! ....| Implicit scheme |............
       if (flagImpl .eq. 1_ip) then

          call runend('MOD_SLD_STRESS_MODEL_154: NOT IMPLEMENTED FOR IMPLICIT!')

       end if
       !...........| Implicit scheme |.....

       ! ---------------------------------------------------
       ! STORING STATE VARIABLES
       svegm_sld(ielem)%a( 1,igaus,ITER_K_STATE) = rFTnew   ! Damage thresholds
       svegm_sld(ielem)%a( 2,igaus,ITER_K_STATE) = rFCnew
       svegm_sld(ielem)%a( 3,igaus,ITER_K_STATE) = rMTnew
       svegm_sld(ielem)%a( 4,igaus,ITER_K_STATE) = rMCnew
       svegm_sld(ielem)%a( 5,igaus,ITER_K_STATE) = d1       ! Damage variable FT/FC
       svegm_sld(ielem)%a( 6,igaus,ITER_K_STATE) = d2       ! Damage varaible MT/MC
       svegm_sld(ielem)%a( 7,igaus,ITER_K_STATE) = d3       !
       svegm_sld(ielem)%a( 8,igaus,ITER_K_STATE) = d4       !
       svegm_sld(ielem)%a( 9,igaus,ITER_K_STATE) = d5       !
       svegm_sld(ielem)%a(10,igaus,ITER_K_STATE) = d6       ! Damage variable shear
       svegm_sld(ielem)%a(11,igaus,ITER_K_STATE) = grePlast
       svegm_sld(ielem)%a(12,igaus,ITER_K_STATE) = greInela

       rmate_sld(ielem)%a(:,:,igaus) = gprot(:,:)           ! Rotating matrix for material orientations
       
    end do material_points

  end subroutine sld_stress_model_154

end module mod_sld_stress_model_154

