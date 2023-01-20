!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!---------------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_interface_element.f90
!> @author  Adria Quintanas (adria.quintanas@udg.edu)
!> @date    January, 2017
!>          - Subroutine written
!> @author  Adria Quintanas and Gerard Guillamet
!> @date    May, 2018
!>          - Includes tangent stiffness calculation methods
!>          - Includes 2-d formulation
!>          - Includes latest formulation Turon et. al 2018
!>          - Includes viscosity regularization
!> @author  Adria Quintanas and Gerard Guillamet
!> @date    July, 2018
!>          - Includes contact law
!>          - Includes analytical tangent Turon et. al 2018
!>          - ...
!>          September, 2021
!> @author  Gerard Guillamet
!>          - Adding dMax as input for the user
!> @brief   Cohesive elements for progressive delamination.
!>
!> @details
!>
!>          Output variables:\n
!>             DCOHE (ISVAR=1)       : Damage variable d
!>             ----- (ISVAR=2)       : Damage threshold r
!>
!>          References:\n
!>
!>          A. Turon, P.P. Camanho, J. Costa, C.G. Davila. A damage model for the
!>          the simulation of delamination in advanced composites under
!>          variable-mode loading. Mechanics of materials, 2006.\n
!>
!>          A. Turon, E.V. Gonzalez, C. Sarrado, G. Guillamet, P. Maimi. Accurate
!>          simulation of delamination under mixed-mode loading using a cohesive
!>          model with a mode-dependent penalty stiffness. Composite
!>          structures, 2018.\n
!>
!>          X. Martinez, S. Oller, F. Rastellini, A.H. Barbat. A numerical
!>          procedure simulation RC structures reinforced with FRP using serial/
!>          parallel mixing theory. Comput Struct, 2008\n
!>          G. Duvaut, J.L. Lions. Les inequations en mecanique et en physique.
!>          Dunod, Paris, 1972\n
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
!>          - Distinguish svegm variable at element and gauss level.
!>          - Posprocess normal vector, mixed mode ratio and ratio dissipated energy.
!>
!> @}
!---------------------------------------------------------------------------------

module mod_sld_interface_element
  !
  ! ==============================================================================
  ! INIT
  !
  use def_kintyp, only                     :  &
      ip, rp, lg
  use def_solidz, only                     :  &
      densi_sld,                              & ! Density
      kfl_fixno_sld,                          & ! Flag nodes with Dirichlet BCs
      lmate_sld,                              & ! List of materials
      lawta_sld,                              & ! Material tangent computation options
      lawch_sld,                              & ! Cohesive law identification
      parch_sld,                              & ! Cohesive law parameters
      parcf_sld,                              & ! Contact model parameters
      nsvar_sld,                              & ! Number of state variables
      svegm_sld                                 ! State variables at gauss point per element
  use def_master, only                     :  &
      zeror
  use def_domain, only                     :  &
      lnods
  ! -----------------------------------------------------------------------------
  implicit none
  ! -----------------------------------------------------------------------------
  real(rp), parameter                      :: &
    tolPe = 1.0e-12_rp                          ! Tolerance for penetration (overclosure)
  integer(ip), parameter                   :: &
    nprop = 9                                   ! Number of material properties
  integer(ip), parameter                   :: &
    ELINT_TANGENT_EXPLICIT   = 0,             & ! VAR tangent computation
    ELINT_TANGENT_ANALYTICAL = 1,             & ! VAR ""
    ELINT_TANGENT_NUMERICAL  = 2,             & ! VAR ""
    ELINT_TANGENT_SECANT     = 3,             & ! VAR ""
    ELINT_LAW_TURON          = 904,           & ! VAR cohesive law identification
    ELINT_LAW_TURON_2018     = 905,           & ! VAR ""
    ELINT_LAW_CONTACT        = 800              ! VAR contact law identification
  integer(ip), parameter                   :: &
    closing    = 0,                           & ! VAR loading status
    opening    = 1,                           & ! VAR ""
    elastic    = 0,                           & ! VAR damage status
    damage     = 1,                           & ! VAR ""
    broken     = 2                              ! VAR ""
  !
  !=============================================================| init |=========
  !==============================================================================
  abstract interface
    ! ----------------------------------
    ! TURON COHESIVE LAW
    subroutine  coheLaw( &
      itang, pdime, dt, properties, jumps, stateOld, stateNew, t, C)
      import                               :: &
        ip, rp, lg
      integer(ip), intent(inout)           :: &
        itang
      integer(ip), intent(in)              :: &
        pdime
      real(rp),    intent(in)              :: &
        dt,                                   &
        properties(:),                        &
        jumps(:)
      real(rp),    intent(inout)           :: &
        stateOld(:),                          &
        stateNew(:)
      real(rp),    intent(out)             :: &
        t(:),                                 &
        C(:,:)
    end subroutine coheLaw
    ! ----------------------------------
    end interface
    procedure(coheLaw),  pointer           :: &
      cohesiveLaw
  !==============================================================================
  !==============================================================================
  ! PUBLIC
  !
  public                                   :: &
    ELINT_TANGENT_EXPLICIT,                   &
    ELINT_TANGENT_ANALYTICAL,                 &
    ELINT_TANGENT_NUMERICAL,                  &
    ELINT_LAW_TURON,                          &
    ELINT_LAW_TURON_2018,                     &
    ELINT_LAW_CONTACT,                        &
    ELINT_elemental_operations,               &
    ELINT_stable_time_increment
    
  !=============================================================| public |=======
  !==============================================================================
  ! PRIVATE
  !
  private                                  :: &
    evaluate_cohesive_law,                    & ! FUN Cohesive law evaluation
    LAW_COH_TURON_ORIGINAL,                   & ! FUN Cohesive law Turon et al. 2006
    LAW_COH_TURON_CURRENT,                    & ! FUN Cohesive law Turon et al. 2018
    LAW_COH_TURON_CURRENT_Cuu,                & ! FUN Cuu analytical Turon et al. 2018
    CONTACT_LAW_UNILATERAL,                   & ! FUN Contact law
    midcoord_and_djump,                       &
    shape_functions,                          &
    rotation_matrix,                          &
    b_matrix_half_positive,                   &
    displacement_jump_local_csys,             &
    internal_force_vector,                    &
    tangent_stiffness_matrix,                 &
    nodes_dirichlet_correction,               &
    nodes_dirichlet_reaction_forces
  !
  !=============================================================| private |======
  !==============================================================================
  ! CONTAINS
  !
  contains
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !------------------------------------------------------------------------------
  !> @author  Adria Quintanas (adria.quintanas@udg.edu)
  !> @date    13/01/2017
  !> @brief   Elemental operations
  !> @details
  !------------------------------------------------------------------------------
  subroutine ELINT_elemental_operations( &
    itask, dt, ielem, pdime, nsvar, pnode, pevat, elcoor, eldisp, elrhs, elmat, &
    elfrx )
    use def_master,                    only : ITER_K_STATE, TIME_N_STATE
    use def_solidz,                    only : SLD_IMPLICIT_SCHEME
    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    integer(ip),    intent(in)             :: &
      itask,                                  & !< Itask
      ielem,                                  & !< Element number
      pdime,                                  & !< Number of dimensions
      nsvar,                                  & !< Number of state variables
      pnode,                                  & !< Number of nodes
      pevat                                     !< Number of degree of freedom {ndime*node}
    real(rp),       intent(in)             :: &
      dt,                                     & !< Time increment
      elcoor(pdime,pnode),                    & !< Nodal element coordinates {dime,node}
      eldisp(pdime,pnode)                       !< Nodal element displacement {dime,node}
    real(rp),       intent(out)            :: &
      elrhs(pevat),                           & !< Elemental RHS vector {dof}
      elmat(pevat,pevat),                     & !< Elemental tangent matrix {dof,dof}
      elfrx(pdime,pnode)                        !< Elemental reaction forces vector {dime,node}
    ! ---------------------------------------------------------------------------
    integer(ip)                            :: &
      imate,                                  & ! Material number
      ilaw,                                   & ! Cohesive law number
      itang,                                  & ! Cohesive law material tangent computation
      sgaus,                                  & ! Mid-surface number of gauss points
      snode,                                  & ! Mid-surface number of nodes
      sdime,                                  & ! Mid-surface number of dimensions
      pevah, peva1,                           & ! Auxiliar variable, half of the total dofs
      igaus                                     ! Loop index
    real(rp)                               :: &
      elfint(pevat/2),                        & ! Positive part of the elemental internal force vector (node*dof/2)
      elstif(pevat/2,pevat/2),                & ! Positive part of the elemental tangent matrix  (node*dof/2, node*dof/2)
      sfcoor(pdime,pnode/2),                  & ! Midsurface coordinates (dime, mid_node)
      sfjump(pevat/2),                        & ! Midsurface relative displacement (mid_node*dof/2)
      gprota(pdime,pdime),                    & ! Integration point midsurface rotation matrix R (dime, dime)
      gpshap(pnode/2),                        & ! Integration point midsurface shape function N (mid_gauss)
      gpderi(pdime-1,pnode/2),                & ! Integration point midsurface deriv. of the shape function dN/dn (mid_dime,mid_gau)
      gpbmat(pdime,pevat/2),                  & ! Integration point midsurface relative displacement matrix (dime, mid_gauss*dof/2)
      gpjump(pdime),                          & ! Integration point midsurface displacement jump (dime)
      gptrac(pdime),                          & ! Integration point midsurface tractions / tau (local coordinate system) (dime)
      gpdtdj(pdime,pdime),                    & ! Integration point midsurface derivative of tractions wrt displ. jump (local csys) 
      gpweig,                                 & ! Integration point weight
      gparea,                                 & ! Integration point mapped area
      gpSVold(nsvar), gpSVnew(nsvar),         & ! Integration point state variable
      properties(nprop)
    !
    ! ---------------------------------------------------------------------------
    !
    ! Read material properties and options in function of the material number
    !
    imate = lmate_sld(ielem)
    ilaw  = lawch_sld(imate)
    !
    if( itask == SLD_IMPLICIT_SCHEME ) then
      itang = lawta_sld(imate)
    else
      itang = ELINT_TANGENT_EXPLICIT
    endif
    !
    if ( ilaw == ELINT_LAW_CONTACT ) then
       ! Contact law parameters
       properties(:) = parcf_sld(1:nprop,imate)
    else
       ! Cohesive law parameters
       properties(:) = parch_sld(1:nprop,imate)
    end if
    !
    ! Auxiliar dimensions
    !
    sgaus = int(pnode/2,ip) ! Mid-surface gauss points (It has to be equal to half of the number of nodes)
    snode = int(pnode/2,ip) ! Mid-surface number of nodes
    sdime = int(pdime-1,ip) ! Mid-surface number of dimensions
    pevah = int(pevat/2,ip) ! Mid-surface number of degree of freedom
    !
    ! Midsurface coordinates and displacement jump
    !
    call midcoord_and_djump(pdime, pnode, snode, pevah, eldisp, elcoor, sfcoor, sfjump)
    !
    ! Loop over the interpolating points
    !
    elfint(:)   = 0.0_rp
    elstif(:,:) = 0.0_rp
    elfrx(:,:)  = 0.0_rp
    !
    gausspoints: do igaus = 1, sgaus
      !
      ! Get gp state old variable
      !
      gpSVold(:) = svegm_sld(ielem)%a(:,igaus,TIME_N_STATE)
      !
      ! Interpolating function and its derivative (N and dN/dn)
      !
      call shape_functions( igaus, pdime, sdime, snode, gpweig, gpshap(:), gpderi(:,:) )
      !
      ! Rotation matrix and area diferential(R and dA)
      !
      call rotation_matrix( pdime, sdime, snode, sfcoor(:,:), gpderi(:,:), gprota(:,:), gparea )
      !
      ! B matrix (B)
      !
      call b_matrix_half_positive( pdime, snode, pevah, gpshap(:), gprota(:,:), gpbmat(:,:) )
      !
      ! Relative displacements / displacement jump (delta)
      !
      call displacement_jump_local_csys( pdime, pevah, gpbmat(:,:), sfjump(:), gpjump(:) )
      !
      ! Constitutive relationship tensors (Traction and dTraction/dJump)
      !
      call evaluate_cohesive_law( ilaw, itang, pdime, nsvar, dt, properties(:), gpjump(:), &
           gpSVold(:), gpSVnew(:), gptrac(:), gpdtdj(:,:) )
      !
      ! Residual vector (RHS) (add the current integration point contribution)
      !
      call internal_force_vector( pdime, pevah, gpbmat(:,:), gptrac(:), gpweig, gparea, elfint(:) )
      !
      ! Tangent stiffness matrix (K) (add the current integration point contribution)
      !
      if( itask == SLD_IMPLICIT_SCHEME ) then
         call tangent_stiffness_matrix( pdime, pevah, gpbmat(:,:), gpdtdj(:,:), gpweig, gparea, elstif(:,:))
      end if
      !
      ! Set gp state new variable
      !
      svegm_sld(ielem)%a(:,igaus,ITER_K_STATE) = gpSVnew(:)
      !
    end do gausspoints
    !
    ! Complete vectors and matrix
    !
    elrhs(:)   = 0.0_rp
    elmat(:,:) = 0.0_rp
    peva1 = pevah + 1
    if( pdime == 2 ) then
       !  - residual vector (RHS = Fext - Fint)
       elrhs(1:2) = elfint(1:2)
       elrhs(3:4) = elfint(3:4)
       elrhs(5:6) = -elrhs(3:4)
       elrhs(7:8) = -elrhs(1:2)
       if( itask == SLD_IMPLICIT_SCHEME ) then
          !  - tangent stiffnes matrix
          elmat(1:2,7:8) = -elstif(1:2,1:2) ! I quadrant (-)
          elmat(3:4,5:6) = -elstif(3:4,3:4) ! I quadrant (-)
          !
          elmat(1:2,1:2) = elstif(1:2,1:2)  ! II quadrant (+)
          elmat(3:4,3:4) = elstif(3:4,3:4)  ! II quadrant (+)
          ! Symetric
          elmat(5:6,3:4) = elmat(3:4,5:6)
          elmat(7:8,1:2) = elmat(1:2,7:8)
          elmat(5:6,5:6) = elmat(3:4,3:4)
          elmat(7:8,7:8) = elmat(1:2,1:2)
       end if
    else
       !  - residual vector (RHS = Fext - Fint)
       elrhs(1:pevah)     =  elfint(1:pevah)
       elrhs(peva1:pevat) = -elfint(1:pevah)
       if( itask == SLD_IMPLICIT_SCHEME ) then
          !  - tangent stiffnes matrix
          elmat(1:pevah,1:pevah)         =  elstif(1:pevah,1:pevah)
          elmat(peva1:pevat,peva1:pevat) =  elstif(1:pevah,1:pevah)
          elmat(1:pevah,peva1:pevat)     = -elstif(1:pevah,1:pevah)
          elmat(peva1:pevat,1:pevah)     = -elstif(1:pevah,1:pevah)
       end if
    end if
    !
    ! Reaction forces at Dirichlet
    !
    call nodes_dirichlet_reaction_forces( ielem, pdime, pnode, pevat, elrhs(:), elfrx(:,:) )
    !
    ! Dirichlet BC correction
    !
    if ( itask == SLD_IMPLICIT_SCHEME ) call nodes_dirichlet_correction( ielem, pdime, pnode, pevat, elmat(:,:), elrhs(:) )

  end subroutine ELINT_elemental_operations

  ! #####################################################################
  !
  ! ---------------------------------------------------------------------
  !> @author  Adria Quintanas
  !> @brief   Evaluate cohesive law
  !> @details Evaluate cohesive law
  ! ---------------------------------------------------------------------
  subroutine evaluate_cohesive_law( &
    ilaw, itang, pdime, nsvar, dt, properties, jumps, stateOld, stateNew, &
    t, C )
    ! ----------------------------------
    implicit none
    ! ----------------------------------
    integer(ip),  intent(in)          :: &
        ilaw,                            & ! ID-law
        itang,                           & ! ID-tangent stiffness computation
        pdime,                           & ! Number of dimensions
        nsvar                              ! Number of state variables
    real(rp), intent(in)              :: &
        dt,                              & ! Time increment
        properties(:),                   & ! Interface properties vector
        jumps(:)                           ! Displacement jump vector
    real(rp),  intent(inout)          :: &
        stateOld(:),                     & ! State variables OLD
        stateNew(:),                     & ! State variables NEW
        t(:),                            & ! Traction vector
        C(:,:)                             ! Tangent stiffnes
    ! ----------------------------------
    integer(ip)                       :: &
        idime, jtang
    real(rp)                          :: &
        dJ,                              & !
        Jnum(pdime),                     & !
        tNum(pdime),                     &
        CNum(pdime,pdime),               &
        Caux(pdime,pdime),               &
        stateNum(nsvar)
    ! ----------------------------------
    !
    t(:) = 0.0_rp
    C(:,:) = 0.0_rp
    !
    nullify(cohesiveLaw)
    !
    ! Selecting law
    !
    select case( ilaw )

        case( ELINT_LAW_TURON )

          cohesiveLaw => LAW_COH_TURON_ORIGINAL

        case( ELINT_LAW_TURON_2018 )

          cohesiveLaw => LAW_COH_TURON_CURRENT

        case( ELINT_LAW_CONTACT )

          cohesiveLaw => CONTACT_LAW_UNILATERAL

        case default

          call runend('MOD_SLD_INTERFACE_ELEMENT: COHESIVE LAW NOT ASSOCIATED')

    end select
    !
    ! Set a tangent type
    !
    jtang = itang
    !
    ! Compute law
    !
    call cohesiveLaw( jtang, pdime, dt, properties(:), jumps(:), stateOld(:), stateNew(:),  &
      t(:), C(:,:) )
    !
    if( jtang == ELINT_TANGENT_NUMERICAL ) then
      !
      ! - Cuu
      CNum(:,:) = 0.0_rp
      !
      do idime = 1, pdime
        ! dJ
        Jnum(:) = jumps(:)
        dJ = get_variation(Jnum(idime)) ! Get ficticious perturbation
        Jnum(idime) = Jnum(idime) + dJ  ! Adds ficticious perturbation
        ! dTu
        call cohesiveLaw( jtang, pdime, dt, properties(:), Jnum(:), stateOld(:), stateNum(:),  &
          tNum(:), Caux(:,:) )
        ! Contribution
        CNum(:,idime) = (tNum(:) - t(:))/dJ
      enddo
      C(:,:) = Cnum(:,:)
    endif
    !
  contains
    !
    ! ---------------------------------------------------------------------
    !> @author  Adria Quintanas
    !> @date    May, 2018
    !> @brief   Ficticious perturbation in the displacement jump
    !> @details Get variation of the displacement jump. It is important
    !>          to take into account the sign.
    ! ---------------------------------------------------------------------
    function get_variation(x) result(dx)
      real(rp), intent(in)      :: x
      real(rp)                  :: dx
      real(rp), parameter       :: p = 1.0E-4_rp
      if( x > 1.0E-19_rp ) then
         dx =  1.0e-8_rp ! max(x*p, 1.0E-14_rp)
      else
         dx = -1.0e-8_rp ! min(x*p,-1.0E-14_rp)
      endif
    end function get_variation

  end subroutine evaluate_cohesive_law

  !------------------------------------------------------------------------------
  !> @author  Adria Quintanas (adria.quintanas@udg.edu)
  !> @date    10/01/2017
  !> @brief   Computation of Tractions and dTractions/dJump
  !> @details
  !------------------------------------------------------------------------------
  subroutine LAW_COH_TURON_ORIGINAL( &
       itang, pdime, dt, properties, gpjump, stateOld, stateNew, gptrac, gpdtdj )
    ! ---------------------------------------
    implicit none
    ! ---------------------------------------
    integer(ip),    intent(in)             :: &
         pdime                                  !< Number of dimensions of the problem
    integer(ip),    intent(inout)          :: &
         itang                                  !< ID element
    real(rp),       intent(in)             :: &
         dt,                                  & !< Time increment
         properties(:),                       & !< Cohesive law properties
         gpjump(:)                              !< Displacement jump (local coordinate system) (dime)
    real(rp),       intent(inout)          :: &
         stateOld(:),                         & !< State variables : old (for possible initialize prupouses)
         stateNew(:)                            !< State variables : new
    real(rp),       intent(out)            :: &
         gptrac(:),                           & !< Tractions / tau (local coordinate system) (dime)
         gpdtdj(:,:)                            !< Tangent stiffness tensor (dime, dime)
    ! ---------------------------------------
    logical(lg)                            :: &
         fload, fpene,                        & ! Loading and penetration flag
         flagVisco                              ! Viscosity flag
    integer(ip)                            :: &
         idime, jdime
    real(rp)                               :: &
         tauI, tauII, kpen, gIc, gIIc, eta,   & ! Interface properties
         rho,                                 & ! Viscosity
         dMax,                                & ! Max. damage variable
         beta, grat,                          & ! Onset criterion parameters
         openI, openII,                       & ! Critical open displacement
         dispI, dispII,                       & ! Mode I and II displacement jump
         dini, dfin,                          & ! Initial and final displacement
         jeqv,                                & ! Equivalent displacement jump
         r, d, dhist,                         & ! Internal damage and state damage
         gpstif(pdime,pdime),                 & ! Stiffness tensor (dime, dime)
         auxs1
    ! ---------------------------------------
    !
    ! Reading material properties
    !
    gIc   = properties(1)
    gIIc  = properties(2)
    tauI  = properties(3)
    tauII = tauI*sqrt(gIIc/gIc)
    kpen  = properties(4)
    eta   = properties(5)
    rho   = properties(6)
    if( rho > 0.0_rp ) then
       flagVisco = .true.
    else
       flagVisco = .false.
    end if
    dMax  = properties(8)
    ! Defaults
    if( abs(dMax) < zeror ) then
      dMax = 1.0_rp
    end if
    !
    ! Initial variables
    !
    fload = .false.
    fpene = .false.
    !
    ! Mode I and II displacement jeqv (dI and dII)
    !
    dispI  = gpjump(pdime)
    !
    dispII = 0.0_rp
    do idime = 1, pdime - 1
       dispII = dispII + gpjump(idime)**2
    end do
    dispII = sqrt(dispII)
    !
    ! Mixed mode ratios
    !
    if (dispI <= tolPe) then
       beta = 1.0_rp
       jeqv = dispII
    else
      beta = dispII/(dispI + dispII)
      jeqv = sqrt(dispI**2 + dispII**2)
    end if
    !
    ! Mixed mode onset and final displacement
    !
    openI  = tauI/kpen
    openII = tauII/kpen
    if (eta >= 0.0_rp) then
      !   - B-K criterion
      grat = (beta**2)/(1.0_rp + 2.0_rp*(beta**2) - 2.0_rp*beta)
      dini = sqrt((openI**2) + (openII**2 - openI**2)*(grat**eta))
      dfin = 2.0_rp*(gIc + (gIIc - gIc)*(grat**eta))/(kpen*dini)
    else
      !  - linear criterion
      grat = (beta**2)/(1.0_rp + 2.0_rp*(beta**2) - 2.0_rp*beta)
      dini = openI*openII/sqrt((1.0_rp - grat)*openI**2 + grat*(openII**2))
      dfin = 2.0_rp*((gIc*gIIc)/((1.0_rp - grat)*gIIc + grat*gIc))/(kpen*dini)
    end if
    !
    ! Mixed mode damage threeshold
    !
    if( properties(7) > zeror ) then
      dhist = properties(7)
    else
      dhist = stateOld(1)
    endif
    !
    r = dini*dfin/(dfin - dhist*(dfin - dini))
    !
    ! Viscous regularization (Duvaut and Lions)
    !
    if (flagVisco) then
       r = r*dt/(rho + dt) + stateOld(2)*rho/(rho + dt)
    end if
    !
    ! Internal variables & loading state flags
    !
    if (jeqv > r) then
       r = jeqv
       fload = .true.   ! loading keyflag
    else
       fload = .false.  ! unloading keyflag
    end if
    !
    if (dispI < tolPe) then
       fpene = .true.   ! penetration keyflag
    else
       fpene = .false.
    end if
    !
    ! Damage (d)
    !
    if( r > 0.0_rp ) then
       d = dfin*(r - dini)/(r*(dfin - dini))
    else
       d = 0.0_rp
    end if
    !  - Is it fully damaged (fully open)?
    if( d >= dMax ) then
       d = dMax
       fload = .false.
    end if
    !
    ! Stiffness tensor (D)
    !
    gpstif(:,:) = 0.0_rp
    do idime = 1, pdime
       gpstif(idime,idime) = (1.0_rp - d)*kpen
    end do
    !
    if (fpene) then
       gpstif(pdime, pdime) = kpen
    end if
    !
    ! Tractions
    !
    gptrac(:) = 0.0_rp
    do idime = 1, pdime
       do jdime = 1, pdime
          gptrac(idime) = gptrac(idime) + gpstif(idime,jdime)*gpjump(jdime)
       end do
    end do
    !
    ! Tangent stiffness (dTdj)
    !
    if(     itang == ELINT_TANGENT_NUMERICAL ) then
       if( d > 0.995_rp ) then
          itang = ELINT_TANGENT_SECANT
       endif
    endif
    !
    if(     itang == ELINT_TANGENT_SECANT ) then
       !
       gpdtdj(:,:) = gpstif(:,:)

    elseif( itang == ELINT_TANGENT_ANALYTICAL) then
       !
       gpdtdj(:,:) = gpstif(:,:)
       !   - Loading
       if (fload) then
          if( r > 0.0_rp ) then
             auxs1 =  (dfin*dini*kpen)/(r**3*(dfin - dini))
          else
             auxs1 = 0.0_rp
          end if
          do idime = 1, pdime
             do jdime = 1, pdime
                gpdtdj(idime,jdime) = gpdtdj(idime,jdime) - auxs1*gpjump(idime)*gpjump(jdime)
             end do
          end do
          !   - Penetration
          if (fpene) then
             do idime = 1, pdime
                gpdtdj(pdime,idime) = gpstif(pdime,idime)
                gpdtdj(idime,pdime) = gpstif(idime,pdime)
             end do
          end if
       endif

    else! ELINT_TANGENT_EXPLICIT or ELINT_TANGENT_NUMERICA
       gpdtdj(:,:) = 0.0_rp

    end if
    !
    ! Update state variables
    !
    stateNew(1) = d
    !
    ! ---------------------------------------
  end subroutine LAW_COH_TURON_ORIGINAL

  !------------------------------------------------------------------------------
  !> @author  Adria Quintanas & Gerard GUillamet
  !> @date    May, 2018
  !> @brief   Cohesive law Turon 2018
  !> @details
  !------------------------------------------------------------------------------
  subroutine LAW_COH_TURON_CURRENT( &
       itang, pdime, dt, properties, jumps, stateOld, stateNew, tu, Cuu )
    ! ----------------------------------
    implicit none
    ! ----------------------------------
    integer(ip),        intent(inout)  :: &
         itang                              !< Tangent calculation method
    integer(ip),        intent(in)     :: &
         pdime                              !< Number of dimensions of the problem
    real(rp),           intent(in)     :: &
         dt,                              & !< Time increment
         properties(:),                   & !< Cohesive law properties
         jumps(:)                           !< Displacement jump (local coordinate system) (dime)
    real(rp),           intent(inout)  :: &
         stateOld(:),                     & !< State variables : old (for possible initialize purpouses)
         stateNew(:)                        !< State variables : new
    real(rp),           intent(out)    :: &
         tu(:),                           & !< Tractions / tau (local coordinate system) (dime)
         Cuu(:,:)                           !< Tangent stiffness tensor (dime, dime)
    ! ----------------------------------
    logical(lg)                        :: &
         flagVisco                          ! Viscosity flag
    integer(ip)                        :: &
         idime,                           &
         kfLoad, kfDamage                   ! Keyflags for damage and loading state
    real(rp)                           :: &
         tauI, tauII,                     & ! Properties: strength
         GI, GII,                         & ! Properties: toughness
         eta,                             & ! Properties: mode-mix
         rho,                             & ! Viscosity
         dMax,                            & ! Max. damage variable
         K11, K22, K33, KCONT,            & ! Properties: penalty stiffness
         JI0, JII0, JIC, JIIC,            & ! Properties: pure mode softening (0) and propagation (C) jumps
         J1, J2, J3, J3p, J3n,            & ! Current jumps + effective opening jump
         Ep,                              & ! Positive elastic energy
         B, KB, KS,                       & ! Mixed mode ratio & Mode depenedent penalty stiffness
         h, h0, hc,                       & ! Equivalent law jumps
         r, d,                            & ! Damage variables
         denom, numer                       ! Auxiliary varaible
    ! ----------------------------------
    !
    ! Init flags
    kfLoad   = opening
    kfDamage = elastic
    !
    ! Material properties
    GI    = properties(1)
    GII   = properties(2)
    tauI  = properties(3)
    tauII = properties(4)
    eta   = properties(5)
    K33   = properties(6)
    KCONT = properties(6) ! Assumed equal to the user-defined
    rho   = properties(7)
    if ( rho > 0.0_rp) then
       flagVisco = .true.
    else
       flagVisco = .false.
    end if
    if (properties(8) > 0.0_rp) then
       stateOld(1) = 1.0_rp
       stateOld(2) = 1.0_rp
    end if
    dMax  = properties(9)
    ! Defaults
    if( abs(dMax) < zeror ) then
      dMax = 1.0_rp
    end if
    !
    ! Effective opening jumps
    J1 = jumps(1)
    if (pdime == 2_ip) then  ! 2-d
       J2  = 0.0_rp
    else                     ! 3-d
       J2  = jumps(2)
    end if
    J3 = jumps(pdime)
    ! Numerical correction (values close to 0)
    ! Positive (J3p) and negative (J3n) normal jumps
    if( J3 > tolPe ) then
       J3p = J3
       J3n = 0.0_rp
       kfLoad = opening
    else
       J3p = 0.0_rp
       J3n =  min(J3, 0.0_rp)
       kfLoad = closing
    endif
    !
    ! Penalty stiffness (mode-II / shear)
    K11 = K33*(gI/gII)*(tauII/tauI)**2
    K22 = K11            ! 2-d (K22 not exists); 3-d (K22 = K11 assumption)
    !
    ! Pure mode openings
    JI0  = tauI/K33
    JII0 = tauII/K11
    JIC  = 2.0_rp*GI/tauI
    JIIC = 2.0_rp*GII/tauII
    !
    ! Mix-mode ratio
    denom = (K11*(J1**2) + K22*(J2**2) + K33*(J3p**2))
    numer = (K11*(J1**2) + K22*(J2**2))
    if( denom /= 0.0_rp ) then
       B = numer/denom
    else
       B = 0.0_rp
    endif
    !
    ! Mode dependent penalty stiffness
    KS = K11
    KB = (1.0_rp - B)*K33 + B*KS
    !
    ! Equivalent jumps
    ! - softening
    h0 = 0.0_rp
    if( KB > 0.0_rp ) then
       h0 = sqrt((K33*(JI0**2) + (KS*(JII0**2) - K33*(JI0**2))*(B**eta))/KB)
    else
       call runend('MOD_SLD_INTERFACE_ELEMENT: COHESIVE TURON HYP A: KB < 0.0')
    end if
    ! - propagation
    hc = 0.0_rp
    if( KB > 0.0_rp .and. h0 > 0.0_rp ) then
       hc = (K33*(JI0*JIC) + (KS*(JII0*JIIC) - K33*(JI0*JIC))*(B**eta))/(KB*h0)
    else
       call runend('MOD_SLD_INTERFACE_ELEMENT: COHESIVE TURON HYP A: h0 < 0.0')
    endif
    ! - checking
    if( hc - h0 < 0.0_rp )then
       call runend('MOD_SLD_INTERFACE_ELEMENT: COHESIVE TURON HYP A: hc < h0')
    endif
    ! - current jump
    denom = sqrt((K11**2)*(J1**2) + (K22**2)*(J2**2) + (K33**2)*(J3p**2))
    numer = (K11*(J1**2) + K22*(J2**2) + K33*(J3p**2))
    if (denom /= 0.0_rp) then
       h = numer/denom
    else
       h = 0.0_rp
    endif
    !
    ! Damage state
    ! - current threshold
    if( h < h0 )then
        r = 0.0_rp
        kfDamage = elastic
    elseif( h0 <= h .and. h < hc  )then
        ! - damage threshold
        r = (h - h0)/(hc - h0)
        ! - Viscous regularization (Duvaut and Lions)
        if (flagVisco) then
           r = r*dt/(rho + dt) + stateOld(2)*rho/(rho + dt)
        end if
        kfDamage = damage
    else
        r = 1.0_rp
        kfDamage = broken
    endif
    ! - historic threshold
    if( r < stateOld(2) )then
       r = stateOld(2)
       kfDamage = elastic
    endif
    ! - state
    denom = r*hc + (1.0_rp - r)*h0
    if( denom > 0.0_rp ) then
       d = min((r*hc)/denom,dMax  )
    else
       d = 0.0_rp
    end if
    !
    ! Tractions
    do idime=1,pdime - 1
       tu(idime) = (1.0_rp - d)*K11*jumps(idime)  ! K11 = K22 assumption
    end do
    if( kfLoad == closing )then
       tu(pdime) = KCONT*J3n
    else
       tu(pdime) =  (1.0_rp - d)*K33*J3p
    end if
    !
    ! Elastic Energy
    Ep = 0.5_rp*(K11*(J1**2) + K22*(J2**2) + K33*(J3p**2))
    !
    ! Material tangent
    if(     itang == ELINT_TANGENT_NUMERICAL ) then
       if( r > 0.995_rp ) then
          itang = ELINT_TANGENT_SECANT
       endif
    endif
    !
    if(     itang == ELINT_TANGENT_SECANT ) then
       !
       Cuu(:,:) = 0.0_rp
       do idime=1,pdime-1
          Cuu(idime,idime) = (1.0_rp - d)*K11  ! K11 = K22 assumption
       end do
       if( kfLoad == opening ) then
          ! - Loading
          Cuu(pdime,pdime) = (1.0_rp - d)*K33
       else
          ! - Penetration
          Cuu(pdime,pdime) = KCONT
       endif

    elseif( itang == ELINT_TANGENT_ANALYTICAL) then

       call LAW_COH_TURON_CURRENT_Cuu( pdime, J1, J2, J3p, K11, K22, K33, KCONT, B, JI0, JIC, JII0, JIIC, &
            eta, KS, KB, d, hc, h0, h, r, kfLoad, kfDamage, Cuu(:,:) )

    else
      ! ELINT_TANGENT_EXPLICIT or ELINT_TANGENT_NUMERICAL
      Cuu(:,:) = 0.0_rp

    end if
    !
    ! Update state variables
    !
    stateNew(1) = d
    stateNew(2) = r
    !
    ! ---------------------------------------
  end subroutine LAW_COH_TURON_CURRENT

  ! ---------------------------------------------------------------------
  !> @author  Adria Quintanas and Gerard Guillamet
  !> @date
  !> @brief   Tangent sitffness matrix Turon et al. 2018
  !> @details Compatibility with 2d and 3d
  !> @todo
  ! ---------------------------------------------------------------------
  subroutine LAW_COH_TURON_CURRENT_Cuu( &
       pdime, J1, J2, J3, K11, K22, K33, KCONT, B, J30, J3C, JS0, JSC, eta, KS, KB, &
       d, hc, h0, h, r, kfLoad, kfDamage, Cuu )
    ! ----------------------------------
    implicit none
    ! ----------------------------------
    integer(ip), intent(in)            :: &
         pdime, kfLoad, kfDamage
    real(rp),    intent(in)            :: &
         J1, J2, J3,                      &
         K11, K22, K33, KCONT, KB, KS,    &
         B, eta,                          &
         JS0, JSC, J30, J3C,              &
         h0, hc, h,                       &
         d, r
    ! ----------------------------------
    real(rp),    intent(inout)         :: &
         Cuu(:,:)
    ! ----------------------------------
    integer(ip)                        :: &
         idime
    ! ----------------------------------
    real(rp)                           :: &
         drdh,  drdh0, drdhc,             &
         dhdJ1, dhdJ2, dhdJ3,             &
         dBdJ1, dBdJ2, dBdJ3,             &
         dKBdB,                           &
         dh0dB, dh0dKB,                   &
         dh0dJ1, dh0dJ2, dh0dJ3,          &
         dhcdB, dhcdKB, dhcdh0,           &
         dhcdJ1, dhcdJ2, dhcdJ3,          &
         dddh0, dddhc, dddr,              &
         dddJ1, dddJ2, dddJ3,             &
         denom
    ! ----------------------------------
    real(rp), parameter                :: &
         zeroerr = 0.0_rp,                &
         oneerr  = 1.0_rp,                &
         Kres    = 0.0_rp
    ! ----------------------------------
    !
    if(     kfDamage == elastic )then
       !
       ! Elastic stage
       !
       Cuu(:,:) = 0.0_rp
       do idime=1,pdime-1
          Cuu(idime,idime) = (1.0_rp - d)*K11  ! K11 = K22 assumption
       end do
       if(     kfLoad == opening )then
          ! - Loading
          Cuu(pdime,pdime) = (1.0_rp - d)*K33
       elseif( kfLoad == closing )then
          ! - Unloading
          Cuu(pdime,pdime) = KCONT
       endif

    elseif( kfDamage == damage  )then
       !
       ! Damage stage
       !
       ! dr/dh, dr/dh0, dh/dhc
       drdh  = 1.0_rp/(hc - h0)
       drdh0 = (h - hc)/(hc - h0)**2
       drdhc = (h0 - h)/(hc - h0)**2
       ! dd/dr, dd/dh0, dd/dhc
       if( d > zeroerr .and. d < oneerr ) then
          dddr  =  h0*hc/(hc*r - h0*(r - 1.0_rp))**2
          dddh0 =  hc*r*(r - 1.0_rp)/(r*(hc - h0) + h0)**2
          dddhc = -h0*r*(r - 1.0_rp)/(r*(hc - h0) + h0)**2
       else
          dddr  = 0.0_rp
          dddh0 = 0.0_rp
          dddhc = 0.0_rp
       endif
       ! dh/dJ
       denom = sqrt( ((K11**2)*(J1**2) + (K22**2)*(J2**2) + (K33**2)*(J3**2))**3 )
       if( denom > zeroerr ) then
          if( J3 > 0.0_rp ) then
             dhdJ1 = K11*J1*((K11**2)*(J1**2) - K11*(K22*(J2**2) + K33*J3**2) + 2.0_rp*((K22**2)*(J2**2) + (K33**2)*(J3**2)))/denom
             dhdJ2 = K11*J2*((K22**2)*(J2**2) - K22*(K11*(J1**2) + K33*J3**2) + 2.0_rp*((K11**2)*(J1**2) + (K33**2)*(J3**2)))/denom
             dhdJ3 = K33*J3*((K33**2)*(J3**2) - K33*(K11*(J1**2) + K22*J2**2) + 2.0_rp*((K11**2)*(J1**2) + (K22**2)*(J2**2)))/denom
          else
             dhdJ1 = sign(1.0_rp, J1)
             dhdJ2 = sign(1.0_rp, J2)
             dhdJ3 = 0.0_rp
          endif
       else
          dhdJ1 = 0.0_rp
          dhdJ2 = 0.0_rp
          dhdJ3 = 0.0_rp
       endif
       ! dB/dJ
       denom = (K11*(J1**2) + K22*(J2**2) +  K33*(J3**2))**2
       if( denom > zeroerr ) then
          if( J3 > 0.0_rp ) then
             dBdJ1 =  2.0_rp*K11*K33*J1*(J3**2)/denom
             dBdJ2 =  2.0_rp*K22*K33*J2*(J3**2)/denom
             dBdJ3 = -2.0_rp*K33*J3*(K11*(J1**2) + K22*(J2**2))/denom
          else
             dBdJ1 = 0.0_rp
             dBdJ2 = 0.0_rp
             dBdJ3 = 0.0_rp
          endif
       else
          dBdJ1 = 0.0_rp
          dBdJ2 = 0.0_rp
          dBdJ3 = 0.0_rp
       endif
       ! dKB/dB
       dKBdB  = KS - K33
       ! dh0/dB, dh0/dKB, dh0/dJ
       dh0dB  = ((KS*(JS0**2) - K33*(J30**2))*eta*(B**eta))/(2.0_rp*h0*B*KB)
       dh0dKB = -h0/(2.0_rp*KB)
       dh0dJ1 = (dh0dB + dh0dKB*dKBdB)*dBdJ1
       dh0dJ2 = (dh0dB + dh0dKB*dKBdB)*dBdJ2
       dh0dJ3 = (dh0dB + dh0dKB*dKBdB)*dBdJ3
       ! dhc/dB, dhc/dKB, dhc/dh0, dhc/dJ
       dhcdB  = ((KS*JS0*JSC - K33*J30*J3C)*eta*(B**eta))/(h0*B*KB)
       dhcdKB = -hc/KB
       dhcdh0 = -hc/h0
       dhcdJ1 = (dhcdB + dhcdKB*dKBdB)*dBdJ1 + dhcdh0*dh0dJ1
       dhcdJ2 = (dhcdB + dhcdKB*dKBdB)*dBdJ2 + dhcdh0*dh0dJ2
       dhcdJ3 = (dhcdB + dhcdKB*dKBdB)*dBdJ3 + dhcdh0*dh0dJ3
       ! dd/dJ
       dddJ1 = dddr*(drdh*dhdJ1 + drdh0*dh0dJ1 + drdhc*dhcdJ1) + dddh0*dh0dJ1 + dddhc*dhcdJ1
       dddJ2 = dddr*(drdh*dhdJ2 + drdh0*dh0dJ2 + drdhc*dhcdJ2) + dddh0*dh0dJ2 + dddhc*dhcdJ2
       dddJ3 = dddr*(drdh*dhdJ3 + drdh0*dh0dJ3 + drdhc*dhcdJ3) + dddh0*dh0dJ3 + dddhc*dhcdJ3
       ! Cuu
       Cuu(:,:) = 0.0_rp
       if(     kfLoad == opening )then
          ! - Loading
          if ( pdime == 2_ip ) then
             Cuu(1,1) = - dddJ1*K11*J1 + (1.0_rp - d)*K11
             Cuu(1,2) = - dddJ3*K11*J1
             Cuu(2,1) = - dddJ1*K33*J3
             Cuu(2,2) = - dddJ3*K33*J3 + (1.0_rp - d)*K33
          else if (pdime == 3_ip) then
             Cuu(1,1) = - dddJ1*K11*J1 + (1.0_rp - d)*K11
             Cuu(1,2) = - dddJ2*K11*J1
             Cuu(1,3) = - dddJ3*K11*J1
             Cuu(2,1) = - dddJ1*K22*J2
             Cuu(2,2) = - dddJ2*K22*J2 + (1.0_rp - d)*K22
             Cuu(2,3) = - dddJ3*K22*J2
             Cuu(3,1) = - dddJ1*K33*J3
             Cuu(3,2) = - dddJ2*K33*J3
             Cuu(3,3) = - dddJ3*K33*J3 + (1.0_rp - d)*K33
          end if
       elseif( kfLoad == closing ) then
          ! - Unloading
          if (pdime == 2_ip) then
             Cuu(1,1) = - dddJ1*K11*J1  + (1.0_rp - d)*K11
             Cuu(2,2) = KCONT
          else if (pdime == 3_ip) then
             Cuu(1,1) = - dddJ1*K11*J1  + (1.0_rp - d)*K11
             Cuu(2,2) = - dddJ2*K22*J2  + (1.0_rp - d)*K22
             Cuu(3,3) = KCONT
          end if
       endif

    elseif( kfDamage == broken )then
       !
       ! Broken (fully damaged)
       !
       Cuu(:,:) = 0.0_rp
       do idime=1,pdime-1
          Cuu(idime,idime) = Kres
       end do
       if(     kfLoad == opening ) then
          ! - Loading
          Cuu(pdime,pdime) = Kres
       elseif( kfLoad == closing ) then
          ! - Penetration
          Cuu(pdime,pdime) = KCONT
       endif

    endif
    !
    ! ----------------------------------
  end subroutine LAW_COH_TURON_CURRENT_Cuu

  !------------------------------------------------------------------------------
  !> @author  Adria Quintanas and Gerard Guillamet
  !> @date    July, 2018
  !> @brief   Contact law
  !> @details
  !> @note    Kres corresponds to a residual stiffnes in order to overcome
  !>          numerical issues for the global system matrix of the problem
  !>          Testing is required.
  !------------------------------------------------------------------------------
  subroutine CONTACT_LAW_UNILATERAL( &
       itang, pdime, dt, properties, jumps, stateOld, stateNew, tu, Cuu )
    ! ----------------------------------
    implicit none
    ! ----------------------------------
    integer(ip), intent(inout)        :: &
         itang                             ! Tangent calculation method
    integer(ip), intent(in)           :: &
         pdime                             ! Number of dimensions of the problem
    real(rp),    intent(in)           :: &
         dt,                             & ! Time increment
         properties(:),                  & ! Contact law properties
         jumps(:)                          ! Displacement jump (local coordinate system) (dime)
    real(rp),    intent(inout)        :: &
         stateOld(:),                    &
         stateNew(:)
    real(rp),    intent(out)          :: &
         tu(:),                          & ! Tractions / tau (local coordinate system) (dime)
         Cuu(:,:)                          ! Tangent stiffness tensor (dime, dime)
    ! ----------------------------------
    integer(ip)                       :: &
         idime,                          & ! Loading keyflag
         kfLoad
    real(rp)                          :: &
         KCONT, Kres,                    & ! Properties: penalty stiffness
         J3, J3n,                        & ! Effective opening jump
         dtim
    ! ----------------------------------
    !
    !
    ! Variables not used
    stateOld(1:2) = 0.0_rp
    stateNew(1:2) = 0.0_rp
    dtim = dt
    !
    ! Init flags
    !
    kfLoad = opening
    !
    ! Contact properties
    !
    KCONT = properties(1)
    Kres  = properties(2)
    !
    ! Effective opening jump
    J3 = jumps(pdime)
    ! Numerical correction (values close to 0)
    ! Negative (J3n) normal jump
    if( J3 >= tolPe ) then
       J3n = 0.0_rp
       kfLoad = opening
    else
       J3n = min(J3, -tolPe)
       kfLoad = closing
    endif
    !
    ! Tractions (contact)
    tu(:) = 0.0_rp
    if( kfLoad == closing )then
       tu(pdime) = KCONT*J3n
    endif
    !
    ! Tangent stiffness
    if( itang == ELINT_TANGENT_ANALYTICAL .or. & ! Analytical/Secant
         itang == ELINT_TANGENT_SECANT ) then
       !
       Cuu(:,:) = 0.0_rp
       do idime=1,pdime-1
          Cuu(idime,idime) = Kres   ! Leave a meaninful value
       end do
       if( kfLoad == opening ) then
          ! - Loading
          Cuu(pdime,pdime) = Kres   ! Leave a meaninful value
       else
          ! - Penetration
          Cuu(pdime,pdime) = KCONT
       endif

    else                                        ! Explicit/Numerical
       Cuu(:,:) = 0.0_rp

    end if
    !
    ! ----------------------------------
  end subroutine CONTACT_LAW_UNILATERAL

  !------------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    23/04/2018
  !> @brief   Midsurface coordinates and displacement jump for 2d and 3d problems
  !> @details
  !------------------------------------------------------------------------------
  subroutine midcoord_and_djump(&
       pdime, pnode, snode, pevah, eldisp, elcoor, sfcoor, sfjump)
    ! ---------------------------------------
    implicit none
    ! ---------------------------------------
    integer(ip),    intent(in)             :: &
         pdime,                               & ! Dimension
         pnode,                               & ! Number of nodes
         snode,                               & ! Mid surface number of nodes
         pevah                                  ! Midsurface number degrees of freedom
    real(rp),       intent(in)             :: &
         elcoor(pdime,pnode),                 & ! Nodal element coordinates
         eldisp(pdime,pnode)                    ! Nodal element displacement
    real(rp),       intent(out)            :: &
         sfcoor(pdime,snode),                 & ! Midsurface coordinates
         sfjump(pevah)                          ! Midsurface relative displacement (displ jump)
    !----------------------------------------
    integer(ip)                            :: &
         inode, jnode, idime, idof
    ! ---------------------------------------
    !
    sfjump(:)   = 0.0_rp
    sfcoor(:,:) = 0.0_rp
    if ( pdime == 2_ip ) then       ! 2-d
       !
       ! Midsurface definition: micod
       !   - coordinates
       sfcoor(1:pdime,1) = 0.5_rp*(elcoor(1:pdime,1) + elcoor(1:pdime,4) &
            + eldisp(1:pdime,1) + eldisp(1:pdime,4))
       sfcoor(1:pdime,2) = 0.5_rp*(elcoor(1:pdime,2) + elcoor(1:pdime,3) &
            + eldisp(1:pdime,2) + eldisp(1:pdime,3))
       !   - relative displacement / displacement jump
       sfjump(1) = eldisp(1,4) - eldisp(1,1)
       sfjump(2) = eldisp(2,4) - eldisp(2,1)
       sfjump(3) = eldisp(1,3) - eldisp(1,2)
       sfjump(4) = eldisp(2,3) - eldisp(2,2)

    else                            ! 3-d
       !
       ! Midsurface definition: micod
       !   - coordinates
       do inode = 1, snode
          jnode = inode + snode
          do idime = 1, pdime
             sfcoor(idime,inode) = 0.5_rp*(elcoor(idime,inode) + elcoor(idime,jnode) + &
                  eldisp(idime,inode) + eldisp(idime,jnode))
          end do
       end do
       !   - relative displacement / displacement jump
       do inode = 1, snode
          jnode = inode + snode
          do idime = 1, pdime
             idof = pdime*(inode - 1) + idime
             sfjump(idof) = eldisp(idime,jnode) - eldisp(idime,inode)
          end do
       end do

    end if
    !
    !----------------------------------------
  end subroutine midcoord_and_djump

  !------------------------------------------------------------------------------
  !> @author  Adria Quintanas (adria.quintanas@udg.edu)
  !> @date    10/01/2017
  !> @brief   Compute shape function
  !> @details
  !------------------------------------------------------------------------------
  subroutine shape_functions( &
    igaus, pdime, sdime, snode, gpweig, gpshap, gpderi)
    ! ---------------------------------------
    implicit none
    ! ---------------------------------------
    integer(ip),    intent(in)             :: &
      igaus,                                  & ! Integration point identification
      pdime,                                  & ! Dimension
      sdime,                                  & ! Dimension number of the mid-surface
      snode                                     ! Number of nodes of the mid-surface
    real(rp),       intent(out)            :: &
      gpweig,                                 &
      gpshap(snode),                          & ! GP Shape function
      gpderi(sdime,snode)                       ! GP Derivative shape function
    !----------------------------------------
    real(rp)                               :: &
      t, s,                                   & ! Natural coordinates (s, t)
      gpposi(sdime,snode)                       ! Natural coordintaes of the integration points
    ! ---------------------------------------
    !
    ! Natural coordiantes of the integration points position
    !            s                         s
    !  2-d       |            3-d         /
    !            |                4 X ---/--- X 3
    !            |                 /    /    /
    !     1 X- - ----X 2 ---t     /    /_ _ /_ _ t
    !                            /         /
    !                         1 X ------- X 2
    !
    ! Newton Cotes
    !
    ! GP Weigth
    gpweig = 1.0_rp
    !
    gpposi(:,:) = 0.0_rp
    gpshap(:)   = 0.0_rp
    gpderi(:,:) = 0.0_rp
    !
    if ( pdime == 2_ip ) then             ! 2-d
       gpposi(1,1) = -1.0_rp
       gpposi(1,2) =  1.0_rp
       t = gpposi(1,igaus)
       !
       ! GP Shape functions (N)
       gpshap(1) = 0.5_rp*(1.0_rp - t)
       gpshap(2) = 0.5_rp*(1.0_rp + t)
       ! GP Derivative Shape functions w.r.t. natural coordinates (dN/dn)
       gpderi(1,1) = -0.5_rp
       gpderi(1,2) =  0.5_rp

    else if( pdime == 3_ip ) then                                  ! 3-d
       gpposi(1,1) = -1.0_rp
       gpposi(2,1) = -1.0_rp
       gpposi(1,2) =  1.0_rp
       gpposi(2,2) = -1.0_rp
       gpposi(1,3) =  1.0_rp
       gpposi(2,3) =  1.0_rp
       gpposi(1,4) = -1.0_rp
       gpposi(2,4) =  1.0_rp
       t = gpposi(1,igaus)
       s = gpposi(2,igaus)
       !
       ! GP Shape functions (N)
       gpshap(1) = 0.25_rp*(1.0_rp + s*t - t - s)
       gpshap(2) = 0.25_rp*(1.0_rp - s*t + t - s)
       gpshap(3) = 0.25_rp*(1.0_rp + s*t + t + s)
       gpshap(4) = 0.25_rp*(1.0_rp - s*t - t + s)
       ! GP Derivative Shape function w.r.t. natural coordinates (dN/dn)
       gpderi(1,1) = 0.25_rp*(-1.0_rp + s)
       gpderi(1,2) = 0.25_rp*( 1.0_rp - s)
       gpderi(1,3) = 0.25_rp*( 1.0_rp + s)
       gpderi(1,4) = 0.25_rp*(-1.0_rp - s)
       !
       gpderi(2,1) = 0.25_rp*(-1.0_rp + t)
       gpderi(2,2) = 0.25_rp*(-1.0_rp - t)
       gpderi(2,3) = 0.25_rp*( 1.0_rp + t)
       gpderi(2,4) = 0.25_rp*( 1.0_rp - t)
    end if
    !
    !----------------------------------------
  end subroutine shape_functions

  !------------------------------------------------------------------------------
  !> @author  Adria Quintanas (adria.quintanas@udg.edu)
  !> @date    10/01/2017
  !> @brief   Rotation matrix
  !> @details Maths:
  !>            x  = (X^{top} + X^{bot} + q^{top} + q^{bot})
  !>            r1 = 1/2*(dN/dn1)*x
  !>            rt = 1/2*(dN/dn2)*x
  !>            r3 = r1 x rt
  !>            r2 = r3 x r1
  !>            R  = {r1/|r1| , r2/|r2|, r3/|r3|}
  !>          where:
  !>            x are the nodal coordinates at the current configuration
  !>            X are the nodal coordinates at the reference configuration
  !>            q are the nodal displacements at the reference configuration
  !>            N are the shape functions of the 2D interface element, i.e. for
  !>              a 3D 8-nodes cohesive element -> a 2D 4-nodes element.
  !>            (*)^{bot} are the bottom nodes, i.e. 1 - 2 - 3 - 4 (3D 8-nodes)
  !>            (*)^{top} are the top nodes,    i.e. 5 - 6 - 7 - 8 (3D 8-nodes)
  !------------------------------------------------------------------------------
  subroutine rotation_matrix( &
    pdime, sdime, snode, sfcoor, gpderi, gprota, gparea )
    ! ---------------------------------------
    implicit none
    ! ---------------------------------------
    integer(ip),    intent(in)             :: &
      pdime,                                  &  ! Number of dimension
      sdime,                                  &  ! Mid-surface number of dimension
      snode                                      ! Mid-surface number of nodes
    real(rp),       intent(in)             :: &
      sfcoor(pdime,snode),                    &  ! Mid-surface coordinates
      gpderi(sdime,snode)                        ! GP derivative shape function
    real(rp),       intent(out)            :: &
      gprota(pdime,pdime),                    &  ! GP rotation operator
      gparea                                     ! GP area
    !----------------------------------------
    integer(ip)                            :: &
      idime, inode
    real(rp)                               :: &
      r1(pdime), r2(pdime), r3(pdime), rM
    ! ---------------------------------------
    !
    ! Initialize variables
    !
    r1(:) = 0.0_rp
    r2(:) = 0.0_rp
    r3(:) = 0.0_rp
    gprota(:,:) = 0.0_rp
    !
    if (pdime == 2_ip ) then      ! 2-d
       !
       ! Gradient matrix (g)
       do inode = 1, snode
          do idime = 1, pdime
             r1(idime) = r1(idime) + gpderi(1,inode)*sfcoor(idime,inode)
          end do
       end do
       ! Rotation matrix
       gparea = sqrt(r1(1)**2 + r1(2)**2)
       ! Tangential direction
       gprota(:,1) = r1(:)/gparea
       ! Normal direction
       gprota(1,2) = -gprota(2,1)
       gprota(2,2) =  gprota(1,1)
       !
    else                          ! 3-d
       !
       ! Gradient matrix (g)
       do inode = 1, snode
          do idime = 1, pdime
             r1(idime) = r1(idime) + gpderi(1,inode)*sfcoor(idime,inode)
             r2(idime) = r2(idime) + gpderi(2,inode)*sfcoor(idime,inode)
          end do
       end do
       ! Compute normal direction vector
       r3(1) =  r1(2)*r2(3) - r1(3)*r2(2)
       r3(2) =  r1(3)*r2(1) - r1(1)*r2(3)
       r3(3) =  r1(1)*r2(2) - r1(2)*r2(1)
       ! Rotation matrix: 1s tangential and normal direction
       rM = sqrt(r1(1)**2 + r1(2)**2 + r1(3)**2)
       gprota(:,:) = 0.0_rp
       gprota(1,1) = r1(1)/rM
       gprota(2,1) = r1(2)/rM
       gprota(3,1) = r1(3)/rM
       !
       gparea = sqrt(r3(1)**2 + r3(2)**2 + r3(3)**2)
       gprota(1,3) = r3(1)/gparea
       gprota(2,3) = r3(2)/gparea
       gprota(3,3) = r3(3)/gparea
       ! Rotation matrix: 2nd tangential
       gprota(1,2) =  gprota(2,3)*gprota(3,1) - gprota(3,3)*gprota(2,1)
       gprota(2,2) =  gprota(3,3)*gprota(1,1) - gprota(1,3)*gprota(3,1)
       gprota(3,2) =  gprota(1,3)*gprota(2,1) - gprota(2,3)*gprota(1,1)
       !
    end if

  end subroutine rotation_matrix

  !------------------------------------------------------------------------------
  !> @author  Adria Quintanas (adria.quintanas@udg.edu)
  !> @date    12/04/2018
  !> @brief   Compute B matrix
  !> @details
  !------------------------------------------------------------------------------
  subroutine b_matrix_half_positive( &
    pdime, snode, pevah, gpshap, gprota, gpbmat )
    ! ---------------------------------------
    implicit none
    ! ---------------------------------------
    integer(ip),    intent(in)             :: &
      pdime,                                  &  ! Dimension of the problem
      snode,                                  &  ! MS number of nodes
      pevah                                      ! MS number of degree of freedom
    real(rp),       intent(in)             :: &
      gpshap(snode),                          &  ! GP Shape function
      gprota(pdime,pdime)                        ! GP Derivative shape function
    real(rp),       intent(out)            :: &
      gpbmat(pdime,pevah)                        ! GP B-matrix operator ( only positive part )
    ! ---------------------------------------
    integer(ip)                            :: &
      inode, idime, jdime, idofn
    ! ---------------------------------------
    !
    ! Half B matrix
    !
    gpbmat(:,:) = 0.0_rp
    do inode = 1, snode
       do jdime = 1, pdime
          idofn = pdime*(inode - 1) + jdime
          do idime = 1, pdime
             gpbmat(idime,idofn) = gprota(jdime,idime)*gpshap(inode)
          end do
       end do
    end do
    ! ---------------------------------------
    end subroutine b_matrix_half_positive

  !------------------------------------------------------------------------------
  !> @author  Adria Quintanas (adria.quintanas@udg.edu)
  !> @date    10/01/2017
  !> @brief   Compute displacement jump w.r.t local coordinate system
  !> @details Math:
  !>            delta = R^{T} [-N N] q
  !>            delta = B d
  !>          where
  !>             B is the relative displacement matrix. Due to its antysimetry
  !>               it can be composed by the positive and negative part [-B B]
  !>             q is the nodal displacement.
  !>             d is the midsurface nodal relative displacement (d1 = q5 - q1,
  !>               d2 = q6 - q2, d3 = q7 - q3, d4 = q6 - q4)
  !------------------------------------------------------------------------------
  subroutine displacement_jump_local_csys( &
    pdime, pevah, gpbmat, sfjump, gpjump )
    ! ---------------------------------------
    implicit none
    ! ---------------------------------------
    integer(ip),    intent(in)             :: &
      pdime,                                  &  ! Number of dimension
      pevah                                      ! MS degree of freedom
    real(rp),       intent(in)             :: &
      gpbmat(pdime,pevah),                    &  ! GP B-matrix positive operator
      sfjump(pevah)                              ! MS displacement jump wrt global coordinate system
    real(rp),       intent(out)            :: &
      gpjump(pdime)                              ! GP displacement jump wrt local coordinate system
    ! ---------------------------------------
    integer(ip)                            :: &
      idime, idof
    ! ---------------------------------------
    !
    ! Half B matrix
    !
    gpjump(:) = 0.0_rp
    do idime = 1, pdime
      do idof = 1, pevah
        gpjump(idime) = gpjump(idime) + gpbmat(idime,idof)*sfjump(idof)
      end do
    end do
    !
    end subroutine displacement_jump_local_csys

  !------------------------------------------------------------------------------
  !> @author  Adria Quintanas (adria.quintanas@udg.edu)
  !> @date    10/01/2017
  !> @brief   Add the integration point contribution to the internal force vector
  !> @details
  !------------------------------------------------------------------------------
  subroutine internal_force_vector( &
    pdime, pevah, gpbmat, gptrac, gpweig, gparea, elfint )
    ! ---------------------------------------
    implicit none
    ! ---------------------------------------
    integer(ip),    intent(in)             :: &
      pdime,                                  &  ! Dimension
      pevah                                      ! MS degree of freedom
    real(rp),       intent(in)             :: &
      gpbmat(pdime,pevah),                    &  ! GP Bmatrix positive operator
      gptrac(pdime),                          &  ! GP cohesive tractions
      gpweig,                                 &  ! GP weight
      gparea                                     ! GP mapped area (|Jacobian|)
    real(rp),       intent(inout)          :: &
      elfint(pevah)                              ! ELEMENT internal force
    ! ---------------------------------------
    integer(ip)                            :: &
      idime, idof
    ! ---------------------------------------
    !
    ! Internal force vector
    !
    do idof = 1, pevah
      do idime = 1, pdime
        elfint(idof) = elfint(idof) + gpbmat(idime,idof)*gptrac(idime)*gpweig*gparea
      end do
    end do
    !
    ! ---------------------------------------
    end subroutine internal_force_vector

  !------------------------------------------------------------------------------
  !> @author  Adria Quintanas (adria.quintanas@udg.edu)
  !> @date    10/01/2017
  !> @brief   Add the integration point contribution to the tangent stiffness matrix
  !> @details
  !------------------------------------------------------------------------------
  subroutine tangent_stiffness_matrix( &
    pdime, pevah, gpbmat, gpdtdj, gpweig, gparea, elstif)
    ! ---------------------------------------
    implicit none
    ! ---------------------------------------
    integer(ip),    intent(in)             :: &
      pdime,                                  &  ! Dimension
      pevah                                      ! MS degree of freedom
    real(rp),       intent(in)             :: &
      gpbmat(pdime,pevah),                    &  ! GP B-matrix positive operator
      gpdtdj(pdime,pdime),                    &  ! GP tractions vecotr
      gpweig,                                 &  ! GP weight
      gparea                                     ! GP mapped area (|Jacobian|)
    real(rp),       intent(inout)          :: &
      elstif(pevah,pevah)                        ! ELEMENT stiffness matrix
    ! ---------------------------------------
    integer(ip)                            :: &
      idime, jdime, idof, jdof
    real(rp)                               :: &
      auxs1
    ! ---------------------------------------
    !
    ! Tangent stiffness matrix
    !
    auxs1 = gpweig*gparea
    do idof = 1, pevah
      do jdof = 1, pevah
        do idime = 1, pdime
          do jdime = 1, pdime
            elstif(idof,jdof) = elstif(idof,jdof) + gpbmat(idime,jdof)*&
                                gpbmat(jdime,idof)*gpdtdj(idime,jdime)*auxs1
          end do
        end do
      end do
    end do
    ! ---------------------------------------
    end subroutine tangent_stiffness_matrix

    !------------------------------------------------------------------------------
    !> @author  Gerard Guillamet and Adria Quintanas (adria.quintanas@udg.edu)
    !> @date    10/01/2017
    !> @brief   Dirichlet boundary condition correction
    !> @details
    !------------------------------------------------------------------------------
    subroutine nodes_dirichlet_correction( &
         ielem, pdime, pnode, pevat, elmat, elrhs )
      ! ---------------------------------------
      use def_solver, only                   :  &
           solve_sol
      ! ---------------------------------------
      implicit none
      ! ---------------------------------------
      integer(ip), intent(in)                :: &
           ielem,                               & !< ID element
           pdime,                               & !< Number of dimensions
           pnode,                               & !< Number of nodes
           pevat                                  !< Number of degree of freedom
      real(rp), intent(inout)                :: &
           elmat(pevat,pevat),                  & !< Elemental matrix
           elrhs(pevat)                           !< Elemental rhs
      ! ---------------------------------------
      integer(ip)                            :: &
           ipoin, idime, jdime, idof, jdof,     &
           inode, jnode
      real(rp)                               :: &
           adiag
      ! ---------------------------------------
      !
      ! Correction elmat & elrhs according to the Dirichlet boundary conditions
      ! Correction is only done when no iffix (Fixity not imposed by solver)
      !
      if( solve_sol(1) % kfl_iffix == 0 ) then
         do inode = 1,pnode
            ipoin = lnods(inode, ielem)
            do idime = 1, pdime
               if( kfl_fixno_sld(idime,ipoin) > 0 ) then
                  idof = (inode - 1)*pdime + idime
                  adiag = elmat(idof,idof)
                  do jnode = 1, pnode
                     do jdime = 1, pdime
                        jdof = (jnode - 1)*pdime + jdime
                        elmat(idof,jdof) = 0.0_rp
                        elmat(jdof,idof) = 0.0_rp
                     end do
                  end do
                  elmat(idof,idof) = adiag
                  elrhs(idof)      = 0.0_rp
               end if
            end do
         end do
      end if
      !
      ! ---------------------------------------
    end subroutine nodes_dirichlet_correction

    !------------------------------------------------------------------------------
    !> @author  Gerard Guillamet and Adria Quintanas (adria.quintanas@udg.edu)
    !> @date    10/01/2017
    !> @brief   Compute reaction forces at Diritchlet BC
    !> @details
    !------------------------------------------------------------------------------
    subroutine nodes_dirichlet_reaction_forces( &
         ielem, pdime, pnode, pevat, elrhs, elfrx )
      ! ---------------------------------------
      implicit none
      ! ---------------------------------------
      integer(ip), intent(in)                :: &
           ielem,                               & !< ID element
           pdime,                               & !< Number of dimensions
           pnode,                               & !< Number of nodes
           pevat                                  !< Number of degree of freedom
      real(rp),    intent(in)                :: &
           elrhs(pevat)                           !< Elemental right hand side
      real(rp),    intent(out)               :: &
           elfrx(pdime,pnode)                     !< Elemental rection force
      ! ---------------------------------------
      integer(ip)                            :: &
           ipoin, inode, idime, idof
      ! ---------------------------------------
      !
      ! Reaction force at Dirichlet BC
      !
      do inode = 1, pnode
         ipoin = lnods(inode, ielem)
         do idime = 1, pdime
            if( kfl_fixno_sld(idime,ipoin) > 0 ) then
               idof = (inode-1)*pdime + idime
               elfrx(idime,inode) = elrhs(idof)
            end if
         end do
      end do
      ! ---------------------------------------
    end subroutine nodes_dirichlet_reaction_forces

    !------------------------------------------------------------------------------
    !> @author  Gerard Guillamet and Adria Quintanas (adria.quintanas@udg.edu)
    !> @date    10/01/2017
    !> @brief   Compute the stable time increment
    !> @details The stable time increment is:
    !>             dt = sqrt(rho_c/K_c)
    !>          where
    !>             density is a superficial density (mass per unite area)
    !>             penalty_stiffness is the bulk_stiffness/Thickness
    !>          Reference: Abaqus 6.14 Manual, Section 32.5
    !------------------------------------------------------------------------------
    subroutine ELINT_stable_time_increment( &
         ielem, stable_time_inc )
      ! ---------------------------------------
      implicit none
      ! ---------------------------------------
      integer(ip), intent(in)                :: &
           ielem
      real(rp),    intent(out)               :: &
           stable_time_inc                        ! Stable time increment
      ! ---------------------------------------
      integer(ip)                            :: &
           pmate
      real(rp)                               :: &
           penalty_stiffness,                   & ! Penalty stiffness
           density                                ! Superficial density
      ! ---------------------------------------
      !
      ! Read information
      !
      pmate = lmate_sld(ielem)
      density = densi_sld(1,pmate)
      penalty_stiffness = densi_sld(2,pmate)
      !
      ! Stable time increment
      !
      stable_time_inc = sqrt( density/(penalty_stiffness+zeror) )
      ! ---------------------------------------
    end subroutine ELINT_stable_time_increment

 !=============================================================| contains |=====

end module mod_sld_interface_element
