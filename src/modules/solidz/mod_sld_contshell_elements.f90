!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!---------------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_contshell_elements.f90
!> @author  Eva Casoni
!> @date    December, 2017
!>          - Subroutine written
!> @author  Gerard Guillamet
!> @date    January, 2019
!>          - Compativility for explicit and implicit transient problems
!>          - Adds flags for EAS and ANS
!>          - Fix bug material oritentations ortho law.
!>          - Refactoring material law manager and material laws.
!> @date    March, 2019
!>          - Boundary conditions for local axes compatibility in Explicit
!>            and Explicit analysis.
!> @date    March, 2022
!>          - Hyperelastic Neo-Hooken model
!>
!> @brief   Continuum shell elements.
!>
!> @details References:\n
!>          J. Reinoso\n
!>
!> @todo    To do list:\n
!>          - Mass matrix should be calculated following the cshel integral.
!>          - Surface and body force vectors (external forces)
!>          - Cauchy stresses for postrprocess purposes.
!>
!> @}
!---------------------------------------------------------------------------------

module mod_sld_contshell_elements
  ! ==============================================================================
  ! INIT
  !
  use def_kintyp, only                     :  &
       ip, rp, lg
  use def_master, only                     :  &
       ITER_K, leinv_loc
  use def_solidz, only                     :  &
       SLD_IMPLICIT_SCHEME,                   &
       SLD_EXPLICIT_SCHEME,                   &
       SLD_DYNAMIC_PROBLEM,                   &
       SLD_CSHEL_EAS,                         &
       SLD_CSHEL_ANS_SHEAR,                   &
       SLD_CSHEL_ANS_TRAPEZOIDAL,             &
       kfl_cshel_sld,                         &
       kfl_timei_sld,                         & ! Transient or equilibrium problem
       kfl_fixno_sld,                         & ! Dirichelt nodes
       densi_sld,                             & ! Density
       accel_sld,                             & ! Nodal acceleration vector
       parco_sld,                             & ! Cohesive law parameters
       svegm_sld,                             & ! State variables at gauss point per element
       lawco_sld,                             & ! Constitutive law
       dunkn_sld,                             & ! Displacement correction from N-R
       lmate_sld                                ! List of materials
  use def_domain, only                     :  &
       lnods
  use def_solver, only                     :  &
       solve_sol
  use mod_sld_stress_model_comput               ! Include functions to pass to voig notation

  ! -----------------------------------------------------------------------------
  implicit none
  !
  integer(ip),    parameter                   :: &
       pdime = 3,                                & ! Dimensions (2D / 3D)
       pdimh = 2,                                & ! Dimensions surface (2D or 1D)
       pnode = 8,                                & ! Element: number of nodes
       pnodh = 4,                                & ! Mid-surface: number of nodes (elem_node/2)
       pgausT = 2,                               & ! Thickness number of integration points
       pgaus = 2,                                & ! Mid-surface: number of integration points in each direction (elem_gauss/2)
       pdofs = 24,                               & ! Element: number of degrees of freedom (elem_node*dof_node)
       pdofh = 12                                  ! Mid-surface: number of degrees of freedom (elem_dof/2)
  !
  !=============================================================| init |=========
  ! PUBLIC
  !
  public                                   :: &
       SHELL_elemental_operations
  !=============================================================| public |=======
  !==============================================================================
  ! PRIVATE
  !
  private                                  :: &
       shape_functions,                          &
       a_midsurface_basis,                       &
       g_shell_basis,                            &
       g_shell_basis2,                           &
       b_matrix,                                 &
       s8vthv,                                   &
       material_law,                             &
       integrate_material_law,                   &
       cshell_stress_model_100,                  &
       cshell_stress_model_151,                  &
       compute_matrices,                         &
       compute_EAS_force,                        &
       computeStiffMatrix,                       &
       compute_mass_matrix,                      &
       internalForceVector,                      &
       inertialContribution,                     &
       EASM_local,                               &
       EASM_transform,                           &
       compute_enhanced_strains,                 &
       ANScolocationpointsShearQ,                &
       ANScolocationpointsShearT,                &
       ANS_bmatrix_shearLocking,                 &
       ANS_bmatrix_trapezoidalLocking,           &
       K_element_reorder,                        &
       b_matrixANS_shearLocking,                 &
       b_matrixANS_trapezoidalLocking,           &
       computeGeomStiffMatrixGlobal,             &
       ForceVectorReord,                         &
       s9Creordering,                            &
       s81CurvTransf81,                          &
       nodes_dirichlet_correction,               &
       nodes_dirichlet_reaction_forces,          &
       sld_local_axes_rhs_and_amatr

  !
  !=============================================================| private |======
  !==============================================================================
  ! CONTAINS
  !
contains
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !------------------------------------------------------------------------------
  !> @author  Eva Casoni (eva.casoni@bsc.es)
  !> @date    218/09/2016
  !> @brief   Elemental operations
  !> @details
  !------------------------------------------------------------------------------
  subroutine SHELL_elemental_operations(itask, ielem, elcoor, eldisp, &
       elmas, elmuu, elfext, elfint, elfine, elrhs, elfrx, elmat, &
       pgausMass, gpvol, gpsha) !< PROVISIONAL: Eliminar el dia que calculem la mass matrix ben fet

    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    integer(ip),    intent(in)                :: &
         itask,                                  & !< Itask for time treatment
         ielem,                                  & !< Element number
         pgausMass                                 !>>> Eliminar el dia que calculem la mass matrix ben fet
    real(rp),       intent(in)                :: &
         elcoor(pdime,pnode),                    & !< Nodal element coordinates (dime, node)
         eldisp(pdime,pnode),                    & !< Nodal element displacement (dime, node)
         gpvol(pgausMass),                       & !< PROVISIONAL: Eliminar el dia que calculem la mass matrix ben fet
         gpsha(pnode,pgausMass)                    !< PROVISIONAL: Eliminar el dia que calculem la mass matrix ben fet
    real(rp),       intent(out)               :: &
         elfext(pdofs),                          & !< Vector of external forces
         elfint(pdofs),                          & !< Vector of internal forces
         elfine(pdofs),                          & !< Vector of inertial forces
         elmas(pnode,pnode),                     & !< Consistent mass matrix
         elmuu(pnode),                           & !< Lumped mass matrix
         elrhs(pdofs),                           & !< Elemental RHS vector (node*dof)
         elfrx(pdime,pnode),                     & !< Elemental reaction forces vector (node*dof)
         elmat(pdofs,pdofs)                        !< Elemental tangent matrix (node*dof, node*dof)
    ! ---------------------------------------------------------------------------
    integer(ip)                               :: &
         ipoin, jpoin, idime, igaus, jgaus,      &  !
         igausT,                                 &  ! counters
         i, j, inode, itott,                     &  !
         pmate,                                  &
         flag_EAS,                               &  ! Flag for EAS method
         flag_SL, flag_TL                           ! Flags for ANS method (shear and trapezoidal)
    real(rp)                                  :: &
         posgp(pgaus),                           &  ! Gauss points in the midplane RS
         posgpT(pgausT),                         &  ! Gauss points in the thickness
         wgp(pgaus),                             &  ! Weight of gauss points
         wgpT(pgausT)                               ! Weight of the thickness

    real(rp)                               :: &
                                ! Geometrical variables
         elddisp(pdime,pnode),                   & ! Increment of displacement at element level
         elcoorc(pdime,pnode),                   & ! Element current coordinates
         sfref(pdime,pnodh),                     & ! Midsurface coordinates in reference configuration (dime, mid_node)
         sfcur(pdime,pnodh),                     & ! Midsurface coordinates in current configuration
         a3r(pdime,pnodh),                       & ! Normal vector to the midplane in reference coordinates
         a3c(pdime,pnodh),                       & ! Normal vector to the midplane in current coordinates
                                ! Variables at the element center (evaluated at (0,0))
         shap0(pnodh),                           & ! Shape functions evaluated at midpoint of the shell
         dshap0(pdimh,pnodh),                    & ! Derivatives of shape functions evaluated at midpoint of the shell
         a0rkov(pdime,pdime),                    & !
         a0rkon(pdime,pdime),                    & !
         am0rkov(pdime,pdime),                   & !
         am0rkon(pdime,pdime),                   & !
         det0r,                                  & !
                                ! Variables of the midsurface
         gpshap(pnodh),                          & ! Shape functions of the midsurface nodes
         gpderi(pdimh,pnodh),                    & ! Derivative of the shape function dN/dn at nodes dNi/dxi1, dNi/dxi2
         arkov(pdime,pdime),            & ! Curvilinear basis (covariant vectors) on the shell midsurface in reference configuration
         ackov(pdime,pdime),            & ! Curvilinear basis (covariant vectors) on the shell midsurface in current configuration
         a3rkovder(pdime,2),            & ! Derivatives of a3 in reference coordinates (of third verctor of curvilinear basis)
         a3ckovder(pdime,2),            & ! Derivatives of a3 in current coordinates (of third verctor of curvilinear basis)
         detra,                                  & ! Determinant of curvilinear basis in reference configuration
         detca,                                  & ! Determinant of curvilinear basis in current configuration
                                ! Variables of the shell (with curvadure)
         grkov(pdime,pdime),                     & ! Curvilinear basis (covariant vectors) on the shell in reference configuration
         grkon(pdime,pdime),                     & !
         gckov(pdime,pdime),                     & ! Curvilinear basis (covariant vectors) on the shell in current configuration
         gckon(pdime,pdime),                     & !
         gmrkov(pdime,pdime),                    & ! Metric of shell in reference configuration
         gmckov(pdime,pdime),                    & ! Metric of shell in current configuration
         gmrkon(pdime,pdime),                    & ! Inverse of gmrkov
         gmckon(pdime,pdime),                    & ! Inverse of gmckon
         detrg,                                  & ! Determinant of curvilinear basis in reference configuration
         detcg,                                  & ! Determinant of curvilinear basis in current configuration
                                ! Material law variables
         vodds(6,6),                             & !
         stress(6),                              & !
         gpddsVoig(12,12),                       & ! Constitutive tensor in Voig notation
         gpstrVoig(12),                          & ! 2PK in Voig notation
                                ! Magnitudes due to EAS method
         Tmat(12,12),                            & ! Transformation matrix s.t: Mglo = T*Mloc
         Mloc(12,22),                            & ! Transformation matrix in local coordinates
         Mglo(12,22),                            & ! Matrix for enchanced strains: Etilde = Mglo*alpha
         DEAS(22,22),                            & ! DEAS = M^t*gpddsVo*M = k_alphaalpha
         DEASinv(22,22),                         & ! Inverse of DEAS
         LEAS(22,24),                            & ! LEAS = M^t*gpddsVo*B = k_alphau
         LEASt(24,22),                           & ! Transpose of LEAS
         fEAS(22),                               & ! fEAS = M* t*gpstrVo
         Etilde(12),                             & ! Enhanced strains: Etilde = Mglo*alpha
         alpha(22),                              & ! Enhancing strain vector
         sumepsilon,                             & ! Auxiliar
                                ! Magnitudes due to ANS method
                                ! alpha13
         gpshap1q(2,1,pnodh),                    &
         gpderi1q(2,pdimh,pnodh),                &
         ackov1q(2,pdime,pdime),                 & ! Current
         a3kvpc1q(2,pdime,pdime),                &
                                ! alpha 23
         gpshap2q(2,1,pnodh),                    &
         gpderi2q(2,pdimh,pnodh),                &
         ackov2q(2,pdime,pdime),                 & ! Current
         a3kvpc2q(2,pdime,pdime),                &
                                ! alpha 33 (curvature thickness locking)
         gpshapT(4,1,pnodh),                     &
         ackovT(4,pdime,pdime),                  & ! Current

                                ! B matrix variables
         frQ(2),                                 &
         fsQ(2),                                 &
         frT(pnodh),                             &
                                ! Other magnitudes
         stiff(pdofs,pdofs),                     & !
         stiffGeom(pdofs,pdofs),                 & !
         stiffTotal(pdofs,pdofs),                & !
         elstif(pdofs,pdofs),                    &
         fint(pdofs),                            & !
         fTotal(pdofs),                          & !
         auxMat(pdofs,pdofs),                    & !
         auxK(22,pdofs),                         & !
         auxF(22),                               &
         auxFor(pdofs),                          & !
                                ! Other magnitudes
         gpbmat(pdofh,pdofs),                 & ! B matrix at integration point: E^u = Bu  (E^u compatible strains, u: displacement)
         xi1,xi2,xi3,                         & ! Quadrature
         weight1,weight2,wT,                  & ! Quadrature
         hi(pdime),                           & !
         hnormi
    !
    ! ---------------------------------------------------------------------------
    !
    ! Initialize variables
    !
    elrhs(:)   = 0.0_rp
    elfint(:)  = 0.0_rp
    elfext(:)  = 0.0_rp
    elfine(:)  = 0.0_rp
    elstif(:,:)= 0.0_rp
    elmat(:,:) = 0.0_rp
    elmuu(:)   = 0.0_rp
    elmas(:,:) = 0.0_rp
    !
    pmate    = lmate_sld(ielem)
    flag_EAS = kfl_cshel_sld(1)
    flag_SL  = kfl_cshel_sld(2)
    flag_TL  = kfl_cshel_sld(3)
    !
    ! Compute current coordinates
    !
    elcoorc(:,:) = 0.0_rp
    do ipoin = 1,pnode
       do idime = 1,pdime
          elcoorc(idime,ipoin) = elcoor(idime,ipoin) + eldisp(idime,ipoin)
       end do
    end do

    !
    ! Midsurface definition in reference and updated configuration
    !   - coordinates
    sfref(:,:) = 0.0_rp
    sfcur(:,:) = 0.0_rp
    do ipoin = 1, pnodh
       jpoin = ipoin + pnodh
       do idime = 1, pdime
          sfref(idime,ipoin) = 0.5_rp*(elcoor(idime,ipoin) + elcoor(idime,jpoin))
          sfcur(idime,ipoin) = 0.5_rp*(elcoorc(idime,ipoin) + elcoorc(idime,jpoin))
       end do
    end do

    !
    ! Compute normal vector to the midsurface in reference and updated configuration
    !
    a3r(:,:) = 0.0_rp
    a3c(:,:) = 0.0_rp
    do ipoin = 1,pnodh
       do idime = 1,pdime
          a3r(idime,ipoin) = 0.5_rp*(elcoor(idime,ipoin+pnodh) - elcoor(idime,ipoin))
          a3c(idime,ipoin) = 0.5_rp*(elcoorc(idime,ipoin+pnodh) - elcoorc(idime,ipoin))
       end do
    end do
    !
    ! EAS method
    !
    if ( flag_EAS == SLD_CSHEL_EAS .and. itask == SLD_IMPLICIT_SCHEME ) then
       !
       ! Get nodal increment displacement
       !
       do inode = 1,pnode
          ipoin = lnods(inode,ielem)
          itott = (ipoin-1)*pdime
          do idime = 1, pdime
             itott = itott + 1
             elddisp(idime,inode) = dunkn_sld(itott)
          end do
       end do
       !
       ! Computations over the midpoint of the shell (gausspoint = (0,0)): curvilinear basis in reference
       ! and current configuration
       !
       call shape_functions(0.0_rp, 0.0_rp, shap0(:), dshap0(:,:))
       call g_shell_basis(sfref(:,:), a3r(:,:), shap0(:), dshap0(:,:), 0.0_rp, &
            a0rkov(:,:), a0rkon(:,:), am0rkov(:,:), am0rkon(:,:), det0r, ielem)
       !
       ! Comute enhanced strain (alpha) previous to the iteration i:
       ! only solve part of the system: not the full static condensation system
       !
       call compute_enhanced_strains(ielem,flag_SL,flag_TL,elcoor,elcoorc,elddisp,alpha)
    end if
    !
    ! ANS Shear Locking
    !
    if ( flag_SL == SLD_CSHEL_ANS_SHEAR ) then
       call ANScolocationpointsShearQ(ielem,sfcur(:,:),a3c(:,:),            &
            gpshap1q(:,:,:),gpderi1q(:,:,:),ackov1q(:,:,:),a3kvpc1q(:,:,:), &
            gpshap2q(:,:,:),gpderi2q(:,:,:),ackov2q(:,:,:),a3kvpc2q(:,:,:))
    end if
    !
    ! ANS Trapezoidal Locking
    !
    if ( flag_TL == SLD_CSHEL_ANS_TRAPEZOIDAL ) then
       call ANScolocationpointsShearT(ielem,sfcur(:,:),a3c(:,:), &
            gpshapT(:,:,:),ackovT(:,:,:))
    end if
    !
    ! Quadrature
    !
    call quadrature_contshell(pgaus,pgausT,posgp,wgp,posgpT,wgpT)
    !
    ! Initializations
    !
    DEAS = 0.0_rp
    LEAS = 0.0_rp
    fEAS = 0.0_rp
    fTotal = 0.0_rp
    stiffTotal = 0.0_rp
    !
    ! Loop on gauss points of the surface
    !
    do igaus = 1,pgaus
       xi1 = posgp(igaus)
       weight1 = wgp(igaus)
       do jgaus = 1,pgaus
          xi2 = posgp(jgaus)
          weight2 = wgp(jgaus)
          !
          ! Interpolating function at the middle plane and its derivative (N and dN/dn)
          !
          call shape_functions(xi1, xi2,  gpshap(:), gpderi(:,:))
          !
          ! Basis in the midsurface of the shell in reference and current configuration
          !
          call a_midsurface_basis(sfref(:,:), gpshap(:) ,gpderi(:,:), a3r(:,:), &
               arkov(:,:), a3rkovder(:,:), detra, ielem)
          call a_midsurface_basis(sfcur(:,:), gpshap(:), gpderi(:,:), a3c(:,:), &
               ackov(:,:), a3ckovder(:,:), detca, ielem)
          !
          ! Compute H norm
          !
          hnormi = 0.0_rp
          hi(1) = arkov(2,1)*arkov(3,2) - arkov(3,1)*arkov(2,2)
          hi(2) = arkov(3,1)*arkov(1,2) - arkov(1,1)*arkov(3,2)
          hi(3) = arkov(1,1)*arkov(2,2) - arkov(2,1)*arkov(1,2)
          do idime=1, pdime
             hnormi = hnormi + hi(idime)*hi(idime)
          end do
          hnormi = sqrt(hnormi)
          if (hnormi > 0.0_rp) then
             hi = hi/hnormi
          else
             write(*,*) 'Norm of the thickness is 0'
          end if

          if ( flag_EAS == SLD_CSHEL_EAS .and. itask == SLD_IMPLICIT_SCHEME ) then
            !
            ! M matrix (EAS method) in the local coordinate system (M_xi in order to obtain enhanced strains)
            !
            call EASM_local(xi1, xi2, Mloc(:,:))
            !
            ! Transforming to the global coordinate system (M = (detJ0/detJ)*(T0^-t)*M_xi
            !
            call EASM_transform(Mloc(:,:), Mglo(:,:), arkov(:,:), a0rkon(:,:), detra, det0r, Tmat(:,:))
            !
            ! Compute compatible strains: Etilde = M*alpha
            !
            do i = 1, 12
               sumepsilon = 0.0_rp
               do j = 1, 22
                  sumepsilon = sumepsilon + Mglo(i,j)*alpha(j)
               end do
               Etilde(i) = sumepsilon
            end do
          end if
          !
          ! B matrix (B)
          !
          call b_matrix(gpshap(:), gpderi(:,:), ackov(:,:), a3ckovder(:,:), gpbmat(:,:))
          !
          ! Modify B matrix Shear and Trapezoidal Lockings
          !
          if ( flag_SL == SLD_CSHEL_ANS_SHEAR ) then
             call ANS_bmatrix_shearLocking(xi1,xi2,frQ(:),fsQ(:))
             call b_matrixANS_shearLocking(&
                  gpshap1q(:,:,:),gpderi1q(:,:,:),ackov1q(:,:,:), &
                  gpshap2q(:,:,:),gpderi2q(:,:,:),ackov2q(:,:,:), &
                  frQ,fsQ,gpbmat(:,:))
          end if
          if ( flag_TL == SLD_CSHEL_ANS_TRAPEZOIDAL ) then
             call ANS_bmatrix_trapezoidalLocking(xi1,xi2,frT(:))
             call b_matrixANS_trapezoidalLocking(gpshapT(:,:,:),ackovT(:,:,:),&
                  frT(:),gpbmat(:,:))
          end if
          !
          ! Loop over the interpolating points along the thickness
          !
          gpddsVoig = 0.0_rp
          gpstrVoig = 0.0_rp
          gausspointsthickness: do igausT = 1, pgausT
             xi3 = posgpT(igausT)
             wT = wgpT(igausT)
             !
             ! Curvilinear basis in the shell
             !
             call g_shell_basis(sfref(:,:), a3r(:,:), gpshap(:), gpderi(:,:), xi3, &
                  grkov(:,:), grkon(:,:), gmrkov(:,:), gmrkon(:,:), detrg, ielem)
             call g_shell_basis(sfcur(:,:), a3c(:,:), gpshap(:), gpderi(:,:), xi3, &
                  gckov(:,:), gckon(:,:), gmckov(:,:), gmckon(:,:), detcg, ielem)
             !
             ! Change the metrics due to EAS
             !
             if ( flag_EAS == SLD_CSHEL_EAS .and. itask == SLD_IMPLICIT_SCHEME ) call s8vthv(gmckov(:,:), Etilde, xi3)
             !
             ! Constitutive relationship tensors
             !
             call material_law(ielem, gmrkov(:,:), gmckov(:,:), gmrkon(:,:), gmckon(:,:), &
                  grkon(:,:), detrg, detcg, vodds(:,:), stress(:))
             !
             ! Integrate through the thickness
             !
             wT = wT*detrg/hnormi
             call integrate_material_law(vodds(:,:),stress(:),wT,xi3,gpddsVoig(:,:),gpstrVoig(:))

          end do gausspointsthickness
          !
          ! Internal force vector
          !
          call internalForceVector(gpbmat(:,:), gpstrVoig(:), weight1, weight2, hnormi, fint(:))
          fTotal = fTotal + fint
          !
          ! Stiffness matrix
          !
          if ( itask == SLD_IMPLICIT_SCHEME ) then
             if ( flag_EAS == SLD_CSHEL_EAS ) then
                !
                ! Computation of auxiliar matrices: k_alphaalpha=D=M^t*C*M , k_alphau=L=M^t*C*B
                !
                call compute_matrices(gpddsVoig(:,:),gpbmat(:,:),Mglo(:,:),weight1,weight2,hnormi,DEAS(:,:),LEAS(:,:))
                !
                ! Compute external forces: f_EAS = int(M^t*S)
                !
                call compute_EAS_force(gpstrVoig(:),Mglo(:,:),weight1,weight2,hnormi,fEAS)
             end if
             !
             ! Compute stiffness matrix
             !
             call computeStiffMatrix(gpbmat(:,:), gpddsVoig(:,:), weight1, weight2, hnormi, stiff(:,:))
             !
             ! Geometrical stiffness with/without EAS and ANS
             !
             call computeGeomStiffMatrixGlobal(gpstrVoig(:),gpshap(:),gpderi(:,:),weight1,weight2,hnormi, &
                  gpshap1q(:,:,:),gpderi1q(:,:,:),gpshap2q(:,:,:),gpderi2q(:,:,:),gpshapT(:,:,:),         &
                  frQ,fsQ,frT(:),flag_SL,flag_TL,stiffGeom(:,:))

             stiffTotal =  stiffTotal + stiff + stiffGeom
          end if

       end do
    end do  ! end do gauss points of the surface

    !--------------------------------------------------------------------------
    !
    ! Static condensation due to EAS
    !
    if ( flag_EAS == SLD_CSHEL_EAS .and. itask == SLD_IMPLICIT_SCHEME ) then
       !
       ! Compute DEAS^(-1) = k_alphaalpha^-1 , size 22x22
       !
       ! Invert matrix DEAS and transpose LEAS
       DEASinv = DEAS
       call invert(DEASinv,22_ip,22_ip)
       LEASt = transpose(LEAS)
       !
       ! Auxiliar matrix and vector due to static condensation
       auxK   = matmul(DEASinv,LEAS)
       auxF   = matmul(DEASinv,fEAS)
       auxMat = matmul(LEASt,auxK)
       auxFor = matmul(LEASt,auxF)
       !
       ! Construct stiffness matrix and force vector
       ! fint = fint - k_ualpha*(k_alphalpha^-1*f_EAS)
       ! k_uu = k_uu - k_ualpha*(k_alphalpha^-1*k_alphau)
       fTotal(:)       = fTotal(:) - auxFor(:)
       stiffTotal(:,:) = stiffTotal(:,:) - auxMat(:,:)

    end if
    !
    ! Reordering K and f to be compatible with ALYA
    !
    if ( itask == SLD_IMPLICIT_SCHEME ) call K_element_reorder(stiffTotal,elstif)
    call ForceVectorReord(fTotal,elfint)
    !
    ! Mass matrices and inertial contributions
    !
    if ( kfl_timei_sld == SLD_DYNAMIC_PROBLEM ) then
       !
       ! Mass matrices
       !
       call compute_mass_matrix(pmate,pnode,pgausMass,gpvol,gpsha,elmas,elmuu)

       if ( itask == SLD_IMPLICIT_SCHEME ) then
          !
          ! Inertial contribution
          !
          call inertialContribution(pdime,pnode,ielem,elmas,elfine,elstif)

       end if
    end if
    !
    ! Complete vectors and matrix:
    !
    ! - residual vector (RHS = Fext - Fint - Ma)
    elrhs = elfext - elfint - elfine
    !
    ! - tangent stiffnes matrix
    ! Make it symmetric
    if ( itask == SLD_IMPLICIT_SCHEME ) elmat = 0.5_rp*(elstif + transpose(elstif))
    !
    ! Local Boundary conditions
    !
    call sld_local_axes_rhs_and_amatr(itask, ielem, pnode, pdime, elrhs, elmat)
    !
    ! Save reaction force
    !
    call nodes_dirichlet_reaction_forces(ielem, elrhs(:), elfrx(:,:))
    !
    ! Dirichlet BC correction
    !
    if ( itask == SLD_IMPLICIT_SCHEME ) call nodes_dirichlet_correction(ielem, elmat(:,:), elrhs(:))

  end subroutine SHELL_elemental_operations

  !------------------------------------------------------------------------------
  !> @author  Eva Casoni (eva.casoni@bsc.es)
  !> @date    13/09/2016
  !> @brief   Compute shape functions at middle surface (2D)
  !> @details
  !------------------------------------------------------------------------------
  subroutine shape_functions(r, s, gpshap, gpderi)

    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    real(rp),    intent(in)                :: &
         r,s                                     ! Integration point coordinates
    real(rp),       intent(out)            :: &
         gpshap(pnodh),                       &  ! Shape function at integration point (gauss)
         gpderi(pdimh,pnodh)                     ! Derivative shape function at integration point (nat_dime, node)

    !
    ! ---------------------------------------------------------------------------
    ! Natural coordiantes of the integration points position
    !                 s
    !                /
    !        2 X ---/--- X 1
    !        /     /    /
    !       /     /_ _ / _ _ r
    !      /          /
    !  3 X -------- X 4
    !
    gpshap(1) = 0.25_rp*(1.0_rp + s*r - r - s)
    gpshap(2) = 0.25_rp*(1.0_rp - s*r + r - s)
    gpshap(3) = 0.25_rp*(1.0_rp + s*r + r + s)
    gpshap(4) = 0.25_rp*(1.0_rp - s*r - r + s)
    !
    ! Derivative of shape function w.r.t. natural coordinates (dN/dn)
    !
    gpderi(1,1) = 0.25_rp*(-1.0_rp + s)
    gpderi(1,2) = 0.25_rp*(1.0_rp - s)
    gpderi(1,3) = 0.25_rp*(1.0_rp + s)
    gpderi(1,4) = 0.25_rp*(-1.0_rp - s)
    !
    gpderi(2,1) = 0.25_rp*(-1.0_rp + r)
    gpderi(2,2) = 0.25_rp*(-1.0_rp - r)
    gpderi(2,3) = 0.25_rp*(1.0_rp + r)
    gpderi(2,4) = 0.25_rp*(1.0_rp - r)
    !
  end subroutine shape_functions

  !------------------------------------------------------------------------------
  !> @author  Eva Casoni (eva.casoni@bsc.es)
  !> @date    13/09/2016
  !> @brief   Compute curvilinear basis
  !> @details
  !------------------------------------------------------------------------------
  subroutine a_midsurface_basis(sfcoo,gpsha,gpderi,a3,akov,a3kovder,detakov,jelem)

    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    real(rp),       intent(in)             :: &
         sfcoo(pdime,pnodh),                     &  ! Midsurface coordinates
         gpsha(pnodh),                           &  ! Shape function at integration point (gauss)
         gpderi(pdimh,pnodh),                    &  ! Derivative of shape function at integration point
         a3(pdime,pnodh)
    integer(ip),    intent(in)             :: &
         jelem
    real(rp),       intent(out)            :: &
         akov(pdime,pdime),                      &  ! Covariant basis of the midsurface
         a3kovder(pdime,2),                      &  ! Derivatitves of ackov(:,3)=a3
         detakov                                    ! Determinant of contravariant basis
    ! ---------------------------------------------------------------------------
    integer(ip)                            :: &
         idime, inode

    !
    ! ---------------------------------------------------------------------------
    ! Initialization
    !
    akov = 0.0_rp
    detakov = 0.0_rp
    a3kovder = 0.0_rp

    !
    ! Covariant basis at the shell middle surface in the plane rs
    !
    do idime = 1,pdime
       do inode = 1,pnodh
          akov(idime,1) = akov(idime,1) + gpderi(1,inode)*sfcoo(idime,inode)
          akov(idime,2) = akov(idime,2) + gpderi(2,inode)*sfcoo(idime,inode)
       end do
    end do

    do idime = 1,pdime
       do inode=1,pnodh
          akov(idime,3) = akov(idime,3) + gpsha(inode)*a3(idime,inode)
       end do
    end do

    !
    ! Determinant of the covariant basis
    !
    detakov = akov(1,1)*akov(2,2)*akov(3,3) + akov(1,3)*akov(2,1)*akov(3,2) + &
         akov(3,1)*akov(1,2)*akov(2,3) - akov(3,1)*akov(2,2)*akov(1,3) - &
         akov(3,3)*akov(1,2)*akov(2,1) - akov(1,1)*akov(2,3)*akov(3,2)

    if (detakov < 1.0e-15_rp)  then
       print*, '----> Shell element error: negative determinant in a_midsurface_element in element', jelem
       print*, '----> Value: ', detakov
       call runend("MOD_SLD_CONTSHELL_ELEMENTS: NEGATIVE DETERMINANT IN A_MIDSURFACE_ELEMENT")
    end if

    !
    ! Derivatives of a3 in current coordinates
    !
    do idime = 1,pdime
       do inode = 1,pnodh
          a3kovder(idime,1) = a3kovder(idime,1) + gpderi(1,inode)*a3(idime,inode)
          a3kovder(idime,2) = a3kovder(idime,2) + gpderi(2,inode)*a3(idime,inode)
       end do
    end do

  end subroutine a_midsurface_basis

  !------------------------------------------------------------------------------
  !> @author  Eva Casoni (eva.casoni@bsc.es)
  !> @date    13/09/2016
  !> @brief   Compute curvilinear basis in continuum shell (3D), not in the midsurface
  !> @details
  !------------------------------------------------------------------------------
  subroutine g_shell_basis(sfcoo,a3,gpsha,gpderi,gpT,gkov,gkon,gmkov,gmkon,detg,jelem)

    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    real(rp),       intent(in)             :: &
         sfcoo(pdime,pnodh),                     &  ! Coordinates of midsurface
         a3(pdime,pnodh),                        &  ! Normal to midplane
         gpsha(pnodh),                           &  ! Shape function at integration point (gauss)
         gpderi(pdimh,pnodh),                    &  ! Derivative of shape function at integration point
         gpT                                        ! Coordinate of the gauss point along thickness (xi3)
    integer(ip),    intent(in)             :: &
         jelem
    real(rp),       intent(out)            :: &
         gkov(pdime,pdime),                      &  ! Covariant basis
         gkon(pdime,pdime),                      &  ! Inverse of gkov
         gmkov(pdime,pdime),                     &  ! Metrics of covariant basis
         gmkon(pdime,pdime),                     &  ! Inverse of gmkov
         detg                                       ! Determinant of gkov
    ! ---------------------------------------------------------------------------
    integer(ip)                            :: &
         inode, idime, jdime, kdime
    real(rp)                               :: &
         detgmkov,                            &
         hmlay,thickness,zeta,thicklay
    !
    !---------------------------------------------------------------------------------
    !
    ! Initialization
    !
    gkov = 0.0_rp
    gkon = 0.0_rp
    gmkov = 0.0_rp
    gmkon = 0.0_rp
    detg = 0.0_rp

    hmlay = 0.131_rp      ! height actual thickness
    thickness = 0.131_rp ! is the total element thickness
    thicklay = 0.131_rp  ! ply thickness

    !zeta  = -1.0_rp + (-thicklay*(1.0_rp - gpT) + 2.0_rp*hmlay)/thickness
    zeta = gpT

    !
    ! Covariant basis at the midsurface plane of the shell
    !
    do idime = 1,pdime
       do inode = 1,pnodh
          gkov(idime,1) = gkov(idime,1) + gpderi(1,inode)*sfcoo(idime,inode) + zeta*gpderi(1,inode)*a3(idime,inode)
          gkov(idime,2) = gkov(idime,2) + gpderi(2,inode)*sfcoo(idime,inode) + zeta*gpderi(2,inode)*a3(idime,inode)
       end do
    end do

    do idime = 1,pdime
       do inode=1,pnodh
          gkov(idime,3) = gkov(idime,3) + gpsha(inode)*a3(idime,inode)
       end do
    end do
    !
    ! Determinant of the covariant basis
    !
    detg  = gkov(1,1)*gkov(2,2)*gkov(3,3) + gkov(1,3)*gkov(2,1)*gkov(3,2) + &
            gkov(3,1)*gkov(1,2)*gkov(2,3) - gkov(3,1)*gkov(2,2)*gkov(1,3) - &
            gkov(3,3)*gkov(1,2)*gkov(2,1) - gkov(1,1)*gkov(2,3)*gkov(3,2)

    if (detg < 1.0e-15_rp)  then
       print*, '----> Shell element error: negative determinant in g_shell_basis in element', jelem
       print*, '----> Value: ', detg
       call runend("MOD_SLD_CONTSHELL_ELEMENTS: NEGATIVE DETERMINANT IN G_SHELL_BASIS")
    end if
    !
    ! Inverse of gkov
    !
    if (detg > 0.0_rp) then
       gkon(1,1) =  (gkov(2,2) *gkov(3,3) -gkov(3,2) *gkov(2,3))/detg
       gkon(1,2) = -(gkov(1,2) *gkov(3,3) - gkov(1,3) *gkov(3,2))/detg
       gkon(1,3) =  (gkov(1,2) *gkov(2,3) - gkov(2,2) *gkov(1,3))/detg
       gkon(2,1) = -(gkov(2,1) *gkov(3,3) - gkov(3,1) *gkov(2,3))/detg
       gkon(2,2) =  (gkov(1,1) *gkov(3,3) - gkov(1,3) *gkov(3,1))/detg
       gkon(2,3) = -(gkov(1,1) *gkov(2,3) - gkov(2,1) *gkov(1,3))/detg
       gkon(3,1) =  (gkov(2,1) *gkov(3,2) - gkov(3,1) *gkov(2,2))/detg
       gkon(3,2) = -(gkov(1,1) *gkov(3,2) - gkov(3,1) *gkov(1,2))/detg
       gkon(3,3) =  (gkov(1,1) *gkov(2,2) - gkov(1,2) *gkov(2,1))/detg
    end if

    gkon = Transpose(gkon)

    !
    ! Metric coefficients
    !
    do idime = 1,pdime
       do jdime = 1,pdime
          gmkov(idime,jdime) = 0.0_rp
          do kdime = 1,pdime
             gmkov(idime,jdime) = gmkov(idime,jdime) + gkov(kdime,idime)*gkov(kdime,jdime)
          end do
       end do
    end do
    ! make it symmetric
    gmkov = 0.5_rp*(gmkov + transpose(gmkov))
    !
    ! Compute the inverse gmkon
    !
    detgmkov = gmkov(1,1)*gmkov(2,2)*gmkov(3,3) + gmkov(1,3)*gmkov(2,1)*gmkov(3,2) + &
               gmkov(3,1)*gmkov(1,2)*gmkov(2,3) - gmkov(3,1)*gmkov(2,2)*gmkov(1,3) - &
               gmkov(3,3)*gmkov(1,2)*gmkov(2,1) - gmkov(1,1)*gmkov(2,3)*gmkov(3,2)

    if (detgmkov > 0.0_rp) then
       gmkon(1,1) = (gmkov(2,2)*gmkov(3,3) - gmkov(3,2)*gmkov(2,3))/detgmkov
       gmkon(1,2) = -(gmkov(1,2)*gmkov(3,3) - gmkov(1,3)*gmkov(3,2))/detgmkov
       gmkon(1,3) = (gmkov(1,2)*gmkov(2,3) - gmkov(2,2)*gmkov(1,3))/detgmkov
       gmkon(2,1) = -(gmkov(2,1)*gmkov(3,3) - gmkov(3,1)*gmkov(2,3))/detgmkov
       gmkon(2,2) = (gmkov(1,1)*gmkov(3,3) - gmkov(1,3)*gmkov(3,1))/detgmkov
       gmkon(2,3) = -(gmkov(1,1)*gmkov(2,3) - gmkov(2,1)*gmkov(1,3))/detgmkov
       gmkon(3,1) = (gmkov(2,1)*gmkov(3,2) - gmkov(3,1)*gmkov(2,2))/detgmkov
       gmkon(3,2) = -(gmkov(1,1)*gmkov(3,2) - gmkov(3,1)*gmkov(1,2))/detgmkov
       gmkon(3,3) = (gmkov(1,1)*gmkov(2,2) - gmkov(1,2)*gmkov(2,1))/detgmkov
    else
       write(*,*) 'Warning: element',jelem, 'has negative jacobian gmkov in function g_shell_basis'
    end if

  end subroutine g_shell_basis

  !------------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    January 2019
  !> @brief   Compute curvilinear basis in continuum shell (3D), not in the midsurface
  !> @details
  !------------------------------------------------------------------------------
  subroutine g_shell_basis2(sfcoo,a3,gpsha,gpderi,gpT,gkov,gkon,gmkov,gmkon,detgkov,jelem)

    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    real(rp),       intent(in)             :: &
         sfcoo(pdime,pnodh),                  &  ! Coordinates of midsurface
         a3(pdime,pnodh),                     &  ! Normal to midplane
         gpsha(pnodh),                        &  ! Shape function at integration point (gauss)
         gpderi(pdimh,pnodh),                 &  ! Derivative of shape function at integration point
         gpT                                     ! Coordinate of the gauss point along thickness (xi3)
    integer(ip),    intent(in)             :: &
         jelem
    real(rp),       intent(out)            :: &
         gkov(pdime,pdime),                   &  ! Covariant basis
         gkon(pdime,pdime),                   &  ! Contravariant basis (Inverse of gkov)
         gmkov(pdime,pdime),                  &  ! Metrics of covariant basis
         gmkon(pdime,pdime),                  &  ! Metrics of contravariant basis (Inverse of gmkov)
         detgkov                                 ! Determinant of gkov

    integer(ip)                         :: inode,idime,jdime,kdime
    real(rp)                            :: detgmkov
    !
    ! Covariant basis vector
    !
    gkov = 0.0_rp
    do idime = 1,pdime
       do inode = 1,pnodh
          gkov(idime,1) = gkov(idime,1) + gpderi(1,inode)*sfcoo(idime,inode) + gpT*gpderi(1,inode)*a3(idime,inode)
          gkov(idime,2) = gkov(idime,2) + gpderi(2,inode)*sfcoo(idime,inode) + gpT*gpderi(2,inode)*a3(idime,inode)
          gkov(idime,3) = gkov(idime,3) + gpsha(inode)*a3(idime,inode)
       end do
    end do
    !
    ! Contravariant basis vector
    !
    gkon = 0.0_rp
    call invmtx(gkov, gkon, detgkov, pdime)
    if ( detgkov < 1.0e-15_rp )  then
       print*,'ELEMENT:',leinv_loc(jelem),"VALUE:",detgkov
       call runend("MOD_SLD_CONTSHELL_ELEMENTS: NEGATIVE DETERMINANT IN G_SHELL_BASIS FOR ELEMENT")
    end if
    !
    ! Metrics of Covariant
    !
    do idime = 1,pdime
       do jdime = 1,pdime
          gmkov(idime,jdime) = 0.0_rp
          do kdime = 1,pdime
             gmkov(idime,jdime) = gmkov(idime,jdime) + gkov(kdime,idime)*gkov(kdime,jdime)
          end do
       end do
    end do
    ! Make it symmetric
    gmkov = 0.5_rp*(gmkov + transpose(gmkov))
    !
    ! Metrics of Contravariant
    !
    call invmtx(gmkov, gmkon, detgmkov, pdime)

  end subroutine g_shell_basis2

  !------------------------------------------------------------------------------
  !> @author  Eva Casoni (eva.casoni@bsc.es)
  !> @date    13/09/2016
  !> @brief   Compute B matrix
  !> @details
  !------------------------------------------------------------------------------
  subroutine b_matrix(gpsha, gpderi, ackov, a3ckovder, gpbmat)

    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    real(rp),       intent(in)             :: &
         gpsha(pnodh),                          &  ! Shape function at integration point (gauss)
         gpderi(2,pnodh),                        &  ! Derivative of shape function at integration point
         ackov(pdime,pdime),                     &  ! Curvilienar basis in current configuration
         a3ckovder(pdime,2)                         ! Derivativies of a3 in current coordinates
    real(rp),       intent(out)            :: &
         gpbmat(pdofh,pdofs)                        ! B matrix
    ! ---------------------------------------------------------------------------
    integer(ip)                            :: &
         inode, node_start
    real(rp)                               :: &
         a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z, &
         a31x, a31y, a31z, a32x, a32y, a32z
    !
    ! ---------------------------------------------------------------------------
    !
    ! Defined only for simplicity to construct the matrix
    !
    a1x = ackov(1,1)
    a1y = ackov(2,1)
    a1z = ackov(3,1)

    a2x = ackov(1,2)
    a2y = ackov(2,2)
    a2z = ackov(3,2)

    a3x = ackov(1,3)
    a3y = ackov(2,3)
    a3z = ackov(3,3)

    a31x = a3ckovder(1,1)
    a31y = a3ckovder(2,1)
    a31z = a3ckovder(3,1)

    a32x = a3ckovder(1,2)
    a32y = a3ckovder(2,2)
    a32z = a3ckovder(3,2)
    !
    ! B matrix
    !
    gpbmat(:,:) = 0.0_rp
    node_start = 1

    do inode = 1, pnodh
       gpbmat(1,node_start) = 0.5_rp*gpderi(1,inode)*a1x
       gpbmat(1,node_start+1) = 0.5_rp*gpderi(1,inode)*a1y
       gpbmat(1,node_start+2) = 0.5_rp*gpderi(1,inode)*a1z
       gpbmat(1,node_start+3) = 0.5_rp*gpderi(1,inode)*a1x
       gpbmat(1,node_start+4) = 0.5_rp*gpderi(1,inode)*a1y
       gpbmat(1,node_start+5) = 0.5_rp*gpderi(1,inode)*a1z

       gpbmat(2,node_start) = 0.5_rp*(gpderi(2,inode)*a1x + gpderi(1,inode)*a2x)
       gpbmat(2,node_start+1) = 0.5_rp*(gpderi(2,inode)*a1y + gpderi(1,inode)*a2y)
       gpbmat(2,node_start+2) = 0.5_rp*(gpderi(2,inode)*a1z + gpderi(1,inode)*a2z)
       gpbmat(2,node_start+3) = 0.5_rp*(gpderi(2,inode)*a1x + gpderi(1,inode)*a2x)
       gpbmat(2,node_start+4) = 0.5_rp*(gpderi(2,inode)*a1y + gpderi(1,inode)*a2y)
       gpbmat(2,node_start+5) = 0.5_rp*(gpderi(2,inode)*a1z + gpderi(1,inode)*a2z)

       gpbmat(3,node_start) = 0.5_rp*(gpderi(1,inode)*a3x + gpsha(inode)*a1x)
       gpbmat(3,node_start+1) = 0.5_rp*(gpderi(1,inode)*a3y + gpsha(inode)*a1y)
       gpbmat(3,node_start+2) = 0.5_rp*(gpderi(1,inode)*a3z + gpsha(inode)*a1z)
       gpbmat(3,node_start+3) = 0.5_rp*(gpderi(1,inode)*a3x - gpsha(inode)*a1x)
       gpbmat(3,node_start+4) = 0.5_rp*(gpderi(1,inode)*a3y - gpsha(inode)*a1y)
       gpbmat(3,node_start+5) = 0.5_rp*(gpderi(1,inode)*a3z - gpsha(inode)*a1z)

       gpbmat(4,node_start) = 0.5_rp*gpderi(2,inode)*a2x
       gpbmat(4,node_start+1) = 0.5_rp*gpderi(2,inode)*a2y
       gpbmat(4,node_start+2) = 0.5_rp*gpderi(2,inode)*a2z
       gpbmat(4,node_start+3) = 0.5_rp*gpderi(2,inode)*a2x
       gpbmat(4,node_start+4) = 0.5_rp*gpderi(2,inode)*a2y
       gpbmat(4,node_start+5) = 0.5_rp*gpderi(2,inode)*a2z

       gpbmat(5,node_start) = 0.5_rp*(gpderi(2,inode)*a3x + gpsha(inode)*a2x)
       gpbmat(5,node_start+1) = 0.5_rp*(gpderi(2,inode)*a3y + gpsha(inode)*a2y)
       gpbmat(5,node_start+2) = 0.5_rp*(gpderi(2,inode)*a3z + gpsha(inode)*a2z)
       gpbmat(5,node_start+3) = 0.5_rp*(gpderi(2,inode)*a3x - gpsha(inode)*a2x)
       gpbmat(5,node_start+4) = 0.5_rp*(gpderi(2,inode)*a3y - gpsha(inode)*a2y)
       gpbmat(5,node_start+5) = 0.5_rp*(gpderi(2,inode)*a3z - gpsha(inode)*a2z)

       gpbmat(6,node_start) = 0.5_rp*gpsha(inode)*a3x
       gpbmat(6,node_start+1) = 0.5_rp*gpsha(inode)*a3y
       gpbmat(6,node_start+2) = 0.5_rp*gpsha(inode)*a3z
       gpbmat(6,node_start+3) = -0.5_rp*gpsha(inode)*a3x
       gpbmat(6,node_start+4) = -0.5_rp*gpsha(inode)*a3y
       gpbmat(6,node_start+5) = -0.5_rp*gpsha(inode)*a3z

       gpbmat(7,node_start) = 0.5_rp*(gpderi(1,inode)*a31x + gpderi(1,inode)*a1x)
       gpbmat(7,node_start+1) = 0.5_rp*(gpderi(1,inode)*a31y + gpderi(1,inode)*a1y)
       gpbmat(7,node_start+2) = 0.5_rp*(gpderi(1,inode)*a31z + gpderi(1,inode)*a1z)
       gpbmat(7,node_start+3) = 0.5_rp*(gpderi(1,inode)*a31x - gpderi(1,inode)*a1x)
       gpbmat(7,node_start+4) = 0.5_rp*(gpderi(1,inode)*a31y - gpderi(1,inode)*a1y)
       gpbmat(7,node_start+5) = 0.5_rp*(gpderi(1,inode)*a31z - gpderi(1,inode)*a1z)

       gpbmat(8,node_start) = 0.5_rp*(gpderi(1,inode)*a32x + gpderi(2,inode)*a1x + gpderi(2,inode)*a31x + gpderi(1,inode)*a2x)
       gpbmat(8,node_start+1) = 0.5_rp*(gpderi(1,inode)*a32y + gpderi(2,inode)*a1y + gpderi(2,inode)*a31y + gpderi(1,inode)*a2y)
       gpbmat(8,node_start+2) = 0.5_rp*(gpderi(1,inode)*a32z + gpderi(2,inode)*a1z + gpderi(2,inode)*a31z + gpderi(1,inode)*a2z)
       gpbmat(8,node_start+3) = 0.5_rp*(gpderi(1,inode)*a32x - gpderi(2,inode)*a1x + gpderi(2,inode)*a31x - gpderi(1,inode)*a2x)
       gpbmat(8,node_start+4) = 0.5_rp*(gpderi(1,inode)*a32y - gpderi(2,inode)*a1y + gpderi(2,inode)*a31y - gpderi(1,inode)*a2y)
       gpbmat(8,node_start+5) = 0.5_rp*(gpderi(1,inode)*a32z - gpderi(2,inode)*a1z + gpderi(2,inode)*a31z - gpderi(1,inode)*a2z)

       gpbmat(9,node_start) = 0.5_rp*(gpsha(inode)*a31x + gpderi(1,inode)*a3x)
       gpbmat(9,node_start+1) = 0.5_rp*(gpsha(inode)*a31y + gpderi(1,inode)*a3y)
       gpbmat(9,node_start+2) = 0.5_rp*(gpsha(inode)*a31z + gpderi(1,inode)*a3z)
       gpbmat(9,node_start+3) = -0.5_rp*(gpsha(inode)*a31x + gpderi(1,inode)*a3x)
       gpbmat(9,node_start+4) = -0.5_rp*(gpsha(inode)*a31y + gpderi(1,inode)*a3y)
       gpbmat(9,node_start+5) = -0.5_rp*(gpsha(inode)*a31z + gpderi(1,inode)*a3z)

       gpbmat(10,node_start) = 0.5_rp*(gpderi(2,inode)*a32x + gpderi(2,inode)*a2x)
       gpbmat(10,node_start+1) = 0.5_rp*(gpderi(2,inode)*a32y + gpderi(2,inode)*a2y)
       gpbmat(10,node_start+2) = 0.5_rp*(gpderi(2,inode)*a32z + gpderi(2,inode)*a2z)
       gpbmat(10,node_start+3) = 0.5_rp*(gpderi(2,inode)*a32x - gpderi(2,inode)*a2x)
       gpbmat(10,node_start+4) = 0.5_rp*(gpderi(2,inode)*a32y - gpderi(2,inode)*a2y)
       gpbmat(10,node_start+5) = 0.5_rp*(gpderi(2,inode)*a32z - gpderi(2,inode)*a2z)

       gpbmat(11,node_start) = 0.5_rp*(gpsha(inode)*a32x + gpderi(2,inode)*a3x)
       gpbmat(11,node_start+1) = 0.5_rp*(gpsha(inode)*a32y + gpderi(2,inode)*a3y)
       gpbmat(11,node_start+2) = 0.5_rp*(gpsha(inode)*a32z + gpderi(2,inode)*a3z)
       gpbmat(11,node_start+3) = -0.5_rp*(gpsha(inode)*a32x + gpderi(2,inode)*a3x)
       gpbmat(11,node_start+4) = -0.5_rp*(gpsha(inode)*a32y + gpderi(2,inode)*a3y)
       gpbmat(11,node_start+5) = -0.5_rp*(gpsha(inode)*a32z + gpderi(2,inode)*a3z)

       gpbmat(12,node_start) =  0.0_rp
       gpbmat(12,node_start+1) =  0.0_rp
       gpbmat(12,node_start+2) =  0.0_rp
       gpbmat(12,node_start+3) =  0.0_rp
       gpbmat(12,node_start+4) =  0.0_rp
       gpbmat(12,node_start+5) =  0.0_rp

       node_start = node_start + 6
    end do

  end subroutine b_matrix

  !------------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    January 2019
  !> @brief   Material law manager
  !> @details Material law manager
  !------------------------------------------------------------------------------

  subroutine material_law(ielem,gmrkov,gmckov,gmrkon,gmckon,grkon,detrg,detcg, &
       stiff,stress)

    implicit none

    integer(ip), intent(in)  :: ielem               !< Element number
    real(rp),    intent(in)  :: gmrkov(pdime,pdime) !< Metrics of covariant basis in reference configuration
    real(rp),    intent(in)  :: gmckov(pdime,pdime) !< Metrics of covariant basis in current configuration
    real(rp),    intent(in)  :: gmrkon(pdime,pdime) !< Metrics of contravariant basis in reference configuration
    real(rp),    intent(in)  :: gmckon(pdime,pdime) !< Metrics of contravariant basis in current configuration
    real(rp),    intent(in)  :: grkon(pdime,pdime)  !< Contravariant basis in reference configuration
    real(rp),    intent(in)  :: detrg               !< 
    real(rp),    intent(in)  :: detcg               !< 
    real(rp),    intent(out) :: stress(6)           !<
    real(rp),    intent(out) :: stiff(6,6)          !<

    integer(ip)              :: i, j
    integer(ip)              :: pmate, pmlaw
    real(rp)                 :: gpgre(pdime,pdime)
    real(rp)                 :: strain(6)

    !
    ! Material code and law
    !
    pmate = lmate_sld(ielem)
    pmlaw = lawco_sld(pmate)
    !
    ! Compute Green-Lagrange strain tensor
    !
    gpgre(:,:) = 0.0_rp
    do i = 1, pdime
       do j = 1, pdime
          gpgre(i,j) = 0.5_rp*(gmckov(i,j) - gmrkov(i,j))
       end do
    end do
    strain(1) = gpgre(1,1)
    strain(2) = gpgre(1,2)*2.0_rp
    strain(3) = gpgre(1,3)*2.0_rp
    strain(4) = gpgre(2,2)
    strain(5) = gpgre(2,3)*2.0_rp
    strain(6) = gpgre(3,3)
    !
    ! Material law selector
    !
    if ( pmlaw == 1_ip ) then
       !
       ! Isotropic linear elastic model
       !
       call cshell_stress_model_100(pmate,gmrkon,strain,stress,stiff)

    else if ( pmlaw == 2_ip) then
       !
       ! Orthotropic linear elastic model
       !
       call cshell_stress_model_151(ielem,pmate,grkon,strain,stress,stiff)

    else if ( pmlaw == 3_ip ) then
       !
       ! Orthotropic damage Maimi (2018)
       !
       call runend('MATERIAL_LAW_NOT_IMPLEMENTED')
       
    else if ( pmlaw == 4_ip ) then
       !
       ! Hyperelastic Neo-Hookean
       !                                              G^{ij},g^{ij},     S,    C
       call cshell_stress_model_101(pmate,detrg,detcg,gmrkon,gmckon,stress,stiff)
   
    end if

  end subroutine material_law

  !------------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    January 2019
  !> @brief   Isotropic linear elastic model
  !> @details
  !------------------------------------------------------------------------------

  subroutine cshell_stress_model_100(pmate,gmrkon,strain,stress,stiff)

    implicit none

    integer(ip), intent(in)  :: pmate
    real(rp),    intent(in)  :: gmrkon(pdime,pdime) !< Inverse of gmrkov
    real(rp),    intent(in)  :: strain(6)           !<
    real(rp),    intent(out) :: stress(6)           !< 2PK stress tensor
    real(rp),    intent(out) :: stiff(6,6)          !<

    integer(ip)              :: i,j,k,l
    real(rp)                 :: E, nu, lame1, lame2
    real(rp)                 :: C(pdime,pdime,pdime,pdime)
    real(rp)                 :: vodds_glo(6,6)
    !
    ! Material properties
    !
    E  = parco_sld(1,pmate)
    nu = parco_sld(2,pmate)
    lame2 = E/(2.0_rp*(1.0_rp+nu))
    lame1 = nu*E/((1.0_rp+nu)*(1.0_rp-2.0_rp*nu))
    !
    ! Compute 4th order stiffness tensor
    !
    C(:,:,:,:) = 0.0_rp
    do i = 1,pdime
       do j = 1,pdime
          do k = 1,pdime
             do l = 1,pdime
                C(i,j,k,l)= lame1 * gmrkon(i,j)*gmrkon(k,l)+ &
                     lame2 * ( gmrkon(i,k)*gmrkon(j,l) + gmrkon(i,l)*gmrkon(k,j) )
             end do
          end do
       end do
    end do
    !
    ! Tensor to Voigt notation
    !
    call SM_tensor_to_voigt_fourth(pdime, 0_ip, C(:,:,:,:), vodds_glo(:,:))
    !
    ! Stiffness tensor (Reordering to Shells Voigt notation)
    !
    ! Original order      -->  New order (Shells)
    ! [11 22 33 23 13 12] -->  [11 12 13 22 23 33]
    !
    call s9Creordering(vodds_glo(:,:),stiff(:,:))
    !
    ! Stress (2PK)
    !
    call SM_stress_tensor(0_ip, strain, stiff, stress)

  end subroutine cshell_stress_model_100

  !------------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    January 2019
  !> @brief   Orthotropic linear elastic model
  !> @details
  !------------------------------------------------------------------------------

  subroutine cshell_stress_model_151(ielem,pmate,grkon,strain,stress,stiff)

    use def_solidz, only : stiff0_sld,rmate_sld

    implicit none

    integer(ip), intent(in)  :: ielem
    integer(ip), intent(in)  :: pmate
    real(rp),    intent(in)  :: grkon(pdime,pdime)  !<
    real(rp),    intent(in)  :: strain(6)           !<
    real(rp),    intent(out) :: stress(6)           !< 2PK stress tensor
    real(rp),    intent(out) :: stiff(6,6)          !<

    integer(ip)              :: idime,igaus
    real(rp)                 :: rotmat(pdime,pdime)
    real(rp)                 :: C(pdime,pdime,pdime,pdime)
    real(rp)                 :: vodds_glo(6,6)
    real(rp)                 :: vodds_loc(6,6)
    real(rp)                 :: vodds_reo(6,6)
    real(rp)                 :: F(pdime,pdime)
    !
    ! Stiffness matrix
    !
    vodds_loc(:,:) = stiff0_sld(:,:,pmate)
    !
    ! Rotate from material to global coordinate system (cartesian CSYS)
    !
    F = 0.0_rp
    do idime = 1,pdime
       F(idime,idime) = 1.0_rp
    end do
    call SM_rotate_basis_creation(ielem, F, rotmat)
    call SM_rotate_voigt_fourth(2_ip,rotmat,vodds_glo,vodds_loc)
    !
    ! Reordening the Cij components
    !
    ! Original order      -->  New order (Shells)
    ! [11 22 33 23 13 12] -->  [11 12 13 22 23 33]
    !
    call s9Creordering(vodds_glo(:,:),vodds_reo(:,:))
    !
    ! Compute 4th order stiffness tensor (Reads Voigt Shell ordering)
    !
    call s81CurvTransf81(grkon(:,:),vodds_reo(:,:),C(:,:,:,:))
    !
    ! Tensor to Voig notation
    !
    call SM_tensor_to_voigt_fourth(pdime, 0_ip, C(:,:,:,:), vodds_glo(:,:))
    !
    ! Sriffness tensor (Reordering to Shells Voigt notation)
    !
    ! Original order      -->  New order (Shells)
    ! [11 22 33 23 13 12] -->  [11 12 13 22 23 33]
    !
    call s9Creordering(vodds_glo(:,:),stiff(:,:))
    !
    ! Stress (2PK)
    !
    call SM_stress_tensor(0_ip, strain, stiff, stress)
    !
    ! Store variables
    ! GGU: Provisional we need to fill this variable for each gauss point 
    do igaus = 1,8
       rmate_sld(ielem) % a(:,:,igaus) = rotmat(:,:)
    end do
    
  end subroutine cshell_stress_model_151
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillamet
  !> @date    2022-03-07
  !> @brief   Hyperelastic Neo-Hookean model for continnuu-shells
  !> @details 
  !>          Ref: Kiendl, J. et al. Isogeometric Kirchhoff–Love
  !>               shell formulations for general hyperelastic materials
  !>  
  !-----------------------------------------------------------------------
  
  subroutine cshell_stress_model_101(pmate,detrg,detcg,gmrkon,gmckon,stress,stiff)

    implicit none

    integer(ip), intent(in)  :: pmate       !< Material code
    real(rp),    intent(in)  :: detrg       !< det(G_{ij})
    real(rp),    intent(in)  :: detcg       !< det(g_{ij})
    real(rp),    intent(in)  :: gmrkon(3,3) !< G^{ij}
    real(rp),    intent(in)  :: gmckon(3,3) !< g^{ij}
    real(rp),    intent(out) :: stress(6)   !< 2PK stress tensor
    real(rp),    intent(out) :: stiff(6,6)  !< Elasticity tensor

    integer(ip)              :: i,j,k,l
    real(rp)                 :: E, nu, lambda, mu
    real(rp)                 :: J_0
    real(rp)                 :: C_SE(3,3,3,3),C_SEvo(6,6)
    real(rp)                 :: S(3,3)
    !
    ! Material properties
    !
    E  = parco_sld(1,pmate)
    nu = parco_sld(2,pmate)
    lambda = nu*E/((1.0_rp+nu)*(1.0_rp-2.0_rp*nu))
    mu     = E/(2.0_rp*(1.0_rp+nu))
    J_0    = sqrt(detcg/detrg)  ! In-plane Jacobian determinant
    !
    ! 4th order elasticity tensor (tangent moduli)
    !
    C_SE(:,:,:,:) = 0.0_rp
    do i = 1,3
       do j = 1,3
          do k = 1,3
             do l = 1,3
                C_SE(i,j,k,l) = mu*((1.0_rp/J_0)**2)*( 2.0_rp*gmckon(i,j)*gmckon(k,l) + &
                     gmckon(i,k)*gmckon(j,l) + gmckon(i,l)*gmckon(j,k) )
             end do
          end do
       end do
    end do
    !
    ! Tensor to Voigt notation
    !
    call SM_tensor_to_voigt_fourth(3_ip, 0_ip, C_SE(:,:,:,:), C_SEvo(:,:))
    !
    ! Stiffness tensor (Reordering to Shells Voigt notation)
    !
    ! Original order      -->  New order (Shells)
    ! [11 22 33 23 13 12] -->  [11 12 13 22 23 33]
    !
    call s9Creordering(C_SEvo(:,:),stiff(:,:))
    !
    ! Stress (2PK)
    !
    S(:,:) = mu*( gmrkon(:,:) - ((1.0_rp/J_0)**2)*gmckon(:,:) )
    !
    ! Tensor to Voigt and reordering
    !
    ! Original order      -->  New order (Shells)
    ! [11 22 33 23 13 12] -->  [11 12 13 22 23 33]
    !
    stress(1) = S(1,1)
    stress(2) = S(1,2)
    stress(3) = S(1,3)
    stress(4) = S(2,2)
    stress(5) = S(2,3)
    stress(6) = S(3,3)

  end subroutine cshell_stress_model_101

  !------------------------------------------------------------------------------
  !> @author  Eva Casoni (eva.casoni@bsc.es)
  !> @date    13/09/2016
  !> @brief   Integrate material law along thickness of the shell
  !> @details
  !------------------------------------------------------------------------------
  subroutine integrate_material_law(vodds_aux,stress,wT,T,gpddsVoig,gpstrVoig)

    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    real(rp),       intent(in)            :: &
         wT,                                 & ! weight for the gauss point of the thickness
         T,                                  & ! gauss point of the thickness
         vodds_aux(6,6),                     &
         stress(6)                             ! input stress in voigt notation

    real(rp),       intent(inout)         :: &
         gpstrVoig(12),                      & ! stress resultant
         gpddsVoig(12,12)

    ! ---------------------------------------------------------------------------
    integer(ip)                           :: &
         i,j,k,l
    !
    !---------------------------------------------------------------------------------
    !
    ! Integrate along thickness
    !
    do i = 1,6
       k = i+6_ip
       gpstrVoig(i) = gpstrVoig(i) + stress(i)*wT
       gpstrVoig(k) = gpstrVoig(k) + stress(i)*wT*T
       do j = 1,6
          l = j+6_ip
          gpddsVoig(i,j) = gpddsVoig(i,j) + vodds_aux(i,j)*wT
          gpddsVoig(k,j) = gpddsVoig(k,j) + vodds_aux(i,j)*wT*T
          gpddsVoig(k,l) = gpddsVoig(k,l) + vodds_aux(i,j)*wT*T*T
       end do
    end do
    !
    ! Force symmetry
    gpddsVoig = 0.5_rp*(gpddsVoig + transpose(gpddsVoig))

  end subroutine integrate_material_law

  !------------------------------------------------------------------------------
  !> @author  Eva Casoni (eva.casoni@bsc.es)
  !> @date    17/09/2016
  !> @brief   Compute linear part of k_uu : int(B^t*C*B)
  !> @details
  !------------------------------------------------------------------------------
  subroutine computeStiffMatrix(Bmat,ddsdde,w1,w2,hnormi,stiff)

    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    real(rp),       intent(in)             :: &
         Bmat(12,24),                            &    ! Bla bla
         ddsdde(12,12),                          &    ! Tangent moduli
         w1,w2,hnormi
    real(rp),       intent(out)            :: &
         stiff(24,24)
    real(rp)                               :: &
         aux(12,24),                          &
         BmatT(24,12),                        &    ! Transpose of Bmat
         factor                                       ! w1*w2*hnorm

    !
    ! Initialization
    !
    stiff = 0.0_rp
    factor = w1*w2*hnormi
    !
    ! B^t*(C*B)
    BmatT = transpose(Bmat)
    aux   = matmul(ddsdde,Bmat)
    stiff = matmul(BmatT,aux)*factor

  end subroutine computeStiffMatrix

  !------------------------------------------------------------------------------
  !> @author  Eva Casoni (eva.casoni@bsc.es)
  !> @date    17/09/2016
  !> @brief   Compute internal force vector: int(B*stress)
  !> @details
  !------------------------------------------------------------------------------
  subroutine  internalForceVector(Bmat,stress,w1,w2,hnormi,fint)

    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    real(rp),       intent(in)             :: &
         Bmat(12,24),                            &    ! Bla bla
         stress(12),                             &    ! Tangent moduli
         w1,w2,hnormi
    real(rp),       intent(out)            :: &
         fint(24)                                     ! Internal force
    integer(ip)                            :: &
         i, j
    real(rp)                               :: &
         factor

    fint = 0.0_rp
    factor = w1*w2*hnormi

    do j = 1,24
       do i= 1,12
          fint(j) = fint(j) + Bmat(i,j)*stress(i)*factor
       end do
    end do

  end subroutine internalForceVector

  !------------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    January 2019
  !> @brief   Compute inertial contributions
  !> @details
  !>
  !------------------------------------------------------------------------------

  subroutine inertialContribution(pdime,pnode,ielem,elmas,elfine,elstif)

    use def_master, only : dtime
    use def_solidz, only : tifac_sld

    implicit none

    integer(ip), intent(in)    :: pdime               !< Dimensions
    integer(ip), intent(in)    :: pnode               !< Number of element nodes
    integer(ip), intent(in)    :: ielem               !< Element number
    real(rp),    intent(in)    :: elmas(pnode,pnode)  !< Consistent mass matrix
    real(rp),    intent(inout) :: elfine(pnode*pdime) !< Inertial force vector
    real(rp),    intent(inout) :: elstif(pdofs,pdofs) !< Inertial term

    integer(ip)                :: idime,inode,ipoin,jnode,ievat,jevat
    real(rp)                   :: betanewmark

    !
    ! Newmark parameters
    !
    betanewmark = 0.0_rp
    if( tifac_sld(1) > 0.0_rp) betanewmark = 1.0_rp / (tifac_sld(1)*dtime*dtime)
    !
    ! Inertial force vector
    !
    do idime = 1, pdime
       do inode = 1, pnode
          ievat = (inode-1)*pdime + idime
          do jnode = 1, pnode
             ipoin = lnods(jnode,ielem)
             !
             ! Inertial forces
             !   f_kb += M_ba * accel_ak
             elfine(ievat) = elfine(ievat) + elmas(inode,jnode)*accel_sld(idime,ipoin,ITER_K)
          end do
       end do
    end do
    !
    ! Inertial term for tangent stiffness matrix
    !
    do idime = 1, pdime
       do inode = 1, pnode
          ievat = (inode-1)*pdime + idime
          do jnode = 1, pnode
             jevat = (jnode-1)*pdime + idime
             !
             ! Inertial term contribution
             !   K_iakb += 1/(beta*dtime^2)*M
             elstif(ievat,jevat) = elstif(ievat,jevat) + elmas(inode,jnode)*betanewmark
          end do
       end do
    end do

  end subroutine inertialContribution

  !------------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    January 2019
  !> @brief   Compute mass matrices
  !> @details De moment s'estant calculant amb els shape functions d'un element solid
  !-----------------------------------------------------------------------------

  subroutine compute_mass_matrix(pmate,pnode,pgausMass,gpvol,gpsha,elmas,elmuu)

    implicit none
    
    integer(ip), intent(in)    :: pmate                  !< Material code
    integer(ip), intent(in)    :: pnode                  !< Number of element nodes
    integer(ip), intent(in)    :: pgausMass              !< Number of gauss points
    real(rp),    intent(in)    :: gpvol(pgausMass)       !< PROVISIONAL
    real(rp),    intent(in)    :: gpsha(pnode,pgausMass) !< PROVISIONAL
    real(rp),    intent(inout) :: elmas(pnode,pnode)     !< Consistent mass matrix
    real(rp),    intent(inout) :: elmuu(pnode)           !< Lumped mass matrix

    integer(ip)                :: igaus,inode,jnode
    real(rp)                   :: volux

    do igaus = 1,pgausMass  !> GGU Gauss de fora!!
       volux = gpvol(igaus)
       do inode = 1,pnode
          do jnode = 1,pnode
             !
             ! Consistent mass matrix
             !   M_ab = int_\Omega /rho * Na * Nb d\Omega
             elmas(inode,jnode) = elmas(inode,jnode) + densi_sld(1,pmate)*volux*gpsha(inode,igaus)*gpsha(jnode,igaus)
             !
             ! Lumped mass matrix
             !
             elmuu(inode) = elmuu(inode) + densi_sld(1,pmate)*volux*gpsha(inode,igaus)*gpsha(jnode,igaus)
          end do
       end do
    end do

  end subroutine compute_mass_matrix

  !------------------------------------------------------------------------------------
  !
  ! Subroutines for EAS method
  !
  !-------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !> @author  Eva Casoni (eva.casoni@bsc.es)
  !> @date    17/09/2016
  !> @brief   Compute matrices for the stiffness matrix of the system with the enhanced strains
  !> @details
  !------------------------------------------------------------------------------
  subroutine compute_matrices(gpddsVoig,Bmat,Mglo,w1,w2,hnormi,DEAS,LEAS)

    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    real(rp),       intent(in)            :: &
         gpddsVoig(12,12),                      &  ! tangent moduli
         Bmat(12,24),                           &  ! B matrix
         Mglo(12,22),                           &  ! M EAS matrix
         w1,w2,hnormi                              ! weights of gauss points of the plane

    real(rp),       intent(inout)         :: &
         DEAS(22,22),                           &  !
         LEAS(22,24)                               !

    integer(ip)                           :: &
         i, j
    real(rp)                              :: &
         MgloT(22,12),                          &  ! Transpose of Mglo
         D(22,22),                              &  ! Auxiliar matrix: D=M^t*C*M
         auxD(12,22),                           &
         L(22,24),                              &  ! Auxiliar matrix: L=M^t*C*B
         auxL(12,24),                           &
         facti

    !
    ! Integration factor
    !
    facti = w1*w2*hnormi !this can go outside the function
    !
    ! Construct k_alphaalpha , k_ualpha, k_alphau, k_uu
    !
    !      k_alphaalpha = int(D)      D = (M^T)*C*M
    !      k_alphau = int(L)          L = (M^T)*C*B
    !      k_ualpha = int(L^T)
    !      k_uu = int((B^T)*C*B + (dB/du)^T*Stress)
    !
    MgloT = transpose(Mglo)
    auxD  = matmul(gpddsVoig,Mglo)
    auxL  = matmul(gpddsVoig,Bmat)
    ! D = M^t * C * M
    ! L = M^t * C * B
    D = matmul(MgloT,auxD)
    L = matmul(MgloT,auxL)

    do i = 1,22
       do j = 1,22
          DEAS(i,j) = DEAS(i,j) + D(i,j)*facti
       end do
    end do

    do i = 1,22
       do j = 1,24
          LEAS(i,j) = LEAS(i,j) + L(i,j)*facti
       end do
    end do

  end subroutine compute_matrices

  !------------------------------------------------------------------------------
  !> @author  Eva Casoni (eva.casoni@bsc.es)
  !> @date    07/04/2017
  !> @brief   Compute external force (f_EAS)
  !> @details
  !------------------------------------------------------------------------------

  subroutine compute_EAS_force(gpstrVoig,Mglo,w1,w2,hnormi,fEAS)

    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    real(rp),       intent(in)            :: &
         gpstrVoig(12),                      &  ! 2PK
         Mglo(12,22),                        &  ! M EAS matrix
         w1,w2,hnormi                           ! weights of gauss points of the plane
    real(rp),       intent(inout)         :: &
         fEAS(22)                               ! external force due to enhanced strains
    integer(ip)                           :: &
         i, j
    real(rp)                              :: &
         MgloT(22,12),facti
    !
    ! Integration factor
    !
    facti = w1*w2*hnormi !this can go outside the function
    !
    ! Global M tranpose  matrix
    !
    MgloT = transpose(Mglo)
    !
    ! f_EAS = int(Mt*stress)
    !
    do i = 1,22
       do j = 1, 12
          fEAS(i) = fEAS(i) + MgloT(i,j)*gpstrVoig(j)*facti
       end do
    end do

  end subroutine compute_EAS_force

  !------------------------------------------------------------------------------
  !> @author  Eva Casoni (eva.casoni@bsc.es)
  !> @date    10/10/2016
  !> @brief   Compute M_xi matrix for enhanced strains (EAS method)
  !> @details Math:
  !>            M = [Mp Mq]^T = [M_7  0_7  0_4      0_4
  !>                             0_7  M_7  M_4^q33  M_4]
  !>  such that:
  !>             E_enhanced = M_xi*alpha
  !>             alpha: enhancing strain vector
  !>
  !------------------------------------------------------------------------------
  subroutine EASM_local(r,s,Mloc)

    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    real(rp),       intent(in)             :: &
         r,s

    integer(ip)                            :: &  ! Local coordinates xi1,xi2
         rr,rs,qr,ss,qs,qt,                   &
         mrr,mrs,mqr,mss,mqs,mqt,             &
         indexcol
    real(rp),       intent(out)            :: &
         Mloc(12,22)

    rr = 1_ip
    rs = 2_ip
    qr = 3_ip
    ss = 4_ip
    qs = 5_ip
    qt = 6_ip

    mrr = 7_ip
    mrs = 8_ip
    mqr = 9_ip
    mss = 10_ip
    mqs = 11_ip
    mqt = 12_ip
    !
    ! Initialization
    !
    Mloc = 0.0_rp
    indexcol = 1_ip
    !
    ! Construct M_7
    !
    Mloc(rr,indexcol)   = r
    Mloc(ss,indexcol+1) = s
    Mloc(rs,indexcol+2) = r
    Mloc(rs,indexcol+3) = s
    Mloc(rr,indexcol+4) = r*s
    Mloc(ss,indexcol+5) = r*s
    Mloc(rs,indexcol+6) = r*s
    !
    ! Add M_7 to second block row
    !
    indexcol = indexcol + 7

    Mloc(mrr,indexcol)   = r
    Mloc(mss,indexcol+1) = s
    Mloc(mrs,indexcol+2) = r
    Mloc(mrs,indexcol+3) = s
    Mloc(mrr,indexcol+4) = r*s
    Mloc(mss,indexcol+5) = r*s
    Mloc(mrs,indexcol+6) = r*s
    !
    ! Construct M_4^q33
    !
    indexcol = indexcol+7

    Mloc(mqt,indexcol)   = 1.0_rp
    Mloc(mqt,indexcol+1) = r
    Mloc(mqt,indexcol+2) = s
    Mloc(mqt,indexcol+3) = r*s
    !
    ! Construct M_4
    !
    indexcol = indexcol + 4

    Mloc(mqr,indexcol)   = r
    Mloc(mqr,indexcol+1) = r*s
    Mloc(mqs,indexcol+2) = s
    Mloc(mqs,indexcol+3) = r*s

  end subroutine EASM_local

  !------------------------------------------------------------------------------
  !> @author  Eva Casoni (eva.casoni@bsc.es)
  !> @date    10/10/2016
  !> @brief   Compute M  matrix for enhanced strains in global coordinates  (EAS method)
  !> @details Math:
  !>            M = (detJ0/detJ)*T0^(-t)*M_xi
  !>
  !>  where:
  !>             detJ0 = det(akov) in midpoint shell
  !>             detJ = det(akov) in shell in reference configuration
  !>             T0:
  !>
  !------------------------------------------------------------------------------

  subroutine  EASM_transform(Mloc, Mglo, arkov, a0rkov, detra, det0r, T)

    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    real(rp),       intent(in)             :: &
         Mloc(12,22),                         & ! M_xi
         arkov(pdime,pdime),                  & ! curvilinear basis of the midplane shell in reference coordinates
         a0rkov(pdime,pdime),                 & ! curvilinear basis of the midpoint of midplane shell in reference coordinates
         detra,                               & ! detJ: determinant of arkov
         det0r                                     ! detJ0: determinant of a0rkov
    real(rp),       intent(out)            :: &
         Mglo(12,22)                               ! M in global coordinates
    real(rp)                               :: &
         T(12,12),                            & ! Transformation matrix
         factor0,                             & ! detJ0/detJ
         t11,t22,t33,t12,t13,t21,t23,t31,t32
    integer(ip)                            :: &
         i

    !
    ! Initialization
    !
    T = 0.0_rp
    Mglo = 0.0_rp

    t11 = 0.0_rp
    t12 = 0.0_rp
    t13 = 0.0_rp

    t21 = 0.0_rp
    t22 = 0.0_rp
    t23 = 0.0_rp

    t31 = 0.0_rp
    t32 = 0.0_rp
    t33 = 1.0_rp
    !
    ! Construct transformation matrix
    !
    do i=1,3
       t11 = t11 + arkov(i,1)*a0rkov(i,1)
       t12 = t12 + arkov(i,1)*a0rkov(i,2)
       t13 = t13 + arkov(i,1)*a0rkov(i,3)

       t21 = t21 + arkov(i,2)*a0rkov(i,1)
       t22 = t22 + arkov(i,2)*a0rkov(i,2)
       !t23 = t23 + arkov(i,2)*a0rkov(i,3)

       !t31 = t31 + arkov(i,3)*a0rkov(i,1)
       !t32 = t32 + arkov(i,3)*a0rkov(i,2)
       !t33 = t33 + arkov(i,3)*a0rkov(i,3)
    end do

    factor0 = 0.0_rp
    if (abs(detra) > 1.0e-15_rp) then
       factor0 = det0r/detra
    else
       write(*,*) '--- detra = 0 in cmoputing EASTmatrix'
    end if

    T(1,1) = factor0*t11*t11
    T(2,1) = factor0*2*t11*t21
    T(3,1) = factor0*2*t11*t31

    T(4,1) = factor0*t21*t21
    T(5,1) = factor0*2*t21*t31
    T(6,1) = factor0*t31*t31

    T(1,2) = factor0*t11*t12
    T(2,2) = factor0*(t11*t22 + t21*t12)
    T(3,2) = factor0*(t11*t32 + t31*t12)

    T(4,2) = factor0*t21*t22
    T(5,2) = factor0*(t21*t32 + t31*t22)
    T(6,2) = factor0*t31*t32

    T(1,3) = factor0*t11*t13
    T(2,3) = factor0*(t11*t23 + t21*t13)
    T(3,3) = factor0*(t11*t33 + t31*t13)

    T(4,3) = factor0*t21*t23
    T(5,3) = factor0*(t21*t33 + t31*t23)
    T(6,3) = factor0*t31*t33

    T(1,4) = factor0*t12*t12
    T(2,4) = factor0*2*t12*t22
    T(3,4) = factor0*2*t12*t32

    T(4,4) = factor0*t22*t22
    T(5,4) = factor0*2*t22*t32
    T(6,4) = factor0*t32*t32

    T(1,5) = factor0*t12*t13
    T(2,5) = factor0*(t12*t23 + t22*t13)
    T(3,5) = factor0*(t12*t33 + t32*t13)

    T(4,5) = factor0*t22*t23
    T(5,5) = factor0*(t22*t33 + t32*t23)
    T(6,5) = factor0*t32*t33

    T(1,6) = factor0*t13*t13
    T(2,6) = factor0*2*t13*t23
    T(3,6) = factor0*2*t13*t33

    T(4,6) = factor0*t23*t23
    T(5,6) = factor0*2*t23*t33
    T(6,6) = factor0*t33*t33

    T(7,7) = factor0 *t11*t11
    T(8,7) = factor0*2*t11*t21
    T(9,7) = factor0*2*t11*t31

    T(10,7) = factor0*t21*t21
    T(11,7) = factor0*2*t21*t31
    T(12,7) = factor0*t31*t31

    T(7,8) = factor0*t11*t12
    T(8,8) = factor0*(t11*t22 + t21*t12)
    T(9,8) = factor0*(t11*t32 + t31*t12)

    T(10,8) = factor0*t21*t22
    T(11,8) = factor0*(t21*t32 + t31*t22)
    T(12,8) = factor0*t31*t32

    T(7,9) = factor0*t11*t13
    T(8,9) = factor0*(t11*t23 + t21*t13)
    T(9,9) = factor0*(t11*t33 + t31*t13)

    T(10,9) = factor0*t21*t23
    T(11,9) = factor0*(t21*t33 + t31*t23)
    T(12,9) = factor0*t31*t33

    T(7,10) = factor0*t12*t12
    T(8,10) = factor0*2*t12*t22
    T(9,10) = factor0*2*t12*t32

    T(10,10) = factor0*t22*t22
    T(11,10) = factor0*2*t22*t32
    T(12,10) = factor0*t32*t32

    T(7,11) = factor0*t12*t13
    T(8,11) = factor0*(t12*t23 + t22*t13)
    T(9,11) = factor0*(t12*t33 + t32*t13)

    T(10,11) = factor0*t22*t23
    T(11,11) = factor0*(t22*t33 + t32*t23)
    T(12,11) = factor0*t32*t33

    T(7,12) = factor0*t13*t13
    T(8,12) = factor0*2*t13*t23
    T(9,12) = factor0*2*t13*t33

    T(10,12) = factor0*t23*t23
    T(11,12) = factor0*2*t23*t33
    T(12,12) = factor0*t33*t33

    Mglo = matmul(T,Mloc)

  end subroutine EASM_transform

  !-----------------------------------------------------------------------------
  !> @author  Eva Casoni (eva.casoni@bsc.es)
  !> @date    17/09/2016
  !> @brief   Modify metric due to EAS
  !> @details
  !-------------------------------------------------------------------------------
  subroutine s8vthv(gmckov, Etilde, T)

    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    real(rp),       intent(in)             :: &
         Etilde(12),                          &   ! gpderi(1,pnodh) = dN/dxi1(inode)
         T                                        ! Integration point of the thickness
    real(rp),       intent(inout)          :: &
         gmckov(pdime,pdime)

    gmckov(1,1) = gmckov(1,1) + 2.0_rp*(Etilde(1) + T*Etilde(7))
    gmckov(2,1) = gmckov(2,1) + Etilde(2) + T*Etilde(8)
    gmckov(3,1) = gmckov(3,1) + Etilde(3) + T*Etilde(9)
    gmckov(2,2) = gmckov(2,2) + 2.0_rp*(Etilde(4) + T*Etilde(10))
    gmckov(3,2) = gmckov(3,2) + Etilde(5) + T*Etilde(11)
    gmckov(3,3) = gmckov(3,3) + 2.0_rp*(Etilde(6) + T*Etilde(12))

    ! making gmckov symmetric
    gmckov = 0.5_rp*(gmckov + transpose(gmckov))

  end subroutine s8vthv

  !-----------------------------------------------------------------------------------
  !
  ! Quadrature for the shell element
  !
  !----------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !> @author  Eva Casoni (eva.casoni@bsc.es)
  !> @date    13/10/2016
  !> @brief   Gauss points and weights for the integration of continuum shell
  !> @details
  !------------------------------------------------------------------------------
  subroutine quadrature_contshell(pgaus, pgausT, posgp, wgp, posgpT, wgpT)
    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    integer(ip),    intent(in)             :: &
         pgaus, pgausT
    real(rp),       intent(out)            :: &
         posgp(pgaus),                        &  ! Gauss points in the midplane RS
         posgpT(pgausT),                      &  ! Gauss points in the thickness
         wgp(pgaus),                          &  ! Weight of gauss points
         wgpT(pgausT)                            ! Weight of the thickness

    ! ---------------------------------------------------------------------------

    posgp(1) = -sqrt(1.0_rp/3.0_rp)
    posgp(2) = sqrt(1.0_rp/3.0_rp)

    wgp(1) = 1.0_rp
    wgp(2) = 1.0_rp

    posgpT(1) = -sqrt(1.0_rp/3.0_rp)
    posgpT(2) = sqrt(1.0_rp/3.0_rp)

    wgpT(1) = 1.0_rp
    wgpT(2) = 1.0_rp

  end subroutine quadrature_contshell

  !***************************************************************************************

  !------------------------------------------------------------------------------
  !> @author  Eva Casoni (eva.casoni@bsc.es)
  !> @date    10/04/2017
  !> @brief   Function to compute enhanced strains for being inputs to iteration i
  !> @details
  !-------------------------------------------------------------------------------
  subroutine compute_enhanced_strains(ielem,flag_SL,flag_TL,elcoor,elcoorc,elddisp,alpha)

    use def_master,            only :  ITER_K_STATE, TIME_N_STATE

    implicit none

    integer(ip), intent(in)        :: & !
         ielem,                       & !
         flag_SL, flag_TL
    real(rp), intent(in)           :: & !
         elcoor(pdime,pnode),         & !
         elcoorc(pdime,pnode),        & !
         elddisp(pdime,pnode)           ! Nodal displacement increment
    real(rp), intent(out)          :: &
         alpha(22)
    integer(ip)                    :: &
         idime, ipoin, jpoin,         &
         igaus, jgaus,                &
         igausT,                      &
         idof,i,j,inode
    real(rp)                       :: &
         dalpha(22),                  &
         ddispl(pdofs),               &
         duReor(pdofs),               &
                                ! Quadrature variables
         posgp(pgaus),                &
         posgpT(pgausT),              &
         wgp(pgaus),                  &
         wgpT(pgausT),                &
                                ! coordinates and current configuration
         elcoorci(pdime,pnode),          &
         sfref(pdime,pnode),             &
         sfcuri(pdime,pnode),            &
         a3r(pdime,pnode),               &
         a3ci(pdime,pnode),              &
         shap0i(pnodh),                  &
         dshap0i(pdimh,pnodh),           &
                                ! Variables at midsurface
         gpshapi(pnodh),                 &
         gpderii(pdimh,pnodh),           &
         arkovi(pdime,pdime),        & ! Curvilinear basis (covariant vectors) on the shell midsurface in reference configuration
         ackovi(pdime,pdime),        & ! Curvilinear basis (covariant vectors) on the shell midsurface in current configuration
         a3rkovderi(pdime,2),        & ! Derivatives of a3 in reference coordinates (of third verctor of curvilinear basis)
         a3ckovderi(pdime,2),        & ! Derivatives of a3 in current coordinates (of third verctor of curvilinear basis)
         detrai,                                  & ! Determinant of curvilinear basis in reference configuration
         detcai,                                  & ! Determinant of curvilinear basis in current configuration
                                ! metric at midsurface
         a0rkovi(pdime,pdime),                    & !
         a0rkoni(pdime,pdime),                    & !
         am0rkovi(pdime,pdime),                   & !
         am0rkoni(pdime,pdime),                   & !
         det0ri,                                  & !
                                ! metric: curvilinear
         grkovi(pdime,pdime),                     & ! Curvilinear basis (covariant vectors) on the shell in reference configuration
         grkoni(pdime,pdime),                     & !
         gckovi(pdime,pdime),                     & ! Curvilinear basis (covariant vectors) on the shell in current configuration
         gckoni(pdime,pdime),                     & !
         gmrkovi(pdime,pdime),                    & ! Metric of shell in reference configuration
         gmckovi(pdime,pdime),                    & ! Metric of shell in current configuration
         gmrkoni(pdime,pdime),                    & ! Inverse of gmrkov
         gmckoni(pdime,pdime),                    & ! Inverse of gmckon
         detrgi,                                  & ! Determinant of curvilinear basis in reference configuration
         detcgi,                                  & ! Determinant of curvilinear basis in current configuration
                                ! Magnitudes due to ANS method
                                ! alpha13 (shear locking)
         gpshap1q(2,1,pnodh),         &
         gpderi1q(2,pdimh,pnodh),     &
         ackov1q(2,pdime,pdime),      & ! Current
         a3kvpc1q(2,pdime,pdime),     &
                                ! alpha 23 (shear locking)
         gpshap2q(2,1,pnodh),         &
         gpderi2q(2,pdimh,pnodh),     &
         ackov2q(2,pdime,pdime),      & ! Current
         a3kvpc2q(2,pdime,pdime),     &
                                ! alpha 33 (trapezoidal locking)
         gpshapT(4,1,pnodh),          &
         ackovT(4,pdime,pdime),       & ! Current
                                ! B matrix variables
         fsQ(2),                      & !
         frQ(2),                      & !
         frT(pnodh),                  & !
                                ! EAS variables
         Tmati(12,12),                            & ! Transformation matrix s.t: Mglo = T*Mloc
         Mloci(12,22),                            & ! Transformation matrix in local coordinates
         Mgloi(12,22),                            & ! Matrix for enchanced strains: Etilde = Mglo*alpha
         DEASi(22,22),                            & ! DEAS = M^t*gpddsVo*M = k_alphaalpha
         DEASinvi(22,22),                         & ! Inverse of DEAS
         LEASi(22,24),                            & ! LEAS = M^t*gpddsVo*B = k_alphau
         fEASi(22),                               & ! fEAS = M* t*gpstrVo
         fEASauxi(22),                            &
         Etildei(12),                             & ! Enhanced strains: Etilde = Mglo*alpha
                                ! Matrices due to FEA
         gpbmati(12,24),                          & !
                                ! Constitutive relationship
         gpddsVoig(12,12),                        & ! Constitutive tensor resultant (integrated along thickness)
         vodds(6,6),                              & !
         gpstrVoig(12),                           & ! 2PK resultant
         stress(6),                               & ! 2PD in Voigt notation
                                ! Other Magnitudes
         xi1,xi2,xi3,                             & !
         weight1,weight2,wT,                      &
         hi(pdime), hnormi,                       &
         sumepsilon
    !
    ! Configuration at state i-1
    !
    !*************************************************************************************************
    !
    ! Compute current coordinates at state i-1
    !
    do inode = 1, pnode
       do idime = 1, pdime
          elcoorci(idime,inode) = elcoorc(idime,inode) - elddisp(idime,inode)
       end do
    end do
    idof = 0
    do inode = 1, pnode
      do idime = 1, pdime
        idof = idof + 1
        ddispl(idof) = elddisp(idime,inode)
      end do
    end do

    !
    ! Incremental displacement: OJO!!! AQUI HAY QUE PONER EL ORDEN DEL ELEMENTO SHELL!!!
    !        * elddisp has the standar ordering (1234-bottom nodes  5678-top nodes)
    !        * elddispReor has to be reorder s.t:
    !
    !
    ! Reorder elddispReor to solve the system at the end
    !
    ! Standard node (pdofs) ordering:
    !
    ! Bottom face: 1----------4        (10-11-12)---------(7-8-9)
    !              |          |            |                 |
    !              |          |            |                 |
    !              2----------3         (1-2-3)------------(4-5-6)
    !
    ! Top face:    5----------8        (22-23-24)--------(19-20-21)
    !              |          |            |                 |
    !              |          |            |                 |
    !              6----------7        (13-14-15)--------(16-17-18)
    !
    ! Shell (pdofs) ordering:
    !
    ! Bottom face: 2----------8        (16-17-18)--------(22-23-24)
    !              |          |            |                 |
    !              |          |            |                 |
    !              4----------6         (4-5-6)----------(10-11-12)
    !
    ! Top face:    1----------7        (13-14-15)--------(19-20-21)
    !              |          |            |                 |
    !              |          |            |                 |
    !              3----------5         (1-2-3)-----------(7-8-9)
    !

    duReor(1) = ddispl(13)
    duReor(2) = ddispl(14)
    duReor(3) = ddispl(15)

    duReor(4) = ddispl(1)
    duReor(5) = ddispl(2)
    duReor(6) = ddispl(3)

    duReor(7) = ddispl(16)
    duReor(8) = ddispl(17)
    duReor(9) = ddispl(18)

    duReor(10) = ddispl(4)
    duReor(11) = ddispl(5)
    duReor(12) = ddispl(6)

    duReor(13) = ddispl(19)
    duReor(14) = ddispl(20)
    duReor(15) = ddispl(21)

    duReor(16) = ddispl(7)
    duReor(17) = ddispl(8)
    duReor(18) = ddispl(9)

    duReor(19) = ddispl(22)
    duReor(20) = ddispl(23)
    duReor(21) = ddispl(24)

    duReor(22) = ddispl(10)
    duReor(23) = ddispl(11)
    duReor(24) = ddispl(12)

    !

    ! Midsurface definition in reference and updated configuration
    !   - coordinates
    !<GGU> Aixo ja esta calculat a fora
    sfref(:,:) = 0.0_rp
    do ipoin = 1, pnodh
       jpoin = ipoin + pnodh
       do idime = 1, pdime
          sfref(idime,ipoin) = 0.5_rp*(elcoor(idime,ipoin) + elcoor(idime,jpoin))
          a3r(idime,ipoin) = 0.5_rp*(elcoor(idime,ipoin+pnodh) - elcoor(idime,ipoin))
       end do
    end do
    !
    ! Compute midsurface coordinates and normal vector to midsurface
    !
    sfcuri(:,:) = 0.0_rp
    do ipoin = 1, pnodh
       jpoin = ipoin + pnodh
       do idime = 1, pdime
          sfcuri(idime,ipoin) = 0.5_rp*(elcoorci(idime,ipoin)       + elcoorci(idime,jpoin))
          a3ci(idime,ipoin)   = 0.5_rp*(elcoorci(idime,ipoin+pnodh) - elcoorci(idime,ipoin))
       end do
    end do
    !
    ! Computations over the midpoint of the shell (gausspoint = (0,0)):
    ! curvilinear basis in reference configuration
    !
    call shape_functions(0.0_rp, 0.0_rp, shap0i(:), dshap0i(:,:))
    !
    ! Curvilinear shell basis only for the reference configuration
    !
    call g_shell_basis(sfref(:,:), a3r(:,:), shap0i(:), dshap0i(:,:), 0.0_rp, &
         a0rkovi(:,:), a0rkoni(:,:), am0rkovi(:,:), am0rkoni(:,:), det0ri, ielem)
    !
    ! Metrics for Shear and Trapezoidal Lockings
    !
    ! Shear Locking
    !
    if ( flag_SL == SLD_CSHEL_ANS_SHEAR ) then
       call ANScolocationpointsShearQ(ielem,sfcuri(:,:),a3ci(:,:),          &
            gpshap1q(:,:,:),gpderi1q(:,:,:),ackov1q(:,:,:),a3kvpc1q(:,:,:), &
            gpshap2q(:,:,:),gpderi2q(:,:,:),ackov2q(:,:,:),a3kvpc2q(:,:,:))
    end if
    !
    ! Trapezoidal Locking
    !
    if ( flag_TL == SLD_CSHEL_ANS_TRAPEZOIDAL ) then
       call ANScolocationpointsShearT(ielem,sfcuri(:,:),a3ci(:,:), &
            gpshapT(:,:,:),ackovT(:,:,:))
    end if
    !
    ! Recover state variables from last increment
    !
    do i = 1,22
       alpha(i) = svegm_sld(ielem)%a(i,1_ip,TIME_N_STATE)
    end do
    !
    ! Quadrature
    !
    call quadrature_contshell(pgaus,pgausT,posgp,wgp,posgpT,wgpT)
    !
    ! Initialize matrices
    !
    DEASi = 0.0_rp
    LEASi = 0.0_rp
    fEASi = 0.0_rp
    !
    ! Loop on gauss points of the surface
    !
    do igaus = 1,2
       xi1 = posgp(igaus)
       weight1 = wgp(igaus)
       do jgaus = 1,2
          xi2 = posgp(jgaus)
          weight2 = wgp(jgaus)
          !
          ! Interpolating function at the middle plane and its derivative (N and dN/dn)
          !
          call shape_functions(xi1, xi2,  gpshapi(:), gpderii(:,:))
          !
          ! Basis in the midsurface of the shell in reference and current configuration
          !
          call a_midsurface_basis(sfref(:,:), gpshapi(:) ,gpderii(:,:), a3r(:,:),   &
               arkovi(:,:), a3rkovderi(:,:), detrai, ielem)
          call a_midsurface_basis(sfcuri(:,:), gpshapi(:), gpderii(:,:), a3ci(:,:), &
               ackovi(:,:), a3ckovderi(:,:), detcai, ielem)
          !
          ! Compute H norm
          !
          hnormi = 0.0_rp
          hi(1) = arkovi(2,1)*arkovi(3,2) - arkovi(3,1)*arkovi(2,2)
          hi(2) = arkovi(3,1)*arkovi(1,2) - arkovi(1,1)*arkovi(3,2)
          hi(3) = arkovi(1,1)*arkovi(2,2) - arkovi(2,1)*arkovi(1,2)
          do idime=1, pdime
             hnormi = hnormi + hi(idime)*hi(idime)
          end do
          hnormi = sqrt(hnormi)
          if (hnormi > 0.0_rp) then
             hi = hi/hnormi
          else
             write(*,*) 'Norm of the thickness is 0'
          end if
          !
          ! M matrix (EAS method) in the local coordinate system (M_xi in order to obtain enhanced strains)
          !
          call EASM_local(xi1, xi2, Mloci(:,:))
          !
          ! Transforming to the global coordinate system (M = (detJ0/detJ)*(T0^-t)*M_xi
          !
          call EASM_transform(Mloci(:,:), Mgloi(:,:), arkovi(:,:), a0rkoni(:,:), detrai, det0ri, Tmati(:,:))
          !
          ! Compute compatible strains: Etilde = M*alpha
          !
          do i = 1, 12
             sumepsilon = 0.0_rp
             do j = 1, 22
                sumepsilon = sumepsilon + Mgloi(i,j)*alpha(j)
             end do
             Etildei(i) = sumepsilon
          end do
          !
          ! B matrix (B)
          !
          call b_matrix(gpshapi(:), gpderii(:,:), ackovi(:,:), a3ckovderi(:,:), gpbmati(:,:))
          !
          ! Modify B matrix due to ANS Shear Locking or Trapezoidal Locking
          !
          if (flag_SL == SLD_CSHEL_ANS_SHEAR ) then
             call ANS_bmatrix_shearLocking(xi1,xi2,frQ(:),fsQ(:))
             call b_matrixANS_shearLocking(&
                  gpshap1q(:,:,:),gpderi1q(:,:,:),ackov1q(:,:,:), &
                  gpshap2q(:,:,:),gpderi2q(:,:,:),ackov2q(:,:,:), &
                  frQ,fsQ,gpbmati(:,:))
          end if

          if ( flag_TL == SLD_CSHEL_ANS_TRAPEZOIDAL ) then
             call ANS_bmatrix_trapezoidalLocking(xi1,xi2,frT(:))
             call b_matrixANS_trapezoidalLocking(gpshapT(:,:,:),ackovT(:,:,:),frT(:),gpbmati(:,:))
          end if
          !
          ! Initializations for thickness integration
          !
          gpddsVoig = 0.0_rp
          gpstrVoig = 0.0_rp
          !
          ! Loop over the interpolating points along the thickness
          !
          gausspointsthicknessi: do igausT = 1, pgausT
             xi3 = posgpT(igausT)
             wT = wgpT(igausT)
             !
             ! Curvilinear basis in the shell
             !
             call g_shell_basis(sfref(:,:), a3r(:,:), gpshapi(:), gpderii(:,:), xi3,   &
                  grkovi(:,:), grkoni(:,:), gmrkovi(:,:), gmrkoni(:,:), detrgi, ielem)
             call g_shell_basis(sfcuri(:,:), a3ci(:,:), gpshapi(:), gpderii(:,:), xi3, &
                  gckovi(:,:), gckoni(:,:), gmckovi(:,:), gmckoni(:,:),detcgi, ielem)
             !
             ! Change the metrics due to EAS
             !
             call s8vthv(gmckovi(:,:), Etildei, xi3)
             !
             ! Constitutive relationship tensors
             !
             call material_law(ielem, gmrkovi(:,:), gmckovi(:,:), gmrkoni(:,:), gmckoni(:,:), &
                  grkoni(:,:), detrgi, detcgi, vodds(:,:), stress(:))
             !
             ! Integrate material law
             !
             wT = wT*detrgi/hnormi
             call integrate_material_law(vodds(:,:),stress(:),wT,xi3,gpddsVoig(:,:),gpstrVoig(:))

          end do gausspointsthicknessi
          !
          ! Computation of auxiliar matrices: k_alphaalpha=D=M^t*C*M , k_alphau=L=M^t*C*B
          !
          call compute_matrices(gpddsVoig(:,:),gpbmati(:,:),Mgloi(:,:),weight1,weight2,hnormi,DEASi(:,:),LEASi(:,:))
          !
          ! Compute external forces: f_EAS = int(M^t*S)
          !
          call compute_EAS_force(gpstrVoig(:),Mgloi(:,:),weight1,weight2,hnormi,fEASi)
          !
       end do
    end do  ! end do gauss points of the surface
    !
    ! Solve the system for intermediate enhaning strains:
    ! Dalpha = -(D^-1)*(fEAS + L*du) of size 22x1
    !
    DEASinvi = DEASi
    call invert(DEASinvi,22_ip,22_ip)
    !
    ! Update alpha
    !
    fEASauxi = matmul(LEASi,duReor)
    fEASauxi = fEASauxi + fEASi
    dalpha   = -matmul(DEASinvi,fEASauxi)
    !
    ! Update alpha
    !
    alpha = alpha + dalpha
    !
    ! Store enhanced strains
    !
    svegm_sld(ielem)%a(:,1_ip,ITER_K_STATE) = alpha

  end subroutine compute_enhanced_strains

  !************************************************************************************
  !
  ! ANS FUNCTIONS
  !
  !***********************************************************************************
  !------------------------------------------------------------------------------
  !> @author  Eva Casoni and Gerard Guillamet (eva.casoni@bsc.es)
  !> @date    10/04/2017
  !> @brief   Compute metrics ANS shear functions
  !> @details Math: modify shear components of Green-Lagrange strain tensor
  !------------------------------------------------------------------------------
  subroutine ANScolocationpointsShearQ(ielem,sfcur,a3c,gpshap1q,gpderi1q, &
       ackov1q,a3kvpc1q,gpshap2q,gpderi2q,ackov2q,a3kvpc2q)
    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    real(rp),       intent(in)             ::    &
         sfcur(pdime,pnodh),                     &
         a3c(pdime,pnodh)
    integer(ip),    intent(in)             ::    &
         ielem
    real(rp),       intent(out)            ::    &
                                ! alpha 13
         gpshap1q(2,1,pnodh),                    &
         gpderi1q(2,pdimh,pnodh),                &
         ackov1q(2,pdime,pdime),                 & ! Current
         a3kvpc1q(2,pdime,pdime),                &
                                ! alpha 23
         gpshap2q(2,1,pnodh),                    &
         gpderi2q(2,pdimh,pnodh),                &
         ackov2q(2,pdime,pdime),                 & ! Current
         a3kvpc2q(2,pdime,pdime)
    integer(ip)                            ::    &
         iq,j,k
    real(rp)                               ::    &
         xr1(2),xs1(2),xr2(2),xs2(2),            &
                                ! Auxiliary variables
         gpshap1Qaux(pnodh),                     &
         gpderi1Qaux(pdimh,pnodh),               &
         gpshap2Qaux(pnodh),                     &
         gpderi2Qaux(pdimh,pnodh),               &
                                ! alpha 13
         ackov1Qaux(pdime,pdime),                & ! Current
         a3ckov1Qaux(pdime,2),                   &
         detc1Qaux,                              &
                                ! alpha 23
         ackov2Qaux(pdime,pdime),                & ! Current
         a3ckov2Qaux(pdime,2),                   &
         detc2Qaux
    !
    ! Coordinates of the collocation points for ANS transverse shear
    !
    ! Points (0,-1), (0,1) alpha13
    xr1(1) =  0.0_rp
    xs1(1) = -1.0_rp
    !
    xr1(2) =  0.0_rp
    xs1(2) =  1.0_rp
    !
    !Points (-1,0), (1,0) alpha23
    xr2(1) = -1.0_rp
    xs2(1) =  0.0_rp
    !
    xr2(2) =  1.0_rp
    xs2(2) =  0.0_rp
    !
    ! Loop over collocation points (alpha 13)
    !
    do iq=1,2
       !
       ! Evaluation of the shape functions and their derivatives
       ! at calculation points
       !
       call shape_functions(xr1(iq),xs1(iq),gpshap1Qaux(:),gpderi1Qaux(:,:))
       !
       ! Auxiliar shape function for each point
       !
       do j=1,pnodh
          gpshap1q(iq,1,j) = gpshap1Qaux(j)
       end do
       do j=1,pnodh
          do k=1,2
             gpderi1q(iq,k,j) =  gpderi1Qaux(k,j)
          end do
       end do
       !
       ! Metrics in the current configuration
       !
       call a_midsurface_basis(sfcur(:,:), gpshap1Qaux(:), gpderi1Qaux(:,:), a3c(:,:), &
            ackov1Qaux(:,:), a3ckov1Qaux(:,:), detc1Qaux, ielem)
       !
       ! Restoring the magnitudes
       !
       do j=1,pdime
          do k=1,pdime
             ackov1q(iq,j,k)  =  ackov1Qaux(j,k)
          end do
       end do
       do j=1,pdime
          do k=1,2
             a3kvpc1q(iq,j,k) =  a3ckov1Qaux(j,k)
          end do
       end do

    end do
    !
    ! Loop over collocation points (alpha 23)
    !
    do iq=1,2
       !
       ! Evaluation of the shape functions and their derivatives
       ! at calculation points
       !
       call shape_functions(xr2(iq),xs2(iq),gpshap2Qaux(:),gpderi2Qaux(:,:))

       ! Restoring the magnitdes
       do j=1,pnodh
          gpshap2q(iq,1,j) = gpshap2Qaux(j)
       end do
       do j=1,pnodh
          do k=1,2
             gpderi2q(iq,k,j) =  gpderi2Qaux(k,j)
          end do
       end do
       !
       ! Metrics in the current configuration
       !
       call a_midsurface_basis(sfcur(:,:), gpshap2Qaux(:), gpderi2Qaux(:,:), a3c(:,:), &
            ackov2Qaux(:,:), a3ckov2Qaux(:,:), detc2Qaux, ielem)
       !
       ! Restoring the magnitudes
       !
       do j=1,pdime
          do k=1,pdime
             ackov2q(iq,j,k)  =  ackov2Qaux(j,k)
          end do
       end do
       do j=1,pdime
          do k=1,2
             a3kvpc2q(iq,j,k) =  a3ckov2Qaux(j,k)
          end do
       end do

    end do

  end subroutine ANScolocationpointsShearQ

  !------------------------------------------------------------------------------
  !> @author  Eva Casoni (eva.casoni@bsc.es)
  !> @date    6/04/2017
  !> @brief   Compute colocation points ANS Trapezoidal locking
  !> @details Math: modify shear components of Green-Lagrange strain tensor
  !------------------------------------------------------------------------------
  subroutine ANScolocationpointsShearT(ielem,sfcur,a3c,gpshapT,ackovT)
    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    real(rp),       intent(in)                :: &
         sfcur(pdime,pnodh),                     &
         a3c(pdime,pnodh)
    integer(ip),    intent(in)                :: &
         ielem
    real(rp),       intent(out)               :: &
         gpshapT(4,1,pnodh),                     & ! alpha33 (curvature thickness locking)
         ackovT(4,pdime,pdime)                     ! Current
    real(rp)                                  :: &
         gpshapTaux(4),                          &
         gpderiTaux(2,4),                        &
         ackovTaux(pdime,pdime),                 &
         a3ckovTaux(pdime,2),                    &
         detcTaux
    integer(ip)                               :: &
         qt,j,k
    real(rp)                                  :: &
         xrT(4),xsT(4)
    !
    ! Coordinates of the collocation points for ANS curvature locking
    !
    ! Points (1,1), (-1,1) alpha33
    !
    xrT(1) =  1.0_rp
    xrT(2) = -1.0_rp
    !
    xsT(1) =  1.0_rp
    xsT(2) =  1.0_rp
    !
    ! Points (-1,-1), (1,-1) alpha23
    !
    xrT(3) = -1.0_rp
    xrT(4) =  1.0_rp
    !
    xsT(3) = -1.0_rp
    xsT(4) = -1.0_rp
    !
    ! Loop over the collocation points
    !
    do qt=1,4
       !
       ! Shape functions
       !
       call shape_functions(xrT(qt),xsT(qt),gpshapTaux(:),gpderiTaux(:,:))
       !
       ! Auxiliar shape function for each point
       !
       do j=1,pnodh
          gpshapT(qt,1,j) = gpshapTaux(j)
       end do
       !
       ! Metrics in the current configuration
       !
       call a_midsurface_basis(sfcur(:,:), gpshapTaux(:), gpderiTaux(:,:), a3c(:,:), &
            ackovTaux(:,:), a3ckovTaux(:,:), detcTaux, ielem)
       !
       ! Auxiliar metric for each point
       !
       do j=1,pdime
          do k=1,pdime
             ackovT(qt,j,k) = ackovTaux(j,k)
          end do
       end do

    end do

  end subroutine ANScolocationpointsShearT

  !------------------------------------------------------------------------------
  !> @author  Eva Casoni (eva.casoni@bsc.es)
  !> @date    10/10/2016
  !> @brief   Compute modified Bmatrix with the modified ANS shear functions
  !> @details Math: modify shear components of Green-Lagrange strain tensor
  !------------------------------------------------------------------------------
  subroutine ANS_bmatrix_shearLocking(r,s,frQ,fsQ)

    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    real(rp),       intent(in)             :: &
         r,s                                    ! Local coordinates xi1,xi2
    real(rp),       intent(out)            :: &
         frQ(2), fsQ(2)

    frQ(1) = 0.5_rp*(1.0_rp - s)
    frQ(2) = 0.5_rp*(1.0_rp + s)

    fsQ(1) = 0.5_rp*(1.0_rp - r)
    fsQ(2) = 0.5_rp*(1.0_rp + r)

  end subroutine ANS_bmatrix_shearLocking

  !------------------------------------------------------------------------------
  !> @author  Eva Casoni (eva.casoni@bsc.es)
  !> @date    10/10/2016
  !> @brief   Compute modified Bmatrix with the modified ANS shear functions
  !> @details Math: modify shear components of Green-Lagrange strain tensor
  !------------------------------------------------------------------------------
  subroutine ANS_bmatrix_trapezoidalLocking(r,s,frT)

    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    real(rp),       intent(in)             :: &
         r,s                                     ! Local coordinates xi1,xi2
    real(rp),       intent(out)            :: &
         frT(pnodh)

    frT(1) = ((1.0_rp + r)*(1.0_rp + s))/4.0_rp
    frT(2) = ((1.0_rp - r)*(1.0_rp + s))/4.0_rp
    frT(3) = ((1.0_rp - r)*(1.0_rp - s))/4.0_rp
    frT(4) = ((1.0_rp + r)*(1.0_rp - s))/4.0_rp

  end subroutine ANS_bmatrix_trapezoidalLocking

  !-------------------------------------------------------------------------------
  !> @author Gerard Guillamet, Eva CAsoni
  !> @date    10/04/2017
  !> @brief
  !> @details Math: modify B operator due to ANS shear method
  !------------------------------------------------------------------------------
  subroutine b_matrixANS_shearLocking(shapef1q,dshape1q,akovc1q,  &
       shapef2q,dshape2q,akovc2q,frq,fsq,Bmat)
    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    real(rp), intent(in)        :: &
                                     ! alpha 13
         shapef1q(2,1,4),          &
         dshape1q(2,2,4),          &
         akovc1q(2,3,3),           &
                                     ! alpha 23
         shapef2q(2,1,4),          &
         dshape2q(2,2,4),          &
         akovc2q(2,3,3),           &
         frq(2),                   &
         fsq(2)
    real(rp), intent(inout)     :: &
         Bmat(12,24)
    integer(ip)                 :: &
         inode,jq,node_start
    real(rp)                    :: &
         a1x1q, a1y1q, a1z1q,      &
         a3x1q, a3y1q, a3z1q,      &
         a2x2q, a2y2q, a2z2q,      &
         a3x2q, a3y2q, a3z2q,      &
         frqj,fsqj,                &
         pk1q, pk2q,               &
         pk1, pk2
    !
    ! Modigy Bmat due to ANS Shear Locking
    !
    node_start = 1_ip
    do inode=1,4
       !
       ! Initialization
       !
       ! ANS alpha13
       Bmat(3,node_start)   = 0.0_rp
       Bmat(3,node_start+1) = 0.0_rp
       Bmat(3,node_start+2) = 0.0_rp
       Bmat(3,node_start+3) = 0.0_rp
       Bmat(3,node_start+4) = 0.0_rp
       Bmat(3,node_start+5) = 0.0_rp
       ! ANS alpha23
       Bmat(5,node_start) =   0.0_rp
       Bmat(5,node_start+1) = 0.0_rp
       Bmat(5,node_start+2) = 0.0_rp
       Bmat(5,node_start+3) = 0.0_rp
       Bmat(5,node_start+4) = 0.0_rp
       Bmat(5,node_start+5) = 0.0_rp

       do jq =1,2
          a1x1q = akovc1q(jq,1,1)
          a1y1q = akovc1q(jq,2,1)
          a1z1q = akovc1q(jq,3,1)

          a3x1q = akovc1q(jq,1,3)
          a3y1q = akovc1q(jq,2,3)
          a3z1q = akovc1q(jq,3,3)

          a2x2q = akovc2q(jq,1,2)
          a2y2q = akovc2q(jq,2,2)
          a2z2q = akovc2q(jq,3,2)

          a3x2q = akovc2q(jq,1,3)
          a3y2q = akovc2q(jq,2,3)
          a3z2q = akovc2q(jq,3,3)

          frqj = frq(jq)
          fsqj = fsq(jq)
          pk1q = shapef1q(jq,1,inode)
          pk2q = shapef2q(jq,1,inode)
          pk1  = dshape1q(jq,1,inode)
          pk2  = dshape2q(jq,2,inode)
          !
          ! ANS alpha 13
          !
          Bmat(3,node_start)   = Bmat(3,node_start)   + 0.5_rp*(pk1*a3x1q*frqj + pk1q*a1x1q*frqj)
          Bmat(3,node_start+1) = Bmat(3,node_start+1) + 0.5_rp*(pk1*a3y1q*frqj + pk1q*a1y1q*frqj)
          Bmat(3,node_start+2) = Bmat(3,node_start+2) + 0.5_rp*(pk1*a3z1q*frqj + pk1q*a1z1q*frqj)
          Bmat(3,node_start+3) = Bmat(3,node_start+3) + 0.5_rp*(pk1*a3x1q*frqj - pk1q*a1x1q*frqj)
          Bmat(3,node_start+4) = Bmat(3,node_start+4) + 0.5_rp*(pk1*a3y1q*frqj - pk1q*a1y1q*frqj)
          Bmat(3,node_start+5) = Bmat(3,node_start+5) + 0.5_rp*(pk1*a3z1q*frqj - pk1q*a1z1q*frqj)
          !
          ! ANS alpha 23
          !
          Bmat(5,node_start)   = Bmat(5,node_start)   + 0.5_rp*(pk2*a3x2q*fsqj + pk2q*a2x2q*fsqj)
          Bmat(5,node_start+1) = Bmat(5,node_start+1) + 0.5_rp*(pk2*a3y2q*fsqj + pk2q*a2y2q*fsqj)
          Bmat(5,node_start+2) = Bmat(5,node_start+2) + 0.5_rp*(pk2*a3z2q*fsqj + pk2q*a2z2q*fsqj)
          Bmat(5,node_start+3) = Bmat(5,node_start+3) + 0.5_rp*(pk2*a3x2q*fsqj - pk2q*a2x2q*fsqj)
          Bmat(5,node_start+4) = Bmat(5,node_start+4) + 0.5_rp*(pk2*a3y2q*fsqj - pk2q*a2y2q*fsqj)
          Bmat(5,node_start+5) = Bmat(5,node_start+5) + 0.5_rp*(pk2*a3z2q*fsqj - pk2q*a2z2q*fsqj)

       end do

       node_start = node_start + 6_ip

    end do

  end subroutine b_matrixANS_shearLocking

  !------------------------------------------------------------------------------
  !> @author  Eva Casoni
  !> @date    06/04/2017
  !> @brief   Compute modified B matrix due to ANS trapezoidal locking
  !> @details
  !-------------------------------------------------------------------------------
  subroutine b_matrixANS_trapezoidalLocking(gpshapT,ackovT,frT,Bmat)

    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    real(rp), intent(in)        :: &
         gpshapT(4,1,4),           &
         ackovT(4,pdime,pdime),    &
         frT(pnodh)
    real(rp), intent(inout)     :: &
         Bmat(12,24)
    integer(ip)                 :: &
         inode,jt,node_start
    real(rp)                    :: &
         a3xT,a3yT,a3zT,           &
         pkT

    node_start = 1_ip

    do inode=1,4
       !
       ! ANS alpha33
       !
       Bmat(6,node_start)   = 0.0_rp
       Bmat(6,node_start+1) = 0.0_rp
       Bmat(6,node_start+2) = 0.0_rp
       Bmat(6,node_start+3) = 0.0_rp
       Bmat(6,node_start+4) = 0.0_rp
       Bmat(6,node_start+5) = 0.0_rp

       do jt=1,4

          a3xT = ackovT(jt,1,3)
          a3yT = ackovT(jt,2,3)
          a3zT = ackovT(jt,3,3)

          pkT =  gpshapT(jt,1,inode)

          Bmat(6,node_start)   = Bmat(6,node_start)   + 0.5_rp*(pkT*a3xT*frT(jt))
          Bmat(6,node_start+1) = Bmat(6,node_start+1) + 0.5_rp*(pkT*a3yT*frT(jt))
          Bmat(6,node_start+2) = Bmat(6,node_start+2) + 0.5_rp*(pkT*a3zT*frT(jt))
          Bmat(6,node_start+3) = Bmat(6,node_start+3) - 0.5_rp*(pkT*a3xT*frT(jt))
          Bmat(6,node_start+4) = Bmat(6,node_start+4) - 0.5_rp*(pkT*a3yT*frT(jt))
          Bmat(6,node_start+5) = Bmat(6,node_start+5) - 0.5_rp*(pkT*a3zT*frT(jt))

       end do

       node_start = node_start + 6_ip

    end do

  end subroutine b_matrixANS_trapezoidalLocking

  !------------------------------------------------------------------------------
  !> @author  Gerard Guillamet and Eva Casoni
  !> @date    05/04/2017
  !> @brief   Compute geometrical part of k_uu: int(dB/du * stress)
  !> @details ANS Shear and Trapezoidal lockings available
  !-------------------------------------------------------------------------------
  subroutine computeGeomStiffMatrixGlobal(stress,gpshap,gpderi,w1,w2,hnormi, &
       gpshap13SL,gpderi13SL,gpshap23SL,gpderi23SL,gpshap33TL,frq,fsq,frT,   &
       flag_SL,flag_TL,stiff)

    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    integer(ip), intent(in)                ::    &
         flag_SL,                                &    ! Flag ANS Shear Locking (SL)
         flag_TL                                      ! Flag ANS Trapezoidal Locking (TL)
    real(rp),       intent(in)             ::    &
         stress(12),                             &    ! 2PK
         gpshap(pnodh),                          &    ! Shape functions at nodes of midshell
         gpderi(pdimh,pnodh),                    &    ! gpderi(1,pnodh) = dN/dxi1(inode)
         gpshap13SL(2,1,pnodh),                  &    ! Shape functions SL 13 component
         gpderi13SL(2,2,pnodh),                  &    ! Derivatives shape functions SL 13 component
         gpshap23SL(2,1,pnodh),                  &    ! Shape functions SL 23 component
         gpderi23SL(2,2,pnodh),                  &    ! Derivatives shape functions SL 23 component
         gpshap33TL(pnodh,1,pnodh),              &    ! Shape functions TL 33 component
         frT(pnodh),                             &    !
         frq(2),                                 &    !
         fsq(2),                                 &    !
         w1,w2,hnormi
    real(rp),       intent(out)            ::    &    !
         stiff(24,24)
    real(rp)                               ::    &
         Ni, Nj, dNidr, dNjdr, dNids, dNjds,     &    ! Shape functions and derivatives
         dNidrNj_SL, NidNjdr_SL, dNidsNj_SL,     &    ! Shape functions and derivatives affecting SL and TL
         NidNjds_SL, NiNj_TL,                    &
         n11,n12,n13,n22,n23,n33,                &    ! Stress constant components
         m11,m12,m13,m22,m23,m33,                &    ! Stress linear components
         aIK, bIK, cIK, dIK,                     &    ! Auxiliar terms of the matrix
         factor                                       ! w1*w2*hnormi
    integer(ip)                            ::    &
         inodh, jnodh, irow, jcol, ic, inode
    !
    ! Stress components
    ![n m] = [n11 n22 ... | m11 m22 ...]
    !         constant     linear
    n11 = stress(1)
    n12 = stress(2)
    n13 = stress(3)
    n22 = stress(4)
    n23 = stress(5)
    n33 = stress(6)
    m11 = stress(7)
    m12 = stress(8)
    m13 = stress(9)
    m22 = stress(10)
    m23 = stress(11)
    m33 = stress(12)
    !
    ! Initialization
    !
    stiff  = 0.0_rp
    factor = w1*w2*hnormi

    irow = 1
    do inodh = 1, pnodh
       jcol = 1
       do jnodh = 1, inodh
          !
          ! Shape functions and their derivatives
          !
          Ni    = gpshap(inodh)
          Nj    = gpshap(jnodh)
          dNidr = gpderi(1,inodh)
          dNjdr = gpderi(1,jnodh)
          dNids = gpderi(2,inodh)
          dNjds = gpderi(2,jnodh)
          !
          ! Transverse Shear Locking (SL)
          ! The interpolation scheme due to ANS method only affects to n13 and n23
          !
          if (flag_SL .ne. 0_ip) then
             dNidrNj_SL = 0.0_rp
             NidNjdr_SL = 0.0_rp
             dNidsNj_SL = 0.0_rp
             NidNjds_SL = 0.0_rp
             do ic=1,2
                dNidrNj_SL = dNidrNj_SL + gpshap13SL(ic,1,jnodh)*gpderi13SL(ic,1,inodh)*frq(ic)
                NidNjdr_SL = NidNjdr_SL + gpshap13SL(ic,1,inodh)*gpderi13SL(ic,1,jnodh)*frq(ic)
                dNidsNj_SL = dNidsNj_SL + gpshap23SL(ic,1,jnodh)*gpderi23SL(ic,2,inodh)*fsq(ic)
                NidNjds_SL = NidNjds_SL + gpshap23SL(ic,1,inodh)*gpderi23SL(ic,2,jnodh)*fsq(ic)
             end do
          else
             dNidrNj_SL = Nj*dNidr  ! n13
             NidNjdr_SL = Ni*dNjdr  ! n13
             dNidsNj_SL = Nj*dNids  ! n23
             NidNjds_SL = Ni*dNjds  ! n23
          end if
          !
          ! Trapezoidal Locking (TL)
          ! The interpolaticon scheme due to ANS method only affects n33
          !
          if (flag_TL .ne. 0_ip) then
             NiNj_TL = 0.0_rp
             do inode=1,pnodh
                NiNj_TL = NiNj_TL + gpshap33TL(inode,1,inodh)*gpshap33TL(inode,1,jnodh)*frT(inode)
             end do
          else
             NiNj_TL = Ni*Nj
          end if
          !
          ! Auxiliar terms of the geometric stiffnes matix
          !
          aIK = 0.25_rp*(dNidr*dNjdr)*n11              + 0.25_rp*(dNids*dNjds)*n22               &
               + 0.25_rp*(NiNj_TL)*n33                 + 0.25_rp*(dNidsNj_SL + NidNjds_SL)*n23   &
               + 0.25_rp*(dNidrNj_SL + NidNjdr_SL)*n13 + 0.25_rp*(dNids*dNjdr + dNidr*dNjds)*n12 &
               + 0.50_rp*(dNidr*dNjdr)*m11             + 0.50_rp*(dNids*dNjds)*m22               &
               + 0.25_rp*(Ni*dNjds + dNids*Nj)*m23     + 0.25_rp*(Ni*dNjdr + dNidr*Nj)*m13       &
               + 0.50_rp*(dNids*dNjdr + dNidr*dNjds)*m12

          cIK = 0.25_rp*(dNidr*dNjdr)*n11              + 0.25_rp*(dNids*dNjds)*n22               &
               - 0.25_rp*(NiNj_TL)*n33                 + 0.25_rp*(dNidsNj_SL - NidNjds_SL)*n23   &
               + 0.25_rp*(dNidrNj_SL - NidNjdr_SL)*n13 + 0.25_rp*(dNids*dNjdr + dNidr*dNjds)*n12 &
               - 0.25_rp*(Ni*dNjds + dNids*Nj)*m23     - 0.25_rp*(Ni*dNjdr + dNidr*Nj)*m13

          bIK = 0.25_rp*(dNidr*dNjdr)*n11              + 0.25_rp*(dNids*dNjds)*n22               &
               - 0.25_rp*(NiNj_TL)*n33                 - 0.25_rp*(dNidsNj_SL - NidNjds_SL)*n23   &
               - 0.25_rp*(dNidrNj_SL - NidNjdr_SL)*n13 + 0.25_rp*(dNids*dNjdr + dNidr*dNjds)*n12 &
               - 0.25_rp*(Ni*dNjds + dNids*Nj)*m23     - 0.25_rp*(Ni*dNjdr + dNidr*Nj)*m13

          dIK = 0.25_rp*(dNidr*dNjdr)*n11              + 0.25_rp*(dNids*dNjds)*n22               &
               + 0.25_rp*(NiNj_TL)*n33                 - 0.25_rp*(dNidsNj_SL + NidNjds_SL)*n23   &
               - 0.25_rp*(dNidrNj_SL + NidNjdr_SL)*n13 + 0.25_rp*(dNids*dNjdr + dNidr*dNjds)*n12 &
               - 0.50_rp*(dNidr*dNjdr)*m11             - 0.50_rp*(dNids*dNjds)*m22               &
               + 0.25_rp*(Ni*dNjds + dNids*Nj)*m23     + 0.25_rp*(Ni*dNjdr + dNidr*Nj)*m13       &
               - 0.50_rp*(dNidr*dNjds + dNids*dNjdr)*m12

          stiff(irow,jcol)     = aIK*factor    ! First row first column
          stiff(irow+1,jcol+1) = aIK*factor
          stiff(irow+2,jcol+2) = aIK*factor

          stiff(irow,jcol+3)   = bIK*factor    ! First row second column
          stiff(irow+1,jcol+4) = bIK*factor
          stiff(irow+2,jcol+5) = bIK*factor

          stiff(irow+3,jcol)   = cIK*factor    ! Second row first column
          stiff(irow+4,jcol+1) = cIK*factor
          stiff(irow+5,jcol+2) = cIK*factor

          stiff(irow+3,jcol+3) = dIK*factor    ! Second row second column
          stiff(irow+4,jcol+4) = dIK*factor
          stiff(irow+5,jcol+5) = dIK*factor

          if (inodh.ne.jnodh) then             ! Non-diagonal terms
             stiff(jcol,irow)     = aIK*factor
             stiff(jcol+1,irow+1) = aIK*factor
             stiff(jcol+2,irow+2) = aIK*factor

             stiff(jcol+3,irow)   = bIK*factor
             stiff(jcol+4,irow+1) = bIK*factor
             stiff(jcol+5,irow+2) = bIK*factor

             stiff(jcol,irow+3)   = cIK*factor
             stiff(jcol+1,irow+4) = cIK*factor
             stiff(jcol+2,irow+5) = cIK*factor

             stiff(jcol+3,irow+3) = dIK*factor
             stiff(jcol+4,irow+4) = dIK*factor
             stiff(jcol+5,irow+5) = dIK*factor
          end if

          jcol = jcol + 6
       end do
       irow = irow + 6
    end do

  end subroutine computeGeomStiffMatrixGlobal

  !**********************************************************************************

  !-----------------------------------------------------------------------------
  !> @author  Eva Casoni (eva.casoni@bsc.es)
  !> @date    17/09/2016
  !> @brief   Reorder matrix
  !> @details
  !
  !   Kmatrix                  Kmatrix_reorder
  !
  !   1     3  (top)              3      4
  !   |     |                     |      |
  !   |     |                     |      |
  !   2     4  (bottom)           1      2
  !
  !-------------------------------------------------------------------------------
  subroutine K_element_reorder(Kmat,Kreorder)

    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    real(rp),       intent(in)             :: &
         Kmat(pdofs,pdofs)                          ! Input matrix (ordered in shell formulation)
    real(rp),       intent(out)            :: &
         Kreorder(pdofs,pdofs)                      ! Output matrix (ordered as Alya:
    !                bottom nodes (1234) top nodes(5678)
    integer(ip)                            :: &
         j, k, m,                             &  ! counters
         icontBot, icontTop,                     &  !
         indexBot, indexTop,                     &  !
         icontBotCol, icontTopCol,               &  !
         indexBotCol, indexTopCol                   !

    real(rp)                              :: &  !
         auxMatrix(pdofs,pdofs),                 &  !
         auxVector(pdofs)                           !

    !
    ! Initialization
    !
    auxMatrix = 0.0_rp
    auxVector = 0.0_rp
    Kreorder = 0.0_rp

    icontBot = 4 !counter for the bottom  dof
    icontTop = 1 !counter for th top dof

    !Reordening the rows

    indexBot =1
    indexTop =13

    do m=1,4
       do j=1,pdofs
          auxMatrix(indexBot,j)   = Kmat(icontBot,j)
          auxMatrix(indexBot+1,j) = Kmat(icontBot+1,j)
          auxMatrix(indexBot+2,j) = Kmat(icontBot+2,j)

          auxMatrix(indexTop,j)   = Kmat(icontTop,j)
          auxMatrix(indexTop+1,j) = Kmat(icontTop+1,j)
          auxMatrix(indexTop+2,j) = Kmat(icontTop+2,j)
       end do

       icontBot  = icontBot + 6
       indexBot = indexBot + 3

       icontTop  = icontTop + 6
       indexTop = indexTop + 3
    end do

    icontBotCol =4
    indexBotCol =1

    icontTopCol =1
    indexTopCol =13

    do m =1,4
       do j=1,24
          do k=1,24
             auxVector(k) = auxMatrix(j,k)  !auxliar vector
          end do

          Kreorder(j,indexBotCol)   = auxVector(icontBotCol)
          Kreorder(j,indexBotCol+1) = auxVector(icontBotCol+1)
          Kreorder(j,indexBotCol+2) = auxVector(icontBotCol+2)

          Kreorder(j,indexTopCol)   = auxVector(icontTopCol)
          Kreorder(j,indexTopCol+1) = auxVector(icontTopCol+1)
          Kreorder(j,indexTopCol+2) = auxVector(icontTopCol+2)
       end do

       icontBotCol = icontBotCol + 6
       indexBotCol = indexBotCol + 3

       icontTopCol = icontTopCol + 6
       indexTopCol = indexTopCol + 3
    end do

  end subroutine K_element_reorder

  !------------------------------------------------------------------------------
  !> @author  Eva Casoni (eva.casoni@bsc.es)
  !> @date    17/09/2016
  !> @brief   Reorder rhs vector to be compatible with standard node numbering
  !> @details
  !-------------------------------------------------------------------------------

  subroutine ForceVectorReord(force,forceRe)

    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    real(rp),       intent(in)             :: &
         force(pdofs)                               ! Input vector  (ordered in shell formulation)
    real(rp),       intent(out)            :: &
         forceRe(pdofs)                             ! Output matrix (ordered as Alya:
    !                bottom nodes (1234) top nodes(5678)

    !initialization
    forceRe = 0.0_rp

    forceRe(1) = force(4)
    forceRe(2) = force(5)
    forceRe(3) = force(6)

    forceRe(4) = force(10)
    forceRe(5) = force(11)
    forceRe(6) = force(12)

    forceRe(7) = force(16)
    forceRe(8) = force(17)
    forceRe(9) = force(18)

    forceRe(10) = force(22)
    forceRe(11) = force(23)
    forceRe(12) = force(24)

    forceRe(13) = force(1)
    forceRe(14) = force(2)
    forceRe(15) = force(3)

    forceRe(16) = force(7)
    forceRe(17) = force(8)
    forceRe(18) = force(9)

    forceRe(19) = force(13)
    forceRe(20) = force(14)
    forceRe(21) = force(15)

    forceRe(22) = force(19)
    forceRe(23) = force(20)
    forceRe(24) = force(21)

  end subroutine ForceVectorReord

  subroutine nodes_dirichlet_reaction_forces(ielem, elrhs, elfrx)

    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    real(rp), intent(in)                   :: &
         elrhs(pdofs)
    real(rp), intent(out)                  :: &
         elfrx(pdime,pnode)
    integer(ip), intent(in)                :: &
         ielem
    ! ---------------------------------------------------------------------------
    integer(ip)                            :: &
         ipoin, inode, idime, idof
    !
    ! ---------------------------------------------------------------------------
    !
    ! Reaction force at Dirichlet BC
    !
    elfrx(:,:) = 0.0_rp
    do inode = 1, pnode
       ipoin = lnods(inode, ielem)
       do idime = 1, pdime
          if( kfl_fixno_sld(idime,ipoin) > 0 ) then
             idof = (inode-1)*pdime + idime
             elfrx(idime,inode) = elrhs(idof)
          end if
       end do
    end do

  end subroutine nodes_dirichlet_reaction_forces

  subroutine nodes_dirichlet_correction(ielem, elmat, elrhs)

    ! ---------------------------------------------------------------------------
    implicit none
    ! ---------------------------------------------------------------------------
    real(rp), intent(inout)                :: &
         elmat(pdofs,pdofs), elrhs(pdofs)
    integer(ip), intent(in)                :: &
         ielem
    ! ---------------------------------------------------------------------------
    integer(ip)                            :: &
         ipoin, idime, jdime, idof, jdof,        &
         inode, jnode
    real(rp)                               :: &
         adiag
    !
    ! ---------------------------------------------------------------------------
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

  end subroutine nodes_dirichlet_correction

  !-----------------------------------------------------------------------
  !>
  !> @author  Gerard Guillamet
  !> @date    March 2019
  !> @brief   Local axes rotation for rhs and amatr
  !> @details Local axes rotation for rhs and amatr
  !>
  !-----------------------------------------------------------------------

  subroutine sld_local_axes_rhs_and_amatr(itask, ielem, pnode, ndofn, elrhs, elmat)

    use def_domain,   only : lpoty
    use def_solidz,   only : kfl_fixno_sld, kfl_fixrs_sld
    use def_solidz,   only : jacrot_du_dq_sld, jacrot_dq_du_sld
    use mod_sld_csys, only : sld_csys_rotuni

    implicit none

    integer(ip), intent(in)    :: itask              !< What to do
    integer(ip), intent(in)    :: ielem              !< Element number
    integer(ip), intent(in)    :: pnode              !< Element number of nodes
    integer(ip), intent(in)    :: ndofn              !< Dimensions
    real(rp),    intent(inout) :: elrhs(pdofs)       !< Matrix
    real(rp),    intent(inout) :: elmat(pdofs,pdofs) !< rhs

    integer(ip)                :: inode, ipoin, ibopo, itott

    select case (itask)

    case(1_ip)
       !
       ! Explicit scheme
       !
       do inode = 1, pnode
          ipoin = lnods(inode,ielem)
          ibopo = lpoty(ipoin)
          itott = (inode - 1) * ndofn
          if ( ibopo > 0 ) then
             if ( kfl_fixno_sld(1,ipoin) == 3_ip .or. kfl_fixrs_sld(ipoin) /= 0_ip ) then

                call sld_csys_rotuni(1_ip,ndofn,ipoin,elrhs(itott+1:itott+ndofn))

             end if
          end if
       end do

    case(2_ip)
       !
       ! Implicit scheme
       !
       do inode = 1, pnode
          ipoin = lnods(inode,ielem)
          ibopo = lpoty(ipoin)
          itott = (inode - 1) * ndofn
          if ( ibopo > 0 ) then
             if ( kfl_fixno_sld(1,ipoin) == 3_ip .or. kfl_fixrs_sld(ipoin) /= 0_ip ) then

                call sld_rotsys(1_ip,inode,pnode,ndofn,ndofn*pnode,&
                     elmat,elrhs,jacrot_du_dq_sld(1,1,ipoin),jacrot_dq_du_sld(1,1,ipoin))

             end if
          end if
       end do

    end select

  end subroutine sld_local_axes_rhs_and_amatr

  !------------------------------------------------------------------------------
  !> @author  Eva Casoni (eva.casoni@bsc.es)
  !> @date    17/09/2016
  !> @brief   Function just to debug matrices and print in a friendly way
  !> @details
  !-------------------------------------------------------------------------------
  subroutine s81CurvTransf81(gkon,CC,C)

    implicit none

    real(rp), intent(in)          :: &
         gkon(3,3),CC(6,6)
    real(rp), intent(out)         :: &
         C(3,3,3,3)
    integer(ip)                   :: &
         i,j,k,m
    !
    !	Reassigning the 4th order tensor
    !    !material tensor in curvilinear coordinates coord
    !
    do i=1,3
       do j=1,3
          do k=1,3
             do m=1,3
                C(i,j,k,m) = &
                     CC(1,1)*gkon(1,i)*gkon(1,j)*gkon(1,k)*gkon(1,m) + &
                     CC(1,2)*gkon(1,i)*gkon(1,j)*gkon(1,k)*gkon(2,m) +  &
                     CC(1,3)*gkon(1,i)*gkon(1,j)*gkon(1,k)*gkon(3,m) + &
                     CC(1,2)*gkon(1,i)*gkon(1,j)*gkon(2,k)*gkon(1,m) + &
                     CC(1,4)*gkon(1,i)*gkon(1,j)*gkon(2,k)*gkon(2,m) + &
                     CC(1,5)*gkon(1,i)*gkon(1,j)*gkon(2,k)*gkon(3,m) +&
                     CC(1,3)*gkon(1,i)*gkon(1,j)*gkon(3,k)*gkon(1,m) +&
                     CC(1,5)*gkon(1,i)*gkon(1,j)*gkon(3,k)*gkon(2,m) +&
                     CC(1,6)*gkon(1,i)*gkon(1,j)*gkon(3,k)*gkon(3,m) +     &
                     CC(2,1)*gkon(1,i)*gkon(2,j)*gkon(1,k)*gkon(1,m) + &
                     CC(2,2)*gkon(1,i)*gkon(2,j)*gkon(1,k)*gkon(2,m) + &
                     CC(2,3)*gkon(1,i)*gkon(2,j)*gkon(1,k)*gkon(3,m) +&
                     CC(2,2)*gkon(1,i)*gkon(2,j)*gkon(2,k)*gkon(1,m) +&
                     CC(2,4)*gkon(1,i)*gkon(2,j)*gkon(2,k)*gkon(2,m) +&
                     CC(2,5)*gkon(1,i)*gkon(2,j)*gkon(2,k)*gkon(3,m) +&
                     CC(2,3)*gkon(1,i)*gkon(2,j)*gkon(3,k)*gkon(1,m) +&
                     CC(2,5)*gkon(1,i)*gkon(2,j)*gkon(3,k)*gkon(2,m) +&
                     CC(2,6)*gkon(1,i)*gkon(2,j)*gkon(3,k)*gkon(3,m) +   &
                     CC(3,1)*gkon(1,i)*gkon(3,j)*gkon(1,k)*gkon(1,m) +&
                     CC(3,2)*gkon(1,i)*gkon(3,j)*gkon(1,k)*gkon(2,m) +&
                     CC(3,3)*gkon(1,i)*gkon(3,j)*gkon(1,k)*gkon(3,m) +&
                     CC(3,2)*gkon(1,i)*gkon(3,j)*gkon(2,k)*gkon(1,m) +&
                     CC(3,4)*gkon(1,i)*gkon(3,j)*gkon(2,k)*gkon(2,m) +&
                     CC(3,5)*gkon(1,i)*gkon(3,j)*gkon(2,k)*gkon(3,m) +&
                     CC(3,3)*gkon(1,i)*gkon(3,j)*gkon(3,k)*gkon(1,m) +&
                     CC(3,5)*gkon(1,i)*gkon(3,j)*gkon(3,k)*gkon(2,m) +&
                     CC(3,6)*gkon(1,i)*gkon(3,j)*gkon(3,k)*gkon(3,m) +  &
                     CC(2,1)*gkon(2,i)*gkon(1,j)*gkon(1,k)*gkon(1,m) +&
                     CC(2,2)*gkon(2,i)*gkon(1,j)*gkon(1,k)*gkon(2,m) +&
                     CC(2,3)*gkon(2,i)*gkon(1,j)*gkon(1,k)*gkon(3,m)+&
                     CC(2,2)*gkon(2,i)*gkon(1,j)*gkon(2,k)*gkon(1,m)+&
                     CC(2,4)*gkon(2,i)*gkon(1,j)*gkon(2,k)*gkon(2,m)+&
                     CC(2,5)*gkon(2,i)*gkon(1,j)*gkon(2,k)*gkon(3,m)+&
                     CC(2,3)*gkon(2,i)*gkon(1,j)*gkon(3,k)*gkon(1,m)+&
                     CC(2,5)*gkon(2,i)*gkon(1,j)*gkon(3,k)*gkon(2,m)+&
                     CC(2,6)*gkon(2,i)*gkon(1,j)*gkon(3,k)*gkon(3,m)+ &
                     CC(4,1)*gkon(2,i)*gkon(2,j)*gkon(1,k)*gkon(1,m)+&
                     CC(4,2)*gkon(2,i)*gkon(2,j)*gkon(1,k)*gkon(2,m)+&
                     CC(4,3)*gkon(2,i)*gkon(2,j)*gkon(1,k)*gkon(3,m)+&
                     CC(4,2)*gkon(2,i)*gkon(2,j)*gkon(2,k)*gkon(1,m)+&
                     CC(4,4)*gkon(2,i)*gkon(2,j)*gkon(2,k)*gkon(2,m)+&
                     CC(4,5)*gkon(2,i)*gkon(2,j)*gkon(2,k)*gkon(3,m)+&
                     CC(4,3)*gkon(2,i)*gkon(2,j)*gkon(3,k)*gkon(1,m)+&
                     CC(4,5)*gkon(2,i)*gkon(2,j)*gkon(3,k)*gkon(2,m)+&
                     CC(4,6)*gkon(2,i)*gkon(2,j)*gkon(3,k)*gkon(3,m)+ &
                     CC(5,1)*gkon(2,i)*gkon(3,j)*gkon(1,k)*gkon(1,m)+&
                     CC(5,2)*gkon(2,i)*gkon(3,j)*gkon(1,k)*gkon(2,m)+&
                     CC(5,3)*gkon(2,i)*gkon(3,j)*gkon(1,k)*gkon(3,m)+&
                     CC(5,2)*gkon(2,i)*gkon(3,j)*gkon(2,k)*gkon(1,m)+&
                     CC(5,4)*gkon(2,i)*gkon(3,j)*gkon(2,k)*gkon(2,m)+&
                     CC(5,5)*gkon(2,i)*gkon(3,j)*gkon(2,k)*gkon(3,m)+&
                     CC(5,3)*gkon(2,i)*gkon(3,j)*gkon(3,k)*gkon(1,m)+&
                     CC(5,5)*gkon(2,i)*gkon(3,j)*gkon(3,k)*gkon(2,m)+&
                     CC(5,6)*gkon(2,i)*gkon(3,j)*gkon(3,k)*gkon(3,m)+ &
                     CC(3,1)*gkon(3,i)*gkon(1,j)*gkon(1,k)*gkon(1,m)+&
                     CC(3,2)*gkon(3,i)*gkon(1,j)*gkon(1,k)*gkon(2,m)+&
                     CC(3,3)*gkon(3,i)*gkon(1,j)*gkon(1,k)*gkon(3,m)+&
                     CC(3,2)*gkon(3,i)*gkon(1,j)*gkon(2,k)*gkon(1,m)+&
                     CC(3,4)*gkon(3,i)*gkon(1,j)*gkon(2,k)*gkon(2,m)+&
                     CC(3,5)*gkon(3,i)*gkon(1,j)*gkon(2,k)*gkon(3,m)+&
                     CC(3,3)*gkon(3,i)*gkon(1,j)*gkon(3,k)*gkon(1,m)+&
                     CC(3,5)*gkon(3,i)*gkon(1,j)*gkon(3,k)*gkon(2,m)+&
                     CC(3,6)*gkon(3,i)*gkon(1,j)*gkon(3,k)*gkon(3,m)+ &
                     CC(5,1)*gkon(3,i)*gkon(2,j)*gkon(1,k)*gkon(1,m)+&
                     CC(5,2)*gkon(3,i)*gkon(2,j)*gkon(1,k)*gkon(2,m)+&
                     CC(5,3)*gkon(3,i)*gkon(2,j)*gkon(1,k)*gkon(3,m)+&
                     CC(5,2)*gkon(3,i)*gkon(2,j)*gkon(2,k)*gkon(1,m)+&
                     CC(5,4)*gkon(3,i)*gkon(2,j)*gkon(2,k)*gkon(2,m)+&
                     CC(5,5)*gkon(3,i)*gkon(2,j)*gkon(2,k)*gkon(3,m)+&
                     CC(5,3)*gkon(3,i)*gkon(2,j)*gkon(3,k)*gkon(1,m)+&
                     CC(5,5)*gkon(3,i)*gkon(2,j)*gkon(3,k)*gkon(2,m)+&
                     CC(5,6)*gkon(3,i)*gkon(2,j)*gkon(3,k)*gkon(3,m)+ &
                     CC(6,1)*gkon(3,i)*gkon(3,j)*gkon(1,k)*gkon(1,m)+&
                     CC(6,2)*gkon(3,i)*gkon(3,j)*gkon(1,k)*gkon(2,m)+&
                     CC(6,3)*gkon(3,i)*gkon(3,j)*gkon(1,k)*gkon(3,m)+&
                     CC(6,2)*gkon(3,i)*gkon(3,j)*gkon(2,k)*gkon(1,m)+&
                     CC(6,4)*gkon(3,i)*gkon(3,j)*gkon(2,k)*gkon(2,m)+&
                     CC(6,5)*gkon(3,i)*gkon(3,j)*gkon(2,k)*gkon(3,m)+&
                     CC(6,3)*gkon(3,i)*gkon(3,j)*gkon(3,k)*gkon(1,m)+&
                     CC(6,5)*gkon(3,i)*gkon(3,j)*gkon(3,k)*gkon(2,m)+&
                     CC(6,6)*gkon(3,i)*gkon(3,j)*gkon(3,k)*gkon(3,m)
             end do
          end do
       end do
    end do

  end subroutine s81CurvTransf81

  subroutine s9Creordering(Caux,CC)

    implicit none

    real(rp), intent(in)          :: &
         Caux(6,6)
    real(rp), intent(out)         :: &
         CC(6,6)
    !
    ! rearranges the cosntitutive matrix
    ! from [11, 22, 33,  12, 23, 13] in the book Altenbach to
    !      [11, 12, 13,  22, 23, 33]
    CC(1,1) = Caux(1,1) !e11
    CC(1,2) = Caux(1,6) !gamma 12
    CC(1,3) = Caux(1,5) !gamma 13
    CC(1,4) = Caux(1,2) !e22
    CC(1,5) = Caux(1,4) !gamma 23
    CC(1,6) = Caux(1,3) !e33

    CC(2,1) = Caux(6,1) !e11
    CC(2,2) = Caux(6,6) !gamma 12
    CC(2,3) = Caux(6,5) !gamma 13
    CC(2,4) = Caux(6,2) !e22
    CC(2,5) = Caux(6,4) !gamma 23
    CC(2,6) = Caux(6,3) !e33

    CC(3,1) = Caux(5,1) !e11
    CC(3,2) = Caux(5,6) !gamma 12
    CC(3,3) = Caux(5,5) !gamma 13
    CC(3,4) = Caux(5,2) !e22
    CC(3,5) = Caux(5,4) !gamma 23
    CC(3,6) = Caux(5,3) !e33

    CC(4,1) = Caux(2,1) !e11
    CC(4,2) = Caux(2,6) !gamma 12
    CC(4,3) = Caux(2,5) !gamma 13
    CC(4,4) = Caux(2,2) !e22
    CC(4,5) = Caux(2,4) !gamma 23
    CC(4,6) = Caux(2,3) !e33

    CC(5,1) = Caux(4,1) !e11
    CC(5,2) = Caux(4,6) !gamma 12
    CC(5,3) = Caux(4,5) !gamma 13
    CC(5,4) = Caux(4,2) !e22
    CC(5,5) = Caux(4,4) !gamma 23
    CC(5,6) = Caux(4,3) !e33

    CC(6,1) = Caux(3,1) !e11
    CC(6,2) = Caux(3,6) !gamma 12
    CC(6,3) = Caux(3,5) !gamma 13
    CC(6,4) = Caux(3,2) !e22
    CC(6,5) = Caux(3,4) !gamma 23
    CC(6,6) = Caux(3,3) !e33

  end subroutine s9Creordering

end module mod_sld_contshell_elements
