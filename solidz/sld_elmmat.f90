!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!--------------------------------------------------------------------------
!> @author  <MRI> Matias Rivero (matias.rivero@bsc.es)
!> @author  <AQC> Adria Quintanas (adria.quintanas@udg.edu)
!> @author  <ECR> Eva Casoni (eva.casoni@bsc.es)
!> @author  <GGU> Gerard Guillamet (gerard.guillamet@bsc.es)
!>
!> @date    22/02/2017
!> @brief   Right Hand Side (RHS) and Matrix computations
!> @details Scheme:
!>          Status: (1) 10/02/2017 : AQC  : TODO : Rayleight damping revision
!>                  (2) 10/02/2017 : AQC  : TODO : Extract of BC correction and the
!>                                                 reaction force computation
!> @param[in] itask
!> @param[in] pgaus
!> @param[in] pmate
!> @param[in] pnode
!> @param[in] ielem
!> @param[in] lnods
!>
!--------------------------------------------------------------------------
subroutine sld_elmmat(itask,&
     ielem,pgaus,pmate,pnode,lnods,gpgdi,gppio,gppio_eps,gpstr,gpvol,gpvel,gpacc,gptmo,&
     gprat,gpsha,gpcar,gphea,gpdis,gpdds,gprestre,eldis,elcod,elrhs,elmat,elmuu,&
     elmaa,elfrx,elfint,elfext,elfine,eldip,elfco)
  use def_kintyp, only       :  ip,rp,lg
  use def_elmtyp, only       :  ELEXT
  use def_master, only       :  dtime,solve, ITER_K
  use def_domain, only       :  ndime,mnode,mgaus,lelch,lpoty
  use def_solidz, only       :  ndofn_sld,densi_sld,grnor_sld
  use def_solidz, only       :  gravi_sld,dampi_sld,kfl_fixno_sld
  use def_solidz, only       :  bodyf_sld,ncomp_sld,tifac_sld,accel_sld
  use def_solidz, only       :  kfl_timei_sld
  use def_solidz, only       :  kfl_exacs_sld
  use def_solidz, only       :  kfl_fixrs_sld, jacrot_du_dq_sld,jacrot_dq_du_sld,thick_sld
  use def_solidz, only       :  SLD_EXPLICIT_SCHEME, SLD_IMPLICIT_SCHEME, SLD_DYNAMIC_PROBLEM, SLD_STATIC_PROBLEM
  use mod_sld_bclocal, only  :  sld_bclocal_elmat_and_elrhs

  implicit none

  integer(ip), intent(in)    :: itask
  integer(ip), intent(in)    :: pgaus
  integer(ip), intent(in)    :: pmate
  integer(ip), intent(in)    :: pnode
  integer(ip), intent(in)    :: ielem
  integer(ip), intent(in)    :: lnods(pnode)                           ! List of nodes
  real(rp),    intent(in)    :: gppio(ndime,ndime,pgaus)               ! First Piola Kirchoff stress tensor at the gauss point
  real(rp),    intent(in)    :: gppio_eps(ndime,ndime,pgaus)           ! Perturbed First Piola Kirchoff stress tensor at the gauss point
  real(rp),    intent(in)    :: gpstr(ndime,ndime,pgaus,2)             ! Second Piola Kirchoff stress tensor at the gauss point
  real(rp),    intent(in)    :: gpvol(pgaus)                           ! Weight of gauss point
  real(rp),    intent(in)    :: gphea(pnode,mgaus)                     ! <GGU> Not used. Shifted heaviside (the dimensions is mgaus)
  real(rp),    intent(in)    :: gprestre(ndime,ndime,pgaus)            ! Residual stresses (from a prestressed field)
  real(rp),    intent(in)    :: gptmo(ndime,ndime,ndime,ndime,pgaus)   !
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)                     ! Shape functions at gauss points
  real(rp),    intent(in)    :: gpcar(ndime,mnode,mgaus)               ! Cartesian derivatives of the shape functions at gauss points
  real(rp),    intent(in)    :: gpgdi(ndime,ndime,mgaus)               ! Displacement gradient at gauss points
  real(rp),    intent(in)    :: gpvel(ndime,pgaus)                     ! Velocity at the gauss point
  real(rp),    intent(in)    :: gpacc(ndime,pgaus,3)                   ! Velocity at the gauss point (first ncomp_nsa)
  real(rp),    intent(in)    :: gprat(ndime,ndime,mgaus)               ! Rate of deformation Fdot at time n at gauss points
  real(rp),    intent(in)    :: gpdis(ndime,pgaus,*)                   ! Displacement at gauss points
  real(rp),    intent(in)    :: gpdds(ndime,ndime,ndime,ndime,mgaus)   ! Material tangent (dP/dF)
  real(rp),    intent(in)    :: eldis(ndime,pnode,ncomp_sld)           ! Elemental displacement
  real(rp),    intent(in)    :: elcod(ndime,pnode)                     ! Elemental node coordinates
  real(rp),    intent(in)    :: eldip(ndime,mnode)                     ! Inexact Newton
  real(rp),    intent(out)   :: elrhs(ndofn_sld*pnode)                 ! Elemental right hand side (RHS)
  real(rp),    intent(out)   :: elmat(ndofn_sld*pnode,ndofn_sld*pnode) ! Elemental matrix (AMAT)
  real(rp),    intent(out)   :: elmuu(pnode)                           !
  real(rp),    intent(out)   :: elmaa(pnode)                           !
  real(rp),    intent(out)   :: elfrx(ndofn_sld,pnode)                 ! Elemental reaction forces (RHS where Dirichlet BC are applied)
  real(rp),    intent(out)   :: elfco(ndofn_sld,pnode)                 ! Elemental reaction forces from contact
  real(rp),    intent(out)   :: elfint(ndofn_sld*pnode)                ! Internal force vector
  real(rp),    intent(out)   :: elfext(ndofn_sld*pnode)                ! External force vector
  real(rp),    intent(out)   :: elfine(ndofn_sld*pnode)                ! Inertial force vector
  real(rp)                   :: elrayl_alpha(ndofn_sld*pnode)          ! Rayleigh damping - alpha term (mass matrix) contribution
  real(rp)                   :: elfint_eps(ndofn_sld*pnode)            ! Auxiliar internal force vector for ninex

  integer(ip)                :: igaus,inode,idime,jdime,kdime,ldime,ibopo,ievat,jevat,pevat,itott
  integer(ip)                :: jnode,ipoin
  real(rp)                   :: elmas(pnode,pnode)
  real(rp)                   :: adiag,dtfa1,dtfa2,volux,betanewmark,dummr,grale
  real(rp)                   :: betar,alphr

  !
  ! Rayleigh coefficients (used in the implicit part)
  !
  alphr = (dampi_sld(1,pmate) / dtime)             ! Rayleigh damping: alpha/dt
  betar = (1.0_rp + dampi_sld(2,pmate) / dtime)    ! Rayleigh damping: 1 + beta/dt

  !
  ! Elemental gravity
  !
  grale = grnor_sld

  ! Elemental matrix/vector size
  pevat = ndofn_sld*pnode

  !
  ! Initialization
  !
  elrhs(:) = 0.0_rp
  elmuu(:) = 0.0_rp
  elmaa(:) = 0.0_rp
  elfrx(:,:) = 0.0_rp
  elfco(:,:) = 0.0_rp
  elmat(:,:) = 0.0_rp

  if( itask == SLD_EXPLICIT_SCHEME ) then

     !-------------------------------------------------------------------
     !
     ! EXPLICIT SCHEME
     !     RHS  = Fext - Fint
     !         According to the central differences scheme
     !
     !-------------------------------------------------------------------

     !
     ! Initialization
     !
     elfint(:) = 0.0_rp
     elfext(:) = 0.0_rp
     elfine(:) = 0.0_rp
     elrayl_alpha(:) = 0.0_rp

     !
     ! RHS: Internal & External foce contribution
     !     RHS = Fext - Fint
     !
     do igaus = 1,pgaus
        volux = gpvol(igaus)
        do inode = 1,pnode
           ipoin = lnods(inode)
           ievat = (inode-1) * ndofn_sld
           do idime = 1,ndime
              ievat = ievat + 1
              do jdime = 1,ndime
                 !
                 ! Internal force contribution
                 !    fint_kb = int_\Omega P_kL * dNb/dxL d\Omega
                 elfint(ievat) = elfint(ievat) + (gppio(idime,jdime,igaus) + gprestre(idime,jdime,igaus)) &
                      *gpcar(jdime,inode,igaus)*volux
              end do
              !
              ! External force contribution
              !    fext_kb = int_\Omega b_0*N d\Omega + int_\dOmega tract*N d\S
              !        where b_0 are bodyforces and  gravity
              !              tract are traction forces (Neumann)
              elfext(ievat) = elfext(ievat) + (gpsha(inode,igaus)*densi_sld(1,pmate)*volux) &
                   *(grale*gravi_sld(idime) + bodyf_sld(idime)) ! < AQC > In the elmmat_viejo the term bodyf_sld is negative
              !
              ! Rayleigh alpha contribution < AQC > Rayleigh damping revision
              if (dampi_sld(1,pmate) > 0.0_rp) then
                 elrayl_alpha(ievat) = elrayl_alpha(ievat) &
                      + dampi_sld(1,pmate)*gpvel(idime,igaus)*gpsha(inode,igaus)*densi_sld(1,pmate)*volux
              end if
           end do
        end do
     end do
     !
     ! Add contributions to the RHSID
     !
     elrhs(1:pevat) = elfext(1:pevat) - elfint(1:pevat) - elrayl_alpha(1:pevat)

     ! RHS stiffness matrix term of the rayleigh damping
     !    < AQC > Raylegih revision required
     if (dampi_sld(2,pmate) > 0.0_rp) then
        do igaus=1,pgaus
           volux= gpvol(igaus)
           do inode=1,pnode
              ievat= (inode-1)*ndime
              do idime=1,ndime
                 ievat=ievat+1
                 do jdime=1,ndime
                    do kdime=1,ndime
                       elrhs(ievat)= elrhs(ievat) &
                            + dampi_sld(2,pmate) &
                            * gpcar(jdime,inode,igaus) &
                            * gprat(kdime,jdime,igaus) &
                            * gpstr(jdime,kdime,igaus,1) &
                            * volux
                    end do
                 end do
              end do
           end do
        end do
     end if

     !
     ! Lumped mass matrix
     !
     do igaus = 1,pgaus
        volux = gpvol(igaus)
        do inode = 1,pnode
           do jnode = 1,pnode
              !
              ! Mass matrix
              !   M_ab = int_\Omega /rho * Na * Nb d\Omega
              elmuu(inode) = elmuu(inode) + densi_sld(1,pmate)*volux*gpsha(inode,igaus)*gpsha(jnode,igaus)
           end do
        end do
     end do

     !
     ! Exact solution
     !
     if( kfl_exacs_sld /= 0 ) then
        call sld_elmexa(&
             1_ip,pgaus,pnode,elcod,gpsha,gpcar,gpdds,gpvol,&
             eldis,gpdis,gpgdi,elrhs)
     end if
     !
     ! Boundary conditions: Dirichlet (essential B.Cs)
     !
     ! Local bases: Global --> Local
     !
     call sld_bclocal_elmat_and_elrhs(ndofn_sld,pnode,lnods,elrhs)
     !
     ! Save reaction forces
     !
     do inode = 1,pnode
        ipoin = lnods(inode)
        itott = (inode - 1) * ndofn_sld
        !
        ! Reactions (all)
        do idime = 1,ndime
           if( kfl_fixno_sld(idime,ipoin) > 0 ) then
              ievat = (inode-1)*ndofn_sld + idime
              elfrx(idime,inode) = elrhs(ievat)*thick_sld
           end if
        end do
        !
        ! Contact force
        if( kfl_fixno_sld(1,ipoin) == 3_ip ) then
           elfco(1:ndofn_sld,inode) = elrhs(itott+1:itott+ndofn_sld)
        end if
     end do

  else if( itask == SLD_IMPLICIT_SCHEME ) then

     !-------------------------------------------------------------------
     !
     ! IMPLICIT SCHEME:
     !     Static  : RHS  = Fext - Fint
     !               Amat = dFint/dU
     !     Dynamic : RHS  = Fext - Fint - Fine = Fext - Fint - M·a
     !               Amat = 1/(beta*dtime^2)*M + dFint/dU
     !               According to the Beta Newmark scheme
     !
     !-------------------------------------------------------------------
     !
     ! Initialization
     !
     elfint(:) = 0.0_rp
     elfint_eps(:) = 0.0_rp
     elfext(:) = 0.0_rp
     elfine(:) = 0.0_rp
     elrayl_alpha(:) = 0.0_rp

     !
     ! Newmark parameters
     !
     betanewmark = 0.0_rp
     if( tifac_sld(1) > 0.0_rp) betanewmark = 1.0_rp / (tifac_sld(1)*dtime*dtime)
     if( kfl_timei_sld == SLD_STATIC_PROBLEM ) betanewmark= 0.0_rp

     dtfa1 = 1.0_rp - tifac_sld(3)
     dtfa2 = tifac_sld(3)

     !
     ! RHS: Internal & External foce contribution
     !     RHS = Fext - Fint
     !
     do igaus = 1,pgaus
        volux = gpvol(igaus)
        do inode = 1,pnode
           ipoin = lnods(inode)
           ievat = (inode-1)*ndofn_sld
           do idime = 1,ndime
              ievat = ievat + 1
              do jdime = 1,ndime
                 !
                 ! Internal force contribution
                 !    fint_kb = int_\Omega P_kL * dNb/dxL d\Omega
                 elfint(ievat) = elfint(ievat) + (gppio(idime,jdime,igaus) + gprestre(idime,jdime,igaus)) &
                      *gpcar(jdime,inode,igaus)*volux
              end do
              !
              ! External force contribution
              !    fext_kb = int_\Omega b_0*N d\Omega + int_\dOmega tract*N d\S
              !        where b_0 are bodyforces and  gravity
              !              tract are traction forces (Neumann)
              elfext(ievat) = elfext(ievat) + (gpsha(inode,igaus)*densi_sld(1,pmate)*volux) &
                   *(grale*gravi_sld(idime) - bodyf_sld(idime))
              !
              ! Rayleigh alpha contribution
              !   < AQC > Alpha rayleigh shouldn't be proportional to the mass matrix?
              if (dampi_sld(1,pmate) > 0.0_rp) then
                 elrayl_alpha(ievat) = elrayl_alpha(ievat) + dampi_sld(1,pmate)*gpvel(idime,igaus)*gpsha(inode,igaus)*densi_sld(1,pmate)*volux
              end if
           end do
        end do
     end do

     ! < MR > Look for a way to optimize this...
     elrhs(1:pevat) = elfext(1:pevat) - elfint(1:pevat) - elrayl_alpha(1:pevat)

     !
     ! ELMAT: Internal force linearitzation contribution
     !     A = d Fint / d U
     do igaus = 1,pgaus
        volux= gpvol(igaus)
        do inode = 1,pnode
           do idime = 1,ndime
              ievat = ( inode - 1 ) * ndofn_sld + idime
              do jnode = 1,pnode
                 do kdime = 1,ndime
                    jevat = ( jnode - 1 ) * ndofn_sld + kdime
                    do jdime = 1,ndime
                       do ldime = 1,ndime
                          !
                          ! Internal force linearitzation and beta rayleigh damping contribution
                          !    K_iakb = int_\Omega (dP_iJ / dF_kL * dNa/dxJ * dNb/dxL * (1 + beta_rayleigh / dt) d\Omega
                          !    < AQC > The rayleigh contribution should be verified.
                          elmat(ievat,jevat) = elmat(ievat,jevat) + (gptmo(idime,jdime,kdime,ldime,igaus) &
                               *gpcar(jdime,inode,igaus)*gpcar(ldime,jnode,igaus)*betar)*volux
                          !elmat(ievat,jevat) = elmat(ievat,jevat) + (gptmo(idime,jdime,kdime,ldime,igaus) &
                          !     *gpcar(jdime,inode,igaus)*gpcar(ldime,jnode,igaus))*volux
                       end do
                    end do
                    !
                    ! Alpha rayleigh damping  contribution
                    !    K_iakb += int_\Omega (alpha_rayleigh / dt) Na * Nb  ) d\Omega
                    !    < AQC > The rayleigh contribution should be verified.
                    !elmat(ievat,jevat) = elmat(ievat,jevat) + densi_sld(1,pmate)*alphr*volux &
                    !     *gpsha(inode,igaus)*gpsha(jnode,igaus)
                 end do
              end do
           end do
        end do
     end do

     !
     ! RHS & ELMAT: Inertial terms contribution
     !
     if( kfl_timei_sld == SLD_DYNAMIC_PROBLEM ) then
        !
        ! Initialization
        !
        elmas(:,:) = 0.0_rp

        !
        ! Mass matrices
        !
        do igaus = 1,pgaus
           volux = gpvol(igaus)
           do inode = 1,pnode
              do jnode = 1,pnode
                 !
                 ! Consistent mass matrix
                 !   M_ab = int_\Omega /rho * Na * Nb d\Omega
                 elmas(inode,jnode) = elmas(inode,jnode) + volux*gpsha(inode,igaus)*gpsha(jnode,igaus)
                 !
                 ! Lumped mass matrix
                 !
                 elmuu(inode) = elmuu(inode) + densi_sld(1,pmate)*volux*gpsha(inode,igaus)*gpsha(jnode,igaus)
              end do
           end do
        end do

        !
        ! RHS: Inertial forces contribution
        !    RHS = Fint - Fext - Ma
        !
        do idime = 1,ndime
           do inode = 1,pnode
              ievat = (inode-1)*ndofn_sld + idime
              do jnode = 1,pnode
                 ipoin = lnods(jnode)
                 !
                 ! Inertial forces
                 !    f_kb += M_ba * accel_ak
                 elfine(ievat) = elfine(ievat) +  densi_sld(1,pmate)*elmas(inode,jnode)*accel_sld(idime,ipoin,1)
              end do
           end do
        end do

        !
        ! RHS inertial contribution
        !
        do ievat = 1,(ndofn_sld*pnode)
           elrhs(ievat) = elrhs(ievat) - elfine(ievat)
        end do

        !
        ! ELMAT: Inertial contribution
        !    A = 1/(beta*dtime^2)*M + dFint/dU
        do idime = 1,ndime
           do inode = 1,pnode
              ievat = ( inode - 1 )*ndofn_sld + idime
              do jnode = 1,pnode
                 jevat = ( jnode - 1 ) *ndofn_sld + idime
                 !
                 ! Inertial term contribution
                 !    K_iakb += 1/(beta*dtime^2)*M
                 elmat(ievat,jevat) = elmat(ievat,jevat) + elmas(inode,jnode)*(densi_sld(1,pmate)*betanewmark)
              end do
           end do
        end do

     end if ! ----- | END : DYNAMIC | -----

     !
     ! Exact solution
     !
     if( kfl_exacs_sld /= 0 ) then
        call sld_elmexa(1_ip,pgaus,pnode,elcod,gpsha,gpcar,gpdds,gpvol,eldis,gpdis(1,1,ITER_K),gpgdi,elrhs)
     end if

     !
     ! Extension elements (Dodeme)
     !
     if( lelch(ielem) == ELEXT ) then
        call elmext(4_ip,ndime,pnode,dummr,dummr,dummr,dummr,elmat,dummr,elrhs,dummr)
     end if
     !
     ! Boundary conditions: Dirichlet (essential B.Cs)
     !
     ! Local bases
     !
     do inode = 1,pnode
        ipoin = lnods(inode)
        ibopo = lpoty(ipoin)
        ievat = (inode-1)*ndofn_sld
        if (ibopo > 0) then
           if( kfl_fixno_sld(1,ipoin) == 2_ip .or. &
               kfl_fixno_sld(1,ipoin) == 3_ip .or. &
               kfl_fixrs_sld(ipoin)   /= 0_ip ) then
              call sld_rotsys(1_ip, &
                   inode,pnode,ndofn_sld,pevat,&
                   elmat,elrhs,&
                   jacrot_du_dq_sld(1,1,ipoin),jacrot_dq_du_sld(1,1,ipoin))
           end if
        end if
     end do
     !
     ! Save reaction forces
     !
     do inode = 1,pnode
        ipoin = lnods(inode)
        itott = (inode - 1) * ndofn_sld
        !
        ! Reactions (all)
        do idime = 1,ndime
           if( kfl_fixno_sld(idime,ipoin) > 0 ) then
              ievat = (inode-1)*ndofn_sld + idime
              elfrx(idime,inode) = elrhs(ievat)*thick_sld
           end if
        end do
        !
        ! Contact force
        if( kfl_fixno_sld(1,ipoin) == 3_ip ) then
           elfco(1:ndofn_sld,inode) = elrhs(itott+1:itott+ndofn_sld)
        end if
     end do
     !
     ! All, when no iffix
     !
     if( solve(1) % kfl_iffix == 0 ) then
        do inode = 1,pnode
           ipoin = lnods(inode)
           do idime = 1,ndime
              if( kfl_fixno_sld(idime,ipoin) > 0 ) then
                 ievat = (inode-1)*ndofn_sld + idime
                 adiag = elmat(ievat,ievat)
                 do jnode = 1,pnode
                    do jdime = 1,ndime
                       jevat = (jnode-1)*ndofn_sld + jdime
                       elmat(ievat,jevat) = 0.0_rp
                       elmat(jevat,ievat) = 0.0_rp
                    end do
                 end do
                 elmat(ievat,ievat) = adiag
                 elrhs(ievat)       = 0.0_rp ! We work with increments
              end if
           end do
        end do
     end if
  end if

end subroutine sld_elmmat
