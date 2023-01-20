!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup NastinMatrixAssembly
!> @{
!> @file    nsi_elmope_omp.f90
!> @author  Guillaume Houzeaux
!> @brief   Navier-Stokes system element assembly and other element
!>          calculations
!> @details Elemental operations, according to ITASK:
!>
!>          \verbatim
!>
!>          1 ........ Element calculations and assembly of global system:
!>                     A <= A^(e): matrix ........................... AMATR
!>                     b <= b^(e): RHS .............................. RHSID
!>                     Q <= Q^(e): Schur complement precond ......... LAPLA_NSI
!>                     M <= M^(e): Consistent Mass Matrix ........... CMAMA_NSI
!>
!>                     Resulting system is A u* = b. For the pressure Schur
!>                     complement-based solver, the following blocks are
!>                     assembled:
!>
!>                     +-        +  +- -+     +-  -+
!>                     | Auu Aup |  | u |     | bu |
!>                     |         |  |   |  =  |    |
!>                     | Apu App |  | p |     | bp |
!>                     +-       -+  +- -+     +-  -+
!>
!>          4 ........ Compute SGS .................................. VESGS_NSI
!>          5 ........ Assemble limiter to be L2 projected .......... RHSID
!>          6 ........ Assemble pressure Schur complement precond ... LAPLA_NSI
!>          10->19 ... Assemble properties to be L2 projected ....... PROPE_NSI
!>          30     ... Assemble Coordinate derivatives of R and F ....... 
!>
!>          \endverbatim
!>
!> @}
!------------------------------------------------------------------------
subroutine nsi_elmope_omp(itask)
  use def_elmtyp,            only : ELCUT
  use def_kintyp,            only : ip,rp
  use def_master,            only : amatr,rhsid,cutim,kfl_timco,zeror
  use def_master,            only : tesgs,lumma,vesgs,rgene
  use def_master,            only : IEMPTY
  use def_kermod,            only : kfl_kxmod_ker
  use def_kermod,            only : kfl_cos_opt,kfl_adj_prob,costf
  use def_kermod,            only : kfl_duatss
  use def_domain,            only : ltype,nnode
  use def_domain,            only : ngaus,llapl,lorde,ltopo
  use def_domain,            only : ngaus,ndime,lnods,nmate
  use def_domain,            only : lelch,elmar,mnode,mgaus
  use def_domain,            only : lmate,hnatu,ntens
  use mod_elmgeo,            only : elmgeo_bubble
  use mod_elmgeo,            only : FREE_SURFACE_BUBBLE
  use def_elmtyp,            only : ELEXT
  use def_nastin,            only : ncomp_nsi,dtmax_nsi
  use def_nastin,            only : kfl_stabi_nsi,nbdfp_nsi
  use def_nastin,            only : kfl_sgste_nsi,kfl_advec_nsi
  use def_nastin,            only : tamin_nsi,tamax_nsi
  use def_nastin,            only : dtinv_nsi,safet_nsi, safma_nsi
  use def_nastin,            only : dtcri_nsi,saflo_nsi
  use def_nastin,            only : dtsgs_nsi,kfl_stead_nsi
  use def_nastin,            only : kfl_timei_nsi,NSI_MONOLITHIC
  use def_nastin,            only : NSI_SCHUR_COMPLEMENT
  use def_nastin,            only : NSI_FRACTIONAL_STEP
  use def_nastin,            only : kfl_grvir_nsi
  use def_nastin,            only : prope_nsi
  use def_nastin,            only : lapla_nsi,kfl_sgscp_nsi
  use def_nastin,            only : cmama_nsi,kfl_corre_nsi
  use def_nastin,            only : kfl_force_nsi,kfl_shock_nsi
  use def_nastin,            only : kfl_matdi_nsi,poauu_nsi
  use def_nastin,            only : poaup_nsi,poapp_nsi
  use def_nastin,            only : poapu_nsi,ndbgs_nsi
  use def_nastin,            only : kfl_ellen_nsi
  use def_nastin,            only : kfl_ellsh_nsi,lforc_material_nsi
  use def_nastin,            only : xforc_material_nsi
  use def_nastin,            only : resis_nsi,itsta_nsi
  use def_nastin,            only : rmsgs_nsi,resgs_nsi
  use def_kermod,            only : sens_mesh
  use def_nastin,            only : kfl_bubbl_nsi
  use def_nastin,            only : dt_rho_nsi
  use def_nastin,            only : mass_rho_nsi
  use def_nastin,            only : tau_nsi
  use mod_nsi_subgrid_scale, only : nsi_subgrid_scale_gather
  use mod_nsi_subgrid_scale, only : nsi_subgrid_scale_residual_and_update
  use mod_nsi_element_operations, only : nsi_rhodt_rhotau_nu
  use mod_nsi_assembly_global_system, only : nsi_assembly_fractional_step_scalar
  use mod_nsi_bubble,        only : nsi_bubble_assembly
  use mod_ker_proper,        only : ker_proper
  use mod_parall,            only : par_omp_num_colors
  use mod_parall,            only : par_omp_ia_colors
  use mod_parall,            only : par_omp_ja_colors
#if !defined(OMPSS) && defined(_OPENMP)
  use def_kermod,            only : kfl_prope
  use def_nastin,            only : gradv_nsi
  use def_master,            only : current_zone
  use def_nastin,            only : bubble_aqq_nsi,bubble_aqu_nsi
  use def_nastin,            only : bubble_aqp_nsi,bubble_bq_nsi
  use def_nastin,            only : kfl_predi_nsi
  use mod_parall,            only : par_omp_nelem_chunk
#endif
#if !defined(OMPSS) && defined(_OPENMP) || defined(NO_COLORING)
  use def_domain,            only : nelem
#endif

  implicit none

  integer(ip), intent(in) :: itask                     !< What to do
  !
  ! Element matrices and vectors (stiffness and preconditioner)
  !
  real(rp)    :: elauu(mnode*ndime,mnode*ndime)        ! Auu
  real(rp)    :: elaup(mnode*ndime,mnode)              ! Aup
  real(rp)    :: elapp(mnode,mnode)                    ! App
  real(rp)    :: elapu(mnode,mnode*ndime)              ! Apu
  real(rp)    :: elrbu(ndime,mnode)                    ! bu
  real(rp)    :: elrbp(mnode)                          ! bp
  real(rp)    :: elrhs(6*mnode)                        ! Generic RHS
  real(rp)    :: elmap(mnode,mnode)                    ! Q
  real(rp)    :: elcmm(mnode*ndime,mnode*ndime)        ! Consistent mass matrix
  !
  ! Bubble matrices
  !
  real(rp)    :: elauq(mnode*ndime,1)
  real(rp)    :: elapq(mnode,1)
  real(rp)    :: elaqu(1,mnode*ndime)
  real(rp)    :: elaqp(1,mnode)
  real(rp)    :: elaqq(1,1)
  real(rp)    :: elrbq(1)
  !
  ! Gather
  !
  real(rp)    :: elvel(ndime,mnode,ncomp_nsi)          ! u
  real(rp)    :: elpre(mnode,ncomp_nsi-1)              ! p
  real(rp)    :: elfle(mnode)                          ! phi
  real(rp)    :: elcod(ndime,mnode)                    ! x
  real(rp)    :: elvep(ndime,mnode)                    ! Pi(momentum)
  real(rp)    :: elprp(mnode)                          ! Pi(div(u))
  real(rp)    :: elgrp(ndime,mnode)                    ! Pi(grad(p))
  real(rp)    :: eltem(mnode,ncomp_nsi)                ! T
  real(rp)    :: elwmean(mnode,ncomp_nsi)              ! W mean
  real(rp)    :: elmsh(ndime,mnode)                    ! u mesh
  real(rp)    :: elnor(ndime,mnode)                    ! Normal to the Level Set interface
  real(rp)    :: elcur(mnode)                          ! Level Set interface curvature
  real(rp)    :: ellum(mnode)                          ! Lumped mass matrix
  real(rp)    :: elbub(mnode)                          ! Element bubble
  real(rp)    :: elrhodt(mnode)                        ! Element bubble
  real(rp)    :: elmurho(mnode)                        ! Element bubble
  real(rp)    :: elrhotau(mnode)                       ! Element bubble
  !
  ! Indices and dimensions
  !
  integer(ip) :: ielem,inode
  integer(ip) :: pnode,pgaus,pevat,kelem,dummi
  integer(ip) :: pelty,plapl,porde,pmate,ptopo
  integer(ip) :: ipoin,icolo,pbubl
  !
  ! Gauss point values
  !
  real(rp)    :: gpsha(mnode,mgaus)                    ! N
  real(rp)    :: gpder(ndime,mnode,mgaus)              ! dN/dsi
  real(rp)    :: gpcar(ndime,mnode,mgaus)              ! dN/dxi
  real(rp)    :: gphes(ntens,mnode,mgaus)              ! d2N/dxidxj
  real(rp)    :: gplap(mnode,mgaus)                    ! Lapl(N)
  real(rp)    :: gpvol(mgaus)                          ! w*|J|, |J|
  real(rp)    :: gpvis(mgaus)                          ! Viscosity
  real(rp)    :: gpgvi(ndime,mgaus)                    ! Viscosity gradients
  real(rp)    :: gpnut(mgaus)                          ! mut
  real(rp)    :: grvis(ndime,mgaus)                    ! grad(mut)
  real(rp)    :: gppor(mgaus)                          ! Porosity
  real(rp)    :: gpden(mgaus)                          ! Density
  real(rp)    :: gpfle(mgaus)                          ! Level set function
  real(rp)    :: gpst1(mgaus)                          ! tau1
  real(rp)    :: gpst2(mgaus)                          ! tau2
  real(rp)    :: gpsp1(mgaus)                          ! tau1'
  real(rp)    :: gpsp2(mgaus)                          ! tau2'
  real(rp)    :: gptt1(mgaus)                          ! tau1'/tau1
  real(rp)    :: gptt2(mgaus)                          ! tau2'/tau2
  real(rp)    :: gpadv(ndime,mgaus)                    ! u+u'
  real(rp)    :: gprhs(ndime,mgaus)                    ! RHS
  real(rp)    :: gprhs_sgs(ndime,mgaus)                ! RHS due to subscales
  real(rp)    :: gprhc(mgaus)                          ! RHS for the continuity equation (Low Mach+penalization)
  real(rp)    :: gprh2(mgaus)                          ! RHS for the residual of continuity equation (Low Mach)
  real(rp)    :: gpsgs(ndime,mgaus,2)                  ! u'
  real(rp)    :: gpsgt(mgaus)                          ! T'
  real(rp)    :: gppre(mgaus,ncomp_nsi-1)              ! p
  real(rp)    :: gpvel(ndime,mgaus,ncomp_nsi-1)        ! u
  real(rp)    :: gpgve(ndime,ndime,mgaus)              ! grad(u)
  real(rp)    :: gpgpr(ndime,mgaus,2)                  ! grad(p)
  real(rp)    :: gpgde(ndime,mgaus)                    ! grad(den)
  real(rp)    :: gpdde(mnode,mgaus)                    ! Density derivatives w.r.t nodal temperature
  real(rp)    :: gpgdd(ndime,mnode,mgaus)              ! Density derivatives w.r.t nodal temperature and coordinates
  real(rp)    :: gpdvi(mnode,mgaus)                    ! Viscosity derivatives w.r.t nodal temperature
  real(rp)    :: gpgdv(ndime,mnode,mgaus)              ! Viscosity derivatives w.r.t nodal temperature and coordinates
  real(rp)    :: gptem(mgaus,ncomp_nsi)                ! T
  real(rp)    :: gpsgi(ndime,mgaus)                    ! SGS (work array)
  real(rp)    :: gpvep(ndime,mgaus)                    ! -tau1'*R(u) or tau1*rho*(a.grad)u
  real(rp)    :: gpprp(mgaus)                          ! tau2*div(u)
  real(rp)    :: gpgrp(ndime,mgaus)                    ! tau1'*( grad(p) - rho*f )
#if !defined(OMPSS) && defined(_OPENMP)
  real(rp)    :: gpfli(mgaus)                          ! phy_hyd
#endif
  real(rp)    :: gphyd(mgaus)                          ! rho_hyd
  real(rp)    :: gpmsh(ndime,mgaus)                    ! u_mesh
  real(rp)    :: gpnor(ndime,mgaus)                    ! LS normal
  real(rp)    :: gpcur(mgaus)                          ! LS curvature
  !
  ! Enrichement
  !
  real(rp)    :: gpsha_bub(mgaus)                      ! Ne
  real(rp)    :: gpcar_bub(ndime,mgaus)                ! dNe/dxi
  !
  ! Element characteristics (to compute tau1 and tau2)
  !
  real(rp)    :: tragl(ndime,ndime)
  real(rp)    :: chave(ndime,2)
  real(rp)    :: chale(2)
  real(rp)    :: hleng(ndime)
  real(rp)    :: dummr(2),dtcri
  real(rp)    :: dtinv_loc,dtsgs_loc
  !
  ! Perturbation and residuals
  !
  real(rp)    :: rmomu(mnode,mgaus)                    ! Residual velocity in momentum
  real(rp)    :: rmom2(ndime,ndime,mnode,mgaus)        ! Residual velocity in momentum
  real(rp)    :: rcont(ndime,mnode,mgaus)              ! Residual velocity in continuity
  real(rp)    :: wgrgr(mnode,mnode,mgaus)              ! grad(Ni).grad(Nj)
  real(rp)    :: wgrvi(mnode,mgaus)                    ! grad(mu).grad(Ni)
  real(rp)    :: agrau(mnode,mgaus)                    ! a.grad(Ni)
  real(rp)    :: p1vec(mnode,mgaus)                    ! Test funct. velocity in momentum
  real(rp)    :: p1ve2(ndime,ndime,mnode,mgaus)        ! Test funct. velocity in momentum
  real(rp)    :: p2vec(ndime,mnode,mgaus)              ! Test funct. velocity in continuity
  real(rp)    :: p2sca(mnode,mgaus)                    ! Test function pressure in continuity
  !
  ! Exact linearization and adjoint
  !
  real(rp)    :: gpstrm(ndime,mgaus)                   ! Momentum strong residual
  real(rp)    :: gpstrc(mgaus)                         ! Continuity strong residual
  real(rp)    :: elaut(mnode*ndime,mnode)
  real(rp)    :: elapt(mnode,mnode)
  real(rp)    :: elaputrans(mnode*ndime,mnode)         ! Transpose of the Apu
  real(rp)    :: elauptrans(mnode,mnode*ndime)         ! Transpose of the Aup
  real(rp)    :: dgpmut_dvel(ndime,mnode,mgaus)        ! Turbulence viscosity derivatives w.r.t. nodal velocity
  real(rp)    :: dgpmut_dtur(1,mnode,mgaus)            ! Turbulence viscosity derivatives w.r.t. turbulence unk
  real(rp)    :: gpvol_der(ndime,mnode,mgaus)          ! Derivative of w*|J|, |J| w.r.t. coordinate nodes
  real(rp)    :: gpcar_der(ndime,mnode,ndime,mnode,mgaus) ! Derivative of gpcar w.r.t. coordinate nodes
  real(rp)    :: hleng_der(ndime,mnode,ndime)          ! Derivative of hleng w.r.t. coordinate nodes
  real(rp)    :: chale_der(ndime,mnode,2)              ! Derivative of chale w.r.t. coordinate nodes
  real(rp)    :: densi(mgaus,nbdfp_nsi)
  real(rp)     :: timea,timeb, facto
#if !defined(OMPSS) && defined(_OPENMP)
  integer(4)   :: idime,igaus
#endif
  !
  ! Initialization
  !
  call cputim(timea)
  if( IEMPTY ) return
  
  dtmax_nsi = -1.0_rp

#ifndef NO_COLORING
  colors: do icolo = 1,par_omp_num_colors
#endif

#ifndef OMPSS
     !
     ! OpenMP declarations
     !$OMP  PARALLEL DO                                                            &
     !$OMP  SCHEDULE     ( DYNAMIC , par_omp_nelem_chunk )                         &
     !$OMP  DEFAULT      ( NONE )                                                  &
     !$OMP  PRIVATE      ( elauu, elaup, elapp, elapu, elrbu, elrhs, elrbp, elmap, &
     !$OMP                 elvel, elpre, elfle, elcod, elvep, elprp, elgrp, gpder, &
     !$OMP                 eltem, ielem, kelem, igaus, pnode, pgaus, pevat,        &
     !$OMP                 pelty, plapl, porde, pmate, ptopo, gpcar, gphes, gplap, &
     !$OMP                 gpvol, gpnut, grvis, gpfle, gpst1, gpst2, gpsp1, gpsp2, &
     !$OMP                 gptt1, gptt2, gpadv, gpsgt, gppre, gpvel, gpgve, gpgpr, &
     !$OMP                 gptem, gpsgi, gpvep, gpprp, gpgrp, gpfli, gphyd, tragl, &
     !$OMP                 chave, chale, hleng, dummr, rmomu, rcont, wgrgr, wgrvi, &
     !$OMP                 p1vec, p2vec, p2sca, elmsh, elnor, elcur, dtcri, gpsha, &
     !$OMP                 idime, gprhc, gpnor, agrau, dummi, gpcur, ellum, inode, &
     !$OMP                 ipoin, gprhs_sgs, elwmean, gpgde , gpden, gpvis, gppor, &
     !$OMP                 gprhs, gprh2, gpgvi, gpsgs, gpmsh, rmom2, p1ve2,        &
     !$OMP                 gpdde,gpgdd,gpdvi, gpgdv, gpsha_bub, gpcar_bub, pbubl,  &
     !$OMP                 dgpmut_dvel,dgpmut_dtur, elrhodt,elrhotau,elmurho,      &
     !$OMP                 gpvol_der,gpcar_der,hleng_der,chale_der,densi,          &
     !$OMP                 gpstrm, gpstrc,elaut, elapt, elaputrans, elauptrans,    &
     !$OMP                 dtinv_loc, dtsgs_loc ,elcmm, elauq, elapq, elaqu,       &
     !$OMP                 elaqp, elaqq, elrbq, elbub, facto )                     &
     !$OMP  SHARED       ( kfl_ellen_nsi, elmar, hnatu, kfl_advec_nsi, lorde,      &
     !$OMP                 itask, llapl, lmate, current_zone, lumma,               &
     !$OMP                 ltype, nnode,  ngaus, kfl_sgste_nsi, tesgs, lnods,      &
     !$OMP                 nmate, ltopo, lelch, kfl_bubbl_nsi, nelem,              &
     !$OMP                 poauu_nsi, poaup_nsi, poapu_nsi,         &
     !$OMP                 poapp_nsi, ndbgs_nsi, kfl_ellsh_nsi, kfl_predi_nsi,     &
     !$OMP                 safet_nsi, saflo_nsi, safma_nsi, kfl_stead_nsi, kfl_timei_nsi,     &
     !$OMP                 kfl_kxmod_ker, kfl_grvir_nsi, mgaus, cutim,             &
     !$OMP                 lforc_material_nsi,xforc_material_nsi, kfl_matdi_nsi,   &
     !$OMP                 kfl_stabi_nsi, kfl_timco, dtcri_nsi, kfl_prope,         &
     !$OMP                 kfl_force_nsi, kfl_sgscp_nsi, kfl_shock_nsi,            &
     !$OMP                 par_omp_ja_colors, par_omp_ia_colors, icolo, mnode,     &
     !$OMP                 par_omp_num_colors, NSI_MONOLITHIC, kfl_duatss,         & 
     !$OMP                 NSI_SCHUR_COMPLEMENT,dtinv_nsi,dtsgs_nsi,kfl_corre_nsi, &
     !$OMP                 kfl_cos_opt,kfl_adj_prob,costf,sens_mesh,               &
     !$OMP                 NSI_FRACTIONAL_STEP,par_omp_nelem_chunk,                &
#ifndef NDIMEPAR
     !$OMP                 ndime,                                                  &
#endif
     !$OMP                 vesgs, prope_nsi, gradv_nsi, rhsid, amatr, lapla_nsi,   &
     !$OMP                 cmama_nsi, bubble_aqq_nsi,bubble_aqu_nsi,               &
     !$OMP                 bubble_aqp_nsi,bubble_bq_nsi,dt_rho_nsi,tau_nsi,        &
     !$OMP                 mass_rho_nsi)                                             &
     !$OMP REDUCTION       ( +:resis_nsi,itsta_nsi,resgs_nsi )                     &
     !$OMP REDUCTION       ( MIN:tamin_nsi )                                       &
     !$OMP REDUCTION       ( MAX:rmsgs_nsi,tamax_nsi,dtmax_nsi )
     !
     ! Loop over elements
     !
#endif

#ifndef NO_COLORING
     elements: do kelem = par_omp_ia_colors(icolo),par_omp_ia_colors(icolo+1)-1
        ielem = par_omp_ja_colors(kelem)
#else
     elements: do ielem = 1,nelem
#endif
        pelty = ltype(ielem)
        !
        ! Element properties and dimensions
        !
        if( pelty > 0 ) then
           pnode = nnode(pelty)
           pgaus = ngaus(pelty)
           plapl = llapl(pelty)
           porde = lorde(pelty)
           ptopo = ltopo(pelty)
           pevat = ndime * pnode
           pbubl = 0
           if( kfl_stabi_nsi == 2 ) plapl = 0
           !
           ! Check if element is a solid
           !
           if( nmate > 1 ) then
              pmate = lmate(ielem)
           else
              pmate = 1
           end if
           !
           ! Initializations
           !
           gpden     =  0.0_rp
           gpvis     =  0.0_rp
           gppor     =  0.0_rp
           gprhc     =  0.0_rp
           gprh2     =  0.0_rp
           gpgvi     =  0.0_rp
           gpsgs     =  0.0_rp
           gpmsh     =  0.0_rp
           rmom2     =  0.0_rp
           p1ve2     =  0.0_rp
           gpgde     =  0.0_rp
           gpdde     =  0.0_rp
           gpgdd     =  0.0_rp
           gpdvi     =  0.0_rp
           gpgdv     =  0.0_rp
           gpnut     =  0.0_rp
           dgpmut_dvel = 0.0_rp
           dgpmut_dtur = 0.0_rp
           grvis     = 0.0_rp
           dtinv_loc = dtinv_nsi
           dtsgs_loc = dtsgs_nsi
           !
           ! Initializations of the subgrid scales
           !
           call nsi_subgrid_scale_gather(ndime,pgaus,ielem,vesgs,gpsgs)
           if( kfl_sgste_nsi == 1 ) then
              gpsgt(1:pgaus) = tesgs(ielem) % a(1,1:pgaus,1)
           else
              gpsgt(1:pgaus) = 0.0_rp
           end if
           
           !
           ! Gather
           !
           call nsi_elmga3(&
                pnode,ielem,lnods(1,ielem),elcod,elpre,elvel,elfle,&
                elvep,elprp,elgrp,eltem,elmsh,elnor,elcur,         &
                elwmean,elbub)
           !
           ! HLENG and TRAGL at center of gravity
           !
           call elmlen(&
                ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                hnatu(pelty),hleng)
           !
           ! Derivatives of HLENG w.r.t. nodal coordinates
           !
           if(itask == 30) & 
             call elmlen_der(&
                ndime,pnode,elmar(pelty)%dercg,elcod,&
                hnatu(pelty),hleng_der)
           
           !
           ! Compute the characteristic length: CHALE
           !
           call elmchl(&
                tragl,hleng,elcod,elvel,chave,chale,pelty,pnode,&
                porde,hnatu(pelty),kfl_advec_nsi,kfl_ellen_nsi)
           !
           ! Derivatives of the characteristic length: CHALE w.r.t. nodal coordinates
           !
           if(itask == 30) & 
             call elmchl_der(&
                tragl,hleng_der,elcod,elvel,chale_der,pnode,&
                porde,hnatu(pelty),kfl_advec_nsi,kfl_ellen_nsi)
           
           !
           ! Local time step: DTINV_LOC
           !
           ! Maximum time step between that given by the global safety factor saflo_nsi,
           ! and local safet_nsi
           !
           if( kfl_timco == 2 ) then
              call nsi_elmtss(&
                   pelty,pmate,pnode,lnods(1,ielem),ielem,elcod,elvel,&
                   gpcar,chale,hleng,dtcri)
              facto = safet_nsi/safma_nsi 
!              facto = 1.0_rp
              dtinv_loc = min(1.0_rp / (dtcri*safet_nsi), 1.0_rp/(dtcri_nsi*saflo_nsi*sqrt(facto)))
              dtsgs_loc = min(1.0_rp / (dtcri*safet_nsi), 1.0_rp/(dtcri_nsi*saflo_nsi*sqrt(facto)))
              if( kfl_stead_nsi == 1 ) dtinv_loc = 0.0_rp
              if( kfl_timei_nsi == 0 ) dtinv_loc = 0.0_rp
              dtmax_nsi = max(dtmax_nsi,1.0_rp/(dtinv_loc+zeror))  ! Stores maximum time step              
           end if
           !
           ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
           !
           call elmca2(&
                pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
                elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpsha,&
                gpder,gpcar,gphes,ielem)
           !
           ! Derivatives of PGVOL w.r.t. nodal coordinates
           !
           if(itask == 30) & 
             call elmca2_der(&
                pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%deriv,elcod,gpvol_der,gpcar_der,ielem)
           
           !
           ! Enrichement by bubble
           !
           if( kfl_bubbl_nsi /= 0 ) then
              if(    kfl_bubbl_nsi /= FREE_SURFACE_BUBBLE .or. &
                   ( kfl_bubbl_nsi == FREE_SURFACE_BUBBLE .and. lelch(ielem) == ELCUT ) .or. &
                   ( kfl_bubbl_nsi == FREE_SURFACE_BUBBLE .and. minval(elfle(1:pnode)) * maxval(elfle(1:pnode)) < 0.0_rp ) ) then
                 pbubl = 1
                 call elmgeo_bubble(&
                      kfl_bubbl_nsi,ndime,mnode,pnode,pgaus,elcod,gpsha,elmar(pelty)%deriv,gpcar,&
                      gpsha_bub,gpcar_bub,elfle)
              end if
           end if
           
           !
           ! Properties
           !
           call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,gpsha,gpcar)
           call ker_proper('GRDEN','PGAUS',dummi,ielem,gpgde,pnode,pgaus,gpsha,gpcar)
           call ker_proper('VISCO','PGAUS',dummi,ielem,gpvis,pnode,pgaus,gpsha,gpcar)
           ! BIFLU ignores laminar viscosity gradient - BIFL2 includes it
           call ker_proper('POROS','PGAUS',dummi,ielem,gppor,pnode,pgaus,gpsha,gpcar)
           call ker_proper('TURBU','PGAUS',dummi,ielem,gpnut,pnode,pgaus,gpsha,gpcar)

           if( kfl_grvir_nsi == 1 ) then
               call ker_proper('GRTUR','PGAUS',dummi,ielem,grvis,pnode,pgaus,gpsha,gpcar)
               call ker_proper('GRVIS','PGAUS',dummi,ielem,gpgvi,pnode,pgaus,gpsha,gpcar)
           endif

           if( kfl_adj_prob == 1 ) then
!              call ker_proper('DRDEN','PGAUS',dummi,ielem,gpdde,pnode,pgaus,gpsha,gpcar)
!              call ker_proper('GDDEN','PGAUS',dummi,ielem,gpgdd,pnode,pgaus,gpsha,gpcar)
!              call ker_proper('DRVIS','PGAUS',dummi,ielem,gpdvi,pnode,pgaus,gpsha,gpcar)
!              call ker_proper('GDVIS','PGAUS',dummi,ielem,gpgdv,pnode,pgaus,gpsha,gpcar)           
             call ker_proper('TDRTU','PGAUS',dummi,ielem,dgpmut_dtur,pnode,pgaus,gpsha,gpcar)
!              call ker_proper('VDRTU','PGAUS',dummi,ielem,dgpmut_dvel,pnode,pgaus,gpsha,gpcar)
           end if           
           ! Adds gpmut to gpvis 
           call nsi_turbul(&
                itask,1_ip,pnode,pgaus,1_ip,pgaus,gpsha,gpcar,elvel,&
                gpden,gpvis,gpnut,gpgvi,grvis,gpgve,ielem,kfl_kxmod_ker)
          
           if( itask >= 10 .and. itask < 20 ) then
              !
              ! Assemble properties
              !
              call asspro(&
                   itask,pnode,2_ip,pgaus,lnods(1,ielem),lelch(ielem),gpden,gpvis,&
                   gpvol,gpsha,elrhs,prope_nsi)

           else if( itask == 6 ) then

              !-------------------------------------------------------------
              !
              ! Assemble pressure equation only
              !
              !-------------------------------------------------------------

              call nsi_elmpri(&
                   pnode,pgaus,lnods(1,ielem),gpcar,gpvol,gpden,&
                   gpvis,gppor,gpsha,chale,dtinv_loc,elmap)
              !               call assrhs(&
              !                   solve(2) % ndofn,pnode,lnods(1,ielem),elrhs,rhsid)
              !               call assmat(&
              !                   solve(2) % ndofn,pnode,pnode,npoin,solve(2) % kfl_algso,&
              !                   ielem,lnods(1,ielem),elmap,lapla_nsi)

              call nsi_assmat(&
                   -1_ip,pnode,pnode,lnods(1,ielem),elmap,dummr,dummr,&
                   dummr,lapla_nsi)

           else
              !
              ! Residual and RHS
              !
              call nsi_elmres(                                             &
                   pnode,pgaus,plapl,gpsha,gpcar,gphes,gpgvi,gpden,gpvis,  &
                   gppor,gptem,gpsgs,elvel,elpre,elvep,elprp,elgrp,        &
                   eltem,elmsh,elcod,elnor,elcur,elbub,elwmean,hleng,chale,&
                   gpvel,gpgpr,rmomu,rmom2,rcont,gprhs,gprhc,gplap,gpadv,  &
                   gpvep,gpprp,gpgrp,gphyd,gpmsh,gpgve,gpnor,gpcur,gpfle,  &
                   ielem,gprh2,gppre,gprhs_sgs,dtinv_loc,gpgde,gpsha_bub,  &
                   densi)
              !
              ! Exact solution: GPRHS
              !
              call nsi_elmexa(                                            &
                   pgaus,pnode,gpsha,elcod,gpden,gpvis,gppor,gpgvi,cutim, &
                   dummr,gprhs,gprhc,gprh2)
              !
              ! External force: GPRHS
              !
              if( kfl_force_nsi == 1 ) then
                 call nsi_elmexf(                                  &
                      ndime,pgaus,pnode,lforc_material_nsi(pmate),gpden, &
                      xforc_material_nsi(1,pmate),gprhs,gpvel,elcod,gpsha, pmate)
              end if
              !
              ! Compute and assemble subgrid scale: VESGS
              !
              if( kfl_sgscp_nsi == 1 ) then
                 call nsi_updsgs(                                            &
                      pgaus,pnode,ndime,ielem,chale,elvel,gpadv,gpvis,gpden, &
                      rmomu,rmom2,gprhs,gpgpr,gpvel,gpcar,gpsp1,gpsgs,gpsgi, &
                      gpgve,gpvep,gpgrp,gpst1,gprhs_sgs,dtsgs_loc,resis_nsi, &
                      itsta_nsi,rmsgs_nsi,resgs_nsi,gppor)
                 call nsi_subgrid_scale_residual_and_update(                 &
                      ndime,pgaus,ielem,gpsgs,vesgs,resgs_nsi)
              end if

              if( itask == 4 ) then

                 !-------------------------------------------------------------
                 !
                 ! Compute and assemble subgrid scale: VESGS, RESGS_NSI
                 !
                 !-------------------------------------------------------------

                 call nsi_updsgs(                                            &
                      pgaus,pnode,ndime,ielem,chale,elvel,gpadv,gpvis,gpden, &
                      rmomu,rmom2,gprhs,gpgpr,gpvel,gpcar,gpsp1,gpsgs,gpsgi, &
                      gpgve,gpvep,gpgrp,gpst1,gprhs_sgs,dtsgs_loc,resis_nsi, &
                      itsta_nsi,rmsgs_nsi,resgs_nsi,gppor)
                 call nsi_subgrid_scale_residual_and_update(                 &
                      ndime,pgaus,ielem,gpsgs,vesgs,resgs_nsi)

                 if( kfl_stabi_nsi == 1 .or. kfl_stabi_nsi == 2 ) then
                    call nsi_elmsgs(                                            &
                         pgaus,pnode,chale,hleng,gpadv,gpvis,gpden,gpcar,gpst1, &
                         gpst2,gpsp1,gpsp2,gptt1,gptt2,rmomu,gppor,dtsgs_loc,   &
                         dtinv_loc,tamin_nsi,tamax_nsi, elvel,gprhs,gpgpr,rmom2,&
                         gpgve)
                    call nsi_elmort(                                            &
                         ielem,pgaus,pnode,ndime,elvel,elpre,rmomu,rmom2,gprhs, &
                         gpgpr,gpsha,gpvol,gpden,gpadv,gpcar,gpsp1,gpsp2,gpst1  )
                 end if

              else if( itask == 5 ) then

                 !-------------------------------------------------------------
                 !
                 ! Limiter
                 !
                 !-------------------------------------------------------------

                 call nsi_elmsgs(                                               &
                      pgaus,pnode,chale,hleng,gpadv,gpvis,gpden,gpcar,gpst1,    &
                      gpst2,gpsp1,gpsp2,gptt1,gptt2,rmomu,gppor,dtsgs_loc,      &
                      dtinv_loc,tamin_nsi,tamax_nsi, elvel,gprhs,gpgpr,rmom2,   &
                      gpgve)
                 call nsi_asslim(                                               &
                      pnode,pgaus,lnods(1,ielem),gpden,gpsp1,gpadv,gpvep,gpvol, &
                      elvel,gpsha,gpcar,wgrvi,elrhs,rhsid)

              else if( itask /= 4 ) then

                 !-------------------------------------------------------------
                 !
                 ! Assemble equations
                 !
                 !-------------------------------------------------------------
                 !
                 ! Stabilization parameters
                 !
                 call nsi_elmsgs(                                              &
                      pgaus,pnode,chale,hleng,gpadv,gpvis,gpden,gpcar,gpst1,   &
                      gpst2,gpsp1,gpsp2,gptt1,gptt2,rmomu,gppor,dtsgs_loc,     &
                      dtinv_loc,tamin_nsi,tamax_nsi, elvel,gprhs,gpgpr,rmom2,  &
                      gpgve)
                 !
                 ! Weighted mass matrix
                 !
                 if( NSI_FRACTIONAL_STEP ) then
                    call nsi_rhodt_rhotau_nu(                                      &
                         pnode,pgaus,porde,gpsha,gpvol,gpden,gpvis,gppor,gpadv, &
                         chale,dtinv_loc,densi,elrhodt,elrhotau,elmurho)
                 end if
                 !
                 ! Assembly
                 !
                 if(kfl_stabi_nsi == 2 .or. kfl_stabi_nsi == -1 ) then
                    call nsi_elmma4(                                            &
                         pnode,pgaus,pevat,gpden,gpvis,gppor,gpsp1,gpsp2,gpvol, &
                         gpsha,gpcar,gpadv,gpvep,gpprp,gpgrp,gprhs,gprhc,gpvel, &
                         gpsgs,wgrgr,agrau,elvel,elpre,elbub,elauu,elaup,elapp, &
                         elapu,elrbu,elrbp,dtinv_loc,dtsgs_loc,pbubl,gpsha_bub, &
                         gpcar_bub,elauq,elapq,elaqu,elaqp,elaqq,elrbq,densi)
                 else
                    call nsi_elmmat(                                      &
                         pnode,pgaus,pevat,gpden,gpvis,gppor,gpgvi,gpsp1, &
                         gptt1,gpsp2,gptt2,gpvol,gpsha,gpcar,gplap,gphes, &
                         gpadv,gpvep,gprhs,gprhc,rmomu,rcont,p1vec,p2vec, &
                         p2sca,wgrgr,wgrvi,elauu,elaup,elapp,elapu,elrbu, &
                         elrbp,rmom2,p1ve2,gpst1,gpgve,gprh2,gprhs_sgs,   &
                         elvel,ellum,dtinv_loc,pbubl,gpsha_bub,gpcar_bub, &
                         gppre,elauq,elapq,elaqu,elaqp,elaqq,elrbq,elpre, &
                         elbub,densi)

                    !
                    !  calculate cost function F
                    !
                    if (kfl_cos_opt == 1) call nsi_elmcost_all(elvel,pnode,pgaus,gpsha,gpvol,lnods(1,ielem),costf)

                    if(kfl_adj_prob == 1 ) then
                      !
                      ! Calculate strong residual & 
                      !
                      call nsi_elmstr(                                              &
                           pgaus,pnode,ndime,ielem,lnods(1,ielem),chale,elvel,gpadv,gpvis,gpden, &
                           rmomu,rmom2,gprhs,gpvel,gpcar,gpsp1,gpstrm, gpstrc, &
                           gpvep,gpgrp,gpst1,gprhs_sgs,gpsha,gpgve)
                      !
                      ! Modify elmat and elrhs for exact linearization
                      ! Send/receive adjoint rhs vectors to/from other modules
                      !
                      if (itask /= 30) &
                         call nsi_elmmat_der_all(&
                           pnode,pgaus,pevat,gpden,gpvis,gppor,gpgvi,gpsp1, &
                           gptt1,gpsp2,gptt2,gpvol,gpsha,gpcar,gplap,gphes, &
                           gpadv,gpvep,gpprp,gpgrp,elauu,elaup,elapp,elapu,gpst1,&
                           gpstrm,gpstrc,elrbu,elrbp,elvel,ielem,gptem,gpgde,&
                           gpdvi,elaut,elapt,elpre,gpgve,p1vec,gpgdv,gpdde,gpgdd,chale,&
                           dgpmut_dvel,dgpmut_dtur)
                      !
                      !  calculate cost function derivatives w.r.t. unknowns dF/dU
                      !
                      call nsi_elmdcost_all(elvel,pnode,pgaus,gpsha,gpvol,lnods(1:pnode,ielem),elrbu,elrbp)
                      !
                      !  calculate residal derivatives (partial) w.r.t. unknowns dR/dD
                      !
                      if (itask == 30) &
                         call nsi_elmresdiff_all(pnode,pgaus,pevat,gpvol,gpvol_der,gpsha,lnods(1,ielem),&
                           p1vec,p2vec,p2sca,elvel,elpre,wgrgr,chale(2),chale_der(:,:,2),gpsp1,gpsp2,rmomu,gpden,gpadv,gpvis,&
                           gpcar,gpcar_der,rcont,sens_mesh,ielem)
                           
                    endif

                 end if
                 !
                 ! Consistent mass matrix
                 !
                 if( kfl_corre_nsi == 3)  call nsi_elmma5(&
                      pnode,pgaus,pevat,gpden, &
                      gpvol,gpsha,dtinv_loc,elcmm)
                 !
                 ! Shock capturing
                 !
                 if( kfl_shock_nsi /= 0 ) then
                    call nsi_elmshc(                                      &
                         pnode,pgaus,ptopo,pevat,ndime,gpden,gpvel,gprhs, &
                         gpsp1,gpsgs,gpvol,elvel,gpcar,chale,rmomu,rmom2, &
                         elauu)
                 end if
                 !
                 ! Extension elements
                 !
                 if( lelch(ielem) == ELEXT ) then
                    call nsi_elmext(&
                         1_ip,pnode,elauu,elaup,elapu,elapp,elmap,elrbu,elrbp)
                 end if
                 !
                 ! Prescribe Dirichlet boundary conditions
                 !
                 if( kfl_matdi_nsi == 0 ) &
                      call nsi_elmdi3(&
                      pnode,pevat,lnods(1,ielem),&
                      elauu,elaup,elapp,elapu,elrbu,elrbp,elcmm)
                 !
                 ! Transpose of the jacobian matrix for adjoint
                 !
                 if(kfl_adj_prob == 1) then
                   elauu      = transpose(elauu)
                   elauptrans = transpose(elaup)
                   elaputrans = transpose(elapu)
                   elapu      = elauptrans
                   elaup      = elaputrans
                   elapp      = transpose(elapp)
                 endif
                 !
                 ! Schur complement preconditioner: LAPLA_NSI
                 !
                 if( NSI_SCHUR_COMPLEMENT .or. NSI_FRACTIONAL_STEP ) then

                    call elmchl(&
                         tragl,hleng,elcod,elvel,chave,chale,pelty,pnode,&
                         porde,hnatu(pelty),kfl_advec_nsi,kfl_ellsh_nsi)
                    call nsi_elmsch(&
                         pnode,pgaus,lnods(1,ielem),gpcar,gpvol,gpden,&
                         gpvis,gppor,gpsha,elvel,chale,gpsp1,&
                         elmap,dtinv_loc)
                    !
                    ! Extension elements
                    !
                    if( lelch(ielem) == ELEXT ) then
                       call nsi_elmext(&
                            2_ip,pnode,elauu,elaup,elapu,elapp,elmap,elrbu,elrbp)
                    end if
                 end if

                 !
                 ! Assembly: AMATR and RHSID
                 !
                 if( NSI_MONOLITHIC ) then

                    call nsi_assemble_monolithic(&
                         pnode,pevat,lnods(1,ielem),elauu,elaup,elapp,elapu,&
                         elrbu,elrbp,amatr,rhsid)

                 else if( NSI_SCHUR_COMPLEMENT ) then

                    call nsi_assemble_schur(&
                         1_ip,pnode,pevat,ielem,lnods(1,ielem),elauu,elaup,elapp,elapu,&
                         elrbu,elrbp,amatr(poauu_nsi),amatr(poaup_nsi),amatr(poapp_nsi),&
                         amatr(poapu_nsi),rhsid,rhsid(ndbgs_nsi+1))
                    call nsi_assemble_schur(&
                         2_ip,pnode,pevat,ielem,lnods(1,ielem),dummr,dummr,elmap,dummr,&
                         dummr,dummr,dummr,dummr,lapla_nsi,&
                         dummr,dummr,dummr)

                 else if( NSI_FRACTIONAL_STEP ) then

                    call nsi_assembly_fractional_step_scalar(&
                         pnode,pevat,ielem,lnods(:,ielem),elvel,elpre,&
                         elauu,elaup,elapp,elapu,elmap,elrbu,elrbp,&
                         amatr(poaup_nsi:),amatr(poapu_nsi:),&
                         rhsid,rhsid(ndbgs_nsi+1:),lapla_nsi)

                 end if
                 
                 !
                 ! Lumped matrices
                 !
                 ! DT_RHO_NSI  = rho / dt  * M
                 ! MU_RHO_NSI  = mu / rho  * M
                 ! TAU_NSI = rho / tau * M
                 !
                 if( NSI_FRACTIONAL_STEP ) then
                    do inode = 1,pnode
                       ipoin = lnods(inode,ielem)
#ifdef NO_COLORING
                       !$OMP ATOMIC
#endif
                       dt_rho_nsi(ipoin) = dt_rho_nsi(ipoin) + elrhodt(inode)
#ifdef NO_COLORING
                       !$OMP ATOMIC
#endif
                       mass_rho_nsi(ipoin,1) = mass_rho_nsi(ipoin,1) + elmurho(inode)
                       tau_nsi(ipoin)        = tau_nsi(ipoin)        + elrhotau(inode)
                    end do
                 end if
                 !
                 ! Dual time step preconditioner
                 !
                 if( kfl_duatss == 1 ) then
                    do inode = 1,pnode
                       ipoin = lnods(inode,ielem)
#ifdef NO_COLORING
                       !$OMP ATOMIC
#endif
                       lumma(ipoin) = lumma(ipoin) + ellum(inode)
                    end do
                 end if

                 if( kfl_corre_nsi == 3 ) then   ! consistent mass matrix for the mass correction

                    call nsi_assemble_schur(&
                         3_ip,pnode,pevat,ielem,lnods(1,ielem),elcmm,dummr,dummr,dummr,&
                         dummr,dummr,cmama_nsi,dummr,dummr,&
                         dummr,dummr,dummr)

                 end if
                 !
                 ! Pressure bubble assembly
                 !
                 if( kfl_bubbl_nsi /= 0 ) then
                    call nsi_bubble_assembly(&
                         pnode,ielem,elaqq,elaqu,elaqp,elrbq)
                 end if

              end if
           end if

        end if

     end do elements
#ifndef OMPSS
        !$OMP END PARALLEL DO
#endif

#ifndef NO_COLORING
  end do colors
#endif

  call cputim(timeb)  ! CONYO
  rgene = timeb-timea

  !print*,'OLD TOTAL TIME=',timeb-timea

end subroutine nsi_elmope_omp

