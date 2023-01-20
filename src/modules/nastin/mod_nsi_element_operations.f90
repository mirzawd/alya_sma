!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup NastinMatrixAssembly
!> @{
!> @file    mod_nsi_element_operations.f90
!> @author  Guillaume Houzeaux
!> @brief   Navier-Stokes system element assembly and other element
!>          calculations
!> @details Elemental operations, according to ITASK using vectorized subroutines:
!>
!>          \verbatim
!>
!>          1 ........ Element calculations and assembly of global system:
!>                     A <= A^(e): matrix ........................... AMATR
!>                     b <= b^(e): RHS .............................. RHSID
!>                     Q <= Q^(e): Schur complement precond ......... LAPLA_NSI
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
!>
!>          \endverbatim
!>
!>          CORRESPONDANCE OLD TO NEW SUBROUTINES:
!>          --------------------------------------
!>
!>          nsi_elmma4         <= nsi_element_assembly_split_oss
!>          nsi_elmga3         <= nsi_element_operations_gather
!>          nsi_elmlen         <= elmgeo_element_characteristic_length
!>          elmchl             <= elmgeo_element_length
!>          elmca2             <= element_shape_function_derivatives_jacobian
!>          nsi_elmtss         <= nsi_element_time_step
!>          nsi_elmres         <= nsi_element_residual
!>          nsi_updsgs         <= nsi_element_subgrid_scale
!>          nsi_elmsgs         <= nsi_element_stabilization
!>          nsi_elmort         <= nsi_assembly_projections
!>          nsi_elmope_omp     <= nsi_element_operations
!>          nsi_elmmat         <= nsi_element_assembly_asgs_oss
!>                                nsi_element_assembly_asgs_oss_old
!>          nsi_elmdi3         <= nsi_element_dirichlet
!>                                nsi_element_schur
!>          nsi_assemble_schur <= nsi_assembly_schur_method
!>          nsi_elmext         <= nsi_element_extension
!>          nsi_elmexa         <= nsi_element_manufactured_solution
!>          nsi_elmexf         <= nsi_element_external_force
!>          nsi_elmshc         <= nsi_element_shock_capturing
!>
!>          TO BE VECTORIZED:
!>          -----------------
!>
!>          elmgeo_bubble
!>          All optimization subroutines (Mohammad)
!>             nsi_elmcost_all
!>             nsi_element_assembly_der_all
!>             nsi_elmdcost_all
!>             nsi_elmresdiff_all
!>          nsi_elmpri
!>          nsi_asslim
!>          nsi_assemble_monolithic
!>
!------------------------------------------------------------------------

module mod_nsi_element_operations

#include "def_vector_size.inc"
  use def_kintyp,                     only : ip,rp
  use def_master,                     only : solve,kfl_paral,ittim
  use def_master,                     only : INOTMASTER
  use def_domain,                     only : ndime
  use def_elmtyp,                     only : ELFEM
  use def_elmtyp,                     only : ELEXT
  use def_solver,                     only : SOL_MATRIX_HAS_CHANGED
  use def_kermod,                     only : kfl_element_to_csr

  use mod_ker_space_time_function,    only : ker_space_time_function

  use mod_nsi_subgrid_scale,          only : nsi_subgrid_scale_gather
  use mod_nsi_subgrid_scale,          only : nsi_subgrid_scale_residual_and_update
  use mod_element_integration,        only : element_shape_function_derivatives_jacobian
  use mod_nsi_element_assembly,       only : nsi_element_assembly_asgs_oss
  use mod_nsi_element_assembly,       only : nsi_element_assembly_split_oss
  use mod_nsi_assembly_global_system, only : nsi_assembly_monolithic
  use mod_nsi_assembly_global_system, only : nsi_assembly_schur_method
  use mod_nsi_assembly_global_system, only : nsi_assembly_fractional_step
  use mod_nsi_assembly_global_system, only : nsi_assembly_fractional_step_rhs
  use mod_nsi_assembly_global_system, only : nsi_assembly_algebraic_split_oss
  use mod_nsi_assembly_global_system, only : nsi_assembly_dt_rho_tau_nu

  use def_nastin,                     only : kfl_anipo_nsi 
  use def_nastin,                     only : kfl_exacs_nsi
  use def_nastin,                     only : kfl_grad_div_nsi
  use def_nastin,                     only : NSI_GALERKIN
  use def_nastin,                     only : NSI_ASGS
  use def_nastin,                     only : NSI_OSS
  use def_nastin,                     only : NSI_SPLIT_OSS
  use def_nastin,                     only : NSI_ALGEBRAIC_SPLIT_OSS
  use def_nastin,                     only : NSI_SCHUR_COMPLEMENT
  use def_nastin,                     only : NSI_FRACTIONAL_STEP
  use def_nastin,                     only : NSI_DIRICHLET_ELEMENT

  use def_kermod,                     only : kfl_noslw_ker
  use def_kermod,                     only : kfl_nswel_ker
  use def_kermod,                     only : normal_nsw_ker
  use def_domain,                     only : lpoty
  
  implicit none
  private

  real(rp), parameter :: zeror = epsilon(1.0_rp)

  interface nsi_element_operations_gather
     module procedure nsi_element_operations_gather_scalar,&
          &           nsi_element_operations_gather_vector
  end interface nsi_element_operations_gather

  interface nsi_rhodt_rhotau_nu
     module procedure nsi_rhodt_rhotau_nu_scalar,&
          &           nsi_rhodt_rhotau_nu_vector
  end interface nsi_rhodt_rhotau_nu

  public :: nsi_element_operations_gather   ! solucion trucha cambiar antes de commitear 
  public :: nsi_element_operations
  public :: nsi_rhodt_rhotau_nu

contains




#if defined PNODE_VALUE && defined PGAUS_VALUE 
  subroutine nsi_element_operations(&
       itask,qnode,qgaus,list_elements,resis_nsi,itsta_nsi,resgs_nsi,&
       tamin_nsi,rmsgs_nsi,tamax_nsi,dtmax_nsi,time1)
#else
  subroutine nsi_element_operations(&
       itask,pnode,pgaus,list_elements,resis_nsi,itsta_nsi,resgs_nsi,&
       tamin_nsi,rmsgs_nsi,tamax_nsi,dtmax_nsi,time1)
#endif
    
    use def_elmtyp,            only : ELCUT
    use def_kintyp,            only : ip,rp
    use def_master,            only : amatr,rhsid,cutim,kfl_timco       ! CONYO
    use def_master,            only : lumma,vesgs
    use def_kermod,            only : kfl_cos_opt,kfl_adj_prob,costf
    use def_kermod,            only : kfl_duatss
    use def_domain,            only : ltype
    use def_domain,            only : llapl,lorde,ltopo
    use def_domain,            only : lnods
    use def_domain,            only : lelch,elmar,mnode
    use def_domain,            only : lmate,ntens
    use mod_elmgeo,            only : element_type
    use mod_elmgeo,            only : elmgeo_bubble
    use mod_elmgeo,            only : FREE_SURFACE_BUBBLE
    use mod_elmgeo,            only : QUADRATIC_BUBBLE
    use mod_elmgeo,            only : FLAT_BUBBLE
    use mod_elmgeo,            only : elmgeo_element_length
    use mod_elmgeo,            only : elmgeo_element_characteristic_length
    use def_nastin,            only : ncomp_nsi,nbdfp_nsi
    use def_nastin,            only : kfl_stabi_nsi
    use def_nastin,            only : kfl_advec_nsi
    use def_nastin,            only : dtinv_nsi,safet_nsi, safma_nsi
    use def_nastin,            only : dtcri_nsi,saflo_nsi
    use def_nastin,            only : dtsgs_nsi,kfl_stead_nsi
    use def_nastin,            only : kfl_timei_nsi,NSI_MONOLITHIC
    use def_nastin,            only : kfl_grvir_nsi
    use def_nastin,            only : lapla_nsi,kfl_sgscp_nsi
    use def_nastin,            only : kfl_force_nsi,kfl_shock_nsi, kfl_vegeo_time_nsi
    use def_nastin,            only : kfl_matdi_nsi,poauu_nsi
    use def_nastin,            only : poaup_nsi,poapp_nsi
    use def_nastin,            only : poapu_nsi,ndbgs_nsi
    use def_nastin,            only : kfl_predi_nsi,kfl_ellen_nsi
    use def_nastin,            only : kfl_ellsh_nsi
!    use def_nastin,            only : resdiff_nsi
    use def_nastin,            only : kfl_bubbl_nsi
    use def_nastin,            only : dt_rho_nsi
    use def_nastin,            only : mass_rho_nsi
    use def_nastin,            only : tau_nsi
    use def_nastin,            only : kfl_immer_nsi
    use mod_nsi_bubble,        only : nsi_bubble_assembly
    use mod_nsi_bubble,        only : nsi_eliminate_bubble
    use mod_ker_proper,        only : ker_proper
    use mod_ker_nsw_visc,      only : ker_nsw_visc

    !
    ! AB: densi
    !
    use mod_ker_tendencies,  only : kfl_tendencies_ker
    implicit none
    !
    ! Input and output variables
    !
    integer(ip), intent(in)          :: itask                        !< What to do
#if defined PNODE_VALUE && defined PGAUS_VALUE 
    integer(ip), intent(in)          :: qnode                        !< Number of nodes
    integer(ip), intent(in)          :: qgaus                        !< Number of Gauss points
    integer(ip), parameter           :: pnode = PNODE_VALUE
    integer(ip), parameter           :: pgaus = PGAUS_VALUE 
#else
    integer(ip), intent(in)          :: pnode                        !< Number of nodes
    integer(ip), intent(in)          :: pgaus                        !< Number of Gauss points
#endif
    integer(ip), intent(in)          :: list_elements(VECTOR_SIZE)   !< List of elements

    real(rp),    intent(inout)       :: resis_nsi(*)
    integer(ip), intent(inout)       :: itsta_nsi(*)
    real(rp),    intent(inout)       :: resgs_nsi(2)
    real(rp),    intent(inout)       :: tamin_nsi                    ! Minimum tau
    real(rp),    intent(inout)       :: rmsgs_nsi
    real(rp),    intent(inout)       :: tamax_nsi                    ! Maximum tau
    real(rp),    intent(inout)       :: dtmax_nsi                    ! Maximum local time step
    real(rp),    intent(inout)       :: time1(10)                    ! Timings
    !
    ! Element matrices and vectors (stiffness and preconditioner)
    !
    real(rp)    :: elauu(VECTOR_SIZE,pnode*ndime,pnode*ndime)        ! Auu
    real(rp)    :: elaup(VECTOR_SIZE,pnode*ndime,pnode)              ! Aup
    real(rp)    :: elapp(VECTOR_SIZE,pnode,pnode)                    ! App
    real(rp)    :: elapu(VECTOR_SIZE,pnode,pnode*ndime)              ! Apu
    real(rp)    :: elrbu(VECTOR_SIZE,ndime,pnode)                    ! bu
    real(rp)    :: elporrbu(VECTOR_SIZE,ndime,pnode)                 ! the part of bu due to anisotropic porosity - it has alreday been added to elrbu but I need it to calculate porous force 
    real(rp)    :: elrbp(VECTOR_SIZE,pnode)                          ! bp
!    real(rp)    :: elrhs(VECTOR_SIZE,6*pnode)                        ! Generic RHS
    real(rp)    :: elmap(VECTOR_SIZE,pnode,pnode)                    ! Q
    real(rp)    :: ellum(VECTOR_SIZE,pnode)                          ! Lumped mass matrix
    real(rp)    :: elcmm(VECTOR_SIZE,pnode*ndime,pnode*ndime)        ! Consistent mass matrix
    real(rp)    :: eldtrho(VECTOR_SIZE,pnode)                        ! Projection of rho/dt
    real(rp)    :: elmurho(VECTOR_SIZE,pnode)                        ! Projection of mu/rho
    real(rp)    :: eltau(VECTOR_SIZE,pnode)                          ! Projection of rho/tau
    !
    ! Bubble matrices
    !
    real(rp)    :: elauq(VECTOR_SIZE,pnode*ndime,1)
    real(rp)    :: elapq(VECTOR_SIZE,pnode,1)
    real(rp)    :: elaqu(VECTOR_SIZE,1,pnode*ndime)
    real(rp)    :: elaqp(VECTOR_SIZE,1,pnode)
    real(rp)    :: elaqq(VECTOR_SIZE,1,1)
    real(rp)    :: elrbq(VECTOR_SIZE,1)
    !
    ! Gather
    !
    real(rp)    :: elvel(VECTOR_SIZE,ndime,pnode,ncomp_nsi+1)        ! u
    real(rp)    :: elpre(VECTOR_SIZE,pnode,ncomp_nsi-1)              ! p
    real(rp)    :: elfle(VECTOR_SIZE,pnode)                          ! phi
    real(rp)    :: elcod(VECTOR_SIZE,ndime,pnode)                    ! x
    real(rp)    :: elvep(VECTOR_SIZE,ndime,pnode)                    ! Pi(momentum)
    real(rp)    :: elprp(VECTOR_SIZE,pnode)                          ! Pi(div(u))
    real(rp)    :: elgrp(VECTOR_SIZE,ndime,pnode)                    ! Pi(grad(p))
    real(rp)    :: eltem(VECTOR_SIZE,pnode,nbdfp_nsi)                ! T
    real(rp)    :: elwmean(VECTOR_SIZE,pnode,nbdfp_nsi)              ! W mean
    real(rp)    :: elmsh(VECTOR_SIZE,ndime,pnode)                    ! u mesh
    real(rp)    :: elnor(VECTOR_SIZE,ndime,pnode)                    ! Normal to the Level Set interface
    real(rp)    :: elcur(VECTOR_SIZE,pnode)                          ! Level Set interface curvature
    real(rp)    :: elbub(VECTOR_SIZE)                                ! Element bubble
    real(rp)    :: elden(VECTOR_SIZE,pnode)                                ! Element density (for low Mach and compressible)
    real(rp)    :: ellag(VECTOR_SIZE,ndime,ndime,pnode)              ! Lagrange multiplier for IBM
    !
    ! Gather No slip wall law
    !
    real(rp)    :: elibopo(VECTOR_SIZE,pnode)
    real(rp)    :: elnnsw(VECTOR_SIZE,ndime)
    real(rp)    :: elavv(VECTOR_SIZE,ndime,pnode)
    real(rp)    :: elywal(VECTOR_SIZE)
    !
    ! Indices and dimensions
    !
    integer(ip) :: ielem,inode,ivect
    integer(ip) :: pevat,dummi,ptopo
    integer(ip) :: pelty,plapl,porde
    integer(ip) :: ipoin
    integer(ip) :: pbubl(VECTOR_SIZE)
    integer(ip) :: lmate_loc(VECTOR_SIZE)
    integer(ip) :: lnods_loc(VECTOR_SIZE,pnode)
    integer(ip) :: lelch_loc(VECTOR_SIZE)
    integer(ip) :: list_elements_p(VECTOR_SIZE)                      ! List of elements (always positive)
    !
    ! Gauss point values
    !
    real(rp)    :: gpsha(VECTOR_SIZE,pnode,pgaus)                    ! N
    real(rp)    :: gpder(VECTOR_SIZE,ndime,pnode,pgaus)              ! dN/dsi
    real(rp)    :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)              ! dN/dxi
    real(rp)    :: gphes(VECTOR_SIZE,ntens,mnode,pgaus)              ! d2N/dxidxj
    real(rp)    :: gplap(VECTOR_SIZE,pnode,pgaus)                    ! Lapl(N)
    real(rp)    :: gpvol(VECTOR_SIZE,pgaus)                          ! w*|J|, |J|
    real(rp)    :: gpvis(VECTOR_SIZE,pgaus)                          ! Viscosity
    real(rp)    :: gpvis_nsw(VECTOR_SIZE,pgaus)                      ! Viscosity for no slip wall
    real(rp)    :: gpgvi(VECTOR_SIZE,ndime,pgaus)                    ! Viscosity gradients
    real(rp)    :: gpmut(VECTOR_SIZE,pgaus)                          ! mut
    real(rp)    :: gpnut(VECTOR_SIZE,pgaus)                           ! nut
    real(rp)    :: grvis(VECTOR_SIZE,ndime,pgaus)                    ! grad(mut)
    real(rp)    :: gppor(VECTOR_SIZE,pgaus)                          ! Porosity
    real(rp)    :: gpapo(VECTOR_SIZE,ndime,ndime,pgaus)              ! Porosity
    real(rp)    :: gpden(VECTOR_SIZE,pgaus)                          ! Density
    real(rp)    :: gpfle(VECTOR_SIZE,pgaus)                          ! Level set function
    real(rp)    :: gpmix(VECTOR_SIZE,pgaus)                          ! Mixing function
    real(rp)    :: gpst1(VECTOR_SIZE,pgaus)                          ! tau1
    real(rp)    :: gpst2(VECTOR_SIZE,pgaus)                          ! tau2
    real(rp)    :: gpsp1(VECTOR_SIZE,pgaus)                          ! tau1'
    real(rp)    :: gpsp2(VECTOR_SIZE,pgaus)                          ! tau2'
    real(rp)    :: gptt1(VECTOR_SIZE,pgaus)                          ! tau1'/tau1
    real(rp)    :: gptt2(VECTOR_SIZE,pgaus)                          ! tau2'/tau2
    real(rp)    :: gpadv(VECTOR_SIZE,ndime,pgaus)                    ! u+u'
    real(rp)    :: gprhs(VECTOR_SIZE,ndime,pgaus)                    ! RHS
    real(rp)    :: gprhs_sgs(VECTOR_SIZE,ndime,pgaus)                ! RHS due to subscales
    real(rp)    :: gprhc(VECTOR_SIZE,pgaus)                          ! RHS for the continuity equation (Low Mach)
    real(rp)    :: gprh2(VECTOR_SIZE,pgaus)                          ! RHS for the residual of continuity equation (Low Mach)
    real(rp)    :: gpsgs(VECTOR_SIZE,ndime,pgaus,2)                  ! u'
    real(rp)    :: gpsgt(VECTOR_SIZE,pgaus,nbdfp_nsi)                ! T'
    real(rp)    :: gppre(VECTOR_SIZE,pgaus,ncomp_nsi-1)              ! p
    real(rp)    :: gpvel(VECTOR_SIZE,ndime,pgaus,ncomp_nsi-1)        ! u
    real(rp)    :: gpgve(VECTOR_SIZE,ndime,ndime,pgaus)              ! grad(u)
    real(rp)    :: gpgpr(VECTOR_SIZE,ndime,pgaus,2)                  ! grad(p)
    real(rp)    :: gpgde(VECTOR_SIZE,ndime,pgaus)                    ! grad(den)
    real(rp)    :: gpdde(VECTOR_SIZE,pnode,pgaus)                    ! Density   derivatives w.r.t nodal temperature
    real(rp)    :: gpgdd(VECTOR_SIZE,ndime,pnode,pgaus)              ! Density   derivatives w.r.t nodal temperature and coordinates
    real(rp)    :: gpdvi(VECTOR_SIZE,pnode,pgaus)                    ! Viscosity derivatives w.r.t nodal temperature
    real(rp)    :: gpgdv(VECTOR_SIZE,ndime,pnode,pgaus)              ! Viscosity derivatives w.r.t nodal temperature and coordinates
    real(rp)    :: gptem(VECTOR_SIZE,pgaus,ncomp_nsi)                ! T
    real(rp)    :: gpsgi(VECTOR_SIZE,ndime,pgaus)                    ! SGS (work array)
    real(rp)    :: gpvep(VECTOR_SIZE,ndime,pgaus)                    ! -tau1'*R(u) or tau1*rho*(a.grad)u
    real(rp)    :: gpprp(VECTOR_SIZE,pgaus)                          ! tau2*div(u)
    real(rp)    :: gpgrp(VECTOR_SIZE,ndime,pgaus)                    ! tau1'*( grad(p) - rho*f )
    real(rp)    :: gphyd(VECTOR_SIZE,pgaus)                          ! rho_hyd
    real(rp)    :: gpmsh(VECTOR_SIZE,ndime,pgaus)                    ! u_mesh
    real(rp)    :: gpnor(VECTOR_SIZE,ndime,pgaus)                    ! LS normal
    real(rp)    :: gpcur(VECTOR_SIZE,pgaus)                          ! LS curvature
    real(rp)    :: densi(VECTOR_SIZE,pgaus,nbdfp_nsi)                ! Density at previous time-steps
    real(rp)    :: gplag(VECTOR_SIZE,ndime,ndime,pgaus)              ! Lagrange multiplier
    !
    ! Enrichement
    !
    real(rp)    :: gpsha_bub(VECTOR_SIZE,pgaus)                      ! Ne
    real(rp)    :: gpcar_bub(VECTOR_SIZE,ndime,pgaus)                ! dNe/dxi
    !
    ! Element characteristics (to compute tau1 and tau2)
    !
    real(rp)    :: tragl(VECTOR_SIZE,ndime,ndime)
    real(rp)    :: chave(VECTOR_SIZE,ndime,2)
    real(rp)    :: chale(VECTOR_SIZE,2)
    real(rp)    :: hleng(VECTOR_SIZE,ndime)
    real(rp)    :: dtcri(VECTOR_SIZE)
    real(rp)    :: dtinv_loc(VECTOR_SIZE)
    real(rp)    :: dtsgs_loc(VECTOR_SIZE)
    real(rp)    :: baloc(VECTOR_SIZE,ndime)
    real(rp)    :: dummr(2),dumm2(2,2),dumm3(2,2,2)
    !
    ! Perturbation and residuals
    !
    real(rp)    :: rmomu(VECTOR_SIZE,pnode,pgaus)                    ! Residual velocity in momentum
    real(rp)    :: rmom2(VECTOR_SIZE,ndime,ndime,pnode,pgaus)        ! Residual velocity in momentum
    real(rp)    :: rcont(VECTOR_SIZE,ndime,pnode,pgaus)              ! Residual velocity in continuity
    !
    ! Exact linearization and adjoint
    !
!    real(rp)    :: gpstrm(VECTOR_SIZE,ndime,pgaus)                   ! Momentum strong residual
!    real(rp)    :: gpstrc(VECTOR_SIZE,pgaus)                         ! Continuity strong residual
!    real(rp)    :: elaut(VECTOR_SIZE,pnode*ndime,pnode)
!    real(rp)    :: elapt(VECTOR_SIZE,pnode,pnode)
    real(rp)    :: elaputrans(VECTOR_SIZE,pnode*ndime,pnode)         ! Transpose of the Apu
    real(rp)    :: elauptrans(VECTOR_SIZE,pnode,pnode*ndime)         ! Transpose of the Aup
    real(rp)    :: dgpmut_dvel(VECTOR_SIZE,pnode,ndime,pgaus)        ! Turbulence viscosity derivatives w.r.t. nodal velocity
    real(rp)    :: dgpmut_dtur1(VECTOR_SIZE,pnode,pgaus)             ! Turbulence viscosity derivatives w.r.t first turbulence unk
    real(rp)    :: dgpmut_dtur2(VECTOR_SIZE,pnode,pgaus)             ! Turbulence viscosity derivatives w.r.t second turbulence unk

    integer(ip) :: nsize
    real(rp)    :: timea,timeb, facto

    !--------------------------------------------------------------------
    !
    ! Gather: global to local
    !
    !--------------------------------------------------------------------

    ielem = list_elements(1)
    pelty = abs(ltype(ielem))
    plapl = llapl(pelty)
    porde = lorde(pelty)
    ptopo = ltopo(pelty)
    pevat = ndime * pnode
    nsize = VECTOR_SIZE
    elcmm = 0.0_rp

    if(  kfl_stabi_nsi == NSI_GALERKIN .or. &
         kfl_stabi_nsi == NSI_ALGEBRAIC_SPLIT_OSS ) plapl = 0

    call cputim(timea)

    !--------------------------------------------------------------------
    !
    ! List of elements
    !
    ! +----+----+----+----+
    ! | 23 | 24 | 25 | 0  | <= list_elements
    ! +----+----+----+----+
    !
    ! +----+----+----+----+
    ! | 23 | 24 | 25 | 23 | <= list_elements_p
    ! +----+----+----+----+
    !
    !--------------------------------------------------------------------

    list_elements_p = list_elements
    do ivect = 1,VECTOR_SIZE
       ielem = abs(list_elements(ivect))
       if( ielem /= 0 ) then
          lnods_loc(ivect,1:pnode) = lnods(1:pnode,ielem)
          lelch_loc(ivect)         = lelch(ielem)
          lmate_loc(ivect)         = lmate(ielem)
       else
          list_elements_p(ivect)   = list_elements(1)
          lnods_loc(ivect,1:pnode) = 0
          lelch_loc(ivect)         = ELFEM
          lmate_loc(ivect)         = 1
       end if
    end do

    call nsi_element_operations_gather(&
         pnode,pgaus,list_elements,lnods_loc,elcod,elpre,elvel,&
         elfle,elvep,elprp,elgrp,eltem,elmsh,elnor,elcur,&
         elwmean,gphyd,gpsgt,gpsgs,elbub,ellag,elden,elavv,elibopo,&
         elnnsw,elywal)


    call cputim(timeb)
    time1(1) = time1(1) + timeb - timea

    !--------------------------------------------------------------------
    !
    ! Element shape functions and derivatives
    !
    !--------------------------------------------------------------------

    call cputim(timea)
    call element_shape_function_derivatives_jacobian(&
         pnode,pgaus,plapl,elmar(pelty) % weigp,elmar(pelty) % shape,&
         elmar(pelty) % deriv,elmar(pelty) % heslo,&
         elcod,gpvol,gpsha,gpder,gpcar,gphes,&
         list_elements)

    !--------------------------------------------------------------------
    !
    ! ELement characteristic length: HLENG, CHALE, CHAVE
    !
    !--------------------------------------------------------------------

    hleng = 0.0_rp
    chave = 0.0_rp
    chale = 0.0_rp

    call elmgeo_element_characteristic_length(&
         ndime,pnode,elmar(pelty) % dercg(:,:),elcod,hleng,element_type(pelty) % natural_length,tragl)
    call elmgeo_element_length(&
         ndime,pnode,porde,tragl,hleng,elcod,elvel(:,:,:,1),chave,chale,&
         element_type(pelty) % natural_length,kfl_advec_nsi,kfl_ellen_nsi)

    !--------------------------------------------------------------------
    !
    ! Enrichement by bubble
    !
    !--------------------------------------------------------------------

    gpsha_bub  = 0.0_rp
    gpcar_bub  = 0.0_rp
    pbubl      = 0

    if( kfl_bubbl_nsi /= 0 ) then
       do ivect = 1,VECTOR_SIZE
          ielem = list_elements(ivect)
          if( ielem > 0 ) then
             pelty = ltype(ielem)
             if( pelty > 0 ) then
                if(     kfl_bubbl_nsi == FREE_SURFACE_BUBBLE ) then
                   if( lelch_loc(ivect) == ELCUT  .or. &
                        ( minval(elfle(ivect,1:pnode)) * maxval(elfle(ivect,1:pnode)) < 0.0_rp ) ) then
                      pbubl(ivect) = 1
                   end if
                else if( kfl_bubbl_nsi == FLAT_BUBBLE .or. &
                     &   kfl_bubbl_nsi == QUADRATIC_BUBBLE ) then
                   pbubl(ivect) = 1
                end if
                if( pbubl(ivect) == 1 ) then
                   call elmgeo_bubble(&
                        kfl_bubbl_nsi,ndime,mnode,pnode,pgaus,elcod(ivect,:,:),gpsha(ivect,:,:),elmar(pelty) % deriv,gpcar(ivect,:,:,:),&
                        gpsha_bub(ivect,:),gpcar_bub(ivect,:,:),elfle(ivect,:))
                end if
             end if
          end if
       end do
    end if

    call cputim(timeb)
    time1(2) = time1(2) + timeb - timea
    call cputim(timea)

    !--------------------------------------------------------------------
    !
    ! Properties and local time step
    !
    !--------------------------------------------------------------------

    call ker_proper('DENSI','PGAUS',dummi,list_elements_p,gpden,pnode,pgaus,gpsha,gpcar)      ! rho
    call ker_proper('GRDEN','PGAUS',dummi,list_elements_p,gpgde,pnode,pgaus,gpsha,gpcar)      ! grad(rho)
    call ker_proper('VISCO','PGAUS',dummi,list_elements_p,gpvis,pnode,pgaus,gpsha,gpcar)      ! mu
    call ker_proper('POROS','PGAUS',dummi,list_elements_p,gppor,pnode,pgaus,gpsha,gpcar)      ! sig
    call ker_proper('TURBU','PGAUS',dummi,list_elements_p,gpnut,pnode,pgaus,gpsha,gpcar)      ! mut
    if( kfl_anipo_nsi /= 0 ) &
         call ker_proper('ANIPO','PGAUS',dummi,list_elements_p,gpapo,pnode,pgaus,gpsha,gpcar) ! sig

    gpmut = gpden * gpnut
    if ( kfl_noslw_ker /= 0_ip) call ker_nsw_visc(ndime,pnode,pgaus,list_elements,elavv,gpcar,elnnsw,gpvis,&
         gpden,gpmut,elibopo,elywal,elcod,gpvis_nsw)  ! obtains gpvis_nsw   
    
    gpvis = gpvis + gpmut   ! Effective viscosity <= mu+mut
    !
    ! Viscosity gradient
    !
    if( kfl_grvir_nsi == 0 ) then
       gpgvi = 0.0_rp
       grvis = 0.0_rp
    else
       call ker_proper('GRVIS','PGAUS',dummi,list_elements_p,gpgvi,pnode,pgaus,gpsha,gpcar)  ! grad(mu)
       call ker_proper('GRTUR','PGAUS',dummi,list_elements_p,grvis,pnode,pgaus,gpsha,gpcar)  ! grad(mut)
       gpgvi = gpgvi + grvis                                                                 ! Effective viscosity gradient <= grad(mu)+grad(mut)
    end if
    !
    ! Adjoint problem
    !
    if( kfl_adj_prob == 1 ) then
       call ker_proper('DRDEN','PGAUS',dummi,list_elements_p,gpdde,pnode,pgaus,gpsha,gpcar)
       call ker_proper('GDDEN','PGAUS',dummi,list_elements_p,gpgdd,pnode,pgaus,gpsha,gpcar)
       call ker_proper('DRVIS','PGAUS',dummi,list_elements_p,gpdvi,pnode,pgaus,gpsha,gpcar)
       call ker_proper('GDVIS','PGAUS',dummi,list_elements_p,gpgdv,pnode,pgaus,gpsha,gpcar)
    end if
    !
    ! Immersed boundary mixing
    !
    if( kfl_immer_nsi /= 0 ) then
       call ker_proper('MIXIN','PGAUS',dummi,list_elements_p,gpmix,pnode,pgaus,gpsha,gpcar)  ! Mixing function
    end if

    call cputim(timeb)
    time1(3) = time1(3) + timeb - timea
    !
    ! Local time steps: DTINV_LOC and DTSGS_LOC
    !
    if( kfl_timco == 2 ) then ! local time step
       !
       ! Maximum time step between that given by the global safety factor saflo_nsi,
       ! and local safet_nsi
       !
       call nsi_element_time_step(&
            pnode,pgaus,list_elements,elcod,elvel,&
            chale,hleng,dtcri,gpden,gpvis,gpmut,gppor)
       
       facto = safet_nsi/safma_nsi
       !facto =1.0_rp
       dtinv_loc(1:VECTOR_SIZE) = min(1.0_rp / (dtcri(1:VECTOR_SIZE)*safet_nsi), 1.0_rp/(dtcri_nsi*saflo_nsi*sqrt(facto)))       
       
       
       if( kfl_stead_nsi == 1 ) dtinv_loc(1:VECTOR_SIZE) = 0.0_rp
       if( kfl_timei_nsi == 0 ) dtinv_loc(1:VECTOR_SIZE) = 0.0_rp
       dtsgs_loc(1:VECTOR_SIZE) = dtinv_loc(1:VECTOR_SIZE)
       dtmax_nsi = max(dtmax_nsi,1.0_rp/maxval(dtinv_loc(1:VECTOR_SIZE)))  ! Stores maximum time step

    else
       !
       ! Use global time step
       !
       dtinv_loc(1:VECTOR_SIZE) = dtinv_nsi
       dtsgs_loc(1:VECTOR_SIZE) = dtsgs_nsi

    end if

    !--------------------------------------------------------------------
    !
    ! Element residual
    !
    !--------------------------------------------------------------------

    call cputim(timea)

    gptem        = 0.0_rp
    gpgpr        = 0.0_rp
    rmomu        = 0.0_rp
    rmom2        = 0.0_rp
    rcont        = 0.0_rp
    gprhs        = 0.0_rp
    gprhc        = 0.0_rp
    gplap        = 0.0_rp
    gpadv        = 0.0_rp
    gpvep        = 0.0_rp
    gpprp        = 0.0_rp
    gpgrp        = 0.0_rp
    gpmsh        = 0.0_rp
    gpgve        = 0.0_rp
    gpnor        = 0.0_rp
    gpcur        = 0.0_rp
    gprh2        = 0.0_rp
    gppre        = 0.0_rp
    gprhs_sgs    = 0.0_rp

    dgpmut_dvel  = 0.0_rp
    dgpmut_dtur1 = 0.0_rp
    dgpmut_dtur2 = 0.0_rp

    call nsi_element_residual(&
         pnode,pgaus,plapl,gpsha,gpcar,gphes,gpgvi,gpden,gpvis,    &
         gppor,gpapo,gptem,gpsgs,elvel,elpre,elvep,elprp,elgrp,ellag, &
         eltem,elmsh,elcod,elnor,elcur,elbub,elwmean,hleng,chale,  &
         gpvel,gpgpr,rmomu,rmom2,rcont,gprhs,gprhc,gplap,gpadv,    &
         gpvep,gpprp,gpgrp,gphyd,gpmsh,gpgve,gpnor,gpcur,gpfle,    &
         ielem,gprh2,gppre,gprhs_sgs,dtinv_loc,gpgde,gpsgt,        &
         gpsha_bub,cutim,densi,gplag,elden,list_elements)



    call cputim(timeb)
    time1(4) = time1(4) + timeb - timea

    !--------------------------------------------------------------------
    !
    ! Exact solution and external force
    !
    !--------------------------------------------------------------------
    !
    ! External force: GPRHS
    !
    if( kfl_force_nsi == 1 .or.kfl_tendencies_ker.or.kfl_vegeo_time_nsi.gt.0) then
       call nsi_element_external_force(&
            pnode, pgaus,list_elements,lnods_loc, gpsha,lmate_loc,gpden,gprhs,gpvel,elcod, gppor, gpapo)
    end if
    !
    ! Exact solution: GPRHS
    !
    if( kfl_exacs_nsi > 0 ) then
       call nsi_element_manufactured_solution(&
            pgaus,pnode,gpsha,elcod,gpden,gpvis,gppor,gpgvi,&
            cutim,baloc,gprhs,gprhc,gprh2)
    end if

    if( itask == 4 ) then

       !--------------------------------------------------------------------
       !
       ! Coupled grid scale and subgrid scale ... ITASK = 1, KFL_SGSCP_NSI = 1
       ! Update subgrid scale ................... ITASK = 4
       ! Orthogonol SGS or Split OSS ............ ITASK = 4, KFL_STABI_NSI = 1,2
       !
       !--------------------------------------------------------------------

       call cputim(timea)

       call nsi_element_subgrid_scale(                             &
            pgaus,pnode,ndime,ielem,chale,elvel,gpadv,gpvis,gpden, &
            rmomu,rmom2,gprhs,gpgpr,gpvel,gpcar,gpsp1,gpsgs,gpsgi, &
            gpgve,gpvep,gpgrp,gpst1,gprhs_sgs,dtsgs_loc,resis_nsi, &
            itsta_nsi,rmsgs_nsi,resgs_nsi,gppor)
       call nsi_subgrid_scale_residual_and_update(                 &
            nsize,ndime,pgaus,list_elements,gpsgs,vesgs,resgs_nsi)

       if( kfl_stabi_nsi == NSI_OSS .or. kfl_stabi_nsi == NSI_SPLIT_OSS ) then
          call nsi_element_stabilization(                               &
               pgaus,pnode,chale,hleng,gpadv,gpvis,gpden,gpcar,gpst1,   &
               gpst2,gpsp1,gpsp2,gptt1,gptt2,rmomu,gppor,dtsgs_loc,     &
               dtinv_loc,tamin_nsi,tamax_nsi)
          call nsi_assembly_projections(&
               list_elements,pgaus,pnode,ndime,lnods_loc,elvel,elpre,&
               rmomu,rmom2,gprhs,gpgpr,gpsha,gpvol,gpden,gpadv,gpcar,&
               gpsp1,gpsp2,gpst1)
       end if

       call cputim(timeb)
       time1(5) = time1(5) + timeb - timea

    else if( itask == 5 ) then

       !--------------------------------------------------------------------
       !
       ! Limiter
       !
       !--------------------------------------------------------------------

       call runend('ASSEMBLY NOT CODED FOR LIMITER')
       call nsi_element_stabilization(                               &
            pgaus,pnode,chale,hleng,gpadv,gpvis,gpden,gpcar,gpst1,   &
            gpst2,gpsp1,gpsp2,gptt1,gptt2,rmomu,gppor,dtsgs_loc,     &
            dtinv_loc,tamin_nsi,tamax_nsi)
       !call nsi_asslim(                                               &
       !     pnode,pgaus,lnods(1,ielem),gpden,gpsp1,gpadv,gpvep,gpvol, &
       !     elvel,gpsha,gpcar,wgrvi,elrhs,rhsid)

    else if( itask == 6 ) then

       !-----------------------------------------------------------------
       !
       ! Assemble pressure equation only
       !
       !-----------------------------------------------------------------

       call nsi_element_laplacian(&
            pnode,pgaus,gpcar,gpden,gpvol,elmap)
       call nsi_assembly_schur_method(&
            2_ip,pnode,pevat,list_elements,lnods_loc,dumm3,dumm3,elmap,dumm3,&
            dumm3,dumm2,dummr,dummr,lapla_nsi,&
            dummr,dummr,dummr)

    else if( itask == 1 ) then

       !-----------------------------------------------------------------
       !
       ! Element matrices and RHS and assembly in global system
       !
       !-----------------------------------------------------------------

       call cputim(timea)

       if( kfl_sgscp_nsi == 1 ) then
          call nsi_element_subgrid_scale(                             &
               pgaus,pnode,ndime,ielem,chale,elvel,gpadv,gpvis,gpden, &
               rmomu,rmom2,gprhs,gpgpr,gpvel,gpcar,gpsp1,gpsgs,gpsgi, &
               gpgve,gpvep,gpgrp,gpst1,gprhs_sgs,dtsgs_loc,resis_nsi, &
               itsta_nsi,rmsgs_nsi,resgs_nsi,gppor)
          call nsi_subgrid_scale_residual_and_update(                 &
               nsize,ndime,pgaus,list_elements,gpsgs,vesgs,           &
               resgs_nsi)
       end if
       
       call nsi_element_stabilization(                               &
            pgaus,pnode,chale,hleng,gpadv,gpvis,gpden,gpcar,gpst1,   &
            gpst2,gpsp1,gpsp2,gptt1,gptt2,rmomu,gppor,dtsgs_loc,     &
            dtinv_loc,tamin_nsi,tamax_nsi)

       call cputim(timeb)
       time1(6) = time1(6) + timeb - timea
       !
       ! Projections of rho/dt, mu/rho and rho/tau
       !
       if( kfl_predi_nsi == 7 .or. kfl_predi_nsi == 8 .or. kfl_predi_nsi == 9 ) then
          call nsi_rhodt_rhotau_nu(                                      &
               pnode,pgaus,porde,gpsha,gpvol,gpden,gpvis,gppor,gpadv, &
               chale,dtinv_loc,densi,eldtrho,eltau,elmurho)
       end if
       !
       ! Element matrices and RHS
       !
       call cputim(timea)

       if( itask == 1 ) then
          if(  kfl_stabi_nsi == NSI_SPLIT_OSS           .or. &
               kfl_stabi_nsi == NSI_GALERKIN            .or. &
               kfl_stabi_nsi == NSI_ALGEBRAIC_SPLIT_OSS      ) then
             call nsi_element_assembly_split_oss(                  &
                  pnode,pgaus,gpden,gpvis,gpvis_nsw,gppor,gpsp1,gpsp2,gpvol, &
                  gpsha,gpcar,gpadv,gpvep,gpgrp,gprhs,gprhc,gpvel, &
                  gpgve,gpsgs,elvel,elpre,elbub,elauu,elaup,elapp, &
                  elapu,elrbu,elrbp,dtinv_loc,dtsgs_loc,pbubl,     &
                  gpsha_bub,gpcar_bub,elauq,elapq,elaqu,elaqp,     &
                  elaqq,elrbq,densi,elavv,elporrbu)
          else
             call nsi_element_assembly_asgs_oss(                   &
                  pnode,pgaus,gpden,gpvis,gppor,gpgvi,gpsp1,gptt1, &
                  gpsp2,gptt2,gpvol,gpsha,gpcar,gplap,gphes,gpadv, &
                  gpvep,gprhs,gprhc,rmomu,rcont,elpre,elbub,elauu, &
                  elaup,elapp,elapu,elrbu,elrbp,rmom2,gpst1,gpgve, &
                  gprh2,gprhs_sgs,elvel,ellum,dtinv_loc,pbubl,     &
                  gpsha_bub,gpcar_bub,gppre,elauq,elapq,elaqu,     &
                  elaqp,elaqq,elrbq,densi,gplag,gpmix)
             !call nsi_element_assembly_asgs_oss_old(               &
             !     pnode,pgaus,pevat,gpden,gpvis,gppor,gpgvi,gpsp1, &
             !     gptt1,gpsp2,gptt2,gpvol,gpsha,gpcar,gplap,gphes, &
             !     gpadv,gpvep,gprhs,gprhc,rmomu,rcont,p1vec,p2vec, &
             !     p2sca,wgrgr,wgrvi,elauu,elaup,elapp,elapu,elrbu, &
             !     elrbp,rmom2,p1ve2,gpst1,gpgve,gprh2,gprhs_sgs,   &
             !     elvel,ellum,dtinv_loc,pbubl,gpsha_bub,gpcar_bub, &
             !     gppre)
          end if
          !
          ! Eliminate bubble
          !
          if( maxval(pbubl) == 1 ) then
             call nsi_eliminate_bubble(&
                  pnode,pevat,pbubl,elauu,elaup,elapu,elapp,elrbu,elrbp,&
                  elauq,elapq,elaqu,elaqp,elaqq,elrbq)
          end if

       end if

       call cputim(timeb)
       time1(7) = time1(7) + timeb - timea
       !
       ! Assembly in global system
       !
       if( kfl_cos_opt == 1 ) then
          do ivect = 1,VECTOR_SIZE
             ielem = list_elements(ivect)
             if( ielem > 0 ) then
                pelty = ltype(ielem)
                !
                ! Element properties and dimensions
                !
                if( pelty > 0 ) then
                   !
                   !  calculate cost function F
                   !
                   call nsi_elmcost_all(elvel(ivect,:,:,:),pnode,pgaus,gpsha(ivect,:,:),gpvol(ivect,:),lnods_loc(ivect,:),costf)

                   !if(kfl_adj_prob == 1 ) then
                   !   !
                   !   ! Calculate strong residual &
                   !   !
                   !   call nsi_elmstr(                                              &
                   !        pgaus,pnode,ndime,ielem,lnods(1,ielem),chale,elvel,gpadv,gpvis,gpden, &
                   !        rmomu,rmom2,gprhs,gpvel,gpcar,gpsp1,gpstrm, gpstrc, &
                   !        gpvep,gpgrp,gpst1,gprhs_sgs,gpsha,gpgve)
                   !   !
                   !   ! Modify elmat and elrhs for exact linearization
                   !   ! Send/receive adjoint rhs vectors to/from other modules
                   !   !
                   !   call nsi_element_assembly_der_all(&
                   !        pnode,pgaus,pevat,gpden,gpvis,gppor,gpgvi,gpsp1, &
                   !        gptt1,gpsp2,gptt2,gpvol,gpsha,gpcar,gplap,gphes, &
                   !        gpadv,gpvep,gpprp,gpgrp,elauu,elaup,elapp,elapu,gpst1,&
                   !        gpstrm,gpstrc,elrbu,elrbp,elvel,ielem,gptem,gpgde,&
                   !        gpdvi,elaut,elapt,elpre,gpgve,p1vec,gpgdv,gpdde,gpgdd,chale,&
                   !        dgpmut_dvel,dgpmut_dtur1,dgpmut_dtur2)
                   !   !
                   !   !  calculate cost function derivatives w.r.t. unknowns dF/dU
                   !   !
                   !   call nsi_elmdcost_all(elvel,pnode,pgaus,gpsha,gpvol,lnods(1:pnode,ielem),elrbu,elrbp)
                   !   !
                   !   !  calculate residal derivatives (partial) w.r.t. unknowns dR/dD
                   !   !
                   !   call nsi_elmresdiff_all(pnode,pgaus,pevat,gpvol,gpsha,lnods(1,ielem),&
                   !        p1vec,p2vec,elvel,wgrgr,chale(2),gpsp1,gpsp2,rmomu,gpden,gpadv,gpcar,rcont,resdiff_nsi)
                   !endif
                end if
             end if
          end do
       end if
       !
       ! Shock capturing
       !
       if( kfl_shock_nsi /= 0 ) then
          call nsi_element_shock_capturing(&
               pnode,pgaus,ptopo,pevat,ndime,gpden,gpvel,gprhs,gpsp1,&
               gpvol,elvel,gpcar,chale,rmomu,rmom2,elauu)
       end if
       !
       ! Extension elements
       !
       call nsi_element_extension(&
            1_ip,list_elements,lelch_loc,pnode,elauu,elaup,elapu,elapp,&
            elmap,elrbu,elrbp)
       !
       ! Prescribe Dirichlet boundary conditions
       !
       if( kfl_matdi_nsi == NSI_DIRICHLET_ELEMENT ) then
          call nsi_element_dirichlet(&
               pnode,pevat,list_elements,lnods_loc,&
               elauu,elaup,elapp,elapu,elrbu,elrbp,elcmm)
       end if
       !
       ! Transpose of the jacobian matrix for adjoint
       !
       if( kfl_adj_prob == 1 ) then
          do ivect = 1,VECTOR_SIZE
             elauu(ivect,:,:)      = transpose(elauu(ivect,:,:))
             elauptrans(ivect,:,:) = transpose(elaup(ivect,:,:))
             elaputrans(ivect,:,:) = transpose(elapu(ivect,:,:))
             elapu(ivect,:,:)      = elauptrans(ivect,:,:)
             elaup(ivect,:,:)      = elaputrans(ivect,:,:)
             elapp(ivect,:,:)      = transpose(elapp(ivect,:,:))
          end do
       end if
       !
       ! Schur complement preconditioner element matrix
       !
       if( NSI_SCHUR_COMPLEMENT .or. NSI_FRACTIONAL_STEP ) then
          if( kfl_grad_div_nsi == 0 ) then
             call elmgeo_element_length(&
                  ndime,pnode,porde,tragl,hleng,elcod,elvel(:,:,:,1),chave,chale,&
                  element_type(pelty) % natural_length,kfl_advec_nsi,kfl_ellsh_nsi)
             call nsi_element_schur(&
                  pnode,pgaus,list_elements,lnods_loc,gpcar,gpvol,gpden,gpvis,gppor,&
                  gpsha,elvel,chale,gpsp1,elapp,elmap,dtinv_loc)
          end if
       end if
       !
       ! Assembly in global system
       !
       call cputim(timea)

       if( NSI_MONOLITHIC ) then
          !
          ! Monolithic scheme: A,b
          !
          call nsi_assembly_monolithic(&
               pnode,pevat,list_elements,lnods_loc,elauu,elaup,elapp,elapu,&
               elrbu,elrbp,amatr,rhsid)

       else if( NSI_SCHUR_COMPLEMENT ) then
          !
          ! Monolithic scheme: Auu,Aup,Apu,Aup,bu,bp,Q
          !
          call nsi_assembly_schur_method(&
               1_ip,pnode,pevat,list_elements,lnods_loc,elauu,elaup,elapp,elapu,&
               elrbu,elrbp,amatr(poauu_nsi:),amatr(poaup_nsi:),amatr(poapp_nsi:),&
               amatr(poapu_nsi:),rhsid,rhsid(ndbgs_nsi+1:))
          call nsi_assembly_schur_method(&
               2_ip,pnode,pevat,list_elements,lnods_loc,dumm3,dumm3,elmap,dumm3,&
               dumm3,dumm2,dummr,dummr,lapla_nsi,&
               dummr,dummr,dummr)

       else if( NSI_FRACTIONAL_STEP ) then
          !
          ! Fractional step scheme: Aup,Apu,Aup,bu,bp,Q
          !
          if( kfl_grad_div_nsi /= 0 ) then
             call nsi_assembly_fractional_step_rhs(&
                  pnode,pevat,list_elements,lnods_loc,elvel,elpre,elauu,elaup,elapp,&
                  elapu,elmap,elrbu,elrbp,rhsid,rhsid(ndbgs_nsi+1:))
             
          else
             call nsi_assembly_fractional_step(&
                  pnode,pevat,list_elements,lnods_loc,elvel,elpre,elauu,elaup,elapp,&
                  elapu,elmap,elrbu,elrbp,amatr(poaup_nsi:),amatr(poapu_nsi:),&
                  rhsid,rhsid(ndbgs_nsi+1:),lapla_nsi)
          end if
       end if
       !
       ! Pressure bubble assembly
       !
       if( kfl_bubbl_nsi /= 0 ) then
          call nsi_bubble_assembly(&
               pnode,list_elements,elaqq,elaqu,elaqp,elrbq)
       end if
       !
       ! Lumped matrices
       !
       ! DT_RHO_NSI   = rho / dt  * M
       ! MASS_RHO_NSI = rho * M
       ! TAU_NSI      = tau * M
       !
       if( kfl_predi_nsi == 7 .or. kfl_predi_nsi == 8 .or. kfl_predi_nsi == 9 ) then
          call nsi_assembly_dt_rho_tau_nu(&
               pnode,list_elements,lnods_loc,eldtrho,eltau,elmurho,&
               dt_rho_nsi,tau_nsi,mass_rho_nsi)
       end if
       !
       ! Dual time step preconditioner
       !
       if( kfl_duatss == 1 ) then
          do ivect = 1,VECTOR_SIZE
             ielem = list_elements(ivect)
             if( ielem > 0 ) then
                do inode = 1,pnode
                   ipoin = lnods_loc(ivect,inode)
#ifdef NO_COLORING
                   !$OMP ATOMIC
#endif
                   lumma(ipoin) = lumma(ipoin) + ellum(ivect,inode)
                end do
             end if
          end do
       end if

       call cputim(timeb)
       time1(8) = time1(8) + timeb - timea

    end if

  end subroutine nsi_element_operations

  subroutine nsi_element_operations_gather_scalar(&
       pnode,pgaus,ielem,lnods,elcod,elpre,elvel,elfle,elvep,elprp,&
       elgrp,eltem,elmsh,elnor,elcur,elwmean,gphyd,&
       gpsgt,gpsgs,elbub,elden,elavv,elibopo,elnnsw,elywal)
    !-----------------------------------------------------------------------
    !****f* Nastin/nsi_elmga3
    ! NAME
    !    nsi_elmga3
    ! DESCRIPTION
    !    Compute some variables at the Gauss points
    !    ELVEL, ELCOD, ELPRE
    ! OUTPUT
    ! USES
    ! USED BY
    !***
    !-----------------------------------------------------------------------
    use def_kintyp, only       :  ip,rp
    use def_kermod, only       :  kfl_adj_prob,gasco
    use def_master, only       :  veloc,press,tempe,fleve,&
         &                        velom,kfl_coupl,wmean,&
         &                        ID_NASTIN,ID_ALEFOR,ID_TEMPER,&
         &                        ID_LEVELS,ID_CHEMIC,veloc_forw,&
         &                        tempe_forw,tesgs,vesgs,prthe
    use def_domain, only       :  ndime,coord,ywale
    use def_kermod, only       :  kfl_noslw_ker,avupo_ker
    use def_nastin, only       :  kfl_timei_nsi,kfl_regim_nsi,nbdfp_nsi,&
         &                        curle_nsi,norle_nsi,kfl_colev_nsi,&
         &                        kfl_stabi_nsi,vepro_nsi,prpro_nsi,&
         &                        grpro_nsi,kfl_cotem_nsi,&
         &                        nbdfp_nsi,kfl_surte_nsi,kfl_sgste_nsi,&
         &                        hydro_density_nsi,kfl_hydro_gravity_nsi,&
         &                        kfl_bubbl_nsi,bubble_nsi
    implicit none
    integer(ip), intent(in)    :: pnode
    integer(ip), intent(in)    :: pgaus
    integer(ip), intent(in)    :: ielem
    integer(ip), intent(in)    :: lnods(pnode)
    real(rp),    intent(out)   :: elcod(ndime,pnode)
    real(rp),    intent(out)   :: elpre(pnode,*)
    real(rp),    intent(out)   :: elvel(ndime,pnode,*)
    real(rp),    intent(out)   :: elfle(pnode)
    real(rp),    intent(out)   :: elvep(ndime,pnode)
    real(rp),    intent(out)   :: elprp(pnode,*)
    real(rp),    intent(out)   :: elgrp(ndime,pnode)
    real(rp),    intent(out)   :: eltem(pnode,nbdfp_nsi)
    real(rp),    intent(out)   :: elwmean(pnode,nbdfp_nsi)
    real(rp),    intent(out)   :: elmsh(ndime,pnode)
    real(rp),    intent(out)   :: elnor(ndime,pnode)
    real(rp),    intent(out)   :: elcur(pnode)
    real(rp),    intent(out)   :: gphyd(pgaus)
    real(rp),    intent(out)   :: gpsgt(pgaus,nbdfp_nsi)
    real(rp),    intent(out)   :: gpsgs(ndime,pgaus,2)
    real(rp),    intent(out)   :: elbub
    real(rp),    intent(out)   :: elden(pnode)
    real(rp),    intent(out)   :: elavv(ndime,pnode)
    real(rp),    intent(out)   :: elibopo(pnode)
    real(rp),    intent(out)   :: elnnsw(ndime)
    real(rp),    intent(out)   :: elywal

    integer(ip)                :: inode,idime,ipoin,itime

    if( kfl_timei_nsi == 0 ) then
       !
       ! Stationary
       !
       do inode = 1,pnode
          ipoin = lnods(inode)
          do idime = 1,ndime
             if (kfl_adj_prob == 0) then
                elvel(idime,inode,1) = veloc(idime,ipoin,1)
                elvel(idime,inode,2) = 0.0_rp
             else
                elvel(idime,inode,1) = veloc_forw(idime,ipoin,1)
                elvel(idime,inode,2) = 0.0_rp
             endif
             elcod(idime,inode)   = coord(idime,ipoin)
          end do
          elpre(inode,1) = press(ipoin,1)
       end do
    else
       !
       ! Transient
       !
       do inode = 1,pnode
          ipoin = lnods(inode)
          do idime = 1,ndime
             if (kfl_adj_prob == 0) then
                elvel(idime,inode,1) = veloc(idime,ipoin,1)
                elvel(idime,inode,2) = veloc(idime,ipoin,3)
                do itime=3,nbdfp_nsi
                   elvel(idime,inode,itime) = veloc(idime,ipoin,itime+1)
                end do
             else
                elvel(idime,inode,1) = veloc_forw(idime,ipoin,1)
                elvel(idime,inode,2) = veloc(idime,ipoin,3)
                do itime=3,nbdfp_nsi
                   elvel(idime,inode,itime) = veloc(idime,ipoin,itime+1)
                end do
             endif
             elcod(idime,inode)   = coord(idime,ipoin)
          end do
          elpre(inode,1) = press(ipoin,1)
       end do
    end if
    !
    ! Surface tension
    !
    if( kfl_colev_nsi /= 0 ) then
       elfle(1:pnode) = fleve(lnods(1:pnode),1)
    end if
    !
    ! Surface tension
    !
    if( kfl_surte_nsi /= 0 ) then
       do inode = 1,pnode
          ipoin = lnods(inode)
          elcur(inode) = curle_nsi(ipoin)
          do idime =1,ndime
             elnor(1:ndime,inode) = norle_nsi(1:ndime,ipoin)
          end do
       end do
    end if
    !
    ! Projections
    !
    if( kfl_stabi_nsi > 0 ) then
       do inode = 1,pnode
          ipoin = lnods(inode)
          do idime = 1,ndime
             elvep(idime,inode) = vepro_nsi(idime,ipoin)
          end do
          elprp(inode,1) = prpro_nsi(ipoin)
       end do
       if( kfl_stabi_nsi == 2 ) then
          do inode = 1,pnode
             ipoin = lnods(inode)
             do idime = 1,ndime
                elgrp(idime,inode) = grpro_nsi(idime,ipoin)
             end do
          end do
       end if
    end if
    !
    ! Coupling with Temper
    !
    if( ( kfl_coupl(ID_NASTIN,ID_TEMPER) /= 0 .or. kfl_cotem_nsi == 1 .or. kfl_regim_nsi==3 ) ) then
       do inode = 1,pnode
          ipoin = lnods(inode)
          if( kfl_coupl(ID_NASTIN,ID_CHEMIC) /= 0) then
             elwmean(inode,1) = wmean(ipoin,1)
          else
             elwmean(inode,1) = 1.0_rp
          endif
          if (kfl_adj_prob == 0) then
             eltem(inode,1) = tempe(ipoin,1)
          else
             eltem(inode,1) = tempe_forw(ipoin,1)
          endif
       end do
       if (kfl_regim_nsi == 3) then ! Low Mach needs time information
          do itime=2,nbdfp_nsi
             do inode = 1,pnode
                ipoin = lnods(inode)
                if( kfl_coupl(ID_NASTIN,ID_CHEMIC) /= 0) then
                   elwmean(inode,itime) = wmean(ipoin,itime+1)
                else
                   elwmean(inode,itime) = 1.0_rp
                endif
                if (kfl_adj_prob == 0) then
                   eltem(inode,itime) = tempe(ipoin,itime+1)
                else
                   eltem(inode,itime) = tempe_forw(ipoin,1)
                endif
             end do
          end do

          ! Oriol: test for lowMach
          do inode = 1,pnode
             elden(inode) = (prthe(1)/gasco) * (elwmean(inode,1)/eltem(inode,1))
          end do
       endif
    end if
    !
    ! Mesh velocity
    !
    if( associated(velom) ) then ! TESTEO
       do inode = 1,pnode
          ipoin = lnods(inode)
          elmsh(1:ndime,inode) = velom(1:ndime,ipoin)
       end do
    end if
    !
    ! rho of the hydrostatic state when couplign with level set
    !
    if( kfl_hydro_gravity_nsi /= 0 ) then
       gphyd(1:pgaus) = hydro_density_nsi(ielem) % a(1:pgaus)
    else
       gphyd = 0.0_rp
    end if
    !
    ! Subgrid scales
    !
    call nsi_subgrid_scale_gather(ndime,pgaus,ielem,vesgs,gpsgs(:,:,:))
    if( kfl_regim_nsi == 3 .and. associated(tesgs) ) then
       gpsgt(1:pgaus,1:nbdfp_nsi) = tesgs(ielem) % a(1,1:pgaus,1:nbdfp_nsi)
    else
       if( kfl_sgste_nsi == 1 ) then
          gpsgt(1:pgaus,1) = tesgs(ielem) % a(1,1:pgaus,1)
       else
          gpsgt(1:pgaus,1) = 0.0_rp
       end if
    end if
    !
    ! Bubble
    !
    if( kfl_bubbl_nsi /= 0 ) then
       elbub = bubble_nsi(ielem)
    end if

    !
    ! no slip wall law - - I could add flag so that i does it only on those elements where it is needed kfl_nswel_ker(ielem)
    !
    if ( kfl_noslw_ker /= 0_ip ) then
       do inode = 1,pnode
          ipoin = lnods(inode)
          if (lpoty(ipoin) > 0 ) then
             elibopo(inode) = 1.0_rp
          else
             elibopo(inode) = 0.0_rp
          end if
          do idime = 1,ndime
             elavv(idime,inode) = avupo_ker(idime,ipoin)
          end do
       end do
       if (kfl_nswel_ker(ielem) > 0_ip ) then
          do idime = 1,ndime
             elnnsw(idime) = normal_nsw_ker(idime,kfl_nswel_ker(ielem)) 
          end do
       else
          elnnsw = 0.0_rp
       end if
       elywal = ywale(ielem)

    end if

  end subroutine nsi_element_operations_gather_scalar

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    6/10/2016
  !> @brief   Gather
  !> @details Gather arrays from nesh to element
  !>
  !-----------------------------------------------------------------------
  
  subroutine nsi_element_operations_gather_vector(&
       pnode,pgaus,list_elements,lnods,elcod,elpre,elvel,elfle,&
       elvep,elprp,elgrp,eltem,elmsh,elnor,elcur,elwmean,&
       gphyd,gpsgt,gpsgs,elbub,ellag,elden,elavv,elibopo,elnnsw,&
       elywal)

    use def_kintyp, only       :  ip,rp
    use def_kermod, only       :  kfl_adj_prob
!    use def_kermod, only       :  gasco
    use def_master, only       :  veloc,press,tempe,fleve,&
         &                        velom,kfl_coupl,wmean,&
         &                        ID_NASTIN,ID_ALEFOR,ID_TEMPER,&
         &                        ID_LEVELS,ID_CHEMIC,veloc_forw,&
         &                        tempe_forw,tesgs,vesgs
!    use def_master, only       :  prthe
    use def_domain, only       :  ndime,coord,ywale
    use def_kermod, only       :  kfl_noslw_ker,avupo_ker, kfl_delta
    use def_nastin, only       :  kfl_timei_nsi,kfl_regim_nsi,nbdfp_nsi,&
         &                        curle_nsi,norle_nsi,kfl_colev_nsi,&
         &                        kfl_stabi_nsi,vepro_nsi,prpro_nsi,&
         &                        grpro_nsi,kfl_cotem_nsi,&
         &                        nbdfp_nsi,kfl_surte_nsi,kfl_sgste_nsi,&
         &                        hydro_density_nsi,kfl_hydro_gravity_nsi,&
         &                        ncomp_nsi,bubble_nsi,kfl_bubbl_nsi,&
         &                        kfl_immer_nsi,lagra_nsi,tauib_nsi,&
         &                        kfl_lookg_nsi
    use mod_ker_proper, only   :  ker_proper

    use def_kermod,     only   :  delta_dom
    implicit none
    integer(ip), intent(in)    :: pnode
    integer(ip), intent(in)    :: pgaus
    integer(ip), intent(in)    :: list_elements(VECTOR_SIZE)
    integer(ip), intent(in)    :: lnods(VECTOR_SIZE,pnode)
    real(rp),    intent(out)   :: elcod(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(out)   :: elpre(VECTOR_SIZE,pnode,ncomp_nsi-1)
    real(rp),    intent(out)   :: elvel(VECTOR_SIZE,ndime,pnode,ncomp_nsi+1)
    real(rp),    intent(out)   :: elfle(VECTOR_SIZE,pnode)
    real(rp),    intent(out)   :: elvep(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(out)   :: elprp(VECTOR_SIZE,pnode)
    real(rp),    intent(out)   :: elgrp(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(out)   :: eltem(VECTOR_SIZE,pnode,nbdfp_nsi)
    real(rp),    intent(out)   :: elwmean(VECTOR_SIZE,pnode,nbdfp_nsi)
    real(rp),    intent(out)   :: elmsh(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(out)   :: elnor(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(out)   :: elcur(VECTOR_SIZE,pnode)
    real(rp),    intent(out)   :: gphyd(VECTOR_SIZE,pgaus)
    real(rp),    intent(out)   :: gpsgt(VECTOR_SIZE,pgaus,nbdfp_nsi)
    real(rp),    intent(out)   :: gpsgs(VECTOR_SIZE,ndime,pgaus,2)
    real(rp),    intent(out)   :: elbub(VECTOR_SIZE)
    real(rp),    intent(out)   :: ellag(VECTOR_SIZE,ndime,ndime,pnode)
    real(rp),    intent(out)   :: elden(VECTOR_SIZE,pnode)
    real(rp),    intent(out)   :: elavv(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(out)   :: elibopo(VECTOR_SIZE,pnode)
    real(rp),    intent(out)   :: elnnsw(VECTOR_SIZE,ndime)
    real(rp),    intent(out)   :: elywal(VECTOR_SIZE)
    
    integer(ip)                :: inode,idime,ipoin,itime,ivect,ielem,dummi

    do ivect = 1,VECTOR_SIZE
       ielem = list_elements(ivect)
       if( ielem > 0 ) then
          if( kfl_timei_nsi == 0 ) then
             !
             ! Stationary
             !
             do inode = 1,pnode
                ipoin = lnods(ivect,inode)
                do idime = 1,ndime
                   if (kfl_adj_prob == 0) then
                      elvel(ivect,idime,inode,1) = veloc(idime,ipoin,1)
                      elvel(ivect,idime,inode,2) = 0.0_rp
                   else
                      elvel(ivect,idime,inode,1) = veloc_forw(idime,ipoin,1)
                      elvel(ivect,idime,inode,2) = 0.0_rp
                   endif
                   elcod(ivect,idime,inode)   = coord(idime,ipoin)
                end do
                elpre(ivect,inode,1) = press(ipoin,1)
             end do
          else
             !
             ! Transient
             !
             do inode = 1,pnode
                ipoin = lnods(ivect,inode)
                do idime = 1,ndime
                   if (kfl_adj_prob == 0) then
                      elvel(ivect,idime,inode,1) = veloc(idime,ipoin,1)
                      elvel(ivect,idime,inode,2) = veloc(idime,ipoin,3)
                      do itime=3,nbdfp_nsi
                         elvel(ivect,idime,inode,itime) = veloc(idime,ipoin,itime+1)
                      end do
                      if  (kfl_stabi_nsi == NSI_ASGS.and.NSI_FRACTIONAL_STEP) &
                           ! need and extra velocity for the derivative in the residual
                           elvel(ivect,idime,inode,3) = veloc(idime,ipoin,4)
                   else
                      elvel(ivect,idime,inode,1) = veloc_forw(idime,ipoin,1)
                      elvel(ivect,idime,inode,2) = veloc(idime,ipoin,3)
                      do itime=3,nbdfp_nsi
                         elvel(ivect,idime,inode,itime) = veloc(idime,ipoin,itime+1)
                      end do
                   endif
                   elcod(ivect,idime,inode)   = coord(idime,ipoin)
                end do
                elpre(ivect,inode,1) = press(ipoin,1)
             end do
          end if
          !
          ! Surface tension
          !
          if( kfl_colev_nsi /= 0 ) then
             elfle(ivect,1:pnode) = fleve(lnods(ivect,1:pnode),1)
          else
             elfle(ivect,1:pnode) = 0.0_rp
          end if
          !
          ! Surface tension
          !
          if( kfl_surte_nsi /= 0 ) then
             do inode = 1,pnode
                ipoin = lnods(ivect,inode)
                elcur(ivect,inode) = curle_nsi(ipoin)
                do idime =1,ndime
                   elnor(ivect,1:ndime,inode) = norle_nsi(1:ndime,ipoin)
                end do
             end do
          end if
          !
          ! Projections
          !
          if( kfl_stabi_nsi > 0 ) then
             do inode = 1,pnode
                ipoin = lnods(ivect,inode)
                elvep(ivect,1:ndime,inode) = vepro_nsi(1:ndime,ipoin)
                elprp(ivect,inode)         = prpro_nsi(ipoin)
             end do
             if( kfl_stabi_nsi == 2 ) then
                do inode = 1,pnode
                   ipoin = lnods(ivect,inode)
                   elgrp(ivect,1:ndime,inode) = grpro_nsi(1:ndime,ipoin)
                end do
             end if
          end if
          !
          ! Coupling with Temper
          !
          if( ( kfl_coupl(ID_NASTIN,ID_TEMPER) /= 0 .or. kfl_cotem_nsi == 1 .or. kfl_regim_nsi==3 ) ) then

             if(kfl_surte_nsi == 2_ip) then
                call ker_proper('DENSI','PNODE',dummi,ielem,elden(ivect,:),pnode)
             else
                eltem(ivect,:,:)   = 0.0_rp
                elwmean(ivect,:,:) = 0.0_rp

                do inode = 1,pnode
                   ipoin = lnods(ivect,inode)
                   if( kfl_coupl(ID_NASTIN,ID_CHEMIC) /= 0 .and. associated(wmean)) then
                      elwmean(ivect,inode,1) = wmean(ipoin,1)
                   else
                      elwmean(ivect,inode,1) = 1.0_rp
                   endif
                   if (associated(tempe)) then
                      if (kfl_adj_prob == 0) then
                         eltem(ivect,inode,1) = tempe(ipoin,1)
                      else
                         eltem(ivect,inode,1) = tempe_forw(ipoin,1)
                      endif
                   endif
                end do

                if (kfl_regim_nsi == 3 ) then ! (not in for fractional step)

                   !do itime=2,nbdfp_nsi
                   !   do inode = 1,pnode
                   !      ipoin = lnods(ivect,inode)
                   !      if( kfl_coupl(ID_NASTIN,ID_CHEMIC) /= 0) then
                   !         elwmean(ivect,inode,itime) = wmean(ipoin,itime+1)
                   !      else
                   !         elwmean(ivect,inode,itime) = 1.0_rp
                   !      endif
                   !      if (kfl_adj_prob == 0) then
                   !         eltem(ivect,inode,itime) = tempe(ipoin,itime+1)
                   !      else
                   !         eltem(ivect,inode,itime) = tempe_forw(ipoin,1)
                   !      endif
                   !   end do
                   !end do


                   if (kfl_lookg_nsi >0) then
                      call ker_proper('DENSI','PNODE',dummi,ielem,elden(ivect,:),pnode)
                   else
                      call ker_proper('DENSI','PNODE',dummi,ielem,elden(ivect,:),pnode)
                      !! Oriol: test for lowMach
                      !do inode = 1,pnode
                      !   elden(ivect,inode) = (prthe(1)/gasco) * (elwmean(ivect,inode,1)/eltem(ivect,inode,1))
                      !end do
                   end if
                endif
             end if
          end if
          !
          ! Mesh velocity
          !
          if( associated(velom) ) then ! TESTEO
             do inode = 1,pnode
                ipoin = lnods(ivect,inode)
                elmsh(ivect,1:ndime,inode) = velom(1:ndime,ipoin)
             end do
          end if
          !
          ! rho of the hydrostatic state when couplign with level set
          !
          if( kfl_hydro_gravity_nsi /= 0 ) then
             gphyd(ivect,1:pgaus) = hydro_density_nsi(ielem) % a(1:pgaus)
          else
             gphyd(ivect,:) = 0.0_rp
          end if
          !
          ! Subgrid scales
          !
          call nsi_subgrid_scale_gather(ndime,pgaus,ielem,vesgs,gpsgs(ivect,:,:,:))
          if( kfl_regim_nsi == 3 .and. associated(tesgs) ) then
             gpsgt(ivect,1:pgaus,1:nbdfp_nsi) = tesgs(ielem) % a(1,1:pgaus,1:nbdfp_nsi)
          else
             if( kfl_sgste_nsi == 1 ) then
                gpsgt(ivect,1:pgaus,1) = tesgs(ielem) % a(1,1:pgaus,1)
             else
                gpsgt(ivect,1:pgaus,1) = 0.0_rp
             end if
          end if
          !
          ! Bubble
          !
          if( kfl_bubbl_nsi /= 0 ) then
             elbub(ivect) = bubble_nsi(ielem)
          else
             elbub(ivect) = 0.0_rp
          end if
          !
          ! Immersed boundary
          !
          if(    kfl_immer_nsi == 1 ) then
             do inode = 1,pnode
                ipoin = lnods(ivect,inode)
                ellag(ivect,1:ndime,1,inode) = lagra_nsi(1:ndime,ipoin,1)
             end do
          else if( kfl_immer_nsi == 2 ) then
             do inode = 1,pnode
                ipoin = lnods(ivect,inode)
                ellag(ivect,1:ndime,1:ndime,inode) = tauib_nsi(1:ndime,1:ndime,ipoin)
             end do
          end if
          !
          ! no slip wall law - I could add flag so that i does it only on those elements where it is needed kfl_nswel_ker(ielem)
          !
          if ( kfl_noslw_ker /= 0_ip ) then
             do inode = 1,pnode
                ipoin = lnods(ivect,inode)
                if (lpoty(ipoin) > 0 ) then
                   elibopo(ivect,inode) = 1.0_rp
                else
                   elibopo(ivect,inode) = 0.0_rp
                end if
                do idime = 1,ndime
                   elavv(ivect,idime,inode) = avupo_ker(idime,ipoin)
                end do
             end do



             
             if (kfl_nswel_ker(ielem) > 0_ip ) then
                do idime = 1,ndime
                   elnnsw(ivect,idime) = normal_nsw_ker(idime,kfl_nswel_ker(ielem)) 
                end do
             else
                elnnsw(ivect,:) = 0.0_rp
             end if
             if (kfl_delta == 0) then
                elywal(ivect) = delta_dom
             else
                elywal(ivect) = ywale(ielem)
             end if
          end if
       else
          !
          ! Element number is null
          !
          elcod(ivect,:,:)   = 0.0_rp
          elpre(ivect,:,:)   = 0.0_rp
          elvel(ivect,:,:,:) = 0.0_rp
          elfle(ivect,:)     = 0.0_rp
          elvep(ivect,:,:)   = 0.0_rp
          elprp(ivect,:)     = 0.0_rp
          elgrp(ivect,:,:)   = 0.0_rp
          eltem(ivect,:,:)   = 0.0_rp
          elmsh(ivect,:,:)   = 0.0_rp
          elnor(ivect,:,:)   = 0.0_rp
          elcur(ivect,:)     = 0.0_rp
          elwmean(ivect,:,:) = 0.0_rp
          gphyd(ivect,:)     = 0.0_rp
          gpsgt(ivect,:,:)   = 0.0_rp
          gpsgs(ivect,:,:,:) = 0.0_rp
          elbub(ivect)       = 0.0_rp
          ellag(ivect,:,:,:) = 0.0_rp
          elden(ivect,:)     = 0.0_rp
          elavv(ivect,:,:)   = 0.0_rp
          elnnsw(ivect,:)    = 0.0_rp
          elibopo(ivect,:)   = 0.0_rp
       end if

    end do

  end subroutine nsi_element_operations_gather_vector

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    6/10/2016
  !> @brief   Compute the residual of the momentum and continuity equations
  !> @details Compute the residual of the momentum and continuity equations:
  !>          rho/(dt*theta)*u + sig*u - mu*Lap(u)
  !>          rho*(uc.grad)u is computed in nsi_elmsgs as uc can depend on the
  !>          subgrid scale
  !>
  !>          Three forms of the viscous term are considered.
  !>          1. Laplacian form:  div[mu*grad(u)]
  !>          2. Divergence form: div[2*mu*eps(u)]
  !>          3. Complete form:   div[2*mu*epsi(u)]
  !>          eps(u)  = 1/2[grad(u)+grad(u)^t]
  !>          eps'(u) = 1/2[grad(u)+grad(u)^t]-1/3*div(u)I
  !>
  !>          <------ Laplacian
  !>          <------------------------ Divergence
  !>          <------------------------------------------------- Complete
  !>          -mu*(d^2ui/dxj^2)-mu*(d^2uj/dxj*dxi)+2/3*mu*d(div(u))/dx
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_element_residual(                              &
       pnode,pgaus,plapl,gpsha,gpcar,gphes,gpgvi,gpden,gpvis,   &
       gppor,gpapo,gptem,gpsgs,elvel,elpre,elvep,elprp,elgrp,ellag,   &
       eltem,elmsh,elcod,elnor,elcur,elbub,elwmean,hleng,chale, &
       gpvel,gpgpr,rmomu,rmom2,rcont,gprhs,gprhc,gplap,gpadv,   &
       gpvep,gpprp,gpgrp,gphyd,gpmsh,gpgve,gpnor,gpcur,gpfle,   &
       ielem,gprh2,gppre,gprhs_sgs,dtinv_loc,gpgde,gpsgt,       &
       gpsha_bub,cutim,densi,gplag,elden,list_elements)

    use def_kintyp, only     :  ip,rp
    use def_parame, only     :  pi
    use def_master, only     :  kfl_coupl,&
         &                      ID_NASTIN,ID_ALEFOR,ID_TEMPER,   &
         &                      ID_CHEMIC,tesgs,kfl_lumped,            &
         &                      prthe,velom
!    use def_master, only     :  densi_gp
    use def_kermod, only     :  gasco,kfl_adj_prob
!    use def_kermod, only     :  thicl
    use def_domain, only     :  ndime,ntens,mnode
    use def_nastin, only     :  grnor_nsi,gravi_nsi,     &
         &                      kfl_sgsti_nsi,bougr_nsi,boube_nsi,     &
         &                      boutr_nsi,kfl_cotem_nsi,kfl_advec_nsi, &
         &                      dtsgs_nsi,kfl_stabi_nsi,kfl_sgsco_nsi, &
         &                      kfl_hydro_gravity_nsi,kfl_force_nsi,   &
         &                      nbdfp_nsi,pabdf_nsi,                   &
         &                      corio_nsi,facca_nsi,fvela_nsi,         &
         &                      kfl_regim_nsi,                         &
         &                      kfl_prthe_nsi,                         &
         &                      faccl_nsi,frotc_nsi,fvins_nsi,         &
         &                      kfl_surte_nsi,surte_nsi,kfl_rmom2_nsi, &
         &                      centr_nsi,gravb_nsi,                   &
         &                      kfl_bubbl_nsi,                         &
         &                      kfl_immer_nsi,                         &
         &                      kfl_convection_type_nsi,               &
         &                      NSI_CONVECTION_EMAC,                   &
         &                      NSI_CONVECTION_SKEW,                   &
         &                      kfl_lookg_nsi

!    use def_nastin, only     :  lowtr_nsi
    !
    ! AB: densi
    !
    use def_master, only     :  wmean_gp,tempe_gp,prthe
    use mod_ker_proper, only  : ker_proper
    

    implicit none
    integer(ip), intent(in)    :: pnode,pgaus,plapl,ielem
    integer(ip), intent(in)    :: list_elements(VECTOR_SIZE)
    real(rp),    intent(in)    :: gpsha(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(in)    :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)
    real(rp),    intent(in)    :: gphes(VECTOR_SIZE,ntens,mnode,pgaus)
    real(rp),    intent(in)    :: gpgvi(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: gpden(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpvis(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gppor(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpapo(VECTOR_SIZE,ndime,ndime,pgaus)
    real(rp),    intent(out)   :: gptem(VECTOR_SIZE,pgaus,nbdfp_nsi)
    real(rp),    intent(in)    :: gpsgs(VECTOR_SIZE,ndime,pgaus,*)
    real(rp),    intent(in)    :: gpgde(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: elvel(VECTOR_SIZE,ndime,pnode,nbdfp_nsi)
    real(rp),    intent(in)    :: elpre(VECTOR_SIZE,pnode)
    real(rp),    intent(in)    :: elvep(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(in)    :: elprp(VECTOR_SIZE,pnode)
    real(rp),    intent(in)    :: elgrp(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(in)    :: ellag(VECTOR_SIZE,ndime,ndime,pnode)
    real(rp),    intent(in)    :: eltem(VECTOR_SIZE,pnode,nbdfp_nsi)
    real(rp),    intent(in)    :: elwmean(VECTOR_SIZE,pnode,nbdfp_nsi)
    real(rp),    intent(in)    :: elmsh(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(in)    :: elcod(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(in)    :: elnor(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(in)    :: elcur(VECTOR_SIZE,pnode)
    real(rp),    intent(in)    :: elden(VECTOR_SIZE,pnode)
    real(rp),    intent(in)    :: elbub(VECTOR_SIZE)
    real(rp),    intent(in)    :: hleng(VECTOR_SIZE,ndime)
    real(rp),    intent(in)    :: chale(VECTOR_SIZE,2)
    real(rp),    intent(in)    :: dtinv_loc(VECTOR_SIZE)
    real(rp),    intent(out)   :: gpvel(VECTOR_SIZE,ndime,pgaus,nbdfp_nsi)
    real(rp),    intent(out)   :: gpgpr(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(out)   :: rmomu(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(out)   :: rmom2(VECTOR_SIZE,ndime,ndime,pnode,pgaus)
    real(rp),    intent(out)   :: rcont(VECTOR_SIZE,ndime,pnode,pgaus)
    real(rp),    intent(out)   :: gprhs(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(out)   :: gprhs_sgs(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(out)   :: gprhc(VECTOR_SIZE,pgaus)
    real(rp),    intent(out)   :: gprh2(VECTOR_SIZE,pgaus)
    real(rp),    intent(out)   :: gplap(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(out)   :: gpadv(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(out)   :: gpvep(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(out)   :: gpprp(VECTOR_SIZE,pgaus)
    real(rp),    intent(out)   :: gpgrp(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(inout) :: gphyd(VECTOR_SIZE,pgaus)
    real(rp),    intent(out)   :: gpmsh(VECTOR_SIZE,ndime,pgaus)           ! Mehs velocity
    real(rp),    intent(out)   :: gpgve(VECTOR_SIZE,ndime,ndime,pgaus)     ! Velocity gradient
    real(rp),    intent(out)   :: gpnor(VECTOR_SIZE,ndime,pgaus)           ! LS normal
    real(rp),    intent(out)   :: gpcur(VECTOR_SIZE,pgaus)                 ! LS curvature
    real(rp),    intent(in)    :: gpfle(VECTOR_SIZE,pgaus)                 ! LS - already calculated in nsi_proper
    real(rp),    intent(out)   :: gppre(VECTOR_SIZE,pgaus)                 ! pressure
    real(rp),    intent(in)    :: gpsgt(VECTOR_SIZE,pgaus,*)               ! LS - already calculated in nsi_proper
    real(rp),    intent(out)   :: gpsha_bub(VECTOR_SIZE,pgaus)             ! bubble shape function
    real(rp),    intent(in)    :: cutim                                    ! Current time
    real(rp),    intent(out)   :: densi(VECTOR_SIZE,pgaus,nbdfp_nsi)
    real(rp),    intent(out)   :: gplag(VECTOR_SIZE,ndime,ndime,pgaus)     ! Lagrange multiplier
    integer(ip)                :: idime,inode,igaus,jdime
    integer(ip)                :: itime,ivect
    real(rp)                   :: fact0(1:VECTOR_SIZE),fact1(1:VECTOR_SIZE)
    real(rp)                   :: fact2(1:VECTOR_SIZE),fact3(1:VECTOR_SIZE)
    real(rp)                   :: w1
    real(rp)                   :: gpext(1:VECTOR_SIZE,3)                   ! External force
    real(rp)                   :: gpcod(1:VECTOR_SIZE,3)                   ! Axes motion
    real(rp)                   :: alpha(1:VECTOR_SIZE,3)                   ! Axes motion
    real(rp)                   :: dummr(1:VECTOR_SIZE,3)                   ! Axes motion
    real(rp)                   :: centf(1:VECTOR_SIZE,3)                   ! Axes motion
    real(rp)                   :: gpgrt(1:VECTOR_SIZE,pgaus)               ! Low mach intermmediate vars
    real(rp)                   :: gpdet(1:VECTOR_SIZE,pgaus)               ! Low mach intermmediate vars
#ifdef matiaslma
    real(rp)                   :: gpgte(1:VECTOR_SIZE,ndime,pgaus)         ! Low mach intermmediate vars
#endif
    real(rp)                   :: xvis2
    real(rp)                   :: gpwmean(1:VECTOR_SIZE,pgaus,nbdfp_nsi)
    real(rp)                   :: dtinv_res(VECTOR_SIZE),  dtinv_mod(VECTOR_SIZE)

    real(rp)                         :: gpden2(pgaus)

    integer(ip) :: ivect_loc, ielem_loc, dummi

#ifdef GABLS1TREF
    real(rp)    :: boutr_aux(VECTOR_SIZE)             ! reference temperature
    real(rp)    :: twall,ttop
#endif

#ifdef OPENACC
#define DEF_VECT ivect
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif

    !----------------------------------------------------------------------
    !
    ! Time step
    !
    !----------------------------------------------------------------------

    if( kfl_lumped == 2 ) then
       dtinv_res = 0.0_rp
    else
       dtinv_res = dtinv_loc
    end if

    gprhs     = 0.0_rp
    gprhs_sgs = 0.0_rp
    gprhc     = 0.0_rp
    gprh2     = 0.0_rp

#ifdef OPENACC
    do ivect = 1,VECTOR_SIZE
#endif

       !----------------------------------------------------------------------
       !
       ! Gauss point values
       !
       !----------------------------------------------------------------------
       !
       ! GPVEL, GPGPR: u, grad(p)
       !
       gpvel(DEF_VECT,:,:,:) = 0.0_rp
       gpgpr(DEF_VECT,:,:)   = 0.0_rp

       do igaus = 1,pgaus
          do inode = 1,pnode
             do idime = 1,ndime
                gpvel(DEF_VECT,idime,igaus,1) = gpvel(DEF_VECT,idime,igaus,1) + elvel(DEF_VECT,idime,inode,1) * gpsha(DEF_VECT,inode,igaus)
                gpvel(DEF_VECT,idime,igaus,2) = gpvel(DEF_VECT,idime,igaus,2) + elvel(DEF_VECT,idime,inode,2) * gpsha(DEF_VECT,inode,igaus)
                do itime=3,nbdfp_nsi
                   gpvel(DEF_VECT,idime,igaus,itime) = gpvel(DEF_VECT,idime,igaus,itime) + elvel(DEF_VECT,idime,inode,itime) * gpsha(DEF_VECT,inode,igaus)
                end do
                gpgpr(DEF_VECT,idime,igaus) = gpgpr(DEF_VECT,idime,igaus) + elpre(DEF_VECT,inode) * gpcar(DEF_VECT,idime,inode,igaus)
             end do
          end do
       end do
       !
       ! Pressure
       !
       if( kfl_regim_nsi == 3 ) then
          do igaus = 1,pgaus
             gppre(DEF_VECT,igaus) = 0.0_rp
             do inode = 1,pnode
                gppre(DEF_VECT,igaus) = gppre(DEF_VECT,igaus) + elpre(DEF_VECT,inode) * gpsha(DEF_VECT,inode,igaus)
             end do
          end do
          if( kfl_bubbl_nsi /= 0 ) then
             do igaus = 1,pgaus
                gppre(DEF_VECT,igaus) = gppre(DEF_VECT,igaus) + elbub(DEF_VECT) * gpsha_bub(DEF_VECT,igaus)
             end do
          end if
       end if
       !
       ! GPADV: Advection velocity
       !
       if( kfl_advec_nsi /= 0 ) then
          gpadv(:,:,:) = gpvel(:,:,:,1)
       else
          gpadv = 0.0_rp
       end if
       !
       ! ALE
       !
       if( associated(velom) ) then ! TESTEO
          do igaus = 1,pgaus
             gpmsh(DEF_VECT,1:ndime,igaus) = 0.0_rp
             do inode = 1,pnode
                do idime = 1,ndime
                   gpmsh(DEF_VECT,idime,igaus) = gpmsh(DEF_VECT,idime,igaus) + elmsh(DEF_VECT,idime,inode) * gpsha(DEF_VECT,inode,igaus)
                end do
             end do
             do idime = 1,ndime
                gpadv(DEF_VECT,idime,igaus) = gpadv(DEF_VECT,idime,igaus) - gpmsh(DEF_VECT,idime,igaus)
             end do
          end do
       end if
       !
       ! Convection tracking: uc = u + u'
       !
       if( kfl_sgsco_nsi >= 1 ) then
          do igaus = 1,pgaus
             gpadv(DEF_VECT,1:ndime,igaus) = gpadv(DEF_VECT,1:ndime,igaus) + gpsgs(DEF_VECT,1:ndime,igaus,1)
          end do
       end if
       !
       ! Temperature and wmean
       !
       if( (kfl_coupl(ID_NASTIN,ID_TEMPER) /= 0 .or. kfl_cotem_nsi == 1 .or. kfl_regim_nsi==3 ) ) then
          if (kfl_surte_nsi == 2_ip) then
             do ivect_loc = 1,VECTOR_SIZE
                ielem_loc = abs(list_elements(ivect_loc))
                if (ielem_loc > 0) then
                   do itime = 1,nbdfp_nsi
                      call ker_proper('DENSI','PGAUS',dummi,ielem_loc,densi(ivect_loc,1:pgaus,itime),pnode,pgaus,gpsha(ivect_loc,1:mnode,1:pgaus),gpcar(ivect_loc,1:ndime,1:mnode,1:pgaus))     
                   end do
                end if
             end do
          else

             if (kfl_lookg_nsi > 0) then

                !
                ! AB: Vectorize this!
                !

                if (kfl_coupl(ID_NASTIN,ID_CHEMIC) /= 0 .and. kfl_coupl(ID_NASTIN,ID_TEMPER) == 0) then  !!! FOR CMC
                   do ivect_loc = 1,VECTOR_SIZE
                      ielem_loc = abs(list_elements(ivect_loc))
                      if (ielem_loc > 0) then
                         call ker_proper('DENSI','PGAUS',dummi,ielem_loc,gpden2,pnode,pgaus,gpsha(ivect_loc,1:pnode,1:pgaus),gpcar(ivect_loc,1:ndime,1:pnode,1:pgaus))
                         do igaus=1,pgaus
                            do itime = 1,nbdfp_nsi
                               !!!densi(ivect_loc,igaus,itime) = prthe(itime)*densi_gp(ielem_loc) % a(igaus,1,1) / prthe(1)
                               densi(ivect_loc,igaus,itime) = prthe(itime)*gpden2(igaus) / prthe(1)
                            enddo
                         enddo
                      endif
                   enddo

                else
                   do ivect_loc = 1,VECTOR_SIZE
                      ielem_loc = abs(list_elements(ivect_loc))
                      if (ielem_loc > 0) then
                         do igaus=1,pgaus
                            do itime = 1,nbdfp_nsi
                               densi(ivect_loc,igaus,itime) = prthe(itime)*wmean_gp(ielem_loc) % a(igaus,1,1)/(gasco*tempe_gp(ielem_loc) % a(igaus,1,1))
                            enddo
                         enddo
                      endif
                   enddo
                end if

             else
                gptem    = 0.0_rp
                gpwmean  = 0.0_rp

                do igaus = 1,pgaus
                   do itime = 1,nbdfp_nsi
                      do inode = 1,pnode
                         gptem(DEF_VECT,igaus,itime)   = gptem(DEF_VECT,igaus,itime)   + eltem(DEF_VECT,inode,itime)   * gpsha(DEF_VECT,inode,igaus)
                         gpwmean(DEF_VECT,igaus,itime) = gpwmean(DEF_VECT,igaus,itime) + elwmean(DEF_VECT,inode,itime) * gpsha(DEF_VECT,inode,igaus)
                      end do
                   end do
                end do

                if (kfl_regim_nsi == 3) then   ! loads density in last time steps
                   if(associated(tesgs)) then
                      do igaus = 1,pgaus
                         do itime =1,nbdfp_nsi
                            densi(DEF_VECT,igaus, itime) =  prthe(itime) * gpwmean(DEF_VECT,igaus,itime) &
                               /gasco / ( gptem(DEF_VECT,igaus,itime) + gpsgt(DEF_VECT,igaus,itime))
                         end do
                      end do
                   else
                      do igaus = 1,pgaus
                         do itime =1,nbdfp_nsi
                            densi(DEF_VECT,igaus, itime) =  prthe(itime) * gpwmean(DEF_VECT,igaus,itime) / ( gasco * ( gptem(DEF_VECT,igaus,itime) + 1.0e-6_rp) )
                         end do
                      end do
                   end if
                end if
             endif
          endif


       end if
       !
       ! Surface tension
       !
       if( kfl_surte_nsi /= 0 ) then
          do igaus = 1,pgaus
             gpcur(DEF_VECT,igaus) = 0.0_rp
             do idime = 1,ndime
                gpnor(DEF_VECT,idime,igaus) = 0.0_rp
             end do
             do inode = 1,pnode
                gpcur(DEF_VECT,igaus) = gpcur(DEF_VECT,igaus) + elcur(DEF_VECT,inode) * gpsha(DEF_VECT,inode,igaus)
                do idime = 1,ndime
                   gpnor(DEF_VECT,idime,igaus) = gpnor(DEF_VECT,idime,igaus) + elnor(DEF_VECT,idime,inode) * gpsha(DEF_VECT,inode,igaus)
                end do
             end do
          end do
       end if
       !
       ! Lagrange multiplier -(eps(u):eps(v))
       ! GPLAG(I,J) = 1/2 * ( dw_j / dx_i + dw_i / dx_j )
       !
       if( kfl_immer_nsi == 1 ) then

          gplag(DEF_VECT,:,:,:) = 0.0_rp
          do igaus = 1,pgaus
             do jdime = 1,ndime
                do idime = 1,ndime
                   do inode = 1,pnode
                      gplag(DEF_VECT,idime,jdime,igaus) = gplag(DEF_VECT,idime,jdime,igaus) &
                           + 0.5_rp * ( ellag(DEF_VECT,jdime,1,inode) * gpcar(DEF_VECT,idime,inode,igaus) &
                           &          + ellag(DEF_VECT,idime,1,inode) * gpcar(DEF_VECT,jdime,inode,igaus) )
                   end do
                end do
             end do
          end do

       else if( kfl_immer_nsi == 2 ) then

          gplag(DEF_VECT,:,:,:) = 0.0_rp
          do igaus = 1,pgaus
             do jdime = 1,ndime
                do idime = 1,ndime
                   do inode = 1,pnode
                      gplag(DEF_VECT,idime,jdime,igaus) = gplag(DEF_VECT,idime,jdime,igaus) &
                           + ellag(DEF_VECT,idime,jdime,inode) * gpsha(DEF_VECT,inode,igaus)
                   end do
                end do
             end do
          end do
       end if

       !----------------------------------------------------------------------
       !
       ! Projections at Gauss points
       !
       !----------------------------------------------------------------------

       if( kfl_stabi_nsi /= NSI_ASGS .and. kfl_stabi_nsi /= NSI_GALERKIN .and. kfl_stabi_nsi /= NSI_ALGEBRAIC_SPLIT_OSS ) then

          do igaus = 1,pgaus
             do inode = 1,pnode
                do idime = 1,ndime
                   gpvep(DEF_VECT,idime,igaus) = gpvep(DEF_VECT,idime,igaus) + gpsha(DEF_VECT,inode,igaus) * elvep(DEF_VECT,idime,inode)
                end do
                gpprp(DEF_VECT,igaus) = gpprp(DEF_VECT,igaus) + gpsha(DEF_VECT,inode,igaus) * elprp(DEF_VECT,inode)
             end do
          end do

          if( kfl_stabi_nsi == NSI_SPLIT_OSS ) then
             gpgrp = 0.0_rp
             do igaus = 1,pgaus
                do inode = 1,pnode
                   do idime = 1,ndime
                      gpgrp(DEF_VECT,idime,igaus) = gpgrp(DEF_VECT,idime,igaus) + gpsha(DEF_VECT,inode,igaus) * elgrp(DEF_VECT,idime,inode)
                   end do
                end do
             end do
          end if

       end if

       if( kfl_stabi_nsi == NSI_GALERKIN .or. kfl_stabi_nsi == NSI_ALGEBRAIC_SPLIT_OSS ) then

          !----------------------------------------------------------------------
          !
          ! No stabilization
          !
          !----------------------------------------------------------------------

          continue

       else if( kfl_stabi_nsi == NSI_SPLIT_OSS ) then

          !----------------------------------------------------------------------
          !
          ! Split OSS: only initialize RMOMU
          !
          !----------------------------------------------------------------------

          rmomu(DEF_VECT,:,:) = 0.0_rp

       else

          !----------------------------------------------------------------------
          !
          ! ASGS and full OSS: residual Gauss point values
          !
          !----------------------------------------------------------------------

          if( plapl == 1 .and. ndime == 2 ) then
             do igaus = 1,pgaus
                do inode = 1,pnode
                   gplap(DEF_VECT,inode,igaus) =        &
                        + gphes(DEF_VECT,1,inode,igaus) &
                        + gphes(DEF_VECT,2,inode,igaus)
                end do
             end do
          else if( plapl == 1 .and. ndime == 3 ) then
             do igaus = 1,pgaus
                do inode = 1,pnode
                   gplap(DEF_VECT,inode,igaus) =        &
                        + gphes(DEF_VECT,1,inode,igaus) &
                        + gphes(DEF_VECT,2,inode,igaus) &
                        + gphes(DEF_VECT,3,inode,igaus)
                end do
             end do
          else
             do igaus = 1,pgaus
                do inode = 1,pnode
                   gplap(DEF_VECT,inode,igaus) = 0.0_rp
                end do
             end do
          end if

          !----------------------------------------------------------------------
          !
          ! RMOMU, RCONT: Momentum and continuity residuals
          !
          !----------------------------------------------------------------------
          !
          ! rho/(dt*theta)*u + sig*u - mu*d^2u/dxk^2
          !
          if (NSI_FRACTIONAL_STEP) then
             dtinv_mod = 0.0_rp 
          else
             dtinv_mod = dtinv_res
          end if
          do igaus = 1,pgaus
             fact1(DEF_VECT) = dtinv_mod(DEF_VECT) * pabdf_nsi(1) * gpden(DEF_VECT,igaus) + gppor(DEF_VECT,igaus)
             do inode = 1,pnode
                rmomu(DEF_VECT,inode,igaus) = fact1(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) - gpvis(DEF_VECT,igaus) * gplap(DEF_VECT,inode,igaus)
             end do
          end do
          !
          ! Continuity and Momentum turbulent viscosity: -grad(mu).grad(u)
          !
          do igaus = 1,pgaus
             do inode = 1,pnode
                rcont(DEF_VECT,1:ndime,inode,igaus) = gpcar(DEF_VECT,1:ndime,inode,igaus)
                do idime = 1,ndime
                   rmomu(DEF_VECT,inode,igaus) = rmomu(DEF_VECT,inode,igaus) - gpgvi(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
                end do
             end do
          end do

       end if

       !----------------------------------------------------------------------
       !
       ! GPRHS
       !
       ! 1. ASGS + Full OSS:
       !    RHS = rho*dt*u^n + rho*g - rho*beta*g*(T-Tr)
       !    RHS_SGS = rho/dt*u'n + proje
       ! 2. Split OSS: only force term without temporal one
       !    RHS = rho*g - rho*beta*g*(T-Tr)
       !
       !----------------------------------------------------------------------
       !
       ! GPHYD: Hydrostatic density
       !
       if( kfl_regim_nsi == 3 .and. kfl_prthe_nsi == 0 ) then
          !
          ! Low-Mach regime with open flow, constant thermodynamic pressure
          !
          if (kfl_coupl(ID_TEMPER,ID_CHEMIC) /= 0 .or. kfl_coupl(ID_NASTIN,ID_CHEMIC) /= 0 ) then
             w1 = 0.0289_rp ! Molecular weight of air for reference
          else
             w1 = 1.0_rp
          endif
          ! gphyd(1:pgaus) =  prthe(1)*w1/(lowtr_nsi * gasco)
          gphyd = 0.0_rp

       else if( kfl_hydro_gravity_nsi == 0 ) then

          gphyd = 0.0_rp

       end if
       !
       ! GPGVE: grad(u)
       !
       if (kfl_regim_nsi==3 .and. kfl_convection_type_nsi == NSI_CONVECTION_EMAC ) then
          do igaus = 1,pgaus
             do idime = 1,ndime
                do jdime = 1,ndime
                   gpgve(DEF_VECT,jdime,idime,igaus) = 0.0_rp
                end do
             end do
             do inode = 1,pnode
                do idime = 1,ndime
                   do jdime = 1,ndime
                      gpgve(DEF_VECT,jdime,idime,igaus) = gpgve(DEF_VECT,jdime,idime,igaus) &
                         + gpcar(DEF_VECT,jdime,inode,igaus) * elvel(DEF_VECT,idime,inode,1) * &
                         & elden(DEF_VECT,inode) / sqrt(elden(DEF_VECT,inode))
                   end do
                end do
             end do
          end do
       elseif (kfl_regim_nsi==3 .and. kfl_convection_type_nsi == NSI_CONVECTION_SKEW ) then
          do igaus = 1,pgaus
             do idime = 1,ndime
                do jdime = 1,ndime
                   gpgve(DEF_VECT,jdime,idime,igaus) = 0.0_rp
                end do
             end do
             do inode = 1,pnode
                do idime = 1,ndime
                   do jdime = 1,ndime
                      gpgve(DEF_VECT,jdime,idime,igaus) = gpgve(DEF_VECT,jdime,idime,igaus) &
                         + gpcar(DEF_VECT,jdime,inode,igaus) * elvel(DEF_VECT,idime,inode,1) * &
                         & elden(DEF_VECT,inode) 
                   end do
                end do
             end do
          end do
       else
          do igaus = 1,pgaus
             do idime = 1,ndime
                do jdime = 1,ndime
                   gpgve(DEF_VECT,jdime,idime,igaus) = 0.0_rp
                end do
             end do
             do inode = 1,pnode
                do idime = 1,ndime
                   do jdime = 1,ndime
                      gpgve(DEF_VECT,jdime,idime,igaus) = gpgve(DEF_VECT,jdime,idime,igaus) &
                         + gpcar(DEF_VECT,jdime,inode,igaus) * elvel(DEF_VECT,idime,inode,1)
                   end do
                end do
             end do
          end do
       end if
       !
       ! Time derivative: rho/dt*u^n
       !
       if( kfl_stabi_nsi /= NSI_SPLIT_OSS .and. kfl_stabi_nsi /= NSI_GALERKIN .and. kfl_stabi_nsi /= NSI_ALGEBRAIC_SPLIT_OSS) then
          do itime = 2, nbdfp_nsi
             do igaus = 1,pgaus
                fact1(DEF_VECT) = gpden(DEF_VECT,igaus) * dtinv_mod(DEF_VECT) * pabdf_nsi(itime)
                do idime = 1,ndime
                   gprhs(DEF_VECT,idime,igaus) = gprhs(DEF_VECT,idime,igaus) - fact1(DEF_VECT) * gpvel(DEF_VECT,idime,igaus,itime)
                end do
             end do
          end do
       end if
       !
       ! Gravity: rho*g
       !
       do igaus = 1,pgaus
          fact2(DEF_VECT) = ( gpden(DEF_VECT,igaus) - gphyd(DEF_VECT,igaus) ) * grnor_nsi
          do idime = 1,ndime
             gprhs(DEF_VECT,idime,igaus) = gprhs(DEF_VECT,idime,igaus) + fact2(DEF_VECT) * gravi_nsi(idime)
          end do
       end do
       !
       ! Boussinesq coupling: -rho*beta*g*(T-Tr)
       !
       if( kfl_cotem_nsi == 1) then
          fact1 = bougr_nsi * boube_nsi
          do igaus = 1,pgaus
#ifdef GABLS1TREF
             !
             ! Here I could just define the y coordinate
             !
             gpcod = 0.0_rp
             do idime = 1,ndime
                do inode = 1,pnode
                   gpcod(DEF_VECT,idime) = gpcod(DEF_VECT,idime) + gpsha(DEF_VECT,inode,igaus) * elcod(DEF_VECT,idime,inode)
                end do
             end do

             twall     = 265.0_rp - 6.9444e-05 * cutim
             ttop      = 268.0_rp
             boutr_aux = twall + (ttop-twall) * gpcod(DEF_VECT,2) / 400.0_rp
             fact2 = gpden(DEF_VECT,igaus) * fact1(DEF_VECT) * ( gptem(DEF_VECT,igaus,1) - boutr_aux )
#else
             fact2 = gpden(DEF_VECT,igaus) * fact1(DEF_VECT) * ( gptem(DEF_VECT,igaus,1) - boutr_nsi )
#endif
             do idime = 1,ndime
                gprhs(DEF_VECT,idime,igaus) = gprhs(DEF_VECT,idime,igaus) - fact2(DEF_VECT) * gravb_nsi(idime)
             end do
          end do
       end if
       !
       ! Multiphase flow problem
       !
       if( kfl_surte_nsi == 1 ) then
          !
          ! Classical FEM level set momentum conservation
          !
          call runend('VECTOR THIS FOR KFL_SURTE_NSI')
          !   do igaus = 1,pgaus
          !      do ivect = 1,VECTOR_SIZE
          !         if (abs(gpfle(ivect,igaus)) < thicl) then
          !            fact1(ivect)               = (0.5_rp/thicl) * (1.0_rp + cos ( pi * gpfle(ivect,igaus) / thicl ) )
          !            fact2(ivect)               = surte_nsi * gpcur(ivect,igaus) * fact1(ivect)
          !            do idime = 1,ndime
          !               gprhs(ivect,idime,igaus) = gprhs(ivect,idime,igaus) - fact2(ivect) * gpnor(ivect,idime,igaus)
          !               !gprhs(ivect,2,igaus) = gprhs(ivect,2,igaus) - fact2(ivect) * gpnor(2,igaus)
          !            end do
          !         end if
          !      end do
          !   end do

       else if ( kfl_surte_nsi == 2 ) then
          !
          ! ELSA model momentum conservation
          !
          do igaus = 1,pgaus
             fact1(DEF_VECT) = surte_nsi * gpcur(DEF_VECT,igaus)/gpden(DEF_VECT,igaus)
             do idime = 1,ndime
                gprhs(DEF_VECT,idime,igaus) = gprhs(DEF_VECT,idime,igaus) - fact1(DEF_VECT) * gpnor(DEF_VECT,idime,igaus)
             end do
          end do

       end if
       !
       ! External force
       !
       if( kfl_force_nsi < 0 ) then
          do igaus = 1,pgaus
             gpcod = 0.0_rp
             do idime = 1,ndime
                do inode = 1,pnode
                   gpcod(DEF_VECT,idime) = gpcod(DEF_VECT,idime) + gpsha(DEF_VECT,inode,igaus) * elcod(DEF_VECT,idime,inode)
                end do
             end do
             do ivect = 1,VECTOR_SIZE
                call ker_space_time_function(&
                     -kfl_force_nsi,gpcod(ivect,1),gpcod(ivect,2),gpcod(ivect,ndime),cutim,gpext(ivect,:))
             end do
             gprhs(DEF_VECT,1:ndime,igaus) = gprhs(DEF_VECT,1:ndime,igaus) + gpext(DEF_VECT,1:ndime)
          end do
       end if
       !
       ! Subgrid scale time derivative: rho/dt*u'^n
       !
       if(  kfl_sgsti_nsi == 1                       .and. &
            kfl_stabi_nsi /= NSI_SPLIT_OSS           .and. &
            kfl_stabi_nsi /= NSI_GALERKIN            .and. &
            kfl_stabi_nsi /= NSI_ALGEBRAIC_SPLIT_OSS       ) then
          if( kfl_regim_nsi == 3 ) then
             do itime = 2,2
                do igaus = 1,pgaus
                   do idime = 1,ndime
                      gprhs_sgs(DEF_VECT,idime,igaus) = &
                           gprhs_sgs(DEF_VECT,idime,igaus) + densi(DEF_VECT,igaus,itime) * gpsgs(DEF_VECT,idime,igaus,itime) * dtsgs_nsi
                   end do
                end do
             end do
          else
             do igaus = 1,pgaus
                fact1(DEF_VECT) = gpden(DEF_VECT,igaus) * dtsgs_nsi
                do idime = 1,ndime
                   gprhs_sgs(DEF_VECT,idime,igaus) = gprhs_sgs(DEF_VECT,idime,igaus) &
                        + fact1(DEF_VECT) * gpsgs(DEF_VECT,idime,igaus,2) ! + gpvep(1,igaus) / gpst1(igaus)
                end do
             end do
          end if
       end if

       !
       !!FER Note: Variable viscosity is now assembled on the left side
       !

!!$    else
!!$       !
!!$       ! Time derivative: rho/dt*u^n
!!$       !
!!$       if( kfl_stabi_nsi /= 2 ) then
!!$          do itime = 2, nbdfp_nsi
!!$             do igaus = 1,pgaus
!!$                fact1 = gpden(igaus) * dtinv_res * pabdf_nsi(itime)
!!$                gprhs(1,igaus) = gprhs(1,igaus) - fact1 * gpvel(1,igaus,itime)
!!$                gprhs(2,igaus) = gprhs(2,igaus) - fact1 * gpvel(2,igaus,itime)
!!$                gprhs(3,igaus) = gprhs(3,igaus) - fact1 * gpvel(3,igaus,itime)
!!$             end do
!!$          end do
!!$       end if
!!$       !
!!$       ! Gravity: rho*g
!!$       !
!!$       do igaus = 1,pgaus
!!$          fact2 = ( gpden(igaus) - gphyd(igaus) ) * grnor_nsi
!!$          gprhs(1,igaus) = gprhs(1,igaus) + fact2 * gravi_nsi(1)
!!$          gprhs(2,igaus) = gprhs(2,igaus) + fact2 * gravi_nsi(2)
!!$          gprhs(3,igaus) = gprhs(3,igaus) + fact2 * gravi_nsi(3)
!!$       end do
!!$       !
!!$       ! Boussinesq coupling: -rho*beta*g*(T-Tr)
!!$       !
!!$       if( kfl_cotem_nsi == 1) then
!!$          fact1 = bougr_nsi * boube_nsi
!!$          do igaus = 1,pgaus
!!$             fact2 = gpden(igaus) * fact1 * ( gptem(igaus,1) - boutr_nsi )
!!$             gprhs(1,igaus) = gprhs(1,igaus) - fact2 * gravb_nsi(1)
!!$             gprhs(2,igaus) = gprhs(2,igaus) - fact2 * gravb_nsi(2)
!!$             gprhs(3,igaus) = gprhs(3,igaus) - fact2 * gravb_nsi(3)
!!$          end do
!!$       end if
!!$       !
!!$       ! Surface tension: sigma*curvature*regularized_delta*normal
!!$       !
!!$       if( kfl_surte_nsi == 1) then
!!$          do igaus = 1,pgaus
!!$             if ( abs(gpfle(igaus)) < thicl ) then
!!$                fact1 = (0.5_rp/thicl) * (1.0_rp + cos ( pi * gpfle(igaus) / thicl ) )
!!$                fact2 = surte_nsi * gpcur(igaus) * fact1
!!$                gprhs(1,igaus) = gprhs(1,igaus) - fact2 * gpnor(1,igaus)
!!$                gprhs(2,igaus) = gprhs(2,igaus) - fact2 * gpnor(2,igaus)
!!$                gprhs(3,igaus) = gprhs(3,igaus) - fact2 * gpnor(3,igaus)
!!$             end if
!!$          end do
!!$       end if
!!$       !
!!$       ! Subgrid scale time derivative: rho/dt*u'n
!!$       !
!!$       if( kfl_sgsti_nsi == 1 .and. kfl_stabi_nsi /= 2 ) then
!!$          if (kfl_regim_nsi==3) then
!!$             do itime =2, 2
!!$                do igaus = 1,pgaus
!!$                   gprhs_sgs(1,igaus) = gprhs_sgs(1,igaus) + densi(igaus, itime)*gpsgs(1,igaus,itime)*dtsgs_nsi
!!$                   gprhs_sgs(2,igaus) = gprhs_sgs(2,igaus) + densi(igaus, itime)*gpsgs(2,igaus,itime)*dtsgs_nsi
!!$                   gprhs_sgs(3,igaus) = gprhs_sgs(3,igaus) + densi(igaus, itime)*gpsgs(3,igaus,itime)*dtsgs_nsi
!!$                end do
!!$             end do
!!$          else
!!$             do igaus = 1,pgaus
!!$                fact1 = gpden(igaus) * dtsgs_nsi
!!$                gprhs_sgs(1,igaus) = gprhs_sgs(1,igaus) + fact1 * gpsgs(1,igaus,2) !+ gpvep(1,igaus)/ gpst1(igaus)
!!$                gprhs_sgs(2,igaus) = gprhs_sgs(2,igaus) + fact1 * gpsgs(2,igaus,2) !+ gpvep(2,igaus)/ gpst1(igaus)
!!$                gprhs_sgs(3,igaus) = gprhs_sgs(3,igaus) + fact1 * gpsgs(3,igaus,2) !+ gpvep(3,igaus)/ gpst1(igaus)
!!$             end do
!!$          end if
!!$       end if
!!$    end if

       !----------------------------------------------------------------------
       !
       ! RMOM2
       !
       ! Off-diagonal momentum operator
       !
       !----------------------------------------------------------------------

       if( kfl_stabi_nsi /= NSI_GALERKIN  ) then

          if( kfl_rmom2_nsi /= 0 ) rmom2 = 0.0_rp
          !
          ! Coriolis force
          !
          ! 2*rho*(w x u)
          ! x-equation: w x u = wy*uz - wz*uy
          ! y-equation: w x u = wz*ux - wx*uz
          ! z-equation: w x u = wx*uy - wy*ux
          !
          !
          if( corio_nsi > 1.0e-12_rp ) then
             if( ndime == 2 ) then
                do igaus = 1,pgaus
                   fact3(DEF_VECT) = 2.0_rp * gpden(DEF_VECT,igaus) * fvela_nsi(3)
                   do inode = 1,pnode
                      rmom2(DEF_VECT,1,2,inode,igaus) = -fact3(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)     ! -wz*uy
                      rmom2(DEF_VECT,2,1,inode,igaus) =  fact3(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)     !  wz*ux
                   end do
                end do
             else
                do igaus = 1,pgaus
                   fact0(DEF_VECT) = 2.0_rp * gpden(DEF_VECT,igaus)
                   fact1(DEF_VECT) = fact0(DEF_VECT) * fvela_nsi(1)
                   fact2(DEF_VECT) = fact0(DEF_VECT) * fvela_nsi(2)
                   fact3(DEF_VECT) = fact0(DEF_VECT) * fvela_nsi(3)
                   do inode = 1,pnode
                      rmom2(DEF_VECT,1,2,inode,igaus) = - fact3(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! -wz*uy
                      rmom2(DEF_VECT,1,3,inode,igaus) =   fact2(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  !  wy*uz
                      rmom2(DEF_VECT,2,1,inode,igaus) =   fact3(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  !  wz*ux
                      rmom2(DEF_VECT,2,3,inode,igaus) = - fact1(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! -wx*uz
                      rmom2(DEF_VECT,3,1,inode,igaus) = - fact2(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! -wy*ux
                      rmom2(DEF_VECT,3,2,inode,igaus) =   fact1(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  !  wx*uy
                   end do
                end do
             end if
          end if
          !
          ! Anisotropic porosity
          !
          if( kfl_anipo_nsi /= 0 ) then
             do igaus = 1,pgaus
                do inode = 1,pnode
                   do idime = 1,ndime
                      do jdime = 1,ndime
                         rmom2(DEF_VECT,idime,jdime,inode,igaus) = rmom2(DEF_VECT,idime,jdime,inode,igaus) &
                              + gpapo(DEF_VECT,jdime,idime,igaus) * gpsha(DEF_VECT,inode,igaus)
                      end do
                   end do
                end do
             end do
          end if
          !
          !!FER CAREFUL THIS DOES NOT TAKE INTO ACCOUNT CHEMIC
          ! Viscous term: only terms from divergence and complete forms
          !
          ! - ( div[-2 mu eps(u)] , v ) = - d/dxj ( -2*mu* ( dui/dxj + duj/dxi - 2/3*mu div(u) delta_ij )
          !
          ! Laplacian  form: A=0, B=0, eps(u) = 1/2 grad(u)
          ! Divergence form: A=1, B=0, eps(u) = 1/2 ( grad(u) + grad(u)^t )
          ! Complete   form: A=1, B=1, eps(u) = 1/2 ( grad(u) + grad(u)^t ) - 1/3 div(u) I
          !
          ! + ( mu d^2ui/dxj^2 , vi )               (6)        already computed
          ! + ( dmu/dxj dui/dxj , vi )              (7)        already computed
          ! + A * ( mu  d^2uj/dxidxj , vi )         (8)        divergence
          ! + A * ( dmu/dxj duj/dxi , vi )          (9)        divergence
          ! - 2/3 * B * ( dmu/dxi (div u) , vi )   (10)        complete
          ! - 2/3 * B * ( mu d(div u)/dxi , vi )   (11)        complete
          !
          if( fvins_nsi > 0.9_rp ) then

             if( fvins_nsi > 1.9_rp ) then
                xvis2 = 2.0_rp / 3.0_rp
             else
                xvis2 = 0.0_rp
             end if

             if( ndime == 2 ) then
                do igaus = 1,pgaus
                   fact0(DEF_VECT) = xvis2 * gpvis(DEF_VECT,igaus) - gpvis(DEF_VECT,igaus)
                   fact1(DEF_VECT) = xvis2 * gpgvi(DEF_VECT,1,igaus)
                   fact2(DEF_VECT) = xvis2 * gpgvi(DEF_VECT,2,igaus)
                   do inode = 1,pnode
                      rmom2(DEF_VECT,1,1,inode,igaus) = rmom2(DEF_VECT,1,1,inode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,1,inode,igaus)            ! (8) + (11)  - mu * d^2ux/dx^2 + 2/3 * mu * d^2ux/dx^2
                      rmom2(DEF_VECT,1,2,inode,igaus) = rmom2(DEF_VECT,1,2,inode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,3,inode,igaus)            ! (8) + (11)  - mu * d^2uy/dxdy + 2/3 * mu * d^2uy/dxdy
                      rmom2(DEF_VECT,2,1,inode,igaus) = rmom2(DEF_VECT,2,1,inode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,3,inode,igaus)            ! (8) + (11)  - mu * d^2ux/dxdy + 2/3 * mu * d^2ux/dxdy
                      rmom2(DEF_VECT,2,2,inode,igaus) = rmom2(DEF_VECT,2,2,inode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,2,inode,igaus)            ! (8) + (11)  - mu * d^2uy/dy^2 + 2/3 * mu * d^2uy/dy^2

                      rmom2(DEF_VECT,1,1,inode,igaus) = rmom2(DEF_VECT,1,1,inode,igaus) - gpgvi(DEF_VECT,1,igaus) * gpcar(DEF_VECT,1,inode,igaus)   ! (9)         - dmu/dx * dux/dx
                      rmom2(DEF_VECT,1,2,inode,igaus) = rmom2(DEF_VECT,1,2,inode,igaus) - gpgvi(DEF_VECT,2,igaus) * gpcar(DEF_VECT,1,inode,igaus)   ! (9)         - dmu/dy * duy/dx
                      rmom2(DEF_VECT,2,1,inode,igaus) = rmom2(DEF_VECT,2,1,inode,igaus) - gpgvi(DEF_VECT,1,igaus) * gpcar(DEF_VECT,2,inode,igaus)   ! (9)         - dmu/dx * dux/dy
                      rmom2(DEF_VECT,2,2,inode,igaus) = rmom2(DEF_VECT,2,2,inode,igaus) - gpgvi(DEF_VECT,2,igaus) * gpcar(DEF_VECT,2,inode,igaus)   ! (9)         - dmu/dy * dvy/dy

                      rmom2(DEF_VECT,1,1,inode,igaus) = rmom2(DEF_VECT,1,1,inode,igaus) + fact1(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus)            ! (10)        + 2/3 * dmu/dx * dux/dx
                      rmom2(DEF_VECT,1,2,inode,igaus) = rmom2(DEF_VECT,1,2,inode,igaus) + fact1(DEF_VECT) * gpcar(DEF_VECT,2,inode,igaus)            ! (10)        + 2/3 * dmu/dx * duy/dy
                      rmom2(DEF_VECT,2,1,inode,igaus) = rmom2(DEF_VECT,2,1,inode,igaus) + fact2(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus)            ! (10)        + 2/3 * dmu/dy * dux/dx
                      rmom2(DEF_VECT,2,2,inode,igaus) = rmom2(DEF_VECT,2,2,inode,igaus) + fact2(DEF_VECT) * gpcar(DEF_VECT,2,inode,igaus)            ! (10)        + 2/3 * dmu/dy * duy/dy
                   end do
                end do

             else if( ndime == 3 ) then
                do igaus = 1,pgaus
                   fact0(DEF_VECT) = xvis2 * gpvis(DEF_VECT,igaus) - gpvis(DEF_VECT,igaus)
                   fact1(DEF_VECT) = xvis2 * gpgvi(DEF_VECT,1,igaus)
                   fact2(DEF_VECT) = xvis2 * gpgvi(DEF_VECT,2,igaus)
                   fact3(DEF_VECT) = xvis2 * gpgvi(DEF_VECT,3,igaus)
                   do inode = 1,pnode
                      rmom2(DEF_VECT,1,1,inode,igaus) = rmom2(DEF_VECT,1,1,inode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,1,inode,igaus)            ! (8) + (11)  - mu * d^2ux/dx^2 + 2/3 * mu * d^2ux/dx^2
                      rmom2(DEF_VECT,1,2,inode,igaus) = rmom2(DEF_VECT,1,2,inode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,4,inode,igaus)            ! (8) + (11)  - mu * d^2uy/dxdy + 2/3 * mu * d^2uy/dxdy
                      rmom2(DEF_VECT,1,3,inode,igaus) = rmom2(DEF_VECT,1,3,inode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,5,inode,igaus)            ! (8) + (11)  - mu * d^2uz/dxdz + 2/3 * mu * d^2uz/dxdz

                      rmom2(DEF_VECT,2,1,inode,igaus) = rmom2(DEF_VECT,2,1,inode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,4,inode,igaus)            ! (8) + (11)  - mu * d^2ux/dxdy + 2/3 * mu * d^2ux/dxdy
                      rmom2(DEF_VECT,2,2,inode,igaus) = rmom2(DEF_VECT,2,2,inode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,2,inode,igaus)            ! (8) + (11)  - mu * d^2uy/dy^2 + 2/3 * mu * d^2uy/dy^2
                      rmom2(DEF_VECT,2,3,inode,igaus) = rmom2(DEF_VECT,2,3,inode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,6,inode,igaus)            ! (8) + (11)  - mu * d^2uz/dydz + 2/3 * mu * d^2uz/dydz

                      rmom2(DEF_VECT,3,1,inode,igaus) = rmom2(DEF_VECT,3,1,inode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,5,inode,igaus)            ! (8) + (11)  - mu * d^2ux/dxdz + 2/3 * mu * d^2ux/dxdz
                      rmom2(DEF_VECT,3,2,inode,igaus) = rmom2(DEF_VECT,3,2,inode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,6,inode,igaus)            ! (8) + (11)  - mu * d^2uy/dydz + 2/3 * mu * d^2uy/dydz
                      rmom2(DEF_VECT,3,3,inode,igaus) = rmom2(DEF_VECT,3,3,inode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,3,inode,igaus)            ! (8) + (11)  - mu * d^2uz/dz^2 + 2/3 * mu * d^2uz/dz^2

                      rmom2(DEF_VECT,1,1,inode,igaus) = rmom2(DEF_VECT,1,1,inode,igaus) - gpgvi(DEF_VECT,1,igaus) * gpcar(DEF_VECT,1,inode,igaus)   ! (9)         - dmu/dx * dux/dx
                      rmom2(DEF_VECT,1,2,inode,igaus) = rmom2(DEF_VECT,1,2,inode,igaus) - gpgvi(DEF_VECT,2,igaus) * gpcar(DEF_VECT,1,inode,igaus)   ! (9)         - dmu/dy * duy/dx
                      rmom2(DEF_VECT,1,3,inode,igaus) = rmom2(DEF_VECT,1,3,inode,igaus) - gpgvi(DEF_VECT,3,igaus) * gpcar(DEF_VECT,1,inode,igaus)   ! (9)         - dmu/dz * duz/dx

                      rmom2(DEF_VECT,2,1,inode,igaus) = rmom2(DEF_VECT,2,1,inode,igaus) - gpgvi(DEF_VECT,1,igaus) * gpcar(DEF_VECT,2,inode,igaus)   ! (9)         - dmu/dx * dux/dy
                      rmom2(DEF_VECT,2,2,inode,igaus) = rmom2(DEF_VECT,2,2,inode,igaus) - gpgvi(DEF_VECT,2,igaus) * gpcar(DEF_VECT,2,inode,igaus)   ! (9)         - dmu/dy * duy/dy
                      rmom2(DEF_VECT,2,3,inode,igaus) = rmom2(DEF_VECT,2,3,inode,igaus) - gpgvi(DEF_VECT,3,igaus) * gpcar(DEF_VECT,2,inode,igaus)   ! (9)         - dmu/dz * duz/dy

                      rmom2(DEF_VECT,3,1,inode,igaus) = rmom2(DEF_VECT,3,1,inode,igaus) - gpgvi(DEF_VECT,1,igaus) * gpcar(DEF_VECT,3,inode,igaus)   ! (9)         - dmu/dx * dux/dz
                      rmom2(DEF_VECT,3,2,inode,igaus) = rmom2(DEF_VECT,3,2,inode,igaus) - gpgvi(DEF_VECT,2,igaus) * gpcar(DEF_VECT,3,inode,igaus)   ! (9)         - dmu/dy * duy/dz
                      rmom2(DEF_VECT,3,3,inode,igaus) = rmom2(DEF_VECT,3,3,inode,igaus) - gpgvi(DEF_VECT,3,igaus) * gpcar(DEF_VECT,3,inode,igaus)   ! (9)         - dmu/dz * duz/dz

                      rmom2(DEF_VECT,1,1,inode,igaus) = rmom2(DEF_VECT,1,1,inode,igaus) + fact1(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus)            ! (10)        + 2/3 * dmu/dx * dux/dx
                      rmom2(DEF_VECT,1,2,inode,igaus) = rmom2(DEF_VECT,1,2,inode,igaus) + fact1(DEF_VECT) * gpcar(DEF_VECT,2,inode,igaus)            ! (10)        + 2/3 * dmu/dx * duy/dy
                      rmom2(DEF_VECT,1,3,inode,igaus) = rmom2(DEF_VECT,1,3,inode,igaus) + fact1(DEF_VECT) * gpcar(DEF_VECT,3,inode,igaus)            ! (10)        + 2/3 * dmu/dx * duz/dz

                      rmom2(DEF_VECT,2,1,inode,igaus) = rmom2(DEF_VECT,2,1,inode,igaus) + fact2(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus)            ! (10)        + 2/3 * dmu/dy * dux/dx
                      rmom2(DEF_VECT,2,2,inode,igaus) = rmom2(DEF_VECT,2,2,inode,igaus) + fact2(DEF_VECT) * gpcar(DEF_VECT,2,inode,igaus)            ! (10)        + 2/3 * dmu/dy * duy/dy
                      rmom2(DEF_VECT,2,3,inode,igaus) = rmom2(DEF_VECT,2,3,inode,igaus) + fact2(DEF_VECT) * gpcar(DEF_VECT,3,inode,igaus)            ! (10)        + 2/3 * dmu/dy * duz/dz

                      rmom2(DEF_VECT,3,1,inode,igaus) = rmom2(DEF_VECT,3,1,inode,igaus) + fact3(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus)            ! (10)        + 2/3 * dmu/dz * dux/dx
                      rmom2(DEF_VECT,3,2,inode,igaus) = rmom2(DEF_VECT,3,2,inode,igaus) + fact3(DEF_VECT) * gpcar(DEF_VECT,2,inode,igaus)            ! (10)        + 2/3 * dmu/dz * duy/dy
                      rmom2(DEF_VECT,3,3,inode,igaus) = rmom2(DEF_VECT,3,3,inode,igaus) + fact3(DEF_VECT) * gpcar(DEF_VECT,3,inode,igaus)            ! (10)        + 2/3 * dmu/dz * duz/dz

                   end do
                end do

             end if

          end if

          !
          ! Adjoint Newton term: rho*( uj d(adv)/d(xj) , v )
          !
          if( kfl_adj_prob == 1 ) then
             if( ndime == 2 ) then
                do igaus = 1,pgaus
                   do inode = 1,pnode
                      fact0(DEF_VECT) = gpden(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)
                      rmom2(DEF_VECT,1,1,inode,igaus) = rmom2(DEF_VECT,1,1,inode,igaus) + fact0(DEF_VECT) * gpgve(DEF_VECT,1,1,igaus)   ! rho * ux * d(adv_x)/dx
                      rmom2(DEF_VECT,1,2,inode,igaus) = rmom2(DEF_VECT,1,2,inode,igaus) + fact0(DEF_VECT) * gpgve(DEF_VECT,2,1,igaus)   ! rho * uy * d(adv_x)/dy
                      rmom2(DEF_VECT,2,1,inode,igaus) = rmom2(DEF_VECT,2,1,inode,igaus) + fact0(DEF_VECT) * gpgve(DEF_VECT,1,2,igaus)   ! rho * ux * d(adv_y)/dx
                      rmom2(DEF_VECT,2,2,inode,igaus) = rmom2(DEF_VECT,2,2,inode,igaus) + fact0(DEF_VECT) * gpgve(DEF_VECT,2,2,igaus)   ! rho * uy * d(adv_y)/dy
                   end do
                end do
             else
                do igaus = 1,pgaus
                   do inode = 1,pnode
                      fact0(DEF_VECT) = gpden(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)
                      rmom2(DEF_VECT,1,1,inode,igaus) = rmom2(DEF_VECT,1,1,inode,igaus) + fact0(DEF_VECT) * gpgve(DEF_VECT,1,1,igaus)   ! rho * ux * d(adv_x)/dx
                      rmom2(DEF_VECT,1,2,inode,igaus) = rmom2(DEF_VECT,1,2,inode,igaus) + fact0(DEF_VECT) * gpgve(DEF_VECT,2,1,igaus)   ! rho * uy * d(adv_x)/dy
                      rmom2(DEF_VECT,1,3,inode,igaus) = rmom2(DEF_VECT,1,3,inode,igaus) + fact0(DEF_VECT) * gpgve(DEF_VECT,3,1,igaus)   ! rho * uz * d(adv_x)/dz
                      rmom2(DEF_VECT,2,1,inode,igaus) = rmom2(DEF_VECT,2,1,inode,igaus) + fact0(DEF_VECT) * gpgve(DEF_VECT,1,2,igaus)   ! rho * ux * d(adv_y)/dx
                      rmom2(DEF_VECT,2,2,inode,igaus) = rmom2(DEF_VECT,2,2,inode,igaus) + fact0(DEF_VECT) * gpgve(DEF_VECT,2,2,igaus)   ! rho * uy * d(adv_y)/dy
                      rmom2(DEF_VECT,2,3,inode,igaus) = rmom2(DEF_VECT,2,3,inode,igaus) + fact0(DEF_VECT) * gpgve(DEF_VECT,3,2,igaus)   ! rho * uz * d(adv_y)/dz
                      rmom2(DEF_VECT,3,1,inode,igaus) = rmom2(DEF_VECT,3,1,inode,igaus) + fact0(DEF_VECT) * gpgve(DEF_VECT,1,3,igaus)   ! rho * ux * d(adv_z)/dx
                      rmom2(DEF_VECT,3,2,inode,igaus) = rmom2(DEF_VECT,3,2,inode,igaus) + fact0(DEF_VECT) * gpgve(DEF_VECT,2,3,igaus)   ! rho * uy * d(adv_z)/dy
                      rmom2(DEF_VECT,3,3,inode,igaus) = rmom2(DEF_VECT,3,3,inode,igaus) + fact0(DEF_VECT) * gpgve(DEF_VECT,3,3,igaus)   ! rho * uz * d(adv_z)/dz
                   end do
                end do
             end if
          end if

       end if ! not galerkin

       !----------------------------------------------------------------------
       !
       ! GPRHS
       !
       ! Rotation term:            f = f - rho * ( w x w x r + dw/dt x r )
       ! Linear acceleration term: f = f - rho * a
       !
       !----------------------------------------------------------------------

       if( corio_nsi > 1.0e-12_rp ) then

          do igaus = 1,pgaus
             !
             ! Gauss point coordinates wrt center of rotation
             !
             gpcod = 0.0_rp
             do idime = 1,ndime
                do inode = 1,pnode
                   gpcod(DEF_VECT,idime) = gpcod(DEF_VECT,idime) + gpsha(DEF_VECT,inode,igaus) * elcod(DEF_VECT,idime,inode)
                end do
                gpcod(DEF_VECT,idime) = gpcod(DEF_VECT,idime) - frotc_nsi(idime)
             end do
             !
             ! Angular acceleration dw/dt x r
             !
             if( ndime == 2 ) then
                alpha(DEF_VECT,1) =-facca_nsi(3) * gpcod(DEF_VECT,2)
                alpha(DEF_VECT,2) = facca_nsi(3) * gpcod(DEF_VECT,1)
             else if( ndime == 3 ) then
                alpha(DEF_VECT,1) = facca_nsi(2) * gpcod(DEF_VECT,3) - facca_nsi(3) * gpcod(DEF_VECT,2)
                alpha(DEF_VECT,2) = facca_nsi(3) * gpcod(DEF_VECT,1) - facca_nsi(1) * gpcod(DEF_VECT,3)
                alpha(DEF_VECT,3) = facca_nsi(1) * gpcod(DEF_VECT,2) - facca_nsi(2) * gpcod(DEF_VECT,1)
             end if
             !
             ! Centrifugal force w x (w x r)
             !
             if( frotc_nsi(1) < 1.0e10_rp ) then
                if( ndime == 2 ) then
                   dummr(DEF_VECT,1) =-fvela_nsi(3) * gpcod(DEF_VECT,2)
                   dummr(DEF_VECT,2) = fvela_nsi(3) * gpcod(DEF_VECT,1)
                   centf(DEF_VECT,1) =-fvela_nsi(3) * dummr(DEF_VECT,2)
                   centf(DEF_VECT,2) = fvela_nsi(3) * dummr(DEF_VECT,1)
                else if( ndime ==3 ) then
                   dummr(DEF_VECT,1) = fvela_nsi(2) * gpcod(DEF_VECT,3) - fvela_nsi(3) * gpcod(DEF_VECT,2)
                   dummr(DEF_VECT,2) = fvela_nsi(3) * gpcod(DEF_VECT,1) - fvela_nsi(1) * gpcod(DEF_VECT,3)
                   dummr(DEF_VECT,3) = fvela_nsi(1) * gpcod(DEF_VECT,2) - fvela_nsi(2) * gpcod(DEF_VECT,1)
                   centf(DEF_VECT,1) = fvela_nsi(2) * dummr(DEF_VECT,3) - fvela_nsi(3) * dummr(DEF_VECT,2)
                   centf(DEF_VECT,2) = fvela_nsi(3) * dummr(DEF_VECT,1) - fvela_nsi(1) * dummr(DEF_VECT,3)
                   centf(DEF_VECT,3) = fvela_nsi(1) * dummr(DEF_VECT,2) - fvela_nsi(2) * dummr(DEF_VECT,1)
                end if
                centf(DEF_VECT,1:3) = centf(DEF_VECT,1:3) * centr_nsi
             else
                centf = 0.0_rp
             end if
             !
             ! Total force: rho * [ - w x (w x r) - dw/dt x r - a ]
             !
             do idime = 1,ndime
                gprhs(DEF_VECT,idime,igaus) = gprhs(DEF_VECT,idime,igaus) - gpden(DEF_VECT,igaus) &
                     &               * (   centf(DEF_VECT,idime)                &     ! w x (w x r)
                     &                   + alpha(DEF_VECT,idime)                &     ! dw/dt x r
                     &                   + faccl_nsi(idime)                     )     ! a
             end do
          end do

       end if

       !----------------------------------------------------------------------
       !
       ! Low-Mach regime: We calculate RHS of continuity equation
       !
       !----------------------------------------------------------------------

       !
       ! Integrate by parts mass conservation equation
       !
       !if( kfl_regim_nsi == 3 .and. (.not.NSI_FRACTIONAL_STEP)) then
       if( kfl_regim_nsi == 3 ) then 
          !
          ! rho = p^th / RT
          !
#ifdef matiaslma
          ! integrate by pats mass conservation equation
          !  gprhc2 = (vel grad rho + dtrho )/rho
          ! as the divergence term is assembled integrating by parts, we need to
          ! only take account for time derivative of density
          do igaus = 1,pgaus
             ! add time derivative of density
             fact0(DEF_VECT) = gpden(DEF_VECT,igaus) * pabdf_nsi(1)
             do itime = 2, nbdfp_nsi
                fact0(DEF_VECT) = fact0(DEF_VECT) + pabdf_nsi(itime) * densi(DEF_VECT,igaus,itime)
             end do
             fact0(DEF_VECT) = fact0(DEF_VECT) * dtinv_res(DEF_VECT)

             gprhc(DEF_VECT,igaus) = - fact0(DEF_VECT) / gpden(DEF_VECT,igaus)

             do idime = 1,ndime
                gpgte(DEF_VECT,idime, igaus) =0.0_rp
                do inode = 1, pnode
                   gpgte(DEF_VECT,idime,igaus)  = gpgte(DEF_VECT,idime,igaus)   + eltem(inode,1) * gpcar(DEF_VECT,idime,inode,igaus)
                end do
                gprh2(DEF_VECT,igaus) = gprh2(DEF_VECT,igaus) +  gpadv(DEF_VECT,idime,igaus) * gpgte(DEF_VECT,idime,igaus)
             end do
             gprh2(DEF_VECT,igaus) = gprh2(DEF_VECT,igaus) / gptem(DEF_VECT,igaus,1) - fact0(DEF_VECT) / gpden(DEF_VECT,igaus)

          end do

#else
          if(.not. NSI_FRACTIONAL_STEP) then

             do igaus = 1,pgaus

                gpdet(DEF_VECT,igaus) = 0.0_rp ! time derivative of density
                gpgrt(DEF_VECT,igaus) = 0.0_rp ! vel grad rho
                fact0(DEF_VECT) = gpden(DEF_VECT,igaus) * pabdf_nsi(1)
                do itime = 2,nbdfp_nsi
                   fact0(DEF_VECT) = fact0(DEF_VECT) + pabdf_nsi(itime) * densi(DEF_VECT,igaus,itime)
                end do
                !
                ! gpdet = Delta rho / Delta t = rho^n+1 - rho^n / dt
                !
                gpdet(DEF_VECT,igaus) = fact0(DEF_VECT) * dtinv_res(DEF_VECT)
                do idime = 1,ndime
                   gpgrt(DEF_VECT,igaus)  = gpgrt(DEF_VECT,igaus) + gpadv(DEF_VECT,idime,igaus) * gpgde(DEF_VECT,idime,igaus)
                end do
                gprhc(DEF_VECT,igaus) = ( - gpdet(DEF_VECT,igaus) - gpgrt(DEF_VECT,igaus) ) / gpden(DEF_VECT,igaus)
                gprh2(DEF_VECT,igaus) = gprhc(DEF_VECT,igaus)

             end do

          end if
#endif
       endif

#ifdef OPENACC
    end do
#endif

  end subroutine nsi_element_residual

  subroutine nsi_element_subgrid_scale(&
       pgaus,pnode,ndofn,ielem,chale,elvel,gpadv,gpvis,&
       gpden,rmom1,rmom2,gprhs,gpgpr,gpvel,gpcar,gpsp1,&
       gpsgs,gpsgi,gpgve,gpvep,gpgrp,gpst1,gprhs_sgs,  &
       dtsgs,resis_nsi,itsta_nsi,rmsgs_nsi,resgs_nsi,  &
       gppor)
    !------------------------------------------------------------------------
    !****f* Nastin/nsi_updsgs
    ! NAME
    !    nsi_updsgs
    ! DESCRIPTION
    !    This subroutine updates the subgrid scale
    ! USES
    ! USED BY
    !***
    !------------------------------------------------------------------------
    use def_kintyp, only       :  ip,rp
    use def_domain, only       :  ndime,mnode
    use def_nastin, only       :  staco_nsi,kfl_sgsti_nsi,&
         &                        kfl_sgsco_nsi,misgs_nsi,tosgs_nsi,&
         &                        kfl_taust_nsi,kfl_stabi_nsi,corio_nsi,&
         &                        kfl_sgsli_nsi,kfl_rmom2_nsi,relsg_nsi,&
         &                        kfl_stabi_nsi,mmsgs_nsi
    use mod_tauadr, only       :  tauadr
    implicit none
    integer(ip), intent(in)    :: pgaus,pnode,ndofn,ielem
    real(rp),    intent(in)    :: chale(VECTOR_SIZE,2)
    real(rp),    intent(in)    :: elvel(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(in)    :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)
    real(rp),    intent(in)    :: gprhs(VECTOR_SIZE,ndofn,pgaus)
    real(rp),    intent(in)    :: gprhs_sgs(VECTOR_SIZE,ndofn,pgaus)
    real(rp),    intent(in)    :: gpgpr(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: gpvel(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(inout) :: gpadv(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: gpvis(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpden(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: rmom1(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(in)    :: rmom2(VECTOR_SIZE,ndime,ndime,pnode,pgaus)
    real(rp),    intent(out)   :: gpsp1(VECTOR_SIZE,pgaus)
    real(rp),    intent(inout) :: gpsgs(VECTOR_SIZE,ndime,pgaus,*)
    real(rp),    intent(out)   :: gpsgi(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(out)   :: gpgve(VECTOR_SIZE,ndime,ndime,pgaus)
    real(rp),    intent(in)    :: gpvep(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: gpgrp(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(out)   :: gpst1(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: dtsgs(VECTOR_SIZE)
    real(rp),    intent(inout) :: resis_nsi(2,*)
    integer(ip), intent(inout) :: itsta_nsi(*)
    real(rp),    intent(out)   :: rmsgs_nsi
    real(rp),    intent(out)   :: resgs_nsi(*)
    real(rp),    intent(in)    :: gppor(VECTOR_SIZE,pgaus)
    integer(ip)                :: idime,igaus,inode,itsgs,jdime,ivect,jtsgs
    real(rp)                   :: resgs,gpnew(3),gpnor,gpnum,gpdnm,rels1
    real(rp)                   :: adv,dif,rea,gpvno,dummr
    real(rp)                   :: dj_fi(ndime,ndime)       ! d(f_i(u'))/du'_j   ! for NEWTON-R (kfl_sgsli_nsi==2)
    real(rp)                   :: djinv(ndime,ndime)       ! its inverse
    real(rp)                   :: funcu(ndime)             ! f(u')
    real(rp)                   :: deltu(ndime)             ! delta u'
    real(rp)                   :: auxi1,auxi2

#ifdef OPENACC
#define DEF_VECT 1:VECTOR_SIZE
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif

    if( kfl_stabi_nsi == NSI_GALERKIN .or. kfl_stabi_nsi == NSI_ALGEBRAIC_SPLIT_OSS ) return
    if( kfl_sgsti_nsi == 0 .and. kfl_sgsco_nsi == 0 ) return

    rels1 = 1.0_rp - relsg_nsi

    !----------------------------------------------------------------------
    !
    ! GPGVE: Velocity gradients
    !
    !----------------------------------------------------------------------

    gpgve = 0.0_rp
    do igaus = 1,pgaus
       do jdime = 1,ndime
          do idime = 1,ndime
             do inode = 1,pnode
                gpgve(DEF_VECT,idime,jdime,igaus) = gpgve(DEF_VECT,idime,jdime,igaus)&
                     + gpcar(DEF_VECT,idime,inode,igaus) * elvel(DEF_VECT,jdime,inode)
             end do
          end do
       end do
    end do

    !----------------------------------------------------------------------
    !
    ! GPSGI: Never changing residual (all without PICARD advection)
    !        Ri(u) = f - rho*u/dt - sig*u - grad(p) + div[2*mu*eps(u)]
    !        where f = (body forces) + rho*u^n/dt + rho*u'^n/dt
    !        For Split OSS: add time term of the SGS
    !
    !----------------------------------------------------------------------

    do igaus = 1,pgaus
       do idime = 1,ndime
          gpsgi(1:VECTOR_SIZE,idime,igaus) = gprhs(1:VECTOR_SIZE,idime,igaus) + gprhs_sgs(1:VECTOR_SIZE,idime,igaus) - gpgpr(1:VECTOR_SIZE,idime,igaus)

          do inode = 1,pnode
             gpsgi(1:VECTOR_SIZE,idime,igaus) = gpsgi(1:VECTOR_SIZE,idime,igaus)&
                  - rmom1(1:VECTOR_SIZE,inode,igaus) * elvel(1:VECTOR_SIZE,idime,inode)
          end do
       end do
    end do
    if( kfl_rmom2_nsi /= 0 ) then
       do igaus = 1,pgaus
          do idime = 1,ndime
             do inode = 1,pnode
                do jdime = 1,ndime
                   gpsgi(1:VECTOR_SIZE,idime,igaus) = gpsgi(1:VECTOR_SIZE,idime,igaus)&
                        - rmom2(1:VECTOR_SIZE,idime,jdime,inode,igaus) * elvel(1:VECTOR_SIZE,jdime,inode)
                end do
             end do
          end do
       end do
    end if

    if( kfl_stabi_nsi == 1 ) then
       !
       ! OSS: rho * u'^n/dt was not included in nsi_elmre3
       ! u ' <= u ' + u'^n / dt
       !
!!$     if( kfl_sgsti_nsi /= 0 ) then
!!$        do igaus = 1,pgaus
!!$           dummr = gpden(igaus) * dtsgs
!!$           do idime = 1,ndime
!!$              gpsgi(idime,igaus) = gpsgi(idime,igaus) + dummr * gpsgs(idime,igaus,2)
!!$           end do
!!$        end do
!!$     end if
    end if

    !----------------------------------------------------------------------
    !
    ! Solve subgrid scale equation usig PICARD (1) or NEWTON RAPHSON (2)
    !
    !----------------------------------------------------------------------

    if( kfl_sgsli_nsi == 0 ) then

       do igaus = 1,pgaus
          !
          ! tau1 = 1 / ( 4*mu/h^2 + 2*rho*u/h + rho*|w| + sig )
          !
          do ivect = 1,VECTOR_SIZE
             gpvno = sqrt(dot_product(gpvel(ivect,1:ndime,igaus),gpvel(ivect,1:ndime,igaus)))
             adv   = gpden(ivect,igaus)*gpvno                                 ! Convective term: rho*|u|
             dif   = gpvis(ivect,igaus)                                       ! Viscous term:    mu
             rea   = gpden(ivect,igaus)*corio_nsi + abs(gppor(ivect,igaus))   ! Coriolis: w + Porosity: sig
             call tauadr(&
                  kfl_taust_nsi,staco_nsi,adv,dif,rea,&
                  chale(ivect,1),chale(ivect,2),gpst1(ivect,igaus))
          end do

          do idime = 1,ndime
             do jdime = 1,ndime
                gpsgi(1:VECTOR_SIZE,idime,igaus) = gpsgi(1:VECTOR_SIZE,idime,igaus) - gpden(1:VECTOR_SIZE,igaus) &
                     * gpvel(1:VECTOR_SIZE,jdime,igaus) * gpgve(1:VECTOR_SIZE,jdime,idime,igaus)
             end do
          end do

       end do

       if( kfl_sgsti_nsi == 1 ) then
          !
          ! tau1' = 1 / ( rho/dt + 1/tau1 )
          !
          do igaus = 1,pgaus
             gpsp1(1:VECTOR_SIZE,igaus) = 1.0_rp / ( gpden(1:VECTOR_SIZE,igaus) * dtsgs(1:VECTOR_SIZE) + 1.0_rp / gpst1(1:VECTOR_SIZE,igaus) )
             do idime = 1,ndime
                gpsgs(1:VECTOR_SIZE,idime,igaus,1) = gpsp1(1:VECTOR_SIZE,igaus) * gpsgi(1:VECTOR_SIZE,idime,igaus)
             end do
          end do
       else
          do igaus = 1,pgaus
             do idime = 1,ndime
                gpsgs(1:VECTOR_SIZE,idime,igaus,1) = gpst1(1:VECTOR_SIZE,igaus) * gpsgi(1:VECTOR_SIZE,idime,igaus)
             end do
          end do
       end if

       gpadv(1:VECTOR_SIZE,1:ndime,1:pgaus) = gpvel(1:VECTOR_SIZE,1:ndime,1:pgaus) + gpsgs(1:VECTOR_SIZE,1:ndime,1:pgaus,1)


    else if( kfl_sgsli_nsi == 1 ) then

       do ivect = 1,VECTOR_SIZE

          do igaus = 1,pgaus

             itsgs = 0
             resgs = 1.e06_rp

             do while( itsgs < misgs_nsi .and. resgs > tosgs_nsi )
                itsgs = itsgs + 1
                jtsgs = min(mmsgs_nsi,itsgs)
                !
                ! tau1 = 1 / ( 4*mu/h^2 + 2*rho*uc/h + rho*|w| + sig )
                ! tau2 = 1 / [ (eps + 1/(h^2*4*mu/h^2 + 2*rho*uc/h + rho*|w|)) ]
                !
                gpvno = sqrt(dot_product(gpadv(ivect,1:ndime,igaus),gpadv(ivect,1:ndime,igaus)))
                adv   = gpden(ivect,igaus)*gpvno                                 ! Convective term: rho*|u+u'|
                dif   = gpvis(ivect,igaus)                                       ! Viscous term:    mu
                rea   = gpden(ivect,igaus)*corio_nsi + abs(gppor(ivect,igaus))   ! Coriolis: w + Porosity: sig
                call tauadr(&
                     kfl_taust_nsi,staco_nsi,adv,dif,rea,&
                     chale(ivect,1),chale(ivect,2),gpst1(ivect,igaus))
                !
                ! tau1' = 1 / ( rho/dt + 1/tau1 )
                !
                gpsp1(ivect,igaus) = 1.0_rp / ( gpden(ivect,igaus) * dtsgs(ivect) + 1.0_rp / gpst1(ivect,igaus) )
                !
                ! GPNEW: new subgrid scale
                ! 0. ASGS ........ u' = tau1' * [ Ri(u) - rho*(uc.grad)u ]
                ! 1. Full OSS .... u' = tau1' * [ Ri(u) - rho*(uc.grad)u ] + P
                ! 2. Split OSS ... u' = tau1' * [ - rho*(uc.grad)u  + rho*u'/dt ] + P
                !
                do idime = 1,ndime
                   gpnew(idime) = gpsgi(ivect,idime,igaus)
                   do jdime = 1,ndime
                      gpnew(idime) = gpnew(idime) - gpden(ivect,igaus)&
                           * gpadv(ivect,jdime,igaus) * gpgve(ivect,jdime,idime,igaus)
                   end do
                   gpnew(idime) = gpsp1(ivect,igaus) * gpnew(idime)
                end do
                if( kfl_stabi_nsi == 2 ) then
                   call runend('NSI_UPDSGS: NOT CODED')
!!$           else if( kfl_stabi_nsi /= 0 ) then
!!$              do idime = 1,ndime
!!$                 gpnew(idime) = gpnew(idime) + gpsp1(igaus) * gpvep(idime,igaus) / gpst1(igaus)
!!$              end do
                end if
                !
                ! RESGS: Residual and update subgrid scale
                ! GPSGS: New subgrid scale=GPNEW
                !
                gpnum = 0.0_rp
                gpdnm = 0.0_rp
                do idime = 1,ndime
                   gpnew(idime)               = relsg_nsi * gpnew(idime) + rels1 * gpsgs(ivect,idime,igaus,1)
                   dummr                      = gpsgs(ivect,idime,igaus,1) - gpnew(idime)
                   gpnum                      = gpnum + dummr * dummr
                   gpdnm                      = gpdnm + gpnew(idime) * gpnew(idime)
                   gpsgs(ivect,idime,igaus,1) = gpnew(idime)
                end do
                if( gpdnm /= 0.0_rp ) then
                   resgs = sqrt( gpnum/gpdnm )
                else
                   resgs = sqrt( gpnum )
                end if

                if( kfl_sgsco_nsi >= 1 ) then
                   !
                   ! GPADV: Update advection uc = u + u'
                   !
                   gpadv(ivect,1:ndime,igaus) = gpvel(ivect,1:ndime,igaus) + gpsgs(ivect,1:ndime,igaus,1)
                   !
                   ! RESIS_NSI: Inner residual
                   !
                   gpnor = sqrt(dot_product(gpnew(1:ndime),gpnew(1:ndime)))
                   resis_nsi(1,jtsgs) = resis_nsi(1,jtsgs) + gpnum * gpnum
                   resis_nsi(2,jtsgs) = resis_nsi(2,jtsgs) + gpnor * gpnor
                end if

             end do

             rmsgs_nsi = max(rmsgs_nsi,resgs)
             !
             ! ITSTA_NSI: Subgrid scale statistics
             !
             if( kfl_sgsco_nsi >= 1 ) itsta_nsi(jtsgs) = itsta_nsi(jtsgs) + 1

          end do

       end do

    else ! Newton Raphson

       do ivect = 1,VECTOR_SIZE

          do igaus = 1,pgaus
             itsgs = 0
             resgs = 1.e06_rp

             do while( itsgs < misgs_nsi .and. resgs > tosgs_nsi )
                itsgs = itsgs + 1
                jtsgs = min(mmsgs_nsi,itsgs)
                !
                ! tau1 = 1 / ( 4*mu/h^2 + 2*rho*uc/h + rho*|w| + sig )
                ! tau2 = 1 / [ (eps + 1/(h^2*4*mu/h^2 + 2*rho*uc/h + rho*|w|)) ]
                !
                gpvno = sqrt(dot_product(gpadv(ivect,1:ndime,igaus),gpadv(ivect,1:ndime,igaus)))
                adv   = gpden(ivect,igaus)*gpvno                           ! Convective term: rho*|u+u'|
                dif   = gpvis(ivect,igaus)                                 ! Viscous term:    mu
                rea   = abs(gppor(ivect,igaus))                            ! Porosity term:   sig
                call tauadr(&
                     kfl_taust_nsi,staco_nsi,adv,dif,rea,&
                     chale(ivect,1),chale(ivect,2),gpst1(ivect,igaus))
                !
                ! tau1' = 1 / ( rho/dt + 1/tau1 )
                !
                gpsp1(ivect,igaus) = 1.0_rp / ( gpden(ivect,igaus) * dtsgs(ivect) + 1.0_rp / gpst1(ivect,igaus) )
                !
                ! obtain d(f_i(u'))/du'_j and its inverse
                !
                if ( gpvno>1.0d-8 ) then
                   auxi1 = 2.0_rp * gpden(ivect,igaus) / ( chale(ivect,1) * gpvno )     ! C2 * rho / ( h * | u + u'| )
                else
                   auxi1 =  0.0_rp
                end if
                do idime = 1,ndime
                   do jdime = 1,ndime
                      dj_fi(idime,jdime) = &
                           + auxi1 * gpadv(ivect,jdime,igaus) * gpsgs(ivect,idime,igaus,1) &  ! auxi1 * ( u_jd + u'_jd ) * u'_id
                           + gpden(ivect,igaus) * gpgve(ivect,jdime,idime,igaus)              ! rho * du_id / dx_jd
                   end do
                end do

                auxi1 = 1.0_rp / gpsp1(ivect,igaus)
                do idime=1,ndime
                   dj_fi(idime,idime) = dj_fi(idime,idime) + auxi1
                end do
                call invmtx(dj_fi,djinv,auxi2,ndime)
                !
                ! obtain f(u')
                !
                do idime=1,ndime
                   funcu(idime) = - gpsgi(ivect,idime,igaus)  &  ! - non changing part of the residual - rho u'_n /(tita*dt) -- See before
                        + auxi1 * gpsgs(ivect,idime,igaus,1)     ! + ( inv(tau1') * u'
                   do jdime=1,ndime
                      funcu(idime) = funcu(idime) + gpden(ivect,igaus) * gpadv(ivect,jdime,igaus) * gpgve(ivect,jdime,idime,igaus)    ! rho * ( u_jd + u'_jd ) * dv_i / dx_j (changing part of residual)
                   end do
                end do
                if( kfl_stabi_nsi == 2 ) then
                   call runend('NSI_UPDSGS: NOT CODED')
!!$           else if( kfl_stabi_nsi /= 0 ) then
!!$              do idime = 1,ndime
!!$                 funcu(idime) = funcu(idime) - gpvep(idime,igaus) / gpst1(igaus)
!!$              end do
                end if
                !
                ! from A_ij * deltau'_j =  - f(u')_i ; (where A_ij = d(f_i(u'))/du'_j )
                ! obtain deltau'_j = - invA_ji *  f(u')_i
                !
                do jdime=1,ndime
                   deltu(jdime) = 0.0_rp
                   do idime=1,ndime
                      deltu(jdime) = deltu(jdime) - djinv(jdime,idime) * funcu(idime)
                   end do
                end do
                !
                ! RESGS: Residual and update subgrid scale
                ! GPSGS: New subgrid scale=GPNEW
                !
                gpnum = 0.0_rp
                gpdnm = 0.0_rp
                do idime = 1,ndime
                   gpnew(idime)               = gpsgs(ivect,idime,igaus,1) + relsg_nsi * deltu(idime)
                   dummr                      = deltu(idime)
                   gpnum                      = gpnum + dummr * dummr
                   gpdnm                      = gpdnm + gpnew(idime) * gpnew(idime)
                   gpsgs(ivect,idime,igaus,1) = gpnew(idime)
                end do
                if( gpdnm /= 0.0_rp ) then
                   resgs = sqrt( gpnum/gpdnm )
                else
                   resgs = sqrt( gpnum )
                end if

                if( kfl_sgsco_nsi >= 1 ) then
                   !
                   ! GPADV: Update advection uc = u + u'
                   !
                   do idime = 1,ndime
                      gpadv(ivect,idime,igaus) = gpvel(ivect,idime,igaus) + gpsgs(ivect,idime,igaus,1)
                   end do
                   !
                   ! RESIS_NSI: Inner residual
                   !
                   gpnor = sqrt(dot_product(gpnew(1:ndime),gpnew(1:ndime)))
                   resis_nsi(1,jtsgs) = resis_nsi(1,jtsgs) + gpnum * gpnum
                   resis_nsi(2,jtsgs) = resis_nsi(2,jtsgs) + gpnor * gpnor
                end if

             end do

             rmsgs_nsi = max(rmsgs_nsi,resgs)
             !
             ! ITSTA_NSI: Subgrid scale statistics
             !
             if( kfl_sgsco_nsi >= 1 ) itsta_nsi(jtsgs) = itsta_nsi(jtsgs) + 1

          end do

       end do

    end if


  end subroutine nsi_element_subgrid_scale

  subroutine nsi_element_stabilization(&
       pgaus,pnode,chale,hleng,gpadv,gpvis,gpden,gpcar,gpst1,&
       gpst2,gpsp1,gpsp2,gptt1,gptt2,rmom1,gppor,dtsgs,&
       dtinv_loc,tamin_loc,tamax_loc)
    !------------------------------------------------------------------------
    !****f* Nastin/nsi_elmsgs
    ! NAME
    !    nsi_elmsgs
    ! DESCRIPTION
    !    1. Compute stability parameters
    !       - tau1 [ L/(rho*U) ]
    !       - tau2 [ rho*L*U ]
    !    3. Update momentum residual
    ! USES
    ! USED BY
    !***
    !------------------------------------------------------------------------
    use def_kintyp, only       :  ip,rp
    use def_domain, only       :  ndime,mnode
    use def_nastin, only       :  staco_nsi,kfl_sgsti_nsi,&
         &                        kfl_taust_nsi,kfl_stabi_nsi,&
         &                        kfl_ellen_nsi,corio_nsi,&
         &                        kfl_stabi_nsi
    use mod_tauadr, only       :  tauadr
    implicit none
    integer(ip), intent(in)    :: pgaus,pnode
    real(rp),    intent(in)    :: chale(VECTOR_SIZE,2)
    real(rp),    intent(in)    :: hleng(VECTOR_SIZE,ndime)
    real(rp),    intent(in)    :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)
    real(rp),    intent(in)    :: gpadv(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: gpvis(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpden(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gppor(VECTOR_SIZE,pgaus)
    real(rp),    intent(out)   :: gpst1(VECTOR_SIZE,pgaus),gpst2(VECTOR_SIZE,pgaus)
    real(rp),    intent(out)   :: gpsp1(VECTOR_SIZE,pgaus),gpsp2(VECTOR_SIZE,pgaus)
    real(rp),    intent(out)   :: gptt1(VECTOR_SIZE,pgaus),gptt2(VECTOR_SIZE,pgaus)
    real(rp),    intent(out)   :: rmom1(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(in)    :: dtsgs(VECTOR_SIZE)
    real(rp),    intent(in)    :: dtinv_loc(VECTOR_SIZE)
    real(rp),    intent(inout) :: tamin_loc
    real(rp),    intent(inout) :: tamax_loc
    integer(ip)                :: idime,igaus,inode,ivect
    real(rp)                   :: adv,dif,rea,h2,gpvno

    if( kfl_stabi_nsi == NSI_GALERKIN ) then

       gpst1 = 0.0_rp
       gpst2 = 0.0_rp
       gpsp1 = 0.0_rp
       gpsp2 = 0.0_rp
       gptt1 = 0.0_rp
       gptt2 = 0.0_rp

    else

       !----------------------------------------------------------------------
       !
       ! TAU1 and TAU2
       !
       !----------------------------------------------------------------------

       do igaus = 1,pgaus
          do ivect = 1,VECTOR_SIZE
             !
             ! tau1 = 1 / ( 4*mu/h^2 + 2*rho*uc/h + rho*|w| + sig )
             ! tau2 = 1 / [ (eps + 1/(h^2*4*mu/h^2 + 2*rho*uc/h + rho*|w|)) ]
             !
             gpvno = sqrt(dot_product(gpadv(ivect,1:ndime,igaus),gpadv(ivect,1:ndime,igaus)))
             adv   = gpden(ivect,igaus)*gpvno                                         ! Convective term: rho*|u+u'|
             dif   = gpvis(ivect,igaus)                                               ! Viscous term:    mu
             rea   = gpden(ivect,igaus)*corio_nsi + abs(gppor(ivect,igaus))           ! Coriolis: w + Porosity: sig
             h2    = chale(ivect,2) * chale(ivect,2)
             call tauadr(&
                  kfl_taust_nsi,staco_nsi,adv,dif,rea,&
                  chale(ivect,1),chale(ivect,2),gpst1(ivect,igaus),gpden(ivect,igaus)*dtinv_loc(ivect))
             if( gpst1(ivect,igaus) /= 0.0_rp ) then
                if( kfl_ellen_nsi == 6 ) then   ! MIXLE
                   gpst2(ivect,igaus) = staco_nsi(4)*(4.0_rp*dif + 2.0_rp*adv*hleng(ivect,1))    ! Use maximum elem length only for tau2
                else
                   gpst2(ivect,igaus) = staco_nsi(4)*h2/gpst1(ivect,igaus)
                end if
             else
                gpst2(ivect,igaus) = 0.0_rp
             end if
          end do
       end do

       do igaus = 1,pgaus
          !
          ! (tau1',tau2') and (tau1'/tau1,tau2'/tau2)
          !
          if( kfl_sgsti_nsi == 1 ) then
             gpsp1(1:VECTOR_SIZE,igaus) = 1.0_rp / ( gpden(1:VECTOR_SIZE,igaus) * dtsgs(1:VECTOR_SIZE) + 1.0_rp / gpst1(1:VECTOR_SIZE,igaus) )
             gpsp2(1:VECTOR_SIZE,igaus) = gpst2(1:VECTOR_SIZE,igaus)
             gptt1(1:VECTOR_SIZE,igaus) = gpsp1(1:VECTOR_SIZE,igaus) / gpst1(1:VECTOR_SIZE,igaus)
             gptt2(1:VECTOR_SIZE,igaus) = 1.0_rp
          else
             gpsp1(1:VECTOR_SIZE,igaus) = gpst1(1:VECTOR_SIZE,igaus)
             gpsp2(1:VECTOR_SIZE,igaus) = gpst2(1:VECTOR_SIZE,igaus)
             gptt1(1:VECTOR_SIZE,igaus) = 1.0_rp
             gptt2(1:VECTOR_SIZE,igaus) = 1.0_rp
          end if
          if( kfl_stabi_nsi == 1 ) gptt1(1:VECTOR_SIZE,igaus) = 1.0_rp

       end do

       !----------------------------------------------------------------------
       !
       ! RMOMU: Add advection to residual; split OSS needs only advection
       !
       !----------------------------------------------------------------------

       do igaus = 1,pgaus
          do inode = 1,pnode
             do idime = 1,ndime
                rmom1(1:VECTOR_SIZE,inode,igaus) = rmom1(1:VECTOR_SIZE,inode,igaus)&
                     + gpden(1:VECTOR_SIZE,igaus) * gpadv(1:VECTOR_SIZE,idime,igaus)&
                     * gpcar(1:VECTOR_SIZE,idime,inode,igaus)
             end do
          end do
       end do

       !----------------------------------------------------------------------
       !
       ! TAMIN_NSI, TAMAX_NSI: Minimum and maximum tau
       !
       !----------------------------------------------------------------------

       tamax_loc = max(tamax_loc,maxval(gpsp1))
       tamin_loc = min(tamin_loc,minval(gpsp1))

    end if

  end subroutine nsi_element_stabilization

  subroutine nsi_element_dirichlet(&
       pnode,pevat,list_elements,lnods,elauu,elaup,&
       elapp,elapu,elrbu,elrbp,elcmm)
    !-----------------------------------------------------------------------
    !****f* Nastin/nsi_elmdi3
    ! NAME
    !    nsi_elmdi3
    ! DESCRIPTION
    ! USES
    !    nsi_rotma3
    ! USED BY
    !    nsi_elmop3
    !***
    !-----------------------------------------------------------------------
    use def_kintyp, only       :  ip,rp
    use def_domain, only       :  ndime,lpoty
    use def_nastin, only       :  kfl_confi_nsi,nodpr_nsi,kfl_local_nsi,&
         &                        kfl_fixno_nsi,kfl_fixrs_nsi,bvess_nsi,&
         &                        valpr_nsi,&
         &                        kfl_imppr_nsi,kfl_fixpr_nsi,kfl_stabi_nsi

    use def_kermod, only       :  kfl_adj_prob
    use mod_local_basis, only  :  local_basis_matrix
    implicit none
    integer(ip), intent(in)    :: pnode
    integer(ip), intent(in)    :: pevat
    integer(ip), intent(in)    :: list_elements(VECTOR_SIZE)
    integer(ip), intent(in)    :: lnods(VECTOR_SIZE,pnode)
    real(rp),    intent(inout) :: elauu(VECTOR_SIZE,pevat,pevat)
    real(rp),    intent(inout) :: elaup(VECTOR_SIZE,pevat,pnode)
    real(rp),    intent(inout) :: elapp(VECTOR_SIZE,pnode,pnode)
    real(rp),    intent(inout) :: elapu(VECTOR_SIZE,pnode,pevat)
    real(rp),    intent(inout) :: elrbu(VECTOR_SIZE,pevat)
    real(rp),    intent(inout) :: elrbp(VECTOR_SIZE,pnode)
    real(rp),    intent(inout) :: elcmm(VECTOR_SIZE,pevat,pevat)
    real(rp)                   :: adiag,worma(9),vimpr_nsi,adia2
    real(rp)                   :: rotma(ndime,ndime)
    integer(ip)                :: inode,ipoin,ievav,idime
    integer(ip)                :: jevav,iroty,ibopo,jnode
    integer(ip)                :: ivect
    !
    ! Prescribe one pressure degree of freedom if the flow is confined
    !
    if( kfl_confi_nsi >= 0 .and. nodpr_nsi > 0 ) then

       do ivect = 1,VECTOR_SIZE
          if( list_elements(ivect) > 0 ) then
             do inode = 1,pnode
                ipoin = lnods(ivect,inode)
                if( ipoin == nodpr_nsi ) then
                   if( kfl_stabi_nsi == NSI_GALERKIN .or. kfl_stabi_nsi == NSI_ALGEBRAIC_SPLIT_OSS ) then
                      adiag = 1.0_rp
                   else
                      adiag = elapp(ivect,inode,inode)
                   end if
                   do jevav = 1,pnode * ndime
                      elrbu(ivect,jevav)       = elrbu(ivect,jevav)-elaup(ivect,jevav,inode) * valpr_nsi
                      elaup(ivect,jevav,inode) = 0.0_rp
                      elapu(ivect,inode,jevav) = 0.0_rp
                   end do
                   do jnode = 1,pnode
                      elrbp(ivect,jnode)       = elrbp(ivect,jnode)-elapp(ivect,jnode,inode) * valpr_nsi
                      elapp(ivect,jnode,inode) = 0.0_rp
                      elapp(ivect,inode,jnode) = 0.0_rp
                   end do
                   elapp(ivect,inode,inode) = adiag
                   elrbp(ivect,inode)       = valpr_nsi*adiag
                end if
             end do
          end if
       end do
    end if
    !
    ! Impose pressure at nodes with kfl_fixpr_nsi - We are modifying the monolithic problem
    !
    if( kfl_imppr_nsi > 0 ) then
       vimpr_nsi = 0.0_rp
       do ivect = 1,VECTOR_SIZE
          if( list_elements(ivect) > 0 ) then
             do inode = 1,pnode
                ipoin = lnods(ivect,inode)
                if ( kfl_fixpr_nsi(1,ipoin) > 0 ) then
                   if( kfl_stabi_nsi == NSI_GALERKIN .or. kfl_stabi_nsi == NSI_ALGEBRAIC_SPLIT_OSS ) then
                      adiag = elapp(ivect,inode,inode)
                   else
                      adiag = 1.0_rp
                   end if
                   do jevav = 1,pnode * ndime
                      elrbu(ivect,jevav)       = elrbu(ivect,jevav)-elaup(ivect,jevav,inode) * vimpr_nsi
                      elaup(ivect,jevav,inode) = 0.0_rp
                      elapu(ivect,inode,jevav) = 0.0_rp
                   end do
                   do jnode = 1,pnode
                      elrbp(ivect,jnode)       = elrbp(ivect,jnode)-elapp(ivect,jnode,inode) * vimpr_nsi
                      elapp(ivect,jnode,inode) = 0.0_rp
                      elapp(ivect,inode,jnode) = 0.0_rp
                   end do
                   elapp(ivect,inode,inode) = adiag
                   elrbp(ivect,inode)       = vimpr_nsi*adiag
                end if
             end do
          end if
       end do
    end if
    !
    ! Rotate the nodal matrices, the element nodal velocities and RHS to
    ! prescribe boundary conditions in a skew-system, either the tangent
    ! one or another prescribed by the user. Also, periodical boundary
    ! conditions are accounted for.
    !
    if( kfl_local_nsi == 1 ) then
       do ivect = 1,VECTOR_SIZE
          if( list_elements(ivect) > 0 ) then
             do inode = 1,pnode
                ipoin = lnods(ivect,inode)
                iroty = kfl_fixrs_nsi(ipoin)
                ibopo = lpoty(ipoin)
                if( iroty /= 0 ) then
                   call local_basis_matrix(ipoin,ibopo,iroty,rotma)
                   call nsi_rotma3(&
                        &      inode,pnode,pevat,elauu(ivect,:,:),elaup(ivect,:,:),&
                        &      elapu(ivect,:,:),elrbu(ivect,:),elcmm(ivect,:,:),rotma,worma)
                end if
             end do
          end if
       end do
    end if
    !
    ! Velocity Dirichlet boundary conditions
    !
    do ivect = 1,VECTOR_SIZE
       if( list_elements(ivect) > 0 ) then

          do inode = 1,pnode
             ievav = (inode-1) * ndime
             ipoin = lnods(ivect,inode)
             do idime = 1,ndime
                ievav = ievav+1
                if(      kfl_fixno_nsi(idime,ipoin) ==  1 &
                     .or.kfl_fixno_nsi(idime,ipoin) ==  3 &
                     .or.kfl_fixno_nsi(idime,ipoin) ==  8 &
                     .or.kfl_fixno_nsi(idime,ipoin) ==  9 &
                     .or.kfl_fixno_nsi(idime,ipoin) ==  5 &
                     .or.kfl_fixno_nsi(idime,ipoin) ==  6 &
                     .or.kfl_fixno_nsi(idime,ipoin) ==  7 &
                     .or.kfl_fixno_nsi(idime,ipoin) == 11 ) then
                   adiag = elauu(ivect,ievav,ievav)
                   adia2 = elcmm(ivect,ievav,ievav)
                   do jevav = 1,pevat
                      elauu(ivect,ievav,jevav) = 0.0_rp
                      elcmm(ivect,ievav,jevav) = 0.0_rp
                   end do
                   do jevav = 1,pevat
                      elrbu(ivect,jevav)       = elrbu(ivect,jevav) - elauu(ivect,jevav,ievav) * bvess_nsi(idime,ipoin,1)
                      elauu(ivect,jevav,ievav) = 0.0_rp
                      ! Since nothing is sent to the rhs - This implies that the Consistent mass matrix is only ready for cases were
                      ! the value to be prescribed is 0. This is what happens for the end of step velocity since one solves for an increment
                      elcmm(ivect,ievav,jevav) = 0.0_rp
                   end do
                   do jnode = 1,pnode
                      elaup(ivect,ievav,jnode) = 0.0_rp
                      if (kfl_adj_prob == 1) elapu(ivect,jnode,ievav) = 0.0_rp     ! added for the adjoint problem
                   end do
                   elauu(ivect,ievav,ievav) = adiag
                   elrbu(ivect,ievav)       = adiag * bvess_nsi(idime,ipoin,1)
                   elcmm(ivect,ievav,ievav) = adiag
                end if
             end do
          end do
       end if
    end do

  end subroutine nsi_element_dirichlet

  subroutine nsi_element_schur(&
       pnode,pgaus,list_elements,lnods,gpcar,gpvol,gpden,&
       gpvis,gppor,gpsha,elvel,chale,gpsp1,elapp,elmap,&
       dtinv_loc)
    !----------------------------------------------------------------------
    !****f* Nastin/nsi_elmsch
    ! NAME
    !    nsi_elmsch
    ! DESCRIPTION
    !    Compute the Schur complement preconditioner
    ! USES
    ! USED BY
    !***
    !----------------------------------------------------------------------
    use def_kintyp, only     :  ip,rp
    use def_nastin, only     :  nodpr_nsi,&
         &                      kfl_confi_nsi,kfl_fixpr_nsi,&
         &                      staco_nsi,kfl_taush_nsi,&
         &                      kfl_predi_nsi,kfl_matdi_nsi,&
         &                      pabdf_nsi,kfl_regim_nsi,kfl_confi_nsi
    use def_domain, only     :  mnode,ndime,lpoty
    use mod_tauadr, only     :  tauadr
    implicit none
    integer(ip), intent(in)  :: pnode
    integer(ip), intent(in)  :: pgaus
    integer(ip), intent(in)  :: list_elements(VECTOR_SIZE)
    integer(ip), intent(in)  :: lnods(VECTOR_SIZE,pnode)
    real(rp),    intent(in)  :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)
    real(rp),    intent(in)  :: gpvol(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)  :: gpden(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)  :: gpvis(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)  :: gppor(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)  :: gpsha(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(in)  :: elvel(VECTOR_SIZE,ndime,mnode)
    real(rp),    intent(in)  :: chale(VECTOR_SIZE,2)
    real(rp),    intent(in)  :: gpsp1(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)  :: elapp(VECTOR_SIZE,pnode,pnode)
    real(rp),    intent(out) :: elmap(VECTOR_SIZE,pnode,pnode)
    real(rp),    intent(in)  :: dtinv_loc(VECTOR_SIZE)
    integer(ip)              :: inode,jnode,kdime,igaus,ipoin,ibopo
    integer(ip)              :: idime,ivect
!    integer(ip)              :: ihang
    real(rp)                 :: fact1(VECTOR_SIZE),fact2(VECTOR_SIZE)
    real(rp)                 :: fact3(VECTOR_SIZE),fact4(VECTOR_SIZE),fact0(VECTOR_SIZE)
    real(rp)                 :: tau(VECTOR_SIZE),penal(VECTOR_SIZE,pgaus)
    real(rp)                 :: adv(VECTOR_SIZE),dif(VECTOR_SIZE),rea(VECTOR_SIZE)
    real(rp)                 :: gpadv(VECTOR_SIZE,3),gpvno(VECTOR_SIZE),facts
    real(rp)                 :: tauinv(VECTOR_SIZE)

    !----------------------------------------------------------------------
    !
    ! Initialization
    !
    !----------------------------------------------------------------------

    elmap = 0.0_rp
    penal = 0.0_rp

    !!DMM-LM #ifdef matiaslma
    if (kfl_regim_nsi==3 .and. kfl_confi_nsi == 1) then
       penal(1:VECTOR_SIZE,1:pgaus) = 1.0e-4_rp*gpden(1:VECTOR_SIZE,1:pgaus)/gpvis(1:VECTOR_SIZE,1:pgaus)
    end if
    !!DMM-LM #endif

    !----------------------------------------------------------------------
    !
    ! Elemental matrix
    !
    !----------------------------------------------------------------------

    if( kfl_predi_nsi == 2 ) then
       !
       ! P (+ App) = ( (dt/ (rho*pabdf_nsi(1)) )*grad p , grad q ) (+ App)
       !
       do ivect = 1,VECTOR_SIZE
          if( dtinv_loc(ivect) /= 0.0_rp ) then
             fact2(ivect) = 1.0_rp/(dtinv_loc(ivect)*pabdf_nsi(1))
          else
             fact2(ivect) = 0.0_rp
          end if
       end do

       do igaus = 1,pgaus
          fact3(1:VECTOR_SIZE) = gpvol(1:VECTOR_SIZE,igaus)*( fact2(1:VECTOR_SIZE) / gpden(1:VECTOR_SIZE,igaus) + gpsp1(1:VECTOR_SIZE,igaus) )
          do inode = 1,pnode
             fact0 = penal(1:VECTOR_SIZE,igaus)*gpsha(1:VECTOR_SIZE,inode,igaus)*gpvol(1:VECTOR_SIZE,igaus)
             do jnode = inode+1,pnode
                fact1(1:VECTOR_SIZE) = 0.0_rp
                do kdime = 1,ndime
                   fact1(1:VECTOR_SIZE) = fact1(1:VECTOR_SIZE) + gpcar(1:VECTOR_SIZE,kdime,inode,igaus)&
                        &        * gpcar(1:VECTOR_SIZE,kdime,jnode,igaus)
                end do
                fact1(1:VECTOR_SIZE) = fact1(1:VECTOR_SIZE)*fact3(1:VECTOR_SIZE)
                elmap(1:VECTOR_SIZE,inode,jnode) = elmap(1:VECTOR_SIZE,inode,jnode)+fact1 + fact0(1:VECTOR_SIZE) * gpsha(1:VECTOR_SIZE,jnode,igaus)
                elmap(1:VECTOR_SIZE,jnode,inode) = elmap(1:VECTOR_SIZE,jnode,inode)+fact1 + fact0(1:VECTOR_SIZE) * gpsha(1:VECTOR_SIZE,jnode,igaus)
             end do
             fact1(1:VECTOR_SIZE) = 0.0_rp
             do kdime = 1,ndime
                fact1(1:VECTOR_SIZE) = fact1(1:VECTOR_SIZE) + gpcar(1:VECTOR_SIZE,kdime,inode,igaus)&
                     &        * gpcar(1:VECTOR_SIZE,kdime,inode,igaus)
             end do
             fact1(1:VECTOR_SIZE) = fact1(1:VECTOR_SIZE) * fact3(1:VECTOR_SIZE)
             elmap(1:VECTOR_SIZE,inode,inode) = elmap(1:VECTOR_SIZE,inode,inode) + fact1(1:VECTOR_SIZE) + fact0(1:VECTOR_SIZE) * gpsha(1:VECTOR_SIZE,inode,igaus)
          end do
       end do

    else if( kfl_predi_nsi == 3 ) then
       !
       ! P (+ App) = ( tau'*grad p , grad q ) (+ App)
       !
       do igaus = 1,pgaus

          gpadv(1:VECTOR_SIZE,:) = 0.0_rp
          do inode = 1,pnode
             do idime = 1,ndime
                gpadv(1:VECTOR_SIZE,idime) = gpadv(1:VECTOR_SIZE,idime) + gpsha(1:VECTOR_SIZE,inode,igaus) * elvel(1:VECTOR_SIZE,idime,inode)
             end do
          end do
          gpvno(1:VECTOR_SIZE) = 0.0_rp
          do idime = 1,ndime
             gpvno(1:VECTOR_SIZE) = gpvno(1:VECTOR_SIZE) + gpadv(1:VECTOR_SIZE,idime) * gpadv(1:VECTOR_SIZE,idime)
          end do
          gpvno(1:VECTOR_SIZE) = sqrt(gpvno(1:VECTOR_SIZE))
          adv(1:VECTOR_SIZE)   = gpvno(1:VECTOR_SIZE) * gpden(1:VECTOR_SIZE,igaus)    ! Convective term rho * u
          dif(1:VECTOR_SIZE)   = gpvis(1:VECTOR_SIZE,igaus)                           ! Viscous term mu
          rea(1:VECTOR_SIZE)   = gppor(1:VECTOR_SIZE,igaus)                           ! Reaction term sig
          !
          !             1
          ! M^-1 = -----------
          !        rho     1
          !        --- +  ---
          !         dt    tau
          !
          do ivect = 1,VECTOR_SIZE
             call tauadr(&
                  kfl_taush_nsi,staco_nsi,adv(ivect),dif(ivect),rea(ivect),&
                  chale(ivect,1),chale(ivect,2),tau(ivect))
          end do
          tauinv(1:VECTOR_SIZE) = 1.0_rp / max(zeror,tau(1:VECTOR_SIZE))
          fact2(1:VECTOR_SIZE)  = 1.0_rp / ( gpden(1:VECTOR_SIZE,igaus) * dtinv_loc(1:VECTOR_SIZE) * pabdf_nsi(1) + tauinv(1:VECTOR_SIZE) )
          !
          ! P + App <= tau + M^-1
          !
          fact3(1:VECTOR_SIZE) = ( fact2(1:VECTOR_SIZE) + gpsp1(1:VECTOR_SIZE,igaus) ) * gpvol(1:VECTOR_SIZE,igaus)

          do inode = 1,pnode
             fact0(1:VECTOR_SIZE) = penal(1:VECTOR_SIZE,igaus) * gpsha(1:VECTOR_SIZE,inode,igaus) * gpvol(1:VECTOR_SIZE,igaus)
             do jnode = inode+1,pnode
                fact1(1:VECTOR_SIZE) = 0.0_rp
                do kdime = 1,ndime
                   fact1(1:VECTOR_SIZE) = fact1(1:VECTOR_SIZE) + gpcar(1:VECTOR_SIZE,kdime,inode,igaus)&
                        &                                      * gpcar(1:VECTOR_SIZE,kdime,jnode,igaus)
                end do
                fact1(1:VECTOR_SIZE) = fact1(1:VECTOR_SIZE) * fact3(1:VECTOR_SIZE)
                elmap(1:VECTOR_SIZE,inode,jnode) = elmap(1:VECTOR_SIZE,inode,jnode) + fact1(1:VECTOR_SIZE)  + fact0(1:VECTOR_SIZE) * gpsha(1:VECTOR_SIZE,jnode,igaus)
                elmap(1:VECTOR_SIZE,jnode,inode) = elmap(1:VECTOR_SIZE,jnode,inode) + fact1(1:VECTOR_SIZE)  + fact0(1:VECTOR_SIZE) * gpsha(1:VECTOR_SIZE,jnode,igaus)
             end do
             fact1(1:VECTOR_SIZE) = 0.0_rp
             do kdime = 1,ndime
                fact1(1:VECTOR_SIZE) = fact1(1:VECTOR_SIZE) + gpcar(1:VECTOR_SIZE,kdime,inode,igaus) * gpcar(1:VECTOR_SIZE,kdime,inode,igaus)
             end do
             fact1(1:VECTOR_SIZE) = fact1(1:VECTOR_SIZE) * fact3(1:VECTOR_SIZE)
             elmap(1:VECTOR_SIZE,inode,inode) = elmap(1:VECTOR_SIZE,inode,inode) + fact1(1:VECTOR_SIZE)  + fact0(1:VECTOR_SIZE) * gpsha(1:VECTOR_SIZE,inode,igaus)
          end do
       end do
       
    else if( kfl_predi_nsi == 4 ) then
       !
       ! P (+ App) = mu*M (+ App)
       !
       do igaus = 1,pgaus

          fact3(1:VECTOR_SIZE) = gpsp1(1:VECTOR_SIZE,igaus)*gpvol(1:VECTOR_SIZE,igaus)
          fact4(1:VECTOR_SIZE) = gpvol(1:VECTOR_SIZE,igaus)/gpvis(1:VECTOR_SIZE,igaus)

          do inode=1,pnode
             fact0(1:VECTOR_SIZE) = penal(1:VECTOR_SIZE,igaus)*gpsha(1:VECTOR_SIZE,inode,igaus)*gpvol(1:VECTOR_SIZE,igaus)
             do jnode=inode+1,pnode
                fact1(1:VECTOR_SIZE)=0.0_rp
                do kdime=1,ndime
                   fact1(1:VECTOR_SIZE) = fact1(1:VECTOR_SIZE) + gpcar(1:VECTOR_SIZE,kdime,inode,igaus)&
                        &                                      * gpcar(1:VECTOR_SIZE,kdime,jnode,igaus)
                end do
                fact1(1:VECTOR_SIZE) = fact1(1:VECTOR_SIZE) * fact3(1:VECTOR_SIZE) + gpsha(1:VECTOR_SIZE,inode,igaus) * gpsha(1:VECTOR_SIZE,jnode,igaus) * fact4(1:VECTOR_SIZE)
                elmap(1:VECTOR_SIZE,inode,jnode) = elmap(1:VECTOR_SIZE,inode,jnode) + fact1(1:VECTOR_SIZE) + fact0(1:VECTOR_SIZE) * gpsha(1:VECTOR_SIZE,jnode,igaus)
                elmap(1:VECTOR_SIZE,jnode,inode) = elmap(1:VECTOR_SIZE,jnode,inode) + fact1(1:VECTOR_SIZE) + fact0(1:VECTOR_SIZE) * gpsha(1:VECTOR_SIZE,jnode,igaus)
             end do
             fact1(1:VECTOR_SIZE) = 0.0_rp
             do kdime = 1,ndime
                fact1(1:VECTOR_SIZE) = fact1(1:VECTOR_SIZE) + gpcar(1:VECTOR_SIZE,kdime,inode,igaus)&
                     &                                      * gpcar(1:VECTOR_SIZE,kdime,inode,igaus)
             end do
             fact1(1:VECTOR_SIZE) = fact1(1:VECTOR_SIZE) * fact3(1:VECTOR_SIZE) + gpsha(1:VECTOR_SIZE,inode,igaus) * gpsha(1:VECTOR_SIZE,inode,igaus) * fact4(1:VECTOR_SIZE)
             elmap(1:VECTOR_SIZE,inode,inode) = elmap(1:VECTOR_SIZE,inode,inode) + fact1(1:VECTOR_SIZE) + fact0(1:VECTOR_SIZE) * gpsha(1:VECTOR_SIZE,inode,igaus)
          end do
       end do

    else if( kfl_predi_nsi == 7 .or. kfl_predi_nsi == 8 .or. kfl_predi_nsi == 9 ) then
       !
       ! "Pure" Laplacian (to be scaled later on)
       !
       do igaus = 1,pgaus
          do jnode = 1,pnode
             do inode = 1,pnode
                do idime = 1,ndime
                   elmap(1:VECTOR_SIZE,inode,jnode) = elmap(1:VECTOR_SIZE,inode,jnode) &
                        + gpvol(1:VECTOR_SIZE,igaus) &
                        * gpcar(1:VECTOR_SIZE,idime,inode,igaus) &
                        * gpcar(1:VECTOR_SIZE,idime,jnode,igaus)
                end do
             end do
          end do
       end do

    end if

    !----------------------------------------------------------------------
    !
    ! Prescribe boundary conditions
    !
    !----------------------------------------------------------------------

    if( kfl_matdi_nsi <= 0 .or. kfl_matdi_nsi == 2 ) then

       do ivect = 1,VECTOR_SIZE

          if( list_elements(ivect) > 0 ) then

             do inode = 1,pnode
                !
                ! Neumann nodes
                !
                ipoin = lnods(ivect,inode)
                ibopo = lpoty(ipoin)
                if( ibopo /= 0 ) then
                   if( kfl_fixpr_nsi(1,ipoin) > 0 ) then
                      facts = elmap(ivect,inode,inode)
                      do jnode = 1,pnode
                         elmap(ivect,jnode,inode) = 0.0_rp    ! Beware perhaps this line should be commented out - at least that is what
                         ! I do for Takuji. This leads to a non sym matrix but at least values are
                         ! not wiped out without sending them tp RHS - THIS is something we must anlyse
                         ! further with guillaume
                         elmap(ivect,inode,jnode) = 0.0_rp
                      end do
                      elmap(ivect,inode,inode) = facts
                   end if
                end if
             end do

             if(kfl_confi_nsi==1.or.kfl_confi_nsi==-3) then
                !
                ! Confined flow and periodicity
                !
                do inode=1,pnode
                   ipoin=lnods(ivect,inode)
                   if(ipoin==nodpr_nsi) then
                      facts=elmap(ivect,inode,inode)
                      do jnode=1,pnode
                         elmap(ivect,jnode,inode)=0.0_rp
                         elmap(ivect,inode,jnode)=0.0_rp
                      end do
                      elmap(ivect,inode,inode)=facts
                   end if
                end do

             end if

             !if(nhang<0) then
             !   !
             !   ! Hanging nodes
             !   !
             !   do inode=1,pnode
             !      ipoin=lnods(inode)
             !      ihang=0
             !      do while(ihang<nhang)
             !         ihang=ihang+1
             !         if(lhang(1,ihang)==ipoin) then
             !            ihang=nhang
             !            facts=elmap(inode,inode)
             !            do jnode=1,pnode
             !               elmap(jnode,inode)=0.0_rp
             !               elmap(inode,jnode)=0.0_rp
             !            end do
             !            elmap(inode,inode)=facts
             !         end if
             !      end do
             !   end do
             !end if

          end if

       end do

    end if

  end subroutine nsi_element_schur

  subroutine nsi_assembly_projections(&
       list_elements,pgaus,pnode,ndofn,lnods_loc,elvel,elpre,rmom1,rmom2,gprhs,&
       gpgpr,gpsha,gpvol,gpden,gpadv,gpcar,gpsp1,gpsp2,gpst1)
    !------------------------------------------------------------------------
    !****f* Nastin/nsi_elmort
    ! NAME
    !    nsi_elmort
    ! DESCRIPTION
    !    This subroutine assemble the projected residuals
    ! USES
    ! USED BY
    !***
    !------------------------------------------------------------------------
    use def_kintyp, only       :  ip,rp
    use def_elmtyp, only       :  ELEXT
    use def_domain, only       :  ndime,mnode
!    use def_domain, only       :  lelch
    use def_nastin, only       :  kfl_stabi_nsi,vepr2_nsi,prpr2_nsi, kfl_regim_nsi
    use def_nastin, only       :  grpr2_nsi,kfl_rmom2_nsi
    implicit none
    integer(ip), intent(in)    :: list_elements(VECTOR_SIZE)
    integer(ip), intent(in)    :: pgaus
    integer(ip), intent(in)    :: pnode
    integer(ip), intent(in)    :: ndofn
    integer(ip), intent(in)    :: lnods_loc(VECTOR_SIZE,pnode)
    real(rp),    intent(in)    :: elvel(VECTOR_SIZE,ndime,pnode,*)
    real(rp),    intent(in)    :: elpre(VECTOR_SIZE,pnode,*)
    real(rp),    intent(inout) :: rmom1(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(in)    :: rmom2(VECTOR_SIZE,ndime,ndime,pnode,pgaus)
    real(rp),    intent(in)    :: gprhs(VECTOR_SIZE,ndofn,pgaus)
    real(rp),    intent(in)    :: gpgpr(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: gpsha(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(in)    :: gpvol(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpden(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpadv(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)
    real(rp),    intent(in)    :: gpsp1(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpsp2(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpst1(VECTOR_SIZE,pgaus)
    integer(ip)                :: idime,igaus,inode,jdime,ielem,ivect
    real(rp)                   :: gpvep(VECTOR_SIZE,ndime,pgaus)
    real(rp)                   :: gpprp(VECTOR_SIZE,pgaus)
    real(rp)                   :: gpgrp(VECTOR_SIZE,ndime,pgaus)
    real(rp)                   :: elvep(VECTOR_SIZE,ndime,pnode)
    real(rp)                   :: elprp(VECTOR_SIZE,pnode)
    real(rp)                   :: elgrp(VECTOR_SIZE,ndime,pnode)
    real(rp)                   :: taupr(VECTOR_SIZE,pgaus)
    real(rp)                   :: fact1(1:VECTOR_SIZE)
!    real(rp)                   :: dummr

    if( kfl_stabi_nsi == 0 ) then

       return

    else if( kfl_stabi_nsi == 1 ) then

       !-------------------------------------------------------------------
       !
       ! Full OSS
       !
       ! GPVEP: Projection of momentum residual
       ! = - tau1 * R(u)
       ! =   tau1 * [ du/dt + rho*(a.grad) u + grad(p) - div[2*mu*eps(u)] - f ]
       !
       !-------------------------------------------------------------------

       if( kfl_regim_nsi == 3 ) then
          taupr = 1.0_rp
       else
          taupr = gpst1
       end if

       do igaus = 1,pgaus
          do idime = 1,ndime
             gpvep(1:VECTOR_SIZE,idime,igaus) = gprhs(1:VECTOR_SIZE,idime,igaus) - gpgpr(1:VECTOR_SIZE,idime,igaus)
             do inode=1,pnode
                gpvep(1:VECTOR_SIZE,idime,igaus) = gpvep(1:VECTOR_SIZE,idime,igaus)&
                     - rmom1(1:VECTOR_SIZE,inode,igaus) * elvel(1:VECTOR_SIZE,idime,inode,1)
             end do
             gpvep(1:VECTOR_SIZE,idime,igaus) = - taupr(1:VECTOR_SIZE,igaus) * gpvep(1:VECTOR_SIZE,idime,igaus) * gpvol(1:VECTOR_SIZE,igaus)
          end do
       end do

       if( kfl_rmom2_nsi /= 0 ) then
          do igaus = 1,pgaus
             fact1 =  taupr(1:VECTOR_SIZE,igaus)* gpvol(1:VECTOR_SIZE,igaus)
             do idime = 1,ndime
                do inode = 1,pnode
                   do jdime = 1,ndime
                      gpvep(1:VECTOR_SIZE,idime,igaus) = gpvep(1:VECTOR_SIZE,idime,igaus) &
                           + rmom2(1:VECTOR_SIZE,idime,jdime,inode,igaus) * elvel(1:VECTOR_SIZE,jdime,inode,1) * fact1(1:VECTOR_SIZE)
                   end do
                end do
             end do
          end do
       end if

    else if( kfl_stabi_nsi == 2 ) then

       !-------------------------------------------------------------------
       !
       ! Split OSS
       !
       ! GPVEP: Projection of tau1' * rho * (a.grad) u
       ! GPGRP: Projection of tau1' * [ grad(p) - f ]
       !
       !-------------------------------------------------------------------

       do igaus = 1,pgaus

          do idime = 1,ndime
             gpvep(1:VECTOR_SIZE,idime,igaus) = 0.0_rp
             do inode = 1,pnode
                gpvep(1:VECTOR_SIZE,idime,igaus) = gpvep(1:VECTOR_SIZE,idime,igaus)&
                     + rmom1(1:VECTOR_SIZE,inode,igaus) * elvel(1:VECTOR_SIZE,idime,inode,1)
             end do
             fact1(1:VECTOR_SIZE)             = gpsp1(1:VECTOR_SIZE,igaus) * gpvol(1:VECTOR_SIZE,igaus)
             gpvep(1:VECTOR_SIZE,idime,igaus) = fact1(1:VECTOR_SIZE) * gpvep(1:VECTOR_SIZE,idime,igaus)
             gpgrp(1:VECTOR_SIZE,idime,igaus) = fact1(1:VECTOR_SIZE) * ( gpgpr(1:VECTOR_SIZE,idime,igaus) - gprhs(1:VECTOR_SIZE,idime,igaus) )
          end do

       end do

    end if

    !----------------------------------------------------------------------
    !
    ! Projection of tau2 * div(u)
    !
    !----------------------------------------------------------------------

    do igaus = 1,pgaus
       gpprp(1:VECTOR_SIZE,igaus) = 0.0_rp
       do idime = 1,ndime
          do inode=1,pnode
             gpprp(1:VECTOR_SIZE,igaus) = gpprp(1:VECTOR_SIZE,igaus) &
                  + gpcar(1:VECTOR_SIZE,idime,inode,igaus) * elvel(1:VECTOR_SIZE,idime,inode,1)
          end do
       end do
       gpprp(1:VECTOR_SIZE,igaus) = gpprp(1:VECTOR_SIZE,igaus) * gpsp2(1:VECTOR_SIZE,igaus) * gpvol(1:VECTOR_SIZE,igaus)
    end do

    !----------------------------------------------------------------------
    !
    ! Assemble RHS
    !
    !----------------------------------------------------------------------

    elvep = 0.0_rp
    elprp = 0.0_rp

    do igaus = 1,pgaus
       do inode = 1,pnode
          do idime = 1,ndime
             elvep(1:VECTOR_SIZE,idime,inode) = elvep(1:VECTOR_SIZE,idime,inode) + gpvep(1:VECTOR_SIZE,idime,igaus) * gpsha(1:VECTOR_SIZE,inode,igaus)
             elprp(1:VECTOR_SIZE,inode)       = elprp(1:VECTOR_SIZE,inode)       + gpprp(1:VECTOR_SIZE,igaus)   * gpsha(1:VECTOR_SIZE,inode,igaus)
          end do
       end do
    end do
    if( kfl_stabi_nsi == 2 ) then
       elgrp = 0.0_rp
       do igaus = 1,pgaus
          do inode = 1,pnode
             do idime = 1,ndime
                elgrp(1:VECTOR_SIZE,idime,inode) = elgrp(1:VECTOR_SIZE,idime,inode) + gpgrp(1:VECTOR_SIZE,idime,igaus) * gpsha(1:VECTOR_SIZE,inode,igaus)
             end do
          end do
       end do
    end if
    !
    ! Extension elements (Dodeme)
    !
    !if( lelch(ielem) == ELEXT ) then
    !   call elmext(&
    !        3_ip,ndime,pnode,dummr,dummr,dummr,dummr,dummr,&
    !        dummr,elvep,dummr)
    !   call elmext(&
    !        3_ip, 1_ip,pnode,dummr,dummr,dummr,dummr,dummr,&
    !        dummr,elprp,dummr)
    !   if( kfl_stabi_nsi == 2 ) then
    !      call elmext(&
    !           3_ip,ndime,pnode,dummr,dummr,dummr,dummr,dummr,&
    !           dummr,elgrp,dummr)
    !   end if
    !end if

    do ivect = 1,VECTOR_SIZE
       ielem = list_elements(ivect)
       if( ielem > 0 ) then
          call assrhs(ndime,pnode,lnods_loc(ivect,:),elvep(ivect,:,:),vepr2_nsi)
          call assrhs(1_ip, pnode,lnods_loc(ivect,:),elprp(ivect,:),prpr2_nsi)
          if( kfl_stabi_nsi == 2 ) then
             call assrhs(ndime,pnode,lnods_loc(ivect,:),elgrp(ivect,:,:),grpr2_nsi)
          end if
       end if
    end do

  end subroutine nsi_assembly_projections

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    21/09/2016
  !> @brief   This routine computes the element time step
  !> @details If properties are given as an input, then an average is
  !>          computed. If not, properties are computed at the center
  !>          of gravity
  !>
  !-----------------------------------------------------------------------


  subroutine nsi_element_time_step(&
       pnode,pgaus,list_elements,elcod,elvel,&
       chale,hleng,dtcri,gpden,gpvis,gpmut,gppor)
    !-----------------------------------------------------------------------
    !****f* Nastin/nsi_elmtss
    ! NAME
    !    nsi_elmtss
    ! DESCRIPTION
    !    This routine computes the element time step
    !    If properties are given as an input, then an average is
    !    computed. If not, properties are computed at the center
    !    of gravity
    ! USED BY
    !    nsi_updtss
    !***
    !-----------------------------------------------------------------------

    use def_domain,     only   :  ndime
    use def_nastin,     only   :  kfl_taust_nsi,staco_nsi,corio_nsi,&
         &                        kfl_advec_nsi
    use mod_ker_proper, only   :  ker_proper
    use mod_tauadr, only       :  tauadr

    implicit none
    integer(ip), intent(in)           :: pnode
    integer(ip), intent(in)           :: pgaus
    integer(ip), intent(in)           :: list_elements(VECTOR_SIZE)
    real(rp),    intent(in)           :: elcod(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(in)           :: elvel(VECTOR_SIZE,ndime,pnode,*)
    real(rp),    intent(in)           :: chale(VECTOR_SIZE,2)
    real(rp),    intent(in)           :: hleng(VECTOR_SIZE,ndime)
    real(rp),    intent(out)          :: dtcri(VECTOR_SIZE)

    real(rp),    intent(in), optional :: gpden(VECTOR_SIZE,pgaus)
    real(rp),    intent(in), optional :: gpvis(VECTOR_SIZE,pgaus)
    real(rp),    intent(in), optional :: gpmut(VECTOR_SIZE,pgaus)
    real(rp),    intent(in), optional :: gppor(VECTOR_SIZE,pgaus)

    integer(ip)                       :: idime,inode,igaus,ivect
    real(rp)                          :: rnode,rgaus

    real(rp)                          :: adv(VECTOR_SIZE)
    real(rp)                          :: dif(VECTOR_SIZE)
    real(rp)                          :: rea(VECTOR_SIZE)
    real(rp)                          :: elvno(VECTOR_SIZE)
    real(rp)                          :: gpvno_ave(VECTOR_SIZE)
    real(rp)                          :: gpvel_ave(VECTOR_SIZE, ndime)
    real(rp)                          :: gpden_ave(VECTOR_SIZE)
    real(rp)                          :: gpvis_ave(VECTOR_SIZE)
    real(rp)                          :: gppor_ave(VECTOR_SIZE)
!    real(rp)                          :: gpmut_ave(VECTOR_SIZE)
    !
    ! Initialization
    !
    elvno     = 0.0_rp
    rnode     = 1.0_rp / real(pnode,rp)
    rgaus     = 1.0_rp / real(pgaus,rp)
    gpvno_ave = 0.0_rp
    gpden_ave = 0.0_rp
    gpvis_ave = 0.0_rp
    gppor_ave = 0.0_rp
!    gpmut_ave = 0.0_rp
    !
    ! Average properties
    !
    if( present(gpden) ) then
       do igaus = 1,pgaus
          gpden_ave(1:VECTOR_SIZE) = gpden_ave(1:VECTOR_SIZE) + gpden(1:VECTOR_SIZE,igaus)
          gpvis_ave(1:VECTOR_SIZE) = gpvis_ave(1:VECTOR_SIZE) + gpvis(1:VECTOR_SIZE,igaus)
          gppor_ave(1:VECTOR_SIZE) = gppor_ave(1:VECTOR_SIZE) + gppor(1:VECTOR_SIZE,igaus)
!          gpmut_ave(1:VECTOR_SIZE) = gpmut_ave(1:VECTOR_SIZE) + gpmut(1:VECTOR_SIZE,igaus)
       end do
       gpden_ave(1:VECTOR_SIZE) = gpden_ave(1:VECTOR_SIZE) * rgaus
       gpvis_ave(1:VECTOR_SIZE) = gpvis_ave(1:VECTOR_SIZE) * rgaus
       gppor_ave(1:VECTOR_SIZE) = gppor_ave(1:VECTOR_SIZE) * rgaus
!       gpmut_ave(1:VECTOR_SIZE) = gpmut_ave(1:VECTOR_SIZE) * rgaus
    else
       call ker_proper('DENSI','COG  ',1_ip,list_elements,gpden_ave)
       call ker_proper('VISCO','COG  ',1_ip,list_elements,gpvis_ave)
       call ker_proper('POROS','COG  ',1_ip,list_elements,gppor_ave)
!       call ker_proper('TURBU','COG  ',1_ip,list_elements,gpmut_ave)
    end if
    !
    ! ELVNO_AVE = average velocity
    ! GPVNO_AVE = average velocity norm
    !
    if( kfl_advec_nsi /= 0 ) then
       gpvel_ave(1:VECTOR_SIZE, 1:ndime) = 0.0_rp
       do inode = 1,pnode
          gpvel_ave(1:VECTOR_SIZE, 1:ndime) =  gpvel_ave(1:VECTOR_SIZE, 1:ndime) + elvel(1:VECTOR_SIZE,1:ndime,inode,1)
       end do
       gpvel_ave = rnode*gpvel_ave
       gpvno_ave (1:VECTOR_SIZE) = 0.0_rp
       do idime =1, ndime
          gpvno_ave (1:VECTOR_SIZE) =  gpvno_ave (1:VECTOR_SIZE) +  gpvel_ave(1:VECTOR_SIZE, idime)*gpvel_ave(1:VECTOR_SIZE, idime)
       end do
       gpvno_ave (1:VECTOR_SIZE) = sqrt( gpvno_ave (1:VECTOR_SIZE)  )
    end if
    !
    ! DTCRI: Critical time step
    !
    adv(1:VECTOR_SIZE)   = gpden_ave(1:VECTOR_SIZE) * gpvno_ave(1:VECTOR_SIZE)                   ! Convective term ... rho*|u|
    dif(1:VECTOR_SIZE)   = gpvis_ave(1:VECTOR_SIZE) ! + gpmut_ave(1:VECTOR_SIZE)                   ! Viscous term ...... mu + mut ( gpvis accounts for gpvis+gptur)
    rea(1:VECTOR_SIZE)   = gpden_ave(1:VECTOR_SIZE) * corio_nsi + abs(gppor_ave(1:VECTOR_SIZE))  ! Coriolis .......... w + Porosity: sig

    do ivect = 1,VECTOR_SIZE
       call tauadr(&
            kfl_taust_nsi,staco_nsi,adv(ivect),dif(ivect),rea(ivect),&
            chale(ivect,1),chale(ivect,2),dtcri(ivect))       
    end do
    dtcri(1:VECTOR_SIZE) = gpden_ave(1:VECTOR_SIZE) * dtcri(1:VECTOR_SIZE)
    
   
  end subroutine nsi_element_time_step

  subroutine nsi_element_extension(&
       itask,list_elements,lelch_loc,pnode,elauu,elaup,elapu,elapp,&
       elmap,elrbu,elrbp)
    !-----------------------------------------------------------------------
    !****f* Nastin/nsi_elmext
    ! NAME
    !    nsi_elmext
    ! DESCRIPTION
    !    Modify element matrix when extension elements are used
    !    Only equation of the first node should be assembled as
    !    it corresponds to its extension test function
    !
    ! USES
    ! USED BY
    !    nsi_elmext
    !***
    !-----------------------------------------------------------------------

    integer(ip), intent(in)    :: itask
    integer(ip), intent(in)    :: pnode
    integer(ip), intent(in)    :: list_elements(VECTOR_SIZE)
    integer(ip), intent(in)    :: lelch_loc(VECTOR_SIZE)
    real(rp),    intent(out)   :: elauu(VECTOR_SIZE,pnode*ndime,pnode*ndime)
    real(rp),    intent(out)   :: elaup(VECTOR_SIZE,pnode*ndime,pnode)
    real(rp),    intent(out)   :: elapu(VECTOR_SIZE,pnode,pnode*ndime)
    real(rp),    intent(out)   :: elapp(VECTOR_SIZE,pnode,pnode)
    real(rp),    intent(out)   :: elmap(VECTOR_SIZE,pnode,pnode)
    real(rp),    intent(out)   :: elrbu(VECTOR_SIZE,pnode*ndime)
    real(rp),    intent(out)   :: elrbp(VECTOR_SIZE,pnode)
    integer(ip)                :: inode,jnode,idofn,jdofn
    integer(ip)                :: ielem,iwhat,ivect

    if( itask == 1 ) then
       !
       ! Auu, Aup, Apu, App, bu, bp
       !
       do ivect = 1,VECTOR_SIZE
          ielem = list_elements(ivect)
          if( ielem > 0 ) then
             if( lelch_loc(ivect) == ELEXT ) then
                do idofn = ndime+1,ndime*pnode
                   do jdofn = 1,ndime*pnode
                      elauu(ivect,idofn,jdofn) = 0.0_rp
                   end do
                   do jnode = 1,pnode
                      elaup(ivect,idofn,jnode) = 0.0_rp
                   end do
                   elrbu(ivect,idofn) = 0.0_rp
                end do
                do inode = 2,pnode
                   do jdofn = 1,ndime*pnode
                      elapu(ivect,inode,jdofn) = 0.0_rp
                   end do
                   do jnode = 1,pnode
                      elapp(ivect,inode,jnode) = 0.0_rp
                   end do
                   elrbp(ivect,inode) = 0.0_rp
                end do
             end if
          end if
       end do

    else if( itask == 2 ) then
       !
       ! Q
       !
       iwhat = 1

       if( iwhat == 1 ) then
          !
          ! Full matrix
          !
          return

       else if( iwhat == 2 ) then
          !
          ! All zero
          !
          do ivect = 1,VECTOR_SIZE
             ielem = list_elements(ivect)
             if( ielem > 0 ) then
                if( lelch_loc(ivect) == ELEXT ) then
                   do inode = 1,pnode
                      do jnode = 1,pnode
                         elmap(ivect,inode,jnode) = 0.0_rp
                      end do
                   end do
                end if
             end if
          end do

       else  if( iwhat == 3 ) then
          !
          ! Diagonal
          !
          do ivect = 1,VECTOR_SIZE
             ielem = list_elements(ivect)
             if( ielem > 0 ) then
                if( lelch_loc(ivect) == ELEXT ) then
                   do inode = 1,pnode
                      do jnode = 1,pnode
                         if( inode /= jnode ) elmap(ivect,inode,jnode) = 0.0_rp
                      end do
                   end do
                end if
             end if
          end do

       else if( iwhat == 4 ) then
          !
          ! Symmetrized
          !
          do ivect = 1,VECTOR_SIZE
             ielem = list_elements(ivect)
             if( ielem > 0 ) then
                if( lelch_loc(ivect) == ELEXT ) then
                   do inode = 2,pnode
                      do jnode = 2,pnode
                         elmap(ivect,inode,jnode) = 0.0_rp
                      end do
                   end do
                   do inode = 1,pnode
                      do jnode = 1,pnode
                         if( abs(elmap(ivect,inode,jnode)-elmap(ivect,jnode,inode)) > 1.0e-8_rp ) print*,'merde'
                      end do
                   end do
                end if
             end if
          end do

       end if

    end if

  end subroutine nsi_element_extension

  subroutine nsi_element_manufactured_solution(&
       pgaus,pnode,gpsha,elcod,gpden,gpvis,gppor,gpgvi,&
       cutim,baloc,gprhs,gprhc,gprh2)
    !-----------------------------------------------------------------------
    !****f* Nastin/nsi_elmexa
    ! NAME
    !    nsi_elmexa
    ! DESCRIPTION
    !    Compute RHS for exact solution
    ! USES
    !    nsi_exacso
    ! USED BY
    !***
    !-----------------------------------------------------------------------

    integer(ip), intent(in)  :: pgaus
    integer(ip), intent(in)  :: pnode
    real(rp),    intent(in)  :: gpsha(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(in)  :: elcod(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(in)  :: gpden(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)  :: gpvis(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)  :: gppor(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)  :: gpgvi(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)  :: cutim
    real(rp),    intent(in)  :: baloc(VECTOR_SIZE,ndime)
    real(rp),    intent(out) :: gprhs(VECTOR_SIZE,ndime,*)
    real(rp),    intent(out) :: gprhc(VECTOR_SIZE,*)
    real(rp),    intent(out) :: gprh2(VECTOR_SIZE,*)
    integer(ip)              :: igaus,inode,idime,ivect
    real(rp)                 :: dummr
    real(rp)                 :: gpcod(VECTOR_SIZE,3)

    if( kfl_exacs_nsi /= 0 ) then

       if( pgaus == -1 ) then
          !
          ! Contribution on a a particular boundary Gauss point
          !
          igaus = 1
          gpcod = 0.0_rp
          do inode = 1,pnode
             do idime = 1,ndime
                gpcod(1:VECTOR_SIZE,idime) = gpcod(1:VECTOR_SIZE,idime) + elcod(1:VECTOR_SIZE,idime,inode) * gpsha(1:VECTOR_SIZE,inode,1)
             end do
          end do
          do ivect = 1,VECTOR_SIZE
             call nsi_exacso(&
                  3_ip,gpcod(ivect,:),gpden(ivect,igaus),gpvis(ivect,igaus),&
                  gppor(ivect,igaus),gpgvi(ivect,:,igaus),&
                  dummr,dummr,dummr,dummr,baloc(ivect,:),&
                  gprhs(ivect,:,igaus),gprhc(ivect,igaus),gprh2(ivect,igaus))
          end do

       else
          !
          ! Contribution on an element
          !
          do igaus = 1,pgaus
             gpcod = 0.0_rp
             do inode = 1,pnode
                do idime = 1,ndime
                   gpcod(1:VECTOR_SIZE,idime) = gpcod(1:VECTOR_SIZE,idime) + elcod(1:VECTOR_SIZE,idime,inode) * gpsha(1:VECTOR_SIZE,inode,igaus)
                end do
             end do
             do ivect = 1,VECTOR_SIZE
                call nsi_exacso(&
                     2_ip,gpcod(ivect,:),gpden(ivect,igaus),gpvis(ivect,igaus),&
                     gppor(ivect,igaus),gpgvi(ivect,:,igaus),&
                     dummr,dummr,dummr,dummr,baloc(ivect,:),&
                     gprhs(ivect,:,igaus),gprhc(ivect,igaus),gprh2(ivect,igaus))
             end do
          end do

       end if

    end if

  end subroutine nsi_element_manufactured_solution

  !-----------------------------------------------------------------------
  !> @addtogroup Nastin
  !> @{
  !> @file    nsi_elmexf.f90
  !> @author  Guillaume Houzeaux
  !> @brief   Add a force term to the momentum equations
  !> @details
  !> @}
  !-----------------------------------------------------------------------

  subroutine nsi_element_external_force(&
       pnode,pgaus,list_elements,lnods_loc,gpsha,lmate_loc,gpden,gprhs,gpvel, elcod, gppor, gpapo)
    
    use def_nastin,          only : lforc_material_nsi
    use def_nastin,          only : xforc_material_nsi, ntabr_nsi, fvela_nsi, kfl_force_nsi
    use def_nastin,          only : kfl_vegeo_time_nsi
    use mod_ker_tendencies,  only : kfl_tendencies_ker, get_tendencies_u, get_tendencies_uwrf
    use def_master,          only : cutim,FUNCTION_DISCRETE
    use def_domain,          only : walld
    use def_kermod,          only : kfl_dampi_ker !  anipo_ker, poros_ker, 
    use mod_ker_functions         ,  only : ker_functions
    use mod_ker_discrete_function,   only : ker_discrete_function
    integer(ip), intent(in)  :: pgaus, pnode
    integer(ip), intent(in)  :: list_elements(VECTOR_SIZE)
    integer(ip), intent(in)  :: lmate_loc(VECTOR_SIZE)
    real(rp),    intent(in)  :: gpden(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)  :: gpvel(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(out) :: gprhs(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)  :: elcod(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(in)  :: gpsha(VECTOR_SIZE,pnode,pgaus)
    integer(ip), intent(in)  :: lnods_loc(VECTOR_SIZE, pnode)
    real(rp),    intent(in)  :: gppor(VECTOR_SIZE,pgaus), gpapo(VECTOR_SIZE, ndime, ndime, pgaus)
    integer(ip)              :: idime,igaus,ivect,pmate,ielem, inode
    real(rp)                 :: n_disk(3),Ct,uref2,uref3,V_disk
    real(rp)                 :: Force_normal, d_disk, veinf, Cp, fcori
    real(rp)                 :: gphei, elhei(pnode)
    real(rp)                 :: gpten_geo(2), gpten_adv(2), gpvel_ref(2), ugeos(2)

    !
    ! ADD TENDENCIES to the RHS of Momentum EQ
    !
    if (kfl_tendencies_ker) then
       fcori = 2.0_rp*sqrt(fvela_nsi(1)*fvela_nsi(1) + &
            fvela_nsi(2)*fvela_nsi(2) + &
            fvela_nsi(3)*fvela_nsi(3))
       loop_ivect1:    do ivect = 1,VECTOR_SIZE
          ielem = list_elements(ivect)
          if( ielem > 0 ) then
             pmate = lmate_loc(ivect)
!           interpolate tendencies to present height (over terrain?)
             elhei(1:pnode) = walld(lnods_loc(ivect, 1:pnode))  !elcod(ivect, 3, 1:pnode)  ! by the moment height is z coordinate
             do igaus = 1,pgaus
                gphei =0.0_rp
                do inode =1, pnode
                   gphei = gphei + gpsha(ivect,inode, igaus)*elhei(inode)
                end do
                
                call get_tendencies_u(gphei,cutim, gpten_geo, gpten_adv) ! get geostrophic and advection forces
                gprhs(ivect,1, igaus) =  gprhs(ivect,1, igaus) + gpden(ivect,igaus)*fcori*(-gpten_geo(2) + gpten_adv(1))
                gprhs(ivect,2, igaus) =  gprhs(ivect,2, igaus) + gpden(ivect,igaus)*fcori*( gpten_geo(1) + gpten_adv(2))
                !
                ! RHS DAMPING  damp*vel (not possible to differentiate between damping and forest)
                ! 

                if (kfl_dampi_ker(pmate).gt.0) then !.and.kfl_anipo_nsi == 0) then
                   call get_tendencies_uwrf(gphei, cutim, gpvel_ref) 
                   if (kfl_anipo_nsi==0) then
                      gprhs(ivect,1:2, igaus) =  gprhs(ivect,1:2, igaus) + gppor(ivect, igaus)*gpvel_ref(1:2)
                   else
                      gprhs(ivect,1:2, igaus) =  gprhs(ivect,1:2, igaus) + gpapo(ivect, 1,1, igaus)*gpvel_ref(1:2)
                   end if

                end if
             end do  
            
          end if
       end do loop_ivect1
    else if (kfl_vegeo_time_nsi.gt.0) then  !Time dependent pressure gradient (acts like a time dependent and uniform tendency)
       fcori = 2.0_rp*sqrt(fvela_nsi(1)*fvela_nsi(1) + &
            fvela_nsi(2)*fvela_nsi(2) + &
            fvela_nsi(3)*fvela_nsi(3))
       loop_ivect2:    do ivect = 1,VECTOR_SIZE
          ielem = list_elements(ivect)
          if( ielem > 0 ) then
             pmate = lmate_loc(ivect)
!            Time interpolation of geostrophic velocity      
             call ker_discrete_function(kfl_vegeo_time_nsi,cutim,ugeos)
             gprhs(ivect,1, 1:pgaus) =  gprhs(ivect,1, 1:pgaus) -  gpden(ivect,1:pgaus)*fcori*ugeos(2)
             gprhs(ivect,2, 1:pgaus) =  gprhs(ivect,2, 1:pgaus) +  gpden(ivect,1:pgaus)*fcori*ugeos(1)
             ! add RHS damping
             if (kfl_dampi_ker(pmate).gt.0) then
                if (kfl_anipo_nsi == 0) then              
                   gprhs(ivect,1, 1:pgaus) =  gprhs(ivect,1, 1:pgaus) + gppor(ivect, 1:pgaus)*ugeos(1)
                   gprhs(ivect,2, 1:pgaus) =  gprhs(ivect,2, 1:pgaus) + gppor(ivect, 1:pgaus)*ugeos(2)
                else
                   gprhs(ivect,1, 1:pgaus) =  gprhs(ivect,1, 1:pgaus) + gpapo(ivect, 1,1,1:pgaus)*ugeos(1)                   
                   gprhs(ivect,2, 1:pgaus) =  gprhs(ivect,2, 1:pgaus) + gpapo(ivect, 1,1,1:pgaus)*ugeos(2)                   
                end if
             end if
             
          end if
       end do loop_ivect2
    end if

    if (kfl_force_nsi==1) then ! actuator disk force
       loop_ivect:    do ivect = 1,VECTOR_SIZE
          ielem = list_elements(ivect)
          if( ielem > 0 ) then
             pmate = lmate_loc(ivect)


             if( lforc_material_nsi(pmate) == 1 ) then
                !
                ! Wake model
                !
                n_disk(1) = xforc_material_nsi( 3,pmate)
                n_disk(2) = xforc_material_nsi( 4,pmate)
                n_disk(3) = xforc_material_nsi( 5,pmate)
                V_disk    = xforc_material_nsi(15,pmate)
                d_disk    = xforc_material_nsi( 1,pmate)
                Veinf     = xforc_material_nsi(18, pmate) ! infinite velocity
                !
                ! traction coefficient
                !
                Ct        =  xforc_material_nsi(16,pmate)
                Cp        =  xforc_material_nsi(17,pmate)
                uref2     =  veinf*veinf
                uref3     =  veinf*veinf*veinf
                !
                !  modification of Ct to use local velocity
                !  assuming theoretical relationship Ct=4a*(1-a)
                !
                !             Ct         = 4.0_rp*Ct/(2.0_rp+2.0_rp*sqrt(1.0_rp-Ct)-Ct)
                if (ntabr_nsi(pmate)==0) then !Uniform loaded disc, non normal and tangential forces distibution

                   do igaus = 1,pgaus

                      Force_normal = 0.5_rp * gpden(ivect, igaus) * Ct * uref2/d_disk
                      ! assembly to rhs
                      gprhs(ivect, 1:ndime,igaus) = gprhs(ivect,1:ndime,igaus) + n_disk(1:ndime) * Force_normal
                   end do
                else
                   call runend('nsi_element_external_force: force with rotational force not coded in the new Alya, Im lazy')
                end if

             else if ( lforc_material_nsi(pmate) == 2 ) then
                !
                ! Constant model
                !
                do igaus = 1,pgaus
                   do idime= 1,ndime
                      gprhs(ivect,idime,igaus) =  gprhs(ivect,idime,igaus) +  xforc_material_nsi(idime,pmate)
                   end do
                end do

             end if


          end if ! ielem > 0

       end do loop_ivect
    end if

  end subroutine nsi_element_external_force

  !------------------------------------------------------------------------
  ! NOTES
  !
  ! Shock capturing for the advection-diffusion-reaction equation
  !
  !     du
  ! rho -- + rho*a.grad(u) - div[k*grad(u)] + s*u = f
  !     dt
  !
  ! k    = Diffusion                  [M/(L*T)]
  ! R    = Residual of the equation   [M*U/(L^3*T)]
  ! C    = Shock capturing constant
  ! tau  = Stabilization parameter:
  !        Its units are h/(rho*a)=   [L^3*T/M]
  !        so that rho*tau is in [T]
  ! kiso = Isotropic SC diffusion     [M/(L*T)]
  ! k'   = Anisotropic SC diffusion   [M/(L*T)]
  !
  !        1    |R|    h
  ! Pe   = - --------- -
  !        2 |grad(u)| k
  !
  !            +          2k   +           R
  ! ac   = max | 0 , C - ----- | , a*= ----------- grad(u), therefore
  !            +         |a*|h +       |grad(u)|^2
  !
  !            +           1   +
  ! ac   = max | 0 , C - ----- |
  !            +          Pe   +
  !        1          |R|
  ! kiso = - ac*h  --------- , k'=(rho*tau)*rho*a^2
  !        2       |grad(u)|
  !
  ! Isotropic:     kiso*grad(u).grad(v)
  !                                         a x a
  ! Anisotropic:   (<kiso-k'>-kiso)*grad(u).-----.grad(v)
  !                                          a^2
  !***
  !------------------------------------------------------------------------

  subroutine nsi_element_shock_capturing(&
       pnode,pgaus,ptopo,pevat,ndofn,gpden,gpvel,gprhs,gpsp1,&
       gpvol,elvel,gpcar,chale,rmomu,rmom2,elauu)

    use def_kintyp, only       :  ip,rp
    use def_domain, only       :  ndime,mnode
    use def_nastin, only       :  kfl_shock_nsi,shock_nsi,kfl_rmom2_nsi
    implicit none

    integer(ip), intent(in)    :: pnode
    integer(ip), intent(in)    :: pgaus
    integer(ip), intent(in)    :: ptopo
    integer(ip), intent(in)    :: pevat
    integer(ip), intent(in)    :: ndofn
    real(rp),    intent(in)    :: gpvel(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)
    real(rp),    intent(in)    :: elvel(VECTOR_SIZE,ndime,pnode,1)
    real(rp),    intent(in)    :: gpden(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gprhs(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: gpvol(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: chale(VECTOR_SIZE,2)
    real(rp),    intent(in)    :: gpsp1(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: rmomu(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(in)    :: rmom2(VECTOR_SIZE,ndime,ndime,pnode,pgaus)
    real(rp),    intent(inout) :: elauu(VECTOR_SIZE,pevat,pevat)

    real(rp)                   :: scdif,sldif,vepar,fact1,addit
    real(rp)                   :: factt,facta
    integer(ip)                :: idime,inode,jnode,kdime,ldime,igaus
    integer(ip)                :: idofn,jdofn,jdime,ivect

    real(rp)                   :: gpres(VECTOR_SIZE,pgaus)
    real(rp)                   :: gpve2(VECTOR_SIZE)
    real(rp)                   :: gpgrd(VECTOR_SIZE,ndime,ndime)
    real(rp)                   :: grnor(VECTOR_SIZE,ndime)
    !
    ! Factor according to element topology
    !
    if( ptopo == 0 ) then
       factt = 1.5_rp             ! Quadrilateral/Hexahedra
    else if( ptopo == 1 ) then
       factt = 0.7_rp             ! Triangles/Tetrahedra
    else
       factt = 1.0_rp             ! Other elements
    end if
    !
    ! Isotropic/anisotropic shock capturing
    !
    if( kfl_shock_nsi == 1 ) then
       facta = 0.0_rp             ! Isotropic SC
    else
       facta = 1.0_rp             ! Anisotropic SC
    end if
    !
    ! Coarse grid residual |R(u)|
    ! |R(u)| = |f - rho*u/(theta*dt)-rho*a.grad(u)+div[k*grad(u)]-s*u|
    ! f = Q + rho*u^n/(theta*dt)
    !
    do igaus = 1,pgaus
       gpres(1:VECTOR_SIZE,igaus) = 0.0_rp
       do inode = 1,pnode
          do idime = 1,ndime
             gpres(1:VECTOR_SIZE,igaus) = gpres(1:VECTOR_SIZE,igaus)&
                  - rmomu(1:VECTOR_SIZE,inode,igaus) * elvel(1:VECTOR_SIZE,idime,inode,1)
          end do
       end do
    end do
    if( kfl_rmom2_nsi /= 0 ) then
       do igaus = 1,pgaus
          do idime = 1,ndime
             do inode = 1,pnode
                do jdime = 1,ndime
                   gpres(1:VECTOR_SIZE,igaus) = gpres(1:VECTOR_SIZE,igaus)&
                        - rmom2(1:VECTOR_SIZE,idime,jdime,inode,igaus) * elvel(1:VECTOR_SIZE,jdime,inode,1)
                end do
             end do
          end do
       end do
    end if
    gpres(1:VECTOR_SIZE,1:pgaus) = abs(gpres(1:VECTOR_SIZE,1:pgaus))
    !
    ! Compute the contribution to the element stiffness matrix
    !
    do igaus = 1,pgaus
       !
       ! Square velocity norm
       !
       gpve2(1:VECTOR_SIZE) = 0.0_rp
       do idime = 1,ndime
          gpve2(1:VECTOR_SIZE) = gpve2(1:VECTOR_SIZE) + gpvel(VECTOR_SIZE,idime,igaus) * gpvel(VECTOR_SIZE,idime,igaus)
       end do
       !
       ! Velocity gradient
       !
       gpgrd(VECTOR_SIZE,:,:) = 0.0_rp
       do inode = 1,pnode
          do jdime = 1,ndime
             do idime=1,ndime
                gpgrd(1:VECTOR_SIZE,idime,jdime) = gpgrd(1:VECTOR_SIZE,idime,jdime)&
                     + gpcar(1:VECTOR_SIZE,idime,inode,igaus) * elvel(1:VECTOR_SIZE,jdime,inode,1)
             end do
          end do
       end do
       grnor(VECTOR_SIZE,:) = 0.0_rp
       do idime = 1,ndime
          do jdime = 1,ndime
             grnor(1:VECTOR_SIZE,idime) = grnor(1:VECTOR_SIZE,idime) + gpgrd(1:VECTOR_SIZE,jdime,idime) * gpgrd(1:VECTOR_SIZE,jdime,idime)
          end do
       end do
       grnor = sqrt(grnor)

       do ivect = 1,VECTOR_SIZE

          if( grnor(ivect,1) > 1.0e-10_rp ) then
             vepar = gpres(ivect,igaus) / grnor(ivect,1)
             scdif = factt * 0.5_rp * shock_nsi * chale(ivect,1) * vepar
             if( scdif > 0.0_rp ) then
                !
                ! Compute diffusion introduced along the streamlines
                !
                sldif = max(0.0_rp,scdif-facta*gpsp1(ivect,igaus)*gpden(ivect,igaus)*gpve2(ivect))
                if( gpve2(ivect) > 1.0e-10_rp ) then
                   fact1 = (sldif-scdif) / gpve2(ivect)
                else
                   fact1 = 0.0_rp
                end if
                do inode = 1,pnode
                   idofn = inode * ndofn
                   do jnode = 1,pnode
                      jdofn = jnode * ndofn
                      addit = 0.0_rp
                      do kdime = 1,ndime
                         addit = addit + scdif * gpcar(ivect,kdime,inode,igaus)&   ! kiso*grad(Ni).grad(Nj)
                              &                * gpcar(ivect,kdime,jnode,igaus)
                         do ldime = 1,ndime                                        ! (<kiso-k'>-kiso)/u^2*
                            addit = addit + fact1 &                                ! grad(Ni).(u x u).grad(Nj)
                                 * gpcar(ivect,kdime,inode,igaus)&
                                 * gpvel(ivect,kdime,igaus)&
                                 * gpvel(ivect,ldime,igaus)&
                                 * gpcar(ivect,ldime,jnode,igaus)
                         end do
                      end do
                      elauu(ivect,idofn,jdofn) = elauu(ivect,idofn,jdofn) &
                           + addit * gpvol(ivect,igaus)
                   end do
                end do
             end if
          end if

       end do

    end do

  end subroutine nsi_element_shock_capturing

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   Lumped matrices
  !> @details Lumped matrices
  !>
  !>          DT_RHO_NSI = Projection of dt / rho
  !>          TAU_NSI    = Projection of tau
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_rhodt_rhotau_nu_vector(&
       pnode,pgaus,porde,gpsha,gpvol,gpden,gpvis,gppor,gpadv,&
       chale,dtinv_loc,densi,eldtrho,eltau,elmurho)

    use def_nastin, only : corio_nsi
    use def_nastin, only : kfl_taust_nsi
    use def_nastin, only : kfl_regim_nsi,nbdfp_nsi
    use def_nastin, only : staco_nsi
    use mod_tauadr, only : tauadr

    integer(ip), intent(in)    :: pnode
    integer(ip), intent(in)    :: pgaus
    integer(ip), intent(in)    :: porde
    real(rp),    intent(in)    :: gpsha(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(in)    :: gpvol(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpden(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpvis(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gppor(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpadv(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: chale(VECTOR_SIZE,2)
    real(rp),    intent(in)    :: dtinv_loc(VECTOR_SIZE)
    real(rp),    intent(in)    :: densi(VECTOR_SIZE,pgaus,nbdfp_nsi)
    real(rp),    intent(out)   :: eldtrho(VECTOR_SIZE,pnode)
    real(rp),    intent(out)   :: elmurho(VECTOR_SIZE,pnode)
    real(rp),    intent(out)   :: eltau(VECTOR_SIZE,pnode)

    integer(ip)                :: inode,jnode,igaus,idime,ivect
    real(rp)                   :: gpsp1(VECTOR_SIZE,pgaus)
    real(rp)                   :: fact2(VECTOR_SIZE)
    real(rp)                   :: fact3(VECTOR_SIZE)
    real(rp)                   :: fact4(VECTOR_SIZE)
    real(rp)                   :: gpvno(VECTOR_SIZE)
    real(rp)                   :: adv(VECTOR_SIZE)
    real(rp)                   :: dif(VECTOR_SIZE)
    real(rp)                   :: rea(VECTOR_SIZE)
    real(rp)                   :: h2(VECTOR_SIZE)
!    real(rp)                   :: d_dtrho(VECTOR_SIZE)
!    real(rp)                   :: d_murho(VECTOR_SIZE)
!    real(rp)                   :: d_tau(VECTOR_SIZE)
!    real(rp)                   :: T_dtrho(VECTOR_SIZE)
!    real(rp)                   :: T_murho(VECTOR_SIZE)
!    real(rp)                   :: T_tau(VECTOR_SIZE)

    real(rp)                   :: elmat1(VECTOR_SIZE,pnode,pnode)
    real(rp)                   :: trace1(VECTOR_SIZE), elmass1(VECTOR_SIZE)
    real(rp)                   :: elmat2(VECTOR_SIZE,pnode,pnode)
    real(rp)                   :: trace2(VECTOR_SIZE), elmass2(VECTOR_SIZE)

    !
    ! tau1 = 1 / ( 4*mu/h^2 + 2*rho*uc/h + rho*|w| + sig )
    !
    do igaus = 1,pgaus

       gpvno(1:VECTOR_SIZE) = 0.0_rp
       do idime = 1,ndime
          gpvno(1:VECTOR_SIZE) = gpvno(1:VECTOR_SIZE) + gpadv(1:VECTOR_SIZE,idime,igaus) * gpadv(1:VECTOR_SIZE,idime,igaus)
       end do
       gpvno(1:VECTOR_SIZE) = sqrt(gpvno(1:VECTOR_SIZE))
       adv(1:VECTOR_SIZE)   = gpden(1:VECTOR_SIZE,igaus) * gpvno(1:VECTOR_SIZE)                          ! Convective term: rho*|u+u'|
       dif(1:VECTOR_SIZE)   = gpvis(1:VECTOR_SIZE,igaus)                                                 ! Viscous term:    mu
       rea(1:VECTOR_SIZE)   = gpden(1:VECTOR_SIZE,igaus) * corio_nsi + abs(gppor(1:VECTOR_SIZE,igaus))   ! Coriolis: w + Porosity: sig
       h2(1:VECTOR_SIZE)    = chale(1:VECTOR_SIZE,2) * chale(1:VECTOR_SIZE,2)

       do ivect = 1,VECTOR_SIZE
          call tauadr(&
               kfl_taust_nsi,staco_nsi,adv(ivect),dif(ivect),rea(ivect),&
               chale(ivect,1),chale(ivect,2),gpsp1(ivect,igaus), gpden(ivect,igaus)*dtinv_loc(ivect))
       end do

    end do
    !
    ! Element assembly
    !
    eldtrho(1:VECTOR_SIZE,:) = 0.0_rp
    eltau(1:VECTOR_SIZE,:)   = 0.0_rp
    elmurho(1:VECTOR_SIZE,:) = 0.0_rp
    if( porde == 1 ) then
       if(kfl_regim_nsi == 3)  then
          do igaus = 1,pgaus
             fact2(1:VECTOR_SIZE) = gpvol(1:VECTOR_SIZE,igaus) / ( gpden(1:VECTOR_SIZE,igaus) * dtinv_loc(1:VECTOR_SIZE) )
             fact3(1:VECTOR_SIZE) = gpvol(1:VECTOR_SIZE,igaus) * gpsp1(1:VECTOR_SIZE,igaus)
             fact4(1:VECTOR_SIZE) =  gpvol(1:VECTOR_SIZE,igaus) * densi(DEF_VECT, igaus, 1)
             do inode = 1,pnode
                eldtrho(1:VECTOR_SIZE,inode) = eldtrho(1:VECTOR_SIZE,inode) + gpsha(1:VECTOR_SIZE,inode,igaus) * fact2(1:VECTOR_SIZE)
                eltau(1:VECTOR_SIZE,inode)   = eltau(1:VECTOR_SIZE,inode)   + gpsha(1:VECTOR_SIZE,inode,igaus) * fact3(1:VECTOR_SIZE)
                elmurho(1:VECTOR_SIZE,inode) = elmurho(1:VECTOR_SIZE,inode) + gpsha(1:VECTOR_SIZE,inode,igaus) * fact4(1:VECTOR_SIZE)
             end do
          end do
       else  
          do igaus = 1,pgaus
             fact2(1:VECTOR_SIZE) = gpvol(1:VECTOR_SIZE,igaus) / ( gpden(1:VECTOR_SIZE,igaus) * dtinv_loc(1:VECTOR_SIZE) )
             fact3(1:VECTOR_SIZE) = gpvol(1:VECTOR_SIZE,igaus) * gpsp1(1:VECTOR_SIZE,igaus)
             !fact4(1:VECTOR_SIZE) = gpvol(1:VECTOR_SIZE,igaus) * gpvis(1:VECTOR_SIZE,igaus) / ( gpden(1:VECTOR_SIZE,igaus)  )
             fact4(1:VECTOR_SIZE) =  gpvol(1:VECTOR_SIZE,igaus) * ( gpden(1:VECTOR_SIZE,igaus)  )
             do inode = 1,pnode
                eldtrho(1:VECTOR_SIZE,inode) = eldtrho(1:VECTOR_SIZE,inode) + gpsha(1:VECTOR_SIZE,inode,igaus) * fact2(1:VECTOR_SIZE)
                eltau(1:VECTOR_SIZE,inode)   = eltau(1:VECTOR_SIZE,inode)   + gpsha(1:VECTOR_SIZE,inode,igaus) * fact3(1:VECTOR_SIZE)
                elmurho(1:VECTOR_SIZE,inode) = elmurho(1:VECTOR_SIZE,inode) + gpsha(1:VECTOR_SIZE,inode,igaus) * fact4(1:VECTOR_SIZE)
             end do
          end do
       end if
    else
      ! T_dtrho(1:VECTOR_SIZE) = 0.0_rp
      ! d_dtrho(1:VECTOR_SIZE) = 0.0_rp
      ! T_tau(1:VECTOR_SIZE)   = 0.0_rp
      ! d_tau(1:VECTOR_SIZE)   = 0.0_rp
      ! T_murho(1:VECTOR_SIZE) = 0.0_rp
      ! d_murho(1:VECTOR_SIZE) = 0.0_rp
      ! do igaus = 1,pgaus
      !    fact2(1:VECTOR_SIZE) = gpvol(1:VECTOR_SIZE,igaus) / ( gpden(1:VECTOR_SIZE,igaus) * dtinv_loc(1:VECTOR_SIZE) )
      !    fact3(1:VECTOR_SIZE) = gpvol(1:VECTOR_SIZE,igaus) * gpsp1(1:VECTOR_SIZE,igaus)
      !    !fact4(1:VECTOR_SIZE) = gpvol(1:VECTOR_SIZE,igaus) * gpvis(1:VECTOR_SIZE,igaus) / ( gpden(1:VECTOR_SIZE,igaus)  )
      !    fact4(1:VECTOR_SIZE) = gpvol(1:VECTOR_SIZE,igaus) * ( gpden(1:VECTOR_SIZE,igaus)  )
      !    T_dtrho(1:VECTOR_SIZE) = T_dtrho(1:VECTOR_SIZE) + fact2(1:VECTOR_SIZE)
      !    T_tau(1:VECTOR_SIZE)   = T_tau(1:VECTOR_SIZE)   + fact3(1:VECTOR_SIZE)
      !    T_murho(1:VECTOR_SIZE) = T_murho(1:VECTOR_SIZE) + fact4(1:VECTOR_SIZE)
      !    do inode = 1,pnode
      !       eldtrho(1:VECTOR_SIZE,inode) = eldtrho(1:VECTOR_SIZE,inode) + gpsha(1:VECTOR_SIZE,inode,igaus)**2 * fact2(1:VECTOR_SIZE)
      !       eltau(1:VECTOR_SIZE,inode)   = eltau(1:VECTOR_SIZE,inode)   + gpsha(1:VECTOR_SIZE,inode,igaus)**2 * fact3(1:VECTOR_SIZE)
      !       elmurho(1:VECTOR_SIZE,inode) = elmurho(1:VECTOR_SIZE,inode) + gpsha(1:VECTOR_SIZE,inode,igaus)**2 * fact4(1:VECTOR_SIZE)
      !    end do
      ! end do
      ! do inode = 1,pnode
      !    d_dtrho(1:VECTOR_SIZE) = d_dtrho(1:VECTOR_SIZE) + eldtrho(1:VECTOR_SIZE,inode)
      !    d_tau(1:VECTOR_SIZE)   = d_tau(1:VECTOR_SIZE)   + eltau(1:VECTOR_SIZE,inode)
      !    d_murho(1:VECTOR_SIZE) = d_murho(1:VECTOR_SIZE) + elmurho(1:VECTOR_SIZE,inode)
      ! end do

      ! do inode = 1,pnode
      !    eldtrho(1:VECTOR_SIZE,inode) = eldtrho(1:VECTOR_SIZE,inode) * T_dtrho(1:VECTOR_SIZE) / d_dtrho(1:VECTOR_SIZE)
      !    eltau(1:VECTOR_SIZE,inode)   = eltau(1:VECTOR_SIZE,inode)   * T_tau(1:VECTOR_SIZE) / d_tau(1:VECTOR_SIZE)
      !    elmurho(1:VECTOR_SIZE,inode) = elmurho(1:VECTOR_SIZE,inode) * T_murho(1:VECTOR_SIZE) / d_murho(1:VECTOR_SIZE)
      ! end do

      do inode=1,pnode
         do jnode=1,pnode
            elmat1(1:VECTOR_SIZE,inode,jnode)=0.0_rp
            elmat2(1:VECTOR_SIZE,inode,jnode)=0.0_rp
         end do
      end do

      do igaus=1,pgaus
         do inode=1,pnode
            fact2(1:VECTOR_SIZE) = gpvol(1:VECTOR_SIZE,igaus)*gpsha(1:VECTOR_SIZE,inode,igaus)/(gpden(1:VECTOR_SIZE,igaus)*dtinv_loc(1:VECTOR_SIZE))
            fact3(1:VECTOR_SIZE) = gpvol(1:VECTOR_SIZE,igaus)*gpsha(1:VECTOR_SIZE,inode,igaus)*(gpden(1:VECTOR_SIZE,igaus))
            
            
            do jnode=1,pnode
               elmat1(1:VECTOR_SIZE,inode,jnode)=elmat1(1:VECTOR_SIZE,inode,jnode) +fact2(1:VECTOR_SIZE)*gpsha(1:VECTOR_SIZE,jnode,igaus)
               elmat2(1:VECTOR_SIZE,inode,jnode)=elmat2(1:VECTOR_SIZE,inode,jnode) +fact3(1:VECTOR_SIZE)*gpsha(1:VECTOR_SIZE,jnode,igaus)
            end do     
         end do
      end do



      trace1  = 0.0_rp
      trace2  = 0.0_rp
      elmass1 = 0.0_rp
      elmass2 = 0.0_rp
      do inode = 1,pnode                       
         trace1 = trace1 + elmat1(1:VECTOR_SIZE,inode,inode)
         trace2 = trace2 + elmat2(1:VECTOR_SIZE,inode,inode)
         do jnode = 1,pnode                       
            elmass1(1:VECTOR_SIZE) = elmass1(1:VECTOR_SIZE) + elmat1(1:VECTOR_SIZE,inode,jnode)
            elmass2(1:VECTOR_SIZE) = elmass2(1:VECTOR_SIZE) + elmat2(1:VECTOR_SIZE,inode,jnode)
         end do
      end do

       do inode = 1,pnode
          eldtrho(1:VECTOR_SIZE,inode) = elmat1(1:VECTOR_SIZE,inode,inode)*(elmass1(1:VECTOR_SIZE)/trace1(1:VECTOR_SIZE))
          elmurho(1:VECTOR_SIZE,inode) = elmat2(1:VECTOR_SIZE,inode,inode)*(elmass2(1:VECTOR_SIZE)/trace2(1:VECTOR_SIZE))
       end do
       ! AB TODO do inode = 1,pnode
       ! AB TODO    eldtrho(1:VECTOR_SIZE,inode) = eldtrho(1:VECTOR_SIZE,inode) * T_dtrho(1:VECTOR_SIZE) / d_dtrho(1:VECTOR_SIZE)
       ! AB TODO    eltau(1:VECTOR_SIZE,inode)   = eltau(1:VECTOR_SIZE,inode)   * T_tau(1:VECTOR_SIZE)   / d_tau(1:VECTOR_SIZE)
       ! AB TODO    elmurho(1:VECTOR_SIZE,inode) = elmurho(1:VECTOR_SIZE,inode) * T_murho(1:VECTOR_SIZE) / d_murho(1:VECTOR_SIZE)
       ! AB TODO end do
    end if

  end subroutine nsi_rhodt_rhotau_nu_vector

  subroutine nsi_rhodt_rhotau_nu_scalar(&
       pnode,pgaus,porde,gpsha,gpvol,gpden,gpvis,gppor,gpadv,&
       chale,dtinv_loc,densi,eldtrho,eltau,elmurho)

    use def_nastin, only : corio_nsi
    use def_nastin, only : kfl_taust_nsi
    use def_nastin, only : nbdfp_nsi
    use def_nastin, only : staco_nsi
    use mod_tauadr, only : tauadr

    integer(ip), intent(in)    :: pnode
    integer(ip), intent(in)    :: pgaus
    integer(ip), intent(in)    :: porde
    real(rp),    intent(in)    :: gpsha(pnode,pgaus)
    real(rp),    intent(in)    :: gpvol(pgaus)
    real(rp),    intent(in)    :: gpden(pgaus)
    real(rp),    intent(in)    :: gpvis(pgaus)
    real(rp),    intent(in)    :: gppor(pgaus)
    real(rp),    intent(in)    :: gpadv(ndime,pgaus)
    real(rp),    intent(in)    :: chale(2)
    real(rp),    intent(in)    :: dtinv_loc
    real(rp),    intent(in)    :: densi(pgaus,nbdfp_nsi)
    real(rp),    intent(out)   :: eldtrho(pnode)
    real(rp),    intent(out)   :: eltau(pnode)
    real(rp),    intent(out)   :: elmurho(pnode)

    integer(ip)                :: inode,igaus,idime
    real(rp)                   :: gpsp1(pgaus)
    real(rp)                   :: fact2
    real(rp)                   :: fact3
    real(rp)                   :: fact4
    real(rp)                   :: gpvno
    real(rp)                   :: adv
    real(rp)                   :: dif
    real(rp)                   :: rea
    real(rp)                   :: h2
    real(rp)                   :: d_dtrho
    real(rp)                   :: d_tau
    real(rp)                   :: d_murho
    real(rp)                   :: T_dtrho
    real(rp)                   :: T_tau
    real(rp)                   :: T_murho
    !
    ! tau1 = 1 / ( 4*mu/h^2 + 2*rho*uc/h + rho*|w| + sig )
    !
    do igaus = 1,pgaus

       gpvno = 0.0_rp
       do idime = 1,ndime
          gpvno = gpvno + gpadv(idime,igaus) * gpadv(idime,igaus)
       end do
       gpvno = sqrt(gpvno)
       adv   = gpden(igaus) * gpvno                          ! Convective term: rho*|u+u'|
       dif   = gpvis(igaus)                                  ! Viscous term: mu
       rea   = gpden(igaus) * corio_nsi + abs(gppor(igaus))  ! Coriolis: w + Porosity: sig
       h2    = chale(2) * chale(2)

       call tauadr(&
            kfl_taust_nsi,staco_nsi,adv,dif,rea,&
            chale(1),chale(2),gpsp1(igaus),dtinv_loc)

    end do
    !
    ! Element assembly
    !
    eldtrho = 0.0_rp
    eltau   = 0.0_rp
    elmurho = 0.0_rp
    if( porde == 1 ) then
       do igaus = 1,pgaus
          fact2 = gpvol(igaus) / ( gpden(igaus) * dtinv_loc )
          fact3 = gpvol(igaus) * gpsp1(igaus)
          fact4 = gpvol(igaus) * gpvis(igaus) / gpden(igaus)
          do inode = 1,pnode
             eldtrho(inode) = eldtrho(inode) + gpsha(inode,igaus) * fact2
             eltau(inode)   = eltau(inode)   + gpsha(inode,igaus) * fact3
             elmurho(inode) = elmurho(inode) + gpsha(inode,igaus) * fact4
          end do
       end do
    else
       T_dtrho = 0.0_rp
       d_dtrho = 0.0_rp
       T_tau   = 0.0_rp
       d_tau   = 0.0_rp
       T_murho = 0.0_rp
       d_murho = 0.0_rp
       do igaus = 1,pgaus
          fact2 = gpvol(igaus) / ( gpden(igaus) * dtinv_loc )
          fact3 = gpvol(igaus) * gpsp1(igaus)
          fact4 = gpvol(igaus) * gpvis(igaus) / gpden(igaus)
          T_dtrho = T_dtrho + fact2
          T_tau   = T_tau   + fact3
          T_murho = T_murho + fact4
          do inode = 1,pnode
             eldtrho(inode) = eldtrho(inode) + gpsha(inode,igaus)**2 * fact2
             eltau(inode)   = eltau(inode)   + gpsha(inode,igaus)**2 * fact3
             elmurho(inode) = elmurho(inode) + gpsha(inode,igaus)**2 * fact4
          end do
       end do
       do inode = 1,pnode
          d_dtrho = d_dtrho + eldtrho(inode)
          d_tau   = d_tau   + eltau(inode)
          d_murho = d_murho + elmurho(inode)
       end do
       do inode = 1,pnode
          eldtrho(inode) = eldtrho(inode) * T_dtrho / d_dtrho
          eltau(inode)   = eltau(inode)   * T_tau / d_tau
          elmurho(inode) = elmurho(inode) * T_murho / d_murho
       end do
    end if

   end subroutine nsi_rhodt_rhotau_nu_scalar

   !-----------------------------------------------------------------------
   !>
   !> @author  Guillaume Houzeaux
   !> @date    6/10/2016
   !> @brief   Laplacian
   !> @details Assemble a Laplacian matrix ( grad p , grad q )
   !>
   !-----------------------------------------------------------------------

   subroutine nsi_element_laplacian(&
        pnode,pgaus,gpcar,gpden,gpvol,elmap)

     use def_domain, only : mnode
     use def_nastin, only : KFL_SURTE_NSI

     integer(ip), intent(in)  :: pnode
     integer(ip), intent(in)  :: pgaus
     real(rp),    intent(in)  :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)
     real(rp),    intent(in)  :: gpden(VECTOR_SIZE,pgaus)
     real(rp),    intent(in)  :: gpvol(VECTOR_SIZE,pgaus)
     real(rp),    intent(out) :: elmap(VECTOR_SIZE,pnode,pnode)
     integer(ip)              :: inode,jnode,idime,igaus

     elmap = 0.0_rp

     if(kfl_surte_nsi == 2_ip) then
        do igaus = 1,pgaus
           do inode = 1,pnode
              do jnode = 1,pnode
                 do idime = 1,ndime
                    elmap(1:VECTOR_SIZE,inode,jnode) = elmap(1:VECTOR_SIZE,inode,jnode) &
                       & + (1.0_rp/gpden(1:VECTOR_SIZE,igaus))*gpcar(1:VECTOR_SIZE,idime,inode,igaus) &
                       & * gpcar(1:VECTOR_SIZE,idime,jnode,igaus) &
                       & * gpvol(1:VECTOR_SIZE,igaus) 
                 end do
              end do
           end do
        end do
     else
        do igaus = 1,pgaus
           do inode = 1,pnode
              do jnode = 1,pnode
                 do idime = 1,ndime
                    elmap(1:VECTOR_SIZE,inode,jnode) = elmap(1:VECTOR_SIZE,inode,jnode) &
                       & + gpcar(1:VECTOR_SIZE,idime,inode,igaus) &
                       & * gpcar(1:VECTOR_SIZE,idime,jnode,igaus) &
                       & * gpvol(1:VECTOR_SIZE,igaus) 
                 end do
              end do
           end do
        end do
     end if

   end subroutine nsi_element_laplacian

end module mod_nsi_element_operations
!> @}
