!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_chm_finiteRate_fast

  use def_kintyp, only   : ip,rp
  use def_domain, only   : ndime,npoin,coord,mgaus,mnode
  use def_chemic, only   : nclas_chm
  use def_chemic, only   : nspec_chm
  use def_master, only   : INOTMASTER
  use mod_memory, only   : memory_alloca, memory_deallo
  use def_master, only   : mem_modul,modul

  implicit none

  private

!@  public :: chm_element_operations_finiteRate_fast
!@  public :: chm_elmgac_finiteRate
!@  public :: chm_gatherProp_finiteRate
!@  public :: chm_elmpre_finiteRate
!@  public :: chm_turbul_finiteRate
!@  public :: chm_elmprc_finiteRate
!@  public :: chm_updtcc_finiteRate
!@  public :: chm_getProp_finiteRate
!@  public :: get_rates_pfa
!@  public :: get_reduced_reaction_rates
!@  public :: chm_IntegrateSource_finiteRate
!@  public :: chm_calc_enthalpy_transport_finiteRate
!@  public :: chm_calc_div_enthalpy_transport_finiteRate
!@  public :: get_correlation
!@  public :: get_mixture_fraction
!@  public :: get_prog_var
!@  public :: chm_heatRelease_finiteRate
contains

!@  subroutine chm_calc_hk_grad_Yk_others
!@    !
!@    !  Compute terms needed for enthalpy flux.
!@    !  o Simplified way:
!@    !  q_j  = - rho * D * grad( h )
!@    !      ~= - rho * D * grad( cp * T  + sum( hk^0 Y_k ) ) Poinsot's interpretation
!@    !         (with constant cp)
!@    !      ~= - lambda * grad( T ) - rho * D * hk^0 * sum( grad( Y_k ) )
!@    !
!@    !  In the simplified way we do:
!@    !
!@    !  q_j = - rho * D * grad( h )
!@    !      = - rho * D * grad( sum( hk Y_k ) )
!@    !      = - rho * sum( D_k * grad( hk Y_k ) )
!@    !      = - rho * sum( D_k * grad( hk Y_k ) ) + rho*D_therm*grad(h) - rho*D_therm*grad(h)
!@    !                           \____________/                 \_____/   \_________________/
!@    !                                  |                          |               |
!@    !                             grad_Yk_hk                   grad_h         in TEMPER
!@    !
!@    !
!@    !  q_j*  = - rho * sum( D_k * grad( hk Y_k ) ) +rho*D_therm*grad(h)
!@    !  -1 *  q_j  is calculated in chm_calc_enthalpy_transport_finiteRate
!@    !
!@    !
!@    !  o Detailed way:
!@    !  q_j  = - lambda * grad( T ) - rho * sum ( D_k * h_k * grad( Y_k ) )
!@    !
!@    !  h_k is the sensible plus chemical enthalpy for species k
!@    !  When computing q_j in this way make thermal diffusion zero in temper
!@    !
!@    !  chm_calc_enthalpy_transport_finiteRate calculates -q_j
!@    !  chm_calc_div_enthalpy_transport_finiteRate calculates -div(q_j)
!@    !
!@    use def_master,                 only : conce,tempe,prthe,therm,kfl_htran
!@#ifdef CANTERA
!@    use def_chemic,                 only : gas_chm
!@#endif
!@    use def_chemic,                 only : grad_Yk_hk,aux_nodes,grad_k_single,&
!@                                           grad_Yk, grad_T, entha_chm, grad_h, W_k
!@    use mod_gradie
!@    use def_kermod,                 only : gasco
!@
!@    implicit none
!@    integer(ip)                         :: iclas,dummi,ipoin
!@    real(rp)                            :: hkn(nspec_chm, npoin)
!@    real(rp)                            :: hk(nspec_chm)
!@    real(rp)                            :: tem
!@
!@
!@    if(kfl_htran == 0) then
!@
!@#ifdef CANTERA
!@
!@       do ipoin = 1,npoin
!@          !
!@          ! Set Yk,H,P at node
!@          !
!@          call setMassFractions(gas_chm,conce(ipoin,1:nspec_chm,1))
!@          call setState_HP(gas_chm,therm(ipoin,1),prthe(1))
!@
!@          !
!@          ! Get temperature (might be slightly different from tempe(ipoin,1)
!@          ! Get h_k_molar / (Ru*T) = a1 + a2 T /2 + a3 T^2 /3 + a4 T^3 /4 + a5 T^4 /5 + a6/T
!@          ! Cantera uses NASA polynomials of species
!@          !
!@          tem  = temperature(gas_chm)
!@          call getenthalpies_rt(gas_chm, hk(1:nspec_chm) )
!@
!@          !
!@          ! Get h_k_mass = h_k_molar / (Ru*T) * Ru * T / Wk
!@          !
!@          hk(1:nspec_chm) = hk(1:nspec_chm) * tem * gasco / W_k(1:nspec_chm)
!@
!@          !
!@          ! Get enthalpy fraction associated with each species h_k_mass * Y_k
!@          !
!@          do iclas=1,nspec_chm
!@             hkn(iclas,ipoin) = hk(iclas) * conce(ipoin,iclas,1)
!@          end do
!@       enddo
!@#endif
!@
!@       !
!@       ! Get gradients of h_k_mass * Y_k
!@       !
!@       do iclas=1,nspec_chm
!@          call gradie(hkn(iclas,1:npoin),grad_k_single)
!@          grad_Yk_hk(iclas,1:ndime,1:npoin) = grad_k_single(1:ndime,1:npoin)
!@       enddo
!@
!@       !
!@       ! Get gradient of -H
!@       !
!@       aux_nodes(1:npoin) = -1.0_rp * therm(1:npoin,1)
!@       call gradie(aux_nodes,grad_h)
!@
!@    else
!@       if (INOTMASTER) then
!@
!@          !
!@          ! Get gradients of Y_k and T
!@          !
!@          do iclas=1,nspec_chm
!@             call gradie(conce(1:npoin,iclas,1),grad_k_single)
!@             grad_Yk(iclas,1:ndime,1:npoin) = grad_k_single(1:ndime,1:npoin)
!@          enddo
!@
!@          !
!@          ! Get gradient of T
!@          !
!@          call gradie(tempe,grad_T)
!@
!@          !
!@          ! Get enthalpies
!@          !
!@#ifdef CANTERA
!@
!@          do ipoin = 1,npoin
!@             !
!@             ! Set Yk,H,P at node
!@             !
!@             call setMassFractions(gas_chm,conce(ipoin,1:nspec_chm,1))
!@             call setState_HP(gas_chm,therm(ipoin,1),prthe(1))
!@
!@             !
!@             ! Get temperature (might be slightly different from tempe(ipoin,1)
!@             ! Get h_k_molar / (Ru*T) = a1 + a2 T /2 + a3 T^2 /3 + a4 T^3 /4 + a5 T^4 /5 + a6/T
!@             ! Cantera uses NASA polynomials of species
!@             !
!@             tem  = temperature(gas_chm)
!@             !! entha_chm also contains soot enthalpies hs (nclas_chm,npoin)
!@             call getenthalpies_rt(gas_chm, entha_chm(1:nspec_chm,ipoin) )
!@
!@             !
!@             ! Get h_k_mass = h_k_molar / (Ru*T) * Ru * T / Wk
!@             !
!@             entha_chm(1:nspec_chm,ipoin) = entha_chm(1:nspec_chm,ipoin) * tem * gasco / W_k(1:nspec_chm)
!@
!@          enddo
!@#endif
!@
!@       end if
!@    endif
!@
!@  end subroutine chm_calc_hk_grad_Yk_others
!@
!@
  !@subroutine chm_element_operations_finiteRate_fast(order,pnode,pgaus,list_elements)

!@    use def_parame
!@    use def_master
!@    use def_domain
!@    use def_chemic
!@    use mod_ker_proper
!@    use def_kermod
!@    use mod_matrix
!@    use def_solver
!@    use mod_chm_finiteRate
!@
!@    use mod_ADR,    only : ADR_element_assembly
!@    use mod_ADR,    only : ADR_projections_and_sgs_assembly
!@    use mod_ADR,    only : ADR_add_sgs_or_bubble
!@    use mod_ADR,    only : ELEMENT_ASSEMBLY                 ! 1
!@    use mod_ADR,    only : PROJECTIONS_AND_SGS_ASSEMBLY     ! 4
!@    use mod_ADR,    only : BUBBLE_ASSEMBLY                  ! 5
!@    use mod_ADR,    only : mreac_adr
!@    use mod_chm_spray, only : chm_elmprc_spray
!@    use mod_chm_spray, only : chm_rhodt_spray
!@    use mod_chm_spray, only : chm_elmpre_spray
!@    use mod_chm_entropy, only : chm_entropy_viscosity
!@
!@    implicit none
!@
!@    integer(ip), intent(in)          :: order                                    ! =1 defaul or =2 compute SGS only
!@    integer(ip), intent(in)          :: pnode                                    !< Number of nodes
!@    integer(ip), intent(in)          :: pgaus                                    !< Number of Gauss points
!@    integer(ip), intent(in), pointer :: list_elements(:)                         !< List of elements
!@
!@    integer(ip)                      :: ielem,iclas,igaus
!@    integer(ip)                      :: izmat,izrhs,pelty
!@    integer(ip)                      :: porde,ptopo,plapl
!@    integer(ip)                      :: dummi,kelem
!@
!@    real(rp)                         :: elmat(mnode,mnode)
!@    real(rp)                         :: elrhs(mnode)
!@    real(rp)                         :: elcon(mnode,nclas_chm)! <=> conce
!@    real(rp)                         :: elcod(ndime,mnode)                       ! <=> coord
!@    real(rp)                         :: elvel(ndime,mnode)
!@    real(rp)                         :: gpvol(mgaus)                             ! |J|*w
!@    real(rp)                         :: gphco(mgaus)                             ! heat conductivity
!@    real(rp)                         :: gpdit(mgaus)                             ! thermal diffusivity
!@    real(rp)                         :: gpsph(mgaus)                             ! specific heat
!@    real(rp)                         :: gpcon(mgaus,nclas_chm,ncomp_chm)         ! <=> conce
!@    real(rp)                         :: gprea(mgaus,mreac_adr)                   ! r
!@    real(rp)                         :: gpvel(ndime,mgaus)                       ! u
!@    real(rp)                         :: gpdif(mgaus,nclas_chm)                   ! D_k
!@    real(rp)                         :: del_gpdif(mgaus,nclas_chm)               ! change in diffusivity from entropy stable
!@                                                                                 !  stabilization
!@    real(rp)                         :: gpgrd(ndime,mgaus)                       ! grad(k) = grad(D_k)
!@    real(rp)                         :: gprhs(mgaus,nclas_chm)                   ! f (all terms)
!@    real(rp)                         :: gpden(mgaus)                             ! fake rho for elmadr
!@    real(rp)                         :: gptur(mgaus)                             ! turbulent viscosity
!@    real(rp)                         :: gpdiv(mgaus)                             ! Divergence of convection
!@    real(rp)                         :: gpcar(ndime,mnode,mgaus)                 ! dNk/dxj
!@    real(rp)                         :: gphes(ntens,mnode,mgaus)                 ! dNk/dxidxj
!@    real(rp)                         :: gpsigma(mgaus)                           ! RHS of liquid gas surface density equation
!@    real(rp)                         :: dummr(mgaus*ndime)
!@    real(rp)                         :: chale(3),chave(3),hleng(3),tragl(9)
!@    real(rp)                         :: grad_Yk_hk_node(mnode,nspec_chm,ndime)
!@    real(rp)                         :: grad_Yk_hk_gp(mgaus,nspec_chm,ndime)
!@    real(rp)                         :: grad_h_node(mnode,ndime)
!@    real(rp)                         :: grad_h_gp(mgaus,ndime)
!@    real(rp)                         :: grad_Yk_node(mnode,nspec_chm,ndime)
!@    real(rp)                         :: entha_node(mnode,nclas_chm)
!@    real(rp)                         :: hk_grad_Yk_gp(mgaus,nspec_chm,ndime)
!@    real(rp)                         :: grad_T_node(mnode,ndime)
!@    real(rp)                         :: grad_T_gp(mgaus,ndime)
!@
!@    !
!@    ! Loop over elements
!@    !
!@    elements: do kelem = 1,size(list_elements)
!@
!@       ielem = list_elements(kelem)
!@
!@       if( ielem > 0 ) then
!@          !
!@          ! Element dimensions
!@          !
!@          pelty = ltype(ielem)
!@          porde = lorde(pelty)
!@          ptopo = ltopo(pelty)
!@
!@          !
!@          ! Initialization variables
!@          !
!@          gpdiv         = 0.0_rp
!@          gpdif         = 0.0_rp
!@          gpdit         = 0.0_rp
!@          del_gpdif     = 0.0_rp
!@          gprea         = 0.0_rp
!@          gptur         = 0.0_rp
!@          gpcon         = 0.0_rp
!@          gpsigma       = 0.0_rp
!@          gpgrd         = 0.0_rp
!@          gpvel         = 0.0_rp
!@          gprhs         = 0.0_rp
!@          gphes         = 0.0_rp
!@          grad_Yk_hk_gp = 0.0_rp
!@          grad_h_gp     = 0.0_rp
!@          hk_grad_Yk_gp = 0.0_rp
!@          grad_T_gp     = 0.0_rp
!@          plapl         = 0
!@
!@          !
!@          ! Gather all
!@          ! elQ = Y_k,U,h, +
!@          !        (grad_Y_k*grad_h),grad_Y_k,grad_h,grad_T)
!@          !
!@          call chm_elmgac_finiteRate(&
!@                   pnode,lnods(1:pnode,ielem),elcod,elcon(1:pnode,1:nclas_chm),&
!@                   elvel,grad_Yk_hk_node,grad_h_node,grad_Yk_node,entha_node,grad_T_node)
!@
!@
!@          !
!@          ! CHALE, HLENG and TRAGL
!@          !
!@          if( kfl_taust_chm /= 0 ) then
!@             call elmlen(&
!@                  ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),hleng)
!@             call elmchl(&
!@                  tragl,hleng,elcod,dummr,chave,chale,pelty,pnode,porde,hnatu(pelty),&
!@                  kfl_advec_chm,kfl_ellen_chm)
!@          else
!@             plapl = 0
!@          end if
!@
!@          !
!@          ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
!@          !
!@          call elmcar(&
!@               pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
!@               elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
!@               gphes,ielem)
!@
!@          call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,elmar(pelty)%shape,gpcar)
!@          call ker_proper('CONDU','PGAUS',dummi,ielem,gphco,pnode,pgaus,elmar(pelty)%shape,gpcar)
!@          call ker_proper('SPHEA','PGAUS',dummi,ielem,gpsph,pnode,pgaus,elmar(pelty)%shape,gpcar)
!@          call ker_proper('TURBU','PGAUS',dummi,ielem,gptur,pnode,pgaus,elmar(pelty)%shape,gpcar)
!@          gpdit(1:pgaus) = gphco(1:pgaus) / gpsph(1:pgaus)
!@
!@          !
!@          ! Send quantities to gauss points
!@          ! gpQ = gpY_k,gpU, +
!@          !       gp_(h_grad_Y_k),gp_(grad_Y_k*grad_h),gp_(grad_Y_k),gp_(grad_h),gp_(grad_T), +
!@          !       elY_k,elU
!@          !
!@          call chm_elmpre_finiteRate(&
!@                   pnode,pgaus,elcon(1:pnode,1:nclas_chm),elvel,elmar(pelty)%shape,&
!@                   gpcon(1:pgaus,1:nclas_chm,:),gpvel,grad_Yk_hk_node,grad_h_node,grad_Yk_hk_gp,grad_h_gp, &
!@                   grad_Yk_node,entha_node,grad_T_node,hk_grad_Yk_gp,grad_T_gp)
!@
!@          !
!@          ! Compute molecular transport
!@          ! TODO: Probably should be after entropy viscosity
!@          !
!@          call chm_elmprc_finiteRate(&
!@                   pgaus,gphco,gpsph,gpdif) ! At this point gpdif only contains the laminar diffusivity
!@
!@          !
!@          ! Compute enthalpy transport
!@          !
!@          call chm_calc_enthalpy_transport_finiteRate(&
!@                   ielem,pgaus,gpdif,gpdit,gphco,grad_Yk_hk_gp,grad_h_gp,hk_grad_Yk_gp,grad_T_gp)
!@
!@          !
!@          ! Entropy stable viscosity
!@          !
!@          if(kfl_entropy_chm == 1_ip) then
!@             del_gpdif = 0.0_rp
!@             do iclas = iclai_chm,iclaf_chm
!@                call chm_entropy_viscosity(ielem,pnode,pgaus,1_ip,pgaus,iclas,elvel,gpden,hleng,del_gpdif)
!@             enddo
!@
!@             do iclas = iclai_chm,iclaf_chm
!@                gpdif(1:pgaus,iclas) = gpdif(1:pgaus,iclas) + del_gpdif(1:pgaus,iclas)
!@             enddo ! At this point gpdif contains the laminar and entropy viscosity contributions
!@          end if
!@
!@          ! At this point gpdif contains the laminar, turbulent and entropy visc. contributions
!@          call chm_turbul_finiteRate(pgaus,gpdif,gpden,gptur)
!@          !
!@          ! RHS calculation of liquid surface density at gauss points
!@          !
!@          if ( kfl_spray_chm /= 0_ip ) &
!@              call chm_elmpre_spray(ielem,pnode,pgaus,elvel,gpcar,gpcon(1:pgaus,:,1),&
!@                   gpden,gptur,hleng,gpsigma)
!@
!@          !
!@          ! Projections of rho/dt and 1/dt
!@          !
!@          if ( kfl_spray_chm == 0_ip )then
!@
!@             !
!@             ! Gas phase only
!@             !
!@             call chm_rhodt(  &
!@                  pnode,pgaus,porde,lnods(:,ielem),elmar(pelty)%shape,gpvol,gpden,dt_rho_chm)
!@
!@          else
!@
!@             !
!@             ! Liquid and gas phase for ELSA model
!@             !
!@             call chm_rhodt_spray(  &
!@                    pnode,pgaus,porde,lnods(:,ielem),elmar(pelty)%shape,gpvol,gpden,dt_rho_chm,dt_chm)
!@          endif
!@
!@
!@          !
!@          ! Calculate RHS soot source terms
!@          !
!@          !! call chm_assembly_soot_sources_ssm(ielem,pgaus,gpcon,gptem,gpden,gprhs,gpsoot)
!@
!@          !
!@          ! Assemble matrix
!@          !
!@          izmat = 1
!@          izrhs = 1
!@
!@          ASSEMBLY_ICLAS: do iclas = iclai_chm,iclaf_chm
!@
!@             !
!@             ! Assembly spray terms
!@             !
!@             if (kfl_spray_chm /= 0_ip .and. iclas > (nclas_chm - 2_ip) ) then
!@                call chm_elmprc_spray(&
!@                         iclas,pgaus,gptur,gpsigma,gpden,gpdif(1:pgaus,iclas),gprhs(1:pgaus,iclas))
!@
!@             end if
!@
!@             !
!@             ! Assemble equation for iclas
!@             !
!@             if( order == ELEMENT_ASSEMBLY ) then
!@
!@                !DMM: Revise communicating variables style
!@                call ADR_element_assembly(&
!@                     ielem,pnode,pgaus,elcod,elmar(pelty)%shape,gpcar,elmar(pelty)%deriv,gphes,gpvol,chale,&
!@                     elmar(pelty)%shape_bub,elmar(pelty)%deriv_bub,ADR_chm(iclas),cutim,gpden,gpvel,gpdif(1:pgaus,iclas),&
!@                     gpgrd,gprea,gprhs(1:pgaus,iclas),gpcon(1:pgaus,iclas:iclas,:),elcon(1:pnode,iclas:iclas),elmat,elrhs)
!@
!@
!@                call matrix_assexp(solve(1)%ndofn,1_ip,pnode,npoin,lnods(1:pnode,ielem),elrhs,elmat, &
!@                                   elcon(1:pnode,iclas:iclas),rhsid,iclas)
!@
!@                izrhs = izrhs + npoin                                       !solve(1)%nzrhs
!@                izmat = izmat + solve(1)%nzmat/nclas_chm**2
!@
!@             end if
!@
!@          end do ASSEMBLY_ICLAS
!@
!@       end if
!@
!@    end do elements

!@  end subroutine chm_element_operations_finiteRate_fast
!@
!@   subroutine chm_updtcc_finiteRate(dtmin)
!@     !-----------------------------------------------------------------------
!@     !****f* Chemic/chm_updtcc_finiteRate
!@     ! NAME
!@     !    chm_updtcc_finiteRate
!@     ! DESCRIPTION
!@     !    This routine computes the critical time step size in flamelet models
!@     ! USED BY
!@     !    chm_updtss
!@     !***
!@     !-----------------------------------------------------------------------
!@     use def_parame
!@     use def_master
!@     use def_domain
!@     use def_chemic
!@     use mod_ker_proper
!@     use def_kermod
!@     use mod_ADR,            only : ADR_critical_time_step
!@     use mod_ADR,            only : mreac_adr
!@     use mod_ADR,            only : FROM_CRITICAL
!@     use mod_communications, only : PAR_MIN
!@
!@     use mod_chm_spray, only : chm_elmprc_spray
!@     use mod_chm_spray, only : chm_rhodt_spray
!@     use mod_chm_spray, only : chm_elmpre_spray
!@
!@     implicit none
!@     real(rp),   intent(inout) :: dtmin
!@     real(rp)                  :: dtcri(2)
!@     integer(ip) :: ielem,iclas               ! Indices and dimensions
!@     integer(ip) :: pelty,pnode
!@     integer(ip) :: pgaus,plapl,porde,ptopo,dummi
!@
!@     real(rp)    :: elcon(mnode,nclas_chm)                   ! <=> Y_k
!@     real(rp)    :: elcod(ndime,mnode)                       ! <=> coord
!@     real(rp)    :: elvel(ndime,mnode)
!@     real(rp)    :: gpvol(mgaus)                             ! |J|*w
!@     real(rp)    :: gpcon(mgaus,nclas_chm,ncomp_chm)         ! <=> conce
!@     real(rp)    :: gprea(mgaus,mreac_adr)                   ! r
!@     real(rp)    :: gpvel(ndime,mgaus)                       ! u
!@     real(rp)    :: gpdif(mgaus,nclas_chm)                   ! D_k
!@     real(rp)    :: gpgrd(ndime,mgaus)                       ! grad(k) = grad(D_k)
!@     real(rp)    :: gprhs(mgaus,nclas_chm)                   ! f (all terms)
!@     real(rp)    :: gpden(mgaus)                             ! rho
!@     real(rp)    :: gpdiv(mgaus)                             ! Divergence of convection
!@     real(rp)    :: gpcar(ndime,mnode,mgaus)                 ! dNk/dxj
!@     real(rp)    :: gphes(ntens,mnode,mgaus)                 ! dNk/dxidxj
!@     real(rp)    :: gphco(mgaus)                             ! heat conductivity
!@     real(rp)    :: gpsph(mgaus)                             ! specific heat
!@     real(rp)    :: gptur(mgaus)
!@     real(rp)    :: gpsigma(mgaus)                           ! RHS of liquid gas surface density equation
!@
!@     real(rp)    :: dummr_node(mnode,nspec_chm,ndime)        ! For grad(hk*Yk) at nodes
!@     real(rp)    :: dumr2_node(mnode,ndime)                  ! For grad(h) at nodes
!@     real(rp)    :: dumr3_node(mnode,nspec_chm,ndime)        ! For grad(Yk) at nodes
!@     real(rp)    :: dumr4_node(mnode,nclas_chm)              ! For enthalpy at nodes
!@     real(rp)    :: dumr5_node(mnode,ndime)                  ! For grad(T) at nodes
!@     real(rp)    :: dummr_gp(mgaus,nspec_chm,ndime)          ! For grad(hk*Yk) at Gaussian points
!@     real(rp)    :: dumr2_gp(mgaus,ndime)                    ! For grad(h) at Gaussian points
!@     real(rp)    :: dumr3_gp(mgaus,nspec_chm,ndime)          ! For hk*grad(Yk) at Gaussian points
!@     real(rp)    :: dumr4_gp(mgaus,ndime)                    ! For grad(T) at Gaussian points
!@     real(rp)    :: dummr(mgaus*ndime)
!@     real(rp)    :: chale(3),chave(3),hleng(3),tragl(9)
!@
!@
!@     if( INOTMASTER ) then
!@
!@        do ielem = 1,nelem
!@
!@           !
!@           ! Element dimensions
!@           !
!@           pelty = ltype(ielem)
!@           pnode = nnode(pelty)
!@           pgaus = ngaus(pelty)
!@           porde = lorde(pelty)
!@           ptopo = ltopo(pelty)
!@
!@           !
!@           ! Initialization
!@           !
!@           gpdiv      = 0.0_rp
!@           gpdif      = 0.0_rp
!@           gprea      = 0.0_rp
!@           gptur      = 0.0_rp
!@           gpcon      = 0.0_rp
!@           gpsigma    = 0.0_rp
!@           gpgrd      = 0.0_rp
!@           gpvel      = 0.0_rp
!@           gprhs      = 0.0_rp
!@           gphes      = 0.0_rp
!@           dummr_gp   = 0.0_rp
!@           dumr2_gp   = 0.0_rp
!@           dumr3_gp   = 0.0_rp
!@           dumr4_gp   = 0.0_rp
!@
!@           !
!@           ! Gather all
!@           !
!@           call chm_elmgac_finiteRate(&
!@                    pnode,lnods(1:pnode,ielem),elcod,elcon(1:pnode,1:nclas_chm),&
!@                    elvel,dummr_node,dumr2_node,dumr3_node,dumr4_node,dumr5_node)
!@
!@           !
!@           ! CHALE, HLENG and TRAGL
!@           !
!@           call elmlen(&
!@                ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),hleng)
!@
!@           call elmchl(&
!@                tragl,hleng,elcod,dummr,chave,chale,pelty,pnode,porde,hnatu(pelty),&
!@                kfl_advec_chm,kfl_ellen_chm)
!@
!@           !
!@           ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
!@           !
!@           call elmcar(&
!@                pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
!@                elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
!@                gphes,ielem)
!@
!@           !
!@           ! Transport properties and density
!@           !
!@           call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,elmar(pelty) % shape,gpcar)
!@           call ker_proper('CONDU','PGAUS',dummi,ielem,gphco,pnode,pgaus,elmar(pelty) % shape,gpcar)
!@           call ker_proper('SPHEA','PGAUS',dummi,ielem,gpsph,pnode,pgaus,elmar(pelty) % shape,gpcar)
!@           call ker_proper('TURBU','PGAUS',dummi,ielem,gptur,pnode,pgaus,elmar(pelty) % shape,gpcar)
!@
!@           !
!@           ! Send quantities to gauss points
!@           !
!@           call chm_elmpre_finiteRate(&
!@                    pnode,pgaus,elcon(1:pnode,1:nclas_chm),elvel,elmar(pelty)%shape,&
!@                    gpcon(1:pgaus,1:nclas_chm,:),gpvel,dummr_node,dumr2_node,dummr_gp,dumr2_gp, &
!@                    dumr3_node,dumr4_node,dumr5_node,dumr3_gp,dumr4_gp)
!@
!@
!@           call chm_elmprc_finiteRate(&
!@                    pgaus,gphco,gpsph,gpdif)
!@
!@           !
!@           ! RHS calculation of liquid surface density at gauss points
!@           !
!@           if ( kfl_spray_chm /= 0_ip ) &
!@               call chm_elmpre_spray(ielem,pnode,pgaus,elvel,gpcar,gpcon(1:pgaus,:,1),&
!@                        gpden,gptur,hleng,gpsigma)
!@
!@           !
!@           ! Loop over variables
!@           !
!@           do iclas = iclai_chm,iclaf_chm
!@
!@              !
!@              ! Assembly spray terms
!@              !
!@              if (kfl_spray_chm /= 0_ip .and. iclas > (nclas_chm - 2_ip) ) then
!@                 call chm_elmprc_spray(iclas,pgaus,gptur,gpsigma,gpden,gpdif(1:pgaus,iclas),gprhs(1:pgaus,iclas))
!@
!@              end if
!@
!@              !
!@              ! Compute time-step
!@              !
!@              call ADR_critical_time_step(ADR_chm(iclas),gpden,gpvel,gpdif(1:pgaus,iclas),gprea,dtcri,chale(1),chale(2))
!@
!@              !
!@              ! Take minimum
!@              !
!@              dtmin = min(dtmin,dtcri(1))
!@
!@           end do
!@
!@        end do
!@
!@     end if
!@     !
!@     ! Look for minimum over subdomains
!@     !
!@     call PAR_MIN(dtmin,'IN MY CODE')
!@
!@  end subroutine chm_updtcc_finiteRate
!@
!@
!@  subroutine chm_elmgac_finiteRate( &
!@                 pnode,lnods,elcod,elcon,elvel,grad_Yk_hk_node,grad_h_node, &
!@                 grad_Yk_node,entha_node,grad_T_node)
!@     !------------------------------------------------------------------------
!@     !****f* Chemic/mod_chm_element_operations/chm_elmgac_finiteRate
!@     ! NAME
!@     !    chm_elmgac_finiteRate
!@     ! DESCRIPTION
!@     !    Gather operations for the combustion models
!@     ! USES
!@     ! USED BY
!@     !    chm_element_operations_finiteRate
!@     !***
!@     !------------------------------------------------------------------------
!@     use def_kintyp, only     :  ip,rp
!@     use def_master, only     :  conce,advec,kfl_htran
!@     use def_chemic, only     :  kfl_advec_chm
!@     use def_chemic, only     :  grad_Yk_hk, grad_h, grad_Yk, grad_T, &
!@                                 entha_chm
!@
!@     implicit none
!@     integer(ip), intent(in)  :: pnode
!@     integer(ip), intent(in)  :: lnods(pnode)
!@     real(rp),    intent(out) :: elcod(ndime,pnode)
!@     real(rp),    intent(out) :: elcon(pnode,nclas_chm)
!@     real(rp),    intent(out) :: elvel(ndime,pnode)
!@     real(rp),    intent(out) :: grad_Yk_hk_node(pnode,nspec_chm,ndime)
!@     real(rp),    intent(out) :: grad_h_node(pnode,ndime)
!@     real(rp),    intent(out) :: grad_Yk_node(pnode,nspec_chm,ndime)
!@     real(rp),    intent(out) :: grad_T_node(pnode,ndime)
!@     real(rp),    intent(out) :: entha_node(pnode,nspec_chm)
!@
!@     integer(ip)              :: inode,ipoin,iclas
!@
!@     !
!@     ! Initialization
!@     !
!@     elvel            = 0.0_rp
!@     elcod            = 0.0_rp
!@     elcon            = 0.0_rp
!@     grad_Yk_hk_node  = 0.0_rp
!@     grad_h_node      = 0.0_rp
!@     grad_Yk_node     = 0.0_rp
!@     grad_T_node      = 0.0_rp
!@     entha_node       = 0.0_rp
!@
!@     !
!@     ! Species mass fractions Yk, grad Yk and coordinates
!@     !
!@     if (kfl_htran == 0) then
!@        do inode=1,pnode
!@           ipoin=lnods(inode)
!@           do iclas=1,nspec_chm
!@              elcon(inode,iclas)                   = conce(ipoin,iclas,1)
!@              grad_Yk_hk_node(inode,iclas,1:ndime) = grad_Yk_hk(iclas,1:ndime,ipoin)
!@           end do
!@           grad_h_node(inode,1:ndime) = grad_h(1:ndime,ipoin)
!@           elcod(1:ndime,inode)                 = coord(1:ndime,ipoin)
!@        end do
!@
!@     else
!@        do inode=1,pnode
!@           ipoin=lnods(inode)
!@           do iclas=1,nspec_chm
!@              elcon(inode,iclas)                 = conce(ipoin,iclas,1)
!@              grad_Yk_node(inode,iclas,1:ndime)  = grad_Yk(iclas,1:ndime,ipoin)
!@              entha_node(inode,iclas)            = entha_chm(iclas,ipoin)
!@           end do
!@           grad_T_node(inode,1:ndime) = grad_T(1:ndime,ipoin)
!@           elcod(1:ndime,inode)       = coord(1:ndime,ipoin)
!@        end do
!@     end if
!@
!@     !
!@     ! Advection
!@     !
!@     if( kfl_advec_chm /= 0 ) then
!@        do inode = 1,pnode
!@           ipoin = lnods(inode)
!@           elvel(1:ndime,inode) = advec(1:ndime,ipoin,1)
!@        end do
!@     end if
!@
!@  end subroutine chm_elmgac_finiteRate
!@
!@  subroutine chm_elmpre_finiteRate(&
!@                 pnode,pgaus,elcon,elvel,gpsha,&
!@                 gpcon,gpvel,grad_Yk_hk_node,grad_h_node,grad_Yk_hk_gp,grad_h_gp, &
!@                 grad_Yk_node,entha_node,grad_T_node,hk_grad_Yk_gp,grad_T_gp)
!@
!@     !------------------------------------------------------------------------
!@     !****f* Chemic/mod_chm_element_operations/chm_elmpre_finiteRate
!@     ! NAME
!@     !    chm_elmpre_finiteRate
!@     ! DESCRIPTION
!@     !    Gather operations for the combustion models
!@     ! USES
!@     ! USED BY
!@     !    chm_element_operations_finiteRate
!@     !***
!@     !------------------------------------------------------------------------
!@     use def_kintyp, only      :  ip,rp
!@     use def_master, only      :  kfl_htran
!@
!@     implicit none
!@     integer(ip),  intent(in)  :: pnode,pgaus
!@     real(rp),     intent(in)  :: elcon(pnode,nclas_chm)
!@     real(rp),     intent(in)  :: elvel(ndime,pnode)
!@     real(rp),     intent(in)  :: gpsha(pnode,pgaus)
!@     real(rp),     intent(in)  :: grad_Yk_hk_node(pnode,nspec_chm,ndime)
!@     real(rp),     intent(in)  :: grad_h_node(pnode,ndime)
!@     real(rp),     intent(in)  :: grad_Yk_node(pnode,nspec_chm,ndime)
!@     real(rp),     intent(in)  :: grad_T_node(pnode,ndime)
!@     real(rp),     intent(in)  :: entha_node(pnode,nclas_chm)
!@     real(rp),     intent(out) :: gpcon(pgaus,nclas_chm,*)
!@     real(rp),     intent(out) :: grad_Yk_hk_gp(pgaus,nspec_chm,ndime)
!@     real(rp),     intent(out) :: grad_h_gp(pgaus,ndime)
!@     real(rp),     intent(out) :: hk_grad_Yk_gp(pgaus,nspec_chm,ndime) ! Product of enthalpy and Yk gradient
!@     real(rp),     intent(out) :: grad_T_gp(pgaus,ndime)
!@     real(rp),     intent(out) :: gpvel(ndime,pgaus)
!@
!@     integer(ip)               :: iclas,igaus,inode,idime
!@     real(rp)                  :: grad_Yk_gp(pgaus,nspec_chm,ndime), entha_k_gp(pgaus,nspec_chm)
!@
!@     if (kfl_htran == 0) then
!@        !
!@        ! Concentration
!@        !
!@        grad_Yk_hk_gp = 0.0_rp
!@        do iclas = 1,nspec_chm
!@           do igaus = 1,pgaus
!@              gpcon(igaus,iclas,1) = 0.0_rp
!@              do inode = 1,pnode
!@                 gpcon(igaus,iclas,1)               = gpcon(igaus,iclas,1) &
!@                                                      + gpsha(inode,igaus) * elcon(inode,iclas)
!@
!@                 grad_Yk_hk_gp(igaus,iclas,1:ndime) = grad_Yk_hk_gp(igaus,iclas,1:ndime) &
!@                                                      + gpsha(inode,igaus) * grad_Yk_hk_node(inode,iclas,1:ndime)
!@
!@              end do
!@           end do
!@        end do
!@
!@        !
!@        ! Gradient of enthalpy
!@        !
!@        grad_h_gp = 0.0_rp
!@        do igaus = 1,pgaus
!@           do inode = 1,pnode
!@              grad_h_gp(igaus,1:ndime) = grad_h_gp(igaus,1:ndime) &
!@                                              + gpsha(inode,igaus) * grad_h_node(inode,1:ndime)
!@           end do
!@        end do
!@
!@     else
!@        !
!@        ! Concentration
!@        !
!@        grad_Yk_gp = 0.0_rp
!@        entha_k_gp = 0.0_rp
!@        do iclas = 1,nspec_chm
!@           do igaus = 1,pgaus
!@              gpcon(igaus,iclas,1) = 0.0_rp
!@              do inode = 1,pnode
!@                 gpcon(igaus,iclas,1)            = gpcon(igaus,iclas,1) &
!@                                                 + gpsha(inode,igaus) * elcon(inode,iclas)
!@
!@                 grad_Yk_gp(igaus,iclas,1:ndime) = grad_Yk_gp(igaus,iclas,1:ndime) &
!@                                                 + gpsha(inode,igaus) * grad_Yk_node(inode,iclas,1:ndime)
!@
!@                 entha_k_gp(igaus,iclas)         = entha_k_gp(igaus,iclas) &
!@                                                 + gpsha(inode,igaus) * entha_node(inode,iclas)
!@
!@              end do
!@           end do
!@        end do
!@
!@        do idime = 1,ndime
!@           hk_grad_Yk_gp(1:pgaus,1:nspec_chm,idime) = entha_k_gp(1:pgaus,1:nspec_chm) * grad_Yk_gp(1:pgaus,1:nspec_chm,idime)
!@        end do
!@
!@        !
!@        ! Gradient of temperature
!@        !
!@        grad_T_gp = 0.0_rp
!@        do igaus = 1,pgaus
!@           do inode = 1,pnode
!@              grad_T_gp(igaus,1:ndime) = grad_T_gp(igaus,1:ndime) &
!@                                              + gpsha(inode,igaus) * grad_T_node(inode,1:ndime)
!@           end do
!@        end do
!@     end if
!@
!@     !
!@     ! Fluid velocity
!@     !
!@     gpvel = 0.0_rp
!@     do igaus = 1, pgaus
!@        do inode = 1,pnode
!@           do idime = 1,ndime
!@              gpvel(idime,igaus) = gpvel(idime,igaus) &
!@                   +  gpsha(inode,igaus) * elvel(idime,inode)
!@           end do
!@        enddo
!@     enddo
!@
!@
!@  end subroutine chm_elmpre_finiteRate
!@
!@
!@  subroutine chm_turbul_finiteRate(pgaus,gpdif,gpden,gptur)
!@    !-----------------------------------------------------------------------
!@    !****f* Chemic/mod_chm_element_operations/chm_turbul
!@    ! NAME
!@    !    chm_turbul_finiteRate
!@    ! DESCRIPTION
!@    !    Add subgrid contribution to diffusion term
!@    ! USES
!@    ! USED BY
!@    !    chm_turbul_finiteRate
!@    !***
!@    !-----------------------------------------------------------------------
!@    use def_chemic, only      : diffu_chm
!@
!@    implicit none
!@
!@    integer(ip),  intent(in)  :: pgaus
!@    real(rp),     intent(in)  :: gpden(pgaus)
!@    real(rp),     intent(in)  :: gptur(pgaus)
!@    real(rp),     intent(out) :: gpdif(pgaus,nclas_chm)
!@    integer(ip)               :: iclas
!@
!@    !
!@    ! Subgrid contribution for LES
!@    !
!@    if(turmu_ker % kfl_exist /= 0_ip) then
!@       do iclas = 1,nclas_chm
!@          gpdif(1:pgaus,iclas) = gpdif(1:pgaus,iclas) + gptur(1:pgaus) * gpden(1:pgaus) / diffu_chm(1,1)
!@       enddo
!@    end if
!@
!@  end  subroutine chm_turbul_finiteRate
!@
!@
!@  subroutine chm_elmprc_finiteRate(&
!@                 pgaus,gphco,gpsph,gpdif)
!@    !-----------------------------------------------------------------------
!@    !****f* Chemic/mod_chm_element_operations/chm_elmprc_finiteRate
!@    ! NAME
!@    !    chm_elmprc_finiteRate
!@    ! DESCRIPTION
!@    !    Compute laminar diffusivity mass coef. for each species ADR equation
!@    ! USES
!@    ! USED BY
!@    !    chm_element_operations_finiteRate
!@    !***
!@    !-----------------------------------------------------------------------
!@    use def_chemic, only      : diffu_chm, Le_k
!@
!@    implicit none
!@
!@    integer(ip),  intent(in)  :: pgaus
!@    real(rp),     intent(in)  :: gphco(pgaus)
!@    real(rp),     intent(in)  :: gpsph(pgaus)
!@    real(rp),     intent(out) :: gpdif(pgaus,nclas_chm)
!@
!@    integer(ip)               :: iclas
!@    real(rp)                  :: diff(pgaus)
!@
!@    !
!@    ! Diffusion coefficient
!@    !
!@    diff = gphco(1:pgaus) / gpsph(1:pgaus)
!@
!@    do iclas = 1,nclas_chm
!@       gpdif(1:pgaus,iclas) = diff(1:pgaus) / Le_k(iclas)
!@    enddo
!@
!@  end subroutine chm_elmprc_finiteRate
!@
!@
!@  subroutine chm_getProp_finiteRate
!@    !-----------------------------------------------------------------------
!@    !****f* chemic/getProp_finiteRate
!@    ! NAME
!@    !    getProp_finiteRate
!@    ! DESCRIPTION
!@    !    Compute properties for finite rate combustion model
!@    ! USES
!@    ! USED BY
!@    !    chm_iniunk: initialize flow field
!@    !    chm_endite:
!@    !    chm_endste: properties available for all modules at the end of time-step
!@    !                (there is no update of properties during chemic iteration)
!@    !-----------------------------------------------------------------------
!@    use def_kintyp, only      : ip,rp
!@    use def_master, only      : condu_gp,visco_gp,wmean_gp,prthe,&
!@                                sphec_gp_ht,sphec_gp_lt,tempe,&
!@                                conce,sphec, sphec_gp
!@    use def_kermod, only      : gasco
!@
!@    use def_domain, only      : nelem,ltype,nnode,ltypb,nboun,lnodb,&
!@                                ngaus,llapl,lorde,ltopo,elmar,lnods
!@    use def_chemic, only      : kfl_transport_chm,hk_gp,coeff_cp_k,W_k,dummy_enthalpy, &
!@                                ittot_chm, kfl_model_chm
!@
!@#ifdef CANTERA
!@    use def_chemic, only      : gas_chm
!@#endif
!@    use mod_physics, only     : physics_T_2_HCp
!@
!@    implicit none
!@
!@    integer(ip) :: iboun,pblty,inodb,pnodb,ipoin
!@    integer(ip) :: ielem,igaus,inode,iclas,ivalu
!@    integer(ip) :: pelty,pnode,lnods_loc(mnode)
!@    integer(ip) :: pgaus,plapl,porde,ptopo
!@    real(rp)    :: elcon(mnode,nclas_chm)
!@    real(rp)    :: elcod(ndime,mnode)
!@    real(rp)    :: elh(mnode)
!@    real(rp)    :: eltem(mnode)
!@    real(rp)    :: gpcon(mgaus,nclas_chm)
!@    real(rp)    :: gph(mgaus)
!@    real(rp)    :: gptem(mgaus)
!@    real(rp)    :: dummr
!@    real(rp)    :: cploc(6,2)
!@    real(rp)    :: aux_h
!@
!@    if (INOTMASTER) then
!@       !
!@       ! Loop over elements
!@       !
!@
!@       elements: do ielem = 1,nelem
!@          !
!@          ! Element dimensions
!@          !
!@          pelty = ltype(ielem)
!@
!@          if( pelty > 0 ) then
!@              pnode = nnode(pelty)
!@              pgaus = ngaus(pelty)
!@              plapl = llapl(pelty)
!@              porde = lorde(pelty)
!@              ptopo = ltopo(pelty)
!@
!@              !
!@              ! Gather all
!@              !
!@              lnods_loc(1:pnode) = lnods(1:pnode,ielem)
!@              call chm_gatherProp_finiteRate( &
!@                       pnode,lnods_loc,elcod,elcon(1:pnode,1:nclas_chm),elh,eltem)
!@
!@              !
!@              ! Initialization variables
!@              !
!@              gpcon = 0.0_rp
!@              gph   = 0.0_rp
!@              gptem = 0.0_rp
!@              condu_gp(ielem) % a               = 0.0_rp
!@              hk_gp(ielem) % a                  = 0.0_rp
!@              visco_gp(ielem) % a               = 0.0_rp
!@              wmean_gp(ielem) % a (1:pgaus,1,1) = 0.0_rp
!@              sphec_gp_ht(ielem) % a            = -666.0_rp
!@              sphec_gp_lt(ielem) % a            = -777.0_rp
!@              sphec_gp(ielem) % a(1:pgaus,1,1)  = 0.0_rp
!@
!@              !
!@              ! Species mass fraction Y_k at gauss points
!@              !
!@              do iclas = 1,nclas_chm
!@                 do igaus = 1,pgaus
!@                    do inode = 1,pnode
!@                       gpcon(igaus,iclas) = gpcon(igaus,iclas)&
!@                                            + elmar(pelty)%shape(inode,igaus) * elcon(inode,iclas)
!@                    end do
!@                 end do
!@              end do
!@
!@              !
!@              ! Enthalpy at gauss points
!@              !
!@              do igaus = 1,pgaus
!@                 do inode = 1,pnode
!@                    gph(igaus)   = gph(igaus) &
!@                                         + elmar(pelty)%shape(inode,igaus) * elh(inode)
!@                    gptem(igaus) = gptem(igaus) &
!@                                         + elmar(pelty)%shape(inode,igaus) * eltem(inode)
!@                 end do
!@              end do
!@
!@              !
!@              ! Compute transport properties
!@              !
!@              if (kfl_transport_chm == 1) then
!@
!@                 !
!@                 ! Cantera properties
!@                 !
!@#ifdef CANTERA
!@                 do igaus = 1,pgaus
!@                    call setState_TPX(gas_chm,gptem(igaus),prthe(1),gpcon(igaus,1:nspec_chm))
!@                    call setMassFractions(gas_chm,gpcon(igaus,1:nspec_chm))
!@
!@                    call cantera_alya_cp(nspec_chm,coeff_cp_k,gpcon(igaus,1:nspec_chm), &
!@                                      W_k,sphec_gp_lt(ielem) % a(igaus,:,1), sphec_gp_ht(ielem) % a(igaus,:,1))
!@
!@                    if (kfl_model_chm == 4) then  ! Compute cp for the mixture for CMC model
!@                       do ivalu = 1,6
!@                          cploc(ivalu,1) = sphec_gp_lt(ielem) % a(igaus,ivalu,1)
!@                          cploc(ivalu,2) = sphec_gp_ht(ielem) % a(igaus,ivalu,1)
!@                        end do
!@                        call physics_T_2_HCp(gptem(igaus), cploc, aux_h, sphec_gp(ielem) % a(igaus,1,1))
!@                    end if
!@
!@                    wmean_gp(ielem) % a(igaus,1,1)   = meanmolecularweight(gas_chm) / 1000.0_rp
!@                    visco_gp(ielem) % a(igaus,1,1)   = viscosity(gas_chm)
!@                    condu_gp(ielem) % a(igaus,1,1)   = thermalConductivity(gas_chm)
!@                    !
!@                    ! h_k = (h_k/RT)_cantera * Ru/W*T
!@                    !
!@                    call setState_TPX(gas_chm,298.15_rp,prthe(1),gpcon(igaus,1:nspec_chm))
!@                    call setMassFractions(gas_chm,gpcon(igaus,1:nspec_chm))
!@                    call getenthalpies_rt(gas_chm,hk_gp (ielem) % a (1:nspec_chm,igaus,1) )
!@                    hk_gp (ielem) % a (1:nspec_chm,igaus,1) = hk_gp (ielem) % a (1:nspec_chm,igaus,1)&
!@                        * gasco / W_k(1:nspec_chm) * 298.15_rp
!@
!@                 end do
!@#endif
!@              elseif (kfl_transport_chm == 2) then
!@                 !
!@                 ! User-defined properties
!@                 !
!@                 call runend('In chm_getProp_finiteRate: User-defined properties not implemented yet.')
!@              else
!@                 !
!@                 ! Unrecognised
!@                 !
!@                 call runend('In chm_getProp_finiteRate: kfl_transport_chm is not recognised.')
!@              endif
!@
!@          end if
!@       end do elements
!@
!@       !
!@       ! Cp coefficients on boundary for BC's
!@       !
!@
!@#ifdef CANTERA
!@       boundaries: do iboun = 1,nboun
!@          pblty = ltypb(iboun)
!@          pnodb = nnode(pblty)
!@          do inodb = 1,pnodb
!@             ipoin = lnodb(inodb,iboun)
!@
!@
!@             call setState_TPX(gas_chm, max(200.0_rp,min(3000.0_rp,tempe(ipoin,1))), &
!@                               prthe(1), conce(ipoin,1:nspec_chm,1))
!@             call setMassFractions(gas_chm,conce(ipoin,1:nspec_chm,1))
!@
!@             call cantera_alya_cp(nspec_chm,coeff_cp_k,conce(ipoin,1:nspec_chm,1), &
!@                                      W_k,sphec(ipoin,1:6,1),sphec(ipoin,1:6,2))
!@          end do
!@       end do boundaries
!@#endif
!@
!@    end if
!@
!@  end subroutine chm_getProp_finiteRate
!@
!@  subroutine chm_gatherProp_finiteRate( &
!@                 pnode,lnods,elcod,elcon,elh,eltem)
!@     !------------------------------------------------------------------------
!@     !****f* Chemic/mod_chm_element_operations/chm_gatherProp_finiteRate
!@     ! NAME
!@     !    chm_gatherProp_finiteRate
!@     ! DESCRIPTION
!@     !    Gather operations for the combustion models
!@     ! USES
!@     ! USED BY
!@     !    chm_getProp_finiteRate
!@     !***
!@     !------------------------------------------------------------------------
!@     use def_kintyp, only     :  ip,rp
!@     use def_master, only     :  conce,therm,tempe
!@
!@     implicit none
!@     integer(ip), intent(in)  :: pnode
!@     integer(ip), intent(in)  :: lnods(pnode)
!@     real(rp),    intent(out) :: elcod(ndime,pnode)
!@     real(rp),    intent(out) :: elcon(pnode,nclas_chm)
!@     real(rp),    intent(out) :: elh(pnode)
!@     real(rp),    intent(out) :: eltem(pnode)
!@
!@     integer(ip)              :: inode,ipoin,iclas,idime
!@
!@     !
!@     ! Initialization
!@     !
!@     elh   = 0.0_rp
!@     eltem = 0.0_rp
!@     elcod = 0.0_rp
!@     elcon = 0.0_rp
!@
!@     !
!@     ! Concentration and coordinates
!@     !
!@     do inode=1,pnode
!@        ipoin=lnods(inode)
!@        do iclas=1,nclas_chm
!@           elcon(inode,iclas) = conce(ipoin,iclas,1)
!@        end do
!@
!@        do idime=1,ndime
!@           elcod(idime,inode)   = coord(idime,ipoin)
!@        end do
!@
!@        elh(inode)   = therm(ipoin,1)
!@        eltem(inode) = max(200.0_rp, min(3000.0_rp, tempe(ipoin,1)))
!@     end do
!@
!@  end subroutine chm_gatherProp_finiteRate
!@
!@  subroutine chm_calc_enthalpy_transport_finiteRate(ielem,pgaus,gpdif,gpdit,gphco, &
!@                     grad_Yk_hk_gp,grad_h_gp,hk_grad_Yk_gp,grad_T_gp)
!@    !-----------------------------------------------------------------------
!@    !****f* chemic/chm_calc_enthalpy_transport_finiteRate
!@    ! NAME
!@    !    chm_calc_enthalpy_transport_finiteRate
!@    ! DESCRIPTION
!@    !    Compute enthalpy transport term to be used as RHS in temper
!@    ! USES
!@    ! USED BY
!@    !-----------------------------------------------------------------------
!@    use def_kintyp, only      : ip,rp
!@    use def_master, only      : enthalpy_transport, kfl_htran
!@
!@    implicit none
!@    integer(ip), intent(in) :: ielem
!@    integer(ip), intent(in) :: pgaus
!@    real(rp), intent(in)    :: gpdif(pgaus,nclas_chm)
!@    real(rp), intent(in)    :: gpdit(pgaus)
!@    real(rp), intent(in)    :: gphco(pgaus)
!@    real(rp), intent(in)    :: grad_Yk_hk_gp(pgaus,nspec_chm,ndime)
!@    real(rp), intent(in)    :: grad_h_gp(pgaus,ndime)
!@    real(rp), intent(in)    :: hk_grad_Yk_gp(pgaus,nspec_chm,ndime)
!@    real(rp), intent(in)    :: grad_T_gp(pgaus,ndime)
!@
!@    integer(ip) :: igaus,iclas
!@
!@    if (kfl_htran == 0) then
!@       !
!@       ! Calculate transport of enthalpy at gauss in the simplified way
!@       !
!@       do igaus=1,pgaus
!@          enthalpy_transport(ielem)%a(igaus,1:ndime,1) = 0.0_rp
!@          do iclas=1,nspec_chm
!@             enthalpy_transport(ielem)%a(igaus,1:ndime,1) = enthalpy_transport(ielem)%a(igaus,1:ndime,1) + gpdif(igaus,iclas)&
!@                 * grad_Yk_hk_gp(igaus,iclas,1:ndime)
!@          end do
!@          enthalpy_transport(ielem)%a(igaus,1:ndime,1) = enthalpy_transport(ielem)%a(igaus,1:ndime,1) + gpdit(igaus)&
!@              * grad_h_gp(igaus,1:ndime)
!@
!@       end do
!@
!@    else
!@       !
!@       ! Calculate transport of enthalpy at gauss in the detailed way:
!@       ! - rho * sum ( h_k * Y_k * V_k ) + k * grad(T) = rho * sum ( D_k * h_k * grad (Y_k) ) + k * grad(T)
!@       ! where D_k is the laminar mass difusivity
!@       ! Then, add the contribution from Fourier law
!@       !
!@       do igaus=1,pgaus
!@          enthalpy_transport(ielem)%a(igaus,1:ndime,1) = 0.0_rp
!@          do iclas=1,nspec_chm
!@             enthalpy_transport(ielem)%a(igaus,1:ndime,1) = enthalpy_transport(ielem)%a(igaus,1:ndime,1) + &
!@                                                            gpdif(igaus,iclas) * hk_grad_Yk_gp(igaus,iclas,1:ndime)
!@          end do
!@          enthalpy_transport(ielem)%a(igaus,1:ndime,1) = enthalpy_transport(ielem)%a(igaus,1:ndime,1) + gphco(igaus)&
!@              * grad_T_gp(igaus,1:ndime)
!@
!@       end do
!@
!@    end if
!@
!@  end  subroutine chm_calc_enthalpy_transport_finiteRate
!@
!@
!@  subroutine chm_calc_div_enthalpy_transport_finiteRate
!@    !-----------------------------------------------------------------------
!@    !****f* chemic/chm_calc_div_enthalpy_transport_finiteRate
!@    ! NAME
!@    !    chm_calc_div_enthalpy_transport_finiteRate
!@    ! DESCRIPTION
!@    !    Compute divergence of enthalpy transport by diffusion
!@    ! USES
!@    ! USED BY
!@    !-----------------------------------------------------------------------
!@    use def_kintyp, only      : ip,rp,r2p
!@    use def_master, only      : div_enthalpy_transport
!@    use def_master, only      : enthalpy_transport
!@
!@    use def_domain, only      : nelem,nnode,ntens,&
!@                                ngaus,llapl,lorde,ltopo,elmar,lnods,&
!@                                ltype
!@    use def_chemic, only      : enthalpy_transport_nodes
!@
!@    implicit none
!@    integer(ip) :: ielem,igaus,inode,ipoin,idime
!@    integer(ip) :: pelty,pnode
!@    integer(ip) :: pgaus,plapl,porde,ptopo
!@    real(rp)    :: el_h_transport(ndime,mnode)
!@    real(rp)    :: elcod(ndime,mnode)
!@    real(rp)    :: gpcar(ndime,mnode,mgaus)                 ! dNk/dxj
!@    real(rp)    :: gphes(ntens,mnode,mgaus)                 ! dNk/dxidxj
!@    real(rp)    :: gpvol(mgaus)
!@    type(r2p),pointer :: aux_r2p(:)
!@
!@    !
!@    ! Initialization
!@    !
!@    enthalpy_transport_nodes = 0.0_rp
!@
!@    !
!@    ! Smooth enthalpy transport to nodes
!@    !
!@    nullify(aux_r2p)
!@    call memory_alloca(mem_modul(1:2,modul),'AUX_R2P','chm_calc_div_enthalpy_transport_finiteRate',aux_r2p,max(1_ip,nelem))
!@    do ielem = 1,nelem
!@       pelty = ltype(ielem)
!@       pgaus = ngaus(pelty)
!@       call memory_alloca(mem_modul(1:2,modul),'AUX_R2P % A','chm_calc_div_enthalpy_transport_finiteRate',aux_r2p(ielem)%a,pgaus,&
!@           ndime)
!@       aux_r2p(ielem) % a = enthalpy_transport(ielem) % a(:,1:ndime,1)
!@    end do
!@    call smoot5(aux_r2p,enthalpy_transport_nodes,ndime)
!@    call memory_deallo(mem_modul(1:2,modul),'AUX_R2P','chm_calc_div_enthalpy_transport_finiteRate',aux_r2p)
!@
!@
!@    if (INOTMASTER) then
!@       !
!@       ! Loop over elements
!@       !
!@       elements: do ielem = 1,nelem
!@          !
!@          ! Element dimensions
!@          !
!@          pelty = ltype(ielem)
!@
!@          if( pelty > 0 ) then
!@              pnode = nnode(pelty)
!@              pgaus = ngaus(pelty)
!@              plapl = llapl(pelty)
!@              porde = lorde(pelty)
!@              ptopo = ltopo(pelty)
!@
!@              !
!@              ! Initialization
!@              !
!@              gpvol = 0.0_rp
!@              gphes = 0.0_rp
!@
!@              !
!@              ! Gather for concentration and coordinates
!@              !
!@              do inode=1,pnode
!@                 ipoin=lnods(inode,ielem)
!@
!@                 el_h_transport(1:ndime,inode) = enthalpy_transport_nodes(1:ndime,ipoin)
!@                 elcod(1:ndime,inode)          = coord(1:ndime,ipoin)
!@              end do
!@
!@              !
!@              ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
!@              !
!@              call elmcar(&
!@                   pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
!@                   elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
!@                   gphes,ielem)
!@
!@              !
!@              ! Divergence of enthalpy transport at gauss points
!@              !
!@              do igaus = 1,pgaus
!@                 div_enthalpy_transport(ielem) % a (igaus,1,1) = 0.0_rp
!@                 do inode = 1,pnode
!@                    do idime = 1,ndime
!@                       div_enthalpy_transport(ielem) % a (igaus,1,1) = div_enthalpy_transport(ielem) % a (igaus,1,1) + &
!@                                              gpcar(idime,inode,igaus) * el_h_transport(idime,inode) !* gpvol(igaus)
!@                    end do
!@                 end do
!@              end do
!@
!@          end if
!@       end do elements
!@
!@    end if
!@
!@  end subroutine chm_calc_div_enthalpy_transport_finiteRate
!@
!@  subroutine get_rates_pfa(temp,press,y_k,gen,r_mat)
!@    !-----------------------------------------------------------------------
!@    !****f* chemic/get_rates_pfa
!@    ! NAME
!@    !    get_rates_pfa
!@    ! DESCRIPTION
!@    !    Computes the Realtion matrix by PFA(A path flux analysis method for the
!@    !    reduction of detailed chemical kinetic mechanisms by Wenting Sun et.al)
!@    ! USES Mass Fractions, Temperature, Pressure and genration of reduction (1
!@    ! or 2)
!@    ! USED BY called by chm_IntegrateSource_finiteRate
!@    !-----------------------------------------------------------------------
!@
!@      use def_chemic,                 only : nreac_chm
!@      use def_master,                 only : zeror
!@
!@#ifdef CANTERA
!@      use def_chemic,                 only : gas_chm
!@#endif
!@      implicit none
!@      real(rp)                          :: y_k(nspec_chm)
!@      real(rp)                          :: temp,press
!@      integer(ip), intent(in)           :: gen
!@      integer                           :: i , j , k, m
!@      double precision                  :: div, addp, addc, dummp,dummc
!@      double precision                  :: ca(nspec_chm), pa(nspec_chm),pab(nspec_chm,nspec_chm), &
!@                                           cab(nspec_chm,nspec_chm), rp1(nspec_chm,nspec_chm),    &
!@                                           rp2(nspec_chm,nspec_chm), rc1(nspec_chm,nspec_chm),    &
!@                                           rc2(nspec_chm,nspec_chm), rp3(nspec_chm,nspec_chm),    &
!@                                           rc3(nspec_chm,nspec_chm), react(nspec_chm,nreac_chm),  &
!@                                           prod(nspec_chm,nreac_chm),nrx(nreac_chm)
!@      double precision, intent(inout)   :: r_mat(nspec_chm,nspec_chm)
!@
!@
!@! Initialize Some Stuff
!@
!@      addp  = 0.0_rp
!@      addc  = 0.0_rp
!@      nrx   = 0.0_rp
!@      pa    = 0.0_rp
!@      ca    = 0.0_rp
!@      pab   = 0.0_rp
!@      cab   = 0.0_rp
!@      rp1   = 0.0_rp
!@      rp2   = 0.0_rp
!@      rc1   = 0.0_rp
!@      rc2   = 0.0_rp
!@      rc1   = 0.0_rp
!@      rc2   = 0.0_rp
!@      r_mat = 0.0_rp
!@
!@#ifdef CANTERA
!@      do j=1,nreac_chm
!@        do i=1,nspec_chm
!@            react(i,j)=  reactantstoichcoeff(gas_chm,i,j)
!@            prod(i,j)=   productstoichcoeff(gas_chm,i,j)
!@        end do
!@      end do
!@#endif
!@
!@      !
!@      ! Set gas State to get production rates of species
!@      !
!@#ifdef CANTERA
!@      call setState_TPX(gas_chm,temp,press,y_k)
!@      call setMassFractions(gas_chm,y_k)
!@      !
!@      ! Get Production Rates of Species
!@      !
!@      call getnetratesofprogress(gas_chm,nrx)
!@#endif
!@      do i=1,nspec_chm
!@        do k = 1,nreac_chm
!@            if (nrx(k) > 0.0_rp) then
!@                pa(i) = pa(i) + abs(prod(i,k)*nrx(k))
!@                ca(i) = ca(i) + abs(react(i,k)*nrx(k))
!@            end if
!@            if (nrx(k) < 0.0_rp) then
!@                pa(i) = pa(i) + abs(react(i,k)*nrx(k))
!@                ca(i) = ca(i) + abs(prod(i,k)*nrx(k))
!@            end if
!@        end do
!@      end do
!@
!@      do i=1,nspec_chm
!@          do j=1,nspec_chm
!@              if (i /= j) then
!@                  pab(i,j) = 0.0_rp
!@                  cab(i,j) = 0.0_rp
!@                  do k=1,nreac_chm
!@                      if((prod(j,k) /= 0.0_rp) .or. (react(j,k) /= 0.0_rp)) then
!@                          addp = 0.0_rp
!@                          addc = 0.0_rp
!@                          if (nrx(k) > 0.0_rp)  then
!@                              addp =  abs(prod(i,k)*nrx(k))
!@                              pab(i,j) = pab(i,j) + addp
!@                              addc =  abs(react(i,k)*nrx(k))
!@                              cab(i,j) = cab(i,j) + addc
!@                          else if (nrx(k) < 0.0_rp) then
!@                              addp = abs(react(i,k)*nrx(k))
!@                              pab(i,j) = pab(i,j) + addp
!@                              addc = abs(prod(i,k)*nrx(k))
!@                              cab(i,j) = cab(i,j) + addc
!@                          end if
!@                      end if
!@                  end do
!@              end if
!@          end do
!@      end do
!@
!@      do i=1,nspec_chm
!@          div = max(pa(i),ca(i))
!@          do j=1,nspec_chm
!@              if (div == 0.0_rp) then
!@                 rp1(i,j) = 0.0_rp
!@                 rc1(i,j) = 0.0_rp
!@              else
!@#ifdef CANTERA
!@                 if (i /= j) then
!@                    call cantera_divide(pab(i,j),div,rp1(i,j))
!@                    call cantera_divide(cab(i,j),div,rc1(i,j))
!@                 end if
!@#endif
!@              end if
!@          end do
!@      end do
!@
!@      do i=1,nspec_chm
!@          do j=1,nspec_chm
!@              if (i /= j) then
!@                  do k = 1,nspec_chm
!@                      if ((k /= i ) .and. (k /= j)) then
!@                          rp2(i,j) = rp2(i,j) + (rp1(i,k) * rp1(k,j))
!@                          rc2(i,j) = rc2(i,j) + (rc1(i,k) * rc1(k,j))
!@                      end if
!@                  end do
!@              end if
!@          end do
!@      end do
!@
!@      !
!@      ! If Second Generation path fluxes are to be included
!@      !
!@
!@      if (gen == 2) then
!@          do j=1,nspec_chm
!@            do i=1,nspec_chm
!@                if (i /= j) then
!@                    do k = 1,nspec_chm
!@                        do m =1,nspec_chm
!@                            if ((k /= i ) .and. (k /= j) .and. (k /= m) .and. (m /= i) .and. (m /= j)) then
!@                                rp3(i,j) = rp3(i,j) + (rp1(i,m) * rp1(m,k) * rp1(k,j))
!@                                rc3(i,j) = rc3(i,j) + (rc1(i,m) * rc1(m,k) * rc1(k,j))
!@                            end if
!@                        end do
!@                    end do
!@                end if
!@            end do
!@          end do
!@      end if
!@
!@
!@      if (gen == 2) then
!@          do j=1,nspec_chm
!@            do i=1,nspec_chm
!@                if (i /= j) then
!@                    r_mat(i,j) = rp1(i,j) + rp2(i,j) + rc1(i,j) + rc2(i,j) + rp3(i,j)+ rc3(i,j)
!@                end if
!@            end do
!@          end do
!@      else
!@          do i=1,nspec_chm
!@            do j=1,nspec_chm
!@                if (i /= j) then
!@                    r_mat(i,j) = rp1(i,j) + rp2(i,j) + rc1(i,j) + rc2(i,j)
!@                end if
!@            end do
!@          end do
!@      end if
!@
!@  end subroutine get_rates_pfa
!@
!@
!@  subroutine get_reduced_reaction_rates(r_mat,key_spec,react)
!@    !-----------------------------------------------------------------------
!@    !****f* chemic/get_reduced_reaction_rates
!@    ! NAME
!@    !    get_reduced_reaction_rates
!@    ! DESCRIPTION
!@    !    Uses the Realtion matrix obtained by get_rates_pfa to generate an array
!@    !    of 1 and 0's which are the reation coeffcient multipliers
!@    ! USES Rmat, crit
!@    ! USED BY called by chm_IntegrateSource_finiteRate after get_rates_pfa
!@    !-----------------------------------------------------------------------
!@
!@      use def_chemic,                 only : kfl_key_chm
!@      use def_chemic,                 only : dac_crit_chm
!@      use def_chemic,                 only : nreac_chm
!@      use def_master,                 only : zeror
!@
!@#ifdef CANTERA
!@      use def_chemic,                 only : gas_chm
!@#endif
!@
!@      implicit none
!@
!@      integer                                     :: i,j
!@      double precision, intent(inout)             :: r_mat(nspec_chm,nspec_chm)
!@      double precision                            :: r_mat_temp(nspec_chm,nspec_chm)
!@      double precision                            :: crit, maxv, div
!@      integer                                     :: spec(nspec_chm)
!@      integer,intent(inout)                       :: react(nreac_chm)
!@      integer(ip), intent(in)                     :: key_spec(kfl_key_chm)
!@      spec       = 1_ip
!@      react      = 1_ip
!@      r_mat_temp = 0.0_rp
!@
!@      do i=1,size(key_spec)
!@         r_mat_temp(key_spec(i),:) = r_mat(key_spec(i),:)
!@      end do
!@
!@      do i=1,size(key_spec)
!@         do j=1,nspec_chm
!@#ifdef CANTERA
!@            call cantera_divide(r_mat_temp(i,j),zeror,div)
!@            if (div /= 0.0_rp) then
!@                call cantera_log(div)
!@            else
!@                div = 0.0_rp
!@            end if
!@#endif
!@            r_mat_temp(i,j) = div
!@         end do
!@      end do
!@      maxv = maxval(r_mat_temp)
!@
!@      do i=1,size(key_spec)
!@         do j=1,nspec_chm
!@            if (maxv < 1e-16_rp) then
!@               r_mat_temp(i,j) = 0.0_rp
!@            else
!@               r_mat_temp(i,j) = r_mat_temp(i,j)/maxv
!@            end if
!@         end do
!@      end do
!@
!@      do i=1,nspec_chm
!@         if (maxval(r_mat_temp(:,i)) < (dac_crit_chm)) then
!@            spec(i) = 0_ip
!@         end if
!@      end do
!@
!@      do i=1,size(key_spec)
!@         spec(key_spec(i)) = 1_ip
!@      end do
!@
!@      do i=1,nspec_chm
!@         if (spec(i) <  1_ip) then
!@            do j=1,nreac_chm
!@#ifdef CANTERA
!@               if ((reactantstoichcoeff(gas_chm,i,j) > 0.001_rp) .or. (productstoichcoeff(gas_chm,i,j) > 0.001_rp)) then
!@                  react(j) = 0_ip
!@               end if
!@#endif
!@            end do
!@         end if
!@      end do
!@
!@  end subroutine get_reduced_reaction_rates
!@
!@  subroutine get_prog_var(y_k,prog_var)
!@    !-----------------------------------------------------------------------
!@    !****f* chemic/get_prog_var
!@    ! NAME
!@    !    get_prog_var
!@    ! DESCRIPTION
!@    ! Get Progress Varible for TDAC implementation
!@    ! USES y_k
!@    ! USED BY called by chm_IntegrateSource_finiteRate
!@    !-----------------------------------------------------------------------
!@
!@#ifdef CANTERA
!@      use def_chemic, only       : gas_chm
!@#endif
!@
!@      implicit none
!@
!@      integer                                     :: i
!@      double precision, intent(in)                :: y_k(nspec_chm)
!@      double precision, intent(inout)             :: prog_var
!@      character (LEN=100)                         :: key_species(4)
!@      character(len=:), allocatable               :: key_species_trimmed
!@      integer                                     :: index_key(4)
!@
!@      key_species(1) = 'CO2'
!@      key_species(2) = 'H2O'
!@      key_species(3) = 'CO'
!@      key_species(4) = 'H2'
!@
!@#ifdef CANTERA
!@      do i=1,4
!@         key_species_trimmed = trim(adjustl(key_species(i)))
!@         index_key(i) = speciesIndex(gas_chm, key_species_trimmed)
!@      end do
!@#endif
!@      !prog_var = y_k(index_key(1)) * 4.0/44.0 + y_k(index_key(2)) * 2.0/18.0 +  y_k(index_key(3)) * 1/28.0 + y_k(index_key(4))&
!@           * 0.5/2.0
!@      !prog_var = y_k(index_key(1)) * 0.125 + y_k(index_key(2)) * 0.5 +  y_k(index_key(3)) * 0.125 + y_k(index_key(4)) * 0.25
!@      prog_var = y_k(index_key(1)) * 0.1_rp + y_k(index_key(2)) * 0.4_rp +  y_k(index_key(3)) * 0.10_rp + y_k(index_key(4))&
!@          * 0.2_rp
!@  end subroutine  get_prog_var
!@
!@  subroutine get_correlation(y_k_old,y_k_new,tempe_old,tempe_new,corr)
!@    !-----------------------------------------------------------------------
!@    !****f* chemic/get_correlation
!@    ! NAME
!@    !    get_correlation
!@    ! DESCRIPTION
!@    ! Get Correlation for CO-DAC Implementation
!@    ! USES y_k_old,y_k_new,tempe_old,tempe_new,error
!@    ! USED BY called by chm_IntegrateSource_finiteRate
!@    !-----------------------------------------------------------------------
!@
!@#ifdef CANTERA
!@      use def_chemic, only       : gas_chm
!@#endif
!@      use def_master, only       : zeror
!@      use def_chemic, only       : dac_cor_chm
!@
!@      implicit none
!@
!@      integer                                     :: i,j
!@      double precision, intent(in)                :: y_k_old(nspec_chm), y_k_new(nspec_chm)
!@      double precision                            :: error(5)
!@      double precision, intent(in)                :: tempe_old, tempe_new
!@      character (LEN=100)                         :: key_species(4)
!@      character(len=:), allocatable               :: key_species_trimmed
!@      integer                                     :: index_key(4)
!@      integer                                     :: corr
!@      double precision                            :: corr_coeff(5)
!@      double precision                            :: max_corr_coeff, num, den
!@
!@      key_species(1) = 'CH4'
!@      key_species(2) = 'OH'
!@      key_species(3) = 'CH2O'
!@      key_species(4) = 'HO2'
!@
!@      do i=1,5
!@         error(i) = dac_cor_chm
!@      end do
!@#ifdef CANTERA
!@      do i=1,4
!@         key_species_trimmed = trim(adjustl(key_species(i)))
!@         index_key(i) = speciesIndex(gas_chm, key_species_trimmed)
!@      end do
!@#endif
!@      !corr_coeff(1) = abs(tempe_new - tempe_old + zeror)/error(1)
!@      corr_coeff(1) = abs(tempe_new - tempe_old + zeror)/0.05_rp
!@      !num = tempe_new - tempe_old + zeror
!@      !den = error(1)
!@      !print *, num, den
!@      !call cantera_divide(num,corr_coeff(1),den)
!@
!@      do i =1,4
!@         j = 2
!@         if ((y_k_new(index_key(i)) - y_k_old(index_key(i))) < 1e-16_rp) then
!@             corr_coeff(j) = 0.0_rp
!@         else
!@             !corr_coeff(j) = (abs(log(y_k_new(index_key(i))) - log(y_k_old(index_key(i)))) + zeror)/error(j)
!@             corr_coeff(j) = (abs(log(y_k_new(index_key(i))) - log(y_k_old(index_key(i)))) + zeror)/0.2_rp
!@         end if
!@        j = j+1
!@      end do
!@
!@      max_corr_coeff = maxval(corr_coeff)
!@      if(max_corr_coeff < 1.0_rp) then
!@          corr = 1_ip
!@      else
!@          corr = 0_ip
!@      end if
!@
!@  end subroutine  get_correlation
!@
!@  subroutine get_mixture_fraction(y_k,z)
!@    !-----------------------------------------------------------------------
!@    !****f* chemic/get_correlation
!@    ! NAME
!@    !    get_mixture_fraction
!@    ! DESCRIPTION
!@    ! Get Mixture fraction to be used by TDAC model
!@    ! USES y_k, b_f_chm and b_o_chm (from chm_inivar)
!@    ! USED BY called by chm_IntegrateSource_finiteRate
!@    !-----------------------------------------------------------------------
!@
!@#ifdef CANTERA
!@      use def_chemic, only       : gas_chm
!@#endif
!@      use def_chemic,            only : bf_fuel_chm,bo_oxy_chm
!@
!@      implicit none
!@
!@      integer                                     :: i,j
!@      double precision                            :: dummi_h, dummi_c, dummi_o, b_chm
!@      double precision, intent(in)                :: y_k(nspec_chm)
!@      double precision, intent(inout)             :: z
!@
!@
!@#ifdef CANTERA
!@        !
!@        ! B_O Calculation (For Mixture fraction)
!@        !
!@        call cantera_elemh(y_k,dummi_h)
!@        call cantera_elemo(y_k,dummi_o)
!@        call cantera_elemc(y_k,dummi_c)
!@        b_chm = ( 2.0_rp* dummi_c/12.0107_rp )  + ( 0.5_rp * dummi_h/1.00794 ) - dummi_o/15.999_rp
!@#endif
!@        z = (b_chm - bo_oxy_chm) / (bf_fuel_chm - bo_oxy_chm)
!@
!@  end subroutine  get_mixture_fraction
!@
!@
!@  subroutine chm_IntegrateSource_finiteRate(dt)
!@
!@    use def_domain,                 only : ndime
!@    use def_domain,                 only : npoin
!@    use def_master,                 only : conce,tempe,prthe,dtinv
!@    use def_master,                 only : therm
!@    use def_master,                 only : zeror
!@    use def_master,                 only : mitim
!@    use def_chemic,                 only : nreac_chm
!@    use def_chemic,                 only : mechanism_path
!@    use def_chemic,                 only : Red_spec
!@    use def_chemic,                 only : kfl_pfa_chm
!@    use def_chemic,                 only : Y_k_n
!@    use def_chemic,                 only : React_ind
!@    use def_chemic,                 only : Corr_chm
!@    use def_chemic,                 only : nsize_red
!@    use def_chemic,                 only : kfl_key_chm
!@    use def_chemic,                 only : dac_crit_chm
!@    use def_chemic,                 only : mixfr_chm
!@    use def_chemic,                 only : prog_var_chm
!@    use def_chemic,                 only : kfl_freq_chm
!@    use def_chemic,                 only : kfl_z_chm
!@    use def_chemic,                 only : ittot_chm
!@    use def_chemic,                 only : src_chm
!@    use def_chemic,                 only : table_fw
!@    use def_chemic,                 only : kfl_tdac_write_chm
!@    use def_chemic,                 only : table_coords
!@    use def_chemic,                 only : kfl_cont_chm
!@    use mod_interp_tab,             only : tab_interp
!@    use mod_interp_tab,             only : tab_load_file
!@    use mod_interp_tab,             only : i2inds
!@    use mod_gradie
!@    use mod_interp_tab,             only : tab_par_exchange
!@#ifdef CANTERA
!@    use def_chemic,                 only : gas_chm
!@#endif
!@    use mod_ker_proper
!@    implicit none
!@    real(rp), intent(in)                :: dt
!@
!@    integer(ip)                         :: ipoin ,gen, corr, dummi, irow, ii, jj, nrow
!@    integer(ip)                         :: shap(6), inds(6)
!@    real(rp)                            :: test_output(6)
!@    real(rp), pointer                   :: Y_k_0(:,:)
!@    real(rp), pointer                   :: T_0(:)
!@    integer(ip)                         :: all_spec_index(nspec_chm)
!@    double precision                    :: r_mat(nspec_chm,nspec_chm)
!@    double precision                    :: temp_lim
!@    integer(ip)                         :: key_spec(kfl_key_chm)
!@    real(rp), pointer                   :: tab_scale_control(:)
!@    real(rp), pointer                   :: retva(:)
!@    real(rp), pointer                   :: prope_tmp(:)
!@
!@    if (INOTMASTER) then
!@
!@       nullify(Y_k_0)
!@       nullify(T_0)
!@
!@       call memory_alloca(mem_modul(1:2,modul),'Y_K_0','chm_IntegrateSource_finiteRate',Y_k_0,npoin,nspec_chm)
!@       call memory_alloca(mem_modul(1:2,modul),'T_0'  ,'chm_IntegrateSource_finiteRate',T_0,npoin)
!@
!@       Y_k_0 = 0.0_rp
!@       T_0   = 0.0_rp
!@
!@       all_spec_index = 1_ip
!@       gen            = 1_ip
!@       corr           = 1_ip
!@       key_spec       = 1_ip
!@       temp_lim       = 600.99_rp ! check this
!@       !
!@       ! Integration source term using CVODE from Cantera
!@       !
!@#ifdef CANTERA
!@
!@       !
!@       ! Complete mechanism
!@       !
!@       if ( kfl_pfa_chm == -1_ip ) then
!@
!@            do ipoin=1,npoin
!@               Y_k_0(ipoin,1:nspec_chm) = conce(ipoin,1:nspec_chm,1)
!@               T_0(ipoin)               = tempe(ipoin,1)
!@               if (T_0(ipoin) > temp_lim) call cantera_integrate(T_0(ipoin),prthe(1),Y_k_0(ipoin,:),dt)
!@            end do
!@       !
!@       ! Equilibrium in first timestep and reduced mechanims in the next
!@       !
!@       else if (kfl_pfa_chm == 0_ip) then
!@
!@           if (ittot_chm == 1_ip) then
!@               do ipoin=1,npoin
!@                  T_0(ipoin)                 = tempe(ipoin,1)
!@                  Y_k_n(ipoin,1:nspec_chm,1) = conce(ipoin,1:nspec_chm,1)
!@                  call cantera_equilibrate(T_0(ipoin),prthe(1),Y_k_n(ipoin,:,1),Red_spec,all_spec_index)
!@               end do
!@           endif
!@
!@           !
!@           ! Integration based on Key_Species with mass fractions of non_key kept
!@           ! constant
!@           !
!@           do ipoin=1,npoin
!@              T_0(ipoin)               = tempe(ipoin,1)
!@              Y_k_0(ipoin,1:nspec_chm) = conce(ipoin,1:nspec_chm,1)
!@              call cantera_equi_integrate(T_0(ipoin),prthe(1),Y_k_0(ipoin,:),Y_k_n(ipoin,:,1),dt,all_spec_index)
!@           end do
!@       !
!@       ! Static PFA
!@       !
!@       else if (kfl_pfa_chm == 1_ip) then
!@            if (ittot_chm == 100) then
!@                do ipoin=1,npoin
!@                   Y_k_0(ipoin,1:nspec_chm) = conce(ipoin,1:nspec_chm,1)
!@                   T_0(ipoin)               = tempe(ipoin,1)
!@                   call cantera_trim(Red_spec,nsize_red,key_spec)
!@                   if (T_0(ipoin) > temp_lim) then
!@                       call get_rates_pfa(T_0(ipoin),prthe(1),Y_k_0(ipoin,:),gen,r_mat)
!@                       call get_reduced_reaction_rates(r_mat,key_spec,React_ind(ipoin,:))
!@                       call cantera_dac_integrate(T_0(ipoin),prthe(1),Y_k_0(ipoin,:),React_ind(ipoin,:),dt)
!@                   else
!@                       React_ind(ipoin,:) = 0_ip
!@                       call cantera_dac_integrate(T_0(ipoin),prthe(1),Y_k_0(ipoin,:),React_ind(ipoin,:),dt)
!@                   end if
!@                end do
!@            else
!@                do ipoin=1,npoin
!@                   Y_k_0(ipoin,1:nspec_chm) = conce(ipoin,1:nspec_chm,1)
!@                   T_0(ipoin)               = tempe(ipoin,1)
!@                   if (T_0(ipoin) > temp_lim) then
!@                       call cantera_dac_integrate(T_0(ipoin),prthe(1),Y_k_0(ipoin,:),React_ind(ipoin,:),dt)
!@                   else
!@                       React_ind(ipoin,:) = 0_ip
!@                       call cantera_dac_integrate(T_0(ipoin),prthe(1),Y_k_0(ipoin,:),React_ind(ipoin,:),dt)
!@                   end if
!@                end do
!@            end if
!@
!@       !
!@       ! DAC
!@       !
!@       else if (kfl_pfa_chm == 2_ip) then
!@
!@            !
!@            ! DAC implementation based on Iterations
!@            !
!@          if (ittot_chm == 0) then
!@              do ipoin=1,npoin
!@                 Y_k_0(ipoin,1:nspec_chm) = conce(ipoin,1:nspec_chm,1)
!@                 T_0(ipoin)               = tempe(ipoin,1)
!@                 React_ind(ipoin,:)       = 1_ip
!@                 call cantera_dac_integrate(T_0(ipoin),prthe(1),Y_k_0(ipoin,:),React_ind(ipoin,:),dt)
!@              end do
!@          else
!@             if (mod(ittot_chm,kfl_freq_chm) == 0) then
!@                 do ipoin=1,npoin
!@                    Y_k_0(ipoin,1:nspec_chm) = conce(ipoin,1:nspec_chm,1)
!@                    T_0(ipoin)               = tempe(ipoin,1)
!@                    call cantera_trim(Red_spec,nsize_red,key_spec)
!@                    if (T_0(ipoin) > temp_lim) then
!@                        call get_rates_pfa(T_0(ipoin),prthe(1),Y_k_0(ipoin,:),gen,r_mat)
!@                        call get_reduced_reaction_rates(r_mat,key_spec,React_ind(ipoin,:))
!@                        call cantera_dac_integrate(T_0(ipoin),prthe(1),Y_k_0(ipoin,:),React_ind(ipoin,:),dt)
!@                    else
!@                        React_ind(ipoin,:) = 0_ip
!@                        call cantera_dac_integrate(T_0(ipoin),prthe(1),Y_k_0(ipoin,:),React_ind(ipoin,:),dt)
!@                    end if
!@                 end do
!@             else
!@                 do ipoin=1,npoin
!@                    Y_k_0(ipoin,1:nspec_chm) = conce(ipoin,1:nspec_chm,1)
!@                    T_0(ipoin)               = tempe(ipoin,1)
!@                    if (T_0(ipoin) > temp_lim) then
!@                        call cantera_dac_integrate(T_0(ipoin),prthe(1),Y_k_0(ipoin,:),React_ind(ipoin,:),dt)
!@                    else
!@                        React_ind(ipoin,:) = 0_ip
!@                        call cantera_dac_integrate(T_0(ipoin),prthe(1),Y_k_0(ipoin,:),React_ind(ipoin,:),dt)
!@                    end if
!@                 end do
!@             end if
!@          end if
!@       !
!@       ! CODAC
!@       !
!@       else if (kfl_pfa_chm == 3_ip) then
!@           if(ittot_chm == 1_ip) then
!@                do ipoin=1,npoin
!@                   Y_k_0(ipoin,1:nspec_chm) = conce(ipoin,1:nspec_chm,1)
!@                   T_0(ipoin)               = tempe(ipoin,1)
!@                   React_ind(ipoin,:) = 1_ip
!@                   call cantera_trim(Red_spec,nsize_red,key_spec)
!@                   call get_rates_pfa(T_0(ipoin),prthe(1),Y_k_0(ipoin,:),gen,r_mat)
!@                   call get_reduced_reaction_rates(r_mat,key_spec,React_ind(ipoin,:))
!@                   call cantera_dac_integrate(T_0(ipoin),prthe(1),Y_k_0(ipoin,:),React_ind(ipoin,:),dt)
!@                end do
!@           else
!@                do ipoin=1,npoin
!@                   Y_k_0(ipoin,1:nspec_chm) = conce(ipoin,1:nspec_chm,1)
!@                   T_0(ipoin)               = tempe(ipoin,1)
!@                   call get_correlation(conce(ipoin,1:nspec_chm,3),conce(ipoin,1:nspec_chm,1), &
!@                                        tempe(ipoin,3),tempe(ipoin,1),Corr_chm(ipoin))
!@                   if (Corr_chm(ipoin) < 1_ip) then
!@                           call cantera_trim(Red_spec,nsize_red,key_spec)
!@                           call get_rates_pfa(T_0(ipoin),prthe(1),Y_k_0(ipoin,:),gen,r_mat)
!@                           call get_reduced_reaction_rates(r_mat,key_spec,React_ind(ipoin,:))
!@                   end if
!@                   call cantera_dac_integrate(T_0(ipoin),prthe(1),Y_k_0(ipoin,:),React_ind(ipoin,:),dt)
!@                end do
!@           end if
!@       !
!@       ! TDAC
!@       !
!@       else if (kfl_pfa_chm == 4_ip) then
!@          if(ittot_chm == 1_ip) then
!@             do ipoin=1,npoin
!@                Y_k_0(ipoin,1:nspec_chm) = conce(ipoin,1:nspec_chm,1)
!@                T_0(ipoin)               = tempe(ipoin,1)
!@                React_ind(ipoin,:) = 1_ip
!@                call cantera_dac_integrate(T_0(ipoin),prthe(1),Y_k_0(ipoin,:),React_ind(ipoin,:),dt)
!@             end do
!@          else
!@             nullify(tab_scale_control)
!@             nullify(retva)
!@             call memory_alloca(mem_modul(1:2,modul),'TAB_SCALE_CONTROL'  ,'chm_IntegrateSource_finiteRate', tab_scale_control,&
!@                 table_fw % main_table % ndim )
!@             call memory_alloca(mem_modul(1:2,modul),'RETVA'              ,'chm_IntegrateSource_finiteRate', retva            ,&
!@                 table_fw % main_table % nvar )
!@
!@             do ipoin=1,npoin
!@                Y_k_0(ipoin,1:nspec_chm) = conce(ipoin,1:nspec_chm,1)
!@                T_0(ipoin)               = tempe(ipoin,1)
!@                !
!@                ! Get Mixture Fraction
!@                !
!@                call get_mixture_fraction(Y_k_0(ipoin,1:nspec_chm),mixfr_chm(ipoin))
!@                call get_prog_var(Y_k_0(ipoin,1:nspec_chm),prog_var_chm(ipoin))
!@                tab_scale_control(1) = mixfr_chm(ipoin)
!@                tab_scale_control(2) = prog_var_chm(ipoin)
!@                if (kfl_cont_chm > 2_ip) tab_scale_control(3) = T_0(ipoin)
!@                !
!@                ! Look up table
!@                !
!@                call tab_interp(table_fw % main_table, tab_scale_control, retva, inds, snap_on_nearest=.true., irow=irow)
!@                if (retva(1) /= 1.0_rp) then
!@                   React_ind(ipoin,:) = 1_ip
!@                   call cantera_trim(Red_spec,nsize_red,key_spec)
!@                   call get_rates_pfa(T_0(ipoin),prthe(1),Y_k_0(ipoin,:),gen,r_mat)
!@                   call get_reduced_reaction_rates(r_mat,key_spec,React_ind(ipoin,:))
!@                   table_fw % main_table % tab(irow,2:(nreac_chm+1)) = React_ind(ipoin,:)
!@                   table_fw % main_table % tab(irow,1) = 1.0_rp
!@                else
!@                   React_ind(ipoin,:) = table_fw % main_table % tab(irow,2:(nreac_chm+1))
!@                endif
!@                call cantera_dac_integrate(T_0(ipoin),prthe(1),Y_k_0(ipoin,:),React_ind(ipoin,:),dt)
!@             end do
!@
!@             call memory_deallo(mem_modul(1:2,modul),'TAB_SCALE_CONTROL'  ,'chm_IntegrateSource_finiteRate', tab_scale_control)
!@             call memory_deallo(mem_modul(1:2,modul),'RETVA'              ,'chm_IntegrateSource_finiteRate', retva            )
!@          end if
!@       else
!@           call runend('In chm_IntegrateSource_finiteRate: Invalid option for chemistry reduction')
!@       end if
!@#endif
!@       !
!@       ! Compute source terms for post
!@       !
!@       nullify ( prope_tmp )
!@       allocate( prope_tmp(npoin) )
!@       call ker_proper('DENSI','NPOIN',dummi,dummi,prope_tmp)
!@
!@       do ipoin=1,npoin
!@          src_chm(ipoin,:)  = (Y_k_0(ipoin,1:nspec_chm) - conce(ipoin,1:nspec_chm,1)) * dtinv * prope_tmp(ipoin)
!@       end do
!@
!@       deallocate( prope_tmp )
!@       !
!@       ! Update conce with new mass fractions
!@       !
!@       do ipoin=1,npoin
!@          conce(ipoin,1:nspec_chm,1) = Y_k_0(ipoin,1:nspec_chm)
!@       end do
!@
!@       call memory_deallo(mem_modul(1:2,modul),'Y_K_0'  ,'chm_IntegrateSource_finiteRate',Y_k_0)
!@       call memory_deallo(mem_modul(1:2,modul),'T_0'    ,'chm_IntegrateSource_finiteRate',T_0)
!@
!@       if (kfl_tdac_write_chm == 1_ip) then
!@          if (kfl_pfa_chm /= 4_ip) &
!@             call runend('In chm_IntegrateSource_finiteRate: Invalid Reduction option for Table')
!@          if (ittot_chm == mitim ) then
!@               call tab_par_exchange( table_coords, table_fw % main_table )
!@          end if
!@       end if
!@
!@       !
!@       ! Write TDAC Reduced Table
!@       !
!@       if (kfl_tdac_write_chm < 15_ip .and. kfl_tdac_write_chm > 0_ip ) then
!@          if (kfl_pfa_chm /= 4_ip) &
!@             call runend('In chm_IntegrateSource_finiteRate: Invalid Reduction option for Table')
!@          if (ittot_chm == mitim ) then
!@             print *, "Writting Reduced Reaction Table"
!@             nrow = 1_ip
!@             do ii = 1, table_fw % main_table % ndim
!@                nrow = nrow * table_coords(ii) % leng
!@             enddo
!@             shap = 0_ip
!@             do ii = 1, table_fw % main_table % ndim
!@                 shap(ii) = table_coords(ii) % leng
!@             enddo
!@             print'(A,300(I5,1X))', '--| ALYA         TABLE: tab % ndim ',table_fw % main_table % ndim
!@             print'(A,300(I5,1X))', '--| ALYA         TABLE: tab % nvar ', table_fw % main_table % nvar
!@             print'(A,300(I5,1X))', '--| ALYA         TABLE: shap ', shap(1:table_fw % main_table % ndim)
!@             print'(A,I8)', '--| ALYA         TABLE: nrow ', nrow
!@             do ii = 1, nrow
!@                 call i2inds(ii, table_fw % main_table % ndim ,shap,inds)
!@                 do jj = 1, table_fw % main_table % ndim
!@                     test_output(jj) = table_fw % main_table % coords(jj) % x(inds(jj))
!@                 enddo
!@                 write(660+table_fw % main_table % ndim,'(500(2X,E20.8))') test_output(1:table_fw % main_table % ndim),  &
!@                         table_fw % main_table % tab(ii,1:table_fw % main_table % nvar)
!@             enddo
!@          end if
!@       end if
!@
!@    end if
!@
!@  end subroutine chm_IntegrateSource_finiteRate
!@
!@
!@  subroutine chm_heatRelease_finiteRate
!@    use def_domain !,                 only : npoin
!@    use def_master !,                 only : conce, tempe, prthe, speci
!@    use def_chemic,                 only : W_k
!@    use def_chemic,                 only : src_chm, hrr_chm, hrr_int_chm
!@#ifdef CANTERA
!@    use def_chemic,                 only : gas_chm
!@#endif
!@
!@    implicit none
!@    integer(ip)                         :: iclas, ipoin, ielem, igaus, inode
!@    integer(ip)                         :: pnode, pgaus, pelty
!@    real(rp)                            :: gpvol, gpdet, gphre, par_mh(nspec_chm) ! Partial molar enthalpy
!@    real(rp)                            :: elhre(mnode), elcod(ndime,mnode)
!@    real(rp)                            :: gpcar(ndime,mnode,mgaus)
!@    real(rp)                            :: xjaci(ndime,ndime),xjacm(ndime,ndime)
!@
!@
!@    !
!@    ! Compute heat release field
!@    ! heat_release = sum_j=1^N(par_mh_j * omega^molar_j)
!@    ! omega^molar_j = omega^mass_j / W_j = dY_j/dt * rho * volume. Substituting:
!@    ! heat_release = density * volume * sum_j=1^N(par_mh_j * dY_j/dt / W_j)
!@    ! However, the heat release per volume unit is obtained at each cell so volume is omitted.
!@    !
!@
!@    if (INOTMASTER) then
!@
!@       do ipoin=1,npoin
!@          hrr_chm(ipoin) = 0.0_rp
!@       end do
!@#ifdef CANTERA
!@       do ipoin=1,npoin
!@          call cantera_partial_molar_enthalpies(nspec_chm,tempe(ipoin,1),prthe(1),conce(ipoin,:,1),par_mh(:))
!@          do iclas=1,nspec_chm
!@             hrr_chm(ipoin) = hrr_chm(ipoin) - par_mh(iclas) * src_chm(ipoin,iclas) / W_k(iclas)
!@          end do
!@       end do
!@#endif
!@
!@       !
!@       ! Compute heat release integrated over the whole domain
!@       ! Transfer heat release field from nodes to Gaussian points and do the volumetric integral
!@       hrr_int_chm = 0.0_rp
!@       elements: do ielem = 1, nelem
!@          pelty = ltype(ielem)
!@          if (pelty > 0) then
!@             pnode = nnode(pelty)
!@             pgaus = ngaus(pelty)
!@
!@             !
!@             ! Gather operations
!@             !
!@             do inode = 1,pnode
!@                ipoin = lnods(inode,ielem)
!@                elhre(inode) = hrr_chm(ipoin)
!@                elcod(1:ndime,inode) = coord(1:ndime,ipoin)
!@             end do
!@
!@             !
!@             ! 1st and 2nd order Cartesian derivatives, and dV:=GPVOL=|J|*wg
!@             !
!@             do igaus = 1, pgaus
!@                call elmder(&
!@                      pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&        ! Cartesian derivative
!@                      elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)          ! and Jacobian
!@                 gpvol = elmar(pelty)%weigp(igaus) * gpdet               ! |J|*wg
!@                 gphre = 0.0_rp
!@                 do inode = 1, pnode
!@                    gphre = gphre + elmar(pelty) % shape(inode,igaus) * elhre(inode)
!@                 end do
!@                 hrr_int_chm = hrr_int_chm + gphre * gpvol
!@             end do
!@          end if
!@       end do elements
!@
!@    else
!@       hrr_int_chm = 0.0_rp
!@    end if
!@
!@    call pararr('SUM',0_ip,1_ip,hrr_int_chm)
!@
!@  end subroutine chm_heatRelease_finiteRate

end module mod_chm_finiteRate_fast

