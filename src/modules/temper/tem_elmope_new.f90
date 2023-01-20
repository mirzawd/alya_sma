!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_elmope_new(order)
  !------------------------------------------------------------------------
  !****f* Temper/tem_elmop2
  ! NAME 
  !    tem_elmop2
  ! DESCRIPTION
  !    ORDER=1:
  !      Temperature equation, elemental operations:
  !      1. Compute elemental matrix and RHS 
  !      2. Impose Dirichlet boundary conditions
  !      3. Assemble them
  !    ORDER=4:
  !      Update the subgrid scale
  ! USES
  ! USED BY
  !    tem_matrix
  !------------------------------------------------------------------------
  use def_parame
  use def_elmtyp
  use def_master, only : kfl_timco, cutim, solve, solve_sol, rhsid, amatr, nspec
  use def_domain
  use def_kermod
  use mod_ker_proper 
  use def_temper
  use mod_tem_entropy, only : tem_entropy_viscosity
  use mod_ADR,    only : ADR_element_assembly
  use mod_ADR,    only : ADR_bubble_assembly
  use mod_ADR,    only : ADR_projections_and_sgs_assembly
  use mod_ADR,    only : ADR_add_sgs_or_bubble
  use mod_ADR,    only : ELEMENT_ASSEMBLY             ! 1
  use mod_ADR,    only : PROJECTIONS_AND_SGS_ASSEMBLY ! 4
  use mod_ADR,    only : BUBBLE_ASSEMBLY              ! 5
  use mod_ADR,    only : mreac_adr
  use mod_solver, only : solver_assemble_element_matrix_scalar
  use mod_matrix, only : matrix_assemble_element_RHS
  use mod_matrix, only : matrix_assexp
  use mod_tem_turbul
  implicit none

  integer(ip), intent(in) :: order                     ! =2: compute SGS only

  real(rp)    :: elmat(mnode,mnode),elrhs(mnode)
  integer(ip) :: ielem,igaus                     ! Indices and dimensions
  integer(ip) :: pelty,pmate,pnode
  integer(ip) :: pgaus,plapl,porde,ptopo
  integer(ip) :: kfl_advec_old

  real(rp)    :: eltem(mnode,ADR_tem % ntime)          ! Gather 
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: elvel(ndime,mnode)
  real(rp)    :: elmsh(ndime,mnode)

  real(rp)    :: tragl(ndime,ndime),chave(ndime,2)     ! Stabilization
  real(rp)    :: chale(2),hleng(ndime)

  real(rp)    :: gpdtc                                 ! Values at Gauss points
  real(rp)    :: gpvol(mgaus)                          ! |J|*w
  real(rp)    :: gprea(mgaus,mreac_adr)                ! r
  real(rp)    :: gptke(mgaus)                          ! tke-sgs 
  real(rp)    :: gpvel(ndime,mgaus)                    ! a
  real(rp)    :: gpcon(mgaus),gpcod(ndime,mgaus)       ! k
  real(rp)    :: gpdif(mgaus),gpgrd(ndime,mgaus)       ! k+kt, grad(k+kt)
  real(rp)    :: gpdiv(mgaus)                          ! Divergence of convection
  real(rp)    :: gprhs(mgaus)                          ! f (all terms)
  real(rp)    :: gpden(mgaus)                          ! rho and then rho*cp
  real(rp)    :: gpsph(mgaus)                          ! cp
  real(rp)    :: gptem(mgaus,ADR_tem % ntime)          ! T
  real(rp)    :: gpgrt(ndime,mgaus)                    ! grad(T)
  real(rp)    :: gpsgv(ndime,mgaus)                    ! u'
  real(rp)    :: gpsha(mnode,mgaus)                    ! N
  real(rp)    :: gpder(ndime,mnode,mgaus)              ! dNk/dsj
  real(rp)    :: gpcar(ndime,mnode,mgaus)              ! dNk/dxj
  real(rp)    :: gphes(ntens,mnode,mgaus)              ! dNk/dxidxj
  real(rp)    :: gpsta(mgaus)                          ! tau
  real(rp)    :: gptur(mgaus)                          ! Turbulent viscosity
  real(rp)    :: gpmsh(ndime,mgaus)                    ! Mesh velocity
  integer(ip) :: dummi
  real(rp)    :: dtmin

  
  real(rp)    :: dsou_dtem(mnode,mgaus)                ! d(hW)/dT
  real(rp)    :: dsou_dcon(nspec,mnode,mgaus)          ! d(hW)/dYk
  real(rp)    :: gpdde(mnode,mgaus)                    ! Density derivatives w.r.t nodal temperature
  real(rp)    :: eltem_forw(mnode)                     ! tempe_forw at elemental nodes
  real(rp)    :: densi(mgaus)
#ifdef EVENT
  call mpitrace_user_function(1)
#endif

  !
  ! Initialization
  !
  gpdif = 0.0_rp
  gpsph = 0.0_rp
  gpden = 0.0_rp
  gprea = 0.0_rp
  gpsta = 0.0_rp
  gpdiv = 0.0_rp
  gprhs = 0.0_rp
  gptur = 0.0_rp
  gptem = 0.0_rp
  gpgrd = 0.0_rp
  gpvel = 0.0_rp
  gpmsh = 0.0_rp

  gpdde = 0.0_rp
  dsou_dtem = 0.0_rp
  dsou_dcon = 0.0_rp
  eltem_forw = 0.0_rp
  !
  ! Loop over elements
  !  
  dtmin = 1.0e6_rp
  elements: do ielem = 1,nelem
     !
     ! Element dimensions
     !
     pelty = ltype(ielem)

     if( pelty > 0 ) then
        pnode = lnnod(ielem)
        pgaus = ngaus(pelty)
        plapl = llapl(pelty)
        porde = lorde(pelty)
        ptopo = ltopo(pelty)
        pmate = lmate(ielem)
        kfl_advec_old = kfl_advec_tem
        !
        ! Gather operations 
        !
        call tem_elmgat(&
             ielem,pnode,lnods(1,ielem),eltem,elvel,elcod,elmsh)
        !
        ! hleng and tragl at center of gravity
        !
        call elmlen(ndime,pnode,elmar(pelty) % dercg,tragl,elcod,&
             hnatu(pelty),hleng)
        !
        ! Compute the characteristic length CHALE
        !
        if( kfl_ellen_tem == -1 ) then 
           call tem_elmchl(&
                ielem,pelty,pnode,plapl,pmate,lnods(1,ielem),elcod,eltem,&
                elvel,gpcar,gphes,chale)
        else
           call elmchl(tragl,hleng,elcod,elvel,chave,chale,pelty,pnode,&
                porde,hnatu(pelty),kfl_advec_tem,kfl_ellen_tem)
        end if
        !
        ! Local time step DTINV_TEM
        !
        if(kfl_timco==2) then 
           call tem_elmtss(&
                ielem,pelty,pnode,pmate,elcod,elvel,&
                eltem,gpcar,chale,gpdtc)        
           dtinv_tem = 1.0_rp/(gpdtc*safet_tem)
           if( kfl_stead_tem == 1 ) dtinv_tem = 0.0_rp
           if( kfl_timei_tem == 0 ) dtinv_tem = 0.0_rp
        end if
        !
        ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
        !
        if( plapl == 0 ) gphes = 0.0_rp
        call elmca2(&
             pnode,pgaus,plapl,elmar(pelty) % weigp,elmar(pelty) % shape,&
             elmar(pelty) % deriv,elmar(pelty) % heslo,elcod,gpvol,gpsha,&
             gpder,gpcar,gphes,ielem)
        !
        ! Temperature: GPTEM
        !
        call gather(&
             2_ip,pgaus,pnode,1_ip,dummi,gpsha,eltem,gptem)   
        !
        ! Add SGS or bubble to current unknown
        !
        call ADR_add_sgs_or_bubble(&
             ielem,pgaus,elmar(pelty) % shape_bub,ADR_tem,gptem)
        !
        ! Properties: GPDEN, GPDIF, GPGRD and GPREA
        !
        call ker_proper('DENSI','PGAUS',dummi,ielem,densi,pnode,pgaus,gpsha,gpcar)
        call ker_proper('CONDU','PGAUS',dummi,ielem,gpcon,pnode,pgaus,gpsha,gpcar)
        call ker_proper('SPHEA','PGAUS',dummi,ielem,gpsph,pnode,pgaus,gpsha,gpcar)
        call ker_proper('TURBU','PGAUS',dummi,ielem,gptur,pnode,pgaus,gpsha,gpcar)

        if ( kfl_entpred_tem == 1_ip ) then
           gpcon = 0.0_rp
        end if

        if(kfl_grdif_tem /= 0)  call ker_proper('GRCON','PGAUS',dummi,ielem,gpgrd,pnode,&
             pgaus,gpsha,gpcar)        
        !
        ! Reaction term
        !
        call tem_elmrea( &
             1_ip,pnode,pgaus,1_ip,pgaus,elvel,densi,gpcar,&
             gprea)
        !
        ! Equation coefficients in GP
        !
        call tem_elmpre(&   ! Care! Here densi is density, and gpden can be rhocp, rhocv, ... 
             ielem,pnode,pgaus,pmate,densi,gpden,gpsph,gpsgv,gpsha,gpcar,gphes,elvel,eltem,&
             elcod,elmsh,gpvel,gptem,gprhs,gpcod,gpgrt,lnods(1,ielem),gpmsh,gprea, gptke)        
        !
        ! Coupling with turbul
        !
        call tem_turbul(&
             pnode,pgaus,1_ip,pgaus,gpcon,gpsph,gpdif,gpgrd,densi,gptur, gpgrt, gptke, hleng)
        
        if ( kfl_entpred_tem == 1_ip ) gpdif = 0.0_rp
        !
        ! Coupling with RADIAT
        !   
        call tem_radiat(&
             ielem,pnode,pgaus,lnods(1,ielem),gpsha,gprhs)
        !
        ! Coupling with CHEMIC
        !
        call tem_chemic(ielem,pgaus,gprhs)      
        !
        ! Entropy stable viscosity
        !
        if ((kfl_entropy_tem == 1_ip)) then   !!!!.and.(kfl_diven_tem == 0)
           call tem_entropy_viscosity(ielem,pnode,pgaus,1_ip,pgaus,gpsha,gpcar,elvel,gpden,hleng,gpvol,gpdif)
           if ( kfl_entpred_tem == 1_ip ) gpdif = 0.0_rp
        end if

        !
        ! Scale source terms with rho (enthalpy) or rho*cp (temperature)
        !
        if ( kfl_rhs_scal_tem > 0 ) then
           do igaus=1,pgaus
              gprea(igaus,:) = gprea(igaus,:) / gpden(igaus)
              gprhs(igaus)   = gprhs(igaus)   / gpden(igaus)
           end do
        end if

        if( order == ELEMENT_ASSEMBLY ) then
           ! 
           ! Assemble equation
           !
           call ADR_element_assembly(&
                ielem,pnode,pgaus,elcod,gpsha,gpcar,elmar(pelty) % deriv,gphes,gpvol,chale,&
                elmar(pelty) % shape_bub,elmar(pelty) % deriv_bub, ADR_tem,&
                cutim,gpden,gpvel,gpdif,gpgrd,gprea,gprhs,gptem,eltem,elmat,elrhs, eladv=elvel)
          
           !
           !  calculate cost function F
           !
           if (kfl_cos_opt == 1) call tem_elmcost_all(eltem,pnode,pgaus,gpsha,gpvol,lnods(1,ielem))
           !
           !  calculate terms related to the adjoint generally: dF/dU, dR/dD and additional terms related to dR/dU
           !	              
           if(kfl_adj_prob == 1 ) then

              !
              !  Calculate derivatives of rho,Cp,K,W at each gauss point respect to T,U,Y at each node
              !
              call tem_der_aux_all(ielem,pgaus,pnode,gpsha,gpcar,dummi,dsou_dtem,dsou_dcon,gpdde)
              !
              !  Calculate matrix dRt/dT and vectors [dRt/dYs] [Lambda_t] and [dRt/dV] [Lambda_t]
              !
              call tem_elmmul_der_all(&
                 pnode     ,pgaus    ,gpvel ,gpdif    ,gprea, &
                 gptem, gpgrd    ,gpden ,gpsha    ,gpcar,gpvol,&
                 elmat, elrhs, chale    ,gphes ,eltem,&
                 gpsph, ielem,lnods(1,ielem),gpdde,dsou_dtem,dsou_dcon,eltem_forw)
              !
              !  calculate cost function derivatives w.r.t. unknowns dF/dU
              !
              call tem_elmdcost_all(pnode,pgaus,gpsha,gpvol,elrhs,eltem_forw)
              !
              !  calculate residal derivatives (partial) w.r.t. unknowns dR/dD
              !
              call tem_resdiff_all(pnode,pgaus,gpden,gpsha,gpdif,gpcar,lnods(1,ielem),ielem,gprea, &
                 chale,gpvel,gpcod,gpvol,gptem,gprhs,gphes,eltem_forw)
           endif


           !
           ! Projections of rho*cp/dt
           !
           if(kfl_explicit_tem == 1) then
              call tem_rhocpdt(                                            &
                 ielem,pnode,pgaus,porde,lnods(:,ielem),gpsha,gpvol,gpden, &
                 dt_rho_cp_tem)
           end if
           !
           ! Prescribe Dirichlet boundary conditions
           !
           if( solve(1) % kfl_iffix == 0 ) &
              call tem_elmdir(&
              pnode,lnods(1,ielem),elmat,elrhs,ielem)
           
           if( solve_sol(1)  %  kfl_algso == 9 ) then
              call tem_assdia(&
                 pnode,elmat,eltem,elrhs)
              call matrix_assemble_element_RHS(&
                 1_ip,1_ip,pnode,lnods(:,ielem),elrhs,rhsid)  
           else
              !
              ! elmat = Transpose[elmat] for adjoint
              !
              if (kfl_adj_prob == 1) elmat = transpose(elmat)
              !
              ! Assembly
              !
              if(kfl_explicit_tem ==1 ) then
                 call matrix_assexp(1_ip,1_ip,pnode,npoin,lnods(1:pnode,ielem),elrhs,elmat,eltem(1:pnode,1),rhsid) 

              else
                 call matrix_assemble_element_RHS(&
                      1_ip,1_ip,pnode,lnods(:,ielem),elrhs,rhsid)
                 call solver_assemble_element_matrix_scalar(&
                      solve(1),1_ip,pnode,pnode,ielem,lnods(:,ielem),&
                      elmat,amatr)
              end if
           end if
           
        else if ( order == PROJECTIONS_AND_SGS_ASSEMBLY ) then
           !
           ! Assemble residual projection
           !
           call ADR_projections_and_sgs_assembly(&
              ielem,pnode,pgaus,elcod,gpsha,gpcar,gphes,gpvol,chale,ADR_tem,&
              cutim,gpden,gpvel,gpdif,gpgrd,gprea,&
              gprhs,gptem,eltem)

        else if ( order == BUBBLE_ASSEMBLY ) then
           !
           ! Update bubble
           !
           call ADR_bubble_assembly(&
              ielem,pnode,pgaus,elcod,gpsha,gpcar,elmar(pelty) % deriv,gphes,gpvol,chale,&
              elmar(pelty) % shape_bub,elmar(pelty) % deriv_bub,ADR_tem,&
              cutim,gpden,gpvel,gpdif,gpgrd,gprea,gprhs,&
              gptem,eltem,elmat,elrhs)

        end if

        kfl_advec_tem=kfl_advec_old
     end if

  end do elements

#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine tem_elmope_new

