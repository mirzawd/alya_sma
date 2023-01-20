!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_elmop2(itask)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_elmop2
  ! NAME
  !   tur_elmop2
  ! DESCRIPTION
  !    Compute elemental matrix and RHS for the Spalart-Allmaras model.
  ! USES
  !    tur_elmgat
  ! USED BY
  !    tur_matrix
  !***
  !-----------------------------------------------------------------------
  use def_elmtyp
  use def_master
  use def_domain
  use def_turbul
  use mod_ker_proper 
  use def_kermod 
  use mod_ADR,    only : ADR_element_assembly
  use mod_ADR,    only : ADR_bubble_assembly
  use mod_ADR,    only : ADR_projections_and_sgs_assembly
  use mod_ADR,    only : ADR_add_sgs_or_bubble
  use mod_ADR,    only : ELEMENT_ASSEMBLY             ! 1
  use mod_ADR,    only : PROJECTIONS_AND_SGS_ASSEMBLY ! 4
  use mod_ADR,    only : BUBBLE_ASSEMBLY              ! 5
  use mod_ADR,    only : mreac_adr
  use mod_matrix, only : matrix_assemble_element_matrix_to_CSR
  use mod_matrix, only : matrix_assemble_element_RHS
  implicit none

  integer(ip), intent(in) :: itask
  real(rp)    :: elmat(mnode,mnode)                       ! Element matrices
  real(rp)    :: elrhs(2*mnode)                           ! For properties
  real(rp)    :: elrpr(2*mnode)

  real(rp)    :: eledd(mnode,2)                           ! Gather 
  real(rp)    :: elcod(ndime,mnode)                       ! Coordinates x
  real(rp)    :: elvel(ndime,mnode)                       ! Velocity u
  real(rp)    :: eltem(mnode)                             ! Temperature T
  real(rp)    :: elwal(mnode)                             ! Distance to wall y
  real(rp)    :: elwai_coo(ndime,mnode)                   ! Coordinates of the nearest wall node
  real(rp)    :: eltur(nturb_tur,mnode,nbdfp_tur+1)       ! Turb. variables
  real(rp)    :: elunk(mnode,2)
  real(rp)    :: elust(mnode)                             ! U*
  real(rp)    :: elrg2(mnode)                             ! (d^2 ui) / ( dxm dxn )
  real(rp)    :: elsqk(mnode)                             ! sqrt(k)
  real(rp)    :: elgrp(ndime,mnode)                       ! grad(phi)
  real(rp)    :: elpro(mnode)                             ! Projection
  real(rp)    :: elprr(mnode)                             ! Projection of reaction (split oss)
  real(rp)    :: elpgr(ndime, mnode)
  real(rp)    :: ellmax(mnode)
  real(rp)    :: elfle(mnode)                             ! Level set
  real(rp)    :: elmsh(ndime,mnode)                       ! u mesh
  real(rp)    :: elprd(mnode)                             ! Smoothed production term

  real(rp)    :: tragl(ndime,ndime),chave(ndime,2)        ! Element length
  real(rp)    :: chale(2),hleng(ndime) 
  real(rp)    :: gpstt(mgaus)                             ! Laplacian

  integer(ip) :: ielem,igaus,ipoin,dummi                  ! Indices and dimensions
  integer(ip) :: pelty,pnode,pgaus,porde,pevat
  integer(ip) :: plapl,pmate,ptopo
  integer(ip) :: elwao(mnode)                             ! Order of the Number of the nearest wall node
  
  real(rp)    :: dtcri                                    ! Critical time step  
  real(rp)    :: fddes(mgaus)                             ! DDES blending factor
  real(rp)    :: gddes(mgaus)                             ! DDES blending function
  real(rp)    :: sstf1(mgaus)                             ! SST blending function
  real(rp)    :: sasso(mgaus)                             ! SAS source term
  real(rp)    :: gpvol(mgaus)                             ! |J|*w,|J|
  real(rp)    :: gprea(mgaus)                             ! r
  real(rp)    :: gpvel(ndime,mgaus), conve(ndime, mgaus)  ! velocity and convetion term
  real(rp)    :: gppro(mgaus)                             ! Weighted residual L2-projection
  real(rp)    :: gpprr(mgaus)                             ! Weighted residual L2-projection
  real(rp)    :: gppgr(ndime, mgaus)                      ! L2-projection turbulence gradient
  real(rp)    :: grtup(ndime, mgaus)
  real(rp)    :: gpdiv(mgaus)                             ! Divergence of convection
  real(rp)    :: gpdif(mgaus)                             ! k
  real(rp)    :: gpgrd(ndime,mgaus)                       ! k, grad(k)
  real(rp)    :: gprhs(mgaus)                             ! f
  real(rp)    :: gprec(mgaus)                             ! reactive term
  real(rp)    :: gpcon(mgaus)                             ! convective term
  real(rp)    :: gpres(mgaus)                             ! residual term
  real(rp)    :: gpgrv(ndime,ndime,mgaus)                 ! grad(a)
  real(rp)    :: gpden(mgaus),gpvis(mgaus)                ! rho,mu
  real(rp)    :: gptur(nturb_tur,nbdfp_tur+1,mgaus)       ! Turb. variables
  real(rp)    :: gpsha(mnode,mgaus)                       ! N
  real(rp)    :: gpder(ndime,mnode,mgaus)                 ! dNk/dsj
  real(rp)    :: gpcar(ndime,mnode,mgaus)                 ! dNk/dxj
  real(rp)    :: gphes(ntens,mnode,mgaus)                 ! dNk/dxidxj
  real(rp)    :: gpcan(mgaus)                             ! canopy flow ( = cd*LAD*|vel|)
  real(rp)    :: sreac(mgaus)                             ! reaction term in perturbation of test function

  real(rp)    :: dgpmut_dtur(nturb_tur,mnode,mgaus)       ! Turbulence viscousity derivatives w.r.t. nodal turbulence unknown
  real(rp)    :: dgpmut_dvel(ndime,mnode,mgaus)           ! Turbulence viscousity derivatives w.r.t. nodal velocity
  real(rp)    :: drea_dtur(nturb_tur,mnode,mgaus)         ! derivatives of reaction term w.r.t turbulent unknown
  real(rp)    :: dsou_dtur(nturb_tur,mnode,mgaus)         ! derivatives of source term w.r.t turbulent unknown
  real(rp)    :: dsou_dvel(ndime,mnode,mgaus)             ! derivatives of source term w.r.t velocity

  real(rp)    :: gprea_diff(ndime,mnode,mgaus)             ! derivatives of reaction term w.r.t nodal coordinate
  real(rp)    :: gprhs_diff(ndime,mnode,mgaus)             ! derivatives of gprhs term w.r.t nodal coordinate
  real(rp)    :: gprhswal_diff(ndime,mnode,mgaus)         ! derivatives of gprhs wall term w.r.t nearest wall nodal coordinate
  real(rp)    :: gpreawal_diff(ndime,mnode,mgaus)         ! derivatives of gprea wall term w.r.t nearest wall nodal coordinate
  
  real(rp)    :: gpvol_der(ndime,mnode,mgaus)          ! Derivative of w*|J|, |J| w.r.t. coordinate nodes
  real(rp)    :: gpcar_der(ndime,mnode,ndime,mnode,mgaus) ! Derivative of gpcar w.r.t. coordinate nodes
  real(rp)    :: hleng_der(ndime,mnode,ndime)          ! Derivative of hleng w.r.t. coordinate nodes
  real(rp)    :: chale_der(ndime,mnode,2)              ! Derivative of chale w.r.t. coordinate nodes
  
  
  real(rp)    :: gptvi(mgaus)                             ! turbulent viscosity
  real(rp)    :: facto
  integer(ip), save :: ipass=0
  
  !
  ! Initialization
  !
  do igaus = 1,mgaus
     gpstt(igaus) = 1.0_rp
     gpdiv(igaus) = 0.0_rp
     gppro(igaus) = 0.0_rp
     gpprr(igaus) = 0.0_rp
     gppgr(1:ndime, igaus) =0.0_rp
     gpcan(igaus) = 0.0_rp
     gpres(igaus) = 0.0_rp
     gprec(igaus) = 0.0_rp
     gpcon(igaus) = 0.0_rp
     grtup(1:ndime, igaus) =0.0_rp
     gptur(:,1:nbdfp_tur+1,igaus) = 0.0_rp
     dgpmut_dtur(1:nturb_tur,1:mnode,igaus) = 0.0_rp
     dgpmut_dvel(1:ndime,1:mnode,igaus) = 0.0_rp
     gpgrd(1:ndime,igaus) = 0.0_rp
     gptvi(igaus) = 0.0_rp
     dsou_dvel(:,:,igaus) = 0.0_rp
     dsou_dtur(:,:,igaus) = 0.0_rp
     drea_dtur(:,:,igaus) = 0.0_rp
     gprea_diff(:,:,igaus) = 0.0_rp
     gprhs_diff(:,:,igaus) = 0.0_rp
     gpreawal_diff(:,:,igaus) = 0.0_rp
     gprhswal_diff(:,:,igaus) = 0.0_rp
  end do
  dtmax_tur = -1.0_rp
  !
  ! Initialize postprocessing values
  !
  if( kfl_ddesm_tur >= 1 .and.  iunkn_tur == 1 ) then
     do ipoin = 1,npoin
        fddes_tur(ipoin) = 0.0_rp
        gddes_tur(ipoin) = 0.0_rp
     end do
  end if
  if( TUR_SST_K_OMEGA ) then
     do ipoin = 1,npoin
        sstf1_tur(ipoin) = 0.0_rp
     end do
     if ( kfl_sasim_tur == 1 ) then
        do ipoin = 1,npoin
           sasso_tur(ipoin) = 0.0_rp
        end do
     end if
  end if
  !
  ! Loop over elements
  !
  elements: do ielem = 1,nelem

     if( lelch(ielem) /= ELHOL ) then
        !
        ! Initialize
        !
        pelty = ltype(ielem)
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)
        porde = lorde(pelty)
        plapl = llapl(pelty)
        ptopo = ltopo(pelty)
        pevat = pnode
        !
        ! Check if element is a solid
        !
        pmate = 1
        if( nmate > 1 ) pmate = lmate(ielem) 
        !
        ! Gather operations
        !
        call tur_elmgat(&
             pnode,lnods(1,ielem),eltur,elvel,elcod,elwal,eledd,&
             elust,eltem,elrg2,elsqk,elgrp,elpro,elunk,elfle,elmsh,&
             elprd, elprr, elpgr, ellmax)
        !
        ! Gather operations on the wall
        !
        if(itask == 30 .and. kfl_walld == 2 .or. kfl_walld == 3) & 
           call tur_elmgat_wall(pnode,lnods(1,ielem),elwao,elwai_coo)
        !
        ! HLENG and TRAGL at center of gravity
        !
        call elmlen(& 
             ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),&
             hleng)
        !
        ! Derivatives of HLENG w.r.t. nodal coordinates
        !
        if(itask == 30) & 
           call elmlen_der(&
                ndime,pnode,elmar(pelty)%dercg,elcod,&
                hnatu(pelty),hleng_der)
        !
        ! Compute the characteristic length CHALE
        !
        call elmchl(&
             tragl,hleng,elcod,elvel,chave,chale,pelty,pnode,&
             porde,hnatu(pelty),kfl_advec_tur,kfl_ellen_tur)
        !
        ! Derivatives of the characteristic length: CHALE w.r.t. nodal coordinates
        ! 
        if(itask == 30) & 
           call elmchl_der(&
                tragl,hleng_der,elcod,elvel,chale_der,pnode,&
                porde,hnatu(pelty),kfl_advec_tur,kfl_ellen_tur)
        !
        ! Local time step
        !
        if( kfl_timco == 2 .and. solve(iunkn_tur) % kfl_algso /= -2) then 
           call tur_elmtss(&
                pelty,pnode,ielem,elvel,eledd,dtcri,eltur,elfle,&
                lnods(:,ielem),chale,hleng)
           ! Maximum time step between that given by the global safety factor saflo_nsi, 
           ! and local safet_nsi
           facto = safet_tur/safma_tur
!           facto = 1.0_rp
           dtinv_tur = min(1.0_rp / (dtcri*safet_tur), 1.0_rp/(dtcri_tur*saflo_tur*sqrt(facto)))
           if( kfl_stead_tur == 1 ) dtinv_tur = 0.0_rp
           if( kfl_timei_tur == 0 ) dtinv_tur = 0.0_rp
           ! Stores maximum time step
           if( kfl_adj_prob == 0  ) dtmax_tur = max(dtmax_tur, 1.0_rp/dtinv_tur)           
        end if
        !
        ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
        !
        gphes = 0.0_rp
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
             
        !call elmcar(&
        !     pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
        !     elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
        !     gphes,ielem)
        ! 
        !  Properties: GPDEN and GPVIS
        !
        call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,gpsha,gpcar)
        call ker_proper('VISCO','PGAUS',dummi,ielem,gpvis,pnode,pgaus,gpsha,gpcar)
        ! porosity for canopy flow
        call ker_proper('POROS','PGAUS',dummi,ielem,gpcan,pnode,pgaus,gpsha,gpcar)        
        call ker_proper('TURBU','PGAUS',dummi,ielem,gptvi,pnode,pgaus,gpsha,gpcar)  
        do igaus =1, pgaus
           gptvi(igaus) = gpden(igaus) *gptvi(igaus)
        end do
!         call ker_proper('GRTUR','PGAUS',dummi,ielem,gpgrd,pnode,pgaus,gpsha,gpcar)         
        if (kfl_adj_prob == 1) then
          call ker_proper('TDRTU','PGAUS',dummi,ielem,dgpmut_dtur,pnode,pgaus,gpsha,gpcar)
        endif
            
        if( itask >= 10 .and. itask <= 20 ) then
           !
           ! Properties
           !
           call asspro(&
                itask,pnode,2_ip,pgaus,lnods(:,ielem),lelch(ielem),gpden,&
                gpvis,gpvol,gpsha,elrpr,rhsid)

        else
           !
           ! Equation coefficients: GPDIF, GPREA, GPRHS, GPGRD and GPVEL(NDIME)
           !
           do igaus = 1,pgaus     
              call tur_elmco2(&
                   pnode,gpvis(igaus),gpden(igaus),elmar(pelty)%shape(1,igaus),  &
                   gpcar(1,1,igaus),elvel,elmsh,eltur,eledd,elwal,elust,eltem,   &
                   elrg2,elsqk,elgrp,elpro,elprd,gpvol(igaus),chale,gpdif(igaus),&
                   gpgrd(1,igaus),gprea(igaus),gpvel(1,igaus),gprhs(igaus),      &
                   gpgrv(1,1,igaus),gptur(1,1,igaus),gppro(igaus),hleng,fddes(igaus),&
                   gddes(igaus),sstf1(igaus),sasso(igaus), gpcan(igaus), sreac(igaus), &
                   gpprr(igaus), elprr,gppgr(1, igaus),elpgr, pmate, turvi_tur(1, igaus, ielem), &
                   ipass, conve(1, igaus), gptvi(igaus), ellmax)
           end do
           if ( kfl_ddesm_tur >= 1 .and.  iunkn_tur == 1) then    ! DDES model flag
              call tur_pprdes(lnods(:,ielem),pnode,fddes,gddes,pgaus,gpvol,gpsha)
           end if
           if ( TUR_SST_K_OMEGA ) then
              call tur_pprsst(1_ip,lnods(:,ielem),pnode,sstf1,pgaus,gpvol,gpsha)
           end if
           if ( kfl_sasim_tur == 1 ) then
              call tur_pprsst(2_ip,lnods(:,ielem),pnode,sasso,pgaus,gpvol,gpsha)
           end if
           !
           ! Assembly 
           !
           if( TUR_K_EPS_STD.or.TUR_TKE_SGS ) then
              !
              ! kfl_ortho_tur= -3 for supg
              ! kfl_ortho_tur=  0 for ASGS
              ! kfl_shock_tur=  0 , 1, 2  no, isotropic, and  anisotropic shock captu
              ! shock_tur    = shock capturing parameter 
              !
              ! SUPG 
              if( kfl_ortho_tur == -3 )   sreac(1:pgaus) = 0.0_rp 
           else 
              !
              ! Other models
              !
              if( kfl_ortho_tur == -3 ) then !SUPG
                 sreac(1:pgaus) = 0.0_rp                
              else                       
                 sreac(1:pgaus) = gprea(1:pgaus)
              end if
           end if

           call tur_elmmsu(&
                kfl_ortho_tur,kfl_shock_tur,shock_tur,pnode,plapl,pgaus,gpvel,gpdif,gprea,&
                gptur,gpgrd,gprhs,gpden,gpsha,gpcar,gpvol,elmat,elrhs,chale,&
                elunk,gphes,sreac,itask,gpres,gprec,gpcon,gppro,gpprr, gppgr,grtup, conve)
           
    
           if (kfl_cos_opt == 1) call tur_elmcost_all(eltur,pnode,pgaus,gpsha,gpvol,lnods(:,ielem),costf)
           
           !
           !  calculate terms related to the adjoint generally: dF/dU, dR/dD and additional terms related to dR/dU
           !              
           if(kfl_adj_prob == 1 ) then
             !
             !  Calculate derivatives of reaction term, source term, miu_tur at each gauss point respect to K, E, U and coordinates at each node
             !
             do igaus = 1,pgaus
               call tur_elmco2_der_all( &
                  itask,pnode,gpvis(igaus),gpden(igaus),elmar(pelty)%shape(1,igaus),gpcar(1,1,igaus),elvel,elmsh,eltur,eledd,&
                  elwal,elust,eltem,elrg2,elsqk,elgrp,elpro,elprd,gpvol(igaus),&
                  chale,gpdif(igaus),gpvel(1,igaus),gpgrv(1,1,igaus),gptur(1,1,igaus),&
                  hleng,gpcan(igaus), elprr, pmate, gpcar_der(1,1,1,1,igaus),&
                  drea_dtur(1,1,igaus),dsou_dtur(1,1,igaus), &
                  dsou_dvel(1,1,igaus),gptvi(igaus),dgpmut_dtur(1,1,igaus),&
                  dgpmut_dvel(1,1,igaus),gprea_diff(1,1,igaus),gprhs_diff(1,1,igaus),gpreawal_diff(1,1,igaus),gprhswal_diff(1,1,igaus),&
                  elcod,elwai_coo)
             enddo
             !
             !  Calculate matrix dRk/dK and dRw/dW as well as vectors [dRk/dV][Lambda_k] and [dRw/dV][Lambda_w]
             !
             if(itask /= 30) &
               call tur_elmmsu_der_all(&
                  kfl_ortho_tur,kfl_shock_tur,shock_tur,pnode,plapl,pgaus,gpvel,gpdif,gprea,&
                  drea_dtur,dsou_dtur,dsou_dvel,gptvi,dgpmut_dtur,dgpmut_dvel,&
                  gptur,gpgrd,gprhs,gpden,gpsha,gpcar,gpvol,elmat,elrhs,chale,&
                  gphes,sreac,eltur,ielem)
             
             !
             !  calculate cost function derivatives w.r.t. unknowns dF/dU
             !
             call tur_elmdcost_all(eltur,pnode,pgaus,gpsha,gpvol,lnods(:,ielem),elrhs)
             !
             !  calculate residal derivatives (partial) w.r.t. unknowns dR/dD
             !
             call tur_resdiff_all( &
                  kfl_ortho_tur, kfl_shock_tur,shock_tur, pnode,plapl,pgaus,gpvel,gpdif,&
                  gprea,gptur,gpgrd,gprhs,gpden,gpsha,gpcar,gpvol,chale, &
                  elunk, gphes, sreac, itask, gppro, gpprr, gppgr,gptvi,lnods(:,ielem),eltur,elwal,gpvis, &
                  gprea_diff,gprhs_diff, gpreawal_diff,gprhswal_diff,gpvol_der, gpcar_der,chale_der(:,:,2),elwao,sens_mesh,sens_wall,ielem)
                    
           endif !kfl_adj_prob
                      
           if( itask == 4 ) then 
              !
              ! Projection
              ! 
              call tur_elmort(gpres,gpcon,gprec,grtup,gprhs,ielem,pgaus,pnode,gpsha,gpvol)

           else if( itask == 1 ) then 
              !
              ! Assembly
              !
              if( solve(iunkn_tur) % kfl_iffix == 0 ) &
                   call tur_elmdir(&
                   1_ip,pnode,pevat,lnods(:,ielem),elmat,elrhs)
              !
              ! elmat = Transpose[elmat] for adjoint
              !
              if (kfl_adj_prob == 1) elmat = transpose(elmat)
              
              call matrix_assemble_element_RHS(&
                   solve(1) % ndofn,solve(1) % ndofn,pnode,lnods(:,ielem),elrhs,rhsid)
              call matrix_assemble_element_matrix_to_CSR(&
                   kfl_element_to_csr,1_ip,pnode,pnode,&
                   ielem,lnods(:,ielem),elmat,r_dom,c_dom,amatr,lezdo)
           end if
           
        end if

     end if

  end do elements
  ipass=1

end subroutine tur_elmop2
!-----------------------------------------------------------------------
! NOTES
! 
! Governing equation:
! L(T) = rho*u/dt + rho*cp*a.grad(u) - div[k*grad(T)] + r*u = f
! The equation is stabilized using the ASGS model:
!  +-            +-                        +-             
!  | L(u)*v dw - | tau*L'(v)*(L(u)-f) dw = | f*v dw 
! -+            -+                        -+
! 
! where L'(v) = -a.grad(v) - k*Lapl(v) - grad(k).grad(v) + r*v
! Let gpadj = -tau*L'(v)
! The weak form is:
!                                        Gal. ASGS
!                                         |    |
!  +-                                     |    |
!  | ( rho*u/dt + rho*a.grad(u) + r*u )*( v + gpadj ) dw
! -+ 
!     +-                            +-
!   + | -div[k*grad(u)]*gpadj dw +  | k*grad(u).grad(v) dw 
!    -+                            -+
!     +-                              
!   = | ( f + rho*u^(n-1)/dt )*( v + gpadj ) dw 
!    -+                   
!
!***
!-----------------------------------------------------------------------




