!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_elmcor.f90
!> @author  Guillaume Houzeaux
!> @brief   Velocity correction
!> @details Correct velocity in order to obtain a continuity 
!>          conservation scheme
!>          KFL_PREDI_NSI = 2 ... u = u* - (dt/rho)*[ grad(p^i+1) - grad(p^i) ]
!>                        = 3 ... u = u* - tau'1   *[ grad(p^i+1) - grad(p^i) ]
!>                                with tau' = 1/( rho/dt + 1/tau )
!> @} 
!-----------------------------------------------------------------------

subroutine nsi_elmcor() 
  use def_kintyp,               only : ip,rp
  use def_master,               only : unkno,INOTEMPTY,amatr
  use def_master,               only : kfl_timco
  use def_master,               only : times
  use def_elmtyp,               only : ELHOL
  use def_kermod,               only : kfl_kxmod_ker
  use def_domain,               only : mnode,mgaus,ndime
  use def_domain,               only : lnods,ltype,elmar,lelch
  use def_domain,               only : hnatu,ltypb,lnodb
  use def_domain,               only : ngaus,nnode,lorde,nboun
  use def_domain,               only : npoin,lmate,ntens
  use mod_parall,               only : par_omp_ia_colors,par_omp_ja_colors
  use mod_parall,               only : par_omp_num_colors
  use def_nastin,               only : kfl_fixno_nsi,kfl_fixbo_nsi,bvess_nsi
  use def_nastin,               only : saflo_nsi,safet_nsi,kfl_timei_nsi
  use def_nastin,               only : kfl_stead_nsi,kfl_predi_nsi
  use def_nastin,               only : kfl_grvir_nsi,kfl_ellsh_nsi
  use def_nastin,               only : kfl_corre_nsi
  use def_nastin,               only : kfl_advec_nsi,dtinv_nsi
  use def_nastin,               only : dtcri_nsi,ndbgs_nsi
  use def_nastin,               only : kfl_fixrs_nsi
  use mod_nsi_schur_complement, only : xx
  use mod_ker_proper,           only : ker_proper
  use mod_local_basis,          only : local_basis_global_to_local
  use mod_local_basis,          only : local_basis_local_to_global

  implicit none
  integer(ip)       :: ielem,dummi,inodb,iboun           ! Indices and dimensions
  integer(ip)       :: pnode,pgaus,pelty,pmate,porde
  integer(ip)       :: idime,idofn,ipoin,icolo,kelem
  real(rp)          :: dummr, dtcri
  real(rp)          :: elcod(ndime,mnode)                ! x
  real(rp)          :: elvel(ndime,mnode)                ! u
  real(rp)          :: elpre(mnode)                      ! p
  real(rp)          :: eltem(mnode)                      ! T
  real(rp)          :: gpsha(mnode,mgaus)                ! N
  real(rp)          :: gpder(ndime,mnode,mgaus)          ! dN/dsi         
  real(rp)          :: gpcar(ndime,mnode,mgaus)          ! dN/dxi
  real(rp)          :: gphes(ntens,mnode,mgaus)          ! d^2N/dxidxj
  real(rp)          :: gpvol(mgaus)                      ! w*|J|, |J|
  real(rp)          :: gpgpr(ndime,mgaus)                ! grad(p)
  real(rp)          :: gppor(mgaus)                      ! sig
  real(rp)          :: gpgvi(ndime,mgaus)                ! grad(mu+mut)
  real(rp)          :: gpgve(9,mgaus)                    ! grad(u)
  real(rp)          :: gpmut(mgaus)                      ! mut
  real(rp)          :: grvis(ndime,mgaus)                ! grad(mut)
  real(rp)          :: gpden(mgaus)                      ! rho
  real(rp)          :: gpvis(mgaus)                      ! mu
  real(rp)          :: chale(3)                          ! Characteristic lengths
  real(rp)          :: hleng(3)                          ! Element lengths (max to min)
  real(rp)          :: chave(3)                          ! Characteristic velocity
  real(rp)          :: dtinv_loc                         ! Time step
  real(rp)          :: tragl(9)
  real(rp), pointer :: my_gpsha(:,:),my_gpder(:,:,:),gpwei(:)

  if( kfl_corre_nsi == 0 ) return 
  call times(4) % ini() 

  if ( INOTEMPTY ) then

     !-------------------------------------------------------------------
     !
     ! Assemble RHS: \int_V tau1'*[ grad(p^i+1) - grad(p^i) ].v dx
     ! 
     !-------------------------------------------------------------------

     do idofn = 1,ndbgs_nsi
        xx(idofn) = 0.0_rp                               ! Initialization
     end do
     dtinv_loc = dtinv_nsi

     colors: do icolo = 1,par_omp_num_colors   
        ! 
        ! Loop over elements
        !     
        !$OMP  PARALLEL DO SCHEDULE (STATIC)                                        & 
        !$OMP  DEFAULT      ( NONE )                                                &
        !$OMP  FIRSTPRIVATE ( dtinv_loc )                                           &
        !$OMP  PRIVATE      ( kelem, pelty, pnode, porde, gpsha, my_gpsha,          &
        !$OMP                 gpder, gpwei, pgaus, pmate, elcod, elpre, elvel,      &
        !$OMP                 eltem, gpvol, gpcar, ielem, dummr,             &
        !$OMP                 tragl, hleng, chave, chale, dtcri, gpden, gpvis,      &
        !$OMP                 gppor, gpgvi, grvis, gpgve, gpmut, dummi,             &
        !$OMP                 gpgpr, gphes, my_gpder )                              &
        !$OMP  SHARED       ( icolo, par_omp_ia_colors, par_omp_ja_colors, ltype,   &
        !$OMP                 kfl_corre_nsi, elmar, ngaus, lmate,                   &
        !$OMP                 lelch, lnods, kfl_predi_nsi, kfl_timco,               &
        !$OMP                 hnatu, kfl_advec_nsi, kfl_ellsh_nsi,  &
        !$OMP                 kfl_kxmod_ker, kfl_grvir_nsi, safet_nsi,              &
        !$OMP                 dtcri_nsi, kfl_timei_nsi, saflo_nsi,                  &
#ifndef NDIMEPAR
        !$OMP                 ndime,                                                &
#endif
        !$OMP                 kfl_stead_nsi, nnode, lorde, xx )                                               

        elements: do kelem = par_omp_ia_colors(icolo),par_omp_ia_colors(icolo+1)-1
           ielem = par_omp_ja_colors(kelem) 
           pelty = ltype(ielem)                             ! Element properties and dimensions
           if( pelty > 0 ) then
              pnode = nnode(pelty)
              porde = lorde(pelty)
              pmate = lmate(ielem)

              if( kfl_corre_nsi == 1 ) then                 ! Close rule for RHS
                 my_gpsha => elmar(pelty) % shapc
                 my_gpder => elmar(pelty) % deric
                 gpwei    => elmar(pelty) % weigc 
                 pgaus    =  nnode(pelty) 
              else                                          ! Open rule for RHS
                 my_gpsha => elmar(pelty) % shape
                 my_gpder => elmar(pelty) % deriv
                 gpwei    => elmar(pelty) % weigp
                 pgaus    =  ngaus(pelty)
              end if

              if( lelch(ielem) /= ELHOL ) then
                 !
                 ! Gather
                 !
                 call nsi_elmgap(&              
                      pnode,pmate,lnods(1,ielem),elcod,elpre,&
                      elvel,eltem)
                 !
                 ! GPVOL, GPCAR: Jabobian, Cartesian derivatives
                 ! 
                 call elmca2(&
                      pnode,pgaus,0_ip,gpwei,my_gpsha,&
                      my_gpder,elmar(pelty)%heslo,elcod,gpvol,gpsha,&
                      gpder,gpcar,gphes,ielem)
                 !
                 ! CHALE: characteristic length
                 !
                 if( kfl_predi_nsi == 3 .or. kfl_timco == 2 ) then   
                    call elmlen(&
                         ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                         hnatu(pelty),hleng)
                    call elmchl(&
                         tragl,hleng,elcod,elvel,chave,chale,pelty,pnode,&
                         porde,hnatu(pelty),kfl_advec_nsi,kfl_ellsh_nsi)
                 end if
                 !
                 ! Local time step: DTINV_LOC
                 ! 
                 if( kfl_timco == 2 ) then 
                    call nsi_elmtss(&
                         pelty,pmate,pnode,lnods(1,ielem),ielem,elcod,elvel,&
                         gpcar,chale,hleng,dtcri)
                    dtinv_loc = min(1.0_rp / (dtcri*safet_nsi), 1.0_rp/(dtcri_nsi*saflo_nsi))
                    if( kfl_stead_nsi == 1 ) dtinv_loc = 0.0_rp
                    if( kfl_timei_nsi == 0 ) dtinv_loc = 0.0_rp
                 end if

                 call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,gpsha,gpcar)
                 call ker_proper('VISCO','PGAUS',dummi,ielem,gpvis,pnode,pgaus,gpsha,gpcar)
                 call ker_proper('POROS','PGAUS',dummi,ielem,gppor,pnode,pgaus,gpsha,gpcar)
                 call ker_proper('GRVIS','PGAUS',dummi,ielem,gpgvi,pnode,pgaus,gpsha,gpcar)
                 call ker_proper('TURBU','PGAUS',dummi,ielem,gpmut,pnode,pgaus,gpsha,gpcar)  
                 call ker_proper('GRTUR','PGAUS',dummi,ielem,grvis,pnode,pgaus,gpsha,gpcar) 

                 call nsi_turbul(&
                      1_ip,0_ip,pnode,pgaus,1_ip,pgaus,&
                      gpsha,gpcar,elvel,gpden,gpvis,gpmut,&
                      gpgvi,grvis,gpgve,ielem,kfl_kxmod_ker)
                 
                 if( kfl_grvir_nsi == 0 ) gpgvi = 0.0_rp ! zero viscosity gradient

                 call nsi_elmrhc(&                    
                      ielem,pgaus,pnode,lnods(1,ielem),gpvol,gpsha,&
                      gpcar,elvel,chale,gpden,gpvis,gppor,gpgpr,&
                      dtinv_loc,xx)

              end if

           end if

        end do elements
        !$OMP  END PARALLEL DO 

     end do colors

  end if

  !-------------------------------------------------------------------
  !
  ! Do not modify weak form no-penetration
  ! 
  !-------------------------------------------------------------------

  do iboun = 1,nboun
     if( kfl_fixbo_nsi(iboun) == 18 ) then
        do inodb = 1,nnode(abs(ltypb(iboun)))
           ipoin = lnodb(inodb,iboun)
           idofn = (ipoin-1)*ndime
           do idime = 1,ndime
              idofn = idofn + 1
              !xx(idofn) = 0.0_rp
           end do
        end do
     end if
  end do

  !-------------------------------------------------------------------
  !
  ! Impose boundary conditions and solve system using diagonal solver
  ! 
  !-------------------------------------------------------------------

  call nsi_inisol(5_ip)
  !
  ! Beware as the solver is richardson : unkno(k+1) = unkno(k) + inv(Ml) * rhsid 
  !
  call solver(xx,unkno,amatr,dummr) 

  if( INOTEMPTY ) then
     call local_basis_global_to_local(kfl_fixrs_nsi,unkno) ! Transform global to local  
     idofn = 0
     do ipoin = 1,npoin                                    ! Prescribe bc 
        do idime = 1,ndime
           idofn = idofn + 1
           if( kfl_fixno_nsi(idime,ipoin) > 0 ) &
                unkno(idofn) = bvess_nsi(idime,ipoin,1)
        end do
     end do
     call local_basis_local_to_global(kfl_fixrs_nsi,unkno) ! Transform local to global  
  end if

  call times(4) % add() 

end subroutine nsi_elmcor
