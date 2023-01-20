!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_updtcc_flamLet(dtmin)
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_updtcc_flamLet
  ! NAME
  !    chm_updtcc_flamLet
  ! DESCRIPTION
  !    This routine computes the critical time step size in flamelet models
  ! USED BY
  !    chm_updtss
  !***
  !-----------------------------------------------------------------------
  use def_master,         only : INOTMASTER
  use def_kintyp,         only : ip, rp
  use def_domain,         only : mnode, ndime, mgaus, ntens, elmar, nelem, ltype, nnode, ngaus, llapl, lorde, ltopo, hnatu, lnods
  use def_chemic,         only : nclas_chm, ADR_chm, iclaf_chm, iclai_chm, kfl_advec_chm, kfl_ellen_chm, kfl_spray_chm, ncomp_chm
  use mod_ker_proper,     only : ker_proper
  use mod_ADR,            only : ADR_critical_time_step
  use mod_ADR,            only : mreac_adr
  use mod_communications, only : PAR_MIN

  use mod_chm_spray,      only : chm_elmprc_spray
  use mod_chm_spray,      only : chm_elmpre_spray

  implicit none
  real(rp),   intent(inout) :: dtmin
  real(rp)                  :: dtcri(2)
  integer(ip) :: ielem,iclas,iclas_start               ! Indices and dimensions
  integer(ip) :: pelty,pnode
  integer(ip) :: pgaus,plapl,porde,ptopo,dummi

  real(rp)    :: elcon(mnode,nclas_chm,ADR_chm(1)%ntime)  ! <=> conce
  real(rp)    :: elcod(ndime,mnode)                       ! <=> coord
  real(rp)    :: elmas(mnode,nclas_chm)                   ! Mass source terms
  real(rp)    :: elvel(ndime,mnode)
  real(rp)    :: gpvol(mgaus)                             ! |J|*w
  real(rp)    :: gpcon(mgaus,nclas_chm,ncomp_chm)         ! <=> conce
  real(rp)    :: gprea(mgaus,mreac_adr)                   ! r
  real(rp)    :: gpvel(ndime,mgaus)                       ! u
  real(rp)    :: gpdif(mgaus,nclas_chm)                   ! D_k
  real(rp)    :: gpgrd(ndime,mgaus)                       ! grad(k) = grad(D_k)
  real(rp)    :: gprhs(mgaus)                             ! f (all terms)
  real(rp)    :: gpden(mgaus)                             ! rho
  real(rp)    :: gpdiv(mgaus)                             ! Divergence of convection
  real(rp)    :: gpmas(mgaus,nclas_chm)                   ! Mass realease of each reaction
  real(rp)    :: gpcar(ndime,mnode,mgaus)                 ! dNk/dxj
  real(rp)    :: gphes(ntens,mnode,mgaus)                 ! dNk/dxidxj
  real(rp)    :: gphco(mgaus)                             ! heat conductivity
  real(rp)    :: gpsph(mgaus)                             ! specific heat
  real(rp)    :: gpprd(mgaus,nclas_chm)                   ! Production term of c and f equations in the Flamelet model
  real(rp)    :: gpdis(mgaus,nclas_chm)
  real(rp)    :: gptur(mgaus)
  real(rp)    :: gpsigma(mgaus)                           ! RHS of liquid gas surface density equation

  real(rp)    :: dummr(mgaus*ndime*mnode)
  real(rp)    :: chale(3),chave(3),hleng(3),tragl(9)

  external    :: elmlen
  external    :: elmchl
  external    :: elmcar
  external    :: chm_elmgac_flamLet
  external    :: chm_elmpre_flamLet
  external    :: chm_elmprc_flamLet

  if( INOTMASTER ) then

     !
     ! Assembly only surface density when spray with level set
     !
     if ( kfl_spray_chm == 2_ip ) then
        iclas_start = iclaf_chm
     else
        iclas_start = iclai_chm
     endif

     do ielem = 1,nelem
        !
        ! Initialization
        !
        gpdiv = 0.0_rp
        gprhs = 0.0_rp
        gpdif = 0.0_rp
        gprea = 0.0_rp
        gptur = 0.0_rp
        gpgrd = 0.0_rp
        gpcon = 0.0_rp
        gpsigma = 0.0_rp

        !
        ! Element dimensions
        !
        pelty = ltype(ielem)
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)
        plapl = llapl(pelty)
        porde = lorde(pelty)
        ptopo = ltopo(pelty)

        !
        ! Gather all
        !
        call chm_elmgac_flamLet(pnode,lnods(1:pnode,ielem),elcod,elcon(1:pnode,:,:),elvel,elmas)

        !
        ! CHALE, HLENG and TRAGL
        !
        call elmlen(&
             ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),hleng)

        call elmchl(&
             tragl,hleng,elcod,dummr,chave,chale,pelty,pnode,porde,hnatu(pelty),&
             kfl_advec_chm,kfl_ellen_chm)

        !
        ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
        !
        call elmcar(&
             pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
             elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
             gphes,ielem)

        !
        ! Transport properties and density
        !
        call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,elmar(pelty) % shape,gpcar)
        call ker_proper('CONDU','PGAUS',dummi,ielem,gphco,pnode,pgaus,elmar(pelty) % shape,gpcar)
        call ker_proper('SPHEA','PGAUS',dummi,ielem,gpsph,pnode,pgaus,elmar(pelty) % shape,gpcar)
        call ker_proper('TURBU','PGAUS',dummi,ielem,gptur,pnode,pgaus,elmar(pelty) % shape,gpcar)
        !
        ! Send quantities to gauss points
        !
        call chm_elmpre_flamLet(&
             ielem,pnode,pgaus,elcon(1:pnode,:,:),elvel,elmas,elmar(pelty)%shape,gpcar, &
             gpcon(1:pgaus,:,:),gpvel,gpmas,hleng,gpdiv,gpdis(1:pgaus,:),gpprd(1:pgaus,:),&
             gptur,gpden,gphco,gpsph)

        !
        ! RHS calculation of liquid surface density at gauss points
        !
        if ( kfl_spray_chm /= 0_ip ) &
              call chm_elmpre_spray(ielem,pnode,pgaus,elvel,gpcar,gpcon(1:pgaus,:,1),&
                   gpden,gptur,hleng,gpsigma)

        !
        ! Loop over variables
        !
        do iclas = iclas_start,iclaf_chm

           !
           ! Assembly spray terms
           !
           if (kfl_spray_chm /= 0_ip .and. iclas > (nclas_chm - 2_ip) ) then
              call chm_elmprc_spray(iclas,pgaus,gptur,gpsigma,gpden,gpdif,gprhs)

           !
           ! Assembly gas phase terms
           !
           else
              call chm_elmprc_flamLet(iclas,pgaus,gpden,gpmas,gphco,gpsph,gptur,gpdis(1:pgaus,:),gpprd(1:pgaus,:),gpdif,gprhs)

           end if

           !
           ! Compute time-step
           !
           call ADR_critical_time_step(ADR_chm(iclas),gpden,gpvel,gpdif(1:pgaus,iclas),gprea,dtcri,chale(1),chale(2))

           !
           ! Take minimum
           !
           dtmin = min(dtmin,dtcri(1))

        end do

     end do

  end if
  !
  ! Look for minimum over subdomains
  !
  call PAR_MIN(dtmin,'IN MY CODE')

end subroutine chm_updtcc_flamLet
