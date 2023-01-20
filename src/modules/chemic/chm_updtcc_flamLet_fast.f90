!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_updtcc_flamLet_fast(VECTOR_DIM, pnode, pgaus, list_elements, dtmin)
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
  use def_master,                      only : INOTMASTER
  use def_domain,                      only : elmar, ltype, llapl, lorde, ltopo, ndime, mnode, mgaus, ntens, lnods, hnatu
  use def_chemic,                      only : nclas_chm, DtRho_gp, ADR_chm, iclaf_chm, iclai_chm, kfl_advec_chm,&
                                              kfl_DtRho_tab_index_chm, kfl_ellen_chm, kfl_spray_chm, ncomp_chm
  use mod_ker_proper,                  only : ker_proper
  use mod_ADR,                         only : ADR_critical_time_step
  use mod_ADR,                         only : mreac_adr
  use def_kintyp,                      only : ip,rp

  use mod_element_integration,         only : element_shape_function_derivatives_jacobian
  use mod_chm_element_operations_fast, only : chm_elmlen_fast, chm_elmchl_fast, chm_elmpre_flamlet_fast, chm_elmpre_spray_fast,&
                                              chm_elmprc_spray_fast, chm_elmprc_flamlet_fast, chm_elmgac_flamLet_fast

  implicit none
  integer(ip), intent(in)    :: VECTOR_DIM
  integer(ip), intent(in)    :: pnode
  integer(ip), intent(in)    :: pgaus
  integer(ip), intent(in)    :: list_elements(VECTOR_DIM)
  real(rp),    intent(inout) :: dtmin

  real(rp)     :: dtcri(2)

  integer(ip) :: ielem,iclas,iclas_start                             ! Indices and dimensions
  integer(ip) :: pelty
  integer(ip) :: plapl,porde,ptopo,dummi
  real(rp)    :: gpsha(VECTOR_DIM,pnode,pgaus)                       ! N
  real(rp)    :: gpder(VECTOR_DIM,ndime,mnode,pgaus)                 ! dNk/dsj
  real(rp)    :: elcon(VECTOR_DIM,pnode,nclas_chm,ADR_chm(1)%ntime)  ! <=> conce
  real(rp)    :: elcod(VECTOR_DIM,ndime,pnode)                       ! <=> coord
  real(rp)    :: elmas(VECTOR_DIM,pnode,nclas_chm)                   ! Mass source terms
  real(rp)    :: elvel(VECTOR_DIM,ndime,pnode)                       ! Velocity
  real(rp)    :: gpvol(VECTOR_DIM,pgaus)                             ! |J|*w
  real(rp)    :: gpcon(VECTOR_DIM,pgaus,nclas_chm,ncomp_chm)         ! <=> conce
  real(rp)    :: gprea(VECTOR_DIM,pgaus,mreac_adr,nclas_chm)         ! r
  real(rp)    :: gpvel(VECTOR_DIM,ndime,mgaus)                       ! u
  real(rp)    :: gpdif(VECTOR_DIM,pgaus,nclas_chm)                   ! D_k
  real(rp)    :: gpgrd(VECTOR_DIM,ndime,pgaus)                       ! grad(k) = grad(D_k)
  real(rp)    :: gprhs(VECTOR_DIM,pgaus)                             ! f (all terms)
  real(rp)    :: gpden(VECTOR_DIM,pgaus)                             ! rho
  real(rp)    :: gpdiv(VECTOR_DIM,pgaus)                             ! Divergence of convection
  real(rp)    :: gpmas(VECTOR_DIM,pgaus,nclas_chm)                   ! Mass realease of each reaction
  real(rp)    :: gpmascon(VECTOR_DIM,pgaus,nclas_chm)                ! Mass consumption / Y_k
  real(rp)    :: gpcar(VECTOR_DIM,ndime,mnode,mgaus)                 ! dNk/dxj
  real(rp)    :: gphes(VECTOR_DIM,ntens,mnode,mgaus)                 ! dNk/dxidxj
  real(rp)    :: gphco(VECTOR_DIM,pgaus)                             ! heat conductivity
  real(rp)    :: gpsph(VECTOR_DIM,pgaus)                             ! specific heat
  real(rp)    :: gpdtr(VECTOR_DIM,pgaus)                             ! Dt*rho
  real(rp)    :: gpprd(VECTOR_DIM,pgaus,nclas_chm)                   ! Production term of c and f equations in the Flamelet model
  real(rp)    :: gpdis(VECTOR_DIM,pgaus,nclas_chm)
  real(rp)    :: gptur(VECTOR_DIM,pgaus)
  real(rp)    :: gpsigma(VECTOR_DIM,pgaus)                           ! RHS of liquid gas surface density equation

  real(rp)    :: dummr(mgaus*ndime*mnode)
  real(rp)    :: chale(VECTOR_DIM,3),chave(VECTOR_DIM,3),hleng(VECTOR_DIM,3),tragl(VECTOR_DIM,9)

  integer(ip) :: ivect, ielemone
  integer(ip) :: list_elements_p(VECTOR_DIM)                         ! List of elements (always positive)
  integer(ip) :: lnods_loc(VECTOR_DIM,pnode)



  if( INOTMASTER ) then
     !
     ! Assembly only surface density when spray with level set
     !
     if ( kfl_spray_chm == 2_ip ) then
        iclas_start = iclaf_chm
     else
        iclas_start = iclai_chm
     endif


     !
     ! Clean list of elements with repeated indeces if necessary for padding
     !
     do ivect = 1,VECTOR_DIM
         ielem = abs(list_elements(ivect))
         if( ielem /= 0 ) then
             list_elements_p(ivect)   = list_elements(ivect)
             lnods_loc(ivect,1:pnode) = lnods(1:pnode,ielem)
         else
             list_elements_p(ivect)   = list_elements(1)
             lnods_loc(ivect,1:pnode) = lnods(1:pnode,list_elements(1))
         end if
     end do

     ielemone = list_elements_p(1)

     !
     ! Element dimensions (same for all in a vector)
     !
     pelty = ltype(ielemone)
     plapl = llapl(pelty)
     porde = lorde(pelty)
     ptopo = ltopo(pelty)

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
     gpdtr = 0.0_rp

     !
     ! Gather nodal unkown, velocity, and mass source term if available on the nodes.
     !
     call chm_elmgac_flamLet_fast(&
             VECTOR_DIM,pnode,lnods_loc,elcod,elcon,elvel,elmas)

     !
     ! Get length scales
     ! CHALE, HLENG and TRAGL
     !
     call chm_elmlen_fast(&
       VECTOR_DIM,ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),hleng)

     call chm_elmchl_fast(&
          VECTOR_DIM,tragl,hleng,elcod,dummr,chave,chale,pnode,&
          porde,hnatu(pelty),kfl_advec_chm,kfl_ellen_chm)

     !
     ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
     !
     call element_shape_function_derivatives_jacobian(&
     pnode,pgaus,plapl,elmar(pelty) % weigp,elmar(pelty) % shape,&
     elmar(pelty) % deriv,elmar(pelty) % heslo,&
     elcod,gpvol,gpsha,gpder,gpcar,gphes,list_elements=list_elements_p)

     !
     ! Density and turbulent viscosity
     !
     call ker_proper('DENSI','PGAUS',dummi,list_elements_p,gpden,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM)
     call ker_proper('TURBU','PGAUS',dummi,list_elements_p,gptur,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM)

     !
     ! Dt*rho from tabulation if available, or from Dt*rho = k/cp if it is not available
     !
     if (kfl_DtRho_tab_index_chm > 0) then
        do ivect = 1,VECTOR_DIM
           ielem = list_elements_p(ivect)
           if(ivect== 1 .or. (ivect > 1 .and. ielem /=ielemone)) then
              gpdtr(ivect,1:pgaus) = DtRho_gp(ielem) % a(1:pgaus,1,1)
           end if
        end do
     else
        call ker_proper('CONDU','PGAUS',dummi,list_elements_p,gphco,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM)
        call ker_proper('SPHEA','PGAUS',dummi,list_elements_p,gpsph,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM)
        gpdtr(VECTOR_DIM,1:pgaus) = gphco(VECTOR_DIM,1:pgaus) / gpsph(VECTOR_DIM,1:pgaus)
     endif

     !
     ! Calculate source terms on Gaussian integration points
     !
     call chm_elmpre_flamLet_fast(&
          VECTOR_DIM,pnode,pgaus,elcon,elvel,elmar(pelty)%shape,&
          gpcar,gpcon,gpvel,gpmas,hleng,gpdiv,gpdis,&
          gpprd,gptur,gpden,gpdtr,gpmascon,list_elements_p)

     !
     ! RHS calculation of liquid surface density at gauss points
     !
     if ( kfl_spray_chm /= 0_ip ) then
         call chm_elmpre_spray_fast(VECTOR_DIM,pnode,pgaus,elvel,gpcar,gpcon,&
                                    gpden,gptur,hleng,gpsigma,list_elements_p)
     end if

     !
     ! Loop over each variable
     !
     do iclas = iclas_start,iclaf_chm

        if (kfl_spray_chm /= 0_ip .and. iclas > (nclas_chm - 2_ip) ) then
           !
           ! Assembly spray terms
           !
           call chm_elmprc_spray_fast(VECTOR_DIM,iclas,pgaus,gptur,gpsigma,gpden,gpdif,gprhs(:,:))

        else
           !
           ! Assembly gas phase terms
           !
           call chm_elmprc_flamLet_fast(&
             VECTOR_DIM,iclas,pgaus,gpden,gpmas,gpmascon,gpdtr,gptur,gpdis,&
             gpprd,gpdif,gprhs(:,:),gprea(:,:,:,iclas:iclas))

        end if



        do ivect = 1,VECTOR_DIM
           ielem = list_elements_p(ivect)
           if( ivect == 1 .or. (ivect>1 .and. ielem /= ielemone)) then
              !
              ! Compute time-step
              !
              call ADR_critical_time_step(ADR_chm(iclas),gpden(ivect,:),gpvel(ivect,:,:),gpdif(ivect,1:pgaus,iclas),&
                                          gprea(ivect,:,:,iclas:iclas),dtcri,chale(ivect,1),chale(ivect,2))

              !
              ! Take minimum
              !
              dtmin = min(dtmin,dtcri(1))
           end if
        end do
     end do
  end if
end subroutine chm_updtcc_flamLet_fast
