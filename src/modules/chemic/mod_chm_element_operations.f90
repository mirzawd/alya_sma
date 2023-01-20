!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Chemic
!> @{
!> @file    mod_chm_element_operations.f90
!> @author  houzeaux
!> @date    2018-07-20
!> @brief   Element operations
!> @details Element operations
!-----------------------------------------------------------------------

module mod_chm_element_operations
    use def_kintyp,                   only : ip,rp
  implicit none

  private

  public :: chm_element_operations

contains

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-07-20
  !> @brief   Elemental operations for the Flamelet combustion model
  !> @details Elemental operations for the Flamelet combustion model
  !>
  !-----------------------------------------------------------------------

  subroutine chm_element_operations(order,pnode,pgaus,list_elements)

    use def_master,      only : cutim, solve, rhsid
    use def_domain,      only : mnode, ndime, mgaus, ntens, elmar, npoin, ltype, llapl, lorde, ltopo, lnods, hnatu
    use def_chemic,      only : nclas_chm, dt_chm, dt_rho_chm, iclaf_chm, iclai_chm, kfl_advec_chm, kfl_ellen_chm, kfl_entropy_chm,&
                                kfl_spray_chm, kfl_taust_chm, ADR_chm, ncomp_chm
    use mod_ker_proper,  only : ker_proper
    use mod_matrix,      only : matrix_assexp

    use mod_ADR,         only : ADR_element_assembly
    use mod_ADR,         only : ADR_add_sgs_or_bubble
    use mod_ADR,         only : ELEMENT_ASSEMBLY                 ! 1
    use mod_ADR,         only : mreac_adr
    use mod_chm_spray,   only : chm_elmprc_spray
    use mod_chm_spray,   only : chm_rhodt_spray
    use mod_chm_spray,   only : chm_elmpre_spray
    use mod_chm_entropy, only : chm_entropy_viscosity

    implicit none

    integer(ip), intent(in)          :: order                                    ! =1 defaul or =2 compute SGS only
    integer(ip), intent(in)          :: pnode                                    !< Number of nodes
    integer(ip), intent(in)          :: pgaus                                    !< Number of Gauss points
    integer(ip), intent(in), pointer :: list_elements(:)                         !< List of elements

    real(rp)                         :: elmat(mnode,mnode)
    real(rp)                         :: elrhs(mnode)
    integer(ip)                      :: ielem,igaus,idime,iclas,iclas_start                  ! Indices and dimensions
    integer(ip)                      :: izmat,izrhs,pelty
    integer(ip)                      :: plapl,porde,ptopo
    integer(ip)                      :: dummi,inode,kelem

    real(rp)                         :: elcon(mnode,nclas_chm,ADR_chm(1) % ntime)! <=> conce
    real(rp)                         :: elcod(ndime,mnode)                       ! <=> coord
    real(rp)                         :: elmas(mnode,nclas_chm)                   ! Mass source terms
    real(rp)                         :: elvel(ndime,mnode)
    real(rp)                         :: gpvol(mgaus)                             ! |J|*w
    real(rp)                         :: gphco(mgaus)                             ! heat conductivity
    real(rp)                         :: gpsph(mgaus)                             ! specific heat
    real(rp)                         :: gpcon(mgaus,nclas_chm,ncomp_chm)         ! <=> conce
    real(rp)                         :: gprea(mgaus,mreac_adr)                   ! r
    real(rp)                         :: gpvel(ndime,mgaus)                       ! u
    real(rp)                         :: gpdif(mgaus,nclas_chm)                   ! D_k
    real(rp)                         :: del_gpdif(mgaus,nclas_chm)               ! change in diffusivity from entropy stable
                                                                                 !  stabilization
    real(rp)                         :: gpgrd(ndime,mgaus)                       ! grad(k) = grad(D_k)
    real(rp)                         :: gprhs(mgaus)                             ! f (all terms)
    real(rp)                         :: gpden(mgaus)                             ! fake rho for elmadr
    real(rp)                         :: gptur(mgaus)                             ! turbulent viscosity
    real(rp)                         :: gppro(mgaus)                             ! Weighted residual L2-projection
    real(rp)                         :: gpdiv(mgaus)                             ! Divergence of convection
    real(rp)                         :: gpmas(mgaus,nclas_chm)                   ! Mass realease of each reaction
    real(rp)                         :: gpcar(ndime,mnode,mgaus)                 ! dNk/dxj
    real(rp)                         :: gphes(ntens,mnode,mgaus)                 ! dNk/dxidxj
    real(rp)                         :: gplap(mgaus,mnode)                       ! Laplacian
    real(rp)                         :: gpdis(mgaus,nclas_chm)                   ! Dissipation rate
    real(rp)                         :: gpprd(mgaus,nclas_chm)                   ! Production term
    real(rp)                         :: gpsigma(mgaus)                           ! RHS of liquid gas surface density equation

    real(rp)                         :: dummr(mgaus*ndime)
    real(rp)                         :: gptau(mgaus)
    real(rp)                         :: chale(3),chave(3),hleng(3),tragl(9)

    real(rp)                         :: dtmin

    external                         :: elmlen
    external                         :: elmchl
    external                         :: elmcar
    external                         :: chm_rhodt
    external                         :: chm_elmpre_flamLet
    external                         :: chm_elmprc_flamLet
    external                         :: chm_elmgac_flamLet

    !
    ! Loop over elements
    !
    dtmin = 1.0e6_rp

    !
    ! Assembly only surface density when spray with level set
    !
    if ( kfl_spray_chm == 2_ip ) then
       iclas_start = iclaf_chm
    else
       iclas_start = iclai_chm
    endif

    elements: do kelem = 1_ip,size(list_elements,kind=ip)

       ielem = list_elements(kelem)

       if( ielem > 0 ) then
          !
          ! Element dimensions
          !
          pelty = ltype(ielem)
          plapl = llapl(pelty)
          porde = lorde(pelty)
          ptopo = ltopo(pelty)

          !
          ! Initialization variables
          !
          gptau = 0.0_rp
          gpdiv = 0.0_rp
          gpdif = 0.0_rp
          gprea = 0.0_rp
          gptur = 0.0_rp
          gppro = 0.0_rp
          gpgrd = 0.0_rp
          gpcon = 0.0_rp
          gpsigma = 0.0_rp

          !
          ! Gather all
          !
          call chm_elmgac_flamLet(pnode,lnods(1,ielem),elcod,elcon(1:pnode,:,:),elvel,elmas)

          !
          ! CHALE, HLENG and TRAGL
          !
          if( kfl_taust_chm /= 0 ) then
             call elmlen(&
                  ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),&
                  hleng)
             call elmchl(&
                  tragl,hleng,elcod,dummr,chave,chale,pelty,pnode,&
                  porde,hnatu(pelty),kfl_advec_chm,kfl_ellen_chm)
          else
             plapl = 0
          end if

          !
          ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
          !
          call elmcar(&
               pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
               elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
               gphes,ielem)

          !
          ! Compute laplacian
          !
          if (plapl /= 0 ) then
             do igaus = 1,pgaus
                do inode = 1,pnode
                   gplap(igaus,inode)=0.0_rp
                   do idime=1,ndime
                      gplap(igaus,inode) = gplap(igaus,inode) + gphes(idime,inode,igaus)
                   enddo
                enddo
             enddo
          else
             gphes = 0.0_rp
             gplap = 0.0_rp
          endif

          call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,elmar(pelty)%shape,gpcar)
          call ker_proper('CONDU','PGAUS',dummi,ielem,gphco,pnode,pgaus,elmar(pelty)%shape,gpcar)
          call ker_proper('SPHEA','PGAUS',dummi,ielem,gpsph,pnode,pgaus,elmar(pelty)%shape,gpcar)
          call ker_proper('TURBU','PGAUS',dummi,ielem,gptur,pnode,pgaus,elmar(pelty)%shape,gpcar)

          !
          ! Send quantities to gauss points
          !
          call chm_elmpre_flamLet(&
               ielem,pnode,pgaus,elcon(1:pnode,:,:),elvel,elmas,elmar(pelty)%shape,&
               gpcar,gpcon(1:pgaus,:,:),gpvel,gpmas,hleng,gpdiv,gpdis(1:pgaus,:),&
               gpprd(1:pgaus,:),gptur,gpden,gphco,gpsph)

          !
          ! RHS calculation of liquid surface density at gauss points
          !
          if ( kfl_spray_chm /= 0_ip ) &
              call chm_elmpre_spray(ielem,pnode,pgaus,elvel,gpcar,gpcon(1:pgaus,:,1),&
                   gpden,gptur,hleng,gpsigma)

          !
          ! Projections of rho/dt and 1/dt for explicit
          !     Gas phase only
          !
          if ( kfl_spray_chm == 0_ip )then
               call chm_rhodt(  &
                    pnode,pgaus,porde,lnods(:,ielem),elmar(pelty)%shape,gpvol,gpden,dt_rho_chm)

          !
          ! Projections of rho/dt and 1/dt for explicit
          !     Liquid and gas phase for ELSA model
          !
          else
               call chm_rhodt_spray(  &
                    pnode,pgaus,porde,lnods(:,ielem),elmar(pelty)%shape,gpvol,gpden,dt_rho_chm,dt_chm)

          endif

          !
          ! Assemble matrix
          !
          izmat = 1
          izrhs = 1


          if(kfl_entropy_chm == 1_ip) then
             del_gpdif = 0.0_rp
             do iclas = iclai_chm,iclaf_chm
                call chm_entropy_viscosity(ielem,pnode,pgaus,1_ip,pgaus,iclas,elvel,gpden,hleng,del_gpdif)
             enddo

             !!do igaus = 1,pgaus
             !!    del_gpdif(igaus,1:nclas_chm) = maxval( del_gpdif(igaus,1:nclas_chm) )
             !!enddo
             !!print'(A,100(E18.8,X))','del_gpdif(1,1:nclas_chm)', del_gpdif(1,1:nclas_chm)

          end if






          ASSEMBLY_ICLAS: do iclas = iclas_start,iclaf_chm

             !
             ! Initialization RHS for assembly Monolithic
             !
             gprhs = 0.0_rp

             !
             ! Add sgs or bubble
             !
             call ADR_add_sgs_or_bubble(&
                  ielem,pgaus,elmar(pelty)%shape_bub,ADR_chm(iclas),gpcon(1:pgaus,iclas:iclas,:))

             !
             ! Assembly spray terms
             !
             if (kfl_spray_chm /= 0_ip .and. iclas > (nclas_chm - 2_ip) ) then
                call chm_elmprc_spray(iclas,pgaus,gptur,gpsigma,gpden,gpdif(1:pgaus,1:nclas_chm),gprhs)

             !
             ! Assembly gas phase terms
             !
             else
                call chm_elmprc_flamLet(&
                     iclas,pgaus,gpden,gpmas,gphco,gpsph,gptur,gpdis(1:pgaus,:),&
                     gpprd(1:pgaus,:),gpdif(1:pgaus,1:nclas_chm),gprhs)

             end if

             !
             ! Entropy stable viscosity
             !
             if(kfl_entropy_chm == 1_ip) then
             !!   call chm_entropy_viscosity(ielem,pnode,pgaus,1_ip,pgaus,iclas,elvel,gpden,hleng,gpdif)
                 gpdif(1:pgaus,iclas:iclas) = gpdif(1:pgaus,iclas:iclas) + del_gpdif(1:pgaus,iclas:iclas)
             end if

             !
             ! Assemble equation for iclas
             !
             if( order == ELEMENT_ASSEMBLY ) then

                call ADR_element_assembly(&
                     ielem,pnode,pgaus,elcod,elmar(pelty)%shape,gpcar,elmar(pelty)%deriv,gphes,gpvol,chale,&
                     elmar(pelty)%shape_bub,elmar(pelty)%deriv_bub,ADR_chm(iclas),cutim,gpden,gpvel,gpdif(1:pgaus,iclas:iclas),&
                     gpgrd,gprea,gprhs,gpcon(1:pgaus,iclas:iclas,:),elcon(1:pnode,iclas:iclas,:),elmat,&
                     elrhs)


                call matrix_assexp(solve(1)%ndofn,1_ip,pnode,npoin,lnods(1:pnode,ielem),elrhs,elmat, &
                                   elcon(1:pnode,iclas:iclas,1),rhsid,iclas)

                izrhs = izrhs + npoin                                       !solve(1)%nzrhs
                izmat = izmat + solve(1)%nzmat/nclas_chm**2

             end if

          end do ASSEMBLY_ICLAS

       end if

    end do elements

  end subroutine chm_element_operations

end module mod_chm_element_operations
!> @}
