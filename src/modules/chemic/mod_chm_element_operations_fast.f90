!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Chemic
!> @{
!> @file    mod_chm_element_operations_fast.f90
!> @author  goyarzun
!> @date    2020-08-25
!> @brief   Vectorized element operations
!> @details Vectorized element operations
!-----------------------------------------------------------------------

module mod_chm_element_operations_fast
#include "def_vector_size.inc"

#define DEF_VECT 1:VECTOR_SIZE

    use def_kintyp,                   only : ip,rp
  implicit none

  private

  public :: chm_element_operations_fast
  public :: chm_elmgac_flamLet_fast
  public :: chm_elmprc_flamlet_fast
  public :: chm_elmpre_flamlet_fast
  public :: chm_elmchl_fast
  public :: chm_elmlen_fast
  public :: chm_elmpre_spray_fast
  public :: chm_elmprc_spray_fast
contains

  !-----------------------------------------------------------------------
  !>
  !> @author  goyarzun
  !> @date    2020-08-25
  !> @brief   Vectorized elemental operations for the Flamelet combustion model
  !> @details Vectorized elemental operations for the Flamelet combustion model
  !>
  !-----------------------------------------------------------------------

  subroutine chm_element_operations_fast(order,VECTOR_DIM,pnode,pgaus,list_elements)

    use def_master,                         only : tempe_gp, solve, rhsid
    use def_domain,                         only : mnode, ndime, mgaus, ntens, lnods, elmar, npoin, ltype, llapl, lorde, ltopo,&
                                                   hnatu
    use def_chemic,                         only : nclas_chm, ADR_chm, mixedEq_groups_chm, dt_chm, dt_rho_chm, iclaf_chm,&
                                                   iclai_chm, kfl_advec_chm, kfl_DtRho_tab_index_chm, kfl_ellen_chm,&
                                                   kfl_entropy_chm, kfl_soot_chm, kfl_spray_chm, kfl_taust_chm, nclas_chm,&
                                                   ngrou_chm, ncomp_chm, DtRho_gp, Yk_ssm_gp
    use mod_ker_proper,                     only : ker_proper
    use def_kermod,                         only : kfl_soot_vect

    use mod_ADR,                            only : ELEMENT_ASSEMBLY                 ! 1
    use mod_ADR,                            only : mreac_adr
    use mod_chm_spray,                      only : chm_elmprc_spray
    use mod_chm_spray,                      only : chm_rhodt_spray
    use mod_chm_spray,                      only : chm_elmpre_spray
    use mod_chm_entropy,                    only : chm_entropy_viscosity
    use mod_element_integration,            only : element_shape_function_derivatives_jacobian
    use mod_chm_sectional_soot_model_fast,  only : assembly_soot_sources_ssm
    use mod_chm_thermophoretic,             only : chm_thermophoretic_gather
    use mod_chm_thermophoretic,             only : chm_thermophoretic_elmpre
    use mod_chm_mixedEq,                    only : CHM_GR_SECTIONAL
    use mod_chm_sectional_soot_model_fast,  only : assembly_soot_sources_ssm_fast

    implicit none

    integer(ip), intent(in)          :: order                                    ! =1 defaul or =2 compute SGS only
    integer(ip), intent(in)          :: VECTOR_DIM
    integer(ip), intent(in)          :: pnode                                    !< Number of nodes
    integer(ip), intent(in)          :: pgaus                                    !< Number of Gauss points
    integer(ip), intent(in), pointer :: list_elements(:)                         !< List of elements

    real(rp)                         :: elmat(VECTOR_DIM,mnode,mnode)
    real(rp)                         :: elrhs(VECTOR_DIM,mnode)
    integer(ip)                      :: ielem,igaus,idime,iclas,iclas_start,igrou           ! Indices and dimensions
    integer(ip)                      :: izmat,izrhs,pelty
    integer(ip)                      :: plapl,porde,ptopo
    integer(ip)                      :: dummi,inode,kelem

    real(rp)                         :: elcon(VECTOR_DIM,pnode,nclas_chm,ADR_chm(1) % ntime)! <=> conce
    real(rp)                         :: elcod(VECTOR_DIM,ndime,pnode)                       ! <=> coord
    real(rp)                         :: elmas(VECTOR_DIM,pnode,nclas_chm)                   ! Mass source terms
    real(rp)                         :: elvel(VECTOR_DIM,ndime,pnode)
    real(rp)                         :: elVtherm(VECTOR_DIM,pnode,nclas_chm)                ! Thermophoretic velocity term for
                                                                                            !  sectional method
    real(rp)                         :: gpvol(VECTOR_DIM,pgaus)                             ! |J|*w
    real(rp)                         :: gphco(VECTOR_DIM,pgaus)                             ! heat conductivity
    real(rp)                         :: gpsph(VECTOR_DIM,pgaus)                             ! specific heat
    real(rp)                         :: gpdtr(VECTOR_DIM,pgaus)                             ! Dt*rho
    real(rp)                         :: gpcon(VECTOR_DIM,pgaus,nclas_chm,ncomp_chm)         ! <=> conce
    real(rp)                         :: gprea(VECTOR_DIM,pgaus,mreac_adr,nclas_chm)         ! r
    real(rp)                         :: gpvel(VECTOR_DIM,ndime,mgaus)                       ! u
    real(rp)                         :: gpdif(VECTOR_DIM,pgaus,nclas_chm)                   ! D_k
    real(rp)                         :: del_gpdif(VECTOR_DIM,pgaus,nclas_chm)               ! change in diffusivity from entropy
                                                                                            !  stable stabilization
    real(rp)                         :: gpgrd(VECTOR_DIM,ndime,pgaus)                       ! grad(k) = grad(D_k)
    real(rp)                         :: gprhs(VECTOR_DIM,pgaus,nclas_chm)                   ! f (all terms)
    real(rp)                         :: gpden(VECTOR_DIM,pgaus)                             ! fake rho for elmadr
    real(rp)                         :: gptur(VECTOR_DIM,pgaus)                             ! turbulent viscosity
    real(rp)                         :: gppro(VECTOR_DIM,mgaus)                             ! Weighted residual L2-projection ERASE?
    real(rp)                         :: gpdiv(VECTOR_DIM,pgaus)                             ! Divergence of convection
    real(rp)                         :: gpmas(VECTOR_DIM,pgaus,nclas_chm)                   ! Mass realease of each reaction
    real(rp)                         :: gpmascon(VECTOR_DIM,pgaus,nclas_chm)                ! Mass consumption term / Y_k
    real(rp)                         :: gpcar(VECTOR_DIM,ndime,mnode,mgaus)                 ! dNk/dxj
    real(rp)                         :: gphes(VECTOR_DIM,ntens,mnode,mgaus)                 ! dNk/dxidxj
    real(rp)                         :: gplap(VECTOR_DIM,pgaus,pnode)                       ! Laplacian
    real(rp)                         :: gpdis(VECTOR_DIM,pgaus,nclas_chm)                   ! Dissipation rate
    real(rp)                         :: gpprd(VECTOR_DIM,pgaus,nclas_chm)                   ! Production term
    real(rp)                         :: gpsigma(VECTOR_DIM,pgaus)                           ! RHS of liquid gas surface density
                                                                                            !  equation
    real(rp)                         :: gpdivVt(VECTOR_DIM,pgaus,nclas_chm)                 ! Thermophoretic velocity term for
                                                                                            !  sectional method


    real(rp)                         :: tempe_gp_vec(VECTOR_DIM, pgaus)
    real(rp)                         :: dummr(VECTOR_DIM,mgaus*ndime) !To erase?
    real(rp)                         :: gptau(VECTOR_DIM,mgaus)       !(To erase?)
    real(rp)                         :: chale(VECTOR_DIM,3),chave(VECTOR_DIM,3),hleng(VECTOR_DIM,ndime)
    real(rp)                         :: tragl(VECTOR_DIM,ndime,ndime)

    real(rp)                         :: dtmin

    integer(ip) :: list_elements_p(VECTOR_DIM)           ! List of elements (always positive)
    integer(ip) :: lnods_loc(VECTOR_DIM,pnode)
    integer(ip) :: ivect , ielemone, iltest
    real(rp)    :: gpsha(VECTOR_DIM,pnode,pgaus)                    ! N
    real(rp)    :: gpder(VECTOR_DIM,ndime,mnode,pgaus)              ! dNk/dsj

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
    !
    ! Element dimensions
    !
    ielem = list_elements(1)

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
    gptur = 0.0_rp
    gppro = 0.0_rp
    gpgrd = 0.0_rp
    gpcon = 0.0_rp
    gpsigma = 0.0_rp
    gpdivVt = 0.0_rp

    ielemone= ielem !! Just to avoid repeating terms in pseudo vectorized version

    do ivect = 1,VECTOR_DIM
        ielem = abs(list_elements(ivect))
        if( ielem /= 0 ) then
            list_elements_p(ivect)   = list_elements(ivect)
        else
            list_elements_p(ivect)   = list_elements(1)
        end if
    end do

    do ivect = 1, VECTOR_DIM
      ielem = abs(list_elements(ivect))
      if( ielem /= 0 ) then
         lnods_loc(ivect,1:pnode) = lnods(1:pnode,ielem)
         ielem                    = list_elements(ivect)
      else
         lnods_loc(ivect,1:pnode) = lnods(1:pnode,list_elements(1))
         ielem                    = list_elements(1)
      end if
    end do

    !
    ! Gather all
    !
    call chm_elmgac_flamLet_fast(&
              VECTOR_DIM,pnode,lnods_loc,elcod,elcon,elvel,elmas)

    !
    ! Gather thermophoretic velocity term
    !
    !
    ! TODO PLEASE VECTORIZE THIS
    !
    do igrou = 1,ngrou_chm
       if (mixedEq_groups_chm(igrou) % kfl_therm_phor /= 0) then
          do ivect = 1, VECTOR_DIM
             ielem = abs(list_elements(ivect))
             if( ielem /= 0 ) then
                ielem                    = list_elements(ivect)
             else
                ielem                    = list_elements(1)
             end if
             call chm_thermophoretic_gather(pnode,                               &
                 &                          mixedEq_groups_chm(igrou) % nequa,   &
                 &                          mixedEq_groups_chm(igrou) % i_start, &
                 &                          lnods(1:pnode,ielem),                &
                 &                          elVtherm(ivect,1:pnode,mixedEq_groups_chm(igrou)%i_start:&
                                            mixedEq_groups_chm(igrou)%i_end))
          end do
       endif
    enddo

    !
    ! CHALE, HLENG and TRAGL
    !
    if( kfl_taust_chm /= 0 ) then
        call chm_elmlen_fast(&
          VECTOR_DIM,ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),hleng)

        call chm_elmchl_fast(&
             VECTOR_DIM,tragl,hleng,elcod,dummr,chave,chale,pnode,&
             porde,hnatu(pelty),kfl_advec_chm,kfl_ellen_chm)
    else
        plapl = 0
    end if
    !
    ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
    !
    call element_shape_function_derivatives_jacobian(&
       pnode,pgaus,plapl,elmar(pelty) % weigp,elmar(pelty) % shape,&
       elmar(pelty) % deriv,elmar(pelty) % heslo,&
       elcod,gpvol,gpsha,gpder,gpcar,gphes,list_elements=list_elements_p)

    !
    ! Compute laplacian
    !
    if (plapl /= 0 ) then
       do igaus = 1,pgaus
          do inode = 1,pnode
             gplap(:,igaus,inode)=0.0_rp
             do idime=1,ndime
                gplap(:,igaus,inode) = gplap(:,igaus,inode) + gphes(:,idime,inode,igaus)
             enddo
          enddo
       enddo
    else
       gphes = 0.0_rp
       gplap = 0.0_rp
    endif

    call ker_proper('DENSI','PGAUS',dummi,list_elements_p,gpden,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM)
    call ker_proper('TURBU','PGAUS',dummi,list_elements_p,gptur,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM)

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
       gpdtr(:,1:pgaus) = gphco(:,1:pgaus) / gpsph(:,1:pgaus)
    endif

    !
    ! Send quantities to Gauss points
    !
    call chm_elmpre_flamLet_fast(&
         VECTOR_DIM,pnode,pgaus,elcon,elvel,elmar(pelty)%shape,&
         gpcar,gpcon,gpvel,gpmas,hleng,gpdiv,gpdis,&
         gpprd,gptur,gpden,gpdtr,gpmascon,list_elements_p)

    !
    ! Thermophoretic velocity term of Gauss points
    !
    !
    ! TODO PLEASE VECTORIZE THIS
    !
    do igrou = 1,ngrou_chm
       if (mixedEq_groups_chm(igrou) % kfl_therm_phor /= 0) then
          do ivect = 1, VECTOR_DIM
            ielem = abs(list_elements(ivect))
            if( ielem /= 0 ) then
               ielem                    = list_elements(ivect)
            else
               ielem                    = list_elements(1)
            end if
            call chm_thermophoretic_elmpre(pnode,                               &
                &                          pgaus,                               &
                &                          mixedEq_groups_chm(igrou) % nequa,   &
                &                          elmar(pelty)%shape,                  &
                &                          elVtherm(ivect,1:pnode,mixedEq_groups_chm(igrou)%i_start:&
                                           mixedEq_groups_chm(igrou)%i_end),gpdivVt(ivect,1:pgaus,mixedEq_groups_chm(&
                                           igrou)%i_start:mixedEq_groups_chm(igrou)%i_end))
          end do
       endif
    enddo

    !
    ! RHS calculation of liquid surface density at gauss points
    !
    if ( kfl_spray_chm /= 0_ip ) &
        call chm_elmpre_spray_fast(VECTOR_DIM,pnode,pgaus,elvel,gpcar,gpcon,&
             gpden,gptur,hleng,gpsigma,list_elements_p)

    !
    ! Projections of rho/dt and 1/dt for explicit
    !     Gas phase only
    !
    if ( kfl_spray_chm == 0_ip )then
         call chm_rhodt_fast(  &
              VECTOR_DIM,pnode,pgaus,porde,lnods_loc,elmar(pelty)%shape,gpvol,gpden,dt_rho_chm,list_elements_p)
    else
    !
    ! Projections of rho/dt and 1/dt for explicit
    !     Liquid and gas phase for ELSA model
    !
       call chm_rhodt_spray_fast(  &
              VECTOR_DIM,pnode,pgaus,porde,lnods_loc,elmar(pelty)%shape,gpvol,gpden,dt_rho_chm,dt_chm,list_elements_p)
    end if



    testing4: do kelem = 1_ip,size(list_elements,kind=ip)

     ielem = list_elements(kelem)
     iltest = list_elements(kelem)

     if( ielem > 0 ) then

        ielem = kelem
        !
        ! Assemble matrix
        !
        izmat = 1
        izrhs = 1


        if(kfl_entropy_chm == 1_ip) then
           del_gpdif = 0.0_rp
           do iclas = iclai_chm,iclaf_chm
              call chm_entropy_viscosity(iltest,pnode,pgaus,1_ip,pgaus,iclas,elvel(ielem,:,:),gpden(ielem,:),hleng(ielem,:),&
                  del_gpdif(ielem,:,:))
           enddo
        end if
      end if
    end do testing4
    !
    ! Initialization RHS for assembly Monolithic
    !
    gprhs = 0.0_rp
    gprea = 0.0_rp


    if (kfl_soot_chm /= 0) then
       !
       ! Group based source treatment
       !
       do igrou = 1,ngrou_chm
          !
          ! SECTIONAL SOOT MODEL
          !
          if (mixedEq_groups_chm(igrou) % kfl_grtype == CHM_GR_SECTIONAL) then
             !
             ! TODO PLEASE VECTORIZE THIS
             !
             tempe_gp_vec = 0.0d0
!             Yk_ssm_gp_vec = 0.0d0

             do ivect = 1, VECTOR_DIM
               ielem = abs(list_elements(ivect))
               if( ielem /= 0 ) then
                  ielem                    = list_elements(ivect)
               else
                  ielem                    = list_elements(1)
               end if

               tempe_gp_vec(ivect, 1:pgaus) = tempe_gp(ielem) % a(1:pgaus,1,1)
!               Yk_ssm_gp_vec(ivect,1:pgaus,:) = Yk_ssm_gp(ielem) % a(1:pgaus,:,1)

             end do

             if( kfl_soot_vect /= 0_ip) then

                 call assembly_soot_sources_ssm_fast(list_elements_p,pgaus,&
                                                  gpcon(:,1:pgaus,1:nclas_chm,:),&
                                                  tempe_gp_vec(:,1:pgaus),&
                                                  gpden(:,1:pgaus),&
                                                  gprhs(:,1:pgaus,1:nclas_chm))
             else
                 do ivect = 1, VECTOR_DIM
                   ielem = abs(list_elements(ivect))
                   if( ielem /= 0 ) then
                      ielem                    = list_elements(ivect)
                   else
                      ielem                    = list_elements(1)
                   end if

                   call assembly_soot_sources_ssm(ielem,pgaus,&
                                                  gpcon(ivect,1:pgaus,1:nclas_chm,:),&
                                                  tempe_gp(ielem) % a(1:pgaus,1,1),&
                                                  gpden(ivect,1:pgaus),&
                                                  gprhs(ivect,1:pgaus,1:nclas_chm),&
                                                  Yk_ssm_gp(ielem) % a(1:pgaus,:,1))
                 end do
             end if
          endif
       enddo
    endif


    do igrou = 1,ngrou_chm
       if (mixedEq_groups_chm(igrou) % kfl_therm_phor /= 0) then
          !
          !  Assembly thermophoresis term on the RHS
          !
          gprhs(:,1:pgaus,mixedEq_groups_chm(igrou)%i_start:mixedEq_groups_chm(igrou)%i_end) =        &
              &   gprhs(:,1:pgaus,mixedEq_groups_chm(igrou)%i_start:mixedEq_groups_chm(igrou)%i_end)  &
              & - gpdivVt(:,1:pgaus,mixedEq_groups_chm(igrou)%i_start:mixedEq_groups_chm(igrou)%i_end)
       endif
    enddo


    ASSEMBLY_ICLAS0: do iclas = iclas_start,iclaf_chm

       !
       ! Add sgs or bubble
       !
       call chm_ADR_add_sgs_or_bubble_fast(&
            VECTOR_DIM,iclas,pgaus,elmar(pelty)%shape_bub,ADR_chm(iclas),gpcon,list_elements_p)

       !
       ! Assembly spray terms
       !
       if (kfl_spray_chm /= 0_ip .and. iclas > (nclas_chm - 2_ip) ) then
          call chm_elmprc_spray_fast(VECTOR_DIM,iclas,pgaus,gptur,gpsigma,gpden,gpdif,gprhs(:,:,iclas:iclas))
       !
       ! Assembly gas phase terms
       !
       else
          call chm_elmprc_flamLet_fast(&
               VECTOR_DIM,iclas,pgaus,gpden,gpmas,gpmascon,gpdtr,gptur,gpdis,&
               gpprd,gpdif,gprhs(:,:,iclas:iclas),gprea(:,:,:,iclas:iclas))
       end if

       !
       ! Entropy stable viscosity
       !
       if(kfl_entropy_chm == 1_ip) then
       !! call chm_entropy_viscosity(ielem,pnode,pgaus,1_ip,pgaus,iclas,elvel,gpden(ielem,:),hleng,gpdif)
         gpdif(:,1:pgaus,iclas:iclas) = gpdif(:,1:pgaus,iclas:iclas) + del_gpdif(:,1:pgaus,iclas:iclas)
       end if

    end do ASSEMBLY_ICLAS0


       !
       ! Assemble equation for iclas
       !
       if( order == ELEMENT_ASSEMBLY ) then
       ASSEMBLY_ICLAS1: do iclas = iclas_start,iclaf_chm


           call chm_element_assembly(&
                VECTOR_DIM,pnode,pgaus,list_elements_p,elmar(pelty)%shape,gpcar,&
                gphes,gpvol,chale,ADR_chm(iclas),&
                gpden,gpvel,gpdif(:,:,iclas:iclas),gpgrd,gprea(:,:,:,iclas),gprhs(:,:,iclas),&
                gpcon(:,:,iclas:iclas,:),elmat,elrhs)

!@      end do ASSEMBLY_ICLAS1


!@      ASSEMBLY_ICLAS: do iclas = iclas_start,iclaf_chm


          call chm_matrix_assexp_fast(VECTOR_DIM,solve(1)%ndofn,1_ip,pnode,lnods_loc,elrhs,elmat, &
                                    elcon(:,:,iclas:iclas,1),rhsid,iclas,list_elements_p)

          elements: do kelem = 1_ip,size(list_elements,kind=ip)
            ielem = list_elements(kelem)
            if( ielem > 0 ) then
                  izrhs = izrhs + npoin                                       !solve(1)%nzrhs
                  izmat = izmat + solve(1)%nzmat/nclas_chm**2
              end if
          end do elements


!        end do ASSEMBLY_ICLAS
        end do ASSEMBLY_ICLAS1

       end if


  end subroutine chm_element_operations_fast

  subroutine chm_elmgac_flamLet_fast(&
     VECTOR_DIM,pnode,lnods_loc,elcod,elcon,elvel,elmas)
  !------------------------------------------------------------------------
  !****f* Chemic/chm_elmgac_cfi
  ! NAME
  !    chm_elmgac_cfi
  ! DESCRIPTION
  !    Gather operations for the combustion models
  ! USES
  ! USED BY
  !    chm_elmcfi
  !***
  !------------------------------------------------------------------------
  use def_domain, only     :  ndime,coord
  use def_master, only     :  conce,massk,advec
  use def_chemic, only     :  kfl_advec_chm,&
       &                      nclas_chm,&
       &                      ADR_chm,&
       &                      kfl_lookg_chm
  use mod_ADR,    only     :  BDF

  implicit none
  integer(ip), intent(in)  :: VECTOR_DIM,pnode
  integer(ip), intent(in)  :: lnods_loc(VECTOR_DIM,pnode)
  real(rp),    intent(out) :: elcod(VECTOR_DIM,ndime,pnode)
  real(rp),    intent(out) :: elcon(VECTOR_DIM,pnode,nclas_chm,ADR_chm(1) % ntime)
  real(rp),    intent(out) :: elvel(VECTOR_DIM,ndime,pnode)
  real(rp),    intent(out) :: elmas(VECTOR_DIM,pnode,nclas_chm)                ! Mass source terms
  integer(ip)              :: inode,ipoin,idime,itime,iclas
  integer(ip)              :: ivect
  !
  ! Concentration and coordinates
  !

  do ivect=1, VECTOR_DIM
     do inode=1,pnode
         ipoin=lnods_loc(ivect,inode)
        do iclas=1,nclas_chm
            elcon(ivect,inode,iclas,1) = conce(ipoin,iclas,1)
        end do
        do idime=1,ndime
            elcod(ivect,idime,inode)   = coord(idime,ipoin)
        end do
    end do
  end do

  !
  ! Time integration
  !
  if( ADR_chm(1) % kfl_time_integration == 1 ) then

     do ivect = 1,VECTOR_DIM
        do iclas = 1,nclas_chm
            do inode = 1,pnode
            ipoin = lnods_loc(ivect,inode)
            elcon(ivect,inode,iclas,2) = conce(ipoin,iclas,3)
            end do
        end do
     end do

     if( ADR_chm(1) % kfl_time_scheme == BDF ) then

      do ivect = 1,VECTOR_DIM
        do iclas = 1,nclas_chm
           do itime = 3,ADR_chm(iclas) % kfl_time_order + 1
              do inode = 1,pnode
                 ipoin = lnods_loc(ivect,inode)
                 elcon(ivect,inode,iclas,itime) = conce(ipoin,iclas,itime+1)
              end do
           end do
        end do
      end do
     end if

  end if

  !
  ! Advection
  !
  if( kfl_advec_chm /= 0 ) then

    do ivect =1,VECTOR_DIM
     do inode = 1,pnode
        ipoin = lnods_loc(ivect,inode)
        do idime = 1,ndime
           elvel(ivect,idime,inode) = advec(idime,ipoin,1)
        end do
     end do
    end do

  else
     elvel = 0.0_rp
  end if


  !
  ! Mass source terms coefficients
  !
  elmas = 0.0_rp
  if (kfl_lookg_chm == 0) then
     do ivect = 1, VECTOR_DIM
        do inode = 1,pnode
          ipoin = lnods_loc(ivect,inode)
           do iclas = 1, nclas_chm
             elmas(ivect,inode,iclas) = massk(ipoin,iclas)
           enddo
       end do
     end do
  endif

end subroutine chm_elmgac_flamLet_fast


    subroutine chm_elmlen_fast(VECTOR_DIM,ndime,pnode,deriv,tragl,elcod,hnatu,hleng)
      implicit none
      integer(ip), intent(in)  :: VECTOR_DIM                                       !< Number of nodes
      integer(ip), intent(in)  :: ndime                                            !<
      integer(ip), intent(in)  :: pnode                                            !< Number of nodes
      real(rp),    intent(out) :: tragl(VECTOR_DIM,ndime,ndime)
      real(rp),    intent(in)  :: hnatu,elcod(VECTOR_DIM,ndime,pnode)
      real(rp),    intent(in)  :: deriv(ndime,pnode)
      real(rp),    intent(out) :: hleng(VECTOR_DIM,ndime)


!  integer(ip), intent(in)  :: ndime,pnode
!  real(rp),    intent(out) :: tragl(VECTOR_DIM,ndime,ndime)
!  real(rp),    intent(in)  :: hnatu,elcod(VECTOR_DIM,ndime,pnode)
!  real(rp),    intent(in)  :: deriv(ndime,pnode)
!  real(rp),    intent(out) :: hleng(VECTOR_DIM,ndime)

  integer(ip)              :: k,ivect
  real(rp)                 :: enor0(VECTOR_DIM),h_tem,gpdet(VECTOR_DIM),denom(VECTOR_DIM)
  real(rp)                 :: xjacm(VECTOR_DIM,ndime,ndime),t1(VECTOR_DIM),t2(VECTOR_DIM),t3(VECTOR_DIM)

!        call elmlen(ndime,pnode,elmar(pelty) % dercg,tragl,elcod,&
!             hnatu(pelty),hleng)

      if( ndime == 2 ) then

         xjacm(DEF_VECT,1,1) = 0.0_rp
         xjacm(DEF_VECT,1,2) = 0.0_rp
         xjacm(DEF_VECT,2,1) = 0.0_rp
         xjacm(DEF_VECT,2,2) = 0.0_rp
         do k = 1,pnode
            xjacm(DEF_VECT,1,1) = xjacm(DEF_VECT,1,1) + elcod(DEF_VECT,1,k) * deriv(1,k)
            xjacm(DEF_VECT,1,2) = xjacm(DEF_VECT,1,2) + elcod(DEF_VECT,1,k) * deriv(2,k)
            xjacm(DEF_VECT,2,1) = xjacm(DEF_VECT,2,1) + elcod(DEF_VECT,2,k) * deriv(1,k)
            xjacm(DEF_VECT,2,2) = xjacm(DEF_VECT,2,2) + elcod(DEF_VECT,2,k) * deriv(2,k)
         end do

         gpdet(DEF_VECT)      =  xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,2,2) - xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,1,2)
         denom(DEF_VECT)      =  1.0_rp/gpdet(DEF_VECT)
         tragl(DEF_VECT,1,1) =  xjacm(DEF_VECT,2,2) * denom(DEF_VECT)
         tragl(DEF_VECT,2,2) =  xjacm(DEF_VECT,1,1) * denom(DEF_VECT)
         tragl(DEF_VECT,2,1) = -xjacm(DEF_VECT,2,1) * denom(DEF_VECT)
         tragl(DEF_VECT,1,2) = -xjacm(DEF_VECT,1,2) * denom(DEF_VECT)

         enor0(DEF_VECT)      =  tragl(DEF_VECT,1,1) * tragl(DEF_VECT,1,1) + tragl(DEF_VECT,1,2) * tragl(DEF_VECT,1,2)
         hleng(DEF_VECT,1)   =  hnatu/sqrt(enor0(DEF_VECT))
         enor0(DEF_VECT)      =  tragl(DEF_VECT,2,1) * tragl(DEF_VECT,2,1) + tragl(DEF_VECT,2,2) * tragl(DEF_VECT,2,2)
         hleng(DEF_VECT,2)   =  hnatu/sqrt(enor0(DEF_VECT))

         do ivect = 1,VECTOR_DIM
          if( hleng(ivect,2) > hleng(ivect,1) ) then
             h_tem    = hleng(ivect,2)
             hleng(ivect,2) = hleng(ivect,1)
             hleng(ivect,1) = h_tem
          end if
         end do
        else if( ndime == 3 ) then

          xjacm(DEF_VECT,1,1) = 0.0_rp ! xjacm = elcod * deriv^t
          xjacm(DEF_VECT,1,2) = 0.0_rp ! tragl = xjacm^-1
          xjacm(DEF_VECT,1,3) = 0.0_rp
          xjacm(DEF_VECT,2,1) = 0.0_rp
          xjacm(DEF_VECT,2,2) = 0.0_rp
          xjacm(DEF_VECT,2,3) = 0.0_rp
          xjacm(DEF_VECT,3,1) = 0.0_rp
          xjacm(DEF_VECT,3,2) = 0.0_rp
          xjacm(DEF_VECT,3,3) = 0.0_rp
          do k = 1,pnode
             xjacm(DEF_VECT,1,1) = xjacm(DEF_VECT,1,1) + elcod(DEF_VECT,1,k) * deriv(1,k)
             xjacm(DEF_VECT,1,2) = xjacm(DEF_VECT,1,2) + elcod(DEF_VECT,1,k) * deriv(2,k)
             xjacm(DEF_VECT,1,3) = xjacm(DEF_VECT,1,3) + elcod(DEF_VECT,1,k) * deriv(3,k)
             xjacm(DEF_VECT,2,1) = xjacm(DEF_VECT,2,1) + elcod(DEF_VECT,2,k) * deriv(1,k)
             xjacm(DEF_VECT,2,2) = xjacm(DEF_VECT,2,2) + elcod(DEF_VECT,2,k) * deriv(2,k)
             xjacm(DEF_VECT,2,3) = xjacm(DEF_VECT,2,3) + elcod(DEF_VECT,2,k) * deriv(3,k)
             xjacm(DEF_VECT,3,1) = xjacm(DEF_VECT,3,1) + elcod(DEF_VECT,3,k) * deriv(1,k)
             xjacm(DEF_VECT,3,2) = xjacm(DEF_VECT,3,2) + elcod(DEF_VECT,3,k) * deriv(2,k)
             xjacm(DEF_VECT,3,3) = xjacm(DEF_VECT,3,3) + elcod(DEF_VECT,3,k) * deriv(3,k)
          end do

          t1(DEF_VECT)         =  xjacm(DEF_VECT,2,2) * xjacm(DEF_VECT,3,3) - xjacm(DEF_VECT,3,2) * xjacm(DEF_VECT,2,3)
          t2(DEF_VECT)         = -xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,3,3) + xjacm(DEF_VECT,3,1) * xjacm(DEF_VECT,2,3)
          t3(DEF_VECT)         =  xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,3,2) - xjacm(DEF_VECT,3,1) * xjacm(DEF_VECT,2,2)

          gpdet(DEF_VECT)      =  xjacm(DEF_VECT,1,1) * t1(DEF_VECT) + xjacm(DEF_VECT,1,2) * t2(DEF_VECT) + xjacm(DEF_VECT,1,3) * t3(DEF_VECT)
          denom(DEF_VECT)      =  1.0_rp / gpdet(DEF_VECT)
          tragl(DEF_VECT,1,1) =  t1(DEF_VECT) * denom(DEF_VECT)
          tragl(DEF_VECT,2,1) =  t2(DEF_VECT) * denom(DEF_VECT)
          tragl(DEF_VECT,3,1) =  t3(DEF_VECT) * denom(DEF_VECT)
          tragl(DEF_VECT,2,2) =  ( xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,3,3) - xjacm(DEF_VECT,3,1) * xjacm(DEF_VECT,1,3)) * denom(DEF_VECT)
          tragl(DEF_VECT,3,2) =  (-xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,3,2) + xjacm(DEF_VECT,1,2) * xjacm(DEF_VECT,3,1)) * denom(DEF_VECT)
          tragl(DEF_VECT,3,3) =  ( xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,2,2) - xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,1,2)) * denom(DEF_VECT)
          tragl(DEF_VECT,1,2) =  (-xjacm(DEF_VECT,1,2) * xjacm(DEF_VECT,3,3) + xjacm(DEF_VECT,3,2) * xjacm(DEF_VECT,1,3)) * denom(DEF_VECT)
          tragl(DEF_VECT,1,3) =  ( xjacm(DEF_VECT,1,2) * xjacm(DEF_VECT,2,3) - xjacm(DEF_VECT,2,2) * xjacm(DEF_VECT,1,3)) * denom(DEF_VECT)
          tragl(DEF_VECT,2,3) =  (-xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,2,3) + xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,1,3)) * denom(DEF_VECT)
           !
           ! Element length HLENG
           !
           enor0(DEF_VECT)    = tragl(DEF_VECT,1,1) * tragl(DEF_VECT,1,1) &
                 &     + tragl(DEF_VECT,1,2) * tragl(DEF_VECT,1,2) &
                 &     + tragl(DEF_VECT,1,3) * tragl(DEF_VECT,1,3)
           hleng(DEF_VECT,1)  = hnatu/sqrt(enor0(DEF_VECT))
           enor0(DEF_VECT)    = tragl(DEF_VECT,2,1) * tragl(DEF_VECT,2,1) &
                 &     + tragl(DEF_VECT,2,2) * tragl(DEF_VECT,2,2) &
                 &     + tragl(DEF_VECT,2,3) * tragl(DEF_VECT,2,3)
           hleng(DEF_VECT,2)  = hnatu/sqrt(enor0(DEF_VECT))
           enor0(DEF_VECT)    = tragl(DEF_VECT,3,1) * tragl(DEF_VECT,3,1) &
                 &     + tragl(DEF_VECT,3,2) * tragl(DEF_VECT,3,2) &
                 &     + tragl(DEF_VECT,3,3) * tragl(DEF_VECT,3,3)
           hleng(DEF_VECT,3)  = hnatu/sqrt(enor0(DEF_VECT))
           !
           ! Sort hleng: hleng(1)=max; hleng(ndime)=min
           !
          do ivect = 1,VECTOR_DIM

           if( hleng(ivect,2) > hleng(ivect,1) ) then
              h_tem    = hleng(ivect,2)
              hleng(ivect,2) = hleng(ivect,1)
              hleng(ivect,1) = h_tem
           end if
           if( hleng(ivect,3) > hleng(ivect,1) ) then
              h_tem    = hleng(ivect,3)
              hleng(ivect,3) = hleng(ivect,1)
              hleng(ivect,1) = h_tem
           end if
           if( hleng(ivect,3) > hleng(ivect,2) ) then
              h_tem    = hleng(ivect,3)
              hleng(ivect,3) = hleng(ivect,2)
              hleng(ivect,2) = h_tem
           end if
         end do

       end if


    end subroutine chm_elmlen_fast



    subroutine chm_elmchl_fast(VECTOR_DIM, &
     tragl,hleng,elcod,elvel,chave,chale,pnode,porde,&
     hnatu,kfl_advec,kfl_ellen)
      !-----------------------------------------------------------------------
      !****f* Domain/elmchl
      ! NAME
      !   elmchl
      ! DESCRIPTION
      !   This routine computes the characteristic element lengths CHALE
      !   according to a given strategy. CHALE is divided by two for
      !   quadratic elements:
      !   KFL_ELLEN = 0 ... CHALE(1) = Minimum element length
      !                 ... CHALE(2) = Minimum element length
      !   KFL_ELLEN = 1 ... CHALE(1) = Maximum element length
      !                 ... CHALE(2) = Maximum element length
      !   KFL_ELLEN = 2 ... CHALE(1) = Average element length
      !                 ... CHALE(2) = Average element length
      !   KFL_ELLEN = 3 ... IF KFL_ADVEC = 1:
      !                     CHALE(1) = Flow direction
      !                     CHALE(2) = Flow direction
      !                     ELSE IF KFL_ADVEC =0:
      !                     CHALE(1) = Minimum element length
      !                     CHALE(2) = Minimum element length
      !   KFL_ELLEN = 4 ... CHALE(1) = Approx. diameter=sqrt(hmin*hmax)
      !                 ... CHALE(2) = Approx. diameter=sqrt(hmin*hmax)
      !   KFL_ELLEN = 5 ... CHALE(1) = Length in flow direction
      !                 ... CHALE(2) = Minimum element kength
      ! OUTPUT
      !   CHALE
      ! USES
      ! USED BY
      !***
      !-----------------------------------------------------------------------
      use def_domain, only     :  ndime
      implicit none
      integer(ip), intent(in)  :: VECTOR_DIM,pnode,porde,kfl_advec,kfl_ellen
      real(rp),    intent(in)  :: hnatu
      real(rp),    intent(out) :: chale(VECTOR_DIM,2)
      real(rp),    intent(in)  :: tragl(VECTOR_DIM,ndime,ndime),hleng(VECTOR_DIM,ndime)
      real(rp),    intent(in)  :: elcod(VECTOR_DIM,ndime,pnode)
      real(rp),    intent(in)  :: elvel(VECTOR_DIM,ndime,pnode)
      real(rp),    intent(out) :: chave(VECTOR_DIM,ndime,2)
      integer(ip)              :: idime,inode,ivect
      real(rp)                 :: elno1,elno2

      external                 :: mbvab1
      external                 :: velchl

    if(kfl_ellen==0) then
       !
       ! Minimum element length
       !
       chale(DEF_VECT,1)=hleng(DEF_VECT,ndime)
       chale(DEF_VECT,2)=chale(DEF_VECT,1)

    else if(kfl_ellen==1) then
       !
       ! Maximum element length
       !
       chale(DEF_VECT,1)=hleng(DEF_VECT,1)
       chale(DEF_VECT,2)=chale(DEF_VECT,1)

    else if(kfl_ellen==2) then
       !
       ! Average length
       !
       chale(DEF_VECT,1)=0.0_rp
       do idime=1,ndime
          chale(DEF_VECT,1)=chale(DEF_VECT,1)+hleng(DEF_VECT,idime)
       end do
       chale(DEF_VECT,1)=chale(DEF_VECT,1)/real(ndime,rp)
       chale(DEF_VECT,2)=chale(DEF_VECT,1)

    else if(kfl_ellen==3) then
       !
       ! Length in flow direction
       !
       if(kfl_advec/=0) then
          !
          ! Characteristic element velocity (average)
          !
          chave=0.0_rp
          do idime=1,ndime
             do inode=1,pnode
                chave(DEF_VECT,idime,1)=chave(DEF_VECT,idime,1)+elvel(DEF_VECT,idime,inode)
             end do
             chave(DEF_VECT,idime,1)=chave(DEF_VECT,idime,1)/real(pnode,rp)
          end do
          !
          ! Characteristic element length u^l = J^(-t) u^g
          !
          do ivect=1,VECTOR_DIM
            call mbvab1(chave(ivect,1,2),tragl(ivect,:,:),chave(ivect,1,1),ndime,ndime,elno2,elno1)
            if(elno2>1.0e-16_rp.and.elno1>1.0e-16_rp) then
               chale(ivect,1)=hnatu*elno1/elno2
            else
               chale(ivect,1)=hleng(ivect,ndime)
            end if
          end do

          chale(DEF_VECT,2)=chale(DEF_VECT,1)
          chale(DEF_VECT,2)=hleng(DEF_VECT,ndime)

          if (ndime ==3 ) then
             chale(DEF_VECT,2)=(hleng(DEF_VECT,ndime)*hleng(DEF_VECT,2)*hleng(DEF_VECT,1))**(1.0_rp/3.0_rp)
          else if (ndime==2) then
             chale(DEF_VECT,2)=sqrt(hleng(DEF_VECT,2)*hleng(DEF_VECT,1))
          end if
       else
          chale(DEF_VECT,1)=hleng(DEF_VECT,ndime)
          chale(DEF_VECT,2)=chale(DEF_VECT,1)
       end if

    else if(kfl_ellen==4) then
       !
       ! sqrt(hmin*hmax)
       !
       chale(DEF_VECT,1)=sqrt(hleng(DEF_VECT,1)*hleng(DEF_VECT,ndime))
       chale(DEF_VECT,2)=chale(DEF_VECT,1)

    else if(kfl_ellen==5) then
       !
       ! Along velocity direction
       !
       do ivect=1,VECTOR_DIM
         call velchl(pnode,elcod(ivect,:,:),elvel(ivect,:,:),chale(ivect,:),hleng(ivect,:))
       end do

    else if(kfl_ellen==6) then
       !
       ! Mixed element length - hmin for tau1, hmax for tau2 - here we only obtain the values for tau1 - tau2 directly in nsi_elmsgs
       !
       chale(DEF_VECT,1)=hleng(DEF_VECT,ndime)
       chale(DEF_VECT,2)=chale(DEF_VECT,1)

    end if
    !
    ! Divide h by 2 for quadratic elements and 3 for cubic elements
    !
    chale(DEF_VECT,1) = chale(DEF_VECT,1)/real(porde,rp)
    chale(DEF_VECT,2) = chale(DEF_VECT,2)/real(porde,rp)

  end subroutine chm_elmchl_fast

  subroutine chm_elmpre_flamLet_fast(&
     VECTOR_DIM,pnode,pgaus,elcon,elvel,gpsha,gpcar,&
     gpcon,gpvel,gpmas,hleng,gpdiv,gpdis,gpprd,gptur,gpden,gpdtr,&
     gpmascon,list_elements_p)

  !-----------------------------------------------------------------------
  !****f* Chemic/chm_elmpre_flamLet
  ! NAME
  !    chm_elmpre
  ! DESCRIPTION
  !    Compute quantities at Gauss points for flamelet combustion model
  ! USES
  ! USED BY
  !    mod_chm_element_operations
  !***
  !-----------------------------------------------------------------------
  use def_domain, only      :  ndime,mnode
  use def_chemic, only      :  nclas_chm,&
                               diffu_chm,&
                               ADR_chm,mass_gp,massConsumption_gp
  use def_chemic, only      :  ncomp_chm
  use def_chemic, only      :  mixedEq_eqs_chm
  use mod_ADR,    only      :  BDF
  use mod_chm_mixedEq, only : CHM_EQ_ZVAR
  use mod_chm_mixedEq, only : CHM_EQ_ZZ
  use mod_chm_mixedEq, only : CHM_EQ_YCVAR
  use mod_chm_mixedEq, only : CHM_EQ_YCYC

  implicit none
  integer(ip),  intent(in)  :: VECTOR_DIM,pnode,pgaus
  real(rp),     intent(in)  :: elcon(VECTOR_DIM,pnode,nclas_chm,ADR_chm(1) % ntime)
  real(rp),     intent(in)  :: elvel(VECTOR_DIM,ndime,pnode)
  real(rp),     intent(in)  :: gpsha(pnode,pgaus)
  real(rp),     intent(in)  :: gpcar(VECTOR_DIM,ndime,mnode,pgaus)
  real(rp),     intent(in)  :: hleng(VECTOR_DIM,3)
  real(rp),     intent(in)  :: gpden(VECTOR_DIM,pgaus)                  ! Density at gauss points
  real(rp),     intent(in)  :: gpdtr(VECTOR_DIM,pgaus)                  ! Dt*rho at gauss points
  real(rp),     intent(in)  :: gptur(VECTOR_DIM,pgaus)                  ! Turbulence viscosity at gauss points
  real(rp),     intent(out) :: gpcon(VECTOR_DIM,pgaus,nclas_chm,ncomp_chm)
  real(rp),     intent(out) :: gpvel(VECTOR_DIM,ndime,pgaus)
  real(rp),     intent(out) :: gpmas(VECTOR_DIM,pgaus,nclas_chm)        ! Mass source term
  real(rp),     intent(out) :: gpdiv(VECTOR_DIM,pgaus)                  ! Velocity divergence
  real(rp),     intent(out) :: gpdis(VECTOR_DIM,pgaus,nclas_chm)        ! Dissipation rate for combustion model
  real(rp),     intent(out) :: gpprd(VECTOR_DIM,pgaus,nclas_chm)        ! Production term of the variance of c and z in the
                                                                        !  combustion model
  real(rp),     intent(out) :: gpmascon(VECTOR_DIM,pgaus,nclas_chm)     ! Mass consumption term / Y_k
  integer(ip),  intent(in)  :: list_elements_p(VECTOR_DIM)              ! List of elements (always positive)


  real(rp)                  :: delta(VECTOR_DIM),delta2(VECTOR_DIM),rdelta2(VECTOR_DIM),seci4(VECTOR_DIM)
  real(rp)                  :: gpgve(VECTOR_DIM,ndime,ndime,pgaus)
  real(rp)                  :: rtau_turb(VECTOR_DIM)                     ! Inverse of turbulent time scale

  !
  ! Generic variance
  !
  real(rp)                  :: sgs_Chi(VECTOR_DIM),res_Chi(VECTOR_DIM)
  real(rp)                  :: varia(VECTOR_DIM)
  real(rp)                  :: grad_Mean(VECTOR_DIM,pgaus,ndime),mod_grad_Mean(VECTOR_DIM,pgaus)

  integer(ip)               :: igaus,iclas,inode,idime,jdime,itime
  integer(ip)               :: ielem,ielemone,ivect

  real(rp)                  :: auxprod(VECTOR_DIM)
  !
  ! Initialization
  !
  gpmas = 0.0_rp
  gpvel = 0.0_rp
  gpdis = 0.0_rp
  gpprd = 0.0_rp
  gpmascon = 0.0_rp


  ielemone = list_elements_p(1)


  !
  ! Concentration
  !
  do iclas = 1,nclas_chm
     do igaus = 1,pgaus
        gpcon(DEF_VECT,igaus,iclas,1) = 0.0_rp
        do inode = 1,pnode
           gpcon(DEF_VECT,igaus,iclas,1) = gpcon(DEF_VECT,igaus,iclas,1)&
                                + gpsha(inode,igaus) * elcon(DEF_VECT,inode,iclas,1)
        end do
     end do
  end do

  !
  ! Time integration
  !
  if( ADR_chm(1) % kfl_time_integration == 1 ) then
     do iclas = 1,nclas_chm
        do igaus = 1,pgaus
           gpcon(DEF_VECT,igaus,iclas,2) = 0.0_rp
           do inode = 1,pnode
              gpcon(DEF_VECT,igaus,iclas,2) = gpcon(DEF_VECT,igaus,iclas,2)&
                                   + gpsha(inode,igaus) * elcon(DEF_VECT,inode,iclas,2)
           end do
        end do
     end do

     if( ADR_chm(1) % kfl_time_scheme == BDF ) then
        do iclas = 1,nclas_chm
           do itime = 3,ADR_chm(iclas) % kfl_time_order + 1
              do igaus = 1,pgaus
                 gpcon(DEF_VECT,igaus,iclas,itime) = 0.0_rp
                 do inode = 1,pnode
                    gpcon(DEF_VECT,igaus,iclas,itime) = gpcon(DEF_VECT,igaus,iclas,itime)&
                                             + gpsha(inode,igaus) * elcon(DEF_VECT,inode,iclas,itime)
                 end do
              end do
           end do
        end do
     end if
  end if

  !
  ! Fluid velocity
  !
  do igaus = 1, pgaus
     do inode = 1,pnode
        do idime = 1,ndime
           gpvel(DEF_VECT,idime,igaus) = gpvel(DEF_VECT,idime,igaus) &
                +  gpsha(inode,igaus) * elvel(DEF_VECT,idime,inode)
        end do
     enddo
  enddo

  !
  ! Fluid velocity divergence
  !
  do igaus = 1,pgaus
     do inode = 1,pnode
        do idime = 1,ndime
           gpdiv(DEF_VECT,igaus) = gpdiv(DEF_VECT,igaus) + gpcar(DEF_VECT,idime,inode,igaus) * elvel(DEF_VECT,idime,inode)
        end do
     end do
  end do


  !
  ! Mass source and consumption terms
  !
  do ivect=1,VECTOR_DIM
     ielem = list_elements_p(ivect)
     if(ivect== 1 .or. (ivect > 1 .and. ielem /=ielemone)) then
        gpmas(ivect,:,:)    = mass_gp(ielem) % a(:,:,1)
        gpmascon(ivect,:,:) = massConsumption_gp(ielem) % a(:,:,1)
     end if
  end do

  !
  ! Special production and dissipation terms that depend on other variables
  !
  do iclas = 1,nclas_chm
     select case(mixedEq_eqs_chm(iclas) % kfl_eqtype)
     case ( CHM_EQ_ZVAR,  &
          & CHM_EQ_ZZ,    &
          & CHM_EQ_YCVAR, &
          & CHM_EQ_YCYC)
         !
         ! Variance of passive scalar
         !
         grad_Mean     = 0.0_rp
         mod_grad_Mean = 0.0_rp

         !
         ! Gradient of mean
         !
         do igaus = 1,pgaus
            do inode = 1,pnode
               do idime = 1,ndime
                  grad_Mean(DEF_VECT,igaus,idime) = grad_Mean(DEF_VECT,igaus,idime) +  &
                   &  gpcar(DEF_VECT,idime,inode,igaus) * elcon(DEF_VECT,inode,mixedEq_eqs_chm(iclas) % kfl_ieq_mean,1)
               end do
            end do
         end do

         !
         ! Length scales
         !
         delta(DEF_VECT)   = ( hleng(DEF_VECT,1) * hleng(DEF_VECT,2) * hleng(DEF_VECT,3) )** 0.3333333_rp
         delta2(DEF_VECT)  = delta(DEF_VECT)*delta(DEF_VECT)
         rdelta2(DEF_VECT) = 1.0_rp / delta2(DEF_VECT)

         !
         ! Loop over GPs
         !
         gauss_loop: do igaus = 1,pgaus
            !
            ! |grad Mean|^2
            !
            auxprod = 0.0_rp
            do idime=1,ndime
              auxprod(DEF_VECT) = auxprod(DEF_VECT) +  grad_Mean(DEF_VECT,igaus,idime)*grad_Mean(DEF_VECT,igaus,idime)
            end do
            mod_grad_Mean(DEF_VECT,igaus) = auxprod(DEF_VECT)

            !
            ! Compute strain rate square |S_ij|*|S_ij|
            !
            gpgve = 0.0_rp
            do inode = 1,pnode
               do idime = 1,ndime
                  do jdime = 1,ndime
                     gpgve(DEF_VECT,jdime,idime,igaus) = gpgve(DEF_VECT,jdime,idime,igaus) + elvel(DEF_VECT,idime,inode)&
                         * gpcar(DEF_VECT,jdime,inode,igaus)
                  end do
               end do
            end do

            seci4(DEF_VECT) = 0.0_rp
            do idime = 1,ndime                     ! |S_ij|*|S_ij|
               do jdime = 1,ndime
                  seci4(DEF_VECT) = seci4(DEF_VECT) + 0.5_rp * gpgve(DEF_VECT,idime,jdime,igaus)&
                      * (gpgve(DEF_VECT,idime,jdime,igaus) + gpgve(DEF_VECT,jdime,idime,igaus))
               end do
            end do

            !
            ! Compute inverse of turbulent time scale: 1 / tau = ( C_eps **2 * nu_t * |S|**2 / Delta**2 ) **(1/3)
            !
            rtau_turb(DEF_VECT) = ( 3.24_rp * gptur(DEF_VECT,igaus) * seci4(DEF_VECT) * rdelta2(DEF_VECT) ) ** (0.33333_rp)


            !
            ! Specific terms in case of non-zero source terms
            !
            if (mixedEq_eqs_chm(iclas) % kfl_eqtype == CHM_EQ_YCVAR ) then
               !
               ! Transport of reaction rate fluctuations: W_c  = 2 * ( {w_c * c} - {w_c}*{c} )
               !                                          W_Yc = 2 * ( {w_Yc * Yc} - {w_Yc}*{Yc} )
               !
               gpmas(DEF_VECT,igaus,iclas) = 2.0_rp * (gpmas(DEF_VECT,igaus,iclas) - gpmas(DEF_VECT,igaus,&
                   mixedEq_eqs_chm(iclas) % kfl_ieq_mean) * gpcon(DEF_VECT,igaus,mixedEq_eqs_chm(iclas) % kfl_ieq_mean,1))

            elseif (mixedEq_eqs_chm(iclas) % kfl_eqtype == CHM_EQ_YCYC ) then
               !
               ! Transport of reaction rate fluctuations: W_Y_c = 2 * {w_Y_c * Y_c}
               !
               gpmas(DEF_VECT,igaus,iclas) = 2.0_rp * gpmas(DEF_VECT,igaus,iclas)
            endif




            !
            ! Common terms for all scalar variances
            !
            if (mixedEq_eqs_chm(iclas) % kfl_eqtype == CHM_EQ_ZVAR .or. mixedEq_eqs_chm(iclas) % kfl_eqtype == CHM_EQ_YCVAR ) then
               !
               ! Production term: P_vari = 2 * D_t * |grad Mean|^2
               !
               gpprd(DEF_VECT,igaus,iclas) = 2.0_rp * gptur(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * mod_grad_Mean(DEF_VECT,igaus)&
                   / diffu_chm(1,1)

               !
               ! Dissipation term: D_vari = 2 * rho * 1 / tau * vari
               !
               gpdis(DEF_VECT,igaus,iclas) = 2.0_rp * rtau_turb(DEF_VECT) * gpden(DEF_VECT,igaus) * gpcon(DEF_VECT,igaus,iclas,1)

            elseif (mixedEq_eqs_chm(iclas) % kfl_eqtype == CHM_EQ_ZZ .or. mixedEq_eqs_chm(iclas) % kfl_eqtype == CHM_EQ_YCYC ) then
               !
               ! Production term: P_mean*Mean = 0
               !
               gpprd(DEF_VECT,igaus,iclas) = 0.0_rp
               !
               ! Dissipation term: D_MM = 2 D |grad Mean|^2 + X_sgs,   X_sgs = 2 * rho * 1 / tau * vari
               !
               varia(DEF_VECT) = gpcon(DEF_VECT,igaus,iclas,1) - gpcon(DEF_VECT,igaus,mixedEq_eqs_chm(iclas) % kfl_ieq_mean,1)&
                   * gpcon(DEF_VECT,igaus,mixedEq_eqs_chm(iclas) % kfl_ieq_mean,1)
               res_Chi(DEF_VECT) = 2.0_rp * gpdtr(DEF_VECT,igaus) / mixedEq_eqs_chm(mixedEq_eqs_chm(iclas) % kfl_ieq_mean) % Lewis&
                   * mod_grad_Mean(DEF_VECT,igaus)
               sgs_Chi(DEF_VECT) = 2.0_rp * rtau_turb(DEF_VECT) * gpden(DEF_VECT,igaus) * varia(DEF_VECT)

               gpdis(DEF_VECT,igaus,iclas) = res_Chi(DEF_VECT) + sgs_Chi(DEF_VECT)
            endif

         enddo gauss_loop

     end select
  enddo

end subroutine chm_elmpre_flamLet_fast


subroutine chm_elmpre_spray_fast(VECTOR_DIM,pnode,pgaus,elvel,gpcar,gpcon,gpden,gptur,hleng,gpsigma,list_elements_p)
!-----------------------------------------------------------------------
!****f* Chemic/chm_elmpre_spray
! NAME
!    chm_elmpre_spray
! DESCRIPTION
!    Compute gauss point terms for spray model
! USES
! USED BY
!    chm_elmcfi
!***
!-----------------------------------------------------------------------
use def_domain, only      :  ndime,mnode
use def_chemic, only      :  nclas_chm,surf_tension_chm,ncomp_chm
use def_chemic, only      :  sigma_gp_chm,sigma0_gp_chm,d32_gp_chm
use def_kermod, only      :  turmu_ker

implicit none
integer(ip), intent(in)   :: VECTOR_DIM
integer(ip),  intent(in)  :: pnode
integer(ip),  intent(in)  :: pgaus
real(rp),     intent(in)  :: elvel(VECTOR_DIM,ndime,pnode)
real(rp),     intent(in)  :: gpcar(VECTOR_DIM,ndime,mnode,pgaus)
real(rp),     intent(in)  :: gpcon(VECTOR_DIM,pgaus,ncomp_chm)
real(rp),     intent(in)  :: gpden(VECTOR_DIM,pgaus)
real(rp),     intent(in)  :: gptur(VECTOR_DIM,pgaus)
real(rp),     intent(in)  :: hleng(VECTOR_DIM,3)
real(rp),     intent(out) :: gpsigma(VECTOR_DIM,pgaus)
integer(ip),   intent(in)  :: list_elements_p(VECTOR_DIM)           ! List of elements (always positive)

real(rp)                  :: gpgve(VECTOR_DIM,ndime,ndime,pgaus)
real(rp)                  :: delta(VECTOR_DIM),delta2(VECTOR_DIM),rdelta(VECTOR_DIM),rdelta2(VECTOR_DIM)
real(rp)                  :: seci4(VECTOR_DIM),rtau_turb(VECTOR_DIM),k_sgs(VECTOR_DIM),sigma_crit(VECTOR_DIM)
real(rp)                  :: phi_clip(VECTOR_DIM)
integer(ip)               :: igaus,idime,jdime,inode
integer(ip)               :: iclas_phi
integer(ip)               :: iclas_sigma
integer(ip)               :: ielem, ivect, ielemone

real(rp)                  :: vec_sigma_gp_chm(VECTOR_SIZE,pgaus)
real(rp)                  :: vec_sigma0_gp_chm(VECTOR_SIZE,pgaus)
real(rp)                  :: vec_d32_gp_chm(VECTOR_SIZE,pgaus)

ielemone = list_elements_p(1)

vec_sigma_gp_chm  = 0.0_rp
vec_sigma0_gp_chm = 0.0_rp
vec_d32_gp_chm    = 0.0_rp

do ivect =1, VECTOR_DIM

   ielem = list_elements_p(ivect)
   if(ivect==1 .or. (ivect>1 .and. ielem /= ielemone)) then

   !
   ! Define liquid volume fraction and surface density variables
   !
   iclas_phi   = nclas_chm - 1
   iclas_sigma = nclas_chm

   !
   ! Initialization
   !


  ! sigma_gp_chm(ielem) % a  = 0.0_rp
  ! sigma0_gp_chm(ielem) % a = 0.0_rp
   !d32_gp_chm(ielem) % a    = 0.0_rp

     do igaus = 1,pgaus

        phi_clip(ivect) = min( 1.0_rp,max(0.0_rp,gpcon(ivect,igaus,iclas_phi)) )

        !
        ! Compute length scale from filter size
        !
        if (ndime == 2 ) then
            delta(ivect)   = ( hleng(ivect,1) * hleng(ivect,2) )**(0.5_rp)
        else
            delta(ivect)   = ( hleng(ivect,1) * hleng(ivect,2) * hleng(ivect,3) )**(0.3333333_rp)
        endif

         delta2(ivect)  = delta(ivect)*delta(ivect)
         rdelta(ivect)  = 1.0_rp / delta(ivect)
         rdelta2(ivect) = 1.0_rp / delta2(ivect)

         !
         ! Compute Sigma_min = alpha / Delta * sqrt (phi * (1 - phi))
         !    Chesnel et al., (2009) Int. J. Multifase Flow
         !
         vec_sigma0_gp_chm(ivect,igaus) = 2.4_rp * rdelta(ivect) * sqrt( min(0.25_rp,max(0.0_rp, (phi_clip(ivect)&
             * (1.0_rp - phi_clip(ivect))) )) )
       !  sigma0_gp_chm(ielem) % a(igaus,1,1) = 2.4_rp * rdelta(ivect) * sqrt( min(0.25_rp,max(0.0_rp, (phi_clip(ivect)
       !      * (1.0_rp - phi_clip(ivect))) )) )

         !
         ! Compute Sigma = Sigma_min + Sigma'
         !
         vec_sigma_gp_chm(ivect,igaus) = max(0.0_rp,vec_sigma0_gp_chm(ivect,igaus) + gpcon(ivect,igaus,iclas_sigma))
 !        sigma_gp_chm(ielem) % a(igaus,1,1) = max(0.0_rp,sigma0_gp_chm(ielem) % a(igaus,1,1) + gpcon(ivect,igaus,iclas_sigma))

         !
         ! Compute Sauter mean diameter d32 ( scaled by phi_L*(1-phi_L)/0.25 )
         !
         if ( vec_sigma_gp_chm(ivect,igaus) > 0.0_rp ) then
             vec_d32_gp_chm(ivect,igaus) = 4.0_rp * phi_clip(ivect) * (1.0_rp - phi_clip(ivect)) * max(0.0_rp,6.0_rp&
                 * phi_clip(ivect) / vec_sigma_gp_chm(ivect,igaus))
!            d32_gp_chm(ielem) % a(igaus,1,1) = 4.0_rp * phi_clip(ivect) * (1.0_rp - phi_clip(ivect)) * max(0.0_rp,6.0_rp&
!                * phi_clip(ivect) / sigma_gp_chm(ielem) % a(igaus,1,1))
         else
             vec_d32_gp_chm(ivect,igaus) = 0.0_rp
 !            d32_gp_chm(ielem) % a(igaus,1,1) = 0.0_rp
         end if

         !
         ! LES
         !
         if(turmu_ker % kfl_exist /= 0_ip) then
           !
           ! Compute strain rate square |S_ij|*|S_ij|
           !
           gpgve = 0.0_rp
           do inode = 1,pnode
               do idime = 1,ndime
                   do jdime = 1,ndime
                       gpgve(ivect,jdime,idime,igaus) = gpgve(ivect,jdime,idime,igaus) + elvel(ivect,idime,inode)&
                           * gpcar(ivect,jdime,inode,igaus)
                    end do
                end do
            end do

            seci4 = 0.0_rp
            do idime = 1,ndime                     ! |S_ij|*|S_ij|
                do jdime = 1,ndime
                    seci4(ivect) = seci4(ivect) + 0.5_rp * gpgve(ivect,idime,jdime,igaus) * (gpgve(ivect,idime,jdime,igaus)&
                        + gpgve(ivect,jdime,idime,igaus))
                 end do
             end do

            !
            ! Compute inverse of turbulent time scale (tau_b ~= tau_t): 1 / tau = ( C_eps **2 * nu_t * |S|**2 / Delta**2 ) **(1/3)
            !
            rtau_turb(ivect) = ( 3.24_rp * gptur(ivect,igaus) * seci4(ivect) * rdelta2(ivect) ) ** (0.33333_rp)

            !
            ! Compute turbulent kinetic energy, k = eps_sgs * tau_t
            !
            if ( rtau_turb(ivect) == 0.0_rp ) then
                gpsigma(ivect,igaus) = 0.0_rp

            else
                k_sgs(ivect)          = gptur(ivect,igaus) * abs( seci4(ivect) ) / rtau_turb(ivect)

                !
                ! Compute production/destruction of surface area by mean and shear flow, turbulence and interactions
                !
                sigma_crit(ivect) = max(0.0_rp,vec_sigma0_gp_chm(ivect,igaus) + phi_clip(ivect) * (1.0_rp - phi_clip(ivect))&
                    * gpden(ivect,igaus) * k_sgs(ivect) / surf_tension_chm)
!               sigma_crit(ivect) = max(0.0_rp,sigma0_gp_chm(ielem) % a(igaus,1,1) + phi_clip(ivect) * (1.0_rp - phi_clip(ivect))&
!                   * gpden(ivect,igaus) * k_sgs(ivect) / surf_tension_chm)

                if ( sigma_crit(ivect) == 0.0_rp ) then
                    gpsigma(ivect,igaus) = 0.0_rp
                else
                    gpsigma(ivect,igaus) = vec_sigma_gp_chm(ivect,igaus) * rtau_turb(ivect) * ( 1.0_rp&
                        - vec_sigma_gp_chm(ivect,igaus) / sigma_crit(ivect) )
!                   gpsigma(ivect,igaus) = sigma_gp_chm(ielem) % a(igaus,1,1) * rtau_turb(ivect) * ( 1.0_rp&
!                       - sigma_gp_chm(ielem) % a(igaus,1,1) / sigma_crit(ivect) )
                 endif

             end if
         !
         ! Laminar
         !
         else
             gpsigma(ivect,igaus) = vec_sigma0_gp_chm(ivect,igaus)
 !            gpsigma(ivect,igaus) = sigma0_gp_chm(ielem) % a(igaus,1,1)
         end if

    end do
  end if
end do


do ivect =1, VECTOR_DIM

   ielem = list_elements_p(ivect)
   if(ivect==1 .or. (ivect>1 .and. ielem /= ielemone)) then
     do igaus = 1,pgaus

        sigma0_gp_chm(ielem) % a(igaus,1,1) = vec_sigma0_gp_chm(ivect,igaus)
        sigma_gp_chm(ielem) % a(igaus,1,1)  = vec_sigma_gp_chm(ivect,igaus)
        d32_gp_chm(ielem) % a(igaus,1,1)    = vec_d32_gp_chm(ivect,igaus)

     end do

   end if
end do

end subroutine chm_elmpre_spray_fast

subroutine chm_rhodt_fast(VECTOR_DIM,pnode,pgaus,porde,lnods_loc,gpsha,gpvol,gpden,dt_rho_chm_loc,list_elements_p)
  !------------------------------------------------------------------------
  ! NAME
  !    chm_rhodt
  ! DESCRIPTION
  !    Projection of rho/dt for explicit time step
  ! USES
  ! USED BY
  !    chm_elmcfi
  !***
  !------------------------------------------------------------------------

  use def_master,      only : dtinv

  implicit none
  integer(ip), intent(in)  :: VECTOR_DIM
  integer(ip), intent(in)  :: pnode
  integer(ip), intent(in)  :: pgaus
  integer(ip), intent(in)  :: porde
  integer(ip), intent(in)  :: lnods_loc(VECTOR_DIM,pnode)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: gpvol(VECTOR_DIM,pgaus)
  real(rp),    intent(in)  :: gpden(VECTOR_DIM,pgaus)
  real(rp),    intent(out) :: dt_rho_chm_loc(*)
  integer(ip),   intent(in)  :: list_elements_p(VECTOR_DIM)           ! List of elements (always positive)
  integer(ip) :: inode,jnode,ipoin,igaus
  real(rp)     :: eldtrho(VECTOR_DIM,pnode), fact(VECTOR_DIM)
  real(rp)                :: elmat(VECTOR_DIM,pnode,pnode)
  real(rp)                :: trace(VECTOR_DIM), elmass(VECTOR_DIM)

  integer(ip)             :: ielem,ivect,ielemone

  eldtrho = 0.0_rp
  elmat   = 0.0_rp
  trace   = 0.0_rp
  elmass  = 0.0_rp

  ielemone  = list_elements_p(1)



       if( porde == 1 ) then
          !
          ! Element assembly
          !
          do igaus = 1,pgaus
             fact(DEF_VECT) = gpvol(DEF_VECT,igaus) / ( gpden(DEF_VECT,igaus)  * dtinv )
             do inode = 1,pnode
                eldtrho(DEF_VECT,inode) = eldtrho(DEF_VECT,inode) + gpsha(inode,igaus) * fact(DEF_VECT)
             end do
          end do

       else
          !
          ! Element assembly
          !
     !@     eldtrho = 0.0_rp



     !@     do inode=1,pnode
     !@        do jnode=1,pnode
     !@           elmat(ivect,inode,jnode)=0.0_rp
     !@        end do
     !@     end do

          do igaus=1,pgaus
             do inode=1,pnode
                fact(DEF_VECT)=gpvol(DEF_VECT,igaus)*gpsha(inode,igaus)/(gpden(DEF_VECT,igaus)*dtinv)
                do jnode=1,pnode
                   elmat(DEF_VECT,inode,jnode)=elmat(DEF_VECT,inode,jnode) +fact(DEF_VECT)*gpsha(jnode,igaus)
                end do
             end do
          end do

     !@     trace  = 0.0_rp
     !@     elmass = 0.0_rp
          do inode = 1,pnode
             trace(DEF_VECT) = trace(DEF_VECT) + elmat(DEF_VECT,inode,inode)
             do jnode = 1,pnode
                elmass(DEF_VECT) = elmass(DEF_VECT) + elmat(DEF_VECT,inode,jnode)
             end do
          end do

  end if

!  end if

!end do


 !
 ! Nodal projection
 !

  do ivect=1, VECTOR_DIM
      ielem = list_elements_p(ivect)
      if( ivect == 1 .or. (ivect>1 .and. ielem /= ielemone)) then
          if( porde == 1 ) then
              do inode = 1,pnode
                  ipoin = lnods_loc(ivect,inode)
                  dt_rho_chm_loc(ipoin) = dt_rho_chm_loc(ipoin) + eldtrho(ivect,inode)
              end do

          else
             do inode = 1,pnode
                ipoin = lnods_loc(ivect,inode)
                dt_rho_chm_loc(ipoin) = dt_rho_chm_loc(ipoin) + elmat(ivect,inode,inode)*(elmass(ivect)/trace(ivect))
            end do
          end if

     end if
end do
end subroutine chm_rhodt_fast

   subroutine chm_rhodt_spray_fast(VECTOR_DIM,pnode,pgaus,porde,lnods_loc,gpsha,gpvol,gpden,&
       dt_rho_chm_loc,dt_chm_loc,list_elements_p)
     !------------------------------------------------------------------------
     ! NAME
     !    chm_rhodt_spray
     ! DESCRIPTION
     !    Projection of rho/dt for explicit time step for the gas phase
     !    and 1/dt for the liquid phase
     ! USES
     ! USED BY
     !    chm_element_operations
     !***
     !------------------------------------------------------------------------

     use def_master,      only : dtinv

     implicit none

     integer(ip), intent(in)  :: VECTOR_DIM
     integer(ip), intent(in)  :: pnode
     integer(ip), intent(in)  :: pgaus
     integer(ip), intent(in)  :: porde
     integer(ip), intent(in)  :: lnods_loc(VECTOR_DIM,pnode)
     real(rp),    intent(in)  :: gpsha(pnode,pgaus)
     real(rp),    intent(in)  :: gpvol(VECTOR_DIM,pgaus)
     real(rp),    intent(in)  :: gpden(VECTOR_DIM,pgaus)
     real(rp),    intent(out) :: dt_rho_chm_loc(*)
     real(rp),    intent(out) :: dt_chm_loc(*)
     integer(ip),   intent(in)  :: list_elements_p(VECTOR_DIM)           ! List of elements (always positive)

     integer(ip)             :: inode,jnode,ipoin,igaus
     real(rp)                :: eldtrho(VECTOR_DIM,pnode),eldt(VECTOR_DIM,pnode)
     real(rp)                :: fact1(VECTOR_DIM),fact2(VECTOR_DIM)
     real(rp)                :: elmat1(VECTOR_DIM,pnode,pnode),elmat2(VECTOR_DIM,pnode,pnode)
     real(rp)                :: trace1(VECTOR_DIM),elmass1(VECTOR_DIM)
     real(rp)                :: trace2(VECTOR_DIM),elmass2(VECTOR_DIM)
     integer(ip)             :: ielem,ivect,ielemone

     eldtrho = 0.0_rp
     elmat1 = 0.0_rp
     elmat2 = 0.0_rp

     trace1   = 0.0_rp
     elmass1  = 0.0_rp
     trace2   = 0.0_rp
     elmass2  = 0.0_rp

     ielemone  = list_elements_p(1)

     do ivect=1, VECTOR_DIM
         ielem = list_elements_p(ivect)

         if( ivect == 1 .or. (ivect>1 .and. ielem /= ielemone)) then



             if( porde == 1 ) then
                !
                ! Element assembly
                !
                !@eldtrho = 0.0_rp
                !@eldt    = 0.0_rp
                do igaus = 1,pgaus
                   fact1(ivect) = gpvol(ivect,igaus) / ( gpden(ivect,igaus)  * dtinv )
                   fact2(ivect) = gpvol(ivect,igaus) / dtinv
                   do inode = 1,pnode
                      eldtrho(ivect,inode) = eldtrho(ivect,inode) + gpsha(inode,igaus) * fact1(ivect)
                      eldt(ivect,inode)    = eldt(ivect,inode)    + gpsha(inode,igaus) * fact2(ivect)
                   end do
                end do
             else
                !
                ! Element assembly
                !
 !@               eldtrho = 0.0_rp
 !@               eldt    = 0.0_rp
 !@
 !@               elmat1 = 0.0_rp
 !@               elmat2 = 0.0_rp

                do igaus=1,pgaus
                   do inode=1,pnode
                      fact1(ivect) = gpvol(ivect,igaus)*gpsha(inode,igaus)/(gpden(ivect,igaus)*dtinv)
                      fact2(ivect) = gpvol(ivect,igaus)*gpsha(inode,igaus)/dtinv
                      do jnode=1,pnode
                         elmat1(ivect,inode,jnode)=elmat1(ivect,inode,jnode) + fact1(ivect)*gpsha(jnode,igaus)
                         elmat2(ivect,inode,jnode)=elmat2(ivect,inode,jnode) + fact2(ivect)*gpsha(jnode,igaus)
                      end do
                   end do
                end do

!@                trace1  = 0.0_rp
!@                trace2  = 0.0_rp
!@                elmass1 = 0.0_rp
!@                elmass2 = 0.0_rp

                do inode = 1,pnode
                   trace1(ivect) = trace1(ivect) + elmat1(ivect,inode,inode)
                   trace2(ivect) = trace2(ivect) + elmat2(ivect,inode,inode)
                   do jnode = 1,pnode
                      elmass1(ivect) = elmass1(ivect) + elmat1(ivect,inode,jnode)
                      elmass2(ivect) = elmass2(ivect) + elmat2(ivect,inode,jnode)
                   end do
                end do

             end if

         end if
     end do

    !
    ! Nodal projection
    !

     do ivect=1, VECTOR_DIM
         ielem = list_elements_p(ivect)

         if( ivect == 1 .or. (ivect>1 .and. ielem /= ielemone)) then



             if( porde == 1 ) then
                do inode = 1,pnode
                   ipoin = lnods_loc(ivect,inode)
                   dt_rho_chm_loc(ipoin) = dt_rho_chm_loc(ipoin) + eldtrho(ivect,inode)
                   dt_chm_loc(ipoin)     = dt_chm_loc(ipoin)     + eldt(ivect,inode)
                end do
             else
                do inode = 1,pnode
                   ipoin = lnods_loc(ivect,inode)
                   dt_rho_chm_loc(ipoin) = dt_rho_chm_loc(ipoin) + elmat1(ivect,inode,inode)*(elmass1(ivect)/trace1(ivect))
                   dt_chm_loc(ipoin)     = dt_chm_loc(ipoin)     + elmat2(ivect,inode,inode)*(elmass2(ivect)/trace2(ivect))
                end do
             end if

         end if
     end do
   end subroutine chm_rhodt_spray_fast

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux and Matias Avila
  !> @brief   Add SGS or Bubble to unknown
  !> @details u <= u + u'
  !
  !-----------------------------------------------------------------------

  subroutine chm_ADR_add_sgs_or_bubble_fast(VECTOR_DIM,iclas,pgaus,gpsha_bub,ADR,gpcon,list_elements_p)
    use def_chemic,          only : nclas_chm,ncomp_chm
    use mod_ADR,             only : ADR_typ,BUBBLE

    integer(ip),   intent(in)    :: VECTOR_DIM,iclas            !< VECTOR_DIM
    integer(ip),   intent(in)    :: pgaus            !< Number of Gauss points
    real(rp),      intent(in)    :: gpsha_bub(pgaus) !< Bubble shape function
    type(ADR_typ), intent(in)    :: ADR              !< ADR type
    real(rp),      intent(inout) :: gpcon(VECTOR_DIM,pgaus,nclas_chm,ncomp_chm)   !< Unknown
    integer(ip),   intent(in)    :: list_elements_p(VECTOR_DIM)           ! List of elements (always positive)
    integer(ip)                  :: ivect,ielem,ielemone

    ielemone  = list_elements_p(1)

     do ivect=1, VECTOR_DIM
         ielem = list_elements_p(ivect)

         if( ivect == 1 .or. (ivect>1 .and. ielem /= ielemone)) then

            if( ADR % kfl_stabilization == BUBBLE ) then
               !
               !      n
               ! u = sum Ni*ui + Ne*ue
               !     i=1
               !
               gpcon(ivect,1:pgaus,iclas,1) = gpcon(ivect,1:pgaus,iclas,1) + ADR % bubble(ielem,1) * gpsha_bub(1:pgaus)

            else if( ADR % kfl_nonlinear_sgs /= 0 ) then
               !
               !      n
               ! u = sum Ni*ui + u'
               !     i=1
               !
               gpcon(ivect,1:pgaus,iclas,1) = gpcon(ivect,1:pgaus,iclas,1) + ADR % sgs(ielem) % a(1,1:pgaus,1)

            end if

        end if
    end do
  end subroutine chm_ADR_add_sgs_or_bubble_fast

subroutine chm_elmprc_flamLet_fast(&
  VECTOR_DIM,iclas,pgaus,gpden,gpmas,gpmascon,gpdtr,gptur,gpdis,gpprd, &
  gpdif,gprhs,gprea)

  !-----------------------------------------------------------------------
  !****f* chemic/chm_elmprc
  ! NAME
  !    chm_elmprc
  ! DESCRIPTION
  !    Compute terms for each species ADR equation
  ! USES
  ! USED BY
  !    chm_element_operations
  !***
  !-----------------------------------------------------------------------
  use def_chemic, only      :  nclas_chm, mixedEq_eqs_chm
  use def_chemic, only      :  diffu_chm
  use def_kermod, only      :  turmu_ker
  use mod_ADR,    only      :  mreac_adr

  implicit none
  integer(ip),  intent(in)  :: VECTOR_DIM
  integer(ip),  intent(in)  :: iclas
  integer(ip),  intent(in)  :: pgaus
  real(rp),     intent(in)  :: gpden(VECTOR_DIM,pgaus)
  real(rp),     intent(in)  :: gpmas(VECTOR_DIM,pgaus,nclas_chm)                ! Mass source term (reaction rate)
  real(rp),     intent(in)  :: gpmascon(VECTOR_DIM,pgaus,nclas_chm)             ! Mass consumption rate / Y_k
  real(rp),     intent(in)  :: gpdtr(VECTOR_DIM,pgaus)                          ! Dt*rho
  real(rp),     intent(in)  :: gptur(VECTOR_DIM,pgaus)                          ! turbulent viscosity
  real(rp),     intent(in)  :: gpdis(VECTOR_DIM,pgaus,nclas_chm)                ! Dissipation rate term in the Flamelet model
  real(rp),     intent(in)  :: gpprd(VECTOR_DIM,pgaus,nclas_chm)                ! Production term in the Flamelet model

  real(rp),     intent(out) :: gpdif(VECTOR_DIM,pgaus,nclas_chm)
  real(rp),     intent(out) :: gprhs(VECTOR_DIM,pgaus)
  real(rp),     intent(out) :: gprea(VECTOR_DIM,pgaus,mreac_adr)

  integer(ip)               :: igaus

  do igaus = 1,pgaus
     !
     ! RHS terms: Source + Production term + Dissipation term
     !
     gprhs(DEF_VECT,igaus) = gprhs(DEF_VECT,igaus) + gpmas(DEF_VECT,igaus,iclas) + gpprd(DEF_VECT,igaus,iclas) - gpdis(DEF_VECT,igaus,iclas)

     !
     ! Reac term: put mass consumption rate here, since it is a linear function of the unknown
     !
     !gprhs(DEF_VECT,igaus) = gprhs(DEF_VECT,igaus) + gpmascon(DEF_VECT,igaus,iclas) * gpcon(DEF_VECT,igaus,iclas)
     gprea(DEF_VECT,igaus,1) = gprea(DEF_VECT,igaus,1) - gpmascon(DEF_VECT,igaus,iclas)


     if (mixedEq_eqs_chm(iclas) % kfl_fix_diffusion == 1) then
        !
        ! Diffusion coefficient for soot (rho*Dsoot), or other predifined constant
        !
        gpdif(DEF_VECT,igaus,iclas) = gpden(DEF_VECT,igaus) * mixedEq_eqs_chm(iclas) % diffusivity
     else
        !
        ! Diffusion coefficient for conventional scalar transport
        !
        gpdif(DEF_VECT,igaus,iclas) = gpdtr(DEF_VECT,igaus) / mixedEq_eqs_chm(iclas) % Lewis
     endif


     !
     ! Adding turbulent part
     !
     if( turmu_ker % kfl_exist /= 0_ip ) then
        gpdif(DEF_VECT,igaus,iclas) = gpdif(DEF_VECT,igaus,iclas) + gptur(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) / diffu_chm(1,1)
     end if

 enddo

end subroutine chm_elmprc_flamLet_fast


   subroutine chm_elmprc_spray_fast(&
         VECTOR_DIM,iclas,pgaus,gptur,gpsigma,gpden,gpdif,gprhs)

     !-----------------------------------------------------------------------
     !****f* Chemic/chm_elmprc_spray
     ! NAME
     !    chm_elmprc_spray
     ! DESCRIPTION
     !    Compute spray terms for assembly in ADR equation
     ! USES
     ! USED BY
     !    chm_element_operations
     !***
     !-----------------------------------------------------------------------
     use def_chemic, only      :  kfl_spray_chm
     use def_chemic, only      :  diffu_chm
     use def_chemic, only      :  nclas_chm
     use def_kermod, only      :  turmu_ker

     implicit none

     integer(ip),  intent(in)  :: VECTOR_DIM
     integer(ip),  intent(in)  :: iclas
     integer(ip),  intent(in)  :: pgaus
     real(rp),     intent(in)  :: gptur(VECTOR_DIM,pgaus)
     real(rp),     intent(in)  :: gpsigma(VECTOR_DIM,pgaus)
     real(rp),     intent(out) :: gpden(VECTOR_DIM,pgaus)
     real(rp),     intent(out) :: gpdif(VECTOR_DIM,pgaus,nclas_chm)
     real(rp),     intent(out) :: gprhs(VECTOR_DIM,pgaus)

     !
     ! Solve for non-density weighted transport
     !
     if (iclas == nclas_chm ) then
        gpden(DEF_VECT,1:pgaus) = 1.0_rp
     end if

     !
     ! LES calculation
     !
     if (turmu_ker % kfl_exist /= 0_ip .and. kfl_spray_chm == 1) then
        gpdif(DEF_VECT,1:pgaus,iclas) = gptur(DEF_VECT,1:pgaus) / diffu_chm(1,1)

     !
     ! Level set calculation or laminar
     !
     else
        gpdif(DEF_VECT,1:pgaus,iclas) = 0.0_rp
     end if

     !
     !  Add production/destruction of surface area by interactions Sigma_int to the
     !  liquid-gas interface density
     !
     if (iclas == nclas_chm ) &
         gprhs(DEF_VECT,1:pgaus) = gprhs(DEF_VECT,1:pgaus) + gpsigma(DEF_VECT,1:pgaus)


   end subroutine chm_elmprc_spray_fast


   subroutine chm_element_assembly(&
            VECTOR_DIM,pnode,pgaus,list_elements_p,&
            gpsha,gpcar,gphes,gpvol,chale,ADR_chm, &
            gpden,gpvel,gpdif,gpgrd,gprea,gprhs,gpcon,elmat,elrhs)
      ! -------------------------------------------------------------------------------------------------------------------------- !
      ! DESCRIPTION                                                                                                                !   
      ! Element assembly operation in Chemic, it takes care of computing the linear problem matrix 'elmat' and the right hand side !
      ! vector 'elrhs'. For simplicity, we only use the following stabilization methods:                                           !
      ! 1) Galerkin or OFF (kfl_stabilization == -2)                                                                               !
      ! 2) ASGS (kfl_stabilization == 0)                                                                                           !
      !                                                                                                                            !   
      ! USED BY:                                                                                                                   !
      ! chm_element_operations_fast                                                                                                !
      ! -------------------------------------------------------------------------------------------------------------------------- !

      use def_domain, only     : mnode,ndime,ntens
      use mod_ADR,    only     : ADR_typ
      use def_chemic, only     : ncomp_chm

      ! Element dimensions
      integer(ip), intent(in)            :: VECTOR_DIM                               !< Current element
      integer(ip), intent(in)            :: pnode                                    !< # nodes
      integer(ip), intent(in)            :: pgaus                                    !< # Gauss points
      integer(ip), intent(in)            :: list_elements_p(VECTOR_DIM)              !< List of elements

      ! Element characteristics at Gauss point
      real(rp),    intent(in)            :: gpsha(pnode,pgaus)                       !< Shape function Nk
      real(rp),    intent(in)            :: gpcar(VECTOR_DIM,ndime,mnode,pgaus)      !< Shape function Cartesian derivatives dNk/dxi
      real(rp),    intent(in)            :: gphes(VECTOR_DIM,ntens,mnode,pgaus)      !< Hessian dNk/dxidxj
      real(rp),    intent(in)            :: gpvol(VECTOR_DIM,pgaus)                  !< Element Jacobian
      real(rp),    intent(in)            :: chale(VECTOR_DIM,2)                      !< Element characteristic length

      ! Numerical strategy
      type(ADR_typ), intent(inout)       :: ADR_chm                                  !< ADR type

      ! Equation coefficients
      real(rp),    intent(in)            :: gpden(VECTOR_DIM,pgaus)                  !< Density
      real(rp),    intent(in)            :: gpvel(VECTOR_DIM,ndime,pgaus)            !< Advection vector
      real(rp),    intent(in)            :: gpdif(VECTOR_DIM,pgaus)                  !< Diffusion
      real(rp),    intent(in)            :: gpgrd(VECTOR_DIM,ndime,pgaus)            !< Diffusion gradient
      real(rp),    intent(in)            :: gprea(VECTOR_DIM,pgaus,4)                !< Reaction
      real(rp),    intent(in)            :: gprhs(VECTOR_DIM,pgaus)                  !< RHS
      real(rp),    intent(in)            :: gpcon(VECTOR_DIM,pgaus,ncomp_chm)        !< Unknown at Gauss point

      ! Output
      real(rp),    intent(out)           :: elmat(VECTOR_DIM,pnode,pnode)            !< Element matrix
      real(rp),    intent(out)           :: elrhs(VECTOR_DIM,pnode)                  !< Element RHS

      ! To assembly of the matrix
      integer(ip) :: ielem, igaus, inode, jnode, ivect, ielemone, itime, idime       !< Indices and dimensions
      real(rp)    :: react(VECTOR_DIM,pgaus), sreac(VECTOR_DIM,pgaus), sgs(VECTOR_DIM,pgaus)
      real(rp)    :: rhsit(VECTOR_DIM, pgaus), gptau(VECTOR_DIM,pgaus), gptau_time(VECTOR_DIM,pgaus)
      real(rp)    :: gpnve(VECTOR_DIM), ka(VECTOR_DIM), aa(VECTOR_DIM), sa(VECTOR_DIM), tau(VECTOR_DIM)
      real(rp)    :: freq1(VECTOR_DIM), freq2(VECTOR_DIM), freq3(VECTOR_DIM)         !< Codina
      real(rp)    :: dtinv_elem
      real(rp)    :: resi1(VECTOR_DIM,pnode), resi2(VECTOR_DIM,pnode), gpad1(VECTOR_DIM,pnode)
      real(rp)    :: gpadv(VECTOR_DIM,pnode), gppe1(VECTOR_DIM,pnode), gppe2(VECTOR_DIM,pnode)
      real(rp)    :: grvgr(VECTOR_DIM), gplap(VECTOR_DIM), dd(VECTOR_DIM)
      real(rp)    :: xmuit(VECTOR_DIM)
      real(rp)    :: gppr1(VECTOR_DIM,pgaus)
      real(rp)    :: galer(VECTOR_DIM,pnode)

      ! To perform auxiliar math operations
      real(rp)    :: fact1(VECTOR_DIM), fact2(VECTOR_DIM), fact3(VECTOR_DIM), fact4(VECTOR_DIM)
      real(rp)    :: fact5(VECTOR_DIM), fact6(VECTOR_DIM), fact7(VECTOR_DIM)

      external    :: runend

      ielemone = list_elements_p(1) ! To avoid repeating terms in pseudo vectorized version

      ! Initialization
      sgs(DEF_VECT,1:pgaus) = 0.0_rp
      elrhs(DEF_VECT,1:pnode) = 0.0_rp
      elmat(DEF_VECT,1:pnode,1:pnode) = 0.0_rp
      gppr1(DEF_VECT,1:pgaus) = 0.0_rp
      dtinv_elem = ADR_chm % dtinv

      rhsit(DEF_VECT,1:pgaus) = gprhs(DEF_VECT,1:pgaus)

      do igaus = 1,pgaus
         fact1(DEF_VECT) = gpden(DEF_VECT,igaus) * dtinv_elem
         do itime = 2,ADR_chm % ntime
            rhsit(DEF_VECT,igaus) = rhsit(DEF_VECT,igaus) - fact1(DEF_VECT) * (ADR_chm % time_parameters(itime))*gpcon(DEF_VECT,igaus,itime)
         end do
      end do

      react(DEF_VECT,1:pgaus) = gprea(DEF_VECT,1:pgaus,1) + 2.0_rp * gprea(DEF_VECT,1:pgaus,2) * gpcon(DEF_VECT,1:pgaus,1)
      sreac(DEF_VECT,1:pgaus) = gprea(DEF_VECT,1:pgaus,1) + gprea(DEF_VECT,1:pgaus,2) * gpcon(DEF_VECT,1:pgaus,1)
      rhsit(DEF_VECT,1:pgaus) = rhsit(DEF_VECT,1:pgaus)   + gprea(DEF_VECT,1:pgaus,2) * (gpcon(DEF_VECT,1:pgaus,1)) ** 2


      if( ADR_chm%kfl_stabilization == -2 ) then  ! GALERKIN
         !-------------------------------------------------------------------
         !                                                                   
         ! Galerkin                                                              
         !
         ! k (grad(u),grad(v)) + (1-b)*(rho*a.grad(u),v) -b*(u,rho*a.grad(v))
         !
         !-------------------------------------------------------------------


         do igaus = 1,pgaus
            !
            ! Coefficients
            !
            fact2(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpdif(DEF_VECT,igaus)               ! k
            fact3(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus)               ! rho
            fact4(DEF_VECT) = gpvol(DEF_VECT,igaus) * rhsit(DEF_VECT,igaus)               ! f
            fact5(DEF_VECT) = gpvol(DEF_VECT,igaus) * react(DEF_VECT,igaus)               ! r
            fact6(DEF_VECT) = fact3(DEF_VECT) * dtinv_elem * ADR_chm % time_parameters(1) ! rho/dt
            fact7(DEF_VECT) = ADR_chm % bemol * fact3(DEF_VECT)
            
            !
            ! Advections = a.grad(Ni)
            !
            gpadv(DEF_VECT,1:pnode) = 0.0_rp
            do idime = 1,ndime
               do inode=1, pnode
                  gpadv(DEF_VECT,inode) = gpadv(DEF_VECT,inode) + gpvel(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
               end do
            end do

            !
            ! Galerkin operator = Nj/dt + r*Nj + rho*a.grad(Nj)
            !
            do inode=1, pnode
               galer(DEF_VECT,inode) = ( fact6(DEF_VECT) + fact5(DEF_VECT) ) * gpsha(inode,igaus) + (fact3(DEF_VECT) &
                     &* ( 1.0_rp - ADR_chm%bemol ) * gpadv(DEF_VECT,inode))
            end do

            !
            ! Galerkin operator + diffusion + convection = Galerkin operator + k grad(Ni).grad(Nj) - b*( u, rho*a.grad(v) )
            !
            do inode = 1,pnode
               do jnode = 1,pnode
                  elmat(DEF_VECT,inode,jnode) = elmat(DEF_VECT,inode,jnode) + gpsha(inode,igaus) * galer(DEF_VECT,jnode)
                  elmat(DEF_VECT,inode,jnode) = elmat(DEF_VECT,inode,jnode) - fact7(DEF_VECT) * gpadv(DEF_VECT,inode)&
                        &* gpsha(jnode,igaus)
                  do idime = 1,ndime
                     elmat(DEF_VECT,inode,jnode) = elmat(DEF_VECT,inode,jnode) + fact2(DEF_VECT)&
                        &* gpcar(DEF_VECT,idime,inode,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
                  end do
               end do
               elrhs(DEF_VECT,inode) = elrhs(DEF_VECT,inode) + gpsha(inode,igaus) * fact4(DEF_VECT)
            end do
         end do  ! Galerkin stabilization

      else if( ADR_chm%kfl_stabilization == 0) then ! ASGS stabilization

         ! Tau calculations
         elements44: do ivect = 1,VECTOR_DIM
            ielem=list_elements_p(ivect)
            if (ivect==1 .or. (ivect > 1 .and. ielem /= ielemone ))  then
               ! Element characteristic length
               if( chale(ivect,1) == 0.0_rp .and. chale(ivect,2) == 0.0_rp ) then
                  gptau(ivect,:) = 0.0_rp
               end if
            end if
         end do elements44

         do igaus = 1,pgaus
            gpnve(DEF_VECT) = 0.0_rp
   
            do idime=1,ndime
               gpnve(DEF_VECT) = gpnve(DEF_VECT) + gpvel(DEF_VECT,idime,igaus) * gpvel(DEF_VECT,idime,igaus)
            end do
   
            gpnve(DEF_VECT) = gpden(DEF_VECT,igaus) * sqrt( gpnve(DEF_VECT) + epsilon(1.0_rp) )
   
            ka(DEF_VECT) = ADR_chm % tau_parameters(1) * gpdif(DEF_VECT,igaus)  ! Diffusion
            aa(DEF_VECT) = ADR_chm % tau_parameters(2) * gpnve(DEF_VECT)        ! Advection
            sa(DEF_VECT) = ADR_chm % tau_parameters(3) * sreac(DEF_VECT,igaus)  ! Reaction
   
            !
            ! Codina
            !
            freq1(DEF_VECT) = 4.0_rp*ka(DEF_VECT) / (chale(DEF_VECT,2) * chale(DEF_VECT,2))
            freq2(DEF_VECT) = 2.0_rp*aa(DEF_VECT) / chale(DEF_VECT,1)
            freq3(DEF_VECT) = abs(sa(DEF_VECT))
            tau(DEF_VECT)   = freq1(DEF_VECT) + freq2(DEF_VECT) + freq3(DEF_VECT)
   
            elements45: do ivect = 1,VECTOR_DIM
               ielem=list_elements_p(ivect)
               if (ivect==1 .or. (ivect > 1 .and. ielem /= ielemone ))  then
                  if(tau(ivect)/=0.0_rp) tau(ivect)=1.0_rp/tau(ivect)
               end if
            end do elements45
   
            gptau(DEF_VECT,igaus) = tau(DEF_VECT)
         end do
         gptau_time(DEF_VECT,:) = gptau(DEF_VECT,:)

         ! Main loop
         do igaus = 1,pgaus
            !
            ! Calculus of residual resid and perturbation function Pi=gppre
            !
            ! RESI1 = r1(u) =  rho/dt*Nj + rho*a.grad(Nj) + r1*Nj
            ! RESI2 = r2(u) =  -grad(k).grad(Nj) - k*Lap(Nj)
            ! GPPE1 = p1(v) =  [ Ni*(1-tau*s) + tau*rho*a.grad(Ni) ]
            ! GPPE2 = p2(v) =  [ -Ni*tau*s  + tau*rho*a.grad(Ni) ] = p1(v) - v * |dv|
            !
            tau(DEF_VECT) = gptau_time(DEF_VECT,igaus)
            fact1(DEF_VECT) = gpden(DEF_VECT,igaus) * dtinv_elem
            dd(DEF_VECT) = fact1(DEF_VECT) * sgs(DEF_VECT,igaus)  ! rho*u'n/dt

            do inode = 1,pnode
               resi1(DEF_VECT,inode) = fact1(DEF_VECT) * ADR_chm % time_parameters(1) * gpsha(inode,igaus)
               gpad1(DEF_VECT,inode) = 0.0_rp
               grvgr(DEF_VECT) = 0.0_rp
               gplap(DEF_VECT) = 0.0_rp

               do idime = 1,ndime
                  gpad1(DEF_VECT,inode)  = gpad1(DEF_VECT,inode) + gpvel(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
                  grvgr(DEF_VECT) = grvgr(DEF_VECT) + gpgrd(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
                  gplap(DEF_VECT) = gplap(DEF_VECT) + gphes(DEF_VECT,idime,inode,igaus)
               end do

               gpadv(DEF_VECT,inode) = gpad1(DEF_VECT,inode) * gpden(DEF_VECT,igaus)
               resi1(DEF_VECT,inode) = resi1(DEF_VECT,inode) + gpadv(DEF_VECT,inode) + react(DEF_VECT,igaus) * gpsha(inode,igaus)
               resi2(DEF_VECT,inode) = - grvgr(DEF_VECT) - gplap(DEF_VECT) * gpdif(DEF_VECT,igaus)
               gppe1(DEF_VECT,inode) = ( gpsha(inode,igaus)*(1.0_rp-tau(DEF_VECT)*sreac(DEF_VECT,igaus)) + tau(DEF_VECT)&
                  * gpadv(DEF_VECT,inode) ) * gpvol(DEF_VECT,igaus)
               gppe2(DEF_VECT,inode) = gppe1(DEF_VECT,inode) - gpsha(inode,igaus) * gpvol(DEF_VECT,igaus)
            end do

            !
            ! Diffusion term
            !
            fact2(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpdif(DEF_VECT,igaus)
            do inode = 1,pnode
               do jnode = 1,inode-1
                  xmuit(DEF_VECT) = 0.0_rp
                  do idime = 1,ndime
                     xmuit(DEF_VECT) = xmuit(DEF_VECT) + gpcar(DEF_VECT,idime,jnode,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
                  end do
                  xmuit(DEF_VECT) = xmuit(DEF_VECT) * fact2(DEF_VECT)
                  elmat(DEF_VECT,inode,jnode) = elmat(DEF_VECT,inode,jnode) + xmuit(DEF_VECT)
                  elmat(DEF_VECT,jnode,inode) = elmat(DEF_VECT,jnode,inode) + xmuit(DEF_VECT)
               end do
               xmuit(DEF_VECT) = 0.0_rp
               do idime = 1,ndime
                  xmuit(DEF_VECT) = xmuit(DEF_VECT) + gpcar(DEF_VECT,idime,inode,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
               end do
               elmat(DEF_VECT,inode,inode) = elmat(DEF_VECT,inode,inode) + xmuit(DEF_VECT) * fact2(DEF_VECT)
            end do

            !
            ! bemol
            !
            do inode = 1,pnode
               fact1(DEF_VECT) = gpsha(inode,igaus) * ADR_chm % bemol * gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus)
               do jnode = 1,pnode
                  fact2(DEF_VECT) = gpadv(DEF_VECT,jnode) * fact1(DEF_VECT)
                  elmat(DEF_VECT,inode,jnode) = elmat(DEF_VECT,inode,jnode) - fact2(DEF_VECT)
                  elmat(DEF_VECT,jnode,inode) = elmat(DEF_VECT,jnode,inode) - fact2(DEF_VECT)
               end do
            end do

            !
            ! Assembly of the matrix and rhs
            !
            ! GPPR1 <=> f - sum_i ri*u^i
            ! GPPR1 <=> - [ rho * a - grad(k) ] . grad(u)
            ! GPPR1 <=> f - rho*(u - u^n)/dt - sum_i ri*u^i - [rho*u - grad(k)].grad(u) - k*lapl(u)
            !
            do inode = 1,pnode
               do jnode = 1,pnode
                  elmat(DEF_VECT,inode,jnode) = elmat(DEF_VECT,inode,jnode) + resi1(DEF_VECT,jnode) * gppe1(DEF_VECT,inode)&
                                                &+ resi2(DEF_VECT,jnode) * gppe2(DEF_VECT,inode)
               end do
               elrhs(:,inode) = elrhs(DEF_VECT,inode) + rhsit(DEF_VECT,igaus) * gppe1(DEF_VECT,inode)&
                                &+ ( dd(DEF_VECT) - gppr1(DEF_VECT,igaus) ) * gppe2(DEF_VECT,inode)
            end do

            if( ADR_chm % kfl_time_sgs /= 0 ) then ! Time tracking of subgrid scale
               call runend('CHEMIC DOES NOT USE TIME TRACKING OF SUBGRID SCALE')
            end if

            if( ADR_chm % kfl_time_lumped == 1 ) then  ! Lumped mass evolution matrix
               call runend('CHEMIC DOES NOT USE LUMPED MASS EVOLUTION MATRIX')
            end if
         end do
      else
         call runend('CHEMIC ONLY USES GALERKIN OR ASGS STABILIZATION METHODS')
      end if
   end subroutine chm_element_assembly


  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Assemble a RHS
  !> @details Assemble a RHS
  !
  !----------------------------------------------------------------------

  subroutine chm_matrix_assexp_fast(VECTOR_DIM,ndofn_glo,ndofn_ele,pnode,lnods_loc,elrhs,elmat,elunk,rhsid,kdofn,list_elements_p)
    implicit none
    integer(ip),  intent(in)           :: VECTOR_DIM
    integer(ip),  intent(in)           :: ndofn_glo                              !< Number of dof of global RHS
    integer(ip),  intent(in)           :: ndofn_ele                              !< Number of dof of element matrix and RHS
    integer(ip),  intent(in)           :: pnode                                  !< Number of nodes of element
    integer(ip),  intent(in)           :: lnods_loc(VECTOR_DIM,pnode)                           !< Element connectivity
    real(rp),     intent(in)           :: elrhs(VECTOR_DIM,pnode)                               !< Element RHS
    real(rp),     intent(in)           :: elmat(VECTOR_DIM,pnode*ndofn_ele,pnode*ndofn_ele) !< Element matrix
    real(rp),     intent(in)           :: elunk(VECTOR_DIM,pnode*ndofn_ele)                 !< Element unknwon
    real(rp),     intent(inout)        :: rhsid(*)                               !< Global RHS
    integer(ip),  intent(in), optional :: kdofn                                  !< Local dof to assemble
    integer(ip),  intent(in)           :: list_elements_p(VECTOR_DIM)           ! List of elements (always positive)
    integer(ip)                        :: inode,ipoin,idofg
    integer(ip)                        :: jnode
    integer(ip)                        :: ielem,ivect,ielemone

    external                           :: runend

    ielemone  = list_elements_p(1)

    if( ndofn_glo > 1 .and. ndofn_ele == 1 .and. ( .not. present(kdofn) ) ) then
       !
       ! NDOFN_GLO DOF's but only KDOFN is assembled
       !
       call runend('MATRIX_ASSEXP: KDOFN IS MISSING')

    else

       if( ndofn_glo == 1 ) then
          !
          ! 1 DOF
          !
          do ivect=1, VECTOR_DIM
              ielem = list_elements_p(ivect)

              if( ivect == 1 .or. (ivect>1 .and. ielem /= ielemone)) then

                  do inode = 1,pnode
                     ipoin = lnods_loc(ivect,inode)
                     rhsid(ipoin) = rhsid(ipoin) + elrhs(ivect,inode)
                     do jnode = 1,pnode
#ifdef NO_COLORING
                        !$OMP ATOMIC
#endif
                        rhsid(ipoin) = rhsid(ipoin) - elmat(ivect,inode,jnode) * elunk(ivect,jnode)
                     end do
                  end do
               endif
          end do


       else if( ndofn_glo > 1 ) then
          !
          ! DOF > 1
          !
          if( ndofn_ele == 1 ) then


              do ivect=1, VECTOR_DIM
                  ielem = list_elements_p(ivect)

                  if( ivect == 1 .or. (ivect>1 .and. ielem /= ielemone)) then




                      do inode = 1,pnode
                        ipoin = lnods_loc(ivect,inode)
                        idofg = (ipoin-1)*ndofn_glo+kdofn
                        rhsid(idofg) = rhsid(idofg) + elrhs(ivect,inode)
                        do jnode = 1,pnode
#ifdef NO_COLORING
                           !$OMP ATOMIC
#endif
                           rhsid(idofg) = rhsid(idofg) - elmat(ivect,inode,jnode) * elunk(ivect,jnode)
                        end do
                     end do
                  end if
              end do
          else
             call runend('NOT CODED, ESPECIE DE VAGO')
          end if

       end if
    end if

  end subroutine chm_matrix_assexp_fast



end module mod_chm_element_operations_fast
!> @}
