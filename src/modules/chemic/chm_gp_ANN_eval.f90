!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_gp_ANN_eval(itask)
  !-----------------------------------------------------------------------
  !****f* chemic/chm_gp_ANN_eval
  ! NAME
  !    chm_reatab
  ! DESCRIPTION
  !    Evaluate ANN frameworks for Flamelet-like combustion model
  ! USES
  ! USED BY
  !    chm_iniunk: initialize properties
  !    chm_begste: update source terms
  !    chm_endite: update properties
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only          : ip,rp
  use def_master, only          : ITASK_BEGSTE,ITASK_ENDITE
  use def_master, only          : therm
  use def_master, only          : tempe_gp

  use def_chemic, only          : mass_gp,mixedEq_eqs_chm

  use def_domain, only          : nnode,nelem,ltype,&
                                  ngaus,llapl,lorde,ltopo,elmar,hnatu,lnods,&
                                  mgaus,ndime,mnode,ntens

  use def_chemic, only          : kfl_ufpv_chm, &
                                  kfl_advec_chm,kfl_ellen_chm,ADR_chm,nclas_chm,&
                                  kfl_min_src_annfw_chm,kfl_max_src_annfw_chm,&
                                  kfl_annfw_has_sources, kfl_annfw_src_equa_list
  use def_chemic,      only     : kfl_max_nvar_ann_in_chm
  use def_chemic,      only     : kfl_max_nvar_ann_out_chm
  use mod_ker_proper,  only     : ker_proper

  use def_kermod,      only     : ann_fw

  implicit none
  integer(ip), intent(in)   :: itask
  integer(ip)               :: ipoin,iclas,iequa,ii,ifw

  integer(ip) :: ielem,igaus,inode
  integer(ip) :: dummi
  integer(ip) :: pelty,pnode
  integer(ip) :: pgaus,plapl,porde,ptopo
  real(rp)    :: gphco(mgaus)
  real(rp)    :: gpsph(mgaus)
  real(rp)    :: gp_res_Chi_Yc(mgaus)
  real(rp)    :: gp_sgs_Chi_Yc(mgaus)
  real(rp)    :: gp_res_Chi_z(mgaus)
  real(rp)    :: gp_sgs_Chi_z(mgaus)
  real(rp)    :: gpcar(ndime,mnode,mgaus)
  real(rp)    :: gphes(ntens,mnode,mgaus)
  real(rp)    :: gpcon(mgaus,nclas_chm)
  real(rp)    :: gpthe(mgaus)
  real(rp)    :: gpden(mgaus)
  real(rp)    :: gptur(mgaus)
  real(rp)    :: gpvol(mgaus)
  real(rp)    :: chale(3),chave(3),hleng(3),tragl(9)
  real(rp)    :: elcon(mnode,nclas_chm,ADR_chm(1) % ntime)
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: dummr(ndime,mnode)
  integer(ip) :: lnods_loc(mnode)

  real(rp)    :: input(mgaus, kfl_max_nvar_ann_in_chm)
  real(rp)    :: output(mgaus, kfl_max_nvar_ann_out_chm)

  external    :: chm_post_gather
  external    :: elmlen
  external    :: elmchl
  external    :: elmcar
  external    :: chm_calc_scalar_dissip

  !
  ! Exit if not needed
  !
  select case(itask)
  case(ITASK_BEGSTE)
      !
      ! Return if ANNs are not used in begste
      !
      if (all( kfl_annfw_has_sources == 0 )) then
          return
      endif
  case(ITASK_ENDITE)
      !
      ! return? Be careful with this
      !
      return
  end select


  !
  ! Loop over elements
  !
  elements: do ielem = 1,nelem
     !
     ! Element dimensions
     !
     pelty = ltype(ielem)
     if( pelty > 0 ) then
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)
        plapl = llapl(pelty)
        porde = lorde(pelty)
        ptopo = ltopo(pelty)

        !
        ! Gather all
        !
        lnods_loc(1:pnode) = lnods(1:pnode,ielem)
        call chm_post_gather(&
             pnode,lnods_loc,elcon(1:pnode,:,:),elcod,dummr)


        if (kfl_ufpv_chm > 0) then
           !
           ! CHALE, HLENG and TRAGL
           !
           call elmlen(&
                ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),&
                hleng)
           call elmchl(&
                tragl,hleng,elcod,dummr,chave,chale,pelty,pnode,&
                porde,hnatu(pelty),kfl_advec_chm,kfl_ellen_chm)

           !
           ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, GPVOL
           !
           call elmcar(&
                pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,elmar(pelty)%deriv, &
                elmar(pelty)%heslo,elcod,gpvol,gpcar,gphes,ielem)


           call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,elmar(pelty)%shape,gpcar)
           call ker_proper('CONDU','PGAUS',dummi,ielem,gphco,pnode,pgaus,elmar(pelty)%shape,gpcar)
           call ker_proper('SPHEA','PGAUS',dummi,ielem,gpsph,pnode,pgaus,elmar(pelty)%shape,gpcar)
           call ker_proper('TURBU','PGAUS',dummi,ielem,gptur,pnode,pgaus,elmar(pelty)%shape,gpcar)

           !
           ! Compute scalar dissipation rate at gauss points
           !
           call chm_calc_scalar_dissip(&
                pnode,pgaus,elcon(1:pnode,:,:),dummr,elmar(pelty)%shape,gpcar,hleng,gptur,&
                gphco,gpsph,gpden,gp_res_Chi_Yc,gp_sgs_Chi_Yc,gp_res_Chi_z,gp_sgs_Chi_z)
        endif

        !
        ! Initialization variables
        !
        gpcon = 0.0_rp
        gpthe = 0.0_rp
        select case(itask)
        case(ITASK_BEGSTE)
           do ifw = kfl_min_src_annfw_chm,kfl_max_src_annfw_chm
              if (kfl_annfw_has_sources(ifw) > 0) then
                  do ii = 1,kfl_annfw_has_sources(ifw)
                     iequa = kfl_annfw_src_equa_list(ifw,ii)
                     mass_gp(ielem) % a(1:pgaus,iequa,1)   = 0.0_rp
                  enddo
              endif
           enddo
!        case(ITASK_ENDITE)
        end select

        !
        ! Concentration
        !
        do iclas = 1,nclas_chm
           do igaus = 1,pgaus
              do inode = 1,pnode
                 gpcon(igaus,iclas) = gpcon(igaus,iclas)&
                      + elmar(pelty)%shape(inode,igaus) * elcon(inode,iclas,1)
              end do
           end do
        end do

        if (associated(therm)) then
           do igaus = 1,pgaus
              do inode = 1,pnode
                 ipoin=lnods_loc(inode)
                 gpthe(igaus) = gpthe(igaus)&
                      +  elmar(pelty)%shape(inode,igaus) * therm(ipoin,1)
              end do
           end do
        endif


        !
        ! Evaluate frameworks
        !
        input = 0.0_rp
        output = 0.0_rp


        select case(itask)
        case(ITASK_BEGSTE)
           !
           ! Loop though all frameworks potentially involved in source terms
           !
           do ifw = kfl_min_src_annfw_chm,kfl_max_src_annfw_chm
              !
              ! If the number of equations depending on the framework is non-zero
              ! then do the lookup
              !
              if (kfl_annfw_has_sources(ifw) > 0) then

                  !
                  ! Assemble input
                  !
                  do ii = 1, ann_fw(ifw) % scaling_in % n_prod_shape
                     if (ann_fw(ifw) % scaling_in % iaux1(ii) > 0) then
                        !
                        ! >0: one of the conces
                        !
                        input(1:pgaus,ii) = gpcon(1:pgaus,ann_fw(ifw) % scaling_in % iaux1(ii))
                     else
                        if (ann_fw(ifw) % scaling_in % iaux1(ii) == -1) then
                           !
                           ! -1: enthalpy
                           !
                           input(1:pgaus,ii) = gpthe(1:pgaus)
                        elseif (ann_fw(ifw) % scaling_in % iaux1(ii) == -2) then
                           !
                           ! -2: scalar dissipation rate
                           !
                           input(1:pgaus,ii) = gp_res_Chi_z(1:pgaus) + gp_sgs_Chi_z(1:pgaus)
                        elseif (ann_fw(ifw) % scaling_in % iaux1(ii) == -3) then
                           !
                           ! -3: temperature
                           !
                           input(1:pgaus,ii) = tempe_gp(ielem) % a(1:pgaus,1,1)
                        endif
                     endif
                  enddo

                  !
                  ! Forward pass of neural network
                  !
                  call ann_fw(ifw) % eval(nvec   = pgaus,                                                     &
                      &                   input  = input(1:pgaus,1:ann_fw(ifw) % scaling_in % n_prod_shape),  &
                      &                   output = output(1:pgaus,1:ann_fw(ifw) % scaling_out % n_prod_shape))

                  !
                  ! Go through all the equations related to the framework
                  !
                  do ii = 1,kfl_annfw_has_sources(ifw)
                     iequa = kfl_annfw_src_equa_list(ifw,ii)
                     mass_gp(ielem) % a(1:pgaus,iequa,1)   = output(1:pgaus, mixedEq_eqs_chm(iequa) % kfl_source_col)
                  enddo

              endif
           enddo

        case(ITASK_ENDITE)
           !
           ! Lookup properties for kernel
           !
           !call runend('Chemic chm_gp_ANN_eval: Property evaluation is not implemented.')
        end select

     endif
  end do elements

end subroutine chm_gp_ANN_eval
