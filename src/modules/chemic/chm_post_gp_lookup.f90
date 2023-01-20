!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_post_gp_lookup(aux_r2p,fw)
  !------------------------------------------------------------------------
  ! lookup from postprocessing table on Gauss points
  !------------------------------------------------------------------------
  use def_chemic,     only  : ADR_chm, nclas_chm, kfl_table_ind_chm
  use def_master,     only  : therm
  use def_kintyp,     only  : ip,rp,r2p
  use def_domain,     only  : nelem,mgaus,mnode
  use def_domain,     only  : lnods,elmar,ndime
  use def_domain,     only  : ltype,ngaus,nnode
  use mod_interp_tab, only  : fw_lookup,typ_lookup_framework
  implicit none

  type(r2p),                   intent(inout) :: aux_r2p(nelem)
  type(typ_lookup_framework),  intent(in)    :: fw

  integer(ip)              :: ipoin,ipostvar,iclas,idimt,ielem,igaus,inode
  integer(ip)              :: pelty,pgaus,pnode
  real(rp)                 :: retva(fw % main_table % nvar)
  real(rp)                 :: control(fw % main_table % ndim)
  integer(ip)              :: ind(fw % main_table % ndim)
  real(rp)                 :: scale_control(fw % main_table % ndim)
  integer(ip)              :: lnods_loc(mnode)
  real(rp)                 :: elcon(mnode,nclas_chm,ADR_chm(1) % ntime)
  real(rp)                 :: gpcon(mgaus,nclas_chm)
  real(rp)                 :: gpthe(mgaus)
  real(rp)                 :: elcod(ndime,mnode)
  real(rp)                 :: dummr(ndime,mnode)

  external                 :: chm_post_gather
  external                 :: runend

  !
  ! Set variable to zero
  !
  do ielem = 1,nelem
     pelty = ltype(ielem)
     if( pelty > 0 ) then
        pgaus = ngaus(pelty)
        do ipostvar=1,fw % main_table % nvar
           do igaus = 1,pgaus
              aux_r2p(ielem)%a(igaus,ipostvar) = 0.0_rp
           enddo
        enddo
     endif
  enddo

  !
  ! Loop over elements
  !
  ind = 1_ip
  elements: do ielem = 1,nelem
     pelty = ltype(ielem)
     if( pelty > 0 ) then
        pgaus = ngaus(pelty)
        pnode = nnode(pelty)
        !
        ! Gather all
        !
        lnods_loc(1:pnode) = lnods(1:pnode,ielem)
        call chm_post_gather(&
             pnode,lnods_loc,elcon(1:pnode,:,:),elcod,dummr)

        !
        ! Initialization variables
        !
        gpcon = 0.0_rp
        gpthe = 0.0_rp

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

        if (fw % kfl_needs_enthalpy /= 0) then
           do igaus = 1,pgaus
              do inode = 1,pnode
                 ipoin=lnods_loc(inode)
                 gpthe(igaus) = gpthe(igaus)&
                      +  elmar(pelty)%shape(inode,igaus) * therm(ipoin,1)
              end do
           end do
        endif

        !
        ! Lookup
        !
        do igaus = 1,pgaus
           control = 0.0_rp
           do idimt = 1, fw % main_table % ndim
              if (fw % kfl_chm_control(idimt) > 0) then
                 !
                 ! >0: one of the conces
                 !
                 control(idimt) = gpcon(igaus,fw % kfl_chm_control(idimt))
              else
                 if (fw % kfl_chm_control(idimt) == -1) then
                    !
                    ! -1: enthalpy
                    !
                    control(idimt) = gpthe(igaus)
                 elseif (fw % kfl_chm_control(idimt) == -2) then
                    !
                    ! -2: scalar dissipation rate
                    !
                    call runend('chm_post_gp_lookup: table post-processing is not implemented for UFPV model.')
                 endif
              endif
           enddo

           call fw_lookup( control, scale_control, fw, retva, kfl_table_ind_chm(1:fw%main_table%ndim,ielem) )

           do ipostvar=1,fw % main_table % nvar
              aux_r2p(ielem)%a(igaus,ipostvar) = retva(ipostvar)
           enddo
        enddo
     endif
  enddo elements


end subroutine chm_post_gp_lookup

