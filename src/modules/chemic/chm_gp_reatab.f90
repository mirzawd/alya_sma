!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_gp_reatab(itask)
  !-----------------------------------------------------------------------
  !****f* chemic/chm_reatab
  ! NAME
  !    chm_reatab
  ! DESCRIPTION
  !    Read table properties for Flamelet combustion model
  ! USES
  ! USED BY
  !    chm_iniunk: initialize properties
  !    chm_begste: update source terms
  !    chm_endite: update properties
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only          : ip,rp
  use def_master, only          : conce,sphec,therm,ITASK_BEGSTE,ITASK_ENDITE
  use def_master, only          : condu_gp,sphec_gp,visco_gp,wmean_gp,&
                                  sphec_gp_ht,sphec_gp_lt,tempe_gp
  use def_chemic, only          : mass_gp,massConsumption_gp,mixedEq_eqs_chm,&
                                  Yk_ssm_gp,kfl_yk_fw_ssm

  use def_domain, only          : npoin,nnode,nelem,ltype,&
                                  ngaus,llapl,lorde,ltopo,elmar,hnatu,lnods,&
                                  mgaus,ndime,mnode,ntens,lpoty
  use def_chemic, only          : kfl_radia_chm, &
                                  kfl_ufpv_chm, xZr_chm, xZs_chm, &
                                  kfl_advec_chm,kfl_ellen_chm,ADR_chm,nclas_chm,&
                                  kfl_table_ind_chm,&
                                  kfl_min_srcfw_chm,kfl_max_srcfw_chm,&
                                  kfl_fw_has_sources, kfl_fw_src_equa_list
  use def_chemic,      only     : kfl_W_tab_index_chm
  use def_chemic,      only     : kfl_k_tab_index_chm
  use def_chemic,      only     : kfl_mu_tab_index_chm
  use def_chemic,      only     : kfl_cpCoefLT_tab_index_chm
  use def_chemic,      only     : kfl_cpCoefHT_tab_index_chm
  use def_chemic,      only     : kfl_cpCoefLT_end_tab_index_chm
  use def_chemic,      only     : kfl_cpCoefHT_end_tab_index_chm
  use def_chemic,      only     : kfl_T_tab_index_chm
  use def_chemic,      only     : kfl_DtRho_tab_index_chm
  use def_chemic,      only     : kfl_max_nvar_lookup_chm
  use def_chemic,      only     : table_fw
  use def_chemic,      only     : DtRho_gp
  use mod_ker_proper,  only     : ker_proper
  use mod_interp_tab,  only     : fw_lookup
  use mod_interp_tab,  only     : tab_interp
  use mod_interp_tab,  only     : fw_scale_cont_var

  use def_kermod,          only : lookup_fw

  implicit none
  integer(ip), intent(in)   :: itask
  integer(ip)               :: ipoin,iclas,idimt,iequa,ii,ifw
  real(rp)                  :: retva(kfl_max_nvar_lookup_chm)                   ! Values read in from flamelet table:
  !                                                                               (1) S_c (2) Wmean (3) Lambda (4) Mu (5-10) cp_low
  !                                                                               (11-16) cp_high (17) S_c*c
  real(rp)                  :: tab_scale_control(table_fw % main_table % ndim)  ! Concentration for each variable at each node
  real(rp)                  :: lim_control(table_fw % main_table % ndim,2_ip)   ! Limiting values of control variables
  real(rp)                  :: control(table_fw % main_table % ndim)            ! input of table lookup function
  integer(ip)               :: ind(table_fw % main_table % ndim)

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

  external    :: chm_post_gather
  external    :: elmlen
  external    :: elmchl
  external    :: elmcar
  external    :: chm_calc_scalar_dissip
  external    :: runend

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
           mass_gp(ielem) % a             = 0.0_rp
           massConsumption_gp(ielem) % a  = 0.0_rp
        case(ITASK_ENDITE)
           if (kfl_k_tab_index_chm > 0)        condu_gp(ielem) % a      = 0.0_rp
           if (kfl_mu_tab_index_chm > 0)       visco_gp(ielem) % a      = 0.0_rp
           if (kfl_W_tab_index_chm > 0)        wmean_gp(ielem) % a (1:pgaus,1,1)     = 0.0_rp ! Only put current value to 0
                                                                                              !  We need to save old steps
           if (kfl_cpCoefHT_tab_index_chm > 0) sphec_gp_ht(ielem) % a   = 0.0_rp
           if (kfl_cpCoefLT_tab_index_chm > 0) sphec_gp_lt(ielem) % a   = 0.0_rp
           if (kfl_cpCoefHT_tab_index_chm > 0 .or. kfl_cpCoefLT_tab_index_chm > 0)   sphec_gp(ielem) % a = 0.0_rp
           if (kfl_T_tab_index_chm > 0)        tempe_gp(ielem) % a (1:pgaus,1,1)     = 0.0_rp ! Only put current value to 0
                                                                                              !  We need to save old steps
           if (kfl_DtRho_tab_index_chm > 0)    DtRho_gp(ielem) % a      = 0.0_rp
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

        if (table_fw % kfl_needs_enthalpy /= 0) then
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
           do idimt = 1, table_fw % main_table % ndim
              if (table_fw % kfl_chm_control(idimt) > 0) then
                 !
                 ! >0: one of the conces
                 !
                 control(idimt) = gpcon(igaus,table_fw % kfl_chm_control(idimt))
              else
                 if (table_fw % kfl_chm_control(idimt) == -1) then
                    !
                    ! -1: enthalpy
                    !
                    control(idimt) = gpthe(igaus)
                 elseif (table_fw % kfl_chm_control(idimt) == -2) then
                    !
                    ! -2: scalar dissipation rate
                    !
                    control(idimt) = gp_res_Chi_z(igaus) + gp_sgs_Chi_z(igaus)
                 endif
              endif
           enddo

           !
           ! Sale control variables
           !
           call fw_scale_cont_var( control, tab_scale_control, lim_control, &
               table_fw, kfl_table_ind_chm(1:table_fw%main_table%ndim,ielem) )

           select case(itask)

           case(ITASK_BEGSTE)
              !
              ! Loop though all frameworks potentially involved in source terms
              !
              do ifw = kfl_min_srcfw_chm,kfl_max_srcfw_chm
                 !
                 ! If the number of equations depending on the framework is non-zero
                 ! then do the lookup
                 !
                 if (kfl_fw_has_sources(ifw) > 0) then
                     call tab_interp( lookup_fw(ifw) % main_table, tab_scale_control, retva(1:lookup_fw(ifw)%main_table%nvar),&
                         kfl_table_ind_chm(1:table_fw%main_table%ndim,ielem) )

                     !
                     ! Go through all the equations related to the framework
                     !
                     do ii = 1,kfl_fw_has_sources(ifw)
                        iequa = kfl_fw_src_equa_list(ifw,ii)
                        mass_gp(ielem) % a(igaus,iequa,1)   = retva(mixedEq_eqs_chm(iequa) % kfl_source_col)
                        if ( mixedEq_eqs_chm(iequa) % kfl_source_split == 1) &
                           massConsumption_gp(ielem) % a(igaus,iequa,1)  = retva(mixedEq_eqs_chm(iequa) % kfl_consum_col)
                     enddo

                 endif
              enddo

              !
              ! Look up tabulated species if necessary
              !
              if (kfl_yk_fw_ssm > 0) then
                 call tab_interp( lookup_fw(kfl_yk_fw_ssm) % main_table, tab_scale_control, Yk_ssm_gp(ielem) % a(igaus,&
                     1:lookup_fw(kfl_yk_fw_ssm)%main_table%nvar,1), kfl_table_ind_chm(1:table_fw%main_table%ndim,ielem) )
              endif


           case(ITASK_ENDITE)
              !
              ! Lookup properties for kernel
              !
              call tab_interp( table_fw % main_table, tab_scale_control, retva(1:table_fw%main_table%nvar),&
                  kfl_table_ind_chm(1:table_fw%main_table%ndim,ielem) )

              if (kfl_W_tab_index_chm > 0)        wmean_gp(ielem) % a(igaus,1,1)      = retva(kfl_W_tab_index_chm)
              if (kfl_k_tab_index_chm > 0)        condu_gp(ielem) % a(igaus,1,1)      = retva(kfl_k_tab_index_chm)
              if (kfl_mu_tab_index_chm > 0)       visco_gp(ielem) % a(igaus,1,1)      = retva(kfl_mu_tab_index_chm)

              if (kfl_cpCoefLT_tab_index_chm > 0) sphec_gp_lt(ielem) % a(igaus,1:6,1) = retva(kfl_cpCoefLT_tab_index_chm:&
                  kfl_cpCoefLT_end_tab_index_chm)
              if (kfl_cpCoefHT_tab_index_chm > 0) sphec_gp_ht(ielem) % a(igaus,1:6,1) = retva(kfl_cpCoefHT_tab_index_chm:&
                  kfl_cpCoefHT_end_tab_index_chm)

              if (kfl_T_tab_index_chm > 0)        tempe_gp(ielem) % a(igaus,1,1)      = retva(kfl_T_tab_index_chm)
              if (kfl_DtRho_tab_index_chm > 0)    DtRho_gp(ielem) % a(igaus,1,1)      = retva(kfl_DtRho_tab_index_chm)

              if (kfl_radia_chm > 0) then
                 call runend('Chemic gp_reatab: Optically thin radiation model is not implemented.')
              end if
           end select
        enddo   !igaus loop

     endif
  end do elements

  !
  ! LOOKUP Cp coefficients on boundary for temperature BC:
  !
  ind = 1_ip
  bou_nodes: do ipoin = 1,npoin
     !
     ! Check if boundary node
     !
     if (lpoty(ipoin) /= 0) then
        control = 0.0_rp

        do idimt = 1, table_fw % main_table % ndim
           if (table_fw % kfl_chm_control(idimt) > 0) then
              !
              ! >0: one of the conces
              !
              control(idimt) = conce(ipoin,table_fw % kfl_chm_control(idimt),1)
           else
              if (table_fw % kfl_chm_control(idimt) == -1) then
                 !
                 ! -1: enthalpy
                 !
                 control(idimt) = therm(ipoin,1)
              elseif (table_fw % kfl_chm_control(idimt) == -2) then
                 !
                 ! -2: scalar dissipation rate
                 !
                 control(idimt) = xZr_chm(ipoin) + xZs_chm(ipoin)
              endif
           endif
        enddo

        call fw_lookup( control, tab_scale_control, table_fw, retva(1:table_fw%main_table%nvar), ind )

        if (kfl_cpCoefLT_tab_index_chm > 0) sphec(ipoin,1:6,1) = retva(kfl_cpCoefLT_tab_index_chm:kfl_cpCoefLT_end_tab_index_chm)
        if (kfl_cpCoefHT_tab_index_chm > 0) sphec(ipoin,1:6,2) = retva(kfl_cpCoefHT_tab_index_chm:kfl_cpCoefHT_end_tab_index_chm)
     endif
  enddo bou_nodes

end subroutine chm_gp_reatab
