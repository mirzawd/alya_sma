!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_gp_multi_reatab(itask)
  !-----------------------------------------------------------------------
  !****f* chemic/chm_gp_multi_reatab
  ! NAME 
  !    chm_reatab
  ! DESCRIPTION
  !    Read table properties for Flamelet Multi-regime combustion model
  ! USES
  ! USED BY
  !    chm_iniunk: initialize properties 
  !    chm_begste: update source terms
  !    chm_endite: update properties
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only          : ip,rp
  use def_master, only          : conce,sphec,therm,ITASK_BEGSTE,ITASK_ENDITE
  use def_kermod, only          : lookup_fw
  use def_master, only          : condu_gp,sphec_gp,visco_gp,wmean_gp,&
                                  sphec_gp_ht,sphec_gp_lt,tempe_gp
  use def_chemic, only          : mass_gp,massConsumption_gp,mixedEq_eqs_chm,&
                                  kfl_izmean_chm,kfl_icmean_chm,&
                                  flame_index_gp,zgrad_gp, ittot_chm

  use def_domain, only          : npoin,nnode,nelem,ltype,&
                                  ngaus,llapl,lorde,ltopo,elmar,lnods,&
                                  mgaus,ndime,mnode,ntens,lpoty
  use def_chemic, only          : kfl_radia_chm, &
                                  xZr_chm, xZs_chm, &
                                  ADR_chm,nclas_chm,&
                                  kfl_min_srcfw_chm,kfl_max_srcfw_chm,&
                                  kfl_fw_src_equa_list, &
                                  kfl_tab_fw_chm_diff, kfl_tab_fw_chm_prem, &
                                  kfl_fw_has_hybrid_sources,kfl_zg_fw_chm
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
  use def_chemic,      only     : DtRho_gp
  use def_chemic,      only     : flame_index_gp
  use mod_ker_proper,  only     : ker_proper
  use mod_interp_tab,  only     : fw_lookup 
  use mod_interp_tab,  only     : tab_interp 
  use mod_interp_tab,  only     : tab_interpt 
  use mod_interp_tab,  only     : fw_scale_cont_var

  implicit none
  integer(ip), intent(in)   :: itask
  integer(ip)               :: ipoin,iclas,idimt,iequa,ii,ifw
  real(rp)                  :: retva1(kfl_max_nvar_lookup_chm)                   ! Values read in from flamelet table: 
  real(rp)                  :: retva2(kfl_max_nvar_lookup_chm)                   ! Values read in from flamelet table: 
  real(rp)                  :: retvapost(lookup_fw(kfl_zg_fw_chm) &
                               % main_table % nvar)                               ! Values read in from flamelet table: 
  !                                                                               (1) S_c (2) Wmean (3) Lambda (4) Mu (5-10) cp_low
  !                                                                               (11-16) cp_high (17) S_c*c
  real(rp)                  :: tab_scale_control1(lookup_fw(kfl_tab_fw_chm_diff) & 
                               % main_table % ndim) ! Concentration for each variable at each node
  real(rp)                  :: control(lookup_fw(kfl_tab_fw_chm_diff) &
                               % main_table % ndim)  ! input of table lookup function 
  integer(ip)               :: ind(lookup_fw(kfl_tab_fw_chm_diff) % main_table % ndim) 

  integer(ip) :: ielem,igaus,inode
  integer(ip) :: pelty,pnode
  integer(ip) :: pgaus,plapl,porde,ptopo
  real(rp)    :: gpcar(ndime,mnode,mgaus)                 
  real(rp)    :: gphes(ntens,mnode,mgaus)
  real(rp)    :: gpcon(mgaus,nclas_chm)
  real(rp)    :: gpthe(mgaus)
  real(rp)    :: gpvol(mgaus)
  real(rp)    :: elcon(mnode,nclas_chm,ADR_chm(1) % ntime)
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: dummr(ndime,mnode)
  real(rp)    :: wDiff,wPrem
  integer(ip) :: lnods_loc(mnode)

  external    :: chm_post_gather
  external    :: elmlen
  external    :: elmchl
  external    :: elmcar
  external    :: chm_calc_scalar_dissip
  external    :: chm_calc_flame_index
  external    :: runend

  !
  ! Loop over elements
  !
  ind = 1_ip
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

        !
        ! Initialization variables
        !
        gpcon = 0.0_rp 
        gpthe = 0.0_rp 
        select case(itask)
        case(ITASK_BEGSTE)
           do iequa = 1,nclas_chm 
              if (mixedEq_eqs_chm(iequa) % kfl_diffsource_fw /=0_ip .or. mixedEq_eqs_chm(iequa) % kfl_premsource_fw /=0_ip) then
                 mass_gp(ielem) % a(:,iequa,:)             = 0.0_rp
                 massConsumption_gp(ielem) % a(:,iequa,:)  = 0.0_rp
              endif
           end do
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

        if (lookup_fw(kfl_tab_fw_chm_diff) % kfl_needs_enthalpy /= 0) then
           do igaus = 1,pgaus
              do inode = 1,pnode
                 ipoin=lnods_loc(inode)
                 gpthe(igaus) = gpthe(igaus)&
                      +  elmar(pelty)%shape(inode,igaus) * therm(ipoin,1)
              end do
           end do
        endif

        !
        ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, GPVOL
        !
        call elmcar(&
             pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,elmar(pelty)%deriv, &
             elmar(pelty)%heslo,elcod,gpvol,gpcar,gphes,ielem)
        
        if (ittot_chm > 0_ip) then 
           call chm_calc_flame_index(pnode,pgaus,elcon(1:pnode,kfl_izmean_chm,1),&
                elcon(1:pnode,kfl_icmean_chm,1),gpcar,flame_index_gp(ielem)%a&
                (1:pgaus,1),zgrad_gp(ielem)%a(1:pgaus,1),gpcon(1:pgaus,kfl_izmean_chm)&
                ,gpcon(1:pgaus,kfl_icmean_chm))
        else
           flame_index_gp(ielem)%a(1:pgaus,1) = 0.0_rp
        end if 
        
        !
        ! Lookup
        !
        do igaus = 1,pgaus
           !
           ! Weights for lookup in premixed and diffusion tables
           !
           wPrem = flame_index_gp(ielem)%a(igaus,1) 
           wDiff = 1.0_rp - flame_index_gp(ielem)%a(igaus,1)
           control = 0.0_rp
           do idimt = 1, lookup_fw(kfl_tab_fw_chm_diff) % main_table % ndim
              if (lookup_fw(kfl_tab_fw_chm_diff) % kfl_chm_control(idimt) > 0) then
                 !
                 ! >0: one of the conces
                 !
                 control(idimt) = gpcon(igaus,lookup_fw(kfl_tab_fw_chm_diff) % kfl_chm_control(idimt))
              else
                 if (lookup_fw(kfl_tab_fw_chm_diff) % kfl_chm_control(idimt) == -1) then
                    !
                    ! -1: enthalpy
                    !
                    control(idimt) = gpthe(igaus)
                 endif
              endif
           enddo

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
                 if (kfl_fw_has_hybrid_sources(ifw) > 0) then
                     
                     call fw_lookup( control, tab_scale_control1, lookup_fw(ifw),&
                          retva1(1:lookup_fw(ifw)%main_table%nvar),ind)
                     !
                     ! Go through all the equations related to the framework
                     !
                     do ii = 1,kfl_fw_has_hybrid_sources(ifw)
                        iequa = kfl_fw_src_equa_list(ifw,ii)
                        if (mixedEq_eqs_chm(iequa) % kfl_diffsource_fw ==ifw) then
                           mass_gp(ielem) % a(igaus,iequa,1) = mass_gp(ielem) % a(igaus,iequa,1) + &
                           retva1(mixedEq_eqs_chm(iequa) % kfl_diffsource_col) * wdiff 
                        endif
                        if (mixedEq_eqs_chm(iequa) % kfl_premsource_fw ==ifw) then
                           mass_gp(ielem) % a(igaus,iequa,1) = mass_gp(ielem) % a(igaus,iequa,1) + &
                           retva1(mixedEq_eqs_chm(iequa) % kfl_premsource_col) * wprem
                        endif
                     enddo
                 endif
              enddo
              ! 
              ! Look up Zgrad
              !
              call fw_lookup(control,tab_scale_control1,lookup_fw(kfl_zg_fw_chm),retvapost,ind)

              zgrad_gp(ielem)%a(igaus,1) = retvapost(1)
  
           case(ITASK_ENDITE)
              !
              ! Lookup properties for kernel with both the tables 
              !
              call fw_lookup(control,tab_scale_control1,lookup_fw(kfl_tab_fw_chm_diff),&
                             retva1(1:lookup_fw(kfl_tab_fw_chm_diff)%main_table%nvar),ind)
              call fw_lookup(control,tab_scale_control1,lookup_fw(kfl_tab_fw_chm_prem),&
                             retva2(1:lookup_fw(kfl_tab_fw_chm_prem)%main_table%nvar),ind)

              if (kfl_W_tab_index_chm > 0)        wmean_gp(ielem) % a(igaus,1,1)      = retva2(kfl_W_tab_index_chm) * wPrem &
                 + retva1(kfl_W_tab_index_chm) * wDiff
              if (kfl_k_tab_index_chm > 0)        condu_gp(ielem) % a(igaus,1,1)      = retva2(kfl_k_tab_index_chm) * wPrem &
                 + retva1(kfl_k_tab_index_chm) * wDiff
              if (kfl_mu_tab_index_chm > 0)       visco_gp(ielem) % a(igaus,1,1)      = retva2(kfl_mu_tab_index_chm) * wPrem &
                 + retva1(kfl_mu_tab_index_chm) * wDiff
              if (kfl_cpCoefLT_tab_index_chm > 0) sphec_gp_lt(ielem) % a(igaus,1:6,1) = &
                  retva2(kfl_cpCoefLT_tab_index_chm:kfl_cpCoefLT_end_tab_index_chm) * wPrem + &
                  retva1(kfl_cpCoefLT_tab_index_chm:kfl_cpCoefLT_end_tab_index_chm) * wDiff
              if (kfl_cpCoefHT_tab_index_chm > 0) sphec_gp_ht(ielem) % a(igaus,1:6,1) = &
                  retva2(kfl_cpCoefHT_tab_index_chm:kfl_cpCoefHT_end_tab_index_chm) * wPrem + &
                  retva1(kfl_cpCoefHT_tab_index_chm:kfl_cpCoefHT_end_tab_index_chm) * wDiff 
              if (kfl_T_tab_index_chm > 0)        tempe_gp(ielem) % a(igaus,1,1)      = retva2(kfl_T_tab_index_chm) * wPrem &
                 + retva1(kfl_T_tab_index_chm) * wDiff
              if (kfl_DtRho_tab_index_chm > 0)    DtRho_gp(ielem) % a(igaus,1,1)      = retva2(kfl_DtRho_tab_index_chm) * wPrem &
                 + retva1(kfl_DtRho_tab_index_chm) * wDiff


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

        do idimt = 1, lookup_fw(kfl_tab_fw_chm_diff) % main_table % ndim
           if (lookup_fw(kfl_tab_fw_chm_diff) % kfl_chm_control(idimt) > 0) then
              !
              ! >0: one of the conces
              !
              control(idimt) = conce(ipoin,lookup_fw(kfl_tab_fw_chm_diff) % kfl_chm_control(idimt),1)
           else
              if (lookup_fw(kfl_tab_fw_chm_diff) % kfl_chm_control(idimt) == -1) then
                 !
                 ! -1: enthalpy
                 !
                 control(idimt) = therm(ipoin,1)
              elseif (lookup_fw(kfl_tab_fw_chm_diff) % kfl_chm_control(idimt) == -2) then
                 !
                 ! -2: scalar dissipation rate
                 !
                 control(idimt) = xZr_chm(ipoin) + xZs_chm(ipoin)
              endif
           endif
        enddo

        call fw_lookup(control, tab_scale_control1, lookup_fw(kfl_tab_fw_chm_diff), &
                       retva1(1:lookup_fw(kfl_tab_fw_chm_diff)%main_table%nvar),ind)
 
        if (kfl_cpCoefLT_tab_index_chm > 0) sphec(ipoin,1:6,1) = retva1(kfl_cpCoefLT_tab_index_chm:kfl_cpCoefLT_end_tab_index_chm)
        if (kfl_cpCoefHT_tab_index_chm > 0) sphec(ipoin,1:6,2) = retva1(kfl_cpCoefHT_tab_index_chm:kfl_cpCoefHT_end_tab_index_chm)
     endif
  enddo bou_nodes

end subroutine chm_gp_multi_reatab
