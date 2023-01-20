!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_outvar(ivari,imesh)
  !------------------------------------------------------------------------
  !****f* partis/chm_output
  ! NAME
  !    chm_output
  ! DESCRIPTION
  !    Output a postprocess variable
  ! USES
  !    postpr
  !    memgen
  ! USED BY
  !    chm_output
  !***
  !------------------------------------------------------------------------
  use def_parame,                         only : zero
  use def_master,                         only : INOTEMPTY, INOTMASTER, ITASK_BEGSTE, postp, div_enthalpy_transport, wmean_gp,&
                                                 tempe_gp, unkno, radiative_heat, enthalpy_transport, intost, mem_modul, modul,&
                                                 cutim, ger3p, gesca, gevec, ittim, mass_sink, tempe, conce, therm, wmean, massk
  use def_domain,                         only : mnode, ndime, mgaus, nelem, npoin, ltype, ngaus, elmar, nnode, lnods
  use def_kintyp,                         only : ip, rp, r1p, r2p
  use def_chemic,                         only : nclas_chm, ADR_chm, kfl_icmean_chm, mass_gp, mixedEq_eqs_chm, table_fw,&
                                                 posttable_fw, Yk_ssm_gp, DtRho_gp, aux_val_CMC_chm, avtim_chm, hrr_chm,&
                                                 kfl_DtRho_tab_index_chm, kfl_lookg_chm, kfl_max_srcfw_chm, kfl_model_chm,&
                                                 kfl_post_gp_CMC_chm, kfl_radia_chm, kfl_solve_cond_CMC_chm,&
                                                 kfl_solve_enth_CMC_chm, kfl_soot_chm, kfl_T_tab_index_chm,&
                                                 kfl_transfer_condField_CMC_chm, kfl_yk_fw_ssm, matr_scal_CMC_chm, nclas_chm,&
                                                 nspec_chm, nspec_cond_write_CMC_chm, nspec_transf_CMC_chm,&
                                                 nspec_uncond_write_CMC_chm, nZ_write_CMC_chm, posttable_fw, transf_entha_CMC_chm,&
                                                 transf_spec_CMC_chm, write_uncond_spec_CMC_chm, entha_chm, avY_chm, avYv_chm,&
                                                 avZv_chm, avZ_chm, avZ2_chm, rspec_chm, zgradmax_chm, phi_chm, xYr_chm, xZr_chm,&
                                                 xYs_chm, xZs_chm, avY2_chm, avL_chm, avL2_chm, avS_chm, avS0_chm, avd32_chm,&
                                                 avden_chm, Sigma_chm, Sigm0_chm, d32_chm,&
                                                 massConsumption_gp, hrr_avg_chm, avmsk_chm, avposttab_chm, write_iZ_spec_CMC_chm,&
                                                 write_cond_spec_CMC_chm, Zavg_CMC_chm, av_Z_flux_chm, av_Yk_CMC_chm,&
                                                 av_Yk_int_CMC_chm, av_enthalp_CMC_chm, av_enthalp_int_CMC_chm, av_temp_CMC_chm,&
                                                 av_temp_int_CMC_chm, av_visco_lam_int_CMC_chm, av_Named_Unk_chm, Zvar_CMC_chm,&
                                                 grad_Yk, React_ind, Yk_CMC_chm, Yk_int_CMC_chm, enthalp_CMC_chm,&
                                                 enthalp_int_CMC_chm, temp_CMC_chm, temp_int_CMC_chm, src_Yk_CMC_chm,&
                                                 t_chem_integ_CMC_chm, src_Yk_int_CMC_chm, flame_index_gp, zgrad_gp,&
                                                 kfl_multimod_chm,kfl_tab_fw_chm_diff
  use mod_postpr,                         only : postpr
  use mod_ker_proper,                     only : ker_proper
  use def_kermod,                         only : max_lookup_dim, lookup_fw, kfl_chemic_vect
  use mod_memory,                         only : memory_alloca, memory_deallo
  use mod_arrays,                         only : arrays
  use mod_solver,                         only : solver_lumped_mass_system
  use mod_chm_entropy,                    only : chm_entropy_postprocess
  use mod_interp_tab,                     only : fw_scale_cont_var
  use mod_interp_tab,                     only : fw_lookup
  use mod_chm_finiteRate,                 only : get_mixture_fraction
  use mod_chm_thermophoretic,             only : Vthermoph_nodes
  use mod_chm_sectional_soot_model_fast,  only : Qnucl_gp_chm
  use mod_chm_sectional_soot_model_fast,  only : Qcoag_gp_chm
  use mod_chm_sectional_soot_model_fast,  only : Qcond_gp_chm
  use mod_chm_sectional_soot_model_fast,  only : Qsurf_gp_chm
  use mod_chm_sectional_soot_model_fast,  only : Qtot_gp_chm
  use mod_chm_sectional_soot_model_fast,  only : RhoC_ssm
  use mod_chm_sectional_soot_model_fast,  only : nsect_ssm
  use mod_chm_sectional_soot_model_fast,  only : idsm_0_ssm
  use mod_chm_sectional_soot_model_fast,  only : Qsoot_gp_chm, NDsoot_gp_chm
  use mod_outvar,                         only : outvar
  use mod_chm_mixedEq,                    only : CHM_SRC_OFF,CHM_SRC_TABLE,CHM_SRC_DSM,CHM_SRC_ANN
  use mod_chm_operations_CMC,             only : chm_PDF_calc_CMC, chm_smooth_field_CMC
  use mod_arrays,                         only : arrays_name
  use mod_arrays,                         only : arrays_number
#ifdef CANTERA
  use def_chemic,                         only : src_chm, gas_chm
  use cantera,                            only : getSpeciesName
#endif
  implicit none
  integer(ip), intent(in) :: ivari
  integer(ip), intent(in) :: imesh   !< Mesh to postprocess on
  integer(ip)             :: ipoin,iclas,ipostvar,dummi,ielem,pelty,pgaus,imixf
  integer(ip)             :: idime
  integer(ip)             :: idsm
  integer(ip)             :: idimt
  real(rp)                :: rutim,auxi,retva(100_ip)
  character(5)            :: wopos(3)
  character(5)            :: wopos_aux
  character(2)            :: wclas
  type(r1p),pointer       :: aux_r1p(:)
  type(r2p),pointer       :: aux_r2p(:)
  real(rp), pointer       :: auxvar(:,:)

  real(rp)                  :: control(max_lookup_dim)   ! input of table lookup function
  real(rp)                  :: scale_control(max_lookup_dim)
  real(rp)                  :: lim_control(max_lookup_dim,2_ip)
  integer(ip), pointer      :: aux_name_ispec(:)
  integer(ip)               :: ind(max_lookup_dim), aux_nspec

  character(len=10)         :: valZ, val_spec

  !
  ! Variables to get unknown on Gaussian Integration Points
  !
  integer(ip) :: igaus, inode
  integer(ip) :: pnode
  integer(ip) :: lnods_loc(mnode)
  real(rp)    :: elcon(mnode,nclas_chm,ADR_chm(1) % ntime)
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: dummr(ndime,mnode)
  real(rp)    :: gpcon(mgaus,nclas_chm)

  external    :: memgen
  external    :: smooth
  external    :: smoot5
  external    :: runend
  external    :: chm_reatab
  external    :: chm_gp_reatab
  external    :: chm_gp_multi_reatab
  external    :: chm_post_scalar_dissipation_rate
  external    :: chm_post_scalar_dist
  external    :: chm_post_gather
  external    :: chm_post_gp_lookup
  external    :: chm_gp_ANN_eval
  external    :: cantera_elemh
  external    :: cantera_elemo
  external    :: cantera_elemc
  external    :: chm_post_flame_index

  if( ivari == 0 ) return
  !
  ! Define postprocess variable
  !
  rutim=cutim
  nullify(gesca)
  nullify(gevec)
  nullify(ger3p)

  select case ( arrays_name(ivari) )  


  case( 'CONCE' )
     !
     ! Class concentration
     !
     call arrays(ivari,'POSTPROCESS',conce)
     return

  case( 'SOURC' )
     !
     ! Source term SOURC
     !

     !
     ! Update based on current unknown
     !
     if ((kfl_model_chm == 1 .or. kfl_model_chm == 2) .and. kfl_max_srcfw_chm > 0) then
        !
        ! Taulated sources
        !
        if (kfl_lookg_chm > 0) then
            if (kfl_multimod_chm > 0_ip) then
               call chm_gp_multi_reatab(ITASK_BEGSTE)
               call chm_gp_multi_reatab(ITASK_BEGSTE)
            else 
               call chm_gp_reatab(ITASK_BEGSTE)
               call chm_gp_reatab(ITASK_BEGSTE)
            end if 
            call chm_gp_ANN_eval(ITASK_BEGSTE)
        else
            call chm_reatab()
        endif
     endif

     if (kfl_model_chm == 1 ) then
        !
        ! Flamelet model, source only for class 1
        !
        if (kfl_lookg_chm > 0) then
           nullify(aux_r1p)
           call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p,max(1_ip,nelem))
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              call memory_alloca(mem_modul(1:2,modul),'AUX_R1P % A','chm_outvar',aux_r1p(ielem)%a,pgaus)
              aux_r1p(ielem) % a = mass_gp(ielem) % a(:,kfl_icmean_chm,1)
           end do
           call memgen(zero,max(1_ip,npoin),zero)
           call smooth (aux_r1p, gesca)
           call memory_deallo(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p)
        else
            if(INOTMASTER ) gesca => massk(:,kfl_icmean_chm)
        endif

     else if ( kfl_model_chm == 2  ) then
        !
        ! Mixed equation model, output sources that are non-passive
        !
        wopos(2)=postp(1)%wopos(2,arrays_number('SOURC'))
        wopos(3)=postp(1)%wopos(3,arrays_number('SOURC'))

        if (kfl_lookg_chm > 0) then
           if( INOTEMPTY ) then
              nullify(auxvar)
              call memory_alloca(mem_modul(1:2,modul),'AUXVAR', 'chm_outvar',auxvar, nclas_chm, npoin)
              nullify(aux_r2p)
              call memory_alloca(mem_modul(1:2,modul),'AUX_R2P','chm_outvar',aux_r2p,nelem)
              do ielem = 1,nelem
                 pelty = ltype(ielem)
                 pgaus = ngaus(pelty)
                 call memory_alloca(mem_modul(1:2,modul),'AUX_R2P % A','chm_outvar',aux_r2p(ielem)%a,pgaus,nclas_chm)
              end do
              !
              ! Get sources
              !
              do iclas = 1,nclas_chm
                 if ( mixedEq_eqs_chm (iclas) % kfl_source_type ==  CHM_SRC_TABLE .or. &
                    & mixedEq_eqs_chm (iclas) % kfl_source_type ==  CHM_SRC_ANN) then
                    !
                    ! Tabulated source terms
                    !
                    do ielem = 1,nelem
                       pelty = ltype(ielem)
                       pgaus = ngaus(pelty)
                       aux_r2p(ielem) % a (1:pgaus,iclas) = mass_gp(ielem) % a(1:pgaus,iclas,1)
                    enddo
                 endif
              enddo

              !
              ! Add contributon of split source term
              !
              do ielem = 1,nelem
                 !
                 ! Get elemental quantities
                 !
                 pelty = ltype(ielem)
                 pnode = nnode(pelty)
                 pgaus = ngaus(pelty)
                 lnods_loc(1:pnode) = lnods(1:pnode,ielem)
                 call chm_post_gather(&
                      pnode,lnods_loc,elcon(1:pnode,:,:),elcod,dummr)

                 !
                 ! Mass fractions on Gaussian integration points
                 !
                 gpcon = 0.0_rp
                 do iclas = 1,nclas_chm
                    do igaus = 1,pgaus
                       do inode = 1,pnode
                          gpcon(igaus,iclas) = gpcon(igaus,iclas)&
                               + elmar(pelty)%shape(inode,igaus) * elcon(inode,iclas,1)
                       end do
                    end do
                 end do

                 do iclas = 1,nclas_chm
                    if ( mixedEq_eqs_chm (iclas) % kfl_source_type ==  CHM_SRC_TABLE ) then
                       if ( mixedEq_eqs_chm(iclas) % kfl_source_split == 1) then
                           !
                           ! Add contribution of first order dependence (consumption)
                           !
                           aux_r2p(ielem) % a (1:pgaus,iclas) = aux_r2p(ielem) % a (1:pgaus,iclas) + gpcon(1:pgaus,iclas)&
                               * massConsumption_gp(ielem) % a(1:pgaus,iclas,1)
                       endif
                    endif
                 enddo
              enddo

              if (kfl_soot_chm /= 0) then
                 !
                 ! Source terms from discrete sectional method
                 !
                 do idsm = 1,nsect_ssm
                    iclas = idsm + idsm_0_ssm

                    if ( mixedEq_eqs_chm (iclas) % kfl_source_type == CHM_SRC_DSM) then
                       do ielem = 1,nelem
                          pelty = ltype(ielem)
                          pgaus = ngaus(pelty)
                          aux_r2p(ielem) % a (1:pgaus,iclas) = Qtot_gp_chm(ielem) % a (1:pgaus,idsm,1)
                       enddo
                    end if
                 end do
              end if
              call smoot5(aux_r2p,auxvar,nclas_chm)
              call memory_deallo(mem_modul(1:2,modul),'AUX_R2P','chm_outvar',aux_r2p)
           endif

        endif

        do iclas = 1,nclas_chm
           if ( mixedEq_eqs_chm (iclas) % kfl_source_type /=  CHM_SRC_OFF ) then
              if( INOTEMPTY ) then
                 if (kfl_lookg_chm > 0) then
                    gesca => auxvar(iclas,:)
                 endif
              endif

              !
              ! Get name
              !
              wclas =  trim(intost(iclas))
              if(iclas<10) then
                 wopos(1)='SRC'//'0'//wclas(1:1)
              else
                 wopos(1)='SRC'//wclas
              end if

              !
              ! Write postprocessed data
              !
              call postpr(gesca,wopos,ittim,rutim)
           endif
        end do
        if( INOTMASTER ) nullify(gesca)

        if (kfl_lookg_chm > 0) then
           if( INOTEMPTY ) then
              call memory_deallo(mem_modul(1:2,modul),'AUXVAR','chm_outvar',auxvar)
           endif
        endif

        return

     else if (kfl_model_chm == 3 ) then
        !
        ! Finite Rate Chemistry - Chemical Sources
        !
        wopos(2)=postp(1)%wopos(2,arrays_number('SOURC'))
        wopos(3)=postp(1)%wopos(3,arrays_number('SOURC'))

#ifdef CANTERA
        do iclas=1,nclas_chm
           if( INOTMASTER ) gesca => src_chm(1:npoin,iclas)
           if (iclas < nspec_chm ) then
              call getSpeciesName(gas_chm,iclas,wopos(1))
              wopos(1) = 'S_'//wopos(1)
           else
              wopos(1) = 'Q_'
           end if
           call postpr(gesca,wopos,ittim,rutim)
        end do
#endif
        if( INOTMASTER ) nullify(gesca)
        return
     end if

  case( 'VISCO' )
    !
    ! Viscosity
    !
    if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       call ker_proper('VISCO','NPOIN',dummi,dummi,gesca)
    end if


  case( 'SPHEA' )
    !
    ! Specific heat
    !
    if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       call ker_proper('SPHEA','NPOIN',dummi,dummi,gesca)
    end if

  case( 'CONDU' )
    !
    ! Heat conductivity
     if(INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call ker_proper('CONDU','NPOIN',dummi,dummi,gesca)
     end if

  case( 'ENTHA' )
     !
     ! Enthalpy
     !
     if(INOTMASTER ) then
        call memgen(zero,npoin,zero)
        gesca=0.0_rp
        do iclas=1,nspec_chm
           do ipoin=1,npoin
              gesca(ipoin) = gesca(ipoin) + entha_chm(ipoin,iclas) * conce(ipoin,iclas,1)
           enddo
        enddo
     end if


  case( 'DIVEN' )
     !
     ! Enthalpy heat source term
     !
     if(INOTMASTER ) then
        nullify(aux_r1p)
        call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'AUX_R1P % A','chm_outvar',aux_r1p(ielem)%a,pgaus)
           aux_r1p(ielem) % a = div_enthalpy_transport(ielem) % a(:,1,1)
        end do
        call memgen(zero,npoin,zero)
        call smooth (aux_r1p, gesca)
        call memory_deallo(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p)
     end if

  case( 'SUMCO' )
     !
     ! Sum of concentration
     !
     if(INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin=1,npoin
           gesca(ipoin) = 0.0_rp
           do iclas=1,nspec_chm
              gesca(ipoin) = gesca(ipoin) + conce(ipoin,iclas,1)
           end do
        end do
     end if

  case( 'MOLEC' )
     !
     ! Molecular weight
     !
     if(INOTMASTER ) then
        if (kfl_model_chm /= 4) then
           if (kfl_lookg_chm > 0) then
               nullify(aux_r1p)
               call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p,nelem)
               do ielem = 1,nelem
                  pelty = ltype(ielem)
                  pgaus = ngaus(pelty)
                  call memory_alloca(mem_modul(1:2,modul),'AUX_R1P % A','chm_outvar',aux_r1p(ielem)%a,pgaus)
                  aux_r1p(ielem) % a = wmean_gp(ielem) % a(:,1,1)
               end do
               call memgen(zero,npoin,zero)
               call smooth (aux_r1p, gesca)
               call memory_deallo(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p)
           else
               gesca => wmean(:,1)
           endif

        else
           if (kfl_solve_cond_CMC_chm == 1) then
              gesca => wmean(:,1)
           end if
        end if
     end if

  case( 'TEMPE' )
    !
    ! Temperature
    !
    if(INOTMASTER ) then
        if(associated(tempe)) then
            gesca => tempe(:,1)
        elseif(associated(tempe_gp) .and. kfl_lookg_chm > 0) then
           nullify(aux_r1p)
           call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p,nelem)
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              call memory_alloca(mem_modul(1:2,modul),'AUX_R1P % A','chm_outvar',aux_r1p(ielem)%a,pgaus)
              aux_r1p(ielem) % a = tempe_gp(ielem) % a(:,1,1)
           end do
           call memgen(zero,npoin,zero)
           call smooth (aux_r1p, gesca)
           call memory_deallo(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p)
        end if
    end if

  case( 'AVY  ' )
    !
    ! AV Yc or C
    !
    if (cutim > avtim_chm) then
       auxi = cutim - avtim_chm
       if ( INOTMASTER ) then
          do ipoin=1,npoin
             unkno(ipoin)     = avY_chm(ipoin) / auxi
             avY_chm(ipoin) = 0.0_rp
          end do
          gesca => unkno
       endif
    else
       return
    endif

  case( 'AVYV ' )
    !
    ! AVYv
    !
    if (cutim > avtim_chm) then
       auxi = cutim - avtim_chm
       if ( INOTMASTER ) then
          do ipoin=1,npoin
             unkno(ipoin)     = avYv_chm(ipoin) / auxi
             avYv_chm(ipoin) = 0.0_rp
          end do
          gesca => unkno
       endif
    else
       return
    endif

  case( 'AVZV ' )
    !
    ! AV variance of Z
    !
    if (cutim > avtim_chm) then
       auxi = cutim - avtim_chm
       if ( INOTMASTER ) then
          do ipoin=1,npoin
             unkno(ipoin)     = avZv_chm(ipoin) / auxi
             avZv_chm(ipoin) = 0.0_rp
          end do
          gesca => unkno
       endif
    else
       return
    endif

  case( 'AVZ  ' )
    !
    ! AVZ
    !
    if (cutim > avtim_chm) then
       auxi = cutim - avtim_chm
       if ( INOTMASTER ) then
          do ipoin=1,npoin
             unkno(ipoin)     = avZ_chm(ipoin) / auxi
             avZ_chm(ipoin) = 0.0_rp
          end do
          gesca => unkno
       endif
    else
       return
    endif

  case( 'AVZ2 ' )
    !
    ! AVZ2
    !
    if (cutim > avtim_chm) then
       auxi = cutim - avtim_chm
       if ( INOTMASTER ) then
          do ipoin=1,npoin
             unkno(ipoin)     = avZ2_chm(ipoin) / auxi
             avZ2_chm(ipoin) = 0.0_rp
          end do
          gesca => unkno
       endif
    else
       return
    endif

  case( 'SPEC1' )
    !
    !  Species post-processing for Flamelet, radiation
    !
    if (kfl_radia_chm > 0) then
      if(INOTMASTER ) then
         call memgen(zero,npoin,zero)
         do ipoin=1,npoin
            gesca(ipoin) = rspec_chm(1,ipoin)
         end do
      end if
    end if

  case( 'SPEC2' )
    !
    !  Species post-processing for Flamelet, radiation
    !
    if (kfl_radia_chm > 0) then
      if(INOTMASTER ) then
         call memgen(zero,npoin,zero)
         do ipoin=1,npoin
            gesca(ipoin) = rspec_chm(2,ipoin)
         end do
      end if
    end if

  case( 'RADIA' )
     !
     ! Radiative source term for heat equation, Flamelet model
     !
    if (kfl_radia_chm > 0) then
     if(INOTMASTER ) then
        nullify(aux_r1p)
        call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'AUX_R1P % A','chm_outvar',aux_r1p(ielem)%a,pgaus)
           aux_r1p(ielem) % a = radiative_heat(ielem) % a(:,1,1)
        end do
        call memgen(zero,npoin,zero)
        call smooth (aux_r1p, gesca)
        call memory_deallo(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p)
     end if
    end if

  case( 'ZGRMA' )
     !
     ! Maximum gradient of the mixture fraction from diffusion flamelets
     !
     if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       do ipoin=1,npoin
          gesca(ipoin) = zgradmax_chm(ipoin)
       end do
     end if
  case( 'PHI  ' )
     !
     ! Weighting factor between premixed and diffusion flamelets
     !
     if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       do ipoin=1,npoin
          gesca(ipoin) = phi_chm(ipoin)
       end do
     end if

  case( 'XYR  ' )
     !
     ! Scalar dissipation rate of the reaction progress variable (resolved)
     !
     if(kfl_chemic_vect /= 1_ip) then
        call chm_post_scalar_dissipation_rate(23_ip)
     else
        call chm_post_scalar_dist(23_ip)
     end if
     if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       do ipoin=1,npoin
          gesca(ipoin) = xYr_chm(ipoin)
       end do
     end if

  case( 'XZR  ' )
     !
     ! Scalar dissipation rate of the mixture fraction (resolved)
     !
     if(kfl_chemic_vect /= 1_ip) then
        call chm_post_scalar_dissipation_rate(24_ip)
     else
        call chm_post_scalar_dist(24_ip)
     end if
     if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       do ipoin=1,npoin
          gesca(ipoin) = xZr_chm(ipoin)
       end do
     end if

  case( 'XYS  ' )
     !
     ! Scalar dissipation rate of the reaction progress variable (subgrid)
     !
     if(kfl_chemic_vect /= 1_ip) then
        call chm_post_scalar_dissipation_rate(25_ip)
     else
        call chm_post_scalar_dist(25_ip)
     end if
     if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       do ipoin=1,npoin
          gesca(ipoin) = xYs_chm(ipoin)
       end do
     end if

  case( 'XZS  ' )
     !
     ! Scalar dissipation rate of the mixture fraction (subgrid)
     !
     if(kfl_chemic_vect /= 1_ip) then
        call chm_post_scalar_dissipation_rate(26_ip)
     else
        call chm_post_scalar_dist(26_ip)
     end if
     if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       do ipoin=1,npoin
          gesca(ipoin) = xZs_chm(ipoin)
       end do
     end if

  case('AVXYR','AVXZR','AVXYS', 'AVXZS')
      call runend('chm_outvar: average scalar dissipation rates are not coded')

  case( 'AVY2 ' )
     !
     ! Average progress variable squared
     !
     if (cutim > avtim_chm) then
        auxi = cutim - avtim_chm
        if ( INOTMASTER ) then
           do ipoin=1,npoin
              unkno(ipoin)     = avY2_chm(ipoin) / auxi
              avY2_chm(ipoin) = 0.0_rp
           end do
           gesca => unkno
        endif
     else
        return
     endif

  case( 'AVL  ' )
     !
     ! Average liquid volume fraction
     !
     if (cutim > avtim_chm) then
        auxi = cutim - avtim_chm

        if ( INOTMASTER ) then
           do ipoin=1,npoin
              unkno(ipoin)     = avL_chm(ipoin) / auxi
              avL_chm(ipoin) = 0.0_rp
           end do
           gesca => unkno
        endif
     else
        return
     endif

  case( 'AVL2 ' )
     !
     ! Average liquid volume fraction squared
     !
     if (cutim > avtim_chm) then
        auxi = cutim - avtim_chm
        if ( INOTMASTER ) then
           do ipoin=1,npoin
              unkno(ipoin)     = avL2_chm(ipoin) / auxi
              avL2_chm(ipoin) = 0.0_rp
           end do
           gesca => unkno
        endif
     else
        return
     endif

  case( 'AVS  ' )
     !
     ! Average interface surface density Sigma
     !
     if (cutim > avtim_chm) then
        auxi = cutim - avtim_chm
        if ( INOTMASTER ) then
           do ipoin=1,npoin
              unkno(ipoin)     = avS_chm(ipoin) / auxi
              avS_chm(ipoin) = 0.0_rp
           end do
           gesca => unkno
        endif
     else
        return
     endif

  case( 'AVS0 ' )
     !
     ! Average interface surface density Sigma_0
     !
     if (cutim > avtim_chm) then
        auxi = cutim - avtim_chm
        if ( INOTMASTER ) then
           do ipoin=1,npoin
              unkno(ipoin)     = avS0_chm(ipoin) / auxi
              avS0_chm(ipoin) = 0.0_rp
           end do
           gesca => unkno
        endif
     else
        return
     endif

  case( 'AVD32' )
     !
     ! Average Sauter mean diameter
     !
     if (cutim > avtim_chm) then
        auxi = cutim - avtim_chm
        if ( INOTMASTER ) then
           do ipoin=1,npoin
              unkno(ipoin)     = avd32_chm(ipoin) / auxi
              avd32_chm(ipoin) = 0.0_rp
           end do
           gesca => unkno
        endif
     else
        return
     endif

  case( 'AVDEN' )
     !
     ! Average density
     !
     if (cutim > avtim_chm) then
        auxi = cutim - avtim_chm
        if ( INOTMASTER ) then
           do ipoin=1,npoin
              unkno(ipoin)     = avden_chm(ipoin) / auxi
              avden_chm(ipoin) = 0.0_rp
           end do
           gesca => unkno
        endif
     else
        return
     endif

  case( 'SIGMA' )
     !
     ! Interface surface density Sigma
     !
     if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       do ipoin=1,npoin
          gesca(ipoin) = Sigma_chm(ipoin)
       end do
     end if

  case( 'SIGM0' )
     !
     ! Interface surface density Sigma_0 or Sigma_min
     !
     if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       do ipoin=1,npoin
          gesca(ipoin) = Sigm0_chm(ipoin)
       end do
     end if

  case( 'D32  ' )
     !
     ! Sauter mean diameter
     !
     if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       do ipoin=1,npoin
          gesca(ipoin) = d32_chm(ipoin)
       end do
     end if

  case( 'GRADY' )
     !
     ! Gradient mass fractions
     !
     if( INOTMASTER )then
        call memgen(zero,ndime,npoin)
        do ipoin=1,npoin
           gevec(1:ndime,ipoin) = grad_Yk(1,1:ndime,ipoin)
        end do
     end if

  case( 'HTRAN' )
     !
     ! Enthalpy transport by diffusion
     !
     if (kfl_model_chm /= 4) then
        if( INOTMASTER )then
           nullify(aux_r2p)
           call memory_alloca(mem_modul(1:2,modul),'AUX_R2P','chm_outvar',aux_r2p,nelem)
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              call memory_alloca(mem_modul(1:2,modul),'AUX_R2P % A','chm_outvar',aux_r2p(ielem)%a,pgaus,ndime)
              aux_r2p(ielem) % a = enthalpy_transport(ielem) % a(:,1:ndime,1)
           end do
           call memgen(zero,ndime,npoin)
           call smoot5(aux_r2p,gevec,ndime)
           call memory_deallo(mem_modul(1:2,modul),'AUX_R2P','chm_outvar',aux_r2p)
        end if
     end if

  case( 'ELEMH' )
       !
       ! Elemental Mass Fraction - H
       !
       if( INOTMASTER )then
          call memgen(zero,npoin,zero)
#ifdef CANTERA
          do ipoin=1,npoin
             call cantera_elemh(conce(ipoin,:,1),gesca(ipoin))
          end do
#endif
       end if

  case( 'ELEMO' )
       !
       ! Elemental Mass Fraction - O
       !
       if( INOTMASTER )then
          call memgen(zero,npoin,zero)
#ifdef CANTERA
          do ipoin=1,npoin
             call cantera_elemo(conce(ipoin,:,1),gesca(ipoin))
          end do
#endif
       end if

  case( 'ELEMC' )
       !
       ! Elemental Mass Fraction - C
       !
       if( INOTMASTER )then
          call memgen(zero,npoin,zero)
#ifdef CANTERA
          do ipoin=1,npoin
             call cantera_elemc(conce(ipoin,:,1),gesca(ipoin))
          end do
#endif
       end if


  case( 'MASSK' )
     !
     ! Mass source from partis
     !
     call memgen(zero,max(1_ip,npoin),zero)
     if (associated(mass_sink)) then
        do ipoin=1,npoin
           gesca(ipoin) = mass_sink(ipoin)
        enddo
        call solver_lumped_mass_system(1_ip,gesca,EXCHANGE=.false.)
     else
        do ipoin=1,npoin
           gesca(ipoin) = 0.0_rp
        end do
     endif

  case( 'ENVIC' )
     !
     ! Entropy viscosity maximum among species
     !
     call memgen(zero,max(npoin,1_ip),zero)
     call chm_entropy_postprocess(0_ip,gesca)

  case( 'SCONC' )
     !
     ! Scaled progress variable
     !
     wopos(2)=postp(1)%wopos(2,arrays_number('SCONC'))
     wopos(3)=postp(1)%wopos(3,arrays_number('SCONC'))

     if(kfl_multimod_chm == 1_ip) then
        if( INOTMASTER ) then
           nullify(auxvar)
           call memory_alloca(mem_modul(1:2,modul),'AUXVAR', 'chm_outvar',auxvar, npoin,&
              lookup_fw(kfl_tab_fw_chm_diff) % main_table % ndim)
        endif

        !
        ! Do lookup
        !
        ind = 1_ip
        do ipoin=1,npoin
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

            call fw_scale_cont_var( control, scale_control, lim_control, lookup_fw(kfl_tab_fw_chm_diff), ind)
            do idimt=1,lookup_fw(kfl_tab_fw_chm_diff) % main_table % ndim
               auxvar(ipoin,idimt) = scale_control(idimt)
            enddo

        end do

        do idimt=1,lookup_fw(kfl_tab_fw_chm_diff) % main_table % ndim
           if( INOTMASTER ) then
               gesca => auxvar(:,idimt)
           endif

           wclas =  trim(intost(idimt))
           if(idimt<10) then
              wopos(1)=postp(1)%wopos(1,arrays_number('SCONC'))(1:3)//'0'//wclas(1:1)
           else
              wopos_aux=postp(1)%wopos(1,arrays_number('SCONC'))(1:3) ! Asign in two steps to avoid bug in ibm xlf
              wopos(1)=trim(wopos_aux)//wclas
           end if
           call postpr(gesca,wopos,ittim,rutim)
        end do
        if( INOTMASTER ) nullify(gesca)

        if( INOTMASTER ) then
           call memory_deallo(mem_modul(1:2,modul),'AUXVAR','chm_outvar',auxvar)
        endif
     else
        if( INOTMASTER ) then
           nullify(auxvar)
           call memory_alloca(mem_modul(1:2,modul),'AUXVAR', 'chm_outvar',auxvar, npoin, table_fw % main_table % ndim)
        endif

        !
        ! Do lookup
        !
        ind = 1_ip
        do ipoin=1,npoin
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

            call fw_scale_cont_var( control, scale_control, lim_control, table_fw, ind)
            do idimt=1,table_fw % main_table % ndim
               auxvar(ipoin,idimt) = scale_control(idimt)
            enddo

        end do

        do idimt=1,table_fw % main_table % ndim
           if( INOTMASTER ) then
               gesca => auxvar(:,idimt)
           endif

           wclas =  trim(intost(idimt))
           if(idimt<10) then
              wopos(1)=postp(1)%wopos(1,arrays_number('SCONC'))(1:3)//'0'//wclas(1:1)
           else
              wopos_aux=postp(1)%wopos(1,arrays_number('SCONC'))(1:3) ! Asign in two steps to avoid bug in ibm xlf
              wopos(1)=trim(wopos_aux)//wclas
           end if
           call postpr(gesca,wopos,ittim,rutim)
        end do
        if( INOTMASTER ) nullify(gesca)

        if( INOTMASTER ) then
           call memory_deallo(mem_modul(1:2,modul),'AUXVAR','chm_outvar',auxvar)
        endif
     end if 
     return

  case( 'NREAC' )
     !
     ! Sum of Reactions (Finite Rate)
     !
     if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       do ipoin=1,npoin
          gesca(ipoin) = real(sum(React_ind(ipoin,:)),rp)
       end do
     end if

  case( 'MIXFR' )
     !
     ! Mixture Fraction
     !
     if(INOTMASTER ) then
       call memgen(zero,npoin,zero)
       do ipoin=1,npoin
          call get_mixture_fraction(conce(ipoin,:,1),gesca(ipoin))
       end do
     end if

  case( 'HRR  ' )
        !
        ! Finite Rate Chemistry + CMC - Instantaneous heat release
        !
        if( INOTMASTER ) gesca => hrr_chm

  case( 'AVHRR' )
        !
        ! Finite Rate Chemistry + CMC - Averaged heat release
        !
        if (cutim > avtim_chm) then
           auxi = cutim - avtim_chm
           if ( INOTMASTER ) then
              do ipoin=1,npoin
                 unkno(ipoin)           = hrr_avg_chm(ipoin) / auxi
                 hrr_avg_chm(ipoin) = 0.0_rp
              end do
              gesca => unkno
           endif
        else
           return
        endif

  case( 'AVMSK' )
     !
     ! Average mass source from spray
     !
     if (cutim > avtim_chm) then
        auxi = cutim - avtim_chm
        do ipoin=1,npoin
           unkno(ipoin)     = avmsk_chm(ipoin) / auxi
           avmsk_chm(ipoin) = 0.0_rp
        end do
        gesca => unkno
     else
        return
     endif

  case( 'POSTT' )
     !
     ! Post-processing lookup variables
     !
     wopos(2)=postp(1)%wopos(2,arrays_number('POSTT'))
     wopos(3)=postp(1)%wopos(3,arrays_number('POSTT'))

     if( INOTEMPTY ) then
        nullify(auxvar)
        call memory_alloca(mem_modul(1:2,modul),'AUXVAR', 'chm_outvar',auxvar,  posttable_fw % main_table % nvar, npoin)
        nullify(aux_r2p)
        call memory_alloca(mem_modul(1:2,modul),'AUX_R2P','chm_outvar',aux_r2p,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'AUX_R2P % A','chm_outvar',aux_r2p(ielem)%a,pgaus,&
               posttable_fw % main_table % nvar)
        end do
        !
        ! Lookup from postprocessing table
        !
        call chm_post_gp_lookup(aux_r2p,posttable_fw)
        call smoot5(aux_r2p,auxvar,posttable_fw % main_table % nvar)
        call memory_deallo(mem_modul(1:2,modul),'AUX_R2P','chm_outvar',aux_r2p)
     endif


     do ipostvar=1,posttable_fw % main_table % nvar
        if( INOTEMPTY ) then
            gesca => auxvar(ipostvar,:)
        endif
        wopos(1)=posttable_fw % main_table % varname(ipostvar)
        call postpr(gesca,wopos,ittim,rutim)
     end do
     if( INOTMASTER ) nullify(gesca)

     if( INOTEMPTY ) then
        call memory_deallo(mem_modul(1:2,modul),'AUXVAR','chm_outvar',auxvar)
     endif

     return

  case( 'AVPOT' )
    !
    ! AVPOT:  Average posttab properties
    !
    wopos(2)=postp(1)%wopos(2,arrays_number('AVPOT'))
    wopos(3)=postp(1)%wopos(3,arrays_number('AVPOT'))
    if (cutim > avtim_chm) then
       auxi = cutim - avtim_chm

       do ipostvar=1,posttable_fw % main_table % nvar
          do ipoin=1,npoin
             unkno(ipoin)                  = avposttab_chm(ipoin,ipostvar) / auxi
             avposttab_chm(ipoin,ipostvar) = 0.0_rp
          end do
          if( INOTEMPTY ) then
              gesca => unkno
          endif
          if (posttable_fw % main_table % varname(ipostvar)(1:2) == 'YK') then
             wopos(1)='AV'//trim(posttable_fw % main_table % varname(ipostvar)(3:5))
          else
             wopos(1)='AV'//trim(posttable_fw % main_table % varname(ipostvar)(1:3))
          endif
          call postpr(gesca,wopos,ittim,rutim)
       end do
       if( INOTMASTER ) nullify(gesca)
    endif

    return

  case( 'MASCN' )
      !
      ! Conditional mass fractions for CMC model
      !
      if ( kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1 ) then
         wopos(2) = 'SCALA'  !! Yk_CMC_chm is MATRIX so wopos(2) has to be changed
         wopos(3) = postp(1)%wopos(3,arrays_number('MASCN'))

         do imixf = 1, nZ_write_CMC_chm
            do iclas = 1, nspec_cond_write_CMC_chm
               if ( INOTEMPTY ) then
                  if (kfl_post_gp_CMC_chm == 1_ip) then
                     nullify(matr_scal_CMC_chm)
                     matr_scal_CMC_chm => Yk_CMC_chm(write_iZ_spec_CMC_chm(imixf),1:npoin,write_cond_spec_CMC_chm(iclas))
                     call chm_smooth_field_CMC
                     gesca => aux_val_CMC_chm
                  else
                     gesca => Yk_CMC_chm(write_iZ_spec_CMC_chm(imixf),1:npoin,write_cond_spec_CMC_chm(iclas))
                  end if
               end if
               write(valZ,'(I2.2)') write_iZ_spec_CMC_chm(imixf)
               if (write_cond_spec_CMC_chm(iclas) < 100) then
                  write(val_spec,'(I2.2)') write_cond_spec_CMC_chm(iclas)
                  wopos(1) = trim(val_spec) // '_' // valZ
               else
                  write(val_spec,'(I3.3)') write_cond_spec_CMC_chm(iclas)
                  wopos(1) = trim(val_spec) // valZ
               end if
               call postpr(gesca,wopos,ittim,rutim)
            end do
         end do
         if( INOTEMPTY ) nullify(gesca)
      end if


  case( 'MASUN' )
      !
      ! Unconditional mass fractions for CMC model
      !
      if ( kfl_model_chm == 4 ) then
         wopos(2)=postp(1)%wopos(2,arrays_number('MASUN'))
         wopos(3)=postp(1)%wopos(3,arrays_number('MASUN'))

         nullify(aux_name_ispec)
         if (kfl_solve_cond_CMC_chm == 1 ) then
            aux_nspec = nspec_uncond_write_CMC_chm
            aux_name_ispec => write_uncond_spec_CMC_chm
         else
            if (kfl_transfer_condField_CMC_chm == 1_ip .and. nspec_transf_CMC_chm > 0_ip) then
               aux_nspec = nspec_transf_CMC_chm
               aux_name_ispec => transf_spec_CMC_chm
            else
               aux_nspec = 0_ip
            end if
         end if

         do iclas = 1, aux_nspec
            if( INOTMASTER ) gesca => Yk_int_CMC_chm(1:npoin,write_uncond_spec_CMC_chm(iclas))
            if (aux_name_ispec(iclas) < 100) then
               write(val_spec,'(I2.2)') aux_name_ispec(iclas)
            else
               write(val_spec,'(I3.3)') aux_name_ispec(iclas)
            end if
            wopos(1) = 'Y' // trim(val_spec)
            call postpr(gesca,wopos,ittim,rutim)
         end do
         if( INOTMASTER ) nullify(gesca)
      end if


  case( 'ENTCN' )
      !
      ! Conditional enthalpy for CMC model
      !
      if ( kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1 ) then
         wopos(2) = 'SCALA'  !! enthalp_CMC_chm is MATRIX so wopos(2) has to be changed
         wopos(3)=postp(1)%wopos(3,arrays_number('ENTCN'))

         if (kfl_solve_enth_CMC_chm /= 0) then
            do imixf = 1, nZ_write_CMC_chm
               if( INOTEMPTY ) then
                  if (kfl_post_gp_CMC_chm == 1_ip) then
                     nullify(matr_scal_CMC_chm)
                     matr_scal_CMC_chm => enthalp_CMC_chm(write_iZ_spec_CMC_chm(imixf),1:npoin)
                     call chm_smooth_field_CMC
                     gesca => aux_val_CMC_chm
                  else
                     gesca => enthalp_CMC_chm(write_iZ_spec_CMC_chm(imixf),1:npoin)
                  end if
               end if
               write(valZ,'(I2.2)') write_iZ_spec_CMC_chm(imixf)
               wopos(1) = 'h_' // trim(valZ)
               call postpr(gesca,wopos,ittim,rutim)
            end do
            if( INOTEMPTY )   nullify(gesca)
         end if
      end if


  case( 'ENTUN' )
      !
      ! Unconditional enthalpy for CMC model
      !
      if ( kfl_model_chm == 4 ) then
         if (kfl_solve_cond_CMC_chm == 1 .or. &
             (kfl_solve_cond_CMC_chm == 0 .and. kfl_transfer_condField_CMC_chm == 1_ip .and. &
             transf_entha_CMC_chm > 0_ip)) then
            if( INOTMASTER ) gesca => enthalp_int_CMC_chm(1:npoin)
         end if
      end if

  case( 'TEMCN' )
      !
      ! Conditional temperature for CMC model
      !
      if ( kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1 ) then
         wopos(2) = 'SCALA'  !! temp_CMC_chm is MATRIX so wopos(2) has to be changed
         wopos(3)=postp(1)%wopos(3,arrays_number('TEMCN'))

         do imixf = 1, nZ_write_CMC_chm
            if ( INOTMASTER ) then
               if (kfl_post_gp_CMC_chm == 1_ip) then
                  nullify(matr_scal_CMC_chm)
                  matr_scal_CMC_chm => temp_CMC_chm(write_iZ_spec_CMC_chm(imixf),1:npoin)
                  call chm_smooth_field_CMC
                  gesca => aux_val_CMC_chm
               else
                  gesca => temp_CMC_chm(write_iZ_spec_CMC_chm(imixf),1:npoin)
               end if
            end if
            write(valZ,'(I2.2)') write_iZ_spec_CMC_chm(imixf)
            wopos(1) = 'T_' // trim(valZ)
            call postpr(gesca,wopos,ittim,rutim)
         end do
         if( INOTMASTER ) nullify(gesca)
      end if

  case( 'TEMUN' )
      !
      ! Unconditional temperature for CMC model
      !
      if ( kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1 ) then
         if( INOTMASTER ) gesca => temp_int_CMC_chm(1:npoin)
      end if


  case( 'SRCCN' )
      !
      ! Conditional mass fractions source terms for CMC model
      !
      if ( kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1 ) then
         wopos(2)=postp(1)%wopos(2,arrays_number('SRCCN'))
         wopos(3)=postp(1)%wopos(3,arrays_number('SRCCN'))

         do imixf = 1, nZ_write_CMC_chm
            do iclas = 1, nspec_cond_write_CMC_chm
               if( INOTMASTER ) then
                  if (kfl_post_gp_CMC_chm == 1_ip) then
                     nullify(matr_scal_CMC_chm)
                     matr_scal_CMC_chm => src_Yk_CMC_chm(write_iZ_spec_CMC_chm(imixf),1:npoin,write_cond_spec_CMC_chm(iclas))
                     call chm_smooth_field_CMC
                     gesca => aux_val_CMC_chm
                  else
                     gesca => src_Yk_CMC_chm(write_iZ_spec_CMC_chm(imixf),1:npoin,write_cond_spec_CMC_chm(iclas))
                  end if
               end if
               write(valZ,'(I2.2)') write_iZ_spec_CMC_chm(imixf)
               write(val_spec,'(I2.2)') write_cond_spec_CMC_chm(iclas)
               wopos(1) = 'w' // trim(val_spec)  // valZ
               call postpr(gesca,wopos,ittim,rutim)
            end do
         end do
         if( INOTMASTER ) nullify(gesca)
      end if

  case( 'SRCUN' )
      !
      ! Unconditional mass fractions source terms for CMC model
      !
      if ( kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1 ) then
         wopos(2)=postp(1)%wopos(2,arrays_number('SRCUN'))
         wopos(3)=postp(1)%wopos(3,arrays_number('SRCUN'))

         do iclas = 1, nspec_uncond_write_CMC_chm
            if( INOTMASTER ) gesca => src_Yk_int_CMC_chm(1:npoin,write_uncond_spec_CMC_chm(iclas))
            if (write_uncond_spec_CMC_chm(iclas) < 100) then
               write(val_spec,'(I2.2)') write_uncond_spec_CMC_chm(iclas)
            else
               write(val_spec,'(I3.3)') write_uncond_spec_CMC_chm(iclas)
            end if
            wopos(1) = 'wU' // trim(val_spec)
            call postpr(gesca,wopos,ittim,rutim)
         end do
         if( INOTMASTER ) nullify(gesca)
      end if

  case( 'VTHER' )
     !
     ! Thermophoretic velocity: VTHER
     !
     if( INOTMASTER )then
        call memgen(zero,ndime,npoin)
        do ipoin=1,npoin
           gevec(1:ndime,ipoin) = Vthermoph_nodes(1:ndime,ipoin)
        end do
     end if

  case( 'TTABL' )
     !
     ! Temperature from table
     !
     ind = 1_ip
     call memgen(zero,max(1_ip,npoin),zero)
     if (kfl_T_tab_index_chm > 0) then
         do ipoin=1,npoin
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

             call fw_lookup( control, scale_control, table_fw, retva, ind )
             gesca(ipoin) = retva(kfl_T_tab_index_chm)
         enddo
     else
         do ipoin=1,npoin
            gesca(ipoin) = 0.0_rp
         enddo
     endif

  case( 'NAMED' )
     !
     ! NAMED unknown of mixed equation model
     !
     if ( kfl_model_chm == 2  ) then
        !
        ! Mixed equation model output unknowns
        !
        wopos(2)=postp(1)%wopos(2,arrays_number('NAMED'))
        wopos(3)=postp(1)%wopos(3,arrays_number('NAMED'))

        do iclas = 1,nclas_chm
           if ( mixedEq_eqs_chm (iclas) % name /=  '' .and.  mixedEq_eqs_chm (iclas) % kfl_do_post /= 0_ip  ) then
              call memgen(zero,npoin,zero)
              if( INOTEMPTY ) then
                  do ipoin = 1,npoin
                     gesca(ipoin) = conce(ipoin,iclas,1)
                  enddo
              endif
              wopos(1) = mixedEq_eqs_chm (iclas) % name
              call postpr(gesca,wopos,ittim,rutim)
              call memgen(2_ip,npoin,zero)
           endif
        end do
        if( INOTMASTER ) nullify(gesca)

        return
     endif

  case( 'PDF  ' )
      !
      ! PDF values for each mixture fraction and node for CMC model
      !
      if ( kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1 ) then
         wopos(2)=postp(1)%wopos(2,arrays_number('PDF  '))
         wopos(3)=postp(1)%wopos(3,arrays_number('PDF  '))

         do imixf = 1, nZ_write_CMC_chm
            ! Compute PDF values at nodes for each mixutre fraction and node
            ! At Z=0 and Z=Zs PDF values 0 or infinite -> leave 0
            aux_val_CMC_chm = 0.0_rp
            do ipoin = 1, npoin
               call chm_PDF_calc_CMC(Zavg_CMC_chm(ipoin), Zvar_CMC_chm(ipoin), &
                       write_iZ_spec_CMC_chm(imixf), aux_val_CMC_chm(ipoin))
            end do
            if( INOTMASTER ) gesca => aux_val_CMC_chm
            write(valZ,'(I2.2)') write_iZ_spec_CMC_chm(imixf)
            wopos(1) = 'PDF' // trim(valZ)
            call postpr(gesca,wopos,ittim,rutim)
         end do
         if( INOTMASTER ) nullify(gesca)
      end if

  case( 'ALPHA' )
      !
      ! Alpha values for boundary and initial conditions
      !

!!!      if ( kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1 ) then
!!!         wopos(2) = 'SCALA'  !! alpha_val_spec_CMC_chm is MATRIX so wopos(2) has to be changed
!!!         wopos(3)=postp(1)%wopos(3,68)

!!!#ifdef CANTERA
!!!         do iclas = 1, nclas_chm
!!!            if( INOTMASTER ) gesca => alpha_val_spec_CMC_chm(1,1:npoin,iclas)
!!!            call getSpeciesName(gas_chm,iclas,aux_name)
!!!            wopos(1) = 'AL' // trim(aux_name)
!!!            call postpr(gesca,wopos,ittim,rutim)
!!!         end do
!!!#endif
!!!         if( INOTMASTER ) nullify(gesca)
!!!
!!!      end if

  case( 'YKSSM' )
      !
      ! Species mass fractions to look up in sectional soot model
      !
      wopos(2)=postp(1)%wopos(2,arrays_number('YKSSM'))
      wopos(3)=postp(1)%wopos(3,arrays_number('YKSSM'))

      if (kfl_lookg_chm > 0 .and. kfl_yk_fw_ssm > 0) then
         if( INOTEMPTY ) then
            nullify(auxvar)
            call memory_alloca(mem_modul(1:2,modul),'AUXVAR', 'chm_outvar',auxvar, lookup_fw(kfl_yk_fw_ssm)%main_table%nvar, npoin)
            nullify(aux_r2p)
            call memory_alloca(mem_modul(1:2,modul),'AUX_R2P','chm_outvar',aux_r2p,nelem)
            do ielem = 1,nelem
               pelty = ltype(ielem)
               pgaus = ngaus(pelty)
               call memory_alloca(mem_modul(1:2,modul),'AUX_R2P','chm_outvar',aux_r2p(ielem)%a,pgaus,&
                   lookup_fw(kfl_yk_fw_ssm)%main_table%nvar)
            end do
            !
            ! Get species
            !
            do iclas = 1,lookup_fw(kfl_yk_fw_ssm)%main_table%nvar
               do ielem = 1,nelem
                  pelty = ltype(ielem)
                  pgaus = ngaus(pelty)
                  aux_r2p(ielem) % a (1:pgaus,iclas) = Yk_ssm_gp(ielem) % a(1:pgaus,iclas,1)
               enddo
            enddo
            call memgen(zero,lookup_fw(kfl_yk_fw_ssm)%main_table%nvar,npoin)
            call smoot5(aux_r2p,auxvar,lookup_fw(kfl_yk_fw_ssm)%main_table%nvar)
            call memory_deallo(mem_modul(1:2,modul),'AUX_R2P','chm_outvar',aux_r2p)
         endif
      endif

      do iclas = 1,lookup_fw(kfl_yk_fw_ssm)%main_table%nvar
         if( INOTEMPTY ) then
            gesca => auxvar(iclas,:)
         endif
         wclas =  trim(intost(iclas))
         if(iclas<10) then
            wopos(1)=postp(1)%wopos(1,arrays_number('YKSSM'))(1:3)//'0'//wclas(1:1)
         else
            wopos_aux=postp(1)%wopos(1,arrays_number('YKSSM'))(1:3) ! Asign in two steps to avoid bug in ibm xlf
            wopos(1)=trim(wopos_aux)//wclas
         end if
         call postpr(gesca,wopos,ittim,rutim)
      end do

      nullify(gesca)

      if (kfl_lookg_chm > 0 .and. kfl_yk_fw_ssm > 0) then
         if( INOTEMPTY ) then
            call memory_deallo(mem_modul(1:2,modul),'AUXVAR','chm_outvar',auxvar)
         endif
      endif

      return

  case( 'QNUCT' )
      !
      ! Discrete sectional method: total soot source due to nucleation:     QNUCT
      !
      if(INOTEMPTY .and. kfl_soot_chm /= 0 ) then
         nullify(aux_r1p)
         call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p,nelem)
         do ielem = 1,nelem
            pelty = ltype(ielem)
            pgaus = ngaus(pelty)
            call memory_alloca(mem_modul(1:2,modul),'AUX_R1P % A','chm_outvar',aux_r1p(ielem)%a,pgaus)
            !
            ! Sum source terms of sections
            !
            do idsm = 1,nsect_ssm
               aux_r1p(ielem) % a = aux_r1p(ielem) % a + RhoC_ssm * Qnucl_gp_chm(ielem) % a (1:pgaus,idsm,1)
            enddo
         end do
         call memgen(zero,npoin,zero)
         call smooth (aux_r1p, gesca)
         call memory_deallo(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p)
      end if
      
  case( 'QCOAT' )
      !
      ! Discrete sectional method: total soot source due to coagulation:    QCOAT
      !
      if(INOTEMPTY .and. kfl_soot_chm /= 0 ) then
         nullify(aux_r1p)
         call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p,nelem)
         do ielem = 1,nelem
            pelty = ltype(ielem)
            pgaus = ngaus(pelty)
            call memory_alloca(mem_modul(1:2,modul),'AUX_R1P % A','chm_outvar',aux_r1p(ielem)%a,pgaus)
            !
            ! Sum source terms of sections
            !
            do idsm = 1,nsect_ssm
               aux_r1p(ielem) % a = aux_r1p(ielem) % a + RhoC_ssm * Qcoag_gp_chm(ielem) % a (1:pgaus,idsm,1)
            enddo
         end do
         call memgen(zero,npoin,zero)
         call smooth (aux_r1p, gesca)
         call memory_deallo(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p)
      end if
  case( 'QCONT' )
      !
      ! Discrete sectional method: total soot source due to condensation:   QCONT
      !
      if(INOTEMPTY .and. kfl_soot_chm /= 0 ) then
         nullify(aux_r1p)
         call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p,nelem)
         do ielem = 1,nelem
            pelty = ltype(ielem)
            pgaus = ngaus(pelty)
            call memory_alloca(mem_modul(1:2,modul),'AUX_R1P % A','chm_outvar',aux_r1p(ielem)%a,pgaus)
            !
            ! Sum source terms of sections
            !
            do idsm = 1,nsect_ssm
               aux_r1p(ielem) % a = aux_r1p(ielem) % a + RhoC_ssm * Qcond_gp_chm(ielem) % a (1:pgaus,idsm,1)
            enddo
         end do
         call memgen(zero,npoin,zero)
         call smooth (aux_r1p, gesca)
         call memory_deallo(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p)
      end if

  case( 'QSURT' )
      !
      ! Discrete sectional method: total soot source due to surface growth: QSURT
      !
      if(INOTEMPTY .and. kfl_soot_chm /= 0 ) then
         nullify(aux_r1p)
         call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p,nelem)
         do ielem = 1,nelem
            pelty = ltype(ielem)
            pgaus = ngaus(pelty)
            call memory_alloca(mem_modul(1:2,modul),'AUX_R1P % A','chm_outvar',aux_r1p(ielem)%a,pgaus)
            !
            ! Sum source terms of sections
            !
            do idsm = 1,nsect_ssm
               aux_r1p(ielem) % a = aux_r1p(ielem) % a + RhoC_ssm * Qsurf_gp_chm(ielem) % a (1:pgaus,idsm,1)
            enddo
         end do
         call memgen(zero,npoin,zero)
         call smooth (aux_r1p, gesca)
         call memory_deallo(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p)
      end if

  case( 'QTOTT' )
      !
      ! Discrete sectional method: total soot source due to all effects: QTOTT
      !
      if(INOTEMPTY .and. kfl_soot_chm /= 0 ) then
         nullify(aux_r1p)
         call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p,nelem)
         do ielem = 1,nelem
            pelty = ltype(ielem)
            pgaus = ngaus(pelty)
            call memory_alloca(mem_modul(1:2,modul),'AUX_R1P % A','chm_outvar',aux_r1p(ielem)%a,pgaus)
            !
            ! Sum source terms of sections
            !
            do idsm = 1,nsect_ssm
               aux_r1p(ielem) % a = aux_r1p(ielem) % a + Qtot_gp_chm(ielem) % a (1:pgaus,idsm,1)
            enddo
         end do
         call memgen(zero,npoin,zero)
         call smooth (aux_r1p, gesca)
         call memory_deallo(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p)
      end if

  case( 'DTRHO' )
    !
    ! DTRHO: Dt * rho
    !
    if(INOTMASTER ) then
        if( kfl_DtRho_tab_index_chm > 0 .and. kfl_lookg_chm > 0) then
           nullify(aux_r1p)
           call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p,nelem)
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              call memory_alloca(mem_modul(1:2,modul),'AUX_R1P % A','chm_outvar',aux_r1p(ielem)%a,pgaus)
              aux_r1p(ielem) % a = DtRho_gp(ielem) % a(:,1,1)
           end do
           call memgen(zero,npoin,zero)
           call smooth (aux_r1p, gesca)
           call memory_deallo(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p)
        end if
    end if

  case( 'AVZFL' )
     !
     ! Average Z flux
     !
     if (cutim > avtim_chm) then
        auxi = cutim - avtim_chm
        if ( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin=1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin)         = av_Z_flux_chm(idime,ipoin) / auxi
                 av_Z_flux_chm(idime,ipoin) = 0.0_rp
              enddo
           end do
        endif
     else
        return
     endif


  case( 'AVYCN' )
     !
     ! Temporal average conditional mass fractions for CMC model
     !
     if ( kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1 ) then
        wopos(2)=postp(1)%wopos(2,arrays_number('AVYCN'))
        wopos(3)=postp(1)%wopos(3,arrays_number('AVYCN'))
        if (cutim > avtim_chm) then
           auxi = cutim - avtim_chm
           do iclas = 1, nspec_cond_write_CMC_chm
              do imixf = 1, nZ_write_CMC_chm
                 if ( INOTMASTER ) then
                    nullify(gesca)
                    do ipoin=1, npoin
                       unkno(ipoin) = av_Yk_CMC_chm(imixf,ipoin,iclas) / auxi
                       av_Yk_CMC_chm(imixf,ipoin,iclas) = 0.0_rp
                    end do
                    if (kfl_post_gp_CMC_chm == 1_ip) then
                       nullify(matr_scal_CMC_chm)
                       matr_scal_CMC_chm => unkno
                       call chm_smooth_field_CMC
                       gesca => aux_val_CMC_chm
                    else
                       gesca => unkno
                    end if
                 endif
                 write(valZ,'(I2.2)') write_iZ_spec_CMC_chm(imixf)
                 write(val_spec,'(I2.2)') write_cond_spec_CMC_chm(iclas)
                 wopos(1) = 'A' // trim(val_spec) // valZ
                 call postpr(gesca,wopos,ittim,rutim)
              end do
           end do
        else
           return
        endif
     endif


  case( 'AVYUN' )
     !
     ! Temporal average unconditional mass fractions for CMC model
     !
     if ( kfl_model_chm == 4 ) then
        wopos(2)=postp(1)%wopos(2,arrays_number('AVYUN'))
        wopos(3)=postp(1)%wopos(3,arrays_number('AVYUN'))
        if (cutim > avtim_chm) then
           nullify(aux_name_ispec)
           if (kfl_solve_cond_CMC_chm == 1 ) then
              aux_name_ispec => write_uncond_spec_CMC_chm
           else
              if (kfl_transfer_condField_CMC_chm == 1_ip .and. nspec_transf_CMC_chm > 0_ip) then
                 aux_name_ispec => transf_spec_CMC_chm
              end if
           end if

           auxi = cutim - avtim_chm
           do iclas = 1, nspec_uncond_write_CMC_chm
              if ( INOTMASTER ) then
                 do ipoin=1, npoin
                    unkno(ipoin) = av_Yk_int_CMC_chm(ipoin,iclas) / auxi
                    av_Yk_int_CMC_chm(ipoin,iclas) = 0.0_rp
                 end do
                 nullify(gesca)
                 gesca => unkno
              endif
              if (aux_name_ispec(iclas) < 100) then
                 write(val_spec,'(I2.2)') aux_name_ispec(iclas)
              else
                 write(val_spec,'(I3.3)') aux_name_ispec(iclas)
              end if
              !!wopos(1) = 'AYU' // trim(val_spec)
              wopos(1) = 'AU' // trim(val_spec)
              call postpr(gesca,wopos,ittim,rutim)
           end do
        else
           return
        endif
     endif


  case( 'AVHCN' )
     !
     ! Temporal average conditional enthalpy for CMC model
     !
     if ( kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1 ) then
        if (kfl_solve_enth_CMC_chm == 1) then
           wopos(2)=postp(1)%wopos(2,arrays_number('AVHCN'))
           wopos(3)=postp(1)%wopos(3,arrays_number('AVHCN'))
           if (cutim > avtim_chm) then
              auxi = cutim - avtim_chm
              do imixf = 1, nZ_write_CMC_chm
                 if ( INOTMASTER ) then
                    nullify(gesca)
                    do ipoin=1, npoin
                       unkno(ipoin) = av_enthalp_CMC_chm(imixf,ipoin) / auxi
                       av_enthalp_CMC_chm(imixf,ipoin) = 0.0_rp
                    end do
                    if (kfl_post_gp_CMC_chm == 1_ip) then
                       nullify(matr_scal_CMC_chm)
                       matr_scal_CMC_chm => unkno
                       call chm_smooth_field_CMC
                       gesca => aux_val_CMC_chm
                    else
                       gesca => unkno
                    end if
                 endif
                 write(valZ,'(I2.2)') write_iZ_spec_CMC_chm(imixf)
                 wopos(1) = 'AHC_' // trim(valZ)
                 call postpr(gesca,wopos,ittim,rutim)
              end do
           else
              return
           endif
        endif
     endif


  case( 'AVHUN' )
     !
     ! Temporal average unconditional enthalpy for CMC model
     !
     if ( kfl_model_chm == 4 ) then
        if (cutim > avtim_chm) then
           auxi = cutim - avtim_chm
           if ( INOTMASTER ) then
              nullify(gesca)
              do ipoin=1, npoin
                 unkno(ipoin) = av_enthalp_int_CMC_chm(ipoin) / auxi
                 av_enthalp_int_CMC_chm(ipoin) = 0.0_rp
              end do
              gesca => unkno
           endif
        else
           return
        endif
     endif


  case( 'AVTCN' )
     !
     ! Temporal average conditional temperature for CMC model
     !
     if ( kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1 ) then
        wopos(2)=postp(1)%wopos(2,arrays_number('AVTCN'))
        wopos(3)=postp(1)%wopos(3,arrays_number('AVTCN'))
        if (cutim > avtim_chm) then
           auxi = cutim - avtim_chm
           do imixf = 1, nZ_write_CMC_chm
              if ( INOTMASTER ) then
                 nullify(gesca)
                 do ipoin=1, npoin
                    unkno(ipoin) = av_temp_CMC_chm(imixf,ipoin) / auxi
                    av_temp_CMC_chm(imixf,ipoin) = 0.0_rp
                 end do
                 if (kfl_post_gp_CMC_chm == 1_ip) then
                    nullify(matr_scal_CMC_chm)
                    matr_scal_CMC_chm => unkno
                    call chm_smooth_field_CMC
                    gesca => aux_val_CMC_chm
                 else
                    gesca => unkno
                 end if
              endif
              write(valZ,'(I2.2)') write_iZ_spec_CMC_chm(imixf)
              wopos(1) = 'ATC' // trim(valZ)
              call postpr(gesca,wopos,ittim,rutim)
           end do
        else
           return
        endif
     endif


  case( 'AVTUN' )
     !
     ! Temporal average unconditional temperature for CMC model
     !
     if ( kfl_model_chm == 4 ) then
        if (cutim > avtim_chm) then
           auxi = cutim - avtim_chm
           if ( INOTMASTER ) then
              nullify(gesca)
              do ipoin=1, npoin
                 unkno(ipoin) = av_temp_int_CMC_chm(ipoin) / auxi
                 av_temp_int_CMC_chm(ipoin) = 0.0_rp
              end do
              gesca => unkno
           endif
        else
           return
        endif
     endif


  case( 'AVVIS' )
     !
     ! Temporal average unconditional laminar viscosity for CMC model
     !
     if ( kfl_model_chm == 4 ) then
        if (cutim > avtim_chm) then
           auxi = cutim - avtim_chm
           if ( INOTMASTER ) then
              nullify(gesca)
              do ipoin=1, npoin
                 unkno(ipoin) = av_visco_lam_int_CMC_chm(ipoin) / auxi
                 av_visco_lam_int_CMC_chm(ipoin) = 0.0_rp
              end do
              gesca => unkno
           endif
        else
           return
        endif
     endif


  case( 'SVF  ' )
     !
     ! Discrete sectional method: total soot volume fraction: SVF
     !
     if(INOTEMPTY ) then
        nullify(aux_r1p)
        call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'AUX_R1P % A','chm_outvar',aux_r1p(ielem)%a,pgaus)
           !
           ! Soot volume fraction
           !
           do idsm = 1,nsect_ssm
              aux_r1p(ielem) % a = aux_r1p(ielem) % a + Qsoot_gp_chm(ielem) % a (1:pgaus,idsm,1)
           enddo
        end do
        call memgen(zero,npoin,zero)
        call smooth (aux_r1p, gesca)
        call memory_deallo(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p)
     end if

  case( 'NDEN ' )
     !
     ! Discrete sectional method: total soot volume fraction: SVF
     !
     if(INOTEMPTY ) then
        nullify(aux_r1p)
        call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'AUX_R1P % A','chm_outvar',aux_r1p(ielem)%a,pgaus)
           !
           ! Soot volume fraction
           !
           do idsm = 1,nsect_ssm
              aux_r1p(ielem) % a = aux_r1p(ielem) % a + NDsoot_gp_chm(ielem) % a (1:pgaus,idsm,1)
           enddo
        end do
        call memgen(zero,npoin,zero)
        call smooth (aux_r1p, gesca)
        call memory_deallo(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p)
     end if

  case( 'AVNAM' )
     !
     ! Averaged NAMED unknown of mixed equation model
     !
     if ( kfl_model_chm == 2  ) then

       if (cutim > avtim_chm) then

          !
          ! Time interval
          !
          auxi = cutim - avtim_chm

          !
          ! Mixed equation model output unknowns
          !
          wopos(2)=postp(1)%wopos(2,arrays_number('AVNAM'))
          wopos(3)=postp(1)%wopos(3,arrays_number('AVNAM'))

          do iclas = 1,nclas_chm
             if ( mixedEq_eqs_chm (iclas) % name /=  '' .and.  mixedEq_eqs_chm (iclas) % kfl_do_post /= 0_ip  ) then
                call memgen(zero,npoin,zero)
                if( INOTEMPTY ) then
                    do ipoin = 1,npoin
                       gesca(ipoin) = av_Named_Unk_chm(ipoin,iclas) / auxi
                       av_Named_Unk_chm(ipoin,iclas) = 0.0_rp
                    enddo
                endif
                wclas =  trim(intost(iclas))

                if(iclas<10) then
                   wopos(1)='AV'//'00'//wclas(1:1)
                else if(iclas<100) then
                   wopos(1)='AV'//'0'//wclas
                else
                   wopos(1)='AV'//wclas
                end if
                call postpr(gesca,wopos,ittim,rutim)
                call memgen(2_ip,npoin,zero)
             endif
          end do

          if( INOTMASTER ) nullify(gesca)

          return

       else
          return
       endif

     endif

  case( 'TCHMI' )
     !
     ! Time spent in chemical integration for CMC model
     !
     if ( kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1 ) then
         wopos(2) = postp(1)%wopos(2,arrays_number('SVF  '))
         wopos(3) = postp(1)%wopos(3,arrays_number('SVF  '))

        do imixf = 1, nZ_write_CMC_chm
           if ( INOTMASTER ) then
              nullify(gesca)
              gesca => t_chem_integ_CMC_chm(write_iZ_spec_CMC_chm(imixf),1:npoin)
           endif
           write(valZ,'(I2.2)') write_iZ_spec_CMC_chm(imixf)
           wopos(1) = 'CI_' // trim(valZ)
           call postpr(gesca,wopos,ittim,rutim)
        end do
        if( INOTMASTER ) nullify(gesca)
     endif

  case( 'PRODU' )
     !
     ! PRODU, Production part of reactive source term
     !
     !
     ! Update based on current unknown
     !
     if (kfl_model_chm == 2 .and. kfl_max_srcfw_chm > 0) then
        !
        ! Taulated sources
        !
        if (kfl_lookg_chm > 0) then
            if (kfl_multimod_chm > 0_ip) then
               call chm_gp_multi_reatab(ITASK_BEGSTE)
            else
               call chm_gp_reatab(ITASK_BEGSTE)
            end if
        else
            call chm_reatab()
        endif

        !
        ! Mixed equation model, output sources that are non-passive
        !
        wopos(2)=postp(1)%wopos(2,arrays_number('SOURC'))
        wopos(3)=postp(1)%wopos(3,arrays_number('SOURC'))

        if (kfl_lookg_chm > 0) then
           if( INOTEMPTY ) then
              nullify(auxvar)
              call memory_alloca(mem_modul(1:2,modul),'AUXVAR', 'chm_outvar',auxvar, nclas_chm, npoin)
              nullify(aux_r2p)
              call memory_alloca(mem_modul(1:2,modul),'AUX_R2P','chm_outvar',aux_r2p,nelem)
              do ielem = 1,nelem
                 pelty = ltype(ielem)
                 pgaus = ngaus(pelty)
                 call memory_alloca(mem_modul(1:2,modul),'AUX_R2P % A','chm_outvar',aux_r2p(ielem)%a,pgaus,nclas_chm)
              end do
              !
              ! Get sources
              !
              do iclas = 1,nclas_chm
                 if ( mixedEq_eqs_chm (iclas) % kfl_source_type ==  CHM_SRC_TABLE ) then
                    !
                    ! Tabulated source terms
                    !
                    do ielem = 1,nelem
                       pelty = ltype(ielem)
                       pgaus = ngaus(pelty)
                       aux_r2p(ielem) % a (1:pgaus,iclas) = mass_gp(ielem) % a(1:pgaus,iclas,1)
                    enddo
                 endif
              enddo

              call smoot5(aux_r2p,auxvar,nclas_chm)
              call memory_deallo(mem_modul(1:2,modul),'AUX_R2P','chm_outvar',aux_r2p)
           endif

        endif

        do iclas = 1,nclas_chm
           if ( mixedEq_eqs_chm(iclas) % kfl_source_split == 1 ) then
              if( INOTEMPTY ) then
                 if (kfl_lookg_chm > 0) then
                    gesca => auxvar(iclas,:)
                 endif
              endif

              !
              ! Get name
              !
              wclas =  trim(intost(iclas))
              if(iclas<10) then
                 wopos(1)='PSR'//'0'//wclas(1:1)
              else
                 wopos(1)='PSR'//wclas
              end if

              !
              ! Write postprocessed data
              !
              call postpr(gesca,wopos,ittim,rutim)
           endif
        end do
        if( INOTMASTER ) nullify(gesca)

        if (kfl_lookg_chm > 0) then
           if( INOTEMPTY ) then
              call memory_deallo(mem_modul(1:2,modul),'AUXVAR','chm_outvar',auxvar)
           endif
        endif

        return
     endif


  case( 'CONSU' )
     !
     ! CONSU, Consumption part of reactive source term
     !
     !
     ! Update based on current unknown
     !
     if (kfl_model_chm == 2 .and. kfl_max_srcfw_chm > 0) then
        !
        ! Taulated sources
        !
        if (kfl_lookg_chm > 0) then
            if (kfl_multimod_chm > 0_ip) then
              call chm_gp_multi_reatab(ITASK_BEGSTE)
            else
              call chm_gp_reatab(ITASK_BEGSTE)
            end if
        else
            call chm_reatab()
        endif

        !
        ! Mixed equation model, output sources that are non-passive
        !
        wopos(2)=postp(1)%wopos(2,arrays_number('SOURC'))
        wopos(3)=postp(1)%wopos(3,arrays_number('SOURC'))

        if (kfl_lookg_chm > 0) then
           if( INOTEMPTY ) then
              nullify(auxvar)
              call memory_alloca(mem_modul(1:2,modul),'AUXVAR', 'chm_outvar',auxvar, nclas_chm, npoin)
              nullify(aux_r2p)
              call memory_alloca(mem_modul(1:2,modul),'AUX_R2P','chm_outvar',aux_r2p,nelem)
              do ielem = 1,nelem
                 pelty = ltype(ielem)
                 pgaus = ngaus(pelty)
                 call memory_alloca(mem_modul(1:2,modul),'AUX_R2P % A','chm_outvar',aux_r2p(ielem)%a,pgaus,nclas_chm)
              end do

              !
              ! Consumption is the contribution of the split source term
              !
              do ielem = 1,nelem
                 !
                 ! Get elemental quantities
                 !
                 pelty = ltype(ielem)
                 pnode = nnode(pelty)
                 pgaus = ngaus(pelty)
                 lnods_loc(1:pnode) = lnods(1:pnode,ielem)
                 call chm_post_gather(&
                      pnode,lnods_loc,elcon(1:pnode,:,:),elcod,dummr)

                 !
                 ! Mass fractions on Gaussian integration points
                 !
                 gpcon = 0.0_rp
                 do iclas = 1,nclas_chm
                    do igaus = 1,pgaus
                       do inode = 1,pnode
                          gpcon(igaus,iclas) = gpcon(igaus,iclas)&
                               + elmar(pelty)%shape(inode,igaus) * elcon(inode,iclas,1)
                       end do
                    end do
                 end do

                 do iclas = 1,nclas_chm
                    if ( mixedEq_eqs_chm (iclas) % kfl_source_type ==  CHM_SRC_TABLE ) then
                       if ( mixedEq_eqs_chm(iclas) % kfl_source_split == 1) then
                           !
                           ! Add contribution of first order dependence (consumption)
                           !
                           aux_r2p(ielem) % a (1:pgaus,iclas) = gpcon(1:pgaus,iclas) * massConsumption_gp(ielem) % a(1:pgaus,iclas,1)
                       endif
                    endif
                 enddo
              enddo

              call smoot5(aux_r2p,auxvar,nclas_chm)
              call memory_deallo(mem_modul(1:2,modul),'AUX_R2P','chm_outvar',aux_r2p)
           endif

        endif

        do iclas = 1,nclas_chm
           if ( mixedEq_eqs_chm(iclas) % kfl_source_split == 1 ) then
              if( INOTEMPTY ) then
                 if (kfl_lookg_chm > 0) then
                    gesca => auxvar(iclas,:)
                 endif
              endif

              !
              ! Get name
              !
              wclas =  trim(intost(iclas))
              if(iclas<10) then
                 wopos(1)='CSR'//'0'//wclas(1:1)
              else
                 wopos(1)='CSR'//wclas
              end if

              !
              ! Write postprocessed data
              !
              call postpr(gesca,wopos,ittim,rutim)
           endif
        end do
        if( INOTMASTER ) nullify(gesca)

        if (kfl_lookg_chm > 0) then
           if( INOTEMPTY ) then
              call memory_deallo(mem_modul(1:2,modul),'AUXVAR','chm_outvar',auxvar)
           endif
        endif

        return
     endif

  case( 'FLAME' )
    !
    ! Flame index
    !
    if(INOTMASTER ) then
       nullify(aux_r1p)
       call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p,nelem)
       do ielem = 1,nelem
          pelty = ltype(ielem)
          pgaus = ngaus(pelty)
          call memory_alloca(mem_modul(1:2,modul),'AUX_R1P % A','chm_outvar',aux_r1p(ielem)%a,pgaus)
          aux_r1p(ielem) % a = flame_index_gp(ielem) % a(:,1)
       end do
       call memgen(zero,npoin,zero)
       call smooth (aux_r1p, gesca)
       call memory_deallo(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p)
    end if
  case( 'ZGRAD' )
    !
    ! ZGRAD - looked up from diffusion table
    !
    if(INOTMASTER ) then
       nullify(aux_r1p)
       call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p,nelem)
       do ielem = 1,nelem
          pelty = ltype(ielem)
          pgaus = ngaus(pelty)
          call memory_alloca(mem_modul(1:2,modul),'AUX_R1P % A','chm_outvar',aux_r1p(ielem)%a,pgaus)
          aux_r1p(ielem) % a = zgrad_gp(ielem) % a(:,1)
       end do
       call memgen(zero,npoin,zero)
       call smooth (aux_r1p, gesca)
       call memory_deallo(mem_modul(1:2,modul),'AUX_R1P','chm_outvar',aux_r1p)
    end if
  end select

  call outvar(&
       ivari,&
       ittim,rutim,postp(1)%wopos(:,ivari),MESH_ID=imesh)

end subroutine chm_outvar
