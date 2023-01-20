!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Chemic
!> @{
!> @file    chm_averag.f90
!> @author  Daniel Mira
!> @date    13/06/204
!> @brief   Average variables: temperature
!> @details Average variables: temperature
!> @}
!-----------------------------------------------------------------------
subroutine chm_averag()
  use def_master,             only : cutim, dtime, INOTEMPTY, INOTMASTER, mass_sink, conce, advec, densi, tempe, visco, &
                                     mem_modul, modul
  use def_domain,             only : ndime, nelem, npoin, ltype, ngaus
  use def_chemic,             only : mixedeq_eqs_chm, posttable_fw, avtim_chm, kfl_icmean_chm, kfl_icvar_chm, kfl_izmean_chm,&
                                     kfl_izvar_chm, kfl_model_chm, kfl_solve_cond_CMC_chm, kfl_solve_enth_CMC_chm, kfl_spray_chm,&
                                     nclas_chm, nspec_cond_write_CMC_chm, nspec_uncond_write_CMC_chm, nZ_write_CMC_chm, avY_chm,&
                                     avYv_chm, avZv_chm, avZ_chm, avZ2_chm, avY2_chm, avL_chm, avL2_chm, avS_chm, Sigma_chm,&
                                     avS0_chm, Sigm0_chm, avd32_chm, d32_chm, densi_int_CMC_chm, avden_chm,&
                                     av_Z_flux_chm, hrr_avg_chm, hrr_chm, avmsk_chm, avposttab_chm, write_cond_spec_CMC_chm,&
                                     av_Yk_CMC_chm, write_iZ_spec_CMC_chm, write_uncond_spec_CMC_chm, av_enthalp_CMC_chm,&
                                     av_enthalp_int_CMC_chm, enthalp_int_CMC_chm, av_temp_CMC_chm, av_temp_int_CMC_chm,&
                                     temp_int_CMC_chm, av_visco_lam_int_CMC_chm, Yk_CMC_chm, av_Yk_int_CMC_chm, Yk_int_CMC_chm,&
                                     enthalp_CMC_chm, visco_lam_int_CMC_chm, temp_CMC_chm
  use mod_ker_proper,         only : ker_proper
  use def_kintyp,             only : ip,rp,r2p
  use mod_memory,             only : memory_alloca, memory_deallo
  use mod_solver,             only : solver_lumped_mass_system
  use mod_output_postprocess, only : output_postprocess_check_variable_postprocess
  use mod_chm_mixedEq,        only : CHM_EQ_ZVAR
  use mod_chm_mixedEq,        only : CHM_EQ_ZZ
  use mod_chm_mixedEq,        only : CHM_EQ_YCVAR
  use mod_chm_mixedEq,        only : CHM_EQ_YCYC
  implicit none
  integer(ip) :: ipoin,iclas_phi,dummi,ielem,pelty,pgaus,ipostvar,idime
  integer(ip) :: iclas, imixf, pos_clas
  real(rp)    :: zechm
  real(rp), pointer         :: prope_tmp(:)
  real(rp), pointer         :: auxvar(:,:)
  type(r2p),pointer         :: aux_r2p(:)

  external                  :: runend
  external                  :: chm_post_gp_lookup
  external                  :: smoot5

  zechm = epsilon(0.0_rp)

  if( cutim > avtim_chm ) then

     if( INOTMASTER ) then

        !
        ! AVY: average reaction progress Yc or C
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVY  ') .and. kfl_icmean_chm > 0_ip  ) then
          do ipoin=1,npoin
             avY_chm(ipoin) = avY_chm(ipoin)&
                                    + conce(ipoin,kfl_icmean_chm,1) * dtime
          end do
        end if

        !
        ! AVYv: average variance of reaction progress Yc
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVYV ') .and. kfl_icvar_chm > 0_ip ) then
          !
          ! Variance is transported directly
          !
          if ( mixedEq_eqs_chm(kfl_icvar_chm) % kfl_eqtype == CHM_EQ_YCVAR ) then
             do ipoin=1,npoin
                avYv_chm(ipoin) = avYv_chm(ipoin)&
                                         + conce(ipoin,kfl_icvar_chm,1) * dtime
             end do
          !
          ! Yc^2 is transported and variance must be computed
          !
          else if ( mixedEq_eqs_chm(kfl_icvar_chm) % kfl_eqtype == CHM_EQ_YCYC ) then
             do ipoin=1,npoin
                avYv_chm(ipoin) = avYv_chm(ipoin)&
                                + ( conce(ipoin,kfl_icvar_chm,1) - conce(ipoin,kfl_icmean_chm,1)*conce(ipoin,kfl_icmean_chm,1) ) &
                                * dtime
             end do
          end if

        end if

        !
        ! AVZv: average variance of mixture fraction
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVZV ') .and. kfl_izvar_chm>0_ip ) then
          !
          ! Variance is transported directly
          !
          if ( mixedEq_eqs_chm(kfl_izvar_chm) % kfl_eqtype == CHM_EQ_ZVAR ) then
             do ipoin=1,npoin
                avZv_chm(ipoin) = avZv_chm(ipoin)&
                                         + conce(ipoin,kfl_izvar_chm,1) * dtime
             end do
          !
          ! Z^2 is transported and variance must be computed
          !
          else if ( mixedEq_eqs_chm(kfl_izvar_chm) % kfl_eqtype == CHM_EQ_ZZ ) then
             do ipoin=1,npoin
                avZv_chm(ipoin) = avZv_chm(ipoin)&
                                + ( conce(ipoin,kfl_izvar_chm,1) - conce(ipoin,kfl_izmean_chm,1)*conce(ipoin,kfl_izmean_chm,1) ) &
                                * dtime
             end do
          end if

        end if

        !
        ! AVZ: average mixture fraction Z
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVZ  ') .and. kfl_izmean_chm>0_ip ) then
          do ipoin=1,npoin
             avZ_chm(ipoin) = avZ_chm(ipoin)&
                                    + conce(ipoin,kfl_izmean_chm,1) * dtime
          end do
        end if

        !
        ! AVZ2: average squared of mixture fraction Z*Z
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVZ2 ') .and. kfl_izmean_chm>0_ip ) then
          do ipoin=1,npoin
             avZ2_chm(ipoin) = avZ2_chm(ipoin)&
                                      + conce(ipoin,kfl_izmean_chm,1)*conce(ipoin,kfl_izmean_chm,1) * dtime
          end do
        end if

        !
        ! AVY2: average squared of progress variable Yc*Yc or C*C
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVY2 ') .and. kfl_icmean_chm>0_ip ) then
          do ipoin=1,npoin
             avY2_chm(ipoin) = avY2_chm(ipoin)&
                                      + conce(ipoin,kfl_icmean_chm,1)*conce(ipoin,kfl_icmean_chm,1) * dtime
          end do
        end if

        !
        ! AVL: average liquid volume fraction phi_L
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVL  ') ) then

          if ( kfl_spray_chm /= 0 ) then
             iclas_phi = nclas_chm - 1
             do ipoin=1,npoin
                   avL_chm(ipoin) = avL_chm(ipoin)&
                                      + conce(ipoin,iclas_phi,1) * dtime
             end do
          else
             call runend('CHEMIC CHM_AVERAG: Liquid volume fraction only valid with spray model')
          end if
        end if

        !
        ! AVL2: average liquid volume fraction squared phi_L*phi_L
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVL2 ') ) then

          if ( kfl_spray_chm /= 0 ) then
             iclas_phi = nclas_chm - 1
             do ipoin=1,npoin
                   avL2_chm(ipoin) = avL2_chm(ipoin)&
                                            + conce(ipoin,iclas_phi,1)*conce(ipoin,iclas_phi,1) * dtime
             end do
          else
             call runend('CHEMIC CHM_AVERAG: Squared of liquid volume fraction only valid with spray model')
          end if
        end if

        !
        ! AVS: Average interface surface density Sigma
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVS  ') ) then

          do ipoin=1,npoin
             avS_chm(ipoin) = avS_chm(ipoin)&
                                      + dtime*Sigma_chm(ipoin)
          end do

        end if

        !
        ! AVS0: Average interface surface density Sigma_0
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVS0 ') ) then

          do ipoin=1,npoin
             avS0_chm(ipoin) = avS0_chm(ipoin)&
                                      + dtime*Sigm0_chm(ipoin)
          end do

        end if

        !
        ! AVD32: Average Sauter mean diameter
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVD32') ) then

          do ipoin=1,npoin
             avd32_chm(ipoin) = avd32_chm(ipoin)&
                                      + dtime*d32_chm(ipoin)
          end do

        end if

        !
        ! AVDEN: average density
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVDEN') ) then

           if (kfl_model_chm == 4) then
              if (kfl_solve_cond_CMC_chm == 1) then
                 do ipoin=1,npoin
                    avden_chm(ipoin) = avden_chm(ipoin) &
                         + densi_int_CMC_chm(ipoin) * dtime
                 end do
              else
                 do ipoin=1,npoin
                    avden_chm(ipoin) = avden_chm(ipoin) &
                         + densi(ipoin,1) * dtime
                 end do
              end if
           else
              nullify ( prope_tmp )
              call memory_alloca(mem_modul(1:2,modul),'PROPE_TMP', 'chm_averag',prope_tmp, npoin)
              call ker_proper('DENSI','NPOIN',dummi,dummi,prope_tmp)

              do ipoin=1,npoin
                 avden_chm(ipoin) = avden_chm(ipoin) + prope_tmp(ipoin) * dtime
              end do

              call memory_deallo(mem_modul(1:2,modul),'PROPE_TMP', 'chm_averag',prope_tmp)
           end if

        end if

        !
        ! AVZFL: average Z flux
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVZFL') ) then

           nullify ( prope_tmp )
           call memory_alloca(mem_modul(1:2,modul),'PROPE_TMP', 'chm_averag',prope_tmp, npoin)
           call ker_proper('DENSI','NPOIN',dummi,dummi,prope_tmp)

           do ipoin=1,npoin
              do idime = 1, ndime
                 av_Z_flux_chm(idime,ipoin) = av_Z_flux_chm(idime,ipoin) &
                     &                        + prope_tmp(ipoin) * conce(ipoin,kfl_izmean_chm,1) * advec(idime,ipoin,1) * dtime
              enddo
           enddo

           call memory_deallo(mem_modul(1:2,modul),'PROPE_TMP', 'chm_averag',prope_tmp)

        end if

        !
        ! AVHRR: average heat release
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVHRR') ) then
           do ipoin=1,npoin
              hrr_avg_chm(ipoin) = hrr_avg_chm(ipoin) + hrr_chm(ipoin) * dtime
           end do
        end if

        !
        ! AVMSK: average mass source term from spray
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVMSK') ) then

           nullify ( prope_tmp )
           call memory_alloca(mem_modul(1:2,modul),'PROPE_TMP','chm_averag',prope_tmp,npoin)
           if (associated(mass_sink)) then
              do ipoin=1,npoin
                 prope_tmp(ipoin) = mass_sink(ipoin)
              enddo
              call solver_lumped_mass_system(1_ip,prope_tmp,EXCHANGE=.false.)
           else
              do ipoin=1,npoin
                 prope_tmp(ipoin) = 0.0_rp
              end do
           endif

           do ipoin=1,npoin
              avmsk_chm(ipoin) = avmsk_chm(ipoin) + prope_tmp(ipoin) * dtime
           end do
           call memory_deallo(mem_modul(1:2,modul),'PROPE_TMP','chm_averag',prope_tmp)

        end if


        !
        ! AVPOT: average postprocessing table quantities
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVPOT') ) then
           if( INOTEMPTY ) then
              nullify(auxvar)
              call memory_alloca(mem_modul(1:2,modul),'AUXVAR', 'chm_averag',auxvar, posttable_fw % main_table % nvar,  npoin)
              nullify(aux_r2p)
              call memory_alloca(mem_modul(1:2,modul),'AUX_R2P','chm_averag',aux_r2p,nelem)
              do ielem = 1,nelem
                 pelty = ltype(ielem)
                 pgaus = ngaus(pelty)
                 call memory_alloca(mem_modul(1:2,modul),'AUX_R2P % A','chm_averag',aux_r2p(ielem)%a,pgaus,&
                     posttable_fw % main_table % nvar)
              end do
              !
              ! Lookup from postprocessing table
              !
              call chm_post_gp_lookup(aux_r2p,posttable_fw)
              call smoot5(aux_r2p,auxvar,posttable_fw % main_table % nvar)
              call memory_deallo(mem_modul(1:2,modul),'AUX_R2P','chm_averag',aux_r2p)
           endif

           do ipostvar=1,posttable_fw % main_table % nvar
              do ipoin=1,npoin
                 avposttab_chm(ipoin,ipostvar) = avposttab_chm(ipoin,ipostvar)&
                                    + auxvar(ipostvar,ipoin) * dtime
              end do
           enddo
           if( INOTEMPTY ) then
              call memory_deallo(mem_modul(1:2,modul),'AUXVAR','chm_averag',auxvar)
           endif
        end if


        !
        ! AVYCN: temp. average conditional mass fractions for CMC model
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVYCN') ) then
           do iclas = 1, nspec_cond_write_CMC_chm
              pos_clas = write_cond_spec_CMC_chm(iclas)
              do ipoin=1,npoin
                 do imixf = 1, nZ_write_CMC_chm
                    av_Yk_CMC_chm(imixf,ipoin,iclas) = av_Yk_CMC_chm(imixf,ipoin,iclas) &
                         + Yk_CMC_chm(write_iZ_spec_CMC_chm(imixf),ipoin,pos_clas) * dtime
                 end do
              end do
           end do
        end if


        !
        ! AVYUN: temp. average unconditional mass fractions for CMC model
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVYUN') ) then
           do iclas = 1, nspec_uncond_write_CMC_chm
              pos_clas = write_uncond_spec_CMC_chm(iclas)
              do ipoin=1,npoin
                 av_Yk_int_CMC_chm(ipoin,iclas) = av_Yk_int_CMC_chm(ipoin,iclas) &
                      + Yk_int_CMC_chm(ipoin,pos_clas) * dtime
              end do
           end do
        end if


        !
        ! AVHCN: temp. average conditional enthalpy for CMC model
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVHCN') .and. kfl_solve_enth_CMC_chm == 1) then
           do ipoin=1,npoin
              do imixf = 1, nZ_write_CMC_chm
                 av_enthalp_CMC_chm(imixf,ipoin) = av_enthalp_CMC_chm(imixf,ipoin) &
                      + enthalp_CMC_chm(write_iZ_spec_CMC_chm(imixf),ipoin) * dtime
              end do
           end do
        end if


        !
        ! AVHUN: temp. average unconditional enthalpy for CMC model
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVHUN') ) then
           do ipoin=1,npoin
              av_enthalp_int_CMC_chm(ipoin) = av_enthalp_int_CMC_chm(ipoin) &
                   + enthalp_int_CMC_chm(ipoin) * dtime
           end do
        end if


        !
        ! AVTCN: temp. average conditional temperature for CMC model
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVTCN') ) then
           do ipoin=1,npoin
              do imixf = 1, nZ_write_CMC_chm
                 av_temp_CMC_chm(imixf,ipoin) = av_temp_CMC_chm(imixf,ipoin) &
                      + temp_CMC_chm(write_iZ_spec_CMC_chm(imixf),ipoin) * dtime
              end do
           end do
        end if


        !
        ! AVTUN: temp. average unconditional temperature for CMC model
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVTUN') ) then
           if (kfl_solve_cond_CMC_chm == 1_ip) then
              do ipoin=1,npoin
                 av_temp_int_CMC_chm(ipoin) = av_temp_int_CMC_chm(ipoin) &
                      + temp_int_CMC_chm(ipoin) * dtime
              end do
           else
              do ipoin=1,npoin
                 av_temp_int_CMC_chm(ipoin) = av_temp_int_CMC_chm(ipoin) &
                      + tempe(ipoin,1) * dtime
              end do

           end if
        end if


        !
        ! AVVIS: temp. average unconditional viscosity for CMC model
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVVIS') ) then
           if (kfl_solve_cond_CMC_chm == 1_ip) then
              do ipoin=1,npoin
                 av_visco_lam_int_CMC_chm(ipoin) = av_visco_lam_int_CMC_chm(ipoin) &
                      + visco_lam_int_CMC_chm(ipoin) * dtime
              end do
           else
              do ipoin=1,npoin
                 av_visco_lam_int_CMC_chm(ipoin) = av_visco_lam_int_CMC_chm(ipoin) &
                      + visco(ipoin,1) * dtime
              end do
           end if
        end if


     end if

  end if

end subroutine chm_averag
