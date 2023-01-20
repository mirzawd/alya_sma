!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Chemic
!> @{
!> @file    chm_averag_fast.f90
!> @author  Daniel Mira y Guillermo Oyarzun
!> @date    11/02/2021
!> @brief   Average variables: temperature
!> @details Average variables: temperature , vectorization friendly
!> @}
!-----------------------------------------------------------------------
subroutine chm_averag_fast()
  use def_master,             only : conce, cutim, dtime, INOTEMPTY, INOTMASTER, mass_sink, mem_modul, modul, advec
  use def_domain,             only : ndime, nelem, npoin, ltype, ngaus
  use def_chemic,             only : avY_chm, mixedEq_eqs_chm, posttable_fw, avtim_chm, kfl_icmean_chm, kfl_icvar_chm,&
                                     kfl_izmean_chm, kfl_izvar_chm, kfl_spray_chm, nclas_chm, avYv_chm, avZv_chm, avZ_chm,&
                                     avZ2_chm, avY2_chm, avL_chm, avL2_chm, avS_chm, Sigma_chm, avS0_chm, Sigm0_chm, avd32_chm,&
                                     d32_chm, avden_chm, av_Z_flux_chm, hrr_avg_chm, hrr_chm, avmsk_chm, av_Named_Unk_chm, &
                                     avposttab_chm
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
  integer(ip)               :: ipoin,iclas_phi,dummi,ielem,pelty,pgaus,ipostvar,idime
  real(rp)                  :: zechm
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
          avY_chm(1:npoin) = avY_chm(1:npoin)&
                                    + conce(1:npoin,kfl_icmean_chm,1) * dtime
        end if

        !
        ! AVYv: average variance of reaction progress Yc
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVYV ') .and. kfl_icvar_chm > 0_ip ) then
          !
          ! Variance is transported directly
          !
          if ( mixedEq_eqs_chm(kfl_icvar_chm) % kfl_eqtype == CHM_EQ_YCVAR ) then
             avYv_chm(1:npoin) = avYv_chm(1:npoin)&
                                         + conce(1:npoin,kfl_icvar_chm,1) * dtime

          !
          ! Yc^2 is transported and variance must be computed
          !
          else if ( mixedEq_eqs_chm(kfl_icvar_chm) % kfl_eqtype == CHM_EQ_YCYC ) then
             avYv_chm(1:npoin) = avYv_chm(1:npoin)&
                                + ( conce(1:npoin,kfl_icvar_chm,1) - conce(1:npoin,kfl_icmean_chm,1)&
                                *conce(1:npoin,kfl_icmean_chm,1) ) * dtime

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
             avZv_chm(1:npoin) = avZv_chm(1:npoin)&
                                         + conce(1:npoin,kfl_izvar_chm,1) * dtime

          !
          ! Z^2 is transported and variance must be computed
          !
          else if ( mixedEq_eqs_chm(kfl_izvar_chm) % kfl_eqtype == CHM_EQ_ZZ ) then
            avZv_chm(1:npoin) = avZv_chm(1:npoin)&
                                + (conce(1:npoin,kfl_izvar_chm,1) - conce(1:npoin,kfl_izmean_chm,1)*conce(1:npoin,kfl_izmean_chm,1))&
                                * dtime
          end if

        end if

        !
        ! AVZ: average mixture fraction Z
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVZ  ') .and. kfl_izmean_chm>0_ip ) then
          avZ_chm(1:npoin) = avZ_chm(1:npoin)&
                                    + conce(1:npoin,kfl_izmean_chm,1) * dtime

        end if

        !
        ! AVZ2: average squared of mixture fraction Z*Z
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVZ2 ') .and. kfl_izmean_chm>0_ip ) then
          avZ2_chm(1:npoin) = avZ2_chm(1:npoin)&
                                      + conce(1:npoin,kfl_izmean_chm,1)*conce(1:npoin,kfl_izmean_chm,1) * dtime

        end if

        !
        ! AVY2: average squared of progress variable Yc*Yc or C*C
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVY2 ') .and. kfl_icmean_chm>0_ip ) then
             avY2_chm(1:npoin) = avY2_chm(1:npoin)&
                                      + conce(1:npoin,kfl_icmean_chm,1)*conce(1:npoin,kfl_icmean_chm,1) * dtime

        end if

        !
        ! AVL: average liquid volume fraction phi_L
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVL  ') ) then

          if ( kfl_spray_chm /= 0 ) then
             iclas_phi = nclas_chm - 1
             avL_chm(1:npoin) = avL_chm(1:npoin)&
                                      + conce(1:npoin,iclas_phi,1) * dtime

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
             avL2_chm(1:npoin) = avL2_chm(1:npoin)&
                                            + conce(1:npoin,iclas_phi,1)*conce(1:npoin,iclas_phi,1) * dtime

          else
             call runend('CHEMIC CHM_AVERAG: Squared of liquid volume fraction only valid with spray model')
          end if
        end if

        !
        ! AVS: Average interface surface density Sigma
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVS  ') ) then
          avS_chm(1:npoin) = avS_chm(1:npoin)&
                                      + dtime*Sigma_chm(1:npoin)


        end if

        !
        ! AVS0: Average interface surface density Sigma_0
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVS0 ') ) then

          avS0_chm(1:npoin) = avS0_chm(1:npoin)&
                                      + dtime*Sigm0_chm(1:npoin)


        end if

        !
        ! AVD32: Average Sauter mean diameter
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVD32') ) then

          avd32_chm(1:npoin) = avd32_chm(1:npoin)&
                                      + dtime*d32_chm(1:npoin)


        end if

        !
        ! AVDEN: average density
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVDEN') ) then

           nullify ( prope_tmp )
           call memory_alloca(mem_modul(1:2,modul),'PROPE_TMP', 'chm_averag',prope_tmp, npoin)
           call ker_proper('DENSI','NPOIN',dummi,dummi,prope_tmp)

           avden_chm(1:npoin) = avden_chm(1:npoin) + prope_tmp(1:npoin) * dtime

           call memory_deallo(mem_modul(1:2,modul),'PROPE_TMP', 'chm_averag',prope_tmp)

        end if

        !
        ! AVZFL: average Z flux
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVZFL') ) then

           nullify ( prope_tmp )
           call memory_alloca(mem_modul(1:2,modul),'PROPE_TMP', 'chm_averag',prope_tmp, npoin)
           call ker_proper('DENSI','NPOIN',dummi,dummi,prope_tmp)

           do idime = 1, ndime
              av_Z_flux_chm(idime,1:npoin) = av_Z_flux_chm(idime,1:npoin) &
                  &                        + prope_tmp(1:npoin) * conce(1:npoin,kfl_izmean_chm,1) * advec(idime,1:npoin,1) * dtime
           enddo

           call memory_deallo(mem_modul(1:2,modul),'PROPE_TMP', 'chm_averag',prope_tmp)

        end if

        !
        ! AVHRR: average heat release
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVHRR') ) then

           hrr_avg_chm(1:npoin) = hrr_avg_chm(1:npoin) + hrr_chm(1:npoin) * dtime
        end if

        !
        ! AVMSK: average mass source term from spray
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVMSK') ) then

           nullify ( prope_tmp )
           call memory_alloca(mem_modul(1:2,modul),'PROPE_TMP','chm_averag',prope_tmp,npoin)
           if (associated(mass_sink)) then

              prope_tmp(1:npoin) = mass_sink(1:npoin)

              call solver_lumped_mass_system(1_ip,prope_tmp,EXCHANGE=.false.)
           else

              prope_tmp(1:npoin) = 0.0_rp

           endif

           avmsk_chm(1:npoin) = avmsk_chm(1:npoin) + prope_tmp(1:npoin) * dtime

           call memory_deallo(mem_modul(1:2,modul),'PROPE_TMP','chm_averag',prope_tmp)

        end if

        !
        ! AVNAM: average heat release
        !
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVNAM') ) then
           do ipoin=1,npoin
              av_Named_Unk_chm(ipoin,1:nclas_chm) = av_Named_Unk_chm(ipoin,1:nclas_chm)     &
                                                   + conce(ipoin,1:nclas_chm,1) * dtime
           end do
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
              avposttab_chm(1:npoin,ipostvar) = avposttab_chm(1:npoin,ipostvar)&
                                    + auxvar(ipostvar,1:npoin) * dtime
           enddo
           if( INOTEMPTY ) then
              call memory_deallo(mem_modul(1:2,modul),'AUXVAR','chm_averag',auxvar)
           endif
        end if



     end if

  end if

end subroutine chm_averag_fast
