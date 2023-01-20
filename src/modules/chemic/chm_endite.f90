!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_endite(itask)
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_endite
  ! NAME
  !    chm_endite
  ! DESCRIPTION
  !    This routine checks convergence and performs updates of the
  !    temperature  at:
  !    - itask=1 The end of an internal iteration
  !    - itask=2 The end of the internal loop iteration
  ! USES
  !    chm_cvgunk
  !    chm_updunk
  ! USED BY
  !    chm_doiter
  !***
  !-----------------------------------------------------------------------
  use def_master,                 only : ITASK_ENDINN, ITASK_ENDITE, ITASK_INNITE, INOTMASTER, INOTSLAVE, postp, dtinv, ittim,&
                                         modul, npoin, itinn
  use def_kintyp,                 only : ip, rp
  use def_chemic,                 only : hrr_int_chm, kfl_hrr_fw_chm, kfl_lookg_chm, kfl_model_chm, kfl_negat_chm, kfl_overs_chm,&
                                         kfl_posit_chm, kfl_post_gp_CMC_chm, kfl_premix_chm, kfl_solve_cond_CMC_chm, kfl_soot_chm,&
                                         kfl_split_CFD_CMC, kfl_split_chm, kfl_spray_chm, kfl_start_CMC_chm, kfl_tab_fw_chm,&
                                         kfl_trans_mxt_spc_CMC_chm, kfl_trans_phs_spc_CMC_chm, kfl_transfer_condField_CMC_chm,&
                                         kfl_ufpv_chm, kfl_under_chm, nZ_CMC_chm, hrr_chm, densi_int_CMC_chm, kfl_multimod_chm
  use mod_ker_updpro,             only : ker_updpro
  use mod_chm_finiteRate,         only : chm_getProp_finiteRate
  use mod_chm_finiteRate,         only : chm_IntegrateSource_finiteRate
  use mod_chm_finiteRate,         only : chm_heatRelease_finiteRate
  use mod_chm_operations_CMC,     only : chm_integrate_chem_source_CMC, &
                                         chm_integrate_flow_var_points_CMC, &
                                         chm_integrate_flow_var_gauss_CMC, &
                                         chm_nodal_projection_CMC, &
                                         chm_limit_Yk_whole_domain_CMC, chm_solve_mixt_frac_diff_CMC, &
                                         chm_solve_mixt_frac_convec_CMC, &
                                         chm_global2local_CMC, chm_calc_temp_CMC, &
                                         chm_heatRelease_field_CMC, &
                                         chm_heatRelease_integral_CMC, &
                                         chm_calc_diff_condVar_mixfraction_CMC, &
                                         chm_update_properties_CMC, chm_integrate_var_CFD_CMC
  use mod_messages,               only : livinf
  use def_kermod,                 only : kfl_chemic_vect
  use mod_communications_global,  only : PAR_SUM
  use mod_arrays,                 only : arrays_number
  implicit none
  integer(ip) :: itask, imixf, ipoin
  real(rp)    :: dt_split, dt

  external    :: chm_clippi
  external    :: chm_heatRelease_lookup
  external    :: chm_cvgunk
  external    :: chm_updunk
  external    :: parari
  external    :: chm_post_scalar_dissipation_rate
  external    :: chm_post_scalar_dist
  external    :: chm_gp_reatab
  external    :: chm_gp_multi_reatab
  external    :: chm_reatab
  external    :: chm_upwmea
  external    :: chm_gp_ANN_eval

  select case(itask)

  case( ITASK_ENDINN )

     !
     ! Cut off over/under shoots
     !
     call chm_clippi()

     !
     ! Calculate total heat release rate in domain
     !
     if (kfl_model_chm == 1 .or. kfl_model_chm == 2_ip ) then
        hrr_int_chm = -666.0_rp
        if (kfl_hrr_fw_chm > 0_ip) then
            call chm_heatRelease_lookup()
        endif
     endif

     !
     ! Compute residual: ||UNKNO(:)-CONCE(:,ICLAS_CHM,1)||
     !
     call chm_cvgunk(ITASK_ENDINN)

     !
     ! Update unknown: CONCE(:,ICLAS_CHM,1)=UNKNO
     !
     call chm_updunk(ITASK_ENDINN)


  case(ITASK_ENDITE)
     !
     !  Compute residual of coupling iteration + update unknown
     !  This step is called after the module has converged. It is the final step before leaving the module
     !
     call livinf(16_ip,' ',itinn(modul))
     call chm_cvgunk(ITASK_ENDITE) ! Residual: ||CONCE(:,:,2)-CONCE(:,:,1)||
     call chm_updunk(ITASK_ENDITE) ! Update:   CONCE(:,:,2) = CONCE(:,:,1)

     call PAR_SUM(kfl_under_chm)
     if( kfl_under_chm > 0 .and. INOTSLAVE ) then
        call livinf(-9_ip,'UNDERSHOOTS HAVE BEEN FOUND= ',kfl_under_chm)
     endif

     call PAR_SUM(kfl_overs_chm)
     if( kfl_overs_chm > 0 .and. INOTSLAVE ) then
        call livinf(-9_ip,'OVERSHOOTS HAVE BEEN FOUND= ',kfl_overs_chm)
     endif


     if (kfl_model_chm == 1 .or. kfl_model_chm == 2_ip ) then
       !
       ! Read flamelet table for gas phase
       !
       if ((kfl_spray_chm == 0 .or. ( kfl_spray_chm /= 0 .and. kfl_premix_chm == 0)) .and. kfl_tab_fw_chm>=0_ip) then

          !
          ! Calculate scalar dissipation rate for UFPV model
          !
          if (kfl_ufpv_chm > 0) then
              if(kfl_chemic_vect /= 1_ip) then
                  call chm_post_scalar_dissipation_rate(24_ip)
              else
                  call chm_post_scalar_dist(24_ip)
              end if
          endif

          !
          ! Table lookup on Gauss points OR on nodes
          !
          if (kfl_lookg_chm > 0) then
              if (kfl_multimod_chm > 0_ip) then
                 call chm_gp_multi_reatab(ITASK_ENDITE)
              else
                 call chm_gp_reatab(ITASK_ENDITE)
              end if 
              call chm_gp_ANN_eval(ITASK_ENDITE)
          else
              call chm_reatab()
          endif

       end if
       call chm_upwmea(ITASK_ENDITE)                         ! wmean(ipoin,1) ==> wmean(ipoin,2)


     elseif (kfl_model_chm == 3) then

       !
       ! Integrate source terms
       !
       if ( kfl_split_chm == 1_ip .or. kfl_split_chm == 3_ip ) then
          dt_split = 1.0_rp / dtinv
          call chm_IntegrateSource_finiteRate(dt_split)
          call chm_heatRelease_finiteRate
       end if
       if ( kfl_split_chm == 2_ip ) then
          dt_split = (1.0_rp / dtinv)/0.5_rp
          call chm_IntegrateSource_finiteRate(dt_split)
          call chm_heatRelease_finiteRate
       end if

       !
       ! Calculate transport properties
       !
       call chm_getProp_finiteRate()
       call chm_upwmea(ITASK_ENDITE)                         ! wmean(ipoin,1) ==> wmean(ipoin,2)


     elseif (kfl_model_chm == 4) then

        if (kfl_solve_cond_CMC_chm == 1) then  ! Conditional variables
           !
           ! Diffusion in mixture fraction
           !

           if (kfl_trans_mxt_spc_CMC_chm == 1_ip) then
              dt = 1.0_rp / dtinv
              call chm_solve_mixt_frac_diff_CMC(dt)
              if (kfl_soot_chm > 0_ip)  call chm_solve_mixt_frac_convec_CMC(dt)
           end if

           if (kfl_trans_phs_spc_CMC_chm == 1_ip .or. kfl_trans_mxt_spc_CMC_chm == 1_ip) then
              ! Update temperature
              do imixf = 2, nZ_CMC_chm-1_ip
                 call chm_global2local_CMC(imixf)
                 call chm_calc_temp_CMC(imixf)
              end do
           end if

           !
           ! Integrate source terms for CMC
           !
           if ( kfl_split_chm > 0_ip  ) then
              dt_split = 1.0_rp / dtinv / dble(kfl_split_chm)
              call chm_integrate_chem_source_CMC(dt_split)
              if( INOTMASTER )   call chm_heatRelease_field_CMC   ! Compute heat release at nodes
           end if

           call chm_limit_Yk_whole_domain_CMC

           ! Update properties with the current composition
           call chm_update_properties_CMC

           ! Do integrations
           if (kfl_post_gp_CMC_chm == 0_ip)   call chm_integrate_flow_var_points_CMC
           call chm_integrate_flow_var_gauss_CMC
           call chm_nodal_projection_CMC

           ! Include density in heat release field (heat release per volume unit)
           do ipoin = 1, npoin
              hrr_chm(ipoin) = hrr_chm(ipoin) * densi_int_CMC_chm(ipoin)
           end do
           call chm_heatRelease_integral_CMC   ! Compute the integral of heat release in the whole domain

           if ( kfl_split_CFD_CMC == 0_ip ) then
              call ker_updpro   ! Update properties with the unconditional values
           else
              if (postp(1) % npp_stepi(arrays_number('MASUN'),0) /= 0) then
                 if( mod(ittim, postp(1) % npp_stepi(arrays_number('MASUN'),0) ) == 0_ip ) then
                    call ker_updpro   ! Update properties with the unconditional values for later writing
                 end if
              end if
           end if

           if (postp(1) % npp_stepi(arrays_number('MASUN'),0) /= 0) then
              if( mod(ittim, postp(1) % npp_stepi(arrays_number('MASUN'),0) ) == 0_ip ) then
                 call chm_calc_diff_condVar_mixfraction_CMC   ! Compute second deriv. in mixt. frac.
              end if
           end if

        else  ! Unconditional mixing variables
           if(kfl_chemic_vect /= 1_ip) then
              call chm_post_scalar_dissipation_rate(24_ip)  ! Xres
              call chm_post_scalar_dissipation_rate(26_ip)  ! Xsgs
           else
              call chm_post_scalar_dist(24_ip)
              call chm_post_scalar_dist(26_ip)
           end if

           if (kfl_transfer_condField_CMC_chm > 0_ip ) then
              call chm_integrate_var_CFD_CMC(6_ip)
              call chm_heatRelease_integral_CMC
           end if

        end if

        if (kfl_start_CMC_chm == 1)    kfl_start_CMC_chm = 0
     endif


  case(ITASK_INNITE)
      !
      ! After Runge-Kutta substeps:
      ! * Clip under and overshoots
      ! * CONCE(:,ICLAS_CHM,1)=UNKNO
      !
     if ( kfl_negat_chm == 1 .or. kfl_posit_chm == 1) &
          call chm_clippi

      call chm_updunk(ITASK_INNITE)

  end select

end subroutine chm_endite
