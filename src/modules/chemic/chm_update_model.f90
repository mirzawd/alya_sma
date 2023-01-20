!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_update_model()

  use def_master,                 only : INOTMASTER, ITASK_BEGSTE, ITASK_ENDITE, kfl_rstar
  use def_kintyp,                 only : ip
  use def_chemic,                 only : kfl_bc_init_method_CMC_chm, kfl_fields_scale_chm, kfl_lookg_chm, kfl_model_chm,&
                                         kfl_premix_chm, kfl_solve_cond_CMC_chm, kfl_spray_chm, kfl_tab_fw_chm, kfl_ufpv_chm,&
                                         kfl_multimod_chm
  use mod_chm_finiteRate,         only : chm_getProp_finiteRate
  use mod_chm_operations_CMC,     only : chm_save_unconditional_fields_bc_CMC, &
                                         chm_bc_type_species_CMC, &
                                         chm_compute_initial_fields_CMC, &
                                         chm_initialization_domain_CMC, &
                                         chm_preprocessing_RK_CMC, &
                                         chm_construct_initialSol_from_chemistry_CMC
  use def_kermod,                 only : kfl_chemic_vect
  use mod_ker_updpro,             only : ker_updpro

#ifdef CANTERA
  use def_chemic,                 only : mechanism_path, nsize_mech_name, coeff_cp_k, W_k, kfl_pfa_chm, Red_spec
#endif

  implicit none
  integer(ip) :: ii

  external    :: chm_scale_Yc
  external    :: chm_gp_reatab
  external    :: chm_post_scalar_dissipation_rate
  external    :: chm_post_scalar_dist
  external    :: chm_reatab
  external    :: chm_fields
  external    :: chm_gp_ANN_eval
  external    :: chm_gp_multi_reatab
#ifdef CANTERA
  external    :: cantera_initialization
  external    :: cantera_coeff
  external    :: cantera_reduced
#endif
  if (kfl_model_chm == 1_ip .or. kfl_model_chm == 2_ip ) then

     !=========================!
     ! FLAMELET or MIXED MODEL !
     !=========================!

     !
     ! Scaling reaction progress Yc from c for initialization
     !
     if ( kfl_premix_chm == 0_ip .and. kfl_fields_scale_chm /= 0_ip ) call chm_scale_Yc()

     !
     ! Read flamelet table for gas phase
     !
     if ((kfl_spray_chm == 0 .or. ( kfl_spray_chm /= 0 .and. kfl_premix_chm == 0)) .and. kfl_tab_fw_chm>=0_ip) then
        if (kfl_lookg_chm > 0) then
           if (kfl_multimod_chm > 0_ip) then
             call chm_gp_multi_reatab(ITASK_BEGSTE)
             call chm_gp_multi_reatab(ITASK_ENDITE)
           else 
              call chm_gp_reatab(ITASK_BEGSTE)
              call chm_gp_reatab(ITASK_ENDITE)
           end if 
           call chm_gp_ANN_eval(ITASK_BEGSTE)
           call chm_gp_ANN_eval(ITASK_ENDITE)
        else
           call chm_reatab()
        endif

        !
        ! Update a few times to get the scalar dissipation right
        !
        if (kfl_ufpv_chm > 0) then
           do ii=1,4
              call ker_updpro()
              if(kfl_chemic_vect /= 1_ip) then
                  call chm_post_scalar_dissipation_rate(24_ip)
                  call chm_post_scalar_dissipation_rate(26_ip)
              else
                  call chm_post_scalar_dist(24_ip)
                  call chm_post_scalar_dist(26_ip)
              end if


              if ( kfl_fields_scale_chm == 1_ip ) call chm_scale_Yc()
              if (kfl_lookg_chm > 0) then
                 if (kfl_multimod_chm > 0_ip) then
                   call chm_gp_multi_reatab(ITASK_BEGSTE)
                   call chm_gp_multi_reatab(ITASK_ENDITE)
                 else 
                    call chm_gp_reatab(ITASK_BEGSTE)
                    call chm_gp_reatab(ITASK_ENDITE)
                 end if 
                 call chm_gp_ANN_eval(ITASK_BEGSTE)
                 call chm_gp_ANN_eval(ITASK_ENDITE)
              else
                 call chm_reatab()
              endif
           enddo
        endif
     end if

  elseif (kfl_model_chm == 3_ip) then

     !=============================!
     ! FINITE RATE CHEMISTRY MODEL !
     !=============================!

#ifdef CANTERA
     !
     ! NASA polinomial coefficients and molecular weights for individual species
     !
     call cantera_initialization(1_ip,mechanism_path,nsize_mech_name)
     call cantera_coeff(mechanism_path,nsize_mech_name,coeff_cp_k, W_k)
     if (kfl_pfa_chm == 1_ip) then
        call cantera_reduced(Red_spec)
     end if
#endif
     call chm_getProp_finiteRate()


  elseif (kfl_model_chm == 4_ip) then
        !==================================!
        ! CONDITIONAL MOMENT CLOSURE MODEL !
        !==================================!

     if( INOTMASTER ) then
        if (kfl_solve_cond_CMC_chm == 1) then ! Conditional species (and enthalpy)
#ifdef CANTERA
           !
           ! NASA polinomial coefficients and molecular weights for individual species
           !
           call cantera_initialization(1_ip,mechanism_path,nsize_mech_name)
           call cantera_coeff(mechanism_path,nsize_mech_name,coeff_cp_k, W_k)
#endif
           call chm_bc_type_species_CMC
           call chm_save_unconditional_fields_bc_CMC
           if (kfl_rstar == 0) then  ! Starting from scratch
              call chm_fields  ! Initialization by fields activated
           end if
           if (kfl_bc_init_method_CMC_chm == 2_ip) then
              call chm_construct_initialSol_from_chemistry_CMC
           end if
           call chm_compute_initial_fields_CMC
           call chm_initialization_domain_CMC
           call chm_preprocessing_RK_CMC

        else ! Unconditional variables
           if(kfl_chemic_vect /= 1_ip) then
               call chm_post_scalar_dissipation_rate(24_ip)
               call chm_post_scalar_dissipation_rate(26_ip)
           else
               call chm_post_scalar_dist(24_ip)
               call chm_post_scalar_dist(26_ip)
           end if
        end if
     end if
  endif

end subroutine chm_update_model
