!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_chm_operations_CMC
!---------------------------------------------
! Module for the CMC subroutines.
! The module agglutinates subroutines in blocks
! with the orientative names:
! a) Initial subroutines
! b) Memory allocation
! c) Boundary conditions
! d) Subroutines for communication
! e) Transfer local <--> global
! f) Terms modelling and assembly
! g) Points transfer
! h) Terms modelling
! i) Diffusion in mixture fraction (equation)
! j) AMC model
! k) Mixture fraction diffusion (term)
! l) Compute time step
! m) Compute variables of interest
! n) Integrations
! o) Variables of interest (others)
! p) Chemical calculation
! q) Clippings
! r) Other functions
! s) Mathematical functions
!---------------------------------------------

  use def_kintyp,              only   : ip,rp
  use def_domain,              only   : ndime,npoin,coord,mgaus,mnode
  use def_chemic,              only   : nclas_chm, nspec_chm
  use mod_memory,              only   : memory_alloca, memory_deallo
  use def_master,              only   : INOTMASTER
  use def_master,              only   : mem_modul,modul
  use mod_communications_global, only : PAR_MAX
  use mod_communications_global, only : PAR_MIN
  use mod_arrays,                only : arrays_number
  
  implicit none

  private

  public :: chm_inert_eq_CMC
  public :: chm_initial_actions_reaphy_CMC
  public :: chm_initial_actions_reanut_CMC
  public :: chm_chemistry_limits_CMC
  public :: chm_construct_vector_chem_integ_CMC
  public :: chm_initial_actions_reaous_CMC
  public :: chm_initial_actions_reabcs_CMC
  public :: chm_compute_initial_fields_CMC
  public :: chm_initialization_domain_CMC
  public :: chm_construct_initialSol_from_chemistry_CMC
  public :: chm_allocate_memory_CMC
  public :: chm_bc_type_species_CMC
  public :: chm_save_unconditional_fields_bc_CMC
  public :: chm_get_cond_fields_bc_Dir_CMC
  public :: chm_find_master_CMC
  public :: chm_transfer_info_CMC_to_CFD
  public :: chm_transfer_info_CFD_to_CMC
  public :: chm_write_info_transfer_CMC
  public :: chm_data_CFD_to_CMC_domain_CMC
  public :: chm_data_CMC_to_CFD_domain_CMC
  public :: chm_local2global_CMC
  public :: chm_global2local_CMC
  public :: chm_element_operations_CMC
  public :: chm_solve_mixt_frac_diff_CMC
  public :: chm_solve_mixt_frac_convec_CMC
  public :: chm_preprocessing_RK_CMC
  public :: chm_AMC_generate_Z_S_vectors_CMC
  public :: chm_AMC_integrals_CMC
  public :: compute_Xnormalized_profile_CMC
  public :: chm_compute_interp_fields_mesh_CMC
  public :: chm_calc_diff_condVar_mixfraction_CMC
  public :: chm_calc_temp_CMC
  public :: chm_update_properties_CMC
  public :: chm_integrate_var_CFD_CMC
  public :: chm_integrate_flow_var_points_CMC
  public :: chm_integrate_flow_var_gauss_CMC
  public :: chm_nodal_projection_CMC
  public :: chm_soot_transfer_source_terms_CMC
  public :: chm_smooth_field_CMC
  public :: chm_integrate_chem_source_CMC
  public :: chm_PDF_calc_CMC
  public :: min_PDF_value_CMC
  public :: chm_updtcc_CMC
  public :: chm_limit_Yk_whole_domain_CMC
  public :: chm_limit_Yk_mixt_fraction_slice_CMC
  public :: clipping_var_from_CFD_CMC
  public :: chm_limit_extremes_PDF_CMC
  public :: chm_heatRelease_field_CMC
  public :: chm_heatRelease_integral_CMC
  public :: chm_values_cvgunk_CMC
contains


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! I N I T I A L   S U B R O U T I N E S !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine chm_inert_eq_CMC

    !-----------------------------------------------------------------------
    !****f* Chemic/mod_chm_operations_CMC/chm_inert_eq_CMC
    ! NAME
    !    chm_inert_eq_CMC
    ! DESCRIPTION
    !    Compute the inert and equilbrium to initialize CMC fields for the
    !    first time.
    ! USES
    ! USED BY
    !    chm_reaphy
    !    chm_initialization_domain_CMC
    !***
    !-----------------------------------------------------------------------

    use def_master,     only :  speci
    use def_chemic,     only :  nZ_CMC_chm, react_scal_bc_CMC_chm, rscal_inert_CMC_chm, &
                                rscal_equil_CMC_chm, Z_CMC_chm, Zs_CMC_chm, &
                                temp_inert_CMC_chm, temp_equil_CMC_chm, &
                                Zstq_CMC_chm, pos_Zstq_CMC_chm, &
                                max_diff_Yk_inert_eq_CMC_chm
    use mod_physics,    only :  physics_H_2_TCp, physics_T_2_HCp
#ifdef CANTERA
    use def_chemic,     only :  mechanism_path, nsize_mech_name, W_k, T_bc_CMC_chm
    use def_master,     only :  prthe
#endif

    implicit none
    integer(ip)              :: imixf, iclas
    real(rp)                 :: temp_chem_equil(nZ_CMC_chm), rscal_chem_equil(nZ_CMC_chm,nspec_chm+1)
    real(rp)                 :: ratio_mf
    character(len=20)        :: names_rscal(nspec_chm+3), number_str
#ifdef CANTERA
    integer(ip)              :: iter, ivalu
    real(rp)                 :: h_mixt, Tini, Taux
    real(rp)                 :: coeff_cp_k_aux(nspec_chm,15) ! Memory still not assigned to coeff_cp_k
    real(rp)                 :: cploc(6,2), aux_cp_lt(8), aux_cp_ht(8)
    real(rp)                 :: aux_cp

    external                 :: cantera_initialization
    external                 :: cantera_coeff
    external                 :: cantera_alya_cp
    external                 :: cantera_get_temperature_from_hy
    external                 :: cantera_equilibrium_from_hy
#endif

#ifdef CANTERA
    call cantera_initialization(1_ip,mechanism_path,nsize_mech_name) ! Create gas object to be used in next functions
    call cantera_coeff(mechanism_path,nsize_mech_name,coeff_cp_k_aux, W_k)

    ! Get enthalpies at boundaries from temperature and composition
    do iter = 1,2
       ! calc_h_cp_from_TY cannot be called since coeff_cp_k has no memory assigned yet
       call cantera_alya_cp(nspec_chm, coeff_cp_k_aux, react_scal_bc_CMC_chm(2_ip:nspec_chm+1_ip,iter), &
                 W_k, aux_cp_lt, aux_cp_ht)
       do ivalu = 1, 6
          cploc(ivalu,1) = aux_cp_lt(ivalu)
          cploc(ivalu,2) = aux_cp_ht(ivalu)
       end do
       call physics_T_2_HCp(T_bc_CMC_chm(iter), cploc, h_mixt, aux_cp)

       react_scal_bc_CMC_chm(1,iter) = h_mixt
    end do
#endif


    ! Get inert profiles
    do iclas = 1,nspec_chm+1_ip
       rscal_inert_CMC_chm(1,iclas)          = react_scal_bc_CMC_chm(iclas,1)
       rscal_inert_CMC_chm(nZ_CMC_chm,iclas) = react_scal_bc_CMC_chm(iclas,2)
    end do

    do imixf = 2,nZ_CMC_chm-1_ip
       ratio_mf = Z_CMC_chm(imixf) / Zs_CMC_chm
       do iclas = 1,nspec_chm+1_ip
          rscal_inert_CMC_chm(imixf,iclas) = react_scal_bc_CMC_chm(iclas,1) + &
               (react_scal_bc_CMC_chm(iclas,2) - react_scal_bc_CMC_chm(iclas,1)) * ratio_mf
       end do
    end do

    ! Get inert temperature from enthalpy and mass fractions
#ifdef CANTERA
    do imixf  = 1,nZ_CMC_chm
       ! Cantera is called to get an initial temperature to start iterating from
       ! Cantera and Alya values are very close but for consistency we take as the good
       ! one the one from Alya
       call cantera_get_temperature_from_hy(prthe(1), rscal_inert_CMC_chm(imixf,1), &
               rscal_inert_CMC_chm(imixf,2_ip:nspec_chm+1_ip), Tini)

       call cantera_alya_cp(nspec_chm, coeff_cp_k_aux, rscal_inert_CMC_chm(imixf,2_ip:nspec_chm+1_ip), &
              W_k, aux_cp_lt, aux_cp_ht)
       do ivalu = 1, 6
          cploc(ivalu,1) = aux_cp_lt(ivalu)
          cploc(ivalu,2) = aux_cp_ht(ivalu)
       end do
       Taux = Tini
       call physics_H_2_TCp(rscal_inert_CMC_chm(imixf,1), cploc, Taux, aux_cp)

       temp_inert_CMC_chm(imixf) = Taux
    end do
#endif

    ! Get equilibrium profiles (only computed to write the corresponding file)
    temp_chem_equil  = temp_inert_CMC_chm
    rscal_chem_equil = rscal_inert_CMC_chm  ! Enthalpy is the same in inert and equilibrium conditions. Species mass fractions
                                               ! are overwritten in the following
#ifdef CANTERA
    do imixf = 2,nZ_CMC_chm-1_ip
       ! Get chemical equilibrium
       call cantera_equilibrium_from_hy(prthe(1), rscal_chem_equil(imixf,1), &
               rscal_chem_equil(imixf,2_ip:nspec_chm+1_ip), Tini)

       ! Cantera and Alya values are very close but for consistency we take as the good
       ! one the one from Alya
       call cantera_alya_cp(nspec_chm, coeff_cp_k_aux, rscal_chem_equil(imixf,2_ip:nspec_chm+1_ip), &
              W_k, aux_cp_lt, aux_cp_ht)
       do ivalu = 1, 6
          cploc(ivalu,1) = aux_cp_lt(ivalu)
          cploc(ivalu,2) = aux_cp_ht(ivalu)
       end do
       Taux = Tini
       call physics_H_2_TCp(rscal_chem_equil(imixf,1), cploc, Taux, aux_cp)

       temp_chem_equil(imixf) = Taux
    end do
#endif


    ! Get Burke-Schumann solution
    temp_equil_CMC_chm  = temp_chem_equil
    rscal_equil_CMC_chm = rscal_chem_equil  ! Enthalpy is the same in inert and equilibrium conditions. Species mass fractions
                                               ! are overwritten in the following

    do imixf = 2,pos_Zstq_CMC_chm-1_ip
       ratio_mf = Z_CMC_chm(imixf) / Zstq_CMC_chm
       do iclas = 1, nspec_chm+1_ip
          rscal_equil_CMC_chm(imixf,iclas) = rscal_equil_CMC_chm(1,iclas) + &
               (rscal_equil_CMC_chm(pos_Zstq_CMC_chm,iclas) - rscal_equil_CMC_chm(1,iclas)) * ratio_mf
       end do
       temp_equil_CMC_chm(imixf) = temp_equil_CMC_chm(1) + &      ! Initialization for later calculation
           (temp_equil_CMC_chm(pos_Zstq_CMC_chm) - temp_equil_CMC_chm(1)) * ratio_mf
    end do

    do imixf = pos_Zstq_CMC_chm+1_ip, nZ_CMC_chm-1_ip
       ratio_mf = (Z_CMC_chm(imixf) - Zstq_CMC_chm) / (Zs_CMC_chm - Zstq_CMC_chm)
       do iclas = 1, nspec_chm+1_ip
          rscal_equil_CMC_chm(imixf,iclas) = rscal_equil_CMC_chm(pos_Zstq_CMC_chm,iclas) + &
               (rscal_equil_CMC_chm(nZ_CMC_chm,iclas) - rscal_equil_CMC_chm(pos_Zstq_CMC_chm,iclas)) * ratio_mf
       end do
       temp_equil_CMC_chm(imixf) = temp_equil_CMC_chm(pos_Zstq_CMC_chm) + &   ! Initialization for later calculation
           (temp_equil_CMC_chm(nZ_CMC_chm) - temp_equil_CMC_chm(pos_Zstq_CMC_chm)) * ratio_mf
    end do

#ifdef CANTERA
    do imixf = 2,nZ_CMC_chm-1_ip   ! Find real temperature for Burke-Schumann
       call cantera_alya_cp(nspec_chm, coeff_cp_k_aux, rscal_equil_CMC_chm(imixf,2_ip:nspec_chm+1_ip), &
              W_k, aux_cp_lt, aux_cp_ht)
       do ivalu = 1, 6
          cploc(ivalu,1) = aux_cp_lt(ivalu)
          cploc(ivalu,2) = aux_cp_ht(ivalu)
       end do
       Taux = temp_equil_CMC_chm(imixf)
       call physics_H_2_TCp(rscal_equil_CMC_chm(imixf,1), cploc, Taux, aux_cp)

       temp_equil_CMC_chm(imixf) = Taux
    end do
#endif

    ! Find maximum values and maximum differences between inert and equilibrium
    do iclas = 1, nspec_chm
       max_diff_Yk_inert_eq_CMC_chm(iclas) = abs(rscal_equil_CMC_chm(pos_Zstq_CMC_chm,iclas+1_ip) - &
          rscal_inert_CMC_chm(pos_Zstq_CMC_chm,iclas+1_ip))

       do imixf = 2,nZ_CMC_chm-1_ip
          max_diff_Yk_inert_eq_CMC_chm(iclas) = max(max_diff_Yk_inert_eq_CMC_chm(iclas), &
             abs(rscal_chem_equil(imixf,iclas+1_ip)-rscal_inert_CMC_chm(imixf,iclas+1_ip)))
       end do
    end do

    !
    ! Write files for inert and equilibrium solutions
    !

    ! Get the names of the reactive scalars in a vector for the header
    names_rscal(1) = 'MIXT_FRAC [-]:1'
    names_rscal(2) = 'TEMPERATURE [K]:2'
    names_rscal(3) = 'ENTHALPY [J/kg]:3'
    do iclas = 1, nspec_chm
       write(number_str,fmt="(i4)") (iclas+3)
       names_rscal(iclas+3) = 'Y' // trim(speci(iclas)%name) // ' [-]:' // number_str
    end do

    ! Write inert mixture file
    open(unit=1, file='inert_solution.log', status="replace", action="write")
    write(unit=1,fmt="(300(a20))") names_rscal
    do imixf = 1,nZ_CMC_chm
       write(unit=1,fmt="(300(e13.6,8x))") Z_CMC_chm(imixf), temp_inert_CMC_chm(imixf), &
            rscal_inert_CMC_chm(imixf,1_ip:nspec_chm+1_ip)
    end do
    close(unit=1)
    print*, '--| ALYA     CHEMIC: INERT SOLUTION WRITTEN FOR CMC MODEL'

    ! Write equilibrium mixture file
    open(unit=1, file='equilibrium_solution.log', status="replace", action="write")
    write(unit=1,fmt="(300(a20))") names_rscal
    do imixf = 1,nZ_CMC_chm
       write(unit=1,fmt="(300(e13.6,8x))") Z_CMC_chm(imixf), temp_chem_equil(imixf), &
            rscal_chem_equil(imixf,1_ip:nspec_chm+1_ip)
    end do
    close(unit=1)
    print*, '--| ALYA     CHEMIC: EQUILIBRIUM SOLUTION WRITTEN FOR CMC MODEL'

    ! Write Burke-Schumann mixture file
    open(unit=1, file='BK_solution.log', status="replace", action="write")
    write(unit=1,fmt="(300(a20))") names_rscal
    do imixf = 1,nZ_CMC_chm
       write(unit=1,fmt="(300(e13.6,8x))") Z_CMC_chm(imixf), temp_equil_CMC_chm(imixf), &
            rscal_equil_CMC_chm(imixf,1_ip:nspec_chm+1_ip)
    end do
    close(unit=1)
    print*, '--| ALYA     CHEMIC: BURKE-SCHUMANN SOLUTION WRITTEN FOR CMC MODEL'

  end subroutine chm_inert_eq_CMC


  subroutine chm_initial_actions_reaphy_CMC
     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_initial_actions_reaphy_CMC
     ! NAME
     !    chm_initial_actions_reaphy_CMC
     ! DESCRIPTION
     !    Do some initial actions for the treatment of the data read in
     !    in chm_reaphy.
     ! USED BY
     !    chm_reaphy
     !***
     !-----------------------------------------------------------------------

     use def_master,             only : momod, modul, kfl_rstar, speci
     use def_kintyp,             only : lg
     use def_chemic,             only : kfl_bc_alpha_CMC_chm, kfl_bc_init_method_CMC_chm, kfl_incl_PDF_trans_CMC_chm,&
                                        kfl_lookg_chm, kfl_solve_cond_CMC_chm, kfl_solve_enth_CMC_chm, kfl_soot_chm,&
                                        kfl_split_CFD_CMC, kfl_trans_mxt_spc_CMC_chm, kfl_trans_phs_spc_CMC_chm,&
                                        kfl_weigh_in_eq_CMC_chm, nS_AMC_CMC_chm, nvar_CMC_chm, nvar_therm_CMC_chm, nZ_AMC_CMC_chm,&
                                        nZ_CMC_chm, pos_Zstq_CMC_chm, react_scal_bc_CMC_chm, size_tncod_CMC_chm, Zs_CMC_chm,&
                                        Zstq_CMC_chm, Z_AMC_CMC_chm, S_AMC_CMC_chm, Xnormalized_prof_CMC_chm, extremes_Z_CMC_chm,&
                                        Z_CMC_chm, PDF_min_CMC_chm, kfl_field_chm, T_bc_CMC_chm, index_N2

     implicit none
     integer(ip)                   :: iclas, index_zbc, line_f
     real(rp)                      :: sum_Yk
     logical(lg)                   :: found

     external                      :: runend
     external                      :: chm_memphy

     if (kfl_solve_cond_CMC_chm == 1) then   ! Conditional species (and enthalpy)

        write(momod(modul) % lun_outpu,*)
        write(momod(modul) % lun_outpu,*) '---------------------------------------------------------'
        write(momod(modul) % lun_outpu,*) ' C M C   P H Y S I C A L   P R O B L E M   S E C T I O N '
        write(momod(modul) % lun_outpu,*) '---------------------------------------------------------'
        write(momod(modul) % lun_outpu,*)

        !!! BORRAR
        !!!write(momod(modul) % lun_outpu,*) 'CHEMIC REAPHY: solving species (and enthalpy) conditional CMC transport equations'

        if (nZ_CMC_chm == 1) then
           call runend('CHEMIC REAPHY: mixture fraction vector cannot contain one single value')
        end if

        if(kfl_solve_enth_CMC_chm == 0) then
           nvar_CMC_chm       = nclas_chm ! CMC only solves species and soot
           nvar_therm_CMC_chm = nspec_chm

        else
           nvar_CMC_chm       = nclas_chm + 1_ip  ! CMC solves species + soot + enthalpy
           nvar_therm_CMC_chm = nspec_chm + 1_ip
        end if

        ! AMC model for scalar dissipation rate
        call chm_memphy(5_ip)
        call chm_AMC_generate_Z_S_vectors_CMC

        write(momod(modul) % lun_outpu,*)
        write(momod(modul) % lun_outpu,*) 'Mixture fraction (', nZ_AMC_CMC_chm, ' entries) for AMC model:'
        do line_f = 1,nZ_AMC_CMC_chm
           write(momod(modul) % lun_outpu,*) line_f, ': ', Z_AMC_CMC_chm(line_f)
        end do
        write(momod(modul) % lun_outpu,*)

        write(momod(modul) % lun_outpu,*)
        write(momod(modul) % lun_outpu,*) 'Segregation factor (', nS_AMC_CMC_chm, ' entries) for AMC model:'
        do line_f = 1,nS_AMC_CMC_chm
           write(momod(modul) % lun_outpu,*) line_f, ': ', S_AMC_CMC_chm(line_f)
        end do
        write(momod(modul) % lun_outpu,*)

        call chm_memphy(6_ip)
        call compute_Xnormalized_profile_CMC   ! Compute scalar dissipation rate normalized profile
        write(momod(modul) % lun_outpu,*)
        write(momod(modul) % lun_outpu,*) 'Normalized scalar dissipation rate profile (', nZ_CMC_chm, ' entries) for AMC model:'
        do line_f = 1,nZ_CMC_chm
           write(momod(modul) % lun_outpu,*) line_f, ': ', Xnormalized_prof_CMC_chm(line_f)
        end do
        write(momod(modul) % lun_outpu,*)
        call chm_AMC_integrals_CMC   ! Compute integrals for AMC model

        call chm_limit_extremes_PDF_CMC
        write(momod(modul) % lun_outpu,*)
        write(momod(modul) % lun_outpu,*) 'Mixture fraction limits for PDF considerations are ', extremes_Z_CMC_chm(1), ' and ',&
            extremes_Z_CMC_chm(2)

        call chm_memphy(8_ip)
        call min_PDF_value_CMC
        write(momod(modul) % lun_outpu,*)
        write(momod(modul) % lun_outpu,*) 'Minimum PDF values considered for each mixture fraction for CMC PDEs are '
        do line_f = 1,nZ_CMC_chm
           write(momod(modul) % lun_outpu,*) Z_CMC_chm(line_f), ': ', PDF_min_CMC_chm(line_f)
        end do
        write(momod(modul) % lun_outpu,*)


        ! Find position for stoichiometric mixture fraction
        if (Zstq_CMC_chm <= 0.0_rp)   call runend('CHEMIC REAPHY: Stoichiometric mixture fraction cannot be less than 0')
        if (Zstq_CMC_chm >= Zs_CMC_chm)   call runend('CHEMIC REAPHY: Stoichiometric mixture fraction cannot be higher than&
            & saturation mixture fraction')
        pos_Zstq_CMC_chm = 2_ip
        do while (Z_CMC_chm(pos_Zstq_CMC_chm) < Zstq_CMC_chm)
           pos_Zstq_CMC_chm = pos_Zstq_CMC_chm + 1_ip
        end do
        if ((Z_CMC_chm(pos_Zstq_CMC_chm+1)-Zstq_CMC_chm) < (Zstq_CMC_chm-Z_CMC_chm(pos_Zstq_CMC_chm))) then
           pos_Zstq_CMC_chm = pos_Zstq_CMC_chm + 1_ip
        end if
        Zstq_CMC_chm = Z_CMC_chm(pos_Zstq_CMC_chm)
        write(momod(modul) % lun_outpu,*)
        write(momod(modul) % lun_outpu,*) 'Stoichiometric mixture fraction (for calculations): ', Zstq_CMC_chm
        write(momod(modul) % lun_outpu,*)


        call assign_val_clipZvar

        if (kfl_trans_phs_spc_CMC_chm == 1_ip) then
           write(momod(modul) % lun_outpu,*)
           write(momod(modul) % lun_outpu,*) 'Transport in physical space included in CMC equations'
           write(momod(modul) % lun_outpu,*)
        else
           write(momod(modul) % lun_outpu,*)
           write(momod(modul) % lun_outpu,*) 'Transport in physical space not inlcuded in CMC equations'
           write(momod(modul) % lun_outpu,*)
        end if

        if (kfl_trans_mxt_spc_CMC_chm == 1_ip) then
           write(momod(modul) % lun_outpu,*)
           write(momod(modul) % lun_outpu,*) 'Transport in mixture fraction space included in CMC equations'
           write(momod(modul) % lun_outpu,*)
        else
           write(momod(modul) % lun_outpu,*)
           write(momod(modul) % lun_outpu,*) 'Transport in mixture fraction space not inlcuded in CMC equations'
           write(momod(modul) % lun_outpu,*)
        end if


        if (kfl_incl_PDF_trans_CMC_chm == 1_ip) then
           write(momod(modul) % lun_outpu,*)
           write(momod(modul) % lun_outpu,*) 'PDF included in CMC transport equations'
           write(momod(modul) % lun_outpu,*)
        else
           write(momod(modul) % lun_outpu,*)
           write(momod(modul) % lun_outpu,*) 'PDF not included in CMC transport equations'
           write(momod(modul) % lun_outpu,*)
        end if

        ! Check we have all the information to compute initial fields and read BCs
        if ( kfl_rstar == 0 .and. kfl_weigh_in_eq_CMC_chm == 1 ) then
           if (kfl_field_chm(3) == 0)  call runend('CHEMIC REAPHY: alpha field not provided')
        end if

        ! If enthalpy is solved kfl_field_chm(2) is needed for BCs and may be used for initial cond.
        if (kfl_field_chm(2) == 0 .and. kfl_solve_enth_CMC_chm /= 0) &
           call runend('CHEMIC REAPHY: temperature not provided but required to solve enthalpy')


        ! Compute inert and equilibrium and do some additional checks
        do index_zbc = 1,2
           ! Do some checks of the mass fractions
           react_scal_bc_CMC_chm = abs(react_scal_bc_CMC_chm)
           sum_Yk = 0.0_rp
           do iclas = 2, nspec_chm+1
              sum_Yk = sum_Yk + react_scal_bc_CMC_chm(iclas,index_zbc)
           end do
           if (sum_Yk /= 0.0_rp)  then  ! If Yk do not sum up 1.0 => normalize
              react_scal_bc_CMC_chm(2:nspec_chm+1,index_zbc) = react_scal_bc_CMC_chm(2:nspec_chm+1,index_zbc) / sum_Yk
           else
              call runend('CHEMIC REAPHY: mass fractions at the mixture fractions boundaries not provided')
           end if

           write(momod(modul) % lun_outpu,*) 'Reactive scalar values at the boundary conditions'
           if (index_zbc == 1) then
              write(momod(modul) % lun_outpu,*) 'OXIDIZER:'
           else
              write(momod(modul) % lun_outpu,*) 'FUEL:'
           end if
           write(momod(modul) % lun_outpu,*) 'Temperature:', T_bc_CMC_chm(index_zbc), ' K'
           write(momod(modul) % lun_outpu,*) 'Mass fractions'
           do iclas = 1, nspec_chm
              write(momod(modul) % lun_outpu,*) speci(iclas)%name, ': ', react_scal_bc_CMC_chm(iclas+1,index_zbc)
           end do
           write(momod(modul) % lun_outpu,*)
           write(momod(modul) % lun_outpu,*)
        end do

        call chm_memphy(7_ip)
        call chm_inert_eq_CMC

        if (kfl_solve_enth_CMC_chm /= 0) then
           iclas = 1
           found = .false.
           do while( (iclas <= nspec_chm) .and. (.not. found))
              if ('N2' == speci(iclas)%name) then
                 index_N2 = iclas
                 found = .true.
              end if
              iclas = iclas + 1
           end do
           if( .not. found )   index_N2 = 0
        end if

        if ( kfl_bc_init_method_CMC_chm == 1_ip ) then
           write(momod(modul) % lun_outpu,*)
           write(momod(modul) % lun_outpu,*) 'Initial and boundary conditions generated with linear relationships'
           write(momod(modul) % lun_outpu,*)
        !!!elseif ( kfl_bc_init_method_CMC_chm == 2_ip ) then
        !!!   inv_exp_alpha_grid_CMC_chm = 1.0_rp / exp_alpha_grid_CMC_chm
        !!!   call chm_memphy(11_ip)
        !!!   write(momod(modul) % lun_outpu,*)
        !!!   write(momod(modul) % lun_outpu,*) 'Initial and boundary conditions generated from homogeneous reactor evolutions'
        !!!   write(momod(modul) % lun_outpu,*)
        end if

        ! Definition of variables for boundary conditions
        if (kfl_bc_alpha_CMC_chm == -1_ip)  then
           ! For additional comments see boundary conditions section
           call runend('CHEMIC REAPHY: the type of boundary conditions (species or alpha) has to be specified')
        else if (kfl_bc_alpha_CMC_chm == 0_ip) then
           if (kfl_soot_chm > 0)  call runend('CHEMIC REAPHY: not compatible unconditional species BCs and soot model')
           if (kfl_field_chm(2) > 0) then
              ! Boundary conditions for temperature are written in the file even if they are not used
              size_tncod_CMC_chm = nspec_chm+1_ip
           else
              ! Nothing related to temperature is written in the file
              size_tncod_CMC_chm = nspec_chm
           end if
        else
           if (kfl_field_chm(2) > 0) then
              ! Boundary conditions for temperature are written in the file even if they are not used
              size_tncod_CMC_chm = 2_ip
           else
              ! Nothing related to temperature is written in the file
              size_tncod_CMC_chm = 1_ip
           end if
        end if


     else  ! Unconditional mixing variables

        !!! BORRAR
        !!!write(momod(modul) % lun_outpu,*) 'CHEMIC REAPHY: solving mixture fraction, its variance, etc. unconditional transport equations'

        !!! BORRAR
        !!!nclas_chm = 4 ! Even if in this case only mixture fraction and its
        !!!              ! variance transport equations are solved, progress
        !!!              ! variable and its variance are considered but they
        !!!              ! are meaningless in this context

        if (kfl_split_CFD_CMC==0) call runend('CHEMIC REAPHY: no splitting and solving unconditional variables is not valid in&
            & CMC model')
     end if

     if (kfl_lookg_chm == 0) then
        call runend('CHEMIC REAPHY: CMC model not compatible with calculation at nodes')
     end if

     if (kfl_split_CFD_CMC==0) then
        write(momod(modul) % lun_outpu,*) 'CHEMIC REAPHY: whole simulation in one single execution'
     else
        write(momod(modul) % lun_outpu,*) 'CHEMIC REAPHY: split simulation coupled by two executions'
     end if

  end subroutine chm_initial_actions_reaphy_CMC


  subroutine chm_initial_actions_reanut_CMC
     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_initial_actions_reanut_CMC
     ! NAME
     !    chm_initial_actions_reanut_CMC
     ! DESCRIPTION
     !    Do some checks in the variables defined in NUMERICAL TREATMENT.
     ! USED BY
     !    chm_reanut
     !***
     !-----------------------------------------------------------------------
     use def_chemic,       only :  kfl_mesh_interp_CMC_chm, kfl_split_CFD_CMC, &
                                   kfl_avg_cond_CMC_chm, kfl_solve_cond_CMC_chm, &
                                   nZ_chm_int_CMC_chm, posZ_chem_integr_CMC_chm, &
                                   chem_int_iZ_CMC_chm, Z_CMC_chm, kfl_transfer_condField_CMC_chm
     use def_coupli,       only :  ON_FLOATING_POINTS, coupling_type, mcoup
     use def_master,       only :  momod

     implicit none
     integer(ip)                :: icoup, ii

     external                   :: runend

     write(momod(modul) % lun_outpu,*)
     write(momod(modul) % lun_outpu,*) '---------------------------------------------------------------'
     write(momod(modul) % lun_outpu,*) ' C M C   N U M E R I C A L   T R E A T M E N T   S E C T I O N '
     write(momod(modul) % lun_outpu,*) '---------------------------------------------------------------'
     write(momod(modul) % lun_outpu,*)

     if (kfl_mesh_interp_CMC_chm == 1_ip .and. kfl_split_CFD_CMC == 0_ip) then
        call runend('CHEMIC REANUT: not possible to interpolate between meshes and run CFD and CMC in the same instance')
     end if
     do icoup = 1, mcoup
        if (coupling_type(icoup) % variable == 'CF2CM') then
           if (coupling_type(icoup) % where_type == ON_FLOATING_POINTS .and. &
               kfl_mesh_interp_CMC_chm == 0_ip) then
              call runend('CHEMIC REANUT: not possible to define coupling with FLOAT but not doing interpolations')
           else if (coupling_type(icoup) % where_type /= ON_FLOATING_POINTS .and. &
               kfl_mesh_interp_CMC_chm == 1_ip) then
              call runend('CHEMIC REANUT: not possible to not define coupling with FLOAT but doing interpolations')
           end if
        end if
     end do

     write(momod(modul) % lun_outpu,*)''
     if (kfl_mesh_interp_CMC_chm == 0_ip) then
        write(momod(modul) % lun_outpu,*) 'Same mesh for CMC and CFD'
     else
        write(momod(modul) % lun_outpu,*) 'Different meshes for CMC and CFD. Fields will be interpolated'
        if (kfl_solve_cond_CMC_chm == 1_ip) then
           if (kfl_avg_cond_CMC_chm == 1_ip) then
              write(momod(modul) % lun_outpu,*) 'Volumetric averages computed from the conditional fields weighing with the PDF'
           else
              write(momod(modul) % lun_outpu,*) 'Volumetric averages computed from the unconditional fields'
           end if
        end if
     end if
     write(momod(modul) % lun_outpu,*)''


     ! Compute between which positions do chemical integration

     write(momod(modul) % lun_outpu,*)
     if (nZ_chm_int_CMC_chm == posZ_chem_integr_CMC_chm(2) - posZ_chem_integr_CMC_chm(1) + 1_ip) then
        write(momod(modul) % lun_outpu,*) 'Chemistry integrated between (if done)'
        write(momod(modul) % lun_outpu,*) 'Zmin: ', chem_int_iZ_CMC_chm(1), ' Z = ', &
           Z_CMC_chm(chem_int_iZ_CMC_chm(1))
        write(momod(modul) % lun_outpu,*) 'Zmax: ', chem_int_iZ_CMC_chm(nZ_chm_int_CMC_chm), ' Z = ', &
           Z_CMC_chm(chem_int_iZ_CMC_chm(nZ_chm_int_CMC_chm))
     else
        write(momod(modul) % lun_outpu,*) 'Chemistry integrated at mixtures (if done)'
        do ii = 1, nZ_chm_int_CMC_chm
           write(momod(modul) % lun_outpu,*) ii, ': ', chem_int_iZ_CMC_chm(ii), ' Z = ', &
              Z_CMC_chm(chem_int_iZ_CMC_chm(ii))
        end do

     end if
     write(momod(modul) % lun_outpu,*)

     if (kfl_transfer_condField_CMC_chm == 1_ip .and. kfl_mesh_interp_CMC_chm == 0_ip) call runend('CHEMIC REANUT: transfer&
        & activated but meshes are coincident')

     if (kfl_transfer_condField_CMC_chm == 1_ip .and. kfl_split_CFD_CMC == 0_ip)  call runend('CHEMIC REANUT: transfer activated&
        & but only one single execution')

  end subroutine chm_initial_actions_reanut_CMC


  subroutine chm_chemistry_limits_CMC
     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_chemistry_limits_CMC
     ! NAME
     !    chm_chemistry_limits_CMC
     ! DESCRIPTION
     !    Do some checks for mixtures to be integrated
     ! USED BY
     !    chm_reanut
     !***
     !-----------------------------------------------------------------------
     use def_chemic,       only:  extr_Z_chem_integr_CMC_chm, posZ_chem_integr_CMC_chm, &
                                  Z_CMC_chm, Zs_CMC_chm, nZ_CMC_chm

     implicit none
     integer(ip)                :: posZ

     external                   :: runend

     ! Compute between which positions to do chemical integration
     if (extr_Z_chem_integr_CMC_chm(1) > extr_Z_chem_integr_CMC_chm(2)) call runend('CHEMIC REAPHY: Zmin > Zmax for chemical&
         & integration')

     ! Minimum value
     if (extr_Z_chem_integr_CMC_chm(1) == 0.0_rp) then
        ! 1 corresponds to Z=0 but Z=0 is not integrated -> integrate from the first solved mixture fraction
        posZ_chem_integr_CMC_chm(1) = 2
     else
        posZ = 2
        do while (Z_CMC_chm(posZ) < extr_Z_chem_integr_CMC_chm(1))
           posZ = posZ + 1
        end do
        posZ_chem_integr_CMC_chm(1) = posZ
     end if

     ! Maximum value
     if (extr_Z_chem_integr_CMC_chm(2) >= Zs_CMC_chm) then
        ! nZ_CMC_chm corresponds to Z=Zs but Z=Zs is not integrated -> integrate till the last solved mixture fraction
        posZ_chem_integr_CMC_chm(2) = nZ_CMC_chm - 1
     else
        posZ = nZ_CMC_chm - 1
        do while (extr_Z_chem_integr_CMC_chm(2) < Z_CMC_chm(posZ))
           posZ = posZ - 1
        end do
        posZ_chem_integr_CMC_chm(2) = posZ
     end if

  end subroutine chm_chemistry_limits_CMC


  subroutine chm_construct_vector_chem_integ_CMC(aux_nZ, aux_Z, kfl_chm_int)
     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_construct_vector_chem_integ_CMC
     ! NAME
     !    chm_construct_vector_chem_integ_CMC
     ! DESCRIPTION
     !    Create vector of mixture fractions to be chemically integrated.
     ! USED BY
     !    chm_reanut
     !***
     !-----------------------------------------------------------------------

     use def_chemic,       only : nZ_chm_int_CMC_chm, posZ_chem_integr_CMC_chm, &
                                  chem_int_iZ_CMC_chm, nZ_no_chm_int_CMC_chm
     use def_kintyp,       only : lg

     implicit none
     integer(ip), intent(in)    :: aux_nZ, aux_Z(aux_nZ), kfl_chm_int
     integer(ip)                :: ii, jj, kk
     logical(lg)                :: found

     external                   :: chm_memnut
     external                   :: runend

     ! Check if vector ordered
     do jj = 1, nZ_chm_int_CMC_chm - 1_ip
        if (aux_Z(jj) >= aux_Z(jj+1_ip))  call runend('CHEMIC REANUT: vector for mixtures to be chemically integrated not&
            & monotically increasing')
     end do

     !
     ! Construct mixture fraction vector for chemical integration
     !
     if (kfl_chm_int == 1_ip) then

        ! The given positions correspond to integrable mixtures
        nZ_chm_int_CMC_chm = aux_nZ

        ! Count --> look for the intersection [extremes_Z_CMC_chm(1), extremes_Z_CMC_chm(2)] and
        ! mixture fraction positions to keep for integration
        do jj = 1, nZ_chm_int_CMC_chm
           if (aux_Z(jj) < posZ_chem_integr_CMC_chm(1) .or. &
              aux_Z(jj) > posZ_chem_integr_CMC_chm(2))  then
              nZ_chm_int_CMC_chm = nZ_chm_int_CMC_chm - 1_ip
           end if
        end do
        call chm_memnut(1_ip)

        ! Assign
        ii = 1_ip
        do jj = 1, aux_nZ
           if (aux_Z(jj) >= posZ_chem_integr_CMC_chm(1) .and. &
              aux_Z(jj) <= posZ_chem_integr_CMC_chm(2))  then
              chem_int_iZ_CMC_chm(ii) = aux_Z(jj)
              ii = ii + 1_ip
           end if
        end do

        ! Find number of mixtures to not do integration
        nZ_no_chm_int_CMC_chm = posZ_chem_integr_CMC_chm(2) - posZ_chem_integr_CMC_chm(1) &
           + 1_ip - nZ_chm_int_CMC_chm

     else

        ! The given positions do not correspond to integrable mixtures
        nZ_chm_int_CMC_chm = posZ_chem_integr_CMC_chm(2) - posZ_chem_integr_CMC_chm(1) &
           + 1_ip - aux_nZ

        ! Count --> look for the intersection [extremes_Z_CMC_chm(1), extremes_Z_CMC_chm(2)] and
        ! mixture fraction positions to keep for integration
        do jj = 1, aux_nZ
           if (aux_Z(jj) < posZ_chem_integr_CMC_chm(1) .or. &
              aux_Z(jj) > posZ_chem_integr_CMC_chm(2))  then
              nZ_chm_int_CMC_chm = nZ_chm_int_CMC_chm + 1_ip
           end if
        end do
        call chm_memnut(1_ip)

        ! Assign
        kk = 1
        do jj = posZ_chem_integr_CMC_chm(1), posZ_chem_integr_CMC_chm(2)
           ii = 1_ip
           found = .false.
           do while ((.not. found) .and. ii <= aux_nZ)
              if (jj == aux_Z(ii))  found = .true.
              ii = ii + 1_ip
           end do
           if (.not. found)  then
              chem_int_iZ_CMC_chm(kk) = jj
              kk = kk + 1_ip
           end if
        end do

        ! Find number of mixtures to not do integration
        nZ_no_chm_int_CMC_chm = posZ_chem_integr_CMC_chm(2) - posZ_chem_integr_CMC_chm(1) &
           + 1_ip - nZ_chm_int_CMC_chm

     end if

  end subroutine chm_construct_vector_chem_integ_CMC


  subroutine chm_initial_actions_reaous_CMC
     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_initial_actions_reaous_CMC
     ! NAME
     !    chm_initial_actions_reaous_CMC
     ! DESCRIPTION
     !    Do some checks in the variables defined in OUTPUT & POSTPROCESS.
     ! USED BY
     !    chm_reaous
     !***
     !-----------------------------------------------------------------------
     use def_chemic,       only :  nspec_uncond_write_CMC_chm, nspec_cond_write_CMC_chm, &
                                   nZ_write_CMC_chm, write_uncond_spec_CMC_chm, &
                                   write_cond_spec_CMC_chm, write_iZ_spec_CMC_chm, &
                                   kfl_post_gp_CMC_chm, kfl_dt_calc_CMC_chm, Z_CMC_chm, &
                                   nZ_CMC_chm, kfl_soot_chm
     use def_master,       only :  momod, speci

     implicit none
     integer(ip)                :: ii


     write(momod(modul) % lun_outpu,*)
     write(momod(modul) % lun_outpu,*) '-------------------------------------------------------'
     write(momod(modul) % lun_outpu,*) ' C M C   P O S T - P R O C E S S I N G   S E C T I O N '
     write(momod(modul) % lun_outpu,*) '-------------------------------------------------------'
     write(momod(modul) % lun_outpu,*)

     ! Unconditional species mass fractions
     write(momod(modul) % lun_outpu,*)''
     write(momod(modul) % lun_outpu,*)''
     if (nspec_uncond_write_CMC_chm < nclas_chm) then
        if ( kfl_soot_chm /= 0 ) then
           write(momod(modul) % lun_outpu,*) 'Writing ', nspec_uncond_write_CMC_chm, ' unconditional species and soot section mass&
               & fractions'
        else
           write(momod(modul) % lun_outpu,*) 'Writing ', nspec_uncond_write_CMC_chm, ' unconditional species mass fractions'
        end if

        do ii = 1, nspec_uncond_write_CMC_chm
           if (write_uncond_spec_CMC_chm(ii) <= nspec_chm) then
              write(momod(modul) % lun_outpu,*) 'Species: ', write_uncond_spec_CMC_chm(ii), ' - ',&
                  speci(write_uncond_spec_CMC_chm(ii))%name
           else
              write(momod(modul) % lun_outpu,*) 'Soot section: ', write_uncond_spec_CMC_chm(ii) - nspec_chm, ' class: ',&
                  write_uncond_spec_CMC_chm(ii)
           end if
        end do

     else
        if (kfl_soot_chm /= 0) then
           write(momod(modul) % lun_outpu,*)'Post-processing all the unconditional species and soot section mass fractions'
        else
           write(momod(modul) % lun_outpu,*)'Post-processing all the unconditional species mass fractions'
        end if
     end if


     ! Conditional species mass fractions --> species (and soot)
     write(momod(modul) % lun_outpu,*)''
     write(momod(modul) % lun_outpu,*)''
     if (nspec_cond_write_CMC_chm < nclas_chm) then
        if ( kfl_soot_chm /= 0 ) then
           write(momod(modul) % lun_outpu,*)'Writing ', nspec_cond_write_CMC_chm, ' species and soot section for conditional mass&
               & fractions'
        else
           write(momod(modul) % lun_outpu,*)'Writing ', nspec_cond_write_CMC_chm, ' species for conditional mass fractions'
        end if

        do ii = 1, nspec_cond_write_CMC_chm
           if (write_cond_spec_CMC_chm(ii) <= nspec_chm) then
              write(momod(modul) % lun_outpu,*) 'Species: ', write_cond_spec_CMC_chm(ii), ' - ',&
                  speci(write_cond_spec_CMC_chm(ii))%name
           else
              write(momod(modul) % lun_outpu,*) 'Soot section: ', write_cond_spec_CMC_chm(ii) - nspec_chm, ' class: ',&
                  write_uncond_spec_CMC_chm(ii)
           end if

        end do
     else
        if ( kfl_soot_chm /= 0 ) then
           write(momod(modul) % lun_outpu,*)'Post-processing all the species and soot sections for conditional mass fractions'
        else
           write(momod(modul) % lun_outpu,*)'Post-processing all the species for conditional mass fractions'
        end if
     end if


     ! Conditional species mass fractions --> mixture fractions
     write(momod(modul) % lun_outpu,*)''
     write(momod(modul) % lun_outpu,*)''
     if (nZ_write_CMC_chm < nZ_CMC_chm) then
        write(momod(modul) % lun_outpu,*)'Writing ', nZ_write_CMC_chm, ' mixture fractions for conditional species mass fractions'
        write(momod(modul) % lun_outpu,*) 'Position: mixture fraction'
        do ii = 1, nZ_write_CMC_chm
           write(momod(modul) % lun_outpu,*) write_iZ_spec_CMC_chm(ii), ': ', Z_CMC_chm(write_iZ_spec_CMC_chm(ii))
        end do
     else
        write(momod(modul) % lun_outpu,*)'Post-processing all the mixture fractions for conditional mass fractions'
     end if

     write(momod(modul) % lun_outpu,*)''
     write(momod(modul) % lun_outpu,*)''

     write(momod(modul) % lun_outpu,*)''
     if (kfl_post_gp_CMC_chm == 1_ip) then
        write(momod(modul) % lun_outpu,*) 'Postprocessing variables on Gauss points'
     else
        write(momod(modul) % lun_outpu,*) 'Postprocessing variables on nodes'
     end if
     write(momod(modul) % lun_outpu,*)''

     write(momod(modul) % lun_outpu,*)''
     if (kfl_dt_calc_CMC_chm == 1_ip) then
        write(momod(modul) % lun_outpu,*) 'Critical time step computed'
     else
        write(momod(modul) % lun_outpu,*) 'Critical time step not computed'
     end if
     write(momod(modul) % lun_outpu,*)''

  end subroutine chm_initial_actions_reaous_CMC


  subroutine chm_initial_actions_reabcs_CMC
     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_initial_actions_reabcs_CMC
     ! NAME
     !    chm_initial_actions_reabcs_CMC
     ! DESCRIPTION
     !    Do some checks in the variables defined in BOUNDARY_CONDITIONS.
     ! USED BY
     !    chm_reabcs
     !***
     !-----------------------------------------------------------------------
     use def_chemic,       only :  kfl_bc_alpha_CMC_chm, size_tncod_CMC_chm, kfl_soot_chm
     use def_master,       only :  momod

     implicit none


     write(momod(modul) % lun_outpu,*)
     write(momod(modul) % lun_outpu,*) '---------------------------------------------------------------'
     write(momod(modul) % lun_outpu,*) ' C M C   B O U N D A R Y   C O N D I T I O N S   S E C T I O N '
     write(momod(modul) % lun_outpu,*) '---------------------------------------------------------------'
     write(momod(modul) % lun_outpu,*)


     write(momod(modul) % lun_outpu,*)''
     write(momod(modul) % lun_outpu,*)''
     if (kfl_bc_alpha_CMC_chm == 1_ip) then
        write(momod(modul) % lun_outpu,*) 'Boundary conditions provided through alpha values'
     else
        write(momod(modul) % lun_outpu,*) 'Boundary conditions provided through unconditional mass fractions for species'
     end if
     write(momod(modul) % lun_outpu,*)''
     write(momod(modul) % lun_outpu,*) 'The number of fields for boundary conditions is ', size_tncod_CMC_chm
     write(momod(modul) % lun_outpu,*)''
     write(momod(modul) % lun_outpu,*)''

     if ( kfl_soot_chm /= 0 ) then
        write(momod(modul) % lun_outpu,*)''
        write(momod(modul) % lun_outpu,*)'-------------'
        write(momod(modul) % lun_outpu,*)'C A U T I O N'
        write(momod(modul) % lun_outpu,*)'-------------'
        write(momod(modul) % lun_outpu,*) 'Setting initial conditions and Dirichlet boundary conditions for conditional soot&
            & section mass fractions to 0'
        write(momod(modul) % lun_outpu,*)''
     end if

  end subroutine chm_initial_actions_reabcs_CMC


  subroutine chm_compute_initial_fields_CMC
     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_compute_initial_fields_CMC
     ! NAME
     !    chm_compute_initial_fields_CMC
     ! DESCRIPTION
     !    CMC fields initialization for bvess_CMC_chm
     ! USED BY
     !    chm_iniunk
     !***
     !-----------------------------------------------------------------------
     use def_master,      only : kfl_rstar
     use def_chemic,      only : nZ_CMC_chm, rscal_inert_CMC_chm, rscal_equil_CMC_chm, &
                                 kfl_weigh_in_eq_CMC_chm, bvess_chm, bvess_CMC_chm, &
                                 kfl_solve_enth_CMC_chm, nvar_CMC_chm, Zavg_CMC_chm, &
                                 Zvar_CMC_chm, Yk_CMC_chm, enthalp_CMC_chm, &
                                 kfl_bc_init_method_CMC_chm

     implicit none
     integer(ip)              :: ipoin, imixf, iclas
     real(rp)                 :: aux_dif
     real(rp)                 :: Yk_poin(nZ_CMC_chm,nspec_chm)

     external                 :: runend

     ! Initialization of CMC domain based on a coefficient field and inert and equilibrium solutions

     if ( kfl_rstar == 0_ip) then ! Start from scratch
        if (kfl_weigh_in_eq_CMC_chm == 1_ip) then  ! Initial fields defined by a coefficient alpha
           ! At this point the only field that is filled in bvess_chm is the first
           ! one with the alpha values
           if (kfl_bc_init_method_CMC_chm == 1_ip) then  ! Linear relationships

              ! Fill the inner nodes of bvess_CMC_chm matrix (physical + mixt. frac. spaces)
              do ipoin = 1,npoin
                 do iclas = 1, nspec_chm  ! Mass fractions for soot set to 0
                    ! Mixture fraction boundaries
                    bvess_CMC_chm(1,ipoin,iclas)          = rscal_inert_CMC_chm(1,iclas+1_ip)
                    bvess_CMC_chm(nZ_CMC_chm,ipoin,iclas) = rscal_inert_CMC_chm(nZ_CMC_chm,iclas+1_ip)

                    ! Intermediate mixture fractions
                    do imixf = 2,nZ_CMC_chm-1_ip
                       aux_dif = rscal_equil_CMC_chm(imixf,iclas+1_ip) - rscal_inert_CMC_chm(imixf,iclas+1_ip)
                       bvess_CMC_chm(imixf,ipoin,iclas) = rscal_inert_CMC_chm(imixf,iclas+1_ip) + &
                           bvess_chm(nvar_CMC_chm+1_ip,ipoin) * aux_dif
                    end do
                 end do
              end do

           !!!else if (kfl_bc_init_method_CMC_chm == 2_ip) then  ! Homogeneous reactors

           !!!   ! Fill the inner nodes of bvess_CMC_chm matrix (physical + mixt. frac. spaces)
           !!!   do ipoin = 1,npoin
           !!!      aux = bvess_chm(nvar_CMC_chm+1_ip,ipoin)**inv_exp_alpha_grid_CMC_chm
           !!!      pos = 1_ip + floor(aux*(n_alpha_grid_CMC_chm-1_ip))
           !!!      if (pos == n_alpha_grid_CMC_chm) then  ! bvess_chm(nvar_CMC_chm+1,ipoin)=1.0
           !!!         pos = pos - 1_ip
           !!!      end if
           !!!      ratio = (bvess_chm(nvar_CMC_chm+1_ip,ipoin) - alpha_grid_CMC_chm(pos)) / (alpha_grid_CMC_chm(pos+1_ip)&
           !!!          - alpha_grid_CMC_chm(pos))

           !!!      do iclas = 1,nclas_chm
           !!!         ! Mixture fraction boundaries
           !!!         bvess_CMC_chm(1,ipoin,iclas)          = rscal_inert_CMC_chm(1,iclas+1_ip)
           !!!         bvess_CMC_chm(nZ_CMC_chm,ipoin,iclas) = rscal_inert_CMC_chm(nZ_CMC_chm,iclas+1_ip)

           !!!         ! Intermediate mixture fractions
           !!!         do imixf = 2,nZ_CMC_chm-1_ip
           !!!            aux_dif = homog_reactors_CMC_chm(imixf,pos+1_ip,iclas) - homog_reactors_CMC_chm(imixf,pos,iclas)
           !!!            bvess_CMC_chm(imixf,ipoin,iclas) = homog_reactors_CMC_chm(imixf,pos,iclas) + &
           !!!                ratio * aux_dif
           !!!         end do
           !!!      end do
           !!!   end do

           !!!   call chm_initialize_alpha_CMC  ! Initialize alpha matrix

           end if

           if (kfl_solve_enth_CMC_chm /= 0_ip) then
              do ipoin = 1,npoin
                 do imixf = 1, nZ_CMC_chm
                    do iclas = 1, nspec_chm
                       Yk_poin(imixf,iclas) = bvess_CMC_chm(imixf,ipoin,iclas)
                    end do
                 end do

                 call chm_get_cond_enthalpy_average_CMC(bvess_chm(nvar_CMC_chm,ipoin), &
                         Yk_poin(1:nZ_CMC_chm,:), Zavg_CMC_chm(ipoin), &
                         Zvar_CMC_chm(ipoin), bvess_chm(nvar_CMC_chm+1_ip,ipoin), &
                         bvess_CMC_chm(1:nZ_CMC_chm,ipoin,nvar_CMC_chm))
              end do
           end if

        else
           call runend('CHEMIC: Reading initial fields when starting from scratch not implemented yet')
        end if

     else ! Restart

        do imixf = 1,nZ_CMC_chm
           do ipoin = 1,npoin
              do iclas = 1,nclas_chm
                 bvess_CMC_chm(imixf,ipoin,iclas) = Yk_CMC_chm(imixf,ipoin,iclas)
              end do
           end do
        end do

        if (kfl_solve_enth_CMC_chm /= 0_ip) then
           do imixf = 1,nZ_CMC_chm
              do ipoin = 1,npoin
                 bvess_CMC_chm(imixf,ipoin,nvar_CMC_chm) = enthalp_CMC_chm(imixf,ipoin)
              end do
           end do
        else
           do imixf = 1,nZ_CMC_chm
              enthalp_CMC_chm(imixf,1) = rscal_inert_CMC_chm(imixf,1)
           end do
        end if
     end if

  end subroutine chm_compute_initial_fields_CMC


  subroutine chm_initialization_domain_CMC
     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_initialization_domain_CMC
     ! NAME
     !    chm_initialization_domain_CMC
     ! DESCRIPTION
     !    CMC fields initialization
     ! USED BY
     !    chm_iniunk
     !***
     !-----------------------------------------------------------------------
     use def_master,      only : kfl_rstar, therm, conce
     use def_chemic,      only : nZ_CMC_chm, rscal_inert_CMC_chm, &
                                 enthalp_CMC_chm, Yk_CMC_chm, temp_CMC_chm, &
                                 bvess_CMC_chm, kfl_solve_enth_CMC_chm, nvar_CMC_chm, &
                                 temp_inert_CMC_chm, temp_equil_CMC_chm, bvess_chm, &
                                 kfl_split_chm, kfl_post_gp_CMC_chm

     implicit none
     integer(ip)              :: ipoin,imixf,iclas
     real(rp)                 :: aux_dif

     if (kfl_rstar == 0) then
        ! Fill matrices Yk_CMC_chm and enthalp_CMC_chm
        do imixf = 1,nZ_CMC_chm
           do ipoin = 1,npoin
              do iclas = 1,nclas_chm
                 Yk_CMC_chm(imixf,ipoin,iclas) = bvess_CMC_chm(imixf,ipoin,iclas)
              end do
           end do
        end do

        if (kfl_solve_enth_CMC_chm /= 0_ip) then
           do imixf = 1,nZ_CMC_chm
              do ipoin = 1,npoin
                 enthalp_CMC_chm(imixf,ipoin) = bvess_CMC_chm(imixf,ipoin,nvar_CMC_chm)
              end do
           end do
        else
           do imixf = 1,nZ_CMC_chm
              enthalp_CMC_chm(imixf,1) = rscal_inert_CMC_chm(imixf,1)
           end do
        end if

        ! Compute conditional temperature at nodes
        do imixf = 1,nZ_CMC_chm

           aux_dif = temp_equil_CMC_chm(imixf) - temp_inert_CMC_chm(imixf)
           ! Assign to conce and therm
           do ipoin = 1,npoin
              if (kfl_solve_enth_CMC_chm /= 0_ip) then
                 therm(ipoin,1) = enthalp_CMC_chm(imixf,ipoin)
              else
                 therm(ipoin,1) = enthalp_CMC_chm(imixf,1)
              end if

              do iclas = 1,nclas_chm
                 conce(ipoin,iclas,1) = Yk_CMC_chm(imixf,ipoin,iclas)
              end do

              ! Just an initial temperature to start iterating from
              temp_CMC_chm(imixf,ipoin) = temp_inert_CMC_chm(imixf) + &
                 bvess_chm(nvar_CMC_chm+1_ip,ipoin) * aux_dif
           end do

           call chm_calc_temp_CMC(imixf)
        end do
     end if

     ! Compute conditional properties
     call chm_update_properties_CMC

     if( INOTMASTER ) then
        if (kfl_split_chm > 0_ip)   call chm_heatRelease_field_CMC   ! Compute heat release at nodes
     end if

     ! Integrations
     if (kfl_post_gp_CMC_chm == 0_ip)   call chm_integrate_flow_var_points_CMC
     call chm_integrate_flow_var_gauss_CMC
     call chm_nodal_projection_CMC

  end subroutine chm_initialization_domain_CMC


  !!!subroutine chm_initialize_alpha_CMC
  !!!   !-----------------------------------------------------------------------
  !!!   !****f* Chemic/chm_initialize_alpha_CMC
  !!!   ! NAME
  !!!   !    chm_initialize_alpha_CMC
  !!!   ! DESCRIPTION
  !!!   !    Initialize alpha values when using homogeneous reactor solutions for
  !!!   !    initial and boundary conditions.
  !!!   ! USED BY
  !!!   !
  !!!   !***
  !!!   !-----------------------------------------------------------------------

  !!!   use def_chemic,             only :  alpha_val_spec_CMC_chm, bvess_chm, &
  !!!                                       kfl_fixno_chm, nvar_CMC_chm

  !!!   implicit none
  !!!   integer(ip)                      :: ipoin, iclas

  !!!   alpha_val_spec_CMC_chm = 0.0_rp

  !!!   do ipoin = 1, npoin
  !!!      do iclas = 1, nclas_chm
  !!!         if( kfl_fixno_chm(1,ipoin) > 0_ip ) then
  !!!            alpha_val_spec_CMC_chm(1_ip,ipoin,iclas) = bvess_chm(nvar_CMC_chm+1,ipoin)
  !!!         end if
  !!!      end do
  !!!   end do

  !!!end subroutine chm_initialize_alpha_CMC


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!! M E M O R Y   A L L O C A T I O N !!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine chm_allocate_memory_CMC
     !-----------------------------------------------------------------------
     !****f* Chemic/chm_allocate_memory_CMC
     ! NAME
     !    chm_allocate_memory_CMC
     ! DESCRIPTION
     !    Allocate memory for most of the variables used by CMC.
     ! USED BY
     !    chm_solmem
     !***
     !-----------------------------------------------------------------------

!     use def_parame
!     use def_inpout
     use def_master,             only : condu_gp, sphec_gp, visco_gp, densi_gp, condk, densi, lescl, sphec_gp, sphek, tempe, therm,&
                                        visco
     use def_domain,             only : nelem, ltype, ngaus, gp_total_cloud_gp
!     use def_solver
     use def_chemic,             only : condu_gp_CMC_chm, sphec_gp_CMC_chm, visco_gp_CMC_chm, spvol_gp_CMC_chm,&
                                        react_scalars_gp_CMC_chm, mass_gp, aux_cond_fields_CMC_chm, aux_interp_fields_CMC_chm,&
                                        aux_val_CMC_chm, av_enthalp_CMC_chm, av_enthalp_int_CMC_chm,&
                                        av_temp_CMC_chm, av_temp_int_CMC_chm, av_visco_lam_int_CMC_chm, av_Yk_CMC_chm,&
                                        av_Yk_int_CMC_chm, coeff_cp_k, condu_gp_CMC_chm, condu_int_CMC_chm, densi_int_CMC_chm,&
                                        deriv2_enthalp_CMC_chm, deriv2_Yk_CMC_chm, enthalp_CMC_chm, enthalp_int_CMC_chm,&
                                        fact_BC_RK_CMC_chm, fact_diag_low_RK_CMC_chm, fact_diag_RK_CMC_chm,&
                                        fact_diag_up_RK_CMC_chm, freq_CFD_coup_CMC_chm, grad_Zavg_CMC_chm, homog_reactors_CMC_chm,&
                                        hrr_avg_chm, hrr_chm, hrr_mass_cond_CMC_chm, kfl_avg_cond_CMC_chm,&
                                        kfl_bc_type_spec_CMC_chm, kfl_mesh_interp_CMC_chm, kfl_post_gp_CMC_chm,&
                                        kfl_solve_cond_CMC_chm, kfl_solve_enth_CMC_chm, kfl_split_CFD_CMC,&
                                        kfl_transfer_condField_CMC_chm, mass_gp, matr_scal_CMC_chm, normal_CMC_mesh_CMC_chm,&
                                        nspec_cond_write_CMC_chm, nspec_transf_CMC_chm, nspec_uncond_write_CMC_chm, nZ_CMC_chm,&
                                        nZ_write_CMC_chm, order_RK_CMC_chm, react_scalars_gp_CMC_chm, rec_CFD_from_CMC_chm,&
                                        rec_CFD_from_CMC_post_chm, rec_CFD_old_from_CMC_chm, rec_CMC_from_CFD_chm,&
                                        send_CFD_to_CMC_chm, send_CMC_to_CFD_chm, send_CMC_to_CFD_post_chm, sphea_int_CMC_chm,&
                                        sphec_gp_CMC_chm, spvol_gp_CMC_chm, src_Yk_CMC_chm, src_Yk_int_CMC_chm,&
                                        sum_Yk_ext_cond_CMC_chm, t_chem_integ_CMC_chm, temp_int_CMC_chm, Text_cond_CMC_chm,&
                                        transf_entha_CMC_chm, turb_kin_visc_CMC_chm, veloc_CMC_chm, visco_gp_CMC_chm,&
                                        visco_lam_int_CMC_chm, w_k, weights_ext_RK_CMC_chm, weights_int_RK_CMC_chm, Xtot_CMC_chm,&
                                        Xtot_const_CMC_chm, Xtot_whole_mesh_CMC_chm, xZr_chm, xZs_chm, Yk_int_CMC_chm,&
                                        Zavg_CMC_chm, Zavg_const_CMC_chm, Zvar_CMC_chm, Zvar_const_CMC_chm, avZ_chm, avZV_chm
     use mod_memory,             only : memory_alloca
     use mod_output_postprocess, only : output_postprocess_check_variable_postprocess

     implicit none
     integer(ip) :: ielem, pelty, pgaus
     integer(ip) :: nvar_int

     if (kfl_solve_cond_CMC_chm == 1) then
        nullify(src_Yk_CMC_chm)
        nullify(Yk_int_CMC_chm)
        nullify(densi_int_CMC_chm)
        nullify(visco_lam_int_CMC_chm)
        nullify(condu_int_CMC_chm)
        nullify(sphea_int_CMC_chm)
        nullify(enthalp_int_CMC_chm)
        nullify(temp_int_CMC_chm)
        nullify(src_Yk_int_CMC_chm)
        nullify(hrr_mass_cond_CMC_chm)
        nullify(deriv2_Yk_CMC_chm)
        nullify(veloc_CMC_chm)
        nullify(Zavg_CMC_chm)
        nullify(Zvar_CMC_chm)
        nullify(Xtot_CMC_chm)
        nullify(grad_Zavg_CMC_chm)
        nullify(turb_kin_visc_CMC_chm)
        nullify(kfl_bc_type_spec_CMC_chm)
        nullify(W_k)
        nullify(coeff_cp_k)
        nullify(hrr_chm)
        nullify(hrr_avg_chm)
        nullify(send_CMC_to_CFD_chm)
        nullify(send_CMC_to_CFD_post_chm)
        nullify(rec_CMC_from_CFD_chm)
        nullify(weights_int_RK_CMC_chm)
        nullify(weights_ext_RK_CMC_chm)
        nullify(fact_diag_low_RK_CMC_chm)
        nullify(fact_diag_RK_CMC_chm)
        nullify(fact_diag_up_RK_CMC_chm)
        nullify(homog_reactors_CMC_chm)
        nullify(aux_val_CMC_chm)
        nullify(matr_scal_CMC_chm)
        nullify(av_Yk_CMC_chm)
        nullify(av_enthalp_CMC_chm)
        nullify(av_temp_CMC_chm)
        nullify(av_Yk_int_CMC_chm)
        nullify(av_enthalp_int_CMC_chm)
        nullify(av_temp_int_CMC_chm)
        nullify(av_visco_lam_int_CMC_chm)
        nullify(t_chem_integ_CMC_chm)

        if (kfl_mesh_interp_CMC_chm == 1_ip) then
           nullify(normal_CMC_mesh_CMC_chm)
           nullify(aux_cond_fields_CMC_chm)
           nullify(aux_interp_fields_CMC_chm)
           if (kfl_avg_cond_CMC_chm == 1_ip) then
              nullify(Xtot_whole_mesh_CMC_chm)
           end if
        end if

        if (kfl_solve_enth_CMC_chm == 0) then
           nullify(enthalp_CMC_chm)
           call memory_alloca(mem_modul(1:2,modul),'CONDITIONAL_ENTHALPY_CMC_CHM','chm_solmem',enthalp_CMC_chm,nZ_CMC_chm,1_ip)
           enthalp_CMC_chm  = 0.0_rp
        end if
        call memory_alloca(mem_modul(1:2,modul),'CONDITIONAL_OMEGA_MASS_FRACT_CMC_CHM','chm_solmem',src_Yk_CMC_chm,nZ_CMC_chm,&
            npoin,nclas_chm)
        call memory_alloca(mem_modul(1:2,modul),'UNCONDITIONAL_MASS_FRACTIONS_CMC_CHM','chm_solmem',Yk_int_CMC_chm,npoin,nclas_chm)
        call memory_alloca(mem_modul(1:2,modul),'UNCONDITIONAL_DENSITY_CMC_CHM','chm_solmem',densi_int_CMC_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'UNCONDITIONAL_LAMINAR_VISC_CMC_CHM','chm_solmem',visco_lam_int_CMC_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'UNCONDITIONAL_THERMAL_CONDUCTIVITY_CMC_CHM','chm_solmem',condu_int_CMC_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'UNCONDITIONAL_SPECIFIC_HEAT_CMC_CHM','chm_solmem',sphea_int_CMC_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'UNCONDITIONAL_ENTHALPY_CMC_CHM','chm_solmem',enthalp_int_CMC_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'UNCONDITIONAL_TEMPERATURE_CMC_CHM','chm_solmem',temp_int_CMC_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'UNCONDITIONAL_OMEGA_MASS_FRACT_CMC_CHM','chm_solmem',src_Yk_int_CMC_chm,npoin,&
            nclas_chm)
        call memory_alloca(mem_modul(1:2,modul),'HEAT_RELEASE_CMC_CHM','chm_solmem',hrr_mass_cond_CMC_chm,nZ_CMC_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'DERIV2_Yk_CMC_CHM','chm_solmem',deriv2_Yk_CMC_chm,nZ_CMC_chm,npoin,nclas_chm)
        call memory_alloca(mem_modul(1:2,modul),'VELOCITY_CFD_CHM','chm_solmem',veloc_CMC_chm,ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'MIXT_FRAC_CFD_CHM','chm_solmem',Zavg_CMC_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'MIXT_FRAC_VAR_CFD_CHM','chm_solmem',Zvar_CMC_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'TOTAL_SCAL_DISSIP_RATE_CFD_CHM','chm_solmem',Xtot_CMC_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'MIXTURE_FRACTION_GRADIENT_CFD_CHM','chm_solmem',grad_Zavg_CMC_chm,ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'MASS_TURBULENT_DIFFUSION_COEFF_CFD_CHM','chm_solmem',turb_kin_visc_CMC_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'TYPE_BC_SPECIES_CMC_CHM','chm_solmem',kfl_bc_type_spec_CMC_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'W_K','chm_solmem',W_k,nspec_chm)
        call memory_alloca(mem_modul(1:2,modul),'COEFF_CP_K','chm_solmem',coeff_cp_k,nspec_chm,15_ip)
        call memory_alloca(mem_modul(1:2,modul),'HRR_CHM','chm_solmem',hrr_chm,npoin)
        if (order_RK_CMC_chm >= 2_ip) then
           call memory_alloca(mem_modul(1:2,modul),'INTERNAL_WEIGHTS_RK_CHM','chm_solmem',weights_int_RK_CMC_chm,&
               order_RK_CMC_chm-1_ip, order_RK_CMC_chm-1_ip)
        end if
        call memory_alloca(mem_modul(1:2,modul),'EXTERNAL_WEIGHTS_RK_CHM','chm_solmem',weights_ext_RK_CMC_chm,order_RK_CMC_chm)
        call memory_alloca(mem_modul(1:2,modul),'FACTOR_LOWER_ELEM_DIAG_RK_CHM','chm_solmem',fact_diag_low_RK_CMC_chm,&
            nZ_CMC_chm-3_ip)
        call memory_alloca(mem_modul(1:2,modul),'FACTOR_ELEM_DIAG_RK_CHM','chm_solmem',fact_diag_RK_CMC_chm,nZ_CMC_chm-2_ip)
        call memory_alloca(mem_modul(1:2,modul),'FACTOR_UPPER_ELEM_DIAG_RK_CHM','chm_solmem',fact_diag_up_RK_CMC_chm,&
            nZ_CMC_chm-3_ip)
        call memory_alloca(mem_modul(1:2,modul),'MIN_MAX_COND_TEMP_CHM','chm_solmem',Text_cond_CMC_chm,2_ip*nZ_CMC_chm)
        call memory_alloca(mem_modul(1:2,modul),'MIN_MAX_COND_SUM_YK_CHM','chm_solmem',sum_Yk_ext_cond_CMC_chm,2_ip*nZ_CMC_chm)
        call memory_alloca(mem_modul(1:2,modul),'AUXILIAR_MATRIX_CHM','chm_solmem',aux_val_CMC_chm,npoin)

        !!!if (kfl_bc_init_method_CMC_chm == 2_ip) then
        !!!   call memory_alloca(mem_modul(1:2,modul),'HOMOGENEOUS_REACTORS_CHM','chm_solmem',homog_reactors_CMC_chm,nZ_CMC_chm,&
        !!!       n_alpha_grid_CMC_chm, nclas_chm)
        !!!   homog_reactors_CMC_chm = 0.0_rp
        !!!end if

        if (kfl_mesh_interp_CMC_chm == 1_ip) then
            call memory_alloca(mem_modul(1:2,modul),'NORMALIZATION_INTEG_CMC_MESH','chm_plugin',normal_CMC_mesh_CMC_chm,npoin)
        end if

        if (kfl_split_CFD_CMC == 1_ip) then

           if (kfl_transfer_condField_CMC_chm == 1_ip) then
              call memory_alloca(mem_modul(1:2,modul),'SEND_CMC_TO_CFD','chm_plugin',send_CMC_to_CFD_chm,6_ip,nZ_CMC_chm,npoin)
              call memory_alloca(mem_modul(1:2,modul),'SEND_CMC_TO_CFD_POST','chm_plugin',send_CMC_to_CFD_post_chm,nZ_CMC_chm,npoin)
           else
              call memory_alloca(mem_modul(1:2,modul),'SEND_CMC_TO_CFD','chm_plugin',send_CMC_to_CFD_chm,5_ip,1_ip,npoin)
           end if

           if (kfl_mesh_interp_CMC_chm == 1_ip) then
              ! Information received at Gauss points of CMC mesh for interpolation (CMC mesh /= CFD mesh)
              call memory_alloca(mem_modul(1:2,modul),'RECEIVE_CMC_FROM_CFD','chm_plugin',rec_CMC_from_CFD_chm,4_ip+2_ip*ndime,&
                  gp_total_cloud_gp)

              if (kfl_avg_cond_CMC_chm == 1_ip) then
                 call memory_alloca(mem_modul(1:2,modul),'AUX_COND_FIELDS_CMC_MESH','chm_solmem',aux_cond_fields_CMC_chm,2_ip+2_ip&
                     * ndime,gp_total_cloud_gp)
                 call memory_alloca(mem_modul(1:2,modul),'AUX_INTERP_FIELDS_CMC_MESH','chm_solmem',aux_interp_fields_CMC_chm,npoin,&
                     2_ip+2_ip*ndime)
                 call memory_alloca(mem_modul(1:2,modul),'SCALAR_DISSIP_RATE_ON_WHOLE_MESH','chm_solmem',Xtot_whole_mesh_CMC_chm,&
                     nZ_CMC_chm,npoin)

              else
                 call memory_alloca(mem_modul(1:2,modul),'AUX_INTERP_FIELDS_CMC_MESH','chm_solmem',aux_interp_fields_CMC_chm,npoin,&
                     4_ip+2_ip*ndime)
              end if

           else
              ! Information received at nodes of CMC mesh (CMC mesh = CFD mesh)
              call memory_alloca(mem_modul(1:2,modul),'RECEIVE_CMC_FROM_CFD','chm_plugin',rec_CMC_from_CFD_chm,4_ip+2_ip*ndime,&
                  npoin)
           end if

        end if

        if ( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVYCN') ) then
           call memory_alloca(mem_modul(1:2,modul),'AVG_CONDITIONAL_MASS_FRAC_CMC','chm_solmem',av_Yk_CMC_chm,nZ_write_CMC_chm,&
               npoin,nspec_cond_write_CMC_chm)
           av_Yk_CMC_chm = 0.0_rp
        end if

        if ( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVYUN') ) then
           call memory_alloca(mem_modul(1:2,modul),'AVG_UNCONDITIONAL_MASS_FRAC_CMC','chm_solmem',av_Yk_int_CMC_chm,npoin,&
               nspec_uncond_write_CMC_chm)
           av_Yk_int_CMC_chm = 0.0_rp
        end if

        if ( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVHCN') .and. kfl_solve_enth_CMC_chm == 1 ) then
           call memory_alloca(mem_modul(1:2,modul),'AVG_CONDITIONAL_ENTHALPY_CMC','chm_solmem',av_enthalp_CMC_chm,nZ_write_CMC_chm,&
               npoin)
           av_enthalp_CMC_chm = 0.0_rp
        end if

        if ( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVHUN') ) then
           call memory_alloca(mem_modul(1:2,modul),'AVG_UNCONDITIONAL_ENTHALPY_CMC','chm_solmem',av_enthalp_int_CMC_chm,npoin)
           av_enthalp_int_CMC_chm = 0.0_rp
        end if

        if ( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVTCN') ) then
           call memory_alloca(mem_modul(1:2,modul),'AVG_CONDITIONAL_TEMPERATURE_CMC','chm_solmem',av_temp_CMC_chm,nZ_write_CMC_chm,&
               npoin)
           av_temp_CMC_chm = 0.0_rp
        end if

        if ( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVTUN') ) then
           call memory_alloca(mem_modul(1:2,modul),'AVG_UNCONDITIONAL_TEMPERATURE_CMC','chm_solmem',av_temp_int_CMC_chm,npoin)
           av_temp_int_CMC_chm = 0.0_rp
        end if

        if ( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVVIS') ) then
           call memory_alloca(mem_modul(1:2,modul),'AVG_UNCONDITIONAL_VISCOSITY_CMC','chm_solmem',av_visco_lam_int_CMC_chm,npoin)
           av_visco_lam_int_CMC_chm = 0.0_rp
        end if

        if ( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVHRR') ) then
           call memory_alloca(mem_modul(1:2,modul),'HRR_AVG_CHM','chm_solmem',hrr_avg_chm,npoin)
           hrr_avg_chm = 0.0_rp
        end if

        call memory_alloca(mem_modul(1:2,modul),'TIME_CHEMICAL_INTEGTATION_CMC','chm_solmem',t_chem_integ_CMC_chm,nZ_CMC_chm,npoin)

        ! Initialize previous matrices
        src_Yk_CMC_chm            = 0.0_rp
        Yk_int_CMC_chm            = 0.0_rp
        densi_int_CMC_chm         = 1.0_rp
        visco_lam_int_CMC_chm     = 0.0_rp
        condu_int_CMC_chm         = 0.0_rp
        sphea_int_CMC_chm         = 0.0_rp
        enthalp_int_CMC_chm       = 0.0_rp
        temp_int_CMC_chm          = 300.0_rp
        src_Yk_int_CMC_chm        = 0.0_rp
        hrr_mass_cond_CMC_chm     = 0.0_rp
        deriv2_Yk_CMC_chm         = 0.0_rp
        veloc_CMC_chm             = 0.0_rp
        Zavg_CMC_chm              = Zavg_const_CMC_chm
        Zvar_CMC_chm              = Zvar_const_CMC_chm
        Xtot_CMC_chm              = Xtot_const_CMC_chm
        grad_Zavg_CMC_chm         = 0.0_rp
        turb_kin_visc_CMC_chm     = 0.0_rp
        kfl_bc_type_spec_CMC_chm  = 1_ip
        W_k                       = 0.0_rp
        coeff_cp_k                = 0.0_rp
        hrr_chm                   = 0.0_rp
        if (order_RK_CMC_chm >= 2_ip)  weights_int_RK_CMC_chm = 0.0_rp
        weights_ext_RK_CMC_chm    = 0.0_rp
        fact_diag_low_RK_CMC_chm  = 0.0_rp
        fact_diag_RK_CMC_chm      = 0.0_rp
        fact_diag_up_RK_CMC_chm   = 0.0_rp
        fact_BC_RK_CMC_chm        = 0.0_rp
        Text_cond_CMC_chm         = 0.0_rp
        sum_Yk_ext_cond_CMC_chm   = 1.0_rp
        aux_val_CMC_chm           = 0.0_rp
        t_chem_integ_CMC_chm      = 0.0_rp


        if (.not. associated(therm)) then
           nullify(therm)
           call memory_alloca(mem_modul(1:2,modul),'THERM_CMC_CHM','chm_solmem',therm,npoin,1_ip)
           therm = 0.0_rp
        end if

        if (.not. associated(tempe)) then
           nullify(tempe)
           call memory_alloca(mem_modul(1:2,modul),'THERM_CMC_CHM','chm_solmem',tempe,npoin,1_ip)
           tempe = 300.0_rp
        end if

        if (kfl_solve_enth_CMC_chm /= 0) then
           nullify(deriv2_enthalp_CMC_chm)
           call memory_alloca(mem_modul(1:2,modul),'DERIV2_ENTHALP_CMC_CHM','chm_solmem',deriv2_enthalp_CMC_chm,nZ_CMC_chm,npoin)
           deriv2_enthalp_CMC_chm = 0.0_rp
        end if

        call memory_alloca(mem_modul(1:2,modul),'CONDU_GP_CMC_CHM','chm_solmem',condu_gp_CMC_chm,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'CONDU_GP_CMC_CHM % A','chm_solmem',condu_gp_CMC_chm(ielem)%a,pgaus,&
               nZ_CMC_chm,1_ip)
           condu_gp_CMC_chm(ielem)%a=0.0_rp
        end do

        call memory_alloca(mem_modul(1:2,modul),'CONDU_GP','chm_solmem',condu_gp,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'CONDU_GP % A','chm_solmem',condu_gp(ielem)%a,pgaus,1_ip,1_ip)
           condu_gp(ielem)%a=1.0_rp
        end do

        call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP_CMC_CHM','chm_solmem',sphec_gp_CMC_chm,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP_CMC_CHM % A','chm_solmem',sphec_gp_CMC_chm(ielem)%a,pgaus,&
               nZ_CMC_chm,1_ip)
           sphec_gp_CMC_chm(ielem)%a=1.0_rp
        end do

        call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP','chm_solmem',sphec_gp,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP % A','chm_solmem',sphec_gp(ielem)%a,pgaus,1_ip,1_ip)
           sphec_gp(ielem)%a=1.0_rp
        end do

        call memory_alloca(mem_modul(1:2,modul),'VISCO_GP_CMC_CHM','chm_solmem',visco_gp_CMC_chm,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'VISCO_GP_CMC_CHM % A','chm_solmem',visco_gp_CMC_chm(ielem)%a,pgaus,&
               nZ_CMC_chm,1_ip)
           visco_gp_CMC_chm(ielem)%a=0.0_rp
        end do

        call memory_alloca(mem_modul(1:2,modul),'VISCO_GP','chm_solmem',visco_gp,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'VISCO_GP % A','chm_solmem',visco_gp(ielem)%a,pgaus,1_ip,1_ip)
           visco_gp(ielem)%a=0.0_rp
        end do

        call memory_alloca(mem_modul(1:2,modul),'SPEC_VOL_GP_CMC_CHM','chm_solmem',spvol_gp_CMC_chm,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'SPEC_VOL_GP_CMC_CHM % A','chm_solmem',spvol_gp_CMC_chm(ielem)%a,pgaus,&
               nZ_CMC_chm,1_ip)
           spvol_gp_CMC_chm(ielem)%a=1.0_rp
        end do

        call memory_alloca(mem_modul(1:2,modul),'DENSI_GP','chm_solmem',densi_gp,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'DENSI_GP % A','chm_solmem',densi_gp(ielem)%a,pgaus,1_ip,1_ip)
           densi_gp(ielem)%a=1.0_rp
        end do

        if (kfl_post_gp_CMC_chm == 1_ip) then
           nvar_int  = 2_ip * nclas_chm + 4_ip
           ! It contains: <T>, <MW>, <hrr>, <Yk>, <wYk>, <h>
           call memory_alloca(mem_modul(1:2,modul),'REACTING_SCALARS_GP','chm_solmem',react_scalars_gp_CMC_chm,nelem)
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              call memory_alloca(mem_modul(1:2,modul),'REACTING_SCALARS_GP % A','chm_solmem',react_scalars_gp_CMC_chm(ielem)%a,&
                  nvar_int,pgaus,1_ip)
              react_scalars_gp_CMC_chm(ielem)%a=0.0_rp
           end do
        end if

        if (kfl_transfer_condField_CMC_chm == 1_ip) then
           nullify(densi)
           nullify(visco)
           nullify(condk)
           nullify(sphek)
           call memory_alloca(mem_modul(1:2,modul),'DENSITY_CHM','chm_solmem',densi,npoin,1_ip)
           call memory_alloca(mem_modul(1:2,modul),'VISCOSITY_CHM','chm_solmem',visco,npoin,1_ip)
           call memory_alloca(mem_modul(1:2,modul),'CONDUCTIVITY_CHM','chm_solmem',condk,npoin,1_ip)
           call memory_alloca(mem_modul(1:2,modul),'SPECIFIC_HEAT_CHM','chm_solmem',sphek,npoin,1_ip)
           densi   = 1.0_rp
           visco   = 0.0_rp
           condk   = 0.0_rp
           sphek   = 0.0_rp
        end if

     else  ! Solve mixture fraction, its variance, etc. for its coupling with CMC
        nullify(xZr_chm)
        nullify(xZs_chm)
        nullify(lescl)
        nullify(densi)
        nullify(visco)
        nullify(condk)
        nullify(sphek)
        nullify(tempe)
        nullify(send_CFD_to_CMC_chm)
        nullify(rec_CFD_from_CMC_chm)
        nullify(rec_CFD_old_from_CMC_chm)
        nullify(rec_CFD_from_CMC_post_chm)
        nullify(enthalp_int_CMC_chm)
        nullify(avZ_chm)
        nullify(avZV_chm)
        nullify(hrr_avg_chm)
        nullify(av_Yk_int_CMC_chm)
        nullify(av_enthalp_int_CMC_chm)
        nullify(av_temp_int_CMC_chm)
        nullify(av_visco_lam_int_CMC_chm)
        !!!! ¿DESCOMENTAR?
        !!!!nullify(turb_kin_visco_CFD_chm)
        !!!! nullify(grad_Zavg_CFD_chm)

        call memory_alloca(mem_modul(1:2,modul),'XZR_CHM','chm_solmem',xZr_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'XZS_CHM','chm_solmem',xZs_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'LESCL_CHM','chm_solmem',lescl,npoin)
        call memory_alloca(mem_modul(1:2,modul),'DENSITY_CHM','chm_solmem',densi,npoin,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'VISCOSITY_CHM','chm_solmem',visco,npoin,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'CONDUCTIVITY_CHM','chm_solmem',condk,npoin,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'SPECIFIC_HEAT_CHM','chm_solmem',sphek,npoin,1_ip)

        if (kfl_split_CFD_CMC == 1) then
           call memory_alloca(mem_modul(1:2,modul),'SEND_CFD_TO_CMC','chm_solmem',send_CFD_to_CMC_chm,4_ip+2_ip*ndime,npoin)
           if (kfl_transfer_condField_CMC_chm == 1_ip) then
              if (freq_CFD_coup_CMC_chm > 1_ip) then
                 call memory_alloca(mem_modul(1:2,modul),'RECEIVE_CFD_OLD_FROM_CMC','chm_solmem',rec_CFD_old_from_CMC_chm,6_ip,&
                     nZ_CMC_chm,npoin)
              end if
              call memory_alloca(mem_modul(1:2,modul),'RECEIVE_CFD_FROM_CMC','chm_solmem',rec_CFD_from_CMC_chm,6_ip,nZ_CMC_chm,&
                  npoin)
              call memory_alloca(mem_modul(1:2,modul),'RECEIVE_CFD_FROM_CMC_POST','chm_solmem',rec_CFD_from_CMC_post_chm,&
                  nZ_CMC_chm,npoin)
              call memory_alloca(mem_modul(1:2,modul),'UNCONDITIONAL_HEAT_RELEASE_CMC_CHM','chm_solmem',hrr_chm,npoin)
              if (nspec_transf_CMC_chm > 0_ip) &
                 call memory_alloca(mem_modul(1:2,modul),'UNCONDITIONAL_MASS_FRACTIONS_CMC_CHM','chm_solmem',Yk_int_CMC_chm,npoin,&
                     nspec_transf_CMC_chm)
              if (transf_entha_CMC_chm > 0_ip) &
                 call memory_alloca(mem_modul(1:2,modul),'UNCONDITIONAL_ENTHALPY_CMC_CHM','chm_solmem',enthalp_int_CMC_chm,npoin)
           else
              call memory_alloca(mem_modul(1:2,modul),'RECEIVE_CFD_FROM_CMC','chm_solmem',rec_CFD_from_CMC_chm,5_ip,1_ip,npoin)
           end if

           if ( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVZ  ') ) then
              call memory_alloca(mem_modul(1:2,modul),'AVZ_CHM','chm_solmem',avZ_chm, npoin)
              avZ_chm = 0.0_rp
           end if

           if ( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVZV ') ) then
              call memory_alloca(mem_modul(1:2,modul),'AVZV_CHM','chm_solmem',avZV_chm, npoin)
              avZV_chm = 0.0_rp
           end if

           if ( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVHRR') ) then
              call memory_alloca(mem_modul(1:2,modul),'AVHRR_CHM','chm_solmem',hrr_avg_chm, npoin)
              hrr_avg_chm = 0.0_rp
           end if

           if ( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVYUN') ) then
              call memory_alloca(mem_modul(1:2,modul),'AVYUN_CHM','chm_solmem',av_Yk_int_CMC_chm, npoin, nspec_transf_CMC_chm)
              av_Yk_int_CMC_chm = 0.0_rp
           end if

           if ( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVHUN') ) then
              call memory_alloca(mem_modul(1:2,modul),'AVHUN_CHM','chm_solmem',av_enthalp_int_CMC_chm, npoin)
              av_enthalp_int_CMC_chm = 0.0_rp
           end if

           if ( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVTUN') ) then
              call memory_alloca(mem_modul(1:2,modul),'AVTUN_CHM','chm_solmem',av_temp_int_CMC_chm, npoin)
              av_temp_int_CMC_chm = 0.0_rp
           end if

           if ( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVVIS') ) then
              call memory_alloca(mem_modul(1:2,modul),'AVVIS_CHM','chm_solmem',av_visco_lam_int_CMC_chm, npoin)
              av_visco_lam_int_CMC_chm = 0.0_rp
           end if

           !!!!!!!! ¿DESCOMENTAR?
           !!!!call memory_alloca(mem_modul(1:2,modul),'GRADIENT_MIXT_FRACT_CFD','chm_solmem',grad_Zavg_CFD_chm,ndime,npoin)
           !!!!allocate(grad_Zavg_CFD_chm(ndime,npoin))
           !!!!grad_Zavg_CFD_chm(1:ndime,1:npoin) = 0.0_rp
           !!!!if (turmu_ker % kfl_exist /= 0_ip) then
           !!!!   allocate(turb_kin_visco_CFD_chm(npoin))
           !!!!end if

        end if


        xZr_chm = 0.0_rp
        xZs_chm = 0.0_rp
        lescl   = 0.0_rp
        densi   = 1.0_rp
        visco   = 0.0_rp
        condk   = 0.0_rp
        sphek   = 0.0_rp


        if (.not. associated(tempe)) then  ! Used only by some subroutines from Nastin
           nullify(tempe)
           call memory_alloca(mem_modul(1:2,modul),'TEMPERATURE_CHM','chm_solmem',tempe,npoin,1_ip)
           tempe = 300.0_rp
        end if


        ! To use flamelet routines we need to allocate mass_gp
        call memory_alloca(mem_modul(1:2,modul),'MASS_GP','chm_solmem',mass_gp,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'MASS_GP % A','chm_solmem',mass_gp(ielem)%a,pgaus,nclas_chm,1_ip)
           mass_gp(ielem)%a=0.0_rp
        end do
     end if
  end subroutine chm_allocate_memory_CMC



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! B O U N D A R Y   C O N D I T I O N S !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine chm_bc_type_species_CMC
     !-----------------------------------------------------------------------
     !****f* Chemic/chm_bc_type_species_CMC
     ! NAME
     !    chm_bc_type_species_CMC
     ! DESCRIPTION
     !    Idenitify if all the species have the same type of boundary conditions
     !    for the physical points of the mesh
     ! USED BY
     !    chm_iniunk
     !***
     !-----------------------------------------------------------------------

     use def_chemic,             only :  kfl_bc_type_spec_CMC_chm, kfl_fixno_chm

     implicit none
     integer(ip)                      :: ipoin, iclas, type_clas1

     do ipoin = 1, npoin
        type_clas1 = kfl_fixno_chm(1,ipoin)
        do iclas = 2, nclas_chm  ! Since we do not explicitly give Dirichlet BCs for soot
           if (kfl_fixno_chm(iclas,ipoin) /= type_clas1) then
              kfl_bc_type_spec_CMC_chm(ipoin) = 0
           end if
        end do
     end do

  end subroutine chm_bc_type_species_CMC


  subroutine chm_save_unconditional_fields_bc_CMC
     !-----------------------------------------------------------------------
     !****f* Chemic/chm_save_unconditional_fields_bc_CMC
     ! NAME
     !    chm_save_unconditional_fields_bc_CMC
     ! DESCRIPTION
     !    Assign bvess_chm values to bvess_ufield_CMC_chm when starting
     !    the simulation.
     ! USED BY
     !    chm_iniunk
     !***
     !-----------------------------------------------------------------------

     use def_chemic,             only :  bvess_chm, bvess_ufield_CMC_chm, &
                                         nvar_CMC_chm, kfl_bc_alpha_CMC_chm, &
                                         kfl_solve_enth_CMC_chm

     implicit none
     integer(ip)                      :: ipoin, iclas

     if (kfl_bc_alpha_CMC_chm == 1_ip) then
        bvess_ufield_CMC_chm(:,1) = bvess_chm(1,:)
        if (kfl_solve_enth_CMC_chm /= 0_ip) then
           bvess_ufield_CMC_chm(:,2) = bvess_chm(nvar_CMC_chm,:)
        end if
     else
        do ipoin = 1, npoin
           do iclas = 1, nvar_CMC_chm
              bvess_ufield_CMC_chm(ipoin,iclas) = bvess_chm(iclas,ipoin)
           end do
        end do
     end if

  end subroutine chm_save_unconditional_fields_bc_CMC


  subroutine chm_get_cond_species_all_average_linear_CMC(ipoin,Yk_avg,Zavg,Zvar,Yk_prof)
     !-----------------------------------------------------------------------
     !****f* Chemic/chm_get_cond_species_all_average_linear_CMC
     ! NAME
     !    chm_get_cond_species_all_average_linear_CMC
     ! DESCRIPTION
     !    Find conditional species profiles along mixture fraction from average
     !    values with linear profile in alfa. Do this operation for all the
     !    species at the same time.
     ! USED BY
     !    chm_get_cond_fields_bc_Dir_linear_CMC
     !***
     !-----------------------------------------------------------------------

     use def_chemic,             only :  nZ_CMC_chm, rscal_inert_CMC_chm, &
                                         rscal_equil_CMC_chm, Yk_CMC_chm, &
                                         max_diff_Yk_inert_eq_CMC_chm

     implicit none
     real(rp), parameter              :: small = 5.0e-4_rp
     integer(ip), intent(in)          :: ipoin
     real(rp), intent(in)             :: Yk_avg(nspec_chm)
     real(rp), intent(in)             :: Zavg
     real(rp), intent(in)             :: Zvar
     real(rp), intent(out)            :: Yk_prof(nZ_CMC_chm,nspec_chm)
     integer(ip)                      :: iclas, imixf
     real(rp)                         :: inert_eq_int(2_ip*nspec_chm)
     real(rp)                         :: alpha
     real(rp)                         :: Yk_sum
     real(rp)                         :: diff_aux

     Yk_prof = 0.0_rp

     ! Integrate inert and equilibrium profiles
     call chm_mxt_fr_integr_previous_steps_CMC(2_ip*nspec_chm, &
            (/rscal_inert_CMC_chm(1:nZ_CMC_chm,2:nspec_chm+1_ip), &
            rscal_equil_CMC_chm(1:nZ_CMC_chm,2:nspec_chm+1_ip)/), &
            Zavg, Zvar, inert_eq_int(1:2_ip*nspec_chm))

     do iclas = 1, nspec_chm
        diff_aux = inert_eq_int(iclas+nspec_chm) - inert_eq_int(iclas)
        if (abs(diff_aux) < small * max_diff_Yk_inert_eq_CMC_chm(iclas)) then
           ! Take the profile from previous time step
           Yk_prof(1:nZ_CMC_chm,iclas) = Yk_CMC_chm(1:nZ_CMC_chm,ipoin,iclas)
        else
           alpha = (Yk_avg(iclas) - inert_eq_int(iclas)) / diff_aux
           alpha = min(1.0_rp, max(0.0_rp, alpha))
           Yk_prof(1:nZ_CMC_chm,iclas) = rscal_inert_CMC_chm(1:nZ_CMC_chm,iclas+1_ip) + &
              alpha * (rscal_equil_CMC_chm(1:nZ_CMC_chm,iclas+1_ip) - &
              rscal_inert_CMC_chm(1:nZ_CMC_chm,iclas+1_ip))
        end if
     end do

     ! Normalize in case mass fractions do not sum 1
     do imixf = 2, nZ_CMC_chm-1_ip
        Yk_sum = 0.0_rp
        ! Notice that soot mass fractions are set to 0 at Dirichlet
        do iclas = 1, nspec_chm
           Yk_sum = Yk_sum + Yk_prof(imixf,iclas)
        end do
        if (Yk_sum /= 0.0_rp) then
           Yk_prof(imixf,1:nspec_chm) = Yk_prof(imixf,1:nspec_chm) / Yk_sum
        end if
     end do

  end subroutine chm_get_cond_species_all_average_linear_CMC


  subroutine chm_get_cond_species_average_linear_CMC(iclas,ipoin,Yk_avg,Zavg,Zvar,Yk_prof)
     !-----------------------------------------------------------------------
     !****f* Chemic/chm_get_cond_field_average_CMC
     ! NAME
     !    chm_get_cond_field_average_CMC
     ! DESCRIPTION
     !    Find conditional species profiles along mixture fraction from average
     !    values with linear profile in alfa. Do this operation for just one
     !    single species.
     ! USED BY
     !    chm_get_cond_fields_bc_Dir_linear_CMC
     !***
     !-----------------------------------------------------------------------

     use def_chemic,             only :  nZ_CMC_chm, rscal_inert_CMC_chm, &
                                         rscal_equil_CMC_chm, Yk_CMC_chm, &
                                         max_diff_Yk_inert_eq_CMC_chm

     implicit none
     real(rp), parameter              :: small = 5.0e-4_rp
     integer(ip), intent(in)          :: iclas
     integer(ip), intent(in)          :: ipoin
     real(rp), intent(in)             :: Yk_avg
     real(rp), intent(in)             :: Zavg
     real(rp), intent(in)             :: Zvar
     real(rp), intent(out)            :: Yk_prof(nZ_CMC_chm)
     real(rp)                         :: inert_eq_int(2)
     real(rp)                         :: alpha
     real(rp)                         :: diff_aux

     Yk_prof = 0.0_rp

     ! Integrate inert and equilibrium profiles
     call chm_mxt_fr_integr_previous_steps_CMC(2_ip, &
            (/rscal_inert_CMC_chm(:,iclas+1_ip), rscal_equil_CMC_chm(:,iclas+1_ip)/), &  ! Enthalpy is the first entry
            Zavg, Zvar, inert_eq_int)

     diff_aux = inert_eq_int(2) - inert_eq_int(1)
     if (abs(diff_aux) < small * max_diff_Yk_inert_eq_CMC_chm(iclas)) then
        ! Take the profile from previous time step
        Yk_prof(1:nZ_CMC_chm) = Yk_CMC_chm(1:nZ_CMC_chm,ipoin,iclas)
     else
        alpha = (Yk_avg - inert_eq_int(1)) / diff_aux
        alpha = min(1.0_rp, max(0.0_rp, alpha))
        Yk_prof(1:nZ_CMC_chm) = rscal_inert_CMC_chm(1:nZ_CMC_chm,iclas+1_ip) + &
           alpha * (rscal_equil_CMC_chm(1:nZ_CMC_chm,iclas+1_ip) - rscal_inert_CMC_chm(1:nZ_CMC_chm,iclas+1_ip))
     end if

  end subroutine chm_get_cond_species_average_linear_CMC


  !!!subroutine chm_get_cond_species_average_HR_CMC(iclas,ipoin,Yk_avg,Zavg,Zvar,Yk_prof)
  !!!   !-----------------------------------------------------------------------
  !!!   !****f* Chemic/chm_get_cond_field_average_CMC
  !!!   ! NAME
  !!!   !    chm_get_cond_field_average_CMC
  !!!   ! DESCRIPTION
  !!!   !    Find conditional species profiles along mixture fraction from average
  !!!   !    values using profiles from homogeneous reactors. Do this operation for
  !!!   !    just one single species.
  !!!   ! USED BY
  !!!   !    chm_get_cond_fields_bc_Dir_linear_CMC
  !!!   !***
  !!!   !-----------------------------------------------------------------------

  !!!   use def_chemic,             only :  nZ_CMC_chm, rscal_inert_CMC_chm, &
  !!!                                       rscal_equil_CMC_chm, Yk_CMC_chm, &
  !!!                                       alpha_grid_CMC_chm, n_alpha_grid_CMC_chm, &
  !!!                                       inv_exp_alpha_grid_CMC_chm, &
  !!!                                       homog_reactors_CMC_chm, alpha_val_spec_CMC_chm

  !!!   implicit none
  !!!   integer(ip), intent(in)          :: iclas
  !!!   integer(ip), intent(in)          :: ipoin
  !!!   real(rp), intent(in)             :: Yk_avg
  !!!   real(rp), intent(in)             :: Zavg
  !!!   real(rp), intent(in)             :: Zvar
  !!!   real(rp), intent(out)            :: Yk_prof(nZ_CMC_chm)

  !!!   integer(ip)                      :: pos_tentat, pos, pos_int
  !!!   real(rp)                         :: aux, val_int(4), val_dif_aux(4)
  !!!   real(rp)                         :: val_int_l, val_int_r, ratio

  !!!   Yk_prof = 0.0_rp

  !!!   ! Find appropriate alpha from homogeneous reactors table
  !!!   aux = alpha_val_spec_CMC_chm(1_ip,ipoin,iclas) ** inv_exp_alpha_grid_CMC_chm
  !!!   pos_tentat = 1_ip + floor(aux*(n_alpha_grid_CMC_chm-1_ip))
  !!!   if (pos_tentat == n_alpha_grid_CMC_chm) then  ! alpha_val_spec_CMC_chm(1,ipoin,iclas)=1.0
  !!!      pos_tentat = pos_tentat - 1_ip
  !!!   end if

  !!!   if (pos_tentat == 1_ip) then  ! Extreme interval at the left
  !!!      ! Integrate profiles
  !!!      call chm_mxt_fr_integr_previous_steps_CMC(3_ip, &
  !!!          homog_reactors_CMC_chm(:,1:3,iclas), &
  !!!          Zavg, Zvar, val_int(1:3))
  !!!      val_dif_aux(1:3) = val_int(1:3) - Yk_avg
  !!!      if (val_dif_aux(1)*val_dif_aux(2) <= 0.0_rp) then
  !!!         pos     = 1_ip
  !!!         pos_int = 1_ip
  !!!      else if (val_dif_aux(2)*val_dif_aux(3) <= 0.0_rp) then
  !!!         pos     = 2_ip
  !!!         pos_int = 2_ip
  !!!      else
  !!!         call runend('CHEMIC CMC: not possible to obtain conditional profiles at boundaries')
  !!!      end if
  !!!
  !!!   else if (pos_tentat == n_alpha_grid_CMC_chm-1_ip) then  ! Extreme interval at the right
  !!!      ! Integrate profiles
  !!!      call chm_mxt_fr_integr_previous_steps_CMC(3_ip, &
  !!!          homog_reactors_CMC_chm(:,nZ_CMC_chm-2_ip:nZ_CMC_chm,iclas), &
  !!!          Zavg, Zvar, val_int(1:3))
  !!!      val_dif_aux(1:3) = val_int(1:3) - Yk_avg
  !!!      if (val_dif_aux(2)*val_dif_aux(3) <= 0.0_rp) then
  !!!         pos     = nZ_CMC_chm - 1_ip
  !!!         pos_int = 2_ip
  !!!      else if (val_dif_aux(1)*val_dif_aux(2) <= 0.0_rp) then
  !!!         pos     = nZ_CMC_chm - 2_ip
  !!!         pos_int = 1_ip
  !!!      else
  !!!         call runend('CHEMIC CMC: not possible to obtain conditional profiles at boundaries')
  !!!      end if

  !!!   else
  !!!      ! Integrate profiles
  !!!      call chm_mxt_fr_integr_previous_steps_CMC(4_ip, &
  !!!          homog_reactors_CMC_chm(:,pos_tentat-1_ip:pos_tentat+2_ip,iclas), &
  !!!          Zavg, Zvar, val_int)
  !!!      val_dif_aux = val_int - Yk_avg
  !!!      if (val_dif_aux(2)*val_dif_aux(3) <= 0.0_rp) then
  !!!         pos     = pos_tentat
  !!!         pos_int = 2_ip
  !!!      else
  !!!         if (val_dif_aux(1)*val_dif_aux(2) <= 0.0_rp) then
  !!!            pos     = pos_tentat - 1_ip
  !!!            pos_int = 1_ip
  !!!         else if (val_dif_aux(3)*val_dif_aux(4) <= 0.0_rp) then
  !!!            pos     = pos_tentat + 1_ip
  !!!            pos_int = 3_ip
  !!!         else
  !!!            call runend('CHEMIC CMC: not possible to obtain conditional profiles at boundaries')
  !!!         end if
  !!!      end if

  !!!   end if

  !!!   val_int_l = val_int(pos_int)
  !!!   val_int_r = val_int(pos_int+1_ip)

  !!!   if (val_int_l /= val_int_r) then
  !!!      ratio = (Yk_avg - val_int_l) / (val_int_r - val_int_l)

  !!!      alpha_val_spec_CMC_chm(1_ip,ipoin,iclas) = alpha_grid_CMC_chm(pos) + &
  !!!         (alpha_grid_CMC_chm(pos+1) - alpha_grid_CMC_chm(pos)) * ratio

  !!!      ! Conditional profile
  !!!      ratio = (alpha_val_spec_CMC_chm(1_ip,ipoin,iclas) - alpha_grid_CMC_chm(pos)) / &
  !!!                 (alpha_grid_CMC_chm(pos+1) - alpha_grid_CMC_chm(pos))
  !!!      Yk_prof = homog_reactors_CMC_chm(:,pos,iclas) + ratio * &
  !!!                    (homog_reactors_CMC_chm(:,pos+1_ip,iclas) - homog_reactors_CMC_chm(:,pos,iclas))

  !!!   else
  !!!      ! Take the profile from previous time step
  !!!      Yk_prof(1:nZ_CMC_chm) = Yk_CMC_chm(1:nZ_CMC_chm,ipoin,iclas)
  !!!   end if

  !!!end subroutine chm_get_cond_species_average_HR_CMC


  subroutine chm_get_cond_enthalpy_average_CMC(temp_uncond,Yk_avg,Zavg,Zvar,alpha,enthalp_prof)
     !-----------------------------------------------------------------------
     !****f* Chemic/chm_get_cond_enthalpy_average_CMC
     ! NAME
     !    chm_get_cond_enthalpy_average_CMC
     ! DESCRIPTION
     !    Find conditional enthalpy profile along mixture fraction from average
     !    temperature and conditional mass fractions. It can be applied to both
     !    the use of alfa with linear relationships and HRs.
     ! USED BY
     !    chm_compute_initial_fields_CMC
     !***
     !-----------------------------------------------------------------------

     use def_chemic,                 only :  nZ_CMC_chm, temp_inert_CMC_chm, &
                                             temp_equil_CMC_chm

     implicit none
     real(rp), parameter                  :: Tlim = 200.0_rp
     real(rp), intent(in)                 :: temp_uncond
     real(rp), intent(in)                 :: Yk_avg(nZ_CMC_chm,nclas_chm)
     real(rp), intent(in)                 :: Zavg
     real(rp), intent(in)                 :: Zvar
     real(rp), intent(in)                 :: alpha
     real(rp), intent(out)                :: enthalp_prof(nZ_CMC_chm)
     integer(ip)                          :: imixf
     real(rp)                             :: temp_prof(nZ_CMC_chm)
     real(rp)                             :: temp_adbiab_cond(nZ_CMC_chm,1)
     real(rp)                             :: temp_int(1)
     real(rp)                             :: factor_T
     real(rp)                             :: h_aux
     real(rp)                             :: dummy

     external                             :: runend

     dummy = 0.0_rp

     temp_adbiab_cond(1:nZ_CMC_chm,1) = temp_inert_CMC_chm(1:nZ_CMC_chm) + &
         alpha * (temp_equil_CMC_chm(1:nZ_CMC_chm) - temp_inert_CMC_chm(1:nZ_CMC_chm))


     ! Integrate temperature profile with adiabatic condition
     call chm_mxt_fr_integr_previous_steps_CMC(1_ip, &
            temp_adbiab_cond, Zavg, Zvar, temp_int)

     factor_T = temp_uncond / temp_int(1)

     temp_prof(1:nZ_CMC_chm) = factor_T * temp_adbiab_cond(1:nZ_CMC_chm,1)

     do imixf = 1, nZ_CMC_chm
        if (temp_prof(imixf) < Tlim) &
           call runend('CHEMIC OPERATIONS MOD. CMC: Conditional temperature lower than the minimum possible temperature')
     end do

     ! Get enthalpy
     do imixf = 1, nZ_CMC_chm
        call calc_h_cp_from_TY(Yk_avg(imixf,1:nclas_chm), temp_prof(imixf), &
                h_aux, dummy)
        enthalp_prof(imixf) = h_aux
     end do

  end subroutine chm_get_cond_enthalpy_average_CMC


  subroutine chm_compute_average_alpha_enthalpy_linear_CMC(Yk_cond, &
                      Zavg, Zvar, alpha_avg)
     !-----------------------------------------------------------------------
     !****f* Chemic/chm_compute_average_alpha_enthalpy_linear_CMC
     ! NAME
     !    chm_compute_average_alpha_enthalpy_linear_CMC
     ! DESCRIPTION
     !    Compute an average value for alpha from species for being used for
     !    enthalpy profile determination at Dirichlet BCs.
     ! USED BY
     !    chm_get_cond_fields_bc_Dir_CMC
     !***
     !-----------------------------------------------------------------------
     use def_chemic,           only :  rscal_inert_CMC_chm, rscal_equil_CMC_chm, &
                                       index_N2, nZ_CMC_chm, max_diff_Yk_inert_eq_CMC_chm

     implicit none
     real(rp), parameter            :: small = 5.0e-4_rp
     real(rp), intent(in)           :: Yk_cond(nZ_CMC_chm,nspec_chm)
     real(rp), intent(in)           :: Zavg
     real(rp), intent(in)           :: Zvar
     real(rp), intent(out)          :: alpha_avg

     integer(ip)                    :: iclas
     real(rp)                       :: inert_interm_eq_int(3_ip*nspec_chm)
     real(rp)                       :: Yk_sum_aux, diff_aux

     alpha_avg  = 0.0_rp
     Yk_sum_aux = 0.0_rp

     ! Get a parameter that weighes inert and equilibrium solutions
     call chm_mxt_fr_integr_previous_steps_CMC(3_ip*nspec_chm, &
          (/rscal_inert_CMC_chm(:,2:nspec_chm+1_ip), Yk_cond, &
          rscal_equil_CMC_chm(:,2:nspec_chm+1_ip)/), &
          Zavg, Zvar, inert_interm_eq_int(1:3_ip*nspec_chm))

     do iclas = 1, nspec_chm
        if (iclas /= index_N2) then  ! N2 is the same in inert and equilibrium conditions except
           ! for nitrogen oxides which are very low in mass. Better not to account for N2 in the
           ! following average since YN2 is very large and can lead to spurious results
           diff_aux = inert_interm_eq_int(iclas+2_ip*nspec_chm) - inert_interm_eq_int(iclas)

           if (abs(diff_aux) >= small * max_diff_Yk_inert_eq_CMC_chm(iclas)) then
              alpha_avg = alpha_avg + inert_interm_eq_int(iclas+nspec_chm) * &
                 (inert_interm_eq_int(iclas+nspec_chm) - inert_interm_eq_int(iclas)) / diff_aux
              Yk_sum_aux = Yk_sum_aux + inert_interm_eq_int(iclas+nspec_chm)
           end if

        end if
     end do

     alpha_avg = alpha_avg / Yk_sum_aux

  end subroutine chm_compute_average_alpha_enthalpy_linear_CMC


  !!!subroutine chm_compute_average_alpha_enthalpy_HR_CMC(ipoin, Yk_cond, alpha_avg)
  !!!   !-----------------------------------------------------------------------
  !!!   !****f* Chemic/chm_compute_average_alpha_enthalpy_HR_CMC
  !!!   ! NAME
  !!!   !    chm_compute_average_alpha_enthalpy_HR_CMC
  !!!   ! DESCRIPTION
  !!!   !    Compute an average value for alpha from species for being used for
  !!!   !    enthalpy profile determination at Dirichlet BCs.
  !!!   ! USED BY
  !!!   !    chm_get_cond_fields_bc_Dir_CMC
  !!!   !***
  !!!   !-----------------------------------------------------------------------
  !!!   use def_chemic,           only :  index_N2, bvess_ufield_CMC_chm, &
  !!!                                     alpha_val_spec_CMC_chm, nZ_CMC_chm

  !!!   implicit none
  !!!   integer(ip), intent(in)        :: ipoin
  !!!   real(rp), intent(in)           :: Yk_cond(nZ_CMC_chm,nclas_chm)
  !!!   real(rp), intent(out)          :: alpha_avg

  !!!   integer(ip)                    :: iclas
  !!!   real(rp)                       :: Yk_sum_aux

  !!!   alpha_avg  = 0.0_rp
  !!!   Yk_sum_aux = 0.0_rp

  !!!   do iclas = 1, nclas_chm
  !!!      if (iclas /= index_N2) then  ! N2 is the same in inert and equilibrium conditions except
  !!!         ! for nitrogen oxides which are very low in mass. Better not to account for N2 in the
  !!!         ! following average since YN2 is very large and can lead to spurious results

  !!!         alpha_avg = alpha_avg + bvess_ufield_CMC_chm(ipoin,iclas) * alpha_val_spec_CMC_chm(1_ip,ipoin,iclas)
  !!!         Yk_sum_aux = Yk_sum_aux + bvess_ufield_CMC_chm(ipoin,iclas)

  !!!      end if
  !!!   end do

  !!!   alpha_avg = alpha_avg / Yk_sum_aux

  !!!end subroutine chm_compute_average_alpha_enthalpy_HR_CMC


  subroutine chm_get_cond_fields_bc_Dir_CMC
     !-----------------------------------------------------------------------
     !****f* Chemic/chm_get_cond_fields_bc_Dir_CMC
     ! NAME
     !    chm_get_cond_fields_bc_Dir_CMC
     ! DESCRIPTION
     !    Compute conditional fields for points that belong to the boundary
     !    with Dirichlet boundary condition.
     ! USED BY
     !    chm_updbcs
     !***
     !-----------------------------------------------------------------------

     use def_chemic,                only :  kfl_bc_type_spec_CMC_chm, kfl_fixno_chm, &
                                            bvess_ufield_CMC_chm, Zavg_CMC_chm, Zvar_CMC_chm, &
                                            bvess_CMC_chm, Yk_CMC_chm, enthalp_CMC_chm, &
                                            kfl_solve_enth_CMC_chm, nZ_CMC_chm, nvar_CMC_chm, &
                                            kfl_bc_init_method_CMC_chm, &
                                            kfl_bc_alpha_CMC_chm, rscal_inert_CMC_chm, &
                                            rscal_equil_CMC_chm, kfl_soot_chm

     implicit none
     integer(ip)                         :: ipoin, iclas, imixf, kfl_dir
     real(rp)                            :: Yk_sum, alpha_avg
     real(rp)                            :: Yk_cond_bc_Dir(nZ_CMC_chm,nclas_chm)
     real(rp)                            :: Yk_poin(nZ_CMC_chm,nclas_chm)
     real(rp)                            :: aux_dif(nspec_chm)

     if (kfl_bc_alpha_CMC_chm == 1_ip) then
        if (kfl_bc_init_method_CMC_chm == 1_ip) then
           do ipoin = 1,npoin
              if ( kfl_fixno_chm(1,ipoin) > 0_ip ) then
                 do imixf = 2, nZ_CMC_chm-1_ip
                    aux_dif = rscal_equil_CMC_chm(imixf,2:nspec_chm+1_ip) - rscal_inert_CMC_chm(imixf,2:nspec_chm+1_ip)
                    bvess_CMC_chm(imixf,ipoin,1:nspec_chm) = rscal_inert_CMC_chm(imixf,2:nspec_chm+1_ip) + &
                       bvess_ufield_CMC_chm(ipoin,1) * aux_dif
                    Yk_CMC_chm(imixf,ipoin,1:nspec_chm) = bvess_CMC_chm(imixf,ipoin,1:nspec_chm)
                 end do
                 ! Soot mass fractions set to 0 for Dirichlet
                 if ( kfl_soot_chm > 0 ) then
                    bvess_CMC_chm(:,ipoin,nspec_chm+1:nclas_chm) = 0.0_rp
                    Yk_CMC_chm(:,ipoin,nspec_chm+1:nclas_chm)    = 0.0_rp
                 end if
              end if
           end do
        end if

     else
        do ipoin = 1,npoin

           Yk_cond_bc_Dir = 0.0_rp

           if (kfl_bc_type_spec_CMC_chm(ipoin) == 1_ip .and. kfl_bc_init_method_CMC_chm == 1_ip) then
             !BCs for all species of the same type and linear relationships
              if( kfl_fixno_chm(1,ipoin) > 0_ip ) then
                 call chm_get_cond_species_all_average_linear_CMC(ipoin,bvess_ufield_CMC_chm(ipoin,1:nspec_chm), &
                        Zavg_CMC_chm(ipoin), Zvar_CMC_chm(ipoin), Yk_cond_bc_Dir(1:nZ_CMC_chm,1:nspec_chm))
                 do imixf = 1,nZ_CMC_chm
                    do iclas = 1, nspec_chm
                       bvess_CMC_chm(imixf,ipoin,iclas) = Yk_cond_bc_Dir(imixf,iclas)
                       Yk_CMC_chm(imixf,ipoin,iclas)    = bvess_CMC_chm(imixf,ipoin,iclas)
                    end do
                 end do
                 ! Soot mass fractions set to 0 for Dirichlet
                 if ( kfl_soot_chm > 0 ) then
                    bvess_CMC_chm(:,ipoin,nspec_chm+1:nclas_chm) = 0.0_rp
                    Yk_CMC_chm(:,ipoin,nspec_chm+1:nclas_chm)    = 0.0_rp
                 end if
              end if

           else  ! Probably never used
              kfl_dir = 0_ip
              do iclas = 1, nspec_chm
                 if( kfl_fixno_chm(iclas,ipoin) > 0_ip ) then
                    if (kfl_bc_init_method_CMC_chm == 1_ip) then
                       call chm_get_cond_species_average_linear_CMC(iclas,ipoin,bvess_ufield_CMC_chm(ipoin,iclas), &
                               Zavg_CMC_chm(ipoin), Zvar_CMC_chm(ipoin), Yk_cond_bc_Dir(1:nZ_CMC_chm,iclas))
                    !!!else if (kfl_bc_init_method_CMC_chm == 2_ip) then
                    !!!   call chm_get_cond_species_average_HR_CMC(iclas,ipoin,bvess_ufield_CMC_chm(ipoin,iclas), &
                    !!!           Zavg_CMC_chm(ipoin), Zvar_CMC_chm(ipoin), Yk_cond_bc_Dir(1:nZ_CMC_chm,iclas))
                    end if
                    kfl_dir = 1_ip
                 end if
              end do

              if (kfl_dir == 1_ip) then
                 do iclas = 1, nclas_chm
                    if( kfl_fixno_chm(iclas,ipoin) == 0_ip ) then
                       ! In case for some species we have Dirichlet and for other species Neumann
                       ! (kfl_dir). For Neumann assign values from previous time step
                       Yk_cond_bc_Dir(1:nZ_CMC_chm,iclas) = Yk_CMC_chm(1:nZ_CMC_chm,ipoin,iclas)
                    end if
                 end do

                 ! Normalize in case mass fractions do not sum 1
                 do imixf = 2, nZ_CMC_chm-1_ip
                    Yk_sum = 0.0_rp
                    do iclas = 1, nclas_chm
                       Yk_sum = Yk_sum + Yk_cond_bc_Dir(imixf,iclas)
                    end do
                    if (Yk_sum /= 0.0_rp) then
                       Yk_cond_bc_Dir(imixf,1:nclas_chm) = Yk_cond_bc_Dir(imixf,1:nclas_chm) / Yk_sum
                    end if
                 end do

                 ! Assign to bvess_CMC_chm
                 do imixf = 1,nZ_CMC_chm
                    do iclas = 1, nclas_chm
                       bvess_CMC_chm(imixf,ipoin,iclas) = Yk_cond_bc_Dir(imixf,iclas)
                       Yk_CMC_chm(imixf,ipoin,iclas)    = bvess_CMC_chm(imixf,ipoin,iclas)
                    end do
                 end do
              end if

           end if
        end do

     end if

     if (kfl_solve_enth_CMC_chm /= 0_ip) then
        do ipoin = 1, npoin
           if( kfl_fixno_chm(nvar_CMC_chm,ipoin) > 0_ip ) then

              do imixf = 1, nZ_CMC_chm
                 do iclas = 1, nclas_chm
                    Yk_poin(imixf,iclas) = Yk_CMC_chm(imixf,ipoin,iclas)
                 end do
              end do

              if (kfl_bc_alpha_CMC_chm == 1_ip) then
                 call chm_get_cond_enthalpy_average_CMC(bvess_ufield_CMC_chm(ipoin,2), &
                          Yk_poin(1:nZ_CMC_chm,1:nclas_chm), Zavg_CMC_chm(ipoin), &
                          Zvar_CMC_chm(ipoin), bvess_ufield_CMC_chm(ipoin,1), &
                          bvess_CMC_chm(1:nZ_CMC_chm,ipoin,nvar_CMC_chm))

              else
                 if (kfl_bc_init_method_CMC_chm == 1_ip) then
                    call chm_compute_average_alpha_enthalpy_linear_CMC(Yk_poin, &
                            Zavg_CMC_chm(ipoin), Zvar_CMC_chm(ipoin), alpha_avg)
                 !!!else if (kfl_bc_init_method_CMC_chm == 2_ip) then
                 !!!   call chm_compute_average_alpha_enthalpy_HR_CMC(ipoin, Yk_poin, alpha_avg)
                 end if

                 call chm_get_cond_enthalpy_average_CMC(bvess_ufield_CMC_chm(ipoin,nvar_CMC_chm), &
                          Yk_poin(1:nZ_CMC_chm,1:nclas_chm), Zavg_CMC_chm(ipoin), &
                          Zvar_CMC_chm(ipoin), alpha_avg, &
                          bvess_CMC_chm(1:nZ_CMC_chm,ipoin,nvar_CMC_chm))
              end if

              do imixf = 1, nZ_CMC_chm
                 enthalp_CMC_chm(imixf,ipoin) = bvess_CMC_chm(imixf,ipoin,nvar_CMC_chm)
              end do
           end if
        end do
     end if

  end subroutine chm_get_cond_fields_bc_Dir_CMC


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! S U B R O U T I N E S   F O R   C O M M U N I C A T I O N !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine chm_find_master_CMC(who_I_am, rank_master)
     !-----------------------------------------------------------------------
     !****f* Chemic/chm_find_master_CMC
     ! NAME
     !    chm_find_master_CMC
     ! DESCRIPTION
     !    Find rank for master
     ! USED BY
     !    chm_plugin
     !-----------------------------------------------------------------------
     use def_master,             only :  kfl_paral
     use mod_parall,             only :  PAR_MY_WORLD_RANK
     use mod_communications,     only :  PAR_SUM
     use def_chemic,             only :  kfl_solve_cond_CMC_chm

     implicit none
     integer(ip), intent(in)          :: who_I_am
     integer(ip), intent(inout)       :: rank_master

     rank_master = 0_ip
     if (kfl_paral == 0 .and. kfl_solve_cond_CMC_chm == who_I_am) then
        rank_master = PAR_MY_WORLD_RANK
     end if
     call PAR_SUM(rank_master,'IN THE UNIVERSE')

  end subroutine chm_find_master_CMC


  subroutine chm_transfer_info_CMC_to_CFD(rank_master_CMC)
     !-----------------------------------------------------------------------
     !****f* Chemic/chm_transfer_info_CMC_to_CFD
     ! NAME
     !    chm_transfer_info_CMC_to_CFD
     ! DESCRIPTION
     !    Send information for the transfer from CMC to CFD
     ! USED BY
     !    chm_parall
     !-----------------------------------------------------------------------
     use def_chemic,            only :  nZ_CMC_chm, Z_CMC_chm, &
                                        kfl_transfer_condField_CMC_chm, &
                                        nspec_transf_CMC_chm, transf_spec_CMC_chm, &
                                        kfl_solve_cond_CMC_chm, &
                                        nvar_trans_CMC_chm, transf_entha_CMC_chm, &
                                        nspec_uncond_write_CMC_chm, write_uncond_spec_CMC_chm
     use mod_communications,    only :  PAR_BROADCAST

     implicit none
     integer(ip), intent(in)         :: rank_master_CMC
     integer(ip)                     :: ii

     external                        :: chm_memphy
     external                        :: chm_memnut
     external                        :: chm_memous

     call PAR_BROADCAST(kfl_transfer_condField_CMC_chm, 'IN THE UNIVERSE', root_rank = rank_master_CMC)
     if (kfl_transfer_condField_CMC_chm == 1_ip) then

        ! Transfer information about the mixture fraction vector
        call PAR_BROADCAST(nZ_CMC_chm, 'IN THE UNIVERSE', root_rank = rank_master_CMC)
        if (kfl_solve_cond_CMC_chm == 0_ip) call chm_memphy(3_ip)
        do ii = 1, nZ_CMC_chm
           call PAR_BROADCAST(Z_CMC_chm(ii), 'IN THE UNIVERSE', root_rank = rank_master_CMC)
        end do

        ! Transfer species
        call PAR_BROADCAST(nspec_transf_CMC_chm, 'IN THE UNIVERSE', root_rank = rank_master_CMC)
        if (nspec_transf_CMC_chm /= 0_ip) then
           if (kfl_solve_cond_CMC_chm == 0_ip) then
              call chm_memnut(2_ip)
              nspec_uncond_write_CMC_chm = nspec_transf_CMC_chm
              call chm_memous(1_ip)
              write_uncond_spec_CMC_chm = (/( ii, ii=1_ip, nspec_uncond_write_CMC_chm )/)
           end if
           if (kfl_solve_cond_CMC_chm == 1_ip .and. INOTMASTER)  call chm_memnut(2_ip)
           do ii = 1, nspec_transf_CMC_chm
              call PAR_BROADCAST(transf_spec_CMC_chm(ii), 'IN THE UNIVERSE', root_rank = rank_master_CMC)
           end do
        end if

        ! Transfer enthalpy
        call PAR_BROADCAST(transf_entha_CMC_chm, 'IN THE UNIVERSE', root_rank = rank_master_CMC)

        nvar_trans_CMC_chm = 6_ip + nspec_transf_CMC_chm + transf_entha_CMC_chm
     end if

  end subroutine chm_transfer_info_CMC_to_CFD


  subroutine chm_transfer_info_CFD_to_CMC(rank_master_CFD)
     !-----------------------------------------------------------------------
     !****f* Chemic/chm_transfer_info_CFD_to_CMC
     ! NAME
     !    chm_transfer_info_CFD_to_CMC
     ! DESCRIPTION
     !    Send information for the transfer from CFD to CMC
     ! USED BY
     !    chm_parall
     !-----------------------------------------------------------------------
     use def_chemic,             only :  nwrite_sol_CFD_CMC_chm, kfl_solve_cond_CMC_chm, &
                                         kfl_transfer_condField_CMC_chm, dt_CFD_CMC_chm, &
                                         freq_CFD_coup_CMC_chm, last_time_step_CFD_CMC_chm, &
                                         last_time_step_CMC_CMC_chm, kfl_av_species_CMC_chm, &
                                         kfl_av_enthalp_CMC_chm
     use mod_communications,     only :  PAR_BROADCAST
     use def_master,             only :  postp, dtinv, IMASTER, ittim
     use def_coupli,             only :  coupling_type
     use mod_output_postprocess, only :  output_postprocess_check_variable_postprocess

     implicit none
     integer(ip), intent(in)         :: rank_master_CFD
     integer(ip)                     :: icoup, aux2
     real(rp)                        :: aux

     external                        :: runend

     nwrite_sol_CFD_CMC_chm = 0_ip

     if (kfl_solve_cond_CMC_chm == 0_ip) then
        dt_CFD_CMC_chm = 1.0_rp / dtinv
        last_time_step_CFD_CMC_chm = ittim
        ! Averages
        if ( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVYUN') )  &
           kfl_av_species_CMC_chm = 1_ip
        if ( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVHUN') )  &
           kfl_av_enthalp_CMC_chm = 1_ip
     end if

     call PAR_BROADCAST(dt_CFD_CMC_chm, 'IN THE UNIVERSE', root_rank = rank_master_CFD)
     call PAR_BROADCAST(last_time_step_CFD_CMC_chm, 'IN THE UNIVERSE', root_rank = rank_master_CFD)
     call PAR_BROADCAST(kfl_av_species_CMC_chm, 'IN THE UNIVERSE', root_rank = rank_master_CFD)
     call PAR_BROADCAST(kfl_av_enthalp_CMC_chm, 'IN THE UNIVERSE', root_rank = rank_master_CFD)

     if (kfl_transfer_condField_CMC_chm == 1_ip) then
        if (kfl_solve_cond_CMC_chm == 0_ip) then
            if (postp(1) % npp_stepi(arrays_number('MASUN'),0) /= 0) then  ! MASUN
               nwrite_sol_CFD_CMC_chm = postp(1) % npp_stepi(arrays_number('MASUN'),0)
            else if (postp(1) % npp_stepi(arrays_number('HRR  '),0) /= 0) then  ! HRR
               nwrite_sol_CFD_CMC_chm = postp(1) % npp_stepi(arrays_number('HRR  '),0)
            end if
        end if

        ! Checkings for the coupling
        if (kfl_solve_cond_CMC_chm == 1_ip) then
           if (IMASTER) then
              do icoup = 1,2
                 if (coupling_type(icoup) % variable == 'CM2CF') then
                    if (coupling_type(icoup) % frequ_send /= 1_ip) call runend('Coupling: In the coupling from CMC to CFD, CMC has&
                        & to send at each time step')
                 else if (coupling_type(icoup) % variable == 'CF2CM') then
                    if (coupling_type(icoup) % frequ_recv /= 1_ip) call runend('Coupling: In the coupling from CFD to CMC, CMC has&
                        & to receive at each time step')
                    freq_CFD_coup_CMC_chm = coupling_type(icoup) % frequ_send
                 end if
              end do
              if (coupling_type(1) % frequ_send /= coupling_type(2) % frequ_recv) call runend('Coupling: mismatch between the&
                  & frequencies from CMC and CFD')
              if (coupling_type(2) % frequ_send /= coupling_type(1) % frequ_recv) call runend('Coupling: mismatch between the&
                  & frequencies from CMC and CFD')

              ! Check time step
              aux = 1.0_rp / (dt_CFD_CMC_chm*dtinv)
              if ( abs(aux - real(freq_CFD_coup_CMC_chm,rp)) > 1.0e-3_rp ) call runend('Coupling: The relation between the time&
                  & steps does not match with the frequency to exchange information')
           end if
           last_time_step_CMC_CMC_chm = ittim
        end if

        do icoup = 1,2
           if (coupling_type(icoup) % variable == 'CF2CM')  freq_CFD_coup_CMC_chm = coupling_type(icoup) % frequ_send
        end do

        if (kfl_solve_cond_CMC_chm == 0_ip .and. IMASTER) then
           aux2 = nwrite_sol_CFD_CMC_chm / freq_CFD_coup_CMC_chm
           if (aux2 == 0_ip)  aux2 = 1_ip
           nwrite_sol_CFD_CMC_chm = aux2 * freq_CFD_coup_CMC_chm
        end if
        call PAR_BROADCAST(nwrite_sol_CFD_CMC_chm, 'IN THE UNIVERSE', root_rank = rank_master_CFD)

     else
        if (kfl_solve_cond_CMC_chm == 1_ip) then
           if (IMASTER) then
              aux = dt_CFD_CMC_chm*dtinv
              if ( abs(aux - 1.0_rp) > 1.0e-3_rp ) call runend('Coupling: Time steps in CFD and CMC are not the same')
           end if
        end if
     end if

     !! NOTA: quizás falta alguna comprobación del número de pasos del restart

  end subroutine chm_transfer_info_CFD_to_CMC


  subroutine chm_write_info_transfer_CMC
     !-----------------------------------------------------------------------
     !****f* Chemic/chm_write_info_transfer_CMC
     ! NAME
     !    chm_write_info_transfer_CMC
     ! DESCRIPTION
     !    Write useful information about the transfer
     ! USED BY
     !    chm_parall
     !-----------------------------------------------------------------------
     use def_master,             only : momod, modul, IMASTER
     use def_chemic,             only : nvar_trans_CMC_chm, nspec_transf_CMC_chm, &
                                        transf_spec_CMC_chm, transf_entha_CMC_chm, &
                                        nZ_CMC_chm, Z_CMC_chm, nwrite_sol_CFD_CMC_chm, &
                                        kfl_transfer_condField_CMC_chm, kfl_solve_cond_CMC_chm, &
                                        kfl_av_species_CMC_chm, kfl_av_enthalp_CMC_chm
     use mod_output_postprocess, only :  output_postprocess_check_variable_postprocess

     implicit none
     external                        :: runend
     integer(ip)                     :: ii

     if (IMASTER .and. kfl_transfer_condField_CMC_chm == 0_ip .and. kfl_solve_cond_CMC_chm == 0_ip) then
        if (     output_postprocess_check_variable_postprocess(VARIABLE_NAME='MASUN')  &
            .or. output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENTUN')  &
            .or. output_postprocess_check_variable_postprocess(VARIABLE_NAME='HRR  ')  &
            .or. output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVHRR')  &
            .or. output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVYUN')  &
            .or. output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVHUN')) &
           call runend('CHEMIC TRANSFER: some of the postprocessing variables cannot be written since the transfer from CMC to CFD&
              & is not activated')
     end if

     if (IMASTER .and. kfl_transfer_condField_CMC_chm == 1_ip) then

        if ((     output_postprocess_check_variable_postprocess(VARIABLE_NAME='MASUN') &
            .or. output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVYUN')) &
            .and. nspec_transf_CMC_chm == 0_ip) &
           call runend('CHEMIC TRANSFER: postprocessing for species cannot be written since the transfer from CMC to CFD for&
              & species is not activated')

        if ((     output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENTUN') &
            .or. output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVHUN')) &
            .and. transf_entha_CMC_chm == 0_ip) &
           call runend('CHEMIC TRANSFER: postprocessing for enthalpy cannot be written since the transfer from CMC to CFD for&
              & enthalpy is not activated')

        write(momod(modul) % lun_outpu,*)
        write(momod(modul) % lun_outpu,*) '---------------------------------------------------------'
        write(momod(modul) % lun_outpu,*) ' T R A N S F E R   I N F O R M A T I O N   S E C T I O N'
        write(momod(modul) % lun_outpu,*) '---------------------------------------------------------'
        write(momod(modul) % lun_outpu,*) 'Total number of fields to be transferred: ', nvar_trans_CMC_chm
        write(momod(modul) % lun_outpu,*)
        if (nspec_transf_CMC_chm > 0_ip) then
           write(momod(modul) % lun_outpu,*) 'Total number of species to be transferred: ', nspec_transf_CMC_chm
           do ii = 1, nspec_transf_CMC_chm
              write(momod(modul) % lun_outpu,*) ii, ': ', transf_spec_CMC_chm(ii)
           end do
           write(momod(modul) % lun_outpu,*)
        end if
        write(momod(modul) % lun_outpu,*)
        if (transf_entha_CMC_chm > 0_ip) write(momod(modul) % lun_outpu,*) 'Enthalpy transferred'
        write(momod(modul) % lun_outpu,*)
        write(momod(modul) % lun_outpu,*) 'Mixture fraction vector: number of entries: ', nZ_CMC_chm
        do ii = 1, nZ_CMC_chm
           write(momod(modul) % lun_outpu,*) ii, ': ', Z_CMC_chm(ii)
        end do
        write(momod(modul) % lun_outpu,*)
        write(momod(modul) % lun_outpu,*) 'Fequency for CFD post-processing: ', nwrite_sol_CFD_CMC_chm
        write(momod(modul) % lun_outpu,*)
        write(momod(modul) % lun_outpu,*)
        if (kfl_av_species_CMC_chm == 1_ip) then
           write(momod(modul) % lun_outpu,*)
           write(momod(modul) % lun_outpu,*) 'Averages for species will be computed in CFD'
           write(momod(modul) % lun_outpu,*) 'CAUTION: conditional mass fractions for species will be communicated from CMC to CFD&
             & at each time step'
           write(momod(modul) % lun_outpu,*)
        end if
        if (kfl_av_enthalp_CMC_chm == 1_ip) then
           write(momod(modul) % lun_outpu,*)
           write(momod(modul) % lun_outpu,*) 'Average for enthalpy will be computed in CFD'
           write(momod(modul) % lun_outpu,*)
        end if
     end if

  end subroutine chm_write_info_transfer_CMC


  subroutine chm_data_CFD_to_CMC_domain_CMC(icoup)
     !-----------------------------------------------------------------------
     !****f* Chemic/chm_data_CFD_to_CMC_domain_CMC
     ! NAME
     !    chm_data_CFD_to_CMC_domain_CMC
     ! DESCRIPTION
     !    Send information from CFD to CMC
     !    Variables: Z, Zv, grad_Z, nu_t, u, N
     ! USED BY
     !    chm_plugin
     !-----------------------------------------------------------------------

     use def_master,         only : conce, advec
     use def_kermod,         only : turmu_ker
     use def_coupli,         only : coupling_type
     use mod_couplings,      only : COU_INTERPOLATE_NODAL_VALUES
     use mod_interpolation,  only : COU_GET_INTERPOLATE_POINTS_VALUES
     use mod_ker_proper,     only : ker_proper
     use mod_gradie,         only : gradie
     use def_chemic,         only : kfl_model_chm, kfl_solve_cond_CMC_chm, &
                                    Zavg_CMC_chm, Zvar_CMC_chm, Xtot_CMC_chm, &
                                    grad_Zavg_CMC_chm, grad_Zavg_CFD_chm, &
                                    turb_kin_visc_CMC_chm, veloc_CMC_chm, &
                                    xZr_chm, xZs_chm, kfl_start_CMC_chm, &
                                    send_CFD_to_CMC_chm, rec_CMC_from_CFD_chm, &
                                    turb_kin_visco_CFD_chm, kfl_mesh_interp_CMC_chm, &
                                    kfl_avg_cond_CMC_chm, aux_interp_fields_CMC_chm

     implicit none
     integer(ip), intent(in) :: icoup
     real(rp), pointer       :: dummr2(:,:)

     integer(ip)             :: ipoin, dummi, irow

     nullify(dummr2)

     !
     ! Send information from CFD to CMC: mixing variables and velocity
     !
     if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 0) then

        if (INOTMASTER) then
           if (kfl_start_CMC_chm == 1) then
              allocate(grad_Zavg_CFD_chm(ndime,npoin))
              allocate(turb_kin_visco_CFD_chm(npoin))
           end if

           call gradie(conce(1:npoin,3,1),grad_Zavg_CFD_chm)

           do ipoin = 1, npoin
              send_CFD_to_CMC_chm(1,ipoin) = conce(ipoin,3,1)
              send_CFD_to_CMC_chm(2,ipoin) = conce(ipoin,4,1)
              send_CFD_to_CMC_chm(3,ipoin) = xZr_chm(ipoin) + xZs_chm(ipoin)
              send_CFD_to_CMC_chm(5:4+ndime,ipoin) = grad_Zavg_CFD_chm(1:ndime,ipoin)
              send_CFD_to_CMC_chm(5+ndime:4+2*ndime,ipoin) = advec(1:ndime,ipoin,1)
           end do


           ! Turbulent viscosity
           if (turmu_ker % kfl_exist /= 0_ip) then  ! Turbulent case
              call ker_proper('TURBU','NPOIN',dummi,dummi,turb_kin_visco_CFD_chm) ! TURBU returns turbulent kinematic viscosity
              do ipoin = 1, npoin
                 send_CFD_to_CMC_chm(4,ipoin) = turb_kin_visco_CFD_chm(ipoin)
              end do

           else ! Laminar case
              do ipoin = 1, npoin
                 send_CFD_to_CMC_chm(4,ipoin) = 0.0_rp
              end do
           end if
        end if

        if (kfl_mesh_interp_CMC_chm == 1_ip) then
           call COU_GET_INTERPOLATE_POINTS_VALUES(send_CFD_to_CMC_chm,dummr2,coupling_type(icoup))
        else
           call COU_INTERPOLATE_NODAL_VALUES(icoup,4_ip+2_ip*ndime,dummr2,send_CFD_to_CMC_chm)
        end if

     end if

     !
     ! Receive information by CMC from CFD: mixing variables and velocity
     !
     if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1) then

        if (kfl_mesh_interp_CMC_chm == 1_ip) then
           call COU_GET_INTERPOLATE_POINTS_VALUES(dummr2,rec_CMC_from_CFD_chm,coupling_type(icoup))

           if (INOTMASTER) then

              if (kfl_avg_cond_CMC_chm == 1_ip) then
                 call chm_average_field_mesh_CMC(2_ip, &
                        (/rec_CMC_from_CFD_chm(1,:), rec_CMC_from_CFD_chm(2,:)/), &
                        aux_interp_fields_CMC_chm(:,1:2))

                 do ipoin = 1, npoin
                    Zavg_CMC_chm(ipoin) = aux_interp_fields_CMC_chm(ipoin,1)
                    Zvar_CMC_chm(ipoin) = aux_interp_fields_CMC_chm(ipoin,2)
                 end do

              else
                 ! rec_CMC_from_CFD_chm has to be passed as
                 ! (/(rec_CMC_from_CFD_chm(irow,:), irow = 1_ip, 4_ip+2_ip*ndime)/)
                 ! Just giving rec_CMC_from_CFD_chm the subroutine will not work
                 ! as desired.

                 call chm_average_field_mesh_CMC(4_ip+2_ip*ndime, &
                         (/(rec_CMC_from_CFD_chm(irow,:), irow = 1_ip, 4_ip+2_ip*ndime)/), &
                         aux_interp_fields_CMC_chm)

                 do ipoin = 1, npoin
                    Zavg_CMC_chm(ipoin) = aux_interp_fields_CMC_chm(ipoin,1)
                    Zvar_CMC_chm(ipoin) = aux_interp_fields_CMC_chm(ipoin,2)
                    Xtot_CMC_chm(ipoin) = aux_interp_fields_CMC_chm(ipoin,3)
                    turb_kin_visc_CMC_chm(ipoin) = aux_interp_fields_CMC_chm(ipoin,4)
                    grad_Zavg_CMC_chm(1:ndime,ipoin) = aux_interp_fields_CMC_chm(ipoin,5:4+ndime)
                    veloc_CMC_chm(1:ndime,ipoin) = aux_interp_fields_CMC_chm(ipoin,5+ndime:4+2*ndime)
                 end do

              end if
           end if

        else
           call COU_INTERPOLATE_NODAL_VALUES(icoup,4_ip+2_ip*ndime,rec_CMC_from_CFD_chm)

           if (INOTMASTER) then
              do ipoin = 1, npoin
                 Zavg_CMC_chm(ipoin) = rec_CMC_from_CFD_chm(1,ipoin)
                 Zvar_CMC_chm(ipoin) = rec_CMC_from_CFD_chm(2,ipoin)
                 Xtot_CMC_chm(ipoin) = rec_CMC_from_CFD_chm(3,ipoin)
                 turb_kin_visc_CMC_chm(ipoin) = rec_CMC_from_CFD_chm(4,ipoin)
                 grad_Zavg_CMC_chm(1:ndime,ipoin) = rec_CMC_from_CFD_chm(5:4+ndime,ipoin)
                 veloc_CMC_chm(1:ndime,ipoin) = rec_CMC_from_CFD_chm(5+ndime:4+2*ndime,ipoin)
              end do
           end if
        end if

        call clipping_var_from_CFD_CMC   ! Clipping of variables from CFD

     end if

  end subroutine chm_data_CFD_to_CMC_domain_CMC


  subroutine chm_data_CMC_to_CFD_domain_CMC(icoup)
     !-----------------------------------------------------------------------
     !****f* Chemic/chm_data_CMC_to_CFD_domain_CMC
     ! NAME
     !    chm_data_CMC_to_CFD_domain_CMC
     ! DESCRIPTION
     !    Send information from CMC to CFD
     !    Variables: density, dyncamic laminar viscosity, conductivity,
     !    specific heat, temperature
     ! USED BY
     !    chm_plugin
     !-----------------------------------------------------------------------

     use def_master,          only : ittim
     use def_master,          only : densi, conce
     use def_master,          only : visco, condk, sphek, tempe
     use mod_couplings,       only : COU_INTERPOLATE_NODAL_VALUES
     use def_chemic,          only : freq_CFD_coup_CMC_chm, kfl_model_chm, kfl_post_gp_CMC_chm, kfl_solve_cond_CMC_chm,&
                                     kfl_solve_enth_CMC_chm, kfl_split_chm, kfl_start_CMC_chm, kfl_transfer_condField_CMC_chm,&
                                     last_time_step_CFD_CMC_chm, last_time_step_CMC_CMC_chm, nspec_transf_CMC_chm,&
                                     nvar_trans_CMC_chm, nwrite_sol_CFD_CMC_chm, nZ_CMC_chm, rec_CFD_from_CMC_post_chm,&
                                     transf_entha_CMC_chm, send_CMC_to_CFD_chm, temp_CMC_chm, hrr_mass_cond_CMC_chm,&
                                     transf_spec_CMC_chm, enthalp_CMC_chm, densi_int_CMC_chm, visco_lam_int_CMC_chm,&
                                     condu_int_CMC_chm, sphea_int_CMC_chm, temp_int_CMC_chm, Yk_int_CMC_chm, enthalp_int_CMC_chm,&
                                     send_CMC_to_CFD_post_chm, Yk_CMC_chm, rec_CFD_from_CMC_chm, rec_CFD_old_from_CMC_chm, &
                                     kfl_av_enthalp_CMC_chm, kfl_av_species_CMC_chm
     use mod_ker_updpro,      only : ker_updpro
     use def_kintyp,          only : lg

     implicit none
     integer(ip), intent(in)      :: icoup
     real(rp), pointer            :: dummr2(:,:,:), aux_var_uncond(:), aux_var_cond(:,:)
     logical(lg)                  :: cond_time_step

     integer(ip)                  :: ipoin, imixf, iclas, nvar, aux_t

     nullify(dummr2)

     !
     ! Send information from CMC to CFD: density, laminar viscosity,
     ! conductivity, cp and temperature
     !
     if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1) then

        if (kfl_start_CMC_chm == 1) then
           if( INOTMASTER ) then
              if (kfl_split_chm > 0_ip)   call chm_heatRelease_field_CMC   ! Compute heat release at nodes
           end if
           if (kfl_post_gp_CMC_chm == 0_ip)   call chm_integrate_flow_var_points_CMC
           call chm_integrate_flow_var_gauss_CMC
           call chm_nodal_projection_CMC
        end if

        ! Transference of information
        if (kfl_transfer_condField_CMC_chm == 1_ip) then
           do imixf = 1, nZ_CMC_chm
              call chm_nodal_projection_fundamental_CMC(imixf)
              do ipoin = 1, npoin
                 send_CMC_to_CFD_chm(1,imixf,ipoin) = densi(ipoin,1)  ! CAUTION: From CMC the specific volume is sent
                 send_CMC_to_CFD_chm(2,imixf,ipoin) = visco(ipoin,1)
                 send_CMC_to_CFD_chm(3,imixf,ipoin) = condk(ipoin,1)
                 send_CMC_to_CFD_chm(4,imixf,ipoin) = sphek(ipoin,1)
                 send_CMC_to_CFD_chm(5,imixf,ipoin) = temp_CMC_chm(imixf,ipoin)
                 send_CMC_to_CFD_chm(6,imixf,ipoin) = hrr_mass_cond_CMC_chm(imixf,ipoin)
              end do
           end do

           call COU_INTERPOLATE_NODAL_VALUES(icoup,6_ip*nZ_CMC_chm,dummr2,send_CMC_to_CFD_chm)

           ! Data only for post-processing
           if (nvar_trans_CMC_chm > 6_ip) then
              aux_t = last_time_step_CFD_CMC_chm + (ittim - last_time_step_CMC_CMC_chm + 1_ip) * freq_CFD_coup_CMC_chm
              cond_time_step = .false.
              if (nwrite_sol_CFD_CMC_chm /= 0_ip) &
                 cond_time_step = mod(aux_t, nwrite_sol_CFD_CMC_chm ) == 0_ip

              if (cond_time_step .or. kfl_av_species_CMC_chm == 1_ip) then

                 if (nspec_transf_CMC_chm > 0_ip) then
                    do iclas = 1, nspec_transf_CMC_chm
                       do ipoin = 1, npoin
                          do imixf = 1, nZ_CMC_chm
                             send_CMC_to_CFD_post_chm(imixf,ipoin) = &
                                Yk_CMC_chm(imixf,ipoin,transf_spec_CMC_chm(iclas))
                          end do
                       end do
                       call COU_INTERPOLATE_NODAL_VALUES(icoup,nZ_CMC_chm,dummr2,send_CMC_to_CFD_post_chm)
                    end do
                 end if
              end if

              if (cond_time_step .or. kfl_av_enthalp_CMC_chm == 1_ip) then

                 if (transf_entha_CMC_chm > 0_ip) then
                    if (kfl_solve_enth_CMC_chm /= 0_ip) then
                       do ipoin = 1, npoin
                          do imixf = 1, nZ_CMC_chm
                             send_CMC_to_CFD_post_chm(imixf,ipoin) = &
                                enthalp_CMC_chm(imixf,ipoin)
                          end do
                       end do
                    else
                       do ipoin = 1, npoin
                          do imixf = 1, nZ_CMC_chm
                             send_CMC_to_CFD_post_chm(imixf,ipoin) = &
                                  enthalp_CMC_chm(imixf,1)
                          end do
                       end do
                    end if
                    call COU_INTERPOLATE_NODAL_VALUES(icoup,nZ_CMC_chm,dummr2,send_CMC_to_CFD_post_chm)
                 end if
              end if
           end if

        else
           do ipoin = 1, npoin
              send_CMC_to_CFD_chm(1,1,ipoin) = densi_int_CMC_chm(ipoin)
              send_CMC_to_CFD_chm(2,1,ipoin) = visco_lam_int_CMC_chm(ipoin)
              send_CMC_to_CFD_chm(3,1,ipoin) = condu_int_CMC_chm(ipoin)
              send_CMC_to_CFD_chm(4,1,ipoin) = sphea_int_CMC_chm(ipoin)
              send_CMC_to_CFD_chm(5,1,ipoin) = temp_int_CMC_chm(ipoin)
           end do

           call COU_INTERPOLATE_NODAL_VALUES(icoup,5_ip,dummr2,send_CMC_to_CFD_chm)

        end if

     end if

     !
     ! Receive information by CFD from CMC: density, laminar viscosity,
     ! conductivity, cp and temperature
     !
     if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 0) then

        ! Transference of information
        if (kfl_transfer_condField_CMC_chm == 1_ip) then

           if (freq_CFD_coup_CMC_chm  > 1_ip .and. kfl_start_CMC_chm == 0_ip) then
              do ipoin = 1, npoin
                 rec_CFD_old_from_CMC_chm(:,:,ipoin) = rec_CFD_from_CMC_chm(:,:,ipoin)
              end do
           end if

           call COU_INTERPOLATE_NODAL_VALUES(icoup,6_ip*nZ_CMC_chm,rec_CFD_from_CMC_chm)

           if (kfl_start_CMC_chm == 1)  then
              if (freq_CFD_coup_CMC_chm  > 1_ip) then
                 do ipoin = 1, npoin
                    rec_CFD_old_from_CMC_chm(:,:,ipoin) = rec_CFD_from_CMC_chm(:,:,ipoin)
                 end do
              end if
              call chm_integrate_var_CFD_CMC(6_ip)
           end if

           if (nvar_trans_CMC_chm > 6_ip) then

              cond_time_step = .false.
              if (nwrite_sol_CFD_CMC_chm /= 0_ip) &
                 cond_time_step = mod(ittim+freq_CFD_coup_CMC_chm, nwrite_sol_CFD_CMC_chm ) == 0_ip

              if (nspec_transf_CMC_chm > 0_ip) then
                 if (cond_time_step .or. kfl_av_species_CMC_chm == 1_ip) then

                    nvar = 1_ip
                    allocate(aux_var_cond(nZ_CMC_chm,nvar))
                    allocate(aux_var_uncond(nvar))

                    do iclas = 1, nspec_transf_CMC_chm
                       call COU_INTERPOLATE_NODAL_VALUES(icoup,nZ_CMC_chm,rec_CFD_from_CMC_post_chm)

                       do ipoin = 1, npoin
                          aux_var_cond(1:nZ_CMC_chm,1) = rec_CFD_from_CMC_post_chm(1:nZ_CMC_chm,ipoin)

                          ! Integrations
                          call chm_mxt_fr_integr_previous_steps_CMC(1_ip, aux_var_cond, &
                                  conce(ipoin,3,1), conce(ipoin,4,1), &
                                  aux_var_uncond)

                          Yk_int_CMC_chm(ipoin,iclas) = aux_var_uncond(1)

                       end do
                    end do

                    deallocate(aux_var_cond)
                    deallocate(aux_var_uncond)
                 end if
              end if

              if (transf_entha_CMC_chm > 0_ip) then
                 if (cond_time_step .or. kfl_av_enthalp_CMC_chm == 1_ip) then

                    nvar = 1_ip
                    allocate(aux_var_cond(nZ_CMC_chm,nvar))
                    allocate(aux_var_uncond(nvar))

                    call COU_INTERPOLATE_NODAL_VALUES(icoup,nZ_CMC_chm,rec_CFD_from_CMC_post_chm)

                    do ipoin = 1, npoin
                       aux_var_cond(1:nZ_CMC_chm,1) = rec_CFD_from_CMC_post_chm(1:nZ_CMC_chm,ipoin)

                       ! Integrations
                       call chm_mxt_fr_integr_previous_steps_CMC(1_ip, aux_var_cond, &
                               conce(ipoin,3,1), conce(ipoin,4,1), &
                               aux_var_uncond)

                       enthalp_int_CMC_chm(ipoin) = aux_var_uncond(1)

                    end do

                    deallocate(aux_var_cond)
                    deallocate(aux_var_uncond)
                 end if

              end if
           end if

        else
           call COU_INTERPOLATE_NODAL_VALUES(icoup,5_ip,rec_CFD_from_CMC_chm)

           do ipoin = 1, npoin
              densi(ipoin,1) = rec_CFD_from_CMC_chm(1,1,ipoin)
              visco(ipoin,1) = rec_CFD_from_CMC_chm(2,1,ipoin)
              condk(ipoin,1) = rec_CFD_from_CMC_chm(3,1,ipoin)
              sphek(ipoin,1) = rec_CFD_from_CMC_chm(4,1,ipoin)
              tempe(ipoin,1) = rec_CFD_from_CMC_chm(5,1,ipoin)
           end do
        end if

        call ker_updpro()  ! Update properties in kernel for CFD

     end if

  end subroutine chm_data_CMC_to_CFD_domain_CMC


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! T R A N S F E R   L O C A L   <--->   G L O B A L !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine chm_global2local_CMC(imixf)
     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_global2local_CMC
     ! NAME
     !    chm_global2local_CMC
     ! DESCRIPTION
     !    This routine transfers information from Yk_CMC_chm, enthalp_CMC_chm,
     !    etc. to conce, therm, etc.
     ! USED BY
     !***
     !-----------------------------------------------------------------------

     use def_master,         only : therm, tempe, conce
     use def_chemic,         only : enthalp_CMC_chm, temp_CMC_chm, Yk_CMC_chm, &
                                    kfl_solve_enth_CMC_chm

     implicit none
     integer(ip), intent(in)     :: imixf
     integer(ip)                 :: ipoin, iclas

     if( INOTMASTER ) then
        ! Assign conce, therm and tempe
        do ipoin = 1,npoin
           tempe(ipoin,1) = temp_CMC_chm(imixf,ipoin)  ! Used to obtain properties
           do iclas = 1,nclas_chm
              conce(ipoin,iclas,1) = Yk_CMC_chm(imixf,ipoin,iclas)
           end do
        end do

        if (kfl_solve_enth_CMC_chm == 0_ip) then
           therm(1:npoin,1) = enthalp_CMC_chm(imixf,1)
        else
           do ipoin = 1,npoin
              therm(ipoin,1) = enthalp_CMC_chm(imixf,ipoin)
           end do
        end if
     end if

  end subroutine chm_global2local_CMC


  subroutine chm_local2global_CMC(imixf)
     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_local2global_CMC
     ! NAME
     !    chm_local2global_CMC
     ! DESCRIPTION
     !    This routine transfers information from conce, therm, etc. to Yk_CMC_chm,
     !    enthalp_CMC_chm, etc.
     ! USED BY
     !***
     !-----------------------------------------------------------------------

     use def_master,         only : therm, conce
     use def_chemic,         only : enthalp_CMC_chm, Yk_CMC_chm, &
                                    kfl_solve_enth_CMC_chm

     implicit none
     integer(ip), intent(in)     :: imixf
     integer(ip)                 :: ipoin, iclas

     if( INOTMASTER ) then
        ! Assign conce, therm and tempe
        do ipoin = 1,npoin
           do iclas = 1,nclas_chm
              Yk_CMC_chm(imixf,ipoin,iclas) = conce(ipoin,iclas,1)
           end do
        end do

        if (kfl_solve_enth_CMC_chm /= 0_ip) then
           do ipoin = 1,npoin
              enthalp_CMC_chm(imixf,ipoin) = therm(ipoin,1)
           end do
        end if
     end if

  end subroutine chm_local2global_CMC


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! T E R M S   M O D E L L I N G   A N D   A S S E M B L Y !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine chm_element_operations_CMC(order,pnode,pgaus,list_elements,imixf)
     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_element_operations_CMC
     ! NAME
     !    chm_element_operations_CMC
     ! DESCRIPTION
     !    This routine computes the conditioned variables not solved in CMC,
     !    that is, the terms that appear in the CMC transport equations and
     !    have to be modelled for a given mixture fraction level. Then, it
     !    assembles the terms in order to solve CMC transport
     !    equations for all the species and one mixture fraction level.
     ! USED BY
     !    chm_elmope_all
     !***
     !-----------------------------------------------------------------------
     use def_master,                        only : densi_gp, solve, cutim, rhsid
     use def_domain,                        only : ntens, elmar, ltype, lorde, ltopo, hnatu, lnods
     use def_chemic,                        only : Z_CMC_chm, Zs_CMC_chm, &
                                                   condu_gp_CMC_chm, sphec_gp_CMC_chm, &
                                                   spvol_gp_CMC_chm, &
                                                   ADR_chm, kfl_taust_chm, kfl_advec_chm, &
                                                   kfl_ellen_chm, kfl_entropy_chm, nvar_CMC_chm, &
                                                   kfl_solve_enth_CMC_chm, Le_k, &
                                                   dt_rho_chm, kfl_mesh_interp_CMC_chm, &
                                                   kfl_avg_cond_CMC_chm, kfl_incl_PDF_trans_CMC_chm, &
                                                   kfl_soot_chm
     use mod_chm_sectional_soot_model_fast, only : assembly_soot_sources_ssm
     use mod_matrix,                        only : matrix_assexp

     use mod_ADR,                           only : ADR_element_assembly
     use mod_ADR,                           only : ELEMENT_ASSEMBLY                 ! 1
     use mod_ADR,                           only : mreac_adr
     use mod_chm_entropy,                   only : chm_entropy_viscosity

     implicit none
     integer(ip), intent(in)          :: order             ! =1 defaul or =2 compute SGS only
     integer(ip), intent(in)          :: pnode             ! Number of nodes
     integer(ip), intent(in)          :: pgaus             ! Number of Gauss points
     integer(ip), intent(in), pointer :: list_elements(:)  ! List of elements
     integer(ip), intent(in)          :: imixf             ! Level of mixture fraction

     integer(ip)              :: iclas,kelem,ielem,igaus
     integer(ip)              :: pelty,porde,ptopo,plapl,izmat,izrhs

     real(rp)                 :: chale(3),chave(3),hleng(3),tragl(9)
     real(rp)                 :: dummr(ndime,mnode), kcp_la

     ! Matrices for elements (values at nodes)
     real(rp)                 :: elmat(mnode,mnode)                     ! Elemental matrix contribution
     real(rp)                 :: elrhs(mnode)                           ! Right hand side
     real(rp)                 :: elcod(ndime,mnode)                     ! Coordinates
     real(rp)                 :: elvel_CFD(ndime,mnode)                 ! Velocity from CFD
     real(rp)                 :: elZavg_CFD(mnode)                      ! Average mixture fraction from CFD
     real(rp)                 :: elZvar_CFD(mnode)                      ! Mixture fraction variance from CFD
     real(rp)                 :: elZgrad_CFD(ndime,mnode)               ! Average mixture fraction gradient from CFD
     real(rp)                 :: elturb_dif_CFD(mnode)                  ! Mass turbulent diffusion coefficient from CFD
     real(rp)                 :: eltem(mnode)                           ! Conditional temperature
     real(rp)                 :: elunkno_CMC(mnode,nvar_CMC_chm)        ! Yk and enthalpy at iZ and nodes

     ! Matrices for elements (values at Gaussian points)
     real(rp)                 :: gpvol(mgaus)                           ! |J|*w
     real(rp)                 :: gpcar(ndime,mnode,mgaus)               ! dNk/dxj
     real(rp)                 :: gphes(ntens,mnode,mgaus)               ! d2Nk/dxidxj
     real(rp)                 :: gpdif(mgaus,nvar_CMC_chm)              ! D_k
     real(rp)                 :: gprea(mgaus,mreac_adr)                 ! r
     real(rp)                 :: gpgrd(ndime,mgaus)                     ! grad(k) = grad(D_k)
     real(rp)                 :: gprhs(mgaus,nclas_chm)                 ! f (all terms)
     real(rp)                 :: gpvel_CFD(ndime,mgaus)                 ! Velocity
     real(rp)                 :: gpZavg_CFD(mgaus)                      ! Average mixture fraction
     real(rp)                 :: gpZvar_CFD(mgaus)                      ! Mixture fraction variance
     real(rp)                 :: gpZgrad_CFD(ndime,mgaus)               ! Average mixture fraction gradient
     real(rp)                 :: gpturb_dif_CFD(mgaus)                  ! Mass turbulent diffusion coefficient
     real(rp)                 :: del_gpdif(mgaus,nvar_CMC_chm)          ! change in diffusivity from entropy stable stabilization
     real(rp)                 :: PDF_val(mgaus)
     real(rp)                 :: gp_densi(mgaus)                        ! Uncond. rho
     real(rp)                 :: gp_densi_PDF_CMC(mgaus)                ! Uncond. rho * probability density function
     real(rp)                 :: gptem(mgaus)                           ! Conditional emperature
     real(rp)                 :: gp_veloc_CMC_chm(ndime,mgaus)          ! Conditional velocity
     real(rp)                 :: gp_diff_phys_spc(mgaus,nvar_CMC_chm)   ! Diffusion in physical space
     real(rp)                 :: gp_diff_mf(mgaus,nvar_CMC_chm)         ! Diffusion in mixture fraction direction
     real(rp)                 :: gpXcond_CMC                            ! Conditional scalar dissipation rate
     real(rp)                 :: gpPDF_param(7)                         ! PDF parameters at Gaussian points
     real(rp)                 :: factor0_vel(ndime)                     ! Constant in the conditional velocity model
     real(rp)                 :: factor1_vel(ndime)                     ! Slope in the conditional velocity model
     real(rp)                 :: gpX0_CMC                               ! Modulator scalar dissipation rate for AMC model
     real(rp)                 :: gpunkno_CMC(mgaus,nvar_CMC_chm)        ! Yk and enthalpy at iZ and Gaussian points

     external                 :: elmlen
     external                 :: elmchl
     external                 :: elmcar
     external                 :: chm_rhodt

     ! Loop over all the elements
     elements: do kelem = 1_ip, size(list_elements, kind=ip)
        ielem = list_elements(kelem)
        if( ielem > 0 ) then
           !
           ! Element dimensions
           !
           pelty = ltype(ielem)
           porde = lorde(pelty)
           ptopo = ltopo(pelty)

           !
           ! Initialization of variables in Gaussian points
           !
           elmat              = 0.0_rp
           elrhs              = 0.0_rp
           elcod              = 0.0_rp
           elvel_CFD          = 0.0_rp
           elZavg_CFD         = 0.0_rp
           elZvar_CFD         = 0.0_rp
           elZgrad_CFD        = 0.0_rp
           elturb_dif_CFD     = 0.0_rp
           eltem              = 0.0_rp
           elunkno_CMC        = 0.0_rp
           gpvol              = 0.0_rp
           gpcar              = 0.0_rp
           gphes              = 0.0_rp
           gpvel_CFD          = 0.0_rp
           gpZavg_CFD         = 0.0_rp
           gpZvar_CFD         = 0.0_rp
           gpZgrad_CFD        = 0.0_rp
           gpturb_dif_CFD     = 0.0_rp
           PDF_val            = 0.0_rp
           gp_densi           = 0.0_rp
           gp_densi_PDF_CMC   = 0.0_rp
           gptem              = 0.0_rp
           gp_veloc_CMC_chm   = 0.0_rp
           gpdif              = 0.0_rp
           gp_diff_mf         = 0.0_rp
           gp_diff_phys_spc   = 0.0_rp
           gpXcond_CMC        = 0.0_rp
           gpX0_CMC           = 0.0_rp
           gpPDF_param        = 0.0_rp
           gpPDF_param(4)     = Zs_CMC_chm
           gprea              = 0.0_rp
           gpgrd              = 0.0_rp
           gprhs              = 0.0_rp
           gpunkno_CMC        = 0.0_rp
           del_gpdif          = 0.0_rp
           factor0_vel        = 0.0_rp
           factor1_vel        = 0.0_rp
           dummr              = 0.0_rp
           plapl              = 0_ip


           !
           ! Gather values at the element. The turbulent diffusion coef. in elturb_dif_CFD is for mass
           !
           call chm_elmgac_CMC(&
                   pnode, lnods(1:pnode,ielem), elcod(:,1:pnode), elvel_CFD(:,1:pnode), &
                   elZavg_CFD(1:pnode), elZvar_CFD(1:pnode), elZgrad_CFD(:,1:pnode), &
                   elturb_dif_CFD(1:pnode), eltem(1:pnode), elunkno_CMC(1:pnode,:), imixf)
           ! FALTA TERMOFORESIS -> paso a matriz elemental

           !
           ! CHALE, HLENG and TRAGL
           !
           if( kfl_taust_chm /= 0 ) then
              call elmlen(&
                   ndime,pnode,elmar(pelty)%dercg,tragl,elcod(:,1:pnode),hnatu(pelty),hleng)
              call elmchl(&
                   tragl,hleng,elcod(:,1:pnode),dummr(:,1:pnode),chave,chale,pelty,pnode,porde,hnatu(pelty),&
                   kfl_advec_chm,kfl_ellen_chm)
           else
              plapl = 0
           end if

           !
           ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
           !
           call elmcar(&
               pnode, pgaus, plapl, elmar(pelty)%weigp, elmar(pelty)%shape,&
               elmar(pelty)%deriv, elmar(pelty)%heslo, elcod(:,1:pnode), gpvol(1:pgaus), &
               gpcar(:,:,1:pgaus), gphes(:,:,1:pgaus), ielem)

           !
           ! Send quantities to Gaussian points
           !
           call chm_elmpre_CMC(&
                   pnode, pgaus, elmar(pelty)%shape, elvel_CFD(:,1:pnode), elZavg_CFD(1:pnode), &
                   elZvar_CFD(1:pnode), elZgrad_CFD(:,1:pnode), &
                   elturb_dif_CFD(1:pnode), eltem(1:pnode), elunkno_CMC(1:pnode,:), &
                   gpvel_CFD(:,1:pgaus), gpZavg_CFD(1:pgaus), &
                   gpZvar_CFD(1:pgaus), gpZgrad_CFD(:,1:pgaus), &
                   gpturb_dif_CFD(1:pgaus), gptem(1:pgaus), gpunkno_CMC(1:pgaus,:))

           ! FALTA TERMOFORESIS -> enviar a puntos de Gauss


           !
           ! Compute the laminar diffusion coefficient D at Gaussian points
           !
           do igaus = 1, pgaus
              kcp_la = condu_gp_CMC_chm(ielem)%a(igaus,imixf,1) * spvol_gp_CMC_chm(ielem)%a(igaus,imixf,1) / &
                                 sphec_gp_CMC_chm(ielem) % a(igaus,imixf,1)

              do iclas = 1,nspec_chm
                 gpdif(igaus,iclas) = kcp_la / Le_k(iclas)  ! It contains rho*D and not only D
              end do
              if (kfl_solve_enth_CMC_chm /= 0_ip) then  ! In case enthalpy is transported
                 gpdif(igaus,nvar_CMC_chm) = kcp_la
              end if
           end do

           !
           ! Entropy stable viscosity. !!!!!!!!!!!!!!!!NEED TO BE INCLUDED?
           !
           !!!!!!!!!!!!!!!!!! CHECK THE FOLLOWING LINES
           if(kfl_entropy_chm == 1_ip) then
              del_gpdif = 0.0_rp
              do iclas = 1,nvar_CMC_chm
              !!!! AHORA NO SE PUEDE USAR PORQUE LOS TAMAÑOS DE LOS VECTORES NO SON COHERENTES CON LOS DE chm_entropy_viscosity
              ! call chm_entropy_viscosity(ielem,pnode,pgaus,1_ip,pgaus,iclas,elvel,gpden,hleng,del_gpdif)
              enddo

              !   do iclas = iclai_chm,iclaf_chm
              !      gpdif(1:pgaus,iclas) = gpdif(1:pgaus,iclas) + del_gpdif(1:pgaus,iclas)
              !   enddo ! At this point gpdif contains the laminar and entropy viscosity contributions
           end if

           !
           ! Find values at Gauss points
           !
           do igaus = 1, pgaus
              gp_densi(igaus)         = densi_gp(ielem) % a(igaus,1,1)
              gp_densi_PDF_CMC(igaus) = densi_gp(ielem) % a(igaus,1,1)
           end do

           if (kfl_incl_PDF_trans_CMC_chm == 1_ip) then
              do igaus = 1, pgaus
                 ! Compute PDF value
                 call chm_PDF_calc_CMC(gpZavg_CFD(igaus), gpZvar_CFD(igaus), imixf, PDF_val(igaus))

                 ! Find density times PDF
                 gp_densi_PDF_CMC(igaus) = gp_densi_PDF_CMC(igaus) * PDF_val(igaus)

              end do
           end if

           ! Diffusion in physical space
           do iclas = 1,nvar_CMC_chm
              gp_diff_phys_spc(1:pgaus,iclas) = gp_densi_PDF_CMC(1:pgaus) * &
                                               (gpdif(1:pgaus,iclas) + gpturb_dif_CFD(1:pgaus))
           end do


           ! Conditional velocity
           if (kfl_mesh_interp_CMC_chm == 1_ip .and. kfl_avg_cond_CMC_chm == 1_ip) then
              gp_veloc_CMC_chm(1:ndime,1:pgaus) = gpvel_CFD(1:ndime,1:pgaus)
           else
              do igaus = 1, pgaus
                 call chm_find_veloc_parameters_CMC(gpvel_CFD(1:ndime,igaus), &
                           factor0_vel(1:ndime), factor1_vel(1:ndime))
                 gp_veloc_CMC_chm(1:ndime,igaus) = factor0_vel(1:ndime) + factor1_vel(1:ndime) * Z_CMC_chm(imixf)
              end do
           end if

           ! Once values at all Gaussian points are computed, do assembly, if required
           !
           ! Projections of rho/dt and 1/dt
           !
           call chm_rhodt(  &
                pnode,pgaus,porde,lnods(1:pnode,ielem),elmar(pelty)%shape,gpvol(1:pgaus), &
                gp_densi_PDF_CMC(1:pgaus),dt_rho_chm)

           !
           ! Calculate RHS soot source terms
           !
           if (kfl_soot_chm > 0 ) then
              call assembly_soot_sources_ssm(ielem,pgaus,gpunkno_CMC(1:pgaus,1:nclas_chm), &
                     gptem(1:pgaus),gp_densi(1:pgaus),gprhs(1:pgaus,1:nclas_chm))
              if (kfl_incl_PDF_trans_CMC_chm == 1_ip) then
                 do igaus = 1, pgaus
                    gprhs(igaus,:) = gprhs(igaus,:) * PDF_val(igaus)
                 end do
              end if
           end if

           ! FALTA TERMOFORESIS -> COMPLETAR

           !
           ! Assemble matrix
           !
           izmat = 1
           izrhs = 1

           ASSEMBLY_ICLAS: do iclas = 1,nvar_CMC_chm
              !
              ! Assemble equation for iclas
              !
              if( order == ELEMENT_ASSEMBLY ) then

                 call ADR_element_assembly(&
                      ielem,pnode,pgaus,elcod(:,1:pnode),elmar(pelty)%shape, &
                      gpcar(:,:,1:pgaus),elmar(pelty)%deriv,gphes(:,:,1:pgaus), &
                      gpvol(1:pgaus),chale,&
                      elmar(pelty)%shape_bub,elmar(pelty)%deriv_bub,ADR_chm(iclas),cutim, &
                      gp_densi_PDF_CMC(1:pgaus),gp_veloc_CMC_chm(:,1:pgaus), &
                      gp_diff_phys_spc(1:pgaus,iclas), &
                      gpgrd(:,1:pgaus), gprea(1:pgaus,:), gprhs(1:pgaus,iclas), &
                      gpunkno_CMC(1:pgaus,iclas), &
                      elunkno_CMC(1:pnode,iclas),elmat(1:pnode,1:pnode),elrhs(1:pnode))


                 !
                 ! Call solver
                 !
                 call matrix_assexp(solve(1)%ndofn,1_ip,pnode,npoin,lnods(1:pnode,ielem), &
                           elrhs(1:pnode),elmat(1:pnode,1:pnode), &
                           elunkno_CMC(1:pnode,iclas),rhsid,iclas)

                 izrhs = izrhs + npoin                                       !solve(1)%nzrhs
                 izmat = izmat + solve(1)%nzmat/nvar_CMC_chm**2_ip
              end if

           end do ASSEMBLY_ICLAS

        end if
     end do elements

  end subroutine chm_element_operations_CMC


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! P O I N T S   T R A N S F E R !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine chm_elmgac_CMC(&
                 pnode, lnods, elcod, elvel_CFD, elZavg_CFD, elZvar_CFD,&
                 elZgrad_CFD, elturb_dif_CFD, eltem, elunkno_CMC, imixf)

     !------------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_elmgac_CMC
     ! NAME
     !    chm_elmgac_CMC
     ! DESCRIPTION
     !    Gather operations for the set of not solved variables that do not
     !    depend on the mixture fraction for CMC model.
     ! USES
     ! USED BY
     !    chm_element_operations_CMC
     !***
     !------------------------------------------------------------------------

     use def_chemic,      only : veloc_CMC_chm, Zavg_CMC_chm, Zvar_CMC_chm, &
                                 grad_Zavg_CMC_chm, turb_kin_visc_CMC_chm, Yk_CMC_chm, &
                                 enthalp_CMC_chm, nvar_CMC_chm, kfl_solve_enth_CMC_chm, &
                                 diffu_chm, temp_CMC_chm
     use mod_matrix,      only : matrix_assexp

     implicit none
     integer(ip), intent(in)  :: pnode
     integer(ip), intent(in)  :: lnods(pnode)
     integer(ip), intent(in)  :: imixf
     real(rp),    intent(out) :: elcod(ndime,pnode)
     real(rp),    intent(out) :: elvel_CFD(ndime,pnode)
     real(rp),    intent(out) :: elZavg_CFD(pnode)
     real(rp),    intent(out) :: elZvar_CFD(pnode)
     real(rp),    intent(out) :: elZgrad_CFD(ndime,pnode)
     real(rp),    intent(out) :: elturb_dif_CFD(pnode)
     real(rp),    intent(out) :: eltem(pnode)
     real(rp),    intent(out) :: elunkno_CMC(pnode,nvar_CMC_chm)

     integer(ip)              :: inode, ipoin

     !
     ! Initialization
     !
     elcod              = 0.0_rp
     elvel_CFD          = 0.0_rp
     elZavg_CFD         = 0.0_rp
     elZvar_CFD         = 0.0_rp
     elZgrad_CFD        = 0.0_rp
     elturb_dif_CFD     = 0.0_rp
     eltem              = 0.0_rp
     elunkno_CMC        = 0.0_rp

     !
     ! Values transference to elemental matrices
     !
     do inode = 1,pnode
        ipoin                                  = lnods(inode)
        elcod(1:ndime,inode)                   = coord(1:ndime,ipoin)
        elvel_CFD(1:ndime,inode)               = veloc_CMC_chm(1:ndime,ipoin)
        elZavg_CFD(inode)                      = Zavg_CMC_chm(ipoin)
        elZvar_CFD(inode)                      = Zvar_CMC_chm(ipoin)
        elZgrad_CFD(1:ndime,inode)             = grad_Zavg_CMC_chm(1:ndime,ipoin)
        elturb_dif_CFD(inode)                  = turb_kin_visc_CMC_chm(ipoin) / diffu_chm(1,1)
        eltem(inode)                           = temp_CMC_chm(imixf,ipoin)
        elunkno_CMC(inode,1:nclas_chm)         = Yk_CMC_chm(imixf,ipoin,1:nclas_chm)
     end do

     if (kfl_solve_enth_CMC_chm /= 0_ip) then   ! In case enthalpy is transported
        do inode = 1,pnode
           ipoin                                  = lnods(inode)
           elunkno_CMC(inode,nvar_CMC_chm)        = enthalp_CMC_chm(imixf,ipoin)
        end do
     end if

  end subroutine chm_elmgac_CMC


  subroutine chm_elmgac_coords(pnode,lnods,elcod)

     !------------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_elmgac_coords
     ! NAME
     !    chm_elmgac_coords
     ! DESCRIPTION
     !    Gather coordinates.
     ! USES
     ! USED BY
     !
     !***
     !------------------------------------------------------------------------

     implicit none
     integer(ip), intent(in)  :: pnode
     integer(ip), intent(in)  :: lnods(pnode)
     real(rp),    intent(out) :: elcod(ndime,pnode)

     integer(ip)              :: inode, ipoin

     !
     ! Initialization
     !
     elcod  = 0.0_rp

     !
     ! Values transference to elemental matrices
     !
     do inode = 1,pnode
        ipoin                = lnods(inode)
        elcod(1:ndime,inode) = coord(1:ndime,ipoin)
     end do

  end subroutine chm_elmgac_coords


  subroutine chm_elmpre_CMC(&
                   pnode, pgaus, gpsha, elvel_CFD, elZavg_CFD, elZvar_CFD, &
                   elZgrad_CFD, elturb_dif_CFD, eltem, elunkno_CMC, &
                   gpvel_CFD, gpZavg_CFD, gpZvar_CFD, gpZgrad_CFD, &
                   gpturb_dif_CFD, gptem, gpunkno_CMC)

     !------------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_elmpre_CMC
     ! NAME
     !    chm_elmpre_CMC
     ! DESCRIPTION
     !    Transfer elemental values to Gaussian points for the set of not
     !    solved variables that do not depend on the mixture fraction for
     !    CMC model.
     ! USES
     ! USED BY
     !    chm_element_operations_CMC
     !***
     !------------------------------------------------------------------------

     use def_chemic,     only :  nvar_CMC_chm

     implicit none
     integer(ip), intent(in)  :: pnode, pgaus
     real(rp),    intent(in)  :: gpsha(pnode,pgaus)
     real(rp),    intent(in)  :: elvel_CFD(ndime,pnode)
     real(rp),    intent(in)  :: elZavg_CFD(pnode)
     real(rp),    intent(in)  :: elZvar_CFD(pnode)
     real(rp),    intent(in)  :: elZgrad_CFD(ndime,pnode)
     real(rp),    intent(in)  :: elturb_dif_CFD(pnode)
     real(rp),    intent(in)  :: eltem(pnode)
     real(rp),    intent(in)  :: elunkno_CMC(pnode,nvar_CMC_chm)
     real(rp),    intent(out) :: gpvel_CFD(ndime,pgaus)
     real(rp),    intent(out) :: gpZavg_CFD(pgaus)
     real(rp),    intent(out) :: gpZvar_CFD(pgaus)
     real(rp),    intent(out) :: gpZgrad_CFD(ndime,pgaus)
     real(rp),    intent(out) :: gpturb_dif_CFD(pgaus)
     real(rp),    intent(out) :: gptem(pgaus)
     real(rp),    intent(out) :: gpunkno_CMC(pgaus,nvar_CMC_chm)

     integer(ip)              :: inode, igaus

     do igaus = 1,pgaus
        do inode = 1,pnode
           gpvel_CFD(1:ndime,igaus)                 = gpvel_CFD(1:ndime,igaus) &
                                                       + gpsha(inode,igaus) * elvel_CFD(1:ndime,inode)
           gpZavg_CFD(igaus)                        = gpZavg_CFD(igaus) &
                                                       + gpsha(inode,igaus) * elZavg_CFD(inode)
           gpZvar_CFD(igaus)                        = gpZvar_CFD(igaus) &
                                                       + gpsha(inode,igaus) * elZvar_CFD(inode)
           gpZgrad_CFD(1:ndime,igaus)               = gpZgrad_CFD(1:ndime,igaus) &
                                                       + gpsha(inode,igaus) * elZgrad_CFD(1:ndime,inode)
           gpturb_dif_CFD(igaus)                    = gpturb_dif_CFD(igaus) &
                                                       + gpsha(inode,igaus) * elturb_dif_CFD(inode)
           gptem(igaus)                             = gptem(igaus) &
                                                       + gpsha(inode,igaus) * eltem(inode)
           gpunkno_CMC(igaus,1:nvar_CMC_chm)        = gpunkno_CMC(igaus,1:nvar_CMC_chm) &
                                                       + gpsha(inode,igaus) * elunkno_CMC(inode,1:nvar_CMC_chm)
        end do
     end do

  end subroutine chm_elmpre_CMC



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! T E R M S   M O D E L L I N G !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine chm_find_veloc_parameters_CMC(gpvel_CFD, factor0, factor1)

     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_calc_conditional_veloc_CMC
     ! NAME
     !    chm_calc_conditional_veloc_CMC
     ! DESCRIPTION
     !    It computes the conditional velocity assuming a normal joint distribution
     !    between mixture fraction and velocity.
     ! USED BY
     !    chm_element_operations_CMC
     !***
     !-----------------------------------------------------------------------

     implicit none
     real(rp),    intent(in)  :: gpvel_CFD(ndime)
     real(rp),    intent(out) :: factor0(ndime), factor1(ndime)


     !!!! ¿DESCOMENTAR?
     !if (gpZavg_CFD <= extremes_Z_CMC_chm(1) .or. gpZavg_CFD >= extremes_Z_CMC_chm(2)) then
        factor0(1:ndime) = gpvel_CFD(1:ndime)
        factor1(1:ndime) = 0.0_rp
     !else
     !   S = gpZvar_CFD / (gpZavg_CFD*(Zs_CMC_chm-gpZavg_CFD))
     !   if (S <= S_threshold) then
     !      ! In this case we take <v|mixf> = <v> because if Zvar is very small grad(Z) will be
     !      ! very small but the ratio can be misleading
     !      factor0(1:ndime) = gpvel_CFD(1:ndime)
     !      factor1(1:ndime) = 0.0_rp
     !   else
     !      aux(1:ndime)     = gpturb_dif_CFD * gpZgrad_CFD(1:ndime) / gpZvar_CFD
     !      factor0(1:ndime) = gpvel_CFD(1:ndime) + aux(1:ndime) *  gpZavg_CFD
     !      factor1(1:ndime) = - aux(1:ndime)
     !   end if
     !end if

  end subroutine chm_find_veloc_parameters_CMC


  subroutine chm_find_scalar_dissip_rate_X0_CMC(Zavg, Zvar, Xtot, X0)

     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_calc_scalar_dissip_rate_CMC
     ! NAME
     !    chm_calc_scalar_dissip_rate_CMC
     ! DESCRIPTION
     !    It computes the scalar dissipation rate at Gaussian points.
     ! USED BY
     !    chm_element_operations_CMC
     !***
     !-----------------------------------------------------------------------

     use def_chemic,      only: Zs_CMC_chm, Z_AMC_CMC_chm, S_AMC_CMC_chm, &
                                nZ_AMC_CMC_chm, nS_AMC_CMC_chm, Smax_AMC_CMC_chm, &
                                Xintegrated_table_AMC_CMC_chm, extremes_Z_CMC_chm

     implicit none
     real(rp), parameter      :: small = 1.0e-4_rp

     real(rp),    intent(in)  :: Zavg
     real(rp),    intent(in)  :: Zvar
     real(rp),    intent(in)  :: Xtot
     real(rp),    intent(out) :: X0

     integer(ip)              :: index_mf1, index_mf2, index_S1, index_S2
     real(rp)                 :: S, denominator

     if (Zavg <= extremes_Z_CMC_chm(1) .or. Zavg >= extremes_Z_CMC_chm(2)) then
        ! Both Xtot and denominator will be very small so to avoid numerical noise 0 is assigned
        X0 = 0.0_rp
     else
        S = Zvar / (Zavg * (Zs_CMC_chm - Zavg))

        call get_index_vector(0.0_rp, Zs_CMC_chm, nZ_AMC_CMC_chm, 1.0_rp, Zavg, index_mf1)
        index_mf2 = index_mf1 + 1_ip

        if (S >= Smax_AMC_CMC_chm) then
           index_S2 = nS_AMC_CMC_chm
           index_S1 = index_S2 - 1_ip
        else
           call get_index_vector(0.0_rp, Smax_AMC_CMC_chm, nS_AMC_CMC_chm, 1.0_rp, S, index_S1)
           index_S2 = index_S1 + 1_ip
        end if

       ! Bilinear interpolation for AMC denominator
        call bilinear_intepolation((/Z_AMC_CMC_chm(index_mf1), Z_AMC_CMC_chm(index_mf2), S_AMC_CMC_chm(index_S1),&
            S_AMC_CMC_chm(index_S2)/), (/Xintegrated_table_AMC_CMC_chm(index_mf1,index_S1), Xintegrated_table_AMC_CMC_chm(index_mf1,&
            index_S2), Xintegrated_table_AMC_CMC_chm(index_mf2,index_S1), Xintegrated_table_AMC_CMC_chm(index_mf2,index_S2)/), Zavg,&
            S, denominator)

        if (denominator > small) then  ! The maximum value for the normalized X profile is 1,
                                       ! then small is a good estimator to discern non-meaningful values
           !!!!!! VERIFICAR SI HAY QUE INCLUIR EL 2
           X0 = Xtot / (2.0_rp * denominator)
        else
           X0 = 0.0_rp
        end if
     end if

  end subroutine chm_find_scalar_dissip_rate_X0_CMC


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! D I F F U S I O N   A N D   C O N V E C T I O N   I N   M I X T U R E   F R A C T I O N !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine chm_solve_mixt_frac_diff_CMC(dt)

     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_solve_mixt_frac_diff_CMC
     ! NAME
     !    chm_solve_mixt_frac_diff_CMC
     ! DESCRIPTION
     !    Solve mixture fraction diffusion by splitting this contribution from
     !    CMC transport equation. More particularly, it solves for all the
     !    physical points:
     !
     !     \partial Q_i                \partial^2 Q_i
     !    -------------- = <N | \eta> ----------------    i =1,..., nvar_therm_CMC_chm
     !      \partial t                  \partial Z^2
     !
     ! USED BY
     !
     !***
     !-----------------------------------------------------------------------

     use def_chemic,       only :  nvar_therm_CMC_chm, nZ_CMC_chm, &
                                   order_RK_CMC_chm, weights_int_RK_CMC_chm, &
                                   weights_ext_RK_CMC_chm, fact_diag_low_RK_CMC_chm, &
                                   fact_diag_RK_CMC_chm, fact_diag_up_RK_CMC_chm, &
                                   fact_BC_RK_CMC_chm, Yk_CMC_chm, enthalp_CMC_chm, &
                                   kfl_solve_enth_CMC_chm, Zavg_CMC_chm, Zvar_CMC_chm, &
                                   Xtot_CMC_chm, Xnormalized_prof_CMC_chm, &
                                   Xtot_whole_mesh_CMC_chm, kfl_mesh_interp_CMC_chm, &
                                   kfl_avg_cond_CMC_chm

     implicit none
     real(rp), parameter        :: small_X = 1.0e-4_rp

     real(rp), intent(in)       :: dt                    ! Time step
     real(rp), allocatable      :: diffu_coeff(:)        ! Diffusion coefficient (N+2 points)
     real(rp), allocatable      :: diag_low(:)           ! Lower diagonal
     real(rp), allocatable      :: diag(:)               ! Diagonal
     real(rp), allocatable      :: diag_up(:)            ! Upper diagonal
     real(rp), allocatable      :: B(:,:)                ! Independent term
     real(rp), allocatable      :: Y(:,:)                ! Unknowns

     integer(ip)                :: ipoin, imixf
     integer(ip)                :: compute
     real(rp)                   :: X0

     allocate(diffu_coeff(nZ_CMC_chm))  ! It includes contour points
     allocate(diag_low(nZ_CMC_chm-3))
     allocate(diag(nZ_CMC_chm-2))
     allocate(diag_up(nZ_CMC_chm-3))
     allocate(B(2,nvar_therm_CMC_chm))
     allocate(Y(nZ_CMC_chm,nvar_therm_CMC_chm))

     diffu_coeff = 0.0_rp
     diag_low    = 0.0_rp
     diag        = 0.0_rp
     diag_up     = 0.0_rp
     B           = 0.0_rp
     Y           = 0.0_rp

     if( INOTMASTER ) then
        do ipoin = 1, npoin
           compute = 0_ip

           if (kfl_mesh_interp_CMC_chm == 1_ip .and. kfl_avg_cond_CMC_chm == 1_ip) then
              diffu_coeff = Xtot_whole_mesh_CMC_chm(:,ipoin)
              if (maxval(diffu_coeff) > small_X)   compute = 1_ip
           else
              ! Find <N | \eta>
              call chm_find_scalar_dissip_rate_X0_CMC(Zavg_CMC_chm(ipoin), Zvar_CMC_chm(ipoin), Xtot_CMC_chm(ipoin), X0)
              if (X0 > small_X) then
                 diffu_coeff = X0 * Xnormalized_prof_CMC_chm
                 compute = 1_ip
              end if
           end if

           if (compute == 1_ip) then
              ! Fill Y
              do imixf = 1, nZ_CMC_chm
                 Y(imixf,1:nspec_chm) = Yk_CMC_chm(imixf,ipoin,1:nspec_chm)
              end do

              if (kfl_solve_enth_CMC_chm /= 0) then
                 do imixf = 1, nZ_CMC_chm
                    Y(imixf,nvar_therm_CMC_chm) = enthalp_CMC_chm(imixf,ipoin)
                 end do
              end if

              ! Compute tridiagonal matrix
              call RK_compute_diag_diff(nZ_CMC_chm-2_ip, diffu_coeff, fact_diag_low_RK_CMC_chm, &
                       fact_diag_RK_CMC_chm, fact_diag_up_RK_CMC_chm, &
                       diag_low, diag, diag_up)

              ! Compute independent terms
              call RK_compute_independent_term_diff(nvar_therm_CMC_chm, fact_BC_RK_CMC_chm, &
                      diffu_coeff(2), diffu_coeff(nZ_CMC_chm-1), Y(1,:), Y(nZ_CMC_chm,:), B)

              ! Solve PDE
              call RK_ODE_tridiagonal_system(nZ_CMC_chm-2_ip, nvar_therm_CMC_chm, order_RK_CMC_chm, dt, diag_low, diag, &
                     diag_up, B, weights_int_RK_CMC_chm, weights_ext_RK_CMC_chm, Y(2:nZ_CMC_chm-1,:))

              ! Assign result to CMC matrices
              do imixf = 2, nZ_CMC_chm-1
                 Yk_CMC_chm(imixf,ipoin,1:nspec_chm) = Y(imixf,1:nspec_chm)
              end do

              if (kfl_solve_enth_CMC_chm /= 0) then
                 do imixf = 2, nZ_CMC_chm-1
                    enthalp_CMC_chm(imixf,ipoin) = Y(imixf,nvar_therm_CMC_chm)
                 end do
              end if
           end if

        end do
     end if

  end subroutine chm_solve_mixt_frac_diff_CMC


  subroutine chm_solve_mixt_frac_convec_CMC(dt)

     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_solve_mixt_frac_convec_CMC
     ! NAME
     !    chm_solve_mixt_frac_convec_CMC
     ! DESCRIPTION
     !    Solve mixture fraction convection by splitting this contribution from
     !    CMC transport equation. More particularly, it solves for all the
     !    physical points:
     !
     !     \partial Q_i                       \partial Q_i
     !    -------------- + <veloc_mixt_frac> -------------- = 0    i =nspec_chm+1, nclas_chm
     !      \partial t                         \partial Z
     !
     ! USED BY
     !
     !***
     !-----------------------------------------------------------------------

     use def_chemic,       only :  nZ_CMC_chm, diff_Z_CMC_chm, &
                                   order_RK_CMC_chm, weights_int_RK_CMC_chm, &
                                   weights_ext_RK_CMC_chm, Yk_CMC_chm, &
                                   Zavg_CMC_chm, Zvar_CMC_chm, &
                                   Xtot_CMC_chm, Xnormalized_prof_CMC_chm, &
                                   Xtot_whole_mesh_CMC_chm, kfl_mesh_interp_CMC_chm, &
                                   kfl_avg_cond_CMC_chm, kfl_incl_PDF_trans_CMC_chm

     implicit none
     real(rp), parameter        :: small_X = 1.0e-4_rp

     real(rp), intent(in)       :: dt                    ! Time step
     real(rp), allocatable      :: diffu_coeff(:)        ! Diffusion coefficient (N+2 points)
     real(rp), allocatable      :: diag_low(:)           ! Lower diagonal
     real(rp), allocatable      :: diag(:)               ! Diagonal
     real(rp), allocatable      :: diag_up(:)            ! Upper diagonal
     real(rp), allocatable      :: B(:,:)                ! Independent term
     real(rp), allocatable      :: Y(:,:)                ! Unknowns
     real(rp), allocatable      :: PDF_val(:)            ! PDF
     real(rp), allocatable      :: veloc_mixt_frac(:)    ! Velocity in mixture fraction direction

     integer(ip)                :: ipoin, imixf
     integer(ip)                :: compute
     real(rp)                   :: X0

     allocate(diffu_coeff(nZ_CMC_chm))  ! It includes contour points
     allocate(diag_low(nZ_CMC_chm-3))
     allocate(diag(nZ_CMC_chm-2))
     allocate(diag_up(nZ_CMC_chm-3))
     allocate(B(2,nclas_chm-nspec_chm))
     allocate(Y(nZ_CMC_chm,nclas_chm-nspec_chm))
     allocate(PDF_val(nZ_CMC_chm-2))
     allocate(veloc_mixt_frac(nZ_CMC_chm-2))

     diffu_coeff     = 0.0_rp
     diag_low        = 0.0_rp
     diag            = 0.0_rp
     diag_up         = 0.0_rp
     B               = 0.0_rp
     Y               = 0.0_rp
     PDF_val         = 0.0_rp
     veloc_mixt_frac = 0.0_rp

     if( INOTMASTER ) then
        do ipoin = 1, npoin
           compute = 0_ip

           if (kfl_mesh_interp_CMC_chm == 1_ip .and. kfl_avg_cond_CMC_chm == 1_ip) then
              diffu_coeff = Xtot_whole_mesh_CMC_chm(:,ipoin)
              if (maxval(diffu_coeff) > small_X)   compute = 1_ip
           else
              ! Find <N | \eta>
              call chm_find_scalar_dissip_rate_X0_CMC(Zavg_CMC_chm(ipoin), Zvar_CMC_chm(ipoin), Xtot_CMC_chm(ipoin), X0)
              if (X0 > small_X) then
                 diffu_coeff = X0 * Xnormalized_prof_CMC_chm
                 compute = 1_ip
              end if
           end if

           if (compute == 1_ip) then

              if (kfl_incl_PDF_trans_CMC_chm == 0_ip) then
                 call compute_1st_deriv_2nd_order(nZ_CMC_chm, diff_Z_CMC_chm, diffu_coeff, veloc_mixt_frac)
              else
                 ! The beta pdf takes 0 or infinity at extremes -> to compute the derivative
                 ! impose 0 at extremes regardless the case
                 do imixf = 2, nZ_CMC_chm-1
                    call chm_PDF_calc_CMC(Zavg_CMC_chm(ipoin), Zvar_CMC_chm(ipoin), imixf, PDF_val(imixf-1))
                    diffu_coeff(imixf) = diffu_coeff(imixf) * PDF_val(imixf-1)
                 end do
                 call compute_1st_deriv_2nd_order(nZ_CMC_chm, diff_Z_CMC_chm, diffu_coeff, veloc_mixt_frac)
                 do imixf = 1, nZ_CMC_chm-2
                    veloc_mixt_frac(imixf) = veloc_mixt_frac(imixf) / PDF_val(imixf)
                 end do
              end if

              ! Fill Y
              do imixf = 1, nZ_CMC_chm
                 Y(imixf,:) = Yk_CMC_chm(imixf,ipoin,nspec_chm+1:nclas_chm)
              end do

              ! Compute tridiagonal matrix
              call RK_compute_diag_convec(nZ_CMC_chm-2_ip, veloc_mixt_frac, diff_Z_CMC_chm, &
                       diag_low, diag, diag_up)

              ! Compute independent terms
              call RK_compute_independent_term_convec(nclas_chm-nspec_chm, diff_Z_CMC_chm(1), &
                      diff_Z_CMC_chm(nZ_CMC_chm-1), veloc_mixt_frac(1), &
                      veloc_mixt_frac(nZ_CMC_chm-2), Y(1,:), Y(nZ_CMC_chm,:), B)

              ! Till the moment the equation has been posed as if it was dY/dt+A*Y+B=0 but
              ! RK_ODE_tridiagonal_system solves dY/dt=A*Y+B
              diag_low = - diag_low
              diag     = - diag
              diag_up  = - diag_up
              B        = - B

              ! Solve PDE
              call RK_ODE_tridiagonal_system(nZ_CMC_chm-2_ip, nclas_chm-nspec_chm, order_RK_CMC_chm, dt, diag_low, diag, &
                     diag_up, B, weights_int_RK_CMC_chm, weights_ext_RK_CMC_chm, Y(2:nZ_CMC_chm-1,:))

              ! Assign result to CMC matrices
              do imixf = 2, nZ_CMC_chm-1
                 Yk_CMC_chm(imixf,ipoin,nspec_chm+1:nclas_chm) = Y(imixf,:)
              end do

           end if

        end do
     end if

  end subroutine chm_solve_mixt_frac_convec_CMC


  subroutine RK_ODE_tridiagonal_system(N, nvar, order, dt, diag_low, diag, &
                 diag_up, B, weights_int, weights_ext, Y)

     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/RK_ODE_tridiagonal_system
     ! NAME
     !    RK_ODE_tridiagonal_system
     ! DESCRIPTION
     !    Solve heat equation PDE with a Runge-Kutta for order
     !    from 1 to 4. Thermal diffusivity is assumed constant
     !    during the time step (not updated during the sub-steps)
     !    and a second (if equispaced) or first (if not equispaced)
     !    order spatial derivative is applied to discretize the second
     !    derivative in space.
     !    In this way the equation to be solved is reduced to
     !    the ODE:
     !    d{y}/dt = [A]*{y} + {b}
     !    where [A] is a tridiagonal matrix and {b} is a null vector
     !    except at the first and last entries (BCs).
     !    The program is specifically coded to directly solve
     !    several variables (nvar) which share the same matrix [A]
     !    (vector {b} may be different). The set of all vectors
     !    {b} conform matrix [B] with shape (2,nvar). In this way
     !    we define [Y] = [{y1} {y2} ... {y_nvar}] and it is solved
     !    d[Y]/dt = [A]*[Y] + [B]
     ! USED BY
     !    chm_solve_mixt_frac_diff_CMC
     !***
     !-----------------------------------------------------------

     implicit none
     integer(ip), intent(in)     :: N                               ! Number of points in physical space (whitout considering
                                                                    !  boundary conditions)
     integer(ip), intent(in)     :: nvar                            ! Number of unknown fields
     integer(ip), intent(in)     :: order
     real(rp),    intent(in)     :: dt                              ! Time step
     real(rp),    intent(in)     :: diag_low(N-1)                   ! Lower diagonal
     real(rp),    intent(in)     :: diag(N)                         ! Diagonal
     real(rp),    intent(in)     :: diag_up(N-1)                    ! Upper diagonal
     real(rp),    intent(in)     :: B(2,nvar)                       ! Independent term matrix (2,nvar)
     real(rp),    intent(in)     :: weights_int(order-1,order-1)    ! RK coefficients to determine where evaluate A*y+b
     real(rp),    intent(in)     :: weights_ext(order)              ! RK coefficients to wheigh K
     real(rp),    intent(inout)  :: Y(N,nvar)                       ! Unknown matrix (N,nvar)

     integer(ip)                 :: nsteps, isteps
     integer(ip)                 :: ii, jj, nn
     real(rp)                    :: dt_sub                          ! Time substep
     real(rp)                    :: K(N,nvar,order)                 ! RK matrix substeps (N,nvar,order)
     real(rp)                    :: Z(N,nvar)                       ! Intermediate arguments where evaluate RK (N,nvar)

     K = 0.0_rp
     Z = 0.0_rp

     call RK_compute_dt_sub(N, diag, dt_sub)

     nsteps = floor(dt/dt_sub, kind=ip) + 1_ip
     dt_sub = dt / real(nsteps,rp)

     do isteps = 1, nsteps  ! Temporal substeps to assure stabilitiy
        ! Compute K1
        K(1,:,1) = diag(1) * Y(1,:) + diag_up(1) * Y(2,:)
        do nn = 2, N-1
           K(nn,:,1) = diag_low(nn-1) * Y(nn-1,:) + diag(nn) * Y(nn,:) + diag_up(nn) * Y(nn+1,:)
        end do
        K(N,:,1) = diag_low(N-1) * Y(N-1,:) + diag(N) * Y(N,:)
        K(1,:,1) = K(1,:,1) + B(1,:)  ! Independent term
        K(N,:,1) = K(N,:,1) + B(2,:)
        K(:,:,1) = K(:,:,1) * dt_sub

        ! Compute Ki i>1
        do ii = 2, order
           ! Compute internal argument
           Z = Y
           do jj = 1, ii-1
              if (weights_int(ii-1,jj) /= 0.0_rp) then
                 Z = Z + weights_int(ii-1,jj) * K(:,:,jj)
              end if
           end do
           K(1,:,ii) = diag(1) * Z(1,:) + diag_up(1) * Z(2,:)
           do nn = 2, N-1
              K(nn,:,ii) = diag_low(nn-1) * Z(nn-1,:) + diag(nn) * Z(nn,:) + diag_up(nn) * Z(nn+1,:)
           end do
           K(N,:,ii) = diag_low(N-1) * Z(N-1,:) + diag(N) * Z(N,:)
           K(1,:,ii) = K(1,:,ii) + B(1,:)  ! Independent term
           K(N,:,ii) = K(N,:,ii) + B(2,:)
           K(:,:,ii) = K(:,:,ii) * dt_sub
        end do

        ! Sum all the contributions
        do ii = 1, order
           Y = Y + weights_ext(ii) * K(:,:,ii)
        end do
     end do

  end subroutine RK_ODE_tridiagonal_system


  subroutine RK_compute_dt_sub(N, D, dt)

     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/RK_compute_dt_sub
     ! NAME
     !    RK_compute_dt_sub
     ! DESCRIPTION
     !    Compute delta_t that assures stability in heat equation when applying
     !    explicit RK method.
     ! USED BY
     !    RK_ODE_tridiagonal_system
     !***
     !-----------------------------------------------------------------------

     implicit none
     integer(ip), intent(in)    :: N
     real(rp), intent(in)       :: D(N)
     real(rp), intent(out)      :: dt
     real(rp)                   :: D_max

     D_max = maxval(abs(D))
     dt = 1.0_rp / D_max

  end subroutine RK_compute_dt_sub


  subroutine RK_compute_weights(order, weights_int, weights_ext)

     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/RK_compute_weights
     ! NAME
     !    RK_compute_weights
     ! DESCRIPTION
     !    Save coefficients for RK.
     ! USED BY
     !    chm_preprocessing_RK_CMC
     !***
     !-----------------------------------------------------------------------

    implicit none
    integer(ip), intent(in)    :: order
    real(rp),    intent(out)   :: weights_int(order-1, order-1)
    real(rp),    intent(out)   :: weights_ext(order)

    select case (order)
       case(1_ip)   ! Euler
          weights_ext(1) = 1.0_rp
       case(2_ip)   ! Heun
          weights_int(1,1) = 1.0_rp
          weights_ext(1) = 0.5_rp
          weights_ext(2) = 0.5_rp
       case(3_ip)   ! RK3
          weights_int    = 0.0_rp
          weights_int(1,1) = 1.0_rp
          weights_int(2,1) = -1.0_rp
          weights_int(2,2) = 2.0_rp
          weights_ext(1) = 1.0_rp/6.0_rp
          weights_ext(2) = 2.0_rp/3.0_rp
          weights_ext(3) = 1.0_rp/6.0_rp
       case(4_ip)   ! RK4
          weights_int    = 0.0_rp
          weights_int(1,1) = 0.5_rp
          weights_int(2,2) = 0.5_rp
          weights_int(3,3) = 1.0_rp
          weights_ext(1) = 1.0_rp/6.0_rp
          weights_ext(2) = 1.0_rp/3.0_rp
          weights_ext(3) = 1.0_rp/3.0_rp
          weights_ext(4) = 1.0_rp/6.0_rp
       end select
  end subroutine RK_compute_weights


  subroutine RK_compute_mesh_factor_diag(N, x_dif, fact_diag_low, fact_diag, fact_diag_up)

     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/RK_compute_mesh_factor_diag
     ! NAME
     !    RK_compute_mesh_factor_diag
     ! DESCRIPTION
     !    Compute factors related to the mesh for tridiagonal matrix [A].
     ! USED BY
     !    chm_preprocessing_RK_CMC
     !***
     !-----------------------------------------------------------------------

     implicit none
     integer(ip), intent(in)    :: N
     real(rp),    intent(in)    :: x_dif(N+1)
     real(rp),    intent(out)   :: fact_diag_low(N-1)
     real(rp),    intent(out)   :: fact_diag(N)
     real(rp),    intent(out)   :: fact_diag_up(N-1)

     integer(ip)                :: ii

     ! Lower diagonal
     do ii = 1,N-1
        fact_diag_low(ii) = 2.0_rp / (x_dif(ii+1) * (x_dif(ii+1) + x_dif(ii+2)))
     end do

     ! Diagonal
     do ii = 1,N
        fact_diag(ii) = -2.0_rp / (x_dif(ii) * x_dif(ii+1))
     end do

     ! Upper diagonal
     do ii = 1,N-1
        fact_diag_up(ii) = 2.0_rp / (x_dif(ii+1) * (x_dif(ii) + x_dif(ii+1)))
     end do

  end subroutine RK_compute_mesh_factor_diag


  subroutine RK_compute_diag_diff(N, diffu_coeff, fact_diag_low, fact_diag, fact_diag_up, &
                     diag_low, diag, diag_up)

     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/RK_compute_diag_diff
     ! NAME
     !    RK_compute_diag_diff
     ! DESCRIPTION
     !    Compute non-null values for tridiagonal matrix [A].
     ! USED BY
     !    chm_solve_mixt_frac_diff_CMC
     !***
     !-----------------------------------------------------------------------

     implicit none
     integer(ip), intent(in)    :: N
     real(rp),    intent(in)    :: diffu_coeff(N+2)
     real(rp),    intent(in)    :: fact_diag_low(N-1)
     real(rp),    intent(in)    :: fact_diag(N)
     real(rp),    intent(in)    :: fact_diag_up(N-1)
     real(rp),    intent(out)   :: diag_low(N-1)
     real(rp),    intent(out)   :: diag(N)
     real(rp),    intent(out)   :: diag_up(N-1)

     integer(ip)                :: ii

     ! Lower diagonal
     do ii = 1,N-1
        diag_low(ii) = fact_diag_low(ii) * diffu_coeff(ii+2)
     end do

     ! Diagonal
     do ii = 1,N
        diag(ii) = fact_diag(ii) * diffu_coeff(ii+1)
     end do

     ! Upper diagonal
     do ii = 1,N-1
        diag_up(ii) = fact_diag_up(ii) * diffu_coeff(ii+1)
     end do

  end subroutine RK_compute_diag_diff


  subroutine RK_compute_diag_convec(N, veloc, incr_x, &
                     diag_low, diag, diag_up)

     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/RK_compute_diag_convec
     ! NAME
     !    RK_compute_diag_convec
     ! DESCRIPTION
     !    Compute non-null values for tridiagonal matrix [A] for convective
     !    term assuming an up-wind scheme.
     ! USED BY
     !    chm_solve_mixt_frac_convec_CMC
     !***
     !-----------------------------------------------------------------------

     implicit none
     integer(ip), intent(in)    :: N
     real(rp),    intent(in)    :: veloc(N)
     real(rp),    intent(in)    :: incr_x(N+1)
     real(rp),    intent(out)   :: diag_low(N-1)
     real(rp),    intent(out)   :: diag(N)
     real(rp),    intent(out)   :: diag_up(N-1)

     integer(ip)                :: ii
     real(rp)                   :: aux

     diag_low = 0.0_rp
     diag     = 0.0_rp
     diag_up  = 0.0_rp

     if (veloc(1) > 0.0_rp) then
        diag(1) = veloc(1) / incr_x(1)
     else
        aux = veloc(1) / incr_x(2)
        diag(1) = - aux
        diag_up(1) = aux
     end if

     do ii = 2, N-1
        if (veloc(ii) > 0.0_rp) then
           aux = veloc(ii) / incr_x(ii)
           diag_low(ii-1) = - aux
           diag(ii) = aux
        else
           aux = veloc(ii) / incr_x(ii+1)
           diag(ii) = - aux
           diag_up(ii) = aux
        end if
     end do

     if (veloc(N) > 0.0_rp) then
        aux = veloc(N) / incr_x(N)
        diag_low(N-1) = - aux
        diag(N) = aux
     else
        diag(N) = - veloc(N) / incr_x(N+1)
     end if

  end subroutine RK_compute_diag_convec


  subroutine RK_compute_mesh_factor_BC(N, x_dif, fact_BC)

     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/RK_compute_mesh_factor_BC
     ! NAME
     !    RK_compute_mesh_factor_BC
     ! DESCRIPTION
     !    Compute factors related to the mesh for the BCs.
     ! USED BY
     !    chm_preprocessing_RK_CMC
     !***
     !-----------------------------------------------------------------------

     implicit none
     integer(ip), intent(in)    :: N
     real(rp),    intent(in)    :: x_dif(N+1)
     real(rp),    intent(out)   :: fact_BC(2)

     fact_BC(1) = 2.0_rp / (x_dif(1) * (x_dif(1) + x_dif(2)))
     fact_BC(2) = 2.0_rp / (x_dif(N+1) * (x_dif(N+1) + x_dif(N)))

  end subroutine RK_compute_mesh_factor_BC


  subroutine RK_compute_independent_term_diff(nvar, fact_BC, diffu_coeff_l, &
               diffu_coeff_r, Y_l, Y_r, B)

     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/RK_compute_independent_term_diff
     ! NAME
     !    RK_compute_independent_term_diff
     ! DESCRIPTION
     !    Compute independent term in heat equation discretization.
     ! USED BY
     !    chm_solve_mixt_frac_diff_CMC
     !***
     !-----------------------------------------------------------------------

     implicit none
     integer(ip), intent(in)    :: nvar
     real(rp),    intent(in)    :: fact_BC(2)
     real(rp),    intent(in)    :: diffu_coeff_l
     real(rp),    intent(in)    :: diffu_coeff_r
     real(rp),    intent(in)    :: Y_l(nvar)
     real(rp),    intent(in)    :: Y_r(nvar)
     real(rp),    intent(out)   :: B(2,nvar)

     real(rp)                   :: aux

     ! Independent term from left BC
     aux = fact_BC(1) * diffu_coeff_l
     B(1,:) = aux * Y_l

     ! Independent term from right BC
     aux = fact_BC(2) * diffu_coeff_r
     B(2,:) = aux * Y_r

  end subroutine RK_compute_independent_term_diff


  subroutine RK_compute_independent_term_convec(nvar, incr_l, incr_r, vel_l, &
               vel_r, Y_l, Y_r, B)

     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/RK_compute_independent_term_convec
     ! NAME
     !    RK_compute_independent_term_convec
     ! DESCRIPTION
     !    Compute independent term in convective equation discretization.
     ! USED BY
     !    chm_solve_mixt_frac_convec_CMC
     !***
     !-----------------------------------------------------------------------

     implicit none
     integer(ip), intent(in)    :: nvar
     real(rp),    intent(in)    :: incr_l
     real(rp),    intent(in)    :: incr_r
     real(rp),    intent(in)    :: vel_l
     real(rp),    intent(in)    :: vel_r
     real(rp),    intent(in)    :: Y_l(nvar)
     real(rp),    intent(in)    :: Y_r(nvar)
     real(rp),    intent(out)   :: B(2,nvar)


     B = 0.0_rp
     if (vel_l > 0.0_rp)  B(1,:) = - vel_l * Y_l / incr_l
     if (vel_r < 0.0_rp)  B(2,:) =   vel_r * Y_r / incr_r

  end subroutine RK_compute_independent_term_convec


  subroutine chm_preprocessing_RK_CMC

     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_preprocessing_RK_CMC
     ! NAME
     !    chm_preprocessing_RK_CMC
     ! DESCRIPTION
     !    Compute and save information used in RK for mixture fraction diffusion.
     ! USED BY
     !    chm_update_model
     !***
     !-----------------------------------------------------------------------

     use def_chemic,       only :  order_RK_CMC_chm, weights_int_RK_CMC_chm, &
                                   weights_ext_RK_CMC_chm, fact_diag_low_RK_CMC_chm, &
                                   fact_diag_RK_CMC_chm, fact_diag_up_RK_CMC_chm, &
                                   fact_BC_RK_CMC_chm, nZ_CMC_chm, diff_Z_CMC_chm

     implicit none

     if (INOTMASTER) then
        ! Compute weights for RK
        call RK_compute_weights(order_RK_CMC_chm, weights_int_RK_CMC_chm, &
                  weights_ext_RK_CMC_chm)

        ! Compute diagonal, upper and lower elements to the diagonal mesh factors
        call RK_compute_mesh_factor_diag(nZ_CMC_chm-2_ip, diff_Z_CMC_chm, &
                  fact_diag_low_RK_CMC_chm, fact_diag_RK_CMC_chm, fact_diag_up_RK_CMC_chm)

        ! Compute independent term mesh factors
        call RK_compute_mesh_factor_BC(nZ_CMC_chm-2_ip, diff_Z_CMC_chm, fact_BC_RK_CMC_chm)
     end if

  end subroutine chm_preprocessing_RK_CMC



  !!!!!!!!!!!!!!!!!!!!!!!!!
  !!! A M C   M O D E L !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine chm_AMC_generate_Z_S_vectors_CMC

     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_AMC_generate_Z_S_vectors_CMC
     ! NAME
     !    chm_AMC_generate_Z_S_vectors_CMC
     ! DESCRIPTION
     !    It generates the vectors for mixture fraction and segregation factor
     !    that define the table where the integrals for AMC model are saved.
     ! USED BY
     !
     !***
     !-----------------------------------------------------------------------

     use def_chemic,     only:  nZ_AMC_CMC_chm, nS_AMC_CMC_chm, Smax_AMC_CMC_chm, &
                                Zs_CMC_chm, Z_AMC_CMC_chm, S_AMC_CMC_chm

     implicit none
     real(rp)                :: delta_Z, delta_S
     integer(ip)             :: iZ, iS

     delta_Z = Zs_CMC_chm / real(nZ_AMC_CMC_chm-1,rp)
     Z_AMC_CMC_chm(nZ_AMC_CMC_chm) = Zs_CMC_chm
     do iZ = 2,nZ_AMC_CMC_chm-1_ip
        Z_AMC_CMC_chm(iZ) = real((iZ-1_ip),rp) * delta_Z
     end do

     delta_S = Smax_AMC_CMC_chm / real((nS_AMC_CMC_chm-1),rp)
     S_AMC_CMC_chm(nS_AMC_CMC_chm) = Smax_AMC_CMC_chm
     do iS = 2,nS_AMC_CMC_chm-1_ip
        S_AMC_CMC_chm(iS) = real((iS-1_ip),rp) * delta_S
     end do

  end subroutine chm_AMC_generate_Z_S_vectors_CMC


  subroutine chm_AMC_integrals_CMC

     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_AMC_integrals_CMC
     ! NAME
     !    chm_AMC_integrals_CMC
     ! DESCRIPTION
     !    It generates a matrix that contains the integrals for AMC model.
     ! USED BY
     !
     !***
     !-----------------------------------------------------------------------

     use def_chemic,     only:  nZ_AMC_CMC_chm, nS_AMC_CMC_chm, Zs_CMC_chm, &
                                Xintegrated_table_AMC_CMC_chm, Xnormalized_prof_CMC_chm, &
                                Z_AMC_CMC_chm, S_AMC_CMC_chm, S_threshold

     implicit none
     integer(ip)             :: iZ, iS
     real(rp)                :: Zvar, PDF_param(7), aux(1)
     character(len=20)       :: names_rscal(3)

     PDF_param(4) = Zs_CMC_chm

     do iZ = 2, nZ_AMC_CMC_chm-1_ip
        do iS = 1, nS_AMC_CMC_chm
           if (S_AMC_CMC_chm(iS) <= S_threshold) then
              ! Direct interpolation
              call chm_mixture_fraction_interpolation_CMC(Xnormalized_prof_CMC_chm, 1_ip, &
                               Z_AMC_CMC_chm(iZ), aux(1))
           else
              ! Perform integration
              Zvar = S_AMC_CMC_chm(iS) * Z_AMC_CMC_chm(iZ) * (Zs_CMC_chm - Z_AMC_CMC_chm(iZ))
              call chm_find_PDF_parameters_CMC(Z_AMC_CMC_chm(iZ), Zvar, PDF_param)
              call chm_mixtFraction_integration_CMC(Xnormalized_prof_CMC_chm, 1_ip, &
                               PDF_param, aux(1))
           end if
           Xintegrated_table_AMC_CMC_chm(iZ,iS) = aux(1)
        end do
     end do

     ! Write file with integrated normalized profiles

     ! Get the names in a vector for the header
     names_rscal(1) = 'MF_AVG [-]:1'
     names_rscal(2) = 'S [-]:2'
     names_rscal(3) = 'X_INTEGRAL [-]:3'

     ! Write inert mixture file
     open(unit=1, file='Xintegrals.log', status="replace", action="write")
     write(unit=1,fmt="(300(a20))") names_rscal
     do iZ = 1, nZ_AMC_CMC_chm
        do iS = 1, nS_AMC_CMC_chm
           write(unit=1,fmt="(3(e13.6,8x))") Z_AMC_CMC_chm(iZ), S_AMC_CMC_chm(iS), Xintegrated_table_AMC_CMC_chm(iZ,iS)
        end do
     end do
     close(unit=1)

  end subroutine chm_AMC_integrals_CMC


  subroutine get_index_vector(bound_min, bound_max, n_vec, expon, val, index_v)

     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/get_index_vector
     ! NAME
     !    get_index_vector
     ! DESCRIPTION
     !    Given a value it finds the rounded down index for which the corresponding
     !    value in a vector defined by a potential distribution is closest.
     ! USED BY
     !    chm_calc_scalar_dissip_rate_CMC
     !***
     !-----------------------------------------------------------------------

     implicit none
     integer(ip), intent(in)          :: n_vec
     integer(ip), intent(out)         :: index_v
     real(rp), intent(in)             :: bound_min, bound_max, expon, val

     if (expon == 1.0_rp) then  ! Separate from others just to avoid numerical noise when doing 1.0_rp/1.0_rp
        index_v = floor(1.0_rp + real(n_vec-1_ip,rp) * ((val - bound_min) / (bound_max - bound_min)), kind=ip)
     else
        index_v = floor(1.0_rp + real(n_vec-1,rp) * ((val - bound_min) / (bound_max - bound_min))**(1.0_rp/expon), kind=ip)
     end if

  end subroutine get_index_vector


  subroutine bilinear_intepolation(coordinates, heights, x, y, z)

     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/bilinear_intepolation
     ! NAME
     !    bilinear_intepolation
     ! DESCRIPTION
     !    Do a bilinear interpolation.
     ! USED BY
     !    chm_calc_scalar_dissip_rate_CMC
     !***
     !-----------------------------------------------------------------------

     implicit none
     real(rp), intent(in)       :: coordinates(4), heights(4), x, y
     real(rp), intent(out)      :: z
     real(rp)                   :: z1, z2

     z1 = heights(1) + (heights(3) - heights(1)) * (x-coordinates(1)) / (coordinates(2) - coordinates(1))
     z2 = heights(2) + (heights(4) - heights(2)) * (x-coordinates(1)) / (coordinates(2) - coordinates(1))
     z  = z1 + (z2 -z1) * (y-coordinates(3)) / (coordinates(4) - coordinates(3))

  end subroutine bilinear_intepolation


  subroutine compute_Xnormalized_profile_CMC

     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/compute_Xnormalized_profile_CMC
     ! NAME
     !    compute_Xnormalized_profile_CMC
     ! DESCRIPTION
     !    Compute normalized scalar dissipation rate profile for AMC model from
     !    CMC model.
     ! USED BY
     !    chm_reaphy
     !***
     !-----------------------------------------------------------------------

     use def_chemic,         only:  Xnormalized_prof_CMC_chm, Z_CMC_chm, nZ_CMC_chm, &
                                    Zs_CMC_chm

     implicit none
     real(rp), parameter         :: small = 1.0e-8_rp

     integer(ip)                 :: iZ
     real(rp)                    :: aux(nZ_CMC_chm), aux2(nZ_CMC_chm)

     aux = 2.0_rp * Z_CMC_chm/Zs_CMC_chm - 1.0_rp
     aux2(1:nZ_CMC_chm) = 0.0_rp

     do iZ = 2,nZ_CMC_chm-1
        call compute_erfinv(aux(iZ),aux2(iZ))
        Xnormalized_prof_CMC_chm(iZ) = exp(-2.0_rp * aux2(iZ)**2.0_rp) * Zs_CMC_chm**2.0_rp
        if (Xnormalized_prof_CMC_chm(iZ) < small) then
           Xnormalized_prof_CMC_chm(iZ) = 0.0_rp
        end if
     end do

  end subroutine compute_Xnormalized_profile_CMC


  subroutine chm_compute_interp_fields_mesh_CMC(imixf)
     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_compute_interp_fields_mesh_CMC
     ! NAME
     !    chm_compute_interp_fields_mesh_CMC
     ! DESCRIPTION
     !    Compute interpolated fields when using different meshes in CMC and
     !    CFD and the volumetric integrals are weighed by the PDF.
     ! USED BY
     !    mod_chm_rk_explicit
     !***
     !-----------------------------------------------------------------------

     use def_chemic,          only :  aux_cond_fields_CMC_chm, &
                                      turb_kin_visc_CMC_chm, &
                                      grad_Zavg_CMC_chm, veloc_CMC_chm, &
                                      aux_interp_fields_CMC_chm, rec_CMC_from_CFD_chm, &
                                      Z_CMC_chm, Xnormalized_prof_CMC_chm, &
                                      Xtot_whole_mesh_CMC_chm
     use def_domain,          only :  gp_total_cloud_gp

     implicit none
     integer(ip), intent(in)       :: imixf
     integer(ip)                   :: irow, ipoin, igaus
     real(rp)                      :: X0, factor0_vel(ndime), factor1_vel(ndime)

     do igaus = 1, gp_total_cloud_gp
        ! Find <N | \eta>
        call chm_find_scalar_dissip_rate_X0_CMC(rec_CMC_from_CFD_chm(1,igaus), &
                rec_CMC_from_CFD_chm(2,igaus), rec_CMC_from_CFD_chm(3,igaus), X0)
        aux_cond_fields_CMC_chm(1,igaus) = X0 * Xnormalized_prof_CMC_chm(imixf)

        ! Turbulent viscosity
        aux_cond_fields_CMC_chm(2,igaus) = rec_CMC_from_CFD_chm(4,igaus)

        ! Mixture fraction gradient
        aux_cond_fields_CMC_chm(3:2+ndime,igaus) = rec_CMC_from_CFD_chm(5:4+ndime,igaus)

        ! Velocity
        call chm_find_veloc_parameters_CMC(rec_CMC_from_CFD_chm(5+ndime:4+2*ndime,igaus), &
                factor0_vel(1:ndime), factor1_vel(1:ndime))
        aux_cond_fields_CMC_chm(3+ndime:2+2*ndime,igaus) = &
            factor0_vel(1:ndime) + factor1_vel(1:ndime) * Z_CMC_chm(imixf)
     end do

     ! In a similar way rec_CMC_from_CFD_chm has to be passed as
     ! (/(rec_CMC_from_CFD_chm(irow,:), irow = 1_ip, 4_ip+2_ip*ndime)/)
     ! Just giving rec_CMC_from_CFD_chm the subroutine will not work
     ! as desired.

     call chm_average_field_mesh_CMC(2_ip+2_ip*ndime, &
             (/(aux_cond_fields_CMC_chm(irow,:), irow = 1_ip, 2_ip+2_ip*ndime)/), &
             aux_interp_fields_CMC_chm, rec_CMC_from_CFD_chm(1,:), rec_CMC_from_CFD_chm(2,:), imixf)

     do ipoin = 1, npoin
        ! Diff. in mixt. frac. solved after RK -> save Xtot in a matrix for the whole mesh
        Xtot_whole_mesh_CMC_chm(imixf,ipoin) = aux_interp_fields_CMC_chm(ipoin,1)
        turb_kin_visc_CMC_chm(ipoin)         = aux_interp_fields_CMC_chm(ipoin,2)
        ! ¿ELIMINAR? NO SE ESTÁ USANDO EL GRAD_Z CON INTERP. COND.
        grad_Zavg_CMC_chm(1:ndime,ipoin)     = aux_interp_fields_CMC_chm(ipoin,3:2+ndime)
        veloc_CMC_chm(1:ndime,ipoin)         = aux_interp_fields_CMC_chm(ipoin,3+ndime:2+2*ndime)
     end do

  end subroutine chm_compute_interp_fields_mesh_CMC


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! M I X T.   F R A C T.  D I F F U S I O N !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine chm_calc_diff_condVar_mixfraction_CMC

     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_calc_diff_condVar_mixfraction_CMC
     ! NAME
     !    chm_calc_diff_condVar_mixfraction_CMC
     ! DESCRIPTION
     !    It computes the second derivative of the conditioned variables to be
     !    solved.
     ! USED BY
     !    chm_element_operations_CMC
     !***
     !-----------------------------------------------------------------------
     use def_chemic,         only:  nZ_CMC_chm, diff_Z_CMC_chm, Yk_CMC_chm, &
                                    enthalp_CMC_chm, deriv2_Yk_CMC_chm, deriv2_enthalp_CMC_chm, &
                                    kfl_solve_enth_CMC_chm

     implicit none
     integer(ip)                 :: imixf, ipoin, iclas
     real(rp)                    :: numerator, denominator_inv(nZ_CMC_chm-2_ip), aux_var

     !
     ! This second derivative is of second order if the nodes are equally spaced and
     ! first order otherwise.
     !

     if( INOTMASTER ) then
        do imixf = 1, nZ_CMC_chm-2_ip
           aux_var = diff_Z_CMC_chm(imixf) * diff_Z_CMC_chm(imixf+1_ip) &
                   * (diff_Z_CMC_chm(imixf) + diff_Z_CMC_chm(imixf+1_ip))
           denominator_inv(imixf) = 1.0_rp / aux_var
        end do

        do imixf = 2, nZ_CMC_chm-1_ip
           do ipoin = 1, npoin
              do iclas = 1, nclas_chm
                 numerator = Yk_CMC_chm(imixf+1_ip,ipoin,iclas) * diff_Z_CMC_chm(imixf-1_ip) &
                            + Yk_CMC_chm(imixf-1_ip,ipoin,iclas) * diff_Z_CMC_chm(imixf) &
                            - Yk_CMC_chm(imixf,ipoin,iclas) * (diff_Z_CMC_chm(imixf-1)+diff_Z_CMC_chm(imixf))
                 deriv2_Yk_CMC_chm(imixf,ipoin,iclas) = 2.0_rp * numerator * denominator_inv(imixf-1_ip)
              end do
           end do
        end do

        if (kfl_solve_enth_CMC_chm /= 0_ip) then
           do imixf = 2, nZ_CMC_chm-1_ip
              do ipoin = 1, npoin
                 numerator = enthalp_CMC_chm(imixf+1_ip,ipoin) * diff_Z_CMC_chm(imixf-1_ip) &
                               + enthalp_CMC_chm(imixf-1_ip,ipoin) * diff_Z_CMC_chm(imixf) &
                               - enthalp_CMC_chm(imixf,ipoin) * (diff_Z_CMC_chm(imixf-1_ip)+diff_Z_CMC_chm(imixf))
                 deriv2_enthalp_CMC_chm(imixf,ipoin) = 2.0_rp * numerator * denominator_inv(imixf-1_ip)
              end do
           end do
        end if
     end if
  end subroutine chm_calc_diff_condVar_mixfraction_CMC



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! C O M P U T E   T I M E   S T E P !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine chm_updtcc_CMC(dtmin)
     !-----------------------------------------------------------------------
     !****f* Chemic/chm_updtcc_CMC
     ! NAME
     !    chm_updtcc_CMC
     ! DESCRIPTION
     !    This routine computes the critical time step for CMC model.
     ! USED BY
     !    chm_updtcc_CMC
     !***
     !-----------------------------------------------------------------------
     use def_master,         only : densi_gp
     use def_domain,         only : ntens, elmar, nelem, ltype, nnode, ngaus, lorde, ltopo, hnatu, lnods
     use def_chemic,         only : nZ_CMC_chm, kfl_advec_chm, kfl_ellen_chm, &
                                    condu_gp_CMC_chm, sphec_gp_CMC_chm, &
                                    spvol_gp_CMC_chm, &
                                    kfl_mesh_interp_CMC_chm, kfl_avg_cond_CMC_chm, &
                                    kfl_incl_PDF_trans_CMC_chm, &
                                    ADR_chm, Zs_CMC_chm, Z_CMC_chm, nvar_CMC_chm, &
                                    kfl_solve_enth_CMC_chm, Le_k
     use mod_ADR,            only : ADR_critical_time_step
     use mod_ADR,            only : mreac_adr
     use mod_communications, only : PAR_MIN

     implicit none
     real(rp),   intent(inout)    :: dtmin

     real(rp)                     :: dtcri(2)
     integer(ip)                  :: ielem,iclas,imixf,igaus                   ! Indices and dimensions
     integer(ip)                  :: pelty,pnode
     integer(ip)                  :: pgaus,plapl,porde,ptopo

     real(rp)                     :: PDF_val
     real(rp)                     :: dummr(ndime,mnode), kcp_la
     real(rp)                     :: chale(3),chave(3),hleng(3),tragl(9)

     ! Matrices for elements (values at nodes)
     real(rp)                     :: elcod(ndime,mnode)                     ! Coordinates
     real(rp)                     :: elvel_CFD(ndime,mnode)                 ! Velocity from CFD
     real(rp)                     :: elZavg_CFD(mnode)                      ! Average mixture fraction from CFD
     real(rp)                     :: elZvar_CFD(mnode)                      ! Mixture fraction variance from CFD
     real(rp)                     :: elZgrad_CFD(ndime,mnode)               ! Average mixture fraction gradient from CFD
     real(rp)                     :: elturb_dif_CFD(mnode)                  ! Mass turbulent diffusion coefficient from CFD
     real(rp)                     :: eldumm_CMC(mnode)
     real(rp)                     :: eldumm2_CMC(mnode,nvar_CMC_chm)

     ! Matrices for elements (values at Gaussian points)
     real(rp)                     :: gpvol(mgaus)                           ! |J|*w
     real(rp)                     :: gprea(mgaus,mreac_adr)                 ! r
     real(rp)                     :: gpcar(ndime,mnode,mgaus)               ! dNk/dxj
     real(rp)                     :: gphes(ntens,mnode,mgaus)               ! dNk/dxidxj
     real(rp)                     :: gpdif(mgaus,nvar_CMC_chm)              ! D_k
     real(rp)                     :: gpvel_CFD(ndime,mgaus)                 ! Velocity
     real(rp)                     :: gpZavg_CFD(mgaus)                      ! Average mixture fraction
     real(rp)                     :: gpZvar_CFD(mgaus)                      ! Mixture fraction variance
     real(rp)                     :: gpZgrad_CFD(ndime,mgaus)               ! Average mixture fraction gradient
     real(rp)                     :: gpturb_dif_CFD(mgaus)                  ! Mass turbulent diffusion coefficient
     real(rp)                     :: gp_densi_PDF_CMC(mgaus)                ! Uncond. rho * probability density function
     real(rp)                     :: gp_veloc_CMC_chm(ndime,mgaus)          ! Conditional velocity
     real(rp)                     :: gp_diff_phys_spc(mgaus,nvar_CMC_chm)   ! Diffusion in physical space
     real(rp)                     :: gpdumm_CMC(mgaus)
     real(rp)                     :: gpdumm2_CMC(mgaus,nvar_CMC_chm)
     real(rp)                     :: gpPDF_param(7)                         ! PDF parameters at Gaussian points
     real(rp)                     :: factor0_vel(ndime)                     ! Constant in the conditional velocity model
     real(rp)                     :: factor1_vel(ndime)                     ! Slope in the conditional velocity model

     external                     :: elmlen
     external                     :: elmchl
     external                     :: elmcar

     if( INOTMASTER ) then
        mixt_fr: do imixf = 2, nZ_CMC_chm-1_ip

           call chm_global2local_CMC(imixf)

           elements: do ielem = 1,nelem
              !
              ! Element dimensions
              !
              pelty = ltype(ielem)
              pnode = nnode(pelty)
              pgaus = ngaus(pelty)
              porde = lorde(pelty)
              ptopo = ltopo(pelty)

              !
              ! Initialization
              !
              elcod            = 0.0_rp
              elvel_CFD        = 0.0_rp
              elZavg_CFD       = 0.0_rp
              elZvar_CFD       = 0.0_rp
              elZgrad_CFD      = 0.0_rp
              elturb_dif_CFD   = 0.0_rp
              eldumm_CMC      = 0.0_rp
              eldumm2_CMC      = 0.0_rp
              gpvol            = 0.0_rp
              gprea            = 0.0_rp
              gpcar            = 0.0_rp
              gphes            = 0.0_rp
              gpdif            = 0.0_rp
              gp_diff_phys_spc = 0.0_rp
              gpdumm_CMC       = 0.0_rp
              gpdumm2_CMC      = 0.0_rp
              gpvel_CFD        = 0.0_rp
              gpZavg_CFD       = 0.0_rp
              gpZvar_CFD       = 0.0_rp
              gpZgrad_CFD      = 0.0_rp
              gpturb_dif_CFD   = 0.0_rp
              gp_densi_PDF_CMC = 0.0_rp
              gp_veloc_CMC_chm = 0.0_rp
              gpPDF_param      = 0.0_rp
              gpPDF_param(4)   = Zs_CMC_chm
              factor0_vel      = 0.0_rp
              factor1_vel      = 0.0_rp

              !
              ! Gather values at the element
              !
               call chm_elmgac_CMC(&
                      pnode, lnods(1:pnode,ielem), elcod(:,1:pnode), elvel_CFD(:,1:pnode), &
                      elZavg_CFD(1:pnode), elZvar_CFD(1:pnode), elZgrad_CFD(:,1:pnode), &
                      elturb_dif_CFD(1:pnode), eldumm_CMC(1:pnode), eldumm2_CMC(1:pnode,:), imixf)

              !
              ! CHALE, HLENG and TRAGL
              !
              call elmlen(&
                   ndime,pnode,elmar(pelty)%dercg,tragl,elcod(:,1:pnode),hnatu(pelty),hleng)

              call elmchl(&
                   tragl,hleng,elcod(:,1:pnode),dummr(:,1:pnode),chave,chale,pelty,pnode, &
                   porde,hnatu(pelty),kfl_advec_chm,kfl_ellen_chm)

              !
              ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
              !
              call elmcar(&
                   pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
                   elmar(pelty)%deriv,elmar(pelty)%heslo,elcod(:,1:pnode), &
                   gpvol(1:pgaus),gpcar(:,:,1:pgaus),gphes(:,:,1:pgaus),ielem)

              !
              ! Send quantities to Gaussian points
              !
              call chm_elmpre_CMC(&
                      pnode, pgaus, elmar(pelty)%shape, elvel_CFD(:,1:pnode), &
                      elZavg_CFD(1:pnode), elZvar_CFD(1:pnode), &
                      elZgrad_CFD(:,1:pnode), elturb_dif_CFD(1:pnode), &
                      eldumm_CMC(1:pnode), eldumm2_CMC(1:pnode,:), gpvel_CFD(:,1:pgaus), &
                      gpZavg_CFD(1:pgaus), gpZvar_CFD(1:pgaus), &
                      gpZgrad_CFD(:,1:pgaus), gpturb_dif_CFD(1:pgaus), &
                      gpdumm_CMC(1:pgaus), gpdumm2_CMC(1:pgaus,:))

              ! Compute the laminar diffusion coefficient D at Gaussian points
              do igaus = 1, pgaus
                 kcp_la = condu_gp_CMC_chm(ielem)%a(igaus,imixf,1) * spvol_gp_CMC_chm(ielem)%a(igaus,imixf,1) / &
                                 sphec_gp_CMC_chm(ielem) % a(igaus,imixf,1)
                 do iclas = 1,nspec_chm
                    gpdif(igaus,iclas) = kcp_la / Le_k(iclas)  ! It contains rho*D and not only D
                 end do
                 if (kfl_solve_enth_CMC_chm /= 0_ip) then  ! In case enthalpy is transported
                    gpdif(igaus,nvar_CMC_chm) = kcp_la
                 end if
              end do


              !
              ! Find values at Gauss points
              !
              do igaus = 1, pgaus
                 gp_densi_PDF_CMC(igaus) = densi_gp(ielem) % a(igaus,1,1)
              end do

              if (kfl_incl_PDF_trans_CMC_chm == 1_ip) then
                 do igaus = 1, pgaus
                    ! Compute PDF value
                    call chm_PDF_calc_CMC(gpZavg_CFD(igaus), gpZvar_CFD(igaus), imixf, PDF_val)

                    ! Find density times PDF
                    gp_densi_PDF_CMC(igaus) = gp_densi_PDF_CMC(igaus) * PDF_val

                 end do
              end if

              ! Diffusion in physical space
              do iclas = 1,nvar_CMC_chm
                 gp_diff_phys_spc(1:pgaus,iclas) = gp_densi_PDF_CMC(1:pgaus) * &
                                                  (gpdif(1:pgaus,iclas) + gpturb_dif_CFD(1:pgaus))
              end do

              ! Conditional velocity
              if (kfl_mesh_interp_CMC_chm == 1_ip .and. kfl_avg_cond_CMC_chm == 1_ip) then
                 gp_veloc_CMC_chm(1:ndime,1:pgaus) = gpvel_CFD(1:ndime,1:pgaus)
              else
                 do igaus = 1, pgaus
                    call chm_find_veloc_parameters_CMC(gpvel_CFD(1:ndime,igaus), &
                              factor0_vel(1:ndime), factor1_vel(1:ndime))
                    gp_veloc_CMC_chm(1:ndime,igaus) = factor0_vel(1:ndime) + factor1_vel(1:ndime) * Z_CMC_chm(imixf)
                 end do
              end if

              do iclas = 1, nvar_CMC_chm
                 ! Compute time-step

                 call ADR_critical_time_step(ADR_chm(iclas),gp_densi_PDF_CMC(1:pgaus), &
                        gp_veloc_CMC_chm(1:ndime,1:pgaus),gp_diff_phys_spc(1:pgaus,iclas), &
                        gprea(1:pgaus,:),dtcri,chale(1),chale(2))
                 ! Take minimum time step
                 dtmin = min(dtmin,dtcri(1))
              end do

           end do elements

        end do mixt_fr

     end if
     !
     ! Look for minimum over subdomains
     !
     call PAR_MIN(dtmin,'IN MY CODE')

  end subroutine chm_updtcc_CMC



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! C O M P U T E   V A R I A B L E S   O F   I N T E R E S T !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine chm_calc_temp_CMC(imixf)
    !-----------------------------------------------------------------------
    !****f* chemic/chm_calc_temp_CMC
    ! NAME
    !    chm_calc_temp_CMC
    ! DESCRIPTION
    !    It computes conditional values for temperature at nodes.
    ! USES
    ! USED BY
    !-----------------------------------------------------------------------
    use def_master,      only : therm, conce
    use def_chemic,      only : temp_CMC_chm

    implicit none

    integer(ip), intent(in)  :: imixf
    integer(ip)              :: ipoin
    real(rp)                 :: T, Tini, h_species

    if (INOTMASTER) then
       !
       ! Loop over points
       !

       points: do ipoin = 1, npoin
          ! Find conditional temperature
          Tini = temp_CMC_chm(imixf,ipoin)  ! Initial temperature from which start iteration
          T    = Tini                       ! Assignement to not leave an empty value
          h_species = therm(ipoin,1)
          !!!if (kfl_soot_chm /= 0) then
          !!!   sum_Yk_soot = 0.0_rp
          !!!   do iclas = nspec_chm + 1_ip, nclas_chm
          !!!      sum_Yk_soot = sum_Yk_soot + conce(ipoin,iclas,1)
          !!!   end do
          !!!   h_species = h_species - sum_Yk_soot * h_soot  !!!!! QUITAR
          !!!end if
          call calc_T_from_hY(conce(ipoin,1:nclas_chm,1), h_species, &
                 Tini, T)
          temp_CMC_chm(imixf,ipoin) = T
       end do points
    end if
  end subroutine chm_calc_temp_CMC



  subroutine calc_T_from_hY(Yk,h,Tini,T)
    !-----------------------------------------------------------------------
    !****f* chemic/calc_T_from_hY
    ! NAME
    !    calc_T_from_hY
    ! DESCRIPTION
    !    It computes temperature from enthalpy and species mass fractions.
    ! USES
    ! USED BY
    !-----------------------------------------------------------------------
#ifdef CANTERA
    use def_chemic,      only : coeff_cp_k, W_k
#endif
    use mod_physics,     only : physics_H_2_TCp

    implicit none

    real(rp), intent(in)     :: Yk(nclas_chm), h, Tini
    real(rp), intent(out)    :: T
#ifdef CANTERA
    integer(ip)              :: ivalu
    real(rp)                 :: cploc(6,2), aux_cp_lt(8), aux_cp_ht(8)
    real(rp)                 :: aux_cp

    external                 :: cantera_alya_cp
#endif
    external                 :: runend

#ifndef CANTERA
    ! Early exit until alternatives are added
    call runend('calc_T_from_hY: Please, build with support for cantera')
    if (.false.) print *, Yk,h,Tini ! avoid compiler warning (unused-dummy-argument)
    if (.false.) T = 0.0_rp         ! avoid intel compiler warning (#6843)
#endif

#ifdef CANTERA
    ! Find temperature from enthalpy and species mass fractions
    !!! ¿DESCOMENTAR?
    !!!alpha = 0.0_rp
    !!!do ii = 1, nspec_chm
    !!!   alpha = alpha + Yk(ii)
    !!!end do
    !!!Yk_aux = Yk(1:nspec_chm) / alpha
    !!!call cantera_alya_cp(nspec_chm,coeff_cp_k, Yk_aux, &
    call cantera_alya_cp(nspec_chm,coeff_cp_k, Yk, &
              W_k, aux_cp_lt, aux_cp_ht)
    do ivalu = 1, 6
       cploc(ivalu,1) = aux_cp_lt(ivalu)
       cploc(ivalu,2) = aux_cp_ht(ivalu)
    end do
    T = Tini
    call physics_H_2_TCp(h, cploc, T, aux_cp)
#endif
  end subroutine calc_T_from_hY


  subroutine calc_h_cp_from_TY(Yk,T,h,cp)
    !-----------------------------------------------------------------------
    !****f* chemic/calc_h_cp_from_TY
    ! NAME
    !    calc_h_cp_from_TY
    ! DESCRIPTION
    !    It computes enthalpy and specific heat from temperature and species
    !    mass fractions.
    ! USES
    ! USED BY
    !-----------------------------------------------------------------------
#ifdef CANTERA
    use def_chemic,            only :  coeff_cp_k, W_k
#endif
    use mod_physics,           only :  physics_T_2_HCp

    implicit none

    real(rp), intent(in)            :: Yk(nclas_chm), T
    real(rp), intent(out)           :: h
    real(rp), intent(out)           :: cp
#ifdef CANTERA
    real(rp)                        :: cploc(6,2), aux_cp_lt(8), aux_cp_ht(8)
    integer(ip)                     :: ivalu

    external                        :: cantera_alya_cp
#endif
    external                        :: runend

#ifndef CANTERA
    ! Early exit until alternatives are added
    call runend('calc_h_cp_from_TY: Please, build with support for cantera')
    if (.false.) print *, Yk,T            ! avoid compiler warning (unused-dummy-argument)
    if (.false.) cp = 0.0_rp; h = 0.0_rp  ! avoid intel compiler warning (#6843)
#endif

#ifdef CANTERA
    ! Find enthalpy from temperature and species mass fractions
    !!! ¿DESCOMENTAR?
    !!!alpha = 0.0_rp
    !!!do ii = 1, nspec_chm
    !!!   alpha = alpha + Yk(ii)
    !!!end do
    !!!Yk_aux = Yk(1:nspec_chm) / alpha
    !!!call cantera_alya_cp(nspec_chm,coeff_cp_k, Yk_aux, &
    call cantera_alya_cp(nspec_chm,coeff_cp_k, Yk, &
              W_k, aux_cp_lt, aux_cp_ht)
    do ivalu = 1, 6
       cploc(ivalu,1) = aux_cp_lt(ivalu)
       cploc(ivalu,2) = aux_cp_ht(ivalu)
    end do
    call physics_T_2_HCp(T, cploc, h, cp)
#endif

    !!!if (kfl_soot_chm /= 0) then
    !!!   do ii = nspec_chm+1, nclas_chm
    !!!      h  = h + Yk(ii) * h_soot   !!!QUITAR
    !!!      cp = cp + Yk(ii) * cp_soot !!!QUITAR
    !!!   end do
    !!!end if
  end subroutine calc_h_cp_from_TY


  subroutine chm_update_properties_CMC
    !-----------------------------------------------------------------------
    !****f* chemic/chm_update_properties_CMC
    ! NAME
    !    chm_update_properties_CMC
    ! DESCRIPTION
    !    Get properties at Gauss points and update through the kernel the
    !    properties.
    ! USES
    ! USED BY
    !-----------------------------------------------------------------------
    use def_chemic,       only : nZ_CMC_chm, spvol_gp_CMC_chm, visco_gp_CMC_chm, &
                                 sphec_gp_CMC_chm, condu_gp_CMC_chm
    use def_domain,       only : nelem, ltype, ngaus, nnode, mgaus, mnode
    use mod_ker_updpro,   only : ker_updpro
    use mod_ker_proper,   only : ker_proper

    implicit none
    integer(ip)               :: imixf, ielem, igaus
    integer(ip)               :: pelty, pgaus, pnode
    integer(ip)               :: dumm1
    real(rp)                  :: dumm2(mnode,mgaus)
    real(rp)                  :: dumm3(ndime,mnode,mgaus)
    real(rp)                  :: gpden(mgaus)
    real(rp)                  :: gpvis(mgaus)
    real(rp)                  :: gpsph(mgaus)
    real(rp)                  :: gpcon(mgaus)

    mixt_fr: do imixf = 1, nZ_CMC_chm
       !
       ! Get conditional properties according to chemical composition
       !
       call chm_global2local_CMC(imixf)
       call chm_calc_properties_gauss_CMC

       !
       ! Update properties through the kernel
       !
       call ker_updpro

       !
       ! Transfer properties to CMC matrices (values may be updated according to
       ! the material)
       !

       ! Loop over elements
       elements: do ielem = 1, nelem

          ! Element dimensions
          pelty = ltype(ielem)
          if( pelty > 0 ) then
              pgaus = ngaus(pelty)
              pnode = nnode(pelty)

              ! Get values from the kernel. The matrices in CMC conditional properties are
              ! treated in Gaussian points. The following is just a transfer of data and neither
              ! shape function nor derivatives are required
              call ker_proper('DENSI','PGAUS',dumm1,ielem,gpden,pnode,pgaus,dumm2,dumm3)
              call ker_proper('VISCO','PGAUS',dumm1,ielem,gpvis,pnode,pgaus,dumm2,dumm3)
              call ker_proper('SPHEA','PGAUS',dumm1,ielem,gpsph,pnode,pgaus,dumm2,dumm3)
              call ker_proper('CONDU','PGAUS',dumm1,ielem,gpcon,pnode,pgaus,dumm2,dumm3)

              do igaus = 1, pgaus
                 spvol_gp_CMC_chm(ielem) % a(igaus,imixf,1) = 1.0_rp / gpden(igaus)
                 visco_gp_CMC_chm(ielem) % a(igaus,imixf,1) = gpvis(igaus)
                 sphec_gp_CMC_chm(ielem) % a(igaus,imixf,1) = gpsph(igaus)
                 condu_gp_CMC_chm(ielem) % a(igaus,imixf,1) = gpcon(igaus)
              end do
          end if

       end do elements

    end do mixt_fr

  end subroutine chm_update_properties_CMC


  subroutine chm_calc_properties_gauss_CMC
    !-----------------------------------------------------------------------
    !****f* chemic/chm_calc_properties_gauss_CMC
    ! NAME
    !    chm_calc_properties_gauss_CMC
    ! DESCRIPTION
    !    It computes properties according to the chemical composition.
    ! USES
    !    chm_update_properties_CMC
    ! USED BY
    !-----------------------------------------------------------------------
#ifdef CANTERA
    use def_master,      only : prthe, visco_gp, densi_gp, condu_gp, sphec_gp
    use def_kermod,      only : gasco
    use def_chemic,      only : gas_chm
    use def_domain,      only : nelem,ltype,nnode,ngaus,lnods,elmar
    use cantera,         only : viscosity, thermalConductivity, setState_TPX, setMassFractions
#endif

    implicit none

#ifdef CANTERA
    integer(ip)             :: ielem,igaus,inode,iclas
    integer(ip)             :: pgaus
    integer(ip)             :: pelty,pnode
    real(rp)                :: gpcon(mgaus,nclas_chm)
    real(rp)                :: gph(mgaus)
    real(rp)                :: gptem(mgaus)
    real(rp)                :: elcon(mnode,nclas_chm)
    real(rp)                :: elcod(ndime,mnode)
    real(rp)                :: elh(mnode)
    real(rp)                :: eltem(mnode)
    real(rp)                :: aux_wmean_mixf
    real(rp)                :: dummy
#endif

    if (INOTMASTER) then
#ifdef CANTERA
       !
       ! Loop over elements
       !
       elements: do ielem = 1, nelem
          !
          ! Element dimensions
          !
          pelty = ltype(ielem)

          if( pelty > 0 ) then
              pnode = nnode(pelty)
              pgaus = ngaus(pelty)

              !
              ! Initialization variables
              !
              gpcon = 0.0_rp
              gph   = 0.0_rp
              gptem = 0.0_rp

              !
              ! Gather all
              !
              call chm_gatherProp_CMC( &
                       pnode,lnods(1:pnode,ielem),elcod(:,1:pnode),elcon(1:pnode,:), &
                       elh(1:pnode),eltem(1:pnode))

              !
              ! Species mass fraction Y_k at Gauss points
              !
              do iclas = 1,nclas_chm
                 do igaus = 1,pgaus
                    do inode = 1,pnode
                       gpcon(igaus,iclas) = gpcon(igaus,iclas)&
                                            + elmar(pelty)%shape(inode,igaus) * elcon(inode,iclas)
                    end do
                 end do
              end do

              !
              ! Enthalpy and temperature at Gauss points
              !
              do igaus = 1,pgaus
                 do inode = 1,pnode
                    gph(igaus)   = gph(igaus) &
                                         + elmar(pelty)%shape(inode,igaus) * elh(inode)
                    gptem(igaus) = gptem(igaus) &
                                         + elmar(pelty)%shape(inode,igaus) * eltem(inode)
                 end do
              end do

              !
              ! Compute transport properties
              ! Cantera properties
              !
              gauss: do igaus = 1, pgaus

                 ! Molecular weight
                 call compute_molWeight_mixt(gpcon(igaus,1:nspec_chm),aux_wmean_mixf)

                 ! Density
                 densi_gp(ielem) % a(igaus,1,1) = prthe(1) * aux_wmean_mixf / &
                                                     gasco / gptem(igaus)

                 ! Specific heat
                 call calc_h_cp_from_TY(gpcon(igaus,:),gptem(igaus),dummy, &
                         sphec_gp(ielem) % a(igaus,1,1))

                 ! Cantera functions
                 !!! ¿DESCOMENTAR?
                 !!!alpha = 0.0_rp
                 !!!do ii = 1, nspec_chm
                 !!!   alpha = alpha + gpcon(igaus,ii)
                 !!!end do
                 !!!gpcon(igaus,:) = gpcon(igaus,:) / alpha
                 call setState_TPX(gas_chm,gptem(igaus),prthe(1),gpcon(igaus,1:nspec_chm))
                 call setMassFractions(gas_chm,gpcon(igaus,1:nspec_chm))

                 ! Viscosity
                 visco_gp(ielem) % a(igaus,1,1) = viscosity(gas_chm)

                 ! Conductivity
                 condu_gp(ielem) % a(igaus,1,1) = thermalConductivity(gas_chm)

              end do gauss

          end if
       end do elements
#endif

    end if

  end subroutine chm_calc_properties_gauss_CMC


#ifdef CANTERA
  subroutine chm_gatherProp_CMC( &
                 pnode,lnods,elcod,elcon,elh,eltem)
     !------------------------------------------------------------------------
     !****f* Chemic/mod_chm_element_operations/chm_gatherProp_CMC
     ! NAME
     !    chm_gatherProp_CMC
     ! DESCRIPTION
     !    Gather operations for the combustion models
     ! USES
     ! USED BY
     !    chm_calc_properties_gauss_CMC
     !***
     !------------------------------------------------------------------------
     use def_master, only     :  conce,therm,tempe

     implicit none
     real(rp), parameter      :: Tmin = 200.0_rp, Tmax = 3000.0_rp
     integer(ip), intent(in)  :: pnode
     integer(ip), intent(in)  :: lnods(pnode)
     real(rp),    intent(out) :: elcod(ndime,pnode)
     real(rp),    intent(out) :: elcon(pnode,nclas_chm)
     real(rp),    intent(out) :: elh(pnode)
     real(rp),    intent(out) :: eltem(pnode)

     integer(ip)              :: inode,ipoin,iclas,idime

     !
     ! Initialization
     !
     elh    = 0.0_rp
     eltem  = 0.0_rp
     elcod  = 0.0_rp
     elcon  = 0.0_rp

     !
     ! Concentration and coordinates
     !
     do inode=1,pnode
        ipoin=lnods(inode)
        do iclas=1,nclas_chm
           elcon(inode,iclas) = conce(ipoin,iclas,1)
        end do

        do idime=1,ndime
           elcod(idime,inode)   = coord(idime,ipoin)
        end do

        elh(inode)   = therm(ipoin,1)
        eltem(inode) = max(Tmin, min(Tmax, tempe(ipoin,1)))
     end do

  end subroutine chm_gatherProp_CMC
#endif


  subroutine compute_molWeight_mixt(Yk,wmean_mixf)
     !-----------------------------------------------------------------------
     !****f* Chemic/compute_molWeight_mixt
     ! NAME
     !    compute_molWeight_mixt
     ! DESCRIPTION
     !    This routine computes the molecular weight for all the mixture
     !    fraction levels for a given physical point.
     ! USED BY
     !
     !***
     !-----------------------------------------------------------------------

     use def_chemic,             only :  W_k

     implicit none
     real(rp), intent(in)             :: Yk(nspec_chm)
     real(rp), intent(out)            :: wmean_mixf
     integer(ip)                      :: iclas

     wmean_mixf = 0.0_rp
     do iclas = 1, nspec_chm
        wmean_mixf =  wmean_mixf + Yk(iclas) / W_k(iclas)
     end do
     wmean_mixf = 1.0_rp / wmean_mixf

  end subroutine compute_molWeight_mixt



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! I N T E G R A T I O N S   P D F s !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine chm_integrate_var_CFD_CMC(nvar)

     !-----------------------------------------------------------------------
     !****f* Chemic/chm_integrate_var_CFD_CMC
     ! NAME
     !    chm_integrate_var_CFD_CMC
     ! DESCRIPTION
     !    Integration of density, viscosity, etc. when conditional variables
     !    are transferred to CFD from CMC.
     ! USED BY
     !
     !***
     !-----------------------------------------------------------------------

     use def_chemic,         only :  rec_CFD_from_CMC_chm, hrr_chm, nZ_CMC_chm, &
                                     rec_CFD_old_from_CMC_chm, freq_CFD_coup_CMC_chm
     use def_master,         only :  densi, visco, condk, sphek, tempe, conce, ittim

     implicit none
     integer(ip), intent(in)      :: nvar
     integer(ip)                  :: ipoin, imixf
     real(rp)                     :: aux_var_cond(nZ_CMC_chm,nvar), aux_var_cond_old(nZ_CMC_chm,nvar)
     real(rp)                     :: aux_var_uncond(nvar), alpha

     aux_var_cond   = 0.0_rp
     aux_var_uncond = 0.0_rp

     if (freq_CFD_coup_CMC_chm > 1_ip) then
        if (mod(ittim, freq_CFD_coup_CMC_chm) == 0_ip) then
           alpha = 1.0_rp
        else
           alpha = real(mod(ittim, freq_CFD_coup_CMC_chm),rp) / real(freq_CFD_coup_CMC_chm,rp)
        end if
     end if

     do ipoin = 1, npoin
        do imixf = 1, nZ_CMC_chm
           aux_var_cond(imixf,1:nvar) = rec_CFD_from_CMC_chm(1:nvar,imixf,ipoin)
        end do

        if (freq_CFD_coup_CMC_chm > 1_ip) then
           do imixf = 1, nZ_CMC_chm
              aux_var_cond_old(imixf,1:nvar) = rec_CFD_old_from_CMC_chm(1:nvar,imixf,ipoin)
           end do

           aux_var_cond = aux_var_cond_old + alpha * (aux_var_cond - aux_var_cond_old)
        end if

        ! Integrations
        call chm_mxt_fr_integr_previous_steps_CMC(nvar, aux_var_cond, &
                conce(ipoin,3,1), conce(ipoin,4,1), &
                aux_var_uncond)

        ! Assignation
        densi(ipoin,1) = 1.0_rp / aux_var_uncond(1)  ! CAUTION: the specific volume is sent from CMC
        visco(ipoin,1) = aux_var_uncond(2)
        condk(ipoin,1) = aux_var_uncond(3)
        sphek(ipoin,1) = aux_var_uncond(4)
        tempe(ipoin,1) = aux_var_uncond(5)
        hrr_chm(ipoin) = aux_var_uncond(6) * densi(ipoin,1)  ! The heat release per mass unit is obtained from CMC

     end do

  end subroutine chm_integrate_var_CFD_CMC


  subroutine chm_integrate_flow_var_points_CMC

     !-----------------------------------------------------------------------
     !****f* Chemic/chm_integrate_flow_var_points_CMC
     ! NAME
     !    chm_integrate_flow_var_points_CMC
     ! DESCRIPTION
     !    This routine integrates density, viscosity, heat release and if required
     !    mass fractions, enthalpy and temperature at nodes.
     ! USED BY
     !
     !***
     !-----------------------------------------------------------------------

     use def_master,              only :  postp, wmean, ittim
     use def_chemic,              only :  hrr_mass_cond_CMC_chm, nZ_CMC_chm, temp_int_CMC_chm, &
                                          Yk_CMC_chm, enthalp_CMC_chm, temp_CMC_chm, &
                                          hrr_chm, Yk_int_CMC_chm, enthalp_int_CMC_chm, &
                                          src_Yk_CMC_chm, src_Yk_int_CMC_chm, &
                                          Zavg_CMC_chm, Zvar_CMC_chm, &
                                          kfl_solve_enth_CMC_chm

     implicit none
     integer(ip)                       :: ipoin, imixf, iclas
     integer(ip)                       :: nvar_int, integ_all
     real(rp), allocatable             :: aux_integrals_CMC_chm(:)
     real(rp), allocatable             :: aux_var(:,:)


     if( INOTMASTER ) then

        !!!!!!! ÉSTE ES EL CORRECTO
        !if( mod(ittim, postp(1) % npp_stepi(arrays_number('MASCN') ) == 0 ) then   ! 56 is the position for conditional mass fractions

        !!!!!!! PROVISIONAL -> BORRAR
        nvar_int  = 3_ip
        integ_all = 0_ip
        if(postp(1) % npp_stepi(arrays_number('MASUN'),0) /= 0) then
           if ( mod(ittim, postp(1) % npp_stepi(arrays_number('MASUN'),0) ) == 0_ip ) then   ! arrays_number('MASUN') is for
                                                                                             ! unconditional mass fractions
              ! When writing integrate all variables
              nvar_int  = 2_ip*nclas_chm+4_ip
              integ_all = 1_ip
           end if
        end if

        allocate(aux_integrals_CMC_chm(nvar_int))
        allocate(aux_var(nZ_CMC_chm,nvar_int))

        do ipoin = 1,npoin
           aux_var(1:nZ_CMC_chm,1_ip) = temp_CMC_chm(1:nZ_CMC_chm,ipoin)
           ! Compute cp, conductivity and molecular weight at nodes
           do imixf = 1, nZ_CMC_chm
              call compute_molWeight_mixt(Yk_CMC_chm(imixf,ipoin,1:nspec_chm),aux_var(imixf,2_ip))
           end do
           aux_var(1:nZ_CMC_chm,3_ip) = hrr_mass_cond_CMC_chm(1:nZ_CMC_chm,ipoin)

          if (integ_all == 1_ip) then
             ! Integrate all the variables
             do iclas = 1, nclas_chm
                aux_var(1:nZ_CMC_chm,3_ip+iclas)           = Yk_CMC_chm(1:nZ_CMC_chm,ipoin,iclas)
                aux_var(1:nZ_CMC_chm,3_ip+iclas+nclas_chm) = src_Yk_CMC_chm(1:nZ_CMC_chm,ipoin,iclas)
             end do

             if (kfl_solve_enth_CMC_chm == 0) then
                aux_var(1:nZ_CMC_chm,2_ip*nclas_chm+4_ip) = enthalp_CMC_chm(1:nZ_CMC_chm,1)
             else
                aux_var(1:nZ_CMC_chm,2_ip*nclas_chm+4_ip) = enthalp_CMC_chm(1:nZ_CMC_chm,ipoin)
             end if
          end if

          ! Integrations
          call chm_mxt_fr_integr_previous_steps_CMC(nvar_int, aux_var, &
                 Zavg_CMC_chm(ipoin), Zvar_CMC_chm(ipoin), &
                 aux_integrals_CMC_chm)

          ! Assign values
          temp_int_CMC_chm(ipoin) = aux_integrals_CMC_chm(1)
          wmean(ipoin,1)          = aux_integrals_CMC_chm(2)
          hrr_chm(ipoin)          = aux_integrals_CMC_chm(3)

          if (integ_all == 1_ip) then
             Yk_int_CMC_chm(ipoin,1:nclas_chm)     = aux_integrals_CMC_chm(4:nclas_chm+3_ip)
             src_Yk_int_CMC_chm(ipoin,1:nclas_chm) = aux_integrals_CMC_chm(nclas_chm+4:2_ip*nclas_chm+3_ip)
             enthalp_int_CMC_chm(ipoin)            = aux_integrals_CMC_chm(2_ip*nclas_chm+4_ip)
          end if

        end do

        deallocate(aux_integrals_CMC_chm)
        deallocate(aux_var)

     end if

  end subroutine chm_integrate_flow_var_points_CMC



  subroutine chm_integrate_flow_var_gauss_CMC

     !-----------------------------------------------------------------------
     !****f* Chemic/chm_integrate_flow_var_gauss_CMC
     ! NAME
     !    chm_integrate_flow_var_gauss_CMC
     ! DESCRIPTION
     !    This routine integrates density, viscosity, heat release and if required
     !    mass fractions and source terms, enthalpy and temperature at Gauss points.
     ! USED BY
     !
     !***
     !-----------------------------------------------------------------------

     use def_master,              only :  densi_gp, visco_gp, condu_gp, sphec_gp, &
                                          postp, ittim
     use def_domain,              only :  nelem,ltype,nnode,ngaus,elmar,lnods
     use def_chemic,              only :  nZ_CMC_chm, spvol_gp_CMC_chm, &
                                          visco_gp_CMC_chm, condu_gp_CMC_chm, &
                                          sphec_gp_CMC_chm, kfl_post_gp_CMC_chm, &
                                          react_scalars_gp_CMC_chm

     implicit none
     integer(ip)                       :: ielem, pelty, pnode, pgaus
     integer(ip)                       :: imixf, igaus, ivar
     integer(ip)                       :: nvar_int, integ_all
     real(rp)                          :: elZavg(mnode)
     real(rp)                          :: elZvar(mnode)

     real(rp)                          :: gpZavg(mgaus)
     real(rp)                          :: gpZvar(mgaus)
     real(rp), allocatable             :: gpvar_int(:,:,:)
     real(rp), allocatable             :: gp_aux(:)


     if (INOTMASTER) then
        if (kfl_post_gp_CMC_chm == 1_ip) then
           nvar_int  = 7_ip
           integ_all = 0_ip
           if(postp(1) % npp_stepi(arrays_number('MASUN'),0) /= 0) then
              if ( mod(ittim, postp(1) % npp_stepi(arrays_number('MASUN'),0) ) == 0_ip ) then   ! arrays_number('MASUN') is for
                                                                                                ! unconditional mass fractions
                 nvar_int  = 2_ip * nclas_chm + 8_ip  ! 8 instead of 6 since properties are included
                 integ_all = 1_ip
              end if
           end if
        else
           nvar_int  = 4_ip  ! Density, viscosity, specific heat, conductivity
           integ_all = 0_ip
        end if

        allocate(gpvar_int(nZ_CMC_chm,nvar_int,mgaus))
        allocate(gp_aux(nvar_int))
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

               ! Initialization
               gpZavg    = 0.0_rp
               gpZvar    = 0.0_rp
               gpvar_int = 0.0_rp
               gp_aux    = 0.0_rp

               ! Gather mixing variables
               call chm_gatherIntegral_mixt_var(pnode,lnods(1:pnode,ielem), &
                       elZavg(1:pnode),elZvar(1:pnode))

               ! Send mixing variables to Gaussian points
               call chm_elmpreIntegral_mixt_var(pnode,pgaus,elmar(pelty)%shape, &
                       elZavg(1:pnode),elZvar(1:pnode),gpZavg(1:pgaus),gpZvar(1:pgaus))

               ! Gather reacting variables
               do imixf = 1,nZ_CMC_chm
                  do igaus = 1,pgaus
                     gpvar_int(imixf,1_ip,igaus) = spvol_gp_CMC_chm(ielem) % a(igaus,imixf,1)
                     gpvar_int(imixf,2_ip,igaus) = visco_gp_CMC_chm(ielem) % a(igaus,imixf,1)
                     gpvar_int(imixf,3_ip,igaus) = sphec_gp_CMC_chm(ielem) % a(igaus,imixf,1)
                     gpvar_int(imixf,4_ip,igaus) = condu_gp_CMC_chm(ielem) % a(igaus,imixf,1)
                  end do
               end do

               if(kfl_post_gp_CMC_chm == 1_ip) then
                  call chm_elmpreIntegral_react_var(nvar_int-4_ip, pgaus, &
                          pnode, lnods(1:pnode,ielem), elmar(pelty)%shape, &
                          integ_all, gpvar_int(:,5:nvar_int,1:pgaus))
               end if

               do igaus = 1,pgaus
                  ! Integrations in mixture fraction
                  call chm_mxt_fr_integr_previous_steps_CMC(nvar_int, &
                          gpvar_int(:,:,igaus), gpZavg(igaus), gpZvar(igaus), gp_aux)
                  densi_gp(ielem) % a(igaus,1,1) = 1.0_rp / gp_aux(1)
                  visco_gp(ielem) % a(igaus,1,1) = gp_aux(2)
                  sphec_gp(ielem) % a(igaus,1,1) = gp_aux(3)
                  condu_gp(ielem) % a(igaus,1,1) = gp_aux(4)

                  if(kfl_post_gp_CMC_chm == 1_ip) then
                     do ivar = 5, nvar_int
                        react_scalars_gp_CMC_chm(ielem) % a(ivar-4_ip,igaus,1) = gp_aux(ivar)
                     end do
                  end if
               end do

           end if
        end do elements

        deallocate(gpvar_int)
        deallocate(gp_aux)

     end if

  end subroutine chm_integrate_flow_var_gauss_CMC


  subroutine chm_gatherIntegral_mixt_var(pnode,lnods,elZavg,elZvar)

     !-----------------------------------------------------------------------
     !****f* Chemic/chm_gatherIntegral_mixt_var
     ! NAME
     !    chm_gatherIntegral_mixt_var
     ! DESCRIPTION
     !    Gather values for chm_integrate_flow_var_gauss_CMC.
     ! USED BY
     !    chm_integrate_flow_var_gauss_CMC
     !***
     !-----------------------------------------------------------------------

     use def_chemic,              only :  Zavg_CMC_chm, Zvar_CMC_chm

     implicit none
     integer(ip), intent(in)  :: pnode, lnods(pnode)
     real(rp),    intent(out) :: elZavg(pnode)
     real(rp),    intent(out) :: elZvar(pnode)

     integer(ip)              :: inode,ipoin

     !
     ! Initialization
     !
     elZavg = 0.0_rp
     elZvar = 0.0_rp

     do inode = 1,pnode
        ipoin = lnods(inode)
        elZavg(inode) = Zavg_CMC_chm(ipoin)
        elZvar(inode) = Zvar_CMC_chm(ipoin)
     end do

  end subroutine chm_gatherIntegral_mixt_var


  subroutine chm_elmpreIntegral_mixt_var(pnode,pgaus,gpsha,elZavg,elZvar,gpZavg,gpZvar)

     !------------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_elmpreIntegral_mixt_var
     ! NAME
     !    chm_elmpreIntegral_mixt_var
     ! DESCRIPTION
     !    Transfer elemental values to Gaussian points for the set of not
     !    solved variables that do not depend on the mixture fraction.
     ! USES
     ! USED BY
     !    chm_integrate_flow_var_gauss_CMC
     !***
     !------------------------------------------------------------------------

     implicit none
     integer(ip), intent(in)  :: pnode, pgaus
     real(rp),    intent(in)  :: gpsha(pnode,pgaus)
     real(rp),    intent(in)  :: elZavg(pnode)
     real(rp),    intent(in)  :: elZvar(pnode)
     real(rp),    intent(out) :: gpZavg(pgaus)
     real(rp),    intent(out) :: gpZvar(pgaus)

     integer(ip)              :: igaus, inode


     do igaus = 1,pgaus
        do inode = 1,pnode
           gpZavg(igaus)  = gpZavg(igaus) + gpsha(inode,igaus) * elZavg(inode)
           gpZvar(igaus)  = gpZvar(igaus) + gpsha(inode,igaus) * elZvar(inode)
        end do
     end do

  end subroutine chm_elmpreIntegral_mixt_var


  subroutine chm_elmpreIntegral_react_var(nvar, pgaus, pnode, lnods, gpsha, integ_all, &
                gpvar_int)
     !------------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_elmpreIntegral_react_var
     ! NAME
     !    chm_elmpreIntegral_react_var
     ! DESCRIPTION
     !    Obtain values at Gaussian points.
     ! USES
     ! USED BY
     !    chm_integrate_flow_var_gauss_CMC
     !***
     !------------------------------------------------------------------------
     use def_chemic,        only :  nZ_CMC_chm, Yk_CMC_chm, enthalp_CMC_chm, &
                                    temp_CMC_chm, hrr_mass_cond_CMC_chm, src_Yk_CMC_chm, &
                                    kfl_solve_enth_CMC_chm

     implicit none
     integer(ip), intent(in)     :: nvar
     integer(ip), intent(in)     :: pgaus
     integer(ip), intent(in)     :: pnode
     integer(ip), intent(in)     :: lnods(pnode)
     real(rp),    intent(in)     :: gpsha(pnode,pgaus)
     integer(ip), intent(in)     :: integ_all
     real(rp),    intent(out)    :: gpvar_int(nZ_CMC_chm, nvar, pgaus)
     real(rp)                    :: elm_Yk(nclas_chm,pnode), elm_Yk_src(nclas_chm,pnode)
     real(rp)                    :: elm_temp(pnode), elm_enthalp(pnode), elm_hrr(pnode)
     real(rp)                    :: gp_Yk(nclas_chm,pgaus), gp_Yk_src(nclas_chm,pgaus)
     real(rp)                    :: gp_temp(pgaus), gp_enthalp(pgaus), gp_hrr(pgaus)
     integer(ip)                 :: imixf, inode, ipoin, igaus

     do imixf = 1, nZ_CMC_chm

        ! Transfer information to nodes
        elm_Yk      = 0.0_rp
        elm_Yk_src  = 0.0_rp
        elm_temp    = 0.0_rp
        elm_enthalp = 0.0_rp
        elm_hrr     = 0.0_rp
        gp_Yk       = 0.0_rp
        gp_Yk_src   = 0.0_rp
        gp_temp     = 0.0_rp
        gp_enthalp  = 0.0_rp
        gp_hrr      = 0.0_rp

        do inode = 1,pnode
           ipoin = lnods(inode)
           elm_temp(inode)     = temp_CMC_chm(imixf,ipoin)
           elm_Yk(:,inode)     = Yk_CMC_chm(imixf,ipoin,:)
           elm_hrr(inode)      = hrr_mass_cond_CMC_chm(imixf,ipoin)
        end do

        if (integ_all == 1_ip) then
           do inode = 1,pnode
              ipoin = lnods(inode)
              if (kfl_solve_enth_CMC_chm == 0_ip) then
                 elm_enthalp(inode) = enthalp_CMC_chm(imixf,1)
              else
                 elm_enthalp(inode) = enthalp_CMC_chm(imixf,ipoin)
              end if
              elm_Yk_src(:,inode) = src_Yk_CMC_chm(imixf,ipoin,:)
           end do
        end if

        ! Transfer information to Gaussian points and compute variables
        do igaus = 1,pgaus
           do inode = 1,pnode
              gp_temp(igaus) = gp_temp(igaus) + gpsha(inode,igaus) * elm_temp(inode)
              gp_Yk(:,igaus) = gp_Yk(:,igaus) + gpsha(inode,igaus) * elm_Yk(:,inode)
              gp_hrr(igaus)  = gp_hrr(igaus)  + gpsha(inode,igaus) * elm_hrr(inode)
           end do
           gpvar_int(imixf, 1_ip, igaus) = gp_temp(igaus)
           call compute_molWeight_mixt(gp_Yk(:,igaus), gpvar_int(imixf,2_ip,igaus))
           gpvar_int(imixf, 3_ip, igaus) = gp_hrr(igaus)

           ! Assign values
           if (integ_all == 1_ip) then
              ! All the variables
              do inode = 1,pnode
                 gp_enthalp(igaus)   = gp_enthalp(igaus)  + gpsha(inode,igaus) * elm_enthalp(inode)
                 gp_Yk_src(:,igaus)  = gp_Yk_src(:,igaus) + gpsha(inode,igaus) * elm_Yk_src(:,inode)
              end do
              gpvar_int(imixf,4_ip:nclas_chm+3_ip,igaus)                = gp_Yk(:,igaus)
              gpvar_int(imixf,nclas_chm+4_ip:2_ip*nclas_chm+3_ip,igaus) = gp_Yk_src(:,igaus)
              gpvar_int(imixf,2_ip*nclas_chm+4_ip,igaus)                = gp_enthalp(igaus)
           end if
        end do
     end do

  end subroutine chm_elmpreIntegral_react_var


  subroutine chm_nodal_projection_CMC
  !------------------------------------------------------------------------
  ! NAME
  !****f* Chemic/mod_chm_operations_CMC/chm_nodal_projection_CMC
  ! DESCRIPTION
  !    Projection of unconditional properties and other postprocessing
  !    variables if necessary.
  ! USES
  ! USED BY
  !    chm_initialization_domain_CMC, chm_begite, chm_endite
  !***
  !------------------------------------------------------------------------

  use def_master,              only : postp, densi_gp, visco_gp, sphec_gp, condu_gp, wmean, ittim
  use mod_solver,              only : solver_lumped_mass_system
  use def_domain,              only : nelem,ltype,nnode,ngaus,lorde,elmar,lnods, &
                                      ntens
  use def_chemic,              only : densi_int_CMC_chm, visco_lam_int_CMC_chm, &
                                      condu_int_CMC_chm, sphea_int_CMC_chm, &
                                      temp_int_CMC_chm, Yk_int_CMC_chm, hrr_chm, &
                                      src_Yk_int_CMC_chm, enthalp_int_CMC_chm, &
                                      aux_val_CMC_chm, react_scalars_gp_CMC_chm, &
                                      kfl_post_gp_CMC_chm
  
  implicit none
  integer(ip)                      :: ielem,inode,ipoin,igaus, ivar, iclas
  integer(ip)                      :: pnode, pgaus, porde, aux_plapl,pelty
  integer(ip)                      :: lnods_loc(mnode)
  real(rp)                         :: elcod(ndime,mnode)
  real(rp)                         :: aux_gpcar(ndime,mnode,mgaus)    ! dNk/dxj
  real(rp)                         :: aux_gphes(ntens,mnode,mgaus)    ! dNk/dxidxj
  real(rp)                         :: gpvol(mgaus)
  real(rp)                         :: fact_densi, fact_visco
  real(rp)                         :: fact_sphec, fact_condu
  real(rp)                         :: fact_var
  real(rp)                         :: eldensi(mnode)
  real(rp)                         :: elvisco(mnode)
  real(rp)                         :: elsphec(mnode)
  real(rp)                         :: elcondu(mnode)
  real(rp), allocatable            :: elvar(:,:)

  integer(ip)                      :: nvar, integ_all

  external                         :: elmcar
  external                         :: runend

  if (INOTMASTER) then

     if (kfl_post_gp_CMC_chm == 1_ip) then
        nvar      = 3_ip
        integ_all = 0_ip
        if( postp(1) % npp_stepi(arrays_number('MASUN'),0) /= 0 ) then
           if ( mod(ittim, postp(1) % npp_stepi(arrays_number('MASUN'),0) ) == 0_ip ) then   ! arrays_number('MASUN') is for
                                                                                             ! unconditional mass fractions
              nvar      = 2_ip * nclas_chm + 4_ip
              integ_all = 1_ip
           end if
        end if
     else
        nvar      = 0_ip
        integ_all = 0_ip
     end if

     densi_int_CMC_chm     = 0.0_rp
     visco_lam_int_CMC_chm = 0.0_rp
     condu_int_CMC_chm     = 0.0_rp
     sphea_int_CMC_chm     = 0.0_rp
     if (kfl_post_gp_CMC_chm == 1_ip) then
        temp_int_CMC_chm = 0.0_rp
        wmean            = 0.0_rp
        hrr_chm          = 0.0_rp
        if (integ_all == 1_ip) then
           Yk_int_CMC_chm      = 0.0_rp
           src_Yk_int_CMC_chm  = 0.0_rp
           enthalp_int_CMC_chm = 0.0_rp
        end if
        allocate(elvar(mnode,nvar))
     end if


     !
     ! Loop over elements
     !
     elements: do ielem = 1, nelem
        !
        ! Element dimensions
        !
        pelty = ltype(ielem)

        if( pelty > 0 ) then
            pnode = nnode(pelty)
            pgaus = ngaus(pelty)
            porde = lorde(pelty)

            !
            ! Initialization
            !
            elcod     = 0.0_rp
            aux_plapl = 0_ip
            gpvol     = 0.0_rp
            aux_gpcar = 0.0_rp
            aux_gphes = 0.0_rp
            eldensi   = 0.0_rp
            elvisco   = 0.0_rp
            elsphec   = 0.0_rp
            elcondu   = 0.0_rp
            if (kfl_post_gp_CMC_chm == 1_ip)   elvar = 0.0_rp

            !
            ! Gather coordinates at the element
            !
            call chm_elmgac_coords(pnode, lnods(1:pnode,ielem), elcod(:,1:pnode))

            !
            ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
            !
            call elmcar(&
                 pnode,pgaus,aux_plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
                 elmar(pelty)%deriv,elmar(pelty)%heslo,elcod(:,1:pnode), &
                 gpvol(1:pgaus),aux_gpcar(:,:,1:pgaus),&
                 aux_gphes(:,:,1:pgaus),ielem)


            if( porde == 1 ) then
               !
               ! Element assembly
               !
               do igaus = 1, pgaus
                  fact_densi = gpvol(igaus) * densi_gp(ielem) % a(igaus,1,1)
                  fact_visco = gpvol(igaus) * visco_gp(ielem) % a(igaus,1,1)
                  fact_sphec = gpvol(igaus) * sphec_gp(ielem) % a(igaus,1,1)
                  fact_condu = gpvol(igaus) * condu_gp(ielem) % a(igaus,1,1)
                  do inode = 1,pnode
                     eldensi(inode) = eldensi(inode) + &
                          elmar(pelty)%shape(inode,igaus) * fact_densi
                     elvisco(inode) = elvisco(inode) + &
                          elmar(pelty)%shape(inode,igaus) * fact_visco
                     elsphec(inode) = elsphec(inode) + &
                          elmar(pelty)%shape(inode,igaus) * fact_sphec
                     elcondu(inode) = elcondu(inode) + &
                          elmar(pelty)%shape(inode,igaus) * fact_condu
                  end do

                  if (kfl_post_gp_CMC_chm == 1_ip) then
                     do ivar = 1, nvar
                        fact_var = gpvol(igaus) * react_scalars_gp_CMC_chm(ielem) % a(ivar,igaus,1)
                        do inode = 1,pnode
                           elvar(inode,ivar) = elvar(inode,ivar) + &
                              elmar(pelty)%shape(inode,igaus) * fact_var
                        end do
                     end do
                  end if
               end do

               !
               ! Nodal projection
               !
               lnods_loc(1:pnode) = lnods(1:pnode,ielem)
               do inode = 1,pnode
                  ipoin = lnods_loc(inode)
                  densi_int_CMC_chm(ipoin)     = densi_int_CMC_chm(ipoin)     + eldensi(inode)
                  visco_lam_int_CMC_chm(ipoin) = visco_lam_int_CMC_chm(ipoin) + elvisco(inode)
                  condu_int_CMC_chm(ipoin)     = condu_int_CMC_chm(ipoin)     + elcondu(inode)
                  sphea_int_CMC_chm(ipoin)     = sphea_int_CMC_chm(ipoin)     + elsphec(inode)

                  if (kfl_post_gp_CMC_chm == 1_ip) then
                     temp_int_CMC_chm(ipoin) = temp_int_CMC_chm(ipoin) + elvar(inode,1)
                     wmean(ipoin,1)          = wmean(ipoin,1)          + elvar(inode,2)
                     hrr_chm(ipoin)          = hrr_chm(ipoin)          + elvar(inode,3)
                     if (integ_all == 1_ip) then
                        do iclas = 1, nclas_chm
                           Yk_int_CMC_chm(ipoin,iclas)     = Yk_int_CMC_chm(ipoin,iclas) + &
                              elvar(inode,3_ip+iclas)
                           src_Yk_int_CMC_chm(ipoin,iclas) = src_Yk_int_CMC_chm(ipoin,iclas) + &
                              elvar(inode,3_ip+nclas_chm+iclas)
                        end do
                        enthalp_int_CMC_chm(ipoin) = enthalp_int_CMC_chm(ipoin) + &
                           elvar(inode,2_ip*nclas_chm+4_ip)
                     end if
                  end if
               end do
            else
               call runend('CHEMIC OPERATIONS MOD. CMC: projections for density and viscosity not implemented for orders higher&
                   & than 1')
            end if

        end if
     end do elements

     call solver_lumped_mass_system(1_ip,densi_int_CMC_chm)
     call solver_lumped_mass_system(1_ip,visco_lam_int_CMC_chm)
     call solver_lumped_mass_system(1_ip,condu_int_CMC_chm)
     call solver_lumped_mass_system(1_ip,sphea_int_CMC_chm)
     if (kfl_post_gp_CMC_chm == 1_ip) then
        call solver_lumped_mass_system(1_ip,temp_int_CMC_chm)
        call solver_lumped_mass_system(1_ip,wmean)
        call solver_lumped_mass_system(1_ip,hrr_chm)
        if (integ_all == 1_ip) then
           do iclas = 1, nclas_chm
              ! Species mass fractions
              do ipoin = 1, npoin
                 aux_val_CMC_chm(ipoin) = Yk_int_CMC_chm(ipoin,iclas)
              end do
              call solver_lumped_mass_system(1_ip,aux_val_CMC_chm)
              do ipoin = 1, npoin
                 Yk_int_CMC_chm(ipoin,iclas) = aux_val_CMC_chm(ipoin)
              end do

              ! Chemical source terms for species
              do ipoin = 1, npoin
                 aux_val_CMC_chm(ipoin) = src_Yk_int_CMC_chm(ipoin,iclas)
              end do
              call solver_lumped_mass_system(1_ip,aux_val_CMC_chm)
              do ipoin = 1, npoin
                 src_Yk_int_CMC_chm(ipoin,iclas) = aux_val_CMC_chm(ipoin)
              end do
           end do
           call solver_lumped_mass_system(1_ip,enthalp_int_CMC_chm)
        end if
        deallocate(elvar)
     end if

  end if

  end subroutine chm_nodal_projection_CMC


  subroutine chm_nodal_projection_fundamental_CMC(imixf)
  !------------------------------------------------------------------------
  ! NAME
  !****f* Chemic/mod_chm_operations_CMC/chm_nodal_projection_fundamental_CMC
  ! DESCRIPTION
  !    Projection of conditional specific volume and transport properties.
  !    CAUTION: the specific volume is saved in densi
  ! USES
  ! USED BY
  !    chm_initialization_domain_CMC, chm_begite, chm_endite
  !***
  !------------------------------------------------------------------------

  use def_master,              only : condk, densi, sphek, visco
  use mod_solver,              only : solver_lumped_mass_system
  use def_domain,              only : nelem,ltype,nnode,ngaus,lorde,elmar,lnods, &
                                      ntens
  use def_chemic,              only : spvol_gp_CMC_chm, visco_gp_CMC_chm, &
                                      sphec_gp_CMC_chm, condu_gp_CMC_chm

  implicit none
  integer(ip), intent(in)          :: imixf
  integer(ip)                      :: ielem,inode,ipoin,igaus
  integer(ip)                      :: pnode, pgaus, porde, aux_plapl,pelty
  integer(ip)                      :: lnods_loc(mnode)
  real(rp)                         :: elcod(ndime,mnode)
  real(rp)                         :: aux_gpcar(ndime,mnode,mgaus)    ! dNk/dxj
  real(rp)                         :: aux_gphes(ntens,mnode,mgaus)    ! dNk/dxidxj
  real(rp)                         :: gpvol(mgaus)
  real(rp)                         :: fact_densi, fact_visco
  real(rp)                         :: fact_sphec, fact_condu
  real(rp)                         :: eldensi(mnode)
  real(rp)                         :: elvisco(mnode)
  real(rp)                         :: elsphec(mnode)
  real(rp)                         :: elcondu(mnode)

  external                         :: elmcar
  external                         :: runend

  if (INOTMASTER) then


     densi = 0.0_rp
     visco = 0.0_rp
     condk = 0.0_rp
     sphek = 0.0_rp

     !
     ! Loop over elements
     !
     elements: do ielem = 1, nelem
        !
        ! Element dimensions
        !
        pelty = ltype(ielem)

        if( pelty > 0 ) then
            pnode = nnode(pelty)
            pgaus = ngaus(pelty)
            porde = lorde(pelty)

            !
            ! Initialization
            !
            elcod     = 0.0_rp
            aux_plapl = 0_ip
            gpvol     = 0.0_rp
            aux_gpcar = 0.0_rp
            aux_gphes = 0.0_rp
            eldensi   = 0.0_rp
            elvisco   = 0.0_rp
            elsphec   = 0.0_rp
            elcondu   = 0.0_rp

            !
            ! Gather coordinates at the element
            !
            call chm_elmgac_coords(pnode, lnods(1:pnode,ielem), elcod(:,1:pnode))

            !
            ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
            !
            call elmcar(&
                 pnode,pgaus,aux_plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
                 elmar(pelty)%deriv,elmar(pelty)%heslo,elcod(:,1:pnode), &
                 gpvol(1:pgaus),aux_gpcar(:,:,1:pgaus),&
                 aux_gphes(:,:,1:pgaus),ielem)


            if( porde == 1 ) then
               !
               ! Element assembly
               !
               do igaus = 1, pgaus
                  fact_densi = gpvol(igaus) * spvol_gp_CMC_chm(ielem) % a(igaus,imixf,1)
                  fact_visco = gpvol(igaus) * visco_gp_CMC_chm(ielem) % a(igaus,imixf,1)
                  fact_sphec = gpvol(igaus) * sphec_gp_CMC_chm(ielem) % a(igaus,imixf,1)
                  fact_condu = gpvol(igaus) * condu_gp_CMC_chm(ielem) % a(igaus,imixf,1)
                  do inode = 1,pnode
                     eldensi(inode) = eldensi(inode) + &
                          elmar(pelty)%shape(inode,igaus) * fact_densi
                     elvisco(inode) = elvisco(inode) + &
                          elmar(pelty)%shape(inode,igaus) * fact_visco
                     elsphec(inode) = elsphec(inode) + &
                          elmar(pelty)%shape(inode,igaus) * fact_sphec
                     elcondu(inode) = elcondu(inode) + &
                          elmar(pelty)%shape(inode,igaus) * fact_condu
                  end do

               end do

               !
               ! Nodal projection
               !
               lnods_loc(1:pnode) = lnods(1:pnode,ielem)
               do inode = 1,pnode
                  ipoin = lnods_loc(inode)
                  densi(ipoin,1) = densi(ipoin,1) + eldensi(inode)
                  visco(ipoin,1) = visco(ipoin,1) + elvisco(inode)
                  condk(ipoin,1) = condk(ipoin,1) + elcondu(inode)
                  sphek(ipoin,1) = sphek(ipoin,1) + elsphec(inode)
               end do
            else
               call runend('CHEMIC OPERATIONS MOD. CMC: projections for density and viscosity not implemented for orders higher&
                   & than 1')
            end if

        end if
     end do elements

     call solver_lumped_mass_system(1_ip,densi)
     call solver_lumped_mass_system(1_ip,visco)
     call solver_lumped_mass_system(1_ip,condk)
     call solver_lumped_mass_system(1_ip,sphek)

  end if

  end subroutine chm_nodal_projection_fundamental_CMC


  subroutine chm_soot_transfer_source_terms_CMC(imixf_rk)
  !------------------------------------------------------------------------
  ! NAME
  !****f* Chemic/mod_chm_operations_CMC/chm_soot_transfer_source_terms_CMC
  ! DESCRIPTION
  !    Project source terms for species affected by soot and soot into nodes.
  ! USES
  ! USED BY
  !    mod_chm_rk_explicit
  !***
  !------------------------------------------------------------------------

  use mod_solver,                        only :  solver_lumped_mass_system
  use def_domain,                        only :  nelem,ltype,nnode,ngaus,lorde,elmar,lnods, &
                                                 ntens
  use mod_chm_sectional_soot_model_fast, only :  Qtot_gp_chm, SSTtot_gp_chm
  use def_chemic,                        only :  src_Yk_CMC_chm, aux_val_CMC_chm, densi_int_CMC_chm

  implicit none
  integer(ip), intent(in)                :: imixf_rk
  integer(ip)                            :: ielem,inode,ipoin,igaus, ivar, iclas
  integer(ip)                            :: pnode, pgaus, porde, aux_plapl,pelty
  integer(ip)                            :: lnods_loc(mnode)
  real(rp)                               :: elcod(ndime,mnode)
  real(rp)                               :: aux_gpcar(ndime,mnode,mgaus)    ! dNk/dxj
  real(rp)                               :: aux_gphes(ntens,mnode,mgaus)    ! dNk/dxidxj
  real(rp)                               :: gpvol(mgaus)
  real(rp)                               :: fact_var
  real(rp), allocatable                  :: elvar(:,:)

  external                               :: elmcar
  external                               :: runend

  ! Project source terms from Gauss points into nodes

  if (INOTMASTER) then

     allocate(elvar(mnode,nclas_chm))

     !
     ! Loop over elements
     !
     elements: do ielem = 1, nelem
        !
        ! Element dimensions
        !
        pelty = ltype(ielem)

        if( pelty > 0 ) then
            pnode = nnode(pelty)
            pgaus = ngaus(pelty)
            porde = lorde(pelty)

            !
            ! Initialization
            !
            elcod     = 0.0_rp
            aux_plapl = 0_ip
            gpvol     = 0.0_rp
            aux_gpcar = 0.0_rp
            aux_gphes = 0.0_rp
            elvar     = 0.0_rp

            !
            ! Gather coordinates at the element
            !
            call chm_elmgac_coords(pnode, lnods(1:pnode,ielem), elcod(:,1:pnode))

            !
            ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
            !
            call elmcar(&
                 pnode,pgaus,aux_plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
                 elmar(pelty)%deriv,elmar(pelty)%heslo,elcod(:,1:pnode), &
                 gpvol(1:pgaus),aux_gpcar(:,:,1:pgaus),&
                 aux_gphes(:,:,1:pgaus),ielem)


            if( porde == 1 ) then
               !
               ! Element assembly
               !
               do igaus = 1, pgaus
                  do ivar = 1, nspec_chm
                     fact_var = gpvol(igaus) * SSTtot_gp_chm(ielem) % a(igaus,ivar,1)
                     do inode = 1,pnode
                        elvar(inode,ivar) = elvar(inode,ivar) + &
                           elmar(pelty)%shape(inode,igaus) * fact_var
                     end do
                  end do
                  do ivar = 1, nclas_chm - nspec_chm
                     fact_var = gpvol(igaus) * Qtot_gp_chm(ielem) % a(igaus,ivar,1)
                     do inode = 1,pnode
                        elvar(inode,ivar+nspec_chm) = elvar(inode,ivar+nspec_chm) + &
                           elmar(pelty)%shape(inode,igaus) * fact_var
                     end do
                  end do
               end do

               !
               ! Nodal projection
               !
               lnods_loc(1:pnode) = lnods(1:pnode,ielem)
               do inode = 1,pnode
                  ipoin = lnods_loc(inode)
                  do iclas = 1, nclas_chm
                     src_Yk_CMC_chm(imixf_rk,ipoin,iclas) = src_Yk_CMC_chm(imixf_rk,ipoin,iclas) + &
                        elvar(inode,iclas)
                  end do
               end do
            else
               call runend('CHEMIC OPERATIONS MOD. CMC: projections for density and viscosity not implemented for orders higher&
                   & than 1')
            end if

        end if
     end do elements

     do iclas = 1, nclas_chm
        ! Chemical source terms for species
        do ipoin = 1, npoin
           aux_val_CMC_chm(ipoin) = src_Yk_CMC_chm(imixf_rk,ipoin,iclas)
        end do
        call solver_lumped_mass_system(1_ip,aux_val_CMC_chm)
        do ipoin = 1, npoin
           src_Yk_CMC_chm(imixf_rk,ipoin,iclas) = aux_val_CMC_chm(ipoin)
        end do
     end do
     deallocate(elvar)

     ! Remove density

     !!!!! VERIFICAR SI ESTO ES NECESARIO PARA LAS ESPECIES
     do iclas = 1, nclas_chm
        do ipoin = 1, npoin
           src_Yk_CMC_chm(imixf_rk,ipoin,iclas) = src_Yk_CMC_chm(imixf_rk,ipoin,iclas) / &
              densi_int_CMC_chm(ipoin)
        end do
     end do
  end if

  end subroutine chm_soot_transfer_source_terms_CMC


  subroutine chm_smooth_field_CMC
  !------------------------------------------------------------------------
  ! NAME
  !****f* Chemic/mod_chm_operations_CMC/chm_smooth_field_CMC
  ! DESCRIPTION
  !    Field smoothing by averaging with Gauss points.
  ! USES
  ! USED BY
  !    chm_outvar
  !***
  !------------------------------------------------------------------------

  use mod_solver,             only :  solver_lumped_mass_system
  use def_domain,             only :  nelem,ltype,nnode,ngaus,lorde,elmar,lnods, &
                                      ntens
  use def_chemic,             only :  matr_scal_CMC_chm, aux_val_CMC_chm

  implicit none
  integer(ip)                      :: ielem,inode,ipoin,igaus
  integer(ip)                      :: pnode, pgaus, porde, aux_plapl,pelty
  real(rp)                         :: elcod(ndime,mnode)
  real(rp)                         :: aux_gpcar(ndime,mnode,mgaus)    ! dNk/dxj
  real(rp)                         :: aux_gphes(ntens,mnode,mgaus)    ! dNk/dxidxj
  real(rp)                         :: gpvol(mgaus)
  real(rp)                         :: gpaux(mgaus)
  real(rp)                         :: fact_aux
  real(rp)                         :: elaux(mnode)
  real(rp)                         :: aux_smooth_nodes(mnode)

  external                         :: elmcar
  external                         :: runend


  if (INOTMASTER) then

     aux_val_CMC_chm = 0.0_rp

     !
     ! Loop over elements
     !
     elements: do ielem = 1, nelem
        !
        ! Element dimensions
        !
        pelty = ltype(ielem)

        if( pelty > 0 ) then
            pnode = nnode(pelty)
            pgaus = ngaus(pelty)
            porde = lorde(pelty)

            !
            ! Initialization
            !
            elcod            = 0.0_rp
            aux_plapl        = 0_ip
            gpvol            = 0.0_rp
            gpaux            = 0.0_rp
            aux_gpcar        = 0.0_rp
            aux_gphes        = 0.0_rp
            elaux            = 0.0_rp
            aux_smooth_nodes = 0.0_rp

            !
            ! Gather coordinates at the element
            !
            call chm_elmgac_coords(pnode, lnods(1:pnode,ielem), elcod(:,1:pnode))

            !
            ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
            !
            call elmcar(&
                 pnode,pgaus,aux_plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
                 elmar(pelty)%deriv,elmar(pelty)%heslo,elcod(:,1:pnode), &
                 gpvol(1:pgaus),aux_gpcar(:,:,1:pgaus),&
                 aux_gphes(:,:,1:pgaus),ielem)

            ! Nodal values
            do inode = 1,pnode
               ipoin        = lnods(inode,ielem)
               elaux(inode) = matr_scal_CMC_chm(ipoin)
            end do

            ! Values at Gauss points
            do igaus = 1,pgaus
               do inode = 1,pnode
                  gpaux(igaus)  = gpaux(igaus) + elmar(pelty)%shape(inode,igaus) * elaux(inode)
               end do
            end do

            if( porde == 1 ) then
               !
               ! Element assembly
               !
               do igaus = 1, pgaus
                  fact_aux = gpvol(igaus) * gpaux(igaus)
                  do inode = 1,pnode
                     aux_smooth_nodes(inode) = aux_smooth_nodes(inode) + &
                          elmar(pelty)%shape(inode,igaus) * fact_aux
                  end do
               end do

               !
               ! Nodal projection
               !
               do inode = 1,pnode
                  ipoin = lnods(inode,ielem)
                  aux_val_CMC_chm(ipoin) = aux_val_CMC_chm(ipoin) + aux_smooth_nodes(inode)
               end do
            else
               call runend('CHEMIC OPERATIONS MOD. CMC: projections for density and viscosity not implemented for orders higher&
                & than 1')
            end if

        end if
     end do elements

     call solver_lumped_mass_system(1_ip,aux_val_CMC_chm)

  end if

  end subroutine chm_smooth_field_CMC


  subroutine chm_mxt_fr_integr_previous_steps_CMC(nfield, fields_Z, &
                            Zavg, Zvar, int_fields)

     !-----------------------------------------------------------------------
     !****f* Chemic/chm_mxt_fr_integr_previous_steps_CMC
     ! NAME
     !    chm_mxt_fr_integr_previous_steps_CMC
     ! DESCRIPTION
     !    This routine takes the fields in mixture fraction
     !    for a given domain, does some preliminary actions and performs the
     !    corresponding integrals for each field.
     !
     !    - nfield: number of fields to be integrated.
     !    - fields_Z(nZ_CMC_chm,nfield): matrix containing all the fields for each
     !    mixture fraction.
     !    - Zavg: field of averaged/filtered mixture fraction.
     !    - Zvar: field of averaged/filtered mixture fraction variance.
     !    - int_fields(nfield): integrated fields.
     ! USED BY
     !    chm_integrate_flow_var_points_CMC
     !***
     !-----------------------------------------------------------------------

     use def_chemic,          only :  nZ_CMC_chm, S_threshold, Zs_CMC_chm, &
                                      extremes_Z_CMC_chm

     implicit none
     integer(ip), intent(in)       :: nfield
     real(rp),    intent(in)       :: fields_Z(nZ_CMC_chm, nfield), Zavg, Zvar
     real(rp),    intent(out)      :: int_fields(nfield)

     real(rp)                      :: PDF_param(7)
     real(rp)                      :: S

     PDF_param(4) = Zs_CMC_chm

     if (Zavg <= extremes_Z_CMC_chm(1) .or. Zavg >= extremes_Z_CMC_chm(2)) then
        ! Direct interpolation
        call chm_mixture_fraction_interpolation_CMC(fields_Z(1:nZ_CMC_chm, 1:nfield), &
                nfield, Zavg, int_fields(1:nfield))

     else
        S = Zvar / (Zavg*(Zs_CMC_chm - Zavg))
        if (S <= S_threshold) then
           ! Direct interpolation
           call chm_mixture_fraction_interpolation_CMC(fields_Z(1:nZ_CMC_chm, 1:nfield), &
                   nfield, Zavg, int_fields(1:nfield))
        else
           ! Perform integration
           call chm_find_PDF_parameters_CMC(Zavg, Zvar, PDF_param)
           call chm_mixtFraction_integration_CMC(fields_Z(1:nZ_CMC_chm, 1:nfield), &
                   nfield, PDF_param, int_fields(1:nfield))
        end if
     end if

  end subroutine chm_mxt_fr_integr_previous_steps_CMC



  subroutine chm_mixtFraction_integration_CMC(fields_Z, nfield, PDF_param, int_fields)

     !-----------------------------------------------------------------------
     !****f* Chemic/chm_mixtFraction_integration_CMC
     ! NAME
     !    chm_mixtFraction_integration_CMC
     ! DESCRIPTION
     !    This routine computes the integrals for the PDF/FPDF times each
     !    variable profile along mixture fraction. It only computes such
     !    integral for one pair of averaged/filtered mixture fraction and variance
     !    but for all the fields.
     !
     !    - fields_Z(nmixt,nfield): fields along mixture fraction (point fixed
     !    in physical space).
     !    - nfield: number of fields to be integrated.
     !    - Z_CMC_chm: vector of mixture fractions.
     !    - nZ_CMC_chm: number of mixture fractions.
     !    - PDF_param: defining parameters for the PDF/FPDF.
     !    - p_cut: vector with the positions of Z_CMC_chm where divide the integral.
     !    - int_fields(nfield): integrated fields.
     !
     !    Actions:
     !    1. Evaluation of the regularized incomplete beta function for coefficients
     !    pairs (alfa, beta) and (alfa+1, beta).
     !    2. Compute the integral from the field values and the incomplete beta
     !    functions considering the constant and linear contributions from the fields
     !    (piecewise linear interpolation).
     !    The algorithm is not iterative.
     ! USED BY
     !    chm_mxt_fr_integr_previous_steps_CMC
     !***
     !-----------------------------------------------------------------------

     use def_chemic,          only :  nZ_CMC_chm, Z_CMC_chm

     implicit none
     integer(ip), intent(in)       :: nfield
     real(rp),    intent(in)       :: fields_Z(nZ_CMC_chm,nfield), PDF_param(7)
     real(rp),    intent(out)      :: int_fields(nfield)

     integer(ip)                   :: imixf, ifield
     real(rp)                      :: incompl_beta(nZ_CMC_chm,2), ratio_f, diffZ_inv
     real(rp)                      :: factor_0, factor_1

     int_fields(1:nfield) = 0.0_rp

     incompl_beta(1,1:2) = 0.0_rp       ! Since Z_CMC_chm(1) = 0
     incompl_beta(nZ_CMC_chm,1:2) = 1.0_rp  ! This is B(x,y)/B(x,y)=1
     do imixf=2, nZ_CMC_chm-1
        call incob(PDF_param(1), PDF_param(2), Z_CMC_chm(imixf)/PDF_param(4), PDF_param(6),incompl_beta(imixf,1))
        call incob(PDF_param(1)+1.0_rp, PDF_param(2), Z_CMC_chm(imixf)/PDF_param(4), PDF_param(7), incompl_beta(imixf,2))
     end do

     do imixf = 1, nZ_CMC_chm-1
        diffZ_inv = 1.0_rp / (Z_CMC_chm(imixf+1) - Z_CMC_chm(imixf))
        factor_0 = incompl_beta(imixf+1,1) - incompl_beta(imixf,1)
        factor_1 = PDF_param(5) * (incompl_beta(imixf+1,2) - incompl_beta(imixf,2))

        do ifield = 1, nfield
           ratio_f = (fields_Z(imixf+1,ifield) - fields_Z(imixf,ifield)) * diffZ_inv

           ! Contribution from the constant part
           int_fields(ifield) = int_fields(ifield) + (fields_Z(imixf,ifield) - ratio_f*Z_CMC_chm(imixf)) * factor_0

           ! Contribution from the linear part
           int_fields(ifield) = int_fields(ifield) + ratio_f * factor_1
        end do
     end do


  end subroutine chm_mixtFraction_integration_CMC



  subroutine chm_mixture_fraction_interpolation_CMC(fields_Z, nfield, Zavg, int_fields)

     !-----------------------------------------------------------------------
     !****f* Chemic/chm_mixture_fraction_interpolation_CMC
     ! NAME
     !    chm_mixture_fraction_interpolation_CMC
     ! DESCRIPTION
     !    This routine interpolates the fields at a given mixture fraction value.
     !
     !    - fields_Z(nmixt,nfield): fields along mixture fraction (point fixed
     !    in physical space).
     !    - nfield: number of fields to be integrated.
     !    - Z_CMC_chm: vector of mixture fractions.
     !    - nZ_CMC_chm: number of mixture fractions.
     !    - Zavg: averaged/filtered mixture fraction.
     !    - int_fields(nfield): interpolated fields.
     ! USED BY
     !    chm_mxt_fr_integr_previous_steps_CMC
     !***
     !-----------------------------------------------------------------------

     use def_chemic,          only :  nZ_CMC_chm, Z_CMC_chm

     implicit none
     integer(ip), intent(in)       :: nfield
     real(rp),    intent(in)       :: fields_Z(nZ_CMC_chm,nfield), Zavg
     real(rp),    intent(out)      :: int_fields(nfield)

     integer(ip)                   :: pos_aux1, pos_aux2, ifield
     real(rp)                      :: fact, diff_Z_int, diff_Z


     call find_pos_min_diff(Z_CMC_chm, nZ_CMC_chm, Zavg, pos_aux1)
     if (Zavg == Z_CMC_chm(pos_aux1)) then
        int_fields(1:nfield) = fields_Z(pos_aux1, 1:nfield)
     else
        if(Zavg >= Z_CMC_chm(pos_aux1)) then
           pos_aux2 = pos_aux1 + 1_ip
        else
           pos_aux2 = pos_aux1
           pos_aux1 = pos_aux1 - 1_ip
        end if

        if (pos_aux2 == 1) then  ! It may happen if Zavg=-1e-30
           pos_aux1 = 1
           pos_aux2 = 2
        else if (pos_aux1 == nZ_CMC_chm) then ! It may happen if Z_CMC_chm(pos_aux1) = Zs
           pos_aux1 = nZ_CMC_chm - 1
           pos_aux2 = nZ_CMC_chm
        end if

        diff_Z_int = Z_CMC_chm(pos_aux2) - Z_CMC_chm(pos_aux1)
        diff_Z = Zavg - Z_CMC_chm(pos_aux1)
        do ifield = 1, nfield
           fact = (fields_Z(pos_aux2, ifield) - fields_Z(pos_aux1, ifield)) / diff_Z_int
           int_fields(ifield) = fields_Z(pos_aux1, ifield) + diff_Z * fact
        end do
     end if

  end subroutine chm_mixture_fraction_interpolation_CMC



  subroutine chm_find_PDF_parameters_CMC(Zavg, Zvar_avg, PDF_param)

     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/chm_find_PDF_parameters_CMC
     ! NAME
     !    chm_find_PDF_parameters_CMC
     ! DESCRIPTION
     !    It computes the PDF parameters for a beta function: alpha, beta,
     !    Zs^(1-alfa-beta)/B(alfa,beta), Zs, Zs*B(alfa+1,beta)/B(alfa,beta),
     !    B(alfa,beta), B(alfa+1,beta)
     ! USED BY
     !    chm_mxt_fr_integr_previous_steps_CMC, chm_elmpre_pdf_CMC
     !***
     !-----------------------------------------------------------------------

     implicit none

     real(rp), parameter           :: fact_Zv = 0.95_rp

     real(rp), intent(in)          :: Zavg, Zvar_avg
     real(rp), intent(inout)       :: PDF_param(7)
     real(rp)                      :: Zvar_aux

     Zvar_aux = min(Zvar_avg, fact_Zv*Zavg*(PDF_param(4)-Zavg))

     PDF_param(1) = ( Zavg*(PDF_param(4)-Zavg)/Zvar_aux - 1.0_rp ) * Zavg / PDF_param(4)
     PDF_param(2) = PDF_param(1) * (PDF_param(4)/Zavg - 1.0_rp)
     call beta_func(PDF_param(1), PDF_param(2), PDF_param(6))
     PDF_param(3) = PDF_param(4)**(1.0_rp - PDF_param(1) - PDF_param(2)) / PDF_param(6)
     call beta_func(PDF_param(1)+1.0_rp, PDF_param(2), PDF_param(7))
     PDF_param(5) = PDF_param(4) * PDF_param(7) / PDF_param(6)

  end subroutine chm_find_PDF_parameters_CMC



  subroutine find_pos_min_diff(vector, N, val, pos)

     !-----------------------------------------------------------------------
     !****f* Chemic/find_pos_min_diff
     ! NAME
     !    find_pos_min_diff
     ! DESCRIPTION
     !    Find the integer pos such that minimizes abs(vector(i)-val) i=1,...,N.
     !    It there are several indexes for which the difference is equal and
     !    minimum it takes the lowest index.
     ! USED BY
     !    chm_mxt_fr_integr_previous_steps_CMC, chm_mixture_fraction_interpolation_CMC
     !***
     !-----------------------------------------------------------------------

     implicit none
     integer(ip), intent(in)       :: N
     real(rp),    intent(in)       :: vector(N), val
     integer(ip), intent(out)      :: pos

     integer(ip)                   :: i
     real(rp)                      :: v_aux(N)

     v_aux(1:N) = abs(vector(1:N) - val)
     pos = 1_ip
     do i = 2, N
        if ( v_aux(i) < v_aux(pos) )  pos = i
     end do

  end subroutine find_pos_min_diff


  subroutine chm_PDF_calc_CMC(Zavg, Zvar, imixf, PDF_val)
     !-----------------------------------------------------------------------
     !****f* Chemic/chm_PDF_calc_CMC
     ! NAME
     !    chm_PDF_calc_CMC
     ! DESCRIPTION
     !    Find PDF value for given Zaverage and Zvariance.
     ! USED BY
     !    chm_integrate_flow_var_points_CMC
     !***
     !-----------------------------------------------------------------------
     use def_chemic,       only :  extremes_Z_CMC_chm, PDF_min_CMC_chm, &
                                   S_threshold, Z_CMC_chm, Zs_CMC_chm

     implicit none
     real(rp), intent(in)       :: Zavg
     real(rp), intent(in)       :: Zvar
     integer(ip), intent(in)    :: imixf
     real(rp), intent(out)      :: PDF_val
     real(rp)                   :: S, aux_Zvar, PDF_param(7)

     PDF_param(4) = Zs_CMC_chm

     if (Zavg <= extremes_Z_CMC_chm(1) .or. Zavg >= extremes_Z_CMC_chm(2)) then
        ! We have a Dirac's delta at one extreme: we assume a PDF that is zero
        ! everywhere except from the extreme to the next/previous level of
        ! mixture fraction. As the extremes are not solved PDF=0 for the rest
        ! of mixture fraction levels
        PDF_val = PDF_min_CMC_chm(imixf)

     else
        S = Zvar / (Zavg*(Zs_CMC_chm-Zavg))

        if (S < S_threshold) then
           S = S_threshold
           aux_Zvar = S * Zavg*(Zs_CMC_chm-Zavg)
        else
           aux_Zvar = Zvar
        end if

        ! Find pdf value
        call chm_find_PDF_parameters_CMC(Zavg, aux_Zvar, PDF_param(1:7))
        PDF_val = PDF_param(3) * Z_CMC_chm(imixf)**(PDF_param(1)-1.0_rp) * &
                   (PDF_param(4) - Z_CMC_chm(imixf))**(PDF_param(2)-1.0_rp)

        PDF_val = max(PDF_val, PDF_min_CMC_chm(imixf))
     end if

  end subroutine chm_PDF_calc_CMC


  subroutine min_PDF_value_CMC
     !-----------------------------------------------------------------------
     !****f* Chemic/min_PDF_value_CMC
     ! NAME
     !    min_PDF_value_CMC
     ! DESCRIPTION
     !    Find minima values for the PDF at each mixture fraction to avoid
     !    jumps in physical space for the PDF. The maximum of all possible
     !    values is taken for Zaverage equal to the extremes found close to
     !    Z=0 and Z=Zs with minimum segregation factor (it provides the maximum
     !    PDF value for fixed Zaverage and Z).
     ! USED BY
     !    chm_reaphy
     !***
     !-----------------------------------------------------------------------

     use def_chemic,       only : Zs_CMC_chm, S_threshold, PDF_min_CMC_chm, &
                                   extremes_Z_CMC_chm, nZ_CMC_chm, Z_CMC_chm
     implicit none

     real(rp), parameter        :: small_PDF = 0.1_rp

     integer(ip)                :: iZ
     real(rp)                   :: Zvar_Z0, Zvar_Zs, PDF_param_Z0(7), &
                                   PDF_param_Zs(7),  PDF_Z0, PDF_Zs

     PDF_param_Z0(4) = Zs_CMC_chm
     PDF_param_Zs(4) = Zs_CMC_chm

     Zvar_Z0 = S_threshold * extremes_Z_CMC_chm(1)  * (Zs_CMC_chm - extremes_Z_CMC_chm(1))
     call chm_find_PDF_parameters_CMC(extremes_Z_CMC_chm(1), Zvar_Z0, PDF_param_Z0)

     Zvar_Zs = S_threshold * extremes_Z_CMC_chm(2)  * (Zs_CMC_chm - extremes_Z_CMC_chm(2))
     call chm_find_PDF_parameters_CMC(extremes_Z_CMC_chm(2), Zvar_Zs, PDF_param_Zs)

     do iZ = 2, nZ_CMC_chm-1
        PDF_Z0 = PDF_param_Z0(3) * Z_CMC_chm(iZ) ** (PDF_param_Z0(1)-1.0_rp) &
                     * (Zs_CMC_chm - Z_CMC_chm(iZ)) ** (PDF_param_Z0(2)-1.0_rp)
        PDF_Zs = PDF_param_Zs(3) * Z_CMC_chm(iZ) ** (PDF_param_Zs(1)-1.0_rp) &
                     * (Zs_CMC_chm - Z_CMC_chm(iZ)) ** (PDF_param_Zs(2)-1.0_rp)
        PDF_min_CMC_chm(iZ) = max(PDF_Z0, PDF_Zs)
        PDF_min_CMC_chm(iZ) = max(PDF_min_CMC_chm(iZ), small_PDF)
     end do

  end subroutine min_PDF_value_CMC


  subroutine chm_average_field_mesh_CMC(n_field,fields,fields_int,field_Z,field_Zv,imixf)
    !-----------------------------------------------------------------------
    !****f* chemic/chm_average_field_mesh_CMC
    ! NAME
    !    chm_average_field_mesh_CMC
    ! DESCRIPTION
    !    Compute the volumetric integrals for passing information from CFD to
    !    CMC mesh in case meshes are different. Integrals can be done in 2 ways:
    !    a) Weighed by the PDF: in this case there is need of fields for mixture
    !       fraction and its variance and field_Z and field_Zv are compulsory.
    !    b) Not weighed by the PDF: do not provide neither field_Z nor field_Zv.
    ! USES
    ! USED BY
    !
    !-----------------------------------------------------------------------

    use def_domain,               only :  lnods, ltype, nnode, mgaus_cloud_gp, &
                                          gp_total_cloud_gp, ngaus_cloud_gp, &
                                          elmar_cloud_gp, nelem
    use def_chemic,               only :  normal_CMC_mesh_CMC_chm

    implicit none
    integer(ip), intent(in)            :: n_field   ! Number of fields
    real(rp),    intent(in)            :: fields(gp_total_cloud_gp,n_field)  ! Fields to be integrated
    real(rp),    intent(out)           :: fields_int(npoin,n_field)  ! Integrated fields
    real(rp),    intent(in), optional  :: field_Z(gp_total_cloud_gp)   ! Mixture fraction field
    real(rp),    intent(in), optional  :: field_Zv(gp_total_cloud_gp)  ! Mixture fraction variance field
    integer(ip), intent(in), optional  :: imixf  ! Level of mixture fraction

    integer(ip)                        :: ielem, igaus, inode, ipoin, ifield
    integer(ip)                        :: pelty, pnode, pgaus, count_gaus, cond_fields
    integer(ip)                        :: ipoin_elem(mnode)
    real(rp)                           :: fields_elem_g(mgaus_cloud_gp,n_field)
    real(rp)                           :: normal_elem_g(mgaus_cloud_gp)
    real(rp)                           :: integr_field_elem(n_field), elcod(ndime,mnode)
    real(rp)                           :: xjacm(ndime,ndime)
    real(rp)                           :: gpdet, integr_normal_elem, aux_vol

    external                           :: jacdet

    if (INOTMASTER) then
       fields_int              = 0.0_rp
       normal_CMC_mesh_CMC_chm = 0.0_rp
       count_gaus              = 1_ip
       cond_fields             = 0_ip
       normal_elem_g           = 1.0_rp

       if (present(field_Zv)) then
          cond_fields = 1_ip
       end if

       elements: do ielem = 1, nelem
          pelty = ltype(ielem)
          if (pelty > 0) then
             integr_field_elem  = 0.0_rp
             integr_normal_elem = 0.0_rp

             pnode = nnode(pelty)
             pgaus = ngaus_cloud_gp(pelty)

             do inode = 1,pnode
                ipoin_elem(inode)    = lnods(inode,ielem)
                elcod(1:ndime,inode) = coord(1:ndime,ipoin_elem(inode))
             end do

             ! Find integrand and transfer to local matrix
             do ifield = 1, n_field
                fields_elem_g(1:pgaus,ifield) = fields(count_gaus:count_gaus+pgaus-1_ip,ifield)
             end do

             if (cond_fields == 1_ip) then
                ! Compute PDF value
                normal_elem_g = 1.0_rp
                do igaus = 1, pgaus
                   call chm_PDF_calc_CMC(field_Z(count_gaus+igaus-1_ip), field_Zv(count_gaus+igaus-1_ip), &
                           imixf, normal_elem_g(igaus))
                end do
             end if

             do igaus = 1, pgaus
                ! Find jacobian
                call jacdet(ndime, pnode, elcod, elmar_cloud_gp(pelty)%deriv(1,1,igaus), &
                        xjacm,gpdet)
                aux_vol = elmar_cloud_gp(pelty)%weigp(igaus) * gpdet

                integr_field_elem  = integr_field_elem  + fields_elem_g(igaus,:) * normal_elem_g(igaus) * aux_vol
                integr_normal_elem = integr_normal_elem + normal_elem_g(igaus)   * aux_vol
             end do

             do inode = 1,pnode
                fields_int(ipoin_elem(inode),:) = fields_int(ipoin_elem(inode),:) + &
                   integr_field_elem
                normal_CMC_mesh_CMC_chm(ipoin_elem(inode)) = &
                   normal_CMC_mesh_CMC_chm(ipoin_elem(inode)) + integr_normal_elem
             end do

             count_gaus = count_gaus + pgaus
          end if
       end do elements

       do ipoin = 1, npoin
          fields_int(ipoin,:) = fields_int(ipoin,:) / normal_CMC_mesh_CMC_chm(ipoin)
       end do
    end if

  end subroutine chm_average_field_mesh_CMC



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! C H E M I C A L   C A L C U L A T I O N S !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine chm_integrate_chem_source_CMC(dt)

     !-----------------------------------------------------------------------
     !****f* Chemic/chm_integrate_chem_source_CMC
     ! NAME
     !    chm_integrate_chem_source_CMC
     ! DESCRIPTION
     !    Compute chemical source terms for CMC transport equations.
     ! USED BY
     !
     !***
     !-----------------------------------------------------------------------

     use def_chemic,             only :  src_Yk_CMC_chm, Z_CMC_chm, chem_int_iZ_CMC_chm, &
                                         nZ_no_chm_int_CMC_chm, nZ_chm_int_CMC_chm
     use mod_physics,            only :  physics_T_2_HCp
#ifdef CANTERA
     use def_master,             only :  prthe
     use def_chemic,             only :  Yk_CMC_chm, t_chem_integ_CMC_chm, temp_CMC_chm
#endif

     implicit none
     real(rp), intent(in)          :: dt
     integer(ip)                   :: ii, jj, nZ_inf, nZ_sup
     real(rp)                      :: inv_dif_Z, fact_Z
#ifdef CANTERA
     real(rp), parameter           :: temp_min = 500.0_rp
     integer(ip)                   :: imixf, ipoin, imixf_chem
     real(rp)                      :: dt_inv
     real(rp)                      :: T_next, Yk_next(nspec_chm)
     real(rp)                      :: t_start, t_end

     external                      :: cantera_integrate
#endif
     external                      :: runend

#ifndef CANTERA
    ! Early exit until alternatives are added
    call runend('chm_integrate_chem_source_CMC: Please, build with support for cantera')
    if (.false.) print *, dt ! avoid compiler warning (unused-dummy-argument)
#endif

     if( INOTMASTER ) then

#ifdef CANTERA
        t_chem_integ_CMC_chm = 0.0_rp
        dt_inv = 1.0_rp / dt
        do imixf = 1, nZ_chm_int_CMC_chm
           imixf_chem = chem_int_iZ_CMC_chm(imixf)
           do ipoin = 1, npoin

              if (temp_CMC_chm(imixf_chem,ipoin) > temp_min)  then
                 T_next  = temp_CMC_chm(imixf_chem,ipoin)
                 Yk_next = Yk_CMC_chm(imixf_chem,ipoin,1:nspec_chm)
                 !!! ¿DESCOMENTAR?
                 !!!alpha = 0.0_rp
                 !!!do kk = 1, nspec_chm
                 !!!   alpha = alpha + Yk_next(kk)
                 !!!end do
                 !!!Yk_next_chm = Yk_next / alpha
                 call cpu_time(t_start)
                 !!! ¿DESCOMENTAR?
                 !!!call cantera_integrate(T_next, prthe(1), Yk_next_chm, dt)
                 call cantera_integrate(T_next, prthe(1), Yk_next, dt)
                 call cpu_time(t_end)
                 !!! ¿DESCOMENTAR?
                 !!!Yk_next = Yk_next_chm * alpha

                 ! Compute chemical source terms and update mass fractions (note: consider
                 ! contributions from soot formation/oxidation for species chemical source terms)
                 src_Yk_CMC_chm(imixf_chem,ipoin,1:nspec_chm) = &
                    src_Yk_CMC_chm(imixf_chem,ipoin,1:nspec_chm) + &
                    (Yk_next(1:nspec_chm) - Yk_CMC_chm(imixf_chem,ipoin,1:nspec_chm)) * dt_inv
                 temp_CMC_chm(imixf_chem,ipoin)           = T_next
                 Yk_CMC_chm(imixf_chem,ipoin,1:nspec_chm) = Yk_next(1:nspec_chm)
                 t_chem_integ_CMC_chm(imixf_chem,ipoin)   = t_end - t_start

              end if

           end do

        end do
#endif

        ! Interpolations if required
        if (nZ_no_chm_int_CMC_chm > 0_ip) then
           do ii = 1, nZ_chm_int_CMC_chm - 1_ip
              nZ_inf = chem_int_iZ_CMC_chm(ii)
              nZ_sup = chem_int_iZ_CMC_chm(ii+1_ip)
              inv_dif_Z = 1.0_rp / (Z_CMC_chm(nZ_sup) - Z_CMC_chm(nZ_inf))
              do jj = nZ_inf + 1_ip, nZ_sup - 1_ip
                 fact_Z = (Z_CMC_chm(jj) - Z_CMC_chm(nZ_inf)) * inv_dif_Z
                 src_Yk_CMC_chm(jj,:,:) = src_Yk_CMC_chm(nZ_inf,:,:) * (1.0_rp - fact_Z) + &
                    src_Yk_CMC_chm(nZ_sup,:,:) * fact_Z
              end do
           end do
        end if
     end if

  end subroutine chm_integrate_chem_source_CMC


  subroutine chm_construct_initialSol_from_chemistry_CMC
     !-----------------------------------------------------------------------
     !****f* Chemic/chm_construct_initialSol_from_chemistry_CMC
     ! NAME
     !    chm_construct_initialSol_from_chemistry_CMC
     ! DESCRIPTION
     !    Compute homogeneous reactors solutions and save in a table
     ! USED BY
     !
     !***
     !-----------------------------------------------------------------------

     !!!use def_chemic,               only :  alpha_grid_CMC_chm, n_alpha_grid_CMC_chm, &
     !!!                                      posZ_chem_integr_CMC_chm, Z_CMC_chm, &
     !!!                                      nZ_CMC_chm, homog_reactors_CMC_chm, &
     !!!                                      rscal_inert_CMC_chm, rscal_equil_CMC_chm, &
     !!!                                      temp_inert_CMC_chm, temp_equil_CMC_chm, &
     !!!                                      Zs_CMC_chm

     !!!implicit none
     !!!
     !!!integer(ip)                        :: iZ, ialfa, iclas
     !!!real(rp)                           :: Yk_start(nclas_chm), Yk_end(nclas_chm)
     !!!real(rp)                           :: T_start, T_end
     !!!real(rp)                           :: Yk_on_normGrid(n_alpha_grid_CMC_chm,nclas_chm)
     !!!real(rp)                           :: Yk_left(n_alpha_grid_CMC_chm,nclas_chm)
     !!!real(rp)                           :: Yk_right(n_alpha_grid_CMC_chm,nclas_chm)
     !!!real(rp), allocatable              :: Yk_bottom_l(:,:)
     !!!real(rp), allocatable              :: Yk_top_l(:,:)
     !!!real(rp), allocatable              :: Yk_bottom_r(:,:)
     !!!real(rp), allocatable              :: Yk_top_r(:,:)

     !!!! Boundary conditions
     !!!do iclas = 1, nclas_chm
     !!!   homog_reactors_CMC_chm(1,1:n_alpha_grid_CMC_chm,iclas) = rscal_inert_CMC_chm(1,iclas+1_ip)
     !!!   homog_reactors_CMC_chm(nZ_CMC_chm,1:n_alpha_grid_CMC_chm,iclas) = rscal_inert_CMC_chm(nZ_CMC_chm,iclas+1_ip)
     !!!end do

     !!!! Compute homogeneous reactor solutions
     !!!do iZ = posZ_chem_integr_CMC_chm(1), posZ_chem_integr_CMC_chm(2)

     !!!   Yk_start = rscal_inert_CMC_chm(iZ,2:nclas_chm+1_ip)
     !!!   T_start  = temp_inert_CMC_chm(iZ)
     !!!   Yk_end   = rscal_equil_CMC_chm(iZ,2:nclas_chm+1_ip)
     !!!   T_end    = temp_equil_CMC_chm(iZ)
     !!!   call chm_homogeneous_reactor_CMC(Yk_start, T_start, Yk_end, T_end, Yk_on_normGrid)

     !!!   ! Transfer information to homog_reactors_CMC_chm
     !!!   do ialfa = 1, n_alpha_grid_CMC_chm
     !!!      homog_reactors_CMC_chm(iZ,ialfa,1:nclas_chm) = Yk_on_normGrid(ialfa,1:nclas_chm)
     !!!   end do
     !!!end do

     !!!! Do interpolations
     !!!if (posZ_chem_integr_CMC_chm(1) >= 3_ip) then
     !!!   do ialfa = 1, n_alpha_grid_CMC_chm
     !!!      Yk_left(ialfa,1:nclas_chm)  = homog_reactors_CMC_chm(1,ialfa,1:nclas_chm)
     !!!      Yk_right(ialfa,1:nclas_chm) = homog_reactors_CMC_chm(posZ_chem_integr_CMC_chm(1),ialfa,1:nclas_chm)
     !!!   end do

     !!!   allocate(Yk_bottom_l(posZ_chem_integr_CMC_chm(1),nclas_chm))
     !!!   allocate(Yk_top_l(posZ_chem_integr_CMC_chm(1),nclas_chm))

     !!!   do iZ = 1, posZ_chem_integr_CMC_chm(1)
     !!!      Yk_bottom_l(iZ,1:nclas_chm) = rscal_inert_CMC_chm(iZ,2:nclas_chm+1_ip)
     !!!      Yk_top_l(iZ,1:nclas_chm)    = rscal_equil_CMC_chm(iZ,2:nclas_chm+1_ip)
     !!!   end do

     !!!   call chm_interpolate_HR_CMC(1_ip, posZ_chem_integr_CMC_chm(1), Z_CMC_chm(1), &
     !!!           Z_CMC_chm(posZ_chem_integr_CMC_chm(1)), &
     !!!           Yk_left, Yk_right, Yk_bottom_l, Yk_top_l, &
     !!!           homog_reactors_CMC_chm(2:posZ_chem_integr_CMC_chm(1)-1_ip,:,:))
     !!!end if

     !!!if (posZ_chem_integr_CMC_chm(2) <= nZ_CMC_chm-2_ip) then
     !!!   do ialfa = 1, n_alpha_grid_CMC_chm
     !!!      Yk_left(ialfa,1:nclas_chm)  = homog_reactors_CMC_chm(posZ_chem_integr_CMC_chm(2),ialfa,1:nclas_chm)
     !!!      Yk_right(ialfa,1:nclas_chm) = homog_reactors_CMC_chm(nZ_CMC_chm,ialfa,1:nclas_chm)
     !!!   end do

     !!!   allocate(Yk_bottom_r(nZ_CMC_chm - posZ_chem_integr_CMC_chm(2) + 1_ip,nclas_chm))
     !!!   allocate(Yk_top_r(nZ_CMC_chm - posZ_chem_integr_CMC_chm(2) + 1_ip,nclas_chm))

     !!!   do iZ = posZ_chem_integr_CMC_chm(2), nZ_CMC_chm
     !!!      Yk_bottom_r(iZ-posZ_chem_integr_CMC_chm(2)+1_ip,1:nclas_chm) = rscal_inert_CMC_chm(iZ,2:nclas_chm+1_ip)
     !!!      Yk_top_r(iZ-posZ_chem_integr_CMC_chm(2)+1_ip,1:nclas_chm)    = rscal_equil_CMC_chm(iZ,2:nclas_chm+1_ip)
     !!!   end do

     !!!   call chm_interpolate_HR_CMC(posZ_chem_integr_CMC_chm(2), &
     !!!           nZ_CMC_chm-posZ_chem_integr_CMC_chm(2)+1_ip, &
     !!!           Z_CMC_chm(posZ_chem_integr_CMC_chm(2)), Zs_CMC_chm, &
     !!!           Yk_left, Yk_right, Yk_bottom_r, Yk_top_r, &
     !!!           homog_reactors_CMC_chm(posZ_chem_integr_CMC_chm(2)+1_ip:nZ_CMC_chm-1_ip,:,:))
     !!!end if

  end subroutine chm_construct_initialSol_from_chemistry_CMC


  !!!subroutine chm_homogeneous_reactor_CMC(Yk_start, T_start, Yk_end, T_end, Yk_on_normGrid)
  !!!   !-----------------------------------------------------------------------
  !!!   !****f* Chemic/chm_homogeneous_reactor_CMC
  !!!   ! NAME
  !!!   !    chm_homogeneous_reactor_CMC
  !!!   ! DESCRIPTION
  !!!   !    Compute the evolution of a homogeneous reactor for given initial
  !!!   !    conditions and interpolate on normalized temperature grid.
  !!!   ! USED BY
  !!!   !
  !!!   !***
  !!!   !-----------------------------------------------------------------------

  !!!   use def_master,               only :  prthe
  !!!   use def_chemic,               only :  alpha_grid_CMC_chm, n_alpha_grid_CMC_chm

  !!!   implicit none
  !!!   real(rp), parameter                :: dt = 1.0e-6_rp

  !!!   real(rp), intent(in)               :: Yk_start(nclas_chm)
  !!!   real(rp), intent(in)               :: T_start
  !!!   real(rp), intent(in)               :: Yk_end(nclas_chm)
  !!!   real(rp), intent(in)               :: T_end
  !!!   real(rp), intent(out)              :: Yk_on_normGrid(n_alpha_grid_CMC_chm,nclas_chm)

  !!!   integer(ip)                        :: pos_grid
  !!!   real(rp)                           :: T_grid(n_alpha_grid_CMC_chm)
  !!!   real(rp)                           :: Yk_next(nclas_chm)
  !!!   real(rp)                           :: T_next
  !!!   real(rp)                           :: Yk_prev(nclas_chm)
  !!!   real(rp)                           :: T_prev
  !!!   real(rp)                           :: alfa

  !!!   ! Compute temperatures in which interpolate
  !!!   T_grid = T_start + alpha_grid_CMC_chm * (T_end-T_start)
  !!!
  !!!   ! Integrate chemical evolution between T_start and T_end
  !!!   Yk_prev = Yk_start
  !!!   T_prev  = T_start
  !!!   pos_grid = 2_ip
  !!!   Yk_next = Yk_prev
  !!!   T_next  = T_prev

  !!!   do while (T_prev <= T_grid(n_alpha_grid_CMC_chm-1_ip))
  !!!      do while (T_next < T_grid(pos_grid))
  !!!         Yk_prev = Yk_next
  !!!         T_prev  = T_next
  !!!         call cantera_integrate(T_next, prthe(1), Yk_next, dt)
  !!!      end do
  !!!      do while (T_next >= T_grid(pos_grid))
  !!!         ! Do linear interpolation
  !!!         alfa = (T_grid(pos_grid)-T_prev) / (T_next-T_prev)
  !!!         Yk_on_normGrid(pos_grid,1:nclas_chm) = Yk_prev + alfa * (Yk_next-Yk_prev)

  !!!         ! Update position
  !!!         pos_grid = pos_grid + 1_ip
  !!!      end do

  !!!      ! Update Yk and T
  !!!      Yk_prev = Yk_next
  !!!      T_prev  = T_next
  !!!   end do

  !!!   ! Fill first and last entries
  !!!   Yk_on_normGrid(1,1:nclas_chm)                     = Yk_start
  !!!   Yk_on_normGrid(n_alpha_grid_CMC_chm,1:nclas_chm)  = Yk_end
  !!!
  !!!end subroutine chm_homogeneous_reactor_CMC


  !!!subroutine chm_interpolate_HR_CMC(posZ_left, nZ_lr, Z_left, Z_right, &
  !!!                 Yk_left, Yk_right, Yk_bottom, Yk_top, Yk_interp)
  !!!   !-----------------------------------------------------------------------
  !!!   !****f* Chemic/chm_interpolate_HR_CMC
  !!!   ! NAME
  !!!   !    chm_interpolate_HR_CMC
  !!!   ! DESCRIPTION
  !!!   !    Compute the profiles for species by interpolating
  !!!   !
  !!!   !                     Yk_top
  !!!   !          x--x----x---x-----x---x----x
  !!!   !          |                          |
  !!!   !          |                          |
  !!!   !          x       o Interpolation    x
  !!!   !          |         points           |
  !!!   !          |                          |
  !!!   !          |                          |
  !!!   ! Yk_left  x       o                  x   Yk_right
  !!!   !          |                          |
  !!!   !          |                          |
  !!!   !          |                          |
  !!!   !          x       o                  x
  !!!   !          |                          |
  !!!   !          |                          |
  !!!   !          x--x----x---x-----x---x----x
  !!!   !                    Yk_bottom
  !!!   ! USED BY
  !!!   !
  !!!   !***
  !!!   !-----------------------------------------------------------------------

  !!!   use def_chemic,               only :  alpha_grid_CMC_chm, n_alpha_grid_CMC_chm, &
  !!!                                         Z_CMC_chm

  !!!   implicit none
  !!!   integer(ip), intent(in)            :: posZ_left
  !!!   integer(ip), intent(in)            :: nZ_lr
  !!!   real(rp),    intent(in)            :: Z_left
  !!!   real(rp),    intent(in)            :: Z_right
  !!!   real(rp),    intent(in)            :: Yk_left(n_alpha_grid_CMC_chm,nclas_chm)
  !!!   real(rp),    intent(in)            :: Yk_right(n_alpha_grid_CMC_chm,nclas_chm)
  !!!   real(rp),    intent(in)            :: Yk_bottom(nZ_lr,nclas_chm)
  !!!   real(rp),    intent(in)            :: Yk_top(nZ_lr,nclas_chm)
  !!!   real(rp),    intent(out)           :: Yk_interp(nZ_lr-2_ip,n_alpha_grid_CMC_chm,nclas_chm) ! Yk_left and Yk_right are known
  !!!
  !!!   integer(ip)                        :: iZ, ialfa, iclas
  !!!   real(rp)                           :: values(4), aux

  !!!   do iZ = 1, nZ_lr-2_ip
  !!!      ! Fill bottom and top
  !!!      Yk_interp(iZ,1,1:nclas_chm)                    = Yk_bottom(iZ+1_ip,1:nclas_chm)
  !!!      Yk_interp(iZ,n_alpha_grid_CMC_chm,1:nclas_chm) = Yk_top(iZ+1_ip,1:nclas_chm)

  !!!      ! Do interpolations
  !!!      do ialfa = 2, n_alpha_grid_CMC_chm-1_ip
  !!!         do iclas = 1, nclas_chm
  !!!            values(1) = (Yk_left(ialfa,iclas)  + Yk_bottom(iZ,iclas)) / 2.0_rp
  !!!            values(2) = (Yk_left(ialfa,iclas)  + Yk_top(iZ,iclas))    / 2.0_rp
  !!!            values(3) = (Yk_right(ialfa,iclas) + Yk_bottom(iZ,iclas)) / 2.0_rp
  !!!            values(4) = (Yk_right(ialfa,iclas) + Yk_top(iZ,iclas))    / 2.0_rp
  !!!            call bilinear_intepolation((/Z_left, Z_right, alpha_grid_CMC_chm(1), &
  !!!                   alpha_grid_CMC_chm(n_alpha_grid_CMC_chm)/), values, &
  !!!                   Z_CMC_chm(posZ_left+iZ), alpha_grid_CMC_chm(ialfa), aux)
  !!!            Yk_interp(iZ,ialfa,iclas) = aux
  !!!         end do
  !!!      end do
  !!!   end do

  !!!end subroutine chm_interpolate_HR_CMC


  subroutine chm_heatRelease_field_CMC

     !-----------------------------------------------------------------------
     !****f* Chemic/chm_integrate_chem_source_CMC
     ! NAME
     !    chm_integrate_chem_source_CMC
     ! DESCRIPTION
     !    Compute heat release per mass unit field at each mixture fraction level.
     !    heat_release = sum_j=1^N(par_mh_j * omega^molar_j)
     !    omega^molar_j = omega^mass_j / W_j = dY_j/dt * rho * volume. Substituting:
     !    heat_release = density * volume * sum_j=1^N(par_mh_j * dY_j/dt / W_j)
     !    However, the heat release per mass unit is obtained at each cell so volume
     !    and density are omitted.
     ! USED BY
     !
     !***
     !-----------------------------------------------------------------------

     use def_chemic,                 only : hrr_mass_cond_CMC_chm
#ifdef CANTERA
     use def_master,                 only : prthe
     use def_chemic,                 only : temp_CMC_chm,Yk_CMC_chm,src_Yk_CMC_chm,W_k,posZ_chem_integr_CMC_chm
#endif

     implicit none

#ifdef CANTERA
     integer(ip)                         :: iclas, ipoin, iZ
     real(rp)                            :: par_mh(nspec_chm) ! Partial molar enthalpy

     external                            :: cantera_partial_molar_enthalpies
#endif

     if (INOTMASTER) then

        hrr_mass_cond_CMC_chm = 0.0_rp
#ifdef CANTERA
        do iZ = posZ_chem_integr_CMC_chm(1), posZ_chem_integr_CMC_chm(2)
           do ipoin= 1, npoin
              ! Get molar enthalpies
              call cantera_partial_molar_enthalpies(nspec_chm,temp_CMC_chm(iZ,ipoin),prthe(1), &
                        Yk_CMC_chm(iZ,ipoin,1:nspec_chm),par_mh(1:nspec_chm))

              do iclas= 1, nspec_chm
                 hrr_mass_cond_CMC_chm(iZ,ipoin) = hrr_mass_cond_CMC_chm(iZ,ipoin) - &
                      par_mh(iclas) * src_Yk_CMC_chm(iZ,ipoin,iclas) / W_k(iclas)
              end do
           end do
        end do
#endif
     end if

  end subroutine chm_heatRelease_field_CMC



  subroutine chm_heatRelease_integral_CMC

     !-----------------------------------------------------------------------
     !****f* Chemic/chm_integrate_chem_source_CMC
     ! NAME
     !    chm_integrate_chem_source_CMC
     ! DESCRIPTION
     !    Compute heat release in the whole CMC domain.
     !    heat_release = sum_j=1^N(par_mh_j * omega^molar_j)
     !    omega^molar_j = omega^mass_j / W_j = dY_j/dt * rho * volume. Substituting:
     !    heat_release = density * volume * sum_j=1^N(par_mh_j * dY_j/dt / W_j)
     !    However, the heat release per volume unit is obtained at each cell so volume is omitted.
     ! USED BY
     !
     !***
     !-----------------------------------------------------------------------

     use def_domain,                 only : elmar, nelem, ltype, nnode, lnods, ngaus
     use mod_communications_global,  only : PAR_SUM
     use def_chemic,                 only : hrr_chm, hrr_int_chm

     implicit none
     integer(ip)                         :: ipoin, ielem, igaus, inode
     integer(ip)                         :: pnode, pgaus, pelty
     real(rp)                            :: gpvol, gpdet, gphre
     real(rp)                            :: elhre(mnode), elcod(ndime,mnode)
     real(rp)                            :: gpcar(ndime,mnode,mgaus)
     real(rp)                            :: xjaci(ndime,ndime),xjacm(ndime,ndime)

     external                            :: elmder
     external                            :: pararr


     if (INOTMASTER) then
        !
        ! Compute heat release integrated over the whole domain
        ! Transfer heat release field from nodes to Gaussian points and do the volumetric integral
        hrr_int_chm = 0.0_rp
        elements: do ielem = 1, nelem
           pelty = ltype(ielem)
           if (pelty > 0) then
              pnode = nnode(pelty)
              pgaus = ngaus(pelty)

              !
              ! Gather operations
              !
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 elhre(inode) = hrr_chm(ipoin)
                 elcod(1:ndime,inode) = coord(1:ndime,ipoin)
              end do

              !
              ! 1st and 2nd order Cartesian derivatives, and dV:=GPVOL=|J|*wg
              !
              do igaus = 1, pgaus
                 call elmder(&
                       pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&        ! Cartesian derivative
                       elcod(:,1:pnode),gpcar(1,1,igaus),gpdet,xjacm,xjaci)          ! and Jacobian
                  gpvol = elmar(pelty)%weigp(igaus) * gpdet               ! |J|*wg
                  gphre = 0.0_rp
                  do inode = 1, pnode
                     gphre = gphre + elmar(pelty) % shape(inode,igaus) * elhre(inode)
                  end do
                  hrr_int_chm = hrr_int_chm + gphre * gpvol
              end do
           end if
        end do elements

     else
        hrr_int_chm = 0.0_rp
     end if

     call PAR_SUM(hrr_int_chm)

  end subroutine chm_heatRelease_integral_CMC


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!  C L  I P P I N G S  !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine chm_limit_Yk_whole_domain_CMC
     !-----------------------------------------------------------------------
     !****f* Chemic/chm_limit_Yk_whole_domain_CMC
     ! NAME
     !    chm_limit_Yk_whole_domain_CMC
     ! DESCRIPTION
     !    Limit the values of mass fractions for the whole domain.
     ! USED BY
     !
     !***
     !-----------------------------------------------------------------------

     use def_chemic,             only :  nZ_CMC_chm, nclas_chm, Yk_CMC_chm

     implicit none
     integer(ip)                      :: imixf, ipoin

     if (INOTMASTER) then
        do imixf = 1, nZ_CMC_chm
           do ipoin = 1, npoin
              call chm_limit_Yk_CMC(nclas_chm, Yk_CMC_chm(imixf,ipoin,1:nclas_chm) )
           end do
        end do
     end if

  end subroutine chm_limit_Yk_whole_domain_CMC


  subroutine chm_limit_Yk_mixt_fraction_slice_CMC
     !-----------------------------------------------------------------------
     !****f* Chemic/chm_limit_Yk_mixt_fraction_slice_CMC
     ! NAME
     !    chm_limit_Yk_mixt_fraction_slice_CMC
     ! DESCRIPTION
     !    Limit the values of mass fractions for a mixture fraction slice.
     ! USED BY
     !
     !***
     !-----------------------------------------------------------------------

     use def_master,             only :  conce
     use def_chemic,             only :  nclas_chm

     implicit none
     integer(ip)                      :: ipoin

     if (INOTMASTER) then
        do ipoin = 1, npoin
           call chm_limit_Yk_CMC(nclas_chm, conce(ipoin,1:nclas_chm,1) )
        end do
     end if

  end subroutine chm_limit_Yk_mixt_fraction_slice_CMC


  subroutine chm_limit_Yk_CMC(nclas, Yk)
     !-----------------------------------------------------------------------
     !****f* Chemic/chm_limit_Yk_CMC
     ! NAME
     !    chm_limit_Yk_CMC
     ! DESCRIPTION
     !    Limit the values of mass fractions.
     ! USED BY
     !    chm_limit_Yk_whole_domain_CMC
     !***
     !-----------------------------------------------------------------------

     implicit none
     integer(ip), intent(in)          :: nclas
     real(rp),    intent(inout)       :: Yk(nclas)
     integer(ip)                      :: iclas
     real(rp)                         :: alpha, inv_alpha

     alpha = 0.0_rp
     do iclas = 1, nclas
        Yk(iclas) = max(0.0_rp, min(Yk(iclas),1.0_rp))
        alpha = alpha + Yk(iclas)
     end do
     inv_alpha = 1.0_rp / alpha
     Yk(1:nclas_chm) = inv_alpha * Yk(1:nclas_chm)

  end subroutine chm_limit_Yk_CMC


  subroutine clipping_var_from_CFD_CMC
     !-----------------------------------------------------------------------
     !****f* Chemic/clipping_var_from_CFD_CMC
     ! NAME
     !    clipping_var_from_CFD_CMC
     ! DESCRIPTION
     !    Clipping of the mixing variables coming from CFD. In principle such
     !    variables are clipped in CFD but due to a possible interpolation
     !    between CFD and CMC it is recommendable to clip at CMC again.
     ! USED BY
     !    chm_plugin
     !***
     !-----------------------------------------------------------------------

     use def_chemic,             only :  Zs_CMC_chm, Zavg_CMC_chm, Zvar_CMC_chm, &
                                         Xtot_CMC_chm

     implicit none
     ! To avoid too small time steps, scalar dissipation rate
     ! from CFD has to be upperly bounded
     real(rp), parameter              :: Xmax = 1.0e5_rp !!!!! PROVISIONAL
     real(rp), parameter              :: factor_Zv = 0.98_rp
     integer(ip)                      :: ipoin


     if (INOTMASTER) then
        do ipoin = 1, npoin
           Zavg_CMC_chm(ipoin) = max(min(Zavg_CMC_chm(ipoin), Zs_CMC_chm), 0.0_rp)
           Zvar_CMC_chm(ipoin) = max(min(Zvar_CMC_chm(ipoin), &
                                    factor_Zv*Zavg_CMC_chm(ipoin)*(Zs_CMC_chm-Zavg_CMC_chm(ipoin))), 0.0_rp)
           call clip_Zvar(Zavg_CMC_chm(ipoin), Zvar_CMC_chm(ipoin))
           Xtot_CMC_chm(ipoin) = max(min(Xtot_CMC_chm(ipoin), Xmax), 0.0_rp)
        end do
     end if

  end subroutine clipping_var_from_CFD_CMC


  subroutine clip_Zvar(Zavg, Zvar)
     !-----------------------------------------------------------------------
     !****f* Chemic/clip_Zvar
     ! NAME
     !    clip_Zvar
     ! DESCRIPTION
     !    Clip Zvar at the extremes to avoid false behaviours.
     ! USED BY
     !    clipping_var_from_CFD_CMC
     !***
     !-----------------------------------------------------------------------
     use def_chemic,        only :  Zlim_clipZvar_CMC_chm, Slim_clipZvar_CMC_chm, &
                                    Zs_CMC_chm

     implicit none
     real(rp), intent(in)        :: Zavg
     real(rp), intent(inout)     :: Zvar

     real(rp)                    :: S

     if (Zavg <= Zlim_clipZvar_CMC_chm(1) .or. Zavg >= Zlim_clipZvar_CMC_chm(2)) then
        S = Slim_clipZvar_CMC_chm(1)
     else
        S = Slim_clipZvar_CMC_chm(2)
     end if

     Zvar = min(Zvar, S*Zavg*(Zs_CMC_chm-Zavg))

  end subroutine clip_Zvar


  subroutine assign_val_clipZvar
     !-----------------------------------------------------------------------
     !****f* Chemic/assign_val_clipZvar
     ! NAME
     !    assign_val_clipZvar
     ! DESCRIPTION
     !    Assign mixt. frac. limits and segretation factor for Zvariance clipping.
     ! USED BY
     !    chm_initial_actions_reaphy_CMC
     !***
     !-----------------------------------------------------------------------
     use def_chemic,          only:  Zlim_clipZvar_CMC_chm, Slim_clipZvar_CMC_chm, &
                                     Zs_CMC_chm

     implicit none

     Zlim_clipZvar_CMC_chm(1) = 5.0e-04_rp
     Zlim_clipZvar_CMC_chm(2) = 0.99_rp * Zs_CMC_chm
     Slim_clipZvar_CMC_chm(1) = 0.5_rp        !!!!0.98_rp * S_threshold   ! Interpolation
     Slim_clipZvar_CMC_chm(2) = 0.5_rp        !!!0.3_rp

  end subroutine assign_val_clipZvar


  subroutine chm_limit_extremes_PDF_CMC
     !-----------------------------------------------------------------------
     !****f* Chemic/chm_limit_extremes_PDF_CMC
     ! NAME
     !    chm_limit_extremes_PDF_CMC
     ! DESCRIPTION
     !    Estimate limits at mixture fraction extremes (Z=0 and Z=Zs) to know
     !    when to apply the models
     ! USED BY
     !    chm_element_operations_CMC
     !    chm_updtcc_CMC
     !***
     !-----------------------------------------------------------------------

     use def_chemic,             only : extremes_Z_CMC_chm, Zs_CMC_chm

     implicit none
     !!!!!real(rp), parameter              :: fact = 0.8_rp, eps = 5.0e-3_rp
     !!!!! fact has to be higher than 0.5 to avoid errors when dealing with the extremes
     !!!!! and the PDF when deciding if the models have to be applied

     !!!!!extremes_Z_CMC_chm(1) = max(fact * Z_CMC_chm(2), eps * Zs_CMC_chm)
     !!!!!extremes_Z_CMC_chm(2) = min(fact*Z_CMC_chm(nZ_CMC_chm-1) + (1.0_rp - fact)*Zs_CMC_chm, &
     !!!!!                            (1.0_rp - eps) * Zs_CMC_chm)

     extremes_Z_CMC_chm(1) = 5.0e-4_rp
     extremes_Z_CMC_chm(2) = Zs_CMC_chm - 5.0e-4_rp

  end subroutine chm_limit_extremes_PDF_CMC


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!  O T H E R  F U N C T I O N S  !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine chm_values_cvgunk_CMC
     !-----------------------------------------------------------------------
     !****f* Chemic/chm_values_cvgunk_CMC
     ! NAME
     !    chm_values_cvgunk_CMC
     ! DESCRIPTION
     !    Obtain some values of interest for the convergence file
     ! USED BY
     !    chm_endste
     !***
     !-----------------------------------------------------------------------
     use def_chemic,      only :  Text_cond_CMC_chm, Text_uncond_CMC_chm, &
                                  temp_int_CMC_chm, temp_CMC_chm, nZ_CMC_chm, &
                                  temp_inert_CMC_chm, sum_Yk_ext_cond_CMC_chm, &
                                  Yk_CMC_chm

     implicit none
     integer(ip)               :: iZ, ipoin, iclas
     real(rp)                  :: alpha

     external                  :: pararr

     sum_Yk_ext_cond_CMC_chm = 1.0_rp

     if (INOTMASTER) then
        ! Find extreme values for unconditional temperature
        Text_uncond_CMC_chm(1) = minval(temp_int_CMC_chm)
        Text_uncond_CMC_chm(2) = maxval(temp_int_CMC_chm)

        ! Find extreme values for conditional temperature
        do iZ = 1, nZ_CMC_chm
           Text_cond_CMC_chm(2*iZ-1) = minval(temp_CMC_chm(iZ,:))
           Text_cond_CMC_chm(2*iZ)   = maxval(temp_CMC_chm(iZ,:))
        end do

        do ipoin = 1, npoin
           do iZ = 2, nZ_CMC_chm-1
              alpha = 0.0_rp
              do iclas = 1, nclas_chm
                 alpha = alpha + Yk_CMC_chm(iZ,ipoin,iclas)
              end do
              if (alpha < sum_Yk_ext_cond_CMC_chm(2*iZ-1))  &
                 sum_Yk_ext_cond_CMC_chm(2*iZ-1) = alpha
              if (alpha > sum_Yk_ext_cond_CMC_chm(2*iZ))    &
                 sum_Yk_ext_cond_CMC_chm(2*iZ) = alpha
           end do
        end do
     else
        Text_uncond_CMC_chm(1) = temp_inert_CMC_chm(1)
        Text_uncond_CMC_chm(2) = temp_inert_CMC_chm(1)

        do iZ = 1, nZ_CMC_chm
           Text_cond_CMC_chm(2*iZ-1) = temp_inert_CMC_chm(iZ)
           Text_cond_CMC_chm(2*iZ)   = temp_inert_CMC_chm(iZ)
        end do

     end if

     ! Look for maximum over subdomains
     
     call PAR_MIN(Text_uncond_CMC_chm(1))
     call PAR_MAX(Text_uncond_CMC_chm(2))

     do iZ = 1, nZ_CMC_chm
        call PAR_MIN(Text_cond_CMC_chm(2*iZ-1))
        call PAR_MAX(Text_cond_CMC_chm(2*iZ))
        call PAR_MIN(sum_Yk_ext_cond_CMC_chm(2*iZ-1))
        call PAR_MAX(sum_Yk_ext_cond_CMC_chm(2*iZ))
     end do

  end subroutine chm_values_cvgunk_CMC


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! M A T H E M A T I C A L   F U N C T I O N S !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!! SHOULDN'T THESE FUNCTIONS BE IN THE KERNEL?


  subroutine incob ( a, b, x, bt, bix )

     !*****************************************************************************
     !
     !  INCOB computes the regularized incomplete beta function Ix(a,b). Hence, if we
     !  denote the incomplete beta function as
     !
     !  B(x; a, b) = int_0^x (t^(a-1) * (1-t)^(b-1) * dt)
     !
     !  and the beta function as
     !
     !  B(a,b) = int_0^1 (t^(a-1) * (1-t)^(b-1) * dt)
     !
     !  then the regularized incomplete beta function Ix(a,b) is
     !
     !  Ix(a,b) = B(x; a, b) / B(a,b)
     !
     !
     !  Licensing:
     !
     !    This routine is copyrighted by Shanjie Zhang and Jianming Jin. However,
     !    they give permission to incorporate this routine into a user program
     !    provided that the copyright is acknowledged.
     !
     !  Modified:
     !
     !    22 July 2012
     !
     !  Author:
     !
     !    Shanjie Zhang, Jianming Jin.
     !
     !  Reference:
     !
     !    Shanjie Zhang, Jianming Jin,
     !    Computation of Special Functions,
     !    Wiley, 1996,
     !    ISBN: 0-471-11963-6,
     !    LC: QA351.C45.
     !
     !  Parameters:
     !    a, b: parameters of the beta function
     !    bt: value of B(a,b)
     !    x: value for which evaluate the incomplete beta function
     !    bix: value of the incomplete beta function

     implicit none
     real(rp), intent(in)     :: a, b, x, bt
     real(rp), intent(out)    :: bix

     integer(ip)              :: k
     real(rp)                 :: dk(51), fk(51), s0
     real(rp)                 :: t1, t2, ta, tb


     ! ¿CAMBIAR LA PRECISIÓN DE LOS DOUBLE?

     s0 = ( a + 1.0D+00 ) / ( a + b + 2.0D+00 )

     if ( x <= s0 ) then

       do k = 1, 20
         dk(2*k) = real(k,rp) * ( b - real(k,rp) ) * x / &
           ( a + 2.0D+00 * real(k,rp) - 1.0D+00 ) / ( a + 2.0D+00 * real(k,rp) )
       end do

       do k = 0, 20
         dk(2*k+1) = - ( a + real(k,rp) ) * ( a + b + real(k,rp) ) * x &
           / ( a + 2.0D+00 * real(k,rp) ) / ( a + 2.0D+00 * real(k,rp) + 1.0D+00 )
       end do

       t1 = 0.0D+00
       do k = 20, 1, -1
         t1 = dk(k) / ( 1.0D+00 + t1 )
       end do
       ta = 1.0D+00 / ( 1.0D+00 + t1 )
       bix = x ** a * ( 1.0D+00 - x ) ** b / ( a * bt ) * ta

     else

       do k = 1, 20
         fk(2*k) = real(k,rp) * ( a - real(k,rp) ) * ( 1.0D+00 - x ) &
           / ( b + 2.0D+00 * real(k,rp) - 1.0D+00 ) / ( b + 2.0D+00 * real(k,rp) )
       end do

       do k = 0,20
         fk(2*k+1) = - ( b + real(k,rp) ) * ( a + b + real(k,rp) ) * ( 1.0D+00 - x ) &
           / ( b + 2.0D+00 * real(k,rp) ) / ( b + 2.0D+00 * real(k,rp) + 1.0D+00 )
       end do

       t2 = 0.0D+00
       do k = 20, 1, -1
         t2 = fk(k) / ( 1.0D+00 + t2 )
       end do
       tb = 1.0D+00 / ( 1.0D+00 + t2 )
       bix = 1.0D+00 - x ** a * ( 1.0D+00 - x ) ** b / ( b * bt ) * tb

     end if

     return

  end subroutine incob



  subroutine beta_func(x, y, beta_f)

     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/beta_func
     ! NAME
     !    beta_func
     ! DESCRIPTION
     !    It computes beta function. All the values are assumed to be strictly
     !    positive.
     ! USED BY
     !    chm_find_PDF_parameters_CMC
     !***
     !-----------------------------------------------------------------------

     implicit none
     real(rp), intent(in)            :: x, y
     real(rp), intent(out)           :: beta_f

     beta_f = gamma(x) * gamma(y) / gamma(x+y)

  end subroutine beta_func



  subroutine compute_erfinv(y0,x)

     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/compute_erfinv
     ! NAME
     !    compute_erfinv
     ! DESCRIPTION
     !    Compute inverse error function by a bisection method.
     ! USED BY
     !    compute_Xnormalized_profile_CMC
     !***
     !-----------------------------------------------------------------------

     implicit none
     real(rp), parameter             :: small=1.0e-8_rp, big=1.0e20_rp

     real(rp), intent(in)            :: y0
     real(rp), intent(out)           :: x

     real(rp)                        :: y, xmiddle, x_aux, residual1, residual2, xmin, xmax
     logical                         :: negative, convergence

     xmin        = 0.0_rp
     xmax        = 4.0_rp
     negative    = .false.
     convergence = .false.

     if (y0 < 0.0_rp) then
        negative = .true.
        y = -y0
     else
        y = y0
     end if

     if (y >= erf(xmax)) then
        x = big
     else if (y==0.0_rp) then
        x = 0.0_rp
     else
        do while (.not. convergence)
           xmiddle = (xmin + xmax) / 2.0_rp
           residual1 =  y - erf(xmiddle)

           if (residual1 == 0.0_rp) then
              x = xmiddle
              convergence = .true.
           else
              if (residual1 < 0.0_rp) then
                 xmax = xmiddle
                 residual2 = y-erf(xmin)
                 x_aux = xmin
              else
                 xmin = xmiddle
                 residual2 = y-erf(xmax)
                 x_aux = xmax
              end if

              if (abs(residual1) < small) then
                 x = xmiddle
                 convergence = .true.
              else if (abs(residual2) < small) then
                 x = x_aux
                 convergence = .true.
              end if
           end if

        end do
     end if

     if (negative) then
        x = -x
     end if

  end subroutine compute_erfinv



  subroutine compute_1st_deriv_2nd_order(n_x, incr_x, y, y_prime)

     !-----------------------------------------------------------------------
     !****f* Chemic/mod_chm_operations_CMC/compute_1st_deriv_2nd_order
     ! NAME
     !    compute_1st_deriv_2nd_order
     ! DESCRIPTION
     !    Compute first derivative with second order excluding the extremes.
     ! USED BY
     !
     !***
     !-----------------------------------------------------------------------

     implicit none
     integer(ip), intent(in)     :: n_x             ! Total number of points in the mesh
     real(rp),    intent(in)     :: incr_x(n_x-1)   ! Step between points
     real(rp),    intent(in)     :: y(n_x)          ! Ordinates of the function
     real(rp),    intent(out)    :: y_prime(n_x-2)  ! First derivative

     integer(ip)                 :: i
     real(rp)                    :: coeff_left, coeff_cent, coeff_right
     real(rp)                    :: inv_left, inv_cent, inv_right

     do i = 2, n_x-1
        inv_left     = 1.0_rp / incr_x(i-1)
        inv_cent     = 1.0_rp / (incr_x(i-1) + incr_x(i))
        inv_right    = 1.0_rp / incr_x(i)
        coeff_left   = -incr_x(i) * inv_left * inv_cent
        coeff_cent   = (incr_x(i) - incr_x(i-1)) * inv_left * inv_right
        coeff_right  = incr_x(i-1) * inv_cent * inv_right
        y_prime(i-1) = y(i-1) * coeff_left + y(i) * coeff_cent + y(i+1) * coeff_right
     end do

  end subroutine compute_1st_deriv_2nd_order

end module mod_chm_operations_CMC
