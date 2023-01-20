!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_memphy(itask)
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_memphy
  ! NAME
  !    chm_memphy
  ! DESCRIPTION
  !    Allocate memory for the physical problem
  ! OUTPUT
  ! USES
  ! USED BY
  !    chm_reaphy
  !***
  !-----------------------------------------------------------------------
  use def_master,       only : mem_modul, modul, speci
  use def_chemic,       only : alpha_grid_CMC_chm, diff_Z_CMC_chm, diffu_chm, kfl_annfw_src_equa_list, kfl_fw_src_equa_list, Le_k,&
                               max_ann_fw, max_diff_Yk_inert_eq_CMC_chm, max_lookup_fw, mixedEq_eqs_chm, mixedEq_groups_chm,&
                               n_alpha_grid_CMC_chm, nclas_chm, ngrou_chm, nS_AMC_CMC_chm, nspec_chm, nZ_AMC_CMC_chm, nZ_CMC_chm,&
                               PDF_min_CMC_chm, react_scal_bc_CMC_chm, rscal_equil_CMC_chm, rscal_inert_CMC_chm, S_AMC_CMC_chm,&
                               T_bc_CMC_chm, temp_equil_CMC_chm, temp_inert_CMC_chm, W_k, Xintegrated_table_AMC_CMC_chm,&
                               Xnormalized_prof_CMC_chm, Z_AMC_CMC_chm, Z_CMC_chm, exp_alpha_grid_CMC_chm
  use def_kintyp,       only : ip, rp
  use def_parame,       only : two
  use mod_memory,       only : memory_alloca
  use mod_memory,       only : memory_deallo
  use mod_chm_mixedEq,  only : chm_mixedEq_initialize
  implicit none
  integer(ip), intent(in) :: itask
  integer(4)              :: istat
  integer(ip)             :: ii

  external                :: memerr

  select case(itask)

  case(1_ip)
     !
     ! Allocate properties
     !
     call memory_alloca(mem_modul(1:2,modul),'DIFFU_CHM','chm_memphy',diffu_chm,2_ip,nclas_chm)

     allocate(speci(nclas_chm)) !!**

     call memory_alloca(mem_modul(1:2,modul),'LE_K','chm_memphy',Le_k,nclas_chm)
     if (nclas_chm > 0_ip) Le_k = 1.0_rp

  case(2_ip)
     !
     ! Deallocate properties
     !
     call memory_deallo(mem_modul(1:2,modul),'DIFFU_CHM','chm_memphy',diffu_chm)

     deallocate(speci,stat=istat) !!**
     if(istat/=0) call memerr(two,'speci','chm_memphy',0_ip)

     call memory_deallo(mem_modul(1:2,modul),'LE_K','chm_memphy',Le_k)

  case(3_ip)
     !
     ! Allocate mixture fraction vector for CMC model
     !
     nullify(Z_CMC_chm)
     nullify(diff_Z_CMC_chm)
     call memory_alloca(mem_modul(1:2,modul),'MIXT_FRAC_VEC_CMC_CHM','chm_memphy',Z_CMC_chm,nZ_CMC_chm)
     call memory_alloca(mem_modul(1:2,modul),'MIXT_FRAC_DIFFERENCE_VEC_CMC_CHM','chm_memphy',diff_Z_CMC_chm,nZ_CMC_chm-1_ip)
     Z_CMC_chm(1:nZ_CMC_chm) = 0.0_rp
     diff_Z_CMC_chm(1:nZ_CMC_chm-1) = 0.0_rp


  case(4_ip)
     !
     ! Allocate memory for reactive scalars (temperature, enthalpy and mass fractions)
     !  at boundaries (mixture fraction space) for CMC model
     !
     nullify(T_bc_CMC_chm)
     nullify(react_scal_bc_CMC_chm)
     call memory_alloca(mem_modul(1:2,modul),'TEMPERATURE_BC_CMC_CHM','chm_memphy',T_bc_CMC_chm,2_ip)
     call memory_alloca(mem_modul(1:2,modul),'REACTIVE_SCALARS_BC_CMC_CHM','chm_memphy',react_scal_bc_CMC_chm,nspec_chm+1_ip,2_ip)
     T_bc_CMC_chm(1:2) = 298.0_rp   ! Temperature
     react_scal_bc_CMC_chm(1:nspec_chm+1,1:2) = 0.0_rp    ! Enthalpy and mass fractions


  case(5_ip)
     !
     ! Allocate memory for mixture fraction and segregation factor for AMC model for CMC model
     !
     nullify(Z_AMC_CMC_chm)
     nullify(S_AMC_CMC_chm)
     nullify(Xintegrated_table_AMC_CMC_chm)
     call memory_alloca(mem_modul(1:2,modul),'MIXTURE_FRAC_AMC_MODEL_CMC_CHM','chm_memphy',Z_AMC_CMC_chm,nZ_AMC_CMC_chm)
     call memory_alloca(mem_modul(1:2,modul),'SEGREGATION_FACTOR_AMC_MODEL_CMC_CHM','chm_memphy',S_AMC_CMC_chm,nS_AMC_CMC_chm)
     call memory_alloca(mem_modul(1:2,modul),'SCAL_DISSIP_RATE_INTEGRATED_PROFILE_CMC_CHM','chm_memphy',&
         Xintegrated_table_AMC_CMC_chm,nZ_AMC_CMC_chm,nS_AMC_CMC_chm)

     Z_AMC_CMC_chm(1:nZ_AMC_CMC_chm) = 0.0_rp
     S_AMC_CMC_chm(1:nS_AMC_CMC_chm) = 0.0_rp
     Xintegrated_table_AMC_CMC_chm(1:nZ_AMC_CMC_chm,1:nS_AMC_CMC_chm) = 0.0_rp


  case(6_ip)
     !
     ! Allocate memory for scalar dissipation rate normalized profile
     !
     nullify(Xnormalized_prof_CMC_chm)
     call memory_alloca(mem_modul(1:2,modul),'SCALAR_DISSIP_RATE_NORMALIZED_PROF_CMC_CHM','chm_memphy',Xnormalized_prof_CMC_chm,&
         nZ_CMC_chm)
     Xnormalized_prof_CMC_chm(1:nZ_CMC_chm) = 0.0_rp

  case(7_ip)
     !
     ! Allocate memory for inert and equilibrium profiles and molecular weight.
     ! Note that W_k is only allocated for the master while later in chm_memall
     ! is allocated for the slaves
     !
     nullify(W_k)
     nullify(rscal_inert_CMC_chm)
     nullify(rscal_equil_CMC_chm)
     nullify(temp_inert_CMC_chm)
     nullify(temp_equil_CMC_chm)
     nullify(max_diff_Yk_inert_eq_CMC_chm)

     call memory_alloca(mem_modul(1:2,modul),'MOLECULAR_WEIGHT_CMC_CHM','chm_memphy',W_k,nspec_chm)
     call memory_alloca(mem_modul(1:2,modul),'INERT_PROFILE_CMC_CHM','chm_memphy',rscal_inert_CMC_chm,nZ_CMC_chm,nspec_chm+1_ip)
     call memory_alloca(mem_modul(1:2,modul),'EQUILIBRIUM_PROFILE_CMC_CHM','chm_memphy',rscal_equil_CMC_chm,nZ_CMC_chm,&
         nspec_chm+1_ip)
     call memory_alloca(mem_modul(1:2,modul),'INERT_TEMPERATURE_CMC_CHM','chm_memphy',temp_inert_CMC_chm,nZ_CMC_chm)
     call memory_alloca(mem_modul(1:2,modul),'EQUILIBRIUM_TEMPERATURE_CMC_CHM','chm_memphy',temp_equil_CMC_chm,nZ_CMC_chm)
     call memory_alloca(mem_modul(1:2,modul),'MAX_DIFF_YK_INERT_EQUILIBRIUM_CMC_CHM','chm_memphy',max_diff_Yk_inert_eq_CMC_chm,&
         nspec_chm)

     W_k                          = 0.0_rp
     rscal_inert_CMC_chm          = 0.0_rp
     rscal_equil_CMC_chm          = 0.0_rp
     temp_inert_CMC_chm           = 0.0_rp
     temp_equil_CMC_chm           = 0.0_rp
     max_diff_Yk_inert_eq_CMC_chm = 0.0_rp


  case(8_ip)
     !
     ! Allocate memory for PDF_min_CMC_chm
     !
     nullify(PDF_min_CMC_chm)
     call memory_alloca(mem_modul(1:2,modul),'MINIMUM_PDF_VALUE_CMC_CHM','chm_memall',PDF_min_CMC_chm,nZ_CMC_chm)
     PDF_min_CMC_chm(1:nZ_CMC_chm) = 0.0_rp

  case(9_ip)
     !
     ! Allocate equation structures
     !
     allocate(mixedEq_groups_chm(ngrou_chm))
     allocate(mixedEq_eqs_chm(nclas_chm))
     call chm_mixedEq_initialize(ngrou_chm,nclas_chm,mixedEq_groups_chm,mixedEq_eqs_chm)

     call memory_alloca(mem_modul(1:2,modul),'KFL_FW_SRC_EQUA_LIST','chm_memphy',kfl_fw_src_equa_list,max_lookup_fw,nclas_chm)
     call memory_alloca(mem_modul(1:2,modul),'KFL_ANNFW_SRC_EQUA_LIST','chm_memphy',kfl_annfw_src_equa_list,max_ann_fw,nclas_chm)

  case(11_ip)
     !
     ! Allocate memory for alpha_grid_CMC_chm
     !
     nullify(alpha_grid_CMC_chm)
     call memory_alloca(mem_modul(1:2,modul),'ALFA_DISTRIBUTION_HR_CMC_CHM','chm_memall',alpha_grid_CMC_chm,n_alpha_grid_CMC_chm)
     alpha_grid_CMC_chm(1) = 0.0_rp
     alpha_grid_CMC_chm(n_alpha_grid_CMC_chm) = 1.0_rp
     do ii = 2, n_alpha_grid_CMC_chm-1_ip
        alpha_grid_CMC_chm(ii) = ( real(ii-1_ip,rp) / real(n_alpha_grid_CMC_chm-1_ip,rp) ) ** exp_alpha_grid_CMC_chm
     end do

  end select

end subroutine chm_memphy

