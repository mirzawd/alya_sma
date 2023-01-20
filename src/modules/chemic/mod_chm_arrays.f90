!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Chemic
!> @{
!> @file    mod_chm_arrays.f90
!> @author  houzeaux
!> @date    2019-11-16
!> @brief   Nastin arrays
!> @details Nastin arrays
!-----------------------------------------------------------------------

module mod_chm_arrays

  use def_master,              only : mem_modul, modul, conce
  use def_domain,              only : ndime, npoin
  use def_chemic,              only : av_Z_flux_chm, enthalp_CMC_chm, kfl_model_chm, kfl_solve_cond_CMC_chm,&
                                      kfl_solve_enth_CMC_chm, nclas_chm, ncomp_chm, nZ_CMC_chm, temp_CMC_chm, Yk_CMC_chm, ADR_chm
  use def_kintyp,              only : ip
  use mod_output_postprocess,  only : output_postprocess_check_variable_postprocess
  use mod_arrays,              only : arrays
  use mod_arrays,              only : arrays_number
  use mod_memory,              only : memory_alloca
  use mod_ADR,                 only : ADR_check_and_compute_data
  use mod_ADR,                 only : ADR_arrays

  implicit none

  private

  public :: chm_arrays

contains

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2019-11-16
  !> @brief   Nastin arrays
  !> @details Do what you have to do with nastin arrays
  !>
  !-----------------------------------------------------------------------

  subroutine chm_arrays(wtask)

    character(len=*), intent(in) :: wtask
    integer(ip)                  :: iclas
    !
    ! CONCE
    !
    call arrays(arrays_number('CONCE'),wtask,conce,npoin,nclas_chm,ncomp_chm)
    if( wtask == 'ALLOCATE' .and. npoin == 0 ) then
       call memory_alloca(mem_modul(1:2,modul),'CONCE','chm_memall',conce,1_ip,nclas_chm,ncomp_chm)
    end if
    !
    ! ADR type
    !
    do iclas=1,nclas_chm
       call ADR_check_and_compute_data(ADR_chm(iclas))
       !call ADR_arrays(wtask,ADR_chm(iclas),'BUBBT','PROJ1','PROJ2','TESG2',TAG1=iclas)
    end do

    !
    ! CMC model
    !
    if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1) then
       call arrays(arrays_number('MASCN'),wtask,Yk_CMC_chm,nZ_CMC_chm,npoin,nclas_chm)
       call arrays(arrays_number('TEMCN'),wtask,temp_CMC_chm,nZ_CMC_chm,npoin)
       if( wtask == 'ALLOCATE' .and. npoin == 0 ) then
          call memory_alloca(mem_modul(1:2,modul),'CONDITIONAL_MASS_FRACTIONS_CMC_CHM','chm_solmem',Yk_CMC_chm,nZ_CMC_chm,1_ip,&
              nclas_chm)
          call memory_alloca(mem_modul(1:2,modul),'CONDITIONAL_TEMPERATURE_CMC_CHM','chm_solmem',temp_CMC_chm,nZ_CMC_chm,1_ip)
       end if

       if (kfl_solve_enth_CMC_chm /= 0) then
          call arrays(arrays_number('ENTCN'),wtask,enthalp_CMC_chm,nZ_CMC_chm,npoin)
          if( wtask == 'ALLOCATE' .and. npoin == 0 ) then
             call memory_alloca(mem_modul(1:2,modul),'CONDITIONAL_ENTHALPY_CMC_CHM','chm_solmem',enthalp_CMC_chm,nZ_CMC_chm,1_ip)
          end if
       end if

       !!!if (kfl_bc_init_method_CMC_chm == 2) then
       !!!   call arrays(arrays_number('ALPHA'),wtask,alpha_val_spec_CMC_chm,1_ip,npoin, nclas_chm)
       !!!   if( wtask == 'ALLOCATE' .and. npoin == 0 ) then
       !!!      call memory_alloca(mem_modul(1:2,modul),'ALPHA_VALUES_SPECIES_CHM','chm_solmem',alpha_val_spec_CMC_chm,1_ip,1_ip,nclas_chm)
       !!!   end if
       !!!end if

    end if
    !
    ! AVZFL: average Z flux
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVZFL') ) then
       call arrays(arrays_number('AVZFL'),trim(wtask),av_Z_flux_chm,ndime,npoin)
    end if

  end subroutine chm_arrays

end module mod_chm_arrays
!> @}
