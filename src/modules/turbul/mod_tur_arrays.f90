!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Turbul
!> @{
!> @file    mod_tur_arrays.f90
!> @author  houzeaux
!> @date    2022-09-20
!> @brief   Turbul arrays
!> @details Turbul arrays
!-----------------------------------------------------------------------

module mod_tur_arrays

  use def_master
  use def_domain 
  use def_turbul
  use def_kermod
  use def_solver
  use mod_output_postprocess,  only : output_postprocess_check_variable_postprocess
  use mod_arrays,              only : arrays
  use mod_arrays,              only : arrays_number
  use mod_memory,              only : memory_alloca

  implicit none

  private

  public :: tur_arrays
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-16
  !> @brief   Turbul arrays
  !> @details Do what you have to do with turbul arrays
  !>          THERM is the module unknown, in principle the turbulature.
  !>          When kfl_regim_tem = 4, THERM is the
  !>          enthalpy and we need to allocate TEMPE as well
  !> 
  !-----------------------------------------------------------------------
  
  subroutine tur_arrays(wtask) 

    character(len=*), intent(in) :: wtask
    !
    ! TURBUL: Turbulent unknown and viscosity
    !
    call arrays(arrays_number('UNTUR'),wtask,untur,nturb_tur,npoin,ncomp_tur)
    call arrays(arrays_number('TURMU'),wtask,turmu,npoin)
    !
    ! Averages
    !
    if( output_postprocess_check_variable_postprocess(arrays_number('AVKEY')) ) then
       call arrays(arrays_number('AVKEY'),wtask,avkey_tur,npoin)
     end if
     if( output_postprocess_check_variable_postprocess(arrays_number('AVOME')) ) then
       call arrays(arrays_number('AVOME'),wtask,avome_tur,npoin)
    end if
     !
     ! averaged values
     !
    if( output_postprocess_check_variable_postprocess(arrays_number('AVTVI')) ) then
       call arrays(arrays_number('AVTVI'),wtask,avtvi_tur,npoin)
       call arrays(arrays_number('OLDED'),wtask,avtvi_tur,npoin)
    end if
    !
    ! UNPRO_TUR: Residual projections for OSS methods
    !
    if( kfl_ortho_tur >= 1 ) then
       call arrays(arrays_number('UNPRO'),wtask,unpro_tur,nturb_tur,npoin)
       if( kfl_ortho_tur == 2 ) then  ! split oss
          call arrays(arrays_number('UNPRR'),wtask,unprr_tur,nturb_tur,npoin)
       end if
    end if

  end subroutine tur_arrays
   
end module mod_tur_arrays
!> @}
