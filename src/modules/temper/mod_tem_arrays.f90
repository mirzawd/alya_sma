!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Temper
!> @{
!> @file    mod_tem_arrays.f90
!> @author  houzeaux
!> @date    2019-11-16
!> @brief   Temper arrays
!> @details Temper arrays
!-----------------------------------------------------------------------

module mod_tem_arrays

  use def_master
  use def_domain 
  use def_temper
  use def_kermod
  use def_solver
  use mod_output_postprocess,  only : output_postprocess_check_variable_postprocess
  use mod_arrays,              only : arrays
  use mod_arrays,              only : arrays_number
  use mod_memory,              only : memory_alloca
  use mod_ADR,                 only : ADR_initialize_type
  use mod_ADR,                 only : ADR_check_and_compute_data
  use mod_ADR,                 only : ADR_allocate_projections_bubble_sgs
  use mod_ADR,                 only : ADR_arrays

  implicit none

  private

  public :: tem_arrays
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-16
  !> @brief   Temper arrays
  !> @details Do what you have to do with temper arrays
  !>          THERM is the module unknown, in principle the temperature.
  !>          When kfl_regim_tem = 4, THERM is the
  !>          enthalpy and we need to allocate TEMPE as well
  !> 
  !-----------------------------------------------------------------------
  
  subroutine tem_arrays(wtask)

    character(len=*), intent(in) :: wtask
    !
    ! Finite element vs finite volume
    !
    if(      kfl_discr_tem == 0 ) then
       solve(1) % kfl_where = SOL_NODES
       nunkn_tem            = npoin
    else if( kfl_discr_tem == 1 ) then
       solve(1) % kfl_where = SOL_ELEMENTS
       nunkn_tem            = nelem
    end if
    !
    ! ADR type
    !
    call ADR_check_and_compute_data(ADR_tem)
    call ADR_arrays(wtask,ADR_tem,'BUBBT','PROJ1','PROJ2','TESG2')
    !
    ! TEMPER: Temperature unknown 
    !
    call arrays(arrays_number('THERM'),wtask,therm,nunkn_tem,ncomp_tem)
    if (kfl_regim_tem == 4) then
       call arrays(arrays_number('TEMPE'),wtask,tempe,nunkn_tem,ncomp_tem)       
    else
       tempe => therm
    end if
    !
    ! TEOLD_TEM: Old temperaure
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='RESID') ) then
       call arrays(arrays_number('RESID'),wtask,teold_tem,nunkn_tem)
    end if
    !
    ! AVTEM_TEM: Average temperature
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVTEM') ) then     
       call arrays(arrays_number('AVTEM'),wtask,avtem_tem,nunkn_tem)
    end if
    !
    ! FVTEM_TEM: Favre-average temperature
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='FVTEM') ) then     
       call arrays(arrays_number('FVTEM'),wtask,fvtem_tem,nunkn_tem)
    end if
    !
    ! FVTE2_TEM: Favre-average temperature squared
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='FVTE2') ) then     
       call arrays(arrays_number('FVTE2'),wtask,fvte2_tem,nunkn_tem)
    end if
    !
    ! AVTE2_TEM: Average tempe**2
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVTE2') ) then     
       call arrays(arrays_number('AVTE2'),wtask,avte2_tem,nunkn_tem)
    end if
    !
    ! AVTEV_TEM: Average tempe*veloc
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVTEV') ) then     
       call arrays(arrays_number('AVTEV'),wtask,avtev_tem,ndime, nunkn_tem)
    end if
    !
    ! AVDEN_TEM: Average density
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVDEN') ) then     
       call arrays(arrays_number('AVDEN'),wtask,avden_tem,nunkn_tem)
    end if
    !
    ! FVVEL_TEM: Average rho*veloc
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='FVVEL') ) then  
       call arrays(arrays_number('FVVEL'),wtask,fvvel_tem,ndime, nunkn_tem)
    end if
    !
    ! AVRES_TEM: Average residual heat flux
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVRES') ) then  
       call arrays(arrays_number('AVRES'),wtask,avres_tem,nunkn_tem)
    end if
    !
    ! AVHSK_TEM: Average spray heat source
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVHSK') ) then  
       call arrays(arrays_number('AVHSK'),wtask,avhsk_tem,nunkn_tem)
    end if

  end subroutine tem_arrays
   
end module mod_tem_arrays
!> @}
