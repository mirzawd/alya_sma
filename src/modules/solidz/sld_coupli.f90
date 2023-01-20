!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_coupli.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   This routine manages solidz coupling with other modules
!> @details This routine manages solidz coupling with other modules
!> @} 
!-----------------------------------------------------------------------

subroutine sld_coupli(itask)

  use def_kintyp,  only : ip
  use def_master,  only : ITASK_INIUNK,ITASK_MATRIX
  use def_master,  only : ITASK_BEGITE,ITASK_ENDITE
  use mod_sld_atm, only : kfl_therm_sld
  use mod_sld_atm, only : sld_atm_set_initial_temperature_at_POINTS
  
  implicit none

  integer(ip), intent(in) :: itask   !< who is calling sld_coupli

  select case( itask )

  case ( ITASK_INIUNK )

     !-------------------------------------------------------------------
     !
     ! INIUNK
     !
     !-------------------------------------------------------------------
     !
     ! Thermomechanical coupling
     !
     if( kfl_therm_sld ) then
        call sld_atm_set_initial_temperature_at_POINTS()
     endif

  case ( ITASK_BEGITE )

     !-------------------------------------------------------------------
     !
     ! BEGITE
     !
     !-------------------------------------------------------------------

  case ( ITASK_ENDITE )

     !-------------------------------------------------------------------
     !
     ! ENDITE
     !
     !-------------------------------------------------------------------

  case ( ITASK_MATRIX )

     !-------------------------------------------------------------------
     !
     ! MATRIX: After assembly
     !
     !-------------------------------------------------------------------

  end select

end subroutine sld_coupli
