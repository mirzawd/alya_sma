!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Coupling
!> @{
!> @file    Endzon.f90
!> @author  Guillaume Houzeaux
!> @date    30/09/2014
!> @brief   Check block convergence if coupling
!> @details Check block convergence if coupling
!> @} 
!-----------------------------------------------------------------------

subroutine Endzon()
  
  use def_master,    only : iblok
  use def_master,    only : kfl_gocou
  use def_kintyp,    only : ip
  use def_master,    only : ITASK_ENDZON
  use def_coupli,    only : mcoup
  use def_coupli,    only : kfl_gozon
  use def_coupli,    only : coupling_driver_iteration
  use def_coupli,    only : coupling_driver_number_couplings
  use mod_couplings, only : COU_CHECK_CONVERGENCE
  use mod_moduls,    only : moduls 
  use mod_immersed,  only : cou_update_couplings
  use def_coupli,    only : kfl_efect
  
  implicit none
  
  !
  ! Check convergence if we have coupling
  !
  if( mcoup > 0_ip ) then

     !
     ! Call modules to exchange data
     !
     call moduls(ITASK_ENDZON)
     !
     ! Update couplings when using Embedded Finite Element Coupling Method (EFECT)
     !
     if( kfl_efect .and. coupling_driver_iteration(iblok) > 0_ip ) call cou_update_couplings()
     !
     ! We are currently in block IBLOK
     !  
     if( coupling_driver_number_couplings(iblok) /= 0_ip ) then
        call COU_CHECK_CONVERGENCE(iblok,kfl_gozon)
        if( kfl_gozon == 1_ip ) kfl_gocou = 1_ip
     else
        kfl_gozon = 0_ip 
     end if
     !
     ! Output convergence
     !
     call cou_cvgunk()   
     
  else
     !
     ! Go to next block
     !
     kfl_gozon = 0_ip
     
  end if

#ifdef COMMDOM  
  call moduls(ITASK_ENDZON) 
#endif 

end subroutine Endzon
