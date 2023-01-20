!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Coupling
!> @{
!> @file    Begzon.f90
!> @author  Guillaume Houzeaux
!> @date    30/09/2014
!> @brief   Start a zonal coupling
!> @details Initialize some variables when zonal coupling is activated
!>          for current block IBLOCK
!> @} 
!-----------------------------------------------------------------------

subroutine Begzon()
  
  use def_kintyp,    only : ip
  use def_master,    only : iblok
  use def_master,    only : mmodu
  use def_master,    only : lmord
  use def_master,    only : itinn
  use def_master,    only : ITASK_BEGZON
  use def_coupli,    only : mcoup
  use def_coupli,    only : coupling_driver_iteration
  use def_coupli,    only : coupling_driver_number_couplings
  use mod_messages,  only : livinf
  use mod_moduls,    only : moduls 
  implicit none
  integer(ip) :: iorde,imodu

  if( mcoup > 0 ) then

     if( coupling_driver_number_couplings(iblok) /= 0 .and. coupling_driver_iteration(iblok) == 0 ) then
        call livinf(-6_ip,'ZONAL COUPLING FOR BLOCK ',iblok)
     end if
     !
     ! Put inner iterations to zero
     !
     do iorde = 1,mmodu
        imodu = lmord(iorde,iblok)
        itinn(imodu) = 0
     end do
     !
     ! Iteration counter
     !
     coupling_driver_iteration(iblok) = coupling_driver_iteration(iblok) + 1
     !
     ! Call modules to exchange data
     !
     call moduls(ITASK_BEGZON)

  end if

#ifdef COMMDOM
     call moduls(ITASK_BEGZON)
#endif 

end subroutine Begzon
