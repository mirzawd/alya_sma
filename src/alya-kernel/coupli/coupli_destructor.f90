!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Coupling
!> @{
!> @file    coupling_destructor.f90
!> @author  houzeaux
!> @date    2019-06-11
!> @brief   Destroy coupling
!> @details Reinitialize all coupling arrays
!> @} 
!-----------------------------------------------------------------------

subroutine coupli_destructor()

  use def_kintyp,          only : ip
  use mod_coupling_memory, only : COU_DEALLOCATE_SINGLE_COUPLING
  use def_coupli,          only : mcoup
  use def_coupli,          only : coupling_type
  implicit none

  integer(ip) :: icoup 

  do icoup = 1,mcoup
     call COU_DEALLOCATE_SINGLE_COUPLING(coupling_type(icoup)) 
  end do
  
end subroutine coupli_destructor
