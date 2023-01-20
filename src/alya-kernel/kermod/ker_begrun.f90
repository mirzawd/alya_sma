!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    ker_begrun.f90
!> @date    14/06/2019
!> @author  Guillaume Houzeaux
!> @brief   Beginning the run... 
!> @details Beginning the run... we can compute matrices!
!> @}
!-----------------------------------------------------------------------

subroutine ker_begrun()

  use def_master,        only : ITASK_BEGRUN
  use mod_ker_noslwa,    only : ker_noslwa
  use mod_wall_exchange, only : ker_waexlo
  
  implicit none
  !
  ! Wall exchange strategy
  ! 
  call ker_waexlo()
  !
  ! NO Slip Wall law - wall law adding extra viscosity 
  !
  call ker_noslwa() 
  !
  ! Velocity and temperature field, etc.
  ! Some field may be needed at this point
  !
  call ker_velfun(ITASK_BEGRUN)
  call ker_temfun(ITASK_BEGRUN)
  call ker_confun(ITASK_BEGRUN)
  call ker_disfun(ITASK_BEGRUN)
  call ker_arefun(ITASK_BEGRUN)

end subroutine ker_begrun
