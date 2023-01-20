!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    gus_timste.f90
!> @author  houzeaux
!> @date    2020-10-21
!> @brief   Compute time step
!> @details Compute time step
!> @} 
!-----------------------------------------------------------------------

subroutine gus_timste()

  use def_parame
  use def_master
  use def_domain
  use def_gusano
  implicit none
  ! 
  ! Time step size
  !
  if( kfl_stead_gus /= 1 ) call gus_updtss()

end subroutine gus_timste

