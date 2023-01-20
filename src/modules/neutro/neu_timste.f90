!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Neutro 
!> @{
!> @file    neu_timste.f90
!> @author  Guillaume Houzeaux
!> @brief   Initialize time step
!> @details Initialize time step
!> @} 
!------------------------------------------------------------------------

subroutine neu_timste
  use def_parame
  use def_master
  use def_domain
  use def_neutro
  implicit none
  !
  ! Actualizes safety factor
  ! ittim is last time step number
  ! Do not modify safet at first time step
  ! 
  !if( ittim > 0 ) safet_nsi = min(safeo_nsi*((safex_nsi)**ittim), safma_nsi)
  ! 
  ! Time step size
  !
  !if( kfl_stead_nsi /= 1 ) call neu_updtss()

end subroutine neu_timste
