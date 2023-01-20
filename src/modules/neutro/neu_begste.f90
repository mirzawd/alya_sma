!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



 !------------------------------------------------------------------------
!> @addtogroup Neutro 
!> @{
!> @file    neu_begste.f90
!> @date    31/03/2016
!> @author  Guillaume Houzeaux
!> @brief   Begin time step
!> @details begin time step
!> @} 
!------------------------------------------------------------------------

subroutine neu_begste()

  use def_kintyp
  use def_neutro
  implicit none
  external ::  neu_updunk

  ! if(kfl_stead_neu/=1) then     
      !
      ! Initial guess fo the velocity: u(n,0,*) <-- u(n-1,*,*).
      !
      call neu_updunk(1_ip)
  ! endif

end subroutine neu_begste
