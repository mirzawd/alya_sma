!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Magnet
!> @{
!> @file    mag_begste.f90
!> @author  houzeaux
!> @date    2020-04-14
!> @brief   Begin a time step
!> @details Begin a time step
!> @} 
!-----------------------------------------------------------------------

subroutine mag_begste()
  !
  ! In Sources/alya/core/Begste.f90 replace the existing reset code with:
  !
  ! if (kfl_reset == 1) then
  !    call iniste(2_ip)
  !    cutim  = oltim !cutim - dtime
  !    call setgts(ITASK_TIMSTE)
  !    call livinf(201_ip, ' ',1_ip)
  !    kfl_reset = -1
  ! endif
  !
  use def_master, only: ITASK_BEGSTE
  use def_magnet, only: He_mag, Hp_mag

  implicit none
  !
  ! Store previous time step
  !
  Hp_mag = He_mag
  !
  call mag_updbdf(ITASK_BEGSTE)
  !
  call mag_updcon(ITASK_BEGSTE)
  !
  call mag_updbcs(ITASK_BEGSTE)
  !  
end subroutine mag_begste
