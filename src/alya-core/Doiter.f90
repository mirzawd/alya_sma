!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine Doiter()
!-----------------------------------------------------------------------
!****f* master/Doiter
! NAME
!    Doiter
! DESCRIPTION
!    This routine calls the different problems to be solved within
!    one iteration      
! USES
!    Nastin
!    Temper
!    Codire
!    Alefor
! USED BY
!    Alya
!***
!-----------------------------------------------------------------------
  use def_kintyp,        only : ip
  use def_master,        only : itinn
  use def_master,        only : ittim
  use def_master,        only : ITASK_DOITER
  use mod_ker_detection, only : ker_detection_doiter
  use mod_messages,      only : livinf
  use mod_moduls,        only : moduls 
  use mod_alya2talp,     only : alya2talp_MonitoringRegionStart
  use mod_alya2talp,     only : alya2talp_MonitoringRegionStop
  use mod_reset,         only : reset_check
  implicit none
  
  call livinf(5_ip,' ',0_ip)
  call livinf(6_ip,' ',0_ip)

  itinn(0) = ittim

  call moduls(ITASK_DOITER)  
  call Kermod(ITASK_DOITER)
  !
  ! Check reset
  !
  call reset_check()

end subroutine Doiter
