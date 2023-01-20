!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_turnon.f90
!> @author  Solidz
!> @date
!> @brief   Turn on Solidz module
!> @details
!>          \verbatim
!>          This routine performs the following tasks:
!>           - Initialize variables
!>           - Read data for the solid mechanics.
!>           - Allocate memory
!>           - Write some info
!>           - Check if there are warnings or errors
!>          \endverbatim
!>
!> @}
!------------------------------------------------------------------------

subroutine sld_turnon()

  use def_kintyp,   only : ip
  use def_master,   only : ITASK_TURNON
  use def_solidz,   only : kfl_rigid_sld
  use mod_sld_rbo,  only : sld_rbo_updbcs
  use mod_sld_fe2,  only : fe2_micropp_create

  implicit none
  !
  ! Initial variables
  !
  call sld_inivar(0_ip)
  !
  ! Read the physical problem
  !
  call sld_reaphy()
  !
  ! Read the numerical treatment
  !
  call sld_reanut()
  !
  ! Read the output strategy
  !
  call sld_reaous()
  !
  ! Read the boundary conditions
  !
  call sld_reabcs()
  !
  ! Service: Parall
  !
  call sld_sendat()  
  !
  ! Initialize boundary conditions
  !
  call sld_inibcs()
  if( kfl_rigid_sld == 0 ) then
     call sld_updbcs(ITASK_TURNON)
  else
     call sld_rbo_updbcs(ITASK_TURNON)
  end if
  !
  ! Initial variables
  !
  call sld_inivar(1_ip)
  !
  ! Allocate memory
  !
  call sld_memall()
  !
  ! Warnings and Errors
  !
  call sld_outerr()
  !
  ! Open additional files
  !
  call sld_openfi(1_ip)
  call sld_outinf(1_ip)
  !
  ! Activate MicroPP for FE2
  !
  call fe2_micropp_create( )
  !
  ! Initialize sysnet (if needed)
  !
  call sld_inivar(2_ip)

end subroutine sld_turnon
