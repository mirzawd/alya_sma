!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    gus_turnon.f90
!> @author  Guillaume Houzeaux
!> @brief   Turn on Nastin module
!> @details Read data and allocate memory
!> @}
!------------------------------------------------------------------------
subroutine gus_turnon()

  use def_kintyp_basic, only : ip
  use def_master,       only : ITASK_TURNON
  
  implicit none
  !
  ! Initial variables
  !
  call gus_inivar(0_ip)
  !
  ! Read the physical problem
  !
  call gus_reaphy()
  !
  ! Read the numerical treatment
  !
  call gus_reanut()
  !
  ! Read the output strategy
  !
  call gus_reaous()
  !
  ! Read the boundary conditions
  !
  call gus_reabcs()
  !
  ! Service: Parall
  !
  call gus_parall()
  !
  ! Impose and modify boundary conditions
  !
  call gus_inibcs()
  call gus_updbcs(ITASK_TURNON)
  !
  ! Initial variables
  !
  call gus_inivar(1_ip)
  !
  ! Allocate memory
  !
  call gus_memall()
!!$  !
!!$  ! Warnings and errors
!!$  !
!!$  call gus_outerr()
!!$  !
!!$  ! Open additional files
!!$  !
!!$  call gus_openfi(2_ip)
  
end subroutine gus_turnon

