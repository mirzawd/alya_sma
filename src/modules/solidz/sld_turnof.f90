!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @defgroup Solidz
!> @{
!> @file    sld_turnof.f90
!> @date    2006
!> @author  Guillaume Houzeaux
!> @brief   This routine closes the run for the solidz module
!> @details This routine closes the run for the solidz module
!> @}
!------------------------------------------------------------------------

subroutine sld_turnof()

  use def_kintyp, only : ip
  
  implicit none

  external :: sld_openfi
  !
  ! Close files
  !
  call sld_openfi(2_ip)
  
end subroutine sld_turnof

