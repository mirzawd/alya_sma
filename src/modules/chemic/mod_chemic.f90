!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Chemic
!> @{
!> @file    mod_chemic.f90
!> @author  aboth
!> @date    2020-05-27
!> @brief   Module for Chemic
!> @details Main chemic module with some useful tools
!-----------------------------------------------------------------------

module mod_chemic

  use def_kintyp, only : ip
  implicit none
  private

  public :: chemic_main

contains

  subroutine chemic_main(order)
    implicit none
    integer(ip), intent(in) :: order
    external                :: Chemic

    call Chemic(order)
  end subroutine chemic_main


end module mod_chemic
!> @}
