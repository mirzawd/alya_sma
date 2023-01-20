!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Magnet
!> @{
!> @file    mod_magnet.f90
!> @author  houzeaux
!> @date    2020-04-22
!> @brief   Module for Magnet
!> @details Main magnet module with some useful tools
!-----------------------------------------------------------------------

module mod_magnet

  use def_master
  use def_magnet
  implicit none
  private

  public :: magnet_main

contains

  subroutine magnet_main(order)
    implicit none
    integer(ip), intent(in) :: order
    call Magnet(order)
  end subroutine magnet_main

end module mod_magnet
!> @}
