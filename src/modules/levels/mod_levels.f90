!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_levels

  use def_master
  use def_levels
  implicit none
  private

  public :: levels_main

contains

  subroutine levels_main(order)
    implicit none
    integer(ip), intent(in) :: order
    call Levels(order)
  end subroutine levels_main

end module mod_levels
