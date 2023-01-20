!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_neutro

  use def_master
  use def_neutro
  implicit none
  private

  public :: neutro_main

contains

  subroutine neutro_main(order)
    implicit none
    integer(ip), intent(in) :: order
    external Neutro
    call Neutro(order)
  end subroutine neutro_main

end module mod_neutro
