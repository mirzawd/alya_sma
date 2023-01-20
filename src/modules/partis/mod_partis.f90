!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    mod_partis.f90
!> @author  houzeaux
!> @date    2020-04-02
!> @brief   Partis main module
!> @details Manage memory of partis, etc.
!-----------------------------------------------------------------------

module mod_partis

  use def_master
  use def_domain
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_deallo
  use def_partis

  private

  public :: partis_main
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-02
  !> @brief   IMain subroutine
  !> @details Main subroutine
  !> 
  !-----------------------------------------------------------------------

  subroutine partis_main(order)
    implicit none
    integer(ip), intent(in) :: order
    call Partis(order)
  end subroutine partis_main

end module mod_partis
!> @}
