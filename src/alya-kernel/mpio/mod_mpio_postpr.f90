!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup mpio
!> @{
!> @file    mod_mpio_postpr.f90
!> @author  Damien Dosimont
!> @date    27/09/2017
!> @brief   MPI-IO post process
!> @details This module is a bridge that redirects to the sequential or the parallel
!>          versions of the post process parallel I/O operations
!-----------------------------------------------------------------------

module mod_mpio_postpr

  use def_master
  use mod_mpio_seq_postpr
  use mod_mpio_par_postpr

  implicit none

  private

  public :: posmpio_real_v, posmpio_real_m, posmpio_int_v, posmpio_int_m
  
contains

  subroutine posmpio_real_v()
    if (ISEQUEN) then
       call seq_posmpio_real_v()
    else
       call par_posmpio_real_v()
    end if
  end subroutine posmpio_real_v

  subroutine posmpio_real_m()
    if (ISEQUEN) then
       call seq_posmpio_real_m()
    else
       call par_posmpio_real_m()
    end if
  end subroutine posmpio_real_m

  subroutine posmpio_int_v()
    if (ISEQUEN) then
       call seq_posmpio_int_v()
    else
       call par_posmpio_int_v()
    end if
  end subroutine posmpio_int_v

  subroutine posmpio_int_m()
    if (ISEQUEN) then
       call seq_posmpio_int_m()
    else
       call par_posmpio_int_m()
    end if
  end subroutine posmpio_int_m

end module mod_mpio_postpr
!> @}
