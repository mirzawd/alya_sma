!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    gus_memall.f90
!> @date    06/06/1966
!> @author  Guillaume Houzeaux
!> @brief   Allocate memory
!> @details Allocate memory. Here are allocated variables which should
!>          redistributed. For example, main variables like the velocity
!>          but also postprcess variables like averaging.
!> @}
!-----------------------------------------------------------------------
subroutine gus_memall()

  use mod_gus_arrays,   only : gus_arrays
  
  implicit none

  call gus_arrays('ALLOCATE')

end subroutine gus_memall
 
