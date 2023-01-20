!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Alefor
!> @{
!> @file    ale_memall.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966
!> @brief   Arrays allocation subroutine
!> @details Arrays allocation subroutine
!> @} 
!-----------------------------------------------------------------------
subroutine ale_memall()
  use      mod_ale_arrays
  implicit none

  call ale_arrays('ALLOCATE')

end subroutine ale_memall
      
