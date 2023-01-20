!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    def_kermod_generic.f90
!> @author  houzeaux
!> @date    2020-11-20
!> @brief   Definition for properties
!> @details Definitions for generic calculations of properties
!-----------------------------------------------------------------------

module def_vectorization

  use def_kintyp_basic, only : ip
  implicit none
  public
  
  integer(ip) :: ivect

end module def_vectorization
!> @}
