!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @defgroup Maths_Toolbox
!> Toolbox for mathematical functions and subroutines
!> @{
!> @name    ToolBox for mathematics operations
!> @file    def_maths.f90
!> @author  Guillaume Houzeaux
!> @brief   Variables
!> @details Variables fot mod_maths.f90
!
!-----------------------------------------------------------------------
module def_maths
  
  use def_kintyp_basic,      only : ip,rp,lg,i1p,qp
  use mod_memory_basic,      only : memory_alloca
  use mod_memory_basic,      only : memory_deallo
  use mod_memory_basic,      only : memory_resize
  use mod_memory_basic,      only : memory_size
  use mod_memory_basic,      only : memory_copy
  use mod_optional_argument, only : optional_argument
  use mod_std                                  ! defintion of qp and count
  implicit none
  
  integer(8)              :: memor(2)
  real(rp),    parameter  :: epsil = epsilon(1.0_rp)
  real(rp),    parameter  :: pi    = 3.141592653589793238462643383279502884197_rp
  
end module def_maths
!> @}

