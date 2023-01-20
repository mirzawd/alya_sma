!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup kermod
!> @{
!> @file    ker_interp.f90
!> @author  houzeaux
!> @date    2019-06-17
!> @brief   Array redistribution
!> @details Array redistribution and reallocation in case values are
!>          not needed
!> @} 
!-----------------------------------------------------------------------

subroutine ker_interp()
  
  use def_master
  use def_domain
  use def_kermod
  use mod_parall,        only : par_memor
  use mod_ker_proper,    only : ker_proper_allocate_properties
  use mod_ker_arrays,    only : ker_arrays
  use mod_ker_updpro,    only : ker_updpro
  use def_AMR
  use mod_AMR_interpolate
  implicit none
  !
  ! Boundary conditions
  !
  call ker_inibcs()
  call AMR_interpolate(walld,'NPOIN',MEMOR=par_memor,VARIABLE_NAME='WALLD')
  call AMR_interpolate(walln,'NPOIN',MEMOR=par_memor,VARIABLE_NAME='WALLN')
  !
  ! Arrays
  !
  call ker_arrays('INTERPOLATE')
  !
  ! Reallocate properties and recompute them
  !
  call ker_proper_allocate_properties()
  call ker_updpro()

end subroutine ker_interp

