!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Magnet
!> @{
!> @file    mag_membcs.f90
!> @author  houzeaux
!> @date    2020-04-14
!> @brief   Boundary conditions
!> @details Boundary conditions
!> @} 
!-----------------------------------------------------------------------

subroutine mag_membcs()

  use def_master
  use def_domain
  use def_magnet
  use def_kermod
  use mod_memory
  implicit none
  
  call memory_alloca(mem_modul(1:2,modul), 'KFL_FIXNO_MAG', 'mag_membcs' , kfl_fixno_mag, ndime, meshe(ndivi) % nedge)
  call memory_alloca(mem_modul(1:2,modul), 'BVESS_MAG'    , 'nmag_membcs', bvess_mag    , ndime, meshe(ndivi) % nedge, 1_ip)

end subroutine mag_membcs
