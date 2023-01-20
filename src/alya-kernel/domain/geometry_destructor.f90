!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    geometry_destructor.f90
!> @author  houzeaux
!> @date    2019-06-11
!> @brief   Destroy domain
!> @details Reinitialize all domain arrays
!> @} 
!-----------------------------------------------------------------------

subroutine geometry_destructor()

  use def_kintyp
  use def_master
  use mod_domain
  implicit none

  call domain_memory_deallocate('GEOMETRY' )
  call domain_memory_deallocate('LGROU_DOM')
  call domain_memory_deallocate('LNINV_LOC')
  call domain_memory_deallocate('LEINV_LOC')
  call domain_memory_deallocate('LBINV_LOC')
  call domain_memory_deallocate('LESET'    )
  call domain_memory_deallocate('LBSET'    )
  call domain_memory_deallocate('LNSET'    )
  call domain_memory_deallocate('KFL_CODNO')
  call domain_memory_deallocate('KFL_CODBO')
  call domain_memory_deallocate('XFIEL % A')

end subroutine geometry_destructor
