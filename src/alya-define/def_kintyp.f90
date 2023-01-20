!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @defgroup Kinds_and_types
!> Kinds ands types of Alya
!> @{
!> @file    def_kintyp.f90
!> @author  houzeaux
!> @date    2018-12-28
!> @brief   Definition of kinds and types.
!> @details Definition of kinds and types.
!>
!-----------------------------------------------------------------------

module def_kintyp
  
  use def_kintyp_basic               ! Symbolc names for integers, reals and logicals
  use def_kintyp_boundary_conditions ! Boundary conditions
  use def_kintyp_functions           ! Functions
  use def_kintyp_solvers             ! Solvers
  use def_kintyp_postprocess         ! Postprocess
  use def_kintyp_domain              ! Mesh, finite element
  use def_kintyp_module              ! Module
  use def_kintyp_physics             ! Physics
  
end module def_kintyp
!> @}
