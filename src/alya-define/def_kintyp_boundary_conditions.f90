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
!> @details Definition of kinds and types for boundary conditions
!>
!-----------------------------------------------------------------------

module def_kintyp_boundary_conditions

  use def_kintyp_basic

  type bccodes
     integer(ip)          :: varna
     real(rp)             :: xcent(3)
     integer(ip)          :: ldofr(5)
     integer(ip)          :: funty
     real(rp)             :: param(6)
  end type bccodes
  !
  ! Node codes
  !
  type bc_nodes1
     integer(ip),     pointer :: lcode(:)
     integer(ip)              :: kfl_fixno
     integer(ip)              :: kfl_value
     real(rp),        pointer :: bvess(:)
     integer(ip)              :: kfl_funtyp ! 0: no function, 1: space_time, 2: time, 3: windkessel, 4: fields
     integer(ip)              :: kfl_funno  ! Old time functions
     integer(ip)              :: kfl_fixrs
     character(5)             :: tag
     character(5)             :: fname
  end type bc_nodes1
  type bc_nodes
     integer(ip)              :: kfl_ibopo
     integer(ip)              :: ndofn
     integer(ip)              :: ncode
     type(bc_nodes1), pointer :: l(:)
  end type bc_nodes
  !
  ! Boundary codes
  !
  type bc_bound1
     integer(ip)              :: lcode
     integer(ip)              :: kfl_fixbo
     integer(ip)              :: kfl_value
     integer(ip)              :: kfl_funtyp ! 0: no function, 1: space_time, 2: time, 3: windkessel, 4: fields, 5: pump
     integer(ip)              :: kfl_funbo
     character(5)             :: tag
     real(rp),        pointer :: bvnat(:)
     character(5)             :: fname
  end type bc_bound1
  type bc_bound
     integer(ip)              :: ndofn
     integer(ip)              :: ncode
     type(bc_bound1), pointer :: l(:)
  end type bc_bound

end module def_kintyp_boundary_conditions
!> @}
