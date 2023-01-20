!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @defgroup Magnet
!> @{
!> @file    Magnet.f90
!> @author  
!> @date    2019-10-09
!> @brief   Magnet
!> @details This module solves the H-formulation of Maxwell's equation
!>	    for the magnetic field using linear Edge Finite Elements
!>	    (EFEM)
!> @} 
!-----------------------------------------------------------------------
subroutine Magnet(order)
  !
  ! IMPORTANT
  ! Please check the tips in the following Magnet files before compiling:
  ! * mag_begste.f90
  !
  use def_master

  implicit none

  integer(ip), intent(in) :: order
  !
  ! OBSOLETE:
  ! Remember to set variable kfl_edge_elements = 1 in kernel
  ! In kernel/kermod/ker_readat.f90
  ! Otherwise edge data will not be computed!
  !
  ! UPDATE: no need to set kfl_edge_elements = 1 in kernel
  ! In case.ker.dat add in NUMERICAL_TREATMENT / MESH the following:
  ! EDGE_ELEMENT:       On
  !
  ! This will set kfl_edge_elements to 1 in kernel/kermod/ker_readat.f90
  !
  select case (order)
    case(ITASK_TURNON)
      ! Launch module
      call mag_turnon()
    case(ITASK_INIUNK)
      ! Launch initial condition
      call mag_iniunk()
    case(ITASK_TIMSTE) 
      ! Time step
      call mag_timste()
    case(ITASK_BEGSTE)
      call mag_begste()
    case(ITASK_DOITER)
      call mag_doiter()
    case(ITASK_CONCOU)
      ! call mag_concou()
    case(ITASK_CONBLK)
      ! call mag_conblk()
    case(ITASK_ENDSTE)
      call mag_endste()
    case(ITASK_OUTPUT)
      call mag_output()
    case(ITASK_TURNOF)
      call mag_turnof()
  end select

end subroutine Magnet
