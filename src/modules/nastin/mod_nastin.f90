!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    mod_nastin.f90
!> @author  houzeaux
!> @date    2020-04-22
!> @brief   Module for Nastin
!> @details Main nastin module with some useful tools
!-----------------------------------------------------------------------

module mod_nastin

  use def_master
  use def_nastin
  implicit none
  private

  public :: nastin_main, nastin_solution_strategy

contains

  subroutine nastin_main(order)
    implicit none
    integer(ip), intent(in) :: order
    call Nastin(order)
  end subroutine nastin_main
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-22
  !> @brief   Numerical strategy
  !> @details Set some parameters determining the solution strategy
  !> 
  !-----------------------------------------------------------------------

  subroutine nastin_solution_strategy()

    if(      kfl_algor_nsi == 1 ) then
       !
       ! Monolithic
       !
       NSI_MONOLITHIC          = .true.
       NSI_SCHUR_COMPLEMENT    = .false.
       NSI_FRACTIONAL_STEP     = .false.
       NSI_SEMI_IMPLICIT       = .false.
       NSI_ASSEMBLY_CONVECTIVE = .true.
       NSI_ASSEMBLY_VISCOUS    = .true.
       
    else if( kfl_algor_nsi == 5 ) then
       !
       ! Schur complement
       !
       NSI_MONOLITHIC          = .false.
       NSI_SCHUR_COMPLEMENT    = .true.
       NSI_FRACTIONAL_STEP     = .false.
       NSI_SEMI_IMPLICIT       = .false.
       NSI_ASSEMBLY_CONVECTIVE = .true.
       NSI_ASSEMBLY_VISCOUS    = .true.
       
    else if( kfl_algor_nsi == 3 ) then
       !
       ! Semi implicit 
       !
       NSI_MONOLITHIC          = .false.
       NSI_SCHUR_COMPLEMENT    = .false.
       NSI_FRACTIONAL_STEP     = .false.
       NSI_SEMI_IMPLICIT       = .true.
       NSI_ASSEMBLY_CONVECTIVE = .false.
       NSI_ASSEMBLY_VISCOUS    = .true.
       
    else if( kfl_algor_nsi == 2 ) then
       !
       ! Fractional step
       !
       NSI_MONOLITHIC          = .false.
       NSI_SCHUR_COMPLEMENT    = .false.
       NSI_FRACTIONAL_STEP     = .true.
       NSI_SEMI_IMPLICIT       = .false.
       NSI_ASSEMBLY_CONVECTIVE = .false.
       NSI_ASSEMBLY_VISCOUS    = .false.
              
    end if

  end subroutine nastin_solution_strategy

end module mod_nastin
!> @}
