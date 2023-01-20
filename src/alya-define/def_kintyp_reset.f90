!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kinds_and_types
!> @{
!> @file    def_kintyp_reset.f90
!> @author  houzeaux
!> @date    2020-04-04
!> @brief   Reset
!> @details Reset common criterion for module
!-----------------------------------------------------------------------

module def_kintyp_reset

  use def_kintyp_basic, only : ip,rp
  implicit none
  private
  
  integer(ip), parameter :: max_reset                  = 5_ip
  integer(ip), parameter :: RESET_NULL                 = 0_ip
  integer(ip), parameter :: RESET_SOLVER_NOT_CONVERGED = 1_ip
  integer(ip), parameter :: RESET_INNER_NOT_CONVERGED  = 2_ip
  integer(ip), parameter :: RESET_TIME                 = 3_ip
  
  type restyp
     integer(ip)                   :: num_criteria
     integer(ip)                   :: kfl_criterion(max_reset)
     real(rp)                      :: param(10,max_reset)
   contains
     procedure,        pass        :: init => init_reset
  end type restyp

  public :: max_reset
  public :: restyp
  public :: RESET_NULL                 
  public :: RESET_SOLVER_NOT_CONVERGED 
  public :: RESET_INNER_NOT_CONVERGED  
  public :: RESET_TIME
  
contains

  subroutine init_reset(self)
    class(restyp), intent(out) :: self

    self % num_criteria  = 0
    self % kfl_criterion = RESET_NULL
    self % param         = 0.0_rp
    
  end subroutine init_reset
        
end module def_kintyp_reset
!> @}
