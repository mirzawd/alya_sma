!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kinds_and_types
!> @{
!> @file    def_kintyp_functions.g90
!> @author  houzeaux
!> @date    2020-04-04
!> @brief   Functions
!> @details Differents space and time functions
!-----------------------------------------------------------------------

module def_kintyp_functions

  use def_kintyp_basic
  !
  ! Space time function
  !
  type typ_space_time_function
     integer(ip)             :: ndime         ! Number of dimensions
     integer(ip)             :: nexpr         ! Size of the expression
     integer(ip)             :: numfield      ! Field number, when the function is a field
     character(5)            :: name          ! Name
     character(400), pointer :: expression(:) ! Expression to be parsed
  end type typ_space_time_function
  !
  ! Time function
  !
  type typ_time_function
     integer(ip)             :: kfl_type      ! Function type
     integer(ip)             :: npara         ! Number of parameters
     real(rp),       pointer :: parameters(:) ! Parameters
     character(5)            :: name          ! Name
  end type typ_time_function
  !
  ! Windkessel functions
  !
  type typ_windk_system
     integer(ip)                           :: sysid             ! ID of the system
     character(5)                          :: name              ! Name
     integer(ip)                           :: wdks_model        ! Code for the windkessel model
     integer(ip)                           :: nparam            ! Number of parameters
     real(rp), dimension(:), pointer       :: params            ! parameters
     integer(ip)                           :: ID_IN             ! module ID for the input
     integer(ip)                           :: ID_OUT            ! module ID for the output
     integer(ip)                           :: tag_in            ! tag got the input
     integer(ip)                           :: tag_out           ! tag for the output
     integer(ip)                           :: ndxs              ! Number of derivatives
     real(rp), dimension(:), pointer       :: xprev             ! previous converged values of the input
     real(rp), dimension(:), pointer       :: yprev             ! previous converged values of the output
     real(rp), dimension(2)                :: yrelaxed          ! relaxed values of Y for the iteration
     real(rp), dimension(2)                :: yunrelaxed        ! unrelaxed values of Y for the iteration
     real(rp), dimension(2)                :: w                 ! relaxation coefficients computed
     real(rp)                              :: x_in              ! read variable
     real(rp)                              :: y_out             ! computed output
     real(rp)                              :: stored_time_step  ! computed output
     integer(ip)                           :: iflow_nsi         ! tag for the Bazilev BC
     integer(ip)                           :: discrete_function ! Time function to change the parameters, 1-4 time dependent parameters. =0 if no function is associated
  endtype typ_windk_system
  !
  ! Windkessel functions
  !
  type typ_pump_curve
     character(5)                          :: name             ! Name
     integer(ip)                           :: model            ! Code for the pump curve model
     integer(ip)                           :: nparam           ! Number of parameters
     real(rp), dimension(:), pointer       :: params           ! parameters
     character(5)                          :: vhvad            ! variable speed
 end type 
end module def_kintyp_functions
!> @}
