!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_day_of_week
  
  implicit none

  public :: day_of_week                   ! Day of the week
contains

  function day_of_week(d, m, y)

    integer(8), intent(in) :: d
    integer(8), intent(in) :: m
    integer(8), intent(in) :: y
    character(9) :: day_of_week
    integer(8) :: day_of_week_i
    integer(8), dimension(12) :: t = (/ 0_8, 3_8, 2_8, 5_8, 0_8, 3_8, 5_8, 1_8, 4_8, 6_8, 2_8, 4_8 /)
    integer(8) :: y_new

    y_new = y
    if (m < 3) y_new = y_new - 1
    day_of_week_i = mod(y_new + y_new/4_8 - y_new*3_8/400_8 - y_new/4000_8 + t(m) + d,7_8)
    if (day_of_week_i == 0) day_of_week_i = 7

    select case (day_of_week_i)
    case (1_8) ; day_of_week = 'Monday'
    case (2_8) ; day_of_week = 'Tuesday'
    case (3_8) ; day_of_week = 'Wednesday'
    case (4_8) ; day_of_week = 'Thursday'
    case (5_8) ; day_of_week = 'Friday'
    case (6_8) ; day_of_week = 'Saturday'
    case (7_8) ; day_of_week = 'Sunday'
    case default ; day_of_week = 'Unknown'
    end select

  end function day_of_week
      
end module mod_day_of_week
