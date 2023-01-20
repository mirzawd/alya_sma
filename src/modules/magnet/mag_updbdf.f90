!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_updbdf(itask)

  use def_master
  use def_magnet

  implicit none

  integer(ip), intent(in) :: itask

  select case (itask)

  case (ITASK_BEGSTE)
    !
    ! Construct Hpi and Hpe
    !
    if (.not. kfl_reset_mag) then
      !
      ! Correct BDF order:
      ! order s cannot be applied until time step s
      !
      bdfode_mag % sc = min( bdfode_mag % s, timeStep_mag)
      call bdfodeUpdate(bdfode_mag, dt_mag, Hp_mag)
      !
    end if
   
  end select

end subroutine mag_updbdf
