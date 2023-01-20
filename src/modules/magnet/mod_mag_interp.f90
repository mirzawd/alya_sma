!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_mag_interp

  use def_kintyp, only: ip, rp, lg

  implicit none

  type stcntp
    integer(ip) :: ni
    real(rp), allocatable :: xi(:)
    real(rp), allocatable :: yi(:)
    integer :: status = 0
    logical(lg) :: basc = .false., bdes = .false.
  end type stcntp

contains

  !##############################################################
  function interpReadat(datfil, iounit) result(interp)

    implicit none

    character(20), intent(in) :: datfil

    integer, intent(in) :: iounit

    integer(ip) :: irow

    type (stcntp) :: interp

    real(rp) :: x, y

    integer :: status

    open(unit=iounit, file=datfil, action='read', iostat=status)
    if (status /= 0) then
      call runend('interpReadat: interpolation data file could not be opened')
      return
    end if

    interp % ni = 0

    do
      read(iounit, *, iostat=status) x, y
      if (status /= 0) then
        exit
      end if

      interp % ni = interp % ni + 1_ip

    end do

    close(iounit)


    allocate(interp % xi(interp % ni), interp % yi(interp % ni))


    open(unit=iounit, file=datfil, action='read', iostat=status)
    if (status /= 0) then
      call runend('interpReadat: interpolation data file could not be opened')
      return
    end if

!    irow = 1_ip

    do irow = 1, interp % ni
      read(iounit, *, iostat=status) interp % xi(irow), interp % yi(irow)
      if (status /= 0) then
        exit
      end if

      if (irow > 1_ip) then
        if (interp % xi(irow) == interp % xi(irow-1)) then
          interp % status = -1
          exit
        else if (interp % xi(irow) >= interp % xi(irow-1)) then
          interp % basc = .true.
        else if (interp % xi(irow) <= interp % xi(irow-1)) then
          interp % bdes = .true.
        end if
      end if

      if (interp % basc .and. interp % bdes) then
        interp % status = -2
        exit
      end if

!      irow = irow + 1_ip

    end do

    close(iounit)

    if (interp % status /= 0) then
      call runend('interpReadat: x vector is not sorted')
    end if

  end function interpReadat
  !##############################################################

  !##############################################################
  function interpLinear1(interp, x) result(y)

    implicit none

    type (stcntp), intent(in) :: interp

    real(rp), intent(in) :: x

    real(rp) :: y

    integer(ip) :: ix

    y = 0.0_rp

    if ((interp % basc .and. interp % bdes) .or. .not. (interp % basc .or. interp % bdes)) return

    if (interp % basc .and. x <= interp % xi(1)) then
      y = interp % yi(1)
      return
    else if (interp % basc .and. x >= interp % xi(interp % ni)) then
      y = interp % yi(interp % ni)
      return
    else if (interp % bdes .and. x >= interp % xi(1)) then
      y = interp % yi(1)
      return
    else if (interp % bdes .and. x <= interp % xi(interp % ni)) then
      y = interp % yi(interp % ni)
      return
    end if

    do ix = 1, interp % ni - 1_ip

      if (x >= interp % xi (ix) .and. x <= interp % xi (ix + 1_ip)) then
        y = interp % yi(ix) + (interp % yi(ix + 1_ip) - interp % yi(ix)) / (interp % xi(ix + 1_ip) - interp % xi(ix)) * (x - interp % xi(ix))
        return
      end if

      if (x <= interp % xi (ix) .and. x >= interp % xi (ix + 1_ip)) then
        y = interp % yi(ix) + (interp % yi(ix + 1_ip) - interp % yi(ix)) / (interp % xi(ix + 1_ip) - interp % xi(ix)) * (x - interp % xi(ix))
        return
      end if

    end do

  end function interpLinear1
  !##############################################################
 
end module mod_mag_interp
