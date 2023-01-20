!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_selfie()
  !-----------------------------------------------------------------------
  !****f* Magnet/mag_selfie
  ! NAME
  !    mag_selfie
  ! DESCRIPTION
  !    This routine sets Self-Field boundary conditions
  ! USES
  !    mag_biosav
  ! USED BY
  !    mag_solite
  !***
  !-----------------------------------------------------------------------

  use def_magnet
  use def_master, only: INOTMASTER, ISEQUEN
  use mod_memory, only: memory_size
  use mod_communications, only: PAR_SUM

  implicit none

  integer(ip) ::    &
    iedge,    &
    jedge,    &
    curproc,    &
    procacc
  !
  ! Contribution of local elements to the whole domain boundary
  !
  if (INOTMASTER) then
    !
    ! Self field computed by Biot-Savart on the domain boundary: Hbs_mag
    !
    call mag_biosav_mod()
    !
    do iedge = 1_ip, memory_size(globoucen_mag, 2_ip)
      !
      Hsfg_mag(iedge) = dot_product(Hbs_mag(:, iedge), globoutan_mag(:, iedge))
      !
    end do
    !
  end if
  !
  ! Contribution of all elements to the local domain boundary
  !
  if (ISEQUEN) then
    !
    Hsfl_mag = Hsfg_mag
    !
  else
    !
    curproc = -1_ip
    procacc = 0_ip
    jedge = 0_ip
    !
    do iedge = 1_ip, memory_size(Hsfg_mag)
      !
      call PAR_SUM(Hsfg_mag(iedge))
      !
      do while (iedge > procacc)
        curproc = curproc + 1_ip
        if (curproc >= nproc_par) call runend('mag_selfie: unexpected error')

        procacc = procacc + proc_nbedg(curproc + 1)
      end do
      !
      if (curproc == iproc_par) then
        jedge = jedge + 1_ip
        Hsfl_mag(jedge) = Hsfg_mag(iedge)
      end if
      !
    end do
    !
    if (jedge /= proc_nbedg(iproc_par+1)) call runend('mag_selfie: unexpected error')
    !
  end if
  !
end subroutine mag_selfie
