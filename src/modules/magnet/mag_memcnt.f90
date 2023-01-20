!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_memcnt()

  use def_magnet

  implicit none

  integer(ip) :: &
    iconstr

  if (constr_total > 0) then
    allocate(constr_mag(constr_total))

    do iconstr = 1, constr_total
      nullify(constr_mag(iconstr) % cnstr)
      nullify(constr_mag(iconstr) % unkno)
    end do

    allocate(cntmat_mag(constr_total, constr_total))
    allocate(cntrhs_mag(constr_total), cntunk_mag(constr_total))
  end if

  allocate(refres_mag(1_ip + constr_total))
  allocate(residu_mag(1_ip + constr_total))
  allocate(refval_mag(1_ip + constr_total))
  allocate(delval_mag(1_ip + constr_total))

end subroutine mag_memcnt
