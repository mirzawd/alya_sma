!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_tistep()
  !-----------------------------------------------------------------------
  !****f* partis/chm_tistep
  ! NAME
  !    chm_tittim
  ! DESCRIPTION
  !    This routine sets the time step
  ! USES
  ! USED BY
  !    chm_begite
  !***
  !-----------------------------------------------------------------------
  use def_master,      only : momod, modul, dtinv, dtinv_old, ittim, lun_outpu, routp, ioutp
  use def_kintyp,      only : ip, rp
  use def_chemic,      only : dtcri_chm, kfl_timei_chm, nclas_chm, ADR_chm
  use mod_ADR,         only : ADR_time_strategy
  use mod_outfor,      only : outfor
  implicit none
  integer(ip)              :: iclas

  do iclas = 1, nclas_chm
     call ADR_time_strategy(ittim,dtinv,dtinv_old,ADR_chm(iclas))
  end do

  routp(1) = dtcri_chm
  routp(2) = 0.0_rp
  routp(3) = 0.0_rp
  ioutp(1) = kfl_timei_chm
  ioutp(2) = momod(modul) % kfl_stead
  call outfor(8_ip,lun_outpu,' ')

end subroutine chm_tistep
