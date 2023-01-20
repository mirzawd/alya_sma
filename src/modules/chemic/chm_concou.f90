!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_concou()
  !-----------------------------------------------------------------------
  !****f* partis/chm_concou
  ! NAME
  !    chm_concou
  ! DESCRIPTION
  !    This routine checks the ADS convergence of the run.
  ! USED BY
  !    Temper
  !***
  !-----------------------------------------------------------------------
!  use def_parame
  use def_master, only : kfl_gocou, lun_outpu, modul, coutp, routp, kfl_conve, glres
  use def_chemic, only : cotol_chm, resid_chm
  use def_kintyp, only : ip
  use mod_outfor, only : outfor
  implicit none
  !
  ! Check convergence
  !
  if(kfl_conve(modul)==1) then
     if(resid_chm>cotol_chm) kfl_gocou = 1
  end if
  glres(modul) = resid_chm
  !
  ! Output residuals
  !
  coutp(1)='CONCENTRATION'
  routp(1)=resid_chm
  call outfor(9_ip,lun_outpu,' ')

end subroutine chm_concou
