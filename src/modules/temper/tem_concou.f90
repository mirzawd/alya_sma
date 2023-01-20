!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_concou()
!-----------------------------------------------------------------------
!****f* Temper/tem_concou
! NAME 
!    tem_concou
! DESCRIPTION
!    This routine checks the temperature convergence of the run.
! USED BY
!    Temper
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_temper
  use mod_outfor, only : outfor
  implicit none
  !
  ! Check convergence
  !
  if(kfl_conve(modul)==1) then
     if(resid_tem>cotol_tem) kfl_gocou = 1
  end if
  glres(modul) = resid_tem
  !
  ! Output residuals
  !
  coutp(1)='TEMPERATURE'
  routp(1)=resid_tem
  call outfor(9_ip,lun_outpu,' ')

  call tem_coupli(ITASK_CONCOU)

end subroutine tem_concou
