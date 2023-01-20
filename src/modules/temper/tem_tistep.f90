!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_tistep
  !-----------------------------------------------------------------------
  !****f* Temper/tem_tistep
  ! NAME 
  !    tem_tittim
  ! DESCRIPTION
  !    This routine sets the time step
  ! USES
  ! USED BY
  !    tem_begite
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod, only : kfl_adj_prob
  use def_temper
  use mod_ADR,    only : ADR_time_strategy
  use mod_outfor, only : outfor
  implicit none

  !
  ! Actualize time integration parameters: above is obsolete
  !
  if(kfl_adj_prob == 1 ) dtinv = 0.0_rp  
  call ADR_time_strategy(ittim,dtinv,dtinv_old,ADR_tem)

  routp(1) = dtcri_tem
  routp(2) = 0.0_rp
  routp(3) = 0.0_rp
  ioutp(1) = kfl_timei_tem
  ioutp(2) = kfl_stead_tem
  call outfor(8_ip,lun_outpu,' ')

end subroutine tem_tistep
