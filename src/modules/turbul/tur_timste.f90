!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_timste()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_timste
  ! NAME 
  !    tur_timste
  ! DESCRIPTION
  !    This routine computes the time step
  ! USES
  !    tur_iniunk
  !    tur_updtss
  !    tur_updbcs
  !    tur_updunk
  !    tur_radvuf
  ! USED BY
  !    Temper
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_turbul
  implicit none
  !
  ! Actualizes safety factor
  !
  if(ittim > 0)  safet_tur = min(safeo_tur*((safex_tur)**ittim), safma_tur)
  !
  ! Time step size
  !
  if(kfl_stead_tur/=1) call tur_updtss()

end subroutine tur_timste

