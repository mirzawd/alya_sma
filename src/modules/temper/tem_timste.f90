!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_timste()
  !-----------------------------------------------------------------------
  !****f* Temper/tem_timste
  ! NAME 
  !    tem_begste
  ! DESCRIPTION
  !    This routine prepares for a new time step of the temperature
  !    equation      
  ! USES
  !    tem_iniunk
  !    tem_updtss
  !    tem_updbcs
  !    tem_updunk
  !    tem_radvuf
  ! USED BY
  !    Temper
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  implicit none
  !
  ! Time step size 
  !
  if( kfl_stead_tem /= 1 ) call tem_updtss()
  
end subroutine tem_timste

