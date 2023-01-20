!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_turnof()
  !-----------------------------------------------------------------------
  !****f* Temper/tem_turnof
  ! NAME 
  !    tem_turnof
  ! DESCRIPTION
  !    This routine closes the run for the temperature equation
  ! USES
  !    tem_outcpu
  !    tem_output
  ! USED BY
  !    Temper
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_solver
  use def_temper
  implicit none
  !-------------------------------------------------------------------
  !
  ! Interpolation from coarse to fine mesh
  !
  !-------------------------------------------------------------------
  !
  if(kfl_meshi_tem /= 0_ip) call tem_coarfine(2_ip)
 
end subroutine tem_turnof

