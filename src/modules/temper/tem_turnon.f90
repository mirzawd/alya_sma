!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_turnon()
  !-----------------------------------------------------------------------
  !****f* Temper/tem_turnon
  ! NAME 
  !    tem_turnon
  ! DESCRIPTION
  !    This routine performs the following tasks:
  !    - Gets file names and open them.
  !    - Read data for the temperature equation.
  !    - Write some info
  !    - Check if there are errrors
  !    - Allocate memory
  ! USES
  !    tem_openfi
  !    tem_reaphy
  !    tem_reabcs
  !    tem_reanut
  !    tem_reaous
  !    tem_outinf
  !    tem_memall
  !    tem_restar
  ! USED BY
  !    Temper
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_temper
  implicit none
  !
  ! Initial variables
  !
  call tem_inivar(zero)
  !
  ! Read the physical problem
  !
  call tem_reaphy()
  !
  ! Read the numerical treatment
  !
  call tem_reanut()
  !
  ! Read the output strategy
  !
  call tem_reaous()
  !
  ! Read the boundary conditions
  !
  call tem_reabcs()
  !
  ! Parall service
  !
  call tem_sendat()
  !
  ! Initial variables
  !
  call tem_inibcs()
  !
  ! Initial variables
  !
  call tem_inivar(1_ip)
  !
  ! Write info
  !
  call tem_outinf()
  !
  ! Warnings and errors
  !
  call tem_outerr()
  !
  ! Allocate memory
  !
  call tem_memall()
  !
  ! Open additional files
  !
  call tem_openfi(two)

end subroutine tem_turnon
