!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_turnon
  !-----------------------------------------------------------------------
  !****f* Levels/lev_turnon
  ! NAME 
  !    lev_turnon
  ! DESCRIPTION
  !    This routine performs the following tasks:
  !    - Gets file names and open them.
  !    - Read data for the temperature equation.
  !    - Write some info
  !    - Check if there are errrors
  !    - Allocate memory
  ! USES
  !    lev_openfi
  !    lev_reaphy
  !    lev_reabcs
  !    lev_reanut
  !    lev_reaous
  !    lev_outinf
  !    lev_memall
  !    lev_restar
  ! USED BY
  !    Levels
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_levels
  implicit none
  !
  ! Initial variables
  !
  call lev_inivar(1_ip)
  !
  ! Read the physical problem
  !
  call lev_reaphy()
  !
  ! Read the numerical treatment
  !
  call lev_reanut()
  !
  ! Read the output strategy
  !
  call lev_reaous()
  !
  ! Read the boundary conditions
  !
  call lev_reabcs()
  !
  ! Service: Parall
  !
  call lev_parall(1_ip)
  !
  ! Initial variables
  !
  call lev_inibcs()
  !
  ! Initial variables
  !
  call lev_inivar(2_ip)
  !
  ! Write info
  !
  call lev_outinf()
  !
  ! Warnings and errors
  !
  call lev_outerr()
  !
  ! Allocate memory
  !
  call lev_memall(1_ip)
  !
  ! Open additional files
  !
  call lev_openfi(2_ip)
  !
  ! Initial values for fleve
  !
  call lev_iniun0()

end subroutine lev_turnon
