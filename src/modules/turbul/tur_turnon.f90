!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_turnon()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_turnon
  ! NAME 
  !    tur_turnon
  ! DESCRIPTION
  !    This routine performs the following tasks:
  !    - Gets file names and open them.
  !    - Read data for the temperature equation.
  !    - Write some info
  !    - Check if there are errrors
  !    - Allocate memory
  ! USES
  !    tur_reaphy
  !    tur_reabcs
  !    tur_reanut
  !    tur_reaous
  !    tur_outinf
  !    tur_memall
  !    tur_restar
  ! USED BY
  !    Turbul
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_turbul
  use mod_iofile
  implicit none
  !
  ! Initial variables
  !
  call tur_inivar(0_ip)
  !
  ! Read the physical problem
  !
  call tur_reaphy()
  !
  ! Read the numerical treatment
  !
  call tur_reanut()
  !
  ! Read the output strategy
  !
  call tur_reaous()
  !
  ! Parall service
  !
  call tur_parall(1_ip)
  !
  ! Initial variables
  !
  call tur_inivar(3_ip)
  !
  ! Read the boundary conditions
  !
  call tur_reabcs()
  !
  ! Parall service
  !
  call tur_parall(2_ip)
  !
  ! Initial conditions
  !
  call tur_inibcs()
  !
  ! Initial variables
  !
  call tur_inivar(2_ip)
  !
  ! Compute additional arrays needed by the model
  !    
  call tur_addarr()
  !
  ! Write info
  !
  call tur_outinf(1_ip)
  !
  ! Warnings and errors
  !
  call tur_outerr()
  !
  ! Initial variables
  !
  call tur_inivar(1_ip)
  !
  ! Allocate memory
  !
  call tur_memall()

end subroutine tur_turnon
