!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine pts_turnon()
  !-----------------------------------------------------------------------
  !****f* Partis/pts_turnon
  ! NAME 
  !    pts_turnon
  ! DESCRIPTION
  !    This routine performs the following tasks:
  !    - Gets file names and open them.
  !    - Read data for the incompressible NS equations.
  !    - Write some info
  !    - Check errors
  !    - Allocate memory
  ! USES
  !    pts_openfi
  !    pts_reaphy
  !    pts_reabcs
  !    pts_reanut
  !    pts_reaous
  !    pts_outinf
  !    pts_outerr
  !    pts_memall
  ! USED BY
  !    Nastin
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_master
  use def_domain
  use def_partis
  implicit none
  !
  ! Initial variable
  !
  call pts_inivar(0_ip)
  !
  ! Read the physical and numerical problems
  !
  call pts_reaphy()
  call pts_reanut()
  !
  ! Default postprocess variables
  !
  call pts_inivar(3_ip)
  !
  ! Read output and boundary conditions
  !
  call pts_reaous()
  call pts_reabcs()
  !
  ! Parall service
  !
  call pts_parall(1_ip)  
  !
  ! Impose boundary conditions
  !
  call pts_inibcs()
  !
  ! Initial variable
  !
  call pts_inivar(1_ip)
  !
  ! Allocate memory: can only have MLAGR living at the same time
  !
  call pts_memall()
  !
  ! Initialize variables that needed allocation
  !
  call pts_inivar(2_ip)
  !
  ! Open files
  !
  call pts_openfi(1_ip)
  !
  ! Check errors
  !
  call pts_outerr()

end subroutine pts_turnon
