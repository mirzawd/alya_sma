!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine exm_turnon
  !-----------------------------------------------------------------------
  !****f* Exmedi/exm_turnon
  ! NAME 
  !    exm_turnon
  ! DESCRIPTION
  !    This routine performs the following tasks:
  !    - Gets file names and open them.
  !    - Read some data for the run
  !    - Write some info
  !    - Allocate memory
  ! USES
  !    _openfi
  !    exm_reaphy
  !    exm_reabcs
  !    exm_reanut
  !    exm_reaous
  !    exm_outinf
  !    exm_memall
  !    exm_massma
  ! USED BY
  !    Exmedi
  !***
  !-----------------------------------------------------------------------
  
  use def_parame
  use def_master
  implicit none
  !
  ! Open set files
  !
  call exm_openfi(1_ip)  
  !
  ! Open most of the files
  !
  call exm_inivar(0_ip)
  !
  ! Read the physical problem (except the material fields)
  !
  call exm_reaphy
  !
  ! Read the numerical treatment
  !
  call exm_reanut
  !
  ! Read the output strategy
  !
  call exm_reaous
  !
  ! Read the field of material properties
  !
!!!!  call exm_reamtr
  !
  ! Read the boundary conditions
  !
  call exm_reabcs
  !
  ! Service: Parall
  !
  call exm_sendat(0_ip)
  !
  ! Write info
  !
  call exm_outinf
  !
  ! Dimensions needed for allocation
  ! 
  call exm_inivar(1_ip)
  !
  ! Service: Parall
  !
  call exm_sendat(1_ip)
  !
  ! Allocate info
  !
  call exm_memall()
  !
  ! Compute fiber fields with a Streeter model
  !
  call exm_street
end subroutine exm_turnon
