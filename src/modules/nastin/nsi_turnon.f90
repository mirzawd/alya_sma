!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_turnon.f90
!> @author  Guillaume Houzeaux
!> @brief   Turn on Nastin module
!> @details Read data and allocate memory
!> @}
!------------------------------------------------------------------------
subroutine nsi_turnon()

  use def_kintyp
  use def_master
  use def_domain
  use def_nastin

  implicit none
#ifdef OPENACCHHH
  integer(ip)  :: gpunum, ngpus
#endif
  !
  ! Initial variables
  !
  call nsi_inivar(0_ip)
  !
  ! Read the physical problem
  !
  call nsi_reaphy()
  !
  ! Read the numerical treatment
  !
  call nsi_reanut()
  !
  ! Read the output strategy
  !
  call nsi_reaous()
  !
  ! Read the boundary conditions
  !
  call nsi_reabcs()
  !
  ! Service: Parall
  !
  call nsi_sendat()
  !
  ! Modify boundary conditions
  !
  call nsi_inibcs()
  call nsi_updbcs(ITASK_TURNON)
  !
  ! Initial variables
  !
  call nsi_inivar(1_ip)
  !
  ! Allocate memory
  !
  call nsi_memall()
  !
  ! Warnings and errors
  !
  call nsi_outerr()
  !
  ! Open additional files
  !
  call nsi_openfi(2_ip)
  !
  ! For openacc data exchange cpu-gpu
  ! Perhaps add contiguous to pointers
  !
#if defined OPENACCHHH && defined _OPENACC
  if (kfl_paral /= 0 ) then
     if( kfl_savda == 0 ) then
        !$acc enter data copyin(coord,ltype,lnods,lnodb,gravi_nsi)
     else
        !$acc enter data copyin(coord,ltype,lnods,lnodb,gravi_nsi,   &
        !$acc                   elmda_gpvol, elmda_gpcar )
     end if
  end if
#endif

end subroutine nsi_turnon

