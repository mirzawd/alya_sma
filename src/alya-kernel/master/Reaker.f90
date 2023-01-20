!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!> @addtogroup Turnon
!> @{
!> @file    Turnon.f90
!> @author  Guillaume Houzeaux
!> @brief   Turnon the run
!> @details Read Kermod to define some important parameters before the
!>          construction of the domain
!>
!> @}
!-----------------------------------------------------------------------
subroutine Reaker()

  use def_kintyp, only : ip
  use def_master, only : ITASK_TURNON
  use def_master, only : kfl_check_data_file

  implicit none
  
  !----------------------------------------------------------------------
  !
  ! Export mesh
  !
  !----------------------------------------------------------------------

  call mpio_export_domain()

  !----------------------------------------------------------------------
  !
  ! Read kermod data... AND
  ! do some calculations required by Kermod= NMATE, NSUBD
  !
  !----------------------------------------------------------------------

  call domvar(0_ip)

  call Kermod(-ITASK_TURNON)
  if( kfl_check_data_file == 1) call runend('O.K.!')
  !
  ! Turn on coupling
  !
  call cou_turnon()
  
  !----------------------------------------------------------------------
  !
  ! Open domain files
  !
  !----------------------------------------------------------------------

  call opfdom(4_ip)

end subroutine Reaker
