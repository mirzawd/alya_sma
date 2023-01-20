!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Partis
!! @{
!> @name    Partis end of a time step
!! @file    pts_endste.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!! @brief   This routine ends a time step
!! @details Restart file
!> @} 
!------------------------------------------------------------------------

subroutine pts_endste()
  use def_master
  use def_kermod
  use def_partis
  use def_domain
  implicit none
  !
  ! Go on in time
  !
  kfl_gotim = 1
  
end subroutine pts_endste

