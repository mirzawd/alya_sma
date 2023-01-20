!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Magnet
!> @{
!> @file    neu_output.f90
!> @date    01/04/2016
!> @author  Guillaume Houzeaux
!> @brief   Finalize coupling iteration
!> @details Check convergence 
!> @} 
!-----------------------------------------------------------------------

subroutine mag_output()

  use def_kintyp,             only: ip
  use def_master,             only: nvarp
  use def_magnet,             only : postce_mag, postev_mag
  use mod_output_postprocess, only : output_postprocess_variables

  implicit none

  external :: mag_outvar
  !
  ! Default
  ! postev_mag = .false.
  ! postce_mag = .false.
  !
  if (postce_mag) then
    !
    ! Evaluate fields at element centroids
    !
    call mag_edgcen()
    !
  else
    !
    ! Integrate fields over elements
    ! Default
    !
    call mag_edgelm()
    !
  end if
  !
  if (postev_mag) then
    !
    ! GAUSS on
    ! Compute Fields at Gauss points in order to:
    ! Extrapolate from element integration points to element nodes
    !
    call mag_gpedge()
    !
  end if
  !
  ! Mesh dependent post-processing
  !
  call output_postprocess_variables(mag_outvar)
  !
end subroutine mag_output
