!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @name    Partis inner iteration
!> @file    pts_doiter.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!> @brief   This routine solves a time step
!> @details This routine solves a time step
!>
!>          call pts_begite
!>          call pts_transport_particles
!>             do communication loops
!>                do particles
!>                  call pts_transport_single_particle
!>                end do
!>                call pts_parallelization_migration
!>             end do
!>          call pts_endite
!> @} 
!------------------------------------------------------------------------

subroutine pts_doiter()
  use def_partis
  use mod_pts_transport
  use def_master,  only : ITASK_ENDITE,ITASK_BEGITE
  use mod_timings, only : timings_ini
  use mod_timings, only : timings_end
  implicit none

  call timings_ini()
  call pts_begite()
  call timings_end(ITASK_BEGITE)
  
  call pts_transport_particles()
  
  call timings_ini()
  call pts_endite()
  call timings_end(ITASK_ENDITE)

end subroutine pts_doiter
