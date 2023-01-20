!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    pts_compute_permutation.f90
!> @author  houzeaux
!> @date    2020-02-05
!> @brief   Permutaiton
!> @details Recompute permutation whenever particles have been created
!>          of have been removed
!> @} 
!-----------------------------------------------------------------------

subroutine pts_compute_permutation()

  use def_kintyp_basic, only : ip 
  use def_master,       only : mem_modul
  use def_master,       only : modul
  use def_master,       only : INOTMASTER
  use mod_memory,       only : memory_alloca
  use mod_memory,       only : memory_deallo
  use mod_memory,       only : memory_size
  use def_partis

  implicit none
  integer(ip) :: ilagr
  
  if( INOTMASTER .and. memory_size(permu_nlagr_pts) < mlagr ) then
     call memory_deallo(mem_modul(1:2,modul),'PERMU_NLAGR','pts_solite',permu_nlagr_pts)
     call memory_alloca(mem_modul(1:2,modul),'PERMU_NLAGR','pts_solite',permu_nlagr_pts,mlagr)
  end if
    
  nlagr_local_pts         = 0
  nlagr_free_pts          = 0
  nlagr_non_migrating_pts = 0

  do ilagr = 1,mlagr

     if(      lagrtyp(ilagr) % kfl_exist ==  0 ) then
        !
        ! Free places
        !
        nlagr_free_pts = nlagr_free_pts + 1

     else if( lagrtyp(ilagr) % kfl_exist == PTS_PARTICLE_EXISTS ) then
        !
        ! Existing particle
        !
        nlagr_local_pts = nlagr_local_pts + 1
        permu_nlagr_pts(nlagr_local_pts) = ilagr

     else if( lagrtyp(ilagr) % kfl_exist <= -2 .and. lagrtyp(ilagr) % kfl_exist /= PTS_PARTICLE_MOVING_MESH ) then
        !
        ! Particles migrating
        !
        nlagr_non_migrating_pts = nlagr_non_migrating_pts + 1
     end if
  end do
  
end subroutine pts_compute_permutation
