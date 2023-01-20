!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



  !------------------------------------------------------------------------
!> @addtogroup Partis
!! @{
!> @name    Partis memory
!! @file    pts_memall.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!! @brief   This routine allocate memory for particles
!! @details Allocate memory and initialize particle type
!> @}
!------------------------------------------------------------------------

subroutine pts_memall()

  use mod_pts_arrays, only : pts_arrays
  use def_partis
  use mod_messages
  implicit none
  integer(8)  :: memorysize_for_particles  ! Type integer and of the system-dependent kind C_SIZE_T (from the ISO_C_BINDING module)
  integer(4)  :: istat
  integer(ip) :: ilagr
  !
  ! Allocate 
  !
  call pts_arrays('ALLOCATE')
  !
  ! Allocate memory: can only have MLAGR living at the same time
  !
  allocate( lagrtyp(mlagr), stat=istat )
  if(istat/=0) then
     memorysize_for_particles = storage_size(lagrtyp,KIND=8)*int(mlagr,8)
     call livinf(-7_ip,"Partis failed to allocate storage for particles (lagrtyp), Number of particles:",mlagr)
     call livinf(-7_ip,"Requested memory in bytes (lagrtyp):",int(memorysize_for_particles,ip))
     call runend("pts_memall: Particle storage allocation failure")
  end if

  do ilagr = 1,mlagr
     call lagrtyp(ilagr) % init()
  end do

end subroutine pts_memall
