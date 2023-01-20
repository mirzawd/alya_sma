!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    Domtra.f90
!> @author  houzeaux
!> @date    2019-06-11
!> @brief   Mesh transformation
!> @details Mesh transformation, multiplication, scaling, etc.
!> @} 
!-----------------------------------------------------------------------

subroutine Domtra()

  use def_kintyp,              only : ip,rp
  use def_master,              only : cpu_start
  use def_master,              only : CPU_MESH_MULTIPLICATION
  use def_master,              only : CPU_CONSTRUCT_DOMAIN
  use def_domain,              only : kfl_chege
  use mod_mesh_type,           only : mesh_type_save_original_mesh
  use mod_mesh_type,           only : mesh_type_allocate_initialize
  use mod_mesh_multiplication, only : mesh_multiplication
  use mod_periodicity,         only : periodicity_setup
  implicit none

  real(rp) :: time1,time2

  call cputim(time1) 
  !
  ! Check mesh if required
  !
  call mescek(2_ip)
  !
  ! Transform geometry (rotation, translation, scaling)
  !
  call domvar(1_ip)                                       ! COORD
  !
  ! Compute some derived parameters
  !
  call domvar(2_ip)                                       ! LNUTY, LTYPF, NNODF
  !
  ! Periodicity
  !
  call periodicity_setup()
  !
  ! Allocate mesh type
  !
  call mesh_type_allocate_initialize()                    ! MESHE
  !
  ! Save original mesh juste in case
  !
  call mesh_type_save_original_mesh()                     ! MESHE(0)
  !
  ! Timing
  !
  call cputim(time2) 
  cpu_start(CPU_CONSTRUCT_DOMAIN) = cpu_start(CPU_CONSTRUCT_DOMAIN) + time2 - time1
  !
  ! Mesh multiplication
  !
  call cputim(time1) 
  call mesh_multiplication()
  !
  ! Compute shape functions & derivatives
  !
  call cshder(3_ip)                                       ! SHAPE, DERIV, HESLO, WEIGP...116
  !
  ! Modify mesh for chebyshev interpolation functions
  !
  call chebyshev_coordinates()
  !
  ! Global mesh dimensions
  !
  call par_mesh_dimensions()                              ! NPOIN_TOTAL, NELEM_TOTAL, NBOUN_TOTAL
  
  call cputim(time2) 
  cpu_start(CPU_MESH_MULTIPLICATION) = cpu_start(CPU_MESH_MULTIPLICATION) + time2 - time1

  
end subroutine Domtra
