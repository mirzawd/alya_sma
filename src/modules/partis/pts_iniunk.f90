!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Partis
!! @{
!> @name    Partis initial solution
!! @file    pts_iniunk.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!! @brief   This routine determines the initial solution
!! @details Initial solution, allocate memory and read restart file
!>          \verbatim
!>          NTYLA_PTS ... Maximum particle type number
!>          DEPOE_PTS ... Number of deposited particle per element
!>          LEDEP_PTS ... .TRUE. if an eement hosts a deposited particle 
!>          DEPOS_PTS ... Smoothed number of particles per node 
!>                        (computed from DEPOE_PTS in pts_endite)
!>          \endverbatim
!> @} 
!------------------------------------------------------------------------

subroutine pts_iniunk()
  use def_parame
  use def_domain
  use def_master 
  use def_kermod
  use def_partis
  use mod_memory
  use mod_elmgeo, only : elmgeo_element_characteristic_length
  use mod_messages
  implicit none
  
end subroutine pts_iniunk
