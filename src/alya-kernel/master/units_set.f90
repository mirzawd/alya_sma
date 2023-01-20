!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Initia
!> @{
!> @file    units_set.f90
!> @author  Constantine Butakoff
!> @brief   initialize global file handles
!> @details initialize global file handles
!> @}
!-----------------------------------------------------------------------
subroutine units_set()
  use def_master
  use def_domain
  use mod_memory, only : lun_varcount, lun_memor
  use def_coupli, only : lun_coupl_res, lun_coupl_dat, lun_coupl_cvg

  implicit none

  !
  ! Units: General
  !
  lun_pdata     = 11             ! Data file unitm
  lun_outpu     = 12             ! Output (log) file unit
  lun_memor     = 13             ! Memory file unit
  lun_conve     = 14             ! Convergence file unit
  lun_rstar     = 17             ! Restart file unit
  lun_latex     = 18             ! Latex file unit: text
  lun_gnupl     = 19             ! Latex file unit: gnuplot
  lun_commu     = 20             ! Communication with Alya
  lun_binar     = 24             ! Geometry binary file
  lun_pdata_dom = 21             ! Domain data file unit
  lun_elsta_dom = 25             ! Elsest statistics
  lun_elmsh_dom = 26             ! Elsest mesh
  lun_elres_dom = 27             ! Elsest results
  lun_syste     = 28             ! System info
  lun_tempo     = 29             ! Temporary unit
  lun_rstib     = 45             ! IB Restart file unit
  lun_detec     = 33             ! Automatic detection file
  lun_timeline  = 42             ! Timeline file
  lun_memory    = 34             ! Memory evolution
  lun_varcount  = 50             ! Memory variable counter file unit
  lun_state     = 53             ! State file unit
  !  lun_time      = 50             ! Time values file unit
  !
  ! Units: Domain postprocess
  !
  lun_outpu_dom = 22             ! Output domain file unit
  lun_postp     = 15             ! Postprocess domain unit
  lun_posvx     = 16             ! Postprocess voxel unit
  lun_pos00     = 30             ! Additional output file (VU)
  lun_pos01     = 31             ! Additional output file (VU)
  lun_pos02     = 32             ! Additional output file (VU)
  lun_mesh      = 51             ! Unit for mesh output
  lun_perf      = 52             ! Unit for performance
  !
  ! Units 54-55-56-57 RESERVED BY TUBES !
  !
  !kfl_oumes     =  0             ! Do not output
  !kfl_oumes(1)  =  1             ! Output mesh
  !
  ! Units: Set postprocess
  !
  lun_pos09     = 43             ! Additional output file (VU) (Filter)
  lun_pos10     = 44             ! Additional output file (VU) (Filter)
  !
  ! IB
  !
  lun_mshib     = 46             ! IB mesh
  lun_resib     = 47             ! IB results
  lun_mshi2     = 48             ! IB mesh (2)
  lun_resi2     = 49             ! IB results (2)
  !
  ! Lagrangian particles
  !
  lun_rstla     = 39
  lun_posla     = 40
  lun_cvgla     = 41
  !
  ! Coupling
  !
  lun_coupl_dat = 36
  lun_coupl_res = 37
  lun_coupl_cvg = 38
  !
  ! Read/write units
  !
  lun_pdata_dom = 21             ! Domain data file unit
  lun_outpu_dom = 22             ! Output domain file unit
end subroutine