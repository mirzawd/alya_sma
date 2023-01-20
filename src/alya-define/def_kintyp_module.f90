!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kinds_and_types
!> @{
!> @file    def_kintyp_module.g90
!> @author  houzeaux
!> @date    2020-04-04
!> @brief   Functions
!> @details Communications
!-----------------------------------------------------------------------

module def_kintyp_module

  use def_kintyp_basic
  use def_kintyp_boundary_conditions 
  use def_kintyp_solvers
  use def_kintyp_postprocess
  use def_kintyp_reset
  use def_kintyp_performance

  !----------------------------------------------------------------------
  !
  ! Module
  !
  !----------------------------------------------------------------------

  type tymod
     type(typos),         pointer :: postp(:)           ! Postprocess
     type(soltyp),        pointer :: solve(:)           ! Algebraic solver
     type(soltyp),        pointer :: solad(:)           ! Adjoint algebraic solver
     type(eigtyp),        pointer :: eigen(:)           ! Eigenvalue solver
     type(bc_nodes),      pointer :: tncod(:)           ! Node code type
     type(bc_nodes),      pointer :: tgcod(:)           ! Geometrical node code type
     type(bc_bound),      pointer :: tbcod(:)           ! Boundary code type
     type(restyp),        pointer :: reset              ! Reset type
     type(perf),          pointer :: times(:)           ! Timings
     integer(ip)                  :: nvarn_bcs          ! Number of bc variables (nodes)
     integer(ip)                  :: nvarg_bcs          ! Number of bc variables (geometrical)
     integer(ip)                  :: nvarb_bcs          ! Number of bc variables (boundaries)
     integer(ip)                  :: kfl_modul          ! Existence of module
     integer(ip)                  :: kfl_delay          ! Delay module
     integer(ip)                  :: kfl_conve          ! Cconvergence required
     integer(ip)                  :: kfl_solve          ! When to solve module
     integer(ip)                  :: kfl_goite          ! Keep on iterating
     integer(ip)                  :: kfl_timei          ! Problem is transient
     integer(ip)                  :: kfl_stead          ! Problem is steady
     integer(ip)                  :: ndela              ! Steps to delay module
     integer(ip)                  :: itinn              ! Module inner iteration
     integer(ip)                  :: ittot              ! Total number of iteration
     integer(ip)                  :: miinn              ! Total number of inner iteration
     integer(8)                   :: mem_modul(2)       ! Module memory
     real(rp)                     :: cpu_modul(30)      ! Module CPU time
     real(rp)                     :: dtcri              ! Module critical time
     real(rp)                     :: glres              ! Problem residuals
     character(6)                 :: namod              ! Module name
     character(3)                 :: exmod              ! Module extension
     integer(ip)                  :: lun_pdata          ! File units
     integer(ip)                  :: lun_outpu          ! ...
     integer(ip)                  :: lun_conve
     integer(ip)                  :: lun_rstpo
     integer(ip)                  :: lun_rstar
     integer(ip)                  :: lun_timin
     character(150)               :: fil_pdata          ! File names
     character(150)               :: fil_outpu          ! ...
     character(150)               :: fil_conve
     character(150)               :: fil_rstar
     character(150)               :: fil_timin
  end type tymod

end module def_kintyp_module
!> @}
