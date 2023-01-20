!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Moduls
!> @{
!> @file    mod_moduls_conf.f90
!> @author  houzeaux
!> @date    2019-06-19
!> @brief   Calling subroutines to the modules of Alya
!> @details Warning: IBLOK and MODUL are global variables
!-----------------------------------------------------------------------

module mod_moduls_conf

  use def_kintyp_basic,       only : ip,rp
  use def_master
  use def_kermod
  use def_domain
  use def_solver
  use def_inpout
  use mod_ker_proper
  use mod_ker_timeline
  use def_kintyp_performance, only : perf_alloca
  use mod_parall,             only : par_code_zone_subd_to_color
  use mod_parall,             only : PAR_COMM_COLOR_ARRAY
  use mod_parall,             only : PAR_COMM_WORLD,commd
  use mod_communications,     only : PAR_BARRIER
  use mod_communications,     only : PAR_MAX
  use mod_messages,           only : livinf
  use mod_timings,            only : timings_doiter
  use mod_alya2talp,          only : alya2talp_MonitoringRegionStart
  use mod_alya2talp,          only : alya2talp_MonitoringRegionStop
  implicit none
  private

  public :: moduls_set_current_module  ! Set current module
  public :: moduls_solver              ! Solver pointer
  public :: moduls_set_current_solver  ! Set current solver of a module
  public :: moduls_allocate_timing     ! Allocate timing structure
  public :: moduls_in_block            ! Module in a block

contains
  
  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2020-04-09
  !> @brief   Allocate timers
  !> @details Allocate timers
  !>
  !-----------------------------------------------------------------------

  subroutine moduls_allocate_timing(nn,CURRENT_MODULE)

    integer(ip),           intent(in) :: nn
    integer(ip), optional, intent(in) :: CURRENT_MODULE
    integer(ip)                       :: kmodu,ii

    kmodu = optional_argument(modul,CURRENT_MODULE)
    
    if( nn > 0 ) then
       call perf_alloca(momod(kmodu) % times,nn)
       times(1:) => momod(kmodu) % times(1:)      ! Lower bound (1:) is required by PGI
       do ii = 1,nn
          call momod(kmodu) % times(ii) % init()
       end do
    end if

  end subroutine moduls_allocate_timing
  
  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2020-04-09
  !> @brief   Update module
  !> @details Pointer to current module structures
  !>
  !-----------------------------------------------------------------------

  subroutine moduls_set_current_module(imodu)

    integer(ip), optional, intent(in) :: imodu
    integer(ip)                       :: kmodu

    kmodu = optional_argument(modul,imodu)
    if( present(imodu) ) modul = imodu

    postp     => momod(kmodu) % postp
    solve     => momod(kmodu) % solve
    eigeg     => momod(kmodu) % eigen
    solve_sol => momod(kmodu) % solve
    eigen_sol => momod(kmodu) % eigen
    times(1:) => momod(kmodu) % times(1:)  ! Lower bound (1:) is required by PGI
    tncod     => momod(kmodu) % tncod
    tgcod     => momod(kmodu) % tgcod
    tbcod     => momod(kmodu) % tbcod
    reset     => momod(kmodu) % reset
    if( associated(postp) ) then
       veset     => postp(1)     % veset
       vbset     => postp(1)     % vbset
       vnset     => postp(1)     % vnset
       witne     => postp(1)     % witne
       witng     => postp(1)     % witng
    end if
    
    if( kmodu > 0 ) then
       current_zone = lzone(kmodu)
    else
       current_zone = 0
    end if

  end subroutine moduls_set_current_module

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2020-04-21
  !> @brief   Point to a solver
  !> @details Point to a solver
  !>
  !-----------------------------------------------------------------------

  function moduls_solver(imodu,vanam) result(solverp)

    integer(ip),      optional, intent(in) :: imodu
    character(len=*), optional, intent(in) :: vanam
    integer(ip)                            :: ivari
    type(soltyp),     pointer              :: solverp(:)

    do ivari = 1,size(momod(imodu) % solve)
       if( trim(momod(imodu) % solve(ivari) % wprob) == vanam ) then
          solverp => momod(imodu) % solve(ivari:)
          return
       end if
    end do
    call runend('MOD_MODULS: SOLVER NOT FOUND')

  end function moduls_solver

  subroutine moduls_set_current_solver(imodu,vanam)

    integer(ip),      optional, intent(in) :: imodu
    character(len=*), optional, intent(in) :: vanam
    integer(ip)                            :: ivari

    do ivari = 1,size(momod(imodu) % solve)
       if( trim(momod(imodu) % solve(ivari) % wprob) == vanam ) then
          solve_sol => momod(imodu) % solve(ivari:)
          return
       end if
    end do
    call runend('MOD_MODULS: SOLVER NOT FOUND')

  end subroutine moduls_set_current_solver

  !-----------------------------------------------------------------------
  !> 
  !> @author  bsc21943
  !> @date    2022-01-13
  !> @brief   Tell if a module is in a block
  !> @details Tell if a module is in a block
  !> 
  !-----------------------------------------------------------------------

  pure function moduls_in_block(kblok,kmodu) result(in)

    integer(ip), intent(in) :: kblok
    integer(ip), intent(in) :: kmodu
    logical(lg)             :: in
    integer(ip)             :: ii,nn,kk

    nn = size(lmord,1)
    ii = 0
    loop_kk: do kk = 1,nn
       if( lmord(kk,kblok) == kmodu ) then
          ii = kk
          exit loop_kk
       end if
    end do loop_kk
    if( ii /= 0 ) then
       in = .true.
    else
       in = .false.
    end if
    
  end function moduls_in_block

end module mod_moduls_conf
!> @}
