!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kinds_and_types
!> @{
!> @file    def_kintyp_postprocess.g90
!> @author  houzeaux
!> @date    2020-04-04
!> @brief   Functions
!> @details Communications
!-----------------------------------------------------------------------

module def_kintyp_postprocess

  use def_kintyp_basic, only : ip,rp,lg
  
  !----------------------------------------------------------------------
  !
  ! Postprocess
  !
  !----------------------------------------------------------------------

  integer(ip), parameter ::     nvars     = 200                 ! # set variables
  integer(ip), parameter ::     nvart     = 10                  ! # times for postprocess
#ifdef ALYA_NVARW
  nvarw = ALYA_NVARW
#else
  integer(ip), parameter ::     nvarw     = 60                  ! # witness points
#endif  
#ifdef ALYA_NVARG
  nvarg = ALYA_NVARG
#else  
  integer(ip), parameter ::     nvarg     = 60                  ! # witness geometries
#endif  
  integer(ip), parameter ::     nvarp     = 300                 ! # postprocess variables
  integer(ip), parameter ::     nvati     = 10                  ! # max number of time steps
  integer(ip), parameter ::     nvarm     = nvarg               ! # max number of meshes
  integer(ip), parameter ::     ncomp_max = 20                  ! # max number of components to postprocess

  type typos_mesh
     ! Applies to all meshes
     real(rp)               :: pos_tinit                        ! Postprocess initial time
     integer(ip)            :: npp_inits                        ! Postprocess initial step
     integer(ip)            :: npp_iniso                        ! Postprocess initial condition
     logical(lg)            :: rst_time(nvati,nvarp)            ! Time steps to postprocess
     integer(ip)            :: kfl_oonce(nvarp)                 ! Postprocess only once
     ! Mesh dependent variables
     integer(ip)            :: npp_stepi(nvarp,0:nvarm)         ! Postprocess step interval for u,p, etc.
     integer(ip)            :: vox_stepi(nvarp,0:nvarm)         ! Postprocess step interval for u,p, etc.
     integer(ip)            :: pos_alrea(nvarp,0:nvarm)         ! Already postprocessed
     integer(ip)            :: vox_alrea(nvarp,0:nvarm)         ! Already postprocessed
     integer(ip)            :: lcomp(ncomp_max,nvarp,0:nvarm)   ! List of components to postprocess
     real(rp)               :: pos_times(nvart,nvarp,0:nvarm)   ! Postprocess times for u,p, etc.
     real(rp)               :: pos_perio(nvarp,0:nvarm)         ! Postprocess time period for u,p, etc.
     real(rp)               :: vox_times(nvart,nvarp,0:nvarm)   ! Postprocess times for u,p, etc.
     real(rp)               :: vox_perio(nvarp,0:nvarm)         ! Postprocess time period for u,p, etc.
  end type typos_mesh
  
  type, extends(typos_mesh) :: typos
     integer(ip)            :: npp_stepelset            ! Postprocess element set interval
     integer(ip)            :: npp_stepnoset            ! Postprocess node set interval
     integer(ip)            :: npp_stepboset            ! Postprocess boundary set interval
     integer(ip)            :: npp_stepw                ! Postprocess witness point interval
     integer(ip)            :: npp_stepg                ! Postprocess witness geometry interval
     integer(ip)            :: npp_witne(nvarw)         ! Postprocess witness points
     integer(ip)            :: npp_setse(nvars)         ! Postprocess element sets calculation
     integer(ip)            :: npp_setsb(nvars)         ! Postprocess boundary sets calculation
     integer(ip)            :: npp_setsn(nvars)         ! Postprocess node sets calculation
     integer(ip)            :: per_setse(nvars)         ! Postprocess element sets calculation
     integer(ip)            :: per_setsb(nvars)         ! Postprocess boundary sets calculation
     integer(ip)            :: per_setsn(nvars)         ! Postprocess node sets calculation
     integer(ip)            :: npp_witng(nvarg)         ! Postprocess witness geometries
     integer(ip)            :: nvaes                    ! Element set variables
     integer(ip)            :: nvabs                    ! Boundary set variables
     integer(ip)            :: nvans                    ! Node set variables
     integer(ip)            :: per_nvaes                ! Element set variables to exchange
     integer(ip)            :: per_nvabs                ! Boundary set variables to exchange
     integer(ip)            :: per_nvans                ! Node set variables to exchange
     integer(ip)            :: nvawi                    ! Witness point variable number
     integer(ip)            :: nvawg                    ! Witness geometry variable number
     integer(ip)            :: ipass                    ! Set memory allocated and header
     integer(ip)            :: lun_setse                ! Element set unit imodu*10+6
     integer(ip)            :: lun_setsb                ! Boundary set unit imodu*10+7
     integer(ip)            :: lun_setsn                ! Node set unit imodu*10+8
     integer(ip)            :: lun_setsi                ! Immersed bounday set
     integer(ip)            :: lun_witne                ! Witness points
     integer(ip)            :: lun_witng                ! Witness geometries
     character(150)         :: fil_setse
     character(150)         :: fil_setsb
     character(150)         :: fil_setsn
     character(150)         :: fil_setsi
     character(150)         :: fil_witne
     character(150)         :: fil_witng
     real(rp)               :: paese(5,nvars)           ! Element set parameters
     real(rp)               :: pabse(5,nvars)           ! Boundary set parameters
     real(rp)               :: panse(5,nvars)           ! Node set parameters
     character(5)           :: woese(nvars)             ! Name of the element set variables
     character(5)           :: wobse(nvars)             ! Name of the boundary set variables
     character(5)           :: wonse(nvars)             ! Name of the node set variables
     character(5)           :: wowit(nvarw)             ! Name of the witness points
     character(5)           :: wowig(nvarw)             ! Name of the witness geometries
     character(5)           :: wopos(5,nvarp)           ! Name and character of postprocess variable
     integer(ip)            :: enti_posit(nvarp)        ! Position of the entity in the array
     integer(ip)            :: comp_posit(nvarp)        ! Position of the component in the array
     integer(ip)            :: time_posit(nvarp)        ! Position of time index in the array
     integer(ip)            :: time_num(nvarp)          ! Number of time components
     integer(ip)            :: comp_num(nvarp)          ! Number of components
     integer(ip)            :: dime_num(nvarp)          ! Number of dimension if this is a vector
     integer(ip)            :: array_allocated(nvarp)   ! If variable has been allocated
     integer(ip)            :: array_used(nvarp)        ! If variable is used
     integer(ip)            :: array_registered(nvarp)  ! If variable has been registered
     real(rp),      pointer :: veset(:,:)               ! Set element values
     real(rp),      pointer :: vbset(:,:)               ! Set boundary values
     real(rp),      pointer :: vnset(:,:)               ! Set node values
     real(rp),      pointer :: viset(:,:)               ! Set IB values
     real(rp),      pointer :: witne(:,:)               ! Witness point values
     real(rp),      pointer :: witng(:,:)               ! Witness geometry values
     integer(ip)            :: witne_kfl_time(nvarg)    ! Witness points:   Time strategy
     real(rp)               :: witne_period(nvarg)      ! Witness points: Averaging period
     real(rp)               :: witne_average(nvarg)     ! Witness points: Averaging time
     real(rp)               :: witne_dt(nvarg)          ! Witness points: Averaging time step
     integer(ip)            :: witng_kfl_time(nvarg)    ! Witness geometry: Time strategy
     real(rp)               :: witng_perio_max(nvarg)   ! Witness geometry: Averaging: max of accummulated period at wich we restart
     real(rp),      pointer :: witng_perio_count(:,:)   ! Witness geometry: Averaging: accummulated period
     real(rp),      pointer :: witng_denom(:,:)         ! Witness geometry: Averaging: accummulated denominator
     real(rp),      pointer :: witng_deldenom(:,:)      ! Witness geometry: Averaging: change in denominator
     real(rp),      pointer :: witng_delperio(:,:)      ! Witness geometry: Averaging: change in period
  end type typos
  
  type witness_geo
     integer(ip)            :: kfl_geometry                 ! Geometry type 
     real(rp)               :: param(9)                     ! Geometry parameters
   contains
     procedure,  pass       :: init   => init_witness_geo   ! Initialization
     procedure,  pass       :: deallo => deallo_witness_geo ! Initialization
  end type witness_geo
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-27
  !> @brief   Initialization
  !> @details Initialization of class WITNESS_GEO
  !> 
  !-----------------------------------------------------------------------

  subroutine init_witness_geo(self)
    class(witness_geo), intent(inout) :: self
    
    self % kfl_geometry = 0                      ! No type 
    self % param        = 0.0_rp                 ! Geometry parameters

  end subroutine init_witness_geo
    
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-27
  !> @brief   Deallocate
  !> @details Deallocate class WITNESS_GEO
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo_witness_geo(self)
    class(witness_geo), intent(inout) :: self
    
  end subroutine deallo_witness_geo
    
end module def_kintyp_postprocess
!> @}
