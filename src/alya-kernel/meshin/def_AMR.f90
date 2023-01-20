!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    mod_AMR.f90
!> @author  houzeaux
!> @date    2020-03-07
!> @brief   Adaptive mesh refinement
!> @details All tools for adaptive mesh refinement
!-----------------------------------------------------------------------

module def_AMR

  use def_kintyp_basic,         only : ip,rp,lg
  use def_coupli,               only : typ_color_coupling
  use def_search_method,        only : search_method
  use def_interpolation_method, only : interpolation
  use def_search_strategy,      only : search
  implicit none

  private

  !----------------------------------------------------------------------
  !
  ! Interpolation strategies
  !
  !----------------------------------------------------------------------

  integer(ip), parameter   :: AMR_ON_NODE_ELEMENT_INTERPOLAITON = 0
  integer(ip), parameter   :: AMR_ON_NODE_NEAREST_NODE          = 1
  integer(ip), parameter   :: AMR_ON_ELEMENT                    = 2
  integer(ip), parameter   :: AMR_BOUNDARY                      = 3
  
  !----------------------------------------------------------------------
  !
  ! Remeshing strategies
  !
  !----------------------------------------------------------------------

  integer(ip), parameter   :: AMR_COPY_MESH                     = 0
  integer(ip), parameter   :: AMR_GMSH                          = 1
  integer(ip), parameter   :: AMR_ALYA                          = 2
  
  !----------------------------------------------------------------------
  !
  ! Error estimate
  !
  !----------------------------------------------------------------------

  integer(ip), parameter   :: AMR_SHARP_LAPLACIAN               = 0
  integer(ip), parameter   :: AMR_LAPLACIAN                     = 1
  integer(ip), parameter   :: AMR_DISTANCE_POINT                = 2
  integer(ip), parameter   :: AMR_HESSIAN                       = 3
  
  !----------------------------------------------------------------------
  !
  ! Other (to classify in some category)
  !
  !----------------------------------------------------------------------
  ! 
  integer(ip), parameter   :: AMR_GIVEN_VARIABLE               = -1
  integer(ip), parameter   :: AMR_QUALITY_VARIABLE             = -2
  real(rp),    parameter   :: AMR_factorMaxSize                = 3.0!5.0     ! factor to determine maximum allowed mesh size 

  
  !----------------------------------------------------------------------
  !
  ! Input variables
  !
  !----------------------------------------------------------------------

  integer(ip)              :: kfl_amr                 ! AMR
  integer(ip)              :: kfl_amr_post            ! Postprocess technique for AMR
  integer(ip)              :: kfl_amr_freq            ! Frequency
  integer(ip)              :: nelem_amr               ! Target number of elements
  integer(ip)              :: npoin_amr               ! Target number of nodes
  integer(ip)              :: maxit_amr               ! Maximum number of iterations
  integer(ip)              :: kfl_amr_varia           ! Variable
  integer(ip)              :: kfl_amr_remesh          ! Remeshing method
  integer(ip)              :: kfl_amr_repart          ! Repartitioining method
  integer(ip)              :: kfl_size_amr            ! Size strategy
  integer(ip)              :: kfl_mesh_size_amr       ! Mesh size interpolaiton strategy  
  integer(ip)              :: kfl_error_amr           ! Error estimate
  real(rp)                 :: min_size_amr            ! Min size fot eh mesh
  real(rp)                 :: max_size_amr            ! Max size for the mesh
  real(rp)                 :: amrp0(3)                ! Center of the amr refinement
  real(rp)                 :: mesh_size_param_amr(5)  ! Parameters for mesh size interpolation
  
  real(rp)                 :: AMR_qualityTarget       ! Target quality to attain (sets the region to perform adaptivity)
  real(rp)                 :: AMR_qualityActivation   ! Quality to activate adaptation (below this, it is not activated)
  logical(lg)              :: AMR_bouSmoothSize       ! Do smooth target mesh size (from sol) to match boundary mesh size (fixed)
  logical(lg)              :: AMR_anisotropic         ! Do smooth target mesh size (from sol) to match boundary mesh size (fixed)
  integer(ip)              :: AMR_LPnorm              ! Lp norm to set hessian metric
  integer(ip)              :: AMR_iniTimeStep         ! Time step (ittim) to start AMR
  integer(ip)              :: AMR_numVariables        ! Number of variables to adapt to
  integer(ip)              :: kfl_amr_varia_multi(10) ! Multiple variables (maximum 10... can be increased if necessary)
  

  !----------------------------------------------------------------------
  !
  ! Local variables
  !
  !----------------------------------------------------------------------

  class(interpolation),      pointer :: interp_AMR_npoin
  class(interpolation),      pointer :: interp_AMR_nelem
  class(interpolation),      pointer :: interp_AMR_nboun
  type(search),              target  :: search_vol_seq
  type(search),              target  :: search_vol_par
  type(search),              target  :: search_bou_seq
  type(search),              target  :: search_bou_par

  class(search_method),    pointer :: search_amr                 ! Search method for mesh size

  integer(ip)              :: num_amr                    ! Number of mesh refinement steps
  integer(ip)              :: npoin_new
  integer(ip)              :: nelem_new
  integer(ip)              :: nboun_new

  logical(lg)              :: is_alloca_hmesh                  = .false.
  real(rp), pointer        :: hmesh(:)

  public :: AMR_SHARP_LAPLACIAN
  public :: AMR_LAPLACIAN
  public :: AMR_DISTANCE_POINT
  public :: AMR_HESSIAN
  
  public :: AMR_COPY_MESH
  public :: AMR_GMSH
  public :: AMR_ALYA
  
  public :: AMR_GIVEN_VARIABLE
  public :: AMR_QUALITY_VARIABLE
  public :: AMR_factorMaxSize
  
  public :: AMR_qualityTarget
  public :: AMR_qualityActivation
  public :: AMR_bouSmoothSize
  public :: AMR_anisotropic
  public :: AMR_LPnorm
  public :: is_alloca_hmesh
  public :: hmesh
  public :: AMR_iniTimeStep
  public :: AMR_numVariables
  public :: kfl_amr_varia_multi

  public :: num_amr
  public :: kfl_amr
  public :: kfl_amr_post
  public :: kfl_amr_freq
  public :: nelem_amr
  public :: npoin_amr
  public :: maxit_amr
  public :: kfl_amr_varia
  public :: kfl_amr_remesh
  public :: kfl_amr_repart
  public :: kfl_size_amr 
  public :: kfl_error_amr 
  public :: min_size_amr 
  public :: max_size_amr 
  public :: kfl_mesh_size_amr 
  public :: mesh_size_param_amr 

  public :: interp_AMR_npoin
  public :: interp_AMR_nelem
  public :: interp_AMR_nboun
  public :: search_amr
  public :: search_vol_seq
  public :: search_vol_par
  public :: search_bou_seq
  public :: search_bou_par
  
  public :: npoin_new
  public :: nelem_new
  public :: nboun_new
  public :: amrp0
  
end module def_AMR
!> @}

