!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    def_gusano.f90
!> @author  houzeaux
!> @date    2020-04-22
!> @brief   Gusano definitions
!> @details All definitions of Gusano
!-----------------------------------------------------------------------
module def_gusano

  use def_kintyp_basic,               only : ip,rp,i1p
  use def_kintyp_boundary_conditions, only : bc_nodes
  use def_kintyp_boundary_conditions, only : bc_bound
  implicit none
  !
  ! Parameters
  !
  integer(ip), parameter :: GUS_ASGS                = 1
  integer(ip), parameter :: GUS_OSS                 = 2
  integer(ip), parameter :: GUS_SPLIT_OSS           = 3
  
  integer(ip), parameter :: GUS_MONOLITHIC          = 0
  integer(ip), parameter :: GUS_SCHUR_COMPLEMENT    = 1
  
  integer(ip), parameter :: GUS_LAMINAR_PIPE        = 0
  integer(ip), parameter :: GUS_TURBULENT_PIPE      = 1
  integer(ip), parameter :: GUS_BEND_OFF            = 0
  integer(ip), parameter :: GUS_BEND_LAMINAR_PIPE   = 1
  integer(ip), parameter :: GUS_BEND_TURBULENT_PIPE = 2
  !
  ! Input
  !
  type(bc_nodes), pointer  :: &     
       tncod_gus(:)                  ! Node code type
  type(bc_bound), pointer  :: &     
       tbcod_gus(:)                  ! Boundary code type
  !
  ! Physical poblem
  !
  integer(ip)              :: &
       kfl_timei_gus,         &      ! Time integration off
       kfl_conve_gus,         &      ! Convective term
       kfl_regim_gus,         &      ! Regime to determine velocity profile and friction
       kfl_bendi_gus                 ! Bending model
  !
  ! Numerical treatment
  !
  integer(ip)              :: &
       kfl_stabi_gus,         &      ! Stabilization strategy
       kfl_algor_gus                 ! Solution algorithm
  
  real(rp)                 :: &
       cotol_gus,             &      ! Non-linearity tolerance
       sstol_gus,             &      ! Steady state tolerance
       safet_gus                     ! Safet factor for time step
  !
  !
  ! Boundary conditions
  !
  integer(ip), pointer     :: &
       kfl_fixno_gus(:,:),    &      ! Nodal fixity 
       kfl_fixbo_gus(:),      &      ! Function
       kfl_funno_gus(:),      &      ! Function
       kfl_funbo_gus(:),      &      ! Function
       kfl_funtn_gus(:),      &      ! Function
       kfl_funtb_gus(:)              ! Function
  real(rp),    pointer     :: &
       bvess_gus(:,:,:),      &      ! Essential velocity bc values
       bvnat_gus(:,:,:)
  !
  ! Time
  !
  integer(ip)              :: &       
       kfl_stead_gus,         &
       kfl_goite_gus
  real(rp)                 :: &
       dtcri_gus,             &
       dtinv_gus,             &
       dtmax_gus,             &
       saflo_gus
  !
  ! Stabilization
  !
  real(rp)                  :: &
       xoss_gus,               &
       xsoss_gus
  !
  ! Arrays
  !
  integer(ip), pointer     :: &
       neuman_gus(:)                 ! Neumann condition
  type(i1p),  pointer      :: &
       dirich_gus(:)                 ! Dirichlet condition
  real(rp),   pointer      :: &
       projm_gus(:,:),        &      ! Momentum residual projection (pressure for split oss)
       projc_gus(:,:),        &      ! Continuity residual projection
       exn1d_gus(:),          &      ! 1D normal
       angle_gus(:),          &      ! Angles for bending model
       bendi_gus(:,:),        &      ! Bending array (radius of curvature)
       densi_gus(:),          &      ! Density
       visco_gus(:),          &      ! Viscosity
       vel3d_gus(:,:),        &      ! 3D velocity
       pre3d_gus(:)                  ! 3D pressure
  !
  ! Schur complement pointers for matrices
  !
  integer(ip), pointer     :: &
       kfl_fixsc_gus(:,:)            ! Pressure schur complement prescription
  real(rp),    contiguous,    &
       &       pointer     :: &
       schur_gus(:)                  ! Pressure Schur complement matrix
  
end module def_gusano
!> @}
