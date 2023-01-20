!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    def_nastin_parameters.f90
!> @author  houzeaux
!> @date    2020-04-22
!> @brief   Parameters
!> @details parameters used in Nastin
!-----------------------------------------------------------------------

module def_nastin_parameters

  use def_kintyp_basic
  !
  ! Global parameters
  !
  integer(ip),   parameter :: &
       NSI_INCREMENTAL_PROJECTION          =  2, &
       NSI_PREDICTOR_CORRECTOR             =  3, &
       NSI_BLOCK_GAUSS_SEIDEL              =  4, &
       NSI_MOMENTUM                        =  2, &
       NSI_CONTINUITY                      =  3, &
       NSI_MOMENTUM_AND_CONTINUITY         =  1, &
       NSI_INCOMPRESSIBLE                  =  0, &
       NSI_COMPRESSIBLE                    =  1, &
       NSI_LOW_MACH                        =  3, &
       NSI_ANALYTICAL_HYDROSTATIC_PRESSURE =  1, &
       NSI_PDE_HYDROSTATIC_PRESSURE        =  2, &
       NSI_GALERKIN                        = -1, &
       NSI_ASGS                            =  0, &
       NSI_OSS                             =  1, &
       NSI_SPLIT_OSS                       =  2, &
       NSI_ALGEBRAIC_SPLIT_OSS             =  3, &
       NSI_SUPG                            =  4, &
       NSI_CONVECTION_NON_CONSERVATIVE     =  0, &
       NSI_CONVECTION_CONSERVATIVE         =  1, &
       NSI_CONVECTION_SKEW                 =  2, &
       NSI_CONVECTION_EMAC                 =  3, &
       NSI_CONVECTION_EMAC2                =  4, &
       NSI_DIRICHLET_ELEMENT               =  0, &
       NSI_DIRICHLET_MATRIX                =  1, &
       NSI_DIRICHLET_ALGORITHM             =  2, &
       NSI_GLOBAL_TO_LOCAL                 =  1, &
       NSI_LOCAL_TO_GLOBAL                 =  2, &
       NSI_LUMPED_MASS                     =  0, &
       NSI_CONSISTENT_MASS                 =  1
!
! density and vsicosity to be use in mod_nsi_element_operations_hh*  &  mod_nsi_multi_step_fs.f90 for kfl_assem_nsi >39_ip
!
  real(rp),parameter                   :: densi_aux = 1.229000D+00
  real(rp),parameter                   :: visco_aux = 0.1730000D-05
  !
  ! Solvers
  !
  integer(ip),   parameter :: &
       NSI_SOLVER_NAVER_STOKES             = 1, & ! For monolithic
       NSI_SOLVER_MOMENTUM                 = 1, & ! For non-monolithic
       NSI_SOLVER_CONTINUITY               = 2, & ! For non-monolithic
       NSI_SOLVER_VISCOUS_TERM             = 9, & ! For semi-implicit
       NSI_SOLVER_BOUNDARY_CONDITIONS      = 3, & ! 
       NSI_SOLVER_MASS_CORRECTION          = 4, & ! 
       NSI_SOLVER_HYDROSTATIC_PRESSURE     = 5, & ! 
       NSI_SOLVER_ZERO_DIVERGENCE          = 6, & ! 
       NSI_SOLVER_NORMAL_EXTENSIONS        = 7, & ! 
       NSI_SOLVER_CONSISTENT_MASS          = 8    ! 
  !
  ! File units
  !
  integer(ip),   parameter :: &
       lun_bound_nsi = 110,   &
       lun_stasg_nsi = 112,   &
       lun_cvgsg_nsi = 113,   &
       lun_psmat_nsi = 114,   &
       lun_refer_nsi = 117,   &
       lun_recvg_nsi = 118,   &
       lun_lmach_nsi = 120,   &
       lun_dynin_nsi = 122,   &
       lun_dynou_nsi = 123,   &
       lun_dynlo_nsi = 124,   &
       lun_dynre_nsi = 125,   &
       lun_windt_nsi = 1100
  !
  ! Variable dimensions
  !
  integer(ip),   parameter :: &
       nvars_nsi=25,          &      ! # set variables
       nvart_nsi=10,          &      ! # times for postprocess
       nvarw_nsi=10,          &      ! # witness point
       nvarp_nsi=40,          &      ! # postprocess variables
       ncoef_nsi=10,          &      ! # coefficient for properties
       mtabl_nsi=500,         &      ! maximum number of tabulated ct cp coefs for disk properties
       mforc_material_nsi=21, &      ! Maximum number of material force parameters 
       ivari_nsi_mom=1,       &      ! Momentum 
       ivari_nsi_cont=2,      &      ! Continuity 
       ivari_nsi_corre=4,     &      ! end of step correction 
       ivari_nsi_hydro=5,     &      ! hydrostatic state
       ivari_nsi_divfree=6,   &      ! Divergence Free
       mmsgs_nsi=50                  ! Max. number of iterations for SGS convergence statistics

  real(rp),      parameter :: &
       zensi = epsilon(1.0_rp)

end module def_nastin_parameters
