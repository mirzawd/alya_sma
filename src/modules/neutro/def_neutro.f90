!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module def_neutro

  !------------------------------------------------------------------------
  !    
  ! Heading for the NEUTRO routines
  !
  !------------------------------------------------------------------------

  use def_kintyp, only : ip,rp,i1p,i2p,r3p
  use def_kintyp, only : bc_nodes
  use def_kintyp, only : bc_bound
  use mod_ADR,    only : ADR_typ

  !------------------------------------------------------------------------
  !
  ! Parameters
  !
  !------------------------------------------------------------------------

  character(150)           :: &
       fil_rstar_neu,         &
       fil_outCUR_neu='Current.csv',         &
       fil_outFLX_neu='flux.csv',         &
       fil_outRAD_neu='rad_neu.csv'
  character(150),allocatable           :: &
       fil_totalXS_neu(:),         &
       fil_scattXS_neu(:),         &
       fil_fisionXS_neu(:),         &
       fil_source_neu(:),          &
       fil_kerma_neu(:)
       
  integer :: unit_curr_neu=3112  
  integer :: unit_flux_neu=3113
  integer :: unit_radd_neu=3114
  integer :: unit_totXs_neu=3115  
  integer :: unit_SctXs_neu=3116
  integer :: unit_fisXs_neu=3117
  integer :: unit_source_neu=3118
  integer :: unit_kerma_neu=3119

  real(rp),      parameter :: &
       zeneu = epsilon(1.0_rp)

  !--BEGIN REA GROUP
  !------------------------------------------------------------------------
  !
  ! Physical problem: read in neu_reaphy
  !
  !------------------------------------------------------------------------
 
  integer(ip)               :: &
       num_energies_neu,       &       ! Number of energy groups
       num_directions_neu,     &       ! Number of directions
       num_sources_neu,       &        ! Number of boundary neutron sources
       num_legendre_neu,     &         ! Number of order of legendre
       num_legendre_lee,     &         ! Number of order of legendre que va a leer
       num_materials_neu,   &          ! number of materials
       num_isotopes_max_neu=0_ip,  &      ! maximum number of isotopes for effective material
       max_val_acum_neu, &              ! auxiliar para matriz de grupos sparse
       kfl_icosa_neu,          &       !
       kfl_snord_neu,           &                !
       kfl_stead_neu                    ! Steady-state has been reached 

  integer(ip) :: output_level_neu = 0_ip ! Verbosity level for neutro

  real(rp)                  :: & 
        albedo_neu                     ! albedo coef for reflective boundaries

 integer(ip) :: units_factor_neu=2_ip !0        ! 0 es nada, 1 es cm, 2 es metros
 real(rp) :: uma_a_kg = 1.6605389e-27_rp ! kg
 real(rp) :: barns_metros = 1.e-28_rp ! m^2
 real(rp) :: uma_a_g = 1.6605389e-24_rp ! g
 real(rp) :: barns_cm = 1.e-24_rp ! cm^2
 real(rp) :: ev_a_joule = 1.602176634e-19_rp ! eV to Joule

   integer(ip),  pointer    ::     & !,allocatable
        efectivos_neu_in(:),    & ! minimo grupo no nulo
        efectivos_neu_out(:),   & ! maximo grupo no nulo
        num_isotopes_neu(:), &  ! cantidad de isotopos por material
        efectivos_neu_sparse(:,:), &   ! grupos no nulos por grupo de energia por material (simil spasrse)
        n_efectivos_neu_sparse(:,:), &     ! cantidad de grupos no nulos por grupo de energia por material
        n_efectivos_neu_sparse_acum(:,:)

   real(rp),pointer ::     &
       source_bound_neu(:,:),    &       ! sources in boundaries por grupo
       absor_neu(:,:),   &       ! absorbcion coefficient by material and by energy group
       scatt_neu(:,:,:,:), &     ! scattering matrix coefficient by material and by energy group and legendre order
       fiss_neu(:,:,:), &        ! fission matrix coefficient by material and by energy group
       funsource(:),  &          ! fuente by material
       At_weight(:),  &         ! peso atomico por material  Uma
       Densidad_(:),  &         ! densidad por material  ! Kg/m3
       Isotope(:), &              ! isotopo por material
       aniso_neu(:),&              ! Linear anisotropic coefficient
       absor_neu_cte(:),   & ! absorcion cte sino se lee
       scatt_neu_cte(:),  &     ! scattering cte sino se lee
       grupo_energias(:,:), &    ! intervalos de enrgias utilizados
       phi_neu(:),      &   ! phi de cada direccion
       tita_neu(:),     &    ! tita de cada direccion 
       At_weight_isotope(:,:),  &         ! peso atomico por isotopo de material efectivo  Uma
       mass_percentage_isotope_neu(:,:), &     ! porcentaje en masa de cada isotopo de material efectivo [0-1]
       atom_percentage_isotope_neu(:,:), &     ! porcentaje atomico de cada isotopo de material efectivo [0-1]
       kerma_neu(:,:),   &         ! kerma factors by material and energy group
       kerma_poin_neu(:,:)          ! kerma factor for each node and energy group


  !> Variables only used in neu_readXS, but defined here to avoid "may be used uninitialized" warning
  real(rp),allocatable :: grupo_energias_isotope(:,:,:), absor_neu_isotope(:,:,:), &
                          scatt_neu_isotope(:,:,:,:,:), kerma_isotope_neu(:,:,:)

  !> Variables only used in neu_reaphy, but defined here to avoid "may be used uninitialized" warning
  real(rp),allocatable :: At_weight_isotope_aux(:,:), mass_percentage_isotope_neu_aux(:,:), &
                          atom_percentage_isotope_neu_aux(:,:) 
  !------------------------------------------------------------------------
  !
  ! Numerical problem: read in neu_reanut
  !
  !------------------------------------------------------------------------

  integer(ip)               :: &
       miinn_neu,              &       ! Maximum inner iterations
       kfl_smobo_neu                   ! B.c. smoothing
  real(rp)                  :: & 
       cotol_neu,              &       ! Tolerance inner iterations       
       relax_neu,              &       ! Relaxation
       nitsche_neu                     ! Nitsche coefficient

  !------------------------------------------------------------------------
  !
  ! Boundary conditions: read in neu_reabcs
  !
  !------------------------------------------------------------------------

  type(bc_nodes), pointer  :: &     
       tncod_neu(:)                  ! Node code type
  type(bc_bound), pointer  :: &     
       tbcod_neu(:)                  ! Boundary code type

  !------------------------------------------------------------------------
  !
  ! Output and Postprocess: read in neu_reaous
  !
  !------------------------------------------------------------------------

  !--END REA GROUP
  !------------------------------------------------------------------------
  !
  ! Others
  !
  !------------------------------------------------------------------------
  !
  ! Dimensions, etc.
  !
  integer(ip)              :: &
       nunkn_neu,             &      ! Number of unknowns
       ncomp_neu,             &      ! Number of components
       nprev_neu,             &      ! Last time step or global iteration
       current_energy_neu,    &      ! Current energy being solved
       current_direction_neu, &      ! Current direction being solved
       kfl_goite_neu                 ! Continue inner iterations
  real(rp)                 :: &
       resid_neu,             &      ! Residual of outer iterations
       residual_neu

  real(rp), pointer :: resid_energy_group_neu(:)

  !
  ! Directions
  !
  real(rp),  pointer       :: &
       direc_neu(:,:),        &      ! Directions
       weigd_neu(:),          &      ! Weights of directions
       ener_weigd_neu(:),     &      ! Weights of energy
       scattering_neu(:,:)          ! Scattering coefficient


  !
  ! Boundary conditions
  !
  type(i2p),   pointer     :: &
       kfl_fixno_neu(:,:)            ! Nodal fixity   CODES= 1 ::> Vacuum ; 2 ::> Albedo ; 3 ::> Reflex
  type(i1p),   pointer     :: &
       kfl_fixbo_neu(:,:)            ! Element boundary fixity
  type(i1p),   pointer     :: &
       kfl_funbo_neu(:,:), &            ! Element boundary function number (for NEUTRO, which source when multiple)
       kfl_funtb_neu(:,:)               ! function type (not used, but needs to be allocated if using funbo)
  type(r3p),   pointer     :: &
       bvess_neu(:,:),        &      ! Essential velocity bc values
       bvnat_neu(:,:)                ! Natural bc values
  !
  ! ADR type
  !
  type(ADR_typ)            :: &
       ADR_NEU                       ! ADR type


  interface
     subroutine runend(message) 
       implicit none
       character(*)          :: message
     end subroutine runend
  end interface

end module def_neutro
