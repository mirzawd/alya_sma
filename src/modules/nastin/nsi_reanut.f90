!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup NastinInput
!> @{
!> @file    nsi_reanut.f90
!> @author  Guillaume Houzeaux
!> @brief   Read numerical data
!> @details Read numerical data
!> @}
!-----------------------------------------------------------------------
subroutine nsi_reanut()
  !.md<module>nastin
  !.md<input>case.nsi.dat
  !.md<pos>1
  !.md<sec>
  use def_parame
  use def_inpout
  use def_master
  use def_nastin
  use def_domain
  use def_solver
  use mod_elmgeo, only : FLAT_BUBBLE
  use mod_elmgeo, only : QUADRATIC_BUBBLE
  use mod_elmgeo, only : FREE_SURFACE_BUBBLE
  use mod_ecoute, only : ecoute

  implicit none
  integer(ip) :: imome,icont
  integer(ip) :: kfl_nota2,kfl_ta2ma  ! using these two auxiliary flags I amke it independent of the order in which it is read

  if( INOTSLAVE ) then
     !
     !  Initializations (defaults)
     !
     staco_nsi(1)            = 1.0_rp                               ! Viscous part of tau1
     staco_nsi(2)            = 1.0_rp                               ! Convective part of tau1
     staco_nsi(3)            = 1.0_rp                               ! Rotation part of tau1
     staco_nsi(4)            = 1.0_rp                               ! Coefficient of tau2
     staco_corio_nsi         = 0.02_rp                              ! Coriolis stabilization paarmeter - recomended value
     tosgs_nsi               = 1.0e-2_rp                            ! Subgrid scale tolerance
     kfl_ellen_nsi           = 0                                    ! Minimum element length
     kfl_relax_nsi           = 0                                    ! Velocity: no relaxation
     kfl_relap_nsi           = 0                                    ! Pressure: no relaxation
     kfl_dttyp_nsi           = 1                                    ! Module time step strategy
     kfl_sgsco_nsi           = 0                                    ! Stabilization strategy
     kfl_sgsti_nsi           = 0                                    ! Stabilization strategy
     kfl_sgsac_nsi           = -1                                   ! SGS time accuracy
     kfl_sgsli_nsi           = 1                                    ! SGS convection linearization PICARD
     kfl_sgscp_nsi           = 0                                    ! Coupling of SGS with grid scale
     kfl_algor_nsi           = 5                                    ! Schur
     kfl_predi_nsi           = 3                                    ! Predictor-corrector like term
     kfl_taush_nsi           = 1                                    ! Schur complement. Tau strategy
     kfl_ellsh_nsi           = 0                                    ! Schur complement.           =0,1 for min/max element length
     kfl_updpr_nsi           = 0                                    ! Update pressure for Crank-Nicolson
     kfl_intpr_nsi           = 1                                    ! =1 pressure term integrated by parts
     kfl_assem_nsi           = 3                                    ! Assembly type (3,5)
     kfl_asbou_nsi           = 3                                    ! Assembly type (3,5)
     kfl_taust_nsi           = 1                                    ! Tau strategy
     kfl_stabi_nsi           = 0                                    ! Stabilization method
     kfl_limit_nsi           = 0                                    ! No limiter
     kfl_trres_nsi           = 1                                    ! Transient residual
     kfl_prtre_nsi           = 0                                    ! Pressure treatment (explicit/implicit)
     kfl_matdi_nsi           = 0                                    ! Impose bc at element level
     kfl_intfo_nsi           = 0                                    ! Internal force residual
     kfl_press_nsi           = 1                                    ! Integrate pressure term in momentum equations
     kfl_bubbl_nsi           = 0                                    ! Pressure enrichment
     kfl_fsgrb_nsi           = 0                                    ! Correction for Fractional Step with gravity or Boussinesq
     kfl_stab_corio_nsi      = 0                                    ! Coriolis stabilization
     fsrot_nsi               = 0.0_rp                               ! FSM in rotational form factor
     kfl_fscon_nsi           = 0                                    ! FSM consistent (solve Poisson equation each step)
     relax_nsi               = 1.0_rp                               ! Momentum relaxation
     relap_nsi               = 1.0_rp                               ! Pressure relaxation
     relsg_nsi               = 1.0_rp                               ! Subgrid scale relaxation
     kfl_shock_nsi           = 0                                    ! Shock capturing off
     shock_nsi               = 0.0_rp                               ! SC parameter
     kfl_tiacc_nsi           = 1                                    ! First order time integ.
     neule_nsi               = 1                                    ! Number of Euler time steps
     sstol_nsi               = 1.0e-5_rp                            ! Steady-satate tolerance
     safet_nsi               = 1.0e10_rp                            ! Safety factor
     bemol_nsi               = 0.0_rp                               ! Integration of convective term by parts
     kfl_normc_nsi           = 2                                    ! L2 norm for convergence
     kfl_refer_nsi           = 0                                    ! Error with reference solution
     momod(modul) % miinn    = 1                                    ! One internal iteration
     kfl_stain_nsi           = 0                                    ! Time step to start inner iterations
     kfl_immer_nsi           = 0                                    ! Immersed boundary method
     kfl_grad_div_nsi        = 0                                    ! Use grad and div matices for Aup and Apu
     misgs_nsi               = 1                                    ! Subgrid scale iterations
     cotol_nsi               = 0.1_rp                               ! Internal tolerance
     kfl_penal_nsi           = 0                                    ! No penalization
     kfl_prepe_nsi           = 0                                    ! No pressure penalization
     penal_nsi               = 0.0_rp                               ! Penalization factor
     prepe_nsi               = 0.0_rp                               ! Pressure penalization factor
     kfl_linea_nsi           = 1                                    ! Linearization (RHS   =0, Picard   =1, Newton   =2)
     npica_nsi               = 1                                    ! Number of Picard iteration (Newton's lin.)
     kfl_meshi_nsi           = 0                                    ! Mesh interpolator activation (OFF   =0,ON   =1)
     kfl_tisch_nsi           = 1                                    ! Time integration scheme
     kfl_savco_nsi           = 0                                    ! Save constant matrix
     kfl_corre_nsi           = 2                                    ! Do not correct momentum
     kfl_sosch_nsi           = 2                                    ! Momentum preserving Orthomin(1)
     kfl_modfi_nsi           = 0                                    ! Do not modify kfl_fixno_nsi
     kfl_expco_nsi           = 0                                    ! Treat the convective term explicitly, that is, assemble the matrix only in the first ortomin iteration
     kfl_addpr_nsi           = 0                                    ! Add contribution due to pressure in matrix side (do nothing BC') on wall law boundaries
     solve_sol               => solve(1:)                           ! Solver type
     xfree_nsi               = 10000000.0_rp                        ! X coordinate of the plane where to free
     kfl_grvir_nsi           = 0                                    ! Viscosity (laminar + turbulent) gradient in residual (default set to zero)
     kfl_hydro_nsi           = 0                                    ! How to compute hydristatic pressure
     kfl_update_hydro_nsi    = ITASK_INIUNK                         ! When do we update the hydrostatic pressure
     kfl_hydro_interface_nsi = 0                                    ! Constant interface height
     mitri_nsi               = 1                                    ! Number of Richardson iterations and FS + Tau solver
     kfl_adres_nsi           = 0                                    ! FS + Tau solver: adaptive residual
     adres_nsi               = 0.1_rp                               ! FS + Tau solver: adaptive residual factor
     toler_nsi               = 1.0e-6_rp                            ! FS + Tau solver: tolerance
     kfl_incre_nsi           = 1                                    ! Solve pressure in incremental form
     safex_nsi               = 1.0_rp                               ! Time function parameter for safety factor
     safma_nsi               = 1.0e9_rp                             ! Maximum safety factor
     safeo_nsi               = safet_nsi                            ! Initial safety factor
     saflo_nsi               = safet_nsi                            ! Global safety factor for local time step
     gamma_nsi               = 0.0_rp                               ! Pressure factor in momentum equations
     kfl_nota2               = 0                                    ! Internal use
     kfl_ta2ma               = 0                                    ! Internal use
     kfl_nota1_nsi           = 0                                    ! stabilize conv, react, corio in momentum
     kfl_ini_ts_guess_order_nsi = 1                                 ! First order extrapolation - u_guess = u_n
     kfl_vector_nsi          = 1                                    ! Vectorized version of assembly
     kfl_press_stab_nsi      = 0                                    ! Pressure stabilization for fractional step
     kfl_stop_by_wit_nsi     = 0                                    ! Stop by convergence of witness points - uses velocity
     kfl_massm_nsi           = NSI_LUMPED_MASS                      ! Mass matrix inf FS methods
     kfl_dampi_nsi           = 0                                    ! no damping
     top_r_damp_nsi          = 2000.0_rp                            ! top rayleigh damping
     top_v_damp_nsi          = 2000.0_rp                            ! top viscous  damping
     bot_r_damp_nsi          = 1600.0_rp                            ! bottom rayleigh damping
     bot_v_damp_nsi          = 1600.0_rp                            ! bottom viscous  damping
     val_r_damp_nsi          = 0.0001_rp                            ! value rayleigh damping
     mul_v_damp_nsi          = 1000.0_rp                            ! value to multiply laminar viscosity
     v_geo_damp_nsi(1)       = 10.0_rp                              ! x component of geostrophic velocity
     v_geo_damp_nsi(2)       = 0.0_rp                               ! y component of geostrophic velocity
     v_geo_damp_nsi(3)       = 0.0_rp                               ! z component of geostrophic velocity
     k_dim_damp_nsi          = 2_ip                                 ! vertical dimension 
     !
     ! Reach the section
     !
     call ecoute('nsi_reanut')
     do while(words(1)/='NUMER')
        call ecoute('nsi_reanut')
     end do
     !
     !.md# Numerical Treatment Definition
     !.md<code>
     !.md<0>NUMERICAL_TREATMENT
     !
     ivari_nsi = ivari_nsi_mom
     imome     = 0
     icont     = 0
     
     do while(words(1)/='ENDNU')
        call ecoute('nsi_reanut')

        if( words(1) == 'VECTO' ) then
           !
           ! Vectorized version of Alya
           !
           if( option('VECTO') ) then
              kfl_vector_nsi = 1
           else
              kfl_vector_nsi = 0
           end if

        else if( words(1) == 'TAUST' ) then
           !
           !.md<1>TAU_STRATEGY:               CODINA | SHAKIB | EXACT  | INCLUDE_TIME |TIME_STEP                                 $ Equation for stabilization parameter tau
           !.md<field>TAU_STRATEGY
           !.md<com> Equation for stabilization parameter tau1.
           !.md<com>
           !.md<com>    -  CODINA : tau1 = 1 / [ 4 mu/h^2 + 2 |u|/h + |w| ] 
           !.md<com>    -  SHAKIB : tau1 = 1 / [ 9*(4 mu/h^2)^2 + (2|u|/h)^2 + |w|^2 ]^1/2 
           !.md<com>    -  INCLUDE: tau1 = 1/ [ 1/dt +  4 mu/h^2 + 2 |u|/h + |w|]
           !.md<com>    -  TIMES     : tau1 = dt   
           !
           call reatau(kfl_taust_nsi)

        else if( words(1) == 'STABI' ) then
           !
           !.md<1>STABILIZATION:              ASGS | FULL_OSS | SPLIT_OSS[,FIRST_ORDER | SOTO | NO_LIMITER | GRADV] $ Stabilization strategy
           !.md<field>STABILIZATION
           !.md<com> Determines the stabilization technique of the Navier-Stokes equations.
           !.md<com>Both ASGS and FULL_OSS adds to the momentum and continuity equations a stabilization term.
           !.md<com>The difference if the equation of the subgrid scale.
           !.md<com>
           !.md<com>    -  ASGS:     Algebraic subgrid scale 
           !.md<com>    -  FULL_OSS: Orthogonal subgrid scale 
           !.md<com>    -  SPLIT:    Split orthogonal subgrid scale.
           !.md<com>The ASGS is generally more robust but less accurate. For the SPLIT_OSS, an additional option is required to set the
           !.md<com>limiter.
           !.md<com>    -  GRADV:       Force viscosity (laminar + turbulent) gradient in residual
           !.md<com>    -  NO_LIMITER:  No limiter is used 
           !.md<com>    -  FIRST_ORDER: SUPG method (only first order) 
           !.md<com>    -  SOTO:        Soto's limiter 
           !.md<com>
           !
           if( words(2) == 'OFF  ' .or. words(2) == 'GALER' ) then
              kfl_stabi_nsi = NSI_GALERKIN                               ! OFF

           else if( words(2) == 'ALGEB' ) then
              kfl_stabi_nsi = NSI_ALGEBRAIC_SPLIT_OSS                    ! Algebraic Split OSS (only pressure)

           else if( words(2) == 'ASGS ' ) then
              kfl_stabi_nsi = NSI_ASGS                                   ! ASGS
              if(exists('GRADV')) kfl_grvir_nsi = 1                      ! Force viscosity (laminar + turbulent) gradient in residual

           else if( words(2) == 'SUPG ' ) then
              kfl_stabi_nsi = NSI_SUPG                                   ! SUPG

           else if( words(2) == 'FULLO' ) then
              kfl_stabi_nsi = NSI_OSS                                    ! Full OSS
              if(exists('GRADV')) kfl_grvir_nsi = 1                      ! Force viscosity (laminar + turbulent) gradient in residual

           else if( words(2) == 'SPLIT' .or. words(2) == 'OSS  ' ) then
              kfl_stabi_nsi = NSI_SPLIT_OSS                              ! Split OSS

              if( exists('NOLIM') ) kfl_limit_nsi =  0                   ! No limiter
              if( exists('SOTO ') ) kfl_limit_nsi =  1                   ! Soto limiter
              if( exists('DIFFU') ) kfl_limit_nsi =  2                   ! Very diffusive limiter
              if( exists('FIRST') ) kfl_limit_nsi = -1                   ! First order

           else if (words(2) == 'CONST' ) then                           ! stability constants
              staco_nsi(1) = param(2)                                    ! c1/4 diff
              staco_nsi(2) = param(3)                                    ! c2/2 conv
              staco_nsi(3) = param(4)                                    ! c3   react
              staco_nsi(4) = param(5)                                    ! c4   pressu
              if(exists('TA2MA'))   kfl_ta2ma = 1                        ! applied at the end of the subroutine
              !              if(exists('NOTA2'))   kfl_nota2 = 1                        ! to avoid different behaviours depending on the order of the lines
              if ( nnpar /= 4 ) call runend ('NSI_REANUT : if stability constants are given all 4 must be given')

           end if
           if(exists('NOTA1'))   kfl_nota1_nsi = 1                           ! do not stabilize convection, reaction, coriolis in momentum eq.
           if(exists('NOTA2'))   kfl_nota2 = 1
        else if( words(1) == 'DAMPI' ) then
           kfl_dampi_nsi = 1
           top_r_damp_nsi          = param(1)                            ! top rayleigh damping
           top_v_damp_nsi          = param(2)                            ! top viscous  damping
           bot_r_damp_nsi          = param(3)                            ! bottom rayleigh damping
           bot_v_damp_nsi          = param(4)                            ! bottom viscous  damping
           val_r_damp_nsi          = param(5)                            ! value rayleigh damping
           mul_v_damp_nsi          = param(6)                            ! value to multiply laminar viscosity
           v_geo_damp_nsi(1)       = param(7)                            ! x component of geostrophic velocity
           v_geo_damp_nsi(2)       = param(8)                            ! y component of geostrophic velocity
           v_geo_damp_nsi(3)       = param(9)                            ! z component of geostrophic velocity
           k_dim_damp_nsi          = nint(param(10))                     ! vertical dimension
        else if( words(1) == 'TRACK' ) then
           !
           !.md<1>TRACKING_SUBGRID_SCALE:     TIME[, ORDER= 1 | 2 ], COUPLED | DECOUPLED \
           !.md<6>CONVECTION[, PICARD | NEWTON, ITERA= int1 , TOLER= real1 , RELAX= real2 ]  $ Subgrid scale tracking
           !.md<field>TRACKING_SUBGRID_SCALE
           !.md<com>Decide if the subgrid scale is tracked in convection and/or in time.
           !.md<com>CONVECTION for convection tracking. In the case, the subgrid scale is added to the convection velocity
           !.md<com>which appers both in the Galerkin term and in the equation for the subgrid scale. The scheme is therefore non-linear
           !.md<com>and requires a linearization procedure (PICARD/NEWTON) and convergence control parameters (ITERATIONS, TOLERANCE, RELAXATION)
           !.md<com>It subgrid scale is tracked in time (TIME option), then the order of the time discretization should be given: 1 or 2.
           !.md<com>COUPLED: SGS are computed at the same time as grid scale.
           !
           if(exists('CONVE') ) then
              kfl_sgsco_nsi = 1
              misgs_nsi     = getint('ITERA',1_ip,     '#Subgrid scale iterations')
              tosgs_nsi     = getrea('TOLER',1.0e-2_rp,'#Subgrid scale Tolerance')
              relsg_nsi     = getrea('RELAX',1.0_rp,   '#Subgrid scale Relaxation')
           end if
           if(exists('TIME ')) kfl_sgsti_nsi = 1
           kfl_sgsac_nsi = getint('ORDER',1_ip,'#Time integration order')
           if(exists('FIRST')) kfl_sgsac_nsi = 1
           if(exists('SECON')) kfl_sgsac_nsi = 2
           if(exists('PICAR')) kfl_sgsli_nsi = 1
           if(exists('NEWTO')) kfl_sgsli_nsi = 2
           if(exists('COUPL')) kfl_sgscp_nsi = 1
           if(exists('DECOU')) kfl_sgscp_nsi = 0
           if(exists('LINEA')) then ! Linear subscales model with convection tracking
              ! Very simple model that satisfies global momentum conservation
              kfl_sgsli_nsi = 0
              kfl_sgsco_nsi = 0 ! turned 0 because subscales are calculated in nsi_elmsgs when matrix assembly
           end if

        else if( words(1) == 'TIMES' ) then
           !
           ! Time step strategy: in development
           !
           if( words(2) == 'LOCAL' ) then
              kfl_dttyp_nsi = 1
           else if( words(2) == 'TAU  ' ) then
              kfl_dttyp_nsi = 3
           else if( words(2) == 'EIGEN' ) then
              kfl_dttyp_nsi = 4
           else if( words(2) == 'ADAPT' ) then
              kfl_dttyp_nsi = 2
              strec_nsi = getrea('STRET',2.0_rp,  '#Adaptive dt: Stretching factor')
              dampi_nsi = getrea('DAMPI',2.0_rp,  '#Adaptive dt: damping')
              epsht_nsi = getrea('EPSR ',0.025_rp,'#Adaptive dt: eps_R')
              epstr_nsi = getrea('EPSA ',0.025_rp,'#Adaptive dt: eps_A')
           end if


        else if( words(1) == 'DIRIC' ) then
           !
           !.md<1>DIRICHLET_CONDITION:  MATRIX | ELEMENT                                                    $ WHere Dirichlet bc are applied
           !.md<field>DIRICHLET_CONDITION
           !.md<com>Dirichlet boundary conditions can be applied at the element level, or directly
           !.md<com>on the matrix. The second option should be selected when using PERIODIC nodes.
           !
           if(      words(2) == 'MATRI' ) then
              kfl_matdi_nsi =  NSI_DIRICHLET_MATRIX
           else if( words(2) == 'ELEME' ) then
              kfl_matdi_nsi =  NSI_DIRICHLET_ELEMENT
           else if( words(2) == 'OFF  ' ) then
              kfl_matdi_nsi = -1
           else if( words(2) == 'ALGOR' ) then
              kfl_matdi_nsi =  NSI_DIRICHLET_MATRIX
              !kfl_matdi_nsi =  NSI_DIRICHLET_ALGORITHM ! This is obsolete
           end if

        else if( words(1) == 'INTER' ) then
           !
           !.md<1>INTERNAL_FORCES:      RESIDUAL | INTEGRAL                                                 $ How to compute forces on boundaries
           !.md<field>INTERNAL_FORCES
           !.md<com>The force of the fluid on the boundaries can be computed in two ways. On the one hand,
           !.md<com>by integrating the traction over the boundaries. On the other hand, by computing the nodal residual
           !.md<com>of the momentum equations where a Dirichlet condition is prescibed: f^i = buu^i - (Auu.u)^i - (Aup.p)^i
           !
           if( words(2) == 'RESID' ) then
              kfl_intfo_nsi = 1
           else if( words(2) == 'INTEG' ) then
              kfl_intfo_nsi = 0
           end if

           ! ************************************************
           !            kfl_intfo_nsi = 1_ip ! Added for FSI coupling test
           ! ************************************************

        else if( words(1) == 'ELEME' ) then
           !
           !.md<1>ELEMENT_LENGTH:       MINIMUM | MAXIMUM                                                   $ Element length h for stablization parameter
           !.md<field>ELEMENT_LENGTH
           !.md<com>Element length used to compute the stabilization parameters. MINIMUM leads to a less
           !.md<com>diffusive algorithm. MAXIMUM length is more diffusive, and more non-linear as the stabilization
           !.md<com>is a non-linear term.
           !
           call realen(kfl_ellen_nsi)

        else if( words(1) == 'ASSEM' ) then
           !
           ! Assembly type 
           !
           !.md<1>ASSEMBLY:             GPU2(FAST)/STABLE                                                   $ Assembly of NS equations
           !.md<field>ASSEMBLY
           !.md<com>Mainly three Navier-Stokes assemblies are available:
           !.md<com>    - FULL:       same assembly as the implicit method: ASGS technique available
           !.md<com>    - GPU2(FAST): assembly exclusive for fractional step methods: only Galerkin available
           !.md<com>    - STABLE:     assembly exclusive for fractional step method with reduced number of option. Galerkin
           !.md<com>                  and SUPG are available. Element length is taken as minimum.
           !
           if( exists('BOUND' ) ) then
              if(      words(2) == 'FAST ' ) then
                 kfl_asbou_nsi= 5
              else if( words(2) == 'GPU2 ' ) then        ! Usar mod_nsi_boundary_operations_fast 
                 kfl_asbou_nsi= 5
              else
                 call runend('NSI_REANUT: ASSEMBLY OPTION DEPRECTAED')                 
              end if
           else
              if(      words(2) == 'FAST ' ) then
                 kfl_assem_nsi= 5
              else if( words(2) == 'FULL ' ) then
                 kfl_assem_nsi= 0
              else if( words(2) == 'STABL' ) then
                 kfl_assem_nsi= 2
                 if( exists('SIMD ' ) )  kfl_assem_nsi = 6
              else if( words(2) == 'GPU2 ' ) then        ! Usar mod_nsi_element_operations_fast 
                 kfl_assem_nsi= 5
              else if( words(2) == 'GPU71' ) then        ! Usar mod_nsi_element_operations_hh71
                 kfl_assem_nsi= 71
              else if( words(2) == 'GPU80' ) then        ! Usar mod_nsi_element_operations_ 
                 kfl_assem_nsi= 80
              else if( words(2) == 'GPU90' ) then        ! Usar mod_nsi_element_operations_ 
                 kfl_assem_nsi= 90
              else if( words(2) == 'GPU91' ) then        ! Usar mod_nsi_element_operations_ 
                 kfl_assem_nsi= 91
              else
                 call runend('NSI_REANUT: ASSEMBLY OPTION DEPRECTAED')
              end if
           end if

        else if( words(1) == 'PENAL' ) then
           !
           ! Obsolete option
           !
           if(exists('ON   ')) kfl_penal_nsi = 1
           if(exists('CLASS')) kfl_penal_nsi = 1
           if( kfl_penal_nsi == 1 ) call runend('CLASSICAL PENALIZATION NOT CODED')
           if(exists('ITERA')) kfl_penal_nsi = 2
           if(exists('ARTIF')) kfl_penal_nsi = 3
           if(kfl_penal_nsi>=1 ) then
              penal_nsi=getrea('VALUE',0.0_rp,'#Penalty parameter')
              if(exists('AUTOM')) kfl_penal_nsi = 4
           end if

        else if( words(1) == 'PREPE' ) then
           !
           !.md<1>PRE_PENALIZATION:     VALUE=real                                                          $ Pressure penalization
           !.md<field>PRE_PENALIZATION
           !.md<com>Penalize pressure in time. To be used when converging to a steady state when
           !.md<com>the pressure presents an odd-even decoupling in time. This option helps converging
           !.md<com>by adding a damping effect.
           !
           kfl_prepe_nsi = 1
           prepe_nsi = 1.0_rp
           if(kfl_prepe_nsi>=1 ) then
              prepe_nsi=getrea('VALUE',1.0_rp,'#Pressure Penalty parameter')
           end if

        else if( words(1) == 'SHOCK' ) then
           !
           !.md<1>SHOCK_CAPTURING:      ISOTROPIC | ANISOTROPIC, VALUE=real                                 $ Shock capturing for the momentum equations
           !.md<field>SHOCK_CAPTURING
           !.md<com>Shock capturing for the momentum equations. It enables to obtain
           !.md<com>a "monotonic" scheme. The algorithm can be ISOTROPIC (very diffusive)
           !.md<com>or ANISOTROPIC (by adding diffusion only in the crosswind direction).
           !.md<com>The added term can be weighted by the VALUE=real. 0 means nothing added
           !.md<com>and 1 means all the term is added.
           !
           if(exists('ISOTR').or.exists('ON   ') ) then
              kfl_shock_nsi = 1
              shock_nsi     = getrea('VALUE',0.0_rp,'#Shock capturing parameter')
           else if(exists('ANISO') ) then
              kfl_shock_nsi = 2
              shock_nsi     = getrea('VALUE',0.0_rp,'#Shock capturing parameter')
           end if

        else if( words(1) == 'UPDAT' ) then
           !
           !.md<1>UPDATE_PRESSURE:      ON                                                                  $ Pressure extrapolation at end of time step
           !.md<field>UPDATE_PRESSURE
           !.md<com>If this option is activated, the pressure is extrapolated at the end of
           !.md<com>a time step: p^{n+1} = 2 p^{n+1} - p^n
           !
           if( words(2) == 'ON   ') kfl_updpr_nsi=1

        else if( words(1) == 'EULER' ) then
           !
           ! Retro-compatibility
           !
           neule_nsi=getint('EULER',0_ip,'#EULER TIME STEPS')

        else if( words(1) == 'TIMEI' ) then
           !
           !.md<1>TIME_INTEGRATION:     TRAPEZOIDAL | BDF , ORDER=int1 , EULER=int2
           !.md<field>TIME_INTEGRATION
           !.md<com>Time integration scheme and options. Two schemes are available: trapezoidal rule
           !.md<com>of order 1 or 2, and backward finite difference BDF of order from 1 to 6. For order 2, BDF
           !.md<com>is generally more stable than the trapezoidal rule. In addition, for the sake of time stability,
           !.md<com>one can carry out some few EULER=int2 iterations before triggering a higher order method.
           !
           if(exists('TRAPE')) kfl_tisch_nsi=1
           if(exists('BDF  ')) kfl_tisch_nsi=2
           if(exists('ADAMS')) kfl_tisch_nsi=3 ! Adams-Bashforth explicit
           if(exists('RUNGE')) kfl_tisch_nsi=4 ! Runge-Kutta explicit
           kfl_tiacc_nsi = getint('ORDER',1_ip,'#Time integration order')
           neule_nsi     = getint('EULER',0_ip,'#EULER TIME STEPS')

        else if( words(1) == 'SAFET' ) then
           !
           !.md<1>SAFETY_FACTOR=        real
           !.md<field>SAFETY_FACTOR
           !.md<com>When using an explicit scheme, the CFL condition determines a maximum time step
           !.md<com>to obtain a stable scheme. When using
           !
           safet_nsi = param(1)
           safeo_nsi = safet_nsi ! initial safety factor

           safex_nsi  = getrea('EXPON',1.0_rp,'#Safety Factor Exponential time function')

           safma_nsi  = getrea('MAXIM',1.0e9_rp,'#Maximum safety factor function')

           saflo_nsi  = getrea('MINGL',safet_nsi,'#Minimum global safety factor when using local time step. Fixes minimum time step over elements')
           !
           ! if uniform safety factor
           !
           if (abs(safex_nsi-1.0_rp).lt.1.0e-7_rp) safma_nsi=safeo_nsi
        else if( words(1) == 'BEMOL' ) then
           !
           !.md<1>BEMOL=                real                                                                $ Put 1 to integrate convection by parts
           !.md<field>BEMOL
           !.md<com>When BEMOL=0, the convection term is not integrated by parts. When BEMOL=1, it is fully
           !.md<com>integrated by parts. Intermediate values are also possible. BEMOL=1 can have some
           !.md<com>robustness advantages in some cases.
           !
           bemol_nsi = param(1)

        else if( words(1) == 'STEAD' ) then
           !
           !.md<1>STEADY_STATE_TOLERANCE=real                                                               $ Tolerance for steady state
           !.md<field>STEADY_STATE_TOLERANCE
           !.md<com>Tolerance to decide when the solution of the module should stop because the steady
           !.md<com>state have been reached.
           !
           sstol_nsi = param(1)

        else if( words(1) == 'MODIF' ) then
           !
           ! Users defined boundary conditions
           !
           if( words(2) == 'FREE ' ) then
              kfl_modfi_nsi = 1
              xfree_nsi = getrea('XFREE',0.0_rp,'#X Coordinate of the plane where to free')
           end if
           if( words(2) == 'PRNAC') kfl_modfi_nsi=2
           if( words(2) == 'PRNA3') kfl_modfi_nsi=3
           if( words(2) == 'PRNA4') kfl_modfi_nsi=4
           if( words(2) == 'PRNA5') kfl_modfi_nsi=5
           if( words(2) == 'PRNA6') kfl_modfi_nsi=6

        else if( words(1) == 'NORMO' ) then
           !
           ! Almost obsolete
           !
           if(exists('BACKW').or.exists('LAGGE').or.exists('ALGEB') ) then
              kfl_normc_nsi = 3
           else if(exists('L2   ') ) then
              kfl_normc_nsi = 2
           else if(exists('L1   ') ) then
              kfl_normc_nsi = 1
           else if(exists('L-inf') ) then
              kfl_normc_nsi = 0
           else if(exists('REFER') ) then
              kfl_normc_nsi = 4
              if( exists('FIELD') ) then
                 kfl_refer_nsi=getint('FIELD',1_ip,'#Field of the reference solution')
              else
                 call runend('NSI_REAOUS: REFERENCE SOLUTION NOT CODED')
              end if

           end if
        else if( words(1) == 'GRADD' ) then
           !
           ! Immersed boundary method
           !
           if( exists('GRADD') ) kfl_grad_div_nsi = 1

        else if( words(1) == 'IMMER' ) then
           !
           ! Immersed boundary method
           !
           if( exists('ON   ') ) kfl_immer_nsi = 1
           if( exists('VELOC') ) kfl_immer_nsi = 1
           if( exists('DEFOR') ) kfl_immer_nsi = 2

        else if( words(1) == 'MAXIM' ) then
           !-->                <inputLine>
           !-->                    <inputLineName>MAXIMUM_ITERATIONS</inputLineName>
           !-->                    <inputLineHelp>MAXIMUM_ITERATIONS:
           !-->                        Maximum number of inner iterations (non-linearity+Orthomin(1)). START_AT=int2 enables on to use only one
           !-->                        iterations until time step int2 is reached. This can be used to go fast when
           !-->                        flow is establishing. The inner iterations are also controlled by the option CONVERGENCE_TOLERANCE.</inputLineHelp>
           !-->                    <inputElement>
           !-->                        <inputElementType>edit2</inputElementType>
           !-->                <inputElementValueType>INT</inputElementValueType>
           !-->                        <inputLineEditValue>5</inputLineEditValue>
           !-->                    </inputElement>
           !-->                    <inputElement>
           !-->                        <inputElementType>edit</inputElementType>
           !-->                <inputElementValueType>INT</inputElementValueType>
           !-->                        <inputLineEditName>START_AT</inputLineEditName>
           !-->                        <inputLineEditValue>5</inputLineEditValue>
           !-->                    </inputElement>
           !-->                </inputLine>
           !
           !.md<1>MAXIMUM_ITERATIONS=   int1 [START_AT=int2]                                                $ Max number of inner iterations
           !.md<field>MAXIMUM_ITERATIONS
           !.md<com>Maximum number of inner iterations (non-linearity+Orthomin(1)). START_AT=int2 enables on to use only one
           !.md<com>iterations until time step int2 is reached. This can be used to go fast when
           !.md<com>flow is establishing. The inner iterations are also controlled by the option CONVERGENCE_TOLERANCE.
           !
           momod(modul) % miinn = int(param(1))
           if( exists('START') ) kfl_stain_nsi = getint('START',0_ip,'#Time step to start')

        else if( words(1) == 'CONVE' ) then
           !-->                <inputLine>
           !-->                    <inputLineName>CONVERGENCE_TOLERANCE</inputLineName>
           !-->                    <inputLineHelp>CONVERGENCE_TOLERANCE:
           !-->                       Convergence tolerance to control the inner iterations (non-linearity+Orthomin(1)). The
           !-->                       residual is based on the maximum lagged algebraic residual of the momentum and continuity
           !-->                       equations. The inner iterations are also controlled by the option MAXIMUM_ITERATIONS.
           !-->                    </inputLineHelp>
           !-->                    <inputElement>
           !-->                        <inputElementType>edit2</inputElementType>
           !-->                <inputElementValueType>REAL</inputElementValueType>
           !-->                        <inputLineEditValue>5</inputLineEditValue>
           !-->                    </inputElement>
           !-->                </inputLine>
           !
           !.md<1>CONVERGENCE_TOLERANCE=real                                                                $ Tolerance for inner iterations
           !.md<field>CONVERGENCE_TOLERANCE
           !.md<com>Convergence tolerance to control the inner iterations (non-linearity+Orthomin(1)). The
           !.md<com>residual is based on the maximum lagged algebraic residual of the momentum and continuity
           !.md<com>equations. The inner iterations are also controlled by the option MAXIMUM_ITERATIONS.
           !
           cotol_nsi = param(1)

        else if( words(1) == 'LINEA' ) then
           !-->                <inputLine>
           !-->                    <inputLineName>LINEARIZATION</inputLineName>
           !-->                    <inputLineHelp>LINEARIZATION:
           !-->                        Linearization strategy for the convective term: PICARD or NEWTON. When NEWTON is
           !-->                        selected, an additional option PICARD=int enables to carry out some int Picard
           !-->                        iterations at each inner iterations in order to have a more robust convergence. Tipically,
           !-->                        one carries out 1 or 2 Picard iterations. </inputLineHelp>
           !-->                    <inputElement>
           !-->                        <inputElementType>combo</inputElementType>
           !-->                        <item><itemName>PICARD</itemName></item><item><itemDependence>0</itemDependence><itemName>NEWTON</itemName></item>
           !-->                    </inputElement>
           !-->                    <inputElement>
           !-->                        <inputElementGroup>0</inputElementGroup>
           !-->                        <inputElementType>edit</inputElementType>
           !-->                <inputElementValueType>INT</inputElementValueType>
           !-->                        <inputLineEditName>PICARD</inputLineEditName>
           !-->                        <inputLineEditValue>5</inputLineEditValue>
           !-->                    </inputElement>
           !-->                </inputLine>
           !
           !.md<1>LINEARIZATION:        PICARD | NEWTON[PICARD=int]                                         $ Convective term linearization
           !.md<field>LINEARIZATION
           !.md<com>Linearization strategy for the convective term: PICARD or NEWTON. When NEWTON is
           !.md<com>selected, an additional option PICARD=int enables to carry out some int Picard
           !.md<com>iterations at each inner iterations in order to have a more robust convergence. Tipically,
           !.md<com>one carries out 1 or 2 Picard iterations.
           !
           if( words(2) == 'PICAR' ) then
              kfl_linea_nsi = 1
           else if( words(2) == 'NEWTO' ) then
              kfl_linea_nsi = 2
              npica_nsi     = getint('PICAR',0_ip,'#Number of Picard iterations')
           end if

        else if( words(1) == 'SAVEC' ) then
           !
           ! Obsolete option
           !
           if( words(2) == 'YES  ') kfl_savco_nsi=1

        else if( words(1) == 'MOMEN' ) then
           !
           ! Current problem is momentum equations (used to define solver)
           !
           ivari_nsi = ivari_nsi_mom
           imome     = 1

        else if( words(1) == 'MASSC' ) then
           !
           ! Current problem is mass correction  (used to define solver - needed with consistent mass matrix)
           !
           ivari_nsi = ivari_nsi_corre
           if (kfl_corre_nsi /= 3) call runend('nsi_neanut: setting solver information if the Correction is not Consistent will lead to an error')

        else if( words(1) == 'DIVER' ) then
           !
           ! Current problem is divergence free correction equation (used to define solver)
           !
           ivari_nsi = ivari_nsi_divfree

        else if( words(1) == 'NORMA' ) then
           !
           ! Current problem is normal extension equation (used to define solver)
           !
           ivari_nsi = 7

        else if( words(1) == 'VISCO' ) then
           !
           ! Current problem is normal extension equation (used to define solver)
           !
           ivari_nsi = 9

        else if( words(1) == 'ENDMO' ) then
           !
           ! Current problem is no longer momentum equations (used to define solver)
           !
           imome     = 0

        else if( words(1) == 'CONTI' ) then
           !
           ! Current problem is continuity equation (used to define solver)
           !
           ivari_nsi = ivari_nsi_cont
           icont     = 1

        else if( words(1) == 'ENDCO' ) then
           !
           ! Current problem is no longer continuity equation (used to define solver)
           !
           icont     = 0

        else if( words(1) == 'ENDMA' ) then
           !
           ! Current problem is no longer mass correction (used to define solver)
           !
           ivari_nsi = ivari_nsi_mom

        else if( words(1) == 'ENDHY' ) then
           !
           ! Current problem is no longer hydrostatic state equation (used to define solver)
           !
           ivari_nsi = ivari_nsi_mom

        else if( words(1) == 'CONSI' ) then
           !
           ! Consistent matrix
           !
           ivari_nsi = 8

        else if( words(1) == 'HYDRO' ) then
           !
           ! Current problem is hydrostatic state equation (used to define solver)
           !
           do while( words(1) /= 'ENDHY' )

              if( words(1) == 'METHO' ) then

                 if( words(2) == 'ANALY' ) then
                    kfl_hydro_nsi = NSI_ANALYTICAL_HYDROSTATIC_PRESSURE
                 else if( words(2) == 'PDE  ' ) then
                    kfl_hydro_nsi = NSI_PDE_HYDROSTATIC_PRESSURE
                 else
                    call runend('NON EXISTING SOLUTION METHOD FOR HYDROSTATIC PRESSURE')
                 end if

              else if( words(1) == 'ALGEB' ) then
                 ivari_nsi =  ivari_nsi_hydro
                 solve_sol => solve(ivari_nsi:)
                 call reasol(1_ip)

              else if( words(1) == 'UPDAT' ) then

                 if( words(2) == 'ON   ') then
                    kfl_update_hydro_nsi = ITASK_BEGITE
                 else if( words(2) == 'OFF  ') then
                    kfl_update_hydro_nsi = ITASK_INIUNK
                 end if

              else if( words(1) == 'INTER' ) then

                 if( words(2) == 'CONST') then
                    kfl_hydro_interface_nsi = 0
                    heihy_nsi = getrea('VALUE',0.0_rp,'#Height of free surface')
                 else
                    kfl_hydro_interface_nsi = 1
                 end if

              end if
              call ecoute('NSI_REANUT')
           end do

        else if( words(1) == 'INTEG' ) then
           !
           ! Obsolete option
           !
           if( words(2) == 'OFF  ' ) then
              kfl_intpr_nsi=0
           else
              kfl_intpr_nsi=1
           end if
           !-->            </subGroup>
        else if( words(1) == 'ALGEB' ) then
           !-->            <subGroup>
           !-->                <subGroupName>MOMENTUM</subGroupName>
           !-->                <subGroupCheckeable>true</subGroupCheckeable>
           !-->                <subGroupType>ALGEBRAIC_SOLVER</subGroupType>
           !-->            </subGroup>
           !-->            <subGroup>
           !-->                <subGroupName>CONTINUITY</subGroupName>
           !-->                <subGroupCheckeable>true</subGroupCheckeable>
           !-->                <subGroupType>ALGEBRAIC_SOLVER</subGroupType>
           !-->            </subGroup>
           !-->            <subGroup>
           !-->                <subGroupName>HYDROSTATIC_STATE</subGroupName>
           !-->                <subGroupCheckeable>true</subGroupCheckeable>
           !-->                <subGroupType>ALGEBRAIC_SOLVER</subGroupType>
           !-->            </subGroup>
           !
           !.md<1>XXX_EQUATION                                                                              $ Define solver for equation XXX
           !.md<2>  ALGEBRAIC_SOLVER
           !.md<2>    INCLUDE ./sample-solver.dat
           !.md<2>  END_ALGEBRAIC_SOLVER
           !.md<1>END_XXX_EQUATION
           !.md<field>XXX_EQUATION
           !.md<com>Define the algebraic solver options for equation XXX.
           !.md<com>XXX: MOMENTUM | CONTINUITY | HYDROSTATIC_STATE
           !
           solve_sol => solve(ivari_nsi:)
           call reasol(1_ip)

        else if( words(1) == 'PRECO' ) then
           !
           ! See previous option
           !
           solve_sol => solve(ivari_nsi:)
           call reasol(2_ip)

        else if( words(1) == 'PRESS' ) then
           !
           ! Pressure integration by parts
           !
           if( option('PRESS') ) then
              kfl_press_nsi = 1
           else
              kfl_press_nsi = 0
           end if

        else if( words(1) == 'ENRIC' ) then
           !
           ! Pressure enrichement
           !
           if(      words(2) == 'FLAT ' ) then
              kfl_bubbl_nsi = FLAT_BUBBLE
           else if( words(2) == 'QUADR' ) then
              kfl_bubbl_nsi = QUADRATIC_BUBBLE
           else if( words(2) == 'FREES' ) then
              kfl_bubbl_nsi = FREE_SURFACE_BUBBLE
           end if

        else if( words(1) == 'FSGRA' ) then
           !
           ! Fractional Step with Gravity or Boussinesq correction
           !    
           kfl_fsgrb_nsi = 1

        else if( words(1) == 'STABC' ) then
           !
           ! Coriolis Stabilization
           !    
           kfl_stab_corio_nsi = 1_ip
           staco_corio_nsi = param(1)
           
        else if( words(1) == 'RELAX' ) then
           !
           ! Obsolete option: only available for monolithic scheme
           !
           relax_nsi            = param(1)
           relap_nsi            = param(1)

           if(exists('CONST') ) then
              if(exists('VELOC') ) then
                 kfl_relax_nsi=1
              end if
              if(exists('PRESS') ) then
                 kfl_relap_nsi=1
              end if
              if(exists('BOTH ') ) then
                 kfl_relax_nsi=1
                 kfl_relap_nsi=1
              end if

           else if(exists('AITKE') ) then
              if(exists('VELOC') ) then
                 kfl_relax_nsi=2
              end if
              if(exists('PRESS') ) then
                 kfl_relap_nsi=2
              end if
              if(exists('BOTH ') ) then
                 kfl_relax_nsi=2
                 kfl_relap_nsi=2
              end if
           end if

        else if( words(1) == 'NOTA2' ) then
           !
           ! This option is being tested: do not stabilize continuity
           !
           kfl_nota2 = 1

        else if( words(1) == 'TA2MA' ) then
           !
           ! This option is being tested: do not stabilize continuity
           !
           kfl_ta2ma = 1

        else if( words(1) == 'EXPCO' ) then
           !
           ! This option is being tested: Treat the convective term explicitly,
           ! that is, assemble the matrix only in the first ortomin iteration
           !
           kfl_expco_nsi = 1

        else if( words(1) == 'INITI' ) then   ! INITIAL_GUESS_ORDER = 1 or 2
           !
           ! Order of the initial guess extrapolation at each time step  - for the moment only velocity.
           !
           kfl_ini_ts_guess_order_nsi    = nint(param(1))

        else if( words(1) == 'ADDPR' ) then
           !
           ! This option is being tested
           !
           kfl_addpr_nsi = 1

        else if( words(1) == 'STOPB' ) then
           !
           ! Stop by convergence of witness points
           !
           kfl_stop_by_wit_nsi = 1

        else if( words(1) == 'MESHI' ) then
           !--><inputLine>
           !-->            <inputLineName>MESH_INTERPOLATION</inputLineName>
           !-->            <inputLineHelp><![CDATA[MESH_INTERPOLATION:
           !-->                        This tool allows you to interpolate a solution from an original grid or a particular setup
           !-->                        to a different mesh or setup. It requires the correct set of the coupling files and allows
           !-->                        to use different number of processors.
           !-->            <inputElement>
           !-->                <inputElementType>combo</inputElementType>
           !-->                <item><itemName>MINIMUM</itemName></item><item><itemName>MAXIMUM</itemName></item>
           !-->            </inputElement>
           !--> </inputLine>
           !
           !.md<1>MESH_INTERPOLATION:       OFF | ON
           !.md<field>MESH_INTERPOLATION
           !.md<com>OFF: no interpolation, ON: interpolation activation.
           !
           if(exists('ON   ') ) kfl_meshi_nsi = 1

        else if( words(1) == 'ALGOR' ) then
           !-->            <subGroup>
           !-->                <subGroupName>ALGORITHM</subGroupName>
           !-->                <!--este valor se imprime cuando se crea el fichero de alya así: ALGORITHM: SCHUR_COMPLEMENT-->
           !-->                <subGroupValue>SCHUR_COMPLEMENT</subGroupValue>
           !-->                <inputLine>
           !-->                    <inputLineName>Solver</inputLineName>
           !-->                    <inputLineHelp><![CDATA[ALGORITHM: SCHUR_COMPLEMENT
           !-->                      Only the Schur complement based solver is available. This field defines
           !-->                      all the options of the solver. Four solvers are available:
           !-->                      
           !-->                             -  ORTHOMIN,   MOMENTUM_PRESERVING;
           !-->                             -  ORTHOMIN,   CONTINUITY_PRESERVING;
           !-->                             -  RICHARDSON, MOMENTUM_PRESERVING;
           !-->                             -  RICHARDSON, CONTINUITY_PRESERVING.
           !-->                      ]]></inputLineHelp>
           !-->                    <inputElement>
           !-->                        <inputElementType>combo</inputElementType>
           !-->                        <item>
           !-->                            <itemName>ORTHOMIN</itemName>
           !-->                        </item>
           !-->                        <item>
           !-->                            <itemName>RICHARDSON</itemName>
           !-->                        </item>
           !-->                    </inputElement>
           !-->                    <inputElement>
           !-->                        <inputElementType>combo</inputElementType>
           !-->                        <item>
           !-->                            <itemName>MOMENTUM_PRESERVING</itemName>
           !-->                        </item>
           !-->                        <item>
           !-->                            <itemName>CONTINUITY_PRESERVING</itemName>
           !-->                        </item>
           !-->                    </inputElement>
           !-->                </inputLine>
           !-->                <inputLine>
           !-->                    <inputLineName>Preconditioner</inputLineName>
           !-->                    <inputLineHelp><![CDATA[Additional information is required. The definition of the pressure Schur complement
           !-->                       preconditioner, and in particular, the approximation of the Uzawa operator. Two options
           !-->                       are available: DT or TAU.]]></inputLineHelp>
           !-->                    <inputElement>
           !-->                        <inputElementType>combo</inputElementType>
           !-->                        <item>
           !-->                            <itemName>DT</itemName>
           !-->                        </item>
           !-->                        <item>
           !-->                            <itemName>TAU</itemName>
           !-->                        </item>
           !-->                    </inputElement>
           !-->                </inputLine>
           !-->                <inputLine>
           !-->                    <inputLineName>TAU_STRATEGY</inputLineName>
           !-->                    <inputLineHelp> If TAU is selected one then have to choose if TAU is approximated
           !-->                       using the CODINA, SHAKIB or EXACT solution based tau.</inputLineHelp>
           !-->                    <inputElement>
           !-->                        <inputElementType>combo</inputElementType>
           !-->                        <item>
           !-->                            <itemName>CODINA</itemName>
           !-->                        </item>
           !-->                        <item>
           !-->                            <itemName>SHAKIB</itemName>
           !-->                        </item>
           !-->                        <item>
           !-->                            <itemName>EXACT</itemName>
           !-->                        </item>
           !-->                    </inputElement>
           !-->                </inputLine>
           !-->                <inputLine>
           !-->                    <inputLineName>ELEMENT_LENGTH</inputLineName>
           !-->                    <inputLineHelp>These expressions require also
           !-->                    the definition of the element length: MINIMUM or MAXIMUM.</inputLineHelp>
           !-->                    <inputElement>
           !-->                        <inputElementType>combo</inputElementType>
           !-->                        <item>
           !-->                            <itemName>MINIMUM</itemName>
           !-->                        </item>
           !-->                        <item>
           !-->                            <itemName>MAXIMUM</itemName>
           !-->                        </item>
           !-->                    </inputElement>
           !-->                </inputLine>
           !-->                <inputLine>
           !-->                    <inputLineName>CORRECTION</inputLineName>
           !-->                    <inputLineHelp><![CDATA[For CONTINUITY_PRESERVING version
           !-->                         one has to finally define the rule to define the correction matrix.
           !-->                         Usually, the following set of options exhibits good convergence:
           !-->                    
           !-->                           -  ORTHOMIN, MOMENTUM_PRESERVING; TAU; CODINA; MINIMUM
           !-->                           -  ORTHOMIN, CONTINUITY_PRESERVING; DT; CODINA; MINIMUM
           !-->                    
           !-->                    The first set requires one less solve of the continuity equation, but the second
           !-->                    one is usually more robust.]]></inputLineHelp>
           !-->                    <inputElement>
           !-->                        <inputElementType>combo</inputElementType>
           !-->                        <item>
           !-->                            <itemName>OFF</itemName>
           !-->                        </item>
           !-->                        <item>
           !-->                            <itemName>CLOSE</itemName>
           !-->                        </item>
           !-->                        <item>
           !-->                            <itemName>OPEN</itemName>
           !-->                        </item>
           !-->                    </inputElement>
           !-->                </inputLine>
           !-->            </subGroup>
           !
           !.md<1>ALGORITHM: SCHUR_COMPLEMENT
           !.md<2>  SOLVER:         ORTHOMIN | RICHARDSON , MOMENTUM_PRESERVING | CONTINUITY_PRESERVING     $ Schur complement solver
           !.md<2>  PRECONDITIONER: DT | TAU                                                                $ Dt or tau for preconditioner
           !.md<2>  TAU_STRATEGY:   CODINA | SHAKIB | EXACT                                                 $ tau approximation
           !.md<2>  ELEMENT_LENGTH: MINIMUM | MAXIMUM                                                       $ h in tau
           !.md<2>  CORRECTION:     OFF | CLOSE | OPEN                                                      $ Rule for mass correction in continuity preserving
           !.md<1>END_ALGORITHM
           !.md<field>ALGORITHM: SCHUR_COMPLEMENT
           !.md<com>Only the Schur complement based solver is available. This field defines
           !.md<com>all the options of the solver. Four solvers are available:
           !.md<com>    -  ORTHOMIN,   MOMENTUM_PRESERVING;
           !.md<com>    -  ORTHOMIN,   CONTINUITY_PRESERVING;
           !.md<com>    -  RICHARDSON, MOMENTUM_PRESERVING;
           !.md<com>    -  RICHARDSON, CONTINUITY_PRESERVING.
           !.md<com>
           !.md<com>Additional information is required. The definition of the pressure Schur complement
           !.md<com>preconditioner, and in particular, the approximation of the Uzawa operator. Two options
           !.md<com>are available: DT or TAU. If TAU is selected one then have to choose if TAU is approximated
           !.md<com>using the CODINA, SHAKIB or EXACT solution based tau. These expressions require also
           !.md<com>the definition of the element length: MINIMUM or MAXIMUM. For CONTINUITY_PRESERVING version
           !.md<com>one has to finally define the rule to define the correction matrix.
           !.md<com>
           !.md<com>Usually, the following set of options exhibits good convergence:
           !.md<com>    -  ORTHOMIN, MOMENTUM_PRESERVING; TAU; CODINA; MINIMUM
           !.md<com>    -  ORTHOMIN, CONTINUITY_PRESERVING; DT; CODINA; MINIMUM
           !.md<com>The first set requires one less solve of the continuity equation, but the second
           !.md<com>one is usually more robust.
           !
           if(exists('MONOL')) then
              !
              ! Monolithic
              !
              kfl_algor_nsi = 1

           else if(exists('SEMII')) then
              !
              ! Semi implicit method
              !
              kfl_algor_nsi = 3
              kfl_predi_nsi = 9
              call ecoute('nsi_reanut')
              do while(words(1)/='ENDAL')
              end do

           else if(exists('FRACT')) then
              !
              ! Fractional step
              !
              kfl_algor_nsi = 2 ! Fraction step method
              kfl_predi_nsi = 7 ! Pure Laplacian for Q

              call ecoute('nsi_reanut')
              do while(words(1)/='ENDAL')

                 if( words(1) == 'MASSM' ) then
                    !
                    ! Mass matrix
                    !
                    if( words(2) == 'LUMPE' .or. words(2) == 'DIAGO' ) then
                       kfl_massm_nsi = NSI_LUMPED_MASS
                    else if( words(2) == 'CONSI' ) then
                       kfl_massm_nsi = NSI_CONSISTENT_MASS
                    end if

                 else if( words(1) == 'GAMMA' ) then
                    !
                    ! Gamma: shoudl be =1 if ASGS is used, and 0 for Galerking
                    !
                    gamma_nsi = getrea('GAMMA',1.0_rp,'#Gamma for pressure')
                 else if( words(1) == 'ROTAT' ) then

                    fsrot_nsi = getrea('ROTAT',1.0_rp,'#Use the rotational form')

                 else if( words(1) == 'CONSI' ) then

                    kfl_fscon_nsi = 1_ip

                 else if( words(1) == 'PRESS' ) then
                    !
                    ! Pressure stabilization
                    !
                    if(      words(2) == 'DT   ' ) then
                       kfl_press_stab_nsi = 0
                    else if( words(2) == 'TAU  ' ) then
                       if( exists('ORTHO') ) then
                          kfl_press_stab_nsi = 2
                       else if( exists('CG   ') ) then
                          kfl_press_stab_nsi = 3
                       else if( exists('BICGS') ) then
                          kfl_press_stab_nsi = 4
                       else
                          kfl_press_stab_nsi = 1
                       end if
                    else
                       call runend('NSI_REANUT: UNKNOWN PRESSURE STABILIZATION')
                    end if

                 else if( words(1) == 'SOLVE' ) then
                    !
                    ! Max number of iterations used with Tau pressure stabilization
                    !
                    if(      words(2) == 'ORTHO' ) then
                       kfl_press_stab_nsi = 2
                    else if( words(2) == 'CG   ' ) then
                       kfl_press_stab_nsi = 3
                    else if( words(2) == 'BICGS' ) then
                       kfl_press_stab_nsi = 4
                    else
                       kfl_press_stab_nsi = 1
                    end if

                 else if( words(1) == 'ITERA' ) then
                    !
                    ! Max number of iterations used with Tau pressure stabilization
                    !
                    mitri_nsi = getint('ITERA',1_ip,'#Number of Richardson iterations')

                 else if( words(1) == 'TOLER' ) then
                    !
                    ! Tolerance used with Tau pressure stabilization
                    !
                    toler_nsi = getrea('TOLER',1.0e-6_rp,'#Number of Tau solver iterations')

                    if( exists('ADAPT') ) then
                       !
                       ! Tolerance used with Tau pressure stabilization
                       !
                       kfl_adres_nsi = 1
                       if(exists('RATIO')) adres_nsi = getrea('RATIO',0.1_rp, '#Ratio wr to initial residual') ! Adpative residual
                    end if

                 else if( words(1) == 'PRECO' ) then
                    !
                    ! Preconditioner
                    !
                    if(exists('DT   ')) kfl_predi_nsi = 2
                    if(exists('TAU  ')) kfl_predi_nsi = 3
                    if(exists('MASS ')) kfl_predi_nsi = 4
                    if(exists('LUMPE')) kfl_predi_nsi = 5
                    if(exists('LAPLA')) kfl_predi_nsi = 7

                 end if

                 call ecoute('nsi_reanut')
              end do

           else if(exists('SCHUR') ) then
              !
              ! Schur complement
              !
              kfl_algor_nsi = 5
              kfl_taush_nsi = kfl_taust_nsi  ! Tau strategy for preconditioner
              kfl_predi_nsi = 3              ! Preconditioner=tau
              kfl_ellsh_nsi = 0              ! Minimum length
              kfl_corre_nsi = 2              ! Open correction
              kfl_sosch_nsi = 2              ! Momentum preserving Orthomin
              mitri_nsi     = 1              ! Number of Richardson iterations
              call ecoute('nsi_reanut')
              do while(words(1)/='ENDAL')
                 if( words(1) == 'SOLVE' ) then
                    if( words(2) == 'ORTHO' ) then
                       kfl_sosch_nsi = 1
                    else if( words(2) == 'RICHA' ) then
                       kfl_sosch_nsi = 4
                    else
                       call runend('WRONG SCHUR COMPLEMENT SOLVER')
                    end if
                    if(exists('MOMEN') ) then
                       kfl_sosch_nsi=kfl_sosch_nsi+1
                    else if(exists('CONTI') ) then
                       kfl_sosch_nsi=kfl_sosch_nsi+2
                    end if

                 else if( words(1) == 'RICHA' ) then
                    mitri_nsi = getint('RICHA',1_ip,'#Number of Richardson iterations')

                 else if( words(1) == 'INCRE' ) then
                    if( option('INCRE') ) then
                       kfl_incre_nsi = 1
                    else
                       kfl_incre_nsi = 0
                    end if

                 else if( words(1) == 'PRECO' ) then
                    if(exists('DT   ')) kfl_predi_nsi = 2 ! \int_V Tau grad(p).grad(q) dV
                    if(exists('TAU  ')) kfl_predi_nsi = 3 ! \int_V Tau grad(p).grad(q) dV
                    if(exists('MASS ')) kfl_predi_nsi = 4
                    if(exists('LUMPE')) kfl_predi_nsi = 5
                    if(exists('LAPLA')) kfl_predi_nsi = 7

                    if(exists('DIAGT')) kfl_predi_nsi = 8 ! Diagonally scaled Laplacian with Tau Q = Tim L
                    if(exists('DIAGD')) kfl_predi_nsi = 9 ! Diagonally scaled with dt  Q = dt/rho L

                 else if( words(1) == 'TAUST' ) then
                    call reatau(kfl_taush_nsi)

                 else if( words(1) == 'ELEME' ) then
                    call realen(kfl_ellsh_nsi)

                 else if( words(1) == 'CORRE' ) then
                    if(exists('OFF  '))                     kfl_corre_nsi = 0
                    if(exists('CLOSE'))                     kfl_corre_nsi = 1
                    if(exists('OPENR').or.exists('OPEN ') ) kfl_corre_nsi = 2
                    if(exists('CONSI'))                     kfl_corre_nsi = 3  ! Use consistent(non-diagonal) mass matrix
                 end if
                 call ecoute('nsi_reanut')
              end do

           end if

        end if
     end do
     !-->        </group>
     !
     !.md<0>END_NUMERICAL_TREATMENT
     !.md</code>
     !
     if ( kfl_nota2 /= 0 .and. kfl_ta2ma /=0 ) call runend ('NSI_REANUT: NOTA2 and TA2MA must not be used at the same time')
     if ( kfl_nota2 == 1 ) staco_nsi(4) = 0.0_rp
     if ( kfl_ta2ma == 1 ) staco_nsi(4) = 0.25_rp / staco_nsi(1)
     if ( kfl_nota1_nsi ==1 ) then ! if no stabilization of momentum equation, then subscales are static
        kfl_sgsco_nsi = 0
        kfl_sgsac_nsi = -1
        kfl_sgsti_nsi = 0
     end if
  end if

end subroutine nsi_reanut

