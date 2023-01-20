!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_inivar.f90
!> @author  Mariano Vazquez
!> @date    August, 2006
!>          - Subroutine written
!> @brief   Initialize some variables
!>
!> @details
!>          \verbatim
!>          Initialize some variables:
!>          ITASK = 0 ... When starting the run (from Turnon)
!>                  1 ... When starting the run (from Turnon)
!>                  2 ... Correct fibers when needed
!>          \endverbatim
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_inivar(itask)

  use def_parame,         only : pi
  use def_kintyp,         only : ip,rp
  use def_master,         only : postp,intost,retost
  use def_master,         only : kfl_timei,coupling,solve,ITASK_INITIA
  use def_master,         only : displ, INOTSLAVE, TIME_N, INOTMASTER
  use def_master,         only : times
  use def_master,         only : donna_gp
  use mod_moduls_conf,    only : moduls_allocate_timing
  use def_domain,         only : xfiel, lmate, ndime, nmate, nelem
  use def_solver,         only : SOL_LOCAL_MASS_MATRIX,SOL_SOLVER_RICHARDSON
  use mod_communications, only : PAR_MAX
#if defined(_OPENMP)
  use mod_parall,         only : par_omp_nelem_chunk
#endif
  use mod_arrays,         only : arrays_register
  use def_solidz,         only : lmate_sld, nmate_sld, nvoig_sld, celen_sld
  use def_solidz,         only : kfl_fixbo_sld, kfl_fixno_sld, kfl_fixrs_sld
  use def_solidz,         only : kfl_funbo_sld, kfl_funno_sld, kfl_funtn_sld, kfl_funtb_sld
  use def_solidz,         only : tbcod_sld, tncod_sld, tgcod_sld
  use def_solidz,         only : jacrot_du_dq_sld, jacrot_dq_du_sld
  use def_solidz,         only : bvess_sld, bvnat_sld
  use def_solidz,         only : veloc_sld, accel_sld, dunkn_sld, ddisp_sld
  use def_solidz,         only : ndofn_sld, dtcri_sld
  use def_solidz,         only : fint2_sld, fext2_sld, fine2_sld, fnatu_sld
  use def_solidz,         only : finte_sld, fexte_sld, macce_sld, frxid_sld
  use def_solidz,         only : fintt_sld, fextt_sld
  use def_solidz,         only : allie_sld, allwk_sld, allke_sld, etota_sld
  use def_solidz,         only : epsel_sld, lepse_sld
  use def_solidz,         only : green_sld, caust_sld, seqvm_sld
  use def_solidz,         only : caunp_sld, caulo_sld, cause_sld
  use def_solidz,         only : svegm_sld, srpro_sld
  use def_solidz,         only : isoch_sld, iswav_sld
  use def_solidz,         only : nvgij_inv_sld, nvgij_sld, nsvar_sld
  use def_solidz,         only : rorig_sld
  use def_solidz,         only : vmass_sld
  use def_solidz,         only : kfl_xfeme_sld, vxfem_sld, dxfem_sld, axfem_sld
  use def_solidz,         only : cranx_sld, crapx_sld, lnenr_sld, leenr_sld
  use def_solidz,         only : lecoh_sld, treff_sld, cohnx_sld, nocoh_sld
  use def_solidz,         only : dcmax_sld, dceff_sld, dslip_sld
  use def_solidz,         only : nopio_sld
  use def_solidz,         only : iwave_sld, kfl_sdvar_sld
  use def_solidz,         only : lawch_sld, lawst_sld, lawpl_sld
  use def_solidz,         only : relpa_sld, resid_sld
  use def_solidz,         only : kfl_restr_sld, restr_sld
  use def_solidz,         only : kfl_timei_sld, kfl_stead_sld, kfl_tisch_sld
  use def_solidz,         only : ncomp_sld, nprev_sld, nderi_sld, dtinv_sld
  use def_solidz,         only : kfl_moduf_sld, fiemo_sld
  use def_solidz,         only : axis1_sld, axis2_sld, axis3_sld, orien_sld, rmate_sld
  use def_solidz,         only : dttau_sld, ddism_sld
  use def_solidz,         only : cpu_ass_sol_sld
  use def_solidz,         only : SLD_STATIC_PROBLEM
  use def_solidz,         only : fcont_sld, saved_rhsid
  use def_solidz,         only : ebfil_sld, sbfil_sld
  use def_solidz,         only : kfl_immer_sld 
  use def_solidz,         only : water_sld, enden_sld
  use mod_sld_csys,       only : sld_csys_assign_material_axes
  use mod_sld_atm,        only : sld_atm_set_flags
  use mod_sld_cardiac_cycle, only : sld_cardiac_cycle_initialization 
  use mod_messages,       only : livinf
  use mod_communications, only : PAR_BROADCAST
  use mod_sld_strain,     only : strain_basis
#ifndef PROPER_ELEM_PRIVATE_OFF
  use mod_ker_sysnet, only: kfl_sysnet, sysnet_compute_pvloop, sysnet_compute_volumes, sysnet_initialize_guess
  use mod_ker_sysnet, only: sysnet_alya_couple, sysnet_update_phase, cavities_sysnet, in_startup_cycles
  use mod_ker_sysnet, only: in_rmax, in_rmin
#endif

  implicit none

  integer(ip), intent(in) :: itask       !< when this subroutine is called
  integer(ip)             :: ielem,pmate
  character(300)          :: messa
  real(rp)                :: prehg, premks
  
  select case(itask)
     !
     ! IMPORTANT: nvarp is the maximum number of wopos. it is defined in def_kintyp.
     ! check it before adding new woposes.
     !
  case(0_ip)
     !
     ! Modules initialization
     !
     call sld_atm_set_flags()
     call sld_cardiac_cycle_initialization()
     !
     ! Postprocess Variable names and types
     !
     call arrays_register((/'DISPL','VECTO','NPOIN','PRIMA'/),     displ, ENTITY_POSITION= 2_ip, TIME_POSITION= 3_ip )
     call arrays_register((/'VELOC','VECTO','NPOIN','PRIMA'/), veloc_sld, ENTITY_POSITION= 2_ip, TIME_POSITION= 3_ip )
     call arrays_register((/'ACCEL','VECTO','NPOIN','PRIMA'/), accel_sld, ENTITY_POSITION= 2_ip, TIME_POSITION= 3_ip )
     call arrays_register((/'SIGMA','TENSO','NPOIN','SECON'/) )
     call arrays_register((/'EPSIL','TENSO','NPOIN','SECON'/) )
     call arrays_register((/'LNEPS','TENSO','NPOIN','SECON'/) )
     call arrays_register((/'SEQVM','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'FIBER','VECTO','NPOIN','SECON'/) )
     call arrays_register((/'BVESS','VECTO','NPOIN','SECON'/) )
     call arrays_register((/'FIBEL','R3PVE','NPOIN','SECON'/) )
     call arrays_register((/'DXFEM','VECTO','NPOIN','SECON'/) )
     call arrays_register((/'VXFEM','VECTO','NPOIN','SECON'/) )
     call arrays_register((/'AXFEM','VECTO','NPOIN','SECON'/) )
     call arrays_register((/'LINEL','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'NSIGN','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'XYROT','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'ZXROT','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'DCOHE','SCALA','NELEM','SECON'/) )
     call arrays_register((/'LNENR','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'ERROR','VECTO','NPOIN','SECON'/) )
     call arrays_register((/'CRACK','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'GROUP','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'SIGEI','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'SDV1E','SCALA','NELEM','SECON'/) )
     call arrays_register((/'SDV2E','SCALA','NELEM','SECON'/) )
     call arrays_register((/'SDV3E','SCALA','NELEM','SECON'/) )
     call arrays_register((/'SDV4E','SCALA','NELEM','SECON'/) )
     call arrays_register((/'FIXRS','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'SRPRO','MULTI','NELEM','SECON'/) )
     call arrays_register((/'SDVAR','MULTI','NPOIN','SECON'/) )
     call arrays_register((/'NOSET','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'NELEM','SCALA','NELEM','SECON'/) )
     call arrays_register((/'PELCH','SCALA','NELEM','SECON'/) )
     call arrays_register((/'FIXNO','VECTO','NPOIN','SECON'/) )
     call arrays_register((/'EXAFO','VECTO','NPOIN','SECON'/) )
     call arrays_register((/'FORCF','VECTO','NPOIN','SECON'/) )
     call arrays_register((/'SEGMA','R3PVE','NELEM','SECON'/) )
     call arrays_register((/'FRXID','VECTO','NPOIN','SECON'/) )
     call arrays_register((/'SIGNP','R3PVE','NPOIN','SECON'/) )
     call arrays_register((/'SIGLO','R3PVE','NPOIN','SECON'/) )
     call arrays_register((/'AXIS1','VECTO','NELEM','SECON'/) )
     call arrays_register((/'AXIS2','VECTO','NELEM','SECON'/) )
     call arrays_register((/'AXIS3','VECTO','NELEM','SECON'/) )
     call arrays_register((/'ORIEN','SCALA','NELEM','SECON'/) )
     call arrays_register((/'SREAC','VECTO','NPOIN','SECON'/) )
     call arrays_register((/'PMATE','SCALA','NELEM','SECON'/) )
     call arrays_register((/'RMAT1','VECTO','NELEM','SECON'/) )
     call arrays_register((/'RMAT2','VECTO','NELEM','SECON'/) )
     call arrays_register((/'RMAT3','VECTO','NELEM','SECON'/) )
     call arrays_register((/'DIISO','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'VMISO','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'BOSET','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'EPSIS','TENSO','NPOIN','SECON'/) )
     call arrays_register((/'PK2  ','TENSO','NPOIN','SECON'/) )
     call arrays_register((/'INTLI','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'SVEGM','R3PVE','NELEM','PRIMA'/), svegm_sld, ENTITY_POSITION= 1_ip, TIME_POSITION= 3_ip )
     call arrays_register((/'CELEN','SCALA','NELEM','SECON'/) )
     call arrays_register((/'FIXBO','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'BOCOD','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'NOCOD','VECTO','NPOIN','SECON'/) )
     call arrays_register((/'ELNOR','VECTO','NELEM','SECON'/) )
     call arrays_register((/'FEXTE','VECTO','NPOIN','SECON'/) )
     call arrays_register((/'WETNO','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'SBVNA','VECTO','NPOIN','SECON'/) )
     call arrays_register((/'ELSET','SCALA','NELEM','SECON'/) )
     call arrays_register((/'PARTI','SCALA','NELEM','SECON'/) )
     call arrays_register((/'ELEPS','TENSO','NPOIN','SECON'/) )
     call arrays_register((/'ROTM1','VECTO','NPOIN','SECON'/) )
     call arrays_register((/'ROTM2','VECTO','NPOIN','SECON'/) )
     call arrays_register((/'ROTM3','VECTO','NPOIN','SECON'/) )
     call arrays_register((/'FCONT','VECTO','NPOIN','PRIMA'/), fcont_sld, ENTITY_POSITION= 2_ip)
     call arrays_register((/'DAMAG','MULTI','NELEM','SECON'/) )
     call arrays_register((/'FFLUI','VECTO','NPOIN','SECON'/) )
     call arrays_register((/'MICNL','SCALA','NELEM','SECON'/) )
     call arrays_register((/'MICCO','SCALA','NELEM','SECON'/) )
     call arrays_register((/'MICCV','SCALA','NELEM','SECON'/) )
     call arrays_register((/'EPRIN','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'CALCI','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'EBFIL','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'SBFIL','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'STACK','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'BVNAT','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'TEMPE','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'PRESS','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'NPOIN','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'EPSUL','SCALA','NPOIN','SECON'/) ) !strain in user supplied long direction
     call arrays_register((/'EPSUR','SCALA','NPOIN','SECON'/) ) !strain in user supplied radial direction
     call arrays_register((/'EPSUC','SCALA','NPOIN','SECON'/) ) !strain in user supplied circumferential direction
     call arrays_register((/'PELTY','SCALA','NELEM','SECON'/) )
     call arrays_register((/'ENDEN','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'WATER','SCALA','NPOIN','SECON'/) )
     call arrays_register((/'DONNA','SCALA','NPOIN','SECON'/) )
     !
     ! Sets variables
     !
     ! Boundary sets
     postp(1) % wobse (1)     = 'FORCE'  ! Force
     postp(1) % wobse (2)     = 'F_y'
     postp(1) % wobse (3)     = 'F_z'
     postp(1) % wobse (4)     = 'NORMA'  ! Normal Force
     postp(1) % wobse (5)     = 'DIBOX'  ! Averaged displacement in x-direction
     postp(1) % wobse (6)     = 'DIBOY'  ! Averaged displacement in y-direction
     postp(1) % wobse (7)     = 'DIBOZ'  ! Averaged displacement in z-direction
     postp(1) % wobse (8)     = 'DIBOU'  ! Averaged displacement (Magnitude)
     postp(1) % wobse (9)     = 'FRBOX'  ! Sum of reaction force in x-direction
     postp(1) % wobse (10)    = 'FRBOY'  ! Sum of reaction force in y-direction
     postp(1) % wobse (11)    = 'FRBOZ'  ! Sum of reaction force in z-direction
     postp(1) % wobse (12)    = 'FRBOU'  ! Sum of reaction force (Magnitude)
     postp(1) % wobse (13)    = 'FCONX'  ! Sum of contact force in x-direction
     postp(1) % wobse (14)    = 'FCONY'  ! Sum of contact force in y-direction
     postp(1) % wobse (15)    = 'FCONZ'  ! Sum of contact force in z-direction
     postp(1) % wobse (16)    = 'FCONT'  ! Sum of contact force (Magnitude)
     postp(1) % wobse (17)    = 'PRESS'  ! Pressure
     postp(1) % wobse (18)    = 'TRNOR'  ! Normal traction
     postp(1) % wobse (19)    = 'FCONO'  ! Sum of contact normal force 
     postp(1) % wobse (20)    = 'FCOT1'  ! Sum of contact tangent force 1 
     postp(1) % wobse (21)    = 'FCOT2'  ! Sum of contact tangent force 2 
     postp(1) % wobse (22)    = 'FOBOX'  ! Sum of forces in x-direction
     postp(1) % wobse (23)    = 'FOBOY'  ! Sum of forces in y-direction
     postp(1) % wobse (24)    = 'FOBOZ'  ! Sum of forces in z-direction
     postp(1) % wobse (25)    = 'FOBOU'  ! Sum of forces (Magnitude)
     ! Node sets
     postp(1) % wonse (1)     = 'DISPX'  ! x-displacement
     postp(1) % wonse (2)     = 'DISPY'  ! y-displacement
     postp(1) % wonse (3)     = 'DISPZ'  ! z-displacement
     postp(1) % wonse (4)     = 'VELOX'  ! x-velocity
     postp(1) % wonse (5)     = 'VELOY'  ! y-velocity
     postp(1) % wonse (6)     = 'VELOZ'  ! z-velocity
     postp(1) % wonse (7)     = 'ACCEX'  ! x-acceleration
     postp(1) % wonse (8)     = 'ACCEY'  ! y-acceleration
     postp(1) % wonse (9)     = 'ACCEZ'  ! z-acceleration
     postp(1) % wonse (10)    = 'FRXIX'  ! x-reaction force
     postp(1) % wonse (11)    = 'FRXIY'  ! y-reaction force
     postp(1) % wonse (12)    = 'FRXIZ'  ! z-reaction force
     postp(1) % wonse (13)    = 'SIGXX'  ! Stress-XX
     postp(1) % wonse (14)    = 'SIGYY'  ! Stress-YY
     postp(1) % wonse (15)    = 'SIGZZ'  ! Stress-ZZ
     postp(1) % wonse (16)    = 'SIGYZ'  ! Stress-YZ
     postp(1) % wonse (17)    = 'SIGXZ'  ! Stress-XZ
     postp(1) % wonse (18)    = 'SIGXY'  ! Stress-XY
     postp(1) % wonse (19)    = 'EPSXX'  ! Strain-XX
     postp(1) % wonse (20)    = 'EPSYY'  ! Strain-YY
     postp(1) % wonse (21)    = 'EPSZZ'  ! Strain-ZZ
     postp(1) % wonse (22)    = 'EPSYZ'  ! Strain-YZ
     postp(1) % wonse (23)    = 'EPSXZ'  ! Strain-XZ
     postp(1) % wonse (24)    = 'EPSXY'  ! Strain-XY
     postp(1) % wonse (25)    = 'LEPXX'  ! Logarithmic Strain-XX
     postp(1) % wonse (26)    = 'LEPYY'  ! Logarithmic Strain-YY
     postp(1) % wonse (27)    = 'LEPZZ'  ! Logarithmic Strain-ZZ
     postp(1) % wonse (28)    = 'LEPYZ'  ! Logarithmic Strain-YZ
     postp(1) % wonse (29)    = 'LEPXZ'  ! Logarithmic Strain-XZ
     postp(1) % wonse (30)    = 'LEPXY'  ! Logarithmic Strain-XY
     postp(1) % wonse (31)    = 'COORX'  ! x-coordinate
     postp(1) % wonse (32)    = 'COORY'  ! y-coordinate
     postp(1) % wonse (33)    = 'COORZ'  ! z-coordinate
     postp(1) % wonse (34)    = 'SEQVM'  ! Von Mises Stress
     postp(1) % wonse (35)    = 'FEXTX'  ! External force x-component
     postp(1) % wonse (36)    = 'FEXTY'  ! External force y-component
     postp(1) % wonse (37)    = 'FEXTZ'  ! External force z-component
     postp(1) % wonse (38)    = 'FCONX'  ! Contact force x-component
     postp(1) % wonse (39)    = 'FCONY'  ! Contact force y-component
     postp(1) % wonse (40)    = 'FCONZ'  ! Contact force z-component
     postp(1) % wonse (41)    = 'EPRIN'  ! Principal stretches
     postp(1) % wonse (42)    = 'BVESX'  ! Dirichlet x-component
     postp(1) % wonse (43)    = 'BVESY'  ! Dirichlet y-component
     postp(1) % wonse (44)    = 'BVESZ'  ! Dirichlet z-component
     postp(1) % wonse (45)    = 'SBVNX'  ! Solver Neumaan value x-component
     postp(1) % wonse (46)    = 'SBVNY'  ! Solver Neumaan value y-component
     postp(1) % wonse (47)    = 'SBVNZ'  ! Solver Neumaan value z-component
     postp(1) % wonse (48)    = 'TEMPE'  ! Temperature
     ! Element sets                            
     postp(1) % woese (1)     = 'INV1E'  ! 1st Invariant of green-lagrange strain tensor
     postp(1) % woese (2)     = 'INV2E'  ! 2nd Invariant of green-lagrange strain tensor
     postp(1) % woese (3)     = 'INV3E'  ! 3rd Invariant of green-lagrange strain tensor
     postp(1) % woese (4)     = 'GRL11'  ! component of green-lagrange strain tensor
     postp(1) % woese (5)     = 'GRL12'  ! component of green-lagrange strain tensor
     postp(1) % woese (6)     = 'GRL13'  ! component of green-lagrange strain tensor
     postp(1) % woese (7)     = 'GRL21'  ! component of green-lagrange strain tensor
     postp(1) % woese (8)     = 'GRL22'  ! component of green-lagrange strain tensor
     postp(1) % woese (9)     = 'GRL23'  ! component of green-lagrange strain tensor
     postp(1) % woese (10)    = 'GRL31'  ! component of green-lagrange strain tensor
     postp(1) % woese (11)    = 'GRL32'  ! component of green-lagrange strain tensor
     postp(1) % woese (12)    = 'GRL33'  ! component of green-lagrange strain tensor
     postp(1) % woese (13)    = 'EPSRT'  ! Averaged green strain in polar cylindrical (eps_rr + eps_tt)
     postp(1) % woese (14)    = 'LEPRT'  ! Averaged log strain in polar cylindrical (eps_rr + eps_tt)
     postp(1) % woese (15)    = 'VELOC'  ! Mean velocity
     postp(1) % woese (16)    = 'TMASS'  ! Total mass
     postp(1) % woese (17)    = 'ALLKE'  ! All kinetic energy
     postp(1) % woese (18)    = 'EBFIL'  ! Mean Strain projection to biofiber longitudinal
     postp(1) % woese (19)    = 'SBFIL'  ! Mean Stress projection to biofiber longitudinal
     postp(1) % woese(20)     = 'EPSUL'  ! Mean Strain along user supplied longitudinal direction
     postp(1) % woese(21)     = 'EPSUR'  ! Mean Strain along user supplied radial direction
     postp(1) % woese(22)     = 'EPSUC'  ! Mean Strain along user supplied circumferential direction
     postp(1) % woese(23)     = 'ENDEN'  ! Energy density from material model
     !
     ! Witness variables
     !
     postp(1) % wowit (1)     = 'DISPX'
     postp(1) % wowit (2)     = 'DISPY'
     postp(1) % wowit (3)     = 'DISPZ'
     postp(1) % wowit (4)     = 'VELOX'
     postp(1) % wowit (5)     = 'VELOY'
     postp(1) % wowit (6)     = 'VELOZ'
     postp(1) % wowit (7)     = 'ACCEX'
     postp(1) % wowit (8)     = 'ACCEY'
     postp(1) % wowit (9)     = 'ACCEZ'
     postp(1) % wowit (10)    = 'FRXIX'
     postp(1) % wowit (11)    = 'FRXIY'
     postp(1) % wowit (12)    = 'FRXIZ'
     postp(1) % wowit (13)    = 'SIGXX'
     postp(1) % wowit (14)    = 'SIGYY'
     postp(1) % wowit (15)    = 'SIGZZ'
     postp(1) % wowit (16)    = 'SIGYZ'
     postp(1) % wowit (17)    = 'SIGXZ'
     postp(1) % wowit (18)    = 'SIGXY'
     postp(1) % wowit (19)    = 'EPSXX'
     postp(1) % wowit (20)    = 'EPSYY'
     postp(1) % wowit (21)    = 'EPSZZ'
     postp(1) % wowit (22)    = 'EPSYZ'
     postp(1) % wowit (23)    = 'EPSXZ'
     postp(1) % wowit (24)    = 'EPSXY'
     postp(1) % wowit (25)    = 'LEPXX'
     postp(1) % wowit (26)    = 'LEPYY'
     postp(1) % wowit (27)    = 'LEPZZ'
     postp(1) % wowit (28)    = 'LEPYZ'
     postp(1) % wowit (29)    = 'LEPXZ'
     postp(1) % wowit (30)    = 'LEPXY'
     postp(1) % wowit (31)    = 'COORX'
     postp(1) % wowit (32)    = 'COORY'
     postp(1) % wowit (33)    = 'COORZ'
     postp(1) % wowit (34)    = 'SEQVM'
     postp(1) % wowit (35)    = 'SIGFI'
     postp(1) % wowit (36)    = 'MICRO'
     postp(1) % wowit (37)    = 'EBFIL'
     postp(1) % wowit (38)    = 'SBFIL'
     postp(1) % wowit (39)    = 'TEMPE'
     !
     ! Solvers
     !
     call soldef(-2_ip)                            ! Allocate memory for solvers
     solve(1) % kfl_solve = 1                      ! Output flag
     solve(1) % wprob     = 'DISPLACEMENT'
     solve(1) % ndofn     = ndime

     solve(1) % kfl_algso = SOL_SOLVER_RICHARDSON
     solve(1) % kfl_preco = SOL_LOCAL_MASS_MATRIX  ! explicit is the default

     solve(1) % kfl_force = 1                      ! Force solution continuity after calling the solver

     solve(2) % kfl_solve = 1                      ! Output flag ! OJO CAMBIAR
     solve(2) % wprob     = 'DISP-ROT'
     solve(2) % ndofn     = 6
     !solve(3) % kfl_solve = 1                     ! Output flag
     !
     ! Materials
     !
     nmate_sld =  nmate
     lmate_sld => lmate
     !
     ! Postprocess rotation origin
     !
     rorig_sld(1) = 0.0_rp
     rorig_sld(2) = 0.0_rp
     rorig_sld(3) = 0.0_rp
     !
     ! Force vector norms
     !
     fext2_sld = 0.0_rp
     fint2_sld = 0.0_rp
     fine2_sld = 0.0_rp
     fnatu_sld = 0.0_rp
     !
     ! Relaxation parameters
     !
     relpa_sld = 0.0_rp ! for solvers
     !
     ! Others
     !
     cpu_ass_sol_sld = 0.0_rp
     dtcri_sld       = 0.0_rp
     resid_sld       = 0.0_rp
     !
     !
     ! Nullify pointers
     !
     ! Fixities and bc codes
     nullify(kfl_fixbo_sld)
     nullify(kfl_fixno_sld)
     nullify(kfl_fixrs_sld)
     nullify(kfl_funbo_sld)
     nullify(kfl_funno_sld)
     nullify(kfl_funtn_sld)
     nullify(kfl_funtb_sld)
     nullify(kfl_immer_sld)
     nullify(tbcod_sld)
     nullify(tncod_sld)
     nullify(tgcod_sld)
     nullify(bvess_sld)
     nullify(bvnat_sld)
     nullify(jacrot_du_dq_sld)
     nullify(jacrot_dq_du_sld)
     ! Motion variables and mass
     nullify(dunkn_sld)
     nullify(ddisp_sld)
     nullify(vmass_sld)
     nullify(veloc_sld)
     nullify(accel_sld)
     ! Force vectors
     nullify(finte_sld)
     nullify(fexte_sld)
     nullify(macce_sld)
     nullify(fintt_sld)
     nullify(fextt_sld)
     nullify(frxid_sld)
     ! Energies
     nullify(allie_sld)
     nullify(allwk_sld)
     nullify(allke_sld)
     nullify(etota_sld)
     ! Stresses and strains
     nullify(green_sld)
     nullify(caust_sld)
     nullify(seqvm_sld)
     nullify(epsel_sld)
     nullify(lepse_sld)
     nullify(ebfil_sld)
     nullify(sbfil_sld)
     ! Contact
     nullify(fcont_sld)
     nullify(saved_rhsid)
     ! Other
     nullify(celen_sld)
     nullify(ddism_sld)
     nullify(dttau_sld)
     nullify(cause_sld)
     ! X-FEM
     nullify(dxfem_sld)
     nullify(vxfem_sld)
     nullify(axfem_sld)
     nullify(crapx_sld)
     nullify(cranx_sld)
     nullify(lnenr_sld)
     nullify(leenr_sld)
     nullify(lecoh_sld)

     nullify(dslip_sld)
     nullify(dceff_sld)
     nullify(dcmax_sld)
     nullify(nopio_sld)
     nullify(treff_sld)
     nullify(nocoh_sld)
     nullify(cohnx_sld)
     !
     ! Composite stuff
     !
     nullify(axis1_sld)
     nullify(axis2_sld)
     nullify(axis3_sld)
     nullify(rmate_sld)
     nullify(orien_sld)
     nullify(svegm_sld)
     nullify(srpro_sld)
     ! 
     nullify(caunp_sld)
     nullify(caulo_sld)
     !
     nullify(isoch_sld)
     nullify(iswav_sld)
     !
     ! Timings
     !
     call moduls_allocate_timing(10_ip)
     call times(1) % register('time step'         ,'step beginning' )
     call times(2) % register('energy'            ,'step beginning' )
     
     call times(3) % register('gather'            ,'element assembly' )
     call times(4) % register('computation'       ,'element assembly' )
     call times(5) % register('electro-mechanical','element assembly' )
     call times(6) % register('scatter'           ,'element assembly' )
     call times(7) % register('element matrix'    ,'element assembly' )
     !
     ! This initizlization is required to activate timing for the master
     !
     call times(3) % ini() ; call times(3) % add() 
     call times(4) % ini() ; call times(4) % add()
     call times(5) % ini() ; call times(5) % add()
     call times(6) % ini() ; call times(6) % add()
     call times(7) % ini() ; call times(7) % add()
     !
     ! Poroelasticity
     !
     nullify(enden_sld)
     nullify(water_sld)
     nullify(donna_gp)

  case(1_ip)
     !
     ! Voigt dimensions
     !
     nvoig_sld = (ndime + 1) * ndime / 2
     !
     ! Voigt rule conversion
     !
     if ( ndime == 2_ip) then
        !
        ! 2D Voigt rule
        ! \sigma_{ij} \sigma_{a}
        !         ij          a
        !         11          1
        !         22          2
        !         12          3
        !
        nvgij_inv_sld(1,1) = 1_ip
        nvgij_inv_sld(2,2) = 2_ip
        nvgij_inv_sld(1,2) = 3_ip
        ! Symmetric
        nvgij_inv_sld(2,1) = nvgij_inv_sld(1,2)

        nvgij_sld(1,1) = 1
        nvgij_sld(1,2) = 1
        nvgij_sld(2,1) = 2
        nvgij_sld(2,2) = 2
        nvgij_sld(3,1) = 1
        nvgij_sld(3,2) = 2

     else if ( ndime == 3_ip) then
        !
        ! 3D Voigt rule
        ! \sigma_{ij} \sigma_{a}
        !         ij          a
        !         11          1
        !         22          2
        !         33          3
        !         23          4
        !         13          5
        !         12          6
        !
        nvgij_inv_sld(1,1) = 1_ip
        nvgij_inv_sld(2,2) = 2_ip
        nvgij_inv_sld(3,3) = 3_ip
        nvgij_inv_sld(2,3) = 4_ip
        nvgij_inv_sld(1,3) = 5_ip
        nvgij_inv_sld(1,2) = 6_ip
        ! Symmetric
        nvgij_inv_sld(2,1) = nvgij_inv_sld(1,2)
        nvgij_inv_sld(3,1) = nvgij_inv_sld(1,3)
        nvgij_inv_sld(3,2) = nvgij_inv_sld(2,3)

        nvgij_sld(1,1) = 1
        nvgij_sld(1,2) = 1
        nvgij_sld(2,1) = 2
        nvgij_sld(2,2) = 2
        nvgij_sld(3,1) = 3
        nvgij_sld(3,2) = 3
        nvgij_sld(4,1) = 2
        nvgij_sld(4,2) = 3
        nvgij_sld(5,1) = 1
        nvgij_sld(5,2) = 3
        nvgij_sld(6,1) = 1
        nvgij_sld(6,2) = 2

     end if
     !
     ! Dimensions
     !
     ndofn_sld = ndime
     if (kfl_xfeme_sld == 1) ndofn_sld = 2*ndime
     !
     ! Number of temporal variables
     !
     if(      kfl_tisch_sld == 1 ) then
        !
        ! Newmark-Beta Scheme/Central differences
        !
        ncomp_sld = 4_ip
     
     else if( kfl_tisch_sld == 2 ) then
        !
        ! Dissipative Tchamwa - Wielgosz scheme
        !
        ncomp_sld = 3_ip

     else if( kfl_tisch_sld == 3 ) then
        !
        ! Runge-Kutta Scheme
        !
        ncomp_sld = 3_ip
        nderi_sld = ndime*2

     else
        ncomp_sld = 3_ip
        nderi_sld = ndime*2

     end if
     nprev_sld = min(3_ip,ncomp_sld)  ! Last time step or global iteration
     !
     ! Time variables
     !
     dtinv_sld = 0.0_rp
     kfl_timei = 1_ip
     if ( kfl_timei_sld == SLD_STATIC_PROBLEM ) dtinv_sld = 1.0_rp
     kfl_stead_sld = 0_ip
     !
     ! Dimensions state dependent variables
     !
     if( kfl_sdvar_sld == 1_ip ) then
        if( INOTMASTER ) then
           nsvar_sld = 0
           !------------------------------------------------------------
           !$OMP PARALLEL  DO                                          &
           !$OMP SCHEDULE  ( DYNAMIC , par_omp_nelem_chunk )           &
           !$OMP DEFAULT   (NONE)                                      &
           !$OMP PRIVATE   (ielem,pmate)                               &
           !$OMP SHARED    (nelem,lmate,lawst_sld,lawch_sld,           &
           !$OMP            lawpl_sld,                                 & 
           !$OMP            par_omp_nelem_chunk)                       &
           !$OMP REDUCTION (MAX:nsvar_sld)
           !------------------------------------------------------------
           do ielem = 1,nelem
              pmate = lmate(ielem)
              if(      lawch_sld(pmate) == 800_ip .or. &
                       lawch_sld(pmate) == 904_ip .or. & 
                       lawch_sld(pmate) == 905_ip ) then
                 ! Cohesive elements: 904 and 905 and contact law 800
                 nsvar_sld = max(nsvar_sld,2_ip)
              else if( lawst_sld(pmate) == 152_ip ) then
                 ! sm152
                 nsvar_sld = max(nsvar_sld,7_ip)
              else if( lawst_sld(pmate) == 153_ip ) then
                 ! sm153
                 nsvar_sld = max(nsvar_sld,3_ip)
              else if( lawst_sld(pmate) == 160_ip ) then
                 ! nitinol
                 nsvar_sld = max(nsvar_sld,3_ip)
              else if( lawst_sld(pmate) == 154_ip ) then
                 ! sm154
                 nsvar_sld = max(nsvar_sld,20_ip)
              else if( lawst_sld(pmate) == 0_ip .and. &
                       lawpl_sld(pmate) == 1_ip ) then
                 ! smXXX
                 nsvar_sld = max(nsvar_sld,20_ip)
              else if( lawst_sld(pmate) == 200_ip ) then
                 ! sm200
                 nsvar_sld = max(nsvar_sld,22_ip)
              end if
           end do
           !$OMP END PARALLEL DO
           !------------------------------------------------------------
        end if
        !
        ! Parall: look for the maximum over all subdomains
        !
        call PAR_MAX(nsvar_sld,'IN MY CODE')
     end if
     !
     ! Stresses field
     !
     if( kfl_restr_sld < 0 ) then
        restr_sld => xfiel(-kfl_restr_sld) % a(:,:,1)
     end if
     !
     ! Rigidity modulator field
     !
     if( kfl_moduf_sld(1) > 0 ) then
        fiemo_sld => xfiel(kfl_moduf_sld(1)) % a(:,:,1)
     end if
     !
     ! Polarization wave index
     !
     iwave_sld = 0
     !
     ! Strain basis
     !
     if (strain_basis%on) then
        call strain_basis%assign_fields()
     end if

  case(2_ip)
     !
     ! Sysnet
     !
#ifndef PROPER_ELEM_PRIVATE_OFF
     if (kfl_sysnet) then
        messa = 'SOLIDZ: SYSNET VASCULAR SYSTEM NETWORK ENABLED'
        call sysnet_initialize_guess()
        call livinf(0_ip,messa,1_ip)
        messa = '        SYSNET: STARTUP PHASE...'
        call livinf(0_ip,messa,1_ip)
        messa = &
             '           VALVE RESISTANCES  (OPEN - CLOSE in mmHg s / mL): '
        call livinf(0_ip,messa,1_ip)
        messa = &
             '               MITRAL:    '//trim(retost(in_rmin(1)))//' - '//trim(retost(in_rmax(1)))
        call livinf(0_ip,messa,1_ip)
        messa = &
             '               AORTIC:    '//trim(retost(in_rmin(2)))//' - '//trim(retost(in_rmax(2)))
        call livinf(0_ip,messa,1_ip)
        messa = &
             '               TRICUSPID: '//trim(retost(in_rmin(3)))//' - '//trim(retost(in_rmax(3)))
        call livinf(0_ip,messa,1_ip)
        messa = &
             '               PULMONIC:  '//trim(retost(in_rmin(4)))//' - '//trim(retost(in_rmax(4)))
        call livinf(0_ip,messa,1_ip)
        ! compute cavity initial volumes to define a synthetic volume function ONLY for the startup phase
        call sysnet_compute_volumes(ITASK_INITIA)   
        ! update phase
        call sysnet_update_phase(ITASK_INITIA)
        if (INOTSLAVE) then
           !
           ! STARTUP LOOP: perform initialization loop, only the master
           !
           call sysnet_compute_pvloop(ITASK_INITIA)
           !
           ! Once the startup is finished, update alya pressures coming from sysnet, converting from mm/hg to CGS
           !
           call sysnet_alya_couple(2_ip)
        end if
        ! Broadcast pressures to all the slaves
        if (INOTMASTER) then
           cavities_sysnet(1) % pres(TIME_N) = 0.0_rp
           cavities_sysnet(2) % pres(TIME_N) = 0.0_rp
        end if
        call PAR_BROADCAST(cavities_sysnet(1) % pres(TIME_N))
        call PAR_BROADCAST(cavities_sysnet(2) % pres(TIME_N))
        messa = &
             '           STARTUP CYCLES COMPLETED: '//trim(intost(in_startup_cycles))
        call livinf(0_ip,messa,1_ip)
        messa = &
             '           RESULTS WRITTEN IN *-startup-* FILES '
        call livinf(0_ip,messa,1_ip)
        premks= cavities_sysnet(1) % pres(TIME_N) 
        prehg = cavities_sysnet(1) % pres(TIME_N) /  1333.22_rp
        messa = &
             '           POST-STARTUP (PRESTRESS) LEFT VENTRICLE PRESSURE  : '//trim(retost(premks))//' ('//trim(retost(prehg))//' mmhg)'
        call livinf(0_ip,messa,1_ip)
        premks= cavities_sysnet(2) % pres(TIME_N) 
        prehg = cavities_sysnet(2) % pres(TIME_N) /  1333.22_rp
        messa = &
             '           POST-STARTUP (PRESTRESS) RIGHT VENTRICLE PRESSURE : '//trim(retost(premks))//' ('//trim(retost(prehg))//' mmhg)'
        call livinf(0_ip,messa,1_ip)
        messa = &
             '        SYSNET: STARTUP PHASE... END.'
        call livinf(0_ip,messa,1_ip)
     end if
#endif

  end select

end subroutine sld_inivar

