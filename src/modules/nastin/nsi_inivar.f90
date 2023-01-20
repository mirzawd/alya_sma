!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_inivar.f90
!> @author  Guillaume Houzeaux
!> @brief   Initialize some variables
!> @details Initialize some variables\n
!!          ITASK=1 ... When starting the run (from Turnon)\n
!!          ITASK=2 ... First time step. This is needed as some variables
!!                      are not initialized before\n
!!          ITASK=3 ... When starting a time step (from nsi_begste)\n
!> @}
!------------------------------------------------------------------------
subroutine nsi_inivar(itask)
  use def_parame
  use def_master
  use def_nastin
  use def_domain
  use def_solver
  use mod_memchk
  use def_kintyp_performance
  use def_kermod,               only : gasco, kfl_adj_prob
  use mod_nsi_schur_complement, only : nsi_schur_complement_initialization
  use mod_nsi_multi_step_fs,    only : nsi_multi_step_fs_initialization
  use mod_nsi_fractional_step,  only : nsi_fractional_step_initialization
  use mod_arrays,               only : arrays_register
  use mod_arrays,               only : arrays_number
  use mod_moduls_conf,          only : moduls_allocate_timing
  use mod_output_postprocess,   only : output_postprocess_check_variable_postprocess
  
#ifndef NASTIN_PRIVATE_OFF
  use mod_nsi_actuator_disc,    only : nsi_cfdwak
#endif

  implicit none
  integer(ip), intent(in) :: itask

  select case(itask)

  case(0_ip)
     !
     ! Modules initialization
     !
     call nsi_schur_complement_initialization()
     call nsi_multi_step_fs_initialization()
     call nsi_fractional_step_initialization()
     !
     ! Postprocess Variable names and types
     !
     call arrays_register((/'VELOC','VECTO','NPOIN','PRIMA'/),veloc,ENTITY_POSITION=2_ip,TIME_POSITION=3_ip,EXCLUDE_RST=(/2_ip,3_ip/))
     call arrays_register((/'PRESS','SCALA','NPOIN','PRIMA'/),press,ENTITY_POSITION=1_ip,TIME_POSITION=2_ip,EXCLUDE_RST=(/2_ip,3_ip/))
     call arrays_register((/'STREA','SCALA','NPOIN','SECON'/))      
     call arrays_register((/'RESID','SCALA','NPOIN','SECON'/))      
     call arrays_register((/'VESGS','R3PVE','NELEM','PRIMA'/),vesgs,ENTITY_POSITION=1_ip,TIME_POSITION=3_ip,EXCLUDE_RST=(/2_ip/))
     call arrays_register((/'BOUND','VECTO','NPOIN','SECON'/))      
     call arrays_register((/'DENSI','SCALA','NPOIN','PRIMA'/),densi,ENTITY_POSITION=1_ip,TIME_POSITION=2_ip)
     call arrays_register((/'PMV  ','SCALA','NPOIN','SECON'/))
     call arrays_register((/'PPD  ','SCALA','NPOIN','SECON'/))
     call arrays_register((/'MACHN','SCALA','NPOIN','SECON'/))
     call arrays_register((/'TANGE','VECTO','NPOIN','SECON'/))
     call arrays_register((/'VORTI','VECTO','NPOIN','SECON'/))
     call arrays_register((/'MODVO','SCALA','NPOIN','SECON'/))
     call arrays_register((/'LAMB2','SCALA','NPOIN','SECON'/))
     call arrays_register((/'GROUP','SCALA','NPOIN','SECON'/))
     call arrays_register((/'LINEL','SCALA','NPOIN','SECON'/))
     call arrays_register((/'VISCO','SCALA','NPOIN','SECON'/))
     call arrays_register((/'FIXPR','SCALA','NPOIN','SECON'/))
     call arrays_register((/'FIXNO','VECTO','NPOIN','SECON'/))
     call arrays_register((/'TAU  ','SCALA','NPOIN','SECON'/))                                        ! FREE POSITION
     call arrays_register((/'AVVEL','VECTO','NPOIN','PRIMA'/),avvel_nsi,ENTITY_POSITION=2_ip)
     call arrays_register((/'SCHUR','SCALA','NPOIN','SECON'/))         
     call arrays_register((/'VEPRO','VECTO','NPOIN','PRIMA'/),vepro_nsi,ENTITY_POSITION=2_ip)
     call arrays_register((/'PRPRO','SCALA','NPOIN','PRIMA'/),prpro_nsi,ENTITY_POSITION=1_ip)
     call arrays_register((/'YPLUS','SCALA','NPOIN','SECON'/))         
     call arrays_register((/'LIMIT','SCALA','NPOIN','SECON'/))         
     call arrays_register((/'GRPRO','VECTO','NPOIN','PRIMA'/),grpro_nsi,ENTITY_POSITION=2_ip)
     call arrays_register((/'AVPRE','SCALA','NPOIN','PRIMA'/),avpre_nsi,ENTITY_POSITION=1_ip)
     call arrays_register((/'VEERR','SCALA','NPOIN','SECON'/))
     call arrays_register((/'PECLE','SCALA','NPOIN','SECON'/))
     call arrays_register((/'HYDRO','SCALA','NPOIN','SECON'/))
     call arrays_register((/'DUDT ','VECTO','NPOIN','SECON'/))
     call arrays_register((/'D2UDT','VECTO','NPOIN','SECON'/))
     call arrays_register((/'DT1  ','SCALA','NPOIN','SECON'/))
     call arrays_register((/'DT2  ','SCALA','NPOIN','SECON'/))
     call arrays_register((/'PRERR','SCALA','NPOIN','SECON'/))
     call arrays_register((/'LINVE','SCALA','NPOIN','SECON'/))
     call arrays_register((/'DEFOR','VECTO','NPOIN','SECON'/))
     call arrays_register((/'NODPR','SCALA','NPOIN','SECON'/))
     call arrays_register((/'-----','VECTO','NPOIN','SECON'/))                                       ! total traction                              
     call arrays_register((/'VELO2','VECTO','NPOIN','SECON'/))                                       ! For BDF restart   
     call arrays_register((/'PRES2','SCALA','NPOIN','SECON'/))                                       ! For BDF restart   
     call arrays_register((/'VELO3','VECTO','NPOIN','SECON'/))                                       ! For BDF restart   
     call arrays_register((/'PRES3','SCALA','NPOIN','SECON'/))                                       ! For BDF restart   
     call arrays_register((/'MODVE','SCALA','NPOIN','SECON'/))                                                                                           
     call arrays_register((/'AVTAN','VECTO','NPOIN','PRIMA'/),avtan_nsi,ENTITY_POSITION=2_ip)        ! Average tangent component of boundary traction
     call arrays_register((/'VELOM','VECTO','NPOIN','SECON'/))                                                                                           
     call arrays_register((/'FORCF','VECTO','NPOIN','PRIMA'/),forcf,ENTITY_POSITION=2_ip)                                                                      
     call arrays_register((/'INTFO','VECTO','NPOIN','SECON'/))                                                                                           
     call arrays_register((/'USTAR','SCALA','NPOIN','SECON'/))                                                                                           
     call arrays_register((/'TRACT','VECTO','NPOIN','SECON'/))                                       ! traction                                                        
     call arrays_register((/'FIXRS','SCALA','NPOIN','SECON'/))                                                                                 
     call arrays_register((/'FVELO','VECTO','NPOIN','SECON'/))                                       ! new velocity for filter output             
     call arrays_register((/'FPRES','SCALA','NPOIN','SECON'/))                                       ! new pressure for filter output             
     call arrays_register((/'FTANG','VECTO','NPOIN','SECON'/))                                       ! new wall shear stress for filter output    
     call arrays_register((/'AVVE2','VECTO','NPOIN','PRIMA'/),avve2_nsi,ENTITY_POSITION=2_ip)        ! Average Vi**2                              
     call arrays_register((/'AVVXY','VECTO','NPOIN','PRIMA'/),avvxy_nsi,ENTITY_POSITION=2_ip)        ! Average Vx*Vy, Vy*Vz, Vz*Vx                
     call arrays_register((/'AVPR2','SCALA','NPOIN','PRIMA'/),avpr2_nsi,ENTITY_POSITION=1_ip)        ! Average Pr**2                              
     call arrays_register((/'MESHR','VECTO','NPOIN','SECON'/))                                       ! Mesh Rotation                              
     call arrays_register((/'NODEF','SCALA','NPOIN','SECON'/))                                       ! If reaction force are computed on nodes    
     call arrays_register((/'MOMEN','VECTO','NPOIN','SECON'/))                                       ! Momentum equaiton residual                 
     call arrays_register((/'FIXPP','SCALA','NPOIN','SECON'/))                                       ! Pressure fixity                            
     call arrays_register((/'BVNAT','SCALA','NPOIN','SECON'/))                                       ! Algebraic Neumann condition                
     call arrays_register((/'-----','VECTO','NPOIN','SECON'/))                                       ! Algebraic reaction (b-Ax)                  
     call arrays_register((/'-----','VECTO','NPOIN','SECON'/))                                                                             
     call arrays_register((/'-----','VECTO','NPOIN','SECON'/))                                                                             
     call arrays_register((/'-----','VECTO','NPOIN','SECON'/))                                                                             
     call arrays_register((/'-----','VECTO','NPOIN','SECON'/))                                                                             
     call arrays_register((/'-----','VECTO','NPOIN','SECON'/))                                                                             
     call arrays_register((/'WETNO','SCALA','NPOIN','SECON'/))                                                                             
     call arrays_register((/'AVVRE','VECTO','WHATE','SECON'/))                                       ! For time-averaging restart              
     call arrays_register((/'TURMU','SCALA','NPOIN','SECON'/))                                       ! LES turbulent viscosity                 
     call arrays_register((/'AVMUT','SCALA','NPOIN','PRIMA'/),avmut_nsi,ENTITY_POSITION=1_ip)        ! Average turbulent viscosity             
     call arrays_register((/'AVSTX','VECTO','NPOIN','PRIMA'/),avstx_nsi,ENTITY_POSITION=2_ip)        ! Average stress (mu_t * grad(u))         
     call arrays_register((/'AVSTY','VECTO','NPOIN','PRIMA'/),avsty_nsi,ENTITY_POSITION=2_ip)        ! Average stress (mu_t * grad(v))         
     call arrays_register((/'AVSTZ','VECTO','NPOIN','PRIMA'/),avstz_nsi,ENTITY_POSITION=2_ip)        ! Average stress (mu_t * grad(w))         
     call arrays_register((/'BUBBL','SCALA','NELEM','SECON'/))                                       ! Pressure bubble                         
     call arrays_register((/'MFCFO','SCALA','NPOIN','SECON'/))                                       ! For mass flow control restart           
     call arrays_register((/'UBPRE','SCALA','NPOIN','SECON'/))                                       ! For mass flow control restart           
     call arrays_register((/'NORMA','SCALA','NPOIN','SECON'/))                                       ! Normals                                 
     call arrays_register((/'SENSM','VECTO','NPOIN','SECON'/))                                       ! Mesh sensitivities                      
     call arrays_register((/'LAGRA','VECTO','NPOIN','PRIMA'/),lagra_nsi,ENTITY_POSITION=2_ip)        ! Lagrange multiplier                     
     call arrays_register((/'ENVEL','VECTO','NPOIN','PRIMA'/),envel_nsi,ENTITY_POSITION=2_ip)                                                                   
     call arrays_register((/'ENPRE','SCALA','NPOIN','PRIMA'/),enpre_nsi,ENTITY_POSITION=1_ip)                                                                   
     call arrays_register((/'ENVE2','VECTO','NPOIN','PRIMA'/),enve2_nsi,ENTITY_POSITION=2_ip)        ! Ensemble Vi**2                                               
     call arrays_register((/'ENVXY','VECTO','NPOIN','PRIMA'/),envxy_nsi,ENTITY_POSITION=2_ip)        ! Ensemble Vx*Vy, Vy*Vz, Vz*Vx                                 
     call arrays_register((/'ENPR2','SCALA','NPOIN','PRIMA'/),enpr2_nsi,ENTITY_POSITION=1_ip)        ! Ensemble Pr**2                                               
     call arrays_register((/'ENTAN','VECTO','NPOIN','PRIMA'/),entan_nsi,ENTITY_POSITION=2_ip)        ! Ensemble tau                                                 
     call arrays_register((/'ENMUT','SCALA','NPOIN','PRIMA'/),enmut_nsi,ENTITY_POSITION=1_ip)        ! Ensemble turbulent viscosity                                 
     call arrays_register((/'ENSTX','VECTO','NPOIN','PRIMA'/),enstx_nsi,ENTITY_POSITION=2_ip)        ! Ensemble stress (mu_t * grad(u))                             
     call arrays_register((/'ENSTY','VECTO','NPOIN','PRIMA'/),ensty_nsi,ENTITY_POSITION=2_ip)        ! Ensemble stress (mu_t * grad(v))                             
     call arrays_register((/'ENSTZ','VECTO','NPOIN','PRIMA'/),enstz_nsi,ENTITY_POSITION=2_ip)        ! Ensemble stress (mu_t * grad(w))                             
     call arrays_register((/'BTRAC','VECTO','NPOIN','PRIMA'/),btrac_nsi,ENTITY_POSITION=2_ip)        ! For two-layer model restart - LES part                       
     call arrays_register((/'TRACR','VECTO','NPOIN','PRIMA'/),tluav_nsi,ENTITY_POSITION=2_ip)        ! For two-layer model restart - RANS part                      
     call arrays_register((/'TLUAV','VECTO','NPOIN','PRIMA'/),tracr_nsi,ENTITY_POSITION=2_ip)        ! For two-layer model restart - Average velocity               
     call arrays_register((/'AVVTA','VECTO','WHATE','SECON'/))                                       ! Average viscous part of tangential stres
     call arrays_register((/'VAFOR','VECTO','NPOIN','PRIMA'/),vafor_nsi,ENTITY_POSITION=2_ip)        ! VARIATIONAL FORCE                                             
     call arrays_register((/'AVVAF','VECTO','NPOIN','PRIMA'/),avvaf_nsi,ENTITY_POSITION=2_ip)        ! Average VARIATIONAL FORCE                                     
     call arrays_register((/'NOTRA','VECTO','NPOIN','PRIMA'/),notra_nsi,ENTITY_POSITION=2_ip)        ! VARIATIONAL Tangential traction                               
     call arrays_register((/'AVNTR','VECTO','NPOIN','PRIMA'/),avntr_nsi,ENTITY_POSITION=2_ip)        ! Average VARIATIONAL Tangential traction  - beware running average             
     call arrays_register((/'AVGTR','VECTO','NPOIN','PRIMA'/),avgtr_nsi,ENTITY_POSITION=2_ip)        ! Average gradient based Tangential traction  - beware running average
!     call arrays_register,(/'FANSW','SCALA','NPOIN','SECON'/))                                       ! Factor for no slip wall law - 4 restart    ! NO longer needed - hope this is ok
     call arrays_register((/'PORFO','VECTO','NPOIN','SECON'/))                                       ! Porous force  
     call arrays_register((/'MASRH','SCALA','NPOIN','PRIMA'/),mass_rho_nsi,ENTITY_POSITION=1_ip,TIME_POSITION=2_ip) ! M * rho for low-Mach                                     
     call arrays_register((/'MOMSK','VECTO','NPOIN','SECON'/))                                       ! Momentum source from partis               
     call arrays_register((/'DRHOD','SCALA','NPOIN','SECON'/))                                       ! unlumped drho/dt for restarting low mach  
     call arrays_register((/'VEOLD','VECTO','NPOIN','PRIMA'/),veold_nsi,     ENTITY_POSITION=2_ip)   ! Old velocity        
     call arrays_register((/'TAUIB','VECTO','NPOIN','PRIMA'/),tauib_nsi,     ENTITY_POSITION=3_ip)   ! tauib_nsi          
     call arrays_register((/'DUNKN','VECTO','NPOIN','PRIMA'/),dunkn_nsi,     ENTITY_POSITION=2_ip)   ! dunkn_nsi    
     call arrays_register((/'DUNKP','SCALA','NPOIN','PRIMA'/),dunkp_nsi,     ENTITY_POSITION=1_ip)   ! dunkp_nsi                                                               
     call arrays_register((/'BUBBL','SCALA','NELEM','PRIMA'/),bubble_nsi    ,ENTITY_POSITION=1_ip)   ! Bubble
     call arrays_register((/'BUAQQ','SCALA','NELEM','PRIMA'/),bubble_aqq_nsi,ENTITY_POSITION=1_ip)   ! Bubble 
     call arrays_register((/'BUAQU','MATRI','NELEM','PRIMA'/),bubble_aqu_nsi,ENTITY_POSITION=2_ip)   ! Bubble 
     call arrays_register((/'BUAQP','MATRI','NELEM','PRIMA'/),bubble_aqp_nsi,ENTITY_POSITION=2_ip)   ! Bubble 
     call arrays_register((/'BUBQ ','SCALA','NELEM','PRIMA'/),bubble_bq_nsi ,ENTITY_POSITION=1_ip)   ! Bubble
     call arrays_register((/'VEPAV','VECTO','NPOIN','SECON'/))                                       ! Velocity at boundaries projected on nodes
     call arrays_register((/'QCRIT','SCALA','NPOIN','SECON'/))                                       ! Q criterion
     call arrays_register((/'AVMOS','VECTO','NPOIN','PRIMA'/),avmos_nsi,ENTITY_POSITION=2_ip)        ! Average momentum source from spray
     call arrays_register((/'EIGEN','SCALA','NPOIN','SECON'/))                                       ! Projection of eigenvalue time step   
     call arrays_register((/'AVMFL','VECTO','NPOIN','PRIMA'/),av_mass_flux_nsi,ENTITY_POSITION=2_ip) ! Average mass flux 
     call arrays_register((/'AVRUU','VECTO','NPOIN','PRIMA'/),av_mom_flux_diag_nsi,ENTITY_POSITION=2_ip) ! Average rho*Vx*Vx, rho*Vy*Vy, rho*Vz*Vz
     call arrays_register((/'AVRUV','VECTO','NPOIN','PRIMA'/),av_mom_flux_off_nsi,ENTITY_POSITION=2_ip)  ! Average rho*Vx*Vy, rho*Vy*Vz, rho*Vz*Vx
     call arrays_register((/'DISTO','SCALA','NPOIN','SECON'/))
     call arrays_register((/'VMAXP','VECTO','NPOIN','PRIMA'/), vmaxp_nsi, ENTITY_POSITION=2_ip)          ! Maximum node-wise modulus ov velocity
     call arrays_register((/'VMINP','VECTO','NPOIN','SECON'/), vminp_nsi, ENTITY_POSITION=2_ip)          ! Minimum node-wise modulus ov velocity
     call arrays_register((/'VAVEP','VECTO','NPOIN','SECON'/), vavep_nsi, ENTITY_POSITION=2_ip)          ! Average node-wise modulus ov velocity
     call arrays_register((/'PINDE','SCALA','NPOIN','SECON'/), pinde_nsi, ENTITY_POSITION=2_ip)          ! Pulsatility index PI=(v_max-v_min)/v_ave
     call arrays_register((/'SBVNA','VECTO','NPOIN','SECON'/))                                       ! Algebraic Neumann boundary conditions 
     call arrays_register((/'FSIFO','VECTO','NPOIN','SECON'/), fsifo_nsi)                            ! FSI body force received from solid in IBcoupling for deformable bodies
     call arrays_register((/'DISTA','SCALA','NPOIN','SECON'/))                                       ! Distance
     call arrays_register((/'VELDI','VECTO','NPOIN','SECON'/))                                       ! Velocity-mesh velocity
     !                                                         
     ! Set variables
     !
     postp(1) % woese (1)     = 'VELOC'
     postp(1) % woese (2)     = 'VORTI'
     postp(1) % woese (3)     = 'KINET'
     postp(1) % woese (4)     = 'DIVU2'
     postp(1) % woese (5:7)   = 'MOMEN'
     postp(1) % woese (8)     = 'MEANP'
     postp(1) % woese (9)     = 'DISTO'
     postp(1) % woese (10)    = 'VMAXP'
     postp(1) % woese (11)    = 'VMINP'
     postp(1) % woese (12)    = 'VAVEP'
     postp(1) % woese (13)    = 'PINDE'

     postp(1) % wobse (1)     = 'MEANP'   ! Mean pressure
     postp(1) % wobse (2)     = 'MASS '   ! Mass rho*u.n
     postp(1) % wobse (3)     = 'FORCE'   ! Viscous force
     postp(1) % wobse (4)     = 'F_v_y'
     postp(1) % wobse (5)     = 'F_v_z'
     postp(1) % wobse (6)     = 'F_p_x'   ! Pressure force
     postp(1) % wobse (7)     = 'F_p_y'
     postp(1) % wobse (8)     = 'F_p_z'
     postp(1) % wobse (9)     = 'TORQU'   ! Viscous torque
     postp(1) % wobse (10)    = 'T_v_y'
     postp(1) % wobse (11)    = 'T_v_z'
     postp(1) % wobse (12)    = 'T_p_x'   ! Pressure torque
     postp(1) % wobse (13)    = 'T_p_y'
     postp(1) % wobse (14)    = 'T_p_z'
     postp(1) % wobse (15)    = 'MEANY'   ! Mean y+
     postp(1) % wobse (16)    = 'MEANV'   ! Mean velocity
     postp(1) % wobse (17)    = 'WETFO'   ! Viscous wet force
     postp(1) % wobse (18)    = 'Fwv_y'
     postp(1) % wobse (19)    = 'Fwv_z'
     postp(1) % wobse (20)    = 'Fwp_x'   ! Pressure wet force
     postp(1) % wobse (21)    = 'Fwp_y'
     postp(1) % wobse (22)    = 'Fwp_z'
     postp(1) % wobse (23)    = 'WETSU'
     postp(1) % wobse (24)    = 'INTFX'   ! Internal force
     postp(1) % wobse (25)    = 'INTFY'
     postp(1) % wobse (26)    = 'INTFZ'
     postp(1) % wobse (27)    = 'REATT'   ! Reattachment
     postp(1) % wobse (28)    = 'REATY'
     postp(1) % wobse (29)    = 'REATZ'
     postp(1) % wobse (30)    = 'REATM'
     postp(1) % wobse (31)    = 'UDOTN'   ! u.n
     postp(1) % wobse (32)    = 'FPORX'   ! Porous force -- ojo no meti como si fuera boundary pero en realidad es volumetrico
     postp(1) % wobse (33)    = 'FPORY'
     postp(1) % wobse (34)    = 'FPORZ'
     postp(1) % wobse (35)    = 'ACCEL'   ! Acceleration 
     postp(1) % wobse (36)    = 'MOMEN'   ! Linear momentum 
     postp(1) % wobse (37)    = 'MOMEY'   ! Linear momentum 
     postp(1) % wobse (38)    = 'MOMEZ'   ! Linear momentum 

     postp(1) % wonse (1)     = 'VELOX'   ! x-velocity
     postp(1) % wonse (2)     = 'VELOY'   ! y-velocity
     postp(1) % wonse (3)     = 'VELOZ'   ! z-velocity
     postp(1) % wonse (4)     = 'PRESS'   ! Pressure
     postp(1) % wonse (5)     = 'YPLUS'   ! y+
     postp(1) % wonse (6)     = 'PMV  '   ! PMV - comfort index
     postp(1) % wonse (7)     = 'PPD  '   ! PPD - comfort index
     !
     ! Witness variables
     !
     postp(1) % wowit ( 1)    = 'VELOX'   ! x-velocity
     postp(1) % wowit ( 2)    = 'VELOY'   ! y-velocity
     postp(1) % wowit ( 3)    = 'VELOZ'   ! z-velocity
     postp(1) % wowit ( 4)    = 'PRESS'   ! Pressure
     postp(1) % wowit ( 5)    = 'S11  '   ! Strain rate
     postp(1) % wowit ( 6)    = 'S22  '   ! Strain rate
     postp(1) % wowit ( 7)    = 'S12  '   ! Strain rate
     postp(1) % wowit ( 8)    = 'S33  '   ! Strain rate
     postp(1) % wowit ( 9)    = 'S13  '   ! Strain rate
     postp(1) % wowit (10)    = 'S23  '   ! Strain rate
     postp(1) % wowit (11)    = 'COORX'   ! x-coordinate
     postp(1) % wowit (12)    = 'COORY'   ! y-coordinate
     postp(1) % wowit (13)    = 'COORZ'   ! z-coordinate
     postp(1) % wowit (14)    = 'TURBU'   ! Turbulent viscosity
     postp(1) % wowit (15)    = 'VELMX'   ! x-velocity mesh
     postp(1) % wowit (16)    = 'VELMY'   ! y-velocity mesh
     postp(1) % wowit (17)    = 'VELMZ'   ! z-velocity mesh
     !
     ! Solvers
     !
     call soldef(-9_ip)                   ! Allocate memory
     solve(1) % kfl_solve = 1             ! Momentum:   Output flag
     solve(2) % kfl_solve = 1             ! Continuity: Output flag
     !
     ! Pointers: matrices used in Schur complement methods
     !
     Auu_nsi => nulir
     Aup_nsi => nulir
     Apu_nsi => nulir
     App_nsi => nulir
     !
     ! Nullify pointers
     !
     nullify(lforc_material_nsi)
     nullify(xforc_material_nsi)
     nullify(velta_nsi)
     nullify(thrta_nsi)
     nullify(powta_nsi)
     nullify(veave_nsi)
     nullify(radiu_nsi)
     nullify(forcn_nsi)
     nullify(forct_nsi)
     nullify(kfl_fixno_nsi)
     nullify(kfl_fixpr_nsi)
     nullify(kfl_fixpp_nsi)
     nullify(kfl_fixbo_nsi)
     nullify(kfl_fixrs_nsi)
     nullify(kfl_funno_nsi)
     nullify(kfl_funbo_nsi)
     nullify(kfl_funtn_nsi)
     nullify(kfl_funtb_nsi)
     nullify(kfl_wlawf_nsi)
     nullify(bvess_nsi)
     nullify(bvnat_nsi)
     nullify(skcos_nsi)
     nullify(tncod_nsi)
     nullify(tgcod_nsi)
     nullify(tbcod_nsi)
     nullify(tscod_nsi)
     nullify(veold_nsi)
     nullify(gradv_nsi)
     nullify(unk2n_nsi)
     nullify(bpess_nsi)
     nullify(dunkn_nsi)
     nullify(dunkp_nsi)
     nullify(avvel_nsi)
     nullify(avve2_nsi)
     nullify(avvxy_nsi)
     nullify(avpre_nsi)
     nullify(avpr2_nsi)
     nullify(avtan_nsi)
     nullify(resch_nsi)
     nullify(prope_nsi)
     nullify(norle_nsi)
     nullify(curle_nsi)
     nullify(outflow_mass)
     nullify(intfo_nsi)
     nullify(iboun_nsi)
     nullify(Auu_nsi)
     nullify(Aup_nsi)
     nullify(Apu_nsi)
     nullify(App_nsi)
     nullify(amatr_nsi)
     nullify(lapla_nsi)
     nullify(visco_nsi)
     nullify(cmama_nsi)
     nullify(deltp_nsi)
     nullify(dt_rho_nsi)
     nullify(mass_rho_nsi)
     nullify(tau_nsi)
     nullify(bubble_nsi)
     nullify(bubble_aqq_nsi)
     nullify(bubble_aqu_nsi)
     nullify(bubble_aqp_nsi)
     nullify(bubble_bq_nsi)
     nullify(vepro_nsi)
     nullify(grpro_nsi)
     nullify(prpro_nsi)
     nullify(vepr2_nsi)
     nullify(grpr2_nsi)
     nullify(prpr2_nsi)
     nullify(kfl_fixno_div_nsi)
     nullify(notra_nsi)
     nullify(massb_nsi)
     nullify(avmut_nsi)
     nullify(avstx_nsi)
     nullify(avsty_nsi)
     nullify(avstz_nsi)
     nullify(envel_nsi)
     nullify(enve2_nsi)
     nullify(envxy_nsi)
     nullify(enpre_nsi)
     nullify(enpr2_nsi)
     nullify(entan_nsi)
     nullify(enmut_nsi)
     nullify(enstx_nsi)
     nullify(ensty_nsi)
     nullify(enstz_nsi)
     nullify(btrac_nsi)
     nullify(tracr_nsi)
     nullify(tluav_nsi)
     nullify(vafor_nsi)
     nullify(avvaf_nsi)
     nullify(bupor_nsi)
     nullify(avntr_nsi)
     nullify(avgtr_nsi)
     nullify(drhodt_nsi)
     nullify(prdivcor_nsi)
     nullify(avmos_nsi)
     nullify(av_mass_flux_nsi)
     nullify(av_mom_flux_diag_nsi)
     nullify(av_mom_flux_off_nsi)
     nullify(fsifo_nsi)
     nullify(vefix)

     nullify(resdiff_nsi)
     nullify(dcost_dx_nsi)
     nullify(costdiff_nsi)
     
     nullify(lbpse)
     !
     ! Flow rates
     !
     nullify(kfl_flow_rate_codes_nsi)  ! Flow rates codes
     nullify(kfl_flow_rate_set_nsi)    ! Flow rates set to compute pressure
     nullify(kfl_flow_rate_normal_nsi) ! Flow rates normals 
     nullify(kfl_flow_rate_stfn_nsi)   ! Flow rates space time function number       
     nullify(kfl_flow_rate_tfn_nsi)    ! Flow rates time function number       
     nullify(kfl_flow_rate_pfn_nsi)    ! Flow rates discrete function number       
     nullify(kfl_flow_rate_alg_nsi)    ! Flow rates algorithm
     nullify(flow_rate_values_nsi)     ! Flow rates values
     nullify(flow_rate_press_nsi)      ! Flow rates imposed pressure values
     nullify(flow_rate_relax_nsi)      ! Flow rates temporal_relaxation
     nullify(flow_rate_normals_nsi)    ! Flow rates values
     mflow_nsi = ncodb
     if( mflow_nsi > 0 ) then
        allocate(kfl_flow_rate_codes_nsi (mflow_nsi))
        allocate(kfl_flow_rate_set_nsi   (mflow_nsi))
        allocate(kfl_flow_rate_normal_nsi(mflow_nsi))
        allocate(kfl_flow_rate_stfn_nsi  (mflow_nsi))
        allocate(kfl_flow_rate_tfn_nsi   (mflow_nsi))
        allocate(kfl_flow_rate_pfn_nsi   (mflow_nsi))
        allocate(kfl_flow_rate_alg_nsi   (mflow_nsi))
        allocate(flow_rate_values_nsi    (mflow_nsi))
        allocate(flow_rate_press_nsi     (mflow_nsi)) 
        allocate(flow_rate_relax_nsi     (mflow_nsi)) 
        allocate(flow_rate_normals_nsi   (ndime,mflow_nsi))
     end if
     !
     ! Others
     !
     nodpr_nsi           = 0      ! Node where to impose pressure
     kfl_exist_fixi7_nsi = 0      ! exists nodal fixity of type 7
     kfl_exist_fib20_nsi = 0      ! exists boundary fixity of type 20
     pabdf_nsi           = 0.0_rp
     dtinv_nsi           = 0.0_rp
     press_lev_nsi       = 0.0_rp ! Pressure level for coupling and postprocess
     
  case(1_ip)
     !
     ! Dimensions
     !
     ivari_nsi = ivari_nsi_mom
     !
     ! Number of velocity components
     !
     if(kfl_timei_nsi==1) then
        if(kfl_tisch_nsi==1) then
           !
           ! Trapezoidal rule
           !
           ncomp_nsi=3
        else if(kfl_tisch_nsi==2) then
           !
           ! BDF scheme
           !
           ncomp_nsi=2+kfl_tiacc_nsi
        else if(kfl_tisch_nsi==3 .or. kfl_tisch_nsi==4) then !=3 adams, =4 RK
           !
           ! Explicit time step
           !
           ncomp_nsi = 4
           ! AB WHAT IS THIS TODO TODO TODO if(kfl_tisch_nsi==4 .and. kfl_fscon_nsi == 0) then !RK extrapolating pressure
           ! AB WHAT IS THIS TODO TODO TODO    ncomp_nsi=4   ! only needed for the pressure unknown 
           ! AB WHAT IS THIS TODO TODO TODO else
           ! AB WHAT IS THIS TODO TODO TODO    ncomp_nsi=3 
           ! AB WHAT IS THIS TODO TODO TODO end if
        end if
        if ( kfl_ini_ts_guess_order_nsi > 1 ) then
           ncomp_nsi = max (ncomp_nsi,kfl_ini_ts_guess_order_nsi + 2 )
        end if
     else
        ncomp_nsi=2
     end if
     nprev_nsi = min(3_ip,ncomp_nsi)  ! Last time step or global iteration
     nbdfp_nsi = 2
     !
     ! First or second order velocity derivatives
     !
     if( kfl_dttyp_nsi == 2 .or. output_postprocess_check_variable_postprocess(arrays_number('DUDT ')) .or.&
          &                      output_postprocess_check_variable_postprocess(arrays_number('D2UDT')) ) then
        ncomp_nsi = max(5_ip,ncomp_nsi)
     end if
     !
     ! Time variables
     !
     if( kfl_timei_nsi == 0 ) then
        dtinv_nsi = 1.0_rp
     else
        kfl_timei = 1
     end if
     kfl_stead_nsi = 0
     !
     ! Dimensions
     !
     if( NSI_MONOLITHIC ) then
        solve(1) % ndofn = ndime + 1
        solve(1) % wprob = 'NAVIER_STOKES'
     else
        solve(1) % ndofn = ndime
        solve(2) % ndofn = 1
        solve(1) % wprob = 'MOMENTUM'
        solve(2) % wprob = 'CONTINUITY'
     end if
     if( NSI_SEMI_IMPLICIT ) then
        solve(9) % kfl_solve = 1
        solve(9) % ndofn     = ndime
        solve(9) % wprob     = 'VISCOUS_TERM'        
     end if
     !
     ! Multiple solves using the same matrix
     !
     if( NSI_SCHUR_COMPLEMENT ) then
        solve(1) % num_multiple_solves = 2
        if( kfl_sosch_nsi == 3 ) then
           solve(2) % num_multiple_solves = 2
        end if
     end if
     !
     ! Off-diagonal part of momentum equations
     !
     kfl_rmom2_nsi = 0
     kfl_p1ve2_nsi = 0

     if( fvnoa_nsi     >  0.0_rp )                          kfl_rmom2_nsi = 1
     if( kfl_visco_nsi == 1 .and. fvins_nsi     >= 0.9_rp ) kfl_rmom2_nsi = 1
     if( kfl_advec_nsi == 1 .and. kfl_linea_nsi == 2      ) kfl_rmom2_nsi = 1
     if( kfl_advec_nsi == 1 .and. kfl_adj_prob  == 1      ) kfl_rmom2_nsi = 1
     if( fvnoa_nsi     >  0.0_rp )                          kfl_p1ve2_nsi = 1
     if( kfl_anipo_nsi /= 0      )                          kfl_rmom2_nsi = 1
     
     ndbgs_nsi = ndime * npoin

     if( kfl_initi_nsi == 5 ) then
        solve(3)             = solve(1)                ! Check boundary conditions (same solver as Momentum)
        solve(3) % wprob     = 'BOUNDARY_CONDITIONS'
        solve(3) % ndofn     = ndime
        solve(3) % kfl_solve = 0
        if( kfl_initi_nsi == 5 ) then
           solve(3) % kfl_solve = 1
           solve(3) % kfl_cvgso = 1
        end if
     end if
     !
     ! Mass correction solver
     !
     solve(4) % wprob     = 'MASS_CORRECTION'
     if( kfl_corre_nsi == 3 ) then
        solve(4) % ndofn     = ndime
        solve(4) % kfl_solve = 1                        ! Output flag
     else
        solve(4) % kfl_algso = SOL_SOLVER_RICHARDSON    ! Mass correction: Diagonal solver
        solve(4) % ndofn     = ndime
        solve(4) % xdiag     = 1.0_rp
        if( kfl_corre_nsi == 1 ) then
           solve(4) % kfl_preco = SOL_CLOSE_MASS_MATRIX ! Uses VMASC
        else
           solve(4) % kfl_preco = SOL_MASS_MATRIX       ! Uses VMASS
        end if
        solve(4) % kfl_solve = 0
        solve(4) % kfl_cvgso = 0
     end if

     if( kfl_inipr_nsi == NSI_PDE_HYDROSTATIC_PRESSURE ) then
        solve(5) % wprob     = 'HYDROSTATIC_PRESSURE'
        solve(5) % ndofn     = 1
        solve(5) % kfl_solve = 1         
     end if

     if( kfl_divcorrec_nsi /= 0 ) then
        solve(6) % wprob     = 'ZERO_DIVERGENCE'
        solve(6) % ndofn     = 1
        solve(6) % kfl_solve = 1        
     end if

     solve(7) % kfl_solve = 0
     if( 1 == 0 ) then
        solve(7) % wprob     = 'NORMAL_EXTENSIONS'
        solve(7) % ndofn     = 1
        solve(7) % kfl_solve = 1 
     end if
     
     solve(8) % wprob     = 'CONSISTENT'                  ! Consistent matrix
     solve(8) % kfl_solve = 1
     solve(8) % kfl_iffix = 2
     solve(8) % ndofn     = ndime

     nevat_nsi          = (ndime+1)*mnode
     nschu_nsi          = 0                               ! # Schur complement solves
     nmome_nsi          = 0                               ! # Momentum solves

     kfl_perip_nsi      = 0                               ! Periodicity: if pressure prescribed on periodic nodes
     pcoef_nsi          = 1.0_rp-gasco/sphea_nsi
     gamth_nsi          = sphea_nsi/(sphea_nsi-gasco)     ! Low-Mach: gamma=Cp/(Cp-R)

     prthe(1)           = prthe_nsi                       ! Low-Mach: p0
     prthe(2)           = prthe_nsi                       ! Low-Mach: p0^{n}
     prthe(3)           = prthe_nsi                       ! Low-Mach: p0^{n-1}
     prthe(4)           = prthe_nsi                       ! Low-Mach: p0^{0}

     imass              = tmass_nsi                       ! initial density
     kfl_prthe          = kfl_prthe_nsi                   ! type of thpressure calculation
     tflux              = 0.0_rp                          ! Low-Mach: heat flux
     dpthe              = 0.0_rp                          ! Low-Mach: dp0/dt
     xmass_nsi          = 0.0_rp                          ! Low-Mach: mass
     cputi_assembly_nsi = 0.0_rp
     cpu_ass_sol_nsi    = 0.0_rp
     iteqn_nsi(1)       = 0
     iteqn_nsi(2)       = 0
     ittot_nsi          = 0
     resin_nsi(1)       = 1.0_rp                          ! Algebraic inner residual
     resin_nsi(2)       = 1.0_rp
     resou_nsi(1)       = 1.0_rp                          ! Algebraic outer residual
     resou_nsi(2)       = 1.0_rp
     resss_nsi(1)       = 1.0_rp                          ! Algebraic Steady state residual
     resss_nsi(2)       = 1.0_rp
     reinf_nsi(1)       = 1.0_rp                          ! Linf residual
     reinf_nsi(2)       = 1.0_rp
     relpa_nsi(1)       = 0.0_rp                          ! Relaxation
     relpa_nsi(2)       = 0.0_rp
     corio_nsi          = 0.0_rp                          ! Coriolis force module
     kfl_tiaor_nsi      = kfl_tiacc_nsi                   ! Time accuracy: save original value
     resgs_nsi(1)       = 0.0_rp
     resgs_nsi(2)       = 0.0_rp
     rmsgs_nsi          = 0.0_rp
     dtsgs_nsi          = 0.0_rp
     kfl_wlare_nsi      = 0
     porfo_nsi          = 0.0_rp
     dtmax_nsi          = 0.0_rp                          ! otherwise full ckeck stops because it is unitialized when doing max in nsi_tistep
     actav_nsi          = 0.0_rp                          ! Accumulated time for averaging
     if(kfl_sgsco_nsi==0) misgs_nsi=1                     ! Subgrid scale number of iterations
     if( kfl_sgsac_nsi == -1 ) then                       ! SGS time accuracy default
        if(kfl_tisch_nsi/=2) then
           kfl_sgsac_nsi = kfl_tiacc_nsi
        else
           kfl_sgsac_nsi = 1_ip                           ! For BDF use order 1 for the subscale by default
        end if
     end if
     !
     ! Performance timings
     !
     call moduls_allocate_timing(11_ip)
     call times( 1) % register('schur complement'   ,'iteration other')
     call times( 2) % register('wall exchange'      ,'iteration other')
     call times( 3) % register('subgrid scale'      ,'iteration other')
     call times( 4) % register('velocity correction','iteration other')
     call times( 5) % register('multistep step'     ,'iteration other')
     call times( 6) % register('semi-implicit'      ,'iteration other')
     call times( 7) % register('time step'          ,'step beginning' )
     call times( 8) % register('FS matrices'        ,'step beginning' )

     call times( 9) % register('gather'             ,'element assembly' )
     call times(10) % register('computation'        ,'element assembly' )
     call times(11) % register('scatter'            ,'element assembly' )
     
     call times( 9) % ini() ; call times( 9) % add()
     call times(10) % ini() ; call times(10) % add()
     call times(11) % ini() ; call times(11) % add()

     !
     ! CFD wake
     !
#ifndef NASTIN_PRIVATE_OFF
     call nsi_cfdwak(1_ip)
#endif
#ifdef matiaslma
     !
     ! for low Mach, never solve pressure
     !
     if (kfl_regim_nsi==3) then
        kfl_confi_nsi = -1   ! do not fix pressure!!!
        fvins_nsi = 2.0_rp     ! Complete form
     end if
     !  pressure stabilization coefficient
     staco_nsi(4) = min(staco_nsi(4), 0.25_rp/staco_nsi(1))
#endif
     
  case(2_ip)
     !
     ! If temperature subgrid scale should be considered
     !
     if( associated(tesgs) .and. (kfl_cotem_nsi == 1 .or. kfl_regim_nsi == 3) ) then !!FER Check if SGS of temper is needed for low mach
        kfl_sgste_nsi = 1
     else
        kfl_sgste_nsi = 0
     end if
     !
     ! Activate flag for Low Mach in combustion (properties on GP or Node)
     !
     kfl_lookg_nsi=0
     if (kfl_regim_nsi==3 .and. (associated(tempe_gp) .or. kfl_coupl(ID_NASTIN,ID_CHEMIC) /= 0)) kfl_lookg_nsi = 1

  case(3_ip)
     !
     ! Before starting a time step
     !
     relax_nsi = 1.0_rp
     relap_nsi = 1.0_rp

  end select

end subroutine nsi_inivar
