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
!> @}
!-----------------------------------------------------------------------
subroutine gus_inivar(itask)
  
  use def_kintyp_basic, only : ip,rp
  use def_master,       only : solve
  use def_master,       only : postp
  use def_master,       only : veloc
  use def_master,       only : press
  use def_master,       only : vel1d
  use def_master,       only : areas
  use def_master,       only : flowr
  use mod_arrays,       only : arrays_register
  use def_gusano
  
  implicit none
  
  integer(ip), intent(in) :: itask

  select case ( itask )

  case ( 0_ip )
     !
     ! Module varianle
     !
     call arrays_register((/'VEL1D','SCALA','NPOIN','SECON'/),vel1d    ,ENTITY_POSITION=1_ip,TIME_POSITION=2_ip)
     call arrays_register((/'VELOC','VECTO','NPOIN','SECON'/),veloc    ,ENTITY_POSITION=2_ip,TIME_POSITION=3_ip)
     call arrays_register((/'PRESS','SCALA','NPOIN','PRIMA'/),press    ,ENTITY_POSITION=1_ip,TIME_POSITION=2_ip,EXCLUDE_RST=(/2_ip,3_ip/))
     call arrays_register((/'PROJM','SCALA','NPOIN','SECON'/),projm_gus,ENTITY_POSITION=1_ip)
     call arrays_register((/'ANGLE','SCALA','NPOIN','SECON'/),angle_gus,ENTITY_POSITION=1_ip)
     call arrays_register((/'AREAS','SCALA','NPOIN','SECON'/),areas    ,ENTITY_POSITION=1_ip,TIME_POSITION=2_ip,EXCLUDE_RST=(/2_ip,3_ip/))
     call arrays_register((/'VEL3D','VECTO','NPOIN','SECON'/),vel3d_gus,ENTITY_POSITION=2_ip)
     call arrays_register((/'PRE3D','SCALA','NPOIN','SECON'/),pre3d_gus,ENTITY_POSITION=1_ip)
     call arrays_register((/'FLOWR','SCALA','NPOIN','PRIMA'/),flowr    ,ENTITY_POSITION=1_ip,TIME_POSITION=2_ip,EXCLUDE_RST=(/2_ip,3_ip/))
     call arrays_register((/'PROJC','SCALA','NPOIN','SECON'/),projc_gus,ENTITY_POSITION=1_ip)
     call arrays_register((/'LINEL','SCALA','NPOIN','SECON'/),          ENTITY_POSITION=1_ip)
      !
     ! Sets
     !
     postp(1) % wobse (1) = 'PRESS'   ! Pressure
     postp(1) % wobse (2) = 'FLOWR'   ! Flow rate A*u.n     
     !
     ! Nullify
     !
     nullify(kfl_fixno_gus)
     nullify(kfl_fixbo_gus)
     nullify(kfl_funno_gus)
     nullify(kfl_funbo_gus)
     nullify(kfl_funtn_gus)
     nullify(kfl_funtb_gus)
     nullify(kfl_fixsc_gus)
     nullify(bvess_gus)
     nullify(bvnat_gus)
     nullify(projm_gus)
     nullify(projc_gus)
     nullify(neuman_gus)
     nullify(dirich_gus)
     nullify(exn1d_gus)
     nullify(angle_gus)
     nullify(bendi_gus)
     nullify(densi_gus)
     nullify(visco_gus)
     nullify(vel3d_gus)
     nullify(pre3d_gus)
     nullify(schur_gus)
     !
     ! Initialize
     !
     dtinv_gus     = 0.0_rp
     dtmax_gus     = 0.0_rp
     dtcri_gus     = 0.0_rp
     kfl_stead_gus = 0
     !
     ! Solvers
     !
     call soldef(-3_ip)                   ! Allocate memory
     solve(1) % kfl_solve = 1             ! Momentum+continuity:   Output flag
     
  case ( 1_ip )
     !
     ! Counters
     !
     kfl_stead_gus = 0
     kfl_goite_gus = 0
     dtinv_gus     = 0.0_rp
     !
     ! Stabilization
     !
     select case ( kfl_stabi_gus )
     case ( GUS_ASGS      ) ; xsoss_gus = 0.0_rp ; xoss_gus  = 0.0_rp
     case ( GUS_OSS       ) ; xsoss_gus = 1.0_rp ; xoss_gus  = 1.0_rp
     case ( GUS_SPLIT_OSS ) ; xsoss_gus = 1.0_rp ; xoss_gus  = 0.0_rp
     case default           ; call runend('GUS_INIVAR: UNKNOWN STABILIZATION STRATEGY')
     end select
     !
     ! Solvers
     !
     solve(1) % ndofn = 2
     solve(1) % wprob = '1D-NAVIER-STOKES'        
     if( kfl_algor_gus == GUS_SCHUR_COMPLEMENT ) then
        solve(2) % kfl_solve = 1             ! Momentum:   Output flag
        solve(3) % kfl_solve = 1             ! Continuity: Output flag
        solve(2) % ndofn     = 1
        solve(2) % wprob     = 'MOMENTUM'        
        solve(3) % ndofn     = 1
        solve(3) % wprob     = 'CONTINUITY'        
     end if
     
  end select
  
end subroutine gus_inivar
