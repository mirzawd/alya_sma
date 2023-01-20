!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ker_memory(itask)
  !-----------------------------------------------------------------------
  !****f* Kermod/ker_memory
  ! NAME 
  !    ker_memcbs
  ! DESCRIPTION
  !    This routine allocates memory for the boundary conditions arrays
  ! USES
  !    ecoute
  ! USED BY
  !    ker_memory
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_kermod
  use def_domain
  use mod_ker_memory, only : memory_alloca
  use mod_ker_memory, only : memory_deallo
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ifunc,idime

  select case (itask)

  case( 1_ip)
     !
     ! KFL_FIXNO_WALLD_KER: Allocate kfl_fixno_walld_ker
     !
     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXNO_WALLD_KER','ker_memory' , kfl_fixno_walld_ker , 1_ip , npoin )

  case(-1_ip)
     !
     ! KFL_FIXNO_WALLD_KER: Dellocate kfl_fixno_walld_ker
     !
     call memory_deallo(mem_modul(1:2,modul),'KFL_FIXNO_WALLD_KER','ker_memory' , kfl_fixno_walld_ker )

  case( 2_ip)
     !
     ! KFL_FIXBO_WALLD_KER: Allocate kfl_fixbo_walld_ker
     !
     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXBO_WALLD_KER','ker_memory' , kfl_fixbo_walld_ker , nboun ) 

  case(-2_ip)
     !
     ! KFL_FIXBO_WALLD_KER: Dellocate kfl_fixbo_walld_ker
     !
     call memory_deallo(mem_modul(1:2,modul),'KFL_FIXBO_WALLD_KER','ker_memory' , kfl_fixbo_walld_ker )

  case( 3_ip)
     !
     ! KFL_FIXNO_ROUGH_KER: Allocate kfl_fixno_rough_ker
     !
     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXNO_ROUGH_KER','ker_memory' , kfl_fixno_rough_ker , 1_ip , npoin )

  case(-3_ip)
     !
     ! KFL_FIXNO_ROUGH_KER: Dellocate kfl_fixno_rough_ker
     !
     call memory_deallo(mem_modul(1:2,modul),'KFL_FIXNO_ROUGH_KER','ker_memory' , kfl_fixno_rough_ker )

  case( 4_ip)
     !
     ! KFL_FIXBO_ROUGH_KER: Allocate kfl_fixbo_rough_ker
     !
     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXBO_ROUGH_KER','ker_memory' , kfl_fixbo_rough_ker , nboun )

  case(-4_ip)
     !
     ! KFL_FIXBO_ROUGH_KER: Dellocate kfl_fixbo_rough_ker
     !
     call memory_deallo(mem_modul(1:2,modul),'KFL_FIXBO_ROUGH_KER','ker_memory' , kfl_fixbo_rough_ker )

  case( 5_ip)
     !
     ! Support geometry
     !
     call memory_alloca(mem_modul(1:2,modul),'COORD_MM','ker_memory' , coord_mm , ndime    , npoin_mm )
     call memory_alloca(mem_modul(1:2,modul),'LNODB_MM','ker_memory' , lnodb_mm , mnodb_mm , nboun_mm )
     call memory_alloca(mem_modul(1:2,modul),'LTYPB_MM','ker_memory' , ltypb_mm , nboun_mm )

  case(-5_ip)
     !
     ! Support geometry
     !
     call memory_deallo(mem_modul(1:2,modul),'COORD_MM','ker_memory' , coord_mm )
     call memory_deallo(mem_modul(1:2,modul),'LNODB_MM','ker_memory' , lnodb_mm )
     call memory_deallo(mem_modul(1:2,modul),'LTYPB_MM','ker_memory' , ltypb_mm )

  case( 6_ip)
     !
     ! Support geometry boundary conditions: KFL_FIXNO_SUPPO_KER and BVESS_SUPPO_KER
     !
     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXNO_SUPPO_KER','ker_memory' , kfl_fixno_suppo_ker , ndime, npoin )
     call memory_alloca(mem_modul(1:2,modul),'BVESS_SUPPO_KER'    ,'ker_memory' , bvess_suppo_ker     , ndime, npoin )

  case(-6_ip)
     !
     ! Support geometry boundary conditions: KFL_FIXNO_SUPPO_KER and BVESS_SUPPO_KER
     !
     call memory_deallo(mem_modul(1:2,modul),'KFL_FIXNO_SUPPO_KER','ker_memory' , kfl_fixno_suppo_ker )
     call memory_deallo(mem_modul(1:2,modul),'BVESS_SUPPO_KER'    ,'ker_memory' , bvess_suppo_ker     )

  case( 7_ip)
     !
     ! Space time functions
     !
     ifunc = igene
     idime = space_time_function(ifunc) % ndime
     call memory_alloca(mem_modul(1:2,modul),'SPACE_TIME_FUNCTION','ker_memory' , space_time_function(ifunc) % expression,idime )

  case(-7_ip)
     !
     ! Space time functions
     !
     ifunc = igene
     idime = space_time_function(ifunc) % ndime
     call memory_deallo(mem_modul(1:2,modul),'SPACE_TIME_FUNCTION','ker_memory' , space_time_function(ifunc) % expression )

  case( 8_ip)
     !
     ! PARAMETERS_TIME_FUNCTION
     !
     ifunc = igene
     call memory_alloca(mem_modul(1:2,modul),'PARAMETERS_TIME_FUNCTION','ker_memory',&
          time_function(ifunc) % parameters,time_function(ifunc) % npara)

  case(-8_ip)
     !
     ! PARAMETERS_TIME_FUNCTION
     !
     ifunc = igene
     call memory_deallo(mem_modul(1:2,modul),'PARAMETERS_TIME_FUNCTION','ker_memory',&
          time_function(ifunc) % parameters)

  case( 9_ip)
     !
     ! KFL_FIXNO_DEFOR_KER
     ! BVESS_DEFOR_KER
     !
     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXNO_DEFOR_KER','ker_memory' , kfl_fixno_defor_ker , ndime , npoin )
     call memory_alloca(mem_modul(1:2,modul),'BVESS_DEFOR_KER'    ,'ker_memory' , bvess_defor_ker     , ndime , npoin )

  case(10_ip)
     !
     ! KFL_FIXNO_WALLN_KER: Allocate kfl_fixno_walln_ker
     !
     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXNO_WALLN_KER','ker_memory' , kfl_fixno_walln_ker , 1_ip , npoin )
     call memory_alloca(mem_modul(1:2,modul),'BVESS_WALLN_KER'    ,'ker_memory' , bvess_walln_ker     , 1_ip , npoin )

  case(-10_ip)
     !
     ! KFL_FIXNO_WALLN_KER: Dellocate kfl_fixno_walln_ker
     !
     call memory_deallo(mem_modul(1:2,modul),'KFL_FIXNO_WALLN_KER','ker_memory' , kfl_fixno_walln_ker )
     call memory_deallo(mem_modul(1:2,modul),'BVESS_WALLN_KER'    ,'ker_memory' , bvess_walln_ker )

  case(11_ip)
     !
     ! KFL_FIXBO_WALLN_KER: Allocate kfl_fixbo_walln_ker
     !
     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXBO_WALLN_KER','ker_memory' , kfl_fixbo_walln_ker , nboun ) 

  case(-11_ip)
     !
     ! KFL_FIXBO_WALLN_KER: Dellocate kfl_fixbo_walln_ker
     !
     call memory_deallo(mem_modul(1:2,modul),'KFL_FIXBO_WALLN_KER','ker_memory' , kfl_fixbo_walln_ker )
     
  case(12_ip)
     !
     ! KFL_FIXBO_NSW_KER and LNSW_EXCH: Allocate kfl_fixbo_nsw_ker
     !
     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXBO_NSW_KER','ker_memory' , kfl_fixbo_nsw_ker , max(1_ip,nboun) )
     call memory_alloca(mem_modul(1:2,modul),'LNSW_EXCH','ker_inibcs'         , lnsw_exch,max(1_ip,nelem))

  case(-12_ip)
     !
     ! KFL_FIXBO_NSW_KER and LNSW_EXCH: Dellocate kfl_fixbo_walln_ker
     !
     call memory_deallo(mem_modul(1:2,modul),'KFL_FIXBO_NSW_KER','ker_memory' , kfl_fixbo_nsw_ker )
     call memory_deallo(mem_modul(1:2,modul),'LNSW_EXCH','ker_inibcs'         , lnsw_exch)
     
  case( 13_ip)
     !
     ! PARAMETERS_WINDK_SYSTEMS
     !
     ifunc = igene
     call memory_alloca(mem_modul(1:2,modul),'PARAMETERS_WINDK_SYSTEMS','ker_memory',&
          windk_systems(ifunc) % params,windk_systems(ifunc) % nparam)
     call memory_alloca(mem_modul(1:2,modul),'XPREV_WINDK_SYSTEMS','ker_memory',&
          windk_systems(ifunc) % xprev,windk_systems(ifunc) % ndxs)
     call memory_alloca(mem_modul(1:2,modul),'YPREV_WINDK_SYSTEMS','ker_memory',&
          windk_systems(ifunc) % yprev,windk_systems(ifunc) % ndxs)

          windk_systems(ifunc) % params = 0.0_rp
          windk_systems(ifunc) % xprev = 0.0_rp
          windk_systems(ifunc) % yprev = 0.0_rp

  case(-13_ip)
     !
     ! PARAMETERS_WINDK_SYSTEMS
     !
     ifunc = igene
     ! params should stay allocated
     !call memory_deallo(mem_modul(1:2,modul),'PARAMETERS_WINDK_SYSTEMS','ker_memory',&
     !     windk_systems(ifunc) % params)
     call memory_deallo(mem_modul(1:2,modul),'XPREV_WINDK_SYSTEMS','ker_memory',&
          windk_systems(ifunc) % xprev)
     call memory_deallo(mem_modul(1:2,modul),'YPREV_WINDK_SYSTEMS','ker_memory',&
          windk_systems(ifunc) % yprev)


  case( 14_ip)
     !
     ! PARAMETERS_PUMP_CURVES
     !
     ifunc = igene
     call memory_alloca(mem_modul(1:2,modul),'PARAMETERS_PUMP_CURVES','ker_memory',&
          pump_curve(ifunc) % params,pump_curve(ifunc) % nparam)
          pump_curve(ifunc) % params = 0.0_rp

  case(-14_ip)
     !
     ! PARAMETERS_BOMB_CURVES
     !
     ifunc = igene
     ! params should stay allocated
     !call memory_deallo(mem_modul(1:2,modul),'PARAMETERS_BOMB_CURVES','ker_memory',&
     !     pump_curve(ifunc) % params)

  end select

end subroutine ker_memory
