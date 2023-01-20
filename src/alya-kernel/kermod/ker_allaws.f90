!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ker_allaws()

  !-----------------------------------------------------------------------
  !****f* kermod/ker_allaws
  ! NAME
  !   ker_allaws
  !   lresp = -1 ... Update property always
  ! DESCRIPTION
  !   Define laws
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_master
  use def_kermod
  use mod_ker_polynomial, only : ker_polynomial_allaws !< 2016JUN29
  implicit none
  !
  ! Names of properties
  !
  densi_ker % name = 'DENSI'
  poros_ker % name = 'POROS'
  visco_ker % name = 'VISCO'
  condu_ker % name = 'CONDU'
  sphea_ker % name = 'SPHEA'
  dummy_ker % name = 'DUMMY'
  turmu_ker % name = 'TURMU'
  absor_ker % name = 'ABSOR'
  scatt_ker % name = 'SCATT'
  mixin_ker % name = 'MIXIN'
  anipo_ker % name = 'ANIPO'  
  walvi_ker % name = 'WALLV'  
  !
  ! Law 1 is constant for all properties
  !
  densi_ker % llaws(1) % wname = 'CONST'
  poros_ker % llaws(1) % wname = 'CONST'
  visco_ker % llaws(1) % wname = 'CONST'
  condu_ker % llaws(1) % wname = 'CONST'
  sphea_ker % llaws(1) % wname = 'CONST'
  dummy_ker % llaws(1) % wname = 'CONST'
  turmu_ker % llaws(1) % wname = 'CONST'
  absor_ker % llaws(1) % wname = 'CONST'
  scatt_ker % llaws(1) % wname = 'CONST'
  mixin_ker % llaws(1) % wname = 'CONST'
  anipo_ker % llaws(1) % wname = 'CONST'
  walvi_ker % llaws(1) % wname = 'CONST'

  densi_ker % llaws(1) % lresp = -2
  poros_ker % llaws(1) % lresp = -2
  visco_ker % llaws(1) % lresp = -2
  condu_ker % llaws(1) % lresp = -2
  sphea_ker % llaws(1) % lresp = -2
  dummy_ker % llaws(1) % lresp = -2
  turmu_ker % llaws(1) % lresp = -2
  absor_ker % llaws(1) % lresp = -2
  scatt_ker % llaws(1) % lresp = -2
  mixin_ker % llaws(1) % lresp = -2
  anipo_ker % llaws(1) % lresp = -2
  walvi_ker % llaws(1) % lresp = -2

  densi_ker % llaws(1) % where = 'CONST'
  poros_ker % llaws(1) % where = 'CONST'
  visco_ker % llaws(1) % where = 'CONST'
  condu_ker % llaws(1) % where = 'CONST'
  sphea_ker % llaws(1) % where = 'CONST'
  dummy_ker % llaws(1) % where = 'CONST'
  turmu_ker % llaws(1) % where = 'CONST'
  absor_ker % llaws(1) % where = 'CONST'
  scatt_ker % llaws(1) % where = 'CONST'
  mixin_ker % llaws(1) % where = 'CONST'
  anipo_ker % llaws(1) % where = 'CONST'
  walvi_ker % llaws(1) % where = 'IELEM'

  !----------------------------------------------------------------------
  !
  ! Density
  !
  !----------------------------------------------------------------------

  densi_ker % llaws(2) % wname     = 'BIFLU'       ! Bifluid
  densi_ker % llaws(2) % lresp(1)  =  ID_LEVELS  
  densi_ker % llaws(2) % where     = 'IELEM'
  densi_ker % llaws(2) % kfl_gradi = 0             ! Despite there is a gradient we never use it for Biflid flow 

  densi_ker % llaws(3) % wname     = 'LOWMA'       ! Low-Mach rho = p/RT
  densi_ker % llaws(3) % lresp(1)  =  ID_TEMPER
  densi_ker % llaws(3) % lresp(2)  =  ID_NASTIN
  densi_ker % llaws(3) % where     = 'IPOIN'
  densi_ker % llaws(3) % kfl_gradi = 1

  densi_ker % llaws(4) % wname     = 'KLOWM'       ! Low Mach with mixture rho = pW/RT
  densi_ker % llaws(4) % lresp(1)  =  ID_TEMPER
  densi_ker % llaws(4) % lresp(2)  =  ID_NASTIN
  densi_ker % llaws(4) % lresp(3)  =  ID_CHEMIC
  densi_ker % llaws(4) % where     = 'IPOIN'
  densi_ker % llaws(4) % kfl_gradi = 1
 
  densi_ker % llaws(5) % wname     = 'MIXTU'       ! Mixture of constant density fluids
  densi_ker % llaws(5) % lresp(1)  =  ID_CHEMIC
  densi_ker % llaws(5) % where     = 'IPOIN'
  densi_ker % llaws(5) % kfl_gradi = 1

  densi_ker % llaws(6) % wname     = 'BIPHA'       ! Mixture of constant density fluids in two phases
  densi_ker % llaws(6) % lresp(1)  =  ID_CHEMIC
  densi_ker % llaws(6) % lresp(2)  =  ID_LEVELS
  densi_ker % llaws(6) % where     = 'IPOIN'
  densi_ker % llaws(6) % kfl_gradi = 1

  densi_ker % llaws(7) % wname     = 'TBIPH'       ! Mixture of temperature dependent density fluids in two phases
  densi_ker % llaws(7) % lresp(1)  =  ID_CHEMIC
  densi_ker % llaws(7) % lresp(2)  =  ID_LEVELS
  densi_ker % llaws(7) % lresp(3)  =  ID_TEMPER
  densi_ker % llaws(7) % where     = 'IPOIN'
  densi_ker % llaws(7) % kfl_gradi = 1
 
  densi_ker % llaws(8) % wname     = 'DNAST'       ! Test1
  densi_ker % llaws(8) % lresp(1)  =  ID_NASTAL
  densi_ker % llaws(8) % where     = 'IPOIN'
  densi_ker % llaws(8) % kfl_gradi = 1

  densi_ker % llaws(9) % wname     = 'TEST1'       ! Test1
  densi_ker % llaws(9) % lresp(1)  =  -1
  densi_ker % llaws(9) % where     = 'IELEM'
  densi_ker % llaws(9) % kfl_gradi = 1

  densi_ker % llaws(10) % wname     = 'TEST2'       ! Test2
  densi_ker % llaws(10) % lresp(1)  =  -1
  densi_ker % llaws(10) % where     = 'IPOIN'
  densi_ker % llaws(10) % kfl_gradi = 0
 
  densi_ker % llaws(11) % wname     = 'TEST3'       ! Test3 
  densi_ker % llaws(11) % lresp(1)  =  -1
  densi_ker % llaws(11) % where     = 'IPOIN'
  densi_ker % llaws(11) % kfl_gradi = 1

  densi_ker % llaws(12) % wname     = 'LOWMG'       ! Low mach in Gauss Points
  densi_ker % llaws(12) % lresp(1)  =  ID_TEMPER
  densi_ker % llaws(12) % where     = 'IELEM'
  densi_ker % llaws(12) % kfl_gradi = 0

  densi_ker % llaws(13) % wname     = 'TLOWM'       ! Low-Mach for syncronized CFI combustion model
  densi_ker % llaws(13) % lresp(1)  =  ID_TEMPER    ! density is only updated by temper
  densi_ker % llaws(13) % where     = 'IPOIN'
  densi_ker % llaws(13) % kfl_gradi = 1

  densi_ker % llaws(14) % wname     = 'GKLOW'       ! Low Mach with mixture rho = pW/RT at Gauss Points
  densi_ker % llaws(14) % lresp(1)  =  ID_TEMPER
  densi_ker % llaws(14) % where     = 'IELEM'
  densi_ker % llaws(14) % kfl_gradi = 1
  densi_ker % llaws(14) % kfl_deriv = 1
  densi_ker % llaws(14) % kfl_grder = 1

  call ker_polynomial_allaws(densi_ker, 15_ip, ID_TEMPER, 'IELEM', 0_ip)

  densi_ker % llaws(16) % wname     = 'SPRAY'       ! Liquid -gas density: 1/rho = Y_L / rho_L + (1-Y_L)/rho_g at Gauss Points
  densi_ker % llaws(16) % lresp(1)  =  ID_CHEMIC
  densi_ker % llaws(16) % where     = 'IELEM'
  densi_ker % llaws(16) % kfl_gradi = 0
  densi_ker % llaws(16) % kfl_deriv = 0
  densi_ker % llaws(16) % kfl_grder = 0

  densi_ker % llaws(17) % wname     = 'LMGPL'       ! Low Mach with mixture rho = pW/RT at Gauss Points (also lookup on GP)
  densi_ker % llaws(17) % lresp(1)  =  ID_CHEMIC
  densi_ker % llaws(17) % where     = 'IELEM'
  densi_ker % llaws(17) % kfl_gradi = 0
  densi_ker % llaws(17) % kfl_deriv = 0
  densi_ker % llaws(17) % kfl_grder = 0

  densi_ker % llaws(18) % wname     = 'DENGP'       ! Density at Gauss points
  densi_ker % llaws(18) % lresp(1)  =  ID_CHEMIC
  densi_ker % llaws(18) % where     = 'IELEM'
  densi_ker % llaws(18) % kfl_gradi = 0
  densi_ker % llaws(18) % kfl_deriv = 0
  densi_ker % llaws(18) % kfl_grder = 0  

  densi_ker % llaws(19) % wname     = 'DENNP'       ! Density at nodes
  densi_ker % llaws(19) % lresp(1)  =  ID_CHEMIC
  densi_ker % llaws(19) % lresp(2)  =  ID_NASTIN
  densi_ker % llaws(19) % where     = 'IPOIN'
  !!!!densi_ker % llaws(19) % kfl_gradi = 1
  !!!!densi_ker % llaws(19) % kfl_deriv = 0
  !!!!densi_ker % llaws(19) % kfl_grder = 0

  !----------------------------------------------------------------------
  !
  ! Viscosity
  !
  !----------------------------------------------------------------------

  visco_ker % llaws(2) % wname     = 'BIFLU'       ! Bifluid
  visco_ker % llaws(2) % lresp(1)  =  ID_LEVELS
  visco_ker % llaws(2) % where     = 'IELEM'
  visco_ker % llaws(2) % kfl_gradi = 0             ! Despite there is a gradient we never use it for Bifulid flow - see if in nsi_elmop3 grvis

  visco_ker % llaws(3) % wname     = 'SUTHE'       ! Sutherland
  visco_ker % llaws(3) % lresp(1)  =  ID_NASTAL
  visco_ker % llaws(3) % lresp(2)  =  ID_TEMPER
  visco_ker % llaws(3) % where     = 'IPOIN'
  visco_ker % llaws(3) % kfl_gradi = 1
  
  visco_ker % llaws(4) % wname     = 'MUMIX'       ! Mixture of species
  visco_ker % llaws(4) % lresp(1)  =  ID_CHEMIC
  visco_ker % llaws(4) % where     = 'IPOIN'
  visco_ker % llaws(4) % kfl_gradi = 1

  visco_ker % llaws(5) % wname     = 'TEST4'       ! Linear viscosity
  visco_ker % llaws(5) % lresp(1)  =  -1
  visco_ker % llaws(5) % where     = 'IPOIN'
  visco_ker % llaws(5) % kfl_gradi = 1

  visco_ker % llaws(6) % wname     = 'BIFL2'       ! Bifluid
  visco_ker % llaws(6) % lresp(1)  =  ID_LEVELS
  visco_ker % llaws(6) % where     = 'IELEM'
  visco_ker % llaws(6) % kfl_gradi = 1             ! BIFL2 includes gradient

  visco_ker % llaws(7) % wname     = 'GPSUT'       ! Sutherland in Gauss points
  visco_ker % llaws(7) % lresp(1)  =  ID_NASTAL
  visco_ker % llaws(7) % lresp(2)  =  ID_TEMPER
  visco_ker % llaws(7) % where     = 'IELEM'
  visco_ker % llaws(7) % kfl_gradi = 1
  visco_ker % llaws(7) % kfl_deriv = 1
  visco_ker % llaws(7) % kfl_grder = 1

  visco_ker % llaws(8) % wname     = 'ABL  '       ! ABL at Gauss points: mu = rho*kap/Cmu*(z+z0)*U_*
  visco_ker % llaws(8) % lresp(1)  =  ID_NASTIN
  visco_ker % llaws(8) % where     = 'IELEM'
  visco_ker % llaws(8) % kfl_gradi = 1

  visco_ker % llaws(9) % wname     = 'MUTAB'       ! from table (Flamelet combustion model)
  visco_ker % llaws(9) % lresp(1)  =  ID_CHEMIC
  visco_ker % llaws(9) % where     = 'IPOIN'
  visco_ker % llaws(9) % kfl_gradi = 1

  visco_ker % llaws(10) % wname     = 'TMUTA'      ! from table (Flamelet combustion model)
  visco_ker % llaws(10) % lresp(1)  =  ID_TEMPER   ! update only in temper
  visco_ker % llaws(10) % where     = 'IPOIN'
  visco_ker % llaws(10) % kfl_gradi = 1

  call ker_polynomial_allaws(visco_ker, 11_ip, ID_TEMPER, 'IELEM', 0_ip)

  visco_ker % llaws(12) % wname     = 'GAUSS'      ! Gaussian filter for generating turbulent fluctuations 
  visco_ker % llaws(12) % lresp(1)  =  ID_NASTIN   !  
  visco_ker % llaws(12) % where     = 'IELEM'
  visco_ker % llaws(12) % kfl_gradi = 1

  visco_ker % llaws(13) % wname     = 'SPRAY'      ! Viscosity of liquid/gas mixture constant 
  visco_ker % llaws(13) % lresp(1)  =  ID_CHEMIC   !  
  visco_ker % llaws(13) % where     = 'IELEM'
  visco_ker % llaws(13) % kfl_gradi = 0

  visco_ker % llaws(14) % wname     = 'MUGPT'       ! from table on Gauss points
  visco_ker % llaws(14) % lresp(1)  =  ID_CHEMIC
  visco_ker % llaws(14) % where     = 'IELEM'
  visco_ker % llaws(14) % kfl_gradi = 0

  visco_ker % llaws(15) % wname     = 'NAHME'      ! Nahme's approximation
  visco_ker % llaws(15) % lresp(1)  =  ID_TEMPER   
  visco_ker % llaws(15) % where     = 'IELEM'
  visco_ker % llaws(15) % kfl_gradi = 0

  visco_ker % llaws(16) % wname     = 'VISNP'       ! Viscosity at nodes
  visco_ker % llaws(16) % lresp(1)  =  ID_CHEMIC
  visco_ker % llaws(16) % where     = 'IPOIN'
  visco_ker % llaws(16) % kfl_gradi = 0
   
  !----------------------------------------------------------------------
  !
  ! Conductivity
  !
  !----------------------------------------------------------------------

  condu_ker % llaws(2) % wname    = 'KMIXT'        ! Mixture of species
  condu_ker % llaws(2) % lresp(1) =  ID_CHEMIC 
  condu_ker % llaws(2) % where     = 'IPOIN'
  condu_ker % llaws(2) % kfl_gradi = 1

  condu_ker % llaws(3) % wname    = 'SUTHE'        ! Sutherland
  condu_ker % llaws(3) % lresp(1) =  ID_TEMPER
  condu_ker % llaws(3) % where     = 'IPOIN'
  condu_ker % llaws(3) % kfl_gradi = 1

  condu_ker % llaws(4) % wname     = 'GPSUT'       ! Sutherland in Gauss points
  condu_ker % llaws(4) % lresp(1)  =  ID_NASTAL
  condu_ker % llaws(4) % lresp(2)  =  ID_TEMPER
  condu_ker % llaws(4) % where     = 'IELEM'
  condu_ker % llaws(4) % kfl_gradi =  1

  condu_ker % llaws(5) % wname    = 'KTABL'        ! from table (Flamelet combustion model)
  condu_ker % llaws(5) % lresp(1) =  ID_CHEMIC
  condu_ker % llaws(5) % where     = 'IPOIN'
  condu_ker % llaws(5) % kfl_gradi = 1

  condu_ker % llaws(6) % wname    = 'TKTAB'        ! from table (Flamelet combustion model)
  condu_ker % llaws(6) % lresp(1) =  ID_TEMPER     ! update only in temper
  condu_ker % llaws(6) % where     = 'IPOIN'
  condu_ker % llaws(6) % kfl_gradi = 1

  condu_ker % llaws(7) % wname    = 'KQUAR'        ! quartz glass, from polynomial function, for CHT
  condu_ker % llaws(7) % lresp(1) =  ID_TEMPER     ! update only in temper
  condu_ker % llaws(7) % where     = 'IPOIN'
  condu_ker % llaws(7) % kfl_gradi = 1

  condu_ker % llaws(8) % wname    = 'KSTAI'        ! stainless steel, from polynomial function, for CHT
  condu_ker % llaws(8) % lresp(1) =  ID_TEMPER     ! update only in temper
  condu_ker % llaws(8) % where     = 'IPOIN'
  condu_ker % llaws(8) % kfl_gradi = 1

  call ker_polynomial_allaws(condu_ker, 9_ip, ID_TEMPER, 'IELEM', 0_ip)

  condu_ker % llaws(10) % wname    = 'KGPTA'        ! from table on Gauss points (Flamelet combustion model)
  condu_ker % llaws(10) % lresp(1) =  ID_CHEMIC
  condu_ker % llaws(10) % where     = 'IELEM'
  condu_ker % llaws(10) % kfl_gradi = 0

  condu_ker % llaws(11) % wname    = 'TKNPO'        ! Thermal conductivity (unconditional if CMC) at nodes
  condu_ker % llaws(11) % lresp(1) =  ID_CHEMIC
  condu_ker % llaws(11) % where     = 'IPOIN'
  condu_ker % llaws(11) % kfl_gradi = 0
  
  !----------------------------------------------------------------------
  !
  ! Specific Heat
  !
  !----------------------------------------------------------------------

  sphea_ker % llaws(2) % wname    = 'CPMIX'        ! Mixture of species
  sphea_ker % llaws(2) % lresp(1) =  ID_CHEMIC 
  sphea_ker % llaws(2) % where    = 'IPOIN'

  sphea_ker % llaws(3) % wname    = 'CPTAB'        ! from table (Flamelet combustion model)
  sphea_ker % llaws(3) % lresp(2) =  ID_TEMPER
  sphea_ker % llaws(3) % where    = 'IPOIN'

  sphea_ker % llaws(4) % wname    = 'TCPTA'        ! from table (Flamelet combustion model)
  sphea_ker % llaws(4) % lresp(1) =  ID_TEMPER     ! update only in temper
  sphea_ker % llaws(4) % where    = 'IPOIN'

  sphea_ker % llaws(5) % wname    = 'CPQUA'        ! quartz glass, from polynomial function, for CHT
  sphea_ker % llaws(5) % lresp(1) =  ID_TEMPER     ! update only in temper
  sphea_ker % llaws(5) % where    = 'IPOIN'

  sphea_ker % llaws(6) % wname    = 'CPSTA'        ! stainless steel, from polynomial function, for CHT 
  sphea_ker % llaws(6) % lresp(1) =  ID_TEMPER     ! update only in temper
  sphea_ker % llaws(6) % where    = 'IPOIN'

  call ker_polynomial_allaws(sphea_ker, 7_ip, ID_TEMPER, 'IELEM', 0_ip) ! 

  sphea_ker % llaws(8) % wname    = 'CPGPT'        ! from table (Flamelet combustion model)
  sphea_ker % llaws(8) % lresp(1) =  ID_TEMPER     ! update in temper
  sphea_ker % llaws(8) % lresp(2) =  ID_CHEMIC     ! update also in chemic
  sphea_ker % llaws(8) % where    = 'IELEM'

  sphea_ker % llaws(9) % wname    = 'CPNPO'        ! Specific heat (unconditional if CMC) at nodes
  sphea_ker % llaws(9) % lresp(1) =  ID_CHEMIC
  sphea_ker % llaws(9) % where    = 'IPOIN'

  ! CMC is run without temper module --> CPGPT cannot be used
  sphea_ker % llaws(10) % wname    = 'CPGPO'        ! Specific heat (unconditional if CMC) at Gauss points
  sphea_ker % llaws(10) % lresp(1) =  ID_CHEMIC
  sphea_ker % llaws(10) % where    = 'IELEM'

  !----------------------------------------------------------------------
  !
  ! Dummy variable can be used for whatever
  !
  !----------------------------------------------------------------------

  dummy_ker % llaws(2) % wname     = 'BIFL2'       ! Bifluid
  dummy_ker % llaws(2) % lresp(1)  =  ID_LEVELS
  dummy_ker % llaws(2) % where     = 'IELEM'
  dummy_ker % llaws(2) % kfl_gradi = 1  ! BIFL2 includes gradient

  !----------------------------------------------------------------------
  !
  ! Turbulent viscosity
  !
  !----------------------------------------------------------------------

  turmu_ker % llaws( 2) % wname     = 'SMAGO'       ! Turbulent viscosity mu_t from Smagorinsky
  turmu_ker % llaws( 2) % lresp(1)  =  ID_NASTIN
  turmu_ker % llaws( 2) % where     = 'IELEM'
  turmu_ker % llaws( 2) % kfl_gradi = 1

  turmu_ker % llaws( 3) % wname     = 'WALE '       ! Turbulent viscosity mu_t from WALE
  turmu_ker % llaws( 3) % lresp(1)  =  ID_NASTIN
  turmu_ker % llaws( 3) % where     = 'IELEM'
  turmu_ker % llaws( 3) % kfl_gradi = 1 

  turmu_ker % llaws( 4) % wname     = 'SSTKO'       ! Turbulent viscosity mu_t from RANS K-W SST Model
  turmu_ker % llaws( 4) % lresp(1)  =  ID_TURBUL
  turmu_ker % llaws( 4) % where     = 'IELEM'
  turmu_ker % llaws( 4) % kfl_gradi = 1
  turmu_ker % llaws( 4) % kfl_deriv_tur = 1
  turmu_ker % llaws( 4) % kfl_deriv_vel = 1

  turmu_ker % llaws( 5) % wname     = 'STDKE'       ! Turbulent viscosity mu_t from RANS K-EPS Model
  turmu_ker % llaws( 5) % lresp(1)  =  ID_TURBUL
  turmu_ker % llaws( 5) % lresp(2)  =  ID_NASTIN    ! for ke-fp and realizable models
  turmu_ker % llaws( 5) % where     = 'IELEM'
  turmu_ker % llaws( 5) % kfl_gradi = 1

  turmu_ker % llaws( 6) % wname     = 'SPALA'       ! Turbulent viscosity mu_t from RANS Spalart-Almaras Model
  turmu_ker % llaws( 6) % lresp(1)  =  ID_TURBUL
  turmu_ker % llaws( 6) % where     = 'IELEM'
  turmu_ker % llaws( 6) % kfl_deriv_tur = 1

  turmu_ker % llaws( 7) % wname     = 'VRMAN'       ! Turbulent viscosity mu_t from Vreman
  turmu_ker % llaws( 7) % lresp(1)  =  ID_NASTIN
  turmu_ker % llaws( 7) % where     = 'IELEM'
  turmu_ker % llaws( 7) % kfl_gradi = 0

  turmu_ker % llaws( 8) % wname     = 'ILSA'       ! Turbulent viscosity mu_t from ILSA
  turmu_ker % llaws( 8) % lresp(1)  =  ID_NASTIN
  turmu_ker % llaws( 8) % where     = 'IELEM'
  turmu_ker % llaws( 8) % kfl_gradi = 1

  turmu_ker % llaws( 9) % wname     = 'MIXIN'       ! Turbulent viscosity mu_t from algebraic mixing-length model
  turmu_ker % llaws( 9) % lresp(1)  =  ID_NASTIN
  turmu_ker % llaws( 9) % where     = 'IELEM'
  turmu_ker % llaws( 9) % kfl_gradi = 1

  turmu_ker % llaws(10) % wname     = 'SIGMA'       ! Turbulent viscosity mu_t from Vreman
  turmu_ker % llaws(10) % lresp(1)  =  ID_NASTIN
  turmu_ker % llaws(10) % where     = 'IELEM'
  turmu_ker % llaws(10) % kfl_gradi = 1

  turmu_ker % llaws(11) % wname     = 'KOMEG'       ! Turbulent viscosity mu_t from RANS K-OMEGA model
  turmu_ker % llaws(11) % lresp(1)  =  ID_TURBUL
  turmu_ker % llaws(11) % where     = 'IELEM'
  turmu_ker % llaws(11) % kfl_gradi = 1

  turmu_ker % llaws(12) % wname     = 'TEST '       ! Turbulent viscosity mu_t from RANS K-OMEGA model
  turmu_ker % llaws(12) % lresp(1)  =  ID_NASTIN
  turmu_ker % llaws(12) % where     = 'IELEM'
  turmu_ker % llaws(12) % kfl_gradi = 0

  turmu_ker % llaws(13) % wname     = 'TKESG'       ! Turbulent viscosity mu_t from a TKE-SGS transport equation
  turmu_ker % llaws(13) % lresp(1)  =  ID_TURBUL
  turmu_ker % llaws(13) % lresp(1)  =  ID_NASTIN
  turmu_ker % llaws(13) % where     = 'IELEM'
  turmu_ker % llaws(13) % kfl_gradi = 0 
  
  turmu_ker % llaws(14) % wname     = 'FIELD'       ! Turbulent viscosity mu_t given by a field
  turmu_ker % llaws(14) % lresp(1)  =  ID_KERMOD
  turmu_ker % llaws(14) % where     = 'IPOIN'
  turmu_ker % llaws(14) % kfl_gradi = 0 
  
  !----------------------------------------------------------------------
  !
  ! Porosity
  !
  !----------------------------------------------------------------------

  poros_ker % llaws(2) % wname     = 'CANOP'     ! Canopy Porosity
  poros_ker % llaws(2) % lresp(1)  =  ID_NASTIN
  poros_ker % llaws(2) % where     = 'IELEM'
  poros_ker % llaws(2) % kfl_gradi = 0

  poros_ker % llaws(3) % wname     = 'VALVE'     ! Canopy Porosity
  poros_ker % llaws(3) % lresp(1)  =  ID_NASTIN
  poros_ker % llaws(3) % where     = 'IELEM'
  poros_ker % llaws(3) % kfl_gradi = 0

  poros_ker % llaws(4) % wname     = 'DAFOR' ! Darcy-Forchheimer porosity
  poros_ker % llaws(4) % lresp(1)  =  ID_NASTIN
  poros_ker % llaws(4) % where     = 'IELEM'
  poros_ker % llaws(4) % kfl_gradi = 0


  poros_ker % llaws(5) % wname     = 'FUNCT'     ! Time function porosity
  poros_ker % llaws(5) % lresp(1)  =  ID_NASTIN
  poros_ker % llaws(5) % where     = 'IELEM'
  poros_ker % llaws(5) % kfl_gradi = 0

  poros_ker % llaws(6) % wname     = 'DAMPI'     !Damping term for ABL
  poros_ker % llaws(6) % lresp(1)  =  ID_NASTIN
  poros_ker % llaws(6) % where     = 'IELEM'
  poros_ker % llaws(6) % kfl_gradi = 0


  !----------------------------------------------------------------------
  !
  ! Mixing function
  !
  !----------------------------------------------------------------------

  mixin_ker % llaws(2) % wname     = 'IMMER'       ! Bifluid
  mixin_ker % llaws(2) % lresp(1)  =  ID_LEVELS
  mixin_ker % llaws(2) % where     = 'IELEM'
  mixin_ker % llaws(2) % kfl_gradi = 0             ! Despite there is a gradient we never use it for Bifulid flow - see if in nsi_elmop3 grvis

  !----------------------------------------------------------------------
  !
  ! Anisotropic porosity
  !
  !----------------------------------------------------------------------

  anipo_ker % llaws(2) % wname     = 'MATRI'       ! Matrix porosity? really used?
  anipo_ker % llaws(2) % lresp(1)  =  ID_NASTIN  
  anipo_ker % llaws(2) % where     = 'IELEM'
  anipo_ker % llaws(2) % kfl_gradi = 0            

  anipo_ker % llaws(3) % wname     = 'DAMPI'       ! anisotropic damping term for ABL
  anipo_ker % llaws(3) % lresp(1)  =  ID_NASTIN
  anipo_ker % llaws(3) % where     = 'IELEM'
  anipo_ker % llaws(3) % kfl_gradi = 0

  anipo_ker % llaws(4) % wname     = 'FIELD'       ! Given by a field
  anipo_ker % llaws(4) % lresp(1)  =  ID_KERMOD  
  anipo_ker % llaws(4) % where     = 'IELEM'
  anipo_ker % llaws(4) % kfl_gradi = 0         
  
  !----------------------------------------------------------------------
  !
  ! Wall viscosity
  !
  !----------------------------------------------------------------------

  walvi_ker % llaws(2) % wname     = 'WVSTD' 
  walvi_ker % llaws(2) % lresp(1)  =  ID_NASTIN  
  walvi_ker % llaws(2) % where     = 'IELEM'
  walvi_ker % llaws(2) % kfl_gradi = 0      

end subroutine ker_allaws
