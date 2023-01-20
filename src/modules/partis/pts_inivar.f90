!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine pts_inivar(itask)
  !-----------------------------------------------------------------------
  !****f* Partis/pts_inivar
  ! NAME 
  !    pts_inivar
  ! DESCRIPTION
  !    Initial variables
  ! USES
  ! USED BY
  !    pts_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_partis
  use def_solver 
  use mod_memory,              only : memory_alloca
  use mod_memory,              only : memory_deallo
  use mod_communications,      only : PAR_SUM
  use mod_pts_injection,       only : pts_injection_initialization
  use mod_pts_parallelization, only : pts_parallelization_initialization
  use mod_physics,             only : physics_initialize_liquid
  use mod_physics,             only : physics_set_liquid_temperature
  use mod_physics,             only : physics_mole_2_mass
  use mod_physics,             only : universal_gas_constant 
  use mod_interp_tab,          only : fw_lookup 
  use mod_arrays,              only : arrays_register
  use mod_maths,               only : maths_local_orthonormal_basis
      
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: itype,ii,ivari,ivar2,iinj
  real(rp), pointer       :: retva(:)
  real(rp), pointer       :: control(:)
  real(rp), pointer       :: tab_conce(:)
  integer(ip), pointer    :: ind(:) 
  real(rp)                :: T_print

  select case( itask )

  case( 0_ip )
     !
     ! Postprocess
     !
     call arrays_register((/'MOMSK','VECTO','NPOIN','PRIMA'/),momentum_sink,ENTITY_POSITION=2_ip)
     call arrays_register((/'HEASK','SCALA','NPOIN','PRIMA'/),heat_sink    ,ENTITY_POSITION=1_ip)
     call arrays_register((/'MASSK','SCALA','NPOIN','PRIMA'/),mass_sink    ,ENTITY_POSITION=1_ip)
     
     call arrays_register((/'DEPOE','SCALA','NELEM','PRIMA'/),depoe_pts,    ENTITY_POSITION=2_ip)
     call arrays_register((/'DEPOB','SCALA','NBOUN','PRIMA'/),depob_pts,    ENTITY_POSITION=2_ip)     
     call arrays_register((/'RESID','SCALA','NELEM','PRIMA'/),resid_pts,    ENTITY_POSITION=2_ip)

     call arrays_register((/'NRESI','SCALA','NPOIN','SECON'/),              ENTITY_POSITION=2_ip,COMPONENT_POSITION=1_ip)
     call arrays_register((/'PARTI','SCALA','NPOIN','SECON'/),              ENTITY_POSITION=2_ip,COMPONENT_POSITION=1_ip)

     call arrays_register((/'AVMIE','SCALA','NPOIN','PRIMA'/),avg_mie_pts  ,ENTITY_POSITION=1_ip)

     call arrays_register((/'DEPOS','MULTI','NPOIN','SECON'/))
     call arrays_register((/'SLIPW','SCALA','NPOIN','SECON'/))
     call arrays_register((/'FRICT','SCALA','NPOIN','SECON'/))
     call arrays_register((/'DEPOT','SCALA','NPOIN','SECON'/))
     call arrays_register((/'BOUNC','SCALA','NPOIN','SECON'/))
     call arrays_register((/'TEMPE','SCALA','NPOIN','SECON'/))

     !
     ! Witness geometry variables
     !
     postp(1) % wowig  ( 1)  = 'NUMBE' ! Number of particles     
     postp(1) % wowig  ( 2)  = 'DIA10' ! Diameter 
     postp(1) % wowig  ( 3)  = 'DIA32' ! SMD
     postp(1) % wowig  ( 4)  = 'TEMPE' ! Temperature
     postp(1) % wowig  (5:7) = 'VELOC' ! Velocity
     postp(1) % wowig  ( 8)  = 'DIA20' ! Diameter square
     postp(1) % wowig  ( 9)  = 'TEMP2' ! Temperature square
     postp(1) % wowig (10:12)= 'VELO2' ! Velocity square
     postp(1) % wowig  (13)  = 'MASS ' ! Mass
     postp(1) % wowig (14:16)= 'MFLUX' ! Mass flux

     !
     ! Others
     !
     nlagr                 = 0       ! Absolute number of lagrangian particles
     nlacc_pts             = 0       ! Number of accumulated particles
     kfl_depos_pts         = 0       ! Deposition map on mesh is not postprocess
     kfl_timei             = 1       ! Problem is always transient
     particles_sent        = 0       ! # of particles received
     particles_recv        = 0       ! # of particles sent
     comm_loops_pts        = 0       ! # communication loops 
     kfl_resid_pts         = 0       ! Residence time not required
     kfl_injec             = 0       ! No injection
     nlagr_local_pts       = 0       ! Position taken in type
     nlagr_free_pts        = 0       ! Free positions taken in type
     kfl_slip_wall_pts     = 0       ! Slip wall
     kfl_bouncing_wall_pts = 0       ! Bouncing wall
     nvarp_pts             = 0       ! Number of deposition variables
     nvard_pts             = 0       ! Number of deposition variables
 
     nullify(permu_nlagr_pts)        ! Permutation
     nullify(depoe_pts)              ! Deposition map over elements
     nullify(depob_pts)              ! Deposition map over boundaries
     nullify(defor_pts)              ! Deformation tensor
     nullify(hleng_pts)              ! Element characteristics length
     nullify(lboue_pts)
     nullify(kfl_fixbo_pts)     
     nullify(bvnat_pts)  
     nullify(tbcod_pts)
     nullify(leleboun_pts) 
     nullify(bouno_pts)  
     nullify(walld_slip_pts)
     nullify(kfl_fixno_walld_slip_pts)
     nullify(friction_pts)
     nullify(walld_bouncing_pts)
     nullify(kfl_fixno_walld_bouncing_pts)
     nullify(resid_pts)
     nullify(postprocess_list_pts)
     nullify(deposition_list_pts)
     nullify(avg_mie_pts)
     !
     ! Postprocess variables
     !     
     postprocess_name_pts     = 'NULL'
     
     postprocess_name_pts( 1) = 'T'       ! t              
     postprocess_name_pts( 2) = 'ILAGR'   ! ilagr          
     postprocess_name_pts( 3) = 'ITYPE'   ! itype          
     postprocess_name_pts( 4) = 'IELEM'   ! ielem          
     postprocess_name_pts( 5) = 'EXIST'   ! kfl_exist      
     postprocess_name_pts( 6) = 'ITTIM'   ! ittim          
     postprocess_name_pts( 7) = 'SET  '   ! boundary_set     
     postprocess_name_pts( 8) = 'COORX'   ! coord(1)       
     postprocess_name_pts( 9) = 'COORY'   ! coord(2)       
     postprocess_name_pts(10) = 'COORZ'   ! coord(3)       
     postprocess_name_pts(11) = 'VELOX'   ! veloc(1)       
     postprocess_name_pts(12) = 'VELOY'   ! veloc(2)       
     postprocess_name_pts(13) = 'VELOZ'   ! veloc(3)       
     postprocess_name_pts(14) = 'ACCEX'   ! accel(1)       
     postprocess_name_pts(15) = 'ACCEY'   ! accel(2)       
     postprocess_name_pts(16) = 'ACCEZ'   ! accel(3)       
     postprocess_name_pts(17) = 'DTK  '   ! dt_k           
     postprocess_name_pts(18) = 'CD   '   ! Cd             
     postprocess_name_pts(19) = 'STK1 '   ! Stk(1)         
     postprocess_name_pts(20) = 'STK2 '   ! Stk(2)         
     postprocess_name_pts(21) = 'VF1  '   ! v_fluid_k(1)   
     postprocess_name_pts(22) = 'VF1  '   ! v_fluid_k(2)   
     postprocess_name_pts(23) = 'VF1  '   ! v_fluid_k(3)   
     postprocess_name_pts(24) = 'ACCDX'   ! acced(1)       
     postprocess_name_pts(25) = 'ACCDY'   ! acced(2)       
     postprocess_name_pts(26) = 'ACCDZ'   ! acced(3)       
     postprocess_name_pts(27) = 'ACCFX'   ! accee(1)       
     postprocess_name_pts(28) = 'ACCFY'   ! accee(2)       
     postprocess_name_pts(29) = 'ACCFZ'   ! accee(3)       
     postprocess_name_pts(30) = 'ACCGX'   ! acceg(1)       
     postprocess_name_pts(31) = 'ACCGY'   ! acceg(2)       
     postprocess_name_pts(32) = 'ACCGZ'   ! acceg(3)       
     postprocess_name_pts(33) = 'STRET'   ! stret          
     postprocess_name_pts(34) = 'DTKM1'   ! dt_km1         
     postprocess_name_pts(35) = 'DTKM2'   ! dt_km2         
     postprocess_name_pts(36) = 'DTG  '   ! dtg            
     postprocess_name_pts(37) = 'VFM1X'   ! v_fluid_km1(1) 
     postprocess_name_pts(38) = 'VFM1Y'   ! v_fluid_km1(2) 
     postprocess_name_pts(39) = 'VFM1Z'   ! v_fluid_km1(3) 
     postprocess_name_pts(40) = 'VFM2X'   ! v_fluid_km2(1) 
     postprocess_name_pts(41) = 'VFM2Y'   ! v_fluid_km2(2) 
     postprocess_name_pts(42) = 'VFM2Z'   ! v_fluid_km2(3) 
     postprocess_name_pts(43) = 'TINJE'   ! t_inject       
     postprocess_name_pts(44) = 'COKX '   ! coord_k(1)     
     postprocess_name_pts(45) = 'COKY '   ! coord_k(2)     
     postprocess_name_pts(46) = 'COKZ '   ! coord_k(3)     
     postprocess_name_pts(47) = 'COK1X'   ! coord_km1(1)   
     postprocess_name_pts(48) = 'COK1Y'   ! coord_km1(2)   
     postprocess_name_pts(49) = 'COK1Z'   ! coord_km1(3)     
     postprocess_name_pts(50) = 'DISTA'   ! dista          
     postprocess_name_pts(51) = 'COO1D'   ! coord1d        
     postprocess_name_pts(52) = 'SIGN '   ! sign           
     postprocess_name_pts(53) = 'TEMPK'   ! tempe_k           
     postprocess_name_pts(54) = 'TEKM1'   ! tempe_km1     
     postprocess_name_pts(55) = 'MASSK'   ! mass_k      
     postprocess_name_pts(56) = 'MAKM1'   ! mass_km1
     postprocess_name_pts(57) = 'MPIRA'   ! mpi_rank
     postprocess_name_pts(58) = 'DIAMK'   ! diam_k      
     postprocess_name_pts(59) = 'TEMPF'   ! Temp_fluid_k
     postprocess_name_pts(60) = 'YVAPF'   ! Yvap_fluid_k
     postprocess_name_pts(61) = 'DIAM0'   ! diam_0
     postprocess_name_pts(62) = 'MASS0'   ! mass_0
     postprocess_name_pts(63) = 'BM   '   ! BM
     !                    
     ! Solver
     !
     call soldef(-2_ip)                     ! Allocate memory
     solve(1) % kfl_solve = 1               ! Slip wall solver
     solve(1) % wprob     = 'SLIP_WALL'     ! Name
     solve(1) % kfl_iffix = 2               ! Fix Dirichlet with value given by initial solution

     solve(2) % kfl_solve = 1               ! Bouncing wall solver
     solve(2) % wprob     = 'BOUNCING_WALL' ! Name
     solve(2) % kfl_iffix = 2               ! Fix Dirichlet with value given by initial solution
     !
     ! Others
     ! 
     nlagr_existing_pts      = 0
     nlagr_non_migrating_pts = 0
     nlagr_going_out_pts     = 0
     nlagr_zero_time_pts     = 0
     nlagr_deposited_pts     = 0
     nlagr_hits_wall_pts     = 0
     
     !
     ! Injection
     !  
     call pts_injection_initialization()
   
  case( 1_ip )
     !
     ! After reading the data
     !
     !
     ! Number of used types
     !
     ntyla_pts = 0
     number_types_pts = 0
     do itype = 1,mtyla
        if( parttyp(itype) % kfl_exist == 1 ) then
           ntyla_pts = max(ntyla_pts,itype)
           number_types_pts = number_types_pts + 1           
        end if
     end do

     !
     ! Initialize for droplet size associated with injector
     !
     do iinj = 1, kfl_imax_injector_pts
        itype = 1_ip
        if (injection_pts(iinj) % kfl_particle_type > 0 ) then 
           itype = injection_pts(iinj) % kfl_particle_type
        endif
        if (parttyp(itype) % kfl_exist == 1) then
           if ( injection_pts(iinj) % size_diame == 0.0_rp) then
              injection_pts(iinj) % size_diame = parttyp(itype) % diame 
           endif
        endif
     enddo

     !
     ! Thermodynamic pressure
     !
     if (prthe_pts /= 0) then
        prthe(1) = prthe_pts
        prthe(2) = prthe_pts
        prthe(3) = prthe_pts
        prthe(4) = prthe_pts
     endif


     !
     ! Initialize liquid structure
     !
     do itype = 1,mtyla
         if (( parttyp(itype) % kfl_exist == 1 ) .and. &
             ( parttyp(itype) % kfl_therm /= 0 )) then
             
             if (parttyp(itype) % liq % name == 'USER ') then
                 call physics_initialize_liquid(parttyp(itype) % liq, 298.15_rp, prthe(1),   & 
                                     parttyp(itype) % L_vapor, parttyp(itype) % T_boiling,   &
                                     parttyp(itype) % Cp,      parttyp(itype) % denpa,       & 
                                     parttyp(itype) % w,       parttyp(itype) % T_crit,      &
                                     parttyp(itype) % P_boiling)
             else
                 !
                 ! Using hardcoded liquids
                 !
                 call physics_initialize_liquid(parttyp(itype) % liq, 298.15_rp, prthe(1))
             endif 
             


             if ( parttyp(itype) % kfl_therm == 2 ) then
                 !
                 ! Check if old table is used:
                 !
                 if (parttyp(itype) % kfl_spr_fw <=0) then
                    parttyp(itype) % kfl_spr_fw =  parttyp(itype) % kfl_tab_fw 
                    parttyp(itype) % spr_tab_fw => parttyp(itype) % table_fw
                 endif


                 !
                 ! DETERMINE TABLE INDECES
                 !
                 if (parttyp(itype) % kfl_tab_fw > 0)  then
                    parttyp(itype) % kfl_W_tab_index             = 0_ip 
                    parttyp(itype) % kfl_k_tab_index             = 0_ip 
                    parttyp(itype) % kfl_mu_tab_index            = 0_ip 
                    parttyp(itype) % kfl_cpCoefLT_tab_index      = 0_ip 
                    parttyp(itype) % kfl_cpCoefHT_tab_index      = 0_ip 
                    parttyp(itype) % kfl_cpCoefLT_end_tab_index  = 0_ip 
                    parttyp(itype) % kfl_cpCoefHT_end_tab_index  = 0_ip 
                    do ivari=1,parttyp(itype) % table_fw % main_table % nvar
                       select case(parttyp(itype) % table_fw % main_table % varname(ivari)) 
                       case('W','WSTAR')
                           parttyp(itype) % kfl_W_tab_index       = ivari
                       case('K')
                           parttyp(itype) % kfl_k_tab_index       = ivari
                       case('MU')
                           parttyp(itype) % kfl_mu_tab_index      = ivari
                       case('CP11')
                           parttyp(itype) % kfl_cpCoefLT_tab_index     = ivari
                           parttyp(itype) % kfl_cpCoefLT_end_tab_index = ivari+5
                       case('CP21')
                           parttyp(itype) % kfl_cpCoefHT_tab_index     = ivari
                           parttyp(itype) % kfl_cpCoefHT_end_tab_index = ivari+5
                       end select
                    enddo
                 endif

                 if (parttyp(itype) % kfl_spr_fw > 0) then
                    parttyp(itype) % kfl_Yfuel_spr_index         = 0_ip 
                    parttyp(itype) % kfl_Dfuel_spr_index         = 0_ip 
                    parttyp(itype) % kfl_dkdT_spr_index          = 0_ip 
                    parttyp(itype) % kfl_dmudT_spr_index         = 0_ip 
                    parttyp(itype) % kfl_dDfueldT_spr_index      = 0_ip 
                    parttyp(itype) % kfl_T_spr_index             = 0_ip 
                    do ivari=1,parttyp(itype) % spr_tab_fw % main_table % nvar
                       select case(parttyp(itype) % spr_tab_fw % main_table % varname(ivari)) 
                       case('YFUEL')
                           parttyp(itype) % kfl_Yfuel_spr_index     = ivari
                       case('DFUEL')
                           parttyp(itype) % kfl_Dfuel_spr_index     = ivari
                       case('DKDT')
                           parttyp(itype) % kfl_dkdT_spr_index      = ivari
                       case('DMUDT')
                           parttyp(itype) % kfl_dmudT_spr_index     = ivari
                       case('DDFUE')
                           parttyp(itype) % kfl_dDfueldT_spr_index  = ivari
                       case('T','TSTAR')
                           parttyp(itype) % kfl_T_spr_index         = ivari
                       end select
                    enddo
                 endif


                 !
                 ! INITIALIZE CP COEFFICENTS OF FUEL
                 !
                 nullify(control) 
                 nullify(retva) 
                 nullify(tab_conce) 
                 nullify(ind) 

                 call memory_alloca(mem_modul(1:2,modul),'CONTROL'  ,'pts_inivar',control  ,parttyp(itype) % table_fw % main_table % ndim)
                 call memory_alloca(mem_modul(1:2,modul),'RETVA'    ,'pts_inivar',retva    ,parttyp(itype) % table_fw % main_table % nvar)
                 call memory_alloca(mem_modul(1:2,modul),'TAB_CONCE','pts_inivar',tab_conce,parttyp(itype) % table_fw % main_table % ndim)
                 call memory_alloca(mem_modul(1:2,modul),'IND'      ,'pts_inivar',ind      ,parttyp(itype) % table_fw % main_table % ndim)

                 !
                 ! Lookup seen state ASSUMES ADAIBATIC FOR NOW
                 !        
                 control = 0.0_rp
                 ind     = 1_ip
                 do ii = 1, parttyp(itype) % table_fw % main_table % ndim
                    select case (parttyp(itype) % table_fw % main_table % coords(ii) % name)
                    case ('CMEAN','C    ')
                        control(ii) = 0.0_rp
                    case ('CVAR ')
                        control(ii) = 0.0_rp
                    case ('CHIST')
                        control(ii) = 0.0_rp 
                    case ('ZMEAN','Z    ')
                        control(ii) = 1.0_rp
                    case ('ZVAR ')
                        control(ii) = 0.0_rp
                    case ('IMEAN','I    ')
                        !
                        ! At this point the enthalpy shouldn't matter, 
                        ! because the Cp coefficeints of the fule should not change
                        ! with enthalpy
                        !
                        control(ii) = 0.0_rp
                    end select
                 enddo
                 call fw_lookup( control, tab_conce, parttyp(itype) % table_fw, retva, ind )
                 
                 !
                 ! First is low temperature, second is high temperature in table
                 !
                 parttyp(itype) % cpcoef_v_chm(1:6,1) = retva(parttyp(itype)%kfl_cpCoefLT_tab_index:parttyp(itype)%kfl_cpCoefLT_end_tab_index)
                 parttyp(itype) % cpcoef_v_chm(1:6,2) = retva(parttyp(itype)%kfl_cpCoefHT_tab_index:parttyp(itype)%kfl_cpCoefHT_end_tab_index)
                 
                 call memory_deallo(mem_modul(1:2,modul),'CONTROL'  ,'pts_inivar',control  )
                 call memory_deallo(mem_modul(1:2,modul),'RETVA'    ,'pts_inivar',retva    )
                 call memory_deallo(mem_modul(1:2,modul),'TAB_CONCE','pts_inivar',tab_conce)
                 call memory_deallo(mem_modul(1:2,modul),'IND'      ,'pts_inivar',ind      )
              endif




              if (INOTSLAVE) then
                !
                ! Print data of table lookup
                !
                if (parttyp(itype) % kfl_tab_fw > 0)  then
                    !
                    ! Common tabulated terms
                    !
                    write(momod(modul) % lun_outpu,'(A)')'---------------------------------------------'
                    write(momod(modul) % lun_outpu,'(A,I3)')        'GAS PROPERTIES OF TYPE:',itype
                    write(momod(modul) % lun_outpu,'(2X,A,1X,I4)')  'COLUMN INDECES OF GAS PROPERTIES FROM FRAMEWORK:', parttyp(itype) % kfl_tab_fw
                    write(momod(modul) % lun_outpu,'(5X,A9,1X,I4)')       'W        ', parttyp(itype) % kfl_W_tab_index         
                    write(momod(modul) % lun_outpu,'(5X,A9,1X,I4)')       'K        ', parttyp(itype) % kfl_k_tab_index         
                    write(momod(modul) % lun_outpu,'(5X,A9,1X,I4)')       'MU       ', parttyp(itype) % kfl_mu_tab_index         
                    write(momod(modul) % lun_outpu,'(5X,A9,1X,I4,1X,I4)') 'CPijL    ', parttyp(itype) % kfl_cpCoefLT_tab_index, parttyp(itype) % kfl_cpCoefLT_end_tab_index    
                    write(momod(modul) % lun_outpu,'(5X,A9,1X,I4,1X,I4)') 'CPijH    ', parttyp(itype) % kfl_cpCoefHT_tab_index, parttyp(itype) % kfl_cpCoefHT_end_tab_index    
                endif

                if (parttyp(itype) % kfl_spr_fw > 0)  then
                    !
                    ! Spray specific tabulated terms
                    !
                    write(momod(modul) % lun_outpu,'(2X,A,1X,I4)')  'COLUMN INDECES OF SPRAY SPECIFIC PROPERTIES FROM FRAMEWORK:', parttyp(itype) % kfl_spr_fw
                    write(momod(modul) % lun_outpu,'(5X,A9,1X,I4)')       'Y_VAP    ', parttyp(itype) % kfl_Yfuel_spr_index         
                    write(momod(modul) % lun_outpu,'(5X,A9,1X,I4)')       'D_VAP    ', parttyp(itype) % kfl_Dfuel_spr_index         
                    write(momod(modul) % lun_outpu,'(5X,A9,1X,I4)')       'dK/dT    ', parttyp(itype) % kfl_dkdT_spr_index       
                    write(momod(modul) % lun_outpu,'(5X,A9,1X,I4)')       'dMU/dT   ', parttyp(itype) % kfl_dmudT_spr_index         
                    write(momod(modul) % lun_outpu,'(5X,A9,1X,I4)')       'dD_VAP/dT', parttyp(itype) % kfl_dDfueldT_spr_index         
                    write(momod(modul) % lun_outpu,'(5X,A9,1X,I4)')       'T        ', parttyp(itype) % kfl_T_spr_index         
                endif


                !
                ! Print some data about the LIQUID that we use
                !
                write(momod(modul) % lun_outpu,'(A)')'---------------------------------------------'
                write(momod(modul) % lun_outpu,'(A,I3)')           'LIQUID PROPERTIES OF TYPE:',itype
                write(momod(modul) % lun_outpu,'(A,A)')            'NAME:                   ',parttyp(itype) % liq % name
                write(momod(modul) % lun_outpu,'(A,f16.8,1X,A)')   'CRITICAL PRESSURE:     ',parttyp(itype) % liq % Pc / 1.0e5_rp, 'bar'
                write(momod(modul) % lun_outpu,'(A,f16.8,1X,A)')   'CRITICAL TEMPERATURE:  ',parttyp(itype) % liq % Tc , 'K'
                write(momod(modul) % lun_outpu,'(A,f16.8,1X,A)')   'SATURATION TEMPERATURE:',parttyp(itype) % liq % Tsat , 'K'
                if ( parttyp(itype) % kfl_therm == 2 ) then
                    !
                    ! NASA Polynomials 
                    !
                    write(momod(modul) % lun_outpu,'(A,6(e16.8,1X))')   'LOW  TEMPERATURE CP:  ',parttyp(itype) % cpcoef_v_chm(:,1)
                    write(momod(modul) % lun_outpu,'(A,6(e16.8,1X))')   'HIGH TEMPERATURE CP:  ',parttyp(itype) % cpcoef_v_chm(:,2)

                endif
                
                !
                ! Output some properties to .pts.log file
                !
                T_print = 240.0_rp
                do while (T_print < parttyp(itype) % liq % Tsat) 
                   T_print = T_print + 10.0_rp
                   T_print = min(T_print,parttyp(itype) % liq % Tsat)
                   call physics_set_liquid_temperature(parttyp(itype) % liq, T_print)
                   write(momod(modul) % lun_outpu,'(A,f16.8,1X,A)')   '@',parttyp(itype) % liq % T , 'K:'
                   write(momod(modul) % lun_outpu,'(A,e16.8E3,1X,A)') '   RHO =               ',parttyp(itype) % liq % rho, 'kg/m^3'
                   write(momod(modul) % lun_outpu,'(A,e16.8E3,1X,A)') '   CP =                ',parttyp(itype) % liq % cp, 'J/kgK'
                   write(momod(modul) % lun_outpu,'(A,e16.8E3,1X,A)') '   P_SAT =             ',parttyp(itype) % liq % psat  / 1.0e5_rp, 'bar'
                   write(momod(modul) % lun_outpu,'(A,e16.8E3,1X,A)') '   X = P_SAT/P =       ',parttyp(itype) % liq % psat / parttyp(itype) % liq % P, ''
                   write(momod(modul) % lun_outpu,'(A,e16.8E3,1X,A)') '   Y in N2 bath gas =  ',physics_mole_2_mass(parttyp(itype) % liq % psat / parttyp(itype) % liq % P, parttyp(itype) % liq % W, 28.0e-3_rp), ''
                   write(momod(modul) % lun_outpu,'(A,e16.8E3,1X,A)') '   LV =                ',parttyp(itype) % liq % Lv, 'J/kg'
                   write(momod(modul) % lun_outpu,'(A)')'---------------------------------------------'
                enddo

                !
                ! Output PROPERTIES TO FILE
                !
                open(1766, file=adjustl(trim(namda))//'-'//trim(parttyp(itype) % liq % name)//'-prop.'//exmod(modul)//'.res',status='replace')

                !
                ! Add header
                !
                write(1766,'(7(A16,1x))')  '#          T [K]', 'rho [kg/m^3]', 'cp [J/kgK]', 'Psat [Pa]', 'Lv [J/kg]', 'Xsat', 'Ysat \w N2'

                !
                ! Subcooled conditions
                !
                T_print = 249.0_rp
                do while (T_print < parttyp(itype) % liq % Tsat) 
                   T_print = T_print + 1.0_rp
                   T_print = min(T_print,parttyp(itype) % liq % Tsat)
                   call physics_set_liquid_temperature(parttyp(itype) % liq, T_print)
                   write(1766,'(7(e16.8E3,1x))')          &
                      &   parttyp(itype) % liq % T,         &
                      &   parttyp(itype) % liq % rho,       &
                      &   parttyp(itype) % liq % cp,        &
                      &   parttyp(itype) % liq % psat,      &
                      &   parttyp(itype) % liq % Lv,        &
                      &   parttyp(itype) % liq % psat / parttyp(itype) % liq % P, &
                      &   physics_mole_2_mass(parttyp(itype) % liq % psat / parttyp(itype) % liq % P, parttyp(itype) % liq % W, 28.0e-3_rp)

                enddo

                !
                ! Go a bit beyond saturation temperature (superheated conditions):
                !
                do ii = 1,50
                   T_print = parttyp(itype) % liq % Tsat + real(ii,KIND=rp)
                   call physics_set_liquid_temperature(parttyp(itype) % liq, T_print)
                   write(1766,'(5(e16.8E3,1x),2(a,1x))')          &
                      &   parttyp(itype) % liq % T,         &
                      &   parttyp(itype) % liq % rho,       &
                      &   parttyp(itype) % liq % cp,        &
                      &   parttyp(itype) % liq % psat,      &
                      &   parttyp(itype) % liq % Lv,        &
                      &   'NaN', &
                      &   'NaN'
                     
                enddo
 
                close(1766)                    
              endif


         endif
     enddo
     !
     ! Postprocess
     !
     nvarp_pts = count(postprocess_var_pts)
     call memory_alloca(mem_modul(1:2,modul),'POSTPROCESS_LIST_PTS','pts_output',postprocess_list_pts,nvarp_pts)
     ivar2 = 0
     do ivari = 1,mvarp_pts
        if( postprocess_var_pts(ivari) ) then
           ivar2 = ivar2 + 1
           postprocess_list_pts(ivar2) = ivari
        end if
     end do
     nvard_pts = count(deposition_var_pts)
     call memory_alloca(mem_modul(1:2,modul),'DEPOSITION_LIST_PTS','pts_output',deposition_list_pts,nvard_pts)
     ivar2 = 0
     do ivari = 1,mvarp_pts
        if( deposition_var_pts(ivari) ) then
           ivar2 = ivar2 + 1
           deposition_list_pts(ivar2) = ivari
        end if
     end do

  case( 2_ip )
     !
     ! Injection through the boundary
     !
     if( kfl_boundary_injection == 1_ip) call pts_injbou()

     !
     ! Initialize the coordinate basis asssociated with injector
     !
     do iinj = 1, kfl_imax_injector_pts
        if ( dot_product(injection_pts(iinj) % geo_normal,injection_pts(iinj) % geo_normal) > 1.0e-16_rp ) then
           injection_pts(iinj) % geo_basis(1:ndime,1) = injection_pts(iinj) % geo_normal(1:ndime)
           call maths_local_orthonormal_basis(ndime, injection_pts(iinj) % geo_basis)
        endif
     enddo

     !
     ! Initialize time management of injector 
     !
     do iinj = 1, kfl_imax_injector_pts
        if (injection_pts(iinj) % time_initial < 0) injection_pts(iinj) % time_initial = tinla_pts
        if (injection_pts(iinj) % time_period  < 0) injection_pts(iinj) % time_period  = tpela_pts
        if (injection_pts(iinj) % time_final   < 0) injection_pts(iinj) % time_final   = tfila_pts
        injection_pts(iinj) % time_cumulative = 1e12_rp
     enddo


     !
     ! Injection and parallelization
     !
     call pts_parallelization_initialization()


     !
     ! Pre-compute used geometrical witness properties
     !
     kfl_max_prop_gwit_pts = 0_ip
     do ii = 1,nvarg
        if (postp(1) % npp_witng(ii) > 0) kfl_max_prop_gwit_pts = ii
     enddo  

   case ( 3_ip )
      !
     ! Default postprocess
     !
     postprocess_var_pts(pts_name_to_variable_number(    'T')) = .true.
     postprocess_var_pts(pts_name_to_variable_number('ILAGR')) = .true.
     postprocess_var_pts(pts_name_to_variable_number('ITYPE')) = .true.
     postprocess_var_pts(pts_name_to_variable_number('EXIST')) = .true.
     postprocess_var_pts(pts_name_to_variable_number('COORX')) = .true.
     postprocess_var_pts(pts_name_to_variable_number('COORY')) = .true.
     postprocess_var_pts(pts_name_to_variable_number('COORZ')) = .true.
     postprocess_var_pts(pts_name_to_variable_number('VELOX')) = .true.
     postprocess_var_pts(pts_name_to_variable_number('VELOY')) = .true.
     postprocess_var_pts(pts_name_to_variable_number('VELOZ')) = .true.
     postprocess_var_pts(pts_name_to_variable_number(  'DTK')) = .true.
     postprocess_var_pts(pts_name_to_variable_number(   'CD')) = .true.

     if( kfl_thermo_pts /= 0 ) then
        postprocess_var_pts(pts_name_to_variable_number('TEMPK')) = .true.
        postprocess_var_pts(pts_name_to_variable_number('MASSK')) = .true.
     end if

     deposition_var_pts( pts_name_to_variable_number(    'T')) = .true.
     deposition_var_pts( pts_name_to_variable_number('ILAGR')) = .true.
     deposition_var_pts( pts_name_to_variable_number('ITYPE')) = .true.
     deposition_var_pts( pts_name_to_variable_number('COORX')) = .true.
     deposition_var_pts( pts_name_to_variable_number('COORY')) = .true.
     deposition_var_pts( pts_name_to_variable_number('COORZ')) = .true.
     deposition_var_pts( pts_name_to_variable_number('EXIST')) = .true.
     deposition_var_pts( pts_name_to_variable_number(  'SET')) = .true.
    
     !postprocess_var_pts( 1: 9) = .true.    ! t,ilagr,itype,coord,veloc,dt_k
     !postprocess_var_pts(14:16) = .true.    ! dt_k,kfl_exist,Cd     
     !deposition_var_pts(  1: 6) = .true.    ! t,ilagr,itype,coord
     !deposition_var_pts(    15) = .true.    ! kfl_exist
     !deposition_var_pts(    31) = .true.    ! boundary_set

  end select

end subroutine pts_inivar
