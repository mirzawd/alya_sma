!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ker_inivar(itask)
  !-----------------------------------------------------------------------
  !****f* Kermod/ker_inivar
  ! NAME
  !    ker_inicar
  ! DESCRIPTION
  !    This routine initializes some variables
  !    ITASK=1 ... When starting the run (from Turnon)
  !    ITASK=2 ... First time step. This is needed as some variables
  !                are not initialized before
  !    ITASK=3 ... When starting a time step (from ker_begste)
  ! USES
  ! USED BY
  !    ker_turnon
  !***
  !-----------------------------------------------------------------------
#include "def_vector_size.inc"
  use def_parame
  use def_elmtyp
  use def_master
  use def_kermod
  use def_domain
  use def_solver
  use mod_cutele
  use mod_ker_space_time_function
  use mod_elsest
  use mod_ker_detection,          only : ker_events_directory_name
  use mod_ker_subdomain,          only : ker_subdomain_initialization
  use mod_ker_proper,             only : ker_proper_initialization
  use mod_ker_proper,             only : ker_proper_default_values
  use mod_ker_proper,             only : ker_proper_allocate_material_laws
  use mod_ker_proper,             only : ker_proper_smoothing
  use mod_ker_discrete_function,  only : ker_discrete_function_initialization
  use mod_arrays,                 only : arrays_register  
  use mod_interp_tab,             only : tab_init_tab 
  use mod_interp_tab,             only : tab_init_coord
  use mod_interp_tab,             only : tab_init_fw
  use mod_eccoupling,             only : kfl_exmsld_ecc, eccou_manage_arrays, eccou_initialise_flags, eccou_allocate_memory
  use mod_biofibers,              only : biofibers
  use mod_random,                 only : random_initialization
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ifunc,jfunc

  select case ( itask )

  case ( 0_ip )
     !
     ! Initialize modules
     !
     call ker_proper_initialization()
     call ker_subdomain_initialization()
     call ker_discrete_function_initialization()
     !
     ! Postprocess initial solution
     !
     postp(1) % npp_iniso = 1
     !
     ! Register primary variables
     !
     call arrays_register((/'VELAV','MATRI','NBOUN','PRIMA'/),velav_ker,           ENTITY_POSITION=3_ip)                    ! Velocity at boundaries                                       
     call arrays_register((/'AVUPO','VECTO','NPOIN','PRIMA'/),avupo_ker,           ENTITY_POSITION=2_ip)                    ! For time-averaging restart              
     call arrays_register((/'AVTA1','VECTO','NELEM','PRIMA'/),avta1_nsw_ker,       ENTITY_POSITION=2_ip)                    ! Average viscous(lam+tur) part of tangential stress
     
     call arrays_register((/'ADVEC','VECTO','NPOIN','PRIMA'/),advec,               ENTITY_POSITION=2_ip,TIME_POSITION=3_ip)
     call arrays_register((/'DISPM','VECTO','NPOIN','PRIMA'/),dispm,               ENTITY_POSITION=2_ip,TIME_POSITION=3_ip)
     call arrays_register((/'CONCE','SCALA','NPOIN','PRIMA'/),conce,               ENTITY_POSITION=1_ip,TIME_POSITION=3_ip)
     call arrays_register((/'TEMPE','SCALA','NPOIN','PRIMA'/),tempe,               ENTITY_POSITION=1_ip,TIME_POSITION=2_ip)
     call arrays_register((/'AREAS','SCALA','NPOIN','PRIMA'/),areas,               ENTITY_POSITION=1_ip,TIME_POSITION=2_ip)
     !
     ! Register arrays forofibers
     !
     call arrays_register((/'BFIBL','VECTO','NPOIN','PRIMA'/),biofibers % nod_lng, ENTITY_POSITION=2_ip,TIME_POSITION=3_ip)
     call arrays_register((/'BFIBS','VECTO','NPOIN','PRIMA'/),biofibers % nod_sht, ENTITY_POSITION=2_ip,TIME_POSITION=3_ip)
     call arrays_register((/'BFIBN','VECTO','NPOIN','PRIMA'/),biofibers % nod_nrm, ENTITY_POSITION=2_ip,TIME_POSITION=3_ip)
     !
     ! Register arrays for the electro-mechanical coupling {exmedi-solidz}
     ! 
     call eccou_manage_arrays(0_ip)    
     !
     ! Secondary variables
     !
     call arrays_register((/'EXNOR','VECTO','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'ERRNO','SCALA','NPOIN','SECON'/))
     call arrays_register((/'LPOIN','SCALA','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'SKCO1','VECTO','NPOIN','SECON'/))
     call arrays_register((/'SKCO2','VECTO','NPOIN','SECON'/))
     call arrays_register((/'SKCO3','VECTO','NPOIN','SECON'/))
     call arrays_register((/'LPOTY','SCALA','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'PONUM','SCALA','NPOIN','SECON'/))
     call arrays_register((/'CODNO','VECTO','NPOIN','SECON'/))
     call arrays_register((/'YWALP','SCALA','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'LNTIB','SCALA','NPOIN','SECON'/))
     call arrays_register((/'HOLES','SCALA','NPOIN','SECON'/))
     call arrays_register((/'DENSI','SCALA','NPOIN','SECON'/))
     call arrays_register((/'VISCO','SCALA','NPOIN','SECON'/))
     call arrays_register((/'POROS','SCALA','NPOIN','SECON'/))
     call arrays_register((/'CONDU','SCALA','NPOIN','SECON'/))
     call arrays_register((/'SPECI','SCALA','NPOIN','SECON'/))
     call arrays_register((/'GROUP','SCALA','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'MASSM','SCALA','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'MASSC','SCALA','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'GEONO','SCALA','NPOIN','SECON'/))
     call arrays_register((/'SUBDO','SCALA','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'WALLD','SCALA','NPOIN','SECON'/))
     call arrays_register((/'ROUGH','SCALA','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'KEKET','SCALA','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'CODBO','SCALA','NBOUN','SECON'/))
     call arrays_register((/'MATER','SCALA','NELEM','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'LMATN','SCALA','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'COMMU','MATRI','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'COMMI','MATRI','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'LELEV','SCALA','NELEM','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'ELNUM','SCALA','NELEM','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'CODBB','SCALA','NPOIN','SECON'/))
     call arrays_register((/'LETIB','SCALA','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'LTYP2','SCALA','NPOIN','SECON'/))
     call arrays_register((/'LELC2','SCALA','NPOIN','SECON'/))
     call arrays_register((/'CONTA','SCALA','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'RDOM ','SCALA','NPOIN','SECON'/))
     call arrays_register((/'VORTX','VECTO','NPOIN','SECON'/))
     call arrays_register((/'LESET','SCALA','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'DISMM','VECTO','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'VELOC','VECTO','NPOIN','SECON'/))
     call arrays_register((/'LBSET','SCALA','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'LMESH','SCALA','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'TURBU','SCALA','NPOIN','SECON'/))
     call arrays_register((/'LNSUB','SCALA','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'WETNO','SCALA','NPOIN','SECON'/))
     call arrays_register((/'LESUB','SCALA','NELEM','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'CANOP','SCALA','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'HEIGH','SCALA','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'BEATR','SCALA','NELEM','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'COLOR','SCALA','NELEM','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'WALLN','VECTO','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'LOCAL','SCALA','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'RENPO','SCALA','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'OMPSS','SCALA','NELEM','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'OPENM','SCALA','NELEM','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'GAUSS','SCALA','NELEM','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'WALLO','SCALA','NPOIN','SECON'/))
     call arrays_register((/'MIXIN','SCALA','NPOIN','SECON'/))
     call arrays_register((/'FIELD','SCALA','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)    
     call arrays_register((/'VELOM','VECTO','NPOIN','SECON'/))
     call arrays_register((/'OMPSB','SCALA','NBOUN','SECON'/))
     call arrays_register((/'NSWVI','SCALA','NPOIN','SECON'/))
     call arrays_register((/'NUMBE','SCALA','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'MPIRA','SCALA','NELEM','SECON'/),POST_ONLY_ONCE=.true.)   
     call arrays_register((/'VELAE','VECTO','NELEM','SECON'/))
     call arrays_register((/'LAD  ','SCALA','NELEM','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'INTNO','SCALA','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)     
     call arrays_register((/'VELNO','VECTO','NPOIN','SECON'/))
     call arrays_register((/'FRING','SCALA','NPOIN','SECON'/))
     call arrays_register((/'DAVID','SCALA','NPOIN','SECON'/))
     call arrays_register((/'GRWAL','VECTO','NPOIN','SECON'/),POST_ONLY_ONCE=.true.)
     call arrays_register((/'ANIPO','VECT2','NELEM','SECON'/),POST_ONLY_ONCE=.true.)
     !
     ! Nullify pointers
     !
     nullify(lnodb_mm)
     nullify(coord_mm)
     nullify(tncod_ker)
     nullify(tgcod_ker)
     nullify(tbcod_ker)
     nullify(cowit)
     nullify(cowit_origi)
     nullify(uwall_ker)
     nullify(uwal2_ker)
     nullify(shwit)
     nullify(dewit)
     nullify(displ_ker)
     nullify(lewit)
     nullify(kfl_funno_walld_ker)
     nullify(kfl_funbo_walld_ker)
     nullify(kfl_fixno_walld_ker)
     nullify(kfl_fixbo_walld_ker)
     nullify(kfl_funty_walld_ker)
     nullify(kfl_fixno_walln_ker)
     nullify(kfl_fixbo_walln_ker)
     nullify(funpa_walld_ker)
     nullify(bvess_walld_ker)
     nullify(bvnat_walld_ker)
     nullify(bvess_walln_ker)
     nullify(kfl_funno_defor_ker)
     nullify(kfl_fixno_defor_ker)
     nullify(kfl_funty_defor_ker)
     nullify(bvess_defor_ker)
     nullify(kfl_funno_rough_ker)
     nullify(kfl_funbo_rough_ker)
     nullify(kfl_fixno_rough_ker)
     nullify(kfl_fixbo_rough_ker)
     nullify(kfl_fixbo_nsw_ker)
     nullify(kfl_funty_rough_ker)
     nullify(funpa_rough_ker)
     nullify(bvess_rough_ker)
     nullify(bvnat_rough_ker)
     nullify(kfl_fixno_suppo_ker)
     nullify(bvess_suppo_ker)
     nullify(avta1_nsw_ker)
     nullify(velav_ker)
     nullify(avupo_ker)
     nullify(lnsw_exch)
     nullify(velel_ker)
     nullify(temel_ker)
     nullify(lexlo_ker)
     nullify(shape_waexl)
     nullify(kfl_boexc_ker)
     nullify(kfl_dampi_ker)
     nullify(subdomain)
     nullify(gewit)
     nullify(witness_mesh)
     nullify(witness_mesh)
     call search_waexlo_seq % init()
     call search_waexlo_par % init()
     call interp_waexlo     % init()
     call search_elsest_seq % init()
     call search_elsest_par % init()
     call interp_elsest     % init()
     !
     ! Initialize modules
     !
     call ker_proper_initialization()
     call ker_subdomain_initialization()     

     do ifunc = 1,max_space_time_function
        nullify( space_time_function(ifunc) % expression )
     end do
     do ifunc = 1,max_time_function
        nullify( time_function(ifunc) % parameters )
     end do
     do ifunc = 1,max_windk_systems
        nullify( windk_systems(ifunc) % params )
        nullify( windk_systems(ifunc) % xprev )
        nullify( windk_systems(ifunc) % yprev )
         windk_systems(ifunc) % yrelaxed   = 0.0_rp
         windk_systems(ifunc) % yunrelaxed = 0.0_rp
         windk_systems(ifunc) % w          = 0.0_rp
     end do
     do ifunc = 1,max_lookup_tab
        call tab_init_tab(lookup_tab(ifunc))
        do jfunc = 1,max_lookup_dim
           call tab_init_coord(lookup_coords(jfunc,ifunc))
        end do
     end do
     do ifunc = 1,max_lookup_fw
        call tab_init_fw(lookup_fw(ifunc))
     end do
     if( .not. associated(ann_fw) ) then
       allocate(ann_fw(max_ann_fw))
     end if
     do ifunc = 1,max_ann_fw
        call ann_fw(ifunc) % init() 
     end do
#ifdef TORCH
     call initialize_neural_networks(max_ann_fw)
#endif
     !
     ! Solvers - kernel problems
     !
     call moddef( 9_ip)
     call soldef(-6_ip)   ! Allocate memory for 7 solvers

     solve(1) % ndofn     = 1
     solve(1) % kfl_solve = 1
     solve(1) % wprob     = 'ROUGHNESS'

     solve(2) % ndofn     = 1
     solve(2) % kfl_solve = 1
     solve(2) % wprob     = 'WALL_DISTANCE'

     solve(3) % ndofn     = ndime
     solve(3) % kfl_solve = 1
     solve(3) % wprob     = 'SUPPORT_GEOMETRY'

     solve(4) % ndofn     = ndime
     solve(4) % kfl_solve = 1
     solve(4) % wprob     = 'MESH_DEFORMATION'

     solve(5) % ndofn     = 1
     solve(5) % kfl_solve = 1
     solve(5) % wprob     = 'WALL_NORMAL'

     solve(6) % ndofn     = 1
     solve(6) % kfl_solve = 1
     solve(6) % wprob     = 'GROUPS'
     solve(6) % kfl_iffix = 1     
     !
     ! Laws for properties
     !
     call ker_proper_allocate_material_laws()
     call ker_proper_default_values()
     call ker_allaws()
     !
     ! Others
     !
     pres_set_history=0.0_rp
     number_event = 0
     !
     ! Get events directory name
     !
     call ker_events_directory_name()

  case ( 1_ip )
     !
     ! Redefine MGAUS if there are cut elements
     !
     if( kfl_cutel == 1 ) then
        if( ndime == 2 ) then
           if( lexis(TRI03) == 0 ) call runend('CDERDA: WHEN USING CUT ELEMENTS, DECLARE TRI03 ELEMENTS')
           mgaus = max(mgaus,9*ngaus(TRI03))
        else
           if( lexis(HEX08) /= 0 ) then
              mgaus = max(mgaus,36*ngaus(TET04))
           else
              if( lexis(TET04) == 0 ) call runend('CDERDA: WHEN USING CUT ELEMENTS, DECLARE TET04 ELEMENTS')
              mgaus = max(mgaus,6*ngaus(TET04))
           end if
        end if
     end if
     !
     ! Initialize space/time functions
     !
     call ker_init_space_time_function()
     !
     ! Electro-mechanical coupling
     !
     if( kfl_exmsld_ecc )then
        call eccou_initialise_flags()
        call eccou_allocate_memory(4_ip)
     endif
     !
     ! Read neural networks on all MPI processes
     !
     do ifunc = 1,max_ann_fw
        call ann_fw(ifunc) % read_file() 
     end do
     !
     ! Random number initialization
     !
     call random_initialization(random_from_file=random_from_file)

  case ( 2_ip )
     !
     ! Allocate memory for cut elements
     !
     if( kfl_cutel == 1 ) then
        if( INOTMASTER ) call buicut(1_ip)
        mnoga = max(mnode,mgaus)
     end if

  end select

end subroutine ker_inivar
