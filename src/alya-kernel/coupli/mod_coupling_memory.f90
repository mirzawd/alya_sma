!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Coupling
!> Manage the memory coupling
!> @{
!> @file    mod_coupling_memory.f90
!> @author  houzeaux
!> @date    2018-03-27
!> @brief   Coupling memory
!> @details Manages the memory of the coupling structure
!>          
!-----------------------------------------------------------------------

module mod_coupling_memory

  use def_kintyp,            only : ip,rp,lg
  use def_spmat,             only : spmat
  use def_coupli,            only : typ_color_coupling
  use def_coupli,            only : typ_coupling_geometry
  use def_coupli,            only : typ_coupling_wet
  use def_master,            only : AT_BEGINNING
  use def_coupli,            only : memor_cou
  use def_coupli,            only : RELAXATION_SCHEME
  use def_coupli,            only : VALUES_ON_NODES
  use def_coupli,            only : UNKNOWN
  use mod_kdtree,            only : kdtree_deallocate
  use mod_kdtree,            only : kdtree_initialize
  use mod_memory,            only : memory_alloca
  use mod_memory,            only : memory_deallo
  use mod_optional_argument, only : optional_argument
  implicit none
  private

  interface cou_initialization
     module procedure cou_initialization_array,&
          &           cou_initialization_scalar 
  end interface cou_initialization
 
  interface cou_deallocate
     module procedure cou_deallocate_array,&
          &           cou_deallocate_scalar 
  end interface cou_deallocate
 
  public :: cou_initialization
  public :: cou_deallocate
  public :: COU_DEALLOCATE_GEOMETRY   
  public :: COU_DEALLOCATE_TRANSMISSION          
  public :: COU_DEALLOCATE_SINGLE_COUPLING
  public :: COU_INITIALIZATION_SINGLE_COUPLING
  public :: COU_REINITIALIZATION_ALL_COUPLINGS
  
contains
  
  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    12/03/2018
  !> @brief   Initialize coupling  
  !> @details Initialize coupling variables of all couplings
  !>
  !----------------------------------------------------------------------
 
  subroutine cou_initialization_array(coupling)
 
    type(typ_color_coupling), intent(inout), pointer :: coupling(:)
    integer(ip) :: icoup
    
    do icoup = 1,size(coupling)
       call COU_INITIALIZATION_SINGLE_COUPLING(coupling(icoup))
    end do

  end subroutine cou_initialization_array

  subroutine cou_initialization_scalar(coupling)
 
    type(typ_color_coupling),           intent(inout) :: coupling
   
    call COU_INITIALIZATION_SINGLE_COUPLING(coupling)

  end subroutine cou_initialization_scalar

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    12/03/2018
  !> @brief   Initialize coupling  
  !> @details Deallocate coupling structure
  !>
  !----------------------------------------------------------------------
 
  subroutine cou_deallocate_array(coupling,COUPLING_NAME)
 
    type(typ_color_coupling),           intent(inout), pointer :: coupling(:)
    character(LEN=*),         optional, intent(in)             :: COUPLING_NAME
    integer(ip)                                                :: icoup
    
    do icoup = 1,size(coupling)
       call COU_DEALLOCATE_SINGLE_COUPLING(coupling(icoup),COUPLING_NAME=COUPLING_NAME)
    end do

  end subroutine cou_deallocate_array

  subroutine cou_deallocate_scalar(coupling,COUPLING_NAME)
 
    type(typ_color_coupling),           intent(inout) :: coupling
    character(LEN=*),         optional, intent(in)    :: COUPLING_NAME
    
    call COU_DEALLOCATE_SINGLE_COUPLING(coupling,COUPLING_NAME=COUPLING_NAME)

  end subroutine cou_deallocate_scalar

  !---------------------------------------------------------------------- 
  !
  !> @author  Matias Rivero
  !> @date    16/12/2014
  !> @brief   Deallocate coupling structure (geome)
  !> @details Deallocate coupling structure (geome)
  !>
  !----------------------------------------------------------------------

  subroutine COU_DEALLOCATE_WET(wet,SAVE_HISTORIC)

    type(typ_coupling_wet), intent(inout) :: wet
    logical(lg), optional,  intent(in)    :: SAVE_HISTORIC
    integer(ip)                           :: iboun_wet

    wet % number_wet_points = 0
    wet % npoin_wet         = 0
    wet % nboun_wet         = 0
    
    call memory_deallo(memor_cou,'COUPLING % WET % VMASB_WET'     ,'DEALLOCATE_COUPLING_STRUCTURE_WET',wet % vmasb_wet      )
    call memory_deallo(memor_cou,'COUPLING % WET % LBOUN_WET'     ,'DEALLOCATE_COUPLING_STRUCTURE_WET',wet % lboun_wet      )
    call memory_deallo(memor_cou,'COUPLING % WET % LSUBD'         ,'DEALLOCATE_COUPLING_STRUCTURE_WET',wet % lsubd          )
    call memory_deallo(memor_cou,'COUPLING % WET % LPOIN_WET'     ,'DEALLOCATE_COUPLING_STRUCTURE_WET',wet % lpoin_wet      )
    call memory_deallo(memor_cou,'COUPLING % WET % COORD_WET'     ,'DEALLOCATE_COUPLING_STRUCTURE_WET',wet % coord_wet      )
    call memory_deallo(memor_cou,'COUPLING % WET % COORD_WET_INI' ,'DEALLOCATE_COUPLING_STRUCTURE_WET',wet % coord_wet_ini  )
    call memory_deallo(memor_cou,'COUPLING % WET % WEIGHT_WET'    ,'DEALLOCATE_COUPLING_STRUCTURE_WET',wet % weight_wet     )
    call memory_deallo(memor_cou,'COUPLING % WET % WEIGHT_WET_IMP','DEALLOCATE_COUPLING_STRUCTURE_WET',wet % weight_wet_imp )
    call memory_deallo(memor_cou,'COUPLING % WET % MASS_RATIO_WET','DEALLOCATE_COUPLING_STRUCTURE_WET',wet % mass_ratio_wet )

    if( associated(wet % proje_target) ) then
       do iboun_wet = 1,size(wet % proje_target)
          if( associated(wet % proje_target(iboun_wet) % permu) ) deallocate(wet % proje_target(iboun_wet) % permu)
          if( associated(wet % proje_target(iboun_wet) % shapb) ) deallocate(wet % proje_target(iboun_wet) % shapb)
          if( associated(wet % proje_target(iboun_wet) % gbsur) ) deallocate(wet % proje_target(iboun_wet) % gbsur)
       end do
       deallocate(wet % proje_target)
    end if

  end subroutine COU_DEALLOCATE_WET

  !---------------------------------------------------------------------- 
  !
  !> @author  Matias Rivero
  !> @date    16/12/2014
  !> @brief   Deallocate coupling structure (geome)
  !> @details Deallocate coupling structure (geome)
  !>
  !----------------------------------------------------------------------

  subroutine COU_DEALLOCATE_GEOMETRY(geome)

    type(typ_coupling_geometry), intent(inout) :: geome

    geome % ndime     = 0

    geome % nelem_source = 0
    geome % nboun_source = 0
    geome % npoin_source = 0

    call memory_deallo(memor_cou,'COUPLING % GEOME % LELEM_SOURCE', 'DEALLOCATE_COUPLING_STRUCTURE_GEOME',geome % lelem_source )
    call memory_deallo(memor_cou,'COUPLING % GEOME % LBOUN_SOURCE', 'DEALLOCATE_COUPLING_STRUCTURE_GEOME',geome % lboun_source )
    call memory_deallo(memor_cou,'COUPLING % GEOME % LPOIN_SOURCE', 'DEALLOCATE_COUPLING_STRUCTURE_GEOME',geome % lpoin_source )
    call memory_deallo(memor_cou,'COUPLING % GEOME % STATUS'      , 'DEALLOCATE_COUPLING_STRUCTURE_GEOME',geome % status       )
    call memory_deallo(memor_cou,'COUPLING % GEOME % SHAPF'       , 'DEALLOCATE_COUPLING_STRUCTURE_GEOME',geome % shapf        )    
    call memory_deallo(memor_cou,'COUPLING % GEOME % VMASB'       , 'DEALLOCATE_COUPLING_STRUCTURE_GEOME',geome % vmasb        )
    call memory_deallo(memor_cou,'COUPLING % GEOME % SCHED_PERM'  , 'DEALLOCATE_COUPLING_STRUCTURE_GEOME',geome % sched_perm   )

    call kdtree_deallocate(geome % kdtree)

  end subroutine COU_DEALLOCATE_GEOMETRY

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/04/2014
  !> @brief   Initialize one single coupling 
  !> @details Initialize one single coupling 
  !>
  !----------------------------------------------------------------------
 
  subroutine COU_INITIALIZATION_SINGLE_COUPLING(coupling,WET,INPUT)

    use def_coupli,      only : NODE_TARGET_ENTITY

    type(typ_color_coupling), intent(inout)        :: coupling
    logical(lg),              intent(in), optional :: WET
    logical(lg),              intent(in), optional :: INPUT

    call coupling % init(WET,INPUT)
    
  end subroutine COU_INITIALIZATION_SINGLE_COUPLING

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/04/2014
  !> @brief   Initialize coupling 
  !> @details Initialize coupling transmission matrices 
  !>
  !----------------------------------------------------------------------

  subroutine COU_INITIALIZATION_TRANSMISSION(ltransmat)
    
    type(spmat), pointer, intent(inout) :: ltransmat(:) ! Transmission matrix
!    integer(ip)                         :: ii

    nullify(ltransmat)
    !if( associated(ltransmat) ) then
    !   do ii = 1,size(ltransmat)
    !      ltransmat(ii) % ndof  = 0
    !      ltransmat(ii) % nrows = 0
    !      ltransmat(ii) % ncols = 0
    !      
    !      nullify(ltransmat(ii) % iA)
    !      nullify(ltransmat(ii) % jA)
    !      nullify(ltransmat(ii) % vA)
    !   end do
    !end if
    
  end subroutine COU_INITIALIZATION_TRANSMISSION
  
  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/04/2014
  !> @brief   Initialize coupling 
  !> @details Initialize coupling geometry
  !>
  !----------------------------------------------------------------------

  subroutine COU_INITIALIZATION_WET(wet)

    use def_coupli,      only : NODE_WET_POINT
    type(typ_coupling_wet), intent(inout) :: wet

    call wet % init()

    !wet % number_wet_points = 0  
    !wet % npoin_wet         = 0   
    !wet % nboun_wet         = 0
    !wet % point_type        = NODE_WET_POINT

    !nullify(wet % lboun_wet)               
    !nullify(wet % lpoin_wet)               
    !nullify(wet % coord_wet)             
    !nullify(wet % weight_wet)              
    !nullify(wet % weight_wet_imp)              
    !nullify(wet % vmasb_wet)               
    !nullify(wet % mass_ratio_wet)          
    !nullify(wet % proje_target)

  end subroutine COU_INITIALIZATION_WET
  
  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/04/2014
  !> @brief   Initialize coupling 
  !> @details Initialize coupling geometry
  !>
  !----------------------------------------------------------------------

  subroutine COU_INITIALIZATION_GEOMETRY(geome)

    use def_coupli, only : typ_coupling_geometry
    type(typ_coupling_geometry), intent(inout) :: geome

    call geome % init()
    
    !geome % ndime             = 0     
    !geome % nelem_source      = 0      
    !geome % nboun_source      = 0         
    !geome % npoin_source      = 0 
    
    !nullify(geome % lelem_source)            
    !nullify(geome % lboun_source)            
    !nullify(geome % lpoin_source)            
    !nullify(geome % status)                  
    !nullify(geome % shapf)               
    !nullify(geome % vmasb)
    !nullify(geome % sched_perm) 
    
    !call kdtree_initialize(geome % kdtree) 

  end subroutine COU_INITIALIZATION_GEOMETRY

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/04/2014
  !> @brief   Initialize coupling 
  !> @details Initialize coupling transmission matrices 
  !>
  !----------------------------------------------------------------------

  subroutine COU_DEALLOCATE_TRANSMISSION(ltransmat)

    type(spmat), pointer, intent(inout) :: ltransmat(:)     ! Transmission matrix

    if( associated(ltransmat) ) then
       ltransmat(:) % ndof1 = 0
       ltransmat(:) % ndof2 = 0
       ltransmat(:) % nrows = 0
       ltransmat(:) % ncols = 0  
       call memory_deallo(memor_cou,'LTRANSMAT_TARGET','COU_DEALLOCATE_SINGLE_COUPLING',ltransmat)
    end if
    
  end subroutine COU_DEALLOCATE_TRANSMISSION
  
  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/04/2014
  !> @brief   Initialize one single coupling 
  !> @details Initialize one single coupling 
  !>
  !----------------------------------------------------------------------
 
  subroutine COU_DEALLOCATE_SINGLE_COUPLING(coupling,WET,COUPLING_NAME)

    type(typ_color_coupling),           intent(inout)        :: coupling
    logical(lg),                        intent(in), optional :: WET
    character(LEN=*),         optional, intent(in)           :: COUPLING_NAME
    logical(lg)                                              :: ifwet
    character(LEN=:),         allocatable                    :: my_name

    my_name = optional_argument('coupling',COUPLING_NAME)
    ifwet   = optional_argument(.true.,WET)
    
    call memory_deallo(memor_cou,my_name//' % values'                  ,'COU_DEALLOCATE_SINGLE_COUPLING',coupling % values          )
    call memory_deallo(memor_cou,my_name//' % values_frequ'            ,'COU_DEALLOCATE_SINGLE_COUPLING',coupling % values_frequ    )     
    call memory_deallo(memor_cou,my_name//' % values_predicted'        ,'COU_DEALLOCATE_SINGLE_COUPLING',coupling % values_predicted)   
    call memory_deallo(memor_cou,my_name//' % values_converged'        ,'COU_DEALLOCATE_SINGLE_COUPLING',coupling % values_converged) 
    call memory_deallo(memor_cou,my_name//' % jacobian_inverse'        ,'COU_DEALLOCATE_SINGLE_COUPLING',coupling % jacobian_inverse) 
    call memory_deallo(memor_cou,my_name//' % dincr_predicted'         ,'COU_DEALLOCATE_SINGLE_COUPLING',coupling % dincr_predicted )      
    call memory_deallo(memor_cou,my_name//' % residues_iqnls'          ,'COU_DEALLOCATE_SINGLE_COUPLING',coupling % residues_iqnls  )    
    call memory_deallo(memor_cou,my_name//' % relaxed_iqnls'           ,'COU_DEALLOCATE_SINGLE_COUPLING',coupling % relaxed_iqnls   )     
    call memory_deallo(memor_cou,my_name//' % unrelaxed_iqnls'         ,'COU_DEALLOCATE_SINGLE_COUPLING',coupling % unrelaxed_iqnls )   
    call memory_deallo(memor_cou,my_name//' % valincr_iqnls'           ,'COU_DEALLOCATE_SINGLE_COUPLING',coupling % valincr_iqnls   )     
    call memory_deallo(memor_cou,my_name//' % residincr_iqnls'         ,'COU_DEALLOCATE_SINGLE_COUPLING',coupling % residincr_iqnls )    
    call memory_deallo(memor_cou,my_name//' % valincr_history_iqnls'   ,'COU_DEALLOCATE_SINGLE_COUPLING',coupling % valincr_history_iqnls )    
    call memory_deallo(memor_cou,my_name//' % residincr_history_iqnls' ,'COU_DEALLOCATE_SINGLE_COUPLING',coupling % residincr_history_iqnls )    
    call memory_deallo(memor_cou,my_name//' % V_current_history_iqnls' ,'COU_DEALLOCATE_SINGLE_COUPLING',coupling % V_current_history_iqnls )    
    call memory_deallo(memor_cou,my_name//' % W_current_history_iqnls' ,'COU_DEALLOCATE_SINGLE_COUPLING',coupling % W_current_history_iqnls )    
    call memory_deallo(memor_cou,my_name//' % history_tracking_iqnls'  ,'COU_DEALLOCATE_SINGLE_COUPLING',coupling % history_tracking_iqnls )    
    
    call COU_DEALLOCATE_GEOMETRY           (coupling % geome)
    if( ifwet ) call COU_DEALLOCATE_WET    (coupling % wet)

    call coupling % commd % deallo(memor_cou,PAR_COMM_OPT=.false.,COMM_NAME='COUPLING % COMMD',INITIALIZE=.true.)
    
    !call PAR_DEALLOCATE_COMMUNICATION_ARRAY(coupling % commd,memor_cou,PAR_COMM_OPT=.false.,COMM_NAME='COUPLING % COMMD',INITIALIZE=.true.)
    call COU_DEALLOCATE_TRANSMISSION       (coupling % ltransmat_target)   

    if( allocated(my_name) ) deallocate(my_name)
    
  end subroutine COU_DEALLOCATE_SINGLE_COUPLING

  !----------------------------------------------------------------------
  !
  !> @author  David Oks
  !> @date    23/03/2022
  !> @brief   Renitialize all couplings
  !> @details Renitialize all couplings
  !>
  !----------------------------------------------------------------------

  subroutine COU_REINITIALIZATION_ALL_COUPLINGS(WET,INPUT)

    use def_coupli, only : coupling_type, mcoup

    implicit none

    logical(lg), intent(in), optional :: WET
    logical(lg), intent(in), optional :: INPUT
    integer(ip)                       :: icoup
    logical(lg)                       :: ifwet
    logical(lg)                       :: ifinput

    ifwet   = optional_argument(.true.,WET)
    ifinput = optional_argument(.true.,INPUT)

    do icoup = 1,mcoup
       call COU_DEALLOCATE_SINGLE_COUPLING(coupling_type(icoup),WET=ifwet)
       call COU_INITIALIZATION_SINGLE_COUPLING(coupling_type(icoup),WET=ifwet,INPUT=ifinput)
    end do

  end subroutine COU_REINITIALIZATION_ALL_COUPLINGS

end module mod_coupling_memory
!> @}
