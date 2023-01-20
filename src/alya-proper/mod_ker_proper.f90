!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @defgroup Properties_Toolbox
!> Toolbox for creating, updating and computing properties
!> @{
!> @name    ToolBox for property creations
!> @file    mod_ker_proper.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for properties
!> @details To create a new property:
!>          - Change variable NUMBER_OF_PROPERTIES
!>          - Add your property ****_ker in def_kermod.f90
!>          - Look for densi_ker... and add yours
!>          - ker_readat.f90: Add the reading of the property
!>          - Add laws in ker_allaws.f90
!>          _ Add property to ker_parall.f90
!
!------------------------------------------------------------------------

module mod_ker_proper

#include "def_vector_size.inc"
  use def_kintyp,             only : ip,rp,r1p
  use def_parame,             only : pi,zero,one
  use def_master,             only : itinn, ittim
  use mod_parall,             only : PAR_MY_CODE_RANK
  use def_master,             only : INOTMASTER
  use def_master,             only : modul,mem_modul
  use def_master,             only : ITASK_DOITER
  use def_master,             only : ITASK_ENDSTE
  use def_master,             only : ITASK_INIUNK
  use def_master,             only : INOTMASTER
  use def_master,             only : condu_gp,sphec_gp,visco_gp
  use def_master,             only : mmodu,kfl_modul
  use def_master,             only : namod
  use def_master,             only : ID_NASTIN
  use def_master,             only : ID_TEMPER
  use def_master,             only : ID_PARTIS
  use def_master,             only : nturb
  use def_domain,             only : nmate,nelem,npoin,ltypb
  use def_domain,             only : ltype,lmate,ngaus,nboun
  use def_domain,             only : lgaus,ndime,mnode,lnods
  use def_domain,             only : ndime,nnode,kfl_naxis
  use def_domain,             only : mgaus,elmar,coord,lelch
  use def_domain,             only : lnnod,lnodb,nmatn,lelbo
  use def_domain,             only : lmatn,hnatu,ntens,canhe
  use def_domain,             only : canla,walld,heiov,rough
  use def_domain,             only : mnoga
  use def_elmtyp,             only : ELEXT,ELCUT
  use def_elmtyp,             only : TET04
  use def_elmtyp,             only : TRI03
  use mod_parall,             only : par_omp_nelem_chunk
  use mod_ker_polynomial,     only : ker_polynomial_proper
  use mod_ker_polynomial,     only : ker_polynomial_name 
  use mod_memory,             only : memory_alloca
  use mod_memory,             only : memory_deallo
  use mod_ker_regularization, only : regul_k, regul_e, kfl_regularization
  use mod_interpolation,      only : COU_GET_INTERPOLATE_POINTS_VALUES
  use mod_communications,     only : PAR_EXCHANGE
  use mod_ker_proper_generic, only : ker_proper_element_operations
  use mod_exchange,           only : exchange_add
  use def_ker_proper
  use def_master
  use def_kermod
  implicit none
  private

  interface ker_proper
     module procedure ker_proper_scalar_000, &
          &           ker_proper_scalar_00,  &
          &           ker_proper_scalar_0,   &
          &           ker_proper_scalar_1,   &
          &           ker_proper_scalar_2,   &
          &           ker_proper_vector_0,   &
          &           ker_proper_vector_1,   &
          &           ker_proper_vector_2,   &
          &           ker_proper_vector_3
  end interface ker_proper

  interface ker_proper_pointer
     module procedure ker_proper_pointer_w, &
          &           ker_proper_pointer_ip
  end interface ker_proper_pointer

  public :: ker_proper_initialization          
  public :: ker_proper_default_values
  public :: ker_proper_allocate_material_laws
  public :: ker_proper_allocate_properties
  public :: ker_proper_property_array_exists
  public :: ker_proper
  public :: ker_proper_updated_property 
  public :: ker_proper_exists                  ! Return if a property exists
  public :: ker_proper_check_properties        ! Check if appropriate properties have been defined
  public :: ker_proper_parall                  ! Boradcast properties inputs
  public :: ker_proper_updpro_time
  public :: ker_proper_smoothing               ! Determine if smoothing is required
  public :: ker_proper_number_properties       ! Number of existing properties
  public :: ker_proper_on_the_fly              ! Set global on-the-fly option
  public :: ker_proper_pointer                 ! Pointer to a property
  
  public :: NUMBER_OF_PROPERTIES
  public :: SCALAR_PROPERTY
  public :: VECTOR_PROPERTY
  public :: MATRIX_PROPERTY
  public :: UPDATE_PROPERTIES
  public :: cpu_prope

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-22
  !> @brief   Smoothing
  !> @details Determine if smoothing is required
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_proper_smoothing(prope_ker) 

    type(typ_valpr_ker), intent(inout) :: prope_ker
!    integer(ip)                        :: ilaws,imate

    !if( nmate > 1 ) prope_ker % kfl_nedsm = 1
    !do imate = 1,nmate
    !   ilaws = prope_ker % ilaws(imate)
    !   if( prope_ker % llaws(ilaws) % where == 'IELEM' ) prope_ker % kfl_nedsm = 1
    !end do

  end subroutine ker_proper_smoothing
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-22
  !> @brief   Get timings
  !> @details Get timings of ker_proper
  !> 
  !-----------------------------------------------------------------------

  real(rp) function ker_proper_updpro_time() 

    ker_proper_updpro_time = cpu_prope(1)
    
  end function ker_proper_updpro_time
    
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-17
  !> @brief   Check property array needed
  !> @details Return if a given array is required according to the law
  !> 
  !-----------------------------------------------------------------------

  logical(lg) function ker_proper_property_array_exists(warray,prope_ker)

    character(len=*),             intent(in)    :: warray
    type(typ_valpr_ker), pointer, intent(inout) :: prope_ker     
    integer(ip)                                 :: imate,ilaws
    integer(ip)                                 :: kfl_ipoin,kfl_ielem,kfl_const
    integer(ip)                                 :: kfl_grele,kfl_grpoi,kfl_drele,kfl_gdele
    integer(ip)                                 :: kfl_drele_tur,kfl_drele_vel

    kfl_ipoin     = 0
    kfl_ielem     = 0
    kfl_const     = 0
    kfl_grele     = 0
    kfl_grpoi     = 0
    kfl_drele     = 0
    kfl_gdele     = 0
    kfl_drele_tur = 0
    kfl_drele_vel = 0

    do imate = 1,nmate
       ilaws = prope_ker % ilaws(imate)
       if(      prope_ker % llaws(ilaws) % where == 'IPOIN' ) then
          kfl_ipoin = 1
          if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then
             kfl_grpoi = 1
          end if
       else if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then
          kfl_ielem = 1
          if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then
             kfl_grele = 1
          end if
          if( prope_ker % llaws(ilaws) % kfl_deriv == 1 ) then
             kfl_drele = 1
          end if
          if( prope_ker % llaws(ilaws) % kfl_deriv_tur == 1 ) then
             kfl_drele_tur = 1
          end if
          if( prope_ker % llaws(ilaws) % kfl_deriv_vel == 1 ) then
             kfl_drele_vel = 1
          end if
          if( prope_ker % llaws(ilaws) % kfl_grder == 1 ) then
             kfl_gdele = 1
          end if
       else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then
          kfl_const = 1
       end if
    end do

    ker_proper_property_array_exists = .false.

    select case ( trim(warray) )

    case ( 'VALUE_IPOIN'     )  ; if( kfl_ipoin == 1 .or.  kfl_ielem     == 1 ) ker_proper_property_array_exists = .true.
    case ( 'GRVAL_IPOIN'     )  ; if( kfl_ipoin == 1 .and. kfl_grpoi     == 1 ) ker_proper_property_array_exists = .true.
    case ( 'VALUE_IELEM'     )  ; if( kfl_ielem == 1                          ) ker_proper_property_array_exists = .true.
    case ( 'GRVAL_IELEM'     )  ; if( kfl_ielem == 1 .and. kfl_grele     == 1 ) ker_proper_property_array_exists = .true.
    case ( 'VALUE_IBOUN'     )  ; if( kfl_ielem == 1                          ) ker_proper_property_array_exists = .true.
    case ( 'VALUE_CONST'     )  ; if( kfl_const == 1                          ) ker_proper_property_array_exists = .true.
    case ( 'DRVAL_IELEM'     )  ; if( kfl_ielem == 1 .and. kfl_drele     == 1 ) ker_proper_property_array_exists = .true.
    case ( 'GDVAL_IELEM'     )  ; if( kfl_ielem == 1 .and. kfl_gdele     == 1 ) ker_proper_property_array_exists = .true.
    case ( 'DRVAL_TUR_IELEM' )  ; if( kfl_ielem == 1 .and. kfl_drele_tur == 1 ) ker_proper_property_array_exists = .true.
    case ( 'DRVAL_VEL_IELEM' )  ; if( kfl_ielem == 1 .and. kfl_drele_vel == 1 ) ker_proper_property_array_exists = .true.
    case default
       call runend('ker_proper_property_array_exists: UNKNOWN PROPERTY ARRAY')
    end select

  end function ker_proper_property_array_exists

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-17
  !> @brief   Point to a specific property
  !> @details Point to a specific property
  !> 
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-01-02
  !> @brief   Number of properties
  !> @details Number of existing properties 
  !> 
  !-----------------------------------------------------------------------
  
  function ker_proper_pointer_ip(iprop) result(pp)
    integer(ip),                 intent(in) :: iprop
    type(typ_valpr_ker), pointer            :: pp

    select case ( iprop )
       
    case (    1_ip ) ; pp => densi_ker 
    case (    2_ip ) ; pp => visco_ker 
    case (    3_ip ) ; pp => poros_ker 
    case (    4_ip ) ; pp => condu_ker 
    case (    5_ip ) ; pp => sphea_ker 
    case (    6_ip ) ; pp => dummy_ker 
    case (    7_ip ) ; pp => turmu_ker 
    case (    8_ip ) ; pp => absor_ker 
    case (    9_ip ) ; pp => scatt_ker 
    case (   10_ip ) ; pp => mixin_ker 
    case (   11_ip ) ; pp => anipo_ker 
    case (   12_ip ) ; pp => walvi_ker 
    case   default   ; call runend('KER_PROPER_LIST_PROPERTIES: UNKNOWN PROPERTY NUMBER')
       
    end select

  end function ker_proper_pointer_ip
  
  function ker_proper_pointer_w(wprop) result(pp)
    character(len=*),            intent(in) :: wprop
    type(typ_valpr_ker), pointer            :: pp

    select case ( wprop(1:5) )
       
    case ( 'DENSI' ) ; pp => densi_ker 
    case ( 'VISCO' ) ; pp => visco_ker 
    case ( 'POROS' ) ; pp => poros_ker 
    case ( 'CONDU' ) ; pp => condu_ker 
    case ( 'SPHEA' ) ; pp => sphea_ker 
    case ( 'DUMMY' ) ; pp => dummy_ker 
    case ( 'TURMU' ) ; pp => turmu_ker 
    case ( 'ABSOR' ) ; pp => absor_ker 
    case ( 'SCATT' ) ; pp => scatt_ker 
    case ( 'MIXIN' ) ; pp => mixin_ker 
    case ( 'ANIPO' ) ; pp => anipo_ker 
    case ( 'WALLV' ) ; pp => walvi_ker 

    case   default   ; call runend('KER_PROPER_LIST_PROPERTIES: UNKNOWN PROPERTY NUMBER')
       
    end select

  end function ker_proper_pointer_w

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-10-08
  !> @brief   Allocate memory
  !> @details Allocate memory for properties
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_proper_allocate_properties()

    integer(ip)                   :: iprop
    type(typ_valpr_ker), pointer  :: prope_ker

    !if( INOTMASTER ) then
    do iprop = 1,NUMBER_OF_PROPERTIES
       prope_ker => ker_proper_pointer(iprop)
       if( prope_ker % kfl_exist == 1 ) then
          call ker_wiprop(prope_ker)
          call ker_allpro(prope_ker)
          call ker_proper_smoothing(prope_ker)
       end if
    end do
    !end if

  end subroutine ker_proper_allocate_properties

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-10-08
  !> @brief   Allocate memory
  !> @details Allocate memory for material laws
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_proper_allocate_material_laws() 

    integer(ip)                   :: iprop
    type(typ_valpr_ker), pointer  :: prope_ker

    do iprop = 1,NUMBER_OF_PROPERTIES
       prope_ker => ker_proper_pointer(iprop)
       allocate( prope_ker % comp (nmate) )
       allocate( prope_ker % on_the_fly(nmate) )
       allocate( prope_ker % wlaws(nmate) )
       allocate( prope_ker % ilaws(nmate) )
       allocate( prope_ker % rlaws(mlapa_ker,nmate) )
       allocate( prope_ker % value_default(nmate) )
       allocate( prope_ker % update(2,nmate) )
       allocate( prope_ker % time_function(nmate) )
    end do

  end subroutine ker_proper_allocate_material_laws

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-10-08
  !> @brief   Initialization
  !> @details Give properties default values
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_proper_initialization() 

    integer(ip)                   :: iprop
    type(typ_valpr_ker), pointer  :: prope_ker

    cpu_prope = 0.0_rp
    
    do iprop = 1,NUMBER_OF_PROPERTIES
       prope_ker => ker_proper_pointer(iprop)
       nullify( prope_ker % comp  )
       nullify( prope_ker % on_the_fly )
       nullify( prope_ker % wlaws )
       nullify( prope_ker % ilaws )
       nullify( prope_ker % rlaws )
       nullify( prope_ker % value_default )
       nullify( prope_ker % update )
       nullify( prope_ker % time_function )

       nullify( prope_ker % value_ipoin     ) 
       nullify( prope_ker % grval_ipoin     )
       nullify( prope_ker % value_ielem     )
       nullify( prope_ker % grval_ielem     )
       nullify( prope_ker % value_iboun     )
       nullify( prope_ker % value_const     )
       nullify( prope_ker % drval_ielem     )
       nullify( prope_ker % gdval_ielem     )
       nullify( prope_ker % drval_tur_ielem )
       nullify( prope_ker % drval_vel_ielem )
    end do

  end subroutine ker_proper_initialization

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-10-08
  !> @brief   Initialization
  !> @details Give properties default values
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_proper_default_values() 

    integer(ip)                   :: iprop,ilaws
    type(typ_valpr_ker), pointer  :: prope_ker

    do iprop = 1,NUMBER_OF_PROPERTIES
       prope_ker => ker_proper_pointer(iprop)      
       !
       ! Material dependent 
       !
       prope_ker % comp          = 1
       prope_ker % on_the_fly    = 0
       prope_ker % wlaws         = 'CONST'
       prope_ker % ilaws         = 0
       prope_ker % rlaws         = 0.0_rp
       prope_ker % value_default = 0.0_rp
       prope_ker % kfl_exist     = 0
       prope_ker % kfl_nedsm     = 0
       prope_ker % kfl_type      = SCALAR_PROPERTY
       prope_ker % dim           = 1
       prope_ker % update(1,:)   = ITASK_DOITER       
       prope_ker % update(2,:)   = ITASK_ENDSTE
       prope_ker % time_function = 'NULL '
       !
       ! Law dependent
       !
       do ilaws = 1,mlaws_ker
          prope_ker % llaws(ilaws) % wname         = ''
          prope_ker % llaws(ilaws) % lresp         = -2
          prope_ker % llaws(ilaws) % where         = ''
          prope_ker % llaws(ilaws) % kfl_gradi     =  0
          prope_ker % llaws(ilaws) % kfl_deriv     =  0
          prope_ker % llaws(ilaws) % kfl_grder     =  0
          prope_ker % llaws(ilaws) % kfl_deriv_tur =  0
          prope_ker % llaws(ilaws) % kfl_deriv_vel =  0
       end do
    end do

  end subroutine ker_proper_default_values

  subroutine ker_allpro(prope_ker)
    !------------------------------------------------------------------------
    !****f* Master/ker_allpro
    ! NAME
    !    ker_allpro
    ! DESCRIPTION
    !    Allocate memory
    ! USES
    !    postpr
    !    memgen
    ! USED BY
    !    ker_output
    !***
    !------------------------------------------------------------------------

    type(typ_valpr_ker), intent(inout) :: prope_ker
    integer(ip)                        :: kfl_ipoin,kfl_ielem,kfl_const
    integer(ip)                        :: ielem,imate,ilaws,pgaus,pgaub
    integer(ip)                        :: iboun,kfl_grele,kfl_grpoi,kfl_drele,kfl_gdele,pelty,pnode
    integer(ip)                        :: kfl_drele_tur,kfl_drele_vel,ifact
    !
    ! Define some flags
    !
    kfl_ipoin = 0
    kfl_ielem = 0
    kfl_const = 0

    kfl_grele = 0
    kfl_grpoi = 0
    kfl_drele = 0
    kfl_drele_tur = 0
    kfl_drele_vel = 0
    kfl_gdele = 0

    do imate = 1,nmate
       ilaws = prope_ker % ilaws(imate)
       if(      prope_ker % llaws(ilaws) % where == 'IPOIN' ) then
          kfl_ipoin = 1
          if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then
             kfl_grpoi = 1
          end if
       else if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then
          kfl_ielem = 1
          if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then
             kfl_grele = 1
          end if
          if( prope_ker % llaws(ilaws) % kfl_deriv == 1 ) then
             kfl_drele = 1
          end if
          if( prope_ker % llaws(ilaws) % kfl_deriv_tur == 1 ) then
             kfl_drele_tur = 1
          end if
          if( prope_ker % llaws(ilaws) % kfl_deriv_vel == 1 ) then
             kfl_drele_vel = 1
          end if
          if( prope_ker % llaws(ilaws) % kfl_grder == 1 ) then
             kfl_gdele = 1
          end if
       else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then
          kfl_const = 1
       end if
    end do
    !
    ! Type of property
    !
    if( prope_ker % kfl_type == SCALAR_PROPERTY ) then
       ifact = 1
    else if( prope_ker % kfl_type == VECTOR_PROPERTY ) then
       ifact = ndime
    else if( prope_ker % kfl_type == MATRIX_PROPERTY ) then
       ifact = ndime * ndime
    end if
    prope_ker % dim = ifact
    !
    !
    ! Allocate memory
    !
    if ( (kfl_ipoin == 1 ) .or. ( kfl_ielem == 1 ) ) then
       !
       ! On nodes
       !
       ! even if the law is stored at ielem it might be smoothed if at some moment it is need at ipoin on npoin, so we allocate mem
       ! in order to save memory one could find a way to check if it will need to be smoothed -- needs to be thought
       !
       call memory_alloca(mem_modul(1:2,modul),'value_ipoin','ker_allpro',prope_ker % value_ipoin,npoin)
    end if
    if( ( kfl_ipoin == 1 ) .and. ( kfl_grpoi == 1 ) ) then 
       call memory_alloca(mem_modul(1:2,modul),'grval_ipoin','ker_allpro',prope_ker % grval_ipoin,ndime,npoin)
    end if

    if( kfl_ielem == 1 ) then
       !
       ! On elements
       !
       call memory_alloca(mem_modul(1:2,modul),'value_ielem','ker_allpro',prope_ker % value_ielem,nelem)
       do ielem = 1,nelem
          pgaus = ngaus(abs(ltype(ielem)))
          if( kfl_cutel /= 0 ) pgaus = mgaus
          call memory_alloca(mem_modul(1:2,modul),'VALUE_IELEM % A','ker_allpro',&
               prope_ker % value_ielem(ielem)%a,pgaus*ifact)
       end do
       if( kfl_grele == 1 ) then
          call memory_alloca(mem_modul(1:2,modul),'GRVAL_IELEM','ker_allpro',prope_ker % grval_ielem,nelem)
          do ielem = 1,nelem
             pgaus = ngaus(abs(ltype(ielem)))
             if( kfl_cutel /= 0 ) pgaus = mgaus
             call memory_alloca(mem_modul(1:2,modul),'GRVAL_IELEM % A','ker_allpro',&
                  prope_ker % grval_ielem(ielem)%a,ndime,pgaus)
          end do
       end if


       if( kfl_drele == 1 ) then
          call memory_alloca(mem_modul(1:2,modul),'DRVAL_IELEM','ker_allpro',prope_ker % drval_ielem,nelem)
          do ielem = 1,nelem
             pgaus = ngaus(abs(ltype(ielem)))
             pelty = ltype(ielem)
             pnode = nnode(pelty)
             call memory_alloca(mem_modul(1:2,modul),'DRVAL_IELEM % A','ker_allpro',&
                  prope_ker % drval_ielem(ielem)%a,pnode,pgaus)
          end do
       end if

       if( kfl_drele_tur == 1 ) then
          call memory_alloca(mem_modul(1:2,modul),'DRVAL_TUR_IELEM','ker_allpro',prope_ker % drval_tur_ielem,nelem)
          nturb = 1
          do ielem = 1,nelem
             pgaus = ngaus(abs(ltype(ielem)))
             pelty = ltype(ielem)
             pnode = nnode(pelty)
             call memory_alloca(mem_modul(1:2,modul),'DRVAL_TUR_IELEM % A','ker_allpro',&
                  prope_ker % drval_tur_ielem(ielem)%a,nturb,pnode,pgaus)
          end do
       end if

       if( kfl_drele_vel == 1 ) then
          call memory_alloca(mem_modul(1:2,modul),'DRVAL_VEL_IELEM','ker_allpro',prope_ker % drval_vel_ielem,nelem)
          do ielem = 1,nelem
             pgaus = ngaus(abs(ltype(ielem)))
             pelty = ltype(ielem)
             pnode = nnode(pelty)
             call memory_alloca(mem_modul(1:2,modul),'DRVAL_VEL_IELEM % A','ker_allpro',&
                  prope_ker % drval_vel_ielem(ielem)%a,ndime,pnode,pgaus)
          end do
       end if

       if( kfl_gdele == 1 ) then
          call memory_alloca(mem_modul(1:2,modul),'GDVAL_IELEM','ker_allpro',prope_ker % gdval_ielem,nelem)
          do ielem = 1,nelem
             pgaus = ngaus(abs(ltype(ielem)))
             pelty = ltype(ielem)
             pnode = nnode(pelty)
             call memory_alloca(mem_modul(1:2,modul),'GDVAL_IELEM % A','ker_allpro',&
                  prope_ker % gdval_ielem(ielem)%a,ndime,pnode,pgaus)
          end do
       end if

       !
       ! On boundaries
       !
       call memory_alloca(mem_modul(1:2,modul),'VALUE_IBOUN','ker_allpro',prope_ker % value_iboun,nboun)
       do iboun = 1,nboun
          pgaub = ngaus(abs(ltypb(iboun)))
          call memory_alloca(mem_modul(1:2,modul),'VALUE_IBOUN % A','ker_allpro',&
               prope_ker % value_iboun(iboun)%a,pgaub)
       end do

    end if
    !
    ! Constant properties
    !
    if( kfl_const == 1 ) then
       call memory_alloca(mem_modul(1:2,modul),'VALUE_CONST','ker_allpro',prope_ker % value_const,ifact,nmate)
    end if

  end subroutine ker_allpro

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-18
  !> @brief   Name to number of property
  !> @details Correspondance number vs name of properties
  !> 
  !-----------------------------------------------------------------------

  integer(ip) function ker_proper_law_name_to_number(wprop,wlaws)

    character(len=*),            intent(in) :: wprop
    character(len=*),            intent(in) :: wlaws
    type(typ_valpr_ker), pointer            :: prope_ker
    integer(ip)                             :: imate,ilaws

    prope_ker => ker_proper_pointer(wprop)

    do imate = 1,nmate
       ilaws = 1
       do while( ilaws <=  mlaws_ker )
          if( prope_ker % llaws(ilaws) % wname == wlaws ) then
             ker_proper_law_name_to_number = ilaws
             ilaws = mlaws_ker + 1
          end if
          ilaws = ilaws + 1
       end do
       if( ilaws /= mlaws_ker + 2 ) then
          write(*,*) 'KER_PROPER_LAW_NAME_TO_NUMBER: PROPERTY LAW DOES NOT EXIST:imate,trim(prope_ker % wlaws(imate))',&
               imate,trim(prope_ker % wlaws(imate))
          ! this can happen if in your geo file you have more materias than those you define in ker.dat
          call runend('KER_PROPER_LAW_NAME_TO_NUMBER: PROPERTY LAW DOES NOT EXIST')
       end if
    end do

  end function ker_proper_law_name_to_number

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-18
  !> @brief   Name to number of property
  !> @details Correspondance number vs name of properties
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_wiprop(prope_ker)

    type(typ_valpr_ker), intent(inout) :: prope_ker
    integer(ip)                        :: imate,ilaws

    do imate = 1,nmate
       ilaws = 1
       do while( ilaws <=  mlaws_ker )
          if( prope_ker % llaws(ilaws) % wname == prope_ker % wlaws(imate) ) then
             prope_ker % ilaws(imate) = ilaws
             ilaws = mlaws_ker + 1
          end if
          ilaws = ilaws + 1
       end do
       if( ilaws /= mlaws_ker + 2 ) then
          write(*,*) 'KER_WIPROP: PROPERTY LAW DOES NOT EXIST:imate,trim(prope_ker % wlaws(imate))',&
               imate,trim(prope_ker % wlaws(imate))
          ! this can happen if in your geo file you have more materias than those you define in ker.dat
          call runend('KER_WIPROP: PROPERTY LAW DOES NOT EXIST')
       end if
    end do

  end subroutine ker_wiprop

  subroutine ker_smopro(prope_ker,xvalu)
    !------------------------------------------------------------------------
    !****f* Master/ker_smopro
    ! NAME
    !    ker_smopro
    ! DESCRIPTION
    !    Smooth a property
    ! USES
    !    postpr
    !    memgen
    ! USED BY
    !    ker_output
    !***
    !------------------------------------------------------------------------
   
    type(typ_valpr_ker), pointer      :: prope_ker
    real(rp),     intent(out), target :: xvalu(*)
    integer(ip)                       :: ipoin,idime,inode,ielem,igaus
    integer(ip)                       :: pnode,pelty,pgaus,ilaws,imate
    integer(ip)                       :: jnode,knode
    real(rp)                          :: detjm,gpvol,gpcar(ndime,mnode)
    real(rp)                          :: elcod(ndime,mnode)
    real(rp)                          :: xjaci(9),xjacm(9),xfact
    integer(ip),  pointer             :: lmele(:)
    real(rp),     pointer             :: vmass_ker(:)

    if( INOTMASTER ) then
       !
       ! Smooth property
       !
       nullify( lmele )
       nullify( vmass_ker )
       call memory_alloca(mem_modul(1:2,modul),'LMELE'    ,'ker_smopro',lmele,nelem)
       call memory_alloca(mem_modul(1:2,modul),'VMASS_KER','ker_smopro',vmass_ker,npoin)
       do imate = 1,nmate
          ilaws = prope_ker % ilaws(imate)
          if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then
             do ielem = 1,nelem
                if( lmate(ielem) == imate ) lmele(ielem) = 1
             end do
          end if
       end do
       !
       ! Initialization
       !
       do ipoin = 1,npoin
          vmass_ker(ipoin) = 0.0_rp
          xvalu(ipoin)     = 0.0_rp
       end do
       !
       ! Loop over elements
       !
       elements: do ielem = 1,nelem
          pelty = ltype(ielem)

          if( pelty > 0 .and. lmele(ielem) == 1 ) then
             pnode = nnode(pelty)
             pgaus = ngaus(pelty)
             !
             ! Gather vectors
             !
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                do idime = 1,ndime
                   elcod(idime,inode) = coord(idime,ipoin)
                end do
             end do
             !
             ! Loop over Gauss points
             !
             gauss_points: do igaus = 1,pgaus
                call elmder(&
                     pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                     elcod,gpcar,detjm,xjacm,xjaci)
                gpvol = elmar(pelty) % weigp(igaus) * detjm
                if( kfl_naxis == 1 ) then
                   call runend('MOD_GRADIE: NOT CODED')
                end if
                !
                ! Extension
                !
                if( lelch(ielem) == ELEXT ) then
                   knode = 1
                else
                   knode = pnode
                end if
                !
                ! Assemble
                !
                do inode = 1,knode
                   ipoin        = lnods(inode,ielem)
                   xfact        = gpvol * elmar(pelty) % shape(inode,igaus)
                   xvalu(ipoin) = xvalu(ipoin) + xfact * prope_ker % value_ielem(ielem) % a(igaus)
                   do jnode = 1,pnode
                      vmass_ker(ipoin) = vmass_ker(ipoin) + xfact * elmar(pelty) % shape(jnode,igaus)
                   end do
                end do

             end do gauss_points
          end if
       end do elements
       !
       ! Parallelization
       !
       call rhsmod(1_ip,xvalu)
       call rhsmod(1_ip,vmass_ker)
       !memory_alloca
       ! Solve diagonal system
       !
       do ipoin = 1,npoin
          if( vmass_ker(ipoin) /= 0.0_rp ) &
               xvalu(ipoin) = xvalu(ipoin) / vmass_ker(ipoin)
       end do
       !
       ! Deallocate memory
       !
       call memory_deallo(mem_modul(1:2,modul),'VMASS_KER','ker_smopro',vmass_ker)
       call memory_deallo(mem_modul(1:2,modul),'LMELE'    ,'ker_smopro',lmele)

    end if

  end subroutine ker_smopro

  subroutine ker_grapro(prope_ker,xvalu)
    !------------------------------------------------------------------------
    !****f* Master/ker_grapro
    ! NAME
    !    ker_grapro
    ! DESCRIPTION
    !    Smooth the gradient of a property
    ! USES
    !    postpr
    !    memgen
    ! USED BY
    !    ker_output
    !***
    !------------------------------------------------------------------------

    type(typ_valpr_ker), pointer      :: prope_ker
    real(rp),     intent(out), target :: xvalu(ndime,npoin)
    integer(ip)                       :: ipoin,idime,inode,ielem,igaus
    integer(ip)                       :: pnode,pelty,pgaus,ilaws,imate
    integer(ip)                       :: jnode,knode
    real(rp)                          :: detjm,gpvol,gpcar(ndime,mnode)
    real(rp)                          :: elcod(ndime,mnode)
    real(rp)                          :: grunk(ndime)
    real(rp)                          :: xjaci(9),xjacm(9),xfact
    integer(ip),  pointer             :: lmele(:)
    real(rp),     pointer             :: vmass_ker(:)

    if( INOTMASTER ) then
       !
       ! Smooth property
       !
       nullify( lmele )
       nullify( vmass_ker )
       call memory_alloca(mem_modul(1:2,modul),'LMELE'    ,'ker_grapro',lmele,nelem)
       call memory_alloca(mem_modul(1:2,modul),'VMASS_KER','ker_grapro',vmass_ker,npoin)
       do imate = 1,nmate
          ilaws = prope_ker % ilaws(imate)
          if(    ( prope_ker % llaws(ilaws) % where == 'IELEM' .and.&
               &   prope_ker % llaws(ilaws) % kfl_gradi == 1 ) .or. &
               & ( prope_ker % llaws(ilaws) % where == 'IPOIN' .and.&
               &   prope_ker % llaws(ilaws) % kfl_gradi == 0 ) ) then
             do ielem = 1,nelem
                if( lmate(ielem) == imate ) lmele(ielem) = ilaws
             end do
          end if
       end do
       !
       ! Initialization
       !
       do ipoin = 1,npoin
          vmass_ker(ipoin) = 0.0_rp
       end do
       do ipoin = 1,npoin
          do idime = 1,ndime
             xvalu(idime,ipoin) = 0.0_rp
          end do
       end do
       !
       ! Loop over elements
       !
       elements: do ielem = 1,nelem
          pelty = ltype(ielem)

          if( pelty > 0 .and. lmele(ielem) /= 0 ) then
             pnode = nnode(pelty)
             pgaus = ngaus(pelty)
             ilaws = lmele(ielem)
             !
             ! Gather vectors
             !
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                do idime = 1,ndime
                   elcod(idime,inode) = coord(idime,ipoin)
                end do
             end do
             !
             ! Loop over Gauss points
             !
             gauss_points: do igaus = 1,pgaus
                call elmder(&
                     pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                     elcod,gpcar,detjm,xjacm,xjaci)
                gpvol = elmar(pelty) % weigp(igaus) * detjm

                if( kfl_naxis == 1 ) then
                   call runend('MOD_GRADIE: NOT CODED')
                end if
                !
                ! Extension
                !
                if( lelch(ielem) == ELEXT ) then
                   knode = 1
                else
                   knode = pnode
                end if
                !
                ! Assemble
                !
                if(        prope_ker % llaws(ilaws) % where == 'IELEM' &
                     .and. prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then

                   do inode = 1,knode
                      ipoin = lnods(inode,ielem)
                      xfact = gpvol * elmar(pelty) % shape(inode,igaus)
                      do idime = 1,ndime
                         xvalu(idime,ipoin) = xvalu(idime,ipoin) &
                              + xfact * prope_ker % grval_ielem(ielem) % a(idime,igaus)
                      end do
                      do jnode = 1,knode
                         vmass_ker(ipoin) = vmass_ker(ipoin) + xfact * elmar(pelty) % shape(jnode,igaus)
                      end do
                   end do

                else if(   prope_ker % llaws(ilaws) % where == 'IPOIN' &
                     .and. prope_ker % llaws(ilaws) % kfl_gradi == 0 ) then

                   grunk = 0.0_rp
                   do inode = 1,knode
                      ipoin = lnods(inode,ielem)
                      do idime = 1,ndime
                         grunk(idime) = grunk(idime) + prope_ker % value_ipoin(ipoin) * gpcar(idime,inode)
                      end do
                   end do
                   do inode = 1,knode
                      ipoin = lnods(inode,ielem)
                      xfact = gpvol * elmar(pelty) % shape(inode,igaus)
                      do idime = 1,ndime
                         xvalu(idime,ipoin) = xvalu(idime,ipoin) + xfact * grunk(idime)
                      end do
                      do jnode = 1,pnode
                         vmass_ker(ipoin) = vmass_ker(ipoin) + xfact * elmar(pelty) % shape(jnode,igaus)
                      end do
                   end do

                end if

             end do gauss_points
          end if
       end do elements
       !
       ! Parallelization
       !
       call rhsmod(ndime,xvalu)
       call rhsmod( 1_ip,vmass_ker)
       !
       ! Solve diagonal system
       !
       do ipoin = 1,npoin
          if( vmass_ker(ipoin) /= 0.0_rp ) then
             do idime = 1,ndime
                xvalu(idime,ipoin) = xvalu(idime,ipoin) / vmass_ker(ipoin)
             end do
          end if
       end do
       !
       ! Deallocate memory
       !
       call memory_deallo(mem_modul(1:2,modul),'VMASS_KER','ker_grapro',vmass_ker)
       call memory_deallo(mem_modul(1:2,modul),'LMELE'    ,'ker_grapro',lmele)

    end if

  end subroutine ker_grapro

  subroutine ker_proper_scalar(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor)
    !-----------------------------------------------------------------------
    !****f* Kermod/ker_proper
    ! NAME
    !    ker_proper
    ! DESCRIPTION
    !    Update properties
    !
    !    Property can be defined with 3 formats
    !    It can be required at 8 different places
    !
    !       WHEREIN   IPOSI  IELBO
    !       --------------------
    !    1. IPOIN   IPOIN  NONSENSE
    !    2. NPOIN   DUMMI  Do smoothing => Loop over IELEM
    !    3. IGAUS   IGAUS  IELEM
    !    4. PGAUS   DUMMI  IELEM
    !    5. IGAUB   IGAUB  IBOUN        => IELEM = LBOEL(,IBOUN)
    !    6. PGAUB   DUMMI  IBOUN        => IELEM = LBOEL(,IBOUN)
    !    7. COG     DUMMI  IELEM
    !    8. PNODE   DUMMI  IELEM
    !
    !    GPSHA is dimensioned as a vector in order to use it as optional
    !    argument
    !
    ! USED BY
    !    many subroutines
    !***
    !-----------------------------------------------------------------------

    use def_kintyp,      only                 :  ip,rp
    use def_elmtyp,      only                 :  ELCUT
    use def_master,      only                 :  mem_modul,modul,nturb
    use def_domain,      only                 :  nmate,npoin,mnode,lelch
    use def_domain,      only                 :  lnods,lnnod,lmate,lmatn,nmatn
    use def_domain,      only                 :  ngaus,ltype,elmar,ndime
    use def_domain,      only                 :  lelbo,ltypb,lnodb,nnode
    use def_kermod,      only                 :  typ_valpr_ker
    use def_kermod,      only                 :  densi_ker,visco_ker
    use def_kermod,      only                 :  poros_ker,condu_ker
    use def_kermod,      only                 :  sphea_ker,dummy_ker
    use def_kermod,      only                 :  absor_ker,scatt_ker
    use def_kermod,      only                 :  turmu_ker,mixin_ker
    use def_kermod,      only                 :  anipo_ker,walvi_ker
    use mod_memory
    implicit none
    character(5),        intent(in)           :: wname,wherein
    integer(ip),         intent(in)           :: iposi,ielbo
    real(rp),            intent(out)          :: xvalu(*)
    integer(ip),         intent(in), optional :: qnode
    integer(ip),         intent(in), optional :: qgaus
    real(rp),            intent(in), optional :: gpsha(*)
    real(rp),            intent(in), optional :: gpcar(ndime,mnode,*)
    real(rp),            intent(in), optional :: gpcor(ndime,*)
    type(typ_valpr_ker), pointer              :: prope_ker
    integer(ip)                               :: imate,ielem,inode,ipoin,igaus,kauxi,iturb
    integer(ip)                               :: kpoin,ilaws,pelty,pgaus,pnode,jpoin,kk
    integer(ip)                               :: pnodb,igaub,inodb,pblty,iboun
    integer(ip)                               :: pgaub,kposi,kfl_gradi,idime,kfl_deriv,kfl_grder
    integer(ip)                               :: kfl_deriv_tur,kfl_deriv_vel
    integer(ip)                               :: lposi
    real(rp)                                  :: elpro(mnode)
    real(rp)                                  :: elgrp(ndime,mnoga)
    real(rp),   pointer                       :: auxva(:)

    nullify(auxva)

    !----------------------------------------------------------------------
    !
    ! Some definitions
    !
    !----------------------------------------------------------------------

    kfl_gradi     = 0
    kfl_deriv     = 0
    kfl_grder     = 0
    kfl_deriv_tur = 0
    kfl_deriv_vel = 0

    if(      wname(1:2) == 'GR'  ) then
       kfl_gradi = 1
    else if( wname(1:2) == 'DR'  ) then
       kfl_deriv = 1
    else if( wname(1:2) == 'GD'  ) then
       kfl_grder = 1
    else if( wname(1:3) == 'TDR' ) then
       kfl_deriv_tur = 1
    else if( wname(1:3) == 'VDR' ) then
       kfl_deriv_vel = 1
    end if

    if(      wname == 'DENSI' .or. wname == 'GRDEN' .or. wname == 'DRDEN' .or. wname == 'GDDEN' ) then
       prope_ker => densi_ker
    else if( wname == 'VISCO' .or. wname == 'GRVIS' .or. wname == 'DRVIS' .or. wname == 'GDVIS' ) then
       prope_ker => visco_ker
    else if( wname == 'MIXIN' ) then
       prope_ker => mixin_ker
    else if( wname == 'POROS' .or. wname == 'GRPOS' ) then
       prope_ker => poros_ker
    else if( wname == 'CONDU' .or. wname == 'GRCON' ) then
       prope_ker => condu_ker
    else if( wname == 'SPHEA' .or. wname == 'GRSPE' ) then
       prope_ker => sphea_ker
    else if( wname == 'DUMMY' .or. wname == 'GRDUM' ) then
       prope_ker => dummy_ker
    else if( wname == 'TURBU' .or. wname == 'TURVI' .or. wname == 'GRTUR' .or. wname == 'TDRTU' .or. wname == 'VDRTU' ) then
       prope_ker => turmu_ker
    else if( wname == 'ABSOR' ) then
       prope_ker => absor_ker
    else if( wname == 'SCATT' ) then
       prope_ker => scatt_ker
    else if( wname == 'ANIPO' ) then
       prope_ker => anipo_ker
    else if( wname == 'WALLV' ) then
       prope_ker => walvi_ker
    else
       print *,'--- I DONT UNDERSTAND ',wname
       call runend('KER_PROPER CALLED WITH WRONG PROPERTY') ! This avoids the bug of calling wrongly the ker_proper routine
    end if

    if( prope_ker % kfl_exist == 0 ) then
       !
       ! Put value to zero
       !
       if( wherein == 'IPOIN' .or. wherein == 'IGAUS' .or. wherein == 'IGAUB' .or. wherein == 'COG' ) then
          kposi = 1
       else if( wherein == 'NPOIN' ) then
          kposi = npoin
       else if( wherein == 'PNODE' ) then
          kposi = lnnod(ielbo)
       else if( wherein == 'PGAUS' ) then
          if( present(gpsha) ) then
             if( present(qgaus) ) then
                kposi = qgaus
             else
                call runend('MOD_KER_PROPER: ARGUMENT IS MISSING')
             end if
          else
             kposi = ngaus(abs(ltype(ielbo)))
          end if
       else if( wherein == 'PGAUB' ) then
          if( present(gpsha) ) then
             if( present(qgaus) ) then
                kposi = qgaus
             else
                call runend('MOD_KER_PROPER: ARGUMENT IS MISSING')
             end if
          else
             kposi = ngaus(abs(ltypb(ielbo)))
          end if
       end if
       if( kfl_gradi == 1 ) kposi = kposi * ndime
       do lposi = 1,kposi
          xvalu(lposi) = 0.0_rp
       end do

    else

       !----------------------------------------------------------------------
       !
       ! 1. Property required on IPOIN
       !
       !----------------------------------------------------------------------

       if( wherein == 'IPOIN' ) then

          if( kfl_gradi == 1 ) then
             call runend('KER_PROPER: DONT KNOW WHAT TO DO WITH THIS OPTION')
          end if

          ipoin = iposi
          !
          ! Smooth value - If in some material it is stored by ielem it is smoothed
          ! In the materials where it is not stored by ielem it will step over with erroneous values
          ! Therefore if is then recalculated in materials where it is stored by IPOIN or CONST
          !
          kauxi = 0
          if ( prope_ker % kfl_nedsm == 1 ) then
             imat2_loop: do imate = 1,nmate
                ilaws = prope_ker % ilaws(imate)
                if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then
                   call memory_alloca(mem_modul(1:2,modul),'auxva','ker_proper',auxva,npoin)
                   call ker_smopro(prope_ker,auxva)
                   prope_ker % kfl_nedsm = 0
                   kauxi = 1
                    exit imat2_loop
                end if
             end do imat2_loop
          end if

          if (kauxi == 1) then  ! correct the value (auxva) in points where it is calculated by ipoin or const
             ! and has been steped over when smoothing due to materials whith law at ielem
             do imate = 1,nmate
                ilaws = prope_ker % ilaws(imate)

                if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

                   do kpoin = 1,nmatn(imate)
                      jpoin = lmatn(imate) % l(kpoin)
                      auxva(jpoin) = prope_ker % value_ipoin(jpoin)
                   end do

                else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

                   do kpoin = 1,nmatn(imate)
                      jpoin = lmatn(imate) % l(kpoin)
                      auxva(jpoin) = prope_ker % value_const(1,imate)
                   end do

                end if
             end do
             !
             ! Save values obtained in auxva to prope_ker % value_ipoin(jpoin)
             !
             do kpoin =1,npoin
                prope_ker % value_ipoin(kpoin) = auxva(kpoin)
             end do
             call memory_deallo(mem_modul(1:2,modul),'auxva','ker_proper',auxva)
          end if


          do imate = 1,nmate
             ilaws = prope_ker % ilaws(imate)

             if( ( prope_ker % llaws(ilaws) % where == 'IPOIN' ) .or. ( prope_ker % llaws(ilaws) % where == 'IELEM' ) )then

                xvalu(1) = prope_ker % value_ipoin(ipoin)
                
             else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

                xvalu(1) = prope_ker % value_const(1,imate)

             end if
          end do
       end if

       !----------------------------------------------------------------------
       !
       ! 2. Property required on NPOIN
       !
       !----------------------------------------------------------------------

       if( wherein == 'NPOIN' ) then

          if( kfl_gradi == 1 ) then
             !
             ! Gradient
             !
             imat0_loop: do imate = 1,nmate
                ilaws = prope_ker % ilaws(imate)
                if(    ( prope_ker % llaws(ilaws) % where == 'IELEM' .and.&
                     &   prope_ker % llaws(ilaws) % kfl_gradi == 1 ) .or. &
                     & ( prope_ker % llaws(ilaws) % where == 'IPOIN' .and.&
                     &   prope_ker % llaws(ilaws) % kfl_gradi == 0 ) ) then
                   call ker_grapro(prope_ker,xvalu)
                   exit imat0_loop
                end if
             end do imat0_loop

             do imate = 1,nmate
                ilaws = prope_ker % ilaws(imate)

                if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

                   if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then

                      do kpoin = 1,nmatn(imate)
                         ipoin = lmatn(imate) % l(kpoin)
                         kposi = (ipoin-1) * ndime
                         do idime = 1,ndime
                            kposi = kposi + 1
                            xvalu(kposi) = prope_ker % grval_ipoin(idime,ipoin)
                         end do
                      end do

                   end if

                else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      kposi = (ipoin-1) * ndime
                      do idime = 1,ndime
                         kposi = kposi + 1
                         xvalu(kposi) = 0.0_rp
                      end do
                   end do

                end if
             end do

          else
             !
             ! Value - If in some material it is stored by ielem it is smoothed
             ! In the materials where it is not stored by ielem it will step over with erroneous values
             ! Therefore if is then recalculated in materials where it is stored by IPOIN or CONST
             !
             kauxi = 0
             if ( prope_ker % kfl_nedsm == 1 ) then
                imat1_loop: do imate = 1,nmate
                   ilaws = prope_ker % ilaws(imate)
                   if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then
                      call memory_alloca(mem_modul(1:2,modul),'auxva','ker_proper',auxva,npoin)
                      call ker_smopro(prope_ker,auxva)
                      prope_ker % kfl_nedsm = 0
                      kauxi = 1
                      exit imat1_loop
                   end if
                end do imat1_loop
             end if

             if (kauxi == 1) then  ! correct the value (auxva) in points where it is calculated by ipoin or const
                ! and has been steped over when smoothing due to materials whith law at ielem
                do imate = 1,nmate
                   ilaws = prope_ker % ilaws(imate)

                   if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

                      do kpoin = 1,nmatn(imate)
                         jpoin = lmatn(imate) % l(kpoin)
                         auxva(jpoin) = prope_ker % value_ipoin(jpoin)
                      end do

                   else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

                      do kpoin = 1,nmatn(imate)
                         jpoin = lmatn(imate) % l(kpoin)
                         auxva(jpoin) = prope_ker % value_const(1,imate)
                      end do

                   end if
                end do
                !
                ! Save values obtained in auxva to prope_ker % value_ipoin(jpoin)
                !
                do kpoin =1,npoin
                   prope_ker % value_ipoin(kpoin) = auxva(kpoin)
                end do
                call memory_deallo(mem_modul(1:2,modul),'auxva','ker_proper',auxva)
             end if

             do imate = 1,nmate
                ilaws = prope_ker % ilaws(imate)

                if( ( prope_ker % llaws(ilaws) % where == 'IPOIN' ) .or. ( prope_ker % llaws(ilaws) % where == 'IELEM' ) )then

                   do kpoin = 1,nmatn(imate)
                      jpoin = lmatn(imate) % l(kpoin)
                      xvalu(jpoin) = prope_ker % value_ipoin(jpoin)
                   end do

                else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

                   do kpoin = 1,nmatn(imate)
                      jpoin = lmatn(imate) % l(kpoin)
                      xvalu(jpoin) = prope_ker % value_const(1,imate)
                   end do

                end if
             end do
          end if

       end if

       !----------------------------------------------------------------------
       !
       ! 3. Property required on IGAUS
       !
       !----------------------------------------------------------------------

       if( wherein == 'IGAUS' ) then

          igaus = iposi
          ielem = ielbo
          imate = lmate(ielem)
          ilaws = prope_ker % ilaws(imate)

          if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then

             if( kfl_gradi == 1 ) then
                !
                ! Gradient
                !
                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then
                   do idime = 1,ndime
                      xvalu(idime) = prope_ker % grval_ielem(ielem) % a(idime,igaus)
                   end do
                else
                   xvalu(1:ndime) = 0.0_rp
                end if
             else
                !
                ! Value
                !
                if( present(gpsha) ) then
                   pelty = abs(ltype(ielem))
                   pnode = lnnod(ielem)
                   pgaus = ngaus(pelty)
                   xvalu(1) = 0.0_rp
                   do inode = 1,pnode
                      do igaus = 1,pgaus
                         xvalu(1) = xvalu(1) + gpsha(inode) * elmar(pelty) % shaga(igaus,inode) * prope_ker % value_ielem(ielem) % a(igaus)
                      end do
                   end do
                else
                   xvalu(1) = prope_ker % value_ielem(ielem) % a(igaus)
                end if
             end if

          else if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

             pelty = abs(ltype(ielem))
             pnode = lnnod(ielem)
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                elpro(inode) = prope_ker % value_ipoin(ipoin)
             end do

             if( kfl_gradi == 1 ) then
                !
                ! Gradient
                !
                xvalu(1:ndime) = 0.0_rp

                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then

                   do inode = 1,pnode
                      ipoin = lnods(inode,ielem)
                      do idime = 1,ndime
                         elgrp(idime,inode) = prope_ker % grval_ipoin(idime,ipoin)
                      end do
                   end do
                   if( present(gpsha) ) then
                      kposi = (igaus-1) * qnode
                      do inode = 1,pnode
                         kposi = kposi + 1
                         do idime = 1,ndime
                            xvalu(idime) = xvalu(idime) + elgrp(idime,inode) * gpsha(kposi)
                         end do
                      end do
                   else
                      do inode = 1,pnode
                         do idime = 1,ndime
                            xvalu(idime) = xvalu(idime) + elgrp(idime,inode) * elmar(pelty) % shape(idime,igaus)
                         end do
                      end do
                   end if

                else

                   if( present(gpcar) ) then
                      do inode = 1,pnode
                         do idime = 1,ndime
                            xvalu(idime) = xvalu(idime) + elpro(inode) * gpcar(idime,inode,1)
                         end do
                      end do
                   else
                      call runend('MOD_KER_PROPER: GPCAR NEEDED')
                   end if
                end if

             else
                !
                ! Value
                !
                xvalu(1) = 0.0_rp
                if( present(gpsha) ) then
                   kposi    = (igaus-1) * qnode
                   do inode = 1,pnode
                      kposi    = kposi + 1
                      xvalu(1) = xvalu(1) + elpro(inode) * gpsha(kposi)
                   end do
                else
                   do inode = 1,pnode
                      xvalu(1) = xvalu(1) + elpro(inode) * elmar(pelty) % shape(inode,igaus)
                   end do
                end if
             end if

          else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

             if( kfl_gradi == 1 ) then
                xvalu(1:ndime) = 0.0_rp
             else
                xvalu(1) = prope_ker % value_const(1,imate)
             end if
          end if

       end if

       !----------------------------------------------------------------------
       !
       ! 4. Property required on IGAUS=1,PGAUS
       !
       !----------------------------------------------------------------------

       if( wherein == 'PGAUS' ) then

          ielem = ielbo
          imate = lmate(ielem)
          ilaws = prope_ker % ilaws(imate)
          pelty = abs(ltype(ielem))
          pgaus = ngaus(pelty)
          if( lelch(ielem) == ELCUT ) pgaus = lgaus(ielem)

          if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then

             if( kfl_gradi == 1 ) then

                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then

                   if( present(gpsha) ) then
                      if( qgaus /= pgaus .and. kfl_cutel == 0 ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                      do igaus = 1,qgaus
                         kposi = (igaus - 1) * ndime
                         do idime = 1,ndime
                            kposi = kposi + 1
                            xvalu(kposi) = prope_ker % grval_ielem(ielem) % a(idime,igaus)
                         end do
                      end do
                   else
                      if( qgaus /= pgaus .and. kfl_cutel == 0 ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                      do igaus = 1,qgaus
                         kposi = (igaus - 1) * ndime
                         do idime = 1,ndime
                            kposi = kposi + 1
                            xvalu(kposi) = prope_ker % grval_ielem(ielem) % a(idime,igaus)
                         end do
                      end do
                   end if
                else
                   do kposi = 1,pgaus*ndime
                      xvalu(kposi) = 0.0_rp
                   end do
                end if

             else if( kfl_deriv == 1 ) then

                if( prope_ker % llaws(ilaws) % kfl_deriv == 1 ) then

                   if( present(gpsha) ) then
                      if( qgaus /= pgaus ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                      do igaus = 1,qgaus
                         kposi = (igaus - 1) * qnode
                         do inode = 1,qnode
                            kposi = kposi + 1
                            xvalu(kposi) = prope_ker % drval_ielem(ielem) % a(inode,igaus)
                         end do
                      end do
                   else
                      do igaus = 1,pgaus
                         kposi = (igaus - 1) * qnode
                         do inode = 1,qnode
                            kposi = kposi + 1
                            xvalu(kposi) = prope_ker % drval_ielem(ielem) % a(inode,igaus)
                         end do
                      end do
                   end if
                else
                   do kposi = 1,pgaus*qnode
                      xvalu(kposi) = 0.0_rp
                   end do
                end if

             else if( kfl_grder == 1 ) then

                if( prope_ker % llaws(ilaws) % kfl_grder == 1 ) then

                   if( present(gpsha) ) then
                      if( qgaus /= pgaus ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                      do igaus = 1,qgaus
                         kposi = (igaus - 1)* qnode*ndime
                         do inode = 1,qnode
                            do idime = 1, ndime
                               kposi = kposi + 1
                               xvalu(kposi) = prope_ker % gdval_ielem(ielem) % a(idime,inode,igaus)
                            enddo
                         end do
                      end do
                   else
                      do igaus = 1,pgaus
                         do inode = 1,qnode
                            kposi = (igaus - 1)*(inode-1) * qnode*ndime
                            do idime = 1, ndime
                               kposi = kposi + 1
                               xvalu(kposi) = prope_ker % gdval_ielem(ielem) % a(idime,inode,igaus)
                            enddo
                         end do
                      end do
                   end if
                else
                   do kposi = 1,pgaus*qnode*ndime
                      xvalu(kposi) = 0.0_rp
                   end do
                end if

             else if( kfl_deriv_tur == 1 ) then

                nturb = 1
                if( prope_ker % llaws(ilaws) % kfl_deriv_tur == 1 ) then

                   if( present(gpsha) ) then
                      if( qgaus /= pgaus ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                      do igaus = 1,qgaus
                         kposi = (igaus - 1)* qnode*nturb
                         do inode = 1,qnode
                            do iturb = 1, nturb
                               kposi = kposi + 1
                               xvalu(kposi) = prope_ker % drval_tur_ielem(ielem) % a(iturb,inode,igaus)
                            enddo
                         end do
                      end do
                   else
                      do igaus = 1,pgaus
                         do inode = 1,qnode
                            kposi = (igaus - 1)*(inode-1) * qnode*nturb
                            do iturb = 1, nturb
                               kposi = kposi + 1
                               xvalu(kposi) = prope_ker % drval_tur_ielem(ielem) % a(iturb,inode,igaus)
                            enddo
                         end do
                      end do
                   end if
                else
                   do kposi = 1,pgaus*qnode*nturb
                      xvalu(kposi) = 0.0_rp
                   end do
                end if

             else if( kfl_deriv_vel == 1 ) then

                if( prope_ker % llaws(ilaws) % kfl_deriv_vel == 1 ) then

                   if( present(gpsha) ) then
                      if( qgaus /= pgaus ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                      do igaus = 1,qgaus
                         kposi = (igaus - 1)* qnode*ndime
                         do inode = 1,qnode
                            do idime = 1, ndime
                               kposi = kposi + 1
                               xvalu(kposi) = prope_ker % drval_vel_ielem(ielem) % a(idime,inode,igaus)
                            enddo
                         end do
                      end do
                   else
                      do igaus = 1,pgaus
                         do inode = 1,qnode
                            kposi = (igaus - 1)*(inode-1) * qnode*ndime
                            do idime = 1, ndime
                               kposi = kposi + 1
                               xvalu(kposi) = prope_ker % drval_vel_ielem(ielem) % a(idime,inode,igaus)
                            enddo
                         end do
                      end do
                   end if
                else
                   do kposi = 1,pgaus*qnode*ndime
                      xvalu(kposi) = 0.0_rp
                   end do
                end if

             else

                if( present(gpsha) .and. kfl_cutel == 0 ) then
                   if( qgaus /= pgaus ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                   do igaus = 1,qgaus
                      xvalu(igaus) = prope_ker % value_ielem(ielem) % a(igaus)
                   end do
                else
                   do igaus = 1,pgaus
                      xvalu(igaus) = prope_ker % value_ielem(ielem) % a(igaus)
                   end do
                end if
             end if

          else if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

             pelty = abs(ltype(ielem))
             pnode = lnnod(ielem)
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                elpro(inode) = prope_ker % value_ipoin(ipoin)
             end do

             if( kfl_gradi == 1 ) then
                !
                ! Gradient
                !
                do kposi = 1,pgaus*ndime
                   xvalu(kposi) = 0.0_rp
                end do

                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then

                   do inode = 1,pnode
                      ipoin = lnods(inode,ielem)
                      do idime = 1,ndime
                         elgrp(idime,inode) = prope_ker % grval_ipoin(idime,ipoin)
                      end do
                   end do
                   if( present(gpsha) ) then
                      do igaus = 1,qgaus
                         kposi = (igaus-1) * qnode
                         do inode = 1,pnode
                            kposi = kposi + 1
                            lposi = (igaus-1) * ndime
                            do idime = 1,ndime
                               lposi = lposi + 1
                               xvalu(lposi) = xvalu(lposi) + elgrp(idime,inode) * gpsha(kposi)
                            end do
                         end do
                      end do
                   else
                      do igaus = 1,qgaus
                         do inode = 1,pnode
                            lposi = (igaus-1) * ndime
                            do idime = 1,ndime
                               lposi = lposi + 1
                               xvalu(lposi) = xvalu(lposi) + elgrp(idime,inode) * elmar(pelty) % shape(inode,igaus)
                            end do
                         end do
                      end do
                   end if

                else

                   if( present(gpcar) ) then
                      do igaus = 1,qgaus
                         kposi = (igaus-1) * ndime
                         do idime = 1,ndime
                            kposi = kposi + 1
                            do inode = 1,pnode
                               xvalu(kposi) = xvalu(kposi) + elpro(inode) * gpcar(idime,inode,igaus)
                            end do
                         end do
                      end do
                   end if

                end if

             else
                !
                ! Value
                !
                if( present(gpsha) ) then
                   do igaus = 1,qgaus
                      xvalu(igaus) = 0.0_rp
                      kposi        = (igaus-1) * qnode
                      do inode = 1,pnode
                         kposi = kposi + 1
                         xvalu(igaus) = xvalu(igaus) + elpro(inode) * gpsha(kposi)
                      end do
                   end do
                else
                   do igaus = 1,ngaus(pelty)
                      xvalu(igaus) = 0.0_rp
                      do inode = 1,pnode
                         xvalu(igaus) = xvalu(igaus) + elpro(inode) * elmar(pelty) % shape(inode,igaus)
                      end do
                   end do
                end if
             end if

          else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

             if( kfl_gradi == 1 ) then
                !
                ! Gradient
                !
                do kposi = 1,pgaus*ndime
                   xvalu(kposi) = 0.0_rp
                end do

             else
                !
                ! Value
                !
                if( prope_ker % kfl_type == SCALAR_PROPERTY ) then
                   if( present(gpsha) ) then
                      do igaus = 1,qgaus
                         xvalu(igaus) = prope_ker % value_const(1,imate)
                      end do
                   else
                      do igaus = 1,ngaus(abs(ltype(ielem)))
                         xvalu(igaus) = prope_ker % value_const(1,imate)
                      end do
                   end if
                else
                   kk = 0
                   if( present(gpsha) ) then
                      do igaus = 1,qgaus
                         do idime = 1,prope_ker % dim
                            kk = kk + 1
                            xvalu(kk) = prope_ker % value_const(idime,imate)
                         end do
                      end do
                   else
                      do igaus = 1,ngaus(abs(ltype(ielem)))
                         do idime = 1,prope_ker % dim
                            kk = kk + 1
                            xvalu(kk) = prope_ker % value_const(idime,imate)
                         end do
                      end do
                   end if                   
                end if
             end if

          end if

       end if

       !----------------------------------------------------------------------
       !
       ! 5. Property required on IGAUB
       !
       !----------------------------------------------------------------------

       if( wherein == 'IGAUB' ) then

          if( kfl_gradi == 1 )  call runend('MOD_KER_PROPER: GRADIENT NOT AVILABLE ON BOUNDARY')
          igaub = iposi
          iboun = ielbo
          pblty = abs(ltypb(iboun))
          pnodb = nnode(pblty)
          ielem = lelbo(iboun)
          imate = lmate(ielem)
          ilaws = prope_ker % ilaws(imate)

          if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then

             if( present(gpsha) ) then
                if( qgaus /= ngaus(pblty) ) call runend('MOD_KER_PROPER: IGAUB NOT CODED')
             end if
             xvalu(1) = prope_ker % value_iboun(iboun) % a(igaub)

          else if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

             do inodb = 1,pnodb
                ipoin = lnodb(inodb,iboun)
                elpro(inodb) = prope_ker % value_ipoin(ipoin)
             end do
             xvalu(1) = 0.0_rp
             if( present(gpsha) ) then
                kposi = (igaub-1) * qnode
                do inodb = 1,pnodb
                   kposi = kposi + 1
                   xvalu(1) = xvalu(1) + elpro(inodb) * gpsha(kposi)
                end do
             else
                do inodb = 1,pnodb
                   xvalu(1) = xvalu(1) + elpro(inodb) * elmar(pblty) % shape(inodb,igaub)
                end do
             end if

          else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

             xvalu(1) = prope_ker % value_const(1,imate)

          end if

       end if

       !----------------------------------------------------------------------
       !
       ! 6. Property required on IGAUB=1,PGAUB
       !
       !----------------------------------------------------------------------

       if( wherein == 'PGAUB' ) then

          if( kfl_gradi == 1 ) call runend('MOD_KER_PROPER: GRADIENT NOT AVILABLE ON BOUNDARY')

          iboun = ielbo
          pblty = abs(ltypb(iboun))
          pnodb = nnode(pblty)
          ielem = lelbo(iboun)
          imate = lmate(ielem)
          ilaws = prope_ker % ilaws(imate)

          if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then

             if( present(gpsha) ) then
                if( qgaus /= ngaus(pblty) ) call runend('MOD_KER_PROPER: IGAUB NOT CODED')
             end if
             do igaub = 1,ngaus(pblty)
                xvalu(igaub) = prope_ker % value_iboun(iboun) % a(igaub)
             end do

          else if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

             do inodb = 1,pnodb
                ipoin = lnodb(inodb,iboun)
                elpro(inodb) = prope_ker % value_ipoin(ipoin)
             end do

             if( present(gpsha) ) then
                do igaub = 1,qgaus
                   xvalu(igaub) = 0.0_rp
                   kposi = (igaub-1) * qnode
                   do inodb = 1,pnodb
                      kposi        = kposi + 1
                      xvalu(igaub) = xvalu(igaub) + elpro(inodb) * gpsha(kposi)
                   end do
                end do
             else
                pgaub = ngaus(pblty)
                do igaub = 1,pgaub
                   xvalu(igaub) = 0.0_rp
                   do inodb = 1,pnodb
                      xvalu(igaub) = xvalu(igaub) + elpro(inodb) * elmar(pblty) % shape(inodb,igaub)
                   end do
                end do
             end if

          else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

             do igaub = 1,ngaus(abs(ltypb(iboun)))
                xvalu(igaub) = prope_ker % value_const(1,imate)
             end do

          end if

       else if( wherein(1:3) == 'COG' ) then

          !----------------------------------------------------------------------
          !
          ! 7. Property required on C.O.G.
          !
          !----------------------------------------------------------------------

          ielem = ielbo
          imate = lmate(ielem)
          ilaws = prope_ker % ilaws(imate)
          pelty = abs(ltype(ielem))
          pgaus = ngaus(pelty)
          pnode = nnode(pelty) 
          if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then

             if( kfl_gradi == 1 )  then
                !
                ! Gradient
                !
                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then
                   xvalu(1:ndime) = 0.0_rp
                   do igaus = 1,pgaus
                      do idime = 1,ndime
                         xvalu(idime) = xvalu(idime) + prope_ker % grval_ielem(ielem) % a(idime,igaus)
                      end do
                   end do
                   do idime = 1,ndime
                      xvalu(idime) = xvalu(idime) / real(pgaus,rp)
                   end do
                else
                   call runend('MOD_KER_PROPER: GRDAIENT ON COG NOT CODED')
                end if

             else
                !
                ! Value
                !
                xvalu(1) = 0.0_rp
                if( present(gpsha) ) then
                   do igaus = 1,qgaus
                      xvalu(1) = xvalu(1) + prope_ker % value_ielem(ielem) % a(igaus)
                   end do
                   xvalu(1) = xvalu(1) / real(qgaus,rp)
                else
                   do igaus = 1,pgaus
                      xvalu(1) = xvalu(1) + prope_ker % value_ielem(ielem) % a(igaus)
                   end do
                   xvalu(1) = xvalu(1) / real(pgaus,rp)
                end if
             end if

          else if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

             if( kfl_gradi == 1 )  then
                !
                ! Gradient
                !
                xvalu(1:ndime) = 0.0_rp

                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then

                   do inode = 1,pnode
                      ipoin = lnods(inode,ielem)
                      do idime = 1,ndime
                         xvalu(idime) = xvalu(idime) + prope_ker % grval_ipoin(idime,ipoin)
                      end do
                   end do
                   do idime = 1,ndime
                      xvalu(idime) = xvalu(idime) / real(pnode,rp)
                   end do

                else

                   if( present(gpcar) ) then
                      do igaus = 1,qgaus
                         do inode = 1,pnode
                            ipoin = lnods(inode,ielem)
                            do idime = 1,ndime
                               xvalu(idime) = xvalu(idime) &
                                    + prope_ker % value_ipoin(ipoin) * gpcar(idime,inode,1)
                            end do
                         end do
                      end do
                      do idime = 1,ndime
                         xvalu(idime) = xvalu(idime) / real(qgaus,rp)
                      end do
                   else
                      call runend('MOD_KER_PROPER: GRADIENT NOT AVAILABLE AT COG')
                   end if

                end if

             else
                !
                ! Value
                !
                xvalu(1) = 0.0_rp
                do inode = 1,pnode
                   ipoin    = lnods(inode,ielem)
                   xvalu(1) = xvalu(1) + prope_ker % value_ipoin(ipoin)
                end do
                xvalu(1) = xvalu(1) / real(pnode,rp)
             end if

          else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

             if( kfl_gradi == 1 )  then
                xvalu(1:ndime) = 0.0_rp
             else
                !if( .not. associated(prope_ker) ) call runend('PROB1')
                !if( .not. associated(prope_ker % value_const) ) call runend('PROB2')
                xvalu(1) = prope_ker % value_const(1,imate)
             end if

          end if

       end if

       !----------------------------------------------------------------------
       !
       ! 8. Property required on PNODE
       !
       !----------------------------------------------------------------------

       if( wherein == 'PNODE' ) then

          ielem = ielbo
          pelty = abs(ltype(ielem))
          imate = lmate(ielem)
          ilaws = prope_ker % ilaws(imate)
          pgaus = ngaus(pelty)
          pnode = nnode(pelty)

          if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then

             if( present(gpsha) ) then
                if( pgaus /= qgaus ) call runend('MOD_KER_PROPER: COG NOT CODED')
             end if

             if( kfl_gradi == 1 ) then
                !
                ! Gradient
                !
                do inode = 1,ndime*pnode
                   xvalu(inode) = 0.0_rp
                end do

                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then

                   do igaus = 1,pgaus
                      do inode = 1,pnode
                         kposi = (inode-1) * ndime
                         do idime = 1,ndime
                            kposi = kposi + 1
                            xvalu(kposi) = xvalu(kposi) + prope_ker % grval_ielem(ielem) % a(idime,igaus) &
                                 * elmar(pelty) % shaga(igaus,inode)
                         end do
                      end do
                   end do

                else

                   do inode = 1,pnode
                      ipoin = lnods(inode,ielem)
                      elpro(inode) = 0.0_rp
                      do igaus = 1,pgaus
                         elpro(inode) = elpro(inode) + prope_ker % value_ielem(ielem) % a(igaus) &
                              * elmar(pelty) % shaga(igaus,inode)
                      end do
                   end do
                   do igaus = 1,pgaus
                      do inode = 1,pnode
                         do idime = 1,ndime
                            elgrp(idime,igaus) = elgrp(idime,igaus) + elpro(inode) * elmar(pelty) % shape(inode,igaus)
                         end do
                      end do
                   end do
                   do igaus = 1,pgaus
                      do inode = 1,pnode
                         kposi = (inode-1) * ndime
                         do idime = 1,ndime
                            kposi = kposi + 1
                            xvalu(kposi) = xvalu(kposi) + elgrp(idime,igaus) * elmar(pelty) % shaga(igaus,inode)
                         end do
                      end do
                   end do

                end if

             else
                !
                ! Value
                !
                do inode = 1,pnode
                   xvalu(inode) = 0.0_rp
                   do igaus = 1,pgaus
                      xvalu(inode) = xvalu(inode) + prope_ker % value_ielem(ielem) % a(igaus) &
                           * elmar(pelty) % shaga(igaus,inode)
                   end do
                end do

             end if

          else if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

             if( kfl_gradi == 1 ) then
                !
                ! Gradient
                !
                do inode = 1,ndime*pnode
                   xvalu(inode) = 0.0_rp
                end do
                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then
                   do inode = 1,pnode
                      kposi = (inode-1) * ndime
                      do idime = 1,ndime
                         kposi = kposi + 1
                         xvalu(kposi) = prope_ker % grval_ipoin(idime,ipoin)
                      end do
                   end do
                else
                   if( .not. present(gpcar) ) call runend('MOD_KER_PROPER: GPCAR MISSING')
                   do igaus = 1,pgaus
                      do inode = 1,pnode
                         ipoin = lnods(inode,ielem)
                         do idime = 1,ndime
                            elgrp(idime,igaus) = elgrp(idime,igaus) + prope_ker % value_ipoin(ipoin) * gpcar(idime,inode,igaus)
                         end do
                      end do
                   end do
                   do igaus = 1,pgaus
                      do inode = 1,pnode
                         kposi = (inode-1) * ndime
                         do idime = 1,ndime
                            kposi = kposi + 1
                            xvalu(kposi) = xvalu(kposi) + elgrp(idime,igaus) * elmar(pelty) % shaga(igaus,inode)
                         end do
                      end do
                   end do
                end if
             else
                !
                ! Value
                !
                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   xvalu(inode) = prope_ker % value_ipoin(ipoin)
                end do
             end if

          else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

             if( kfl_gradi == 1 ) then
                do inode = 1,ndime*pnode
                   xvalu(inode) = 0.0_rp
                end do
             else
                do inode = 1,pnode
                   xvalu(inode) = prope_ker % value_const(1,imate)
                end do
             end if

          end if

       end if

       !----------------------------------------------------------------------
       !
       ! 9. Property required Anywhere
       !
       !----------------------------------------------------------------------

       if( wherein == 'ANYWH' ) then

          igaus = iposi
          ielem = ielbo
          imate = lmate(ielem)
          ilaws = prope_ker % ilaws(imate)

          if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then

             call runend('MOD_KER_PROPER: TEST THIS WITH LAGRANGE')

          else if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

             pelty = abs(ltype(ielem))
             pnode = lnnod(ielem)
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                elpro(inode) = prope_ker % value_ipoin(ipoin)
             end do
             xvalu(1) = 0.0_rp
             if( present(gpsha) ) then
                kposi    = (igaus-1) * qnode
                do inode = 1,pnode
                   kposi    = kposi + 1
                   xvalu(1) = xvalu(1) + elpro(inode) * gpsha(kposi)
                end do
             else
                call runend('MOD_KER_PROPER: NOT CODED')
             end if

          else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

             xvalu(1) = prope_ker % value_const(1,imate)

          end if

       end if

    end if
    !nullify(prope_ker) ! This avoids the bug of calling wrongly the ker_proper routine

  end subroutine ker_proper_scalar

  subroutine ker_proper_vector_0(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor,VECTOR_DIM)
    character(5),        intent(in)           :: wname
    character(5),        intent(in)           :: wherein
    integer(ip),         intent(in)           :: iposi
    integer(ip),         intent(in)           :: ielbo(VECTOR_SIZE)
    real(rp),            intent(out)          :: xvalu(VECTOR_SIZE)
    integer(ip),         intent(in), optional :: qnode
    integer(ip),         intent(in), optional :: qgaus
    real(rp),            intent(in), optional :: gpsha(VECTOR_SIZE,1,*)
    real(rp),            intent(in), optional :: gpcar(VECTOR_SIZE,1,1,*)
    real(rp),            intent(in), optional :: gpcor(VECTOR_SIZE,1,*)
    integer(ip),                     optional :: VECTOR_DIM
    integer(ip)                               :: VECTOR_SIZE_LOC

    if( present(VECTOR_DIM) ) then
       VECTOR_SIZE_LOC = VECTOR_DIM
    else
       VECTOR_SIZE_LOC = VECTOR_SIZE       
    end if

    call ker_proper_vector(&
         wname,wherein,iposi,ielbo,xvalu,&
         qnode,qgaus,gpsha,gpcar,gpcor,VECTOR_SIZE_LOC)

  end subroutine ker_proper_vector_0

  subroutine ker_proper_vector_1(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor,VECTOR_DIM)
    character(5),        intent(in)           :: wname
    character(5),        intent(in)           :: wherein
    integer(ip),         intent(in)           :: iposi
    integer(ip),         intent(in)           :: ielbo(VECTOR_SIZE)
    real(rp),            intent(out)          :: xvalu(VECTOR_SIZE,*)
    integer(ip),         intent(in)           :: qnode
    integer(ip),         intent(in)           :: qgaus
    real(rp),            intent(in), optional :: gpsha(VECTOR_SIZE,1,*)
    real(rp),            intent(in), optional :: gpcar(VECTOR_SIZE,1,1,*)
    real(rp),            intent(in), optional :: gpcor(VECTOR_SIZE,1,*)
    integer(ip),                     optional :: VECTOR_DIM
    integer(ip)                               :: VECTOR_SIZE_LOC

    if( present(VECTOR_DIM) ) then
       VECTOR_SIZE_LOC = VECTOR_DIM
    else
       VECTOR_SIZE_LOC = VECTOR_SIZE       
    end if

    call ker_proper_vector(&
         wname,wherein,iposi,ielbo,xvalu,&
         qnode,qgaus,gpsha,gpcar,gpcor,VECTOR_SIZE_LOC)

  end subroutine ker_proper_vector_1

  subroutine ker_proper_vector_2(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor,VECTOR_DIM)
    character(5),        intent(in)           :: wname
    character(5),        intent(in)           :: wherein
    integer(ip),         intent(in)           :: iposi
    integer(ip),         intent(in)           :: ielbo(VECTOR_SIZE)
    real(rp),            intent(out)          :: xvalu(VECTOR_SIZE,1,*)
    integer(ip),         intent(in)           :: qnode
    integer(ip),         intent(in)           :: qgaus
    real(rp),            intent(in), optional :: gpsha(VECTOR_SIZE,1,*)
    real(rp),            intent(in), optional :: gpcar(VECTOR_SIZE,1,1,*)
    real(rp),            intent(in), optional :: gpcor(VECTOR_SIZE,1,*)
    integer(ip),                     optional :: VECTOR_DIM
    integer(ip)                               :: VECTOR_SIZE_LOC

    if( present(VECTOR_DIM) ) then
       VECTOR_SIZE_LOC = VECTOR_DIM
    else
       VECTOR_SIZE_LOC = VECTOR_SIZE       
    end if

    call ker_proper_vector(&
         wname,wherein,iposi,ielbo,xvalu,&
         qnode,qgaus,gpsha,gpcar,gpcor,VECTOR_SIZE_LOC)

  end subroutine ker_proper_vector_2

  subroutine ker_proper_vector_3(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor,VECTOR_DIM)

    character(5),        intent(in)           :: wname
    character(5),        intent(in)           :: wherein
    integer(ip),         intent(in)           :: iposi
    integer(ip),         intent(in)           :: ielbo(VECTOR_SIZE)
    real(rp),            intent(out)          :: xvalu(VECTOR_SIZE,1,1,*)
    integer(ip),         intent(in)           :: qnode
    integer(ip),         intent(in)           :: qgaus
    real(rp),            intent(in), optional :: gpsha(VECTOR_SIZE,1,*)
    real(rp),            intent(in), optional :: gpcar(VECTOR_SIZE,1,1,*)
    real(rp),            intent(in), optional :: gpcor(VECTOR_SIZE,1,*)
    integer(ip),                     optional :: VECTOR_DIM
    integer(ip)                               :: VECTOR_SIZE_LOC

    if( present(VECTOR_DIM) ) then
       VECTOR_SIZE_LOC = VECTOR_DIM
    else
       VECTOR_SIZE_LOC = VECTOR_SIZE       
    end if

    call ker_proper_vector(&
         wname,wherein,iposi,ielbo,xvalu,&
         qnode,qgaus,gpsha,gpcar,gpcor,VECTOR_SIZE_LOC)

  end subroutine ker_proper_vector_3

  subroutine ker_proper_scalar_00(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor)
    ! For example: only on one Gauss point (IGAUS)
    character(5),        intent(in)           :: wname
    character(5),        intent(in)           :: wherein
    integer(ip),         intent(in)           :: iposi
    integer(ip),         intent(in)           :: ielbo
    real(rp),            intent(out)          :: xvalu(*)
    integer(ip),         intent(in)           :: qnode
    integer(ip),         intent(in)           :: qgaus
    real(rp),            intent(in)           :: gpsha(*)
    real(rp),            intent(in), optional :: gpcar(1,*)
    real(rp),            intent(in), optional :: gpcor(*)
    call ker_proper_scalar(&
         wname,wherein,iposi,ielbo,xvalu,&
         qnode,qgaus,gpsha,gpcar,gpcor)
  end subroutine ker_proper_scalar_00

  subroutine ker_proper_scalar_000(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor)
    ! For example: only on one Gauss point (IGAUS)
    character(5),        intent(in)           :: wname
    character(5),        intent(in)           :: wherein
    integer(ip),         intent(in)           :: iposi
    integer(ip),         intent(in)           :: ielbo
    real(rp),            intent(out)          :: xvalu
    integer(ip),         intent(in)           :: qnode
    integer(ip),         intent(in)           :: qgaus
    real(rp),            intent(in)           :: gpsha(*)
    real(rp),            intent(in), optional :: gpcar(1,*)
    real(rp),            intent(in), optional :: gpcor(*)
    real(rp)                                  :: xval2(2)
    call ker_proper_scalar(&
         wname,wherein,iposi,ielbo,xval2,&
         qnode,qgaus,gpsha,gpcar,gpcor)
    xvalu = xval2(1)
  end subroutine ker_proper_scalar_000

  subroutine ker_proper_scalar_0(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor)
    character(5),        intent(in)           :: wname
    character(5),        intent(in)           :: wherein
    integer(ip),         intent(in)           :: iposi
    integer(ip),         intent(in)           :: ielbo
    real(rp),            intent(out)          :: xvalu(*)
    integer(ip),         intent(in), optional :: qnode
    integer(ip),         intent(in), optional :: qgaus
    real(rp),            intent(in), optional :: gpsha(1,*)
    real(rp),            intent(in), optional :: gpcar(1,1,*)
    real(rp),            intent(in), optional :: gpcor(1,*)

    if(  .not. present(qnode) .and. &
         .not. present(qgaus) .and. &
         .not. present(gpsha) .and. &
         .not. present(gpcar) .and. &
         .not. present(gpcor) ) then
       call ker_proper_scalar(&
            wname,wherein,iposi,ielbo,xvalu)
    else
       call ker_proper_scalar(&
            wname,wherein,iposi,ielbo,xvalu,&
            qnode,qgaus,gpsha,gpcar,gpcor)
    end if
  end subroutine ker_proper_scalar_0

  subroutine ker_proper_scalar_1(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor)
    character(5),        intent(in)           :: wname
    character(5),        intent(in)           :: wherein
    integer(ip),         intent(in)           :: iposi
    integer(ip),         intent(in)           :: ielbo
    real(rp),            intent(out)          :: xvalu(1,*)
    integer(ip),         intent(in), optional :: qnode
    integer(ip),         intent(in), optional :: qgaus
    real(rp),            intent(in), optional :: gpsha(1,*)
    real(rp),            intent(in), optional :: gpcar(1,1,*)
    real(rp),            intent(in), optional :: gpcor(1,*)
    call ker_proper_scalar(&
         wname,wherein,iposi,ielbo,xvalu,&
         qnode,qgaus,gpsha,gpcar,gpcor)
  end subroutine ker_proper_scalar_1

  subroutine ker_proper_scalar_2(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor)
    character(5),        intent(in)           :: wname
    character(5),        intent(in)           :: wherein
    integer(ip),         intent(in)           :: iposi
    integer(ip),         intent(in)           :: ielbo
    real(rp),            intent(out)          :: xvalu(1,1,*)
    integer(ip),         intent(in), optional :: qnode
    integer(ip),         intent(in), optional :: qgaus
    real(rp),            intent(in), optional :: gpsha(1,*)
    real(rp),            intent(in), optional :: gpcar(1,1,*)
    real(rp),            intent(in), optional :: gpcor(1,*)
    call ker_proper_scalar(&
         wname,wherein,iposi,ielbo,xvalu,&
         qnode,qgaus,gpsha,gpcar,gpcor)
  end subroutine ker_proper_scalar_2

  subroutine ker_proper_vector(&
       wname,wherein,iposi,ielbo,xvalu,&
       qnode,qgaus,gpsha,gpcar,gpcor,VECTOR_SIZE_LOC)
    !-----------------------------------------------------------------------
    !****f* Kermod/ker_proper
    ! NAME
    !    ker_proper
    ! DESCRIPTION
    !    Update properties
    !
    !    Property can be defined with 3 formats
    !    It can be required at 8 different places
    !
    !       WHEREIN   IPOSI  IELBO
    !       --------------------
    !    1. IPOIN   IPOIN  NONSENSE
    !    2. NPOIN   DUMMI  Do smoothing => Loop over IELEM
    !    3. IGAUS   IGAUS  IELEM
    !    4. PGAUS   DUMMI  IELEM
    !    5. IGAUB   IGAUB  IBOUN        => IELEM = LBOEL(,IBOUN)
    !    6. PGAUB   DUMMI  IBOUN        => IELEM = LBOEL(,IBOUN)
    !    7. COG     DUMMI  IELEM
    !    8. PNODE   DUMMI  IELEM
    !
    !    GPSHA is dimensioned as a vector in order to use it as optional
    !    argument
    !
    ! USED BY
    !    many subroutines
    !***
    !-----------------------------------------------------------------------

    use def_kintyp,      only                 :  ip,rp
    use def_elmtyp,      only                 :  ELCUT
    use def_master,      only                 :  nturb
    use def_domain,      only                 :  npoin,mnode,lelch
    use def_domain,      only                 :  lnods,lnnod,lmate,lnnod
    use def_domain,      only                 :  ngaus,ltype,elmar,ndime
    use def_domain,      only                 :  nnode
    use def_kermod,      only                 :  typ_valpr_ker
    use def_kermod,      only                 :  densi_ker,visco_ker
    use def_kermod,      only                 :  poros_ker,condu_ker
    use def_kermod,      only                 :  sphea_ker,dummy_ker
    use def_kermod,      only                 :  turmu_ker,mixin_ker
    use def_kermod,      only                 :  anipo_ker,walvi_ker
    use mod_memory

#ifdef OPENACCHHH
    use openacc
#endif
    implicit none
    integer(ip),         intent(in)           :: VECTOR_SIZE_LOC
    character(5),        intent(in)           :: wname
    character(5),        intent(in)           :: wherein
    integer(ip),         intent(in)           :: iposi
    integer(ip),         intent(in)           :: ielbo(VECTOR_SIZE_LOC)
    real(rp),            intent(out)          :: xvalu(VECTOR_SIZE_LOC,*)
    integer(ip),         intent(in)           :: qnode
    integer(ip),         intent(in)           :: qgaus
    real(rp),            intent(in), optional :: gpsha(VECTOR_SIZE_LOC,*) !!DMM ...,qnode)
    real(rp),            intent(in), optional :: gpcar(VECTOR_SIZE_LOC,ndime,mnode,*)
    real(rp),            intent(in), optional :: gpcor(VECTOR_SIZE_LOC,ndime,*)
    type(typ_valpr_ker), pointer              :: prope_ker
    integer(ip)                               :: imate
    integer(ip)                               :: ielem(VECTOR_SIZE_LOC)
    integer(ip)                               :: pgaus
    integer(ip)                               :: inode,ipoin,igaus,iturb,kelem
    integer(ip)                               :: ilaws,pelty,pnode,kk
    integer(ip)                               :: ivect
    integer(ip)                               :: kposi,kfl_gradi,idime,kfl_deriv,kfl_grder
    integer(ip)                               :: kfl_deriv_tur,kfl_deriv_vel,jposi
    integer(ip)                               :: lposi
    real(rp)                                  :: elpro(VECTOR_SIZE_LOC,mnode)
    real(rp)                                  :: elgrp(VECTOR_SIZE_LOC,ndime,mnoga)
    real(rp),   pointer                       :: auxva(:)
    integer(ip)                               :: ngaus_size
#ifdef OPENACCHHH
    logical                                   :: lpalloc
!    logical                                   :: lpalloc_ind
!    logical                                   :: lpalloc_stru
#endif



#ifdef OPENACCHHH
#define DEF_VECT ivect
#else
#define DEF_VECT 1:VECTOR_SIZE_LOC
#endif

    nullify(auxva)

    !----------------------------------------------------------------------
    !
    ! What to compute
    !
    !----------------------------------------------------------------------

    if( wname(1:2) == 'GR' ) then
       kfl_gradi = 1
    else
       kfl_gradi = 0
    end if

    if( wname(1:2) == 'DR' ) then
       kfl_deriv = 1
    else
       kfl_deriv = 0
    end if

    if( wname(1:2) == 'GD' ) then
       kfl_grder = 1
    else
       kfl_grder = 0
    end if

    if( wname(1:3) == 'TDR' ) then
       kfl_deriv_tur = 1
    else
       kfl_deriv_tur = 0
    end if

    if( wname(1:3) == 'VDR' ) then
       kfl_deriv_vel = 1
    else
       kfl_deriv_vel = 0
    end if

    !----------------------------------------------------------------------
    !
    ! Which property to compute
    !
    !----------------------------------------------------------------------

    if(      wname == 'DENSI' .or. wname == 'GRDEN' .or. wname == 'DRDEN' .or. wname == 'GDDEN') then
       prope_ker => densi_ker
    else if( wname == 'VISCO' .or. wname == 'GRVIS' .or. wname == 'DRVIS' .or. wname == 'GDVIS') then
       prope_ker => visco_ker
   else if( wname == 'MIXIN' ) then
       prope_ker => mixin_ker
    else if( wname == 'POROS' .or. wname == 'GRPOS' ) then
       prope_ker => poros_ker
    else if( wname == 'CONDU' .or. wname == 'GRCON' ) then
       prope_ker => condu_ker
    else if( wname == 'SPHEA' .or. wname == 'GRSPE' ) then
       prope_ker => sphea_ker
    else if( wname == 'DUMMY' .or. wname == 'GRDUM' ) then
       prope_ker => dummy_ker
    else if( wname == 'TURBU' .or. wname == 'TURVI' .or. wname == 'GRTUR' .or. wname == 'TDRTU' .or. wname == 'VDRTU' ) then
       prope_ker => turmu_ker
    else if( wname == 'ABSOR' ) then
       prope_ker => absor_ker
    else if( wname == 'SCATT' ) then
       prope_ker => scatt_ker
    else if( wname == 'ANIPO' ) then
       prope_ker => anipo_ker
    else if( wname == 'WALLV' ) then
       prope_ker => walvi_ker
    else
       print *,'--- I DONT UNDERSTAND ',wname
       call runend('KER_PROPER CALLED WITH WRONG PROPERTY') ! This avoids the bug of calling wrongly the ker_proper routine
    end if

#ifdef OPENACCHHH
    if (associated(prope_ker) .and. associated(prope_ker % value_const)) then
      !$acc enter data copyin(prope_ker, prope_ker%value_const)
    else if (associated(prope_ker)) then
      !$acc enter data copyin(prope_ker)
    end if
#endif

    if( prope_ker % kfl_exist == 0 ) then
       !
       ! Put value to zero
       !
       if( wherein == 'IPOIN' .or. wherein == 'IGAUS' .or. wherein == 'IGAUB' .or. wherein == 'COG' ) then
          kposi = 1
       else if( wherein == 'NPOIN' ) then
          kposi = npoin
       else if( wherein == 'PNODE' ) then
          kposi = lnnod(ielbo(1))
       else if( wherein == 'PGAUS' ) then
          if( present(gpsha) ) then
             kposi = qgaus
          else
             kposi = ngaus(abs(ltype(ielbo(1))))
          end if
       else
          call runend('KER_PROPER_VECTOR: NOT CODED 2')
       end if
       if( kfl_gradi == 1 ) kposi = kposi * ndime
       do kelem = 1,kposi
          xvalu(1:VECTOR_SIZE_LOC,kelem) = 0.0_rp
       end do

    else

       if( wherein == 'IGAUS' ) then

          !----------------------------------------------------------------------
          !
          ! 3. Property required on IGAUS
          !    Note: it has not been validated 03/03/2017
          !----------------------------------------------------------------------

          igaus = iposi
          ielem = ielbo(1)
          imate = lmate(ielem(1))
          ilaws = prope_ker % ilaws(imate)

          if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then

             if( kfl_gradi == 1 ) then
                !
                ! Gradient
                !
                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then
                   do ivect = 1,VECTOR_SIZE_LOC
                      if( ielem(ivect) > 0 ) then
                         do idime = 1,ndime
                            xvalu(ivect,idime) = prope_ker % grval_ielem(ielem(ivect)) % a(idime,igaus)
                         end do
                      end if
                   end do
                else
                   xvalu(1:VECTOR_SIZE_LOC,1:ndime) = 0.0_rp
                end if

             else
                !
                ! Value
                !
                if( present(gpsha) ) then
                   do ivect = 1,VECTOR_SIZE_LOC
                      kelem = ielem(ivect)
                      if( kelem > 0 ) then
                         pelty = abs(kelem)
                         pnode = lnnod(kelem)
                         pgaus = ngaus(pelty)
                         xvalu(ivect,1) = 0.0_rp
                         do inode = 1,pnode
                            do igaus = 1,pgaus
                               xvalu(ivect,1) = xvalu(ivect,1) + gpsha(ivect,inode) * elmar(pelty) % shaga(igaus,inode) &
                                    * prope_ker % value_ielem(kelem) % a(igaus)
                            end do
                         end do
                      end if
                   end do
                else
                   do ivect = 1,VECTOR_SIZE_LOC
                      if( ielem(ivect) > 0 ) then
                         xvalu(ivect,1) = prope_ker % value_ielem(ielem(ivect)) % a(igaus)
                      end if
                   end do
                end if

             end if

          else if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

             call runend('KER_PROPER: NOT AVAILABLE IGAUS 2')

!!$             pelty = abs(ltype(ielem))
!!$             pnode = lnnod(ielem)
!!$             do inode = 1,pnode
!!$                ipoin = lnods(inode,ielem)
!!$                elpro(inode) = prope_ker % value_ipoin(ipoin)
!!$             end do
!!$
!!$             if( kfl_gradi == 1 ) then
!!$                !
!!$                ! Gradient
!!$                !
!!$                xvalu(1:ndime) = 0.0_rp
!!$
!!$                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then
!!$
!!$                   do inode = 1,pnode
!!$                      ipoin = lnods(inode,ielem)
!!$                      do idime = 1,ndime
!!$                         elgrp(idime,inode) = prope_ker % grval_ipoin(idime,ipoin)
!!$                      end do
!!$                   end do
!!$                   if( present(gpsha) ) then
!!$                      kposi = (igaus-1) * qnode
!!$                      do inode = 1,pnode
!!$                         kposi = kposi + 1
!!$                         do idime = 1,ndime
!!$                            xvalu(idime) = xvalu(idime) + elgrp(idime,inode) * gpsha(kposi)
!!$                         end do
!!$                      end do
!!$                   else
!!$                      do inode = 1,pnode
!!$                         do idime = 1,ndime
!!$                            xvalu(idime) = xvalu(idime) + elgrp(idime,inode) * elmar(pelty) % shape(idime,igaus)
!!$                         end do
!!$                      end do
!!$                   end if
!!$
!!$                else
!!$
!!$                   if( present(gpcar) ) then
!!$                      do inode = 1,pnode
!!$                         do idime = 1,ndime
!!$                            xvalu(idime) = xvalu(idime) + elpro(inode) * gpcar(idime,inode,1)
!!$                         end do
!!$                      end do
!!$                   else
!!$                      call runend('MOD_KER_PROPER: GPCAR NEEDED')
!!$                   end if
!!$                end if
!!$
!!$          else
!!$             !
!!$             ! Value
!!$             !
!!$             call runend('KER_PROPER: NOT AVAILABLE IGAUS 3')
!!$
!!$                xvalu(1) = 0.0_rp
!!$                if( present(gpsha) ) then
!!$                   kposi    = (igaus-1) * qnode
!!$                   do inode = 1,pnode
!!$                      kposi    = kposi + 1
!!$                      xvalu(1) = xvalu(1) + elpro(inode) * gpsha(kposi)
!!$                   end do
!!$                else
!!$                   do inode = 1,pnode
!!$                      xvalu(1) = xvalu(1) + elpro(inode) * elmar(pelty) % shape(inode,igaus)
!!$                   end do
!!$                end if

          else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

             if( kfl_gradi == 1 ) then
                xvalu(1:VECTOR_SIZE_LOC,1:ndime) = 0.0_rp
             else
                xvalu(1:VECTOR_SIZE_LOC,1) = prope_ker % value_const(1,imate)
             end if

          end if

       else if( wherein == 'PGAUS' ) then

          !----------------------------------------------------------------------
          !
          ! 4. Property required on IGAUS=1,PGAUS
          !    Note: it has been validated 03/03/2017
          !----------------------------------------------------------------------

          ielem(1:VECTOR_SIZE_LOC) = ielbo(1:VECTOR_SIZE_LOC)
          imate                = lmate(ielem(1))
          ilaws                = prope_ker % ilaws(imate)
          pelty                = abs(ltype(ielem(1)))
          pgaus                = ngaus(pelty)
          if( lelch(ielem(1)) == ELCUT ) pgaus = lgaus(ielem(1))

          if( prope_ker % llaws(ilaws) % where == 'CONST' ) then
             !
             ! Property is constant
             !

             if( kfl_gradi == 1 ) then
                !
                ! Gradient
                !
                do kposi = 1,pgaus*ndime*prope_ker % dim
                   xvalu(1:VECTOR_SIZE_LOC,kposi) = 0.0_rp
                end do

             else
                !
                ! Value
                !
                if( prope_ker % kfl_type == SCALAR_PROPERTY ) then

                   if( present(gpsha) ) then
#ifdef OPENACCHHH
#ifdef _OPENACC
                    lpalloc = acc_is_present(xvalu) 
                    if ( .not. lpalloc) then
                        !$acc enter data copyin (xvalu(1:VECTOR_SIZE_LOC,1:qgaus))
                    end if 
#endif
                    !$acc parallel loop gang vector default(present) 
                    do ivect = 1, VECTOR_SIZE_LOC
#endif                        
                      do igaus = 1,qgaus
                         xvalu(DEF_VECT,igaus) = prope_ker % value_const(1,imate)
                      end do

#ifdef OPENACCHHH                      
                    end do
                    !$acc end parallel loop

                    if ( .not. lpalloc) then
                        !$acc update self (xvalu(1:VECTOR_SIZE_LOC,1:qgaus))
                        !$acc exit data delete (xvalu(1:VECTOR_SIZE_LOC,1:qgaus)) 
                    end if 
#endif                 
                   else

                   ngaus_size = ngaus(abs(ltype(ielem(1))))


#ifdef OPENACCHHH                       
#ifdef _OPENACC
                   lpalloc = acc_is_present(xvalu)
#endif
                    if ( .not. lpalloc) then
                        !$acc enter data copyin (xvalu(1:VECTOR_SIZE_LOC,1:ngaus_size))
                    end if 

                   !$acc parallel loop gang vector default(present) 
                   do ivect = 1, VECTOR_SIZE_LOC
#endif                       
                      do igaus = 1,ngaus_size 
                         xvalu(DEF_VECT,igaus) = prope_ker % value_const(1,imate)
                      end do
#ifdef OPENACCHHH                       
                   end do
                   !$acc end parallel loop

                   if ( .not. lpalloc) then
                        !$acc update self (xvalu(1:VECTOR_SIZE_LOC,1:ngaus_size))
                        !$acc exit data delete (xvalu(1:VECTOR_SIZE_LOC,1:ngaus_size)) 
                   end if 


#endif
                   end if
                else

                   kk = 0
                   if( present(gpsha) ) then
                      do igaus = 1,qgaus
                         do idime = 1,prope_ker % dim
                            kk = kk + 1
                            xvalu(1:VECTOR_SIZE_LOC,kk) = prope_ker % value_const(idime,imate)
                         end do
                      end do

                   else
                      do igaus = 1,ngaus(abs(ltype(ielem(1))))
                         do idime = 1,prope_ker % dim
                            kk = kk + 1
                            xvalu(1:VECTOR_SIZE_LOC,kk) = prope_ker % value_const(idime,imate)
                         end do
                      end do

                   end if
                end if

             end if
             
          else if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then
             !
             ! Property is computed element-wise
             !

             if( kfl_gradi == 1 ) then

                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then

                   if( present(gpsha) ) then
                      if( qgaus /= pgaus .and. kfl_cutel == 0 ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                      do ivect = 1,VECTOR_SIZE_LOC
                         do igaus = 1,qgaus
                            jposi = (igaus - 1) * ndime
                            do idime = 1,ndime
                               kposi = jposi + idime
                               xvalu(ivect,kposi) = prope_ker % grval_ielem(ielem(ivect)) % a(idime,igaus)
                            end do
                         end do
                      end do
                   else
                      if( qgaus /= pgaus .and. kfl_cutel == 0 ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                      do ivect = 1,VECTOR_SIZE_LOC
                         do igaus = 1,qgaus
                            jposi = (igaus - 1) * ndime
                            do idime = 1,ndime
                               kposi = jposi + idime
                               xvalu(ivect,kposi) = prope_ker % grval_ielem(ielem(ivect)) % a(idime,igaus)
                            end do
                         end do
                      end do
                   end if
                else
                   do kposi = 1,pgaus*ndime
                      xvalu(1:VECTOR_SIZE_LOC,kposi) = 0.0_rp
                   end do
                end if

             else if( kfl_deriv == 1 ) then

                call runend('KER_PROPER_VECTOR: NOT CODED 1')

                if( prope_ker % llaws(ilaws) % kfl_deriv == 1 ) then

                   if( present(gpsha) ) then
                      if( qgaus /= pgaus ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                      do ivect = 1,VECTOR_SIZE_LOC
                         do igaus = 1,qgaus
                            kposi = (igaus - 1) * qnode
                            do inode = 1,qnode
                               kposi = kposi + 1
                               xvalu(ivect,kposi) = prope_ker % drval_ielem(ielem(ivect)) % a(inode,igaus)
                            end do
                         end do
                      end do
                   else
                      do ivect = 1,VECTOR_SIZE_LOC
                         do igaus = 1,pgaus
                            kposi = (igaus - 1) * qnode
                            do inode = 1,qnode
                               kposi = kposi + 1
                               xvalu(ivect,kposi) = prope_ker % drval_ielem(ielem(ivect)) % a(inode,igaus)
                            end do
                         end do
                      end do
                   end if
                else
                   do kposi = 1,pgaus*qnode
                      xvalu(1:VECTOR_SIZE_LOC,kposi) = 0.0_rp
                   end do
                end if

             else if( kfl_grder == 1 ) then

                if( prope_ker % llaws(ilaws) % kfl_grder == 1 ) then

                   if( present(gpsha) ) then
                      if( qgaus /= pgaus ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                      do ivect = 1,VECTOR_SIZE_LOC
                         do igaus = 1,qgaus
                            kposi = (igaus - 1)* qnode*ndime
                            do inode = 1,qnode
                               do idime = 1, ndime
                                  kposi = kposi + 1
                                  xvalu(ivect,kposi) = prope_ker % gdval_ielem(ielem(ivect)) % a(idime,inode,igaus)
                               end do
                            end do
                         end do
                      end do
                   else
                      do ivect = 1,VECTOR_SIZE_LOC
                         do igaus = 1,pgaus
                            do inode = 1,qnode
                               kposi = (igaus - 1)*(inode-1) * qnode*ndime
                               do idime = 1, ndime
                                  kposi = kposi + 1
                                  xvalu(ivect,kposi) = prope_ker % gdval_ielem(ielem(ivect)) % a(idime,inode,igaus)
                               end do
                            end do
                         end do
                      end do
                   end if
                else
                   do kposi = 1,pgaus*qnode*ndime
                      xvalu(1:VECTOR_SIZE_LOC,kposi) = 0.0_rp
                   end do
                end if

             else if( kfl_deriv_tur == 1 ) then

                if( prope_ker % llaws(ilaws) % kfl_deriv_tur == 1 ) then

                   if( present(gpsha) ) then
                      nturb = 2
                      if( qgaus /= pgaus ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                      do ivect = 1,VECTOR_SIZE_LOC
                         do igaus = 1,qgaus
                            kposi = (igaus - 1)* qnode*ndime
                            do inode = 1,qnode
                               do iturb = 1, nturb
                                  kposi = kposi + 1
                                  xvalu(ivect,kposi) = prope_ker % drval_tur_ielem(ielem(ivect)) % a(iturb,inode,igaus)
                               end do
                            end do
                         end do
                      end do
                   else
                      do ivect = 1,VECTOR_SIZE_LOC
                         do igaus = 1,pgaus
                            do inode = 1,qnode
                               kposi = (igaus - 1)*(inode-1) * qnode*ndime
                               do iturb = 1, nturb
                                  kposi = kposi + 1
                                  xvalu(ivect,kposi) = prope_ker % drval_tur_ielem(ielem(ivect)) % a(iturb,inode,igaus)
                               end do
                            end do
                         end do
                      end do
                   end if
                else
                   do kposi = 1,pgaus*qnode*nturb
                      xvalu(1:VECTOR_SIZE_LOC,kposi) = 0.0_rp
                   end do
                end if

             else if( kfl_deriv_vel == 1 ) then

                if( prope_ker % llaws(ilaws) % kfl_deriv_vel == 1 ) then

                   if( present(gpsha) ) then
                      if( qgaus /= pgaus ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                      do ivect = 1,VECTOR_SIZE_LOC
                         do igaus = 1,qgaus
                            kposi = (igaus - 1)* qnode*ndime
                            do inode = 1,qnode
                               do idime = 1, ndime
                                  kposi = kposi + 1
                                  xvalu(ivect,kposi) = prope_ker % drval_vel_ielem(ielem(ivect)) % a(idime,inode,igaus)
                               end do
                            end do
                         end do
                      end do
                   else
                      do ivect = 1,VECTOR_SIZE_LOC
                         do igaus = 1,pgaus
                            do inode = 1,qnode
                               kposi = (igaus - 1)*(inode-1) * qnode*ndime
                               do idime = 1, ndime
                                  kposi = kposi + 1
                                  xvalu(ivect,kposi) = prope_ker % drval_vel_ielem(ielem(ivect)) % a(idime,inode,igaus)
                               enddo
                            end do
                         end do
                      end do
                   end if
                else
                   do kposi = 1,pgaus*qnode*ndime
                      xvalu(1:VECTOR_SIZE_LOC,kposi) = 0.0_rp
                   end do
                end if

             else                
                if( present(gpsha) .and. kfl_cutel == 0 ) then
                   if( qgaus /= pgaus ) call runend('MOD_KER_PROPER: PGAUS NOT CODED')
                 
!@@                   lpalloc = acc_is_present(xvalu)
!@@                   if( .not. lpalloc) then 
!@@                        !$acc enter data copyin(xvalu(1:VECTOR_SIZE_LOC,1:qgaus))
!@@                   end if 
!@@
!@@                   lpalloc_ind = acc_is_present(ielbo)
!@@                   if(lpalloc_ind == .false.) then 
!@@                        !$acc enter data copyin(ielbo(1:VECTOR_SIZE_LOC))
!@@                   end if 
!@@
!@@                   lpalloc_stru = acc_is_present(prope_ker%value_ielem)
!@@                   if(lpalloc_stru == .false.) then 
!@@                        !$acc enter data copyin(prope_ker)
!@@                        !$acc enter data copyin(prope_ker%value_ielem)
!@@                        do ivect = 1,VECTOR_SIZE_LOC
!@@                            !$acc enter data copyin(prope_ker%value_ielem(ielbo(ivect))%a)
!@@                        end do
!@@                 
!@@                   end if 
!@@                   !$acc parallel loop gang vector default(present) 
                   do ivect = 1,VECTOR_SIZE_LOC
                      do igaus = 1,qgaus * prope_ker % dim
                         xvalu(ivect,igaus) = prope_ker % value_ielem(ielbo(ivect)) % a(igaus)
                        !@ xvalu(ivect,igaus) = prope_ker % value_ielem(ielem(ivect)) % a(igaus)                    
                      end do
                   end do
!@@                   !$acc end parallel loop
!@@
!@@                   if( .not. lpalloc) then 
!@@                        !$acc update self (xvalu(1:VECTOR_SIZE_LOC,1:qgaus))
!@@                        !$acc exit data delete(xvalu(1:VECTOR_SIZE_LOC,1:qgaus))
!@@                   end if 
!@@ 
!@@                   if(lpalloc_ind == .false.) then 
!@@                        !$acc exit data delete(ielbo(1:VECTOR_SIZE_LOC))
!@@                   end if 
!@@
!@@                   if(lpalloc_stru == .false.) then 
!@@                        !$acc exit data delete(prope_ker)
!@@                        !$acc exit data delete(prope_ker%value_ielem)
!@@                        do ivect = 1,VECTOR_SIZE_LOC
!@@                            !$acc exit data delete(prope_ker%value_ielem(ielbo(ivect))%a)
!@@                        end do
!@@ 
!@@                   end if 

                else
                   do ivect = 1,VECTOR_SIZE_LOC
                      do igaus = 1,pgaus * prope_ker % dim
                         xvalu(ivect,igaus) = prope_ker % value_ielem(ielem(ivect)) % a(igaus)
                      end do
                   end do
                end if
             end if

          else if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then
             !
             ! Property is computed node-wise
             ! Note: it has not been validated 03/03/2017
             do ivect = 1,VECTOR_SIZE_LOC
                pelty = abs(ltype(ielem(ivect)))
                pnode = nnode(pelty)
                do inode = 1,pnode
                   ipoin = lnods(inode,ielem(ivect))
                   elpro(ivect,inode) = prope_ker % value_ipoin(ipoin)
                end do
             end do

             if( kfl_gradi == 1 ) then
                !
                ! Gradient
                !
                do kposi = 1,pgaus*ndime
                   xvalu(1:VECTOR_SIZE_LOC,kposi) = 0.0_rp
                end do

                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then

                   do ivect = 1,VECTOR_SIZE_LOC
                      do inode = 1,pnode
                         ipoin = lnods(inode,ielem(ivect))
                         do idime = 1,ndime
                            elgrp(ivect,idime,inode) = prope_ker % grval_ipoin(idime,ipoin)
                         end do
                      end do
                   end do

                   if( present(gpsha) ) then
                      do igaus = 1,qgaus
                         kposi = (igaus-1) * qnode
                         do inode = 1,pnode
                            kposi = kposi + 1
                            lposi = (igaus-1) * ndime
                            do idime = 1,ndime
                               lposi = lposi + 1
                               xvalu(1:VECTOR_SIZE_LOC,lposi) = xvalu(1:VECTOR_SIZE_LOC,lposi) + elgrp(1:VECTOR_SIZE_LOC,idime,inode) * gpsha(1:VECTOR_SIZE_LOC,kposi)
                            end do
                         end do
                      end do
                   else
                      do igaus = 1,qgaus
                         do inode = 1,pnode
                            lposi = (igaus-1) * ndime
                            do idime = 1,ndime
                               lposi = lposi + 1
                               xvalu(1:VECTOR_SIZE_LOC,lposi) = xvalu(1:VECTOR_SIZE_LOC,lposi) + elgrp(1:VECTOR_SIZE_LOC,idime,inode) * elmar(pelty) % shape(inode,igaus)
                            end do
                         end do
                      end do
                   end if

                else

                   if( present(gpcar) ) then
                      do igaus = 1,qgaus
                         kposi = (igaus-1) * ndime
                         do idime = 1,ndime
                            kposi = kposi + 1
                            do inode = 1,pnode
                               xvalu(1:VECTOR_SIZE_LOC,kposi) = xvalu(1:VECTOR_SIZE_LOC,kposi) + elpro(1:VECTOR_SIZE_LOC,inode) * gpcar(1:VECTOR_SIZE_LOC,idime,inode,igaus)
                            end do
                         end do
                      end do
                   end if

                end if

             else
                !
                ! Value
                !
                xvalu(1:VECTOR_SIZE_LOC,1:qgaus) = 0.0_rp

                if( present(gpsha) ) then
                   do igaus = 1,qgaus
                      kposi        = (igaus-1) * qnode
                      do inode = 1,pnode
                         kposi = kposi + 1
                         xvalu(1:VECTOR_SIZE_LOC,igaus) = xvalu(1:VECTOR_SIZE_LOC,igaus) + elpro(1:VECTOR_SIZE_LOC,inode) * gpsha(1:VECTOR_SIZE_LOC,kposi)
                      end do
                   end do
                else
                   do igaus = 1,ngaus(pelty)
                      do inode = 1,pnode
                         xvalu(1:VECTOR_SIZE_LOC,igaus) = xvalu(1:VECTOR_SIZE_LOC,igaus) + elpro(1:VECTOR_SIZE_LOC,inode) * elmar(pelty) % shape(inode,igaus)
                      end do
                   end do
                end if

             end if

          end if

          !----------------------------------------------------------------------
          !
          ! 7. Property required on C.O.G.
          !    Note: it has not been validated 03/03/2017
          !----------------------------------------------------------------------

       else if( wherein(1:3) == 'COG' ) then

          ielem(1:VECTOR_SIZE_LOC) = ielbo(1:VECTOR_SIZE_LOC)
          imate = lmate(ielem(1))
          ilaws = prope_ker % ilaws(imate)
          pelty = abs(ltype(ielem(1)))
          pgaus = ngaus(pelty)
          pnode = nnode(pelty)

          if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then

             if( kfl_gradi == 1 )  then
                !
                ! Gradient
                !
                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then
                   xvalu(1:VECTOR_SIZE_LOC,1:ndime) = 0.0_rp
                   do ivect = 1,VECTOR_SIZE_LOC
                      if( ielem(ivect) > 0 ) then
                         do igaus = 1,pgaus
                            do idime = 1,ndime
                               xvalu(ivect,idime) = xvalu(ivect,idime) + prope_ker % grval_ielem(ielem(ivect)) % a(idime,igaus)
                            end do
                         end do
                      end if
                   end do
                   do idime = 1,ndime
                      do ivect = 1,VECTOR_SIZE_LOC
                         xvalu(ivect,idime) = xvalu(ivect,idime) / real(pgaus,rp)
                      end do
                   end do
                else
                   call runend('MOD_KER_PROPER: GRDAIENT ON COG NOT CODED')
                end if

             else
                !
                ! Value
                !
                if (ielem(1) < 0 ) call runend('ielem < 0_rp should not have happened ')   ! before version 13823  there was an abs(ielem(1)  when calculating xvalu
                ! Guillaume said that should not happen. I leave this runend just in case - It can be removed after some time    
                xvalu(1:VECTOR_SIZE_LOC,1) = 0.0_rp
                if( present(gpsha) ) then
                   do ivect = 1,VECTOR_SIZE_LOC
                      do igaus = 1,qgaus
                         xvalu(ivect,1) = xvalu(ivect,1) + prope_ker % value_ielem(ielem(ivect)) % a(igaus)
                      end do
                   end do
                   xvalu(1:VECTOR_SIZE_LOC,1) = xvalu(1:VECTOR_SIZE_LOC,1) / real(qgaus,rp)
                else
                   do ivect = 1,VECTOR_SIZE_LOC
                      do igaus = 1,pgaus
                         xvalu(ivect,1) = xvalu(ivect,1) + prope_ker % value_ielem(ielem(ivect)) % a(igaus)
                      end do
                   end do
                   xvalu(1:VECTOR_SIZE_LOC,1) = xvalu(1:VECTOR_SIZE_LOC,1) / real(pgaus,rp)
                end if
             end if

          else if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

             if( kfl_gradi == 1 )  then
                !
                ! Gradient
                !
                xvalu(1:VECTOR_SIZE_LOC,1:ndime) = 0.0_rp

                if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then

                   do ivect = 1,VECTOR_SIZE_LOC
                      if( ielem(ivect) > 0 ) then
                         do inode = 1,pnode
                            ipoin = lnods(inode,ielem(ivect))
                            do idime = 1,ndime
                               xvalu(ivect,idime) = xvalu(ivect,idime) + prope_ker % grval_ipoin(idime,ipoin)
                            end do
                         end do
                      end if
                   end do
                   do idime = 1,ndime
                      xvalu(1:VECTOR_SIZE_LOC,idime) = xvalu(1:VECTOR_SIZE_LOC,idime) / real(pnode,rp)
                   end do

                else

                   if( present(gpcar) ) then
                      do ivect = 1,VECTOR_SIZE_LOC
                         if( ielem(ivect) > 0 ) then
                            do igaus = 1,qgaus
                               do inode = 1,pnode
                                  ipoin = lnods(inode,ielem(ivect))
                                  do idime = 1,ndime
                                     xvalu(ivect,idime) = xvalu(ivect,idime) &
                                          + prope_ker % value_ipoin(ipoin) * gpcar(ivect,idime,inode,1)
                                  end do
                               end do
                            end do
                         end if
                      end do
                      do idime = 1,ndime
                         xvalu(1:VECTOR_SIZE_LOC,idime) = xvalu(1:VECTOR_SIZE_LOC,idime) / real(qgaus,rp)
                      end do
                   else
                      call runend('MOD_KER_PROPER: GRADIENT NOT AVAILABLE AT COG')
                   end if

                end if

             else
                !
                ! Value
                !
                xvalu(1:VECTOR_SIZE_LOC,1) = 0.0_rp
                do ivect = 1,VECTOR_SIZE_LOC
                   if( ielem(ivect) > 0 ) then
                      do inode = 1,pnode
                         ipoin          = lnods(inode,ielem(ivect))
                         xvalu(ivect,1) = xvalu(ivect,1) + prope_ker % value_ipoin(ipoin)
                      end do
                   end if
                end do
                xvalu(1:VECTOR_SIZE_LOC,1) = xvalu(1:VECTOR_SIZE_LOC,1) / real(pnode,rp)
             end if

          else if( prope_ker % llaws(ilaws) % where == 'CONST' ) then

             if( kfl_gradi == 1 )  then
                xvalu(1:VECTOR_SIZE_LOC,1:ndime) = 0.0_rp
             else
                xvalu(1:VECTOR_SIZE_LOC,1) = prope_ker % value_const(1,imate)
             end if

          end if

       else

          call runend('KER_PROPER_VECTOR: NOT CODED 3')

       end if

    end if
    !
    ! This fails with OpenMP
    !nullify(prope_ker) ! This avoids the bug of calling wrongly the ker_proper routine

    !
    ! This part deallocates the prope_ker from the GPU
#ifdef OPENACCHHH
    if (associated(prope_ker) .and. associated(prope_ker % value_const)) then
     !$acc exit data delete(prope_ker, prope_ker%value_const)
    else if (associated(prope_ker)) then
     !$acc exit data delete(prope_ker)
    end if
#endif

  end subroutine ker_proper_vector

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-28
  !> @brief   Check if a property has just been updated
  !> @details Check if a property has just been updated
  !> 
  !-----------------------------------------------------------------------

  logical(lg) function ker_proper_updated_property(wprop,itask)
    character(len=*),            intent(in) :: wprop
    integer(ip),                 intent(in) :: itask
    integer(ip)                             :: iprop,imate,ilaws,iresp
    type(typ_valpr_ker), pointer            :: prope_ker

    prope_ker => ker_proper_pointer(wprop)
    
    ker_proper_updated_property = .false.
    if( itask == ITASK_INIUNK ) then
       if( prope_ker % kfl_exist == 1 ) then
          do imate = 1,nmate
             if( prope_ker % wlaws(imate) == 'CONST' ) &
                  ker_proper_updated_property = .true.
          end do
       end if
    end if

    do iprop = 1,NUMBER_OF_PROPERTIES
       if( prope_ker % kfl_exist == 1 ) then
          do imate = 1,nmate
             if( prope_ker % wlaws(imate) /= 'CONST' ) then
                ilaws = prope_ker % ilaws(imate)
                do iresp = 1,mresp_ker
                   if( prope_ker % llaws(ilaws) % lresp(iresp) == lastm_ker&
                        .and.(prope_ker % update(1, imate)==itask&
                        .or.prope_ker % update(2, imate)==itask)) then
                      ker_proper_updated_property = .true.
                   else if ( prope_ker % llaws(ilaws) % lresp(iresp) == -1 ) then
                      ker_proper_updated_property = .true.
                   end if
                end do
             end if
          end do
       end if
    end do

  end function ker_proper_updated_property

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-11-12
  !> @brief   Parallelization of properties
  !> @details Brodcast of properties in arallel
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_proper_parall()

    integer(ip)                  :: iprop
    type(typ_valpr_ker), pointer :: prope_ker

    do iprop = 1,NUMBER_OF_PROPERTIES
       prope_ker => ker_proper_pointer(iprop)
       call exchange_add(prope_ker % kfl_exist)
       call exchange_add(prope_ker % kfl_type)
       call exchange_add(prope_ker % rlaws)
       call exchange_add(prope_ker % value_default)
       call exchange_add(prope_ker % wlaws)
       call exchange_add(prope_ker % comp)
       call exchange_add(prope_ker % on_the_fly)
       call exchange_add(prope_ker % update)
       call exchange_add(prope_ker % time_function)
    end do

  end subroutine ker_proper_parall

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-01-02
  !> @brief   Return if a property exists
  !> @details retrun if a properties have been created
  !> 
  !-----------------------------------------------------------------------

  logical(lg) function ker_proper_exists(wprop) 

    character(len=*),            intent(in) :: wprop
    type(typ_valpr_ker), pointer            :: prope_ker

    prope_ker => ker_proper_pointer(wprop)

    if( prope_ker % kfl_exist == 0 ) then
       ker_proper_exists = .false.
    else
       ker_proper_exists = .true.
    end if

  end function ker_proper_exists

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-01-02
  !> @brief   Check properties
  !> @details Check that properties have beed defined according to the
  !>          Alya modules activated. Check that responsible modules
  !>          are also on.
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_proper_check_properties()

    integer(ip)                   :: imodu,imate,iprop,iresp,ilaws
    type(typ_valpr_ker), pointer  :: prope_ker
    !
    ! Required properties
    !
    do imodu = 1,mmodu
       if( kfl_modul(imodu) /= 0 ) then
          select case ( imodu )

          case ( ID_NASTIN )

             if( .not. ker_proper_exists('VISCO') ) call runend(namod(imodu)//': VISCOSITY SHOULD BE DEFINED AS A PROPERTY IN *.KER.DAT FILE')
             if( .not. ker_proper_exists('DENSI') ) call runend(namod(imodu)//': DENSITY SHOULD BE DEFINED AS A PROPERTY IN *.KER.DAT FILE')

          case ( ID_TEMPER )

             if( .not. ker_proper_exists('SPHEA') ) call runend(namod(imodu)//': SPECIFIC HEAT SHOULD BE DEFINED AS A PROPERTY IN *.KER.DAT FILE')
             if( .not. ker_proper_exists('CONDU') ) call runend(namod(imodu)//': CONDUCTIVITY SHOULD BE DEFINED AS A PROPERTY IN *.KER.DAT FILE')
             if( .not. ker_proper_exists('DENSI') ) call runend(namod(imodu)//': DENSITY SHOULD BE DEFINED AS A PROPERTY IN *.KER.DAT FILE')

          case ( ID_PARTIS )

             if( .not. ker_proper_exists('VISCO') ) call runend(namod(imodu)//': VISCOSITY SHOULD BE DEFINED AS A PROPERTY IN *.KER.DAT FILE')
             if( .not. ker_proper_exists('DENSI') ) call runend(namod(imodu)//': DENSITY SHOULD BE DEFINED AS A PROPERTY IN *.KER.DAT FILE')

          end select
       end if
    end do
    !
    ! Check if responsible is on
    !
    do iprop = 1,NUMBER_OF_PROPERTIES
       prope_ker => ker_proper_pointer(iprop)
       if( prope_ker % kfl_exist == 1 ) then
          do imate = 1,nmate
             ilaws = prope_ker % ilaws(imate)
             do iresp = 1,size(prope_ker % llaws(ilaws) % lresp)
                imodu = prope_ker % llaws(ilaws) % lresp(iresp)
                if( imodu > 0 ) then
                   if( kfl_tefun == 0 .and. kfl_vefun == 0 .and. kfl_cofun == 0 ) then 
                      if( kfl_modul(imodu) == 0 ) then
                         call runend(&
                              namod(imodu)//' MODULE IS OFF AND IS RESPONSIBLE FOR MODEL '//&
                              trim(prope_ker % llaws(ilaws) % wname)//' OF PROPERTY '//&
                              trim(prope_ker % name))
                      end if
                   end if
                end if
             end do
          end do
       end if
    end do

  end subroutine ker_proper_check_properties

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-01-02
  !> @brief   Number of properties
  !> @details Number of existing properties 
  !> 
  !-----------------------------------------------------------------------
  
  function ker_proper_number_properties() result(num)
    integer(ip)                   :: num
    type(typ_valpr_ker), pointer  :: prope_ker
    integer(ip)                   :: iprop
    
    num = 0
    do iprop = 1,NUMBER_OF_PROPERTIES
       prope_ker => ker_proper_pointer(iprop)
       num = num + prope_ker % kfl_exist
    end do
    
  end function ker_proper_number_properties
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-01-02
  !> @brief   Number of properties
  !> @details Number of existing properties 
  !> 
  !-----------------------------------------------------------------------
  
  subroutine ker_proper_on_the_fly(kfl_on_the_fly) 
    integer(ip),                  intent(in) :: kfl_on_the_fly
    type(typ_valpr_ker), pointer             :: prope_ker
    integer(ip)                              :: iprop

    if( kfl_on_the_fly /= 0 ) then
       do iprop = 1,NUMBER_OF_PROPERTIES
          prope_ker => ker_proper_pointer(iprop)
          if( prope_ker % kfl_exist /= 0 ) prope_ker % on_the_fly = max(0_ip,kfl_on_the_fly)
       end do
    end if
    
  end subroutine ker_proper_on_the_fly
  
end module mod_ker_proper
!> @}
