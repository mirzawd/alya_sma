!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Master
!> @{
!> @file    mod_arrays.f90
!> @author  houzeaux
!> @date    2019-11-13
!> @brief   Manipulate arrays
!> @details Manipulate arrays:
!>          1. Register
!>          2. Allocate
!>          3. Deallocate
!>          4. Redistribute
!>          5. Write in restart file
!>          6. Read in restart file
!-----------------------------------------------------------------------

module mod_arrays

  use def_master
  use def_kintyp_basic,      only : ip,rp,lg,r3p
  use def_kintyp_mesh_basic, only : mesh_type_basic
  use def_kermod,            only : witness_mesh
  use def_domain,            only : nelem,nboun,npoin
  use mod_memory_basic,      only : memory_alloca
  use mod_memory_basic,      only : memory_deallo
  use mod_memory_basic,      only : memory_size
  use mod_memory_basic,      only : memory_copy
  use mod_memory_basic,      only : memory_copy_r3p_1
  use mod_redistribution,    only : redistribution_array
  use mod_postpr,            only : postpr_read_restart
  use mod_postpr,            only : postpr_write_restart
  use mod_postpr,            only : postpr_postprocess
  use mod_communications,    only : PAR_MAX
  use mod_messages,          only : messages_live
  use mod_AMR_interpolate,   only : AMR_interpolate
  use mod_optional_argument, only : optional_argument
  use mod_strings,           only : integer_to_string
  
  implicit none

  character(10),  parameter :: vacal = 'mod_arrays'
  integer(8)                :: memor_loc(2)
  integer(ip)               :: enti_posit_loc
  integer(ip)               :: nenti_loc
  character(200)            :: variable_name_loc
  character(200)            :: subroutine_name_loc
  integer(ip)               :: modul_loc
  character(5)              :: wopos_loc(5)
  integer(ip)               :: kfl_reawr_loc
  integer(ip)               :: minmem
  logical(lg),    parameter :: verbose = .true.
  integer(8)                :: memor_arrays(2)
  integer(ip)               :: array_total_number(mmodu)
  integer(ip)               :: array_counter(mmodu)
  !
  ! Primary variables register
  !
  type arrays_all_types
     integer(ip)                     :: ivari
     real(rp),           pointer     :: xr1(:)
     integer(ip),        pointer     :: xi1(:)
     logical(lg),        pointer     :: xl1(:)
     type(r3p),          pointer     :: xp1(:)
     real(rp),           pointer     :: xr2(:,:)
     integer(ip),        pointer     :: xi2(:,:)
     logical(lg),        pointer     :: xl2(:,:)
     real(rp),           pointer     :: xr3(:,:,:)
     integer(ip),        pointer     :: xi3(:,:,:)
     logical(lg),        pointer     :: xl3(:,:,:)
   contains
     procedure,          pass        :: init => init_arrays_all_types
  end type arrays_all_types
  type arrays_reset_module
     type(arrays_all_types), pointer :: a(:)
  end type arrays_reset_module
  integer(ip)                        :: num_primary_variables(0:mmodu)
  integer(ip)                        :: num_secondary_variables(0:mmodu)
  type(arrays_reset_module)          :: xx_reg(mmodu)
  type(arrays_reset_module)          :: xx_sav(mmodu)
  
  private

  interface arrays
     module procedure &
          &           arrays_RP_1,arrays_RP_2,arrays_RP_3,arrays_R3P_1,&
          &           arrays_IP_1
  end interface arrays

  interface arrays_used
     module procedure &
          &           arrays_used_RP_1,arrays_used_RP_2,&
          &           arrays_used_RP_3,arrays_used_R3P_1
  end interface arrays_used

  interface arrays_register
     module procedure &
          &           arrays_register_0,&
          &           arrays_register_RP_1,arrays_register_RP_2,&
          &           arrays_register_RP_3,arrays_register_R3P_1,&
          &           arrays_register_IP_1,&
          &           arrays_register_IP_2,&
          &           arrays_register_0b,&
          &           arrays_register_RP_1b,arrays_register_RP_2b,&
          &           arrays_register_RP_3b,arrays_register_R3P_1b,&
          &           arrays_register_IP_1b,&
          &           arrays_register_IP_2b
  end interface arrays_register

  interface arrays_save_pointer
     module procedure &
          arrays_save_pointer_1,&
          arrays_save_pointer_2,&
          arrays_save_pointer_3
  end interface arrays_save_pointer
  
  public :: arrays
  public :: arrays_initialization
  public :: arrays_register          ! Register a variable
  public :: arrays_number
  public :: arrays_name
  public :: arrays_allocated
  public :: arrays_tag               ! wopos of a module variable
  public :: arrays_allocate_primary  ! Allocate memory for reset
  public :: arrays_save_primary      ! Save primary variables
  public :: arrays_recover_primary   ! Recover primary variables
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-13
  !> @brief   Initialize arrays
  !> @details Initialize arrays
  !> 
  !-----------------------------------------------------------------------
  
  subroutine init_arrays_all_types(self)
    class(arrays_all_types), intent(inout) :: self

     nullify(self % xr1)
     nullify(self % xi1)
     nullify(self % xl1)
     nullify(self % xp1)
     nullify(self % xr2)
     nullify(self % xi2)
     nullify(self % xl2)
     nullify(self % xr3)
     nullify(self % xi3)
     nullify(self % xl3)
    
  end subroutine init_arrays_all_types

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-13
  !> @brief   Free position
  !> @details Find free position 
  !> 
  !-----------------------------------------------------------------------
  
  function arrays_primary_position(MODULE_NUMBER) result(ivari)
    
    integer(ip), optional, intent(in) :: MODULE_NUMBER
    integer(ip)                       :: ivari,ii,imodu
    
    ivari = 0
    imodu = optional_argument(modul,MODULE_NUMBER)    

    if( associated(xx_reg(imodu) % a) ) then 
       imodu = optional_argument(modul,MODULE_NUMBER)
       loop_ii: do ii = 1,size(xx_reg(imodu) % a,KIND=ip)
          if( xx_reg(imodu) % a(ii) % ivari == 0 ) then
             ivari = ii
             exit loop_ii
          end if
       end do loop_ii
    end if
   
  end function arrays_primary_position

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-13
  !> @brief   Allocate memory for reset
  !> @details Allocate memory for reset
  !> 
  !-----------------------------------------------------------------------
  
  subroutine arrays_save_pointer_1(ivari,xx,MODULE_NUMBER)

    integer(ip),                    intent(in) :: ivari
    class(*),              target,  intent(in) :: xx(:)
    integer(ip), optional,          intent(in) :: MODULE_NUMBER
    integer(ip)                                :: kvari,imodu
    
    imodu = optional_argument(modul,MODULE_NUMBER)

    if( momod(imodu) % postp(1) % wopos(4,ivari) == 'PRIMA' ) then
       
       call arrays_allocate_primary(MODULE_NUMBER)
       kvari = arrays_primary_position(imodu)
       if( kvari == 0 ) call runend('ARRAYS_SAVE_POINTER: ARRAY POSITION NOT FOUND X(:)')

       if( verbose ) call messages_live(namod(imodu)//': REGISTERING PRIMARY VARIABLE '//momod(imodu) % postp(1) % wopos(1,ivari))
     
       xx_reg(imodu) % a(kvari) % ivari = ivari
       
       select type ( xx )
          
       type is ( integer (kind=ip ) ) ; xx_reg(imodu) % a(kvari) % xi1 => xx
       type is ( real    (kind=rp ) ) ; xx_reg(imodu) % a(kvari) % xr1 => xx       
       type is ( logical (kind=lg ) ) ; xx_reg(imodu) % a(kvari) % xl1 => xx
       type is ( r3p                ) ; xx_reg(imodu) % a(kvari) % xp1 => xx
          
       end select

    end if
    
  end subroutine arrays_save_pointer_1
  
  subroutine arrays_save_pointer_2(ivari,xx,MODULE_NUMBER)

    integer(ip),                    intent(in) :: ivari
    class(*),              target,  intent(in) :: xx(:,:)
    integer(ip), optional,          intent(in) :: MODULE_NUMBER
    integer(ip)                                :: kvari,imodu

    imodu = optional_argument(modul,MODULE_NUMBER)
    
    if( momod(imodu) % postp(1) % wopos(4,ivari) == 'PRIMA' ) then

       call arrays_allocate_primary(MODULE_NUMBER)
       kvari = arrays_primary_position(imodu)
       if( kvari == 0 ) call runend('ARRAYS_SAVE_POINTER: ARRAY POSITION NOT FOUND X(:,:)')
       
       if( verbose ) call messages_live(namod(imodu)//': REGISTERING PRIMARY VARIABLE '//momod(imodu) % postp(1) % wopos(1,ivari))

       xx_reg(imodu) % a(kvari) % ivari = ivari
       
       select type ( xx )
          
       type is ( integer (kind=ip) ) ; xx_reg(imodu) % a(kvari) % xi2 => xx
       type is ( real    (kind=rp) ) ; xx_reg(imodu) % a(kvari) % xr2 => xx       
       type is ( logical (kind=lg) ) ; xx_reg(imodu) % a(kvari) % xl2 => xx
          
       end select

    end if
    
  end subroutine arrays_save_pointer_2
  
  subroutine arrays_save_pointer_3(ivari,xx,MODULE_NUMBER)

    integer(ip),                    intent(in) :: ivari
    class(*),              target,  intent(in) :: xx(:,:,:)
    integer(ip), optional,          intent(in) :: MODULE_NUMBER
    integer(ip)                                :: kvari,imodu

    imodu = optional_argument(modul,MODULE_NUMBER)
    
    if( momod(imodu) % postp(1) % wopos(4,ivari) == 'PRIMA' ) then
       
       call arrays_allocate_primary(MODULE_NUMBER)
       kvari = arrays_primary_position(imodu)
       if( kvari == 0 ) call runend('ARRAYS_SAVE_POINTER: ARRAY POSITION NOT FOUND X(:,:,:)')
       
       if( verbose ) call messages_live(namod(imodu)//': REGISTERING PRIMARY VARIABLE '//momod(imodu) % postp(1) % wopos(1,ivari))

       xx_reg(imodu) % a(kvari) % ivari = ivari
       
       select type ( xx )
          
       type is ( integer (kind=ip) ) ; xx_reg(imodu) % a(kvari) % xi3 => xx
       type is ( real    (kind=rp) ) ; xx_reg(imodu) % a(kvari) % xr3 => xx       
       type is ( logical (kind=lg) ) ; xx_reg(imodu) % a(kvari) % xl3 => xx
          
       end select

    end if
    
  end subroutine arrays_save_pointer_3
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-13
  !> @brief   Allocate memory for reset
  !> @details Allocate memory for reset
  !> 
  !-----------------------------------------------------------------------
  
  subroutine arrays_allocate_primary(MODULE_NUMBER)

    integer(ip), optional, intent(in) :: MODULE_NUMBER
    integer(ip)                       :: imodu,nvari,ivari

    imodu = optional_argument(modul,MODULE_NUMBER)

    !if( verbose ) call messages_live('REGISTERING PRIMARY VARIABLES','START SECTION')
    if( .not. associated(xx_reg(imodu) % a) ) then
       nvari = num_primary_variables(imodu) !arrays_number_primary_arrays(imodu)
       !call messages_live(namod(imodu)//': '//integer_to_string(nvari)//' POSSIBLE VARIABLES IN TOTAL')
       if( nvari > 0 ) then
          allocate(xx_reg(imodu) % a(nvari))
          allocate(xx_sav(imodu) % a(nvari))
          do ivari = 1,nvari
             xx_reg(imodu) % a(ivari) % ivari = 0
             xx_sav(imodu) % a(ivari) % ivari = 0
             call xx_reg(imodu) % a(ivari) % init()
             call xx_sav(imodu) % a(ivari) % init()
          end do
       end if
    end if
    !if( verbose ) call messages_live('REGISTERING PRIMARY VARIABLES','END SECTION')
    
  end subroutine arrays_allocate_primary
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-13
  !> @brief   Save primary variables
  !> @details Save primary variables for further use
  !> 
  !-----------------------------------------------------------------------
  
  subroutine arrays_save_primary()

    integer(ip) :: imodu,nvari,kvari,ivari

    do imodu = 1,mmodu
       if( kfl_modul(imodu) /= 0 .and. associated( xx_reg(imodu) % a) ) then
          nvari = arrays_number_primary_arrays(imodu)
          if( nvari > 0 ) then
             do kvari = 1,nvari
                ivari = xx_reg(imodu) % a(kvari) % ivari
                if( ivari /= 0 ) then
                   !print*,'saving=',momod(imodu) % postp(1) % wopos(1,ivari)
                   if(      associated(xx_reg(imodu) % a(kvari) % xr1) ) then
                      call memory_copy      (mem_modul(1:2,imodu),'XR1',vacal,xx_reg(imodu) % a(kvari) % xr1,xx_sav(imodu) % a(kvari) % xr1,'DO_NOT_DEALLOCATE')
                   else if( associated(xx_reg(imodu) % a(kvari) % xi1) ) then
                      call memory_copy      (mem_modul(1:2,imodu),'XI1',vacal,xx_reg(imodu) % a(kvari) % xi1,xx_sav(imodu) % a(kvari) % xi1,'DO_NOT_DEALLOCATE')
                   else if( associated(xx_reg(imodu) % a(kvari) % xp1) ) then
                      call memory_copy_r3p_1(mem_modul(1:2,imodu),'XP1',vacal,xx_reg(imodu) % a(kvari) % xp1,xx_sav(imodu) % a(kvari) % xp1,'DO_NOT_DEALLOCATE')
                   else if( associated(xx_reg(imodu) % a(kvari) % xr2) ) then
                      call memory_copy      (mem_modul(1:2,imodu),'XR2',vacal,xx_reg(imodu) % a(kvari) % xr2,xx_sav(imodu) % a(kvari) % xr2,'DO_NOT_DEALLOCATE')
                   else if( associated(xx_reg(imodu) % a(kvari) % xi2) ) then
                      call memory_copy      (mem_modul(1:2,imodu),'XI2',vacal,xx_reg(imodu) % a(kvari) % xi2,xx_sav(imodu) % a(kvari) % xi2,'DO_NOT_DEALLOCATE')
                   else if( associated(xx_reg(imodu) % a(kvari) % xr3) ) then
                      call memory_copy      (mem_modul(1:2,imodu),'XR3',vacal,xx_reg(imodu) % a(kvari) % xr3,xx_sav(imodu) % a(kvari) % xr3,'DO_NOT_DEALLOCATE')
                   else if( associated(xx_reg(imodu) % a(kvari) % xi3) ) then
                      call memory_copy      (mem_modul(1:2,imodu),'XI3',vacal,xx_reg(imodu) % a(kvari) % xi3,xx_sav(imodu) % a(kvari) % xi3,'DO_NOT_DEALLOCATE')                   
                   end if
                end if
             end do
          end if
       end if
    end do


  end subroutine arrays_save_primary
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-13
  !> @brief   Recover primary
  !> @details Recover primary variables from last saved
  !> 
  !-----------------------------------------------------------------------
  
  subroutine arrays_recover_primary()

    integer(ip) :: imodu,nvari,kvari,ivari

    do imodu = 1,mmodu
       if( kfl_modul(imodu) /= 0 .and. associated( xx_reg(imodu) % a) ) then
          nvari = arrays_number_primary_arrays(imodu) 
          if( nvari > 0 ) then
             do kvari = 1,nvari
                ivari = xx_reg(imodu) % a(kvari) % ivari
                if( ivari /= 0 ) then
                   !print*,'recover=',ivari,momod(imodu) % postp(1) % wopos(1,ivari)
                   if(      associated(xx_reg(imodu) % a(kvari) % xr1) ) then
                      call memory_copy      (mem_modul(1:2,imodu),'XR1',vacal,xx_sav(imodu) % a(kvari) % xr1,xx_reg(imodu) % a(kvari) % xr1,'DO_NOT_DEALLOCATE')
                   else if( associated(xx_reg(imodu) % a(kvari) % xi1) ) then
                      call memory_copy      (mem_modul(1:2,imodu),'XI1',vacal,xx_sav(imodu) % a(kvari) % xi1,xx_reg(imodu) % a(kvari) % xi1,'DO_NOT_DEALLOCATE')
                   else if( associated(xx_reg(imodu) % a(kvari) % xp1) ) then
                      call memory_copy_r3p_1(mem_modul(1:2,imodu),'XP1',vacal,xx_sav(imodu) % a(kvari) % xp1,xx_reg(imodu) % a(kvari) % xp1,'DO_NOT_DEALLOCATE')
                   else if( associated(xx_reg(imodu) % a(kvari) % xr2) ) then
                      call memory_copy      (mem_modul(1:2,imodu),'XR2',vacal,xx_sav(imodu) % a(kvari) % xr2,xx_reg(imodu) % a(kvari) % xr2,'DO_NOT_DEALLOCATE')
                   else if( associated(xx_reg(imodu) % a(kvari) % xi2) ) then
                      call memory_copy      (mem_modul(1:2,imodu),'XI2',vacal,xx_sav(imodu) % a(kvari) % xi2,xx_reg(imodu) % a(kvari) % xi2,'DO_NOT_DEALLOCATE')
                   else if( associated(xx_reg(imodu) % a(kvari) % xr3) ) then
                      call memory_copy      (mem_modul(1:2,imodu),'XR3',vacal,xx_sav(imodu) % a(kvari) % xr3,xx_reg(imodu) % a(kvari) % xr3,'DO_NOT_DEALLOCATE')
                   else if( associated(xx_reg(imodu) % a(kvari) % xi3) ) then
                      call memory_copy      (mem_modul(1:2,imodu),'XI3',vacal,xx_sav(imodu) % a(kvari) % xi3,xx_reg(imodu) % a(kvari) % xi3,'DO_NOT_DEALLOCATE')
                   end if
                end if
             end do
          end if
       end if
    end do

  end subroutine arrays_recover_primary
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-13
  !> @brief   Number of primary variables
  !> @details Number of primary variables of a module
  !> 
  !-----------------------------------------------------------------------
  
  integer(ip) function arrays_number_primary_arrays(MODULE_NUMBER) result(nn)

    integer(ip), optional, intent(in) :: MODULE_NUMBER
    integer(ip)                       :: imodu,ivari

    imodu = optional_argument(modul,MODULE_NUMBER)

    if( imodu > 0 ) then
       if( num_primary_variables(imodu) == 0 ) then
          nn = 0     
          do ivari = 1,int(size(postp(1) % wopos,2),ip)
             if(  momod(imodu) % postp(1) % wopos (4,ivari) == 'PRIMA' .and. &
                  momod(imodu) % postp(1) % wopos (1,ivari) /= 'NULL ' ) nn = nn + 1
          end do
          num_primary_variables(imodu) = nn
       else
          nn = num_primary_variables(imodu) 
       end if
    end if
    
  end function arrays_number_primary_arrays
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-13
  !> @brief   Variable name
  !> @details Find the name of a variable given its number
  !> 
  !-----------------------------------------------------------------------
  
  pure function arrays_name(ivari,MODULE_NUMBER) result(wname)

    integer(ip),                   intent(in)    :: ivari
    integer(ip),         optional, intent(in)    :: MODULE_NUMBER
    character(len=5)                             :: wname
    integer(ip)                                  :: imodu

    if( ivari <= 0 ) then
       wname = 'NULL '
    else
       imodu = optional_argument(modul,MODULE_NUMBER)
       wname = momod(imodu) % postp(1) % wopos (1,ivari)
    end if
    
  end function arrays_name
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-13
  !> @brief   Variable number
  !> @details Find the number of a variable
  !> 
  !-----------------------------------------------------------------------
  
  integer(ip) pure function arrays_number(wname,MODULE_NUMBER)

    character(len=*),              intent(in)    :: wname
    integer(ip),         optional, intent(in)    :: MODULE_NUMBER
    integer(ip)                                  :: ivari,imodu
    
    imodu         = optional_argument(modul,MODULE_NUMBER)
    arrays_number = 0_ip
    do ivari = 1,nvarp
       if( momod(imodu) % postp(1) % wopos (1,ivari) == trim(wname) ) then
          arrays_number = ivari
          return
       end if
    end do
       
  end function arrays_number
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-13
  !> @brief   Initialization
  !> @details Initialization of the module
  !> 
  !-----------------------------------------------------------------------
  
  subroutine arrays_initialization()

    integer(ip) :: imodu
    
    memor_arrays       = 0_8
    array_total_number = 0
    array_counter      = 0
    minmem             = 0

    num_primary_variables   = 0
    num_secondary_variables = 0
    do imodu = 1,mmodu
       nullify(xx_reg(imodu) % a)
       nullify(xx_sav(imodu) % a)
    end do
    
  end subroutine arrays_initialization
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-13
  !> @brief   Register arrays
  !> @details Register an array:
  !>
  !>          ARRAY_NUMBER ...... Number of the array, must be unique
  !>          WOPOS(1) .......... Name 
  !>          WOPOS(2) .......... Dimension=                  SCALA/VECTO/MATR/R3PVE
  !>          WOPOS(3) .......... Entity=                     NPOIN/NELEM/NBOUN
  !>          WOPOS(4) .......... Primary or secondary array= PRIMA/SECON
  !>          WOPOS(5) .......... Kind of the array=          REA/INTEG
  !>
  !>          ENTITY_POSITION ... Where the entity dimension is located
  !>          TIME_POSITION ..... For transient variables, where time is
  !>
  !>          Primary variables are those required for a restart, a repartitioning
  !>          or a remeshing. Secondary are basically used for postprocess
  !>
  !-----------------------------------------------------------------------
  
  subroutine arrays_register_0(ARRAY_NUMBER,wopos,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,POST_ONLY_ONCE)

    integer(ip),                   intent(in)    :: ARRAY_NUMBER
    character(len=*),              intent(in)    :: wopos(:)
    integer(ip),         optional, intent(in)    :: ENTITY_POSITION
    integer(ip),         optional, intent(in)    :: TIME_POSITION
    integer(ip),         optional, intent(in)    :: COMPONENT_POSITION
    integer(ip),         optional, intent(in)    :: MODULE_NUMBER
    logical(lg),         optional, intent(in)    :: POST_ONLY_ONCE

    call arrays_register_go(wopos,ARRAY_NUMBER,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,POST_ONLY_ONCE)

  end subroutine arrays_register_0
  
  subroutine arrays_register_IP_1(ARRAY_NUMBER,wopos,xx,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,POST_ONLY_ONCE)

    integer(ip),                   intent(in)    :: ARRAY_NUMBER
    character(len=*),              intent(in)    :: wopos(:)
    integer(ip),         pointer,  intent(inout) :: xx(:)
    integer(ip),         optional, intent(in)    :: ENTITY_POSITION
    integer(ip),         optional, intent(in)    :: TIME_POSITION
    integer(ip),         optional, intent(in)    :: COMPONENT_POSITION
    integer(ip),         optional, intent(in)    :: MODULE_NUMBER
    logical(lg),         optional, intent(in)    :: POST_ONLY_ONCE

    nullify(xx)
    call arrays_register_go(wopos,ARRAY_NUMBER,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,VAR_INT=.true.,POST_ONLY_ONCE=POST_ONLY_ONCE)
    
  end subroutine arrays_register_IP_1
  
  subroutine arrays_register_IP_2(ARRAY_NUMBER,wopos,xx,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,POST_ONLY_ONCE)

    integer(ip),                   intent(in)    :: ARRAY_NUMBER
    character(len=*),              intent(in)    :: wopos(:)
    integer(ip),         pointer,  intent(inout) :: xx(:,:)
    integer(ip),         optional, intent(in)    :: ENTITY_POSITION
    integer(ip),         optional, intent(in)    :: TIME_POSITION
    integer(ip),         optional, intent(in)    :: COMPONENT_POSITION
    integer(ip),         optional, intent(in)    :: MODULE_NUMBER
    logical(lg),         optional, intent(in)    :: POST_ONLY_ONCE

    if( present(ENTITY_POSITION) ) then
       if( ENTITY_POSITION > 2 ) call runend('ARRAYS_REGISTER: WRONG POSITION')
    end if
    nullify(xx)
    call arrays_register_go(wopos,ARRAY_NUMBER,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,VAR_INT=.true.,POST_ONLY_ONCE=POST_ONLY_ONCE)
    !
    ! Guess component position
    !
    if( .not. present(COMPONENT_POSITION) ) then
       if( wopos(2) == 'SCALA' ) then
          if( .not. present(TIME_POSITION) ) then
             if( momod(modul_loc) % postp(1) % time_posit(ARRAY_NUMBER) == 2 ) then
                momod(modul_loc) % postp(1) % comp_posit(ARRAY_NUMBER) = 1
             else if( momod(modul_loc) % postp(1) % time_posit(ARRAY_NUMBER) == 1 ) then
                momod(modul_loc) % postp(1) % comp_posit(ARRAY_NUMBER) = 2
             end if
          end if
       end if
    end if
    
  end subroutine arrays_register_IP_2
  
  subroutine arrays_register_RP_1(ARRAY_NUMBER,wopos,xx,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,POST_ONLY_ONCE)

    integer(ip),                   intent(in)    :: ARRAY_NUMBER
    character(len=*),              intent(in)    :: wopos(:)
    real(rp),            pointer,  intent(inout) :: xx(:)
    integer(ip),         optional, intent(in)    :: ENTITY_POSITION
    integer(ip),         optional, intent(in)    :: TIME_POSITION
    integer(ip),         optional, intent(in)    :: COMPONENT_POSITION
    integer(ip),         optional, intent(in)    :: MODULE_NUMBER
    logical(lg),         optional, intent(in)    :: POST_ONLY_ONCE

    nullify(xx)
    call arrays_register_go(wopos,ARRAY_NUMBER,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,POST_ONLY_ONCE)
    
  end subroutine arrays_register_RP_1
  
  subroutine arrays_register_RP_2(ARRAY_NUMBER,wopos,xx,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,EXCLUDE_RST,POST_ONLY_ONCE)

    integer(ip),                             intent(in)    :: ARRAY_NUMBER
    character(len=*),                        intent(in)    :: wopos(:)
    real(rp),                      pointer,  intent(inout) :: xx(:,:)
    integer(ip),         optional,           intent(in)    :: ENTITY_POSITION
    integer(ip),         optional,           intent(in)    :: TIME_POSITION
    integer(ip),         optional,           intent(in)    :: COMPONENT_POSITION
    integer(ip),         optional,           intent(in)    :: EXCLUDE_RST(:)
    integer(ip),         optional,           intent(in)    :: MODULE_NUMBER
     logical(lg),        optional,           intent(in)    :: POST_ONLY_ONCE
   
    if( present(ENTITY_POSITION) ) then
       if( ENTITY_POSITION > 2 ) call runend('ARRAYS_REGISTER: WRONG POSITION')
    end if
    nullify(xx)
    call arrays_register_go(wopos,ARRAY_NUMBER,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,EXCLUDE_RST=EXCLUDE_RST,POST_ONLY_ONCE=POST_ONLY_ONCE)
    !
    ! Guess component position
    !
    if( .not. present(COMPONENT_POSITION) ) then
       if( wopos(2) == 'SCALA' ) then
          if( .not. present(TIME_POSITION) ) then
             if( momod(modul_loc) % postp(1) % time_posit(ARRAY_NUMBER) == 2 ) then
                momod(modul_loc) % postp(1) % comp_posit(ARRAY_NUMBER) = 1
             else if( momod(modul_loc) % postp(1) % time_posit(ARRAY_NUMBER) == 1 ) then
                momod(modul_loc) % postp(1) % comp_posit(ARRAY_NUMBER) = 2
             end if
          end if
       end if
    end if
    
  end subroutine arrays_register_RP_2
  
  subroutine arrays_register_RP_3(ARRAY_NUMBER,wopos,xx,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,EXCLUDE_RST,POST_ONLY_ONCE)

    integer(ip),                           intent(in)    :: ARRAY_NUMBER
    character(len=*),                      intent(in)    :: wopos(:)
    real(rp),                    pointer,  intent(inout) :: xx(:,:,:)
    integer(ip),       optional,           intent(in)    :: ENTITY_POSITION
    integer(ip),       optional,           intent(in)    :: TIME_POSITION
    integer(ip),       optional,           intent(in)    :: COMPONENT_POSITION
    integer(ip),       optional,           intent(in)    :: MODULE_NUMBER
    integer(ip),       optional,           intent(in)    :: EXCLUDE_RST(:)
    logical(lg),       optional,           intent(in)    :: POST_ONLY_ONCE

    if( present(ENTITY_POSITION) ) then
       if( ENTITY_POSITION > 3 ) call runend('ARRAYS_REGISTER: WRONG POSITION')
    end if

    nullify(xx)
    call arrays_register_go(wopos,ARRAY_NUMBER,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,EXCLUDE_RST=EXCLUDE_RST,POST_ONLY_ONCE=POST_ONLY_ONCE)
    !
    ! Guess component position
    !
    if( .not. present(COMPONENT_POSITION) ) then
       if( wopos(2) == 'SCALA' ) then
          if( present(ENTITY_POSITION) ) then
             if( momod(modul_loc) % postp(1) % time_posit(ARRAY_NUMBER) == 3 ) then
                if(      ENTITY_POSITION == 2 ) then
                   momod(modul_loc) % postp(1) % comp_posit(ARRAY_NUMBER) = 1
                else if( ENTITY_POSITION == 1 ) then
                   momod(modul_loc) % postp(1) % comp_posit(ARRAY_NUMBER) = 2
                end if
             else
                call runend('ARRAYS_REGISTER: CANNOT GUESS')
             end if
          end if
       end if
    end if

  end subroutine arrays_register_RP_3
  
  subroutine arrays_register_R3P_1(ARRAY_NUMBER,wopos,xx,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,EXCLUDE_RST,POST_ONLY_ONCE)

    integer(ip),                             intent(in)    :: ARRAY_NUMBER
    character(len=*),                        intent(in)    :: wopos(:)
    type(r3p),                     pointer,  intent(inout) :: xx(:)
    integer(ip),         optional,           intent(in)    :: ENTITY_POSITION
    integer(ip),         optional,           intent(in)    :: MODULE_NUMBER
    integer(ip),         optional,           intent(in)    :: TIME_POSITION
    integer(ip),         optional,           intent(in)    :: COMPONENT_POSITION
    integer(ip),         optional,           intent(in)    :: EXCLUDE_RST(:)
    logical(lg),         optional,           intent(in)    :: POST_ONLY_ONCE

    if( present(ENTITY_POSITION) ) then
       if( ENTITY_POSITION > 1 ) call runend('ARRAYS_REGISTER: WRONG POSITION')
    end if
    nullify(xx)
    call arrays_register_go(wopos,ARRAY_NUMBER,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,EXCLUDE_RST=EXCLUDE_RST,POST_ONLY_ONCE=POST_ONLY_ONCE)
    
  end subroutine arrays_register_R3P_1

  subroutine arrays_register_0b(wopos,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,POST_ONLY_ONCE)

    
    character(len=*),              intent(in)    :: wopos(:)
    integer(ip),         optional, intent(in)    :: ENTITY_POSITION
    integer(ip),         optional, intent(in)    :: TIME_POSITION
    integer(ip),         optional, intent(in)    :: COMPONENT_POSITION
    integer(ip),         optional, intent(in)    :: MODULE_NUMBER
    integer(ip)                                  :: imodu
    logical(lg),         optional, intent(in)    :: POST_ONLY_ONCE
    
    imodu = optional_argument(modul,MODULE_NUMBER)
    array_counter(imodu) = array_counter(imodu) + 1
    call arrays_register_go(wopos,array_counter(modul),ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,POST_ONLY_ONCE=POST_ONLY_ONCE)

  end subroutine arrays_register_0b
  
  subroutine arrays_register_IP_1b(wopos,xx,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,POST_ONLY_ONCE)

    
    character(len=*),              intent(in)    :: wopos(:)
    integer(ip),         pointer,  intent(inout) :: xx(:)
    integer(ip),         optional, intent(in)    :: ENTITY_POSITION
    integer(ip),         optional, intent(in)    :: TIME_POSITION
    integer(ip),         optional, intent(in)    :: COMPONENT_POSITION
    integer(ip),         optional, intent(in)    :: MODULE_NUMBER
    integer(ip)                                  :: imodu
    logical(lg),         optional, intent(in)    :: POST_ONLY_ONCE
   
    imodu = optional_argument(modul,MODULE_NUMBER)
    array_counter(imodu) = array_counter(imodu) + 1

    nullify(xx)
    call arrays_register_go(wopos,array_counter(modul),ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,VAR_INT=.true.,POST_ONLY_ONCE=POST_ONLY_ONCE)
    
  end subroutine arrays_register_IP_1b
  
  subroutine arrays_register_IP_2b(wopos,xx,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,POST_ONLY_ONCE)

    
    character(len=*),              intent(in)    :: wopos(:)
    integer(ip),         pointer,  intent(inout) :: xx(:,:)
    integer(ip),         optional, intent(in)    :: ENTITY_POSITION
    integer(ip),         optional, intent(in)    :: TIME_POSITION
    integer(ip),         optional, intent(in)    :: COMPONENT_POSITION
    integer(ip),         optional, intent(in)    :: MODULE_NUMBER
    integer(ip)                                  :: imodu,iarra
    logical(lg),         optional, intent(in)    :: POST_ONLY_ONCE
    
    imodu = optional_argument(modul,MODULE_NUMBER)
    array_counter(imodu) = array_counter(imodu) + 1
    iarra = array_counter(imodu)
    
    if( present(ENTITY_POSITION) ) then
       if( ENTITY_POSITION > 2 ) call runend('ARRAYS_REGISTER: WRONG POSITION')
    end if
    nullify(xx)
    call arrays_register_go(wopos,array_counter(modul),ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,VAR_INT=.true.,POST_ONLY_ONCE=POST_ONLY_ONCE)
    !
    ! Guess component position
    !
    if( .not. present(COMPONENT_POSITION) ) then
       if( wopos(2) == 'SCALA' ) then
          if( .not. present(TIME_POSITION) ) then
             if( momod(modul_loc) % postp(1) % time_posit(iarra) == 2 ) then
                momod(modul_loc) % postp(1) % comp_posit(iarra) = 1
             else if( momod(modul_loc) % postp(1) % time_posit(iarra) == 1 ) then
                momod(modul_loc) % postp(1) % comp_posit(iarra) = 2
             end if
          end if
       end if
    end if
    
  end subroutine arrays_register_IP_2b
  
  subroutine arrays_register_RP_1b(wopos,xx,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,POST_ONLY_ONCE)

    
    character(len=*),              intent(in)    :: wopos(:)
    real(rp),            pointer,  intent(inout) :: xx(:)
    integer(ip),         optional, intent(in)    :: ENTITY_POSITION
    integer(ip),         optional, intent(in)    :: TIME_POSITION
    integer(ip),         optional, intent(in)    :: COMPONENT_POSITION
    integer(ip),         optional, intent(in)    :: MODULE_NUMBER
    logical(lg),         optional, intent(in)    :: POST_ONLY_ONCE
    integer(ip)                                  :: imodu
    
    imodu = optional_argument(modul,MODULE_NUMBER)
    array_counter(imodu) = array_counter(imodu) + 1

    nullify(xx)
    call arrays_register_go(wopos,array_counter(modul),ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,POST_ONLY_ONCE=POST_ONLY_ONCE)
    
  end subroutine arrays_register_RP_1b
  
  subroutine arrays_register_RP_2b(wopos,xx,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,EXCLUDE_RST,POST_ONLY_ONCE)

    character(len=*),                        intent(in)    :: wopos(:)
    real(rp),                      pointer,  intent(inout) :: xx(:,:)
    integer(ip),         optional,           intent(in)    :: ENTITY_POSITION
    integer(ip),         optional,           intent(in)    :: TIME_POSITION
    integer(ip),         optional,           intent(in)    :: COMPONENT_POSITION
    integer(ip),         optional,           intent(in)    :: EXCLUDE_RST(:)
    integer(ip),         optional,           intent(in)    :: MODULE_NUMBER
    logical(lg),         optional,           intent(in)    :: POST_ONLY_ONCE
    integer(ip)                                            :: imodu,iarra
   
    imodu = optional_argument(modul,MODULE_NUMBER)
    array_counter(imodu) = array_counter(imodu) + 1
    iarra = array_counter(imodu)
    
    if( present(ENTITY_POSITION) ) then
       if( ENTITY_POSITION > 2 ) call runend('ARRAYS_REGISTER: WRONG POSITION')
    end if
    nullify(xx)
    call arrays_register_go(wopos,array_counter(modul),ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,EXCLUDE_RST=EXCLUDE_RST,POST_ONLY_ONCE=POST_ONLY_ONCE)
    !
    ! Guess component position
    !
    if( .not. present(COMPONENT_POSITION) ) then
       if( wopos(2) == 'SCALA' ) then
          if( .not. present(TIME_POSITION) ) then
             if( momod(modul_loc) % postp(1) % time_posit(iarra) == 2 ) then
                momod(modul_loc) % postp(1) % comp_posit(iarra) = 1
             else if( momod(modul_loc) % postp(1) % time_posit(iarra) == 1 ) then
                momod(modul_loc) % postp(1) % comp_posit(iarra) = 2
             end if
          end if
       end if
    end if
    
  end subroutine arrays_register_RP_2b
  
  subroutine arrays_register_RP_3b(wopos,xx,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,EXCLUDE_RST,POST_ONLY_ONCE)

    character(len=*),                      intent(in)    :: wopos(:)
    real(rp),                    pointer,  intent(inout) :: xx(:,:,:)
    integer(ip),       optional,           intent(in)    :: ENTITY_POSITION
    integer(ip),       optional,           intent(in)    :: TIME_POSITION
    integer(ip),       optional,           intent(in)    :: COMPONENT_POSITION
    integer(ip),       optional,           intent(in)    :: MODULE_NUMBER
    integer(ip),       optional,           intent(in)    :: EXCLUDE_RST(:)
    logical(lg),       optional,           intent(in)    :: POST_ONLY_ONCE
    integer(ip)                                          :: imodu,iarra
    
    imodu = optional_argument(modul,MODULE_NUMBER)
    array_counter(imodu) = array_counter(imodu) + 1
    iarra = array_counter(imodu)

    if( present(ENTITY_POSITION) ) then
       if( ENTITY_POSITION > 3 ) call runend('ARRAYS_REGISTER: WRONG POSITION')
    end if

    nullify(xx)
    call arrays_register_go(wopos,array_counter(modul),ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,EXCLUDE_RST=EXCLUDE_RST,POST_ONLY_ONCE=POST_ONLY_ONCE)
    !
    ! Guess component position
    !
    if( .not. present(COMPONENT_POSITION) ) then
       if( wopos(2) == 'SCALA' ) then
          if( present(ENTITY_POSITION) ) then
             if( momod(modul_loc) % postp(1) % time_posit(iarra) == 3 ) then
                if(      ENTITY_POSITION == 2 ) then
                   momod(modul_loc) % postp(1) % comp_posit(iarra) = 1
                else if( ENTITY_POSITION == 1 ) then
                   momod(modul_loc) % postp(1) % comp_posit(iarra) = 2
                end if
             else
                call runend('ARRAYS_REGISTER: CANNOT GUESS')
             end if
          end if
       end if
    end if

  end subroutine arrays_register_RP_3b
  
  subroutine arrays_register_R3P_1b(wopos,xx,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,EXCLUDE_RST,POST_ONLY_ONCE)

    character(len=*),                        intent(in)    :: wopos(:)
    type(r3p),                     pointer,  intent(inout) :: xx(:)
    integer(ip),         optional,           intent(in)    :: ENTITY_POSITION
    integer(ip),         optional,           intent(in)    :: MODULE_NUMBER
    integer(ip),         optional,           intent(in)    :: TIME_POSITION
    integer(ip),         optional,           intent(in)    :: COMPONENT_POSITION
    integer(ip),         optional,           intent(in)    :: EXCLUDE_RST(:)
    integer(ip)                                            :: imodu
    logical(lg),         optional,           intent(in)    :: POST_ONLY_ONCE
    
    imodu = optional_argument(modul,MODULE_NUMBER)
    array_counter(imodu) = array_counter(imodu) + 1

    if( present(ENTITY_POSITION) ) then
       if( ENTITY_POSITION > 1 ) call runend('ARRAYS_REGISTER: WRONG POSITION')
    end if
    nullify(xx)
    call arrays_register_go(wopos,array_counter(modul),ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,EXCLUDE_RST=EXCLUDE_RST,POST_ONLY_ONCE=POST_ONLY_ONCE)
    
  end subroutine arrays_register_R3P_1b
 
  subroutine arrays_register_go(wopos,ARRAY_NUMBER,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,VAR_INT,EXCLUDE_RST,POST_ONLY_ONCE)

    character(len=*),              intent(in)    :: wopos(:)
    integer(ip),         optional, intent(in)    :: ARRAY_NUMBER
    integer(ip),         optional, intent(in)    :: ENTITY_POSITION
    integer(ip),         optional, intent(in)    :: TIME_POSITION
    integer(ip),         optional, intent(in)    :: COMPONENT_POSITION
    integer(ip),         optional, intent(in)    :: MODULE_NUMBER
    logical(lg),         optional, intent(in)    :: VAR_INT
    integer(ip),         optional, intent(in)    :: EXCLUDE_RST(:)
    logical(lg),         optional, intent(in)    :: POST_ONLY_ONCE
    integer(ip)                                  :: ii,itime,ivari
    logical(lg)                                  :: if_var_int,only_once

    modul_loc = optional_argument(modul,MODULE_NUMBER)
    only_once = optional_argument(.false.,POST_ONLY_ONCE)
    
    if( present(ARRAY_NUMBER) ) then
       ivari = ARRAY_NUMBER
    else
       array_total_number(modul_loc) = array_total_number(modul_loc) + 1
       ivari = array_total_number(modul_loc)
    end if
    if( present(VAR_INT) ) then
       if_var_int = VAR_INT
    else
       if_var_int = .false.       
    end if
    if( present(EXCLUDE_RST) ) then
       do ii = 1,size(EXCLUDE_RST,KIND=ip)
          itime = EXCLUDE_RST(ii)
          momod(modul_loc) % postp(1) % rst_time(itime,ivari) = .false. 
       end do
    end if
    
    if( momod(modul_loc) % postp(1) % array_registered(ivari) == 1 ) then
       call runend('ARRAY_REGISTER: ARRAY '//trim(wopos(1))//' HAS ALREADY BEEN REGISTERED')
    else
       momod(modul_loc) % postp(1) % array_registered(ivari) = 1
       momod(modul_loc) % postp(1) % wopos (1,ivari)         = trim(wopos(1))
       momod(modul_loc) % postp(1) % wopos (2,ivari)         = trim(wopos(2))
       momod(modul_loc) % postp(1) % wopos (3,ivari)         = trim(wopos(3))
       momod(modul_loc) % postp(1) % wopos (4,ivari)         = trim(wopos(4))
       if( size(wopos) < 5 ) then
          if( if_var_int ) then
             momod(modul_loc) % postp(1) % wopos (5,ivari)      = 'INT'
          else
             momod(modul_loc) % postp(1) % wopos (5,ivari)      = 'REAL'
          end if
       else
          momod(modul_loc) % postp(1) % wopos (5,ivari)      = trim(wopos(5))
       end if
       
       if( present(ENTITY_POSITION) ) then
          momod(modul_loc) % postp(1) % enti_posit (ivari) = ENTITY_POSITION
       else       
          momod(modul_loc) % postp(1) % enti_posit (ivari) = 1
       end if
       
       if( present(TIME_POSITION) ) then
          momod(modul_loc) % postp(1) % time_posit (ivari) = TIME_POSITION
       else       
          momod(modul_loc) % postp(1) % time_posit (ivari) = 0
       end if
       
       if( present(COMPONENT_POSITION) ) then
          momod(modul_loc) % postp(1) % comp_posit (ivari) = COMPONENT_POSITION
       else       
          momod(modul_loc) % postp(1) % comp_posit (ivari) = 0
       end if
    end if

    if( only_once ) then
       momod(modul_loc) % postp(1) % kfl_oonce(ivari) = 1
    else
       momod(modul_loc) % postp(1) % kfl_oonce(ivari) = 0
    end if
    
    if( momod(modul_loc) % postp(1) % wopos (4,ivari) == 'PRIMA' ) then
       num_primary_variables(modul_loc) = num_primary_variables(modul_loc) + 1
    else if( momod(modul_loc) % postp(1) % wopos (4,ivari) == 'SECON' ) then
       num_secondary_variables(modul_loc) = num_secondary_variables(modul_loc) + 1
    else
       call runend('ARRAYS_REGISTER_GO: VARIABLE KIND PRIMA/SECON NOT CORRECT')
    end if
    
  end subroutine arrays_register_go

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-13
  !> @brief   Do what we have to do with arrays
  !> @details Do what we have to do with arrays
  !> 
  !-----------------------------------------------------------------------

  subroutine arrays_RP_1(ivari,wtask,xx,ndim1,ENTITY_POSITION,TIME_POSITION,TAG1,FORCE,MESH,MESH_ID)

    integer(ip),                                intent(in)    :: ivari
    character(len=*),                           intent(in)    :: wtask
    real(rp),                         pointer,  intent(inout) :: xx(:)
    integer(ip),            optional,           intent(in)    :: ndim1
    integer(ip),            optional,           intent(in)    :: ENTITY_POSITION
    integer(ip),            optional,           intent(in)    :: TIME_POSITION
    integer(ip),            optional,           intent(in)    :: TAG1
    logical(lg),            optional,           intent(in)    :: FORCE
    type(mesh_type_basic),  optional,           intent(in)    :: MESH
    integer(ip),            optional,           intent(in)    :: MESH_ID
    
    call arrays_options_start(ivari)
    
    select case ( trim(wtask) )
       
    case ( 'ALLOCATE' )
       !
       ! Allocate
       !
       if( present(ndim1) ) then
          if( enti_posit_loc == 1 .and. ndim1 < nenti_loc ) call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' NOT REGISTERED CORRECTLY')
          if( enti_posit_loc /= 1 )                          call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' NOT REGISTERED CORRECTLY')
          !if( trim(wopos_loc(4)) /= 'PRIMA' )           call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' SHOULD BE REGISTERED AS PRIMARY')
          call memory_alloca(memor_loc,trim(variable_name_loc),trim(subroutine_name_loc),xx,max(minmem,ndim1))
          if( associated(xx) ) call arrays_save_pointer(ivari,xx,modul)
          postp(1) % array_used(ivari) = 1
          if( memory_size(xx) > 0 ) postp(1) % array_allocated(ivari) = 1 
       else
          call runend('ARRAYS: DIMENSIONS ARE MISSING')
       end if
       
    case ( 'DEALLOCATE' )
       !
       ! deallocate
       !
       call memory_deallo(memor_loc,trim(variable_name_loc),trim(subroutine_name_loc),xx)
       postp(1) % array_used(ivari) = 0
       postp(1) % array_allocated(ivari) = 0
       
    case ( 'REDISTRIBUTE' )
       !
       ! Redistribute
       !
       call redistribution_array(xx,wopos_loc(3),POSIT=enti_posit_loc,MEMOR=memor_loc,VARIABLE_NAME=variable_name_loc)
       
    case ( 'INTERPOLATE' )
       ! 
       ! Interpolate
       !
       call AMR_interpolate(xx,wopos_loc(3),POSIT=enti_posit_loc,MEMOR=memor_loc,VARIABLE_NAME=variable_name_loc)
       
    case ( 'WRITE RESTART' )
       !
       ! Write restart
       !
       if( postp(1) % array_used(ivari) == 1 ) then
          call messages_live('VARIABLE: '//wopos_loc(1)//' FROM '//namod(modul),'START SECTION')
          call postpr_write_restart(xx,ivari,ittim,cutim,TAG1=TAG1)   
          call messages_live('VARIABLE','END SECTION')
       end if
         
    case ( 'READ RESTART' )
       !
       ! Write restart
       !  
       if( postp(1) % array_used(ivari) == 1 ) then
          call messages_live('VARIABLE: '//wopos_loc(1)//' FROM '//namod(modul),'START SECTION')
          call postpr_read_restart(xx,ivari,ittim,cutim,TAG1=TAG1)   
          call messages_live('VARIABLE','END SECTION')
      end if
       
    case ( 'POSTPROCESS' )
       !
       ! Postprocess
       !  
       if( postp(1) % array_used(ivari) == 1 .or. optional_argument(.false.,FORCE) ) then
          if( present(MESH) ) then
             call postpr_postprocess(xx,ivari,ittim,cutim,MESH=MESH)             
          else if( present(MESH_ID) ) then
             if(     MESH_ID == 0 ) then
                call postpr_postprocess(xx,ivari,ittim,cutim)
             else if( MESH_ID > 0 ) then
                call postpr_postprocess(xx,ivari,ittim,cutim,MESH=witness_mesh(MESH_ID) % mesh)
             end if
          else
             call postpr_postprocess(xx,ivari,ittim,cutim)
          end if
      end if
       
    case default
       !
       ! Do not know what to do
       !
       call runend('ARRAYS: DO NOT KNOW WHAT TO DO '//trim(wtask)//' WITH '//trim(wopos_loc(1)))
       
    end select

    call arrays_options_end()

  end subroutine arrays_RP_1
  
  subroutine arrays_RP_2(ivari,wtask,xx,ndim1,ndim2,ENTITY_POSITION,TIME_POSITION,TAG1,FORCE,MESH,MESH_ID)

    integer(ip),                                intent(in)    :: ivari
    character(len=*),                           intent(in)    :: wtask
    real(rp),                         pointer,  intent(inout) :: xx(:,:)
    integer(ip),            optional,           intent(in)    :: ndim1
    integer(ip),            optional,           intent(in)    :: ndim2
    integer(ip),            optional,           intent(in)    :: ENTITY_POSITION
    integer(ip),            optional,           intent(in)    :: TIME_POSITION
    integer(ip),            optional,           intent(in)    :: TAG1
    logical(lg),            optional,           intent(in)    :: FORCE
    type(mesh_type_basic),  optional,           intent(in)    :: MESH
    integer(ip),            optional,           intent(in)    :: MESH_ID
    
    call arrays_options_start(ivari)

    select case ( trim(wtask) )

    case ( 'ALLOCATE' )
       !
       ! Allocate
       !
       if( present(ndim1) .and. present(ndim2) ) then
          if( enti_posit_loc == 1 .and. ndim1 < nenti_loc ) call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' NOT REGISTERED CORRECTLY')
          if( enti_posit_loc == 2 .and. ndim2 < nenti_loc ) call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' NOT REGISTERED CORRECTLY')
          !if( trim(wopos_loc(4)) /= 'PRIMA' )           call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' SHOULD BE REGISTERED AS PRIMARY')
          call memory_alloca(memor_loc,trim(variable_name_loc),trim(subroutine_name_loc),xx,max(minmem,ndim1),max(minmem,ndim2))
          if( associated(xx) ) call arrays_save_pointer(ivari,xx,modul)
          postp(1) % array_used(ivari) = 1
          if( memory_size(xx) > 0 ) postp(1) % array_allocated(ivari) = 1 
          if(      postp(1) % time_posit(ivari) == 1 ) then
             postp(1) % time_num(ivari) = ndim1
          else if( postp(1) % time_posit(ivari) == 2 ) then
             postp(1) % time_num(ivari) = ndim2
          end if
          if(      postp(1) % comp_posit(ivari) == 1 ) then
             postp(1) % comp_num(ivari) = ndim1
          else if( postp(1) % comp_posit(ivari) == 2 ) then
             postp(1) % comp_num(ivari) = ndim2
          end if
          if( postp(1) % wopos(2,ivari) == 'VECTO' ) then
             if(      postp(1) % comp_posit(ivari) /= 1 .and. postp(1) % time_posit(ivari) /= 1 ) then
                postp(1) % dime_num(ivari) = ndim1
             else 
                postp(1) % dime_num(ivari) = ndim2
             end if
          end if
          
       else
          call runend('ARRAYS: DIMENSIONS ARE MISSING')
       end if
       
    case ( 'DEALLOCATE' )
       !
       ! deallocate
       !
       call memory_deallo(memor_loc,trim(variable_name_loc),trim(subroutine_name_loc),xx)
       postp(1) % array_used(ivari) = 0
       postp(1) % array_allocated(ivari) = 0
       
    case ( 'REDISTRIBUTE' )
       !
       ! Redistribute
       !
       call redistribution_array(xx,wopos_loc(3),POSIT=enti_posit_loc,MEMOR=memor_loc,VARIABLE_NAME=variable_name_loc)

    case ( 'INTERPOLATE' )
       ! 
       ! Interpolate
       !
       call AMR_interpolate(xx,wopos_loc(3),POSIT=enti_posit_loc,MEMOR=memor_loc,VARIABLE_NAME=variable_name_loc)
       
    case ( 'WRITE RESTART' )
       !
       ! Write restart
       !
       if( postp(1) % array_used(ivari) == 1 ) then
          call messages_live('VARIABLE: '//wopos_loc(1)//' FROM '//namod(modul),'START SECTION')
          call postpr_write_restart(xx,ivari,ittim,cutim,posit=enti_posit_loc,TAG1=TAG1)   
          call messages_live('VARIABLE','END SECTION')
       end if
       
    case ( 'READ RESTART' )
       !
       ! Read restart
       !  
       if( postp(1) % array_used(ivari) == 1 ) then
          call messages_live('VARIABLE: '//wopos_loc(1)//' FROM '//namod(modul),'START SECTION')
          call postpr_read_restart(xx,ivari,ittim,cutim,posit=enti_posit_loc,TAG1=TAG1)   
          call messages_live('VARIABLE','END SECTION')
       end if
       
    case ( 'POSTPROCESS' )
       !
       ! Postprocess
       !
       if( postp(1) % array_used(ivari) == 1 .or. optional_argument(.false.,FORCE) ) then
          if( present(MESH) ) then
             call postpr_postprocess(xx,ivari,ittim,cutim,MESH=MESH)             
          else if( present(MESH_ID) ) then
             if(     MESH_ID == 0 ) then
                call postpr_postprocess(xx,ivari,ittim,cutim)
             else if( MESH_ID > 0 ) then
                call postpr_postprocess(xx,ivari,ittim,cutim,MESH=witness_mesh(MESH_ID) % mesh)
             end if
          else
             call postpr_postprocess(xx,ivari,ittim,cutim)
          end if
      end if

    case default
       !
       ! Do not know what to do
       !
       call runend('ARRAYS: DO NOT KNOW WHAT TO DO '//trim(wtask)//' WITH '//trim(wopos_loc(1)))
       
    end select

    call arrays_options_end()

  end subroutine arrays_RP_2
  
  subroutine arrays_RP_3(ivari,wtask,xx,ndim1,ndim2,ndim3,ENTITY_POSITION,TIME_POSITION,TAG1,FORCE,MESH,MESH_ID)

    integer(ip),                                intent(in)    :: ivari    
    character(len=*),                           intent(in)    :: wtask
    real(rp),                         pointer,  intent(inout) :: xx(:,:,:)
    integer(ip),            optional,           intent(in)    :: ndim1
    integer(ip),            optional,           intent(in)    :: ndim2
    integer(ip),            optional,           intent(in)    :: ndim3
    integer(ip),            optional,           intent(in)    :: ENTITY_POSITION
    integer(ip),            optional,           intent(in)    :: TIME_POSITION
    integer(ip),            optional,           intent(in)    :: TAG1
    logical(lg),            optional,           intent(in)    :: FORCE
    type(mesh_type_basic),  optional,           intent(in)    :: MESH
    integer(ip),            optional,           intent(in)    :: MESH_ID

    call arrays_options_start(ivari)

    select case ( trim(wtask) )

    case ( 'ALLOCATE' )
       !
       ! Allocate
       !
       if( present(ndim1) .and. present(ndim2) .and. present(ndim3) ) then
          if( enti_posit_loc == 1 .and. ndim1 /= nenti_loc ) call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' NOT REGISTERED CORRECTLY')
          if( enti_posit_loc == 2 .and. ndim2 /= nenti_loc ) call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' NOT REGISTERED CORRECTLY')
          if( enti_posit_loc == 3 .and. ndim3 /= nenti_loc ) call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' NOT REGISTERED CORRECTLY')
          !if( trim(wopos_loc(4)) /= 'PRIMA' )           call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' SHOULD BE REGISTERED AS PRIMARY')
          call memory_alloca(memor_loc,trim(variable_name_loc),trim(subroutine_name_loc),xx,max(minmem,ndim1),max(minmem,ndim2),max(minmem,ndim3))
          if( associated(xx) ) call arrays_save_pointer(ivari,xx,modul)
          postp(1) % array_used(ivari) = 1
          if( memory_size(xx) > 0 ) postp(1) % array_allocated(ivari) = 1 
          if(      postp(1) % time_posit(ivari) == 1 ) then
             postp(1) % time_num(ivari) = ndim1
          else if( postp(1) % time_posit(ivari) == 2 ) then
             postp(1) % time_num(ivari) = ndim2
          else if( postp(1) % time_posit(ivari) == 3 ) then
             postp(1) % time_num(ivari) = ndim3
          end if
          if(      postp(1) % comp_posit(ivari) == 1 ) then
             postp(1) % comp_num(ivari) = ndim1
          else if( postp(1) % comp_posit(ivari) == 2 ) then
             postp(1) % comp_num(ivari) = ndim2
          else if( postp(1) % comp_posit(ivari) == 3 ) then
             postp(1) % comp_num(ivari) = ndim3
          end if
          if( postp(1) % wopos(2,ivari) == 'VECTO' ) then
             if(      postp(1) % comp_posit(ivari) /= 1 .and. postp(1) % time_posit(ivari) /= 1 ) then
                postp(1) % dime_num(ivari) = ndim1
             else if( postp(1) % comp_posit(ivari) /= 2 .and. postp(1) % time_posit(ivari) /= 2 ) then
                postp(1) % dime_num(ivari) = ndim2
             else
                postp(1) % dime_num(ivari) = ndim3
             end if
          end if          
       else
          call runend('ARRAYS: DIMENSIONS ARE MISSING')
       end if

    case ( 'DEALLOCATE' )
       !
       ! deallocate
       !
       call memory_deallo(memor_loc,trim(variable_name_loc),trim(subroutine_name_loc),xx)
       postp(1) % array_used(ivari) = 0
       postp(1) % array_allocated(ivari) = 0

    case ( 'REDISTRIBUTE' )
       !
       ! Redistribute
       !
       call redistribution_array(xx,wopos_loc(3),POSIT=enti_posit_loc,MEMOR=memor_loc,VARIABLE_NAME=variable_name_loc)

    case ( 'INTERPOLATE' )
       ! 
       ! Interpolate
       !
       call AMR_interpolate(xx,wopos_loc(3),POSIT=enti_posit_loc,MEMOR=memor_loc,VARIABLE_NAME=variable_name_loc)

    case ( 'WRITE RESTART' )
       !
       ! Write restart
       !  
       if( postp(1) % array_used(ivari) == 1 ) then
          call messages_live('VARIABLE: '//wopos_loc(1)//' FROM '//namod(modul),'START SECTION')          
          call postpr_write_restart(xx,ivari,ittim,cutim,TAG1=TAG1)
          call messages_live('VARIABLE','END SECTION')
       end if

    case ( 'READ RESTART' )
       !
       ! Read restart
       !  
       if( postp(1) % array_used(ivari) == 1 ) then
          call messages_live('VARIABLE: '//wopos_loc(1)//' FROM '//namod(modul),'START SECTION')          
          call postpr_read_restart(xx,ivari,ittim,cutim,TAG1=TAG1)   
          call messages_live('VARIABLE','END SECTION')
       end if
       
    case ( 'POSTPROCESS' )
       !
       ! Postprocess
       !
       if( postp(1) % array_used(ivari) == 1 .or. optional_argument(.false.,FORCE) ) then
          if( present(MESH) ) then
             call postpr_postprocess(xx,ivari,ittim,cutim,MESH=MESH)             
          else if( present(MESH_ID) ) then
             if(     MESH_ID == 0 ) then
                call postpr_postprocess(xx,ivari,ittim,cutim)
             else if( MESH_ID > 0 ) then
                call postpr_postprocess(xx,ivari,ittim,cutim,MESH=witness_mesh(MESH_ID) % mesh)
             end if
          else
             call postpr_postprocess(xx,ivari,ittim,cutim)
          end if
      end if

    case default
       !
       ! Do not know what to do
       !
       call runend('ARRAYS: DO NOT KNOW WHAT TO DO '//trim(wtask)//' WITH '//trim(wopos_loc(1)))

    end select

    call arrays_options_end()

  end subroutine arrays_RP_3
  
  subroutine arrays_IP_1(ivari,wtask,xx,ndim1,ENTITY_POSITION,TIME_POSITION,TAG1,FORCE,MESH,MESH_ID)

    integer(ip),                                intent(in)    :: ivari
    character(len=*),                           intent(in)    :: wtask
    integer(ip),                      pointer,  intent(inout) :: xx(:)
    integer(ip),            optional,           intent(in)    :: ndim1
    integer(ip),            optional,           intent(in)    :: ENTITY_POSITION
    integer(ip),            optional,           intent(in)    :: TIME_POSITION
    integer(ip),            optional,           intent(in)    :: TAG1
    logical(lg),            optional,           intent(in)    :: FORCE
    type(mesh_type_basic),  optional,           intent(in)    :: MESH
    integer(ip),            optional,           intent(in)    :: MESH_ID

    call arrays_options_start(ivari)
    
    select case ( trim(wtask) )
       
    case ( 'ALLOCATE' )
       !
       ! Allocate
       !
       if( present(ndim1) ) then
          if( enti_posit_loc == 1 .and. ndim1 < nenti_loc ) call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' NOT REGISTERED CORRECTLY')
          if( enti_posit_loc /= 1 )                         call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' NOT REGISTERED CORRECTLY')
          if( trim(wopos_loc(4)) /= 'PRIMA' )               call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' SHOULD BE REGISTERED AS PRIMARY')
          call memory_alloca(memor_loc,trim(variable_name_loc),trim(subroutine_name_loc),xx,max(minmem,ndim1))
          if( associated(xx) ) call arrays_save_pointer(ivari,xx,modul)
          postp(1) % array_used(ivari) = 1
          if( memory_size(xx) > 0 ) postp(1) % array_allocated(ivari) = 1 
       else
          call runend('ARRAYS: DIMENSIONS ARE MISSING')
       end if
       
    case ( 'DEALLOCATE' )
       !
       ! deallocate
       !
       call memory_deallo(memor_loc,trim(variable_name_loc),trim(subroutine_name_loc),xx)
       postp(1) % array_used(ivari) = 0
       postp(1) % array_allocated(ivari) = 0
       
    case ( 'REDISTRIBUTE' )
       !
       ! Redistribute
       !
       call redistribution_array(xx,wopos_loc(3),POSIT=enti_posit_loc,MEMOR=memor_loc,VARIABLE_NAME=variable_name_loc)
       
    case ( 'INTERPOLATE' )
       ! 
       ! Interpolate
       !
       call AMR_interpolate(xx,wopos_loc(3),POSIT=enti_posit_loc,MEMOR=memor_loc,VARIABLE_NAME=variable_name_loc)

    case ( 'WRITE RESTART' )
       !
       ! Write restart
       !
       if( postp(1) % array_used(ivari) == 1 ) then
          call messages_live('VARIABLE: '//wopos_loc(1)//' FROM '//namod(modul),'START SECTION')
          call postpr_write_restart(xx,ivari,ittim,cutim,TAG1=TAG1)   
          call messages_live('VARIABLE','END SECTION')
       end if
         
    case ( 'READ RESTART' )
       !
       ! Write restart
       !  
       if( postp(1) % array_used(ivari) == 1 ) then
          call messages_live('VARIABLE: '//wopos_loc(1)//' FROM '//namod(modul),'START SECTION')
          call postpr_read_restart(xx,ivari,ittim,cutim,TAG1=TAG1)   
          call messages_live('VARIABLE','END SECTION')
       end if
       
    case ( 'POSTPROCESS' )
       !
       ! Postprocess
       !  
       if( postp(1) % array_used(ivari) == 1 .or. optional_argument(.false.,FORCE) ) then
          if( present(MESH) ) then
             call postpr_postprocess(xx,ivari,ittim,cutim,MESH=MESH)             
          else if( present(MESH_ID) ) then
             if(     MESH_ID == 0 ) then
                call postpr_postprocess(xx,ivari,ittim,cutim)
             else if( MESH_ID > 0 ) then
                call postpr_postprocess(xx,ivari,ittim,cutim,MESH=witness_mesh(MESH_ID) % mesh)
             end if
          else
             call postpr_postprocess(xx,ivari,ittim,cutim)
          end if
      end if

    case default
       !
       ! Do not know what to do
       !
       call runend('ARRAYS: DO NOT KNOW WHAT TO DO '//trim(wtask)//' WITH '//trim(wopos_loc(1)))
       
    end select

    call arrays_options_end()

  end subroutine arrays_IP_1
    
  subroutine arrays_R3P_1(ivari,wtask,xx,ndim1,ENTITY_POSITION,TIME_POSITION,FORCE,MESH,MESH_ID)

    integer(ip),                                intent(in)    :: ivari
    character(len=*),                           intent(in)    :: wtask
    type(r3p),                        pointer,  intent(inout) :: xx(:)
    integer(ip),            optional,           intent(in)    :: ndim1
    integer(ip),            optional,           intent(in)    :: ENTITY_POSITION
    integer(ip),            optional,           intent(in)    :: TIME_POSITION
    logical(lg),            optional,           intent(in)    :: FORCE
    type(mesh_type_basic),  optional,           intent(in)    :: MESH
    integer(ip),            optional,           intent(in)    :: MESH_ID

    call arrays_options_start(ivari)
    
    select case ( trim(wtask) )
       
    case ( 'ALLOCATE' )
       !
       ! Allocate
       !
       if( present(ndim1) ) then
          if( trim(wopos_loc(4)) /= 'PRIMA' )           call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' SHOULD BE REGISTERED AS PRIMARY')
          call memory_alloca(memor_loc,trim(variable_name_loc),trim(subroutine_name_loc),xx,max(minmem,ndim1))
          if( associated(xx) ) call arrays_save_pointer(ivari,xx,modul)
          postp(1) % array_used(ivari) = 1
          if( memory_size(xx) > 0 ) postp(1) % array_allocated(ivari) = 1 
       else
          call runend('ARRAYS: DIMENSIONS ARE MISSING')
       end if
       
    case ( 'DEALLOCATE' )
       !
       ! deallocate
       !
       call memory_deallo(memor_loc,trim(variable_name_loc),trim(subroutine_name_loc),xx)
       postp(1) % array_used(ivari) = 0
       postp(1) % array_allocated(ivari) = 0
       
    case ( 'REDISTRIBUTE' )
       !
       ! Redistribute
       !
       call redistribution_array(xx,wopos_loc(3),POSIT=enti_posit_loc,MEMOR=memor_loc,VARIABLE_NAME=variable_name_loc)
       
    case ( 'INTERPOLATE' )
       ! 
       ! Interpolate
       !
       call runend('NOT CODED')
       !call AMR_interpolate(xx,wopos_loc(3),POSIT=enti_posit_loc,MEMOR=memor_loc,VARIABLE_NAME=variable_name_loc)
       
    case ( 'WRITE RESTART' )
       !
       ! Write restart
       !
       if( postp(1) % array_used(ivari) == 1 ) then
          call messages_live('VARIABLE: '//wopos_loc(1)//' FROM '//namod(modul),'START SECTION')         
          call postpr_write_restart(xx,ivari,ittim,cutim)   
          call messages_live('VARIABLE','END SECTION')
       end if
       
    case ( 'READ RESTART' )
       !
       ! Write restart
       !  
       if( postp(1) % array_used(ivari) == 1 ) then
          call messages_live('VARIABLE: '//wopos_loc(1)//' FROM '//namod(modul),'START SECTION')
          call postpr_read_restart(xx,ivari,ittim,cutim)   
          call messages_live('VARIABLE','END SECTION')
       end if
       
    case ( 'POSTPROCESS' )
       !
       ! Postprocess
       !  
       if( postp(1) % array_used(ivari) == 1 .or. optional_argument(.false.,FORCE) ) then
          if( present(MESH) ) then
             call runend('arrays_R3P_1: POSTPROCESS NOT CODED 1')                
             !call postpr_postprocess(xx,ivari,ittim,cutim,MESH=MESH)             
          else if( present(MESH_ID) ) then
             if(     MESH_ID == 0 ) then
                call postpr_postprocess(xx,ivari,ittim,cutim)
             else if( MESH_ID > 0 ) then
                call runend('arrays_R3P_1: POSTPROCESS NOT CODED 2')                
                !call postpr_postprocess(xx,ivari,ittim,cutim,MESH=witness_mesh(MESH_ID) % mesh)
             end if
          else
             call postpr_postprocess(xx,ivari,ittim,cutim)
          end if
       end if
       
    case default
       !
       ! Do not know what to do
       !
       call runend('ARRAYS: DO NOT KNOW WHAT TO DO '//trim(wtask)//' WITH '//trim(wopos_loc(1)))
       
    end select

    call arrays_options_end()

  end subroutine arrays_R3P_1
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-13
  !> @brief   Starting operations
  !> @details Starting operations to deal with optional arguments
  !> 
  !-----------------------------------------------------------------------

  subroutine arrays_options_start(ivari)
    
    integer(ip),  intent(in) :: ivari

    if( modul == 0 ) then
       call runend('ARRAYS: YOU ARE TRYING TO REGISTER A VARIABLE OUT OF A MODULE')
    else
       kfl_reawr_loc       = kfl_reawr
       modul_loc           = modul
       variable_name_loc   = momod(modul_loc) % postp(1) % wopos(1,ivari)
       wopos_loc           = momod(modul_loc) % postp(1) % wopos(:,ivari)
       subroutine_name_loc = trim(exmod(modul_loc))//'_ARRAYS'
       memor_loc           = mem_modul(1:2,modul_loc)
       enti_posit_loc      = momod(modul_loc) % postp(1) % enti_posit (ivari)
       
       if(      momod(modul_loc) % postp(1) % wopos(3,ivari) == 'NPOIN' ) then
          nenti_loc = npoin
       else if( momod(modul_loc) % postp(1) % wopos(3,ivari) == 'NELEM' ) then
          nenti_loc = nelem
       else if( momod(modul_loc) % postp(1) % wopos(3,ivari) == 'NBOUN' ) then
          nenti_loc = nboun
       end if
    end if
    
  end subroutine arrays_options_start

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-13
  !> @brief   Ending operations
  !> @details Ending operations to deal with optional arguments
  !> 
  !-----------------------------------------------------------------------

  subroutine arrays_options_end()
        
    mem_modul(1:2,modul_loc) = memor_loc
    kfl_reawr = kfl_reawr_loc
    
  end subroutine arrays_options_end

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-15
  !> @brief   Check if variable is allocated
  !> @details Check if variable is allocated at least in one MPI
  !> 
  !-----------------------------------------------------------------------

  logical(lg) function arrays_allocated(wname,MODULE_NUMBER)
    character(len=*),           intent(in) :: wname
    integer(ip),      optional, intent(in) :: MODULE_NUMBER
    integer(ip)                            :: ivari

    ivari = arrays_number(wname,MODULE_NUMBER)

    if( ivari == 0 ) call runend('MOD_ARRAYS: YOU ARE CHECKING A NON REGISTERED VARIABLE '//trim(wname))
    arrays_allocated = .false.
    if( present(MODULE_NUMBER) ) then
       if( momod(MODULE_NUMBER) % postp(1) % array_allocated(ivari) /= 0 )  arrays_allocated = .true.
    else
       if( postp(1) % array_allocated(ivari) /= 0 ) arrays_allocated = .true.
    end if
    
  end function arrays_allocated

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-14
  !> @brief   Return array name
  !> @details ???
  !> 
  !-----------------------------------------------------------------------

  function arrays_tag(wname,MODULE_NUMBER) result(wopos)

    character(len=*),              intent(in) :: wname
    integer(ip),         optional, intent(in) :: MODULE_NUMBER
    character(5)                              :: wopos(5)
    integer(ip)                               :: ivari,imodu
    
    if( present(MODULE_NUMBER) ) then
       imodu = MODULE_NUMBER
    else
       imodu = modul
    end if

    wopos = ''
    do ivari = 1,nvarp
       if( momod(imodu) % postp(1) % wopos (1,ivari) == trim(wname) ) then
          wopos(1:5) = momod(imodu) % postp(1) % wopos (1:5,ivari)
          return
       end if
    end do
       
  end function arrays_tag
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-15
  !> @brief   Check if variable is allocated
  !> @details Check if variable is allocated at least in one MPI
  !> 
  !-----------------------------------------------------------------------

  logical(lg) function arrays_used_RP_1(xx)

    real(rp),   pointer, intent(in) :: xx(:)
    integer(ip)                     :: if_allocated

    if_allocated = memory_size(xx)
    call PAR_MAX(if_allocated)
    if( if_allocated > 0 ) then
       arrays_used_RP_1 = .true.
    else
       arrays_used_RP_1 = .false.
    end if
    
  end function arrays_used_RP_1
  
  logical(lg) function arrays_used_RP_2(xx)

    real(rp),   pointer, intent(in) :: xx(:,:)
    integer(ip)                     :: if_allocated

    if_allocated = memory_size(xx)
    call PAR_MAX(if_allocated)
    if( if_allocated > 0 ) then
       arrays_used_RP_2 = .true.
    else
       arrays_used_RP_2 = .false.
    end if
    
  end function arrays_used_RP_2
  
  logical(lg) function arrays_used_RP_3(xx)

    real(rp),   pointer, intent(in) :: xx(:,:,:)
    integer(ip)                     :: if_allocated

    if_allocated = memory_size(xx)
    call PAR_MAX(if_allocated)
    if( if_allocated > 0 ) then
       arrays_used_RP_3 = .true.
    else
       arrays_used_RP_3 = .false.
    end if
    
  end function arrays_used_RP_3
  
  logical(lg) function arrays_used_R3P_1(xx)

    type(r3p),  pointer, intent(in) :: xx(:)
    integer(ip)                     :: if_allocated

    if_allocated = memory_size(xx)
    call PAR_MAX(if_allocated)
    if( if_allocated > 0 ) then
       arrays_used_R3P_1 = .true.
    else
       arrays_used_R3P_1 = .false.
    end if
    
  end function arrays_used_R3P_1

end module mod_arrays
!> @}
