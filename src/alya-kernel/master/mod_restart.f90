!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup IO
!> @{
!> @file    mod_restart.f90
!> @author  houzeaux and Eduardo Perez
!> @date    2020-05-20
!> @brief   Restart
!> @details Some tools for restart
!-----------------------------------------------------------------------

module mod_restart

  use def_kintyp_basic,   only : ip,rp,lg
  use def_master,         only : momod
  use def_master,         only : modul,IPARALL
  use def_master,         only : lun_rstar,namod
  use def_master,         only : ISLAVE
  use def_master,         only : INOTSLAVE
  use def_master,         only : ITASK_READ_RESTART
  use def_master,         only : ITASK_WRITE_RESTART
  use mod_communications, only : PAR_BROADCAST
  use mod_messages,       only : messages_live
  use def_linked_list
#ifdef ALYA_FTI
  use mod_alya2fti,       only : alya2fti_RecoverVarInit
  use mod_alya2fti,       only : alya2fti_RecoverVarFinalize
  use mod_alya2fti,       only : alya2fti_InitICP
  use mod_alya2fti,       only : alya2fti_FinalizeICP
  use mod_alya2fti,       only : FTI_st
  use mod_alya2fti,       only : FTI_write_ckpt
  use mod_fti_config,     only : fti_config
  use FTI
#endif
  implicit none
  private

#ifdef __PGI
  interface restart_add
     module procedure &
          restart_add_s_ip,&
          restart_add_1_ip,&
          restart_add_2_ip,&
          restart_add_s_rp,&
          restart_add_1_rp,&
          restart_add_2_rp
  end interface restart_add
#else
  interface restart_add
     module procedure &
          restart_add_s,&
          restart_add_1,&
          restart_add_2
  end interface restart_add
#endif

#ifdef __PGI
  type ptr_inte_pack
     integer(ip), pointer :: p
  end type ptr_inte_pack

  type ptr_real_pack
     real(rp),    pointer :: p
  end type ptr_real_pack

  type ptr_logi_pack
     logical(lg), pointer :: p
  end type ptr_logi_pack

  type ptr_char_pack
     character(len=:), pointer :: p
  end type ptr_char_pack

  type(ptr_inte_pack), allocatable :: inte_pack(:)
  type(ptr_real_pack), allocatable :: real_pack(:)
  type(ptr_logi_pack), allocatable :: logi_pack(:)
  type(ptr_char_pack), allocatable :: char_pack(:)

  type(ptr_inte_pack), allocatable :: ptr_inte_tmp(:)
  type(ptr_real_pack), allocatable :: ptr_real_tmp(:)
  type(ptr_logi_pack), allocatable :: ptr_logi_tmp(:)
  type(ptr_char_pack), allocatable :: ptr_char_tmp(:)

  integer(ip)                      :: inte_size
  integer(ip)                      :: real_size
  integer(ip)                      :: logi_size
  integer(ip)                      :: char_size
#endif
  !
  ! FTI variable
  !
#ifdef ALYA_FTI
  integer(ip),         pointer     :: FTI_Tmp_ptr_scal_inte
  real(rp),            pointer     :: FTI_Tmp_ptr_scal_real
  logical(lg),         pointer     :: FTI_Tmp_ptr_scal_logi
  integer(ip),         pointer     :: FTI_Tmp_ptr_vect_inte(:)
  real(rp),            pointer     :: FTI_Tmp_ptr_vect_real(:)
  logical(lg),         pointer     :: FTI_Tmp_ptr_vect_logi(:)
  integer(ip),         pointer     :: FTI_Tmp_ptr_matr_inte(:,:)
  real(rp),            pointer     :: FTI_Tmp_ptr_matr_real(:,:)
  logical(lg),         pointer     :: FTI_Tmp_ptr_matr_logi(:,:)
  integer(ip)                      :: FTI_id
  integer(ip)                      :: FTI_error
#endif

  type(llist_typ),     pointer     :: llist
  integer(ip)                      :: itask_pgi
  integer(ip)                      :: nunit_loc
  logical(lg),         parameter   :: debug = .false.
  logical(lg),         parameter   :: ascii = .false.

  public :: restart_initialize
  public :: restart_finalize
  public :: restart_ini
  public :: restart_add
  public :: restart_end
  public :: itask_pgi

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-20
  !> @brief   Unit
  !> @details Unit to be used
  !> 
  !-----------------------------------------------------------------------

  integer(ip) function restart_unit(NUNIT)

    integer(ip), optional, intent(in)    :: NUNIT

    if( present(NUNIT) ) then
       restart_unit = NUNIT
    else
       if( modul == 0 ) then
          restart_unit = lun_rstar
       else
          restart_unit = momod(modul) % lun_rstar
       end if
    end if
    
  end function restart_unit

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-12-05
  !> @brief   Initialize
  !> @details Initialize restart by master 
  !> 
  !-----------------------------------------------------------------------

  subroutine restart_initialize(itask)

    integer(ip), optional, intent(in) :: itask
    if(       itask == ITASK_READ_RESTART ) then
#ifdef ALYA_FTI
       call alya2fti_RecoverVarInit()
#endif
    else if( itask == ITASK_WRITE_RESTART ) then
#ifdef ALYA_FTI
       call alya2fti_InitICP()               
#endif
    end if
    !
    ! Fix itask
    !
    itask_pgi = itask
    !
    ! Undefined unit
    !
    nunit_loc = 0
    
  end subroutine restart_initialize

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-12-05
  !> @brief   Initialize
  !> @details Initialize restart by master 
  !> 
  !-----------------------------------------------------------------------

  subroutine restart_finalize(itask)
    integer(ip), optional, intent(in) :: itask
    
#ifdef ALYA_FTI
    if( fti_config%enabled ) then
       if(       itask == ITASK_READ_RESTART ) then
          call alya2fti_RecoverVarFinalize()

       else if( itask == ITASK_WRITE_RESTART ) then
          call alya2fti_FinalizeICP()
       end if
    end if
#endif

    nunit_loc = 0
    
  end subroutine restart_finalize

#ifndef __PGI

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-12-05
  !> @brief   End restart 
  !> @details End restart 
  !> 
  !-----------------------------------------------------------------------

  subroutine restart_end(itask,NUNIT)

    use mod_mpio_config, only: mpio_config

    implicit none

    integer(ip),                 intent(in) :: itask
    integer(ip),       optional, intent(in) :: NUNIT
    class(node_typ),   pointer              :: node
    logical(lg)                             :: fti_enabled

    fti_enabled = .false.

#ifdef ALYA_FTI
    fti_enabled = fti_config%enabled
#endif


    if(      itask == ITASK_WRITE_RESTART ) then
       !
       ! Write restart
       !
       call llist % iterate(restart_write,nunit_loc)

    else if( itask == ITASK_READ_RESTART ) then
       !
       ! Read data
       !
       call restart_read()

       if(.not. fti_enabled) then
          !
          ! Check data
          !
          if( INOTSLAVE ) then
             node => llist % head
             do while( associated(node) )
                if( .not. node % exists() ) then
                   call messages_live('VARIABLE '//node % name//' NOT FOUND IN RESTART FILE, TAKING DEFAULT VALUE','WARNING')
                end if
                node => node % next
             end do
          end if
          !
          ! Broacast
          !
          if( IPARALL ) call llist % iterate(restart_bcast)
       end if
    end if
    !
    ! Deallocate
    !
    call llist % untarget()
    deallocate(llist)

  end subroutine restart_end

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-12-05
  !> @brief   Initialize list
  !> @details Initialize list
  !> 
  !-----------------------------------------------------------------------

  subroutine restart_ini(itask,NUNIT)

    integer(ip), optional, intent(in) :: itask
    integer(ip), optional, intent(in) :: NUNIT
    !
    ! Allocate linked list
    !
    allocate(llist)
    call llist % init()
    !
    ! Unit to be used
    !
    nunit_loc = restart_unit(NUNIT)
    
  end subroutine restart_ini

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-12-05
  !> @brief   Add variable
  !> @details Add a variable
  !> 
  !-----------------------------------------------------------------------

  subroutine restart_add_s(vals,vanam)
    class(*),               target     :: vals
    character(*), optional, intent(in) :: vanam
    class(node_typ),  pointer          :: node

    allocate(node)
    call node % init()

    node % name      =  vanam
    node % rank      =  0

    select type ( vals )

    type is ( integer(kind=ip) )

       node % vals      => vals
       node % type      =  TYPE_INTEGER       

    type is ( real(kind=rp) )

       node % vals      => vals
       node % type      =  TYPE_REAL       

    type is ( logical(kind=lg) )

       node % vals      => vals
       node % type      =  TYPE_LOGICAL       

    type is ( character(len=*) )

       node % vals      => vals
       node % type      =  TYPE_CHARACTER       

    end select

    call llist % add(node)

  end subroutine restart_add_s

  subroutine restart_add_1(val1,vanam)
    class(*),         target,   intent(in) :: val1(:)
    character(*),               intent(in) :: vanam
    class(node_typ),  pointer              :: node

    allocate(node)
    call node % init()

    node % name     =  vanam
    node % rank     =  1

    select type ( val1 )

    type is ( integer(kind=ip) )

       node % val1     => val1
       node % type     =  TYPE_INTEGER       

    type is ( real(kind=rp) )

       node % val1     => val1
       node % type     =  TYPE_REAL       

    type is ( logical(kind=lg) )

       node % val1     => val1
       node % type     =  TYPE_LOGICAL       

    type is ( character(len=*) )

       node % val1     => val1
       node % type     =  TYPE_CHARACTER       

       class default

       print*,'unknown'

    end select

    call llist % add(node)

  end subroutine restart_add_1

  subroutine restart_add_2(val2,vanam)
    class(*),         target,   intent(in) :: val2(:,:)
    character(*),               intent(in) :: vanam
    class(node_typ),  pointer              :: node

    allocate(node)
    call node % init()

    node % name     =  vanam
    node % rank     =  2

    select type ( val2 )

    type is ( integer(kind=ip) )

       node % val2     => val2
       node % type     =  TYPE_INTEGER       

    type is ( real(kind=rp) )

       node % val2     => val2
       node % type     =  TYPE_REAL       

    type is ( logical(kind=lg) )

       node % val2     => val2
       node % type     =  TYPE_LOGICAL       

    type is ( character(len=*) )

       node % val2     => val2
       node % type     =  TYPE_CHARACTER       

       class default

       print*,'unknown'

    end select

    call llist % add(node)

  end subroutine restart_add_2

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-12-05
  !> @brief   Write value
  !> @details Write value
  !> 
  !-----------------------------------------------------------------------

  subroutine restart_write(node,lunit)
    class(node_typ),           intent(inout) :: node
    integer(ip),     optional, intent(in)    :: lunit
    integer(ip)                              :: nn(2)
    integer(4)                               :: lunit4

#ifdef ALYA_FTI          
    if( fti_config%enabled ) then
       if( FTI_write_ckpt == 1 ) then
          call FTI_setIDFromString(node % name, FTI_id)
          if( associated(node % vals) ) then
             select type ( v => node % vals )
             type is ( integer  (kind=ip) ) 
                FTI_Tmp_ptr_scal_inte => v
                call FTI_Protect(FTI_id,FTI_Tmp_ptr_scal_inte,FTI_error) 
             type is ( real     (kind=rp) ) 
                FTI_Tmp_ptr_scal_real => v 
                call FTI_Protect(FTI_id,FTI_Tmp_ptr_scal_real,FTI_error)
             type is ( logical  (kind=lg) ) 
                FTI_Tmp_ptr_scal_logi => v 
                call FTI_Protect(FTI_id,FTI_Tmp_ptr_scal_logi,FTI_error)                   
             type is ( character(len=*)   ) 
                call runend('MOD_RESTART: NOT HANDLED BY FTI')
             end select
          end if
          if( associated(node % val1) ) then
             select type ( v => node % val1 )
             type is ( integer  (kind=ip) ) 
                FTI_Tmp_ptr_vect_inte => v 
                call FTI_Protect(FTI_id,FTI_Tmp_ptr_vect_inte,FTI_error) 
             type is ( real     (kind=rp) ) 
                FTI_Tmp_ptr_vect_real => v 
                call FTI_Protect(FTI_id,FTI_Tmp_ptr_vect_real,FTI_error) 
             type is ( logical  (kind=lg) )                    
                FTI_Tmp_ptr_vect_logi => v 
                call FTI_Protect(FTI_id,FTI_Tmp_ptr_vect_logi,FTI_error) 
             type is ( character(len=*)   )
                call runend('MOD_RESTART: NOT HANDLED BY FTI')
             end select
          end if
          if( associated(node % val2) ) then
             select type ( v => node % val2 )
             type is ( integer  (kind=ip) ) 
                FTI_Tmp_ptr_matr_inte => v 
                call FTI_Protect(FTI_id,FTI_Tmp_ptr_matr_inte,FTI_error)                    
             type is ( real     (kind=rp) ) 
                FTI_Tmp_ptr_matr_real => v 
                call FTI_Protect(FTI_id,FTI_Tmp_ptr_matr_real,FTI_error)                    
             type is ( logical  (kind=lg) ) 
                FTI_Tmp_ptr_matr_logi => v 
                call FTI_Protect(FTI_id,FTI_Tmp_ptr_matr_logi,FTI_error)                    
             type is ( character(len=*)   ) 
                call runend('MOD_RESTART: NOT HANDLED BY FTI') 
             end select
          end if
          call FTI_AddVarICP(FTI_id,FTI_error)
       end if
    else
#endif 
       if( present(lunit) ) then
          lunit4 = abs(int(lunit,4))
       else
          lunit4 = 10_4
       end if

       if( associated(node % vals) .or. associated(node % val1) .or. associated(node % val2) ) then       
          nn = node % getdim()
          if( nn(1) /= 0 .and. INOTSLAVE ) then

             if( debug ) print*,'writing ',node % name,' in unit ',lunit4,' for module ',namod(modul)

             if( ascii ) then
                write(lunit4,*) nn,node % type,node % rank,len(node % name,KIND=ip)
                write(lunit4,*) node % name

                if( associated(node % vals) ) then
                   select type ( v => node % vals )
                   type is ( integer  (kind=ip) ) ; write(lunit4,*) v 
                   type is ( real     (kind=rp) ) ; write(lunit4,*) v 
                   type is ( logical  (kind=lg) ) ; write(lunit4,*) v
                   type is ( character(len=*)   ) ; write(lunit4,*) v
                   end select
                end if
                if( associated(node % val1) ) then
                   select type ( v => node % val1 )
                   type is ( integer  (kind=ip) ) ; write(lunit4,*) v
                   type is ( real     (kind=rp) ) ; write(lunit4,*) v
                   type is ( logical  (kind=lg) ) ; write(lunit4,*) v
                   type is ( character(len=*)   ) ; write(lunit4,*) v
                   end select
                end if
                if( associated(node % val2) ) then
                   select type ( v => node % val2 )
                   type is ( integer  (kind=ip) ) ; write(lunit4,*) v
                   type is ( real     (kind=rp) ) ; write(lunit4,*) v
                   type is ( logical  (kind=lg) ) ; write(lunit4,*) v
                   type is ( character(len=*)   ) ; write(lunit4,*) v
                   end select
                end if
             else
                write(lunit4) nn,node % type,node % rank,len(node % name,KIND=ip)
                write(lunit4) node % name

                if( associated(node % vals) ) then
                   select type ( v => node % vals )
                   type is ( integer  (kind=ip) ) ; write(lunit4) v 
                   type is ( real     (kind=rp) ) ; write(lunit4) v 
                   type is ( logical  (kind=lg) ) ; write(lunit4) v
                   type is ( character(len=*)   ) ; write(lunit4) v
                   end select
                end if
                if( associated(node % val1) ) then
                   select type ( v => node % val1 )
                   type is ( integer  (kind=ip) ) ; write(lunit4) v
                   type is ( real     (kind=rp) ) ; write(lunit4) v
                   type is ( logical  (kind=lg) ) ; write(lunit4) v
                   type is ( character(len=*)   ) ; write(lunit4) v
                   end select
                end if
                if( associated(node % val2) ) then
                   select type ( v => node % val2 )
                   type is ( integer  (kind=ip) ) ; write(lunit4) v
                   type is ( real     (kind=rp) ) ; write(lunit4) v
                   type is ( logical  (kind=lg) ) ; write(lunit4) v
                   type is ( character(len=*)   ) ; write(lunit4) v
                   end select
                end if

             end if

          end if
       end if

#ifdef ALYA_FTI
    end if
#endif

  end subroutine restart_write

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-12-05
  !> @brief   Broadcast value
  !> @details Broadcast value
  !> 
  !-----------------------------------------------------------------------

  subroutine restart_bcast(node,ivari)
    class(node_typ),           intent(inout) :: node
    integer(ip),     optional, intent(in)    :: ivari
    integer(ip)                              :: nn(2)

    nn = node % getdim()
    if( associated(node % vals) ) then
       select type ( v => node % vals )
       type is ( integer  (kind=ip) ) ; call PAR_BROADCAST(v)
       type is ( real     (kind=rp) ) ; call PAR_BROADCAST(v)
       type is ( logical  (kind=lg) ) ; call PAR_BROADCAST(v)
       type is ( character(len=*)   ) ; call PAR_BROADCAST(nn(1),v,'IN MY CODE')
       end select
    else if( associated(node % val1) ) then
       select type ( v => node % val1 )
       type is ( integer  (kind=ip) ) ; call PAR_BROADCAST(nn(1),v)
       type is ( real     (kind=rp) ) ; call PAR_BROADCAST(nn(1),v)
       type is ( logical  (kind=lg) ) ; call PAR_BROADCAST(nn(1),v)
       end select
    else if( associated(node % val2) ) then
       select type ( v => node % val2 )
       type is ( integer  (kind=ip) ) ; call PAR_BROADCAST(nn(1),nn(2),v)
       type is ( real     (kind=rp) ) ; call PAR_BROADCAST(nn(1),nn(2),v)
       type is ( logical  (kind=lg) ) ; call PAR_BROADCAST(nn(1),nn(2),v)
       end select
    else
       return
    end if

  end subroutine restart_bcast

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-12-05
  !> @brief   Read value
  !> @details Read value
  !> 
  !-----------------------------------------------------------------------

  subroutine restart_read()

    integer(ip)                             :: nn(2),char_len
    class(node_typ),   pointer              :: node
    character(LEN=:),  allocatable          :: name
    integer(ip)                             :: typ,ierr
    integer(ip)                             :: rank
    integer(4)                              :: nunit4,istat4
    integer(ip),       pointer              :: vals_ip
    integer(ip),       pointer              :: val1_ip(:)
    integer(ip),       pointer              :: val2_ip(:,:)
    real(rp),          pointer              :: vals_rp
    real(rp),          pointer              :: val1_rp(:)
    real(rp),          pointer              :: val2_rp(:,:)
    logical(lg),       pointer              :: vals_lg
    logical(lg),       pointer              :: val1_lg(:)
    logical(lg),       pointer              :: val2_lg(:,:)
    character(LEN=:),  pointer              :: vals_ch

    if( llist % size <= 0 ) return

    nullify(vals_ip,val1_ip,val2_ip)
    nullify(vals_rp,val1_rp,val2_rp)
    nullify(vals_lg,val1_lg,val2_lg)
    nullify(vals_ch)

    nunit4 = int(nunit_loc,4)


#ifdef ALYA_FTI          
    if( fti_config%enabled ) then
       if( FTI_st /= 0 ) then
          node => llist % head
          do while( associated(node) )
             call FTI_getIDFromString(node % name, FTI_id)

             if( associated(node % vals) ) then
                select type ( v => node % vals )
                type is ( integer  (kind=ip) ) 
                   FTI_Tmp_ptr_scal_inte => v
                   call FTI_Protect(FTI_id,FTI_Tmp_ptr_scal_inte,FTI_error) 
                type is ( real     (kind=rp) ) 
                   FTI_Tmp_ptr_scal_real => v 
                   call FTI_Protect(FTI_id,FTI_Tmp_ptr_scal_real,FTI_error)
                type is ( logical  (kind=lg) ) 
                   FTI_Tmp_ptr_scal_logi => v 
                   call FTI_Protect(FTI_id,FTI_Tmp_ptr_scal_logi,FTI_error)                   
                type is ( character(len=*)   ) 
                   call runend('MOD_RESTART: NOT HANDLED BY FTI')
                end select
             end if
             if( associated(node % val1) ) then
                select type ( v => node % val1 )
                type is ( integer  (kind=ip) ) 
                   FTI_Tmp_ptr_vect_inte => v 
                   call FTI_Protect(FTI_id,FTI_Tmp_ptr_vect_inte,FTI_error) 
                type is ( real     (kind=rp) ) 
                   FTI_Tmp_ptr_vect_real => v 
                   call FTI_Protect(FTI_id,FTI_Tmp_ptr_vect_real,FTI_error) 
                type is ( logical  (kind=lg) )                    
                   FTI_Tmp_ptr_vect_logi => v 
                   call FTI_Protect(FTI_id,FTI_Tmp_ptr_vect_logi,FTI_error) 
                type is ( character(len=*)   )
                   call runend('MOD_RESTART: NOT HANDLED BY FTI')
                end select
             end if
             if( associated(node % val2) ) then
                select type ( v => node % val2 )
                type is ( integer  (kind=ip) ) 
                   FTI_Tmp_ptr_matr_inte => v 
                   call FTI_Protect(FTI_id,FTI_Tmp_ptr_matr_inte,FTI_error)                    
                type is ( real     (kind=rp) ) 
                   FTI_Tmp_ptr_matr_real => v 
                   call FTI_Protect(FTI_id,FTI_Tmp_ptr_matr_real,FTI_error)                    
                type is ( logical  (kind=lg) ) 
                   FTI_Tmp_ptr_matr_logi => v 
                   call FTI_Protect(FTI_id,FTI_Tmp_ptr_matr_logi,FTI_error)                    
                type is ( character(len=*)   ) 
                   call runend('MOD_RESTART: NOT HANDLED BY FTI') 
                end select
             end if

             call FTI_RecoverVar(FTI_id,FTI_error)
             node => node % next
          end do
       end if

    else
#endif

       if( INOTSLAVE ) then

          istat4 = 0
          if( nunit_loc == 0 ) call runend('MOD_RESTART: WRONG UNIT IN RESTART_END')

          do while( istat4 == 0 )

             if( ascii ) then
                read(nunit4,*,iostat=istat4) nn(1:2),typ,rank,char_len

                if( istat4 /= 0 ) goto 100

                allocate( character(char_len) :: name )
                read(nunit4,*) name

                if( debug ) print*,'reading ',name,' in unit ',nunit4,' for module ',namod(modul)

                select case ( typ )
                case( TYPE_INTEGER )
                   select case ( rank )
                   case ( 0_ip ) ; allocate( vals_ip )                     ; read(nunit4,*,iostat=istat4) vals_ip 
                   case ( 1_ip ) ; allocate( val1_ip(nn(1)) )              ; read(nunit4,*,iostat=istat4) val1_ip
                   case ( 2_ip ) ; allocate( val2_ip(nn(1),nn(2)) )        ; read(nunit4,*,iostat=istat4) val2_ip
                   end select
                case( TYPE_REAL )
                   select case ( rank )
                   case ( 0_ip ) ; allocate( vals_rp )                     ; read(nunit4,*,iostat=istat4) vals_rp 
                   case ( 1_ip ) ; allocate( val1_rp(nn(1)) )              ; read(nunit4,*,iostat=istat4) val1_rp
                   case ( 2_ip ) ; allocate( val2_rp(nn(1),nn(2)) )        ; read(nunit4,*,iostat=istat4) val2_rp
                   end select
                case( TYPE_LOGICAL ) 
                   select case ( rank )
                   case ( 0_ip ) ; allocate( vals_lg )                     ; read(nunit4,*,iostat=istat4) vals_lg
                   case ( 1_ip ) ; allocate( val1_lg(nn(1)) )              ; read(nunit4,*,iostat=istat4) val1_lg
                   case ( 2_ip ) ; allocate( val2_lg(nn(1),nn(2)) )        ; read(nunit4,*,iostat=istat4) val2_lg
                   end select
                case( TYPE_CHARACTER ) 
                   select case ( rank )
                   case ( 0_ip ) ; allocate( character(nn(1)) :: vals_ch ) ; read(nunit4,*,iostat=istat4) vals_ch
                   case default  ; call runend('MOD_RESTART: CHARACTERS MUST BE SCALAR')
                   end select
                end select

             else
                read(nunit4,iostat=istat4) nn(1:2),typ,rank,char_len

                if( istat4 /= 0 ) goto 100

                allocate( character(char_len) :: name )
                read(nunit4) name

                if( debug ) print*,'reading ',name,' in unit ',nunit4,' for module ',namod(modul)

                select case ( typ )
                case( TYPE_INTEGER )
                   select case ( rank )
                   case ( 0_ip ) ; allocate( vals_ip )                     ; read(nunit4,iostat=istat4) vals_ip 
                   case ( 1_ip ) ; allocate( val1_ip(nn(1)) )              ; read(nunit4,iostat=istat4) val1_ip
                   case ( 2_ip ) ; allocate( val2_ip(nn(1),nn(2)) )        ; read(nunit4,iostat=istat4) val2_ip
                   end select
                case( TYPE_REAL )
                   select case ( rank )
                   case ( 0_ip ) ; allocate( vals_rp )                     ; read(nunit4,iostat=istat4) vals_rp 
                   case ( 1_ip ) ; allocate( val1_rp(nn(1)) )              ; read(nunit4,iostat=istat4) val1_rp
                   case ( 2_ip ) ; allocate( val2_rp(nn(1),nn(2)) )        ; read(nunit4,iostat=istat4) val2_rp
                   end select
                case( TYPE_LOGICAL ) 
                   select case ( rank )
                   case ( 0_ip ) ; allocate( vals_lg )                     ; read(nunit4,iostat=istat4) vals_lg
                   case ( 1_ip ) ; allocate( val1_lg(nn(1)) )              ; read(nunit4,iostat=istat4) val1_lg
                   case ( 2_ip ) ; allocate( val2_lg(nn(1),nn(2)) )        ; read(nunit4,iostat=istat4) val2_lg
                   end select
                case( TYPE_CHARACTER ) 
                   select case ( rank )
                   case ( 0_ip ) ; allocate( character(nn(1)) :: vals_ch ) ; read(nunit4,iostat=istat4) vals_ch
                   case default  ; call runend('MOD_RESTART: CHARACTERS MUST BE SCALAR')
                   end select
                end select
             end if
             !
             ! Look for variable
             !
             node => llist % nameis(name)
             if( associated(node) ) then
                ierr = 0
                node % assigned = 1
                if( associated(vals_ip) ) call node % set(vals_ip)
                if( associated(val1_ip) ) then
                   call node % set(val1_ip,ierr)
                   if( ierr /= 0 ) call messages_live('VARIABLE '//name//' READ WITH DIFFERENT SIZE','WARNING')
                end if
                if( associated(val2_ip) ) then
                   call node % set(val2_ip,ierr)
                   if( ierr /= 0 ) call messages_live('VARIABLE '//name//' READ WITH DIFFERENT SIZE','WARNING')
                end if
                if( associated(vals_rp) ) call node % set(vals_rp)
                if( associated(val1_rp) ) then
                   call node % set(val1_rp,ierr) 
                   if( ierr /= 0 ) call messages_live('VARIABLE '//name//' READ WITH DIFFERENT SIZE','WARNING')
                end if
                if( associated(val2_rp) ) then
                   call node % set(val2_rp,ierr) 
                   if( ierr /= 0 ) call messages_live('VARIABLE '//name//' READ WITH DIFFERENT SIZE','WARNING')
                end if
                if( associated(vals_lg) ) call node % set(vals_lg)
                if( associated(val1_lg) ) call node % set(val1_lg)
                if( associated(val2_lg) ) call node % set(val2_lg)
                if( associated(vals_ch) ) call node % set(vals_ch)
             else
                call messages_live('VARIABLE '//name//' FOUND IN RESTART FILE BUT NOT USED','WARNING')
             end if
             if( allocated(name)     ) deallocate(name)
             if( associated(vals_ip) ) deallocate(vals_ip)
             if( associated(val1_ip) ) deallocate(val1_ip)
             if( associated(val2_ip) ) deallocate(val2_ip)
             if( associated(vals_rp) ) deallocate(vals_rp)
             if( associated(val1_rp) ) deallocate(val1_rp)
             if( associated(val2_rp) ) deallocate(val2_rp)
             if( associated(vals_lg) ) deallocate(vals_lg)
             if( associated(val1_lg) ) deallocate(val1_lg)
             if( associated(val2_lg) ) deallocate(val2_lg)
             if( associated(vals_ch) ) deallocate(vals_ch)
          end do

       end if

100    continue
#ifdef ALYA_FTI
    end if
#endif

  end subroutine restart_read

#else


  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-20
  !> @brief   Initialize restart
  !> @details Initialize restart
  !> 
  !-----------------------------------------------------------------------

  subroutine restart_ini(itask,NUNIT)

    integer(ip), optional, intent(in) :: itask
    integer(ip), optional, intent(in) :: NUNIT
    !
    ! Initialize variables
    !
    inte_size = 0
    real_size = 0
    logi_size = 0
    char_size = 0

    if( allocated(inte_pack) ) deallocate(inte_pack)
    if( allocated(real_pack) ) deallocate(real_pack)
    if( allocated(logi_pack) ) deallocate(logi_pack)
    if( allocated(char_pack) ) deallocate(char_pack)   
    !
    ! Unit to be used
    !
    nunit_loc = restart_unit(NUNIT)

  end subroutine restart_ini

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-20
  !> @brief   Finalize the task
  !> @details Finalize the task by writing or broadcasting
  !> 
  !-----------------------------------------------------------------------

  subroutine restart_end(itask)

    integer(ip),                intent(in) :: itask
    integer(4)                             :: nunit4
    integer(ip)                            :: ii,nn
    integer(ip),      allocatable          :: buffi(:)
    real(rp),         allocatable          :: buffr(:)
    logical(lg),      allocatable          :: buffl(:)
    character(len=:), allocatable          :: buffc

#ifndef ALYA_FTI

    nunit4 = int(nunit_loc,4)

    if(     itask == ITASK_WRITE_RESTART .and. INOTSLAVE ) then
       !
       ! Write restart
       ! 
       if( nunit_loc == 0 ) call runend('MOD_RESTART: WRONG UNIT')
       if( allocated(inte_pack) ) write(nunit4) ( inte_pack(ii) % p,ii=1,inte_size )   
       if( allocated(real_pack) ) write(nunit4) ( real_pack(ii) % p,ii=1,real_size )       
       if( allocated(logi_pack) ) write(nunit4) ( logi_pack(ii) % p,ii=1,logi_size )
       if( allocated(char_pack) ) write(nunit4) ( char_pack(ii) % p,ii=1,char_size )
       
    else if( itask == ITASK_READ_RESTART ) then
       !
       ! Read restart
       !       
       if( allocated(inte_pack) ) then
          allocate(buffi(inte_size))
          if( debug ) print*,'reading in unit ',nunit4,' for module ',namod(modul)
          if( INOTSLAVE ) read(nunit4) ( buffi(ii),ii=1,inte_size )
          call PAR_BROADCAST(inte_size,buffi)
          do ii = 1,inte_size
             inte_pack(ii) % p = buffi(ii) 
          end do
          deallocate(buffi)           
       end if

       if( allocated(real_pack) ) then
          allocate(buffr(real_size))
          if( INOTSLAVE ) read(nunit4) ( buffr(ii),ii=1,real_size )
          call PAR_BROADCAST(real_size,buffr)
          do ii = 1,real_size
             real_pack(ii) % p = buffr(ii) 
          end do
          deallocate(buffr)           
       end if

       if( allocated(logi_pack) ) then
          allocate(buffl(logi_size))
          if( INOTSLAVE ) read(nunit4) ( buffl(ii),ii=1,logi_size )
          call PAR_BROADCAST(logi_size,buffl)
          do ii = 1,logi_size
             logi_pack(ii) % p = buffl(ii) 
          end do
          deallocate(buffl)           
       end if

       if( allocated(char_pack) ) then
          if( INOTSLAVE ) read(nunit4) ( char_pack(ii) % p,ii=1,char_size )
          do ii = 1,char_size
             nn    = len(char_pack(ii) % p)
             buffc = char_pack(ii) % p
             call PAR_BROADCAST(nn,buffc,'IN MY CODE')
             char_pack(ii) % p = buffc
             deallocate(buffc)           
          end do
       end if

    end if
    !
    ! Reinitialize and deallocate
    !
    call restart_ini()
    !
    ! 
    !
    nunit_loc = 0_ip

#endif

  end subroutine restart_end

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-20
  !> @brief   Pack scalars 
  !> @details Pack scalars
  !> 
  !-----------------------------------------------------------------------

  subroutine restart_add_s_ip(val,vanam)
    integer(ip),  target, intent(in) :: val
    character(*),         intent(in) :: vanam
    integer(ip)                      :: itask

    itask = itask_pgi

#ifndef ALYA_FTI
    if( itask == ITASK_WRITE_RESTART .and. ISLAVE ) return
#endif


#ifdef ALYA_FTI
    FTI_Tmp_ptr_scal_inte => val
    if(      itask == ITASK_READ_RESTART .and. FTI_st /= 0) then
       call FTI_getIDFromString(vanam, FTI_id)
       call FTI_Protect(FTI_id,FTI_Tmp_ptr_scal_inte,FTI_error)
       call FTI_RecoverVar(FTI_id,FTI_error)
    else if( itask == ITASK_WRITE_RESTART .and. FTI_write_ckpt == 1 ) then
       call FTI_setIDFromString(vanam, FTI_id)
       call FTI_Protect(FTI_id,FTI_Tmp_ptr_scal_inte,FTI_error)
       call FTI_AddVarICP(FTI_id,FTI_error)
    end if
#else

    if( .not. allocated(inte_pack) )  then
       allocate(inte_pack(0))
       inte_size = 0
    end if
    inte_size = inte_size + 1

    if( inte_size > size(inte_pack) ) then
       allocate(ptr_inte_tmp(2_ip*inte_size))
       if( inte_size >= 1 ) ptr_inte_tmp(1:inte_size-1) = inte_pack(1:inte_size-1)
       call move_alloc(from=ptr_inte_tmp, to=inte_pack)
    end if

    inte_pack(inte_size) % p => val
#endif

  end subroutine restart_add_s_ip

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-20
  !> @brief   Pack scalars 
  !> @details Pack scalars
  !> 
  !-----------------------------------------------------------------------

  subroutine restart_add_s_rp(val,vanam)
    real(rp),     target, intent(in) :: val
    character(*),         intent(in) :: vanam
    integer(ip)                      :: itask

    itask = itask_pgi

#ifndef ALYA_FTI
    if( itask == ITASK_WRITE_RESTART .and. ISLAVE ) return
#endif


#ifdef ALYA_FTI
    FTI_Tmp_ptr_scal_real => val
    if(      itask == ITASK_READ_RESTART .and. FTI_st /= 0) then
       call FTI_getIDFromString(vanam, FTI_id)
       call FTI_Protect(FTI_id,FTI_Tmp_ptr_scal_real,FTI_error)
       call FTI_RecoverVar(FTI_id,FTI_error)
    else if( itask == ITASK_WRITE_RESTART .and. FTI_write_ckpt == 1 ) then
       call FTI_setIDFromString(vanam, FTI_id)
       call FTI_Protect(FTI_id,FTI_Tmp_ptr_scal_real,FTI_error)
       call FTI_AddVarICP(FTI_id,FTI_error)
    end if
#else
    if( .not. allocated(real_pack) )  then
       allocate(real_pack(0))
       real_size = 0
    end if
    real_size = real_size + 1

    if( real_size > size(real_pack) ) then 
       allocate(ptr_real_tmp(2_ip*real_size))
       if( real_size >= 1 ) ptr_real_tmp(1:real_size-1) = real_pack(1:real_size-1)
       call move_alloc(from=ptr_real_tmp, to=real_pack)
    end if

    real_pack(real_size) % p => val
#endif

  end subroutine restart_add_s_rp

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-20
  !> @brief   Pack arrays(:) 
  !> @details Pack arrays(:) 
  !> 
  !-----------------------------------------------------------------------

  subroutine restart_add_1_ip(val,vanam)
    character(*), optional, intent(in) :: vanam
    integer(ip),    target               :: val(:)
    integer(ip)                        :: kk,ii,itask

    itask = itask_pgi

#ifndef ALYA_FTI
    if( itask == ITASK_WRITE_RESTART .and. ISLAVE ) return
#endif

#ifdef ALYA_FTI
       FTI_Tmp_ptr_vect_inte => val
       if(      itask == ITASK_READ_RESTART .and. FTI_st /= 0) then
          call FTI_getIDFromString(vanam, FTI_id)
          call FTI_Protect(FTI_id,FTI_Tmp_ptr_vect_inte,FTI_error)
          call FTI_RecoverVar(FTI_id,FTI_error)
       else if( itask == ITASK_WRITE_RESTART .and. FTI_write_ckpt == 1 ) then
          call FTI_setIDFromString(vanam, FTI_id)
          call FTI_Protect(FTI_id,FTI_Tmp_ptr_vect_inte,FTI_error)
          call FTI_AddVarICP(FTI_id,FTI_error)
       end if
#else
       kk = size(val)
       if( .not. allocated(inte_pack) )  then
          allocate(inte_pack(0))
          inte_size = 0
       end if
       inte_size = inte_size + kk

       if( inte_size > size(inte_pack) ) then   
          allocate(ptr_inte_tmp(2_ip*inte_size))
          if( inte_size >= 1 ) ptr_inte_tmp(1:inte_size-kk) = inte_pack(1:inte_size-kk)
          call move_alloc(from=ptr_inte_tmp, to=inte_pack)
       end if

       do ii = 1,kk
          inte_pack(inte_size-kk+ii) % p => val(ii)
       end do
#endif

  end subroutine restart_add_1_ip

  subroutine restart_add_1_rp(val,vanam)
    character(*), optional, intent(in) :: vanam
    real(rp),     target               :: val(:)
    integer(ip)                        :: kk,ii,itask

    itask = itask_pgi

#ifndef ALYA_FTI
    if( itask == ITASK_WRITE_RESTART .and. ISLAVE ) return
#endif


#ifdef ALYA_FTI
       FTI_Tmp_ptr_vect_real => val
       if(      itask == ITASK_READ_RESTART .and. FTI_st /= 0) then
          call FTI_getIDFromString(vanam, FTI_id)
          call FTI_Protect(FTI_id,FTI_Tmp_ptr_vect_real,FTI_error)
          call FTI_RecoverVar(FTI_id,FTI_error)
       else if( itask == ITASK_WRITE_RESTART .and. FTI_write_ckpt == 1 ) then
          call FTI_setIDFromString(vanam, FTI_id)
          call FTI_Protect(FTI_id,FTI_Tmp_ptr_vect_real,FTI_error)
          call FTI_AddVarICP(FTI_id,FTI_error)
       end if
#else
       kk = size(val)
       if( .not. allocated(real_pack) )  then
          allocate(real_pack(0))
          real_size = 0
       end if
       real_size = real_size + kk

       if( real_size > size(real_pack) ) then 
          allocate(ptr_real_tmp(2_ip*real_size))
          if( real_size >= 1 ) ptr_real_tmp(1:real_size-kk) = real_pack(1:real_size-kk)
          call move_alloc(from=ptr_real_tmp, to=real_pack)
       end if

       do ii = 1,kk
          real_pack(real_size-kk+ii) % p => val(ii)
       end do
#endif

     end subroutine restart_add_1_rp

  subroutine restart_add_2_ip(val,vanam)
    integer(ip),      target, intent(in) :: val(:,:)
    character(LEN=*),         intent(in) :: vanam
    integer(ip)                          :: kk,ii,jj,ll,itask

    itask = itask_pgi

#ifndef ALYA_FTI
    if( itask == ITASK_WRITE_RESTART .and. ISLAVE ) return
#endif

#ifdef ALYA_FTI
       FTI_Tmp_ptr_matr_inte => val
       if(      itask == ITASK_READ_RESTART .and. FTI_st /= 0) then
          call FTI_getIDFromString(vanam, FTI_id)
          call FTI_Protect(FTI_id,FTI_Tmp_ptr_matr_inte,FTI_error)
          call FTI_RecoverVar(FTI_id,FTI_error)
       else if( itask == ITASK_WRITE_RESTART .and. FTI_write_ckpt == 1 ) then
          call FTI_setIDFromString(vanam, FTI_id)
          call FTI_Protect(FTI_id,FTI_Tmp_ptr_matr_inte,FTI_error)
          call FTI_AddVarICP(FTI_id,FTI_error)
       end if
#else
       kk = size(val)
       if( .not. allocated(inte_pack) )  then
          allocate(inte_pack(0))
          inte_size = 0
       end if
       inte_size = inte_size + kk

       if( inte_size > size(inte_pack) ) then   
          allocate(ptr_inte_tmp(2_ip*inte_size))
          if( inte_size >= 1 ) ptr_inte_tmp(1:inte_size-kk) = inte_pack(1:inte_size-kk)
          call move_alloc(from=ptr_inte_tmp, to=inte_pack)
       end if

       ll = 0
       do jj = 1,size(val,2)
          do ii = 1,size(val,1)
             ll = ll + 1
             inte_pack(inte_size-kk+ll) % p => val(ii,jj)
          end do
       end do
#endif

     end subroutine restart_add_2_ip

  subroutine restart_add_2_rp(val,vanam)
    real(rp),         target, intent(in) :: val(:,:)
    character(LEN=*),         intent(in) :: vanam
    integer(ip)                          :: kk,ii,jj,ll,itask

    itask = itask_pgi

#ifndef ALYA_FTI
    if( itask == ITASK_WRITE_RESTART .and. ISLAVE ) return
#endif

#ifdef ALYA_FTI
       FTI_Tmp_ptr_matr_real => val
       if(      itask == ITASK_READ_RESTART .and. FTI_st /= 0) then
          call FTI_getIDFromString(vanam, FTI_id)
          call FTI_Protect(FTI_id,FTI_Tmp_ptr_matr_real,FTI_error)
          call FTI_RecoverVar(FTI_id,FTI_error)
       else if( itask == ITASK_WRITE_RESTART .and. FTI_write_ckpt == 1 ) then
          call FTI_setIDFromString(vanam, FTI_id)
          call FTI_Protect(FTI_id,FTI_Tmp_ptr_matr_real,FTI_error)
          call FTI_AddVarICP(FTI_id,FTI_error)
       end if
#else
       kk = size(val)
       if( .not. allocated(real_pack) )  then
          allocate(real_pack(0))
          real_size = 0
       end if
       real_size = real_size + kk

       if( real_size > size(real_pack) ) then 
          allocate(ptr_real_tmp(2_ip*real_size))
          if( real_size >= 1 ) ptr_real_tmp(1:real_size-kk) = real_pack(1:real_size-kk)
          call move_alloc(from=ptr_real_tmp, to=real_pack)
       end if

       ll = 0
       do jj = 1,size(val,2)
          do ii = 1,size(val,1)
             ll = ll + 1
             real_pack(real_size-kk+ll) % p => val(ii,jj)
          end do
       end do
#endif

     end subroutine restart_add_2_rp

#endif

end module mod_restart
!> @}
