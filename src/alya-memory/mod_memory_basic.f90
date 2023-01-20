!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!>
!> @defgroup Memory_Toolbox
!> @{
!> @name    ToolBox for memory management
!> @file    mod_memory.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for memory management
!> @details ToolBox for memory management: allocate, deallocate
!>          for basic types defined in def_kintyp_basic
!>
!------------------------------------------------------------------------

module mod_memory_basic

  use def_kintyp_basic,      only : ip,rp,lg,i1p,i2p,i3p,r1p,r2p,r3p,r4p,i1pp
  use mod_optional_argument, only : optional_argument
#ifndef I_AM_NOT_ALYA
  use mod_std
#endif
  use mod_memory_tools
  
  implicit none

  private
  !
  ! MEMORY_ALLOC: allocate
  !
  interface memory_alloca
          !                 (:)                    (:,:)                 (:,:,:)               (:,:,:,:)
          !           --------------------------------------------------------------------------------------------
     module procedure &
          &           memory_alloca_r81    , memory_alloca_r82     , memory_alloca_r83    , memory_alloca_r84  ,  & ! REAL(RP)
          &           memory_alloca_r81_8  ,                                                                      & !
          &           memory_alloca_r41    , memory_alloca_r42     , memory_alloca_r43    , memory_alloca_r44  ,  & ! 
          &           memory_alloca_r41_8  ,                                                                      & !
          &           memory_alloca_i41    , memory_alloca_i42     , memory_alloca_i43    ,                       & ! INTEGER(4)
          &           memory_alloca_i418   ,                                                                      & ! 
          &           memory_alloca_i81    , memory_alloca_i81_mix , memory_alloca_i82    , memory_alloca_i83 ,   & ! INTEGER(8)
          &           memory_alloca_i1p_1  , memory_alloca_i1p_2   , memory_alloca_i1p_3  ,                       & ! TYPE(I1P)
          &           memory_alloca_i2p_1  , memory_alloca_i2p_2   ,                                              & ! TYPE(I2P)            
          &           memory_alloca_i3p_1  ,                                                                      & ! TYPE(I3P)                         
          &           memory_alloca_i1pp_1 ,                                                                      & ! TYPE(I1PP)
          &           memory_alloca_r1p_1  , memory_alloca_r1p_2   ,                                              & ! TYPE(R1P)
          &           memory_alloca_r2p_1  , memory_alloca_r2p_2   ,                                              & ! TYPE(R2P)
          &           memory_alloca_r3p_1  , memory_alloca_r3p_2   ,                                              & ! TYPE(R3P)
          &           memory_alloca_r4p_1  ,                                                                      & ! TYPE(R4P)
          &           memory_alloca_lg1    , memory_alloca_lg2     , memory_alloca_lg3    ,                       & ! LOGICAL(LG)
          &           memory_alloca_cha_4  , memory_alloca_cha_8   , memory_alloca_cha_0  ,                       & ! CHARACTER(:)
          &           memory_alloca_xp1    , memory_alloca_xp2                                                      ! COMPLEX(RP)
  end interface memory_alloca

  interface memory_deallo
          !                 (:)                    (:,:)                 (:,:,:)               (:,:,:,:)
          !           -------------------------------------------------------------------------------------------
     module procedure &
          &           memory_deallo_r81    , memory_deallo_r82     , memory_deallo_r83    , memory_deallo_r84   , & ! REAL(RP)
          &           memory_deallo_r41    , memory_deallo_r42     , memory_deallo_r43    , memory_deallo_r44   , & !
          &           memory_deallo_i41    , memory_deallo_i42     , memory_deallo_i43    ,                       & ! INTEGER(4)
          &           memory_deallo_i81    , memory_deallo_i82     , memory_deallo_i83    ,                       & ! INTEGER(8)
          &           memory_deallo_i1p_1  , memory_deallo_i1p_2   , memory_deallo_i1p_3  ,                       & ! TYPE(I1P)
          &           memory_deallo_i2p_1  , memory_deallo_i2p_2   ,                                              & ! TYPE(I2P)
          &           memory_deallo_i3p_1  ,                                                                      & ! TYPE(I3P)
          &           memory_deallo_i1pp_1 ,                                                                      & ! TYPE(I1PP)
          &           memory_deallo_r1p_1  , memory_deallo_r1p_2   ,                                              & ! TYPE(R1P)
          &           memory_deallo_r2p_1  , memory_deallo_r2p_2   ,                                              & ! TYPE(R2P)
          &           memory_deallo_r3p_1  , memory_deallo_r3p_2   ,                                              & ! TYPE(R3P)
          &           memory_deallo_r4p_1  ,                                                                      & ! TYPE(R4P)
          &           memory_deallo_lg1    , memory_deallo_lg2     , memory_deallo_lg3    ,                       & ! LOGICAL(LG)
          &           memory_deallo_cha    , memory_deallo_cha_0   ,                                              & ! CHARACTER(:)
          &           memory_deallo_xp1    , memory_deallo_xp2                                                      ! COMPLEX(RP)
  end interface memory_deallo
  
  interface memory_size
          !                 (:)                    (:,:)                 (:,:,:)               (:,:,:,:)
          !           --------------------------------------------------------------------------------------------
     module procedure &
          &           memory_size_rp1      , memory_size_rp2       , memory_size_rp3      , memory_size_rp4 ,     & ! REAL(R)
          &           memory_size_ip1_4    , memory_size_ip2_4     , memory_size_ip3_4    ,                       & ! INTEGER(4)
          &           memory_size_ip1_8    , memory_size_ip2_8     , memory_size_ip3_8    ,                       & ! INTEGER(8)
          &           memory_size_lg       ,                                                                      & ! LOGICAL(LG)
          &           memory_size_i1p_1    ,                                                                      & ! TYPE(I1P)
          &           memory_size_r1p_1    ,                                                                      & ! TYPE(R1P)
          &           memory_size_r2p_1    ,                                                                      & ! TYPE(R2P)
          &           memory_size_r3p_1                                                                             ! TYPE(R3P)
  end interface memory_size

  interface memory_initia
          !                 (:)                    (:,:)                 (:,:,:)               (:,:,:,:)
          !           --------------------------------------------------------------------------------------------
     module procedure &
          &           memory_initia_rp1    , memory_initia_rp2     , memory_initia_rp3                              ! REAL(RP)
  end interface memory_initia

  interface memory_alloca_min
          !                 (:)                    (:,:)                 (:,:,:)               (:,:,:,:)
          !           --------------------------------------------------------------------------------------------
     module procedure &
          &           memory_alloca_min_rp1, memory_alloca_min_rp2 , memory_alloca_min_rp3, memory_alloca_min_rp4,& ! REAL(RP)
          &           memory_alloca_min_ip1, memory_alloca_min_ip2 , memory_alloca_min_ip3,                       & ! INTEGER(IP)              
          &           memory_alloca_min_lg1, memory_alloca_min_lg2 , memory_alloca_min_lg3                          ! LOGICAL(LG)
  end interface memory_alloca_min

  interface memory_resize
          !                 (:)                    (:,:)                 (:,:,:)               (:,:,:,:)
          !           --------------------------------------------------------------------------------------------
     module procedure &
          &           memory_resize_rp1    , memory_resize_rp2     ,  memory_resize_rp3   ,                       & ! REAL(RP)
          &           memory_resize_41     , memory_resize_42      ,  memory_resize_43    ,                       & ! INTEGER(4)
          &           memory_resize_81     , memory_resize_82      ,  memory_resize_83                              ! INTEGER(8)
  end interface memory_resize

  interface memory_copy
          !                 (:)                    (:,:)                 (:,:,:)               (:,:,:,:)
          !           --------------------------------------------------------------------------------------------
     module procedure &
          &           memory_copy_rp1      , memory_copy_rp2       , memory_copy_rp3      ,                       & ! REAL(RP)
          &           memory_copy_i41      , memory_copy_i42       , memory_copy_i43      ,                       & ! INTEGER(4)
          &           memory_copy_i81      , memory_copy_i82       , memory_copy_i83      ,                       & ! INTEGER(8)
          &           memory_copy_i1p_1    ,                                                                      & ! TYPE(I1P)
          &           memory_deallo_r3p_1                                                                           ! TYPE(R3P)
  end interface memory_copy

  interface memory_renumber
     module procedure &
          &           memory_renumber_rp1  , memory_renumber_rp2   ,  memory_renumber_rp3  ,                      & ! REAL(RP)
          &           memory_renumber_i41  ,                                                                      & ! INTEGER(4)
          &           memory_renumber_i81                                                                           ! INTEGER(8)
  end interface memory_renumber

  interface memory_append
          !                 (:)                    (:,:)                 (:,:,:)               (:,:,:,:)
          !           --------------------------------------------------------------------------------------------
     module procedure &
          &           memory_append_ip1                                                                             ! INTEGER(IP)
  end interface memory_append

  public :: memory_initia                     ! Allocate memory
  public :: memory_alloca                     ! Allocate memory
  public :: memory_size                       ! Gives the size of a pointer (=0 if not associated)
  public :: memory_alloca_min                 ! Allocate a minimum memory for arrays
  public :: memory_deallo                     ! Deallocate memory
  public :: memory_copy                       ! Copy an array
  public :: memory_renumber                   ! Renumber an array
  public :: memory_resize                     ! Resize an array
  public :: memory_append                     ! Append an array to an array

  public :: memory_copy_r3p_1
  
contains

  subroutine memory_alloca_r81(memor,vanam,vacal,varia,ndim1,wzero,lboun,INIT_VALUE)
    !
    ! Real(rp)(:)
    !
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(4),   intent(in)              :: ndim1
    real(8),      intent(inout), pointer  :: varia(:)
    character(*), intent(in),    optional :: wzero
    integer(4),   intent(in),    optional :: lboun
    real(8),      intent(in),    optional :: INIT_VALUE
    integer(4)                            :: istat
    integer(4)                            :: idim1
    logical(lg)                           :: lzero

    if( ndim1 > 0 ) then

       if( kfl_alloc == 1 )    call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       if( present(lboun) ) then
          allocate( varia(lboun:lboun+ndim1-1) , stat = istat )
       else
          allocate( varia(ndim1) , stat = istat )
       end if
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( present(INIT_VALUE) ) then
             do idim1 = int(lbound(varia,1),4),int(ubound(varia,1),4)
                varia(idim1) = INIT_VALUE
             end do             
          else if( lzero ) then
             do idim1 = int(lbound(varia,1),4),int(ubound(varia,1),4)
                varia(idim1) = 0.0_8
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_r81

  subroutine memory_alloca_r81_8(memor,vanam,vacal,varia,ndim1,wzero,lboun,INIT_VALUE)
    !
    ! Real(rp)(:)
    !
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(8),   intent(in)              :: ndim1
    real(8),      intent(inout), pointer  :: varia(:)
    integer(8),   intent(in),    optional :: lboun
    character(*), intent(in),    optional :: wzero
    real(8),      intent(in),    optional :: INIT_VALUE
    integer(4)                            :: istat
    integer(8)                            :: idim1
    logical(lg)                           :: lzero

    if( ndim1 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       if( present(lboun) ) then
          allocate( varia(lboun:lboun+ndim1-1) , stat = istat )
       else
          allocate( varia(ndim1) , stat = istat )
       end if
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( present(INIT_VALUE) ) then
             do idim1 = int(lbound(varia,1),8),int(ubound(varia,1),8)
                varia(idim1) = INIT_VALUE
             end do             
          else if( lzero ) then
             do idim1 = int(lbound(varia,1),8),int(ubound(varia,1),8)
                varia(idim1) = 0.0_8
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_r81_8

  subroutine memory_alloca_r82(memor,vanam,vacal,varia,ndim1,ndim2,wzero,lboun1,lboun2)
    !
    ! Real(rp)(:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    integer(ip),  intent(in)            :: ndim2
    real(8),      intent(inout), pointer :: varia(:,:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun1
    integer(ip),  intent(in),  optional :: lboun2
    integer(ip)                         :: lboun1_loc
    integer(ip)                         :: lboun2_loc
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2
    logical(lg)                         :: lzero

    if( ndim1*ndim2 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)
       lboun1_loc = optional_argument(1_ip,lboun1)
       lboun2_loc = optional_argument(1_ip,lboun2)
       
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       allocate( varia(lboun1_loc:lboun1_loc+ndim1-1,lboun2_loc:lboun2_loc+ndim2-1) , stat = istat )

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim2 = int(lbound(varia,2),ip),int(ubound(varia,2),ip)
                do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
                   varia(idim1,idim2) = 0.0_8
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_r82

  subroutine memory_alloca_r83(memor,vanam,vacal,varia,ndim1,ndim2,ndim3,wzero,lboun1,lboun2,lboun3)
    !
    ! Real(rp)(:,:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    integer(ip),  intent(in)            :: ndim2
    integer(ip),  intent(in)            :: ndim3
    real(8),     intent(inout), pointer  :: varia(:,:,:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun1
    integer(ip),  intent(in),  optional :: lboun2
    integer(ip),  intent(in),  optional :: lboun3
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2
    integer(ip)                         :: idim3
    logical(lg)                         :: lzero

    if( ndim1*ndim2*ndim3 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun1) .and. present(lboun2) .and. present(lboun3) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1,lboun3:lboun3+ndim3-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2,ndim3) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim3 = int(lbound(varia,3),ip),int(ubound(varia,3),ip)
                do idim2 = int(lbound(varia,2),ip),int(ubound(varia,2),ip)
                   do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
                      varia(idim1,idim2,idim3) = 0_8
                   end do
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if
  end subroutine memory_alloca_r83

  subroutine memory_alloca_r84(memor,vanam,vacal,varia,ndim1,ndim2,ndim3,ndim4,wzero,lboun1,lboun2,lboun3,lboun4)
    !
    ! Real(rp)(:,:,:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    integer(ip),  intent(in)            :: ndim2
    integer(ip),  intent(in)            :: ndim3
    integer(ip),  intent(in)            :: ndim4
    real(8),     intent(inout), pointer :: varia(:,:,:,:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun1
    integer(ip),  intent(in),  optional :: lboun2
    integer(ip),  intent(in),  optional :: lboun3
    integer(ip),  intent(in),  optional :: lboun4
    integer(4)                          :: istat
    logical(lg)                         :: lzero

    if( ndim1*ndim2*ndim3*ndim4 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun1) .and. present(lboun2) .and. present(lboun3) .and. present(lboun4) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1,lboun3:lboun3+ndim3-1,lboun4:lboun4+ndim4-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2,ndim3,ndim4) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) varia = 0.0_8
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_r84

  subroutine memory_alloca_r41(memor,vanam,vacal,varia,ndim1,wzero,lboun,INIT_VALUE)
    !
    ! Real(rp)(:)
    !
    
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(4),   intent(in)              :: ndim1
    real(4),     intent(inout), pointer   :: varia(:)
    character(*), intent(in),   optional  :: wzero
    integer(4),   intent(in),    optional :: lboun
    real(4),      intent(in),   optional  :: INIT_VALUE
    integer(4)                            :: istat
    integer(4)                            :: idim1
    logical(lg)                           :: lzero

    if( ndim1 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       if( present(lboun) ) then
          allocate( varia(lboun:lboun+ndim1-1) , stat = istat )
       else
          allocate( varia(ndim1) , stat = istat )
       end if
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( present(INIT_VALUE) ) then
             do idim1 = int(lbound(varia,1),4),int(ubound(varia,1),4)
                varia(idim1) = INIT_VALUE
             end do             
          else if( lzero ) then
             do idim1 = int(lbound(varia,1),4),int(ubound(varia,1),4)
                varia(idim1) = 0.0_4
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_r41

  subroutine memory_alloca_r41_8(memor,vanam,vacal,varia,ndim1,wzero)
    !
    ! Real(rp)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   intent(in)            :: ndim1
    real(4),     intent(inout), pointer :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(4)                          :: istat
    integer(8)                          :: idim1
    logical(lg)                         :: lzero

    if( ndim1 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       allocate( varia(ndim1) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = 1,ndim1
                varia(idim1) = 0.0_4
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_r41_8

  subroutine memory_alloca_r42(memor,vanam,vacal,varia,ndim1,ndim2,wzero,lboun1,lboun2)
    !
    ! Real(rp)(:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    integer(ip),  intent(in)            :: ndim2
    real(4),     intent(inout), pointer  :: varia(:,:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun1
    integer(ip),  intent(in),  optional :: lboun2
    integer(ip)                         :: lboun1_loc
    integer(ip)                         :: lboun2_loc
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2
    logical(lg)                         :: lzero

    if( ndim1*ndim2 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)
       lboun1_loc = optional_argument(1_ip,lboun1)
       lboun2_loc = optional_argument(1_ip,lboun2)

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       allocate( varia(lboun1_loc:lboun1_loc+ndim1-1,lboun2_loc:lboun2_loc+ndim2-1) , stat = istat )

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim2 = int(lbound(varia,2),ip),int(ubound(varia,2),ip)
                do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
                   varia(idim1,idim2) = 0.0_4
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_r42

  subroutine memory_alloca_r43(memor,vanam,vacal,varia,ndim1,ndim2,ndim3,wzero,lboun1,lboun2,lboun3)
    !
    ! Real(rp)(:,:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    integer(ip),  intent(in)            :: ndim2
    integer(ip),  intent(in)            :: ndim3
    real(4),     intent(inout), pointer  :: varia(:,:,:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun1
    integer(ip),  intent(in),  optional :: lboun2
    integer(ip),  intent(in),  optional :: lboun3
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2
    integer(ip)                         :: idim3
    logical(lg)                         :: lzero

    if( ndim1*ndim2*ndim3 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun1) .and. present(lboun2) .and. present(lboun3) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1,lboun3:lboun3+ndim3-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2,ndim3) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim3 = int(lbound(varia,3),ip),int(ubound(varia,3),ip)
                do idim2 = int(lbound(varia,2),ip),int(ubound(varia,2),ip)
                   do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
                      varia(idim1,idim2,idim3) = 0_4
                   end do
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if
  end subroutine memory_alloca_r43

  subroutine memory_alloca_r44(memor,vanam,vacal,varia,ndim1,ndim2,ndim3,ndim4,wzero,lboun1,lboun2,lboun3,lboun4)
    !
    ! Real(rp)(:,:,:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    integer(ip),  intent(in)            :: ndim2
    integer(ip),  intent(in)            :: ndim3
    integer(ip),  intent(in)            :: ndim4
    real(4),     intent(inout), pointer :: varia(:,:,:,:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun1
    integer(ip),  intent(in),  optional :: lboun2
    integer(ip),  intent(in),  optional :: lboun3
    integer(ip),  intent(in),  optional :: lboun4
    integer(4)                          :: istat
    logical(lg)                         :: lzero

    if( ndim1*ndim2*ndim3*ndim4 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun1) .and. present(lboun2) .and. present(lboun3) .and. present(lboun4) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1,lboun3:lboun3+ndim3-1,lboun4:lboun4+ndim4-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2,ndim3,ndim4) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) varia = 0.0_4
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_r44

  subroutine memory_alloca_i41(memor,vanam,vacal,varia,ndim1,wzero,lboun,REALLOCATE,INIT_VALUE)
    !
    ! Integer(4)(:)
    !
    
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(4),   intent(in)              :: ndim1
    integer(4),   intent(inout), pointer  :: varia(:)
    character(*), intent(in),    optional :: wzero
    integer(4),   intent(in),    optional :: lboun
    logical(lg),  intent(in),    optional :: REALLOCATE
    integer(4),                  optional :: INIT_VALUE
    integer(4)                            :: istat
    integer(ip)                           :: idim1
    logical(lg)                           :: lzero
    logical(lg)                           :: linde
    integer(4)                            :: xvalu

    if( ndim1 > 0 ) then

       if( present(REALLOCATE) ) then
          if( REALLOCATE ) then
             call memory_deallo(memor,vanam,vacal,varia)
          end if
       end if

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       lzero = .true.
       linde = .false.
       xvalu = 0_4
       if( present(wzero) ) then
          if(      trim(wzero) == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          else if( trim(wzero) == 'HUGE') then
             xvalu = huge(0_4)
          else if( trim(wzero) == 'IDENTITY') then
             lzero = .false.
             linde = .true.
          end if
       end if

       if( present(lboun) ) then
          allocate( varia(lboun:lboun+ndim1-1) , stat = istat )
       else
          allocate( varia(ndim1) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( present(INIT_VALUE) ) then
             do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
                varia(idim1) = INIT_VALUE
             end do             
          else if( lzero ) then
             do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
                varia(idim1) = xvalu
             end do
          else if( linde ) then
             do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
                varia(idim1) = int(idim1,4)
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'integer(4)')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_i41



  subroutine memory_alloca_i418(memor,vanam,vacal,varia,ndim1,wzero,INIT_VALUE)
    !
    ! Integer(4)(:): wrt memory_alloca_i41, the dimension is integer(8)
    !
    
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(8),   intent(in)              :: ndim1
    integer(4),   intent(inout), pointer  :: varia(:)
    character(*), intent(in),    optional :: wzero
    integer(4),   intent(in),    optional :: INIT_VALUE
    integer(4)                            :: istat
    integer(ip)                           :: idim1
    logical(lg)                           :: lzero
    logical(lg)                           :: linde
    integer(4)                            :: xvalu

    if( ndim1 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       lzero = .true.
       linde = .false.
       xvalu = 0_4
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          else if( wzero == 'HUGE') then
             xvalu = huge(0_4)
          else if( wzero == 'IDENTITY') then
             lzero = .false.
             linde = .true.
          end if
       end if

       allocate( varia(ndim1) , stat = istat )

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( present(INIT_VALUE) ) then
             do idim1 = 1,int(ndim1,ip)
                varia(idim1) = INIT_VALUE
             end do
          else if( lzero ) then
             do idim1 = 1,int(ndim1,ip)
                varia(idim1) = xvalu
             end do
          else if( linde ) then
             do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
                varia(idim1) = int(idim1,4)
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'integer(4)')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_i418

  subroutine memory_alloca_i42(memor,vanam,vacal,varia,ndim1,ndim2,wzero,lboun1,lboun2,INIT_VALUE)
    !
    ! Integer(4)(:,:)
    !
    
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(4),   intent(in)              :: ndim1
    integer(4),   intent(in)              :: ndim2
    integer(4),   intent(inout), pointer  :: varia(:,:)
    character(*), intent(in),    optional :: wzero
    integer(4),   intent(in),    optional :: lboun1
    integer(4),   intent(in),    optional :: lboun2
    integer(4),                  optional :: INIT_VALUE
    integer(4)                            :: istat
    integer(4)                            :: idim1
    integer(4)                            :: idim2
    logical(lg)                           :: lzero

    if( ndim1*ndim2 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun1) .and. present(lboun2) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( present(INIT_VALUE) ) then
             do idim2 = int(lbound(varia,2),4),int(ubound(varia,2),4)
                do idim1 = int(lbound(varia,1),4),int(ubound(varia,1),4)
                   varia(idim1,idim2) = INIT_VALUE
                end do
             end do
          else if( lzero ) then
             do idim2 = int(lbound(varia,2),4),int(ubound(varia,2),4)
                do idim1 = int(lbound(varia,1),4),int(ubound(varia,1),4)
                   varia(idim1,idim2) = 0_ip
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'integer(4)')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_i42

  subroutine memory_alloca_i43(memor,vanam,vacal,varia,ndim1,ndim2,ndim3,wzero,lboun1,lboun2,lboun3)
    !
    ! Integer(4)(:,:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(4),   intent(in)            :: ndim1
    integer(4),   intent(in)            :: ndim2
    integer(4),   intent(in)            :: ndim3
    integer(4),   intent(inout), pointer  :: varia(:,:,:)
    character(*), intent(in),  optional :: wzero
    integer(4),   intent(in),  optional :: lboun1
    integer(4),   intent(in),  optional :: lboun2
    integer(4),   intent(in) , optional :: lboun3
    integer(4)                          :: istat
    integer(8)                          :: idim1
    integer(8)                          :: idim2
    integer(8)                          :: idim3
    logical(lg)                         :: lzero

    if( ndim1*ndim2*ndim3 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun1) .and. present(lboun2) .and. present(lboun3) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1,lboun3:lboun3+ndim3-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2,ndim3) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim3 = int(lbound(varia,3),8),int(ubound(varia,3),8)
                do idim2 = int(lbound(varia,2),8),int(ubound(varia,2),8)
                   do idim1 = int(lbound(varia,1),8),int(ubound(varia,1),8)
                      varia(idim1,idim2,idim3) = 0_4
                   end do
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'integer(4)')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_i43

  subroutine memory_alloca_lg1(memor,vanam,vacal,varia,ndim1,wzero,lboun)
    !
    ! Logical(lg)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    logical(lg),  intent(inout), pointer  :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun
    integer(4)                          :: istat
    integer(8)                          :: idim1
    logical(lg)                         :: lzero

    if( ndim1 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun) ) then
          allocate( varia(lboun:lboun+ndim1-1) , stat = istat )
       else
          allocate( varia(ndim1) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = int(lbound(varia,1),8),int(ubound(varia,1),8)
                varia(idim1) = .false.
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'logical')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_lg1

  subroutine memory_alloca_lg2(memor,vanam,vacal,varia,ndim1,ndim2,wzero,lboun1,lboun2)
    !
    ! Logical(lg)(:,:)
    !
    
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(ip),  intent(in)              :: ndim1
    integer(ip),  intent(in)              :: ndim2
    logical(lg),  intent(inout), pointer  :: varia(:,:)
    character(*), intent(in),    optional :: wzero
    integer(ip),  intent(in),    optional :: lboun1
    integer(ip),  intent(in),    optional :: lboun2
    integer(4)                            :: istat
    logical(lg)                           :: lzero

    if( ndim1*ndim2 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun1) .and. present(lboun2) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             varia = .false.
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'logical')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_lg2

  subroutine memory_alloca_lg3(memor,vanam,vacal,varia,ndim1,ndim2,ndim3,wzero,lboun1,lboun2,lboun3)
    !
    ! Logical(lg)(:,:,:)
    !
    
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(ip),  intent(in)              :: ndim1
    integer(ip),  intent(in)              :: ndim2
    integer(ip),  intent(in)              :: ndim3
    logical(lg),  intent(inout), pointer  :: varia(:,:,:)
    character(*), intent(in),    optional :: wzero
    integer(ip),  intent(in),    optional :: lboun1
    integer(ip),  intent(in),    optional :: lboun2
    integer(ip),  intent(in),    optional :: lboun3
    integer(4)                            :: istat
    logical(lg)                           :: lzero

    if( ndim1*ndim2*ndim3 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun1) .and. present(lboun2) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1,lboun3:lboun3+ndim3-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2,ndim3) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             varia = .false.
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'logical')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_lg3

  subroutine memory_alloca_r1p_1(memor,vanam,vacal,varia,ndim1,wzero,lboun)
    !
    ! Type(r1p)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    type(r1p),    intent(inout), pointer  :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun
    integer(ip)                         :: idim1
    integer(4)                          :: istat
    logical(lg)                         :: lzero

    if( ndim1 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       if( present(lboun) ) then
          allocate( varia(lboun:lboun+ndim1-1) , stat = istat )
       else
          allocate( varia(ndim1) , stat = istat )
       end if
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*size_r1p_type
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
                nullify(varia(idim1) % a)
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(r1p)')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_r1p_1

  subroutine memory_alloca_r1p_2(memor,vanam,vacal,varia,ndim1,ndim2,wzero)
    !
    ! Type(r1p)(:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    integer(ip),  intent(in)            :: ndim2
    type(r1p),    intent(inout), pointer  :: varia(:,:)
    character(*), intent(in),  optional :: wzero
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2
    integer(4)                          :: istat
    logical(lg)                         :: lzero

    if( ndim1*ndim2 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       allocate( varia(ndim1,ndim2) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*size_r1p_type
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim2 = 1,ndim2
                do idim1 = 1,ndim1
                   nullify(varia(idim1,idim2) % a)
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(r1p)')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_r1p_2

  subroutine memory_alloca_r2p_1(memor,vanam,vacal,varia,ndim1,wzero,lboun)
    !
    ! Type(r2p)(:)
    !
    
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(ip),  intent(in)              :: ndim1
    type(r2p),    intent(inout), pointer  :: varia(:)
    character(*), intent(in),    optional :: wzero
    integer(ip),  intent(in),    optional :: lboun
    integer(ip)                           :: idim1
    integer(4)                            :: istat
    logical(lg)                           :: lzero

    if( ndim1 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       if( present(lboun) ) then
          allocate( varia(lboun:lboun+ndim1-1) , stat = istat )
       else
          allocate( varia(ndim1) , stat = istat )
       end if

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*size_r2p_type
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
                nullify(varia(idim1) % a)
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(r2p)')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_r2p_1

  subroutine memory_alloca_r2p_2(memor,vanam,vacal,varia,ndim1,ndim2,wzero)
    !
    ! Type(r2p)(:,:)
    !
    
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(ip),  intent(in)              :: ndim1
    integer(ip),  intent(in)              :: ndim2
    type(r2p),    intent(inout), pointer  :: varia(:,:)
    character(*), intent(in),    optional :: wzero
    integer(ip)                           :: idim1
    integer(ip)                           :: idim2
    integer(4)                            :: istat
    logical(lg)                           :: lzero

    if( ndim1*ndim2 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       allocate( varia(ndim1,ndim2) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*size_r2p_type
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim2 = 1,ndim2
                do idim1 = 1,ndim1
                   nullify(varia(idim1,idim2) % a)
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(r2p)')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_r2p_2

  subroutine memory_alloca_r3p_1(memor,vanam,vacal,varia,ndim1,wzero,lboun)
    !
    ! Type(r3p)(:)
    !
    
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(ip),  intent(in)              :: ndim1
    type(r3p),    intent(inout), pointer  :: varia(:)
    character(*), intent(in),    optional :: wzero
    integer(ip),  intent(in),    optional :: lboun
    integer(ip)                           :: idim1
    integer(4)                            :: istat
    logical(lg)                           :: lzero

    if( ndim1 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       if( present(lboun) ) then
          allocate( varia(lboun:lboun+ndim1-1) , stat = istat )
       else
          allocate( varia(ndim1) , stat = istat )
       end if
       
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*size_r3p_type
          !lbytm = int(sizeof(varia),8) ! ndim1*ip
          if( lzero ) then
             do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
                nullify(varia(idim1) % a)
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(r3p)')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_r3p_1

  subroutine memory_alloca_r4p_1(memor,vanam,vacal,varia,ndim1,wzero)
    !
    ! Type(r4p)(:)
    !
    
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(ip),  intent(in)              :: ndim1
    type(r4p),    intent(inout), pointer  :: varia(:)
    character(*), intent(in),    optional :: wzero
    integer(ip)                           :: idim1
    integer(4)                            :: istat
    logical(lg)                           :: lzero

    if( ndim1 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       allocate( varia(ndim1) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*size_r4p_type
          !lbytm = int(sizeof(varia),8) ! ndim1*ip
          if( lzero ) then
             do idim1 = 1,ndim1
                nullify(varia(idim1) % a)
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(r4p)')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_r4p_1

  subroutine memory_alloca_r3p_2(memor,vanam,vacal,varia,ndim1,ndim2,wzero)
    !
    ! Type(r3p)(:,:)
    !
    
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(ip),  intent(in)              :: ndim1
    integer(ip),  intent(in)              :: ndim2
    type(r3p),    intent(inout), pointer  :: varia(:,:)
    character(*), intent(in),    optional :: wzero
    integer(ip)                           :: idim1,idim2
    integer(4)                            :: istat
    logical(lg)                           :: lzero

    if( ndim1*ndim2 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       allocate( varia(ndim1,ndim2) , stat = istat )
       if( istat == 0 ) then
          lbytm = int(ndim1,KIND=8)*int(ndim2,KIND=8)*int(ip,KIND=8)
          if( lzero ) then
             do idim2 = 1,ndim2
                do idim1 = 1,ndim1
                   nullify(varia(idim1,idim2) % a)
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(r3p)')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_r3p_2

  subroutine memory_alloca_i1p_2(memor,vanam,vacal,varia,ndim1,ndim2,wzero,lboun1,lboun2)
    !
    ! Type(i1p)(:,:)
    !
    
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(ip),  intent(in)              :: ndim1
    integer(ip),  intent(in)              :: ndim2
    type(i1p),    intent(inout), pointer  :: varia(:,:)
    character(*), intent(in),    optional :: wzero
    integer(ip),  intent(in),    optional :: lboun1
    integer(ip),  intent(in),    optional :: lboun2
    integer(ip)                           :: idim1
    integer(ip)                           :: idim2
    integer(4)                            :: istat
    logical(lg)                           :: lzero

    if( ndim1*ndim2 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun1) .and. present(lboun2) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)
          if( lzero ) then
             do idim2 = int(lbound(varia,2),ip),int(ubound(varia,2),ip)
                do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
                   nullify(varia(idim1,idim2) % l)
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(i1p)')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_i1p_2

  subroutine memory_alloca_i1p_3(memor,vanam,vacal,varia,ndim1,ndim2,ndim3,wzero,lboun1,lboun2,lboun3)
    !
    ! Type(i1p)(:,:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    integer(ip),  intent(in)            :: ndim2
    integer(ip),  intent(in)            :: ndim3
    type(i1p),    intent(inout), pointer  :: varia(:,:,:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun1
    integer(ip),  intent(in),  optional :: lboun2
    integer(ip),  intent(in),  optional :: lboun3
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2
    integer(ip)                         :: idim3
    integer(4)                          :: istat
    logical(lg)                         :: lzero

    if( ndim1*ndim2*ndim3 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       if( present(lboun1) .and. present(lboun2) .and. present(lboun3) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1,lboun3:lboun3+ndim3-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2,ndim3) , stat = istat )
       end if
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*size_i1p_type
          if( lzero ) then
             do idim3 = int(lbound(varia,3),ip),int(ubound(varia,3),ip)
                do idim2 = int(lbound(varia,2),ip),int(ubound(varia,2),ip)
                   do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
                      nullify(varia(idim1,idim2,idim3) % l)
                   end do
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(i1p)')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_i1p_3

  subroutine memory_alloca_i1p_1(memor,vanam,vacal,varia,ndim1,wzero,lboun)
    !
    ! Type(i1p)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    type(i1p),    intent(inout), pointer  :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun
    integer(ip)                         :: idim1
    integer(4)                          :: istat
    logical(lg)                         :: lzero

    if( ndim1 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       if( present(lboun) ) then
          allocate( varia(lboun:lboun+ndim1-1) , stat = istat )
       else
          allocate( varia(ndim1) , stat = istat )
       end if
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*size_i1p_type
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
                nullify(varia(idim1) % l)
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(i1p)')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_i1p_1

  subroutine memory_alloca_i1pp_1(memor,vanam,vacal,varia,ndim1,wzero)
    !
    ! Type(i1p)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    type(i1pp),   intent(inout), pointer  :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(ip)                         :: idim1
    integer(4)                          :: istat
    logical(lg)                         :: lzero

    if( ndim1 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       allocate( varia(ndim1) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*size_i1pp_type
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = 1,ndim1
                varia(idim1) % n = 0
                nullify(varia(idim1) % l)
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(i1p)')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_i1pp_1

  subroutine memory_alloca_i2p_1(memor,vanam,vacal,varia,ndim1,wzero,lboun)
    !
    ! Type(i2p)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    type(i2p),    intent(inout), pointer  :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun
    integer(ip)                         :: idim1
    integer(4)                          :: istat
    logical(lg)                         :: lzero

    if( ndim1 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun) ) then
          allocate( varia(lboun:lboun+ndim1-1) , stat = istat )
       else
          allocate( varia(ndim1) , stat = istat )
       end if

      if( istat == 0 ) then
          lbytm = size(varia,kind=8)*size_i2p_type
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
                nullify(varia(idim1) % l)
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(i2p)')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_i2p_1

  subroutine memory_alloca_i2p_2(memor,vanam,vacal,varia,ndim1,ndim2,wzero)
    !
    ! Type(i2p)(:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    integer(ip),  intent(in)            :: ndim2
    type(i2p),    intent(inout), pointer  :: varia(:,:)
    character(*), intent(in),  optional :: wzero
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2
    integer(4)                          :: istat
    logical(lg)                         :: lzero

    if( ndim1*ndim2 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       allocate( varia(ndim1,ndim2) , stat = istat )
       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*size_i2p_type
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = 1,ndim1
                do idim2 = 1,ndim2
                   nullify(varia(idim1,idim2) % l)
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(i2p)')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_i2p_2

  subroutine memory_alloca_i3p_1(memor,vanam,vacal,varia,ndim1,wzero)
    !
    ! Type(i3p)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    type(i3p),    intent(inout), pointer  :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(ip)                         :: idim1
    integer(4)                          :: istat
    logical(lg)                         :: lzero

    if( ndim1 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       allocate( varia(ndim1) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*size_i3p_type
          !lbytm = int(sizeof(varia),8) ! ndim1*ip
          if( lzero ) then
             do idim1 = 1,ndim1
                nullify(varia(idim1) % l)
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(i3p)')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_i3p_1

  subroutine memory_deallo_i41(memor,vanam,vacal,varia)
    !
    ! Integer(4)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(4),   pointer               :: varia(:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8)

       deallocate( varia , stat = istat )

       nullify(varia)

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)

       end if
       call memory_info(memor,vanam,vacal,'integer(4)')

    else

       lbytm = 0
       
    end if

  end subroutine memory_deallo_i41

  subroutine memory_deallo_i42(memor,vanam,vacal,varia)
    !
    ! Integer(4)(:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(4),   pointer               :: varia(:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'integer(4)')

    else

       lbytm = 0
      
    end if

  end subroutine memory_deallo_i42

  subroutine memory_deallo_i43(memor,vanam,vacal,varia)
    !
    ! Integer(4)(:,:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(4),   pointer               :: varia(:,:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8) !-size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'integer(4)')

    else

       lbytm = 0
      
    end if

  end subroutine memory_deallo_i43

  subroutine memory_deallo_r81(memor,vanam,vacal,varia)
    !
    ! Real(rp)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(8),     pointer                :: varia(:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8) ! -size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')

    else

       lbytm = 0
      
    end if

  end subroutine memory_deallo_r81

  subroutine memory_deallo_r82(memor,vanam,vacal,varia)
    !
    ! Real(rp)(:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(8),     pointer                :: varia(:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8) ! -size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')
       
    else

       lbytm = 0

    end if

  end subroutine memory_deallo_r82

  subroutine memory_deallo_r83(memor,vanam,vacal,varia)
    !
    ! Real(rp)(:,:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(8),     pointer                :: varia(:,:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')
    end if

  end subroutine memory_deallo_r83

  subroutine memory_deallo_r84(memor,vanam,vacal,varia)
    !
    ! Real(rp)(:,:,:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(8),     pointer                :: varia(:,:,:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8) ! -size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')
    end if

  end subroutine memory_deallo_r84

    subroutine memory_deallo_r41(memor,vanam,vacal,varia)
    !
    ! Real(rp)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(4),     pointer                :: varia(:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8) ! -size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')

    end if

  end subroutine memory_deallo_r41

  subroutine memory_deallo_r42(memor,vanam,vacal,varia)
    !
    ! Real(rp)(:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(4),     pointer                :: varia(:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8) ! -size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')

    end if

  end subroutine memory_deallo_r42

  subroutine memory_deallo_r43(memor,vanam,vacal,varia)
    !
    ! Real(rp)(:,:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(4),     pointer                :: varia(:,:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')
    end if

  end subroutine memory_deallo_r43

  subroutine memory_deallo_r44(memor,vanam,vacal,varia)
    !
    ! Real(rp)(:,:,:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(4),     pointer                :: varia(:,:,:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8) ! -size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')
    end if

  end subroutine memory_deallo_r44


  subroutine memory_deallo_lg1(memor,vanam,vacal,varia)
    !
    ! Logical(lg)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    logical(lg),  pointer               :: varia(:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'logical')

    end if

  end subroutine memory_deallo_lg1

  subroutine memory_deallo_lg2(memor,vanam,vacal,varia)
    !
    ! Logical(lg)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    logical(lg),  pointer               :: varia(:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')

    end if

  end subroutine memory_deallo_lg2

  subroutine memory_deallo_lg3(memor,vanam,vacal,varia)
    !
    ! Logical(lg)(:,:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    logical(lg),  pointer               :: varia(:,:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')

    end if

  end subroutine memory_deallo_lg3

  subroutine memory_deallo_r1p_1(memor,vanam,vacal,varia)
    !
    ! Type(r1p)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(r1p),    pointer               :: varia(:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1

    if( associated(varia) ) then

       do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
          call memory_deallo(memor,vanam//' % A',vacal,varia(idim1) % a)
       end do
       lbytm = -size(varia,kind=8)*size_r1p_type
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(r1p)')
    end if

  end subroutine memory_deallo_r1p_1

  subroutine memory_deallo_r1p_2(memor,vanam,vacal,varia)
    !
    ! Type(r1p)(:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(r1p),    pointer               :: varia(:,:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2

    if( associated(varia) ) then

       do idim2 = int(lbound(varia,2),ip),int(ubound(varia,2),ip)
          do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
             call memory_deallo(memor,vanam//' % A',vacal,varia(idim1,idim2) % a)
          end do 
       end do
       lbytm = -size(varia,kind=8)*size_r1p_type
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(r1p)')
    end if

  end subroutine memory_deallo_r1p_2

  subroutine memory_deallo_r2p_2(memor,vanam,vacal,varia)
    !
    ! Type(r2p)(:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(r2p),    pointer               :: varia(:,:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1,idim2

    if( associated(varia) ) then

       do idim2 = int(lbound(varia,2),ip),int(ubound(varia,2),ip)
          do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
             call memory_deallo(memor,vanam//' % A',vacal,varia(idim1,idim2) % a)
          end do
       end do
       lbytm = -size(varia,kind=8)*size_r2p_type
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(r2p)')
    end if

  end subroutine memory_deallo_r2p_2

  subroutine memory_deallo_r2p_1(memor,vanam,vacal,varia)
    !
    ! Type(r2p)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(r2p),    pointer               :: varia(:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1

    if( associated(varia) ) then

       do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
          call memory_deallo(memor,vanam//' % A',vacal,varia(idim1) % a)
       end do

       lbytm = -size(varia,kind=8)*size_r2p_type
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(r2p)')
    end if

  end subroutine memory_deallo_r2p_1

  subroutine memory_deallo_r3p_1(memor,vanam,vacal,varia)
    !
    ! Type(r3p)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(r3p),    pointer               :: varia(:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1

    if( associated(varia) ) then

       do idim1 = int(lbound(varia,1_ip),ip),int(ubound(varia,1_ip),ip)
          call memory_deallo(memor,vanam//' % A',vacal,varia(idim1) % a)
       end do
       lbytm = -size(varia,kind=8)*size_r3p_type
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(r3p)')
    end if

  end subroutine memory_deallo_r3p_1

  subroutine memory_deallo_r4p_1(memor,vanam,vacal,varia)
    !
    ! Type(r3p)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(r4p),    pointer               :: varia(:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1

    if( associated(varia) ) then

       do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
          call memory_deallo(memor,vanam//' % A',vacal,varia(idim1) % a)
       end do
       lbytm = -int(size(varia),8)*size_r4p_type
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(r4p)')
    end if

  end subroutine memory_deallo_r4p_1

  subroutine memory_deallo_r3p_2(memor,vanam,vacal,varia)
    !
    ! Type(r3p)(:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(r3p),    pointer               :: varia(:,:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2

    if( associated(varia) ) then

       do idim2 = int(lbound(varia,2_ip),ip),int(ubound(varia,2_ip),ip)
          do idim1 = int(lbound(varia,1_ip),ip),int(ubound(varia,1_ip),ip)
             call memory_deallo(memor,vanam//' % A',vacal,varia(idim1,idim2) % a)
          end do
       end do

       lbytm = -size(varia,kind=8)*size_r3p_type
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(r3p)')
    end if

  end subroutine memory_deallo_r3p_2

  subroutine memory_deallo_i1p_1(memor,vanam,vacal,varia)
    !
    ! Type(i1p)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(i1p),    pointer               :: varia(:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1

    if( associated(varia) ) then

       do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
          call memory_deallo(memor,vanam//' % L',vacal,varia(idim1) % l)
       end do
       
       lbytm = -size(varia,kind=8)*size_i1p_type
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(i1p)')

    end if

  end subroutine memory_deallo_i1p_1

  subroutine memory_deallo_i1p_2(memor,vanam,vacal,varia)
    !
    ! Type(i1p)(:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(i1p),    pointer               :: varia(:,:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2

    if( associated(varia) ) then

       do idim2 = int(lbound(varia,2),ip),int(ubound(varia,2),ip)
          do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
             call memory_deallo(memor,vanam//' % L',vacal,varia(idim1,idim2) % l)
          end do
       end do

       lbytm = -size(varia,kind=8)*size_i1p_type
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(i1p)')

    end if

  end subroutine memory_deallo_i1p_2

  subroutine memory_deallo_i1p_3(memor,vanam,vacal,varia)
    !
    ! Type(i1p)(:,:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(i1p),    pointer               :: varia(:,:,:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2
    integer(ip)                         :: idim3

    if( associated(varia) ) then

       do idim3 = int(lbound(varia,3),ip),int(ubound(varia,3),ip)
          do idim2 = int(lbound(varia,2),ip),int(ubound(varia,2),ip)
             do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
                call memory_deallo(memor,vanam//' % L',vacal,varia(idim1,idim2,idim3) % l)
             end do
          end do
       end do
       lbytm = -size(varia,kind=8)*size_i1p_type
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(i1p)')

    end if

  end subroutine memory_deallo_i1p_3

  subroutine memory_deallo_i1pp_1(memor,vanam,vacal,varia)
    !
    ! Type(i1pp)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(i1pp),   pointer               :: varia(:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1

    if( associated(varia) ) then

       do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
          call memory_deallo(memor,vanam//' % L',vacal,varia(idim1) % l)
       end do
       lbytm = -size(varia,kind=8)*size_i1pp_type
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(i1p)')

    end if

  end subroutine memory_deallo_i1pp_1

  subroutine memory_deallo_i2p_2(memor,vanam,vacal,varia)
    !
    ! Type(i2p)(:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(i2p),    pointer               :: varia(:,:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1,idim2

    if( associated(varia) ) then

       do idim2 = int(lbound(varia,2),ip),int(ubound(varia,2),ip)
          do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
             call memory_deallo(memor,vanam//' % L',vacal,varia(idim1,idim2) % l)
          end do
       end do

       lbytm = -size(varia,kind=8)*size_i2p_type
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(i2p)')

    end if

  end subroutine memory_deallo_i2p_2

  subroutine memory_deallo_i2p_1(memor,vanam,vacal,varia)
    !
    ! Type(i2p)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(i2p),    pointer               :: varia(:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1

    if( associated(varia) ) then

       do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
          call memory_deallo(memor,vanam//' % L',vacal,varia(idim1) % l)
       end do

       lbytm = -size(varia,kind=8)*size_i2p_type
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(i2p)')

    end if

  end subroutine memory_deallo_i2p_1

  subroutine memory_deallo_i3p_1(memor,vanam,vacal,varia)
    !
    ! Type(i3p)(:)
    !
    
    character(*), intent(in)             :: vanam         !< Variable name
    character(*), intent(in)             :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)          :: memor(2)      !< Memory counter
    type(i3p),    intent(inout), pointer :: varia(:)
    integer(4)                           :: istat
    integer(ip)                          :: idim1

    if( associated(varia) ) then

       do idim1 = int(lbound(varia,1_ip),ip),int(ubound(varia,1_ip),ip)
          call memory_deallo(memor,vanam//' % L',vacal,varia(idim1) % l)
       end do

       lbytm = -size(varia,kind=8)*size_i3p_type
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(i3p)')

    end if

  end subroutine memory_deallo_i3p_1

  subroutine memory_deallo_xp1(memor,vanam,vacal,varia)
    !
    ! complex(rp)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    complex(rp),     pointer            :: varia(:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)*2
       !lbytm = -int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')

    end if

  end subroutine memory_deallo_xp1

  subroutine memory_deallo_xp2(memor,vanam,vacal,varia)
    !
    ! complex(rp)(:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    complex(rp),  pointer               :: varia(:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)*2
       !lbytm = -int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')

    end if

  end subroutine memory_deallo_xp2

  subroutine memory_alloca_i81(memor,vanam,vacal,varia,ndim1,wzero,lboun,REALLOCATE,INIT_VALUE)
    !
    ! Integer(8)(:)
    !
    
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(8),   intent(in)              :: ndim1
    integer(8),   intent(inout), pointer  :: varia(:)
    character(*), intent(in),    optional :: wzero
    integer(ip),  intent(in),    optional :: lboun
    logical(lg),  intent(in),    optional :: REALLOCATE
    integer(8),   intent(in),    optional :: INIT_VALUE
    integer(4)                            :: istat
    integer(8)                            :: idim1
    logical(lg)                           :: lzero
    logical(lg)                           :: linde
    integer(8)                            :: xvalu

    if( ndim1 > 0 ) then

       if( present(REALLOCATE) ) then
          if( REALLOCATE ) then
             call memory_deallo(memor,vanam,vacal,varia)
          end if
       end if

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       if( present(lboun) ) then
          allocate( varia(lboun:lboun+int(ndim1,KIND=ip)-1_ip) , stat = istat )
       else
          allocate( varia(ndim1) , stat = istat )
       end if
       lzero = .true.
       linde = .false.
       xvalu = 0_8
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          else if( wzero == 'HUGE') then
             xvalu = huge(0_8)
          else if( wzero == 'IDENTITY') then
             lzero = .false.
             linde = .true.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8) !int(ndim1,8)*int(kind(varia),8)
          if( present(INIT_VALUE) ) then
             do idim1 = int(lbound(varia,1),8),int(ubound(varia,1),8)
                varia(idim1) = INIT_VALUE
             end do
          else if( lzero ) then
             do idim1 = int(lbound(varia,1),8),int(ubound(varia,1),8)
                varia(idim1) = xvalu
             end do
          else if( linde ) then
             do idim1 = int(lbound(varia,1),8),int(ubound(varia,1),8)
                varia(idim1) = idim1
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'integer(4)')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_i81

  subroutine memory_alloca_i81_mix(memor,vanam,vacal,varia,ndim1,wzero,lboun,REALLOCATE,INIT_VALUE)
    !
    ! Integer(8)(:)
    !
    
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(4),   intent(in)              :: ndim1
    integer(8),   intent(inout), pointer  :: varia(:)
    character(*), intent(in),    optional :: wzero
    integer(4),   intent(in),    optional :: lboun
    logical(lg),  intent(in),    optional :: REALLOCATE
    integer(8),                  optional :: INIT_VALUE
    integer(4)                            :: istat
    integer(8)                            :: idim1
    logical(lg)                           :: lzero
    logical(lg)                           :: linde
    integer(8)                            :: xvalu

    if( ndim1 > 0 ) then

       if( present(REALLOCATE) ) then
          if( REALLOCATE ) then
             call memory_deallo(memor,vanam,vacal,varia)
          end if
       end if

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       if( present(lboun) ) then
          allocate( varia(lboun:lboun+int(ndim1,KIND=4)-1_4) , stat = istat )
       else
          allocate( varia(ndim1) , stat = istat )
       end if
       lzero = .true.
       linde = .false.
       xvalu = 0_8
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          else if( wzero == 'HUGE') then
             xvalu = huge(0_8)
          else if( wzero == 'IDENTITY') then
             lzero = .false.
             linde = .true.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8) !int(ndim1,8)*int(kind(varia),8)
          if( present(INIT_VALUE) ) then
             do idim1 = int(lbound(varia,1),8),int(ubound(varia,1),8)
                varia(idim1) = INIT_VALUE
             end do             
          else if( lzero ) then
             do idim1 = int(lbound(varia,1),8),int(ubound(varia,1),8)
                varia(idim1) = xvalu
             end do
          else if( linde ) then
             do idim1 = int(lbound(varia,1),8),int(ubound(varia,1),8)
                varia(idim1) = idim1
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'integer(4)')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_i81_mix

  subroutine memory_alloca_i82(memor,vanam,vacal,varia,ndim1,ndim2,wzero,lboun1,lboun2,INIT_VALUE)
    !
    ! Integer(8)(:,:)
    !
    
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(8),   intent(in)              :: ndim1
    integer(8),   intent(in)              :: ndim2
    integer(8),   intent(inout), pointer  :: varia(:,:)
    character(*), intent(in),    optional :: wzero
    integer(ip),  intent(in),    optional :: lboun1
    integer(ip),  intent(in),    optional :: lboun2
    integer(8),                  optional :: INIT_VALUE
    integer(4)                            :: istat
    integer(ip)                           :: idim1
    integer(ip)                           :: idim2
    logical(lg)                           :: lzero

    if( ndim1*ndim2 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)

       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       if( present(lboun1) .and. present(lboun2) ) then
          allocate( varia(lboun1:lboun1+int(ndim1,KIND=ip)-1_ip,lboun2:lboun2+int(ndim2,KIND=ip)-1_ip) , stat = istat )
       else
          allocate( varia(ndim1,ndim2) , stat = istat )
       end if
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
         if( present(INIT_VALUE) ) then
             do idim2 = int(lbound(varia,2),ip),int(ubound(varia,2),ip)
                do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
                   varia(idim1,idim2) = INIT_VALUE
                end do
             end do
         else if( lzero ) then
             do idim2 = int(lbound(varia,2),ip),int(ubound(varia,2),ip)
                do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
                   varia(idim1,idim2) = 0_ip
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'integer(4)')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_i82

  subroutine memory_alloca_i83(memor,vanam,vacal,varia,ndim1,ndim2,ndim3,wzero,lboun1,lboun2,lboun3)
    !
    ! Integer(8)(:,:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   intent(in)            :: ndim1
    integer(8),   intent(in)            :: ndim2
    integer(8),   intent(in)            :: ndim3
    integer(8),   intent(inout), pointer  :: varia(:,:,:)
    character(*), intent(in),  optional :: wzero
    integer(8),   intent(in),  optional :: lboun1
    integer(8),   intent(in),  optional :: lboun2
    integer(8),   intent(in),  optional :: lboun3
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2
    integer(ip)                         :: idim3
    logical(lg)                         :: lzero

    if( ndim1*ndim2*ndim3 > 0 ) then

       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       if( present(lboun1) .and. present(lboun2) .and. present(lboun3) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1,lboun3:lboun3+ndim3-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2,ndim3) , stat = istat )
       end if
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim3 = int(lbound(varia,3),ip),int(ubound(varia,3),ip)
                do idim2 = int(lbound(varia,2),ip),int(ubound(varia,2),ip)
                   do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
                      varia(idim1,idim2,idim3) = 0_8
                   end do
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'integer(4)')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_i83

  subroutine memory_deallo_i81(memor,vanam,vacal,varia)
    !
    ! Integer(8)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   pointer               :: varia(:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'integer(4)')

    end if

  end subroutine memory_deallo_i81

  subroutine memory_deallo_i82(memor,vanam,vacal,varia)
    !
    ! Integer(8)(:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   pointer               :: varia(:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'integer(4)')

    end if

  end subroutine memory_deallo_i82

  subroutine memory_deallo_i83(memor,vanam,vacal,varia)
    !
    ! Integer(8)(:,:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   pointer               :: varia(:,:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'integer(4)')

    end if

  end subroutine memory_deallo_i83

  subroutine memory_alloca_cha_0(memor,vanam,vacal,varia,ndim1,wzero)
    !
    ! Character(:)
    !
    integer(8),       intent(inout)           :: memor(2)      !< Memory counter
    character(*),     intent(in)              :: vanam         !< Variable name
    character(*),     intent(in)              :: vacal         !< Calling subroutine name
    integer(ip),      intent(in)              :: ndim1
    character(len=:), intent(inout), pointer  :: varia
    character(*),     intent(in),    optional :: wzero
    integer(4)                                :: istat
    logical(lg)                               :: lzero

    if( ndim1 > 0 ) then

       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       allocate( character(ndim1) :: varia , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = int(ndim1,kind=8)*int(kind(varia),8)
          if( lzero ) then
             varia = ''
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'character')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_cha_0
  
  subroutine memory_alloca_cha_4(memor,vanam,vacal,varia,ndim1,wzero)
    !
    ! Character(*)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(4),   intent(in)            :: ndim1
    character(*), intent(inout), pointer  :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(4)                          :: istat
    integer(4)                          :: idim1
    logical(lg)                         :: lzero

    if( ndim1 > 0 ) then

       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       allocate( varia(ndim1) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = -int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = 1,ndim1
                varia(idim1) = ''
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'character')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_cha_4

  subroutine memory_alloca_cha_8(memor,vanam,vacal,varia,ndim1,wzero)
    !
    ! Character(*)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   intent(in)            :: ndim1
    character(*), intent(inout), pointer  :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(4)                          :: istat
    integer(8)                          :: idim1
    logical(lg)                         :: lzero

    if( ndim1 > 0 ) then

       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       allocate( varia(ndim1) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = -int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = 1,ndim1
                varia(idim1) = ''
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'character')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_cha_8

  subroutine memory_deallo_cha(memor,vanam,vacal,varia)
    !
    ! Character(*)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    character(*), pointer               :: varia(:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'character')

    end if

  end subroutine memory_deallo_cha

  subroutine memory_deallo_cha_0(memor,vanam,vacal,varia)
    !
    ! Character(*)
    !
    
    character(*),     intent(in)            :: vanam         !< Variable name
    character(*),     intent(in)            :: vacal         !< Calling subroutine name
    integer(8),       intent(inout)         :: memor(2)      !< Memory counter
    character(len=:), pointer               :: varia
    integer(4)                              :: istat
    integer(ip)                             :: ndim1
    
    if( associated(varia) ) then

       ndim1 = len(varia,kind=ip)
       lbytm = -int(ndim1,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'character')

    end if

  end subroutine memory_deallo_cha_0
  
  !----------------------------------------------------------------------
  !
  ! Copy arrays
  !
  ! 1. Allocate varia if it is null
  ! 2. varia <= vacpy (COPY_NAME) <= VANAM
  ! 3. Deallocate vacpy if 'DO_NOT_DEALLOCATE' not present
  !
  !----------------------------------------------------------------------

  subroutine memory_copy_i41(memor,vanam,vacal,vacpy,varia,wzero,COPY_NAME)
    !
    ! Integer(4)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(4),   intent(inout), pointer  :: varia(:)
    integer(4),   intent(inout), pointer  :: vacpy(:)
    character(*), intent(in),  optional :: wzero
    character(*), intent(in),  optional :: COPY_NAME     !< Variable name copy
    integer(ip)                         :: idim1
    logical(lg)                         :: lzero
    integer(ip)                         :: ndim1

    if( .not. associated(vacpy) ) then
       ndim1 = 0
    else
       ndim1 = size(vacpy,1_4,kind=ip)
    end if

    if( ndim1 > 0 ) then

       if( .not. associated(varia) ) then
          if( present(COPY_NAME) ) then
             call memory_alloca(memor,COPY_NAME,vacal,varia,ndim1,'DO_NOT_INITIALIZE')
          else
             call memory_alloca(memor,vanam,vacal,varia,ndim1,'DO_NOT_INITIALIZE')
          end if
       else
          ndim1 = min(ndim1,size(varia,kind=ip))
       end if

       do idim1 = 1,ndim1
          varia(idim1) = vacpy(idim1)
       end do

    else

       nullify(varia)

    end if

    lzero = .true.
    if(present(wzero)) then
       if( wzero == 'DO_NOT_DEALLOCATE') lzero = .false.
    end if
    if( lzero ) call memory_deallo(memor,vanam,vacal,vacpy)
    
  end subroutine memory_copy_i41

  subroutine memory_copy_i81(memor,vanam,vacal,vacpy,varia,wzero,COPY_NAME)
    !
    ! Integer(8)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   intent(inout), pointer  :: varia(:)
    integer(8),   intent(inout), pointer  :: vacpy(:)
    character(*), intent(in),  optional :: wzero
    character(*), intent(in),  optional :: COPY_NAME     !< Variable name copy
    integer(8)                          :: idim1
    logical(lg)                         :: lzero
    integer(8)                          :: ndim1

    if( .not. associated(vacpy) ) then
       ndim1 = 0
    else
       ndim1 = size(vacpy,kind=8)
    end if

    if( ndim1 > 0 ) then

       if( .not. associated(varia) ) then
          if( present(COPY_NAME) ) then
             call memory_alloca(memor,COPY_NAME,vacal,varia,ndim1,'DO_NOT_INITIALIZE')
          else
             call memory_alloca(memor,vanam,vacal,varia,ndim1,'DO_NOT_INITIALIZE')
          end if
       else
          ndim1 = min(ndim1,size(varia,kind=8))
       end if

       do idim1 = 1,ndim1
          varia(idim1) = vacpy(idim1)
       end do

    else

       nullify(varia)

    end if

    lzero = .true.
    if(present(wzero)) then
       if( wzero == 'DO_NOT_DEALLOCATE') lzero = .false.
    end if
    if( lzero ) call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_copy_i81

  subroutine memory_copy_i42(memor,vanam,vacal,vacpy,varia,wzero,COPY_NAME)
    !
    ! Integer(4)(:,:)
    !
    
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(4),   intent(inout), pointer  :: varia(:,:)
    integer(4),   intent(inout), pointer  :: vacpy(:,:)
    character(*), intent(in),    optional :: wzero
    character(*), intent(in),    optional :: COPY_NAME     !< Variable name copy
    integer(4)                            :: idim1,idim2
    logical(lg)                           :: lzero
    integer(4)                            :: ndim1,ndim2

    if( .not. associated(vacpy) ) then
       ndim1 = 0
       ndim2 = 0
    else
       ndim1 = size(vacpy,1_ip,kind=4)
       ndim2 = size(vacpy,2_ip,kind=4)
    end if

    if( ndim1 > 0 .and. ndim2 > 0 ) then

       if( .not. associated(varia) ) then
          if( present(COPY_NAME) ) then
             call memory_alloca(memor,COPY_NAME,vacal,varia,ndim1,ndim2,'DO_NOT_INITIALIZE')
          else
             call memory_alloca(memor,vanam,vacal,varia,ndim1,ndim2,'DO_NOT_INITIALIZE')
          end if
       else
          ndim1 = min(ndim1,size(varia,1,kind=4))
          ndim2 = min(ndim2,size(varia,2,kind=4))
       end if

       do idim2 = 1,ndim2
          do idim1 = 1,ndim1
             varia(idim1,idim2) = vacpy(idim1,idim2)
          end do
       end do

    else

       nullify(varia)

    end if

    lzero = .true.
    if(present(wzero)) then
       if( wzero == 'DO_NOT_DEALLOCATE') lzero = .false.
    end if
    if( lzero ) call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_copy_i42

  subroutine memory_copy_i43(memor,vanam,vacal,vacpy,varia,wzero,COPY_NAME)
    !
    ! Integer(4)(:,:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(4),   intent(inout), pointer  :: varia(:,:,:)
    integer(4),   intent(inout), pointer  :: vacpy(:,:,:)
    character(*), intent(in),  optional :: wzero
    character(*), intent(in),  optional :: COPY_NAME     !< Variable name copy
    integer(4)                          :: idim1,idim2,idim3
    logical(lg)                         :: lzero
    integer(4)                          :: ndim1,ndim2,ndim3

    if( .not. associated(vacpy) ) then
       ndim1 = 0
       ndim2 = 0
       ndim3 = 0
    else
       ndim1 = size(vacpy,1_ip,kind=4)
       ndim2 = size(vacpy,2_ip,kind=4)
       ndim3 = size(vacpy,3_ip,kind=4)
    end if

    if( ndim1 > 0 .and. ndim2 > 0 .and. ndim3 > 0 ) then

       if( .not. associated(varia) ) then
          if( present(COPY_NAME) ) then
             call memory_alloca(memor,COPY_NAME,vacal,varia,ndim1,ndim2,ndim3,'DO_NOT_INITIALIZE')
          else
             call memory_alloca(memor,vanam,vacal,varia,ndim1,ndim2,ndim3,'DO_NOT_INITIALIZE')
          end if
       else
          ndim1 = min(ndim1,size(varia,1,kind=4))
          ndim2 = min(ndim2,size(varia,2,kind=4))
          ndim3 = min(ndim3,size(varia,3,kind=4))
       end if

       do idim3 = 1,ndim3
          do idim2 = 1,ndim2
             do idim1 = 1,ndim1
                varia(idim1,idim2,idim3) = vacpy(idim1,idim2,idim3)
             end do
          end do
       end do

    else

       nullify(varia)

    end if

    lzero = .true.
    if(present(wzero)) then
       if( wzero == 'DO_NOT_DEALLOCATE') lzero = .false.
    end if
    if( lzero ) call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_copy_i43

  subroutine memory_copy_i82(memor,vanam,vacal,vacpy,varia,wzero,COPY_NAME)
    !
    ! Integer(8)(:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   intent(inout), pointer  :: varia(:,:)
    integer(8),   intent(inout), pointer  :: vacpy(:,:)
    character(*), intent(in),  optional :: wzero
     character(*), intent(in),  optional :: COPY_NAME     !< Variable name copy
   integer(8)                          :: idim1,idim2
    logical(lg)                         :: lzero
    integer(8)                          :: ndim1,ndim2

    if( .not. associated(vacpy) ) then
       ndim1 = 0
       ndim2 = 0
    else
       ndim1 = size(vacpy,1_ip,kind=8)
       ndim2 = size(vacpy,2_ip,kind=8)
    end if

    if( ndim1 > 0 .and. ndim2 > 0 ) then

       if( .not. associated(varia) ) then
          if( present(COPY_NAME) ) then
             call memory_alloca(memor,COPY_NAME,vacal,varia,ndim1,ndim2,'DO_NOT_INITIALIZE')
          else
             call memory_alloca(memor,vanam,vacal,varia,ndim1,ndim2,'DO_NOT_INITIALIZE')
          end if
       else
          ndim1 = min(ndim1,int(size(varia,1,kind=8),8))
          ndim2 = min(ndim2,int(size(varia,2,kind=8),8))
       end if

       do idim2 = 1,ndim2
          do idim1 = 1,ndim1
             varia(idim1,idim2) = vacpy(idim1,idim2)
          end do
       end do

    else

       nullify(varia)

    end if

    lzero = .true.
    if(present(wzero)) then
       if( wzero == 'DO_NOT_DEALLOCATE') lzero = .false.
    end if
    if( lzero ) call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_copy_i82

  subroutine memory_copy_i83(memor,vanam,vacal,vacpy,varia,wzero,COPY_NAME)
    !
    ! Integer(8)(:,:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   intent(inout), pointer  :: varia(:,:,:)
    integer(8),   intent(inout), pointer  :: vacpy(:,:,:)
    character(*), intent(in),  optional :: wzero
    character(*), intent(in),  optional :: COPY_NAME     !< Variable name copy
    integer(8)                          :: idim1,idim2,idim3
    logical(lg)                         :: lzero
    integer(8)                          :: ndim1,ndim2,ndim3

    if( .not. associated(vacpy) ) then
       ndim1 = 0
       ndim2 = 0
       ndim3 = 0
    else
       ndim1 = size(vacpy,1_ip,kind=8)
       ndim2 = size(vacpy,2_ip,kind=8)
       ndim3 = size(vacpy,3_ip,kind=8)
    end if

    if( ndim1 > 0 .and. ndim2 > 0 .and. ndim3 > 0 ) then

       if( .not. associated(varia) ) then
          if( present(COPY_NAME) ) then
             call memory_alloca(memor,COPY_NAME,vacal,varia,ndim1,ndim2,ndim3,'DO_NOT_INITIALIZE')
          else
             call memory_alloca(memor,vanam,vacal,varia,ndim1,ndim2,ndim3,'DO_NOT_INITIALIZE')
          end if
       else
          ndim1 = min(ndim1,int(size(varia,1,kind=8),8))
          ndim2 = min(ndim2,int(size(varia,2,kind=8),8))
          ndim3 = min(ndim3,int(size(varia,3,kind=8),8))
       end if

       do idim3 = 1,ndim3
          do idim2 = 1,ndim2
             do idim1 = 1,ndim1
                varia(idim1,idim2,idim3) = vacpy(idim1,idim2,idim3)
             end do
          end do
       end do

    else

       nullify(varia)

    end if

    lzero = .true.
    if(present(wzero)) then
       if( wzero == 'DO_NOT_DEALLOCATE') lzero = .false.
    end if
    if( lzero ) call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_copy_i83

  subroutine memory_copy_i1p_1(memor,vanam,vacal,vacpy,varia,wzero,COPY_NAME)
    !
    ! Type(i1p)(4)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(i1p),    intent(inout), pointer  :: varia(:)
    type(i1p),    intent(inout), pointer  :: vacpy(:)
    character(*), intent(in),  optional :: wzero
    character(*), intent(in),  optional :: COPY_NAME     !< Variable name copy
    integer(ip)                         :: idim1,idim2
    logical(lg)                         :: lzero
    integer(ip)                         :: ndim1
    integer(ip)                         :: ndim2

    if( .not. associated(vacpy) ) then
       ndim1 = 0
    else
       ndim1 = size(vacpy,1_ip,kind=ip)
    end if

    if( ndim1 > 0 ) then

       if( .not. associated(varia) ) then
          if( present(COPY_NAME) ) then
             call memory_alloca(memor,COPY_NAME,vacal,varia,ndim1,'DO_NOT_INITIALIZE')
          else
             call memory_alloca(memor,vanam,vacal,varia,ndim1,'DO_NOT_INITIALIZE')
          end if
          do idim1 = 1,ndim1
             ndim2 = size(vacpy(idim1) % l,1_ip,kind=ip)
             nullify(varia(idim1) % l)
             if( present(COPY_NAME) ) then
                call memory_alloca(memor,COPY_NAME//' % L',vacal,varia(idim1) % l,ndim2,'DO_NOT_INITIALIZE')
             else
                call memory_alloca(memor,vanam//' % L',vacal,varia(idim1) % l,ndim2,'DO_NOT_INITIALIZE')
             end if
          end do
       else
          ndim1 = min(ndim1,size(varia,kind=ip))
       end if

       do idim1 = 1,ndim1
          ndim2 = min(size(vacpy(idim1) % l,1_ip,kind=ip),size(varia(idim1) % l,1_ip,kind=ip))
          do idim2 = 1,ndim2
             varia(idim1) % l(idim2) = vacpy(idim1) % l(idim2)
          end do
       end do

    else

       nullify(varia)

    end if

    lzero = .true.
    if(present(wzero)) then
       if( wzero == 'DO_NOT_DEALLOCATE') lzero = .false.
    end if
    if( lzero ) call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_copy_i1p_1

  subroutine memory_copy_rp1(memor,vanam,vacal,vacpy,varia,wzero,COPY_NAME)
    !
    ! Real(rp)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(rp),     intent(inout), pointer  :: varia(:)
    real(rp),     intent(inout), pointer  :: vacpy(:)
    character(*), intent(in),  optional :: wzero
    character(*), intent(in),  optional :: COPY_NAME     !< Variable name copy
    integer(ip)                         :: idim1
    logical(lg)                         :: lzero
    integer(ip)                         :: ndim1

    ndim1 = memory_size(vacpy)

    if( ndim1 > 0 ) then

       if( .not. associated(varia) ) then
          if( present(COPY_NAME) ) then
             call memory_alloca(memor,COPY_NAME,vacal,varia,ndim1,'DO_NOT_INITIALIZE')
          else
             call memory_alloca(memor,vanam,vacal,varia,ndim1,'DO_NOT_INITIALIZE')
          end if
       else
          ndim1 = min(ndim1,memory_size(varia))
       end if

       do idim1 = 1,ndim1
          varia(idim1) = vacpy(idim1)
       end do

    else

       nullify(varia)

    end if

    lzero = .true.
    if(present(wzero)) then
       if( wzero == 'DO_NOT_DEALLOCATE') lzero = .false.
    end if
    if( lzero ) call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_copy_rp1

  subroutine memory_copy_rp2(memor,vanam,vacal,vacpy,varia,wzero,COPY_NAME)
    !
    ! Real(rp)(:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(rp),     intent(inout), pointer  :: varia(:,:)
    real(rp),     intent(inout), pointer  :: vacpy(:,:)
    character(*), intent(in),  optional :: wzero
    character(*), intent(in),  optional :: COPY_NAME     !< Variable name copy
    integer(ip)                         :: idim1,idim2
    logical(lg)                         :: lzero
    integer(ip)                         :: ndim1,ndim2

    if( .not. associated(vacpy) ) then
       ndim1 = 0
       ndim2 = 0
    else
       ndim1 = size(vacpy,1_ip,kind=ip)
       ndim2 = size(vacpy,2_ip,kind=ip)
    end if

    if( ndim1 > 0 .and. ndim2 > 0 ) then

       if( .not. associated(varia) ) then
          if( present(COPY_NAME) ) then
             call memory_alloca(memor,COPY_NAME,vacal,varia,ndim1,ndim2,'DO_NOT_INITIALIZE')
          else
             call memory_alloca(memor,vanam,vacal,varia,ndim1,ndim2,'DO_NOT_INITIALIZE')
          end if
       else
          ndim1 = min(ndim1,size(varia,1,kind=ip))
          ndim2 = min(ndim2,size(varia,2,kind=ip))
       end if

       do idim2 = 1,ndim2
          do idim1 = 1,ndim1
             varia(idim1,idim2) = vacpy(idim1,idim2)
          end do
       end do

    else

       nullify(varia)

    end if

    lzero = .true.
    if(present(wzero)) then
       if( wzero == 'DO_NOT_DEALLOCATE') lzero = .false.
    end if
    if( lzero ) call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_copy_rp2

 subroutine memory_copy_rp3(memor,vanam,vacal,vacpy,varia,wzero,COPY_NAME)
    !
    ! Real(rp)(:,:,:)
    !
    
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    real(rp),     intent(inout), pointer  :: varia(:,:,:)
    real(rp),     intent(inout), pointer  :: vacpy(:,:,:)
    character(*), intent(in),    optional :: wzero
    character(*), intent(in),    optional :: COPY_NAME     !< Variable name copy
    integer(ip)                           :: idim1,idim2,idim3
    logical(lg)                           :: lzero
    integer(ip)                           :: ndim1,ndim2,ndim3

    if( .not. associated(vacpy) ) then
       ndim1 = 0
       ndim2 = 0
       ndim3 = 0
    else
       ndim1 = size(vacpy,1_ip,kind=ip)
       ndim2 = size(vacpy,2_ip,kind=ip)
       ndim3 = size(vacpy,3_ip,kind=ip)
    end if

    if( ndim1 > 0 .and. ndim2 > 0 .and. ndim3 > 0 ) then

       if( .not. associated(varia) ) then
          if( present(COPY_NAME) ) then
             call memory_alloca(memor,COPY_NAME,vacal,varia,ndim1,ndim2,ndim3,'DO_NOT_INITIALIZE')
          else
             call memory_alloca(memor,vanam,vacal,varia,ndim1,ndim2,ndim3,'DO_NOT_INITIALIZE')
          end if
       else
          ndim1 = min(ndim1,size(varia,1,kind=ip))
          ndim2 = min(ndim2,size(varia,2,kind=ip))
          ndim3 = min(ndim3,size(varia,3,kind=ip))
       end if

       do idim3 = 1,ndim3
          do idim2 = 1,ndim2
             do idim1 = 1,ndim1
                varia(idim1,idim2,idim3) = vacpy(idim1,idim2,idim3)
             end do
          end do
       end do

    else

       nullify(varia)

    end if

    lzero = .true.
    if(present(wzero)) then
       if( wzero == 'DO_NOT_DEALLOCATE') lzero = .false.
    end if
    if( lzero ) call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_copy_rp3

  subroutine memory_copy_r3p_1(memor,vanam,vacal,vacpy,varia,wzero)
    !
    ! type(r3p)(:)
    !    
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    type(r3p),    intent(inout), pointer  :: varia(:)
    type(r3p),    intent(inout), pointer  :: vacpy(:)
    character(*), intent(in),    optional :: wzero
    integer(ip)                           :: idim1
    logical(lg)                           :: lzero
    integer(ip)                           :: ndim1
    integer(ip)                           :: n1,n2,n3,i1,i2,i3
     
    ndim1 = memory_size(vacpy)

    if( ndim1 > 0 ) then

       if( .not. associated(varia) ) then
          call memory_alloca(memor,vanam,vacal,varia,ndim1,'DO_NOT_INITIALIZE')
          do idim1 = 1,ndim1
             if( associated(vacpy(idim1) % a) ) then
                n1 = int(size(vacpy(idim1) % a,1),ip)
                n2 = int(size(vacpy(idim1) % a,2),ip)
                n3 = int(size(vacpy(idim1) % a,3),ip) 
                do i3 = 1,n3
                   do i2 = 1,n2
                      do i1 = 1,n1
                         call memory_alloca(memor,vanam//' % A',vacal,varia(idim1) % a,n1,n2,n3,'DO_NOT_INITIALIZE')
                      end do
                   end do
                end do
             end if
          end do                         
       else
          ndim1 = min(ndim1,memory_size(varia))
       end if

       do idim1 = 1,ndim1
          do i3 = 1,n3
             do i2 = 1,n2
                do i1 = 1,n1
                   varia(idim1) % a(i1,i2,i3) = vacpy(idim1) % a(i1,i2,i3)
                end do
             end do
          end do
       end do

    else

       nullify(varia)

    end if

    lzero = .true.
    if(present(wzero)) then
       if( wzero == 'DO_NOT_DEALLOCATE') lzero = .false.
    end if
    if( lzero ) call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_copy_r3p_1

  !----------------------------------------------------------------------
  !
  ! Resize arrays
  !
  ! 1. Nullify(vacpy)
  ! 2. vacpy <= varia, deallocate varia
  ! 3. Allocate varia with new dimension
  ! 4. varia <= vacpy, deallocate vacpy
  !
  !----------------------------------------------------------------------

  subroutine memory_resize_41(memor,vanam,vacal,varia,ndim1)
    !
    ! Integer(4)(:)
    !
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(4),   intent(inout), pointer  :: varia(:)
    integer(4),   intent(in)              :: ndim1
    integer(4),                  pointer  :: vacpy(:)
    integer(4)                            :: mdim1,ii
    integer(4)                            :: lboun1
    
    if( ndim1 == 0_4 ) then
       call memory_deallo(memor,vanam,vacal,varia)
    else
       nullify(vacpy)
       if( associated(varia) ) then
          if( size(varia) /= ndim1 ) then
             mdim1  = min(size(varia,KIND=4),ndim1)
             lboun1 = int(lbound(varia,1),4)
             call memory_alloca(memor,vanam,vacal,vacpy,ndim1,lboun=lboun1)
             do ii = lboun1,lboun1+mdim1-1
                vacpy(ii) = varia(ii)
             end do
             call memory_deallo(memor,vanam,vacal,varia)
             varia => vacpy
          else
             return
          end if
       else
          call memory_alloca(memor,vanam,vacal,varia,ndim1)
       end if
    end if

  end subroutine memory_resize_41
  
  subroutine memory_resize_42(memor,vanam,vacal,varia,ndim1,ndim2)
    !
    ! Integer(4)(:,:)
    !
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(4),   intent(inout), pointer  :: varia(:,:)
    integer(4),   intent(in)              :: ndim1
    integer(4),   intent(in)              :: ndim2
    integer(4),                  pointer  :: vacpy(:,:)
    integer(4)                            :: mdim1,mdim2,ii,jj
    integer(4)                            :: lboun1,lboun2
    
    if( ndim1 == 0_4 ) then
       call memory_deallo(memor,vanam,vacal,varia)
    else
       nullify(vacpy)
       if( associated(varia) ) then
          if( size(varia,1) /= ndim1 .or. size(varia,2) /= ndim2 ) then
             mdim1  = min(size(varia,KIND=4,DIM=1),ndim1)
             mdim2  = min(size(varia,KIND=4,DIM=2),ndim2)
             lboun1 = int(lbound(varia,1),4)
             lboun2 = int(lbound(varia,2),4)
             call memory_alloca(memor,vanam,vacal,vacpy,ndim1,ndim2,lboun1=lboun1,lboun2=lboun2)
             do jj = lboun2,lboun2+mdim2-1
                do ii = lboun1,lboun1+mdim1-1
                   vacpy(ii,jj) = varia(ii,jj)
                end do
             end do
             call memory_deallo(memor,vanam,vacal,varia)
             varia => vacpy
          else
             return
          end if
       else
          call memory_alloca(memor,vanam,vacal,varia,ndim1,ndim2)
       end if
    end if

  end subroutine memory_resize_42

  subroutine memory_resize_43(memor,vanam,vacal,varia,ndim1,ndim2,ndim3)
    !
    ! Integer(4)(:,:,:)
    !
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(4),   intent(inout), pointer  :: varia(:,:,:)
    integer(4),   intent(in)              :: ndim1
    integer(4),   intent(in)              :: ndim2
    integer(4),   intent(in)              :: ndim3
    integer(4),                  pointer  :: vacpy(:,:,:)
    integer(4)                            :: mdim1,mdim2,mdim3,ii,jj,kk
    integer(4)                            :: lboun1,lboun2,lboun3
    
    if( ndim1 == 0_4 ) then
       call memory_deallo(memor,vanam,vacal,varia)
    else
       nullify(vacpy)
       if( associated(varia) ) then
          if( size(varia,1) /= ndim1 .or. size(varia,2) /= ndim2 .or. size(varia,3) /= ndim3 ) then
             mdim1  = min(size(varia,KIND=4,DIM=1),ndim1)
             mdim2  = min(size(varia,KIND=4,DIM=2),ndim2)
             mdim3  = min(size(varia,KIND=4,DIM=3),ndim3)
             lboun1 = int(lbound(varia,1),4)
             lboun2 = int(lbound(varia,2),4)
             lboun3 = int(lbound(varia,3),4)
             call memory_alloca(memor,vanam,vacal,vacpy,ndim1,ndim2,ndim3,lboun1=lboun1,lboun2=lboun2,lboun3=lboun3)
             do kk = lboun3,lboun3+mdim3-1
                do jj = lboun2,lboun2+mdim2-1
                   do ii = lboun1,lboun1+mdim1-1
                      vacpy(ii,jj,kk) = varia(ii,jj,kk)
                   end do
                end do
             end do
             call memory_deallo(memor,vanam,vacal,varia)
             varia => vacpy
          else
             return
          end if
       else
          call memory_alloca(memor,vanam,vacal,varia,ndim1,ndim2,ndim3)
       end if
    end if

  end subroutine memory_resize_43

  subroutine memory_resize_81(memor,vanam,vacal,varia,ndim1)
    !
    ! Integer(8)(:)
    !
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(8),   intent(inout), pointer  :: varia(:)
    integer(8),   intent(in)              :: ndim1
    integer(8),                  pointer  :: vacpy(:)
    integer(8)                            :: mdim1,ii
    integer(8)                            :: lboun1
   
    if( ndim1 == 0_4 ) then
       call memory_deallo(memor,vanam,vacal,varia)
    else
       nullify(vacpy)
       if( associated(varia) ) then
          if( size(varia,1_ip,KIND=8) /= ndim1 ) then
             mdim1 = min(size(varia,KIND=8),ndim1)
             lboun1 = int(lbound(varia,1),8)
             call memory_alloca(memor,vanam,vacal,vacpy,ndim1,lboun=int(lboun1,ip))
             do ii = lboun1,lboun1+mdim1-1
                vacpy(ii) = varia(ii)
             end do
             call memory_deallo(memor,vanam,vacal,varia)
             varia => vacpy
          else
             return
          end if
       else
          call memory_alloca(memor,vanam,vacal,varia,ndim1)
       end if
    end if

  end subroutine memory_resize_81

  subroutine memory_resize_82(memor,vanam,vacal,varia,ndim1,ndim2)
    !
    ! Integer(8)(:,:)
    !
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(8),   intent(inout), pointer  :: varia(:,:)
    integer(8),   intent(in)              :: ndim1
    integer(8),   intent(in)              :: ndim2
    integer(8),                  pointer  :: vacpy(:,:)
    integer(8)                            :: mdim1,mdim2,ii,jj
    integer(8)                            :: lboun1,lboun2
    
    if( ndim1 == 0_4 ) then
       call memory_deallo(memor,vanam,vacal,varia)
    else
       nullify(vacpy)
       if( associated(varia) ) then
          if( size(varia,1_ip,KIND=8) /= ndim1 .or. size(varia,2_ip,KIND=8) /= ndim2 ) then
             mdim1  = min(int(size(varia,KIND=8,DIM=1),8),ndim1)
             mdim2  = min(int(size(varia,KIND=8,DIM=2),8),ndim2)
             lboun1 = int(lbound(varia,1),8)
             lboun2 = int(lbound(varia,2),8)
             call memory_alloca(memor,vanam,vacal,vacpy,ndim1,ndim2,lboun1=int(lboun1,ip),lboun2=int(lboun2,ip))
             do jj = lboun2,lboun2+mdim2-1
                do ii = lboun1,lboun1+mdim1-1
                   vacpy(ii,jj) = varia(ii,jj)
                end do
             end do
             call memory_deallo(memor,vanam,vacal,varia)
             varia => vacpy
          else
             return
          end if
       else
          call memory_alloca(memor,vanam,vacal,varia,ndim1,ndim2)
       end if
    end if

  end subroutine memory_resize_82

  subroutine memory_resize_83(memor,vanam,vacal,varia,ndim1,ndim2,ndim3)
    !
    ! Integer(8)(:,:,:)
    !
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(8),   intent(inout), pointer  :: varia(:,:,:)
    integer(8),   intent(in)              :: ndim1
    integer(8),   intent(in)              :: ndim2
    integer(8),   intent(in)              :: ndim3
    integer(8),                  pointer  :: vacpy(:,:,:)
    integer(8)                            :: mdim1,mdim2,mdim3,ii,jj,kk
    integer(8)                            :: lboun1,lboun2,lboun3
    
    if( ndim1 == 0_4 ) then
       call memory_deallo(memor,vanam,vacal,varia)
    else
       nullify(vacpy)
       if( associated(varia) ) then
          if( size(varia,1_ip,KIND=8) /= ndim1 .or. size(varia,2_ip,KIND=8) /= ndim2 .or. size(varia,3_ip,KIND=8) /= ndim3 ) then
             mdim1  = min(int(size(varia,KIND=8,DIM=1),8),ndim1)
             mdim2  = min(int(size(varia,KIND=8,DIM=2),8),ndim2)
             mdim3  = min(int(size(varia,KIND=8,DIM=3),8),ndim3)
             lboun1 = int(lbound(varia,1),8)
             lboun2 = int(lbound(varia,2),8)
             lboun3 = int(lbound(varia,3),8)
             call memory_alloca(memor,vanam,vacal,vacpy,ndim1,ndim2,ndim3,lboun1=lboun1,lboun2=lboun2,lboun3=lboun3)
             do kk = lboun3,lboun3+mdim3-1
                do jj = lboun2,lboun2+mdim2-1
                   do ii = lboun1,lboun1+mdim1-1
                      vacpy(ii,jj,kk) = varia(ii,jj,kk)
                   end do
                end do
             end do
             call memory_deallo(memor,vanam,vacal,varia)
             varia => vacpy
          else
             return
          end if
       else
          call memory_alloca(memor,vanam,vacal,varia,ndim1,ndim2,ndim3)
       end if
    end if

  end subroutine memory_resize_83

  subroutine memory_resize_rp1(memor,vanam,vacal,varia,ndim1)
    !
    ! Real(rp)(:)
    !
    character(*),  intent(in)              :: vanam         !< Variable name
    character(*),  intent(in)              :: vacal         !< Calling subroutine name
    integer(8),    intent(inout)           :: memor(2)      !< Memory counter
    real(rp),      intent(inout), pointer  :: varia(:)
    integer(ip),   intent(in)              :: ndim1
    real(rp),                     pointer  :: vacpy(:)
    integer(ip)                            :: mdim1,ii
    integer(ip)                            :: lboun1
    
    if( ndim1 == 0_4 ) then
       call memory_deallo(memor,vanam,vacal,varia)
    else
       nullify(vacpy)
       if( associated(varia) ) then
          if( size(varia,1_ip,KIND=ip) /= ndim1 ) then
             mdim1  = min(size(varia,KIND=ip),ndim1)
             lboun1 = int(lbound(varia,1),ip)
             call memory_alloca(memor,vanam,vacal,vacpy,ndim1,lboun=lboun1)
             do ii = lboun1,lboun1+mdim1-1
                vacpy(ii) = varia(ii)
             end do
             call memory_deallo(memor,vanam,vacal,varia)
             varia => vacpy
          else
             return
          end if
       else
          call memory_alloca(memor,vanam,vacal,varia,ndim1)
       end if
    end if

  end subroutine memory_resize_rp1
  
  subroutine memory_resize_rp2(memor,vanam,vacal,varia,ndim1,ndim2)
    !
    ! Real(rp)(:,:)
    !
    character(*),  intent(in)              :: vanam         !< Variable name
    character(*),  intent(in)              :: vacal         !< Calling subroutine name
    integer(8),    intent(inout)           :: memor(2)      !< Memory counter
    real(rp),      intent(inout), pointer  :: varia(:,:)
    integer(ip),   intent(in)              :: ndim1
    integer(ip),   intent(in)              :: ndim2
    real(rp),                     pointer  :: vacpy(:,:)
    integer(ip)                            :: mdim1,mdim2,ii,jj
    integer(ip)                            :: lboun1,lboun2
    
    if( ndim1 == 0_4 ) then
       call memory_deallo(memor,vanam,vacal,varia)
    else
       nullify(vacpy)
       if( associated(varia) ) then
          if( size(varia,1_ip,KIND=ip) /= ndim1 .or. size(varia,2_ip,KIND=ip) /= ndim2 ) then
             mdim1  = min(size(varia,KIND=ip,DIM=1),ndim1)
             mdim2  = min(size(varia,KIND=ip,DIM=2),ndim2)
             lboun1 = int(lbound(varia,1),ip)
             lboun2 = int(lbound(varia,2),ip)
             call memory_alloca(memor,vanam,vacal,vacpy,ndim1,ndim2,lboun1=lboun1,lboun2=lboun2)
             do jj = lboun2,lboun2+mdim2-1
                do ii = lboun1,lboun1+mdim1-1
                   vacpy(ii,jj) = varia(ii,jj)
                end do
             end do
             call memory_deallo(memor,vanam,vacal,varia)
             varia => vacpy
          else
             return
          end if
       else
          call memory_alloca(memor,vanam,vacal,varia,ndim1,ndim2)
       end if
    end if

  end subroutine memory_resize_rp2

  subroutine memory_resize_rp3(memor,vanam,vacal,varia,ndim1,ndim2,ndim3)
    !
    ! Real(rp)(:,:,:)
    !
    character(*),  intent(in)              :: vanam         !< Variable name
    character(*),  intent(in)              :: vacal         !< Calling subroutine name
    integer(8),    intent(inout)           :: memor(2)      !< Memory counter
    real(rp),      intent(inout), pointer  :: varia(:,:,:)
    integer(ip),   intent(in)              :: ndim1
    integer(ip),   intent(in)              :: ndim2
    integer(ip),   intent(in)              :: ndim3
    real(rp),                     pointer  :: vacpy(:,:,:)
    integer(ip)                            :: mdim1,mdim2,mdim3,ii,jj,kk
    integer(ip)                            :: lboun1,lboun2,lboun3
    
    if( ndim1 == 0_4 ) then
       call memory_deallo(memor,vanam,vacal,varia)
    else
       nullify(vacpy)
       if( associated(varia) ) then
          if( size(varia,1_ip,KIND=ip) /= ndim1 .or. size(varia,2_ip,KIND=ip) /= ndim2 .or. size(varia,3_ip,KIND=ip) /= ndim3 ) then
             mdim1  = min(size(varia,KIND=ip,DIM=1),ndim1)
             mdim2  = min(size(varia,KIND=ip,DIM=2),ndim2)
             mdim3  = min(size(varia,KIND=ip,DIM=3),ndim3)
             lboun1 = int(lbound(varia,1),ip)
             lboun2 = int(lbound(varia,2),ip)
             lboun3 = int(lbound(varia,3),ip)
             call memory_alloca(memor,vanam,vacal,vacpy,ndim1,ndim2,ndim3,lboun1=lboun1,lboun2=lboun2,lboun3=lboun3)
             do kk = lboun3,lboun3+mdim3-1
                do jj = lboun2,lboun2+mdim2-1
                   do ii = lboun1,lboun1+mdim1-1
                      vacpy(ii,jj,kk) = varia(ii,jj,kk)
                   end do
                end do
             end do
             call memory_deallo(memor,vanam,vacal,varia)
             varia => vacpy
          else
             return
          end if
       else
          call memory_alloca(memor,vanam,vacal,varia,ndim1,ndim2,ndim3)
       end if
    end if

  end subroutine memory_resize_rp3


  subroutine memory_alloca_xp1(memor,vanam,vacal,varia,ndim1,wzero)
    !
    ! Complex(rp)(:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    complex(rp),  intent(inout), pointer  :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    logical(lg)                         :: lzero

    if( ndim1 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       allocate( varia(ndim1) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)*2
          !lbytm = -int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = 1,ndim1
                varia(idim1) = CMPLX(0.0_rp,0.0_rp,kind=rp)
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_xp1

  subroutine memory_alloca_xp2(memor,vanam,vacal,varia,ndim1,ndim2,wzero)
    !
    ! Complex(rp)(:,:)
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    integer(ip),  intent(in)            :: ndim2
    complex(rp),  intent(inout), pointer  :: varia(:,:)
    character(*), intent(in),  optional :: wzero
    integer(4)                          :: istat
    integer(ip)                         :: idim1,idim2
    logical(lg)                         :: lzero

    if( ndim1*ndim2 > 0 ) then

       if( kfl_alloc == 1 ) call memory_deallo(memor,vanam,vacal,varia)
       if( associated(varia) ) call memory_already_associated(vanam,vacal)

       allocate( varia(ndim1,ndim2) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)*2
          !lbytm = -int(sizeof(varia),8)
          if( lzero ) then
             do idim2 = 1,ndim2
                do idim1 = 1,ndim1
                   varia(idim1,idim2) = CMPLX(0.0_rp,0.0_rp,kind=rp)
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if
       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if

  end subroutine memory_alloca_xp2

  !----------------------------------------------------------------------
  !
  ! Renumber arrays
  !
  ! 1. Nullify vacpy
  ! 2. Copy vacpy <= varia, deallocate varia
  ! 3. Reallocate varia
  ! 4. Renumber varia using vacpy
  ! 5. Deallocate vacpy
  !
  !----------------------------------------------------------------------

  subroutine memory_renumber_i41(memor,vanam,vacal,varia,lrenu)
    !
    ! Integer(4)(:)
    !
    
    character(*), intent(in)           :: vanam         !< Variable name
    character(*), intent(in)           :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)        :: memor(2)      !< Memory counter
    integer(4),   intent(inout), pointer :: varia(:)
    integer(ip),  intent(in),  pointer :: lrenu(:)
    integer(4),                pointer :: vacpy(:)
    integer(ip)                        :: idim1_new,idim1_old
    integer(ip)                        :: ndim1_new,ndim1_old

    nullify(vacpy)

    ndim1_old = memory_size(varia)
    call memory_copy(memor,vanam,vacal,varia,vacpy)
    ndim1_new = 0
    do idim1_old = 1,ndim1_old
       if( lrenu(idim1_old) /= 0 ) ndim1_new = ndim1_new + 1
    end do
    call memory_alloca(memor,vanam,vacal,varia,ndim1_new)
    do idim1_old = 1,ndim1_old
       if( lrenu(idim1_old) /= 0 ) then
          idim1_new        = lrenu(idim1_old)
          varia(idim1_new) = vacpy(idim1_old)
       end if
    end do
    call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_renumber_i41

  subroutine memory_renumber_i81(memor,vanam,vacal,varia,lrenu)
    !
    ! Integer(8)(:)
    !
    
    character(*), intent(in)           :: vanam         !< Variable name
    character(*), intent(in)           :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)        :: memor(2)      !< Memory counter
    integer(8),   intent(inout), pointer :: varia(:)
    integer(8),   intent(in),  pointer :: lrenu(:)
    integer(8),                pointer :: vacpy(:)
    integer(8)                         :: idim1_new,idim1_old
    integer(8)                         :: ndim1_new,ndim1_old

    nullify(vacpy)

    ndim1_old = int(memory_size(varia),8)
    call memory_copy(memor,vanam,vacal,varia,vacpy)
    ndim1_new = 0
    do idim1_old = 1,ndim1_old
       if( lrenu(idim1_old) /= 0 ) ndim1_new = ndim1_new + 1
    end do
    call memory_alloca(memor,vanam,vacal,varia,ndim1_new)
    do idim1_old = 1,ndim1_old
       if( lrenu(idim1_old) /= 0 ) then
          idim1_new        = lrenu(idim1_old)
          varia(idim1_new) = vacpy(idim1_old)
       end if
    end do
    call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_renumber_i81

  subroutine memory_renumber_rp1(memor,vanam,vacal,varia,lrenu)
    !
    ! Real(rp)(:)
    !
    
    character(*), intent(in)           :: vanam         !< Variable name
    character(*), intent(in)           :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)        :: memor(2)      !< Memory counter
    real(rp),     intent(inout), pointer :: varia(:)
    integer(ip),  intent(in),  pointer :: lrenu(:)
    real(rp),                  pointer :: vacpy(:)
    integer(ip)                        :: idim1_new,idim1_old
    integer(ip)                        :: ndim1_new,ndim1_old

    nullify(vacpy)

    ndim1_old = memory_size(varia)
    call memory_copy(memor,vanam,vacal,varia,vacpy)
    ndim1_new = 0
    do idim1_old = 1,ndim1_old
       if( lrenu(idim1_old) /= 0 ) ndim1_new = ndim1_new + 1
    end do
    call memory_alloca(memor,vanam,vacal,varia,ndim1_new)
    do idim1_old = 1,ndim1_old
       if( lrenu(idim1_old) /= 0 ) then
          idim1_new        = lrenu(idim1_old)
          varia(idim1_new) = vacpy(idim1_old)
       end if
    end do
    call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_renumber_rp1

  subroutine memory_renumber_rp2(memor,vanam,vacal,varia,lrenu,idime)
    !
    ! Real(rp)(:,:)
    ! If renumbering LRENU is present, it assumes by default that
    ! the second dimension is the one that should be renumbered
    !
    
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(rp),     intent(inout), pointer  :: varia(:,:)
    integer(ip),  intent(in),  pointer  :: lrenu(:)
    real(rp),                  pointer  :: vacpy(:,:)
    integer(ip),  intent(in),  optional :: idime
    integer(ip)                         :: idim1_new,idim1_old
    integer(ip)                         :: ndim1_new,ndim1_old
    integer(ip)                         :: ndim2,kdime_1,kdime_2

    nullify(vacpy)
    !
    ! If IDIME is not present: KDIME_1=2, KDIME_2=1
    !
    if( present(idime) ) then
       kdime_1 = idime
       call memory_runend('MEMORY_RENUMBER_RP2: NOT CODED')
    else
       kdime_1 = 2
    end if
    if( kdime_1 == 1 ) then
       kdime_2 = 2
    else if( kdime_1 == 2 ) then
       kdime_2 = 1
    else
       call memory_runend('MEMORY_RENUMBER_RP2: WRONG DIMENSION')
    end if
    !
    ! VARIA(KDIME_2,KDIME_1)
    ! NDIM1_OLD <= size(VARIA,2)
    !
    ndim1_old = memory_size(varia,kdime_1)
    ndim2     = memory_size(varia,kdime_2)
    call memory_copy(memor,vanam,vacal,varia,vacpy)
    ndim1_new = 0
    do idim1_old = 1,ndim1_old
       if( lrenu(idim1_old) /= 0 ) ndim1_new = ndim1_new + 1
    end do
    call memory_alloca(memor,vanam,vacal,varia,ndim2,ndim1_new)
    do idim1_old = 1,ndim1_old
       if( lrenu(idim1_old) /= 0 ) then
          idim1_new = lrenu(idim1_old)
          varia(1:ndim2,idim1_new) = vacpy(1:ndim2,idim1_old)
       end if
    end do
    call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_renumber_rp2

  subroutine memory_renumber_rp3(memor,vanam,vacal,varia,lrenu,idime)
    !
    ! Real(rp)(:,:)
    ! If renumbering LRENU is present, it assumes by default that
    ! the second dimension is the one that should be renumbered
    !
    
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    real(rp),     intent(inout), pointer  :: varia(:,:,:)
    integer(ip),  intent(in),    pointer  :: lrenu(:)
    real(rp),                    pointer  :: vacpy(:,:,:)
    integer(ip),  intent(in),    optional :: idime
    integer(ip)                           :: idim1_new,idim1_old
    integer(ip)                           :: ndim1_new,ndim1_old
    integer(ip)                           :: ndim2,kdime_1,kdime_2
    integer(ip)                           :: ndim3

    nullify(vacpy)
    !
    ! If IDIME is not present: KDIME_1=2, KDIME_2=1
    !
    if( present(idime) ) then
       kdime_1 = idime
       call memory_runend('RENPR3 NOT CODED')
    else
       kdime_1 = 2
    end if
    if( kdime_1 == 1 ) then
       kdime_2 = 2
    else if( kdime_1 == 2 ) then
       kdime_2 = 1
    else
       call memory_runend('MEMORY_RENUMBER_RP2: WRONG DIMENSION')
    end if
    !
    ! VARIA(KDIME_2,KDIME_1)
    ! NDIM1_OLD <= size(VARIA,2)
    !
    ndim1_old = size(varia,kdime_1,kind=ip)
    ndim2     = size(varia,kdime_2,kind=ip)
    ndim3     = size(varia,3,kind=ip)
    call memory_copy(memor,vanam,vacal,varia,vacpy)
    ndim1_new = 0
    do idim1_old = 1,ndim1_old
       if( lrenu(idim1_old) /= 0 ) ndim1_new = ndim1_new + 1
    end do
    call memory_alloca(memor,vanam,vacal,varia,ndim2,ndim1_new,ndim3)
    do idim1_old = 1,ndim1_old
       if( lrenu(idim1_old) /= 0 ) then
          idim1_new = lrenu(idim1_old)
          varia(1:ndim2,idim1_new,1:ndim3) = vacpy(1:ndim2,idim1_old,1:ndim3)
       end if
    end do
    call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_renumber_rp3

  !----------------------------------------------------------------------
  !
  ! Allocate minimum size if pointers are not associated
  !
  !----------------------------------------------------------------------
  !
  ! Real(rp)(:)
  !
  subroutine memory_alloca_min_rp1(memor,vanam,vacal,varia)
    
    integer(8),   intent(inout)          :: memor(2)      !< Memory counter
    character(*), intent(in)             :: vanam         !< Variable name
    character(*), intent(in)             :: vacal         !< Calling subroutine name
    real(rp),     intent(inout), pointer :: varia(:)

    if( .not. associated(varia) ) then
       call memory_alloca(memor,vanam,vacal,varia,1_ip)
    end if
    
  end subroutine memory_alloca_min_rp1
  ! 
  ! Real(rp)(:,:)
  !
  subroutine memory_alloca_min_rp2(memor,vanam,vacal,varia)
    
    integer(8),   intent(inout)          :: memor(2)      !< Memory counter
    character(*), intent(in)             :: vanam         !< Variable name
    character(*), intent(in)             :: vacal         !< Calling subroutine name
    real(rp),     intent(inout), pointer :: varia(:,:)
    if( .not. associated(varia) ) then
       call memory_alloca(memor,vanam,vacal,varia,1_ip,1_ip)
    end if
  end subroutine memory_alloca_min_rp2
  !
  ! Real(rp)(:,:,:)
  !
  subroutine memory_alloca_min_rp3(memor,vanam,vacal,varia)
    
    integer(8),   intent(inout)          :: memor(2)      !< Memory counter
    character(*), intent(in)             :: vanam         !< Variable name
    character(*), intent(in)             :: vacal         !< Calling subroutine name
    real(rp),     intent(inout), pointer :: varia(:,:,:)
    if( .not. associated(varia) ) then
       call memory_alloca(memor,vanam,vacal,varia,1_ip,1_ip,1_ip)
    end if
  end subroutine memory_alloca_min_rp3
  !
  ! Real(rp)(:,:,:,:)
  !
  subroutine memory_alloca_min_rp4(memor,vanam,vacal,varia)
    
    integer(8),   intent(inout)          :: memor(2)      !< Memory counter
    character(*), intent(in)             :: vanam         !< Variable name
    character(*), intent(in)             :: vacal         !< Calling subroutine name
    real(rp),     intent(inout), pointer :: varia(:,:,:,:)
    if( .not. associated(varia) ) then
       call memory_alloca(memor,vanam,vacal,varia,1_ip,1_ip,1_ip,1_ip)
    end if
  end subroutine memory_alloca_min_rp4
  !
  ! integer(ip)(:)
  !
  subroutine memory_alloca_min_ip1(memor,vanam,vacal,varia)
    
    integer(8),   intent(inout)          :: memor(2)      !< Memory counter
    character(*), intent(in)             :: vanam         !< Variable name
    character(*), intent(in)             :: vacal         !< Calling subroutine name
    integer(ip),  intent(inout), pointer :: varia(:)
    if( .not. associated(varia) ) then
       call memory_alloca(memor,vanam,vacal,varia,1_ip)
    end if
  end subroutine memory_alloca_min_ip1
  !
  ! integer(ip)(:,:)
  !
  subroutine memory_alloca_min_ip2(memor,vanam,vacal,varia)
    
    integer(8),   intent(inout)          :: memor(2)      !< Memory counter
    character(*), intent(in)             :: vanam         !< Variable name
    character(*), intent(in)             :: vacal         !< Calling subroutine name
    integer(ip),  intent(inout), pointer :: varia(:,:)
    if( .not. associated(varia) ) then
       call memory_alloca(memor,vanam,vacal,varia,1_ip,1_ip)
    end if
  end subroutine memory_alloca_min_ip2
  !
  ! integer(ip)(:,:,:)
  !
  subroutine memory_alloca_min_ip3(memor,vanam,vacal,varia)
    
    integer(8),   intent(inout)          :: memor(2)      !< Memory counter
    character(*), intent(in)             :: vanam         !< Variable name
    character(*), intent(in)             :: vacal         !< Calling subroutine name
    integer(ip),  intent(inout), pointer :: varia(:,:,:)
    if( .not. associated(varia) ) then
       call memory_alloca(memor,vanam,vacal,varia,1_ip,1_ip,1_ip)
    end if
  end subroutine memory_alloca_min_ip3
  !
  ! logical(lg)(:)
  !
  subroutine memory_alloca_min_lg1(memor,vanam,vacal,varia)
    
    integer(8),   intent(inout)          :: memor(2)      !< Memory counter
    character(*), intent(in)             :: vanam         !< Variable name
    character(*), intent(in)             :: vacal         !< Calling subroutine name
    logical(lg),  intent(inout), pointer :: varia(:)
    if( .not. associated(varia) ) then
       call memory_alloca(memor,vanam,vacal,varia,1_ip)
    end if
  end subroutine memory_alloca_min_lg1
  !
  ! logical(lg)(:,:)
  !
  subroutine memory_alloca_min_lg2(memor,vanam,vacal,varia)
    
    integer(8),   intent(inout)          :: memor(2)      !< Memory counter
    character(*), intent(in)             :: vanam         !< Variable name
    character(*), intent(in)             :: vacal         !< Calling subroutine name
    logical(lg),  intent(inout), pointer :: varia(:,:)
    if( .not. associated(varia) ) then
       call memory_alloca(memor,vanam,vacal,varia,1_ip,1_ip)
    end if
  end subroutine memory_alloca_min_lg2
  !
  ! logical(lg)(:,:,:)
  !
  subroutine memory_alloca_min_lg3(memor,vanam,vacal,varia)
    
    integer(8),   intent(inout)          :: memor(2)      !< Memory counter
    character(*), intent(in)             :: vanam         !< Variable name
    character(*), intent(in)             :: vacal         !< Calling subroutine name
    logical(lg),  intent(inout), pointer :: varia(:,:,:)
    if( .not. associated(varia) ) then
       call memory_alloca(memor,vanam,vacal,varia,1_ip,1_ip,1_ip)
    end if
  end subroutine memory_alloca_min_lg3

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    19/11/2015
  !> @brief   Size of some arrays
  !> @details This routine computes the total size of arrays.
  !>          If the pointer is nullified, it returns 0.
  !>
  !-----------------------------------------------------------------------

  pure function memory_size_rp1(varia,ndim1)
    !
    ! Real(rp)(:)
    !
    
    real(rp),     intent(in),    pointer  :: varia(:)
    integer(ip),  intent(in), optional    :: ndim1
    integer(ip)                           :: memory_size_rp1

    if( associated(varia) ) then
       memory_size_rp1 = size(varia,kind=ip)
    else
       memory_size_rp1 = 0_ip
    end if

  end function memory_size_rp1

  pure function memory_size_rp2(varia,ndim1)
    !
    ! Real(rp)(:,:)
    !
    
    real(rp),     intent(in),    pointer  :: varia(:,:)
    integer(ip),  intent(in),    optional :: ndim1
    integer(ip)                           :: memory_size_rp2

    if( associated(varia) ) then
       if( present(ndim1) ) then
          if( ndim1 >= 1 .and. ndim1 <= 2 ) then
             memory_size_rp2 = size(varia,ndim1,kind=ip)
          end if
       else
          memory_size_rp2 = size(varia,kind=ip)
       end if
    else
       memory_size_rp2 = 0_ip
    end if

  end function memory_size_rp2

  pure function memory_size_rp3(varia,ndim1)
    !
    ! Real(rp)(:,:,:)
    !
    
    real(rp),     intent(in),    pointer  :: varia(:,:,:)
    integer(ip),  intent(in),    optional :: ndim1
    integer(ip)                           :: memory_size_rp3

    if( associated(varia) ) then
       if( present(ndim1) ) then
          if( ndim1 >= 1 .and. ndim1 <= 3 ) then
             memory_size_rp3 = size(varia,ndim1,kind=ip)
          end if
       else
          memory_size_rp3 = size(varia,kind=ip)
       end if
    else
       memory_size_rp3 = 0_ip
    end if

  end function memory_size_rp3

  pure function memory_size_rp4(varia,ndim1)
    !
    ! Real(rp)(:,:,:)
    !
    
    real(rp),     intent(in),    pointer  :: varia(:,:,:,:)
    integer(ip),  intent(in),    optional :: ndim1
    integer(ip)                           :: memory_size_rp4

    if( associated(varia) ) then
       if( present(ndim1) ) then
          if( ndim1 >= 1 .and. ndim1 <= 4 ) then
             memory_size_rp4 = size(varia,ndim1,kind=ip)
          else
             memory_size_rp4 = 0_ip
          end if
       else
          memory_size_rp4 = size(varia,kind=ip)
       end if
    else
       memory_size_rp4 = 0_ip
    end if

  end function memory_size_rp4

  pure function memory_size_ip1_4(varia,ndim1)
    !
    ! Int(ip)(:)
    !
    
    integer(4),   intent(in), pointer  :: varia(:)
    integer(4),   intent(in), optional :: ndim1
    integer(ip)                        :: memory_size_ip1_4

    if( associated(varia) ) then
       memory_size_ip1_4 = size(varia,kind=ip)
    else
       memory_size_ip1_4 = 0_ip
    end if

  end function memory_size_ip1_4

  pure function memory_size_ip1_8(varia,ndim1)
    !
    ! Int(ip)(:)
    !
    
    integer(8),   intent(in), pointer  :: varia(:)
    integer(8),   intent(in), optional :: ndim1
    integer(ip)                        :: memory_size_ip1_8

    if( associated(varia) ) then
       memory_size_ip1_8 = size(varia,kind=ip)
    else
       memory_size_ip1_8 = 0_ip
    end if

  end function memory_size_ip1_8

  pure function memory_size_ip2_4(varia,ndim1)
    !
    ! Int(ip)(:,:)
    !
    
    integer(4),   intent(in), pointer  :: varia(:,:)
    integer(ip),  intent(in), optional :: ndim1
    integer(ip)                        :: memory_size_ip2_4

    if( associated(varia) ) then
       if( present(ndim1) ) then
          if( ndim1 >= 1 .and. ndim1 <= 2 ) then
             memory_size_ip2_4 = size(varia,ndim1,kind=ip)
          else
             memory_size_ip2_4 = 0_ip
          end if
       else
          memory_size_ip2_4 = size(varia,kind=ip)
       end if
    else
       memory_size_ip2_4 = 0_ip
    end if

  end function memory_size_ip2_4

  pure function memory_size_ip2_8(varia,ndim1)
    !
    ! Int(ip)(:,:)
    !
    
    integer(8),   intent(in), pointer  :: varia(:,:)
    integer(ip),  intent(in), optional :: ndim1
    integer(ip)                        :: memory_size_ip2_8

    if( associated(varia) ) then
       if( present(ndim1) ) then
          if( ndim1 >= 1 .and. ndim1 <= 2 ) then
             memory_size_ip2_8 = size(varia,ndim1,kind=ip)
          else
             memory_size_ip2_8 = 0_ip             
          end if
       else
          memory_size_ip2_8 = size(varia,kind=ip)
       end if
    else
       memory_size_ip2_8 = 0_ip
    end if

  end function memory_size_ip2_8

  pure function memory_size_ip3_4(varia,ndim1)
    !
    ! Int(ip)(:,:,:)
    !
    
    integer(4),   intent(in), pointer  :: varia(:,:,:)
    integer(ip),  intent(in), optional :: ndim1
    integer(ip)                        :: memory_size_ip3_4

    if( associated(varia) ) then
       if( present(ndim1) ) then
          if( ndim1 >= 1 .and. ndim1 <= 3 ) then
             memory_size_ip3_4 = size(varia,ndim1,kind=ip)
          else
             memory_size_ip3_4 = 0_ip             
          end if
       else
          memory_size_ip3_4 = size(varia,kind=ip)
       end if
    else
       memory_size_ip3_4 = 0_ip
    end if

  end function memory_size_ip3_4

  pure function memory_size_ip3_8(varia,ndim1)
    !
    ! Int(ip)(:,:,:)
    !
    
    integer(8),   intent(in), pointer  :: varia(:,:,:)
    integer(ip),  intent(in), optional :: ndim1
    integer(ip)                        :: memory_size_ip3_8

    if( associated(varia) ) then
       if( present(ndim1) ) then
          if( ndim1 >= 1 .and. ndim1 <= 3 ) then
             memory_size_ip3_8 = size(varia,ndim1,kind=ip)
          else
             memory_size_ip3_8 = 0_ip             
          end if
       else
          memory_size_ip3_8 = size(varia,kind=ip)
       end if
    else
       memory_size_ip3_8 = 0_ip
    end if

  end function memory_size_ip3_8

  pure function memory_size_lg(varia)
    !
    ! Log(lg)(:)
    !
    
    logical(lg),  intent(in), pointer :: varia(:)
    integer(ip)                       :: memory_size_lg

    if( associated(varia) ) then
       memory_size_lg = int(size(varia,1,ip),ip)
    else
       memory_size_lg = 0_ip
    end if

  end function memory_size_lg

  pure function memory_size_i1p_1(varia,ndime)
    !
    ! type(i1p)(:)
    !
    
    type(i1p),     intent(in),   pointer  :: varia(:)
    integer(ip),   intent(in),   optional :: ndime
    integer(ip)                           :: memory_size_i1p_1
    integer(ip)                           :: ii

    if( present(ndime) ) then
       
       memory_size_i1p_1 = 0_ip
       
       if( associated(varia) ) then
          do ii = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
             if( associated(varia(ii) % l) ) &
                  memory_size_i1p_1 = max(memory_size_i1p_1,size(varia(ii) % l,DIM=ndime,KIND=ip))
          end do
       end if
    else
       if( associated(varia) ) then
          memory_size_i1p_1 = size(varia,KIND=ip)
       else
          memory_size_i1p_1 = 0_ip
       end if
    end if
    
  end function memory_size_i1p_1

  pure function memory_size_r1p_1(varia,ndime)
    !
    ! type(r1p)(:)
    !
    
    type(r1p),     intent(in),   pointer  :: varia(:)
    integer(ip),   intent(in),   optional :: ndime
    integer(ip)                           :: memory_size_r1p_1
    integer(ip)                           :: ii

    if( present(ndime) ) then
       
       memory_size_r1p_1 = 0_ip
       
       if( associated(varia) ) then
          do ii = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
             if( associated(varia(ii) % a) ) &
                  memory_size_r1p_1 = max(memory_size_r1p_1,size(varia(ii) % a,DIM=ndime,KIND=ip))
          end do
       end if
    else
       if( associated(varia) ) then
          memory_size_r1p_1 = size(varia,KIND=ip)
       else
          memory_size_r1p_1 = 0_ip
       end if
    end if
    
  end function memory_size_r1p_1

  function memory_size_r2p_1(varia,ndime)
    !
    ! type(r2p)(:)
    !
    
    type(r2p),     intent(in),   pointer  :: varia(:)
    integer(ip),   intent(in),   optional :: ndime
    integer(ip)                           :: memory_size_r2p_1
    integer(ip)                           :: ii

    if( present(ndime) ) then
       
       memory_size_r2p_1 = 0_ip
       
       if( associated(varia) ) then
          do ii = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
             if( associated(varia(ii) % a) ) &
                  memory_size_r2p_1 = max(memory_size_r2p_1,size(varia(ii) % a,DIM=ndime,KIND=ip))
          end do
       end if
    else
       if( associated(varia) ) then
          memory_size_r2p_1 = size(varia,KIND=ip)
       else
          memory_size_r2p_1 = 0_ip
       end if
    end if
    
  end function memory_size_r2p_1

  function memory_size_r3p_1(varia,ndime)
    !
    ! type(r3p)(:)
    !
    
    type(r3p),     intent(in),   pointer  :: varia(:)
    integer(ip),   intent(in),   optional :: ndime
    integer(ip)                           :: memory_size_r3p_1
    integer(ip)                           :: ii

    if( present(ndime) ) then
       
       memory_size_r3p_1 = 0_ip
       
       if( associated(varia) ) then
          do ii = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
             if( associated(varia(ii) % a) ) &
                  memory_size_r3p_1 = max(memory_size_r3p_1,size(varia(ii) % a,DIM=ndime,KIND=ip))
          end do
       end if
    else
       if( associated(varia) ) then
          memory_size_r3p_1 = size(varia,KIND=ip)
       else
          memory_size_r3p_1 = 0_ip
       end if
    end if
    
  end function memory_size_r3p_1
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-08-09
  !> @brief   Variable initialization
  !> @details Initialization of a variable
  !> 
  !-----------------------------------------------------------------------
  
  subroutine memory_initia_rp1(varia,VARIABLE_VALUE,lboun1,uboun1)
    !
    ! Real(rp)(:)
    !
    real(rp),   intent(inout), pointer          :: varia(:)
    real(rp),   intent(in),            optional :: VARIABLE_VALUE
    integer(ip),intent(in),            optional :: lboun1
    integer(ip),intent(in),            optional :: uboun1
    integer(ip)                                 :: idim1
    integer(ip)                                 :: my_lboun1
    integer(ip)                                 :: my_uboun1

    if( associated(varia) ) then
       if( present(lboun1) ) then
          my_lboun1 = lboun1
       else
          my_lboun1 = int(lbound(varia,1),ip)
       end if
       if( present(uboun1) ) then
          my_uboun1 = uboun1
       else
          my_uboun1 = int(ubound(varia,1),ip)
       end if
       if( present(VARIABLE_VALUE) ) then
          do idim1 = my_lboun1,my_uboun1
             varia(idim1) = VARIABLE_VALUE
          end do
       else
          do idim1 =  my_lboun1,my_uboun1
             varia(idim1) = 0.0_rp
          end do
       end if
    end if
    
  end subroutine memory_initia_rp1
  
  subroutine memory_initia_rp2(varia,VARIABLE_VALUE)
    !
    ! Real(rp)(:,:)
    !
    real(rp),   intent(inout), pointer          :: varia(:,:)
    real(rp),   intent(in),            optional :: VARIABLE_VALUE
    integer(ip)                                 :: idim1,idim2
    
    if( associated(varia) ) then
       if( present(VARIABLE_VALUE) ) then
          do idim2 = int(lbound(varia,2),ip),int(ubound(varia,2),ip)
             do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
                varia(idim1,idim2) = VARIABLE_VALUE
             end do
          end do
       else
          do idim2 = int(lbound(varia,2),ip),int(ubound(varia,2),ip)
             do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
                varia(idim1,idim2) = 0.0_rp
             end do
          end do
       end if
    end if
    
  end subroutine memory_initia_rp2
  
  subroutine memory_initia_rp3(varia,VARIABLE_VALUE)
    !
    ! Real(rp)(:,:)
    !
    real(rp),   intent(inout), pointer          :: varia(:,:,:)
    real(rp),   intent(in),            optional :: VARIABLE_VALUE
    integer(ip)                                 :: idim1,idim2,idim3
    integer(ip)                                 :: ldim1,udim1
    integer(ip)                                 :: ldim2,udim2
    integer(ip)                                 :: ldim3,udim3
    real(rp)                                    :: rvalu
    
    if( associated(varia) ) then
       ldim1 = int(lbound(varia,1),ip)
       ldim2 = int(lbound(varia,2),ip)
       ldim3 = int(lbound(varia,3),ip)
       udim1 = int(ubound(varia,1),ip)
       udim2 = int(ubound(varia,2),ip)
       udim3 = int(ubound(varia,3),ip)
       if( present(VARIABLE_VALUE) ) then
          rvalu = VARIABLE_VALUE
       else
          rvalu = 0.0_rp
       end if      
       do idim3 = ldim3,udim3
          do idim2 = ldim2,udim2
             do idim1 = ldim1,udim1
                varia(idim1,idim2,idim3) = rvalu
             end do
          end do
       end do
    end if
    
  end subroutine memory_initia_rp3

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-05-25
  !> @brief   Add
  !> @details Add values to VARIA <= VARIA = VARIA + ADD
  !> 
  !-----------------------------------------------------------------------
  
  subroutine memory_append_ip1(memor,vanam,vacal,varia,vaadd)
    !
    ! Integer(ip)(:)
    !  
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(ip),  intent(inout), pointer  :: varia(:)
    integer(ip),  intent(in),    pointer  :: vaadd(:)
    integer(ip)                           :: ndim1,idim1,jj
    integer(ip),                 pointer  :: vatot(:)

    if( memory_size(vaadd) > 0 ) then
       
       ndim1 = memory_size(varia) + memory_size(vaadd)
       
       if( ndim1 > 0 ) then
          
          nullify(vatot)
          if( associated(varia) ) then
             call memory_alloca(memor,vanam,vacal,vatot,ndim1,'DO_NOT_INITIALIZE',lboun=int(lbound(varia,1),ip))
             do idim1 = int(lbound(varia,1),ip),int(ubound(varia,1),ip)
                vatot(idim1) = varia(idim1)
             end do
             jj = int(ubound(varia,1),ip)
          else
             call memory_alloca(memor,vanam,vacal,vatot,ndim1,'DO_NOT_INITIALIZE',lboun=int(lbound(vaadd,1),ip))
             jj = int(lbound(vaadd,1),ip)-1
          end if
          call memory_deallo(memor,vanam,vacal,varia)
          do idim1 = int(lbound(vaadd,1),ip),int(ubound(vaadd,1),ip)
             jj = jj + 1
             vatot(jj) = vaadd(idim1)
          end do
          call memory_deallo(memor,vanam,vacal,varia)
          varia => vatot
          
       end if
       
    end if
    
  end subroutine memory_append_ip1
  
end module mod_memory_basic
!> @}

