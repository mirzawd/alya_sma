!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @defgroup Hash_Table_Toolbox
!> Toolbox for Hash table operations
!> @{
!> @file    mod_htable.f90
!> @author  Ricard Borrell
!> @date    2018-11-27
!> @brief   Htable
!> @details Toolbox for hash table
!-----------------------------------------------------------------------

module mod_htable

  use def_kintyp_basic, only : ip,lg,i1p
  use mod_memory_basic, only : memory_size
  use mod_memory_basic, only : memory_alloca
  use mod_memory_basic, only : memory_deallo
  use mod_memory_basic, only : memory_resize
  implicit none

  private
  
  integer(ip), parameter :: size_factor=10
  integer(ip), parameter :: NIL=-1
  integer(ip), parameter :: HTablePrimes(125) = &
       (/      &
       &  7,         11,        19,       31,       127,    &
       &  211,       383,       631,      887,      1151,   &
       &  1399,      1663,      1913,     2161,     2423,   &
       &  2687,      2939,      3191,     3449,     3709,   &
       &  3967,      4219,      4463,     4733,     4987,   &
       &  5237,      5503,      5749,     6011,     6271,   &
       &  6521,      6781,      7039,     7283,     7549,   &
       &  7793,      8059,      8317,     8573,     8831,   &
       &  9067,      9343,      9587,     9851,     10111,  &
       &  10357,     10613,     10867,    11131,    11383,  &
       &  11633,     11903,     12157,    12413,    12671,  &
       &  12923,     13183,     13421,    13693,    13933,  &
       &  14207,     14461,     14717,    14969,    15227,  &
       &  15473,     15739,     15991,    16253,    16493,  &
       &  18553,     20599,     22651,    24697,    26737,  &
       &  28793,     30841,     32887,    34939,    36979,  &
       &  39023,     41081,     43133,    45181,    47221,  &
       &  49279,     51307,     53359,    55411,    57467,  &
       &  59513,     61561,     63611,    65657,    82039,  &
       &  98429,     114809,    131171,   147583,   163927, &
       &  180347,    196727,    213119,   229499,   245881, &
       &  262271,    327799,    393331,   458879,   524413, &
       &  589933,    655471,    721013,   786553,   852079, &
       &  917629,    983153,   1159523,  2750159,  4256233, &
       & 7368787,  10570841,  15485863, 20495843, 40003409  /)

  type hash_t
     integer(ip)          :: sizei
     integer(ip)          :: sizet
     integer(ip)          :: nelem
     integer(ip), pointer :: dades(:)
     integer(ip), pointer :: lids(:)
     integer(ip), pointer :: indic(:)
     type(i1p),   pointer :: list(:)
     logical(lg)          :: lidson
     logical(lg)          :: repeated
   contains
     procedure,   pass    :: init
  end type hash_t
  integer(8)              :: memor_ht(2)
  integer(8)              :: memor_loc(2)

  interface htaadd
     module procedure htaadd_0,htaadd_1,htaadd_1b
  end interface htaadd
  
  interface htable_initialization
     module procedure htable_initialization_0,htable_initialization_1
  end interface htable_initialization

  public :: hash_t                   ! Hash table type
  public :: htaini                   ! Initialize and allocate hasth table
  public :: htades                   ! Deallocate hash table
  public :: htades_list              ! Deallocate hash table
  public :: htares                   ! Restart the hash table
  public :: HTableMaxPrimeNumber     ! Give the Max prime number
  public :: htaadd                   ! Add an element to the hash table
  public :: htalid                   ! Returns local ID for an specific global ID
  public :: htalid_list              ! Returns local ID for an specific global ID
  public :: htable_initialization    ! Initialize a hash table
  
contains

  !-----------------------------------------------------------------------
  !
  !> @brief   Initialize hash table
  !> @details The size of ht % dades is evaluated with the subroutine
  !>          HTablePrimeNumber
  !
  !-----------------------------------------------------------------------

  subroutine init(ht)
    class(hash_t), intent(inout) :: ht
    call htable_initialization_1(ht)
    
  end subroutine init
  
  subroutine htable_initialization_0()

    memor_ht = 0_8
   
  end subroutine htable_initialization_0
  
  subroutine htable_initialization_1(ht)
    type(hash_t), intent(inout) :: ht

    ht % sizei    = 0_ip
    ht % sizet    = 0_ip
    ht % nelem    = 0_ip
    ht % lidson   = .false.
    ht % repeated = .false.
    nullify(ht % dades)
    nullify(ht % lids)
    nullify(ht % indic)
    nullify(ht % list)

  end subroutine htable_initialization_1
  
  !-----------------------------------------------------------------------
  !
  !> @brief   Initialize hash table
  !> @details The size of ht % dades is evaluated with the subroutine
  !>          HTablePrimeNumber
  !
  !-----------------------------------------------------------------------
  
  subroutine htaini(ht, tsize, lidson, AUTOMATIC_SIZE,memor_opt, REPEATED_ELEMENTS)

    implicit none
    type(hash_t), intent(inout)           :: ht
    integer(ip),  intent(in)              :: tsize
    logical(lg),  intent(in),    optional :: lidson
    logical(lg),  intent(in),    optional :: AUTOMATIC_SIZE
    integer(8),   intent(inout), optional :: memor_opt(2)
    logical(lg),  intent(in),    optional :: REPEATED_ELEMENTS

    if( present(memor_opt) ) then
       memor_loc = memor_opt
    else
       memor_loc = memor_ht
    end if
    if( present(REPEATED_ELEMENTS) ) then
       ht % repeated = REPEATED_ELEMENTS
    else
       ht % repeated = .false.
    end if
    
    nullify(ht % dades)
    nullify(ht % lids)
    nullify(ht % list)
    nullify(ht % indic)

    ht % sizet = HTablePrimeNumber( tsize, AUTOMATIC_SIZE , ht )
    ht % nelem = 0_ip
    call memory_alloca(memor_loc,'HT % DADES','mod_htable',ht % dades,ht % sizet)
    ht % dades(1:ht%sizet) = NIL
    if(present(lidson)) then
       ht % lidson = lidson
    else
       ht % lidson = .false.       
    end if
    if( ht % lidson ) then
       call memory_alloca(memor_loc,'HT % LIDS','mod_htable',ht % lids,ht % sizet)
    endif
    if( ht % repeated ) then
       call memory_alloca(memor_loc,'HT % LIST','mod_htable',ht % list,ht % sizet)
    end if

    if( present(memor_opt) ) memor_opt = memor_loc

  end subroutine htaini

  !-----------------------------------------------------------------------
  !
  !> @brief   Deallocate hash table
  !> @details Deallocate hash table
  !
  !-----------------------------------------------------------------------
  
  subroutine htades( ht ,memor_opt )
 
    implicit none
    type(hash_t), intent(inout)             :: ht
    integer(8),   intent(inout),   optional :: memor_opt(2)
    
    if( present(memor_opt) ) then
       memor_loc = memor_opt
    else
       memor_loc = memor_ht
    end if

    call memory_deallo(memor_loc,'HT % DADES','mod_htable',ht % dades)
    call memory_deallo(memor_loc,'HT % LIDS' ,'mod_htable',ht % lids)
    call memory_deallo(memor_loc,'HT % LIST' ,'mod_htable',ht % list)
    nullify(ht % indic)

    if( present(memor_opt) ) then
       memor_opt = memor_loc
    else
       memor_ht  = memor_loc
    end if
    
  end subroutine htades

  subroutine htades_list(local_id,memor_opt) 

    integer(8),   intent(inout), optional :: memor_opt(2)
    integer(ip),  intent(inout), pointer  :: local_id(:)

   if( present(memor_opt) ) then
       memor_loc = memor_opt
    else
       memor_loc = memor_ht
    end if

    call memory_deallo(memor_loc,'LOCAL_ID','mod_htable',local_id)
    
    if( present(memor_opt) ) then
       memor_opt = memor_loc
    else
       memor_ht  = memor_loc
    end if
    
  end subroutine htades_list

  !-----------------------------------------------------------------------
  !
  !> @brief   Restart hash table
  !> @details Restart hash table
  !
  !-----------------------------------------------------------------------
  
  subroutine htares( ht, indic )

    implicit none
    type(hash_t), intent(inout)                  :: ht
    integer(ip),  intent(in),  target, optional  :: indic(:)

    ht % nelem = 0_ip
    if(present(indic)) ht % indic => indic
    ht % dades(1:ht%sizet) = NIL
    if(ht % lidson) ht % lids(1:ht%sizet) = 0_ip

  end subroutine htares

  !-----------------------------------------------------------------------
  !
  !> @brief   Add element into the table
  !> @details Add element into the table
  !
  !-----------------------------------------------------------------------
  
  subroutine htaadd_0( ht, valin, lid, isin, memor_opt, icont)

    type(hash_t), intent(inout)                :: ht
    integer(ip),  intent(in)                   :: valin 
    integer(ip),  intent(out),   optional      :: lid 
    logical(lg),  intent(out),   optional      :: isin
    integer(8),   intent(inout), optional      :: memor_opt(2)
    integer(ip),  intent(in),    optional      :: icont 
    logical(lg)                                :: found
    integer(ip)                                :: pos, incr, val, ii, kk
    logical(lg)                                :: noEncontrado

    if( present(memor_opt) ) then
       memor_loc = memor_opt
    else
       memor_loc = memor_ht
    end if

    pos   = MOD(valin,ht%sizet) + 1_ip
    incr  = MOD(valin,ht%sizet  - 1_ip)

    noEncontrado = .true.
    found        = .false.
    ii           = 0
    
    do while(noEncontrado)
       val = ht%dades(pos)
       if (val==NIL) then
          ht % nelem           = ht % nelem + 1
          if( ht % nelem > ht % sizet ) call runend('HTAADD: HASH TABLE IS TOO SMALL')
          ht % dades(pos)      = valin
          if(associated(ht % indic)) ht % indic(ht%nelem) = valin
          if(ht % lidson) ht % lids(pos) = ht % nelem
          noEncontrado = .false.
       else if (val==valin) then
          noEncontrado = .false.
          found        = .true.
       else
          if( ii > ht%sizet ) call runend('YOUR HASH TABLE IS TOO SMALL!')
          pos = MOD(pos+incr,ht%sizet) + 1
          ii  = ii + 1
       endif
    enddo
    if(present(lid) ) lid  = ht % lids(pos)
    if(present(isin)) isin = found
    if( ht%repeated ) then
       if(.not. present(icont)) call runend('HTAADD: NOT PRESENT ICONT FOR REPEATED')
       kk = memory_size(ht % list(pos) % l)
       if( kk == 0 ) then
          call memory_alloca(memor_ht,'MY_LIST_MASTER_NODES','mod_periodicity',ht % list(pos) % l,1_ip)
       else
          call memory_resize(memor_ht,'MY_LIST_MASTER_NODES','mod_periodicity',ht % list(pos) % l,kk+1_ip)
       end if
       ht % list(pos) % l(kk+1) = icont
    end if
    
    if( present(memor_opt) ) then
       memor_opt = memor_loc
    else
       memor_ht  = memor_loc
    end if
    
  end subroutine htaadd_0

  subroutine htreallocate(ht,memor_opt)
    
    type(hash_t), intent(inout)           :: ht
    integer(8),   intent(inout), optional :: memor_opt(2)
    integer(ip)                           :: tsize
    
    tsize = 2_ip * ht % sizei
    call htades(ht, memor_opt)
    call htaini(ht, tsize, lidson=ht % lidson, AUTOMATIC_SIZE=.false.,memor_opt=memor_opt, REPEATED_ELEMENTS=ht % repeated)
    
  end subroutine htreallocate

  !-----------------------------------------------------------------------
  !
  !> @brief   Add array of elements into the table
  !> @details Add array of elements into the table
  !
  !-----------------------------------------------------------------------
  
  subroutine htaadd_1( ht, valin, memor_opt)

    implicit none
    type(hash_t),          intent(inout)           :: ht
    integer(ip),  pointer, intent(in)              :: valin(:)
    integer(8),            intent(inout), optional :: memor_opt(2)
    integer(ip)                                    :: ii

    do ii=1,memory_size(valin)
       call htaadd_0(ht,valin(ii),memor_opt=memor_opt,icont=ii)
    end do
    
  end subroutine htaadd_1
  
  !-----------------------------------------------------------------------
  !
  !> @brief   Add array of elements into the table
  !> @details Add array of elements into the table
  !
  !-----------------------------------------------------------------------
  
  subroutine htaadd_1b( ht, nvalin, valin, memor_opt)

    implicit none
    type(hash_t),  intent(inout)           :: ht
    integer(ip),   intent(in)              :: nvalin
    integer(ip),   intent(in)              :: valin(*)
    integer(8),    intent(inout), optional :: memor_opt(2)
    integer(ip)                            :: ii

    do ii=1,nvalin
       call htaadd_0(ht,valin(ii),memor_opt=memor_opt,icont=ii)
    end do
    
  end subroutine htaadd_1b
  
  !-----------------------------------------------------------------------
  !
  !> @brief   Returns lid for an specific valin
  !> @details Returns lid for an specific valin. If valin is not in 
  !>          the table returns 0. The lis is the position of valin
  !>          in the array idic(:)
  !
  !-----------------------------------------------------------------------

  function htalid_list(ht, valin, memor_opt) result(local_id)

    integer(ip),  pointer                 :: local_id(:)
    type(hash_t), intent(in)              :: ht
    integer(ip),  intent(in)              :: valin 
    integer(8),   intent(inout), optional :: memor_opt(2)
    integer(ip)                           :: pos, incr, val,kk
    logical(lg)                           :: noEncontrado

    if( present(memor_opt) ) then
       memor_loc = memor_opt
    else
       memor_loc = memor_ht
    end if

    if( .not. ht % lidson ) call runend("htalid: generate hast_t with lidson=true")
    pos   = MOD(valin,ht%sizet) + 1_ip
    incr  = MOD(valin,ht%sizet  - 1_ip)
    nullify(local_id)
    
    noEncontrado = .true.
    do while(noEncontrado)
       val = ht%dades(pos)
       if (val==NIL) then
          pos = -1_ip
          noEncontrado = .false.
       else if (val==valin) then
          noEncontrado = .false.
       else
          pos = MOD(pos+incr,ht%sizet) + 1
       endif
    enddo
    if(pos == -1 ) then
       nullify(local_id)
    else
       kk = memory_size(ht % list(pos) % l)       
       call memory_alloca(memor_loc,'HT % LIST','htalid',local_id,kk)
       local_id = ht % list(pos) % l
    endif
    
    if( present(memor_opt) ) then
       memor_opt = memor_loc
    else
       memor_ht  = memor_loc
    end if

  end function htalid_list
  
  function htalid(ht, valin) result(local_id)

    integer(ip)                 :: local_id
    type(hash_t), intent(in)    :: ht
    integer(ip),  intent(in)    :: valin 

    integer(ip)                 :: pos, incr, val
    logical(lg)                 :: noEncontrado

    if( .not. ht % lidson ) call runend("htalid: generate hast_t with lidson=true")
    pos   = MOD(valin,ht%sizet) + 1_ip
    incr  = MOD(valin,ht%sizet  - 1_ip)

    noEncontrado = .true.
    do while(noEncontrado)
       val = ht%dades(pos)
       if (val==NIL) then
          pos = -1_ip
          noEncontrado = .false.
       else if (val==valin) then
          noEncontrado = .false.
       else
          pos = MOD(pos+incr,ht%sizet) + 1
       endif
    enddo
    if(pos == -1 ) then
       local_id = 0_ip
    else
       local_id = ht % lids(pos)
    endif
    
  end function htalid

  !-----------------------------------------------------------------------
  !
  !> @brief   HTablePrimeNumber
  !> @details Return first prime number of the array HTablePrimes higher 
  !>          than sizet, if all elements are lower than sizet,  runends
  !
  !-----------------------------------------------------------------------

  function HTablePrimeNumber( sizet, AUTOMATIC_SIZE , ht )

    implicit none
    integer(ip)                           :: HTablePrimeNumber
    integer(ip),  intent(in)              :: sizet
    logical(lg),  intent(in),    optional :: AUTOMATIC_SIZE
    type(hash_t), intent(inout), optional :: ht
    integer(ip)                           :: ii,iimax
    logical(lg)                           :: if_automatic_size
    integer(ip)                           :: sizet_try
    
    if_automatic_size = .false.
    if( present(AUTOMATIC_SIZE) ) if_automatic_size = AUTOMATIC_SIZE
    
    ii    = size(HTablePrimes,KIND=ip) ! Last prime number position
    iimax = HTablePrimes(ii)           ! Max prime number
    if( sizet > iimax ) call runend("Not possible to use hash table with this sizet")
    
    if( if_automatic_size ) then
       !
       ! Automatic size
       !
       sizet_try = sizet * size_factor
       if( present(ht) ) ht % sizei = sizet_try
       if( sizet_try > iimax ) then
          ii = size(HTablePrimes,KIND=ip)
       else
          ii = 1
          do while(HTablePrimes(ii)<sizet_try)
             ii = ii + 1
          end do
       end if

    else
       !
       ! Look for nearest greater prime number
       !
       if( present(ht) ) ht % sizei = sizet
       ii = 1
       do while(HTablePrimes(ii)<sizet)
          ii = ii + 1
       end do
       
    end if
    
    HTablePrimeNumber = HTablePrimes(ii)

  end function HTablePrimeNumber

  !-----------------------------------------------------------------------
  !
  !> @brief   HTableMaxPrimeNumber
  !> @details Max prime number of the table
  !
  !-----------------------------------------------------------------------

  function HTableMaxPrimeNumber()

    implicit none
    integer(ip) :: HTableMaxPrimeNumber
    integer(ip) :: ii

    ii = size(HTablePrimes,KIND=ip)
    HTableMaxPrimeNumber= HTablePrimes(ii)

    return

  end function HTableMaxPrimeNumber

end module mod_htable
!> @}
