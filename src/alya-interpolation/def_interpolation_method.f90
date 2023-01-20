!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Maths
!> @{
!> @file    def_search_method.f90
!> @author  houzeaux
!> @date    2020-10-02
!> @brief   Abstract class for search methods
!> @details All search methods are extensions of search_method
!>
!>          The different classes:
!>          ----------------------
!>                             search_method
!>                                  ||
!>                           search_method_par
!>                  ||            ||            ||            ||
!>               maths_bin   maths_octree   maths_octbin   typ_kdrtee
!>
!>          Deferred procedures:
!>          --------------------
!>          init() ........ initialize the search strategy
!>          deallo() ...... deallocate the search strategy structure
!>          candidate() ... give a list of entity candidate when doing
!>                          a point search
!>          fill() ........ Fill int eh search structure?
!>
!-----------------------------------------------------------------------

module def_interpolation_method

  use def_kintyp_basic,                   only : ip,rp,lg,i1p,i2p,r1p,r2p,r3p
  use def_kintyp_mesh_basic,              only : mesh_type_basic
  use def_kintyp_mesh_basic,              only : bmsh_type_basic
  use def_kintyp_comm,                    only : comm_data_par_basic
  use def_master,                         only : IPARALL,kfl_paral
  use def_mat_coo,                        only : mat_coo
  use def_search_parall,                  only : search_parall
  use def_search_method,                  only : search_method
  use def_search_method,                  only : CANDIDATE_INSIDE
  use def_search_method,                  only : CANDIDATE_NEAREST
  use def_search_method,                  only : SEARCH_POINTS
  use def_search_method,                  only : SEARCH_BOUNDING_BOXES
  use mod_memory_basic,                   only : memory_alloca
  use mod_memory_basic,                   only : memory_deallo
  use mod_memory_basic,                   only : memory_size
  use mod_memory_basic,                   only : memory_resize
  use mod_memory_basic,                   only : memory_copy
  use mod_memory_tools,                   only : memory_counter_ini
  use mod_memory_tools,                   only : memory_counter_end
  use mod_optional_argument,              only : optional_argument
  use mod_communications_global,          only : PAR_ALLGATHER
  use mod_communications_point_to_point,  only : PAR_SEND_RECEIVE
  use mod_communications_global,          only : PAR_SUM
  use mod_communications_global,          only : PAR_MIN
  use mod_communications_global,          only : PAR_MAX
  use mod_communications_global,          only : PAR_AVERAGE
  use mod_maths_arrays,                   only : maths_maxloc_nonzero
  use mod_elmgeo,                         only : element_type
  use mod_elmgeo,                         only : elmgeo_natural_coordinates
  use mod_elmgeo,                         only : elmgeo_cartesian_derivatives
  use mod_elmgeo,                         only : elmgeo_natural_coordinates_on_boundaries
  use mod_elmgeo,                         only : elmgeo_projection_on_a_face
  use mod_elmgeo,                         only : elmgeo_jacobian_boundary
  use mod_htable,                         only : hash_t
  use mod_htable,                         only : htaini
  use mod_htable,                         only : htaadd
  use mod_htable,                         only : htades
  use mod_htable,                         only : htalid
  use def_mpi
  use mod_std
#include "def_mpi.inc"
  
  implicit none
  private

  real(rp),      parameter :: epsil = epsilon(1.0_rp)
  character(24), parameter :: vacal = 'def_interpolation_method'
  
  ! Mode
  integer(ip),   parameter :: INT_SEQUENTIAL             =  0
  integer(ip),   parameter :: INT_PARALLEL               =  1
  ! Search method
  integer(ip),   parameter :: INT_BIN                    =  0
  integer(ip),   parameter :: INT_OCTREE                 =  1
  integer(ip),   parameter :: INT_OCTBIN                 =  2
  integer(ip),   parameter :: INT_KDTREE                 =  3
  integer(ip),   parameter :: INT_SKDTREE                =  4
  ! Search method name
  ! Interpolation method
  integer(ip),   parameter :: INT_OFF                    = -1
  integer(ip),   parameter :: INT_ELEMENT_INTERPOLATION  =  0
  integer(ip),   parameter :: INT_BOUNDARY_INTERPOLATION =  1
  integer(ip),   parameter :: INT_NEAREST_NODE           =  2
  integer(ip),   parameter :: INT_GLOBAL_NUMBERING       =  3
  integer(ip),   parameter :: INT_ELEMENT_VALUE          =  4
  integer(ip),   parameter :: INT_BOUNDARY_VALUE         =  5

  type           :: interpolation_data
     integer(ip) :: interpolation_method                                ! Element interpolation=1/boundary interpolation=2/nearset boundary=3
     integer(ip) :: search_method_seq                                   ! bin=1/octree=2/octbin=3/kdtree=4
     integer(ip) :: search_method_par                                   ! bin=1/octree=2/octbin=3/kdtree=4
     real(rp)    :: toler_rel                                           ! Relative tolerance
     real(rp)    :: toler_abs                                           ! Absolute tolerance
     integer(ip) :: mode                                                ! Execution mode
     logical(lg) :: deriv                                               ! If derivative should be computed
     logical(lg) :: myself                                              ! If myself (rank) should be included 
     logical(lg) :: force_find                                          ! When point is knonw to be in mesh, force a solution (lowering tol if necessary)
     integer(ip) :: max_it_force                                        ! Orders of magnitude the tolerance can be lowered to force_find
  end type interpolation_data

  type :: interpolation
     character(LEN=:),  allocatable     :: name                           ! Name
     type(interpolation_data)           :: input_data                     ! Input data
     class(search_method),      pointer :: search_method_seq              ! Sequential search
     type(search_parall)                :: parallel_search                ! Parallel search
     type(mat_coo),             pointer :: matrix(:)                      ! Interpolation matrix
     type(mat_coo),             pointer :: matder(:)                      ! Interpolation matrix for derivatives
     logical(lg),               pointer :: found(:)                       ! It point has been found
     real(rp),                  pointer :: dista(:)
     integer(ip),               pointer :: perm(:)                        ! Permutation on target
     integer(ip)                        :: nn                             ! Number of points 
     integer(ip)                        :: nd                             ! Dimension
     integer(ip)                        :: nrank                          ! Number of ranks involved
     integer(8)                         :: memor(2)                       ! Memory counter
     real(rp)                           :: times(10)                      ! Timings
     real(rp)                           :: stats(10)                      ! Statistics
   contains
     procedure,               pass      :: init
     procedure,               pass      :: input
     procedure,               pass      :: deallo
     procedure,               pass      :: check                          ! Check data
     procedure,               pass      :: sparse_matrix                  ! Compute sparse interpolaiton matrices
     procedure,               pass      :: sparse_matrix_permutation      ! Permute columns
     procedure,               pass      :: interpolation_entity           ! Entity interpolation
     procedure,               pass      :: interpolation_nearest_node     ! Nearest node interpolation
     procedure,               pass      :: interpolation_global_numbering ! Global numbering
     procedure,               pass      :: lost                           ! Number of lost points
     procedure,               pass      :: merge                          ! Merge two interpolations
     procedure,               pass      :: output                         ! Output statistics
     procedure,               pass      :: permutation                    ! Permutation
     procedure,               pass      :: preprocess
     procedure,               pass      :: preprocess_single
     procedure,               pass      :: values_10
     procedure,               pass      :: values_11
     procedure,               pass      :: values_22
     procedure,               pass      :: values_10_i
     procedure,               pass      :: values_11_i
     procedure,               pass      :: values_22_i
     procedure,               pass      :: derivatives_12
     procedure,               pass      :: distances_21
     generic                            :: values =>    &
          &                                values_10,   &
          &                                values_11,   &
          &                                values_22,   &
          &                                values_10_i, &
          &                                values_11_i, &
          &                                values_22_i
     generic                            :: derivatives => &
          &                                derivatives_12
     generic                            :: distances =>   &
          &                                distances_21
     
  end type interpolation
  !
  ! Interface for CANDIDATE_NEAREST, CANDIDATE_INSIDE, CANDIDATE_INSIDE interpolations
  !
  abstract interface
     subroutine interpolation_generic(self,search,xx,mesh,lelem,shapf,deriv,dista,lenty,mask,TOLER,MEMORY_COUNTER)
       import                                                 :: interpolation
       import                                                 :: search_method
       import                                                 :: mesh_type_basic
       import                                                 :: rp
       import                                                 :: ip
       import                                                 :: lg
       class(interpolation),                    intent(inout) :: self
       class(search_method),                    intent(inout) :: search
       real(rp),                       pointer, intent(in)    :: xx(:,:)
       class(mesh_type_basic),                  intent(in)    :: mesh
       integer(ip),                    pointer, intent(inout) :: lelem(:)
       real(rp),                       pointer, intent(inout) :: shapf(:,:)
       real(rp),         optional,     pointer, intent(inout) :: deriv(:,:,:)
       real(rp),         optional,     pointer, intent(inout) :: dista(:)
       integer(ip),      optional,     pointer, intent(inout) :: lenty(:,:)
       logical(lg),      optional,     pointer, intent(in)    :: mask(:)
       real(rp),         optional,              intent(in)    :: TOLER
       integer(8),       optional,              intent(inout) :: MEMORY_COUNTER(2) !< Memory counters       
     end subroutine interpolation_generic
  end interface 

  interface
     subroutine cputim(rtime)
       import :: rp
       implicit none
       real(rp), intent(out) :: rtime
     end subroutine cputim
     subroutine runend(message) 
       implicit none
       character(*)          :: message
     end subroutine runend
  end interface

  public :: interpolation

  public :: INT_SEQUENTIAL
  public :: INT_PARALLEL
  public :: INT_BIN                   
  public :: INT_OCTREE                
  public :: INT_OCTBIN                
  public :: INT_KDTREE                
  public :: INT_SKDTREE                
  public :: INT_ELEMENT_INTERPOLATION                       
  public :: INT_ELEMENT_VALUE                       
  public :: INT_BOUNDARY_VALUE                       
  public :: INT_BOUNDARY_INTERPOLATION
  public :: INT_NEAREST_NODE 
  public :: INT_GLOBAL_NUMBERING

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-16
  !> @brief   Initialize
  !> @details Initialize
  !> 
  !-----------------------------------------------------------------------

  subroutine init(self)

    class(interpolation), intent(inout) :: self

    nullify(self % search_method_seq)
    nullify(self % matrix)
    nullify(self % matder)
    nullify(self % found)
    nullify(self % dista)
    nullify(self % perm)
    call self % parallel_search % init()

    self % input_data % interpolation_method = INT_ELEMENT_INTERPOLATION
    self % input_data % search_method_seq    = -1
    self % input_data % search_method_par    = -1
    self % input_data % toler_rel            = epsil
    self % input_data % toler_abs            = 0.0_rp
    self % input_data % mode                 = INT_SEQUENTIAL
    self % input_data % deriv                = .false.
    self % input_data % myself               = .true.
    self % input_data % force_find           = .false.
    self % input_data % max_it_force         = 4

    self % nn                                = 0
    self % nd                                = 0
    self % nrank                             = 0    
    self % memor                             = 0_8
    self % times                             = 0.0_rp
    self % stats                             = 0.0_rp

  end subroutine init

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-16
  !> @brief   Associate a permutation
  !> @details Associate a permutation
  !> 
  !-----------------------------------------------------------------------

  subroutine permutation(self,perm,MEMORY_COUNTER)

    class(interpolation),           intent(inout) :: self
    integer(ip),          pointer,  intent(inout) :: perm(:)
    integer(8),           optional, intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(8)                                    :: memor_loc(2)
    integer(ip)                                   :: ii,jj
    
    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)
    if( .not. associated(self % perm) ) then
       call memory_copy(memor_loc,'SELF % PERM',vacal,perm,self % perm,'DO_NOT_DEALLOCATE')
    else
       do ii = 1,memory_size(self % perm)
          jj = self % perm(ii)
          self % perm(ii) = perm(jj)
       end do
    end if
    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)

  end subroutine permutation
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-16
  !> @brief   Deallocate
  !> @details Deallocate.
  !>          IMPORTANT: Search method are NOT deallocated, only
  !>          pointers to them are nullified
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo(self,MEMORY_COUNTER)
 
    class(interpolation),           intent(inout) :: self
    integer(ip)                                   :: ii
    integer(8),           optional, intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(8)                                    :: memor_loc(2)

    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)
    !
    ! Matrices
    !
    if( associated(self % matrix) ) then
       do ii = int(lbound(self % matrix,1),ip),int(ubound(self % matrix,1),ip)
          call self % matrix(ii) % deallo(MEMORY_COUNTER=memor_loc)
       end do
       deallocate(self % matrix)
    end if
    if( associated(self % matder) ) then
       do ii = int(lbound(self % matder,1),ip),int(ubound(self % matder,1),ip)
          call self % matder(ii) % deallo(MEMORY_COUNTER=memor_loc)
       end do
       deallocate(self % matder)
    end if
    !
    ! Search methods
    !
    nullify(self % search_method_seq)
    if( self % input_data % mode == INT_PARALLEL ) then       
       call self % parallel_search % deallo(MEMORY_COUNTER=memor_loc)
       nullify(self % parallel_search % search_method_par)
    end if
    !
    ! Others
    !
    if( allocated(self % name) ) deallocate(self % name)    
    call memory_deallo(memor_loc,'SELF % FOUND',vacal,self % found)
    call memory_deallo(memor_loc,'SELF % DISTA',vacal,self % dista)
    call memory_deallo(memor_loc,'SELF % PERM' ,vacal,self % perm)

    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)

  end subroutine deallo

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2021-05-12
  !> @brief   Number of lost points
  !> @details Number of lost points
  !> 
  !-----------------------------------------------------------------------

  pure function lost(self) result(num)

    class(interpolation), intent(in) :: self
    integer(ip)                      :: num
    
    if( associated(self % found) ) then
       num = self % nn - int(count(self % found),ip)
    else
       num = 0
    end if

  end function lost
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-16
  !> @brief   Input
  !> @details Input
  !> 
  !-----------------------------------------------------------------------

  subroutine input(self,            &
       search_method_seq,           &
       search_method_par,           &
       COMM,                        &
       MODE,                        &
       INTERPOLATION_METHOD,        &
       NAME,                        &
       DERIVATIVES,                 &
       MYSELF,                      &
       FORCE_FIND,                  &
       MAX_IT_FORCE                 )

    class(interpolation),                   intent(inout) :: self
    class(search_method), optional, target, intent(in)    :: search_method_seq
    class(search_method), optional, target, intent(in)    :: search_method_par
    MY_MPI_COMM,          optional,         intent(in)    :: COMM
    integer(ip),          optional,         intent(in)    :: MODE
    integer(ip),          optional,         intent(in)    :: INTERPOLATION_METHOD
    character(LEN=*),     optional,         intent(in)    :: NAME
    logical(lg),          optional,         intent(in)    :: DERIVATIVES
    logical(lg),          optional,         intent(in)    :: MYSELF
    logical(lg),          optional,         intent(in)    :: FORCE_FIND
    integer(ip),          optional,         intent(in)    :: MAX_IT_FORCE
    !
    ! Sequential and parallel search methods
    !
    if( present(search_method_seq) ) then
       self % search_method_seq => search_method_seq
       self % input_data % mode =  INT_SEQUENTIAL
    end if

    if( present(search_method_par) .and. IPARALL ) then
       self % parallel_search % search_method_par => search_method_par
       self % input_data % mode                   =  INT_PARALLEL
       if( .not. present(COMM) ) then
          call runend('DEF_INTERPOLATION_METHOD: FOR PARALLEL SEARCH, COMMUNICATOR IS COMPULSORY')
       end if
    end if
    !
    ! Interpolation method
    !
    if( present(INTERPOLATION_METHOD) ) then
       self % input_data % interpolation_method = INTERPOLATION_METHOD
    end if
    !
    ! Communicator
    !
    if( present(MODE) ) self % input_data % mode = MODE
    call self % parallel_search % input(COMM)
    !
    ! Number of ranks
    !
    if( self % input_data % mode == INT_PARALLEL ) then
       self % nrank = int(self % parallel_search % comm % size4,ip)
    else
       self % nrank = 1
       self % parallel_search % comm % RANK4 = 0_4 
    end if
    !
    ! Name
    !
    self % name = optional_argument('My interpolation',NAME)
    !
    ! Derivatives and myself
    !
    self % input_data % deriv  = optional_argument(.false.,DERIVATIVES)
    self % input_data % myself = optional_argument(.true. ,MYSELF)
    !
    ! FORCE_FIND
    !
    self % input_data % force_find   = optional_argument(.false.,FORCE_FIND  )
    self % input_data % max_it_force = optional_argument(4_ip   ,MAX_IT_FORCE)
    !
    ! Method for candidate
    !
!!$    if( present(CANDIDATE_SEQUENTIAL) ) then
!!$       
!!$    else
!!$       if(      self % input_data % interpolation_method == INT_ELEMENT_INTERPOLATION ) then
!!$          self % input_data % candidate_seq = CANDIDATE_INSIDE
!!$          self % input_data % candidate_par = CANDIDATE_INSIDE
!!$       else if( self % input_data % interpolation_method == INT_ELEMENT_INTERPOLATION ) then
!!$          self % input_data % candidate_seq = CANDIDATE_INSIDE
!!$          self % input_data % candidate_par = CANDIDATE_INSIDE
!!$       end if
!!$    end if
    !
    ! Check data
    !
    call self % check()
    
  end subroutine input
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-16
  !> @brief   Check input data
  !> @details Check input data
  !> 
  !-----------------------------------------------------------------------

  subroutine check(self)
    
    class(interpolation), intent(inout) :: self

    if(      self % input_data % interpolation_method == INT_NEAREST_NODE ) then
       if( self % search_method_seq % fill_method /= SEARCH_POINTS ) &
            call runend('NEAREST NODE NEEDS SEARCH METHOD TO BE FILLED WITH POINTS')
    else if( self % input_data % interpolation_method == int_ELEMENT_INTERPOLATION ) then
       if( self % search_method_seq % fill_method /= SEARCH_BOUNDING_BOXES ) &
            call runend('ELEMENT INTERPOLATION SEARCH METHOD TO BE FILLED WITH BOUNDING BOXES')
    else if( self % input_data % interpolation_method == int_BOUNDARY_INTERPOLATION ) then
       if( self % search_method_seq % fill_method /= SEARCH_BOUNDING_BOXES ) &
            call runend('BOUNDARY INTERPOLATION SEARCH METHOD TO BE FILLED WITH BOUNDING BOXES')
    end if

  end subroutine check
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-16
  !> @brief   Output stats
  !> @details Output stats
  !> 
  !-----------------------------------------------------------------------

  subroutine output(self,UNIT,FILENAME,HEADER)

    class(interpolation),           intent(in)    :: self
    integer(ip),          optional, intent(in)    :: UNIT
    character(LEN=*),     optional, intent(in)    :: FILENAME
    character(LEN=*),     optional, intent(in)    :: HEADER
    integer(4)                                    :: unit4
    integer(ip),          parameter               :: rn = 7
    real(rp)                                      :: rmin(rn),rmax(rn),rave(rn)
    real(rp)                                      :: rlost,rself,rtest,rnn
    character(len=:),     allocatable             :: my_header
    logical(lg)                                   :: if_write
    character(len=:),     allocatable             :: seq_search
    character(len=:),     allocatable             :: par_search
    
    unit4     = int(optional_argument(6_ip,UNIT),4)
    my_header = optional_argument('DUMMY',HEADER)
    
    rmin = 0.0_rp ; rmax = 0.0_rp ; rave = 0.0_rp
    !
    ! Integer stats
    !
    rlost = self % stats(1)
    rtest = self % stats(2)
    rself = self % stats(3)
    rnn   = real(self % nn,rp)
    rmin  = (/ rnn,rnn-rlost,rlost,rself,rtest,self % times(1),self % times(2)/)
    rmax  = rmin
    rave  = rmin
    call PAR_MIN    (rn,rmin)
    call PAR_MAX    (rn,rmax)
    call PAR_AVERAGE(rn,rave)
    !
    ! Who writes
    !
    if_write = .false.
    if( associated(self % parallel_search % search_method_par) ) then
       if( self % parallel_search % comm % RANK4 == 0 ) if_write = .true.
    end if
    if( associated(self % parallel_search % search_method_par) ) then
       if( allocated(self % parallel_search % search_method_par % name) ) then
          par_search = self % parallel_search % search_method_par % name
       else
          seq_search = 'PARALLEL SEARCH'
       end if
    else
       par_search = 'NONE'
    end if
    if( associated(self % search_method_seq) ) then
       if( allocated(self % search_method_seq % name) ) then
          seq_search = self % search_method_seq % name
       else
          seq_search = 'SEQUENTIAL SEARCH'
       end if
    else
       seq_search = 'NONE'
    end if
    !
    ! Output
    !
    if( if_write ) then
       write(unit4,'(a)') ' '
       write(unit4, 99) my_header , 'Interpolation results for: '//trim(self % name)
       write(unit4, 99) my_header , '-------------------------'
       write(unit4, 99) my_header , 'Sequential search method: '//trim(seq_search)
       write(unit4, 99) my_header , 'Parallel search method:   '//trim(par_search)
       write(unit4, 99) my_header , 'Stats (min, max, average):'
       write(unit4,  1) my_header , int(rmin(1),ip),int(rmax(1),ip),int(rave(1),ip)
       write(unit4,  2) my_header , int(rmin(2),ip),int(rmax(2),ip),int(rave(2),ip)
       write(unit4,  3) my_header , int(rmin(3),ip),int(rmax(3),ip),int(rave(3),ip)
       write(unit4,  4) my_header , int(rmin(4),ip),int(rmax(4),ip),int(rave(4),ip)
       write(unit4,  5) my_header , int(rmin(5),ip),int(rmax(5),ip),int(rave(5),ip)
       write(unit4,  6) my_header , rmin(6),rmax(6),rave(6)
       write(unit4,  7) my_header , rmin(7),rmax(7),rave(7)
       write(unit4,'(a)') ' '
    end if

    if( allocated(my_header)  ) deallocate(my_header)
    if( allocated(seq_search) ) deallocate(seq_search)
    if( allocated(par_search) ) deallocate(par_search)
    
99  format(a,a)
1   format(a,'# of searches=         ' , 1x , i11 , 1x , i11 , 1x , i11 )
2   format(a,'# successful searches= ' , 1x , i11 , 1x , i11 , 1x , i11 )
3   format(a,'# failed searches=     ' , 1x , i11 , 1x , i11 , 1x , i11 )
4   format(a,'# Myself first=        ' , 1x , i11 , 1x , i11 , 1x , i11 )
5   format(a,'# tests per entity=    ' , 1x , i11 , 1x , i11 , 1x , i11 )
6   format(a,'Time search=           ' , f10.3 , ' s', f10.3 , ' s' , f10.3, ' s' )
7   format(a,'Time candidate=        ' , f10.3 , ' s', f10.3 , ' s' , f10.3, ' s' )
    
  end subroutine output
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-16
  !> @brief   Initialize
  !> @details Initialize
  !> 
  !-----------------------------------------------------------------------

  subroutine preprocess(self1,xx,mesh,ll,mask,MEMORY_COUNTER,SECOND_TRY)

    class(interpolation),                    intent(inout) :: self1
    real(rp),                       pointer, intent(in)    :: xx(:,:)
    class(mesh_type_basic),                  intent(in)    :: mesh
    integer(ip),          optional, pointer, intent(in)    :: ll(:)
    logical(lg),          optional, pointer, intent(in)    :: mask(:)
    integer(8),           optional,          intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    class(interpolation), optional,          intent(inout) :: SECOND_TRY

    integer(ip)                                            :: ii,kk
    integer(ip),                    pointer                :: permu(:)
    integer(ip),                    pointer                :: invpr(:)
    real(rp),                       pointer                :: x2(:,:)
    integer(8)                                             :: memor_loc(2)
    integer(8)                                             :: memor_dum(2)

    call memory_counter_ini(memor_loc,self1 % memor,MEMORY_COUNTER)
    memor_dum = 0_8
    !
    ! Search with first try
    !
    call self1 % preprocess_single(xx,mesh,ll,mask,MEMORY_COUNTER=memor_loc) 

    if( present(SECOND_TRY) ) then
       nullify(x2,permu,invpr)
       kk = self1 % lost()
       call memory_alloca(memor_dum,'X2'   ,vacal,x2   ,self1 % nd,kk)
       call memory_alloca(memor_dum,'PERMU',vacal,permu,kk)
       kk = 0
       do ii = 1,self1 % nn
          if( .not. self1 % found(ii) ) then
             kk        = kk + 1
             permu(kk) = ii
             x2(:,kk)  = xx(:,ii)
          end if
       end do
       !
       ! Search with second try
       !
       call SECOND_TRY % preprocess_single(x2,mesh,ll,mask)       
       !
       ! Inverse permutation, only useful in sequential mode
       !
       if( self1 % input_data % mode ==  INT_SEQUENTIAL ) then
          call memory_alloca(memor_dum,'INVPR',vacal,invpr,self1 % nn)
          kk = 0
          do ii = 1,self1 % nn
             if( self1 % found(ii) ) then
                kk        = kk + 1
                invpr(kk) = ii
             end if
          end do
          do ii = 1,self1 % nn
             if( .not. self1 % found(ii) ) then
                kk        = kk + 1
                invpr(kk) = ii
             end if
          end do
       end if
       !
       ! Merge second try into first try
       !
       call self1 % merge(SECOND_TRY,permu,invpr,MEMORY_COUNTER=memor_loc)
       
       call memory_deallo(memor_dum,'X2'   ,vacal,x2)
       call memory_deallo(memor_dum,'PERMU',vacal,permu)
       call memory_deallo(memor_dum,'INVPR',vacal,invpr)
       
    end if
    
    call memory_counter_end(memor_loc,self1 % memor,MEMORY_COUNTER)
 
  end subroutine preprocess
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-16
  !> @brief   Initialize
  !> @details Initialize
  !> 
  !-----------------------------------------------------------------------

  subroutine preprocess_single(self,xx,mesh,ll,mask,MEMORY_COUNTER)

    class(interpolation),                    intent(inout) :: self
    real(rp),                       pointer, intent(in)    :: xx(:,:)
    class(mesh_type_basic),                  intent(in)    :: mesh
    integer(ip),          optional, pointer, intent(in)    :: ll(:)
    logical(lg),          optional, pointer, intent(in)    :: mask(:)
    integer(8),           optional,          intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(ip)                                            :: irank,nrank,nn,nd
    type(r2p),                      pointer                :: shapf(:)
    type(r3p),                      pointer                :: deriv(:)
    type(i1p),                      pointer                :: lelem(:)
    type(i2p),                      pointer                :: lenty(:)
    logical(lg),                    pointer                :: llost(:)
    real(rp)                                               :: toler
    integer(ip)                                            :: ii,nlost
    integer(8)                                             :: memor_loc(2)
   
    !--------------------------------------------------------------------
    !
    ! Errors 
    !
    !--------------------------------------------------------------------
    
    if( .not. associated(self % search_method_seq) ) &
         call runend('DEF_INTERPOLATION_METHOD: SEQUENTIAL SEARCH METHOD NOT DEFINED')
    if( self % input_data % mode == INT_PARALLEL ) then
       if( .not. associated(self % parallel_search % search_method_par) ) &
            call runend('DEF_INTERPOLATION_METHOD: PARALLEL SEARCH METHOD NOT DEFINED')
    end if

    !if( self % search_method_seq % kfl_exists == 0 ) &
    !     call runend('DEF_INTERPOLATION_METHOD: SEARCH METHOD HAS NOT BEEN FILLED')         
    !if( self % input_data % mode == INT_PARALLEL ) then
    !   if( self % parallel_search % search_method_par % kfl_exists == 0 ) &
    !     call runend('DEF_INTERPOLATION_METHOD: PARALLEL SEARCH METHOD HAS NOT BEEN FILLED')         
    !end if
       
    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)

    nullify(shapf)
    nullify(deriv)
    nullify(lelem)
    nullify(lenty)
    nullify(llost)
    
    toler                       = self % input_data % toler_rel
    nn                          = memory_size(xx,2_ip)
    nd                          = max(mesh % ndime,memory_size(xx,1_ip))
    self % nn                   = nn
    self % nd                   = nd
    nrank                       = self % nrank    
    self % parallel_search % nd = nd

    if( nd /= memory_size(xx,1_ip) .and. memory_size(xx,1_ip) > 0 ) &
         call runend('DEF_INTERPOLATION_METHOD: INCOMPATIBLE DIMENSIONS')

    !--------------------------------------------------------------------
    !
    ! Allocate memory and initialize
    ! By default, points are lost => LLOST(II) = .true.
    !
    !-------------------------------------------------------------------- 
    
    call memory_alloca(memor_loc,'SHAPF'       ,vacal,shapf,nrank,LBOUN=0_ip)
    call memory_alloca(memor_loc,'DERIV'       ,vacal,deriv,nrank,LBOUN=0_ip)
    call memory_alloca(memor_loc,'LELEM'       ,vacal,lelem,nrank,LBOUN=0_ip)
    call memory_alloca(memor_loc,'LENTY'       ,vacal,lenty,nrank,LBOUN=0_ip)
    call memory_alloca(memor_loc,'LLOST'       ,vacal,llost,nn)
    call memory_alloca(memor_loc,'SELF % FOUND',vacal,self % dista,nn,INIT_VALUE=-huge(1.0_rp))

    allocate(self % matrix(0:nrank-1))
    do irank = 0,nrank-1
       call self % matrix(irank) % init()
    end do
    if( self % input_data % deriv ) then
       allocate(self % matder(0:nrank-1))
       do irank = 0,nrank-1
          call self % matder(irank) % init()
       end do
    end if
    do ii = 1,nn
       llost(ii) = .true.
    end do
    
    !--------------------------------------------------------------------
    !
    ! Interpolate
    !
    !--------------------------------------------------------------------

   select case ( self % input_data % interpolation_method )
       
    case ( INT_ELEMENT_INTERPOLATION , INT_ELEMENT_VALUE , INT_BOUNDARY_INTERPOLATION , INT_BOUNDARY_VALUE )
       !
       ! Element and boundary interpolations, and element value
       !
       call self % interpolation_entity(xx,mesh,shapf,deriv,lelem,lenty,llost,mask,MEMORY_COUNTER=memor_loc)
       
    case ( INT_NEAREST_NODE )
       !
       ! Nearest node
       !
       call self % interpolation_nearest_node(xx,mesh,shapf,deriv,lelem,lenty,llost,mask,MEMORY_COUNTER=memor_loc)
       
    case ( INT_GLOBAL_NUMBERING )
       !
       ! Global numbering
       !
       if( present(ll) ) then
          call self % interpolation_global_numbering(xx,ll,mesh,shapf,deriv,lelem,lenty,llost,mask,MEMORY_COUNTER=memor_loc)
       else
          call runend('PREPROCESS_SINGLE: GLOBAL NUMBERING ARRAY REQUIRED')
       end if
       
    case default
       !
       ! Nothing
       !
       call runend('DEF_INTERPOLATION_METHODS: DO NOT KNOW WHAT TO DO')
       
    end select
   
    !--------------------------------------------------------------------
    !
    ! Allocate and construct sparse matrix
    !
    !--------------------------------------------------------------------

    
    call self % sparse_matrix(lelem,lenty,shapf,deriv,MEMORY_COUNTER=memor_loc)

    call memory_alloca(memor_loc,'SELF % FOUND',vacal,self % found,nn)
    do ii = 1,nn
       self % found(ii) = .not. llost(ii)
    end do
        
    !--------------------------------------------------------------------
    !
    ! Statistics
    !
    !--------------------------------------------------------------------
    
    nlost = 0
    if( nn > 0 ) nlost = int(count(llost),ip)
    self % stats(1) = real(nlost,rp)

    !--------------------------------------------------------------------
    !
    ! Deallocate 
    !
    !--------------------------------------------------------------------

    call memory_deallo(memor_loc,'LLOST',vacal,llost)
    call memory_deallo(memor_loc,'LENTY',vacal,lenty)
    call memory_deallo(memor_loc,'LELEM',vacal,lelem)
    call memory_deallo(memor_loc,'DERIV',vacal,deriv)
    call memory_deallo(memor_loc,'SHAPF',vacal,shapf)
    
    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)
    
   end subroutine preprocess_single

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-16
  !> @brief   Merge
  !> @details Merge two interpolations:
  !>          self3 = self1 U self2
  !> 
  !-----------------------------------------------------------------------

  subroutine merge(self1,self2,permu,invpr,MEMORY_COUNTER)

    class(interpolation),                     intent(inout) :: self1
    class(interpolation),                     intent(inout) :: self2
    integer(ip),                    pointer,  intent(in)    :: permu(:)
    integer(ip),                    pointer,  intent(in)    :: invpr(:)
    integer(8),           optional,           intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(ip)                                             :: ii,jj
    integer(ip)                                             :: irank
    integer(8)                                              :: memor_loc(2)
    
    call memory_counter_ini(memor_loc,self1 % memor,MEMORY_COUNTER)
    !
    ! Matrix. In parallel, permutation is carried out by lrecv_perm. In sequential,
    ! the matrix should be reordered explicity
    !
    if( self1 % input_data % mode == INT_PARALLEL ) then       
       do irank = 0,self1 % nrank-1
          call self1 % matrix(irank) % merge(self2 %  matrix(irank),MEMORY_COUNTER=memor_loc)          
       end do
    else
       do irank = 0,self1 % nrank-1
          call self1 % matrix(irank) % merge(self2 %  matrix(irank),MEMORY_COUNTER=memor_loc,PERMU=invpr)          
       end do       
    end if
    if( self1 % input_data % deriv ) then 
       if( self1 % input_data % mode == INT_PARALLEL ) then       
          do irank = 0,self1 % nrank-1
             call self1 % matder(irank) % merge(self2 %  matder(irank),MEMORY_COUNTER=memor_loc)          
          end do
       else
          do irank = 0,self1 % nrank-1
             call self1 % matder(irank) % merge(self2 %  matder(irank),MEMORY_COUNTER=memor_loc,PERMU=invpr)          
          end do
       end if
    end if
    !
    ! List of found
    !
    do ii = 1,self2 % nn
       jj = permu(ii)
       self1 % found(jj) = self2 % found(ii)
    end do
    !
    ! Parallel data structure
    !    
    if( self1 % input_data % mode == INT_PARALLEL ) then       
       do jj = 1,self2 % parallel_search % comm % lrecv_dim
          ii = self2 % parallel_search % comm % lrecv_perm(jj)
          self2 % parallel_search % comm % lrecv_perm(jj) = permu(ii)
       end do  
       call self1 % parallel_search % comm % merge(self2 % parallel_search % comm,MEMORY_COUNTER=memor_loc)
    end if
    
    call memory_counter_end(memor_loc,self1 % memor,MEMORY_COUNTER)    
    
  end subroutine merge
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-29
  !> @brief   Nearest node
  !> @details Nearest node interpolation
  !> 
  !-----------------------------------------------------------------------

  subroutine interpolation_nearest_node(&
    self,xx,mesh,shapf,deriv,lelem,lenty,llost,mask,MEMORY_COUNTER,search_seq_in,search_par_in)

    class(interpolation),                    intent(inout) :: self
    real(rp),                       pointer, intent(in)    :: xx(:,:)
    class(mesh_type_basic),                  intent(in)    :: mesh
    type(r2p),                      pointer, intent(inout) :: shapf(:)
    type(r3p),                      pointer, intent(inout) :: deriv(:)
    type(i1p),                      pointer, intent(inout) :: lelem(:)
    type(i2p),                      pointer, intent(inout) :: lenty(:)
    logical(lg),                    pointer, intent(inout) :: llost(:)
    logical(lg),          optional, pointer, intent(in)    :: mask(:)
    integer(8),           optional,          intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    class(search_method), optional, target,  intent(inout) :: search_seq_in     ! Sequential search
    class(search_method), optional, target,  intent(inout) :: search_par_in     ! Sequential search   
    integer(ip)                                            :: irank,nrank,nn
    integer(ip)                                            :: ii,kk,ineig,ipoin
    integer(ip)                                            :: my_rank
    real(rp)                                               :: toler
    integer(8)                                             :: memor_loc(2)
    type(r1p),                      pointer                :: dista_recv(:)    
    type(r2p),                      pointer                :: xx_recv(:)
    integer(ip),                    pointer                :: nn_send(:)
    integer(ip),                    pointer                :: nn_recv(:)
    type(i1p),                      pointer                :: lista_send(:)  
    type(i1p),                      pointer                :: lista_recv(:)
    class(search_method),           pointer                :: search_seq        ! Sequential search
    class(search_method),           pointer                :: search_par        ! Sequential search   

    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)

    nullify(dista_recv)
    nullify(xx_recv)
    nullify(nn_send)
    nullify(nn_recv)
    nullify(lista_send)
    nullify(lista_recv)

    !--------------------------------------------------------------------
    !
    ! Search methods
    !
    !--------------------------------------------------------------------
    
    if( present(search_seq_in) ) then
       search_seq => search_seq_in
    else
       search_seq => self % search_method_seq
    end if
    
    if( present(search_par_in) ) then
       search_par => search_par_in
    else
       search_par => self % parallel_search % search_method_par
    end if
    
    !--------------------------------------------------------------------
    !
    ! Allocate
    !
    !--------------------------------------------------------------------

    nn      = memory_size(xx,2_ip)
    toler   = self % input_data % toler_rel
    nrank   = self % nrank
    my_rank = int(self % parallel_search % comm % RANK4,ip)

    call memory_alloca(memor_loc,'DISTA_RECV',vacal,dista_recv,nrank,LBOUN=0_ip)

    !--------------------------------------------------------------------
    !
    ! Search
    !
    !--------------------------------------------------------------------

    if( self % input_data % mode == INT_PARALLEL ) then

#ifdef __PGI
       call self % parallel_search % send_list_recv_points(&
            XX=xx,NN_SEND=nn_send,NN_RECV=nn_recv,LISTA_SEND=lista_send,XX_RECV=xx_recv,&
            METHOD=CANDIDATE_NEAREST,MASK=std_log_1,MEMORY_COUNTER=memor_loc)
#else
       call self % parallel_search % send_list_recv_points(&
            xx,nn_send,nn_recv,lista_send,xx_recv,&
            METHOD=CANDIDATE_NEAREST,MEMORY_COUNTER=memor_loc)
#endif
       
!!$       block
!!$         use def_master
!!$         integer(ip) :: irank
!!$         do irank = 0,nrank - 1
!!$            if( nn_send(irank) > 0 ) write(6,*)'send near=',kfl_paral,'->',irank,': ',nn_send(irank)
!!$            if( nn_recv(irank) > 0 ) write(6,*)'recv nrea=',kfl_paral,'<-',irank,': ',nn_recv(irank),memory_size(xx_recv(irank) % a)
!!$         end do
!!$         call runend('O.K.!')
!!$       end block
 
       do irank = 0,nrank-1
          if( nn_recv(irank) > 0 ) then
             call nearest_node(&
                  self,                     &
                  search_seq,               &
                  xx_recv(irank)%a,         &
                  mesh,                     &
                  lelem(irank)%l,           &
                  shapf(irank)%a,           &
                  deriv(irank)%a,           &
                  dista_recv(irank) % a,    &
                  lenty(irank)%l,           &
                  mask,                     &
                  toler,                    &
                  memor_loc                 )
          end if
       end do    

       call self % parallel_search % dista_comm(&
            nn,nn_send,nn_recv,lista_send,lista_recv,dista_recv,&
            MEMORY_COUNTER=memor_loc)

       do irank = 0,self % nrank-1
          do ii = 1,memory_size(lista_recv(irank) % l)
             if( lista_recv(irank) % l(ii) == 0 ) then
                lelem(irank) % l(ii) = 0
             end if
          end do
       end do
       kk = 0
       do ineig = 1,self % parallel_search % comm % nneig
          do ii = 1,self % parallel_search % comm % lrecv_size(ineig+1)-self % parallel_search % comm % lrecv_size(ineig)
             kk    = kk + 1
             ipoin = self % parallel_search % comm % lrecv_perm(kk)
             llost(ipoin) = .false.
          end do
       end do

    else

       call nearest_node(             &
            self,                     &
            search_seq,               &
            xx,                       &
            mesh,                     &
            lelem(my_rank)%l,         &
            shapf(my_rank)%a,         &
            deriv(my_rank)%a,         &
            dista_recv(my_rank) % a,  &
            lenty(my_rank)%l,         &
            mask,                     &
            toler,                    &
            memor_loc                 )

       do ii = 1,memory_size(lelem(my_rank) % l)
          if( lelem(my_rank) % l(ii) /= 0 ) llost(ii) = .false.
       end do

    end if
    
!!$    block
!!$      use def_master, only : kfl_paral
!!$      if( self % input_data % mode == INT_PARALLEL ) then
!!$         do irank = 0,nrank-1
!!$            !
!!$            ! Source
!!$            !
!!$            kk = 0
!!$            do ii = 1,memory_size(lista_recv(irank) % l)
!!$               if( lista_recv(irank) % l(ii) /= 0 ) then
!!$                  kk = kk + 1
!!$                  write(6,*)'Nearest, I am source of point=',kfl_paral,'->',irank,':',lista_recv(irank) % l(ii)
!!$               end if
!!$            end do
!!$            !call memory_alloca(memor_loc,'SELF % SOURCE_XX',vacal,self % source_xx(irank) % l,kk)
!!$            !kk = 0
!!$            !do ii = 1,memory_size(lista_recv(irank) % l)
!!$            !   if( lista_recv(irank) % l(ii) /= 0 ) then
!!$            !      kk = kk + 1
!!$            !      self % source_xx(irank) % l(kk) = lista_recv(irank) % l(ii)
!!$            !   end if
!!$            !end do
!!$            !
!!$            ! Target
!!$            !
!!$            kk = 0
!!$            do ii = 1,memory_size(lista_send(irank) % l)
!!$               if( lista_send(irank) % l(ii) /= 0 ) then
!!$                  kk = kk + 1
!!$                  write(6,*)'Nearest, I am target of point=',kfl_paral,'<-',irank,':',lista_send(irank) % l(ii)
!!$               end if
!!$            end do
!!$            !call memory_alloca(memor_loc,'SELF % TARGET_XX',vacal,self % target_xx(irank) % l,kk)
!!$            !kk = 0
!!$            !do ii = 1,memory_size(lista_send(irank) % l)
!!$            !   if( lista_send(irank) % l(ii) /= 0 ) then
!!$            !      kk = kk + 1
!!$            !      self % target_xx(irank) % l(kk) = lista_send(irank) % l(ii)
!!$            !   end if
!!$            !end do
!!$         end do
!!$      end if
!!$
!!$    end block
    
    !--------------------------------------------------------------------
    !
    ! Deallocate
    !
    !--------------------------------------------------------------------

    call memory_deallo(memor_loc,'DISTA_RECV'   ,vacal,dista_recv)
    call memory_deallo(memor_loc,'LISTA_SEND'   ,vacal,lista_send)
    call memory_deallo(memor_loc,'LISTA_RECV'   ,vacal,lista_recv)
    call memory_deallo(memor_loc,'XX_RECV'      ,vacal,xx_recv   )
    call memory_deallo(memor_loc,'NN_SEND'      ,vacal,nn_send   )
    call memory_deallo(memor_loc,'NN_RECV'      ,vacal,nn_recv   )

    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)

  end subroutine interpolation_nearest_node
    
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-29
  !> @brief   Nearest node and global numbering
  !> @details Global numbering
  !> 
  !-----------------------------------------------------------------------

  subroutine interpolation_global_numbering(&
    self,xx,ll,mesh,shapf,deriv,lelem,lenty,llost,mask,MEMORY_COUNTER,search_seq_in,search_par_in)

    class(interpolation),                    intent(inout) :: self
    real(rp),                       pointer, intent(in)    :: xx(:,:)
    integer(ip),                    pointer, intent(in)    :: ll(:)
    class(mesh_type_basic),                  intent(in)    :: mesh
    type(r2p),                      pointer, intent(inout) :: shapf(:)
    type(r3p),                      pointer, intent(inout) :: deriv(:)
    type(i1p),                      pointer, intent(inout) :: lelem(:)
    type(i2p),                      pointer, intent(inout) :: lenty(:)
    logical(lg),                    pointer, intent(inout) :: llost(:)
    logical(lg),          optional, pointer, intent(in)    :: mask(:)
    integer(8),           optional,          intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    class(search_method), optional, target,  intent(inout) :: search_seq_in     ! Sequential search
    class(search_method), optional, target,  intent(inout) :: search_par_in     ! Sequential search   
    integer(ip)                                            :: irank,nrank,nn
    integer(ip)                                            :: ii,kk,ineig,ipoin
    integer(ip)                                            :: my_rank
    real(rp)                                               :: toler
    integer(8)                                             :: memor_loc(2)
    type(i1p),                      pointer                :: ll_recv(:)
    type(r2p),                      pointer                :: xx_recv(:)
    integer(ip),                    pointer                :: nn_send(:)
    integer(ip),                    pointer                :: nn_recv(:)
    type(i1p),                      pointer                :: lista_send(:)  
    type(i1p),                      pointer                :: lista_recv(:)
    class(search_method),           pointer                :: search_seq        ! Sequential search
    class(search_method),           pointer                :: search_par        ! Sequential search   
    type(hash_t)                                           :: ht

    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)

    nullify(xx_recv)
    nullify(ll_recv)
    nullify(nn_send)
    nullify(nn_recv)
    nullify(lista_send)
    nullify(lista_recv)

    !--------------------------------------------------------------------
    !
    ! Search methods
    !
    !--------------------------------------------------------------------
    
    if( present(search_seq_in) ) then
       search_seq => search_seq_in
    else
       search_seq => self % search_method_seq
    end if
    
    if( present(search_par_in) ) then
       search_par => search_par_in
    else
       search_par => self % parallel_search % search_method_par
    end if
    
    !--------------------------------------------------------------------
    !
    ! Initialization
    !
    !--------------------------------------------------------------------

    nn      = memory_size(xx,2_ip)
    toler   = self % input_data % toler_rel
    nrank   = self % nrank
    my_rank = int(self % parallel_search % comm % RANK4,ip)

    !--------------------------------------------------------------------
    !
    ! Hash table
    !
    !--------------------------------------------------------------------

    if( nn > 0 ) then
       call htaini(ht,nn,lidson=.true.,AUTOMATIC_SIZE=.true.)
       if( associated(ll) ) call htaadd(ht,nn,ll)
    end if

    !--------------------------------------------------------------------
    !
    ! Search
    !
    !--------------------------------------------------------------------

    if( self % input_data % mode == INT_PARALLEL ) then

#ifdef __PGI
       call self % parallel_search % send_list_recv_points(&
            XX=xx,LL=ll,NN_SEND=nn_send,NN_RECV=nn_recv,LISTA_SEND=lista_send,XX_RECV=xx_recv,&
            LL_RECV=ll_recv,METHOD=CANDIDATE_NEAREST,MASK=std_log_1,&
            EXCLUDE_MYSELF=.not.self % input_data % myself,MEMORY_COUNTER=memor_loc)
#else
       call self % parallel_search % send_list_recv_points(&
            xx,ll,nn_send,nn_recv,lista_send,xx_recv,ll_recv,&
            METHOD=CANDIDATE_NEAREST,&
            EXCLUDE_MYSELF=.not.self % input_data % myself,MEMORY_COUNTER=memor_loc)
#endif
       do irank = 0,nrank-1
          if( nn_recv(irank) > 0 ) then
             call global_numbering(&
                  self,                     &
                  search_seq,               &
                  xx_recv(irank)%a,         &
                  ll_recv(irank)%l,         &
                  ht,                       &
                  mesh,                     &
                  lelem(irank)%l,           &
                  shapf(irank)%a,           &
                  deriv(irank)%a,           &
                  lenty(irank)%l,           &
                  mask,                     &
                  toler,                    &
                  memor_loc                 )
          end if
       end do   
       !
       ! Communication strategy
       !
       call self % parallel_search % dista_comm(nn,nn_send,nn_recv,lista_send,lista_recv,lelem,MEMORY_COUNTER=memor_loc)
       !
       do irank = 0,self % nrank-1
          do ii = 1,memory_size(lista_recv(irank) % l)
             if( lista_recv(irank) % l(ii) == 0 ) then
                lelem(irank) % l(ii) = 0
             end if
          end do
       end do
       
       kk = 0
       do ineig = 1,self % parallel_search % comm % nneig
          do ii = 1,self % parallel_search % comm % lrecv_size(ineig+1)-self % parallel_search % comm % lrecv_size(ineig)
             kk           = kk + 1
             ipoin        = self % parallel_search % comm % lrecv_perm(kk)
             llost(ipoin) = .false.
          end do
       end do
    else

       call global_numbering(         &
            self,                     &
            search_seq,               &
            xx,                       &
            ll,                       &
            ht,                       &
            mesh,                     &
            lelem(my_rank)%l,         &
            shapf(my_rank)%a,         &
            deriv(my_rank)%a,         &
            lenty(my_rank)%l,         &
            mask,                     &
            toler,                    &
            memor_loc                 )

       do ii = 1,memory_size(lelem(my_rank) % l)
          if( lelem(my_rank) % l(ii) /= 0 ) llost(ii) = .false.
       end do

    end if
        
    !--------------------------------------------------------------------
    !
    ! Deallocate
    !
    !--------------------------------------------------------------------

    if( nn > 0 ) call htades(ht)
    
    call memory_deallo(memor_loc,'LISTA_SEND'   ,vacal,lista_send)
    call memory_deallo(memor_loc,'LISTA_RECV'   ,vacal,lista_recv)
    call memory_deallo(memor_loc,'LL_RECV'      ,vacal,ll_recv   )
    call memory_deallo(memor_loc,'XX_RECV'      ,vacal,xx_recv   )
    call memory_deallo(memor_loc,'NN_SEND'      ,vacal,nn_send   )
    call memory_deallo(memor_loc,'NN_RECV'      ,vacal,nn_recv   )

    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)

  end subroutine interpolation_global_numbering
    
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-29
  !> @brief   Element interpolation
  !> @details Element interpolation
  !> 
  !-----------------------------------------------------------------------

  subroutine interpolation_entity(&
       self,xx,mesh,shapf,deriv,lelem,lenty,llost,mask,MEMORY_COUNTER)

    class(interpolation),                    intent(inout) :: self
    real(rp),                       pointer, intent(in)    :: xx(:,:)
    class(mesh_type_basic),                  intent(in)    :: mesh
    type(r2p),                      pointer, intent(inout) :: shapf(:)
    type(r3p),                      pointer, intent(inout) :: deriv(:)
    type(i1p),                      pointer, intent(inout) :: lelem(:)
    type(i2p),                      pointer, intent(inout) :: lenty(:)
    logical(lg),                    pointer, intent(inout) :: llost(:)
    logical(lg),          optional, pointer, intent(in)    :: mask(:)
    integer(8),           optional,          intent(inout) :: MEMORY_COUNTER(2) !< Memory counters

    type(r1p),                      pointer                :: dista_recv(:)    
    type(r2p),                      pointer                :: xx_recv(:)
    integer(ip),                    pointer                :: nn_send(:)
    integer(ip),                    pointer                :: nn_recv(:)
    type(i1p),                      pointer                :: lista_send(:)  
    type(i1p),                      pointer                :: lista_recv(:)

    real(rp)                                               :: toler
    integer(ip)                                            :: ii,kk,irank
    integer(ip)                                            :: ineig,ipoin
    integer(ip)                                            :: my_rank,nn,nrank
    integer(ip)                                            :: kfl_method
    integer(8)                                             :: memor_loc(2)
    logical(lg)                                            :: myself_first
    procedure(interpolation_generic),      pointer         :: interpolation_what

    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)

    nullify(dista_recv)    
    nullify(xx_recv)    
    nullify(nn_send)
    nullify(nn_recv)
    nullify(lista_send)  
    nullify(lista_recv)

    !--------------------------------------------------------------------
    !
    ! What to search
    !
    !--------------------------------------------------------------------

    if(      self % input_data % interpolation_method == INT_BOUNDARY_INTERPOLATION ) then
       interpolation_what => interpolation_boundary
       kfl_method         =  CANDIDATE_NEAREST
    else if( self % input_data % interpolation_method == INT_BOUNDARY_VALUE         ) then
       interpolation_what => interpolation_boundary
       kfl_method         =  CANDIDATE_NEAREST
    else if( self % input_data % interpolation_method == INT_ELEMENT_INTERPOLATION  ) then
       interpolation_what => interpolation_element
       kfl_method         =  CANDIDATE_INSIDE
    else if( self % input_data % interpolation_method == INT_ELEMENT_VALUE          ) then
       interpolation_what => interpolation_element
       kfl_method         =  CANDIDATE_INSIDE
     else
       call runend('DEF_INTERPOLATION_METHOD: DO NOT KNOW WHERE TO INTERPOLATE FROM')
    end if

    !--------------------------------------------------------------------
    !
    ! Allocate
    !
    !--------------------------------------------------------------------

    nn           = memory_size(xx,2_ip)
    toler        = self % input_data % toler_rel
    nrank        = self % nrank
    my_rank      = int(self % parallel_search % comm % RANK4,ip)
    if( kfl_method == CANDIDATE_INSIDE ) then
       myself_first = .true.                   ! Check in my subdomain first
    else
       myself_first = .false.                  ! For nearest, this does not make sense
    end if
    
    call memory_alloca(memor_loc,'DISTA_RECV',vacal,dista_recv,nrank,LBOUN=0_ip)

    !--------------------------------------------------------------------
    !
    ! Try myself first
    !
    !--------------------------------------------------------------------
    
    if( myself_first .or. self % input_data % mode == INT_SEQUENTIAL ) then
       call interpolation_what(&
            self                    , &
            self % search_method_seq, &
            xx,                       &
            mesh,                     &
            lelem(my_rank)%l,         &
            shapf(my_rank)%a,         &
            deriv(my_rank)%a,         &
            dista_recv(my_rank) % a,  &
            lenty(my_rank)%l,         &
            mask,                     &
            toler,                    &
            memor_loc                 )
       do ii = 1,memory_size(lelem(my_rank) % l)
          if( lelem(my_rank) % l(ii) /= 0 ) then
             llost(ii) = .false.
             self % dista(ii) = dista_recv(my_rank) % a(ii)
          end if
       end do
       if( nn > 0 ) then
          self % stats(3) = real(nn-int(count(llost),ip),rp)
       end if
    end if

    !--------------------------------------------------------------------
    !
    ! Parallel search
    !
    !--------------------------------------------------------------------

    if( self % input_data % mode == INT_PARALLEL ) then
       !
       ! Parallel mode
       !
#ifdef __PGI
       call self % parallel_search % send_list_recv_points(&
            XX=xx,NN_SEND=nn_send,NN_RECV=nn_recv,LISTA_SEND=lista_send,XX_RECV=xx_recv,&
            METHOD=kfl_method,MASK=llost,MEMORY_COUNTER=memor_loc)
#else
       call self % parallel_search % send_list_recv_points(&
            xx,nn_send,nn_recv,lista_send,xx_recv,&
            METHOD=kfl_method,MASK=llost,MEMORY_COUNTER=memor_loc)
#endif

       if( myself_first ) then
          nn_recv(my_rank) = nn
          nn_send(my_rank) = nn
          call memory_deallo(memor_loc,'LISTA_SEND % L',vacal,lista_send(my_rank) % l)
          call memory_alloca(memor_loc,'LISTA_SEND % L',vacal,lista_send(my_rank) % l,nn)
          do ii = 1,nn
             lista_send(my_rank) % l(ii) = ii
          end do
          !nn_recv(my_rank) = 0
          !nn_send(my_rank) = 0
          !call memory_deallo(memor_loc,'LISTA_SEND % L',vacal,lista_send(my_rank) % l)
          !call memory_alloca(memor_loc,'LISTA_SEND % L',vacal,lista_send(my_rank) % l,nn)
          !do ii = 1,nn
          !   lista_send(my_rank) % l(ii) = ii
          !end do
       end if

!!$       block
!!$         use def_master
!!$         integer(ip) :: irank
!!$         do irank = 0,nrank - 1
!!$            if( nn_send(irank) > 0 ) write(6,*)'send=',kfl_paral,'->',irank,': ',nn_send(irank)!,memory_size(dista_recv(irank) % a)
!!$            if( nn_recv(irank) > 0 ) write(6,*)'recv=',kfl_paral,'<-',irank,': ',nn_recv(irank)!,memory_size(dista_send(irank) % a)
!!$         end do
!!$         !call runend('O.K.!')
!!$       end block
       
       !irank = 0
       !block ; use def_master ; write(6,*)'cacacac=',kfl_paral,nn_recv(irank),memory_size(xx_recv(irank)%a); end block
       !  call runend('O.K.!')
       !block
       !  use def_master
       !  do irank = 0,nrank - 1
       !     write(6,*)'dista1=',kfl_paral,my_rank,'->',irank,': ',nn_recv(irank),memory_size(dista_recv(irank)%a)
       !  end do
       !  call runend('O.K.!')
       !end block
       
       do irank = 0,nrank-1

          if( nn_recv(irank) > 0 .and. ( irank /= my_rank .or. (.not. myself_first) ) ) then

             call interpolation_what(&
                  self,                     &
                  self % search_method_seq, &
                  xx_recv(irank)%a,         &
                  mesh,                     &
                  lelem(irank)%l,           &
                  shapf(irank)%a,           &
                  deriv(irank)%a,           &
                  dista_recv(irank) % a,    &
                  lenty(irank)%l,           &
                  mask,                     &
                  toler,                    &
                  memor_loc                 )

          end if
       end do

!!$       block
!!$         use def_master
!!$         do irank = 0,nrank - 1
!!$            write(6,*)'dista2=',kfl_paral,'->',irank,': ',nn_recv(irank),memory_size(dista_recv(irank)%a)
!!$         end do
!!$         call runend('O.K.!')
!!$       end block

       call self % parallel_search % dista_comm(nn,nn_send,nn_recv,lista_send,lista_recv,dista_recv,self % dista,MEMORY_COUNTER=memor_loc) 

       do irank = 0,nrank-1
          do ii = 1,memory_size(lista_recv(irank) % l)
             if( lista_recv(irank) % l(ii) == 0 ) then
                lelem(irank) % l(ii) = 0
             end if
          end do
       end do
       kk = 0
       do ineig = 1,self % parallel_search % comm % nneig
          do ii = 1,self % parallel_search % comm % lrecv_size(ineig+1)-self % parallel_search % comm % lrecv_size(ineig)
             kk    = kk + 1
             ipoin = self % parallel_search % comm % lrecv_perm(kk)
             llost(ipoin) = .false.
          end do
       end do
    end if

    !--------------------------------------------------------------------
    !
    ! Assign target and source points
    !
    !--------------------------------------------------------------------

!!$    call memory_alloca(memor_loc,'SELF % SOURCE_XX',vacal,self % source_xx,nrank,LBOUN=0_ip)
!!$    call memory_alloca(memor_loc,'SELF % TARGET_XX',vacal,self % target_xx,nrank,LBOUN=0_ip)
!!$
!!$    block
!!$      use def_master, only : kfl_paral
!!$      if( self % input_data % mode == INT_PARALLEL ) then
!!$         do irank = 0,nrank-1
!!$            !
!!$            ! Source
!!$            !
!!$            kk = 0
!!$            do ii = 1,memory_size(lista_recv(irank) % l)
!!$               if( lista_recv(irank) % l(ii) /= 0 ) then
!!$                  kk = kk + 1
!!$                  write(6,*)'Interpo, I am source of point=',kfl_paral,'->',irank,':',lista_recv(irank) % l(ii)
!!$               end if
!!$            end do
!!$            !call memory_alloca(memor_loc,'SELF % SOURCE_XX',vacal,self % source_xx(irank) % l,kk)
!!$            !kk = 0
!!$            !do ii = 1,memory_size(lista_recv(irank) % l)
!!$            !   if( lista_recv(irank) % l(ii) /= 0 ) then
!!$            !      kk = kk + 1
!!$            !      self % source_xx(irank) % l(kk) = lista_recv(irank) % l(ii)
!!$            !   end if
!!$            !end do
!!$            !
!!$            ! Target
!!$            !
!!$            kk = 0
!!$            do ii = 1,memory_size(lista_send(irank) % l)
!!$               if( lista_send(irank) % l(ii) /= 0 ) then
!!$                  kk = kk + 1
!!$                  write(6,*)'Interpo, I am target of point=',kfl_paral,'<-',irank,':',lista_send(irank) % l(ii)
!!$               end if
!!$            end do
!!$            !call memory_alloca(memor_loc,'SELF % TARGET_XX',vacal,self % target_xx(irank) % l,kk)
!!$            !kk = 0
!!$            !do ii = 1,memory_size(lista_send(irank) % l)
!!$            !   if( lista_send(irank) % l(ii) /= 0 ) then
!!$            !      kk = kk + 1
!!$            !      self % target_xx(irank) % l(kk) = lista_send(irank) % l(ii)
!!$            !   end if
!!$            !end do
!!$         end do
!!$      end if
!!$     end block

    !--------------------------------------------------------------------
    !
    ! Deallocate
    !
    !--------------------------------------------------------------------

    call memory_deallo(memor_loc,'DISTA_RECV'   ,vacal,dista_recv)
    call memory_deallo(memor_loc,'LISTA_SEND'   ,vacal,lista_send)
    call memory_deallo(memor_loc,'LISTA_RECV'   ,vacal,lista_recv)
    call memory_deallo(memor_loc,'XX_RECV'      ,vacal,xx_recv   )
    call memory_deallo(memor_loc,'NN_SEND'      ,vacal,nn_send   )
    call memory_deallo(memor_loc,'NN_RECV'      ,vacal,nn_recv   )

    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)

  end subroutine interpolation_entity
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-29
  !> @brief   Sparse matrix
  !> @details Compute interpolation sparse matrix
  !>
  !>          x        x  x
  !>          o-----o-----o-----o
  !>          1     2     3     4
  !>
  !>          1.0  0.0  0.0  0.0 
  !>          0.0  0.5  0.5  0.0
  !>          0.0  0.0  1.0  0.0
  !>
  !>          1.0 0.5 0.5 1.0
  !>          xa = (1,2,2,3)
  !>          ya = (1,2,3,3) 
  !>
  !-----------------------------------------------------------------------

  subroutine sparse_matrix(self,lelem,lenty,shapf,deriv,MEMORY_COUNTER)

    class(interpolation),                       intent(inout) :: self
    type(i1p),                        pointer,  intent(in)    :: lelem(:)
    type(i2p),                        pointer,  intent(in)    :: lenty(:)
    type(r2p),                        pointer,  intent(in)    :: shapf(:)
    type(r3p),                        pointer,  intent(in)    :: deriv(:)
    integer(8),             optional,           intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(ip)                                               :: irank,nz,nd
    integer(ip)                                               :: ii,kk,ielem
    integer(ip)                                               :: pnode
    integer(ip)                                               :: ipoin,inode
    integer(8)                                                :: memor_loc(2)

    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)

    nd = self % nd
    
    do irank = 0,self % nrank-1
       !
       ! Compute dimension
       !
       nz = 0
       if( self % input_data % interpolation_method == INT_ELEMENT_VALUE .or. self % input_data % interpolation_method == INT_BOUNDARY_VALUE ) then
          do ii = 1,memory_size(lelem(irank) % l)
             ielem = lelem(irank) % l(ii)
             if( ielem /= 0 ) then
                pnode = 1
                nz    = nz + 1
             end if
          end do
       else
          do ii = 1,memory_size(lelem(irank) % l)
             ielem = lelem(irank) % l(ii)
             if( ielem /= 0 ) then
                pnode = maths_maxloc_nonzero(lenty(irank) % l(:,ii))
                nz    = nz + pnode
             end if
          end do
       end if
       !
       ! Allocate matrix for values and derivatives
       !
       self % matrix(irank) % ndof1 = 1
       self % matrix(irank) % ndof2 = 1
       self % matrix(irank) % nz    = nz       
       call self % matrix(irank) % alloca(MEMORY_COUNTER=memor_loc)

       if( self % input_data % deriv ) then
          self % matder(irank) % ndof1 = 1
          self % matder(irank) % ndof2 = nd
          self % matder(irank) % nz    = nz
          call self % matder(irank) % alloca(MEMORY_COUNTER=memor_loc)
       end if
       !
       ! Fill in matrix
       !
       nz = 0
       kk = 0

       if( self % input_data % interpolation_method == INT_ELEMENT_VALUE .or. self % input_data % interpolation_method == INT_BOUNDARY_VALUE ) then
          do ii = 1,memory_size(lelem(irank) % l)
             ielem = lelem(irank) % l(ii)
             if( ielem /= 0 ) then
                kk                                       = kk + 1
                nz                                       = nz + 1                

                self % matrix(irank) % xA(nz)            = kk
                self % matrix(irank) % yA(nz)            = ielem
                self % matrix(irank) % vA(1,1,nz)        = 1.0_rp
                self % matrix(irank) % nrows             = max(self % matrix(irank) % nrows,kk)
                self % matrix(irank) % ncols             = max(self % matrix(irank) % ncols,ielem)

                
                if( self % input_data % deriv ) then
                   self % matder(irank) % xA(nz)         = kk
                   self % matder(irank) % yA(nz)         = ielem
                   self % matder(irank) % vA(1,1:nd,nz)  = 0.0_rp
                   self % matder(irank) % nrows          = max(self % matder(irank) % nrows,kk)
                   self % matder(irank) % ncols          = max(self % matder(irank) % ncols,ielem)
                end if
             end if
          end do
       else
          do ii = 1,memory_size(lelem(irank) % l)
             ielem = lelem(irank) % l(ii)
             if( ielem /= 0 ) then

                kk    = kk + 1
                pnode = maths_maxloc_nonzero(lenty(irank) % l(:,ii))
                do inode = 1,pnode
                   nz                                       = nz + 1                
                   ipoin                                    = lenty(irank) % l(inode,ii)

                   self % matrix(irank) % xA(nz)            = kk
                   self % matrix(irank) % yA(nz)            = ipoin
                   self % matrix(irank) % vA(1,1,nz)        = shapf(irank) % a(inode,ii)
                   self % matrix(irank) % nrows             = max(self % matrix(irank) % nrows,kk)
                   self % matrix(irank) % ncols             = max(self % matrix(irank) % ncols,ipoin)

                   if( self % input_data % deriv ) then
                      self % matder(irank) % xA(nz)         = kk
                      self % matder(irank) % yA(nz)         = ipoin
                      self % matder(irank) % vA(1,1:nd,nz)  = deriv(irank) % a(1:nd,inode,ii)
                      self % matder(irank) % nrows          = max(self % matder(irank) % nrows,kk)
                      self % matder(irank) % ncols          = max(self % matder(irank) % ncols,ipoin)
                   end if
                end do
             end if
          end do
       end if
    end do

    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)

  end subroutine sparse_matrix
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-29
  !> @brief   Sparse matrix
  !> @details Permute sparse matrix
  !>
  !-----------------------------------------------------------------------

  subroutine sparse_matrix_permutation(self,perm)

    class(interpolation),                       intent(inout) :: self
    integer(ip),                      pointer,  intent(in)    :: perm(:)
    integer(ip)                                               :: irank
    integer(ip)                                               :: iz
    integer(ip)                                               :: ii_old,ii_new

    do irank = 0,self % nrank-1
       do iz = 1,self % matrix(irank) % nz
          ii_old = self % matrix(irank) % yA(iz)
          if( ii_old /= 0 ) then
             ii_new = perm(ii_old) 
             self % matrix(irank) % yA(iz)    = ii_new
             if( self % input_data % deriv ) then
                self % matder(irank) % yA(iz) = ii_new
             end if
          end if
       end do
    end do

  end subroutine sparse_matrix_permutation
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-18
  !> @brief   Compute value
  !> @details Compute value
  !> 
  !-----------------------------------------------------------------------
  
  subroutine values_10(self,xx_in,xx_out,POINT,MEMORY_COUNTER)
    
    class(interpolation),                   intent(inout) :: self
    real(rp),                      pointer, intent(in)    :: xx_in(:)
    real(rp),                               intent(inout) :: xx_out
    integer(ip),          optional,         intent(in)    :: POINT
    integer(8),           optional,         intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    real(rp),                      pointer                :: xx_tmp(:)
    real(rp)                                              :: time1,time2
    
    call cputim(time1)
    
    if( self % input_data % mode == INT_PARALLEL ) then
       !
       ! Parallel mode
       !
    else    
       !
       ! Sequential mode
       !
       allocate(xx_tmp(1))
       
       call self % matrix(0) % mv(xx_in,xx_tmp,N1=POINT,N2=POINT)
       xx_out = xx_tmp(1)
       deallocate(xx_tmp)
    end if
    
    call cputim(time2) ; self % times(3) = self % times(3) + time2-time1 
   
  end subroutine values_10
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-18
  !> @brief   Compute value
  !> @details Compute value
  !> 
  !-----------------------------------------------------------------------
  
  subroutine values_10_i(self,xx_in,xx_out,POINT,MEMORY_COUNTER)
    
    class(interpolation),                   intent(inout) :: self
    integer(ip),                   pointer, intent(in)    :: xx_in(:)
    integer(ip),                            intent(inout) :: xx_out
    integer(ip),          optional,         intent(in)    :: POINT
    integer(8),           optional,         intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    real(rp),                      pointer                :: xx_in_r(:)
    real(rp)                                              :: xx_out_r
    integer(ip)                                           :: ii,nn_in
    integer(8)                                            :: memor_loc(2)
    real(rp)                                              :: time1,time2

    call cputim(time1)
    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)
    
    nullify(xx_in_r)

    nn_in = memory_size(xx_in,1_ip)
    call memory_alloca(memor_loc,'XX_IN_R' ,vacal,xx_in_r ,nn_in)
    
    do ii = 1,nn_in
       xx_in_r(ii) = real(xx_in(ii),rp)
    end do
    call self % values(xx_in_r,xx_out_r,POINT,MEMORY_COUNTER=memor_loc)
    xx_out = int(xx_out_r,ip)
    
    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)   
    call cputim(time2) ; self % times(3) = self % times(3) + time2-time1 
    
  end subroutine values_10_i

  subroutine values_11_i(self,xx_in,xx_out,POINT,INITIALIZATION,MEMORY_COUNTER)
    
    class(interpolation),                   intent(inout) :: self
    integer(ip),                   pointer, intent(in)    :: xx_in(:)
    integer(ip),                   pointer, intent(inout) :: xx_out(:)
    integer(ip),          optional,         intent(in)    :: POINT
    logical(lg),          optional,         intent(in)    :: INITIALIZATION
    integer(8),           optional,         intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(ip)                                           :: ii,nn_in,nn_out
    real(rp),                      pointer                :: xx_in_r(:)
    real(rp),                      pointer                :: xx_out_r(:)
    integer(8)                                            :: memor_loc(2)
    real(rp)                                              :: time1,time2

    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)
    call cputim(time1)

    nullify(xx_in_r,xx_out_r)

    nn_in  = memory_size(xx_in ,1_ip)
    nn_out = memory_size(xx_out,1_ip)
    
    call memory_alloca(memor_loc,'XX_IN_R' ,vacal,xx_in_r ,nn_in)
    call memory_alloca(memor_loc,'XX_OUT_R',vacal,xx_out_r,nn_out)

    do ii = 1,nn_in
       xx_in_r(ii) = real(xx_in(ii),rp)
    end do
    call self % values(xx_in_r,xx_out_r,POINT,MEMORY_COUNTER=memor_loc)
    if( associated(self % perm) ) then
       do ii = 1,nn_out
          xx_out(self % perm(ii)) = int(xx_out_r(ii),ip)
       end do       
    else
       do ii = 1,nn_out
          xx_out(ii) = int(xx_out_r(ii),ip)
       end do
    end if
    
    call memory_deallo(memor_loc,'XX_IN_R' ,vacal,xx_in_r )
    call memory_deallo(memor_loc,'XX_OUT_R',vacal,xx_out_r)
    
    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)
    call cputim(time2) ; self % times(3) = self % times(3) + time2-time1 

  end subroutine values_11_i
  
  subroutine values_11(self,xx_in,xx_out,POINT,INITIALIZATION,MEMORY_COUNTER)
    
    class(interpolation),                   intent(inout) :: self
    real(rp),                      pointer, intent(in)    :: xx_in(:)
    real(rp),                      pointer, intent(inout) :: xx_out(:)
    integer(ip),          optional,         intent(in)    :: POINT
    logical(lg),          optional,         intent(in)    :: INITIALIZATION
    integer(8),           optional,         intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(ip)                                           :: ii,ipoin
    integer(8)                                            :: memor_loc(2)
    integer(ip)                                           :: nrank,irank,ineig,jj
    integer(ip)                                           :: nsend,nrecv
    real(rp),                      pointer                :: xx_send(:)
    real(rp),                      pointer                :: xx_recv(:)
    logical(lg)                                           :: if_initialization
    real(rp)                                              :: time1,time2
   
    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)
    call cputim(time1)
    !
    ! Initialize
    !
    if_initialization = optional_argument(.true.,INITIALIZATION)
    if( if_initialization .and. associated(xx_out) ) then
       if( present(POINT) ) then          
          xx_out(1) = 0
       else
          do ii = 1,memory_size(xx_out,1_ip)
             xx_out(ii) = 0.0_rp
          end do
       end if
    end if
    
    if( self % input_data % mode == INT_PARALLEL ) then
       !
       ! Parallel mode
       !
       nullify(xx_send)
       nullify(xx_recv)
       nrank = self % nrank
       jj    = 0
       do ineig = 1,self % parallel_search % comm % nneig
          nsend = self % parallel_search % comm % lsend_size(ineig+1) - self % parallel_search % comm % lsend_size(ineig)
          nrecv = self % parallel_search % comm % lrecv_size(ineig+1) - self % parallel_search % comm % lrecv_size(ineig)
          if( nsend /= 0 .or. nrecv /= 0 ) then
             irank = self % parallel_search % comm % neights(ineig)
             call memory_alloca(memor_loc,'XX_SEND',vacal,xx_send,nsend)
             call memory_alloca(memor_loc,'XX_RECV',vacal,xx_recv,nrecv,'DO_NOT_INITIALIZE')
             call self % matrix(irank) % mv(xx_in,xx_send)   
             call PAR_SEND_RECEIVE(xx_send,xx_recv,DOM_I=irank,PAR_COMM_IN=self % parallel_search % comm % PAR_COMM_WORLD)
             if( present(POINT) ) then
                do ii = 1,nrecv
                   jj                 = jj + 1
                   ipoin              = self % parallel_search % comm % lrecv_perm(jj)
                   if( ipoin == POINT ) then
                      xx_out(1) = xx_out(1) + xx_recv(ii)
                   end if
                end do
             else
                do ii = 1,nrecv
                   jj            = jj + 1
                   ipoin         = self % parallel_search % comm % lrecv_perm(jj)
                   if( associated(self % perm) ) then
                      xx_out(self % perm(ipoin)) = xx_out(self % perm(ipoin)) + xx_recv(ii)
                   else
                      xx_out(ipoin) = xx_out(ipoin) + xx_recv(ii)
                   end if
                end do
             end if
             call memory_deallo(memor_loc,'XX_SEND',vacal,xx_send)
             call memory_deallo(memor_loc,'XX_RECV',vacal,xx_recv)
          end if
       end do
    else    
       !
       ! Sequential mode
       !       
       if( associated(self % perm) ) then
          call memory_copy(memor_loc,'XX_RECV',vacal,xx_out,xx_recv,'DO_NOT_DEALLOCATE')
          call self % matrix(0) % mv(xx_in,xx_recv,N1=POINT,N2=POINT,INITIALIZATION=INITIALIZATION)
          do ii = 1,memory_size(xx_out)
             xx_out(self % perm(ii)) = xx_recv(ii)
          end do
          call memory_deallo(memor_loc,'XX_RECV',vacal,xx_recv)
       else
          call self % matrix(0) % mv(xx_in,xx_out,N1=POINT,N2=POINT,INITIALIZATION=INITIALIZATION)
       end if
    end if

    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)
    call cputim(time2) ; self % times(3) = self % times(3) + time2-time1 
    
  end subroutine values_11

  subroutine values_22_i(self,xx_in,xx_out,POINT,INITIALIZATION,MEMORY_COUNTER)

    class(interpolation),                   intent(inout) :: self
    integer(ip),                   pointer, intent(in)    :: xx_in(:,:)
    integer(ip),                   pointer, intent(inout) :: xx_out(:,:)
    integer(ip),          optional,         intent(in)    :: POINT
    logical(lg),          optional,         intent(in)    :: INITIALIZATION
    integer(8),           optional,         intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(ip)                                           :: ii,n1_in,n1_out
    integer(ip)                                           :: n2_in,n2_out
    real(rp),                      pointer                :: xx_in_r(:,:)
    real(rp),                      pointer                :: xx_out_r(:,:)
    integer(8)                                            :: memor_loc(2)
    real(rp)                                              :: time1,time2
   
    call cputim(time1)
    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)

    nullify(xx_in_r,xx_out_r)

    n1_in  = memory_size(xx_in ,1_ip)
    n1_out = memory_size(xx_out,1_ip)
    n2_in  = memory_size(xx_in ,2_ip)
    n2_out = memory_size(xx_out,2_ip)

    call memory_alloca(memor_loc,'XX_IN_R' ,vacal,xx_in_r ,n1_in ,n2_in)
    call memory_alloca(memor_loc,'XX_OUT_R',vacal,xx_out_r,n1_out,n2_out)
    do ii = 1,n1_in
       xx_in_r(:,ii) = real(xx_in(:,ii),rp)
    end do
    call self % values(xx_in_r,xx_out_r,POINT,MEMORY_COUNTER=memor_loc)
    if( associated(self % perm) ) then
       do ii = 1,n1_out
          xx_out(:,self % perm(ii)) = int(xx_out_r(:,ii),ip)
       end do
    else
       do ii = 1,n1_out
          xx_out(:,ii) = int(xx_out_r(:,ii),ip)
       end do
    end if
    
    call memory_deallo(memor_loc,'XX_IN_R' ,vacal,xx_in_r )
    call memory_deallo(memor_loc,'XX_OUT_R',vacal,xx_out_r)
    
    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)
    call cputim(time2) ; self % times(3) = self % times(3) + time2-time1 
    
  end subroutine values_22_i
  
  subroutine values_22(self,xx_in,xx_out,POINT,INITIALIZATION,MEMORY_COUNTER)

    class(interpolation),                   intent(inout) :: self
    real(rp),                      pointer, intent(in)    :: xx_in(:,:)
    real(rp),                      pointer, intent(inout) :: xx_out(:,:)
    integer(ip),          optional,         intent(in)    :: POINT
    logical(lg),          optional,         intent(in)    :: INITIALIZATION
    integer(8),           optional,         intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(ip)                                           :: ii,ipoin,nd,nn
    integer(ip)                                           :: nrank,irank,ineig,jj
    integer(ip)                                           :: nsend,nrecv,idime
    integer(ip)                                           :: nn_out,nn_in
    integer(ip)                                           :: nd_out,nd_in
    integer(8)                                            :: memor_loc(2)
    real(rp),                      pointer                :: xx_send(:,:)
    real(rp),                      pointer                :: xx_recv(:,:)
    real(rp),                      pointer                :: x1_send(:)
    real(rp),                      pointer                :: x1_in(:)
    real(rp),                      pointer                :: x1_out(:)
    logical(lg)                                           :: if_initialization
    real(rp)                                              :: time1,time2
    
    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)
    nullify(x1_in)
    nullify(x1_out)
    call cputim(time1)

    nd     = max(memory_size(xx_in,1_ip),memory_size(xx_out,1_ip))
    nn     = self % nn
    nn_in  = memory_size(xx_in, 2_ip)
    nn_out = memory_size(xx_out,2_ip)
    nd_in  = memory_size(xx_in, 1_ip)
    nd_out = memory_size(xx_out,1_ip)

    if( nd_in /= 0 .and. nd_out /= 0 .and. nd_in /= nd_out ) &
         call runend('VALUES_22: WRONG DIMENSION 1 FOR OUTPUT ARRAY')
    if( nn_out /= 0 .and. nn_out < nn ) &
         call runend('VALUES_22: WRONG DIMENSION 2 FOR OUTPUT ARRAY')
    !
    ! Initialize
    !
    if_initialization = optional_argument(.true.,INITIALIZATION)
    if( if_initialization .and. associated(xx_out) ) then
       if( present(POINT) ) then          
          xx_out(:,1) = 0
       else
          do ii = 1,memory_size(xx_out,2_ip)
             xx_out(:,ii) = 0.0_rp
          end do
       end if
    end if
         
    if( self % input_data % mode == INT_PARALLEL ) then
       !
       ! Parallel mode
       !
       nullify(xx_send)
       nullify(xx_recv)
       nullify(x1_send)

       nrank  = self % nrank
       
       jj    = 0
       do ineig = 1,self % parallel_search % comm % nneig
          nsend = self % parallel_search % comm % lsend_size(ineig+1) - self % parallel_search % comm % lsend_size(ineig)
          nrecv = self % parallel_search % comm % lrecv_size(ineig+1) - self % parallel_search % comm % lrecv_size(ineig)
          if( nsend /= 0 .or. nrecv /= 0 ) then
             irank = self % parallel_search % comm % neights(ineig)
             call memory_alloca(memor_loc,'XX_SEND',vacal,xx_send,nd,nsend)
             call memory_alloca(memor_loc,'XX_RECV',vacal,xx_recv,nd,nrecv,'DO_NOT_INITIALIZE')             
             call memory_alloca(memor_loc,'X1_SEND',vacal,x1_send,nsend)
             call memory_alloca(memor_loc,'X1_IN',  vacal,x1_in  ,nn_in)
             do idime = 1,nd
                do ii = 1,nn_in
                   x1_in(ii) = xx_in(idime,ii)
                end do
                call self % matrix(irank) % mv(x1_in,x1_send)
                do ii = 1,nsend
                   xx_send(idime,ii) = x1_send(ii) 
                end do
             end do
             call PAR_SEND_RECEIVE(xx_send,xx_recv,DOM_I=irank,PAR_COMM_IN=self % parallel_search % comm % PAR_COMM_WORLD)
             if( present(POINT) ) then
                do ii = 1,nrecv
                   jj                 = jj + 1
                   ipoin              = self % parallel_search % comm % lrecv_perm(jj)
                   if( ipoin == POINT ) then
                      xx_out(1:nd,1) = xx_out(1:nd,1) + xx_recv(1:nd,ii)
                   end if
                end do
             else
                do ii = 1,nrecv
                   jj                 = jj + 1
                   ipoin              = self % parallel_search % comm % lrecv_perm(jj)
                   if( associated(self % perm) ) then
                      xx_out(1:nd,self % perm(ipoin)) = xx_out(1:nd,self % perm(ipoin)) + xx_recv(1:nd,ii)
                   else
                      xx_out(1:nd,ipoin) = xx_out(1:nd,ipoin) + xx_recv(1:nd,ii)
                   end if
                end do
             end if
             call memory_deallo(memor_loc,'XX_SEND',vacal,xx_send)
             call memory_deallo(memor_loc,'XX_RECV',vacal,xx_recv)
             call memory_deallo(memor_loc,'X1_SEND',vacal,x1_send)
             call memory_deallo(memor_loc,'X1_IN',  vacal,x1_in  )
          end if
       end do

    else
       !
       ! Sequential mode
       !
       call memory_alloca(memor_loc,'X1_OUT',vacal,x1_out,nn_out)
       call memory_alloca(memor_loc,'X1_IN' ,vacal,x1_in ,nn_in )
       do idime = 1,nd
          do ii = 1,nn_in
             x1_in(ii) = xx_in(idime,ii)
          end do
          call self % matrix(0) % mv(x1_in,x1_out,N1=POINT,N2=POINT)
          if( present(POINT) ) then
             xx_out(idime,1) = xx_out(idime,1) + x1_out(POINT)
          else
             if( associated(self % perm) ) then
                do ii = 1,nn_out
                   xx_out(idime,self % perm(ii)) = xx_out(idime,self % perm(ii)) + x1_out(ii)
                end do
             else
                do ii = 1,nn_out
                   xx_out(idime,ii) = xx_out(idime,ii) + x1_out(ii)
                end do
             end if
          end if
       end do
       call memory_deallo(memor_loc,'X1_OUT',vacal,x1_out)
       call memory_deallo(memor_loc,'X1_IN' ,vacal,x1_in )
       
    end if

    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)
    call cputim(time2) ; self % times(3) = self % times(3) + time2-time1 

  end subroutine values_22

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-18
  !> @brief   Compute derivative
  !> @details Compute derivative
  !> 
  !-----------------------------------------------------------------------
  
  subroutine derivatives_12(self,xx_in,xx_out,MEMORY_COUNTER)

    class(interpolation),                   intent(inout) :: self
    real(rp),                      pointer, intent(in)    :: xx_in(:)
    real(rp),                      pointer, intent(inout) :: xx_out(:,:)
    integer(8),           optional,         intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(ip)                                           :: ii,ipoin,nd,nn
    integer(ip)                                           :: nrank,irank,ineig,jj
    integer(ip)                                           :: nsend,nrecv,nn_in
    integer(ip)                                           :: nn_out
    integer(8)                                            :: memor_loc(2)
    real(rp),                      pointer                :: xx_send(:,:)
    real(rp),                      pointer                :: xx_recv(:,:)

    if( .not. self % input_data % deriv ) &
         call runend('DERIVATIVES_12: DERIVATIVES CANNOT BE COMPTUED. TURNON DERIVATIBES OPTIONS')
    
    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)

    nn     = self % nn
    nd     = memory_size(xx_out,1_ip)
    nn_in  = memory_size(xx_in ,2_ip)
    nn_out = memory_size(xx_out,2_ip)

    if( nn_out > 0 .and. self % matder(0) % ndof2 /= nd ) &
         call runend('DERIVATIVES_12: WRONG DIMENSION 1 FOR OUTPUT ARRAY')
    if( nn     > 0 .and. nn_out < nn ) &
         call runend('DERIVATIVES_12: WRONG DIMENSION 2 FOR OUTPUT ARRAY')
    
    if( self % input_data % mode == INT_PARALLEL ) then
       !
       ! Parallel mode
       !
       nullify(xx_send)
       nullify(xx_recv)

       nrank  = self % nrank       
       jj     = 0
       do ineig = 1,self % parallel_search % comm % nneig
          nsend = self % parallel_search % comm % lsend_size(ineig+1) - self % parallel_search % comm % lsend_size(ineig)
          nrecv = self % parallel_search % comm % lrecv_size(ineig+1) - self % parallel_search % comm % lrecv_size(ineig)
          if( nsend /= 0 .or. nrecv /= 0 ) then
             irank = self % parallel_search % comm % neights(ineig)
             call memory_alloca(memor_loc,'XX_SEND',vacal,xx_send,nd,nsend)
             call memory_alloca(memor_loc,'XX_RECV',vacal,xx_recv,nd,nrecv,'DO_NOT_INITIALIZE')  
             call self % matder(irank) % mv(xx_in,xx_send)
             call PAR_SEND_RECEIVE(xx_send,xx_recv,DOM_I=irank,PAR_COMM_IN=self % parallel_search % comm % PAR_COMM_WORLD)
             do ii = 1,nrecv
                jj                 = jj + 1
                ipoin              = self % parallel_search % comm % lrecv_perm(jj)
                xx_out(1:nd,ipoin) = xx_recv(1:nd,ii)
             end do
             call memory_deallo(memor_loc,'XX_SEND',vacal,xx_send)
             call memory_deallo(memor_loc,'XX_RECV',vacal,xx_recv)
          end if
       end do

    else
       !
       ! Sequential mode
       !
       call self % matder(0) % mv(xx_in,xx_out)
       
    end if

    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)

  end subroutine derivatives_12

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-05
  !> @brief   Search element
  !> @details Compute the distance
  !>          
  !-----------------------------------------------------------------------
  
  subroutine distances_21(self,xx_in,xx_out,POINT,MEMORY_COUNTER)

    class(interpolation),                   intent(inout) :: self
    real(rp),                      pointer, intent(in)    :: xx_in(:,:)
    real(rp),                      pointer, intent(inout) :: xx_out(:)
    integer(ip),          optional,         intent(in)    :: POINT
    integer(8),           optional,         intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(ip)                                           :: ii,nd,nn
    integer(8)                                            :: memor_loc(2)
    real(rp),                      pointer                :: xx_tmp(:,:)

    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)

    nullify(xx_tmp)
    nd = memory_size(xx_in,1_ip)
    nn = memory_size(xx_in,2_ip)
    
    call memory_alloca(memor_loc,'XX_TMP',vacal,xx_tmp,nd,nn)
    call self % values(xx_in,xx_tmp,MEMORY_COUNTER=memor_loc)

    do ii = 1,nn
       xx_out(ii) = sqrt(dot_product(xx_in(1:nd,ii)-xx_tmp(1:nd,ii),xx_in(1:nd,ii)-xx_tmp(1:nd,ii)))
    end do
    
    call memory_deallo(memor_loc,'XX_TMP',vacal,xx_tmp)
    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)
    
  end subroutine distances_21
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-05
  !> @brief   Search element
  !> @details Search element for a list of points XX
  !>          Output:
  !>          DISTA(:) ... -1 if point not found
  !>                        0 if point found
  !>          LELEM(:) ...  List of elements for each point
  !>          
  !-----------------------------------------------------------------------
  
  subroutine interpolation_element(self,search,xx,mesh,lelem,shapf,deriv,dista,lenty,mask,TOLER,MEMORY_COUNTER)

    class(interpolation),                    intent(inout) :: self
    class(search_method),                    intent(inout) :: search
    real(rp),                       pointer, intent(in)    :: xx(:,:)
    class(mesh_type_basic),                  intent(in)    :: mesh
    integer(ip),                    pointer, intent(inout) :: lelem(:)
    real(rp),                       pointer, intent(inout) :: shapf(:,:)
    real(rp),         optional,     pointer, intent(inout) :: deriv(:,:,:)
    real(rp),         optional,     pointer, intent(inout) :: dista(:)
    integer(ip),      optional,     pointer, intent(inout) :: lenty(:,:)
    logical(lg),      optional,     pointer, intent(in)    :: mask(:)
    real(rp),         optional,              intent(in)    :: TOLER
    integer(8),       optional,              intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(ip)                                            :: ndime,nn,ii,kk
    integer(ip)                                            :: nelem,mnode,pnode
    integer(ip)                                            :: inode,ielem,pelty
    integer(ip)                                            :: ifoun
    type(i1p),                      pointer                :: list_entities(:)
    integer(8)                                             :: memor_loc(2)
    real(rp)                                               :: toler_loc,time1,time2
    real(rp)                                               :: derit(mesh % ndime,mesh % mnode)
    real(rp)                                               :: elcod(mesh % ndime,mesh % mnode)
    real(rp)                                               :: coloc(3)
    logical(lg)                                            :: if_mask
    
    real(rp)                                               :: toler_loc_force
    integer(ip)                                            :: iter_force

    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)

    toler_loc = optional_argument(epsil,TOLER)
    if_mask   = .true.
    nn        = memory_size(xx,2_ip)
    mnode     = mesh % mnode
    ndime     = self % nd
    nullify(list_entities)

    if( nn > 0 ) then
       !
       ! Allocate if necessary
       !
       call interpolation_allocate(self,nn,ndime,mnode,lelem,shapf,deriv,dista,lenty,MEMORY_COUNTER=memor_loc)
      
     call cputim(time1)
     
#ifdef __PGI
       call search % candidate(XX=xx,LIST_ENTITIES=list_entities,METHOD=CANDIDATE_INSIDE,MASK=std_log_1,MEMORY_COUNTER=memor_loc)       
#else
       call search % candidate(xx,list_entities,METHOD=CANDIDATE_INSIDE,MEMORY_COUNTER=memor_loc)       
#endif
       
       call cputim(time2)
       self % times(1) = self % times(1) + time2-time1 ; time1 = time2

       do ii = 1,nn
          nelem = memory_size(list_entities(ii)%l)
         
          self % stats(2) = self % stats(2) + real(nelem,rp)
          if( present(dista) ) dista(ii) = -huge(1.0_rp)

          if( nelem == 0 ) then
             lelem(ii)   =  0
             shapf(:,ii) =  0.0_rp
             if( present(deriv) ) deriv(:,:,ii) = 0.0_rp
          else
            loop_kk: do kk = 1,nelem
               ielem = list_entities(ii) % l(kk)
               if( ielem <= 0 ) call runend('DEF_INTERPOLATION_METHOD: WRONG ELEMENT')
                if( present(mask) ) if_mask = mask(ielem)
                if( if_mask ) then
                   pelty = mesh % ltype(ielem)
                   pnode = element_type(pelty) % number_nodes
                   do inode = 1,pnode
                      elcod(1:ndime,inode) = mesh % coord(1:ndime,mesh % lnods(inode,ielem))
                   end do
                   call elmgeo_natural_coordinates(&
                        ndime,pelty,pnode,elcod,&
                        shapf(:,ii),derit,xx(:,ii),coloc,&
                        ifoun,toler_loc)
                   if( ifoun /= 0 ) then
                      lelem(ii) = ielem
                      if( present(lenty) ) then
                         lenty(1:pnode,ii) = mesh % lnods(1:pnode,ielem)
                      end if
                      if( present(dista) ) dista(ii) = 0.0_rp
                      if( present(deriv) ) then
                         call elmgeo_cartesian_derivatives(ndime,pnode,elcod,derit,deriv(:,:,ii))
                      end if
                      exit loop_kk
                   end if
                end if
             end do loop_kk
             if(self % input_data % force_find.and.(ifoun == 0)) then
               toler_loc_force = toler_loc
               iter_force = 0_ip
               do while(ifoun==0)
                 iter_force      = iter_force + 1_ip
                 toler_loc_force = toler_loc_force*10_rp
                 if(kfl_paral==-1_ip) then
                    ! SEQUENTIAL CHECK
                   if(iter_force > self % input_data % max_it_force) call runend('Forced find not succeeded: required more iterations than allowed')
                 end if
                 loop_kk2: do kk = 1,nelem
                    ielem = list_entities(ii) % l(kk)
                    if( ielem <= 0 ) call runend('DEF_INTERPOLATION_METHOD: WRONG ELEMENT')
                     if( present(mask) ) if_mask = mask(ielem)
                     if( if_mask ) then
                        pelty = mesh % ltype(ielem)
                        pnode = element_type(pelty) % number_nodes
                        do inode = 1,pnode
                           elcod(1:ndime,inode) = mesh % coord(1:ndime,mesh % lnods(inode,ielem))
                        end do
                        call elmgeo_natural_coordinates(&
                             ndime,pelty,pnode,elcod,&
                             shapf(:,ii),derit,xx(:,ii),coloc,&
                             ifoun,toler_loc_force)
                        if( ifoun /= 0 ) then
                           lelem(ii) = ielem
                           if( present(lenty) ) then
                              lenty(1:pnode,ii) = mesh % lnods(1:pnode,ielem)
                           end if
                           if( present(dista) ) then
                             !dista(ii) = 0.0_rp
                             dista(ii) = toler_loc_force ! i roughly use as distance the relaxed tolerance
                             ! This way, if it is found in two different parallel domains, the one where it has been found
                             ! with a finer tolerance is the one chosen :-)
                             ! TODO perfect: Abel -> round coloc to valid projecting to faces (if outside),and compute distance
                             !               easy in simplices and hexes
                             ! Should find if inside using coloc. Otherwise find how much outside it is.
                             ! This is just to solve any conflict if it is found in two different domain in parallel
                             ! while relaxing the tolerance... hardly difficult for the same relax tolerance..
                           end if
                           if( present(deriv) ) then
                              call elmgeo_cartesian_derivatives(ndime,pnode,elcod,derit,deriv(:,:,ii))
                           end if
                           exit loop_kk2
                        end if
                     end if
                   end do loop_kk2
                 end do !while
             end if
          end if
       end do

       !write(6,*)'c=',kfl_paral,nn,nelem ; flush(6)

       self % stats(2) = self % stats(2) / real(nn,rp)
       call cputim(time2) ; self % times(2) = self % times(2) + time2-time1 ; time1 = time2
       call memory_deallo(memor_loc,'LIST_ENTITIES',vacal,list_entities)       
    end if

    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)

  end subroutine interpolation_element

  subroutine interpolation_boundary(self,search,xx,mesh,lelem,shapf,deriv,dista,lenty,mask,TOLER,MEMORY_COUNTER)

    class(interpolation),                    intent(inout) :: self
    class(search_method),                    intent(inout) :: search
    real(rp),                       pointer, intent(in)    :: xx(:,:)
    class(mesh_type_basic),                  intent(in)    :: mesh
    integer(ip),                    pointer, intent(inout) :: lelem(:)
    real(rp),                       pointer, intent(inout) :: shapf(:,:)
    real(rp),         optional,     pointer, intent(inout) :: deriv(:,:,:)
    real(rp),         optional,     pointer, intent(inout) :: dista(:)
    integer(ip),      optional,     pointer, intent(inout) :: lenty(:,:)
    logical(lg),      optional,     pointer, intent(in)    :: mask(:)
    real(rp),         optional,              intent(in)    :: TOLER
    integer(8),       optional,              intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(ip)                                            :: ndime,nn,ii,kk
    integer(ip)                                            :: nboun,mnode
    integer(ip)                                            :: inodb,iboun,pelty,mnodb
    integer(ip)                                            :: pblty,pnodb,pnode,ielem
    integer(ip)                                            :: ifoun,mnodf,ipoin,inode
    type(i1p),                      pointer                :: list_entities(:)
    integer(8)                                             :: memor_loc(2)
    real(rp)                                               :: coloc(3),proje(3)
    real(rp)                                               :: pnear(3),time1,time2
    real(rp)                                               :: baloc(mesh %ndime,mesh % ndime)
    real(rp)                                               :: xjaci(3,3),gpdet,my_sign
    real(rp)                                               :: toler_loc,dista_loc,dista_min
    real(rp),                       allocatable            :: shapt(:)
    real(rp),                       allocatable            :: derit(:,:)
    real(rp),                       allocatable            :: elcod(:,:)
    real(rp),                       allocatable            :: bocod(:,:)
    integer(ip),                    pointer                :: lnodb(:,:)
    integer(ip),                    pointer                :: ltypb(:)
    real(rp),                       pointer                :: coord(:,:)
    real(rp),                       pointer                :: coorv(:,:)
    integer(ip),                    pointer                :: lnods(:,:)
    integer(ip),                    pointer                :: ltype(:)
    integer(ip),                    pointer                :: lelbo(:)
    logical(lg)                                            :: if_mask

    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)
    toler_loc = optional_argument(epsil,TOLER)
    if_mask   = .true.
    nn        = memory_size(xx,2_ip)
    ndime     = self % nd
    mnode     = -1
    mnodb     = -1
    nullify(list_entities)

    if( nn > 0 ) then
       !
       ! Boundary mesh pointers
       !
       select type ( mesh )
       class is ( bmsh_type_basic )
          mnodb =  mesh % mnode
          ltypb => mesh % ltype
          lnodb => mesh % lnods
          coord => mesh % coord
       class is ( mesh_type_basic )
          if( associated(mesh % boundary) ) then
             mnodb =  mesh % boundary % mnode
             ltypb => mesh % boundary % ltype
             lnodb => mesh % boundary % lnods
             coord => mesh % coord
          else
             call runend('INTERPOLATION_BOUNDARY: NO BOUNDARY MESH ASSOCIATED TO VOLUME MESH')          
          end if
       end select
       !
       ! Volume mesh pointers if derivatives are required
       !
       mnode = 0_ip
       if( present(deriv) .and. self % input_data % deriv ) then
          select type ( mesh )
          class is ( bmsh_type_basic )
             if( associated(mesh % mesh) ) then
                mnode =  mesh % mesh % mnode
                ltype => mesh % mesh % ltype
                lnods => mesh % mesh % lnods
                coorv => mesh % mesh % coord
                lelbo => mesh % lelbo                
             else
                call runend('INTERPOLATION_BOUNDARY: NO VOLUME MESH ASSOCIATED TO BOUNDARY MESH')
             end if
          class is ( mesh_type_basic )
             if( associated(mesh % boundary) ) then
                mnode =  mesh % mnode
                ltype => mesh % ltype
                lnods => mesh % lnods
                coorv => mesh % coord
                lelbo => mesh % boundary % lelbo
             end if
          end select
       end if
       mnodf = max(mnode,mnodb)   
       !
       ! Allocate local
       !
       allocate(elcod(ndime,mnode))
       allocate(bocod(ndime,mnodb))
       allocate(derit(ndime,mnodf))
       allocate(shapt(mnodf))
       !
       ! Allocate if necessary
       !
       call interpolation_allocate(self,nn,ndime,max(mnode,mnodb),lelem,shapf,deriv,dista,lenty,MEMORY_COUNTER=memor_loc)
       
       call cputim(time1)

#ifdef __PGI
       call search % candidate(XX=xx,LIST_ENTITIES=list_entities,METHOD=CANDIDATE_NEAREST,MASK=std_log_1,MEMORY_COUNTER=memor_loc)       
#else
       call search % candidate(xx,list_entities,METHOD=CANDIDATE_NEAREST,MEMORY_COUNTER=memor_loc)       
#endif
       call cputim(time2) ; self % times(1) = self % times(1) + time2-time1 ; time1 = time2
       !
       ! Check boundary mesh
       !
       if( associated(list_entities) ) then
          do ii = 1,nn
             nboun = memory_size(list_entities(ii)%l)
             if( present(dista) ) dista(ii) = -huge(1.0_rp)
             if( nboun == 0 ) then
                lelem(ii)   =  0
                shapf(:,ii) =  0.0_rp
                if( present(deriv) .and. self % input_data % deriv ) deriv(:,:,ii) = 0.0_rp 
             else
                self % stats(2) = self % stats(2) + real(nboun,rp)
                dista_min = huge(1.0_rp)
                loop_kk: do kk = 1,nboun
                   iboun = list_entities(ii) % l(kk)
                   if( present(mask) ) if_mask = mask(iboun)
                   if( iboun > 0 .and. if_mask ) then
                      pblty = ltypb(iboun)                
                      pnodb = element_type(pblty) % number_nodes
                      do inodb = 1,pnodb
                         bocod(1:ndime,inodb) = coord(1:ndime,lnodb(inodb,iboun))
                      end do
                      !
                      ! Compute projection on boundary
                      !
                      call elmgeo_projection_on_a_face(&
                           ndime,pblty,bocod,xx(:,ii),proje)
                      !
                      ! Check of point is on boundary
                      !
                      call elmgeo_natural_coordinates_on_boundaries(&
                           ndime,pblty,pnodb,bocod, &
                           shapt,derit,proje, & 
                           coloc,ifoun,toler_loc,NEAREST_POINT=pnear)

                      call elmgeo_jacobian_boundary(&
                           ndime,pnodb,bocod,derit,gpdet,xjaci,baloc)
                      
                      if( ifoun /= 0 ) then
                         proje(1:ndime) = pnear(1:ndime)
                         my_sign        = dot_product(baloc(1:ndime,ndime),xx(1:ndime,ii)-proje(1:ndime))
                         dista_loc      = sqrt(dot_product(proje(1:ndime)-xx(1:ndime,ii),proje(1:ndime)-xx(1:ndime,ii)))
                         if( dista_loc < dista_min ) then
                            lelem(ii)         = iboun
                            shapf(1:pnodb,ii) = shapt(1:pnodb)
                            dista_min         = dista_loc
                            if( present(lenty) ) then
                               lenty(1:pnodb,ii) = lnodb(1:pnodb,iboun)
                            end if
                            if( present(dista) ) then
                               if( my_sign >= 0.0_rp ) then
                                  dista(ii) =  dista_loc
                               else
                                  dista(ii) = -dista_loc                                  
                               end if
                            end if
                            !
                            ! Derivatives not available
                            !
                            if( present(deriv) .and. self % input_data % deriv ) then                       
                               ielem = lelbo(iboun)
                               pelty = ltype(ielem)
                               pnode = element_type(pelty) % number_nodes
                               do inode = 1,pnode
                                  ipoin                = lnods(inode,ielem)
                                  elcod(1:ndime,inode) = coorv(1:ndime,ipoin)
                               end do
                               call elmgeo_natural_coordinates(&
                                    ndime,pelty,pnode,elcod,&
                                    shapt,derit,xx(:,ii),coloc,&
                                    ifoun,toler_loc)
                               if( ifoun /= 0 ) then
                                  call elmgeo_cartesian_derivatives(ndime,pnode,elcod,derit,deriv(:,:,ii))
                               end if
                            end if
                         end if
                      end if
                   end if
                end do loop_kk
             end if
          end do
       end if

       self % stats(2) = self % stats(2)/ real(nn,rp)
       call cputim(time2) ; self % times(2) = self % times(2) + time2-time1 ; time1 = time2

       if( allocated(shapt) ) deallocate(shapt)
       if( allocated(elcod) ) deallocate(elcod)
       if( allocated(bocod) ) deallocate(bocod)
       if( allocated(derit) ) deallocate(derit)
       call memory_deallo(memor_loc,'LIST_ENTITIES',vacal,list_entities)       
    end if

    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)

  end subroutine interpolation_boundary

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-05
  !> @brief   Nearest node
  !> @details Nearest node
  !> 
  !-----------------------------------------------------------------------
  
  subroutine nearest_node(self,search,xx,mesh,lelem,shapf,deriv,dista,lenty,mask,TOLER,MEMORY_COUNTER)

    class(interpolation),                    intent(inout) :: self
    class(search_method),                    intent(inout) :: search
    real(rp),                       pointer, intent(in)    :: xx(:,:)
    class(mesh_type_basic),                  intent(in)    :: mesh
    integer(ip),                    pointer, intent(inout) :: lelem(:)
    real(rp),                       pointer, intent(inout) :: shapf(:,:)
    real(rp),         optional,     pointer, intent(inout) :: deriv(:,:,:)
    real(rp),         optional,     pointer, intent(inout) :: dista(:)
    integer(ip),      optional,     pointer, intent(inout) :: lenty(:,:)
    logical(lg),      optional,     pointer, intent(in)    :: mask(:)
    real(rp),         optional,              intent(in)    :: TOLER
    integer(8),       optional,              intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(ip)                                            :: ndime,nn,ii,kk,ipoin
    integer(ip)                                            :: mnode,npoin
    type(i1p),                      pointer                :: list_entities(:)
    integer(8)                                             :: memor_loc(2)
    real(rp)                                               :: toler_loc
    real(rp)                                               :: mydis,dimin,time1,time2
    logical(lg)                                            :: if_mask

    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)

    toler_loc = optional_argument(epsil,TOLER)
    if_mask   = .true.
    nn        = memory_size(xx,2_ip)
    mnode     = mesh % mnode
    ndime     = self % nd 
    nullify(list_entities)
    
    if( nn > 0 ) then
       !
       ! Allocate if necessary
       !
       call interpolation_allocate(self,nn,ndime,mnode,lelem,shapf,deriv,dista,lenty,MEMORY_COUNTER=memor_loc)

       call cputim(time1)
#ifdef __PGI
       call search % candidate(XX=xx,LIST_ENTITIES=list_entities,METHOD=CANDIDATE_NEAREST,MASK=std_log_1,MEMORY_COUNTER=memor_loc)       
#else
       call search % candidate(xx,list_entities,METHOD=CANDIDATE_NEAREST,MEMORY_COUNTER=memor_loc)       
#endif
       call cputim(time2)
       self % times(1) = self % times(1) + time2-time1 ; time1 = time2
       
       if( associated(list_entities) ) then
          do ii = 1,nn
             dimin = huge(1.0_rp)
             npoin = memory_size(list_entities(ii)%l)
             self % stats(2) = self % stats(2) + real(npoin,rp)
             if( present(dista) ) dista(ii) = -huge(1.0_rp)
             do kk = 1,npoin
                ipoin = list_entities(ii) % l(kk)
                if( ipoin > 0 ) then
                   if( present(mask) ) if_mask = mask(ipoin)
                   if( if_mask ) then
                      mydis = sqrt(dot_product(mesh % coord(1:ndime,ipoin)-xx(1:ndime,ii),mesh % coord(1:ndime,ipoin)-xx(1:ndime,ii)))
                      if( mydis < dimin ) then
                         dimin       = mydis
                         lelem(ii)   = ipoin
                         shapf(1,ii) = 1.0_rp
                         if( present(deriv) ) deriv(:,:,ii) = 0.0_rp
                         if( present(dista) ) dista(ii)     = dimin
                         if( present(lenty) ) lenty(1,ii)   = ipoin
                      end if
                   end if
                end if
             end do
          end do
       end if

       self % stats(2) = self % stats(2) / real(nn,rp)
       call cputim(time2) ; self % times(2) = self % times(2) + time2-time1 ; time1 = time2

    end if

    call memory_deallo(memor_loc,'LIST_ENTITIES',vacal,list_entities)       
    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)

  end subroutine nearest_node

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-05
  !> @brief   Nearest node
  !> @details Nearest node
  !> 
  !-----------------------------------------------------------------------
  
  subroutine global_numbering(self,search,xx,ll,ht,mesh,lelem,shapf,deriv,lenty,mask,TOLER,MEMORY_COUNTER)

    class(interpolation),                    intent(inout) :: self
    class(search_method),                    intent(inout) :: search
    real(rp),                       pointer, intent(in)    :: xx(:,:)
    integer(ip),                    pointer, intent(in)    :: ll(:)
    type(hash_t),                            intent(in)    :: ht
    class(mesh_type_basic),                  intent(in)    :: mesh
    integer(ip),                    pointer, intent(inout) :: lelem(:)
    real(rp),                       pointer, intent(inout) :: shapf(:,:)
    real(rp),         optional,     pointer, intent(inout) :: deriv(:,:,:)
    integer(ip),      optional,     pointer, intent(inout) :: lenty(:,:)
    logical(lg),      optional,     pointer, intent(in)    :: mask(:)
    real(rp),         optional,              intent(in)    :: TOLER
    integer(8),       optional,              intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(ip)                                            :: ndime,nn,ii
    integer(ip)                                            :: mnode
    integer(ip)                                            :: ipoin_local,ipoin_global
    type(i1p),                      pointer                :: list_entities(:)
    integer(8)                                             :: memor_loc(2)
    real(rp)                                               :: toler_loc
    real(rp)                                               :: time1,time2
    logical(lg)                                            :: if_mask

    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)

    toler_loc = optional_argument(epsil,TOLER)
    if_mask   = .true.
    nn        = memory_size(xx,2_ip)
    mnode     = mesh % mnode
    ndime     = self % nd 
    nullify(list_entities)
    
    if( nn > 0 ) then
       !
       ! Allocate if necessary
       !
       call cputim(time1)
       call interpolation_allocate(self,nn,ndime,mnode,LELEM=lelem,SHAPF=shapf,DERIV=deriv,LENTY=lenty,MEMORY_COUNTER=memor_loc)

       do ii = 1,nn
          ipoin_global = abs(ll(ii))
          ipoin_local  = htalid(ht,ipoin_global)
          !block ; use def_master ;  if(ipoin_global==26) write(6,*)'AQUI=',kfl_paral,ipoin_local ;end block
          if( ipoin_local > 0 ) then
             self % stats(2) = self % stats(2) + 1.0_rp
             lelem(ii)   = ipoin_local
             shapf(1,ii) = 1.0_rp
             if( present(deriv) ) deriv(:,:,ii) = 0.0_rp
             if( present(lenty) ) lenty(1,ii)   = ipoin_local
          else                  
             lelem(ii) = 0
          end if
       end do
       
       self % stats(2) = self % stats(2) / real(nn,rp)
       call cputim(time2) ; self % stats(5) = self % stats(5) + time2-time1 ; time1 = time2

    end if

    call memory_deallo(memor_loc,'LIST_ENTITIES',vacal,list_entities)       
    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)

  end subroutine global_numbering

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-02-02
  !> @brief   Allocate
  !> @details Allocate search arrays
  !> 
  !-----------------------------------------------------------------------
  
  subroutine interpolation_allocate(self,nn,ndime,mnode,lelem,shapf,deriv,dista,lenty,MEMORY_COUNTER)
  
    class(interpolation),                    intent(inout) :: self
    integer(ip),                             intent(in)    :: nn
    integer(ip),                             intent(in)    :: ndime
    integer(ip),                             intent(in)    :: mnode
    integer(ip),                    pointer, intent(inout) :: lelem(:)
    real(rp),                       pointer, intent(inout) :: shapf(:,:)
    real(rp),         optional,     pointer, intent(inout) :: deriv(:,:,:)
    real(rp),         optional,     pointer, intent(inout) :: dista(:)
    integer(ip),      optional,     pointer, intent(inout) :: lenty(:,:)
    integer(8),       optional,              intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(8)                                             :: memor_loc(2)
    
    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)

    if( .not. associated(lelem)    ) call memory_alloca(memor_loc,'LELEM',vacal,lelem,nn)
    if( .not. associated(shapf)    ) call memory_alloca(memor_loc,'SHAPF',vacal,shapf,mnode,nn)
    if( present(deriv) ) then
       if( .not. associated(deriv) ) call memory_alloca(memor_loc,'DERIV',vacal,deriv,ndime,mnode,nn)
    end if
    if( present(dista) ) then
       if( .not. associated(dista) ) call memory_alloca(memor_loc,'DISTA',vacal,dista,nn)
    end if
    if( present(lenty) ) then
       if( .not. associated(lenty) ) call memory_alloca(memor_loc,'LENTY',vacal,lenty,mnode,nn)
    end if
       
    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)

  end subroutine interpolation_allocate
  
end module def_interpolation_method
!> @}
