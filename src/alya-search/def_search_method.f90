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
!>                                  search_method
!>                                       ||
!>                                search_method_par
!>                  ||            ||            ||            ||
!>               maths_bin   maths_octree   maths_octbin   typ_kdtree
!>
!>          Deferred procedures:
!>          --------------------
!>          init() ............ Initialize the search strategy
!>          deallo() .......... Deallocate the search strategy structure
!>          candidate() ....... Give a list of entity candidate when doing
!>                              a point search
!>          fill() ............ Fill int eh search structure?
!>          mesh() ............ Create a basic mesh
!>          results() ......... Get some results on mesh
!>          self % times(1) ... Time for fill
!>          self % times(2) ... Time for search
!>          self % stats(:) ... Statistics
!>
!-----------------------------------------------------------------------

module def_search_method

  use def_kintyp_basic,      only : ip,rp,lg,i1p,r1p,r2p
  use mod_memory_basic,      only : memory_alloca
  use mod_memory_basic,      only : memory_deallo
  use mod_memory_basic,      only : memory_size
  use mod_optional_argument, only : optional_argument
  use mod_std

  implicit none
  private
  
  real(rp),      parameter :: epsil                 = epsilon(1.0_rp)
  real(rp),      parameter :: big                   = 1.0e16_rp
  character(17), parameter :: vacal                 = 'def_search_method'
  integer(ip),   parameter :: SEARCH_POINTS         = 0
  integer(ip),   parameter :: SEARCH_BOUNDING_BOXES = 1
  integer(ip),   parameter :: SEARCH_BIN            = 1
  integer(ip),   parameter :: SEARCH_OCTBIN         = 2
  integer(ip),   parameter :: SEARCH_OCTREE         = 3
  integer(ip),   parameter :: SEARCH_KDTREE         = 4
  integer(ip),   parameter :: SEARCH_SKDTREE        = 5

  integer(ip),   parameter :: SEARCH_FILL           = 1
  integer(ip),   parameter :: SEARCH_CANDIDATE      = 2
  integer(ip),   parameter :: SEARCH_DEALLO         = 3

  integer(ip),   parameter :: CANDIDATE_INSIDE      = 0
  integer(ip),   parameter :: CANDIDATE_NEAREST     = 1

  integer(ip),   parameter :: SEARCH_MESH_LEAVES    = 0
  integer(ip),   parameter :: SEARCH_MESH_FILL      = 1
  integer(ip),   parameter :: SEARCH_MESH_ALL       = 2

  real(rp)                 :: time1,time2
  
  !----------------------------------------------------------------------
  ! 
  ! Sequential search methods
  !
  !----------------------------------------------------------------------
  
  type, abstract :: search_method
     character(LEN=:), allocatable                       :: name                  ! Name of the method
     integer(ip)                                         :: mode                  ! 0: sequential, 1: parallel
     integer(ip)                                         :: dim                   ! Space dimensions
     integer(ip)                                         :: nelem                 ! Number of elements
     integer(ip)                                         :: kfl_filled            ! If method has been filled
     real(rp)                                            :: comin(3)              ! Minimum coordinates
     real(rp)                                            :: comax(3)              ! Maximum coordinates
     real(rp)                                            :: toler_rel             ! Relative tolerance
     real(rp)                                            :: toler_abs             ! Absolute tolerance
     real(rp)                                            :: times(10)             ! Timings
     real(rp)                                            :: stats(10)             ! Statistics
     integer(8)                                          :: memor(2)              ! Memory counter
     integer(ip)                                         :: fill_method           ! Point or bounding boxes
   contains
     procedure (search_method_init),      pass, deferred :: init                  ! Initialize
     procedure (search_method_deallo),    pass, deferred :: deallo                ! Deallocate
     procedure (search_method_candidate), pass, deferred :: candidate             ! List of candidates
     procedure (search_method_fill),      pass, deferred :: fill                  ! Fill the structures
     procedure (search_method_mesh),      pass, deferred :: mesh                  ! Crate a mesh
     procedure (search_method_results),   pass, deferred :: results               ! Get results on mesh
     procedure (search_method_type_name), pass, deferred :: type_name             ! Type name of the search method
     
     procedure,                           pass           :: bounding_box_coord    ! Bounding box of geometry
     procedure,                           pass           :: bounding_box_bobox    ! Bounding box of geometry
     generic                                             :: bounding_box => &
          &                                                 bounding_box_coord,&
          &                                                 bounding_box_bobox
     procedure,                           pass           :: init_all              ! Common initializations
     procedure,                           pass           :: input_all             ! Input data
     procedure,                           pass           :: mem_ini               ! Init memory counters
     procedure,                           pass           :: mem_end               ! End memory counters
     procedure,                           pass           :: tim_ini               ! Init time counters
     procedure,                           pass           :: tim_end               ! End time counters
  end type search_method
  
  abstract interface
          
     subroutine search_method_init(self)
       import                                                 :: search_method
       class(search_method),                    intent(inout) :: self
     end subroutine search_method_init
          
     subroutine search_method_deallo(self,list_entities,MEMORY_COUNTER)
       import                                                 :: search_method,i1p
       class(search_method),                    intent(inout) :: self
       type(i1p),         optional, pointer,    intent(inout) :: list_entities(:)    
       integer(8),        optional,             intent(inout) :: MEMORY_COUNTER(2)
     end subroutine search_method_deallo

     subroutine search_method_candidate(self,xx,list_entities,METHOD,MASK,MEMORY_COUNTER)
       import                                                 :: search_method,ip,rp,i1p,lg
       class(search_method),                    intent(inout) :: self
       real(rp),                                intent(in)    :: xx(:,:)               !< List coordinates
       type(i1p),                      pointer, intent(inout) :: list_entities(:)      !< List of entities found
       integer(ip),      optional,              intent(in)    :: METHOD                !< Method for candidates
       logical(lg),      optional,     pointer, intent(in)    :: MASK(:)               !< Mask to consider or not points
       integer(8),       optional,              intent(inout) :: MEMORY_COUNTER(2)     !< Memory counters
     end subroutine search_method_candidate
     
     subroutine search_method_fill(self,coord,bobox,MEMORY_COUNTER,COORD_MIN,COORD_MAX,PERMU,MASK)
       import                                                 :: search_method,ip,rp,lg
       class(search_method),                    intent(inout) :: self            
       real(rp),    optional,          pointer, intent(in)    :: coord(:,:)            !< Coordinates
       real(rp),    optional,          pointer, intent(in)    :: bobox(:,:,:)          !< Bounding boxes
       integer(8),  optional,                   intent(inout) :: MEMORY_COUNTER(2)     !< Memory counters
       real(rp),    optional,                   intent(in)    :: COORD_MIN(:)          !< Min coordinates
       real(rp),    optional,                   intent(in)    :: COORD_MAX(:)          !< Max coordinates       
       integer(ip), optional,          pointer, intent(in)    :: PERMU(:)              !< Permutation
       logical(lg), optional,          pointer, intent(in)    :: MASK(:)               !< mask
     end subroutine search_method_fill

     subroutine search_method_mesh(self,ndime,mnode,nelem,npoin,lnods,ltype,coord,MEMORY_COUNTER,&
          OFFSET_IELEM,OFFSET_IPOIN,CRITERION,CENTROID,ONLY_DEALLOCATE)
       import                                                 :: search_method,ip,rp,lg       
       class(search_method),                   intent(inout)  :: self
       integer(ip),                            intent(out)    :: ndime
       integer(ip),                            intent(out)    :: mnode
       integer(ip),                            intent(out)    :: nelem
       integer(ip),                            intent(out)    :: npoin
       integer(ip),                   pointer, intent(inout)  :: lnods(:,:)
       integer(ip),                   pointer, intent(inout)  :: ltype(:)
       real(rp),                      pointer, intent(inout)  :: coord(:,:)
       integer(8),          optional,          intent(inout)  :: MEMORY_COUNTER(2)
       integer(ip),         optional,          intent(in)     :: OFFSET_IELEM
       integer(ip),         optional,          intent(in)     :: OFFSET_IPOIN
       integer(ip),         optional,          intent(in)     :: CRITERION     !< If empty bins should be considered
       real(rp),            optional, pointer, intent(inout)  :: CENTROID(:,:)
       logical(lg),         optional,          intent(in)     :: ONLY_DEALLOCATE
     end subroutine search_method_mesh

     subroutine search_method_results(self,xx,names,OFFSET,MEMORY_COUNTER,CRITERION,ONLY_DEALLOCATE)
       import                                                 :: search_method,ip,rp,lg       
       class(search_method),                    intent(inout) :: self
       real(rp),                       pointer, intent(inout) :: xx(:,:)
       character(len=5),               pointer, intent(inout) :: names(:)
       integer(ip), optional,                   intent(in)    :: OFFSET
       integer(8),  optional,                   intent(inout) :: MEMORY_COUNTER(2)
       integer(ip), optional,                   intent(in)    :: CRITERION     
       logical(lg), optional,                   intent(in)    :: ONLY_DEALLOCATE       
     end subroutine search_method_results

     function search_method_type_name(self) result(name)
       import                                                 :: search_method,ip,rp,lg       
       class(search_method),                    intent(inout) :: self
       character(LEN=:), allocatable                          :: name
     end function search_method_type_name
     
  end interface
  
  public :: big

  public :: search_method
  public :: search_method_read_data
  public :: SEARCH_POINTS
  public :: SEARCH_BOUNDING_BOXES
  public :: SEARCH_BIN 
  public :: SEARCH_OCTBIN 
  public :: SEARCH_OCTREE 
  public :: SEARCH_KDTREE 
  public :: SEARCH_SKDTREE 

  public :: SEARCH_FILL 
  public :: SEARCH_CANDIDATE 
  public :: SEARCH_DEALLO 

  public :: CANDIDATE_INSIDE   
  public :: CANDIDATE_NEAREST  

  public :: SEARCH_MESH_LEAVES 
  public :: SEARCH_MESH_FILL 
  public :: SEARCH_MESH_ALL
  
  !----------------------------------------------------------------------
  !
  ! Parallel search methods
  !
  !----------------------------------------------------------------------
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Start timing
  !> @details End timing
  !> 
  !-----------------------------------------------------------------------

  subroutine tim_ini(self)

    class(search_method), intent(in)  :: self
    
    call cputim(time1)

  end subroutine tim_ini
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Start timing
  !> @details End timing
  !> 
  !-----------------------------------------------------------------------

  subroutine tim_end(self,ii)

    class(search_method), intent(inout) :: self
    integer(ip),          intent(in)    :: ii
    
    call cputim(time2)
    self % times(ii) = self % times(ii) + (time2 - time1)
    
  end subroutine tim_end
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Start memory
  !> @details End memory
  !> 
  !-----------------------------------------------------------------------

  subroutine mem_ini(self,memor,MEMORY_COUNTER)

    class(search_method),           intent(in)  :: self
    integer(8),                     intent(out) :: memor(2)
    integer(8),           optional, intent(in)  :: MEMORY_COUNTER(2)
    
    memor = optional_argument(self % memor,MEMORY_COUNTER)

  end subroutine mem_ini
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Start memory
  !> @details End memory
  !> 
  !-----------------------------------------------------------------------

  subroutine mem_end(self,memor,MEMORY_COUNTER)

    class(search_method),           intent(inout) :: self
    integer(8),                     intent(in)    :: memor(2)
    integer(8),           optional, intent(inout) :: MEMORY_COUNTER(2)
    
    if( present(MEMORY_COUNTER) ) then
       MEMORY_COUNTER = memor
    else
       self % memor = memor
    end if

  end subroutine mem_end
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Input
  !> @details Set input data
  !> 
  !-----------------------------------------------------------------------

  subroutine input_all(self,relative_tolerance,absolute_tolerance,fill_method,name)
    
    class(search_method),       intent(inout) :: self
    real(rp),         optional, intent(in)    :: relative_tolerance
    real(rp),         optional, intent(in)    :: absolute_tolerance
    integer(ip),      optional, intent(in)    :: fill_method
    character(len=*), optional, intent(in)    :: name

    if( present(relative_tolerance) ) self % toler_rel   = relative_tolerance
    if( present(absolute_tolerance) ) self % toler_abs   = absolute_tolerance
    if( present(fill_method)        ) self % fill_method = fill_method
    if( present(name)               ) self % name        = name
    
  end subroutine input_all

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-03
  !> @brief   Read data
  !> @details Read data
  !> 
  !-----------------------------------------------------------------------

  subroutine search_method_read_data(words,param,kfl_search_method,search_method_param)

    character(len=5),     intent(in)  :: words(:)
    real(rp),             intent(in)  :: param(:)
    integer(ip),          intent(out) :: kfl_search_method
    real(rp),             intent(out) :: search_method_param(:)
    
    select case( words(2) )

    case ( 'BIN  ' ) 
       
       kfl_search_method = SEARCH_BIN 
       if( words(3) == 'BOXES' ) search_method_param(1) = param(3)
       search_method_param(2:3) = param(3)
       
    case ('OCTTR' , 'OCT  ' , 'OCTRE' )                

       kfl_search_method = SEARCH_OCTREE
       if( words(3) == 'MAXIM' ) search_method_param(1) = param(3)
       if( words(3) == 'LIMIT' ) search_method_param(1) = param(3)

    end select

  end subroutine search_method_read_data
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-03
  !> @brief   Initialize common variables of the class
  !> @details Initialization  
  !> 
  !-----------------------------------------------------------------------

  subroutine init_all(self)
    
    class(search_method), intent(inout) :: self
    !
    ! Input parameters
    !
    self % mode          =  0
    self % toler_rel     =  epsil
    self % toler_abs     =  0.0_rp
    self % fill_method   =  SEARCH_POINTS
    self % kfl_filled    =  0_ip    
    !
    ! Parameters
    !
    self % dim           =  0
    self % nelem         =  0
    self % comin         =  big
    self % comax         = -big
    self % memor         =  0_8
    !
    ! Output
    !
    self % times         =  0.0_rp
    self % stats         =  0.0_rp

  end subroutine init_all
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-03
  !> @brief   Bounding box
  !> @details Bounding box of bounding boxes
  !> 
  !-----------------------------------------------------------------------

  subroutine bounding_box_coord(self,coord,OFFSET,COORD_MIN,COORD_MAX)
    
    class(search_method),           intent(inout) :: self
    real(rp),                       intent(in)    :: coord(:,:)
    real(rp),             optional, intent(in)    :: OFFSET 
    real(rp),             optional, intent(in)    :: COORD_MIN(:) 
    real(rp),             optional, intent(in)    :: COORD_MAX(:)
    real(rp)                                      :: delta(3),offset_loc
    integer(ip)                                   :: ndime,idime

    if( present(COORD_MIN) .and. present(COORD_MAX) ) then
       
       self % comin(1:self % dim) = COORD_MIN(1:self % dim)
       self % comax(1:self % dim) = COORD_MAX(1:self % dim)
       
    else

       if( size(coord) > 0 ) then

          offset_loc = optional_argument(epsil,OFFSET)
          if( self % dim /= 0 ) then
             ndime = self % dim
          else
             ndime = int(size(coord,1),ip)
          end if
          do idime = 1,self % dim
             self % comin(idime) = minval(coord(idime,:))
             self % comax(idime) = maxval(coord(idime,:))
          end do
          delta        = ( self % comax - self % comin ) *  offset_loc  
          self % comin = self % comin - delta - epsil
          self % comax = self % comax + delta + epsil
       end if
       
    end if

  end subroutine bounding_box_coord

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-03
  !> @brief   Bounding box
  !> @details Boundig box of the geometry
  !> 
  !-----------------------------------------------------------------------

  subroutine bounding_box_bobox(self,bobox,OFFSET,COORD_MIN,COORD_MAX)
    
    class(search_method),           intent(inout) :: self
    real(rp),                       intent(in)    :: bobox(:,:,:)
    real(rp),             optional, intent(in)    :: OFFSET 
    real(rp),             optional, intent(in)    :: COORD_MIN(:) 
    real(rp),             optional, intent(in)    :: COORD_MAX(:)
    real(rp)                                      :: delta(3),offset_loc
    integer(ip)                                   :: ndime,idime

    if( present(COORD_MIN) .and. present(COORD_MAX) ) then
       
       self % comin(1:self % dim) = COORD_MIN(1:self % dim)
       self % comax(1:self % dim) = COORD_MAX(1:self % dim)
       
    else

       if( size(bobox) > 0 ) then

          offset_loc = optional_argument(epsil,OFFSET)
          if( self % dim /= 0 ) then
             ndime = self % dim
          else
             ndime = int(size(bobox,2),ip)
          end if
          do idime = 1,self % dim
             self % comin(idime) = minval(bobox(1,idime,:))
             self % comax(idime) = maxval(bobox(2,idime,:))
          end do
          delta        = ( self % comax - self % comin ) *  offset_loc  
          self % comin = self % comin - delta - epsil
          self % comax = self % comax + delta + epsil
       end if
       
    end if

  end subroutine bounding_box_bobox

end module def_search_method
!> @}
