!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @defgroup Maths_Toolbox 
!> Toolbox for bins
!> @{
!> @name    ToolBox for bins
!> @file    def_maths.f90
!> @author  Guillaume Houzeaux
!> @brief   Tree
!> @details Tree based methods 
!>          maths_octree: OCTREE .... always split into 4 or 8 
!>          maths_kdtree: SKDTREE ... split in 2, in x,y,z alternatively
!>
!>          Input
!>          -----
!>          adjust_bb ......... .true.: Adjust boundaing boxes to actual
!>                              bounding box
!>          divmax ............ Maximum number of points in boxes
!>
!>          fill
!>          ----
!>          self % stats(1) ... Total number of leaves
!>          self % stats(2) ... Percentage of addition elements
!>          self % stats(3) ... Maximum depth (level)
!>          self % stats(4) ... Max element in a box
!>          self % times(1) ... Time for fill
!>
!>          candidate
!>          ---------
!>          self % stats(5) ... Saving ratio: average number of candidates
!>                              divided by total number of elements
!>          self % times(2) ... Time for search
!
!-----------------------------------------------------------------------

module def_maths_tree

  use def_kintyp_basic,      only : ip,rp,lg,i1p
  use def_elmtyp,            only : QUA04
  use def_elmtyp,            only : HEX08
  use mod_memory_basic,      only : memory_alloca
  use mod_memory_basic,      only : memory_deallo
  use mod_memory_basic,      only : memory_size
  use mod_memory_basic,      only : memory_resize
  use def_search_method,     only : search_method   
  use def_search_method,     only : SEARCH_FILL
  use def_search_method,     only : SEARCH_CANDIDATE
  use def_search_method,     only : SEARCH_DEALLO
  use def_search_method,     only : CANDIDATE_INSIDE
  use def_search_method,     only : CANDIDATE_NEAREST
  use def_search_method,     only : SEARCH_POINTS
  use def_search_method,     only : SEARCH_BOUNDING_BOXES
  use def_search_method,     only : SEARCH_MESH_LEAVES
  use def_search_method,     only : SEARCH_MESH_FILL
  use def_search_method,     only : SEARCH_MESH_ALL
  use def_search_method,     only : big
  use mod_optional_argument, only : optional_argument
  use mod_maths_arrays,      only : maths_maxloc_nonzero
  use mod_maths_geometry,    only : maths_min_max_box_vertices
  use mod_maths_geometry,    only : maths_min_box_vertices
  use mod_maths_sort,        only : maths_heap_sort
  use mod_std

  character(14), parameter :: vacal          = 'def_maths_tree'  
  real(rp),      parameter :: epsil          = epsilon(1.0_rp)
  real(rp),      parameter :: epsil1         = epsilon(1.0_rp)
  real(rp),      parameter :: epsil2         = 2.0_rp*epsilon(1.0_rp)
  real(rp),      parameter :: epsil100       = 100.0_rp*epsilon(1.0_rp)
  integer(ip),   parameter :: SEARCH_OCTREE  = 3
  integer(ip),   parameter :: SEARCH_SKDTREE = 4
  integer(ip),   parameter :: SEARCH_KDTREE  = 5

  type octpoi
     class(octbox), pointer :: p
  end type octpoi

  type octbox
     integer(ip)               :: id          ! My global ID
     integer(ip)               :: idfilled    ! My global ID in filled bins (If 0 does not have any tree nodes, 1 otherwise)
     integer(ip)               :: level       ! Generation
     integer(ip)               :: npoinbox    ! Number of elements (number of entities in the leafs)
     integer(ip)               :: childid     ! Child ID (1->4 or 1->8)
     integer(ip)               :: whoiam      ! Father or have nodes (If 1 does not have children, 0 otherwise)
     integer(ip),  pointer     :: nodes(:)    ! List of entities in box
     real(rp)                  :: minc(3)     ! Min coordinates
     real(rp)                  :: maxc(3)     ! Max coordinates
     type(octbox), pointer     :: parent      ! Pointer to parent
     type(octbox), pointer     :: children(:) ! Pointer to children
     type(octbox), pointer     :: leaf_next   ! Pointer to next leaf
     type(octbox), pointer     :: fill_next   ! Pointer to next filled
     type(octbox), pointer     :: all_next    ! Pointer to next all
  end type octbox
  type(octbox)   :: octbox_init = octbox(&
       0_ip,&                                 ! id         
       0_ip,&                                 ! id filled         
       0_ip,&                                 ! level      
       0_ip,&                                 ! npoinbox   
       0_ip,&                                 ! childid    
       0_ip,&                                 ! whoiam     
       null(),&                               ! nodes(:)   
       (/0.0_rp,0.0_rp,0.0_rp/),&             ! minc(3)    
       (/0.0_rp,0.0_rp,0.0_rp/),&             ! maxc(3)    
       null(),&                               ! parent     
       null(),&                               ! children(:)
       null(),&                               ! leaf_next
       null(),&                               ! fill_next
       null())                                ! all_next

  type, extends(search_method) :: maths_tree
     type(octbox),  pointer  :: tree_root         ! Tree root
     integer(ip)             :: limit             ! Maximum number of elements per bin (input)     
     integer(ip)             :: nleaves           ! Number of leaves
     integer(ip)             :: nfilled           ! Numer of boxes
     integer(ip)             :: nboxes            ! Number of boxes
     integer(ip)             :: divmax            ! Number of divisions
     integer(ip)             :: method            ! method
     logical(lg)             :: adjust_bb         ! Adjusted boxes bounding boxes /useful for nearest point)
     logical(lg)             :: unsorted_division ! Division of the kdtree in a z-y-x axis order instead of searching the longest axis  
     class(octbox), pointer  :: leaf_root         ! Leaf root
     class(octbox), pointer  :: fill_root         ! Fill root
     class(octbox), pointer  :: all_root          ! All root
   contains
     procedure,         pass :: init          ! Initialize all
     procedure,         pass :: deallo        ! Deallocate
     procedure,         pass :: candidate     ! Candidate
     procedure,         pass :: fill          ! Fill
     procedure,         pass :: input         ! Input parameters          
     procedure,         pass :: results       ! Get results on the mesh
     procedure,         pass :: mesh          ! Create a mesh
     procedure,         pass :: graph         ! Create a tree graph
     procedure,         pass :: mesh_level    ! Return level as an element array
     procedure,         pass :: type_name     ! Type name

     procedure,         pass :: divide        ! Divide box and compute children coordinates
     procedure,         pass :: host_bin      ! Pointer to host bin
     procedure,         pass :: near          ! Nearest candidates
     procedure,         pass :: mesh_dim      ! Get a mesh dimensions
     procedure,         pass :: id            ! Get the bin ID of a point
     procedure,         pass :: centroid      ! Centroid coordinates of the bins
  end type maths_tree

  type, extends(maths_tree) :: maths_octree
  end type maths_octree
  
  type, extends(maths_tree) :: maths_skdtree
  end type maths_skdtree
  
  type, extends(maths_tree) :: maths_kdtree
  end type maths_kdtree
    
  type stack
     type(octbox),  pointer  :: p
  end type stack
  
  private

  public :: maths_octree
  public :: maths_skdtree
  public :: maths_kdtree
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Input
  !> @details Set input data
  !> 
  !-----------------------------------------------------------------------

  subroutine input(&
       self,limit,relative_tolerance,absolute_tolerance,fill_method,dim,&
       ADJUST_BOUNDING_BOXES,UNSORTED_DIVISION,NAME,PARAM,DIVMAX)

    class(maths_tree),          intent(inout) :: self
    integer(ip),      optional, intent(in)    :: limit
    real(rp),         optional, intent(in)    :: relative_tolerance
    real(rp),         optional, intent(in)    :: absolute_tolerance
    integer(ip),      optional, intent(in)    :: fill_method
    integer(ip),      optional, intent(in)    :: dim
    logical(lg),      optional, intent(in)    :: ADJUST_BOUNDING_BOXES
    logical(lg),      optional, intent(in)    :: UNSORTED_DIVISION
    character(len=*), optional, intent(in)    :: name
    real(rp),         optional, intent(in)    :: PARAM(:)
    integer(ip),      optional, intent(in)    :: DIVMAX

    call self % input_all(relative_tolerance,absolute_tolerance,fill_method,name)

    if( present(PARAM) ) then
       self % limit = int(param(1),ip)
    end if
    
    if( present(dim)                   ) self % dim               = dim
    if( present(limit)                 ) self % limit             = limit
    if( present(ADJUST_BOUNDING_BOXES) ) self % adjust_bb         = ADJUST_BOUNDING_BOXES
    if( present(UNSORTED_DIVISION)     ) self % unsorted_division = UNSORTED_DIVISION
    if( present(DIVMAX)                ) self % divmax            = DIVMAX

    select type ( self )
    class is ( maths_skdtree )
       self % adjust_bb = .true.
       self % limit     = 1
       self % method    = SEARCH_SKDTREE
       self % name      = 'SKD-TREE'
    class is ( maths_kdtree )
       self % adjust_bb = .false.
       self % method    = SEARCH_KDTREE
       self % name      = 'KD-TREE' 
    class is ( maths_octree )
       self % method    = SEARCH_OCTREE
       self % name      = 'OCT-TREE'
   end select
    
  end subroutine input

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Octree
  !> @details Fill in a octree from a coordinate field
  !> 
  !-----------------------------------------------------------------------

  subroutine init(self)

    class(maths_tree), intent(inout) :: self

    call self % init_all()

    nullify(self % tree_root)
    nullify(self % leaf_root)
    nullify(self % fill_root)
    nullify(self % all_root)
    self % divmax            = 0
    self % nleaves           = 0
    self % nfilled           = 0
    self % nboxes            = 0
    self % limit             = 0
    self % adjust_bb         = .false.
    self % unsorted_division = .false.
    self % method            = SEARCH_OCTREE
    
  end subroutine init

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Deallocate
  !> @details Deallocate octree structure
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo(self,list_entities,MEMORY_COUNTER)

    class(maths_tree),                       intent(inout) :: self
    type(i1p),            optional, pointer, intent(inout) :: list_entities(:)    
    integer(8),           optional,          intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                                             :: memor_loc(2)
    type(octbox), pointer                                  :: current_o
    logical(lg)                                            :: conti

    call self % tim_ini()
    call self % mem_ini(memor_loc,MEMORY_COUNTER)

    self % kfl_filled =  0_ip
    
    if( associated(self % tree_root) ) then

       current_o => self % tree_root
       conti     =  .true.

       do while ( conti )
          !
          ! First go to deepest level in first branch
          !
          do while( current_o % whoiam == 0 )
             current_o => current_o % children(1)
          end do
          !
          ! Deallocate list of elements
          !
          if( current_o % whoiam > 0 ) then
             call memory_deallo(memor_loc,'CURRENT_O % NODES','deallo',current_o % nodes)
          end if

          if( current_o % childid < self % divmax .and. current_o % childid /= 0 ) then
             !
             ! I'm not the last child neither the Padrino
             !
             current_o => current_o % parent % children(current_o % childid+1)

          else if( current_o % childid == self % divmax ) then
             !
             ! I'm the last child
             !
             current_o => current_o % parent
             deallocate(current_o % children)
             current_o % whoiam = -1

          else if( current_o % id == 0 ) then
             !
             ! I'm the Padrino: end of deallocation
             !
             deallocate(current_o)
             nullify(self % tree_root)
             conti = .false.

          end if

       end do
       
    end if

    if( allocated(self % name) ) deallocate(self % name)
    if( present(list_entities) ) call memory_deallo(memor_loc,'LIST_ENTITIES','deallo',list_entities)

    call self % mem_end(memor_loc,MEMORY_COUNTER)
    call self % tim_end(SEARCH_DEALLO)

  end subroutine deallo

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Candidate
  !> @details Candidate in parallel
  !> 
  !-----------------------------------------------------------------------

  subroutine candidate(self,xx,list_entities,METHOD,MASK,MEMORY_COUNTER)

    class(maths_tree),                   intent(inout) :: self
    real(rp),                            intent(in)    :: xx(:,:)           !< List coordinates
    type(i1p),                  pointer, intent(inout) :: list_entities(:)    
    integer(8),       optional,          intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(ip),      optional,          intent(in)    :: METHOD            !< Method for candidates
    logical(lg),      optional, pointer, intent(in)    :: MASK(:)           !< Mask to consider or not points
    integer(ip)                                        :: ii,jj,nn,kk,ll
    integer(8)                                         :: memor_loc(2)
    type(octbox),     pointer                          :: my_bin
    logical(lg)                                        :: mask_loc
    integer(ip)                                        :: my_method
    integer(ip)                                        :: xloc_num
    type(stack),      pointer                          :: lcheck(:)

    call self % tim_ini()
    call self % mem_ini(memor_loc,MEMORY_COUNTER)

    my_method = optional_argument(CANDIDATE_INSIDE,METHOD)
    nn        = int(size(xx,2_ip),ip) 
    mask_loc  = .true.
    nullify(lcheck)
    
    if( .not. associated(list_entities) ) & 
         call memory_alloca(memor_loc,'LIST_ENTITIES',vacal,list_entities,nn)
    
    if( my_method == CANDIDATE_INSIDE ) then
       !
       ! Inside
       !
       do ii = 1,nn
          if( present(MASK) ) then
             if( associated(MASK) ) mask_loc = MASK(ii)
          end if
          if( mask_loc ) then
             
             my_bin => self % host_bin(xx(:,ii))
             
             if( associated(my_bin) ) then
                ll = memory_size(my_bin % nodes)
                self % stats(5) = self % stats(5) + real(ll,rp)                                
                call memory_alloca(memor_loc,'LIST_ENTITIES % L',vacal,list_entities(ii)%l,ll)
                do kk = 1,ll
                   list_entities(ii) % l(kk) = my_bin % nodes(kk)
                end do
             end if
          end if
       end do

    else if( my_method == CANDIDATE_NEAREST ) then
       !
       ! Nearest: eliminate repeated entities using intrinsic ANY()
       !
       allocate(lcheck(self % nfilled))
       do ii = 1,nn
          if( present(MASK) ) then
             if( associated(MASK) ) mask_loc = MASK(ii)   
          end if
          if( mask_loc ) then
             call self % near(xx(:,ii),lcheck,xloc_num,ll)
             if( xloc_num > 0 ) then
                self % stats(5) = self % stats(5) + real(ll,rp)  
                call memory_alloca(memor_loc,'LIST_ENTITIES % L',vacal,list_entities(ii)%l,ll)
                ll = 0
                do jj = 1,xloc_num             
                   do kk = 1,memory_size(lcheck(jj) % p  % nodes)
                      if( .not. any( list_entities(ii) % l == lcheck(jj) % p % nodes(kk)) ) then
                         ll = ll + 1
                         list_entities(ii) % l(ll) = lcheck(jj) % p % nodes(kk)
                      end if
                   end do
                end do
                call memory_resize(memor_loc,'LIST_ENTITIES % L','maths_tree',list_entities(ii) % l,ll)
             end if
          end if
       end do
       if( associated(lcheck) ) deallocate(lcheck)
    end if

    self % stats(5) = self % stats(5) / (real(max(1_ip,self % nelem),rp)*real(nn,rp))

    call self % mem_end(memor_loc,MEMORY_COUNTER)
    call self % tim_end(SEARCH_CANDIDATE)

  end subroutine candidate

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Octree
  !> @details Fill in a octree from a coordinate (COORD) field or
  !>          bounding box field (BOBOX). In this last case, the stopping
  !>          criterion is based on the nu,ber of elements having the
  !>          centroid inside the bin.
  !> 
  !-----------------------------------------------------------------------

  subroutine fill(self,coord,bobox,MEMORY_COUNTER,COORD_MIN,COORD_MAX,PERMU,MASK)

    class(maths_tree),              intent(inout) :: self
    real(rp),    optional, pointer, intent(in)    :: coord(:,:)
    real(rp),    optional, pointer, intent(in)    :: bobox(:,:,:)      !< Bounding boxes
    integer(8),  optional,          intent(inout) :: MEMORY_COUNTER(2)
    real(rp),    optional,          intent(in)    :: COORD_MIN(:)
    real(rp),    optional,          intent(in)    :: COORD_MAX(:)
    integer(ip), optional, pointer, intent(in)    :: PERMU(:)
    logical(lg), optional, pointer, intent(in)    :: MASK(:)           !< Mask to consider or not points
    integer(ip)                                   :: ii,idime,ndime
    integer(ip)                                   :: ipoin,kpoin,npoin
    integer(ip)                                   :: npoin_ini,npoin_end
    integer(ip)                                   :: id,nsize
    integer(ip)                                   :: counter,divmax,imeth
    integer(ip)                                   :: half    
    integer(ip),           pointer                :: sorted_points(:)        
    integer(8)                                    :: memor_loc(2)
    logical(lg)                                   :: conti,not_null_bin
    real(rp)                                      :: minsize(3)
    type(octbox),          pointer                :: tree_root
    type(octbox),          pointer                :: current_o
    type(octbox),          pointer                :: current_tmp
    type(octbox),          pointer                :: leaf_o
    type(octbox),          pointer                :: fill_o
    type(octbox),          pointer                :: all_o
    type(octbox),          pointer                :: old_pointer
    type(octbox),          pointer                :: tm1_pointer
    type(octbox),          pointer                :: tm2_pointer    
    real(rp),              pointer                :: centr(:,:)
    real(rp),              pointer                :: min_loc(:,:)
    real(rp),              pointer                :: max_loc(:,:)
    logical(lg)                                   :: if_coord
    logical(lg)                                   :: if_bobox

    ! FOR KDTREE: Point ids that will be sorted after: sorted_points
    nullify(centr,min_loc,max_loc,sorted_points)

    call self % tim_ini()
    call self % mem_ini(memor_loc,MEMORY_COUNTER)
    imeth     = 0
    self % kfl_filled = 1_ip

    !--------------------------------------------------------------------
    !
    ! Check errors
    !
    !--------------------------------------------------------------------

    if( self % limit < 1 ) call runend('DEF_MATHS_TREE: MAXIMUM NUMBER OF ENTITIES IN BOXES SHOULD BE >= 1 ')

    !--------------------------------------------------------------------
    !
    ! Coordinates or bounding boxes
    !
    !--------------------------------------------------------------------

    if_coord  = .false.
    if_bobox  = .false.
    npoin_ini =  0
    npoin_end = -1
    if( present(coord) ) then
       self % fill_method = SEARCH_POINTS
       if( associated(coord) ) if_coord = .true.
    end if
    if( present(bobox) ) then
       self % fill_method = SEARCH_BOUNDING_BOXES
       if( associated(bobox) ) if_bobox = .true.
    end if
    if( if_coord ) then
       imeth = 1
    else if( if_bobox ) then
       imeth = 2
    else
       return
    end if

    if( present(coord)      .and. if_coord ) then
       imeth      = 1
       ndime      = memory_size(coord,1_ip)
       npoin_ini  = lbound(coord,DIM=2_ip,KIND=ip)
       npoin_end  = ubound(coord,DIM=2_ip,KIND=ip)
    else if( present(bobox) .and. if_bobox ) then
       imeth      = 2
       ndime      = memory_size(bobox,2_ip)
       npoin_ini  = lbound(bobox,DIM=3_ip,KIND=ip)
       npoin_end  = ubound(bobox,DIM=3_ip,KIND=ip)
    end if

    !--------------------------------------------------------------------
    !
    ! How to divide
    !
    !--------------------------------------------------------------------

    select case ( self % method )
    case ( SEARCH_OCTREE  ) ; divmax = 2**ndime 
    case ( SEARCH_SKDTREE ) ; divmax = 2        
    case ( SEARCH_KDTREE  ) ; divmax = 2        
    end select

    !--------------------------------------------------------------------
    !
    ! Dimensions
    !
    !--------------------------------------------------------------------

    npoin         = npoin_end-npoin_ini+1
    self % dim    = ndime
    self % divmax = divmax
    self % nelem  = npoin
    if( npoin <= 0 ) return

    !--------------------------------------------------------------------
    !
    ! Allocate tree root
    !
    !--------------------------------------------------------------------

    allocate(self % tree_root) 
    self % tree_root = octbox_init

    !--------------------------------------------------------------------
    !
    ! In case we adjust boundaries
    !
    !--------------------------------------------------------------------

    call memory_alloca(memor_loc,'MIN_LOC','fill',min_loc,ndime,divmax) 
    call memory_alloca(memor_loc,'MAX_LOC','fill',max_loc,ndime,divmax) 

    !--------------------------------------------------------------------
    !
    ! Padrino's (root) bounding box
    !
    !--------------------------------------------------------------------

    if(      imeth == 1 .and. present(coord) ) then

       call self % bounding_box(coord,OFFSET=self % toler_rel,COORD_MIN=COORD_MIN,COORD_MAX=COORD_MAX)
       minsize = 0.0_rp

    else if( imeth == 2 .and. present(bobox) ) then
       !
       ! Use centroids only according stopping citerion
       !
       call self % bounding_box(bobox,OFFSET=self % toler_rel,COORD_MIN=COORD_MIN,COORD_MAX=COORD_MAX)
       minsize = big
       do ipoin = npoin_ini,npoin_end
          do idime = 1,ndime
             !
             ! Check if a boundary box dimension is equal to 0
             ! It can happens for a boundary box of a surface that is parallel to an axis
             !
             if( bobox(2,idime,ipoin)-bobox(1,idime,ipoin) <= epsil100 ) then
                bobox(2,idime,ipoin) = bobox(2,idime,ipoin) + epsil100 
                bobox(1,idime,ipoin) = bobox(1,idime,ipoin) - epsil100 
             end if
             minsize(idime) = min(minsize(idime),abs(bobox(2,idime,ipoin)-bobox(1,idime,ipoin)))                
          end do
       end do
    end if

    !--------------------------------------------------------------------
    !
    ! Fill in Padrino
    !
    !--------------------------------------------------------------------

    tree_root                 => self % tree_root
    current_o                 => tree_root
    current_o % id            =  0
    current_o % level         =  0
    current_o % childid       =  0
    current_o % whoiam        =  0
    current_o % minc(1:ndime) =  self % comin(1:ndime)
    current_o % maxc(1:ndime) =  self % comax(1:ndime)

    !--------------------------------------------------------------------
    !
    ! Initial nodes, and remove wrong nodes with bad boundary boxes
    ! (min > max)
    !
    !--------------------------------------------------------------------

    if( present(bobox) ) then
       current_o % npoinbox = 0
       do ipoin = npoin_ini,npoin_end
          if ( (bobox(1,1,ipoin) < sign(10.0_rp,bobox(2,1,ipoin))*bobox(2,1,ipoin)) .and. &
               (bobox(1,2,ipoin) < sign(10.0_rp,bobox(2,2,ipoin))*bobox(2,2,ipoin)) .and. &
               (bobox(1,ndime,ipoin) < sign(10.0_rp,bobox(2,ndime,ipoin))*bobox(2,ndime,ipoin)) ) then          
             current_o % npoinbox = current_o % npoinbox + 1
          end if
       end do
       call memory_alloca(memor_loc,'CURRENT_O % NODES','fill',self % tree_root % nodes,current_o % npoinbox)
       kpoin = 0
       do ipoin = npoin_ini,npoin_end
          if ( (bobox(1,1,ipoin) < sign(10.0_rp,bobox(2,1,ipoin))*bobox(2,1,ipoin)) .and. &
               (bobox(1,2,ipoin) < sign(10.0_rp,bobox(2,2,ipoin))*bobox(2,2,ipoin)) .and. &
               (bobox(1,ndime,ipoin) < sign(10.0_rp,bobox(2,ndime,ipoin))*bobox(2,ndime,ipoin)) ) then          
             kpoin = kpoin + 1
             current_o % nodes(kpoin) = ipoin
          end if
       end do
    else
       call memory_alloca(memor_loc,'CURRENT_O % NODES','fill',self % tree_root % nodes,npoin)!,INIT_VALUE=-1_ip)
       current_o % npoinbox = npoin
       kpoin = npoin_ini-1
       do ipoin = 1,npoin
          kpoin = kpoin + 1
          current_o % nodes(ipoin) = kpoin
       end do
    end if
    nullify( current_o % parent )

    !--------------------------------------------------------------------
    !
    ! Additional arrays required by SKD-TREE
    !
    !--------------------------------------------------------------------
    
    if( imeth == 2 .and. present(bobox) .and. self % method == SEARCH_SKDTREE ) then
       call memory_alloca(memor_loc,'CENTR'        ,'fill',centr,ndime,npoin,LBOUN2=npoin_ini)
       call memory_alloca(memor_loc,'SORTED_POINTS','fill',sorted_points,current_o % npoinbox)                 
       do ipoin = npoin_ini,npoin_end
          centr(1:ndime,ipoin) = (bobox(1,1:ndime,ipoin)+bobox(2,1:ndime,ipoin))*0.5_rp
       end do
    end if
    
    !--------------------------------------------------------------------
    !
    ! Traverse the oct-tree
    !
    !--------------------------------------------------------------------

    counter = 0
    conti   = .true.
    all_o           => self % tree_root
    self % all_root => self % tree_root
    nullify(old_pointer)
    
    do while( conti )
       !
       ! Check null bins or bins that are too small
       !
       not_null_bin = .true.

       if( present(bobox) .and. current_o % id /= 0 ) then
          if( self % method == SEARCH_OCTREE .or. self % method == SEARCH_KDTREE ) then
             minsize      = big
             do ipoin = 1,current_o % npoinbox
                kpoin = current_o % nodes(ipoin)
                do idime = 1,ndime
                   minsize(idime) = min(minsize(idime),bobox(2,idime,kpoin)-bobox(1,idime,kpoin))
                end do
             end do
             do idime = 1,ndime
                if( abs(current_o % maxc(idime)-current_o % minc(idime)) < minsize(idime)+epsil ) not_null_bin = .false.
             end do
             if( .not. not_null_bin .and. current_o % id == 0 ) not_null_bin = .true.
          end if
       end if
       
       if( current_o % npoinbox > self % limit .and. not_null_bin ) then
          !
          ! If maximum number of points inside current box is exceeded, subdivide
          !
          allocate( current_o % children(divmax) )
          current_o % children(1:divmax) = octbox_init
          !
          ! Give birth to my DIVMAX children
          !
          all_o % all_next => current_o % children(1)
          all_o            => all_o % all_next
          do ii = 1,divmax
             counter                             =  counter+1
             current_o % children(ii) % id       =  counter
             current_o % children(ii) % childid  =  ii
             current_o % children(ii) % level    =  current_o % level + 1 
             current_o % children(ii) % whoiam   =  0  
             current_o % children(ii) % npoinbox =  0
             current_o % children(ii) % idfilled =  0
             current_o % children(ii) % parent   => current_o

             if( ii < divmax ) then
                all_o % all_next => current_o % children(ii+1)
                all_o            => all_o % all_next
             end if
             call memory_alloca(memor_loc,'CURRENT_O % NODES','fill',current_o % children(ii) % nodes,current_o % npoinbox)
          end do
          !
          ! Divide and compute children coordinates
          !
          call self % divide(current_o,sorted_points,coord,centr)        
          !
          ! Min sizes
          !
          min_loc =  big
          max_loc = -big
          !
          ! Offer my nodes to my children
          !
          if( imeth == 1 .and. present(coord) ) then
             !
             ! Coordinates
             !
             if( ndime == 2 ) then
                !
                ! COORD 2D
                !
                do ipoin = 1,current_o % npoinbox
                   kpoin = current_o % nodes(ipoin)
                   if( optional_argument(.true.,MASK,kpoin) ) then
                      loop_coord_2d: do ii = 1,divmax
                         if(  &
                              & coord(1,kpoin) >= current_o % children(ii) % minc(1) - epsil .and. &
                              & coord(2,kpoin) >= current_o % children(ii) % minc(2) - epsil .and. &
                              & coord(1,kpoin) <= current_o % children(ii) % maxc(1) + epsil .and. &
                              & coord(2,kpoin) <= current_o % children(ii) % maxc(2) + epsil ) then
                            current_o % children(ii) % npoinbox     = current_o % children(ii) % npoinbox + 1
                            nsize                                   = current_o % children(ii) % npoinbox
                            current_o % children(ii) % nodes(nsize) = kpoin
                            min_loc(1:2,ii)                         = min(min_loc(1:2,ii),coord(1:2,kpoin))
                            max_loc(1:2,ii)                         = max(max_loc(1:2,ii),coord(1:2,kpoin))
                            exit loop_coord_2d
                         end if
                      end do loop_coord_2d
                   end if
                end do
                if( self % adjust_bb ) then
                   do ii = 1,divmax
                      current_o % children(ii) % minc(1:2) = min_loc(1:2,ii) - epsil1  
                      current_o % children(ii) % maxc(1:2) = max_loc(1:2,ii) + epsil1
                   end do
                end if

             else if( ndime == 3 ) then
                !
                ! COORD 3D
                !
                do ipoin = 1,current_o % npoinbox
                   kpoin = current_o % nodes(ipoin)
                   if( optional_argument(.true.,MASK,kpoin) ) then
                      loop_coord_3d: do ii = 1,divmax
                         if(                coord(1,kpoin) >= current_o % children(ii) % minc(1) - epsil ) then
                            if(             coord(2,kpoin) >= current_o % children(ii) % minc(2) - epsil ) then
                               if(          coord(3,kpoin) >= current_o % children(ii) % minc(3) - epsil ) then 
                                  if(       coord(1,kpoin) <= current_o % children(ii) % maxc(1) + epsil ) then 
                                     if(    coord(2,kpoin) <= current_o % children(ii) % maxc(2) + epsil ) then 
                                        if( coord(3,kpoin) <= current_o % children(ii) % maxc(3) + epsil ) then
                                           current_o % children(ii) % npoinbox = current_o % children(ii) % npoinbox + 1
                                           current_o % children(ii) % nodes(current_o % children(ii) % npoinbox) = kpoin
                                           min_loc(1:3,ii) = min(min_loc(1:3,ii),coord(1:3,kpoin))
                                           max_loc(1:3,ii) = max(max_loc(1:3,ii),coord(1:3,kpoin))
                                           exit loop_coord_3d
                                        end if
                                     end if
                                  end if
                               end if
                            end if
                         end if
                      end do loop_coord_3d
                   end if
                end do
                if( self % adjust_bb ) then
                   do ii = 1,divmax
                      current_o % children(ii) % minc(1:3) = min_loc(1:3,ii) - epsil1
                      current_o % children(ii) % maxc(1:3) = max_loc(1:3,ii) + epsil1
                   end do
                end if

             end if

          else if( imeth == 2 .and. present(bobox) ) then
             !
             ! Bounding boxes
             !
             if( self % method == SEARCH_OCTREE .or. self % method == SEARCH_KDTREE ) then

                if( ndime == 2 ) then
                   !
                   ! Count an fill
                   !
                   do ipoin = 1,current_o % npoinbox
                      kpoin = current_o % nodes(ipoin)
                      if( optional_argument(.true.,MASK,kpoin) ) then
                         loop_bobox_2d_1: do ii = 1,divmax
                            if(          bobox(1,1,kpoin) <= current_o % children(ii) % maxc(1) + epsil ) then
                               if(       bobox(2,1,kpoin) >= current_o % children(ii) % minc(1) - epsil ) then
                                  if(    bobox(1,2,kpoin) <= current_o % children(ii) % maxc(2) + epsil ) then
                                     if( bobox(2,2,kpoin) >= current_o % children(ii) % minc(2) - epsil ) then
                                        current_o % children(ii) % npoinbox     = current_o % children(ii) % npoinbox + 1
                                        nsize                                   = current_o % children(ii) % npoinbox
                                        current_o % children(ii) % nodes(nsize) = kpoin
                                        if( self % adjust_bb ) then
                                           min_loc(1:2,ii)                      = min(min_loc(1:2,ii),bobox(1,1:2,kpoin))
                                           max_loc(1:2,ii)                      = max(max_loc(1:2,ii),bobox(2,1:2,kpoin))
                                        end if
                                     end if
                                  end if
                               end if
                            end if
                         end do loop_bobox_2d_1
                      end if
                   end do
                   if( self % adjust_bb ) then
                      do ii = 1,divmax
                         current_o % children(ii) % minc(1:2) = min_loc(1:2,ii) - epsil1
                         current_o % children(ii) % maxc(1:2) = max_loc(1:2,ii) + epsil1
                      end do
                   end if

                else if( ndime == 3 ) then                
                   !
                   ! Count and fill
                   !
                   do ipoin = 1,current_o % npoinbox
                      kpoin = current_o % nodes(ipoin)
                      if( optional_argument(.true.,MASK,kpoin) ) then
                         loop_bobox_3d_1: do ii = 1,divmax
                            if(                   bobox(2,1,kpoin) >= current_o % children(ii) % minc(1) - epsil ) then
                               if(                bobox(2,2,kpoin) >= current_o % children(ii) % minc(2) - epsil ) then
                                  if(             bobox(2,3,kpoin) >= current_o % children(ii) % minc(3) - epsil ) then
                                     if(          bobox(1,1,kpoin) <= current_o % children(ii) % maxc(1) + epsil ) then
                                        if(       bobox(1,2,kpoin) <= current_o % children(ii) % maxc(2) + epsil ) then
                                           if(    bobox(1,3,kpoin) <= current_o % children(ii) % maxc(3) + epsil ) then
                                              current_o % children(ii) % npoinbox     = current_o % children(ii) % npoinbox + 1
                                              nsize                                   = current_o % children(ii) % npoinbox
                                              current_o % children(ii) % nodes(nsize) = kpoin
                                              if( self % adjust_bb ) then
                                                 min_loc(1:3,ii)                      = min(min_loc(1:3,ii),bobox(1,1:3,kpoin))
                                                 max_loc(1:3,ii)                      = max(max_loc(1:3,ii),bobox(2,1:3,kpoin))
                                              end if
                                           end if
                                        end if
                                     end if
                                  end if
                               end if
                            end if
                         end do loop_bobox_3d_1

                      end if
                   end do

                   if( self % adjust_bb ) then
                      do ii = 1,divmax
                         current_o % children(ii) % minc(1:3) = min_loc(1:3,ii) - epsil1
                         current_o % children(ii) % maxc(1:3) = max_loc(1:3,ii) + epsil1
                      end do
                   end if

                end if

             else if( self % method == SEARCH_SKDTREE ) then

                half = current_o % npoinbox/2_ip
                !
                ! Fill first child
                !
                do ipoin = 1,half
                   kpoin = sorted_points(ipoin)
                   if( optional_argument(.true.,MASK,kpoin) ) then                   
                      current_o % children(1) % npoinbox     = current_o % children(1) % npoinbox + 1
                      nsize                                  = current_o % children(1) % npoinbox
                      current_o % children(1) % nodes(nsize) = kpoin
                      min_loc(1:ndime,1)                     = min(min_loc(1:ndime,1),bobox(1,1:ndime,kpoin))
                      max_loc(1:ndime,1)                     = max(max_loc(1:ndime,1),bobox(2,1:ndime,kpoin))                      
                   end if                   
                end do
                ! Determine new BB
                current_o % children(1) % minc(1:ndime) = min_loc(1:ndime,1) - epsil1
                current_o % children(1) % maxc(1:ndime) = max_loc(1:ndime,1) + epsil1
                !
                ! Fill second child
                !
                do ipoin = half+1,current_o % npoinbox
                   kpoin = sorted_points(ipoin)
                   if( optional_argument(.true.,MASK,kpoin) ) then                   
                      current_o % children(2) % npoinbox     = current_o % children(2) % npoinbox + 1
                      nsize                                  = current_o % children(2) % npoinbox
                      current_o % children(2) % nodes(nsize) = kpoin
                      min_loc(1:ndime,2)                     = min(min_loc(1:ndime,2),bobox(1,1:ndime,kpoin))
                      max_loc(1:ndime,2)                     = max(max_loc(1:ndime,2),bobox(2,1:ndime,kpoin))                      
                   end if                   
                end do                
                ! Determine new BB
                current_o % children(2) % minc(1:ndime) = min_loc(1:ndime,2) - epsil1
                current_o % children(2) % maxc(1:ndime) = max_loc(1:ndime,2) + epsil1
                
             end if
          end if
          !
          ! Deallocate parent list of nodes
          !
          call memory_deallo(memor_loc,'CURRENT_O % NODES','fill',current_o % nodes)

          current_o % whoiam   =  0
          current_o % npoinbox =  0

          !if( current_o % children(ii) % id /= 0 ) all_o % all_next => current_o % children(1)

          current_o            => current_o % children(1)

       else if( current_o % id == 0 .and. current_o % npoinbox <= self % limit ) then
          !
          ! If the Padrino has too few elements
          !          
          conti              =  .false.
          self % nleaves     =  self % nleaves + 1
          current_o % whoiam =  self % nleaves
          if( current_o % npoinbox > 0 ) then
             self % nfilled       =  self % nfilled + 1
             current_o % idfilled =  self % nfilled
          end if
          current_o          => old_pointer

       else 
          !
          ! if limit of points inside box is not exceeded, assign nodes
          !
          call memory_resize(memor_loc,'CURRENT_O % NODES','fill',current_o % nodes,current_o % npoinbox)
          self      % nleaves = self % nleaves + 1
          current_o % whoiam  = self % nleaves
          if( current_o % npoinbox > 0 ) then
             self      % nfilled  = self % nfilled + 1
             current_o % idfilled = self % nfilled
          end if

          if( current_o % childid < divmax .and. current_o % id /= 0 ) then
             !
             ! Go to next children
             !
             tm1_pointer => current_o
             tm2_pointer => tm1_pointer % parent % children(tm1_pointer%childid+1)
             current_o   => tm2_pointer
             goto 10

          else if(current_o % childid == 0 ) then  
             !
             ! Padrino
             !
             goto 10

          else if(current_o % childid == divmax ) then
             !
             ! Last children
             !
             noparent: do while( current_o % id > 0 )
                if(current_o % parent % id == 0) then
                   conti = .false.
                   exit noparent
                else
                   if( current_o % parent % childid /= divmax ) then 
                      tm1_pointer => current_o
                      tm2_pointer => tm1_pointer % parent % parent % children(tm1_pointer%parent%childid+1)
                      current_o   => tm2_pointer
                      exit
                   else 
                      current_o   => current_o % parent
                   end if
                end if
             end do noparent

          else 
             !
             ! Wrong child ID
             !
             stop

          end if

       end if

10     continue

       old_pointer => current_o

    end do

    !--------------------------------------------------------------------
    !
    ! Compute leaf and fill pointers, and permute if required
    !
    !--------------------------------------------------------------------

    self % nboxes =  counter
    current_o     => self % tree_root
    leaf_o        => self % tree_root
    fill_o        => self % tree_root
    conti         = .true.

    loop_conti: do while( conti )
       !
       ! First go to deepest level in first branch
       !
       do while( current_o % whoiam == 0 )
          current_tmp => current_o % children(1)
          current_o   => current_tmp
       end do
       !
       ! Statistics
       !
       self % stats(1) = self % stats(1) + 1.0_rp
       self % stats(2) = self % stats(2) + real(current_o % npoinbox,rp)
       self % stats(3) = max(self % stats(3),real(current_o % level,rp))
       self % stats(4) = max(self % stats(4),real(current_o % npoinbox,rp))
       !
       ! Leaf trasverse
       !
       if( .not. associated(self % leaf_root) ) then
          self   % leaf_root => current_o
       else 
          leaf_o % leaf_next => current_o          
       end if
       leaf_o => current_o
       !
       ! Fill trasverse
       !
       if( current_o % idfilled /= 0 ) then
          if( .not. associated(self % fill_root) ) then
             self   % fill_root => current_o
          else 
             fill_o % fill_next => current_o      
          end if
          fill_o => current_o          
       end if
       !
       ! Permute list of nodes if ncessary
       !
       if( present(permu) ) then
          do ii = 1,memory_size(current_o % nodes)
             ipoin = current_o % nodes(ii)
             current_o % nodes(ii) = permu(ipoin)
          end do
       end if

       if(current_o % childid < divmax .and. current_o % childid /=0 ) then
          !
          ! I'm not the last child neither the Padrino
          !
          id          =  current_o % childid
          current_tmp => current_o % parent % children(id+1)
          current_o   => current_tmp

       else if( current_o % childid == divmax ) then
          !
          ! I'm the last child of this generation
          !
          do while( current_o % id > 0 )
             if(current_o % parent % id == 0) then
                conti = .false. 
                exit loop_conti
             else
                if( current_o % parent % childid /= divmax ) then
                   id          =  current_o % parent % childid
                   current_tmp => current_o % parent % parent % children(id+1)
                   current_o   => current_tmp
                   exit
                else 
                   current_tmp => current_o % parent
                   current_o   => current_tmp
                end if
             end if
          end do

       else if( current_o % id == 0 ) then
          !
          ! I'm the Padrino
          !
          conti = .false.

       end if

    end do loop_conti

    self % stats(2) = self % stats(2) / max(real(self % nelem,rp),1.0_rp)
    self % kfl_filled = 2

    !--------------------------------------------------------------------
    !
    ! Deallocate
    !
    !--------------------------------------------------------------------

    call memory_deallo(memor_loc,'CENTR'        ,'fill',centr)    
    call memory_deallo(memor_loc,'MIN_LOC'      ,'fill',min_loc) 
    call memory_deallo(memor_loc,'MAX_LOC'      ,'fill',max_loc)
    call memory_deallo(memor_loc,'SORTED_POINTS','fill',sorted_points)

    call self % mem_end(memor_loc,MEMORY_COUNTER)
    call self % tim_end(SEARCH_FILL)

  end subroutine fill

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   Get mesh dimensions
  !> @details Get the mesh dimensions
  !> 
  !-----------------------------------------------------------------------

  subroutine mesh_dim(octree,ndime,mnode,nelem,npoin,CRITERION)

    class(maths_tree),                      intent(in)    :: octree       !< Octree
    integer(ip),                            intent(out)   :: ndime        !< Space dimension
    integer(ip),                            intent(out)   :: mnode        !< Max number of nodes per element
    integer(ip),                            intent(out)   :: nelem        !< Number of elements
    integer(ip),                            intent(out)   :: npoin        !< Number of nodes
    integer(ip),         optional,          intent(in)    :: CRITERION    !< If empty bins should be considered
    integer(ip)                                           :: my_criterion

    my_criterion = optional_argument(SEARCH_MESH_LEAVES , CRITERION)

    select case ( my_criterion )
    case ( SEARCH_MESH_ALL    ) ; nelem = octree % nboxes + 1
    case ( SEARCH_MESH_LEAVES ) ; nelem = octree % nleaves
    case ( SEARCH_MESH_FILL   ) ; nelem = octree % nfilled
    end select

    ndime =  octree % dim
    mnode =  (ndime-1)*4
    npoin =  nelem * mnode

  end subroutine mesh_dim

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   Get level
  !> @details Get level
  !> 
  !-----------------------------------------------------------------------
  
  subroutine mesh_level(self,level,CRITERION,MEMORY_COUNTER,ONLY_DEALLOCATE)

    class(maths_tree),                      intent(inout) :: self
    integer(ip),                   pointer, intent(inout) :: level(:)
    integer(ip),         optional,          intent(in)    :: CRITERION     !< If empty bins should be considered
    integer(8),          optional,          intent(inout) :: MEMORY_COUNTER(2)
    logical(lg),         optional,          intent(in)    :: ONLY_DEALLOCATE
    integer(ip)                                           :: ndime
    integer(ip)                                           :: mnode
    integer(ip)                                           :: nelem
    integer(ip)                                           :: npoin,ielem
    type(octbox),                  pointer                :: current_o
    integer(8)                                            :: memor_loc(2)
    integer(ip)                                           :: my_criterion

    call self % mem_ini(memor_loc,MEMORY_COUNTER)

    if( optional_argument(.false.,ONLY_DEALLOCATE) ) then
       
       call memory_deallo(memor_loc,'LEVEL','mesh',level)

    else

       call self % mesh_dim(ndime,mnode,nelem,npoin,CRITERION)
       if( .not. associated(level) ) call memory_alloca(memor_loc,'LEVEL','mesh',level,nelem)

       my_criterion = optional_argument(SEARCH_MESH_LEAVES , CRITERION)

       select case ( my_criterion )
       case ( SEARCH_MESH_ALL    ) ; current_o => self % all_root
       case ( SEARCH_MESH_LEAVES ) ; current_o => self % leaf_root
       case ( SEARCH_MESH_FILL   ) ; current_o => self % fill_root
       end select

       do while( associated(current_o) )

          select case ( my_criterion )  
          case ( SEARCH_MESH_ALL    ) ; ielem = current_o % id + 1      
          case ( SEARCH_MESH_LEAVES ) ; ielem = current_o % whoiam  
          case ( SEARCH_MESH_FILL   ) ; ielem = current_o % idfilled 
          end select

          level(ielem) = current_o % level

          select case ( my_criterion )
          case ( SEARCH_MESH_ALL    ) ; current_o => current_o % all_next
          case ( SEARCH_MESH_LEAVES ) ; current_o => current_o % leaf_next
          case ( SEARCH_MESH_FILL   ) ; current_o => current_o % fill_next  
          end select

       end do

    end if

    call self % mem_end(memor_loc,MEMORY_COUNTER)

  end subroutine mesh_level
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   Create a mesh
  !> @details Create a mesh from an octree
  !> 
  !-----------------------------------------------------------------------

  subroutine mesh(self,ndime,mnode,nelem,npoin,lnods,ltype,coord,MEMORY_COUNTER,&
       OFFSET_IELEM,OFFSET_IPOIN,CRITERION,CENTROID,ONLY_DEALLOCATE)

    class(maths_tree),                      intent(inout) :: self
    integer(ip),                            intent(out)   :: ndime
    integer(ip),                            intent(out)   :: mnode
    integer(ip),                            intent(out)   :: nelem
    integer(ip),                            intent(out)   :: npoin
    integer(ip),                   pointer, intent(inout) :: lnods(:,:)
    integer(ip),                   pointer, intent(inout) :: ltype(:)
    real(rp),                      pointer, intent(inout) :: coord(:,:)
    integer(8),          optional,          intent(inout) :: MEMORY_COUNTER(2)
    integer(ip),         optional,          intent(in)    :: OFFSET_IELEM
    integer(ip),         optional,          intent(in)    :: OFFSET_IPOIN
    integer(ip),         optional,          intent(in)    :: CRITERION     !< If empty bins should be considered
    real(rp),            optional, pointer, intent(inout) :: CENTROID(:,:)
    logical(lg),         optional,          intent(in)    :: ONLY_DEALLOCATE
    type(octbox),                  pointer                :: current_o
    integer(ip)                                           :: ipoin,ielem,idime
    logical(lg)                                           :: conti
    integer(8)                                            :: memor_loc(2)
    integer(ip)                                           :: offset_ielem_loc
    integer(ip)                                           :: offset_ipoin_loc
    logical(lg)                                           :: if_deallocate
    integer(ip)                                           :: my_criterion
    
    call self % mem_ini(memor_loc,MEMORY_COUNTER)

    offset_ielem_loc = optional_argument(0_ip               , OFFSET_IELEM)
    offset_ipoin_loc = optional_argument(0_ip               , OFFSET_IPOIN)
    my_criterion     = optional_argument(SEARCH_MESH_LEAVES , CRITERION)
    if_deallocate    = optional_argument(.false.            , ONLY_DEALLOCATE)

    if( if_deallocate ) then
       !
       ! Deallocate
       !
       call memory_deallo(memor_loc,'LNODS','mesh',lnods)
       call memory_deallo(memor_loc,'LTYPE','mesh',ltype)
       call memory_deallo(memor_loc,'COORD','mesh',coord)
       if( present(CENTROID) ) then
          if( .not. associated(CENTROID) ) call memory_deallo(memor_loc,'CENTROID','mesh',CENTROID)
       end if

    else
       !
       ! Allocate and compute mesh
       !
       call self % mesh_dim(ndime,mnode,nelem,npoin,CRITERION)
       if( .not. associated(lnods) ) call memory_alloca(memor_loc,'LNODS','mesh',lnods,mnode,nelem)
       if( .not. associated(ltype) ) call memory_alloca(memor_loc,'LTYPE','mesh',ltype,nelem)
       if( .not. associated(coord) ) call memory_alloca(memor_loc,'COORD','mesh',coord,ndime,npoin)
       if( present(CENTROID) ) then
          if( .not. associated(CENTROID) ) call memory_alloca(memor_loc,'CENTROID','mesh',CENTROID,ndime,nelem)
       end if
       !
       ! Mesh arrays
       !
       ipoin =  1 + offset_ipoin_loc
       conti =  .true.

       select case ( my_criterion )
       case ( SEARCH_MESH_ALL    ) ; current_o => self % all_root
       case ( SEARCH_MESH_LEAVES ) ; current_o => self % leaf_root
       case ( SEARCH_MESH_FILL   ) ; current_o => self % fill_root
       end select

       do while( associated(current_o) )

          select case ( my_criterion )
          case ( SEARCH_MESH_ALL    ) ; ielem = current_o % id       + offset_ielem_loc + 1
          case ( SEARCH_MESH_LEAVES ) ; ielem = current_o % whoiam   + offset_ielem_loc
          case ( SEARCH_MESH_FILL   ) ; ielem = current_o % idfilled + offset_ielem_loc      
          end select
          
          if( ndime == 2 ) then
             coord(1:2,ipoin)   = (/ current_o % minc(1),current_o % minc(2) /)
             coord(1:2,ipoin+1) = (/ current_o % maxc(1),current_o % minc(2) /)
             coord(1:2,ipoin+2) = (/ current_o % maxc(1),current_o % maxc(2) /)
             coord(1:2,ipoin+3) = (/ current_o % minc(1),current_o % maxc(2) /)
             lnods(1,ielem)     = ipoin 
             lnods(2,ielem)     = ipoin + 1
             lnods(3,ielem)     = ipoin + 2
             lnods(4,ielem)     = ipoin + 3
             ltype(ielem)       = QUA04
             if( present(CENTROID) ) then
                do idime = 1,2
                   CENTROID(idime,ielem) = sum(coord(idime,ipoin:ipoin+3))/4.0_rp
                end do
             end if
             ipoin              = ipoin + 4
          else
             coord(1:3,ipoin)   = (/ current_o % minc(1),current_o % minc(2),current_o % minc(3) /)
             coord(1:3,ipoin+1) = (/ current_o % maxc(1),current_o % minc(2),current_o % minc(3) /)
             coord(1:3,ipoin+2) = (/ current_o % maxc(1),current_o % maxc(2),current_o % minc(3) /)
             coord(1:3,ipoin+3) = (/ current_o % minc(1),current_o % maxc(2),current_o % minc(3) /)
             coord(1:3,ipoin+4) = (/ current_o % minc(1),current_o % minc(2),current_o % maxc(3) /)
             coord(1:3,ipoin+5) = (/ current_o % maxc(1),current_o % minc(2),current_o % maxc(3) /)
             coord(1:3,ipoin+6) = (/ current_o % maxc(1),current_o % maxc(2),current_o % maxc(3) /)
             coord(1:3,ipoin+7) = (/ current_o % minc(1),current_o % maxc(2),current_o % maxc(3) /)
             lnods(1,ielem)     = ipoin 
             lnods(2,ielem)     = ipoin + 1
             lnods(3,ielem)     = ipoin + 2
             lnods(4,ielem)     = ipoin + 3
             lnods(5,ielem)     = ipoin + 4
             lnods(6,ielem)     = ipoin + 5
             lnods(7,ielem)     = ipoin + 6
             lnods(8,ielem)     = ipoin + 7
             ltype(ielem)       = HEX08
             if( present(CENTROID) ) then
                do idime = 1,3
                   CENTROID(idime,ielem) = sum(coord(idime,ipoin:ipoin+7))/8.0_rp
                end do
             end if
             ipoin              = ipoin + 8             
          end if

          select case ( my_criterion )
          case ( SEARCH_MESH_ALL    ) ; current_o => current_o % all_next
          case ( SEARCH_MESH_LEAVES ) ; current_o => current_o % leaf_next
          case ( SEARCH_MESH_FILL   ) ; current_o => current_o % fill_next  
          end select

       end do

    end if

    call self % mem_end(memor_loc,MEMORY_COUNTER)

  end subroutine mesh

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   Create a mesh
  !> @details Create a mesh from an octree
  !> 
  !-----------------------------------------------------------------------

  subroutine graph(self,CRITERION)

    class(maths_tree),                      intent(inout) :: self    
    integer(ip),         optional,          intent(in)    :: CRITERION     !< If empty bins should be considered
    type(octbox),                  pointer                :: current_o
    integer(ip)                                           :: ielem,ifath
    integer(ip)                                           :: ii,iw
    integer(ip)                                           :: my_criterion

    my_criterion     = optional_argument(SEARCH_MESH_LEAVES , CRITERION)
    !
    ! Mesh arrays
    !
    select case ( my_criterion )
    case ( SEARCH_MESH_ALL    ) ; current_o => self % all_root
    case ( SEARCH_MESH_LEAVES ) ; current_o => self % leaf_root
    case ( SEARCH_MESH_FILL   ) ; current_o => self % fill_root
    end select

    do while( associated(current_o) )

       select case ( my_criterion )
       case ( SEARCH_MESH_ALL    ) ; ifath = current_o % id       + 1
       case ( SEARCH_MESH_LEAVES ) ; ifath = current_o % whoiam    
       case ( SEARCH_MESH_FILL   ) ; ifath = current_o % idfilled        
       end select
       
    !print*,'aaa=',associated(current_o % children),self % divmax
    
       if( associated(current_o % children) ) then
          do ii = 1,self % divmax
             select case ( my_criterion )
             case ( SEARCH_MESH_ALL    ) ; ielem = current_o % children(ii) % id + 1   ; iw = current_o % children(ii) % npoinbox
             case ( SEARCH_MESH_LEAVES ) ; ielem = current_o % children(ii) % whoiam   ; iw = current_o % children(ii) % npoinbox 
             case ( SEARCH_MESH_FILL   ) ; ielem = current_o % children(ii) % idfilled ; iw = current_o % children(ii) % npoinbox       
             end select
             write(90,*) ifath,ielem,iw ; flush(90) !; print*,'hijoooo'
          end do
       end if

    !stop

       select case ( my_criterion )
       case ( SEARCH_MESH_ALL    ) ; current_o => current_o % all_next
       case ( SEARCH_MESH_LEAVES ) ; current_o => current_o % leaf_next
       case ( SEARCH_MESH_FILL   ) ; current_o => current_o % fill_next  
       end select

    end do
    
  end subroutine graph

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   Return the centroid
  !> @details Compute the centroids of the leaves
  !> 
  !-----------------------------------------------------------------------

  subroutine centroid(octree,coorc,MEMORY_COUNTER,OFFSET_IELEM,CRITERION)

    class(maths_tree),                    intent(in)    :: octree
    real(rp),                      pointer, intent(inout) :: coorc(:,:)
    integer(8),          optional,          intent(inout) :: MEMORY_COUNTER(2)
    integer(ip),         optional,          intent(in)    :: OFFSET_IELEM
    integer(ip),         optional,          intent(in)    :: CRITERION      !< If empty bins should be considered
    type(octbox),                  pointer                :: current_o
    integer(ip)                                           :: ielem,divmax,ndime
    integer(ip)                                           :: nelem,mnode,idime
    integer(ip)                                           :: npoin
    logical(lg)                                           :: conti
    integer(8)                                            :: memor_loc(2)
    integer(ip)                                           :: offset_ielem_loc
    real(rp)                                              :: xx(3,8),rnode
    logical(lg)                                           :: if_empty
    integer(ip)                                           :: my_criterion 
    
    my_criterion = optional_argument(SEARCH_MESH_LEAVES , CRITERION)
    if_empty     = my_criterion == SEARCH_MESH_LEAVES
    
    if( present(MEMORY_COUNTER) ) then
       memor_loc = MEMORY_COUNTER
    else
       memor_loc = 0_8
    end if
    if( present(OFFSET_IELEM) ) then
       offset_ielem_loc = OFFSET_IELEM
    else
       offset_ielem_loc = 0
    end if

    call octree % mesh_dim(ndime,mnode,nelem,npoin,CRITERION)    
    divmax    =  octree % divmax
    current_o => octree % tree_root
    rnode     =  1.0_rp / real(mnode,rp)
    !
    ! Allocate mesh
    !
    if( .not. associated(coorc) ) call memory_alloca(memor_loc,'COORC','mesh',coorc,ndime,nelem)
    !
    ! Mesh arrays
    !    
    conti =  .true.

    loop_conti2: do while( conti )
       !
       ! First go to deepest level in first branch
       !
       do while( current_o % whoiam == 0 )
          current_o => current_o % children(1)
       end do

       if( if_empty .and. current_o % whoiam /= 0 ) then
          ielem = current_o % whoiam + offset_ielem_loc
       else if( (.not. if_empty) .and. current_o % idfilled /= 0 ) then
          ielem = current_o % idfilled + offset_ielem_loc
       else
          ielem = 0
       end if

       if( ielem /= 0 ) then
          !
          ! Current bin is a leaf
          !
          if( ndime == 2 ) then
             xx(1:2,1) = (/ current_o % minc(1),current_o % minc(2) /)
             xx(1:2,2) = (/ current_o % maxc(1),current_o % minc(2) /)
             xx(1:2,3) = (/ current_o % maxc(1),current_o % maxc(2) /)
             xx(1:2,4) = (/ current_o % minc(1),current_o % maxc(2) /)
          else
             xx(1:3,1) = (/ current_o % minc(1),current_o % minc(2),current_o % minc(3) /)
             xx(1:3,2) = (/ current_o % maxc(1),current_o % minc(2),current_o % minc(3) /)
             xx(1:3,3) = (/ current_o % maxc(1),current_o % maxc(2),current_o % minc(3) /)
             xx(1:3,4) = (/ current_o % minc(1),current_o % maxc(2),current_o % minc(3) /)
             xx(1:3,5) = (/ current_o % minc(1),current_o % minc(2),current_o % maxc(3) /)
             xx(1:3,6) = (/ current_o % maxc(1),current_o % minc(2),current_o % maxc(3) /)
             xx(1:3,7) = (/ current_o % maxc(1),current_o % maxc(2),current_o % maxc(3) /)
             xx(1:3,8) = (/ current_o % minc(1),current_o % maxc(2),current_o % maxc(3) /)
          end if
          do idime = 1,ndime
             coorc(idime,ielem) = sum(xx(idime,1:mnode))  
          end do
          coorc(:,ielem) = coorc(:,ielem) * rnode

       end if

       if(current_o % childid < divmax .and. current_o % childid /=0 ) then
          !
          ! I'm not the last child neither the Padrino
          !
          current_o => current_o % parent % children(current_o % childid+1)

       else if( current_o % childid == divmax ) then
          !
          ! I'm the last child of this generation: postprocess 
          !
          do while(current_o % id > 0 )
             if(current_o % parent % id == 0) then
                conti = .false. 
                exit loop_conti2
             else
                if(current_o % parent % childid /=divmax) then
                   current_o => current_o % parent % parent % children(current_o % parent % childid+1)
                   exit
                else 
                   current_o => current_o % parent
                end if
             end if
          end do

       else if( current_o % id == 0 ) then
          !
          ! I'm the Padrino
          !
          conti = .false.

       end if

    end do loop_conti2

    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc

  end subroutine centroid

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-10
  !> @brief   Get the global ID of the bin
  !> @details Given a point, get the global ID of the bin it
  !>          is located in
  !> 
  !-----------------------------------------------------------------------

  integer(ip) function id(octree,coord) 

    class(maths_tree),      intent(in) :: octree
    real(rp),                 intent(in) :: coord(:) !< Coordinate of the test point
    type(octbox),     pointer            :: current_o
    integer(ip)                          :: ichild,divmax

    current_o  => octree % tree_root
    id         =  0_ip
    divmax     =  octree % divmax
    
    if( octree % dim == 3 ) then

       do while( current_o % whoiam == 0 )    
          childloop3: do ichild = 1,divmax         
             if(    coord(1) >= current_o % children(ichild) % minc(1) .and. &
                  & coord(1) <= current_o % children(ichild) % maxc(1) .and. &
                  & coord(2) >= current_o % children(ichild) % minc(2) .and. &
                  & coord(2) <= current_o % children(ichild) % maxc(2) .and. &
                  & coord(3) >= current_o % children(ichild) % minc(3) .and. &
                  & coord(3) <= current_o % children(ichild) % maxc(3) ) then
                current_o => current_o % children(ichild)
                id        =  current_o % whoiam
                exit childloop3
             end if
          end do childloop3
       end do

    else if( octree % dim == 2 ) then

       do while( current_o % whoiam == 0 )  
          childloop4: do ichild = 1,divmax
             if(    coord(1) >= current_o % children(ichild) % minc(1) .and. &
                  & coord(1) <= current_o % children(ichild) % maxc(1) .and. &
                  & coord(2) >= current_o % children(ichild) % minc(2) .and. &
                  & coord(2) <= current_o % children(ichild) % maxc(2)  ) then
                current_o => current_o % children(ichild)
                id        =  current_o % whoiam
                exit childloop4
             end if
          end do childloop4
       end do

    end if

  end function id

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-10
  !> @brief   Get the global ID of the bin
  !> @details Given a point, get the global ID of the bin it
  !>          is located in
  !> 
  !-----------------------------------------------------------------------

  function host_bin(octree,coord) result(bin)

    class(maths_tree),        intent(in) :: octree
    real(rp),                 intent(in) :: coord(:) !< Coordinate of the test point
    type(octbox),     pointer            :: bin
    type(octbox),     pointer            :: bin_tmp
    integer(ip)                          :: ichild,divmax,idime
    logical(lg)                          :: ifoun 

    if( associated(octree % tree_root) ) then

       bin    => octree % tree_root
       ifoun  =  .false.
       divmax =  octree % divmax

       do idime = 1,octree % dim
          if( coord(idime) < bin % minc(idime) .or. coord(idime) > bin % maxc(idime) ) then
             ifoun = .false.
             nullify(bin)
             return
          end if
       end do

       if( bin % whoiam /= 0 ) then
          !
          ! Root is a leaf!
          !
          if( octree % dim == 3 ) then
             if(    coord(1) >= bin % minc(1) .and. &
                  & coord(1) <= bin % maxc(1) .and. &
                  & coord(2) >= bin % minc(2) .and. &
                  & coord(2) <= bin % maxc(2) .and. &
                  & coord(3) >= bin % minc(3) .and. &
                  & coord(3) <= bin % maxc(3) ) then
                ifoun = .true.
             end if
          else
             if(    coord(1) >= bin % minc(1) .and. &
                  & coord(1) <= bin % maxc(1) .and. &
                  & coord(2) >= bin % minc(2) .and. &
                  & coord(2) <= bin % maxc(2) ) then
                ifoun = .true.
             end if
          end if

       else
          !
          ! Go through the tree
          !
          if( octree % dim == 3 ) then

             do while( bin % whoiam == 0 )    
                childloop3: do ichild = 1,divmax   
                   if(    coord(1) >= bin % children(ichild) % minc(1) .and. &
                        & coord(1) <= bin % children(ichild) % maxc(1) .and. &
                        & coord(2) >= bin % children(ichild) % minc(2) .and. &
                        & coord(2) <= bin % children(ichild) % maxc(2) .and. &
                        & coord(3) >= bin % children(ichild) % minc(3) .and. &
                        & coord(3) <= bin % children(ichild) % maxc(3) ) then
                      bin_tmp  => bin % children(ichild)
                      bin      => bin_tmp
                      ifoun    = .true.
                      exit childloop3
                   end if
                end do childloop3
             end do

          else if( octree % dim == 2 ) then

             do while( bin % whoiam == 0 )
                childloop2: do ichild = 1,divmax

                   if(    coord(1) >= bin % children(ichild) % minc(1) .and. &
                        & coord(1) <= bin % children(ichild) % maxc(1) .and. &
                        & coord(2) >= bin % children(ichild) % minc(2) .and. &
                        & coord(2) <= bin % children(ichild) % maxc(2)  ) then
                      bin_tmp  => bin % children(ichild)
                      bin      => bin_tmp
                      ifoun    = .true.
                      exit childloop2
                   end if
                end do childloop2
             end do

          end if

       end if

       if( .not. ifoun ) nullify(bin)

    else

       nullify(bin)
       
    end if

  end function host_bin

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-10
  !> @brief   Results
  !> @details Get some results
  !> 
  !-----------------------------------------------------------------------

  subroutine results(self,xx,names,OFFSET,MEMORY_COUNTER,CRITERION,ONLY_DEALLOCATE)

    class(maths_tree),                   intent(inout) :: self
    real(rp),                   pointer, intent(inout) :: xx(:,:)
    character(len=5),           pointer, intent(inout) :: names(:)
    integer(ip),      optional,          intent(in)    :: OFFSET
    integer(8),       optional,          intent(inout) :: MEMORY_COUNTER(2)
    integer(ip),      optional,          intent(in)    :: CRITERION     !< If empty bins should be considered
    logical(lg),      optional,          intent(in)    :: ONLY_DEALLOCATE
    integer(ip)                                        :: nelem
    integer(ip)                                        :: ndime,mnode,npoin
    integer(ip)                                        :: ielem
    integer(ip)                                        :: ielem_offset
    integer(8)                                         :: memor_loc(2)
    type(octbox),               pointer                :: current_o
    logical(lg)                                        :: if_deallocate
    integer(ip)                                        :: my_criterion

    call self % mem_ini(memor_loc,MEMORY_COUNTER)

    if_deallocate = optional_argument(.false.,ONLY_DEALLOCATE)
    ielem_offset  = optional_argument(0_ip   ,OFFSET)
    my_criterion  = optional_argument(SEARCH_MESH_LEAVES , CRITERION)

    if( if_deallocate ) then
       !
       ! Deallocate
       !
       call memory_deallo(memor_loc,'XX'   ,'results',xx)
       call memory_deallo(memor_loc,'NAMES','results',names)

    else
       call self % mesh_dim(ndime,mnode,nelem,npoin,CRITERION)
       current_o => self % tree_root
       !
       ! Deallocate results and look for next leaf
       !
       if( .not. associated(xx)    ) then
          call memory_alloca(memor_loc,'XX'   ,'results',xx,nelem,3_ip)
       else
          if( size(xx,2) /= 3 ) call runend('DEF_MATHS_TREE: WRONG XX SIZE')
       end if
       
       if( .not. associated(names) ) then
          call memory_alloca(memor_loc,'NAMES','results',names,3_ip)
       else
          if( size(names) /= 3 ) call runend('DEF_MATHS_TREE: WRONG NAMES SIZE')
       end if

       select case ( my_criterion )
       case ( SEARCH_MESH_ALL    ) ; current_o => self % all_root
       case ( SEARCH_MESH_LEAVES ) ; current_o => self % leaf_root
       case ( SEARCH_MESH_FILL   ) ; current_o => self % fill_root
       end select

       do while( associated(current_o) )
         
          select case ( my_criterion )
          case ( SEARCH_MESH_ALL    ) ; ielem = current_o % id       + ielem_offset + 1
          case ( SEARCH_MESH_LEAVES ) ; ielem = current_o % whoiam   + ielem_offset
          case ( SEARCH_MESH_FILL   ) ; ielem = current_o % idfilled + ielem_offset   
          end select

          names(1)    =  'NPOIN'
          names(2)    =  'ID'
          names(3)    =  'LEVEL'
          xx(ielem,1) =  real(current_o % npoinbox,rp)
          xx(ielem,2) =  real(current_o % id      ,rp)
          xx(ielem,3) =  real(current_o % level   ,rp)

          select case ( my_criterion )
          case ( SEARCH_MESH_ALL    ) ; current_o => current_o % all_next
          case ( SEARCH_MESH_LEAVES ) ; current_o => current_o % leaf_next
          case ( SEARCH_MESH_FILL   ) ; current_o => current_o % fill_next  
          end select

       end do

    end if

    call self % mem_end(memor_loc,MEMORY_COUNTER)

  end subroutine results

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Nearest bins
  !> @details Find the nearest bins
  !> 
  !-----------------------------------------------------------------------


  
  subroutine near(self,coord,lcheck,nn,nn_tot)

    class(maths_tree),         intent(in)    :: self
    real(rp),                  intent(in)    :: coord(:)    !< Coordinate of the test point
    type(stack),      pointer, intent(inout) :: lcheck(:)
    integer(ip),               intent(out)   :: nn          !< Number of candidate boxes
    integer(ip),               intent(out)   :: nn_tot      !< Number of candidate boxes
    integer(ip)                              :: nd
    integer(ip)                              :: divmax,istack
    integer(ip)                              :: ncheck
    integer(ip)                              :: nstack,ii,jj,ll
    integer(ip)                              :: nlist,ncount
    real(rp)                                 :: dimax,dd
    real(rp)                                 :: dimin
    real(rp)                                 :: dismm(2)
    real(rp)                                 :: bobox(2,3)
    type(octbox),     pointer                :: current_o
    integer(ip),      pointer                :: list(:)
    real(rp),         pointer                :: mindi(:,:)
    type(stack),      pointer                :: lstack(:)

    nullify(lstack)
    nullify(mindi)
    nullify(list)


    nd        =  self % dim
    nn        =  0
    nn_tot    =  0
    dimin     =  big
    dimax     =  big
    divmax    =  self % divmax
    current_o => self % tree_root
    !
    ! Initialize check list
    !
    ncheck = 0
    nstack = self % nboxes
    if( nstack > 0 ) then
       !
       ! Initialize stack
       !     
       !allocate(lcheck(nstack))
       allocate(lstack(nstack))
       allocate(mindi(2,divmax))
       allocate(list(divmax))
       istack             =  1
       lstack(istack) % p => current_o

       !
       ! Traverse tree for searching BB candidates
       !
       ncount = 0
       do while( istack > 0 ) 
          current_o => lstack(istack) % p
          !
          !
          !
          if( current_o % whoiam == 0 ) then
             !
             ! Update upper bound dimax. Only if the maximum distancce between the given point
             ! and the current BB of the tree node is smaller than the current upper bound
             !
             bobox(1,1:nd)  =  current_o % minc(1:nd)
             bobox(2,1:nd)  =  current_o % maxc(1:nd)
             dismm          =  maths_min_max_box_vertices(coord(1:nd),bobox(1:2,1:nd))

             istack = istack - 1
             if( dismm(1) <= dimax + epsil100) then

                if( dismm(2) < dimax) then
                   dimax =  dismm(2)             
                end if
                !
                ! Determine BB minimum distances of the children
                !
                !istack = istack - 1
                do ii = 1,divmax
                   bobox(1,1:nd)  =  current_o % children(ii) % minc(1:nd)
                   bobox(2,1:nd)  =  current_o % children(ii) % maxc(1:nd)
                   dismm          =  maths_min_max_box_vertices(coord(1:nd),bobox(1:2,1:nd))
                   mindi(:,ii)    =  dismm(:)
                   ncount         =  ncount + 1

                   !if( dismm(1) <= dimin .and. ncheck == 0 ) then ! CRISTOBAL
                   !   dimin       =  dismm(1)
                   !   dimax       =  dismm(2)
                   !end if
                end do
                !
                ! Select only children which minimum distancce between the given point and its BB 
                ! is smaller than the current upper_bound dimax             
                !
                nlist = 0
                do ii = 1,divmax
                   !   if( mindi(1,ii) <= dimax + epsil100 ) then
                   nlist           =  nlist + 1
                   mindi(1,nlist)  =  mindi(1,ii) 
                   list(nlist)     =  ii
                   !   end if
                end do
                !
                ! Sort selected children into the stack according to their minimum distance to the given point
                !
                do ii = 2,nlist
                   dd =  mindi(1,ii)
                   ll =  list(ii)
                   jj =  ii - 1
                   do while( jj >= 1 )
                      if( mindi(1,jj) - epsil100 > dd ) exit 
                      mindi(1,jj+1) = mindi(1,jj)
                      list(jj+1)    = list(jj)
                      jj            = jj - 1
                   end do
                   mindi(1,jj+1) = dd
                   list(jj+1)    = ll
                end do
                !
                ! Update list of candidates
                !
                do ii = 1,nlist
                   istack             =  istack + 1
                   lstack(istack) % p => current_o % children(list(ii))
                end do
             end if
          else if( current_o % idfilled > 0 ) then ! It is a leaf of the tree
             !
             ! Reset min and max distances 
             !
             !if( ncheck == 0 ) then
             !   dimin = big
             !   dimax = big
             !end if
             !
             ! Update upper bound dimax. Only if the maximum distancce between the given point
             ! and the current BB of the tree node is smaller than the current upper bound
             !
             ncount         =  ncount + 1
             bobox(1,1:nd)  =  current_o % minc(1:nd)
             bobox(2,1:nd)  =  current_o % maxc(1:nd)
             dismm          =  maths_min_max_box_vertices(coord(1:nd),bobox(1:2,1:nd))

             if( dismm(1) <= dimax + epsil100) then

                if( dismm(2) < dimax ) then
                   dimax =  dismm(2)
                end if
                !
                ! Save the current BB if its minimum distance between the given point
                ! and its BB is smaller than the current upper_bound dimax
                !                          
                !if( dismm(1) <= dimax + epsil100 ) then               
                ncheck             =  ncheck + 1
                lcheck(ncheck) % p => current_o
                nn                 =  nn + 1
                nn_tot             =  nn_tot + memory_size(lcheck(nn) % p % nodes)                
                !if( dismm(1) < dimin ) then
                !   dimin           =  dismm(1)
                !end if
             end if
             istack = istack - 1
          else
             istack = istack - 1
          end if

       end do

       !print*,'MIN/MAX=',ncheck,dimin,dimax
       !
       ! Percentage saved with respect to inquiring all leaves
       !print*,'NCOUNT=',real(ncount,rp)/real(self % nfilled,rp)*100.0_rp
       !
       ! Redefine candidates using max distance
       !
       !nn = 0
       !do icheck = 1,ncheck
       !   bobox(1,1:nd)  =  lcheck(icheck) % p % minc(1:nd)
       !   bobox(2,1:nd)  =  lcheck(icheck) % p % maxc(1:nd)
       !   dismm(1)       =  maths_min_box_vertices(coord(1:nd),bobox(1:2,1:nd))
       !   if( dismm(1) <= dimax ) then
       !      nn             =  nn + 1
       !      lcheck(nn) % p => lcheck(icheck) % p
       !      nn_tot         =  nn_tot + memory_size(lcheck(nn) % p % nodes)
       !   end if
       !end do
       !
       ! Deallocate
       !    
       if( associated(list)   ) deallocate(list  ) 
       if( associated(mindi)  ) deallocate(mindi ) 
       if( associated(lstack) ) deallocate(lstack)

    end if

  end subroutine near  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Nearest bins
  !> @details Find the nearest bins
  !> 
  !-----------------------------------------------------------------------

  pure subroutine divide(self,current_o,sorted_points,coord,centr)

    class(maths_tree),                     intent(in)    :: self
    type(octbox),                 pointer, intent(inout) :: current_o
    integer(ip),                  pointer, intent(inout) :: sorted_points(:) ! For KD-TREE        
    real(rp),           optional, pointer, intent(in)    :: coord(:,:)
    real(rp),           optional, pointer, intent(in)    :: centr(:,:)        !< Center of gravity
    real(rp)                                             :: delta_x, delta_y, delta_Z
    real(rp)                                             :: values(current_o % npoinbox)        
    integer(ip)                                          :: ipoin,kpoin,kdime
    logical(lg)                                          :: divide_z, divide_y
    
    select case ( self % method )
       
    case ( SEARCH_OCTREE ) 
       !
       ! OCTREE
       !
       current_o % children(1) % minc(1) = current_o % minc(1)
       current_o % children(1) % minc(2) = current_o % minc(2)
       current_o % children(1) % maxc(1) = (current_o % maxc(1) + current_o % minc(1))*0.5_rp
       current_o % children(1) % maxc(2) = (current_o % maxc(2) + current_o % minc(2))*0.5_rp

       current_o % children(2) % minc(1) = (current_o % maxc(1) + current_o % minc(1))*0.5_rp
       current_o % children(2) % minc(2) = current_o % children(1) % minc(2)
       current_o % children(2) % maxc(1) = current_o % maxc(1)     
       current_o % children(2) % maxc(2) = current_o % children(1) % maxc(2)

       current_o % children(3) % minc(1) = current_o % children(1) % minc(1)
       current_o % children(3) % minc(2) = (current_o % minc(2) + current_o % maxc(2))*0.5_rp
       current_o % children(3) % maxc(1) = current_o % children(1) % maxc(1)
       current_o % children(3) % maxc(2) = current_o % maxc(2)

       current_o % children(4) % minc(1) = current_o % children(2) % minc(1)
       current_o % children(4) % minc(2) = current_o % children(3) % minc(2)
       current_o % children(4) % maxc(1) = current_o % children(2) % maxc(1)
       current_o % children(4) % maxc(2) = current_o % children(3) % maxc(2)

       if( self % dim == 3 ) then

          current_o % children(1) % minc(3) = current_o % minc(3)
          current_o % children(1) % maxc(3) = (current_o % maxc(3) + current_o % minc(3))*0.5_rp

          current_o % children(2) % minc(3) = current_o % children(1) % minc(3)
          current_o % children(2) % maxc(3) = current_o % children(1) % maxc(3)
          current_o % children(3) % minc(3) = current_o % children(1) % minc(3)
          current_o % children(3) % maxc(3) = current_o % children(1) % maxc(3)
          current_o % children(4) % minc(3) = current_o % children(1) % minc(3)
          current_o % children(4) % maxc(3) = current_o % children(1) % maxc(3)

          current_o % children(5) % minc(1) = current_o % children(1) % minc(1)
          current_o % children(5) % minc(2) = current_o % children(1) % minc(2)
          current_o % children(5) % minc(3) = (current_o % minc(3) + current_o % maxc(3))*0.5_rp
          current_o % children(5) % maxc(1) = current_o % children(1) % maxc(1)
          current_o % children(5) % maxc(2) = current_o % children(1) % maxc(2)
          current_o % children(5) % maxc(3) = current_o % maxc(3)

          current_o % children(6) % minc(1) = current_o % children(2) % minc(1)
          current_o % children(6) % minc(2) = current_o % children(1) % minc(2)
          current_o % children(6) % minc(3) = current_o % children(5) % minc(3)
          current_o % children(6) % maxc(1) = current_o % children(2) % maxc(1)     
          current_o % children(6) % maxc(2) = current_o % children(1) % maxc(2)
          current_o % children(6) % maxc(3) = current_o % children(5) % maxc(3)

          current_o % children(7) % minc(1) = current_o % children(1) % minc(1)
          current_o % children(7) % minc(2) = current_o % children(3) % minc(2)
          current_o % children(7) % minc(3) = current_o % children(5) % minc(3)
          current_o % children(7) % maxc(1) = current_o % children(1) % maxc(1)
          current_o % children(7) % maxc(2) = current_o % children(3) % maxc(2)
          current_o % children(7) % maxc(3) = current_o % children(5) % maxc(3)

          current_o % children(8) % minc(1) = current_o % children(2) % minc(1)
          current_o % children(8) % minc(2) = current_o % children(3) % minc(2)
          current_o % children(8) % minc(3) = current_o % children(5) % minc(3)
          current_o % children(8) % maxc(1) = current_o % children(2) % maxc(1)
          current_o % children(8) % maxc(2) = current_o % children(3) % maxc(2)
          current_o % children(8) % maxc(3) = current_o % children(5) % maxc(3)

       end if

    case ( SEARCH_KDTREE  )
       !
       ! KDTREE
       !
       kdime = 1              
       if (self % unsorted_division) then 
          if(      mod(current_o % level + 1_ip,3_ip) == 0 .and. self % dim == 3 ) then
             kdime = 3
          else if( mod(current_o % level + 1_ip,2_ip) == 0 ) then
             kdime = 2
          end if
       else
          delta_x = current_o % maxc(1) - current_o % minc(1) 
          delta_y = current_o % maxc(2) - current_o % minc(2) 
          delta_z = current_o % maxc(self % dim) - current_o % minc(self % dim)           
          if ((delta_z >= delta_x) .and. (delta_z >= delta_y) .and. self % dim == 3) then
             kdime    = 3 
          elseif (delta_y >= delta_x) then
             kdime    = 2              
          end if
       end if
       
       select case ( kdime )
          
       case ( 1_ip )
          
          current_o % children(1) % minc(1) = current_o % minc(1)
          current_o % children(1) % minc(2) = current_o % minc(2)
          current_o % children(1) % maxc(1) = (current_o % maxc(1) + current_o % minc(1))*0.5_rp
          current_o % children(1) % maxc(2) = current_o % maxc(2)
          
          current_o % children(2) % minc(1) = current_o % children(1) % maxc(1)
          current_o % children(2) % minc(2) = current_o % minc(2)
          current_o % children(2) % maxc(1) = current_o % maxc(1)
          current_o % children(2) % maxc(2) = current_o % maxc(2)

          if(self % dim == 3 ) then
             current_o % children(1) % minc(3) = current_o % minc(3)
             current_o % children(1) % maxc(3) = current_o % maxc(3)
             current_o % children(2) % minc(3) = current_o % minc(3)
             current_o % children(2) % maxc(3) = current_o % maxc(3)
          end if
          
       case ( 2_ip )
          
          current_o % children(1) % minc(1) = current_o % minc(1)
          current_o % children(1) % minc(2) = current_o % minc(2)
          current_o % children(1) % maxc(1) = current_o % maxc(1)
          current_o % children(1) % maxc(2) = (current_o % maxc(2) + current_o % minc(2))*0.5_rp
          
          current_o % children(2) % minc(1) = current_o % minc(1)
          current_o % children(2) % minc(2) = (current_o % maxc(2) + current_o % minc(2))*0.5_rp
          current_o % children(2) % maxc(1) = current_o % maxc(1)
          current_o % children(2) % maxc(2) = current_o % maxc(2)
          
          if(self % dim == 3 ) then
             current_o % children(1) % minc(3) = current_o % minc(3)
             current_o % children(1) % maxc(3) = current_o % maxc(3)
             current_o % children(2) % minc(3) = current_o % minc(3)
             current_o % children(2) % maxc(3) = current_o % maxc(3)
          end if
          
       case ( 3_ip )
          
          current_o % children(1) % minc(1) = current_o % minc(1)
          current_o % children(1) % minc(2) = current_o % minc(2)
          current_o % children(1) % maxc(1) = current_o % maxc(1)
          current_o % children(1) % maxc(2) = current_o % maxc(2) 
          
          current_o % children(2) % minc(1) = current_o % minc(1)
          current_o % children(2) % minc(2) = current_o % maxc(2) 
          current_o % children(2) % maxc(1) = current_o % maxc(1)
          current_o % children(2) % maxc(2) = current_o % maxc(2)
          
          current_o % children(1) % minc(3) = current_o % minc(3)
          current_o % children(1) % maxc(3) = (current_o % maxc(3) + current_o % minc(3))*0.5_rp
          current_o % children(2) % minc(3) = current_o % children(1) % maxc(3)
          current_o % children(2) % maxc(3) = current_o % maxc(3)          
          
       end select

    
    case ( SEARCH_SKDTREE )
       !
       ! SKD-TREE
       !
       ! 1. Define the axis will be divided 
       !
       divide_z = .false.
       divide_y = .false.
       !
       ! In a z-y-x axis order or in the longest axis
       !
       if (self % unsorted_division) then 
          if(      mod(current_o % level + 1_ip,3_ip) == 0 .and. self % dim == 3 ) then
             divide_z = .true.
          else if( mod(current_o % level + 1_ip,2_ip) == 0 ) then
             divide_y = .true.          
          end if                    
       else
          delta_x = current_o % maxc(1) - current_o % minc(1) 
          delta_y = current_o % maxc(2) - current_o % minc(2) 
          delta_z = current_o % maxc(self % dim) - current_o % minc(self % dim) 

          kdime    = 1              
          if ((delta_z >= delta_x) .and. (delta_z >= delta_y) .and. self % dim == 3) then
             divide_z = .true.
             kdime    = 3 
          elseif (delta_y >= delta_x) then
             divide_y = .true.
             kdime    = 2              
          end if
       end if
       !
       ! 2. Store point ids and their coordinates on the selected division axis  
       !
       if(      present(coord) ) then          
          do ipoin = 1,current_o % npoinbox
             kpoin                = current_o % nodes(ipoin)
             sorted_points(ipoin) = kpoin
             values(ipoin)        = coord(kdime,kpoin)
          end do
       else if( present(centr) ) then
          do ipoin = 1,current_o % npoinbox
             kpoin                = current_o % nodes(ipoin)
             sorted_points(ipoin) = kpoin
             values(ipoin)        = centr(kdime,kpoin)             
          end do
       end if
       !
       ! Sort points using values in increasing order
       !
       if( current_o % npoinbox > 0 ) &
            call maths_heap_sort(2_ip,current_o % npoinbox,values,sorted_points)
       
    end select
    
  end subroutine divide

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    24/01/2022
  !> @brief   Search type name
  !> @details Search type name
  !
  !-----------------------------------------------------------------------

  function type_name(self) result(name)
    class(maths_tree),                       intent(inout) :: self
    character(LEN=:), allocatable                          :: name

    select type ( self )
    class is ( maths_skdtree ) ; name = 'SKD-TREE'
    class is ( maths_kdtree  ) ; name = 'KD-TREE'
    class is ( maths_octree  ) ; name = 'OCT-TREE'
    end select
    
  end function type_name

end module def_maths_tree
!> @}
