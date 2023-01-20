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
!> @file    def_maths_bin.f90
!> @author  Cristobal Samaniego
!> @brief   Variables
!> @details Variables fot mod_maths.f90
!
!-----------------------------------------------------------------------

module def_maths_skdtree

  use def_kintyp_basic,          only : ip,rp,lg,i1p,r2p
  use def_elmtyp,                only : QUA04
  use def_elmtyp,                only : HEX08
  use mod_memory_basic,          only : memory_alloca
  use mod_memory_basic,          only : memory_deallo
  use mod_memory_basic,          only : memory_size
  use mod_maths_geometry,        only : maths_point_in_box
  use mod_maths_sort,            only : maths_heap_sort
  use def_search_method,         only : search_method
  use def_search_method,         only : SEARCH_FILL
  use def_search_method,         only : SEARCH_CANDIDATE
  use def_search_method,         only : SEARCH_BOUNDING_BOXES
  use def_search_method,         only : CANDIDATE_NEAREST 
  use def_search_method,         only : SEARCH_DEALLO
  use mod_optional_argument,     only : optional_argument

  implicit none
  
  real(rp),      parameter :: epsil         = epsilon(1.0_rp)
  real(rp),      parameter :: epsil1        = epsilon(1.0_rp)
  real(rp),      parameter :: epsil100      = 100.0_rp*epsilon(1.0_rp)
  
  !
  ! Skd tree node class
  !  
  type tree_node
     real(rp)                  :: dista_min
     real(rp)                  :: dista_max     
     real(rp)                  :: my_bobox(2,3)   ! Bounding box containing all my boundary boxes
     integer(ip),      pointer :: my_id_bobox(:)  ! Sorted list of the id od my boundary boxes      
     class(tree_node), pointer :: parent          ! Pointer to parent
     class(tree_node), pointer :: child_left      ! Pointer to left children
     class(tree_node), pointer :: child_right     ! Pointer to right children
     
   contains     
     procedure                 :: initialize     
     procedure                 :: determine_boundary_box
     procedure                 :: determine_min_distance
     procedure                 :: determine_max_distance          
     procedure                 :: reorder_id_bobox          
  end type tree_node
  !
  ! Class containing the skd tree
  !    
  type, extends(search_method) :: maths_skd_tree
     class(tree_node), pointer :: tree_root
     integer(ip),      pointer :: all_id_bobox(:) ! Sorted list of the id of all the boundary boxes     
   contains     
     procedure,           pass :: init            ! Initialize all
     procedure,           pass :: deallo          ! Deallocate
     procedure,           pass :: candidate       ! Candidate
     procedure,           pass :: fill            ! Fill
     procedure,           pass :: results         ! Get results on the mesh
     procedure,           pass :: mesh            ! Create a mesh
     procedure,           pass :: input           ! Input parameters
     procedure,           pass :: type_name       ! Type name
  end type maths_skd_tree
  !
  ! Stack node type
  !  
  type stack_node
     class(tree_node), pointer :: data
     type(stack_node), pointer :: next
  end type stack_node
  !
  ! Class containing the stack structure data
  !
  type stack
     type(stack_node), pointer :: first
     integer(ip)               :: len=0     
   contains
     procedure :: pop
     procedure :: push
     procedure :: clear_stack     
  end type stack
  class(stack), pointer  :: my_stack

  public :: maths_skd_tree

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  samaniego
  !> @date    2021-01-31
  !> @brief   Skdtree
  !> @details Set input data
  !> 
  !-----------------------------------------------------------------------
  
  subroutine input(self,limit,relative_tolerance,absolute_tolerance,fill_method)

    class(maths_skd_tree),       intent(inout) :: self
    integer(ip),      optional, intent(in)    :: limit
    real(rp),         optional, intent(in)    :: relative_tolerance
    real(rp),         optional, intent(in)    :: absolute_tolerance
    integer(ip),      optional, intent(in)    :: fill_method

    call self % input_all(relative_tolerance,absolute_tolerance,fill_method)

  end subroutine input

  !-----------------------------------------------------------------------
  !> 
  !> @author  samaniego
  !> @date    2021-01-28
  !> @brief   Skdtree
  !> @details Init tree root node 
  !> 
  !-----------------------------------------------------------------------
  
  subroutine init(self)
    class(maths_skd_tree), intent(inout) :: self

    call self % init_all()
    !
    ! Initialize tree node
    !
    nullify(self % tree_root)
    self % name = 'SKD-TREE'
    
  end subroutine init

  !-----------------------------------------------------------------------
  !> 
  !> @author  samaniego
  !> @date    2021-01-28
  !> @brief   Skdtree
  !> @details Deallocate skdtree memory
  !> 
  !-----------------------------------------------------------------------
  
  subroutine deallo(self,list_entities,MEMORY_COUNTER)

    class(maths_skd_tree),                   intent(inout) :: self
    type(i1p),            optional, pointer, intent(inout) :: list_entities(:)    
    integer(8),           optional,          intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                                             :: memor_loc(2)

    class(tree_node),               pointer                :: current_node
    
    call self % tim_ini()
    call self % mem_ini(memor_loc,MEMORY_COUNTER)

    self % kfl_filled =  0_ip

    if( associated(self % tree_root) ) then
       !
       ! Initialice stack, put root node to start
       !
       if( associated(my_stack) ) & 
            nullify(my_stack)
       allocate(my_stack)
       call my_stack % push(self % tree_root)
       !
       ! Loop until stack be empty
       !
       do while (my_stack % len > 0)
          current_node   => my_stack % pop()
          !
          ! It is NOT a leaf of the tree. Add children to the stack
          !
          if (size(current_node % my_id_bobox) > 1) then
             call my_stack % push(current_node % child_left)
             call my_stack % push(current_node % child_right)             
          end if
          !
          ! Deallocate current node
          !          
          deallocate(current_node)
       end do
       call my_stack % clear_stack
    end if
    nullify(self % tree_root)
    if( associated(self % all_id_bobox) ) deallocate(self % all_id_bobox) ! Do not deallocate before tree_root       
    
    if( present(list_entities) ) &
         call memory_deallo(memor_loc,'LIST_ENTITIES','def_maths_bin',list_entities)

    if( allocated(self % name) ) deallocate(self % name)
    
    call self % mem_end(memor_loc,MEMORY_COUNTER)
    call self % tim_end(SEARCH_DEALLO)

  end subroutine deallo

  !-----------------------------------------------------------------------
  !> 
  !> @author  samaniego
  !> @date    2021-01-28
  !> @brief   Skdtree
  !> @details Fill skd tree from a list of boundary boxes that correspond to
  !>          faces of a surface mesh
  !> 
  !-----------------------------------------------------------------------
  
  subroutine fill(self, coord,bobox,MEMORY_COUNTER,COORD_MIN,COORD_MAX,PERMU,MASK)

    class(maths_skd_tree),           intent(inout) :: self

    real(rp),    optional, pointer, intent(in)    :: coord(:,:)
    real(rp),    optional, pointer, intent(in)    :: bobox(:,:,:)      !< Bounding boxes
    integer(8),  optional,          intent(inout) :: MEMORY_COUNTER(2)
    real(rp),    optional,          intent(in)    :: COORD_MIN(:)
    real(rp),    optional,          intent(in)    :: COORD_MAX(:)
    integer(ip), optional, pointer, intent(in)    :: PERMU(:)
    logical(lg), optional, pointer, intent(in)    :: MASK(:)

    real(rp),              pointer                :: centr(:,:)
    integer(ip)                                   :: iboun, idime, nboun, current_nboun, new_nboun, ndime
    integer(8)                                    :: memor_loc(2)

    class(tree_node),      pointer                :: current_node

    call self % tim_ini()
    call self % mem_ini(memor_loc,MEMORY_COUNTER)

    self % kfl_filled = 1
    
    if( present(bobox) ) then
       
       self % fill_method = SEARCH_BOUNDING_BOXES
       !
       ! 0. Initializations
       !       
       ndime        = memory_size(bobox,2_ip)
       nboun        = memory_size(bobox,3_ip)
       self % dim   = ndime
       self % nelem = nboun

       if( nboun <= 0 ) return
       !
       ! Root tree
       !
       allocate(self % tree_root)    
       call self % tree_root % initialize()
       nullify(self % all_id_bobox)
       
       nullify(centr)  
       call memory_alloca(memor_loc,'CENTR'       ,'fill',centr,ndime,nboun)
       
       if( associated(self % all_id_bobox) ) deallocate(self % all_id_bobox)       
       call memory_alloca(memor_loc,'ALL_ID_BOBOX','fill',self % all_id_bobox,nboun)
       !
       ! 0.1. Calculate the centroids and save the id of all BBs
       !
       do iboun = 1,nboun
          ! Check if a boundary box dimension is equal to 0
          ! It can happens for a boundary box of a surface that is parallel to an axis
          do idime = 1,ndime          
             if (bobox(2,idime,iboun)-bobox(1,idime,iboun) <= epsil100) then
                bobox(1,idime,iboun) = bobox(1,idime,iboun) - epsil100
                bobox(2,idime,iboun) = bobox(2,idime,iboun) + epsil100
             end if
          end do
          self % all_id_bobox(iboun) = iboun           
          centr(1:ndime,iboun) = (bobox(1,1:ndime,iboun)+bobox(2,1:ndime,iboun))*0.5_rp
       end do
       !
       ! 0.2. Initiliaze root node
       !
       self % tree_root % my_id_bobox => self % all_id_bobox(:)
       call self % tree_root % determine_boundary_box(ndime, bobox)
       !
       ! 1. Build tree
       !
       ! 1.1. Initialice stack, put root node to start
       !
       if( associated(my_stack) ) nullify(my_stack)
       allocate(my_stack)
       call my_stack % push(self % tree_root)
       !
       ! 1.2 Loop until stack be empty
       !
       do while (my_stack % len > 0)
          !
          ! 1.3. Pop an item form the stack
          !
          current_node   => my_stack % pop()
          current_nboun  =  memory_size(current_node % my_id_bobox)
          
          if (current_nboun > 1) then             
             !
             ! 1.4. Create and initilize children
             !             
             allocate( current_node % child_left)
             call current_node % child_left  % initialize()
             allocate( current_node % child_right)       
             call current_node % child_right % initialize()
             !
             ! 1.5. Choose larger dimension and divide equally the BBs bewtween the children
             !
             call current_node % reorder_id_bobox(ndime, nboun, current_nboun, centr)
             !new_nboun = int(real(current_nboun,rp)/2.0_rp,ip)
             new_nboun = current_nboun/2_ip
             !
             ! 1.6. Put children into the stack
             !
             current_node      % child_right % parent      => current_node                          
             current_node      % child_right % my_id_bobox => current_node % my_id_bobox(new_nboun+1:current_nboun)
             call current_node % child_right % determine_boundary_box(ndime, bobox)             
             call my_stack     % push(current_node % child_right)
             
             current_node      % child_left % parent      => current_node                          
             current_node      % child_left % my_id_bobox => current_node % my_id_bobox(1:new_nboun)
             call current_node % child_left % determine_boundary_box(ndime, bobox)             
             call my_stack     % push(current_node % child_left)
          end if
       end do
       
       call my_stack % clear_stack
       deallocate(my_stack)

       deallocate(centr)
    else
       call runend('DEF_MATHS_SKD_TREE: SKD TREE DOES NOT WORK WITHOUT BOUNDARY BOXES')
    end if

    self % kfl_filled = 2
    call self % tim_end(SEARCH_FILL)
    call self % mem_end(memor_loc,MEMORY_COUNTER)

  end subroutine fill

  !-----------------------------------------------------------------------
  !> 
  !> @author  samaniego
  !> @date    2021-01-28
  !> @brief   Skdtree
  !> @details Obtain a list boundary box candidates closest to a
  !>          set of given points
  !> 
  !-----------------------------------------------------------------------
  
  subroutine candidate(self,xx,list_entities,METHOD,MASK,MEMORY_COUNTER)

    class(maths_skd_tree),               intent(inout) :: self
    real(rp),                            intent(in)    :: xx(:,:)           !< List coordinates
    type(i1p),                  pointer, intent(inout) :: list_entities(:)    
    integer(8),       optional,          intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(ip),      optional,          intent(in)    :: METHOD            !< Method for candidates
    logical(lg),      optional, pointer, intent(in)    :: MASK(:)           !< Mask to consider or not points
    integer(ip),                pointer                :: current_list(:)
    integer(8)                                         :: memor_loc(2)
    integer(ip)                                        :: inode, ndime, nnode, iboun, nboun, current_ncand    
    real(rp)                                           :: xcoor(3), upper_bound
    class(tree_node),           pointer                :: current_node

    call self % tim_ini()
    call self % mem_ini(memor_loc,MEMORY_COUNTER)

    if( present(METHOD) ) then
       if( METHOD /= CANDIDATE_NEAREST ) then
          call runend('SKDTREE: WRONG METHOD FOR SKD_TREE')
          stop 2
       end if
    end if
    
    ndime = int(size(xx,1_ip),ip)
    nnode = int(size(xx,2_ip),ip)
    nboun = memory_size(self % tree_root % my_id_bobox)    

    if( nboun > 0 ) then
       
       if( .not. associated(list_entities) ) & 
            call memory_alloca(memor_loc,'LIST_ENTITIES','maths_skd_tree',list_entities,nnode)

       allocate(current_list(nboun))

       !
       ! Loop over input points
       !    
       do inode = 1,nnode

          if( associated(my_stack) ) & 
               nullify(my_stack)
          allocate(my_stack)    
          call my_stack % push(self % tree_root)
          
          xcoor          = 0.0_rp       
          xcoor(1:ndime) = xx(1:ndime,inode)

          call  self % tree_root % determine_max_distance(ndime, xcoor)

          upper_bound   = self % tree_root % dista_max
          current_ncand = 0_ip
          !
          ! Traverse tree for searching BB candidates
          !
          do while (my_stack % len > 0)
             current_node  => my_stack % pop()
             call  current_node % determine_max_distance(ndime, xcoor)          
             !
             ! Check if minimum distacnce between a point and the current BB of the tree node
             ! is smaller than the current upper_bound  
             !
             if (current_node % dista_min <= upper_bound + epsil100) then
                !upper_bound = min(upper_bound,current_node % dista_max)

                if (current_node % dista_max < upper_bound) then
                   upper_bound = current_node % dista_max
                end if
                !
                ! It is a leaf of the tree, save the id BB is a candidate
                !
                if (size(current_node % my_id_bobox) == 1) then
                   current_ncand               = current_ncand + 1
                   current_list(current_ncand) = current_node % my_id_bobox(1)                
                else
                   !
                   ! Put children into the stack according to their minimum distance to the given point
                   !
                   call  current_node % child_left  % determine_min_distance(ndime, xcoor)
                   call  current_node % child_right % determine_min_distance(ndime, xcoor)
                   if  (current_node % child_right % dista_min < current_node % child_left % dista_min - epsil100) then
                      call my_stack % push(current_node % child_left)
                      call my_stack % push(current_node % child_right)
                   else
                      call my_stack % push(current_node % child_right)
                      call my_stack % push(current_node % child_left)
                   end if
                end if
             end if

          end do

          call memory_alloca(memor_loc,'LIST_ENTITIES % L','maths_skd_tree',list_entities(inode)%l,current_ncand)
          do iboun = 1,current_ncand
             list_entities(inode) % l(iboun) = current_list(iboun)
          end do
          
          call my_stack % clear_stack
          deallocate(my_stack)          
       end do

       deallocate(current_list)
    end if

    call self % mem_end(memor_loc,MEMORY_COUNTER)
    call self % tim_end(SEARCH_CANDIDATE)

  end subroutine candidate
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  samaniego
  !> @date    2021-01-28
  !> @brief   Skdtree
  !> @details Mesh
  !> 
  !-----------------------------------------------------------------------
  
  subroutine mesh(self,ndime,mnode,nelem,npoin,lnods,ltype,coord,MEMORY_COUNTER,&
       OFFSET_IELEM,OFFSET_IPOIN,CRITERION,CENTROID,ONLY_DEALLOCATE)

    class(maths_skd_tree),                  intent(inout)  :: self
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

  end subroutine mesh

  !-----------------------------------------------------------------------
  !> 
  !> @author  samaniego
  !> @date    2021-01-28
  !> @brief   Skdtree
  !> @details Results
  !> 
  !-----------------------------------------------------------------------
  
  subroutine results(self,xx,names,OFFSET,MEMORY_COUNTER,CRITERION,ONLY_DEALLOCATE)

    class(maths_skd_tree),                   intent(inout) :: self
    real(rp),                       pointer, intent(inout) :: xx(:,:)
    character(len=5),               pointer, intent(inout) :: names(:)
    integer(ip), optional,                   intent(in)    :: OFFSET
    integer(8),  optional,                   intent(inout) :: MEMORY_COUNTER(2)
    integer(ip), optional,                   intent(in)    :: CRITERION       !< If empty bins should be considered
    logical(lg), optional,                   intent(in)    :: ONLY_DEALLOCATE
    
  end subroutine results


  !-----------------------------------------------------------------------
  !> 
  !> @author  samaniego
  !> @date    2021-02-02
  !> @brief   Skdtree node
  !> @details Determine the boundary box for a given node on the skd tree
  !> 
  !-----------------------------------------------------------------------
  
  pure subroutine determine_boundary_box(self, ndime, bobox)
    
    class(tree_node),          intent(inout) :: self
    integer(ip),               intent(in)    :: ndime
    real(rp),         pointer, intent(in)    :: bobox(:,:,:)      !< Bounding boxes    
    integer(ip)                              :: iitem, nitems, idime, iboun

    nitems = memory_size(self % my_id_bobox)
    do iitem = 1,nitems
       iboun = self % my_id_bobox(iitem)
       do idime = 1,ndime       
          self % my_bobox(1,idime) = min( bobox(1,idime,iboun) , self % my_bobox(1,idime) ) - epsil1 
          self % my_bobox(2,idime) = max( bobox(2,idime,iboun) , self % my_bobox(2,idime) ) + epsil1
       end do
    end do
    
  end subroutine determine_boundary_box

  !-----------------------------------------------------------------------
  !> 
  !> @author  samaniego
  !> @date    2021-02-02
  !> @brief   Skdtree node
  !> @details Initialiaze a node of the skd tree
  !> 
  !-----------------------------------------------------------------------
  
  pure subroutine initialize(self)
    
    class(tree_node), intent(inout) :: self
 
    self % my_bobox = reshape((/1.0e10_rp, -1.0e10_rp, 1.0e10_rp, -1.0e10_rp, 1.0e10_rp, -1.0e10_rp/), &
         (/ 2, 3/))
    
    self % dista_min = 0.0_rp
    self % dista_max = 0.0_rp
    
    nullify(self % my_id_bobox)
    nullify(self % parent     )
    nullify(self % child_left )
    nullify(self % child_right)

  end subroutine initialize

  !-----------------------------------------------------------------------
  !> 
  !> @author  samaniego
  !> @date    2021-02-04
  !> @brief   Skdtree node
  !> @details Reorder the id of BBs according to their centroid values
  !> 
  !-----------------------------------------------------------------------
  
  pure subroutine reorder_id_bobox(self, ndime, nboun, current_nboun, centr)
    
    class(tree_node),          intent(inout) :: self
    integer(ip),               intent(in)    :: ndime
    integer(ip),               intent(in)    :: nboun        
    integer(ip),               intent(in)    :: current_nboun    
    real(rp),                  intent(in)    :: centr(ndime,nboun)    
    integer(ip)                              :: iitem, iboun, max_idime
    real(rp)                                 :: values(current_nboun)
    real(rp)                                 :: delta_x, delta_y, delta_Z    
    !
    ! Choose the largest dimension
    !
    delta_x = self % my_bobox(2,1)     - self % my_bobox(1,1) 
    delta_y = self % my_bobox(2,2)     - self % my_bobox(1,2)
    delta_z = self % my_bobox(2,ndime) - self % my_bobox(1,ndime)
   
    max_idime = 1_ip
    if ((delta_z >= delta_x) .and. (delta_z >= delta_y) .and. ndime == 3) then
       max_idime = 3 
    elseif (delta_y >= delta_x) then
       max_idime = 2              
    end if    
    !
    ! Reorder the id of BBs according to their centroid values
    !
    do iitem = 1, current_nboun
       iboun         = self % my_id_bobox(iitem)              
       values(iitem) = centr(max_idime,iboun)
    end do

    call maths_heap_sort(2_ip,current_nboun,values,IVO1=self % my_id_bobox)

  end subroutine reorder_id_bobox

  !-----------------------------------------------------------------------
  !> 
  !> @author  samaniego
  !> @date    2021-02-02
  !> @brief   Skdtree node
  !> @details Find maximum distance from a point to a BB of a tree node
  !> 
  !-----------------------------------------------------------------------
  
  pure subroutine determine_max_distance(self,ndime,xcoor)
    
    class(tree_node),          intent(inout) :: self
    integer(ip),               intent(in)    :: ndime    
    real(rp),                  intent(in)    :: xcoor(3)    
    real(rp)                                 :: temp    
    integer(ip)                              :: idime

    self % dista_max = 0.0_rp
    do idime = 1,ndime
       temp             = max(xcoor(idime) - self % my_bobox(1,idime), self % my_bobox(2,idime) - xcoor(idime))
       self % dista_max = self % dista_max + temp*temp
    end do
    self % dista_max = sqrt(self % dista_max)
    
  end subroutine determine_max_distance

  !-----------------------------------------------------------------------
  !> 
  !> @author  samaniego
  !> @date    2021-02-02
  !> @brief   Skdtree node
  !> @details Find the minimum distance from a point to a BB of a tree node
  !> 
  !-----------------------------------------------------------------------
  
  pure subroutine determine_min_distance(self,ndime,xcoor)
    
    class(tree_node),          intent(inout) :: self
    integer(ip),               intent(in)    :: ndime    
    real(rp),                  intent(in)    :: xcoor(3)    
    real(rp)                                 :: temp    
    integer(ip)                              :: idime

    self % dista_min = 0.0_rp
    do idime = 1,ndime
       temp             = max(0.0_rp, self % my_bobox(1,idime) - xcoor(idime))
       temp             = max(temp, xcoor(idime) - self % my_bobox(2,idime))
       self % dista_min = self % dista_min + temp*temp
    end do
    self % dista_min = sqrt(self % dista_min)
    
  end subroutine determine_min_distance

  !-----------------------------------------------------------------------
  !> 
  !> @author  samaniego
  !> @date    2021-02-04
  !> @brief   Stack
  !> @details Pop an item from the stack
  !> 
  !-----------------------------------------------------------------------
  
  function pop(this) result(x)
    
    class(stack)              :: this
    class(tree_node), pointer :: x
    type(stack_node), pointer :: tmp
    
    if ( this%len == 0 ) then
       print*, "popping from empty stack"
       !stop
    end if
    tmp => this%first
    x   => this%first%data
    this%first => this%first%next
    deallocate(tmp)
    this%len = this%len -1
    
  end function pop

  !-----------------------------------------------------------------------
  !> 
  !> @author  samaniego
  !> @date    2021-02-04
  !> @brief   Stack
  !> @details Push an item into the stack
  !> 
  !-----------------------------------------------------------------------
  
  subroutine push(this, x)
    
    class(tree_node), pointer :: x
    class(stack)              :: this
    type(stack_node), pointer :: new, tmp
    
    allocate(new)
    new%data => x
    if (.not. associated(this%first)) then
       this%first => new
    else
       tmp => this%first
       this%first => new
       this%first%next => tmp
    end if
    this%len = this%len + 1
    
  end subroutine push

  !-----------------------------------------------------------------------
  !> 
  !> @author  samaniego
  !> @date    2021-01-31
  !> @brief   Stack
  !> @details Clear stack
  !> 
  !-----------------------------------------------------------------------
  
  subroutine clear_stack(this)
    
    class(stack)              :: this
    type(stack_node), pointer :: tmp
    integer(ip)               :: i
    
    if ( this%len == 0 ) then
       return
    end if
    do i = 1, this%len
       tmp => this%first
       if ( .not. associated(tmp)) exit
       this%first => this%first%next
       deallocate(tmp)
    end do
    this%len = 0
    
  end subroutine clear_stack
  
  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    24/01/2022
  !> @brief   Search type name
  !> @details Search type name
  !
  !-----------------------------------------------------------------------

  function type_name(self) result(name)
    
    class(maths_skd_tree),                   intent(inout) :: self
    character(LEN=:), allocatable                          :: name

    name = 'SKD-TREE'
    
  end function type_name

end module def_maths_skdtree
