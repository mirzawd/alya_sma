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
!> @brief   Bin
!> @details Binary tree
!>
!>          fill
!>          ----
!>          self % stats(1) ... Average number of elements per bin
!>          self % stats(4) ... Max number of elements in a bin
!>          self % times(1) ... Total fill time
!>
!>          candidate
!>          ---------
!>          self % stats(5) ... Saving ratio: average number of candidates
!>                              divided by total number of elements
!>          self % times(2) ... Total candidate time
!>
!>          others
!>          ------
!>          self % times(3) ... Total deallocate time
!>
!>          INPUT:
!>          ------
!>          if limit is given as input, the the nuber of boxes in computed 
!>          automatically as:
!>          self % boxip(1)   = ( npoin / limit ) ** (1/ndime)
!>
!>
!-----------------------------------------------------------------------

module def_maths_bin

  use def_kintyp_basic,          only : ip,rp,lg,i1p,r2p
  use def_elmtyp,                only : QUA04
  use def_elmtyp,                only : HEX08
  use mod_memory_basic,          only : memory_alloca
  use mod_memory_basic,          only : memory_deallo
  use mod_memory_basic,          only : memory_size
  use mod_memory_basic,          only : memory_resize
  use mod_maths_geometry,        only : maths_point_in_box
  use mod_maths_geometry,        only : maths_min_box_vertices
  use mod_maths_geometry,        only : maths_max_box_vertices
  use mod_maths_geometry,        only : maths_min_max_box_vertices
  use def_search_method,         only : search_method
  use def_search_method,         only : SEARCH_FILL
  use def_search_method,         only : SEARCH_CANDIDATE
  use def_search_method,         only : SEARCH_DEALLO
  use def_search_method,         only : CANDIDATE_INSIDE
  use def_search_method,         only : CANDIDATE_NEAREST
  use def_search_method,         only : SEARCH_POINTS
  use def_search_method,         only : SEARCH_BOUNDING_BOXES
  use def_search_method,         only : SEARCH_MESH_LEAVES
  use def_search_method,         only : SEARCH_MESH_FILL
  use def_search_method,         only : SEARCH_MESH_ALL
  use def_search_method,         only : big
  use mod_optional_argument,     only : optional_argument
  use mod_std

  real(rp),    parameter :: epsil = epsilon(1.0_rp)
  integer(ip), parameter :: LINKED_LIST_BIN = 0_ip
  integer(ip), parameter :: TYPE_BIN        = 1_ip

  type, extends(search_method) ::  maths_bin
     integer(ip)             :: boxip(3)             ! Number of boxes in integer(ip) (input)
     integer(ip)             :: boxes                ! Total number of boxes 
     type(i1p),   pointer    :: list(:,:,:)          ! List of elements
     integer(ip), pointer    :: num(:,:,:)           ! Number of elements
     integer(ip), pointer    :: ia(:)                ! Linked list ia
     integer(ip), pointer    :: ja(:)                ! Linked list ja
     integer(ip)             :: nfilled              ! Number of filled bins
     integer(ip)             :: fill_type            ! Type of fill, linked list or type
     integer(ip)             :: max_bins             ! Maximum number of bins
     real(rp)                :: boxrp(3)             ! Number of boxes in real(rp)
     real(rp)                :: delta(3)             ! Bin size 
     real(rp)                :: delti(3)             ! Inverse of bin size 
   contains
     procedure,         pass :: init                 ! Initialize all
     procedure,         pass :: deallo               ! Deallocate
     procedure,         pass :: candidate            ! Return the list of candidates
     procedure,         pass :: fill                 ! Allocate 
     procedure,         pass :: input                ! Set input data
     procedure,         pass :: results              ! Get results on the mesh
     procedure,         pass :: mesh                 ! Create a bin mesh
     procedure,         pass :: type_name            ! Name 

     procedure,         pass :: find                 ! Find a bin given a coordinate
     procedure,         pass :: near                 ! Find nearset bins
     procedure,         pass :: find_bb              ! Find a bin given a bounding box
     procedure,         pass :: inbin                ! If a point is inside bounding box of a bin
     procedure,         pass :: box_num              ! Box number: 3D to 1D coordiantes
     procedure,         pass :: num_points           ! Number of points in box
     procedure,         pass :: filled               ! If a box is filled or empty
     procedure,         pass :: coord_1              ! Bounding box of a bin (1 argument)
     procedure,         pass :: coord_2              ! Bounding box of a bin (2 arguments)
     generic                 :: coord => coord_1,&
          &                              coord_2
     procedure,         pass :: mesh_dim             ! Dimension of the bin mesh
     procedure,         pass :: mesh_node            ! Node number of the mesh
     procedure,         pass :: mesh_coord           ! Coordinates of bin nodes
     procedure,         pass :: centroid             ! Centroid    
     procedure,         pass :: search               ! Search a series of points  
  end type maths_bin
  
  private

  public :: maths_bin
  public :: LINKED_LIST_BIN
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Start timing
  !> @details End timing
  !> 
  !-----------------------------------------------------------------------

  subroutine fin(self)

    type(maths_bin), intent(inout) :: self
    
    call self % deallo()

  end subroutine fin
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Input
  !> @details Set input data
  !> 
  !-----------------------------------------------------------------------

  subroutine input(self,boxes,relative_tolerance,absolute_tolerance,fill_method,name,limit,max_bins,PARAM)
    
    class(maths_bin),           intent(inout) :: self
    integer(ip),      optional, intent(in)    :: boxes(:)
    real(rp),         optional, intent(in)    :: relative_tolerance
    real(rp),         optional, intent(in)    :: absolute_tolerance
    integer(ip),      optional, intent(in)    :: fill_method
    character(len=*), optional, intent(in)    :: name
    integer(ip),      optional, intent(in)    :: limit
    integer(ip),      optional, intent(in)    :: max_bins
    real(rp),         optional, intent(in)    :: PARAM(:)
    integer(ip)                               :: ii

    if( present(PARAM) ) then
       self % boxip(1:3) = int(param(1:3),ip)
    end if
    
    if( present(boxes) ) then
       do ii = 1,size(boxes,KIND=ip)
          self % boxip(ii) = boxes(ii)
       end do       
    end if

    if( present(limit) ) then
       self % boxip = -limit
    end if

    if( present(max_bins) ) then
       self % max_bins = max_bins 
    end if
    
    call self % input_all(relative_tolerance,absolute_tolerance,fill_method,name)

  end subroutine input

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Search
  !> @details Search in parallel
  !> 
  !-----------------------------------------------------------------------

  subroutine search(self,xx,list_entities,MEMORY_COUNTER)

    class(maths_bin),                    intent(inout) :: self
    real(rp),                            intent(in)    :: xx(:,:)           !< List coordinates
    type(i1p),                  pointer, intent(inout) :: list_entities(:)    
    integer(8),       optional,          intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(ip)                                        :: ii,nn,kk,ll,ib,iz
    integer(8)                                         :: memor_loc(2)
    integer(ip)                                        :: xloc(3)
    
    call self % mem_ini(memor_loc,MEMORY_COUNTER)
    
    nn = size(xx,2_ip,KIND=ip)
 
    if( .not. associated(list_entities) ) & 
         call memory_alloca(memor_loc,'LIST_ENTITIES','maths_bin',list_entities,nn)
    
    do ii = 1,nn
       call self % find(xx(:,ii),xloc)
       ll = self % num_points(xloc)
       call memory_alloca(memor_loc,'LIST_ENTITIES % L','maths_bin',list_entities(ii)%l,ll)
       if( xloc(1) > 0 ) then
          if( self % fill_type == LINKED_LIST_BIN ) then
             ib = self % box_num(xloc)
             kk = 0
             do iz = self % ia(ib),self % ia(ib+1)-1
                kk = kk + 1
                list_entities(ii) % l(kk) = self % ja(iz)
             end do
          else
             do kk = 1,ll
                list_entities(ii) % l(kk) = self % list(xloc(1),xloc(2),xloc(3)) % l(kk)
             end do
          end if
       end if
    end do
    
    call self % mem_end(memor_loc,MEMORY_COUNTER)
    
  end subroutine search
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Search
  !> @details Candidate elements
  !> 
  !-----------------------------------------------------------------------

  subroutine candidate(self,xx,list_entities,METHOD,MASK,MEMORY_COUNTER)

    class(maths_bin),                    intent(inout) :: self
    real(rp),                            intent(in)    :: xx(:,:)           !< List coordinates
    type(i1p),                  pointer, intent(inout) :: list_entities(:)    
    integer(8),       optional,          intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(ip),      optional,          intent(in)    :: METHOD                !< Method for candidates
    logical(lg),      optional, pointer, intent(in)    :: MASK(:)           !< Mask to consider or not points
    integer(ip)                                        :: ii,jj,nn
    integer(ip)                                        :: kk,ll,ienti
    integer(ip)                                        :: iz,ib
    integer(8)                                         :: memor_loc(2)
    integer(ip)                                        :: xloc(3),xloc_num
    integer(ip),                pointer                :: xloc_list(:,:)
    logical(lg)                                        :: mask_loc
    integer(ip)                                        :: my_method

    call self % tim_ini()
    call self % mem_ini(memor_loc,MEMORY_COUNTER)

    my_method = optional_argument(CANDIDATE_INSIDE,METHOD)
    nn        = size(xx,2_ip,KIND=ip)
    mask_loc  = .true.
    nullify(xloc_list)

    if( nn > 0 ) then

       if( .not. associated(list_entities) ) then
          call memory_alloca(memor_loc,'LIST_ENTITIES','maths_bin',list_entities,nn)          
       end if

       if( self % boxes > 0 ) then
            
          if( my_method == CANDIDATE_INSIDE ) then
             !
             ! Inside
             !
             do ii = 1,nn
                if( present(MASK) ) then
                   if( associated(MASK) ) mask_loc = MASK(ii)   
                end if
                if( mask_loc ) then
                   ll = 0
                   call self % find(xx(:,ii),xloc)
                   if( xloc(1) > 0 ) then
                      if( self % fill_type == LINKED_LIST_BIN ) then
                         ib = self % box_num(xloc)
                         ll = self % ia(ib+1)-self % ia(ib)
                         call memory_alloca(memor_loc,'LIST_ENTITIES % L','maths_bin',list_entities(ii)%l,ll)
                         kk = 0
                         do iz = self % ia(ib),self % ia(ib+1)-1
                            kk = kk + 1
                            list_entities(ii) % l(kk) = self % ja(iz)
                         end do
                      else
                         ll = self % num(xloc(1),xloc(2),xloc(3))
                         call memory_alloca(memor_loc,'LIST_ENTITIES % L','maths_bin',list_entities(ii)%l,ll)
                         do kk = 1,ll
                            list_entities(ii) % l(kk) = self % list(xloc(1),xloc(2),xloc(3)) % l(kk)
                         end do
                      end if
                   end if
                   self % stats(5) = self % stats(5) + real(ll,rp)               
                end if
             end do

          else if( my_method == CANDIDATE_NEAREST ) then
             !
             ! Nearest
             !
             call memory_alloca(memor_loc,'XLOC_LIST','maths_bin',xloc_list,3_ip,self % boxes)
             do ii = 1,nn
                if( present(MASK) ) then
                   if( associated(MASK) ) mask_loc = MASK(ii)
                end if
                if( mask_loc ) then
                   call self % near(xx(:,ii),xloc_list,xloc_num)
                   ll  = 0
                   if( self % fill_type == LINKED_LIST_BIN ) then
                      do jj = 1,xloc_num
                         ib = self % box_num(xloc_list(:,jj)) 
                         ll = ll + self % ia(ib+1)-self % ia(ib)
                      end do
                      self % stats(5) = self % stats(5) + real(ll,rp)                
                      call memory_alloca(memor_loc,'LIST_ENTITIES % L','maths_bin',list_entities(ii)%l,ll)
                      ll = 0
                      do jj = 1,xloc_num
                         ib = self % box_num(xloc_list(:,jj)) 
                         do iz = self % ia(ib),self % ia(ib+1)-1
                            ienti = self % ja(iz)
                            if( .not. any( list_entities(ii) % l == ienti ) ) then
                               ll = ll + 1
                               list_entities(ii) % l(ll) = ienti
                            end if
                         end do
                      end do
                   else
                      do jj = 1,xloc_num
                         ll = ll + self % num(xloc_list(1,jj),xloc_list(2,jj),xloc_list(3,jj))
                      end do
                      self % stats(5) = self % stats(5) + real(ll,rp)                
                      call memory_alloca(memor_loc,'LIST_ENTITIES % L','maths_bin',list_entities(ii)%l,ll)
                      ll = 0
                      do jj = 1,xloc_num
                         do kk = 1,self % num(xloc_list(1,jj),xloc_list(2,jj),xloc_list(3,jj))
                            ienti = self % list(xloc_list(1,jj),xloc_list(2,jj),xloc_list(3,jj)) % l(kk)
                            if( .not. any( list_entities(ii) % l == ienti ) ) then
                               ll = ll + 1
                               list_entities(ii) % l(ll) = ienti
                            end if
                         end do
                      end do
                   end if
                   call memory_resize(memor_loc,'LIST_ENTITIES % L','maths_bin',  list_entities(ii) % l,ll)              
                end if
             end do
             call memory_deallo(memor_loc,'XLOC_LIST','maths_bin',xloc_list)
          end if
    
          self % stats(5) = self % stats(5) / (real(max(1_ip,self % nelem),rp)*real(nn,rp))
       end if
    end if

    call self % mem_end(memor_loc,MEMORY_COUNTER)
    call self % tim_end(SEARCH_CANDIDATE)

  end subroutine candidate
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Bin init
  !> @details Initialization of bin
  !> 
  !-----------------------------------------------------------------------

  subroutine init(self)

    class(maths_bin), intent(inout) :: self

    call self % init_all()

    self % name       = 'BIN'
    self % nfilled    =  0
    self % fill_type  =  LINKED_LIST_BIN ! TYPE_BIN
    self % boxip      =  1_ip
    self % boxrp      =  1.0_rp
    self % boxes      =  0_ip
    self % delta      =  1.0_rp
    self % delti      =  1.0_rp
    self % memor      =  0_8
    self % max_bins   =  int(1e7,ip)
    !
    ! Pinters for type-based fill in
    !
    nullify(self % list)    
    nullify(self % num )
    !
    ! Pinters to linked-list-based fill in
    !
    nullify(self % ia  )    
    nullify(self % ja  )    

  end subroutine init

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Bin
  !> @details Fill in a bin from a coordinate field
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo(self,list_entities,MEMORY_COUNTER)

    class(maths_bin),                     intent(inout) :: self
    type(i1p),         optional, pointer, intent(inout) :: list_entities(:)    
    integer(8),        optional,          intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                                          :: memor_loc(2)

    call self % tim_ini()
    call self % mem_ini(memor_loc,MEMORY_COUNTER)
    
    self % kfl_filled =  0_ip

    call memory_deallo(memor_loc,'SELF % NUM' ,'maths_bin',self % num )
    call memory_deallo(memor_loc,'SELF % LIST','maths_bin',self % list)
    call memory_deallo(memor_loc,'SELF % IA'  ,'maths_bin',self % ia  )
    call memory_deallo(memor_loc,'SELF % JA'  ,'maths_bin',self % ja  )
   
    if( allocated(self % name) ) deallocate(self % name)

    if( present(list_entities) ) call memory_deallo(memor_loc,'LIST_ENTITIES','def_maths_bin',list_entities)
    
    call self % mem_end(memor_loc,MEMORY_COUNTER)
    call self % tim_end(SEARCH_DEALLO)

  end subroutine deallo

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   Return the centroid
  !> @details Compute the centroids of the leaves
  !> 
  !-----------------------------------------------------------------------
  
  subroutine centroid(bin,coorc,MEMORY_COUNTER,OFFSET_IELEM,CRITERION)

    class(maths_bin),                       intent(inout) :: bin   
    real(rp),                      pointer, intent(inout) :: coorc(:,:)
    integer(8),          optional,          intent(inout) :: MEMORY_COUNTER(2)
    integer(ip),         optional,          intent(in)    :: OFFSET_IELEM
    integer(ip),         optional,          intent(in)    :: CRITERION     !< If empty bins should be considered
    integer(ip)                                           :: mnode,ndime
    integer(ip)                                           :: nelem,inode
    integer(ip)                                           :: npoin,ielem
    integer(ip)                                           :: ii,jj,kk
    real(rp)                                              :: rnode,xx(3)
    integer(8)                                            :: memor_loc(2)
    real(rp)                                              :: coord(3,8)
    integer(ip)                                           :: my_criterion

    my_criterion = optional_argument(SEARCH_MESH_LEAVES , CRITERION)
    
    call bin % mem_ini(memor_loc,MEMORY_COUNTER)
    !
    ! Dimensions
    !
    call bin % mesh_dim(NDIME=ndime,MNODE=mnode,NELEM=nelem,NPOIN=npoin,CRITERION=my_criterion)
    !
    ! Allocate mesh
    !
    if( .not. associated(coorc) ) call memory_alloca(memor_loc,'COORC','mesh',coorc,ndime,nelem)
    !
    ! Centroid
    !
    rnode = 1.0_rp / real(mnode,rp)
    ielem = 0
    if( my_criterion == SEARCH_MESH_LEAVES ) then
       do kk = 1,bin % boxip(3)
          do jj = 1,bin % boxip(2)
             do ii = 1,bin % boxip(1)
                call bin % mesh_coord((/ii,jj,kk/),coord)
                xx = 0.0_rp
                do inode = 1,mnode
                   xx(1:ndime) = xx(1:ndime) + coord(1:ndime,inode)
                end do
                ielem = ielem + 1
                coorc(1:ndime,ielem) = xx(1:ndime) * rnode
             end do
          end do
       end do
    else
       do kk = 1,bin % boxip(3)
          do jj = 1,bin % boxip(2)
             do ii = 1,bin % boxip(1)
                if( associated(bin % list(ii,jj,kk) % l) ) then
                   call bin % mesh_coord((/ii,jj,kk/),coord)
                   xx = 0.0_rp
                   do inode = 1,mnode
                      xx(1:ndime) = xx(1:ndime) + coord(1:ndime,inode)
                   end do
                   ielem = ielem + 1
                   coorc(1:ndime,ielem) = xx(1:ndime) * rnode
                end if
             end do
          end do
       end do       
    end if

    call bin % mem_end(memor_loc,MEMORY_COUNTER)
    
  end subroutine centroid
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Bin
  !> @details Fill in a bin from a coordinate field
  !> 
  !-----------------------------------------------------------------------

  subroutine fill(self,coord,bobox,MEMORY_COUNTER,COORD_MIN,COORD_MAX,PERMU,MASK)

    class(maths_bin),                        intent(inout) :: self              !< Bin structure
    real(rp),    optional,          pointer, intent(in)    :: coord(:,:)        !< Coordinates
    real(rp),    optional,          pointer, intent(in)    :: bobox(:,:,:)      !< Bounding boxes
    integer(8),  optional,                   intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    real(rp),    optional,                   intent(in)    :: COORD_MIN(:)
    real(rp),    optional,                   intent(in)    :: COORD_MAX(:)
    integer(ip), optional,          pointer, intent(in)    :: PERMU(:)
    logical(lg), optional,          pointer, intent(in)    :: MASK(:)
    integer(ip)                                            :: ii,jj,kk
    integer(ip)                                            :: npoin_ini,npoin_end
    integer(ip)                                            :: ndime,nz,iz,ipoin,np
    integer(ip)                                            :: nn,nb,ib
    real(rp)                                               :: rr
    integer(ip)                                            :: xx(3)
    integer(ip),           pointer                         :: xx_bb(:,:)
    integer(8)                                             :: memor_loc(2)
    integer(ip)                                            :: imeth
    logical(lg)                                            :: if_coord
    logical(lg)                                            :: if_bobox
    integer(ip),           pointer                         :: na(:)
    
    nullify(xx_bb,na)
    self % kfl_filled = 1

    call self % mem_ini(memor_loc,MEMORY_COUNTER)
    call self % tim_ini()
    !
    ! Determine method
    !
    if_coord = .false.
    if_bobox = .false.
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
       np    = memory_size(coord,2_ip)
    else if( if_bobox ) then
       imeth = 2
       np    = memory_size(bobox,3_ip)
    else
       return
    end if
    !
    ! Space dimensions
    !
    npoin_ini =  0
    npoin_end = -1
    if(      present(coord) .and. imeth == 1 ) then
       if( associated(coord) ) then
          npoin_ini  = lbound(coord,DIM=2_ip,KIND=ip)
          npoin_end  = ubound(coord,DIM=2_ip,KIND=ip)
          self % dim = memory_size(coord,1_ip)
       end if
    else if( present(bobox) .and. imeth == 2 ) then
       if( associated(bobox) ) then
          npoin_ini  = lbound(bobox,DIM=3_ip,KIND=ip)
          npoin_end  = ubound(bobox,DIM=3_ip,KIND=ip)
          self % dim = memory_size(bobox,2_ip)
       end if
    end if
    self % nelem      = npoin_end-npoin_ini+1
    ndime             = self % dim
    !
    ! Automatic boxes
    !
    if( self % boxip(1) < 0 ) then
       self % boxip(1)   = int(( real(np,rp) / real(-self % boxip(1),rp) ) ** (1.0_rp/real(ndime,rp)),ip)
       self % boxip(2:3) = self % boxip(1)
    end if
    if( ndime == 2 ) self % boxip(3) = 1
    !
    ! Limit number of bins
    !
    self % boxes = self % boxip(1) * self % boxip(2) * self % boxip(3) 
    if( self % boxes > self % max_bins ) then
       rr           = ( real(self % max_bins,rp) / real(self % boxes,rp) ) ** (1.0_rp/real(ndime,rp))
       self % boxip = int( real(self % boxip,rp) * rr ,ip)
    end if
    !
    ! Some extra variables
    !
    self % boxip(1)   = max (self % boxip(1)  ,1_ip)
    self % boxip(2)   = max (self % boxip(2)  ,1_ip)
    self % boxip(3)   = max (self % boxip(3)  ,1_ip)
    self % boxrp(1:3) = real(self % boxip(1:3),rp  )
    self % boxes      = self % boxip(1)*self % boxip(2)*self % boxip(3)
    
    if( imeth == 1 .and. present(coord) ) then 

       !-----------------------------------------------------------------
       !
       ! Fill in bin with coordinates
       !
       !-----------------------------------------------------------------

       if( .not. associated(self % list) ) then
          !
          ! Min and max coordinates
          !
          call self % bounding_box(coord,self % toler_rel)
          self % delta   = (self % comax-self % comin) / self % boxrp
          self % delti   = self % boxrp / (self % comax-self % comin)
          !
          ! Count bin elements
          !
          if( self % fill_type == LINKED_LIST_BIN ) then
             nb = self % boxes
             nz = 0
             call memory_alloca(memor_loc,'SELF % IA','maths_bin',self % ia,nb+1_ip)
             call memory_alloca(memor_loc,'NA'       ,'maths_bin',na       ,nb)
             do ipoin = npoin_ini,npoin_end
                if( optional_argument(.true.,MASK,IPOIN) ) then
                   call self % find(coord(:,ipoin),xx)
                   if( xx(1) > 0 ) then
                      nz            = nz + 1
                      ib            = self % box_num(xx)
                      self % ia(ib) = self % ia(ib) + 1
                   end if
                end if
             end do
             call size_to_linked_list(self % ia)
             do ii = 1,nb
                na(ii) = self % ia(ii) - 1
             end do
          else
             call memory_alloca(memor_loc,'SELF % NUM' ,'maths_bin',self % num ,self % boxip(1),self % boxip(2),self % boxip(3))
             call memory_alloca(memor_loc,'SELF % LIST','maths_bin',self % list,self % boxip(1),self % boxip(2),self % boxip(3))
             do ipoin = npoin_ini,npoin_end
                if( optional_argument(.true.,MASK,IPOIN) ) then
                   call self % find(coord(:,ipoin),xx)
                   if( xx(1) > 0 ) self % num(xx(1),xx(2),xx(3)) = self % num(xx(1),xx(2),xx(3)) + 1
                end if
             end do
          end if
          !
          ! Allocate bin
          !
          self % nfilled = 0
          if( self % fill_type == LINKED_LIST_BIN ) then
             call memory_alloca(memor_loc,'SELF % JA','maths_bin',self % ja,nz)
             do ib = 1,nb
                self % nfilled = self % nfilled + min(1_ip,self % ia(ib+1)-self % ia(ib))
             end do
          else
             do kk = 1,self % boxip(3)
                do jj = 1,self % boxip(2)
                   do ii = 1,self % boxip(1)
                      if( self % num(ii,jj,kk) /= 0 ) self % nfilled = self % nfilled + 1
                      call memory_alloca(memor_loc,'SELF % LIST % L','maths_bin',self % list(ii,jj,kk) % l,self % num(ii,jj,kk))
                      self % num(ii,jj,kk) = 0
                   end do
                end do
             end do
          end if
          !
          ! Fill in bin
          !
          if( self % fill_type == LINKED_LIST_BIN ) then
             do ipoin = npoin_ini,npoin_end
                if( optional_argument(.true.,MASK,IPOIN) ) then
                   call self % find(coord(:,ipoin),xx)
                   if( xx(1) > 0 ) then
                      iz     = self % box_num(xx)
                      na(iz) = na(iz) + 1
                      if( present(permu) ) then
                         self % ja(na(iz)) = permu(ipoin)                   
                      else
                         self % ja(na(iz)) = ipoin
                      end if
                   end if
                end if
             end do
             call memory_deallo(memor_loc,'NA','maths_bin',na)
          else
             do ipoin = npoin_ini,npoin_end
                if( optional_argument(.true.,MASK,IPOIN) ) then
                   call self % find(coord(:,ipoin),xx)
                   if( xx(1) > 0 ) then
                      self % num(xx(1),xx(2),xx(3)) = self % num(xx(1),xx(2),xx(3)) + 1
                      if( present(permu) ) then
                         self % list(xx(1),xx(2),xx(3)) % l(self % num(xx(1),xx(2),xx(3))) = permu(ipoin)
                      else
                         self % list(xx(1),xx(2),xx(3)) % l(self % num(xx(1),xx(2),xx(3))) = ipoin
                      end if
                   end if
                end if
             end do
          end if
          
       end if

    else if( imeth == 2 .and. present(bobox) ) then 

       !-----------------------------------------------------------------
       !
       ! Fill in bin with bounding boxes
       !
       !-----------------------------------------------------------------       
       !
       ! Min and max coordinates
       !       
       call self % bounding_box(bobox,self % toler_rel)
       self % delta = (self % comax-self % comin) / self % boxrp
       self % delti = self % boxrp / (self % comax-self % comin)

       if( .not. associated(self % list) ) then
          !
          ! Count bin elements
          !
          call memory_alloca(memor_loc,'XX_BB','maths_bin',xx_bb,3_ip,self % boxes+1_ip)
          
          if( self % fill_type == LINKED_LIST_BIN ) then
             nb = self % boxes
             nz = 0
             call memory_alloca(memor_loc,'SELF % IA','maths_bin',self % ia,nb+1_ip)
             call memory_alloca(memor_loc,'NA'       ,'maths_bin',na       ,nb)
             do ipoin = npoin_ini,npoin_end
                if( optional_argument(.true.,MASK,IPOIN) ) then
                   call self % find_bb(bobox(:,:,ipoin),xx_bb,nn)
                   nz = nz + nn
                   do ii = 1,nn 
                      ib            = self % box_num(xx_bb(:,ii))
                      self % ia(ib) = self % ia(ib) + 1
                   end do
                end if
             end do
             call size_to_linked_list(self % ia)
             do ii = 1,nb
                na(ii) = self % ia(ii) - 1
             end do
          else
             call memory_alloca(memor_loc,'SELF % NUM' ,'maths_bin',self % num ,self % boxip(1),self % boxip(2),self % boxip(3))          
             call memory_alloca(memor_loc,'SELF % LIST','maths_bin',self % list,self % boxip(1),self % boxip(2),self % boxip(3))          
             do ipoin = npoin_ini,npoin_end
                if( optional_argument(.true.,MASK,IPOIN) ) then
                   call self % find_bb(bobox(:,:,ipoin),xx_bb,nn)
                   do ii = 1,nn                
                      self % num(xx_bb(1,ii),xx_bb(2,ii),xx_bb(3,ii)) = self % num(xx_bb(1,ii),xx_bb(2,ii),xx_bb(3,ii)) + 1
                   end do
                end if
             end do
          end if
          !
          ! Allocate bin
          !
          self % nfilled = 0
          if( self % fill_type == LINKED_LIST_BIN ) then
             call memory_alloca(memor_loc,'SELF % JA','maths_bin',self % ja,nz)
             do ib = 1,nb
                self % nfilled = self % nfilled + min(1_ip,self % ia(ib+1)-self % ia(ib))
             end do
          else
             do kk = 1,self % boxip(3)
                do jj = 1,self % boxip(2)
                   do ii = 1,self % boxip(1)
                      if( self % num(ii,jj,kk) /= 0 ) self % nfilled = self % nfilled + 1
                      call memory_alloca(memor_loc,'SELF % LIST % L','maths_bin',self % list(ii,jj,kk) % l,self % num(ii,jj,kk))
                      self % num(ii,jj,kk) = 0
                   end do
                end do
             end do
          end if
          !
          ! Fill in bin
          !
          if( self % fill_type == LINKED_LIST_BIN ) then
             do ipoin = npoin_ini,npoin_end
                if( optional_argument(.true.,MASK,IPOIN) ) then
                   call self % find_bb(bobox(:,:,ipoin),xx_bb,nn)                   
                   do ii = 1,nn
                      iz     = self % box_num(xx_bb(:,ii))
                      na(iz) = na(iz) + 1
                      if( present(permu) ) then
                         self % ja(na(iz)) = permu(ipoin)                   
                      else
                         self % ja(na(iz)) = ipoin
                      end if
                   end do
                end if
             end do
             call memory_deallo(memor_loc,'NA','maths_bin',na)
          else
             do ipoin = npoin_ini,npoin_end
                if( optional_argument(.true.,MASK,IPOIN) ) then
                   call self % find_bb(bobox(:,:,ipoin),xx_bb,nn)
                   do ii = 1,nn                
                      self % num(xx_bb(1,ii),xx_bb(2,ii),xx_bb(3,ii)) = self % num(xx_bb(1,ii),xx_bb(2,ii),xx_bb(3,ii)) + 1             
                      if( present(permu) ) then
                         self % list(xx_bb(1,ii),xx_bb(2,ii),xx_bb(3,ii)) % l(self % num(xx_bb(1,ii),xx_bb(2,ii),xx_bb(3,ii))) = permu(ipoin)                   
                      else
                         self % list(xx_bb(1,ii),xx_bb(2,ii),xx_bb(3,ii)) % l(self % num(xx_bb(1,ii),xx_bb(2,ii),xx_bb(3,ii))) = ipoin
                      end if
                   end do
                end if
             end do
          end if
          call memory_deallo(memor_loc,'XX_BB'  ,'maths_bin',xx_bb)

!!$          block
!!$            use def_master, only : kfl_paral
!!$            if( kfl_paral == 0 ) then
!!$               do kk = 1,self % boxip(3)
!!$                  do jj = 1,self % boxip(2)
!!$                     do ii = 1,self % boxip(1)
!!$                        if( associated(self % list(ii,jj,kk) % l) ) &
!!$                             print*,'inside=',self % list(ii,jj,kk) % l
!!$                     end do
!!$                  end do
!!$               end do
!!$            end if
!!$          end block
!!$          
       end if

    end if
    !
    ! Compute some statistics
    !
    if( self % fill_type == LINKED_LIST_BIN ) then
       self % stats(1) = real(nz,rp)
       do ii = 1,self % boxes
          self % stats(4) = max(self % stats(4),real(self % ia(ii+1)-self % ia(ii),rp))
       end do
     else
       do kk = 1,self % boxip(3)
          do jj = 1,self % boxip(2)
             do ii = 1,self % boxip(1)
                 self % stats(1) = self % stats(1)   + real(self % num(ii,jj,kk),rp)
                 self % stats(4) = max(self % stats(4),real(self % num(ii,jj,kk),rp))
             end do
          end do
       end do
    end if
    self % stats(1) = self % stats(1) / max(real(self % boxes,rp),1.0_rp)
    self % kfl_filled = 2

    call self % mem_end(memor_loc,MEMORY_COUNTER)
    call self % tim_end(SEARCH_FILL)

  end subroutine fill

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   Create a mesh
  !> @details Create a mesh from an octbin
  !> 
  !-----------------------------------------------------------------------
  
  subroutine mesh_dim(bin,ndime,mnode,nelem,npoin,CRITERION,permn)

    class(maths_bin),                       intent(inout) :: bin
    integer(ip),                            intent(out)   :: ndime
    integer(ip),                            intent(out)   :: mnode
    integer(ip),                            intent(out)   :: nelem
    integer(ip),                            intent(out)   :: npoin
    integer(ip),         optional,          intent(in)    :: CRITERION     !< If empty bins should be considered
    integer(ip),         optional, pointer, intent(inout) :: permn(:)
    integer(ip)                                           :: idime,inode,ipoin
    integer(ip)                                           :: ii,jj,kk,kpoin
    integer(ip)                                           :: lnods(8)
    integer(ip),                   pointer                :: permn_loc(:)
    integer(ip)                                           :: my_criterion

    nullify(permn_loc)
    
    my_criterion = optional_argument(SEARCH_MESH_LEAVES , CRITERION)

    
    ndime = bin % dim
    mnode = (ndime-1)*4

    if( my_criterion == SEARCH_MESH_LEAVES ) then
       nelem = 1
       npoin = 1
       do idime = 1,ndime
          nelem = nelem *   bin % boxip(idime)
          npoin = npoin * ( bin % boxip(idime) + 1 )
       end do
    else
       nelem = bin % nfilled
       npoin = 1
       do idime = 1,ndime
          npoin = npoin * ( bin % boxip(idime) + 1 )
       end do
       if( present(permn) ) then
          allocate(permn(npoin))
          permn_loc => permn
       else
          allocate(permn_loc(npoin))
       end if
       do ipoin = 1,npoin
          permn_loc(ipoin) = 0
       end do
       kpoin = 0
       do kk = 1,bin % boxip(3)
          do jj = 1,bin % boxip(2)
             do ii = 1,bin % boxip(1)
                if( bin % filled((/ii,jj,kk/)) ) then
                   call bin % mesh_node((/ii,jj,kk/),lnods)
                   do inode = 1,mnode
                      ipoin = lnods(inode)
                      if( permn_loc(ipoin) == 0 ) then
                         kpoin = kpoin + 1
                         permn_loc(ipoin) = kpoin
                      end if
                   end do
                end if
             end do
          end do
       end do
       npoin = kpoin
       if( .not. present(permn) ) deallocate(permn_loc)
    end if
    
  end subroutine mesh_dim
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   Create a mesh
  !> @details Create a mesh from an octbin
  !> 
  !-----------------------------------------------------------------------
  
  pure subroutine mesh_coord(bin,xx,coord)

    class(maths_bin),  intent(inout) :: bin
    integer(ip),       intent(in)    :: xx(:)
    real(rp),          intent(out)   :: coord(:,:)
    integer(ip)                      :: ndime
    real(rp)                         :: comin(3),comax(3)
    real(rp)                         :: delta(3)
    
    ndime          = bin % dim
    delta(1:ndime) = ( bin % comax(1:ndime) - bin % comin(1:ndime) ) / bin % boxrp(1:ndime)
    comin(1:ndime) = bin % comin(1:ndime) + real(xx(1:ndime) - 1,rp) * delta(1:ndime)
    comax(1:ndime) = bin % comin(1:ndime) + real(xx(1:ndime)    ,rp) * delta(1:ndime)
    
    if( bin % dim == 2 ) then

       coord(1:ndime,1) = (/ comin(1) , comin(2) /)
       coord(1:ndime,2) = (/ comax(1) , comin(2) /)
       coord(1:ndime,3) = (/ comax(1) , comax(2) /)
       coord(1:ndime,4) = (/ comin(1) , comax(2) /)
              
    else if( bin % dim == 3 ) then

       coord(1:ndime,1) = (/ comin(1) , comin(2) , comin(3) /)
       coord(1:ndime,2) = (/ comax(1) , comin(2) , comin(3) /)
       coord(1:ndime,3) = (/ comax(1) , comax(2) , comin(3) /)
       coord(1:ndime,4) = (/ comin(1) , comax(2) , comin(3) /)
       coord(1:ndime,5) = (/ comin(1) , comin(2) , comax(3) /)
       coord(1:ndime,6) = (/ comax(1) , comin(2) , comax(3) /)
       coord(1:ndime,7) = (/ comax(1) , comax(2) , comax(3) /)
       coord(1:ndime,8) = (/ comin(1) , comax(2) , comax(3) /)
 
    end if
    
  end subroutine mesh_coord
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   Create a mesh
  !> @details Create a mesh from an octbin
  !> 
  !-----------------------------------------------------------------------
  
  pure subroutine mesh_node(bin,xx,lnods)

    class(maths_bin),  intent(inout) :: bin
    integer(ip),       intent(in)    :: xx(:)
    integer(ip),       intent(out)   :: lnods(:)
    integer(ip)                      :: ipoin,kpoin
    
    if( bin % dim == 2 ) then

       ipoin    = xx(1) + (bin % boxip(1)+1)*(xx(2)-1) 
       lnods(1) = ipoin
       lnods(2) = ipoin+1
       lnods(3) = ipoin+bin % boxip(1)+2
       lnods(4) = ipoin+bin % boxip(1)+1
       
    else if( bin % dim == 3 ) then

       kpoin    = (bin % boxip(1)+1)*(bin % boxip(2)+1)
       ipoin    = xx(1) + (bin % boxip(1)+1)*(xx(2)-1) + (bin % boxip(1)+1)*(bin % boxip(2)+1)*(xx(3)-1)
       lnods(1) = ipoin
       lnods(2) = ipoin + 1
       lnods(3) = ipoin + bin % boxip(1)+2
       lnods(4) = ipoin + bin % boxip(1)+1
       lnods(5) = ipoin                    + kpoin
       lnods(6) = ipoin + 1                + kpoin
       lnods(7) = ipoin + bin % boxip(1)+2 + kpoin
       lnods(8) = ipoin + bin % boxip(1)+1 + kpoin              
       
    end if
    
  end subroutine mesh_node
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   Create a mesh
  !> @details Create a mesh from an octbin
  !> 
  !-----------------------------------------------------------------------
  
  subroutine mesh(self,ndime,mnode,nelem,npoin,lnods,ltype,coord,MEMORY_COUNTER,&
       OFFSET_IELEM,OFFSET_IPOIN,CRITERION,CENTROID,ONLY_DEALLOCATE)

    class(maths_bin),                       intent(inout) :: self
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
    integer(ip)                                           :: ielem,ipoin
    integer(ip)                                           :: ii,jj,kk,kpoin
    integer(ip)                                           :: offset_ielem_loc
    integer(ip)                                           :: offset_ipoin_loc
    integer(8)                                            :: memor_loc(2)
    real(rp)                                              :: comin(3)
    real(rp)                                              :: comax(3)
    real(rp)                                              :: delta(3)
    real(rp)                                              :: yy,zz
    integer(ip),                   pointer                :: permn(:)
    integer(ip)                                           :: lnods_loc(8)
    logical(lg)                                           :: if_deallocate
    logical(lg)                                           :: if_empty
    integer(ip)                                           :: my_criterion

    nullify(permn)

    offset_ielem_loc = optional_argument(0_ip               , OFFSET_IELEM)
    offset_ipoin_loc = optional_argument(0_ip               , OFFSET_IPOIN)
    my_criterion     = optional_argument(SEARCH_MESH_LEAVES , CRITERION)
    if_deallocate    = optional_argument(.false.            , ONLY_DEALLOCATE)
    if_empty         = my_criterion == SEARCH_MESH_LEAVES
    
    call self % mem_ini(memor_loc,MEMORY_COUNTER)
    call self % mesh_dim(ndime,mnode,nelem,npoin,CRITERION,permn)

    if( if_deallocate ) then
       !
       ! Deallocate
       !
       call memory_deallo(memor_loc,'LNODS','mesh',lnods)
       call memory_deallo(memor_loc,'LTYPE','mesh',ltype)
       call memory_deallo(memor_loc,'COORD','mesh',coord)

    else

       comin(1:ndime) = self % comin(1:ndime)
       comax(1:ndime) = self % comax(1:ndime)
       delta(1:ndime) = comax(1:ndime) - comin(1:ndime)
       !
       ! Allocate mesh
       !
       if( .not. associated(lnods) ) call memory_alloca(memor_loc,'LNODS','mesh',lnods,mnode,nelem)
       if( .not. associated(ltype) ) call memory_alloca(memor_loc,'LTYPE','mesh',ltype,nelem)
       if( .not. associated(coord) ) call memory_alloca(memor_loc,'COORD','mesh',coord,ndime,npoin)
       !
       ! Mesh arrays
       !
       ipoin = offset_ipoin_loc
       ielem = offset_ielem_loc
       if( if_empty ) then
          if( ndime == 2 ) then
             do jj = 1,self % boxip(2)+1
                yy = comin(2) + real(jj-1,rp)/self % boxrp(2) * delta(2)
                do ii = 1,self % boxip(1)+1
                   ipoin = ipoin + 1
                   coord(1,ipoin) = comin(1) + real(ii-1,rp)/self % boxrp(1) * delta(1)
                   coord(2,ipoin) = yy
                end do
             end do
             ipoin = 1
             do jj = 1,self % boxip(2)
                do ii = 1,self % boxip(1)
                   ielem          = ielem + 1
                   lnods(1,ielem) = ipoin
                   lnods(2,ielem) = ipoin+1
                   lnods(3,ielem) = ipoin+self % boxip(1)+2
                   lnods(4,ielem) = ipoin+self % boxip(1)+1
                   ltype(ielem)   = QUA04
                   ipoin          = ipoin + 1
                end do
                ipoin = ipoin + 1
             end do
          else
             do kk = 1,self % boxip(3)+1
                zz = comin(3) + real(kk-1,rp)/self % boxrp(3) * delta(3)
                do jj = 1,self % boxip(2)+1
                   yy = comin(2) + real(jj-1,rp)/self % boxrp(2) * delta(2)
                   do ii = 1,self % boxip(1)+1
                      ipoin = ipoin + 1
                      coord(1,ipoin) = comin(1) + real(ii-1,rp)/self % boxrp(1) * delta(1)
                      coord(2,ipoin) = yy
                      coord(3,ipoin) = zz
                   end do
                end do
             end do
             ipoin = 1
             kpoin = (self % boxip(1)+1)*(self % boxip(2)+1)
             do kk = 1,self % boxip(3)
                do jj = 1,self % boxip(2)
                   do ii = 1,self % boxip(1)
                      ielem          = ielem + 1                
                      lnods(1,ielem) = ipoin
                      lnods(2,ielem) = ipoin+1
                      lnods(3,ielem) = ipoin+self % boxip(1)+2
                      lnods(4,ielem) = ipoin+self % boxip(1)+1
                      lnods(5,ielem) = kpoin + ipoin
                      lnods(6,ielem) = kpoin + ipoin+1
                      lnods(7,ielem) = kpoin + ipoin+self % boxip(1)+2
                      lnods(8,ielem) = kpoin + ipoin+self % boxip(1)+1                
                      ltype(ielem)   = HEX08
                      ipoin          = ipoin + 1
                   end do
                   ipoin = ipoin + 1
                end do
                ipoin = kpoin * kk + 1
             end do
          end if

       else

          if( ndime == 2 ) then
             do jj = 1,self % boxip(2)+1
                yy = comin(2) + real(jj-1,rp)/self % boxrp(2) * delta(2)
                do ii = 1,self % boxip(1)+1
                   ipoin          = ipoin + 1
                   kpoin          = permn(ipoin)
                   if( kpoin > 0 ) then
                      coord(1,kpoin) = comin(1) + real(ii-1,rp)/self % boxrp(1) * delta(1)
                      coord(2,kpoin) = yy
                   end if
                end do
             end do
             do jj = 1,self % boxip(2)
                do ii = 1,self % boxip(1)
                   if( self % filled((/ii,jj,1_ip/)) ) then !associated(self % list(ii,jj,1_ip) % l) ) then
                      ielem                = ielem + 1
                      call self % mesh_node((/ii,jj,1_ip/),lnods_loc)
                      lnods(1:mnode,ielem) = permn(lnods_loc(1:mnode))
                      ltype(ielem)         = QUA04
                   end if
                end do
             end do
          else
             do kk = 1,self % boxip(3)+1
                zz = comin(3) + real(kk-1,rp)/self % boxrp(3) * delta(3)
                do jj = 1,self % boxip(2)+1
                   yy = comin(2) + real(jj-1,rp)/self % boxrp(2) * delta(2)
                   do ii = 1,self % boxip(1)+1
                      ipoin          = ipoin + 1
                      kpoin          = permn(ipoin)
                      if( kpoin > 0 ) then
                         coord(1,kpoin) = comin(1) + real(ii-1,rp)/self % boxrp(1) * delta(1)
                         coord(2,kpoin) = yy
                         coord(3,kpoin) = zz
                      end if
                   end do
                end do
             end do
             do kk = 1,self % boxip(3)
                do jj = 1,self % boxip(2)
                   do ii = 1,self % boxip(1)
                      if( self % filled((/ii,jj,kk/)) ) then  !if( associated(self % list(ii,jj,kk) % l) ) then
                         ielem                = ielem + 1
                         call self % mesh_node((/ii,jj,kk/),lnods_loc)
                         lnods(1:mnode,ielem) = permn(lnods_loc(1:mnode))
                         ltype(ielem)         = HEX08
                      end if
                   end do
                end do
             end do
          end if
       end if

       if( associated(permn) ) deallocate(permn)

       if( present(CENTROID) ) then
          call self % centroid(CENTROID,MEMORY_COUNTER=MEMORY_COUNTER,CRITERION=my_criterion)
       end if

    end if
    call self % mem_end(memor_loc,MEMORY_COUNTER)

  end subroutine mesh
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Coordinates of a bin
  !> @details Coordinates of a bin
  !> 
  !-----------------------------------------------------------------------

  pure subroutine coord_2(bin,xx,comin,comax) 

    class(maths_bin), intent(in)  :: bin
    integer(ip),      intent(in)  :: xx(:)    !< Box in x direction
    real(rp),         intent(out) :: comin(:) !< Coordinate of the test point
    real(rp),         intent(out) :: comax(:) !< Coordinate of the test point
    integer(ip)                   :: ndime
    real(rp)                      :: delta(3)
    
    comin          =  big
    comax          = -big
    ndime          =  bin % dim
    delta(1:ndime) = ( bin % comax(1:ndime) - bin % comin(1:ndime) ) / bin % boxrp(1:ndime)
    
    comin(1:ndime) = bin % comin(1:ndime) + real(xx(1:ndime) - 1,rp) * delta(1:ndime)
    comax(1:ndime) = bin % comin(1:ndime) + real(xx(1:ndime)    ,rp) * delta(1:ndime)
    
  end subroutine coord_2
  
  pure subroutine coord_1(bin,xx,bobox) 

    class(maths_bin), intent(in)  :: bin
    integer(ip),      intent(in)  :: xx(:)      !< Box in x direction
    real(rp),         intent(out) :: bobox(:,:) !< Coordinate of the test point
    integer(ip)                   :: ndime
    real(rp)                      :: delta(3)
    
    bobox(1,:)     =  big
    bobox(2,:)     = -big
    ndime          =  bin % dim
    delta(1:ndime) = ( bin % comax(1:ndime) - bin % comin(1:ndime) ) / bin % boxrp(1:ndime)
    
    bobox(1,1:ndime) = bin % comin(1:ndime) + real(xx(1:ndime) - 1,rp) * delta(1:ndime)
    bobox(2,1:ndime) = bin % comin(1:ndime) + real(xx(1:ndime)    ,rp) * delta(1:ndime)
    
  end subroutine coord_1
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Nearest bins
  !> @details Find the nearest bins
  !> 
  !-----------------------------------------------------------------------

  pure subroutine near(bin,coord,xx,nn)

    class(maths_bin),         intent(in)    :: bin
    real(rp),                 intent(in)    :: coord(:)    !< Coordinate of the test point
    integer(ip),              intent(inout) :: xx(:,:)     !< Box in direction
    integer(ip),              intent(out)   :: nn
    integer(ip)                             :: ii,jj,kk,nd
    integer(ip)                             :: i1,j1,k1
    real(rp)                                :: dimax
    real(rp)                                :: dimin
    real(rp)                                :: dismm(2)
    real(rp)                                :: bobox(2,3)
    real(rp),         pointer               :: mindi(:,:,:)

    nd     = bin % dim
    dimin  = big
    nn     = 0

    if( bin % nfilled > 0 ) then
       
       allocate( mindi(bin % boxip(1),bin % boxip(2),bin % boxip(3)) )

       do kk = 1,bin % boxip(3)
          do jj = 1,bin % boxip(2)
             do ii = 1,bin % boxip(1)
                if( bin % num_points( (/ii,jj,kk/) ) > 0 ) then
                   call bin % coord((/ii,jj,kk/),bobox)
                   dismm(1)        = maths_min_box_vertices(coord(1:nd),bobox(1:2,1:nd))
                   mindi(ii,jj,kk) = dismm(1)
                   if( dismm(1) <= dimin ) then
                      dimin = dismm(1)
                      i1    = ii
                      j1    = jj
                      k1    = kk
                   end if
                else
                   mindi(ii,jj,kk) = big
                end if
             end do
          end do
       end do
       
       call bin % coord((/i1,j1,k1/),bobox)
       dimax = maths_max_box_vertices(coord(1:nd),bobox(1:2,1:nd))

       do kk = 1,bin % boxip(3)
          do jj = 1,bin % boxip(2)
             do ii = 1,bin % boxip(1)
                if( mindi(ii,jj,kk) <= dimax + epsil ) then 
                   nn = nn + 1
                   xx(1:3,nn) = (/ ii,jj,kk /)
                end if
             end do
          end do
       end do

       deallocate(mindi)
    end if

  end subroutine near
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Location in a bin
  !> @details Find the location of a point in a bin
  !> 
  !-----------------------------------------------------------------------

  pure subroutine find(bin,coord,xx) 

    class(maths_bin), intent(in)    :: bin
    real(rp),         intent(in)    :: coord(:) !< Coordinate of the test point
    integer(ip),      intent(inout) :: xx(3)    !< Box in x direction
       
    select case ( bin % dim ) 

    case ( 1_ip )

       if(     coord(1) >= bin % comin(1)-epsil .and. coord(1) <= bin % comax(1)+epsil ) then          
          xx(1) = int( ( (coord(1)-bin % comin(1)-epsil) / (bin % comax(1)-bin % comin(1)) ) * bin % boxrp(1), ip ) + 1
          xx(1) = min(max(1_ip,xx(1)),bin % boxip(1))
          xx(2) = 1
          xx(3) = 1
       else
          xx    = -1
       end if

    case ( 2_ip )

       if(     coord(1) >= bin % comin(1)-epsil .and. coord(1) <= bin % comax(1)+epsil .and. &
            &  coord(2) >= bin % comin(2)-epsil .and. coord(2) <= bin % comax(2)+epsil ) then         
          xx(1:2) = int( ( (coord(1:2)-bin % comin(1:2)-epsil) / (bin % comax(1:2)-bin % comin(1:2)) ) * bin % boxrp(1:2), ip ) + 1
          xx(1:2) = min(max(1_ip,xx(1:2)),bin % boxip(1:2))
          xx(3)   = 1
       else
          xx      = -1
       end if
       
    case ( 3_ip )

      if(      coord(1) >= bin % comin(1)-epsil .and. coord(1) <= bin % comax(1)+epsil .and. &
            &  coord(2) >= bin % comin(2)-epsil .and. coord(2) <= bin % comax(2)+epsil .and. &
            &  coord(3) >= bin % comin(3)-epsil .and. coord(3) <= bin % comax(3)+epsil ) then         
         xx(1:3) = int( ( (coord(1:3)-bin % comin(1:3)-epsil) / (bin % comax(1:3)-bin % comin(1:3)) ) * bin % boxrp(1:3), ip ) + 1
         xx(1:3) = min(max(1_ip,xx(1:3)),bin % boxip(1:3))
       else
          xx     = -1
       end if

    end select

  end subroutine find
  
  pure subroutine find_bb(bin,coord,xx,number_boxes) 

    class(maths_bin),           intent(in)  :: bin
    real(rp),                   intent(in)  :: coord(:,:)   !< Coordinate of the test point
    integer(ip),                intent(out) :: xx(:,:)      !< Boxes in each direction
    integer(ip),      optional, intent(out) :: number_boxes !< Number of boxes
    integer(ip)                             :: ii,jj,kk,ll
    integer(ip)                             :: imin,imax
    integer(ip)                             :: jmin,jmax
    integer(ip)                             :: kmin,kmax
    real(rp)                                :: tole(3)

    tole = bin % delta * bin % toler_rel 
    ll   = 0
    
    select case ( bin % dim ) 

    case ( 1_ip )

       if(  coord(1,1) <= bin % comax(1) + tole(1) .and. coord(2,1) >= bin % comin(1) - tole(1) ) then

          imin = int( ( coord(1,1) - bin % comin(1) - tole(1) ) * bin % delti(1) , ip )  + 1
          imax = int( ( coord(2,1) - bin % comin(1) + tole(1) ) * bin % delti(1) , ip )  + 1      

          imin = max(imin,1_ip)
          imax = min(imax,bin % boxip(1))

          ll = 0
          do ii = imin,imax
             ll       = ll + 1
             xx(1,ll) = ii
             xx(2,ll) = 1
             xx(3,ll) = 1
          end do

       end if

    case ( 2_ip )
          
       if(  coord(1,1) <= bin % comax(1) + tole(1) .and. coord(2,1) >= bin % comin(1) - tole(1) .and. &
            coord(1,2) <= bin % comax(2) + tole(2) .and. coord(2,2) >= bin % comin(2) - tole(2)) then

          imin = int( ( coord(1,1) - bin % comin(1) - tole(1) ) * bin % delti(1) , ip )  + 1
          imax = int( ( coord(2,1) - bin % comin(1) + tole(1) ) * bin % delti(1) , ip )  + 1      

          jmin = int( ( coord(1,2) - bin % comin(2) - tole(2) ) * bin % delti(2) , ip )  + 1
          jmax = int( ( coord(2,2) - bin % comin(2) + tole(2) ) * bin % delti(2) , ip )  + 1      

          imin = max(imin,1_ip)
          imax = min(imax,bin % boxip(1))
          jmin = max(jmin,1_ip)
          jmax = min(jmax,bin % boxip(2))

          ll = 0
          do ii = imin,imax
             do jj = jmin,jmax
                ll       = ll + 1                
                xx(1,ll) = ii
                xx(2,ll) = jj
                xx(3,ll) = 1
             end do
          end do

       end if

    case ( 3_ip )

       if(  coord(1,1) <= bin % comax(1) + tole(1) .and. coord(2,1) >= bin % comin(1) - tole(1) .and. &
            coord(1,2) <= bin % comax(2) + tole(2) .and. coord(2,2) >= bin % comin(2) - tole(2) .and. &
            coord(1,3) <= bin % comax(3) + tole(3) .and. coord(2,3) >= bin % comin(3) - tole(3)) then

          imin = int( ( coord(1,1) - bin % comin(1) - tole(1) ) * bin % delti(1) , ip )  + 1
          imax = int( ( coord(2,1) - bin % comin(1) + tole(1) ) * bin % delti(1) , ip )  + 1      

          jmin = int( ( coord(1,2) - bin % comin(2) - tole(2) ) * bin % delti(2) , ip )  + 1
          jmax = int( ( coord(2,2) - bin % comin(2) + tole(2) ) * bin % delti(2) , ip )  + 1      

          kmin = int( ( coord(1,3) - bin % comin(3) - tole(3) ) * bin % delti(3) , ip )  + 1
          kmax = int( ( coord(2,3) - bin % comin(3) + tole(3) ) * bin % delti(3) , ip )  + 1      

          imin = max(imin,1_ip)
          imax = min(imax,bin % boxip(1))
          jmin = max(jmin,1_ip)
          jmax = min(jmax,bin % boxip(2))
          kmin = max(kmin,1_ip)
          kmax = min(kmax,bin % boxip(3))

          ll = 0
          do kk = kmin,kmax
             do jj = jmin,jmax
                do ii = imin,imax
                   ll       = ll + 1
                   xx(1,ll) = ii
                   xx(2,ll) = jj
                   xx(3,ll) = kk
                end do
             end do
          end do

       end if

    end select

    if( present(number_boxes) ) number_boxes = ll

  end subroutine find_bb
  
   !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Location in a bin
  !> @details Find the location of a point in a bin
  !> 
  !-----------------------------------------------------------------------

  pure logical(lg) function inbin(bin,coord) 

    class(maths_bin), intent(in)  :: bin
    real(rp),         intent(in)  :: coord(:) !< Coordinate of the test point

    inbin = .true.

    select case (  bin % dim )

    case ( 1_ip )
       
       if( coord(1) > bin % comax(1) .or.  coord(1) < bin % comin(1) ) inbin = .false.
       
    case ( 2_ip )
       
       if(     coord(1) > bin % comax(1) ) then
          inbin = .false. ; return
       else if( coord(1) < bin % comin(1) ) then
          inbin = .false. ; return
       else if( coord(2) > bin % comax(2) ) then
          inbin = .false. ; return
       else if( coord(2) < bin % comin(2) ) then
          inbin = .false. ; return
       end if
       
     case ( 3_ip )
      
       if(     coord(1) > bin % comax(1) ) then
          inbin = .false. ; return
       else if( coord(1) < bin % comin(1) ) then
          inbin = .false. ; return
       else if( coord(2) > bin % comax(2) ) then
          inbin = .false. ; return
       else if( coord(2) < bin % comin(2) ) then
          inbin = .false. ; return
       else if( coord(3) > bin % comax(3) ) then
          inbin = .false. ; return
       else if( coord(3) < bin % comin(3) ) then
          inbin = .false. ; return
       end if
       
    end select

  end function inbin

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-10
  !> @brief   Results
  !> @details Get some results
  !> 
  !-----------------------------------------------------------------------

  subroutine results(self,xx,names,OFFSET,MEMORY_COUNTER,CRITERION,ONLY_DEALLOCATE)

    class(maths_bin),                    intent(inout) :: self
    real(rp),                   pointer, intent(inout) :: xx(:,:)
    character(len=5),           pointer, intent(inout) :: names(:)
    integer(ip),      optional,          intent(in)    :: OFFSET
    integer(8),       optional,          intent(inout) :: MEMORY_COUNTER(2)
    integer(ip),      optional,          intent(in)    :: CRITERION     !< If empty bins should be considered
    logical(lg),      optional,          intent(in)    :: ONLY_DEALLOCATE
    integer(ip)                                        :: ii,jj,kk
    integer(ip)                                        :: ielem,nelem
    integer(ip)                                        :: ndime,mnode,npoin
    integer(ip)                                        :: ielem_offset
    integer(8)                                         :: memor_loc(2)
    logical(lg)                                        :: if_deallocate
    integer(ip)                                        :: my_criterion

    call self % mem_ini(memor_loc,MEMORY_COUNTER)

    if_deallocate = optional_argument(.false.,ONLY_DEALLOCATE)
    ielem_offset  = optional_argument(0_ip,OFFSET)
    my_criterion  = optional_argument(SEARCH_MESH_LEAVES , CRITERION)
     
    if( if_deallocate ) then
       !
       ! Deallocate
       !
       call memory_deallo(memor_loc,'XX'   ,'results',xx)
       call memory_deallo(memor_loc,'NAMES','results',names)
       
    else
       !
       ! Deallocate results and fill in results
       !
       call self % mesh_dim(ndime,mnode,nelem,npoin,CRITERION)
       if( .not. associated(xx)    ) call memory_alloca(memor_loc,'XX'   ,'results',xx,nelem,1_ip)
       if( .not. associated(names) ) call memory_alloca(memor_loc,'NAMES','results',names,1_ip)

       if( my_criterion == SEARCH_MESH_LEAVES ) then
          ielem = ielem_offset
          do kk = 1,self % boxip(3)
             do jj = 1,self % boxip(2)
                do ii = 1,self % boxip(1)
                   ielem = ielem + 1                
                   names(1)    = 'NELEM'
                   xx(ielem,1) = real(num_points(self,(/ii,jj,kk/)),rp)
                end do
             end do
          end do
       else
          ielem = ielem_offset
          do kk = 1,self % boxip(3)
             do jj = 1,self % boxip(2)
                do ii = 1,self % boxip(1)
                   if( memory_size(self % list(ii,jj,kk) %l) > 0 ) then
                      ielem = ielem + 1                
                      names(1)    = 'NELEM'
                      xx(ielem,1) = real(num_points(self,(/ii,jj,kk/)),rp)
                   end if
                end do
             end do
          end do          
       end if
       
    end if

    call self % mem_end(memor_loc,MEMORY_COUNTER)
    
  end subroutine results

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-10
  !> @brief   3D to 1D
  !> @details 3D coordinates to 1D coordinates
  !> 
  !-----------------------------------------------------------------------

  pure function box_num(self,xx)
    
    class(maths_bin), intent(in) :: self
    integer(ip),      intent(in) :: xx(3)
    integer(ip)                  :: box_num
    
    box_num = xx(1) + self % boxip(1) * ( xx(2)-1 + (xx(3)-1) * self % boxip(2) )
    
  end function box_num

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-10
  !> @brief   If filled
  !> @details If a box is filled
  !> 
  !-----------------------------------------------------------------------

  pure function filled(self,xx)
    
    class(maths_bin), intent(in) :: self
    integer(ip),      intent(in) :: xx(3)
    logical(lg)                  :: filled
    integer(ip)                  :: ib

    filled = .false.
    
    if( self % fill_type == LINKED_LIST_BIN ) then
       ib = self % box_num(xx)
       if( self % ia(ib+1)-self % ia(ib) > 0 ) filled = .true.
    else
       if( associated(self % list(xx(1),xx(2),xx(3)) % l) ) filled = .true.
    end if
    
  end function filled
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-10
  !> @brief   Number of 
  !> @details Number
  !> 
  !-----------------------------------------------------------------------

  pure function num_points(self,xx)
    
    class(maths_bin), intent(in) :: self
    integer(ip),      intent(in) :: xx(3)
    integer(ip)                  :: num_points
    integer(ip)                  :: ib
    
    if( self % fill_type == LINKED_LIST_BIN ) then
       ib         = self % box_num(xx)
       num_points = self % ia(ib+1)-self % ia(ib) 
    else
       num_points = self % num(xx(1),xx(2),xx(3)) 
    end if
    
  end function num_points

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/01/1016
  !> @brief   Transform a sum to a linked list
  !> @details Transform a sum to a linked list
  !
  !-----------------------------------------------------------------------

  subroutine size_to_linked_list(ia)
    
    integer(ip), pointer, intent(inout) :: ia(:)
    integer(ip)                         :: nn,ii,kk,ll

    nn    = memory_size(ia)-1
    if( nn > 0 ) then
       kk    = ia(1)
       ia(1) = 1 
       do ii = 2,nn+1
          ll     = ia(ii)
          ia(ii) = ia(ii-1) + kk
          kk     = ll
       end do
    end if
    
  end subroutine size_to_linked_list

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    24/01/2022
  !> @brief   Search type name
  !> @details Search type name
  !
  !-----------------------------------------------------------------------

  function type_name(self) result(name)
    class(maths_bin),                        intent(inout) :: self
    character(LEN=:), allocatable                          :: name

    name = 'BIN'
    
  end function type_name

end module def_maths_bin
!> @}
