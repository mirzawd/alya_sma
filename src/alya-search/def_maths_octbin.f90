!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @defgroup Maths_Toolbox
!> Toolbox for octbins
!> @{
!> @name    ToolBox for octbins
!> @file    def_maths.f90
!> @author  Guillaume Houzeaux
!> @brief   Variables
!> @details Variables fot mod_maths.f90
!
!-----------------------------------------------------------------------

module def_maths_octbin

  use def_kintyp_basic,      only : ip,rp,lg,i1p
  use mod_memory_basic,      only : memory_alloca
  use mod_memory_basic,      only : memory_deallo
  use mod_memory_basic,      only : memory_size
  use def_elmtyp,            only : HEX08
  use def_elmtyp,            only : QUA04
  use def_search_method,     only : search_method
  use def_search_method,     only : SEARCH_FILL
  use def_search_method,     only : SEARCH_CANDIDATE
  use def_search_method,     only : SEARCH_MESH_LEAVES
  use mod_optional_argument, only : optional_argument
  use def_maths_bin,         only : maths_bin
  use def_maths_bin,         only : LINKED_LIST_BIN
  use def_maths_tree,        only : maths_octree

  real(rp), parameter :: epsil = epsilon(1.0_rp)

  type test
     type(maths_octree) :: octree
  end type test
  
  type, extends(search_method) :: maths_octbin
     integer(ip)                     :: limit      ! Octree limit  (input data)
     logical(lg)                     :: enable_gap ! Stick or not the octree to the bin
     integer(ip)                     :: num_octree
     type(maths_bin)                 :: bin
     type(maths_octree), allocatable :: octree(:,:,:)
   contains
     procedure,         pass :: init          ! Initialize all
     procedure,         pass :: deallo        ! Deallocate
     procedure,         pass :: candidate     ! Candidate entities
     procedure,         pass :: input         ! Input parameters
     procedure,         pass :: results       ! Get results on the mesh
     procedure,         pass :: mesh          ! Mesh
     procedure,         pass :: type_name     ! Type name

     procedure,         pass :: fill          ! Fill the octbin
     procedure,         pass :: mesh_dim      ! Mesh dimensions
     procedure,         pass :: centroid      ! Get centroid
  end type maths_octbin

  private

  public :: maths_octbin

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Input
  !> @details Set input data
  !> 
  !-----------------------------------------------------------------------

  subroutine input(self,boxes,limit,enable_gap,relative_tolerance,absolute_tolerance,fill_method)
    
    class(maths_octbin),        intent(inout) :: self
    integer(ip),      optional, intent(in)    :: boxes(:)
    integer(ip),      optional, intent(in)    :: limit
    logical(lg),      optional                :: enable_gap
    real(rp),         optional, intent(in)    :: relative_tolerance
    real(rp),         optional, intent(in)    :: absolute_tolerance
    integer(ip),      optional, intent(in)    :: fill_method
    integer(ip)                               :: ii
    
    if( present(boxes) ) then
       do ii = 1,int(size(boxes),ip)
          self % bin % boxip(ii) = boxes(ii)
       end do
    end if
    if( present(limit) )      self % limit      = limit
    if( present(enable_gap) ) self % enable_gap = enable_gap
    
    call self % input_all(relative_tolerance,absolute_tolerance,fill_method)
   
  end subroutine input

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Octbin
  !> @details Fill in a octbin from a coordinate field
  !> 
  !-----------------------------------------------------------------------

  subroutine init(self)

    class(maths_octbin), intent(inout) :: self

    call self % init_all()
    call self % bin % init()
    self % num_octree = 0
    self % name = 'OCT-BIN'
    
  end subroutine init

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Deallocate
  !> @details Deallocate and octebin
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo(self,list_entities,MEMORY_COUNTER)

    class(maths_octbin),                     intent(inout) :: self
    type(i1p),            optional, pointer, intent(inout) :: list_entities(:)    
    integer(8),           optional,          intent(inout) :: MEMORY_COUNTER(2)
    integer(ip)                                            :: ii,jj,kk
    integer(8)                                             :: memor_loc(2)

    call self % mem_ini(memor_loc,MEMORY_COUNTER)    

    do kk = 1,self % bin % boxip(3)
       do jj = 1,self % bin % boxip(2)
          do ii = 1,self % bin % boxip(1)
             call self % octree(ii,jj,kk) % deallo(MEMORY_COUNTER=memor_loc)
          end do
       end do
    end do
    if( allocated(self % octree) ) deallocate(self % octree)
    call self % bin % deallo(MEMORY_COUNTER=memor_loc)
    if( allocated(self % name) ) deallocate(self % name)

    if( present(list_entities) ) call memory_deallo(memor_loc,'LIST_ENTITIES','def_maths_octbin',list_entities)
    
    call self % mem_end(memor_loc,MEMORY_COUNTER)    
    
  end subroutine deallo

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Search
  !> @details Search in parallel
  !> 
  !-----------------------------------------------------------------------

  subroutine candidate(self,xx,list_entities,METHOD,MASK,MEMORY_COUNTER)

    class(maths_octbin),                 intent(inout) :: self
    real(rp),                            intent(in)    :: xx(:,:)           !< List coordinates
    type(i1p),                  pointer, intent(inout) :: list_entities(:)    
    integer(8),       optional,          intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(ip),      optional,          intent(in)    :: METHOD                !< Method for candidates
    logical(lg),      optional, pointer, intent(in)    :: MASK(:)           !< Mask to consider or not points
    integer(ip)                                        :: ii,nn,ndime
    integer(ip)                                        :: nenti
    integer(8)                                         :: memor_loc(2)
    integer(ip)                                        :: xloc(3)
    real(rp),                   allocatable            :: xx1(:,:)
    type(i1p),                  pointer                :: my_list_entities(:)    
    logical(lg)                                        :: mask_loc

    nullify(my_list_entities)
    
    call self % tim_ini()
    call self % mem_ini(memor_loc,MEMORY_COUNTER)    

    ndime    = int(size(xx,1_ip),ip)
    nn       = int(size(xx,2_ip),ip)
    mask_loc = .true.
    
    if( .not. associated(list_entities) ) & 
         call memory_alloca(memor_loc,'LIST_ENTITIES','maths_bin',list_entities,nn)
    
    call memory_alloca(memor_loc,'MY_LIST_ENTITIES','maths_bin',my_list_entities,1_ip)
    allocate(xx1(ndime,1))
    do ii = 1,nn
       if( present(MASK) ) then
          if( associated(MASK) ) mask_loc = MASK(ii)
       end if
       if( mask_loc ) then
          xx1(1:ndime,1) = xx(1:ndime,ii)
          call self % bin % find(xx1(:,1),xloc)
          if( minval(xloc) > 0 ) then
             call self % octree(xloc(1),xloc(2),xloc(3)) % candidate(xx1,my_list_entities,MEMORY_COUNTER=memor_loc)
             nenti = memory_size(my_list_entities(1) % l )
             if( nenti > 0 ) then
                call memory_alloca(memor_loc,'LIST_ENTITIES','maths_bin',list_entities(ii) % l,nenti)
                list_entities(ii) % l(1:nenti) = my_list_entities(1) % l(1:nenti) 
                call memory_deallo(memor_loc,'MY_LIST_ENTITIES','maths_bin',my_list_entities(1) % l)
             end if
          end if
       end if
    end do
    deallocate(xx1)
    call memory_deallo(memor_loc,'MY_LIST_ENTITIES','maths_bin',my_list_entities)
    
    call self % mem_end(memor_loc,MEMORY_COUNTER)    
    call self % tim_end(SEARCH_CANDIDATE)

  end subroutine candidate

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Octbin
  !> @details Fill in a octbin from a coordinate field
  !> 
  !-----------------------------------------------------------------------

  subroutine fill(self,coord,bobox,MEMORY_COUNTER,COORD_MIN,COORD_MAX,PERMU,MASK)

    class(maths_octbin),             intent(inout) :: self
    real(rp),    optional,  pointer, intent(in)    :: coord(:,:)
    real(rp),    optional,  pointer, intent(in)    :: bobox(:,:,:)      !< Bounding boxes
    integer(8),  optional,           intent(inout) :: MEMORY_COUNTER(2)
    real(rp),    optional,           intent(in)    :: COORD_MIN(:)
    real(rp),    optional,           intent(in)    :: COORD_MAX(:)
    integer(ip), optional, pointer,  intent(in)    :: PERMU(:)
    logical(lg), optional, pointer,  intent(in)    :: MASK(:)
    integer(ip)                                    :: ii,jj,kk
    integer(ip)                                    :: npoin,ll,ndime
    integer(ip)                                    :: ipoin,nn,ib,iz
    real(rp)                                       :: comin(3)
    real(rp)                                       :: comax(3)
    integer(ip),            pointer                :: permr(:)
    integer(ip),            pointer                :: invpr(:)
    real(rp),               pointer                :: xx(:,:)
    real(rp),               pointer                :: bb(:,:,:)
    integer(8)                                     :: memor_loc(2)
    logical(lg)                                    :: if_gap
    real(rp),               pointer                :: bobox_tmp(:,:,:) ! Just for PGI to work
    real(rp),               pointer                :: xx_tmp(:,:)      ! Just for PGI to work

    nullify(xx_tmp,bobox_tmp)
      
    call self % tim_ini()
    call self % mem_ini(memor_loc,MEMORY_COUNTER)

    if_gap    = self % enable_gap

    nullify(permr)
    nullify(invpr)
    nullify(xx)
    nullify(bb)
    !
    ! Dimension
    !
    if( present(coord) ) then
       npoin      = memory_size(coord,2_ip)
       self % dim = memory_size(coord,1_ip)
       if(associated(coord) ) then
          if( lbound(coord,2) /= 1 ) call runend('DEF_MATHS_OCTBIN: NOT CODED 1')
       end if
    else if( present(bobox) ) then
       npoin      = memory_size(bobox,3_ip)
       self % dim = memory_size(bobox,2_ip)
       if(associated(bobox) ) then
          if( lbound(bobox,3) /= 1 ) call runend('DEF_MATHS_OCTBIN: NOT CODED 2')
       end if
    end if
    self % nelem = npoin 
    
    call self % bin % input (BOXES=self % bin % boxip,RELATIVE_TOLERANCE=self % toler_rel,ABSOLUTE_TOLERANCE=self % toler_abs)
    call self % bin % fill  (coord,bobox,MEMORY_COUNTER=memor_loc,COORD_MIN=COORD_MIN,COORD_MAX=COORD_MAX)

    ndime = self % dim

    allocate(self % octree(self % bin % boxip(1),self % bin % boxip(2),self % bin % boxip(3)))

    do kk = 1,self % bin % boxip(3)
       do jj = 1,self % bin % boxip(2)
          do ii = 1,self % bin % boxip(1)
             
             call self % octree(ii,jj,kk) % init()
             
             !if( associated(self % bin % list(ii,jj,kk) % l) ) then
             if( self % bin % filled( (/ii,jj,kk /) ) ) then
                !
                ! Construct the octree of this bin
                !
                !nn = size( self % bin % list(ii,jj,kk) % l)
                nn = self % bin % num_points((/ii,jj,kk/))
                call self % octree(ii,jj,kk) % input(LIMIT=self % limit,RELATIVE_TOLERANCE=self % toler_rel,ABSOLUTE_TOLERANCE=self % toler_abs,DIM=self % dim)

                if( present(coord) ) then
                   !
                   ! Coordinates
                   !
                   allocate(xx(ndime,nn),permr(nn))
                   if( self % bin % fill_type == LINKED_LIST_BIN ) then
                      if( nn > 0 ) then
                         ib = self % bin % box_num((/ii,jj,kk/))
                         ll = 0
                         do iz = self % bin % ia(ib),self % bin % ia(ib+1)-1
                            ll    = ll + 1
                            ipoin = self % bin % ja(iz)
                            xx(1:ndime,ll) = coord(1:ndime,ipoin)
                            permr(ll)      = ipoin                         
                         end do
                      end if
                   else
                      do ll = 1,nn
                         ipoin          = self % bin % list(ii,jj,kk) % l(ll)
                         xx(1:ndime,ll) = coord(1:ndime,ipoin)
                         permr(ll)      = ipoin
                      end do
                   end if
#ifdef __PGI
                   if( if_gap ) then
                      call self % octree(ii,jj,kk) % fill(xx,bobox_tmp,MEMORY_COUNTER=memor_loc,PERMU=permr)
                   else
                      call self % bin              % coord ((/ii,jj,kk/),comin,comax)
                      call self % octree(ii,jj,kk) % fill  (xx,bobox_tmp,MEMORY_COUNTER=memor_loc,COORD_MIN=comin,COORD_MAX=comax,PERMU=permr)                   
                   end if
#else
                   if( if_gap ) then
                      call self % octree(ii,jj,kk) % fill(xx,MEMORY_COUNTER=memor_loc,PERMU=permr)
                   else
                      call self % bin              % coord ((/ii,jj,kk/),comin,comax)
                      call self % octree(ii,jj,kk) % fill  (xx,MEMORY_COUNTER=memor_loc,COORD_MIN=comin,COORD_MAX=comax,PERMU=permr)                   
                   end if 
#endif
                   self % num_octree = self % num_octree + 1
                   deallocate(xx,permr)
                   
                else if( present(bobox) ) then
                   !
                   ! Bounding boxes
                   !
                   allocate(bb(2,ndime,nn),permr(nn),invpr(npoin))
                   if( self % bin % fill_type == LINKED_LIST_BIN ) then
                      if( nn > 0 ) then
                         ib = self % bin % box_num((/ii,jj,kk/))
                         ll = 0
                         do iz = self % bin % ia(ib),self % bin % ia(ib+1)-1
                            ll                 = ll + 1
                            ipoin              = self % bin % ja(iz)
                            bb(1:2,1:ndime,ll) = bobox(1:2,1:ndime,ipoin)
                            permr(ll)          = ipoin
                            invpr(ipoin)       = ll
                         end do
                      end if
                   else
                      do ll = 1,nn
                         ipoin              = self % bin % list(ii,jj,kk) % l(ll)
                         bb(1:2,1:ndime,ll) = bobox(1:2,1:ndime,ipoin)
                         permr(ll)          = ipoin
                         invpr(ipoin)       = ll
                      end do
                   end if
#ifdef __PGI
                   if( if_gap ) then
                      call self % octree(ii,jj,kk) % fill(xx_tmp,bb,memor_loc,PERMU=permr)
                   else
                      call self % bin              % coord ((/ii,jj,kk/),comin,comax)
                      call self % octree(ii,jj,kk) % fill  (xx_tmp,bb,memor_loc,comin,comax,permr)                   
                   end if
#else
                   if( if_gap ) then
                      call self % octree(ii,jj,kk) % fill(BOBOX=bb,MEMORY_COUNTER=memor_loc,PERMU=permr)
                   else
                      call self % bin              % coord ((/ii,jj,kk/),comin,comax)
                      call self % octree(ii,jj,kk) % fill  (BOBOX=bb,MEMORY_COUNTER=memor_loc,COORD_MIN=comin,COORD_MAX=comax,PERMU=permr)                   
                   end if
#endif
                   self % num_octree = self % num_octree + 1
                   deallocate(bb,permr,invpr)
                end if
                
                self % stats(1)   = self % stats(1) + self % octree(ii,jj,kk) % stats(1)
                self % stats(2)   = self % stats(2) + self % octree(ii,jj,kk) % stats(2)
                self % stats(3)   = max(self % stats(3),self % octree(ii,jj,kk) % stats(3))
                
             end if
          end do
       end do
    end do

    call self % mem_end(memor_loc,MEMORY_COUNTER)    
    call self % tim_end(SEARCH_FILL)

  end subroutine fill

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   Mesh dimension
  !> @details Find mesh dimension
  !> 
  !-----------------------------------------------------------------------
  
  subroutine mesh_dim(octbin,ndime,mnode,nelem,npoin,CRITERION)

    class(maths_octbin),                    intent(inout) :: octbin
    integer(ip),                            intent(out)   :: ndime
    integer(ip),                            intent(out)   :: mnode
    integer(ip),                            intent(out)   :: nelem
    integer(ip),                            intent(out)   :: npoin
    integer(ip),         optional,          intent(in)    :: CRITERION     !< If empty bins should be considered
    integer(ip)                                           :: ielem,ipoin 
    integer(ip)                                           :: ii,jj,kk
    integer(ip)                                           :: my_criterion

    my_criterion = optional_argument(SEARCH_MESH_LEAVES , CRITERION)
    
    nelem = 0
    npoin = 0
    mnode = (octbin % bin % dim-1)*4
    ndime = octbin % bin % dim
    
    do kk = 1,octbin % bin % boxip(3)
       do jj = 1,octbin % bin % boxip(2)
          do ii = 1,octbin % bin % boxip(1)
             if( octbin % bin % filled((/ii,jj,kk/)) ) then
                call octbin % octree(ii,jj,kk) % mesh_dim(ndime,mnode,ielem,ipoin,CRITERION=my_criterion)
                nelem = nelem + ielem
                npoin = npoin + ipoin
             else if( my_criterion == SEARCH_MESH_LEAVES ) then
                nelem = nelem + 1
                npoin = npoin + mnode
             end if
          end do
       end do
    end do
    
  end subroutine mesh_dim
  
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

    class(maths_octbin),                    intent(inout) :: self
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
    integer(ip)                                           :: ielem,ipoin,pelty
    integer(ip)                                           :: offset_ielem_loc
    integer(ip)                                           :: offset_ipoin_loc
    integer(ip)                                           :: ii,jj,kk
    integer(ip)                                           :: idime,inode
    integer(8)                                            :: memor_loc(2)
    logical(lg)                                           :: if_empty
    logical(lg)                                           :: if_deallocate
    real(rp)                                              :: coord_loc(3,8),rnode
    integer(ip)                                           :: my_criterion

    call self % mem_ini(memor_loc,MEMORY_COUNTER)

    offset_ielem_loc = optional_argument(0_ip               , OFFSET_IELEM)
    offset_ipoin_loc = optional_argument(0_ip               , OFFSET_IPOIN)
    my_criterion     = optional_argument(SEARCH_MESH_LEAVES , CRITERION)
    if_deallocate    = optional_argument(.false.            , ONLY_DEALLOCATE)
    if_empty         = my_criterion == SEARCH_MESH_LEAVES

    if( if_deallocate ) then
       !
       ! Deallocate
       !
       call memory_deallo(memor_loc,'LNODS','mesh',lnods)
       call memory_deallo(memor_loc,'LTYPE','mesh',ltype)
       call memory_deallo(memor_loc,'COORD','mesh',coord)       
       if( present(CENTROID) ) then
          call memory_deallo(memor_loc,'CENTROID','mesh',CENTROID)
       end if
    else
       !
       ! Mesh dimensions
       !
       call self % mesh_dim(ndime,mnode,nelem,npoin,CRITERION)
       rnode = 1.0_rp / real(mnode,rp)
       if( ndime == 2 ) then
          pelty = QUA04
       else
          pelty = HEX08
       end if
       !
       ! Allocate mesh
       !
       if( .not. associated(lnods)    ) call memory_alloca(memor_loc,'LNODS','mesh',lnods,mnode,nelem)
       if( .not. associated(ltype)    ) call memory_alloca(memor_loc,'LTYPE','mesh',ltype,nelem)
       if( .not. associated(coord)    ) call memory_alloca(memor_loc,'COORD','mesh',coord,ndime,npoin)

       if( present(CENTROID) ) then
          if( .not. associated(CENTROID) ) call memory_alloca(memor_loc,'CENTROID','mesh',CENTROID,ndime,nelem)
       end if
       !
       ! Load mesh
       !
       nelem = offset_ielem_loc
       npoin = offset_ipoin_loc
       do kk = 1,self % bin % boxip(3)
          do jj = 1,self % bin % boxip(2)
             do ii = 1,self % bin % boxip(1)
                if( self % bin % filled((/ii,jj,kk/)) ) then
                   call self % octree(ii,jj,kk) % mesh(&
                        ndime,mnode,ielem,ipoin,&
                        lnods,ltype,coord,&
                        MEMORY_COUNTER,&
                        OFFSET_IELEM=nelem,&
                        OFFSET_IPOIN=npoin,&
                        CRITERION=SEARCH_MESH_LEAVES,&
                        CENTROID=CENTROID)
                   nelem = nelem + ielem 
                   npoin = npoin + ipoin

                else if( if_empty ) then

                   call self % bin % mesh_coord((/ii,jj,kk/),coord_loc)
                   ipoin = npoin+1
                   ielem = nelem+1
                   do idime = 1,ndime
                      coord(idime,ipoin:ipoin+mnode-1) = coord_loc(idime,1:mnode)                   
                   end do
                   if( present(CENTROID) ) then
                      do idime = 1,ndime
                         CENTROID(idime,ielem) = sum(coord_loc(idime,1:mnode)) * rnode
                      end do
                   end if
                   ipoin = npoin+1
                   ltype(ielem) = pelty
                   do inode = 1,mnode
                      lnods(inode,ielem) = ipoin
                      ipoin = ipoin + 1
                   end do
                   nelem = nelem + 1
                   npoin = npoin + mnode

                end if
             end do
          end do
       end do
    end if

    call self % mem_end(memor_loc,MEMORY_COUNTER)

  end subroutine mesh

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   Centroid
  !> @details Get teh centroid of each box
  !> 
  !-----------------------------------------------------------------------
  
  subroutine centroid(octbin,coorc,MEMORY_COUNTER)

    class(maths_octbin),                    intent(inout) :: octbin
    real(rp),                      pointer, intent(inout) :: coorc(:,:)
    integer(8),          optional,          intent(inout) :: MEMORY_COUNTER(2)
    integer(ip)                                           :: nelem,ndime
    integer(ip)                                           :: ii,jj,kk
    integer(8)                                            :: memor_loc(2)

    call octbin % mem_ini(memor_loc,MEMORY_COUNTER)

    nelem = 0
    ndime = octbin % bin % dim
    
    do kk = 1,octbin % bin % boxip(3)
       do jj = 1,octbin % bin % boxip(2)
          do ii = 1,octbin % bin % boxip(1)
             if( associated(octbin % bin % list(ii,jj,kk) % l) ) then
                nelem = nelem + octbin % octree(ii,jj,kk) % nleaves
             else
                nelem = nelem + 1
             end if
          end do
       end do
    end do

    if( .not. associated(coorc) ) call memory_alloca(memor_loc,'COORC','mesh',coorc,ndime,nelem)

    nelem = 0
    do kk = 1,octbin % bin % boxip(3)
       do jj = 1,octbin % bin % boxip(2)
          do ii = 1,octbin % bin % boxip(1)
             if( octbin % bin % filled((/ii,jj,kk/)) ) then
                call octbin % octree(ii,jj,kk) % centroid(coorc,MEMORY_COUNTER,OFFSET_IELEM=nelem)
                nelem = nelem + octbin % octree(ii,jj,kk) % nleaves
             else
                
             end if
          end do
       end do
    end do

    call octbin % mem_end(memor_loc,MEMORY_COUNTER)

  end subroutine centroid

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-10
  !> @brief   Results
  !> @details Get some results
  !> 
  !-----------------------------------------------------------------------

  subroutine results(self,xx,names,OFFSET,MEMORY_COUNTER,CRITERION,ONLY_DEALLOCATE)

    class(maths_octbin),                 intent(inout) :: self
    real(rp),                   pointer, intent(inout) :: xx(:,:)
    character(len=5),           pointer, intent(inout) :: names(:)
    integer(ip),      optional,          intent(in)    :: OFFSET
    integer(8),       optional,          intent(inout) :: MEMORY_COUNTER(2)
    integer(ip),      optional,          intent(in)    :: CRITERION    !< If empty bins should be considered
    logical(lg),      optional,          intent(in)    :: ONLY_DEALLOCATE
    integer(ip)                                        :: ii,jj,kk
    integer(ip)                                        :: nelem
    integer(ip)                                        :: mnode,npoin,ndime
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
       ! Allocate results and fill in results
       !
       call self % mesh_dim(ndime,mnode,nelem,npoin,CRITERION=my_criterion)

       if( .not. associated(xx)    ) call memory_alloca(memor_loc,'XX'   ,'results',xx,nelem,4_ip)
       if( .not. associated(names) ) call memory_alloca(memor_loc,'NAMES','results',names,4_ip)
       
       nelem = ielem_offset
       do kk = 1,self % bin % boxip(3)
          do jj = 1,self % bin % boxip(2)
             do ii = 1,self % bin % boxip(1)
                if( self % bin % filled((/ii,jj,kk/)) ) then
                   call self % octree(ii,jj,kk) % results(xx,names,OFFSET=nelem,CRITERION=my_criterion)!if_empty)
                   !if( if_empty ) then
                      nelem = nelem + self % octree(ii,jj,kk) % nleaves
                   !else
                   !   nelem = nelem + self % octree(ii,jj,kk) % nfilled
                   !end if
                else if( my_criterion == SEARCH_MESH_LEAVES ) then
                   nelem = nelem + 1
                end if
             end do
          end do
       end do
       
    end if

    call self % mem_end(memor_loc,MEMORY_COUNTER)
    
  end subroutine results

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    24/01/2022
  !> @brief   Search type name
  !> @details Search type name
  !
  !-----------------------------------------------------------------------
  
  function type_name(self) result(name)
    class(maths_octbin),                    intent(inout) :: self
    character(LEN=:), allocatable                         :: name
    
    name = 'OCTBIN'
    
  end function type_name

end module def_maths_octbin
!> @}
