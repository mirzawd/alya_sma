!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!>
!> @addtogroup Element_Search_Strategy_Toolbox
!> @{
!> @name    ToolBox for elemental search
!> @file    mod_elsest.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for elements
!> @details Different functions to find host element
!>          \verbatim
!>
!>          EL      SE     ST
!>            ement   arch   rategies   L I B R A R Y
!>
!>          INPUT
!>          -----
!>          IPARA( 1) = # of bins in x                 (BIN)
!>          IPARA( 2) = # of bins in y                 (BIN)
!>          IPARA( 3) = # of bins in z                 (BIN)
!>          IPARA( 4) = Data format (0=type,1=list)    (BIN)
!>          IPARA( 5) = Max # of background meshes     (BOTH)
!>          IPARA( 7) = output unit (0 dfor no output) (BOTH)
!>          IPARA( 8) = Search strategy (0=bin,1=oct)
!>          IPARA( 9) = Max # of nodes per bin         (QUAD/OCT)
!>          IPARA(10) = Output frequency (not coded)
!>          IPARA(11) = search radius of
!>                    = 0 for no                       (BIN)
!>                    = 1 for yes                      (BIN)
!>          IPARA(12) = Unit for mesh postprocess      (BOTH)
!>          IPARA(13) = Unit for result postprocess    (BOTH)
!>          IPARA(14) = If element flag should be      (BOTH)
!>                      checked during search          
!>          IPARA(15) = Find an element anyway         (BOTH)
!>                      by projection
!>          IPARA(16) = Save element bounding box      (BIN)
!>
!>          RPARA(1)  = Tolerance
!> 
!>          OUTPUT
!>          ------
!>          IFOUN     = 0: Element not found
!>                    > 0: Host element
!>          RPARA(2)  = Memory
!>          RPARA(3)  = Max. Memory
!>         
!>          \endverbatim
!>
!-----------------------------------------------------------------------

module mod_elsest

  use def_kintyp,               only : ip,rp,lg,i1p
  use def_domain,               only : mesh_type
  use mod_elmgeo,               only : elmgeo_natural_coordinates
  use mod_elmgeo,               only : elmgeo_nearest_point_on_element_faces
  use mod_elmgeo,               only : element_type
  use mod_elmgeo,               only : elmgeo_nearest_intersection_point_on_element_faces
  use mod_elmgeo,               only : elmgeo_nearest_element_node
  use mod_elmgeo,               only : elmgeo_inside_element_bounding_box
  use def_elmtyp,               only : BAR02,BAR03,BAR04,TRI03,TRI06,QUA04,QUA08,QUA09,QUA16
  use def_elmtyp,               only : TET04,TET10,PYR05,PYR14,PEN06,PEN15,PEN18,HEX08
  use def_elmtyp,               only : HEX27,HEX64,SHELL,BAR3D
  use mod_memory,               only : memory_alloca
  use mod_memory,               only : memory_deallo
  use mod_memory,               only : memory_size
  use def_search_strategy,      only : search
  use def_search_strategy,      only : SEARCH_BIN
  use def_search_strategy,      only : SEARCH_OCTREE
  use def_interpolation_method, only : interpolation
  implicit none
  private 


  integer(ip) :: method = 0 ! 0=old, 1=new
  
  type octbox
     integer(ip)               :: id          ! My global ID
     integer(ip)               :: level       ! Generation
     integer(ip)               :: npoinbox    ! Number of nodes
     integer(ip)               :: nelembox    ! Number of elements
     integer(ip)               :: childid     ! Child ID (1->4 or 1->8)
     integer(ip)               :: whoiam      ! Father or have nodes
     integer(ip),  pointer     :: nodes(:)    
     integer(ip),  pointer     :: elems(:)    ! List of elements
     real(rp)                  :: minc(3)     ! Min coordinates
     real(rp)                  :: maxc(3)     ! Max coordinates
     type(octbox), pointer     :: parent      ! Pointer to parent
     type(octbox), pointer     :: children(:) ! Pointer to children
  end type octbox
  !
  ! Bin structures 
  !
  type bintype
     integer(ip)           :: nboxe
     integer(ip)           :: nboxx(3)
     integer(ip)           :: dataf
     integer(ip)           :: iallo  
     integer(ip), pointer  :: lboel(:)
     integer(ip), pointer  :: pboel(:)
     type(i1p),   pointer  :: tboel(:)
     integer(ip), pointer  :: kstat(:)    
     real(rp)              :: delta(3)
     real(rp)              :: comin(3)
     real(rp)              :: comax(3)
     real(rp),    pointer  :: cputi(:)    
     real(rp),    pointer  :: element_bb(:,:,:)
  end type bintype
  !
  ! Oct-tree structure
  !
  type octtype
     integer(ip)           :: iallo
     type(octbox),pointer  :: tree_root
     integer(ip), pointer  :: kstat(:)  
     integer(ip)           :: divmax
     real(rp)              :: comin(3)
     real(rp)              :: comax(3)
     real(rp),    pointer  :: cputi(:)    
     real(rp),    pointer  :: element_bb(:,:,:)
  end type octtype
  !
  ! Elsest structures
  !  
  type(bintype),   pointer  :: bin_struc(:)  => null()    
  type(octtype),   pointer  :: oct_struc(:)  => null()

  type(octbox)   :: octbox_init = octbox(&
       0_ip,&                      ! id         
       0_ip,&                      ! level      
       0_ip,&                      ! npoinbox   
       0_ip,&                      ! nelembox   
       0_ip,&                      ! childid    
       0_ip,&                      ! whoiam     
       null(),&                    ! nodes(:)   
       null(),&                    ! elems(:)   
       (/0.0_rp,0.0_rp,0.0_rp/),&  ! minc(3)    
       (/0.0_rp,0.0_rp,0.0_rp/),&  ! maxc(3)    
       null(),&                    ! parent     
       null())                     ! children(:)
  type(bintype),  parameter :: bin_struc_init=bintype(&
       0_ip,&                      ! nboxe         
       (/0_ip,0_ip,0_ip/),&        ! nboxx(3)      
       0_ip,&                      ! dataf         
       0_ip,&                      ! iallo         
       null(),&                    ! lboel(:)      
       null(),&                    ! pboel(:)      
       null(),&                    ! tboel(:)      
       null(),&                    ! kstat(:)      
       (/0.0_rp,0.0_rp,0.0_rp/),&  ! delta(3)      
       (/0.0_rp,0.0_rp,0.0_rp/),&  ! comin(3)      
       (/0.0_rp,0.0_rp,0.0_rp/),&  ! comax(3)  
       null(),&                    ! cputi(:)
       null())                     ! Bounding boxes
  type(octtype), parameter :: oct_struc_init=octtype(&
       0_ip,&                      ! iallo      
       null(),&                    ! tree_root              
       null(),&                    ! kstat(:)             
       0_ip,&                      ! divmax             
       0.0_rp,&                    ! comin(3)        
       0.0_rp,&                    ! comax(3)     
       null(),&                    ! comax(3)     
       null())                     ! Bounding boxes

  integer(ip), parameter :: ELSEST_BIN_STRATEGY              = 0
  integer(ip), parameter :: ELSEST_OCT_TREE_STRATEGY         = 1
  integer(ip), parameter :: ELSEST_KD_TREE_STRATEGY          = 3
  integer(ip), parameter :: ELSEST_A_TODA_COSTA              = 2
  integer(ip), parameter :: ELSEST_ELEMENT_NOT_FOUND         = 0
  integer(ip), parameter :: ELSEST_MASK_HOST_ELEMENT         = 1
  integer(ip), parameter :: ELSEST_TYPE_BIN_STRUCTURE        = 0
  integer(ip), parameter :: ELSEST_LINKED_LIST_BIN_STRUCTURE = 1
  real(rp),    parameter :: zeror                            = epsilon(1.0_rp)

  integer(ip)            :: iunit(3)  ! Output units
  integer(8)             :: memor8(2) ! Memory counter

  integer(ip)            :: num_elem
  real(rp)               :: time1,time2,time_total=0.0_rp

  type(search)           :: search_elsest_seq             ! Elsest: Sequential search
  type(search)           :: search_elsest_par             ! Elsest: Parallel search
  type(interpolation)    :: interp_elsest                 ! Elsest: Wall exchange iterpolation


  public :: time_total
  public :: num_elem

  public :: elsest_initialization
  public :: elsest_preprocess
  public :: elsest_host_element
  public :: elsest_deallocate
  public :: elsest_allocate

  public :: ELSEST_A_TODA_COSTA

  public :: search_elsest_seq
  public :: search_elsest_par
  public :: interp_elsest

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-09
  !> @brief   Elsest deallocate
  !> @details Elsest deallocate
  !> 
  !-----------------------------------------------------------------------

  subroutine elsest_deallocate(ipara)

    integer(ip), intent(in) :: ipara(*)  !< Integer parameters
    integer(ip)             :: imesh

    if( method == 1 ) then

       call search_elsest_seq % deallo ()
       call search_elsest_par % deallo ()

    else

       imesh = 1

       if( ipara(8) == ELSEST_BIN_STRATEGY ) then                           
          !
          ! Bin strategy
          !
          if( associated(bin_struc) ) then
             call elsest_bin_deallocate(imesh)
             deallocate(bin_struc) 
          end if

       else if( ipara(8) == ELSEST_OCT_TREE_STRATEGY ) then
          !
          ! Quad/Oct strategy
          !
          if( associated(oct_struc) ) then
             call elsest_oct_deallocate(imesh)
             deallocate(oct_struc)
          end if

       else if( ipara(8) == ELSEST_KD_TREE_STRATEGY ) then
          !
          ! Kd-tree strategy
          !
       end if
    end if

  end subroutine elsest_deallocate

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-09
  !> @brief   Elsest initialization
  !> @details Elsest initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine elsest_initialization(ielse)

    integer(ip), intent(inout) :: ielse(:)

    ielse(1)      = 100                                    ! nx
    ielse(2)      = 100                                    ! ny
    ielse(3)      = 100                                    ! nz
    ielse(4)      = 1                                      ! data format (0=type,1=linked list)
    ielse(5)      = 10                                     ! Maximum number of possible meshes
    ielse(6)      = 1                                      ! Not used 
    ielse(7)      = 0                                      ! Output unit
    ielse(8)      = 1                                      ! Search strategy (0=bin,1=Quad)
    ielse(9)      = 100                                    ! Points per node for Quad/Octtree
    ielse(10)     = 0                                      ! Result output frequency 
    ielse(11)     = 1                                      ! Neighboring-boxes-search radius, 0: search until all boxes finished
    ielse(12)     = 0                                      ! Postprocess mesh
    ielse(13)     = 0                                      ! Postprocess results
    ielse(14)     = 0                                      ! Dont check
    ielse(15)     = 0                                      ! Dont force
    memor8        = 0_8                                    ! Memory counter

    nullify(bin_struc)
    nullify(oct_struc)

  end subroutine elsest_initialization

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    18/11/2015
  !> @brief   Find a host element using Elsest
  !> @details Find a host element using Elsest
  !>
  !-----------------------------------------------------------------------

  subroutine elsest_host_element(&
       ipara,rpara,imesh,meshe,xcoor,ifoun,shapt,&
       derit,coloc,dista,lchec)

    integer(ip),     intent(in)            :: ipara(*)  !< Integer parameters
    real(rp),        intent(in)            :: rpara(*)  !< Real parameters
    integer(ip),     intent(in)            :: imesh     !< Mesh number
    type(mesh_type), intent(in)            :: meshe     !< Mesh type
    real(rp),        intent(in)            :: xcoor(*)  !< Test point coordinates
    integer(ip),     intent(out)           :: ifoun     !< Host element
    real(rp),        intent(out)           :: shapt(*)  !< Shape funciton in host element
    real(rp),        intent(out)           :: derit(*)  !< Shae funciton derivatives
    real(rp),        intent(out)           :: coloc(*)  !< Parametric coordinates in host element
    real(rp),        intent(out)           :: dista     !< Distance to element
    integer(ip),     intent(in), optional  :: lchec(*)  !< List of elements to be checked

    integer(ip)                            :: inode,ielem,kelem,ilook
    integer(ip)                            :: pnode,pelty,ipoin
    real(rp)                               :: elcod(meshe % ndime,meshe % mnode)
    real(rp)                               :: toler
    real(rp),        pointer               :: xx(:,:)
    type(i1p),       pointer               :: list_entities(:)
    logical(lg)                            :: if_check
    !
    ! Check errors
    !
    if( ipara(14) == ELSEST_MASK_HOST_ELEMENT .and. .not. present(lchec) ) then
       call runend('ELSEST_HOST_ELEMENT: MASK ARRAY IS MISSING')
    end if
    !
    ! Look for host element
    !
    if( method == 1 ) then

       nullify (xx,list_entities)
       allocate(xx(meshe % ndime,1))
       xx(1:meshe % ndime,1) = xcoor(1:meshe % ndime)
       call search_elsest_seq % method % candidate(xx,list_entities)

       if( ipara(14) /= 0 .and. present(lchec) ) then
          if_check = .true.
       else
          if_check = .false.
       end if
       
       ifoun = 0
       dista = huge(1.0_rp)
       toler = abs(rpara(1))

       loop_kelem: do kelem = 1,memory_size(list_entities(1) % l)
          
          ielem = list_entities(1) % l(kelem)
          pelty = meshe % ltype(ielem)
          
          ilook = 1
          if( if_check ) then
             if( lchec(ielem) /= ipara(14) ) ilook = 0
          end if
          
          if( ilook == 1 .and. pelty > 0 ) then
             
             pnode = meshe % lnnod(ielem)
             do inode = 1,pnode
                ipoin= meshe % lnods(inode,ielem)
                elcod(1:meshe % ndime,inode) = meshe % coord(1:meshe % ndime,ipoin)
             end do
             
             call elmgeo_natural_coordinates(&
                  meshe % ndime,pelty,pnode,elcod,&
                  shapt,derit,xx(:,1),coloc,ifoun,toler)
             
             if( ifoun > 0 ) then                   
                ifoun = ielem
                dista = 0.0_rp
                exit loop_kelem                
             end if
             
          end if
             
       end do loop_kelem       
       
       call memory_deallo(search_elsest_seq % method % memor,'LIST_ENTITIES','mod_elsest',list_entities)
       deallocate(xx)

    else
       
       if( ipara(8) == ELSEST_BIN_STRATEGY ) then

          call elsest_bin_host_element(&                                           ! Bin strategy
               imesh,ipara,rpara,&
               meshe % mnode,meshe % ndime,meshe % npoin,&
               meshe % nelem,meshe % lnnod,meshe % lnods,&
               meshe % ltype,meshe % coord,xcoor,ifoun,shapt,&
               derit,coloc,dista,lchec)

       else if( ipara(8) == ELSEST_OCT_TREE_STRATEGY ) then

          call elsest_oct_host_element(&                                           ! Oct-tree strategy
               imesh,ipara,rpara,&
               meshe % mnode,meshe % ndime,meshe % npoin,&
               meshe % nelem,meshe % lnnod,meshe % lnods,&
               meshe % ltype,meshe % coord,xcoor,ifoun,shapt,&
               derit,coloc,dista,lchec)

       else if( ipara(8) == ELSEST_KD_TREE_STRATEGY ) then

          call elsest_kd_host_element(&                                            ! Kd-tree strategy
               imesh,ipara,rpara,&
               meshe % mnode,meshe % ndime,meshe % npoin,&
               meshe % nelem,meshe % lnnod,meshe % lnods,&
               meshe % ltype,meshe % coord,xcoor,ifoun,shapt,&
               derit,coloc,dista,lchec)

       end if

       if( ( ifoun == ELSEST_ELEMENT_NOT_FOUND .and. ipara(15) == 1 ) .or. ipara(8) == ELSEST_A_TODA_COSTA ) then
          call elsest_host_element_a_toda_costa(&
               imesh,ipara,rpara,&
               meshe % mnode,meshe % ndime,meshe % npoin,&
               meshe % nelem,meshe % lnnod,meshe % lnods,&
               meshe % ltype,meshe % coord,xcoor,ifoun,shapt,&
               derit,coloc,dista,lchec)
       end if

    end if

  end subroutine elsest_host_element

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    18/11/2015
  !> @brief   Initialize Elsest
  !> @details Initialize Elsest for NMESH meshes and construct oct or
  !>          bin trees
  !>
  !-----------------------------------------------------------------------

  subroutine elsest_preprocess(ipara,rpara,meshe,CURRENT_MESH)

    integer(ip),               intent(in)    :: ipara(*)        !< Integer parameters
    real(rp),                  intent(inout) :: rpara(*)        !< Real parameters
    type(mesh_type),           intent(in)    :: meshe           !< Mesh type
    integer(ip),     optional, intent(in)    :: CURRENT_MESH    !< Current mesh to allocate
    integer(ip)                              :: imesh
    real(rp),        pointer                 :: bobox(:,:,:)

    if( present(CURRENT_MESH) ) then
       imesh = CURRENT_MESH
    else
       imesh = 1
    end if

    if( method == 1 ) then
       
       nullify(bobox)
       call meshe             % element_bb (bobox,MEMORY_COUNTER=memor8)
       call search_elsest_seq % set        (BOBOX=bobox,NAME='ELSEST')
       call memory_deallo(memor8,'BOBOX','mod_elsest',bobox)

       
    else

       if( meshe % nelem == 0 ) return

       if( ipara(8) == ELSEST_BIN_STRATEGY ) then                           
          !
          ! Bin strategy
          !        
          call elsest_bin_preprocess(&
               ipara,rpara,imesh,&
               meshe % mnode,meshe % ndime,meshe % npoin,meshe % nelem,&
               meshe % lnnod,meshe % lnods,meshe % ltype,meshe % coord)

       else if( ipara(8) == ELSEST_OCT_TREE_STRATEGY ) then
          !
          ! Quad/Oct strategy
          !
          call elsest_oct_preprocess(&      
               ipara,rpara,imesh,&
               meshe % mnode,meshe % ndime,meshe % npoin,meshe % nelem,&
               meshe % lnnod,meshe % lnods,meshe % ltype,meshe % coord)

       else if( ipara(8) == ELSEST_KD_TREE_STRATEGY ) then                           

          call runend('MOD_ELSEST: KDTREE NO LONGER SUPPORTED')

       end if

    end if
    
    rpara(2) = real(memor8(1),rp)
    rpara(3) = real(memor8(2),rp)

  end subroutine elsest_preprocess

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    18/11/2015
  !> @brief   Bin preprocess
  !>
  !-----------------------------------------------------------------------

  subroutine elsest_bin_preprocess(&
       ipara,rpara,imesh,mnode,ndime,npoin,nelem,lnnod,&
       lnods,ltype,coord)

    integer(ip), intent(in)  :: ipara(*)
    real(rp),    intent(in)  :: rpara(*)
    integer(ip), intent(in)  :: imesh
    integer(ip), intent(in)  :: mnode
    integer(ip), intent(in)  :: ndime
    integer(ip), intent(in)  :: npoin
    integer(ip), intent(in)  :: nelem
    integer(ip), intent(in)  :: lnnod(nelem)
    integer(ip), intent(in)  :: lnods(mnode,nelem)
    integer(ip), intent(in)  :: ltype(nelem)
    real(rp),    intent(in)  :: coord(ndime,npoin)
    integer(ip)              :: idime,iboxe,kk,jj,ii,ll
    integer(ip)              :: imin,imax,jmin,jmax,kmin,kmax,ni,nj,nk
    integer(ip)              :: ielem,box_nr,box_nr1,box_nr2,ninj
    real(rp)                 :: time1,time2,time3,time4,rni,rnj,rnk
    real(rp)                 :: deltx,delty,deltz,dni,dnj,dnk

    real(rp),     pointer    :: xmima(:,:,:)
    integer(ip),  pointer    :: nbono(:)

    integer(ip),    pointer  :: nboxe
    integer(ip),    pointer  :: nboxx(:)
    integer(ip),    pointer  :: dataf
    integer(ip),    pointer  :: lboel(:)
    integer(ip),    pointer  :: pboel(:)
    type(i1p),      pointer  :: tboel(:)
    integer(ip),    pointer  :: kstat(:)    
    real(rp),       pointer  :: delta(:)
    real(rp),       pointer  :: comin(:)
    real(rp),       pointer  :: comax(:)
    real(rp),       pointer  :: cputi(:)   

    call elsest_cputim(time1)
    bin_struc(imesh) % iallo = 1
    nullify(xmima)
    nullify(nbono)

    !----------------------------------------------------------------------
    !
    ! Initialize parameters
    !
    !----------------------------------------------------------------------

    call memory_alloca(memor8,'BIN_STRUC % CPUTI'     ,'elsest_bin_preprocess',bin_struc(imesh) % cputi,10_ip)
    call memory_alloca(memor8,'BIN_STRUC % KSTAT'     ,'elsest_bin_preprocess',bin_struc(imesh) % kstat,10_ip)
    call memory_alloca(memor8,'BIN_STRUC % ELEMENT_BB','elsest_bin_preprocess',bin_struc(imesh) % element_bb,ndime,2_ip,nelem)

    !----------------------------------------------------------------------
    !
    ! Point to current mesh (IMESH) structure
    !
    !----------------------------------------------------------------------

    nboxe => bin_struc(imesh) % nboxe
    nboxx => bin_struc(imesh) % nboxx
    dataf => bin_struc(imesh) % dataf
    lboel => bin_struc(imesh) % lboel
    pboel => bin_struc(imesh) % pboel
    tboel => bin_struc(imesh) % tboel
    kstat => bin_struc(imesh) % kstat
    comin => bin_struc(imesh) % comin
    comax => bin_struc(imesh) % comax
    delta => bin_struc(imesh) % delta
    cputi => bin_struc(imesh) % cputi
    xmima => bin_struc(imesh) % element_bb

    cputi    = 0.0_rp
    kstat    = 0_ip
    kstat(1) = huge(1_ip)
    kstat(3) = huge(1_ip)
    nboxx(1) = ipara(1)
    nboxx(2) = ipara(2)
    nboxx(3) = ipara(3)  
    dataf    = ipara(4)
    comin(1) = 0.0_rp
    comin(2) = 0.0_rp 
    comin(3) = 0.0_rp
    comax(1) = 0.0_rp
    comax(2) = 0.0_rp
    comax(3) = 0.0_rp

    !----------------------------------------------------------------------
    !
    ! Total number of boxes
    !
    !----------------------------------------------------------------------

    if( ndime == 1 ) then
       nboxe = nboxx(1)
    else if( ndime == 2 ) then
       nboxe = nboxx(1) * nboxx(2)
    else
       nboxe = nboxx(1) * nboxx(2) * nboxx(3)
    end if
    ni   = nboxx(1)
    nj   = nboxx(2)
    nk   = nboxx(3)
    bin_struc(imesh) % nboxx(1) = ni
    bin_struc(imesh) % nboxx(2) = nj
    bin_struc(imesh) % nboxx(3) = nk
    rni = real(ni,rp)
    rnj = real(nj,rp)
    rnk = real(nk,rp)

    !----------------------------------------------------------------------
    !
    ! Allocate memory for bin structure
    !
    !----------------------------------------------------------------------

    if( dataf == ELSEST_TYPE_BIN_STRUCTURE ) then
       call memory_alloca(memor8,'BIN_STRUC % TBOEL','elsest_bin_preprocess',bin_struc(imesh) % tboel,nboxe)
    else if( dataf == ELSEST_LINKED_LIST_BIN_STRUCTURE ) then
       call memory_alloca(memor8,'BIN_STRUC % PBOEL','elsest_bin_preprocess',bin_struc(imesh) % pboel,nboxe+1_ip)
    end if

    tboel => bin_struc(imesh) % tboel
    pboel => bin_struc(imesh) % pboel
    lboel => bin_struc(imesh) % lboel

    !----------------------------------------------------------------------
    !
    ! Compute bounding box and bin spacing delta
    !
    !----------------------------------------------------------------------

    call elsest_bounding_box(ndime,npoin,coord,comin,comax)
    comin = comin - zeror
    comax = comax + zeror
    !
    ! Bin spacing DELTA
    !
    do idime = 1,ndime
       delta(idime) = (comax(idime)-comin(idime))/real(nboxx(idime),rp)
    end do
    deltx = 1.0_rp / ( comax(1)-comin(1) )
    dni   = real(ni,rp) * deltx
    if(ndime >= 2 ) then
       delty = 1.0_rp / ( comax(2)-comin(2) )
       dnj   = real(nj,rp) * delty
       ninj  = nboxx(1) * nboxx(2)
    end if
    if( ndime == 3 ) then
       deltz = 1.0_rp / ( comax(3)-comin(3) )
       dnk   = real(nk,rp) * deltz
    end if

    call elsest_cputim(time2)
    cputi(1)=time2-time1

    !----------------------------------------------------------------------
    !
    ! Compute element bounding boxes
    !
    !----------------------------------------------------------------------

    call elsest_element_bounding_boxes(&
         mnode,ndime,npoin,nelem,lnnod,lnods,coord,bin_struc(imesh) % element_bb)
    xmima => bin_struc(imesh) % element_bb

    !call memory_alloca(memor8,'XMIMA','elsest_bin_preprocess',xmima,3_ip,2_ip,nelem)
    !call elsest_element_bounding_box(mnode,ndime,npoin,nelem,lnnod,lnods,coord,xmima)

    !----------------------------------------------------------------------
    !
    ! Number of elements per box
    !
    !----------------------------------------------------------------------

    call memory_alloca(memor8,'NBONO','elsest_bin_preprocess',nbono,nboxe)

    if( ndime == 1 ) then

       do ielem = 1,nelem

          imin = int(( ( xmima(1,1,ielem) - comin(1)) * deltx ) * rni , ip ) + 1
          imax = int(( ( xmima(1,2,ielem) - comin(1)) * deltx ) * rni , ip ) + 1  

          imin = max(imin,1_ip)
          imax = min(imax,nboxx(1))

          do ii = imin,imax
             box_nr = ii
             nbono(box_nr) = nbono(box_nr) + 1
          end do
       end do

    else if( ndime == 2 ) then

       do ielem = 1,nelem

          imin = int(( ( xmima(1,1,ielem) - comin(1) ) * deltx ) * rni , ip ) + 1
          imax = int(( ( xmima(1,2,ielem) - comin(1) ) * deltx ) * rni , ip ) + 1      

          jmin = int(( ( xmima(2,1,ielem) - comin(2) ) * delty ) * rnj , ip ) + 1
          jmax = int(( ( xmima(2,2,ielem) - comin(2) ) * delty ) * rnj , ip ) + 1      

          imin = max(imin,1_ip)
          imax = min(imax,nboxx(1))
          jmin = max(jmin,1_ip)
          jmax = min(jmax,nboxx(2))

          do ii = imin,imax
             do jj = jmin,jmax
                box_nr = (jj-1_ip) * ni + ii
                nbono(box_nr) = nbono(box_nr) + 1
             end do
          end do
       end do

    else if( ndime == 3 ) then

       do ielem = 1,nelem

          imin = int( ( xmima(1,1,ielem) - comin(1) ) * dni , ip )  + 1
          imax = int( ( xmima(1,2,ielem) - comin(1) ) * dni , ip )  + 1      

          jmin = int( ( xmima(2,1,ielem) - comin(2) ) * dnj , ip )  + 1
          jmax = int( ( xmima(2,2,ielem) - comin(2) ) * dnj , ip )  + 1      

          kmin = int( ( xmima(3,1,ielem) - comin(3) ) * dnk , ip )  + 1
          kmax = int( ( xmima(3,2,ielem) - comin(3) ) * dnk , ip )  + 1      

          imin = max(imin,1_ip)
          imax = min(imax,nboxx(1))
          jmin = max(jmin,1_ip)
          jmax = min(jmax,nboxx(2))
          kmin = max(kmin,1_ip)
          kmax = min(kmax,nboxx(3))

          do kk = kmin,kmax
             box_nr2 = ninj * (kk-1)
             do jj = jmin,jmax
                box_nr1 = box_nr2 + nboxx(1) * (jj-1)
                do ii = imin,imax
                   box_nr        = box_nr1 + ii
                   nbono(box_nr) = nbono(box_nr) + 1
                end do
             end do
          end do

       end do

    end if

    !----------------------------------------------------------------------
    !
    ! Fill in box element list
    !
    !----------------------------------------------------------------------

    if( dataf == ELSEST_TYPE_BIN_STRUCTURE ) then
       !
       ! Type
       !
       do iboxe = 1,bin_struc(imesh) % nboxe
          call memory_alloca(memor8,'BIN_STRUC % TBOEL % L','elsest_bin_preprocess',bin_struc(imesh) % tboel(iboxe) % l,nbono(iboxe))
          nbono(iboxe) = 0
       end do

       if( ndime == 1 ) then

          do ielem = 1,nelem

             imin = int(( ( xmima(1,1,ielem) - comin(1)) * deltx ) * rni , ip ) + 1
             imax = int(( ( xmima(1,2,ielem) - comin(1)) * deltx ) * rni , ip ) + 1      

             imin = max(imin,1_ip)
             imax = min(imax,nboxx(1))

             do ii = imin,imax
                box_nr        = ii
                nbono(box_nr) = nbono(box_nr) + 1
                bin_struc(imesh) % tboel(box_nr) % l ( nbono(box_nr) ) = ielem
             end do

          end do

       else if( ndime == 2 ) then

          do ielem = 1,nelem

             imin = int(( ( xmima(1,1,ielem) - comin(1)) * deltx ) * rni , ip ) + 1
             imax = int(( ( xmima(1,2,ielem) - comin(1)) * deltx ) * rni , ip ) + 1      

             jmin = int(( ( xmima(2,1,ielem) - comin(2)) * delty ) * rnj , ip ) + 1
             jmax = int(( ( xmima(2,2,ielem) - comin(2)) * delty ) * rnj , ip ) + 1    

             imin = max(imin,1_ip)
             imax = min(imax,nboxx(1))
             jmin = max(jmin,1_ip)
             jmax = min(jmax,nboxx(2))

             do ii = imin,imax
                do jj = jmin,jmax
                   box_nr        = (jj-1_ip) * ni + ii
                   nbono(box_nr) = nbono(box_nr) + 1
                   bin_struc(imesh) % tboel(box_nr) % l ( nbono(box_nr) ) = ielem
                end do

             end do

          end do

       else

          do ielem = 1,nelem

             imin = int( ( xmima(1,1,ielem) - comin(1) ) * dni , ip ) + 1
             imax = int( ( xmima(1,2,ielem) - comin(1) ) * dni , ip ) + 1      

             jmin = int( ( xmima(2,1,ielem) - comin(2) ) * dnj , ip ) + 1
             jmax = int( ( xmima(2,2,ielem) - comin(2) ) * dnj , ip ) + 1      

             kmin = int( ( xmima(3,1,ielem) - comin(3) ) * dnk , ip ) + 1
             kmax = int( ( xmima(3,2,ielem) - comin(3) ) * dnk , ip ) + 1      

             imin = max(imin,1_ip)
             imax = min(imax,nboxx(1))
             jmin = max(jmin,1_ip)
             jmax = min(jmax,nboxx(2))
             kmin = max(kmin,1_ip)
             kmax = min(kmax,nboxx(3))

             do kk = kmin,kmax
                box_nr2 = ninj * (kk-1)
                do jj = jmin,jmax
                   box_nr1 = box_nr2 + nboxx(1) * (jj-1)
                   do ii = imin,imax
                      box_nr        = box_nr1 + ii
                      nbono(box_nr) = nbono(box_nr) + 1
                      bin_struc(imesh) % tboel(box_nr) % l ( nbono(box_nr) ) = ielem
                   end do
                end do

             end do

          end do

       end if

    else
       !
       ! Linked list
       !
       bin_struc(imesh) % pboel(1) = 1

       do iboxe = 1,bin_struc(imesh) % nboxe
          bin_struc(imesh) % pboel(iboxe+1) = bin_struc(imesh) % pboel(iboxe) + nbono(iboxe) 
          nbono(iboxe) = 0
       end do

       call memory_alloca(memor8,'BIN_STRUC % LBOEL','elsest_bin_preprocess',bin_struc(imesh)%lboel,bin_struc(imesh) % pboel(bin_struc(imesh) % nboxe+1_ip))
       lboel => bin_struc(imesh) % lboel

       if( ndime == 1 ) then

          do ielem = 1,nelem

             imin = int( ( xmima(1,1,ielem) - comin(1)) * dni , ip ) + 1
             imax = int( ( xmima(1,2,ielem) - comin(1)) * dni , ip ) + 1      

             imin = max(imin,1_ip)
             imax = min(imax,nboxx(1))

             do ii = imin,imax
                box_nr        = ii
                ll            = bin_struc(imesh) % pboel(box_nr) + nbono(box_nr)
                nbono(box_nr) = nbono(box_nr) + 1
                lboel(ll)     = ielem
             end do
          end do

       else if( ndime == 2 ) then

          do ielem = 1,nelem

             imin = int( ( xmima(1,1,ielem) - comin(1)) * dni , ip ) + 1
             imax = int( ( xmima(1,2,ielem) - comin(1)) * dni , ip ) + 1      

             jmin = int( ( xmima(2,1,ielem) - comin(2)) * dnj , ip ) + 1
             jmax = int( ( xmima(2,2,ielem) - comin(2)) * dnj , ip ) + 1      

             imin = max(imin,1_ip)
             imax = min(imax,nboxx(1))
             jmin = max(jmin,1_ip)
             jmax = min(jmax,nboxx(2))

             do ii = imin,imax
                do jj = jmin,jmax
                   box_nr        = (jj-1_ip) * ni + ii
                   ll            = bin_struc(imesh) % pboel(box_nr) + nbono(box_nr)
                   nbono(box_nr) = nbono(box_nr) + 1
                   lboel(ll)     = ielem
                end do
             end do
          end do

       else if( ndime == 3 ) then

          do ielem = 1,nelem

             imin = int( ( xmima(1,1,ielem) - comin(1)) * dni , ip ) + 1
             imax = int( ( xmima(1,2,ielem) - comin(1)) * dni , ip ) + 1      

             jmin = int( ( xmima(2,1,ielem) - comin(2)) * dnj , ip ) + 1
             jmax = int( ( xmima(2,2,ielem) - comin(2)) * dnj , ip ) + 1      

             kmin = int( ( xmima(3,1,ielem) - comin(3)) * dnk , ip ) + 1
             kmax = int( ( xmima(3,2,ielem) - comin(3)) * dnk , ip ) + 1      

             imin = max(imin,1_ip)
             imax = min(imax,nboxx(1))
             jmin = max(jmin,1_ip)
             jmax = min(jmax,nboxx(2))
             kmin = max(kmin,1_ip)
             kmax = min(kmax,nboxx(3))

             do kk = kmin,kmax
                box_nr2 = ninj * (kk-1)
                do jj = jmin,jmax
                   box_nr1 = box_nr2 + nboxx(1) * (jj-1)
                   do ii = imin,imax
                      box_nr        = box_nr1 + ii
                      ll            = bin_struc(imesh) % pboel(box_nr) + nbono(box_nr)
                      nbono(box_nr) = nbono(box_nr) + 1
                      lboel(ll)     = ielem
                   end do
                end do
             end do
          end do

       end if

    end if
    !
    ! Deallocate memory
    !
    call elsest_cputim(time3)
    if( ipara(16) == 0 ) call memory_deallo(memor8,'BIN_STRUC % ELEMENT_BB','elsest_bin_preprocess',bin_struc(imesh) % element_bb)
    call memory_deallo(memor8,'NBONO','elsest_bin_preprocess',nbono)
    call elsest_cputim(time4)
    cputi(3)=time4-time3
    !
    ! Post-process bin
    !
!!!call elsest_binpos(ithre,ndime,mnode,npoin,nelem,lnods,ltype,coord,nnode,ipara)
    !
    ! Output statistics
    !
!!!if( ipara(7) /= 0 ) call elsest_statis(1_ip,imesh,ipara,ithre)

!!!call elsest_binpoi(imesh)

  end subroutine elsest_bin_preprocess

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    18/11/2015
  !> @brief   Domain bounding box
  !> @details Compute domain bounding box
  !>
  !-----------------------------------------------------------------------

  subroutine elsest_bounding_box(ndime,npoin,coord,comin,comax)
    integer(ip), intent(in)  :: ndime,npoin
    real(rp),    intent(in)  :: coord(ndime,npoin)
    real(rp),    intent(out) :: comin(ndime),comax(ndime)
    integer(ip)              :: idime

    do idime = 1,ndime
       comin(idime) = minval(coord(idime,:))
       comax(idime) = maxval(coord(idime,:))
    end do
    do idime = 1,ndime
       comax(idime) = comax(idime) + 0.000001_rp * ( comax(idime) - comin(idime) )
       comin(idime) = comin(idime) - 0.000001_rp * ( comax(idime) - comin(idime) )
    end do

  end subroutine elsest_bounding_box

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    18/11/2015
  !> @brief   Allocate memory for bin and/or quad/oc structures
  !>
  !-----------------------------------------------------------------------

  subroutine elsest_allocate(ipara,NUMBER_MESHES)

    integer(ip),           intent(in) :: ipara(*)
    integer(ip), optional, intent(in) :: NUMBER_MESHES
    integer(ip)                       :: nmesh
    !
    ! Number of meshes to allocate
    !
    if( present(NUMBER_MESHES) ) then
       nmesh = NUMBER_MESHES
    else
       nmesh  = 1_ip
    end if
    memor8 = 0

    if( ipara(8) == ELSEST_BIN_STRATEGY ) then
       !
       ! Allocate bin structure
       !
       if( associated(bin_struc) ) then
          call runend('ELSEST_ALLOCATE: BIN STRCTURE HAS ALREADY BEEN ALLOCATED')
       else
          allocate( bin_struc(nmesh) )
          bin_struc(1:nmesh) = bin_struc_init
       end if

    else if( ipara(8) == ELSEST_OCT_TREE_STRATEGY ) then
       !
       ! Allocate oct structure
       !
       if( associated(oct_struc) ) then
          call runend('ELSEST_ALLOCATE: OCT STRCTURE HAS ALREADY BEEN ALLOCATED')
       else
          allocate( oct_struc(nmesh) )
          oct_struc(1:nmesh) = oct_struc_init
       end if

    else if( ipara(8) == ELSEST_KD_TREE_STRATEGY ) then
       !
       ! Allocate kd-tree structure
       !
       !if( associated(tree_struc) ) then
       !   call runend('ELSEST_ALLOCATE: OCT STRCTURE HAS ALREADY BEEN ALLOCATED')
       !else
       !   allocate( tree_struc(nmesh) )
       !   tree_struc(1:nmesh) = tree_struc_init
       !end if

    end if

  end subroutine elsest_allocate

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    18/11/2015
  !> @brief   Find host element
  !> @details Find the host element of a test point using the bin structure
  !>
  !-----------------------------------------------------------------------

  subroutine elsest_bin_host_element(&
       imesh,ipara,rpara,mnode,ndime,npoin,nelem,lnnod,&
       lnods,ltype,coord,point_x,ifoun,shapt,derit,coloc,&
       dista,lchec)


    integer(ip), intent(in)           :: imesh
    integer(ip), intent(in)           :: ipara(*)
    real(rp),    intent(in)           :: rpara(*)
    integer(ip), intent(in)           :: mnode
    integer(ip), intent(in)           :: ndime
    integer(ip), intent(in)           :: npoin
    integer(ip), intent(in)           :: nelem
    integer(ip), intent(in)           :: lnnod(nelem)
    integer(ip), intent(in)           :: lnods(mnode,nelem)
    integer(ip), intent(in)           :: ltype(nelem)
    real(rp),    intent(in)           :: coord(ndime,npoin)
    real(rp),    intent(in)           :: point_x(*)
    integer(ip), intent(out)          :: ifoun
    real(rp),    intent(out)          :: shapt(*)
    real(rp),    intent(out)          :: derit(*)
    real(rp),    intent(out)          :: coloc(*)
    real(rp),    intent(out)          :: dista
    integer(ip), intent(in), optional :: lchec(*)
    integer(ip)                       :: curr_box_coor(3),box_nr
    integer(ip)                       :: array_size,ielem,ii,inode,ipoin,ilook,kk
    integer(ip)                       :: pnode,pelty,idime,kelem,ithre
    real(rp)                          :: toler
    real(rp)                          :: elcod(ndime,mnode)
    !
    ! Local pointers
    !
    integer(ip), pointer              :: nboxx(:)
    integer(ip), pointer              :: dataf
    integer(ip), pointer              :: lboel(:)
    integer(ip), pointer              :: pboel(:)
    type(i1p),   pointer              :: tboel(:)
    real(rp),    pointer              :: comin(:)
    real(rp),    pointer              :: comax(:)
    real(rp),    pointer              :: elebb(:,:,:)
    !
    ! If not allocated, create structure of mesh IMESH
    !
    ithre = 1
    if( bin_struc(imesh) % iallo == 0 ) then
       call runend('ELSEST_BIN_HOST_ELEMENT: MEMORY HAS NOT BEEN ALLOCATED')
    end if

    toler = abs(rpara(1))    
    ifoun = 0
    !
    ! Point to current bin structure
    !
    nboxx => bin_struc(imesh) % nboxx
    dataf => bin_struc(imesh) % dataf
    lboel => bin_struc(imesh) % lboel
    pboel => bin_struc(imesh) % pboel
    tboel => bin_struc(imesh) % tboel
    comin => bin_struc(imesh) % comin
    comax => bin_struc(imesh) % comax
    elebb => bin_struc(imesh) % element_bb
    !
    ! Check if point is outside the bounding box
    !
    do idime = 1,ndime
       if( point_x(idime) < comin(idime) .or. point_x(idime) > comax(idime) ) then
          return
       end if
    end do
    !
    ! Determine in which box (i,j,k) the point lies: curr_box_coor
    ! 
    call elsest_bin_box(ndime,nboxx,point_x,curr_box_coor,comin,comax)
    call elsest_bin_number(ndime,nboxx,curr_box_coor,box_nr)
    !
    ! Geometric CHECK
    !
    if( dataf == ELSEST_TYPE_BIN_STRUCTURE ) then
       array_size = memory_size(bin_struc(imesh) % tboel(box_nr) % l)
    else if( dataf == ELSEST_LINKED_LIST_BIN_STRUCTURE ) then
       array_size = bin_struc(imesh) % pboel(box_nr+1) - bin_struc(imesh) % pboel(box_nr)
       kk         = bin_struc(imesh) % pboel(box_nr) - 1
    end if
    !
    ! First try: Loop over elements in box
    !
    ii       = 0
    kelem    = 0
    dista    = huge(1.0_rp)

    if( ipara(16) == 1 ) then

       do while( ifoun == 0 .and. ii < array_size )
          ii = ii + 1
          if( dataf == ELSEST_TYPE_BIN_STRUCTURE ) then
             ielem = bin_struc(imesh) % tboel(box_nr) % l(ii)
          else if( dataf == ELSEST_LINKED_LIST_BIN_STRUCTURE ) then
             ielem = bin_struc(imesh) % lboel(ii+kk)
          end if

          ilook = 1
          if( ipara(14) /= 0 .and. present(lchec) ) then
             if( lchec(ielem) /= ipara(14) ) ilook = 0
          end if
          pelty = ltype(ielem)
          ifoun = 0

          if( ilook == 1 .and. pelty > 0 ) then

             if( elmgeo_inside_element_bounding_box(ndime,elebb(:,1,ielem),elebb(:,2,ielem),point_x,pelty) ) then

                pnode = lnnod(ielem)
                do inode = 1,pnode
                   ipoin= lnods(inode,ielem)
                   elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                end do

                call elmgeo_natural_coordinates(&
                     ndime,pelty,pnode,elcod,&
                     shapt,derit,point_x,coloc,ifoun,toler,BOUNDING_BOX=.false.)

             end if

             if( ifoun > 0 ) then

                ifoun = ielem
                dista = 0.0_rp

             end if
          end if
          call cputim(time2)
          time_total = time_total + (time2-time1)
       end do

    else
       do while( ifoun == 0 .and. ii < array_size )
          ii = ii + 1
          if( dataf == ELSEST_TYPE_BIN_STRUCTURE ) then
             ielem = bin_struc(imesh) % tboel(box_nr) % l(ii)
          else if( dataf == ELSEST_LINKED_LIST_BIN_STRUCTURE ) then
             ielem = bin_struc(imesh) % lboel(ii+kk)
          end if

          ilook = 1
          if( ipara(14) /= 0 .and. present(lchec) ) then
             if( lchec(ielem) /= ipara(14) ) ilook = 0
          end if
          pelty = ltype(ielem)

          if( ilook == 1 .and. pelty > 0 ) then
             pnode = lnnod(ielem)
             do inode = 1,pnode
                ipoin= lnods(inode,ielem)
                elcod(1:ndime,inode) = coord(1:ndime,ipoin)
             end do

             call elmgeo_natural_coordinates(&
                  ndime,pelty,pnode,elcod,&
                  shapt,derit,point_x,coloc,ifoun,toler)

             if( ifoun > 0 ) then

                ifoun = ielem
                dista = 0.0_rp

             end if
          end if

       end do
    end if

  end subroutine elsest_bin_host_element


  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    18/11/2015
  !> @brief   Find an element a toda costa
  !> @details Project
  !>
  !-----------------------------------------------------------------------

  subroutine elsest_host_element_a_toda_costa(&
       imesh,ipara,rpara,mnode,ndime,npoin,nelem,lnnod,&
       lnods,ltype,coord,point_x,ifoun,shapt,derit,coloc,&
       dista,lchec)

    integer(ip), intent(in)           :: imesh
    integer(ip), intent(in)           :: ipara(*)
    real(rp),    intent(in)           :: rpara(*)
    integer(ip), intent(in)           :: mnode
    integer(ip), intent(in)           :: ndime
    integer(ip), intent(in)           :: npoin
    integer(ip), intent(in)           :: nelem
    integer(ip), intent(in)           :: lnnod(nelem)
    integer(ip), intent(in)           :: lnods(mnode,nelem)
    integer(ip), intent(in)           :: ltype(nelem)
    real(rp),    intent(in)           :: coord(ndime,npoin)
    real(rp),    intent(in)           :: point_x(*)
    integer(ip), intent(out)          :: ifoun
    real(rp),    intent(out)          :: shapt(*)
    real(rp),    intent(out)          :: derit(*)
    real(rp),    intent(out)          :: coloc(*)
    real(rp),    intent(out)          :: dista
    integer(ip), intent(in), optional :: lchec(*)
    integer(ip)                       :: ielem,pnode,inode,pelty,jelem
    integer(ip)                       :: idime
    real(rp)                          :: cog(3),xdist,elcod(ndime,mnode)
!    real(rp)                          :: xx_intersection(3)
    !
    ! If an element is really required, look for an element si o si
    !
    dista = huge(1.0_rp)
    jelem = 0
    !
    ! Look for nearest element JELEM
    !                              
    if( ipara(14) /= 0 .and. present(lchec) ) then

       do ielem = 1,nelem
          if( lchec(ielem) == ipara(14) ) then
             if( ltype(ielem) > 0 ) then
                pnode = lnnod(ielem)
                do idime = 1,ndime
                   cog(idime) = sum(coord(idime,lnods(1:pnode,ielem))) / real(pnode,rp)
                end do
                xdist = dot_product(point_x(1:ndime)-cog(1:ndime),point_x(1:ndime)-cog(1:ndime))
                if( xdist <= dista ) then
                   dista = xdist
                   jelem = ielem
                end if
             end if
          end if
       end do

    else

       do ielem = 1,nelem
          if( ltype(ielem) > 0 ) then
             pnode = lnnod(ielem)
             do idime = 1,ndime
                cog(idime) = sum(coord(idime,lnods(1:pnode,ielem))) / real(pnode,rp)
             end do
             xdist = 0.0_rp
             xdist = dot_product(point_x(1:ndime)-cog(1:ndime),point_x(1:ndime)-cog(1:ndime))
             if( xdist <= dista ) then
                dista = xdist
                jelem = ielem
             end if
          end if
       end do

    end if
    !
    ! Compute intersection between cog and the faces and take minimum distance
    !
    if( jelem > 0 ) then
       dista = sqrt(dista)
       pelty = ltype(jelem)     
       pnode = lnnod(jelem)
       do idime = 1,ndime
          elcod(idime,1:pnode) = coord(idime,lnods(1:pnode,jelem))
       end do

       !do idime = 1,ndime
       !   cog(idime) = sum(elcod(idime,1:pnode)) / real(pnode,rp)
       !end do
       !call elmgeo_nearest_intersection_point_on_element_faces(&
       !     ndime,pelty,elcod,point_x,cog,xx_intersection,ifoun)
       !call elmgeo_natural_coordinates(          &
       !     ndime,pelty,pnode,elcod,shapt,derit, &
       !     xx_intersection,coloc,ifoun,abs(rpara(1)))

       !   write(*,'(a,3(1x,e13.6))') 'elsest_host_element_a_toda_costa, coord= ',xx_intersection(1:ndime)
       !   do inode = 1,pnode
       !      write(*,'(a,i7,60(1x,e13.6))') 'elsest_host_element_a_toda_costa, elcod= ',inode,elcod(1:ndime,inode)
       !   end do
       !
       ! Go for nearest node
       !       
       ifoun = jelem
       call elmgeo_nearest_element_node(             &
            ndime,pelty,pnode,elcod,point_x, &
            inode,coloc,shapt,derit)

    end if

  end subroutine elsest_host_element_a_toda_costa

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    18/11/2015
  !> @brief   Find the bin coordinates
  !> @details Find the bion coordinates i,j,k given the coordinates
  !>          of a test point
  !>
  !-----------------------------------------------------------------------

  subroutine elsest_bin_box(ndime,nboxx,point_x,curr_box_coor,comin,comax)

    integer(ip), intent(in)  :: ndime
    integer(ip), intent(in)  :: nboxx(ndime)
    real(rp),    intent(in)  :: point_x(ndime),comin(ndime),comax(ndime)
    integer(ip), intent(out) :: curr_box_coor(ndime)
    integer(ip)              :: idime

    do idime = 1,ndime
       curr_box_coor(idime) = int( ( (point_x(idime) - comin(idime)) / (comax(idime) - comin(idime)) ) * real(nboxx(idime),rp) , ip ) + 1    
       curr_box_coor(idime) = min( curr_box_coor(idime) , nboxx(idime) )
    end do

  end subroutine elsest_bin_box

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    18/11/2015
  !> @brief   Find the bin number
  !> @details Find the bin number as a function of coordinates i,j,k
  !>
  !-----------------------------------------------------------------------

  subroutine elsest_bin_number(ndime,nboxx,box_coord,box_nr)

    integer(ip), intent(in)  :: ndime
    integer(ip), intent(in)  :: nboxx(ndime),box_coord(ndime)
    integer(ip), intent(out) :: box_nr

    if( ndime == 2 ) then
       box_nr = (box_coord(2)-1)*nboxx(1) + box_coord(1)
    else
       box_nr = (box_coord(3)-1)*(nboxx(1)*nboxx(2)) + (box_coord(2)-1)*nboxx(1) + box_coord(1)
    end if

  end subroutine elsest_bin_number

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    18/11/2015
  !> @brief   Find host element
  !> @details Find the host element of a test point using the oct tree 
  !>          structure
  !>
  !-----------------------------------------------------------------------

  subroutine elsest_oct_host_element(&
       imesh,ipara,rpara,mnode,ndime,npoin,nelem,lnnod,&
       lnods,ltype,coord,point_x,ifoun,shapt,derit,coloc,&
       dista,lchec)

    integer(ip), intent(in)           :: imesh
    integer(ip), intent(in)           :: ipara(*)
    real(rp),    intent(in)           :: rpara(*)
    integer(ip), intent(in)           :: mnode
    integer(ip), intent(in)           :: ndime
    integer(ip), intent(in)           :: npoin
    integer(ip), intent(in)           :: nelem
    integer(ip), intent(in)           :: lnnod(nelem)
    integer(ip), intent(in)           :: lnods(mnode,nelem)
    integer(ip), intent(in)           :: ltype(nelem)
    real(rp),    intent(in)           :: coord(ndime,npoin)
    real(rp),    intent(in)           :: point_x(*)
    integer(ip), intent(out)          :: ifoun
    real(rp),    intent(out)          :: shapt(*)
    real(rp),    intent(out)          :: derit(*)
    real(rp),    intent(out)          :: coloc(*)
    real(rp),    intent(out)          :: dista
    integer(ip), intent(in), optional :: lchec(*)

    integer(ip)                       :: array_size,ielem,ii,inode,ipoin
    integer(ip)                       :: pnode,pelty,idime,ilook,kelem
    integer(ip)                       :: ichild
    real(rp)                          :: elcod(ndime,mnode)
    real(rp)                          :: toler
    real(rp),     pointer             :: comin(:)
    real(rp),     pointer             :: comax(:)
    real(rp),    pointer              :: elebb(:,:,:)
    type(octbox), pointer             :: current_o

    if( oct_struc(imesh) % iallo == 0 ) then
       call runend('ELSEST_OCTPRO: STRUCTURE HAS NOT BEEN ALLOCATED')
    end if
    !
    ! This is a new search
    !
    comin     => oct_struc(imesh) % comin
    comax     => oct_struc(imesh) % comax
    elebb     => oct_struc(imesh) % element_bb
    toler     =  abs(rpara(1))
    ifoun     =  0
    !
    ! Check if point is outside the bounding box
    !
    do idime = 1,ndime
       if( point_x(idime) < comin(idime) .or. point_x(idime) > comax(idime) ) then
          return
       end if
    end do
    !
    ! Find the quad/oct where the point lies
    !
    current_o  => oct_struc(imesh) % tree_root

    if( ndime == 3 ) then

       do while( current_o % whoiam == 0 )    
          childloop8: do ichild = 1,8           
             if(    point_x(1) >= current_o % children(ichild) % minc(1) .and. &
                  & point_x(1) <= current_o % children(ichild) % maxc(1) .and. &
                  & point_x(2) >= current_o % children(ichild) % minc(2) .and. &
                  & point_x(2) <= current_o % children(ichild) % maxc(2) .and. &
                  & point_x(3) >= current_o % children(ichild) % minc(3) .and. &
                  & point_x(3) <= current_o % children(ichild) % maxc(3) ) then
                current_o => current_o % children(ichild)
                exit childloop8
             end if
          end do childloop8
       end do

    else if( ndime == 2 ) then

       do while( current_o % whoiam == 0 )  
          childloop4: do ichild = 1,4
             if(    point_x(1) >= current_o % children(ichild) % minc(1) .and. &
                  & point_x(1) <= current_o % children(ichild) % maxc(1) .and. &
                  & point_x(2) >= current_o % children(ichild) % minc(2) .and. &
                  & point_x(2) <= current_o % children(ichild) % maxc(2)  ) then
                current_o => current_o % children(ichild)
                exit childloop4
             end if
          end do childloop4
       end do

    end if
    ! 
    ! Perform search over elements inside current box
    !
    array_size = current_o % nelembox
    ii         = 0
    kelem      = 0
    dista      = huge(1.0_rp)

    do while( ifoun == 0 .and. ii < array_size )
       ii    = ii+1
       ielem = current_o % elems(ii)

       ilook = 1
       if( ipara(14) /= 0 .and. present(lchec) ) then
          if( lchec(ielem) /= ipara(14) ) ilook = 0
       end if

       if( ilook == 1 ) then
          pelty = ltype(ielem)
          if( pelty > 0 ) then
             pnode = lnnod(ielem)
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                elcod(1:ndime,inode) = coord(1:ndime,ipoin)
             end do

             if( ipara(16) == 1 ) then
                if( elmgeo_inside_element_bounding_box(ndime,elebb(:,1,ielem),elebb(:,2,ielem),point_x,pelty) ) then
                   call elmgeo_natural_coordinates(          &
                        ndime,pelty,pnode,elcod,shapt,derit, &
                        point_x,coloc,ifoun,toler,BOUNDING_BOX=.false.)                   
                end if
             else
                call elmgeo_natural_coordinates(          &
                     ndime,pelty,pnode,elcod,shapt,derit, &
                     point_x,coloc,ifoun,toler)
             end if

             if( ifoun > 0 ) then
                ifoun = ielem
                dista = 0.0_rp
             end if

          end if
       end if
    end do

  end subroutine elsest_oct_host_element

  function elsest_intost(integ)
    !-------------------------------------
    !
    !  Convert an integer(ip) to a string
    !
    !-------------------------------------
    implicit none
    integer(ip)   :: integ
    character(20) :: elsest_intost
    character(20) :: intaux

    write(intaux,*) integ
    elsest_intost=adjustl(intaux)

  end function elsest_intost

  subroutine elsest_oct_preprocess(&
       ipara,rpara,imesh,mnode,ndime,npoin,nelem,lnnod,&
       lnods,ltype,coord)

    integer(ip),  intent(in)  :: ipara(*)
    real(rp),     intent(in)  :: rpara(*)
    integer(ip),  intent(in)  :: imesh
    integer(ip),  intent(in)  :: mnode
    integer(ip),  intent(in)  :: ndime
    integer(ip),  intent(in)  :: npoin
    integer(ip),  intent(in)  :: nelem
    integer(ip),  intent(in)  :: lnnod(nelem)
    integer(ip),  intent(in)  :: lnods(mnode,nelem)
    integer(ip),  intent(in)  :: ltype(nelem)
    real(rp),     intent(in)  :: coord(ndime,npoin)

    integer(ip)               :: ipoin
    integer(ip)               :: i,kpoin,ielem,counter
    integer(ip)               :: divmax
    real(rp)                  :: toler,xsize
    real(rp)                  :: time1,time2,time4,time5      
    logical(lg)               :: conti

    real(rp),     pointer     :: comin(:)
    real(rp),     pointer     :: comax(:)
    type(octbox), pointer     :: old_pointer
    type(octbox), pointer     :: tm1_pointer
    type(octbox), pointer     :: tm2_pointer
    type(octbox), pointer     :: tree_root
    real(rp),     pointer     :: xmima(:,:,:)
    integer(ip),  pointer     :: lboel_oct(:)
    integer(ip)               :: limit
    integer(ip),  pointer     :: kstat(:)
    real(rp),     pointer     :: cputi(:)

    type(octbox), pointer     :: current_o
    !
    ! Nullify pointers
    !
    nullify(comin)
    nullify(comax)
    nullify(old_pointer)
    nullify(tm1_pointer)
    nullify(tm2_pointer)
    nullify(tree_root)
    nullify(lboel_oct)
    nullify(xmima)

    call elsest_cputim(time1)
    oct_struc(imesh) % iallo  = 1 
    oct_struc(imesh) % divmax = 2**ndime
    !
    ! Allocate memory
    !
    call memory_alloca(memor8,'OCT_STRUC % CPUTI'     ,'elsest_oct_preprocess',oct_struc(imesh) % cputi,10_ip)
    call memory_alloca(memor8,'OCT_STRUC % KSTAT'     ,'elsest_oct_preprocess',oct_struc(imesh) % kstat,10_ip)
    call memory_alloca(memor8,'OCT_STRUC % ELEMENT_BB','elsest_oct_preprocess',oct_struc(imesh) % element_bb,ndime,2_ip,nelem)
    !
    ! Point to current mesh (IMESH) structure
    !
    tree_root => oct_struc(imesh) % tree_root
    kstat     => oct_struc(imesh) % kstat
    comin     => oct_struc(imesh) % comin
    comax     => oct_struc(imesh) % comax
    cputi     => oct_struc(imesh) % cputi
    !
    ! Initialize parameters
    !
    cputi    = 0.0_rp
    kstat(:) = 0
    kstat(1) = huge(1_ip)
    kstat(3) = huge(1_ip)
    limit    = ipara(9)
    divmax   = 2**ndime
    toler    = 1.0e-6_rp
    !
    ! Compute bounding box
    !
    call elsest_bounding_box(ndime,npoin,coord,comin,comax)
    comin = comin - zeror
    comax = comax + zeror
    !
    ! Allocate memory for tree root
    ! 
    allocate( oct_struc(imesh) % tree_root )
    oct_struc(imesh) % tree_root = octbox_init

    !call memory_alloca(memor8,'OCT_STRUC % TREE_ROOT % NODES','elsest_oct_preprocess',oct_struc(imesh) % tree_root % nodes,npoin)
    call memory_alloca(memor8,'CURRENT_O % NODES','elsest_oct_preprocess',oct_struc(imesh) % tree_root % nodes,npoin)

    !allocate(oct_struc(imesh) % tree_root % nodes(npoin),stat=istat)
    !call elsest_memchk(0_ip,ithre,istat,memor(2),'TREE_ROOT % NODES','elsest_octpre',oct_struc(imesh) % tree_root % nodes)

    !--------------------------------------------------------------------
    !
    ! Init tree_root values
    !
    !--------------------------------------------------------------------

    tree_root            => oct_struc(imesh) % tree_root
    current_o            => tree_root
    current_o % npoinbox =  0
    current_o % id       =  0
    current_o % level    =  0
    current_o % whoiam   =  0
    current_o % childid  =  0
    current_o % nelembox =  0
    current_o % npoinbox =  npoin

    do ipoin = 1,npoin
       current_o % nodes(ipoin) = ipoin
    end do
    current_o % minc(1:ndime) = comin(1:ndime)
    current_o % maxc(1:ndime) = comax(1:ndime)

    nullify(current_o % parent)
    call elsest_cputim(time2)
    cputi(1) = time2-time1 

    !--------------------------------------------------------------------
    !
    ! Generation of quad/oct tree
    !
    !--------------------------------------------------------------------
    !
    ! Compute element bounding boxes
    !
    !call memory_alloca(memor8,'XMIMA','elsest_oct_preprocess',xmima,3_ip,2_ip,nelem)
    !call elsest_element_bounding_box(mnode,ndime,npoin,nelem,lnnod,lnods,coord,xmima)
    call elsest_element_bounding_boxes(&
         mnode,ndime,npoin,nelem,lnnod,lnods,coord,oct_struc(imesh) % element_bb)
    xmima => oct_struc(imesh) % element_bb
    !
    ! Allocate memory
    !
    call memory_alloca(memor8,'LBOEL_OCT','elsest_oct_preprocess',lboel_oct,nelem)

    divmax  = 2**ndime
    conti   = .true.
    counter = 0

    do while( conti )
       !
       ! If maximum number of points inside current box is exceeded, subdivide
       !     
       if( current_o % npoinbox > limit ) then
          call elsest_cputim(time1)
          allocate( current_o % children(divmax) )
          current_o % children(1:divmax) = octbox_init
          !
          ! Give birth to my DIVMAX children
          !
          do i = 1,divmax
             counter                            =  counter+1
             current_o % children(i) % id       =  counter
             current_o % children(i) % childid  =  i
             current_o % children(i) % level    =  current_o % level + 1 
             current_o % children(i) % whoiam   =  0  
             current_o % children(i) % npoinbox =  0
             current_o % children(i) % nelembox =  0
             current_o % children(i) % parent   => current_o

             call memory_alloca(memor8,'CURRENT_O % NODES','elsest_oct_preprocess',current_o % children(i) % nodes,current_o % npoinbox)
          end do
          !
          ! Compute the coordinates of my children
          !
          do i = 1,ndime
             current_o % children(1) % minc(i) = current_o % minc(i)
             current_o % children(1) % maxc(i) = (current_o % maxc(i) + current_o % minc(i))*0.5_rp
          end do
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

          if( ndime == 3 ) then
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
          !
          ! Offer my nodes to my children
          !
          if( ndime == 2 ) then
             do i = 1,4
                do ipoin = 1,current_o % npoinbox
                   kpoin = current_o % nodes(ipoin)
                   if(    coord(1,kpoin) >= current_o % children(i) % minc(1) .and. &
                        & coord(2,kpoin) >= current_o % children(i) % minc(2) .and. &
                        & coord(1,kpoin) <  current_o % children(i) % maxc(1) .and. &
                        & coord(2,kpoin) <  current_o % children(i) % maxc(2) ) then
                      current_o % children(i) % npoinbox = current_o % children(i) % npoinbox + 1
                      current_o % children(i) % nodes(current_o % children(i) % npoinbox) = kpoin
                   end if
                end do
             end do
          else
             do i = 1,8
                do ipoin = 1,current_o % npoinbox
                   kpoin = current_o % nodes(ipoin)
                   !if(    coord(1,kpoin) >= current_o % children(i) % minc(1) .and. &
                   !     & coord(2,kpoin) >= current_o % children(i) % minc(2) .and. &
                   !     & coord(3,kpoin) >= current_o % children(i) % minc(3) .and. &
                   !     & coord(1,kpoin) <  current_o % children(i) % maxc(1) .and. &
                   !     & coord(2,kpoin) <  current_o % children(i) % maxc(2) .and. &
                   !     & coord(3,kpoin) <  current_o % children(i) % maxc(3) ) then
                   if(                coord(1,kpoin) >= current_o % children(i) % minc(1) ) then
                      if(             coord(2,kpoin) >= current_o % children(i) % minc(2) ) then
                         if(          coord(3,kpoin) >= current_o % children(i) % minc(3) ) then 
                            if(       coord(1,kpoin) <  current_o % children(i) % maxc(1) ) then 
                               if(    coord(2,kpoin) <  current_o % children(i) % maxc(2) ) then 
                                  if( coord(3,kpoin) <  current_o % children(i) % maxc(3) ) then
                                     current_o % children(i) % npoinbox = current_o % children(i) % npoinbox + 1
                                     current_o % children(i) % nodes(current_o % children(i) % npoinbox) = kpoin
                                  end if
                               end if
                            end if
                         end if
                      end if
                   end if
                end do
             end do
          end if

          call memory_deallo(memor8,'CURRENT_O % NODES','elsest_oct_preprocess',current_o % nodes)

          current_o % whoiam   =  0
          current_o % npoinbox =  0
          current_o            => current_o % children(1)
          call elsest_cputim(time2)
          cputi(2) = cputi(2) + time2-time1

       else if(current_o % id == 0 .and. current_o % npoinbox <= limit ) then
          !
          ! If the Padrino has too few elements
          !
          call memory_alloca(memor8,'CURRENT_O % ELEMS','elsest_oct_preprocess',current_o % elems,nelem)

          do ielem = 1,nelem
             current_o % elems(ielem) = ielem
          end do
          conti = .false.
          current_o % nelembox =  nelem
          current_o % whoiam   =  1
          current_o            => old_pointer

       else 
          !
          ! if limit of points inside box is not exceeded, assign elements
          !
          call elsest_cputim(time1)
          kstat(6) = kstat(6) + 1
          kstat(8) = kstat(8) + current_o % npoinbox
          if( current_o % npoinbox < kstat(1) ) kstat(1) = current_o % npoinbox
          if( current_o % npoinbox > kstat(2) ) kstat(2) = current_o % npoinbox 

          xsize = maxval(current_o % maxc(1:ndime)-current_o % minc(1:ndime))*toler

          if( ndime == 1 ) then
             do ielem = 1,nelem
                if( xmima(1,1,ielem) <= current_o % maxc(1) + xsize ) then
                   if( xmima(1,2,ielem) >= current_o % minc(1) - xsize ) then
                      current_o % nelembox = current_o % nelembox + 1
                      lboel_oct(current_o % nelembox) = ielem
                   end if
                end if
             end do
          else if( ndime == 2 ) then
             do ielem = 1,nelem
                if( xmima(1,1,ielem) <= current_o % maxc(1) + xsize ) then
                   if( xmima(1,2,ielem) >= current_o % minc(1) - xsize ) then
                      if( xmima(2,1,ielem) <= current_o % maxc(2) + xsize ) then
                         if( xmima(2,2,ielem) >= current_o % minc(2) - xsize ) then
                            current_o % nelembox = current_o % nelembox + 1
                            lboel_oct(current_o % nelembox) = ielem
                         end if
                      end if
                   end if
                end if
             end do
          else if( ndime == 3 ) then
             do ielem = 1,nelem
                if( xmima(1,1,ielem) <= current_o % maxc(1) + xsize ) then
                   if( xmima(1,2,ielem) >= current_o % minc(1) - xsize ) then
                      if( xmima(2,1,ielem) <= current_o % maxc(2) + xsize ) then
                         if( xmima(2,2,ielem) >= current_o % minc(2) - xsize  ) then
                            if( xmima(3,1,ielem) <= current_o % maxc(3) + xsize ) then
                               if( xmima(3,2,ielem) >= current_o % minc(3) - xsize ) then
                                  current_o % nelembox = current_o % nelembox + 1
                                  lboel_oct(current_o % nelembox) = ielem
                               end if
                            end if
                         end if
                      end if
                   end if
                end if
             end do
          end if
          !
          ! Look for elements crossing the Quad
          !
          if(current_o % nelembox < kstat(3)) kstat(3) = current_o % nelembox
          if(current_o % nelembox > kstat(4)) kstat(4) = current_o % nelembox
          kstat(7) = kstat(7)+current_o % nelembox

          call memory_deallo(memor8,'CURRENT_O % NODES','elsest_oct_preprocess',current_o % nodes)
          !call elsest_memchk(2_ip,ithre,istat,memor(2),'CURRENT_O % NODES','elsest_octsub',current_o % nodes)
          !deallocate(current_o % nodes,stat=istat)
          !if(istat/=0) call elsest_memerr(2_ip,'CURRENT%NODES','elsest_octsub',0_ip)     
          !
          ! Here we assign elements
          !
          if( current_o % nelembox /= 0 ) then
             current_o % whoiam  = 1

             call memory_alloca(memor8,'CURRENT_O % ELEMS','elsest_oct_preprocess',current_o % elems,current_o % nelembox)
             !allocate( current_o % elems(current_o % nelembox),stat=istat)
             !call elsest_memchk(0_ip,ithre,istat,memor(2),'CURRENT%ELEMS','elsest_octsub',current_o % elems)
             do ielem = 1,current_o % nelembox
                current_o %elems(ielem) = lboel_oct(ielem)
             end do
          else
             current_o % whoiam  = 2
          end if
          call elsest_cputim(time2)
          cputi(3)=cputi(3)+time2-time1

          if( current_o % childid < divmax .and. current_o % id /= 0 ) then
             !
             ! Go to next children
             !
             tm1_pointer      => current_o
             tm2_pointer      => tm1_pointer%parent%children(tm1_pointer%childid+1)
             current_o => tm2_pointer
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
                if(current_o % parent%id == 0) then
                   conti=.false.
                   exit noparent
                else
                   if(current_o % parent%childid /=divmax) then 
                      tm1_pointer       => current_o
                      tm2_pointer       => tm1_pointer%parent%parent%children(tm1_pointer%parent%childid+1)
                      current_o  => tm2_pointer
                      exit
                   else 
                      current_o => current_o % parent
                   end if
                end if
             end do noparent

          else 
             !
             ! Wrong child ID
             !
             call elsest_runend('WRONG CHILD ID: '//trim(elsest_intost(current_o % childid)))    
          end if

       end if

10     continue
       old_pointer => current_o

    end do
    !
    ! Deallocate memory
    !
    call memory_deallo(memor8,'LBOEL_OCT','elsest_bin_preprocess',lboel_oct)
    if( ipara(16) == 0 ) call memory_deallo(memor8,'OCT_STRUC % ELEMENT_BB','elsest_oct_preprocess',oct_struc(imesh) % element_bb)
    !call memory_deallo(memor8,'XMIMA'    ,'elsest_bin_preprocess',xmima)
    !
    ! Postprocess of quad/oct tree
    !
    call elsest_oct_postprocess(imesh,ndime,mnode,npoin,nelem,lnods,ltype,coord,lnnod,ipara)

    call elsest_cputim(time4)
    call elsest_cputim(time5) 
    cputi(4) = time5-time4

    !if(ipara(7)/=0) call elsest_statis(1_ip,imesh,ipara,ithre)

  end subroutine elsest_oct_preprocess

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    18/11/2015
  !> @brief   Post-process the Oct/Quad-tree mesh
  !>
  !-----------------------------------------------------------------------

  subroutine elsest_oct_postprocess(&
       imesh,ndime,mnode,npoin,nelem,lnods,ltype,coord,lnnod,ipara)

    integer(ip), intent(in)  :: imesh,ndime,mnode,npoin,nelem
    integer(ip), intent(in)  :: lnods(mnode,nelem),ltype(nelem),lnnod(*)
    integer(ip), intent(in)  :: ipara(*)
    real(rp),    intent(in)  :: coord(ndime,npoin)
    integer(ip)              :: ipoin,ielem,idime,kpoin
    type(octbox), pointer    :: current_o

    iunit(1) = ipara( 7)
    iunit(2) = ipara(12)
    iunit(3) = ipara(13)

    if( iunit(2) > 0 ) then
       !
       ! Coordinates
       !
       current_o =>  oct_struc(imesh) % tree_root
       ipoin = 0
       ielem = 0
       if( ndime == 2 ) then
          write(iunit(2),*) 'MESH  ELSEST_QUAD dimension 2 Elemtype Quadrilateral Nnode 4'
       else
          write(iunit(2),*) 'MESH  ELSEST_OCT  dimension 3 Elemtype Hexahedra Nnode 8'
       end if
       write(iunit(2),*) 'coordinates'
       call elsest_oct_recursive(1_ip,ndime,ipoin,ielem,current_o)
       do kpoin = 1,npoin
          write(iunit(2),*) ipoin+kpoin,(coord(idime,kpoin),idime=1,ndime)
       end do
       write(iunit(2),*)  'end coordinates'

       current_o => oct_struc(imesh) % tree_root
       ipoin = 0
       ielem = 0
       write(iunit(2),*) 'elements'
       call elsest_oct_recursive(2_ip,ndime,ipoin,ielem,current_o)
       write(iunit(2),*) 'end elements'
       call elsest_gid_format(ndime,mnode,nelem,ielem,ipoin,lnnod,lnods,ltype)
    end if

    if( iunit(3) > 0 ) then

       write(iunit(3),*) 'GiD Post Results File 1.0'
       write(iunit(3),*) 'GaussPoints GP_QUAD4 Elemtype Quadrilateral'
       write(iunit(3),*) 'Number of Gauss Points: 1'
       write(iunit(3),*) 'Natural Coordinates: Internal'
       write(iunit(3),*) 'End GaussPoints'

       current_o => oct_struc(imesh) % tree_root
       ielem     =  0
       write(iunit(3),*) 'Result NODE_NUMBER ELSEST 0 Scalar OnGaussPoints GP_QUAD4'
       write(iunit(3),*) 'ComponentNames NODE_NUMBER'
       write(iunit(3),*) 'values'
       call elsest_oct_recursive(3_ip,ndime,ipoin,ielem,current_o)
       write(iunit(3),*) 'end values'

       current_o => oct_struc(imesh) % tree_root
       ielem     =  0
       write(iunit(3),*) 'Result ELEMENT_NUMBER ELSEST 0 Scalar OnGaussPoints GP_QUAD4'
       write(iunit(3),*) 'ComponentNames ELEMENT_NUMBER'
       write(iunit(3),*) 'values'
       call elsest_oct_recursive(4_ip,ndime,ipoin,ielem,current_o)
       write(iunit(3),*) 'end values'

       current_o => oct_struc(imesh) % tree_root
       ielem     =  0
       write(iunit(3),*) 'Result LEVEL ELSEST 0 Scalar OnGaussPoints GP_QUAD4'
       write(iunit(3),*) 'ComponentNames LEVEL'
       write(iunit(3),*) 'values'
       call elsest_oct_recursive(5_ip,ndime,ipoin,ielem,current_o)
       write(iunit(3),*) 'end values'

    end if

  end subroutine elsest_oct_postprocess

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    18/11/2015
  !> @brief   GiD format
  !>
  !-----------------------------------------------------------------------

  subroutine elsest_gid_format(ndime,mnode,nelem,kelem,kpoin,lnnod,lnods,ltype)

    integer(ip), intent(in)  :: ndime,mnode,nelem,kelem,kpoin
    integer(ip), intent(in)  :: lnnod(*),lnods(mnode,nelem),ltype(nelem)
    integer(ip)              :: ielem,iesta_dom,iesto_dom
    integer(ip)              :: ielty,inode
    integer(ip), allocatable :: lexis(:)

    allocate(lexis(nelem))
    do ielem=1,nelem
       lexis(abs(ltype(ielem)))=1
    end do

    if(ndime==1) then
       iesta_dom=BAR02
       iesto_dom=BAR03
    else if(ndime==2) then
       iesta_dom=TRI03
       iesto_dom=QUA09
    else
       iesta_dom=TET04
       iesto_dom=HEX27
    end if

    do ielty=iesta_dom,iesto_dom
       if(lexis(ielty)/=0) then
          !
          ! Header
          !
          if(ielty<10) then
             write(iunit(2),10)&
                  'ELSEST_BACKGROUND',ielty,max(2_ip,ndime),&
                  adjustl(trim(element_type(ielty) % nametopo)),element_type(ielty) % number_nodes
          else
             write(iunit(2),11)&
                  'ELSEST_BACKGROUND',ielty,max(2_ip,ndime),&
                  adjustl(trim(element_type(ielty) % nametopo)),element_type(ielty) % number_nodes
          end if
          !
          ! Connectivity
          !
          write(iunit(2),*) 'elements'
          do ielem=1,nelem
             if(abs(ltype(ielem))==ielty) then
                write(iunit(2),4) ielem+kelem,&
                     (lnods(inode,ielem)+kpoin,inode=1,lnnod(ielem)),0_ip
             end if
          end do
          write(iunit(2),*) 'end elements'
       end if
    end do
    deallocate(lexis)

10  format('MESH ',a,i1,' dimension ',i1,' Elemtype ',a,' Nnode ',i2)
11  format('MESH ',a,i2,' dimension ',i1,' Elemtype ',a,' Nnode ',i2)
2   format(a)
3   format(i7, 3(1x,e16.8e3))
4   format(i7,50(1x,i7))

  end subroutine elsest_gid_format

  subroutine elsest_oct_recursive(itask,ndime,ipoin,ielem,current_o)
    !
    ! Descend down the hierarchy for post-process
    !
    integer(ip), intent(in)    :: itask,ndime
    integer(ip), intent(inout) :: ipoin,ielem
    integer(ip)                :: inode,lnods(8),divmax
    type(octbox), pointer      :: current_o

    divmax  = 2**ndime

    do 
       !
       ! First go to deepest level in first branch
       !
       do while( current_o % whoiam == 0 )
          current_o => current_o % children(1)
       end do
       !
       ! Current bin has elements
       !
       if(current_o % whoiam /= 0 ) then
          if(itask==1) then
             if(ndime==2) then
                ipoin = ipoin+1
                write(iunit(2),*) ipoin,current_o % minc(1),current_o % minc(2)
                ipoin = ipoin+1
                write(iunit(2),*) ipoin,current_o % maxc(1),current_o % minc(2)
                ipoin = ipoin+1
                write(iunit(2),*) ipoin,current_o % maxc(1),current_o % maxc(2)
                ipoin = ipoin+1
                write(iunit(2),*) ipoin,current_o % minc(1),current_o % maxc(2)
             end if
          else if(itask==2) then
             if(ndime==2) then
                inode = 0
                ielem = ielem+1
                ipoin = ipoin+1
                inode = inode+1
                lnods(inode)=ipoin
                ipoin = ipoin+1
                inode = inode+1
                lnods(inode)=ipoin
                ipoin = ipoin+1
                inode = inode+1
                lnods(inode)=ipoin
                ipoin = ipoin+1
                inode = inode+1
                lnods(inode)=ipoin
                write(iunit(2),'(6(1x,i9))')  ielem,lnods(1),lnods(2),lnods(3),lnods(4),current_o % level
             end if
          else if(itask==3) then
             ielem = ielem+1
             write(iunit(3),*) ielem,current_o % npoinbox
          else if(itask==4) then
             ielem = ielem+1
             write(iunit(3),*) ielem,current_o % nelembox   
          else if(itask==5) then
             ielem = ielem+1
             write(iunit(3),*) ielem,current_o % level
          end if
       end if

       if(current_o % childid < divmax .and. current_o % childid /=0) then
          !
          ! I'm not the last child neither the Padrino
          !
          current_o => current_o % parent%children(current_o % childid+1)

       else if(current_o % childid==divmax) then
          !
          ! I'm the last child of this generation: postprocess 
          !
          do while(current_o % id > 0 )
             if(current_o % parent % id == 0) then
                return
             else
                if(current_o % parent%childid /=divmax) then
                   current_o => current_o % parent%parent%children(current_o % parent%childid+1)
                   exit
                else 
                   current_o => current_o % parent
                end if
             end if
          end do

       else if(current_o % id==0) then
          !
          ! I'm the Padrino
          !
          return
       end if

    end do

  end subroutine elsest_oct_recursive

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    18/11/2015
  !> @brief   Deallocate a bin structure
  !>
  !-----------------------------------------------------------------------

  subroutine elsest_bin_deallocate(imesh)
    integer(ip), intent(in) :: imesh

    if( bin_struc(imesh) % iallo /= 0 ) then

       bin_struc(imesh) % iallo = 0
       call memory_deallo(memor8,'BIN_STRUC % LBOEL'     ,'elsest_bin_deallocate',bin_struc(imesh) % lboel)
       call memory_deallo(memor8,'BIN_STRUC % PBOEL'     ,'elsest_bin_deallocate',bin_struc(imesh) % pboel)
       call memory_deallo(memor8,'BIN_STRUC % TBOEL'     ,'elsest_bin_deallocate',bin_struc(imesh) % tboel)
       call memory_deallo(memor8,'BIN_STRUC % KSTAT'     ,'elsest_bin_deallocate',bin_struc(imesh) % kstat)
       call memory_deallo(memor8,'BIN_STRUC % CPUTI'     ,'elsest_bin_deallocate',bin_struc(imesh) % cputi)
       call memory_deallo(memor8,'BIN_STRUC % ELEMENT_BB','elsest_bin_deallocate',bin_struc(imesh) % element_bb)

    end if

  end subroutine elsest_bin_deallocate

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    18/11/2015
  !> @brief   Deallocate a bin structure
  !>
  !-----------------------------------------------------------------------

  subroutine elsest_tree_deallocate(imesh)
    integer(ip), intent(in) :: imesh

  end subroutine elsest_tree_deallocate

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    18/11/2015
  !> @brief   Deallocate an oct-tree structure
  !>
  !-----------------------------------------------------------------------

  subroutine elsest_oct_deallocate(imesh)
    implicit none
    integer(ip),  intent(in) :: imesh
    type(octbox), pointer    :: current_o
    integer(ip)              :: divmax
    logical(lg)              :: conti

    if( oct_struc(imesh) % iallo /= 0 ) then

       oct_struc(imesh) % iallo =  0
       divmax                   =  oct_struc(imesh) % divmax
       current_o                => oct_struc(imesh) % tree_root

       conti=.true.
       do while(conti)
          !
          ! First go to deepest level in first branch
          !
          do while( current_o % whoiam == 0 )
             current_o => current_o % children(1)
          end do
          !
          ! Deallocate list of elements
          !
          if( current_o % whoiam == 1 ) then
             call memory_deallo(memor8,'CURRENT_O % ELEMS','elsest_oct_deallocate',current_o % elems)
             call memory_deallo(memor8,'CURRENT_O % NODES','elsest_oct_deallocate',current_o % nodes)
          end if

          if( current_o % childid < divmax .and. current_o % childid /= 0 ) then
             !
             ! I'm not the last child neither the Padrino
             !
             current_o => current_o % parent % children(current_o % childid+1)

          else if( current_o % childid == divmax ) then
             !
             ! I'm the last child
             !
             current_o => current_o % parent 
             deallocate( current_o % children)
             current_o % whoiam = 3

          else if( current_o % id == 0 ) then
             !
             ! I'm the Padrino: end of deallocation             
             deallocate( current_o)
             conti = .false.

          end if

       end do

       call memory_deallo(memor8,'OCT_STRUC % CPUTI'     ,'elsest_oct_deallocate',oct_struc(imesh) % cputi)
       call memory_deallo(memor8,'OCT_STRUC % KSTAT'     ,'elsest_oct_deallocate',oct_struc(imesh) % kstat)
       call memory_deallo(memor8,'OCT_STRUC % ELEMENT_BB','elsest_oct_deallocate',oct_struc(imesh) % element_bb)

    end if

  end subroutine elsest_oct_deallocate

  subroutine elsest_kdtree_preprocess(&
       ipara,rpara,imesh,mnode,ndime,npoin,nelem,lnnod,&
       lnods,ltype,coord)

    integer(ip), intent(in)  :: ipara(*)
    real(rp),    intent(in)  :: rpara(*)
    integer(ip), intent(in)  :: imesh
    integer(ip), intent(in)  :: mnode
    integer(ip), intent(in)  :: ndime
    integer(ip), intent(in)  :: npoin
    integer(ip), intent(in)  :: nelem
    integer(ip), intent(in)  :: lnnod(nelem)
    integer(ip), intent(in)  :: lnods(mnode,nelem)
    integer(ip), intent(in)  :: ltype(nelem)
    real(rp),    intent(in)  :: coord(ndime,npoin)

  end subroutine elsest_kdtree_preprocess

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    18/11/2015
  !> @brief   Find host element
  !> @details Find the host element of a test point using the bin structure
  !>
  !-----------------------------------------------------------------------

  subroutine elsest_kd_host_element(&
       imesh,ipara,rpara,mnode,ndime,npoin,nelem,lnnod,&
       lnods,ltype,coord,point_x,ifoun,shapt,derit,coloc,&
       dista,lchec)

    integer(ip), intent(in)           :: imesh
    integer(ip), intent(in)           :: ipara(*)
    real(rp),    intent(in)           :: rpara(*)
    integer(ip), intent(in)           :: mnode
    integer(ip), intent(in)           :: ndime
    integer(ip), intent(in)           :: npoin
    integer(ip), intent(in)           :: nelem
    integer(ip), intent(in)           :: lnnod(nelem)
    integer(ip), intent(in)           :: lnods(mnode,nelem)
    integer(ip), intent(in)           :: ltype(nelem)
    real(rp),    intent(in)           :: coord(ndime,npoin)
    real(rp),    intent(in)           :: point_x(ndime)
    integer(ip), intent(out)          :: ifoun
    real(rp),    intent(out)          :: shapt(*)
    real(rp),    intent(out)          :: derit(*)
    real(rp),    intent(out)          :: coloc(*)
    real(rp),    intent(out)          :: dista
    integer(ip), intent(in), optional :: lchec(*)
 
  end subroutine elsest_kd_host_element

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    15/06/2018
  !> @brief   Element bounding boxes
  !> @details ELEMENT_BB(1,1:ndime,IELEM) = min
  !>          ELEMENT_BB(2,1:ndime,IELEM) = max
  !>
  !-----------------------------------------------------------------------

  subroutine elsest_element_bounding_boxes(&
       mnode,ndime,npoin,nelem,lnnod,lnods,coord,element_bb)

    integer(ip), intent(in)  :: mnode
    integer(ip), intent(in)  :: ndime
    integer(ip), intent(in)  :: npoin
    integer(ip), intent(in)  :: nelem
    integer(ip), intent(in)  :: lnnod(nelem)
    integer(ip), intent(in)  :: lnods(mnode,nelem)
    real(rp),    intent(in)  :: coord(ndime,npoin)
    real(rp),    intent(out) :: element_bb(ndime,2,nelem)
    real(rp)                 :: elcod(ndime,mnode)
    integer(ip)              :: ielem,inode,ipoin,pnode

    if(      ndime == 1 ) then
       do ielem = 1,nelem
          pnode = lnnod(ielem)
          do inode = 1,pnode
             ipoin = lnods(inode,ielem)
             elcod(1:ndime,inode) = coord(1:ndime,ipoin)
          end do
          element_bb(1,1,ielem) = minval(elcod(1,1:pnode)) - zeror
          element_bb(1,2,ielem) = maxval(elcod(1,1:pnode)) + zeror
       end do
    else if( ndime == 2 ) then
       do ielem = 1,nelem
          pnode = lnnod(ielem)
          do inode = 1,pnode
             ipoin = lnods(inode,ielem)
             elcod(1:ndime,inode) = coord(1:ndime,ipoin)
          end do
          element_bb(1,1,ielem) = minval(elcod(1,1:pnode)) - zeror
          element_bb(2,1,ielem) = minval(elcod(2,1:pnode)) - zeror
          element_bb(1,2,ielem) = maxval(elcod(1,1:pnode)) + zeror
          element_bb(2,2,ielem) = maxval(elcod(2,1:pnode)) + zeror
       end do
    else
       do ielem = 1,nelem
          pnode = lnnod(ielem)
          do inode = 1,pnode
             ipoin = lnods(inode,ielem)
             elcod(1:ndime,inode) = coord(1:ndime,ipoin)
          end do
          element_bb(1,1,ielem) = minval(elcod(1,1:pnode)) - zeror
          element_bb(2,1,ielem) = minval(elcod(2,1:pnode)) - zeror
          element_bb(3,1,ielem) = minval(elcod(3,1:pnode)) - zeror
          element_bb(1,2,ielem) = maxval(elcod(1,1:pnode)) + zeror
          element_bb(2,2,ielem) = maxval(elcod(2,1:pnode)) + zeror
          element_bb(3,2,ielem) = maxval(elcod(3,1:pnode)) + zeror
       end do
    end if

  end subroutine elsest_element_bounding_boxes

!!$subroutine elsest_binpos(&
!!$     ithre,ndime,mnode,npoin,nelem,lnods,ltype,coord,nnode,ipara)
!!$  !
!!$  ! Post-process the bin mesh
!!$  !
!!$  use def_elsest
!!$  use mod_elsest
!!$  implicit none
!!$  integer(ip), intent(in)  :: ithre,ndime,mnode,npoin,nelem
!!$  integer(ip), intent(in)  :: lnods(mnode,nelem),ltype(nelem),nnode(*)
!!$  integer(ip), intent(in)  :: ipara(*)
!!$  real(rp),    intent(in)  :: coord(ndime,npoin)
!!$  integer(ip)              :: ipoin,ielem,idime,kpoin
!!$  integer(ip)              :: box_coord(3),ii,jj,iboxe
!!$ 
!!$  iunit(1) = ipara( 7)
!!$  iunit(2) = ipara(12)
!!$  iunit(3) = ipara(13)
!!$  
!!$  if(iunit(2)>0) then
!!$     !
!!$     ! Coordinates
!!$     !
!!$     ipoin=0
!!$     ielem=0
!!$     if(ndime==2) then
!!$        write(iunit(2),*) 'MESH ELSEST_BIN dimension 2 Elemtype Quadrilateral Nnode 4'
!!$        write(iunit(2),*) 'coordinates'
!!$        do iboxe=1,nboxe
!!$           call elsest_boxcoo(ndime,iboxe,box_coord)
!!$           ii = box_coord(1)
!!$           jj = box_coord(2)
!!$           ipoin=ipoin+1
!!$           write(iunit(2),*) ipoin,comin(1)+real(ii-1)*delta(1),comin(2)+real(jj-1)*delta(2)
!!$           ipoin=ipoin+1
!!$           write(iunit(2),*) ipoin,comin(1)+real(ii)*delta(1),  comin(2)+real(jj-1)*delta(2)
!!$           ipoin=ipoin+1
!!$           write(iunit(2),*) ipoin,comin(1)+real(ii)*delta(1),  comin(2)+real(jj)*delta(2)
!!$           ipoin=ipoin+1
!!$           write(iunit(2),*) ipoin,comin(1)+real(ii-1)*delta(1),comin(2)+real(jj)*delta(2)
!!$        end do
!!$     end if
!!$     do kpoin=1,npoin
!!$        write(iunit(2),*) ipoin+kpoin,(coord(idime,kpoin),idime=1,ndime)
!!$     end do
!!$     write(iunit(2),*)  'end coordinates'
!!$     !
!!$     ! Elements
!!$     !
!!$     ipoin=1
!!$     ielem=0
!!$     write(iunit(2),*) 'elements'
!!$     if(ndime==2) then
!!$        do iboxe=1,nboxe
!!$           write(iunit(2),*) iboxe,ipoin,ipoin+1,ipoin+2,ipoin+3,1
!!$           ipoin=ipoin+4
!!$        end do
!!$     end if
!!$     write(iunit(2),*) 'end elements'
!!$     ipoin=nboxe*4
!!$     call elsest_geogid(ndime,mnode,nelem,nboxe,ipoin,nnode,lnods,ltype)
!!$  end if
!!$
!!$  if(iunit(3)>0) then
!!$     write(iunit(3),*) 'GiD Post Results File 1.0'
!!$     write(iunit(3),*) ' '
!!$     write(iunit(3),*) 'GaussPoints GP_QUAD4 Elemtype Quadrilateral'
!!$     write(iunit(3),*) 'Number of Gauss Points: 1'
!!$     write(iunit(3),*) 'Natural Coordinates: Internal'
!!$     write(iunit(3),*) 'End GaussPoints'
!!$
!!$     write(iunit(3),*) 'Result ELEMENT_NUMBER ELSEST 0 Scalar OnGaussPoints GP_QUAD4'
!!$     write(iunit(3),*) 'ComponentNames ELEMENT_NUMBER'
!!$     write(iunit(3),*) 'values'
!!$     if(dataf==0) then
!!$        do iboxe=1,nboxe
!!$           if( associated(tboel(iboxe)%l) ) then
!!$              write(iunit(3),*) iboxe,size(tboel(iboxe)%l)
!!$           else
!!$              write(iunit(3),*) iboxe,0
!!$           end if
!!$        end do
!!$     else
!!$        do iboxe=1,nboxe
!!$           write(iunit(3),*) iboxe,pboel(iboxe+1)-pboel(iboxe)
!!$        end do
!!$     end if
!!$     write(iunit(3),*) 'end values'
!!$  end if
!!$
!!$end subroutine elsest_binpos

  subroutine elsest_cputim(rtime)
    !-----------------------------------------------------------------------
    !****f* elsest_cputim
    ! NAME
    !    nsi_elmope
    ! DESCRIPTION
    !    Returns the CPU time in seconds
    !    corrections by Alistair Hart <ahart@cray.com> to avoid non standard etime
    ! OUTPUT
    !    rtime
    ! USES
    ! USED BY
    !***
    !-----------------------------------------------------------------------
#ifdef _OPENMP
    use omp_lib
#endif

!!$ CRAY change
!!$   etime() is not part of the Fortran standard
!!$   Some compilers support it, but it would be better to use
!!$   something more standard. 
!!$   In the meantime, we will use MPI_WTIME()
!!$   Macro _CRAYFTN is automatically defined if we are using
!!$   the Cray compiler
#ifdef _CRAYFTN
    use MPI
#endif
!!$ END CRAY change

    implicit none
    real(rp), intent(out) :: rtime

#ifdef _OPENMP
    rtime = OMP_GET_WTIME()
#else
!!$ CRAY change
#ifdef _CRAYFTN
    rtime = MPI_WTIME()
#else
    call cpu_time(rtime)
#endif
!!$ END CRAY change
#endif

  end subroutine elsest_cputim

  subroutine elsest_runend(message)
    !-----------------------------------------------------------------------
    !
    ! This routine stops the run and writes the summary of CPU time.
    !
    !-----------------------------------------------------------------------

    character(*) :: message
    !
    ! Write message and stop the run
    !
    if(iunit(1)/=0) then
       !*OMP CRITICAL(write)
       write(iunit(1),1) 'ELSEST ERROR WAS FOUND: '//trim(message)
       !*OMP END CRITICAL(write)
    else
       !*OMP CRITICAL(write)
       write(6,1)     'ELSEST ERROR WAS FOUND: '//trim(message)
       !*OMP END CRITICAL(write)
    end if

    stop

1   format(//,5x,'>>> ',a)

  end subroutine elsest_runend

end module mod_elsest
!> @}
