!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




!-----------------------------------------------------------------------
!> @addtogroup Maths
!> @{
!> @file    def_coupling_method.f90
!> @author  SSantoso
!> @date    2022-02-09
!> @brief   Class that regroups interpolation and projection. 
!> @details In the interpolation module, the interpolation interp_matrix is 
!>contained in the source process. 
!>In the projection module, the mass matrix is contained in the target process.
!>In coupling we want all interp_matrix to be in the target process.
!>The interpolation interp_matrix has to be sent to the target thanks to the communicator
!>of def_search_parall. However in the coupling, we want the source to communicate
!>only nodal values. Then, a new communicator has to be built in here.
!>The subroutine of the projection remains the same for the moment
!>
!>
!-----------------------------------------------------------------------

module def_coupling_method

  use def_kintyp_basic,                   only : ip,rp,lg,i1p,i2p,r1p,r2p,r3p
  use def_kintyp_mesh_basic,              only : mesh_type_basic,bmsh_type_basic 
  use def_kintyp_comm,                    only : comm_data_par_basic
  use def_mat_coo,                        only : mat_coo
  use def_search_parall     
  use def_search_method,                  only : search_method,CANDIDATE_INSIDE,CANDIDATE_NEAREST,SEARCH_POINTS,&
                                                 SEARCH_BOUNDING_BOXES
  use mod_memory_basic,                   only : memory_alloca,memory_deallo,memory_size,memory_resize
  use mod_memory_tools,                   only : memory_counter_ini,memory_counter_end
  use mod_optional_argument,              only : optional_argument
  use mod_communications_global,          only : PAR_ALLGATHER,PAR_SUM,PAR_MIN,PAR_MAX,PAR_AVERAGE
  use mod_communications_point_to_point,  only : PAR_SEND_RECEIVE
  use mod_communications,                 only : PAR_BARRIER
  use mod_maths_arrays,                   only : maths_maxloc_nonzero
  use mod_elmgeo,                         only : element_type,elmgeo_natural_coordinates,elmgeo_cartesian_derivatives,&
                                                 elmgeo_natural_coordinates_on_boundaries,elmgeo_projection_on_a_face
  use mod_htable,                         only : hash_t,htaini,htaadd,htades,htalid
  
  use def_interpolation_method
  use def_projection_method
  use mod_std
  use mod_communications_tools,          only : PAR_COMM_RANK_AND_SIZE
  use def_domain
  use def_master
  use mod_comm_basic,                     only : par_interface_exchange
  use mod_parall
  use mod_maths_solver, only : maths_conjugate_gradient
  use def_iterative_solvers
  use def_mat_csr
  use def_mat_dia
  use def_direct_solvers
  use mod_driver_solvers
  use def_solver
  use def_solvers
  use def_all_solvers
  use def_master
  use def_communications
  use mod_communications
  use def_preconditioners
  use def_mpi
#include "def_mpi.inc"
  implicit none
  private

  real(rp),      parameter :: epsil = epsilon(1.0_rp)
  character(24), parameter :: vacal = 'def_coupling_method'

  type :: coupling_met
     character(LEN=:),  allocatable     :: name                           ! Name


     type(mat_coo),             pointer :: interp_matrix(:)               ! Interpolation interp_matrix
     type(mat_coo),             pointer :: interp_matder(:)               ! Interpolation interp_matrix for derivatives
     type(mat_coo)                      :: mass_matrix                    ! Interpolation interp_matrix for derivatives
     type(interpolation)                :: my_interpolation               ! Type of interpolation
     integer(ip)                        :: interpolation_method
     type(projection)                   :: my_projection                  ! Type of coupling
     type(comm_data_par_basic)          :: comm                           ! my communicator
     type(i1p),                 pointer :: point_to_send(:)               !
     integer(ip)                        :: nn                             ! Number of points 
     integer(ip)                        :: nd                             ! Dimension
     integer(ip)                        :: nrank                          ! Number of ranks involved
     integer(ip)                        :: mode
     class(search_method),       pointer :: search_method_seq              ! Sequential search
     !class(search_method),      pointer :: search_method_par              ! Parallel search
     type(search_parall)                :: parallel_search                ! Parallel search
     real(rp)                           :: toler_rel                      ! Relative tolerance
     logical(lg)                        :: deriv                          ! If derivative should be computed
     logical(lg)                        :: lg
     logical(lg),               pointer :: found(:)                       ! It point has been found
     integer(8)                         :: memor(2)                       ! Memory counter
     real(rp)                           :: times(10)                      ! Timings
     real(rp)                           :: stats(10)                      ! Statistics
     logical(lg) :: myself   
   contains
     procedure,               pass      :: init
     procedure,               pass      :: input
     procedure,               pass      :: deallo
     procedure,               pass      :: preprocess
     procedure,               pass      :: coupling_entity
     procedure,               pass      :: coupling_boundary
     procedure,               pass      :: coupling_nearest_node
     procedure,               pass      :: coupling_global_numbering
     procedure,               pass      :: communicate_matrixs
     procedure,               pass      :: compute_interp_matrix
     procedure,               pass      :: coupling_sparse_matrix
     procedure,               pass      :: build_comm
     procedure,               pass      :: values_interp1d
     procedure,               pass      :: values_interp2d  
     procedure,               pass      :: values_proj
     generic                            :: values_interp => &
          &                                values_interp1d, &
          &                                values_interp2d    

     
  end type coupling_met

  abstract interface
   subroutine interpolation_generic(self,search,xx,mesh,lelem,shapf,deriv,dista,lenty,mask,TOLER,MEMORY_COUNTER)
   import                                                 :: coupling_met
   import                                                 :: search_method
   import                                                 :: mesh_type_basic
   import                                                 :: rp
   import                                                 :: ip
   import                                                 :: lg
   class(coupling_met),                    intent(inout) :: self
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
  
  public :: coupling_met
  public :: INT_SEQUENTIAL
  public :: INT_PARALLEL
  public :: INT_BIN                   
  public :: INT_OCTREE                
  public :: INT_OCTBIN                
  public :: INT_KDTREE                
  public :: INT_ELEMENT_INTERPOLATION                       
  public :: INT_ELEMENT_VALUE                       
  public :: INT_BOUNDARY_INTERPOLATION
  public :: INT_NEAREST_NODE 
  public :: INT_GLOBAL_NUMBERING
  public :: BOUNDARY_PROJECTION

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

    class(coupling_met), intent(inout) :: self

    nullify(self % interp_matrix)
    nullify(self % interp_matder)
    nullify(self%point_to_send) 
    self % memor                             = 0_8
    self % times                             = 0.0_rp
    self % stats                             = 0.0_rp
    call self%my_interpolation%init ()
    call self%my_projection%init()   
    call self%comm%init()
    call self % parallel_search % init()
  end subroutine init

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-16
  !> @brief   Deallocate
  !> @details Deallocate
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo(self,MEMORY_COUNTER)

    class(coupling_met),           intent(inout) :: self
    integer(ip)                                   :: ii
    integer(8),           optional, intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(8)                                    :: memor_loc(2)

    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)
    
    
    call self%my_interpolation%deallo()
    call self%my_projection%deallo()
    
    if( associated(self % interp_matrix) ) then
      do ii = lbound(self % interp_matrix,1),ubound(self % interp_matrix,1)
        call self % interp_matrix(ii) % deallo(MEMORY_COUNTER=memor_loc)
      end do
      deallocate(self % interp_matrix)
    end if
    

    
    
    
    
    
    if( associated(self % interp_matder) ) then
       do ii = lbound(self % interp_matder,1),ubound(self % interp_matder,1)
          call self % interp_matder(ii) % deallo(MEMORY_COUNTER=memor_loc)
       end do
       deallocate(self % interp_matder)
    end if
    
    call self % mass_matrix % deallo(MEMORY_COUNTER=memor_loc)


    if( allocated(self % name) ) deallocate(self % name)
    call self%comm % deallo()
      
      
    if(associated(self%point_to_send)) then 
      call memory_deallo(memor_loc,'self % point_to_send',vacal,self%point_to_send) 
      !do ii = lbound(self % point_to_send,1),ubound(self % point_to_send,1)
      !  call memory_deallo(memor_loc,'self%point_to_send%l',vacal,self%point_to_send(ii)%l)
      !end do
      !deallocate(self % point_to_send)
    end if
    
    
        
    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)

  end subroutine deallo

  
  !-----------------------------------------------------------------------
  !> 
  !> @author  SSantoso
  !> @date    2022-02-10
  !> @brief   Input of interpolation and projection
  !> @details Input of interpolation and projection
  !> 
  !-----------------------------------------------------------------------

  subroutine input(self,            &
    mesh,                        &
    search_method_seq,           &
    search_method_par,           &
    COMM,                        &
    MODE,                        &
    INTERPOLATION_METHOD,        &
    NAME,                        &
    DERIVATIVES,                 &
    MYSELF                       )

    class(coupling_met),                    intent(inout) :: self
    class(mesh_type_basic),                 intent(in)    :: mesh
    class(search_method), optional, target, intent(in)    :: search_method_seq
    class(search_method), optional, target, intent(in)    :: search_method_par
    MY_MPI_COMM   ,       optional,         intent(in)    :: COMM
    integer(ip),          optional,         intent(in)    :: MODE
    integer(ip),          optional,         intent(in)    :: INTERPOLATION_METHOD
    character(LEN=*),     optional,         intent(in)    :: NAME
    logical(lg),          optional,         intent(in)    :: DERIVATIVES
    logical(lg),          optional,         intent(in)    :: MYSELF
    type(bmsh_type_basic), pointer                        :: mesh_target
    type(bmsh_type_basic), pointer                        :: mesh_source
          
    call self % parallel_search % input(COMM)

    if(present(search_method_par)) then
      self%mode = INT_PARALLEL
   end if

    if( self % mode == INT_PARALLEL ) then
       self % nrank = int(self % parallel_search % comm % size4,ip)
    else
       self % nrank = 1
       self % parallel_search % comm % RANK4 = 0_4
    end if
    
    if(present(search_method_par)) then
      self%search_method_seq => search_method_seq
      self % parallel_search % search_method_par => search_method_par
    else 
      self%search_method_seq => search_method_seq
    end if  
    self % comm % size4          = self % parallel_search % comm % size4
    self % comm % rank4          = self % parallel_search % comm % rank4
    self % comm % PAR_COMM_WORLD = self % parallel_search % comm % PAR_COMM_WORLD
    
    self % interpolation_method = INTERPOLATION_METHOD
    self % deriv = .false.
    !
    ! Boundary mesh
    !
    mesh_source => mesh % boundary
    mesh_target => mesh % boundary
    call self%my_projection%input(BOUNDARY_PROJECTION,self%my_interpolation,elmar,mesh_source,mesh_target,COMM=commd)
    
  end subroutine input
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  ssantoso
  !> @date    2022-02-09
  !> @brief   Initialize
  !> @details Initialize
  !> 
  !-----------------------------------------------------------------------

  subroutine preprocess(self,mesh,ll,mask,MEMORY_COUNTER)
    
    class(coupling_met),                     intent(inout) :: self
    class(mesh_type_basic),                  intent(in)    :: mesh
    integer(ip),          optional, pointer, intent(in)    :: ll(:)
    logical(lg),          optional, pointer, intent(in)    :: mask(:)
    integer(8),           optional,          intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(ip)                                            :: nrank,irank,size_pts,my_rank
    integer(8)                                             :: memor_loc(2)
    integer(8)                                             :: memor_dum(2)
    type(bmsh_type_basic), pointer                         :: mesh_target
    type(mat_coo),         pointer                         :: source_matrix(:)
    
    my_rank =int(self%comm % RANK4,ip) 
    call memory_counter_ini(memor_loc,self% memor,MEMORY_COUNTER)
    memor_dum = 0_8
    nrank        = self % nrank
     
    allocate(self % interp_matrix(0:nrank-1))
    do irank = 0,nrank-1
       call self % interp_matrix(irank) % init()
    end do    
  
    !calculer les points d'intégrations
    mesh_target => mesh % boundary
    call self%my_projection%integration_points(self%my_projection%mesh_target,self%my_projection%elmar,mask)
    !print*,"Simon//",my_rank,self%my_projection %xg 
    !Calculer les matrices d'interpolations côté source
    call self%compute_interp_matrix(source_matrix,self%my_projection %xg,mesh_target) 
    self % nn = self%my_interpolation%nn 
    size_pts = 0
    call self % mass_matrix % init()
    
    if(IPARALL) then 
     !allocate(self % point_to_send(0:nrank-1))
     nullify(self%point_to_send)
     call memory_alloca(memor_loc,'self % point_to_send',vacal,self%point_to_send,nrank,LBOUN=0_ip)
     do irank = 0, nrank-1
       size_pts=0
       size_pts= memory_size(source_matrix(irank)%ya_conc)
       print*,"Simon//",my_rank,irank,size_pts
       nullify(self%point_to_send(irank)%l)
       if(size_pts > 0 ) then
         call memory_alloca(memor_loc,'self%point_to_send%l',vacal,self%point_to_send(irank)%l,size_pts)
         self%point_to_send(irank)%l(1:size_pts) = source_matrix(irank)%ya_conc(1:size_pts)
       end if
     end do
     call self%communicate_matrixs(source_matrix,MEMORY_COUNTER=memor_loc)
     call self%build_comm()
    else
     self % interp_matrix(0) % ndof1 = 1
     self % interp_matrix(0) % ndof2 = 1
     self % interp_matrix(0) % nz    = source_matrix(0)%nz
     call self % interp_matrix(0) % alloca(MEMORY_COUNTER=memor_loc)
     self%interp_matrix(0)%xA(:)=source_matrix(0)%xA(:)
     self%interp_matrix(0)%yA(:)=source_matrix(0)%yA(:)
     self%interp_matrix(0)%vA(1,1,:)=source_matrix(0)%vA(1,1,:)  

  end if
    !
    ! Compute projection matrix on target
    !
    call self%my_projection%target_matrix(mesh_target,elmar,mask,lumped =.false.)

    !On transmet les matrices de masse dans coupling
    if(memory_size(self%my_projection%matrix%xA) >0) then
     self % mass_matrix % ndof1 = 1
     self % mass_matrix % ndof2 = 1
     self % mass_matrix % nz    = self%my_projection%matrix%nz
     self % mass_matrix % nrows = self%my_projection%matrix % nrows
     call self % mass_matrix % alloca(MEMORY_COUNTER=memor_loc)
     self%mass_matrix%xA(:)=self%my_projection%matrix%xA(:)
     self%mass_matrix%yA(:)=self%my_projection%matrix%yA(:)
     self%mass_matrix%vA(1,1,:)=self%my_projection%matrix%vA(1,1,:)
    end if 
    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)

  end subroutine preprocess
  
   !-----------------------------------------------------------------------
  !> 
  !> @author  SSantoso
  !> @date    2020-10-16
  !> @brief   compute_interp_matrix
  !> @details Subroutine made to communicate matrix in the target process
  !> 
  !-----------------------------------------------------------------------

  subroutine compute_interp_matrix(self,source_matrix,xx,mesh,ll,mask,MEMORY_COUNTER) 
  
    class(coupling_met),                     intent(inout) :: self
    type(mat_coo),                  pointer                :: source_matrix(:)
    real(rp),                       pointer, intent(in)    :: xx(:,:)
    class(mesh_type_basic),                  intent(in)    :: mesh
    integer(ip),          optional, pointer, intent(in)    :: ll(:)
    logical(lg),          optional, pointer, intent(in)    :: mask(:)
    integer(8),           optional,          intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(ip)                                            :: irank,nrank,nn,nd
    type(r2p),                      pointer                :: shapf(:)
    type(r3p),                      pointer                :: deriv(:)
    type(i1p),                      pointer                :: lelem(:),lista_recv(:)
    type(i2p),                      pointer                :: lenty(:)
    logical(lg),                    pointer                :: llost(:)
    real(rp)                                               :: toler
    integer(ip)                                            :: ii,nlost
    integer(8)                                             :: memor_loc(2) 
    
    
    
    if( .not. associated(self % search_method_seq) ) &
         call runend('DEF_COUPLING_METHOD: SEQUENTIAL SEARCH METHOD NOT DEFINED')
    if( self % mode == INT_PARALLEL ) then
       if( .not. associated(self % parallel_search%search_method_par) )  &
            call runend('DEF_COUPLING_METHOD: PARALLEL SEARCH METHOD NOT DEFINED')
    end if   
    
    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)
    nullify(shapf)
    nullify(deriv)
    nullify(lelem)
    nullify(lenty)
    nullify(llost)
    nullify(lista_recv)
    
       
    toler     = self % toler_rel
    nn        = memory_size(xx,2_ip)
    nd        = max(mesh % ndime,memory_size(xx,1_ip))
    self % nn = nn
    self % nd = nd
    nrank     = self % nrank
        
    if( self % mode == INT_PARALLEL ) self % parallel_search % nd = nd
    if( nd /= memory_size(xx,1_ip) .and. memory_size(xx,1_ip) > 0 ) &
         call runend('DEF_COUPLING_METHOD: INCOMPATIBLE DIMENSIONS')

    !--------------------------------------------------------------------
    !
    ! Allocate memory and initialize
    ! By default, points are lost => LLOST(II) = .true.
    !
    !-------------------------------------------------------------------- 
    call memory_alloca(memor_loc,'SHAPF'        ,vacal,shapf        ,nrank,LBOUN=0_ip)
    call memory_alloca(memor_loc,'DERIV'        ,vacal,deriv        ,nrank,LBOUN=0_ip)
    call memory_alloca(memor_loc,'LELEM'        ,vacal,lelem        ,nrank,LBOUN=0_ip)
    call memory_alloca(memor_loc,'LENTY'        ,vacal,lenty        ,nrank,LBOUN=0_ip)
    call memory_alloca(memor_loc,'LISTA_RECV'   ,vacal,lista_recv,nrank,LBOUN=0_ip)
    call memory_alloca(memor_loc,'LLOST'        ,vacal,llost        ,nn)   


    allocate(source_matrix(0:nrank-1))
    do irank = 0,nrank-1
       call source_matrix(irank) % init()
    end do
    !if( self % deriv ) then
    !   allocate(self % interp_matder(0:nrank-1))
    !   do irank = 0,nrank-1
    !      call self % interp_matder(irank) % init()
    !   end do
    !end if
    do ii = 1,nn
       llost(ii) = .true.
    end do

    !--------------------------------------------------------------------
    !
    ! Interpolate
    !
    !--------------------------------------------------------------------

   select case ( self % interpolation_method )
       
    case ( INT_ELEMENT_INTERPOLATION , INT_ELEMENT_VALUE , INT_BOUNDARY_INTERPOLATION )
       !
       ! Element and boundary interpolations, and element value
       !
       call self % coupling_entity(xx,mesh,lista_recv,shapf,deriv,lelem,lenty,llost,mask,MEMORY_COUNTER=memor_loc)
       
    case ( INT_NEAREST_NODE )
       !
       ! Nearest node
       !
       call self % coupling_nearest_node(xx,mesh,shapf,deriv,lelem,lenty,llost,mask,MEMORY_COUNTER=memor_loc)
       
    case ( INT_GLOBAL_NUMBERING )
       !
       ! Global numbering
       !
       if( present(ll) ) then
          call self % coupling_global_numbering(xx,ll,mesh,shapf,deriv,lelem,lenty,llost,mask,MEMORY_COUNTER=memor_loc)
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
    !call self % sparse_matrix(lelem,lenty,shapf,deriv,lista_recv,MEMORY_COUNTER=memor_loc)

    call self%coupling_sparse_matrix(source_matrix,lelem,lenty,shapf,deriv,lista_recv,MEMORY_COUNTER=memor_loc)
    
    nullify(self % found)
    call memory_alloca(memor_loc,'SELF % FOUND',vacal,self % found,nn)
    if(nn>0) self % found = .not. llost
   

    !--------------------------------------------------------------------
    !
    ! Statistics
    !
    !--------------------------------------------------------------------
    
    nlost = 0
    if( nn > 0 ) nlost = count(llost)
    self % stats(1) = real(nlost,rp)

    !--------------------------------------------------------------------
    !
    ! Deallocate 
    !
    !--------------------------------------------------------------------
    call memory_deallo(memor_loc,'LLOST'        ,vacal,llost     )
    call memory_deallo(memor_loc,'LISTA_RECV'   ,vacal,lista_recv)
    call memory_deallo(memor_loc,'LENTY'        ,vacal,lenty     )
    call memory_deallo(memor_loc,'LELEM'        ,vacal,lelem     )
    call memory_deallo(memor_loc,'DERIV'        ,vacal,deriv     )
    call memory_deallo(memor_loc,'SHAPF'        ,vacal,shapf     )
    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)    
    
    
    
    
  end subroutine compute_interp_matrix
  
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  ssantoso
  !> @date    2021-01-29
  !> @brief   Element interpolation
  !> @details Element interpolation
  !> 
  !-----------------------------------------------------------------------

  subroutine coupling_entity(&
       self,xx,mesh,lista_recv,shapf,deriv,lelem,lenty,llost,mask,MEMORY_COUNTER)

    class(coupling_met),                     intent(inout) :: self
    real(rp),                       pointer, intent(in)    :: xx(:,:)
    class(mesh_type_basic),                  intent(in)    :: mesh
    type(i1p),                      pointer, intent(inout) :: lista_recv(:)    
    type(r2p),                      pointer                :: shapf(:)
    type(r3p),                      pointer                :: deriv(:)
    type(i1p),                      pointer                :: lelem(:)
    type(i2p),                      pointer                :: lenty(:)
    logical(lg),                    pointer                :: llost(:)
    logical(lg),          optional, pointer, intent(in)    :: mask(:)
    integer(8),           optional,          intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    type(r1p),                      pointer                :: dista_recv(:)    
    type(r2p),                      pointer                :: xx_recv(:)
    integer(ip),                    pointer                :: nn_send(:)
    integer(ip),                    pointer                :: nn_recv(:)
    type(i1p),                      pointer                :: lista_send(:)  
    real(rp)                                               :: toler
    integer(ip)                                            :: ii,kk,irank
    integer(ip)                                            :: ineig,ipoin
    integer(ip)                                            :: my_rank,nn,nrank
    integer(ip)                                            :: kfl_method
    integer(8)                                             :: memor_loc(2)
    logical(lg)                                            :: myself_first
    procedure(interpolation_generic),      pointer         :: coupling_what
    
    

    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)
    nullify(dista_recv)    
    nullify(xx_recv)    
    nullify(nn_send)
    nullify(nn_recv)
    nullify(lista_send)  

    !--------------------------------------------------------------------
    !
    ! What to search
    !
    !--------------------------------------------------------------------

    if(      self % interpolation_method == INT_BOUNDARY_INTERPOLATION ) then
       coupling_what => coupling_boundary
       kfl_method         =  CANDIDATE_NEAREST
    else if( self % interpolation_method == INT_ELEMENT_INTERPOLATION  ) then
       coupling_what => coupling_element
       kfl_method         =  CANDIDATE_INSIDE
    else if( self % interpolation_method == INT_ELEMENT_VALUE          ) then
       coupling_what => coupling_element
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
    toler        = self  % toler_rel
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
    if( myself_first .or. self  % mode == INT_SEQUENTIAL ) then
       call coupling_what(&
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
          if( lelem(my_rank) % l(ii) /= 0 ) llost(ii) = .false.
       end do
       if( nn > 0 ) then
          self % stats(3) = real(nn-count(llost),rp)
       end if
    end if  

    !--------------------------------------------------------------------
    !
    ! Parallel search
    !
    !--------------------------------------------------------------------
    if( self % mode == INT_PARALLEL ) then
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
       end if

       do irank = 0,nrank-1
          if( nn_recv(irank) > 0 .and. ( irank /= my_rank .or. (.not. myself_first) ) ) then
             call coupling_what(&
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
       call self % parallel_search % dista_comm(nn,nn_send,nn_recv,lista_send,lista_recv,dista_recv,MEMORY_COUNTER=memor_loc) 
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
    !--------------------------------------------------------------------
    !
    ! Deallocate
    !
    !--------------------------------------------------------------------
    call memory_deallo(memor_loc,'DISTA_RECV'   ,vacal,dista_recv)
    call memory_deallo(memor_loc,'LISTA_SEND'   ,vacal,lista_send)
    call memory_deallo(memor_loc,'XX_RECV'      ,vacal,xx_recv   )
    call memory_deallo(memor_loc,'NN_SEND'      ,vacal,nn_send   )
    call memory_deallo(memor_loc,'NN_RECV'      ,vacal,nn_recv   )
    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)
    
    
  end subroutine coupling_entity
  
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-29
  !> @brief   Nearest node
  !> @details Nearest node interpolation
  !> 
  !-----------------------------------------------------------------------

  subroutine coupling_nearest_node(&
    self,xx,mesh,shapf,deriv,lelem,lenty,llost,mask,MEMORY_COUNTER,search_seq_in,search_par_in)

    class(coupling_met),                     intent(inout) :: self
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
    toler   = self % toler_rel

    nrank   = self % nrank
    my_rank = int(self % parallel_search % comm % RANK4,ip)

    call memory_alloca(memor_loc,'DISTA_RECV',vacal,dista_recv,nrank,LBOUN=0_ip)
    !--------------------------------------------------------------------
    !
    ! Search
    !
    !--------------------------------------------------------------------
    if( self  % mode == INT_PARALLEL ) then

#ifdef __PGI
       call self % parallel_search % send_list_recv_points(&
            XX=xx,NN_SEND=nn_send,NN_RECV=nn_recv,LISTA_SEND=lista_send,XX_RECV=xx_recv,&
            METHOD=CANDIDATE_NEAREST,MASK=std_log_1,MEMORY_COUNTER=memor_loc)
#else
       call self % parallel_search % send_list_recv_points(&
            xx,nn_send,nn_recv,lista_send,xx_recv,&
            METHOD=CANDIDATE_NEAREST,MEMORY_COUNTER=memor_loc)
#endif
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

  end subroutine coupling_nearest_node
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-29
  !> @brief   Nearest node and global numbering
  !> @details Global numbering
  !> 
  !-----------------------------------------------------------------------

  subroutine coupling_global_numbering(&
    self,xx,ll,mesh,shapf,deriv,lelem,lenty,llost,mask,MEMORY_COUNTER,search_seq_in,search_par_in)

    class(coupling_met),                    intent(inout) :: self
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
!    class(search_method),           pointer                :: search_par        ! Sequential search
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
    
    !if( present(search_seq_in) ) then
    !   search_seq => search_seq_in
    !else
    !   search_seq => self % search_method_seq
    !end if
   ! 
   ! if( present(search_par_in) ) then
   !    search_par => search_par_in
   ! else
   !    search_par => self % parallel_search % search_method_par
   ! end if
    
    !--------------------------------------------------------------------
    !
    ! Initialization
    !
    !--------------------------------------------------------------------

    nn      = memory_size(xx,2_ip)
    toler   = self  % toler_rel
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

    if( self  % mode == INT_PARALLEL ) then

#ifdef __PGI
       call self % parallel_search % send_list_recv_points(&
            XX=xx,LL=ll,NN_SEND=nn_send,NN_RECV=nn_recv,LISTA_SEND=lista_send,XX_RECV=xx_recv,&
            LL_RECV=ll_recv,METHOD=CANDIDATE_NEAREST,MASK=std_log_1,&
            EXCLUDE_MYSELF=.not.self % myself,MEMORY_COUNTER=memor_loc)
#else
       call self % parallel_search % send_list_recv_points(&
            xx,ll,nn_send,nn_recv,lista_send,xx_recv,ll_recv,&
            METHOD=CANDIDATE_NEAREST,&
            EXCLUDE_MYSELF=.not.self % myself,MEMORY_COUNTER=memor_loc)
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

  end subroutine coupling_global_numbering
  
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
  
  subroutine coupling_element(self,search,xx,mesh,lelem,shapf,deriv,dista,lenty,mask,TOLER,MEMORY_COUNTER)

    class(coupling_met),                     intent(inout) :: self
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
       call coupling_allocate(self,nn,ndime,mnode,lelem,shapf,deriv,dista,lenty,MEMORY_COUNTER=memor_loc)

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
          if( present(dista) ) dista(ii) = -1.0_rp

          if( nelem == 0 ) then
             lelem(ii)   =  0
             shapf(:,ii) =  0.0_rp
             if( present(deriv) ) deriv(:,:,ii) = 0.0_rp
          else
            loop_kk: do kk = 1,nelem
               ielem = list_entities(ii) % l(kk)
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
          end if
       end do

       self % stats(2) = self % stats(2) / real(nn,rp)
       call cputim(time2) ; self % times(2) = self % times(2) + time2-time1 ; time1 = time2
       call memory_deallo(memor_loc,'LIST_ENTITIES',vacal,list_entities)       
    end if

    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)
    
    
  end subroutine coupling_element

  subroutine coupling_boundary(self,search,xx,mesh,lelem,shapf,deriv,dista,lenty,mask,TOLER,MEMORY_COUNTER)

    class(coupling_met),                     intent(inout) :: self
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
       if( present(deriv) ) then
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
       call coupling_allocate(self,nn,ndime,max(mnode,mnodb),lelem,shapf,deriv,dista,lenty,MEMORY_COUNTER=memor_loc)

       
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
             if( present(dista) ) dista(ii) = -1.0_rp
             if( nboun == 0 ) then
                lelem(ii)   =  0
                shapf(:,ii) =  0.0_rp
                if( present(deriv) .and. self % deriv ) deriv(:,:,ii) = 0.0_rp 
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

                      if( ifoun /= 0 ) then
                         proje(1:ndime) = pnear(1:ndime)
                         dista_loc      = sqrt(dot_product(proje(1:ndime)-xx(1:ndime,ii),proje(1:ndime)-xx(1:ndime,ii)))

                         if( dista_loc < dista_min ) then
                            lelem(ii)         = iboun
                            shapf(1:pnodb,ii) = shapt(1:pnodb)
                            dista_min         = dista_loc
                            if( present(lenty) ) then
                               lenty(1:pnodb,ii) = lnodb(1:pnodb,iboun)
                            end if
                            if( present(dista) ) dista(ii) = dista_loc
                            !
                            ! Derivatives not available
                            !
                            if( present(deriv) .and. self % deriv ) then                       
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

  end subroutine coupling_boundary

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-05
  !> @brief   Nearest node
  !> @details Nearest node
  !> 
  !-----------------------------------------------------------------------
  
  subroutine nearest_node(self,search,xx,mesh,lelem,shapf,deriv,dista,lenty,mask,TOLER,MEMORY_COUNTER)

    class(coupling_met),                     intent(inout) :: self
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
       call coupling_allocate(self,nn,ndime,mnode,lelem,shapf,deriv,dista,lenty,MEMORY_COUNTER=memor_loc)

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
             if( present(dista) ) dista(ii) = -1.0_rp
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

    class(coupling_met),                     intent(inout) :: self
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
       call coupling_allocate(self,nn,ndime,mnode,LELEM=lelem,SHAPF=shapf,DERIV=deriv,LENTY=lenty,MEMORY_COUNTER=memor_loc)

       do ii = 1,nn
          ipoin_global = abs(ll(ii))
          ipoin_local  = htalid(ht,ipoin_global)
          !block ; use def_master ;  if(ipoin_global==26) print*,'AQUI=',kfl_paral,ipoin_local ;end block
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
  
  subroutine coupling_allocate(self,nn,ndime,mnode,lelem,shapf,deriv,dista,lenty,MEMORY_COUNTER)
  
    class(coupling_met),                     intent(inout) :: self
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
    if( .not. associated(lelem)    ) then 
      call memory_alloca(memor_loc,'LELEM',vacal,lelem,nn)
    end if
    if( .not. associated(shapf)    ) then
      call memory_alloca(memor_loc,'SHAPF',vacal,shapf,mnode,nn)
    end if
    if( present(deriv) ) then
      if( .not. associated(deriv) ) then
        call memory_alloca(memor_loc,'DERIV',vacal,deriv,ndime,mnode,nn)
      end if
    end if
    if( present(dista) ) then
      if( .not. associated(dista) ) then 
        call memory_alloca(memor_loc,'DISTA',vacal,dista,nn)
      end if
    end if
    if( present(lenty) ) then
      if( .not. associated(lenty) ) then 
        call memory_alloca(memor_loc,'LENTY',vacal,lenty,mnode,nn)
      end if
    end if
    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)


  end subroutine coupling_allocate
  
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

  subroutine coupling_sparse_matrix(self,source_matrix,lelem,lenty,shapf,deriv,lista_recv,MEMORY_COUNTER)

    class(coupling_met),                        intent(inout) :: self
    type(mat_coo),                    pointer,  intent(in)    :: source_matrix(:)
    type(i1p),                        pointer,  intent(in)    :: lelem(:)
    type(i2p),                        pointer,  intent(in)    :: lenty(:)
    type(r2p),                        pointer,  intent(in)    :: shapf(:)
    type(r3p),                        pointer,  intent(in)    :: deriv(:)
    type(i1p),                        pointer,  intent(in)    :: lista_recv(:)
    integer(8),             optional,           intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(ip)                                               :: irank,nz,nd,nrank
    integer(ip)                                               :: ii,kk,ielem
    integer(ip)                                               :: pnode
    integer(ip)                                               :: ipoin,inode
    integer(8)                                                :: memor_loc(2)
    
    

    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)
    nrank = self % nrank
    nd = self % nd
    do irank = 0,self % nrank-1
       !
       ! Compute dimension
       !
       nz = 0
       if( self % interpolation_method == INT_ELEMENT_VALUE ) then
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
       source_matrix(irank) % ndof1 = 1
       source_matrix(irank) % ndof2 = 1
       source_matrix(irank) % nz    = nz
       call source_matrix(irank) % alloca(MEMORY_COUNTER=memor_loc)

       !if( self  % deriv ) then
       !   self % matder(irank) % ndof1 = 1
       !   self % matder(irank) % ndof2 = nd
       !   self % matder(irank) % nz    = nz
       !   call self % matder(irank) % alloca(MEMORY_COUNTER=memor_loc)
       !end if
       !
       ! Fill in matrix
       !
       nz = 0
       kk = 0

       if( self  % interpolation_method == INT_ELEMENT_VALUE ) then
          do ii = 1,memory_size(lelem(irank) % l)
             ielem = lelem(irank) % l(ii)
             if( ielem /= 0 ) then
                kk                                       = kk + 1
                nz                                       = nz + 1                

                source_matrix(irank) % xA(nz)            = kk
                source_matrix(irank) % yA(nz)            = ielem
                source_matrix(irank) % vA(1,1,nz)        = 1.0_rp
                source_matrix(irank) % nrows             = max(source_matrix(irank) % nrows,kk)
                source_matrix(irank) % ncols             = max(source_matrix(irank) % ncols,ielem)
                
                !if( self  % deriv ) then
                !   self % matder(irank) % xA(nz)         = kk
                !   self % matder(irank) % yA(nz)         = ielem
                !   self % matder(irank) % vA(1,1:nd,nz)  = 0.0_rp
                !   self % matder(irank) % nrows          = max(self % matder(irank) % nrows,kk)
                !   self % matder(irank) % ncols          = max(self % matder(irank) % ncols,ielem)
                !end if
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

                   !self % matrix(irank) % xA(nz)            = kk
                   if(self  % mode == INT_PARALLEL) then
                     source_matrix(irank) % xA(nz)            = lista_recv(irank)%l(ii)
                   else 
                     source_matrix(irank) % xA(nz)            = ii
                   end if

                   
                   
                   
                   source_matrix(irank) % yA(nz)            = ipoin
                   source_matrix(irank) % vA(1,1,nz)        = shapf(irank) % a(inode,ii) 
                   source_matrix(irank) % nrows             = max(source_matrix(irank) % nrows,kk)
                   source_matrix(irank) % ncols            = max(source_matrix(irank) % ncols,ipoin) !!! Pas sûr de ça
                   
                   
                   
                   !if( self % deriv ) then
                   !   self % matder(irank) % xA(nz)         = kk
                   !   self % matder(irank) % yA(nz)         = ipoin
                   !   self % matder(irank) % vA(1,1:nd,nz)  = deriv(irank) % a(1:nd,inode,ii)
                   !   self % matder(irank) % nrows          = max(self % matder(irank) % nrows,kk)
                   !   self % matder(irank) % ncols          = max(self % matder(irank) % ncols,ipoin)
                   !end if !! present deriv
                end do ! inode
             end if !! elem different of zero)
          end do ! ii 
          if(nz > 0) then
            call source_matrix(irank) % compute_perm_and_conc(MEMORY_COUNTER=memor_loc)
         end if ! nz > 0           
       end if ! interpolation method 
    end do ! irank

    
    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)
    
    
  end subroutine coupling_sparse_matrix
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  SSantoso
  !> @date    2020-10-16
  !> @brief   communicate_matrixs
  !> @details Subroutine made to communicate matrix in the target process
  !> 
  !-----------------------------------------------------------------------

  subroutine communicate_matrixs(self,interp_matrix_to_send,MEMORY_COUNTER)
    
    class(coupling_met),         intent(inout) :: self
    type(mat_coo), pointer,  intent(in)    :: interp_matrix_to_send(:) ! interp matrix to comm
    integer(8),    optional, intent(inout) :: MEMORY_COUNTER(2)        !< Memory counters
    !-----------------------------------------------------------------------------------!
    integer(ip) :: irank,ineig,my_rank,nsend,nrecv,nrank
    integer(ip) :: size_buff1,size_buff2
    integer(ip), pointer :: send_size(:),recv_buff(:),send_buff(:)
    real(rp), pointer :: send_buff_rp(:),recv_buff_rp(:)
    type(i1p), pointer :: recv_size(:)
    integer(8)                                                :: memor_loc(2)
    logical :: debug
    !type(comm_data_par_basic), pointer       :: comm
    
    
    
    
    !comm => self % parallel_search % comm
    nrank        = self % nrank
    my_rank      = int(self%comm % RANK4,ip)     
  
    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)
    
    
    !We send za and ya perm size
    nullify(recv_buff)
    nullify(recv_size)
    nullify(send_size)
    call memory_alloca(memor_loc,'recv_size',vacal,recv_size,nrank,LBOUN=0_ip)
    call memory_alloca(memor_loc,'send_size',vacal,send_size,2)   
    send_size(:) = 0
    do ineig = 1,self % parallel_search % comm % nneig
      nsend = self % parallel_search % comm % lsend_size(ineig+1) - self % parallel_search % comm % lsend_size(ineig)
      nrecv = self % parallel_search % comm % lrecv_size(ineig+1) - self % parallel_search % comm % lrecv_size(ineig) 
      irank = self % parallel_search % comm % neights(ineig) 
      if( nsend /= 0 .or. nrecv /= 0 ) then  
        send_size(1) = memory_size(interp_matrix_to_send(irank)%yA)
        send_size(2) = memory_size(interp_matrix_to_send(irank)%yA_conc) 
        call memory_alloca(memor_loc,'recv_buff',vacal,recv_buff,2,'DO_NOT_INITIALIZE') 
        call PAR_SEND_RECEIVE(send_size,recv_buff,DOM_I=irank,&
             PAR_COMM_IN=self % parallel_search % comm % PAR_COMM_WORLD )
        !Unpack the buffer
        nullify(recv_size(irank)%l)
        call memory_alloca(memor_loc,'recv_size%l',vacal,recv_size(irank)%l,2)
        recv_size(irank)%l(:) = recv_buff(:)   
        call memory_deallo(memor_loc,'recv_buff',vacal,recv_buff)           
      end if
    end do
    call memory_deallo(memor_loc,'send_size',vacal,send_size)   
    
    print*,"Simon//2.memorloc=",my_rank,memor_loc

    debug = .true.
    if(debug) then
     nullify(send_buff)
     nullify(recv_buff)
     nullify(send_buff_rp)
     nullify(recv_buff_rp)
     !We send za and ya perm
     do ineig = 1,self % parallel_search % comm % nneig
       nsend = self % parallel_search % comm % lsend_size(ineig+1) - self % parallel_search % comm % lsend_size(ineig)
       nrecv = self % parallel_search % comm % lrecv_size(ineig+1) - self % parallel_search % comm % lrecv_size(ineig) 
       irank = self % parallel_search % comm % neights(ineig)
       if(nsend /=0) then
         size_buff1 = memory_size(interp_matrix_to_send(irank)%yA)
         size_buff2 = memory_size(interp_matrix_to_send(irank)%yA_conc)
         call memory_alloca(memor_loc,'send_buff',vacal,send_buff, 3*size_buff1 +size_buff2  ,'DO_NOT_INITIALIZE')   
         call memory_alloca(memor_loc,'send_buff_rp',vacal,send_buff_rp, size_buff1 ,'DO_NOT_INITIALIZE')  
         
         send_buff(1:size_buff1 ) = interp_matrix_to_send(irank)%xA(1:size_buff1)
         send_buff(size_buff1+1: 2*size_buff1 ) = interp_matrix_to_send(irank)%yA(1:size_buff1)
         send_buff(2*size_buff1+1: 3*size_buff1 ) = interp_matrix_to_send(irank)%yA_perm(1:size_buff1)
         send_buff( 3*size_buff1+1: 3*size_buff1 +size_buff2  ) = interp_matrix_to_send(irank)%yA_conc(1:size_buff2)
         send_buff_rp(1:size_buff1 ) = interp_matrix_to_send(irank)%va(1,1,1:size_buff1 )
       end if
       
       if(nrecv /= 0) then 
         size_buff1 = recv_size(irank)%l(1)
         size_buff2 = recv_size(irank)%l(2)
         call memory_alloca(memor_loc,'recv_buff',vacal,recv_buff, 3*size_buff1 +size_buff2,'DO_NOT_INITIALIZE')
         call memory_alloca(memor_loc,'recv_buff_rp',vacal,recv_buff_rp,size_buff1,'DO_NOT_INITIALIZE')
       end if
       
       if( nsend /= 0 .or. nrecv /= 0 ) then  
         call PAR_SEND_RECEIVE(send_buff,recv_buff,DOM_I=irank,&
              PAR_COMM_IN=self % parallel_search % comm % PAR_COMM_WORLD)
         call PAR_SEND_RECEIVE(send_buff_rp,recv_buff_rp,DOM_I=irank,&
              PAR_COMM_IN=self % parallel_search % comm % PAR_COMM_WORLD)
       end if 
  
       
       if(nrecv/=0) then
         size_buff1 = recv_size(irank)%l(1)
         size_buff2 = recv_size(irank)%l(2)
         if(size_buff2>0) then 
           self % interp_matrix(irank) % ndof1 = 1
           self % interp_matrix(irank) % ndof2 = 1
           self % interp_matrix(irank) % nz    = size_buff1
           print*,"Simon//allocating=",my_rank, irank
           call self % interp_matrix(irank) % alloca(MEMORY_COUNTER=memor_loc) 
           call memory_alloca(memor_loc,'self%interp_matrix(irank)%yA_perm',vacal,self%interp_matrix(irank)%yA_perm,size_buff1)
           call memory_alloca(memor_loc,'self%interp_matrix(irank)%yA_conc',vacal,self%interp_matrix(irank)%yA_conc,size_buff2)
           !allocate(self%interp_matrix(irank)%yA_perm(1:size_buff1))
           !allocate(self%interp_matrix(irank)%yA_conc(1:size_buff2))
           self%interp_matrix(irank)%xA(1:size_buff1) = recv_buff(1:size_buff1)
           self%interp_matrix(irank)%yA(1:size_buff1) = recv_buff(size_buff1+1:2*size_buff1)
           self%interp_matrix(irank)%yA_perm(1:size_buff1) = recv_buff(2*size_buff1+1:3*size_buff1)
           self%interp_matrix(irank)%yA_conc(1:size_buff2) = recv_buff(3*size_buff1+1:3*size_buff1+ size_buff2)
           self%interp_matrix(irank)%vA(1,1,1:size_buff1) = recv_buff_rp(1:size_buff1)
           if(associated(recv_buff)) call memory_deallo(memor_loc,'recv_buff',vacal,recv_buff)
           if(associated(recv_buff_rp))call memory_deallo(memor_loc,'recv_buff_rp',vacal,recv_buff_rp)         
        end if
       end if
       if(associated(send_buff)) call memory_deallo(memor_loc,'send_buff',vacal,send_buff) 
       if(associated(send_buff_rp)) call memory_deallo(memor_loc,'send_buff_rp',vacal,send_buff_rp) 
     end do
    end if ! debug   
    
    
    
    
    
    
    if(associated(recv_size)) then
    do ineig =  lbound(recv_size,1),ubound(recv_size,1)
      call memory_deallo(memor_loc,'self%recv_size%l',vacal,recv_size(ineig)%l)
    end do
      call memory_deallo(memor_loc,'recv_size',vacal,recv_size)
    end if
    
    
    
    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER) 
    
    print*,"Simon//debug end communicate matrixs"
    
    
  end subroutine communicate_matrixs
  
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  SSantoso
  !> @date    2020-10-05
  !> @brief   Build Communicator
  !> @details Build communicator to communicate nodal values
  !> 
  !-----------------------------------------------------------------------

  subroutine build_comm(self,MEMORY_COUNTER)

    class(coupling_met),                      intent(inout) :: self                !< Parallel structure
    integer(8),       optional,           intent(inout) :: MEMORY_COUNTER(2)   !< Memory counters
    integer(ip)                                         :: nrank,irank,my_rank
    integer(ip)                                         :: ineig,ii,jj,kk
    integer(8)                                          :: memor_loc(2)
    !
    ! Send/receive results self%my_projection
    ! MY_LIST_RECV(IRANK) % L(II) > 0 ... I'm the owner of point II
    !                             = 0 ... I'm not
    !
    
    print*,"Simon//debug build_comm"
 
    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)

    nrank =self % nrank
    my_rank = int(self%comm% RANK4,ip) 
    self % comm % nneig = nrank
    nullify(self % comm % neights)
    nullify(self % comm % lsend_size  )
    nullify(self % comm % lrecv_size )
    call memory_alloca(memor_loc,'SELF % COMM % NEIGHTS'    ,'def_coupling_method',self % comm % neights,self % comm % nneig)
    call memory_alloca(memor_loc,'SELF % COMM % LSEND_SIZE' ,'def_coupling_method',self % comm % lsend_size    ,self % comm % nneig+1_ip)
    call memory_alloca(memor_loc,'SELF % COMM % LRECV_SIZE' ,'def_coupling_method',self % comm % lrecv_size    ,self % comm % nneig+1_ip)
    
    self % comm % lsend_dim = 0
    self % comm % lrecv_dim = 0

    do ineig = 1,nrank
      irank = ineig - 1
      self % comm % neights(ineig) = irank
      if( associated(self%point_to_send) ) then
        if( associated(self%point_to_send(irank) % l) ) then
          self % comm % lsend_size(ineig) = self % comm % lsend_size(ineig) + memory_size(self%point_to_send(irank)%l)
          self % comm % lsend_dim         = self % comm % lsend_dim + self % comm % lsend_size(ineig)
        end if
      end if
      if( associated(self%interp_matrix) ) then       
        if( associated(self%interp_matrix(irank) % ya_conc) ) then
            self % comm % lrecv_size(ineig) = self % comm % lrecv_size(ineig) + memory_size(self%interp_matrix(irank)%ya_conc)
            self % comm % lrecv_dim         = self % comm % lrecv_dim + self % comm % lrecv_size(ineig)
        end if
      end if
    end do
    
    !
    ! Permutation arrays
    ! lrecv_perm(:)
    ! lsend_perm(:)
    !
    if( associated(self%interp_matrix) ) then
       call memory_alloca(memor_loc,'PAR % COMM % LRECV_PERM','def_search_parall',self % comm % lrecv_perm,self % comm % lrecv_dim)
       jj = 0
       do ineig = 1,nrank
          irank = ineig - 1
          if( associated(self%interp_matrix(irank)%ya_conc) ) then 
           do ii = 1,memory_size(self%interp_matrix(irank)%ya_conc)
              kk = self%interp_matrix(irank)%ya_conc(ii)
              if( kk /= 0 ) then
                 jj = jj + 1
                 self % comm % lrecv_perm(jj) = kk
              end if
           end do
         end if
       end do
    end if
    if( associated(self%point_to_send) ) then
       call memory_alloca(memor_loc,'PAR % COMM % LSEND_PERM','def_search_parall',self % comm % lsend_perm,self % comm % lsend_dim)
       do ii = 1,self % comm % lsend_dim
          self % comm % lsend_perm(ii) = ii
       end do
    end if
    !
    ! LSEND_SIZE and LRECV_SIZE has a linked list
    !
    call number_to_linked_list(self % comm % nneig,self % comm % lsend_size)
    call number_to_linked_list(self % comm % nneig,self % comm % lrecv_size)
    
    
    

    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)
    
    print*,"Simon//debug end  build_comm"

  end subroutine build_comm  
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  SSantoso
  !> @date    2020-10-05
  !> @brief   values_interp version 1
  !> @details values_interp version 1 dimension
  !> 
  !-----------------------------------------------------------------------  
  
  
  
  subroutine values_interp1d(self,xx_in,xx_out,INITIALIZATION,MEMORY_COUNTER)

    class(coupling_met),                    intent(inout) :: self
    real(rp),                      pointer, intent(in)    :: xx_in(:)
    real(rp),                      pointer, intent(inout) :: xx_out(:)
    logical(lg),          optional,         intent(in)    :: INITIALIZATION
    integer(8),           optional,         intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    
    integer(ip)                                           :: ii,nd,nn
    integer(ip)                                           :: nrank,irank,ineig,jj
    integer(ip)                                           :: nsend,nrecv
    integer(ip)                                           :: nd_out,nd_in,my_rank
    integer(ip)                                           :: kk,inode
    real(rp),                      pointer                :: xx_send(:)
    real(rp),                      pointer                :: xx_recv(:)
    logical(lg)                                           :: if_initialization
    real(rp)                                              :: time1,time2
    integer(8)                                            :: memor_loc(2)
    integer(ip)                                          :: index_point,index_recv
    
     
    
    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)
    call cputim(time1)
    
    nd     = max(memory_size(xx_in,1_ip),memory_size(xx_out,1_ip))
    nn     = self % nn
    nd_in  = memory_size(xx_in, 1_ip)
    nd_out = memory_size(xx_out,1_ip)
    my_rank = int(self%comm % RANK4,ip)
    !
    ! Initialize
    !
    if_initialization = optional_argument(.true.,INITIALIZATION)
    if( if_initialization .and. associated(xx_out) ) then
      do ii = 1,memory_size(xx_out,2_ip)
        xx_out(ii) = 0.0_rp
      end do
    end if
    
         
    !if( self % my_interpolation % input_data % mode == INT_PARALLEL ) then
      !
      ! Parallel mode
      !
    nullify(xx_send)
    nullify(xx_recv)
    jj = 0
    nrank  = self % nrank
    if(self %  mode == INT_PARALLEL) then
      jj=0
      do ineig = 1,self % comm % nneig
        nsend = self % comm % lsend_size(ineig+1) - self% comm % lsend_size(ineig)
        nrecv = self % comm % lrecv_size(ineig+1) - self% comm % lrecv_size(ineig)   
        irank = self % comm % neights(ineig)   
        if( nsend /= 0 .or. nrecv /= 0 ) then
          call memory_alloca(memor_loc,'XX_SEND',vacal,xx_send,nsend)
          call memory_alloca(memor_loc,'XX_RECV',vacal,xx_recv,nrecv,'DO_NOT_INITIALIZE')  
          kk=0
          do ii = 1,memory_size(self%point_to_send(irank)%l)
            inode = self%point_to_send(irank) % l(ii)
            xx_send(ii) = xx_in(inode) 
         end do
         call PAR_SEND_RECEIVE(xx_send,xx_recv,DOM_I=irank,PAR_COMM_IN=self  % comm % PAR_COMM_WORLD)
         !Do the computation
          do ii = 1,memory_size(self%interp_matrix(irank)%xA)
            index_point = self%interp_matrix(irank)%xA(ii)
            index_recv = self%interp_matrix(irank)%yA_perm(ii)
            xx_out(index_point) = xx_out(index_point) + self%interp_matrix(irank)%vA(1,1,ii)*xx_recv(index_recv)
          end do
        end if
        call memory_deallo(memor_loc,'XX_SEND',vacal,xx_send)
        call memory_deallo(memor_loc,'XX_RECV',vacal,xx_recv)
      end do 
    else
     do ii = 1,memory_size(self%interp_matrix(0)%xA)
       index_point = self%interp_matrix(0)%xA(ii)
       index_recv = self%interp_matrix(0)%yA(ii)
       xx_out(index_point) = xx_out(index_point) + self%interp_matrix(0)%vA(1,1,ii)*xx_in(index_recv)
     end do
    end if
    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)
    call cputim(time2) ; self % times(3) = self % times(3) + time2-time1 
    
    
    
  end subroutine values_interp1d  
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  SSantoso
  !> @date    2020-10-05
  !> @brief   values_interp version 2
  !> @details values_interp version 2 dimension
  !> 
  !-----------------------------------------------------------------------  
  
  
  
  subroutine values_interp2d(self,xx_in,xx_out,INITIALIZATION,MEMORY_COUNTER)

    class(coupling_met),                   intent(inout) :: self
    real(rp),                      pointer, intent(in)    :: xx_in(:,:)
    real(rp),                      pointer, intent(inout) :: xx_out(:,:)
    logical(lg),          optional,         intent(in)    :: INITIALIZATION
    integer(8),           optional,         intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    !--------------------------------------------------------------------------------------------
    integer(ip)                                           :: ii,nd,nn
    integer(ip)                                           :: nrank,irank,ineig,jj
    integer(ip)                                           :: nsend,nrecv
    integer(ip)                                           :: nn_out,nn_in
    integer(ip)                                           :: nd_out,nd_in,my_rank
    integer(ip)                                           :: kk,inode
    real(rp),                      pointer                :: xx_send(:,:)
    real(rp),                      pointer                :: xx_recv(:,:)
    real(rp),                      pointer                :: xx_values(:,:)
    real(rp),                      pointer                :: x1_send(:)
    real(rp),                      pointer                :: x1_in(:)
    real(rp),                      pointer                :: x1_out(:)
    logical(lg)                                           :: if_initialization
    real(rp)                                              :: time1,time2
    integer(ip)                                           :: index_point,index_recv
    integer(8)                                            :: memor_loc(2)  
       
    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)
    nullify(x1_in)
    nullify(x1_out)
    nullify(xx_values)
    call cputim(time1)
    
    irank = -1_ip

    nd     = max(memory_size(xx_in,1_ip),memory_size(xx_out,1_ip))
    nn     = self % nn
    nn_in  = memory_size(xx_in, 2_ip)
    nn_out = memory_size(xx_out,2_ip)
    nd_in  = memory_size(xx_in, 1_ip)
    nd_out = memory_size(xx_out,1_ip)
    my_rank = int(self%comm % RANK4,ip)
    if( nd_in /= 0 .and. nd_out /= 0 .and. nd_in /= nd_out ) &
         call runend('values_interp2d: WRONG DIMENSION 1 FOR OUTPUT ARRAY')
    !
    ! Initialize
    !
    if_initialization = optional_argument(.true.,INITIALIZATION)
    if( if_initialization .and. associated(xx_out) ) then
      do ii = 1,memory_size(xx_out,2_ip)
        xx_out(:,ii) = 0.0_rp
      end do
    end if
    nullify(xx_send)
    nullify(xx_recv)
    nullify(x1_send)
    jj = 0
    nrank  = self % nrank
    if(self % my_interpolation% input_data % mode == INT_PARALLEL) then
      !
      ! Parallel mode
      !
      jj=0
      do ineig = 1,self % comm % nneig
        nsend = self % comm % lsend_size(ineig+1) - self% comm % lsend_size(ineig)
        nrecv = self % comm % lrecv_size(ineig+1) - self% comm % lrecv_size(ineig)   
        irank = self % comm % neights(ineig)      
        if( nsend /= 0 .or. nrecv /= 0 ) then
          call memory_alloca(memor_loc,'XX_SEND',vacal,xx_send,nd,nsend)
          call memory_alloca(memor_loc,'XX_RECV',vacal,xx_recv,nd,nrecv,'DO_NOT_INITIALIZE')  
          kk=0
          do ii = 1,memory_size(self%point_to_send(irank)%l)
            inode = self%point_to_send(irank) % l(ii)
            xx_send(1:nd,ii) = xx_in(1:nd,inode) 
          end do
          call PAR_SEND_RECEIVE(xx_send,xx_recv,DOM_I=irank,PAR_COMM_IN=self  % comm % PAR_COMM_WORLD)
          !Do the computation
          do ii = 1,memory_size(self%interp_matrix(irank)%xA)
            index_point = self%interp_matrix(irank)%xA(ii)
            index_recv = self%interp_matrix(irank)%yA_perm(ii)
            xx_out(1:nd,index_point) = xx_out(1:nd,index_point) + self%interp_matrix(irank)%vA(1,1,ii)*xx_recv(1:nd,index_recv)
          end do
        end if
        call memory_deallo(memor_loc,'XX_SEND',vacal,xx_send)
        call memory_deallo(memor_loc,'XX_RECV',vacal,xx_recv)
      end do 
    else
      !TODO: irank may not be initialized
      do ii = 1,memory_size(self%my_interpolation%matrix(irank)%xA)
        index_point = self%interp_matrix(0)%xA(ii)
        index_recv = self%interp_matrix(0)%yA(ii)
        xx_out(1:nd,index_point) = xx_out(1:nd,index_point) + self%interp_matrix(0)%vA(1,1,ii)*xx_in(1:nd,index_recv)
      end do
    end if
    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)
    call cputim(time2) ; self % times(3) = self % times(3) + time2-time1 
    
    
  end subroutine values_interp2d  
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  SSantoso
  !> @date    2021-01-18
  !> @brief   Compute value
  !> @details Compute value
  !> 
  !-----------------------------------------------------------------------

  subroutine values_proj(self,xx_in,xx_out,INITIALIZATION,MEMORY_COUNTER)

   class(coupling_met),                    intent(inout) :: self
   real(rp),                      pointer, intent(in)    :: xx_in(:)
   real(rp),                      pointer, intent(inout) :: xx_out(:)
   logical(lg),          optional,         intent(in)    :: INITIALIZATION
   integer(8),           optional,         intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
   integer(ip)                                           :: ii,ipoin
   integer(8)                                            :: memor_loc(2)
   real(rp),                      pointer                :: xx_ng(:)
   real(rp),                      pointer                :: xx_nn(:)
   integer(ip)                                           :: my_rank
   class(iterative_solver), pointer   :: sol
   logical :: debug



    print*,"Simon//debug values proj"


    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)
    nullify(xx_ng,xx_nn)
    if(self%my_projection%ng>0) call memory_alloca(memor_loc,'XX_NG',vacal,xx_ng,self %my_projection% ng)
    if(self%my_projection%nn>0) call memory_alloca(memor_loc,'XX_NN',vacal,xx_nn,self %my_projection% nn)
    my_rank = int(self%comm % RANK4,ip)
    !
    ! Initialization
    !
    if( optional_argument(.true.,INITIALIZATION) ) then
       do ii = 1,self %my_projection% nn
          ipoin     = self % my_projection% invp(ii)
          xx_nn(ii) = xx_out(ipoin) 
       end do
    end if
    !
    ! Source nodes (global) => Target Gauss-points
    !
    call self % values_interp(xx_in,xx_ng,INITIALIZATION,MEMORY_COUNTER=memor_loc) 
    !
    ! Target Gauss points => Target nodes (local)
    !
    call self % mass_matrix % mv(xx_ng,xx_nn,INITIALIZATION=INITIALIZATION)

    debug = .false.
     if(debug) then
     ! Use of CG to solve MUt = us
     allocate(cg :: sol)
     call solver_init(sol)
     allocate(sol % comm)
     call sol % comm % init()
     sol % comm % RANK4          = int(self % comm % RANK4,4)
     sol % comm % SIZE4          = int(self % comm % size4,4)
     !sol % comm % PAR_COMM_WORLD = MPI_COMM_WORLD 
     sol % comm % nneig         = 0
     print*,"Simon//",my_rank,sol % comm % nneig
     sol % comm % bound_dim     = 0
     print*,"Simon//",my_rank,self%comm%bound_dim
     call sol % comm % alloca()
     !sol%comm = self % comm

     !sol % comm % bound_perm(1) = 0
     !sol % comm % bound_size(1) = 0
     !sol % comm % bound_size(2) = 0
  
     call solver_dim (sol,self%mass_matrix)
     sol % input % relax     = 0.5_rp
     sol % input % kfl_symm  = 1
     sol % input % kfl_preco = SOL_SOLVER_DIAGONAL
     sol % input % nkryd     = 4
     sol % input % solco     = 1.0e-8_rp
     sol % input % miter     = 1000
     sol % input % kfl_exres = 1
     sol % input % kfl_cvgso = 1
     sol % input % lun_cvgso = 90+ii
     call preconditioner_alloca(sol)
     if( my_rank == 0 ) print*,'solver setup'
     call solver_setup (sol,self%mass_matrix)
     call sol % parallel_exchange(xx_nn)
     if( my_rank == 0 ) print*,'set preconditioner'
     call preconditioner_setup(sol,self%mass_matrix)
     !xx_nn = 1.0
     print*,'xx_ng = ',my_rank,xx_ng
     call sol % solve (xx_ng,xx_nn,self%mass_matrix)
     if( my_rank == 0 ) then
       print*,'resip_init = ',sol % output % resip_init
       print*,'resid_init = ',sol % output % resid_init
       print*,'resip      = ',sol % output % resip_final
       print*,'resid      = ',sol % output % resid_final
       print*,'iters      = ',sol % output % iters
       print*,'||b||      = ',sol % output % bnorm
       print*,'||L^-1b||  = ',sol % output % bnorp
     end if
    end if

    !
    ! Permute Target nodes (local) => Target nodes (global)
    !
    do ii = 1,self % my_projection%nn
       ipoin         = self %my_projection %invp(ii)
       xx_out(ipoin) = xx_nn(ii)
    end do
    !
    ! Parallel exchange
    !
    if( associated(self % my_projection % comm) ) call par_interface_exchange(xx_out,self %my_projection% comm)


    !call sol % deallo()
    !deallocate(sol)

    call memory_deallo(memor_loc,'XX_NN',vacal,xx_nn)
    call memory_deallo(memor_loc,'XX_NG',vacal,xx_ng)
    

    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)


    print*,"Simon//debug end values proj"



  end subroutine values_proj
  
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  SSantoso
  !> @date    2021-01-18
  !> @brief   Compute value
  !> @details Compute value
  !> 
  !-----------------------------------------------------------------------

  subroutine transpose_matrix(self,matrix_in,matrix_out,MEMORY_COUNTER)

    class(coupling_met),           intent(inout) :: self
    type(mat_coo),                 intent(inout)    :: matrix_in
    type(mat_coo),                      intent(inout) :: matrix_out
    integer(8),           optional,         intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(8)                                            :: memor_loc(2)

    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)
    call matrix_out% init()
    matrix_out % ndof1 = 1
    matrix_out % ndof2 = 1
    matrix_out % nz    = matrix_in%nz
    
    matrix_out% nrows = matrix_in % ncols
    matrix_out% ncols = matrix_in % nrows
    !matrix_out% alloca(MEMORY_COUNTER=memor_loc)
    call matrix_out % alloca(MEMORY_COUNTER=memor_loc)
    matrix_out%xA(:)=matrix_in%yA(:)
    matrix_out%yA(:)=matrix_in%xA(:)
    matrix_out%vA(1,1,:)=self%my_projection%matrix%vA(1,1,:)    
    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)

  end subroutine transpose_matrix  
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-05-19
  !> @brief   Number to linked list
  !> @details Number to linked list
  !> 
  !-----------------------------------------------------------------------

  pure subroutine number_to_linked_list(nn,ia)
    
    integer(ip), intent(in)    :: nn
    integer(ip), intent(inout) :: ia(nn+1)
    integer(ip)                :: ii,kk,ll

    kk    = ia(1)
    ia(1) = 1 
    do ii = 2,nn+1
       ll     = ia(ii)
       ia(ii) = ia(ii-1) + kk
       kk     = ll
    end do

  end subroutine number_to_linked_list
  
end module def_coupling_method
