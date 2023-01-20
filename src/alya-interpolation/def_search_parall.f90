!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module def_search_parall

  use def_kintyp_basic,                  only : ip,rp,lg,i1p,r1p,r2p
  use def_kintyp_comm,                   only : comm_data_par_basic
  use mod_memory_basic,                  only : memory_alloca
  use mod_memory_basic,                  only : memory_deallo
  use mod_memory_basic,                  only : memory_size
  use mod_memory_basic,                  only : memory_resize
  use mod_memory_tools,                  only : memory_counter_ini
  use mod_memory_tools,                  only : memory_counter_end
  use mod_optional_argument,             only : optional_argument
  use mod_maths_geometry,                only : maths_point_in_box
  use mod_communications_tools,          only : PAR_COMM_RANK_AND_SIZE
  use mod_communications_global,         only : PAR_ALLGATHER
  use mod_communications_global,         only : PAR_ALLTOALL
  use mod_communications_point_to_point, only : PAR_SEND_RECEIVE
  use mod_communications_point_to_point, only : PAR_SEND_RECEIVE_RP_2
  use def_search_method,                 only : search_method
  use def_search_method,                 only : SEARCH_BOUNDING_BOXES
  use mod_std
  use def_mpi
#include "def_mpi.inc"
  implicit none 
  
  private

  character(17), parameter :: vacal = 'def_search_parall'
  
  type :: search_parall
     class(search_method), pointer   :: search_method_par       ! Search strategy
     type(comm_data_par_basic)       :: comm                    ! Communication arrays
     integer(ip)                     :: nd                      ! Dimension
     integer(8)                      :: memor(2)                ! Memory counter
   contains
     procedure,               pass   :: init                    ! Initialize all
     procedure,               pass   :: deallo                  ! Deallocate all
     procedure,               pass   :: list_rank               ! List of ranks
     procedure,               pass   :: send_recv_points_1      ! Send/receive points to/from neighbors
     procedure,               pass   :: send_recv_points_2      ! Send/receive points to/from neighbors
     procedure,               pass   :: send_recv_dista         ! Send/receive distance
     procedure,               pass   :: send_list_recv_points_1 ! Send/recv points
     procedure,               pass   :: send_list_recv_points_2 ! Send/recv points
     procedure,               pass   :: input                   ! Input data
     procedure,               pass   :: communication           ! Set communication strategy
     procedure,               pass   :: dista_comm_1            ! Set communication strategy
     procedure,               pass   :: dista_comm_2            ! Set communication strategy

     generic                         :: send_recv_points      => &
          &                             send_recv_points_1,      &
          &                             send_recv_points_2
     generic                         :: send_list_recv_points => &
          &                             send_list_recv_points_1, &
          &                             send_list_recv_points_2
     generic                         :: dista_comm            => &
          &                             dista_comm_1,            &
          &                             dista_comm_2
  end type search_parall

  public :: search_parall
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   InitializationPAR
  !> @details Initialization of the parallel strategy
  !> 
  !-----------------------------------------------------------------------
  
  subroutine init(par)

    class(search_parall), intent(inout) :: par
    
    nullify(par % search_method_par)
    par % nd    = 0_ip
    par % memor = 0_8

    call par % comm % init()
    
  end subroutine init
 
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   InitializationPAR
  !> @details Initialization of the parallel strategy
  !> 
  !-----------------------------------------------------------------------
   
  subroutine deallo(par,MEMORY_COUNTER)

    class(search_parall),         intent(inout) :: par
    integer(8),         optional, intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
      
    call par % comm % deallo(MEMORY_COUNTER)
    
  end subroutine deallo
 
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   Parallel search
  !> @details Give the list of ranks containing each point
  !>
  !>          LIST_SEND(IRANK) % L(:) = List of nodes contained by rank IRANK
  !> 
  !-----------------------------------------------------------------------
  
  subroutine list_rank(par,xx,list_send,METHOD,MASK,EXCLUDE_MYSELF,MEMORY_COUNTER)

    class(search_parall),                intent(inout) :: par               !< Parallel structure
    real(rp),                   pointer, intent(in)    :: xx(:,:)           !< List of coordinates
    type(i1p),                  pointer, intent(inout) :: list_send(:)        !< List of ranks containing each point
    integer(ip),      optional,          intent(in)    :: METHOD                !< Method for candidates
    logical(lg),      optional, pointer, intent(in)    :: MASK(:)
    logical(lg),      optional,          intent(in)    :: EXCLUDE_MYSELF
    integer(8),       optional,          intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(ip)                                        :: nrank,ii,kk
    integer(ip)                                        :: irank,nn,ndime
    integer(ip)                                        :: my_rank
    integer(8)                                         :: memor_loc(2)
    integer(ip),                pointer                :: list_tmp(:)
    type(i1p),                  pointer                :: node_to_rank(:)
    logical(lg)                                        :: if_exclude_myself
    
    nullify(list_send)
    nullify(list_tmp)
    nullify(node_to_rank)

    call memory_counter_ini(memor_loc,par % memor,MEMORY_COUNTER)

    my_rank           = int(par % comm % RANK4,ip)
    nrank             = int(par % comm % size4,ip)
    ndime             = par % nd 
    nn                = memory_size(xx,2_ip)
    if_exclude_myself = optional_argument(.false.,EXCLUDE_MYSELF)
    
    call memory_alloca(memor_loc,'LIST_SEND','list_rank',list_send,nrank,LBOUN=0_ip)
    
    if( nn > 0 ) then

       call memory_alloca(memor_loc,'LIST_TMP'    ,'list_rank',list_tmp,nrank,LBOUN=0_ip)
       call memory_alloca(memor_loc,'NODE_TO_RANK','list_rank',node_to_rank,nn)
       call par % search_method_par % candidate(xx,node_to_rank,METHOD=METHOD,MASK=MASK,MEMORY_COUNTER=memor_loc)

!!$       block
!!$         use def_master
!!$         integer(ip) :: irank
!!$         if( par % search_method_par % name == 'BIN PAR' ) then
!!$            do ii = 1,nn
!!$               if( associated(node_to_rank(ii) % l  )) then
!!$                  print*,'size node=',kfl_paral,memory_size(node_to_rank(ii) % l)
!!$                  
!!$                  print*,'send node=',ii,'from ',kfl_paral,'->',node_to_rank(ii) % l
!!$               end if
!!$               !if( nn_send(irank) > 0 ) print*,'send near=',kfl_paral,'->',irank,': ',nn_send(irank)!,memory_size(dista_recv(irank) % a)
!!$               !if( nn_recv(irank) > 0 ) print*,'recv nrea=',kfl_paral,'<-',irank,': ',nn_recv(irank)!,memory_size(dista_send(irank) % a)
!!$            end do
!!$            !  end do
!!$            !call runend('O.K.!')
!!$         end if
!!$       end block
       !
       ! Exclude myself
       !
       if( if_exclude_myself ) then
          do ii = 1,nn
             do kk = 1,memory_size(node_to_rank(ii) % l)
                irank = node_to_rank(ii) % l(kk)
                if( irank == my_rank ) node_to_rank(ii) % l(kk) = -1
             end do
          end do
       end if
       !
       ! Count nu,ber of points to send
       !
       do ii = 1,nn
          do kk = 1,memory_size(node_to_rank(ii) % l)
             irank = node_to_rank(ii) % l(kk)
             if( irank /= -1 ) list_tmp(irank) = list_tmp(irank) + 1
          end do
       end do
       !
       ! LIST_SEND % L: Fill in send list
       !
       do irank = 0,nrank - 1
          call memory_alloca(memor_loc,'LIST_SEND % L','list_rank',list_send(irank) % l,list_tmp(irank))
       end do

       if( nrank > 0 ) then
          do ii = 0,nrank-1
             list_tmp(ii) = 0
          end do
       end if
       
       do ii = 1,nn
          do kk = 1,memory_size(node_to_rank(ii) % l)
             irank = node_to_rank(ii) % l(kk)
             if( irank /= -1 ) then
                list_tmp(irank) = list_tmp(irank) + 1
                list_send(irank) % l(list_tmp(irank))=ii
             end if
          end do
       end do

       call memory_deallo(memor_loc,'NODE_TO_RANK','list_rank',node_to_rank)
       call memory_deallo(memor_loc,'LIST_TMP'    ,'list_rank',list_tmp)
    end if
     
    call memory_counter_end(memor_loc,par % memor,MEMORY_COUNTER)

  end subroutine list_rank

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-05
  !> @brief   Send and receive points
  !> @details Send points to selected ranks according to MY_LIST
  !>          and receive points
  !>
  !>          XX_SEND(IRANK) % L(:) = List of points to send to rank IRANK
  !>          XX_RECV(IRANK) % L(:) = List of points received from rank IRANK
  !> 
  !-----------------------------------------------------------------------

  subroutine send_recv_points_1(par,xx,nn_send,nn_recv,xx_recv,my_list,MEMORY_COUNTER)

    class(search_parall),                 intent(inout) :: par               !< Parallel structure
    real(rp),                    pointer, intent(in)    :: xx(:,:)           !< List of coordinates
    integer(ip),                 pointer, intent(inout) :: nn_send(:)        !< Number of points to send
    integer(ip),                 pointer, intent(inout) :: nn_recv(:)        !< Number of points to receive
    type(r2p),                   pointer, intent(inout) :: xx_recv(:)        !< List of received points for each rank
    type(i1p),                   pointer, intent(inout) :: my_list(:)        !< List of ranks for each point
    integer(8),       optional,           intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    real(rp),                    pointer                :: xx_send(:,:)
    integer(8)                                          :: memor_loc(2)
    integer(ip)                                         :: irank,ii,kk
    integer(ip)                                         :: ndime,nrank,nn
    integer(ip)                                         :: my_rank
    MY_MPI_COMM                                         :: comm4 

    nullify(nn_send)
    nullify(nn_recv)
    nullify(xx_send)
    nullify(xx_recv)

    call memory_counter_ini(memor_loc,par % memor,MEMORY_COUNTER)
    
    ndime     = par % nd 
    nn        = memory_size(xx,2_ip)    
    nrank     = int(par % comm % size4,ip)
    comm4     = par % comm % PAR_COMM_WORLD
    my_rank   = int(par % comm % RANK4,ip)

    !call PAR_MAX(ndime,COMM4,INCLUDE_ROOT=.true.)

    call memory_alloca(memor_loc,'NN_SEND','send_point',nn_send,nrank,LBOUN=0_ip)
    call memory_alloca(memor_loc,'NN_RECV','send_point',nn_recv,nrank,LBOUN=0_ip)

    if( associated(my_list) ) then
       do irank = 0,nrank - 1
          nn_send(irank) = memory_size(my_list(irank) % l)
       end do
    end if

    call PAR_ALLTOALL(1_ip,1_ip,nn_send,nn_recv,PAR_COMM_IN=comm4)
    !
    ! XX_RECV(IRANK) % L(:,:)
    !
    if( .not. associated(xx_recv) ) & 
         call memory_alloca(memor_loc,'XX_RECV','send_point',xx_recv,nrank,LBOUN=0_ip)
    do irank = 0,nrank - 1
       call memory_alloca(memor_loc,'XX_RECV % A','send_point',xx_recv(irank) % a,ndime,nn_recv(irank))
       call memory_alloca(memor_loc,'XX_SEND'    ,'send_point',xx_send,ndime,nn_send(irank))
       do kk = 1,nn_send(irank)
          ii = my_list(irank) % l(kk)
          xx_send(1:ndime,kk) = xx(1:ndime,ii) 
       end do

       if( irank == my_rank ) then
          do kk = 1,nn_send(irank)
             xx_recv(irank) % a(1:ndime,kk) = xx_send(1:ndime,kk)
          end do
       else
          call PAR_SEND_RECEIVE(xx_send,xx_recv(irank) % a,DOM_I=irank,PAR_COMM_IN=comm4)
       end if
       call memory_deallo(memor_loc,'XX_SEND','send_point',xx_send)
    end do
    
    call memory_counter_end(memor_loc,par % memor,MEMORY_COUNTER)

!!$    block
!!$      use def_master
!!$      integer(ip) :: irank
!!$      do irank = 0,nrank - 1
!!$         if( nn_send(irank) > 0 ) print*,'send A=',kfl_paral,'->',irank,': ',nn_send(irank)
!!$         if( nn_recv(irank) > 0 ) print*,'recv A=',kfl_paral,'<-',irank,': ',nn_recv(irank),ndime,memory_size(xx_recv(irank) % a)
!!$      end do
!!$      call runend('O.K.!')
!!$    end block
    
  end subroutine send_recv_points_1
  
  subroutine send_recv_points_2(par,xx,ll,nn_send,nn_recv,xx_recv,ll_recv,my_list,MEMORY_COUNTER)

    class(search_parall),                 intent(inout) :: par               !< Parallel structure
    real(rp),                    pointer, intent(in)    :: xx(:,:)           !< List of coordinates
    integer(ip),                 pointer, intent(in)    :: ll(:)           !< List of coordinates
    integer(ip),                 pointer, intent(inout) :: nn_send(:)        !< Number of points to send
    integer(ip),                 pointer, intent(inout) :: nn_recv(:)        !< Number of points to receive
    type(r2p),                   pointer, intent(inout) :: xx_recv(:)        !< List of received points for each rank
    type(i1p),                   pointer, intent(inout) :: ll_recv(:)        !< List of received points for each rank
    type(i1p),                   pointer, intent(inout) :: my_list(:)        !< List of ranks for each point
    integer(8),       optional,           intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    real(rp),                    pointer                :: xx_send(:,:)
    integer(ip),                 pointer                :: ll_send(:)
    integer(8)                                          :: memor_loc(2)
    integer(ip)                                         :: irank,ii,kk
    integer(ip)                                         :: ndime,nrank,nn
    integer(ip)                                         :: my_rank
    MY_MPI_COMM                                         :: comm4 

    nullify(nn_send)
    nullify(nn_recv)
    nullify(xx_send)
    nullify(xx_recv)
    nullify(ll_send)
    nullify(ll_recv)

    call memory_counter_ini(memor_loc,par % memor,MEMORY_COUNTER)
    
    ndime     = par % nd 
    nn        = memory_size(xx,2_ip)    
    nrank     = int(par % comm % size4,ip)
    comm4     = par % comm % PAR_COMM_WORLD
    my_rank   = int(par % comm % RANK4,ip)

    !call PAR_MAX(ndime,COMM4,INCLUDE_ROOT=.true.)

    call memory_alloca(memor_loc,'NN_SEND','send_point',nn_send,nrank,LBOUN=0_ip)
    call memory_alloca(memor_loc,'NN_RECV','send_point',nn_recv,nrank,LBOUN=0_ip)

    if( associated(my_list) ) then
       do irank = 0,nrank - 1
          nn_send(irank) = memory_size(my_list(irank) % l)
       end do
    end if

    call PAR_ALLTOALL(1_ip,1_ip,nn_send,nn_recv,PAR_COMM_IN=comm4)
    !
    ! XX_RECV(IRANK) % A(:,:)
    ! LL_RECV(IRANK) % L(:)
    !
    if( .not. associated(xx_recv) ) & 
         call memory_alloca(memor_loc,'XX_RECV','send_point',xx_recv,nrank,LBOUN=0_ip)
    if( .not. associated(ll_recv) ) & 
         call memory_alloca(memor_loc,'LL_RECV','send_point',ll_recv,nrank,LBOUN=0_ip)
    do irank = 0,nrank - 1
       call memory_alloca(memor_loc,'XX_RECV % A','send_point',xx_recv(irank) % a,ndime,nn_recv(irank))
       call memory_alloca(memor_loc,'XX_SEND'    ,'send_point',xx_send,ndime,nn_send(irank))
       call memory_alloca(memor_loc,'LL_RECV % L','send_point',ll_recv(irank) % l,nn_recv(irank))
       call memory_alloca(memor_loc,'LL_SEND'    ,'send_point',ll_send,nn_send(irank))
       do kk = 1,nn_send(irank)
          ii                  = my_list(irank) % l(kk)
          xx_send(1:ndime,kk) = xx(1:ndime,ii) 
          ll_send(kk)         = ll(ii) 
       end do

       if( irank == my_rank ) then
          do kk = 1,nn_send(irank)
             xx_recv(irank) % a(1:ndime,kk) = xx_send(1:ndime,kk)
             ll_recv(irank) % l(kk)         = ll_send(kk)
          end do
       else
          call PAR_SEND_RECEIVE(xx_send,xx_recv(irank) % a,DOM_I=irank,PAR_COMM_IN=comm4)
          call PAR_SEND_RECEIVE(ll_send,ll_recv(irank) % l,DOM_I=irank,PAR_COMM_IN=comm4)
       end if
       call memory_deallo(memor_loc,'XX_SEND','send_point',xx_send)
       call memory_deallo(memor_loc,'ll_SEND','send_point',ll_send)
    end do
    
    call memory_counter_end(memor_loc,par % memor,MEMORY_COUNTER)

!!$    block
!!$      use def_master
!!$      integer(ip) :: irank
!!$      do irank = 0,nrank - 1
!!$         if( nn_send(irank) > 0 ) print*,'send A=',kfl_paral,'->',irank,': ',nn_send(irank)
!!$         if( nn_recv(irank) > 0 ) print*,'recv A=',kfl_paral,'<-',irank,': ',nn_recv(irank),ndime,memory_size(xx_recv(irank) % a)
!!$      end do
!!$      call runend('O.K.!')
!!$    end block
    
  end subroutine send_recv_points_2
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-05
  !> @brief   Take decision
  !> @details Receive decision
  !> 
  !-----------------------------------------------------------------------

  subroutine communication(par,nn_send,nn_recv,my_list_send,my_list_recv,MEMORY_COUNTER)

    class(search_parall),                 intent(inout) :: par                 !< Parallel structure
    integer(ip),                 pointer, intent(in)    :: nn_send(:)
    integer(ip),                 pointer, intent(in)    :: nn_recv(:)
    type(i1p),                   pointer, intent(in)    :: my_list_send(:)     !< List of coordinates
    type(i1p),                   pointer, intent(inout) :: my_list_recv(:)     !< List of coordinates
    integer(8),       optional,           intent(inout) :: MEMORY_COUNTER(2)   !< Memory counters
    integer(ip)                                         :: nrank,irank
    integer(ip)                                         :: ineig,ii,jj,kk
    integer(8)                                          :: memor_loc(2)
    integer(ip),                   pointer              :: lrecv_perm_tmp(:)
    integer(ip),                   pointer              :: lrecv_size_tmp(:)
    logical(lg)                                         :: if_already
    !
    ! Send/receive results
    ! MY_LIST_RECV(IRANK) % L(II) > 0 ... I'm the owner of point II
    !                             = 0 ... I'm not
    !
    call memory_counter_ini(memor_loc,par % memor,MEMORY_COUNTER)

    nrank = int(par % comm % size4,ip)

    nullify(lrecv_perm_tmp)
    nullify(lrecv_size_tmp)
    
    if( .not. associated(my_list_recv) ) &
         call memory_alloca(memor_loc,'MY_LIST_RECV','def_search_parall',my_list_recv,nrank,LBOUN=0_ip)

    do irank = 0,nrank-1
       if( .not. associated(my_list_recv(irank) % l) ) &
            call memory_alloca(memor_loc,'MY_LIST_RECV % L','def_search_parall',my_list_recv(irank) % l,nn_recv(irank))          
    end do
    do irank = 0,nrank-1
       if( nn_recv(irank) > 0 .or. nn_send(irank) > 0 ) then
          call PAR_SEND_RECEIVE(my_list_send(irank) % l,my_list_recv(irank) % l,DOM_I=irank,PAR_COMM_IN=par % comm % PAR_COMM_WORLD)
       end if
    end do
    !
    ! Set communication arrays
    ! lsend_size(:)
    ! lsend_dim
    ! lrecv_size(:) 
    ! lrecv_dim
    !
    par % comm % nneig = nrank
    if_already         = associated(par % comm % neights)

    call memory_alloca(memor_loc,'PAR % COMM % NEIGHTS'    ,'def_search_parall',par % comm % neights,par % comm % nneig)
    call memory_alloca(memor_loc,'PAR % COMM % LSEND_SIZE' ,'def_search_parall',par % comm % lsend_size    ,par % comm % nneig+1_ip)
    call memory_alloca(memor_loc,'PAR % COMM % LRECV_SIZE' ,'def_search_parall',par % comm % lrecv_size    ,par % comm % nneig+1_ip)
    
    par % comm % lsend_dim = 0
    par % comm % lrecv_dim = 0
    
    do ineig = 1,nrank
       irank = ineig - 1
       par % comm % neights(ineig) = irank
       if( associated(my_list_recv) ) then
          if( associated(my_list_recv(irank) % l) ) then
             par % comm % lsend_size(ineig) = par % comm % lsend_size(ineig) + count(my_list_recv(irank) % l/=0,KIND=ip)
             par % comm % lsend_dim         = par % comm % lsend_dim + par % comm % lsend_size(ineig)
          end if
       end if
       if( associated(my_list_send) ) then       
          if( associated(my_list_send(irank) % l) ) then
             par % comm % lrecv_size(ineig) = par % comm % lrecv_size(ineig) + count(my_list_send(irank) % l/=0,KIND=ip)
             par % comm % lrecv_dim         = par % comm % lrecv_dim + par % comm % lrecv_size(ineig)
          end if
       end if
    end do
    !
    ! Permutation arrays
    ! lrecv_perm(:)
    ! lsend_perm(:)
    !
    if( associated(my_list_send) ) then
       call memory_alloca(memor_loc,'PAR % COMM % LRECV_PERM','def_search_parall',par % comm % lrecv_perm,par % comm % lrecv_dim)
       jj = 0
       do ineig = 1,nrank
          irank = ineig - 1
          do ii = 1,memory_size(my_list_send(irank) % l)
             kk = my_list_send(irank) % l(ii)
             if( kk /= 0 ) then
                jj = jj + 1
                par % comm % lrecv_perm(jj) = kk
             end if
          end do
       end do
    end if
    if( associated(my_list_recv) ) then
       call memory_alloca(memor_loc,'PAR % COMM % LSEND_PERM','def_search_parall',par % comm % lsend_perm,par % comm % lsend_dim)
       do ii = 1,par % comm % lsend_dim
          par % comm % lsend_perm(ii) = ii
       end do
    end if
    !
    ! LSEND_SIZE and LRECV_SIZE has a linked list
    !
    call number_to_linked_list(par % comm % nneig,par % comm % lsend_size)
    call number_to_linked_list(par % comm % nneig,par % comm % lrecv_size)
    
    call memory_counter_end(memor_loc,par % memor,MEMORY_COUNTER)

  end subroutine communication

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-05
  !> @brief   Take decision
  !> @details Receive decision
  !> 
  !-----------------------------------------------------------------------

  subroutine send_recv_dista(par,nn_send,nn_recv,dista_send,dista_recv,MEMORY_COUNTER)
    
    class(search_parall),                 intent(inout) :: par               !< Parallel structure
    integer(ip),                 pointer, intent(in)    :: nn_send(:)
    integer(ip),                 pointer, intent(in)    :: nn_recv(:)
    type(r1p),                   pointer, intent(in)    :: dista_send(:)     !< List of coordinates
    type(r1p),                   pointer, intent(inout) :: dista_recv(:)     !< List of coordinates
    integer(8),       optional,           intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(ip)                                         :: nrank,irank
    integer(8)                                          :: memor_loc(2)

    call memory_counter_ini(memor_loc,par % memor,MEMORY_COUNTER)

    nrank = int(par % comm % size4,ip)
    
    if( .not. associated(dista_recv) ) &
         call memory_alloca(memor_loc,'DISTA_RECV','def_search_parall',dista_recv,nrank,LBOUN=0_ip)

    do irank = 0,nrank-1
       if( .not. associated(dista_recv(irank) % a) ) &
            call memory_alloca(memor_loc,'DISTA_RECV % A','def_search_parall',dista_recv(irank) % a,nn_send(irank))          
    end do

!!$    block
!!$      use def_master
!!$      integer(ip) :: irank
!!$      do irank = 0,nrank - 1
!!$         if( nn_send(irank) > 0 ) print*,'send=',kfl_paral,'->',irank,': ',nn_send(irank),memory_size(dista_recv(irank) % a)
!!$         if( nn_recv(irank) > 0 ) print*,'recv=',kfl_paral,'<-',irank,': ',nn_recv(irank),memory_size(dista_send(irank) % a)
!!$      end do
!!$    end block
    
    do irank = 0,nrank-1       
       if( nn_recv(irank) > 0 .or. nn_send(irank) > 0 ) then
          call PAR_SEND_RECEIVE(dista_send(irank) % a,dista_recv(irank) % a,DOM_I=irank,PAR_COMM_IN=par % comm % PAR_COMM_WORLD)
       end if
    end do

    call memory_counter_end(memor_loc,par % memor,MEMORY_COUNTER)

  end subroutine send_recv_dista

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-05
  !> @brief   Set
  !> @details Set communicator, rank and size
  !> 
  !-----------------------------------------------------------------------

  subroutine input(par,comm)

    class(search_parall),            intent(inout) :: par
    MY_MPI_COMM,           optional, intent(in)    :: comm

    if( present(comm) ) then
       
       par % comm % PAR_COMM_WORLD = COMM
       call PAR_COMM_RANK_AND_SIZE(COMM,par % comm % rank4,par % comm % size4)

    end if

  end subroutine input

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-17
  !> @brief   Get send list and receive points
  !> @details Compute:
  !>
  !>          LISTA_SEND(IRANK) % L(:) ... List of points contained by rank IRANK
  !>          XX_RECV(IRANK)   % L(:) ... List of points received from rank IRANK
  !>
  !-----------------------------------------------------------------------

  subroutine send_list_recv_points_1(par,xx,nn_send,nn_recv,lista_send,xx_recv,METHOD,MASK,MEMORY_COUNTER)

    class(search_parall),                    intent(inout) :: par
    real(rp),                       pointer, intent(in)    :: xx(:,:)
    integer(ip),                    pointer, intent(inout) :: nn_send(:)
    integer(ip),                    pointer, intent(inout) :: nn_recv(:)
    type(i1p),                      pointer, intent(inout) :: lista_send(:)  
    type(r2p),                      pointer, intent(inout) :: xx_recv(:)
    integer(ip),      optional,              intent(in)    :: METHOD                !< Method for candidates
    logical(lg),      optional,     pointer, intent(in)    :: MASK(:)
    integer(8),       optional,              intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(8)                                             :: memor_loc(2)

    nullify(nn_send)
    nullify(nn_recv)
    nullify(lista_send)
    nullify(xx_recv)

    call memory_counter_ini(memor_loc,par % memor,MEMORY_COUNTER)

    call par % list_rank       (xx,lista_send,METHOD=METHOD,MASK=MASK,MEMORY_COUNTER=memor_loc)    
    call par % send_recv_points(xx,nn_send,nn_recv,xx_recv,lista_send,MEMORY_COUNTER=memor_loc)
 
    call memory_counter_end(memor_loc,par % memor,MEMORY_COUNTER)

  end subroutine send_list_recv_points_1
  
  subroutine send_list_recv_points_2(par,xx,ll,nn_send,nn_recv,lista_send,xx_recv,ll_recv,METHOD,MASK,EXCLUDE_MYSELF,MEMORY_COUNTER)

    class(search_parall),                    intent(inout) :: par
    real(rp),                       pointer, intent(in)    :: xx(:,:)
    integer(ip),                    pointer, intent(in)    :: ll(:)
    integer(ip),                    pointer, intent(inout) :: nn_send(:)
    integer(ip),                    pointer, intent(inout) :: nn_recv(:)
    type(i1p),                      pointer, intent(inout) :: lista_send(:)  
    type(r2p),                      pointer, intent(inout) :: xx_recv(:)
    type(i1p),                      pointer, intent(inout) :: ll_recv(:)  
    integer(ip),      optional,              intent(in)    :: METHOD                !< Method for candidates
    logical(lg),      optional,     pointer, intent(in)    :: MASK(:)
    logical(lg),      optional,              intent(in)    :: EXCLUDE_MYSELF
    integer(8),       optional,              intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(8)                                             :: memor_loc(2)

    nullify(nn_send)
    nullify(nn_recv)
    nullify(lista_send)
    nullify(xx_recv)
    nullify(ll_recv)

    call memory_counter_ini(memor_loc,par % memor,MEMORY_COUNTER)

    call par % list_rank         (xx,lista_send,METHOD=METHOD,MASK=MASK,EXCLUDE_MYSELF=EXCLUDE_MYSELF,MEMORY_COUNTER=memor_loc)    
    call par % send_recv_points_2(xx,ll,nn_send,nn_recv,xx_recv,ll_recv,lista_send,MEMORY_COUNTER=memor_loc)
 
    call memory_counter_end(memor_loc,par % memor,MEMORY_COUNTER)

  end subroutine send_list_recv_points_2

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-17
  !> @brief   Communicator
  !> @details Construct parallel communicator
  !>
  !>          LISTA_SEND(IRANK) % L(:) ... List of points contained by rank IRANK
  !>          XX_RECV(IRANK)   % L(:) ... List of points received from rank IRANK
  !>
  !-----------------------------------------------------------------------

  subroutine dista_comm_1(par,nn,nn_send,nn_recv,lista_send,lista_recv,dista_recv,dista,MEMORY_COUNTER)
    class(search_parall),                    intent(inout) :: par
    integer(ip),                             intent(in)    :: nn
    integer(ip),                    pointer, intent(in)    :: nn_send(:)
    integer(ip),                    pointer, intent(in)    :: nn_recv(:)
    type(i1p),                      pointer, intent(in)    :: lista_send(:)
    type(i1p),                      pointer, intent(inout) :: lista_recv(:)
    type(r1p),                      pointer, intent(in)    :: dista_recv(:)
    real(rp),         optional,     pointer, intent(inout) :: dista(:) 
    integer(8),       optional,              intent(inout) :: MEMORY_COUNTER(2) !< Memory counters

    integer(ip)                                            :: irank,nrank
    integer(ip)                                            :: ii,kk
    type(r1p),                      pointer                :: dista_send(:)
    integer(ip),                    pointer                :: irank_min(:)
    real(rp),                       pointer                :: dista_min(:)
    integer(8)                                             :: memor_loc(2)

    nrank = int(par % comm % size4,ip)
    if( nrank <= 0 ) return
    
    call memory_counter_ini(memor_loc,par % memor,MEMORY_COUNTER)

    nullify(lista_recv)
    nullify(dista_send)
    nullify(irank_min)
    nullify(dista_min)
    ! 
    ! DISTA_RECV(IRANK) % A(NN_SEND(IRANK))
    !
    call par % send_recv_dista(nn_send,nn_recv,dista_recv,dista_send,MEMORY_COUNTER=memor_loc)
    !
    ! Decide who's the owner of point II according to distance
    !
    if( present(dista) ) then
       dista_min => dista
    else
       call memory_alloca(memor_loc,'DISTA_MIN',vacal,dista_min,nn,INIT_VALUE= huge(1.0_rp))
    end if
    call memory_alloca(memor_loc,'IRANK_MIN',vacal,irank_min,nn,INIT_VALUE=-1_ip)
      
    if( associated(lista_send) ) then
       do irank = 0,nrank-1                          
          do kk = 1,memory_size(lista_send(irank) % l)
             ii = lista_send(irank) % l(kk)
             if( abs(dista_send(irank) % a(kk)) <= abs(dista_min(ii)) .and. dista_send(irank) % a(kk) >= -0.1_rp*huge(1.0_rp) ) then
                irank_min(ii) = irank
                dista_min(ii) = dista_send(irank) % a(kk)
             end if
          end do
       end do
       !
       ! Prepare answer
       !
       do irank = 0,nrank-1
          do kk = 1,memory_size(lista_send(irank) % l)
             ii = lista_send(irank) % l(kk)
             if( irank_min(ii) /= irank ) then
                lista_send(irank) % l(kk) = 0
             else
                continue
             end if
          end do
       end do
    end if
    !
    ! Give answer to my neighbors and compute communicator
    !
    call par % communication(nn_send,nn_recv,lista_send,lista_recv,MEMORY_COUNTER=memor_loc)
    !
    ! Deallocate
    !
    if( .not. present(dista) ) call memory_deallo(memor_loc,'DISTA_MIN' ,vacal,dista_min)
    call memory_deallo(memor_loc,'IRANK_MIN' ,vacal,irank_min)
    call memory_deallo(memor_loc,'DISTA_SEND',vacal,dista_send)

    call memory_counter_end(memor_loc,par % memor,MEMORY_COUNTER)

  end subroutine dista_comm_1

  subroutine dista_comm_2(par,nn,nn_send,nn_recv,lista_send,lista_recv,ll_recv,MEMORY_COUNTER)
    class(search_parall),                    intent(inout) :: par
    integer(ip),                             intent(in)    :: nn
    integer(ip),                    pointer, intent(in)    :: nn_send(:)
    integer(ip),                    pointer, intent(in)    :: nn_recv(:)
    type(i1p),                      pointer, intent(in)    :: lista_send(:)
    type(i1p),                      pointer, intent(inout) :: lista_recv(:)
    type(i1p),                      pointer, intent(in)    :: ll_recv(:)
    integer(8),       optional,              intent(inout) :: MEMORY_COUNTER(2) !< Memory counters

    integer(ip)                                            :: irank,nrank
    integer(ip)                                            :: ii,kk
    type(i1p),                      pointer                :: ll_send(:)
    integer(8)                                             :: memor_loc(2)
    MY_MPI_COMM                                            :: comm4
    
    nrank = int(par % comm % size4,ip)
    comm4 = par % comm % PAR_COMM_WORLD
    if( nrank <= 0 ) return
    
    call memory_counter_ini(memor_loc,par % memor,MEMORY_COUNTER)

    nullify(lista_recv)
    nullify(ll_send)
    ! 
    ! DISTA_RECV(IRANK) % A(NN_SEND(IRANK))
    !
    if( .not. associated(ll_send) ) &
         call memory_alloca(memor_loc,'LL_SEND','def_search_parall',ll_send,nrank,LBOUN=0_ip)
    do irank = 0,nrank-1
       if( .not. associated(ll_send(irank) % l) ) &
            call memory_alloca(memor_loc,'LL_SEND % L','def_search_parall',ll_send(irank) % l,nn_send(irank))          
    end do
    do irank = 0,nrank-1       
       if( nn_recv(irank) > 0 .or. nn_send(irank) > 0 ) then
          call PAR_SEND_RECEIVE(ll_recv(irank) % l,ll_send(irank) % l,DOM_I=irank,PAR_COMM_IN=comm4)
       end if
    end do
    !
    ! Decide who's the owner of point II according to distance
    !
    if( associated(lista_send) ) then
       !
       ! Prepare answer
       !
       do irank = 0,nrank-1
          do kk = 1,memory_size(lista_send(irank) % l)
             ii = ll_send(irank) % l(kk)
             if( ii == 0 ) then
                lista_send(irank) % l(kk) = 0
             else
                continue
             end if
          end do
       end do
    end if
    !
    ! Give answer to my neighbors and compute communicator
    !
    call par % communication(nn_send,nn_recv,lista_send,lista_recv,MEMORY_COUNTER=memor_loc)
    !
    ! Deallocate
    !
    call memory_deallo(memor_loc,'LL_SEND'  ,vacal,ll_send)

    call memory_counter_end(memor_loc,par % memor,MEMORY_COUNTER)

  end subroutine dista_comm_2

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

end module def_search_parall
