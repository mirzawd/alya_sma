!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    mod_par_tools.f90
!> @author  houzeaux
!> @date    2018-10-31
!> @brief   Parallel tools
!> @details Tools for doing sequential things in parallel ;o)
!-----------------------------------------------------------------------

module mod_par_tools

  use def_kintyp_basic,      only : ip,rp,lg
  use def_kintyp_comm,       only : comm_data_par
  use def_kintyp_comm,       only : comm_data_par_basic
  use def_master,            only : npart
  use def_master,            only : INOTMASTER
  use def_master,            only : ISEQUEN
  use def_master,            only : ISLAVE
  use def_master,            only : kfl_paral
  use def_master,            only : lninv_loc
  use def_master,            only : THIS_NODE_IS_MINE
  use def_domain,            only : htable_lninv_loc
  use def_kintyp_dims,       only : npoin_own,npoin
  use def_kintyp_dims,       only : npoin,npoin_2
  use def_parall,            only : kfl_global_numbering_par
  use mod_maths,             only : maths_heap_sort
  use mod_communications,    only : PAR_ALLGATHER
  use mod_communications,    only : PAR_ALLGATHERV
  use mod_communications,    only : PAR_MAX
  use mod_communications,    only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications,    only : PAR_COMM_RANK_AND_SIZE
  use mod_communications,    only : PAR_SUM
  use mod_parall,            only : PAR_CODE_SIZE
  use mod_parall,            only : par_memor
  use mod_memory,            only : memory_alloca
  use mod_memory,            only : memory_deallo
  use mod_memory,            only : memory_size
  use mod_maths,             only : maths_merge_ordered_lists
  use mod_messages,          only : messages_live
  use mod_htable,            only : htalid
  use mod_optional_argument, only : optional_argument
  
  use mod_std  

  implicit none

  public :: par_tools_merge_lists    ! Merges lists generated in parallel
  public :: par_tools_gathered_graph ! gather sequential graphs at root
  public :: par_tools_global_numbering
  public :: par_tools_ownership_interfaces
  public :: par_tools_comm_to_array
  public :: par_subdomain_bb
  public :: PAR_GLOBAL_TO_LOCAL_NODE
  public :: PAR_THIS_NODE_IS_MINE
  
contains

  !-----------------------------------------------------------------------
  !>  
  !> @author  houzeaux
  !> @date    2018-10-30
  !> @brief   Merge parallel lists
  !> @details Marge lists in parallel, each one being generated by the slaves
  !> 
  !-----------------------------------------------------------------------

  subroutine par_tools_merge_lists(array,memor,VARIABLE_NAME)

    integer(ip), pointer, intent(inout)          :: array(:)
    integer(8),           intent(inout)          :: memor(2)
    character(len=*),     intent(in),   optional :: VARIABLE_NAME
    integer(ip), pointer                         :: array_size_gat(:)
    integer(ip), pointer                         :: array_gat(:)
    integer(ip)                                  :: array_size
    integer(ip)                                  :: ii,jj,kk,total_size
    integer(ip)                                  :: minvalue

    minvalue   = 0
    array_size = memory_size(array)
    total_size = array_size
    
    call PAR_MAX(total_size)
    
    if( total_size > 0 ) then

       nullify(array_gat)

       allocate(array_size_gat(0:PAR_CODE_SIZE-1))   
       call PAR_ALLGATHER(array_size,array_size_gat)

       total_size = sum(array_size_gat)
       allocate(array_gat(total_size))
       call PAR_ALLGATHERV(array,array_gat,array_size_gat)

       call maths_heap_sort(2_ip,total_size,array_gat)
       
       jj = minvalue
       kk = 0
       do ii = 1,total_size
          if( array_gat(ii) <= jj ) then
             array_gat(ii) = minvalue
          else
             kk            = kk + 1
             array_gat(kk) = array_gat(ii)
             jj            = array_gat(ii)
          end if
       end do
       if( kk < total_size ) array_gat(kk+1:) = minvalue

       total_size = count(array_gat/=minvalue,KIND=ip)
       if( present(VARIABLE_NAME) ) then
          call memory_deallo(memor,trim(VARIABLE_NAME),'par_tools_merge_list',array)
          call memory_alloca(memor,trim(VARIABLE_NAME),'par_tools_merge_list',array,total_size)
       else
          call memory_deallo(memor,'ARRAY','par_tools_merge_list',array)
          call memory_alloca(memor,'ARRAY','par_tools_merge_list',array,total_size)
       end if
       array(1:total_size) = array_gat(1:total_size)

       deallocate(array_size_gat)
       deallocate(array_gat)

    end if

  end subroutine par_tools_merge_lists

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-10-31
  !> @brief   Merge graphs
  !> @details Marge graphs generated in parallel into a sequential one
  !> 
  !-----------------------------------------------------------------------

  subroutine par_tools_gathered_graph(ja_type,lninv_gat,MEMORY_COUNTER)

    use def_kintyp
    use def_domain
    use def_master
    use mod_communications
    use mod_parall
    use mod_memory
    use mod_domain, only : domain_total_node_number
    implicit none
    type(i1p),   pointer,  intent(inout) :: ja_type(:)
    integer(ip), pointer,  intent(inout) :: lninv_gat(:)
    integer(8),  optional, intent(inout) :: MEMORY_COUNTER(2)
    integer(ip), pointer                 :: inul1(:)
    integer(ip), pointer                 :: inul2(:,:)

    integer(ip), pointer                 :: nz_gat(:)
    integer(4),  pointer                 :: nn_gat(:)
    integer(ip), pointer                 :: ia_gat(:)
    integer(ip), pointer                 :: ja_gat(:)
    !integer(ip), pointer                 :: lninv_gat(:)

    integer(ip), pointer                 :: lninv_send(:)
    integer(ip), pointer                 :: ja_send(:)

    integer(ip)                          :: inul0(1)
    integer(ip)                          :: nn,nz,ii,izdom
    integer(ip)                          :: index_ipoin
    integer(ip)                          :: isize,ipoin,apoin,npoin_mesh
    integer(ip)                          :: ipart,nz_typ,imeth

    integer(8)                           :: memor(2)

    if( present(MEMORY_COUNTER) ) then
       memor = MEMORY_COUNTER
    else
       memor = par_memor
    end if

    if( ISEQUEN ) then

       call runend('PAR_TOOLS_GATHERED_GRAPH: IMPOSSIBLE IN SEQUENTIAL')

    else

       nullify(inul1)
       nullify(inul2)
       nullify(nz_gat)
       nullify(nn_gat)
       nullify(ia_gat) 
       nullify(ja_gat)
       nullify(lninv_gat)
       nullify(lninv_send)
       nullify(ja_send)
       nullify(ja_type)

       imeth = 0

       if( IMASTER ) then

          call memory_alloca(memor,'NZ_GAT','par_tools_gathered_graph',nz_gat,PAR_CODE_SIZE,lboun=0_ip)
          call memory_alloca(memor,'NN_GAT','par_tools_gathered_graph',nn_gat,int(PAR_CODE_SIZE,4),lboun=0_4)
          nz              = 0_ip
          nn_gat(0)       = 0
          nn_gat(1:npart) = npoin_par(1:npart)

       else
          !
          ! Slaves mark boudnary node with negative sign
          !
          nz         = nzdom
          nn         = npoin
          npoin_mesh = npoin_own
          call memory_alloca(memor,'LNINV_SEND','par_tools_gathered_graph',lninv_send,npoin,'DO_NOT_INITIALIZE')
          call memory_alloca(memor,'JA_SEND'   ,'par_tools_gathered_graph',ja_send,   nzdom,'DO_NOT_INITIALIZE')
          do ipoin = 1,npoi1
             lninv_send(ipoin) =  lninv_loc(ipoin)
          end do
          do ipoin = npoi1+1,npoin             
             lninv_send(ipoin) = -lninv_loc(ipoin)
          end do
          do izdom = 1,nzdom
             ja_send(izdom)    =  lninv_loc(c_dom(izdom))
          end do
          do ipoin = 1,npoin
             isize = r_dom(ipoin+1)-r_dom(ipoin)
             call maths_heap_sort(2_ip,isize,ja_send(r_dom(ipoin):r_dom(ipoin+1)-1))
          end do
       end if
       !
       ! NPOIN_MESH: number of nodes (excluding duplicated ndoes)
       ! Gather graph sizes NZ_GAT and global numbering LNINV_GAT
       !
       npoin_mesh = domain_total_node_number()
       call PAR_GATHER(nz,nz_gat)

       if( IMASTER ) then
          call memory_alloca(memor,'JA_TYPE'  ,'par_tools_gathered_graph',ja_type  ,npoin_mesh)
          call memory_alloca(memor,'LNINV_GAT','par_tools_gathered_graph',lninv_gat,npoin_total,'DO_NOT_INITIALIZE')
       end if

       call PAR_GATHERV(lninv_send,lninv_gat,nn_gat)
       !
       ! Merge graphs into type JA_TYPE
       !
       if( IMASTER ) then
          index_ipoin = 0
          do ipart = 1,PAR_CODE_SIZE-1

             call memory_alloca(memor,'IA_GAT'   ,'par_tools_gathered_graph',ia_gat   ,npoin_par(ipart)+1_ip,'DO_NOT_INITIALIZE')
             call memory_alloca(memor,'JA_GAT'   ,'par_tools_gathered_graph',ja_gat   ,nz_gat(ipart)        ,'DO_NOT_INITIALIZE')

             call PAR_SEND_RECEIVE(inul1,ia_gat   ,'IN MY CODE',ipart,'SYNCHRONOUS')
             call PAR_SEND_RECEIVE(inul1,ja_gat   ,'IN MY CODE',ipart,'SYNCHRONOUS')

                do ii = 1,npoin_par(ipart)
                   nz    = ia_gat(ii+1) - ia_gat(ii)
                   ipoin = lninv_gat(ii+index_ipoin)
                   apoin = abs(ipoin)
                   if( ipoin > 0 ) then
                      !
                      ! Interior node: graph is necessarily complete
                      !
                      call memory_alloca(memor,'JA_TYPE % L','par_tools_gathered_graph',ja_type(ipoin) % l,nz)
                      ja_type(ipoin) % l(1:nz) = ja_gat(ia_gat(ii):ia_gat(ii+1)-1)
                   else
                      !
                      ! Boundary node: merge list
                      !
                      if( associated(ja_type(apoin) % l) ) then
                         nz_typ = memory_size(ja_type(apoin) % l)
                         call maths_merge_ordered_lists(nz,ja_gat(ia_gat(ii):),nz_typ,ja_type(apoin) % l,MEMORY_COUNTER=memor,LIST_NAME='JA_TYPE % L')  
                      else
                         nz_typ = int(1.5_rp*real(nz,rp))
                         call memory_alloca(memor,'JA_TYPE % L','par_tools_gathered_graph',ja_type(apoin) % l,nz_typ)
                         ja_type(apoin) % l(1:nz) = ja_gat(ia_gat(ii):ia_gat(ii+1)-1)
                      end if
                   end if
                   lninv_gat(ii+index_ipoin) = apoin
                end do
                
             index_ipoin = index_ipoin + npoin_par(ipart)
             call memory_deallo(memor,'JA_GAT'   ,'par_tools_gathered_graph',ja_gat)
             call memory_deallo(memor,'IA_GAT'   ,'par_tools_gathered_graph',ia_gat)
          end do

          do ipoin = 1,npoin_mesh
             nz = 0
             if( associated(ja_type(ipoin) % l) ) then
                ii_loop: do ii = 1,memory_size(ja_type(ipoin) % l)
                   if( ja_type(ipoin) % l(ii) == 0 ) then
                      call memory_resize(memor,'JA_TYPE % L','par_tools_gathered_graph',ja_type(ipoin) % l,ii-1_ip)
                      exit ii_loop
                   end if
                end do ii_loop
             end if
          end do

          !print*,'cacac=',npoin_mesh
          !do ipoin = 1,npoin_mesh
          !   if( associated(ja_type(ipoin)%l) ) then
          !      print*,'ppo=',ipoin,': ',ja_type(ipoin)%l
          !   else
          !      print*,'caca'
          !   end if
          !end do
          !call runend('O.K.!')
          !
          ! Resize graph
          !
          !if( IMASTER ) then
          !   do ipoin = 1,npoin_mesh
          !      nz = size(ja_type(ipoin) % l)
          !      if( ja_type(ipoin) % l(nz) == 0 ) then
          !         do while( ja_type(ipoin) % l(nz) == 0 )
          !            nz = nz - 1
          !         end do
          !         call memory_resize(memor,'JA_TYPE % L','par_tools_gathered_graph',ja_type(ipoin) % l,nz)
          !      end if
          !   end do
          !end if
          
          call memory_deallo(memor,'NZ_GAT','par_tools_gathered_graph',nz_gat)
          call memory_deallo(memor,'NN_GAT','par_tools_gathered_graph',nn_gat)

       else

          call PAR_SEND_RECEIVE(npoin+1_ip,0_ip,r_dom,inul0,'IN MY CODE',0_ip,'SYNCHRONOUS')
          if( nz > 0 ) call PAR_SEND_RECEIVE(nz,0_ip,ja_send,inul0,'IN MY CODE',0_ip,'SYNCHRONOUS')

          call memory_deallo(memor,'LNINV_SEND','par_tools_gathered_graph',lninv_send)          
          call memory_deallo(memor,'JA_SEND'   ,'par_tools_gathered_graph',ja_send)          

       end if

    end if

    if( present(MEMORY_COUNTER) ) then
       MEMORY_COUNTER = memor 
    else
       par_memor = memor
    end if

  end subroutine par_tools_gathered_graph

  !-----------------------------------------------------------------------
  !
  !> @author  houzeaux
  !> @date    2019-01-02
  !> @brief   Global numbering
  !> @details Global numbering
  !
  !-----------------------------------------------------------------------

  subroutine par_tools_global_numbering(NUMBERING,lninv_opt)

    integer(ip), intent(in),             optional :: NUMBERING
    integer(ip), intent(inout), pointer, optional :: lninv_opt(:)
    integer(ip)                                   :: ipoin,ipart,npoin_offset
    integer(ip), pointer                          :: npoin_gat(:)
    integer(ip)                                   :: numbering_loc
    integer(ip), pointer                          :: lninv_use(:)
    
    if( present(NUMBERING) ) then
       numbering_loc = NUMBERING
    else
       numbering_loc = kfl_global_numbering_par
    end if
    if( present(lninv_opt) ) then
       lninv_use => lninv_opt
       if( .not. associated(lninv_use) ) then
          call memory_alloca(par_memor,'LNINV_OPT','par_global_numbering',lninv_use,npoin)   
       end if
    else
       lninv_use => lninv_loc
    end if
  
    nullify(npoin_gat)
    
    if( numbering_loc == 0 ) then
       !
       ! Do nothing, take the global numbering given by the mesh generator
       !
    else if( numbering_loc == 1 ) then
       !
       ! Lexical numbering
       !
       call messages_live('COMPUTE GLOBAL NUMBERING')
       call memory_alloca(par_memor,'NPOIN_GAT','par_global_numbering',npoin_gat,npart+1_ip,LBOUN=0_ip)
       call PAR_ALLGATHER(npoin_own,npoin_gat,1_4,'IN MY CODE')
       npoin_offset = 0
       do ipart = 1,kfl_paral-1
          npoin_offset = npoin_offset + npoin_gat(ipart)
       end do
       if(INOTMASTER ) lninv_use = 0_ip
       do ipoin = 1,npoin_own
          lninv_use(ipoin) = ipoin + npoin_offset
       end do
       call memory_deallo(par_memor,'NPOIN_GAT','par_global_numbering',npoin_gat)
       call PAR_INTERFACE_NODE_EXCHANGE(lninv_use,'SUM','IN MY CODE')
       
    end if

  end subroutine par_tools_global_numbering

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-06
  !> @brief   Dertrmine Ownerhsip on interfaces
  !> @details Determine own nodes on interfaces
  !> 
  !-----------------------------------------------------------------------

  subroutine par_tools_ownership_interfaces(npoin_new,comm,permr)

    integer(ip),                  intent(in)    :: npoin_new
    type(comm_data_par),          intent(inout) :: comm
    integer(ip),         pointer, intent(inout) :: permr(:)
    integer(ip)                                 :: ineig,ipoin,ii
    integer(ip)                                 :: dom_i,isize,lsize
    integer(ip)                                 :: kpoin
    integer(ip),         pointer                :: node_mult(:)
    integer(ip),         pointer                :: node_type(:)
    integer(ip),         pointer                :: size_inte(:)
    integer(ip)                                 :: PAR_CURRENT_RANK

    PAR_CURRENT_RANK = int(comm % RANK4,ip)

    nullify(node_mult)
    nullify(node_type)
    nullify(size_inte)
    call memory_alloca(par_memor,'NODE_MULT','par_tools_ownership_interfaces',node_mult,npoin_new)
    call memory_alloca(par_memor,'NODE_TYPE','par_tools_ownership_interfaces',node_type,npoin_new)
    call memory_alloca(par_memor,'SIZE_INTE','par_tools_ownership_interfaces',size_inte,comm % nneig)
    if( .not. associated(permr) ) then
       call memory_alloca(par_memor,'PERMR','par_tools_ownership_interfaces',permr,npoin_new)
    end if
    !
    ! SIZE_INTE: order nodes on both sides of the interface
    !
    do ineig = 1,comm % nneig
       lsize            = comm % bound_size(ineig+1)-comm % bound_size(ineig)
       size_inte(ineig) = lsize
    end do
    !
    ! NODE_MULT: compute multiplicity
    !
    do ii = 1,comm % bound_dim
       ipoin            = comm % bound_perm(ii)
       node_type(ipoin) = 2
       node_mult(ipoin) = node_mult(ipoin) + 1
    end do
    !
    ! Split boundary between neighbors
    !
    do ineig = 1,comm % nneig
       dom_i = comm % neights(ineig)
       isize = comm % bound_size(ineig+1)-comm % bound_size(ineig)
       if( isize > 0 ) then
          lsize = size_inte(ineig)/2
          if( PAR_CURRENT_RANK < dom_i ) then
             !
             ! Take the first chunk
             !
             do ii = comm % bound_size(ineig),comm % bound_size(ineig)+lsize
                ipoin = comm % bound_perm(ii)
                if( node_mult(ipoin) <= 1 ) node_type(ipoin) = 1
             end do
          else
             !
             ! Take the second chunk
             !
             do ii = comm % bound_size(ineig)+lsize+1,comm % bound_size(ineig+1)-1
                ipoin = comm % bound_perm(ii)
                if( node_mult(ipoin) <= 1 ) node_type(ipoin) = 1
             end do
          end if
       end if
    end do
    !
    ! Mark nodes with multiplicity > 2 (see node * in next example). It will be
    ! given to the lowest rank MPI process
    !
    !       o
    !       |
    !   1   o
    !       |    3
    !  -----*
    !       |
    !   2   o
    !
    do ipoin = 1,npoin_new
       if( node_mult(ipoin) > 1 ) then
          node_mult(ipoin) = PAR_CURRENT_RANK
       else
          node_mult(ipoin) = 0
       end if
    end do
    do ineig = 1,comm % nneig
       dom_i = comm % neights(ineig)
       do ii = comm % bound_size(ineig),comm % bound_size(ineig+1)-1
          ipoin = comm % bound_perm(ii)
          if( node_mult(ipoin) > 0 ) node_mult(ipoin) = min(node_mult(ipoin),dom_i)
       end do
    end do
    do ipoin = 1,npoin_new
       if( node_mult(ipoin) > 0 ) then
          if( PAR_CURRENT_RANK == node_mult(ipoin) ) then
             node_type(ipoin) = 1
          else
             node_type(ipoin) = 2
          end if
       end if
    end do
    !
    ! Compute permutation to irder interior and boundary nodes
    !
    kpoin = 0
    do ipoin = 1,npoin_new
       if( node_type(ipoin) == 0 ) then
          kpoin = kpoin + 1
          permr(ipoin) = kpoin
       end if
    end do
    comm % npoi1 = kpoin
    comm % npoi2 = kpoin+1
    do ipoin = 1,npoin_new
       if( node_type(ipoin) == 1 ) then
          kpoin = kpoin + 1
          permr(ipoin) = kpoin
       end if
    end do
    comm % npoi3 = kpoin
    do ipoin = 1,npoin_new
       if( node_type(ipoin) == 2 ) then
          kpoin = kpoin + 1
          permr(ipoin) = kpoin
       end if
    end do
    !
    ! Deallocate
    !
    call memory_deallo(par_memor,'NODE_MULT','par_tools_ownership_interfaces',node_mult)
    call memory_deallo(par_memor,'NODE_TYPE','par_tools_ownership_interfaces',node_type)
    call memory_deallo(par_memor,'SIZE_INTE','par_tools_ownership_interfaces',size_inte)
    
  end subroutine par_tools_ownership_interfaces

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    06/03/2014
  !> @brief   Returns the local node number of a global number
  !> @details Returns the local node number of a global number
  !>          and zero if I do not own the node
  !
  !----------------------------------------------------------------------

  function PAR_GLOBAL_TO_LOCAL_NODE(ipoin_global,wherein,mask)

    integer(ip),  intent(in)                    :: ipoin_global
    character(*), intent(in), optional          :: wherein
    logical(lg),  intent(in), optional, pointer :: mask(:)
    integer(ip)                                 :: ipoin
    integer(ip)                                 :: PAR_GLOBAL_TO_LOCAL_NODE
    integer(ip)                                 :: npoin_end

    npoin_end = npoin
    if( present(wherein) ) then
       if( trim(wherein) == 'INCLUDING HALOS' ) then
          npoin_end = npoin_2
       end if
    end if

    PAR_GLOBAL_TO_LOCAL_NODE = 0

    if( ISEQUEN ) then
       PAR_GLOBAL_TO_LOCAL_NODE = ipoin_global
    else if( ISLAVE ) then
       if( htable_lninv_loc % sizet /= 0 ) then
          ipoin = htalid(htable_lninv_loc,ipoin_global)
          if( ipoin < 1 .or. ipoin > npoin_end ) ipoin = 0
          PAR_GLOBAL_TO_LOCAL_NODE = ipoin
          if( present(mask) ) then
             if( .not. mask(ipoin) ) PAR_GLOBAL_TO_LOCAL_NODE = 0
          end if
          return
       else
          if( present(mask) ) then
             do ipoin = 1,npoin_end
                if( mask(ipoin) ) then
                   if( lninv_loc(ipoin) == ipoin_global ) then
                      PAR_GLOBAL_TO_LOCAL_NODE = ipoin
                      return
                   end if
                end if
             end do
          else
             do ipoin = 1,npoin_end
                if( lninv_loc(ipoin) == ipoin_global ) then
                   PAR_GLOBAL_TO_LOCAL_NODE = ipoin
                   return
                end if
             end do
          end if
       end if
    end if

  end function PAR_GLOBAL_TO_LOCAL_NODE

 !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    06/03/2014
  !> @brief   Returns if a node is interior or inside my own boundary
  !> @details Returns if a node is interior or inside my own boundary
  !>
  !----------------------------------------------------------------------

  function PAR_THIS_NODE_IS_MINE(ipoin,where)
    integer(ip),  intent(in)           :: ipoin
    character(*), intent(in), optional :: where
    logical(lg)                        :: PAR_THIS_NODE_IS_MINE

    PAR_THIS_NODE_IS_MINE = THIS_NODE_IS_MINE(ipoin)

  end function PAR_THIS_NODE_IS_MINE
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-07
  !> @brief   Communicator to array
  !> @details Convert a communicator into an array
  !> 
  !-----------------------------------------------------------------------
  
  subroutine par_tools_comm_to_array(COMM,nn,xx,NAME,GATHER)
    type(comm_data_par),                     intent(in)    :: comm
    integer(ip),                             intent(in)    :: nn
    integer(ip),                   pointer,  intent(inout) :: xx(:,:)
    character(len=*),                        intent(in)    :: NAME
    logical(lg),       optional,             intent(in)    :: GATHER
    character(200)                                         :: varname
    integer(ip)                                            :: kdime,ii,kk
    integer(ip)                                            :: ineig,dom_i
    integer(ip),                   pointer                 :: num_neigh(:)
    logical(lg)                                            :: my_gather

    nullify(num_neigh)

    varname = optional_argument('XX',NAME)
    my_gather = optional_argument(.true.,GATHER)

    call memory_alloca(par_memor,'NUM_NEIGH','par_tools',num_neigh,nn)   

    if( my_gather ) then
       do kk = 1,comm % bound_dim
          ii = comm % bound_perm(kk)
          num_neigh(ii) = num_neigh(ii) + 1
       end do
    else
       do kk = 1,comm % bound_dim
          ii = comm % bound_invp(kk)
          num_neigh(ii) = num_neigh(ii) + 1
       end do
    end if

    kdime = 0
    if( nn > 0 ) then
       kdime = maxval(num_neigh)
       do ii = 1,nn
          num_neigh(ii) = 0
       end do
    end if

    call PAR_MAX(kdime)

    if( associated(xx) ) then
       call runend('PAR_TOOLS: ARRAY ALREADY ASSOCIATED')
    else
       call memory_alloca(par_memor,trim(varname),'par_tools',xx,kdime,nn) 
    end if

    if( my_gather ) then
       do ineig = 1,comm % nneig
          dom_i = comm % neights(ineig)
          do kk = comm % bound_size(ineig),comm % bound_size(ineig+1)-1
             ii = comm % bound_perm(kk)
             num_neigh(ii)        = num_neigh(ii) + 1
             xx(num_neigh(ii),ii) = dom_i
          end do
       end do
    else
       do ineig = 1,comm % nneig
          dom_i = comm % neights(ineig)
          do kk = comm % bound_size(ineig),comm % bound_size(ineig+1)-1
             ii = comm % bound_invp(kk)
             num_neigh(ii)        = num_neigh(ii) + 1
             xx(num_neigh(ii),ii) = dom_i
          end do
       end do
    end if

  end subroutine par_tools_comm_to_array

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    03/11/2021
  !> @brief   Subdomain bounding boxes
  !> @details Subdomain bounding boxes from own bounding boxes
  !
  !----------------------------------------------------------------------

  subroutine par_subdomain_bb(COMM,bobox,subox,MEMORY_COUNTER,DEALLOCATE_ONLY)

    class(comm_data_par_basic),                   intent(inout) :: COMM               !< Input communicator
    real(rp),                            pointer, intent(in)    :: bobox(:,:,:)
    real(rp),                            pointer, intent(inout) :: subox(:,:,:)
    integer(8),                optional,          intent(inout) :: MEMORY_COUNTER(2)  !< Memory counter
    logical(lg),               optional,          intent(in)    :: DEALLOCATE_ONLY    !< Order
    integer(8)                                                  :: memor(2)
    integer(ip)                                                 :: nrank,irank
    integer(ip)                                                 :: pdime,idime
    logical(lg)                                                 :: if_deallocate

    memor         = optional_argument(par_memor,MEMORY_COUNTER)
    if_deallocate = optional_argument(.false.,DEALLOCATE_ONLY)

    if( if_deallocate ) then
       
       call memory_deallo(memor,'SUBOX','par_tools',subox)
       
    else

       call PAR_COMM_RANK_AND_SIZE(COMM,irank,nrank)
   
       pdime = memory_size(bobox,2_ip)
       
       call PAR_MAX(pdime,COMM,INCLUDE_ROOT=.true.) 
       
       call memory_alloca(memor,'SUBOX','par_tools',subox,2_ip,pdime,nrank,LBOUN1=1_ip,LBOUN2=1_ip,LBOUN3=0_ip) 

       if( memory_size(bobox) > 0 ) then
          do idime = 1,pdime
             subox(1,idime,irank) = minval(bobox(1,idime,:))
             subox(2,idime,irank) = maxval(bobox(2,idime,:))
          end do
       else
          subox(1,:,irank) =  huge(1.0_rp)*0.1_rp
          subox(2,:,irank) = -huge(1.0_rp)*0.1_rp
       end if

       call PAR_SUM(subox,COMM,INCLUDE_ROOT=.true.) 
              
    end if

    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor

  end subroutine par_subdomain_bb
  
end module mod_par_tools
  !> @}
