!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_edge_communication_array.f90
!> @date    04/02/2016
!> @author  Guillaume Houzeaux
!> @brief   Edge communication array 
!> @details Construct edge communication arrays from
!>          node communication arrays:
!>
!>          COMM % BEDGE_DIM ....... Totla size
!>          COMM % BEDGE_SIZE(:) ... Permutation list 
!>          COMM % BEDGE_PERM(:) ... Permutaiton array
!>          COMM % BEDG1 ........... Interior edges
!>          COMM % BEDG2 ........... Start of own boundary edges (=COMM % BEDG1+1)
!>          COMM % BEDG3 ........... End of own boundary edges
!>
!>          Edge renumbering:
!>
!>                 Interior        Own boundary    Other boundary
!>          <--------------------><------------><------------------>
!>
!>          +--+---------------+--+--+------+--+----------------+--+
!>          |  |               |  |  |      |  |                |  |
!>          +--+---------------+--+--+------+--+----------------+--+
!>            1              BEDG1 BEDG2   BEDG3           MESHE % NEDGE
!>
!>          Then renumber edge arrays:
!>
!>          LEDGS(1:PEDGE,1:NELEM) ...... List of global edges for elements
!>          LEDGB(1:PEDGE,1:NBOUN) ...... List of global edges for boundaries
!>          R_EDG(:), C_EDG(:) .......... Edge graph
!>          EDGE_TO_NODE(1:2,1:NEDGE) ... Edge nodes
!>
!> @} 
!-----------------------------------------------------------------------

subroutine par_edge_communication_array(meshe,COMM)
  use def_kintyp_basic,   only : ip,lg,i1p
  use def_kintyp_comm,    only : comm_data_par
  use def_master,         only : intost
  use def_master,         only : nedg1,nedg2,nedg3
  use def_master,         only : kfl_paral,INOTMASTER
  use def_domain,         only : mesh_type,memor_dom
  use mod_parall,         only : par_memor 
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_memory,         only : memory_size
  use mod_graphs,         only : graphs_number_to_linked_list
  use mod_communications, only : PAR_INTERFACE_EDGE_EXCHANGE
  use mod_communications, only : PAR_SEND_RECEIVE
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_SUM
  use mod_maths,          only : maths_array_permutation
  use mod_graphs,         only : graphs_rengra
  implicit none
  type(mesh_type),     intent(inout) :: meshe                !< Mesh type
  type(comm_data_par), intent(inout) :: COMM                 !< Communication array
  integer(ip),         parameter     :: INTERIOR     = 0
  integer(ip),         parameter     :: OWN_BOUNDARY = 1
  integer(ip),         parameter     :: OTH_BOUNDARY = 2
  integer(ip)                        :: jj,ineig,isize
  integer(ip)                        :: ipoin,jpoin,ierro
  integer(ip)                        :: ipoin1,ipoin2
  integer(ip)                        :: iedgg,kk,iedge
  integer(ip)                        :: counti,dom_i
  integer(ip)                        :: max_edges,jj0
  integer(ip)                        :: bsize,nedge_domi
  integer(ip)                        :: npoin_int,npoin_own
  integer(ip)                        :: npoin_oth,pedge
  integer(ip)                        :: ielem,iboun
  logical(lg)                        :: ifoun
  integer(ip),         pointer       :: permu(:)             ! Temporary permutation
  integer(ip),         pointer       :: invpr(:)             ! Temporary inv permutation
  integer(ip),         pointer       :: edge_type(:)         ! Edge type
  integer(ip),         pointer       :: edge_owner(:)        ! Edge type sum
  integer(ip),         pointer       :: number_edges(:)      ! Number of edges per node
  integer(ip),         pointer       :: global_edges(:)      ! Global edge numbers
  integer(ip),         pointer       :: list_send_edges(:,:) ! List of edges per node
  integer(ip),         pointer       :: list_recv_edges(:,:) ! List of edges per node
  integer(ip),         pointer       :: iexchange(:)         ! Dummy array for exchange
  logical(lg),         pointer       :: interface_edge(:)    ! Interface edge
  type(i1p),           pointer       :: list_neigh_edges(:)
  !
  ! Allocate memory
  ! 
  nullify( permu            )
  nullify( invpr            )
  nullify( edge_type        ) 
  nullify( edge_owner       ) 
  nullify( number_edges     ) 
  nullify( global_edges     ) 
  nullify( list_send_edges  ) 
  nullify( list_recv_edges  ) 
  nullify( iexchange        ) 
  nullify( interface_edge   ) 
  nullify( list_neigh_edges )

if( INOTMASTER ) then
     !
     ! Allocate memory
     !
     call memory_alloca(par_memor,'BEDGE_SIZE'      ,'par_edge_communication_array',COMM % bedge_size,COMM % nneig+1_ip)
     call memory_alloca(par_memor,'EDGE_TYPE'       ,'par_edge_communication_array',edge_type        ,meshe % nedge)
     call memory_alloca(par_memor,'EDGE_OWNER'      ,'par_edge_communication_array',edge_owner       ,meshe % nedge)
     call memory_alloca(par_memor,'INTERFACE_EDGES' ,'par_edge_communication_array',interface_edge   ,meshe % nedge)
     call memory_alloca(par_memor,'NUMBER_EDGES'    ,'par_edge_communication_array',number_edges     ,meshe % npoin)
     call memory_alloca(par_memor,'LIST_NEIGH_EDGES','par_edge_communication_array',list_neigh_edges ,COMM % nneig)
     !
     ! Detect possible boundary edges
     !
     do iedgg = 1,meshe % nedge
        ipoin = meshe % edge_to_node(1,iedgg)
        jpoin = meshe % edge_to_node(2,iedgg)    
        if( ipoin > COMM % npoi1 .and. jpoin > COMM % npoi1 ) then
           edge_type(iedgg)    = 1
           number_edges(ipoin) = number_edges(ipoin) + 1
           number_edges(jpoin) = number_edges(jpoin) + 1
        end if
     end do
     max_edges = 0
     do ipoin = COMM % npoi1+1,meshe % npoin
        max_edges = max(max_edges,number_edges(ipoin))
     end do
  end if
  call PAR_MAX(max_edges,'IN MY CODE')
  !
  ! Two boundary nodes are not necessarily an edge
  ! For partition 1: 1-2 is an edge
  ! For partition 2: 1-2 is NOT an edge!
  ! For partition 3: 1-2 is an edge
  !
  !       o-----o
  !       |  1  |
  !       |     |
  !       1-----2
  ! o-----1-----2-----o
  ! |  2  |  3  |  2  |
  ! |     |     |     |
  ! o-----o-----o-----o
  ! |  2  |  2  |  2  |
  ! |     |     |     |
  ! o-----o-----o-----o
  !
  ! BSIZE is an estimate of the max number of edges I share with a given neighbor
  !
  if( INOTMASTER ) then

     COMM % bedge_dim = 0
     do ineig = 1,COMM % nneig 
        dom_i      = COMM % neights(ineig)       
        nedge_domi = 0
        bsize      = ( COMM % bound_size(ineig+1)-COMM % bound_size(ineig) ) * max_edges / 2 + 1
        call memory_alloca(par_memor,'LIST_SEND_EDGES','par_edge_communication_array',list_send_edges,2_ip,bsize)
        call memory_alloca(par_memor,'LIST_RECV_EDGES','par_edge_communication_array',list_recv_edges,2_ip,bsize)
        call memory_alloca(par_memor,'GLOBAL_EDGES'   ,'par_edge_communication_array',global_edges,bsize)
        call memory_alloca(par_memor,'PERMU'          ,'par_edge_communication_array',permu,bsize) 
        jj0 = 0
        do jj = COMM % bound_size(ineig),COMM % bound_size(ineig+1)-1
           ipoin = COMM % bound_perm(jj)
           jj0   = jj0 + 1
           do kk = jj+1,COMM % bound_size(ineig+1)-1
              !
              ! Check if ipoin1-ipoin2 is an edge
              !
              jpoin  = COMM % bound_perm(kk)
              ipoin1 = min(ipoin,jpoin)
              ipoin2 = max(ipoin,jpoin)
              ifoun  = .false.
              loop_edges: do iedgg = 1,meshe % nedge
                 if( edge_type(iedgg) == 1 ) then
                    if( meshe % edge_to_node(1,iedgg) == ipoin1 .and. meshe % edge_to_node(2,iedgg) == ipoin2 ) then
                       ifoun = .true.
                       exit loop_edges
                    end if
                 end if
              end do loop_edges
              !
              ! It is an edge, add it to the communication list
              !
              if( ifoun ) then
                 nedge_domi                    = nedge_domi + 1
                 ipoin1                        = meshe % lninv_loc(ipoin)
                 ipoin2                        = meshe % lninv_loc(jpoin)
                 list_send_edges(1,nedge_domi) = min(ipoin1,ipoin2)
                 list_send_edges(2,nedge_domi) = max(ipoin1,ipoin2)
                 global_edges(nedge_domi)      = iedgg
              end if

           end do
        end do
        !
        ! Send/receive edge list
        !
        call PAR_SEND_RECEIVE(list_send_edges,list_recv_edges,'IN MY CODE',dom_i)
        !
        ! Check edges in common
        !     
        do isize = 1,bsize
           ipoin1 = list_recv_edges(1,isize)
           ipoin2 = list_recv_edges(2,isize)
           loop_jsize: do iedge = 1,nedge_domi
              if( list_send_edges(1,iedge) == ipoin1 .and. list_send_edges(2,iedge) == ipoin2 ) then
                 iedgg                           = global_edges(iedge)
                 COMM % bedge_dim                = COMM % bedge_dim  + 1
                 COMM % bedge_size(ineig)        = COMM % bedge_size(ineig) + 1
                 permu(COMM % bedge_size(ineig)) = iedgg
                 interface_edge(iedgg)           = .true.
                 exit loop_jsize 
              end if
           end do loop_jsize
        end do
        call memory_alloca(par_memor,'LIST_NEIGH_EDGES','par_edge_communication_array',list_neigh_edges(ineig) % l,COMM % bedge_size(ineig))
        do iedge = 1,COMM % bedge_size(ineig)
           list_neigh_edges(ineig) % l(iedge) = permu(iedge) 
        end do
        !
        ! Deallocate
        !
        call memory_deallo(par_memor,'LIST_SEND_EDGES','par_edge_communication_array',list_send_edges)
        call memory_deallo(par_memor,'LIST_RECV_EDGES','par_edge_communication_array',list_recv_edges)
        call memory_deallo(par_memor,'GLOBAL_EDGES'   ,'par_edge_communication_array',global_edges)
        call memory_deallo(par_memor,'PERMU'          ,'par_edge_communication_array',permu)

     end do
     !
     ! Construct edge arrays 
     !
     call memory_alloca(par_memor,'BEDGE_PERM','par_edge_communication_array',COMM % bedge_perm,COMM % bedge_dim)  
     kk = 0
     do ineig = 1,COMM % nneig 
        dom_i = COMM % neights(ineig)  
        do jj = 1,COMM % bedge_size(ineig)
           kk = kk + 1
           COMM % bedge_perm(kk) = list_neigh_edges(ineig) % l(jj)
        end do
     end do
     call graphs_number_to_linked_list(COMM % nneig,COMM % bedge_size)
     call memory_deallo(par_memor,'LIST_NEIGH_EDGES','par_edge_communication_array',list_neigh_edges)
     ! 
     ! Decide who owns the edge
     ! COUNTI = 2 ... both nodes are mine
     ! COUNTI = 1 ... one node is mine
     !
     do iedgg = 1,meshe % nedge     
        counti  = 0
        if( interface_edge(iedgg) ) then
           ipoin1 = meshe % edge_to_node(1,iedgg)
           ipoin2 = meshe % edge_to_node(2,iedgg)
           if( ipoin1 >= COMM % npoi2 .and. ipoin1 <= COMM % npoi3 ) counti = counti + 1
           if( ipoin2 >= COMM % npoi2 .and. ipoin2 <= COMM % npoi3 ) counti = counti + 1
           if( counti == 2 ) then
              edge_owner(iedgg) = huge(1_ip) 
           else
              edge_owner(iedgg) = kfl_paral
           end if
        end if
        edge_type(iedgg) = counti
     end do
     call PAR_INTERFACE_EDGE_EXCHANGE(edge_owner,'MAX','IN MY CODE')
     !call PAR_INTERFACE_EDGE_EXCHANGE(edge_owner,'TAKE MIN','IN MY CODE')

     do iedgg = 1,meshe % nedge     
        if( interface_edge(iedgg) ) then  

           counti = edge_type(iedgg)

           if(      counti == 2 ) then
              !
              ! Both nodes are mine
              !
              edge_type(iedgg) = OWN_BOUNDARY

           else if( edge_owner(iedgg) == kfl_paral ) then
              !
              ! I have one or none
              !
              edge_type(iedgg) = OWN_BOUNDARY

           else

              edge_type(iedgg) = OTH_BOUNDARY

           end if

        end if
     end do
     !
     ! Renumber edges: Interior-own boundary-other boundary
     ! NEW = INVPR(OLD) 
     ! OLD = PERMR(NEW)
     !
     call memory_alloca(par_memor,'PERMU','par_edge_communication_array',permu,meshe % nedge) 
     call memory_alloca(par_memor,'INVPR','par_edge_communication_array',invpr,meshe % nedge) 

     npoin_int = 0
     npoin_own = 0
     npoin_oth = 0
     do iedgg = 1,meshe % nedge
        if( edge_type(iedgg) == INTERIOR ) then
           npoin_int    = npoin_int + 1
           invpr(iedgg) = npoin_int
        else if( edge_type(iedgg) == OWN_BOUNDARY ) then
           npoin_own    = npoin_own + 1
           invpr(iedgg) = npoin_own 
        else
           npoin_oth    = npoin_oth + 1
           invpr(iedgg) = npoin_oth           
        end if
     end do
     do iedgg = 1,meshe % nedge
        if(      edge_type(iedgg) == OWN_BOUNDARY ) then
           invpr(iedgg) = invpr(iedgg) + npoin_int
        else if( edge_type(iedgg) == OTH_BOUNDARY ) then
           invpr(iedgg) = invpr(iedgg) + npoin_int + npoin_own
        end if
        iedge        = invpr(iedgg)
        permu(iedge) = iedgg
     end do
     do kk = 1,COMM % bedge_dim
        iedgg = COMM % bedge_perm(kk)
        COMM % bedge_perm(kk) = invpr(iedgg)  
     end do
     COMM % nedg1 = npoin_int
     COMM % nedg2 = npoin_int + 1
     COMM % nedg3 = npoin_int + npoin_own
     !
     ! Permute edge arrays LEDGS, LEDGB
     !
     do ielem = 1,meshe % nelem
        pedge = meshe % lnned(ielem)
        do iedge = 1,pedge
           iedgg = invpr(meshe % ledgs(iedge,ielem))
           meshe % ledgs(iedge,ielem) = iedgg
        end do
     end do
     do iboun = 1,meshe % nboun
        pedge = meshe % lnneb(iboun)
        do iedge = 1,pedge
           iedgg = invpr(meshe % ledgb(iedge,iboun))
           meshe % ledgb(iedge,iboun) = iedgg
        end do
     end do
     call maths_array_permutation(invpr,meshe % edge_to_node)
     !
     ! Permute edge graph: R_EDG, C_EDG
     !
     call graphs_rengra(meshe % nedge,0_ip,invpr,permu,meshe % r_edg,meshe % c_edg,memor=memor_dom)     
     !
     ! Check owner ship
     !
     permu = 0
     permu(1:COMM % nedg3) = 1
     call PAR_INTERFACE_EDGE_EXCHANGE(permu,'SUM','IN MY CODE')
     do iedgg = 1,meshe % nedge
        if( permu(iedgg) /= 1 ) then
           print*,'EDGE=',kfl_paral,edge_type(iedgg),edge_owner(iedgg),&
                meshe % lninv_loc(meshe % edge_to_node(1,iedgg)),&
                meshe % lninv_loc(meshe % edge_to_node(2,iedgg))
        end if
     end do
     !
     ! Construct multiplicty and owner rank arrays
     !
     call memory_alloca(par_memor,'BEDGE_OWNER_RANK'  ,'par_edge_communication_array',COMM % bedge_owner_rank  , meshe % nedge-COMM % nedg1)  
     call memory_alloca(par_memor,'BEDGE_MULTIPLICITY','par_edge_communication_array',COMM % bedge_multiplicity, meshe % nedge-COMM % nedg1)  
     call memory_alloca(par_memor,'IEXCHANGE'         ,'par_edge_communication_array',iexchange                , meshe % nedge )  

     iexchange(COMM % nedg2:COMM % nedg3) = kfl_paral
     call PAR_INTERFACE_EDGE_EXCHANGE(iexchange,'SUM','IN MY CODE')
     jj = 0
     do iedge = COMM % nedg1+1,meshe % nedge
        jj = jj + 1
        COMM % bedge_owner_rank(jj) = iexchange(iedge)
     end do

     iexchange(COMM % nedg1+1:meshe % nedge) = 1
     call PAR_INTERFACE_EDGE_EXCHANGE(iexchange,'SUM','IN MY CODE')
     jj = 0
     do iedge = COMM % nedg1+1,meshe % nedge
        jj = jj + 1
        COMM % bedge_multiplicity(jj) = iexchange(iedge)
     end do
     !
     ! Deallocate memory
     !
     call memory_deallo(par_memor,'IEXCHANGE'       ,'par_edge_communication_array',iexchange       )  
     call memory_deallo(par_memor,'PERMU'           ,'par_edge_communication_array',permu           ) 
     call memory_deallo(par_memor,'INVPR'           ,'par_edge_communication_array',invpr           )
     call memory_deallo(par_memor,'EDGE_TYPE'       ,'par_edge_communication_array',edge_type       )
     call memory_deallo(par_memor,'EDGE_OWNER'      ,'par_edge_communication_array',edge_owner      )
     call memory_deallo(par_memor,'INTERFACE EDGE'  ,'par_edge_communication_array',interface_edge  )
     call memory_deallo(par_memor,'NUMBER_EDGES'    ,'par_edge_communication_array',number_edges    )
     call memory_deallo(par_memor,'LIST_NEIGH_EDGES','par_edge_communication_array',list_neigh_edges)

  end if
  !
  ! Check errors
  !
  ierro = 0
  call memory_alloca(par_memor,'IEXCHANGE','par_edge_communication_array',iexchange,meshe % nedge)
  do iedge = 1,COMM % nedg3
     iexchange(iedge) = 1
  end do
  call PAR_INTERFACE_EDGE_EXCHANGE(iexchange,'SUM','IN MY CODE')
  do iedge = 1,meshe % nedge
     if( iexchange(iedge) /= 1 ) ierro = ierro + 1
  end do
  call memory_deallo(par_memor,'IEXCHANGE','par_edge_communication_array',iexchange)

  call PAR_SUM(ierro,'IN MY CODE')
  if( ierro > 0 ) then
     call runend('PAR_EDGE_COMMUNICATION_ARRAY: ERROR DETECTED ON EDGE OWNERSHIP: '//intost(ierro))
  end if
  !
  ! Copy to mesh type
  !
  nedg1         = COMM % nedg1
  nedg2         = COMM % nedg2
  nedg3         = COMM % nedg3
  meshe % nedg1 = nedg1
  meshe % nedg2 = nedg2
  meshe % nedg3 = nedg3
  
end subroutine par_edge_communication_array
