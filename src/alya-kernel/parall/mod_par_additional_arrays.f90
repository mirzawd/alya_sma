!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @name    Parallelization toolbox
!> @file    mod_par_additional_arrays.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!> @brief   "The easiest way of making software scalable is to make it sequentially inefficient."
!>          - Bill Gropp, 1999
!> @details ToolBox and definitions for parall.
!>
!>          \endverbatim
!>
!------------------------------------------------------------------------

module mod_par_additional_arrays
  use def_kintyp_basic,      only : ip,rp,lg,i1p
  use def_kintyp_comm,       only : comm_data_par
  use def_master,            only : routp,ioutp
  use def_master,            only : kfl_paral,ISLAVE,INOTMASTER,IPARALL
  use def_master,            only : INOTEMPTY
  use def_master,            only : IMASTER,lninv_loc
  use def_master,            only : npoi1,npoi2,npoi3
  use def_master,            only : kfl_paral
  use def_domain,            only : npoin,nelem,nboun
  use def_domain,            only : npoin_2,nelem_2,nboun_2
  use def_domain,            only : mesh_type,ndime 
  use def_domain,            only : memor_dom
  use def_kermod,            only : kfl_full_rows
  use mod_maths,             only : maths_heap_sort
  use mod_graphs,            only : graphs_number_to_linked_list
  use def_parall,            only : nneig
  use mod_parall,            only : commd  
  use mod_parall,            only : par_memor
  use mod_parall,            only : NODE_IN_NEIGHBOR
  use mod_par_tools,         only : PAR_GLOBAL_TO_LOCAL_NODE
  use mod_parall,            only : PAR_COMM_MY_CODE
  use mod_parall,            only : PAR_COMM_MY_CODE
  use mod_memory,            only : memory_alloca
  use mod_memory,            only : memory_deallo
  use mod_memory,            only : memory_resize
  use mod_memory,            only : memory_copy
  use mod_memory,            only : memory_size
  use mod_outfor,            only : outfor
  use mod_messages,          only : livinf
  use mod_optional_argument, only : optional_argument

  use mod_communications,    only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications,    only : PAR_INTERFACE_OWN_NODE_EXCHANGE
  use mod_communications,    only : PAR_GHOST_NODE_EXCHANGE
  use mod_communications,    only : PAR_SEND_RECEIVE
  use mod_communications,    only : PAR_MAX
  use mod_communications,    only : PAR_SUM

  use mod_communications,    only : PAR_COMM_RANK_AND_SIZE
  use mod_communications,    only : PAR_START_NON_BLOCKING_COMM
  use mod_communications,    only : PAR_SET_NON_BLOCKING_COMM_NUMBER
  use mod_communications,    only : PAR_END_NON_BLOCKING_COMM

  use def_mpi
#include "def_mpi.inc"
  use mod_std
  implicit none
  private

  public :: par_multiplicity_ownership
  public :: par_matrix_exchange_on_interface_nodes
  public :: par_matrix_w_halos_exchange_on_interface_nodes
  public :: par_matrix_computational_halo_exchange
  public :: par_node_number_in_owner
  public :: par_full_row_communication_arrays
  public :: par_ordered_exchange_update
  public :: par_global_variables_arrays
  public :: par_bounding_box
  
contains 

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-10-30
  !> @brief   Addition arrays
  !> @details Compute some useful arrays required by Alya just after the
  !>          Alya communicator has been computed
  !> 
  !-----------------------------------------------------------------------

  subroutine par_global_variables_arrays()
    !
    ! Compute some useful variables
    !
    if( IMASTER ) then
       npoi1   =  0
       npoi2   =  0
       npoi3   = -1
       npoin   =  0
       nelem   =  0
       nboun   =  0
       nneig   = -1
       npoin_2 = -1
       nelem_2 = -1
       nboun_2 = -1
    else
       npoin_2 = npoin
       nelem_2 = nelem
       nboun_2 = nboun
    end if
    
    commd % npoi1 = npoi1
    commd % npoi2 = npoi2
    commd % npoi3 = npoi3
    commd % npoin = npoin
    nneig         = commd % nneig
    !
    ! Compute ordered exchange update
    !
    call par_ordered_exchange_update(commd,COMM_NAME='COMMD')
    
  end subroutine par_global_variables_arrays
  
  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    06/03/2014
  !> @brief   multiplicty and ownership
  !> @details Construct multiplicty and owner rank arrays for node
  !>          communication. Edge multiplicity is computed in
  !>          par_edge_communication_array
  !>                          NPOI1+1  NPOIN         NPOIN_2
  !>          +-----------------+-------+--------------+
  !>          |                 |       |              |
  !>          +-----------------+-------+--------------+
  !>                            <----------------------> BOUND_OWNER_RANK
  !>                            <----------------------> BOUND_MULTIPLICITY
  !>
  !>          Important note: the information comming from ghost nodes
  !>          not necessarily comes from the subdomain who owns that node
  !>          as the owner can be a second-neighbor!
  !>
  !----------------------------------------------------------------------

  subroutine par_multiplicity_ownership(meshe,COMM,where,COMM_NAME)

    type(mesh_type),     intent(inout)          :: meshe         !< Mesh type
    type(comm_data_par), intent(inout)          :: COMM          !< Communication array
    character(*),        intent(in),   optional :: where         !< Where to exchange
    character(len=*),    intent(in),   optional :: COMM_NAME
    integer(ip)                                 :: ipoin,jj
    integer(ip),         pointer                :: iexchange(:)
    character(20)                               :: my_comm_name

    if( present(COMM_NAME) ) then
       my_comm_name = trim(COMM_NAME)
    else
       my_comm_name = 'COMM'
    end if
    
    if( ISLAVE .and. INOTEMPTY ) then

       nullify(iexchange)
       call memory_alloca(par_memor,trim(my_comm_name)//' % BOUND_OWNER_RANK'  ,'par_multiplicity_ownership',COMM % bound_owner_rank  , meshe % npoin_2-meshe % npoi1)
       call memory_alloca(par_memor,trim(my_comm_name)//' % BOUND_MULTIPLICITY','par_multiplicity_ownership',COMM % bound_multiplicity, meshe % npoin_2-meshe % npoi1)
       call memory_alloca(par_memor,'IEXCHANGE'                                ,'par_multiplicity_ownership',iexchange                , meshe % npoin_2 )

       iexchange(            1:meshe % npoi1) = kfl_paral
       iexchange(meshe % npoi2:meshe % npoi3) = kfl_paral
       if( present(where) ) then
          call PAR_INTERFACE_NODE_EXCHANGE(iexchange,'SUM'       ,trim(where))
          call PAR_GHOST_NODE_EXCHANGE    (iexchange,'SUBSTITUTE',trim(where))
       else
          call PAR_INTERFACE_NODE_EXCHANGE(iexchange,'SUM'       ,'IN MY CODE')
          call PAR_GHOST_NODE_EXCHANGE    (iexchange,'SUBSTITUTE','IN MY CODE')
       end if
       jj = 0
       do ipoin = meshe % npoi1+1,meshe % npoin_2
          jj = jj + 1
          COMM % bound_owner_rank(jj) = iexchange(ipoin)
       end do

       iexchange = 0
       iexchange(              1:meshe % npoi1) = 1
       iexchange(meshe % npoi1+1:meshe % npoin) = 1
       if( present(where) ) then
          call PAR_INTERFACE_NODE_EXCHANGE(iexchange,'SUM'       ,trim(where))
          call PAR_GHOST_NODE_EXCHANGE    (iexchange,'SUBSTITUTE',trim(where))
       else
          call PAR_INTERFACE_NODE_EXCHANGE(iexchange,'SUM'       ,'IN MY CODE')
          call PAR_GHOST_NODE_EXCHANGE    (iexchange,'SUBSTITUTE','IN MY CODE')
       end if
       jj = 0
       do ipoin = meshe % npoi1+1,meshe % npoin_2
          jj = jj + 1
          COMM % bound_multiplicity(jj) = iexchange(ipoin)
       end do

       call memory_deallo(par_memor,'IEXCHANGE','par_multiplicity_ownership',iexchange)

    end if

  end subroutine par_multiplicity_ownership

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    14/04/2016
  !> @brief   Interface node matrix exchange
  !> @details Construct communication arrays for matrix exchange
  !>          on interface nodes
  !>
  !----------------------------------------------------------------------

  subroutine par_matrix_exchange_on_interface_nodes(meshe,COMM,COMM_NAME)

    use mod_htable,           only  : HtableMaxPrimeNumber
    use mod_htable,           only  : hash_t
    use mod_htable,           only  : htaini
    use mod_htable,           only  : htaadd
    use mod_htable,           only  : htades
    use mod_htable,           only  : htalid
    use mod_htable,           only  : htares
    
    implicit none
    type(mesh_type),     intent(inout)          :: meshe         !< Mesh type
    type(comm_data_par), intent(inout)          :: COMM          !< Communication array
    character(len=*),    intent(in),   optional :: COMM_NAME
    character(20)                               :: my_comm_name
    integer(ip)                                 :: ipoin,ii,jj,bsize,ll
    integer(ip)                                 :: kpoin,kk,jpoin,dom_i,ineig,iz

    integer(ip),         pointer                :: list_nodes(:)
    integer(ip),         pointer                :: list_nodes_loc(:)
    !
    ! Hash table
    !
    type(hash_t)         :: ht
    
    if( present(COMM_NAME) ) then
       my_comm_name = trim(COMM_NAME)
    else
       my_comm_name = 'COMM'
    end if

    if( ISLAVE ) then

       nullify(list_nodes)
       nullify(list_nodes_loc)
       !
       ! Allocate memory for data type
       !
       if( associated(COMM % bound_matrix) ) then
          do ineig = 1,COMM % nneig
             call memory_deallo(par_memor,trim(my_comm_name)//' % BOUND_MATRIX % JA'      ,'par_matrix_exchange',COMM % bound_matrix(ineig) % ja)
             call memory_deallo(par_memor,trim(my_comm_name)//' % BOUND_MATRIX % NZDOM_II','par_matrix_exchange',COMM % bound_matrix(ineig) % nzdom_ii)
          end do
          deallocate(COMM % bound_matrix)
       end if
       allocate(COMM % bound_matrix(COMM % nneig))
       do ineig = 1,COMM % nneig
          nullify(COMM % bound_matrix(ineig) % ja)
          nullify(COMM % bound_matrix(ineig) % nzdom_ii)
          COMM % bound_matrix(ineig) % nzdom = 0 
       end do
       !
       ! INEIG = neighbor
       ! IPOIN = each interface node in common with neighbor
       ! JPOIN = neighbor connected to IPOIN
       !
       call htaini( ht, COMM % bound_size(COMM % nneig+1), lidson=.true., AUTOMATIC_SIZE=.true. )
       do ineig = 1,COMM % nneig
          !
          ! Restart htable
          !
          call htares( ht )
          do ii = COMM % bound_size(ineig),COMM % bound_size(ineig+1)-1
             call htaadd(ht, COMM % bound_perm(ii))
          end do 
          !
          ! Loop over subdomains INEIG
          !
          dom_i = COMM % neights(ineig)
          bsize = COMM % bound_size(ineig+1) - COMM % bound_size(ineig)   !

          call memory_alloca(par_memor,'LIST_NODES'                                  ,'par_matrix_exchange',list_nodes,bsize)
          call memory_alloca(par_memor,'LIST_NODES_LOC'                              ,'par_matrix_exchange',list_nodes_loc,bsize)
          call memory_alloca(par_memor,trim(my_comm_name)//' % BOUND_MATRIX % JA'      ,'par_matrix_exchange',COMM % bound_matrix(ineig) % ja,bsize)
          call memory_alloca(par_memor,trim(my_comm_name)//' % BOUND_MATRIX % NZDOM_II','par_matrix_exchange',COMM % bound_matrix(ineig) % nzdom_ii,bsize)

          jj = 0
          do ii = COMM % bound_size(ineig),COMM % bound_size(ineig+1)-1
             !
             ! Loop over the nodes IPOIN shared with INEIG
             !
             ipoin = COMM % bound_perm(ii)
             jj    = jj + 1
             kk    = 0
             do iz = meshe % r_dom(ipoin),meshe % r_dom(ipoin+1)-1
                jpoin = meshe % c_dom(iz)
                !
                ! JPOIN is also in INEIG?
                !
                !if( NODE_IN_NEIGHBOR(jpoin,ineig,COMM) ) then  
                if(  htalid(ht,jpoin) > 0_ip ) then
                   kk = kk + 1
                   list_nodes(kk) = jpoin
                end if
             end do
             !
             ! Fill in neighbors of IPOIN in INEIG: KK nodes common
             !
             if( kk > 0 ) then
                call memory_alloca(par_memor,trim(my_comm_name)//' % BOUND_MATRIX % JA % L','par_matrix_exchange',COMM % bound_matrix(ineig) % ja(jj) % l,kk)
                do ll = 1,kk
                   jpoin              = list_nodes(ll)
                   list_nodes_loc(ll) = meshe % lninv_loc(jpoin)
                   iz                 = meshe % r_dom(ipoin)
                   do while( iz <= meshe % r_dom(ipoin+1)-1 )
                      kpoin = meshe % c_dom(iz)
                      if( kpoin == jpoin ) then
                         COMM % bound_matrix(ineig) % ja(jj) % l(ll) = iz
                         iz = meshe % r_dom(ipoin+1)
                      end if
                      iz = iz + 1
                   end do
                   if( iz /= meshe % r_dom(ipoin+1)+1 ) call runend('PAR_MATRIX_EXCHANGE: COULD NOT FIND NODE')
                end do
                call maths_heap_sort(1_ip,kk,list_nodes_loc,' ',COMM % bound_matrix(ineig) % ja(jj) % l)
             else
                nullify(COMM % bound_matrix(ineig) % ja(jj) % l)
             end if
          end do
          !
          ! Count number of coefficients to send
          !
          bsize = 0
          do ii = 1,jj
             if( associated(COMM % bound_matrix(ineig) % ja(ii) % l) ) then
                COMM % bound_matrix(ineig) % nzdom_ii(ii) = size(COMM % bound_matrix(ineig) % ja(ii) % l)
                bsize = bsize + COMM % bound_matrix(ineig) % nzdom_ii(ii)
             end if
          end do
          COMM % bound_matrix(ineig) % nzdom = bsize

          call memory_deallo(par_memor,'LIST_NODES'    ,'par_matrix_exchange',list_nodes)
          call memory_deallo(par_memor,'LIST_NODES_LOC','par_matrix_exchange',list_nodes_loc)

       end do
       call htades(ht)
    
    end if

  end subroutine par_matrix_exchange_on_interface_nodes


  !----------------------------------------------------------------------
  !>
  !> @author  Herbert Owen (borrowing ideas from par_matrix_exchange_on_interface_nodes)
  !> @date    20/04/2016
  !> @brief   Interface node matrix exchange including halos
  !> @details Construct communication arrays for matrix exchange
  !>          All the row will be sent to subdomain who is the owner of that node
  !>
  !----------------------------------------------------------------------

  subroutine par_matrix_w_halos_exchange_on_interface_nodes(meshe,COMM)
    use def_domain,         only : npoin_2,nelem_2     ! I understand this will later come as part of meshe
    use mod_graphs,         only : graphs_poipoi
    use mod_communications, only : PAR_COMM_RANK_AND_SIZE
    use mod_communications, only : PAR_START_NON_BLOCKING_COMM
    use mod_communications, only : PAR_SET_NON_BLOCKING_COMM_NUMBER
    use mod_communications, only : PAR_END_NON_BLOCKING_COMM

    implicit none
    type(mesh_type),     intent(inout) :: meshe         !< Mesh type
    type(comm_data_par), intent(inout) :: COMM          !< Communication array
    integer(ip)                        :: ipoin,ii,jj,bsize,ll
    integer(ip)                        :: kk,jpoin,dom_i,ineig,iz
    integer(ip),parameter              :: idebug = 0

    type(i1p),           pointer       :: connect_in_glob_num_2send(:),connect_in_glob_num_2recv(:)


    ! por otro lado el ksize_connect2send esel equivalenet al COMM % bound_mat_halo_send(ineig) % nzdom
    ! ademas si voy a hacer todo de tamaño de la frontera - no de tamaño de nodos de frontera que pert al receive probablemenet KSIZE_NEDGE2SEND no lo necesite.

    ! con lo cual lo unico nuevo sea connect_in_glob_num_2send(:),connect_in_glob_num_2recv(:)
    ! ahoar el orden que tengo ahoar paar obtener als cosas es el correcto.

    !TO DO
    ! 1) rever si no hay problema que todo sea de tamaño de la frontera
    ! 2) cambiar nedge2send  por COMM % bound_mat_halo_send(ineig) % nzdom_ii y eliminar la aprte donde se volia a calcular
    ! 3) cambiar  ksize_connect2send    por    COMM % bound_mat_halo_send(ineig) % nzdom
    ! 4) ver si KSIZE_NEDGE2SEND me sirve paar algo sino eliminarlo
    ! 5) cambiar  COMM % bound_mat_halo_send(ineig) % ja(jj)    a COMM % bound_mat_halo_send(ineig) % ja(ii)                   por ejemplo en call memory_alloca(par_memor,'JA','par_matrix_exchange',COMM % bound_mat_halo_send(ineig) % ja(jj) % l,kk)

    integer(ip)                        :: iauxi,mm,kumul,ifoun

    integer(ip)                        :: PAR_CURRENT_RANK
    integer(ip)                        :: PAR_CURRENT_SIZE
    integer(ip)                        :: ssize,rsize
    !
    ! Some explanations:
    !
    ! COMM % bound_mat_halo_send(ineig) % nzdom_ii : This is the number of edges to send for each node in the shared boundary
    ! If the node does not belong to the subdomain that receives then nzdom_ii = 0. At some time I had called it nedge2send.
    ! it is a kind of ia but not in the way ia is typically stored. Just the number of connections for each node. Obvioslly it would be very easy to build an ia from it.
    !
    ! COMM % bound_mat_halo_send(ineig) % nzdom is the total. Adding for all nodes in the boundary
    !
    ! COMM % bound_mat_halo_send(ineig) % ja(ii) % l - is the vector that tells where to get from amatr the values it will send to the neighbor
    ! COMM % bound_mat_halo_recv(ineig) % ja(ii) % l - Once it has received them tells where to put them
    !
    ! connect_in_glob_num_2send(:),connect_in_glob_num_2recv(:) -
    ! is used to build bound_mat_halo_send & recive
    ! fora each point in the boundary (that is owned by the subdomain that receives) gives the nodes to which it is connected in global numbering
    !
    !
    ! bound_mat_halo_...  ja and nzdom_ii will be alloacted of size bsize (all the nodes in the boundary)
    ! those that do not belong to the subdomain that receives will have nzdom_ii = 0
    ! A 2nd option would have been to allocate it of a smaller size : only the nodes in the boundary that belong to the subdomain that receives.
    ! Actually the 2nd option was what I had until oct 16 in this subroutine but it was not consistent with what I had in PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP
    ! I decided to change this subroutine (1st option) because the second option would have implied needing to pass meshe to PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP
    ! and adding an additional if
    !
    !

    if( ISLAVE ) then

       !
       !  Obtain r_dom_2 & c_dom_2 to be used later
       !
       call graphs_poipoi(&
            npoin_2,nelem_2,meshe % mnode,meshe % lnods,meshe % lnnod,meshe % ltype,&
            meshe % r_dom_2,meshe % c_dom_2,&
            IA_NAME='MESHE % R_DOM_2',JA_NAME='MESHE % C_DOM_2',memor=memor_dom)
       ! hhh2check  . Cuando yo hablo de grafo extendido en principio solo necesitaria los conectados
       ! a mis own nodes de boundary si salen otros igual no me importará a AGMG le mandare solo los interior y boundary own y listo
       ! reordenados como corresponde

       !
       ! Allocate memory for data type
       !
       allocate(COMM % bound_mat_halo_send(COMM % nneig))
       allocate(COMM % bound_mat_halo_recv(COMM % nneig))
       do ineig = 1,COMM % nneig
          nullify(COMM % bound_mat_halo_send(ineig) % ja)
          nullify(COMM % bound_mat_halo_send(ineig) % nzdom_ii)
          COMM % bound_mat_halo_send(ineig) % nzdom = 0
          nullify(COMM % bound_mat_halo_recv(ineig) % ja)
          nullify(COMM % bound_mat_halo_recv(ineig) % nzdom_ii)
          COMM % bound_mat_halo_recv(ineig) % nzdom = 0
       end do
       !
       ! INEIG = neighbor
       ! IPOIN = each interface node in common with neighbor
       ! JPOIN = neighbor connected to IPOIN
       !
       !
       ! Obtain the number of conections (for each node) that I will send to my neighbor - set to 0 if the node does not belong to my neighbor
       !
       do ineig = 1,comm % nneig
          bsize = COMM % bound_size(ineig+1) - COMM % bound_size(ineig)
          call memory_alloca(par_memor,'COMM % BOUND_MAT_HALO_SEND % NZDOM_II','par_matrix_w_halos_exchange',COMM % bound_mat_halo_send(ineig) % nzdom_ii,bsize)
          call memory_alloca(par_memor,'COMM % BOUND_MAT_HALO_RECV % NZDOM_II','par_matrix_w_halos_exchange',COMM % bound_mat_halo_recv(ineig) % nzdom_ii,bsize)
       end do

       do ineig = 1,COMM % nneig
          !
          ! Loop over subdomains INEIG
          !
!          dom_i = COMM % neights(ineig)
          bsize = COMM % bound_size(ineig+1) - COMM % bound_size(ineig)
          !
          do ii = COMM % bound_size(ineig),COMM % bound_size(ineig+1)-1   ! Just to count
             jj = ii - COMM % bound_size(ineig) + 1
             !
             ! What I have to send? - Loop over the nodes IPOIN shared with INEIG and see if ineig is the owner of that node
             !
             ipoin = COMM % bound_perm(ii)
             kk = 0
             ! Does the point belong to this neighbour so that I have to send him the row
             if ( COMM % bound_owner_rank(ipoin - meshe % npoi1 ) == COMM % neights(ineig) ) then
                iauxi = meshe % r_dom(ipoin+1) - meshe % r_dom(ipoin)
                COMM % bound_mat_halo_send(ineig) % nzdom_ii(jj) = iauxi
             else
                COMM % bound_mat_halo_send(ineig) % nzdom_ii(jj) = 0
             end if
          end do
       end do
       !
       ! Send & receive number of connections (for each node)
       !
       call PAR_COMM_RANK_AND_SIZE(COMM % PAR_COMM_WORLD,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
       call PAR_START_NON_BLOCKING_COMM(1_ip,PAR_CURRENT_SIZE)
       call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
       do ineig = 1,comm % nneig
          dom_i = comm % neights(ineig)
          bsize = COMM % bound_size(ineig+1) - COMM % bound_size(ineig)
          if( bsize > 0 ) call PAR_SEND_RECEIVE(bsize,bsize,COMM % bound_mat_halo_send(ineig) % nzdom_ii,&
               COMM % bound_mat_halo_recv(ineig) % nzdom_ii,'IN MY CODE',dom_i,'NON BLOCKING')
       end do
       call PAR_END_NON_BLOCKING_COMM(1_ip)
       !
       ! Allocate CONNECT_IN_GLOB_NUM_2SEND & 2RECV - also COMM % bound_mat_halo_send(ineig) % ja  , % ja % l and  the same for recv
       !
       nullify(connect_in_glob_num_2send)
       nullify(connect_in_glob_num_2recv)
       call memory_alloca(par_memor,'CONNECT_IN_GLOB_NUM_2SEND','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',connect_in_glob_num_2send,comm % nneig)
       call memory_alloca(par_memor,'CONNECT_IN_GLOB_NUM_2RECV','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',connect_in_glob_num_2recv,comm % nneig)

       do ineig = 1,comm % nneig
          bsize = COMM % bound_size(ineig+1) - COMM % bound_size(ineig)
          !
          ! Allocate bound_mat_halo_.. % ja  & %ja % l
          !
          call memory_alloca(par_memor,'COMM % BOUND_MAT_HALO_SEND % JA','par_matrix_w_halos_exchange',COMM % bound_mat_halo_send(ineig) % ja,bsize)
          call memory_alloca(par_memor,'COMM % BOUND_MAT_HALO_SEND % JA','par_matrix_w_halos_exchange',COMM % bound_mat_halo_recv(ineig) % ja,bsize)

          COMM % bound_mat_halo_send(ineig) % nzdom = 0  ! send size
          COMM % bound_mat_halo_recv(ineig) % nzdom = 0  ! recv size
          do ii = 1, bsize
             COMM % bound_mat_halo_send(ineig) % nzdom = COMM % bound_mat_halo_send(ineig) % nzdom + COMM % bound_mat_halo_send(ineig) % nzdom_ii (ii)
             COMM % bound_mat_halo_recv(ineig) % nzdom = COMM % bound_mat_halo_recv(ineig) % nzdom + COMM % bound_mat_halo_recv(ineig) % nzdom_ii (ii)

             nullify(COMM % bound_mat_halo_send(ineig) % ja(ii) % l)
             nullify(COMM % bound_mat_halo_recv(ineig) % ja(ii) % l)
             if ( COMM % bound_mat_halo_send(ineig) % nzdom_ii (ii) /= 0 ) &
                  call memory_alloca(par_memor,'COMM % BOUND_MAT_HALO_SEND % NZDOM_II','par_matrix_w_halos_exchange',  &
                  COMM % bound_mat_halo_send(ineig) % ja(ii) % l,COMM % bound_mat_halo_send(ineig) % nzdom_ii (ii)  )
             if ( COMM % bound_mat_halo_recv(ineig) % nzdom_ii (ii) /= 0 ) &
                  call memory_alloca(par_memor,'COMM % BOUND_MAT_HALO_RECV % NZDOM_II','par_matrix_w_halos_exchange',  &
                COMM % bound_mat_halo_recv(ineig) % ja(ii) % l,COMM % bound_mat_halo_recv(ineig) % nzdom_ii (ii)  )
          end do
          !
          ! I had to add max else in the send recv it would complaint it is not allocated
          !
          call memory_alloca(par_memor,'CONNECT_IN_GLOB_NUM_2SEND','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',&
               connect_in_glob_num_2send(ineig) % l,max(1_ip,COMM % bound_mat_halo_send(ineig) % nzdom))
          call memory_alloca(par_memor,'CONNECT_IN_GLOB_NUM_2RECV','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',&
               connect_in_glob_num_2recv(ineig) % l,max(1_ip,COMM % bound_mat_halo_recv(ineig) % nzdom))
       end do
       !
       ! Obtain connect_in_glob_num_2send  & COMM % bound_mat_halo_send
       !
       do ineig = 1,COMM % nneig
          !
          ! Loop over subdomains INEIG
          !
          dom_i = COMM % neights(ineig)
          bsize = COMM % bound_size(ineig+1) - COMM % bound_size(ineig)

          if (idebug==1) write(400+kfl_paral,*)ineig,kfl_paral,COMM % neights(ineig)

          mm = 0
          do ii = COMM % bound_size(ineig),COMM % bound_size(ineig+1)-1
             !
             ! Loop over the nodes IPOIN shared with INEIG
             !
             ipoin = COMM % bound_perm(ii)
             jj = ii - COMM % bound_size(ineig) + 1
             !
             ! Obtain values
             !
             if (idebug==1) write (kfl_paral+510,'(20(i9,1x))') ineig,jj,ipoin,meshe % lninv_loc(ipoin),&
                  COMM % bound_owner_rank(ipoin - meshe % npoi1 ) , COMM % neights(ineig)
             if ( COMM % bound_owner_rank(ipoin - meshe % npoi1 ) == COMM % neights(ineig) ) then
                ll = 0
                do iz = meshe % r_dom(ipoin),meshe % r_dom(ipoin+1)-1
                   mm = mm + 1
                   ll = ll + 1
                   jpoin  = meshe % c_dom(iz)
                   connect_in_glob_num_2send(ineig) % l(mm) = meshe % lninv_loc(jpoin)
                   COMM % bound_mat_halo_send(ineig) % ja(jj) % l(ll) = iz
                   if (idebug==1) write (kfl_paral+500,'(20(i9,1x))') ineig,jj,ipoin,meshe % lninv_loc(ipoin),jpoin,meshe % lninv_loc(jpoin),mm,connect_in_glob_num_2send(ineig) % l(mm)
                end do
             end if
          end do
       end do
       !
       ! Send & receive connections (for each node)
       !
       call PAR_COMM_RANK_AND_SIZE(COMM % PAR_COMM_WORLD,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
       call PAR_START_NON_BLOCKING_COMM(1_ip,PAR_CURRENT_SIZE)
       call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
       do ineig = 1,comm % nneig
          dom_i = comm % neights(ineig)
          ssize = COMM % bound_mat_halo_send(ineig) % nzdom
          rsize = COMM % bound_mat_halo_recv(ineig) % nzdom
          if( ssize > 0 .or. rsize > 0 ) call PAR_SEND_RECEIVE(ssize,rsize,connect_in_glob_num_2send(ineig) % l, &
               connect_in_glob_num_2recv(ineig) % l,'IN MY CODE',dom_i,'NON BLOCKING')
       end do
       call PAR_END_NON_BLOCKING_COMM(1_ip)

       !
       ! Obtain COMM % bound_mat_halo_recv
       !
       do ineig = 1,comm % nneig
          kumul = 0
          do ii = COMM % bound_size(ineig),COMM % bound_size(ineig+1)-1
             jj = ii - COMM % bound_size(ineig) + 1
             !
             ! Loop over the nodes IPOIN shared with INEIG
             !
             if (COMM % bound_mat_halo_recv(ineig) % nzdom_ii (jj) /= 0 ) then
                ipoin = COMM % bound_perm(ii)
                do kk = kumul+1, kumul + COMM % bound_mat_halo_recv(ineig) % nzdom_ii(jj)    !  if the node is not mine  nzdom_ii(ii) = 0 and nothing is done
                   ifoun=0
                   all_con: do iz = meshe % r_dom_2(ipoin),meshe % r_dom_2(ipoin+1)-1
                      jpoin = meshe % lninv_loc(meshe % c_dom_2(iz))
                      if(connect_in_glob_num_2recv(ineig) % l(kk) == jpoin) then
                         ifoun=1
                         COMM % bound_mat_halo_recv(ineig) % ja(jj) % l(kk-kumul) = iz  !REVISAR si es ii o jj!!!!!
                      end if
                      if (idebug==1) write (kfl_paral+600,'(20(i9,1x))') ineig,jj,ipoin,meshe % lninv_loc(ipoin),jpoin,kk,connect_in_glob_num_2recv(ineig) % l(kk)
                      if (ifoun==1) exit all_con
                   end do all_con
                   if (ifoun==0) write (kfl_paral+600,*) 'falloooo'
!                   if (ifoun==0) call runend('par_matrix_w_halos_exchange_on_interface_nodes:connectivy not found')
                end do
                kumul = kumul + COMM % bound_mat_halo_recv(ineig) % nzdom_ii(jj)
             end if
          end do
       end do


       do ineig = 1,comm % nneig
          call memory_deallo(par_memor,'CONNECT_IN_GLOB_NUM_2SEND','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',&
               connect_in_glob_num_2send(ineig) % l)
          call memory_deallo(par_memor,'CONNECT_IN_GLOB_NUM_2RECV','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',&
               connect_in_glob_num_2recv(ineig) % l)
       end do
       call memory_deallo(par_memor,'CONNECT_IN_GLOB_NUM_2SEND','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',connect_in_glob_num_2send)
       call memory_deallo(par_memor,'CONNECT_IN_GLOB_NUM_2RECV','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',connect_in_glob_num_2recv)

    end if

  end subroutine par_matrix_w_halos_exchange_on_interface_nodes





  subroutine par_matrix_computational_halo_exchange(meshe,COMM,COMM_NAME)
    use mod_communications, only : PAR_COMM_RANK_AND_SIZE
    use mod_communications, only : PAR_START_NON_BLOCKING_COMM
    use mod_communications, only : PAR_SET_NON_BLOCKING_COMM_NUMBER
    use mod_communications, only : PAR_END_NON_BLOCKING_COMM

    implicit none
    type(mesh_type),     intent(inout)          :: meshe         !< Mesh type
    type(comm_data_par), intent(inout)          :: COMM          !< Communication array
    character(len=*),    intent(in),   optional :: COMM_NAME
    character(20)                               :: my_comm_name
    integer(ip)                                 :: ipoin,ii,jj,bsize,ll
    integer(ip)                                 :: kk,jpoin,dom_i,ineig,iz
    integer(ip),parameter                       :: idebug = 0

    type(i1p),           pointer                :: connect_in_glob_num_2send(:),connect_in_glob_num_2recv(:)


    ! por otro lado el ksize_connect2send esel equivalenet al COMM % bound_mat_halo_send(ineig) % nzdom
    ! ademas si voy a hacer todo de tamaño de la frontera - no de tamaño de nodos de frontera que pert al receive probablemenet KSIZE_NEDGE2SEND no lo necesite.

    ! con lo cual lo unico nuevo sea connect_in_glob_num_2send(:),connect_in_glob_num_2recv(:)
    ! ahoar el orden que tengo ahoar paar obtener als cosas es el correcto.

    !TO DO
    ! 1) rever si no hay problema que todo sea de tamaño de la frontera
    ! 2) cambiar nedge2send  por COMM % bound_mat_halo_send(ineig) % nzdom_ii y eliminar la aprte donde se volia a calcular
    ! 3) cambiar  ksize_connect2send    por    COMM % bound_mat_halo_send(ineig) % nzdom
    ! 4) ver si KSIZE_NEDGE2SEND me sirve paar algo sino eliminarlo
    ! 5) cambiar  COMM % bound_mat_halo_send(ineig) % ja(jj)    a COMM % bound_mat_halo_send(ineig) % ja(ii)                   por ejemplo en call memory_alloca(par_memor,'JA','par_matrix_exchange',COMM % bound_mat_halo_send(ineig) % ja(jj) % l,kk)

    integer(ip)                        :: iauxi,mm,kumul,ifoun

    integer(ip)                        :: PAR_CURRENT_RANK
    integer(ip)                        :: PAR_CURRENT_SIZE
    integer(ip)                        :: ssize,rsize

    !
    ! Some explanations:
    !
    ! COMM % bound_mat_halo_send(ineig) % nzdom_ii : This is the number of edges to send for each node in the shared boundary
    ! If the node does not belong to the subdomain that receives then nzdom_ii = 0. At some time I had called it nedge2send.
    ! it is a kind of ia but not in the way ia is typically stored. Just the number of connections for each node. Obvioslly it would be very easy to build an ia from it.
    !
    ! COMM % bound_mat_halo_send(ineig) % nzdom is the total. Adding for all nodes in the boundary
    !
    ! COMM % bound_mat_halo_send(ineig) % ja(ii) % l - is the vector that tells where to get from amatr the values it will send to the neighbor
    ! COMM % bound_mat_halo_recv(ineig) % ja(ii) % l - Once it has received them tells where to put them
    !
    ! connect_in_glob_num_2send(:),connect_in_glob_num_2recv(:) -
    ! is used to build bound_mat_halo_send & recive
    ! fora each point in the boundary (that is owned by the subdomain that receives) gives the nodes to which it is connected in global numbering
    !
    !
    ! bound_mat_halo_...  ja and nzdom_ii will be alloacted of size bsize (all the nodes in the boundary)
    ! those that do not belong to the subdomain that receives will have nzdom_ii = 0
    ! A 2nd option would have been to allocate it of a smaller size : only the nodes in the boundary that belong to the subdomain that receives.
    ! Actually the 2nd option was what I had until oct 16 in this subroutine but it was not consistent with what I had in PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP
    ! I decided to change this subroutine (1st option) because the second option would have implied needing to pass meshe to PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP
    ! and adding an additional if
    !
    !

    if( kfl_full_rows == 0 ) return
    
    if( present(COMM_NAME) ) then
       my_comm_name = trim(COMM_NAME)
    else
       my_comm_name = 'COMM'
    end if
    
    if( ISLAVE ) then
       !
       ! Allocate memory for data type
       !
       if( associated(COMM % bound_mat_halo_send) ) deallocate(COMM % bound_mat_halo_send)
       if( associated(COMM % bound_mat_halo_recv) ) deallocate(COMM % bound_mat_halo_recv)
       allocate(COMM % bound_mat_halo_send(COMM % nneig))
       allocate(COMM % bound_mat_halo_recv(COMM % nneig))
       do ineig = 1,COMM % nneig
          nullify(COMM % bound_mat_halo_send(ineig) % ja)
          nullify(COMM % bound_mat_halo_send(ineig) % nzdom_ii)
          COMM % bound_mat_halo_send(ineig) % nzdom = 0
          nullify(COMM % bound_mat_halo_recv(ineig) % ja)
          nullify(COMM % bound_mat_halo_recv(ineig) % nzdom_ii)
          COMM % bound_mat_halo_recv(ineig) % nzdom = 0
       end do
       !
       ! INEIG = neighbor
       ! IPOIN = each interface node in common with neighbor
       ! JPOIN = neighbor connected to IPOIN
       !
       !
       ! Obtain the number of conections (for each node) that I will send to my neighbor - set to 0 if the node does not belong to my neighbor
       !
       do ineig = 1,comm % nneig
          bsize = COMM % bound_size(ineig+1) - COMM % bound_size(ineig)
          call memory_alloca(par_memor,trim(my_comm_name)//' % BOUND_MATH_HALO_SEND % NZDOM','par_matrix_computational_halo_exchange',COMM % bound_mat_halo_send(ineig) % nzdom_ii,bsize)
          call memory_alloca(par_memor,trim(my_comm_name)//' % BOUND_MATH_HALO_RECV % NZDOM','par_matrix_computational_halo_exchange',COMM % bound_mat_halo_recv(ineig) % nzdom_ii,bsize)
       end do

       do ineig = 1,COMM % nneig
          !
          ! Loop over subdomains INEIG
          !
!          dom_i = COMM % neights(ineig)
          bsize = COMM % bound_size(ineig+1) - COMM % bound_size(ineig)
          !
          do ii = COMM % bound_size(ineig),COMM % bound_size(ineig+1)-1   ! Just to count
             jj = ii - COMM % bound_size(ineig) + 1
             !
             ! What I have to send? - Loop over the nodes IPOIN shared with INEIG and see if ineig is the owner of that node
             !
             ipoin = COMM % bound_perm(ii)
             kk = 0
             ! Does the point belong to this neighbour so that I have to send him the row
             if ( COMM % bound_owner_rank(ipoin - meshe % npoi1 ) == COMM % neights(ineig) ) then
                iauxi = meshe % r_dom(ipoin+1) - meshe % r_dom(ipoin)
                COMM % bound_mat_halo_send(ineig) % nzdom_ii(jj) = iauxi
             else
                COMM % bound_mat_halo_send(ineig) % nzdom_ii(jj) = 0
             end if
          end do
       end do
       !
       ! Send & receive number of connections (for each node)
       !
       call PAR_COMM_RANK_AND_SIZE(COMM % PAR_COMM_WORLD,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
       call PAR_START_NON_BLOCKING_COMM(1_ip,PAR_CURRENT_SIZE)
       call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
       do ineig = 1,comm % nneig
          dom_i = comm % neights(ineig)
          bsize = COMM % bound_size(ineig+1) - COMM % bound_size(ineig)
          if( bsize > 0 ) call PAR_SEND_RECEIVE(bsize,bsize,COMM % bound_mat_halo_send(ineig) % nzdom_ii,&
               COMM % bound_mat_halo_recv(ineig) % nzdom_ii,'IN MY CODE',dom_i,'NON BLOCKING')
       end do
       call PAR_END_NON_BLOCKING_COMM(1_ip)
       !
       ! Allocate CONNECT_IN_GLOB_NUM_2SEND & 2RECV - also COMM % bound_mat_halo_send(ineig) % ja  , % ja % l and  the same for recv
       !
       nullify(connect_in_glob_num_2send)
       nullify(connect_in_glob_num_2recv)
       call memory_alloca(par_memor,'CONNECT_IN_GLOB_NUM_2SEND','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',connect_in_glob_num_2send,comm % nneig)
       call memory_alloca(par_memor,'CONNECT_IN_GLOB_NUM_2RECV','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',connect_in_glob_num_2recv,comm % nneig)

       do ineig = 1,comm % nneig
          bsize = COMM % bound_size(ineig+1) - COMM % bound_size(ineig)
          !
          ! Allocate bound_mat_halo_.. % ja  & %ja % l
          !
          call memory_alloca(par_memor,trim(my_comm_name)//' % BOUND_MAT_HALO_SEND % JA','par_matrix_w_halos_exchange',COMM % bound_mat_halo_send(ineig) % ja,bsize)
          call memory_alloca(par_memor,trim(my_comm_name)//' % BOUND_MAT_HALO_RECV % JA','par_matrix_w_halos_exchange',COMM % bound_mat_halo_recv(ineig) % ja,bsize)

          COMM % bound_mat_halo_send(ineig) % nzdom = 0  ! send size
          COMM % bound_mat_halo_recv(ineig) % nzdom = 0  ! recv size
          do ii = 1, bsize
             COMM % bound_mat_halo_send(ineig) % nzdom = COMM % bound_mat_halo_send(ineig) % nzdom + COMM % bound_mat_halo_send(ineig) % nzdom_ii (ii)
             COMM % bound_mat_halo_recv(ineig) % nzdom = COMM % bound_mat_halo_recv(ineig) % nzdom + COMM % bound_mat_halo_recv(ineig) % nzdom_ii (ii)

             nullify(COMM % bound_mat_halo_send(ineig) % ja(ii) % l)
             nullify(COMM % bound_mat_halo_recv(ineig) % ja(ii) % l)
             if ( COMM % bound_mat_halo_send(ineig) % nzdom_ii (ii) /= 0 ) call memory_alloca(par_memor,trim(my_comm_name)//' % BOUND_MAT_HALO_SEND % JA % L','par_matrix_computational_halo_exchange',&
                COMM % bound_mat_halo_send(ineig) % ja(ii) % l,COMM % bound_mat_halo_send(ineig) % nzdom_ii (ii)  )
             if ( COMM % bound_mat_halo_recv(ineig) % nzdom_ii (ii) /= 0 ) call memory_alloca(par_memor,trim(my_comm_name)//' % BOUND_MAT_HALO_RECV % JA % L','par_matrix_computational_halo_exchange',&
                COMM % bound_mat_halo_recv(ineig) % ja(ii) % l,COMM % bound_mat_halo_recv(ineig) % nzdom_ii (ii)  )
          end do
          !
          ! I had to add max else in the send recv it would complaint it is not allocated
          !
          call memory_alloca(par_memor,'CONNECT_IN_GLOB_NUM_2SEND','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',&
               connect_in_glob_num_2send(ineig) % l,max(1_ip,COMM % bound_mat_halo_send(ineig) % nzdom))
          call memory_alloca(par_memor,'CONNECT_IN_GLOB_NUM_2RECV','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',&
               connect_in_glob_num_2recv(ineig) % l,max(1_ip,COMM % bound_mat_halo_recv(ineig) % nzdom))
       end do
       !
       ! Obtain connect_in_glob_num_2send  & COMM % bound_mat_halo_send
       !
       do ineig = 1,COMM % nneig
          !
          ! Loop over subdomains INEIG
          !
          dom_i = COMM % neights(ineig)
          bsize = COMM % bound_size(ineig+1) - COMM % bound_size(ineig)

          if (idebug==1) write(400+kfl_paral,*)ineig,kfl_paral,COMM % neights(ineig)

          mm = 0
          do ii = COMM % bound_size(ineig),COMM % bound_size(ineig+1)-1
             !
             ! Loop over the nodes IPOIN shared with INEIG
             !
             ipoin = COMM % bound_perm(ii)
             jj = ii - COMM % bound_size(ineig) + 1
             !
             ! Obtain values
             !
             if (idebug==1) write (kfl_paral+510,'(20(i9,1x))') ineig,jj,ipoin,meshe % lninv_loc(ipoin),&
                  COMM % bound_owner_rank(ipoin - meshe % npoi1 ) , COMM % neights(ineig)
             if ( COMM % bound_owner_rank(ipoin - meshe % npoi1 ) == COMM % neights(ineig) ) then
                ll = 0
                do iz = meshe % r_dom(ipoin),meshe % r_dom(ipoin+1)-1
                   mm = mm + 1
                   ll = ll + 1
                   jpoin  = meshe % c_dom(iz)
                   connect_in_glob_num_2send(ineig) % l(mm) = meshe % lninv_loc(jpoin)
                   COMM % bound_mat_halo_send(ineig) % ja(jj) % l(ll) = iz
                   if (idebug==1) write (kfl_paral+500,'(20(i9,1x))') ineig,jj,ipoin,meshe % lninv_loc(ipoin),jpoin,meshe % lninv_loc(jpoin),mm,connect_in_glob_num_2send(ineig) % l(mm)
                end do
             end if
          end do
       end do
       !
       ! Send & receive connections (for each node)
       !
       call PAR_COMM_RANK_AND_SIZE(COMM % PAR_COMM_WORLD,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
       call PAR_START_NON_BLOCKING_COMM(1_ip,PAR_CURRENT_SIZE)
       call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
       do ineig = 1,comm % nneig
          dom_i = comm % neights(ineig)
          ssize = COMM % bound_mat_halo_send(ineig) % nzdom
          rsize = COMM % bound_mat_halo_recv(ineig) % nzdom
          if( ssize > 0 .or. rsize > 0 ) call PAR_SEND_RECEIVE(ssize,rsize,connect_in_glob_num_2send(ineig) % l, &
               connect_in_glob_num_2recv(ineig) % l,'IN MY CODE',dom_i,'NON BLOCKING')
       end do
       call PAR_END_NON_BLOCKING_COMM(1_ip)

       !
       ! Obtain COMM % bound_mat_halo_recv
       !
       do ineig = 1,comm % nneig
          kumul = 0
          do ii = COMM % bound_size(ineig),COMM % bound_size(ineig+1)-1
             jj = ii - COMM % bound_size(ineig) + 1
             !
             ! Loop over the nodes IPOIN shared with INEIG
             !
             if (COMM % bound_mat_halo_recv(ineig) % nzdom_ii (jj) /= 0 ) then
                ipoin = COMM % bound_perm(ii)
                do kk = kumul+1, kumul + COMM % bound_mat_halo_recv(ineig) % nzdom_ii(jj)    !  if the node is not mine  nzdom_ii(ii) = 0 and nothing is done
                   ifoun=0
                   all_con: do iz = meshe % r_dom_own(ipoin),meshe % r_dom_own(ipoin+1)-1
                      jpoin = meshe % lninv_loc(meshe % c_dom_own(iz))
                      if(connect_in_glob_num_2recv(ineig) % l(kk) == jpoin) then
                         ifoun=1
                         COMM % bound_mat_halo_recv(ineig) % ja(jj) % l(kk-kumul) = iz  !REVISAR si es ii o jj!!!!!
                      end if
                      if (idebug==1) write (kfl_paral+600,'(20(i9,1x))') ineig,jj,ipoin,meshe % lninv_loc(ipoin),jpoin,kk,connect_in_glob_num_2recv(ineig) % l(kk)
                      if (ifoun==1) exit all_con
                   end do all_con
                   if (ifoun==0) write (kfl_paral+600,*) 'falloooo'
!                   if (ifoun==0) call runend('par_matrix_w_halos_exchange_on_interface_nodes:connectivy not found')
                end do
                kumul = kumul + COMM % bound_mat_halo_recv(ineig) % nzdom_ii(jj)
             end if
          end do
       end do


       do ineig = 1,comm % nneig
          call memory_deallo(par_memor,'CONNECT_IN_GLOB_NUM_2SEND','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',&
               connect_in_glob_num_2send(ineig) % l)
          call memory_deallo(par_memor,'CONNECT_IN_GLOB_NUM_2RECV','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',&
               connect_in_glob_num_2recv(ineig) % l)
       end do
       call memory_deallo(par_memor,'CONNECT_IN_GLOB_NUM_2SEND','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',connect_in_glob_num_2send)
       call memory_deallo(par_memor,'CONNECT_IN_GLOB_NUM_2RECV','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',connect_in_glob_num_2recv)

    end if

  end subroutine par_matrix_computational_halo_exchange


  !----------------------------------------------------------------------
  !>
  !> @author  Herbert Owen
  !> @date    10/06/2016
  !> @brief   Obtains the numeration used in th subdomain that is the owner of the poin
  !> @details Ideas borrowed from par_multiplicity_ownership
  !>
  !----------------------------------------------------------------------

  subroutine par_node_number_in_owner(meshe,COMM,COMM_NAME)

    type(mesh_type),     intent(inout)          :: meshe         !< Mesh type
    type(comm_data_par), intent(inout)          :: COMM          !< Communication array
    character(len=*),    intent(in),   optional :: COMM_NAME
    character(20)                               :: my_comm_name
    integer(ip)                                 :: ipoin
    integer(ip),         pointer                :: iexchange(:)

    if( present(COMM_NAME) ) then
       my_comm_name = trim(COMM_NAME)
    else
       my_comm_name = 'COMM'
    end if
    
    if( ISLAVE .and. meshe % npoin_2 > 0 ) then

       nullify(iexchange)
       call memory_alloca(par_memor, trim(my_comm_name)//' % NODE_NUMBER_IN_OWNER', 'par_multiplicity_ownership', COMM % node_number_in_owner  ,&
            meshe % npoin_2 - meshe % npoi1)
       call memory_alloca(par_memor, 'IEXCHANGE',            'par_multiplicity_ownership', iexchange, meshe % npoin_2 )

       iexchange = 0
       do ipoin=1, meshe % npoi1
          iexchange(ipoin) = ipoin
       end do
       do ipoin = meshe % npoi2, meshe % npoi3
          iexchange(ipoin) = ipoin
       end do
       call PAR_INTERFACE_NODE_EXCHANGE(iexchange,'SUM'       ,'IN MY CODE')
       call PAR_GHOST_NODE_EXCHANGE    (iexchange,'SUBSTITUTE','IN MY CODE')
       do ipoin = meshe % npoi1+1, meshe % npoin_2
          COMM % node_number_in_owner(ipoin - meshe % npoi1) = iexchange(ipoin)
       end do

       call memory_deallo(par_memor,'IEXCHANGE','par_multiplicity_ownership',iexchange)

    end if

  end subroutine par_node_number_in_owner

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    27/03/2017
  !> @brief   Communicator for exchanging halos
  !> @details Communicator for exchanging halos
  !>
  !----------------------------------------------------------------------

  subroutine par_full_row_communication_arrays(meshe,COMM,COMM_NAME)

!    use def_domain,           only  : htable_lninv_loc
    use def_domain,           only  : lmast,nperi
    use mod_htable,           only  : HtableMaxPrimeNumber
    use mod_htable,           only  : hash_t
    use mod_htable,           only  : htaadd
    use mod_htable,           only  : htades
    use mod_htable,           only  : htalid

    implicit none

    type(mesh_type),     intent(inout)          :: meshe            !< Mesh type
    type(comm_data_par), intent(inout)          :: COMM             !< Communication array
    character(len=*),    intent(in),   optional :: COMM_NAME
    integer(ip)                                 :: ipoin,ineig,ii
    integer(ip)                                 :: jj,dom_i
    integer(ip)                                 :: ii_send,ii_send_size
    integer(ip)                                 :: ii_recv,ii_recv_size
    integer(ip)                                 :: kneig,jneig,ipoin_global
    integer(ip)                                 :: ierro
    integer(ip)                                 :: num_second,num_first

    integer(ip)                                 :: nneig_second_send
    integer(ip)                                 :: nneig_second_recv
    integer(ip)                                 :: nneig_second
    integer(ip),         pointer                :: lneig(:)
    integer(ip),         pointer                :: lneig_second_send(:)
    type(i1p),           pointer                :: lneig_second_recv(:)
    integer(ip),         pointer                :: lneig_second(:)
    integer(ip),         pointer                :: full_row_send_size(:)
    integer(ip),         pointer                :: full_row_recv_neights(:)

    integer(ip),         target                 :: full_row_send_perm_ok(2)
    integer(ip),         target                 :: full_row_recv_perm_ok(2)

    integer(ip)                                 :: PAR_CURRENT_SIZE
    integer(ip)                                 :: PAR_CURRENT_RANK
    character(20)                               :: my_comm_name
    ! Check everything's ok
    integer(ip),         pointer       :: lcheck(:)
    real(rp)                           :: ratio_second_total
    logical(lg)                        :: debug
    
    if( nperi <= -1 ) return
    
    if( present(COMM_NAME) ) then
       my_comm_name = trim(COMM_NAME)
    else
       my_comm_name = 'COMM'
    end if

    call livinf(0_ip,'PARALL: COMPUTE FULL ROW COMMUNICATOR FOR COMPUTATIONAL HALO EXCHANGE',0_ip)

    nullify(lneig)
    nullify(lneig_second_send)
    nullify(lneig_second_recv)
    nullify(lneig_second)
    nullify(full_row_send_size)
    nullify(full_row_recv_neights)

    nullify(lcheck)
    nneig_second      = 0
    nneig_second_send = 0
    nneig_second_recv = 0
    
    debug = .false.
    call PAR_COMM_RANK_AND_SIZE(COMM % PAR_COMM_WORLD,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
    
    if( ISLAVE ) then
       !
       ! COMM % FULL_ROW_NNEIG: number of neighbors from which I need to receive halo information
       !
       call memory_alloca(par_memor,'LNEIG','par_full_row_communication_arrays',lneig,2_ip*COMM % nneig)
       jj    = 0
       kneig = 0

       do ipoin = meshe % npoi1+1,meshe % npoin_halo
          jj = jj + 1
          if( ipoin > meshe % npoin_own ) then
             dom_i = COMM % bound_owner_rank(jj)
             if( count(lneig==dom_i,KIND=ip) == 0 ) then
                kneig = kneig + 1
                if( kneig > size(lneig) ) call memory_resize(par_memor,'LNEIG','par_full_row_communication_arrays',lneig,2_ip*int(size(lneig),ip))
                lneig(kneig) = dom_i
             end if
          end if
       end do
       COMM % full_row_recv_nneig = kneig
       !
       ! LNEIG_SECOND(1:NNEIG_SECOND): List of neighbors of my neighbors
       !
       call memory_alloca(par_memor,'LNEIG_SECOND_SEND','par_full_row_communication_arrays',lneig_second_send,COMM % nneig)
       call memory_alloca(par_memor,'LNEIG_SECOND_RECV','par_full_row_communication_arrays',lneig_second_recv,COMM % nneig)
       if( COMM % nneig > 0 ) then
          lneig_second_send = COMM % neights
          nneig_second_send = COMM % nneig
       end if

       do ineig = 1,COMM % nneig
          dom_i = COMM % neights(ineig)
          call PAR_SEND_RECEIVE(nneig_second_send,nneig_second_recv,'IN MY CODE',dom_i)
          call memory_alloca(par_memor,'LNEIG_SECOND_RECV % L','par_full_row_communication_arrays',lneig_second_recv(ineig)%l,nneig_second_recv)
          call PAR_SEND_RECEIVE(lneig_second_send,lneig_second_recv(ineig)%l,'IN MY CODE',dom_i)
       end do

       if( COMM % nneig /= 0 ) then
          call memory_alloca(par_memor,'LNEIG_SECOND','par_full_row_communication_arrays',lneig_second,2_ip*COMM % nneig)
          lneig_second(1:COMM % nneig) = COMM % neights(1:COMM % nneig)
          nneig_second = COMM % nneig
          do ineig = 1,COMM % nneig
             do jneig = 1,size(lneig_second_recv(ineig)%l)
                dom_i = lneig_second_recv(ineig)%l(jneig)
                if( dom_i /= kfl_paral ) then
                   if( count(lneig_second==dom_i,KIND=ip) == 0 ) then
                      nneig_second = nneig_second + 1
                      if( nneig_second > size(lneig_second) ) call memory_resize(par_memor,'LNEIG_SECOND','par_full_row_communication_arrays',lneig_second,2_ip*int(size(lneig_second),ip))
                      lneig_second(nneig_second) = dom_i
                   end if
                end if
             end do
          end do
       end if
       !         
       ! COMM % FULL_ROW_RECV_NEIGHTS: list of neighbors
       !
       COMM % full_row_recv_dim = 0
       call memory_alloca(par_memor,trim(my_comm_name)//' % FULL_ROW_RECV_SIZE'   ,'par_full_row_communication_arrays',COMM % full_row_recv_size,   COMM % full_row_recv_nneig+1_ip)
       call memory_alloca(par_memor,trim(my_comm_name)//' % FULL_ROW_RECV_NEIGHTS','par_full_row_communication_arrays',COMM % full_row_recv_neights,COMM % full_row_recv_nneig)

       jj = 0
       kneig = 0
       do ipoin = meshe % npoi1+1,meshe % npoin_halo
          jj = jj + 1
          if( ipoin > meshe % npoin_own ) then
             dom_i = COMM % bound_owner_rank(jj)
             if( count(COMM % full_row_recv_neights==dom_i,KIND=ip) == 0 ) then
                kneig = kneig + 1
                COMM % full_row_recv_neights(kneig) = dom_i
             end if
          end if
       end do
       if( debug) print*,'A) Neighbors of ',kfl_paral,': ',COMM % full_row_recv_neights(:)
       !
       ! FULL_ROW_RECV_NEIGHTS: Reorder list of neighbors to mimic original scheduling
       !
       call memory_copy(par_memor,trim(my_comm_name)//' % FULL_ROW_RECV_NEIGHTS','par_full_row_communication_arrays',&
            COMM % full_row_recv_neights,full_row_recv_neights,'DO_NOT_DEALLOCATE',COPY_NAME='FULL_ROW_RECV_NEIGHTS')

      jneig = 0
       do ineig = 1,COMM % nneig
          dom_i = COMM % neights(ineig)
          kneig = 0
          jj    = 0
          do while( jj < COMM % full_row_recv_nneig )
             jj = jj + 1
             if( COMM % full_row_recv_neights(jj) == dom_i ) then
                kneig = jj
                jj = COMM % full_row_recv_nneig
             end if
          end do
          if( kneig > 0 ) then
             jneig = jneig + 1
             full_row_recv_neights(kneig) = 0
             COMM % full_row_recv_neights(jneig) = dom_i
          end if
       end do

       do ineig = 1,COMM % full_row_recv_nneig
          dom_i = full_row_recv_neights(ineig)
          if( dom_i > 0 ) then
             jneig = jneig + 1
             COMM % full_row_recv_neights(jneig) = dom_i
          end if
       end do
       call memory_deallo(par_memor,'FULL_ROW_RECV_NEIGHTS','par_full_row_communication_arrays',full_row_recv_neights)

       if( debug) print*,'B) Neighbors of ',kfl_paral,': ',COMM % full_row_recv_neights(:)
       !
       ! FULL_ROW_RECV_DIM, FULL_ROW_RECV_SIZE: size of communication arrays
       !
       do ineig = 1,COMM % full_row_recv_nneig
          dom_i = COMM % full_row_recv_neights(ineig)
          jj    = 0
          COMM % full_row_recv_size(ineig) = 0
          do ipoin = meshe % npoi1+1,meshe % npoin_halo
             jj = jj + 1
             if( ipoin > meshe % npoin_own ) then
                if( COMM % bound_owner_rank(jj) == dom_i ) then
                   COMM % full_row_recv_size(ineig) = COMM % full_row_recv_size(ineig) + 1
                   COMM % full_row_recv_dim = COMM % full_row_recv_dim + 1
                end if
             end if
          end do
       end do
       call graphs_number_to_linked_list(COMM % full_row_recv_nneig,COMM % full_row_recv_size)
       !
       ! COMM % FULL_ROW_RECV_PERM: Fill in receive permutation
       !
       call memory_alloca(par_memor,trim(my_comm_name)//' % FULL_ROW_RECV_PERM','par_full_row_communication_arrays',COMM % full_row_recv_perm,COMM % full_row_recv_dim)
       COMM % full_row_recv_dim = 0
       do ineig = 1,COMM % full_row_recv_nneig
          dom_i = COMM % full_row_recv_neights(ineig)
          jj    = 0
          do ipoin = meshe % npoi1+1,meshe % npoin_halo
             jj = jj + 1
             if( ipoin > meshe % npoin_own ) then
                if( COMM % bound_owner_rank(jj) == dom_i ) then
                   COMM % full_row_recv_dim = COMM % full_row_recv_dim + 1
                   COMM % full_row_recv_perm(COMM % full_row_recv_dim) = meshe % lninv_loc(ipoin)
                end if
             end if
          end do
       end do
       !do ineig = 1,COMM % full_row_recv_nneig
       !   dom_i        = COMM % full_row_recv_neights(ineig)
       !   if( kfl_paral == 1 ) print*,'recv=',dom_i,COMM % full_row_recv_perm(COMM % full_row_recv_size(ineig):COMM % full_row_recv_size(ineig+1)-1)
       !end do
       !
       ! Construct send list
       !
       if( debug ) print*,'NNEIG_SECOND= ',kfl_paral,nneig_second
       call memory_alloca(par_memor,'FULL_ROW_SEND_SIZE','par_full_row_communication_arrays',full_row_send_size,nneig_second)

       COMM % full_row_send_nneig = 0
       COMM % full_row_send_dim   = 0
       if( associated(lneig) ) lneig = 0
       if( nneig_second > 0 ) call maths_heap_sort(1_ip,nneig_second,lneig_second)

       do ineig = 1,nneig_second
          dom_i = lneig_second(ineig)
          kneig = 0
          ii    = 0
          do while( ii < COMM % full_row_recv_nneig )
             ii = ii + 1
             if( COMM % full_row_recv_neights(ii) == dom_i ) then
                kneig = ii
                ii = COMM % full_row_recv_nneig
             end if
          end do
          if( kneig == 0 ) then
             nneig_second_recv = 0
          else
             nneig_second_recv = COMM % full_row_recv_size(kneig+1)-COMM % full_row_recv_size(kneig)
          end if

          call PAR_SEND_RECEIVE(nneig_second_recv,nneig_second_send,'IN MY CODE',dom_i)
          if( nneig_second_send == 0 ) then
             !print*,'RARISIMO= ',kfl_paral,dom_i,kneig
             !call runend('RARISIMO')
          end if
          !
          ! Add dom_i to my send list if nneig_second_send > 0
          !
          if( count(lneig==dom_i,KIND=ip) == 0 .and. nneig_second_send > 0 ) then
             COMM % full_row_send_nneig = COMM % full_row_send_nneig + 1
             if( COMM % full_row_send_nneig > size(lneig) ) call memory_resize(par_memor,'LNEIG','par_full_row_communication_arrays',lneig,2_ip*int(size(lneig),ip))
             lneig(COMM % full_row_send_nneig)              = dom_i
             COMM % full_row_send_dim                       = COMM % full_row_send_dim + nneig_second_send
             full_row_send_size(COMM % full_row_send_nneig) = nneig_second_send
          end if

       end do
       !
       ! Size and list of neighbors
       !
       call memory_alloca(par_memor,trim(my_comm_name)//' % FULL_ROW_SEND_SIZE'   ,'par_full_row_communication_arrays',COMM % full_row_send_size,   COMM % full_row_send_nneig+1)
       call memory_alloca(par_memor,trim(my_comm_name)//' % FULL_ROW_SEND_NEIGHTS','par_full_row_communication_arrays',COMM % full_row_send_neights,COMM % full_row_send_nneig)
       do ineig = 1,COMM % full_row_send_nneig
          dom_i = lneig(ineig)
          COMM % full_row_send_size(ineig)    = full_row_send_size(ineig)
          COMM % full_row_send_neights(ineig) = dom_i
       end do
       call graphs_number_to_linked_list(COMM % full_row_send_nneig,COMM % full_row_send_size)
       !
       ! Exchange list of nodes to be sent
       !
       call memory_alloca(par_memor,trim(my_comm_name)//' % FULL_ROW_SEND_PERM','par_full_row_communication_arrays',COMM % full_row_send_perm,COMM % full_row_send_dim)

       ii = COMM % full_row_recv_nneig + COMM % full_row_send_nneig
       call PAR_START_NON_BLOCKING_COMM(1_ip,ii)
       call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
   
       do ineig = 1,COMM % full_row_recv_nneig
          dom_i        = COMM % full_row_recv_neights(ineig)
          ii_recv_size = COMM % full_row_recv_size(ineig+1)-COMM % full_row_recv_size(ineig)
          ii_recv      = COMM % full_row_recv_size(ineig)
          ii_send_size = 0
          if( debug ) print*,kfl_paral,' sends ',ii_recv_size,'  to  ',dom_i
          if( ii_recv > size(COMM % full_row_recv_perm) .or. ii_recv + ii_recv_size - 1 > size(COMM % full_row_recv_perm) ) then
             print*,'RECV=',kfl_paral,dom_i,ii_recv, ii_recv + ii_recv_size - 1,size(COMM % full_row_recv_perm)
             call runend('TROUBLE')
          end if
          call PAR_SEND_RECEIVE(ii_recv_size,ii_send_size,COMM % full_row_recv_perm(ii_recv:),full_row_send_perm_ok,'IN MY CODE',dom_i,'ASYNCHRONOUS')
       end do

       do ineig = 1,COMM % full_row_send_nneig
          dom_i        = COMM % full_row_send_neights(ineig)
          ii_send_size = COMM % full_row_send_size(ineig+1)-COMM % full_row_send_size(ineig)
          ii_send      = COMM % full_row_send_size(ineig)
          ii_recv_size = 0
          if( debug ) print*,kfl_paral,' recvs ',ii_send_size,' from ',dom_i
          if( ii_send > size(COMM % full_row_send_perm) .or. ii_send + ii_send_size - 1 > size(COMM % full_row_send_perm) ) then
             print*,'SEND=',kfl_paral,dom_i,ii_send, ii_send + ii_send_size - 1,size(COMM % full_row_send_perm)
             call runend('TROUBLE')
          end if
          call PAR_SEND_RECEIVE(ii_recv_size,ii_send_size,full_row_recv_perm_ok,COMM % full_row_send_perm(ii_send:),'IN MY CODE',dom_i,'ASYNCHRONOUS')
       end do

       call PAR_END_NON_BLOCKING_COMM(1_ip)
       !
       ! From global to local numbering
       !
       do ii = 1,COMM % full_row_send_dim
          ipoin_global                  = COMM % full_row_send_perm(ii)
          if( 1 == 2 ) then
             do ipoin = 1,npoin_2
                if( lninv_loc(ipoin) == ipoin_global ) then                   
                   COMM % full_row_send_perm(ii) = ipoin
                end if
             end do
          else
             ipoin                         = PAR_GLOBAL_TO_LOCAL_NODE(ipoin_global,'INCLUDING HALOS')
             !ipoin                         = htalid(htable_lninv_loc,ipoin_global)
             COMM % full_row_send_perm(ii) = ipoin
          end if
       end do
       do ii = 1,COMM % full_row_recv_dim
          ipoin_global                  = COMM % full_row_recv_perm(ii)
          if( 1 == 2 ) then
             do ipoin = 1,npoin_2
                if( lninv_loc(ipoin) == ipoin_global ) then                   
                   COMM % full_row_recv_perm(ii) = ipoin
                end if
             end do
          else
             ipoin                         = PAR_GLOBAL_TO_LOCAL_NODE(ipoin_global,'INCLUDING HALOS')
             !ipoin                         = htalid(htable_lninv_loc,ipoin_global)
             COMM % full_row_recv_perm(ii) = ipoin
          end if
      end do
       !
       ! Condense if here are empty message
       !
       jneig = 0
       do ineig = 1,COMM % full_row_send_nneig
          if( COMM % full_row_send_size(ineig+1)-COMM % full_row_send_size(ineig) > 0 ) then
             jneig = jneig + 1
             COMM % full_row_send_neights(jneig) = COMM % full_row_send_neights(ineig)
             COMM % full_row_send_size(jneig)    = COMM % full_row_send_size(ineig+1)-COMM % full_row_send_size(ineig)
          end if
       end do
       COMM % full_row_send_nneig = jneig
       call graphs_number_to_linked_list(COMM % full_row_send_nneig,COMM % full_row_send_size)

       jneig = 0
       do ineig = 1,COMM % full_row_recv_nneig
          if( COMM % full_row_recv_size(ineig+1)-COMM % full_row_recv_size(ineig) > 0 ) then
             jneig = jneig + 1
             COMM % full_row_recv_neights(jneig) = COMM % full_row_recv_neights(ineig)
             COMM % full_row_recv_size(jneig)    = COMM % full_row_recv_size(ineig+1)-COMM % full_row_recv_size(ineig)
          end if
       end do
       COMM % full_row_recv_nneig = jneig
       call graphs_number_to_linked_list(COMM % full_row_recv_nneig,COMM % full_row_recv_size)
       !
       ! Check if there is a second neighbor
       !
       num_second = 0
       num_first  = 0
       loop_ineig : do ineig = 1,COMM % full_row_recv_nneig
          dom_i = COMM % full_row_recv_neights(ineig)
          ii    = 0
          kneig = 0
          do while( ii < COMM % nneig )
             ii = ii + 1
             if( COMM % neights(ii) == dom_i ) then
                kneig = 1
                ii = COMM % nneig
             end if
          end do
          if( kneig == 0 ) then
             num_second = num_second + 1
          else
             num_first  = num_first  + 1
          end if
       end do loop_ineig

       if( num_second+num_first > 0 )then
          ratio_second_total = real(num_second,rp)/real(num_second+num_first,rp)
       else
          ratio_second_total = 0.0_rp
       end if
    end if
    !
    ! Check I don't receive results on my own nodes
    !
    ierro = 0
    if( ISLAVE ) then
       do ii = 1,COMM % full_row_recv_dim
          ipoin = COMM % full_row_recv_perm(ii)
          if( ipoin <= meshe % npoin_own ) ierro = ierro + 1
       end do
    end if
    call PAR_MAX(ierro)
    if( ierro > 0 .and. IMASTER ) then
       !call runend('SOME PARTITIONS RECEIVE ON THEIR OWN NODES')
    end if
    !
    ! Check I only send reseult of my own nodes
    !
    ierro = 0
    if( ISLAVE ) then
       do ii = 1,COMM % full_row_send_dim
          ipoin = COMM % full_row_send_perm(ii)
          if( ipoin > meshe % npoin_own ) ierro = ierro + 1
       end do

    end if
    call PAR_MAX(ierro)
    if( ierro > 0 .and. IMASTER ) then
       !call runend('SOME PARTITIONS SENDS THEIR HALO NODES')
    end if
    !
    ! Check exchange of owners is ok
    !
    ierro = 0
    if( ISLAVE ) then
       call memory_alloca(par_memor,'LCHECK'     ,'par_full_row_communication_arrays',lcheck,meshe % npoin_halo)
       do ipoin = 1,meshe % npoin_own
          lcheck(ipoin) = kfl_paral
       end do
       call PAR_INTERFACE_OWN_NODE_EXCHANGE(lcheck)
       jj = 0   
       do ipoin = meshe % npoi1+1,meshe % npoin_halo
          jj = jj + 1
          dom_i = COMM % bound_owner_rank(jj)
          if( dom_i /= lcheck(ipoin) ) then
             ierro = ierro + 1
             print*,'PROBLEM: ',kfl_paral,dom_i,lcheck(ipoin),lninv_loc(ipoin),lmast(ipoin),ipoin,meshe % npoin_own
          end if
       end do
    end if
    call PAR_MAX(ierro)
    !!if( ierro > 0 ) then
    !!   call runend('WE ARE HAVING A REALLY BAD TIME...')
    !!end if
    !call runend('O.K.!')
    !
    ! Deallocate memory
    !
    if( ISLAVE ) then
       call memory_deallo(par_memor,'LNEIG'                ,'par_full_row_communication_arrays',lneig)
       call memory_deallo(par_memor,'LNEIG_SECOND_SEND'    ,'par_full_row_communication_arrays',lneig_second_send)
       call memory_deallo(par_memor,'LNEIG_SECOND_RECV'    ,'par_full_row_communication_arrays',lneig_second_recv)
       call memory_deallo(par_memor,'LNEIG_SECOND'         ,'par_full_row_communication_arrays',lneig_second)
       call memory_deallo(par_memor,'FULL_ROW_SEND_SIZE'   ,'par_full_row_communication_arrays',full_row_send_size)
       call memory_deallo(par_memor,'FULL_ROW_RECV_NEIGHTS','par_full_row_communication_arrays',full_row_recv_neights)
       call memory_deallo(par_memor,'LCHECK'               ,'par_full_row_communication_arrays',lcheck)
    end if

    !--------------------------------------------------------------------
    !
    ! Output statistics on halos
    !
    !--------------------------------------------------------------------

    routp = 0.0_rp

    if( INOTMASTER ) then
       !
       ! Computational halos
       !
       if( meshe % npoin_own > 0 ) &
            routp(1) = 100.0_rp * abs(real(meshe % npoin_halo-meshe % npoin_own,rp)) &
            / real(meshe % npoin_own,rp)
       if( COMM % full_row_recv_nneig + COMM % full_row_send_nneig > 0 ) &
            routp(2) = 100.0_rp * 2.0_rp * abs(real(COMM % full_row_recv_nneig - COMM % full_row_send_nneig,rp)) &
            / real(COMM % full_row_recv_nneig + COMM % full_row_send_nneig,rp)
       if( COMM % full_row_recv_dim   + COMM % full_row_send_dim > 0 ) &
            routp(3) = 100.0_rp * 2.0_rp * abs(real(COMM % full_row_recv_dim   - COMM % full_row_send_dim,rp)) &
            / real(COMM % full_row_recv_dim   + COMM % full_row_send_dim,rp)
       routp(4) = 100.0_rp * ratio_second_total
       !
       ! Geometrical halos
       !
       routp(5:6) = 0.0_rp
       if( meshe % nelem > 0 ) routp(5) = 100.0_rp * real(meshe % nelem_2 - meshe % nelem,rp)/real(meshe % nelem,rp)
       if( meshe % npoin > 0 ) routp(6) = 100.0_rp * real(meshe % npoin_2 - meshe % npoin,rp)/real(meshe % npoin,rp)
    end if
    call PAR_MAX(6_ip,routp)

    call outfor(82_ip,0_ip,' ')

  end subroutine par_full_row_communication_arrays

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-06
  !> @brief   Order update after exchange
  !> @details Order the way solution is update after exchanging
  !>          arrays between interfaces
  !> 
  !-----------------------------------------------------------------------

  subroutine par_ordered_exchange_update(COMM,COMM_NAME)

    type(comm_data_par), intent(inout)        :: COMM             !< Communication array
    integer(ip)                               :: ii
    logical(lg)                               :: ifoun
    character(len=*),    intent(in), optional :: COMM_NAME
    character(20)                             :: my_comm_name

    if( present(COMM_NAME) ) then
       my_comm_name = trim(COMM_NAME)
    else
       my_comm_name = 'COMM'
    end if
    
    if( INOTMASTER .and. COMM % nneig > 0 ) then
          
       call memory_alloca(par_memor,trim(my_comm_name)//' % NEIGHTS_ORDERED','par_slexch',COMM % neights_ordered,COMM % nneig,'DO_NOT_INITIALIZE')
       call memory_alloca(par_memor,trim(my_comm_name)//' % PERM_ORDERED'   ,'par_slexch',COMM % perm_ordered   ,COMM % nneig,'DO_NOT_INITIALIZE')
       
       COMM % neights_ordered(1:COMM % nneig) = COMM % neights(1:COMM % nneig)

       do ii = 1,COMM % nneig
          COMM % perm_ordered(ii) = ii
       end do
        call maths_heap_sort(2_ip,COMM % nneig,COMM % neights_ordered,'NORMAL',COMM % perm_ordered)

       ifoun = .true.
       COMM % nneig_1 = 0
       do while( ifoun .and. COMM % nneig_1 < COMM % nneig )
          COMM % nneig_1 = COMM % nneig_1 + 1 
          if( COMM % neights_ordered(COMM % nneig_1) > int(COMM % RANK4,ip) ) then
             COMM % nneig_1 = COMM % nneig_1-1
             ifoun = .false.
          end if
       end do
       COMM % nneig_2 = COMM % nneig_1 + 1

    end if

  end subroutine par_ordered_exchange_update

   !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-15
  !> @brief   Get subdomain bounding boxes
  !> @details Get subdomain bounding boxes
  !> 
  !-----------------------------------------------------------------------

  subroutine par_bounding_box(COMM_WORLD,bobox,subox,MEMORY_COUNTER)

    MY_MPI_COMM   ,                 intent(in)    :: COMM_WORLD
    real(rp),              pointer, intent(inout) :: subox(:,:,:)
    real(rp),              pointer, intent(in)    :: bobox(:,:,:)
    real(rp),              pointer                :: bobox_loc(:,:,:)
    integer(8),   optional,         intent(inout) :: MEMORY_COUNTER(2)
    integer(ip)                                   :: my_rank,comm_size
    integer(ip)                                   :: idime,kdime
    integer(8)                                    :: memor_loc(2)
    
    if( IPARALL ) then

       bobox_loc => bobox
       memor_loc = optional_argument((/0_8,0_8/),MEMORY_COUNTER)
       kdime     = memory_size(bobox_loc,2_ip)

       call PAR_MAX(kdime)

       call PAR_COMM_RANK_AND_SIZE(COMM_WORLD,my_rank,comm_size)

       call memory_alloca(memor_loc,'SUBOX','par_bounding_box',subox,2_ip,kdime,comm_size,LBOUN1=1_ip,LBOUN2=1_ip,LBOUN3=0_ip)

       subox = 0.0_rp 
       if( my_rank /= 0 .and. associated(bobox_loc) ) then
          do idime = 1,kdime
             subox(1,idime,my_rank) = minval(bobox_loc(1,idime,:))
             subox(2,idime,my_rank) = maxval(bobox_loc(2,idime,:))
          end do
       else
          subox(1,:,my_rank) =  huge(1.0_rp)*0.1_rp
          subox(2,:,my_rank) = -huge(1.0_rp)*0.1_rp
       end if

       call PAR_SUM(subox,INCLUDE_ROOT=.true.)

       if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc

    end if

  end subroutine par_bounding_box

end module mod_par_additional_arrays
!> @}
