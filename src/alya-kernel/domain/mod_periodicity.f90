!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @addtogroup Domain
!> @{
!> @name    ToolBox for periodicity
!> @file    mod_periodicity.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox periodicity
!> @details ToolBox periodicity
!> @{
!
!-----------------------------------------------------------------------

module mod_periodicity

  use def_kintyp_basic,          only : ip,rp,lg,i1p
  use def_master,                only : lninv_loc
  use def_master,                only : npoi1,npoi2,npoi3
  use def_master,                only : npart
  use def_kintyp_comm,           only : comm_data_par
  use def_master,                only : INOTEMPTY,IEMPTY
  use def_master,                only : ISEQUEN,IPARALL
  use def_master,                only : INOTMASTER,IMASTER
  use def_master,                only : kfl_paral
  use def_master,                only : intost
  use def_domain,                only : npoin,lmast,npoin_own
  use def_domain,                only : ltypb,lnodb,nboun
  use def_domain,                only : memor_dom
  use def_domain,                only : nperi,nnode
  use def_domain,                only : kfl_codbo
  use def_domain,                only : mcodb
  use mod_graphs,                only : graphs_number_to_linked_list
  use mod_htable,                only : htable_initialization
  use mod_htable,                only : hash_t
  use mod_htable,                only : htaini
  use mod_htable,                only : htaadd
  use mod_htable,                only : htalid
  use mod_htable,                only : htalid_list
  use mod_htable,                only : htades
  use mod_htable,                only : htades_list
  use mod_memory,                only : memory_alloca
  use mod_memory,                only : memory_deallo
  use mod_memory,                only : memory_size
  use mod_memory,                only : memory_copy
  use mod_parall,                only : commd
  use mod_parall,                only : PAR_COMM_MY_CODE_ARRAY
  use mod_parall,                only : par_memor
  use mod_parall,                only : PAR_COMM_MY_CODE
  use mod_messages,              only : messages_live
  use mod_renumbering,           only : renumbering_node_arrays
  use mod_renumbering,           only : renumbering_boundary_arrays
  use mod_renumbering,           only : renumbering_update
  use mod_communications,        only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications,        only : PAR_SEND_RECEIVE
  use mod_communications,        only : PAR_ALLGATHER
  use mod_communications,        only : PAR_ALLGATHERV
  use mod_communications,        only : PAR_ALLTOALL
  use mod_communications,        only : PAR_MAX
  use mod_par_additional_arrays, only : par_global_variables_arrays
  
  implicit none

  integer(ip), pointer :: lnods_sav(:,:) 
  integer(ip), pointer :: lninv_sav(:) 
  integer(ip)          :: mnode_sav

  public :: periodicity_setup
  public :: periodicity_sequential
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-29
  !> @brief   Periodicity setup
  !> @details Add periodic nodes to the communicator
  !> 
  !-----------------------------------------------------------------------

  subroutine periodicity_setup()

    integer(ip)          :: ipoin,iperi,kpoin,ii,jj,kk
    integer(ip)          :: ineig,nneig_new,ipart,knodb
    integer(ip)          :: my_num_periodic_nodes,igo(2),kperi
    integer(ip)          :: tot_num_periodic_nodes,iperm,pnodb
    integer(ip)          :: kboun,iboun,inodb,jpoin,ifact
    integer(ip)          :: ipoin_min,ipoin_max
    integer(ip), pointer :: my_list_periodic_nodes(:)
    integer(ip), pointer :: my_list_periodic_loc(:)
    integer(ip), pointer :: list_periodic_nodes(:)
    integer(ip), pointer :: num_periodic_nodes(:)
    integer(4),  pointer :: num_periodic_nodes4(:)
    integer(ip), pointer :: npoin_send(:)
    integer(ip), pointer :: npoin_recv(:)
    integer(ip), pointer :: owner(:)
    integer(ip), pointer :: permr(:)
    logical(lg), pointer :: touch(:)
    integer(ip), pointer :: lninv_sav(:)
    integer(ip), pointer :: lpoin_send(:,:)
    type(i1p),   pointer :: list_send(:)
    type(i1p),   pointer :: list_recv(:)
    type(i1p),   pointer :: list_send_loc(:)
    type(i1p),   pointer :: list_recv_loc(:)
    ! Getting master nodes
    integer(ip)          :: my_num_master_nodes,ipoin_loc
    integer(ip)          :: tot_num_master_nodes
    integer(ip), pointer :: my_list_master_nodes(:)
    integer(ip), pointer :: list_master_nodes(:)
    integer(ip), pointer :: num_master_nodes(:)
    integer(4),  pointer :: num_master_nodes4(:)

    type(comm_data_par)  :: comm_sav
    type(hash_t)         :: htable_lninv 
    integer(ip), pointer :: local_id(:)
    real(rp)             :: time1
!    real(rp)             :: time2
    
    if( nperi > 0 ) then

       nullify(my_list_master_nodes)
       nullify(list_master_nodes)
       nullify(num_master_nodes)
       nullify(num_master_nodes4)
       
       nullify(my_list_periodic_nodes)
       nullify(my_list_periodic_loc)
       nullify(list_periodic_nodes)
       nullify(num_periodic_nodes)
       nullify(num_periodic_nodes4)

       nullify(npoin_send)
       nullify(npoin_recv)

       nullify(owner)
       nullify(permr)
       nullify(touch)

       nullify(lpoin_send)
       nullify(list_send)
       nullify(list_recv)
       nullify(list_send_loc)
       nullify(list_recv_loc)

       nullify(lninv_sav)
       nullify(local_id)
       
       if( ISEQUEN ) then

          do ipoin = 1,npoin
             jpoin = lmast(ipoin)
             if( jpoin > 0 ) then
                lmast(jpoin) = -jpoin
             end if
          end do

       else if( IPARALL ) then

          call messages_live('COMPUTING PERIODICITY')

          call cputim(time1)
          !
          ! Local htable for LNINV_LOC
          !
          if( INOTEMPTY ) then
             call htable_initialization(htable_lninv)
             call htades( htable_lninv, memor_opt=memor_dom  )
             call htaini( htable_lninv, npoin, lidson=.true., AUTOMATIC_SIZE=.true.,memor_opt=memor_dom)
             call htaadd( htable_lninv, lninv_loc, memor_opt=memor_dom)
          end if
          !
          ! Add master to the list
          !
          my_num_master_nodes = 0
          do ipoin = 1,npoin
             if( lmast(ipoin) /= 0 ) then
                my_num_master_nodes = my_num_master_nodes + 1
             end if
          end do
          call memory_alloca(memor_dom,'MY_LIST_MASTER_NODES','periodicity_setup',my_list_master_nodes,my_num_master_nodes)
          call memory_alloca(memor_dom,'NUM_MASTER_NODES'    ,'periodicity_setup',num_master_nodes    ,npart+1_ip,       'INITIALIZE',0_ip)
          call memory_alloca(memor_dom,'NUM_MASTER_NODES4'   ,'periodicity_setup',num_master_nodes4   ,int(npart+1_ip,4),'INITIALIZE',0_4)
          my_num_master_nodes = 0
          do ipoin = 1,npoin
             if( lmast(ipoin) /= 0 ) then
                my_num_master_nodes                       = my_num_master_nodes + 1
                my_list_master_nodes(my_num_master_nodes) = lmast(ipoin)
             end if
          end do
          call PAR_ALLGATHER(my_num_master_nodes,num_master_nodes,1_4)
          tot_num_master_nodes = sum(num_master_nodes)
          call memory_alloca(memor_dom,'LIST_MASTER_NODES','periodicity_setup',list_master_nodes,tot_num_master_nodes)
          do ii = 0,size(num_master_nodes)-1
             num_master_nodes4(ii) = int(num_master_nodes(ii),4)
          end do
          call PAR_ALLGATHERV(my_list_master_nodes,list_master_nodes,num_master_nodes4)
          if( INOTEMPTY ) then
             do ii = 1,memory_size(list_master_nodes)
                ipoin     = list_master_nodes(ii)
                ipoin_loc = htalid(htable_lninv,ipoin)
                if( ipoin_loc /= 0 ) then
                   lmast(ipoin_loc) = -ipoin
                end if
             end do
          end if
          call memory_deallo(memor_dom,'MY_LIST_MASTER_NODES','periodicity_setup',my_list_master_nodes)
          call memory_deallo(memor_dom,'NUM_MASTER_NODES'    ,'periodicity_setup',num_master_nodes    )
          call memory_deallo(memor_dom,'NUM_MASTER_NODES4'   ,'periodicity_setup',num_master_nodes4   )
          call memory_deallo(memor_dom,'LIST_MASTER_NODES'   ,'periodicity_setup',list_master_nodes   )
          !call cputim(time2) ; time1 = time2 - time1 ; call PAR_MAX(time1) ; if( IMASTER ) print*,'a1=',time1 ; time1 = time2
          !
          ! Deallocate htable
          !
          if( INOTEMPTY) call htades( htable_lninv ,memor_dom )
          !
          ! Copy original global numbering
          !
          call memory_copy(memor_dom,'MY_LIST_PERIODIC_NODES','periodicity_setup',lninv_loc,lninv_sav,'DO_NOT_DEALLOCATE')
          !
          ! List of periodic nodes
          !
          my_num_periodic_nodes = 0
          do ipoin = 1,npoin
             if( lmast(ipoin) /= 0 ) then
                my_num_periodic_nodes = my_num_periodic_nodes + 1
             end if
          end do
          call memory_alloca(memor_dom,'MY_LIST_PERIODIC_NODES','periodicity_setup',my_list_periodic_nodes,my_num_periodic_nodes)
          call memory_alloca(memor_dom,'MY_LIST_PERIODIC_LOC'  ,'periodicity_setup',my_list_periodic_loc,  my_num_periodic_nodes)
          call memory_alloca(memor_dom,'NUM_PERIODIC_NODES'    ,'periodicity_setup',num_periodic_nodes    ,npart+1_ip,       'INITIALIZE',0_ip)
          call memory_alloca(memor_dom,'NUM_PERIODIC_NODES4'   ,'periodicity_setup',num_periodic_nodes4   ,int(npart+1_ip,4),'INITIALIZE',0_4)

          my_num_periodic_nodes = 0
          do ipoin = 1,npoin
             if( lmast(ipoin) /= 0 ) then
                my_num_periodic_nodes                         = my_num_periodic_nodes + 1
                my_list_periodic_nodes(my_num_periodic_nodes) = abs(lmast(ipoin))
                my_list_periodic_loc(my_num_periodic_nodes)   = ipoin
             end if
          end do
          !
          ! Change lninv_loc to the numbering of the master
          !
          do ipoin = 1,npoin
             if( lmast(ipoin) /= 0 ) then
                lninv_loc(ipoin) = abs(lmast(ipoin))
             end if
          end do
          if( INOTEMPTY ) then
             ipoin_min = minval(lninv_loc(1:npoin))
             ipoin_max = maxval(lninv_loc(1:npoin))
             call htable_initialization(htable_lninv)
             call htades( htable_lninv, memor_opt=memor_dom  )
             call htaini( htable_lninv, npoin, lidson=.true., AUTOMATIC_SIZE=.true.,memor_opt=memor_dom,REPEATED_ELEMENTS=.true.)
             call htaadd( htable_lninv, lninv_loc, memor_opt=memor_dom)
          end if
          !call cputim(time2) ; time1 = time2 - time1 ; call PAR_MAX(time1) ; if( IMASTER ) print*,'a3=',time1 ; time1 = time2          
          !
          ! Gather list of masters
          !
          call PAR_ALLGATHER(my_num_periodic_nodes,num_periodic_nodes,1_4)
          tot_num_periodic_nodes = sum(num_periodic_nodes)

          call memory_alloca(memor_dom,'LIST_PERIODIC_NODES','periodicity_setup',list_periodic_nodes,tot_num_periodic_nodes)

          do ii = 0,size(num_periodic_nodes)-1
             num_periodic_nodes4(ii) = int(num_periodic_nodes(ii),4)
          end do
          call PAR_ALLGATHERV(my_list_periodic_nodes,list_periodic_nodes,num_periodic_nodes4)
          !call cputim(time2) ; time1 = time2 - time1 ; call PAR_MAX(time1) ; if( IMASTER ) print*,'b=',time1 ; time1 = time2
          !
          ! Look for coinciding nodes with the periodic ones
          ! LIST_SEND(IPART) % L = numbering to which each node is connected
          !
          call memory_alloca(memor_dom,'NPOIN_SEND'   ,'periodicity_setup',npoin_send   ,npart+1_ip,'INITIALIZE',0_ip)
          ifact = 4
10        continue
          call memory_alloca(memor_dom,'LIST_SEND'    ,'periodicity_setup',list_send    ,npart+1_ip,'INITIALIZE',0_ip)
          call memory_alloca(memor_dom,'LIST_SEND_LOC','periodicity_setup',list_send_loc,npart+1_ip,'INITIALIZE',0_ip)
          if( INOTEMPTY ) then
             call memory_alloca(memor_dom,'LPOIN_SEND','periodicity_setup',lpoin_send,2_ip,ifact*npoin)
             iperi = 0
             do ipart = 1,npart
                kperi = 0
                do jj = 1,num_periodic_nodes(ipart)
                   iperi = iperi + 1
                   kperi = kperi + 1
                   kpoin = list_periodic_nodes(iperi)                   
                   if( kpoin >= ipoin_min .and. kpoin <= ipoin_max ) then
                      local_id  => htalid_list(htable_lninv,kpoin, memor_opt=memor_dom)
                      do ii = 1,memory_size(local_id)
                         ipoin = local_id(ii)
                         if( lninv_loc(ipoin) == kpoin ) then
                            npoin_send(ipart)               = npoin_send(ipart) + 1
                            if( npoin_send(ipart) > size(lpoin_send,2) ) then
                               !
                               ! This check is very unlikely, but just in case...
                               !
                               npoin_send = 0
                               call memory_deallo(memor_dom,'LPOIN_SEND'       ,'periodicity_setup',lpoin_send)
                               call memory_deallo(memor_dom,'LIST_SEND % L'    ,'periodicity_setup',list_send)
                               call memory_deallo(memor_dom,'LIST_SEND_LOC % L','periodicity_setup',list_send_loc)
                               ifact = ifact + 1
                               goto 10
                            end if
                            lpoin_send(1,npoin_send(ipart)) = kperi ! in periodic numbering
                            lpoin_send(2,npoin_send(ipart)) = ipoin ! local node
                         end if
                      end do
                      call htades_list(local_id,memor_opt=memor_dom)
                   end if
                end do
                call memory_alloca(memor_dom,'LIST_SEND % L'    ,'periodicity_setup',list_send(ipart)     % l,npoin_send(ipart))
                call memory_alloca(memor_dom,'LIST_SEND_LOC % L','periodicity_setup',list_send_loc(ipart) % l,npoin_send(ipart))
                do ipoin = 1,npoin_send(ipart)
                   list_send(ipart)     % l(ipoin) = lpoin_send(1,ipoin)  ! in periodic numbering 
                   list_send_loc(ipart) % l(ipoin) = lpoin_send(2,ipoin)  ! local node
                end do
             end do
             call memory_deallo(memor_dom,'LPOIN_SEND','periodicity_setup',lpoin_send)
          end if
          if( INOTEMPTY ) call htades( htable_lninv ,memor_dom )
          !call cputim(time2) ; time1 = time2 - time1 ; call PAR_MAX(time1) ; if( IMASTER ) print*,'c=',time1 ; time1 = time2
          !
          ! Receive list
          !
          call memory_alloca(memor_dom,'NPOIN_RECV','periodicity_setup',npoin_recv   ,npart+1_ip,'INITIALIZE',0_ip)
          call PAR_ALLTOALL(1_ip,1_ip,npoin_send,npoin_recv)

          call memory_alloca(memor_dom,'LIST_RECV'    ,'periodicity_setup',list_recv    ,npart+1_ip,'INITIALIZE',0_ip)
          call memory_alloca(memor_dom,'LIST_RECV_LOC','periodicity_setup',list_recv_loc,npart+1_ip,'INITIALIZE',0_ip)
          do ipart = 1,npart
             call memory_alloca(memor_dom,'LIST_RECV % L'    ,'periodicity_setup',list_recv(ipart)     % l,npoin_recv(ipart))
             call memory_alloca(memor_dom,'LIST_RECV_LOC % L','periodicity_setup',list_recv_loc(ipart) % l,npoin_recv(ipart))
          end do
          nneig_new = 0
          do ipart = 1,npart
             call PAR_SEND_RECEIVE(list_send(ipart) % l,list_recv(ipart) % l,'IN MY CODE',ipart)
          end do
          do ipart = 1,npart
             do ii = 1,npoin_recv(ipart)
                iperi = list_recv(ipart) % l(ii)
                list_recv_loc(ipart) % l(ii) = my_list_periodic_loc(iperi)
             end do
          end do
          !
          ! Remove self nodes
          !
          kk = 0
          ipart = kfl_paral
          if( npoin_recv(ipart) /= npoin_send(ipart) ) call runend('MOD_PERIODICITY: SEND AND RECV SHOULD BE EQUAL') 
          do ii = 1,npoin_recv(ipart)
             iperi = list_recv(ipart) % l(ii)
             if( list_send_loc(ipart) % l(ii) /= list_recv_loc(ipart) % l(ii) ) then
                kk = kk + 1
                list_send_loc(ipart) % l(kk) = list_send_loc(ipart) % l(ii)
                list_send(ipart)     % l(kk) = list_send(ipart)     % l(ii)
                list_recv_loc(ipart) % l(kk) = list_recv_loc(ipart) % l(ii)
             end if
          end do
          npoin_send(ipart) = kk
          npoin_recv(ipart) = kk
          !
          ! Remove nodes already accounted for in communicator
          !
          call memory_alloca(memor_dom,'TOUCH','periodicity_setup',touch,npoin)
          do ipart = 1,npart
             do ii = 1,npoin_recv(ipart)
                ipoin = list_recv_loc(ipart) % l(ii)
                touch(ipoin) = .true.
             end do
             do ii = 1,npoin_send(ipart)
                ipoin = list_send_loc(ipart) % l(ii)
                touch(ipoin) = .true.
             end do
          end do

          if( INOTEMPTY ) then
             jj = 0
             do ineig = 1,commd % nneig
                kk = 0
                do ii = commd % bound_size(ineig),commd % bound_size(ineig+1)-1
                   ipoin = commd % bound_perm(ii)
                   if( .not. touch(ipoin) ) then
                      kk = kk + 1
                      jj = jj + 1
                      commd % bound_perm(jj) = ipoin
                      commd % bound_invp(jj) = ipoin
                   end if
                end do
                commd % bound_size(ineig) = kk
             end do
             call graphs_number_to_linked_list(commd % nneig,commd % bound_size)
             commd % bound_dim = commd % bound_size(commd % nneig+1)-1
          end if
          call memory_deallo(memor_dom,'TOUCH','periodicity_setup',touch)
          !
          ! New number of neighbors
          !
          nneig_new = 0
          do ipart = 1,npart
             if( ipart == kfl_paral ) then
                if( npoin_recv(ipart) > 0 .or. npoin_send(ipart) > 0 ) nneig_new = nneig_new + 1
             else
                if( npoin_recv(ipart) > 0 .or. npoin_send(ipart) > 0 ) nneig_new = nneig_new + 1         
             end if
          end do
          !call cputim(time2) ; time1 = time2 - time1 ; call PAR_MAX(time1) ; if( IMASTER ) print*,'d=',time1 ; time1 = time2

          !--------------------------------------------------------------------
          !
          ! Construct new communicator
          !
          !--------------------------------------------------------------------
          !
          ! Copy old communicator
          !
          call comm_sav % init(COMM_NAME='COMM_SAV')
          call comm_sav % copy(PAR_COMM_MY_CODE_ARRAY(1))

          PAR_COMM_MY_CODE_ARRAY(1) % nneig     = comm_sav % nneig     + nneig_new 
          PAR_COMM_MY_CODE_ARRAY(1) % bound_dim = comm_sav % bound_dim + sum(npoin_send) + 0*sum(npoin_recv)

          call PAR_COMM_MY_CODE_ARRAY(1) % deallo()
          call PAR_COMM_MY_CODE_ARRAY(1) % alloca()
          commd => PAR_COMM_MY_CODE_ARRAY(1)
          ineig =  comm_sav % nneig

          if( INOTMASTER ) then
             do ipart = 1,comm_sav % nneig
                commd % neights(ipart) = comm_sav % neights(ipart)
             end do
             do ipart = 1,comm_sav % nneig+1
                commd % bound_size(ipart) = comm_sav % bound_size(ipart)
             end do
             do ii = 1,comm_sav % bound_dim
                commd % bound_perm(ii) = comm_sav % bound_perm(ii)
                commd % bound_invp(ii) = comm_sav % bound_invp(ii)
             end do
             !
             ! Add periodic nodes to communicator
             !
             iperm = comm_sav % bound_dim

             do ipart = 1,npart

                igo = 0

                if( npoin_send(ipart) > 0 .or. npoin_recv(ipart) > 0 ) then

                   ineig = ineig + 1

                   if( kfl_paral < ipart ) then
                      if( npoin_send(ipart) > 0 ) igo(1) = 1
                      if( npoin_recv(ipart) > 0 ) igo(2) = 0
                   else if( kfl_paral == ipart ) then
                      if( npoin_send(ipart) > 0 ) igo(1) = 1
                      if( npoin_recv(ipart) > 0 ) igo(2) = 0
                   else
                      if( npoin_send(ipart) > 0 ) igo(2) = 0
                      if( npoin_recv(ipart) > 0 ) igo(1) = 2
                   end if
                   commd % bound_size(ineig+1) = commd % bound_size(ineig) + npoin_send(ipart) + 0*npoin_recv(ipart)

                   do kk = 1,2

                      select case ( igo(kk) )

                      case ( 1_ip )

                         commd % neights(ineig) = ipart

                         do ii = 1,npoin_send(ipart)
                            iperm = iperm + 1
                            commd % bound_perm(iperm) = list_send_loc(ipart) % l(ii)
                            if( kfl_paral == ipart ) then
                               commd % bound_invp(iperm) = list_recv_loc(ipart) % l(ii)
                            else
                               commd % bound_invp(iperm) = list_send_loc(ipart) % l(ii)
                            end if
                         end do

                      case ( 2_ip )

                         commd % neights(ineig) = ipart
                         do ii = 1,npoin_recv(ipart)
                            iperm = iperm + 1
                            commd % bound_perm(iperm) = list_recv_loc(ipart) % l(ii)
                            if( kfl_paral == ipart ) then
                               commd % bound_invp(iperm) = list_send_loc(ipart) % l(ii)
                            else
                               commd % bound_invp(iperm) = list_recv_loc(ipart) % l(ii)                         
                            end if
                         end do

                      case default

                         continue

                      end select

                   end do
                end if
                !end if
             end do
          end if

          commd % nneig = ineig 

          !call cputim(time2) ; time1 = time2 - time1 ; call PAR_MAX(time1) ; if( IMASTER ) print*,'e=',time1 ; time1 = time2
          if( commd % nneig > 0 ) commd % bound_dim = commd % bound_size(commd % nneig+1)-1

          call comm_sav % deallo()
          call commd    % collapse(par_memor)
          !
          ! For the moment remove scheduling
          !
          !call par_interface_exchange_scheduling(PAR_COMM_MY_CODE_ARRAY(1),ii)

          PAR_COMM_MY_CODE_ARRAY(1) % PAR_COMM_WORLD =  PAR_COMM_MY_CODE
          PAR_COMM_MY_CODE_ARRAY(1) % RANK4          =  int(kfl_paral,4)
          commd                                      => PAR_COMM_MY_CODE_ARRAY(1)

          call par_global_variables_arrays()
          !
          ! Deallocate
          !
          call memory_deallo(memor_dom,'MY_LIST_PERIODIC_NODES','periodicity_setup',my_list_periodic_nodes)
          call memory_deallo(memor_dom,'MY_LIST_PERIODIC_LOC'  ,'periodicity_setup',my_list_periodic_loc  )
          call memory_deallo(memor_dom,'LIST_PERIODIC_NODES'   ,'periodicity_setup',list_periodic_nodes)
          call memory_deallo(memor_dom,'NUM_PERIODIC_NODES'    ,'periodicity_setup',num_periodic_nodes    )
          call memory_deallo(memor_dom,'NUM_PERIODIC_NODES4'   ,'periodicity_setup',num_periodic_nodes4   )
          call memory_deallo(memor_dom,'NPOIN_SEND'            ,'periodicity_setup',npoin_send            )
          call memory_deallo(memor_dom,'NPOIN_SEND'            ,'periodicity_setup',npoin_recv            )
          call memory_deallo(memor_dom,'LPOIN_SEND'            ,'periodicity_setup',lpoin_send            )
          call memory_deallo(memor_dom,'LIST_SEND'             ,'periodicity_setup',list_send             )
          call memory_deallo(memor_dom,'LIST_RECV'             ,'periodicity_setup',list_recv             )
          call memory_deallo(memor_dom,'LIST_SEND_LOC'         ,'periodicity_setup',list_send_loc         )
          call memory_deallo(memor_dom,'LIST_RECV_LOC'         ,'periodicity_setup',list_recv_loc         )

          !--------------------------------------------------------------------
          !
          ! PERMR: nodes with owner=1 should be declared as non owner
          !
          !--------------------------------------------------------------------
          !
          ! Recover original LNINV_LOC
          !
          call memory_copy(memor_dom,'MY_LIST_PERIODIC_NODES','periodicity_setup',lninv_sav,lninv_loc)    
          !
          ! Treat multiple ownership
          !
          call memory_alloca(memor_dom,'OWNER'   ,'mod_periodicity',owner   ,npoin)
          call periodicity_owner(owner)
          !
          ! Permute
          !
          call memory_alloca(par_memor,'permR','mod_periodicity',permR,npoin)

          if( INOTEMPTY ) then
             !
             ! Interior nodes
             !
             ii = 0
             do ipoin = 1,npoin
                if( owner(ipoin) == 1 ) then
                   ii           = ii + 1
                   permR(ipoin) = ii
                end if
             end do
             npoi1 = ii
             !
             ! Own boundary nodes
             !
             do ipoin = 1,npoin
                if( owner(ipoin) == -1 ) then
                   ii           = ii + 1
                   permR(ipoin) = ii
                end if
             end do
             npoi2     = npoi1+1
             npoi3     = ii
             npoin_own = npoi3
             !
             ! Other nodes
             !
             do ipoin = 1,npoin
                if( permR(ipoin) == 0 ) then
                   ii           = ii + 1
                   permR(ipoin) = ii
                end if
             end do

          end if
          call memory_deallo(par_memor,'owner','mod_periodicity',owner)
          !
          ! Renumber node arrays
          !
          call renumbering_node_arrays(permR)
          call memory_deallo(par_memor,'permR','mod_periodicity',permR)
          call renumbering_update()

       end if
       !call cputim(time2) ; time1 = time2 - time1 ; call PAR_MAX(time1) ; if( IMASTER ) print*,'f=',time1 ; time1 = time2

       !--------------------------------------------------------------------
       !
       ! Eliminate periodic boundaries
       !
       !--------------------------------------------------------------------

       call memory_alloca(par_memor,'permR','mod_periodicity',permR,nboun)
       kboun = 0
       do iboun = 1,nboun
          knodb = 0
          pnodb = nnode(abs(ltypb(iboun)))
          do inodb = 1,pnodb
             if( abs(lmast(lnodb(inodb,iboun))) > 0 ) knodb = knodb + 1 
          end do
          if( associated(kfl_codbo) ) then
             if( knodb /= pnodb .or. kfl_codbo(iboun) /= mcodb+1 ) then
                kboun = kboun + 1
                permr(iboun) = kboun
             end if
          else
             if( knodb /= pnodb ) then
                kboun = kboun + 1
                permr(iboun) = kboun
             end if
          end if
       end do
       call renumbering_boundary_arrays(permr)
       call memory_deallo(par_memor,'permR','mod_periodicity',permR)
       !call cputim(time2) ; time1 = time2 - time1 ; call PAR_MAX(time1) ; if( IMASTER ) print*,'g=',time1 ; time1 = time2

       !--------------------------------------------------------------------
       !
       ! Mesh dimensions
       !
       !--------------------------------------------------------------------

       call par_mesh_dimensions()

       !--------------------------------------------------------------------
       !
       ! Some checks
       !
       !--------------------------------------------------------------------

       call periodicity_check()

       !--------------------------------------------------------------------
       !
       ! Treat LNINV_LOC
       ! NPERI = -1 ... Periodic nodes have their own global numbering
       ! NPERI = -2 ... Periodic nodes have the saem global numbering
       !
       !--------------------------------------------------------------------
       !call cputim(time2) ; time1 = time2 - time1 ; call PAR_MAX(time1) ; if( IMASTER ) print*,'h=',time1 ; time1 = time2

       if( 1 == 2 ) then
          nperi = -2
          do ipoin = 1,npoin
             if( lmast(ipoin) /= 0 ) then
                lninv_loc(ipoin) = abs(lmast(ipoin))
             end if
          end do
       else
          nperi = -1
       end if

    end if

  end subroutine periodicity_setup

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-12-09
  !> @brief   Check exchange with periodicity
  !> @details Check that periodicity is ok at the communication level
  !> 
  !-----------------------------------------------------------------------
  
  subroutine periodicity_check()

    integer(ip)          :: ierro,ipoin
    integer(ip), pointer :: multi(:)

    if( ISEQUEN ) return
    
    nullify(multi)
    call memory_alloca(par_memor,'MULTI','mod_periodicity',multi,npoin)
    
    do ipoin = 1,npoin
       multi(ipoin) = 0
    end do
    do ipoin = 1,npoi3
       multi(ipoin) = 1
    end do
    call PAR_INTERFACE_NODE_EXCHANGE(multi,'SUM','IN MY CODE')
    ierro = 0
    do ipoin = 1,npoin
       if( multi(ipoin) /= 1 ) then
          ierro = ierro + 1
          print*,'ERROR WITH OWN NODES 1=',kfl_paral,ipoin,lninv_loc(ipoin),multi(ipoin)
       end if
    end do
    call PAR_MAX(ierro)
    if( ierro > 0 ) call runend('O.K.!')

    do ipoin = 1,npoin
       multi(ipoin) = 0
    end do
    do ipoin = 1,npoi3
       multi(ipoin) = kfl_paral
    end do
    call PAR_INTERFACE_NODE_EXCHANGE(multi,'SUM','IN MY CODE')
    ierro = 0
    do ipoin = 1,npoi3
       if( multi(ipoin) /= kfl_paral ) then
          ierro = ierro + 1
          print*,'ERROR WITH OWN NODES 2=',kfl_paral,ipoin,lninv_loc(ipoin),multi(ipoin)
       end if
    end do
    call PAR_MAX(ierro)
    if( ierro > 0 ) call runend('O.K.!')

!!$    ierro = 0
!!$    do ipoin = 1,npoin
!!$       !multi(ipoin) = lninv_loc(ipoin)
!!$       multi(ipoin) = lmast(ipoin)
!!$    end do
!!$    call PAR_INTERFACE_NODE_EXCHANGE(multi,'DIFF','IN MY CODE')
!!$    do ipoin = 1,npoin
!!$       if( multi(ipoin) /= lninv_loc(ipoin) ) then
!!$          ierro = ierro + 1
!!$          print*,'ERROR WITH GLOBAL NUMBERING=',kfl_paral,ipoin,lninv_loc(ipoin),multi(ipoin)
!!$       end if
!!$    end do
!!$    call PAR_MAX(ierro)
!!$    if( ierro > 0 ) call runend('O.K.!')
         
    call memory_deallo(par_memor,'MULTI','mod_periodicity',multi)

  end subroutine periodicity_check
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-12-02
  !> @brief   Onwership
  !> @details Determine ownership of periodic nodes
  !>           1 is interior
  !>          -1 is own boundary
  !>           0 is other boundary
  !> 
  !-----------------------------------------------------------------------
  
  subroutine periodicity_owner(owner)

    use mod_communications
    integer(ip), pointer, intent(inout) :: owner(:)
    integer(ip), pointer                :: min_nume(:)
    integer(ip), pointer                :: min_rank(:)
    integer(ip)                         :: ipoin
    integer(ip), pointer                :: owner_lmast(:)
    
    if( IEMPTY ) return

    do ipoin = 1,npoi1
       owner(ipoin) =  1
    end do
    do ipoin = npoi2,npoi3
       owner(ipoin) = -1
    end do
    
    nullify(min_nume)
    nullify(min_rank)
    
    call memory_alloca(memor_dom,'MIN_RANK','mod_periodicity',min_rank,npoin)
    call memory_alloca(memor_dom,'MIN_NUME','mod_periodicity',min_nume,npoin)
    do ipoin = 1,npoin
       min_nume(ipoin) = huge(1_ip)
       min_rank(ipoin) = kfl_paral
    end do
    call PAR_INTERFACE_NODE_EXCHANGE(min_rank,'MIN','IN MY CODE')

    nullify(owner_lmast)
    call memory_alloca(memor_dom,'permR','mod_periodicity',owner_lmast,npoin)    
    do ipoin = 1,npoin
       if( lmast(ipoin) /= 0 ) owner_lmast(ipoin) = 1
    end do
    call PAR_INTERFACE_NODE_EXCHANGE(owner_lmast,'SUM','IN MY CODE')

    do ipoin = 1,npoin
       if( owner_lmast(ipoin) == 1 ) then
          call runend('PERIODICITY_OWNER: THIS IS AN IMPOSSIBLE SITUATION FOR NODE '//trim(intost(ipoin)))
       else if( owner_lmast(ipoin) > 1 ) then
          if( kfl_paral == min_rank(ipoin) ) then
             min_nume(ipoin) = lninv_loc(ipoin)
          else
             owner_lmast(ipoin) = 0
          end if
       end if
    end do    
       
    call PAR_INTERFACE_NODE_EXCHANGE(min_nume,'MIN','IN MY CODE')
    
    do ipoin = 1,npoin
       if( owner_lmast(ipoin) > 1 ) then
          if( lninv_loc(ipoin) == min_nume(ipoin) ) then
             owner_lmast(ipoin) = -1
          else
             owner_lmast(ipoin) =  0
          end if
       end if
    end do
    
    do ipoin = 1,npoin
       if( lmast(ipoin) /= 0 ) owner(ipoin) = owner_lmast(ipoin)
    end do
       
    call memory_deallo(memor_dom,'OWNER_LMAST','mod_periodicity',owner_lmast)

  end subroutine periodicity_owner 

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-12-14
  !> @brief   Impose periodicity
  !> @details Impose periodicity
  !> 
  !-----------------------------------------------------------------------

  subroutine periodicity_sequential(ndofn,xx)

    use def_domain
    integer(ip), intent(in)    :: ndofn
    real(rp),    intent(inout) :: xx(ndofn,*)
    integer(ip)                :: ipoin,jpoin,idime

    if( ISEQUEN ) then
       !
       ! Master = master + slaves
       !
       do ipoin = 1,npoin
          jpoin = lmast(ipoin)
          if(jpoin>0) then
             do idime = 1,ndofn
                xx(idime,jpoin) = xx(idime,jpoin) + xx(idime,ipoin) 
             end do
          end if
       end do
       !
       ! Slaves = master
       !
       do ipoin = 1,npoin
          jpoin = lmast(ipoin)
          if(jpoin>0) then
             do idime = 1,ndofn
                xx(idime,ipoin) = xx(idime,jpoin) 
             end do
          end if
       end do

    end if

  end subroutine periodicity_sequential

end module mod_periodicity
!> @}
