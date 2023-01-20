!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{                                                                   
!> @file    par_extgra.f90
!> @author  Guillaume Houzeaux
!> @date    19/11/2013
!> @brief   Extend graph to add connectivity between interface nodes              
!> @details Add the connection between interface nodes in the
!>          original graph of the matrix and reallocate memory
!>          This is necessay to avoid the following possibility,
!>          where nodes 1 and 2 are not in the graph of subdomain 2
!>
!>         o----o----o----o
!>         |    |    |    |  subdomain 1
!>         o----1----2----o
!>              |    |
!>              o----o
!>
!>         o----1    2----o
!>         |    |    |    |  subdomain 2
!>         o----o----o----o
!>               
!> @}                                                                   
!-----------------------------------------------------------------------

subroutine par_extgra()
  use def_kintyp
  use def_domain
  use def_master
  use mod_memory
  use def_kermod,         only : kfl_graph
  use mod_parall,         only : commd
  use mod_parall,         only : NODE_IN_NEIGHBOR
  use mod_parall,         only : par_memor
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_SEND_RECEIVE
  use mod_messages,       only : livinf
  implicit none
  integer(ip)          :: ipoin,ii,jj,bsize,ll,nsize,kpoin
  integer(ip)          :: kk,jpoin,dom_i,ineig,izdom,msize
  integer(ip)          :: ii_bou,nzdom_old,nzdom_send,nzdom_recv
  integer(ip), pointer :: graph_copy(:)
  integer(ip), pointer :: graph_send(:,:)
  integer(ip), pointer :: graph_recv(:,:)
  integer(ip), pointer :: graph_check_send(:)
  integer(ip), pointer :: graph_check_recv(:)
  integer(ip), pointer :: graph_number(:)
  type(i1p),   pointer :: graph(:)
  logical(lg)          :: ifoun

  if( ISEQUEN ) return

  if( kfl_graph == 1 ) then

     call livinf(0_ip,'COMPUTE EXTENDED GRAPH',0_ip)

     nullify(graph_copy)
     nullify(graph_send)
     nullify(graph_recv)
     nullify(graph_check_send)
     nullify(graph_check_recv)
     nullify(graph_number) 
     nullify(graph) 

     !----------------------------------------------------------------------
     !
     ! Copy original graph (R_DOM,C_DOM) to type GRAPH(1:NPOIN) % L(:)
     !
     !----------------------------------------------------------------------

     if( INOTMASTER ) then
        call memory_alloca(par_memor,'GRAPH',       'par_extgra',graph,npoin)
        call memory_alloca(par_memor,'GRAPH_NUMBER','par_extgra',graph_number,npoin)
        msize = 0
        do ipoin = 1,npoin
           nsize = r_dom(ipoin+1)-r_dom(ipoin)
           graph_number(ipoin) = nsize
           if( nsize > 0 ) allocate( graph(ipoin) % l(nsize) )
           do ii = 1,nsize
              graph(ipoin) % l(ii) = 0
           end do
           ii = 0
           nsize = 0
           do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
              jpoin = c_dom(izdom)
              if( ipoin > npoi1 .and. jpoin > npoi1 ) nsize = nsize + 1
              ii = ii + 1
              graph(ipoin) % l(ii) = jpoin
           end do
           msize = max(msize,nsize)
        end do
     end if

     !----------------------------------------------------------------------
     !
     ! Estimate size of communications: MSIZE=max. of connectivities
     ! between interface nodes over the subdomains
     !
     !----------------------------------------------------------------------

     call PAR_MAX(msize,'IN MY CODE')

     !----------------------------------------------------------------------
     !
     ! Exchange cnnectivities between neighboring nodes and add the ones 
     ! missing in my graph
     !
     !----------------------------------------------------------------------

     if( INOTMASTER ) then

        do ineig = 1,commd % nneig

           dom_i = commd % neights(ineig)
           bsize = commd % bound_size(ineig+1) - commd % bound_size(ineig)

           call memory_alloca(par_memor,'GRAPH_SEND','par_extgra',graph_send,msize,bsize)
           call memory_alloca(par_memor,'GRAPH_RECV','par_extgra',graph_recv,msize,bsize)

           ii_bou = 0
           do ii = commd % bound_size(ineig),commd % bound_size(ineig+1)-1        
              !
              ! For each node IPOIN, put the connected nodes JPOIN in global numbering
              !   
              ii_bou = ii_bou + 1
              ipoin  = commd % bound_perm(ii)
              kk     = 0 
              do jj = commd % bound_size(ineig),commd % bound_size(ineig+1)-1
                 jpoin = commd % bound_perm(jj)
                 izdom = r_dom(ipoin)
                 do while( izdom <= r_dom(ipoin+1)-1 )
                    if( c_dom(izdom) == jpoin ) then
                       kk                    = kk + 1
                       if( kk     > msize ) call runend('TROUBLE 1')
                       if( ii_bou > bsize ) call runend('TROUBLE 2')
                       graph_send(kk,ii_bou) = lninv_loc(jpoin)
                       izdom                 = r_dom(ipoin+1)+2
                    end if
                    izdom = izdom + 1
                 end do
              end do
           end do
           !
           ! Send/Receive to my neighbor INEIG
           !
           call PAR_SEND_RECEIVE(graph_send,graph_recv,'IN MY CODE',dom_i)

           ii_bou = 0
           do ii = commd % bound_size(ineig),commd % bound_size(ineig+1)-1        
              ipoin  = commd % bound_perm(ii)  
              ii_bou = ii_bou + 1
              !
              ! IPOIN graph: loop over what my neighbor has sent to me
              !
              kk = 1
              do while( kk <= msize )
                 jpoin = graph_recv(kk,ii_bou)

                 if( jpoin == 0 ) then
                    !
                    ! Get out of the loop
                    !
                    kk = msize  
                 else
                    !
                    ! Check if JPOIN is already in my graph
                    ! JPOIN is in global numbering
                    !
                    ll    = 1
                    ifoun = .false.
                    nsize = graph_number(ipoin)
                    do while( ll <= nsize )
                       if( lninv_loc( graph(ipoin) % l(ll) ) == jpoin ) then
                          ifoun = .true.
                          ll    = nsize
                       end if
                       ll = ll + 1
                    end do
                    !
                    ! JPOIN is not connected to IPOIN in my graph
                    !
                    if( .not. ifoun ) then
                       !
                       ! Get local numbering 
                       !
                       ll = commd % bound_size(ineig)
                       do while( ll <= commd % bound_size(ineig+1)-1 )
                          kpoin = commd % bound_perm(ll)
                          if( lninv_loc(kpoin) == jpoin ) then
                             jpoin = kpoin
                             ll    = commd % bound_size(ineig+1)+2
                          end if
                          ll = ll + 1
                       end do

                       if( ll /= commd % bound_size(ineig+1)+3 ) call runend('TROUBLE 3')
                       nsize = nsize + 1
                       call memory_resize(par_memor,'GRAPH','par_extgra',graph(ipoin) % l,nsize)
                       graph(ipoin) % l(nsize) = kpoin
                       graph_number(ipoin)     = nsize
                    end if
                 end if
                 kk = kk + 1
              end do
           end do
           call memory_deallo(par_memor,'GRAPH_SEND','par_extgra',graph_send)
           call memory_deallo(par_memor,'GRAPH_RECV','par_extgra',graph_recv)

        end do
     end if

     !----------------------------------------------------------------------
     !
     ! Reconstruct graph R_DOM,C_DOM using new connectivities
     !
     !----------------------------------------------------------------------

     if( INOTMASTER ) then

        nzdom_old = nzdom
        nzdom     = 0
        do ipoin = 1,npoin
           nzdom = nzdom + graph_number(ipoin)
        end do
        nzsol = nzdom
        
        if( nzdom > nzdom_old ) then

           call memory_deallo(memor_dom,'C_DOM','par_extgra',c_dom)
           call memory_alloca(memor_dom,'C_DOM','par_extgra',c_dom,nzdom)
           ii = 1
           jj = 0
           r_dom(1) = 1
           do ipoin = 1,npoin
              nsize = graph_number(ipoin)
              r_dom(ipoin+1) = r_dom(ipoin) + nsize
              do kk = 1,nsize
                 jj = jj + 1
                 c_dom(jj) = graph(ipoin) % l(kk)
              end do
           end do

        end if
        !
        ! Deallocate memory
        !
        call memory_deallo(par_memor,'GRAPH',       'par_extgra',graph)
        call memory_deallo(par_memor,'GRAPH_NUMBER','par_extgra',graph_number) 

     end if

     c_sol => c_dom
     r_sol => r_dom
     
     !----------------------------------------------------------------------
     !
     ! Check everything's ok
     !
     !----------------------------------------------------------------------

     return

     if( INOTMASTER ) then
        do ineig = 1,commd % nneig
           do ii = commd % bound_size(ineig),commd % bound_size(ineig+1)-1        
              dom_i = commd % neights(ineig)
              ipoin = commd % bound_perm(ii)
              nzdom_send = 0
              do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                 jpoin = c_dom(izdom)
                 if( NODE_IN_NEIGHBOR(jpoin,ineig,commd) ) then
                    nzdom_send = nzdom_send + 1
                 end if
              end do
              call PAR_SEND_RECEIVE(nzdom_send,nzdom_recv,'IN MY CODE',dom_i)
              if( nzdom_send /= nzdom_recv ) then
                 call runend('TROUBLE 4')
              else
                 allocate( graph_check_send(nzdom_send) )
                 allocate( graph_check_recv(nzdom_send) )
                 nzdom_send = 0
                 do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                    jpoin = c_dom(izdom)
                    if( NODE_IN_NEIGHBOR(jpoin,ineig,commd) ) then
                       nzdom_send = nzdom_send + 1
                       graph_check_send(nzdom_send) = lninv_loc(jpoin)
                    end if
                 end do
                 call heapsorti1(1_ip,nzdom_send,graph_check_send)

                 call PAR_SEND_RECEIVE(graph_check_send,graph_check_recv,'IN MY CODE',dom_i)
                 do jj = 1,nzdom_send
                    if( graph_check_send(jj) /= graph_check_recv(jj) ) call runend('TROUBLE 5')
                 end do
                 deallocate( graph_check_send )
                 deallocate( graph_check_recv )
              end if
           end do
        end do
     end if

  end if

end subroutine par_extgra
