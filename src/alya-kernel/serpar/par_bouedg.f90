!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_bouedg()
  !----------------------------------------------------------------------
  !****f* Parall/par_bouedg
  ! NAME
  !    par_bouedg
  ! DESCRIPTION
  !    The output of this subroutine is:
  !    - The list of boundary edges and 
  !    - For each edge ipoin-jpoin if ipoin and jpoin are own nodes 
  !      or not 
  ! USED BY
  !    Parall
  !***
  !----------------------------------------------------------------------
  use def_master
  use def_parame
  use def_parall
  use def_domain
  use mod_memchk  
  use mod_parall, only : commd
  use mod_parall, only : PAR_COMM_MY_CODE,PAR_INTEGER
  use mod_parall, only : par_memor
  use def_mpi
  implicit none

  integer(ip)          :: ipoin,ii,jj,bsize,iown1,iown2
  integer(4)           :: istat,bsize4,bsizes,bsizer,dom_i,dom_m
  integer(ip)          :: iedgg,iedgb,knode,ineig,kk,lsize
  integer(ip)          :: kowne,jpoin,node1,node2
  integer(ip)          :: ll,ipoi1,ipoi2,qpoin
  type(i1p),   pointer :: lnown(:)
  integer(ip), pointer :: nneis(:),nneir(:)
  type(i1p),   pointer :: lneis(:),lneir(:),lsend(:),lowns(:)

  if( INOTMASTER ) then
     !
     ! Candidate edges: both edge nodes are boundary nodes
     !
     !allocate( lenei(nedgg) , stat = istat )
     do iedgg = 1,nedgg
        !lenei(iedgg) = 0
        jpoin = ledgg(1,iedgg)
        ipoin = ledgg(2,iedgg)
        if( ipoin > npoi1 .and. jpoin > npoi1 ) then
           ledgg(4,iedgg) = 1
        end if
     end do
     !
     ! Count my edges with my neighbors   
     ! My own candidate edges: LNEIS
     ! My neighbor's candidate edges: LNEIR
     !
     allocate( nneis(nneig) , stat = istat )
     allocate( nneir(nneig) , stat = istat )
     allocate( lneis(nneig) , stat = istat )
     allocate( lneir(nneig) , stat = istat )
     allocate( lsend(nneig) , stat = istat )
     allocate( lowns(nneig) , stat = istat )

     do ineig = 1,nneig
        nneir(ineig) = 0
        nneis(ineig) = 0
     end do

     do iedgg = 1,nedgg
        if( ledgg(4,iedgg) > 0 ) then
           ipoin = ledgg(1,iedgg)
           jpoin = ledgg(2,iedgg)
           do ineig = 1,nneig
              knode = 0
              jj    = commd % bound_size(ineig)
              dom_i = commd % neights(ineig)
              do while( jj <= commd%bound_size(ineig+1)-1 .and. knode < 2 )
                 if( commd % bound_perm(jj) == ipoin .or. commd % bound_perm(jj) == jpoin ) then
                    knode = knode + 1
                 end if
                 jj = jj + 1
              end do
              if( knode == 2 ) then
                 nneis(ineig) = nneis(ineig) + 1
              end if
           end do
        end if
     end do
     !
     ! List of candidate common edges with my neighbors
     !    
     do ineig = 1,nneig
        allocate( lneis(ineig) % l(max(1_ip,2*nneis(ineig))) , stat = istat )
        allocate( lsend(ineig) % l(max(1_ip,nneis(ineig)))   , stat = istat )
        allocate( lowns(ineig) % l(max(1_ip,nneis(ineig)))   , stat = istat )
        nneis(ineig) = 0
        nneir(ineig) = 0
     end do

     do iedgg = 1,nedgg
        if( ledgg(4,iedgg) > 0 ) then
           ipoin = ledgg(1,iedgg)
           jpoin = ledgg(2,iedgg)
           node1 = lninv_loc(ipoin)
           node2 = lninv_loc(jpoin)
           do ineig = 1,nneig
              knode = 0
              jj    = commd % bound_size(ineig)
              dom_i = commd % neights(ineig)
              do while( jj <= commd%bound_size(ineig+1)-1 .and. knode < 2 )
                 if( commd % bound_perm(jj) == ipoin .or. commd % bound_perm(jj) == jpoin ) then
                    knode = knode + 1
                 end if
                 jj = jj + 1
              end do
              if( knode == 2 ) then
                 nneis(ineig) = nneis(ineig) + 1
                 lsend(ineig) % l(nneis(ineig)) = iedgg
                 if( node1 > node2 ) then
                    lneis(ineig) % l( (nneis(ineig)-1)*2+1 ) = node2
                    lneis(ineig) % l( (nneis(ineig)-1)*2+2 ) = node1         
                 else
                    lneis(ineig) % l( (nneis(ineig)-1)*2+1 ) = node1
                    lneis(ineig) % l( (nneis(ineig)-1)*2+2 ) = node2
                 end if
              end if
           end do
        end if
     end do

#ifdef MPI_OFF
#else
     !
     ! Send/receive # my common edges to/from my neighbors
     !
     do ineig = 1, nneig

        dom_i  = commd%neights(ineig)
        bsize4 = int(1_ip,4)

        call MPI_Sendrecv( &
             nneis(ineig:), bsize4,        &
             PAR_INTEGER,  dom_i, 0_4,     &
             nneir(ineig:), bsize4,        &
             PAR_INTEGER,  dom_i, 0_4,     &
             PAR_COMM_MY_CODE, status, istat )
     end do
#endif

     do ineig = 1,nneig
        lsize = max(1_ip,nneir(ineig))
        allocate( lneir(ineig) % l(2*lsize) , stat = istat )
     end do

#ifdef MPI_OFF
#else
     !
     ! Send/receive list of candidate common edges to/from my neighbors
     !
     do ineig = 1, nneig

        dom_i  = commd%neights(ineig)
        bsizes = int(nneis(ineig)*2,4)
        bsizer = int(nneir(ineig)*2,4) 

        call MPI_Sendrecv( &
             lneis(ineig) % l(1:), bsizes, &
             PAR_INTEGER,  dom_i, 0_4,     &
             lneir(ineig) % l(1:), bsizer, &
             PAR_INTEGER,  dom_i, 0_4,     &
             PAR_COMM_MY_CODE, status, istat )

     end do
#endif
     !
     ! List of neighbor's owning
     !
     allocate( lnown(nneig) , stat = istat )
     do ineig = 1, nneig
        bsize  = commd%bound_size(ineig+1) - commd%bound_size(ineig)
        allocate( lnown(ineig) % l(bsize) , stat = istat )
        kk = 0
        do jj = commd%bound_size(ineig),commd%bound_size(ineig+1)-1
           kk = kk + 1
           lnown(ineig) % l(kk) = lownr_par(jj)
        end do
     end do
     !
     ! Eliminate edges not shared by any other subdomain: can be an interior edge
     ! NEDGB: Number of boundary edges
     !
     nedgb = 0
     do ineig = 1,nneig
        dom_i = commd%neights(ineig)
        do ii = 1,nneis(ineig)
           iedgg = lsend(ineig) % l(ii)
           if( ledgg(4,iedgg) > 0 ) then
              ipoin = lneis(ineig) % l( (ii-1)*2+1 )
              jpoin = lneis(ineig) % l( (ii-1)*2+2 )
              knode = 0
              do kk = 1,nneir(ineig)
                 node1 = lneir(ineig) % l( (kk-1)*2+1 )
                 node2 = lneir(ineig) % l( (kk-1)*2+2 )
                 if( node1 == ipoin .and. node2 == jpoin ) knode = 2
              end do
              if( knode == 2 ) then
                 !
                 ! Neighbor INEIG has also edge IEDGG
                 !
                 nedgb = nedgb + 1
                 ledgg(4,iedgg) = -nedgb                 
              end if
           end if
        end do
     end do    
     do iedgg = 1,nedgg
        if( ledgg(4,iedgg) > 0 ) then
           ledgg(4,iedgg) =  0              ! Candidate edge is not good
        else
           ledgg(4,iedgg) = -ledgg(4,iedgg) ! Candidate edge is shared by another subdomain
        end if
     end do
     !
     ! Allocate memory for boundary edges
     !
     allocate( ledgc(4+nneig,nedgb), stat=istat )
     call memchk(zero,istat,par_memor,'LEDGC','par_submsh',ledgc) 
     !
     ! 1. LEDGC: Fill in edge table. List of neighbors sharing this edge
     ! 2. LOWNS: Number of nodes owned by my neighbors for my boundary edges
     !
     do ineig = 1,nneig
        dom_i = commd%neights(ineig)
        do ii = 1,nneis(ineig)
           iedgg = lsend(ineig) % l(ii)
           ipoin = lneis(ineig) % l( (ii-1)*2+1 )
           jpoin = lneis(ineig) % l( (ii-1)*2+2 )
           knode = 0
           kk    = 0
           do while( kk < nneir(ineig) )
              kk    = kk + 1
              node1 = lneir(ineig) % l( (kk-1)*2+1 )
              node2 = lneir(ineig) % l( (kk-1)*2+2 )
              if( node1 == ipoin .and. node2 == jpoin ) then
                 knode = 2
                 kk    = nneir(ineig)
              end if
           end do
           if( knode == 2 ) then
              !
              ! INEIG has also my edge II
              !
              iedgb = ledgg(4,iedgg)
              kk    = ledgc(4,iedgb) + 1
              ledgc(4,   iedgb) = kk                     ! A2: # subd. that share this edge
              ledgc(4+kk,iedgb) = ineig                  ! A3: List of subd. that share this edge
              !
              ! LL = Number of nodes owned by neighbor INEIG for edge II
              !
              kk = 0
              ll = 0
              do jj = commd%bound_size(ineig),commd%bound_size(ineig+1)-1
                 kk = kk + 1
                 qpoin = lninv_loc(commd%bound_perm(jj))
                 if( lnown(ineig) % l(kk) == 1 ) then
                    if( qpoin == ipoin .or. qpoin == jpoin ) then
                      ll = ll + 1
                    end if
                 end if                 
              end do 
              lowns(ineig) % l(ii) = ll ! My neighbor INEIG has LL own nodes for this edge (II)
           else
              lowns(ineig) % l(ii) = 0
              lsend(ineig) % l(ii) = 0
           end if
        end do
     end do
 
     do iedgg = 1,nedgg
 
        if( ledgg(4,iedgg) > 0 ) then

           ipoi1          = ledgg(1,iedgg)
           ipoi2          = ledgg(2,iedgg)
           iedgb          = ledgg(4,iedgg)
           ledgc(1,iedgb) = iedgg
           ledgc(2,iedgb) = lninv_loc(ipoi1)          ! A1: global edge node
           ledgc(3,iedgb) = lninv_loc(ipoi2)          ! A1: global edge node
           iown1          = 0
           iown2          = 0
           if( ipoi1 >= npoi2 .and. ipoi1 <= npoi3 ) iown1 = 1
           if( ipoi2 >= npoi2 .and. ipoi2 <= npoi3 ) iown2 = 1

           if( iown1 == 1 .and. iown2 == 1 ) then
              !
              ! Both are mine: node is mine
              !
              kowne = -1_ip

           else 
              ! 
              ! 1--0 / 0--1 / 0--0
              !
              ll    = 0
              dom_m = kfl_paral
              do kk = 1,ledgc(4,iedgb)
                 ineig = ledgc(4+kk,iedgb)
                 dom_i = commd%neights(ineig)
                 ii    = 0
                 do while( ii < nneis(ineig) )
                    ii = ii + 1
                    if( iedgg == lsend(ineig) % l(ii) ) then
                       dom_m = min(dom_i,dom_m)
                       ll = ll + lowns(ineig) % l(ii)
                    end if
                 end do
              end do

              if( iown1 == 1 .or. iown2 == 1 ) then
                 !
                 ! Only one is mine: IPOIN
                 ! I have:     0---------1
                 ! INEIG has:  0---------0  (stored in lnown(ineig) % l(2,kk))
                 !         or  1---------0
                 !         
                 !
                 if( iown1 == 1 ) then
                    ipoin = ipoi1
                    jpoin = ipoi2
                 else if( iown2 == 1 ) then
                    ipoin = ipoi2
                    jpoin = ipoi1
                 end if
                 if( ll == 0 ) then
                    kowne = -1_ip
                 else 
                    if( lninv_loc(ipoin) <= ledgc(2,iedgb) .and. lninv_loc(ipoin) <= ledgc(3,iedgb) ) then
                       kowne = -1_ip
                    else
                       kowne =  1_ip
                    end if
                 end if
              else
                 !
                 ! None is mine
                 !
                 if( ll == 0 ) then
                    if( kfl_paral == dom_m ) then
                       kowne = -1_ip
                    else
                       kowne =  1_ip
                    end if
                 else
                    kowne = 1_ip
                 end if
              end if
           end if

           ledgg(3,iedgg) = abs(ledgg(3,iedgg)) * kowne

        end if
     end do
     !
     ! Deallocate memory
     !
     do ineig = 1, nneig
        deallocate( lnown(ineig) % l , stat = istat )
        deallocate( lsend(ineig) % l , stat = istat )
        deallocate( lneir(ineig) % l , stat = istat )
        deallocate( lneis(ineig) % l , stat = istat )
     end do
 
     deallocate( lnown , stat = istat )
     deallocate( lsend , stat = istat )
     deallocate( lneir , stat = istat )
     deallocate( lneis , stat = istat )
     deallocate( nneis , stat = istat )
     deallocate( nneir , stat = istat )

  end if

end subroutine par_bouedg
