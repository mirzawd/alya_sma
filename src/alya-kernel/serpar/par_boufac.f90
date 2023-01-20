!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_boufac()
  !-----------------------------------------------------------------------
  !****f* Parall/par_boufac
  ! NAME
  !    par_boufac
  ! DESCRIPTION
  !    This subroutine determines neighbors that share the same edge 
  !    IEDGG (IEDGB)
  ! USED BY
  !    Parall
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_parame
  use def_parall
  use def_domain
  use mod_memchk
  use mod_parall,     only : commd
  use mod_parall,     only : PAR_COMM_MY_CODE,PAR_INTEGER
  use mod_parall,     only : par_memor
  use def_elmgeo,     only : element_type
  use mod_maths_sort, only : maths_heap_sort
  use def_mpi
  implicit none

  integer(4)           :: bsizer,bsizes,istat,bsize4,dom_i
  integer(ip)          :: ii,jj,bsize,iowns,iothe,idof
  integer(ip)          :: ifacg,ineig,kk
  integer(ip)          :: t,i,j,kowne
  integer(ip)          :: lsize,knode,ll,mnodm,pflty,pnodf
  integer(ip)          :: lnode(16),nodes(16),inode
  integer(ip), pointer :: nneis(:),nneir(:),giscl(:),lfnei(:)
  type(i1p),   pointer :: lneis(:),lneir(:),lfass(:)

  if( INOTMASTER ) then
     !
     ! Candidate boundary faces
     !
     allocate( lfnei(nfacg) , stat = istat )
     mnodm = max(4_ip,mnodb)

     do ifacg = 1,nfacg
        lfnei(ifacg) = 0
        pflty = lfacg(mnodm+2,ifacg)
        pnodf = element_type(pflty) % number_nodes
        if( all(lfacg(1:pnodf,ifacg)>npoi1) ) then
           lfacg(mnodm+5,ifacg) = 1
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
     allocate( lfass(nneig) , stat = istat )

     do ineig = 1,nneig
        nneir(ineig) = 0
        nneis(ineig) = 0
     end do

     do ifacg = 1,nfacg
        if( lfacg(4,ifacg) > 0 ) then
           pflty          = lfacg(mnodm+2,ifacg)
           pnodf          = element_type(pflty) % number_nodes
           lnode(1:pnodf) = lfacg(1:pnodf,ifacg)  
           do ineig = 1,nneig
              knode = 0
              jj    = commd % bound_size(ineig)
              dom_i = commd % neights(ineig)
              do while( jj <= commd%bound_size(ineig+1)-1 .and. knode < pnodf )
                 if( any(lnode(1:pnodf)==commd % bound_perm(jj)) ) then
                    knode = knode + 1
                 end if
                 jj = jj + 1
              end do
              if( knode == pnodf ) then
                 nneis(ineig) = nneis(ineig) + 1
              end if
           end do
        end if
     end do
     !
     ! List of candidate common faces with my neighbors
     !    
     do ineig = 1,nneig
        allocate( lneis(ineig) % l(max(1_ip,mnodm*nneis(ineig))) , stat = istat )
        allocate( lfass(ineig) % l(max(1_ip,      nneis(ineig))) , stat = istat )
        nneis(ineig) = 0
        nneir(ineig) = 0
     end do

     do ifacg = 1,nfacg
        if( lfacg(4,ifacg) > 0 ) then
           pflty          = lfacg(mnodm+2,ifacg)
           pnodf          = element_type(pflty) % number_nodes
           lnode(1:pnodf) = lfacg(1:pnodf,ifacg)
           nodes(1:pnodf) = lninv_loc(lnode(1:pnodf))
           do ineig = 1,nneig
              knode = 0
              jj    = commd % bound_size(ineig)
              dom_i = commd % neights(ineig)
              do while( jj <= commd%bound_size(ineig+1)-1 .and. knode < pnodf )
                 if( any(lnode(1:pnodf)==commd % bound_perm(jj)) ) then
                    knode = knode + 1
                 end if
                 jj = jj + 1
              end do
              if( knode == pnodf ) then
                 nneis(ineig) = nneis(ineig) + 1
                 lfass(ineig) % l(nneis(ineig)) = ifacg
                 do ii = 1,pnodf
                    lneis(ineig) % l( (nneis(ineig)-1)*mnodm+ii) = nodes(ii)
                 end do  
                 ii = (nneis(ineig)-1)*mnodm+1
                 call sortin(pnodf,lneis(ineig) % l(ii))

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
        allocate( lneir(ineig) % l(mnodm*lsize) , stat = istat )
     end do

#ifdef MPI_OFF
#else
     !
     ! Send/receive list of candidate common faces to/from my neighbors
     !
     do ineig = 1, nneig

        dom_i  = commd%neights(ineig)
        bsizes = int(max(1_ip,nneis(ineig)*mnodm),4)
        bsizer = int(max(1_ip,nneir(ineig)*mnodm),4) 

        call MPI_Sendrecv( &
             lneis(ineig) % l(1:), bsizes, &
             PAR_INTEGER,  dom_i, 0_4,     &
             lneir(ineig) % l(1:), bsizer, &
             PAR_INTEGER,  dom_i, 0_4,     &
             PAR_COMM_MY_CODE, status, istat )

     end do
#endif

     !
     ! Eliminate edges not shared by any other subdomain: can be an interior edge
     ! NFACB: Number of boundary faces
     !
     nfacb = 0
     do ineig = 1,nneig
        dom_i = commd%neights(ineig)
        do ii = 1,nneis(ineig)
           ifacg = lfass(ineig) % l(ii)
           if( lfacg(4,ifacg) > 0 ) then
              pflty = lfacg(mnodm+2,ifacg)
              pnodf = element_type(pflty) % number_nodes
              do inode = 1,pnodf
                 lnode(inode) = lneis(ineig) % l( (ii-1)*mnodm+inode )
              end do 
              knode = 0
              kk    = 0
              do while( kk < nneir(ineig) )
                 kk    = kk + 1
                 do inode = 1,pnodf
                    nodes(inode) = lneir(ineig) % l( (kk-1)*mnodm+inode )
                 end do
                 if( all(lnode(1:pnodf)==nodes(1:pnodf)) ) then
                    knode = pnodf
                    kk    = nneir(ineig)
                 end if
              end do
              if( knode == pnodf ) then
                 !
                 ! Neighbor INEIG has also edge IFACG
                 !
                 nfacb                =  nfacb + 1
                 lfacg(mnodm+5,ifacg) = -nfacb
                 lfnei(ifacg)         =  ineig
              end if
           end if
        end do
     end do
     do ifacg = 1,nfacg
        if( lfacg(mnodm+5,ifacg) > 0 ) then
           lfacg(mnodm+5,ifacg) =  0                    ! Candidate face is not good
        else
           lfacg(mnodm+5,ifacg) = -lfacg(mnodm+5,ifacg) ! Candidate face is shared by another subdomain
        end if
     end do
     !
     ! Allocate memory for boundary edges
     !
     allocate( lfacb(mnodm+2+1,nfacb), stat=istat )
     call memchk(zero,istat,par_memor,'LFACB','par_submsh',lfacb) 
     !
     ! Order boundary nodes
     !
     bsize = 0
     do ineig = 1,nneig
        bsize = bsize + commd%bound_size(ineig+1)-commd%bound_size(ineig)
     end do
     call memgen(1_ip,bsize,0_ip)
     allocate(giscl(bsize))
     ii = 0
     kk = 1
     do ineig = 1, nneig
        do jj = commd%bound_size(ineig),commd%bound_size(ineig+1)-1
           ii = ii + 1
           gisca(ii) = commd % bound_perm(ii)
           giscl(ii) = ii
        end do
        bsize = commd%bound_size(ineig+1)-commd%bound_size(ineig)
        call maths_heap_sort(2_ip,bsize,gisca(kk:),'NOTHING',giscl(kk:)) 
        !call heapsorti2(2_ip,bsize,gisca(kk),giscl(kk))
        kk = kk + bsize 
     end do
     !
     ! Fill in boundary face
     !
     nfacb = 0
     do ifacg = 1,nfacg
        if( lfacg(mnodm+5,ifacg) > 0 ) then
           nfacb                  = nfacb + 1
           pflty                  = lfacg(mnodm+2,ifacg)
           pnodf                  = element_type(pflty) % number_nodes
           lfacb(1,nfacb)         = ifacg
           lfacb(2:pnodf+1,nfacb) = lfacg(1:pnodf,ifacg) ! Local face node
           
           do i = 1,pnodf-1
              do j = i+1,pnodf
                 if( lfacb(i+1,nfacb) > lfacb(j+1,nfacb) ) then
                    t                = lfacb(i+1,nfacb)
                    lfacb(i+1,nfacb) = lfacb(j+1,nfacb)
                    lfacb(j+1,nfacb) = t
                 end if
              end do
           end do
           lnode(1:pnodf) = lfacb(2:pnodf+1,nfacb)
           ineig          = lfnei(ifacg)
           dom_i          = commd % neights(ineig)
           
           iowns          = 0
           iothe          = 0
           
           ii             = commd%bound_size(ineig) - 1  
           kk             = 0

           do while( ii < commd%bound_size(ineig+1)-1 )
              ii = ii + 1
              loop_idof: do idof = 1,pnodf
                 if( kk == idof-1 .and. lnode(idof) == gisca(ii) ) then
                    kk    = kk + 1
                    ll    = giscl(ii)
                    iothe = iothe + lownr_par(ll)
                    iowns = iowns + lowns_par(ll)
                    exit loop_idof
                 end if
              end do loop_idof
           end do
           !
           ! Order boundary face nodes
           !
           lfacb(2:pnodf+1,nfacb) = lninv_loc(lnode(1:pnodf))
           do i = 1,pnodf-1
              do j = i+1,pnodf
                 if( lfacb(i+1,nfacb) > lfacb(j+1,nfacb) ) then
                    t                = lfacb(i+1,nfacb)
                    lfacb(i+1,nfacb) = lfacb(j+1,nfacb)
                    lfacb(j+1,nfacb) = t
                 end if
              end do
           end do
           lfacb(mnodm+3,nfacb) = ineig
           !
           ! Decide who owns the node
           !
           if( iowns > iothe ) then
              kowne = -1
           else if( iothe > iowns ) then
              kowne =  1
           else
              if( dom_i < kfl_paral ) then
                 kowne = -1
              else
                 kowne =  1
              end if
           end if
           if( iowns + iothe > pnodf ) then
              print*,'PAR_BOUFAC: OUPS=',ineig,iowns,iothe
              call runend('PAR_BOUFAC: OUPS')               
           end if
           lfacg(mnodm+1,ifacg) = abs(lfacg(mnodm+1,ifacg)) * kowne
        end if
     end do
     !
     ! Deallocate memory
     !
     deallocate(giscl)
     call memgen(3_ip,bsize,0_ip)

     do ineig = 1, nneig
        deallocate( lfass(ineig) % l , stat = istat )
        deallocate( lneir(ineig) % l , stat = istat )
        deallocate( lneis(ineig) % l , stat = istat )
     end do

     deallocate( lfass , stat = istat )
     deallocate( lneir , stat = istat )
     deallocate( lneis , stat = istat )
     deallocate( nneis , stat = istat )
     deallocate( nneir , stat = istat )
     deallocate( lfnei , stat = istat )

  end if


end subroutine par_boufac
