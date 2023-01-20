!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_submsh()
  !-----------------------------------------------------------------------
  !****f* domain/par_submsh
  ! NAME
  !    domain
  ! DESCRIPTION
  !    Create edge table
  ! OUTPUT
  !    NNEDG ... Number of edges
  !    LEDGG ... Edge table
  !    LEDGC ... Boundary edge table (when Parall is on)
  ! USED BY
  !    Turnon
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_domain
  use def_master
  use def_parall
  use mod_memory
  use mod_parall,         only : commd,PAR_COMM_MY_CODE
  use mod_parall,         only : PAR_INTEGER
  use mod_parall,         only : par_memor
  use mod_renumbering,    only : renumbering_node_arrays
  use mod_communications, only : PAR_ALLGATHER
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_BROADCAST
  use mod_alya2metis,     only : alya2metis_METIS_NodeND
  use def_elmgeo,         only : element_type
  use mod_elmgeo,         only : elmgeo_mm_face_nodes
  use mod_maths_sort,     only : maths_geometrical_sort_using_coordinates
  use def_mpi
  implicit none 

  integer(ip)              :: iedgg,iedgb,ii,ineig,kk,ifacg,kboun,newfa
  integer(ip)              :: bsize,kpoin,kface,ini,iboun,pflty,pnodf
  integer(ip)              :: kedge,ipoin,jpoin,ifacb,ielem,kelem
  integer(ip), pointer     :: ledg1(:)      => null()
  integer(ip), pointer     :: ledg2(:)      => null()
  integer(ip), pointer     :: ledgc_loc(:)  => null()
  integer(ip), pointer     :: lfac1(:)      => null()
  integer(ip), pointer     :: lfac2(:)      => null()
  integer(ip), pointer     :: lfac3(:)      => null()
  integer(ip), pointer     :: lfac4(:)      => null()
  integer(ip), pointer     :: lfaci(:,:)
  integer(ip), pointer     :: lfacb_loc(:)  => null()
  integer(ip), pointer     :: ltypb_loc(:)  => null()
  integer(ip), pointer     :: permR(:)      => null()
  type(i1p),   pointer     :: bperm(:)      => null()
  integer(4)               :: istat,bsize4,dom_i
  real(rp)                 :: time0,time1,time2,time3,time4
  real(rp)                 :: time5,time6,xcoor(3,16)
  ! Renumbering
  integer(ip)              :: isend(2),mnodm,lnode(16),bsizf
  integer(ip)              :: npoi1_old,npoi2_old,npoi3_old
  integer(ip), pointer     :: perm2(:) 
  integer(ip), pointer     :: invp2(:) 
  integer(ip), pointer     :: ia(:)    
  integer(ip), pointer     :: ja(:)    

  if( IPARALL ) then

     nullify(ledg2)
     nullify(ledgc)
     nullify(lfac1)
     nullify(lfac2)
     nullify(lfac3)
     nullify(lfac4)
     nullify(lfaci)
     nullify(lfacb_loc)
     nullify(ltypb_loc)
     nullify(permR)
     nullify(bperm)
     nullify(perm2)
     nullify(invp2)
     nullify(ia)
     nullify(ja)

     !-------------------------------------------------------------------
     !
     ! List of my own nodes:            LOWNS_PAR
     ! List of my neighbor's own nodes: LOWNR_PAR
     !
     !-------------------------------------------------------------------

     npoi1_old = npoi1
     npoi2_old = npoi2
     npoi3_old = npoi3
     mnodm     = max(4_ip,mnodb)

     call cputim(time0)
     if( INOTMASTER ) then

        call memory_alloca(par_memor,'LOWNS_PAR','par_submsh',lowns_par,commd%bound_dim)
        call memory_alloca(par_memor,'LOWNR_PAR','par_submsh',lownr_par,commd%bound_dim)
        do ii = 1,commd % bound_dim
           ipoin = commd % bound_perm(ii)
           if( ipoin >= npoi2 .and. ipoin <= npoi3 ) then
              lowns_par(ii) = 1
           else
              lowns_par(ii) = 0
           end if
        end do

#ifdef MPI_OFF
#else
        do ii = 1, commd % nneig
           dom_i  = commd % neights(ii)
           ini    = commd % bound_size(ii)
           bsize  = commd % bound_size(ii+1) - ini 
           bsize4 = int(bsize,4)
           call MPI_Sendrecv(                 & 
                lowns_par(ini:), bsize4,      &
                PAR_INTEGER,  dom_i, 0_4,     &
                lownr_par(ini:), bsize4,      &
                PAR_INTEGER,  dom_i, 0_4,     &
                PAR_COMM_MY_CODE, status, istat )
        end do
#endif

     end if

     call cputim(time1)

     !-------------------------------------------------------------------
     !
     ! A. Determine neighbors that share the same edge IEDGG (IEDGB)
     !
     !-------------------------------------------------------------------

     call par_bouedg()
     call cputim(time2)

     call par_boufac()
     call cputim(time3)

     if( INOTMASTER ) then

        call memory_deallo(par_memor,'LOWNS_PAR','par_submsh',lowns_par)
        call memory_deallo(par_memor,'LOWNR_PAR','par_submsh',lownr_par)

     end if

     if( INOTMASTER ) then

        !-------------------------------------------------------------------
        !
        ! B. Reorder. faces were already ordered in par_boufac
        !
        !-------------------------------------------------------------------

        if( nedgb > 0 ) then
           call memory_alloca(par_memor,'LEDGC_LOC','par_submsh',ledgc_loc,nedgb)
           call memory_alloca(par_memor,'LEDG1',    'par_submsh',ledg1,    nedgb)
           call memory_alloca(par_memor,'LEDG2',    'par_submsh',ledg2,    nedgb)
        end if

        if( nfacb > 0 ) then
           call memory_alloca(par_memor,'LFACB_LOC','par_submsh',lfacb_loc,nfacb)
           call memory_alloca(par_memor,'LFACB_LOC','par_submsh',ltypb_loc,nfacb)
           call memory_alloca(par_memor,'LFACI',    'par_submsh',lfaci,mnodm,nfacb)
        end if

        do iedgb = 1,nedgb
           iedgg = ledgc(1,iedgb)
           if( ledgc(2,iedgb) > ledgc(3,iedgb) ) then
              ipoin          = ledgc(2,iedgb)            ! B1: reorder boundary edge
              ledgc(2,iedgb) = ledgc(3,iedgb)            ! B1: reorder boundary edge
              ledgc(3,iedgb) = ipoin
              ipoin          = ledgg(1,iedgg)            ! B2: reorder corresponding global edge
              ledgg(1,iedgg) = ledgg(2,iedgg)            ! B2: reorder corresponding global edge
              ledgg(2,iedgg) = ipoin
           end if
        end do


        call memory_alloca(par_memor,'BPERM','par_submsh',bperm,commd % nneig)

        call memgen(1_ip,npoin,0_ip)

        do ineig = 1,commd % nneig

           !----------------------------------------------------------------
           !
           ! C. Fill in new edge node array in common with neighbor INEIG
           !
           !----------------------------------------------------------------

           dom_i = commd%neights(ineig)
           kedge = 0
           do iedgb = 1,nedgb
              do ii = 1,ledgc(4,iedgb)
                 if( ledgc(4+ii,iedgb) ==  ineig ) then
                    iedgg            = ledgc(1,iedgb)
                    kedge            = kedge + 1
                    ledgc_loc(kedge) = abs(ledgg(3,iedgg))
                    ledg1(kedge)     = ledgc(2,iedgb)
                    ledg2(kedge)     = ledgc(3,iedgb)
                 end if
              end do
           end do

           kface = 0
           do ifacb = 1,nfacb
              if( lfacb(mnodm+3,ifacb) == ineig ) then
                 ifacg                = lfacb(1,ifacb)
                 kface                = kface + 1
                 pflty                = lfacg(mnodm+2,ifacg)
                 lfacb_loc(kface)     = abs(lfacg(mnodm+1,ifacg))
                 ltypb_loc(kface)     = pflty
                 pnodf                = element_type(pflty) % number_nodes
                 lfaci(1:pnodf,kface) = lfacb(2:pnodf+1,ifacb)
              end if
           end do

           !----------------------------------------------------------------
           !
           ! D. Order edges with INEIG so that the send and receive correspond
           !
           !    This only works when all faces are of. We should also order
           !    ltypb_loc!
           !
           !----------------------------------------------------------------

           if( nedgb > 0 ) call hsort2    (2_ip,kedge,ledg1,ledg2,ledgc_loc)
           if( nfacb > 0 ) call par_sortii(kface,mnodm,mnodm,lfaci,lfacb_loc) 

           !----------------------------------------------------------------
           !
           ! E. Add edge nodes with INEIG to permutation list
           !
           !----------------------------------------------------------------

           !
           ! E1. Copy old permutation array and reallocate it
           !
           kk = 0
           do ii = commd % bound_size(ineig),commd % bound_size(ineig+1)-1
              kk = kk + 1
              gisca(kk) = commd % bound_perm(ii)
           end do
           bsize = commd % bound_size(ineig+1)-commd % bound_size(ineig) + kedge 
           bsizf = 0
           do ifacb = 1,kface
              pflty = ltypb_loc(ifacb)
              newfa = elmgeo_mm_face_nodes(pflty)
              bsizf = bsizf + newfa
           end do
           bsize = bsize + bsizf
           
           call memory_alloca(par_memor,'BPERM % L','par_submsh',bperm(ineig) % l,bsize)             
           kk = 0
           do ii = commd % bound_size(ineig),commd % bound_size(ineig+1)-1
              kk = kk + 1
              bperm(ineig) % l(kk) = gisca(kk) 
           end do
           
           do iedgb = 1,kedge
              kk = kk + 1 
              bperm(ineig) % l(kk) = ledgc_loc(iedgb)
           end do
           
           do ifacb = 1,kface
              pflty = ltypb_loc(ifacb)
              newfa = elmgeo_mm_face_nodes(pflty)
              do ii = 1,newfa
                 lnode(ii) = lfacb_loc(ifacb) + ii - 1
              end do

              if( newfa > 1 ) then
                 do ii = 1,newfa
                    xcoor(1:ndime,ii) = coord(1:ndime,lnode(ii))
                 end do
                 call maths_geometrical_sort_using_coordinates(2_ip,ndime,newfa,xcoor,lnode)
              end if
              do ii = 1,newfa
                 kk = kk + 1
                 bperm(ineig) % l(kk) = lnode(ii)
              end do
           end do
           commd % bound_dim = commd % bound_dim + kedge + bsizf
        end do
        !
        ! Deallocate memory
        !
        call memory_deallo(par_memor,'COMMD % BOUND_PERM','par_submsh',commd % bound_perm)
        call memory_alloca(par_memor,'COMMD % BOUND_PERM','par_submsh',commd % bound_perm,commd % bound_dim)

        call memory_deallo(par_memor,'COMMD % BOUND_PERM','par_submsh',commd % bound_invp)
        call memory_alloca(par_memor,'COMMD % BOUND_PERM','par_submsh',commd % bound_invp,commd % bound_dim)

        if( nedgb > 0 ) then
           call memory_deallo(par_memor,'LEDG2','par_submsh',ledg2)
           call memory_deallo(par_memor,'LEDG1','par_submsh',ledg1)
        end if
        if( nfacb > 0 ) then
           call memory_deallo(par_memor,'LFACI','par_submsh',lfaci)
        end if

        commd % bound_size(1) = 1
        kk = 0
        do ineig = 1,commd % nneig
           bsize = size( bperm(ineig) % l )
           commd % bound_size(ineig+1) = commd % bound_size(ineig) + bsize
           do ii = 1,bsize
              kk                     = kk + 1
              commd % bound_perm(kk) = bperm(ineig) % l(ii)
              commd % bound_invp(kk) = bperm(ineig) % l(ii)
           end do
           call memory_deallo(par_memor,'BPERM % L','par_submsh', bperm(ineig) % l)
        end do
        call memory_deallo(par_memor,'BPERM','par_submsh', bperm)
        call memgen(3_ip,npoin,0_ip)

     end if

     call cputim(time4)

     if( INOTMASTER ) then

        !-------------------------------------------------------------------
        !
        ! F. Local renumbering
        !
        !-------------------------------------------------------------------

        call memory_alloca(par_memor,'permR','par_submsh',permR,npoin)
        !
        ! PERMR: Identify boundary nodes. Others'=-1, Own=-2
        !
        do ipoin = npoi1+1,npoi2-1
           permR(ipoin) = -1
        end do
        do ipoin = npoi2,npoi3
           permR(ipoin) = -2
        end do

        do ipoin = npoi3+1,npoin_old
           permR(ipoin) = -1
        end do
        !
        ! Renumber edge nodes
        !
        do iedgb = 1,nedgb
           iedgg = ledgc(1,iedgb)
           ipoin = ledgg(3,iedgg)
           if( ipoin < 0 ) then
              permR(-ipoin) = -2
           else
              permR( ipoin) = -1
           end if
        end do
        !
        ! Renumber face nodes
        !
        do ifacb = 1,nfacb
           ifacg = lfacb(1,ifacb)
           ipoin = lfacg(mnodm+1,ifacg)
           pflty = lfacg(mnodm+2,ifacg)
           newfa = elmgeo_mm_face_nodes(pflty)
           if( ipoin < 0 ) then
              do ii = 1,newfa
                 kpoin         = -ipoin + ii - 1                 
                 permR( kpoin) = -2
              end do
           else
              do ii = 1,newfa
                 kpoin         = ipoin + ii - 1                 
                 permR( kpoin) = -1
              end do              
           end if
        end do
        !
        ! Renumber interior nodes
        ! 
        kpoin = 0
        do ipoin = 1,npoin
           if( permR(ipoin) == 0 ) then
              kpoin        = kpoin + 1
              permR(ipoin) = kpoin
           end if
        end do
     end if

     if( INOTMASTER ) then
        !
        ! Re-renumber interior nodes using METIS
        !      
        if( 1 == 1 .and. kpoin > 0 ) then

           call memory_alloca(par_memor,'IA','par_submsh',ia,kpoin+1_ip)

           call par_intgra(kpoin,permR,ia)         ! Graph of interior nodes
           ja => gisca

           call memory_alloca(par_memor,'INVP2','par_submsh',invp2,kpoin)
           call memory_alloca(par_memor,'PERM2','par_submsh',perm2,kpoin)

           call alya2metis_METIS_NodeND(kpoin,ia,ja,perm2)

           kpoin = 0
           do ipoin = 1,npoin
              if( permR(ipoin) > 0 ) then
                 kpoin        = kpoin + 1
                 jpoin        = perm2(kpoin)
                 permR(ipoin) = jpoin
              end if
           end do

           call memory_deallo(par_memor,'invp2','par_submsh',invp2)
           call memory_deallo(par_memor,'perm2','par_submsh',perm2)

           call memgen(3_ip,1_ip,0_ip)        ! Was allocated in par_intgra()
           call memory_deallo(par_memor,'IA','par_submsh',ia)

        end if
        !
        ! Renumber own boundary nodes
        !
        npoi1 = kpoin 
        npoi2 = kpoin + 1
        do ipoin = 1,npoin
           if( permR(ipoin) == -2 ) then
              kpoin        = kpoin + 1
              permR(ipoin) = kpoin ! NEW => OLD
           end if
        end do
        npoi3 = kpoin
        !
        ! Renumber others' boundary nodes
        !
        do ipoin = 1,npoin
           if( permR(ipoin) == -1 ) then
              kpoin = kpoin + 1
              permR(ipoin) = kpoin
           end if
        end do
        !
        ! Reorder nodal arrays
        !     
        call renumbering_node_arrays(permR,npoin)
        !
        ! Deallocate memory
        !
        call memory_deallo(par_memor,'permR','par_submsh',permR)

     end if

     call cputim(time5)

     !-------------------------------------------------------------------
     !
     ! G. LNINV_LOC: Uniquely renumber the nodes. 
     !    LEINV_LOC: Uniquely renumber the element.
     !    LBINV_LOC: Uniquely renumber the boundary.
     !    This is useful for:
     !    - Postprocess
     !    - To call the division recursively
     !
     !-------------------------------------------------------------------

     call memgen(1_ip,npart_par+1,0_ip)

     call par_algath()

     npoin_total = 0
     nelem_total = 0
     nboun_total = 0

     if( IMASTER ) then
        do kfl_desti_par = 1,npart_par
           npoin_par(kfl_desti_par) = npoia_tot(kfl_desti_par+1)
           nelem_par(kfl_desti_par) = nelem_tot(kfl_desti_par+1)
           nboun_par(kfl_desti_par) = nboun_tot(kfl_desti_par+1)
           npoin_total              = npoin_total + npoin_par(kfl_desti_par)
           nelem_total              = nelem_total + nelem_par(kfl_desti_par)
           nboun_total              = nboun_total + nboun_par(kfl_desti_par)
        end do
     end if
     !
     ! Send complete mesh information to slave
     !
     isend(1) = npoin_total
     isend(2) = nboun_total
     call PAR_BROADCAST(2_ip,isend)
     npoin_total = isend(1) 
     nboun_total = isend(2) 

     do ii = 1,npart_par+1
        gisca(ii) = npoin_tot(ii)
     end do
     do ii = 2,npart_par+1
        gisca(ii) = gisca(ii-1) + gisca(ii)
     end do
     !
     ! LNINV_LOC: nodal global and unique numbering
     !
     call par_number_mesh(2_ip)
     !
     ! LEINV_LOC: element global and unique numbering
     !
     do ii = 1,npart_par+1
        gisca(ii) = nelem_tot(ii)
     end do
     do ii = 2,npart_par+1
        gisca(ii) = gisca(ii-1) + gisca(ii)  
     end do

     if( INOTMASTER ) then

        call memory_deallo(par_memor,'LEINV','par_submsh',leinv_loc)
        call memory_alloca(par_memor,'LEINV','par_submsh',leinv_loc,nelem)

        kelem = gisca(kfl_paral)
        do ielem = 1,nelem
           kelem = kelem + 1
           leinv_loc(ielem) = kelem
        end do

     end if
     !
     ! LBINV_LOC: element global and uniquex2 numbering
     !
     if( nboun_total > 0 ) then
        do ii = 1,npart_par+1
           gisca(ii) = nboun_tot(ii)
        end do
        do ii = 2,npart_par+1
           gisca(ii) = gisca(ii-1) + gisca(ii)  
        end do

        if( INOTMASTER ) then

           call memory_deallo(par_memor,'LBINV','par_submsh',lbinv_loc)
           call memory_alloca(par_memor,'LBINV','par_submsh',lbinv_loc,max(1_ip,nboun))

           kboun = gisca(kfl_paral)
           do iboun = 1,nboun
              kboun = kboun + 1
              lbinv_loc(iboun) = kboun
           end do

        end if
     end if

     call memgen(3_ip,npart_par+1,0_ip)

     !----------------------------------------------------------------------
     !
     ! Deallocate memory
     !
     !----------------------------------------------------------------------

     if( INOTMASTER ) then

        if( nedgb > 0 ) then
           call memory_deallo(par_memor,'LEDGC_LOC','par_submsh',ledgc_loc)
           call memory_deallo(par_memor,'LEDGC',    'par_submsh',ledgc)
        end if

        if( nfacb > 0 ) then
           call memory_deallo(par_memor,'LTYPB_LOC','par_submsh',ltypb_loc)
           call memory_deallo(par_memor,'LFACB_LOC','par_submsh',lfacb_loc)
           call memory_deallo(par_memor,'LFACB',    'par_submsh',lfacb)
        end if

     end if

     call cputim(time6)

     cpu_other(4) = cpu_other(4) + time1 - time0 
     cpu_other(5) = cpu_other(5) + time2 - time1 
     cpu_other(6) = cpu_other(6) + time3 - time2 
     cpu_other(7) = cpu_other(7) + time4 - time3 
     cpu_other(8) = cpu_other(8) + time5 - time4 
     cpu_other(9) = cpu_other(9) + time6 - time5 

  end if

end subroutine par_submsh

subroutine par_sorti2(n,a,b,c)

  !-----------------------------------------------------------------------
  !
  ! Sort a vector c according to a and then b
  !
  !-----------------------------------------------------------------------
  use       def_kintyp
  implicit none
  integer(ip) :: n,i,j,t,u,v
  integer(ip) :: a(n),b(n),c(n)

  do i = 1,n-1
     do j = i+1,n
        if( a(i) > a(j) ) then
           t    = a(i)
           a(i) = a(j) 
           a(j) = t
           u    = b(i)
           b(i) = b(j) 
           b(j) = u
           v    = c(i)
           c(i) = c(j) 
           c(j) = v
        end if
     end do
  end do

  do i = 1,n-1
     do j = i+1,n
        if( a(i) == a(j) ) then
           if( b(i) > b(j) ) then
              t    = a(i)
              a(i) = a(j) 
              a(j) = t
              u    = b(i)
              b(i) = b(j) 
              b(j) = u
              v    = c(i)
              c(i) = c(j) 
              c(j) = v
           end if
        end if
     end do
  end do

end subroutine par_sorti2

subroutine par_sortii(n,mdof,dof,ai,c)

  !-----------------------------------------------------------------------
  !
  ! Sort a vector c according to a and then b
  !
  !-----------------------------------------------------------------------
  use       def_kintyp
  implicit none
  integer(ip) :: n,mdof,dof,i,j,v,k,l,t
  integer(ip) :: ai(mdof,n),c(n)

  do i = 1,n-1
     do j = i+1,n
        if( ai(1,i) > ai(1,j) ) then
           do k = 1,dof
              t       = ai(k,i)
              ai(k,i) = ai(k,j) 
              ai(k,j) = t
           end do
           v    = c(i)
           c(i) = c(j) 
           c(j) = v
        end if
     end do
  end do

  do l = 1,dof-1
     do i = 1,n-1
        do j = i+1,n
           if( all(ai(1:l,i)==ai(1:l,j)) ) then
              if( ai(l+1,i) > ai(l+1,j) ) then
                 do k = 1,dof
                    t       = ai(k,i)
                    ai(k,i) = ai(k,j) 
                    ai(k,j) = t
                 end do
                 v    = c(i)
                 c(i) = c(j) 
                 c(j) = v           
              end if
           end if
        end do
     end do
  end do

end subroutine par_sortii
