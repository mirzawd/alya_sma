!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine csrgro()
  !-----------------------------------------------------------------------
  !****f* domain/csrgro
  ! NAME
  !    csrgro
  ! DESCRIPTION
  !    Set up the CSR format for deflated solvers
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_solver
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_communications, only : PAR_ALLGATHERV
  use mod_communications, only : PAR_ALLGATHER
  use mod_htable 
  use mod_graphs,         only : graphs_permut_metis_postordering
  use mod_graphs,         only : graphs_permut_metis_postordering_deallocate
  use mod_memory
  implicit none
  integer(ip)          :: igrou,ndofn,ipoin,ll,nzgro,izdom,jgrou,nlelp
  integer(ip)          :: pnode,ngrou,lsize,kk,ii,ipart,ifoun,izgro
  integer(ip), pointer :: ia(:)           
  integer(ip), pointer :: ja(:)           
  integer(ip), pointer :: lgrou(:)        
  integer(ip), pointer :: recvbuf_1(:)    
  integer(ip), pointer :: recvcounts_1(:) 
  integer(ip), pointer :: sendbuf_1(:)    
  integer(ip)          :: sendcounts_1
  integer(ip), pointer :: displs_1(:)     
  integer(ip), pointer :: recvbuf_2(:)    
  integer(ip), pointer :: recvcounts_2(:) 
  integer(ip), pointer :: sendbuf_2(:)    
  integer(ip)          :: sendcounts_2
  integer(ip), pointer :: displs_2(:)     
  type(i1p),   pointer :: ita(:)          
  logical(lg)          :: debuggingMode

  logical(lg), pointer :: touch(:)        
  integer(ip), pointer :: ngpoi(:)        
  integer(ip), pointer :: lgrpo(:)        
  integer(ip), pointer :: pgrpo(:)        

  integer(ip), pointer :: permr(:)        
  integer(ip), pointer :: invpr(:)

  real(rp),    pointer :: aa(:)

  nullify(ia)           
  nullify(ja)           
  nullify(lgrou)        
  nullify(recvbuf_1)    
  nullify(recvcounts_1) 
  nullify(sendbuf_1)    
  nullify(displs_1)     
  nullify(recvbuf_2)    
  nullify(recvcounts_2) 
  nullify(sendbuf_2)    
  nullify(displs_2)     
  nullify(ita)          
  nullify(touch)        
  nullify(ngpoi)        
  nullify(lgrpo)        
  nullify(pgrpo)        
  nullify(permr)        
  nullify(invpr)   
  nullify(aa)        
  !
  ! Deallocate memory just in case
  !
  call memory_deallo(memor_dom,'SOLVE % IAGRO','csrgro',solve_sol(1) % iagro)
  call memory_deallo(memor_dom,'SOLVE % JAGRO','csrgro',solve_sol(1) % jagro)

  ndofn =  solve_sol(1) % ndofn !TODO: this initialization was after the if, and lgrou was used without initizalitaion. Check if this is fine.
  lgrou => solve_sol(1) % lgrou
  
  if( solve_sol(1) % kfl_defas == 2 ) then
     ngrou = 0
     do ipoin = 1,npoin
        igrou = lgrou(ipoin)
        if( igrou > 0 ) ngrou = ngrou + 1
     end do
  else 
     ngrou =  solve_sol(1) % ngrou 
  end if

  debuggingMode = .false.

  call memory_alloca(memor_dom,'SOLVE % IAGRO','csrgro',solve_sol(1) % iagro,ngrou+1_ip)
  ia => solve_sol(1) % iagro

  !----------------------------------------------------------------------
  !
  ! Compute graph IA, JA
  !
  !----------------------------------------------------------------------

  if( INOTMASTER ) then
     !
     ! NGPOI(IGROU) = number of nodes per group
     ! 
     call memory_alloca(memor_dom,'NGPOI','csrgro',ngpoi,ngrou)
     do ipoin = 1,npoin
        igrou = lgrou(ipoin)
        if( igrou > 0 ) ngpoi(igrou) = ngpoi(igrou) + 1
     end do
     ! 
     ! Allocate memory for PGRPO and compute it
     ! PGRPO(IGROU): linked list, pointer of each group
     !
     call memory_alloca(memor_dom,'PGRPO','csrgro',pgrpo,ngrou+1_ip)
     pgrpo(1) = 1
     do igrou = 1,ngrou
        pgrpo(igrou+1) = pgrpo(igrou) + ngpoi(igrou)
     end do
     !
     ! Allocate memory for LGRPO and construct the list
     ! LGRPO(IGROU): linked list, list of nodes for each group
     !
     !
     nlelp = pgrpo(ngrou+1)
     call memory_alloca(memor_dom,'LGRPO','csrgro',lgrpo,nlelp)
     mpopo = 0        
     do ipoin = 1,npoin
        igrou = lgrou(ipoin)
        if( igrou > 0 ) then
           lgrpo(pgrpo(igrou)) = ipoin
           pgrpo(igrou) = pgrpo(igrou) + 1
        end if
     end do
     !
     ! Recompute PGRPO 
     !
     pgrpo(1) = 1
     do igrou = 1,ngrou
        pgrpo(igrou+1) = pgrpo(igrou) + ngpoi(igrou)
     end do
     call memory_deallo(memor_dom,'NGPOI','csrgro',ngpoi)
     !
     ! Compute graph size NZGRO
     !
     solve_sol(1) % nzgro = 0
     call memory_alloca(memor_dom,'TOUCH','csrgro',touch,ngrou)
     do igrou = 1,ngrou
        do jgrou = 1,ngrou       
           touch(jgrou) = .false.
        end do
        do izgro = pgrpo(igrou),pgrpo(igrou+1)-1
           ipoin = lgrpo(izgro)
           do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
              jgrou = lgrou(c_dom(izdom))
              if( jgrou > 0 ) then
                 if( .not. touch(jgrou) ) then
                    touch(jgrou) = .true.
                    solve_sol(1) % nzgro = solve_sol(1) % nzgro + 1
                 end if
              end if
           end do
        end do
     end do
    !
     ! Construct the array of indexes
     !
     call memory_alloca(memor_dom,'SOLVE % JAGRO','csrgro',solve_sol(1) % jagro,solve_sol(1) % nzgro)
     ja => solve_sol(1) % jagro

     solve_sol(1) % nzgro = 1
     do igrou = 1,ngrou
        do jgrou = 1,ngrou       
           touch(jgrou) = .false.
        end do
        ia(igrou) = solve_sol(1) % nzgro
        do izgro = pgrpo(igrou),pgrpo(igrou+1)-1
           ipoin = lgrpo(izgro)
           do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
              jgrou = lgrou(c_dom(izdom))
              if( jgrou > 0 ) then
                 if( .not. touch(jgrou) ) then
                    touch(jgrou) = .true.
                    ja(solve_sol(1) % nzgro) = jgrou
                    solve_sol(1) % nzgro = solve_sol(1) % nzgro + 1
                 end if
              end if
           end do
        end do
     end do

     solve_sol(1) % nzgro = solve_sol(1) % nzgro - 1
     ia(ngrou+1)          = solve_sol(1) % nzgro + 1

     call memory_deallo(memor_dom,'TOUCH','csrgro',touch)
     call memory_deallo(memor_dom,'PGRPO','csrgro',pgrpo)
     call memory_deallo(memor_dom,'LGRPO','csrgro',lgrpo)
  else

     solve_sol(1) % nzgro = 0
     nzgro = 0
     call memory_alloca(memor_dom,'SOLVE % JAGRO','csrgro',solve_sol(1) % jagro,1_ip)

  end if

  !----------------------------------------------------------------------
  !
  ! Merge the slave graphs in Parallel
  !
  !----------------------------------------------------------------------
 
  if( IPARALL .and. solve_sol(1) % kfl_defas == 1 ) then
     !
     ! Debug
     !
     if( debuggingMode ) then
        if( INOTMASTER ) then
           write(90+kfl_paral,'(a)') 'MY GRAPH:'
           do igrou = 1,ngrou
              write(90+kfl_paral,'(i3,a1,20(1x,i3))') &
                   igrou,'=',(ja(izdom),izdom=ia(igrou),ia(igrou+1)-1) 
           end do
        end if
     end if
     !
     ! PARIS, PARIN
     !
     izdom = 0   
     if( INOTMASTER ) then 
        do igrou = 1,ngrou
           if( ia(igrou+1)-ia(igrou) > 0 ) &
                izdom = izdom + 1
        end do
     end if
     call memory_alloca(memor_dom,'RECVCOUNTS_1','csrgro',recvcounts_1,npart+1_ip)    
     call PAR_ALLGATHER(izdom,recvcounts_1)
     !
     ! Give my graph size through allgatherv: PARIS, PARIN, PARIG, PARI1
     !
     kk = 0
     do ii = 1,npart+1
        recvcounts_1(ii) = recvcounts_1(ii) * 2
        kk               = kk + recvcounts_1(ii)
     end do
     if( kk > 0 ) then
        call memory_alloca(memor_dom,'DISPLS_1','csrgro',displs_1,npart+1_ip) 
        displs_1(1) = 0
        do ii = 2,npart+1
           displs_1(ii) = displs_1(ii-1) + recvcounts_1(ii-1)
        end do
        call memory_alloca(memor_dom,'RECVBUF_1','csrgro',recvbuf_1,kk)
        sendcounts_1 = 2 * izdom
        if( sendcounts_1 == 0 ) then
        !   sendbuf_1 => dummi_par
        else
           call memory_alloca(memor_dom,'SENDBUF_1','csrgro',sendbuf_1,2*izdom)
           ll = 0
           do igrou = 1,ngrou
              if( ia(igrou+1)-ia(igrou) > 0 ) then
                 ll            = ll + 1
                 sendbuf_1(ll) = ia(igrou+1)-ia(igrou)
                 ll            = ll + 1
                 sendbuf_1(ll) = igrou
              end if
           end do
        end if
        call PAR_ALLGATHERV(sendbuf_1,recvbuf_1,recvcounts_1,'IN MY CODE',displs_1)
        !paris => sendbuf_1
        !npasi =  sendcounts_1
        !parin => recvbuf_1  
        !parig => recvcounts_1
        !pari1 => displs_1
        !call Parall(49_ip)
     end if
     !
     ! Give my graph through allgatherv: PARIS, PARIN, PARIG, PARI1
     !
     call memory_alloca(memor_dom,'RECVCOUNTS_2','csrgro',recvcounts_2,npart+1_ip)
     do ipart = 1,npart+1
        recvcounts_2(ipart) = 0
     end do

     ll = 0
     do ipart = 1,npart+1
        do ii = 1,recvcounts_1(ipart)/2
           recvcounts_2(ipart) = recvcounts_2(ipart) + recvbuf_1(ll+1)
           ll                  = ll + 2
        end do
     end do

     if( ll > 0 ) then
        call memory_alloca(memor_dom,'DISPLS_2','csrgro',displs_2,npart+1_ip)    
        displs_2(1) = 0
        kk = 0
        do ii = 2,npart+1
           kk = kk + recvcounts_2(ii)
           displs_2(ii) = displs_2(ii-1) + recvcounts_2(ii-1)
        end do
        call memory_alloca(memor_dom,'RECVBUF_2','csrgro',recvbuf_2,kk)
        sendcounts_2 = solve_sol(1) % nzgro

        if( sendcounts_2 == 0 ) then
           !sendbuf_2 => dummi_par
        else
           call memory_alloca(memor_dom,'SENDBUF_2','csrgro',sendbuf_2,sendcounts_2)
           ll = 0
           do igrou = 1,ngrou
              do izdom = ia(igrou),ia(igrou+1)-1
                 ll = ll + 1
                 sendbuf_2(ll) = ja(izdom)
              end do
           end do
        end if
        call PAR_ALLGATHERV(sendbuf_2,recvbuf_2,recvcounts_2,'IN MY CODE',displs_2)
        !paris => sendbuf_2
        !npasi =  sendcounts_2
        !parin => recvbuf_2
        !parig => recvcounts_2
        !pari1 => displs_2
        !call Parall(49_ip)
     end if
     !
     ! Debugging mode
     !
     if( debuggingMode ) then
        if( IMASTER ) then
           write(90,'(a)') 'DISPL_1'
           write(90,*) displs_1    
           write(90,'(a)') 'DISPL_2'
           write(90,*) displs_2
           write(90,'(a)') 'RECVBUF_1'
           do ipart = 2,npart+1
              write(90,'(a,i3)') 'List subdomain ',ipart-1
              ll = displs_1(ipart)
              ii = 0
              do while( ii < recvcounts_1(ipart) )
                 write(90,'(i3,1x,i3)') recvbuf_1(ll+1),recvbuf_1(ll+2)
                 ii = ii + 2
                 ll = ll + 2
              end do
           end do
        end if
     end if
     !
     ! Merge the graphs
     !
     call memory_alloca(memor_dom,'ITA','csrgro',ita,ngrou)
     call memgen(1_ip,ngrou,0_ip)
     do igrou = 1,ngrou
        lsize = 0
        do izdom = ia(igrou),ia(igrou+1)-1
           lsize = lsize + 1
           gisca(lsize) = ja(izdom)
        end do
        do ipart = 2,npart+1
           if( kfl_paral /= ipart - 1 ) then
              ll    = displs_1(ipart)
              ii    = 0
              kk    = displs_2(ipart)+1
              ifoun = 0
              do while( ifoun == 0 .and. ii < recvcounts_1(ipart) )
                 if( recvbuf_1(ll+2) == igrou ) then
                    ifoun = 1
                 else
                    kk = kk + recvbuf_1(ll+1)
                    ll = ll + 2
                    ii = ii + 2
                 end if
              end do
              if( ifoun == 1 ) then
                 pnode = recvbuf_1(ll+1)
                 if( debuggingMode ) then
                    write(90+kfl_paral,'(a,i3,a,i3)') 'SUBDOMAIN ',ipart-1,&
                         ' gives me the following graph for group ',igrou
                    write(90+kfl_paral,'(a,i3)') 'It gives me # nodes= ',pnode
                    write(90+kfl_paral,'(20(1x,i3))') recvbuf_2(kk:kk+pnode-1)
                 end if
                 call mergli( gisca, lsize, recvbuf_2(kk) , pnode, -1_ip )
              end if
           end if
        end do
        if( debuggingMode ) then
           if( kfl_paral == 1 ) then
              write(90+kfl_paral,'(a,i3)') 'Connectivity for group ',igrou
              write(90+kfl_paral,'(20(1x,i3))') gisca(1:ngrou)
           end if
        end if
        lsize = 1
        giscaloop: do while( gisca(lsize) /= 0 ) 
           lsize = lsize + 1
           if( lsize > ngrou ) exit giscaloop
        end do giscaloop
        lsize = lsize - 1
        call memory_alloca(memor_dom,'ITA % L','csrgro',ita(igrou) % l,lsize)
        do ll = 1,lsize
           ita(igrou) % l(ll) = gisca(ll)
           gisca(ll)          = 0
        end do
     end do
     !
     ! Deallocate memory
     !
     call memgen(3_ip,ngrou,0_ip)
     call memory_deallo(memor_dom,'recvcounts_1','csrgro',recvcounts_1 )
     call memory_deallo(memor_dom,'displs_1    ','csrgro',displs_1     )
     call memory_deallo(memor_dom,'recvbuf_1   ','csrgro',recvbuf_1    )
     call memory_deallo(memor_dom,'sendbuf_1   ','csrgro',sendbuf_1    )
     call memory_deallo(memor_dom,'recvcounts_2','csrgro',recvcounts_2 )
     call memory_deallo(memor_dom,'displs_2    ','csrgro',displs_2     )
     call memory_deallo(memor_dom,'recvbuf_2   ','csrgro',recvbuf_2    )
     call memory_deallo(memor_dom,'sendbuf_2   ','csrgro',sendbuf_2    )
     !
     ! Reconstruct graph: IA, JA
     !
     call memory_deallo(memor_dom,'SOLVE % JAGRO','csrgro',solve_sol(1) % jagro )
     nzgro = 0
     do igrou = 1,ngrou
        nzgro = nzgro + size( ita(igrou) % l , KIND=ip )
     end do
     solve_sol(1) % nzgro = nzgro
     call memory_alloca(memor_dom,'SOLVE % JAGRO','csrgro',solve_sol(1) % jagro,nzgro )
     ja => solve_sol(1) % jagro
     ia(1) = 1
     do igrou = 1,ngrou
        ia(igrou+1) = ia(igrou) + size( ita(igrou) % l , KIND=ip )
     end do
     nzgro = 0
     do igrou = 1,ngrou
        lsize = 0
        do izdom = ia(igrou),ia(igrou+1)-1
           lsize = lsize + 1
           ja(izdom) = ita(igrou) % l(lsize)
        end do
        call heapsorti1(2_ip,lsize,ja(ia(igrou)))
        !deallocate( ita(igrou) % l )
     end do
     call memory_deallo(memor_dom,'ITA','csrgro',ita )

     !if( associated(ita) ) deallocate( ita ) 

     if( debuggingMode ) then
        write(90+kfl_paral,'(a)') 'GRAPH ---------------------------------------' 
        do igrou = 1,ngrou
           write(90+kfl_paral,'(a,i3,a,20(1x,i3))') 'Connectivity for group ',&
                igrou,': ',(ja(izdom),izdom=ia(igrou),ia(igrou+1)-1 )
        end do
        call runend('SWSWSW')
     end if
  end if

  !----------------------------------------------------------------------
  !
  ! Renumber graph using METIS and postordering to minimize filling
  !
  !----------------------------------------------------------------------

  nzgro = solve_sol(1) % nzgro

#ifdef __PGI
#ifdef I8
#define PGI8
#endif
#endif

#ifdef PGI8

  call memory_alloca(memor_dom,'PERMR_POST','csrgro',permr,ngrou)
  call memory_alloca(memor_dom,'INVPR_POST','csrgro',invpr,ngrou)
  do igrou = 1,ngrou
    permr(igrou) = igrou
    invpr(igrou) = igrou
  end do

#else

  call graphs_permut_metis_postordering(ngrou,nzgro,ia,ja,permr,invpr,&
       PERMR_NAME='PERMR_CSRGRO',INVPR_NAME='INVPR_CSRGRO',memor=memor_dom)

#endif

  !----------------------------------------------------------------------
  !
  ! Renumber groups
  !
  !----------------------------------------------------------------------

  if( INOTMASTER ) then
     lgrou => solve_sol(1) % lgrou
     do ipoin = 1,npoin
        jgrou = lgrou(ipoin)
        if( jgrou > 0 ) then
           igrou = permr(jgrou) ! new = permr(old)
           lgrou(ipoin) = igrou
        end if
     end do
  end if
  
  call graphs_permut_metis_postordering_deallocate(permr,invpr,&
       memor_dom,PERMR_NAME='PERMR_CSRGRO',INVPR_NAME='INVPR_CSRGRO')

  !----------------------------------------------------------------------
  !
  ! Postprocess matrix
  !
  !----------------------------------------------------------------------

  !if( INOTSLAVE ) then
  !   allocate(aa(nzgro))
  !   aa = 1.0_rp
  !   call pspltm(&
  !        ngrou,ngrou,1_ip,0_ip,ja,ia,aa,&
  !        trim(title)//': '//namod(modul),0_ip,18.0_rp,'cm',&
  !        0_ip,0_ip,2_ip,99)
  !   deallocate(aa)
  !end if

end subroutine csrgro
