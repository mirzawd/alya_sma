!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_arrays()
  !-------------------------------------------------------------------------------
  !****f* Parall/par_arrays
  ! NAME
  !    par_arrays
  ! DESCRIPTION
  !
  ! INPUT
  !    Element graph
  ! OUTPUT
  !    Partition of the graph
  !    lepar_par
  !    leper_par
  !    leinv_par
  !    lnpar_par
  !    lnper_par
  !    lninv_par
  !    lneig_par
  !    ginde_par 
  !    lcomm_par
  ! USED BY
  !    par_partit
  !***
  !-------------------------------------------------------------------------------
  use def_parame 
  use def_domain 
  use def_parall
  use def_master
  use mod_memory
  use def_coupli,         only : kfl_graph_cou
  use mod_domain,         only : domain_memory_deallocate
  use mod_parall, only : par_memor
  implicit none
  integer(ip)             :: ii,jj,kk,front
  integer(ip)             :: dummi,ndual_par,offsetB,offsetI
  integer(ip)             :: adjsize,nbNodInter,ipart,jpart
  integer(ip)             :: nbNodBound,inter,dsize,dsiz2,dsiz1
  integer(ip)             :: ielem,jelem,iboun,jboun,inode,ierro
  integer(ip)             :: b_ind,ndomi,domai,jnode
  integer(ip)             :: domli(mepoi)
  integer(ip), pointer    :: permR(:)
  integer(ip), pointer    :: invpR(:)
  integer(ip), pointer    :: permI(:)
  integer(ip), pointer    :: invpI(:)
  integer(ip), pointer    :: permB(:)
  integer(ip), pointer    :: invpB(:)
  integer(ip), pointer    :: xadjSubDom(:)
  integer(ip), pointer    :: adjSubDom(:)
  integer(ip), pointer    :: indice_dom(:)
  real(rp)                :: time2
  real(rp)                :: time3,time4
  real(rp)                :: time5,time6
  real(rp)                :: time7,time8

#ifdef EVENT
  call mpitrace_user_function(1)
#endif
  !
  ! Nullify pointers
  !
  nullify(permR)
  nullify(invpR)
  nullify(permI)
  nullify(invpI) 
  nullify(permB)
  nullify(invpB)
  nullify(xadjSubDom) 
  nullify(adjSubDom)
  nullify(indice_dom)
  !
  ! Deallocate memory of PADJA_PAR (PELEL) and LADJA_PAR (LELEL). They were allocated in par_elmgra
  !
  call cputim(time2)
  if(kfl_graph_cou == 0) then
     call domain_memory_deallocate('PELEL')
     call domain_memory_deallocate('LELEL')
  endif
  ! 
  ! Compute communication strategy
  !  
  call par_livinf(4_ip,' ',dummi)
  call par_domgra(inter,front,par_memor)     ! LNPAR_PAR, NEIGHDOM, LNEIG_PAR, ADJDOM, XADJDOM
  call cputim(time3)
  cpu_paral(6) = time3 - time2
  
  call par_duagra(ndual_par,par_memor)       ! NBDUAL, TRANSLDUAL, TRANSL, IADUAL, JADUAL
  call cputim(time4)
  cpu_paral(7) = time4 - time3

  call par_colgra(ndual_par,nbcol,par_memor) ! NBCOLOURS, COLOUR
  call cputim(time5)
  cpu_paral(8)= time5 - time4

  call par_commun(nbcol,ndual_par)                      ! COMMUNSORT, LCOMM_PAR
  call cputim(time6)
  cpu_paral(9)= time6 - time5
  !
  ! Compute boundary partition LBPAR_PAR and NBOUN_PAR
  !
  call par_bounda()
  !
  ! Deallocate memory
  !
  call par_memory(7_ip)
  !
  ! Number interior and boundary elements for each subdomain
  !  
  dsize = 4_ip*(npoin/npart_par)               ! OPTIMIZAR or DE/ALLOCATE EACH TIME
  dsiz1 = dsize + 1_ip                         ! OPTIMIZAR or DE/ALLOCATE EACH TIME
  dsiz2 = (r_dom(npoin+1)/npoin)*dsize*1_ip    ! OPTIMIZAR or DE/ALLOCATE EACH TIME

  call memory_alloca( par_memor , 'PERMI'      , 'par_arrays' , permI      , npoin )
  call memory_alloca( par_memor , 'PERMB'      , 'par_arrays' , permB      , npoin )

10 continue

  call memory_alloca( par_memor , 'INVPB'      , 'par_arrays' , invpB      , dsize )
  call memory_alloca( par_memor , 'INVPI'      , 'par_arrays' , invpI      , dsize )
  call memory_alloca( par_memor , 'XADJSUBDOM' , 'par_arrays' , xadjSubDom , dsiz1 )
  call memory_alloca( par_memor , 'ADJSUBDOM'  , 'par_arrays' , adjSubDom  , dsiz2 )

  gni     = 0_ip
  gnb     = 0_ip
  offsetI = 0_ip
  offsetB = inter
  ierro   = 0_ip
  
  call par_livinf(18_ip,' ',dummi)

  do ipart = 1,npart_par
     !
     ! For each subdomain: separate interior and boundary nodes 
     ! and create subgraph
     !
     call par_subgra( &                                       ! NBNODINTER, NBNODBOUND
          npoin , r_dom, c_dom, ipart, lnpar_par,        &    ! XADJSUBDOM, ADJSUBDOM, ADJSIZE
          nbNodInter, nbNodBound, xadjSubDom, adjSubDom, &    ! PERMI, INVPI, PERMB, INVPB
          permI, invpI, permB, invpB , adjsize ,dsiz1,   &
          dsiz2, dsize, ierro)

     if( ierro > 0 ) then
        dsiz1 = 2_ip * dsiz1
        dsiz2 = 2_ip * dsiz2
        dsize = 2_ip * dsize
        call memory_deallo( par_memor , 'ADJSUBDOM'  , 'par_arrays' , adjSubDom  )
        call memory_deallo( par_memor , 'XADJSUBDOM' , 'par_arrays' , xadjSubDom )
        call memory_deallo( par_memor , 'INVPI'      , 'par_arrays' , invpI      )
        call memory_deallo( par_memor , 'INVPB'      , 'par_arrays' , invpB      )
        goto 10      
     end if
     
     ginde_par(1,ipart) = gni + 1
     ginde_par(2,ipart) = gnb + 1
     ginde_par(3,ipart) = nbNodInter
     ginde_par(4,ipart) = nbNodBound
     gni                = gni + nbNodInter
     gnb                = gnb + nbNodBound

     do ii = 1,xadjSubDom(nbNodInter+1)-1 
        if( adjSubDom(ii) < 1 .or. adjSubDom(ii) > nbNodInter ) &
             call runend('PAR_METIS: adjSubDom index out of range')
     end do
     !
     ! Renumber interior nodes
     ! 
     if( nbNodInter > 0 ) then
        call memory_alloca( par_memor , 'PERMR' , 'par_arrays' , permR , nbNodInter )
        call memory_alloca( par_memor , 'INVPR' , 'par_arrays' , invpR , nbNodInter )        
        do ii = 1,nbNodInter
           permR(ii) = ii
           invpR(ii) = ii
        end do
     end if
     !
     ! Number interior nodes
     !
     do ii = 1,nbNodInter
        kk            = invpI(ii)
        jj            = permR(ii) + offsetI
        lnper_par(kk) = jj
        lninv_par(jj) = kk
     end do
     offsetI = offsetI + nbNodInter 
     !
     ! Number own boundary nodes
     !
     do ii = 1,nbNodBound
        kk            = invpB(ii)
        jj            = ii + offsetB
        lnper_par(kk) = jj
        lninv_par(jj) = kk
     end do
     offsetB = offsetB + nbNodBound

     if( nbNodInter > 0 ) then
        call memory_deallo(par_memor,'PERMR' ,'par_arrays' , permR )
        call memory_deallo(par_memor,'INVPR' ,'par_arrays' , invpR )
     end if

  end do

  call cputim(time7)
  cpu_paral(11) = time7 - time6
  !
  ! End of GINDE_PAR
  !
  ginde_par( 1, npart_par + 1 ) = gni + 1
  ginde_par( 2, npart_par + 1 ) = gnb + 1
  ginde_par( 3, npart_par + 1 ) = 0
  ginde_par( 4, npart_par + 1 ) = 0
  !
  ! Deallocate memory
  !
  call memory_deallo(par_memor,'ADJSUBDOM' ,'par_arrays' , adjSubDom  )
  call memory_deallo(par_memor,'XADJSUBDOM','par_arrays' , xadjSubDom )
  call memory_deallo(par_memor,'INVPB'     ,'par_arrays' , invpB      )
  call memory_deallo(par_memor,'PERMB'     ,'par_arrays' , permB      )
  call memory_deallo(par_memor,'INVPI'     ,'par_arrays' , invpI      )
  call memory_deallo(par_memor,'PERMI'     ,'par_arrays' , permI      )
  !
  ! Compute subdomains dimensions
  !
  do ipart = 1,npart_par
     nelem_par( ipart ) = 0
  end do
  do ielem = 1,nelem
     ipart            = lepar_par(ielem)
     nelem_par(ipart) = nelem_par(ipart) + 1
  end do
  !
  ! Compute permutation vectors for nodes: lninv_loc, xlnin_loc
  !
  call par_memory(2_ip)
  call memory_alloca(par_memor,'INDICE_DOM','par_arrays' ,indice_dom , npart_par+1_ip )

  xlnin_loc(1)  = 1
  indice_dom(1) = xlnin_loc(1)
  do domai = 1,npart_par
     xlnin_loc(domai+1)  = xlnin_loc(domai) + npoin_par(domai)
     indice_dom(domai+1) = xlnin_loc(domai+1)
  end do
  !
  ! Internal points
  !
  inode = 1 
  do ipart= 1, npart_par
     do ii= 1, ginde_par(3,ipart)
        jnode                       = lninv_par(inode)
        lninv_loc(xlnin_loc(ipart)) = jnode
        xlnin_loc(ipart)            = xlnin_loc(ipart) +  1
        inode                       = inode + 1
     end do
  end do
  !
  ! Domains of every boundary nodes: BADJ, BDOM, BPOIN
  ! BADJ is the index vector for boundary nodes (1->GNB) that points to BDOM.
  ! BDOM contains the list of elements to which the boundary nodes belong to.
  ! BPOIN contains the local numbering of the global boundary node in the
  ! gnb_total numbering
  !
  call par_livinf(5_ip,' ',dummi)
  call par_memory(3_ip)
  !
  ! At this step XLNIN_LOC(IPART) indicates where the list of boundary nodes begins
  ! for subdomain IPART in the npoin_total numbering
  !
  badj(1) = 1
  b_ind   = 1
  inode   = gni
  kk      = 1

  do ipart= 1,npart_par
     !
     ! Begining of self boundary
     !
     slfbo_par(ipart) = xlnin_loc(ipart) - indice_dom(ipart) + 1
     do ii= 1,ginde_par(4,ipart)
        inode = inode + 1        ! Interior/boundary numbering (strarting at end of interior)
        jnode = lninv_par(inode) ! Original numbering
        call par_domlis( pelpo, lelpo, jnode, lepar_par, ndomi, domli )
        do jj= 1, ndomi
           jpart                       = domli(jj)
           lninv_loc(xlnin_loc(jpart)) = jnode
           bdom(b_ind)                 = jpart
           bpoin(b_ind)                = xlnin_loc(jpart) - indice_dom(jpart) + 1
           b_ind                       = b_ind + 1
           xlnin_loc(jpart)            = xlnin_loc(jpart) +  1
        end do
        kk       = kk + 1
        badj(kk) = b_ind
     end do
  end do
  do domai = 1,npart_par+1
     xlnin_loc(domai) = indice_dom(domai)
  end do
  !
  ! Element permutation: LEPER_PAR, LEINV_PAR, LEIND_PAR
  !
  indice_dom(1) = 1
  do domai= 1,npart_par
     indice_dom(domai+1)  = indice_dom(domai) + nelem_par(domai)
  end do
  do domai = 1,npart_par+1
     leind_par(domai)  = indice_dom(domai)
  end do
  do ielem = 1,nelem
     domai             = lepar_par(ielem)
     jelem             = indice_dom(domai)
     leinv_par(jelem)  = ielem
     leper_par(ielem)  = jelem 
     indice_dom(domai) = indice_dom(domai) + 1
  end do
  !
  ! Boundary permutation: LBPER_PAR, LBINV_PAR, LBIND_PAR
  !
  indice_dom(1) = 1
  do domai= 1,npart_par
     indice_dom(domai+1) = indice_dom(domai) + nboun_par(domai)
  end do

  do domai = 1,npart_par+1
     lbind_par(domai)  = indice_dom(domai)
  end do
  do iboun = 1,nboun
     domai             = lbpar_par(iboun)
     jboun             = indice_dom(domai)
     lbinv_par(jboun)  = iboun
     lbper_par(iboun)  = jboun
     indice_dom(domai) = indice_dom(domai) + 1
  end do
  !
  ! Deallocate memory
  !
  call par_memory(-16_ip)
  call memory_deallo(par_memor,'INDICE_DOM','par_arrays' ,indice_dom )

  call cputim(time8)
  cpu_paral(12)=time8-time7

#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine par_arrays
