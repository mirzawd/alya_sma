!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_disbou(memor)
  !-------------------------------------------------------------------------------
  !****f* Parall/par_disbou
  ! NAME
  !    par_disbou
  ! DESCRIPTION
  !    This routine distributes the boundary points between adjacent subdomains
  !    When entering, LNPAR_PAR=0 for boundary nodes
  !                            =subdomain #
  ! INPUT
  !    Element graph
  ! OUTPUT
  !
  ! USED BY
  !    par_create_domain_graph
  !***
  !-------------------------------------------------------------------------------
  use def_parame
  use def_kintyp
  use def_domain
  use def_parall
  use mod_memory
  use mod_alya2metis, only : alya2metis_METIS_PartGraph
  implicit none
  integer(8)           :: memor(2)
  integer(ip), pointer :: indice_dom(:) 
  integer(ip), pointer :: xfrontera(:)  
  integer(ip), pointer :: i_front(:)    
  integer(ip), pointer :: d_front(:)    
  integer(ip), pointer :: domli(:)      
  integer(ip), pointer :: permI(:)      
  integer(ip), pointer :: invpI(:)      
  integer(ip), pointer :: permB(:)      
  integer(ip), pointer :: invpB(:)      
  integer(ip), pointer :: xadjSubDom(:) 
  integer(ip), pointer :: adjSubDom(:)  
  integer(ip), pointer :: lnpar_loc(:)  
  integer(ip), pointer :: wvert_loc(:)  
  integer(ip)          :: wedge_loc, weigh_loc
  integer(ip)          :: nbNodInter, nbNodBound, ierro
!  integer(ip)          :: dummi
  integer(ip)          :: ii, jj, kk, dom1, dom2, n_nei, n_nei2
  integer(ip)          :: inode, ndomi, mascara, maxb,max2,adjsize
  integer(ip)          :: dsiz1, dsiz2

  nullify(indice_dom )
  nullify(xfrontera  )
  nullify(i_front    )
  nullify(d_front    )
  nullify(domli      )
  nullify(permI      )
  nullify(invpI      )
  nullify(permB      )
  nullify(invpB      )
  nullify(xadjSubDom )
  nullify(adjSubDom  )
  nullify(lnpar_loc  )
  nullify(wvert_loc  )
  
  call memory_alloca(memor,'INDICE_DOM','par_disbou' , indice_dom , npart_par      )
  call memory_alloca(memor,'XFRONTERA' ,'par_disbou' , xfrontera  , npart_par+1_ip )

  do dom1 = 1, npart_par
     indice_dom(dom1) = 0
  enddo

  maxb = 0
  do dom1 = 1, npart_par
     do dom2 = 1, dom1-1
        kk = (dom1*(dom1-1))/2 + dom2
        n_nei = neighDom(kk)
        if (n_nei> maxb) maxb = n_nei
        if (n_nei/=0) then
           indice_dom(dom1) = indice_dom(dom1) + n_nei
        endif
     enddo
  enddo

  xfrontera(1) = 1
  do dom1 = 1, npart_par
     xfrontera(dom1+1) = xfrontera(dom1) + indice_dom(dom1)
  enddo

  call memory_alloca(memor,'I_FRONT' ,'par_disbou' , i_front , xfrontera(npart_par+1_ip) )
  call memory_alloca(memor,'D_FRONT' ,'par_disbou' , d_front , xfrontera(npart_par+1_ip) )
  call memory_alloca(memor,'DOMLI'   ,'par_disbou' , domli   , mepoi )

  do inode = 1, npoin
     if (lnpar_par(inode)==0) then
        !
        ! INODE is a boundary node
        ! 
        call par_domlis( pelpo, lelpo, inode, lepar_par, ndomi, domli )
        do ii = 1, ndomi
           dom1 = domli(ii)
           do jj = ii+1, ndomi
              dom2 = domli(jj)
              if (dom2> dom1) then
                 kk              = xfrontera(dom2)
                 d_front(kk)     = dom1
                 i_front(kk)     = inode
                 xfrontera(dom2) = kk + 1
              else
                 kk              = xfrontera(dom1)
                 d_front(kk)     = dom2
                 i_front(kk)     = inode
                 xfrontera(dom1) = kk + 1
              endif
           enddo
        enddo
     endif
  enddo

  call memory_deallo(memor,'DOMLI'   ,'par_disbou' , domli )

  xfrontera(1) = 1
  do dom1 = 1, npart_par
     xfrontera(dom1+1) = xfrontera(dom1) + indice_dom(dom1)
  enddo
  max2  = (r_dom(npoin+1)/npoin)*maxb*2
  dsiz1 = maxb+1_ip
  dsiz2 = max2

  call memory_alloca(memor,'PERMI'     ,'par_disbou' , permI      , npoin )
  call memory_alloca(memor,'PERMB'     ,'par_disbou' , permB      , npoin )
  
  10 continue

  call memory_alloca(memor,'INVPB'     ,'par_disbou' , invpB      , maxb  )
  call memory_alloca(memor,'INVPI'     ,'par_disbou' , invpI      , maxb  )
  
  call memory_alloca(memor,'XADJSUBDOM','par_disbou' , xadjSubDom , dsiz1 )
  call memory_alloca(memor,'ADJSUBDOM' ,'par_disbou' , adjSubDom  , dsiz2 )
  !
  ! Partition graph using METIS
  !
  weigh_loc = 0
  wedge_loc = 0
  ierro     = 0
  mascara   = npart_par+1
  ierro     = 0
  
  do dom1 = 1, npart_par
     do dom2 = 1, dom1-1
        kk = (dom1*(dom1-1))/2 + dom2
        n_nei = neighDom(kk)
        if (n_nei/=0) then
           n_nei2 = 0
           do ii = xfrontera(dom1), xfrontera(dom1+1)-1
              if (d_front(ii)==dom2) then
                 inode = i_front(ii)
                 if (lnpar_par(inode)==0) then
                    lnpar_par(inode) = mascara
                    n_nei2 = n_nei2 + 1
                 endif
              endif
           enddo
           if( n_nei2 > 0 ) then
              !
              !create a subgraph
              !
              call par_subgra( &
                   npoin , r_dom, c_dom, mascara, lnpar_par,        &
                   nbNodInter, nbNodBound, xadjSubDom, adjSubDom,   &
                   permI, invpI, permB, invpB , adjsize , dsiz1 ,   &
                   dsiz2, maxb, ierro  )
              
              if( ierro /= 0 ) then
                 dsiz1 = dsiz1 * 2
                 dsiz2 = dsiz2 * 2
                 maxb  = maxb  * 2 
                 call memory_deallo(memor,'XADJSUBDOM','par_disbou' , xadjSubDom )
                 call memory_deallo(memor,'ADJSUBDOM' ,'par_disbou' , adjSubDom  )
                 goto 10
              end if
              
              call memory_alloca(memor,'LNPAR_LOC' ,'par_disbou' , lnpar_loc  , nbNodInter )

              if( nbNodInter > 0 ) then
                 if( kfl_interface_parti_par == 1 ) then
                    allocate(wvert_loc(nbNodInter))
                    wvert_loc = 1_ip
                    call alya2metis_METIS_PartGraph(2_ip,nbNodInter,xadjSubDom,adjSubDom,wvert_loc,lnpar_loc)
                    deallocate(wvert_loc)                    
                    !call par_metis(&
                    !     1_ip , nbNodInter , maxb, xadjSubDom, adjSubDom, &
                    !     wvert_loc, wedge_loc, weigh_loc,  &
                    !     2_ip , lnpar_loc , dummi, dummi, &
                    !     dummi, dummi, dummi, dummi, 1_ip , 1_ip , memor )
                 else if( kfl_interface_parti_par == 2 ) then
                    lnpar_loc(1) = 1
                    if( nbNodInter > 1 ) then
                       lnpar_loc(1:nbNodInter/2) = 1
                       lnpar_loc(nbNodInter/2+1:nbNodInter) = 2                 
                    end if
                 end if
              end if

              do ii= 1, nbNodInter
                 inode = invpI(ii)
                 if (lnpar_loc(ii)==1) then
                    lnpar_par(inode) = -dom1
                 else
                    lnpar_par(inode) = -dom2
                 endif
              enddo
              call memory_deallo(memor,'LNPAR_LOC' ,'par_disbou' , lnpar_loc  )
              mascara =  mascara + 1
           endif
        endif
     enddo
  enddo
  !
  ! Deallocate memory
  !
  call memory_deallo(memor,'ADJSUBDOM' ,'par_disbou' , adjSubDom  )
  call memory_deallo(memor,'XADJSUBDOM','par_disbou' , xadjSubDom )
  call memory_deallo(memor,'INVPB'     ,'par_disbou' , invpB      )
  call memory_deallo(memor,'PERMB'     ,'par_disbou' , permB      )
  call memory_deallo(memor,'INVPI'     ,'par_disbou' , invpI      )
  call memory_deallo(memor,'PERMI'     ,'par_disbou' , permI      )
  call memory_deallo(memor,'I_FRONT'   ,'par_disbou' , i_front    )
  call memory_deallo(memor,'D_FRONT'   ,'par_disbou' , d_front    )
  call memory_deallo(memor,'XFRONTERA' ,'par_disbou' , xfrontera  )
  call memory_deallo(memor,'INDICE_DOM','par_disbou' , indice_dom )

end subroutine par_disbou
