!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_elmgra()
  !------------------------------------------------------------------------
  !****f* domain/elmgra
  ! NAME
  !    elmgra
  ! DESCRIPTION
  !    This routine computes the element graph using:
  !    ITASK=1 ... The nodal connectivity
  !             =2 ... The face connectivity
  !
  !    For example, given the following mesh, element 5 will have the
  !    following neighbors:
  !    +---+---+---+
  !    | 1 | 2 | 3 |   ITASK=1: 1,2,3,4,6,7,8,9
  !    +---+---+---+   ITASK=2: 2,4,6,8
  !    | 4 | 5 | 6 |
  !    +---+---+---+
  !    | 7 | 8 | 9 |
  !    +---+---+---+
  !
  !    Working arrays:
  !    NEPOI:      Number of occurencies of a node in connectivity
  !                Node ipoin belongs to nepoi(ipoin) elements
  !    PELPO:      Pointer to node connectivities arrays lelpo
  !    LELPO:      List of node to element connectivities
  !    MELEL:      Maximum number of neighbors=max(nepoi)*mnode
  !
  ! OUTPUT
  !    NEDGT:  Number of edges
  !    PELEL:  Pointer to the lelel array of size 
  !                pelel(nelem+1):
  !                - Adjacencies of ielem starts at pelel(ielem) 
  !                  included
  !                - Adjacencies of ielem ends at pelel(ielem+1)-1 
  !                  included
  !                - pelel does not include the element itself
  !    LELEL:  List of adjacencies.
  !
  !    Example:
  !                 1---2---3---4---5
  !                 |   |   |   |   |
  !                 6---7---8---9--10
  !                 |   |   |   |   |
  !                11--12--13--14--15
  !
  !    nedgt:  22
  !    element #:  1  2  3  4  5 ...
  !    pelel:  1  3  8  9 13 ...
  !                |  |  |
  !                |  |  +--------------+
  !                |  |                 |
  !                |  +--+              |
  !                |     |              |
  !    lelel:  2  6  1  3  6  7  8  2  4  7  8  9 ...
  !
  ! USED BY
  !    par_partit
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use mod_memory
  use mod_htable
  use def_parall
  use mod_elmgeo,         only : elmgeo_number_nodes
  use mod_elmgeo,         only : element_type
  use mod_graphs,         only : graphs_number_to_linked_list
  use mod_domain,         only : domain_memory_allocate
  use mod_domain,         only : domain_memory_deallocate
  use mod_memory_config,  only : memory_config
  implicit none  
  integer(ip)           :: ielty,ielem,iface,inodb,ilist,kelem
  integer(ip)           :: inode,jelem,jface,jelty,ipoin
  integer(ip)           :: pepoi,ielpo,dummi,pnodi,pnodb
  integer(ip)           :: melel,pelty,pnodj,jpoin,pnode,nedgt
  integer(8)            :: msize
  logical(lg)           :: equal_faces
  type(hash_t)          :: ht
  integer(ip),  pointer :: faces(:,:,:) 
  integer(ip),  pointer :: lista(:) 
  !
  ! Allocate memory for PELEL
  !
  nullify(faces)
  nullify(lista)
  call domain_memory_allocate('PELEL')
  
  select case ( kfl_parti_par )

  case ( 1_ip )

     !-------------------------------------------------------------------
     !
     ! Graph is based on common nodes
     !
     !-------------------------------------------------------------------

     melel = mepoi*mnode ! Optimize: it is overdimensionalized

     if( memory_config%high_memory ) then
        !
        ! First option: we have a lot of memory
        !
        msize = int(melel,KIND=8)
        msize = msize*int(nelem,KIND=8)
        !
        ! Compute Hash table (initialize, reset, add and destroy)
        !
        call memory_alloca(memor_dom,'LISTA','par_elmgra',lista,msize)
        call htaini( ht, melel )
        pelel(1) = 1
        do ielem= 1, nelem
           call htares( ht, lista(pelel(ielem):) )
           pelty = ltype(ielem)
           pnode = elmgeo_number_nodes(pelty,lnods(:,ielem))
           do inode= 1,pnode
           !do inode= 1,lnnod(ielem)
              jpoin = lnods(inode,ielem)
              do jelem = pelpo(jpoin), pelpo(jpoin+1)-1
                 kelem = lelpo(jelem)
                 if (kelem/=ielem) then
                    call htaadd( ht, kelem )
                 end if
              end do
           end do
           pelel(ielem+1) = pelel(ielem) + ht%nelem
        end do
        call htades( ht )
        nedgt = pelel(nelem+1)-1
        !
        ! Allocate memory and compute list of adjacancies LELEL
        !
        call domain_memory_allocate('LELEL',NUMBER1=nedgt)
        do ielem = 1, nedgt
           lelel(ielem) = lista(ielem)
        end do

     else
        !
        ! Second option: we DO NOT have a lot of memory
        !
        msize = int(melel,KIND=8)
        call memory_alloca(memor_dom,'LISTA','par_elmgra',lista,msize)
        call htaini( ht, melel )
        pelel(1) = 1
        do ielem= 1, nelem
           call htares( ht, lista )
           pelty = ltype(ielem)
           pnode = elmgeo_number_nodes(pelty,lnods(:,ielem))
           do inode= 1,pnode
              jpoin = lnods(inode,ielem)
              do jelem = pelpo(jpoin), pelpo(jpoin+1)-1
                 kelem = lelpo(jelem)
                 if (kelem/=ielem) then
                    call htaadd( ht, kelem )
                 end if
              end do
           end do
           pelel(ielem+1) = pelel(ielem) + ht%nelem
        end do
        nedgt = pelel(nelem+1)-1
        call domain_memory_allocate('LELEL',NUMBER1=nedgt)
        do ielem= 1, nelem
           call htares( ht, lelel(pelel(ielem):) )
           pelty = ltype(ielem)
           pnode = elmgeo_number_nodes(pelty,lnods(:,ielem))
           do inode= 1,pnode
              jpoin = lnods(inode,ielem)
              do jelem = pelpo(jpoin), pelpo(jpoin+1)-1
                 kelem = lelpo(jelem)
                 if (kelem/=ielem) then
                    call htaadd( ht, kelem )
                 end if
              end do
           end do
        end do
        call htades( ht )
     end if
     !
     ! Deallocate LISTA
     !
     if(kfl_parti_par==1) then
        call memory_deallo(memor_dom,'LISTA','par_elmgra',lista)
     end if

  case ( 2_ip )

     !-------------------------------------------------------------------
     !
     ! Graph is based on common faces
     !
     !-------------------------------------------------------------------
     nedgt = 0
     !
     ! Faces graph: Allocate memory for FACES
     !
     call memory_alloca(memor_dom,'FACES','par_elmgra',faces,mnodb+1,mface,nelem)
     !
     ! Construct and sort FACES
     !
     !*OMP  PARALLEL DO SCHEDULE (GUIDED)           & 
     !*OMP  DEFAULT (NONE)                          &
     !*OMP  PRIVATE (ielem,ielty,iface,inodb,inode) &
     !*OMP  SHARED  (ltype,faces,lnods,nelem) 
     !
     do ielem = 1,nelem                                         
        ielty = abs(ltype(ielem))
        do iface = 1,element_type(ielty) % number_faces
           pnodb = element_type(ielty) % node_faces(iface)
           do inodb = 1,pnodb
              inode = element_type(ielty) % list_faces(inodb,iface)
              faces(inodb,iface,ielem) = lnods(inode,ielem)
           end do
           call sortin(pnodb,faces(1,iface,ielem))
        end do
     end do
     !
     ! Compute FACES
     !
     do ielem = 1,nelem                                          ! Compare the faces and 
        ielty = abs(ltype(ielem))                                ! eliminate the repited faces
        do iface=1,element_type(ielty) % number_faces
           ipoin = faces(1,iface,ielem)
           if( ipoin /= 0 ) then
              pnodi = element_type(ielty) % node_faces(iface)
              ilist = 1
              pepoi = pelpo(ipoin+1)-pelpo(ipoin)
              ielpo = pelpo(ipoin)-1
              do while( ilist <= pepoi )
                 ielpo = ielpo+1
                 jelem = lelpo(ielpo)
                 if( jelem /= ielem ) then
                    jelty = abs(ltype(jelem))                              ! eliminate the repited faces
                    jface = 0
                    do while( jface /= element_type(jelty) % number_faces )
                       jface = jface + 1
                       if( faces(1,jface,jelem) /= 0 ) then
                          equal_faces = .true.
                          inodb       = 0
                          pnodj       = element_type(jelty) % node_faces(jface)
                          do while( equal_faces .and. inodb /= pnodj )
                             inodb = inodb + 1
                             if( faces(inodb,iface,ielem) /= faces(inodb,jface,jelem) ) &
                                  equal_faces = .false.
                          end do
                          if( equal_faces ) then
                             faces(1,iface,ielem) = 0                         ! IFACE and JFACE
                             faces(1,jface,jelem) = 0                         ! are eliminated
                             faces(pnodj+1,iface,ielem) = jelem               ! IFACE and JFACE
                             faces(pnodi+1,jface,jelem) = ielem               ! are eliminated
                             nedgt = nedgt + 2
                             jface = element_type(jelty) % number_faces       ! Exit JFACE do
                             ilist = pepoi                                    ! Exit JELEM do
                          end if
                       end if
                    end do
                 end if
                 ilist = ilist + 1
              end do
           end if
        end do
     end do
     !
     ! Allocate memoty for adjacancies LELEL
     !
     call domain_memory_allocate('LELEL',NUMBER1=nedgt)     
     !
     ! Compute PELEL and LELEL
     !
     dummi    = 0
     pelel(1) = 1
     !
     !*OMP  PARALLEL DO SCHEDULE (GUIDED)      & 
     !*OMP  DEFAULT (NONE)                     &
     !*OMP  PRIVATE (dummi,ielem,ielty,iface)  &
     !*OMP  SHARED  (pelel,ltype,faces,lelel,nelem) 
     !
     do ielem = 1,nelem
        ielty = abs(ltype(ielem))
        do iface = 1,element_type(ielty) % number_faces
           if( faces(1,iface,ielem) == 0 ) then
              pnodi          = element_type(ielty) % node_faces(iface)
              dummi          = dummi + 1
              pelel(ielem+1) = pelel(ielem+1) + 1
              lelel(dummi)   = faces(pnodi+1,iface,ielem)
           end if
        end do
        pelel(ielem+1) = pelel(ielem) + pelel(ielem+1)
     end do
     !
     ! Deallocate memory of FACES
     !
     call memory_deallo(memor_dom,'FACES','par_elmgra',faces)

  case(3_ip)

     !-------------------------------------------------------------------
     !
     ! Graph is based on common nodes: brute force
     !
     !-------------------------------------------------------------------

     !call runend('NOT PROGRAMMED')
     melel = mepoi*mnode ! Optimize: it is overdimensionalized
     msize = int(melel,KIND=8)
     call memory_alloca(memor_dom,'LISTA','par_elmgra',lista,msize)

     nedgt = 0
     do ielem = 1,nelem
        ielty = ltype(ielem)
        pnodi = elmgeo_number_nodes(pelty,lnods(:,ielem))
        lista = 0
        do inode = 1,pnodi 
           ipoin = lnods(inode,ielem)
           ielpo = pelpo(ipoin)-1
           do while( ielpo < pelpo(ipoin+1)-1 )
              ielpo = ielpo + 1
              jelem = lelpo(ielpo)
              if( jelem /= ielem ) then
                 ilist = 1
                 loop_lista_1: do while( lista(ilist) /= 0 ) 
                    if( lista(ilist) == jelem ) exit loop_lista_1
                    ilist = ilist + 1
                 end do loop_lista_1
                 if( lista(ilist) == 0 ) then                 
                    lista(ilist) = jelem 
                    nedgt        = nedgt + 1
                    pelel(ielem) = pelel(ielem)+1
                 end if
              end if
           end do
        end do
     end do

     call domain_memory_allocate('LELEL',NUMBER1=nedgt)     
     
     call graphs_number_to_linked_list(nelem,pelel)    

     nedgt = 0
     do ielem = 1,nelem
        ielty = ltype(ielem)
        pnodi = elmgeo_number_nodes(pelty,lnods(:,ielem))
        lista = 0
        do inode = 1,pnodi 
           ipoin = lnods(inode,ielem)
           ielpo = pelpo(ipoin)-1
           do while( ielpo < pelpo(ipoin+1)-1 )
              ielpo = ielpo + 1
              jelem = lelpo(ielpo)
              if( jelem /= ielem ) then
                 ilist = 1
                 loop_lista_2: do while( lista(ilist) /= 0 ) 
                    if( lista(ilist) == jelem ) exit loop_lista_2
                    ilist = ilist + 1
                 end do loop_lista_2
                 if( lista(ilist) == 0 ) then                 
                    lista(ilist) = jelem 
                    nedgt        = nedgt + 1
                    lelel(nedgt) = jelem
                    pelel(ielem) = pelel(ielem)+1
                 end if
              end if
           end do
        end do
     end do

  end select
  !
  ! Maximum number of edges in the mesh
  !
  medge=0
  do ielem=1,nelem
     if(pelel(ielem+1)-pelel(ielem)>medge) then
        medge=pelel(ielem+1)-pelel(ielem)
        dummi=ielem
     end if
  end do
  !
  ! Pointers
  !
  padja_par => pelel  
  ladja_par => lelel

end subroutine par_elmgra
