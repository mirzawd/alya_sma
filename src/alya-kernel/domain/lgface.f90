!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lgface()
  !------------------------------------------------------------------------
  !****f* domain/lgface
  ! NAME
  !    lgface
  ! DESCRIPTION
  !    This routine computes the list of global faces:
  !
  !    LFACG(4,NFACG) and LELFA(1:NELEM) % L(1:NFACE(LTYPE(IELEM)))
  !
  !    LFACG(1,IFACG) = IELEM
  !    LFACG(2,IFACG) = JELEM/0 ...... 0 if boundary face
  !    LFACG(3,IFACG) = IFACE ........ Local IELEM face
  !    LFACG(4,IFACG) = JFACE/0/-1 ... Local JELEM face, 0 is a bounday face, 
  !                                    -1 is an extension element face
  !    Caution: 0 < IELEM < NELEM_2 so IELEM can be a neighbor's element
  !    LELFA(IELEM) % L(1:NFACE) ... Global face connectivity for IELEM
  !
  ! USED BY
  !    submsh
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_elmtyp
  use def_kermod
  use def_master
  use def_domain
  use mod_memory
  use mod_htable
  use mod_elmgeo, only : element_type
  implicit none  
  integer(ip) :: ielty,ielem,iface,inodf,ifacg
  integer(ip) :: inode,jelem,jface,jelty,ipoin,pnodf
  integer(ip) :: ielpo,knode
  logical(lg) :: equal_faces  

  !----------------------------------------------------------------------
  !
  ! LELFA: List of global faces
  !
  !----------------------------------------------------------------------

  if( INOTMASTER .and. kfl_lface == 1 ) then
     !
     ! Allocate memory for lelfa
     !
     kfl_lface = 2
     call memory_alloca(memor_dom,'LELFA','lgface',lelfa,nelem_2)
     !
     ! Faces graph: Allocate memory for FACES, CFAEL AND NNODG
     !
     call memory_alloca(memor_dom,'FACEL','lgface',facel,mnodb+1_ip,mface,nelem_2)
     !
     ! Construct and sort FACES
     !
     !*OMP  PARALLEL DO SCHEDULE (GUIDED)           & 
     !*OMP  DEFAULT (NONE)                          & 
     !*OMP  PRIVATE (ielem,ielty,iface,inodf,inode) &
     !*OMP  SHARED  (ltype,cfael,faces,lnods,nnodg,nelem,nface) 
     !
     do ielem = 1,nelem_2                                         
        ielty = abs(ltype(ielem))
        call memory_alloca(memor_dom,'LELFA % L','lgface',lelfa(ielem)%l,element_type(ielty) % number_faces)
        do iface = 1,element_type(ielty) % number_faces
           pnodf = element_type(ielty) % node_faces(iface)
           do inodf = 1,pnodf 
              inode = element_type(ielty) % list_faces(inodf,iface) 
              facel(inodf,iface,ielem) = lnods(inode,ielem)
           end do
           facel(mnodb+1,iface,ielem) = 1
           call sortin(pnodf,facel(1,iface,ielem))
        end do
     end do
     !
     ! Compute FACES
     !
     do ielem = 1,nelem_2                                          ! Compare the faces and 
        ielty = abs(ltype(ielem))                                  ! eliminate the repited faces
        do iface = 1,element_type(ielty) % number_faces
           if( facel(mnodb+1,iface,ielem) > 0 ) then
              nfacg = nfacg + 1
              ipoin = facel(1,iface,ielem)
              ielpo = pelpo_2(ipoin)-1
              do while( ielpo < pelpo_2(ipoin+1)-1 )
                 ielpo = ielpo + 1
                 jelem = lelpo_2(ielpo)
                 if( jelem /= ielem ) then
                    jelty = abs(ltype(jelem))                      ! eliminate the repited faces
                    jface = 0
                    do while( jface < element_type(jelty) % number_faces )
                       jface = jface + 1
                       if( facel(mnodb+1,jface,jelem) > 0 ) then
                          equal_faces = .true.
                          inodf = 0
                          do while( equal_faces .and. inodf < element_type(jelty) % node_faces(jface) )
                             inodf = inodf + 1 
                             if( facel(inodf,iface,ielem) /= facel(inodf,jface,jelem) ) equal_faces = .false.
                          end do
                          if( equal_faces ) then
                             facel(mnodb+1,iface,ielem) =  jelem                              ! Keep IELEM face
                             facel(mnodb+1,jface,jelem) = -ielem                              ! Elminate JELEM face
                             facel(      1,iface,ielem) = -jface                              ! Remember IFACE face
                             jface                      =  element_type(jelty) % number_faces ! Exit JFACE do
                             ielpo                      =  pelpo_2(ipoin+1)                   ! Exit JELEM do  
                          end if
                       end if
                    end do
                 end if
              end do
           end if
        end do
     end do
     !
     ! Allocate memory
     !
     call memory_alloca(memor_dom,'LFACG','lgface',lfacg,4_ip,nfacg)

     nfacg = 0
     do ielem = 1,nelem                                            ! Compare the faces and 
        ielty = abs(ltype(ielem))                                  ! eliminate the repeated faces
        do iface = 1,element_type(ielty) % number_faces
           if( facel(mnodb+1,iface,ielem) > 0 ) then
              nfacg = nfacg + 1
              pnodf = element_type(ielty) % node_faces(iface)
              if( facel(1,iface,ielem) < 0 ) then
                 jelem                   =  facel(mnodb+1,iface,ielem)  ! 
                 jface                   = -facel(      1,iface,ielem)
                 lelfa(ielem) % l(iface) =  nfacg
                 lelfa(jelem) % l(jface) =  nfacg
                 lfacg(1,nfacg)          =  ielem
                 lfacg(2,nfacg)          =  jelem
                 lfacg(3,nfacg)          =  iface
                 lfacg(4,nfacg)          =  jface
              else
                 knode = 0
                 do inodf = 1,pnodf
                    inode = element_type(ielty) % list_faces(inodf,iface)
                    ipoin = lnods(inode,ielem)
                    if( lpoty(ipoin) /= 0 ) knode = knode + 1
                 end do
                 lelfa(ielem) % l(iface) =  nfacg
                 lfacg(1,nfacg)          =  ielem
                 lfacg(2,nfacg)          =  0
                 lfacg(3,nfacg)          =  iface
                 lfacg(4,nfacg)          =  0                
                 !
                 ! This happens with some faces of extension elements
                 ! They are not shared with any other element but are not necessarily 
                 ! boundary faces
                 !
                 if( knode /= pnodf ) lfacg(2,nfacg) = -1 
              end if
           end if
        end do
     end do
     !
     ! Deallocate memory
     !
     call memory_deallo(memor_dom,'FACEL','lgface',facel)
     !
     ! Construct Parall face communication arrays
     !
     call par_lgface(1_ip)  

  end if

  !----------------------------------------------------------------------
  !
  ! LELBF: List of element boundary faces
  !
  !----------------------------------------------------------------------

  if( INOTMASTER .and. kfl_lelbf == 1 ) then

     kfl_lelbf = 2
     call memory_alloca(memor_dom,'LELBF','lgface',lelbf,nelem_2)

     do ielem = 1,nelem_2
        lelbf(ielem) % n = 0
     end do

     do ifacg = 1,nfacg
        if( lfacg(2,ifacg) == 0 ) then
           ielem            = lfacg(1,ifacg)
           lelbf(ielem) % n = lelbf(ielem) % n + 1
        end if
     end do

     do ielem = 1,nelem_2
        if( lelbf(ielem) % n /= 0 ) then
           call memory_alloca(memor_dom,'LELBF % L','lgface',lelbf(ielem) % l,lelbf(ielem) % n)
           lelbf(ielem) % n = 0
        end if
     end do

     do ifacg = 1,nfacg
        if( lfacg(2,ifacg) == 0 ) then
           ielem            = lfacg(1,ifacg)
           lelbf(ielem) % n = lelbf(ielem) % n + 1
           lelbf(ielem) % l(lelbf(ielem) % n) = ifacg
        end if
     end do
 
  end if

end subroutine lgface
