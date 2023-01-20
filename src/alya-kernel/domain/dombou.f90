!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine dombou
!-----------------------------------------------------------------------
!****f* Domain/dombou
! NAME
!    dombou
! DESCRIPTION
!    Compute boundary arrays automatically
! OUTPUT
!    NBOUN .................. # of boundaries
!    LNODB(MNODB,NBOUN) ..... Boundary connectivity
!    LBOEL(MNODB,NBOUN) ..... Boundary/Element connectivity
!    LELBO(NBOUN) ........... Boundary element
!    LTYPB(NBOUN) ........... Type of boundary
! USED BY
!    Domain
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_deallo
  use def_elmtyp, only : BOFEM
  use mod_elmgeo, only : element_type
  implicit none 
  integer(ip), allocatable :: faces(:,:,:),bomsh(:,:)
  integer(ip)              :: ielty,ielem,iboun,iface,inodb,ilist
  integer(ip)              :: inode,jelem,jface,jelty,ipoin,jpoin,iblty
  integer(ip)              :: pepoi,ielpo,pnodb,pnodi
  logical(lg)              :: equal_faces

  if( kfl_autbo == 1 .and. INOTSLAVE .and. .not. READ_AND_RUN() ) then
     !
     ! Allocate memory for FACES, BOMSH
     ! 
     allocate(faces(mnodb,mface,nelem)) 
     allocate(bomsh(mface,nelem))
     !
     ! Construct and sort FACES
     !
     do ielem = 1,nelem                                         
        ielty = ltype(ielem)
        do iface = 1,element_type(ielty) % number_faces
           pnodb = element_type(ielty) % node_faces(iface)           
           do inodb = 1,pnodb
              inode = element_type(ielty) % list_faces(inodb,iface)
              faces(inodb,iface,ielem)=lnods(inode,ielem)
           end do
           call sortin(pnodb,faces(1,iface,ielem))
        end do
     end do
     !
     ! Compute NBOUN and fill in BOMSH
     !
     nboun=0
     do ielem=1,nelem                                            ! Compare the faces and 
        ielty=ltype(ielem)                                       ! eliminate the repited faces
        do iface=1,element_type(ielty) % number_faces
           ipoin=faces(1,iface,ielem)
           if(ipoin/=0) then
              pnodi = element_type(ielty) % node_faces(iface)
              ilist=1
              pepoi=pelpo(ipoin+1)-pelpo(ipoin)
              ielpo=pelpo(ipoin)-1
              do while(ilist<=pepoi)
                 ielpo=ielpo+1
                 jelem=lelpo(ielpo)
                 if(jelem/=ielem) then
                    jelty=ltype(jelem)                              ! eliminate the repited faces
                    jface=0
                    do while(jface/=element_type(jelty) % number_faces)
                       jface=jface+1
                       if(faces(1,jface,jelem)/=0) then
                          equal_faces=.true.
                          inodb=0
                          do while(equal_faces .and. inodb /= element_type(jelty) % node_faces(jface) )
                             inodb=inodb+1
                             if(faces(inodb,iface,ielem) /=faces(inodb,jface,jelem))&
                                  equal_faces=.false.
                          end do
                          if(equal_faces) then
                             faces(1,iface,ielem)=0                   ! IFACE and JFACE
                             faces(1,jface,jelem)=0                   ! are eliminated
                             jface=element_type(jelty) % number_faces ! Exit JFACE do
                             ilist=pepoi                              ! Exit JELEM do
                          end if
                       end if
                    end do
                 end if
                 ilist=ilist+1
              end do
              if(faces(1,iface,ielem)/=0) then
                 nboun=nboun+1
                 do inodb=1,element_type(ielty) % node_faces(iface) 
                    inode=element_type(ielty) % list_faces(inodb,iface)
                    faces(inodb,iface,ielem)=lnods(inode,ielem)
                 end do
                 bomsh(iface,ielem)=nboun
              end if
           end if
        end do
     end do
     !
     ! Allocate memory for LNODB, LBOEL and LTYPB
     !
     call memory_deallo(memor_dom,'LNODB'    ,'dombou',lnodb)
     call memory_deallo(memor_dom,'LBOEL'    ,'dombou',lboel)
     call memory_deallo(memor_dom,'LELBO'    ,'dombou',lelbo)
     call memory_deallo(memor_dom,'LTYPB'    ,'dombou',ltypb)
     call memory_deallo(memor_dom,'LBOCH'    ,'dombou',lboch)
     call memory_deallo(memor_dom,'LBINV_LOC','dombou',lbinv_loc)
     
     call memory_alloca(memor_dom,'LNODB'    ,'dombou',lnodb,mnodb,nboun)
     call memory_alloca(memor_dom,'LBOEL'    ,'dombou',lboel,mnodb,nboun)
     call memory_alloca(memor_dom,'LELBO'    ,'dombou',lelbo,nboun)
     call memory_alloca(memor_dom,'LTYPB'    ,'dombou',ltypb,nboun)
     call memory_alloca(memor_dom,'LBOCH'    ,'dombou',lboch,nboun)
     call memory_alloca(memor_dom,'LBINV_LOC','dombou',lbinv_loc,nboun,'IDENTITY')
     !
     ! Compute LNODB, LBOEL AND LTYPB
     !
     kfl_bouel=1
     do ielem=1,nelem
        ielty=ltype(ielem)
        do iface=1,element_type(ielty) % number_faces
           if(faces(1,iface,ielem)/=0) then                
              call fintyp(ndimb,element_type(ielty) % node_faces(iface),iblty)
              iboun=bomsh(iface,ielem)
              lelbo(iboun)=ielem
              lelbo(iboun)=ielem              
              ltypb(iboun)=iblty
              lboch(iboun)=BOFEM
              do inodb=1,nnode(iblty) 
                 lnodb(inodb,iboun)=faces(inodb,iface,ielem)
              end do
              do inodb=1,nnode(iblty)
                 ipoin=lnodb(inodb,iboun)
                 nodes: do inode=1,nnode(ielty)
                    jpoin=lnods(inode,ielem)
                    if(ipoin==jpoin) then
                       lboel(inodb,iboun)=inode
                       exit nodes
                    end if
                 end do nodes
              end do

           end if
        end do
     end do

!!   ONLY FOR DEBUGGING AND FOR SEQUENTIAL CASES
!!     do iboun= 1,nboun
!! alf        write(758,*) iboun,coord(1:3,lnodb(1,iboun)),coord(1:3,lnodb(2,iboun)),coord(1:3,lnodb(3,iboun))
!! jaz        write(758,1000) iboun,lnodb(1:4,iboun),lboel(5,iboun)
!!     end do
!!     stop

     !
     ! Deallocate memory of FACES, BOMSH
     !
     deallocate(faces)
     deallocate(bomsh)

  end if


!! jaz 1000 format (6(i2,i8))

end subroutine dombou

