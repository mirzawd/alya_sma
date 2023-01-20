!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chkmas()
  !-----------------------------------------------------------------------
  !****f* domain/chkmas
  ! NAME
  !    chkmas
  ! DESCRIPTION
  !    Check projections
  ! USED BY
  !    Turnon
  !*** 
  !-----------------------------------------------------------------------
  use def_domain 
  use def_elmtyp
  use def_kermod, only : kfl_grpro
  use def_master
  use mod_gradie
  use mod_messages, only : messages_live
  implicit none
  integer(ip)         :: inode,igaus,ielem,ipoin,idime,jdime
  integer(ip)         :: kdime,pdime,pelty,pnode,pgaus,itens
  integer(ip)         :: knode
  real(rp)            :: xfact,gpvol,gpdet,xjaci(9),xjacm(9)
  real(rp)            :: xresu(6,6)
  real(rp)            :: elcod(ndime,mnode)
  real(rp)            :: gpcar(ndime,max(mgaus,mnode))
  real(rp),   pointer :: grunk(:,:,:) => null()
  real(rp),   pointer :: grun2(:,:)   => null()
  real(rp),   pointer :: diffu(:)     => null()

  if( IMASTER ) call runend('END PROJECTION TEST')

  call messages_live('CHECK PROJECTION')
 
  do kfl_grpro = 0,1
     write(*,*) '-----------------------------------------'
     if( kfl_grpro == 0 ) write(*,*) 'CHECK OPEN RULE'
     if( kfl_grpro == 1 ) write(*,*) 'CHECK CLOSE RULE'

     inode = 0
     do idime = 1,6
        do jdime = 1,6
           inode = inode + 1
           xresu(jdime,idime) = real(inode,rp) + 0.5_rp
        end do
     end do

     !----------------------------------------------------------------------
     !
     ! Check projection of constant = constant using open rule
     !
     !----------------------------------------------------------------------

     call memgen(0_ip,npoin,0_ip)
     do ielem = 1,nelem
        pelty = ltype(ielem)
        if( pelty > 0 ) then
           pgaus = ngaus(pelty)
           pnode = nnode(pelty)
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              do idime = 1,ndime
                 elcod(idime,inode) = coord(idime,ipoin)
              end do
           end do
           if( lelch(ielem) == ELEXT ) then
              knode = 1
           else
              knode = pnode
           end if
           do igaus = 1,pgaus
              call elmder(&
                   pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                   elcod,gpcar,gpdet,xjacm,xjaci)
              gpvol = elmar(pelty)%weigp(igaus)*gpdet
              do inode = 1,knode
                 ipoin = lnods(inode,ielem)
                 xfact = gpvol * elmar(pelty) % shape(inode,igaus)
                 gesca(ipoin) = gesca(ipoin) + xfact * 2.5_rp
              end do
           end do
        end if
     end do
     call rhsmod(1_ip,gesca)
     do ipoin = 1,npoin
        xfact = 1.0_rp / vmass(ipoin)
        gesca(ipoin) = xfact * gesca(ipoin)
     end do
     do ipoin = 1,npoin
        if( lnoch(ipoin) /= NOHOL ) then
           if( abs(gesca(ipoin) - 2.5_rp) > 1.0e-12_rp ) then
              write(*,*) gesca(ipoin) ,2.5_rp
              call runend('TEST 1: WRONG PROJECTION')
           end if
        end if
     end do
     call memgen(2_ip,npoin,0_ip)
     write(*,*) 'TEST 1 GOOD'

     !----------------------------------------------------------------------
     !
     ! Check projection of constant = constant using closed rule
     !
     !----------------------------------------------------------------------

     call memgen(0_ip,npoin,0_ip)
     do ielem = 1,nelem
        pelty = ltype(ielem)
        if( pelty > 0 ) then
           pgaus = ngaus(pelty)
           pnode = nnode(pelty)
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              do idime = 1,ndime
                 elcod(idime,inode) = coord(idime,ipoin)
              end do
           end do
           if( lelch(ielem) == ELEXT ) then
              knode = 1
           else
              knode = pnode
           end if
           if( pelty == PYR05 ) then
              do igaus = 1,pgaus
                 call elmder(&
                      pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                      elcod,gpcar,gpdet,xjacm,xjaci)
                 gpvol = elmar(pelty)%weigp(igaus)*gpdet
                 do inode = 1,knode
                    ipoin = lnods(inode,ielem)
                    xfact = gpvol * elmar(pelty) % shape(inode,igaus)
                    gesca(ipoin) = gesca(ipoin) + xfact * 2.5_rp
                 end do
              end do
           else
              do inode = 1,knode
                 call elmder(&
                      pnode,ndime,elmar(pelty)%deric(1,1,inode),&
                      elcod,gpcar,gpdet,xjacm,xjaci)
                 gpvol = elmar(pelty)%weigc(inode)*gpdet
                 ipoin = lnods(inode,ielem)
                 gesca(ipoin) = gesca(ipoin) + gpvol * 2.5_rp
              end do
           end if
        end if
     end do
     call rhsmod(1_ip,gesca)
     do ipoin = 1,npoin
        xfact = 1.0_rp / vmasc(ipoin)
        gesca(ipoin) = xfact * gesca(ipoin)
     end do
     do ipoin = 1,npoin
        if( lnoch(ipoin) /= NOHOL ) then
           if( abs(gesca(ipoin) - 2.5_rp) > 1.0e-12_rp ) then
              write(*,*) gesca(ipoin),2.5_rp
              call runend('TEST 2: WRONG PROJECTION')
           end if
        end if
     end do
     call memgen(2_ip,npoin,0_ip)
     write(*,*) 'TEST 2 GOOD'

     !----------------------------------------------------------------------
     !
     ! GRASCA: Check grad(ui)=xresu(i,1)
     !
     !----------------------------------------------------------------------

     call memgen(0_ip,npoin,0_ip)
     call memgen(0_ip,ndime,npoin)
     do ipoin = 1,npoin
        do idime = 1,ndime
           gesca(ipoin) = gesca(ipoin) + xresu(idime,1) * coord(idime,ipoin) 
        end do
     end do
     call gradie(gesca,gevec)
     do ipoin = 1,npoin
        if( lnoch(ipoin) /= NOHOL ) then
           do idime = 1,ndime
              if(   abs( gevec(idime,ipoin) - xresu(idime,1) ) > 1.0e-12_rp ) then
                 do jdime = 1,ndime
                    write(*,*) gevec(jdime,ipoin),xresu(jdime,1)
                 end do
                 call runend('TEST 3: WRONG PROJECTION')        
              end if
           end do
        end if
     end do
     call memgen(2_ip,npoin,0_ip)
     call memgen(2_ip,ndime,npoin)
     write(*,*) 'TEST 3 GOOD'

     !----------------------------------------------------------------------
     !
     ! GRAVEC: Check grad(uij)
     !
     !----------------------------------------------------------------------

     allocate(grunk(ndime,ndime,npoin))
     call memgen(0_ip,ndime,npoin)
     do ipoin = 1,npoin
        do idime = 1,ndime
           do jdime = 1,ndime
              gevec(jdime,ipoin) = gevec(jdime,ipoin) + xresu(idime,jdime) * coord(idime,ipoin)
           end do
        end do
     end do
     call gradie(gevec,grunk)
     do ipoin = 1,npoin
        if( lnoch(ipoin) /= NOHOL ) then
           do idime = 1,ndime
              do jdime = 1,ndime
                 if(   abs( grunk(idime,jdime,ipoin) - xresu(idime,jdime) ) > 1.0e-12_rp ) then
                    do kdime = 1,ndime
                       do pdime = 1,ndime
                          write(*,*) grunk(kdime,pdime,ipoin),xresu(kdime,pdime)
                       end do
                    end do
                    call runend('TEST 4: WRONG PROJECTION')        
                 end if
              end do
           end do
        end if
     end do
     call memgen(2_ip,ndime,npoin)
     deallocate(grunk)
     write(*,*) 'TEST 4 GOOD'

     !----------------------------------------------------------------------
     !
     ! GRATEN: Check grad(uij)+grad^t(uij)
     !
     !----------------------------------------------------------------------

     allocate(grun2(ntens,npoin))
     call memgen(0_ip,ndime,npoin)
     do ipoin = 1,npoin
        do idime = 1,ndime
           do jdime = 1,ndime
              gevec(jdime,ipoin) = gevec(jdime,ipoin) + xresu(idime,jdime) * coord(idime,ipoin)
           end do
        end do
     end do
     call gradie(gevec,grun2)
     do ipoin = 1,npoin       
        if( lnoch(ipoin) /= NOHOL ) then
           do itens = 1,ntens
              if( itens == 1 ) then
                 idime = 1
                 jdime = 1
              else if( itens == 2 ) then
                 idime = 2
                 jdime = 2
              else if( itens == 3 ) then
                 idime = 1
                 jdime = 2
              else if( itens == 4 ) then
                 idime = 3
                 jdime = 3
              else if( itens == 5 ) then
                 idime = 1
                 jdime = 3
              else if( itens == 6 ) then
                 idime = 2
                 jdime = 3
              end if
              if( abs( grun2(itens,ipoin) - (xresu(idime,jdime)+xresu(jdime,idime)) ) > 1.0e-12_rp ) then
                 write(*,*) itens,grun2(itens,ipoin),(xresu(idime,jdime)+xresu(jdime,idime))
                 call runend('TEST 5: WRONG PROJECTION')        
              end if
           end do
        end if
     end do
     call memgen(2_ip,ndime,npoin)
     deallocate(grun2)
     write(*,*) 'TEST 5 GOOD'

     !----------------------------------------------------------------------
     !
     ! GRASCA: Check k*grad(ui)=k*xresu(i,1)
     !
     !----------------------------------------------------------------------

     call memgen(0_ip,npoin,0_ip)
     call memgen(0_ip,ndime,npoin)
     allocate(diffu(npoin))
     do ipoin = 1,npoin
        diffu(ipoin) = 9.0_rp
        do idime = 1,ndime
           gesca(ipoin) = gesca(ipoin) + xresu(idime,1) * coord(idime,ipoin)
        end do
     end do
     call gradie(gesca,gevec,diffu)
     do ipoin = 1,npoin
        if( lnoch(ipoin) /= NOHOL ) then
           do idime = 1,ndime
              if(   abs( gevec(idime,ipoin) - diffu(ipoin)*xresu(idime,1) ) > 1.0e-12_rp ) then
                 do jdime = 1,ndime
                    write(*,*) gevec(jdime,ipoin),diffu(ipoin)*xresu(jdime,1)
                 end do
                 call runend('TEST 6: WRONG PROJECTION')        
              end if
           end do
        end if
     end do
     call memgen(2_ip,npoin,0_ip)
     call memgen(2_ip,ndime,npoin)
     deallocate(diffu)
     write(*,*) 'TEST 6 GOOD'

  end do

  call runend('END PROJECTION TEST')

end subroutine chkmas
