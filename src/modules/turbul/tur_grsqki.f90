!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_grsqki()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_grsqki
  ! NAME
  !   tur_grsqki
  ! DESCRIPTION
  !    Compute grad(sqrt(k))^2
  ! USES
  !    memgen
  ! USED BY
  !    tur_begite
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_elmtyp
  use def_master
  use def_domain
  use def_turbul
  implicit none
  integer(ip) :: ipoin,idime,inode,ielem,jnode,pnode,pelty
  real(rp)    :: detjm,gpvol,cartc(ndime,mnode) 
  real(rp)    :: elgrk(ndime),elkin(mnode),elcod(ndime,mnode)
  real(rp)    :: xjaci(9),xjacm(9)

  if( INOTMASTER .and. kfl_grsqk_tur==1 ) then

     call memgen(zero,ndime,npoin)

     if(ndime==2) then
        do ielem = 1,nelem
           pelty = ltype(ielem) 
           if( pelty == PYR05 ) call runend('TUR_GRSQKI: NOT CODED FOR PYR05')
           pnode = nnode(pelty)
           do inode = 1,pnode
              ipoin          = lnods(inode,ielem)
              elkin(inode)   = sqrt(untur(1,ipoin,1))
              elcod(1,inode) = coord(1,ipoin)
              elcod(2,inode) = coord(2,ipoin)
           end do
           do inode=1,pnode                         ! Loop over Gauss points (which are nodes)
              ipoin=lnods(inode,ielem)
              call elmder(&
                   pnode,ndime,elmar(pelty)%deric(1,1,inode),&
                   elcod,cartc,detjm,xjacm,xjaci)
              gpvol    = elmar(pelty)%weigc(inode)*detjm
              elgrk(1) = 0.0_rp
              elgrk(2) = 0.0_rp
              do jnode=1,pnode
                 elgrk(1) = elgrk(1) + cartc(1,jnode)*elkin(jnode)
                 elgrk(2) = elgrk(2) + cartc(2,jnode)*elkin(jnode)
              end do
              gevec(1,ipoin) = gevec(1,ipoin) + gpvol*elgrk(1)
              gevec(2,ipoin) = gevec(2,ipoin) + gpvol*elgrk(2)
           end do
        end do
     else
        do ielem = 1,nelem
           pelty = ltype(ielem) 
           if( pelty == PYR05 ) call runend('TUR_GRSQKIb: NOT CODED FOR PYR05')
           pnode = nnode(pelty)
           do inode=1,pnode
              ipoin          = lnods(inode,ielem)
              elkin(inode)   = sqrt(untur(1,ipoin,1))
              elcod(1,inode) = coord(1,ipoin)
              elcod(2,inode) = coord(2,ipoin)
              elcod(3,inode) = coord(3,ipoin)
           end do
           do inode=1,pnode                         ! Loop over Gauss points (which are nodes)
              ipoin=lnods(inode,ielem)
              call elmder(&
                   pnode,ndime,elmar(pelty)%deric(1,1,inode),&
                   elcod,cartc,detjm,xjacm,xjaci)
              gpvol    = elmar(pelty)%weigc(inode)*detjm
              elgrk(1) = 0.0_rp
              elgrk(2) = 0.0_rp
              elgrk(3) = 0.0_rp
              do jnode=1,pnode
                 elgrk(1) = elgrk(1) + cartc(1,jnode)*elkin(jnode)
                 elgrk(2) = elgrk(2) + cartc(2,jnode)*elkin(jnode)
                 elgrk(3) = elgrk(3) + cartc(3,jnode)*elkin(jnode)
              end do
              gevec(1,ipoin) = gevec(1,ipoin) + gpvol*elgrk(1)
              gevec(2,ipoin) = gevec(2,ipoin) + gpvol*elgrk(2)
              gevec(3,ipoin) = gevec(3,ipoin) + gpvol*elgrk(3)
           end do
        end do
     end if
     !
     ! Parall service and Periodicity
     !
     call rhsmod(ndime,gevec)
     !
     ! Solve diagonal system
     !
     do ipoin=1,npoin   
        do idime=1,ndime
           gevec(idime,ipoin)=gevec(idime,ipoin)/vmasc(ipoin)
        end do
     end do
     !
     ! GRSQK_TUR= grad[sqrt(k)].grad[sqrt(k)]
     !
     do ipoin=1,npoin
        grsqk_tur(ipoin)=0.0_rp 
        do idime=1,ndime
           grsqk_tur(ipoin)=grsqk_tur(ipoin)&
                +gevec(idime,ipoin)*gevec(idime,ipoin)
        end do
     end do

     call memgen(two,ndime,npoin)

  end if

end subroutine tur_grsqki
