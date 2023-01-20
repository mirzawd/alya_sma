!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine exm_gravec(veinp,grunk)

  !------------------------------------------------------------------------
  !
  ! Case of a vector
  !
  !------------------------------------------------------------------------
  use def_kintyp
  use def_elmtyp   ! type of elements
  use def_domain   ! geometry information
  use def_master   ! geometry information
  implicit none
  real(rp),   intent(in)          :: veinp(ndime,npoin)
  real(rp),   intent(out), target :: grunk(ndime,ndime,npoin)
  integer(ip)                     :: ipoin,idime,inode,ielem,jdime
  integer(ip)                     :: pnode,pelty
  real(rp)                        :: detjm,gpvol,cartc(ndime,mnode),xauxi
  real(rp)                        :: elgra(ndime,ndime)
  real(rp)                        :: elunk(ndime,mnode),elcod(ndime,mnode)
  real(rp)                        :: xjaci(9),xjacm(9)

  !
  ! Initialization
  !
  do ipoin=1,npoin
     do idime=1,ndime
        do jdime=1,ndime
           grunk(jdime,idime,ipoin)=0.0_rp
        end do
     end do
  end do
  !
  ! Loop over elements
  !

  elements: do ielem = 1,nelem

     pelty=ltype(ielem) 
     pnode=nnode(pelty)
     !
     ! Gather vectors
     !
     do inode=1,pnode
        ipoin=lnods(inode,ielem)
        elcod(1:ndime,inode)=coord(1:ndime,ipoin)
        elunk(1:ndime,inode)=veinp(1:ndime,ipoin)
     end do
     !     if() then   !si no queremos detectar gradiente con fibras a 180
     do inode= 2,pnode
        xauxi = elunk(1,1)*elunk(1,inode) + elunk(2,1)*elunk(2,inode)
        if (ndime == 3) xauxi = xauxi + elunk(3,1)*elunk(3,inode)
        if (xauxi < 0.0_rp) then
           elunk(1:ndime,inode) = -1.0_rp * elunk(1:ndime,inode)
        end if
     end do
     !  else
     !
     ! Loop over Gauss points (which are nodes)
     !
     gauss_points: do inode=1,pnode
        ipoin=lnods(inode,ielem)
        call elmder(&
             pnode,ndime,elmar(pelty)%deric(1,1,inode),&
             elcod,cartc,detjm,xjacm,xjaci)
        gpvol=elmar(pelty)%weigc(inode)*detjm
        !
        ! Vector derivatives
        !
        do idime=1,ndime
           do jdime=1,ndime                 
              elgra(idime,jdime)=&
                   dot_product(cartc(idime,1:pnode),elunk(jdime,1:pnode))
           end do
        end do
        do idime=1,ndime
           do jdime=1,ndime
              grunk(jdime,idime,ipoin)=grunk(jdime,idime,ipoin)&
                   +gpvol*elgra(jdime,idime)
           end do
        end do
     end do gauss_points
  end do elements
  !
  ! Periodicity
  !
  call rhsmod(ndime*ndime,grunk) 
  !
  ! Parall service
  !
  if(kfl_paral>=1) then
     call runend('GRAVEC: NOT PROGRAMMED')
     !parr3 => grunk
     !call vocabu(NPOIN_REAL_3DIM,ntens,0_ip)
     !call par_slexch()
  end if
  !
  ! Solve diagonal system
  !
  do ipoin=1,npoin
     do idime=1,ndime
        do jdime=1,ndime
           grunk(jdime,idime,ipoin)=grunk(jdime,idime,ipoin)/vmass(ipoin)
        end do
     end do
  end do


end subroutine exm_gravec
