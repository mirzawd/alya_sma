!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_grave2()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_grave2
  ! NAME
  !   tur_grave2
  ! DESCRIPTION
  !    Compute the 2nd order velocity gradient norm:
  !    (d^2 ui) / ( dx_m dx_n ) needed by Launder-Sharma turbulence model
  !
  !    1. Project du/dx, dv/dy, du/dy and dv/dx 
  !    2. Set f = du/dx + dv/dy + du/dy + dv/dx 
  !    3. Project df/dx and df/dy
  !    4. Set g = df/dx + df/dy
  !        => g =    d/dx ( du/dx + dv/dy + du/dy + dv/dx )
  !                + d/dy ( du/dx + dv/dy + du/dy + dv/dx )
  !             =    d^2u/dx^2 + d^2v/dxdy + d^2u/dxdy + d^2v/dx^2
  !                + d^2u/dxdy + d^2v/dy^2 + d^2u/dy2  + d^2v/dxdy
  !             =    d^2u/dx^2 + d^2v/dx^2 + d^2u/dy2  + d^2v/dy^2
  !             =  + 2 d^2v/dxdy + 2 d^2u/dxdy
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
  integer(ip) :: ipoin,idime,inode,ielem,jdime,jnode,knode
  integer(ip) :: pnode,pelty
  real(rp)    :: detjm,gpvol,cartc(ndime,mnode) 
  real(rp)    :: elgra(ndime,ndime)
  real(rp)    :: elunk(mnode),elcod(ndime,mnode)
  real(rp)    :: elvel(ndime,mnode),elgrs(ndime)
  real(rp)    :: xjaci(9),xjacm(9)

  if(kfl_paral/=0.and.kfl_grve2_tur==1) then

     call memgen(zero,npoin,zero)
     do ipoin=1,npoin
        gesca(ipoin)     = 0.0_rp
        grve2_tur(ipoin) = 0.0_rp
     end do

     do ielem = 1,nelem
        pelty = ltype(ielem) 
        pnode = nnode(pelty)
        if( pelty == PYR05 ) call runend('TUR_GRAVE2: NOT CODED FOR PYR05')
        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
           do idime=1,ndime
              elcod(idime,inode) = coord(idime,ipoin)
              elvel(idime,inode) = advec(idime,ipoin,1)
           end do
        end do
        knode = pnode
        if( lelch(ielem) == ELEXT ) knode = 1
        do inode=1,knode                         ! Loop over Gauss points (which are nodes)
           ipoin = lnods(inode,ielem)
           call elmder(&
                pnode,ndime,elmar(pelty)%deric(1,1,inode),&
                elcod,cartc,detjm,xjacm,xjaci)
           gpvol=elmar(pelty)%weigc(inode)*detjm
           do idime=1,ndime                      ! Velocity gradients
              do jdime=1,ndime
                 elgra(jdime,idime)=0.0_rp
              end do
           end do
           do jnode=1,pnode
              do idime=1,ndime
                 do jdime=1,ndime
                    elgra(idime,jdime) = elgra(idime,jdime)&
                         +cartc(idime,jnode)*elvel(jdime,jnode)
                 end do
              end do
           end do
           !
           ! du/dx + dv/dy + du/dy + dv/dx
           !
           gesca(ipoin)=gesca(ipoin)&
                +gpvol*(elgra(1,1)+elgra(2,2)+elgra(1,2)+elgra(2,1))
           !
           ! + dw/dz + dw/dx + du/dz + dw/dy + dv/dz
           !
           if(ndime==3)&
                gesca(ipoin)=gesca(ipoin)&
                +gpvol*(elgra(3,3)+elgra(1,3)+elgra(3,1)+elgra(2,3)+elgra(3,2))
        end do
     end do

     call rhsmod(1_ip,gesca)                    ! Periodicity
     do ipoin=1,npoin                           ! Solve diagonal system
        gesca(ipoin)=gesca(ipoin)/vmasc(ipoin)
     end do

     do ielem=1,nelem
        pelty=ltype(ielem) 
        pnode=nnode(pelty)
        if( pelty == PYR05 ) call runend('TUR_GRAVE2b: NOT CODED FOR PYR05')
        do inode=1,pnode
           ipoin=lnods(inode,ielem)
           do idime=1,ndime
              elcod(idime,inode)=coord(idime,ipoin)
           end do
           elunk(inode)=gesca(ipoin)
        end do
        knode = pnode
        if( lelch(ielem) == ELEXT ) knode = 1
        do inode=1,knode
           ipoin=lnods(inode,ielem)
           call elmder(&
                pnode,ndime,elmar(pelty)%deric(1,1,inode),&
                elcod,cartc,detjm,xjacm,xjaci)
           gpvol=elmar(pelty)%weigc(inode)*detjm
           if(kfl_naxis==1) then
              if(elcod(1,inode)==0.0_rp) then
                 gpvol=gpvol*twopi*1.0e-12_rp
              else
                 gpvol=gpvol*twopi*elcod(1,inode)
              end if
           end if
           do idime=1,ndime
              elgrs(idime)=0.0_rp
           end do
           do jnode=1,pnode
              do idime=1,ndime
                 elgrs(idime)=elgrs(idime)&
                      +cartc(idime,jnode)*elunk(jnode)
              end do
           end do
           !
           !   d/dx ( du/dx + dv/dy + du/dy + dv/dx )
           ! + d/dy ( du/dx + dv/dy + du/dy + dv/dx )
           !
           grve2_tur(ipoin)=grve2_tur(ipoin)+gpvol*(elgrs(1)+elgrs(2))
           if(ndime==3) grve2_tur(ipoin)=grve2_tur(ipoin)+gpvol*elgrs(3)
        end do
     end do

     call rhsmod(1_ip,grve2_tur)                ! Periodicity
     do ipoin=1,npoin                           ! Solve diagonal system
        grve2_tur(ipoin)=grve2_tur(ipoin)/vmasc(ipoin)
     end do

     call memgen(two,npoin,zero)

  end if

end subroutine tur_grave2
