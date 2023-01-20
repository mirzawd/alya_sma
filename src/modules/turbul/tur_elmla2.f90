!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_elmla2(&
     pnode,pgaus,lnods,gpcar,gpsha,gpvol,elmal,elrhs)
  !----------------------------------------------------------------------
  !****f* Turbul/tur_elmla2
  ! NAME 
  !    tur_elmla2
  ! DESCRIPTION
  !    Compute the Laplacian matrix
  ! USES
  ! USED BY
  !    tur_elmope
  !***
  !----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp 
  use def_turbul, only     :  kfl_fixno_tur,ustar_tur
  use def_domain, only     :  mnode,ndime,lpoty,exnor
  implicit none
  integer(ip), intent(in)  :: pnode,pgaus
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus),gpvol(pgaus)
  real(rp),    intent(out) :: elmal(pnode,pnode),elrhs(ndime,pnode)
  integer(ip)              :: inode,jnode,idime,igaus,ipoin,ibopo
  real(rp)                 :: fact1
  !
  ! Initialization
  !
  do jnode=1,pnode
     do idime=1,ndime
        elrhs(idime,jnode)=0.0_rp
     end do
     do inode=1,pnode
        elmal(inode,jnode)=0.0_rp
     end do
  end do
  !
  ! Laplacian matrix: ( grad p , grad q )
  !
  do igaus=1,pgaus
     do inode=1,pnode
        do jnode=inode+1,pnode
           fact1=0.0_rp
           do idime=1,ndime
              fact1=fact1+gpcar(idime,inode,igaus)&
                   &     *gpcar(idime,jnode,igaus)
           end do
           fact1=fact1*gpvol(igaus)
           elmal(inode,jnode)=elmal(inode,jnode)+fact1
           elmal(jnode,inode)=elmal(jnode,inode)+fact1
        end do
        fact1=0.0_rp
        do idime=1,ndime
           fact1=fact1+gpcar(idime,inode,igaus)&
                &     *gpcar(idime,inode,igaus)
        end do
        fact1=fact1*gpvol(igaus)
        elmal(inode,inode)=elmal(inode,inode)+fact1
        do idime=1,ndime
           elrhs(idime,inode)=elrhs(idime,inode)+gpvol(igaus)*gpsha(inode,igaus)
        end do
     end do
  end do
  !
  ! Prescribe Laplacian on walls
  !
  do inode=1,pnode
     ipoin=lnods(inode)
     if(kfl_fixno_tur(1,ipoin,1)==3.or.kfl_fixno_tur(1,ipoin,1)==4) then
        ibopo=lpoty(ipoin)
        if(ibopo/=0) then
           fact1=elmal(inode,inode)
           do jnode=1,pnode
              elmal(inode,jnode)=0.0_rp
           end do
           elmal(inode,inode)=fact1
           do idime=1,ndime
              elrhs(idime,inode)=-fact1*ustar_tur(ipoin)*exnor(idime,1,ibopo)
           end do
        end if
     end if
  end do

end subroutine tur_elmla2
