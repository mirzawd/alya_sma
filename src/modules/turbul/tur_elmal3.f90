!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_elmla3(&
     pnode,pgaus,lnods,elvel,gpcar,gpsha,gpvol,gpvel,&
     gpdif,elmal,elrhs)
  !----------------------------------------------------------------------
  !****f* Turbul/tur_elmla3
  ! NAME 
  !    tur_elmla3
  ! DESCRIPTION
  !    Compute the Laplacian matrix
  ! USES
  ! USED BY
  !    tur_elmope
  !***
  !----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp 
  use def_turbul, only     :  kfl_fixno_tur,ustar_tur
  use def_domain, only     :  mnode,ndime,lpoty
  implicit none
  integer(ip), intent(in)  :: pnode,pgaus
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus),elvel(ndime,pnode)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus),gpvol(pgaus)
  real(rp),    intent(in)  :: gpdif
  real(rp),    intent(out) :: gpvel(ndime,pgaus)
  real(rp),    intent(out) :: elmal(pnode,pnode),elrhs(pnode)
  integer(ip)              :: inode,jnode,idime,igaus,ipoin,ibopo
  real(rp)                 :: fact1
  !
  ! Initialization
  !
  do jnode=1,pnode
     elrhs(jnode)=0.0_rp
     do inode=1,pnode
        elmal(inode,jnode)=0.0_rp
     end do
  end do
  !
  ! Velocity
  !
  do igaus=1,pgaus
     do idime=1,ndime
        gpvel(idime,igaus)=0.0_rp
     end do
  end do
  do igaus=1,pgaus
     do inode=1,pnode
        do idime=1,ndime
           gpvel(idime,igaus)=gpvel(idime,igaus)&
                +elvel(idime,inode)*gpsha(inode,igaus)
        end do
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
           fact1=fact1*gpvol(igaus)*gpdif
           elmal(inode,jnode)=elmal(inode,jnode)+fact1
           elmal(jnode,inode)=elmal(jnode,inode)+fact1
        end do
        fact1=0.0_rp
        do idime=1,ndime
           fact1=fact1+gpcar(idime,inode,igaus)&
                &     *gpcar(idime,inode,igaus)
        end do
        elmal(inode,inode)=elmal(inode,inode)+fact1*gpvol(igaus)*gpdif
        fact1=0.0_rp
        do idime=1,ndime
           fact1=fact1+gpcar(idime,inode,igaus)*gpvel(idime,igaus)
        end do
        fact1=fact1*gpvol(igaus)
        do jnode=1,pnode
           elmal(jnode,inode)=elmal(jnode,inode)&
                +fact1*gpsha(jnode,igaus)
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
           elrhs(inode)=fact1*ustar_tur(ipoin)
        end if
     end if
  end do

end subroutine tur_elmla3
