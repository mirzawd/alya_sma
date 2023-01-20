!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_elmlap(&
     itask,pnode,pgaus,lnods,gpcar,gpsha,gpvol,elmal,elrhs)
  !----------------------------------------------------------------------
  !****f* Levels/lev_elmlap
  ! NAME 
  !    tur_elmlap
  ! DESCRIPTION
  !    Compute the Laplacian matrix  - nearly identical to tur_elmlap
  !    Itask == 1 the original version for system: Lapl(f)=-1, with f=0 on wall
  !    Itask == 2 use to extend the wall distance from the wall using system: Lapl(g)=0, with g = ywalp on wall
  ! USES
  ! USED BY
  !    lev_distuc
  !***
  !----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp 
  use def_levels, only     :  icupt_lev,flev0_lev
  use def_domain, only     :  mnode,ndime
  implicit none
  integer(ip), intent(in)  :: itask,pnode,pgaus
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus),gpvol(pgaus)
  real(rp),    intent(out) :: elmal(pnode,pnode),elrhs(pnode)
  integer(ip)              :: inode,jnode,kdime,igaus,ipoin
  real(rp)                 :: fact1,bvess
  !
  ! Initialization
  !
  do jnode=1,pnode
     elrhs(jnode)=0.0_rp
     do inode=1,pnode
        elmal(inode,jnode)=0.0_rp
     end do
  end do

  if (itask==1) then
     !
     ! Laplacian matrix: ( grad p , grad q ), and rhs
     !
     do igaus=1,pgaus
        do inode=1,pnode
           do jnode=inode+1,pnode
              fact1=0.0_rp
              do kdime=1,ndime
                 fact1=fact1+gpcar(kdime,inode,igaus)&
                      &     *gpcar(kdime,jnode,igaus)
              end do
              fact1=fact1*gpvol(igaus)
              elmal(inode,jnode)=elmal(inode,jnode)+fact1
              elmal(jnode,inode)=elmal(jnode,inode)+fact1
           end do
           fact1=0.0_rp
           do kdime=1,ndime
              fact1=fact1+gpcar(kdime,inode,igaus)&
                   &     *gpcar(kdime,inode,igaus)
           end do
           fact1=fact1*gpvol(igaus)
           elmal(inode,inode)=elmal(inode,inode)+fact1
           elrhs(inode)=elrhs(inode)+gpvol(igaus)*gpsha(inode,igaus)
        end do
     end do
     !
     ! Prescribe Laplacian on walls
     !
     do inode=1,pnode
        ipoin=lnods(inode)
        if(icupt_lev(ipoin) > 0_ip) then
           fact1=elmal(inode,inode)
           do jnode=1,pnode
              elmal(jnode,inode)=0.0_rp              
              elmal(inode,jnode)=0.0_rp              
           end do
           elmal(inode,inode)=fact1
           elrhs(inode)=0.0_rp
        end if
     end do
  else
     !
     ! Laplacian matrix: ( grad p , grad q ), and rhs
     !
     do igaus=1,pgaus
        do inode=1,pnode
           do jnode=inode+1,pnode
              fact1=0.0_rp
              do kdime=1,ndime
                 fact1=fact1+gpcar(kdime,inode,igaus)&
                      &     *gpcar(kdime,jnode,igaus)
              end do
              fact1=fact1*gpvol(igaus)
              elmal(inode,jnode)=elmal(inode,jnode)+fact1
              elmal(jnode,inode)=elmal(jnode,inode)+fact1
           end do
           fact1=0.0_rp
           do kdime=1,ndime
              fact1=fact1+gpcar(kdime,inode,igaus)&
                   &     *gpcar(kdime,inode,igaus)
           end do
           fact1=fact1*gpvol(igaus)
           elmal(inode,inode)=elmal(inode,inode)+fact1
        end do
     end do
     !
     ! Prescribe Laplacian on walls
     !
     do inode=1,pnode
        ipoin=lnods(inode)
        if(icupt_lev(ipoin) > 0_ip) then
           bvess=abs(flev0_lev(ipoin))
           fact1=elmal(inode,inode)
           do jnode=1,pnode
              elmal(inode,jnode)=0.0_rp
              elrhs(jnode)=elrhs(jnode)&
                   -elmal(jnode,inode)*bvess
              elmal(jnode,inode)=0.0_rp
           end do
           elmal(inode,inode)=fact1
           elrhs(inode)=fact1*bvess
        end if
     end do
  end if

end subroutine lev_elmlap
