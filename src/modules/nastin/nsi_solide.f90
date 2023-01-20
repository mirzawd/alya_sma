!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_solide(&
     pnode,elauu,elaup,elapp,elapu,elrbu,elrbp,elmap)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_solide
  ! NAME 
  !    nsi_solide
  ! DESCRIPTION
  !    Compute element matrix and RHS
  ! USES
  ! USED BY
  !    nsi_elmop3
  !***
  !----------------------------------------------------------------------- 
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
  implicit none
  integer(ip), intent(in)    :: pnode
  real(rp),    intent(out)   :: elauu(pnode*ndime,pnode*ndime)
  real(rp),    intent(out)   :: elaup(pnode*ndime,pnode)
  real(rp),    intent(out)   :: elapp(pnode,pnode)
  real(rp),    intent(out)   :: elapu(pnode,pnode*ndime)
  real(rp),    intent(out)   :: elrbu(ndime*pnode)
  real(rp),    intent(out)   :: elrbp(pnode)
  real(rp),    intent(out)   :: elmap(pnode,pnode)
  real(rp)                   :: adiag,adia2
  integer(ip)                :: inode,jnode

  do inode = 1,pnode*ndime
     adiag = 1.0_rp
     do jnode = 1,pnode*ndime
        elauu(jnode,inode) = 0.0_rp
     end do
     do jnode = 1,pnode
        elaup(inode,jnode) = 0.0_rp
        elapu(jnode,inode) = 0.0_rp
     end do
     elrbu(inode)       = 0.0_rp
     elauu(inode,inode) = adiag
  end do

  do inode = 1,pnode
     adiag = 1.0_rp
     adia2 = 1.0_rp
     do jnode = 1,pnode
        elapp(inode,jnode) = 0.0_rp
        elmap(inode,jnode) = 0.0_rp
     end do
     elapp(inode,inode) = adiag
     elmap(inode,inode) = adia2
     elrbp(inode)       = 0.0_rp
  end do

end subroutine nsi_solide
