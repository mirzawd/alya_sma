!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_elmini(pnode,elmat,elrhs)

  !------------------------------------------------------------------------
  !****f* Levels/lev_elmini
  ! NAME 
  !    lev_elmini
  ! DESCRIPTION
  !    Initialize element matrix and RHS
  ! USES
  ! USED BY
  !    lev_elmope
  !------------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp  
  implicit none
  integer(ip), intent(in)  :: pnode  
  real(rp),    intent(out) :: elmat(pnode,pnode),elrhs(pnode)
  integer(ip)              :: inode,jnode

  do inode=1,pnode
     do jnode=1,pnode
        elmat(jnode,inode)=0.0_rp
     end do
     elrhs(inode)=0.0_rp
  end do

end subroutine lev_elmini
