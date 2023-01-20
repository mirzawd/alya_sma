!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_assdia(&
     pnode,elmat,eltem,elrhs)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmhan
  ! NAME 
  !    nsi_elmhang
  ! DESCRIPTION
  !    Modify element matrix when extension elements are used
  !    Only equation of the first node should be assembled as
  !    it corresponds to its extension test function
  !
  ! USES
  ! USED BY
  !    nsi_elmhang
  !***
  !----------------------------------------------------------------------- 
  use def_kintyp, only         :  ip,rp
  implicit none
  integer(ip), intent(in)      :: pnode
  real(rp),    intent(in)      :: elmat(pnode,pnode)
  real(rp),    intent(in)      :: eltem(pnode,*)
  real(rp),    intent(inout)   :: elrhs(pnode)
  integer(ip)                  :: inode,jnode

  do inode = 1,pnode
     do jnode = 1,pnode
        elrhs(inode) = elrhs(inode) - elmat(inode,jnode) * eltem(jnode,2)
     end do
  end do

end subroutine tem_assdia
