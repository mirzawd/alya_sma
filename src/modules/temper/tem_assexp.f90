!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_assexp(&
     elmat,elrhs,eltem,pnode,mnode,npoin,lnods,rhsid)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_assexp
  ! NAME
  !   tem_assexp
  ! DESCRIPTION
  !   This routine assemble the RHS when using an explicit method
  ! USES
  ! USED BY
  !    tem_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  implicit none
  integer(ip), intent(in)  :: pnode,npoin,mnode
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(in)  :: eltem(pnode,*)
  real(rp),    intent(in)  :: elmat(pnode,pnode),elrhs(pnode)
  real(rp),    intent(out) :: rhsid(npoin)
  integer(ip)              :: inode,jnode,ipoin

  do inode=1,pnode
     ipoin=lnods(inode)
     rhsid(ipoin)=rhsid(ipoin)+elrhs(inode)
     do jnode=1,pnode
        rhsid(ipoin)=rhsid(ipoin)&
             -elmat(inode,jnode)*eltem(jnode,2)
     end do
  end do
  
end subroutine tem_assexp
