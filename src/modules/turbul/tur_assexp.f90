!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_assexp(&
     elmat,elrhs,eltur,pnode,mnode,npoin,nturb_tur,&
     iunkn_tur,lnods,rhsid)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_assexp
  ! NAME
  !   tur_assexp
  ! DESCRIPTION
  !   This routine assemble the RHS for explicit solver
  ! USES
  ! USED BY
  !    tur_elmop1
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  implicit none
  integer(ip), intent(in)  :: pnode,npoin,mnode,nturb_tur,iunkn_tur
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(in)  :: eltur(nturb_tur,pnode)
  real(rp),    intent(in)  :: elmat(mnode,mnode),elrhs(pnode)
  real(rp),    intent(out) :: rhsid(npoin)
  integer(ip)              :: inode,jnode,ipoin

  do inode=1,pnode
     ipoin=lnods(inode)
     rhsid(ipoin)=rhsid(ipoin)+elrhs(inode)
     do jnode=1,pnode
        rhsid(ipoin)=rhsid(ipoin)&
             -elmat(inode,jnode)*eltur(iunkn_tur,jnode)
     end do
  end do
  
end subroutine tur_assexp
