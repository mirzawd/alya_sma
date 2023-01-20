!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_assdia(&
     ndofn,pnode,pevat,lnods,elmat,pmatr)
  !----------------------------------------------------------------------
  !****f* Nastin/nsi_assdia
  ! NAME 
  !    nsi_assdia
  ! DESCRIPTION
  !    Assemble the diagonal diag(A) of A
  ! INPUT
  !    ELMAT ... Elemental matrix ... A_e
  ! OUTPUT
  !    PMATR ... diag(A) = Ass_e diag(A_e)
  ! USES
  ! USED BY
  !    nsi_bouope
  !***
  !----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  implicit none
  integer(ip), intent(in)    :: ndofn,pnode,pevat
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(inout) :: elmat(pevat,pevat)
  real(rp),    intent(inout) :: pmatr(*)
  integer(ip)                :: inode,ipoin,idofl,idofg,idofn

  do inode=1,pnode
     ipoin=lnods(inode)
     idofl=(inode-1)*ndofn
     idofg=(ipoin-1)*ndofn
     do idofn=1,ndofn
        idofl=idofl+1
        idofg=idofg+1
        pmatr(idofg)=pmatr(idofg)+elmat(idofl,idofl)
     end do
  end do

end subroutine nsi_assdia
