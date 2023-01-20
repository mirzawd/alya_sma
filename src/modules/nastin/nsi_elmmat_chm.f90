!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmmat_chm(pnode,ielem,elrbu)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmmat_chm
  ! NAME 
  !    nsi_elmmat_chm
  ! DESCRIPTION
  !    Modify RHS related to the concentrations terms
  ! USES
  ! USED BY
  !***
  !----------------------------------------------------------------------- 
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
  use def_master, only       :  RhsadjNas_chm
  
  implicit none
  integer(ip), intent(in)    :: pnode,ielem
  real(rp),    intent(inout) :: elrbu(ndime,pnode)
  integer(ip)                :: inode,idime
  
  !----------------------------------------------------------------------
  !
  !                 !  ELRHS = ELRHS - [dRs/du] [Lambda_s]  sent by chemic
  !
  !----------------------------------------------------------------------

  do idime = 1,ndime
    do inode = 1,pnode
      elrbu(idime,inode) = elrbu(idime,inode) - RhsadjNas_chm(ielem)%a(idime, inode)
    enddo
  enddo

  
end subroutine nsi_elmmat_chm
