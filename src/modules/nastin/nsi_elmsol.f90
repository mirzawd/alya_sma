!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmsol(&
     pnode,pgaus,lnods,gpcar,gpvol,elmal)
  !----------------------------------------------------------------------
  !****f* Nastin/nsi_elmsol
  ! NAME 
  !    nsi_elmsol
  ! DESCRIPTION
  !    Compute the Laplacian matrix
  ! USES
  ! USED BY
  !***
  !----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp 
  use def_domain, only     :  mnode,ndime
  implicit none
  integer(ip), intent(in)  :: pnode,pgaus
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus),gpvol(pgaus)
  real(rp),    intent(out) :: elmal(pnode,pnode)
  integer(ip)              :: inode,jnode,ipoin

  do inode=1,pnode
     ipoin=lnods(inode)
     !if(lmatn(ipoin)==-1) then
        elmal(inode,jnode)=0.0_rp
     !end if
  end do

end subroutine nsi_elmsol
