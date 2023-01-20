!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_assrhs(elrhs,rhsid,pnode,lnods)

!-----------------------------------------------------------------------
!****f* Turbul/tur_assrhs
! NAME
!   tur_assrhs
! DESCRIPTION
!   This routine assemble the RHS for the turbulence equations.
! USES
! USED BY
!    tur_elmop1
!***
!-----------------------------------------------------------------------
  use      def_kintyp
  implicit none
  integer(ip),  intent(in) :: pnode
  integer(ip),  intent(in) :: lnods(pnode)
  real(rp), intent(in)     :: elrhs(pnode)
  real(rp), intent(inout)  :: rhsid(*)
  integer(ip)              :: inode,ipoin

  do inode=1,pnode
     ipoin=lnods(inode)
     rhsid(ipoin)=rhsid(ipoin)+elrhs(inode)
  end do
  
end subroutine tur_assrhs
