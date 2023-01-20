!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_bougat(&
     ndime,pnodb,npoin,lnodb,bovel,bocod,veloc,coord)
  !-----------------------------------------------------------------------
  !
  ! Gather operations
  !
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  implicit none
  integer(ip), intent(in)  :: ndime,pnodb,npoin
  integer(ip), intent(in)  :: lnodb(pnodb)
  real(rp),    intent(in)  :: veloc(ndime,npoin,*),coord(ndime,npoin)
  real(rp),    intent(out) :: bovel(ndime,pnodb),bocod(ndime,pnodb)
  integer(ip)              :: inodb,ipoin,idime
  !
  ! Current velocity and coordinates
  !
  do inodb=1,pnodb
     ipoin=lnodb(inodb)
     do idime=1,ndime
        bovel(idime,inodb) = veloc(idime,ipoin,1)
        bocod(idime,inodb) = coord(idime,ipoin)
     end do
  end do

end subroutine nsi_bougat
