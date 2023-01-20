!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine extcog(iboun,ielem,pnodb,pnode,exwor)
  !-----------------------------------------------------------------------
  !
  ! This routine calculates the normal baloc at the center of
  ! gravity of the boundary element and assign it to exwor.
  ! The orientation is not a problem as the point belongs to
  ! two adjacent surfaces.
  !
  !   exwor <- baloc
  !
  !    /|\     /|\
  !     |       |
  !     o-------+-------o
  !           iboun
  !
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime,ndimb,mnode,coord,lnods,&
       &                      lnodb,ltypb,elmar
  use mod_bouder
  implicit  none
  integer(ip), intent(in)  :: iboun,ielem,pnodb,pnode
  real(rp),    intent(out) :: exwor(ndime)
  integer(ip)              :: idime,inodb,ipoin,iblty,inode
  real(rp)                 :: xnorm,eucta
  real(rp)                 :: baloc(ndime,ndime)
  real(rp)                 :: bocod(ndime,pnodb),elcod(ndime,mnode)
  !
  ! Identifies the coordinates of the boundary nodes
  !
  iblty=ltypb(iboun)
  do inodb=1,pnodb
     ipoin=lnodb(inodb,iboun)
     do idime=1,ndime
        bocod(idime,inodb)=coord(idime,ipoin)
     end do
  end do
  do inode=1,pnode
     ipoin=lnods(inode,ielem)
     do idime=1,ndime
        elcod(idime,inode)=coord(idime,ipoin)
     end do
  end do
  !
  ! Coordinates of the center of gravity
  !
  call bouder(&
       pnodb,ndime,ndimb,elmar(iblty)%dercg,bocod,baloc,eucta) ! BALOC
  call chenor(pnode,baloc,bocod,elcod)                ! Check BALOC
  !
  ! Set exwor=baloc and normalize it
  !      
  xnorm=0.0_rp
  do idime=1,ndime
     exwor(idime)=baloc(idime,ndime)
     xnorm=xnorm+exwor(idime)*exwor(idime)
  end do
  do idime=1,ndime
     exwor(idime)=exwor(idime)/sqrt(xnorm)
  end do

end subroutine extcog
