!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_elmchl(&
     ielem,pelty,pnode,plapl,pmate,lnods,elcod,eltem,&
     elvel,gpcar,gphes,chale)
  !------------------------------------------------------------------------
  !****f* Temper/tem_elmchl
  ! NAME 
  !    tem_elmchl
  ! DESCRIPTION
  ! USES
  ! USED BY
  !    tem_matrix
  !------------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime,ntens
  implicit none
  integer(ip), intent(in)  :: ielem,pelty,pnode,plapl,pmate
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(in)  :: elcod(ndime,pnode),eltem(pnode)
  real(rp),    intent(in)  :: elvel(ndime,pnode)
  real(rp),    intent(in)  :: gpcar(ndime,pnode),gphes(ntens,pnode)
  real(rp),    intent(out) :: chale(2)
!  integer(ip)              :: inode,knode,jnode,idime,ipoin,klapl
!  real(rp)                 :: gpdir(3),gpcod(3),gpgrt(3),rnode,gpgr2(3)
!  real(rp)                 :: gpvol,gptem,gpden,gpcon,gpdif
  real(rp)                 :: hleng(3)
!  real(rp)                 :: gpsph,gprea
!  real(rp)                 :: vx1,vy1,vx2,vy2,m,x1,y1,x2,y2
!  real(rp)                 :: xmin1,xmin2,ymin1,ymin2,det,dinor

  call velchl(pnode,elcod,elvel,chale,hleng) 
  return

end subroutine tem_elmchl
