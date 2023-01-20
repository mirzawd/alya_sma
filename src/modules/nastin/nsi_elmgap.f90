!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmgap(&
     pnode,pmate,lnods,elcod,elpre,elvel,eltem)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmgap
  ! NAME 
  !    nsi_elmgap
  ! DESCRIPTION
  !    Compute some variables at the Gauss points
  !    ELVEL, ELCOD, ELPRE
  ! OUTPUT
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_master, only       :  veloc,press,tempe,&
       &                        veloc_forw,press_forw
  use def_domain, only       :  ndime,coord
  use def_nastin, only       :  kfl_cotem_nsi,&
       &                        kfl_regim_nsi, kfl_surte_nsi
  use def_kermod, only : kfl_adj_prob
  
  implicit none
  integer(ip), intent(in)    :: pnode,pmate
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(out)   :: elcod(ndime,pnode)
  real(rp),    intent(out)   :: elpre(pnode)
  real(rp),    intent(out)   :: elvel(ndime,pnode)
  real(rp),    intent(out)   :: eltem(pnode)
  integer(ip)                :: inode,idime,ipoin
  !
  ! ELCOD <- COORD
  !
  do inode = 1,pnode
     ipoin = lnods(inode)
     do idime = 1,ndime
        elcod(idime,inode) = coord(idime,ipoin)
     end do
  end do
  !
  ! ELVEL <- VELOC
  !
  do inode = 1,pnode
     ipoin = lnods(inode)
     do idime = 1,ndime
       if (kfl_adj_prob == 0) then
         elvel(idime,inode) = veloc(idime,ipoin,1)
       else
         elvel(idime,inode) = veloc_forw(idime,ipoin,1)
       endif
     end do
  end do
  !
  ! ELPRE <- PRESS
  !
  do inode = 1,pnode
     ipoin = lnods(inode)
     if (kfl_adj_prob == 0) then
       elpre(inode) = press(ipoin,1)
     else
       elpre(inode) = press_forw(ipoin,1)
     endif
  end do
  !
  ! ELTEM <- TEMPE: Coupling with Temper
  !
  if( kfl_cotem_nsi == 1 .or. kfl_regim_nsi == 3 .and. kfl_surte_nsi /= 2_ip .and. associated(tempe)) then 
     do inode = 1,pnode
        ipoin = lnods(inode)
        eltem(inode) = tempe(ipoin,1)
     end do
  end if

end subroutine nsi_elmgap
