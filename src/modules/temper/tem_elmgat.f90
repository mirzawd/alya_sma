!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_elmgat(&
     ielem,pnode,lnods,eltem,elvel,elcod,elmsh)
  !------------------------------------------------------------------------
  !****f* Temper/tem_elmgat
  ! NAME 
  !    tem_elmgat
  ! DESCRIPTION
  !    This routine performs the gather operations
  ! USES
  ! USED BY
  !    tem_elmope
  !    tem_bouope
  !***
  !------------------------------------------------------------------------ 
  use def_kintyp, only     :  ip,rp 
  use def_domain, only     :  ndime,coord
  use def_master, only     :  therm,advec,veloc_forw,velom
  use def_temper, only     :  ADR_tem,kfl_advec_tem
  use mod_ADR,    only     :  BDF
  use def_kermod, only     :  kfl_adj_prob
  use mod_ker_proper
  implicit none
  integer(ip), intent(in)  :: pnode 
  integer(ip), intent(in)  :: ielem
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(out) :: eltem(pnode,*)
  real(rp),    intent(out) :: elcod(ndime,pnode),elvel(ndime,pnode)
  real(rp),    intent(out) :: elmsh(ndime,pnode)
  integer(ip)              :: inode,ipoin,itime
  !
  ! Current temperature and coordinates
  !
  do inode = 1,pnode
     ipoin = lnods(inode)
     eltem(inode,1) = therm(ipoin,1)
     elcod(1:ndime,inode) = coord(1:ndime,ipoin)
  end do

  !
  ! Time integration
  !
  if( ADR_tem % kfl_time_integration /= 0 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        eltem(inode,2) = therm(ipoin,3)
     end do

     if( ADR_tem % kfl_time_scheme == BDF ) then
        do itime = 3,ADR_tem % kfl_time_order + 1
           do inode = 1,pnode 
              ipoin = lnods(inode)
              eltem(inode,itime) = therm(ipoin,itime+1)
           end do
        end do
     end if
  end if
  !
  ! ALE terms
  !
  if( associated(velom) ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        elmsh(1:ndime,inode)=velom(1:ndime,ipoin)
     end do
  else
     elmsh = 0.0_rp 
  end if
  !
  ! Coupling with flow equations
  !
  if( kfl_advec_tem == 1 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        if (kfl_adj_prob == 0_ip) then
           elvel(1:ndime,inode)=advec(1:ndime,ipoin,1)
        else
           elvel(1:ndime,inode)=veloc_forw(1:ndime,ipoin,1)
        endif
     end do

  else if( kfl_advec_tem >= 2 ) then
     call tem_velfun(pnode,elcod,elvel)

  end if

end subroutine tem_elmgat
