!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_elmgah(&
     pnode,lnods,eltem,elvel,elcod)
  !------------------------------------------------------------------------
  !****f* Temper/tem_elmgat
  ! NAME 
  !    tem_elmgat
  ! DESCRIPTION
  !    This routine performs gather operation for temperature 
  !    independently enthalpy of temperature is solved
  ! USES
  ! USED BY
  !    tem_heatfl
  !***
  !------------------------------------------------------------------------ 
  use def_kintyp, only     :  ip,rp 
  use def_domain, only     :  ndime,mnode,coord
  use def_master, only     :  tempe,veloc
  use def_temper, only     :  kfl_advec_tem
  use def_temper, only     :  ADR_tem
  use mod_ADR,    only     :  BDF
  implicit none
  integer(ip), intent(in)  :: pnode
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(out) :: eltem(pnode,*)
  real(rp),    intent(out) :: elcod(ndime,mnode),elvel(ndime,mnode)
  integer(ip)              :: inode,ipoin,idime,itime
  !
  ! Current temperature and coordinates
  !
  do inode=1,pnode
     ipoin=lnods(inode)
     eltem(inode,1)=tempe(ipoin,1)
     do idime=1,ndime
        elcod(idime,inode)=coord(idime,ipoin)
     end do  
  end do
  !
  ! Time integration
  !
  if( ADR_tem % kfl_time_integration /= 0 ) then
     do inode=1,pnode
        ipoin=lnods(inode)
        eltem(inode,2)=tempe(ipoin,3)
     end do
     if( ADR_tem % kfl_time_scheme == BDF ) then
        do itime = 3,ADR_tem % kfl_time_order + 1
           do inode=1,pnode
              ipoin=lnods(inode)
              eltem(inode,itime) = tempe(ipoin,itime+1)
           end do
        end do
     end if
  end if 
  !
  ! Coupling with flow equations
  !
  if(kfl_advec_tem==1) then
     do inode=1,pnode
        ipoin=lnods(inode)
        do idime=1,ndime
           elvel(idime,inode)=veloc(idime,ipoin,1)
        end do
     end do
  else if(kfl_advec_tem>=2) then
     call tem_velfun(pnode,elcod,elvel)
  end if

end subroutine tem_elmgah
