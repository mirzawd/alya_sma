!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_post_gather(pnode,lnods,elcon,elcod,elvel)
  !------------------------------------------------------------------------
  !****f* Chemic/chm_post_gather
  ! NAME
  !    chm_post_gather
  ! DESCRIPTION
  !    Gather operations for the combustion models
  ! USES
  ! USED BY
  !    chm_post_scalar_dissipation_rate
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  conce,advec
  use def_domain, only     :  ndime,coord
  use def_chemic, only     :  nclas_chm, &
                              ADR_chm,kfl_advec_chm
  use mod_ADR,    only     :  BDF

  implicit none
  integer(ip), intent(in)  :: pnode
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(out) :: elcon(pnode,nclas_chm,*)
  real(rp),    intent(out) :: elvel(ndime,pnode)
  real(rp),    intent(out) :: elcod(ndime,pnode)
  integer(ip)              :: inode,ipoin,idime,itime,iclas

  !
  ! Concentration and coordinates
  !
  do inode=1,pnode
     ipoin=lnods(inode)
     do iclas=1,nclas_chm
        elcon(inode,iclas,1) = conce(ipoin,iclas,1)
     end do
     do idime=1,ndime
        elcod(idime,inode)   = coord(idime,ipoin)
     end do
  end do
  !
  ! Time integration
  !
  if( ADR_chm(1) % kfl_time_integration == 1 ) then
     do iclas = 1,nclas_chm
        do inode = 1,pnode
           ipoin = lnods(inode)
           elcon(inode,iclas,2) = conce(ipoin,iclas,3)
        end do
     end do

     if( ADR_chm(1) % kfl_time_scheme == BDF ) then
        do iclas = 1,nclas_chm
           do itime = 3,ADR_chm(iclas) % kfl_time_order + 1
              do inode = 1,pnode
                 ipoin = lnods(inode)
                 elcon(inode,iclas,itime) = conce(ipoin,iclas,itime+1)
              end do
           end do
        end do
     end if
  end if

  !
  ! Advection
  !
  if( kfl_advec_chm /= 0 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        do idime = 1,ndime
           elvel(idime,inode) = advec(idime,ipoin,1)
        end do
     end do
  else
     elvel = 0.0_rp
  end if

end subroutine chm_post_gather
