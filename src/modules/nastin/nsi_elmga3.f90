!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmga3(&
     pnode,ielem,lnods,elcod,elpre,elvel,elfle,elvep,&
     elprp,elgrp,eltem,elmsh,elnor,elcur,elwmean,&
     elbub)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmga3
  ! NAME 
  !    nsi_elmga3
  ! DESCRIPTION
  !    Compute some variables at the Gauss points
  !    ELVEL, ELCOD, ELPRE
  ! OUTPUT
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_kermod, only       :  kfl_adj_prob
  use def_master, only       :  veloc,press,tempe,fleve,&
       &                        velom,kfl_coupl,wmean,&
       &                        ID_NASTIN,ID_ALEFOR,ID_TEMPER,&
       &                        ID_LEVELS,ID_CHEMIC,veloc_forw,tempe_forw
  use def_domain, only       :  ndime,coord
  use def_nastin, only       :  kfl_timei_nsi,kfl_regim_nsi,&
       &                        curle_nsi,norle_nsi,kfl_colev_nsi,&
       &                        kfl_stabi_nsi,vepro_nsi,prpro_nsi,&
       &                        grpro_nsi,kfl_cotem_nsi,&
       &                        nbdfp_nsi,kfl_surte_nsi,bubble_nsi,&
       &                        kfl_bubbl_nsi
  implicit none
  integer(ip), intent(in)    :: pnode
  integer(ip), intent(in)    :: ielem
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(out)   :: elcod(ndime,pnode),elpre(pnode)
  real(rp),    intent(out)   :: elvel(ndime,pnode,*)
  real(rp),    intent(out)   :: elfle(pnode)
  real(rp),    intent(out)   :: elvep(ndime,pnode)
  real(rp),    intent(out)   :: elprp(pnode)
  real(rp),    intent(out)   :: elgrp(ndime,pnode)
  real(rp),    intent(out)   :: eltem(pnode,nbdfp_nsi)
  real(rp),    intent(out)   :: elwmean(pnode,nbdfp_nsi)
  real(rp),    intent(out)   :: elmsh(ndime,pnode)
  real(rp),    intent(out)   :: elnor(ndime,pnode)
  real(rp),    intent(out)   :: elcur(pnode)
  real(rp),    intent(out)   :: elbub
  integer(ip)                :: inode,idime,ipoin,itime

  if( kfl_timei_nsi == 0 ) then
     !
     ! Stationary
     !     
     do inode = 1,pnode
        ipoin = lnods(inode)
        do idime = 1,ndime
           if (kfl_adj_prob == 0) then
	      elvel(idime,inode,1) = veloc(idime,ipoin,1)
	      elvel(idime,inode,2) = 0.0_rp
           else
	      elvel(idime,inode,1) = veloc_forw(idime,ipoin,1)
	      elvel(idime,inode,2) = 0.0_rp
           endif
           elcod(idime,inode)   = coord(idime,ipoin)
        end do
        elpre(inode) = press(ipoin,1)
     end do
  else
     !
     ! Transient
     !     
     do inode = 1,pnode
        ipoin = lnods(inode)
        do idime = 1,ndime
           if (kfl_adj_prob == 0) then
             elvel(idime,inode,1) = veloc(idime,ipoin,1)
             elvel(idime,inode,2) = veloc(idime,ipoin,3)
             do itime=3,nbdfp_nsi
               elvel(idime,inode,itime) = veloc(idime,ipoin,itime+1)
             end do
           else
             elvel(idime,inode,1) = veloc_forw(idime,ipoin,1)
             elvel(idime,inode,2) = veloc(idime,ipoin,3)
             do itime=3,nbdfp_nsi
               elvel(idime,inode,itime) = veloc(idime,ipoin,itime+1)
             end do
           endif
           elcod(idime,inode)   = coord(idime,ipoin)
        end do
        elpre(inode) = press(ipoin,1)
     end do
  end if
  !
  ! Surface tension
  !
  if( kfl_colev_nsi /= 0 ) then
     elfle(1:pnode) = fleve(lnods(1:pnode),1)
  end if
  !
  ! Surface tension
  !
  if( kfl_surte_nsi /= 0 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        elcur(inode) = curle_nsi(ipoin)
        do idime =1,ndime
           elnor(1:ndime,inode) = norle_nsi(1:ndime,ipoin)
        end do
     end do
  end if
  !
  ! Projections
  !
  if( kfl_stabi_nsi > 0 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        do idime = 1,ndime
           elvep(idime,inode) = vepro_nsi(idime,ipoin) 
        end do
        elprp(inode) = prpro_nsi(ipoin) 
     end do
     if( kfl_stabi_nsi == 2 ) then
        do inode = 1,pnode
           ipoin = lnods(inode)
           do idime = 1,ndime
              elgrp(idime,inode) = grpro_nsi(idime,ipoin) 
           end do
        end do
     end if
  end if
  !
  ! Coupling with Temper
  !
  if( kfl_coupl(ID_NASTIN,ID_TEMPER) /= 0 .or. kfl_cotem_nsi == 1 .or. kfl_regim_nsi==3 ) then 
     do inode = 1,pnode
        ipoin = lnods(inode)
        if( kfl_coupl(ID_NASTIN,ID_CHEMIC) /= 0) then
          elwmean(inode,1) = wmean(ipoin,1)
        else
          elwmean(inode,1) = 1.0_rp
        endif
        if (kfl_adj_prob == 0) then
          eltem(inode,1) = tempe(ipoin,1)
        else
          eltem(inode,1) = tempe_forw(ipoin,1)
        endif
     end do

     if (kfl_regim_nsi == 3) then ! Low Mach needs time information
        do itime=2,nbdfp_nsi
           do inode = 1,pnode
              ipoin = lnods(inode)
              if( kfl_coupl(ID_NASTIN,ID_CHEMIC) /= 0) then
                elwmean(inode,itime) = wmean(ipoin,itime+1)
              else
                elwmean(inode,itime) = 1.0_rp
              endif
              if (kfl_adj_prob == 0) then
                eltem(inode,itime) = tempe(ipoin,itime+1)
              else
                eltem(inode,itime) = tempe_forw(ipoin,1)
              endif
           end do
        end do
     endif
  end if
  !
  ! Mesh velocity
  !     
  if( kfl_coupl(ID_NASTIN,ID_ALEFOR) /= 0 ) then ! TESTEO 
     do inode = 1,pnode
        ipoin = lnods(inode)
        do idime = 1,ndime
           elmsh(idime,inode) = velom(idime,ipoin)
        end do
     end do
  end if
  !
  ! Bubble
  !
  if( kfl_bubbl_nsi /= 0 ) then
     elbub = bubble_nsi(ielem)
  end if

end subroutine nsi_elmga3
