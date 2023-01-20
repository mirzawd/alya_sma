!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_elmgat(&
     pnode,lnods,eltur,elvel,elcod,elwal,eledd,elust,&
     eltem,elgr2,elsqk,elgrp,elpro,elunk,elfle,elmsh,&
     elprd, elprr, elpgr, ellmax)

  !-----------------------------------------------------------------------
  !****f* Turbul/tur_elmgat
  ! NAME
  !   tur_elmgat
  ! DESCRIPTION
  !    Gather operations for a generic turbulence model
  ! USES
  ! USED BY
  !    tur_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  untur,turmu,tempe,advec,fleve,velom,kfl_coupl,&
       &                      ID_TURBUL,ID_ALEFOR,veloc_forw,untur_forw
  use def_domain, only     :  ndime,mnode,coord,walld
  use def_turbul, only     :  kfl_ustar_tur,kfl_walld_tur,kfl_timei_tur,&
       &                      kfl_cotem_tur,lwnei_tur,&
       &                      ustar_tur,nturb_tur,grve2_tur,&
       &                      kfl_grve2_tur,kfl_grsqk_tur,grsqk_tur,&
       &                      grphi_tur,TUR_K_EPS_PHI_F,iunkn_tur,&
       &                      kfl_ortho_tur,unpro_tur,kfl_colev_tur,&
       &                      kfl_produ_tur,produ_tur, unprr_tur, kfl_shock_tur, &
       &                      unpgr_tur,nbdfp_tur, tur_max_mixlen, TUR_K_EPS_STD
  use def_kermod, only     :  kfl_adj_prob
  implicit none
  integer(ip), intent(in)  :: pnode
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(out) :: eltur(nturb_tur,pnode,nbdfp_tur+1)
  real(rp),    intent(out) :: elwal(pnode)
  real(rp),    intent(out) :: elcod(ndime,pnode)
  real(rp),    intent(out) :: elvel(ndime,pnode)
  real(rp),    intent(out) :: eledd(pnode)
  real(rp),    intent(out) :: elust(pnode)
  real(rp),    intent(out) :: eltem(pnode)
  real(rp),    intent(out) :: elgr2(mnode)
  real(rp),    intent(out) :: elsqk(mnode)
  real(rp),    intent(out) :: elgrp(ndime,pnode)
  real(rp),    intent(out) :: elpro(pnode)
  real(rp),    intent(out) :: elprr(pnode) ! projection of reaction
  real(rp),    intent(out) :: elpgr(ndime, pnode) ! projection of unknown gradient
  real(rp),    intent(out) :: ellmax(pnode)  ! maximum mixing length
  real(rp),    intent(out) :: elunk(pnode,*)
  real(rp),    intent(out) :: elfle(pnode)
  real(rp),    intent(out) :: elmsh(ndime,pnode)
  real(rp),    intent(out) :: elprd(pnode)
  integer(ip)              :: inode,ipoin,iturb,idime,itime
  !
  ! Turbulence variable, wall distance, velocity and coordinates
  !
  do inode = 1,pnode
     ipoin = lnods(inode)
     eledd(inode)         = turmu(ipoin)
     do idime = 1,ndime
        if (kfl_adj_prob == 0) then
          elvel(idime,inode) = advec(idime,ipoin,1)
        else
          elvel(idime,inode) = veloc_forw(idime,ipoin,1)
        endif
        elcod(idime,inode) = coord(idime,ipoin)
     end do
     do iturb = 1,nturb_tur
        if (kfl_adj_prob == 0) then
          eltur(iturb,inode,1) = untur(iturb,ipoin,1)
          eltur(iturb,inode,2) = untur(iturb,ipoin,2)
        else
          eltur(iturb,inode,1) = untur_forw(iturb,ipoin,1)
          eltur(iturb,inode,2) = untur_forw(iturb,ipoin,2)
        endif
     end do
     elunk(inode,1) = eltur(iunkn_tur,inode,1)
  end do
  !
  ! Temperature
  !
  if( kfl_cotem_tur /= 0 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        eltem(inode) = tempe(ipoin,1)
     end do
  end if
  !
  ! Distance to the wall
  !
  if( kfl_walld_tur == 1 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        elwal(inode) = walld(ipoin)
     end do
  else
     elwal = 0.0_rp
  end if
  !
  ! Friction velocity
  !
  if( kfl_ustar_tur == 2 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        elust(inode) = ustar_tur(lwnei_tur(ipoin))
     end do
  end if
  !
  ! Velocity 2nd order gradients
  !
  if( kfl_grve2_tur /= 0 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        elgr2(inode) = grve2_tur(ipoin)
     end do
  end if
  !
  ! Grad(sqrt(k))
  !
  if( kfl_grsqk_tur /= 0 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        elsqk(inode) = grsqk_tur(ipoin)
     end do     
  end if
  !
  ! Time integration
  !
  if( kfl_timei_tur /= 0 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        do itime = 3,nbdfp_tur+1
           do iturb = 1,nturb_tur        
              if (kfl_adj_prob == 0) then
                eltur(iturb,inode,itime) = untur(iturb,ipoin,itime)
              else
                eltur(iturb,inode,itime) = untur(iturb,ipoin,itime)
              endif
           end do
        end do
        if (kfl_adj_prob == 0) then
          elunk(inode,2) = eltur(iunkn_tur,inode,3)
        else
          elunk(inode,2) = untur(iunkn_tur,ipoin,3)
        endif
     end do
  end if
  !
  ! Grad(phi)
  !
  if( TUR_K_EPS_PHI_F .and. iunkn_tur == 4 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        do idime = 1,ndime
           elgrp(idime,inode) = grphi_tur(idime,ipoin)
        end do
     end do     
  end if
  !
  ! Projection
  !
  if( kfl_ortho_tur == 1 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        elpro(inode) = unpro_tur(iunkn_tur,ipoin)
     end do          
  else if( kfl_ortho_tur == 2 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        elpro(inode) = unpro_tur(iunkn_tur,ipoin)
        elprr(inode) = unprr_tur(iunkn_tur,ipoin)
     end do          
  end if
  if (kfl_shock_tur/=0) then
     do inode =1, pnode
        ipoin = lnods(inode)
        do idime = 1,ndime
           elpgr(idime, inode) = unpgr_tur(iunkn_tur, idime, ipoin)
        end do
     end do
  end if
  !
  ! Level set
  !
  if( kfl_colev_tur /= 0 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        elfle(inode) = fleve(ipoin,kfl_colev_tur)
     end do          
  end if
  !
  ! Mesh velocity
  !     
  if( kfl_coupl(ID_TURBUL,ID_ALEFOR) /= 0 ) then  
     do inode = 1,pnode
        ipoin = lnods(inode)
        do idime = 1,ndime
           elmsh(idime,inode) = velom(idime,ipoin)
        end do
     end do
  end if
  !
  ! Smoothed production
  ! 
  if( kfl_produ_tur == 1 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        elprd(inode) = produ_tur(ipoin)
     end do               
  end if
  !
  ! maximum mixing length
  !
  if ( TUR_K_EPS_STD ) then
     do inode =1, pnode 
        ipoin = lnods(inode)
        ellmax(inode) = tur_max_mixlen(ipoin) 
     end do
  end if
end subroutine tur_elmgat
