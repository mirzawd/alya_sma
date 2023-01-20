!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmstr(&
     pgaus,pnode,ndofn,ielem,lnods, chale,elvel,gpadv,gpvis,&
     gpden,rmomu,rmom2,gprhs,gpvel,gpcar,gpsp1,&
     gpstrm,gpstrc,gpvep,gpgrp,gpst1,gprhs_sgs,gpsha,gpgve)
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_elmstr
  ! NAME 
  !    nsi_elmstr
  ! DESCRIPTION
  !    This subroutine calculates the strong residual
  ! USES
  ! USED BY
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_master, only       :  press,press_forw
  use def_domain, only       :  ndime,mnode
  use def_kermod, only       :  kfl_adj_prob
  use def_nastin, only       :  ncomp_nsi
  implicit none
  integer(ip), intent(in)    :: pgaus,pnode,ndofn,ielem
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(in)    :: chale(2)
  real(rp),    intent(in)    :: elvel(ndime,pnode,ncomp_nsi)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gprhs(ndofn,pgaus)
  real(rp),    intent(in)    :: gprhs_sgs(ndofn,pgaus)
  real(rp),    intent(in)    :: gpvel(ndime,pgaus)
  real(rp),    intent(inout) :: gpadv(ndime,pgaus)
  real(rp),    intent(in)    :: gpvis(pgaus)
  real(rp),    intent(in)    :: gpden(pgaus)
  real(rp),    intent(in)    :: rmomu(pnode,pgaus)
  real(rp),    intent(in)    :: rmom2(ndime,ndime,pnode,pgaus)
  real(rp),    intent(out)   :: gpsp1(pgaus)
  real(rp),    intent(out)   :: gpstrm(ndime,pgaus)
  real(rp),    intent(out)   :: gpstrc(pgaus)
  real(rp),    intent(in)    :: gpvep(ndime,pgaus)
  real(rp),    intent(in)    :: gpgrp(ndime,pgaus)
  real(rp),    intent(out)   :: gpst1(pgaus)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpgve(ndime,ndime,pgaus)
  integer(ip)                :: idime,igaus,inode,jdime,ipoin
  real(rp)                   :: rmomu1(pnode,pgaus),rmomu2(ndime,ndime,pnode,pgaus),elpre_forw(pnode),gpgpr(ndime,pgaus)
  real(rp)                   :: fact0

  !----------------------------------------------------------------------
  !
  ! Initialization
  !
  !----------------------------------------------------------------------
  
  do igaus = 1,pgaus
    gpstrc(igaus) = 0.0_rp
    do idime = 1,ndime
      gpstrm(idime,igaus) = 0.0_rp
    enddo
  enddo
  
  do igaus = 1,pgaus
    do inode = 1,pnode
      rmomu1(inode,igaus) = rmomu(inode,igaus)
      do idime = 1,ndime
	do jdime = 1,ndime
	  rmomu2(idime,jdime, inode,igaus) = rmom2(idime,jdime, inode,igaus)
	enddo
      enddo
    end do
  end do
  
  !
  ! calculate grad(p) here (without considering as a subroutine input) beacase in the case of adjoint the pressure is different
  !
  do inode = 1,pnode
    ipoin = lnods(inode)
    if (kfl_adj_prob == 0) then
      elpre_forw(inode) = press(ipoin,1)
    else
      elpre_forw(inode) = press_forw(ipoin,1)
    endif
  end do
    
  do igaus = 1,pgaus
     do idime = 1,ndime
        gpgpr(idime,igaus)   = 0.0_rp
     end do
     do inode = 1,pnode
        do idime = 1,ndime
           gpgpr(idime,igaus)   = gpgpr(idime,igaus) + elpre_forw(inode) * gpcar(idime,inode,igaus)
        end do
     end do
  end do
  
  !
  !  rmom2 = rmom2 - (terms related to the exact linearization ADJ ) just in order to calculate gpstrm
  !  
  if( kfl_adj_prob == 1 ) then
     if( ndime == 2 ) then
        do igaus = 1,pgaus           
           do inode = 1,pnode
              fact0 = gpden(igaus) * gpsha(inode,igaus)
              rmomu2(1,1,inode,igaus) = rmomu2(1,1,inode,igaus) - fact0 * gpgve(1,1,igaus)   ! rho * ux * d(adv_x)/dx
              rmomu2(1,2,inode,igaus) = rmomu2(1,2,inode,igaus) - fact0 * gpgve(2,1,igaus)   ! rho * uy * d(adv_x)/dy
              rmomu2(2,1,inode,igaus) = rmomu2(2,1,inode,igaus) - fact0 * gpgve(1,2,igaus)   ! rho * ux * d(adv_y)/dx
              rmomu2(2,2,inode,igaus) = rmomu2(2,2,inode,igaus) - fact0 * gpgve(2,2,igaus)   ! rho * uy * d(adv_y)/dy
           end do
        end do
     else
        do igaus = 1,pgaus           
           do inode = 1,pnode
              fact0 = gpden(igaus) * gpsha(inode,igaus)
              rmomu2(1,1,inode,igaus) = rmomu2(1,1,inode,igaus) - fact0 * gpgve(1,1,igaus)   ! rho * ux * d(adv_x)/dx
              rmomu2(1,2,inode,igaus) = rmomu2(1,2,inode,igaus) - fact0 * gpgve(2,1,igaus)   ! rho * uy * d(adv_x)/dy
              rmomu2(1,3,inode,igaus) = rmomu2(1,3,inode,igaus) - fact0 * gpgve(3,1,igaus)   ! rho * uz * d(adv_x)/dz
              rmomu2(2,1,inode,igaus) = rmomu2(2,1,inode,igaus) - fact0 * gpgve(1,2,igaus)   ! rho * ux * d(adv_y)/dx
              rmomu2(2,2,inode,igaus) = rmomu2(2,2,inode,igaus) - fact0 * gpgve(2,2,igaus)   ! rho * uy * d(adv_y)/dy
              rmomu2(2,3,inode,igaus) = rmomu2(2,3,inode,igaus) - fact0 * gpgve(3,2,igaus)   ! rho * uz * d(adv_y)/dz
              rmomu2(3,1,inode,igaus) = rmomu2(3,1,inode,igaus) - fact0 * gpgve(1,3,igaus)   ! rho * ux * d(adv_z)/dx
              rmomu2(3,2,inode,igaus) = rmomu2(3,2,inode,igaus) - fact0 * gpgve(2,3,igaus)   ! rho * uy * d(adv_z)/dy
              rmomu2(3,3,inode,igaus) = rmomu2(3,3,inode,igaus) - fact0 * gpgve(3,3,igaus)   ! rho * uz * d(adv_z)/dz
           end do
        end do
     end if
  end if
  
  !----------------------------------------------------------------------
  !
  ! GPSTRM: 
  ! For Forward:        Res(u) = rho*(u*d/dx)u + sig*u + grad(p) - div[2*mu*eps(u)] + rho*(u^(n+1)-u^(n))/dt
  !                            = rmomu1*vel + grad(p) - rho*u^(n)/dt
  ! 
  ! For Adjoint:        Res(u) = rho*(u*d/dx)u + sig*u + grad(p) - div[2*mu*eps(u)] 
  !                            = rmomu1*vel + grad(p)
  !
  !----------------------------------------------------------------------
  
  do igaus = 1,pgaus
    do idime = 1,ndime
	gpstrm(idime,igaus) = gpgpr(idime,igaus)
	do inode=1,pnode
	  gpstrm(idime,igaus) = gpstrm(idime,igaus)&
				+ rmomu1(inode,igaus) * elvel(idime,inode,1)
	  do jdime = 1,ndime
	    gpstrm(idime,igaus) = gpstrm(idime,igaus)&
				+ rmomu2(idime,jdime,inode,igaus) * elvel(jdime,inode,1)
	  enddo
	end do
    end do
  end do
    
  !----------------------------------------------------------------------
  !
  ! GPSTRC: 
  !        Res(c) = v.grad(u)
  !
  !----------------------------------------------------------------------

  do igaus = 1,pgaus
    do inode=1,pnode
	if( ndime == 2 ) then
	  gpstrc(igaus) = gpstrc(igaus)&
			+ elvel(1,inode,1)*gpcar(1,inode,igaus) + elvel(2,inode,1)*gpcar(2,inode,igaus)
	else
	  gpstrc(igaus) = gpstrc(igaus)&
			+ elvel(1,inode,1)*gpcar(1,inode,igaus) + elvel(2,inode,1)*gpcar(2,inode,igaus)+ &	
			  elvel(3,inode,1)*gpcar(3,inode,igaus)
	endif
    end do
  end do
  
end subroutine nsi_elmstr
