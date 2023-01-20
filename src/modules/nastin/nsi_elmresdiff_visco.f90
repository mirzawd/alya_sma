!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup NastinInput
!> @{
!> @file    nsi_elmresdiff_visco.f90
!> @author  Guillaume Houzeaux
!> @brief   Compute the derivation of GPRES r. t. design variables
!
!-------------------------------------------------------------------
!        +-           -+    +-                      -+   +- -+   +-  -+
!        | elmresudiff  |   | elauudiff    elaupdiff |   | u |   | elrbudiff |
!        |              | = |                        | * |   | - |           |
!        | elmresbdiff  |   | elapudiff    elappdiff |   | p |   | elrbpdiff |
!        +-           -+    +-                      -+   +- -+   +-  -+
!-------------------------------------------------------------------
!
!> @} 
!-----------------------------------------------------------------------
subroutine nsi_elmresdiff_visco(pnode,pgaus,pevat,gpvol,gpsha,lnods, &
			  p1vec,p2vec,elvel,wgrgr,h, &
			  gpsp1,gpsp2,rmomu,gpden,gpadv,gpcar,rcont,elresudiff,elrespdiff)

  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime,mnode
!  use def_nastin, only     :  dtinv_nsi, pabdf_nsi
  use def_parame, only     :  pi
  use def_master

  implicit none
  integer(ip), intent(in)    :: pnode,pgaus,pevat
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(in)    :: p1vec(pnode,pgaus)
  real(rp),    intent(in)    :: elvel(ndime, pnode)
  real(rp),    intent(in)    :: p2vec(ndime,pnode,pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus),gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: wgrgr(pnode,pnode,pgaus)
  real(rp),    intent(in)    :: h
  real(rp),    intent(in)    :: gpsp1(pgaus)                      ! tau1'
  real(rp),    intent(in)    :: gpsp2(pgaus)                      ! tau2'
  real(rp),    intent(in)    :: rmomu(pnode,pgaus)
  real(rp),    intent(in)    :: gpden(pgaus)
  real(rp),    intent(in)    :: gpadv(ndime,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: rcont(ndime,pnode,pgaus)
  real(rp),    intent(out)   :: elresudiff(ndime,pnode)
  real(rp),    intent(out)   :: elrespdiff(pnode)

  integer(ip)                :: idime,inode,jnode,igaus,idofn,jdofn,idofv,ipoin
  integer(ip)                :: idof1,idof2,idof3,jdof1,jdof2,jdof3,ind
  real(rp)                   :: elauudiff(pnode*ndime,pnode*ndime)
  real(rp)                   :: elaupdiff(pnode*ndime,pnode)
  real(rp)                   :: elappdiff(pnode,pnode)
  real(rp)                   :: elapudiff(pnode,pnode*ndime)
  real(rp)                   :: elrbudiff(ndime,pnode)
  real(rp)	             :: elrbpdiff(pnode)
  real(rp)                   :: elraudiff(ndime,pnode)
  real(rp)	             :: elrapdiff(pnode)
  real(rp)	             :: fact0, fact1,fact2,fact3,fact4,fact6,elvel1(ndime*pnode),elpre1(pnode)
  real(rp)                   :: auuprodu(ndime*pnode),aupprodp(ndime*pnode)
  real(rp)                   :: apuprodu(pnode),appprodp(pnode)
  real(rp)                   :: gpsp1diff(pgaus)                      ! d(tau1')/d(mu)
  real(rp)    		     :: gpsp2diff(pgaus)                      ! d(tau2')/d(mu)
  real(rp)                   :: p1vecdiff(pnode,pgaus)
  real(rp)                   :: p2vecdiff(ndime,pnode,pgaus)
  real(rp)                   :: rmomu1(pnode,pgaus)
  
  !----------------------------------------------------------------------
  !
  ! Initialization
  !
  !----------------------------------------------------------------------

  do inode = 1,pnode
    elrbpdiff(inode) = 0.0_rp
    elrespdiff(inode) = 0.0_rp
    do idime = 1,ndime
	elrbudiff(idime,inode) = 0.0_rp
	elresudiff(idime,inode) = 0.0_rp
    end do
    do jnode = 1,pnode
	elappdiff(jnode,inode) = 0.0_rp
    end do
  end do
  do idofn = 1,pevat
    do jdofn = 1,pevat
	elauudiff(jdofn,idofn) = 0.0_rp
    end do
    do jnode = 1,pnode
	elaupdiff(idofn,jnode) = 0.0_rp
	elapudiff(jnode,idofn) = 0.0_rp
    end do
  end do
  
  
  do inode = 1,pnode
    apuprodu(inode) = 0.0_rp
    appprodp(inode) = 0.0_rp
    do idime = 1,ndime
      ind = (inode-1)*ndime + idime
      auuprodu(ind) = 0.0_rp
      aupprodp(ind) = 0.0_rp
    enddo
  enddo
  
  !
  ! elvel ---> elvel1 
  !
  do inode = 1,pnode
    do idime = 1,ndime
      ind = (inode-1)*ndime + idime
      elvel1(ind) = elvel(idime,inode)
    enddo
  enddo
  !
  ! press_nsi ---> elpre1 
  !
  do inode = 1,pnode
    ipoin = lnods(inode)
    elpre1(inode) = press_forw(ipoin,1)
  end do
  !
  !  gpsp1diff and gpsp2diff
  !
  do igaus = 1,pgaus
    gpsp1diff(igaus) = -4.0_rp*gpsp1(igaus)*gpsp1(igaus)/(h*h)
    gpsp2diff(igaus) = 4.0_rp
  enddo
  !
  !  rmomu --->  rmomu1
  !      
  do inode = 1,pnode
    do igaus = 1,pgaus
      rmomu1(inode,igaus) = rmomu(inode,igaus)
    enddo
  enddo
  
  !----------------------------------------------------------------------
  !
  ! Test functions
  !
  !----------------------------------------------------------------------
  !
  ! P1VEC =  v * [ (tau1'/tau1) - tau1' * sig ] + tau1' * rho * (uc.grad)v
  ! P2VEC =  tau1' * grad(q)  ( and * rho if in Low-Mach, but not yet)
  !
  if( ndime == 2 ) then

    do igaus = 1,pgaus
	fact1 = gpsp1diff(igaus)*gpden(igaus)*gpvol(igaus)
	fact3 = gpsp1diff(igaus)*gpvol(igaus)
	do inode = 1,pnode
	  p1vecdiff(inode,igaus)   = fact1 *    & 
		&                 ( gpadv(1,igaus) * gpcar(1,inode,igaus) &
		&                 + gpadv(2,igaus) * gpcar(2,inode,igaus) )
	  p2vecdiff(1,inode,igaus) = fact3 * gpcar(1,inode,igaus)
	  p2vecdiff(2,inode,igaus) = fact3 * gpcar(2,inode,igaus)
	end do
    end do

  else

    do igaus = 1,pgaus
	fact1 = gpsp1diff(igaus)*gpden(igaus)*gpvol(igaus)
	fact3 = gpsp1diff(igaus)*gpvol(igaus)
	do inode = 1,pnode
	  p1vecdiff(inode,igaus)   = fact1 *    &
		&                 ( gpadv(1,igaus) * gpcar(1,inode,igaus) &
		&                 + gpadv(2,igaus) * gpcar(2,inode,igaus) &
		&                 + gpadv(3,igaus) * gpcar(3,inode,igaus) )
	  p2vecdiff(1,inode,igaus) = fact3 * gpcar(1,inode,igaus)
	  p2vecdiff(2,inode,igaus) = fact3 * gpcar(2,inode,igaus)
	  p2vecdiff(3,inode,igaus) = fact3 * gpcar(3,inode,igaus)
	end do
    end do
    
  end if
  !
  ! RMOMU1 = RMOMU1 - rho/(dt*theta)*u (since gprhs is zero we have to cancle the time temrs from RMOMU1)
  !
!   do igaus = 1,pgaus
!     fact1 = dtinv_nsi * pabdf_nsi(1) * gpden(igaus)
!     do inode = 1,pnode
! 	rmomu1(inode,igaus) = rmomu1(inode,igaus) - fact1 * gpsha(inode,igaus)
!     end do
!   end do

  !----------------------------------------------------------------------
  !
  ! elauudiff
  !
  !----------------------------------------------------------------------
  if( ndime == 2 ) then
    do igaus = 1,pgaus
	fact0 = gpsp2diff(igaus) * gpvol(igaus)                        ! (1) + (5)
	fact6 = gpvol(igaus)
	do inode = 1,pnode
	  idof1 = 2*inode-1
	  idof2 = idof1+1
	  fact1 = fact0 * gpcar(1,inode,igaus) 
	  fact2 = fact0 * gpcar(2,inode,igaus)
	  do jnode = 1,pnode              
	      jdof1              =   2*jnode-1
	      jdof2              =   jdof1+1
	      fact4              =   p1vecdiff(inode,igaus) * rmomu1(jnode,igaus) &            ! (2):       ( rmomu1(u) , p1vec(v) ) 
				    +fact6 * wgrgr(inode,jnode,igaus)                          ! (7):       ( mu dui/dxj , dvi/dxj )
	      elauudiff(idof1,jdof1) = elauudiff(idof1,jdof1) + fact1 * rcont(1,jnode,igaus) + fact4  ! Auu_xx
	      elauudiff(idof2,jdof1) = elauudiff(idof2,jdof1) + fact2 * rcont(1,jnode,igaus)          ! Auu_yx
	      elauudiff(idof1,jdof2) = elauudiff(idof1,jdof2) + fact1 * rcont(2,jnode,igaus)          ! Auu_xy
	      elauudiff(idof2,jdof2) = elauudiff(idof2,jdof2) + fact2 * rcont(2,jnode,igaus) + fact4  ! Auu_yy
	  end do
	end do
    end do
  else
    do igaus = 1,pgaus
	fact0 = gpsp2diff(igaus) * gpvol(igaus)                        ! (1) + (5)
	fact6 = gpvol(igaus)
	do inode = 1,pnode
	  idof1 = 3*inode-2
	  idof2 = idof1+1
	  idof3 = idof2+1
	  fact1 = fact0 * gpcar(1,inode,igaus) 
	  fact2 = fact0 * gpcar(2,inode,igaus)
	  fact3 = fact0 * gpcar(3,inode,igaus)
	  do jnode = 1,pnode              
	      jdof1              =  3*jnode-2
	      jdof2              =  jdof1+1
	      jdof3              =  jdof2+1
	      fact4              =  p1vecdiff(inode,igaus) * rmomu1(jnode,igaus) &             ! (2):       ( rmomu1(u) , p1vec(v) ) 
				    +fact6 * wgrgr(inode,jnode,igaus)                         ! (7):       ( mu*dui/dxk , dv/dxk )
	  elauudiff(idof1,jdof1) = elauudiff(idof1,jdof1) + fact1 * rcont(1,jnode,igaus) + fact4  ! Auu_xx
	  elauudiff(idof2,jdof1) = elauudiff(idof2,jdof1) + fact2 * rcont(1,jnode,igaus)          ! Auu_yx
	  elauudiff(idof3,jdof1) = elauudiff(idof3,jdof1) + fact3 * rcont(1,jnode,igaus)          ! Auu_zx
	  elauudiff(idof1,jdof2) = elauudiff(idof1,jdof2) + fact1 * rcont(2,jnode,igaus)          ! Auu_xy
	  elauudiff(idof2,jdof2) = elauudiff(idof2,jdof2) + fact2 * rcont(2,jnode,igaus) + fact4  ! Auu_yy
	  elauudiff(idof3,jdof2) = elauudiff(idof3,jdof2) + fact3 * rcont(2,jnode,igaus)          ! Auu_zy
	  elauudiff(idof1,jdof3) = elauudiff(idof1,jdof3) + fact1 * rcont(3,jnode,igaus)          ! Auu_xz
	  elauudiff(idof2,jdof3) = elauudiff(idof2,jdof3) + fact2 * rcont(3,jnode,igaus)          ! Auu_yz
	  elauudiff(idof3,jdof3) = elauudiff(idof3,jdof3) + fact3 * rcont(3,jnode,igaus) + fact4  ! Auu_zz
	  end do
	end do
    end do
  end if
  
  !----------------------------------------------------------------------
  !
  ! elaupdiff
  !
  !----------------------------------------------------------------------
  !
  ! Pressure: - ( p , div(v) ) + ( grad(p) , p1vec(v)-v )
  !  
  if( ndime == 2 ) then
    do igaus = 1,pgaus
	do jnode = 1,pnode
	  do inode = 1,pnode
	      fact1              = p1vecdiff(inode,igaus)
	      idofv              = (inode-1)*ndime+1
	      elaupdiff(idofv,jnode) = elaupdiff(idofv,jnode) + fact1*gpcar(1,jnode,igaus)
	      idofv              = idofv+1
	      elaupdiff(idofv,jnode) = elaupdiff(idofv,jnode) + fact1*gpcar(2,jnode,igaus)
	  end do
	end do
    end do
  else
    do igaus = 1,pgaus
	do jnode = 1,pnode
	  do inode = 1,pnode
	      fact1              = p1vecdiff(inode,igaus)
	      idofv              = (inode-1)*ndime+1
	      elaupdiff(idofv,jnode) = elaupdiff(idofv,jnode) + fact1 * gpcar(1,jnode,igaus)
	      idofv              = idofv+1
	      elaupdiff(idofv,jnode) = elaupdiff(idofv,jnode) + fact1 * gpcar(2,jnode,igaus)
	      idofv              = idofv+1
	      elaupdiff(idofv,jnode) = elaupdiff(idofv,jnode) + fact1 * gpcar(3,jnode,igaus)
	  end do
	end do
    end do
  end if
  !----------------------------------------------------------------------
  !
  ! elapudiff
  !
  !----------------------------------------------------------------------
  !
  ! ( div(u) , (tau2^{-1}*tau2')*q )   --- or ( div (rho u), q) if Low-Mach 
  ! + ( rho*(uc.grad)u + 2*rho*(w x u) + sig*u -div[2*mu*eps(u)], tau1' grad(q) )  (only laplacian form of div[2 mu eps(u)] )
  !
  if( ndime == 2 ) then
    do igaus = 1,pgaus
! 	    fact1 = dtinv_nsi * pabdf_nsi(1) * gpden(igaus)*(1-tiprst_nsi) !!! for the time pressure stabilization
	do inode = 1,pnode
	  idof1 = 2*inode-1
	  idof2 = idof1+1
	  do jnode = 1,pnode
	      elapudiff(jnode,idof1) = elapudiff(jnode,idof1) &
		  &                                  + rmomu1(inode,igaus)   * p2vecdiff(1,jnode,igaus)
! 		      &                                  - fact1*gpsha(inode,igaus)   * p2vecdiff(1,jnode,igaus)
	      elapudiff(jnode,idof2) = elapudiff(jnode,idof2) &
		  &                                  + rmomu1(inode,igaus)   * p2vecdiff(2,jnode,igaus)
! 		      &                                  - fact1*gpsha(inode,igaus)   * p2vecdiff(2,jnode,igaus)
	  end do
	end do
    end do
  else
    do igaus = 1,pgaus
! 	    fact1 = dtinv_nsi * pabdf_nsi(1) * gpden(igaus)*(1-tiprst_nsi) !!! for the time pressure stabilization
	do inode = 1,pnode
	  idof1 = 3*inode-2
	  idof2 = idof1+1
	  idof3 = idof2+1
	  do jnode = 1,pnode
	      elapudiff(jnode,idof1) = elapudiff(jnode,idof1) &
		  &                                  + rmomu1(inode,igaus)   * p2vecdiff(1,jnode,igaus)
! 		      &                                  - fact1*gpsha(inode,igaus)   * p2vecdiff(1,jnode,igaus)
	      elapudiff(jnode,idof2) = elapudiff(jnode,idof2) &
		  &                                  + rmomu1(inode,igaus)   * p2vecdiff(2,jnode,igaus)
! 		      &                                  - fact1*gpsha(inode,igaus)   * p2vecdiff(2,jnode,igaus)
	      elapudiff(jnode,idof3) = elapudiff(jnode,idof3) &
		  &                                  + rmomu1(inode,igaus)   * p2vecdiff(3,jnode,igaus)
! 		      &                                  - fact1*gpsha(inode,igaus)   * p2vecdiff(3,jnode,igaus)
	  end do
	end do
    end do
  end if
  !----------------------------------------------------------------------
  !
  ! elappdiff
  !
  !----------------------------------------------------------------------
  !
  ! Pressure: ( grad(p) , tau1' grad(q) )
  ! 
  if( ndime == 2 ) then
    do igaus=1,pgaus
	do inode=1,pnode
	  do jnode=inode+1,pnode
	      fact1 =  p2vecdiff(1,jnode,igaus) * gpcar(1,inode,igaus)&
		  & + p2vecdiff(2,jnode,igaus) * gpcar(2,inode,igaus)
	      elappdiff(jnode,inode) = elappdiff(jnode,inode) + fact1
	      elappdiff(inode,jnode) = elappdiff(inode,jnode) + fact1
	  end do
	  fact1 =  p2vecdiff(1,inode,igaus) * gpcar(1,inode,igaus)&
		& + p2vecdiff(2,inode,igaus) * gpcar(2,inode,igaus)
	  elappdiff(inode,inode) = elappdiff(inode,inode) + fact1
	end do
    end do
  else
    do igaus=1,pgaus
	do inode=1,pnode
	  do jnode=inode+1,pnode
	      fact1 =  p2vecdiff(1,jnode,igaus) * gpcar(1,inode,igaus)&
		  & + p2vecdiff(2,jnode,igaus) * gpcar(2,inode,igaus)&
		  & + p2vecdiff(3,jnode,igaus) * gpcar(3,inode,igaus)
	      elappdiff(jnode,inode) = elappdiff(jnode,inode) + fact1
	      elappdiff(inode,jnode) = elappdiff(inode,jnode) + fact1
	  end do
	  fact1 =  p2vecdiff(1,inode,igaus) * gpcar(1,inode,igaus)&
		& + p2vecdiff(2,inode,igaus) * gpcar(2,inode,igaus)&
		& + p2vecdiff(3,inode,igaus) * gpcar(3,inode,igaus)
	  elappdiff(inode,inode) = elappdiff(inode,inode) + fact1
	end do
    end do
  end if
  
  !----------------------------------------------------------------------
  !
  ! elraudiff = elauudiff * u  + elaupdiff * p 
  ! elrapdiff = elapudiff * u  + elappdiff * p 
  !
  !----------------------------------------------------------------------
  
  !elauudiff * u
  call mbvab0(auuprodu,elauudiff,elvel1,pevat,pevat)
  !elaupdiff * p
  call mbvab0(aupprodp,elaupdiff,elpre1,pevat,pnode)
  !elapudiff * u
  call mbvab0(apuprodu,elapudiff,elvel1,pnode,pevat)
  !elappdiff * p
  call mbvab0(appprodp,elappdiff,elpre1,pnode,pnode)
  
  ! auuprodu + aupprodp --> elraudiff
  ! apuprodu + appprodp --> elrapdiff
  do inode = 1,pnode
    elrapdiff(inode) = apuprodu(inode)+appprodp(inode)
    do idime = 1,ndime
      ind = (inode-1)*ndime + idime
      elraudiff(idime,inode) = auuprodu(ind)+aupprodp(ind)
    enddo
  enddo

  
  ! elauudiff --> 0.0 (for assembly at the end)
  ! elaupdiff --> 0.0 (for assembly at the end)
  ! elapudiff --> 0.0 (for assembly at the end)
  ! elappdiff --> 0.0 (for assembly at the end)
  do idofn = 1,pevat
    do jdofn = 1,pevat
	elauudiff(jdofn,idofn) = 0.0_rp
    end do
    do jnode = 1,pnode
	elaupdiff(idofn,jnode) = 0.0_rp
	elapudiff(jnode,idofn) = 0.0_rp
    end do
  end do
  do inode = 1,pnode
    do jnode = 1,pnode
	elappdiff(jnode,inode) = 0.0_rp
    end do
  end do

  !----------------------------------------------------------------------
  !
  ! elrbudiff and elrbpdiff
  !
  !----------------------------------------------------------------------

!   if( ndime == 2 ) then
!     do igaus = 1,pgaus
! 	do inode = 1,pnode
! 	      elrbudiff(1,inode) = elrbudiff(1,inode) + p1vecdiff(inode,igaus)   * gprhs(1,igaus)
! 	      elrbudiff(2,inode) = elrbudiff(2,inode) + p1vecdiff(inode,igaus)   * gprhs(2,igaus)
! 	      elrbpdiff(inode)   = elrbpdiff(inode)   + p2vecdiff(1,inode,igaus) * (gprhs(1,igaus)-gptemp(1,igaus)) &
! 		    &                          + p2vecdiff(2,inode,igaus) * (gprhs(2,igaus)-gptemp(2,igaus))
! 	end do
!     end do
!   else
!     do igaus=1,pgaus
! 	do inode=1,pnode
! 	      elrbudiff(1,inode) = elrbudiff(1,inode) + p1vecdiff(inode,igaus)   * gprhs(1,igaus)
! 	      elrbudiff(2,inode) = elrbudiff(2,inode) + p1vecdiff(inode,igaus)   * gprhs(2,igaus)
! 	      elrbudiff(3,inode) = elrbudiff(3,inode) + p1vecdiff(inode,igaus)   * gprhs(3,igaus)
! 	      elrbpdiff(inode)   = elrbpdiff(inode)   + p2vecdiff(1,inode,igaus) * (gprhh(1,igaus)-gptemp(1,igaus)) &
! 		    &                          + p2vecdiff(2,inode,igaus) * (gprhh(2,igaus)-gptemp(2,igaus)) &
! 		    &                          + p2vecdiff(3,inode,igaus) * (gprhh(3,igaus)-gptemp(3,igaus))
! 	end do
!     end do
!   end if


  !----------------------------------------------------------------------
  !
  ! elresudiff = elraudiff + elrbudiff
  ! elrespdiff = elrapdiff + elrbpdiff
  !
  !----------------------------------------------------------------------
  
  do inode = 1,pnode
    elrespdiff(inode) = elrapdiff(inode) - elrbpdiff(inode)
    do idime = 1,ndime
      elresudiff(idime,inode) = elraudiff(idime,inode) - elrbudiff(idime,inode)
    enddo
  enddo

	
	
end subroutine nsi_elmresdiff_visco

