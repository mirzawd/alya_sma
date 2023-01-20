!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmrhc(&
     ielem,pgaus,pnode,lnods,gpvol,gpsha,gpcar,elvel,&
     chale,gpden,gpvis,gppor,gpgpr,dtinv_loc,rhsid)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmrhc
  ! NAME
  !    nsi_elmrhc
  ! DESCRIPTION
  !    The correction is:
  !    rho/dt*[ grad(p)^{i+1} - grad(p)^{i} ]
  ! USES
  ! USED BY
  !    nsi_elmcor
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_elmtyp, only     :  ELEXT
  use def_domain, only     :  ndime,mnode,npoin,lelch
  use def_nastin, only     :  deltp_nsi,kfl_taush_nsi,&
       &                      staco_nsi,kfl_predi_nsi,&
       &                      kfl_advec_nsi,&
       &                      pabdf_nsi
  use mod_tauadr, only     :  tauadr
  implicit none
  integer(ip), intent(in)  :: ielem
  integer(ip), intent(in)  :: pgaus
  integer(ip), intent(in)  :: pnode
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(in)  :: gpvol(pgaus),gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)  :: elvel(ndime,pnode)
  real(rp),    intent(in)  :: chale(*)
  real(rp),    intent(in)  :: gpden(pgaus),gpvis(pgaus),gppor(pgaus)
  real(rp),    intent(in)  :: dtinv_loc
  real(rp) ,   intent(out) :: gpgpr(ndime,pgaus),rhsid(ndime*npoin)
  integer(ip)              :: igaus,idime,inode,ipoin,idofn,knode
  real(rp)                 :: fact1,fact2,fact3
  real(rp)                 :: gpadv(3),adv,dif,rea,tau,gpvno
  !
  ! Extension elements
  !
  if( lelch(ielem) == ELEXT ) then
     knode = 1
  else
     knode = pnode
  end if
  !
  ! Pressure gradients at gauss points
  !
  do igaus = 1,pgaus
     do idime = 1,ndime
        gpgpr(idime,igaus) = 0.0_rp
     end do
  end do
  do inode = 1,pnode
     ipoin = lnods(inode)
     fact1 = deltp_nsi(ipoin)
     do igaus = 1,pgaus
        do idime = 1,ndime
           gpgpr(idime,igaus) = gpgpr(idime,igaus)&
                +gpcar(idime,inode,igaus)*fact1
        end do
     end do
  end do
  !
  ! Assembly
  !
  if( kfl_predi_nsi == 2 ) then      ! kfl_predi_nsi :  2 DT ; 3 TAU ; 4 MASS ; 5 LUMPED
     !
     ! (dt/rho)*[ grad(p^i+1) - grad(p^i) ]*v
     !
     if( ndime == 2 ) then
        do igaus = 1,pgaus
           fact2 = gpvol(igaus) / ( gpden(igaus) * dtinv_loc * pabdf_nsi(1) )
           do inode = 1,knode
              fact3        = fact2*gpsha(inode,igaus)
              idofn        = 2*lnods(inode)-1
              rhsid(idofn) = rhsid(idofn)-fact3*gpgpr(1,igaus)
              idofn        = idofn+1
              rhsid(idofn) = rhsid(idofn)-fact3*gpgpr(2,igaus)
           end do
        end do
     else
        do igaus = 1,pgaus
           fact2 = gpvol(igaus) / ( gpden(igaus) * dtinv_loc * pabdf_nsi(1) )
           do inode = 1,knode
              fact3        = fact2*gpsha(inode,igaus)
              idofn        = 3*lnods(inode)-2
              rhsid(idofn) = rhsid(idofn)-fact3*gpgpr(1,igaus)
              idofn        = idofn+1
              rhsid(idofn) = rhsid(idofn)-fact3*gpgpr(2,igaus)
              idofn        = idofn+1
              rhsid(idofn) = rhsid(idofn)-fact3*gpgpr(3,igaus)
           end do
        end do
     end if

  else if( kfl_predi_nsi == 3 ) then
     !
     ! tau'1*[ grad(p^i+1) - grad(p^i) ]*v, tau'=1/( rho/dt + 1/tau )
     !
     do igaus = 1,pgaus

        gpadv = 0.0_rp
        if( kfl_advec_nsi /= 0 ) then
           do inode = 1,pnode
              gpadv(1:ndime) = gpadv(1:ndime) + gpsha(inode,igaus) * elvel(1:ndime,inode)
           end do
        end if
        call vecnor(gpadv,ndime,gpvno,2_ip)

        adv = gpvno * gpden(igaus)                   ! Convective term
        dif = gpvis(igaus)                           ! Viscous term
        rea = gppor(igaus)                           ! Reaction term
        call tauadr(&
             kfl_taush_nsi,staco_nsi,adv,dif,rea,&
             chale(1),chale(2),tau)

         if( tau /= 0.0_rp ) then
           fact2 = 1.0_rp / ( gpden(igaus) * dtinv_loc * pabdf_nsi(1) + 1.0_rp / tau )
        else
           fact2 = 1.0_rp / ( gpden(igaus) * dtinv_loc * pabdf_nsi(1) )
        end if

        fact2 = gpvol(igaus) * fact2

        if( ndime == 2 ) then
           do inode=1,knode
              fact3        = fact2*gpsha(inode,igaus)
              idofn        = 2*lnods(inode)-1
              rhsid(idofn) = rhsid(idofn)-fact3*gpgpr(1,igaus)
              idofn        = idofn+1
              rhsid(idofn) = rhsid(idofn)-fact3*gpgpr(2,igaus)
           end do
        else
           do inode=1,knode
              fact3        = fact2*gpsha(inode,igaus)
              idofn        = 3*lnods(inode)-2
              rhsid(idofn) = rhsid(idofn)-fact3*gpgpr(1,igaus)
              idofn        = idofn+1
              rhsid(idofn) = rhsid(idofn)-fact3*gpgpr(2,igaus)
              idofn        = idofn+1
              rhsid(idofn) = rhsid(idofn)-fact3*gpgpr(3,igaus)
           end do
        end if

     end do

  else if( kfl_predi_nsi == 4 ) then
     !
     ! M/mu*[ grad(p^i+1) - grad(p^i) ]*v
     !
     call runend('NSI_ELMRHC: MASS PRECONDTIONER, CORRECTION NOT IMPLEMENTED')   ! if this line is uncommented remember that
                                                                                 ! from vers 772 dtinv_nsi=1/dt , see if
                                                                                 ! pabdf_nsi(1) needds to be added
     if( ndime == 2 ) then
        do igaus=1,pgaus
           fact2=gpvol(igaus)/(gpden(igaus)*dtinv_loc)
           do inode=1,knode
              fact3        = fact2*gpsha(inode,igaus)
              idofn        = 2*lnods(inode)-1
              rhsid(idofn) = rhsid(idofn)-fact3*gpgpr(1,igaus)
              idofn        = idofn+1
              rhsid(idofn) = rhsid(idofn)-fact3*gpgpr(2,igaus)
           end do
        end do
     else
        do igaus=1,pgaus
           fact2=gpvol(igaus)/(gpden(igaus)*dtinv_loc)
           do inode=1,knode
              fact3        = fact2*gpsha(inode,igaus)
              idofn        = 3*lnods(inode)-2
              rhsid(idofn) = rhsid(idofn)-fact3*gpgpr(1,igaus)
              idofn        = idofn+1
              rhsid(idofn) = rhsid(idofn)-fact3*gpgpr(2,igaus)
              idofn        = idofn+1
              rhsid(idofn) = rhsid(idofn)-fact3*gpgpr(3,igaus)
           end do
        end do
     end if

  else if( kfl_predi_nsi == 7 ) then

     call runend('NSI_ELMRHC: NOT CODED')

  end if

end subroutine nsi_elmrhc
