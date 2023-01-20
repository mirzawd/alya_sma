!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_bouopb(&
     lboel,gbsha,gbcar,baloc,wmatr,pnode,&
     pnodb,ndofn,pevat,gbvis,shapp)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_bouopb
  ! NAME 
  !    nsi_bouopb
  ! DESCRIPTION
  !    This routine computes the contribution to WMATR for the Navier-
  !    Stokes equations due to open boundaries. In this case, the term
  !    S.n (S = -p I + 2 mu Sym(grad(u)) being the Cauchy stress tensor) 
  !    is considered unknown.
  !           +-
  !    LHS =  |  sig.n v ds => 
  !          -+
  !
  !           +-                      +-
  !    LHS -  |  mu*grad(u).n v ds =  |  -p.n ds
  !          -+                      -+
  !    or
  !           +-                              +-
  !    LHS -  |  2 mu Sym(grad(u).n.n v ds =  |  -p.n ds, or
  !          -+                              -+
  ! USES
  ! USED BY
  !    nsi_bouope 
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
  use def_nastin, only       :  fvins_nsi
  implicit none
  integer(ip), intent(in)    :: pnode,pnodb,pevat,ndofn
  integer(ip), intent(in)    :: lboel(pnodb)
  real(rp),    intent(in)    :: gbsha(pnodb),shapp(pnode)
  real(rp),    intent(in)    :: gbcar(ndime,pnode),baloc(ndime)
  real(rp),    intent(inout) :: wmatr(pevat,pevat)
  integer(ip)                :: inodb,ldime,ievab,jnode,jevab,kdime
  integer(ip)                :: jnods,jevat,idime
  real(rp)                   :: gbvis,xmuit,prod1
  !
  ! Contribution from the viscous term
  !
  do inodb = 1,pnodb
     xmuit = gbsha(inodb) * gbvis
     do ldime = 1,ndime
        ievab = (lboel(inodb)-1) * ndofn + ldime
        do jnode = 1,pnode
           jevab = (jnode-1) * ndofn + ldime
           prod1 = 0.0_rp
           do idime = 1,ndime
              prod1 = prod1 + gbcar(idime,jnode) * baloc(idime)
           end do
           wmatr(ievab,jevab) = wmatr(ievab,jevab) - xmuit * prod1
        end do
        if( fvins_nsi > 0.9_rp ) then
           do jnode = 1,pnode
              do kdime = 1,ndime
                 jevab = (jnode-1) * ndofn + kdime
                 wmatr(ievab,jevab) = wmatr(ievab,jevab) &
                      -xmuit * gbcar(ldime,jnode) * baloc(kdime)
              end do
           end do
        end if
        if( fvins_nsi > 1.9_rp ) then
           do jnode = 1,pnode
              do kdime = 1,ndime
                 jevab = (jnode-1) * ndofn + kdime
                 wmatr(ievab,jevab) = wmatr(ievab,jevab) &
                      + 2.0_rp / 3.0_rp * xmuit * gbcar(kdime,jnode) * baloc(ldime)
              end do
           end do
        end if
     end do
  end do

  return

  !
  ! Contribution from the pressure term     
  !     
  do ldime=1,ndime
     prod1=baloc(ldime)
     do inodb=1,pnodb
        ievab=(lboel(inodb)-1)*ndofn+ldime
        xmuit=gbsha(inodb)*prod1
        do jnods=1,pnode
           jevat=jnods*ndofn
           wmatr(ievab,jevat)=wmatr(ievab,jevat)&
                +xmuit*shapp(jnods)
        end do
     end do
  end do

end subroutine nsi_bouopb

subroutine nsi_bouopb3(&
     lboel,gbsha,gbcar,baloc,wmatr,pnode,&
     pnodb,ndofn,pevat,gbvis,shapp,elcod)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_bouopb
  ! NAME 
  !    nsi_bouopb
  ! DESCRIPTION
  !    This routine computes the contribution to WMATR for the Navier-
  !    Stokes equations due to open boundaries. In this case, the term
  !    S.n (S = -p I + 2 mu Sym(grad(u)) being the Cauchy stress tensor) 
  !    is considered unknown.
  ! USES
  ! USED BY
  !    nsi_bouope 
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
  implicit none
  integer(ip), intent(in)    :: pnode,pnodb,pevat,ndofn
  integer(ip), intent(in)    :: lboel(pnodb)
  real(rp),    intent(in)    :: gbsha(pnodb),shapp(pnode)
  real(rp),    intent(in)    :: gbcar(ndime,pnode),baloc(ndime),elcod(ndime,*)
  real(rp),    intent(inout) :: wmatr(pevat,pevat)
  integer(ip)                :: idime,inode
  real(rp)                   :: gbvis,xmuit,elfun(pnode)
  !
  ! Contribution from the viscous term: 2 mu Sym(grad(u).n
  !
  do inode = 1,pnode
     elfun(inode) = 2.0_rp * elcod(1,inode) + 3.0_rp * elcod(2,inode) 
  end do
!print*,elcod(1:ndime,1:pnode)
!print*,gbcar(1:ndime,1:pnode)
  xmuit = 0.0_rp
  do idime = 1,ndime
     xmuit = 0.0_rp
     do inode = 1,pnode
        xmuit = xmuit + gbcar(idime,inode)*elfun(inode)
     end do
     print*,'XX=',xmuit
  end do
stop


end subroutine nsi_bouopb3
