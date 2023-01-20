!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmprc(&
     itask,ndime,ndofn,pnode,pevat,pgaus,lnods,gpsha,gpvol,&
     gpden,dtinv_nsi,gpst1,gpst2,elmat,elpat)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmprc
  ! NAME 
  !    nsi_elmprc
  ! DESCRIPTION
  !    This routine compute the elemental preconditioner
  !    (tau*U,V) approx (L^{-1}*U,V) where tau=(tau1*I,tau2)
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp 
  use def_nastin, only     :  kfl_fixno_nsi,pabdf_nsi
  implicit none
  integer(ip), intent(in)  :: itask,ndime,pnode,pevat
  integer(ip), intent(in)  :: pgaus,ndofn
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(in)  :: gpst1(pgaus),gpst2(pgaus)
  real(rp),    intent(in)  :: gpvol(pgaus),gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: gpden(pgaus),dtinv_nsi
  real(rp),    intent(in)  :: elmat(pevat,pevat)
  real(rp),    intent(out) :: elpat(pevat,pevat)
  integer(ip)              :: idime,inode,jnode,igaus,ipoin
  integer(ip)              :: idofv,jdofv,idofp,jdofp,ievat,jevat
  real(rp)                 :: fact1,fact2,adiag
  !
  ! Initialization
  !
  do ievat=1,pevat  
     do jevat=1,pevat
        elpat(jevat,ievat)=0.0_rp
     end do
  end do
  !
  ! Momentum equation
  !
  if(itask==2) then
     do igaus=1,pgaus
        fact2=gpvol(igaus)/(gpden(igaus)*( dtinv_nsi*pabdf_nsi(1) )+1.0_rp/gpst1(igaus))
        do inode=1,pnode
           idofv=(inode-1)*ndofn
           fact1=fact2*gpsha(inode,igaus)
           do idime=1,ndime
              idofv=idofv+1
              do jnode=1,pnode
                 jdofv=(jnode-1)*ndofn+idime
                 elpat(idofv,jdofv)=elpat(idofv,jdofv)&
                      +fact1*gpsha(jnode,igaus)
              end do
           end do
        end do
     end do
     !do inode=1,pevat
     !   elpat(inode,inode)=elmat(inode,inode)
     !end do
     !
     ! Boundary conditions
     !
     do inode=1,pnode
        idofv=(inode-1)*ndofn
        ipoin=lnods(inode)
        do idime=1,ndime
           idofv=idofv+1
           if(      kfl_fixno_nsi(idime,ipoin)==1&
                .or.kfl_fixno_nsi(idime,ipoin)==9) then
              adiag=elpat(idofv,idofv)
              do jdofv=1,pevat
                 elpat(idofv,jdofv)=0.0_rp
                 elpat(jdofv,idofv)=0.0_rp
              end do
              elpat(idofv,idofv)=adiag
           end if
        end do
     end do
  end if

  return
  !
  ! Continuity equation
  !  
  do igaus=1,pgaus
     do inode=1,pnode
        fact1=gpvol(igaus)*gpst2(igaus)*gpsha(inode,igaus)
        idofp=inode*ndofn
        do jnode=1,pnode
           jdofp=jnode*ndofn
           elpat(jdofp,idofp)=elpat(jdofp,idofp)&
                +fact1*gpsha(jnode,igaus)
        end do
     end do
  end do

end subroutine nsi_elmprc
