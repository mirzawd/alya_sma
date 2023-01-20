!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmcon(&
     pnode,pevat,pgaus,ndofn,elvel,gppre,gpden,gprhs,&
     gptem,gpgrt,rmomu,rcont,gpsha,gpcar,p2vec,p2sca,wmatr,&
     wrhsi,wmatc) 
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmcon
  ! NAME 
  !    nsi_elmcon
  ! DESCRIPTION
  !    Compute the continuity equation at Gauss points
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime,mnode
  use def_nastin, only     :  kfl_penal_nsi,&
       &                      kfl_savco_nsi,penal_nsi,&
       &                      kfl_regim_nsi
  use def_kermod, only     :  gasco
  implicit none
  integer(ip), intent(in)  :: pnode,pevat,pgaus,ndofn
  real(rp),    intent(in)  :: gppre(pgaus),elvel(ndime,*)
  real(rp),    intent(in)  :: gpden(pgaus),gprhs(ndime+1,pgaus)
  real(rp),    intent(in)  :: gptem(pgaus),gpgrt(ndime,pgaus)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)  :: p2vec(ndime,pnode,pgaus)
  real(rp),    intent(in)  :: p2sca(pnode,pgaus)
  real(rp),    intent(in)  :: rmomu(ndime,ndime,pnode,pgaus)
  real(rp),    intent(in)  :: rcont(ndime+1,pnode,pgaus)
  real(rp),    intent(out) :: wmatr(pevat,pevat,pgaus)
  real(rp),    intent(out) :: wmatc(pevat,pevat,pgaus)
  real(rp),    intent(out) :: wrhsi(pevat,pgaus)
  integer(ip)              :: idime,inode,jnode,jdime,igaus
  integer(ip)              :: idofv,jdofv,idofp,jdofp
  real(rp)                 :: fact1,fact2,fact3
  !
  ! Velocity: ( div(u) , q )
  !
  if(kfl_regim_nsi/=1) then
     if(kfl_savco_nsi==1) then
        do igaus=1,pgaus
           do inode=1,pnode
              idofv=(inode-1)*ndofn
              do idime=1,ndime
                 idofv=idofv+1
                 fact1=rcont(idime,inode,igaus)
                 do jnode=1,pnode
                    jdofp=jnode*ndofn
                    wmatc(jdofp,idofv,igaus)=wmatc(jdofp,idofv,igaus)&
                         +fact1*p2sca(jnode,igaus)
                 end do
              end do
           end do
        end do
     else if(kfl_savco_nsi/=2) then
        if(ndime==2) then
           do igaus=1,pgaus
              do inode=1,pnode
                 idofv=(inode-1)*ndofn

                 idofv=idofv+1
                 fact1=rcont(1,inode,igaus)
                 do jnode=1,pnode
                    jdofp=jnode*ndofn
                    wmatr(jdofp,idofv,igaus)=wmatr(jdofp,idofv,igaus)&
                         +fact1*p2sca(jnode,igaus)
                 end do

                 idofv=idofv+1
                 fact1=rcont(2,inode,igaus)
                 do jnode=1,pnode
                    jdofp=jnode*ndofn
                    wmatr(jdofp,idofv,igaus)=wmatr(jdofp,idofv,igaus)&
                         +fact1*p2sca(jnode,igaus)
                 end do
              end do
           end do
        else
           do igaus=1,pgaus
              do inode=1,pnode
                 idofv=(inode-1)*ndofn

                 idofv=idofv+1
                 fact1=rcont(1,inode,igaus)
                 do jnode=1,pnode
                    jdofp=jnode*ndofn
                    wmatr(jdofp,idofv,igaus)=wmatr(jdofp,idofv,igaus)&
                         +fact1*p2sca(jnode,igaus)
                 end do

                 idofv=idofv+1
                 fact1=rcont(2,inode,igaus)
                 do jnode=1,pnode
                    jdofp=jnode*ndofn
                    wmatr(jdofp,idofv,igaus)=wmatr(jdofp,idofv,igaus)&
                         +fact1*p2sca(jnode,igaus)
                 end do

                 idofv=idofv+1
                 fact1=rcont(3,inode,igaus)
                 do jnode=1,pnode
                    jdofp=jnode*ndofn
                    wmatr(jdofp,idofv,igaus)=wmatr(jdofp,idofv,igaus)&
                         +fact1*p2sca(jnode,igaus)
                 end do
              end do
           end do
        end if
     end if
  end if

  if(kfl_regim_nsi==2) then
     !
     ! Density: ( R* [ T*grad(rho)+rho*grad(T) ], tau1' grad(q) ); P2VEC is symmmetric
     !  
     if(ndime==2) then
        do igaus=1,pgaus
           fact2=gasco*gptem(igaus)
           do inode=1,pnode
              idofp=inode*ndofn
              do jnode=inode+1,pnode
                 jdofp=jnode*ndofn
                 fact1=p2vec(1,jnode,igaus)*gpcar(1,inode,igaus)&
                      +p2vec(2,jnode,igaus)*gpcar(2,inode,igaus)
                 fact1=fact1*fact2
                 wmatr(jdofp,idofp,igaus)=wmatr(jdofp,idofp,igaus)+fact1
                 wmatr(idofp,jdofp,igaus)=wmatr(idofp,jdofp,igaus)+fact1
              end do
              fact1=p2vec(1,inode,igaus)*gpcar(1,inode,igaus)&
                   +p2vec(2,inode,igaus)*gpcar(2,inode,igaus)
              wmatr(idofp,idofp,igaus)=wmatr(idofp,idofp,igaus)+fact1*fact2
              fact3=p2vec(1,inode,igaus)*gpgrt(1,igaus)&
                   +p2vec(2,inode,igaus)*gpgrt(2,igaus)
              fact3=fact3*gasco
              do jnode=1,pnode
                 jdofp=jnode*ndofn                  
                 wmatr(idofp,jdofp,igaus)=wmatr(idofp,jdofp,igaus)&
                      +fact3*gpsha(jnode,igaus)
              end do
           end do
        end do
     else
        do igaus=1,pgaus
           fact2=gasco*gptem(igaus)
           do inode=1,pnode
              idofp=inode*ndofn
              do jnode=inode+1,pnode
                 jdofp=jnode*ndofn
                 fact1=p2vec(1,jnode,igaus)*gpcar(1,inode,igaus)&
                      +p2vec(2,jnode,igaus)*gpcar(2,inode,igaus)&
                      +p2vec(3,jnode,igaus)*gpcar(3,inode,igaus)
                 fact1=fact1*fact2
                 wmatr(jdofp,idofp,igaus)=wmatr(jdofp,idofp,igaus)+fact1
                 wmatr(idofp,jdofp,igaus)=wmatr(idofp,jdofp,igaus)+fact1
              end do
              fact1=p2vec(1,inode,igaus)*gpcar(1,inode,igaus)&
                   +p2vec(2,inode,igaus)*gpcar(2,inode,igaus)&
                   +p2vec(3,inode,igaus)*gpcar(3,inode,igaus)
              wmatr(idofp,idofp,igaus)=wmatr(idofp,idofp,igaus)+fact1*fact2
              fact3=p2vec(1,inode,igaus)*gpgrt(1,igaus)&
                   +p2vec(2,inode,igaus)*gpgrt(2,igaus)
              fact3=fact3*gasco
              do jnode=1,pnode
                 jdofp=jnode*ndofn                  
                 wmatr(idofp,jdofp,igaus)=wmatr(idofp,jdofp,igaus)&
                      +fact3*gpsha(jnode,igaus)
              end do
           end do
        end do
     end if

  else
     !
     ! Pressure: ( grad(p) , tau1' grad(q) ); P2VEC is symmmetric
     ! 
     if(ndime==2) then
        do igaus=1,pgaus
           do inode=1,pnode
              idofp=inode*ndofn
              do jnode=inode+1,pnode
                 jdofp=jnode*ndofn
                 fact1=p2vec(1,jnode,igaus)*gpcar(1,inode,igaus)&
                      +p2vec(2,jnode,igaus)*gpcar(2,inode,igaus)
                 wmatr(jdofp,idofp,igaus)=wmatr(jdofp,idofp,igaus)+fact1
                 wmatr(idofp,jdofp,igaus)=wmatr(idofp,jdofp,igaus)+fact1
              end do
              fact1=p2vec(1,inode,igaus)*gpcar(1,inode,igaus)&
                   +p2vec(2,inode,igaus)*gpcar(2,inode,igaus)
              wmatr(idofp,idofp,igaus)=wmatr(idofp,idofp,igaus)+fact1
           end do
        end do
     else
        do igaus=1,pgaus
           do inode=1,pnode
              idofp=inode*ndofn
              do jnode=inode+1,pnode
                 jdofp=jnode*ndofn
                 fact1=p2vec(1,jnode,igaus)*gpcar(1,inode,igaus)&
                      +p2vec(2,jnode,igaus)*gpcar(2,inode,igaus)&
                      +p2vec(3,jnode,igaus)*gpcar(3,inode,igaus)
                 wmatr(jdofp,idofp,igaus)=wmatr(jdofp,idofp,igaus)+fact1
                 wmatr(idofp,jdofp,igaus)=wmatr(idofp,jdofp,igaus)+fact1
              end do
              fact1=p2vec(1,inode,igaus)*gpcar(1,inode,igaus)&
                   +p2vec(2,inode,igaus)*gpcar(2,inode,igaus)&
                   +p2vec(3,inode,igaus)*gpcar(3,inode,igaus)
              wmatr(idofp,idofp,igaus)=wmatr(idofp,idofp,igaus)+fact1
           end do
        end do
     end if
  end if
  !
  ! Velocity: ( rho*(uc.grad)u + 2*rho*(w x u) + sig*u -div[2*mu*eps(u)], tau1' grad(q) )
  !
  if(ndime==2) then
     do igaus=1,pgaus
        do inode=1,pnode
           idofp=inode*ndofn
           do jnode=1,pnode
              jdofv=(jnode-1)*ndofn
              do jdime=1,ndime
                 jdofv=jdofv+1
                 wmatr(idofp,jdofv,igaus)=wmatr(idofp,jdofv,igaus)&
                      +p2vec(1,inode,igaus)&
                      *rmomu(1,jdime,jnode,igaus)
                 wmatr(idofp,jdofv,igaus)=wmatr(idofp,jdofv,igaus)&
                      +p2vec(2,inode,igaus)&
                      *rmomu(2,jdime,jnode,igaus)
              end do
           end do
        end do
     end do
  else
     do igaus=1,pgaus
        do inode=1,pnode
           idofp=inode*ndofn
           do jnode=1,pnode
              jdofv=(jnode-1)*ndofn
              do jdime=1,ndime
                 jdofv=jdofv+1
                 wmatr(idofp,jdofv,igaus)=wmatr(idofp,jdofv,igaus)&
                      +p2vec(1,inode,igaus)&
                      *rmomu(1,jdime,jnode,igaus)
                 wmatr(idofp,jdofv,igaus)=wmatr(idofp,jdofv,igaus)&
                      +p2vec(2,inode,igaus)&
                      *rmomu(2,jdime,jnode,igaus)
                 wmatr(idofp,jdofv,igaus)=wmatr(idofp,jdofv,igaus)&
                      +p2vec(3,inode,igaus)&
                      *rmomu(3,jdime,jnode,igaus)
              end do
           end do
        end do
     end do
  end if
  !
  ! Right-hand side: ( rho*f , tau1' grad(q) )
  !
  if(ndime==2) then
     do igaus=1,pgaus  
        do inode=1,pnode
           idofp=inode*ndofn
           wrhsi(idofp,igaus)=wrhsi(idofp,igaus)&
                +gprhs(1,igaus)*p2vec(1,inode,igaus)
           wrhsi(idofp,igaus)=wrhsi(idofp,igaus)&
                +gprhs(2,igaus)*p2vec(2,inode,igaus)
        end do
     end do
  else
     do igaus=1,pgaus  
        do inode=1,pnode
           idofp=inode*ndofn
           wrhsi(idofp,igaus)=wrhsi(idofp,igaus)&
                +gprhs(1,igaus)*p2vec(1,inode,igaus)
           wrhsi(idofp,igaus)=wrhsi(idofp,igaus)&
                +gprhs(2,igaus)*p2vec(2,inode,igaus)
           wrhsi(idofp,igaus)=wrhsi(idofp,igaus)&
                +gprhs(3,igaus)*p2vec(3,inode,igaus)
        end do
     end do
  end if
  !
  ! Penalization: ( eps*p , q )
  !
  if(kfl_penal_nsi>=1) then
     do igaus=1,pgaus
        do inode=1,pnode             
           idofp=inode*ndofn
           fact1=penal_nsi*gpsha(inode,igaus)
           do jnode=1,pnode
              jdofp=jnode*ndofn
              wmatr(idofp,jdofp,igaus)=wmatr(idofp,jdofp,igaus)&
                   +fact1*gpsha(jnode,igaus)
           end do
        end do
     end do
     if(kfl_penal_nsi/=2) then             ! Term ( eps*p^{i-1} , q )
        do igaus=1,pgaus                   ! or   ( eps*p^{n-1} , q )
           fact1=penal_nsi*gppre(igaus)
           do inode=1,pnode                
              idofp=inode*ndofn
              wrhsi(idofp,igaus)=wrhsi(idofp,igaus)&  
                   +fact1*gpsha(inode,igaus)
           end do
        end do
     end if
  end if

  if(kfl_regim_nsi==1.or.kfl_regim_nsi==2) then
     !
     ! Left-hand side:  ( 1/(RT*theta*dt)*p + u.grad(p)   , v )
     !                  ( 1/(theta*dt)*rho  + u.grad(rho) + div(u)*rho , v )
     !
     do igaus=1,pgaus
        do inode=1,pnode
           idofp=inode*ndofn
           do jnode=1,pnode
              jdofp=inode*ndofn
              wmatr(idofp,jdofp,igaus)=wmatr(idofp,jdofp,igaus)&
                   +p2sca(inode,igaus)*rcont(ndime+1,jnode,igaus)
           end do
        end do
     end do
  end if

  if(kfl_regim_nsi>=1) then
     !
     ! Right-hand side /= 0 for non-incompressible models
     !
     do igaus=1,pgaus
        do inode=1,pnode
           idofp=inode*ndofn
           wrhsi(idofp,igaus)=wrhsi(idofp,igaus)&
                +p2sca(inode,igaus)*gprhs(ndime+1,igaus)
        end do
     end do

  end if

end subroutine nsi_elmcon
