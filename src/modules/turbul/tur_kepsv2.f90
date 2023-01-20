!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_kepsv2(&
     ndime,pnode,gptur,gpden,gpvis,gpmut,gpwal,&
     gpgra,eledd,eltur,gpcar,gppro,gppr2,gplap,&
     gprea,gpdif,gprhs,gpgrd)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_kepsv2
  ! NAME
  !   tur_kepsv2
  ! DESCRIPTION
  !    Compute coefficient of the equation of Chien's k-eps model:
  ! OUTPUT 
  !    GPREA .......... r 
  !    GPDIF .......... k 
  !    GPRHS .......... f 
  !    GPGRD(NDIME) ... grad(k) coefficient
  ! USES
  ! USED BY
  !    tur_elmcoe
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_master, only       :  ittim
  use def_turbul, only       :  nturb_tur,iunkn_tur,&
       &                        param_tur,inits_tur
  implicit none
  integer(ip), intent(in)    :: ndime,pnode
  real(rp),    intent(in)    :: eledd(pnode)
  real(rp),    intent(in)    :: eltur(nturb_tur,pnode)
  real(rp),    intent(in)    :: gptur(nturb_tur)
  real(rp),    intent(in)    :: gpgra,gpden,gpvis,gpmut,gpwal,gppro
  real(rp),    intent(in)    :: gppr2,gplap
  real(rp),    intent(in)    :: gpcar(ndime,pnode)
  real(rp),    intent(out)   :: gpdif
  real(rp),    intent(out)   :: gpgrd(ndime)
  real(rp),    intent(inout) :: gprea
  real(rp),    intent(out)   :: gprhs
  integer(ip)                :: idime,inode
  real(rp)                   :: gpkin,gpeps,gpphi,gpfff
  real(rp)                   :: Ce1,Ce2,Ce3,nu,g,T
  real(rp)                   :: CL,Cn,c1,c2,L,f,oneT
  real(rp)                   :: oneL2,sk,onReT,fmu,Ry,maxfg
  real(rp)                   :: grk(3),grp(3),grkgrp,Cmu
  !
  ! Definitions
  !
  nu    = gpvis/gpden                                       ! nu
  gpkin = gptur(1)                                          ! k
  gpeps = gptur(2)                                          ! eps
  gpphi = gptur(3)                                          ! phi (v2/k)
  gpfff = gptur(4)                                          ! f
  Ce2   =   1.90_rp
  Ce3   =   0.80_rp
  c1    =   1.40_rp
  c2    =   0.30_rp
  CL    =   0.25_rp
  Cn    = 110.00_rp
  Ce1   =   1.40_rp*(1.0_rp+0.05_rp*sqrt(1.0_rp/gpphi))
  Cmu   = param_tur(6)
  sk    = param_tur(1)
  !
  ! ONET = 1/T
  !
  if(iunkn_tur==2.or.iunkn_tur==4) then
     f    = gpkin/gpeps
     g    = 6.0_rp*sqrt(nu/gpeps)
     call mixmax(4.0_rp,f,g,T)
     oneT = 1.0_rp/T
  end if
  !
  ! GRKGRP = grad(k)*grad(phi)
  !
  if(iunkn_tur==3.or.iunkn_tur==4) then     
     grk = 0.0_rp
     grp = 0.0_rp
     do inode=1,pnode
        do idime=1,ndime
           grk(idime)=grk(idime)+gpcar(idime,inode)*eltur(1,inode)
           grp(idime)=grp(idime)+gpcar(idime,inode)*eltur(3,inode)
        end do
     end do
     grkgrp=0.0_rp
     do idime=1,ndime
        grkgrp=grkgrp+grk(idime)*grp(idime)    
     end do
  end if
  !
  ! Force and reaction GPRHS and GPREA
  !
  if(iunkn_tur==1) then
      
     if(ittim>=inits_tur) then
        f     = gpkin
        g     = 6.0_rp*sqrt(nu*gpeps)
        call mixmax(4.0_rp,f,g,maxfg)                  ! = max (k, 6*(nu*eps)^{1/2} )
        gprea = gpden*gpden*Cmu*gpphi/gpmut*maxfg
        gprhs = gppro + gpgra
     else
        Ry    = gpden*sqrt(gpkin)*gpwal/gpvis          ! Ry=k^{1/2}*y/nu
        onReT = gpvis*gpeps/(gpden*gpkin*gpkin)
        fmu   = ((1.0_rp-exp(-0.0165_rp*Ry))**2.0_rp)*(1.0_rp+20.5_rp*onReT)
        fmu   = 1.0_rp
        gprea = 0.09_rp*gpden*gpden*fmu*gpkin/gpmut
        gprhs = gppro + gpgra 
     end if 
 
  else if(iunkn_tur==2) then
     gprhs    = Ce1*gppro*oneT + Ce3*gpgra*gpeps/gpkin
     gprea    = Ce2*gpden*oneT

  else if(iunkn_tur==3) then     
     gprhs    = gpden*gpfff + 2.0_rp/(sk*gpkin)*gpmut*grkgrp 
     gprea    = gppro/gpkin 

  else if(iunkn_tur==4) then
     f       = gpkin**1.5_rp/gpeps
     g       = Cn*(nu*nu*nu/gpeps)**0.25_rp
     call mixmax(4.0_rp,f,g,L)
     L       = CL*L
     oneL2   = 1.0_rp/(L*L)
     gprhs   = oneL2*( &
          &         - oneT*(c1-1.0_rp)*(gpphi-2.0_rp/3.0_rp) &
          &         + c2*gppro/(gpden*gpkin) &
          &         + 2.0_rp*nu/gpkin*grkgrp + nu*gplap )                   
     gprea   = oneL2
  end if
  !
  ! GPDIF, GPGRD: diffusion and its gradient
  !
  if(iunkn_tur==4) then   
     gpdif = 1.0_rp
     gpgrd = 0.0_rp
  else
     call tur_elmdif(pnode,gpvis,gpmut,gpcar,eledd,gpdif,gpgrd)
  end if

end subroutine tur_kepsv2
