!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_elmshc(&
     pnode,plapl,pgaus,ptopo,gpden,gpvel,gpdif,gprea,&
     gprhs,gpgrd,gptur,gpsta,gpvol,eltur,gphes,gpcar,&
     hleng,elmat)
  !-----------------------------------------------------------------------
  !
  ! This routine computes the contribution from the shock capturing
  ! term to the turbulence equation      
  !
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,ntens,mnode
  use def_turbul, only       :  nturb_tur,iunkn_tur,kfl_shock_tur,&
       &                        kfl_timei_tur,shock_tur,dtinv_tur,&
       &                        nbdfp_tur,pabdf_tur
  implicit none  
  integer(ip), intent(in)    :: pnode,plapl,pgaus,ptopo
  real(rp),    intent(in)    :: gpvel(ndime,pgaus)
  real(rp),    intent(in)    :: gphes(ntens,mnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: eltur(nturb_tur,pnode) 
  real(rp),    intent(in)    :: gprea(pgaus),gpden(pgaus),gpgrd(ndime,pgaus)                     
  real(rp),    intent(in)    :: gprhs(pgaus),gptur(nturb_tur,nbdfp_tur+1,pgaus)
  real(rp),    intent(in)    :: hleng,gpsta(pgaus),gpdif(pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus)
  real(rp),    intent(inout) :: elmat(pnode,pnode)                             
  real(rp)                   :: advec,diffu,react,timed,scdif,sldif
  real(rp)                   :: resid,gpve2,grnor,vepar,fact1
  real(rp)                   :: addit,gpgrt(3),factt,facta,gplap
  integer(ip)                :: idime,inode,jnode,kdime,ldime,igaus,itime
  !
  ! Factor according to element topology
  !
  if(ptopo==0) then
     factt=1.5_rp             ! Quadrilateral/Hexahedra
  else if(ptopo==1) then
     factt=0.7_rp             ! Triangles/Tetrahedra
  else
     factt=1.0_rp             ! Other elements
  end if
  !
  ! Isotropic/anisotropic shock capturing
  !
  if(kfl_shock_tur==1) then  
     facta=0.0_rp             ! Isotropic SC:
  else
     facta=1.0_rp             ! Anisotropic SC:
  end if

  do igaus=1,pgaus
     !
     ! Square velocity norm
     !
     gpve2=0.0_rp
     do idime=1,ndime
        gpve2=gpve2+gpvel(idime,igaus)*gpvel(idime,igaus)
     end do
     !
     ! Turbulent variable gradient
     !
     do idime=1,ndime
        gpgrt(idime)=0.0_rp
     end do
     do inode=1,pnode
        do idime=1,ndime
           gpgrt(idime)=gpgrt(idime)&
                +gpcar(idime,inode,igaus)*eltur(iunkn_tur,inode)
        end do
     end do
     call vecnor(gpgrt,ndime,grnor,2_ip)
     !
     ! Advection term
     !
     advec=0.0_rp
     do idime=1,ndime
        advec=advec+gpvel(idime,igaus)*gpgrt(idime)
     end do
     advec=gpden(igaus)*advec
     ! 
     ! Reaction term
     !
     react=gprea(igaus)*gptur(iunkn_tur,1,igaus)
     !
     ! Diffusion term
     !   
     diffu=0.0_rp
     do idime=1,ndime
        diffu=diffu+gpgrd(idime,igaus)*gpgrt(idime)
     end do
     if(plapl==1) then
        do inode=1,pnode
           gplap=0.0_rp
           do idime=1,ndime
              gplap=gplap+gphes(idime,inode,igaus)
           end do
           diffu=diffu+gplap*eltur(iunkn_tur,inode)
        end do
     end if
     !
     ! Time derivative term
     !  
     if(kfl_timei_tur==1) then
        timed = gptur(iunkn_tur,1,igaus) * pabdf_tur(1)
        do itime = 3,nbdfp_tur+1
           timed = timed + gptur(iunkn_tur,itime,igaus) * pabdf_tur(itime-1)
        end do
        timed = gpden(igaus) * dtinv_tur * timed
     else
        timed=0.0_rp
     end if
     !
     ! Residual=rho*(u-u^n-1)/(theta*dt)+rho*a.grad(u)-div(k*grad(u))+r*u-f
     !  
     resid=abs(timed+advec-diffu+react-gprhs(igaus))

     if(grnor>1.0d-10) then
        vepar = resid/grnor
        scdif = (0.5_rp*factt*shock_tur*hleng*vepar-gpdif(igaus))
        if(scdif>0.0_rp) then
           !
           ! Compute the diffusion introduced along the streamlines
           ! sldif = <kiso-k'> = max(0,kiso-(rho*tau)*rho*a^2) for anisotropic
           !
           !f=0.0_rp
           !g=scdif-facta*gpden(igaus)*gpsta(igaus)*gpden(igaus)*gpve2
           !call mixmax(4.0_rp,f,g,sldif)

           sldif=max(0.0_rp,scdif-facta*gpden(igaus)*gpsta(igaus)*gpden(igaus)*gpve2)
           !
           ! Compute the contribution to the element stiffness matrix
           ! 
           if(gpve2>1.0d-10) then
              fact1=(sldif-scdif)/gpve2
           else
              fact1=0.0_rp
           end if
           do inode=1,pnode
              do jnode=1,pnode 
                 addit=0.0_rp
                 do kdime=1,ndime
                    addit=addit+scdif*gpcar(kdime,inode,igaus)& ! kiso*grad(Ni).grad(Nj)
                         *gpcar(kdime,jnode,igaus)
                    do ldime=1,ndime                            ! (<kiso-k'>-kiso)/u^2*
                       addit=addit+fact1&                       ! grad(Ni).(u x u).grad(Nj)
                            *gpcar(kdime,inode,igaus)&
                            *gpvel(kdime,igaus)*gpvel(ldime,igaus)&
                            *gpcar(ldime,jnode,igaus)
                    end do
                 end do
                 elmat(inode,jnode)=elmat(inode,jnode)&
                      +gpvol(igaus)*addit
              end do
           end do
        end if
     end if
  end do

end subroutine tur_elmshc
!------------------------------------------------------------------------
! NOTES
!
! Shock capturing for the advection-diffusion-reaction equation
! 
!     du
! rho -- + rho*a.grad(u) - div[k*grad(u)] + r*u = f
!     dt
!
! k    = Diffusion                  [M/(L*T)]    
! R    = Residual of the equation   [M*U/(L^3*T)]
! C    = Shock capturing constant 
! tau  = Stabilization parameter: 
!        Its units are h/(rho*a)=   [L^3*T/M]
!        so that rho*tau is in [T]
! kiso = Isotropic SC diffusion     [M/(L*T)]
! k'   = Anisotropic SC diffusion   [M/(L*T)]
!
!        1    |R|    h  
! Pe   = - --------- - 
!        2 |grad(u)| k 
!
!            +            2k     +          R/rho
! ac   = max | 0 , C - --------- | , a*= ----------- grad(u), therefore
!            +         rho*|a*|h +       |grad(u)|^2
!
!            +           1   +
! ac   = max | 0 , C - ----- |
!            +          Pe   +
!        1          |R|
! kiso = - ac*h  --------- , k'=(rho*tau)*rho*a^2
!        2       |grad(u)|
!
! Isotropic:     kiso*grad(u).grad(v)  
!                                         a x a
! Anisotropic:   (<kiso-k'>-kiso)*grad(u).-----.grad(v)
!                                          a^2
!***
!------------------------------------------------------------------------
