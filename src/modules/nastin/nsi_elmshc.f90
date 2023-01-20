!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmshc(&
     pnode,pgaus,ptopo,pevat,ndofn,gpden,gpvel,gprhs,gpsp1,&
     gpvol,elvel,gpcar,chale,rmomu,rmom2,elauu)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmshc
  ! NAME
  !   nsi_elmshc
  ! DESCRIPTION
  !   This routine computes the contribution from the shock capturing
  ! OUTPUT
  !    ELAUU ... Element matrix
  ! USES
  ! USED BY
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mgaus,mnode
  use def_nastin, only       :  kfl_shock_nsi,shock_nsi,kfl_rmom2_nsi
  implicit none  
  integer(ip), intent(in)    :: pnode,pgaus,ptopo,pevat,ndofn
  real(rp),    intent(in)    :: gpvel(ndime,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: elvel(ndime,pnode)
  real(rp),    intent(in)    :: gpden(pgaus),gprhs(ndime,pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus)
  real(rp),    intent(in)    :: chale,gpsp1(pgaus)
  real(rp),    intent(in)    :: rmomu(pnode,pgaus)
  real(rp),    intent(in)    :: rmom2(ndime,ndime,pnode,pgaus)
  real(rp),    intent(inout) :: elauu(pevat,pevat)                      
  real(rp)                   :: scdif,sldif,gpve2,grnor(3),vepar,facta
  real(rp)                   :: addit,gpgrd(ndime,ndime),fact1,factt
  real(rp)                   :: gpres(mgaus)
  integer(ip)                :: idime,inode,jnode,kdime,ldime,igaus
  integer(ip)                :: idofn,jdofn,jdime
  !
  ! Factor according to element topology
  !
  if( ptopo == 0 ) then
     factt = 1.5_rp             ! Quadrilateral/Hexahedra
  else if( ptopo == 1 ) then
     factt = 0.7_rp             ! Triangles/Tetrahedra
  else
     factt = 1.0_rp             ! Other elements
  end if
  !
  ! Isotropic/anisotropic shock capturing
  !
  if( kfl_shock_nsi == 1 ) then  
     facta = 0.0_rp             ! Isotropic SC
  else
     facta = 1.0_rp             ! Anisotropic SC
  end if
  !
  ! Coarse grid residual |R(u)|
  ! |R(u)| = |f - rho*u/(theta*dt)-rho*a.grad(u)+div[k*grad(u)]-s*u|
  ! f = Q + rho*u^n/(theta*dt)
  !
  do igaus = 1,pgaus
     gpres(igaus) = 0.0_rp
     do inode = 1,pnode
        do idime = 1,ndime
           gpres(igaus) = gpres(igaus)&
                - rmomu(inode,igaus) * elvel(idime,inode)              
        end do
     end do 
  end do
  if( kfl_rmom2_nsi /= 0 ) then
    do igaus = 1,pgaus
        do idime = 1,ndime
           do inode = 1,pnode
              do jdime = 1,ndime
                 gpres(igaus) = gpres(igaus)&
                      - rmom2(idime,jdime,inode,igaus) * elvel(jdime,inode)
              end do
           end do
        end do
     end do
  end if
  do igaus = 1,pgaus
     gpres(igaus) = abs(gpres(igaus))
  end do
  !
  ! Compute the contribution to the element stiffness matrix
  ! 
  do igaus = 1,pgaus
     !
     ! Square velocity norm
     !
     gpve2 = 0.0_rp
     do idime = 1,ndime
        gpve2 = gpve2 + gpvel(idime,igaus) * gpvel(idime,igaus)
     end do
     !
     ! Velocity gradient
     !
     do jdime = 1,ndime
        do idime = 1,ndime
           gpgrd(idime,jdime) = 0.0_rp
        end do
     end do
     do inode = 1,pnode
        do jdime = 1,ndime
           do idime=1,ndime
              gpgrd(idime,jdime)=gpgrd(idime,jdime)&
                   +gpcar(idime,inode,igaus)*elvel(jdime,inode)
           end do
        end do
     end do
     do idime = 1,ndime
        call vecnor(gpgrd(1,idime),ndime,grnor(idime),2_ip)
     end do

     if(grnor(1)>1.0e-10_rp) then
        vepar = gpres(igaus) / grnor(1)
        scdif = factt * 0.5_rp * shock_nsi * chale * vepar
        if( scdif > 0.0_rp ) then
           !
           ! Compute diffusion introduced along the streamlines
           !
           sldif = max(0.0_rp,scdif-facta*gpsp1(igaus)*gpden(igaus)*gpve2)
           if( gpve2 > 1.0e-10_rp ) then
              fact1 = (sldif-scdif) / gpve2
           else
              fact1 = 0.0_rp
           end if
           do inode = 1,pnode
              idofn = inode*ndofn
              do jnode = 1,pnode 
                 jdofn = jnode*ndofn
                 addit = 0.0_rp
                 do kdime = 1,ndime
                    addit = addit + scdif * gpcar(kdime,inode,igaus)&   ! kiso*grad(Ni).grad(Nj)
                         * gpcar(kdime,jnode,igaus)
                    do ldime = 1,ndime                                  ! (<kiso-k'>-kiso)/u^2*
                       addit = addit + fact1 &                          ! grad(Ni).(u x u).grad(Nj)
                            * gpcar(kdime,inode,igaus)&
                            * gpvel(kdime,igaus)&
                            * gpvel(ldime,igaus)&
                            * gpcar(ldime,jnode,igaus)
                    end do
                 end do
                 elauu(idofn,jdofn) = elauu(idofn,jdofn) &
                      + addit * gpvol(igaus)
              end do
           end do
        end if
     end if
  end do

end subroutine nsi_elmshc
!------------------------------------------------------------------------
! NOTES
!
! Shock capturing for the advection-diffusion-reaction equation
! 
!     du
! rho -- + rho*a.grad(u) - div[k*grad(u)] + s*u = f
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
!            +          2k   +           R
! ac   = max | 0 , C - ----- | , a*= ----------- grad(u), therefore
!            +         |a*|h +       |grad(u)|^2
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
