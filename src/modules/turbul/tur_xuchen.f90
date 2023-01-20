!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_xuchen(&
     ndime,pnode,param_tur,gptur,gpden,gpvis,&
     gpmut,gpwal,gpgra,eledd,gpcar,gppro,gprea,&
     gpdif,gprhs,gpgrd)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_xuchen
  ! NAME
  !   tur_xuchen
  ! DESCRIPTION
  !    Compute coefficient of the equation of Xu-Chen-Nieuwstadt k model:
  !
  !    rho*dK/dt + rho*u.grad(K) = tau:grad(u) + div[(mu+mut/sk)*grad(K)]
  !                                -rho*eps + beta*mut/Prt*g.grad(T)
  !    where
  !
  !    eps  = sqrt(vv)*k/leps
  !    mut  = rho*sqrt(vv)*lmu
  !    tau  = -rho*ui'uj' = mut*(grad(u)+grad(u)^t)-2/3*delta_ij*k
  !
  !    with
  !
  !    leps = 8.8*y/[1+10/(yv*)+5.15e-2*(yv*)]
  !    lmu  = 0.544*y/(1+5.025e-04*(yv*)^1.65)
  !    vv/k = 7.19e-3(y*)-4.33e-5*(y*)^2+8.8e-08*(y*)^3
  !    y*   = y*sqrt(k)/nu
  !    yv*  = y*sqrt(vv)/nu
  !
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
  use      def_kintyp
  implicit none
  integer(ip), intent(in)    :: ndime,pnode
  real(rp),    intent(in)    :: param_tur(*),eledd(pnode)
  real(rp),    intent(in)    :: gptur
  real(rp),    intent(in)    :: gpgra,gpden,gpvis,gpmut,gpwal,gppro
  real(rp),    intent(in)    :: gpcar(ndime,pnode)
  real(rp),    intent(inout) :: gpdif,gpgrd(ndime),gprea,gprhs
  real(rp)                   :: gpkin,leps,vv,ystarv,ystar,nu,xeps
  !
  ! Definitions
  !
  gpkin = gptur                                             ! K
  nu    = gpvis/gpden                                       ! nu
  !
  ! XEPS = sqrt(VV)/LEPS
  !
  ystar  = gpwal*sqrt(gpkin)/nu                             ! y* = y*K^{1/2}/nu
  vv     = 7.19e-03_rp*ystar-4.33e-05_rp*ystar**2_rp+8.8e-08_rp*ystar**3_rp
  vv     = gpkin*vv
  ystarv = gpwal*sqrt(vv)/nu                                ! (yv*) = y*vv^{1/2}/nu
  if(ystarv<1.0e-10_rp) then
     xeps  = nu/(0.88_rp*gpwal*gpwal)                       ! vv^{1/2}/leps=nu/(0.88*y^2)
  else
     leps  = 8.8_rp*gpwal/(1.0_rp+10.0_rp/ystarv&           ! leps = 8.8*y/(1+10/(yv*)+0.0515*(yv*)) 
          &  +0.0515_rp*ystarv)
     xeps  = sqrt(vv)/leps
  end if
  !
  ! GPRHS, GPDIF, GPGRD, GPVEL and GPREA
  !
  gprhs = gppro + gpgra                                     ! P+G
  gprea = gpden*xeps                                        ! rho*sqrt(vv)/leps
  gpdif = gpvis+gpmut/param_tur(1)                          ! mu+mu_t/sk

end subroutine tur_xuchen
