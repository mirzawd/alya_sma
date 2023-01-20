!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_twonat(&
     ndime,pnode,nturb_tur,iunkn_tur,param_tur,&
     gptur,gpden,gpvis,gpmut,gpwal,gpgra,eledd,&
     gpcar,gppro,gprea,gpdif,gprhs,gpgrd,gpvel)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_twonat
  ! NAME
  !   tur_twonat
  ! DESCRIPTION
  !    Compute coefficient of the equation of Chien's k-eps model:
  !
  !       dK         dK       dui            d +-          dK  -+
  !    rho-- + rho*ui--- = tij--- - rho*e + ---|(mu+mut/sk)---  |
  !       dt         dxi      dxj           dxi+-          dxi -+
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
  integer(ip), intent(in)    :: ndime,pnode,nturb_tur,iunkn_tur
  real(rp),    intent(in)    :: param_tur(*),eledd(pnode)
  real(rp),    intent(in)    :: gptur(nturb_tur)
  real(rp),    intent(in)    :: gpgra,gpden,gpvis,gpmut,gpwal,gppro
  real(rp),    intent(in)    :: gpcar(ndime,pnode)
  real(rp),    intent(inout) :: gpdif,gpgrd(ndime),gprea,gprhs
  real(rp),    intent(inout) :: gpvel(ndime)
  real(rp)                   :: gpkin,gpeps
  real(rp)                   :: leps,vv,ystarv,ystar,nu
  real(rp)                   :: gpdif_inn,gpgrd_inn(ndime),gprea_inn
  real(rp)                   :: gprhs_inn,gpvel_inn(ndime)
  real(rp)                   :: wemix,xeps
  integer(ip)                :: idime,inode
  !
  ! Definitions
  !
  gpkin = gptur(1)                                          ! K
  gpeps = gptur(2)                                          ! eps
  nu    = gpvis/gpden                                       ! nu
  !
  ! XEPS = sqrt(VV)/LEPS
  !
  ystar  = gpwal*sqrt(gpkin)/nu                             ! y* = y*K^{1/2}/nu
  vv     = 7.19e-3_rp*ystar-4.33e-5_rp*ystar**2_rp+8.8e-8_rp*ystar**3_rp
  vv     = gpkin*vv
  ystarv = gpwal*sqrt(vv)/nu                                ! (yv*) = y*vv^{1/2}/nu
  if(ystarv<1.0e-10_rp) then
     xeps  = nu/(0.88_rp*gpwal*gpwal)                       ! sqrt(vv)/leps -> nu/(0.88*y^2)
  else
     leps  = 8.8_rp*gpwal/(1.0_rp+10.0_rp/ystarv&           ! leps = 8.8*y/(1+10/(yv*)+0.0515*(yv*)) 
          &  +0.0515_rp*ystarv)
     xeps  = sqrt(vv)/leps
  end if
  !
  ! GPRHS, GPDIF, GPGRD, GPVEL and GPREA
  !
  if(iunkn_tur==1) then
     gprhs_inn = gppro + gpgra                              ! P+G
     gprea_inn = gpden*xeps                                 ! rho*sqrt(vv)/leps
     gpdif_inn=gpvis+gpmut/param_tur(1)                     ! mu+mu_t/sk
     do idime=1,ndime                                       ! grad(mu_t)/sk
        gpgrd_inn(idime)=0.0_rp
     end do
     do inode=1,pnode
        do idime=1,ndime
           gpgrd_inn(idime)=gpgrd_inn(idime)&
                +gpcar(idime,inode)*eledd(inode)
        end do
     end do
     do idime=1,ndime 
        gpgrd_inn(idime)=gpgrd_inn(idime)/param_tur(1)
     end do
     gpvel_inn = gpvel
  else
     gprea_inn = 1.0_rp
     gprhs_inn = gprea_inn*xeps*gpkin                        ! eps=K*sqrt(vv)/leps
     gpdif_inn = 1.0e-06_rp*gprea_inn
     gpgrd_inn = 0.0_rp
     gpvel_inn = 0.0_rp
  end if
  !
  ! Mix the inner and outer layers
  !
  call mixing(1_ip,wemix,param_tur(8),param_tur(9),ystar)    ! Mixing function
  gpdif = wemix*gpdif+(1.0_rp-wemix)*gpdif_inn
  gprea = wemix*gprea+(1.0_rp-wemix)*gprea_inn
  gprhs = wemix*gprhs+(1.0_rp-wemix)*gprhs_inn
  gpgrd = wemix*gpgrd+(1.0_rp-wemix)*gpgrd_inn
  gpvel = wemix*gpvel+(1.0_rp-wemix)*gpvel_inn

end subroutine tur_twonat
