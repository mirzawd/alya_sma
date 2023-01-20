!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_tworod(&
     ndime,pnode,nturb_tur,iunkn_tur,param_tur,&
     gptur,gpden,gpvis,gpmut,gpwal,gpgra,eledd,&
     gpust,gpsha,gpcar,gppro,gprea,gpdif,gprhs,&
     gpgrd,gpvel)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_tworod
  ! NAME
  !   tur_tworod
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
  real(rp),    intent(in)    :: gpust,gpsha(pnode),gptur(nturb_tur)
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
  ! leps
  !
  ystar  = gpwal*sqrt(gpkin)/nu                             ! y*  = y*K^{1/2}/nu
  vv     = gpkin*0.4_rp*(1.0_rp-exp(-ystar**2/4200.0_rp))   ! vv  = 0.4*k*[1-exp(-y*^2/4200)]
  ystarv = gpwal*sqrt(vv)/nu                                ! y*v = y*vv^{1/2}/nu
  if(ystarv<1.0e-10_rp) then
     xeps  = 2.12_rp*nu/(1.3_rp*gpwal*gpwal)                ! vv^{1/2}/leps=2.12*nu/(1.3*y^2)
  else
     leps  = 1.3_rp*gpwal/(1.0_rp+2.12_rp/ystarv&           ! leps = 1.3*y/(1+2.12/y*v+0.028*yv*) 
          &  +0.028_rp*ystarv)
     xeps  = sqrt(vv)/leps
  end if
  !
  ! GPRHS, GPDIF, GPGRD, GPVEL and GPREA
  !
  if(iunkn_tur==1) then
     gprhs_inn = gppro + gpgra                              ! P+G
     gprea_inn = gpden*xeps                                 ! rho*vv^{1/2}/leps
     gpdif_inn=gpvis+gpmut/param_tur(iunkn_tur)             ! mu+mu_t/sk
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
        gpgrd_inn(idime)=gpgrd_inn(idime)&
             /param_tur(iunkn_tur)
     end do
     gpvel_inn = gpvel
  else
     gprea_inn = 1.0_rp
     gprhs_inn = gprea_inn*xeps*gpkin                        ! eps=K*vv^{1/2}/leps
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

end subroutine tur_tworod
