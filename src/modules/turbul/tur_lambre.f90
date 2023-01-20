!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_lambre(&
     ndime,pnode,gptur,gpden,gpvis,gpmut,gpwal,&
     gpgra,eledd,gpcar,gppro,gppr2,gprea,gpdif,&
     gprhs,gpgrd)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_lambre
  ! NAME
  !   tur_lambre
  ! DESCRIPTION
  !    Compute coefficient of the equation of  Lam-Bremhorst k-e model
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
  use def_kintyp, only     :  ip,rp
  use def_turbul, only     :  iunkn_tur,param_tur,kfl_algor_tur
  implicit none
  integer(ip), intent(in)  :: ndime,pnode
  real(rp),    intent(in)  :: eledd(pnode),gptur(2)
  real(rp),    intent(in)  :: gpgra,gpden,gpvis,gpmut,gpwal,gppro,gppr2
  real(rp),    intent(in)  :: gpcar(ndime,pnode)
  real(rp),    intent(out) :: gpdif(kfl_algor_tur)
  real(rp),    intent(out) :: gpgrd(kfl_algor_tur,ndime)
  real(rp),    intent(out) :: gprea(kfl_algor_tur,kfl_algor_tur)
  real(rp),    intent(out) :: gprhs(kfl_algor_tur)
  real(rp)                 :: gpkin,gpeps,ReT,Ry,onReT
  real(rp)                 :: f2,f1,fmu,ff,Ce1,Ce2,Cmu,epsok
  !
  ! Definitions
  !
  gpkin = gptur(1)
  gpeps = gptur(2)
  Ce1   = param_tur(3)                                    ! Ce1
  Ce2   = param_tur(4)                                    ! Ce2
  Cmu   = param_tur(6)                                    ! Cmu
  Ry    = gpden*sqrt(gpkin)*gpwal/gpvis                   ! Ry=k^{1/2}*y/nu
  !
  ! F2
  !
  if(gpeps==0.0_rp) then
     ReT = 0.0_rp                                         
     f2  = 1.0_rp
  else
     ReT = gpden*gpkin*gpkin/(gpvis*gpeps)                ! ReT=k^2/(eps*nu)
     if(ReT>100.0_rp) then
        f2  = 1.0_rp
     else if(ReT<0.01_rp) then
        f2  = 0.0_rp     
     else
        f2  = 1.0_rp-exp(-ReT*ReT)
     end if
  end if
  !
  ! FMU, F1
  !
  if(gpkin==0.0_rp) then
     fmu   = 0.0_rp
     ff    = 0.0_rp
  else
     onReT = gpvis*gpeps/(gpden*gpkin*gpkin)
     fmu   = ((1.0_rp-exp(-0.0165_rp*Ry))**2.0_rp)*(1.0_rp+20.5_rp*onReT)
     f1    = 1.0_rp+(0.05_rp/fmu)**3.0_rp
     ff    = (fmu+0.05_rp**3.0_rp/(fmu*fmu))
  end if
  !
  ! Force and reaction GPRHS and GPREA
  !
  if(kfl_algor_tur==1) then
     !
     ! Decoupled
     !
     if(iunkn_tur==1) then  
        gprea(1,1)    = Cmu*gpden*gpden*fmu*gpkin/gpmut
        gprhs(1)      = gppro + gpgra 
        !gprea(1,1)    = Cmu*gpden*gpden*fmu*gpkin/gpmut-Cmu*gpden*fmu*gpkin/gpeps
        !gprhs(1)      = gpgra 
     else 
        if(gpkin==0.0_rp) then
           gprea(1,1) = 0.0_rp
           gprhs(1)   = 0.0_rp
        else
           epsok      = gpeps/gpkin
           gprea(1,1) = gpden*Ce2*f2*epsok                                 ! rho*Ce2*f2*eps/k
           gprhs(1)   = gpden*Ce1*Cmu*ff*gpkin*gppr2 - param_tur(5)*gpgra  ! Ce1*f1*rho*eps/k*P
        end if
     end if

  else
     !
     ! Coupled
     !
     gprea(1,1) = Cmu*gpden*gpden*fmu*gpkin/gpmut
     gprhs(1)   = gppro + gpgra + gpden*gpeps
     gprea(1,2) = gpden
     if(gpkin==0.0_rp) then
        gprea(2,2) = 0.0_rp
        gprhs(2)   = 0.0_rp
     else
        epsok      = gpeps/gpkin
        gprea(2,2) = gpden*Ce2*f2*epsok                                 ! rho*Ce2*f2*eps/k
        gprhs(2)   = gpden*Ce1*Cmu*ff*gpkin*gppr2 - param_tur(5)*gpgra  ! Ce1*f1*rho*eps/k*P
     end if

  end if
  !
  ! GPDIF, GPGRD: diffusion and its gradient
  !
  call tur_elmdif(pnode,gpvis,gpmut,gpcar,eledd,gpdif,gpgrd)

end subroutine tur_lambre
