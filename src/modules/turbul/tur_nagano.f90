!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_nagano(&
     ndime,pnode,gptur,gpden,gpvis,gpmut,gpwal,&
     gpgra,eledd,gpust,gpcar,gppro,gppr2,gprea,&
     gpdif,gprhs,gpgrd)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_nagano
  ! NAME
  !   tur_nagano
  ! DESCRIPTION
  !    Compute coefficient of the equation of Jaw-Hwang's k-eps model:
  !    Y. Nagano and M. Tagawa
  !    An improved k-epsilon model for boundary layer flows
  !    J. Fluid Eng., 112, 33-39 (1990).
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
  use def_turbul, only       :  nturb_tur,iunkn_tur,kfl_algor_tur,&
       &                        param_tur
  implicit none
  integer(ip), intent(in)    :: ndime,pnode
  real(rp),    intent(in)    :: eledd(pnode)
  real(rp),    intent(in)    :: gpust,gptur(nturb_tur)
  real(rp),    intent(in)    :: gpgra,gpden,gpvis,gpmut,gpwal,gppro
  real(rp),    intent(in)    :: gppr2,gpcar(ndime,pnode)
  real(rp),    intent(out)   :: gpdif(kfl_algor_tur)
  real(rp),    intent(out)   :: gpgrd(kfl_algor_tur,ndime)
  real(rp),    intent(inout) :: gprea(kfl_algor_tur,kfl_algor_tur)
  real(rp),    intent(out)   :: gprhs(kfl_algor_tur)
  real(rp)                   :: gpkin,gpeps,f1,f2,Ce1,Ce2,Ce3
  real(rp)                   :: epsok,gpypl,Cmu,fmu,nu,ReT
  !
  ! Definitions 
  !
  gpkin = gptur(1)                                          ! k
  gpeps = gptur(2)                                          ! eps
  nu    = gpvis/gpden                                       ! nu
  gpypl = gpwal*gpust/nu
  Ce1   = param_tur(3)
  Ce2   = param_tur(4)
  Ce3   = param_tur(5)
  Cmu   = param_tur(6)
  ReT   = gpkin*gpkin/(nu*gpeps)
  !
  ! Force and reaction GPRHS and GPREA
  !
  if(kfl_algor_tur==1) then
     !
     ! Decoupled
     !
     if(iunkn_tur==1) then
        fmu        =  (1.0_rp-exp(-gpypl/26.0_rp))**2&
             &       *(1.0_rp+4.1_rp/ReT**0.75_rp)
        gprhs(1)   = gppro + gpgra
        gprea(1,1) = gpden*gpden*fmu*Cmu/gpmut*gpkin    ! eps = [fmu*Cmu*rho*k/mut]*k
     else 
        if(gpkin/=0.0_rp) then
           f1         =  1.0_rp
           f2         =  (1.0_rp-0.3_rp*exp(-(ReT/6.5_rp)**2))&
                &       *(1.0_rp-exp(-gpypl/6.0_rp))**2
           epsok      = gpeps/gpkin
           gprhs(1)   = Ce1*f1*gppro*epsok + Ce3*gpgra*epsok
           gprea(1,1) = Ce2*f2*gpden*epsok
        end if
     end if

  else
     !
     ! Coupled
     !
     fmu        =  (1.0_rp-exp(-gpypl/26.0_rp))**2&
          &       *(1.0_rp+4.1_rp/ReT**0.75_rp)
     gprhs(1)   = gppro + gpgra
     gprea(1,2) = gpden 
     if(gpkin/=0.0_rp) then
        f1         =  1.0_rp
        f2         =  (1.0_rp-0.3_rp*exp(-(ReT/6.5_rp)**2))&
             &       *(1.0_rp-exp(-gpypl/6.0_rp))**2
        epsok      = gpeps/gpkin
        gprhs(2)   = Ce1*f1*gppro*epsok + Ce3*gpgra*epsok
        gprea(2,2) = Ce2*f2*gpden*epsok
     end if

  end if
  !
  ! GPDIF, GPGRD: diffusion and its gradient
  !
  call tur_elmdif(pnode,gpvis,gpmut,gpcar,eledd,gpdif,gpgrd)

end subroutine tur_nagano
