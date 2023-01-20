!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_kjawhw(&
     ndime,pnode,gptur,gpden,gpvis,gpmut,gpwal,&
     gpgra,eledd,gpust,gpcar,gppro,gppr2,gprea,&
     gpdif,gprhs,gpgrd)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_kjawhw
  ! NAME
  !   tur_kjawhw
  ! DESCRIPTION
  !    Compute coefficient of the equation of Jaw-Hwang's k-eps model:
  !    S.-Y. Jaw and R.R. Hwang
  !    A two-scale low-Reynolds number turbulence model
  !    Int. J. Numer. Meth. Fluids, 33, 695-710 (2000).
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
  real(rp)                   :: oneT,epsok,H,f,g,fog
  real(rp)                   :: fmu,y,dummr,nu,Ry
  !
  ! Definitions
  !
  gpkin = gptur(1)                                          ! k
  gpeps = gptur(2)                                          ! eps
  nu    = gpvis/gpden                                       ! nu
  y     = gpwal
  Ce1   = param_tur(3)
  Ce2   = param_tur(4)
  Ce3   = param_tur(5)
  !
  ! Force and reaction GPRHS and GPREA
  !
  if(kfl_algor_tur==1) then
     !
     ! Decoupled
     !
     if(iunkn_tur==1) then
        !gprhs(1)   = gppro + gpgra - gpden*gpeps            ! P+G-rho*eps
        Ry         = sqrt(gpkin)*y/nu
        fmu        = (1.0_rp-exp(-Ry/70.0_rp))**1.75_rp
        f          = gpkin
        g          = 6.0_rp*sqrt(nu*gpeps)
        fog        = f/g
        call mixing(1_ip,H,4.0_rp,1.0_rp,fog)
        dummr      = f*H+(1.0_rp-H)*g
        gprhs(1)   = gppro + gpgra
        gprea(1,1) = gpden*fmu*param_tur(6)/gpmut*dummr     ! eps = [fmu*Cmu*max(k,6(nu*eps)^{1/2})/mut]*k

     else 
        if(gpkin/=0.0_rp) then
           epsok   = gpeps/gpkin
           f       = epsok
           g       = 1.0_rp/6.0_rp*sqrt(gpeps/nu)
           fog     = f/g
           call mixing(1_ip,H,4.0_rp,1.0_rp,fog)            ! T   = max (k/eps, 6*(nu/eps)^{1/2} )
           oneT    = g*H+(1.0_rp-H)*f                       ! 1/T = min (eps/k, 1/6*(eps/nu)^{1/2} )
        else
           epsok   = 0.0_rp 
           oneT    = 0.0_rp
        end if
        f1         = 1.0_rp+0.1_rp*gppro/(gpden*gpeps)      ! f1 = 1 + 0.1*(P/eps)
        f2         = 1.0_rp                                 ! f2 = 1
        gprhs(1)   = Ce1*f1*oneT*gppro + Ce3*gpgra*epsok
        gprea(1,1) = Ce2*f2*oneT*gpden
     end if

  else
     !
     ! Coupled
     !
     gprhs(1)   = gppro + gpgra
     gprea(1,2) = gpden
     if(gpkin/=0.0_rp) then
        epsok   = gpeps/gpkin
        f       = epsok
        g       = 1.0_rp/6.0_rp*sqrt(gpden*gpeps/gpvis)
        fog     = f/g
        call mixing(1_ip,H,4.0_rp,1.0_rp,fog) 
        oneT    = g*H+(1.0_rp-H)*f
     else
        epsok   = 0.0_rp 
        oneT    = 0.0_rp
     end if
     f1         = 1.0_rp+0.1_rp*gppro/(gpden*gpeps)
     f2         = 1.0_rp
     gprhs(2)   = Ce1*f1*oneT*gppro + Ce3*gpgra*epsok
     gprea(2,2) = Ce2*f2*oneT*gpden
  end if
  !
  ! GPDIF, GPGRD: diffusion and its gradient
  !
  call tur_elmdif(pnode,gpvis,gpmut,gpcar,eledd,gpdif,gpgrd)

end subroutine tur_kjawhw
