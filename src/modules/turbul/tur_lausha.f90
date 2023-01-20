!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_lausha(&
     pnode,gptur,gpden,gpvis,gpmut,gpgra,gphev,&
     eledd,gpcar,gppro,gpsqk,gprea,gpdif,gprhs,&
     gpgrd)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_lausha
  ! NAME
  !   tur_lausha
  ! DESCRIPTION
  !    Compute coefficient of the equation of Launder-Sharma-eps model:
  !    B.E. Launder and B.I. Sharma. Applications of the 
  !    Energy-Dissipation model of turbulence to the calculation of
  !    flow near a spinning disc. Letters in Heat and Mass Transf.
  !    Vol 1 (2), 131-138 (1974).
  !    rho*Dk/Dt-div[(mu+mut/sk) grad(k)] = 
  !               P + G - rho*e - D
  !    rho*De/Dt-div[(mu+mut/se) grad(e)] = 
  !              (Ce1*f1*P + Ce3*G - rho*Ce2*f2*e)*e/k + E
  !
  !    mut = fmu*Cmu*rho*k^2/e
  !    ReT = k^2/(nu*eps)
  !
  !    sk  = 1.00
  !    se  = 1.30
  !    Ce1 = 1.44
  !    Ce2 = 1.92
  !    Ce3 = 0.80
  !    Cmu = 0.09
  !
  !    D   = 2*mu*(grad(sqrt(k))^2
  !    E   = 2*nu*mut*[ d^2 ui/(dxn dxm) ]^2
  !    G   = beta*mut/Prt*[g.grad(T)]
  !    P   = mu_t*(dui/dxj+duj/dxi)*dui/dxj
  !    f1  = 1
  !    f2  = 1-0.3*exp(-ReT^2) 
  !    fmu = exp[-3.4/(1+ReT/50)^2]
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
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
  use def_turbul, only       :  iunkn_tur,kfl_algor_tur,param_tur
  implicit none
  integer(ip), intent(in)    :: pnode
  real(rp),    intent(in)    :: eledd(pnode)
  real(rp),    intent(in)    :: gptur(*)
  real(rp),    intent(in)    :: gpgra,gpden,gpvis,gpmut,gppro
  real(rp),    intent(in)    :: gphev,gpsqk
  real(rp),    intent(in)    :: gpcar(ndime,pnode)
  real(rp),    intent(out)   :: gpdif(kfl_algor_tur)
  real(rp),    intent(out)   :: gpgrd(kfl_algor_tur,ndime)
  real(rp),    intent(inout) :: gprea(kfl_algor_tur,kfl_algor_tur)
  real(rp),    intent(out)   :: gprhs(kfl_algor_tur)
  real(rp)                   :: gpkin,gpeps,ReT,f2,f1,E,D,fmu
  real(rp)                   :: Ce1,Ce2,Ce3,epsok,Cmu
  !
  ! Definitions
  !
  gpkin = gptur(1)
  gpeps = gptur(2)
  Ce1   = param_tur(3)
  Ce2   = param_tur(4)
  Ce3   = param_tur(5)
  Cmu   = param_tur(6)
  !
  ! Functions
  !
  if(gpeps/=0.0_rp) then
     ReT = gpden*gpkin*gpkin/(gpvis*gpeps)                     ! k^2/(eps*nu)
     fmu = (1.0_rp+0.02_rp*ReT)
     fmu = fmu*fmu
     fmu = exp(-3.4_rp/fmu)
     f2  = 1.0_rp-0.3_rp*exp(-ReT*ReT)     
  else
     f2  = 1.0_rp
     fmu = 1.0_rp
  end if
  f1 = 1.0_rp
  D  = 2.0_rp*gpvis*gpsqk                                      ! 2*mu*(grad(sqrt(k))^2
  E  = 2.0_rp*(gpvis/gpden)*gpmut*gphev**2                     ! E=2*nu*mut*[ d^2 ui/(dxn dxm) ]^2
  ! 
  ! Force and reaction GPRHS and GPREA
  ! 
  if(kfl_algor_tur==1) then
     !
     ! Decoupled
     !
     if(iunkn_tur==1) then
        gprhs(1)   = gppro + gpgra - D                         ! P + G - 2*nu*grad[sqrt(k)]^2
        if(gpkin/=0.0_rp) then
           gprea(1,1) = gpden*gpden*fmu*Cmu/gpmut*gpkin        ! rho*eps <= rho*[fmu*Cmu*rho*k/mut]*k
        else
           gprhs(1)   = gprhs(1) - gpden*gpeps 
        end if
     else
        if(gpkin/=0.0_rp) then
           epsok      = gpeps/gpkin
           gprhs(1)   = (Ce1*f1*gppro + Ce3*gpgra)*epsok + E   ! (Ce1*f1*P + Ce3*G)*eps/k + E
           gprea(1,1) = Ce2*f2*gpden*epsok                     ! rho*Ce2*f2*eps/k
        end if
     end if

  else
     !
     ! Coupled
     !
     gprhs(1)   = gppro + gpgra - D                            ! P + G - 2*nu*grad[sqrt(k)]^2
     gprea(1,2) = gpden                                        ! rho
     if(gpkin/=0.0_rp) then
        epsok      = gpeps/gpkin
        gprhs(2)   = (Ce1*f1*gppro + Ce3*gpgra)*epsok + E      ! (Ce1*f1*P + Ce3*G)*eps/k + E
        gprea(2,2) = Ce2*f2*gpden*epsok                        ! rho*Ce2*f2*eps/k
     else
        gprhs(2)   = E
     end if
  end if
  !
  ! GPDIF, GPGRD: diffusion and its gradient
  !
  call tur_elmdif(pnode,gpvis,gpmut,gpcar,eledd,gpdif,gpgrd)

end subroutine tur_lausha
