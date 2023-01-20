!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_kchien(&
     ndime,pnode,gptur,gpden,gpvis,gpmut,gpwal,&
     gpgra,eledd,gpust,gpcar,gppro,gppr2,gprea,&
     gpdif,gprhs,gpgrd)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_kchien
  ! NAME
  !   tur_kchien
  ! DESCRIPTION
  !    Compute coefficient of the equation of Chien's k-eps model:
  !    K.-Y. Chien
  !    Predictions of channel and boundary-layer flows with a low-Reynolds
  !    number turbulence model
  !    AIAA J., 20 (1), 33-38 (1982).
  !
  !    Compute coefficient of the equation of Chien's k-eps model:
  !
  !       dK         dK       dui                 d +-          dK  -+
  !    rho-- + rho*ui--- = tij--- - rho*(e+e0) + ---|(mu+mut/sk)---  | + G
  !       dt         dxi      dxj                dxi+-          dxi -+
  !
  !       de         de               e    dui              e^2    
  !    rho-- + rho*ui--- = rho*Ce1*f1*- tij--- - rho*Ce2*f2 --- + rho*E 
  !       dt         dxi              K    dxj               K    
  !                            e      d +-          de  -+
  !                      + ce3*-*G + ---|(mu+mut/se)---  |
  !                            K     dxi+-          dxi -+
  !
  !    fmu  = 1-exp(-0.0115*y+)
  !    f1   = 1
  !    f2   = 1-0.22*exp[-(ReT/6)**2]
  !    e0   = 2*nu*K/y^2
  !    E    = -2*nu*eps/y^2*exp(-y+/2)
  !    G    = beta*mut/Prt*g.grad(T)
  !
  !    Ce1=1.35, Ce2=1.80, Cmu=0.09, sk=1.0, se=1.3
  !
  !    The equation is assembled as follows:
  !
  !       dK         dK        mu                K^2      dui   
  !    rho-- + rho*ui--- + 2*K --- + rho*Cmu*fmu*--- = tij--- + G
  !       dt         dxi       y^2               nut      dxj 
  !                                               d +-          dK  -+
  !                                           +  ---|(mu+mut/sk)---  |
  !                                              dxi+-          dxi -+
  !
  !       de         de               e^2          e                          
  !    rho-- + rho*ui--- + rho*Ce2*f2 --- + 2*mu* --- exp(-y+/2)  
  !       dt         dxi               K          y^2                               
  !                    e              e    dui    d +-          de  -+
  !           =  Ce3*G*- + rho*Ce1*f1*- tij--- + ---|(mu+mut/se)---  |
  !                    K              K    dxj   dxi+-          dxi -+
  !
  !    The seond term of the RHS of K-equation is substituted by
  !    rho*e=rho*Cmu*fmu*k/nut
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
  real(rp)                   :: gpkin,gpeps,ReT,gpypl,epsok,dummr
  real(rp)                   :: f2,f1,fmu,Cmu,Ce1,Ce2,Ce3
  !
  ! Definitions
  !
  gpypl = gpwal*gpust/(gpvis/gpden)                         ! y+  = y*U*/nu
  gpkin = gptur(1)                                          ! k
  gpeps = gptur(2)                                          ! eps
  Ce1   = param_tur(3)                                      ! Ce1
  Ce2   = param_tur(4)                                      ! Ce2
  Ce3   = param_tur(5)                                      ! Ce3
  Cmu   = param_tur(6)                                      ! Cmu
  fmu   = 1.0_rp-exp(-0.0115_rp*gpypl)                      ! fmu = 1-exp(-0.0115*y+)
  f1    = 1.0_rp                                            ! f1  = 1
  if(gpeps==0.0_rp) then
     ReT = 0.0_rp
  else
     ReT = gpkin*gpkin/(gpeps*gpvis/gpden)                  ! ReT = k^2/(eps*nu)
  end if
  if(ReT>=30.0_rp) then
     f2 = 1.0_rp
  else
     f2 = 1.0_rp-0.22_rp*exp(-(ReT/6.0_rp)**2.0_rp)         ! f2  = 1-0.22*exp[-(ReT/6)^2]
  end if
  !
  ! Force and reaction GPRHS and GPREA
  !
  if(kfl_algor_tur==1) then
     !
     ! Decoupled
     !
     if(iunkn_tur==1) then
        gprea(1,1) = 2.0_rp*gpvis/(gpwal*gpwal)               ! 2*mu/(y^2)
        gprhs(1)   = gprhs(1) + gppro + gpgra - gpden*gpeps   ! -rho*eps

        !gprhs(1)   = gppro + gpgra                           ! P + G 
        !gprea(1,1) = 2.0_rp*gpvis/(gpwal*gpwal)              ! 2*nu*grad[sqrt(k)]^2
        !if(gpkin/=0.0_rp) then
        !   gprea(1,1) = gprea(1,1) &       
        !        &       + gpden*gpden*fmu*Cmu/gpmut*gpkin    ! rho*eps <= rho*[fmu*Cmu*rho*k/mut]*k
        !else

        !   gprhs(1)   = gprhs(1) - gpden*gpeps 
        !end if

     else 
        if(gpkin==0.0_rp) then
            epsok  =  0.0_rp
        else
           epsok   =  gpeps/gpkin
        end if
        gprea(1,1) =  gpden*Ce2*f2*epsok&
             &        +2.0_rp*gpvis/(gpwal*gpwal)*exp(-0.5_rp*gpypl)
        gprhs(1)   =  gprhs(1) + Ce3*epsok*gpgra + Ce1*f1*gpden*Cmu*fmu*gpkin*gppr2
     end if

  else
     !
     ! Coupled
     !
     gprhs(1) = gprhs(1) + gpgra
     if(gpeps>1.0e-10_rp) then
        dummr      = 2.0_rp*gpvis/(gpwal*gpwal)-gpden*Cmu*fmu*gpkin/gpeps*gppr2
        if(dummr>0.0_rp) then
           gprea(1,1) = gprea(1,1) + dummr
        else
           gprea(1,1) = 2.0_rp*gpvis/(gpwal*gpwal)
           gprhs(1)   = gprhs(1) + gppro
        end if
     else
        gprhs(1)   = gprhs(1) + gppro                       ! P+G
        gprea(1,1) = 2.0_rp*gpvis/(gpwal*gpwal)             ! 2*mu/(y^2)
     end if
     gprea(1,2) = gpden                                     ! rho
     if(gpkin==0.0_rp) then
        epsok     =  0.0_rp
     else
        epsok      =  gpeps/gpkin
     end if
     gprea(2,2) = gpden*Ce2*f2*epsok&
          &       +2.0_rp*gpvis/(gpwal*gpwal)*exp(-0.5_rp*gpypl)
     gprhs(2)   = gprhs(2) + Ce3*epsok*gpgra + Ce1*f1*gpden*Cmu*fmu*gpkin*gppr2     
  end if
  !
  ! GPDIF, GPGRD: diffusion and its gradient
  !
  call tur_elmdif(pnode,gpvis,gpmut,gpcar,eledd,gpdif,gpgrd)

end subroutine tur_kchien
