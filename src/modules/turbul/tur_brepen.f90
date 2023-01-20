!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_brepen(&
     ndime,pnode,gptur,gpden,gpvis,gpmut,gpgra,&
     eledd,eltur,gpcar,gppro,gppr2,gprea,gpdif,&
     gprhs,gpgrd)
  !------------------------------------------------------------------------
  !****f* Turbul/tur_brepen
  ! NAME
  !   tur_brepen
  ! DESCRIPTION
  !    Compute coefficient of the equation of k-w model
  !    J. Bredberg, S.-H. Peng and L. Davidson
  !    An improved k-w turbulence model applied to recirculating flows
  !    Int. J. Heat and Fluid Flow 23 (2002) 731-743
  ! OUTPUT 
  !    GPREA .......... r 
  !    GPDIF .......... k  
  !    GPRHS .......... f 
  !    GPGRD(NDIME) ... grad(k) coefficient
  ! USES
  ! USED BY
  !    tur_elmcoe
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_turbul, only     :  nturb_tur,iunkn_tur,param_tur
  implicit none
  integer(ip), intent(in)  :: ndime,pnode
  real(rp),    intent(in)  :: eledd(pnode),gptur(nturb_tur)
  real(rp),    intent(in)  :: eltur(nturb_tur,pnode)
  real(rp),    intent(in)  :: gpgra,gpden,gpvis,gpmut,gppro,gppr2
  real(rp),    intent(in)  :: gpcar(ndime,pnode)
  real(rp),    intent(out) :: gpdif(1),gpgrd(ndime),gprea(1)
  real(rp),    intent(out) :: gprhs(1)
  real(rp)                 :: gpkin,gpome,xprod,gpgrk(3),gpgrw(3)
  real(rp)                 :: ReT,fmu,grads
  integer(ip)              :: idime,inode
  !
  ! Definitions
  !
  gpkin=gptur(1)
  gpome=gptur(2)
  !
  ! Force and reaction GPRHS and GPREA
  !
  if(iunkn_tur==1) then
     !
     ! K-equation
     !
     gprhs(1) = gpgra
     gprea(1) = 0.0_rp

     if(gpome==0.0_rp .or. gpkin==0.0_rp ) then
        gprhs(1) = gprhs(1) + gppro               
     else

        ReT   = gpden*gpkin/(gpvis*gpome)
        if(ReT<0.1_rp) then
           fmu   = 0.09_rp
        else if(ReT>100.0_rp) then
           fmu   = 1.0_rp
        else
           fmu   = 0.09_rp + (0.91_rp + 1.0_rp/ReT**3_rp)&
                *(1.0_rp-exp(-(ReT/25.0_rp)**2.75_rp))
        end if        
        xprod = gpden*(param_tur( 6)*fmu*gppr2/gpome-param_tur( 3)*gpome)   
        if(xprod<0.0_rp) then
           gprea = gprea - xprod
        else
           gprhs(1) = gprhs(1) + xprod*gpkin
           gprea(1) = 0.0_rp
        end if

        !gprhs = gprhs + gppro
        !gprea = param_tur(3)*gpome 

     end if

  else if(iunkn_tur==2) then
     !
     ! w-equation
     !
     gprhs(1) = gpgra
     gprea(1) = 0.0_rp

     if(gpkin/=0.0_rp) then

        do idime=1,ndime
           gpgrk(idime)=0.0_rp
           gpgrw(idime)=0.0_rp
        end do
        do inode=1,pnode
           do idime=1,ndime
              gpgrk(idime)=gpgrk(idime)&                 ! grad(k)
                   +gpcar(idime,inode)*eltur(1,inode)
              gpgrw(idime)=gpgrw(idime)&                 ! grad(w)
                   +gpcar(idime,inode)*eltur(2,inode)              
           end do
        end do
        grads = 0.0_rp
        do idime=1,ndime
           grads = grads + gpgrk(idime)*gpgrw(idime)     ! grad(k).grad(w)
        end do
        grads = param_tur(4)/gpkin*(gpvis+gpmut)*grads   ! Cw/k*(nu+nut)*grad(k).grad(w)

        gprhs(1) = gprhs(1) + grads + param_tur(5)/gpkin*gppro*gpome
        gprea(1) = param_tur(7)*gpden*gpome 

     else

        gprhs(1) = gprhs(1) + param_tur(5)*gppr2*gpden*0.09_rp*param_tur(6)
        gprea(1) = gprea(1) + param_tur(7)*gpden*gpome

     end if

  end if
  !
  ! GPDIF, GPGRD: diffusion and its gradient
  !
  call tur_elmdif(pnode,gpvis,gpmut,gpcar,eledd,gpdif,gpgrd)

end subroutine tur_brepen
