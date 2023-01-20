!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_kepsil(&
     pnode,eledd,gpcar,D,E,G,P,rho,mu,mut,&
     Cmu,Ce1,Ce2,Ce3,fmu,f1,f2,k,eps,&
     r,dif,Q,grdif)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_kepsil
  ! NAME
  !   tur_kepsil
  ! DESCRIPTION
  !    rho*Dk/Dt-div[(mu+mut/sk) grad(k)] = 
  !               P + G - rho*e - D
  !    rho*De/Dt-div[(mu+mut/se) grad(e)] = 
  !              (Ce1*f1*P + Ce3*G - rho*Ce2*f2*e)*e/k + E
  !
  !    mut = fmu*Cmu*rho*k^2/e
  ! OUTPUT 
  !    R .............. Reaction
  !    DIF ............ Diffusion
  !    Q .............. RHS 
  !    GRDIF(NDIME) ... grad(dif)
  ! USES
  ! USED BY
  !    tur_elmcoe
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime
  use def_turbul, only     :  iunkn_tur,kfl_algor_tur
  implicit none
  integer(ip), intent(in)  :: pnode
  real(rp),    intent(in)  :: eledd(pnode)
  real(rp),    intent(in)  :: gpcar(ndime,pnode)
  real(rp),    intent(in)  :: D,E,G,P,rho,mu,mut
  real(rp),    intent(in)  :: Cmu,Ce1,Ce2,Ce3,fmu,f1,f2,k,eps
  real(rp),    intent(out) :: r(kfl_algor_tur,kfl_algor_tur)
  real(rp),    intent(out) :: dif(kfl_algor_tur)
  real(rp),    intent(out) :: Q(kfl_algor_tur)
  real(rp),    intent(out) :: grdif(kfl_algor_tur,ndime)
  real(rp)                 :: epsok
  ! 
  ! Force and reaction GPRHS and GPREA
  ! 
  if(kfl_algor_tur==1) then
     !
     ! Decoupled
     !
     if(iunkn_tur==1) then
        Q(1)   = P + G - D        
        if(k/=0.0_rp) then
           r(1,1) = rho*rho*fmu*Cmu/mut*k 
        else
           Q(1)   = Q(1) - rho*eps
        end if
     else
        if( k /= 0.0_rp ) then
           epsok  = eps/k
           Q(1)   = (Ce1*f1*P + Ce3*G)*epsok + E 
           r(1,1) = Ce2*f2*rho*epsok 
        else
           Q(1)   = E 
           r(1,1) = 0.0_rp        
        end if
     end if

  else
     !
     ! Coupled
     !
     Q(1)   = P + G - D
     r(1,2) = rho  
     if(k/=0.0_rp) then
        epsok  = eps/k
        Q(2)   = (Ce1*f1*P + Ce3*G)*epsok + E  
        r(2,2) = Ce2*f2*rho*epsok 
     else
        Q(2)   = E          
        r(2,2) = 0.0_rp
     end if
  end if
  !
  ! GPDIF, GPGRD: diffusion and its gradient
  !
  call tur_elmdif(pnode,mu,mut,gpcar,eledd,dif,grdif)

end subroutine tur_kepsil
