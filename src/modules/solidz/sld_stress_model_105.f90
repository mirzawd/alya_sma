!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_stress_model_105(pgaus,pmate,gpstr,flagt,gpdds,&
             gpigd_eps,gpgdi_eps,gpdet_eps,&
             gpmof)
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_elmcla
  ! NAME
  !    sld_elmcla
  ! DESCRIPTION
  !    Compute second Piola-Kirchoff stress tensor S_{IJ}
  !
  !    GPGDI ... Deformation tensor ...................... F = grad(phi)
  !    GPCAU ... Right Cauchy-Green deformation tensor ... C = F^t x F
  !    GPSTR ... 2nd P-K Stress tensor ................... S
  !    GPPIO ... 1st Piola-Kirchhoff stress tensor ....... P = F.S
  !    GPENE ... Stored energy function .................. W
  !    GPDDS ... Stress tangent moduli ................... dS/dE
  !    FLAGT ... Flag to activate GPDDS (when implicit)
  ! USES
  ! USED BY
  !    sld_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
  use def_solidz, only       :  parco_sld,densi_sld,velas_sld
  use def_master, only       :  ittim
  implicit none
  integer(ip), intent(in)    :: pgaus,pmate,flagt
  real(rp),    intent(out)   :: gpstr(ndime,ndime,pgaus,2)
  real(rp)                   :: gpdds(ndime,ndime,ndime,ndime,pgaus),gpmof(pgaus)
  integer(ip)                :: igaus,idime
  real(rp)                   :: lambda0
  real(rp),    intent(in)        :: gpgdi_eps(ndime,ndime,pgaus)           ! Displacement Gradient F for the perturbed state
  real(rp),    intent(in)        :: gpigd_eps(ndime,ndime,pgaus)           ! Displacement Gradient Inverse F^{-1} for the perturbed state
  real(rp),    intent(in)        :: gpdet_eps(pgaus)                       ! 


  !
  ! Simple constant model, kind of a spring
  ! S    = lambda0
  !

  lambda0 = parco_sld(1,pmate)

  ! Compute sound velocity only at time-step 0 (initialization)
  if (ittim == 0_ip) then
     velas_sld(1,pmate) = sqrt(lambda0/densi_sld(1,pmate))
  end if


  gpstr   = 0.0_rp
  do igaus=1,pgaus
     do idime=1,ndime
        gpstr(idime,idime,igaus,1)= lambda0
     end do
  end do

end subroutine sld_stress_model_105
