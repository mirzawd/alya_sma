!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_stress_model_101(pgaus,pmate,gpcau,gpgdi,gpene,gpstr,gpdet,flagt,gpdds,&
             gpigd_eps,gpgdi_eps,gpdet_eps,&
             gpmof)
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_elmcla
  ! NAME
  !    sld_stress_model_101
  ! DESCRIPTION
  !    Neo-Hookean material, Belytschko's formulation
  !    Compute second Piola-Kirchoff stress tensor S_{IJ}
  !
  !    GPGDI ... Deformation tensor ...................... F = grad(phi)
  !    GPCAU ... Right Cauchy-Green deformation tensor ... C = F^t x F
  !    GPSTR ... 2nd P-K Stress tensor ........................... S
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
  use def_solidz, only       :  parco_sld
!  use def_solidz, only       :  densi_sld,velas_sld
!  use def_master, only       :  ittim
  implicit none
  integer(ip), intent(in)    :: pgaus,pmate,flagt
  real(rp),    intent(in)    :: gpcau(ndime,ndime,pgaus),gpmof(pgaus)
  real(rp),    intent(in)    :: gpgdi(ndime,ndime,pgaus),gpdet(pgaus)
  real(rp),    intent(out)   :: gpene(pgaus)
  real(rp),    intent(out)   :: gpstr(ndime,ndime,pgaus,2)
  real(rp),    intent(out)   :: gpdds(ndime,ndime,ndime,ndime,pgaus)
  real(rp)                   :: gpcin(ndime,ndime),tkron(ndime,ndime)
  integer(ip)                :: igaus,idime,jdime,i,j,k,l
  real(rp)                   :: lambda0,mu0,tracec,logj,bidon
  real(rp),    intent(in)        :: gpgdi_eps(ndime,ndime,pgaus)           ! Displacement Gradient F for the perturbed state
  real(rp),    intent(in)        :: gpigd_eps(ndime,ndime,pgaus)           ! Displacement Gradient Inverse F^{-1} for the perturbed state
  real(rp),    intent(in)        :: gpdet_eps(pgaus)                       ! 


  ! Neo-Hookean's law
  ! w(C) = 1/2*lambda0*(log J)^2 - mu0*(log J) + 1/2*mu0*(trace(C)-3)
  ! S    = (lambda0*(log J)-mu0)*C^{-1} + mu0*I

  lambda0 = parco_sld(1,pmate)  ! also called D_1, it is the bulk modulus / 2
  mu0     = parco_sld(2,pmate)  ! also called C_1, it is the shear modulus / 2

  ! Compute sound velocity only at time-step 0 (initialization)
  !if (ittim == 0_ip) then
  !   velas_sld(1,pmate) = sqrt((lambda0+2.0_rp*mu0)/densi_sld(1,pmate))
  !end if

  gpstr = 0.0_rp
  gpdds = 0.0_rp

  ! Kronecker delta
  tkron = 0.0_rp
  do idime = 1, ndime
     tkron(idime,idime) = 1.0_rp
  end do

  do igaus=1,pgaus
     gpcin = tkron

     !     call invmtx(gpgdi(1,1,igaus),work2(1,1),gpdet(igaus),ndime)
     !     logj = 0.0_rp
     !     if (gpdet(igaus) > 1.0E-10) then
     !        logj   = log(gpdet(igaus)) ! J=|det(F)|
     !     end if

     logj = log(gpdet(igaus)) ! J=|det(F)|
     call invmtx(gpcau(:,:,igaus),gpcin,bidon,ndime)

     tracec = 0.0_rp
     do idime=1,ndime
        tracec = tracec + gpcau(idime,idime,igaus)
!!!!  
!!!!  this line is to test an incompressibility parameter in the model 
!!!        gpstr(idime,idime,igaus,1) = gpstr(idime,idime,igaus,1) + gpdet(igaus)*1.0e5*(1.0_rp-(1.0_rp/gpdet(igaus)))
!!!! 
     end do

     do idime=1,ndime
        do jdime=1,ndime
           gpstr(idime,jdime,igaus,1) = gpstr(idime,jdime,igaus,1) + (lambda0*gpmof(igaus)*logj-mu0)*gpcin(idime,jdime) + mu0*tkron(idime,jdime)
        end do
     end do

     ! when no modulating fields are present, gpmof is 1.0

     gpene(igaus)=0.5_rp*lambda0*gpmof(igaus)*logj*logj+mu0*(0.5_rp*(tracec-3.0_rp)-logj)

     ! Compute the tangent moduli if required (i.e., implicit scheme)
     if (flagt == 1_ip) then
        ! Elastic tangent moduli (Belytschko, 2000)
        ! dSdE_{ijkl} = lambda0*C^{-1}_{ij}*C^{-1}_{kl}
        !               + [mu0 - lambda0*(log J)][C^{-1}_{ik}*C^{-1}_{jl}+
        !                                         C^{-1}_{il}*C^{-1}_{jk}]
        do i=1,ndime; do j=1,ndime; do k=1,ndime; do l=1,ndime
           gpdds(i,j,k,l,igaus)= &
                lambda0*gpmof(igaus)*gpcin(i,j)*gpcin(k,l)+ &
                (mu0-lambda0*gpmof(igaus)*logj)*(gpcin(i,k)*gpcin(j,l)+ &
                gpcin(i,l)*gpcin(j,k))
        enddo; enddo; enddo; enddo
     endif

  end do


end subroutine sld_stress_model_101
