!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_stress_model_121(pgaus,pmate,gpcau,gpgdi,gpene,gpstr,gprat,gpdet)
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_stress_model_121
  ! NAME
  !    sld__stress_model_121
  ! DESCRIPTION
  !    Compute second Piola-Kirchoff stress tensor S_{IJ}
  !
  !    GPGDI ... Deformation tensor ...................... F = grad(phi)
  !    GPRAT ... Rate of Deformation tensor .............. Fdot = grad(phidot)
  !    GPCAU ... Right Cauchy-Green deformation tensor ... C = F^t x F
  !    GPSTR ... 2nd P-K Stress tensor ........................... S
  !    GPPIO ... 1st Piola-Kirchhoff stress tensor ....... P = F.S
  !    GPENE ... Stored energy function .................. W
  ! USES
  ! USED BY
  !    sld_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
  use def_solidz, only       :  parco_sld,densi_sld
  implicit none
  integer(ip), intent(in)    :: pgaus,pmate
  real(rp),    intent(in)    :: gpcau(ndime,ndime,pgaus)
  real(rp),    intent(in)    :: gpgdi(ndime,ndime,pgaus)
  real(rp),    intent(in)    :: gprat(ndime,ndime,pgaus)
  real(rp),    intent(out)   :: gpene(pgaus),gpdet(pgaus)
  real(rp),    intent(out)   :: gpstr(ndime,ndime,pgaus)
  real(rp)                   :: work1(ndime,ndime), work2(ndime,ndime)
  integer(ip)                :: igaus,idime,jdime
  real(rp)                   :: c1, c2, eta, logj, tracec

  !
  ! Hyperelastic model from Chappelle and co.
  ! w(C) = eta/2 Edot : Edot + k ...
  ! S = 2 rho0 (c1 I + c2 (I1 I - C) + eta/2 Cdot  + sigma1D fiber X fiber)
  !
  ! being "X" the tensor product
  !

  c1  = parco_sld(1,pmate)
  c2  = parco_sld(2,pmate)
  eta = parco_sld(3,pmate)
  do igaus=1,pgaus
     do idime=1,ndime
        do jdime=1,ndime
           work1(idime,jdime)      = 0.0_rp
           work2(idime,jdime)      = 0.0_rp
        end do
        gpstr(idime,idime,igaus)   = 2.0_rp * densi_sld(1,pmate) * c1
     end do
     call invmtx(gpgdi(1,1,igaus),work2,gpdet(igaus),ndime)
     logj = 0.0_rp
     if (gpdet(igaus) > 1.0E-10_rp) then
        logj   = log(gpdet(igaus)) ! J=|det(F)|
     end if
     call invmtx(gpcau(1,1,igaus),work1,gpdet(igaus),ndime)     ! compute C

     tracec = 0.0_rp
     do idime=1,ndime
        tracec = tracec + gpcau(idime,idime,igaus)
     end do
     do idime=1,ndime
        gpstr(idime,idime,igaus)   = gpstr(idime,idime,igaus) + 2.0_rp * densi_sld(1,pmate) * c2 * tracec
        do jdime=1,ndime
           gpstr(idime,jdime,igaus)=&
                gpstr(idime,jdime,igaus) - 2.0_rp * densi_sld(1,pmate) * c2 * gpcau(idime,idime,igaus) &
                + 0.5_rp * eta * ( &
                gpstr(idime,jdime,igaus) * gprat(jdime,idime,igaus) + &
                gpstr(jdime,idime,igaus) * gprat(idime,jdime,igaus) )
        end do
     end do

     gpene(igaus) = 0.0_rp

  end do



end subroutine sld_stress_model_121
