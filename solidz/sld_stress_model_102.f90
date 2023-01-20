!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_stress_model_102(pgaus,pmate,gpgdi,gpstr,gpdet,flagt,gpdds,&
             gpigd_eps,gpgdi_eps,gpdet_eps,&
             gpmof)
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_elmcla
  ! NAME
  !    sld_stress_model_102
  ! DESCRIPTION
  !    Neo-Hookean material, Ansys' formulation
  !    Compute second Piola-Kirchoff stress tensor S_{IJ}
  !
  !    GPGDI ... Deformation tensor ...................... F = grad(phi)
  !    GPCAU ... Right Cauchy-Green deformation tensor ... C = F^t x F
  !    GPSTR ... 2nd P-K Stress tensor ........................... S
  !    GPPIO ... 1st Piola-Kirchhoff stress tensor ....... P = F.S
  !    GPDDS ... Stress tangent moduli ................... dS/dE
  !    FLAGT ... Flag to activate GPDDS (when implicit)
  ! USES
  ! USED BY
  !    sld_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
  use def_solidz, only       :  parco_sld,kfl_serei_sld
!  use def_solidz, only       :  densi_sld,velas_sld
!  use def_master, only       :  ittim
  implicit none
  integer(ip), intent(in)    :: pgaus,pmate,flagt
  real(rp)                   :: gpcau(ndime,ndime,pgaus),gpmof(pgaus)
  real(rp),    intent(in)    :: gpgdi(ndime,ndime,pgaus),gpdet(pgaus)  !NOTE: F NOT USED
  real(rp)                   :: gpgdi2(ndime,ndime,pgaus)
  real(rp),    intent(out)   :: gpstr(ndime,ndime,pgaus,2)
  real(rp),    intent(out)   :: gpdds(ndime,ndime,ndime,ndime,pgaus)
  real(rp)                   :: gpcin(ndime,ndime),tkron(ndime,ndime)
  integer(ip)                :: igaus,idime,switch,i,j,k,l,ktmo,ltmo
  real(rp)                   :: K0,mu0,tracec,traceb,bidon,detf0,hglas,hfac1,hfac2,logj

  real(rp)        :: gpfin(ndime,ndime,pgaus)
  real(rp)        :: cauch(ndime,ndime,pgaus)

  real(rp)        :: gpstr_tmo(ndime,ndime),gpcau_tmo(ndime,ndime)
  real(rp),    intent(in)        :: gpgdi_eps(ndime,ndime,pgaus)           ! Displacement Gradient F for the perturbed state
  real(rp),    intent(in)        :: gpigd_eps(ndime,ndime,pgaus)           ! Displacement Gradient Inverse F^{-1} for the perturbed state
  real(rp),    intent(in)        :: gpdet_eps(pgaus)                       ! 

  integer(ip), parameter     :: tmeth = 1_ip ! tmeth = 0 => approximation by Belytschko model
                                             ! tmeth = 1 => numerical method using finite diff
  real(rp), parameter        :: ptmo  = 1.0e-8_rp

  !
  ! Neo-Hookean's law
  ! w(C) =
  ! S    =
  !

  mu0  = parco_sld(1,pmate)
  K0   = parco_sld(2,pmate)

  ! Compute sound velocity only at time-step 0 (initialization)
  !if (ittim == 0_ip) then
  !   velas_sld(1,pmate) = sqrt((K0+4.0_rp*mu0/3.0_rp)/densi_sld(1,pmate))
  !end if

  ! Kronecker delta
  tkron = 0.0_rp
  do idime = 1, ndime
      tkron(idime,idime) = 1.0_rp
  end do

  gpstr = 0.0_rp
  gpdds = 0.0_rp

  switch = 1_ip

  if (switch == 1_ip) then

     if (kfl_serei_sld == 1) then

        detf0=0.0_rp
        do igaus=1,pgaus
           detf0=detf0 + 0.125_rp*gpdet(igaus)
        end do

        ! Selective reduced integration: (see de Souza Neto et al. (1996))
        !  replace F = F * ( mean(J) / J )^{1/ndime}
        !  where mean(J) = J_0 (evaluate at the centroid)
        do igaus=1,pgaus
           hglas = (detf0/gpdet(igaus))**(1.0_rp/real(ndime,rp))
           gpgdi2(:,:,igaus) = hglas*gpgdi(:,:,igaus)
        end do

     else

        gpgdi2 = gpgdi

     end if

     ! GPCAU: Cauchy tensor
     do igaus=1,pgaus
        gpcau(:,:,igaus)= matmul(transpose(gpgdi2(:,:,igaus)),gpgdi2(:,:,igaus))
     end do

     do igaus=1,pgaus

        gpcin = tkron

        call invmtx(gpcau(:,:,igaus),gpcin,bidon,ndime) !C^(-1)

        tracec = 0.0_rp
        do idime=1,ndime
           tracec = tracec + gpcau(idime,idime,igaus)
        end do

        ! GPSTR: 2nd Piola--Kirchhoff stress tensor
        gpstr(:,:,igaus,1)= &
              ((-mu0*gpmof(igaus)/3.0_rp)*gpdet(igaus)**(-2.0_rp/3.0_rp)*tracec + &
              gpdet(igaus)*K0*(gpdet(igaus)-1.0_rp))*gpcin + &
              mu0*gpmof(igaus)*gpdet(igaus)**(-2.0_rp/3.0_rp)*tkron

        ! Compute the tangent moduli if required (i.e., implicit scheme)
        if (flagt == 1_ip) then

           selectcase (tmeth)

           case (0)
           ! Approximated by elastic tangent moduli (Belytschko, 2000)
           ! dSdE_{ijkl} = lambda0*C^{-1}_{ij}*C^{-1}_{kl}
           !               + [mu0 - lambda0*(log J)][C^{-1}_{ik}*C^{-1}_{jl}+
           !                                         C^{-1}_{il}*C^{-1}_{jk}]
           logj = log(gpdet(igaus))
           do i=1,ndime; do j=1,ndime; do k=1,ndime; do l=1,ndime
              gpdds(i,j,k,l,igaus)= &
                           (K0-2.0_rp*mu0*gpmof(igaus)/3.0_rp)*gpcin(i,j)*gpcin(k,l)+ &
                           (mu0*gpmof(igaus)-(K0-2.0_rp*mu0*gpmof(igaus)/3.0_rp)*logj)*(gpcin(i,k)*gpcin(j,l)+ &
                                                              gpcin(i,l)*gpcin(j,k))
           enddo; enddo; enddo; enddo

           case (1)
           ! Approximated by numerical tangent moduli using finite difference method
           do ktmo=1,ndime
           do ltmo=1,ktmo
              gpcau_tmo = gpcau(:,:,igaus)
              if (ltmo == ktmo) then
                 gpcau_tmo(ktmo,ltmo) = gpcau_tmo(ktmo,ltmo) + ptmo
              else
                 gpcau_tmo(ktmo,ltmo) = gpcau_tmo(ktmo,ltmo) + 0.5_rp*ptmo
                 gpcau_tmo(ltmo,ktmo) = gpcau_tmo(ktmo,ltmo)
              endif

              gpcin = tkron
              call invmtx(gpcau_tmo,gpcin,bidon,ndime) !C^(-1)

              tracec = 0.0_rp
              do idime=1,ndime
                 tracec = tracec + gpcau_tmo(idime,idime)
              end do

              ! GPSTR: 2nd Piola--Kirchhoff stress tensor
              gpstr_tmo = 0.0_rp
              gpstr_tmo = ((-mu0*gpmof(igaus)/3.0_rp)*sqrt(bidon)**(-2.0_rp/3.0_rp)*tracec + &
                           sqrt(bidon)*K0*(sqrt(bidon)-1.0_rp))*gpcin + &
                           mu0*gpmof(igaus)*sqrt(bidon)**(-2.0_rp/3.0_rp)*tkron

              gpdds(:,:,ktmo,ltmo,igaus)= 2.0_rp*(gpstr_tmo - gpstr(:,:,igaus,1))/ptmo
              if (ltmo /= ktmo) gpdds(:,:,ltmo,ktmo,igaus)= gpdds(:,:,ktmo,ltmo,igaus)
           enddo
           enddo

           endselect

        endif

     end do !ifgaus

  else if (switch == 2_ip) then
     ! * * * * * * *
     ! Calculation of Cauchy stress
     ! (equation from http://en.wikipedia.org/wiki/Neo-Hookean_solid) for validation
     ! * * * * * * *

     bidon=0.0_rp
     do igaus=1,pgaus
        bidon=bidon + 0.125_rp*gpdet(igaus)
     end do

     ! Selective reduced integration: (see de Souza Neto et al. (1996))
     !  replace F = F * ( mean(J) / J )^{1/ndime}
     !  where mean(J) = J_0 (evaluate at the centroid)
     do igaus=1,pgaus
        hglas = (bidon/gpdet(igaus))**(1.0_rp/real(ndime,rp))
        gpgdi2(:,:,igaus) = hglas*gpgdi(:,:,igaus)
     end do

     do igaus=1,pgaus

        !Calculate b=F F^T
        gpfin(:,:,igaus) = matmul(gpgdi2(:,:,igaus),transpose(gpgdi2(:,:,igaus)))

        traceb = 0.0_rp
        do idime=1,ndime
           traceb = traceb + gpfin(idime,idime,igaus)
        end do

        ! when no modulating fields are present, gpmof is 1.0

        hfac1 = mu0*gpmof(igaus)*gpdet(igaus)**(-2.0_rp/3.0_rp)
        hfac2 = K0*gpdet(igaus)*(gpdet(igaus)-1.0_rp) - &
                (mu0*gpmof(igaus)/3.0_rp)*traceb*gpdet(igaus)**(-2.0_rp/3.0_rp)
        cauch(:,:,igaus) = hfac1*gpfin(:,:,igaus) - hfac2*tkron

        call invmtx(gpgdi2(:,:,igaus),gpcin,bidon,ndime) !F^(-1)

        gpstr(:,:,igaus,1) = matmul(gpcin,matmul(cauch(:,:,igaus),transpose(gpcin)))

     end do !igaus

  end if

end subroutine sld_stress_model_102
