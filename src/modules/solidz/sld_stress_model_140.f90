!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!$-----------------------------------------------------------------------
!> @addtogroup SolidzMaterials
!> @{
!> @file    sld_stress_model_40.f90
!> @author  Eva Casoni 
!> @date    15/09/2019
!> @brief   Law of Holzapfel et Ogden modified for compressible anisotropic (Nolan, 2014)
!> @details Law of HGO-MA: only changes the anisotropic part, which now depends on the no-isochoric invariants
!> @} 
!-----------------------------------------------------------------------
subroutine sld_stress_model_140(&
     pgaus,pmate,gpvol,gpcau,gpgdi,gpigd,gpstr,gpdet,nfibe,ielem,elcod,&
     pnode,lnods,gpsha,flagt,gpdds,&
     gpigd_eps,gpgdi_eps,gpdet_eps,&
     gpmof)

  !-----------------------------------------------------------------------
  !****f* Solidz/sld_stress_model_134ma
  ! NAME
  !    sld__stress_model_134ma
  ! DESCRIPTION
  !    Law of Holzapfel et Ogden (2009)
  !
  !    GPGDI ... Deformation tensor ...................... F = grad(u) + I
  !    GPIGD ... Inverse of Deformation tensor ............F^(-1)
  !    GPCAU ... Right Cauchy-Green deformation tensor ... C = F^t x F
  !    GPCAL ... Left Cauchy-Green deformation tensor ...  b = F x F^T
  !    GPSTR ... 2nd P-K Stress tensor ....................S
  !    CASTR ... Cauchy Stress tensor .................... sigma
  !    GPPIO ... 1st Piola-Kirchhoff stress tensor ....... P = F.S
  !    GPENE ... Stored energy function .................. W
  !    GPLEP ... Log strain in the {f s n} system
  !    GPDDS ... Tangent moduli in terms of 2nd P-K stress ....... dS/dE
  !
  !    Special postproc values for this material:
  !    OUTPUT_&_POST-PROCESS
  !    POSTPROCESS LOCEPSILON,   => Log strain in the fibers CS - ln(lambda)
  !    END_OUTPUT_&_POST_PROCESS
  !
  ! USES
  ! USED BY
  !    sld_elmope
  !***
  !-------------------------------------------- ---------------------------
  use def_kintyp, only : ip,rp
  use def_domain, only : ndime,mnode
  use def_solidz
  use def_master, only : fiber,gpfib
  use mod_eccoupling


  implicit none

  integer(ip), intent(in) :: pgaus,pmate,pnode,lnods(pnode),flagt

  integer(ip) :: igaus,idime,jdime,kdime,ldime,ielem,inode,ipoin

  real(rp), intent(in) :: gpcau(ndime,ndime,pgaus),gpdet(pgaus),gpgdi(ndime,ndime,pgaus),         &
   gpmof(pgaus),gpigd(ndime,ndime,pgaus),gpsha(pnode,pgaus),gpgdi_eps(ndime,ndime,pgaus),         &
   gpigd_eps(ndime,ndime,pgaus),gpdet_eps(pgaus),gpvol(pgaus),elcod(ndime,mnode)   

  real(rp), intent(out) :: gpstr(ndime,ndime,pgaus,2),gpdds(ndime,ndime,ndime,ndime,pgaus)

  real(rp) :: nfibe_length,nfibe(ndime,pgaus),                                        &
   gpcal(ndime,ndime),castr_active(ndime,ndime),     &
   a,b,af,bf,as,bs,afs,bfs,scalf,i1,i4f,i4s,i4n,i8fs,term1,gpcai(ndime,ndime),Kct,lamda(3),bidon, &
   ovlam(3),elfib(3,mnode),elfis(3,mnode),elfin(3,mnode),tkron(3,3),     &
   fporf(ndime,ndime),fporf_0(ndime,ndime),spors(ndime,ndime),spors_0(ndime,ndime),  &
   sporf(ndime,ndime),fpors(ndime,ndime),dev_fporf(ndime,ndime), dev_spors(ndime,ndime),          &
   dev_fpors(ndime,ndime),fpors_0(ndime,ndime),spk_vol(ndime,ndime),spk_iso(ndime,ndime),         &
   spk_4f(ndime,ndime),spk_4s(ndime,ndime),spk_8fs(ndime,ndime),ixinvc(ndime,ndime,ndime,ndime),  &
   invcxi(ndime,ndime,ndime,ndime),invcxinvc(ndime,ndime,ndime,ndime),                            &
   i_invc(ndime,ndime,ndime,ndime),i1xi1(ndime,ndime,ndime,ndime),                                &
   dev_ipori(ndime,ndime),J,J23,J43,J53,                       &
   tan_opt_aniso_4f(ndime,ndime,ndime,ndime),tan_opt_aniso_4s(ndime,ndime,ndime,ndime),           &
   tan_opt_aniso_8fs(ndime,ndime,ndime,ndime),tan_opt_iso(ndime,ndime,ndime,ndime),               &
   tan_opt_vol(ndime,ndime,ndime,ndime),term4f,term4s,term8fs,term_vol,iso_i1,iso_i4f,iso_i4n,    &
   iso_i4s,iso_i8fs,fxfxfxf(ndime,ndime,ndime,ndime),                                 &
   i4fp,i4sp,                             &
   sxsxsxs(ndime,ndime,ndime,ndime),nxnxnxn(ndime,ndime,ndime,ndime),                             &
   sxsxfxf(ndime,ndime,ndime,ndime),nxnxfxf(ndime,ndime,ndime,ndime),nporn_0(ndime,ndime),        &
   fxfxi_fxfxi(ndime,ndime,ndime,ndime),sxsxi_sxsxi(ndime,ndime,ndime,ndime), invJ,               &
   castr_passive(ndime,ndime), cau_vol(ndime,ndime), cau_iso(ndime,ndime),                        &
   cau_i4f(ndime,ndime), cau_i4s(ndime,ndime), cau_i8fs(ndime,ndime)
  ! nxnxinvc(ndime,ndime,ndime,ndime),sxsxinvc(ndime,ndime,ndime,ndime),                           &
  ! fxfxinvc(ndime,ndime,ndime,ndime),fsxfsxinvc(ndime,ndime,ndime,ndime),                         &
!   i8fsxi8fs(ndime,ndime,ndime,ndime),invcxfsxfs(ndime,ndime,ndime,ndime)

  real(rp) :: nfibe0(ndime,pgaus),norma0(ndime,pgaus),nshet0(ndime,pgaus) !f0, s0, n0 (ref. config)
  real(rp) :: norma(ndime,pgaus),nshet(ndime,pgaus)                                             !f,s,n (fF as define in paper) (nfibe is "output")
  real(rp) :: nfibt(ndime,pgaus),normt(ndime,pgaus),nshtt(ndime,pgaus)                          !f,s,n (updated-Unit)

  integer(ip) :: kfl_yeoh
  real(rp) :: k1, k2, k3, mu1, mu2, mu3


  ! Retrieve the parameters
  scalf = 1.0_rp !default value

  Kct = parco_sld(2,pmate) !0.001 N/cm2 , values of Usyk, 2000
  a   = parco_sld(3,pmate)
  b   = parco_sld(4,pmate)
  af  = parco_sld(5,pmate)
  bf  = parco_sld(6,pmate)
  as  = parco_sld(7,pmate)
  bs  = parco_sld(8,pmate)
  afs = parco_sld(9,pmate)
  bfs = parco_sld(10,pmate)
  scalf = parco_sld(11,pmate)

  ! Yeoh model (McEvoy, 2018)
  kfl_yeoh = 0_ip
  if (kfl_yeoh == 1) then
   !porcine
  !  k1 = 18.23_rp 
  !  k2 = 145.8_rp
  !  k3 = 1.6e6_rp
  !  mu1 = 2.44_rp
  !  mu2 = -6.04_rp
  !  mu3 = 14.56_rp
   !human
    k1 = 7.31_rp 
    k2 = 5.12_rp
    k3 = 2.49e6_rp
    mu1 = 0.98_rp
    mu2 = -0.212_rp
    mu3 = 22.68_rp
  end if

  ! Initialise some needed matrices and vectors
  gpcai = 0.0_rp
  gpstr = 0.0_rp
  gpdds = 0.0_rp

  nfibt = 0.0_rp
  nshtt = 0.0_rp
  normt = 0.0_rp

  nfibe = 0.0_rp
  nshet = 0.0_rp
  norma = 0.0_rp

  nfibe0 = 0.0_rp
  norma0 = 0.0_rp
  nshet0 = 0.0_rp

  elfib = 0.0_rp

  tkron = 0.0_rp
  do idime = 1,ndime
     tkron(idime,idime) = 1.0_rp
  end do

  !
  ! Gather
  !
  if (modfi_sld(pmate) .ne. 0) then
    elfis = 0.0_rp           
    elfin = 0.0_rp
    if (kfl_fiber_sld > 0) then
      if (modfi_sld(pmate) < 0) then
        do inode = 1,pnode
          ipoin = lnods(inode)
          elfib(1:ndime,inode) = fiber(1:ndime,ipoin)
        end do
      else if (modfi_sld(pmate) > 0) then
        gfibe_sld = 0.0_rp
        gfibe_sld(modfi_sld(pmate)) = 1.0_rp
        do inode = 1,pnode
          elfib(1:ndime,inode) = gfibe_sld(1:ndime)
        end do
      end if

      if (kfl_fiber_sld == 2 .or. kfl_fiber_sld == 3) then
        do inode = 1,pnode
          ipoin = lnods(inode)

          elfis(1,inode) = fibts_sld(1,ipoin)
          elfis(2,inode) = fibts_sld(2,ipoin)
          elfis(3,inode) = fibts_sld(3,ipoin)    

          elfin(1,inode) = fibtn_sld(1,ipoin)
          elfin(2,inode) = fibtn_sld(2,ipoin)
          elfin(3,inode) = fibtn_sld(3,ipoin)  
        end do
      end if
    end if
  end if

  ! Gauss points loop
  do igaus = 1,pgaus

    !Initialize active, passive and volumetric stresses to zero
    gpcai = 0.0_rp
    castr_active = 0.0_rp

    !Calcul of b = F F^T
    gpcal = 0.0_rp
    do idime = 1,ndime
      do jdime = 1,ndime
        do kdime = 1,ndime
          gpcal(idime,jdime) = gpcal(idime,jdime)+gpgdi(idime,kdime,igaus)*                       &
           gpgdi(jdime,kdime,igaus)
        end do
      end do
    end do

    if( modfi_sld(pmate) .ne. 0 ) then

      !
      ! Fibers interpolated on Gauss point
      !
      nfibe0(1:3,igaus) = 0.0_rp
      do inode = 1,pnode
        nfibe0(1:3,igaus) = nfibe0(1:3,igaus) + gpsha(inode,igaus)*elfib(1:3,inode)
      end do

      !Normalized f_0
      bidon = sqrt(nfibe0(1,igaus)**2.0_rp+nfibe0(2,igaus)**2.0_rp+nfibe0(3,igaus)**2.0_rp)
      if (bidon == 0.0_rp) bidon= 1.0_rp
      do idime = 1,ndime
        nfibe0(idime,igaus) = nfibe0(idime,igaus)/bidon
      end do
    end if

    if (kfl_fiber_sld == 2 .or. kfl_fiber_sld == 3) then

      !
      ! Orthotropic material
      !
      norma0(1:3,igaus) = 0.0_rp
      nshet0(1:3,igaus) = 0.0_rp
      do inode = 1,pnode
        norma0(1:3,igaus) = norma0(1:3,igaus)+gpsha(inode,igaus)*elfin(1:3,inode)
        nshet0(1:3,igaus) = nshet0(1:3,igaus)+gpsha(inode,igaus)*elfis(1:3,inode)
      end do

      !Normalized n_0
      bidon = sqrt(norma0(1,igaus)**2.0_rp+norma0(2,igaus)**2.0_rp+norma0(3,igaus)**2.0_rp)
      if (bidon == 0.0_rp) bidon= 1.0_rp
      do idime=1,ndime
        norma0(idime,igaus)=norma0(idime,igaus)/bidon
      end do

      !Normalized s_0
      bidon = sqrt(nshet0(1,igaus)**2.0_rp+nshet0(2,igaus)**2.0_rp+nshet0(3,igaus)**2.0_rp)
      if (bidon == 0.0_rp) bidon = 1.0_rp
      do idime = 1,ndime
        nshet0(idime,igaus) = nshet0(idime,igaus)/bidon
      end do
    end if

    ! * * * *
    !INVARIANTS
    ! * * * *

    i1 = 0.0_rp
    do idime = 1,ndime
      i1 = i1 + gpcau(idime,idime,igaus)
    end do

    ! Compute the invariants
    i4f = 0.0_rp
    i4s = 0.0_rp
    i4n = 0.0_rp
    i8fs = 0.0_rp

    do idime = 1,ndime
      do jdime = 1,ndime
        i4f = i4f+nfibe0(idime,igaus)*gpcau(idime,jdime,igaus)*nfibe0(jdime,igaus) ! I_4f  = f0 * C_ij * f0 (eq. 5.1)
        i4s = i4s+nshet0(idime,igaus)*gpcau(idime,jdime,igaus)*nshet0(jdime,igaus) ! I_4s  = s0 * C_ij * s0
        i4n = i4n+norma0(idime,igaus)*gpcau(idime,jdime,igaus)*norma0(jdime,igaus) ! I_4n  = n0 * C_ij * n0
        i8fs = i8fs+0.5_rp*(nfibe0(idime,igaus)*gpcau(idime,jdime,igaus)*nshet0(jdime,igaus)+     &
         nshet0(idime,igaus)*gpcau(idime,jdime,igaus)*nfibe0(jdime,igaus)) ! I_8fs = f0 * C_ij * s0 (eq. 5.3)

        !Transform f=F'*f0 ; s=F'*s0 ; n=F'*n0 (PUSH FORWARD)
        nfibe(idime,igaus) = nfibe(idime,igaus)+gpgdi(idime,jdime,igaus)*nfibe0(jdime,igaus)
        nshet(idime,igaus) = nshet(idime,igaus)+gpgdi(idime,jdime,igaus)*nshet0(jdime,igaus)
        norma(idime,igaus) = norma(idime,igaus)+gpgdi(idime,jdime,igaus)*norma0(jdime,igaus)
      end do
    end do

    do idime = 1,ndime
      gpfib(idime,igaus,ielem) = nfibe(idime,igaus)
    end do

    !
    ! f x f , s x s, f x s, s x f (outer products [3x3])
    !
    do idime = 1,ndime
      do jdime = 1,ndime
        fporf_0(idime,jdime) = nfibe0(idime,igaus)*nfibe0(jdime,igaus)
        spors_0(idime,jdime) = nshet0(idime,igaus)*nshet0(jdime,igaus)
        nporn_0(idime,jdime) = norma0(idime,igaus)*norma0(jdime,igaus)
        fpors_0(idime,jdime) = 0.5_rp*(nfibe0(idime,igaus)*nshet0(jdime,igaus)+                   &
         nshet0(idime,igaus)*nfibe0(jdime,igaus))

        fporf(idime,jdime) = nfibe(idime,igaus)*nfibe(jdime,igaus)
        spors(idime,jdime) = nshet(idime,igaus)*nshet(jdime,igaus)
        !nporn(idime,jdime) = norma(idime,igaus)*norma(jdime,igaus)
        fpors(idime,jdime) = 0.5_rp*(nfibe(idime,igaus)*nshet(jdime,igaus)+                   &
         nshet(idime,igaus)*nfibe(jdime,igaus))
      end do
    end do

    ! * * * * * *
    ! Second Piola-Kirchhoff stress tensor
    ! * * * * * *

    ! Different needed variables related to the Deformation Gradient
    J   = gpdet(igaus)
    invJ = 1.0_rp/gpdet(igaus)
    J23 = J**(-(2.0_rp/3.0_rp))
    J43 = J**(-(4.0_rp/3.0_rp))
    J53 = J**(-(5.0_rp/3.0_rp))

    ! Different expressions of the scalar portion of the volumetric stress tensor
    !term_vol = (Kct/2.0_rp)*(J-1.0_rp)*J
    !term_vol = (Kct/4.0_rp)*(J*(J-1.0_rp)+log(J))
    term_vol = (Kct/4.0_rp)*((J**2.0_rp)-1.0_rp)

    ! Compute the isochoric invariants
    iso_i1 = J23*i1
    iso_i4f = J23*i4f
    iso_i4s = J23*i4s
    iso_i4n = J23*i4n
    iso_i8fs = J23*i8fs

    ! Tension/compression asymmetry
    i4fp = 0.5_rp*((iso_i4f-1.0_rp)+abs(iso_i4f-1.0_rp))
    i4sp = 0.5_rp*((iso_i4s-1.0_rp)+abs(iso_i4s-1.0_rp))

    ! exponential isotropic and isochoric term
    term1 = exp(b*(iso_i1-3.0_rp))

    if (kfl_fiber_sld > 0) then

      ! these checks re necessary, because round-off errors from i4f/s can create spurious small forces in the passive 
      if (abs((i4f-1.0_rp)) < 1.0e-12) i4f = 1.0_rp
      if (abs((i4s-1.0_rp)) < 1.0e-12) i4s = 1.0_rp
 
      if (abs((iso_i4f-1.0_rp)) < 1.0e-12) iso_i4f = 1.0_rp
      if (abs((iso_i4s-1.0_rp)) < 1.0e-12) iso_i4s = 1.0_rp

      ! exponential anisotropic and isochoric term
      term4f = exp(bf*(i4f-1.0_rp)**2.0_rp)
      term4s = exp(bs*(i4s-1.0_rp)**2.0_rp)
      term8fs = exp(bfs*(i8fs**2.0_rp))

    end if

    ! Invert Right Cauchy-Green strain tensor C (gpcai)
    call invmtx(gpcau(:,:,igaus),gpcai,bidon,ndime)

    !
    ! Deviatoric part: Dev(Aij) = Aij - (1/3)Akk tkron_ij
    !
    dev_ipori = tkron-(1.0_rp/3.0_rp)*i1*gpcai
    dev_fporf = fporf_0-(1.0_rp/3.0_rp)*i4f*gpcai
    dev_spors = spors_0-(1.0_rp/3.0_rp)*i4s*gpcai
    dev_fpors = fpors_0-(1.0_rp/3.0_rp)*i8fs*gpcai

    ! Second Piola-Kirchhoff stress tensor
    spk_vol = 2.0_rp*term_vol*gpcai   !this term remains equal as in the original formulation
    spk_iso = J23*a*term1*dev_ipori   !this term remains equal as in the original formulation
    spk_4f = 2.0_rp*af*(i4f-1.0_rp)*term4f*fporf
    spk_4s = 2.0_rp*as*(i4s-1.0_rp)*term4s*spors
    spk_8fs = 2.0_rp*afs*(i8fs-1.0_rp)*term8fs*fpors

 !   gpstr(1:3,1:3,igaus,1) = spk_vol(1:3,1:3)+spk_iso(1:3,1:3)+spk_4f(1:3,1:3)+spk_4s(1:3,1:3)+   &
 !    spk_8fs(1:3,1:3)

!------------------------------------------------------------------------------------
    !Calculate Caucy stress
    if (kfl_yeoh == 1) then
      cau_vol = k1*(J-1.0_rp)*tkron + 2.0_rp*k2*((J-1)**3.0_rp)*tkron + 3.0_rp*k3*((J-1)**5.0_rp)*tkron
    else
      !cau_vol = (Kct/4.0_rp)*(J-1.0_rp)*tkron
      cau_vol = Kct*(J-1.0_rp)*tkron
    end if

    if (kfl_yeoh == 1) then
      cau_iso = mu1*J23*(gpcal-(1.0_rp/3.0_rp)*i1*tkron) + 2.0_rp*mu2*J23*(J23*i1-3.0_rp)*(gpcal-(1.0_rp/3.0_rp)*i1*tkron) + &
                3.0_rp*mu3*J23*((J23*i1-3.0_rp)**2.0_rp)*(gpcal-(1.0_rp/3.0_rp)*i1*tkron) 
    else
      !cau_iso = a*exp(b*(i1-3.0_rp))*gpcal
      cau_iso = a*exp(b*(i1-3.0_rp))*(gpcal-(1.0_rp/3.0_rp)*i1*tkron)
    end if

    cau_i4f = 2.0_rp*af*(i4f-1.0_rp)*term4f*fporf
    cau_i4s = 2.0_rp*as*(i4s-1.0_rp)*term4s*spors
    cau_i8fs = afs*i8fs*term8fs*(fpors + sporf)

    castr_passive = cau_iso + cau_i4f + cau_i4s + cau_i8fs    

    !Transform Cauchy stress to S: S = J*F^-1*sigma*F^-T
    do idime = 1,ndime
       do jdime = 1,ndime
          do kdime = 1,ndime
             do ldime = 1,ndime
                gpstr(idime,jdime,igaus,1) = gpstr(idime,jdime,igaus,1)+ &
                     gpdet(igaus) *&
                     gpigd(idime,kdime,igaus)*(castr_passive(kdime,ldime) + cau_vol(kdime,ldime)* gpigd(jdime,ldime,igaus))
             end do
          end do
       end do
    end do


!------------------------------------------------------------------------------------

    ! Calculate the stretch (lambda). In the fiber direction: lambda^2=f_0*C*f_0 = I4f
    lamda(1) = sqrt(i4f)
    lamda(2) = sqrt(i4s)
    lamda(3) = sqrt(i4n)

    !
    ! Updated fiber direction: lambda_f * f_true = F * f_0 = f
    ! => here f is the UNIT fiber direction in the current configuration
    ! used for active stress and postproc
    !
    ovlam(1) = 0.0_rp
    ovlam(2) = 0.0_rp
    ovlam(3) = 0.0_rp
    if( lamda(1) /= 0.0_rp ) then
      !         gplep(1,1,igaus) = log(lamda(1))
      ovlam(1) = 1.0_rp / lamda(1)
    end if
    if( lamda(2) /= 0.0_rp ) then
      !         gplep(2,2,igaus) = log(lamda(2))
      ovlam(2) = 1.0_rp / lamda(2)
    end if
    if( lamda(3) /= 0.0_rp ) then
      !         gplep(3,3,igaus) = log(lamda(3))
      ovlam(3) = 1.0_rp / lamda(3)
    end if

    nfibe_length= 0.0_rp
    do idime = 1,ndime
      nfibt(idime,igaus) = nfibe(idime,igaus)*ovlam(1)
      nshtt(idime,igaus) = nshet(idime,igaus)*ovlam(2)
      normt(idime,igaus) = norma(idime,igaus)*ovlam(3)
      nfibe_length = nfibe_length + nfibe(idime,igaus)*nfibe(idime,igaus) 
    end do

    ! * * * * * * * *
    ! Tangent moduli
    ! * * * * * * * *
    ! ========================================================================

    ! Passive component
    ! ========================================================================
    if (flagt == 1_ip) then
      do idime = 1,ndime
        do jdime = 1,ndime
          do kdime = 1,ndime
            do ldime = 1,ndime
              fxfxfxf(idime,jdime,kdime,ldime) = nfibe0(idime,igaus)*nfibe0(jdime,igaus)*         &
               nfibe0(kdime,igaus)*nfibe0(ldime,igaus)
              sxsxsxs(idime,jdime,kdime,ldime) = nshet0(idime,igaus)*nshet0(jdime,igaus)*         &
               nshet0(kdime,igaus)*nshet0(ldime,igaus)
              nxnxnxn(idime,jdime,kdime,ldime) = norma0(idime,igaus)*norma0(jdime,igaus)*         &
               norma0(kdime,igaus)*norma0(ldime,igaus)
              sxsxfxf(idime,jdime,kdime,ldime) = nshet0(idime,igaus)*nshet0(jdime,igaus)*         &
               nfibe0(kdime,igaus)*nfibe0(ldime,igaus)
              nxnxfxf(idime,jdime,kdime,ldime) = norma0(idime,igaus)*norma0(jdime,igaus)*         &
               nfibe0(kdime,igaus)*nfibe0(ldime,igaus)                          
              fxfxi_fxfxi(idime,jdime,kdime,ldime) = nfibe0(jdime,igaus)*nfibe0(ldime,igaus)*     &
               tkron(idime,kdime) + nfibe0(idime,igaus)*nfibe0(ldime,igaus)*tkron(jdime,kdime)
              sxsxi_sxsxi(idime,jdime,kdime,ldime) = nshet0(jdime,igaus)*nshet0(ldime,igaus)*     &
               tkron(idime,kdime) + nshet0(idime,igaus)*nshet0(ldime,igaus)*tkron(jdime,kdime)
            end do
          end do
        end do
      end do

      tan_opt_iso = 2.0_rp*J23*a*term1*(-(1.0_rp/3.0_rp)*ixinvc-(1.0_rp/3.0_rp)*invcxi+           &
       (1.0_rp/9.0_rp)*i1*invcxinvc+(1.0_rp/6.0_rp)*i1*i_invc)+2.0_rp*J43*a*b*term1*(i1xi1)

      tan_opt_aniso_4f = 2.0_rp*af*invJ*(i4f-1.0_rp)*term4f*fxfxi_fxfxi + &
       4.0_rp*af*invJ*(2.0_rp*bf*(i4f-1)**2.0_rp + 1.0_rp)*term4f*fxfxfxf 

      tan_opt_aniso_4s = 2.0_rp*as*invJ*(i4s-1.0_rp)*term4s*sxsxi_sxsxi + &
       4.0_rp*as*invJ*(2.0_rp*bs*(i4s-1)**2.0_rp + 1.0_rp)*term4s*sxsxsxs 

      !tan_opt_aniso_8fs = 4.0_rp*J23*afs*iso_i8fs*term8fs*(-(1.0_rp/3.0_rp)*fsxfsxinvc-           &
      ! (1.0_rp/3.0_rp)*invcxfsxfs+(1.0_rp/9.0_rp)*i8fs*invcxinvc+(1.0_rp/6.0_rp)*i8fs*i_invc)+    &
      ! 4.0_rp*J43*afs*term8fs*i8fsxi8fs+8.0_rp*J43*afs*bfs*(iso_i8fs**2.0_rp)*term8fs*i8fsxi8fs

      !tan_opt_vol = 2.0_rp*Kct*J*((2.0_rp*J-1.0_rp)*invcxinvc-0.5_rp*(J-1.0_rp)*i_invc)
      !tan_opt_vol = (Kct/2.0_rp)*(((2.0_rp*J-1.0_rp)*J+1.0_rp)*invcxinvc-(((J**2.0_rp)-J)+        &
      ! log(J))*i_invc)
      tan_opt_vol = Kct*((J**2.0_rp)*invcxinvc-0.5_rp*((J**2.0_rp)-1.0_rp)*i_invc)

      gpdds(1:3,1:3,1:3,1:3,igaus) = tan_opt_iso(1:3,1:3,1:3,1:3)+                                &
       tan_opt_aniso_4f(1:3,1:3,1:3,1:3)+tan_opt_aniso_4s(1:3,1:3,1:3,1:3)+                       &
       tan_opt_aniso_8fs(1:3,1:3,1:3,1:3)+tan_opt_vol(1:3,1:3,1:3,1:3)

    end if
    ! ========================================================================

  end do !gauss points

end subroutine sld_stress_model_140

