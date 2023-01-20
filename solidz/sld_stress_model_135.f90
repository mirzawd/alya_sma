!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_stress_model_135(pgaus,pmate,gpcau,gpgdi,gpigd,gpstr,gpdet,pnode,elcod,lnods,ielem,gpsha,nfibe,&
             gpigd_eps,gpgdi_eps,gpdet_eps,&
             gpmof)
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_stress_model_135
  ! NAME
  !    sld__stress_model_134
  ! DESCRIPTION
  !    Law of Usyk and al.
  !
  !    GPGDI ... Deformation tensor ...................... F = grad(u) + I
  !    GPIGD ... Inverse of Deformation tensor ............F^(-1)
  !    GPCAU ... Right Cauchy-Green deformation tensor ... C = F^t x F
  !    GPCAL ... Left Cauchy-Green deformation tensor ...  b = F x F^T
  !    GPSTR ... 2nd P-K Stress tensor ....................S
  !    CASTR ... Cauchy Stress tensor .................... sigma
  !    GPPIO ... 1st Piola-Kirchhoff stress tensor ....... P = F.S
  !    GPENE ... Stored energy function .................. W
  ! USES
  ! USED BY
  !    sld_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mnode
  use def_solidz
  use def_master, only       :  ID_EXMEDI,ITASK_ENDRUN,fiber,gpfib
!  use def_master, only       :  ittim



  implicit none


  integer(ip), intent(in)    :: pgaus,pmate,ielem,pnode,lnods(pnode)
  real(rp),    intent(in)    :: gpcau(ndime,ndime,pgaus),gpdet(pgaus),gpgdi(ndime,ndime,pgaus),&
       gpigd(ndime,ndime,pgaus),gpsha(pnode,pgaus),elcod(ndime,mnode),gpmof(pgaus)
  real(rp),    intent(out)   :: gpstr(ndime,ndime,pgaus,2),nfibe(ndime,pgaus)
  real(rp),    intent(in)        :: gpgdi_eps(ndime,ndime,pgaus)           ! Displacement Gradient F for the perturbed state
  real(rp),    intent(in)        :: gpigd_eps(ndime,ndime,pgaus)           ! Displacement Gradient Inverse F^{-1} for the perturbed state
  real(rp),    intent(in)        :: gpdet_eps(pgaus)                       ! 
  integer(ip)                :: igaus,idime,jdime,kdime,ldime,inode,ipoin
  real(rp)                   :: elfib(3,mnode),lamda(3),ovlam(3),&
       term1,gpgre(ndime,ndime),gpgrf(ndime,ndime)

  real(rp)                   :: nfibe0(ndime,pgaus),norma0(ndime,pgaus),nshet0(ndime,pgaus) !f0, s0, n0 (ref. config)
  real(rp)                   ::                     norma(ndime,pgaus) ,nshet(ndime,pgaus)!f,s,n (fF as define in paper)
  real(rp)                   :: nfibt(ndime,pgaus) ,normt(ndime,pgaus) ,nshtt(ndime,pgaus)  !f,s,n (updated-Unit)
  real(rp)                   :: Kct,cct,bff,bss,bnn,bfs,bfn,bns
  real(rp)                   :: bidon,rotma(3,3),dQdCi(3,3),qfunc,gpcai(3,3)
  real(rp)                   :: i1,i4f,i4s,i4n

  Kct= parco_sld(2,1) !0.001 N/cm² , values of Usyk, 2000
  cct= parco_sld(3,1)
  bff= parco_sld(4,1)
  bss= parco_sld(5,1)
  bnn= parco_sld(6,1)
  bfs= parco_sld(7,1)
  bfn= parco_sld(8,1)
  bns= parco_sld(9,1)

  ! Compute sound velocity only at time-step 0 (initialization)
  !if (ittim == 0_ip) then
  !   velas_sld(1,1) = sqrt(parco_sld(1,1)/densi_sld(1,1))
  !end if

  gpstr = 0.0_rp

  nfibt = 0.0_rp
  nshtt = 0.0_rp
  normt = 0.0_rp

  nfibe = 0.0_rp
  nshet = 0.0_rp
  norma = 0.0_rp

  nfibe0 = 0.0_rp
  norma0 = 0.0_rp
  nshet0 = 0.0_rp

  !
  ! Gather
  !

  if( modfi_sld(pmate) .ne. 0 ) then
     if (modfi_sld(pmate) < 0) then
        do inode = 1,pnode
           ipoin = lnods(inode)
           elfib(1:ndime,inode) = fiber(1:ndime,ipoin)
        end do
     else if (modfi_sld(pmate) > 0) then
        gfibe_sld= 0.0_rp
        gfibe_sld(modfi_sld(pmate))= 1.0_rp
        do inode = 1,pnode
           elfib(1:ndime,inode)= gfibe_sld(1:ndime)
        end do
     end if
  end if

  do igaus=1,pgaus

     ovlam = 0.0_rp
     lamda = 0.0_rp
     gpgre = 0.0_rp
     gpgrf = 0.0_rp
     gpcai = 0.0_rp
     rotma = 0.0_rp
     dQdCi = 0.0_rp

     ! * * * * * *
     ! fibers definition
     ! * * * * * *

     if( modfi_sld(pmate) < 0 ) then
        !
        ! Fibers interpolated on Gauss point from the FIELD (external to this routine)
        !
        call runend("sld_stress_model_135: THIS LAW NEEDS THE THREE VECTORS OF AN ORTHONORMAL BASIS (f,s and n)")
     else !for the prolate spheroid. WARNING - depend on the geom of the spheroid
        call sld_angleh(modfi_sld(pmate),pgaus,pmate,igaus,ielem,elcod,nfibe0,norma0,nshet0,pnode,lnods,gpsha)
     end if

     !E_IJ = 0.5(F_ki*F_kj-Kro_ij) (Holzapfel eq. 2.69) - in the global CS (XYZ)
     do idime=1,ndime
        do jdime=1,ndime
           do kdime=1,ndime
              gpgre(idime,jdime)= gpgre(idime,jdime)+&
                   0.5_rp*(gpgdi(kdime,idime,igaus)*gpgdi(kdime,jdime,igaus))
           end do
        end do
     end do
     do idime=1,ndime
        gpgre(idime,idime)=gpgre(idime,idime)-0.5_rp
     end do

     ! * * * *
     !INVARIANTS
     ! * * * *

     !I_1=Trace(C)
     i1=0.0_rp
     do idime = 1,ndime
        i1 = i1 + gpcau(idime,idime,pgaus)
     end do

     i4f = 0.0_rp
     i4s = 0.0_rp
     i4n = 0.0_rp
     do idime = 1,ndime
        do jdime = 1,ndime
           i4f  =  i4f + nfibe0(idime,igaus) * gpcau(idime,jdime,igaus) * nfibe0(jdime,igaus) ! I_4f  = f0 * C_ij * f0 (eq. 5.1)
           i4s  =  i4s + nshet0(idime,igaus) * gpcau(idime,jdime,igaus) * nshet0(jdime,igaus) ! I_4s  = s0 * C_ij * s0
           i4n  =  i4n + norma0(idime,igaus) * gpcau(idime,jdime,igaus) * norma0(jdime,igaus) ! I_4n  = n0 * C_ij * n0
        end do
     end do

     ! * * * * *

     !Transform f=F'*f0 ; s=F'*s0 ; n=F'*n0

     do idime = 1,ndime
        do jdime = 1,ndime
           nfibe(idime,igaus) = nfibe(idime,igaus) + gpgdi(idime,jdime,igaus) * nfibe0(jdime,igaus)
           nshet(idime,igaus) = nshet(idime,igaus) + gpgdi(idime,jdime,igaus) * nshet0(jdime,igaus)
           norma(idime,igaus) = norma(idime,igaus) + gpgdi(idime,jdime,igaus) * norma0(jdime,igaus)
        end do
     end do


     !Calculate the stretch (lamda). In the fiber direction: lambda^2=f_0*C*f_0 = I_4f

     lamda(1) = sqrt(i4f)
     lamda(2) = sqrt(i4s)
     lamda(3) = sqrt(i4n)

     !      gplep(1,1,igaus) = log(lamda(1))
     !      gplep(2,2,igaus) = log(lamda(2))
     !      gplep(3,3,igaus) = log(lamda(3))

     ovlam(1) = 0.0_rp
     ovlam(2) = 0.0_rp
     ovlam(3) = 0.0_rp
     if( lamda(1) /= 0.0_rp ) ovlam(1) = 1.0_rp / lamda(1)
     if( lamda(2) /= 0.0_rp ) ovlam(2) = 1.0_rp / lamda(2)
     if( lamda(3) /= 0.0_rp ) ovlam(3) = 1.0_rp / lamda(3)
     !
     ! lambda_f * f_true = F * f_0
     ! => here f_true is the UNIT fiber direction in the current configuration
     !
     do idime = 1,ndime
        nfibt(idime,igaus) = nfibe(idime,igaus) * ovlam(1)
        nshtt(idime,igaus) = nshet(idime,igaus) * ovlam(2)
        normt(idime,igaus) = norma(idime,igaus) * ovlam(3)
     end do


     ! * * * * * * * * * * OJO OJO OJO OJO * * * * * * * * * * * *
     !
     ! DETERMIN IF f_0 OR f_t HAS TO BE USED FOR THE ROTATION
     !
     ! * * * * * * * * * * OJO OJO OJO OJO * * * * * * * * * * * *

     !Rotate E_IJ in the (fns) cs system
     do idime=1,ndime
        rotma(1,idime)=nfibt(idime,igaus)
        rotma(2,idime)=nshtt(idime,igaus)
        rotma(3,idime)=normt(idime,igaus)
     end do

     do idime=1,ndime
        do jdime=1,ndime
           do kdime=1,ndime
              do ldime=1,ndime
                 gpgrf(idime,jdime)= gpgrf(idime,jdime)+&
                                !rotma(kdime,idime)*rotma(ldime,jdime)*gpgre(kdime,ldime)
                      rotma(idime,kdime)*rotma(jdime,ldime)*gpgre(kdime,ldime)     !confirm that
              end do
           end do
        end do
     end do

     ! * * * * * *
     ! Passive terms S = 2 * dW/dC
     ! * * * * * *

     qfunc = bff*gpgrf(1,1)**2.0_rp + bss*gpgrf(2,2)**2.0_rp + bnn*gpgrf(3,3)**2.0_rp + &
          bfs*(gpgrf(1,2)**2.0_rp + gpgrf(2,1)**2.0_rp) + &
          bfn*(gpgrf(1,3)**2.0_rp + gpgrf(3,1)**2.0_rp) + &
          bns*(gpgrf(3,2)**2.0_rp + gpgrf(3,2)**2.0_rp)

     !dQ / dC_{ij}
     dQdCi(1,1) = -bff*(0.5_rp - 0.5_rp*gpcau(1,1,igaus))
     dQdCi(2,2) = -bss*(0.5_rp - 0.5_rp*gpcau(2,2,igaus))
     dQdCi(3,3) = -bnn*(0.5_rp - 0.5_rp*gpcau(3,3,igaus))

     dQdCi(1,2) = 0.5_rp*bfs*gpcau(1,2,igaus)
     dQdCi(1,3) = 0.5_rp*bfn*gpcau(1,3,igaus)
     dQdCi(2,3) = 0.5_rp*bns*gpcau(2,3,igaus)

     dQdCi(2,1) = 0.5_rp*bfs*gpcau(2,1,igaus)
     dQdCi(3,1) = 0.5_rp*bfn*gpcau(3,1,igaus)
     dQdCi(3,2) = 0.5_rp*bns*gpcau(3,2,igaus)

     ! S = 2 * dW/dC = cct*exp(Q)*(dQ/dC)
     term1=cct*exp(qfunc)
     do idime = 1,ndime
        do jdime = 1,ndime
           gpstr(idime,jdime,igaus,1)=term1*dQdCi(idime,jdime)
        end do
     end do

     ! * * * * * *
     ! Volumetric term
     ! * * * * * *

     !C^{-1}
     !
     call invmtx(gpcau(1,1,igaus),gpcai(1,1),bidon,ndime)

     term1=Kct*gpdet(igaus)*log(gpdet(igaus))

     do idime = 1,ndime
        do jdime = 1,ndime
           gpstr(idime,jdime,igaus,1) = gpstr(idime,jdime,igaus,1) + &
                term1*gpcai(idime,jdime)
        end do
     end do

     ! Store nfibe for postproc

     do idime = 1,ndime
        gpfib(idime,igaus,ielem) = nfibe(idime,igaus)
     end do

  end do !gauss points

end subroutine sld_stress_model_135
