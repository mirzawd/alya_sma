!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!$-----------------------------------------------------------------------
!> @addtogroup SolidzMaterials
!> @{
!> @file    sld_stress_model_137.f90
!> @author  Waleed Miraza and Eva Casoni 01/07/2017 
!> @date    01/08/2016
!> @brief   Law of Holzapfel et Ogden modified for compressible anisotropic (Nolan, 2014)
!> @details Law of HGO-MA
!> @} 
!-----------------------------------------------------------------------
subroutine sld_stress_model_137hgoc(&
     pgaus,pmate,gpvol,gpcau,gpgdi,gpigd,gpstr,gpdet,&
     nfibe,ielem,elcod,pnode,lnods,gpsha,flagt,gpdds,&
             gpigd_eps,gpgdi_eps,gpdet_eps,&
             gpmof)
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_stress_model_134b
  ! NAME
  !    sld__stress_model_134b
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
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mnode
  use def_solidz
  use def_master, only       :  ITASK_ENDRUN,fiber,gpfib,coupling


  implicit none
  integer(ip), intent(in)    :: pgaus,pmate,pnode,lnods(pnode),flagt
  real(rp),    intent(in)    :: gpcau(ndime,ndime,pgaus),gpdet(pgaus),gpgdi(ndime,ndime,pgaus),gpmof(pgaus)
  real(rp),    intent(in)    :: gpigd(ndime,ndime,pgaus),gpsha(pnode,pgaus)
  real(rp),    intent(out)   :: gpstr(ndime,ndime,pgaus,2)
  real(rp),    intent(out)   :: gpdds(ndime,ndime,ndime,ndime,pgaus)
  real(rp),    intent(in)        :: gpgdi_eps(ndime,ndime,pgaus)           ! Displacement Gradient F for the perturbed state
  real(rp),    intent(in)        :: gpigd_eps(ndime,ndime,pgaus)           ! Displacement Gradient Inverse F^{-1} for the perturbed state
  real(rp),    intent(in)        :: gpdet_eps(pgaus)                       ! 
  real(rp)                   :: nfibe(ndime,pgaus),nfibe_2(ndime,pgaus)

  integer(ip)                :: igaus,idime,jdime,kdime,ldime,ielem,inode,ipoin
  real(rp)                   :: gpcal(ndime,ndime),castr(ndime,ndime)

  real(rp)                   :: nfibe0(ndime,pgaus),nfibe0_2(ndime,pgaus), norma0(ndime,pgaus), nshet0(ndime,pgaus) !f0, s0, n0 (ref. config)
  real(rp)                   ::                     norma(ndime,pgaus) ,nshet(ndime,pgaus)  !f,s,n (fF as define in paper) (nfibe is "output")
  real(rp)                   :: nfibt(ndime,pgaus) ,normt(ndime,pgaus) ,nshtt(ndime,pgaus)  !f,s,n (updated-Unit)

  real(rp)                   :: a,af,bf,flag_2

  real(rp)                   :: i1,i4f,i6f
  real(rp)                   :: term1,term2,term3,term4,term1_2,term2_2,term3_2,term4_2
  real(rp)                   :: gpcai(ndime,ndime),bidon,bidon_2,Kct

  real(rp)                   :: elfib(3,mnode),elfib_2(3,mnode)

  real(rp),    intent(in)    :: gpvol(pgaus)
  real(rp),    intent(in)    :: elcod(ndime,mnode)

  real(rp)                   :: tkron(3,3)

  real(rp)                   :: fporf(ndime,ndime), fporf_0(ndime,ndime), fporf_2(ndime,ndime), fporf_0_2(ndime,ndime)
 
  real(rp)                   :: trace_fporf,trace_fporf_2

  real(rp)                   :: energy_a1,energy_a1_2,energy_a1_a1,energy_a1_a1_2,&
                                bulk_1,bulk_2,bulk_3,&
                                shear_1,shear_2,shear_3,shear_4,shear,C4,Dev(ndime,ndime),Dev_2(ndime,ndime),&
                                invC_bar(ndime,ndime)
  
  real(rp)                   :: dev_fporf(ndime,ndime),dev_fporf_2(ndime,ndime),&
                                cauchy_iso(ndime,ndime),cauchy_aniso1(ndime,ndime),cauchy_aniso2(ndime,ndime),cauchy_vol(ndime,ndime)
 

  !real(rp)                   :: Fpert(ndime,ndime,pgaus), Cpert(ndime,ndime), bpert(ndime,ndime), epsil






  flag_2 = 0_ip

  Kct = parco_sld( 2,pmate) ! bulk modulus  Phi_vol = 0.5*Kct*(J-1)^2
  a   = parco_sld( 3,pmate) ! mu (shear modulus)
  af  = parco_sld( 4,pmate) ! k1 (in pascals)
  bf  = parco_sld( 5,pmate) ! k2 (adimensional)

  gpcai = 0.0_rp
  gpstr = 0.0_rp
  gpdds = 0.0_rp

  nfibt = 0.0_rp
  nshtt = 0.0_rp
  normt = 0.0_rp

  nfibe = 0.0_rp
  nfibe_2 = 0.0_rp
  nshet = 0.0_rp
  norma = 0.0_rp

  nfibe0 = 0.0_rp
  nfibe0_2 = 0.0_rp
  norma0 = 0.0_rp
  nshet0 = 0.0_rp

  elfib = 0.0_rp
  elfib_2 = 0.0_rp

  tkron = 0.0_rp
  do idime=1,ndime
     tkron(idime,idime) = 1.0_rp
  end do



  !
  ! Gather
  !

  if( modfi_sld(pmate) .ne. 0 ) then
     if (kfl_fiber_sld > 0) then
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
  end if

  do igaus = 1,pgaus

     !Initialize active, passive and volumetric stresses to zero
     gpcai = 0.0_rp

     !Calcul of b = F F^T
     gpcal = 0.0_rp
     do idime = 1,ndime
        do jdime = 1,ndime
           do kdime = 1,ndime
              gpcal(idime,jdime) = &
                   gpcal(idime,jdime) + gpgdi(idime,kdime,igaus) * gpgdi(jdime,kdime,igaus)
           end do
        end do
     end do
 
     if( modfi_sld(pmate) .ne. 0 ) then
        !
        ! Fibers interpolated on Gauss point
        !
        nfibe0(1,igaus) = 0.0_rp
        nfibe0(2,igaus) = 0.0_rp
        nfibe0(3,igaus) = 0.0_rp
        nfibe0_2(1,igaus) = 0.0_rp
        nfibe0_2(2,igaus) = 0.0_rp
        nfibe0_2(3,igaus) = 0.0_rp


        do inode = 1,pnode
          nfibe0(1,igaus) = nfibe0(1,igaus) + gpsha(inode,igaus)*elfib(1,inode)
          nfibe0(2,igaus) = nfibe0(2,igaus) + gpsha(inode,igaus)*elfib(2,inode)
          nfibe0(3,igaus) = nfibe0(3,igaus) + gpsha(inode,igaus)*elfib(3,inode)

          nfibe0_2(1,igaus) = nfibe0_2(1,igaus) + gpsha(inode,igaus)*elfib(1,inode)
          nfibe0_2(2,igaus) = nfibe0_2(2,igaus) + gpsha(inode,igaus)*(elfib(2,inode))
          nfibe0_2(3,igaus) = nfibe0_2(3,igaus) + gpsha(inode,igaus)*(-elfib(3,inode))
        end do

        !Normalized f_0
        bidon=sqrt(nfibe0(1,igaus)**2.0_rp+nfibe0(2,igaus)**2.0_rp+nfibe0(3,igaus)**2.0_rp)
        bidon_2=sqrt(nfibe0_2(1,igaus)**2.0_rp+nfibe0_2(2,igaus)**2.0_rp+nfibe0_2(3,igaus)**2.0_rp)
        if (bidon == 0.0_rp) bidon= 1.0_rp
        if (bidon_2 == 0.0_rp) bidon_2 = 1.0_rp
        do idime=1,ndime
           nfibe0(idime,igaus)=nfibe0(idime,igaus)/bidon
           nfibe0_2(idime,igaus)=nfibe0_2(idime,igaus)/bidon_2
        end do

     end if

     ! * * * *
     !INVARIANTS
     ! * * * *

     !I_1=Trace(C)
     i1=0.0_rp
     do idime = 1,ndime
        i1 = i1 + gpcau(idime,idime,pgaus)
     end do

     i4f = 0.0_rp
     i6f = 0.0_rp

     do idime = 1,ndime
        do jdime = 1,ndime
           i4f  =  i4f + nfibe0(idime,igaus) * gpdet(igaus)**(-2.0_rp/3.0_rp)*gpcau(idime,jdime,igaus) * nfibe0(jdime,igaus) ! I_4f  = f0 * C_ij * f0 (eq. 5.1)
           i6f  =  i6f + nfibe0_2(idime,igaus) * gpdet(igaus)**(-2.0_rp/3.0_rp)*gpcau(idime,jdime,igaus) * nfibe0_2(jdime,igaus) ! I_4f  = f0 * C_ij * f0 (eq. 5.1)

           ! * * * * *
           !Transform f=F'*f0  (PUSH FORWARD)
           nfibe(idime,igaus) = nfibe(idime,igaus) + (gpdet(igaus)**(-1.0_rp/3.0_rp))*gpgdi(idime,jdime,igaus) * nfibe0(jdime,igaus)
           nfibe_2(idime,igaus) = nfibe_2(idime,igaus) + (gpdet(igaus)**(-1.0_rp/3.0_rp))*gpgdi(idime,jdime,igaus) * nfibe0_2(jdime,igaus)
        end do
     end do

     do idime = 1,ndime
        gpfib(idime,igaus,ielem) = nfibe(idime,igaus)
     end do

     !
     ! f x f (outer products)
     !
     fporf = 0.0_rp
     fporf_2 = 0.0_rp

     do idime = 1,ndime
        do jdime = 1,ndime
           fporf(idime,jdime) = nfibe(idime,igaus) * nfibe(jdime,igaus)
           fporf_0(idime,jdime) = nfibe0(idime,igaus) * nfibe0(jdime,igaus)
           fporf_0_2(idime,jdime) = nfibe0_2(idime,igaus) * nfibe0_2(jdime,igaus)
           fporf_2(idime,jdime) = nfibe_2(idime,igaus) * nfibe_2(jdime,igaus)
        end do
     end do


     ! * * * * * *
     ! Anisotropy terms (related to both family of fibers)
     ! * * * * * *
     
     term2= 0.0_rp
     term2_2 = 0.0_rp

     if (kfl_fiber_sld > 0) then
        term2 = 2.0_rp*af*(i4f-1.0_rp)*exp(bf *(i4f-1.0_rp)**2.0_rp)/gpdet(igaus)
        term2_2 = 2.0_rp*af*(i6f-1.0_rp)*exp(bf *(i6f-1.0_rp)**2.0_rp)/gpdet(igaus)
     end if

     !
     ! Deviatoric part: Dev(Aij) = Aij - (1/3)Akk tkron_ij
     !
     trace_fporf = 0.0_rp
     trace_fporf_2 = 0.0_rp

     do idime = 1,ndime
        do jdime = 1,ndime
           trace_fporf = trace_fporf + tkron(idime,jdime)*fporf(jdime,idime)
           trace_fporf_2 = trace_fporf_2 + tkron(idime,jdime)*fporf_2(jdime,idime)
        end do
     end do

     dev_fporf = fporf - (1.0_rp/3.0_rp)*trace_fporf*tkron
     dev_fporf_2 = fporf_2 - (1.0_rp/3.0_rp)*trace_fporf_2*tkron


     ! * * * * * * * 
     ! Cauchy stress
     ! * * * * * * *
     cauchy_vol = 0.0_rp
     cauchy_iso = 0.0_rp
     cauchy_aniso1 = 0.0_rp
     cauchy_aniso2 = 0.0_rp

     do idime = 1,ndime
       do jdime = 1,ndime
         cauchy_vol(idime,jdime) = Kct*(gpdet(igaus)-1)*tkron(idime,jdime)
         cauchy_iso(idime,jdime) = a*(gpcal(idime,jdime) - (1.0_rp/3.0_rp)*tkron(idime,jdime)*i1)*gpdet(igaus)**(-5.0_rp/3.0_rp)

         cauchy_aniso1(idime,jdime) = term2*dev_fporf(idime,jdime)
         cauchy_aniso2(idime,jdime) = term2_2*dev_fporf_2(idime,jdime)
       end do
     end do

     castr = cauchy_iso + cauchy_aniso1 + flag_2*cauchy_aniso2 + cauchy_vol

     !
     !Transform Cauchy stress to S: S = J*F^-1*sigma*F^-T
     !
     do idime = 1,ndime
        do jdime = 1,ndime
           do kdime = 1,ndime
              do ldime = 1,ndime
                 gpstr(idime,jdime,igaus,1) = gpstr(idime,jdime,igaus,1)+ &
                      gpdet(igaus) *&
                      gpigd(idime,kdime,igaus)*(castr(kdime,ldime)* gpigd(jdime,ldime,igaus))
              end do
           end do
        end do
     end do



       ! * * * * * * * *
       ! Tangent moduli
       ! * * * * * * * *

       ! C = C_vol + C
       !
       ! C_vol = p_bar *(C^-1 times C^-1) - 2*J*p*I_(C^-1)
       !       
       !         p = Kct*(J-1)
       !         p_bar = J*p + J^2(dp/dJ)
       !         I_(C) = 0.5*(C_ik*C_jl + C_il*C_jk)
       !
       ! C =    

        call invmtx(gpcau(:,:,igaus),gpcai,bidon,ndime)
        invC_bar = gpcai*gpdet(igaus)**(2.0_rp/3.0_rp)

        energy_a1      = af*(i4f-1.0_rp)*exp(bf*(i4f-1.0_rp)**2)
        energy_a1_2    = af*(i6f-1.0_rp)*exp(bf*(i6f-1.0_rp)**2)
        energy_a1_a1   = af*(1.0_rp + 2.0_rp*bf*(i4f-1.0_rp)**2.0_rp)*exp(bf*(i4f-1.0_rp)**2)
        energy_a1_a1_2 = af*(1.0_rp + 2.0_rp*bf*(i6f-1.0_rp)**2.0_rp)*exp(bf*(i6f-1.0_rp)**2)
 
        Dev   = fporf_0 - (1.0_rp/3.0_rp)*i4f*invC_bar
        Dev_2 = fporf_0_2 - (1.0_rp/3.0_rp)*i6f*invC_bar  

        C4 = 0.0_rp
        shear_1 = 0.0_rp
        shear_2 = 0.0_rp
        shear_3 = 0.0_rp
        shear_4 = 0.0_rp
        term1 = 0.0_rp
        term2 = 0.0_rp
        term3 = 0.0_rp
        term4 = 0.0_rp
        term1_2 = 0.0_rp
        term2_2 = 0.0_rp
        term3_2 = 0.0_rp
        term4_2 = 0.0_rp
        bulk_1 = 0.0_rp
        bulk_2 = 0.0_rp
        bulk_3 = 0.0_rp



      if (flagt == 1_ip) then

        do idime = 1,ndime
           do jdime = 1,ndime
              do kdime = 1,ndime
                 do ldime=1,ndime

                    C4 = 0.5_rp*(gpcai(idime,kdime)*gpcai(jdime,ldime)+gpcai(jdime,kdime)*gpcai(idime,ldime)) 
 
                   
                    term1 = 4.0_rp*(gpdet(igaus)**(-4.0_rp/3.0_rp))*energy_a1_a1*Dev(idime,jdime)*Dev(kdime,ldime)

                    term2 = (4.0_rp/3.0_rp)*gpdet(igaus)*energy_a1*i4f*(C4+(1.0_rp/3.0_rp)*gpcai(idime,jdime)*gpcai(kdime,ldime))

                    term3 = -4.0_rp*energy_a1*(1.0_rp/3.0_rp)*gpdet(igaus)**(-2.0_rp/3.0_rp)*fporf_0(idime,jdime)*gpcai(kdime,ldime)

                    term4 = -4.0_rp*energy_a1*(1.0_rp/3.0_rp)*gpdet(igaus)**(-2.0_rp/3.0_rp)*fporf_0(kdime,ldime)*gpcai(idime,jdime)

                 
                    term1_2 = 4.0_rp*(gpdet(igaus)**(-4.0_rp/3.0_rp))*energy_a1_a1_2*Dev_2(idime,jdime)*Dev_2(kdime,ldime)

                    term2_2 = (4.0_rp/3.0_rp)*gpdet(igaus)*energy_a1_2*i6f*(C4+(1.0_rp/3.0_rp)*gpcai(idime,jdime)*gpcai(kdime,ldime))

                    term3_2 = -4.0_rp*energy_a1_2*(1.0_rp/3.0_rp)*gpdet(igaus)**(-2.0_rp/3.0_rp)*fporf_0_2(idime,jdime)*gpcai(kdime,ldime)

                    term4_2 = -4.0_rp*energy_a1_2*(1.0_rp/3.0_rp)*gpdet(igaus)**(-2.0_rp/3.0_rp)*fporf_0_2(kdime,ldime)*gpcai(idime,jdime)

                 
                    bulk_1 = Kct*gpdet(igaus)*(gpdet(igaus)-1.0_rp)*gpcai(idime,jdime)*gpcai(kdime,ldime)

                    bulk_2 = -2.0_rp*Kct*gpdet(igaus)*(gpdet(igaus)-1.0_rp)*C4

                    bulk_3 = Kct*(gpdet(igaus)**2.0_rp)*gpcai(idime,jdime)*gpcai(kdime,ldime)
 

                    shear_1 = (1.0_rp/3.0_rp)*i1*C4

                    shear_2 = -(1.0_rp/3.0_rp)*tkron(idime,jdime)*gpcai(kdime,ldime)

                    shear_3 = -(1.0_rp/3.0_rp)*tkron(kdime,ldime)*gpcai(idime,jdime)

                    shear_4 =  (1.0_rp/9.0_rp)*i1*gpcai(idime,jdime)*gpcai(kdime,ldime)
 
                    shear = shear_1+shear_2+shear_3+shear_4


                    gpdds(idime,jdime,kdime,ldime,igaus) = gpdds(idime,jdime,kdime,ldime,igaus)+&
                                      term1+term2+term3+term4+bulk_1+bulk_2+bulk_3+&
                                      2.0_rp*a*gpdet(igaus)**(-2.0_rp/3.0_rp)*shear+&
                                      (term1_2+term2_2+term3_2+term4_2)*flag_2
                 end do
              end do
           end do
        end do

     end if

  end do !gauss points


end subroutine sld_stress_model_137hgoc
