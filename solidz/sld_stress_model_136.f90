!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!$-----------------------------------------------------------------------
!> @addtogroup SolidzMaterials
!> @{
!> @file    sld_stress_model_136.f90
!> @author  
!> @date    
!> @brief   Guccione
!> @details 
!> @} 
!-----------------------------------------------------------------------
subroutine sld_stress_model_136(&
     pgaus,pmate,gpcau,gpgdi,gpigd,gpstr,gpdet,&
     nfibe,ielem,elcod,pnode,lnods,gpsha,flagt,gpdds,&
             gpigd_eps,gpgdi_eps,gpdet_eps,&
             gpmof)
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_stress_model_136
  ! NAME
  !    sld__stress_model_136
  ! DESCRIPTION
  !    Guccione
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
!  use def_master, only       :  ittim
  
  
  implicit none
  integer(ip), intent(in)    :: pgaus,pmate,pnode,lnods(pnode),flagt
  real(rp),    intent(in)    :: gpcau(ndime,ndime,pgaus),gpdet(pgaus),gpgdi(ndime,ndime,pgaus),gpmof(pgaus)
  real(rp),    intent(in)    :: gpigd(ndime,ndime,pgaus),gpsha(pnode,pgaus)
  real(rp),    intent(out)   :: gpstr(ndime,ndime,pgaus,2),nfibe(ndime,pgaus)
  real(rp),    intent(out)   :: gpdds(ndime,ndime,ndime,ndime,pgaus)
  real(rp),    intent(in)        :: gpgdi_eps(ndime,ndime,pgaus)           ! Displacement Gradient F for the perturbed state
  real(rp),    intent(in)        :: gpigd_eps(ndime,ndime,pgaus)           ! Displacement Gradient Inverse F^{-1} for the perturbed state
  real(rp),    intent(in)        :: gpdet_eps(pgaus)                       ! 
  real(rp)                   :: nfibe_length
  
  integer(ip)                :: igaus,idime,jdime,kdime,ldime,ielem,inode,ipoin
  real(rp)                   :: gpcal(ndime,ndime),castr(ndime,ndime),castv(ndime,ndime),sleci(3,3,3)
  
  real(rp)                   :: nfibe0(ndime,pgaus),norma0(ndime,pgaus),nshet0(ndime,pgaus) !f0, s0, n0 (ref. config)
  real(rp)                   ::                     norma(ndime,pgaus) ,nshet(ndime,pgaus)  !f,s,n (fF as define in paper) (nfibe is "output")
  real(rp)                   :: nfibt(ndime,pgaus) ,normt(ndime,pgaus) ,nshtt(ndime,pgaus)  !f,s,n (updated-Unit)
  
  real(rp)                   :: c1, c2, c3, c4, Q
  
  real(rp)                   :: i1,i4f,i4s,i4n,i8fs,i8sf,i8sn,i8ns,i8nf,i8fn 
  real(rp)                   :: fporf(ndime,ndime),spors(ndime,ndime),nporn(ndime,ndime)
  real(rp)                   :: fpors(ndime,ndime),sporf(ndime,ndime),fporn(ndime,ndime),nporf(ndime,ndime),sporn(ndime,ndime),npors(ndime,ndime)
  real(rp)                   :: term1,term2,term3,term4,term5,term6
!  real(rp)                   :: termK
  real(rp)                   :: gpcai(ndime,ndime),bidon,Kct
  
  real(rp)                   :: lamda(3)
  real(rp)                   :: ovlam(3)
  
  real(rp)                   :: castr_active(ndime,ndime),gpstr_active(ndime,ndime)
  
  real(rp)                   :: elfib(3,mnode),elfis(3,mnode),elfin(3,mnode),elcac(mnode)

  real(rp),    intent(in)    :: elcod(ndime,mnode)
  
  real(rp)                   :: tkron(3,3)


 Kct= parco_sld( 2,pmate)
 c1 = parco_sld( 3,pmate)
 c2 = parco_sld( 4,pmate)
 c3 = parco_sld( 5,pmate)
 c4 = parco_sld( 6,pmate)
 
  
  !Now in sld_velsnd.f90
  ! Compute sound velocity only at time-step 0 (initialization)
  !if (ittim == 0_ip) then
  !   velas_sld(1,1) = sqrt(parco_sld(1,1)/densi_sld(1,1))
  !end if
  
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
  
  tkron = 0.0_rp
  do idime=1,ndime
     tkron(idime,idime) = 1.0_rp
  end do


  ! levi-civita symbol
  sleci= 0.0_rp
  sleci(1,2,3)= 1.0_rp
  sleci(2,3,1)= 1.0_rp
  sleci(3,1,2)= 1.0_rp
  sleci(3,2,1)= 1.0_rp
  sleci(1,3,2)= 1.0_rp
  sleci(2,1,3)= 1.0_rp

  !
  ! Gather
  !
  if( modfi_sld(pmate) .ne. 0 ) then
     elfis = 0.0_rp           
     elfin = 0.0_rp
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
  elcac = 0.0_rp
  
  do igaus = 1,pgaus
     
     !Initialize active, passive and volumetric stresses to zero
     castr = 0.0_rp
     castv = 0.0_rp
     castr_active = 0.0_rp
     gpstr_active = 0.0_rp

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

     
     ! * * * * * *
     ! fiber from file 
     ! * * * * * *
     if( modfi_sld(pmate) == 5 ) then !NOT USE ANYMORE???
        
        nfibe0(1,igaus) = fiber(1,lnods(1))
        nfibe0(2,igaus) = fiber(2,lnods(1))
        nfibe0(3,igaus) = fiber(3,lnods(1))
        norma0 = 0.0_rp
        nshet0 = 0.0_rp
        
     else if( modfi_sld(pmate) < 0 ) then
        !
        ! Fibers interpolated on Gauss point
        !
        nfibe0(1,igaus) = 0.0_rp
        nfibe0(2,igaus) = 0.0_rp
        nfibe0(3,igaus) = 0.0_rp
        do inode = 1,pnode
           nfibe0(1,igaus) = nfibe0(1,igaus) + gpsha(inode,igaus) * elfib(1,inode)
           nfibe0(2,igaus) = nfibe0(2,igaus) + gpsha(inode,igaus) * elfib(2,inode)
           nfibe0(3,igaus) = nfibe0(3,igaus) + gpsha(inode,igaus) * elfib(3,inode)
        end do
        
        
        !Normalized f_0
        bidon=sqrt(nfibe0(1,igaus)**2.0_rp+nfibe0(2,igaus)**2.0_rp+nfibe0(3,igaus)**2.0_rp)
        if (bidon == 0.0_rp) bidon= 1.0_rp
        do idime=1,ndime
           nfibe0(idime,igaus)=nfibe0(idime,igaus)/bidon
        end do
        
        ! * check if unitary 
        !        if (abs(bidon-1.0)>0.001) then
        !           write(*,*) 'bidon ',bidon
        !           write(*,*) '--| ALYA  sld_stress_model_136. Fibers not unitary: ',bidon
        !           stop
        !        end if
        ! *
        !if (ielem==20) write(*,*) bidon      

        
     else !for the prolate spheroid
        
        !call sld_angleh(&
        !     modfi_sld(pmate),pgaus,pmate,igaus,ielem,elcod,nfibe0,norma0,nshet0,pnode,lnods,gpsha)
        
        !To create the fibers file
        call sld_angle2(1,modfi_sld(pmate))
        write(*,*) '--| ALYA  sld_stress_model_136. Fibers produced in fort.7777'
        stop

     end if
     
     if( kfl_fiber_sld == 2 .or. kfl_fiber_sld == 3) then
        
        !
        ! Orthotropic material
        !
        norma0(1,igaus) = 0.0_rp
        norma0(2,igaus) = 0.0_rp
        norma0(3,igaus) = 0.0_rp
        nshet0(1,igaus) = 0.0_rp
        nshet0(2,igaus) = 0.0_rp
        nshet0(3,igaus) = 0.0_rp
        do inode = 1,pnode
           norma0(1,igaus) = norma0(1,igaus) + gpsha(inode,igaus) * elfin(1,inode)
           norma0(2,igaus) = norma0(2,igaus) + gpsha(inode,igaus) * elfin(2,inode)
           norma0(3,igaus) = norma0(3,igaus) + gpsha(inode,igaus) * elfin(3,inode)
           nshet0(1,igaus) = nshet0(1,igaus) + gpsha(inode,igaus) * elfis(1,inode)
           nshet0(2,igaus) = nshet0(2,igaus) + gpsha(inode,igaus) * elfis(2,inode)
           nshet0(3,igaus) = nshet0(3,igaus) + gpsha(inode,igaus) * elfis(3,inode)
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
     
     i4f  = 0.0_rp
     i4s  = 0.0_rp
     i4n  = 0.0_rp
     i8fs = 0.0_rp
     i8sf = 0.0_rp
     i8ns = 0.0_rp
     i8sn = 0.0_rp
     i8fn = 0.0_rp
     i8nf = 0.0_rp
     do idime = 1,ndime
        do jdime = 1,ndime
           i4f  =  i4f + nfibe0(idime,igaus) * gpcau(idime,jdime,igaus) * nfibe0(jdime,igaus) ! I_4f  = f0 * C_ij * f0 (eq. 5.1)
           i4s  =  i4s + nshet0(idime,igaus) * gpcau(idime,jdime,igaus) * nshet0(jdime,igaus) ! I_4s  = s0 * C_ij * s0
           i4n  =  i4n + norma0(idime,igaus) * gpcau(idime,jdime,igaus) * norma0(jdime,igaus) ! I_4n  = n0 * C_ij * n0
           i8fs = i8fs + nfibe0(idime,igaus) * gpcau(idime,jdime,igaus) * nshet0(jdime,igaus) ! I_8fs = f0 * C_ij * s0 (eq. 5.3)
           i8sf = i8sf + nshet0(idime,igaus) * gpcau(idime,jdime,igaus) * nfibe0(jdime,igaus) ! I_8sn = s0 * C_ij * f0
           i8ns = i8ns + norma0(idime,igaus) * gpcau(idime,jdime,igaus) * nshet0(jdime,igaus) ! I_8ns = n0 * C_ij * s0
           i8sn = i8sn + nshet0(idime,igaus) * gpcau(idime,jdime,igaus) * norma0(jdime,igaus) ! I_8sn = s0 * C_ij * n0
           i8fn = i8fn + nfibe0(idime,igaus) * gpcau(idime,jdime,igaus) * norma0(jdime,igaus) ! I_8sn = f0 * C_ij * n0
           i8nf = i8nf + norma0(idime,igaus) * gpcau(idime,jdime,igaus) * nfibe0(jdime,igaus) ! I_8sn = n0 * C_ij * f0
           
           ! * * * * *
           !Transform f=F'*f0 ; s=F'*s0 ; n=F'*n0 (PUSH FORWARD)
           
           nfibe(idime,igaus) = nfibe(idime,igaus) + gpgdi(idime,jdime,igaus) * nfibe0(jdime,igaus)
           nshet(idime,igaus) = nshet(idime,igaus) + gpgdi(idime,jdime,igaus) * nshet0(jdime,igaus)
           norma(idime,igaus) = norma(idime,igaus) + gpgdi(idime,jdime,igaus) * norma0(jdime,igaus)
        end do
     end do
     

     ! Store nfibe for postproc

     do idime = 1,ndime
        gpfib(idime,igaus,ielem) = nfibe(idime,igaus)
     end do

     !
     ! f x f , s x s, f x s, s x f (outer products [3x3])
     !
     fporf = 0.0_rp
     spors = 0.0_rp
     nporn = 0.0_rp
     fpors = 0.0_rp
     sporf = 0.0_rp
     npors = 0.0_rp
     sporn = 0.0_rp
     nporf = 0.0_rp
     fporn = 0.0_rp
     do idime = 1,ndime
        do jdime = 1,ndime
           fporf(idime,jdime) = nfibe(idime,igaus) * nfibe(jdime,igaus)
           spors(idime,jdime) = nshet(idime,igaus) * nshet(jdime,igaus)
           nporn(idime,jdime) = norma(idime,igaus) * norma(jdime,igaus)
           fpors(idime,jdime) = nfibe(idime,igaus) * nshet(jdime,igaus)
           sporf(idime,jdime) = nshet(idime,igaus) * nfibe(jdime,igaus)
           npors(idime,jdime) = norma(idime,igaus) * nshet(jdime,igaus)
           sporn(idime,jdime) = nshet(idime,igaus) * norma(jdime,igaus)
           nporf(idime,jdime) = norma(idime,igaus) * nfibe(jdime,igaus)
           fporn(idime,jdime) = nfibe(idime,igaus) * norma(jdime,igaus)         
        end do
     end do
     
     ! * * * * * *
     ! Passive terms (calculate directy the cauchy stress (castr))
     ! * * * * * *

     
     ! Non-collageous material
     !muscle fibers in tension
     term1 = c2*0.5_rp*(i4f-1.0_rp)        
     !collagen fibers in tension
     term2 = c3*0.5_rp*(i4s-1.0_rp) 
    
     term3 = c3*0.5_rp*(i4n-1.0_rp) 
     !shear
     term4 = c3*(i8sn)
     term5 = c4*0.5_rp*(i8fs+i8sf)
     term6 = c4*0.5_rp*(i8fn+i8nf) 
     
     Q = c2*(0.5_rp*(i4f-1.0_rp))**2.0_rp +c3*(0.5_rp*(i4s-1.0_rp))**2.0_rp + c3*(0.5_rp*(i4n-1.0_rp))**2.0_rp &
          +2.0_rp*c3*(0.5_rp*i8sn)**2.0_rp + 2.0_rp*c4*(0.25_rp*i8sf*i8fs) + 2.0_rp*c4*(0.25_rp*i8nf*i8fn)
     
     do idime = 1,ndime
        do jdime = 1,ndime
           castr(idime,jdime) = castr(idime,jdime)+(1/gpdet(igaus))*2.0_rp*c1*exp(Q)&   
                *(term1 * fporf(idime,jdime)&         
                + term2 * spors(idime,jdime)&         
                + term3 * nporn(idime,jdime)&
                + term4 *(npors(idime,jdime)+fporn(idime, jdime))&
                + term5 * (fpors(idime,jdime)+sporf(idime,jdime)))
        end do
     end do
     
     
     
     ! * * * * * *
     ! Volumetric term
     ! * * * * * *
    
     !term1 = 1/gpdet(igaus)* Kct * gpmof(igaus)*log(gpdet(igaus)) 
     
     !do idime = 1,ndime   
     !   castv(idime,idime) = castv(idime,idime)+term1  
     !end do
    
     !Calculate the stretch (lamda). In the fiber direction: lambda^2=f_0*C*f_0 = I_4f
     
     lamda(1) = sqrt(i4f)
     lamda(2) = sqrt(i4s)
     lamda(3) = sqrt(i4n)
     
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

     !
     ! Updated fiber direction: lambda_f * f_true = F * f_0 = f
     ! => here f is the UNIT fiber direction in the current configuration
     ! used for active stress and postproc
     !
     nfibe_length= 0.0_rp
     do idime = 1,ndime
        nfibt(idime,igaus) = nfibe(idime,igaus) * ovlam(1)
        nshtt(idime,igaus) = nshet(idime,igaus) * ovlam(2)
        normt(idime,igaus) = norma(idime,igaus) * ovlam(3)
        nfibe_length=   nfibe_length + nfibe(idime,igaus)*nfibe(idime,igaus) 
     end do
     nfibe_length= sqrt(nfibe_length) 

     ! Electromechanical coupling
     ! * * * * * *
     ! Active stress T(lambda, [Ca++])m x m  (Hunter p 688)
     ! * * * * * *
     castr_active = 0.0_rp

     !Transform Cauchy stress to S
     do idime = 1,ndime
        do jdime = 1,ndime
           do kdime = 1,ndime
              do ldime = 1,ndime
                 gpstr(idime,jdime,igaus,1) = gpstr(idime,jdime,igaus,1)+ &
                      gpdet(igaus) *&
                      gpigd(idime,kdime,igaus) * &
                      (castr(kdime,ldime) + castv(kdime,ldime) + castr_active(kdime,ldime))*&
                      gpigd(jdime,ldime,igaus)
              end do
           end do
        end do
     end do
     
     if (flagt == 1_ip) then
        ! Tangent moduli d2W/dE_ij*dE_kl
        call invmtx(gpcau(:,:,igaus),gpcai,bidon,ndime)
        do idime = 1,ndime
           do jdime = 1,ndime
              do kdime = 1,ndime
                 do ldime=1,ndime
                    term1 = c1*c2*exp(Q)*(0.5_rp+c1*c2*((i4f-1)*0.5_rp)**2.0_rp)&
                         *fporf(idime,jdime)*fporf(kdime,ldime)
                    term2 = c1*c3*exp(Q)*(0.5_rp+c1*c3*((i4s-1)*0.5_rp)**2.0_rp)&
                         *spors(idime,jdime)*spors(kdime,ldime)
                    term3 = c1**2.0_rp*c2*c3*exp(Q)*(i4f-1.0_rp)*(i4s-1.0_rp)*0.25_rp&
                         *(spors(idime,jdime)*fporf(kdime,ldime)+fporf(idime,jdime)*spors(kdime,ldime)) 
                    term4 =c1*c4*exp(Q)*(1.0_rp+c1*c4*i8fs)&
                         *(sporf(idime,jdime)+fpors(idime,jdime))*(sporf(kdime,ldime)+fpors(kdime,ldime))*4.0_rp 
                    !termK = Kct*gpmof(igaus)*(gpcai(idime,jdime)*gpcai(kdime,ldime)-log(gpdet(igaus))&  
                    !     *(gpcai(idime,kdime)*gpcai(jdime,ldime)+gpcai(idime,ldime)*gpcai(jdime,kdime)))  
                    gpdds(idime,jdime,kdime,ldime,igaus) = &
                         gpdds(idime,jdime,kdime,ldime,igaus)+term1+term2+term3+term4!+termK
                  end do
               end do
           end do
        end do
     end if
  end do !gauss points
  
end subroutine sld_stress_model_136
