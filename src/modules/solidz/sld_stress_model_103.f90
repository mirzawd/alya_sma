!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_stress_model_103(pgaus,pmate,gpgdi,gpstr,gpdet,flagt,gpdds,&
             gpigd_eps,gpgdi_eps,gpdet_eps,&
             gpmof)
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_elmcla
  ! NAME
  !    sld_stress_model_103
  ! DESCRIPTION
  !    Compressible Mooney-Rivlin stress model
  !    Xiao and Belytschko's formulation
  !    Compute second Piola-Kirchoff stress tensor S_{IJ}
  !
  !    GPGDI ... Deformation tensor ...................... F = grad(phi)
  !    GPCAU ... Right Cauchy-Green deformation tensor ... C = F^t x F
  !    GPCIN ... Inverce of Right Caughy-Green tensor .... C^-1
  !    GPSTR ... 2nd P-K Stress tensor ........................... S
  !    GPPIO ... 1st Piola-Kirchhoff stress tensor ....... P = F.S
  !    GPENE ... Stored energy function .................. W
  !    FLAGT ... Flag to activate GPDDS (when implicit)
  ! USES
  ! USED BY
  !    sld_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
  use def_solidz, only       :  lawmo_sld,parco_sld
  implicit none
  integer(ip), intent(in)    :: pgaus,pmate,flagt
  real(rp),    intent(in)    :: gpmof(pgaus)
  real(rp),    intent(in)    :: gpgdi(ndime,ndime,pgaus),gpdet(pgaus)
  real(rp),    intent(out)   :: gpstr(ndime,ndime,pgaus,2)
  real(rp),    intent(out)   :: gpdds(ndime,ndime,ndime,ndime,pgaus)
  real(rp)                   :: gpcau(ndime,ndime,pgaus), gpene(pgaus)
  real(rp)                   :: gpcin(ndime,ndime,pgaus),tkron(ndime,ndime)
  integer(ip)                :: igaus,idime,jdime,kdime,i,j,k,l
  real(rp)                   :: bidon,i1,i2
  real(rp)                   :: gpcal(ndime,ndime,pgaus), dummr, gpcau2(ndime,ndime,pgaus), Iden(ndime,ndime)
!  real(rp)                   :: gpcal(ndime,ndime,pgaus), dummr, gpcau2(ndime,ndime,pgaus), gppre(pgaus), Iden(ndime,ndime)
  real(rp)                   :: cstc1,cstc2,cstd1,cstd2,s0,s1,cstm1,cstm2
  real(rp),    intent(in)        :: gpgdi_eps(ndime,ndime,pgaus)           ! Displacement Gradient F for the perturbed state
  real(rp),    intent(in)        :: gpigd_eps(ndime,ndime,pgaus)           ! Displacement Gradient Inverse F^{-1} for the perturbed state
  real(rp),    intent(in)        :: gpdet_eps(pgaus)                       ! 


  ! Mooney Rivlin's law (not implemented here)
  ! W(C) = C_1*(I_1-3) + C_2*(I_2-3) + D_1*(J-1)^2
  ! Xiao and Belytschko (compressible Mooney-Rivlin)
  ! W(C) = C_1/2*I_1 + C_2/2*I_2 - [(3/2)*C_1*I_3^{1/3}+(3/2)*C_2*I_3^{2/3}-(lambda/4)*(ln I_3)^2]

  cstd1 = parco_sld(2,pmate)  ! also called D_1, here is lambda
  cstc1 = parco_sld(3,pmate)  ! also called C_1
  cstc2 = parco_sld(4,pmate)  ! C_2=0 is the NeoHook material 

  gpstr = 0.0_rp

  do idime = 1,ndime
     do jdime = 1,ndime
        if (idime==jdime) then
           Iden(idime,jdime) = 1.0_rp
        else
           Iden(idime,jdime) = 0.0_rp
        end if
     enddo
  enddo

  ! Kronecker delta
  tkron = 0.0_rp
  do idime = 1, ndime
      tkron(idime,idime) = 1.0_rp
  end do

  !GPCAL: Left Cauchy tensor b = F F^T
  !
  do igaus=1,pgaus
     do idime = 1,ndime
        do jdime = 1,ndime
           gpcal(idime,jdime,igaus) = 0.0_rp
           do kdime = 1,ndime
              gpcal(idime,jdime,igaus) = gpcal(idime,jdime,igaus) + gpgdi(idime,kdime,igaus)*gpgdi(jdime,kdime,igaus)
           end do
        end do
     end do
  end do

  !Alternative computation
  !do igaus=1,pgaus
  !      gpcal(:,:,igaus) = matmul(gpgdi2(:,:,igaus),transpose(gpgdi2(:,:,igaus)))
  !      traceb = 0.0_rp
  !      do idime=1,ndime
  !         traceb = traceb + gpcal(idime,idime,igaus)
  !      end do
  !enddo



  do igaus=1,pgaus

     ! compute C^2 for the second invariant
     gpcau2 = 0.0_rp
     gpcau2(:,:,igaus) = matmul(gpcau(:,:,igaus),gpcau(:,:,igaus))

     call invmtx(gpcau(:,:,igaus),gpcin(:,:,igaus),bidon,ndime)

     ! * * * *
     !INVARIANTS
     ! * * * *

     !I_1=Trace(C)
     i1=0.0_rp
     i2=0.0_rp
     dummr=0.0_rp
     do idime = 1,ndime
        i1 = i1 + gpcau(idime,idime,pgaus)
        dummr = dummr + gpcau2(idime,idime,pgaus)
     end do
     i2 = 0.5_rp*(i1*i1 - dummr);
  

     select case (lawmo_sld(pmate))

       case (1_ip) !Xiao and Belytschko
          cstc1 = parco_sld(2,pmate)  ! also called C_1
          cstc2 = parco_sld(3,pmate)  ! C_2=0 is the NeoHook material 
          cstd1 = parco_sld(4,pmate)  ! also called D_1, here is lambda

         ! Energy density function
         gpene(igaus) = 0.5_rp*cstc1*i1 + 0.5_rp*cstc2*i2&
                    - (1.5_rp*cstc1*(gpdet(igaus)**(1.0_rp/3.0_rp)) + 1.5_rp*cstc2*(i1**(2.0_rp/3.0_rp))&
                    - 0.25_rp*cstd1*((log(gpdet(igaus)))**2.0_rp))

         !gppre = 0.0_rp
         !gppre(igaus) = -2.0_rp*cstd1*(gpdet(igaus)-1.0_rp) !pressure (needed in Belystchko formulation)

         ! Second Pyola-Kirchoff
         ! S = (c_1 + c_2*I_1)*Iden - c_2*C - (c_1*I_3^(1/3)+2*c_2*I_3^(2/3) - lambda*ln(I_3))*C^(-1)     
         do idime=1,ndime
           do jdime=1,ndime
             gpstr(idime,jdime,igaus,1) = (cstc1 + cstc2*i1)*Iden(idime,jdime) - cstc2*gpcau(idime,jdime,igaus)&
                                    -(cstc1*gpdet(igaus)**(1.0_rp/3.0_rp)+2.0_rp*cstc2*gpdet(igaus)**(2.0_rp/3.0_rp)&
                                    - cstd1*log(gpdet(igaus)))*gpcin(idime,jdime,igaus)
           enddo
         enddo

         if (flagt == 1_ip) then
           ! Stress trangent moduli: gpdds dSdE_{ijkl} only when implicit
           ! C_ijkl = 2*c_2*delta_ij*delta_kl - c_2*(delta_ik*delta_jl + delta_il*delta_jk) 
           !         + s_0*(C_ij^(-1)*C_jl^(-1) + C_ij^(-1)*C_jk^(-1)) - 2*s_1*C_ij^(-1)*C_kl^(-1))    
           ! s_0 = c_1*I_3^(1/3) + 2*c_2*I_2^(2/3)
           ! s_1 = ((1/3)*c_1*I_3^(-2/3) + (4/3)*c_2*(_3^(-1/3)-lambda/I_3)*I_3 
           !   
           s0 = cstc1*gpdet(igaus)**(1.0_rp/3.0_rp) + 2.0_rp*cstc2*i2**(2.0_rp/3.0_rp)
           s1 = gpdet(igaus)*((1.0_rp/3.0_rp)*cstc1*gpdet(igaus)**(-2.0_rp/3.0_rp)&
             + (4.0_rp/3.0_rp)*cstc2*gpdet(igaus)**(-1.0_rp/3.0_rp) - cstd1/gpdet(igaus))
           do i=1,ndime
             do j=1,ndime
               do k=1,ndime
                 do l=1,ndime
                   gpdds(i,j,k,l,igaus) = 2.0_rp*cstc2*tkron(i,j)*tkron(k,l)&
                                            - cstc2*(tkron(i,k)*tkron(j,l) + tkron(i,l)*tkron(j,k))&
                                            + s0*(gpcin(i,j,igaus)*gpcin(j,l,igaus) + gpcin(i,j,igaus)*gpcin(j,k,igaus))&
                                            - 2.0_rp*s1*gpcin(i,j,igaus)*gpcin(k,l,igaus)
                 enddo
               enddo
             enddo
           enddo
        end if

       case (2_ip)
         cstc1 = parco_sld(2,pmate)  ! also called C_1
         cstc2 = parco_sld(3,pmate)  ! C_2=0 is the NeoHook material 
         cstd1 = parco_sld(4,pmate)  ! also called D_1, here is lambda
         cstd2 = 2.0_rp*(cstc1+2.0_rp*cstc2)

         ! Energy density function
         ! W = d1*(I3-1)^2 - d2*log(I3) + c1*(I1-3) + c2*(I2-3)
         gpene(igaus) = cstc1*(i1-3.0_rp) + cstc2*(i2-3.0_rp) +cstd1*(gpdet(igaus)-1.0_rp)**2.0_rp - cstd2*(log(gpdet(igaus)))

         ! Second Pyola-Kirchoff
         !  
         ! Belystschko
         ! S = 2*(c1+c2*I1)*delta_ij - 2*c2*C_ij + (2*d1*I3*(I3-1)-d2)*C_ij^(-1) 
         do idime=1,ndime
           do jdime=1,ndime
             gpstr(idime,jdime,igaus,1) = 2.0_rp*(cstc1 + cstc2*i1)*Iden(idime,jdime) - 2.0_rp*cstc2*gpcau(idime,jdime,igaus)&
                                       +2.0_rp*(2.0_rp*cstd1*gpdet(igaus)*(gpdet(igaus)-1.0_rp)-cstd2)*gpcin(idime,jdime,igaus)        
           end do
         end do

         if (flagt == 1_ip) then

         ! Stress trangent moduli: gpdds dSdE_{ijkl} only when implicit
         ! C_ijkl = 2*c_2*(delta_ik*delta_jl + delta_il*delta_jk)
         !         + d1*I3*(2*I3-1)*(C_ij^(-1)*C_jl^(-1) + C_ij^(-1)*C_jk^(-1)) 
         !         - (2*d1*I3*(I3-1)-d2)*(delta_ik*delta_jl + delta_il*delta_jk)   
           do i=1,ndime
             do j=1,ndime
               do k=1,ndime
                 do l=1,ndime
                   gpdds(i,j,k,l,igaus) = 2.0_rp*cstc2*(tkron(i,k)*tkron(j,l)+tkron(i,l)*tkron(j,k))&
                                              + cstd1*gpdet(igaus)*((2.0_rp*gpdet(igaus)-1.0_rp)&
                                               *(gpcin(i,k,igaus)*gpcin(j,l,igaus)+gpcin(i,l,igaus)*gpcin(j,k,igaus))&
                                              - (2.0_rp*cstd1*gpdet(igaus)*(gpdet(igaus)-1.0_rp)-cstd2)*(tkron(i,k)*tkron(j,l)+tkron(i,l)*tkron(j,k)))
                 enddo
               enddo
             enddo
           enddo
         endif


       case (3_ip)  !Modified Mooney Rivlin
         cstc1 = parco_sld(2,pmate)  ! also called C_1
         cstc2 = parco_sld(3,pmate)  ! C_2=0 is the NeoHook material 
         cstm1 = parco_sld(4,pmate)
         cstm2 = parco_sld(5,pmate)
         cstd1 = parco_sld(6,pmate)  ! also called D_1, here is lambda

         ! Energy density function
         ! W = 
         gpene(igaus) = cstc1*(i1-3.0_rp) + cstc2*(i2-3.0_rp) &
                    + cstm1*(exp(cstm2*(i1-3.0_rp))-1.0_rp)
                    !& - 0.25_rp*cstd1*((log(gpdet(igaus)))**2.0_rp))

         ! S = 
         do idime=1,ndime
           do jdime=1,ndime
             gpstr(idime,jdime,igaus,1) = (cstc1 + cstc2*i1)*Iden(idime,jdime) - cstc2*gpcau(idime,jdime,igaus) &
                                     + 2.0_rp*cstm1*cstm2*(exp(cstm2*(i1-3.0_rp)))* Iden(idime,jdime)
                                      ! +(2.0_rp/cstd1)*(gpdet(igaus)-1)*gpdet(igaus)*gpcin(idime,jdime,igaus)
           end do
         end do

         if (flagt == 1_ip) then

           ! Stress trangent moduli: gpdds dSdE_{ijkl} only when implicit
           ! C_ijkl = MODIFICAR NO ES ESTA  2*c_2*delta_ij*delta_kl - c_2*(delta_ik*delta_jl + delta_il*delta_jk) 
           !         + s_0*(C_ij^(-1)*C_jl^(-1) + C_ij^(-1)*C_jk^(-1)) - 2*s_1*C_ij^(-1)*C_kl^(-1))    

           do i=1,ndime
             do j=1,ndime
               do k=1,ndime
                 do l=1,ndime
                   gpdds(i,j,k,l,igaus) =  2.0_rp*cstc2*(tkron(i,k)*tkron(j,l)+tkron(i,l)*tkron(j,k)) &
                      + 2.0_rp*cstm1*cstm2*cstm2*(exp(cstm2*(i1-3.0_rp)))*(tkron(i,k)*tkron(j,l)+tkron(i,l)*tkron(j,k))
                      !+ 2.0_rp*cstd1*gpdet(igaus)*((2.0_rp*gpdet(igaus)-1.0_rp)&
                      !*(gpcin(i,k,igaus)*gpcin(j,l,igaus)+gpcin(i,l,igaus)*gpcin(j,k,igaus))&
                      !- 0.5_rp*(gpdet(igaus)-1.0_rp)*(tkron(i,k)*tkron(j,l)+tkron(i,l)*tkron(j,k)))
                 enddo
               enddo
             enddo
           enddo
         endif

       case default

     end select 
    
  enddo

end subroutine sld_stress_model_103
