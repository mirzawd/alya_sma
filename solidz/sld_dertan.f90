!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!$---------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_dertan.f90
!> @author  Eva Casoni
!> @date    27/11/2014
!> @brief   Numerical tangent moduli
!> @details Computation of the numerical tangent moduli using FD
!> @}
!-----------------------------------------------------------------------
subroutine sld_dertan(pgaus,pmate,temp0,temp1,&
        gpgi0,gpgi1,gpigd,gpcau,gpdet,gprat,gppio,gpstr,&
        gpene,flagt,gptmo,flags,ggsgo,ielem,nfibe,elcod,pnode,lnods,gpsha,gpdds,gpmof)

!-----------------------------------------------------------------------
!****f* Solidz/sld_builtin_materials
! NAME
!    sld_builtin_materials
! DESCRIPTION
!    Interface for calling built-in subroutines
! INPUT
!    PGAUS ... Number of Gauss points
!    PMATE ... Material number
!    TEMP0 ... Previous temperature ............................ temp(n)
!    TEMP1 ... Updated temperature ............................. temp(n+1)
!    GPGI0 ... Previous deformation gradient tensor ............ F(n)
!    GPGI1 ... Updated deformation gradient tensor ............. F(n+1)
!    GPIGD ... Inverse of updated deformation gradient tensor .. F^{-1}(n+1)
!    GPCAU ... Updated right Cauchy-Green tensor ............... C(n+1)
!    GPDET ... Updated jacobian ................................ J = det(F(n+1))
!    GPRAT ... Rate of Deformation tensor ...................... Fdot = grad(phidot)
!    FLAGT ... Flag for tangent moduli calculations ............ 1 if yes; 0 if no
!    FLAGS ... Flag for strain gradient calculations ........... 1 if yes; 0 if no
! INPUT/OUTPUT
!    GPPIO ... 1st Piola-Kirchhoff stress tensor at t(n)/t(n+1)  P(n)/P(n+1)
!    GPSTR ... 2nd Piola-Kirchhoff stress tensor at t(n)/t(n+1)  S(n)/S(n+1)
!    GPIVO ... Internal variable array at t(n)/t(n+1) .......... Q(n)/Q(n+1)
!    GPENE ... Stored energy function .......................... W
! OUTPUT
!    GPTMO ... Tangent moduli at t(n+1) ........................ dP/dF(n+1)
!    GGSGO ... Deformation gradient tensor ............ dF/dX(n+1)
! USES
!    sld_stress_model_xxx ................... constitutive law xxx
!    (other subroutines to be added as necessary)
! USED BY
!    sld_elmcla
!***
!-----------------------------------------------------------------------

  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mnode
  use def_solidz
  use mod_sld_stress_model_comput, only : sm_tensor_to_voigt_fourth

  implicit none
  integer(ip), intent(in)    :: pgaus,pmate,ielem,pnode,lnods(pnode)
  real(rp),    intent(in)    :: temp0(pgaus)
  real(rp),    intent(in)    :: temp1(pgaus)
  real(rp),    intent(in)    :: gpgi0(ndime,ndime,pgaus)
  real(rp),    intent(in)    :: gpgi1(ndime,ndime,pgaus)
  real(rp),    intent(in)    :: gpcau(ndime,ndime,pgaus)
  real(rp),    intent(in)    :: gpdet(pgaus),gpmof(pgaus)
  real(rp),    intent(in)    :: gprat(ndime,ndime,pgaus)
  integer(ip), intent(in)    :: flagt,flags
  real(rp),    intent(in)    :: gpigd(ndime,ndime,pgaus)
  real(rp),    intent(inout) :: gppio(ndime,ndime,pgaus)
  real(rp),    intent(inout) :: gpstr(ndime,ndime,pgaus)
  real(rp),    intent(out)   :: nfibe(ndime,pgaus)
  real(rp),    intent(inout) :: gpene(pgaus)
  real(rp),    intent(out)   :: gptmo(ndime,ndime,ndime,ndime,pgaus)
  real(rp),    intent(out)   :: ggsgo(ndime,ndime)
  real(rp),    intent(out)   :: gpdds(ndime,ndime,ndime,ndime,pgaus)
  integer(ip)                :: igaus,idime,jdime
  real(rp)                   :: tkron(ndime,ndime)
  real(rp),    intent(in)    :: elcod(ndime,mnode)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)

  real(rp)                   :: eps, epsil, dcau(ndime,ndime), pertcau(ndime,ndime,pgaus), a, b, errfb
  real(rp)                   :: eigva(ndime), eigve(ndime,ndime), Up(ndime,ndime), U(ndime,ndime), Uaux(ndime,ndime), Uval(ndime,ndime), R(ndime,ndime), Uinv(ndime,ndime)
  real(rp)                   :: Fpert(ndime,ndime,pgaus), Fpinv(ndime,ndime,pgaus), Fpdet(pgaus), gpddsvo(2*ndime,2*ndime),  pertdds(ndime,ndime,ndime,ndime,pgaus), ednum(2*ndime,2*ndime)
  real(rp)                   :: Spert(ndime,ndime,pgaus), Svoi(2*ndime,pgaus), Spvoi(2*ndime,pgaus), ddnum(2*ndime,2*ndime,pgaus)
  integer(ip)                :: ii,jj,idumm,ivoig
  integer(ip)  :: flgs1,flgs3

  logical(lg)                :: debuggingMode

  debuggingMode = .false.

  ! declaration of flags for additional inputs/outputs

  flgs1 = 0_ip ! flag for 1st P-K Stress tensor from 2nd P-K in output P(n+1)
  flgs3 = 0_ip ! flag for tangent moduli from material stress tangent  dP/dF

  ! touch the values of the "output" so that the compiler does think it is used
  ! (to be erased as soon as there exists one subroutine that uses it)

  ggsgo(1,1)       = 0.0_rp
  gptmo(1,1,1,1,1) = 0.0_rp
  gpddsvo = 0.0_rp

  Fpdet = 0.0_rp   ! Determinant of F perturbate
  Fpert = 0.0_rp   ! F perturbate
  Fpinv = 0.0_rp   ! Inverse od F perturbate

  U = 0.0_rp
  Uaux = 0.0_rp
  Uval = 0.0_rp
  Uinv = 0.0_rp
  R = 0.0_rp

  Spvoi = 0.0_rp
  Svoi = 0.0_rp

  ddnum = 0.0_rp   ! Numerical tangent moduli

  debuggingMode = .false. ! (ielem == 1_ip)

  ! Kronecker delta

  tkron = 0.0_rp
  do idime = 1, ndime
      tkron(idime,idime) = 1.0_rp
  end do

  eps = 0.00000000000001_rp
  epsil = 0.0_rp

  if (lawst_sld(pmate)==100) then                                         ! isolinear elastic
        flgs1 = 1_ip
        flgs3 = 1_ip
  else if (lawst_sld(pmate)==101) then                                    ! Neo-Hookean belytschko (textbook)
        flgs1 = 1_ip
        flgs3 = 1_ip
  else if (lawst_sld(pmate)==103) then                                    ! Mooney-Rivlin
        flgs1 = 1_ip
        flgs3 = 1_ip
  else if (lawst_sld(pmate)==134) then                                    ! Holza & Ogden
        flgs1 = 1_ip
        flgs3 = 1_ip
  end if

  ! call built-in material constitutive laws

  if (lawst_sld(pmate)==100) then                                         ! isolinear elastic
     call sld_stress_model_100(pgaus,pmate,gpgi1,gpstr,gpcau,gpdet,gpigd,ielem,elcod,flagt,gpdds,gpmof)    !with activation
  else if (lawst_sld(pmate)==101) then                                    ! Neo-Hookean belytschko (textbook)
     call sld_stress_model_101(pgaus,pmate,gpcau,gpgi1,gpene,gpstr,gpdet,flagt,gpdds,gpmof)
  else if (lawst_sld(pmate)==103) then
       call sld_stress_model_103(pgaus,pmate,gpgi1,gpstr,gpdet,flagt,gpdds,gpmof) ! Mooney-Rivlin (Belytschko notation)
  end if

  ! Tensor to voigt notation of gpdds  (to print gpdds)
  call SM_tensor_to_voigt_fourth(ndime, 0_ip, gpdds(:,:,:,:,1), gpddsvo(:,:))

 do igaus=1,pgaus

    ! Compute useful quantities (for polar decomposition: F=RU)
    call spcdec(gpcau(:,:,igaus),eigva,eigve,idumm,1_ip,'SLD_DERTAN')
    do idime=1,ndime
           Uval(idime,idime) = sqrt(eigva(idime))
    end do
    Uaux  = matmul(eigve,Uval)
    U = matmul (Uaux, transpose(eigve))

    call invmtx(U(1,1),Uinv(1,1),gpdet(igaus),ndime)
    R = matmul(gpgi1(:,:,igaus),Uinv)

    ! Perturb right Cauchy-Green deformation tensor to generate columns of material stiffness
    do ii=1,6
       dcau = 0.0_rp
       if (ii<4_ip) then
          if (abs(gpcau(ii,ii,igaus))<eps) then
             epsil = eps
          else
             epsil = sqrt(eps)*gpcau(ii,ii,igaus)
          end if
          dcau(ii,ii) = 1.0_rp
       else
          if (abs(gpcau(nvgij_sld(ii,1), nvgij_sld(ii,2), igaus))<eps) then
             epsil = eps
          else
             epsil = sqrt(eps)*gpcau(nvgij_sld(ii,1), nvgij_sld(ii,2),igaus)
          end if
          dcau(nvgij_sld(ii,1), nvgij_sld(ii,2)) = 0.5_rp
          dcau(nvgij_sld(ii,2), nvgij_sld(ii,1)) = 0.5_rp
       end if

       ! Perturbed right Cauchy-Green Deformation tensor
       pertcau(:,:,igaus) = 0.0_rp
       do idime = 1,ndime
          do jdime = 1,ndime
             pertcau(idime,jdime,igaus) = gpcau(idime,jdime,igaus) + epsil*dcau(idime,jdime)
          end do
       end do

       ! Related perturbed kinematic measures (for polar decomposition Fp=RUp)
       call spcdec(pertcau,eigva,eigve,idumm,1_ip,'SLD_DERTAN, PERTURBED KINEMATIC MEASURES')
       do idime=1,ndime
          Uval(idime,idime) = sqrt(eigva(idime))
       end do
       Uaux  = matmul(eigve,Uval)
       Up = matmul (Uaux, transpose(eigve))

       ! Perturbed gradient deformation
       Fpert(:,:,igaus) = matmul(R,Up)
       call invmtx(Fpert(1,1,igaus),Fpinv(1,1,igaus),Fpdet(igaus),ndime)

       ! Compute pertrubed stresses from the model
       if (lawst_sld(pmate)==100) then                                         ! isolinear elastic
          call sld_stress_model_100(pgaus,pmate,Fpert,Spert,pertcau,Fpdet,Fpinv,ielem,elcod,flagt,pertdds,gpmof)    !with activation
       else if (lawst_sld(pmate)==101) then                                    ! Neo-Hookean belytschko (textbook)
          call sld_stress_model_101(pgaus,pmate,pertcau,Fpert,gpene,Spert,Fpdet,flagt,pertdds,gpmof)
       else if (lawst_sld(pmate)==103) then
          call sld_stress_model_103(pgaus,pmate,Fpert,Spert,Fpdet,flagt,pertdds,gpmof) ! Mooney-Rivlin (Belytschko notation)
       end if

       ! Transform from matrix to Voig notation
      do idime = 1, ndime
         do jdime = 1, ndime
            ivoig = nvgij_inv_sld(idime, jdime)
            Spvoi(ivoig,igaus) = Spert(idime,jdime,igaus)
            Svoi(ivoig,igaus) = gpstr(idime,jdime,igaus)
        end do
      end do

     ! Column ii of material stiffness
     do jj=1,6
         ddnum(jj,ii,igaus) = (2.0_rp/epsil)*(Spvoi(jj,igaus) - Svoi(jj,igaus))
     end do

   end do  !end do of perturbations

end do ! fin del bucle ii de las columnas de tangent moduli

if ( debuggingMode ) then
! Check if the numerical tangent is equal to the analytical
  print*, 'GPDDS VOIGT'
  do ii=1,6
    ! print*, gpddsvo(1,ii),gpddsvo(2,ii),gpddsvo(3,ii),gpddsvo(4,ii),gpddsvo(5,ii),gpddsvo(6,ii)
    print*, 'column:', ii
    print*, gpddsvo(1,ii),gpddsvo(2,ii),gpddsvo(3,ii),gpddsvo(4,ii),gpddsvo(5,ii),gpddsvo(6,ii)
  end do
  print*, ' '
  print*, 'GPDDS numerical'
  do ii=1,6
     print*, 'ddnum: column', ii
     print*, ddnum(1,ii,1), ddnum(2,ii,1), ddnum(3,ii,1), ddnum(4,ii,1), ddnum(5,ii,1), ddnum(6,ii,1)
  end do

  print*, 'Absolute error '
  do ii=1,6
    print*, abs(ddnum(:,ii,pgaus)-gpddsvo(:,ii))
  end do

    print*, 'Relative error '
  do ii=1,6
    do jj = 1,6
    if (abs(gpddsvo(jj,ii)- ddnum(jj,ii,pgaus))< eps ) then
        ednum(jj,ii)= 0.0_rp
    else
        ednum(jj,ii)= abs(gpddsvo(jj,ii)- ddnum(jj,ii,pgaus))/abs(gpddsvo(jj,ii))
    end if
    end do
  end do

  do ii=1,6
    print*, ednum(:,ii)
  end do

  do ii = 1,6
    do jj = 1, 6
       a = a + ednum(jj,ii)
       b = b + gpddsvo(jj,ii)
       errfb = a/b
    end do
  end do

  print*, 'Frobenius norm:', errfb
! Transform gpdds computed in the model to voigt notation and compare it with ddnum
end if

end subroutine sld_dertan
