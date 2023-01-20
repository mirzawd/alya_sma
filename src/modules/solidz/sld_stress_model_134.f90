!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup SolidzMaterials
!> @{
!> @file    sld_stress_model_134.f90
!> @author  Eva Casoni 
!> @author  Francesc Levrero-Florencio
!> @author  Adria Quintanas-Corominas 
!> @author  Alfonso Santiago
!> @date    01/08/2017
!> @brief   Law of Holzapfel et Ogden modified for compressible anisotropic (Nolan, 2014)
!> @details 20/04/2020 : Adria  Quintanas-Corominas and Alfonso Santiago move the activa part
!>                       computation outside of the model
!>          25/03/2020 : Adria Quintanas-Corominas and Alfonso Santiago cleaning due to the
!>                       migration of the electro-mechanical coupling.
!>          01/08/2017 : Eva Casoni and Francesc Levrero-Florencio implementation
!> @} 
!-----------------------------------------------------------------------
subroutine sld_stress_model_134( &
    pgaus,pmate,flagt,gpdet,gpcau,gpbasis,gpstr,gpdds)
    ! -----------------------------------
    use def_kintyp,                only :  ip, rp
    use def_solidz,                only :  parco_sld
    use mod_maths_basic,           only :  maths_MULT_VxMxV
    ! -----------------------------------
    implicit none
    ! -----------------------------------
    integer(ip), intent(in)             :: pgaus
    integer(ip), intent(in)             :: pmate
    integer(ip), intent(in)             :: flagt
    real(rp),    intent(in)             :: gpcau(3,3,pgaus)
    real(rp),    intent(in)             :: gpdet(pgaus)
    real(rp),    intent(in)             :: gpbasis(3,3,pgaus)
    real(rp),    intent(out)            :: gpstr(3,3,pgaus)
    real(rp),    intent(out)            :: gpdds(3,3,3,3,pgaus)
    ! -----------------------------------
    integer(ip)                         :: igaus, idime, jdime, kdime, ldime
    real(rp)                            :: Kct, a, b, af, bf, as, bs, afs, bfs
    real(rp)                            :: i1, i1_iso, i4f, i4f_iso, i4s, i4s_iso !, i4n, i4n_iso
    real(rp)                            :: i4fp, i4sp, i8fs, i8fs_iso
    real(rp)                            :: term_i1, term_i4f, term_i4s, term_i8fs, term_vol
    real(rp)                            :: J, J23, J43
    real(rp), dimension(3,3)            :: fxf_0, fxf_dev, sxs_0, sxs_dev, fxs_0, fxs_dev, IxI_dev
    real(rp), dimension(3,3)            :: spk_vol, spk_iso, spk_4f, spk_4s, spk_8fs
    real(rp), dimension(3,3,3,3)        :: ixCi, Cixi, CixCi, I_Ci, i1xi1, fxfxCi
    real(rp), dimension(3,3,3,3)        :: Cixfxf, i4fxi4f, sxsxCi, Cixsxs, i4sxi4s
    real(rp), dimension(3,3,3,3)        :: fsxfsxCi, Cixfsxfs, i8fsxi8fs
    real(rp), dimension(3,3,3,3)        :: dspk_4f, dspk_4s, dspk_8fs, dspk_iso, dspk_vol
    real(rp), dimension(3,3)            :: gpcai, aux_M3x3
    real(rp), dimension(3)              :: fbr0, sht0 !, nrm0
    ! -----------------------------------
    real(rp), parameter, dimension(3,3) :: tkron = reshape( (/             &
    &                                           1.0_rp, 0.0_rp, 0.0_rp   , &
    &                                           0.0_rp, 1.0_rp, 0.0_rp   , &
    &                                           0.0_rp, 0.0_rp, 1.0_rp /), &
    &                                           (/ 3 , 3 /) )
    ! -----------------------------------
    ! Retrieve the parameters
    Kct   = parco_sld(2,pmate)
    a     = parco_sld(3,pmate)
    b     = parco_sld(4,pmate)
    af    = parco_sld(5,pmate)
    bf    = parco_sld(6,pmate)
    as    = parco_sld(7,pmate)
    bs    = parco_sld(8,pmate)
    afs   = parco_sld(9,pmate)
    bfs   = parco_sld(10,pmate)

    ! -----------------------------------
    ! Passive part
    gpstr = 0.0_rp
    gpdds = 0.0_rp
    loop_gp_passive: &
    do igaus = 1, pgaus

        ! Recober reference longitudinal (fiber) and transversal (sheet, norma) direction
        fbr0 = gpbasis(1:3,1,igaus)
        sht0 = gpbasis(1:3,2,igaus)
        !nrm0 = gpbasis(1:3,3,igaus)

        ! Compute the invariants
        i1   =  TRACE(gpcau(:,:,igaus),3_ip,3_ip)
        i4f  =  maths_MULT_VxMxV(fbr0(:),3_ip,gpcau(:,:,igaus),fbr0(:),3_ip)
        i4s  =  maths_MULT_VxMxV(sht0(:),3_ip,gpcau(:,:,igaus),sht0(:),3_ip)
        !i4n  =  maths_MULT_VxMxV(nrm0(:),3_ip,gpcau(:,:,igaus),nrm0(:),3_ip)
        i8fs = (maths_MULT_VxMxV(fbr0(:),3_ip,gpcau(:,:,igaus),sht0(:),3_ip) + &
        &       maths_MULT_VxMxV(sht0(:),3_ip,gpcau(:,:,igaus),fbr0(:),3_ip))*0.5_rp

        ! Different needed variables related to the Deformation Gradient
        J   = gpdet(igaus)
        J23 = J**(-(2.0_rp/3.0_rp))
        J43 = J**(-(4.0_rp/3.0_rp))

        ! Different expressions of the scalar portion of the volumetric stress tensor
        term_vol = (Kct/4.0_rp)*((J**2) - 1.0_rp)

        ! Compute the isochoric invariants
        i1_iso   = J23*i1
        i4f_iso  = J23*i4f
        i4s_iso  = J23*i4s
        !i4n_iso  = J23*i4n
        i8fs_iso = J23*i8fs

        ! Tension/compression asymmetry
        i4fp = 0.5_rp*((i4f_iso - 1.0_rp) + abs(i4f_iso - 1.0_rp))
        i4sp = 0.5_rp*((i4s_iso - 1.0_rp) + abs(i4s_iso - 1.0_rp))

        ! exponential isotropic and isochoric term
        term_i1 = exp(b*(i1_iso - 3.0_rp))

        ! These checks re necessary, because round-off errors from i4f/s 
        ! can create spurious small forces in the passive 
        i4f     = ROUNDOFF(i4f,1.0_rp)
        i4s     = ROUNDOFF(i4s,1.0_rp)
        i4f_iso = ROUNDOFF(i4f_iso,1.0_rp)
        i4s_iso = ROUNDOFF(i4s_iso,1.0_rp)

        ! Compute exponential anisotropic and isochoric terms
        term_i4f  = exp(bf*(i4f_iso - 1.0_rp)**2)
        term_i4s  = exp(bs*(i4s_iso - 1.0_rp)**2)
        term_i8fs = exp(bfs*(i8fs_iso**2))

        ! Invert Right Cauchy-Green strain tensor C (gpcai)
        gpcai(:,:) = INV3x3(gpcau(:,:,igaus))

        ! f x f , s x s, f x s, s x f (outer products [3x3])
        fxf_0(:,:) =  OUTPROD_VxV(fbr0(:),3_ip,fbr0(:),3_ip)
        sxs_0(:,:) =  OUTPROD_VxV(sht0(:),3_ip,sht0(:),3_ip)
        fxs_0(:,:) = (OUTPROD_VxV(fbr0(:),3_ip,sht0(:),3_ip) + &
        &             OUTPROD_VxV(sht0(:),3_ip,fbr0(:),3_ip))*0.5_rp

        ! Deviatoric part: Dev(Aij) = Aij - (1/3)Akk tkron_ij
        aux_M3x3(:,:)  = gpcai(:,:)/3.0_rp
        IxI_dev(:,:) = tkron(:,:) - i1*aux_M3x3(:,:)
        fxf_dev(:,:) = fxf_0(:,:) - i4f*aux_M3x3(:,:)
        sxs_dev(:,:) = sxs_0(:,:) - i4s*aux_M3x3(:,:)
        fxs_dev(:,:) = fxs_0(:,:) - i8fs*aux_M3x3(:,:)

        ! Second Piola-Kirchhoff stress tensor
        spk_vol(:,:) = 2.0_rp*term_vol*gpcai(:,:)
        spk_iso(:,:) = J23*a*term_i1*IxI_dev(:,:)
        spk_4f(:,:)  = 2.0_rp*J23*af*i4fp*term_i4f*fxf_dev(:,:)
        spk_4s(:,:)  = 2.0_rp*J23*as*i4sp*term_i4s*sxs_dev(:,:)
        spk_8fs(:,:) = 2.0_rp*J23*afs*i8fs_iso*term_i8fs*fxs_dev(:,:)

        gpstr(:,:,igaus) = spk_vol(:,:) + spk_iso(:,:) + spk_4f(:,:) + &
        &                  spk_4s(:,:)  + spk_8fs(:,:)

        ! Tangent moduli
        if( flagt == 1_ip )then

            ! Compute some outper products
            IxCi(:,:,:,:)      = OUTPROD_MxM(tkron(:,:),3_ip,3_ip,gpcai(:,:),3_ip,3_ip)
            Cixi(:,:,:,:)      = OUTPROD_MxM(gpcai(:,:),3_ip,3_ip,tkron(:,:),3_ip,3_ip)
            CixCi(:,:,:,:)     = OUTPROD_MxM(gpcai(:,:),3_ip,3_ip,gpcai(:,:),3_ip,3_ip)
            i1xi1(:,:,:,:)     = OUTPROD_MxM(IxI_dev(:,:),3_ip,3_ip,IxI_dev(:,:),3_ip,3_ip)
            fxfxCi(:,:,:,:)    = OUTPROD_MxM(fxf_0(:,:),3_ip,3_ip,gpcai(:,:),3_ip,3_ip)
            Cixfxf(:,:,:,:)    = OUTPROD_MxM(gpcai(:,:),3_ip,3_ip,fxf_0(:,:),3_ip,3_ip)
            i4fxi4f(:,:,:,:)   = OUTPROD_MxM(fxf_dev(:,:),3_ip,3_ip,fxf_dev(:,:),3_ip,3_ip)
            sxsxCi(:,:,:,:)    = OUTPROD_MxM(sxs_0(:,:),3_ip,3_ip,gpcai(:,:),3_ip,3_ip)
            Cixsxs(:,:,:,:)    = OUTPROD_MxM(gpcai(:,:),3_ip,3_ip,sxs_0(:,:),3_ip,3_ip)
            i4sxi4s(:,:,:,:)   = OUTPROD_MxM(sxs_dev(:,:),3_ip,3_ip,sxs_dev(:,:),3_ip,3_ip)
            fsxfsxCi(:,:,:,:)  = OUTPROD_MxM(fxs_0(:,:),3_ip,3_ip,gpcai(:,:),3_ip,3_ip)
            Cixfsxfs(:,:,:,:)  = OUTPROD_MxM(gpcai(:,:),3_ip,3_ip,fxs_0(:,:),3_ip,3_ip)
            i8fsxi8fs(:,:,:,:) = OUTPROD_MxM(fxs_dev(:,:),3_ip,3_ip,fxs_dev(:,:),3_ip,3_ip)

            do idime = 1, 3_ip
                do jdime = 1, 3_ip
                    do kdime = 1, 3_ip
                        do ldime = 1, 3_ip
                            i_Ci(idime,jdime,kdime,ldime) = (gpcai(idime,kdime)*gpcai(jdime,ldime)+ &
                            gpcai(idime,ldime)*gpcai(jdime,kdime))                 
                        end do
                    end do
                end do
            end do

            ! Compute tangent contributions
            ! - isotropic
            dspk_iso = 2.0_rp*J23*a*term_i1*(-IxCi/3.0_rp - Cixi/3.0_rp + i1*CixCi/9.0_rp + i1*i_Ci/6.0_rp) &
            &   + 2.0_rp*J43*a*b*term_i1*i1xi1

            ! - anisotropic
            dspk_4f = 4.0_rp*J23*af*i4fp*term_i4f*(-fxfxCi/3.0_rp - Cixfxf/3.0_rp + i4f*CixCi/9.0_rp + i4f*i_Ci/6.0_rp) + &
            &   4.0_rp*J43*af*term_i4f*i4fxi4f + 8.0_rp*J43*af*bf*(i4fp**2)*term_i4f*i4fxi4f

            dspk_4s = 4.0_rp*J23*as*i4sp*term_i4s*(-sxsxCi/3.0_rp - Cixsxs/3.0_rp + i4s*CixCi/9.0_rp + i4s*i_Ci/6.0_rp) + &
            &   4.0_rp*J43*as*term_i4s*i4sxi4s + 8.0_rp*J43*as*bs*(i4sp**2)*term_i4s*i4sxi4s

            dspk_8fs = 4.0_rp*J23*afs*i8fs_iso*term_i8fs*(-fsxfsxCi/3.0_rp - Cixfsxfs/3.0_rp + i8fs*CixCi/9.0_rp + i8fs*i_Ci/6.0_rp) + &
            &   4.0_rp*J43*afs*term_i8fs*i8fsxi8fs + 8.0_rp*J43*afs*bfs*(i8fs_iso**2)*term_i8fs*i8fsxi8fs

            ! - volumetric
            dspk_vol = Kct*((J**2)*CixCi - 0.5_rp*((J**2) - 1.0_rp)*i_Ci)

            ! Add contributions
            gpdds(:,:,:,:,igaus) = dspk_iso(:,:,:,:) + dspk_4f(:,:,:,:) + dspk_4s(:,:,:,:) +  &
            &   dspk_8fs(:,:,:,:) + dspk_vol(:,:,:,:) 
             
        endif

    end do &
    loop_gp_passive

    ! -----------------------------------
    contains 

        pure function ROUNDOFF(val,tgt) result(res)
            ! -------------------------------
            implicit none
            ! -------------------------------
            real(rp), intent(in)           :: val
            real(rp), intent(in)           :: tgt
            real(rp)                       :: res
            real(rp), parameter            :: thr = 1.0e-12_rp
            ! -------------------------------
            if( abs((val - tgt)) < thr )then
                res = tgt
            else
                res = val
            endif
            ! -------------------------------
        end function ROUNDOFF

        pure function TRACE(M,mi,mj) result(r)
            ! -------------------------------
            implicit none
            ! -------------------------------
            integer(ip), intent(in)                :: mi, mj
            real(rp), dimension(mi,mj), intent(in) :: M
            real(rp)                               :: r
            integer(ip)                            :: i
            ! -------------------------------
            r = 0.0_rp
            do i = 1, mi
                r = r+ M(i,i)
            enddo
            ! -------------------------------
        end function TRACE

        pure function INV3x3(M) result(Mi)
            ! -------------------------------
            implicit none
            ! -------------------------------
            real(rp), dimension(3,3), intent(in) :: M
            real(rp), dimension(3,3)             :: Mi
            real(rp)                             :: det
            ! -------------------------------
            ! Determinant
            det = M(1,1)*M(2,2)*M(3,3) + M(1,3)*M(2,1)*M(3,2) + &
                M(3,1)*M(1,2)*M(2,3) - M(3,1)*M(2,2)*M(1,3) - &
                M(3,3)*M(1,2)*M(2,1) - M(1,1)*M(2,3)*M(3,2)
            ! Invert matrix
            Mi(1,1) =  (M(2,2)*M(3,3) - M(3,2)*M(2,3))/det
            Mi(1,2) = -(M(1,2)*M(3,3) - M(1,3)*M(3,2))/det
            Mi(1,3) =  (M(1,2)*M(2,3) - M(2,2)*M(1,3))/det
            Mi(2,1) = -(M(2,1)*M(3,3) - M(3,1)*M(2,3))/det
            Mi(2,2) =  (M(1,1)*M(3,3) - M(1,3)*M(3,1))/det
            Mi(2,3) = -(M(1,1)*M(2,3) - M(2,1)*M(1,3))/det
            Mi(3,1) =  (M(2,1)*M(3,2) - M(3,1)*M(2,2))/det
            Mi(3,2) = -(M(1,1)*M(3,2) - M(3,1)*M(1,2))/det
            Mi(3,3) =  (M(1,1)*M(2,2) - M(1,2)*M(2,1))/det
            ! -------------------------------
        end function INV3x3

        pure function MULT_MxV(M,mi,mj,v,vj) result(r)
            ! -------------------------------
            implicit none
            ! -------------------------------
            integer(ip), intent(in)                :: mi, mj, vj
            real(rp), dimension(mi,mj), intent(in) :: M
            real(rp), dimension(vj),    intent(in) :: v
            real(rp), dimension(mi)                :: r
            integer(ip)                            :: i, j
            ! -------------------------------
            r(:) = 0.0_rp
            do j = 1, mi
                do i = 1, mj
                    r(i) =  r(i) + M(i,j)*v(j)
                enddo
            enddo
            ! -------------------------------
        end function MULT_MxV

        pure function OUTPROD_VxV(a,ai,b,bi) result(R)
            ! EQUIV to dot_product(a,b)
            !   where A and B are matrices
            ! -------------------------------
            implicit none
            ! -------------------------------
            integer(ip), intent(in)             :: ai, bi
            real(rp), dimension(ai), intent(in) :: a
            real(rp), dimension(bi), intent(in) :: b
            real(rp), dimension(ai,bi)          :: R
            integer(ip)                         :: i, j
            ! -------------------------------
            do i = 1, ai
                do j = 1, bi
                    R(i,j) = a(i)*b(j)
                enddo
            enddo
            ! -------------------------------
        end function OUTPROD_VxV

        pure function OUTPROD_MxM(A,ai,aj,B,bi,bj) result(R)
            ! -------------------------------
            implicit none
            ! -------------------------------
            integer(ip), intent(in)                :: ai, aj, bi, bj
            real(rp), dimension(ai,aj), intent(in) :: a
            real(rp), dimension(bi,bj), intent(in) :: b
            real(rp), dimension(ai,aj,bi,bj)       :: R
            integer(ip)                            :: i, j, k, l
            ! -------------------------------
            do l = 1, bj
                do k = 1, bi
                    do j = 1, aj
                        do i = 1, ai
                            R(i,j,k,l) = A(i,j)*B(k,l)
                        enddo
                    enddo
                enddo
            enddo
            ! -------------------------------
        end function OUTPROD_MxM

    ! -----------------------------------
end subroutine sld_stress_model_134

