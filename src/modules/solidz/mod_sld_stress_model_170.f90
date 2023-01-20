!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-------------------------------------------------------------------------------
!> @authors Adria Quintanas-Corominas : adria.quintanas@bsc.es
!> @authors Gerard Guillaumet         : gerard.guillamet@bsc.es
!> @date    October, 2020
!------------------------------------------------------------------------------
module mod_sld_stress_model_170
    ! ----------------------------------------
    use def_kintyp_basic,               only :  ip, rp
    use def_domain,                     only :  ndime
    use def_master,                     only :  cutim
    ! ----------------------------------------
    implicit none
    ! ----------------------------------------

    public  :: sm170_get_sound_velocity
    public  :: sm170_get_stress
    private :: Identity
    private :: GreenLagrange
    private :: trace

    contains

    subroutine sm170_get_stress( &
        itask, pgaus, props, gpgdi, gptmp, gpstr, gpdds)
        !-------------------------------------
        implicit none
        !-------------------------------------
        integer(ip),   intent(in)           :: itask
        integer(ip),   intent(in)           :: pgaus
        real(rp),      intent(in)           :: props(:)
        real(rp),      intent(in)           :: gpgdi(:,:,:)
        real(rp),      intent(in)           :: gptmp(:,:)
        real(rp),      intent(out)          :: gpstr(:,:,:)
        real(rp),      intent(out)          :: gpdds(:,:,:,:,:)
        integer(ip)                         :: ig, ii, jj, kk, ll
        real(rp)                            :: Em, nu, lmb, mu, kappa, alpha, tmp1, tmp0
        real(rp)                            :: I(ndime,ndime)
        real(rp)                            :: F(ndime,ndime)
        real(rp)                            :: E(ndime,ndime)
        real(rp)                            :: E_tr
        real(rp)                            :: tmp
        real(rp)                            :: S(ndime,ndime)
        !-------------------------------------
        Em = props(1)
        nu = props(2)
        alpha = props(3)

        lmb = Em * nu / ( (1.0_rp + nu) * (1.0_rp - 2.0_rp * nu) )
        mu = Em / (2.0_rp + 2.0_rp * nu)
        kappa = Em / (3.0_rp - 6.0_rp * nu)

        I(:,:) = Identity()

        ! Second Piola-Kirchoff stress tensor
        do ig = 1, pgaus
            F = gpgdi(:,:,ig)
            tmp0 = gptmp(4,ig) ! Initial simulation value
            tmp1 = gptmp(1,ig) ! Current iteration value
            E = GreenLagrange(F)
            E_tr = trace(E)
            S = lmb * E_tr * I + 2.0_rp * mu * E - 3.0_rp * alpha * ( tmp1 - tmp0 ) * kappa * I
            gpstr(:,:,ig) = S
        enddo

        ! Second elasticity tensor (dSdE)
        gpdds(:,:,:,:,:) = 0.0_rp
        if( itask == 1_ip )then

           do ig = 1, pgaus
               do ll = 1, ndime
                   do kk = 1, ndime
                       do jj = 1, ndime
                           do ii = 1, ndime
                               gpdds(ii,jj,kk,ll,ig) = lmb * I(ii,jj) * I(kk,ll) + mu * ( I(ii,kk) * I(jj,ll) + I(ii,ll) * I(jj,kk) )
                           enddo
                       enddo
                   enddo
               enddo
           enddo

        endif
        !-------------------------------------
    end subroutine sm170_get_stress

    subroutine sm170_get_sound_velocity( & 
        props, rho, velsnd)
        !-------------------------------------
        implicit none
        !-------------------------------------
        real(rp),      intent(in)           :: props(:)
        real(rp),      intent(in)           :: rho
        real(rp),      intent(out)          :: velsnd
        real(rp)                            :: Em, nu, lmb, mu, M
        !-------------------------------------
        Em = props(1)
        nu = props(2)

        lmb = Em * nu / ( (1.0_rp + nu) * (1.0_rp - 2.0_rp * nu) )
        mu = Em / (2.0_rp + 2.0_rp * nu)
        M = lmb + 2.0_rp * mu

        velsnd = sqrt( M / rho )
        !-------------------------------------
    end subroutine sm170_get_sound_velocity

    function Identity() result(Id)
        use def_domain, only : ndime
        implicit none
        real(rp)            :: Id(ndime,ndime)
#ifndef NDIMEPAR
        if( ndime == 2_ip )then
            Id = reshape( (/ 1.0_rp, 0.0_rp, &
                             0.0_rp, 1.0_rp  /) , (/2,2/))
        else

            Id = reshape( (/ 1.0_rp, 0.0_rp, 0.0_rp, &
                             0.0_rp, 1.0_rp, 0.0_rp, & 
                             0.0_rp, 0.0_rp, 1.0_rp /) , (/3,3/))
        endif
#else
        Id = reshape( (/ 1.0_rp, 0.0_rp, 0.0_rp, &
                         0.0_rp, 1.0_rp, 0.0_rp, &
                         0.0_rp, 0.0_rp, 1.0_rp /) , (/3,3/))
#endif
    end function Identity

    function GreenLagrange(F) result(E)
        use def_domain, only : ndime
        implicit none
        real(rp), intent(in) :: F(:,:)
        real(rp)             :: E(ndime,ndime)
        integer(ip)          :: ii, jj, kk
        E = 0.0_rp
        do kk = 1, ndime
            do jj = 1, ndime
                do ii = 1, ndime
                    E(ii,jj) = E(ii,jj) + 0.5_rp * F(kk,ii) * F(kk,jj)
                enddo
            enddo
        enddo
        do ii = 1, ndime
            E(ii,ii) = E(ii,ii) - 0.5_rp
        enddo
    end function GreenLagrange

    function trace(T) result(r)
        use def_domain, only : ndime
        implicit none
        real(rp), intent(in) :: T(:,:)
        real(rp)             :: r
        integer(ip)          :: ii
        r = 0.0_rp
        do ii = 1, ndime
            r = r + T(ii,ii)
        enddo
    end function trace

end module mod_sld_stress_model_170
