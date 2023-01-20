!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_sld_stress_model_400

   use def_kintyp, only : ip, rp, lg
   use def_domain, only : ndime
   use def_solidz, only : enden_sld, water_sld
   use def_master, only : donna_gp
   
   implicit none
   

   public ::             &
   sld_stress_model_400

   private ::           &
   sm400_hyper_PK2,     &
   sm400_hyper_dpdds,   &
   sm400_poro_PK2,      &
   sm400_poro_dpdds,    &
   sm400_donnan_PK2,    &
   sm400_donnan2_PK2,   &
   sm400_donnan_dpdds,  &
   sm400_donnan2_dpdds, &
   sm400_fiber_PK2,     &
   sm400_fiber_dpdds
   !!sm400_fib1_dbg,      &
   !!sm400_fib2_dbg

   contains

   subroutine sld_stress_model_400(ielem, pgaus, pnode, pmate, gpgdi, gpstr, gpdet, flagt, gpdds, gpene, gpcau, gpsha)
      ! Is intented to be used with STATIC. DYNAMIC is also implemented, but should only be used if 
      ! no fibers exist, as they are not considered in the calculation of the critical timestep.
      
      ! When using fibers, 2 unit vectors regarding the initial orientation of the fibers must be passed.
      ! They are passed in the same field.
      ! Reference tests: o Pure hyperelastic + porous -> sm400_hyperelastic
      !                  o Hyperelastic + porous + swelling -> sm400_hyperelastic_swelling
      !                  o Hyperelastic + fibers -> sm400_hyperelastic_fibers

      use def_solidz, only : dtinv_sld
      use def_solidz, only : kfl_donna_sld, kfl_fiber_sld
      use def_solidz, only : axis1_sld, axis2_sld!, fibep_sld
      use def_master, only : dtime, cutim
      use mod_maths,  only : maths_norm2
      !!use def_domain, only : nfiel
      !!use mod_poroelasticity, only : kfl_poroelasticity, poro_sld_get_pressure_at_GP

      implicit none

      integer(ip), intent(in) :: ielem, pgaus, pmate, flagt, pnode
      real(rp),    intent(in) :: gpgdi(ndime, ndime, pgaus), gpcau(ndime, ndime, pgaus), gpdet(pgaus), &
                                 gpsha(pnode, pgaus)
      real(rp), intent(inout) :: gpstr(ndime, ndime, pgaus, 2)
      real(rp), intent(out)   :: gpene(pgaus), gpdds(ndime, ndime, ndime, ndime, pgaus)
      
      integer(ip)  :: igaus, ii
      real(rp)     :: gpprs(pgaus), gpdon(pgaus)
      real(rp)     :: gplng1(3, pgaus), gplng2(3, pgaus)
      real(rp)     :: gpwat(pgaus)
      logical(lg)  :: kfl_fib_internal

      ! Init 
      kfl_fib_internal = .False.

      gpstr(:, :, :, 1) = 0.0_rp
      gpene = 0.0_rp
      gpdds = 0.0_rp
      gpprs = 0.0_rp
      gpdon = 0.0_rp

      !!if (kfl_poroelasticity == 1_ip) then
      !!   call poro_sld_get_pressure_at_GP(ielem, pnode, pgaus, gpsha, gpprs)
      !!end if
      if (kfl_fiber_sld == 4_ip) then
         do igaus = 1, pgaus
            gplng1(:, igaus) = axis1_sld(:, ielem)
            gplng2(:, igaus) = axis2_sld(:, ielem)
         end do
         ! This internal flag is used because there is currently no other way to avoid
         ! looping for non-fibrous materials.
         if (maths_norm2(ndime, gplng1(:, 1)) .ge. 1e-5) kfl_fib_internal = .True.
      end if

      do igaus = 1, pgaus
         
         ! Energy only included for the hyperelastic part
         call sm400_hyper_PK2(pmate, gpgdi(:, :, igaus), gpstr(:, :, igaus, 1), gpdet(igaus),&
                              gpene(igaus), gpcau(:, :, igaus), gpwat(igaus))
         enden_sld(ielem) % a(igaus) = gpene(igaus)
         water_sld(ielem) % a(igaus) = gpwat(igaus)
         if (kfl_fib_internal) then 
            call sm400_fiber_PK2(pmate, gpcau(:, :, igaus), gplng1(:, igaus), gpstr(:, :, igaus, 1))
            call sm400_fiber_PK2(pmate, gpcau(:, :, igaus), gplng2(:, igaus), gpstr(:, :, igaus, 1))
         end if
         !!if (kfl_poroelasticity == 1_ip) then
         !!   call sm400_poro_PK2(gpgdi(:, :, igaus), gpdet(igaus), gpstr(:, :, igaus, 1), gpprs(igaus))
         !!end if
         if (kfl_donna_sld == 1_ip) then
            call sm400_donnan2_PK2(pmate, gpgdi(:, :, igaus), gpdet(igaus), gpstr(:, :, igaus, 1), gpdon(igaus))
            !!! For postprocess reasons and to send to def_master
            donna_gp(ielem) % a(igaus) = gpdon(igaus)
         end if

         ! If implicit
         if (flagt == 1_ip) then
            call sm400_hyper_dpdds(pmate, gpcau(:, :, igaus), gpdet(igaus), gpdds(:, :, :, :, igaus))
            if (kfl_fib_internal) then 
               call sm400_fiber_dpdds(pmate, gpcau(:, :, igaus), gplng1(:, igaus), gpdds(:, :, :, :, igaus))
               call sm400_fiber_dpdds(pmate, gpcau(:, :, igaus), gplng2(:, igaus), gpdds(:, :, :, :, igaus))
            end if         
 !!  if (kfl_poroelasticity == 1_ip) then
          !!     call sm400_poro_dpdds(gpcau(:, :, igaus), gpdet(igaus), gpprs(igaus), gpdds(:, :, :, :, igaus))
          !!  end if
            if (kfl_donna_sld == 1_ip) then
               call sm400_donnan2_dpdds(pmate, gpgdi(:, :, igaus), gpdet(igaus), gpdds(:, :, :, :, igaus))
            end if
         end if
      end do
   
   end subroutine sld_stress_model_400

   subroutine sm400_hyper_PK2(pmate, F, PK2, J, W, C, nf)

      use def_solidz, only       :  parco_sld

      implicit none
  
      integer(ip), intent(in)    :: pmate
      real(rp),    intent(in)    :: C(ndime,ndime)
      real(rp),    intent(in)    :: F(ndime,ndime), J
      real(rp),    intent(out)   :: W
      real(rp),    intent(out)   :: nf
      real(rp),    intent(inout) :: PK2(ndime,ndime)
        
      real(rp)                   :: tkron(ndime,ndime)
      real(rp)                   :: B(ndime, ndime), sigma(ndime, ndime)    ! left cauchy green deformation
                                                                            ! tensor and Cauchy stress tensor
      real(rp)                   :: invF(ndime, ndime), invC(ndime, ndime)
      integer(ip)                :: ii, jj, kk, ll
      real(rp)                   :: K0, tracrc, bidon
      real(rp)                   :: a, ns0, Gm ! a(J) = -1 + 3*(J+ns0)/(-J+ns0) + 3*J*log(J)*ns0/(-J + ns0)^2

      ! Modified Neo-Hookean law
      Gm      = parco_sld(1, pmate)  ! Shear modulus
      ns0     = parco_sld(2, pmate)  ! Initial solid fraction 

      sigma = 0.0_rp
      tkron = 0.0_rp
      do ii = 1, ndime
         tkron(ii, ii) = 1.0_rp
      end do

      ! Form left cauchy green deformation tensor
      B = 0.0_rp
      do ii = 1, ndime; do jj = 1, ndime; do kk = 1, ndime
         B(ii, jj) = B(ii, jj) + F(ii, kk) * F(jj, kk)
      end do; end do; end do

      ! Compute the J dependent term, a
      a = -1.0_rp + 3.0_rp * (J + ns0) / (-J + ns0) + 3.0_rp * J * log(J) * ns0 / (- J + ns0) ** 2

      ! Compute the inverse of the deformation gradient tensor
      call invmtx(F, invF, bidon, ndime)

      do ii =1, ndime
         do jj = 1, ndime   
            do kk = 1, ndime
               do ll = 1, ndime
                  ! Compute the cauchy stress tensor
                  sigma(kk, ll) = - 1.0_rp / 6.0_rp * log(J) / J * Gm * tkron(kk, ll) * a &
                  + Gm / J * (B(kk, ll) - (J ** (2.0_rp/3.0_rp)) * tkron(kk, ll))
                  ! Compute the second piola kirchoff stress tensor
                  PK2(ii, jj) = PK2(ii, jj) + J * invF(ii, kk) * sigma(kk, ll) * invF(jj, ll)
               end do
            end do
         end do
      end do

      ! Compute the trace of the right cauchy-green deformation tensor
      tracrc = trace(C, ndime)

      ! Compute the stored energy for each gauss point
      W = 1.0_rp / 3.0_rp * Gm * (1.0_rp + 0.5_rp * ns0 / J) / (1.0_rp - ns0 / J) * log(J) ** 2 &
      + 0.5_rp * Gm * (tracrc - 3 * J ** (2.0_rp / 3.0_rp))
      ! Compute water content
      nf = (J - ns0) / J

   end subroutine sm400_hyper_PK2 

   subroutine sm400_hyper_dpdds(pmate, C, J, dds)

      use def_solidz, only       :  parco_sld

      implicit none
  
      integer(ip), intent(in) :: pmate
      real(rp),    intent(in) :: C(ndime,ndime)
      real(rp),    intent(in) :: J
      real(rp), intent(inout) :: dds(ndime, ndime, ndime, ndime)

      real(rp)    :: detC, bidon
      real(rp)    :: invC(ndime, ndime)
      integer(ip) :: ii, jj, kk, ll
      real(rp)    :: Km, Y1, Y2, ds1, ds2, ds3
      real(rp)    :: K0, Gm, ns0
      
      ! Modified Neo-Hookean law
      Gm      = parco_sld(1, pmate)  ! Shear modulus
      ns0     = parco_sld(2, pmate)  ! Initial solid fraction 

      detC = DET(C, ndime)
      call invmtx(C, invC, bidon, ndime)

      do ii = 1,ndime; do jj = 1,ndime; do kk = 1,ndime; do ll = 1,ndime
         Km  =  2.0_rp / 3.0_rp * Gm * (J + 0.5_rp * ns0) / (J - ns0) 
         Y1  = -0.5_rp * (invC(jj, kk) * invC(ll, ii) + invC(jj, ll) * invC(kk, ii))
         Y2  = -0.5_rp * (invC(ii, kk) * invC(ll, jj) + invC(ii, ll) * invC(kk, jj))

         ds1 = -0.125_rp * Gm * ns0 * (&
              (-ns0 - J) / (J - ns0) ** 3_ip * J * 0.5_rp * log(detC) ** 2_ip * invC(ii, jj) * invC(kk, ll) + &
              J / (J - ns0) ** 2_ip * 2.0_rp * log(detC) * invC(ii, jj) * invC(ll, kk) +&
              J / (J - ns0) ** 2_ip * log(detC) ** 2_ip * Y2) 
         ds2 = -Gm * detC ** (1.0_rp / 3.0_rp) * (1.0_rp / 3.0_rp * invC(jj, ii) * invC(ll, kk) + Y1) 
         ds3 = 0.5_rp * (- Gm * ns0 / (J - ns0) ** 2_ip * J * 0.5_rp * log(detC) * invC(kk, ll) * invC(jj, ii) + &
              Km * invC(ll, kk) * invC(jj, ii) + &
              Km * log(detC) * Y1)

         dds(ii, jj, kk, ll) = dds(ii, jj, kk, ll) + 2.0_rp * (ds1 + ds2 + ds3)
      end do; end do; end do; end do 

   end subroutine sm400_hyper_dpdds

   subroutine sm400_fiber_PK2(pmate, C, v0, PK2)

      use def_solidz, only : parco_sld
      
      implicit none

      integer(ip), intent(in) :: pmate
      real(rp), intent(in)    :: C(ndime, ndime), v0(ndime)

      real(rp), intent(inout) :: PK2(ndime, ndime)

      integer(ip) :: ii, jj
      real(rp)    :: k1, k2, I4

      k1 = parco_sld(8, pmate)
      k2 = parco_sld(9, pmate)

      I4 = invariant4(C, v0)

      do ii = 1, ndime
         do jj = 1, ndime
            PK2(ii, jj) = PK2(ii, jj) + 2.0_rp * k1 * Aaniso(k2, I4) * (I4 - 1.0_rp) * v0(ii) * v0(jj)
         end do
      end do

   end subroutine sm400_fiber_PK2
      
   subroutine sm400_fiber_dpdds(pmate, C, v0, dds)

      use def_solidz, only : parco_sld
      
      implicit none

      integer(ip), intent(in) :: pmate
      real(rp), intent(in)    :: C(ndime, ndime), v0(ndime)

      real(rp), intent(inout) :: dds(ndime, ndime, ndime, ndime)
      
      integer(ip) :: ii, jj, kk, ll
      real(rp)    :: k1, k2, I4, Aa
      real(rp), dimension(ndime, ndime) :: dAdC

      k1 = parco_sld(8, pmate)
      k2 = parco_sld(9, pmate)

      I4 = invariant4(C, v0)
      Aa = Aaniso(k2, I4)
      dAdC = dAaniso(k2, Aa, C, v0)

      do ii = 1, ndime; do jj = 1, ndime; do kk = 1, ndime; do ll = 1, ndime
        dds(ii, jj, kk, ll) = dds(ii, jj, kk, ll) &
                                + 4.0_rp * k1 * Aa * (2 * k2 * (I4 - 1) ** 2_ip + 1) &
                                * v0(ii) * v0(jj) * v0(kk) * v0(ll)
      enddo; enddo; enddo; enddo

   end subroutine sm400_fiber_dpdds

   subroutine sm400_poro_PK2(F, J, PK2, prs)
   
      implicit none
      
      real(rp), intent(in)    :: F(ndime, ndime)
      real(rp), intent(in)    :: J, prs
      real(rp), intent(inout) :: PK2(ndime, ndime)

      integer(ip) :: ii, jj, kk
      real(rp)    :: presu(ndime, ndime), invF(ndime, ndime)
      real(rp)    :: bidon

      call invmtx(F, invF, bidon, ndime)

      presu = 0.0_rp

      do ii = 1, ndime
         do jj = 1, ndime
            do kk = 1, ndime
               presu(ii, jj) = presu(ii, jj) - J * invF(ii, kk) &
                              * prs * invF(jj, kk)
            end do
         end do
      end do
      
      PK2 = presu + PK2

   end subroutine sm400_poro_PK2

   subroutine sm400_poro_dpdds(C, J, prs, dds)
     
      implicit none

      real(rp), intent(in) :: C(ndime, ndime)
      real(rp), intent(in) :: J, prs

      real(rp), intent(inout) :: dds(ndime, ndime, ndime, ndime)

      integer(ip) :: ii, jj, kk, ll
      real(rp)    :: invC(ndime, ndime)
      real(rp)    :: bidon

      call invmtx(C, invC, bidon, ndime)
      do ii = 1, ndime; do jj = 1, ndime; do kk = 1, ndime; do ll = 1, ndime
         dds(ii, jj, kk, ll) = dds(ii, jj, kk, ll) + &
         J * prs  * (invC(kk, ll) * invC(ii, jj) - &
         invC(ii, kk) * invC(ll, jj) - invC(ii, ll) * invC(kk, jj))
      end do; end do; end do; end do
   
   end subroutine sm400_poro_dpdds
   
   subroutine sm400_donnan_PK2(pmate, F, J, PK2, don)
      ! UNUSED MODEL
      use def_solidz, only  :  parco_sld
   
      implicit none
      
      integer(ip), intent(in) :: pmate
      real(rp), intent(in)    :: F(ndime, ndime)
      real(rp), intent(in)    :: J
      real(rp), intent(inout) :: PK2(ndime, ndime)
      real(rp), intent(out)   :: don

      integer(ip) :: ii, jj, kk
      real(rp)    :: presu(ndime, ndime), invF(ndime, ndime)
      real(rp)    :: bidon, donnan
      real(rp)    :: n_f, n_exf, phi_ci, rho_c, c_fexf, c_F0, ns0
      real(rp)    :: phi_int, phi_ext, gamma_int, gamma_ext, c_ext, Temp, c_F, rho_ctot, R_gas
      real(rp)    :: dummy

      call invmtx(F, invF, bidon, ndime)
      ns0       = parco_sld(2, pmate)
      ! Donnan Osmosis {Ref. Ruiz 2016, Simulating the sensitivity of cell nutritive environment \
      ! to composition changes within the intervertebral disc }
      phi_ci    = parco_sld(3, pmate)
      phi_int   = parco_sld(4, pmate)  ! Internal osmotic coefficient
      phi_ext   = parco_sld(5, pmate)  ! External osmotic coefficient
      gamma_int = parco_sld(6, pmate)  ! Internal activity coefficient
      gamma_ext = parco_sld(7, pmate) ! External activity coefficient
      c_ext     = parco_sld(8, pmate) ! External concentration of salt
      Temp      = parco_sld(9, pmate) ! Temperature of fluid
      c_F0      = parco_sld(10, pmate) ! Normal fixed charge density in mEq /millilitre of the total fluid
      rho_ctot  = parco_sld(11, pmate) ! Collagen content wrt total wet weight
      R_gas     = parco_sld(12, pmate) ! Ideal gas constant
      dummy     = parco_sld(13, pmate)

      n_f    = (J - ns0) / J
      n_exf  = n_f - phi_ci * rho_ctot
      c_F    = (1 - ns0) * c_F0 / (J - ns0) 
      c_fexf = n_f * c_F / n_exf
      don    =  R_gas * Temp * (phi_int * sqrt (c_fexf ** 2 + 4 * (gamma_ext / gamma_int) ** 2 * c_ext ** 2) &
                     - 2.0 * phi_ext * c_ext)
      don    = dummy

      presu = 0.0_rp

      do ii = 1, ndime
         do jj = 1, ndime
            do kk = 1, ndime
               presu(ii, jj) = presu(ii, jj) - J * invF(ii, kk) &
                              * don * invF(jj, kk)
            end do
         end do
      end do
      
      PK2 = presu + PK2
   end subroutine sm400_donnan_PK2

   subroutine sm400_donnan2_PK2(pmate, F, J, PK2, don)
      use def_solidz, only  :  parco_sld
      use def_master, only  : cutim
   
      implicit none
      
      integer(ip), intent(in) :: pmate
      real(rp), intent(in)    :: F(ndime, ndime)
      real(rp), intent(in)    :: J
      real(rp), intent(inout) :: PK2(ndime, ndime)
      real(rp), intent(out)   :: don

      integer(ip) :: ii, jj, kk
      real(rp)    :: presu(ndime, ndime), invF(ndime, ndime)
      real(rp)    :: bidon, donnan
      real(rp)    :: Rgas, Temp
      real(rp)    :: phi_0, c_F0, c_b, tramp, ramp_ratio
      real(rp)    :: c_F

      Rgas    = parco_sld(3, pmate)    ! Ideal gas constant
      Temp    = parco_sld(4, pmate)    ! Absolute Temperature [K]
      c_F0    = parco_sld(5, pmate)    ! Initial fixed charge density {positive or negative, I believe [mm^3/mol]}
      c_b     = parco_sld(6, pmate)    ! Osmolarity bath
      tramp   = parco_sld(7, pmate)    ! Time ramp for application of osmotic pressure
      phi_0   = 1.0_rp - parco_sld(2, pmate)    ! Initial water content

      call invmtx(F, invF, bidon, ndime)

      if (cutim < tramp) then
         ramp_ratio = cutim / tramp
      else
         ramp_ratio = 1.0_rp
      end if

      c_F = phi_0 / (J - 1 + phi_0) * c_F0
      don = ramp_ratio * Rgas * Temp * (sqrt(c_F**2 + c_b**2) - c_b)
      
      presu = 0.0_rp

      do ii = 1, ndime
         do jj = 1, ndime
            do kk = 1, ndime
               presu(ii, jj) = presu(ii, jj) - J * invF(ii, kk) &
                              * don * invF(jj, kk)
            end do
         end do
      end do
      
      PK2 = presu + PK2
      
   end subroutine sm400_donnan2_PK2

   subroutine sm400_donnan2_dpdds(pmate, F, J, dds)
      use def_solidz, only  :  parco_sld
   
      implicit none
      
      integer(ip), intent(in) :: pmate
      real(rp), intent(in)    :: F(ndime, ndime)
      real(rp), intent(in)    :: J
      real(rp), intent(inout) :: dds(ndime, ndime, ndime, ndime)

      integer(ip) :: ii, jj, kk, ll, mm, nn, qq, pp
      real(rp)    :: bidon
      real(rp)    :: A4(ndime, ndime, ndime, ndime), A2(ndime, ndime, ndime, ndime)
      real(rp)    :: tkron(ndime, ndime)
      real(rp)    :: invF(ndime, ndime)
      real(rp)    :: Rgas, Temp, tramp
      real(rp)    :: phi_0, c_F0, c_b
      real(rp)    :: c_F

      A2 = 0.0_rp
      A4 = 0.0_rp
      tkron = 0.0_rp
      do ii = 1, ndime
         tkron(ii, ii) = 1.0_rp
      end do

      Rgas    = parco_sld(3, pmate)    ! Ideal gas constant
      Temp    = parco_sld(4, pmate)    ! Absolute Temperature [K]
      c_F0    = parco_sld(5, pmate)    ! Initial fixed charge density {positive or negative, I believe [mm^3/mol]}
      c_b     = parco_sld(6, pmate)    ! Osmolarity bath
      tramp   = parco_sld(7, pmate)   ! Time ramp for application of osmotic pressure
      phi_0   = 1.0_rp - parco_sld(2, pmate)    ! Initial water content

      call invmtx(F, invF, bidon, ndime)

      c_F = phi_0 / (J - 1 + phi_0) * c_F0 ! Initial water content

      do ii = 1,  ndime; do jj = 1, ndime
         do kk = 1, ndime; do ll = 1, ndime
            A4(ii, jj, kk, ll) = A4(ii, jj, kk, ll) + &
                                 Rgas * Temp * (J * c_F ** 2 / ((J - 1.0_rp + phi_0) * sqrt(c_F ** 2 + c_b ** 2)) * &
                                 tkron(ii, jj) * tkron(kk, ll) + &
                                 (sqrt(c_F ** 2 + c_b ** 2) - c_b) * & 
                                 (tkron(ii, kk) * tkron(jj, ll) + tkron(ii, ll) * tkron(jj, kk) - tkron(ii, jj) * tkron(kk, ll)))
         end do; end do
      end do; end do

      do ii = 1, ndime; do jj = 1, ndime; do kk = 1, ndime; do ll = 1, ndime
         do mm = 1, ndime; do nn = 1, ndime; do pp = 1, ndime; do qq = 1, ndime
            A2(mm, nn, pp, qq) = A2(mm, nn, pp, qq) + &
                                 invF(ll, qq) * invF(kk, pp) * invF(jj, nn) * invF(ii, mm) * &
                                 A4(ii, jj, kk, ll)
            
         end do; end do; end do; end do
      end do; end do; end do; end do

      dds = dds + A2
   
   end subroutine sm400_donnan2_dpdds

   subroutine sm400_donnan_dpdds(C, J, prs, dds)
      ! UNUSED MODEL
      implicit none

      real(rp), intent(in) :: C(ndime, ndime)
      real(rp), intent(in) :: J, prs

      real(rp), intent(inout) :: dds(ndime, ndime, ndime, ndime)

      integer(ip) :: ii, jj, kk, ll
      real(rp)    :: invC(ndime, ndime)
      real(rp)    :: bidon

      call invmtx(C, invC, bidon, ndime)
      do ii = 1, ndime; do jj = 1, ndime; do kk = 1, ndime; do ll = 1, ndime
         dds(ii, jj, kk, ll) = dds(ii, jj, kk, ll) + &
         J * prs  * (invC(kk, ll) * invC(ii, jj) - &
         invC(ii, kk) * invC(ll, jj) - invC(ii, ll) * invC(kk, jj))
      end do; end do; end do; end do
   end subroutine sm400_donnan_dpdds

   pure function trace(A, nsize) result(trc)
      integer(ip), intent(in) :: nsize
      real(rp), intent(in)    :: A(nsize, nsize)
      real(rp)                :: trc

      integer(ip) :: ii
   
      trc = 0.0_rp
      do ii = 1, nsize
         trc = trc + A(ii, ii)
      end do
   end function trace

   pure function DET(a, nsize) result(deter)
      implicit none

      integer(ip),                       intent(in) :: nsize
      real(rp), dimension(nsize, nsize), intent(in) :: a

      real(rp)                                      :: deter
      real(rp)                                      :: t1, t2, t3

      select case( nsize )
   
      case(2_ip)
    
         deter = a(1,1)*a(2,2)-a(2,1)*a(1,2)
    
      case(3_ip)
         t1    = a(2,2)*a(3,3) - a(3,2)*a(2,3)
         t2    =-a(2,1)*a(3,3) + a(3,1)*a(2,3)
         t3    = a(2,1)*a(3,2) - a(3,1)*a(2,2)
         deter = a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3
      end select
      
   end function DET

   pure function invariant4(C, v0) result(res)

      implicit none

      real(rp), intent(in), dimension(3)    :: v0
      real(rp), intent(in), dimension(3, 3) :: C

      real(rp) :: aux(3), res

      integer(ip) :: ii, jj

      res = 0.0_rp
      aux = 0.0_rp

      do ii = 1, 3
         do jj = 1, 3
            aux(ii) = aux(ii) + C(ii, jj) * v0(jj)
         end do
      end do

      do ii = 1, 3
         res = res + aux(ii) * v0(ii)
      end do

   end function invariant4

   pure function Aaniso(k2, I4) result(A)
      
      implicit none

      real(rp), intent(in) :: k2, I4
      
      real(rp) :: A

      A = exp(k2 * (I4 - 1.0_rp) ** 2_ip)

   end function Aaniso

   pure function dAaniso(k2, A, C, v0) result(dAdC)
      
      implicit none

      real(rp), intent(in)                  :: k2, A
      real(rp), intent(in), dimension(3)    :: v0
      real(rp), intent(in), dimension(3, 3) :: C

      real(rp), dimension(3, 3) :: dAdC

      integer(ip) :: ii, jj
      real(rp)    :: I4

      I4 = invariant4(C, v0)

      do ii = 1, 3
         do jj = 1, 3
            dAdC(ii, jj) = 2 * k2 * A * (I4 - 1) * v0(ii) * v0(jj)
         end do
      end do
   end function dAaniso

end module mod_sld_stress_model_400
