!------------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_stress_model_160.f90
!> @author  Adria Quintanas (adria.quintanas@udg.edu), Waleed Mirza (waleed.mirza889@gmail.com)
!>          Ref_A : 10.1016/S0045-7825(96)01147-4
!>          Ref_B : 10.1016/S0749-6419(00)00050-4
!>          Ref_D : 10.1007/s00466-017-1518-9
!> @}
!------------------------------------------------------------------------------
module sld_stress_model_160
   ! --------------------------------------------------------------------------
   use def_kintyp, only                :  ip, rp, lg
   use def_solidz

   !-------------------------------------
   implicit none
   !-------------------------------------
   real(rp), parameter                 :: zero = 0.0_rp
   logical(lg)                         :: go_NR
   !-------------------------------------
   !public                              :: eval_sld_stress_model_160
   !public                              :: goNR
   !-------------------------------------
   !private
   ! --------------------------------------------------------------------------

   contains

   ! --------------------------------------------------------------------------
   !> @author  Adria Quintanas
   !> @date
   !> @brief
   !> @details
   ! --------------------------------------------------------------------------
   subroutine eval_sld_stress_model_160( &
      pgaus, props, gpgdi, gptmp, gpsdv, gpstr, gpdds, flagt)
      !-------------------------------------
      implicit none
      !-------------------------------------
      integer(ip),   intent(in)           :: pgaus
      real(rp),      intent(in)           :: props(:)
      real(rp),      intent(in)           :: gpgdi(:,:,:)
      real(rp),      intent(in)           :: gptmp(:)
      real(rp),      intent(inout)        :: gpsdv(:,:,:)
      real(rp),      intent(out)          :: gpstr(:,:,:)  
      real(rp),      intent(out)          :: gpdds(:,:,:,:,:)
      integer(ip),   intent(in)           :: flagt
      !-------------------------------------
      real(rp),      dimension(18)        :: gp_props
      
      logical(lg)                         :: flow_rate_exp
      real(rp),      dimension(3)         :: gp_SDV
      real(rp),      dimension(3, 3)       :: gp_F, gp_F_inv, gp_F_per, aux1, aux2
      real(rp)                            :: gp_tmp
      real(rp),      dimension(3, 3)       :: gp_str, gp_str_per
      real(rp),      dimension(3, 3, 3, 3)   :: gp_dds
      integer(ip)                         :: igaus, i, j, k, l, m, n, p, q
      integer(ip)                         :: a, b, c, d
      !-------------------------------------
      real(rp), parameter                 :: eps = 1.0e-5_rp
      real(rp), parameter, dimension(3, 3):: car_basis = reshape((/           &
      &                                          1.0_rp, 0.0_rp, 0.0_rp,  &
      &                                          0.0_rp, 1.0_rp, 0.0_rp,  &
      &                                          0.0_rp, 0.0_rp, 1.0_rp   &
      &                                          /), (/3, 3/) ) 
      !-------------------------------------
      !
      ! Get properties
      !

       
      

      
      gp_props(1)  = 5.0e+4_rp  !E
      gp_props(2)  = 0.0_rp     !v
      gp_props(3)  = 0.07_rp     !eL
      gp_props(4)  = 500.0_rp    !Ss_As
      gp_props(5)  = 500.0_rp    !Af_As
      gp_props(6)  = 0.0_rp   !CAS
      gp_props(7)  = 00.0_rp  !Tf_AS
      gp_props(8)  = 00.0_rp  !Ts_AS
      gp_props(9)  = 700.0_rp  !no one is using this
      gp_props(10) = 300_rp ! \sigmas_sa
      gp_props(11) = 300_rp ! \sigmaf_sa
      gp_props(12) = 0.0_rp  !C_SA
      gp_props(13) = 0.0_rp  !Tf_SA
      gp_props(14) = 0.0_rp  ! Ts_SA
      gp_props(15) = 0.      !alpha
      gp_props(16) = 1.0_rp  !beta_As
      gp_props(17) = 1.0_rp  !beta_SA
      gp_props(18) = 1.0_rp  !beta_SA
      flow_rate_exp= .False.

      !
      ! Loop over gauss points
      !
      loop_gauss_points: do igaus = 1, pgaus


         ! Assign current gp variables
         gp_F     = gpgdi(:,:,igaus)
         gp_F_inv = INV3x3(gpgdi(:,:,igaus))
         gp_tmp   = 0.0_rp  ! gptmp(igaus)

         ! Get state variables
         gp_SDV(1:3) = gpsdv(1:3, igaus, 2)

         ! Get stress
         call compute_kirchoff_stress(props(:), gp_F(:,:), gp_tmp, gp_SDV(:), gp_str(:,:))

         ! Set state variables
         gpsdv(1:3, igaus, 1) = gp_SDV(1:3)

         ! Get consistent/algorithmic tangent numerically
         ! { Ref_C : Eq. 9-17 }
         if (flagt == 1_ip) then  
            do i = 1, 3
               do j = 1, 3

                  ! Compute perturbated deformation gradient
                  gp_F_per(:,:) = gp_F(:,:) + 0.5_rp*eps*(    &
                  !&     OUTPROD_v3(car_basis(:,i), car_basis(:,j)) + &
                  !&     OUTPROD_v3(car_basis(:,j), car_basis(:,i)) )
                  &     matmul(OUTPROD_v3(car_basis(:,i), car_basis(:,j)), gp_F(:,:)) + &
                  &     matmul(OUTPROD_v3(car_basis(:,j), car_basis(:,i)), gp_F(:,:)) )

                  ! Recorver state variables
                  gp_SDV(1:3) = gpsdv(1:3, igaus, 2)

                  ! Compute perturbated Kirchhoff stress
                  call compute_kirchoff_stress(props(:), gp_F_per(:,:), &
                  &     gp_tmp, gp_SDV(:), gp_str_per(:,:))

                  ! Compute algorithmic tangent approximation
                  gp_dds(i, j, :,:) = (gp_str_per(:,:) - gp_str(:,:))/eps

               enddo
            enddo
            ! correction {Ref_D: Eq. 68}
            do a = 1, 3
               do b = 1, 3
                  do c = 1, 3
                     do d = 1, 3
                        gp_dds(a, b, c, d) = gp_dds(a, b, c, d) - 0.5_rp*(&
                        !&    car_basis(j, l)*gp_str(i, k) + car_basis(i, l)*gp_str(j, k) + &
                        !&    car_basis(j, k)*gp_str(i, l) + car_basis(i, k)*gp_str(j, l))
                        &    car_basis(a, c)*gp_str(d, b) + car_basis(c, b)*gp_str(a, d) + &
                        &    car_basis(a, d)*gp_str(c, b) + car_basis(d, b)*gp_str(a, c))
                     enddo
                  enddo
               enddo
            enddo
         end if 

         ! From kirchoff to Second-Piola
         ! - stress
         gpstr(:,:,igaus) = 0.0_rp
         do i = 1, 3
            do j = 1, 3
               do k = 1, 3
                  do l = 1, 3
                     gpstr(i, j, igaus) = gpstr(i, j, igaus) + gp_F_inv(i, k)*gp_F_inv(j, l)*gp_str(k, l)
                  enddo
               enddo
            enddo
         enddo
         
         if (flagt == 1_ip) then 
            ! - algorithmic tangent
            gpdds(:,:,:,:,igaus) = 0.0_rp
            do i = 1, 3
               do j = 1, 3
                  do k = 1, 3
                     do l = 1, 3
                        do m = 1, 3
                           do n = 1, 3
                              do p = 1, 3
                                 do q = 1, 3
                                   gpdds(i, j, k, l, igaus) = gpdds(i, j, k, l, igaus) + &
                                   &     gp_F_inv(i, m)*gp_F_inv(j, n)*gp_F_inv(k, p)*gp_F_inv(l, q)*gp_dds(m, n, p, q)
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         end if
      end do loop_gauss_points

      !-------------------------------------
   end subroutine eval_sld_stress_model_160

   subroutine is_not_valid(x, messa)
      use def_master, only: isnain
      implicit none
      real(rp), intent(in):: x
      character(len=*), intent(in):: messa
      if( isnain(x) ) call runend("SMA 160: SOME VALUE IS NAN IN "//trim(messa))

   end subroutine is_not_valid

   ! --------------------------------------------------------------------------
   !> @author  Adria Quintanas
   !> @date
   !> @brief
   !> @details
   ! --------------------------------------------------------------------------
   subroutine compute_kirchoff_stress( PROPS, GRADDEF, TMP, SDV, STR)
      use def_master, only: isnain, itinn, modul
      !-------------------------------------
      implicit none
      !-------------------------------------
      real(rp),      intent(in)           :: PROPS(18)
      real(rp),      intent(in)           :: GRADDEF(3, 3)
      real(rp),      intent(in)           :: TMP
      real(rp),      intent(inout)        :: SDV(3)
      real(rp),      intent(out)          :: STR(3, 3)
      !-------------------------------------
      real(rp)                            :: E, v, K, G, eL, alfa
      real(rp)                            :: Ss_AS, Sf_AS, Sn_AS, Ss_SA, Sf_SA
      real(rp)                            :: Ts_AS, Tf_AS, Ts_SA, Tf_SA
      real(rp)                            :: C_AS, C_SA
      real(rp)                            :: Rs_AS, Rf_AS, Rs_SA, Rf_SA
      real(rp)                            :: B_AS, B_SA
      real(rp)                            :: G_1
      real(rp), dimension(3, 3)           :: b
      real(rp), dimension(3, 3)           :: n_A
      real(rp), dimension(3)              :: h_A
      real(rp)                            :: h_1, h_2, h_3
      real(rp), dimension(3)              :: n_1, n_2, n_3
      real(rp)                            :: e_1, e_2, e_3, e_n
      real(rp)                            :: d_1, d_2, d_3
      real(rp)                            :: J_vol, J_dev, theta
      real(rp)                            :: xS, xS_0, A
      real(rp)                            :: F_AS, F_AS_0, F_AS_dot, Ff_AS, Fs_AS, H_AS
      real(rp)                            :: F_SA, F_SA_0, F_SA_dot, Ff_SA, Fs_SA, H_SA
      logical(lg)                         :: kf_AS, kf_SA, kf_SS
      real(rp)                            :: resi
      real(rp)                            :: RHS_SA, dh_SA
      real(rp)                            :: RHS_AS, dh_AS, SN  
      integer(ip)                         :: igaus, iter, i, j
      real(rp), dimension(2, 2)            :: NR_MAT
      real(rp), dimension(2)              :: NR_RHS
      real(rp), dimension(2)              :: NR_SOL
      real(rp)                            :: aux
      logical(lg)                         :: flow_rate_exp  ! TODO modify this
      !-------------------------------------
      integer(ip), parameter              :: iter_NR = 65_ip
      real(rp),    parameter              :: tole_NR = 1.0e-8_rp
      real(rp), dimension(3, 3), parameter:: I_3x3 = reshape((/           &
      &                                          1.0_rp, 0.0_rp, 0.0_rp,  &
      &                                          0.0_rp, 1.0_rp, 0.0_rp,  &
      &                                          0.0_rp, 0.0_rp, 1.0_rp   &
      &                                          /), (/3, 3/) )
      !-------------------------------------
      !
      ! Get material properties
      ! - direct properties
       E     = props(1)
       v     = props(2)
       eL    = props(3)
       Ss_AS = props(4)
       Sf_AS = props(5)
       C_AS  = props(6)
       Tf_AS = props(7)
       Ts_AS = props(8)
       Sn_AS = props(9)
       Ss_SA = props(10)
       Sf_SA = props(11)
       C_SA  = props(12)
       Tf_SA = props(13)
       Ts_SA = props(14)
       alfa  = props(15)
       B_AS  = props(16)
       B_SA  = props(17)
       flow_rate_exp = .FALSE. ! TODO Pass this as input argument. Exponential flow rate has been tests and it works.
      ! - indirect propertes
      K  =  E/(3.0_rp-6.0_rp*v)
      G  =  E/(2.0_rp+2.0_rp*v) 
      aux   = sqrt(2.0_rp/3.0_rp) + alfa
      Rf_AS = sf_AS*aux+C_AS*Tf_AS
      Rs_AS = ss_AS*aux+C_AS*Ts_AS
      Rf_SA = sf_SA*aux+C_SA*Tf_SA
      Rs_SA = ss_SA*aux+C_SA*Ts_SA
      G_1   = (2.0_rp*G + 9.0_rp*(alfa**2)*K)*eL

      ! Trial left Cauchy tensor (b)
      ! {Ref_A:: Eq. 19}
      b = matmul(GRADDEF(:,:), transpose(GRADDEF(:,:)))

      ! Spectral decompositions of b
      ! {Ref_A:: Eq. 21}
      ! - computation of the eigenvalues and eigenvectors
      call DSYEVQ3(b, n_A, h_A)
      ! - sort eigenvalues and eigenvectors in ascending order
      call eigsrt(h_A, n_A)
      ! - principal stretches
      h_1 = sqrt(h_A(1))
      h_2 = sqrt(h_A(2))
      h_3 = sqrt(h_A(3))
      ! - principal directions
      n_1 = n_A(:,1)
      n_2 = n_A(:,2)
      n_3 = n_A(:,3)

      ! Logarithmic strains
      ! {Ref_A:: Eq. 23-26}
      ! - volumetric contribution
      J_vol = h_1*h_2*h_3
      theta = log(J_vol)

      ! - deviatoric contribution
      J_dev = J_vol**(-1.0_rp/3.0_rp)
      e_1 = log(J_dev*h_1)
      e_2 = log(J_dev*h_2)
      e_3 = log(J_dev*h_3)
      e_n = sqrt(e_1**2 + e_2**2 + e_3**2)

      ! - auxiliar variable
      if( e_n > epsilon(0.0_rp) )then
         d_1 = e_1/e_n
         d_2 = e_2/e_n
         d_3 = e_3/e_n
      else
         d_1 = 0.0_rp
         d_2 = 0.0_rp
         d_3 = 0.0_rp
      endif

      ! Recover state variables at the gauss points
      xS_0   = SDV(1)
      F_AS_0 = SDV(2)
      F_SA_0 = SDV(3)

      ! Check transformation phases
      !   { Ref_A & Ref_B, see inside function but and max added for physical consistency
      !    when dealing with large increments }
      F_AS = eval_loading_function(K, G, alfa, C_AS, eL, theta, e_n, tmp, xS_0 )
      !write
      ! - austenite to martensite
      !   { Ref_B : Eq. 4-8 }
      Fs_AS = F_AS-Rs_AS
      Ff_AS = F_AS-Rf_AS
      F_AS_dot = F_AS-F_AS_0
      !   { Ref_B : Eq. 10, but modified for physical consistency
      !    when dealing with large increments }
      if( F_AS > F_AS_0 .and. Fs_AS > 0.0 .and. xS_0 < 1.0 )then
         kf_AS = .True.
         H_AS = 1.0_rp
      else
         kf_AS = .False.
         H_AS = 0.0_rp
      endif
      ! - martensite to austenite
      !   { Ref_A & Ref_B, see inside function but and max added for physical consistency
      !    when dealing with large increments }
      F_SA = eval_loading_function(K, G, alfa, C_SA, eL, theta, e_n, tmp, xS_0)
      !   { Ref_B : Eq. 13-17 }
      Fs_SA = F_SA-Rs_SA
      Ff_SA = F_SA-Rf_SA
      F_SA_dot = F_SA-F_SA_0
      !   { Ref_B : Eq. 11,  but modified for physical consistency
      !    when dealing with large increments }
      if( F_SA < F_SA_0 .and. Fs_SA < 0.0 .and. xS_0 > 0.0 )then
         kf_SA = .TRUE.
         H_SA = 1.0_rp
      else
         kf_SA = .False.
         H_SA = 0.0_rp
      endif

      ! - check consistency
      if( kf_SA .and. kf_AS )then
         write(*,*)'WARNING 1: SOMETHING IS NOT COHERENT CONCERNING kf_SA AND kf_AS'
      endif

      ! Integration of the time-discrete evolutionary equations using a
      ! returning map algorithm employing a Newton-Raphson scheme.
      go_NR = .true. ! TODO: Pass this as input parameter
      !Find the flow parameters numerically
      if (go_NR) then
         integration_martensite_variable: if( kf_SA .or. kf_AS )then

            ! Initialisation
            xS = xS_0
            dh_SA = 0.0_rp
            dh_AS = 0.0_rp
            iter = 0_ip
            if (flow_rate_exp) then
               ! exponential case
               RHS_AS =  H_AS*(xS-xS_0)*Ff_AS**2 - H_AS*B_AS*F_AS_dot*(1.0-xS)
               RHS_AS =  RHS_AS+H_SA *(xS-xS_0)*Ff_SA**2   - H_SA*B_SA*xS*F_SA_dot
            else
               ! Linear case
               RHS_AS = (F_AS-Rf_AS)*dh_AS+H_AS*F_AS_dot*(1.0-xS)
               RHS_SA = (F_SA-Rf_SA)*dh_SA-H_SA*F_SA_dot*xS
            endif

            resi = RHS_AS**2 + RHS_SA**2

            loop_newton_raphson: do while( iter < iter_NR .and. resi > tole_NR )
               ! Increase iteration number
               iter = iter+1_ip
               ! Compute the increment of martensite state variable
               ! - Assemble the Jacobian coefficients of the linearised system of equations
               if (flow_rate_exp) then
                  NR_MAT(1, 1) =  H_AS*(F_AS-Rf_AS)**2-2*G_1*(F_AS-Rf_AS)*(xS-xS_0)*H_AS+H_AS*B_AS*F_AS_dot+G_1*H_AS*B_AS*(1.0-xS)
                  NR_MAT(1, 1) = NR_MAT(1, 1) +  H_SA *(F_SA-Rf_SA)**2 -  2*G_1*H_SA *(xS-xS_0)*(F_SA-Rf_SA)
                  NR_MAT(1, 1) = NR_MAT(1, 1) -  H_SA*B_SA*F_SA_dot+G_1*H_SA*B_SA*xS

                  xS = xS -  RHS_AS/NR_MAT(1, 1)
              else
                  NR_MAT(1, 1) = -G_1*dh_AS - (H_AS*F_AS_dot+H_AS*(1.0_rp-xS)*G_1) + (F_AS-Rf_AS)
                  NR_MAT(1, 2) = -G_1*dh_AS - (H_AS*F_AS_dot+H_AS*(1.0_rp-xS)*G_1)
                  NR_MAT(2, 1) = -G_1*dh_SA - (H_SA*F_SA_dot-H_SA*xS*G_1)
                  NR_MAT(2, 2) = -G_1*dh_SA - (H_SA*F_SA_dot-H_SA*xS*G_1) + (F_SA-Rf_SA)

                  ! - Assemble the RHS
                  NR_RHS(1) = RHS_AS
                  NR_RHS(2) = RHS_SA
                  ! - Solve the system
                  ! TODO: change the matmul for implicit function
                  NR_SOL(:) = matmul(INV2x2(NR_MAT(:,:)), NR_RHS(:))
                  ! - Assign the solution to auxiliar variables
                  dh_AS = dh_AS-NR_SOL(1)
                  dh_SA = dh_SA-NR_SOL(2)
                  ! Update the martensite state variable
                  xS = xS+dh_SA+dh_AS
               endif
               ! Update the loading function variables
               ! - austenite to martensite
               F_AS = eval_loading_function(K, G, alfa, C_AS, eL, theta, e_n, tmp, xS)
               F_AS_dot = F_AS-F_AS_0
               ! - martensite to austenite
               F_SA = eval_loading_function(K, G, alfa, C_SA, eL, theta, e_n, tmp, xS)
               F_SA_dot = F_SA-F_SA_0
               ! Compute RHS of the linearised system of equations
               if (flow_rate_exp) then
                  ! exponential case
                  RHS_AS = H_AS*(xS-xS_0)*(F_AS-Rf_AS)**2-H_AS*B_AS*F_AS_dot*(1.0-xS)
                  RHS_AS = RHS_AS+H_SA*(xS-xS_0)*(F_SA-Rf_SA)**2  - H_SA*B_SA*xS*(F_SA-F_SA_0)
                  RHS_SA = 0
               else
                  ! Linear case
                  RHS_AS = (F_AS-Rf_AS)*dh_AS+H_AS*F_AS_dot*(1.0-xS)
                  RHS_SA = (F_SA-Rf_SA)*dh_SA-H_SA*F_SA_dot*xS
               endif

               ! Compute the norm of the residual for the tolerance
               resi = sqrt(RHS_AS**2 + RHS_SA**2)


            enddo loop_newton_raphson

            if( iter >= iter_NR )then
                print *,'WARNING 2 ::: Max iter_NR achieve'
                xS   = SDV(1)
                F_AS = SDV(2)
                F_SA = SDV(3)
                xS = xS_0
            endif

         else
            ! If there is not a phase transformation, so the previous solution is correct
            xS = xS_0
         endif integration_martensite_variable
      else  
         ! Anlytical expressions of flow parameter from Ref B
         if (kf_AS) then
      	    A=(1-xS_0)/((Rf_AS-F_AS_0) + (1-xS_0)*G_1)
         else if (kf_SA) then
            A=-xS_0/((Rf_SA-F_SA_0) -xS_0*G_1)
         else
            A=0
         end if
         xS = xS_0 +A*(F_AS-F_AS_0)  
      end if 
      
      ! Set the phyisical limit of the martensite variable
      xS = max(min(1.0_rp, xS), 0.0_rp)

      ! Compute nominal Kirchoff stresses
      STR(:,:) = &
      &  K*(theta-3.0_rp*alfa*eL*xS)*I_3x3(:,:) + &
      &  2.0_rp*G*(e_1-eL*xS*d_1)*OUTPROD_v3(n_1(:), n_1(:)) + &
      &  2.0_rp*G*(e_2-eL*xS*d_2)*OUTPROD_v3(n_2(:), n_2(:)) + &
      &  2.0_rp*G*(e_3-eL*xS*d_3)*OUTPROD_v3(n_3(:), n_3(:))

      !
      ! Set the current state variables
      !
      SDV(1) = xS
      SDV(2) = F_AS
      SDV(3) = F_SA
      !-------------------------------------
   end subroutine compute_kirchoff_stress

   pure function left_cauchy_strain_tensor(F) result(b)
      !-------------------------------------
      implicit none
      !-------------------------------------
      real(rp), dimension(3, 3), intent(in)  :: F
      real(rp), dimension(3, 3)              :: b
      !-------------------------------------
      b = matmul(F, transpose(F))
      !-------------------------------------
   end function left_cauchy_strain_tensor

   function eval_loading_function(K, G, alfa, C, eL, theta, enor, tmp, xS) result(F)
      ! Evaluate the loading function according to the failure criteria
      ! { Modification :
      !     Ref_B : Eq.  3 for austenire  to martensite:: AS
      !     Ref_B : Eq. 12 for martensite to austenite  :: SA
      !     Ref_B : Eq. 21 for martensite reorientation:: SS
      !   according to :
      !     Ref_A : Remark 2.5, i.e. Eqs. 33-37 }
      !-------------------------------------
      implicit none
      !-------------------------------------
      real(rp), intent(in)                 :: K, G, alfa, C, eL
      real(rp), intent(in)                 :: theta, enor, tmp
      real(rp), intent(in)                 :: xS
      real(rp)                             :: F
      !-------------------------------------
      real(rp)                             :: pres, trac, beta
      !------------------------------------

      pres = K*(theta-3.0_rp*alfa*eL*xS)
      trac = 2.0_rp*G*(enor-eL*xS)
      beta = eval_beta(eL, enor, xS)
      F = trac*beta+3.0*alfa*pres-C*tmp
      !-------------------------------------
   end function eval_loading_function

   pure function eval_beta(eL, enor, xS) result(beta)
      ! Evaluate the beta coefficient
      ! { Ref_A : Remark 2.5, i.e. Eqs. 48 }
      !-------------------------------------
      implicit none
      !-------------------------------------
      real(rp), intent(in)                 :: eL, enor, xS
      real(rp)                             :: beta
      !-------------------------------------
      if( enor-eL*xS > 0.0 )then
         beta = 1.0_rp
      else
         beta = -1.0_rp
      endif
      
      beta = 1.0_rp !FIXME Beta=-1 gives convergence issues in contradiction with the paper.
      !-------------------------------------
   end function eval_beta

   pure function INV2x2(M) result(R)
      ! TODO: Add some control to check for 0 determinant
      ! --------------------------------------
      implicit none
      ! --------------------------------------
      real(rp),    dimension(2, 2),     intent(in):: M
      real(rp),    dimension(2, 2)                 :: R
      real(rp)                                    :: det
      ! --------------------------------------
      det = M(1, 1)*M(2, 2) - M(1, 2)*M(2, 1)
      R(1, 1) =  M(2, 2)/det
      R(1, 2) = -M(1, 2)/det
      R(2, 1) = -M(2, 1)/det
      R(2, 2) =  M(1, 1)/det
      ! --------------------------------------
   end function INV2x2

   pure function INV3x3(M) result(Mi)
      ! Ref: (Eq. ) from ...
      ! -------------------------------
      implicit none
      ! -------------------------------
      real(rp), dimension(3, 3), intent(in):: M
      real(rp), dimension(3, 3)             :: Mi
      real(rp)                             :: det
      ! -------------------------------
      ! Determinant
      det = M(1, 1)*M(2, 2)*M(3, 3) + M(1, 3)*M(2, 1)*M(3, 2) + &
            M(3, 1)*M(1, 2)*M(2, 3) - M(3, 1)*M(2, 2)*M(1, 3) - &
            M(3, 3)*M(1, 2)*M(2, 1) - M(1, 1)*M(2, 3)*M(3, 2)
      ! Invert matrix
      Mi(1, 1) =  (M(2, 2)*M(3, 3) - M(3, 2)*M(2, 3))/det
      Mi(1, 2) = -(M(1, 2)*M(3, 3) - M(1, 3)*M(3, 2))/det
      Mi(1, 3) =  (M(1, 2)*M(2, 3) - M(2, 2)*M(1, 3))/det
      Mi(2, 1) = -(M(2, 1)*M(3, 3) - M(3, 1)*M(2, 3))/det
      Mi(2, 2) =  (M(1, 1)*M(3, 3) - M(1, 3)*M(3, 1))/det
      Mi(2, 3) = -(M(1, 1)*M(2, 3) - M(2, 1)*M(1, 3))/det
      Mi(3, 1) =  (M(2, 1)*M(3, 2) - M(3, 1)*M(2, 2))/det
      Mi(3, 2) = -(M(1, 1)*M(3, 2) - M(3, 1)*M(1, 2))/det
      Mi(3, 3) =  (M(1, 1)*M(2, 2) - M(1, 2)*M(2, 1))/det
      ! -------------------------------
   end function INV3x3

   pure function OUTPROD_v3(v, u) result(M)
      ! TODO: Add some control to check for 0 determinant
      ! --------------------------------------
      implicit none
      ! --------------------------------------
      real(rp),    dimension(3),       intent(in):: u
      real(rp),    dimension(3),       intent(in):: v
      real(rp),    dimension(3, 3)                 :: M
      ! --------------------------------------
      M(:,1) = (/ u(1)*v(1), u(1)*v(2), u(1)*v(3) /)
      M(:,2) = (/ u(2)*v(1), u(2)*v(2), u(2)*v(3) /)
      M(:,3) = (/ u(3)*v(1), u(3)*v(2), u(3)*v(3) /)
      ! --------------------------------------
   end function OUTPROD_v3

end module sld_stress_model_160









   

