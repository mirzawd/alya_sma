!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_stress_model_comput.f90
!> @author  Adria Quintanas (adria.quintanas@udg.edu)
!> @date    December, 2015
!>          - Subroutine written
!> @brief   ToolBox for stress models
!>
!> @details ToolBox for material models.
!>
!>          \verbatim
!>          This toolBox includes useful subroutines for finite element
!>          implementations in material models:
!>           - Calculation of strain and stress tensors.
!>           - Tensor to Voigt and Voigt to tensor conversions
!>           - Tensor rotation
!>           - Polar decomposition
!>           - Transport operators
!>
!>          In Alya the Voigt rules in 2D and 3D are the following:
!>
!>          2D Voigt rule             3D Voigt rule
!>          \sigma_{ij} \sigma_{a}    \sigma_{ij} \sigma_{a}
!>                  ij          a             ij          a
!>                  11          1             11          1
!>                  22          2             22          2
!>                  12          3             33          3
!>                                            23          4
!>                                            13          5
!>                                            12          6
!>          \endverbatim
!>
!>          References:\n
!>          T. Belytschko, W. K. Liu, B. Moran, K. I. Elkhodary
!>          Nonlinear Finite elements for Continua and Structures\n
!>          E. J. Barbero. Introduction to Composite Materials Design\n
!>
!>
!> @}
!------------------------------------------------------------------------------

module mod_sld_stress_model_comput
   ! ==============================================================================
   ! INIT
   !
   use def_kintyp, only: &
      ip, rp
   ! -----------------------------------------------------------------------------
   implicit none

   integer(ip), parameter :: SM_PULLBACK = 0
   integer(ip), parameter :: SM_PUSHFORWARD = 1

   !
   !=============================================================| init |=========
   !==============================================================================
   ! PUBLIC
   !
   public                                    :: &
      SM_strain_tensor, & !
      SM_stress_tensor, & !
      SM_stress_transport, & !
      SM_stiffness_transport, & !
      SM_tensor_to_voigt_second, & !
      SM_voigt_to_tensor_second, & !
      SM_tensor_to_voigt_fourth, & !
      SM_rotate_basis_creation, & !
      SM_rotate_tensor_second, & !
      SM_rotate_voigt_second, & !
      SM_rotate_voigt_fourth, &
      SM_polar_decomposition, &
      SM_PUSHFORWARD, &
      SM_PULLBACK
   !
   !=============================================================| public |=======
   !==============================================================================
   ! PRIVATE
   !
   private                                 :: &
      priv_transformation_matrix

   !
   !=============================================================| private |======
   !==============================================================================
   ! CONTAINS
   !
contains

   !------------------------------------------------------------------------------
   !> @author  Adria Quintanas (adria.quintanas@udg.edu)
   !> @date
   !> @brief
   !> @details
   !------------------------------------------------------------------------------
   subroutine SM_strain_tensor(task, gpgdi, gpgre)

      ! ---------------------------------------------------------------------------
      implicit none
      ! ---------------------------------------------------------------------------
      integer(ip), intent(in)              :: &
         task
      real(rp), intent(in)              :: &
         gpgdi(3, 3)                              !
      real(rp), intent(out)             :: &
         gpgre(3, 3)                              !
      ! ---------------------------------------------------------------------------
      integer(ip)                             :: &
         i, j, k                                 !
      real(rp)                                :: &
         auxS1, tkron(3, 3)
      ! ---------------------------------------------------------------------------
      !
      tkron = 0.0_rp
      do i = 1, 3
         tkron(i, i) = 1.0_rp
      end do
      !
      ! ---------------------------------------------------------------------------
      select case (task)
         !
      case (0_ip)
         !
         ! Infinitesimal strain tensor
         !
         gpgre = 0.0_rp
         do i = 1, 3
            do j = 1, 3
               gpgre(i, j) = gpgre(i, j) + 0.5_rp*(gpgdi(i, j) + gpgdi(j, i)) - tkron(i, j)
            end do
         end do
         !
      case (1_ip)
         !
         ! Green-Lagrange strain tensor
         !
         gpgre = 0.0_rp
         do i = 1, 3
            do j = 1, 3
               auxS1 = 0.0_rp
               do k = 1, 3
                  auxS1 = auxS1 + gpgdi(k, i)*gpgdi(k, j)
               end do
               gpgre(i, j) = 0.5_rp*(auxS1 - tkron(i, j))
            end do
         end do
         !
      case default
         return
      end select
      !
      return

   end subroutine SM_strain_tensor

   !------------------------------------------------------------------------------
   !> @author  Gerard Guillamet and Adria Quintanas
   !> @date    September, 2021
   !> @brief   Transport operators for stress measure
   !> @details Transport operators for stress measure (Bely's book)
   !>
   !>          Pull-back:    {S}_ij     = J F_{ik}^{-1} sigma_{kl} F_{lj}^{-T}
   !>          Push-forward: {sigma}_ij = J^{-1} F_{ik} S_{kl} F_{lj}^{T}
   !>
   !------------------------------------------------------------------------------
   subroutine SM_stress_transport(itask, ndofn, detF, F, Finv, sigma, stress)

      implicit none

      integer(ip), intent(in) :: itask                !< Task
      integer(ip), intent(in) :: ndofn                !< Dimensions
      real(rp), intent(in)    :: detF                 !< |F|      : Determinant of the deformation gradient tensor
      real(rp), intent(in)    :: F(ndofn, ndofn)      !< F        : Deformation gradient tensor
      real(rp), intent(in)    :: Finv(ndofn, ndofn)   !< F^{-1}   : Inverse of the deformation gradient tensor
      real(rp), intent(inout) :: sigma(ndofn, ndofn)  !< Cauchy stress tensor
      real(rp), intent(inout) :: stress(ndofn, ndofn) !< 2nd Piola-Kirchoff stress tensor

      integer(ip)                :: idime, jdime, kdime, ldime

      select case (itask)

      case (0_ip)
         !
         ! {S}_ij (From Cauchy to 2PK)
         !
         stress = 0.0_rp
         do ldime = 1, ndofn
            do kdime = 1, ndofn
               do jdime = 1, ndofn
                  do idime = 1, ndofn
                     stress(idime, jdime) = stress(idime, jdime) + detF*Finv(idime, kdime)*sigma(kdime, ldime)*Finv(jdime, ldime)
                  end do
               end do
            end do
         end do

      case (1_ip)
         !
         ! {sigma}_ij (From 2PK to Cauchy)
         !
         sigma = 0.0_rp
         do ldime = 1, ndofn
            do kdime = 1, ndofn
               do jdime = 1, ndofn
                  do idime = 1, ndofn
                     sigma(idime, jdime) = sigma(idime, jdime) + detF*F(idime, kdime)*stress(kdime, ldime)*F(jdime, ldime)
                  end do
               end do
            end do
         end do
         
      case (2_ip)
         !
         ! From Cauchy to nominal stress P
         !
         stress = 0.0_rp
         do jdime = 1, ndofn
            do idime = 1, ndofn
               stress(idime, jdime) = stress(idime, jdime) + detF*Finv(idime, jdime)*sigma(idime, jdime)
            end do
         end do
         
      case default
         return
      end select
      !
      return

   end subroutine SM_stress_transport

   !------------------------------------------------------------------------------
   !> @author  Gerard Guillamet and Adria Quintanas
   !> @date    October, 2018
   !> @brief   Transport operators for stiffness tensor
   !> @details Transport operators for stiffness tensor
   !>
   !>          Pull-back:    C_{IJKL} = F_{iI}^{-1} F_{jJ}^{-1} F_{kK}^{-1} F_{lL}^{-1} c_{ijkl}
   !>          Push-forward: c_{ijkl} = F_{iI} F_{jJ} F_{kK} F_{lL} C_{IJKL}
   !>
   !------------------------------------------------------------------------------
   subroutine SM_stiffness_transport(itask, ndofn, F, Finv, &
                                     c_sma, C_cap)

      implicit none

      integer(ip), intent(in)    :: itask                          !< Task
      integer(ip), intent(in)    :: ndofn                          !< Dimensions
      real(rp), intent(in)    :: F(ndofn, ndofn)                 !< F      : Deformation gradient tensor
      real(rp), intent(in)    :: Finv(ndofn, ndofn)              !< F^{-1} : Inverse of the deformation gradient tensor
      real(rp), intent(inout) :: C_cap(ndofn, ndofn, ndofn, ndofn) !< Cijkl  : 2nd elasticity tensor
      real(rp), intent(inout) :: c_sma(ndofn, ndofn, ndofn, ndofn) !< cijkl  : 2nd elasticity tensor

      integer(ip)                :: idime, jdime, kdime, ldime
      integer(ip)                :: I, J, K, L

      select case (itask)

      case (0_ip)
         !
         ! C_{IJKL}
         !
         C_cap = 0.0_rp
         do I = 1, ndofn
            do J = 1, ndofn
               do K = 1, ndofn
                  do L = 1, ndofn
                     do idime = 1, ndofn
                        do jdime = 1, ndofn
                           do kdime = 1, ndofn
                              do ldime = 1, ndofn
                                 C_cap(I, J, K, L) = C_cap(I, J, K, L) + &
                                                     Finv(idime, I)*Finv(jdime, J)* &
                                                     Finv(kdime, K)*Finv(ldime, L)*c_sma(idime, jdime, kdime, ldime)
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do

      case (1_ip)
         !
         ! c_{ijkl}
         !
         c_sma = 0.0_rp
         do idime = 1, ndofn
            do jdime = 1, ndofn
               do kdime = 1, ndofn
                  do ldime = 1, ndofn
                     do I = 1, ndofn
                        do J = 1, ndofn
                           do K = 1, ndofn
                              do L = 1, ndofn
                                 c_sma(idime, jdime, kdime, ldime) = c_sma(idime, jdime, kdime, ldime) + &
                                                                     F(idime, I)*F(jdime, J)* &
                                                                     F(kdime, K)*F(ldime, L)*C_cap(I, J, K, L)
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do

      case default
         return
      end select
      !
      return

   end subroutine SM_stiffness_transport

   !------------------------------------------------------------------------------
   !> @author  Adria Quintanas (adria.quintanas@udg.edu)
   !> @date
   !> @brief
   !> @details
   !------------------------------------------------------------------------------
   subroutine SM_stress_tensor(task, strain, stiff, stress)

      ! ---------------------------------------------------------------------------
      implicit none
      ! ---------------------------------------------------------------------------
      integer(ip), intent(in)              :: &
         task
      real(rp), intent(in)              :: strain(6)
      real(rp), intent(in)              :: stiff(6, 6) !
      real(rp), intent(out)             :: stress(6)
      ! ---------------------------------------------------------------------------
      integer(iP)                             :: &
         i, j                                    !
      !
      ! ---------------------------------------------------------------------------
      !
      select case (task)
         !
      case (0_ip)
         !
         ! 2nd Piola-Kirchhoff stress tensor
         !   {S}_i = [C]_ij*{E}_j   / (1)-(1,1), (2)-(2,2), (3)-(3,3)
         !                            (4)-(2,3), (5)-(1,3), (6)-(1,2)
         !
         stress = 0.0_rp
         do i = 1, 6
            do j = 1, 6
               stress(i) = stress(i) + stiff(i, j)*strain(j)
            end do
         end do
         !
      case default
         return
      end select
      !
      return

   end subroutine SM_stress_tensor

   !------------------------------------------------------------------------------
   !> @author  Adria Quintanas (adria.quintanas@udg.edu)
   !> @date
   !>          - Subrountine written
   !> @author  Gerard Guillamet
   !> @date    June, 2018
   !>          - 2-d compatibility
   !> @brief   Tensor second-order to Voigt form and viceversa
   !> @details
   !>          The Voigt rule depends on whether a tensor is a kinetic quantity
   !>          such as stress, or kinematic quantity, such as strain.
   !>
   !>          Kinetic Voigt rule (Stresses):
   !>          2D Voigt rule                 3D Voigt rule
   !>          +-         -+   +-        -+  +-         -+   +-        -+
   !>          |\sigma_{11}|   |\sigma_{1}|  |\sigma_{11}|   |\sigma_{1}|
   !>          |\sigma_{22}| = |\sigma_{2}|  |\sigma_{22}|   |\sigma_{2}|
   !>          |\sigma_{12}|   |\sigma_{3}|  |\sigma_{33}|   |\sigma_{3}|
   !>          +-         -+   +-        -+  |\sigma_{23}| = |\sigma_{4}|
   !>                                        |\sigma_{13}|   |\sigma_{5}|
   !>                                        |\sigma_{12}|   |\sigma_{6}|
   !>                                        +-              +-        -+
   !>          Kinematic Voigt rule (Strains):
   !>          2D Voigt rule                             3D Voigt rule
   !>          +-                -+   +-             -+  +-                -+   +-             -+
   !>          |  \varepsilon_{11}|   |\varepsilon_{1}|  |  \varepsilon_{11}|   |\varepsilon_{1}|
   !>          |  \varepsilon_{22}| = |\varepsilon_{2}|  |  \varepsilon_{22}|   |\varepsilon_{2}|
   !>          |2*\varepslion_{12}|   |\varepsilon_{3}|  |  \varepsilon_{33}|   |\varepsilon_{3}|
   !>          +-                -+   +-             -+  |2*\varepsilon_{23}| = |\varepsilon_{4}|
   !>                                                    |2*\varepsilon_{13}|   |\varepsilon_{5}|
   !>                                                    |2*\varepsilon_{12}|   |\varepsilon_{6}|
   !>                                                    +-                -+   +-             -+
   !------------------------------------------------------------------------------
   subroutine SM_tensor_to_voigt_second(ndofn, tensor, voigt)

      use def_solidz, only: nvoig_sld, nvgij_inv_sld
      ! ---------------------------------------------------------------------------
      implicit none
      ! ---------------------------------------------------------------------------
      integer(ip), intent(in)    :: ndofn               !< Dimensions
      real(rp), intent(in)       :: tensor(ndofn, ndofn) !< Tensor second-order
      real(rp), intent(inout)    :: voigt(nvoig_sld)    !< Voigt form
      !
      integer(ip)                :: &
         i, j, ivoig

      !
      ! Tensor to Voigt form (Kinematic rule)
      !
      voigt = 0.0_rp
      do i = 1, ndofn
         do j = 1, ndofn
            ivoig = nvgij_inv_sld(i, j)
            voigt(ivoig) = voigt(ivoig) + tensor(i, j)
         end do
      end do

   end subroutine SM_tensor_to_voigt_second

   subroutine SM_voigt_to_tensor_second(ndofn, tensor, voigt)

      use def_solidz, only: nvoig_sld, nvgij_inv_sld
      ! ---------------------------------------------------------------------------
      implicit none
      ! ---------------------------------------------------------------------------
      integer(ip), intent(in)    :: ndofn               !< Dimensions
      real(rp), intent(inout)    :: tensor(ndofn, ndofn) !< Tensor second-order
      real(rp), intent(in)       :: voigt(nvoig_sld)    !< Voigt form
      !   
      integer(ip)                :: &
         i, j, ivoig
      !
      ! ---------------------------------------------------------------------------
      !
      ! Voigt form to tensor (Kinetic rule)
      !
      tensor = 0.0_rp
      do i = 1, ndofn
         do j = 1, ndofn
            ivoig = nvgij_inv_sld(i, j)
            tensor(i, j) = voigt(ivoig)
         end do
      end do

   end subroutine SM_voigt_to_tensor_second

   !------------------------------------------------------------------------------
   !> @author  Adria Quintanas (adria.quintanas@udg.edu)
   !> @date
   !>          - Subrountine written
   !> @author  Gerard Guillamet
   !> @date    June, 2018
   !>          - 2-d compatibility
   !> @brief   Tensor fourth-order to Voigt form and viceversa
   !> @details
   !>          The mapping ij -> a and kl -> b is:
   !>           a = i*delta_{ij} + (1 - delta_{ij})*(vaux - i - j)
   !>           b = k*delta_{kl} + (1 - delta_{kl})*(vaux - k - l)
   !>          where,
   !>           vaux = 6 for 2-d problems
   !>           vaux = 9 for 3-d problems
   !------------------------------------------------------------------------------
   subroutine SM_tensor_to_voigt_fourth(ndofn, task, tensor, voigt)

      use def_solidz, only: nvoig_sld

      ! ---------------------------------------------------------------------------
      implicit none
      ! ---------------------------------------------------------------------------
      integer(ip), intent(in)    :: ndofn                           !< Dimensions
      integer(ip), intent(in)    :: task                            !< Task
      real(rp), intent(inout) :: tensor(ndofn, ndofn, ndofn, ndofn) !< Tensor fourth-order
      real(rp), intent(inout) :: voigt(nvoig_sld, nvoig_sld)      !< Voigt form
      !
      integer(ip)                :: &
         i, j, k, l, p, q, &
         tkron(ndofn, ndofn), &
         vaux
      !
      ! ---------------------------------------------------------------------------
      !
      ! Kronecker delta
      tkron = 0_ip
      do i = 1, ndofn
         tkron(i, i) = 1_ip
      end do
      ! Auxiliar variable
      vaux = 0_ip
      if (ndofn == 2_ip) then       ! 2-d
         vaux = 6_ip
      else if (ndofn == 3_ip) then  ! 3-d
         vaux = 9_ip
      end if
      !
      ! ---------------------------------------------------------------------------
      select case (task)
         !
      case (0_ip)
         !
         ! Tensor to Voigt form
         !
         voigt = 0.0_rp
         do l = 1, ndofn
            do k = 1, ndofn
               do j = 1, ndofn
                  do i = 1, ndofn
                     p = i*tkron(i, j) + (1_ip - tkron(i, j))*(vaux - i - j)
                     q = k*tkron(k, l) + (1_ip - tkron(k, l))*(vaux - k - l)
                     voigt(p, q) = tensor(i, j, k, l)
                  end do
               end do
            end do
         end do

      case (1_ip)
         !
         ! Voigt form to tensor
         !
         tensor = 0.0_rp
         do l = 1, ndofn
            do k = 1, ndofn
               do j = 1, ndofn
                  do i = 1, ndofn
                     p = i*tkron(i, j) + (1_ip - tkron(i, j))*(vaux - i - j)
                     q = k*tkron(k, l) + (1_ip - tkron(k, l))*(vaux - k - l)
                     tensor(i, j, k, l) = voigt(p, q)
                  end do
               end do
            end do
         end do

      case default
         return

      end select

      return

   end subroutine SM_tensor_to_voigt_fourth

   !------------------------------------------------------------------------------
   !> @author  Adria Quintanas (adria.quintanas@udg.edu)
   !> @date
   !> @brief   Material basis creation
   !> @details Material basis creation for 3d cases
   !>          Material axes are filled by rows
   !>          Let's assume,
   !>             Global system: 1-2-3
   !>             Material system: 1'-2'-3'
   !>             The point coordiantes with respect to the reference (Global) are:
   !>
   !>                    ^ x_2
   !>                    |
   !>          x_2' theta|
   !>              \     |     / x_1'
   !>               \    |    /                 [x']= [R]*[x]
   !>                \   |   /
   !>                 \  |  /  theta
   !>                  \ | /
   !>                   \|/                x_1
   !>                     ----------------- >
   !>
   !>               -    -     -                        -   -   -
   !>              | x'_1 |   |  cos(theta) sin(theta) 0 | | x_1 |
   !>              | x'_2 | = | -sin(theta) cos(theta) 0 |·| x_2 |
   !>              | x'_3 |   |          0          0  1 | | x_3 |
   !>               -    -     -                        -   -   -
   !>
   !>
   !>          Barbero, E.J. (1999). Introduction to composite materials design.
   !>          Taylor and Francis, Philadelphia (USA)
   !>          E.J. Barbero - Finite Element Analysis of Composite Materials
   !>
   !------------------------------------------------------------------------------

   subroutine SM_rotate_basis_creation(ielem, gpgdi, rotmat)

      use def_solidz, only: kfl_fiber_sld, kfl_rmate_sld
      use def_solidz, only: axis1_sld, axis2_sld, axis3_sld

      ! ---------------------------------------------------------------------------
      implicit none
      ! ---------------------------------------------------------------------------
      integer(ip), intent(in)              :: &
         ielem
      real(rp), intent(in)              :: &
         gpgdi(3, 3)                              ! Deformation gradient
      real(rp), intent(out)             :: &
         rotmat(3, 3)
      !
      integer(ip)                             :: &
         i, j, k
      real(rp)                                :: &
         rotrig(3, 3), U(3, 3), & ! Solid rigid rotation matrix
         rotmcs(3, 3)                             ! Material coordinate system rotation matrix
      !
      ! ---------------------------------------------------------------------------
      !
      ! Create material coordinate system (rotation matrix)
      !
      rotmcs = 0.0_rp
      rotrig = 0.0_rp
      if (kfl_fiber_sld > 3_ip) then
         do i = 1, 3
            rotmcs(1, i) = real(axis1_sld(i, ielem), rp)
            rotmcs(2, i) = real(axis2_sld(i, ielem), rp)
            rotmcs(3, i) = real(axis3_sld(i, ielem), rp)
         end do
      else
         do i = 1, 3
            rotmcs(i, i) = 1.0_rp
            rotrig(i, i) = 1.0_rp
         end do
      end if
      !
      ! Compute material rotation matrix
      !
      if (kfl_rmate_sld == 0_ip) then
         rotmat(:, :) = rotmcs(:, :)
      else
         !
         ! Accounting for the solid rigid rotation
         !
         call SM_polar_decomposition(1_ip, 3_ip, gpgdi(:, :), rotrig(:, :), U(:, :))
         do i = 1, 3
            do j = 1, 3
               rotmat(i, j) = 0.0_rp
               do k = 1, 3
                  rotmat(i, j) = rotmat(i, j) + rotmcs(i, k)*rotrig(j, k)
               end do
            end do
         end do
      end if

   end subroutine SM_rotate_basis_creation

   !------------------------------------------------------------------------------
   !> @author  Adria Quintanas (adria.quintanas@udg.edu)
   !> @date
   !> @brief
   !> @details
   !------------------------------------------------------------------------------
   subroutine SM_rotate_tensor_second(task, rotmat, tenGCS, tenMCS)

      ! ---------------------------------------------------------------------------
      implicit none
      ! ---------------------------------------------------------------------------
      integer(ip), intent(in)              :: &
         task
      real(rp), intent(in)              :: &
         rotmat(3, 3)
      real(rp), intent(inout)           :: &
         tenGCS(3, 3), tenMCS(3, 3)
      !
      integer(ip)                             :: &
         i, j, p, q
      !
      ! ---------------------------------------------------------------------------
      !
      select case (task)
         !
      case (0_ip)
         !
         ! From global to local coordinate system
         !     Alcs = Q Agcs Q^{T}
         !
         tenMCS = 0.0_rp
         do i = 1, 3
            do j = 1, 3
               do p = 1, 3
                  do q = 1, 3
                     tenMCS(i, j) = tenMCS(i, j) + rotmat(i, p)*rotmat(j, q)*tenGCS(p, q)
                  end do
               end do
            end do
         end do
         !
      case (1_ip)
         !
         ! From local to global coordinate system
         !     Agcs = Q^{T} Alcs Q
         tenGCS = 0.0_rp
         do i = 1, 3
            do j = 1, 3
               do p = 1, 3
                  do q = 1, 3
                     tenGCS(i, j) = tenGCS(i, j) + rotmat(i, p)*rotmat(j, q)*tenMCS(p, q)
                  end do
               end do
            end do
         end do
         !
      case default
         return
      end select

      return

   end subroutine SM_rotate_tensor_second

   !------------------------------------------------------------------------------
   !> @author  Adria Quintanas (adria.quintanas@udg.edu)
   !> @date    November, 2015
   !> @brief   This subroutine performs a transformation of a tensor in the Voigt
   !>          from a old basis to a new basis.
   !> @details This subroutine performs a transformation of a tensor
   !>            Strain from Local CSYS to Global CSYS:
   !>               E _{i} = T _{ij}^T*E'_{j}
   !>            Strain from Global CSYS to Local CSYS:
   !>               E'_{i} = Tg_{ij}  *E _{j}
   !>            Stress  from Local CSYS to Global CSYS:
   !>               S _{i} = Tg_{ij}^T*S'_{j}
   !>            Stress  from Global CSYS to Local CSYS:
   !>               S'_{i} = T _{ij}  *S _{j}
   !>
   !>         Reference: E.J. Barbero - Finite Element Analysis of Composite
   !>            Materials
   !------------------------------------------------------------------------------
   subroutine SM_rotate_voigt_second(task, rotmat, voigtGCS, voigtLCS)

      ! ---------------------------------------------------------------------------
      implicit none
      ! ---------------------------------------------------------------------------
      integer(ip), intent(in)              :: &
         task
      real(rp), intent(in)              :: &
         rotmat(3, 3)
      real(rp), intent(inout)           :: &
         voigtGCS(6), voigtLCS(6)
      !
      integer(ip)                             :: &
         i, j
      real(rp)                                :: &
         rotvoi(6, 6)
      !
      ! ---------------------------------------------------------------------------
      !
      ! Different cases
      !
      select case (task)
         !
      case (0_ip)
         !
         ! Strain from Local CSYS to Global CSYS:
         !    E _{i} = T _{ij}^T*E'_{j}
         !
         call priv_transformation_matrix(0_ip, rotmat, rotvoi)
         !
         voigtGCS = 0.0_rp
         do i = 1, 6
            do j = 1, 6
               voigtGCS(i) = voigtGCS(i) + rotvoi(j, i)*voigtLCS(j)
            end do
         end do
         !
      case (1_ip)
         !
         ! Strain from Global CSYS to Local CSYS:
         !    E'_{i} = Tg_{ij}  *E _{j}
         !
         call priv_transformation_matrix(1_ip, rotmat, rotvoi)
         !
         voigtLCS = 0.0_rp
         do i = 1, 6
            do j = 1, 6
               voigtLCS(i) = voigtLCS(i) + rotvoi(i, j)*voigtGCS(j)
            end do
         end do
         !
      case (2_ip)
         !
         ! Stress  from Local CSYS to Global CSYS:
         !    S _{i} = Tg_{ij}^T*S'_{j}
         !
         call priv_transformation_matrix(1_ip, rotmat, rotvoi)
         !
         voigtGCS = 0.0_rp
         do i = 1, 6
            do j = 1, 6
               voigtGCS(i) = voigtGCS(i) + rotvoi(j, i)*voigtLCS(j)
            end do
         end do
         !
      case (3_ip)
         !
         ! Stress  from Global CSYS to Local CSYS:
         !    S'_{i} = T _{ij}  *S _{j}
         !
         call priv_transformation_matrix(0_ip, rotmat, rotvoi)
         !
         voigtLCS = 0.0_rp
         do i = 1, 6
            do j = 1, 6
               voigtLCS(i) = voigtLCS(i) + rotvoi(i, j)*voigtGCS(j)
            end do
         end do
         !
      case default
         return
      end select

      return

   end subroutine SM_rotate_voigt_second

   !------------------------------------------------------------------------------
   !> @author  Adria Quintanas (adria.quintanas@udg.edu)
   !> @date    November, 2015
   !> @brief   This subroutine performs a transformation of a tensor in the Voigt
   !>          from a old basis to a new basis.
   !> @details This subroutine performs a transformation of a tensor
   !>            Strain from Local CSYS to Global CSYS:
   !>               E _{i} = T _{ij}^T*E'_{j}
   !>            Strain from Global CSYS to Local CSYS:
   !>               E'_{i} = Tg_{ij}  *E _{j}
   !>            Stress  from Local CSYS to Global CSYS:
   !>               S _{i} = Tg_{ij}^T*S'_{j}
   !>            Stress  from Global CSYS to Local CSYS:
   !>               S'_{i} = T _{ij}  *S _{j}
   !>
   !>         Reference: E.J. Barbero - Finite Element Analysis of Composite
   !>            Materials
   !------------------------------------------------------------------------------
   subroutine SM_rotate_voigt_fourth(task, rotmat, voigtGCS, voigtLCS)

      ! ---------------------------------------------------------------------------
      implicit none
      ! ---------------------------------------------------------------------------
      integer(ip), intent(in)              :: &
         task
      real(rp), intent(in)              :: &
         rotmat(3, 3)
      real(rp), intent(inout)           :: &
         voigtGCS(6, 6), voigtLCS(6, 6)
      !
      integer(ip)                             :: &
         i, j, p, q
      real(rp)                                :: &
         rotvoi(6, 6)
      !
      ! ---------------------------------------------------------------------------
      !
      ! Different cases
      !
      select case (task)
         !
      case (0_ip)
         !
         ! Compliance tensor from LCSYS to GCSYS
         !    H_ {ij} = T_{pi}*T_{qj}*H'_{pq}
         !
         call priv_transformation_matrix(0_ip, rotmat, rotvoi)
         !
         voigtGCS = 0.0_rp
         do i = 1, 6
            do j = 1, 6
               do p = 1, 6
                  do q = 1, 6
                     voigtGCS(i, j) = voigtGCS(i, j) + rotvoi(p, i)*rotvoi(q, j)*voigtLCS(p, q)
                  end do
               end do
            end do
         end do
         !
      case (1_ip)
         !
         ! Compliance tensor from GCSYS to LCSYS
         !    H'_{ij} = Tg_{ip}*H _{pq}*Tg_{jq}
         !
         call priv_transformation_matrix(1_ip, rotmat, rotvoi)
         !
         voigtLCS = 0.0_rp
         do i = 1, 6
            do j = 1, 6
               do p = 1, 6
                  do q = 1, 6
                     voigtLCS(i, j) = voigtLCS(i, j) + rotvoi(i, p)*rotvoi(j, q)*voigtGCS(p, q)
                  end do
               end do
            end do
         end do
         !
      case (2_ip)
         !
         ! Stiffness tensor from LCSYS to GCSYS
         !    C _{ij} = Tg_{pi}*C'_{pq}*Tg_{qj}
         !
         call priv_transformation_matrix(1_ip, rotmat, rotvoi)
         !
         voigtGCS = 0.0_rp
         do i = 1, 6
            do j = 1, 6
               do p = 1, 6
                  do q = 1, 6
                     voigtGCS(i, j) = voigtGCS(i, j) + rotvoi(p, i)*rotvoi(q, j)*voigtLCS(p, q)
                  end do
               end do
            end do
         end do
         !
      case (3_ip)
         !
         ! Stiffness tensor from GCSYS to LCSYS
         !    C'_ {ij} = T _{ip}*T _{jq}*C_{pq}
         !
         call priv_transformation_matrix(0_ip, rotmat, rotvoi)
         !
         voigtLCS = 0.0_rp
         do i = 1, 6
            do j = 1, 6
               do p = 1, 6
                  do q = 1, 6
                     voigtLCS(i, j) = voigtLCS(i, j) + rotvoi(p, i)*rotvoi(j, q)*voigtGCS(p, q)
                  end do
               end do
            end do
         end do
         !
      case default
         return
      end select

      return

   end subroutine SM_rotate_voigt_fourth

   ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !------------------------------------------------------------------------------
   !> @author  Adria Quintanas
   !> @date
   !> @brief   Transformation matrix forms for stress and strain tensors
   !> @details Transformation matrix forms for stress and strain tensors
   !>
   !>          E.J. Barbero - Finite Element Analysis of Composite Materials
   !>          Barbero, E.J. (1999). Introduction to composite materials design.
   !>          Taylor and Francis, Philadelphia (USA)
   !------------------------------------------------------------------------------

   subroutine priv_transformation_matrix(task, rotmat, rotvoi)

      ! ---------------------------------------------------------------------------
      implicit none
      ! ---------------------------------------------------------------------------
      integer(ip), intent(in)              :: &
         task
      real(rp), intent(in)              :: &
         rotmat(3, 3)
      real(rp), intent(out)             :: &
         rotvoi(6, 6)
      !
      ! ---------------------------------------------------------------------------
      !
      rotvoi = 0.0_rp
      !
      ! TRANSFORMATION MATRIX
      !    [T] and [Tg] are almost equal, only change some positions, then, the common part
      !    is defined here.

      rotVoi(1, 1) = rotmat(1, 1)*rotmat(1, 1)
      rotVoi(1, 2) = rotmat(1, 2)*rotmat(1, 2)
      rotVoi(1, 3) = rotmat(1, 3)*rotmat(1, 3)
      rotVoi(1, 4) = rotmat(1, 2)*rotmat(1, 3)
      rotVoi(1, 5) = rotmat(1, 1)*rotmat(1, 3)
      rotVoi(1, 6) = rotmat(1, 1)*rotmat(1, 2)
      rotVoi(2, 1) = rotmat(2, 1)*rotmat(2, 1)
      rotVoi(2, 2) = rotmat(2, 2)*rotmat(2, 2)
      rotVoi(2, 3) = rotmat(2, 3)*rotmat(2, 3)
      rotVoi(2, 4) = rotmat(2, 2)*rotmat(2, 3)
      rotVoi(2, 5) = rotmat(2, 1)*rotmat(2, 3)
      rotVoi(2, 6) = rotmat(2, 1)*rotmat(2, 2)
      rotVoi(3, 1) = rotmat(3, 1)*rotmat(3, 1)
      rotVoi(3, 2) = rotmat(3, 2)*rotmat(3, 2)
      rotVoi(3, 3) = rotmat(3, 3)*rotmat(3, 3)
      rotVoi(3, 4) = rotmat(3, 2)*rotmat(3, 3)
      rotVoi(3, 5) = rotmat(3, 1)*rotmat(3, 3)
      rotVoi(3, 6) = rotmat(3, 1)*rotmat(3, 2)
      rotVoi(4, 1) = rotmat(2, 1)*rotmat(3, 1)
      rotVoi(4, 2) = rotmat(2, 2)*rotmat(3, 2)
      rotVoi(4, 3) = rotmat(2, 3)*rotmat(3, 3)
      rotVoi(4, 4) = rotmat(2, 2)*rotmat(3, 3) + rotmat(2, 3)*rotmat(3, 2)
      rotVoi(4, 5) = rotmat(2, 1)*rotmat(3, 3) + rotmat(2, 3)*rotmat(3, 1)
      rotVoi(4, 6) = rotmat(2, 1)*rotmat(3, 2) + rotmat(2, 2)*rotmat(3, 1)
      rotVoi(5, 1) = rotmat(1, 1)*rotmat(3, 1)
      rotVoi(5, 2) = rotmat(1, 2)*rotmat(3, 2)
      rotVoi(5, 3) = rotmat(1, 3)*rotmat(3, 3)
      rotVoi(5, 4) = rotmat(1, 2)*rotmat(3, 3) + rotmat(1, 3)*rotmat(3, 2)
      rotVoi(5, 5) = rotmat(1, 1)*rotmat(3, 3) + rotmat(1, 3)*rotmat(3, 1)
      rotVoi(5, 6) = rotmat(1, 1)*rotmat(3, 2) + rotmat(1, 2)*rotmat(3, 1)
      rotVoi(6, 1) = rotmat(1, 1)*rotmat(2, 1)
      rotVoi(6, 2) = rotmat(1, 2)*rotmat(2, 2)
      rotVoi(6, 3) = rotmat(1, 3)*rotmat(2, 3)
      rotVoi(6, 4) = rotmat(1, 2)*rotmat(2, 3) + rotmat(1, 3)*rotmat(2, 2)
      rotVoi(6, 5) = rotmat(1, 1)*rotmat(2, 3) + rotmat(1, 3)*rotmat(2, 1)
      rotVoi(6, 6) = rotmat(1, 1)*rotmat(2, 2) + rotmat(1, 2)*rotmat(2, 1)
      !
      ! T or Tg
      !
      select case (task)

      case (0_ip)
         !
         ! Transformation matrix [T]
         !
         rotVoi(1, 4) = rotVoi(1, 4)*2.0_rp
         rotVoi(1, 5) = rotVoi(1, 5)*2.0_rp
         rotVoi(1, 6) = rotVoi(1, 6)*2.0_rp
         rotVoi(2, 4) = rotVoi(2, 4)*2.0_rp
         rotVoi(2, 5) = rotVoi(2, 5)*2.0_rp
         rotVoi(2, 6) = rotVoi(2, 6)*2.0_rp
         rotVoi(3, 4) = rotVoi(3, 4)*2.0_rp
         rotVoi(3, 5) = rotVoi(3, 5)*2.0_rp
         rotVoi(3, 6) = rotVoi(3, 6)*2.0_rp

      case (1_ip)
         !
         ! Transformation matrix [Tg]
         !
         rotVoi(4, 1) = rotVoi(4, 1)*2.0_rp
         rotVoi(4, 2) = rotVoi(4, 2)*2.0_rp
         rotVoi(4, 3) = rotVoi(4, 3)*2.0_rp
         rotVoi(5, 1) = rotVoi(5, 1)*2.0_rp
         rotVoi(5, 2) = rotVoi(5, 2)*2.0_rp
         rotVoi(5, 3) = rotVoi(5, 3)*2.0_rp
         rotVoi(6, 1) = rotVoi(6, 1)*2.0_rp
         rotVoi(6, 2) = rotVoi(6, 2)*2.0_rp
         rotVoi(6, 3) = rotVoi(6, 3)*2.0_rp

      end select

   end subroutine priv_transformation_matrix

   !------------------------------------------------------------------------------
   !> @author  Adria Quintanas
   !> @author  Gerard Guillamet
   !> @date
   !> @brief   Polar decomposition
   !> @details This subroutine performs the polar decomposition of the
   !>          deformation gradient F. The deformation gradient can be writter
   !>          either R*U or V*R. In each case, R is the rotation matrix, and
   !>          U and V are symmetric matrices describing the deformations
   !------------------------------------------------------------------------------

   subroutine SM_polar_decomposition(itask, ndofn, F, R, U, V, TOLERANCE)

      use mod_maths, only : maths_invert_matrix
      use mod_maths, only : maths_eigen_symmetric_matrix

      implicit none

      integer(ip), intent(in)            :: itask          !< Choose method
      integer(ip), intent(in)            :: ndofn          !< Dimensions
      real(rp), intent(in)            :: F(ndofn, ndofn) !< Deformation gradient
      real(rp), intent(out)           :: R(ndofn, ndofn) !< Rotation matrix
      real(rp), intent(out)           :: U(ndofn, ndofn) !< Stretch tensor
      real(rp), intent(out), optional :: V(ndofn, ndofn) !< Stretch tensor
      real(rp), intent(in), optional :: TOLERANCE

      integer(ip)              :: idime, jdime, kdime, ldime, nrot
      real(rp)                 :: C(ndofn, ndofn), omega(ndofn), Q(ndofn, ndofn)
      real(rp)                 :: Uinv(ndofn, ndofn), Utilde(ndofn, ndofn), detU
      real(rp)                 :: Rinv(ndofn, ndofn), detR
      real(rp)                 :: Rold(ndofn, ndofn)
      real(rp)                 :: xxnorm2, yynorm2
      real(rp), save           :: err
      real(rp)                 :: my_toler

      select case (itask)

      case (1_ip)

         !-------------------------------------------------------------------
         !
         ! Finding Square Roots of (Symmetric) Matrices
         !
         !-------------------------------------------------------------------
         !
         ! Right Cauchy-Green Deformation tensor, C = F^T * F
         !
         do idime = 1, ndofn
            do jdime = 1, ndofn
               C(idime, jdime) = 0.0_rp
               do kdime = 1, ndofn
                  C(idime, jdime) = C(idime, jdime) + F(kdime, idime)*F(kdime, jdime)
               end do
            end do
         end do
         !
         ! Compute eigenvalues and eigenvectors
         !
         call maths_eigen_symmetric_matrix(ndofn,C(:,:),omega(:),Q(:,:),maxit=100_ip,toler=1.0e-06_rp)
         !
         ! Calculate the principal values of U'
         !
         Utilde = 0.0_rp
         do idime = 1, ndofn
            Utilde(idime, idime) = sqrt(omega(idime))
         end do
         !
         ! Calculate complete tensor U = Q*U*Q^{T}
         !
         do idime = 1, ndofn
            do jdime = 1, ndofn
               U(idime, jdime) = 0.0_rp
               do kdime = 1, ndofn
                  do ldime = 1, ndofn
                     U(idime, jdime) = U(idime, jdime) + Q(idime, kdime)*Utilde(kdime, ldime)*Q(jdime, ldime)
                  end do
               end do
            end do
         end do
         !
         ! Calculate inverse of U
         !
         call maths_invert_matrix(ndofn, U, detU, Uinv)
         !
         ! Calculate R = F*U^{-1}
         !
         do idime = 1, ndofn
            do jdime = 1, ndofn
               R(idime, jdime) = 0.0_rp
               do kdime = 1, ndofn
                  R(idime, jdime) = R(idime, jdime) + F(idime, kdime)*Uinv(kdime, jdime)
               end do
            end do
         end do
         !
         ! Calculate V = F*R^{-1}
         !
         if (present(V)) then
            call maths_invert_matrix(ndofn, R, detR, Rinv)
            do idime = 1, ndofn
               do jdime = 1, ndofn
                  V(idime, jdime) = 0.0_rp
                  do kdime = 1, ndofn
                     V(idime, jdime) = V(idime, jdime) + F(idime, kdime)*Rinv(kdime, jdime)
                  end do
               end do
            end do
         end if

      case (2_ip)

         !-------------------------------------------------------------------
         !
         ! Iterative Algorithm
         !
         !-------------------------------------------------------------------

         if (present(TOLERANCE)) then
            my_toler = TOLERANCE
         else
            my_toler = epsilon(1.0_rp)
         end if

         ! Initial solution
         Rold(:, :) = F(:, :)
         err = huge(1.0_rp)
         do while (err > my_toler)
            call maths_invert_matrix(ndofn, Rold, detR, Rinv)
            R = 0.5_rp*(Rold + transpose(Rinv))
            ! Error
            xxnorm2 = 0.0_rp
            yynorm2 = 0.0_rp
            do idime = 1, ndofn
               do jdime = 1, ndofn
                  xxnorm2 = xxnorm2 + R(idime, jdime)**2
                  yynorm2 = yynorm2 + Rold(idime, jdime)**2
               end do
            end do
            xxnorm2 = sqrt(xxnorm2)
            yynorm2 = sqrt(yynorm2)
            err = abs(xxnorm2 - yynorm2)
            ! Update value
            Rold = R
         end do
         !
         ! Calculate U
         !
         U = 0.0_rp
         !
         ! Calculate V = F*R^{-1}
         !
         if (present(V)) then
            call maths_invert_matrix(ndofn, R, detR, Rinv)
            do idime = 1, ndofn
               do jdime = 1, ndofn
                  V(idime, jdime) = 0.0_rp
                  do kdime = 1, ndofn
                     V(idime, jdime) = V(idime, jdime) + F(idime, kdime)*Rinv(kdime, jdime)
                  end do
               end do
            end do
         end if

      end select

   end subroutine SM_polar_decomposition

   !=============================================================| contains |=====

end module mod_sld_stress_model_comput
