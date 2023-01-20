!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_calculate_pp_tensors.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Compute some tensors required for postprocess
!> @details Compute some tensors required for postprocess
!> @}
!-----------------------------------------------------------------------

subroutine sld_calculate_pp_tensors()

   use def_kintyp, only: ip, rp
   use def_elmtyp
   use def_master
   use def_domain
   use mod_output_postprocess, only: output_postprocess_check_variable_postprocess
   use mod_output_postprocess, only: output_postprocess_check_variable_witness
   use mod_output_postprocess, only: output_postprocess_check_variable_node_sets
   use mod_output_postprocess, only: output_postprocess_check_variable_element_sets
   use def_solidz

   implicit none

   integer(ip)             :: ipoin, ivoig, ibopo, idumm, i, j, ielem
   real(rp)                :: dummr, nx, ny, vdumm(3, 3), G_ij(3, 3), G_val(3)
   real(rp)                :: E_ij(3, 3), B_ij(3, 3), B_val(3), iden(3, 3)

   !-----------------------------------------------------------------
   !
   ! Initializations
   !
   !-----------------------------------------------------------------

   do ipoin = 1, npoin
      caust_sld(1:nvoig_sld, ipoin) = 0.0_rp
      green_sld(1:nvoig_sld, ipoin) = 0.0_rp
   end do

   if (output_postprocess_check_variable_postprocess(VARIABLE_NAME='LNEPS') .or. &
       output_postprocess_check_variable_witness(VARIABLE_NAME='LEPSI', NCHAR=3_ip) .or. &
       output_postprocess_check_variable_node_sets(VARIABLE_NAME='LEPSI', NCHAR=3_ip) .or. &
       output_postprocess_check_variable_element_sets(VARIABLE_NAME='LEPRT')) then
      do ipoin = 1, npoin
         lepsi_sld(1:nvoig_sld, ipoin) = 0.0_rp
      end do
   end if

   if (output_postprocess_check_variable_postprocess(VARIABLE_NAME='SEQVM') .or. &
       output_postprocess_check_variable_witness(VARIABLE_NAME='SEQVM') .or. &
       output_postprocess_check_variable_node_sets(VARIABLE_NAME='SEQVM')) then
      do ipoin = 1, npoin
         seqvm_sld(ipoin) = 0.0_rp
      end do
   end if

   if (output_postprocess_check_variable_element_sets(VARIABLE_NAME='EBFIL')) then
      do ielem = 1, nelem
         ebfil_sld(ielem) = 0.0_rp
      end do
   end if

   if (output_postprocess_check_variable_element_sets(VARIABLE_NAME='SBFIL')) then
      do ielem = 1, nelem
         sbfil_sld(ielem) = 0.0_rp
      end do
   end if

   !-----------------------------------------------------------------
   !
   ! Calculate stress/strain tensors
   !
   !-----------------------------------------------------------------

   call sld_elmope(3_ip)

   !-----------------------------------------------------------------
   !
   ! Projections to nodes
   !
   !-----------------------------------------------------------------
   !
   ! Cauchy Stress / Green strain
   !
   if (INOTEMPTY) then
      call rhsmod(nvoig_sld, green_sld)
      call rhsmod(nvoig_sld, caust_sld)
   end if
   do ipoin = 1, npoin
      dummr = 1.0_rp/vmass(ipoin)
      do ivoig = 1, nvoig_sld
         green_sld(ivoig, ipoin) = dummr*green_sld(ivoig, ipoin)
         caust_sld(ivoig, ipoin) = dummr*caust_sld(ivoig, ipoin)
      end do
   end do
   !
   ! Logarithmic strain
   !
   if (output_postprocess_check_variable_postprocess(VARIABLE_NAME='LNEPS') .or. &
       output_postprocess_check_variable_witness(VARIABLE_NAME='LEPSI', NCHAR=3_ip) .or. &
       output_postprocess_check_variable_node_sets(VARIABLE_NAME='LEPSI', NCHAR=3_ip) .or. &
       output_postprocess_check_variable_element_sets(VARIABLE_NAME='LEPRT')) then
      if (INOTEMPTY) then
         call rhsmod(nvoig_sld, lepsi_sld)
      end if
      do ipoin = 1, npoin
         dummr = 1.0_rp/vmass(ipoin)
         do ivoig = 1, nvoig_sld
            lepsi_sld(ivoig, ipoin) = dummr*lepsi_sld(ivoig, ipoin)
         end do
      end do
   end if
   !
   ! Von Mises
   !
   if (output_postprocess_check_variable_postprocess(VARIABLE_NAME='SEQVM') .or. &
       output_postprocess_check_variable_witness(VARIABLE_NAME='SEQVM') .or. &
       output_postprocess_check_variable_node_sets(VARIABLE_NAME='SEQVM')) then
      if (ndime == 2) then
         do ipoin = 1, npoin
            !General Plane Stress assumption)
            seqvm_sld(ipoin) = sqrt(caust_sld(1, ipoin)**2 + caust_sld(2, ipoin)**2 &
                                    - caust_sld(1, ipoin)*caust_sld(2, ipoin) + 3.0_rp*caust_sld(3, ipoin)**2)
         end do
      else if (ndime == 3) then
         do ipoin = 1, npoin
            seqvm_sld(ipoin) = &
               (sqrt(0.5_rp*((caust_sld(1, ipoin) - caust_sld(2, ipoin))**2 + &
                             (caust_sld(2, ipoin) - caust_sld(3, ipoin))**2 + &
                             (caust_sld(3, ipoin) - caust_sld(1, ipoin))**2 + &
                             6.0_rp*(caust_sld(4, ipoin)**2 + caust_sld(5, ipoin)**2 + caust_sld(6, ipoin)**2))))
         end do
      else
         call runend('VON MISES NOT PROGRAMMED FOR OTHER NDIME DIFFERENT THAN 2D OR 2D')
      end if
   end if
   !
   ! ELEPS, NSIGN
   !
   if (output_postprocess_check_variable_postprocess(VARIABLE_NAME='ELEPS') .or. &
       output_postprocess_check_variable_postprocess(VARIABLE_NAME='NSIGN')) then

      if (kfl_plast_sld == 1_ip) call rhsmod(nvoig_sld, epsee_sld)
      do ipoin = 1, npoin
         dummr = 1.0_rp/vmass(ipoin)
         do ivoig = 1, nvoig_sld
            if (kfl_plast_sld == 1_ip) epsee_sld(ivoig, ipoin) = dummr*epsee_sld(ivoig, ipoin)
         end do
         ibopo = lpoty(ipoin)
         caunn_sld(ipoin) = 0.0_rp
         if (ibopo > 0) then
            if (ndime == 2) then
               nx = exnor(1, 1, ibopo)
               ny = exnor(2, 1, ibopo)
               caunn_sld(ipoin) = &
                  nx*caust_sld(1, ipoin)*nx + &
                  ny*caust_sld(2, ipoin)*ny + &
                  2.0_rp*nx*caust_sld(3, ipoin)*ny
            else if (ndime == 3) then

            end if
         end if
      end do
   end if
   !
   ! Eprin and sigei
   !
   if (output_postprocess_check_variable_postprocess(VARIABLE_NAME='EPRIN') .or. &
       output_postprocess_check_variable_node_sets(VARIABLE_NAME='EPRIN') .or. &
       output_postprocess_check_variable_postprocess(VARIABLE_NAME='SIGEI')) then

      if (ndime == 3_ip) then

         iden(:, :) = 0.0_rp
         do i = 1, 3
            iden(i, i) = 1.0_rp
         end do

         do ipoin = 1, npoin
            if (kfl_rotei_sld == 1) call sld_troloc(0_ip, G_ij)     !compute roloc_sld tensor
            !G_ij(1,1)= caust_sld(1,ipoin)
            !G_ij(2,2)= caust_sld(2,ipoin)
            !G_ij(3,3)= caust_sld(3,ipoin)
            !G_ij(1,2)= caust_sld(4,ipoin)
            !G_ij(1,3)= caust_sld(5,ipoin)
            !G_ij(2,3)= caust_sld(6,ipoin)
            !G_ij(2,1)= G_ij(1,2)
            !G_ij(3,1)= G_ij(1,3)
            !G_ij(3,2)= G_ij(2,3)
            do i = 1, ndime
               do j = 1, ndime
                  ivoig = nvgij_inv_sld(i, j)
                  G_ij(i, j) = caust_sld(ivoig, ipoin)
                  E_ij(i, j) = green_sld(ivoig, ipoin)
               end do
            end do

            if (kfl_rotei_sld == 1) call sld_troloc(ipoin, G_ij)     !rotate and correct tensor
            if (kfl_rotei_sld == 1) call sld_troloc(ipoin, E_ij)     !rotate and correct tensor

            ! eigenvalues sorted in ascending order
            call spcdec(G_ij, G_val, vdumm, idumm, 0_ip, 'sld_calculate_pp_tensors, COMPUTE CAUCHY EIGENVALUES FOR SIGEI')

            G_val(1) = abs(G_val(1))
            G_val(3) = abs(G_val(3))
            if (G_val(3) .gt. G_val(1)) G_val(1) = G_val(3)

            sigei_sld(ipoin) = G_val(1)

            !
            ! Principal stretches
            !
            B_ij = transpose(2.0_ip*E_ij + iden)

            ! eigenvalues sorted in ascending order
            call spcdec(B_ij, B_val, vdumm, idumm, 0_ip, 'sld_calculate_pp_tensors, COMPUTE GREEN-LAGRANGE EIGENVALUES FOR EPRIN')

            B_val(1) = abs(B_val(1))
            B_val(3) = abs(B_val(3))
            if (B_val(3) .gt. B_val(1)) B_val(1) = B_val(3)

            eprin_sld(ipoin) = B_val(1)
         end do

      end if

   end if

end subroutine sld_calculate_pp_tensors
