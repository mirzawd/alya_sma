!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_outset.f90
!> @author  Solidz Team
!> @date    April, 2010
!>          - Subroutine creation
!> @brief   Compute and write results on sets
!> @details Sets can be defined at element, boundary and node levels.
!> @note    The keyflag kfl_foten_sld is only activated when the
!>          user requires the calculation of tensors related to
!>          strains and stress.
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_outset()

   use def_kintyp, only: ip, rp
   use def_master, only: INOTMASTER, INOTSLAVE, ittim, mitim, cutim, timef
   use def_master, only: postp, vbset, veset
   use def_domain, only: nnset, neset, nbset
   use def_domain, only: lnsec, lbsec, lesec
   use def_solidz, only: kfl_foten_sld
   use mod_output_postprocess, only: output_postprocess_node_sets_parall
   use mod_output_postprocess, only: output_postprocess_boundary_sets_parall
   use mod_output_postprocess, only: output_postprocess_element_sets_parall

   implicit none

   integer(ip) :: ieset, ibset, inset, dummi, i
   ! VELOC, EBFIL, SBFIL, EPSUL, EPSUR, EPSUC
   integer(ip), parameter :: sets_to_normalize(6) = [15, 18, 19, 20, 21, 22]

   !----------------------------------------------------------------------
   !
   ! Element sets
   !
   !----------------------------------------------------------------------

   if (maxval(postp(1)%npp_setse) > 0) then
      if ((mod(ittim, postp(1)%npp_stepelset) == 0) .or. (ittim .eq. mitim) .or. (cutim >= timef)) then
         if (INOTMASTER) then
            !
            ! Call fotens if necessary
            if (any(postp(1)%npp_setse(1:14) /= 0_ip) .or. &
                any(postp(1)%npp_setse(18:22) /= 0_ip)) then
               kfl_foten_sld = kfl_foten_sld + 1_ip
               if (kfl_foten_sld == 1_ip) call sld_calculate_pp_tensors
            end if

            do ieset = 1, neset
               call sld_elmset(lesec(ieset), ieset)
            end do
         end if
         !
         ! Parall
         !
         call output_postprocess_element_sets_parall()
         !
         ! Averaged/Norm of some variables
         !
         if (INOTSLAVE) then
            do ieset = 1, neset
               if (veset(postp(1)%nvaes + 1, ieset) > 0.0_rp) then
                  do i = 1, size(sets_to_normalize)
                     if (postp(1)%npp_setse(sets_to_normalize(i)) /= 0) then
                        veset(sets_to_normalize(i), ieset) = &
                           veset(sets_to_normalize(i), ieset)/ &
                           veset(postp(1)%nvaes + 1, ieset)
                     end if
                  end do
               end if
            end do
         end if
      end if
   end if

   !----------------------------------------------------------------------
   !
   ! Boundary sets
   !
   !----------------------------------------------------------------------

   if (maxval(postp(1)%npp_setsb) > 0) then
      if ((mod(ittim, postp(1)%npp_stepboset) == 0) .or. (ittim .eq. mitim) .or. (cutim >= timef)) then
         if (INOTMASTER) then
            do ibset = 1, nbset
               call sld_bouset(lbsec(ibset), ibset)
            end do
         end if
         !
         ! Parall
         !
         call output_postprocess_boundary_sets_parall()
         !
         ! Averaged/Norm of some variables
         !
         if (INOTSLAVE) then
            do ibset = 1, nbset
               ! DIBOX
               if (postp(1)%npp_setsb(5) /= 0 .and. vbset(postp(1)%nvabs + 1, ibset) > 0) then
                  vbset(5, ibset) = vbset(5, ibset)/vbset(postp(1)%nvabs + 1, ibset)
               end if
               ! DIBOY
               if (postp(1)%npp_setsb(6) /= 0 .and. vbset(postp(1)%nvabs + 1, ibset) > 0) then
                  vbset(6, ibset) = vbset(6, ibset)/vbset(postp(1)%nvabs + 1, ibset)
               end if
               ! DIBOZ
               if (postp(1)%npp_setsb(7) /= 0 .and. vbset(postp(1)%nvabs + 1, ibset) > 0) then
                  vbset(7, ibset) = vbset(7, ibset)/vbset(postp(1)%nvabs + 1, ibset)
               end if
               ! DIBOU
               if (postp(1)%npp_setsb(8) /= 0 .and. vbset(postp(1)%nvabs + 1, ibset) > 0) then
                  vbset(8, ibset) = vbset(8, ibset)/vbset(postp(1)%nvabs + 1, ibset)
               end if
               ! FRBOU
               if (postp(1)%npp_setsb(12) /= 0_ip) then
                  vbset(12, ibset) = sqrt(sum(vbset(9:11, ibset)**2))
               end if
               ! FCONT
               if (postp(1)%npp_setsb(16) /= 0) then
                  vbset(16, ibset) = sqrt(sum(vbset(13:15, ibset)**2))
               end if
               ! PRESS
               if (postp(1)%npp_setsb(17) /= 0 .and. vbset(postp(1)%nvabs + 1, ibset) > 0) then
                  vbset(17, ibset) = vbset(17, ibset)/vbset(postp(1)%nvabs + 1, ibset)
               end if
               ! TRACTION
               if (postp(1)%npp_setsb(18) /= 0 .and. vbset(postp(1)%nvabs + 1, ibset) > 0) then
                  vbset(18, ibset) = vbset(18, ibset)/vbset(postp(1)%nvabs + 1, ibset)
               end if
            end do
         end if
      end if
   end if

   !----------------------------------------------------------------------
   !
   ! Node sets
   !
   !----------------------------------------------------------------------

   if (maxval(postp(1)%npp_setsn) > 0) then
      if ((mod(ittim, postp(1)%npp_stepnoset) == 0) .or. (ittim .eq. mitim) .or. (cutim >= timef)) then
         if (INOTMASTER) then
            !
            ! Call fotens if any of these variables are called:
            ! SIGMA, EPSIL, LNEPS, SEQVM
            if (any(postp(1)%npp_setsn(13:30) /= 0_ip) .or. postp(1)%npp_setsn(34) /= 0_ip) then
               kfl_foten_sld = kfl_foten_sld + 1_ip
               if (kfl_foten_sld == 1_ip) call sld_calculate_pp_tensors
            end if

            do inset = 1, nnset
               call sld_nodset(lnsec(inset), inset)
            end do

         end if
         !
         ! Parall
         !
         call output_postprocess_node_sets_parall()
      end if
   end if

end subroutine sld_outset
