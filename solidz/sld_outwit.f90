!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_outwit.f90
!> @author  Solidz Team
!> @date    September, 2011
!>          - Subroutine creation
!> @brief   Compute and write results on witness points
!>
!> @details
!>
!>          \verbatim
!>          Witness point variables:
!>          --------
!>          1-3.   DISPL = Displacement components
!>          4-6.   VELOC = Velocity components
!>          7-9.   ACCEL = Acceleration components
!>          10-12. FRXIX = Reaction forces components
!>          13-18. SIGMA = Stresses (Cauchy)
!>          19-24. EPSIL = Strains (Green)
!>          25-30. LNEPS = Logarithmic Strains
!>          31-33. COORD = Nodal coordinate components
!>          34.    SEQVM = Von Mises Stress
!>          35.    SIGFI =
!>          36.    MICRO =
!>          37.    EBFIL = Strain projected on longitudinal biofiber direction
!>          38.    SBFIL = Stress projected on longitudinal biofiber direction
!>          39.    TEMPE = Temperature
!>          \endverbatim
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_outwit()

  use def_kintyp,      only : ip, rp
  use def_master,      only : postp, witne, INOTMASTER, ittim, leinv_loc
  use def_master,      only : displ, therm
  use def_master,      only : ITER_K
  use def_kermod,      only : nwitn, lewit, shwit
  use def_domain,      only : ndime, nnode, lnods, ltype, coord
  use def_solidz,      only : veloc_sld, accel_sld
  use def_solidz,      only : frxid_sld
  use def_solidz,      only : kfl_foten_sld, caust_sld, green_sld, lepsi_sld
  use def_solidz,      only : seqvm_sld
  use def_solidz,      only : fibde_sld
  use mod_biofibers,   only : biofib_point_nodal_fibers
  use mod_sld_fe2
  use mod_maths_basic, only :  maths_MULT_VxMxV
  use mod_sld_stress_model_comput, only : SM_tensor_to_voigt_second
  use mod_output_postprocess, only : output_postprocess_check_witness_output

  implicit none

  integer(ip)          :: iwitn,pnode,pelty
  integer(ip)          :: ielem,inode,ipoin,idofn
  real(rp)             :: fb(3),st(6),stfib,fmod,rfaux
  real(rp)             :: matrix(ndime, ndime)
  real(rp),   pointer  :: vecptr(:,:)

  !
  ! Check: write witness or not?
  !
  if( output_postprocess_check_witness_output() ) then
     

     if( INOTMASTER ) then
        !
        ! Call sld_calculate_pp_tensors if any of these variables are called for post-process:
        ! SIGMA, EPSIL, LNEPS, SEQVM, SIGFI, EBFIL, SBFIL
        if( any(postp(1) % npp_witne(13:30) /= 0_ip) .or. &
            any(postp(1) % npp_witne(34:35) /= 0_ip) .or. &
            any(postp(1) % npp_witne(37:38) /= 0_ip) ) then
           kfl_foten_sld = kfl_foten_sld + 1_ip
           if( kfl_foten_sld == 1_ip ) call sld_calculate_pp_tensors
        end if
        !
        ! Fidbes
        !
        ! Extract fibers (tensor projections)
        if (any(postp(1)%npp_witne([35, 37, 38]) == 1_ip)) then
           call biofib_point_nodal_fibers(fibde_sld, 'LONGITUDINAL', 'CURRENT')
        end if

        do iwitn = 1, nwitn
           ielem = lewit(iwitn)
           if( ielem > 0 ) then

              pelty = ltype(ielem)
              pnode = nnode(pelty)
              !
              ! DISPX
              !
              if( postp(1) % npp_witne(1) == 1_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(1,iwitn) = witne(1,iwitn) + shwit(inode,iwitn) *  displ(1,ipoin,1) * postp(1) % witne_dt(1)
                 end do
              end if
              !
              ! DISPY
              !
              if( postp(1) % npp_witne(2) == 1_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(2,iwitn) = witne(2,iwitn) + shwit(inode,iwitn) *  displ(2,ipoin,1) * postp(1) % witne_dt(2)
                 end do
              end if
              !
              ! DISPZ
              !
              if( postp(1) % npp_witne(3) == 1_ip .and. ndime == 3_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(3,iwitn) = witne(3,iwitn) + shwit(inode,iwitn) *  displ(ndime,ipoin,1) * postp(1) % witne_dt(3)
                 end do
              end if
              !
              ! VELOX
              !
              if( postp(1) % npp_witne(4) == 1_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(4,iwitn) = witne(4,iwitn) + shwit(inode,iwitn)*veloc_sld(1,ipoin,1) * postp(1) % witne_dt(4)
                 end do
              end if
              !
              ! VELOY
              !
              if( postp(1) % npp_witne(5) == 1_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(5,iwitn) = witne(5,iwitn) + shwit(inode,iwitn)*veloc_sld(2,ipoin,1) * postp(1) % witne_dt(5) 
                 end do
              end if
              !
              ! VELOZ
              !
              if( postp(1) % npp_witne(6) == 1_ip .and. ndime == 3_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(6,iwitn) = witne(6,iwitn) + shwit(inode,iwitn)*veloc_sld(3,ipoin,1) * postp(1) % witne_dt(6) 
                 end do
              end if
              !
              ! ACCEX
              !
              if( postp(1) % npp_witne(7) == 1_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(7,iwitn) = witne(7,iwitn) + shwit(inode,iwitn)*accel_sld(1,ipoin,1) * postp(1) % witne_dt(7) 
                 end do
              end if
              !
              ! ACCEY
              !
              if( postp(1) % npp_witne(8) == 1_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(8,iwitn) = witne(8,iwitn) + shwit(inode,iwitn)*accel_sld(2,ipoin,1) * postp(1) % witne_dt(8) 
                 end do
              end if
              !
              ! ACCEZ
              !
              if( postp(1) % npp_witne(9) == 1_ip .and. ndime == 3_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(9,iwitn) = witne(9,iwitn) + shwit(inode,iwitn)*accel_sld(3,ipoin,1) * postp(1) % witne_dt(9) 
                 end do
              end if
              !
              ! FRXIX
              !
              if( postp(1) % npp_witne(10) == 1_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    idofn = (ipoin-1)*ndime + 1
                    rfaux = -frxid_sld(idofn) ! opposite sign
                    witne(10,iwitn) = witne(10,iwitn) + shwit(inode,iwitn)*rfaux * postp(1) % witne_dt(10) 
                 end do
              end if
              !
              ! FRXIY
              !
              if( postp(1) % npp_witne(11) == 1_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    idofn = (ipoin-1)*ndime + 2
                    rfaux = -frxid_sld(idofn) ! opposite sign
                    witne(11,iwitn) = witne(11,iwitn) + shwit(inode,iwitn)*rfaux * postp(1) % witne_dt(11) 
                 end do
              end if
              !
              ! FRXIZ
              !
              if( postp(1) % npp_witne(12) == 1_ip .and. ndime == 3_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    idofn = (ipoin-1)*ndime + 3
                    rfaux = -frxid_sld(idofn) ! opposite sign
                    witne(12,iwitn) = witne(12,iwitn) + shwit(inode,iwitn)*rfaux * postp(1) % witne_dt(12) 
                 end do
              end if
              !
              ! SIGXX
              !
              if( postp(1) % npp_witne(13) == 1_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(13,iwitn) = witne(13,iwitn) + shwit(inode,iwitn)*caust_sld(1,ipoin) * postp(1) % witne_dt(13) 
                 end do
              end if
              !
              ! SIGYY
              !
              if( postp(1) % npp_witne(14) == 1_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(14,iwitn) = witne(14,iwitn) + shwit(inode,iwitn)*caust_sld(2,ipoin) * postp(1) % witne_dt(14) 
                 end do
              end if
              !
              ! SIGZZ
              !
              if( postp(1) % npp_witne(15) == 1_ip .and. ndime == 3_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(15,iwitn) = witne(15,iwitn) + shwit(inode,iwitn)*caust_sld(3,ipoin) * postp(1) % witne_dt(15) 
                 end do
              end if
              !
              ! SIGYZ
              !
              if( postp(1) % npp_witne(16) == 1_ip .and. ndime == 3_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(16,iwitn) = witne(16,iwitn) + shwit(inode,iwitn)*caust_sld(4,ipoin) * postp(1) % witne_dt(16) 
                 end do
              end if
              !
              ! SIGXZ
              !
              if( postp(1) % npp_witne(17) == 1_ip .and. ndime == 3_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(17,iwitn) = witne(17,iwitn) + shwit(inode,iwitn)*caust_sld(5,ipoin) * postp(1) % witne_dt(17) 
                 end do
              end if
              !
              ! SIGXY
              !
              if( postp(1) % npp_witne(18) == 1_ip ) then
                 if( ndime == 2_ip ) then
                    do inode = 1,pnode
                       ipoin = lnods(inode,ielem)
                       witne(18,iwitn) = witne(18,iwitn) + shwit(inode,iwitn)*caust_sld(3,ipoin) * postp(1) % witne_dt(18) 
                    end do
                 else
                    do inode = 1,pnode
                       ipoin = lnods(inode,ielem)
                       witne(18,iwitn) = witne(18,iwitn) + shwit(inode,iwitn)*caust_sld(6,ipoin) * postp(1) % witne_dt(18) 
                    end do
                 end if
              end if
              !
              ! EPSXX
              !
              if( postp(1) % npp_witne(19) == 1_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(19,iwitn) = witne(19,iwitn) + shwit(inode,iwitn)*green_sld(1,ipoin) * postp(1) % witne_dt(19) 
                 end do
              end if
              !
              ! EPSYY
              !
              if( postp(1) % npp_witne(20) == 1_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(20,iwitn) = witne(20,iwitn) + shwit(inode,iwitn)*green_sld(2,ipoin) * postp(1) % witne_dt(20) 
                 end do
              end if
              !
              ! EPSZZ
              !
              if( postp(1) % npp_witne(21) == 1_ip .and. ndime == 3_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(21,iwitn) = witne(21,iwitn) + shwit(inode,iwitn)*green_sld(3,ipoin) * postp(1) % witne_dt(21) 
                 end do
              end if
              !
              ! EPSYZ
              !
              if( postp(1) % npp_witne(22) == 1_ip .and. ndime == 3_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(22,iwitn) = witne(22,iwitn) + shwit(inode,iwitn)*green_sld(4,ipoin) * postp(1) % witne_dt(22)
                 end do
              end if
              !
              ! EPSXZ
              !
              if( postp(1) % npp_witne(23) == 1_ip .and. ndime == 3_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(23,iwitn) = witne(23,iwitn) + shwit(inode,iwitn)*green_sld(5,ipoin) * postp(1) % witne_dt(23)
                 end do
              end if
              !
              ! EPSXY
              !
              if( postp(1) % npp_witne(24) == 1_ip ) then
                 if (ndime == 2_ip) then
                    do inode = 1,pnode
                       ipoin = lnods(inode,ielem)
                       witne(24,iwitn) = witne(24,iwitn) + shwit(inode,iwitn)*green_sld(3,ipoin) * postp(1) % witne_dt(24)
                    end do
                 else
                    do inode = 1,pnode
                       ipoin = lnods(inode,ielem)
                       witne(24,iwitn) = witne(24,iwitn) + shwit(inode,iwitn)*green_sld(6,ipoin) * postp(1) % witne_dt(24)
                    end do
                 end if
              end if
              !
              ! LEPXX
              !
              if( postp(1) % npp_witne(25) == 1_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(25,iwitn) = witne(25,iwitn) + shwit(inode,iwitn)*lepsi_sld(1,ipoin) * postp(1) % witne_dt(25)
                 end do
              end if
              !
              ! LEPYY
              !
              if( postp(1) % npp_witne(26) == 1_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(26,iwitn) = witne(26,iwitn) + shwit(inode,iwitn)*lepsi_sld(2,ipoin) * postp(1) % witne_dt(26)
                 end do
              end if
              !
              ! LEPZZ
              !
              if( postp(1) % npp_witne(27) == 1_ip .and. ndime == 3_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(27,iwitn) = witne(27,iwitn) + shwit(inode,iwitn)*lepsi_sld(3,ipoin) * postp(1) % witne_dt(27)
                 end do
              end if
              !
              ! LEPYZ
              !
              if( postp(1) % npp_witne(28) == 1_ip .and. ndime == 3_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(28,iwitn) = witne(28,iwitn) + shwit(inode,iwitn)*lepsi_sld(4,ipoin) * postp(1) % witne_dt(28)
                 end do
              end if
              !
              ! LEPXZ
              !
              if( postp(1) % npp_witne(29) == 1_ip .and. ndime == 3_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(29,iwitn) = witne(29,iwitn) + shwit(inode,iwitn)*lepsi_sld(5,ipoin) * postp(1) % witne_dt(29)
                 end do
              end if
              !
              ! LEPXY
              !
              if( postp(1) % npp_witne(30) == 1_ip ) then
                 if (ndime == 2_ip) then
                    do inode = 1,pnode
                       ipoin = lnods(inode,ielem)
                       witne(30,iwitn) = witne(30,iwitn) + shwit(inode,iwitn)*lepsi_sld(3,ipoin) * postp(1) % witne_dt(30)
                    end do
                 else
                    do inode = 1,pnode
                       ipoin = lnods(inode,ielem)
                       witne(30,iwitn) = witne(30,iwitn) + shwit(inode,iwitn)*lepsi_sld(6,ipoin) * postp(1) % witne_dt(30)
                    end do
                 end if
              end if
              !
              ! COORX
              !
              if( postp(1) % npp_witne(31) == 1_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(31,iwitn) = witne(31,iwitn) + shwit(inode,iwitn)*coord(1,ipoin) * postp(1) % witne_dt(31)
                 end do
              end if
              !
              ! COORY
              !
              if( postp(1) % npp_witne(32) == 1_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(32,iwitn) = witne(32,iwitn) + shwit(inode,iwitn)*coord(2,ipoin) * postp(1) % witne_dt(32)
                 end do
              end if
              !
              ! COORZ
              !
              if( postp(1) % npp_witne(33) == 1_ip .and. ndime == 3_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(33,iwitn) = witne(33,iwitn) + shwit(inode,iwitn)*coord(3,ipoin) * postp(1) % witne_dt(33)
                 end do
              end if
              !
              ! SEQVM
              !
              if( postp(1) % npp_witne(34) == 1_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(34,iwitn) = witne(34,iwitn) + shwit(inode,iwitn)*seqvm_sld(ipoin) * postp(1) % witne_dt(34)
                 end do
              end if
              !
              ! SIGFI
              !
              if( postp(1) % npp_witne(35) == 1_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    st(1:6)= caust_sld(1:6,ipoin)
                    fb(1:3)= fibde_sld(1:3,ipoin)
                    fmod= fb(1)*fb(1)+fb(2)*fb(2)+fb(3)*fb(3)
                    stfib= 0.0_rp
                    if (fmod > 1.0e-10_rp ) then
                       stfib= (st(1)*fb(1)+st(6)*fb(2)+st(5)*fb(3))*fb(1)&
                            + (st(6)*fb(1)+st(2)*fb(2)+st(4)*fb(3))*fb(2)&
                            + (st(5)*fb(1)+st(4)*fb(2)+st(3)*fb(3))*fb(3)
                       stfib= stfib/fmod
                    end if
                    witne(35,iwitn) = witne(35,iwitn) + shwit(inode,iwitn) * stfib * postp(1) % witne_dt(35)
                 end do
              end if
              !
              ! MICRO
              !
              if( postp(1) % npp_witne(36) == 1_ip ) then
                 call fe2_write_micropp(ielem, leinv_loc(ielem), ittim)
              end if
              !
              ! EBFIL
              !
              if( postp(1) % npp_witne(37) == 1_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)

                    witne(37,iwitn) = witne(37,iwitn) + shwit(inode,iwitn) * postp(1) % witne_dt(37) * &
                                        tensor_project_normalized(green_sld(:, ipoin), fibde_sld(1:ndime, ipoin))

                 end do
              end if
              !
              ! SBFIL
              !
              if( postp(1) % npp_witne(38) == 1_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    witne(38,iwitn) = witne(38,iwitn) + shwit(inode,iwitn) * postp(1) % witne_dt(38) * &
                                        tensor_project_normalized(caust_sld(:, ipoin), fibde_sld(1:ndime, ipoin))
                 end do
              end if
              !
              ! TEMPE
              !
              if( postp(1) % npp_witne(39) == 1_ip ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)       
                    witne(39,iwitn) = witne(39,iwitn) + shwit(inode,iwitn)*therm(ipoin,ITER_K) * postp(1) % witne_dt(39)
                 end do
              end if

           end if

        end do

     end if

  end if

contains
   real(rp) function tensor_project_normalized(T, v)
      use mod_maths_basic, only: maths_normalize_vector
      use mod_sld_stress_model_comput, only: SM_voigt_to_tensor_second
      implicit none
      real(rp), dimension(6), intent(in) :: T !tensor
      real(rp), dimension(3), intent(in) :: v !vector
      real(rp), dimension(3, 3)          :: matrix
      real(rp), dimension(3)             :: v1

      call SM_voigt_to_tensor_second(3_ip, matrix, T)
      v1(:) = v(:)
      call maths_normalize_vector(3_ip, v1)
      tensor_project_normalized = maths_MULT_VxMxV(v1, 3_ip, matrix, v1, 3_ip)

   end function

end subroutine sld_outwit

