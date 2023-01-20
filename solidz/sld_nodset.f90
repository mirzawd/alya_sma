!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_nodset.f90
!> @author  Gerard Guillamet
!> @date    July, 2018
!>          - Subroutine written
!> @brief   Compute and write results on node sets
!>
!> @details
!>
!>          \verbatim
!>          Node Sets:
!>          ----------
!>          1-3.   DISPL = Displacement components
!>          4-6.   VELOC = Velocity components
!>          7-9.   ACCEL = Acceleration components
!>          10-12. FRXIX = Reaction forces components
!>          13-18. SIGMA = Stresses (Cauchy)
!>          19-24. EPSIL = Strains (Green)
!>          25-30. LNEPS = Logarithmic Strains
!>          31-33. COORD = Nodal coordinate components
!>          34.    SEQVM = Von Mises Stress
!>          35-37. FEXTE = External force
!>          38-40. FCONT = Contact force
!>          41.    EPRIN = Principal stetches
!>          42-44. BVESS = Dirichlet displacement
!>          45-47. SBVNA = Solver Neumann value
!>          48.    TEMPE = Temperature
!>          \endverbatim
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_nodset(insec,inset)

  use def_kintyp, only : ip, rp
  use def_master, only : postp, vnset
  use def_master, only : displ, therm
  use def_master, only : TIME_N, ITER_K
  use def_master, only : solve_sol
  use def_domain, only : ndime
  use def_domain, only : coord
  use def_solidz, only : kfl_conta_sld, kfl_rigid_sld
  use def_solidz, only : veloc_sld, accel_sld
  use def_solidz, only : frxid_sld, fexte_sld, fcont_sld
  use def_solidz, only : caust_sld, green_sld, lepsi_sld
  use def_solidz, only : seqvm_sld, eprin_sld
  use def_solidz, only : bvess_sld

  implicit none

  integer(ip), intent(in) :: insec     !< Node number
  integer(ip), intent(in) :: inset     !< Set code
  integer(ip)             :: ipoin, idofn, jdofn
  real(rp)                :: rfaux(ndime), faux(ndime), rfcon(ndime)
  real(rp)                :: sbvna(ndime)

  !----------------------------------------------------------------------
  !
  ! Node sets
  !
  !----------------------------------------------------------------------

  ipoin = insec
  if (ipoin/=0) then
     idofn = (ipoin-1) * ndime + 1
     jdofn = idofn + ndime - 1
     !
     ! DISPX, DISPY, DISPZ
     !
     if (postp(1) % npp_setsn(1) /= 0) vnset(1,inset)= displ(1,ipoin,TIME_N)
     if (postp(1) % npp_setsn(2) /= 0) vnset(2,inset)= displ(2,ipoin,TIME_N)
     if (postp(1) % npp_setsn(3) /= 0 .and. ndime == 3_ip) vnset(3,inset)= displ(3,ipoin,TIME_N)
     !
     ! VELOX, VELOY, VELOZ
     !
     if (postp(1) % npp_setsn(4) /= 0) vnset(4,inset)= veloc_sld(1,ipoin,TIME_N)
     if (postp(1) % npp_setsn(5) /= 0) vnset(5,inset)= veloc_sld(2,ipoin,TIME_N)
     if (postp(1) % npp_setsn(6) /= 0 .and. ndime == 3_ip) vnset(6,inset)= veloc_sld(3,ipoin,TIME_N)
     !
     ! ACCEX, ACCEY, ACCEZ
     !
     if (postp(1) % npp_setsn(7) /= 0) vnset(7,inset)= accel_sld(1,ipoin,TIME_N)
     if (postp(1) % npp_setsn(8) /= 0) vnset(8,inset)= accel_sld(2,ipoin,TIME_N)
     if (postp(1) % npp_setsn(9) /= 0 .and. ndime == 3_ip) vnset(9,inset)= accel_sld(3,ipoin,TIME_N)
     !
     ! FRXIX, FRXIY, FRXIZ
     !
     rfaux(:) = 0.0_rp
     if (kfl_rigid_sld == 0) rfaux(1:ndime) =  -frxid_sld(idofn:jdofn) ! opposite sign
     if (postp(1) % npp_setsn(10) /= 0) vnset(10,inset) = rfaux(1)
     if (postp(1) % npp_setsn(11) /= 0) vnset(11,inset) = rfaux(2)
     if (postp(1) % npp_setsn(12) /= 0 .and. ndime == 3_ip) vnset(12,inset) = rfaux(3)
     !
     ! SIGXX, SIGYY, SIGZZ, SIGYZ, SIGXZ, SIGXY
     !
     if (postp(1) % npp_setsn(13) /= 0) vnset(13,inset) = caust_sld(1,ipoin)                    !sxx
     if (postp(1) % npp_setsn(14) /= 0) vnset(14,inset) = caust_sld(2,ipoin)                    !syy
     if (postp(1) % npp_setsn(15) /= 0 .and. ndime == 3_ip) vnset(15,inset) = caust_sld(3,ipoin)!szz
     if (postp(1) % npp_setsn(16) /= 0 .and. ndime == 3_ip) vnset(16,inset) = caust_sld(4,ipoin)!syz
     if (postp(1) % npp_setsn(17) /= 0 .and. ndime == 3_ip) vnset(17,inset) = caust_sld(5,ipoin)!sxz
     if (postp(1) % npp_setsn(18) /= 0 ) then                                                   !sxy
        if( ndime == 2_ip ) then
           vnset(18,inset) = caust_sld(3,ipoin)
        else
           vnset(18,inset) = caust_sld(6,ipoin)
        end if
     end if
     !
     ! EPSXX, EPSYY, EPSZZ, EPSYZ, EPSXZ, EPSXY
     !
     if (postp(1) % npp_setsn(19) /= 0) vnset(19,inset) = green_sld(1,ipoin)                    !exx
     if (postp(1) % npp_setsn(20) /= 0) vnset(20,inset) = green_sld(2,ipoin)                    !eyy
     if (postp(1) % npp_setsn(21) /= 0 .and. ndime == 3_ip) vnset(21,inset) = green_sld(3,ipoin)!ezz
     if (postp(1) % npp_setsn(22) /= 0 .and. ndime == 3_ip) vnset(22,inset) = green_sld(4,ipoin)!eyz
     if (postp(1) % npp_setsn(23) /= 0 .and. ndime == 3_ip) vnset(23,inset) = green_sld(5,ipoin)!exz
     if (postp(1) % npp_setsn(24) /= 0 ) then                                                   !exy
        if( ndime == 2_ip ) then
           vnset(24,inset) = green_sld(3,ipoin)
        else
           vnset(24,inset) = green_sld(6,ipoin)
        end if
     end if
     !
     ! LEPXX, LEPYY, LEPZZ, LEPYZ, LEPXZ, LEPXY
     !
     if (postp(1) % npp_setsn(25) /= 0) vnset(25,inset) = lepsi_sld(1,ipoin)                    !exx
     if (postp(1) % npp_setsn(26) /= 0) vnset(26,inset) = lepsi_sld(2,ipoin)                    !eyy
     if (postp(1) % npp_setsn(27) /= 0 .and. ndime == 3_ip) vnset(27,inset) = lepsi_sld(3,ipoin)!ezz
     if (postp(1) % npp_setsn(28) /= 0 .and. ndime == 3_ip) vnset(28,inset) = lepsi_sld(4,ipoin)!eyz
     if (postp(1) % npp_setsn(29) /= 0 .and. ndime == 3_ip) vnset(29,inset) = lepsi_sld(5,ipoin)!exz
     if (postp(1) % npp_setsn(30) /= 0 ) then                                                   !exy
        if( ndime == 2_ip ) then
           vnset(30,inset) = lepsi_sld(3,ipoin)
        else
           vnset(30,inset) = lepsi_sld(6,ipoin)
        end if
     end if
     !
     ! COORX, COORY, COORZ
     !
     if (postp(1) % npp_setsn(31) /= 0) vnset(31,inset) = coord(1,ipoin)
     if (postp(1) % npp_setsn(32) /= 0) vnset(32,inset) = coord(2,ipoin)
     if (postp(1) % npp_setsn(33) /= 0 .and. ndime == 3_ip) vnset(33,inset) = coord(3,ipoin)
     !
     ! SEQVM
     !
     if (postp(1) % npp_setsn(34) /= 0) vnset(34,inset) = seqvm_sld(ipoin)
     !
     ! FEXTX, FEXTY, FEXTZ
     !
     faux(:) = 0.0_rp
     if (kfl_rigid_sld == 0) faux(1:ndime) = fexte_sld(idofn:jdofn)
     if (postp(1) % npp_setsn(35) /= 0) vnset(35,inset) = faux(1)
     if (postp(1) % npp_setsn(36) /= 0) vnset(36,inset) = faux(2)
     if (postp(1) % npp_setsn(37) /= 0 .and. ndime == 3_ip) vnset(37,inset) = faux(3)
     !
     ! FCONX, FCONY, FCONZ
     !
     rfcon(:) = 0.0_rp
     if (kfl_conta_sld /= 0 .and. kfl_rigid_sld == 0) rfcon(1:ndime) = fcont_sld(1:ndime,ipoin)
     if (postp(1) % npp_setsn(38) /= 0) vnset(38,inset) = rfcon(1)
     if (postp(1) % npp_setsn(39) /= 0) vnset(39,inset) = rfcon(2)
     if (postp(1) % npp_setsn(40) /= 0 .and. ndime == 3_ip) vnset(40,inset) = rfcon(3)
     !
     ! EPRIN: principal sretches
     !
     if (postp(1) % npp_setsn(41) /= 0) vnset(41,inset) = eprin_sld(ipoin)
     !
     ! BVESX, BVESY, BVESZ
     !
     if (postp(1) % npp_setsn(42) /= 0) vnset(42,inset) = bvess_sld(1,ipoin,ITER_K)
     if (postp(1) % npp_setsn(43) /= 0) vnset(43,inset) = bvess_sld(2,ipoin,ITER_K)
     if (postp(1) % npp_setsn(44) /= 0 .and. ndime == 3_ip) vnset(44,inset) = bvess_sld(3,ipoin,ITER_K)
     !
     ! SBVNAT
     !
     sbvna(:) = 0.0_rp
     if (solve_sol(1)%kfl_bvnat == 1 .and. kfl_rigid_sld == 0) sbvna(1:ndime) = solve_sol(1)%bvnat(1:ndime,ipoin)
     if (postp(1) % npp_setsn(45) /= 0) vnset(45,inset) = sbvna(1)
     if (postp(1) % npp_setsn(46) /= 0) vnset(46,inset) = sbvna(2)
     if (postp(1) % npp_setsn(47) /= 0 .and. ndime == 3_ip) vnset(47,inset) = sbvna(3)
     !
     ! TEMPE
     !
     if (postp(1) % npp_setsn(48) /= 0) vnset(48,inset) = therm(ipoin,ITER_K)
     
  end if

end subroutine sld_nodset
