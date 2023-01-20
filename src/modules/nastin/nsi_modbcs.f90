!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



  !------------------------------------------------------------------------
  !> @addtogroup Nastin 
  !> @{
  !> @file    nsi_modbcs.f90
  !> @author  Herbert Owen
  !> @brief   Modify the boundary conditions manually
  !> @details Modify the boundary conditions manually
  !> @} 
  !------------------------------------------------------------------------
subroutine nsi_modbcs()
  use def_elmtyp
  use def_master
  use def_kermod
  use def_domain
  use def_nastin
  implicit none
  integer(ip) :: ipoin,ibopo

  if( IMASTER ) return
  if( kfl_modfi_nsi == 2 ) then
     do ipoin = 1,npoin
        if (abs(coord(1,ipoin)-3.0_rp)<1.0d-6) then
           ibopo = lpoty(ipoin)
           kfl_fixpr_nsi(1,ipoin) = 1
        end if
        ! modification to free the trailing edge of a specific naca 
        if ( (abs(coord(1,ipoin)-0.2022230_rp)<1.0d-4) .and. (abs(coord(2,ipoin)+0.0176851_rp)<1.0d-4) ) then
           ibopo = lpoty(ipoin)
           kfl_fixrs_nsi(ipoin) = 0   ! no skew system
           kfl_fixno_nsi(1,ipoin) = 0
           kfl_fixno_nsi(2,ipoin) = 0
           if (ndime==3) kfl_fixno_nsi(3,ipoin) = 1
        end if
     end do
  else if( kfl_modfi_nsi == 3 ) then
     do ipoin = 1,npoin
        if (abs(coord(1,ipoin)-2.25_rp)<1.0d-6) then
           ibopo = lpoty(ipoin)
           kfl_fixpr_nsi(1,ipoin) = 1
        end if
     end do
  else if( kfl_modfi_nsi == 4 ) then
     do ipoin = 1,npoin
        if (abs(coord(1,ipoin)-11.0837_rp)<1.0d-6) then
           ibopo = lpoty(ipoin)
           kfl_fixpr_nsi(1,ipoin) = 1
        end if
        ! modification to free the trailing edge of a specific naca 
        if ( (abs(coord(1,ipoin)-0.4980720_rp)<1.0d-6) .and. (abs(coord(2,ipoin)+0.0435490_rp)<1.0d-6) ) then
           ibopo = lpoty(ipoin)
           kfl_fixrs_nsi(ipoin) = 0   ! no skew system
           kfl_fixno_nsi(1,ipoin) = 0
           kfl_fixno_nsi(2,ipoin) = 0
           if (ndime==3) kfl_fixno_nsi(3,ipoin) = 1
        end if
     end do
  else if( kfl_modfi_nsi == 5 ) then
     do ipoin = 1,npoin
        if (abs(coord(1,ipoin)-5.0_rp)<1.0d-4) then
           ibopo = lpoty(ipoin)
           kfl_fixpr_nsi(1,ipoin) = 1
        end if
        ! modification to free the trailing edge of a specific naca 
        if ( (abs(coord(1,ipoin)-0.101109_rp)<1.0d-4) .and. (abs(coord(2,ipoin)+0.0088404_rp)<1.0d-4) ) then
           ibopo = lpoty(ipoin)
           kfl_fixrs_nsi(ipoin) = 0   ! no skew system
           kfl_fixno_nsi(1,ipoin) = 0
           kfl_fixno_nsi(2,ipoin) = 0
           if (ndime==3) kfl_fixno_nsi(3,ipoin) = 1
        end if
     end do
  else if( kfl_modfi_nsi == 6 ) then  ! same as 5 but also fixpr in upper face , beware it will also affect the top inlet node
     do ipoin = 1,npoin
        if ((abs(coord(1,ipoin)-5.0_rp)<1.0d-4).or.(abs(coord(2,ipoin)-0.46_rp)<1.0d-4)) then
           ibopo = lpoty(ipoin)
           kfl_fixpr_nsi(1,ipoin) = 1
        end if
        ! modification to free the trailing edge of a specific naca 
        if ( (abs(coord(1,ipoin)-0.101109_rp)<1.0d-4) .and. (abs(coord(2,ipoin)+0.0088404_rp)<1.0d-4) ) then
           ibopo = lpoty(ipoin)
           kfl_fixrs_nsi(ipoin) = 0   ! no skew system
           kfl_fixno_nsi(1,ipoin) = 0
           kfl_fixno_nsi(2,ipoin) = 0
           if (ndime==3) kfl_fixno_nsi(3,ipoin) = 1
        end if
     end do
  end if

end subroutine nsi_modbcs
