!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_roback.f90
!> @author  Mariano Vazquez
!> @date    25/04/2016
!> @brief   Rotate back boundary conditions
!> @details Rotate back boundary conditions
!> @}
!-----------------------------------------------------------------------

subroutine sld_roback

  use def_kintyp, only : ip, rp
  use def_master, only : INOTMASTER
  use def_domain, only : npoin, lpoty
  use def_solidz, only : ndofn_sld
  use def_solidz, only : kfl_fixrs_sld, kfl_fixno_sld
  use def_solidz, only : dunkn_sld
  use def_solidz, only : jacrot_du_dq_sld, jacrot_dq_du_sld

  implicit none

  integer(ip)     :: ipoin,ievat,jevat,ibopo
  integer(ip)     :: dummy_matrix(ndofn_sld,ndofn_sld)

  if ( INOTMASTER ) then
     !
     ! Local Axes (Local --> Global)
     !
     do ipoin = 1, npoin
        ibopo = lpoty(ipoin)
        ievat = (ipoin-1)*ndofn_sld + 1
        jevat = (ipoin-1)*ndofn_sld + ndofn_sld
        if ( ibopo > 0 ) then
           if ( kfl_fixno_sld(1,ipoin) == 2_ip .or. &
                kfl_fixno_sld(1,ipoin) == 3_ip .or. &
                kfl_fixrs_sld(ipoin)   /= 0_ip ) then
              call sld_rotsys(2_ip, &
                   1_ip,1_ip,ndofn_sld,ndofn_sld,&
                   dummy_matrix,dunkn_sld(ievat:jevat),&
                   jacrot_du_dq_sld(1,1,ipoin),jacrot_dq_du_sld(1,1,ipoin))
           end if
        end if
     end do

  end if

end subroutine sld_roback
