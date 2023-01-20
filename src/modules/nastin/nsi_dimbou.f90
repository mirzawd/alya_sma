!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_dimbou.f90
!> @author  David Oks 
!> @brief   Immersed boundary method for deformable bodies 
!> @details Immersed boundary method for deformable bodies 
!> @}
!------------------------------------------------------------------------

subroutine nsi_dimbou()

  use def_master, only : rhsid
  use def_kintyp, only : ip, rp 
  use def_domain, only : npoin, ndime 
  use def_coupli, only : kfl_dimbou
  use def_nastin, only : kfl_fixno_nsi, fsifo_nsi

  implicit none

  integer(ip) :: idime, ipoin, idofn

  if( kfl_dimbou ) then
     do ipoin = 1,npoin
        do idime = 1,ndime
           if( kfl_fixno_nsi(idime,ipoin) <= 0 ) then
              idofn = ndime*(ipoin-1) + idime
              rhsid(idofn) = rhsid(idofn) + fsifo_nsi(idime,ipoin)
           end if
        end do
     end do
  end if

end subroutine nsi_dimbou
