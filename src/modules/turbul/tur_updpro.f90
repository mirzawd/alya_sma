!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_updpro()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_updpro
  ! NAME 
  !    tur_updpro
  ! DESCRIPTION
  !    This routine updates density and viscosity
  ! OUTPUT
  !    DETUR_TUR
  !    VITUR_TUR
  ! USED BY
  !    tur_begite
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_turbul
  implicit none
  integer(ip) :: ipoin
  real(rp)    :: dummr

  if( INOTMASTER .and. kfl_colev_tur /= 0 ) then
     !
     ! Assemble RHS
     !
     iunkn_tur = 1
     do ipoin = 1,2*npoin
        rhsid(ipoin) = 0.0_rp
     end do
     call tur_elmop2(12_ip)
     !
     ! Periodicity and Parall service
     !
     call rhsmod(2_ip,rhsid)
     !
     ! Update unknown
     !
     do ipoin=1,npoin
        dummr            = 1.0_rp / vmass(ipoin)
        detur_tur(ipoin) = rhsid( (ipoin-1)*2 + 1 ) * dummr
        vitur_tur(ipoin) = rhsid( (ipoin-1)*2 + 2 ) * dummr
     end do
     
  end if

end subroutine tur_updpro
