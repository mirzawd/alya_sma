!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_fixgeo(icode,ifixx,kfl_value,bvalu)
  !------------------------------------------------------------------------
  !****f* Turbul/tur_fixgeo
  ! NAME 
  !    tur_fixgeo
  ! DESCRIPTION
  ! USES
  ! USED BY
  !    tur_reabcs
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_turbul
  use mod_chktyp, only : check_type
  implicit none
  integer(ip), intent(in) :: icode,ifixx
  integer(ip), intent(in) :: kfl_value
  real(rp),    intent(in) :: bvalu(4)
  integer(ip)             :: ipoin,ibopo,iturb
  
  if( kfl_value == 0 ) then
     do ipoin = 1,npoin
        if( lpoty(ipoin) > 0 ) then
           ibopo = lpoty(ipoin)
           if( kfl_geono(ibopo) == icode ) then
              ibopo = lpoty(ipoin)
              iturb = iunkn_tur
              kfl_fixno_tur(1,ipoin,iturb) = ifixx
              bvess_tur(1,ipoin,iturb)     = bvalu(iturb)
           end if
        end if
     end do
  else if( kfl_value > 0 ) then 
     call check_type(xfiel,kfl_value,1_ip,npoin)
     do ipoin = 1,npoin
        if( lpoty(ipoin) > 0 ) then
           ibopo = lpoty(ipoin)
           if( kfl_geono(ibopo) == icode ) then
              ibopo = lpoty(ipoin)
              iturb = iunkn_tur
              if (ifixx.ne.0)  kfl_fixno_tur(1,ipoin,iturb) = 1
              bvess_tur(1,ipoin,iturb) =  xfiel(kfl_value) % a(1,ipoin,1)
           end if
        end if
     end do    
  else
     call runend('TUR_FIXGEO: WRONG FIELD')
  end if

end subroutine tur_fixgeo
