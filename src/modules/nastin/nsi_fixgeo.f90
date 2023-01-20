!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_fixgeo(icode,ifixx,ifixy,ifixz,ifixr,kfl_value,bvalu)
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_fixgeo
  ! NAME 
  !    nsi_fixgeo
  ! DESCRIPTION
  ! USES
  ! USED BY
  !    nsi_reabcs
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_nastin 
  use mod_chktyp, only : check_type
  implicit none
  integer(ip), intent(in) :: icode,ifixx,ifixy,ifixz,ifixr
  integer(ip), intent(in) :: kfl_value
  real(rp),    intent(in) :: bvalu(3)
  integer(ip)             :: ipoin,ibopo
  real(rp),    pointer    :: xvalu(:,:,:)

  if( kfl_value /= 0 ) then
     xvalu => xfiel(kfl_value) % a
     call check_type(xfiel,kfl_value,ndime,npoin)
  end if

  if( ndime == 2 ) then

     if( kfl_value == 0 ) then
        do ipoin = 1,npoin
           if( lpoty(ipoin) > 0 ) then
              ibopo = lpoty(ipoin)
              if( kfl_geono(ibopo) == icode ) then
                 kfl_fixno_nsi(1,ipoin) = ifixx
                 kfl_fixno_nsi(2,ipoin) = ifixz
                 bvess_nsi(1,ipoin,1)   = bvalu(1)
                 bvess_nsi(2,ipoin,1)   = bvalu(2)
                 if( ibopo /= 0 ) kfl_fixrs_nsi(ipoin) = ifixr
              end if
           end if
        end do
     else
        do ipoin = 1,npoin
           if( lpoty(ipoin) > 0 ) then
              ibopo = lpoty(ipoin)
              if( kfl_geono(ibopo) == icode ) then
                 kfl_fixno_nsi(1,ipoin) = ifixx
                 kfl_fixno_nsi(2,ipoin) = ifixz
                 bvess_nsi(1,ipoin,1)   = xvalu(1,ipoin,1)
                 bvess_nsi(2,ipoin,1)   = xvalu(2,ipoin,1)
                 if( ibopo /= 0 ) kfl_fixrs_nsi(ipoin) = ifixr
              end if
           end if
        end do
     end if

  else

     if( kfl_value == 0 ) then
        if( kfl_modfi_nsi /= 1 ) then
           do ipoin = 1,npoin
              if( lpoty(ipoin) > 0 ) then
                 ibopo = lpoty(ipoin)
                 if( kfl_geono(ibopo) == icode ) then
                    kfl_fixno_nsi(1,ipoin) = ifixx
                    kfl_fixno_nsi(2,ipoin) = ifixy
                    kfl_fixno_nsi(3,ipoin) = ifixz
                    bvess_nsi(1,ipoin,1)   = bvalu(1)
                    bvess_nsi(2,ipoin,1)   = bvalu(2)
                    bvess_nsi(3,ipoin,1)   = bvalu(3)
                    if( ibopo /= 0 ) kfl_fixrs_nsi(ipoin) = ifixr
                 end if
              end if
           end do
        else if( kfl_modfi_nsi == 1 ) then
           do ipoin = 1,npoin
              if( lpoty(ipoin) > 0 ) then
                 ibopo = lpoty(ipoin)
                 if( kfl_geono(ibopo) == icode ) then
                    ibopo                  = lpoty(ipoin)
                    kfl_fixno_nsi(2,ipoin) = ifixy
                    if((ifixx==1).and.(ifixy==0).and.(ifixz==1).and.(abs(coord(1,ipoin)-xfree_nsi)<0.0001_rp)) then
                       kfl_fixno_nsi(1,ipoin) = 0_ip
                       kfl_fixno_nsi(3,ipoin) = 0_ip
                    else
                       kfl_fixno_nsi(1,ipoin) = ifixx
                       kfl_fixno_nsi(3,ipoin) = ifixz
                    end if
                    bvess_nsi(1,ipoin,1)   = bvalu(1)
                    bvess_nsi(2,ipoin,1)   = bvalu(2)
                    bvess_nsi(3,ipoin,1)   = bvalu(3)
                    if( ibopo /= 0 ) kfl_fixrs_nsi(ipoin) = ifixr
                 end if
              end if
           end do
        end if
     else
        if( kfl_modfi_nsi /= 1 ) then
           do ipoin = 1,npoin
              if( lpoty(ipoin) > 0 ) then
                 ibopo = lpoty(ipoin)
                 if( kfl_geono(ibopo) == icode ) then
                    kfl_fixno_nsi(1,ipoin) = ifixx
                    kfl_fixno_nsi(2,ipoin) = ifixy
                    kfl_fixno_nsi(3,ipoin) = ifixz
                    bvess_nsi(1,ipoin,1)   = xvalu(1,ipoin,1)
                    bvess_nsi(2,ipoin,1)   = xvalu(2,ipoin,1)
                    bvess_nsi(3,ipoin,1)   = xvalu(3,ipoin,1)
                    if( ibopo /= 0 ) kfl_fixrs_nsi(ipoin) = ifixr
                 end if
              end if
           end do
        else if( kfl_modfi_nsi == 1 ) then
           do ipoin = 1,npoin
              if( lpoty(ipoin) > 0 ) then
                 ibopo = lpoty(ipoin)
                 if( kfl_geono(ibopo) == icode ) then
                    kfl_fixno_nsi(2,ipoin) = ifixy
                    if((ifixx==1).and.(ifixy==0).and.(ifixz==1).and.(abs(coord(1,ipoin)-xfree_nsi)<0.0001_rp)) then
                       kfl_fixno_nsi(1,ipoin) = 0_ip
                       kfl_fixno_nsi(3,ipoin) = 0_ip
                    else
                       kfl_fixno_nsi(1,ipoin) = ifixx
                       kfl_fixno_nsi(3,ipoin) = ifixz
                    end if
                    bvess_nsi(1,ipoin,1)      = xvalu(1,ipoin,1)
                    bvess_nsi(2,ipoin,1)      = xvalu(2,ipoin,1)
                    bvess_nsi(3,ipoin,1)      = xvalu(3,ipoin,1)
                    if( ibopo /= 0 ) kfl_fixrs_nsi(ipoin) = ifixr
                 end if
              end if
           end do
        end if
     end if

  end if

end subroutine nsi_fixgeo
