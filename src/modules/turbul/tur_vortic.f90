!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_vortic()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_vortic
  ! NAME
  !   tur_vortic
  ! DESCRIPTION
  !    Compute magnitude of vorticity or strain rate
  ! USES
  !    memgen
  ! USED BY
  !    tur_begite
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_turbul
  use mod_gradie
  implicit none
  integer(ip) :: ipoin
  real(rp)    :: xfact

  if( INOTMASTER .and. kfl_vorti_tur /= 0 ) then

     if( kfl_vorti_tur == 1 ) then        
        call tur_memarr(12_ip)
        kfl_vorti_tur = 2
     end if

     if( ipara_tur(1) == 0 ) then
        !
        ! Vorticity: W = sqrt( wij*wij ); wij = ( ui,j - uj,i )
        !
        if( ndime == 2 ) then
           call memgen(zero,1_ip,npoin)
        else
           call memgen(zero,ndime,npoin)
        end if
        vorti => gevec
        call vortic(2_ip)
        if( ndime == 2 ) then
           do ipoin = 1,npoin
              vorti_tur(ipoin) = sqrt(gevec(1,ipoin)*gevec(1,ipoin))
           end do
        else
           do ipoin = 1,npoin
              vorti_tur(ipoin) = gevec(1,ipoin)*gevec(1,ipoin) &
                   &           + gevec(2,ipoin)*gevec(2,ipoin) &
                   &           + gevec(3,ipoin)*gevec(3,ipoin) 
              vorti_tur(ipoin) = sqrt(vorti_tur(ipoin))
           end do
        end if       
        call memgen(two,1_ip,npoin)

     else 
        !
        ! Strain rate: S = sqrt(2.0*sij*sij); sij = 1/2 ( ui,j + uj,i )
        !
        call memgen(zero,ntens,npoin)
        call gradie(advec(1:ndime,1:npoin,1),gevec)
        if( ndime == 2 ) then
           do ipoin = 1,npoin
              xfact =         gevec(1,ipoin)*gevec(1,ipoin) & ! S11.S11
                   +          gevec(2,ipoin)*gevec(2,ipoin) & ! S22.S22
                   +   2.0_rp*gevec(3,ipoin)*gevec(3,ipoin)   ! 2*S12*S12
              xfact = 0.25_rp*xfact
              vorti_tur(ipoin) = sqrt(2.0_rp*xfact)
           end do
        else
           do ipoin = 1,npoin
              xfact =         gevec(1,ipoin)*gevec(1,ipoin) & ! S11.S11
                   +          gevec(2,ipoin)*gevec(2,ipoin) & ! S22.S22
                   +          gevec(4,ipoin)*gevec(4,ipoin) & ! S33.S33
                   +   2.0_rp*gevec(3,ipoin)*gevec(3,ipoin) & ! 2*S12*S12
                   +   2.0_rp*gevec(5,ipoin)*gevec(5,ipoin) & ! 2*S13*S13
                   +   2.0_rp*gevec(6,ipoin)*gevec(6,ipoin)   ! 2*S23*S23
              xfact = 0.25_rp*xfact
              vorti_tur(ipoin) = sqrt(2.0_rp*xfact)
           end do
        end if
        call memgen(two,ntens,npoin)

     end if

  end if
 
end subroutine tur_vortic
