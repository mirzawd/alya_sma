!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_pstres()
  !------------------------------------------------------------------------
  !****f* Solidz/sld_output
  ! NAME 
  !    sld_pstres
  ! DESCRIPTION
  !    Output a postprocess variable
  ! USES
  !    postpr
  !    memgen
  ! USED BY
  !    sld_begste
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_solidz

  implicit none
  integer(ip)             :: ipoin,idime,jdime
  real(rp)                :: dummr

  if(INOTMASTER) then
  
     do ipoin = 1,npoin
        do idime = 1,ndime
           do jdime = 1,ndime
              nopio_sld(jdime+(idime-1)*ndime,ipoin) = 0.0_rp
           end do
        end do
     end do

     call sld_elmope(7_ip)
     call rhsmod(ndime*ndime,nopio_sld)

     do ipoin = 1,npoin
        dummr = 1.0_rp/vmass(ipoin)
        do idime = 1,ndime
           do jdime = 1,ndime
              nopio_sld(jdime+(idime-1)*ndime,ipoin) = dummr * nopio_sld(jdime+(idime-1)*ndime,ipoin) 
           end do
        end do
     end do
     
  end if

end subroutine sld_pstres
