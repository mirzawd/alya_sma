!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_outite()
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_outite
  ! NAME 
  !    nsi_outite
  ! DESCRIPTION
  ! USES
  ! USED BY
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  implicit none
  !
  ! End of an inner iteration
  !
  if( ISEQUEN ) then
     if(itinn(modul)==kfl_psmat_nsi) then
        kfl_psmat_nsi=0
        call pspltm(&
             npoin,npoin,solve(ivari_nsi)%ndofn,0_ip,c_dom,r_dom,amatr,&
             trim(title)//': '//namod(modul),0_ip,18.0_rp,'cm',&
             0_ip,0_ip,2_ip,lun_psmat_nsi)
        call nsi_openfi(8_ip)
     end if
  end if

end subroutine nsi_outite
