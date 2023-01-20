!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_denvis()
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_denvis
  ! NAME 
  !    nsi_denvis
  ! DESCRIPTION
  !    Smoothing of density and viscosity
  ! USES
  !    postpr
  !    memgen
  ! USED BY
  !    nsi_output
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_nastin
  use mod_ker_proper 
  use mod_memory,  only : memory_alloca, memory_deallo
  implicit none
  integer(ip)       :: ipoin,dummi
  real(rp), pointer :: prope_tmp(:) 

  if( INOTMASTER ) then 

     nullify ( prope_tmp )
     call memory_alloca(mem_modul(1:2,modul),'PROPE_TMP','nsi_denvis',prope_tmp,npoin)
     call ker_proper('DENSI','NPOIN',dummi,dummi,prope_tmp)
     do ipoin = 1,npoin
        prope_nsi(1,ipoin) = prope_tmp(ipoin)
     end do
     call ker_proper('VISCO','NPOIN',dummi,dummi,prope_tmp)
     do ipoin = 1,npoin
        prope_nsi(2,ipoin) = prope_tmp(ipoin)
     end do
     call memory_deallo(mem_modul(1:2,modul),'PROPE_TMP','nsi_denvis',prope_tmp)
     
  end if

end subroutine nsi_denvis
