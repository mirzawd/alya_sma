!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_concou
!-----------------------------------------------------------------------
!****f* Nastin/nsi_concou
! NAME 
!    nsi_concou
! DESCRIPTION
!    This routine checks the Nastin convergence of the run and
!    set the general convergence flags.
! USED BY
!    Nastin
!***
!-----------------------------------------------------------------------
  use      def_master
  use      def_nastin
  use mod_outfor, only : outfor
  implicit none
  !
  ! Check convergence
  !
  glres(modul) = max(resid_nsi,resip_nsi)
  if( kfl_conve(modul) == 1 ) then
     if( glres(modul) > cotol_nsi ) kfl_gocou = 1
  end if   
  !
  ! Output residuals
  !
  coutp(1) = 'VELOCITY'
  routp(1) = resid_nsi
  call outfor(9_ip,lun_outpu,' ')
  coutp(1) = 'PRESSURE'
  routp(1) = resip_nsi
  call outfor(9_ip,lun_outpu,' ')

  call nsi_coupli(ITASK_CONCOU) 

end subroutine nsi_concou
