!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_concou
!-----------------------------------------------------------------------
!****f* Temper/tur_concou
! NAME 
!    tur_concou
! DESCRIPTION
!    This routine checks the temperature convergence of the run.
! USED BY
!    Temper
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_turbul
  use mod_outfor, only : outfor
  implicit none
  !
  ! Check convergence
  !
  if(kfl_conve(modul)==1) then
     if(resid_tur>cotol_tur) kfl_gocou = 1
  end if
  glres(modul) = resid_tur
  !
  ! Output residuals
  !
  coutp(1)='TURBULENCE'
  routp(1)=resid_tur
  call outfor(9_ip,lun_outpu,' ')
    
end subroutine tur_concou
