!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_concou
!-----------------------------------------------------------------------
!****f* Wavequ/lev_concou
! NAME 
!    lev_concou
! DESCRIPTION
!    This routine checks the level set convection convergence of the run.
! USED BY
!    Wavequ
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_levels
  use mod_outfor, only : outfor
  implicit none

  !
  ! Check convergence
  !
  if(kfl_conve(modul)==1) then
     if(resid_lev>cotol_lev) kfl_gocou = 1
  end if
  glres(modul) = resid_lev
  !
  ! Output residuals
  !
  coutp(1)='LEVEL SET'
  routp(1)=resid_lev
  call outfor(9_ip,lun_outpu,' ')

end subroutine lev_concou
