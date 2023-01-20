!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine exm_concou
!-----------------------------------------------------------------------
!****f* Exmedi/exm_concou
! NAME 
!    exm_concou
! DESCRIPTION
!    This routine checks the convergence of this run and
!    set the general convergence flags.
! USED BY
!    Exmedi
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master

  use      def_exmedi
  use mod_outfor, only : outfor
  implicit none

  if(kfl_conve(modul)==1) then
     if(resid_exm(4)>cotol_exm) kfl_gocou = 1
  end if
  glres(modul) = resid_exm(4)

  coutp='InPotential'
  routp(1)=resid_exm(1)
  call outfor(9_ip,lun_outpu,' ')
  coutp='ExPotential'
  routp(1)=resid_exm(2)
  call outfor(9_ip,lun_outpu,' ')
  coutp='Global'
  routp(1)=resid_exm(4)
  call outfor(9_ip,lun_outpu,' ')

end subroutine exm_concou
