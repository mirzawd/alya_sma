!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_turnof
  !-----------------------------------------------------------------------
  !****f* Levels/lev_turnof
  ! NAME 
  !    lev_turnof
  ! DESCRIPTION
  !    This routine closes the run for the module levels
  ! USES
  !    lev_outcpu
  !    lev_output
  ! USED BY
  !    Wavequ
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_levels
  use      def_solver
  use mod_outfor, only : outfor
  implicit none
  real(rp)             :: xfact
  !
  ! Output results
  !
  if((kfl_inlev_lev==2).or.(kfl_inlev_lev==4)) call lev_l1norm()

  call lev_openfi(4_ip)
  !
  ! Output cputimes - for the momment only time for ITASK_ENDSTE that is where redist is done - format similar to nsi_outcpu - perhaps create lev_outcpu
  !
  if( INOTSLAVE ) then
     !
     ! Initializations
     !
     routp(1) = cpu_modul(30,modul)
     call outfor(29_ip,momod(modul)%lun_outpu,' ')

     if( cpu_modul(30,modul) > 0.0_rp ) then
        xfact = 100.0_rp / routp(1)
     else
        xfact = 1.0_rp
     end if
     !
     ! Others
     !
     coutp(1)  = 'ITASK_ENDSTE'
     routp(1)  = cpu_modul(ITASK_ENDSTE,modul) 
     routp(2)  = xfact*routp(1)
     call outfor(30_ip,momod(modul)%lun_outpu,' ')

  end if


end subroutine lev_turnof

