!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine exm_doiter
  !-----------------------------------------------------------------------
  !****f* Exmedi/exm_doiter
  ! NAME 
  !    exm_doiter
  ! DESCRIPTION
  !   CORREGIR ESTO
  !    This routine solves an iteration of the linearized incompressible NS
  !    equations.
  !    Fractional step is to be implemented
  ! USES
  !    exm_comapp
  !    exm_iapupd
  !    exm_eapupd
  !    exm_rcpupd
  ! USED BY
  !    Exmedi
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_domain
  use      def_master
  use      def_solver

  use      def_exmedi
  use mod_timings, only : timings_ini
  use mod_timings, only : timings_end

  implicit none
  real(rp)    :: cpu_refe1,cpu_refe2

  call cputim(cpu_refe1)
  
  call timings_ini()
  call exm_begite
  call timings_end(ITASK_BEGITE)

  do while(kfl_goite_exm==1)

     call exm_solite

     call exm_endite(ITASK_ENDINN)      !  u(,ITER_K) <-- unkno  
     

  end do

  call timings_ini()
  call exm_endite(ITASK_ENDITE)     ! elmag(*,2) <-- elmag(*,1) , refhn(*,*,2) <-- refhn(*,*,1)
  call timings_end(ITASK_ENDITE)

  call cputim(cpu_refe2)
  cpu_exmed(1) = cpu_exmed(1) + cpu_refe2 - cpu_refe1 
  cpu_exmed(2) = cpu_exmed(2) + cpu_solve

end subroutine exm_doiter
