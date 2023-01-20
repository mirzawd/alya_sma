!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_matrix()
  !------------------------------------------------------------------------
  !
  ! Compute elemental matrix and RHS
  !
  !  USED BY
  !    tur_solite
  !  USES
  !     tur_inisol
  !     tur_elmop2
  !     tur_bouope
  !------------------------------------------------------------------------
  use def_master
  use def_domain
  use def_turbul
  use def_solver
  implicit none
  real(rp)    :: time1,time2
  !
  ! Initialization
  !
  call inisol()

  call cputim(time1)

  if( INOTMASTER ) then
     call tur_elmop2(1_ip)
     !call tur_elmope(1_ip)
     call tur_bouope()
  end if

  call cputim(time2) 
  cpu_modul(CPU_ASSEMBLY,modul) = cpu_modul(CPU_ASSEMBLY,modul) + time2 - time1
  solve_sol(1) % kfl_assem = SOL_MATRIX_HAS_CHANGED ! Matrix has been assembled

end subroutine tur_matrix
