!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine soleig(amatr,eigva,eigen,bmatr,iter)
  !-----------------------------------------------------------------------
  !****f* master/soleig
  ! NAME 
  !    solver
  ! DESCRIPTION
  !    This routine calls the solvers
  ! USES
  !    memchk
  !    mediso
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  rp,ip
  use def_solver, only     :  cpu_eigen
  implicit none
  real(rp),    intent(in)  :: eigva(*),eigen(*)
  real(rp),    intent(in)  :: amatr(*),bmatr(*) 
  integer(ip), intent(in)  :: iter
  real(rp)                 :: time1,time2

  call cputim(time1)
  
  call runend('SOLEIG: NO MORE EIGEN SOLVER')

  call cputim(time2)
  cpu_eigen = cpu_eigen + (time2-time1)

end subroutine soleig
