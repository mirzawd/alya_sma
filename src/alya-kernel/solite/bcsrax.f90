!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Algebraic_Solver
!> @addtogroup Krylov_Solver
!> @{
!> @file    bcsrax.f90
!> @author  Guillaume Houzeaux
!> @brief   SpMV
!> @details Parallel SpMV
!>          Multiply a matrix by a vector YY = A XX
!>          where A is in CSR, COO or ELL format. This is the parallel
!>          version, including blocking and non-blocking send-receive
!>          MPI functions.
!>
!>          INPUT
!>             NBNODES .... Number of equations
!>             NBVAR ...... Number of variables
!>             AN ......... Matrix
!>             JA ......... List of elements
!>             IA ......... Pointer to list of elements
!>             XX ......... Vector
!>          OUTPUT
!>             YY ......... result vector
!>
!> @}
!------------------------------------------------------------------------

subroutine bcsrax(itask,nbnodes,nbvar,an,ja,ia,xx,yy)

  use def_kintyp,         only     :  ip,rp
  use def_solver,         only     :  solve_sol
  use mod_solver,         only     :  solver_parallel_SpMV
  implicit none
  integer(ip), intent(in)          :: itask,nbnodes,nbvar
  real(rp),    intent(in)          :: an(nbvar,nbvar,*)
  integer(ip), intent(in)          :: ja(*),ia(*)
  real(rp),    intent(inout)       :: xx(nbvar,*)
  real(rp),    intent(out), target :: yy(nbvar,*)

  if( itask == 0 ) then
     call solver_parallel_SpMV(solve_sol(1),an,xx,yy,OPENMP=.true.,MPI=.false.,TIMING=.true.)
  else
     if( solve_sol(1) % omp_interface == 0 ) then
        call solver_parallel_SpMV(solve_sol(1),an,xx,yy,OPENMP=.true.,TIMING=.true.,OPENMP_INTERFACE=.false.)
     else
        call solver_parallel_SpMV(solve_sol(1),an,xx,yy,OPENMP=.true.,TIMING=.true.)
     end if
  end if

end subroutine bcsrax
