!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine steepe( &
     nbnodes, nbvar, idprecon, maxiter, &
     eps, an, pn, kfl_cvgso, lun_cvgso, kfl_solve, &
     lun_outso, ja, ia, bb, xx )
  !-----------------------------------------------------------------------
  !****f* solite/steepe
  ! NAME
  !    steepe
  ! DESCRIPTION
  !    Steepest descent algebraic solver
  ! USES
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp,lg
  use def_master, only       :  IMASTER
  use def_solver, only       :  memit,SOL_NO_PRECOND,SOL_LINELET,&
       &                        SOL_SQUARE,SOL_DIAGONAL,SOL_MATRIX,&
       &                        resi1,resi2,solve_sol
  use mod_memchk
  implicit none
  integer(ip), intent(in)    :: nbnodes, nbvar, idprecon, maxiter
  integer(ip), intent(in)    :: kfl_cvgso, lun_cvgso
  integer(ip), intent(in)    :: kfl_solve, lun_outso
  real(rp),    intent(in)    :: eps
  real(rp),    intent(in)    :: an(nbvar,nbvar,*), pn(*)
  integer(ip), intent(in)    :: ja(*), ia(*)
  real(rp),    intent(in)    :: bb(*)
  real(rp),    intent(inout) :: xx(nbvar*nbnodes)
  integer(ip)                :: nrows, npoin
  integer(4)                 :: istat
  real(rp)                   :: stopcri, resid
  real(rp)                   :: invnb
  real(rp),    pointer       :: rr(:), invdiag(:)

#ifdef EVENT
  call mpitrace_user_function(1)
#endif

  if(  IMASTER ) then
     nrows   = 1                ! Minimum memory for working arrays
     npoin   = 0                ! master does not perform any loop
  else
     nrows   = nbnodes * nbvar
     npoin   = nbnodes
  end if
  !
  ! Sequential and slaves: Working arrays
  !
  allocate(rr(nrows),stat=istat) 
  call memchk(0_ip,istat,memit,'RR','deflcg',rr)

  allocate(invdiag(nrows),stat=istat)
  call memchk(0_ip,istat,memit,'INVDIAG','deflcg',invdiag)

  if( IMASTER ) nrows = 0 ! Master does not perform any loop
  solve_sol(1) % iters = 0

  !-----------------------------------------------------------------
  !
  !  MAIN LOOP
  !
  !-----------------------------------------------------------------

  resid = 1.0_rp !TODO: might need correction, it was not initialized
  invnb = 1.0_rp !TODO: might need correction, it was not initialized
  stopcri = 0.0_rp !TODO: might need correction, it was not initialized
  do while( solve_sol(1) % iters < maxiter .and. resid > stopcri ) 

     !do ipoin = 1,npoin
     !resid  = || rk ||
     resi2  = resi1
     resi1  = resid * invnb
     solve_sol(1) % iters  = solve_sol(1) % iters + 1
     !
     ! Convergence post process:
     ! kk    = iteration number
     ! resi1 = preconditioned residual
     !
     if( kfl_cvgso == 1 ) &
          call outcso(an,bb,xx)

  end do

  !-----------------------------------------------------------------
  !
  ! Deallocate memory
  !
  !-----------------------------------------------------------------

  call memchk(2_ip,istat,memit,'INVDIAG','deflcg',invdiag)
  deallocate(invdiag,stat=istat)
  if(istat/=0) call memerr(2_ip,'INVDIAG','deflcg',0_ip)

  call memchk(2_ip,istat,memit,'RR','deflcg',rr)
  deallocate(rr,stat=istat)
  if(istat/=0) call memerr(2_ip,'RR','deflcg',0_ip)

110 format(i5,18(2x,e12.6))
120 format(&
       & '# Error at iteration ',i6,&
       & 'Dividing by zero: alpha = rho^k / <p^{k+1},q^{k+1}>')

#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine steepe
