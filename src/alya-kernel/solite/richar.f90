!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine richar(&
     nbnodes, nbvar, idprecon, maxiter, eps, an, pn, &
     kfl_cvgso, lun_cvgso, kfl_solve, lun_outso,&
     ja, ia, bb, xx )  
  !-----------------------------------------------------------------------
  ! Objective: Solve    [A] [M]^-1  [M] x = b
  !                        [A']        x' = b'   by the BiCGstab method.
  !
  !            Three preconditioners are possible:
  !            idprecon = 0 => [M]^-1 = [M] = I
  !            idprecon = 2 => [M]^-1 = diag([A])
  !            idprecon = 3 => [M]^-1 = SPAI

  !            The diagonal terms of [A] must be the first in each row.
  !
  !            If Diag. Scaling is selected the preconditioned system is:
  !                    D^-1/2 [A] D^-1/2   D^1/2 x  =  D^-1/2 b
  !                          [A']              x'   =      b'
  !
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp,lg
  use def_master, only       :  IMASTER,INOTSLAVE
  use def_solver, only       :  memit,resi1,resi2,solve_sol
  use mod_solver, only       :  solver_SpMV
  use mod_solver, only       :  solver_parallel_vector_L2norm

  use mod_memchk
  implicit none
  integer(ip), intent(in)    :: nbnodes,nbvar,idprecon,maxiter
  integer(ip), intent(in)    :: kfl_cvgso,lun_cvgso
  integer(ip), intent(in)    :: kfl_solve,lun_outso
  real(rp),    intent(inout) :: eps
  real(rp),    intent(in)    :: an(nbvar,nbvar,*),pn(*)
  integer(ip), intent(in)    :: ja(*),ia(*)
  real(rp),    intent(in)    :: bb(*)
  real(rp),    intent(inout) :: xx(*)
  real(rp),    pointer       :: rr(:),zz(:),ww(:)
  real(rp),    pointer       :: invdiag(:)
  integer(ip)                :: ii,nrows,ierr,npoin
  integer(4)                 :: istat
  real(rp)                   :: stopcri,invnb,resid,dummr,newrho

  if( IMASTER ) then
     nrows = 1                ! Minimum memory for working arrays
     npoin = 0                ! master does not perform any loop
  else
     npoin = nbnodes
     nrows = nbnodes * nbvar
  end if
  !
  ! Allocate memory for working arrays
  !
  allocate(rr(nrows),stat=istat) 
  call memchk(0_ip,istat,memit,'RR','richar',rr)
  allocate(zz(nrows),stat=istat) 
  call memchk(0_ip,istat,memit,'ZZ','richar',zz)
  allocate(ww(nrows),stat=istat) 
  call memchk(0_ip,istat,memit,'WW','richar',ww)
  allocate(invdiag(nrows),stat=istat) 
  call memchk(0_ip,istat,memit,'INVDIAG','bcgpls',invdiag)

  !----------------------------------------------------------------------
  !
  ! Initial computations
  !
  !----------------------------------------------------------------------

  call solope(&
     1_ip, nbvar, idprecon, eps, an, pn, ja, ia, bb, xx , &
     ierr, stopcri, newrho, resid, invnb, rr, zz, dummr, ww, &
     invdiag, dummr )
  if( ierr /= 0 ) goto 10

  !----------------------------------------------------------------------
  !
  ! MAIN LOOP
  !
  !----------------------------------------------------------------------

  do while( solve_sol(1) % iters < maxiter .and. resid > stopcri )
     !
     ! x = x + z
     !
     do ii = 1,nrows
        xx(ii) = xx(ii) + zz(ii)
     end do
     !
     ! z = L^-1 ( b - A x )
     !
     call solver_SpMV(solve_sol(1),an,xx,rr) 

     do ii = 1,nrows
        rr(ii) = bb(ii) - rr(ii)
     end do
     call precon(&
          3_ip,nbvar,nbnodes,nrows,solve_sol(1)%kfl_symme,idprecon,ia,ja,an,&
             pn,invdiag,ww,rr,zz)
     !
     ! || z ||
     !
     call solver_parallel_vector_L2norm(solve_sol(1),zz,resid)
     resi2 = resi1
     resi1 = resid * invnb           
     solve_sol(1) % iters = solve_sol(1) % iters + 1
     !
     ! Convergence output
     !    
     if( kfl_cvgso == 1 ) call outcso(an,bb,xx)

  end do

  !----------------------------------------------------------------------
  !
  ! END MAIN LOOP
  !
  !----------------------------------------------------------------------
 

10 continue
  call solope(&
       2_ip, nbvar, idprecon, dummr, an, dummr, ja, ia, &
       bb,xx ,ierr, dummr, dummr, dummr, invnb, rr, dummr, &
       dummr, dummr, invdiag, dummr )

  if( kfl_solve == 1 .and. INOTSLAVE ) then
     if( ierr == 1 ) write(lun_outso,201) solve_sol(1) % iters
     if( ierr == 2 ) write(lun_outso,202) solve_sol(1) % iters
  end if

  call memchk(2_ip,istat,memit,'INVDIAG','richar',invdiag)
  deallocate(invdiag,stat=istat)
  if(istat/=0) call memerr(2_ip,'INVDIAG','richar',0_ip)

  call memchk(2_ip,istat,memit,'WW','richar',ww)
  deallocate(ww,stat=istat)
  if(istat/=0) call memerr(2_ip,'WW','richar',0_ip)

  call memchk(2_ip,istat,memit,'ZZ','richar',zz)
  deallocate(zz,stat=istat)
  if(istat/=0) call memerr(2_ip,'ZZ','richar',0_ip)

  call memchk(2_ip,istat,memit,'RR','richar',rr)
  deallocate(rr,stat=istat)
  if(istat/=0) call memerr(2_ip,'RR','richar',0_ip)

100 format(i7,1x,e12.6)
110 format(i5,18(2x,e12.6))
201 format(&
       & '# Error at iteration ',i6,&
       & 'Dividing by zero: alpha = beta = (rho^k/rho^{k-1})*(alpha/omega)')
202 format(&
       & '# Error at iteration ',i6,&
       & 'Dividing by zero: alpha = beta = alpha = rho^k / <r0,q^{k+1}>')

end subroutine richar
