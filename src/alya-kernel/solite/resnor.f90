!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine resnor(ndofn,kfl_symme,kfl_full_rows,rhsid,unkno,amatr,xnorm)
  !-----------------------------------------------------------------------
  !****f* master/resnor
  ! NAME 
  !    resnor
  ! DESCRIPTION
  !    This routine calls the resnors
  !    For diagonal solve which uses vmass, amatr must NOT be modified
  ! USES
  !    memchk
  !    mediso
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_master, only       :  NPOIN_REAL_12DI,INOTMASTER
  use def_domain, only       :  npoin,c_dom,r_dom,c_sym,r_sym
  use def_solver, only       :  solve_sol,memma
  use mod_solver, only       :  solver_SpMV
  use mod_solver, only       :  solver_parallel_vector_L2norm
  use mod_memchk
  implicit none
  integer(ip), intent(in)    :: ndofn,kfl_symme,kfl_full_rows
  real(rp),    intent(inout) :: unkno(*)
  real(rp),    intent(in)    :: amatr(*)
  real(rp),    intent(in)    :: rhsid(*)
  real(rp),    intent(out)   :: xnorm
  real(rp),    pointer       :: rr(:),ww(:)
  integer(ip)                :: iunkn,nunkn
  integer(4)                 :: istat
  real(rp)                   :: raux1,raux2

  if( INOTMASTER ) then
     nunkn = npoin * ndofn
  else
     nunkn = 0
  end if
  allocate(rr(max(nunkn,1_ip)),stat=istat)
  call memchk(0_ip,istat,memma,'RR','resnor',rr)
  allocate(ww(max(nunkn,1_ip)),stat=istat)
  call memchk(0_ip,istat,memma,'RR','resnor',rr)
  !
  ! RR = b - A x
  !
  call solver_SpMV(solve_sol(1),amatr,unkno,rr) 
  do iunkn = 1,nunkn
     rr(iunkn) = rhsid(iunkn) - rr(iunkn)
  end do
  !
  ! Diagonal D
  !
  if( kfl_symme == 0 ) then
     call diagon(npoin,ndofn,kfl_symme,kfl_full_rows,r_dom,c_dom,amatr,ww)
  else
     call diagon(npoin,ndofn,kfl_symme,kfl_full_rows,r_sym,c_sym,amatr,ww)
  end if
  do iunkn = 1,nunkn
     ww(iunkn) = 1.0_rp / ww(iunkn) 
  end do
  !
  ! RR = D^-1 ( b - A x )
  !
  do iunkn = 1,nunkn
     rr(iunkn) = rr(iunkn) * ww(iunkn)
  end do
  call solver_parallel_vector_L2norm(solve_sol(1),rr,raux2)
  !
  ! RR = D^-1 ( b - A x )
  !
  do iunkn = 1,nunkn
     rr(iunkn) = rhsid(iunkn) * ww(iunkn)
  end do
  call solver_parallel_vector_L2norm(solve_sol(1),rr,raux1)
  
  if( raux1 /= 0.0_rp ) then
     xnorm = raux2 / raux1
  else
     xnorm = 0.0_rp
  end if
  !
  ! Deallocate memory
  !
  call memchk(2_ip,istat,memma,'RR','resnor',rr)
  deallocate(rr,stat=istat)
  if(istat/=0) call memerr(2_ip,'RR','resnor',0_ip)
  call memchk(2_ip,istat,memma,'WW','resnor',ww)
  deallocate(ww,stat=istat)
  if(istat/=0) call memerr(2_ip,'WW','resnor',0_ip)
  
end subroutine resnor
