!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_agmgsol(ff,xx,amatr,ia,ja,ndof,nbrows)
  !
  ! Sequential solve with agmg - this subroutine is obsolete I already put the sequential case inside mod_agmg2alya
  !
  use def_kintyp
  use def_solver

  implicit none

  integer(ip), intent(in)    :: ndof,nbrows
  integer(ip), intent(in)    :: ia(nbrows+1)
  real(rp),    intent(in)    :: amatr(ndof,ndof,ia(nbrows+1)-1)
  integer(ip), intent(in)    :: ja(ia(nbrows+1)-1)

  real(rp),    intent(inout)    :: xx(ndof*nbrows)
  real(rp),    intent(inout)    :: ff(ndof*nbrows)

  real(rp),allocatable       :: amatr_ndof(:),faux(:),xaux(:)
  integer(ip),allocatable    :: ja_ndof(:)
  integer(ip),allocatable    :: ia_ndof(:)

  integer(ip)          :: j,idof,irow
  integer(ip)          :: nz,neqn,i,nz_save

  integer(ip)          :: iprint,iter,kryl,korder
  integer(ip)          :: iblock(4)

  real(rp)             :: time1,time2,tol

  korder = 1 !1x,1y,2x,2y...
 
  neqn = nbrows*ndof
  nz  = ia(nbrows+1)-1
  nz_save = nz

  allocate(amatr_ndof(ndof*ndof*nz))
  allocate(ja_ndof(ndof*ndof*nz))
  allocate(ia_ndof(ndof*nbrows+1))
  allocate(faux(ndof*nbrows))
  allocate(xaux(ndof*nbrows))

  call matreorder(ff,xx,amatr,ia,ja,ndof,nbrows,korder,amatr_ndof,ja_ndof,ia_ndof,faux,xaux)
  !
  !       unit number for output messages: 6 => standard output
  iprint=666
  ! perhaps use teh one from alya
  !solve_sol(ivari) % lun_solve)    ! Output unit

  !
  ! obtain parameters from what has been defined in .dat
  !
  !kernel/master/soldef.f90:        solve_sol(ivari) % solco     = 1.0e-6_rp             ! Solver tolerance
  !kernel/master/soldef.f90:        solve_sol(ivari) % adres     = 0.1_rp                ! Adaptive residual tolerance
  !kernel/master/soldef.f90:        solve_sol(ivari) % solmi     = 1.0e-6_rp             ! Minimum solver tolerance
  !
  if (solve_sol(1) % wprob == 'CONTINUITY') then
     kryl = 1
  else
     kryl = solve_sol(1) % nkryd
  end if
  tol = solve_sol(1) % adres      ! here I use adaptive ratio instead of solco
  ! to make if comparable with what we typically use in Alya
  ! beware it must only be used for cases where the initial guess is 0
  iter = solve_sol(1) % miter
  !
  ! Beware I do not have a strict equivalence with adaptive in AGMG prepare runs without using adaptive so taht they are comparable.
  ! In pressure it is clear taht the tol must be put idem ratio (initial guess 0 - I discussed with Yvan) Mom see som logical value to use.
  !
  !
  ! SOLVE
  !
  !
  ! Missing see how to add the CPU Times for agmg directly into solve_sol(1) % cputi(??)
  !
  if (ndof >4 ) call runend('nsi_agmgsol: ndof>4')
  iblock(1) = 1
  iblock(2) = 1+nbrows
  iblock(3) = 1+2*nbrows
  iblock(4) = 1+3*nbrows
  call cputim(time1)
  !
#ifdef solve_w_agmg
  call dagmg(ndof*nbrows,amatr_ndof,ja_ndof,ia_ndof,faux,xaux,0,iprint,kryl,iter,tol) !non_block
  ! call dagmg(ndof*nbrows,amatr_ndof,ja_ndof,ia_ndof,faux,xaux,0,iprint,kryl,iter,tol,ndof,iblock) ! block
#endif  
  !
  call cputim(time2)

!  solve_sol(ivari) % cputi = !not sure how to use it

  write(666,'((a),2(i4,1x),(2e10.3,1x))') 'iter,kryl,tol,cputime'//solve_sol(1) % wprob,iter,kryl,tol,time2-time1
  ! you can get you results easily with: grep 'CONTI' -B 15 fort.666 | egrep -B 1 'CONTI|Convergence reached' |egrep  'CONTI|Iter'
  !
  ! now pass unk_agmg to unkno
  !
  if (korder==1) then
     xx(1:ndof*nbrows) = xaux(1:ndof*nbrows)
  else
     do idof = 1,ndof
        do irow = 1,nbrows
           i = idof+(irow-1)*ndof
           j = irow+(idof-1)*nbrows
           xx(i) = xaux(j)
        end do
     end do
  end if
       
  deallocate(amatr_ndof)
  deallocate(ja_ndof)
  deallocate(ia_ndof)
  deallocate(faux)
  deallocate(xaux)

 end subroutine nsi_agmgsol
