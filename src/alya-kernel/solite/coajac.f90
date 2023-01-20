!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine coajac(&
     nbvar,npoin,nrows,ia,ja,an,invdiag,pp,qq)
  !------------------------------------------------------------------------
  !****f* solite/coajac
  ! NAME 
  !    coajac
  ! DESCRIPTION
  ! USES
  ! USED BY 
  !***
  !------------------------------------------------------------------------
  use def_kintyp,         only      :  ip,rp
  use def_master,         only      :  INOTMASTER
  use def_solver,         only      :  solve_sol
  use def_solver,         only      :  SOL_DIAGONAL,SOL_LINELET
  use def_solver,         only      :  SOL_GAUSS_SEIDEL,SOL_BIDIAGONAL
  use mod_preconditioner, only      :  preconditioner_gauss_seidel
  use mod_preconditioner, only      :  preconditioner_bidiagonal_initialize
  use mod_direct_solver,  only      :  direct_solver_solution
  use mod_linelet,        only      :  linelet_solution
  use mod_solver,         only      :  solver_SpMV
  use mod_deflated_cg,    only      :  wvect
  use mod_deflated_cg,    only      :  wtvect
  implicit none
  integer(ip), intent(in)   :: nbvar,npoin,nrows
  integer(ip), intent(in)   :: ia(*),ja(*)
  real(rp),    intent(in)   :: an(*),invdiag(*),pp(*)  
  real(rp),    intent(out)  :: qq(*)
  integer(ip)               :: ii,kk,ngrou,igrou,nskyl,itpre
  real(rp),    pointer      :: rr(:)
  real(rp),    pointer      :: rrs(:)
  real(rp),    pointer      :: rrrhs(:)
  real(rp),    pointer      :: ww(:)

  ngrou = solve_sol(1) % ngrou
  nskyl = solve_sol(1) % nskyl
  itpre = solve_sol(1) % itpre

  !----------------------------------------------------------------------
  !
  ! Alloocate memory
  ! NROWS = 1 for MASTER
  !
  !----------------------------------------------------------------------

  allocate( rr(nrows) ) 
  if( solve_sol(1) % kfl_defpr == SOL_LINELET ) then
     allocate( ww(nrows) ) 
  end if
  allocate( rrs(nbvar*ngrou) )    
  do igrou = 1,nbvar*ngrou
     rrs(igrou)   = 0.0_rp
  end do
  if( solve_sol(1) % kfl_defas /= 0 ) then 
     allocate( rrrhs(nbvar*ngrou) )  
     do igrou = 1,nbvar*ngrou
        rrrhs(igrou) = rrs(igrou)
     end do
  end if

  !----------------------------------------------------------------------
  !
  ! Initial solution QQ and initial residual RR
  !
  !----------------------------------------------------------------------

  !if( solve_sol(1) % kfl_defpr == SOL_DIAGONAL .or. solve_sol(1) % kfl_defpr == SOL_GAUSS_SEIDEL ) then
  if( solve_sol(1) % kfl_defpr == SOL_DIAGONAL ) then
     !
     ! Jacobi with diagonal
     !
     do ii = 1,nrows
        qq(ii) = invdiag(ii) * pp(ii)
     end do

  else if( solve_sol(1) % kfl_defpr == SOL_GAUSS_SEIDEL ) then
     !
     ! SSOR
     !
     do ii = 1,nrows
        qq(ii) = 0.0_rp
     end do
!      call preconditioner_gauss_seidel(1_ip,npoin,nbvar,an,ja,ia,invdiag,qq,pp) !PCP: check with GH
     call bcsrgs(npoin,nbvar,an,ja,ia,invdiag,qq,pp)

  else if( solve_sol(1) % kfl_defpr == SOL_BIDIAGONAL ) then 
     !
     ! Bidiagonal
     !
     do ii = 1,nrows
        qq(ii) = 0.0_rp
     end do

     call preconditioner_bidiagonal_initialize(1_ip,npoin,nbvar,invdiag,qq,pp)
  
  else if( solve_sol(1) % kfl_defpr == SOL_LINELET ) then
     !
     ! Jacobi with linelet
     !
     do ii = 1,nrows
        qq(ii) = pp(ii)
     end do     
     call linelet_solution(nbvar,qq,ww,invdiag)  
  end if

  call solver_SpMV(solve_sol(1),an,qq,rr) 
  do ii = 1,nrows
     rr(ii) = pp(ii) - rr(ii)
  end do

  !----------------------------------------------------------------------
  !
  ! Smoothing:
  !
  ! 1. Jacobi smoothing: q^{i+1} = q^i + M^{-1} ( p - A q^i )
  !    1.1 Diagonal preconditioner: M = D
  !    1.2 Linelet preconditioner:  M = Linelet
  ! 2. Sym. Gauss-Seidel: (L+D) D^-1 (U+D) q = r
  !
  !----------------------------------------------------------------------

  do kk = 1,itpre

     if( solve_sol(1) % kfl_defpr == SOL_DIAGONAL ) then
        !
        ! Jacobi with diagonal
        !
        do ii = 1,nrows
           qq(ii) = qq(ii) + invdiag(ii) * rr(ii)
        end do

     else if( solve_sol(1) % kfl_defpr == SOL_GAUSS_SEIDEL ) then
        !
        ! SSOR
        !
        call bcsrgs(npoin,nbvar,an,ja,ia,invdiag,qq,pp) 

     else if( solve_sol(1) % kfl_defpr == SOL_LINELET ) then
        !
        ! Jacobi with linelet
        !
        call linelet_solution(nbvar,rr,ww,invdiag)  
        do ii = 1,nrows
           qq(ii) = qq(ii) + rr(ii) 
        end do

     end if

     call solver_SpMV(solve_sol(1),an,qq,rr)
     
     do ii = 1,nrows
        rr(ii) = pp(ii) - rr(ii)
     end do

  end do

  !----------------------------------------------------------------------
  !
  ! Coarse system: 
  ! 1. r     =  p - A q
  ! 2  r     => r'
  ! 3. A' e' =  r'
  ! 4. e'    => e
  ! 5. q     =  q + e
  !
  !----------------------------------------------------------------------

  call wtvect(npoin,ngrou,nbvar,rrs,rr) 

  if( INOTMASTER ) then
     do igrou = 1,ngrou*nbvar
        rrrhs(igrou) = rrs(igrou)
     end do
     call direct_solver_solution(solve_sol(1) % direct_solver_AMG,rrrhs,rrs) 
  end if

  call wvect(npoin,nbvar,rrs,rr)   

  do ii = 1,nrows
     qq(ii) = qq(ii) + rr(ii)   
  end do

  !----------------------------------------------------------------------
  !
  ! Smoothing
  !
  !----------------------------------------------------------------------

  do kk = 1,itpre

     call solver_SpMV(solve_sol(1),an,qq,rr)
     do ii = 1,nrows
        rr(ii) = pp(ii) - rr(ii)
     end do

     if( solve_sol(1) % kfl_defpr == SOL_DIAGONAL ) then
        !
        ! Jacobi with diagonal
        !
        do ii = 1,nrows
           qq(ii) = qq(ii) + invdiag(ii) * rr(ii)
        end do

     else if( solve_sol(1) % kfl_defpr == SOL_GAUSS_SEIDEL ) then        
        !
        ! SSOR
        !
        call bcsrgs(npoin,nbvar,an,ja,ia,invdiag,qq,pp) 

     else if( solve_sol(1) % kfl_defpr == SOL_LINELET ) then
        !
        ! Jacobi with linelet
        !
        call linelet_solution(nbvar,rr,ww,invdiag)  
        do ii = 1,nrows
           qq(ii) = qq(ii) + rr(ii) 
        end do

     end if
  end do

  !----------------------------------------------------------------------
  !
  ! Deallocate
  !
  !----------------------------------------------------------------------

  deallocate( rr )
  if( solve_sol(1) % kfl_defpr == SOL_LINELET ) then
     deallocate( ww )
  end if
  deallocate( rrs )
  if( solve_sol(1) % kfl_defas /= 0 ) then 
     deallocate( rrrhs )
  end if

end subroutine coajac
