!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Algebraic_Solver
!> @{
!> @file    solric.f90
!> @date    01/01/2005
!> @author  Guillaume Houzeaux
!> @brief   Richardson with diagonal preconditioning
!> @details Richardson based solvers:
!>
!>          DIAGONAL ................ KFL_ALGSO = SOL_SOLVER_RICHARDSON
!>          MATRIX_BASED_DIAGONAL ... KFL_ALGSO = SOL_SOLVER_MATRIX_RICHARDSON
!>
!>          Some definitions:
!>          M ....... Mass matrix (computed with user-defined integration rule)
!>          Mc ...... Close mass matrix (computed with close rule)
!>          P ....... Diagonal module-defined preconditioner
!>          xdiag ... Module-define constant (basically the minimum time 
!>                    step for explicit schemes)
!>
!>          The different options are:
!>
!>          MASS MATRIX:        x^{i+1} = x^i + xdiag.M^-1.r^i
!>          CLOSE MASS MATRIX:  x^{i+1} = x^i + xdiag.Mc^-1.r^i
!>          LOCAL MASS MATRIX:  x^{i+1} = x^i + P.r^i
!>          LOCAL DIAGONAL:     x^{i+1} = x^i + P.M^-1.r^i
!>          AUGMENTED DIAGONAL: x^{i+1} = (M/xdiag+P)^-1.(M/xdiag.x^i+r^i)
!>
!>          For example, for an explicit cheme, choose xdiag = dt_min and
!>          the MASS MATRIX option.
!> 
!>          According to the solver, the residual is already assembled or
!>          the matrix contribution must be added r^i <= r^i - A.x^i
!> @} 
!-----------------------------------------------------------------------

subroutine solric(nbnodes,ndofn,rhsid,unkno,amatr,pmatr)
  use def_kintyp
  use def_master, only       :  INOTMASTER,IMASTER
  use def_domain, only       :  vmass,vmasc
  use def_solver, only       :  solve_sol
  use def_solver, only       :  resin,resfi,resi1,resi2,memit,&
       &                        SOL_SOLVER_MATRIX_RICHARDSON,&
       &                        SOL_MASS_MATRIX,&
       &                        SOL_CLOSE_MASS_MATRIX,&
       &                        SOL_LOCAL_MASS_MATRIX,&
       &                        SOL_LOCAL_DIAGONAL,&
       &                        SOL_AUGMENTED_DIAGONAL
  use mod_solver, only       :  solver_SpMV
  use mod_solver, only       :  solver_parallel_vector_L2norm
  use mod_solver, only       :  solver_lumped_mass_system
  use mod_solver, only       :  solver_parallel_SpMV
  use mod_memory
  implicit none
  integer(ip), intent(in)            :: nbnodes
  integer(ip), intent(in)            :: ndofn
  real(rp),    intent(inout)         :: unkno(*)
  real(rp),    intent(in)            :: amatr(1,1,*)
  real(rp),    intent(in)            :: pmatr(*)
  real(rp),    intent(inout), target :: rhsid(1:nbnodes*ndofn)
  integer(ip)                        :: itotn,npoin
  integer(ip)                        :: idofn,ipoin
  integer(ip)                        :: maxiter,nrows
  real(rp)                           :: fact1,xdiag,fact2
  real(rp)                           :: bnorm,invnb,eps,resid,stopcri
  real(rp),    pointer               :: rr(:)

  xdiag   = solve_sol(1) % xdiag
  maxiter = solve_sol(1) % miter
  eps     = solve_sol(1) % solmi
  solve_sol(1) % iters   = 0
  resid   = 1.0_rp
  if( IMASTER ) then
     npoin = 0
  else
     npoin = nbnodes
  end if
  nrows = ndofn * npoin

  nullify(rr)
  if( solve_sol(1) % kfl_algso == SOL_SOLVER_MATRIX_RICHARDSON ) then
     call memory_alloca(memit,'RR','solric',rr,ndofn*max(1_ip,npoin),'DO_NOT_INITIALIZE')
  else
     rr => rhsid
  end if 
  !
  ! Stoping criteria. We have two cases:
  ! 1. Matrix is  available: Residual = || r || / || b ||
  ! 2. Matrix not available: Residual = || r || 
  ! 
  if( solve_sol(1) % kfl_algso == SOL_SOLVER_MATRIX_RICHARDSON ) then
     call solver_parallel_vector_L2norm(solve_sol(1),rhsid,bnorm)
     stopcri = eps * bnorm
     invnb   = 1.0_rp / ( epsilon(1.0_rp) + bnorm )
  else
     bnorm   = 1.0_rp
     stopcri = eps
     invnb   = 1.0_rp
  end if
  if( maxiter == 1 ) stopcri = 0.0_rp

  !----------------------------------------------------------------------
  !
  ! SOLVER ITERATIONS
  !
  !----------------------------------------------------------------------

  do while( solve_sol(1) % iters < maxiter .and. resid > stopcri )
     !
     ! Matrix-based diaognal solver: r^i <= r^i - A.x^i
     !
     if( solve_sol(1) % kfl_algso == SOL_SOLVER_MATRIX_RICHARDSON .and. INOTMASTER ) then

        if( solve_sol(1) % omp_interface == 0 ) then
           call solver_parallel_SpMV(solve_sol(1),amatr,unkno,rr,OPENMP=.true.,TIMING=.true.,OPENMP_INTERFACE=.false.)
        else
           call solver_parallel_SpMV(solve_sol(1),amatr,unkno,rr,OPENMP=.true.,TIMING=.true.)
        end if
        !call solver_SpMV(solve_sol(1),amatr,unkno,rr)
        do ipoin = 1,npoin * ndofn
           rr(ipoin) = rhsid(ipoin) - rr(ipoin)
        end do
     end if
     !
     ! Residual norm: || r^i ||
     !
     call solver_parallel_vector_L2norm(solve_sol(1),rr,resid)
     if( solve_sol(1) % iters == 0 ) then
        resi1 = resid * invnb
        resi2 = resi1 
        resin = resi1
     else
        resi2  = resi1
        resi1  = resid * invnb
     end if
     !
     ! Solve system: u^{k+1} = u^k + D^{-1} r^k 
     !
     if( solve_sol(1) % kfl_preco == SOL_MASS_MATRIX ) then

        !----------------------------------------------------------------
        !
        ! MASS MATRIX [MASSM]:
        !
        ! x^{i+1} = x^i + xdiag.M^-1.r^i
        !
        !----------------------------------------------------------------

        if( ndofn == 1 ) then

           !call solver_lumped_mass_system(ndofn,rr,EXCHANGE=.false.)
           !do ipoin = 1,npoin
           !   unkno(ipoin) = unkno(ipoin) + rr(ipoin) * xdiag 
           !end do
           !return
           
           do ipoin = 1,npoin
              unkno(ipoin) = unkno(ipoin) + rr(ipoin) * xdiag / vmass(ipoin)
              !unkno(ipoin) = unkno(ipoin) + rr(ipoin) * 0.1_rp / vmass(ipoin)
           end do

        else if( ndofn == 2 ) then
           do ipoin = 1,npoin
              itotn        = (ipoin-1)*ndofn+1
              fact1        = xdiag / vmass(ipoin)  
              unkno(itotn) = unkno(itotn) + rr(itotn) * fact1
              itotn        = itotn + 1
              unkno(itotn) = unkno(itotn) + rr(itotn) * fact1
           end do

        else if( ndofn == 3 ) then
           do ipoin = 1,npoin
              itotn        = (ipoin-1)*ndofn+1
              fact1        = xdiag / vmass(ipoin)  
              unkno(itotn) = unkno(itotn) + rr(itotn) * fact1
              itotn        = (ipoin-1)*ndofn+2
              unkno(itotn) = unkno(itotn) + rr(itotn) * fact1
              itotn        = (ipoin-1)*ndofn+3
              unkno(itotn) = unkno(itotn) + rr(itotn) * fact1
           end do
 
        else
           do ipoin = 1,npoin
              itotn = (ipoin-1)*ndofn
              fact1 = xdiag / vmass(ipoin) 
              do idofn = 1,ndofn
                 itotn = itotn+1
                 unkno(itotn) = unkno(itotn) + rr(itotn) * fact1
              end do
           end do
        end if

     else if( solve_sol(1) % kfl_preco == SOL_CLOSE_MASS_MATRIX ) then

        !----------------------------------------------------------------
        !
        ! CLOSE MASS MATRIX [CLOSE]:
        !
        ! x^{i+1} = x^i + xdiag.Mc^-1.r^i
        !
        !----------------------------------------------------------------

        if( ndofn == 1 ) then
           do ipoin = 1,npoin
              unkno(ipoin) = unkno(ipoin) + rr(ipoin) * xdiag / vmasc(ipoin)
           end do

        else if( ndofn == 2 ) then
           do ipoin = 1,npoin
              itotn        = (ipoin-1)*ndofn+1
              fact1        = xdiag / vmasc(ipoin)  
              unkno(itotn) = unkno(itotn) + rr(itotn) * fact1
              itotn        = itotn + 1
              unkno(itotn) = unkno(itotn) + rr(itotn) * fact1
           end do

        else if( ndofn == 3 ) then
           do ipoin = 1,npoin
              itotn        = (ipoin-1)*ndofn+1
              fact1        = xdiag / vmasc(ipoin)  
              unkno(itotn) = unkno(itotn) + rr(itotn) * fact1
              itotn        = itotn + 1
              unkno(itotn) = unkno(itotn) + rr(itotn) * fact1
              itotn        = itotn + 1
              unkno(itotn) = unkno(itotn) + rr(itotn) * fact1
           end do

        else
           do ipoin = 1,npoin
              itotn = (ipoin-1)*ndofn
              fact1 = xdiag / vmasc(ipoin)
              do idofn = 1,ndofn
                 itotn = itotn + 1
                 unkno(itotn) = unkno(itotn) + rr(itotn) * fact1
              end do
           end do
        end if

     else if( solve_sol(1) % kfl_preco  == SOL_LOCAL_DIAGONAL ) then

        !----------------------------------------------------------------
        !
        ! LOCAL DIAGONAL [LOCDI]:
        !
        ! x^{i+1} = x^i + P.M^-1.r^i
        !
        !----------------------------------------------------------------

        if( ndofn == 1 ) then
           do ipoin = 1,npoin
              unkno(ipoin) = unkno(ipoin) + rr(ipoin) * pmatr(ipoin) / vmass(ipoin)
           end do

        else if( ndofn == 2 ) then
           do ipoin = 1,npoin
              itotn        = (ipoin-1)*ndofn+1
              fact1        = 1.0_rp / vmass(ipoin)  
              unkno(itotn) = unkno(itotn) + rr(itotn) * pmatr(itotn) * fact1
              itotn        = itotn + 1
              unkno(itotn) = unkno(itotn) + rr(itotn) * pmatr(itotn) * fact1
           end do

        else if( ndofn == 3 ) then
           do ipoin = 1,npoin
              itotn        = (ipoin-1)*ndofn+1
              fact1        = 1.0_rp / vmass(ipoin)  
              unkno(itotn) = unkno(itotn) + rr(itotn) * pmatr(itotn) * fact1
              itotn        = itotn + 1
              unkno(itotn) = unkno(itotn) + rr(itotn) * pmatr(itotn) * fact1
              itotn        = itotn + 1
              unkno(itotn) = unkno(itotn) + rr(itotn) * pmatr(itotn) * fact1
           end do

        else
           do ipoin = 1,npoin
              itotn = (ipoin-1)*ndofn
              fact1 = 1.0_rp / vmass(ipoin)
              do idofn = 1,ndofn
                 itotn        = itotn + 1
                 unkno(itotn) = unkno(itotn) + rr(itotn) * pmatr(itotn) * fact1
              end do
           end do

        end if

     else if( solve_sol(1) % kfl_preco  == SOL_LOCAL_MASS_MATRIX ) then

        !----------------------------------------------------------------
        !
        ! LOCAL MASS MATRIX:
        !
        ! x^{i+1} = x^i + P^-1.r^i
        !
        !----------------------------------------------------------------

        if( ndofn == 1 ) then
           do ipoin = 1,npoin
              fact1        = 1.0_rp / pmatr(ipoin)  
              unkno(ipoin) = unkno(ipoin) + rr(ipoin) * fact1
           end do

        else if( ndofn == 2 ) then
           do ipoin = 1,npoin
              fact1        = 1.0_rp / pmatr(ipoin)  
              itotn        = (ipoin-1)*ndofn+1
              unkno(itotn) = unkno(itotn) + rr(itotn) * fact1
              itotn        = (ipoin-1)*ndofn+2
              unkno(itotn) = unkno(itotn) + rr(itotn) * fact1
           end do

        else if( ndofn == 3 ) then
           do ipoin = 1,npoin
              fact1        = 1.0_rp / pmatr(ipoin)  
              itotn        = (ipoin-1)*ndofn+1
              unkno(itotn) = unkno(itotn) + rr(itotn) * fact1
              itotn        = (ipoin-1)*ndofn+2
              unkno(itotn) = unkno(itotn) + rr(itotn) * fact1
              itotn        = (ipoin-1)*ndofn+3
              unkno(itotn) = unkno(itotn) + rr(itotn) * fact1
           end do

        else

           do ipoin = 1,npoin
              itotn = (ipoin-1)*ndofn
              fact1 = 1.0_rp / pmatr(ipoin)
              do idofn = 1,ndofn
                 itotn        = itotn + 1
                 unkno(itotn) = unkno(itotn) + rr(itotn) * fact1
              end do
           end do

        end if

     else if( solve_sol(1) % kfl_preco == SOL_AUGMENTED_DIAGONAL ) then

        !----------------------------------------------------------------
        !
        ! AUGMENTED DIAGONAL [AUGME]:
        !
        ! x^{i+1} = ( M/xdiag + P )^-1.( M/xdiag. x^i + r^i )
        !
        !----------------------------------------------------------------

        if( ndofn == 1 ) then
           do ipoin = 1,npoin
              fact2        = vmass(ipoin) / xdiag
              fact1        = 1.0_rp / ( fact2 + pmatr(ipoin) )
              unkno(ipoin) = fact1 * ( fact2 * unkno(ipoin) + rr(ipoin) )
           end do

        else if( ndofn == 2 ) then
           do ipoin = 1,npoin
              itotn        = (ipoin-1)*ndofn+1
              fact2        = vmass(ipoin) / xdiag
              fact1        = 1.0_rp / ( fact2 + pmatr(itotn) )
              unkno(itotn) = fact1 * ( fact2 * unkno(itotn) + rr(itotn) )               
              itotn        = itotn + 1
              fact1        = 1.0_rp / ( fact2 + pmatr(itotn) )
              unkno(itotn) = fact1 * ( fact2 * unkno(itotn) + rr(itotn) )               
           end do

        else if( ndofn == 3 ) then
           do ipoin = 1,npoin
              itotn        = (ipoin-1)*ndofn+1
              fact2        = vmass(ipoin) / xdiag
              fact1        = 1.0_rp / ( fact2 + pmatr(itotn) )
              unkno(itotn) = fact1 * ( fact2 * unkno(itotn) + rr(itotn) )               
              itotn        = itotn + 1
              fact1        = 1.0_rp / ( fact2 + pmatr(itotn) )
              unkno(itotn) = fact1 * ( fact2 * unkno(itotn) + rr(itotn) )               
              itotn        = itotn + 1
              fact1        = 1.0_rp / ( fact2 + pmatr(itotn) )
              unkno(itotn) = fact1 * ( fact2 * unkno(itotn) + rr(itotn) )               
           end do

        else
           do ipoin = 1,npoin
              itotn = (ipoin-1)*ndofn
              fact2 = vmass(ipoin) / xdiag
              do idofn = 1,ndofn
                 itotn        = itotn + 1
                 fact1        = 1.0_rp / ( fact2 + pmatr(itotn) )
                 unkno(itotn) = fact1 * ( fact2 * unkno(itotn) + rr(itotn) )
              end do
           end do
        end if

     else

        call runend('SOLRIC: IMPOSSIBLE PRECONDITIONER')

     end if

     solve_sol(1) % iters = solve_sol(1) % iters + 1
     !
     ! Convergence post process
     !
     if( solve_sol(1) % kfl_cvgso == 1 ) &
          call outcso(amatr,rr,unkno)

  end do

  resfi = resi1
  !
  ! Deallocate memory
  !
  if( solve_sol(1) % kfl_algso == SOL_SOLVER_MATRIX_RICHARDSON ) then
     call memory_deallo(memit,'RR','solric',rr)
  end if

100 format(i7,1x,e12.6)
110 format(i5,18(2x,e12.6))

end subroutine solric
