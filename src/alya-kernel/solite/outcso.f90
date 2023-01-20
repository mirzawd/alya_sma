!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Algebraic_Solver
!> @{                                                                   
!> @file    outcso.f90
!> @author  Guillaume Houzeaux
!> @date    15/07/2015
!> @brief   Output solve convergence
!> @details Output solve convergence and timing
!>          Preconditioned residual: || M^-1 (b - Ax) || / || M^-1 b ||
!>          Real residual            || b - Ax || / || b ||
!>          KFL_EXRES = 0 ... Preconditioned residual 
!>                    = 1 ... Preconditioned residual and real residual
!>
!> @}                                                                   
!-----------------------------------------------------------------------

subroutine outcso(an,bb,xx)

  use def_kintyp, only       :  ip,rp
  use def_master, only       :  INOTSLAVE
  use def_solver, only       :  resi1,solve_sol
  use def_solver, only       :  SOL_SOLVER_RICHARDSON
  use def_solver, only       :  SOL_SOLVER_MATRIX_RICHARDSON
  use def_solver, only       :  SOL_RIGHT_PRECONDITIONING
  use mod_solver, only       :  solver_parallel_residual_and_norm
  implicit none
  real(rp),    intent(in)    :: an(*)
  real(rp),    intent(in)    :: bb(*)
  real(rp),    intent(inout) :: xx(*)
  real(rp)                   :: rnorm
  real(rp),    save          :: time1,time2,timei
  real(rp)                   :: tdiff 
 
  if( solve_sol(1) % kfl_cvgso == 1 ) then
     
     if( solve_sol(1) % kfl_exres /= 0 ) then
        !
        ! If required, compute RNORM = || b - Ax || / || b || 
        !
        if(  solve_sol(1) % kfl_algso == SOL_SOLVER_RICHARDSON .or. &
             solve_sol(1) % kfl_algso == SOL_SOLVER_MATRIX_RICHARDSON ) then
           call runend('OUTCSO: OPTION IMPOSSIBLE')
        else           
           if ( solve_sol(1) % kfl_leftr == SOL_RIGHT_PRECONDITIONING ) then
              !
              ! Right preconditioning
              !
              rnorm = resi1
           else
              !
              ! Left preconditioning
              !
              call solver_parallel_residual_and_norm(solve_sol(1),an,xx,bb,RESIDUAL_NORM=rnorm)
              if( solve_sol(1) % bnorm /= 0.0_rp ) rnorm = rnorm / solve_sol(1) % bnorm
           end if
           
        end if
     else
        rnorm = 0.0_rp         
     end if

     call cputim(time2)
     if( solve_sol(1) % iters >= 1 ) then
        tdiff = time2 - time1
     else
        timei = time2
        tdiff = 0.0_rp
     end if
     time1 = time2

     if( INOTSLAVE ) then
        write(solve_sol(1) % lun_cvgso,100) &
             solve_sol(1) % iters,resi1,rnorm,tdiff,time2-timei,&
             solve_sol(1) % xorth
     end if
     
  end if

100 format(i7,20(1x,e16.8E3))

end subroutine outcso

subroutine outress(npoin,xx)
  !------------------------------------------------------------------------
  !****f* solite/outcso
  ! NAME 
  !    outcso
  ! DESCRIPTION
  !    Exact residual
  ! USES
  ! USED BY 
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_master, only       :  intost
  use def_domain, only       :  coord
  use def_solver, only       :  solve_sol
  implicit none
  integer(ip), intent(in)    :: npoin
  real(rp),    intent(in)    :: xx(*)
  integer(ip)                :: ii
  real(rp)                   :: x
  
  if(solve_sol(1) % iters==5.or.solve_sol(1) % iters==10.or.solve_sol(1) % iters==50) then
     open(unit=90,file='dcg-'//trim(intost(solve_sol(1) % iters))//'.txt',status='unknown')
     !call memgen(0_ip,npoin,0_ip)
     do ii = 1,npoin
        x=coord(1,ii)
        write(90,*) x,xx(ii)-(0.5_rp*x*(1.0_rp-x))
     end do
     close(90)
     !wopos(1) = 'RESID'
     !wopos(2) = 'SCALA'
     !wopos(3) = 'NPOIN'     
     !call postpr(gesca,wopos,solve_sol(1) % iters,real(ittim,rp))   
     !call memgen(2_ip,npoin,0_ip)
  end if

end subroutine outress
