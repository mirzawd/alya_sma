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
!> @file    bsyjes.f90
!> @author  Guillaume Houzeaux
!> @date    19/01/2018
!> @brief   Some operations of DCG
!> @details The following operations of DCG are carried out:
!>
!>          1. w      = A z
!>          2. newrho = r.z
!>          3. mu     = W^T.w
!>          4. All reduce of newrho and mu at the same time
!>
!> @}
!------------------------------------------------------------------------

subroutine bsyjes(nbnodes,nbvar,an,zz,ww,rr,ngrou,newrho,mu)

  use def_kintyp,         only       :  ip,rp
  use def_master,         only       :  INOTMASTER,IPARALL
  use def_master,         only       :  icoml,parre
  use def_solver,         only       :  solve_sol
  use mod_solver,         only       :  solver_parallel_SpMV
  use mod_solver,         only       :  solver_parallel_scalar_product
  use mod_communications, only       :  PAR_SUM
  implicit none
  integer(ip), intent(in)            :: nbnodes
  integer(ip), intent(in)            :: nbvar
  integer(ip), intent(in)            :: ngrou
  real(rp),    intent(in)            :: an(nbvar,nbvar,*)
  real(rp),    intent(inout)         :: zz(nbvar*nbnodes)
  real(rp),    intent(inout)         :: rr(nbvar*nbnodes)
  real(rp),    intent(inout), target :: ww(nbvar,nbnodes)
  real(rp),    intent(out),   target :: mu(nbvar*ngrou+1)
  real(rp),    intent(out)           :: newrho
  integer(ip)                        :: ipoin,igrou,ngrou1
  integer(ip)                        :: ii,kk

  !----------------------------------------------------------------
  !
  ! w      = A z
  ! newrho = r.z
  !
  !----------------------------------------------------------------

  call solver_parallel_SpMV(solve_sol(1),an,zz,ww,MPI=.false.,OPENMP=.true.,TIMING=.false.)
  call solver_parallel_scalar_product(solve_sol(1),rr,zz,newrho,MPI=.false.)

  !----------------------------------------------------------------
  !
  ! mu = W^T.w
  !
  !----------------------------------------------------------------

  do igrou = 1,solve_sol(1)%ngrou*nbvar
     mu(igrou) = 0.0_rp
  end do

  if( INOTMASTER ) then

     if( nbvar == 1 )then
        do ipoin = 1,nbnodes
           igrou = solve_sol(1) % lgrou(ipoin)
           if( igrou > 0 ) mu(igrou) = mu(igrou) + ww(1,ipoin)
        end do
     else
        do ipoin = 1,nbnodes
           igrou = solve_sol(1) % lgrou(ipoin)
           if( igrou > 0 ) then
              ii = (igrou-1)*nbvar
              do kk = 1,nbvar
                 ii = ii + 1
                 mu(ii) = mu(ii) + ww(kk,ipoin)
              end do
           end if
        end do
     end if

  end if

  !----------------------------------------------------------------
  !
  ! All reduce rho and mu
  !
  !----------------------------------------------------------------

  if( IPARALL ) then

     if( solve_sol(1) % kfl_gathe == 1 ) then
        !
        ! All reduce for rho and all gatherv for mu
        !
        call PAR_SUM(newrho)
        parre        => mu
        icoml        =  solve_sol(1)%icoml
        call par_gatsen()

     else if( solve_sol(1) % kfl_gathe == 0 ) then
        !
        ! All reduce for mu and rho
        !
        ngrou1     = solve_sol(1)%ngrou*nbvar+1
        mu(ngrou1) = newrho
        call PAR_SUM(ngrou1,mu)
        newrho = mu(ngrou1)

     else if( solve_sol(1) % kfl_gathe == 2 ) then
        !
        ! All reduce for rho and send/receive for mu
        !
        call PAR_SUM(newrho)
        parre        => mu
        call par_cregrp(3_ip) 

     end if

  end if


end subroutine bsyjes
