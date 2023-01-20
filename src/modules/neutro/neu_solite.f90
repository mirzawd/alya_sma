!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Neutro 
!> @{
!> @file    neu_solite.f90
!> @date    29/03/2016
!> @author  Guillaume Houzeaux
!> @brief   Solve an inner iteration 
!> @details Solve an inner iteration for current energy and direction
!> @}
!------------------------------------------------------------------------

subroutine neu_solite()

  use def_parame 
  use def_master
  use def_domain
  use def_neutro
  use mod_solver, only : solver_postprocess
  use mod_solver, only : solver_solve
  use mod_ADR,    only : ELEMENT_ASSEMBLY
  implicit none
  integer(ip) :: ipoin
!   real(rp)    :: invb2
  external :: inisol, neu_updunk, neu_elmope, neu_bouope
  !
  ! Initialize solver
  !
  call inisol() 
  !momod(modul) % solve % kfl_fixno => kfl_fixno_neu(current_energy_neu,current_direction_neu) % l(:,:)
  !
  ! Fill in advection (e.g. can be useful for iteratove solvers)
  !
  if( associated(advec) ) then
     do ipoin = 1,npoin
        advec(1:ndime,ipoin,1) = direc_neu(1:ndime,current_direction_neu)
     end do
  end if
  !
  ! Obtain the initial guess for UNKNO
  !
  call neu_updunk(6_ip) !> Option 6 is initial guess for solver
  ! 
  ! Assemble equations 
  !
  if( INOTMASTER ) then
     call neu_elmope(ELEMENT_ASSEMBLY) !> We assemble the elementary matrix into the system matrix
     call neu_bouope() !> We assemble the boundary element into the system
  end if
  !
  ! Solve equations 
  !
  call solver_solve(momod(modul) % solve,amatr,rhsid,unkno,pmatr) !> Solve the system of equations
 
  !
  ! Solver postprocess
  !
  call solver_postprocess(momod(modul) % solve,amatr,rhsid,unkno) !> Post-process the outcome


end subroutine neu_solite
