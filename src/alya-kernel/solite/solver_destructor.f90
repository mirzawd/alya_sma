!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solver
!> @{
!> @file    solver_destructor.f90
!> @author  houzeaux
!> @date    2019-06-11
!> @brief   Destroy solver
!> @details Reinitialize all solver arrays
!> @} 
!-----------------------------------------------------------------------

subroutine solver_destructor()

  use def_kintyp, only : ip
  use def_master, only : mmodu
  use def_master, only : momod
  use def_master, only : kfl_modul
  use mod_solver, only : solver_deallocate
  use def_solver, only : SOL_NO_SOLVER
  implicit none

  integer(ip) :: ivari,imodu
  
  do imodu = 1,mmodu
     if( kfl_modul(imodu) == 1 ) then
        if( associated(momod(imodu) % solve) ) then
           do ivari = 1,size(momod(imodu) % solve,KIND=ip)
              if( momod(imodu) % solve(ivari) % kfl_algso /= SOL_NO_SOLVER ) then
                 call solver_deallocate(momod(imodu) % solve(ivari),ivari)
                 if( momod(imodu) % solve(ivari) % kfl_clean_precond == 0 ) then
                    momod(imodu) % solve(ivari) % kfl_clean_precond = -1
                 end if
              end if
           end do
        end if
     end if
  end do

end subroutine solver_destructor
