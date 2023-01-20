!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kernel
!> @{
!> @file    what_is_being_done.f90
!> @author  houzeaux
!> @date    2020-03-04
!> @brief   Do things
!> @details According to what is being done, change some flags
!> @} 
!-----------------------------------------------------------------------

subroutine what_is_being_done()

  use def_master
  use def_solver
  implicit none
  integer(ip) :: imodu,ivari
  !
  ! Force assembly os solvers if we have just read a restart or did a repartitioning
  ! Sometime matrices are only assembled when ittim==1 !
  !
  do imodu = 1,mmodu
     if( kfl_modul(imodu) == 1 ) then
        if( associated(momod(imodu) % solve) ) then
           do ivari = 1,size(momod(imodu) % solve,KIND=ip) 
              if( momod(imodu) % solve(ivari) % kfl_algso /= SOL_NO_SOLVER ) then
                 if( do_repartitioning .or. do_amr .or. read_restart ) then
                    momod(imodu) % solve(ivari) % kfl_force_assembly = 1
                 else
                    momod(imodu) % solve(ivari) % kfl_force_assembly = 0
                 end if
              end if
           end do
        end if
     end if
  end do
  do_repartitioning = .false.
  do_amr            = .false.
  read_restart      = .false.
  write_restart     = .false.

end subroutine what_is_being_done
  
  
