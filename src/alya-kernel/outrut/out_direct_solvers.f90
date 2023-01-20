!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kernel
!> @{
!> @file    out_direct_solvers.f90
!> @author  Guillaume Houzeaux
!> @brief   Direct solver statistic
!> @details Output some direct solver statistics
!> @} 
!-----------------------------------------------------------------------

subroutine out_direct_solvers()

  use def_kintyp,         only : ip,rp,direct_solver_typ
  use def_master,         only : kfl_modul,mmodu,momod,IMASTER
  use def_master,         only : namod,ioutp,routp,coutp,intost
  use def_solver,         only : solve_sol,memdi
  use mod_direct_solver,  only : direct_solver_statistics 
  use mod_direct_solver,  only : direct_solver_name
  use mod_memory,         only : memory_unit
  use mod_outfor,         only : outfor
  use mod_communications, only : PAR_SUM
  implicit none
  integer(ip)     :: imodu,ivari
  integer(ip)     :: num_factorizations,num_solutions,num_threads
  integer(ip)     :: num_initializations,ipass,kfl_solver,ii
  integer(8)      :: memor_max,memor_ave
  real(rp)        :: cputi_max(3),cputi_ave(3),load_balance(3)
  real(rp)        :: memor_factor
  character(6)    :: memor_char 
  character(50)   :: wsolver,wname
  real(rp)        :: solver_statistics(20)
  type(direct_solver_typ) :: solve_dummy

  call outfor(74_ip,0_ip,' ')

  do imodu = 1,mmodu
     if( kfl_modul(imodu) == 1 ) then
        if( associated(momod(imodu) % solve) ) then
           do ivari = 1,size(momod(imodu) % solve,KIND=ip)
              solve_sol => momod(imodu) % solve(ivari:)
              if( solve_sol(1) % kfl_algso /= -999 ) then
                 ipass = 1
                 kfl_solver = 0
                 ii = 0
                 do while( ipass <= 6 )
                    if(      ipass == 1 ) then
                       wsolver       =  'DIRECT SOLVER'
                       kfl_solver    =  solve_sol(1) % direct_solver % kfl_solver
                       wname         =  solve_sol(1) % direct_solver % name
                       num_threads   =  solve_sol(1) % direct_solver % num_threads
                       if( kfl_solver /= 0 ) call direct_solver_statistics(solve_sol(1) % direct_solver,solver_statistics)
                       ipass = 2
                    else if( ipass == 2 ) then
                       wsolver       =  'DEFLATION COARSE SOLVER'                       
                       kfl_solver    =  solve_sol(1) % direct_solver_deflation % kfl_solver
                       wname         =  solve_sol(1) % direct_solver_deflation % name
                       num_threads   =  solve_sol(1) % direct_solver_deflation % num_threads
                       if( kfl_solver /= 0 ) call direct_solver_statistics(solve_sol(1) % direct_solver_deflation,solver_statistics)
                       ipass = 3 
                    else if( ipass == 3 ) then
                       wsolver       =  'COARSE GRID CORRECTION'                       
                       kfl_solver    =  solve_sol(1) % direct_solver_coarse % kfl_solver
                       wname         =  solve_sol(1) % direct_solver_coarse % name
                       num_threads   =  solve_sol(1) % direct_solver_coarse % num_threads
                       if( kfl_solver /= 0 ) call direct_solver_statistics(solve_sol(1) % direct_solver_coarse,solver_statistics)
                       ipass = 4
                    else if( ipass == 4 ) then
                       wsolver       =  'BLOCK LU PRECONDITIONER'                       
                       kfl_solver    =  solve_sol(1) % direct_solver_Block_LU % kfl_solver
                       wname         =  solve_sol(1) % direct_solver_Block_LU % name
                       num_threads   =  solve_sol(1) % direct_solver_Block_LU % num_threads
                       if( kfl_solver /= 0 ) call direct_solver_statistics(solve_sol(1) % direct_solver_Block_LU,solver_statistics)
                       ipass = 5
                    else if( ipass == 5 ) then
                       ii = ii + 1
                       if( associated(solve_sol(1) % direct_solver_RAS) ) then
                          if( solve_sol(1) % kfl_block_ras == 0 ) then
                             wsolver =  'RAS PRECONDITIONER'
                          else
                             wsolver = 'BLOCK RAS PRECONDITIONER: BLOCK '//trim(intost(ii))
                          end if
                          kfl_solver    =  solve_sol(1) % direct_solver_RAS(ii) % kfl_solver
                          wname         =  solve_sol(1) % direct_solver_RAS(ii) % name
                          num_threads   =  solve_sol(1) % direct_solver_RAS(ii) % num_threads
                       end if
                       if( kfl_solver /= 0 )then
                          if( associated(solve_sol(1) % direct_solver_RAS) )then
                             call direct_solver_statistics(solve_sol(1) % direct_solver_RAS(ii),solver_statistics)
                          else
                             call direct_solver_statistics(solve_dummy ,solver_statistics)
                          endif
                       endif
                       ipass = 6
                    else if( ipass == 6 ) then
                       wsolver       =  'MULTIGRID PRECONDITIONER'                       
                       kfl_solver    =  solve_sol(1) % direct_solver_AMG % kfl_solver
                       wname         =  solve_sol(1) % direct_solver_AMG % name
                       num_threads   =  solve_sol(1) % direct_solver_AMG % num_threads
                       if( kfl_solver /= 0 ) call direct_solver_statistics(solve_sol(1) % direct_solver_AMG,solver_statistics)
                       ipass = 7
                    end if 
                    if( kfl_solver /= 0 .and. IMASTER ) then
                       cputi_max(1:3)      = solver_statistics(1:3) 
                       cputi_ave(1:3)      = solver_statistics(4:6) 
                       load_balance(1:3)   = solver_statistics(7:9)
                       num_initializations = int(solver_statistics(10),ip)
                       num_factorizations  = int(solver_statistics(11),ip) 
                       num_solutions       = int(solver_statistics(12),ip)  
                       memor_max           = int(solver_statistics(13),8)
                       memor_ave           = int(solver_statistics(14),8)
                       memdi(1:2)          = memdi(1:2) + int(solver_statistics(15),8)
                       call memory_unit(memor_max,memor_char,memor_factor)
                       ioutp(1)   = num_initializations
                       ioutp(2)   = num_factorizations
                       ioutp(3)   = num_solutions
                       ioutp(4)   = num_threads
                       routp(1:3) = cputi_ave(1:3)
                       routp(4:6) = cputi_max(1:3)
                       routp(7:9) = load_balance(1:3)
                       routp(12)  = solver_statistics(16)
                       routp(13)  = solver_statistics(17)
                       routp(10)  = real(memor_ave,rp)*memor_factor 
                       routp(11)  = real(memor_max,rp)*memor_factor
                       coutp(1)   = trim(namod(imodu))
                       coutp(2)   = trim(solve_sol(1) % wprob)
                       coutp(3)   = trim(wsolver)
                       coutp(4)   = trim(wname)
                       coutp(5)   = memor_char
                       call outfor(75_ip,0_ip,' ')
                    end if
                 end do
              end if
           end do
        end if
     end if
  end do

end subroutine out_direct_solvers
