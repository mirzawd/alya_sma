!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!>
!> @defgroup Autotuning
!> @{
!> @name    Subroutine for auto-tuning of some algorithms
!> @file    mod_auto_tuning.f90
!> @author  Guillaume Houzeaux
!> @date    29/11/2017
!> @brief   ToolBox for autotuning
!> @details ToolBox for autotuning
!>
!------------------------------------------------------------------------

module mod_auto_tuning

  use def_kintyp,   only : ip,rp
  use def_kintyp_solvers,   only : soltyp
  use mod_outfor,   only : outfor
  use mod_messages, only : livinf
  use mod_messages, only : messages_live

  implicit none
  private

  public :: auto_tuning_SpMV_OpenMP

contains

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    29/11/2017     
  !> @brief   Auto-tuning for SpMV
  !> @details Look for best chunk size
  !>
  !-----------------------------------------------------------------------

  subroutine auto_tuning_SpMV_OpenMP()

    use def_kintyp,         only : ip,rp,lg
    use def_master,         only : INOTMASTER,solve_sol
    use def_master,         only : amatr,rhsid,unkno,routp,coutp
    use def_master,         only : mmodu,momod,modul,kfl_modul,lun_outpu
    use def_master,         only : ioutp,coutp
    use mod_communications, only : PAR_SUM
    use mod_communications, only : PAR_MIN
    use mod_communications, only : PAR_MAX
    use mod_communications, only : PAR_AVERAGE
    use mod_communications, only : PAR_BARRIER
    use mod_solver,         only : solver_parallel_SpMV
    use mod_solver,         only : solver_initialization
    use mod_solver,         only : solver_initialize_matrix_and_rhs
    use def_solver,         only : SOL_OMP_DYNAMIC
    use def_solver,         only : SOL_OMP_STATIC
    use def_solver,         only : SOL_OMP_GUIDED
    use def_solver,         only : SOL_OMP_OFF
    use def_solver,         only : memit
    use mod_parall,         only : par_hybrid
    use mod_parall,         only : par_omp_num_threads
    use mod_parall,         only : PAR_HYBRID_OFF
    use mod_memory,         only : memory_alloca
    use mod_memory,         only : memory_deallo
    use mod_memory,         only : memory_size
    use mod_maths,          only : maths_time_unit

    integer(ip),  parameter     :: num_tests = 100
    integer(ip)                 :: ii,ii_min,kk,ivari
    real(rp)                    :: times(-3:num_tests)
    integer(ip)                 :: sizes(-3:num_tests)
    real(rp)                    :: time1,time2,delta
    integer(ip)                 :: chunk_size
    integer(ip)                 :: chunk_size_opt
    integer(ip)                 :: chunk_size_min
    integer(ip)                 :: chunk_size_max
    integer(ip)                 :: chunk_size_ave
    real(rp)                    :: xfact
    real(rp)                    :: time_ave
    real(rp)                    :: time_min
    real(rp)                    :: time_max
    integer(ip)                 :: chunk_size_static
    integer(ip)                 :: chunk_step
    integer(ip)                 :: max_stats
    integer(ip)                 :: nunkn,nzmat,nzrhs,ndofn
    integer(ip)                 :: nrows,ncols,ipass
    integer(ip)                 :: num_scheduling(-3:1)
    real(rp),   pointer         :: xx(:)
    logical(lg)                 :: auto_schedule
    !
    ! Hey! This makes sense only if you use OpenMP!
    !
    if( par_hybrid == PAR_HYBRID_OFF ) return
    if( par_omp_num_threads <= 0 )     return

    ipass = 0
    max_stats = 10
    nullify(xx)

    do modul = 1,mmodu
       if( kfl_modul(modul) == 1 ) then
          if( associated(momod(modul) % solve) ) then
             do ivari = 1,size(momod(modul) % solve,KIND=ip)
                solve_sol => momod(modul) % solve(ivari:)
                if(    solve_sol(1) % kfl_algso >  0 .and.   &
                     & solve_sol(1) % kfl_algso /= 9 .and.   &
                     & solve_sol(1) % omp_chunk_size <= -1 ) then

                   if( ipass == 0 ) then
                      call messages_live('AUTO-TUNING OF OPENMP FOR SpMV','START SECTION')
                      ipass = 1
                      call outfor(87_ip,lun_outpu,' ')
                   end if
                   call messages_live('START AUTO-TUNING FOR '//trim(solve_sol(1) % wprob),'START SECTION')

                   coutp(1) = trim(solve_sol(1) % wprob)
                   call outfor(88_ip,lun_outpu,' ')
                   auto_schedule    = .false.
                   if( solve_sol(1) % omp_chunk_size == -2 ) auto_schedule = .true.
                   
                   times          = 0.0_rp
                   chunk_size_opt = 0
                   !
                   ! Solver array sizes
                   !
                   ndofn = solve_sol(1) % ndofn
                   nrows = solve_sol(1) % nunkn
                   ncols = solve_sol(1) % ncols
                   nzrhs = solve_sol(1) % nzrhs
                   nzmat = solve_sol(1) % nzmat
                   nunkn = max(nrows,ncols)

                   if( INOTMASTER ) then
                      !
                      ! Memory
                      !
                      if( solve_sol(1) % kfl_full_rows == 0 ) then
                         xx => unkno
                      else
                         call memory_alloca(memit,'XX','auto_tuning_SpMV_OpenMP',xx,nunkn)
                      end if
                      if( nzmat < solve_sol(1) % nnz ) call runend('MARDE MARDE 1')
                      if( nzmat > memory_size(amatr) ) call runend('MARDE MARDE 2')
                      if( nunkn > memory_size(xx) )    call runend('MARDE MARDE 3')
                      !
                      ! Initialization
                      !
                      do ii = 1,nzrhs
                         rhsid(ii) = 1.0_rp
                      end do
                      do ii = 1,nunkn
                         xx(ii) = 1.0e-3_rp                         
                      end do
                      do ii = 1,nzmat
                         amatr(ii) = 1.0_rp
                      end do
                      !
                      ! Determine numebr of tests and chunk sizes
                      !
                      chunk_size_static = solve_sol(1) % nequa / max(par_omp_num_threads,1_ip)
                      chunk_size_max    = chunk_size_static
                      chunk_size_min    = max(chunk_size_static/500,num_tests)
                      chunk_step        = ( chunk_size_max - chunk_size_min ) / num_tests
                      delta             = real(chunk_size_max-chunk_size_min,rp) / (real(num_tests,rp)-1.0_rp)
                      times             = huge(1.0_rp)
                      !
                      ! DYNAMIC scheduling
                      !
                      ! Try NUM_TESTS different chunks
                      ! Average over MAX_STATS SpMV and do not take into account first shot, just in case
                      !                      
                      chunk_size = 0
                      do ii = 1,num_tests
                         chunk_size = chunk_size_min + int(delta * (real(ii,rp)-1.0_rp),ip)
                         solve_sol(1) % omp_chunk_size = chunk_size                         
                         solve_sol(1) % omp_schedule   = SOL_OMP_DYNAMIC
                         times(ii) = 0.0_rp                         
                         do kk = 1,max_stats
                            call PAR_BARRIER('IN MY CODE WITHOUT MASTER')
                            call cputim(time1)
                            call solver_parallel_SpMV(&
                                 solve_sol(1),amatr,xx,rhsid,&
                                 OPENMP=.true.)
                            call cputim(time2)
                            if( kk /= 1 ) times(ii) = times(ii) + time2-time1
                         end do
                         times(ii) = times(ii) / real(max_stats-1,rp) 
                         sizes(ii) = chunk_size
                      end do

                      if( auto_schedule ) then
                         !
                         ! STATIC scheduling
                         !
                         solve_sol(1) % omp_schedule = SOL_OMP_STATIC
                         times(-2) = 0.0_rp                         
                         do kk = 1,max_stats
                            call PAR_BARRIER('IN MY CODE WITHOUT MASTER')
                            call cputim(time1)
                            call solver_parallel_SpMV(&
                                 solve_sol(1),amatr,xx,rhsid,&
                                 OPENMP=.true.)
                            call cputim(time2)
                            if( kk /= 1 ) times(-2) = times(-2) + time2-time1
                         end do
                         times(-2) = times(-2) / real(max_stats-1,rp) 
                         sizes(-2) = 0    
                         !
                         ! GUIDED scheduling
                         !
                         solve_sol(1) % omp_schedule = SOL_OMP_GUIDED
                         times(-3) = 0.0_rp                         
                         do kk = 1,max_stats
                            call PAR_BARRIER('IN MY CODE WITHOUT MASTER')
                            call cputim(time1)
                            call solver_parallel_SpMV(&
                                 solve_sol(1),amatr,xx,rhsid,&
                                 OPENMP=.true.)
                            call cputim(time2)
                            if( kk /= 1 ) times(-3) = times(-3) + time2-time1
                         end do
                         times(-3) = times(-3) / real(max_stats-1,rp) 
                         sizes(-3) = 0                      
                         !
                         ! NO OPENMP 
                         !
                         solve_sol(1) % omp_schedule = SOL_OMP_OFF
                         times(0) = 0.0_rp                         
                         do kk = 1,max_stats
                            call PAR_BARRIER('IN MY CODE WITHOUT MASTER')
                            call cputim(time1)
                            call solver_parallel_SpMV(&
                                 solve_sol(1),amatr,xx,rhsid,&
                                 OPENMP=.false.)
                            call cputim(time2)
                            if( kk /= 1 ) times(0) = times(0) + time2-time1
                         end do
                         times(0) = times(0) / real(max_stats-1,rp) 
                         sizes(0) = 0
                      end if
                      !
                      ! Deallocate
                      !
                      if( solve_sol(1) % kfl_full_rows /= 0 ) then
                         call memory_deallo(memit,'XX','auto_tuning_SpMV_OpenMP',xx)
                      end if

                   end if
                   !
                   ! Scheduling and optimum chunk size 
                   ! MINLOC does not work here!
                   !
                   time_min = huge(1.0_rp)
                   do ii = -3,num_tests
                      if( times(ii) < time_min ) then
                         ii_min   = ii
                         time_min = times(ii_min)
                      end if
                   end do

                   num_scheduling = 0

                   if(      ii_min == -3 ) then
                      num_scheduling(-3)            = 1
                      solve_sol(1) % omp_schedule   = SOL_OMP_GUIDED
                      solve_sol(1) % omp_chunk_size = 0
                      chunk_size_min                =  huge(1_ip)
                      chunk_size_max                = -huge(1_ip)
                      chunk_size_ave                = 0
                   else if( ii_min == -2 ) then
                      num_scheduling(-2)            = 1
                      solve_sol(1) % omp_schedule   = SOL_OMP_STATIC
                      solve_sol(1) % omp_chunk_size = 0
                      chunk_size_min                =  huge(1_ip)
                      chunk_size_max                = -huge(1_ip)
                      chunk_size_ave                = 0
                   else if( ii_min ==  0 ) then
                      num_scheduling(0)             = 1
                      solve_sol(1) % omp_schedule   = SOL_OMP_OFF
                      solve_sol(1) % omp_chunk_size = 0
                      chunk_size_min                =  huge(1_ip)
                      chunk_size_max                = -huge(1_ip)
                      chunk_size_ave                = 0
                   else if( ii_min >=  1 ) then
                      num_scheduling(1)             = 1
                      chunk_size_opt                = sizes(ii_min)
                      solve_sol(1) % omp_schedule   = SOL_OMP_DYNAMIC
                      solve_sol(1) % omp_chunk_size = chunk_size_opt
                      chunk_size_min                = chunk_size_opt
                      chunk_size_max                = chunk_size_opt
                      chunk_size_ave                = chunk_size_opt
                   else
                      call runend('auto_tuning_SpMV_OpenMP: UNKNOWN BEST SCHEDULING')
                   end if
                   !
                   ! Output
                   !
                   time_min = times(ii_min)
                   time_max = times(ii_min)
                   time_ave = times(ii_min)

                   call PAR_SUM(5_ip,num_scheduling)

                   call PAR_MIN(chunk_size_min)
                   call PAR_MAX(chunk_size_max)
                   call PAR_SUM(chunk_size_ave)

                   call PAR_MIN(time_min)
                   call PAR_MAX(time_max)
                   call PAR_AVERAGE(time_ave)

                   if( num_scheduling(-3) > 0 ) then
                      call maths_time_unit(time_max,coutp(1),xfact)
                      call livinf(-14_ip,'GUIDED SCHEDULING=              ',num_scheduling(-3))                      
                   end if
                   if( num_scheduling(-2) > 0 ) then
                      call maths_time_unit(time_max,coutp(1),xfact)
                      call livinf(-14_ip,'STATIC SCHEDULING=              ',num_scheduling(-2))   
                   end if
                   if( num_scheduling( 0) > 0 ) then
                      call maths_time_unit(time_max,coutp(1),xfact)
                      call livinf(-14_ip,'OPENMP DESACTIVATED=            ',num_scheduling( 0))                      
                   end if
                   if( num_scheduling(1) > 0 ) then
                      chunk_size_ave = chunk_size_ave / num_scheduling(1)
                      call livinf(-14_ip,'DYNAMIC SCHEDULING=             ',num_scheduling(1))
                      call livinf(-14_ip,'  - MINIMUM OPTIMUM CHUNK SIZE= ',chunk_size_min)
                      call livinf(-14_ip,'  - MAXIMUM OPTIMUM CHUNK SIZE= ',chunk_size_max)
                      call livinf(-14_ip,'  - AVERAGE OPTIMUM CHUNK SIZE= ',chunk_size_ave)
                      ioutp(3) = chunk_size_ave
                      ioutp(4) = chunk_size_min
                      ioutp(5) = chunk_size_max
                   end if

                   ioutp(10) = num_scheduling(-3) ! Guided
                   ioutp(11) = num_scheduling(-2) ! Static
                   ioutp(12) = num_scheduling( 0) ! Off
                   ioutp(13) = num_scheduling( 1) ! Dynamic

                   call maths_time_unit(time_max,coutp(1),xfact) 
                   routp(3) = time_ave * xfact
                   routp(4) = time_min * xfact
                   routp(5) = time_max * xfact
                   routp(6) = routp(3) / (routp(5)+epsilon(1.0_rp))

                   call outfor(89_ip,lun_outpu,' ')
                   call messages_live('AUTO-TUNING','END SECTION')

                end if

             end do
          end if
       end if
    end do

    if( ipass == 1 ) call messages_live('AUTO-TUNING OF OPENMP FOR SpMV','END SECTION') 

  end subroutine auto_tuning_SpMV_OpenMP

end module mod_auto_tuning
!> @}

