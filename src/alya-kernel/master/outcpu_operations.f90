!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup CPU_Time
!> @{
!> @file    outcpu_operations.f90
!> @author  Guillaume Houzeaux
!> @date    22/11/2017
!> @brief   Output some CPU info
!> @details CPU info about some operations with and without OpenMP
!>          1. Matrix-vector product
!>          2. Scalar product
!> @}
!-----------------------------------------------------------------------

subroutine outcpu_operations(meshe)
  use def_kintyp,         only : ip,rp,lg
  use def_kintyp_solvers,         only : soltyp
  use def_master,         only : INOTMASTER
  use def_master,         only : INOTEMPTY
  use def_master,         only : amatr,rhsid,unkno,routp,coutp
  use def_domain,         only : mesh_type
  use mod_maths,          only : maths_time_unit
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_AVERAGE
  use mod_communications, only : PAR_BARRIER
#ifdef NINJA
  use def_kermod,         only : kfl_coo
#endif
  use def_kermod,         only : kfl_algebra_operations
  use mod_solver,         only : solver_parallel_SpMV
  use mod_solver,         only : solver_parallel_scalar_product
  use mod_solver,         only : solver_initialization
  use mod_solver,         only : solver_initialize_matrix_and_rhs
  use def_solver,         only : SOL_CSR_FORMAT
  use def_solver,         only : SOL_OMP_STATIC
  use mod_outfor,         only : outfor
  use mod_memory,         only : memory_size
  implicit none
  type(mesh_type), intent(in) :: meshe          !< Mesh type
  integer(ip)                 :: ii,jj
  integer(ip)                 :: ipass
  real(rp)                    :: my_times(18)
  real(rp)                    :: times_max(18)
  real(rp)                    :: times_ave(18)
  real(rp)                    :: load_balance(18)
  real(rp)                    :: alpha,xfact1,xfact2,xfact3,xfact4
  real(rp)                    :: my_timing(3)
  type(soltyp)                :: solve

  if(  INOTEMPTY                           .and. &
       kfl_algebra_operations > 0          .and. &
       memory_size(amatr) >= meshe % nzdom .and. &
       memory_size(rhsid) >= meshe % npoin .and. &
       memory_size(unkno) >= meshe % npoin ) then
     ipass = 1
  else
     ipass = 0
  end if

  call PAR_SUM(ipass,'IN MY CODE')

  my_times     = 0.0_rp
  times_max    = 0.0_rp
  load_balance = 0.0_rp
  my_timing    = 0.0_rp
  
  if( ipass > 0 ) then
     !
     ! Initialization
     !
     call solver_initialization(solve)
     solve % nzmat        =  meshe % nzdom
     solve % nequa        =  meshe % npoin
     solve % ndofn        =  1_ip
     solve % nzrhs        =  meshe % npoin * solve % ndofn
     solve % nequa_own    =  meshe % npoin_own
     solve % ia           => meshe % r_dom
     solve % ja           => meshe % c_dom
     solve % omp_schedule =  SOL_OMP_STATIC
     solve % kfl_format   =  SOL_CSR_FORMAT

     if( INOTMASTER ) then
        alpha = 0.0_rp
        do ii = 1,min(memory_size(rhsid),solve % nzrhs)
           rhsid(ii) = 1.0_rp
        end do
        do ii = 1,min(memory_size(unkno),solve % nzrhs)
           unkno(ii) = 1.0e-3_rp
        end do
        do ii = 1,min(memory_size(amatr),solve % nzmat)
           amatr(ii) = 1.0_rp
        end do
     end if
     !
     ! Matrix-vector product without OpenMP
     !
     if( INOTMASTER ) then
        do ipass = 1,kfl_algebra_operations
           call PAR_BARRIER('IN MY CODE WITHOUT MASTER')
           call solver_parallel_SpMV(&
                solve,amatr,unkno,rhsid,&
                OPENMP=.false.,MY_TIMING=my_timing)
           times_max(1) = times_max(1) + my_timing(1)
           times_max(2) = times_max(2) + my_timing(2)
           times_max(3) = times_max(3) + my_timing(1) + my_timing(2)
        end do
     end if
     !
     ! Matrix-vector product with OpenMP
     !
     if( INOTMASTER ) then
        do ipass = 1,kfl_algebra_operations
           call PAR_BARRIER('IN MY CODE WITHOUT MASTER')
           call solver_parallel_SpMV(&
                solve,amatr,unkno,rhsid,&
                OPENMP=.true.,MY_TIMING=my_timing)
           times_max(4) = times_max(4) + my_timing(1)
           times_max(5) = times_max(5) + my_timing(2)
           times_max(6) = times_max(6) + my_timing(1) + my_timing(2)
        end do
     end if
     !
     ! Scalar-product without OpenMP
     !
     do ipass = 1,kfl_algebra_operations
        call solver_parallel_scalar_product(&
             solve,rhsid,unkno,alpha,&
             MY_TIMING=my_timing,OPENMP=.false.)
        times_max(7) = times_max(7) + my_timing(1)
        times_max(8) = times_max(8) + my_timing(2)
        times_max(9) = times_max(9) + my_timing(1) + my_timing(2)
     end do
     !
     ! Scalar-product with OpenMP
     !
     do ipass = 1,kfl_algebra_operations
        call solver_parallel_scalar_product(&
             solve,rhsid,unkno,alpha,&
             MY_TIMING=my_timing,OPENMP=.true.)
        times_max(10) = times_max(10) + my_timing(1)
        times_max(11) = times_max(11) + my_timing(2)
        times_max(12) = times_max(12) + my_timing(1) + my_timing(2)
     end do

#ifdef NINJA
     !
     ! Matrix-vector and Dot product on GPU
     !
     if( kfl_coo /= 0 ) then
        !call spmv_time(amatr,meshe % coo_rows, meshe % coo_cols , meshe % npoin, rhsid, unkno, times)
     else
        !call spmv_time(amatr,meshe % r_dom, meshe % c_dom, meshe % npoin, rhsid, unkno, times)
     end if

     !call dottime(rhsid,unkno,meshe % npoin,times)
#endif
     !
     ! Normalize operations
     !
     times_max = times_max / real(kfl_algebra_operations,rp)
     !
     ! Output info in appropriate units
     !
     times_max(13) = my_times(14) - my_times(13) ! Allocate and copy memories
     times_max(14) = my_times(15) - my_times(14) ! SMVP computation
     times_max(15) = my_times(15) - my_times(13) ! all

     times_max(16) = my_times(17) - my_times(16) ! Allocate and copy memories
     times_max(17) = my_times(18) - my_times(17) ! Dot product computation
     times_max(18) = my_times(18) - my_times(16) ! all

     times_ave = times_max

     ii = size(times_max,KIND=ip)
     jj = size(times_ave,KIND=ip)
     call PAR_MAX    (ii,times_max,'IN MY CODE')
     call PAR_AVERAGE(jj,times_ave,'IN MY CODE')
     load_balance = times_ave / (times_max+epsilon(1.0_rp))
     coutp(1:4)   = ' '

     call maths_time_unit(maxval(times_max( 1: 6)),coutp(1)(1:3),xfact1) ! SPMV
     call maths_time_unit(maxval(times_max( 7:12)),coutp(2)(1:3),xfact2) ! Dot product
     call maths_time_unit(maxval(times_max(13:15)),coutp(3)(1:3),xfact3) ! SPMV GPU
     call maths_time_unit(maxval(times_max(16:18)),coutp(4)(1:3),xfact4) ! Dot product GPU

     times_ave( 1: 6) = times_ave( 1: 6) * xfact1                        ! SpMV
     times_max( 1: 6) = times_max( 1: 6) * xfact1                        !
     times_max( 7:12) = times_max( 7:12) * xfact2                        ! Dot product
     times_ave( 7:12) = times_ave( 7:12) * xfact2                        !
     times_max(13:15) = times_max(13:15) * xfact3                        ! SPMV with GPU
     times_ave(13:15) = times_ave(13:15) * xfact3                        !
     times_max(16:18) = times_max(16:18) * xfact4                        ! Dot product with GPU
     times_ave(16:18) = times_ave(16:18) * xfact4
     !
     ! SpMV
     !
     routp( 1) = times_ave( 1) ; routp( 2) = times_max( 1) ; routp( 3) = load_balance( 1)
     routp( 4) = times_ave( 2) ; routp( 5) = times_max( 2) ; routp( 6) = load_balance( 2)
     routp( 7) = times_ave( 3) ; routp( 8) = times_max( 3) ; routp( 9) = load_balance( 3)
     routp(10) = times_ave( 4) ; routp(11) = times_max( 4) ; routp(12) = load_balance( 4)
     routp(13) = times_ave( 5) ; routp(14) = times_max( 5) ; routp(15) = load_balance( 5)
     routp(16) = times_ave( 6) ; routp(17) = times_max( 6) ; routp(18) = load_balance( 6)
     !
     ! Dot product
     !
     routp(19) = times_ave( 7) ; routp(20) = times_max( 7) ; routp(21) = load_balance( 7)
     routp(22) = times_ave( 8) ; routp(23) = times_max( 8) ; routp(24) = load_balance( 8)
     routp(25) = times_ave( 9) ; routp(26) = times_max( 9) ; routp(27) = load_balance( 9)
     routp(28) = times_ave(10) ; routp(29) = times_max(10) ; routp(30) = load_balance(10)
     routp(31) = times_ave(11) ; routp(32) = times_max(11) ; routp(33) = load_balance(11)
     routp(34) = times_ave(12) ; routp(35) = times_max(12) ; routp(36) = load_balance(12)
     !
     ! SpMV GPU
     !
     routp(37) = times_ave(13) ; routp(38) = times_max(13) ; routp(39) = load_balance( 13)
     routp(40) = times_ave(14) ; routp(41) = times_max(14) ; routp(42) = load_balance( 14)
     routp(43) = times_ave(15) ; routp(44) = times_max(15) ; routp(45) = load_balance( 15)
     !
     ! Dot product GPU
     !
     routp(46) = times_ave(16) ; routp(47) = times_max(16) ; routp(48) = load_balance( 16)
     routp(49) = times_ave(17) ; routp(50) = times_max(17) ; routp(51) = load_balance( 17)
     routp(52) = times_ave(18) ; routp(53) = times_max(18) ; routp(54) = load_balance( 18)
     !
     ! Output info in appropriate units
     !
     call outfor(76_ip,0_ip,' ')

  end if
 
end subroutine outcpu_operations




