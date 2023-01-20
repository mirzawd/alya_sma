!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!>------------------------------------------------------------------------
!> @defgroup   Direct_Solver
!> Bridge to different direct solvers, like MUMPS, WSMP, PASTIX and Alya sparse direct solver
!> @ingroup    Algebraic_Solver
!> @{
!> @file    mod_direct_solver.f90
!> @author  Guillaume Houzeaux
!> @date    13/06/2014
!> @brief   ToolBox for direct solvers
!> @details ToolBox for direct solvers
!>
!>          Initialization            |
!>          Ordering                  | direct_solver_initialization
!>          Symbolic factorization    |
!>          +--->
!>          |   
!>          | Numerical factorization | direct_solver_factorization
!>          | Solution                | direct_solver_solution
!>          | Partial cleaning        | direct_solver_partialcleaning (if problem is to be solved once more)
!>          |
!>          <---+
!>          Cleaning                  | direct_solver_cleaning
!>
!>          Partial cleaning deallocate temporary matrix, L and U
!>
!>          Input data in DIRECT_SOLVER type:
!>          ---------------------------------
!>          KFL_SOLVER ................ Type of solver
!>          KFL_SYMMETRIC ............. Matrix is symmetric (1) or not (0)
!>          NN ........................ Number of nodes
!>          NDOF ...................... Number of dof per node
!>          NRHS ...................... Number of RHS
!>          NZ,IA(:),JA(:) ............ Matrix CSR format (only for CSR or CSC solvers)
!>          NSKYL,IDIAG(:),ISKYL(:) ... Matrix skyline format (only Cholesky skyline)
!>
!>          Output data in DIRECT_SOLVER type:
!>          ----------------------------------
!>          MEMOR(2) .................. Memory counter (MEMOR(2) is the max memory allocated)
!>          CPUTI(3) .................. CPU time for each phase (init, fact, sol)
!>          NUM_INITIALIZATIONS ....... Number of initializations
!>          NUM_FACTORIZATIONS ........ Number of factorizations
!>          NUM_SOLUTIONS ............. Number of solutions
!>          NAME ...................... Name of the solver
!>
!>          Internal data in DIRECT_SOLVER type:
!>          ------------------------------------
!>          AA(:,:,:) ................. Matrix which can be allocated by this module
!>          IA_NEW(:), JA_NEW(:) ...... Renumbered matrix graph
!>          AA_SKYL(:) ................ ALYA SKYLINE: Factorized matrix in skyline format
!>          IL,JL,IU,JU ............... ALYA CSR: Symbolical factorization of Alya solver
!>          LN,UN ..................... ALYA CSR: Numerical LU for Alya solver
!>          PERMR(:), INVPR(:) ........ Permutation arrays NEW=PERM(OLD), but the opposite for pastix
!>          IPARM(:),DPARM(:) ......... PASTIX: int and real parameters 
!>          PASTIX_DATA ............... PASTIX: identifier 
!>          PASTIX_COMM ............... PASTIX: communicator
!>          NUM_THREADS ............... PASTIX and MUMPS: number of OMP threads to be used
!>
!>                           MUMPS
!>          Matrix assembly  Lesxical ordering: full/partial rows
!>          Matrix format    COO in global lexical numbering
!>          RHS              Centralized
!>
!------------------------------------------------------------------------

module mod_direct_solver

  use def_kintyp_basic,       only : ip,rp,lg
  use def_kintyp_solvers,     only : direct_solver_typ
  use def_master,             only : INOTMASTER,IMASTER
  use def_master,             only : IPARALL
  use def_master,             only : kfl_paral,npart
  use def_solver,             only : SOL_DIRECT_SOLVER_PASTIX
  use def_solver,             only : SOL_DIRECT_SOLVER_WSMP
  use def_solver,             only : SOL_DIRECT_SOLVER_PWSMP
  use def_solver,             only : SOL_DIRECT_SOLVER_ALYA
  use def_solver,             only : SOL_DIRECT_SOLVER_MUMPS
  use def_solver,             only : SOL_DIRECT_SOLVER_SKYLINE_ALYA
  use def_solver,             only : SOL_DIRECT_SOLVER_GPUQR
  use mod_graphs,             only : graphs_permut_metis_postordering
  use mod_graphs,             only : graphs_copyij
  use mod_graphs,             only : graphs_permut_metis_postordering_deallocate
  use mod_graphs,             only : graphs_copyij_deallocate
  use mod_memory_basic,       only : memory_alloca
  use mod_memory_basic,       only : memory_deallo
  use mod_alya_direct_solver, only : alya_Symbolical_CSR_LU
  use mod_alya_direct_solver, only : alya_Numerical_CSR_LU
  use mod_alya_direct_solver, only : alya_CSR_LUSol
  use mod_alya_direct_solver, only : alya_CSR_LUfin
  use mod_alya_direct_solver, only : alya_Numerical_CSR_LU_Deallocate
  use mod_alya_direct_solver, only : alya_cholesky_factorization
  use mod_alya_direct_solver, only : alya_cholesky_solution
  use mod_alya_direct_solver, only : alya_CSR_memor 
  use mod_communications,     only : PAR_MAX,PAR_SUM
  use mod_communications,     only : PAR_AVERAGE
  use mod_parall,             only : PAR_COMM_MY_CODE_WM,PAR_COMM_MY_CODE
  use mod_alya2mumps,         only : direct_solver_initialization_mumps
  use mod_alya2mumps,         only : direct_solver_factorization_mumps
  use mod_alya2mumps,         only : direct_solver_solution_mumps
  use mod_alya2mumps,         only : direct_solver_cleaning_mumps
  
  implicit none 
  private

!------------------------------------------------------------------------  
#ifdef PASTIX
#include "pastix_fortran.h"
  pastix_int_t :: pastix_integer
!#else
!  integer(4)   :: pastix_integer
#endif
!------------------------------------------------------------------------  
  
  real(rp)     :: time1,time2

  interface direct_solver_type_initialization
     module procedure direct_solver_type_initialization_s,&
          &           direct_solver_type_initialization_1
  end interface direct_solver_type_initialization

  public :: direct_solver_factorization                          ! Symbolic factorization
  public :: direct_solver_initialization                         ! Initialization
  public :: direct_solver_solution                               ! Solution of system Ax=b
  public :: direct_solver_cleaning                               ! Complete cleaning of solvers
  public :: direct_solver_factorization_solution_partialcleaning ! One shot solution
  public :: direct_solver_partialcleaning                        ! Partial cleaning (if to be solved once more)
  public :: direct_allocate_temporary_matrix                     ! Allocate tmp maytrix AA
  public :: direct_solver_matrix_size                            ! Get matrix size
  public :: direct_solver_time_initialization                    ! Get initialization CPU time
  public :: direct_solver_time_factorization                     ! Get factorization CPU time
  public :: direct_solver_time_solution                          ! Get solution CPU time
  public :: direct_solver_memory                                 ! Get solver maximum memory allocation
  public :: direct_solver_name                                   ! Get solver name
  public :: direct_solver_statistics                             ! Some CPU time and memory statistics
  public :: direct_solver_type_initialization                    ! Initialize solver type

contains

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    02/03/2016
  !> @brief   Memory
  !> @details Maximum Memory allocated during all the steps
  !>
  !----------------------------------------------------------------------

  function direct_solver_name(direct_solver)
    type(direct_solver_typ), intent(inout) :: direct_solver 
    character(400)                         :: direct_solver_name
    
    if(      direct_solver % kfl_solver == SOL_DIRECT_SOLVER_SKYLINE_ALYA ) then
       direct_solver_name = 'CHOLESKY IN SKYLINE'
    else if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_ALYA ) then
       direct_solver_name = 'DIRECT SPARSE IN CSR '
    else if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_PASTIX ) then
       direct_solver_name = 'PASTIX - CSC'
    else if(direct_solver % kfl_solver == SOL_DIRECT_SOLVER_WSMP ) then
       direct_solver_name = 'WSMP - CSR'
    else if(direct_solver % kfl_solver == SOL_DIRECT_SOLVER_PWSMP ) then
       direct_solver_name = 'PWSMP - CSR'
    else if(direct_solver % kfl_solver == SOL_DIRECT_SOLVER_GPUQR ) then
       direct_solver_name = 'GPU - QR DECOMPOSITION'
    else if(direct_solver % kfl_solver == SOL_DIRECT_SOLVER_MUMPS ) then
       direct_solver_name = 'MUMPS - COO'
    end if
  end function direct_solver_name

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    02/03/2016
  !> @brief   Memory
  !> @details Maximum Memory allocated during all the steps
  !>
  !----------------------------------------------------------------------

  function direct_solver_memory(direct_solver)
    type(direct_solver_typ), intent(in) :: direct_solver 
    integer(8)                          :: direct_solver_memory
    
    direct_solver_memory = direct_solver % memor(2)
    
    if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_PASTIX ) then
#ifdef PASTIX
       direct_solver_memory = direct_solver_memory + int(direct_solver % iparm(IPARM_ALLOCATED_TERMS),8)
#endif
    else if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_WSMP ) then
#ifdef WSMP
       ! direct_solver_memory = direct_solver_memory + int(direct_solver % iparm(IPARM_ALLOCATED_TERMS),8)
#endif
    end if

  end function direct_solver_memory

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    02/03/2016
  !> @brief   Timings
  !> @details Timings
  !>
  !----------------------------------------------------------------------

  function direct_solver_time_initialization(direct_solver)
    type(direct_solver_typ), intent(in) :: direct_solver 
    real(rp)                            :: direct_solver_time_initialization
    direct_solver_time_initialization = direct_solver % cputi(1)
  end function direct_solver_time_initialization
  function direct_solver_time_factorization(direct_solver)
    type(direct_solver_typ), intent(in) :: direct_solver 
    real(rp)                            :: direct_solver_time_factorization
    direct_solver_time_factorization = direct_solver % cputi(2)
  end function direct_solver_time_factorization
  function direct_solver_time_solution(direct_solver)
    type(direct_solver_typ), intent(in) :: direct_solver 
    real(rp)                            :: direct_solver_time_solution
    direct_solver_time_solution = direct_solver % cputi(3)
  end function direct_solver_time_solution

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    02/03/2016
  !> @brief   CPU time and memory statistics
  !> @details Compute some statistics for the direct solver:
  !>          solver_statistics( 1: 3) ... Max. CPU of init/factor/solut over slaves
  !>          solver_statistics( 4: 6) ... Ave. CPU of init/factor/solut over slaves
  !>          solver_statistics( 7: 9) ... Load balance of init/factor/solut
  !>          solver_statistics( 10)   ... Number of initializations
  !>          solver_statistics( 11)   ... Number of factorizations
  !>          solver_statistics( 12)   ... Number of solutions
  !>          solver_statistics( 13)   ... Max. memory over slaves
  !>          solver_statistics( 14)   ... Ave. memory over slaves
  !>          solver_statistics( 15)   ... Memory of current slave
  !>          solver_statistics( 16)   ... Max. fill in
  !>          solver_statistics( 17)   ... Min. fill in
  !>
  !----------------------------------------------------------------------

  subroutine direct_solver_statistics(direct_solver,solver_statistics)
    type(direct_solver_typ), intent(in)  :: direct_solver 
    real(rp),                intent(out) :: solver_statistics(20)
    real(rp)                             :: cputi_max(3),cputi_ave(3),load_balance(3)
    integer(8)                           :: memor_max,memor_ave
    integer(ip)                          :: num_initializations
    integer(ip)                          :: num_factorizations
    integer(ip)                          :: num_solutions,ii
    real(rp)                             :: rmemor_max,rmemor_ave,rmemor_cur
    real(rp)                             :: rfillin_max,rfillin_ave

    solver_statistics     = 0.0_rp
    cputi_max             = 0.0_rp
    cputi_ave             = 0.0_rp
    memor_max             = 0
    memor_ave             = 0
    rfillin_max           = 0.0_rp
    num_initializations   = 0
    num_factorizations    = 0
    num_solutions         = 0
    rfillin_max           = 0.0_rp
    rmemor_ave            = 0.0_rp
    
    if( INOTMASTER ) then
       cputi_max(1)        = direct_solver_time_initialization(direct_solver)
       cputi_max(2)        = direct_solver_time_factorization(direct_solver)
       cputi_max(3)        = direct_solver_time_solution(direct_solver)
       memor_max           = direct_solver_memory(direct_solver)
       rfillin_max         = direct_solver % fillin
       num_initializations = direct_solver % num_initializations
       num_factorizations  = direct_solver % num_factorizations
       num_solutions       = direct_solver % num_solutions
    end if
    cputi_ave              = cputi_max
    memor_ave              = memor_max
    rmemor_max             = real(memor_max,rp)
    rmemor_ave             = real(memor_ave,rp)
    rmemor_cur             = rmemor_ave
    rfillin_ave            = rfillin_max

    call PAR_MAX    (3_ip,cputi_max      , 'IN MY CODE')
    call PAR_AVERAGE(3_ip,cputi_ave      , 'IN MY CODE')
    call PAR_MAX    (rmemor_max          , 'IN MY CODE')
    call PAR_AVERAGE(rmemor_ave          , 'IN MY CODE')
    call PAR_AVERAGE(num_initializations , 'IN MY CODE')
    call PAR_AVERAGE(num_factorizations  , 'IN MY CODE')
    call PAR_AVERAGE(num_solutions       , 'IN MY CODE')
    call PAR_MAX    (rfillin_max         , 'IN MY CODE')
    call PAR_AVERAGE(rfillin_ave         , 'IN MY CODE')

    do ii = 1,3
       if( cputi_max(ii) > epsilon(1.0_rp) ) then 
          load_balance(ii) = cputi_ave(ii) / cputi_max(ii)
       else
          load_balance(ii) = 1.0_rp
       end if
    end do
    
    solver_statistics( 1: 3) = cputi_max(1:3)                  ! Max. CPU of init/factor/solut over slaves
    solver_statistics( 4: 6) = cputi_ave(1:3)                  ! Ave. CPU of init/factor/solut over slaves
    solver_statistics( 7: 9) = load_balance(1:3)               ! Load balance of init/factor/solut
    solver_statistics( 10)   = real(num_initializations,rp)    ! Number of initializations
    solver_statistics( 11)   = real(num_factorizations,rp)     ! Number of factorizations
    solver_statistics( 12)   = real(num_solutions,rp)          ! Number of solutions
    solver_statistics( 13)   = rmemor_max                      ! Max. memory over slaves
    solver_statistics( 14)   = rmemor_ave                      ! Ave. memory over slaves
    solver_statistics( 15)   = rmemor_cur                      ! Memory of current slave
    solver_statistics( 16)   = rfillin_max                     ! Max. fill-in over slaves
    solver_statistics( 17)   = rfillin_ave                     ! Ave. fill-in of current slave

  end subroutine direct_solver_statistics

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    02/03/2016
  !> @brief   Size of the matrix
  !> @details This function gives the size of the matrix 
  !>
  !----------------------------------------------------------------------

  function direct_solver_matrix_size(direct_solver)

    type(direct_solver_typ), intent(inout) :: direct_solver 
    integer(ip)                            :: direct_solver_matrix_size

    if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_SKYLINE_ALYA ) then
       !
       ! Skyline format
       !
       direct_solver_matrix_size = direct_solver % nskyl 
    else 
       !
       ! CSR/CSC format
       !
       direct_solver_matrix_size = direct_solver % nz * direct_solver % ndof * direct_solver % ndof
    end if

  end function direct_solver_matrix_size

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    19/02/2016
  !> @brief   Solve system  
  !> @details Solve system   
  !>
  !----------------------------------------------------------------------

  subroutine direct_solver_factorization_solution_partialcleaning(direct_solver,aa,rhsid,unkno) 

    type(direct_solver_typ), intent(inout)          :: direct_solver 
    real(rp),                intent(in)             :: aa(*)
    real(rp),                intent(in)             :: rhsid(*)
    real(rp),                intent(out)            :: unkno(*)
    !
    ! Factorize matrix
    !
    call direct_solver_factorization(direct_solver,aa,rhsid,unkno) 
    !
    ! Solve system
    !
    call direct_solver_solution(direct_solver,rhsid,unkno) 
    !
    ! Clean solver arrays
    !
    call direct_solver_partialcleaning(direct_solver) 

  end subroutine direct_solver_factorization_solution_partialcleaning

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    19/02/2016
  !> @brief   Clean everything
  !> @details Complete cleaning of the solver
  !>
  !----------------------------------------------------------------------

  subroutine direct_allocate_temporary_matrix(direct_solver) 

    type(direct_solver_typ), intent(inout) :: direct_solver 

    if( .not. associated(direct_solver % aa) ) then
       if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_SKYLINE_ALYA ) then
          call memory_alloca(&
               direct_solver % memor,'DIRECT_SOLVER % AA','direct_solver_partialcleaning',direct_solver % aa,&
               1_ip,1_ip,direct_solver % nskyl)       
       else
          call memory_alloca(&
               direct_solver % memor,'DIRECT_SOLVER % AA','direct_solver_partialcleaning',direct_solver % aa,&
               direct_solver % ndof,direct_solver % ndof,direct_solver % nz)
       end if
    end if

  end subroutine direct_allocate_temporary_matrix

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    19/02/2016
  !> @brief   Clean everything
  !> @details Complete cleaning of the solver
  !>
  !----------------------------------------------------------------------

  subroutine direct_solver_partialcleaning(direct_solver) 

    type(direct_solver_typ), intent(inout) :: direct_solver 
    !
    ! Temporary matrix
    !
    call memory_deallo(direct_solver % memor,'DIRECT_SOLVER % AA','direct_solver_partialcleaning',direct_solver % aa)

    if(      direct_solver % kfl_solver == SOL_DIRECT_SOLVER_SKYLINE_ALYA ) then
       !
       ! Chelsky factorization
       !
       call memory_deallo(direct_solver % memor,'DIRECT_SOLVER % AA_SKYL','direct_solver_partialcleaning',direct_solver % aa_skyl)

    else if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_ALYA ) then
       !
       call alya_Numerical_CSR_LU_Deallocate(&
            direct_solver % ln , direct_solver % un   )
       direct_solver % memor = alya_CSR_memor

    else if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_PASTIX ) then
       ! 
       ! PASTIX: option not avialable
       !
    else if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_WSMP ) then
       !
       ! WSMP:  option not available
       !
    else if(direct_solver % kfl_solver == SOL_DIRECT_SOLVER_PWSMP ) then
       !
       ! PWSMP: Option not available
       !
    else if(direct_solver % kfl_solver == SOL_DIRECT_SOLVER_MUMPS ) then
       !
       !  MUMPS
       !
    end if

  end subroutine direct_solver_partialcleaning

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    19/02/2016
  !> @brief   Clean everything
  !> @details Complete cleaning of the solver
  !>
  !----------------------------------------------------------------------

  subroutine direct_solver_cleaning(direct_solver) 

    type(direct_solver_typ), intent(inout) :: direct_solver 
#ifdef PASTIX
    real(rp)                               :: dummr(2)
#endif

    if(      direct_solver % kfl_solver == SOL_DIRECT_SOLVER_ALYA ) then
       !
       ! Alya
       !
       alya_CSR_memor = direct_solver % memor
       call alya_CSR_LUfin(&
            direct_solver % IL , direct_solver % JL , &
            direct_solver % LN , direct_solver % IU , &
            direct_solver % JU , direct_solver % UN   )
       direct_solver % memor = alya_CSR_memor
      call graphs_copyij_deallocate(                          &
           direct_solver % ia_new, direct_solver % ja_new,    &
           direct_solver % memor                              )       
      call graphs_permut_metis_postordering_deallocate(       &
           direct_solver % permr  , direct_solver % invpr,    &
           direct_solver % memor  ,                           &
           PERMR_NAME = 'PERMR_'//trim(direct_solver % name), &
           INVPR_NAME = 'INVPR_'//trim(direct_solver % name)  )

#ifdef NINJA
    else if (direct_solver % kfl_solver == SOL_DIRECT_SOLVER_GPUQR) then
       !!call gpuqrdestroy()
#endif

    else if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_PASTIX ) then
       ! 
       ! PASTIX
       !
#ifdef PASTIX
       direct_solver % iparm(IPARM_START_TASK) = API_TASK_CLEAN   
       direct_solver % iparm(IPARM_END_TASK)   = API_TASK_CLEAN  
       call d_pastix_fortran(&
            direct_solver % pastix_data,direct_solver % pastix_comm,&
            direct_solver % nn,direct_solver % ia,direct_solver % ja,dummr,&
            direct_solver % invpr,direct_solver % permr,&
            dummr,direct_solver % nrhs,direct_solver % iparm,direct_solver % dparm)
#endif
    else if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_WSMP ) then
       !
       ! WSMP
       !

#ifdef WSMP
       !       
       !   WSMP: Option not available   
       !     
#endif 
    else if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_PWSMP ) then
       !
       ! PWSMP
       !

#ifdef PWSMP
       !       
       !   PWSMP: Option not available   
       !     
#endif 
    else if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_MUMPS) then
       !
       ! MUMPS
       !
#ifdef MUMPS
       !
       ! MUMPS
       !
       call direct_solver_cleaning_mumps(direct_solver)
 
#endif
    end if

  end subroutine direct_solver_cleaning

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    19/02/2016
  !> @brief   Solve system  
  !> @details Solve system   
  !>
  !----------------------------------------------------------------------

  subroutine direct_solver_solution(direct_solver,rhsid,unkno,aa) 

    type(direct_solver_typ), intent(inout)          :: direct_solver 
    real(rp),                intent(in), optional   :: aa(*)
    real(rp),                intent(in)             :: rhsid(*)
    real(rp),                intent(out)            :: unkno(*)
    integer(ip)                                     :: ntotn,ierro
!    real(rp),                pointer                :: resid(:)
!    real(rp)                                        :: numer,denom


    call cputim(time1) 
    direct_solver % num_solutions = direct_solver % num_solutions + 1

    if(      direct_solver % kfl_solver == SOL_DIRECT_SOLVER_SKYLINE_ALYA ) then
       !
       ! Alya Cholesky with skyline format: nothing to do
       ! 
       ntotn = direct_solver % nn * direct_solver % ndof
       unkno(1:ntotn) = rhsid(1:ntotn)
       call alya_cholesky_solution(&
            ntotn,direct_solver % nskyl,direct_solver % iskyl,&
            direct_solver % nrhs,direct_solver % aa_skyl,unkno,&
            ntotn,ierro)    
       if( ierro /= 0 ) call runend('DIRECT_SOLVER_FACTORIZATION: ERROR IN CHOLESKY SOLUTION')

    else if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_ALYA ) then
       !
       ! Alya
       !  
       alya_CSR_memor = direct_solver % memor
       call alya_CSR_LUsol(&
            direct_solver % nn    , direct_solver % ndof  , &
            direct_solver % invpr , direct_solver % invpr , &
            direct_solver % il    , direct_solver % jl    , &
            direct_solver % ln    , direct_solver % iu    , &
            direct_solver % ju    , direct_solver % un    , &
            rhsid                 , unkno                 )
       direct_solver % memor = alya_CSR_memor 

#ifdef NINJA
    else if (direct_solver % kfl_solver == SOL_DIRECT_SOLVER_GPUQR) then
       !!call GPUQRSOLVE(rhsid,unkno)
#endif

    else if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_PASTIX ) then
       !
       ! Pastix
       !
       call direct_solver_solution_pastix(direct_solver,rhsid,unkno)
    else if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_WSMP ) then
       !
       ! WSMP
       !
       !          if (present(aa)) then

       call direct_solver_solution_wsmp(direct_solver,rhsid,unkno)  
       !          else

       !             call direct_solver_solution_wsmp(direct_solver,rhsid,unkno,direct_solver % aa) 
       !          end if
    else if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_PWSMP ) then
       !
       ! PWSMP
       !
       call direct_solver_solution_pwsmp(direct_solver,rhsid,unkno)  

    else if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_MUMPS ) then
       ! 
       !  MUMPS
       !
       call direct_solver_solution_mumps(direct_solver,rhsid,unkno)
       
    else

       call runend('DIRECT_SOLVER_SOLUTION: not coded')

    end if
    !
    ! Residual of the equation: ||b-Ax||/||b||
    !
    !allocate(resid(direct_solver % nn*direct_solver % ndof))
    !call matrix_CSR_SMVP(direct_solver % nn,direct_solver % ndof,direct_solver%ia,direct_solver%ja,direct_solver % aa,unkno,resid)
    !ntotn = direct_solver % nn * direct_solver % ndof
    !resid(1:ntotn) = resid(1:ntotn) - rhsid(1:ntotn)
    !numer = dot_product(resid(1:ntotn),resid(1:ntotn))
    !denom = dot_product(rhsid(1:ntotn),rhsid(1:ntotn))
    !write(kfl_paral+900,*) sqrt(numer/denom)

    call cputim(time2) 
    direct_solver % cputi(3) = direct_solver % cputi(3) + (time2-time1)

  end subroutine direct_solver_solution

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    19/02/2016
  !> @brief   Factorize matrix 
  !> @details LU Factorization  
  !>
  !----------------------------------------------------------------------

  subroutine direct_solver_factorization(direct_solver,aa,rhsid,unkno) 

    type(direct_solver_typ), intent(inout)          :: direct_solver       
    real(rp),               intent(in),   optional  :: aa(*)
    integer(ip)                                     :: ierro,ntotn
!    integer(ip)                                     :: ii,nnz,ndofn
    real(rp),               intent(in),optional     :: rhsid(*)
    real(rp),               intent(out),optional    :: unkno(*)

    integer(ip), pointer                           :: lninv_lex(:)
    nullify(lninv_lex)

    call cputim(time1)
    direct_solver % num_factorizations = direct_solver % num_factorizations + 1

    if(      direct_solver % kfl_solver == SOL_DIRECT_SOLVER_SKYLINE_ALYA ) then
       !
       ! Alya Cholesky with skyline format: nothing to do
       ! 
       call memory_alloca(direct_solver % memor,'DIRECT_SOLVER % AA_SKYL',&
            'solver_direct_solver_factorization',direct_solver % aa_skyl,direct_solver % nskyl)
       if( present(aa) ) then
          direct_solver % aa_skyl(1:direct_solver % nskyl) = aa(1:direct_solver % nskyl)
       else  
          direct_solver % aa_skyl(1:direct_solver % nskyl) = direct_solver % aa(1,1,1:direct_solver % nskyl)
       end if

       ntotn = direct_solver % nn * direct_solver % ndof
       call alya_cholesky_factorization(&
            ntotn,direct_solver % nskyl,&
            direct_solver % iskyl,direct_solver % aa_skyl,ierro)
       if( ierro /= 0 ) call runend('DIRECT_SOLVER_FACTORIZATION: ERROR IN CHOLESKY FACTORIZATION')

    else if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_ALYA ) then
       !
       ! Alya sparse in CSR format
       !            
       alya_CSR_memor = direct_solver % memor
       if( present(aa) ) then

          call alya_Numerical_CSR_LU(&                      
               direct_solver % nn     , direct_solver % ndof   , &
               direct_solver % ia_new , direct_solver % ja_new , &
               aa                     , direct_solver % IL     , &
               direct_solver % JL     , direct_solver % LN     , &
               direct_solver % IU     , direct_solver % JU     , &
               direct_solver % UN     , ierro                  , & 
               direct_solver % invpr  , direct_solver % permr  , &
               direct_solver % ia     , direct_solver % ja     )

       else

          call alya_Numerical_CSR_LU(&                      
               direct_solver % nn     , direct_solver % ndof   , &
               direct_solver % ia_new , direct_solver % ja_new , &
               direct_solver % aa     , direct_solver % IL     , &
               direct_solver % JL     , direct_solver % LN     , &
               direct_solver % IU     , direct_solver % JU     , &
               direct_solver % UN     , ierro                  , & 
               direct_solver % invpr  , direct_solver % permr  , &
               direct_solver % ia     , direct_solver % ja     )

       end if
       if( ierro /= 0 ) call runend('DIRECT_SOLVER_FACTORIZATION: ERROR')
       direct_solver % memor = alya_CSR_memor 

    else if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_PASTIX ) then
       !
       ! Pastix
       !
       if( present(aa) ) then
          call direct_solver_factorization_pastix(direct_solver,aa) 
       else
          call direct_solver_factorization_pastix(direct_solver,direct_solver % aa) 
       end if

    else if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_WSMP ) then
       !
       ! WSMP
       !
       if( present(aa) ) then
          call direct_solver_factorization_wsmp(direct_solver,aa)
       else
          call direct_solver_factorization_wsmp(direct_solver,direct_solver % aa)
       end if

    else if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_PWSMP ) then
       !
       ! PWSMP
       !

       !
       ! Renumber in lexical order
       !
       !       call renumbering_lexical_order(meshe(ndivi),lninv_lex)      ! lo elimino porque no es compatible con mumps - guillaume ya volvera
       !       a meter wsmp desde 0
       !do ii = 1,direct_solver % nz
       !   direct_solver % ja(ii) = lninv_lex(direct_solver % ja(ii))
       !end do

       !write(kfl_paral+90,'(50(1x,i3))') direct_solver % ia(1:direct_solver % nn+1)
       !write(kfl_paral+90,'(50(1x,i3))') direct_solver % ja(1:direct_solver % nz)
       !call matrix_output_dense_format(direct_solver % nn,1_ip,direct_solver % ia,direct_solver % ja,an)

!!$       if( present(aa) ) then
!!$          nnz   = solve_sol(1) % nnz
!!$          ndofn = solve_sol(1) % ndofn
!!$          if( INOTMASTER ) then
!!$             call memory_alloca(direct_solver % memor,'DIRECT_SOLVER % AA','solver_direct_solver_initialize',direct_solver % aa,ndofn,ndofn,nnz)
!!$             call matrix_copy_matrix(solve_sol(1) % nequa,solve_sol(1) % ndofn,solve_sol(1) % ia,solve_sol(1) % ja,aa,direct_solver % aa)
!!$             !call matrix_remove_null_coefficients(solve_sol(1) % nequa,solve_sol(1) % ndofn,solve_sol(1) % nnz,&
!!$             !     solve_sol(1) % ia,solve_sol(1) % ja,direct_solver % aa)
!!$          else
!!$             call memory_alloca(direct_solver % memor,'DIRECT_SOLVER % AA','solver_direct_solver_initialize',direct_solver % aa,1_ip,1_ip,1_ip)
!!$          end if
!!$          call direct_solver_factorization_pwsmp(direct_solver,direct_solver % aa)
!!$          call memory_deallo(direct_solver % memor,'DIRECT_SOLVER % AA','solver_direct_solver_initialize',direct_solver % aa)
!!$       else
!!$          call direct_solver_factorization_pwsmp(direct_solver,direct_solver % aa)
!!$       end if

    else if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_MUMPS ) then
       !
       ! MUMPS
       !
       if( present(aa) ) then
          call direct_solver_factorization_mumps(direct_solver,aa)
       else
          call direct_solver_factorization_mumps(direct_solver,direct_solver % aa)
       end if

#ifdef NINJA
    else if(direct_solver % kfl_solver == SOL_DIRECT_SOLVER_GPUQR) then
       !call GPUQRFACTORIZATION(&
       !     direct_solver % nn    , direct_solver % ndof  , &
       !     direct_solver % aa    , direct_solver % ia    , &
       !     direct_solver % ja    , direct_solver % permr )
#endif
    else

       call runend('DIRECT_SOLVER_FACTORIZATION: not coded')
    end if

    call cputim(time2) 
    direct_solver % cputi(2) = direct_solver % cputi(2) + (time2-time1)

  end subroutine direct_solver_factorization

  !-----------------------------------------------------------------------
  !
  ! This routine solves linear systems using PASTIX.
  !
  !-----------------------------------------------------------------------

  subroutine direct_solver_solution_pastix(direct_solver,rhsid,unkno)
    type(direct_solver_typ), intent(inout) :: direct_solver 
    real(rp),                intent(in)    :: rhsid(*)
    real(rp),                intent(out)   :: unkno(*)
#ifdef PASTIX
    integer(ip)                            :: ii
    real(rp)                               :: dummr(2)

    direct_solver % iparm(IPARM_START_TASK) = API_TASK_SOLVE 
    direct_solver % iparm(IPARM_END_TASK)   = API_TASK_SOLVE

    do ii = 1,direct_solver % nn * direct_solver % ndof
       unkno(ii) = rhsid(ii)
    end do
    call d_pastix_fortran(&
         direct_solver % pastix_data,direct_solver % pastix_comm,&
         direct_solver % nn,direct_solver % ia,direct_solver % ja,dummr,&
         direct_solver % invpr,direct_solver % permr,&
         unkno,direct_solver % nrhs,direct_solver % iparm,direct_solver % dparm)
#endif

  end subroutine direct_solver_solution_pastix

  !----------------------------------------------------------------------
  ! 
  ! This routine solves linear systems using WSMP.
  !
  !----------------------------------------------------------------------

  subroutine direct_solver_solution_wsmp(direct_solver,rhsid,unkno)
    type(direct_solver_typ), intent(inout) :: direct_solver 
    real(rp),                intent(in)    :: rhsid(*)
    !      real(rp),                intent(in)    :: aa(*)
    real(rp),                intent(out)   :: unkno(*)
#ifdef WSMP
    real(rp)                               :: dummr(2)
    integer(ip)                            :: nn_wsmp


    nn_wsmp = direct_solver % nn * direct_solver % ndof    


    direct_solver % iparm(2) = 3 
    direct_solver % iparm(3) = 3


    call wgsmp(&
         nn_wsmp, direct_solver % ia, direct_solver % ja, &
         dummr,unkno,direct_solver % ldb,direct_solver % nrhs,dummr,direct_solver % iparm, direct_solver % dparm)

    !       call wgsmp(&
    !            direct_solver % nn, direct_solver % ia, direct_solver % ja, &
    !            aa,unkno,direct_solver % ldb,direct_solver % nrhs,dummr,direct_solver % iparm, direct_solver % dparm)

#endif

  end subroutine direct_solver_solution_wsmp

  !----------------------------------------------------------------------
  !
  ! This routine solves linear systems using PWSMP.
  !
  !----------------------------------------------------------------------

  subroutine direct_solver_solution_pwsmp(direct_solver,rhsid,unkno)
    type(direct_solver_typ), intent(inout) :: direct_solver
    real(rp),                intent(in)    :: rhsid(*)
    !      real(rp),                intent(in)    :: aa(*)
    real(rp),                intent(out)   :: unkno(*)
#ifdef PWSMP
    real(rp)                               :: dummr(2)
    integer(ip)                            :: nn_wsmp

    if (INOTMASTER) then
       nn_wsmp = direct_solver % nn * direct_solver % ndof


       direct_solver % iparm(2) = 3
       direct_solver % iparm(3) = 3

       call pwgsmp(&
            nn_wsmp,direct_solver % ia, direct_solver % ja, dummr, &
            unkno,direct_solver % nn, direct_solver %nrhs, dummr, direct_solver% iparm, direct_solver % dparm)
    end if
#endif

  end subroutine direct_solver_solution_pwsmp

  !-----------------------------------------------------------------------
  !
  ! This routine factorize a matrix using PASTIX.
  !
  !-----------------------------------------------------------------------

  subroutine direct_solver_factorization_pastix(direct_solver,aa)
    type(direct_solver_typ), intent(inout) :: direct_solver 
    real(rp),                intent(in)    :: aa(*)
#ifdef PASTIX
    real(rp)                               :: dummr(2)

    direct_solver % iparm(IPARM_START_TASK) = API_TASK_NUMFACT 
    direct_solver % iparm(IPARM_END_TASK)   = API_TASK_NUMFACT 
    call d_pastix_fortran(&
         direct_solver % pastix_data,direct_solver % pastix_comm,&
         direct_solver % nn,direct_solver % ia,direct_solver % ja,aa,&
         direct_solver % invpr,direct_solver % permr,&
         dummr,direct_solver % nrhs,direct_solver % iparm,direct_solver % dparm)
#endif

  end subroutine direct_solver_factorization_pastix

  !-----------------------------------------------------------------------
  !
  ! This routine factorize a matrix using WSMP.
  !
  !-----------------------------------------------------------------------

  subroutine direct_solver_factorization_wsmp(direct_solver,aa)
    type(direct_solver_typ), intent(inout) :: direct_solver 
    real(rp),                intent(in)    :: aa(*)
#ifdef WSMP 
    real(rp)                               :: dummr(2)
    integer(ip)                            :: nn_wsmp

    nn_wsmp = direct_solver % nn * direct_solver % ndof
    !
    ! Analysis and reordering
    !
    direct_solver % iparm(2) = 1 
    direct_solver % iparm(3) = 1 
    call wgsmp(&
         nn_wsmp, direct_solver % ia, direct_solver % ja, aa, &
         dummr,direct_solver % ldb,direct_solver % nrhs,dummr,direct_solver % iparm, direct_solver % dparm) 
    !
    ! Factorization
    !
    direct_solver % iparm(2) = 2 
    direct_solver % iparm(3) = 2 
    call wgsmp(&
         nn_wsmp, direct_solver % ia, direct_solver % ja, aa, &
         dummr,direct_solver % ldb,direct_solver % nrhs,dummr,direct_solver % iparm, direct_solver % dparm)
#endif

  end subroutine direct_solver_factorization_wsmp

  !-----------------------------------------------------------------------
  !
  ! This routine factorize a matrix using WSMP.
  !
  !-----------------------------------------------------------------------

  subroutine direct_solver_factorization_pwsmp(direct_solver,aa)
    type(direct_solver_typ), intent(inout) :: direct_solver 
    real(rp),                intent(inout) :: aa(*)
#ifdef PWSMP
    real(rp)                               :: dummr(2)
#endif
    integer(ip)                            :: nn_wsmp!,ii

    !if (IMASTER) then
    !nn_wsmp = 0
    !       call pwgsmp(&
    !            nn_wsmp,dummr, dummr, dummr, &
    !            dummr,direct_solver % ldb, direct_solver %nrhs, dummr, direct_solver % iparm, direct_solver % dparm)
    !else

    if (INOTMASTER) then

       nn_wsmp = direct_solver % nn * direct_solver % ndof

       !do ii = 1, direct_solver % nz 
       !   print*,'aa',aa(ii)
       !end do
       !
       ! Analysis and reordering
       !
       direct_solver % iparm(2) = 1 
       direct_solver % iparm(3) = 1
       !direct_solver % iparm(8) = 1
       !print*,size(aa)
       !
       ! Renumber matrix in increasing order
       !
!!$       call matrix_heap_sort(2_ip,solve_sol(1) % ndofn,solve_sol(1) % nunkn,solve_sol(1) % ia,solve_sol(1) % ja,aa)

#ifdef PWSMP

       print*,'direct_solver % nn',direct_solver % nz,'ja_full', direct_solver % ja
       call pwgsmp(&
            direct_solver % nn, direct_solver % ia, direct_solver % ja, aa, &
            dummr,direct_solver % nn,direct_solver % nrhs,dummr,direct_solver % iparm, direct_solver % dparm)
       !print *,'Analysis complete in time - ',wsmprtc()-waltime
       if (direct_solver % iparm(64) /= 0) then
          print *,'The following ERROR was detected: ',direct_solver % iparm(64)
          !return
       end if
       print *,'Number of expected nonzeros in factors = 1000 X ',direct_solver % iparm(24)
       print *,'Number of expected FLOPS in factorization = ',direct_solver % dparm(24)
#else
       print*,'PWSMP symbolic factorization'
#endif
       !print*,'aa_be4_factor', aa 
       !print*,'ja_factor',direct_solver % ja
       !print*,'ia_factor',direct_solver % ia
       !
       ! Factorization
       !       
       direct_solver % iparm(2) = 2
       direct_solver % iparm(3) = 2

#ifdef PWSMP
       call pwgsmp(&
            direct_solver % nn,direct_solver % ia, direct_solver % ja, aa, &
            dummr, direct_solver % nn, direct_solver %nrhs, dummr, direct_solver % iparm, direct_solver % dparm)
       !print *,'LU complete in time - ',waltime
       print *,'Actual number of nonzeros in factors = 1000 X ',direct_solver % iparm(23)
       print *,'Actual number of FLOPS in factorization = ',direct_solver % dparm(23)
       !print *,'Factorization MegaFlops = ',(direct_solver % dparm(23) * 1.d-6) / waltime
       if (direct_solver % iparm(64) .ne. 0) then
          print *,'The following ERROR was detected: ',direct_solver % iparm(64)
          !return
       end if
#else
       print*,'PWSMP numerical factorization'
#endif
    end if

  end subroutine direct_solver_factorization_pwsmp
  

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    07/10/2014
  !> @brief   Direct solvers
  !> @details Direct solvers:
  !>          1. Ordering/Renumbering:
  !>             Input:  NN, IA, JA
  !>             Output: PERMR, INVPR, IA_NEW, JA_NEW
  !>          2. Symbolical factorization
  !>          3. Pastix: analyze
  !>
  !----------------------------------------------------------------------

  subroutine direct_solver_initialization(direct_solver,message)

    type(direct_solver_typ), intent(inout)         :: direct_solver 
    character(*),            intent(in),  optional :: message
#if defined (PASTIX) || defined(WSMP) || defined(PWSMP)
    real(rp)                                       :: dummr(2)  
#endif
    integer(ip)                                    :: nn_wsmp
#ifdef PASTIX
    integer(ip)                                    :: num_threads 
    integer(ip)                                    :: INCOMM   
    integer(ip)                                    :: istat4,my_rank4,jcolor4,my_rank_wm4,ikey4,ii 
#endif

    integer(ip), pointer                           :: lninv_lex(:)

    nullify(lninv_lex)

    call cputim(time1)
    direct_solver % num_initializations = direct_solver % num_initializations + 1
    direct_solver % name = trim(direct_solver_name(direct_solver))

    if( present(message) ) then
       if( trim(message) == 'NOT MASTER' .and. IMASTER ) return
    end if

    if(       direct_solver % kfl_solver == SOL_DIRECT_SOLVER_SKYLINE_ALYA ) then

       !
       ! Alya Cholesky with skyline format: nothing to do
       ! 
       continue 

    else if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_ALYA ) then      
       !
       ! Alya sparse BCSR format
       !
       call graphs_copyij(                                   &
            direct_solver % nn     , direct_solver % ia,     &
            direct_solver % ja     , direct_solver % ia_new, &
            direct_solver % ja_new , direct_solver % memor   )
       
       if( associated(direct_solver % invpr) ) then
          call runend('DIRECT SOLVER: IMPOSSIBLE OPTION')
       else
          call graphs_permut_metis_postordering(                  &
               direct_solver % nn     , direct_solver % nz,       &
               direct_solver % ia_new , direct_solver % ja_new,   &
               direct_solver % permr  , direct_solver % invpr,    &
               PERMR_NAME = 'PERMR_'//trim(direct_solver % name), &
               INVPR_NAME = 'INVPR_'//trim(direct_solver % name), &
               memor=direct_solver % memor)
       end if
       !
       ! Symbolical factorization
       !
       alya_CSR_memor = direct_solver % memor 
       call alya_Symbolical_CSR_LU(&
            direct_solver % nn     , direct_solver % ia_new,  &
            direct_solver % jA_new , direct_solver % iL,      &
            direct_solver % jL     , direct_solver % iU,      &
            direct_solver % jU     , direct_solver % fillin   )
       direct_solver % memor = alya_CSR_memor

#ifdef NINJA
    else if (direct_solver % kfl_solver == SOL_DIRECT_SOLVER_GPUQR) then

#endif
    else if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_PASTIX ) then       
       !
       ! Pastix
       !
#ifdef PASTIX
       if( rp /= 8 )                    call runend('PASTIX CAN ONLY WORK WITH 8-BYTES REAL')
       if( kind(pastix_integer) /= ip ) call runend('PASTIX CAN ONLY WORK WITH 4-BYTES INTEGER')
       call memory_alloca(direct_solver % memor,'DIRECT_SOLVER % PERMR','solver_direct_solver_initialize',direct_solver % permr,direct_solver % nn)
       call memory_alloca(direct_solver % memor,'DIRECT_SOLVER % INVPR','solver_direct_solver_initialize',direct_solver % invpr,direct_solver % nn)
       call memory_alloca(direct_solver % memor,'DIRECT_SOLVER % IPARM','solver_direct_solver_initialize',direct_solver % iparm,IPARM_SIZE)
       call memory_alloca(direct_solver % memor,'DIRECT_SOLVER % DPARM','solver_direct_solver_initialize',direct_solver % dparm,DPARM_SIZE)
       !
       ! Initialization
       ! pastix_data is modified by pastix
       ! pastix_comm is given by Alya
       !
#ifndef MPI_OFF
       direct_solver % pastix_comm                      = MPI_COMM_SELF
       direct_solver % pastix_comm                      = 0
#else
       direct_solver % pastix_comm                      = 0
#endif
       direct_solver % pastix_data                      = 0
       direct_solver % iparm(IPARM_MODIFY_PARAMETER)    = API_NO        ! Use default parameters
       direct_solver % iparm(IPARM_START_TASK)          = API_TASK_INIT
       direct_solver % iparm(IPARM_END_TASK)            = API_TASK_INIT


       call d_pastix_fortran(&
            direct_solver % pastix_data,direct_solver % pastix_comm,&
            direct_solver % nn,direct_solver % ia,direct_solver % ja,&
            dummr,direct_solver % invpr,direct_solver % permr,&
            dummr,direct_solver % nrhs,direct_solver % iparm,direct_solver % dparm)  
       !
       ! 1. Ordering: uses IA and JA to compute PERMR and INVPR
       ! 2. Symbolic factorization LU
       ! 3. Task mapping and scheduling
       ! ATTENTION: For PASTIX, permutation and inverse permutation have
       ! opposite meaning as in Alya. In Alya, NEW = PERMR(OLD)
       ! This is why we have changed the order of parameters in PASTIX calls.
       ! That is, now PERMR has the same meaning for ALL solvers
       !
       direct_solver % iparm(IPARM_THREAD_NBR)          = direct_solver % num_threads
       direct_solver % iparm(IPARM_VERBOSE)             = API_VERBOSE_NOT


       if( direct_solver % kfl_symmetric == 1 ) then
          direct_solver % iparm(IPARM_SYM)              = API_SYM_YES
          direct_solver % iparm(IPARM_FACTORIZATION)    = API_FACT_LLT
       else
          direct_solver % iparm(IPARM_SYM)              = API_SYM_NO
          direct_solver % iparm(IPARM_FACTORIZATION)    = API_FACT_LU
       end if
       !direct_solver % iparm(IPARM_SYM)                 = API_SYM_NO

       direct_solver % iparm(IPARM_TRANSPOSE_SOLVE)      = API_YES
       direct_solver % iparm(IPARM_MATRIX_VERIFICATION) = API_NO
       direct_solver % iparm(IPARM_DOF_NBR)             = direct_solver % ndof
       direct_solver % iparm(IPARM_RHS_MAKING)          = API_RHS_B
       direct_solver % iparm(IPARM_START_TASK)          = API_TASK_ORDERING 
       direct_solver % iparm(IPARM_END_TASK)            = API_TASK_ANALYSE

       call d_pastix_fortran(&
            direct_solver % pastix_data,direct_solver % pastix_comm,&
            direct_solver % nn,direct_solver % ia,direct_solver % ja,&
            dummr,direct_solver % invpr,direct_solver % permr,&
            dummr,direct_solver % nrhs,direct_solver % iparm,direct_solver % dparm)    
       !
       ! Output info from MAPHYS
       !
       direct_solver % fillin = direct_solver % dparm(DPARM_FILL_IN)
#else
       call runend('MACRO PASTIX IS REQUIRED: RECOMPILE ALYA WITH -DPASTIX AND CORRESPONDING INCLUDE')
#endif   

    else if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_WSMP ) then
       !
       ! WSMP
       !
       call memory_alloca(direct_solver % memor,'DIRECT_SOLVER % IPARM','solver_direct_solver_initialization',direct_solver % iparm,64_ip)
       call memory_alloca(direct_solver % memor,'DIRECT_SOLVER % DPARM','solver_direct_solver_initialization',direct_solver % dparm,64_ip)

       direct_solver % ldb  = direct_solver % ndof * direct_solver % nn
       direct_solver % nrhs = 1
       nn_wsmp = direct_solver % ndof * direct_solver % nn

       direct_solver % iparm(1) = 0         
       direct_solver % iparm(2) = 0
       direct_solver % iparm(3) = 0
       !
       ! Initialization
       !
#ifdef WSMP
       call wgsmp(&
            nn_wsmp, direct_solver % ia, direct_solver % ja, dummr, &
            dummr,direct_solver % ldb,direct_solver % nrhs, dummr, direct_solver % iparm, direct_solver % dparm)    
#else
       print*,'WSMP initialization'
#endif

    else if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_PWSMP ) then

       !call wsetmaxthrds(1)
       !print*,'AA=',kfl_paral
       !call wsetmpicomm(PAR_COMM_MY_CODE)
       !print*,'BB=',kfl_paral
       !call runend('O.K.!')
       !
       ! PWSMP
       !
       !istat4 = 0
       !my_rank_wm4 = 0
       !call MPI_Init(istat4)

       !if (istat4 .ne. 0) then
       ! print*,'MPI_initialization_error'
       ! stop
       !end if

       !PAR_COMM_MY_CODE = MPI_COMM_WORLD
       !!call MPI_COMM_RANK(PAR_COMM_MY_CODE,my_rank4,istat4)

       !!if( my_rank4 == 0 ) then
       !!jcolor4 = 0
       !!ikey4 = 0
       !!else
       !!jcolor4 = 1
       !!ikey4 = my_rank4
       !!end if

       !!call MPI_COMM_SPLIT(PAR_COMM_MY_CODE,jcolor4,ikey4,PAR_COMM_MY_CODE_WM,istat4)
       !!call MPI_COMM_RANK(PAR_COMM_MY_CODE_WM,my_rank_wm4,istat4)

       !print*,'old, new=',my_rank4,my_rank_wm4

       !!if( my_rank4 /= 0 ) then
!       call renumbering_lexical_order(meshe(ndivi),lninv_lex)  !  lo elimino guillauem comenzara wsmp de 0

       if( INOTMASTER ) then
          !
          ! Renumber in lexical order
          !
          !do ii = 1,direct_solver % nz
          !   direct_solver % ja(ii) = lninv_lex(direct_solver % ja(ii))
          !end do
#ifdef PWSMP
          call wsetmaxthrds(1)
#endif
          !print*,'pasa1=',direct_solver % nz,size(direct_solver % ja)
#ifdef PWSMP
          call wsetmpicomm(PAR_COMM_MY_CODE_WM)
#endif
          print*,'pasa2'
          !!end if
          !if (INOTMASTER) then
          call memory_alloca(direct_solver % memor,'DIRECT_SOLVER % IPARM','solver_direct_solver_initialization',direct_solver % iparm,64_ip)
          call memory_alloca(direct_solver % memor,'DIRECT_SOLVER % DPARM','solver_direct_solver_initialization',direct_solver % dparm,64_ip)
          !if (IMASTER) then
          !nn_wsmp = 0
          direct_solver % iparm(1) = 0
          direct_solver % iparm(2) = 0
          direct_solver % iparm(3) = 0
          !direct_solver % iparm(4) = 0
#ifdef PWSMP
          call pwgsmp(&
               direct_solver % nn, direct_solver % ia, direct_solver % ja, dummr, &
               dummr,direct_solver % ldb,direct_solver % nrhs, dummr, direct_solver % iparm, direct_solver % dparm)
#else
          print*,'PWSMP initialization 2'
#endif
          !else
          !call memory_alloca(direct_solver % memor,'DIRECT_SOLVER % IPARM','solver_direct_solver_initialization',direct_solver % iparm,64_ip)
          !call memory_alloca(direct_solver % memor,'DIRECT_SOLVER % DPARM','solver_direct_solver_initialization',direct_solver % dparm,64_ip)
          !print*,'A=',kfl_paral
          !call wsetmaxthrds(1)
          !print*,'B=',kfl_paral,PAR_COMM_MY_CODE_WM

          !call wsetmpicomm(PAR_COMM_MY_CODE_WM)
          !print*,'C=',kfl_paral
          !!          call wsetmaxthrds(1)
          !!print*,'D=',kfl_paral
          !!          direct_solver % ldb  = direct_solver % ndof * direct_solver % nn
          !!          direct_solver % nrhs = 1
          !!          nn_wsmp = direct_solver % ndof * direct_solver % nn

          !!print*,'ini0=',kfl_paral
          !!      direct_solver % iparm(1) = 0
          !!      direct_solver % iparm(2) = 0
          !!      direct_solver % iparm(3) = 0
          !
          ! Initialization
          !
          !print*,'ini1=',kfl_paral,direct_solver % nn,associated(direct_solver % ia),associated(direct_solver % ja),size(direct_solver % ia),size(direct_solver % ja)

          !!       call pwgsmp(&
          !!            nn_wsmp, direct_solver % ia, direct_solver % ja, dummr, &
          !!            dummr,direct_solver % ldb,direct_solver % nrhs, dummr, direct_solver % iparm, direct_solver % dparm)
          !end if
          print*,'ini2'
       end if

    else if( direct_solver % kfl_solver == SOL_DIRECT_SOLVER_MUMPS ) then
       !
       ! MUMPS initialization
       !
       call direct_solver_initialization_mumps(direct_solver)

    end if

    call cputim(time2)
    direct_solver % cputi(1) = direct_solver % cputi(1) + (time2-time1)  


  end subroutine direct_solver_initialization

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    07/10/2014
  !> @brief   Direct solvers
  !> @details Direct solvers:
  !> @param thread_icntl     parameter controls (5) of multi-threadings
  !>       - 1 TOPO_BIND: binding method of threads to cores
  !>                      - 0 no binding
  !>                      - 1 thread to core binding
  !>                      - 2 grouped binding
  !>       - 2 TOPO_N:    number of nodes
  !>       - 3 TOPO_C:    number of cores per node
  !>       - 4 TOPO_T:    number of threads per process
  !>       - 5 TOPO_P:    number of MPI processes (in total)
  !>
  !>       Example: On MN, each node has TOPO_C=16 cores 
  !>       Using TOPO_N=3 nodes of Marenostrum, for a total of TOPO_P=8*3=24 processes 
  !>       and TOPO_T=2 threads each, we have:
  !>       MPI_rank   Core number
  !>         0           0
  !>         1           2
  !>         2           4
  !>         3           6
  !>         4           8
  !>         5          10
  !>         6          12
  !>         7          14
  !>         8           0
  !>         9           2
  !>        10           4
  !>        11           6
  !>        12           8
  !>        13          10
  !>        14          12
  !>        15          14
  !>
  !----------------------------------------------------------------------  

  subroutine direct_solver_pastix_set_multithread(direct_solver)

    type(direct_solver_typ), intent(inout) :: direct_solver
#ifdef PASTIX
    integer(ip)                            :: topo_n,topo_c,topo_t,topo_p
    integer(ip)                            :: topo_bind,i
    integer(ip)                            :: gbrank,gbnp,p_per_n,my_core_beg
    integer(4)                             :: istat4


    topo_bind = 1
    !topo_n    = par_topo_num_nodes
    !topo_c    = par_topo_num_cores_per_node
    topo_t    = direct_solver % num_threads
    topo_p    = npart

    direct_solver % iparm(IPARM_THREAD_NBR) = topo_t

    if( topo_bind == 0 ) then
       !
       ! No binding
       !
       direct_solver % iparm(IPARM_BINDTHRD)  = API_BIND_AUTO

    else if( topo_bind == 2 ) then
       !
       ! Grouped binding
       !
       direct_solver % iparm(IPARM_BINDTHRD)  = API_BIND_AUTO

    else if( topo_bind == 1 ) then
       !
       ! User specify everything case
       !
       ! set the options
       direct_solver % iparm(IPARM_BINDTHRD)  = API_BIND_TAB

       ! obtain MPI informations   
       gbrank = kfl_paral
       gbnp   = npart

       !Call MPI_Comm_rank(MPI_COMM_WORLD, gbrank, info )
       !Call MPI_Comm_size(MPI_COMM_WORLD, gbnp, info )

       ! Compute the bindtab
       p_per_n     = ( topo_p + topo_n - 1 ) / topo_n
       my_core_beg =  mod(gbrank, p_per_n) * (topo_c / p_per_n )

       allocate( direct_solver % bindtab(topo_t),STAT=istat4 )
       if(istat4 < 0) call runend('BINDTAB FAILED')

       do i= 1,topo_t
          direct_solver % bindtab(i) = my_core_beg + i-1
       end do

       !call DMPH_ARITH_pastix_fortran_bindthreads &
       !     (direct_solver % pastix_data, topo_t, direct_solver % bindtab)
       !call DMPH_pastix_check_status(id_pastix, __LINE__, istat4)
       !if( istat4 < 0 ) call runend('COULD NOT BIND THREADS FOR PASTIX')

    end if

#endif

  end subroutine direct_solver_pastix_set_multithread

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-03-20
  !> @brief   Initialize solver type
  !> @details Initialize solver type
  !> 
  !-----------------------------------------------------------------------

  subroutine direct_solver_type_initialization_1(direct_solver)

    type(direct_solver_typ), pointer, intent(inout) :: direct_solver(:)
    integer(ip) :: ii

    do ii = 1,size(direct_solver)
       call direct_solver_type_initialization_s(direct_solver(ii))
    end do

  end subroutine direct_solver_type_initialization_1
  
  subroutine direct_solver_type_initialization_s(direct_solver)

    type(direct_solver_typ), intent(inout) :: direct_solver 

    direct_solver % kfl_solver          = 0_ip      ! Solver type
    direct_solver % kfl_symmetric       = 0_ip      ! Symmetric matrix
    direct_solver % nn                  = 0_ip      ! # nodes
    direct_solver % nn_own              = 0_ip      ! # own nodes
    direct_solver % ndof                = 0_ip      ! # number dof per node
    direct_solver % nrhs                = 1_ip      ! # RHS
    direct_solver % nz                  = 0_ip      ! CSR or CSC: # non-zero edges
    direct_solver % nz_own              = 0_ip      ! CSR or CSC: # non-zero edges
    direct_solver % kfl_paral           = 0_ip      ! Sequential or parallel
    nullify(direct_solver % ia)                     ! CSR or CSC: Matrix graph
    nullify(direct_solver % ja)                     ! CSR or CSC: Matrix graph
    direct_solver % nskyl               = 0_ip      ! Skyline: matrix size
    nullify(direct_solver % idiag)                  ! Skyline: diagonal
    nullify(direct_solver % iskyl)                  ! Skyline: pointer
    direct_solver % memor               = 0_8       ! Memory counter
    direct_solver % cputi               = 0.0_rp    ! CPU time
    direct_solver % num_initializations = 0_ip      ! Number of initializations
    direct_solver % num_solutions       = 0_ip      ! Number of solves
    direct_solver % num_factorizations  = 0_ip      ! Number of factorizations
    direct_solver % name                = ' '       ! Name of the solver
    direct_solver % fillin              = 0_ip      ! Fill-in     
    nullify(direct_solver % aa)                     ! Matrix
    nullify(direct_solver % ia_new)                 ! Alya renumbered graph
    nullify(direct_solver % ja_new)                 ! Alya renumbered graph
    nullify(direct_solver % aa_skyl)                ! Skyline factorized matrix
    nullify(direct_solver % IL)                     ! IL
    nullify(direct_solver % JL)                     ! JL
    nullify(direct_solver % LN)                     ! LN
    nullify(direct_solver % IU)                     ! IU
    nullify(direct_solver % JU)                     ! JU
    nullify(direct_solver % UN)                     ! UN
    nullify(direct_solver % permr)                  ! For Pastix, pastix_int_t
    nullify(direct_solver % invpr)                  ! For Pastix, pastix_int_t     
    nullify(direct_solver % iparm)                  ! For Pastix, pastix_int_t, WSMP and PWSMP
    nullify(direct_solver % dparm)                  ! For Pastix, pastix_float_t,WSMP and PWSMP
    nullify(direct_solver % bindtab)                ! For Pastix, bindtable
    direct_solver % pastix_data         = 0_8       ! For Pastix, pastix_data_ptr_t
    direct_solver % pastix_comm         = 0_4       ! For Pastix, pastix_mpi_int                               
    direct_solver % ldb                 = 1_ip      ! ldb
    direct_solver % num_threads         = 1_ip      ! OMP number of threads
    nullify(direct_solver % avals)                  ! avals, contains the values different from zero in WSMP and PWSMP
    nullify(direct_solver % gatsca_mumps)           ! Gather MUMP

  end subroutine direct_solver_type_initialization_s

end module mod_direct_solver
!> @}
