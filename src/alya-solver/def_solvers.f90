!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




!-----------------------------------------------------------------------
!> @addtogroup Solvers
!> @{
!> @file    def_iterative_solvers.f90
!> @author  houzeaux
!> @date    2022-03-01
!> @brief   Solvers
!> @details Solvers. The different steps are:
!>
!>          1. alloca .................. Allocate and define solver
!>          2. init .................... Initialize all solver arrays
!>          3. dim ..................... Pass matrix dimensions to solvers
!>          4. input ................... Set solver inputs
!>          5. preconditioner_alloca ... Allocate preconditioner
!>          6. set ..................... Pass some arrays coming from Alya to define
!>                                       some local structures. For example: groups for
!>                                       the deflated; list of linelets for linelet
!>                                       preconditioner; symbolical factorization.
!>                                       Matrix is not needed.
!>
!>          5. setup ................... Set up the solver in terms of dimensions
!>                                       If the matrix size does not change.
!>          6. preprocess .............. Just before solving. RHS exchange, impose
!>                                       Dirichlet boundary conditions, etc.
!>          7. solve ................... Solve the system Ax=b
!>          8. postprocess ............. Perform some postprocess operations like
!>                                       computing reaction forces, output some
!>                                       solver info, etc.
!>         
!-----------------------------------------------------------------------

module def_solvers
  
  use def_kintyp_basic,                  only : ip,rp,lg
  use def_kintyp_comm,                   only : comm_data_par_basic
  use def_kintyp_solvers,                only : soldat
  use def_kintyp_solvers,                only : solout
  use def_mat,                           only : mat
  use def_mat_dia,                       only : mat_dia
  use mod_communications_point_to_point, only : PAR_SYMMETRIC_EXCHANGE
  use mod_communications_global,         only : PAR_SUM
  use mod_optional_argument,             only : optional_argument
  use mod_exchange,                      only : exchange_init
  use mod_exchange,                      only : exchange_add
  use mod_exchange,                      only : exchange_end
  use def_solver
  implicit none

  private
  
  type input_arrays
     integer(ip), contiguous, pointer :: fixno(:) ! Dirichlet fixity flag
     real(rp),    contiguous, pointer :: bvess(:) ! Dirichlet value
     real(rp),    contiguous, pointer :: bvnat(:) ! Natural value
   contains
     procedure,      pass :: init => init_arrays  ! Solver input data
  end type input_arrays
  
  type input_data
     integer(ip)          :: kfl_algso          ! Solver tag
     integer(ip)          :: nkryd              ! Krylov dimension
     integer(ip)          :: kfl_exp_method     ! Method for explicit solver
     integer(ip)          :: kfl_ortho          ! Orthogonolization GMRES (0=classical Gram-Schmidt,1=modified)
     integer(ip)          :: kfl_recov          ! Recover original residual after slave exchange
     integer(ip)          :: kfl_exchange       ! If RHS should be exchanged by the solver
     integer(ip)          :: kfl_limit          ! Algebraic limiter
     integer(ip)          :: kfl_adres          ! Adaptive residual (0=no,1=yes)
     integer(ip)          :: nrhss              ! # of simultaneous RHS
     integer(ip)          :: kfl_coarse         ! If a coarse solver is to be used
     integer(ip)          :: kfl_penal          ! Penalization: method
     real(rp)             :: penal              ! Penalization: parameter
     real(rp)             :: relax              ! Relaxation parameter
     integer(ip)          :: kfl_force          ! Force continuity of solution across interface after solver
     integer(ip)          :: kfl_where          ! Unknown is on node, edge, element center
     integer(ip)          :: num_multiple_solves! Number of multiple solves using the same matrix
     integer(ip)          :: kfl_save_krylov    ! If Krylov subspace should be saved
     integer(ip)          :: kfl_roe_correction ! Round off error corrections should be used
     real(rp)             :: threshold          ! Threshold used for direct solver base preconditioners
     integer(ip)          :: kfl_iffix          ! If solver takes care of imposing Dirichlet b.c.
     integer(ip)          :: kfl_bvnat          ! If natural b.c. are imposed in solver
     integer(ip)          :: kfl_dirichlet      ! Way dirichlet b.c. are imposed
     integer(ip)          :: kfl_exp_order      ! Order for explicit solver
     !
     ! Matrix type
     !
     integer(ip)          :: kfl_format         ! CSR (1), COO (2), ELL (3)
     integer(ip)          :: kfl_scpre          ! Schur preconditioner solver
     integer(ip)          :: kfl_scaii          ! Schur matrix solver
     integer(ip)          :: kfl_cmplx          ! Real or complex solver JELENA
     integer(ip)          :: kfl_full_rows      ! Partial or full row matrix
     integer(ip)          :: kfl_symm           ! Symmetric matrix
     !
     ! Output
     !
     real(rp)             :: reaction_level     ! Level for reaction force
     integer(ip)          :: kfl_exres          ! Exact residual output ||b-Ax||/||b||
     integer(ip)          :: kfl_marke          ! Output of matrix in market format
     integer(ip)          :: kfl_kappa          ! If condition number should be computed
     integer(ip)          :: verbose            ! Verbose level
     !
     ! Convergence
     !
     integer(ip)          :: miter              ! Maximum number of iterations
     real(rp)             :: solco              ! Solver tolerance
     real(rp)             :: adres              ! Adaptive residual tolerance (0=no,1=yes)
     real(rp)             :: solmi              ! Minimum solver tolerance (if adaptive tolerance)
     real(rp)             :: bnorm_min          ! Minimum norm of RHS to consider null system
     !
     ! Units and output
     !
     integer(ip)          :: kfl_cvgso          ! Convergence output
     integer(ip)          :: lun_cvgso          ! Convergence file unit
     integer(ip)          :: kfl_solve          ! Solver info
     integer(ip)          :: lun_solve          ! Solver info file unit
     character(50)        :: conf_file          ! Configuration file
     !
     ! OpenMP
     !
     integer(ip)          :: omp_schedule       ! OMP schedule
     integer(ip)          :: omp_chunk_size     ! OMP chunk size
     integer(ip)          :: omp_interface      ! OMP for SpMV on interface
     !
     ! Preconditioning
     !
     integer(ip)          :: kfl_preco          ! Preconditioner tag
     integer(ip)          :: kfl_clean_precond  ! If preconditioner/coarse solver should be always recomputed   
     integer(ip)          :: kfl_leftr          ! Left(0) or right(1)
     !
     ! RAS
     !
     integer(ip)          :: kfl_block_ras      ! Full or block RAS
     integer(ip)          :: kfl_add_schwarz    ! RAS or additive schwarz 
     integer(ip)          :: max_dof_ras        ! Maximum number of dofs of RAS
     integer(ip)          :: num_subd_ras       ! Number of RAS subdomains
     !
     ! Block Gauss-Seidel
     !
     integer(ip)          :: kfl_blogs          ! Block Gauss-Seidel treatment
     integer(ip)          :: nblok              ! Number of blocks
     integer(ip)          :: kfl_renumbered_gs  ! Renumbered Gauss-Seidel permutation
     !
     ! Deflation
     !
     integer(ip)          :: kfl_defas          ! Deflated CG: Assembly of deflated
     integer(ip)          :: kfl_defso          ! Deflated CG: Solver of deflated
     integer(ip)          :: ngrou              ! Number of groups
     !
     ! Direct solvers
     !
     integer(ip)          :: levels_sparsity    ! Number of sparsity levels
   contains
     procedure,      pass :: init => init_input ! Solver input data
  end type input_data

  type :: output_data
     integer(ip)          :: iters               ! # iterations
     integer(ip)          :: itsol(3)            ! Solver statistics
     integer(ip)          :: nsolv               ! # solves
     real(rp)             :: toler               ! Dimensional tolerance x ||b||
     real(rp)             :: resid_init          ! Initial preconditioned residual
     real(rp)             :: resid_final         ! Final preconditioned residual
     real(rp)             :: resip_init          ! Initial preconditioned residual
     real(rp)             :: resip_final         ! Final preconditioned residual
     real(rp)             :: resip_old           ! Last before last preconiditioned residual
     real(rp)             :: cputi(10)           ! CPU time
     real(rp)             :: time_init           ! Initial time
     real(rp)             :: time_final          ! Final time
     real(rp)             :: time_old            ! Old time
     integer(ip)          :: num_spmv            ! Number of SpMV
     integer(ip)          :: num_dot             ! Number of dot products
     real(rp)             :: cpu_spmv(9)         ! CPU time of SMVP
     real(rp)             :: cpu_dot(9)          ! CPU time of dot product
     real(rp)             :: cpu_schur(10)       ! Schur complement solver timing
     real(rp)             :: xorth               ! Orthogonality check = (p^i-1,p^i)
     real(rp)             :: xdiag               ! Coefficient for diagonal solver
     integer(ip)          :: num_solves          ! Number of solves for multiple solves
     real(rp)             :: bnorm               ! RHS L2 norm
     real(rp)             :: bnorp               ! preconditioned RHS L2 norm
     real(rp)             :: bnorm_inv           ! Inverse RHS L2 norm
     real(rp)             :: bnorp_inv           ! Inverse preconditioned RHS L2 norm
     real(rp)             :: kappa               ! Condition number of the matrix
     real(rp)             :: lambda_min          ! Minimum eigenvalue
     real(rp)             :: lambda_max          ! maximum eigenvalue     
   contains
     procedure,      pass :: init => init_output ! Solver input data
  end type output_data

  type, abstract :: solver
     type(input_arrays)                                    :: arrays
     type(input_data)                                      :: input
     type(output_data)                                     :: output
     type(comm_data_par_basic), pointer                    :: comm
     integer(8)                                            :: memor(2)
     logical(lg)                                           :: iwrite
     integer(ip)                                           :: n                       ! # Number of block unknowns 
     integer(ip)                                           :: nown                    ! # Number of own block unknowns 
     integer(ip)                                           :: nint                    ! # Number of interior unknown 
     integer(ip)                                           :: nn                      ! # Number of unknowns = n * ndof 
     integer(ip)                                           :: ndof                    ! # Number of unknowns
     logical(lg)                                           :: comm_is_allocated       ! If matrix is allocated or points
   contains
     procedure (my_solver_init),            pass, deferred :: init                    ! Initialize
     procedure (my_solver_solve),           pass, deferred :: solve                   ! Solve
     procedure (my_solver_setup),           pass, deferred :: setup                   ! Setup
     procedure (my_solver_deallo),          pass, deferred :: deallo                  ! Deallocate
     procedure,                             pass           :: dim                     ! Set solver dimensions
     procedure,                             pass           :: dof                     ! dof 
     procedure,                             pass           :: exchange                ! Broadcast data  
     procedure,                             pass           :: info                    ! Check output
     procedure,                             pass           :: parallel_mv             ! Parallel mv product
     procedure,                             pass           :: parallel_sum_1          ! Parallel sum
     procedure,                             pass           :: parallel_sum_3          ! Parallel sum
     procedure,                             pass           :: parallel_residual       ! Parallel residual
     procedure,                             pass           :: parallel_dot            ! Parallel dot product
     procedure,                             pass           :: parallel_dots           ! Parallel dot product
     procedure,                             pass           :: parallel_L2norm         ! Parallel L2norm
     procedure,                             pass           :: parallel_exchange_rp    ! Parallel exchange
     procedure,                             pass           :: parallel_exchange_mat   ! Parallel exchange
     generic                                               :: parallel_exchange =>  & ! Sparse matrix-vector product    
          &                                                   parallel_exchange_rp, &
          &                                                   parallel_exchange_mat
     generic                                               :: parallel_sum =>       &
          &                                                   parallel_sum_1,       &
          &                                                   parallel_sum_3
  end type solver

    abstract interface
          
     subroutine my_solver_init(self)
       import                                              :: solver
       class(solver),                        intent(inout) :: self
     end subroutine my_solver_init
     
     subroutine my_solver_solve(self,b,x,a)
       import                                                        :: solver
       import                                                        :: mat
       import                                                        :: rp
       class(solver),                                  intent(inout) :: self 
       real(rp),                              pointer, intent(in)    :: b(:)
       real(rp),                              pointer, intent(inout) :: x(:)
       class(mat),                                     intent(in)    :: a       
     end subroutine my_solver_solve
     
     subroutine my_solver_setup(self,a)
       import                                                        :: solver
       import                                                        :: mat
       class(solver),                                  intent(inout) :: self 
       class(mat),                                     intent(inout) :: a       
     end subroutine my_solver_setup
     
     subroutine my_solver_deallo(self)
       import                                                        :: solver
       class(solver),                                  intent(inout) :: self
     end subroutine my_solver_deallo

  end interface


  public :: solver

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2022-03-04
  !> @brief   Input arrays initialization
  !> @details Input arrays initialization
  !> 
  !-----------------------------------------------------------------------
  
  subroutine init_arrays(self)
    class(input_arrays), intent(inout) :: self

    nullify(self % fixno)
    nullify(self % bvess)
    nullify(self % bvnat)

  end subroutine init_arrays
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2022-03-04
  !> @brief   Input initialization
  !> @details Input initialization
  !> 
  !-----------------------------------------------------------------------
  
  subroutine init_input(self)
    class(input_data), intent(inout) :: self

    self % kfl_algso           = -999                      ! No selfr: do not modify (used to check errors)
    self % kfl_recov           = 0                         ! Recover original residual after slave exchange
    self % kfl_exchange        = 1                         ! Do Exchange RHS in the selfr 
    self % kfl_ortho           = 1                         ! Orthogonalization GMRES: Gramm-Schmidt=0, modified=1
    self % kfl_limit           = 0                         ! Algebraic limiter
    self % kfl_adres           = 0                         ! Adaptive residual
    self % nrhss               = 1                         ! # of simultaneous RHS
    self % miter               = 1                         ! Max. number of iterations: should be 1, for Richardson selfrs
    self % nkryd               = 10                        ! Krylov dimensions
    self % kfl_coarse          = 0                         ! No coarse selfr
    self % solco               = 1.0e-6_rp                 ! Selfr tolerance
    self % adres               = 0.1_rp                    ! Adaptive residual tolerance
    self % solmi               = 1.0e-6_rp                 ! Minimum selfr tolerance
    self % kfl_cvgso           = 0                         ! Convergence flag
    self % lun_cvgso           = 0                         ! Convergence unit
    self % kfl_solve           = 0                         ! Output flag
    self % lun_solve           = 0                         ! Output unit
    self % kfl_penal           = 0                         ! Penalization: method
    self % penal               = 1.0_rp                    ! Penalization: parameter
    self % relax               = 1.0_rp                    ! Relaxation parameter
    self % kfl_preco           = 0                         ! Preconditioner
    self % kfl_leftr           = SOL_LEFT_PRECONDITIONING  ! Left preconditioner
    self % kfl_marke           = 0                         ! Output in market format
    self % kfl_force           = 0                         ! Do not force selfr continuity
    self % kfl_format          = SOL_CSR_FORMAT            ! Matrix format
    self % kfl_where           = SOL_NODES                 ! Where are the unknowns
    self % kfl_clean_precond   = 1                         ! Preconditioner/coarse selfr is not always recomputed at each self

    self % kfl_block_ras       = 0                         ! Full RAS
    self % max_dof_ras         = 0                         ! Maximum number of dofs of RAS
    self % num_subd_ras        = 1                         ! Number of RAS subdomains
    self % kfl_add_schwarz     = 0                         ! RAS instead of AS

    self % num_multiple_solves = 1                         ! Multiple selfs with the same matrix
    self % kfl_full_rows       = 0                         ! By default, matrix are partially assembled
    self % kfl_symm            = 0                         ! Unsymmetric matrix
    self % kfl_save_krylov     = 0                         ! Krylov subspace saving
    self % kfl_exres           = 0                         ! Dont output preconditioned residual
    self % kfl_kappa           = 0                         ! Condition number should not be computed
    self % verbose             = 0                         ! Verbose level
    self % conf_file           = ''                        ! Configuration file for selfrs
    self % bnorm_min           = 1.0e-12_rp                ! Minimum norm of RHS to consider null system
    self % kfl_roe_correction  = 0                         ! No round off error corrections (pipeline DCG)
    self % threshold           = -1.0_rp                   ! Threshold for direct selfr
    self % omp_schedule        = SOL_OMP_STATIC            ! OMP schedule
    self % omp_chunk_size      = 1000                      ! OMP chunk size
    self % omp_interface       = 0                         ! Do not use openMP for interface > 0 = chunk
    self % kfl_iffix           = 0                         ! Fixity not imposed by solver
    self % kfl_bvnat           = 0                         ! Natural condition imposed by solver
    self % kfl_dirichlet       = 1                         ! Way to impose Dirichlet boundary conditions 
    self % kfl_exp_method      = -999                      ! Lumped method
    self % kfl_exp_order       = 1                         ! Way to impose Dirichlet boundary conditions 
    !
    ! Deflated CG
    !
    self % kfl_defas           = SOL_CSR                   ! Deflated CG: Assembly (CSR, DENSE, SKYLINE)
    self % kfl_defso           = 0                         ! Deflated CG: Solver of deflated
    self % ngrou               = 0                         ! Deflated CG: # groups

    self % levels_sparsity     = 1                         ! Number of sparsity levels
  end subroutine init_input

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2022-03-04
  !> @brief   Output initialization
  !> @details Output initialization
  !> 
  !-----------------------------------------------------------------------
  
  subroutine init_output(self)
    class(output_data), intent(inout) :: self

    self % iters       = 0               ! # iterations
    self % itsol       = 0               ! Solver statistics
    self % nsolv       = 0               ! # solves
    self % toler       = 0.0_rp          ! Dimensional tolerance
    self % resid_init  = 0.0_rp          ! Initial residual
    self % resid_final = 0.0_rp          ! Final residual
    self % resip_init  = 0.0_rp          ! Initial preconditioned residual
    self % resip_final = 0.0_rp          ! Final preconditioned residual
    self % resip_old   = 0.0_rp          ! Last before last residual preconditioned residual
    self % time_init   = 0.0_rp          ! Initial time
    self % time_final  = 0.0_rp          ! Final time
    self % time_old    = 0.0_rp          ! Old time
    self % cputi       = 0.0_rp          ! CPU time
    self % num_spmv    = 0               ! Number of SpMV
    self % num_dot     = 0               ! Number of dot products
    self % cpu_spmv    = 0.0_rp          ! CPU time of SMVP
    self % cpu_dot     = 0.0_rp          ! CPU time of dot product
    self % cpu_schur   = 0.0_rp          ! Schur complement solver timing
    self % xorth       = 0.0_rp          ! Orthogonality check = (p^i-1,p^i)
    self % xdiag       = 0.0_rp          ! Coefficient for diagonal solver
    self % num_solves  = 0               ! Number of solves for multiple solves
    self % bnorm       = 0.0_rp          ! RHS L2 norm
    self % bnorp       = 0.0_rp          ! RHS L2 norm
    self % bnorm_inv   = 0.0_rp          ! Inverse RHS L2 norm
    self % bnorp_inv   = 0.0_rp          ! Inverse preconditioned RHS L2 norm
    self % kappa       = 0.0_rp          ! Condition number of the matrix
    self % lambda_min  = 0.0_rp          ! Minimum eigenvalue
    self % lambda_max  = 0.0_rp          ! maximum eigenvalue     
    
  end subroutine init_output

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-05
  !> @brief   Parallel exchange
  !> @details Parallel exchange
  !> 
  !-----------------------------------------------------------------------
  
  subroutine parallel_exchange_rp(self,x) 

    class(solver),                      intent(inout) :: self
    real(rp),                            pointer, intent(inout) :: x(:) !< Multiplicand

    if( associated(self % comm) .and. associated(x) ) then
       call PAR_SYMMETRIC_EXCHANGE(x,'SUM',self % comm,NDOF=self % ndof)
    end if
      
  end subroutine parallel_exchange_rp

  subroutine parallel_exchange_mat(self,a)

    class(solver),                                intent(inout) :: self
    class(mat),                                   intent(inout) :: a
    
    if( associated(self % comm) ) then !.and. associated(a) ) then
       select type ( a )
       class is ( mat_dia )
          if( a % nrows > 0 ) call PAR_SYMMETRIC_EXCHANGE(a % vA,'SUM',self % comm)
       end select
    end if

  end subroutine parallel_exchange_mat

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-05
  !> @brief   Parallel MV product
  !> @details Parallel MV product y = A x
  !> 
  !-----------------------------------------------------------------------
  
  subroutine parallel_residual(&
       self,a,x,b,r,&
       ON_INTERFACE,&
       MPI,OPENMP,INITIALIZATION,TIMING,&
       MY_TIMING,ASYNCHRONISM)

    class(solver),         intent(inout) :: self
    class(mat),                      intent(in)    :: a                !< Matrix type
    real(rp),               pointer, intent(in)    :: x(:)             !< Multiplicand
    real(rp),               pointer, intent(in)    :: b(:)             !< Result vector
    real(rp),               pointer, intent(inout) :: r(:)             !< Result vector
    logical(lg),  optional,          intent(in)    :: ON_INTERFACE     !< Parallel residual only on interface
    logical(lg),  optional,          intent(in)    :: MPI              !< If MPI should be used or not
    logical(lg),  optional,          intent(in)    :: OPENMP           !< If OpenMP should be used or not
    logical(lg),  optional,          intent(in)    :: INITIALIZATION   !< If result vector should be initialized
    logical(lg),  optional,          intent(in)    :: TIMING           !< If timing should be activated 
    real(rp),     optional,          intent(out)   :: MY_TIMING(2)     !< Register timing in this array
    logical(lg),  optional,          intent(in)    :: ASYNCHRONISM     !< Blocking or non-blocking communications
    integer(ip)                                    :: i

    call self % parallel_mv(a,x,r,ON_INTERFACE)

    do i = 1,self % nn
       r(i) = b(i) - r(i)
    end do

  end subroutine parallel_residual
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-05
  !> @brief   Parallel MV product
  !> @details Parallel MV product y = A x
  !> 
  !-----------------------------------------------------------------------
  
  subroutine parallel_mv(&
       self,a,x,y,&
       ON_INTERFACE,&
       MPI,OPENMP,INITIALIZATION,TIMING,&
       MY_TIMING,ASYNCHRONISM)

    class(solver),         intent(inout) :: self
    class(mat),                      intent(in)    :: a                !< Matrix type
    real(rp),               pointer, intent(in)    :: x(:)             !< Multiplicand
    real(rp),               pointer, intent(inout) :: y(:)             !< Result vector
    logical(lg),  optional,          intent(in)    :: ON_INTERFACE     !< MV product only on interface
    logical(lg),  optional,          intent(in)    :: MPI              !< If MPI should be used or not
    logical(lg),  optional,          intent(in)    :: OPENMP           !< If OpenMP should be used or not
    logical(lg),  optional,          intent(in)    :: INITIALIZATION   !< If result vector should be initialized
    logical(lg),  optional,          intent(in)    :: TIMING           !< If timing should be activated 
    real(rp),     optional,          intent(out)   :: MY_TIMING(2)     !< Register timing in this array
    logical(lg),  optional,          intent(in)    :: ASYNCHRONISM     !< Blocking or non-blocking communications
    integer(ip)                                    :: n1,n2,n3,n4
    
    if( associated(x) .and. associated(y) ) then

       if( associated(self % comm) .and. optional_argument(.false.,ON_INTERFACE) ) then
          n3     = self % nint+1
          n4     = self % n
          call a % mv(x,y,n3,n4,INITIALIZATION,OPENMP,self % input % omp_chunk_size,self % input % omp_schedule)
          call PAR_SYMMETRIC_EXCHANGE(y,'SUM',self % comm,NDOF=self % ndof)
          
       else if( associated(self % comm) .and. optional_argument(.true.,MPI) ) then
          n1     = 1
          n2     = self % nint
          n3     = self % nint + 1
          n4     = self % n
          call a % mv(x,y,n3,n4,INITIALIZATION,OPENMP,self % input % omp_chunk_size,self % input % omp_schedule)
          call PAR_SYMMETRIC_EXCHANGE(y,'SUM',self % comm,NDOF=self % ndof)
          call a % mv(x,y,n1,n2,INITIALIZATION,OPENMP,self % input % omp_chunk_size,self % input % omp_schedule)
          
       else
          n1     = 1
          n2     = self % n
          call a % mv(x,y,n1,n2,INITIALIZATION,OPENMP,self % input % omp_chunk_size,self % input % omp_schedule)
          
       end if
    end if

  end subroutine parallel_mv

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    07/10/2014
  !> @brief   Scalar product
  !> @details Compute the parallel scalar product XXDOTYY = XX.YY
  !>
  !----------------------------------------------------------------------

  function parallel_dot(self,x,y,MY_TIMING,OPENMP,MPI) result(xdoty)

    class(solver),                             intent(inout) :: self
    real(rp),                         pointer, intent(in)    :: x(:)             !< Multiplicand
    real(rp),                         pointer, intent(in)    :: y(:)             !< Result vector
    real(rp),               optional,          intent(out)   :: MY_TIMING(2)
    logical(lg),            optional,          intent(in)    :: OPENMP
    logical(lg),            optional,          intent(in)    :: MPI
    real(rp)                                                 :: xdoty    
    real(rp)                                                 :: dot    
    integer(ip)                                              :: ii,nn
    real(rp)                                                 :: time1,time2,time3
    logical(lg)                                              :: use_openmp
    logical(lg)                                              :: use_mpi

#ifdef BLAS
    real(rp) :: DDOT
    external DDOT
#endif
    call cputim(time1)

    dot        = 0.0_rp
    use_openmp = optional_argument(.false.,OPENMP)
    use_mpi    = optional_argument(.true. ,MPI)
    nn         = self % nown * self % ndof
    
    if( associated(x) .and. associated(y) .and. nn > 0 ) then

#ifdef BLAS
       dot = DDOT(nn,x,1_ip,y,1_ip)
#else
       if( use_openmp ) then
          !$OMP PARALLEL    DO SCHEDULE (STATIC) &
          !$OMP DEFAULT   ( NONE )               &
          !$OMP PRIVATE   ( ii )                 &
          !$OMP SHARED    ( x, y, nn )           &
          !$OMP REDUCTION ( +:dot )  
          do ii = 1,nn
             dot = dot + x(ii) * y(ii)
          end do
          !$OMP END PARALLEL DO
       else
          do ii = 1,nn
             dot = dot + x(ii) * y(ii)
          end do
       end if
#endif

    end if
    
    call cputim(time2)
    if( associated(self % comm) .and. use_mpi ) &
         call PAR_SUM(dot,COMM=self % comm,INCLUDE_ROOT=.true.)

    xdoty = dot
    
    call cputim(time3)
    if( present(MY_TIMING) ) then
       MY_TIMING(1) = time2 - time1
       MY_TIMING(2) = time3 - time2
    else
       self % output % cpu_dot(1) = self % output % cpu_dot(1) + time2 - time1
       self % output % cpu_dot(2) = self % output % cpu_dot(2) + time3 - time2
       self % output % cpu_dot(3) = self % output % cpu_dot(1) + self % output % cpu_dot(2)
       self % output % num_dot    = self % output % num_dot + 1
    end if

  end function parallel_dot

  function parallel_dots(self,x,y,z,MY_TIMING,OPENMP,MPI) result(xdoty)

    class(solver),                             intent(inout) :: self
    real(rp),                         pointer, intent(in)    :: x(:)             !< x
    real(rp),                         pointer, intent(in)    :: y(:)             !< y
    real(rp),                         pointer, intent(in)    :: z(:)             !< z
    real(rp),               optional,          intent(out)   :: MY_TIMING(2)
    logical(lg),            optional,          intent(in)    :: OPENMP
    logical(lg),            optional,          intent(in)    :: MPI
    real(rp)                                                 :: xdoty(2)    
    real(rp)                                                 :: dot(2)    
    integer(ip)                                              :: ii,nn
    real(rp)                                                 :: time1,time2,time3
    logical(lg)                                              :: use_openmp
    logical(lg)                                              :: use_mpi

#ifdef BLAS
    real(rp) :: DDOT
    external DDOT
#endif
    call cputim(time1)

    dot        = 0.0_rp
    use_openmp = optional_argument(.false.,OPENMP)
    use_mpi    = optional_argument(.true. ,MPI)
    nn         = self % nown * self % ndof

    if( associated(x) .and. associated(y) .and. nn > 0 ) then

#ifdef BLAS
       dot = DDOT(nn,x,1_ip,y,1_ip)
#else
       if( use_openmp ) then
          !$OMP PARALLEL    DO SCHEDULE (STATIC) &
          !$OMP DEFAULT   ( NONE )               &
          !$OMP PRIVATE   ( ii )                 &
          !$OMP SHARED    ( x, y, z, nn )        &
          !$OMP REDUCTION ( +:dot )  
          do ii = 1,nn
             dot(1) = dot(1)  + x(ii) * y(ii)
             dot(2) = dot(2)  + x(ii) * z(ii)
          end do
          !$OMP END PARALLEL DO
       else
          do ii = 1,nn
             dot(1)  = dot(1)  + x(ii) * y(ii)
             dot(2)  = dot(2)  + x(ii) * z(ii)
           end do
       end if
#endif

    end if
    
    call cputim(time2)
    if( associated(self % comm) .and. use_mpi ) &
         call PAR_SUM(2_ip,dot,COMM=self % comm,INCLUDE_ROOT=.true.)

    xdoty = dot
    
    call cputim(time3)
    if( present(MY_TIMING) ) then
       MY_TIMING(1) = time2 - time1
       MY_TIMING(2) = time3 - time2
    else
       self % output % cpu_dot(1) = self % output % cpu_dot(1) + time2 - time1
       self % output % cpu_dot(2) = self % output % cpu_dot(2) + time3 - time2
       self % output % cpu_dot(3) = self % output % cpu_dot(1) + self % output % cpu_dot(2)
       self % output % num_dot    = self % output % num_dot + 1
    end if

  end function parallel_dots

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    07/10/2014
  !> @brief   Scalar product
  !> @details Compute the parallel scalar product XXDOTYY = XX.YY
  !>
  !----------------------------------------------------------------------

  function parallel_L2norm(self,x,MY_TIMING,OPENMP,MPI) result(xdotx)

    class(solver),                   intent(inout) :: self
    real(rp),                         pointer, intent(in)    :: x(:)             !< Multiplicand
    real(rp),               optional,          intent(out)   :: MY_TIMING(2)
    logical(lg),            optional,          intent(in)    :: OPENMP
    logical(lg),            optional,          intent(in)    :: MPI
    real(rp)                                                 :: xdotx

    xdotx = sqrt(parallel_dot(self,x,x,MY_TIMING,OPENMP,MPI))

  end function parallel_L2norm

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    07/10/2014
  !> @brief   Parallel sum
  !> @details Parallel sum
  !>
  !----------------------------------------------------------------------
  
  subroutine parallel_sum_1(self,x,MPI) 

    class(solver),                      intent(in)    :: self
    real(rp),                  pointer, intent(inout) :: x(:)             !< x
    logical(lg),     optional,          intent(in)    :: MPI 

    if( optional_argument(.true.,MPI) .and. associated(self % comm) ) then
       call PAR_SUM(x,COMM=self % comm ,INCLUDE_ROOT=.true.)
    end if

  end subroutine parallel_sum_1

  subroutine parallel_sum_3(self,x,MPI) 

    class(solver),                      intent(in)    :: self
    real(rp),                  pointer, intent(inout) :: x(:,:,:)         !< x
    logical(lg),     optional,          intent(in)    :: MPI 

    if( optional_argument(.true.,MPI) .and. associated(self % comm) ) then
       call PAR_SUM(x,COMM=self % comm ,INCLUDE_ROOT=.true.)
    end if

  end subroutine parallel_sum_3

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    07/10/2014
  !> @brief   dof
  !> @details dof
  !>
  !----------------------------------------------------------------------
  
  elemental integer(ip) function dof(self,i,idof)
    
    class(solver), intent(in) :: self
    integer(ip),   intent(in) :: i
    integer(ip),   intent(in) :: idof

    dof = (i-1)*self % ndof + idof

  end function dof

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-03
  !> @brief   Set dimensions
  !> @details Set dimensions
  !> 
  !-----------------------------------------------------------------------

  subroutine dim(self,a)

    class(solver),              intent(inout) :: self          !< Solver
    class(mat),                 intent(inout) :: a             !< Matrix A
    
    self % n    = a % nrows
    self % ndof = a % ndof1
    self % nn   = self % n * self % ndof
    
    if( associated(self % comm) ) then
       self % nint = self % comm % npoi1
       self % nown = self % comm % npoi3
    else
       self % nint = a % nrows
       self % nown = a % nrows
    end if
    
    if( associated(self % comm) ) then
       if( self % comm % RANK4 <= 0 ) then
          self % iwrite = .true.
       else
          self % iwrite = .false.
       end if
    else
       self % iwrite = .true.
    end if

  end subroutine dim

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-03
  !> @brief   Broadacast data
  !> @details Broadacast data from master to slaves
  !> 
  !-----------------------------------------------------------------------

  subroutine exchange(self)

    class(solver),              intent(inout) :: self          !< Solver
    !
    ! Initialize
    !
    call exchange_init()
    !
    ! Accumulate
    !
    call exchange_add(self % input % kfl_algso)           
    call exchange_add(self % input % kfl_recov)           
    call exchange_add(self % input % kfl_exchange)        
    call exchange_add(self % input % kfl_ortho)           
    call exchange_add(self % input % kfl_limit)           
    call exchange_add(self % input % kfl_adres)           
    call exchange_add(self % input % nrhss)               
    call exchange_add(self % input % miter)               
    call exchange_add(self % input % nkryd)               
    call exchange_add(self % input % kfl_coarse)          
    call exchange_add(self % input % solco)               
    call exchange_add(self % input % adres)               
    call exchange_add(self % input % solmi)               
    call exchange_add(self % input % kfl_cvgso)           
    call exchange_add(self % input % lun_cvgso)           
    call exchange_add(self % input % kfl_solve)           
    call exchange_add(self % input % lun_solve)           
    call exchange_add(self % input % kfl_penal)           
    call exchange_add(self % input % penal)               
    call exchange_add(self % input % relax)               
    call exchange_add(self % input % kfl_preco)            
    call exchange_add(self % input % kfl_leftr)            
    call exchange_add(self % input % kfl_marke)            
    call exchange_add(self % input % kfl_force)            
    call exchange_add(self % input % kfl_format)           
    call exchange_add(self % input % kfl_where)            
    call exchange_add(self % input % kfl_clean_precond)    

    call exchange_add(self % input % kfl_block_ras)        
    call exchange_add(self % input % max_dof_ras)          
    call exchange_add(self % input % num_subd_ras)         
    call exchange_add(self % input % kfl_add_schwarz)      

    call exchange_add(self % input % num_multiple_solves)  
    call exchange_add(self % input % kfl_full_rows)        
    call exchange_add(self % input % kfl_symm)             
    call exchange_add(self % input % kfl_save_krylov)      
    call exchange_add(self % input % kfl_exres)            
    call exchange_add(self % input % kfl_kappa)            
    call exchange_add(self % input % verbose)            
    call exchange_add(self % input % conf_file)            
    call exchange_add(self % input % bnorm_min)            
    call exchange_add(self % input % kfl_roe_correction)   
    call exchange_add(self % input % threshold)            
    call exchange_add(self % input % omp_schedule)         
    call exchange_add(self % input % omp_chunk_size)       
    call exchange_add(self % input % omp_interface)        
    call exchange_add(self % input % kfl_iffix)            
    call exchange_add(self % input % kfl_bvnat)            
    call exchange_add(self % input % kfl_dirichlet)        
    call exchange_add(self % input % kfl_exp_method)       
    call exchange_add(self % input % kfl_exp_order) 
    
    call exchange_add(self % input % kfl_defas)         
    call exchange_add(self % input % kfl_defso)            
    call exchange_add(self % input % ngrou)                

    call exchange_add(self % input % levels_sparsity)           
    !
    ! Exchange
    !
    call exchange_end()           

  end subroutine exchange

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-04-01
  !> @brief   Should we write
  !> @details Should we write
  !> 
  !-----------------------------------------------------------------------

  logical(lg) pure function info(self,level) result(info_out)
    
    class(solver), intent(in) :: self  !< Solver
    integer(ip),   intent(in) :: level !< Level of message

    if( self % iwrite .and. level <= self % input % verbose ) then
       info_out = .true.
    else
       info_out = .false.
    end if
    
  end function info
  
end module def_solvers
!> @}
