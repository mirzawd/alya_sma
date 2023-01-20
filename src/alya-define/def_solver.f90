!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module def_solver
  !-----------------------------------------------------------------------
  !****f* defmod/def_solver
  ! NAME
  !   def_solver
  ! DESCRIPTION
  !    Heading for the solvers
  !***
  !-----------------------------------------------------------------------
  use def_kintyp_basic,   only : ip,rp,lg
  use def_kintyp_solvers, only : soltyp,eigtyp
  use def_kintyp_solvers, only : direct_solver_typ
  !
  ! Does and dont's
  !
  integer(ip),   parameter :: SOL_YES                          =  0
  integer(ip),   parameter :: SOL_NO                           =  1  
  !
  ! Parameters solvers
  !
  integer(ip),   parameter :: SOL_SOLVER_DIRECT                =   0
  integer(ip),   parameter :: SOL_SOLVER_LDU                   =   0
  integer(ip),   parameter :: SOL_SOLVER_MUMPS                 =  -1
  integer(ip),   parameter :: SOL_SOLVER_PASTIX                =  -3
  integer(ip),   parameter :: SOL_SOLVER_WSMP                  =  -5
  integer(ip),   parameter :: SOL_SOLVER_PWSMP                 =  -7
  integer(ip),   parameter :: SOL_SOLVER_GPUQR                 = -100
  integer(ip),   parameter :: SOL_SOLVER_CG                    =   1
  integer(ip),   parameter :: SOL_SOLVER_DEFLATED_CG           =   2
  integer(ip),   parameter :: SOL_SOLVER_A_DEF2                =   3
  integer(ip),   parameter :: SOL_SOLVER_BICGSTAB              =   5
  integer(ip),   parameter :: SOL_SOLVER_GMRES                 =   8
  integer(ip),   parameter :: SOL_SOLVER_DIAGONAL              =   9
  integer(ip),   parameter :: SOL_SOLVER_RICHARDSON            =   9
  integer(ip),   parameter :: SOL_SOLVER_MATRIX_DIAGONAL       =  10
  integer(ip),   parameter :: SOL_SOLVER_MATRIX_RICHARDSON     =  10
  integer(ip),   parameter :: SOL_SOLVER_DEFLATED_BICGSTAB     =  12
  integer(ip),   parameter :: SOL_SOLVER_DEFLATED_GMRES        =  13
  integer(ip),   parameter :: SOL_SOLVER_SPARSE_DIRECT         =  14
  integer(ip),   parameter :: SOL_SOLVER_STEEPEST_DESCENT      =  15
  integer(ip),   parameter :: SOL_SOLVER_JACOBI                =  16
  integer(ip),   parameter :: SOL_SOLVER_PIPELINED_CG          =  18
  integer(ip),   parameter :: SOL_SOLVER_PIPELINED_DEFLATED_CG =  19
  integer(ip),   parameter :: SOL_SOLVER_MAPHYS_UNSYMMETRIC    =  20
  integer(ip),   parameter :: SOL_SOLVER_MAPHYS_SYMMETRIC      =  21
  integer(ip),   parameter :: SOL_SOLVER_MAPHYS_SPD            =  22
  integer(ip),   parameter :: SOL_SOLVER_AGMG                  =  23
  integer(ip),   parameter :: SOL_SOLVER_MINRES                =  24
  integer(ip),   parameter :: SOL_SOLVER_PSBLAS                =  25
  integer(ip),   parameter :: SOL_SOLVER_EXPLICIT              =  26
  integer(ip),   parameter :: SOL_SOLVER_POLYNOMIAL            =  27
  integer(ip),   parameter :: SOL_SOLVER_ORTHOMIN              =  28
  integer(ip),   parameter :: SOL_SOLVER_GAUSS_SEIDEL          =  29
  integer(ip),   parameter :: SOL_SOLVER_PETSC                 =  30
  integer(ip),   parameter :: SOL_SOLVER_SOR                   =  31
  integer(ip),   parameter :: SOL_SOLVER_SSOR                  =  32
  integer(ip),   parameter :: SOL_SOLVER_GAUSS_ELIMINATION     =  33
  integer(ip),   parameter :: SOL_SOLVER_APPROX_INVERSE        =  34
  integer(ip),   parameter :: SOL_SOLVER_LINELET               =  35
  integer(ip),   parameter :: SOL_SOLVER_RAS                   =  36
  !
  integer(ip),   parameter :: SOL_DEFLATED_CG                  =   2
  integer(ip),   parameter :: SOL_GCG                          = 100
  integer(ip),   parameter :: SOL_GGMRS                        = 101
  integer(ip),   parameter :: SOL_GDECG                        = 102
  integer(ip),   parameter :: SOL_GPCG                         = 103
  integer(ip),   parameter :: SOL_GAMGX                        = 104
  integer(ip),   parameter :: SOL_GCGNOPREC                    = 105
  integer(ip),   parameter :: SOL_NO_SOLVER                    = -999
  !
  ! Parameters preconditioners
  !
  integer(ip),   parameter :: SOL_NO_PRECOND                =  0
  integer(ip),   parameter :: SOL_SQUARE                    =  1
  integer(ip),   parameter :: SOL_DIAGONAL                  =  2
  integer(ip),   parameter :: SOL_MATRIX                    =  3
  integer(ip),   parameter :: SOL_LINELET                   =  4
  integer(ip),   parameter :: SOL_MASS_MATRIX               =  5
  integer(ip),   parameter :: SOL_GAUSS_SEIDEL              =  6
  integer(ip),   parameter :: SOL_SOR                       =  6
  integer(ip),   parameter :: SOL_CLOSE_MASS_MATRIX         =  7
  integer(ip),   parameter :: SOL_RICHARDSON                =  8
  integer(ip),   parameter :: SOL_ORTHOMIN                  =  9
  integer(ip),   parameter :: SOL_APPROXIMATE_SCHUR         = 10
  integer(ip),   parameter :: SOL_ABB                       = 11
  integer(ip),   parameter :: SOL_MOD_DIAGONAL              = 12
  integer(ip),   parameter :: SOL_AII                       = 13
  integer(ip),   parameter :: SOL_DEFLATED                  = 14
  integer(ip),   parameter :: SOL_MULTIGRID                 = 15
  integer(ip),   parameter :: SOL_COARSE_AII                = 16
  integer(ip),   parameter :: SOL_NEUMANN                   = 17
  integer(ip),   parameter :: SOL_LOCAL_DIAGONAL            = 18
  integer(ip),   parameter :: SOL_AUGMENTED_DIAGONAL        = 19
  integer(ip),   parameter :: SOL_RAS                       = 20
  integer(ip),   parameter :: SOL_RAS2                      = 21
  integer(ip),   parameter :: SOL_BLOCK_DIAGONAL            = 22
  integer(ip),   parameter :: SOL_LOCAL_MASS_MATRIX         = 23
  integer(ip),   parameter :: SOL_BIDIAGONAL                = 24
  integer(ip),   parameter :: SOL_SYMMETRIC_GAUSS_SEIDEL    = 26
  integer(ip),   parameter :: SOL_SSOR                      = 26
  integer(ip),   parameter :: SOL_APPROX_INVERSE            = 31
  !
  ! Gauss-Seidel options
  !
  integer(ip),   parameter :: SYMMETRIC_GAUSS_SEIDEL        = -1
  integer(ip),   parameter :: STREAMWISE_GAUSS_SEIDEL       =  1
  integer(ip),   parameter :: CLASSICAL_GAUSS_SEIDEL        =  0
  !
  ! Bidiagonal options
  !
  integer(ip),   parameter :: STREAMWISE_BIDIAGONAL         =  1
  integer(ip),   parameter :: BIDIAGONAL                    =  0
  !
  ! Other parameters
  !
  integer(ip),   parameter :: SOL_SKYLINE                  =   0
  integer(ip),   parameter :: SOL_CSR                       =  1
  integer(ip),   parameter :: SOL_DENSE                     =  2
  integer(ip),   parameter :: SOL_DIRECT                    =  1
  integer(ip),   parameter :: SOL_ITERATIVE                 =  0
  integer(ip),   parameter :: SOL_LEFT_PRECONDITIONING      =  0
  integer(ip),   parameter :: SOL_RIGHT_PRECONDITIONING     =  1
  integer(ip),   parameter :: SOL_LEFT_RIGHT_PRECONDITIONING =  2
  !
  ! Direct solver type
  !
  integer(ip),   parameter :: SOL_DIRECT_SOLVER_MUMPS        = -1
  integer(ip),   parameter :: SOL_DIRECT_SOLVER_ALYA         = -2
  integer(ip),   parameter :: SOL_DIRECT_SOLVER_PASTIX       = -3
  integer(ip),   parameter :: SOL_DIRECT_SOLVER_SKYLINE_ALYA = -4
  integer(ip),   parameter :: SOL_DIRECT_SOLVER_WSMP         = -5
  integer(ip),   parameter :: SOL_DIRECT_SOLVER_PWSMP        = -7
  integer(ip),   parameter :: SOL_DIRECT_SOLVER_GPUQR        = -100
  !
  ! Other parameters
  !
  integer(ip),   parameter :: SOL_BCSR_FORMAT                =  1
  integer(ip),   parameter :: SOL_CSR_FORMAT                 =  1
  integer(ip),   parameter :: SOL_COO_FORMAT                 =  2
  integer(ip),   parameter :: SOL_ELL_FORMAT                 =  3
  integer(ip),   parameter :: SOL_DEN_FORMAT                 =  4
  integer(ip),   parameter :: SOL_NODES                      =  0
  integer(ip),   parameter :: SOL_EDGES                      =  1
  integer(ip),   parameter :: SOL_ELEMENTS                   =  2
  integer(ip),   parameter :: SOL_MATRIX_HAS_CHANGED         =  0
  integer(ip),   parameter :: SOL_SYSTEM_HAS_BEEN_SOLVED     =  1
  integer(ip),   parameter :: SOL_EXPLICIT_LUMPED            =  0
  integer(ip),   parameter :: SOL_EXPLICIT_CONSISTENT        =  1
  integer(ip),   parameter :: SOL_EXPLICIT_APPROX_INVERSE    =  2  
  !
  ! Parallelization  
  !
  integer(ip),   parameter :: SOL_OMP_OFF                    =  0
  integer(ip),   parameter :: SOL_OMP_STATIC                 =  2
  integer(ip),   parameter :: SOL_OMP_GUIDED                 =  3
  integer(ip),   parameter :: SOL_OMP_DYNAMIC                =  4
  !   
  ! Direct solver based preconditioners
  !
  type(direct_solver_typ)     ::    &
       typ_direct_solver_RAS,       &
       typ_direct_solver_coarse,    &
       typ_direct_solver_block_LU,  &
       typ_direct_solver_AMG,       &
       typ_direct_solver_Deflation
  !
  ! Other parameters
  !
  integer(ip)          :: & 
       iters                ! # solver iterations
  integer(8)           :: &
       smemo(2),          & ! Eigenvalue solver memory counter
       memit(2),          & ! Iterative solver memory counter
       memdi(2),          & ! Direct solver memory counter
       memma(2)             ! All solver memory counter
  !
  ! Reals
  !
  real(rp) ::&
       cpu_solve,         & ! Algebraic solver
       cpu_eigen,         & ! Eigen solver
       resin,             & ! Initial residual norm
       resfi,             & ! Final residual norm
       resi2,             & ! Last residual norm
       resi1,             & ! Last before the last residual norm
       resal                ! Algebraic residual before enterring solver
  !
  ! Schur things
  !    
  integer(ip)              :: &
       nzdom_aii,             &      ! Aii
       nzdom_aib,             &      ! Aib
       nzdom_abi,             &      ! Abi
       nzdom_abb,             &      ! Abb
       nzdom_prec                    ! P
  !
  ! Pointers
  !
  integer(ip), allocatable :: &
       lpntn(:),&                       ! Pointer to the renumbering of the nodes
       lpont(:)                         ! Pointer to the skyline mesh graph matrix
  !
  ! Precoditioner
  !
  real(rp),    pointer     ::    &
       block_diagonal(:,:,:),    &      ! Block diagonal
       block_invdiagonal(:,:,:)         ! Inverse block diagonal

  type(soltyp), pointer :: solve_sol(:), solad_sol(:)
  type(eigtyp), pointer :: eigen_sol(:)

end module def_solver
