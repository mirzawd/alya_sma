!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kinds_and_types
!> @{
!> @file    def_kintyp_solvers.g90
!> @author  houzeaux
!> @date    2020-04-04
!> @brief   Kind and types for solvers
!> @details Solvers
!-----------------------------------------------------------------------

module def_kintyp_solvers
  
  use def_kintyp_basic, only : ip,rp,lg,i1p

#ifdef MAPHYS
  use dmph_maphys_mod
#endif

#ifdef MUMPS
  INCLUDE 'dmumps_struc.h'
#endif

#ifdef PETSC
  use def_kintyp_petsc, only : PETSC_LINEAR_SOLVER
#endif

#ifndef MAPHYS
  type DMPH_maphys_t
     integer(ip)          :: n 
     integer(ip)          :: nnz
     integer(ip)          :: sym
     integer(ip)          :: job
     integer(4)           :: comm
     real(rp),    pointer :: values(:)
     real(rp),    pointer :: rhs(:)
     real(rp),    pointer :: sol(:)
     real(rp),    pointer :: RCNTL(:)
     real(rp),    pointer :: RINFOG(:)
     real(rp),    pointer :: RINFO(:)
     integer(ip), pointer :: IINFOG(:)
     integer(ip), pointer :: ICNTL(:)
     integer(ip), pointer :: rows(:)
     integer(ip), pointer :: cols(:)
  end type DMPH_maphys_t
#endif

  type lpoin_block_typ
     type(block_matrix_typ), pointer :: block2_num(:,:)
     type(block_rhs_typ),    pointer :: block1_num(:)
  end type lpoin_block_typ
  type block_matrix_typ
     real(rp), pointer :: matrix(:,:,:)
  end type block_matrix_typ
  type block_rhs_typ
     real(rp),    pointer :: reaction(:)
     real(rp),    pointer :: rhs(:)
     real(rp),    pointer :: bvess(:,:)
     integer(ip), pointer :: kfl_fixno(:,:)     ! Fixity array
     real(rp),    pointer :: bvnat(:,:)         ! Neumann value
  end type block_rhs_typ
  type block_typ
     real(rp),    pointer :: bvess(:,:)
     integer(ip), pointer :: kfl_fixno(:,:)     ! Fixity array
     real(rp),    pointer :: bvnat(:,:)         ! Neumann value
  end type block_typ

#ifdef PASTIX
#include "pastix_fortran.h"
#else
!  #define PASTIX_INT_KIND    4
!  #define pastix_int_t       INTEGER(kind=4)
!  #define pastix_uint_t      unsigned INTEGER(kind=4)
!  #define pastix_data_ptr_t  INTEGER(kind=8)
!  #define MPI_PASTIX_INT     MPI_INTEGER4
!  #define pastix_float_t     REAL(kind=8)
!  #define MPI_PASTIX_FLOAT   MPI_REAL8
#endif

  type direct_solver_typ
     ! Input data
     integer(ip)          :: kfl_solver          ! Solver type
     integer(ip)          :: kfl_symmetric       ! Symmetric matrix
     integer(ip)          :: nn                  ! # nodes
     integer(ip)          :: nn_own              ! # own nodes
     integer(ip)          :: ndof                ! # number dof per node
     integer(ip)          :: nrhs                ! # RHS
     integer(ip)          :: nz                  ! CSR or CSC: # non-zero edges
     integer(ip)          :: nz_own              ! CSR or CSC: # non-zero edges
     integer(ip)          :: kfl_paral           ! Should be run in parallel
     integer(ip), pointer :: ia(:)               ! CSR or CSC: Matrix graph
     integer(ip), pointer :: ja(:)               ! CSR or CSC: Matrix graph
     integer(ip)          :: nskyl               ! Skyline: matrix size
     integer(ip), pointer :: idiag(:)            ! Skyline: diagonal
     integer(ip), pointer :: iskyl(:)            ! Skyline: pointer
     ! Output
     integer(8)           :: memor(2)            ! Memory counter
     real(rp)             :: cputi(3)            ! CPU time
     integer(ip)          :: num_initializations ! Number of initializations
     integer(ip)          :: num_solutions       ! Number of solves
     integer(ip)          :: num_factorizations  ! Number of factorizations
     character(30)        :: name                ! Name of the solver
     real(rp)             :: fillin              ! Fill-in
     ! Solver internal arrays
     real(rp),    pointer :: aa(:,:,:)           ! Matrix
     integer(ip), pointer :: ia_new(:)           ! Alya renumbered graph
     integer(ip), pointer :: ja_new(:)           ! Alya renumbered graph
     real(rp),    pointer :: aa_skyl(:)          ! Skyline factorized matrix
     integer(ip), pointer :: IL(:)
     integer(ip), pointer :: JL(:)
     real(rp),    pointer :: LN(:)
     integer(ip), pointer :: IU(:)
     integer(ip), pointer :: JU(:)
     real(rp),    pointer :: UN(:)
     integer(ip), pointer :: permr(:)            ! For Pastix, pastix_int_t
     integer(ip), pointer :: invpr(:)            ! For Pastix, pastix_int_t
     ! Pastix internal arrays
     integer(ip), pointer :: iparm(:)            ! For Pastix, pastix_int_t, WSMP and PWSMP
     real(rp),    pointer :: dparm(:)            ! For Pastix, pastix_float_t,WSMP and PWSMP
     integer(ip), pointer :: bindtab(:)
     integer(8)           :: pastix_data         ! For Pastix, pastix_data_ptr_t
     integer(4)           :: pastix_comm         ! For Pastix, pastix_mpi_int
     integer(ip)          :: ldb
     integer(ip)          :: num_threads
     real(rp),    pointer :: avals(:)
     ! MUMPS
     integer(ip), pointer :: gatsca_mumps(:)     ! Gather MUMP
#ifdef MUMPS
     TYPE (DMUMPS_STRUC) mumps_par
#endif
  end type direct_solver_typ

  !----------------------------------------------------------------------
  !
  ! Solver input data
  !
  !----------------------------------------------------------------------
  
  type soldat
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
     !
     ! Output
     !
     real(rp)             :: reaction_level     ! Level for reaction force
     integer(ip)          :: kfl_exres          ! Exact residual output ||b-Ax||/||b||
     integer(ip)          :: kfl_marke          ! Output of matrix in market format
     integer(ip)          :: kfl_kappa          ! If condition number should be computed
     !
     ! Convergence
     !
     integer(ip)          :: miter              ! Maximum number of iterations
     real(rp)             :: solco              ! Solver tolerance
     real(rp)             :: adres              ! Adaptive residual tolerance (0=no,1=yes)
     real(rp)             :: solmi              ! Minimum solver tolerance (if adaptive tolerance)
     real(rp)             :: bnorm_min          ! Minimum norm of RHS to consider null system
     integer(ip)          :: kfl_normalization  ! Residual normalization strategy
     real(rp)             :: normalization      ! Residual normalization value
     !
     ! Units and output
     !
     integer(ip)          :: kfl_cvgso          ! Convergence output
     integer(ip)          :: lun_cvgso          ! Convergence file unit
     integer(ip)          :: kfl_solve          ! Solver info
     integer(ip)          :: lun_solve          ! Solver info file unit
     integer(ip)          :: lun_exsol          ! External solver info file unit
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
     integer(ip)          :: kfl_linty          ! Linelet: Linlet type
     integer(ip)          :: kfl_defpr          ! Smoother preconditioner
     integer(ip)          :: itpre              ! Jacobi iterations
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
   contains
     procedure,      pass :: init => init_soldat ! Solver input data
  end type soldat

  !----------------------------------------------------------------------
  !
  ! Solver output
  !
  !----------------------------------------------------------------------
  
  type, extends(soldat) :: solout
     !
     ! Output
     !
     integer(ip)          :: iters              ! # iterations
     integer(ip)          :: itsol(3)           ! Solver statistics
     integer(ip)          :: nsolv              ! # solves
     real(rp)             :: resin              ! Initial preconditioned residual
     real(rp)             :: resfi              ! Final preconditioned residual
     real(rp)             :: resi1              ! Last residual = ||b-Ax||/||b||
     real(rp)             :: resi2              ! Last before last residual = ||b-Ax||/||b||
     real(rp)             :: resf2              ! Final non-preconditioned residual = ||b-Ax||/||b||
     real(rp)             :: cputi(10)          ! CPU time
     integer(ip)          :: num_spmv           ! Number of SpMV
     integer(ip)          :: num_dot            ! Number of dot products
     real(rp)             :: cpu_spmv(9)        ! CPU time of SMVP
     real(rp)             :: cpu_dot(9)         ! CPU time of dot product
     real(rp)             :: cpu_schur(10)      ! Schur complement solver timing
     real(rp)             :: xorth              ! Orthogonality check = (p^i-1,p^i)
     real(rp)             :: xdiag              ! Coefficient for diagonal solver
     integer(ip)          :: num_solves         ! Number of solves for multiple solves
     real(rp)             :: bnorm              ! RHS L2 norm
     real(rp)             :: kappa              ! Condition number of the matrix
     real(rp)             :: lambda_min         ! Minimum eigenvalue
     real(rp)             :: lambda_max         ! maximum eigenvalue     
   contains
     procedure,      pass :: init => init_solout ! Solver input data
  end type solout
  
  !----------------------------------------------------------------------
  !
  ! Solver type
  !
  !----------------------------------------------------------------------
  
  type, extends(solout) :: soltyp
     !
     ! Derived parameters
     !
     integer(8)           :: memor(2)           ! Memory counter
     integer(ip)          :: ndofn              ! # dof per nodes
     character(50)        :: wprob              ! Problem name
     character(50)        :: wsolv              ! Solver name
     character(50)        :: wprec              ! Preconditioner name
     integer(ip)          :: nzmat              ! System size
     integer(ip)          :: nzmat_own          ! System size for full row matrix
     integer(ip)          :: nzrhs              ! Size RHS = ndofn * nequa
     integer(ip)          :: nequa              ! # nodes to solve
     integer(ip)          :: nequa_own          ! # nodes to solve (own nodes for scalar product)
     integer(ip)          :: nequa_halo         ! # nodes to solve (own nodes for scalar product)
     integer(ip)          :: nunkn              ! # unknowns to solve nequa*ndofn
     integer(ip)          :: ncols              ! # columns (non square matrix)
     integer(ip)          :: ndof2              ! # dof^2
     integer(ip)          :: nequ1              ! Interior unknowns
     integer(ip)          :: nequ2              ! Start own boundary
     integer(ip)          :: nequ3              ! End own boundary
     integer(ip)          :: kfl_symme          ! Matrix symmetric assembly (0=no,1=yes)
     integer(ip)          :: kfl_symeq          ! Equation symmetry (0=no,1=yes)
     integer(ip)          :: kfl_assem          ! Matrix has been assembled
     integer(ip)          :: kfl_force_assembly ! Force matrix assembly
     integer(ip)          :: kfl_do_assembly    ! If assembly
     integer(ip)          :: kfl_schur          ! Schur solver (0=no)
     integer(ip)          :: kfl_schum          ! If matrices come from a Schur complement
     integer(ip)          :: kfl_version        ! Version of solver (new one incldues pre and post)
     integer(ip)          :: heade              ! Track if header has already been written in file
     character(150)       :: fil_cvgso          ! Convergence file
     character(150)       :: fil_solve          ! Solver file
     integer(ip)          :: nzpre              ! Size of preconditioner matrix
     integer(ip)          :: kfl_update_precond ! If preconditionner should be updated
     !
     ! Explicit solver
     !
     real(rp),    contiguous, pointer :: mass(:)            ! Mass matrix for explicit solvers
     !
     ! Mask for dot product
     !
     integer(ip)          :: kfl_mask           ! If mask
     real(rp),    pointer :: mask(:)            ! The mask
     !
     ! Matrix graph
     !
     integer(ip)          :: nnz                ! Size of graph
     integer(ip)          :: nnz_ell            ! Size of ELL graph
     integer(ip), pointer :: ia(:)              ! CSR format: IA
     integer(ip), pointer :: ja(:)              ! CSR format: JA
     integer(ip), pointer :: ia_full(:)         ! CSR format: IA for full row format
     integer(ip), pointer :: ia_full_end(:)     ! CSR format: IA for full row format, end until own nodes
     integer(ip), pointer :: ia_full_ini(:)     ! CSR format: IA for full row format, start after own nodes
     integer(ip), pointer :: ja_full(:)         ! CSR format: JA for full row format
     integer(ip), pointer :: rows(:)            ! COO format: rows
     integer(ip), pointer :: cols(:)            ! COO format: columns
     integer(ip), pointer :: cols_ell(:,:)      ! ELL format: columns
     !
     ! Block Gauss-Seidel
     !
     integer(ip)          :: ndofn_per_block    ! Number of degress of freedom per block
     integer(ip), pointer :: lperm_block(:,:)   ! Permutation arrays BGS
     integer(ip), pointer :: linvp_block(:,:)   ! Inverse permutation arrays BGS
     !
     ! Fixity and boundary conditions
     !
     real(rp),    pointer :: bvess(:,:)         ! Dirichlet value
     integer(ip), pointer :: kfl_fixno(:,:)     ! Fixity array
     integer(ip), pointer :: kfl_fixrs(:)       ! Rotation axis
     real(rp),    pointer :: bvnat(:,:)         ! Neumann value
     !
     ! Direct solver
     !
     type(direct_solver_typ)          :: direct_solver              ! Direct solver
     type(direct_solver_typ), pointer :: direct_solver_RAS(:)       ! Direct solver for RAS
     type(direct_solver_typ)          :: direct_solver_coarse       ! Direct solver for coarse
     type(direct_solver_typ)          :: direct_solver_block_LU     ! Direct solver for block LU
     type(direct_solver_typ)          :: direct_solver_AMG          ! Direct solver for AMG
     type(direct_solver_typ)          :: direct_solver_Deflation    ! Direct solver for deflation
     !
     ! Deflated CG: deflcg.f90
     !
     integer(ip)          :: ngrou              ! Deflated CG: # groups
     integer(ip)          :: nskyl              ! Deflated CG: Skyline A' matrix size
     integer(ip)          :: ifbop              ! Deflated CG: If boundary condition is imposed on boundary nodes
     integer(ip)          :: kfl_gathe          ! Deflated CG: All reduce / or gather
     integer(ip), pointer :: limpo(:)           ! Deflated CG: List of imposed nodes
     integer(ip), pointer :: lgrou(:)           ! Deflated CG: Group list
     integer(ip), pointer :: iskyl(:)           ! Deflated CG: Skyline index
     integer(ip), pointer :: idiag(:)           ! Deflated CG: Pointer to skyline diagonal
     integer(ip), pointer :: iagro(:)           ! Deflated CG: Sparse graph
     integer(ip), pointer :: jagro(:)           ! Deflated CG: Sparse graph
     integer(ip)          :: nzgro              ! Deflated CG: Sparse graph
     integer(ip)          :: icoml              ! Deflated CG: Communication level (for parallelization)
     integer(ip)          :: nbig               ! Deflated CG: All gather strategy
     integer(ip)          :: nsmall             ! Deflated CG: All gather strategy
     integer(ip), pointer :: lcoun(:)           ! Deflated CG: All gather strategy
     integer(ip), pointer :: lbig(:)            ! Deflated CG: All gather strategy
     integer(ip), pointer :: displ(:)           ! Deflated CG: All gather strategy
     integer(4),  pointer :: disp4(:)           ! Deflated CG: All gather strategy
     integer(4),  pointer :: lcou4(:)           ! Deflated CG: All gather strategy
     real(rp),    pointer :: xsmall(:)          ! Deflated CG: All gather strategy
     real(rp),    pointer :: xbig(:)            ! Deflated CG: All gather strategy
     !
     ! Renumbered Gauss-Seidel
     !
     integer(ip), pointer :: permr_gs(:)        ! Renumbered Gauss-Seidel permutation
     integer(ip)          :: ngrou_gs           ! Renumbered Gauss-Seidel group number
     integer(ip), pointer :: invpr_gs(:)        ! Renumbered Gauss-Seidel inverse permutation
     integer(ip), pointer :: lgrou_gs(:)        ! Renumbered Gauss-Seidel groups
     real(rp),    pointer :: vecto_gs(:,:)      ! Renumbered Gauss-Seidel vetor
     integer(ip), pointer :: idiag1(:,:)        ! Stores the position of the diagonal terms in Streamwise Bidiagonal precon.
     real(rp),    pointer :: adiag1(:,:,:)      ! Stores the subdiagonal matrix coefficents in Streamwise Bidiagonal precon.
     real(rp)             :: angle_stream       ! Minimum angle to determine next node renumbered in Streamwise direction
     !
     ! Linelet
     !
     integer(ip)          :: npntr              ! Linelet: # non-zero in trima
     integer(ip)          :: nlpntr             ! Linelet: # points in linelet
     integer(ip)          :: nline              ! Linelet: # linelets
     integer(ip)          :: nlin1              ! Linelet: # linelets crossed by boundary
     integer(ip)          :: npoin              ! Linelet: % of nodes in linelets
     integer(ip)          :: kfl_factl          ! Linelet: Factorization flag
     integer(ip), pointer :: lpntr(:)           ! Linelet: tridiag. to original sym. matrix
     integer(ip), pointer :: lrenu(:)           ! Linelet: point position in tridiag.
     integer(ip), pointer :: lrenup(:)          ! Linelet: inverse permutation
     integer(ip), pointer :: limli(:)           ! Linelet: list of imposed node
     integer(ip), pointer :: lline(:)           ! Linelet: Pointer for each linelet
     real(rp),    pointer :: trima(:)           ! Linelet: Tridiagonal matrix
     real(rp)             :: toler              ! Linelet: Tolerance aspect ratio
     !
     ! RAS and AS
     !
     integer(ip), pointer :: lgrou_ras(:)       ! RAS groups
     type(i1p),   pointer :: permr_ras(:)
     type(i1p),   pointer :: invpr_ras(:)
     type(i1p),   pointer :: ia_ras(:)
     type(i1p),   pointer :: ja_ras(:)
     !
     ! Schur solver
     !
     integer(ip)          :: poaii              ! Pointer to Aii
     integer(ip)          :: poaib              ! Pointer to Aib
     integer(ip)          :: poabi              ! Pointer to Abi
     integer(ip)          :: poabb              ! Pointer to Abb
     integer(ip), pointer :: r_dom_aii(:)       ! 
     integer(ip), pointer :: c_dom_aii(:)       ! 
     integer(ip), pointer :: r_dom_aib(:)       ! 
     integer(ip), pointer :: c_dom_aib(:)       ! 
     integer(ip), pointer :: r_dom_abi(:)       ! 
     integer(ip), pointer :: c_dom_abi(:)       ! 
     integer(ip), pointer :: r_dom_abb(:)       ! 
     integer(ip), pointer :: c_dom_abb(:)       ! 
     !
     ! Block character of matrix, used to compute Reaction residuals
     ! Save all blocks of matrix
     ! e.g.: block(1,2) % a(:,:,:) = matrix of block 1,2 in CSR format
     !       block_number = 3
     !       block_dimensions(1:3) = 1,ndime,1
     !
     ! +-----+--------+-----+
     ! |     |        |     |
     ! | 1,1 |  1,2   | 1,3 |
     ! +-----+--------+-----+
     ! |     |        |     |
     ! | 2,1 |  2,2   | 2,3 |
     ! |     |        |     |
     ! +-----+--------+-----+
     ! |     |        |     |
     ! | 3,1 |  3,2   | 3,3 |
     ! +-----+--------+-----+
     !
     ! lpoin_block(1:npoin) % block2_num(i,j) % matrix(:,:,:)
     ! lpoin_block(1:npoin) % block2_num(i,j) % rhsi(:)
     ! lpoin_block(1:npoin) % block1_num(i)   % reaction(:)
     !
     integer(ip)                     :: kfl_react              ! If reaction should be computed
     type(lpoin_block_typ), pointer  :: lpoin_block(:)
     logical(lg),           pointer  :: lpoin_reaction(:)
     real(rp),              pointer  :: reaction(:,:)
     integer(ip)                     :: num_blocks
     integer(ip)                     :: kfl_block              ! Block solver
     integer(ip)                     :: kfl_block_schur        ! Schur complement based solver for blocks
     integer(ip)                     :: block_num
     integer(ip)                     :: block_dimensions(10)
     integer(ip)                     :: block_sizes(10,10)
     type(block_typ)                 :: block_array(10)
     integer(ip),           pointer  :: block_matrix(:,:)
     integer(ip),           pointer  :: block_rhs(:)
     type(soltyp),          pointer  :: block_solve(:)
     !
     ! Main matrix comes from a Schur complement
     ! A = A1 + A2*A3^-1*A4
     !
     integer(ip)          :: ndofn_A3
     real(rp),    pointer :: A1(:)
     real(rp),    pointer :: A2(:)
     real(rp),    pointer :: A3(:)
     real(rp),    pointer :: A4(:)
     real(rp),    pointer :: invA3(:)
     !
     ! Krylov subspace
     !
     real(rp),    pointer :: krylov_dir(:,:)
     real(rp),    pointer :: krylov_dot_products(:)
     integer(ip)          :: krylov_size
     !
     ! External solvers
     !
     type(DMPH_maphys_t)  :: mphs                 ! MAPHYS or dummy type
#ifdef PETSC
     type(PETSC_LINEAR_SOLVER) :: petsc
#endif

  end type soltyp

  !----------------------------------------------------------------------
  !
  ! Eigenvalue solver
  !
  !----------------------------------------------------------------------

  type eigtyp
     character(50)        :: wprob              ! Names
     character(50)        :: wsolv
     character(50)        :: wprec
     integer(ip)          :: ndofn              ! D.o.F
     integer(ip)          :: ndof2              ! D.o.F*D.o.F
     integer(ip)          :: neiva              ! number of eigen values requiered
     integer(ip)          :: neige              ! Size of eigenvector
     integer(ip)          :: nzmbt              ! RHS eigen matrix
     integer(ip)          :: nzmat              ! LHS matrix
     integer(ip)          :: nzrhs              ! RHS
     integer(ip)          :: kfl_algso          ! Solver parameter
     integer(ip)          :: kfl_massm          ! Mass matrix type
     integer(ip)          :: miter              ! Solver number of iterations
     integer(ip)          :: kfl_facto          ! Factorization
     integer(ip)          :: itsol(3)           ! Output: Solver statistics
     integer(ip)          :: nsolv              ! Output: # solves
     integer(ip)          :: lun_cvgei
     integer(ip)          :: lun_solei
     integer(ip)          :: kfl_preco          ! Preconditioner
     real(rp)             :: solco              ! Tolerance
     real(rp)             :: shift              ! shift
     integer(ip)          :: diter              ! iter en eigdir
     integer(ip)          :: error              ! error en eigdir
     integer(ip), pointer :: lpdof(:)
  end type eigtyp

contains

  subroutine init_soldat(self)

    class(soldat) :: self

  end subroutine init_soldat
  
  subroutine init_solout(self)

    class(solout) :: self
    
  end subroutine init_solout
  
end module def_kintyp_solvers
!> @}
