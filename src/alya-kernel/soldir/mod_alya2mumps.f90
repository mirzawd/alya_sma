!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!>------------------------------------------------------------------------
!> @addtogroup Direct_Solver
!> @{
!> @file    mod_alya2mumps.f90
!> @author  Adria Quintanas-Corominas
!> @author  Guillaume Houzeaux
!> @date    2021
!> @brief   ToolBox for MUMPS
!> @details ToolBox for MUMPS
!>
!------------------------------------------------------------------------

module mod_alya2mumps

  use def_kintyp_basic,         only : ip, rp, lg
  use def_kintyp_solvers,       only : direct_solver_typ
  use def_master,               only : INOTMASTER, IMASTER, IPARALL
  use mod_communications_tools, only : PAR_COMM_TO_INT
#ifdef MUMPS
  use omp_lib
#endif

  implicit none 
  private

#ifdef MUMPS
  INTERFACE
     SUBROUTINE DMUMPS( id )
       !!!DEC$ ATTRIBUTES C,REFERENCE,NOMIXED_STR_LEN_ARG, ALIAS:'_DMUMPS'   ::
       !!!DMUMPS
       INCLUDE 'dmumps_struc.h'
       TYPE (DMUMPS_STRUC) :: id
     END SUBROUTINE DMUMPS
  END INTERFACE
#endif

  public  :: direct_solver_initialization_mumps  ! Symbolic factorization
  public  :: direct_solver_factorization_mumps   ! Numerical factorization
  public  :: direct_solver_solution_mumps        ! Solution of system Ax=b
  public  :: direct_solver_cleaning_mumps        ! Clean MUMPS structures
  private :: mumps_csr_to_coo

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  Adria Quintanas-Corominas
  !> @author  Guillaume Houzeaux
  !> @date    2021
  !> @brief   Initializaiton of MUMPS
  !> @details Initializaiton of MUMPS
  !> 
  !-----------------------------------------------------------------------
  subroutine direct_solver_initialization_mumps( direct_solver )
#ifdef MUMPS
    use def_master,         only : kfl_paral
    use def_kermod,         only : ndivi
    use def_domain,         only : npoin, npoin_origi, meshe
    use mod_parall,         only : PAR_COMM_MY_CODE,PAR_CODE_SIZE
#endif
    use mod_renumbering,    only : renumbering_lexical_order_type
    use mod_communications, only : PAR_GATHER
    use mod_communications, only : PAR_SCATTERV
    use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE 

    type(direct_solver_typ), intent(inout) :: direct_solver 
#ifdef MUMPS
    integer(ip)                            :: ii
    integer(ip), pointer                   :: lninv_lex(:) => null()

    ! Set  basic options
    direct_solver % mumps_par % SYM  = 0_ip          ! Unsymmetric problem
    direct_solver % mumps_par % COMM = PAR_COMM_TO_INT(PAR_COMM_MY_CODE)
    if( IPARALL )then
       direct_solver % mumps_par % PAR  = 0_ip       ! The host is not involved on factorisation and solve.
    else
       direct_solver % mumps_par % PAR  = 1_ip       ! The host is involved on factorisation and solve.
    endif

    ! Initialise MUMPS structures
    direct_solver % mumps_par % JOB = -1_4
    call DMUMPS(direct_solver % mumps_par)

    if( direct_solver % mumps_par % INFO(1) < 0_ip )then
       call runend('ERROR IN MOD_ALYA2MUMPS: INITIALIZATION PHASE FAILED')   
    endif

    ! Set dimensions and allocate structures
    if( IPARALL ) then
       if( IMASTER )then
          ! Defining order ot he system
          direct_solver % mumps_par % N = npoin_origi * direct_solver % ndof
          allocate( direct_solver % mumps_par % RHS( direct_solver % mumps_par % N ) )
          ! Number of RHS 
          direct_solver % mumps_par % NRHS = 1_ip
       else
          ! Distributed MATRIX
          direct_solver % mumps_par % NZ_LOC = direct_solver % nz * direct_solver % ndof**2
          allocate( direct_solver % mumps_par % A_LOC( direct_solver % mumps_par % NZ_LOC ) )
          allocate( direct_solver % mumps_par % IRN_LOC( direct_solver % mumps_par % NZ_LOC ) )
          allocate( direct_solver % mumps_par % JCN_LOC( direct_solver % mumps_par % NZ_LOC ) )
          ! Distributed RHS
          direct_solver % mumps_par % NLOC_RHS = direct_solver % nn * direct_solver % ndof
          direct_solver % mumps_par % LRHS_LOC = direct_solver % nn * direct_solver % ndof
          allocate( direct_solver % mumps_par % RHS_LOC( direct_solver % mumps_par % NLOC_RHS ) )
          allocate( direct_solver % mumps_par % IRHS_LOC( direct_solver % mumps_par % NLOC_RHS ) )
       endif
  
    else
       ! RHS 
       direct_solver % mumps_par % NRHS = 1_ip ! How many RHS?
       direct_solver % mumps_par % N = npoin * direct_solver % ndof 
       allocate( direct_solver % mumps_par % RHS( direct_solver % mumps_par % N ) )
       ! MATRIX
       direct_solver % mumps_par % NZ = direct_solver % nz * direct_solver % ndof**2 
       allocate( direct_solver % mumps_par % A( direct_solver % mumps_par % NZ ) )
       allocate( direct_solver % mumps_par % IRN( direct_solver % mumps_par % NZ ) )
       allocate( direct_solver % mumps_par % JCN( direct_solver % mumps_par % NZ ) )
    endif

    ! Set options  
    ! - Output stream for global information, collected on the host
    ! TODO: To output the results in this file, it is needed a clean of the alya file
    direct_solver % mumps_par % ICNTL(3) = 230421
    ! - More output
    if( kfl_paral < 2_ip )then
       direct_solver % mumps_par % ICNTL(4) = 2_ip
    else
       direct_solver % mumps_par % ICNTL(4) = 0_ip
    endif
    !direct_solver % mumps_par % ICNTL(1) = 200_ip + (kfl_paral - 1_ip) * 3_ip ! Output stream for error messages
    !direct_solver % mumps_par % ICNTL(2) = 201_ip + (kfl_paral - 1_ip) * 3_ip ! Output stream for diagnostic printing, statistics, and warning messages
    !direct_solver % mumps_par % ICNTL(3) = 202_ip + (kfl_paral - 1_ip) * 3_ip ! Output stream for global information, collected on the host
    ! - CCO MATRIX format 
    direct_solver % mumps_par % ICNTL(5) = 0_ip
    ! - Distributed MATRIX (the usr provides everything)
    if( IPARALL ) direct_solver % mumps_par % ICNTL(18) = 3_ip
    ! - Distributed RHS (the usr provides the order) 
    if( IPARALL ) direct_solver % mumps_par % ICNTL(20) = 10_ip 
    ! - Iterative Refinament
    !direct_solver % mumps_par % ICNTL(10) = 5000_ip ! Number of iterations
    !direct_solver % mumps_par % CNTL(2) = -1.0_rp   ! Tolerance of the refinment
    ! - Error Analysis (only for centralized RHS
    !direct_solver % mumps_par % ICNTL(11) = 2_ip 
    ! - Performance flags suggested by Alfredo and VIncent
    direct_solver % mumps_par % KEEP(370) = 1_ip
    direct_solver % mumps_par % KEEP(371) = 1_ip
    ! - Multithreading
    direct_solver  % num_threads = max(1_ip, int(OMP_GET_MAX_THREADS(),ip)) 
    if( direct_solver  % num_threads > 1_ip )then
       direct_solver % mumps_par % KEEP(401) = 1_ip 
    else
       direct_solver % mumps_par % KEEP(401) = 0_ip
    endif
    ! - BLR options suggested by Alfredo and Vincent
    !direct_solver % mumps_par % ICNTL(15) = -3_ip ! - Compress the analysis graph considering 3x3 blocks
    !direct_solver % mumps_par % ICNTL(35) = 1_ip
    !direct_solver % mumps_par % ICNTL(36) = 1_ip  ! - BLR without pivoting
    !direct_solver % mumps_par % CNTL(7) = 1.0e-10_rp
    ! - Pivoting stuff
    direct_solver % mumps_par % CNTL(1) = 0.0_rp  ! - No pivoting

    ! Uncomment to write MATRIX and RHS. Then move initialisation after 
    ! filling the vectors A, RHS
    !direct_solver % mumps_par % WRITE_PROBLEM = 'alya'

    ! Create lexicograph array (MPI_1 : 1,2,3,4; MPI_2 : 5,6,7,8 ...)
    if( IPARALL )then
       call renumbering_lexical_order_type( &
          meshe(ndivi), &
          lninv_lex, &
          'WITHOUT HALOS', &
          direct_solver % ndof )
    else
       if( .not. associated(lninv_lex) ) allocate(lninv_lex(direct_solver % nn * direct_solver % ndof))
       do ii = 1, direct_solver % nn * direct_solver % ndof
          lninv_lex(ii) = ii
       end do
    endif
    ! From Alya B-CSR formato to MUMPS arrays
    if( IPARALL )then
       if( INOTMASTER )then
          call mumps_csr_to_coo( &
             direct_solver % nn, direct_solver % ndof, &
             lninv_lex, &
             direct_solver % ia, direct_solver % ja, &
             direct_solver % mumps_par % IRN_LOC, & 
             direct_solver % mumps_par % JCN_LOC )
          ! Distributed RHS
          direct_solver % mumps_par % IRHS_LOC(:) = lninv_lex(:)
       endif
    else 
       call mumps_csr_to_coo( &
          direct_solver % nn, direct_solver % ndof, &
          lninv_lex, &
          direct_solver % ia, direct_solver % ja, &
          direct_solver % mumps_par % IRN, &
          direct_solver % mumps_par % JCN )
    endif

    ! Job analysis (it includes symbolic factorisation)
    direct_solver % mumps_par % JOB = 1_4
    call DMUMPS( direct_solver % mumps_par )

    if( direct_solver % mumps_par % INFO(1) < 0_ip )then
       call runend('ERROR IN MOD_ALYA2MUMPS: ANALYSIS PHASE FAILED')
    endif

    ! Number of DoF that each MPI has 
    if( IPARALL )then
       if( IMASTER )then
          allocate( direct_solver % gatsca_mumps(PAR_CODE_SIZE) )
       else
          allocate( direct_solver % gatsca_mumps(1) )
       endif
       call PAR_GATHER( direct_solver % nn_own * direct_solver % ndof, direct_solver % gatsca_mumps,'IN MY CODE')
    endif

    ! Release memmory
    if( associated(lninv_lex) )then
       deallocate(lninv_lex); lninv_lex => null()
    endif

#endif

  end subroutine direct_solver_initialization_mumps 

  !-----------------------------------------------------------------------
  !> 
  !> @author  Adria Quintanas-Corominas
  !> @author  Guillaume Houzeaux
  !> @date    2021
  !> @brief   MUMPS factorization
  !> @details Analysis and factorization with MUMPS
  !> 
  !-----------------------------------------------------------------------

  subroutine direct_solver_factorization_mumps( direct_solver, amatr )

    type(direct_solver_typ), intent(inout) :: direct_solver    
    real(rp), target,        intent(in)    :: amatr( direct_solver % nz * direct_solver % ndof**2 )
#ifdef MUMPS
    real(rp)                               :: nz_own


    ! From Alya to MUMPS structure
    nz_own = direct_solver % nz * direct_solver % ndof**2
    if( IPARALL )then
       if( INOTMASTER )then 
          direct_solver % mumps_par % A_LOC(1:nz_own) = amatr(1:nz_own)
       endif
    else
       direct_solver % mumps_par % A(1:nz_own) = amatr(1:nz_own)
    endif

    ! Job factorisation
    direct_solver % mumps_par % JOB = 2_4
    call DMUMPS( direct_solver % mumps_par )

    if( direct_solver % mumps_par % INFO(1) < 0_ip )then
       call runend('ERROR IN MOD_ALYA2MUMPS: FACTORIZATION PHASE FAILED')
    endif

#endif

  end subroutine direct_solver_factorization_mumps

  !-----------------------------------------------------------------------
  !> 
  !> @author  Adria Quintanas-Corominas
  !> @author  Guillaume Houzeaux
  !> @date    2021
  !> @brief   MUMPS factorization
  !> @details Analysis and factorization with MUMPS
  !> 
  !-----------------------------------------------------------------------

  subroutine direct_solver_solution_mumps( direct_solver, rhsid, unkno )
    use mod_communications, only : PAR_SCATTERV 
    use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
    use mod_communications, only : PAR_GATHERV

    type(direct_solver_typ), intent(inout) :: direct_solver    
    real(rp), target,        intent(in)    :: rhsid( direct_solver % nn * direct_solver % ndof)
    real(rp), target,        intent(out)   :: unkno( direct_solver % nn * direct_solver % ndof)
#ifdef MUMPS
    integer                                :: nn_own


    ! From Alya to MUMPS structure
    nn_own = direct_solver % nn_own * direct_solver % ndof
    if( IPARALL )then
       if( IMASTER )then
          direct_solver % mumps_par % RHS = 0.0_rp
       endif
       ! - Distributed RHS
       if( INOTMASTER )then
          direct_solver % mumps_par % RHS_LOC(:) = 0.0_rp
          direct_solver % mumps_par % RHS_LOC(1:nn_own) = rhsid(1:nn_own) 
       endif
       ! - Centralized RHS
       !call PAR_GATHERV( rhsid, direct_solver % mumps_par % RHS, nn_own, direct_solver % gatsca_mumps, 'IN MY CODE')
    else
       direct_solver % mumps_par % RHS(1:nn_own) = rhsid(1:nn_own) 
    endif

    ! Job solver
    direct_solver % mumps_par % JOB = 3_4
    call DMUMPS(direct_solver % mumps_par)

    if( direct_solver % mumps_par % INFO(1) < 0_ip )then
       call runend('ERROR IN MOD_ALYA2MUMPS: SOLVER PHASE FAILED')
    endif

    ! Scatter the solution
    unkno(:) = 0.0_rp
    if( IPARALL ) then
       call PAR_SCATTERV(direct_solver % mumps_par % RHS, unkno, direct_solver % gatsca_mumps, nn_own, 'IN MY CODE')
       call PAR_INTERFACE_NODE_EXCHANGE(direct_solver % ndof,unkno, 'SUM', 'IN MY CODE')
    else 
        unkno(1:nn_own) = direct_solver % mumps_par % RHS(1:nn_own)
    endif

#endif    

  end subroutine direct_solver_solution_mumps

  !-----------------------------------------------------------------------
  !> 
  !> @author  Adria Quintanas-Corominas
  !> @author  Guillaume Houzeaux
  !> @date    2021
  !> @brief   MUMPS factorization
  !> @details Analysis and factorization with MUMPS
  !> 
  !-----------------------------------------------------------------------

  subroutine direct_solver_cleaning_mumps( direct_solver )
    
    type(direct_solver_typ), intent(inout) :: direct_solver    

#ifdef MUMPS

    if( IPARALL )then
       if( IMASTER )then
          deallocate( direct_solver % mumps_par % RHS )
       else
          ! - Distributed MATRIX
          deallocate( direct_solver % mumps_par % A_LOC   ) 
          deallocate( direct_solver % mumps_par % IRN_LOC ) 
          deallocate( direct_solver % mumps_par % JCN_LOC )
          ! - Distributed RHS
          deallocate( direct_solver % mumps_par % RHS_LOC )
          deallocate( direct_solver % mumps_par % IRHS_LOC )
       endif
       deallocate( direct_solver % gatsca_mumps )
    else
       deallocate( direct_solver % mumps_par % RHS )
       deallocate( direct_solver % mumps_par % A )
       deallocate( direct_solver % mumps_par % IRN )
       deallocate( direct_solver % mumps_par % JCN )
    endif

    direct_solver % mumps_par % JOB = -2_4
    call DMUMPS(direct_solver % mumps_par)

#endif    
    
  end subroutine direct_solver_cleaning_mumps


  !-----------------------------------------------------------------------
  !> 
  !> @author  Adria Quintanas-Corominas
  !> @author  Guillaume Houzeaux
  !> @date    2021
  !> @brief   MUMPS BCSR to global COO format
  !> @details Analysis and factorization with MUMPS
  !-----------------------------------------------------------------------

  subroutine mumps_csr_to_coo( nn, ndofn, ninv, ia, ja, row, col )
    
    integer(ip),  intent(in)              :: nn
    integer(ip),  intent(in)              :: ndofn
    integer(ip),  intent(in)              :: ninv(:) ! From local to global COO
    integer(ip),  intent(in)              :: ia(nn+1)
    integer(ip),  intent(in)              :: ja(*)
    integer(ip),  intent(inout), pointer  :: row(:)
    integer(ip),  intent(inout), pointer  :: col(:)
    integer(ip)                           :: ii, jj, iz, kz, nz, nnz
    integer(ip)                           :: idofn, jdofn, icoo, jcoo

    ! Size of graph
    nnz = ia(nn+1)-1
    nz = nnz * ndofn * ndofn

    ! Compute COO format including symmetric graph
    if( ndofn == 1_ip ) then

       ! One degree of freedom per node
       do ii = 1, nn
          do iz = ia(ii), ia(ii+1) - 1
             jj      = ja(iz)
             row(iz) = ninv(ii)
             col(iz) = ninv(jj)
          end do
       end do

    else if( ndofn > 1 ) then

       ! Explode graph
       kz = 0_ip
       do ii = 1, nn
          do iz = ia(ii), ia(ii+1) - 1
             jj = ja(iz)
             do idofn = 1, ndofn
                icoo = (ii - 1_ip) * ndofn + idofn
                do jdofn = 1, ndofn 
                   jcoo = (jj - 1_ip) * ndofn + jdofn
                   kz       = kz + 1_ip
                   row(kz)  = ninv(icoo)
                   col(kz)  = ninv(jcoo)
                end do
             end do
          end do
       end do
       if( kz /= nz ) then
          print*,'MATRIX SIZE=',kz,nz
          call runend('GRAPHS_CSR_TO_COO: PROBLEM WITH MATRIX SIZE')
       end if

    end if

  end subroutine mumps_csr_to_coo

end module mod_alya2mumps
