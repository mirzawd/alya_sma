!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup WRAPPERS
!> @{
!> @file    def_alya2petsc.f90
!> @author  Adria Quintanas-Corominas
!> @brief   Bridge to PETSc
!-----------------------------------------------------------------------

module def_alya2petsc
    

#ifdef PETSC
#define ZPLS PETSC_LINEAR_SOLVER
#include <petsc/finclude/petscksp.h>

    use def_kintyp_basic, only: lg, ip, rp
    use def_kintyp_petsc, only: PETSC_LINEAR_SOLVER
    use petscksp

    PetscInt,    parameter :: ZERO = 0
    PetscInt,    parameter :: NSOLAPE = 1
    PetscScalar, parameter :: MINUSONER = -1.0_rp
    PetscReal,   parameter :: RELTOL = 1.0e-10_rp

    private

    public :: PETSC_LINEAR_SOLVER
    public :: ZERO
    public :: NSOLAPE
    public :: MINUSONER
    public :: RELTOL
    public :: petscKSP_initialise
    public :: petscKSP_copyAlyaSolver
    public :: petscKSP_finalise
    public :: petscKSP_allocateSystem
    public :: petscKSP_deallocateSystem
    public :: petscKSP_createSolver
    public :: petscKSP_destroySolver

    contains

    ! =================================================================
    ! 
    ! Methods of the PETSC_LINEAR_SOLVER type
    !
    ! > PETSC_LINEAR_SOLVER not contains the methods in the class to 
    ! > avoid circular dependency with soltyp of def_kintyp_solvers
    !
    subroutine petscKSP_nullify(alyapetsc)

        type(ZPLS),   intent(inout) :: alyapetsc

        alyapetsc % I_PARALL = .false.
        alyapetsc % I_WORKER = .false.
        alyapetsc % I_WORKER_NOT_EMPTY = .false.
        alyapetsc % ndofn = 0
        alyapetsc % nequa = 0
        alyapetsc % nequa_own = 0
        alyapetsc % nunkn = 0
        alyapetsc % nunkn_own = 0
        alyapetsc % nzmat = 0
        alyapetsc % nzmat_own = 0
        
        nullify(alyapetsc % ia)
        nullify(alyapetsc % ja)
        nullify(alyapetsc % ia_full)
        nullify(alyapetsc % ja_full)
        nullify(alyapetsc % ia_full_cooLex)
        nullify(alyapetsc % ja_full_cooLex)
        nullify(alyapetsc % permr)
        nullify(alyapetsc % nndiag)
        nullify(alyapetsc % nnoutd)


    end subroutine petscKSP_nullify


    subroutine petscKSP_initialise(alyapetsc)

        type(ZPLS), intent(inout) :: alyapetsc

        call petscKSP_nullify(alyapetsc)

    end subroutine petscKSP_initialise
    

    subroutine petscKSP_finalise(alyapetsc)

        type(ZPLS), intent(inout) :: alyapetsc

        ! Reset dimensions
        alyapetsc % ndofn = 0_ip
        alyapetsc % nequa = 0_ip
        alyapetsc % nequa_own = 0_ip
        alyapetsc % nunkn = 0_ip
        alyapetsc % nunkn_own = 0_ip
        alyapetsc % nzmat = 0_ip
        alyapetsc % nzmat_own = 0_ip

        ! Deallocate & Destroy memmory
        call petscKSP_deallocateSystem(alyapetsc)
        call petscKSP_destroySolver(alyapetsc)
        call petscKSP_nullify(alyapetsc)
 
    end subroutine petscKSP_finalise
    

    subroutine petscKSP_copyAlyaSolver(alyapetsc, ALYA_IPARALL, ALYA_IMASTER, ALYA_IEMPTY, solve)
        use def_kintyp_solvers, only: soltyp

        type(ZPLS),   intent(inout) :: alyapetsc
        logical(lg),  intent(in)    :: ALYA_IPARALL
        logical(lg),  intent(in)    :: ALYA_IMASTER
        logical(lg),  intent(in)    :: ALYA_IEMPTY
        type(soltyp), intent(inout) :: solve

        alyapetsc % I_PARALL = ALYA_IPARALL
        alyapetsc % I_WORKER = .not. (ALYA_IPARALL .and. ALYA_IMASTER)
        alyapetsc % I_WORKER_NOT_EMPTY = alyapetsc % I_WORKER .and. (.not. ALYA_IEMPTY)
        alyapetsc % ndofn = solve % ndofn
        alyapetsc % nequa = solve % nequa
        alyapetsc % nequa_own = solve % nequa_own
        alyapetsc % nunkn = solve % nunkn
        alyapetsc % nunkn_own = solve % nequa_own * solve % ndofn
        alyapetsc % nzmat = solve % nzmat
        alyapetsc % nzmat_own = solve % nzmat_own

        if( alyapetsc % I_WORKER_NOT_EMPTY )then

            ! CSR or CSC: Matrix graph
            ! - rows
            if( .not. associated(alyapetsc % ia) )then
                allocate(alyapetsc % ia(size(solve % ia)))
                alyapetsc % ia = solve % ia
            endif
            ! - columns
            if( .not. associated(alyapetsc % ja) )then
                allocate(alyapetsc % ja(size(solve % ja)))
                alyapetsc % ja = solve % ja
            endif
            
            ! CSR or CSC: Full matrix graph
            ! - rows
            if( .not. associated(alyapetsc % ia_full) )then
                allocate(alyapetsc % ia_full(size(solve % ia_full)))
                alyapetsc % ia_full = solve % ia_full
            endif
            ! - columns
            if( .not. associated(alyapetsc % ja_full) )then
                allocate(alyapetsc % ja_full(size(solve % ja_full)))
                alyapetsc % ja_full = solve % ja_full
            endif

        endif

    end subroutine petscKSP_copyAlyaSolver


    subroutine petscKSP_allocateSystem(alyapetsc)

        type(ZPLS),      intent(inout) :: alyapetsc

        select case(alyapetsc % PETSC_STRATEGY)
        case('mpiaij'); call allocate_MPIAIJ(alyapetsc)
        end select

    end subroutine petscKSP_allocateSystem


    subroutine petscKSP_deallocateSystem(alyapetsc)

        type(ZPLS), intent(inout) :: alyapetsc
        
        select case(alyapetsc % PETSC_STRATEGY)
        case('mpiaij'); call deallocate_MPIAIJ(alyapetsc)
        end select

    end subroutine petscKSP_deallocateSystem


    subroutine petscKSP_createSolver(alyapetsc)

        type(ZPLS),    intent(inout) :: alyapetsc
        PetscErrorCode               :: ierr

        if( alyapetsc % I_WORKER )then
            call KSPCreate(PETSC_COMM_WORLD, alyapetsc % ksp, ierr)
            call KSPSetOperators(alyapetsc % ksp, alyapetsc % A, alyapetsc % A, ierr)
            call KSPSetTolerances(alyapetsc % ksp, &
            RELTOL, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)
            call KSPSetType(alyapetsc % ksp, KSPCG, ierr)
            call KSPGetPC(alyapetsc % ksp, alyapetsc % pc, ierr)
            call PCSetType(alyapetsc % pc, PCJACOBI, ierr)
            call PCASMSetOverlap(alyapetsc % pc, NSOLAPE, ierr)
            call KSPSetFromOptions(alyapetsc % ksp, ierr)
        endif

    end subroutine petscKSP_createSolver


    subroutine petscKSP_destroySolver(alyapetsc)

        type(ZPLS),     intent(inout) :: alyapetsc
        PetscErrorCode                :: ierr

        if( alyapetsc % I_WORKER )then
            ! PC is destroyed along with KSP
            !call PCDestroy(alyapetsc % pc, ierr)
            call KSPDestroy(alyapetsc % ksp, ierr)
        endif

    end subroutine petscKSP_destroySolver


    ! =================================================================
    ! 
    ! MPIAIJ type procedures
    !
    subroutine allocate_MPIAIJ(alyapetsc)

        type(ZPLS),      intent(inout) :: alyapetsc
        PetscErrorCode                 :: ierr

        ! Preliminar sutff (mainly Alya stuff)
        call compute_LexicographicalCOO_from_BCSR(alyapetsc)
        call compute_MPIAIJ_MATRIX_strategy(alyapetsc)

        ! PETSC structures
        if( alyapetsc % I_WORKER )then
            ! Allocate matrix
            call MatCreate(PETSC_COMM_WORLD, alyapetsc % A, ierr)
            call MatSetType(alyapetsc % A, MATMPIAIJ, ierr)
            call MatSetSizes(alyapetsc % A, alyapetsc % nunkn_own, alyapetsc % nunkn_own, PETSC_DETERMINE, PETSC_DETERMINE, ierr)
            call MatMPIAIJSetPreallocation(alyapetsc % A, ZERO, alyapetsc % nndiag, ZERO, alyapetsc % nnoutd, ierr)

            ! Allocate RHS
            call VecCreate(PETSC_COMM_WORLD, alyapetsc % b, ierr)
            call VecSetType(alyapetsc % b, VECMPI, ierr)
            call VecSetSizes(alyapetsc % b, alyapetsc % nunkn_own, PETSC_DETERMINE, ierr)
    
            ! Allocate UNKNO
            call VecDuplicate(alyapetsc % b, alyapetsc % x, ierr)
        endif

    end subroutine allocate_MPIAIJ


    subroutine deallocate_MPIAIJ(alyapetsc)

        type(ZPLS),     intent(inout) :: alyapetsc
        PetscErrorCode                :: ierr
    
        ! Alya structures
        if( alyapetsc % I_WORKER_NOT_EMPTY )then
            deallocate(alyapetsc % ia)
            deallocate(alyapetsc % ja)
            deallocate(alyapetsc % ia_full)
            deallocate(alyapetsc % ja_full)
            deallocate(alyapetsc % ia_full_cooLex)
            deallocate(alyapetsc % ja_full_cooLex)
            deallocate(alyapetsc % permr)
            deallocate(alyapetsc % nndiag)
            deallocate(alyapetsc % nnoutd)
        endif

        ! PETSC structures
        if( alyapetsc % I_WORKER )then
            call MatDestroy(alyapetsc % A, ierr)
            call VecDestroy(alyapetsc % b, ierr)
            call VecDestroy(alyapetsc % x, ierr)
        endif

    end subroutine deallocate_MPIAIJ


    subroutine compute_MPIAIJ_MATRIX_strategy(alyapetsc)
        ! Get number of ZEROs per row in the diagonal and off-diagonal 
        ! for PETSc preallocation purposes
        
        type(ZPLS),   intent(inout) :: alyapetsc
        integer(ip)                 :: ii, iz, jj, id, jd, aa, bb
        integer(ip)                 :: insd, outs, out2
    
        if( alyapetsc % I_WORKER_NOT_EMPTY )then
            if( .not. associated(alyapetsc % nndiag) )then
                allocate(alyapetsc % nndiag(alyapetsc % nunkn_own))
                alyapetsc % nndiag = 0
            endif

            if( .not. associated(alyapetsc % nnoutd) )then
                allocate(alyapetsc % nnoutd(alyapetsc % nunkn_own))
                alyapetsc % nnoutd = 0
            endif
        endif

        ! Loop BCSR
        do ii = 1, alyapetsc % nequa_own
            do iz = alyapetsc % ia_full(ii), alyapetsc % ia_full(ii+1)-1
                jj = alyapetsc % ja_full(iz)
                ! Recover COO rows
                do aa = 1, alyapetsc % ndofn
                    id = (ii - 1_ip) * alyapetsc % ndofn + aa
                    ! Initialise inside / outside counters
                    insd = 0_ip
                    outs = 0_ip 
                    out2 = 0_ip 
                    ! Recover COO columns
                    do bb = 1, alyapetsc % ndofn 
                        jd = (jj - 1_ip) * alyapetsc % ndofn + bb
                        ! Count according if inside (insd) or outside (outs) the diagonal
                        if( id > alyapetsc % nunkn_own .or. jd > alyapetsc % nunkn_own )then
                            outs = outs + 1_ip
                        else
                            insd = insd + 1_ip
                        endif
                    end do
                    alyapetsc % nndiag(id) = alyapetsc % nndiag(id) + insd
                    alyapetsc % nnoutd(id) = alyapetsc % nnoutd(id) + outs
                end do
            end do
        end do
    
    end subroutine compute_MPIAIJ_MATRIX_strategy


    ! =================================================================
    ! 
    ! General procedures
    !
    subroutine compute_LexicographicalCOO_from_BCSR(alyapetsc)
        use def_master, only: kfl_paral

        type(ZPLS),      intent(inout) :: alyapetsc
        logical(lg)                    :: INOTMASTER
        integer(ip)                    :: ii, iz, jj, kz
        integer(ip)                    :: idofn, ii_idofn, jdofn, jj_jdofn
        integer(ip)                    :: nunkn, nequa, nequa_own, ndofn
        integer(ip), pointer           :: permr(:)
    
        ! Allocate memmory
        if( alyapetsc % I_WORKER_NOT_EMPTY )then
            if( .not. associated(alyapetsc % ia_full_cooLex) )then
                allocate(alyapetsc % ia_full_cooLex(alyapetsc % nzmat))
                alyapetsc % ia_full_cooLex = 0
            endif

            if( .not. associated(alyapetsc % ja_full_cooLex) )then
                allocate(alyapetsc % ja_full_cooLex(alyapetsc % nzmat))
                alyapetsc % ja_full_cooLex = 0
            endif

            if( .not. associated(alyapetsc % permr) ) then
                allocate(alyapetsc % permr(alyapetsc % nunkn))
                alyapetsc % permr = 0
            end if
        endif

        ! Get lexicographical coo format
        ! - from PetscInt to AlyaInt
        nunkn = int(alyapetsc % nunkn,ip)
        nequa = int(alyapetsc % nequa,ip)
        nequa_own = int(alyapetsc % nequa_own,ip)
        ndofn = int(alyapetsc % ndofn,ip)
        allocate(permr(nunkn))
        ! > equivalent to call renumbering_lexical_order_type( &
        ! >     meshe, alyapetsc % permr, 'WITHOUT HALOS', alyapetsc % ndofn)
        ! > which is in alya-kernel/domain/mod_renumbering.
        call renumbering_lexical_order_type(&
            nequa, nequa_own, ndofn, permr)
        ! - from AlyaInt to PetscInt
        if( alyapetsc % I_WORKER_NOT_EMPTY )then
            do ii = 1, alyapetsc % nunkn
                alyapetsc % permr(ii) = permr(ii)
            enddo
        endif
        deallocate(permr)

        ! From BSCR to COO (considering the lexicographical order)
        if( alyapetsc % I_WORKER_NOT_EMPTY )then
            kz = 0_ip
            do ii = 1, alyapetsc % nequa
                do iz = alyapetsc % ia(ii), alyapetsc % ia(ii+1)-1
                    jj = alyapetsc % ja(iz)
                    do idofn = 1, alyapetsc % ndofn
                        ii_idofn = (ii - 1_ip) * alyapetsc % ndofn + idofn
                        do jdofn = 1, alyapetsc % ndofn 
                            jj_jdofn = (jj - 1_ip) * alyapetsc % ndofn + jdofn
                            kz       = kz + 1_ip
                            alyapetsc % ia_full_cooLex(kz)  = alyapetsc % permr(ii_idofn) - 1_ip
                            alyapetsc % ja_full_cooLex(kz)  = alyapetsc % permr(jj_jdofn) - 1_ip
                        end do
                    end do
                end do
            end do
        endif

        ! Uncomment for debbugging proposes 
        !if( kz /= ((meshe % r_dom(nzpoi+1) - 1_ip) * (alyapetsc % ndofn **2)) ) then
        !    print*, kz, ((meshe % r_dom(nzpoi+1) - 1_ip) * (alyapetsc % ndofn **2))
        !    call runend('GRAPHS_CSR_TO_COO: PROBLEM WITH MATRIX SIZE')
        !end if

    end subroutine compute_LexicographicalCOO_from_BCSR


    ! ===================================================================
    ! 
    ! Procedures copyed from Kernel 
    ! 
    subroutine renumbering_lexical_order_type(npoin, npoin_own, ndofn, permr)
        use def_domain, only: memor_dom
        use def_master, only: INOTMASTER, kfl_paral, npart
        use mod_memory, only: memory_alloca, memory_deallo
        use mod_communications, only: PAR_ALLGATHER, PAR_INTERFACE_NODE_EXCHANGE

        integer(ip),              intent(in)    :: npoin             !< Nodes
        integer(ip),              intent(in)    :: npoin_own         !< Nodes
        integer(ip),              intent(in)    :: ndofn             !< If dof > 1
        integer(ip),     pointer, intent(inout) :: permr(:)          !< Permutation
        integer(ip),     pointer                :: npoin_own_tot(:)
        integer(ip)                             :: ipoin
        integer(ip)                             :: npoin_offset
        !
        ! Permute
        !
        nullify(npoin_own_tot)
        !
        ! Renumber nodes (without halos)
        !
        !if( .not. associated(permr) .and. INOTMASTER ) then
        !    call memory_alloca(memor_dom,'permr','renumbering_lexical_order',permr,npoin*ndofn)
        !end if
        call memory_alloca(memor_dom, 'npoin_own_tot', 'renumbering_lexical_order', npoin_own_tot, npart+1,' INITIALIZE', 0_ip)
        call PAR_ALLGATHER(npoin_own, npoin_own_tot)
        if( INOTMASTER )then
            npoin_offset = sum(npoin_own_tot(0:kfl_paral-1))*ndofn
            do ipoin = 1, npoin_own * ndofn
                permr(ipoin) = ipoin + npoin_offset
            end do
            call PAR_INTERFACE_NODE_EXCHANGE(ndofn, permr, 'MAX', 'IN MY CODE')
        end if
        call memory_deallo(memor_dom,'npoin_own_tot','renumbering_lexical_order',npoin_own_tot)

    end subroutine renumbering_lexical_order_type


#endif

end module def_alya2petsc 

