!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup WRAPPERS
!> @{
!> @file    mod_alya2petsc.f90
!> @author  Adria Quintanas-Corominas
!> @brief   Bridge to PETSc
!-----------------------------------------------------------------------

module mod_alya2petsc

#ifdef PETSC
#define ZPLS PETSC_LINEAR_SOLVER 
#include <petsc/finclude/petscksp.h>

    use def_kintyp_basic, only : lg, ip, rp
    use def_alya2petsc
    use petscksp
    use iso_c_binding

    implicit none

    logical(lg)    :: kfl_petsc = .false.
    integer(ip)    :: petsc_file_unit = 23042221
    character(150) :: petsc_file_name
    PetscErrorCode :: ierr

    public

    contains

    subroutine alya2petsc_read_configuration()
        use def_inpout
        use def_master, only: lun_pdata 
        use mod_ecoute, only: ecoute, ecoute_set_read_unit
        use mod_messages, only: messages_live

        ! Read *.dat again looking if the PETSC option exists
        rewind(lun_pdata)
        do while( words(1) /= 'PETSC' )
            call ecoute('READ_PETSC',STOP_END_OF_FILE=.false.)
            if( words(1) == 'EOF  ' ) return
        end do

        ! Read petsc options
        if( option_not_off('PETSC') ) then
            kfl_petsc = .true.
            petsc_file_name = trim('petsc_config.in')
            call ecoute('READ_PETSC')
            do while( words(1) /= 'ENDPE' )
                if( words(1) == 'CONFI' )then
                    petsc_file_name = trim(&
                        getcha_long('CONFI','petsc_config.in','#color',len(petsc_file_name,KIND=ip)))
                endif
                call ecoute('READ_PETSC')
            enddo
         endif

    end subroutine alya2petsc_read_configuration


    subroutine alya2petsc_exchange()
        use mod_exchange, only: exchange_init, exchange_add, exchange_end

        call exchange_init()
        call exchange_add(kfl_petsc)
        call exchange_end()

        if( kfl_petsc )then
            call exchange_init()
            call exchange_add(petsc_file_name)
            call exchange_end()
        endif

    end subroutine alya2petsc_exchange


    subroutine alya2petsc_initialise()
        use def_master, only: IPARALL, INOTMASTER
        use mod_iofile, only: iofile
        use mod_parall, only: PAR_COMM_MY_CODE_WM
        use mod_messages, only: messages_live

        if( .not. kfl_petsc ) return

        call messages_live('START PETSC INITIALISATION')

        if( IPARALL )then
           call messages_live('   SET COMMUNICATOR: PAR_COMM_MY_CODE_WM')
#if defined USEMPIF08
           PETSC_COMM_WORLD = PAR_COMM_MY_CODE_WM % MPI_VAL
#else
           PETSC_COMM_WORLD = PAR_COMM_MY_CODE_WM 
#endif
        end if

        call messages_live('   READ CONFIGURATION FILE: '//trim(petsc_file_name))
        if( .not. IPARALL .or. INOTMASTER ) then
            call iofile(0_ip, petsc_file_unit, petsc_file_name, 'PETSc CONFIG FILE','old')
            call PetscInitialize(petsc_file_name,ierr)
            if( ierr > 0 )then
                call runend("ALYA2PETSC_INITIALIZATION: Could not initialize PETSc, aborting ... ")
            endif
            call iofile(2_ip, petsc_file_unit, petsc_file_name, 'PETSc CONFIG FILE')
        end if

        call messages_live('END PETSC INITIALISATION')

    end subroutine alya2petsc_initialise


    subroutine alya2petsc_finalise()
        use def_master, only: IPARALL, INOTMASTER

        if( .not. kfl_petsc ) return

        if( .not. IPARALL .or. INOTMASTER ) then
            call PetscFinalize( ierr )
        end if

    end subroutine alya2petsc_finalise


    subroutine alya2petsc_initialiseLinearSolver(alyapetsc)

        type(ZPLS), intent(inout) :: alyapetsc

        call petscKSP_initialise(alyapetsc)
  
    end subroutine alya2petsc_initialiseLinearSolver


    subroutine alya2petsc_createLinearSolver(alyapetsc, solve)
        use def_kintyp_solvers, only: soltyp
        use def_master, only: IPARALL, IMASTER, IEMPTY
        use mod_messages, only: messages_live

        type(ZPLS),      intent(inout) :: alyapetsc
        type(soltyp),    intent(inout) :: solve
  
        call messages_live('   GET DIMENIONS FROM ALYA SOLVER')
        call petscKSP_copyAlyaSolver(alyapetsc, IPARALL, IMASTER, IEMPTY, solve)
        call messages_live('   CREATE STRUCTURES (MPIAIJ, MPIVEC)')
        call petscKSP_allocateSystem(alyapetsc)
        call messages_live('   CREATE LINEAR SOLVER (KSP, PC)')
        call petscKSP_createSolver(alyapetsc)
  
    end subroutine alya2petsc_createLinearSolver
  
  
    subroutine alya2petsc_destroyLinearSolver(alyapetsc)

        type(ZPLS), intent(inout) :: alyapetsc

        ! Iside the finalise it is called allocateSystem & 
        ! destroySolver
        if( alyapetsc % I_WORKER )then
            call petscKSP_finalise(alyapetsc)
        endif

    end subroutine alya2petsc_destroyLinearSolver
  
  
    subroutine alya2petsc_solution(alyapetsc, amatr, rhsid, unkno)
  
        type(ZPLS), intent(inout)          :: alyapetsc  !< Solver type
        real(rp),   intent(inout), target  :: amatr(*)   !< Matrix
        real(rp),   intent(inout), target  :: rhsid(*)   !< RHS
        real(rp),   intent(inout), target  :: unkno(*)   !< Solution
        integer(ip)                        :: ii, jj
        PetscInt                           :: ai, aj, bj
        PetscScalar                        :: xx, yy
        PetscScalar,               pointer :: xx_v(:)
        PetscErrorCode                     :: ierr

        if( alyapetsc % I_WORKER )then
 
            ! Fill/assemble the PETSc-RHS
            call VecZeroEntries(alyapetsc % b, ierr)
            do jj = 1, alyapetsc % nunkn
                bj = alyapetsc % permr(jj) - 1 
                yy = rhsid(jj)
                call VecSetValue(alyapetsc % b, bj, yy, ADD_VALUES, ierr)
            enddo
            call VecAssemblyBegin(alyapetsc % b, ierr)
            call VecAssemblyEnd(alyapetsc % b, ierr)

            ! Initialize PETSc-MATRIX
            call MatZeroEntries(alyapetsc % A, ierr)
            do ii = 1, alyapetsc % nzmat
                ai = alyapetsc % ia_full_cooLex(ii)
                aj = alyapetsc % ja_full_cooLex(ii)
                xx = amatr(ii)
                call MatSetValue(alyapetsc % A, ai, aj, xx, ADD_VALUES, ierr)
            enddo
            call MatAssemblyBegin(alyapetsc % A, MAT_FINAL_ASSEMBLY, ierr)
            call MatAssemblyEnd(alyapetsc % A, MAT_FINAL_ASSEMBLY, ierr)

            ! Solve
            call VecZeroEntries(alyapetsc % x, ierr)
            call KSPSetUp(alyapetsc % ksp, ierr)
            call KSPSolve(alyapetsc % ksp, alyapetsc % b, alyapetsc % x, ierr)

            ! Recover solution to Alya matrix
            call VecGetArrayF90(alyapetsc % x,xx_v,ierr)
            do ii = 1, alyapetsc % nunkn_own
                unkno(ii) = xx_v(ii)
            enddo
            unkno(alyapetsc % nunkn_own+1:alyapetsc % nunkn) = 0.0_rp
            call VecRestoreArrayF90(alyapetsc % x, xx_v, ierr)

        endif
 
    end subroutine alya2petsc_solution
  
  
    subroutine alya2petsc_getOutInfo(alyapetsc, niters)
  
        type(ZPLS),  intent(inout) :: alyapetsc  !< Solver type
        integer(ip), intent(out)   :: niters
        PetscInt                   :: iters
        PetscErrorCode             :: ierr
  
        if( alyapetsc % I_WORKER )then
            call KSPGetIterationNumber( alyapetsc % ksp, iters, ierr);
            niters = int(iters,ip)
        endif
  
    end subroutine alya2petsc_getOutInfo
  
  
    subroutine alya2petsc_residualnorm(alyapetsc, resid, rnorm)
  
        type(ZPLS),    intent(inout)          :: alyapetsc
        real(rp),      intent(inout), target  :: resid(*)
        real(rp),      intent(out)            :: rnorm
        integer(ip)                           :: ii
        PetscScalar                           :: xnorm
        PetscScalar,                  pointer :: xv(:)
        PetscErrorCode                        :: ierr
        Vec                                   :: r
  
        if( alyapetsc % I_WORKER )then

            ! Allocate
            call VecDuplicate(alyapetsc % b, r, ierr)

            ! r = b - Ax
            call MatMult(alyapetsc % A, alyapetsc % x, r, ierr)
            call VecAXPY(r, MINUSONER, alyapetsc % b, ierr)
    
            ! Recover solution to Alya matrix
            call VecGetArrayF90( alyapetsc % x, xv, ierr)
            do ii = 1, alyapetsc % nequa_own
            resid(ii) = xv(ii)
            enddo
            call VecRestoreArrayF90(alyapetsc % x, xv, ierr)
    
            ! Compute L2 norm
            call VecNorm(r, NORM_2, xnorm, ierr)
            rnorm = real(xnorm,rp)

            ! Deallocate
            call VecDestroy(r, ierr)

        endif
        
    end subroutine alya2petsc_residualnorm
  
  
    subroutine alya2petsc_rhsidnorm(alyapetsc, rnorm)

        type(ZPLS),  intent(inout) :: alyapetsc
        real(rp),    intent(out)   :: rnorm
        PetscScalar                :: xnorm
        PetscErrorCode             :: ierr

        if( alyapetsc % I_WORKER )then

            call VecNorm(alyapetsc % b, NORM_2, xnorm, ierr)
            rnorm = real(xnorm,rp)

        endif
    
    end subroutine alya2petsc_rhsidnorm
  
  
    subroutine alya2petsc_unknonorm(alyapetsc, unorm)
  
        type(ZPLS),  intent(inout) :: alyapetsc
        real(rp),    intent(out)   :: unorm
        PetscScalar                :: xnorm
        PetscErrorCode             :: ierr

        if( alyapetsc % I_WORKER )then
            
            call VecNorm(alyapetsc % x, NORM_2, xnorm, ierr)
            unorm = real(xnorm,rp)
            
        endif
        
    end subroutine alya2petsc_unknonorm
  
  
    subroutine alya2petsc_interfaceNodeExchange(alyapetsc, rhsid, unkno)
        use mod_communications, only: PAR_INTERFACE_NODE_EXCHANGE
  
        type(ZPLS),  intent(inout)         :: alyapetsc
        real(rp),    intent(inout), target :: rhsid(*)
        real(rp),    intent(inout), target :: unkno(*) 
        integer(ip)                        :: ndofn
  
        ndofn = int(alyapetsc % ndofn, ip)
        call PAR_INTERFACE_NODE_EXCHANGE(ndofn, rhsid, 'SUM', 'IN MY CODE')
        call PAR_INTERFACE_NODE_EXCHANGE(ndofn, unkno, 'SUM', 'IN MY CODE')
  
    end subroutine alya2petsc_interfaceNodeExchange


    subroutine alya2petsc_solverInfoExchange(solve)
        use def_kintyp_solvers, only: soltyp
        use mod_communications, only: PAR_BROADCAST

        type(soltyp), intent(inout) :: solve

        call PAR_BROADCAST(solve % iters, wherein='IN MY CODE', root_rank=1_ip)

    end subroutine alya2petsc_solverInfoExchange


    subroutine alya2petsc_checkConsistency()
        use def_kermod, only: kfl_full_rows
        use mod_messages, only: messages_live
  
        if( kfl_full_rows /= 1_ip )then
           call messages_live("------------------------------------------------------------","WARNING")
           call messages_live(" PETSC Linear Solver (KSP) requires FULL_ROW format         ","WARNING") 
           call messages_live(" To active it add the following keywords in ker.dat file:   ","WARNING") 
           call messages_live(" NUMERICAL_TREATMENT","WARNING")
           call messages_live("    MESH","WARNING")
           call messages_live("       FULL_ROWS: ON","WARNING")
           call messages_live("    END_MESH","WARNING")
           call messages_live(" END_NUMERICAL_TREATMENT","WARNING")
           call messages_live("------------------------------------------------------------","WARNING")
           call runend('MOD_ALYA2PETSC: subroutine alya2petsc_checkConsistency')
        endif
  
    end subroutine alya2petsc_checkConsistency 
#endif

end module mod_alya2petsc
!> @}
