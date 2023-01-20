!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kinds_and_types
!> @{
!> @file    def_kintyp_petsc.f90
!> @author  Adria Quintanas-Corominas
!> @brief   Bridge to PETSc
!-----------------------------------------------------------------------

module def_kintyp_petsc 


#ifdef PETSC
#include <petsc/finclude/petscksp.h>

    use def_kintyp_basic, only: lg, ip, rp
    use petscksp

    implicit none

    type PETSC_LINEAR_SOLVER
        character(10)           :: name
        logical(lg)             :: I_PARALL
        logical(lg)             :: I_WORKER
        logical(lg)             :: I_WORKER_NOT_EMPTY
        character(6)            :: PETSC_STRATEGY = 'mpiaij'
        character(6)            :: PETSC_MATRIX_TYPE = 'mpiaij'
        character(6)            :: PETSC_VECTOR_TYPE = 'vecmpi'
        PetscInt                :: ndofn
        PetscInt                :: nequa
        PetscInt                :: nequa_own
        PetscInt                :: nunkn
        PetscInt                :: nunkn_own
        PetscInt                :: nzmat
        PetscInt                :: nzmat_own
        PetscInt,     pointer   :: ia(:)
        PetscInt,     pointer   :: ja(:)
        PetscInt,     pointer   :: ia_full(:)
        PetscInt,     pointer   :: ja_full(:)
        PetscInt,     pointer   :: ia_full_cooLex(:)
        PetscInt,     pointer   :: ja_full_cooLex(:)
        PetscInt,     pointer   :: permr(:)
        PetscInt,     pointer   :: nndiag(:)
        PetscInt,     pointer   :: nnoutd(:)
        Mat                     :: A
        Vec                     :: b
        Vec                     :: x
        KSP                     :: ksp
        PC                      :: pc
    end type PETSC_LINEAR_SOLVER

    private

    public :: PETSC_LINEAR_SOLVER

#endif

end module def_kintyp_petsc 
!> @}
