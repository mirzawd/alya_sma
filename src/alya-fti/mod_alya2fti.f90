!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!

!-----------------------------------------------------------------------
!> @addtogroup Parall
!> Bridge to FTI
!> @{
!> @file    mod_alya2fti.f90
!> @author  jchiva
!> @date    2020-06-19
!> @brief   Bridge to FTI
!> @details Interfaces with FTI
!-----------------------------------------------------------------------

module mod_alya2fti

#ifdef ALYA_FTI

    use def_master, only: kfl_rstar, ittim, mitim
    use def_kintyp, only: ip, lg
    use def_inpout, only: words
    use def_communications, only: MPI_COMM_WORLD
    use mod_messages, only: messages_live
    use FTI

    implicit none

    integer, target :: FTI_comm_world           ! FTI COMMUNICATOR
    integer(ip)     :: FTI_st = 0, &    ! FTI status (if 0 /= FTI_st then a previous FTI file exists and will be read
                       FTI_write_ckpt = 0, &    ! A checkpoint will be written using FTI in this step
                       FTI_DBG                  ! Is FTI debug

    logical(lg)     :: FTI_usage                ! Is FTI enabled

    public :: FTI_comm_world
    public :: FTI_st
    public :: FTI_write_ckpt
    public :: FTI_usage
    public :: FTI_DBG
    public :: alya2fti_initialization
    public :: alya2fti_finalization
    public :: alya2fti_iexcha
    public :: alya2fti_RecoverVarInit
    public :: alya2fti_RecoverVarFinalize
    public :: alya2fti_InitICP
    public :: alya2fti_FinalizeICP

contains

    !-----------------------------------------------------------------------
    !>
    !> @author  bsc21943
    !> @date    2019-07-19
    !> @brief   Initialization
    !> @details Initialization of FTI
    !>
    !-----------------------------------------------------------------------

    subroutine alya2fti_initialization()
        integer, target :: err
#if defined USEMPIF08
        FTI_comm_world = MPI_COMM_WORLD%MPI_VAL
#else
        FTI_comm_world = MPI_COMM_WORLD
#endif
        call FTI_Init("config.fti", FTI_comm_world, err)
        call FTI_status(FTI_st)
    end subroutine alya2fti_initialization

    subroutine alya2fti_finalization()
        integer(4) :: ierror
        ierror = 1
        call FTI_Finalize(ierror) !ierror is output, should be used for something?
    end subroutine alya2fti_finalization

    subroutine alya2fti_iexcha()
        call iexcha(FTI_usage)
        call iexcha(FTI_DBG)
    end subroutine alya2fti_iexcha

    subroutine alya2fti_RecoverVarInit()
        if (FTI_st /= 0) then
            call FTI_RecoverVarInit()
        end if
    end subroutine alya2fti_RecoverVarInit

    subroutine alya2fti_RecoverVarFinalize()
        if (FTI_st /= 0 .and. kfl_rstar == 2) then
            call FTI_RecoverVarFinalize()
        end if
    end subroutine alya2fti_RecoverVarFinalize

    subroutine alya2fti_InitICP()
        integer(ip) :: FTI_error
        if (FTI_usage .and. (FTI_DBG /= 1 .or. ittim < mitim)) then
            FTI_write_ckpt = 1_ip
            call FTI_initICP(ittim, 4, .true., FTI_error)
        end if
    end subroutine alya2fti_InitICP

    subroutine alya2fti_FinalizeICP()
        integer(ip) :: FTI_error
        if (FTI_write_ckpt /= 0) then
            FTI_write_ckpt = 0
            call FTI_FinalizeICP(FTI_error)
            if (FTI_error /= 0) then
                call messages_live('FTI: FAILED TO WRITE CHECKPOINT FILE', 'WARNING')
            end if
        end if
    end subroutine alya2fti_FinalizeICP

#endif

end module mod_alya2fti
!> @}
