!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup mpio
!> @{
!> @file    mod_mpio_par_configure.f90
!> @author  Damien Dosimont
!> @date    27/09/2017
!> @brief   MPI-IO Wrapper
!> @details This module is a wrapper of MPI-IO functions.
!>          It adds error management and makes type management transparent to the user.
!> @}
!-----------------------------------------------------------------------


module mod_mpio_par_configure

    use def_kintyp_basic,               only : ip,rp,lg,r1p
    use def_kintyp_comm,                only : comm_data_par
    use def_master,                     only : intost, retost
    use def_master,                     only : IPARSLAVE, IMASTER, INOTMASTER, IPARALL, ISEQUEN
    use def_master,                     only : lun_pos00
    use def_master,                     only : kfl_ptask
    use mod_memory,                     only : memory_alloca, memory_deallo
    use def_domain,                     only : ndime, npoin, nboun, nelem, npoin_2, nelem_2, nboun_2
    use mod_messages,                   only : livinf
    use def_mpio
    use mod_mpio_files
    use def_mpi
    use mod_mpio_config,                only : mpio_config
    
#ifndef MPI_OFF
    use mod_communications,             only : PAR_DEFINE_COMMUNICATOR
#endif

    implicit none

    private

    character(150)                          ::  wherein_code="IN MY CODE"           !

    integer(ip)                             ::  dboxc_p(3)      ! #bin boxes per direction in the coarse bin



    public                                  ::  mpio_init,&
                                                par_mpio_finalize,&
                                                dboxc_p

    contains

      subroutine mpio_init()
        use mod_communications,     only : PAR_BROADCAST

#ifndef MPI_OFF
        call PAR_BROADCAST(kfl_ptask)
#endif
        mpio_memor=0_8

        if (IPARALL) then
            call par_mpio_init()
        end if
        if (mpio_config%input%enabled) then
            call find_mpio_mesh_name()
        end if

      end subroutine mpio_init

    subroutine par_mpio_init()
#ifndef MPI_OFF
        if (mpio_config%enabled) then
            call livinf(-4_ip,'INIT MPI-IO COMMUNICATORS',0_ip)
            call par_create_communicators_all()
            call livinf(-5_ip,'END INIT MPI-IO COMMUNICATORS',0_ip)
        end if
#endif
    end subroutine

    subroutine par_mpio_finalize()
#ifndef MPI_OFF
        use mod_mpio_par_async_io
        implicit none
        call PAR_ASYNCHRONOUS_WRITE_POP_ALL(0_ip)
#endif
    end subroutine

    subroutine par_create_communicators_all()
#ifndef MPI_OFF
        use mod_communications, only : PAR_COMM_SPLIT, PAR_BROADCAST
        use mod_parall,         only : PAR_CODE_SIZE, PAR_MY_CODE_RANK
        use mod_parall,         only : PAR_COMM_MPIO, PAR_COMM_MPIO_WM, PAR_COMM_MPIO_RANK_WM, PAR_COMM_MPIO_WM_SIZE
        use def_master,         only : IIOSLAVE

        implicit none
        integer(ip)         :: icolor, iproc
        character(100), PARAMETER :: vacal = "par_create_communicators_all"

        PAR_COMM_MPIO_WM_SIZE = int((PAR_CODE_SIZE-1_ip), ip)

        if (PAR_COMM_MPIO_WM_SIZE<2) then
            call livinf(0_ip,"ONLY ONE MPI-IO WORKER: SWITCHING TO SEQUENTIAL VERSION",0_ip)
            return
        end if

        call livinf(0_ip,"NUMBER OF MPI-IO WORKERS (GEOMETRY READING): "//intost(PAR_COMM_MPIO_WM_SIZE),0_ip)
        call livinf(0_ip,"NUMBER OF MPI-IO WORKERS (POST/RST): "//intost(PAR_CODE_SIZE-1),0_ip)
        !
        ! Create the slaves communicator
        !
        IIOSLAVE = .FALSE.
        if (INOTMASTER .and. (PAR_MY_CODE_RANK-1 < PAR_COMM_MPIO_WM_SIZE)) then
            IIOSLAVE = .TRUE.
            !print*, "process: ", PAR_MY_CODE_RANK
        end if

        icolor = 0
        PAR_COMM_MPIO_WM = PAR_COMM_NULL
        PAR_COMM_MPIO    = PAR_COMM_NULL
        
        !
        ! Creates partitioners and master communicator
        !
        icolor = 0
        if( IIOSLAVE .or. IMASTER ) then
           icolor = 1_4
        else
           icolor = 0_4
        end if
        call PAR_COMM_SPLIT(icolor,PAR_COMM_MPIO,iproc,wherein_code)
        
        if( IIOSLAVE ) then
           icolor = 1_4
        else
           icolor = 0_4
        end if
        call PAR_COMM_SPLIT(icolor,PAR_COMM_MPIO_WM,iproc,wherein_code)
        PAR_COMM_MPIO_RANK_WM = int(iproc,ip)
        
#endif
    end subroutine

end module
