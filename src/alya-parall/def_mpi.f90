!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Communication_Toolbox
!> @{
!> @file    def_mpi.f90
!> @author  houzeaux
!> @date    2022-08-30
!> @brief   MPI
!> @details MPI definitions
!-----------------------------------------------------------------------

module def_mpi

  use def_kintyp_basic, only : ip,rp
#include "def_mpi.inc"

  !----------------------------------------------------------------------
  !
  ! Fortran MPI bindings
  !
  !----------------------------------------------------------------------

#ifndef MPI_OFF
#if defined USEMPIF08
  !
  ! MPI module
  !
  use mpi_f08
#elif defined  MPIFH
  !
  ! MPI include
  !
  include 'mpif.h'  
#else
  !
  ! MPI include (default)
  !
  include 'mpif.h'
#endif
#else
  !
  ! No MPI... MPI variables and status
  !
  integer, parameter :: MPI_SUCCESS      =  0
  integer, parameter :: MPI_REAL4        =  0
  integer, parameter :: MPI_REAL8        =  0
  integer, parameter :: MPI_INTEGER4     =  0
  integer, parameter :: MPI_INTEGER8     =  0
  integer, parameter :: MPI_STATUS_SIZE  =  1
  integer, parameter :: MPI_FILE_NULL    = -1
  integer, parameter :: MPI_COMM_NULL    = -1
  integer, parameter :: MPI_WIN_NULL     = -1
  integer, parameter :: MPI_INFO_NULL    = -1
  integer, parameter :: MPI_ADDRESS_KIND =  8
#endif

  !----------------------------------------------------------------------
  !
  ! Null communicator and generic status
  !
  !----------------------------------------------------------------------
  
  MY_MPI_COMM     :: PAR_COMM_NULL
  
#ifdef USEMPIF08
  MY_MPI_STATUS   :: status  
#else
  MY_MPI_STATUS   :: status(MPI_STATUS_SIZE)
#endif
  
end module def_mpi
!> @}




  
