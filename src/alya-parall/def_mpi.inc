!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



#ifndef __ALYA_MPI_DEF_HEADER__
#define __ALYA_MPI_DEF_HEADER__

#ifdef USEMPIF08
#define MY_MPI_COMM     type(MPI_Comm)
#define MY_MPI_STATUS   type(MPI_Status)
#define MY_MPI_FILE     type(MPI_File)
#define MY_MPI_REQUEST  type(MPI_request)
#define MY_MPI_DATATYPE type(MPI_Datatype)
#define MY_MPI_INFO     type(MPI_Info)
#define MY_MPI_WIN      type(MPI_Win)
#else
#define MY_MPI_COMM     integer(4)
#define MY_MPI_STATUS   integer(4)
#define MY_MPI_FILE     integer(4)
#define MY_MPI_REQUEST  integer(4)
#define MY_MPI_DATATYPE integer(4)
#define MY_MPI_INFO     integer(4)
#define MY_MPI_WIN      integer(4)
#endif

#endif
