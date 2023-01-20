!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



#ifndef __ALYA_VECTOR_SIZE_DEF_HEADER__
#define __ALYA_VECTOR_SIZE_DEF_HEADER__


!------VECTOR_SIZE_OPENACC------!
#ifndef VECTOR_SIZE

#ifdef VECTOR_SIZE_OPENACC
#ifdef OPENACCHHH
#define VECTOR_SIZE VECTOR_SIZE_OPENACC
#endif
#endif

#endif

!------VECTOR_SIZE_CPU------!
#ifndef VECTOR_SIZE

#ifdef VECTOR_SIZE_CPU
#define VECTOR_SIZE VECTOR_SIZE_CPU
#else
#define VECTOR_SIZE 1_ip
#define VECTOR_SIZE_CPU VECTOR_SIZE
#endif

#else

#ifndef VECTOR_SIZE_CPU
#define VECTOR_SIZE_CPU VECTOR_SIZE
#endif

#endif

#endif
