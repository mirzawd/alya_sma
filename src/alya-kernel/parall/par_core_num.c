/*-----------------------------------------------------------------*/
/* Copyright 2005 - 2022 Barcelona Supercomputing Center.          */
/* Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE */
/* for nonprofit scientific purposes only.                         */
/* See companion file LICENSE.txt.                                 */
/*-----------------------------------------------------------------*/



#include <stdio.h>
#ifdef __unix__
#include <sched.h>
#endif

void par_core_num( int *core_num ) {
  
  int thread_num;

#ifdef __unix__
#ifdef _OPENMP
#pragma omp parallel private(thread_num) shared(core_num)
  {
    thread_num = omp_get_thread_num();
    core_num[thread_num] = sched_getcpu();
  }
#else
  core_num[0] = sched_getcpu();
#endif
#else 
  core_num[0] = -1;
#endif
}
