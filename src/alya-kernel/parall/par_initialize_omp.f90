!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_initialize_omp.f90
!> @author  Guillaume Houzeaux
!> @date    13/11/2006
!> @brief   Initialize OpenMP
!> @details Prints the number of OpenMP threads
!>          Define the openmp with coloring option by default
!>
!> @}
!------------------------------------------------------------------------

subroutine par_initialize_omp()

  use def_kintyp,   only : ip
  use def_master,   only : intost
#if defined(ALYA_OMPSS) || defined(_OPENMP)
  use mod_parall,   only : par_omp_num_threads
#endif
  use mod_parall,   only : par_hybrid
  use mod_parall,   only : PAR_OPENMP_COLORING
  use mod_parall,   only : PAR_OPENMP_NO_COLORING
  use mod_parall,   only : PAR_OMPSS
  use mod_parall,   only : PAR_HYBRID_OFF
  use mod_messages, only : livinf
#ifdef _OPENMP  
  use omp_lib
#endif 
  implicit none
  !integer(4)  :: OMP_GET_MAX_THREADS
  !integer(4)  :: OMP_GET_NUM_THREADS
  !integer(4)  :: OMP_GET_THREAD_NUM
#if defined(ALYA_OMPSS) || defined(_OPENMP)
  integer(4)  :: num_nthre4,max_nthre4
  integer(ip) :: num_nthre,max_nthre
#endif
  !
  ! Get OMP MAX_THREADS
  ! When using OmpSs, this is the max. number of threads available on the node
  ! When using normal OMP, this is equivalent to OMP_NUM_THREADS
  !
#if defined ALYA_OMPSS
  !$OMP PARALLEL
  !$OMP SINGLE
#ifdef _OPENMP  
  num_nthre4  = OMP_GET_NUM_THREADS()
  max_nthre4  = OMP_GET_MAX_THREADS() 
#else
  num_nthre4  = 1_4
  max_nthre4  = 1_4
#endif  
  num_nthre  = int(num_nthre4,ip)
  max_nthre  = int(max_nthre4,ip)
  par_omp_num_threads = num_nthre
  !$OMP END SINGLE
  !$OMP END PARALLEL
  if( par_omp_num_threads == 0 ) call runend('INIOMP: OMP_NUM_THREADS SHOULD BE /= 0 ')
  call livinf(0_ip,'OMPSS: NUMBER OF THREADS=     '//trim(intost(num_nthre)),0_ip)
  call livinf(0_ip,'OMPSS: MAX NUMBER OF THREADS= '//trim(intost(max_nthre)),0_ip)
#elif defined _OPENMP 
  !$OMP PARALLEL
  !$OMP SINGLE
  num_nthre4  = OMP_GET_NUM_THREADS()
  max_nthre4  = OMP_GET_MAX_THREADS()
  num_nthre   = int(num_nthre4,ip)
  max_nthre   = int(max_nthre4,ip)
  par_omp_num_threads = num_nthre
  !$OMP END SINGLE
  !$OMP END PARALLEL
  if( par_omp_num_threads == 0 ) call runend('INIOMP: OMP_NUM_THREADS SHOULD BE /= 0 ')
  call livinf(0_ip,'OPENMP: NUMBER OF THREADS=     '//trim(intost(num_nthre)),0_ip)
  call livinf(0_ip,'OPENMP: MAX NUMBER OF THREADS= '//trim(intost(max_nthre)),0_ip)
#endif
  !
  ! Define hybrid method
  !
#ifdef ALYA_OMPSS                      /* OMPSS */
  par_hybrid = PAR_OMPSS
#elif _OPENMP && defined NO_COLORING   /* OMP without coloring */
  par_hybrid = PAR_OPENMP_NO_COLORING
#elif _OPENMP                          /* OMP with coloring */
  par_hybrid = PAR_OPENMP_COLORING
#else                                  /* Hybrid off */
  par_hybrid = PAR_HYBRID_OFF
#endif

end subroutine par_initialize_omp

