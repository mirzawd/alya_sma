!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_initialize_mpi.f90
!> @author  houzeaux
!> @date    2020-05-08
!> @brief   Start Parall 
!> @details Turn on parall
!> @} 
!-----------------------------------------------------------------------

subroutine par_initialize_mpi()
  use def_kintyp_basic,         only : rp
  use mod_outfor,               only : outfor
  use mod_communications_tools, only : PAR_COMM_SET_ERRHANDLER
  use def_parall
  use mod_parall
  real(rp) :: time1,time2
  
  call cputim(time1)
  call par_initia()                                      ! Initialize MPI
  call par_errors(1_ip)
  call cputim(time2)
  cpu_paral(1)=time2-time1

  call PAR_COMM_SET_ERRHANDLER()
  
end subroutine par_initialize_mpi
