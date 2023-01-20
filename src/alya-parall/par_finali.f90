!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_finali.f90
!> @author  Guillaume Houzeaux
!> @date    01/02/2014
!> @brief   Finalize MPI after a call to runend
!> @details Finalize MPI. Abort MPI to kill all slaves if code
!>          is not ending properly
!> @} 
!----------------------------------------------------------------------

subroutine par_finali(kfl_endok)
  use def_kintyp, only   :  ip
  use mod_parall, only   :  PAR_COMM_WORLD
  implicit none
  integer(ip)            :: kfl_endok
  integer(4)             :: istat,ierro

#ifndef MPI_OFF
  include  'mpif.h'
#endif

  ierro = 1 

#ifndef MPI_OFF
  if ( kfl_endok == 1_ip ) then
     call MPI_Finalize(istat)
  else
     call MPI_Abort(PAR_COMM_WORLD,ierro,istat)
  end if
#endif

end subroutine par_finali
