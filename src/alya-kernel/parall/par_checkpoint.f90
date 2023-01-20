!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_checkpoint.f90
!> @author  houzeaux
!> @date    2020-05-08
!> @brief   MPI checkpointing
!> @details MPI checkpointing
!> @} 
!-----------------------------------------------------------------------

subroutine par_checkpoint()
  use def_master
  use def_parall
  use mod_parall
  implicit none
  integer(ip)             :: ichec
  
  if(nproc_par>1) then                                ! Checkpoint for communication
     call par_livinf(11_ip,' ',ichec)
     call par_chkpoi(ichec)
     if(ichec==0) then
        call par_livinf(12_ip,' ',ichec)
     else
        call runend('MPI IS NOT WORKING WELL')
     end if
  end if
  
end subroutine par_checkpoint
