!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Timste
!> @{
!> @file    Timste.f90
!> @author  Guillaume Houzeaux
!> @brief   Compute time step
!> @details Compute time step:
!>          - Compute the time step for each module
!>          - Compute the global time step 
!>
!> @}
!-----------------------------------------------------------------------
subroutine Timste()

  use def_kintyp_basic, only : ip
  use def_master,       only : ITASK_TIMSTE
  use def_master,       only : iblok,nblok
  use def_master,       only : dtinv
  use mod_messages,     only : livinf
  use mod_messages,     only : messages_live
  use mod_couplings,    only : couplings_time_step
  use mod_moduls,       only : moduls 
  implicit none
  !
  ! Initializations
  !
  call livinf(4_ip, ' ',0_ip)
  call iniste(1_ip)
  !
  ! Do things if we are repartitioning, doign AMR, etc.
  !
  call what_is_being_done()
  !
  ! Live information
  ! 
  !call messages_live('MODULES COMPUTE CRITICAL TIME STEPS')
  !
  ! Compute the time step for each module
  !     
  call Kermod(ITASK_TIMSTE)
  do iblok = 1,nblok
     call moduls(ITASK_TIMSTE)
  end do
  !
  ! Coupling time step
  !
  call couplings_time_step(dtinv)
  !
  ! Initializations
  !
  call iniste(2_ip)
  !
  ! Computes the global time step
  !
  call setgts(ITASK_TIMSTE)
  call livinf(18_ip,' ',0_ip)

end subroutine Timste

