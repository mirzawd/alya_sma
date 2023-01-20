!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kernel
!> @{
!> @file    mod_alya2signal.f90
!> @author  houzeaux
!> @date    2020-02-25
!> @brief   Signal handling
!> @details Module to treat signal handling... which is not part of
!>          the standard. See https://en.wikipedia.org/wiki/Signal_(IPC)
!-----------------------------------------------------------------------

module mod_alya2signal

  use def_master, only : ITASK_WRITE_RESTART
  use def_master, only : kfl_timei
  use def_master, only : kfl_stop
  use mod_run_config, only: run_config
#ifdef ALYA_SIGNAL
#ifdef __INTEL_COMPILER 
  USE IFPORT
#endif  
#endif
  implicit none
 

#ifndef __INTEL_COMPILER
  integer(4), parameter :: SIGINT  =  2_4 ! Interrupt from keyboard Ctrl-C 
  integer(4), parameter :: SIGTERM = 15_4 ! Termination signal             
#endif 
  integer(4), parameter :: SIGUSR1 = 10_4 ! SIGUSR1
  integer(4), parameter :: SIGUSR2 = 12_4 ! SIGUSR2
  
contains 

  subroutine alya2signal

#if defined(ALYA_SIGNAL) && (defined(__INTEL_COMPILER) || defined(__GNUC__))
    integer(4) :: status
#endif
    integer(4) :: flag

    flag=-1_4

#ifdef ALYA_SIGNAL
#if defined __INTEL_COMPILER 
    status = signal(SIGINT  , alya2signal_intel  , flag )
    status = signal(SIGTERM , alya2signal_intel  , flag ) 
    status = signal(SIGUSR1 , alya2signal_intel  , flag )
    status = signal(SIGUSR2 , alya2signal_intel  , flag )
#elif defined __PGI

#elif defined __ibmxl__
    call signal(SIGINT  , alya2signal_sigint  )     
    call signal(SIGTERM , alya2signal_sigterm )   
#elif defined __GNUC__
    status = signal(SIGINT  , alya2signal_sigint  )
    status = signal(SIGTERM , alya2signal_sigterm )
    status = signal(SIGUSR1 , alya2signal_sigusr1 )
    status = signal(SIGUSR2 , alya2signal_sigusr2 )
#endif
#endif

  end subroutine alya2signal

  integer(4) function alya2signal_intel(sig_num)

    integer(4) :: sig_num

    alya2signal_intel = 1

#ifdef ALYA_SIGNAL
#ifdef __INTEL_COMPILER 
    select case ( sig_num )
    case ( SIGINT  ) ; call alya2signal_sigint()
    case ( SIGTERM ) ; call alya2signal_sigterm()
    case ( SIGUSR1 ) ; call alya2signal_sigusr1()
    case ( SIGUSR2 ) ; call alya2signal_sigusr2()
    end select 
#endif
#endif
    
  end  function alya2signal_intel

  subroutine alya2signal_sigint()

    print *,'SIGNAL SIGINT'
    run_config%restart%preliminary = .true.
    kfl_stop  = 1

  end subroutine alya2signal_sigint

  subroutine alya2signal_sigterm()

    print *,'SIGNAL SIGTERM'
    run_config%restart%preliminary = .true.
    kfl_stop  = 1

  end subroutine alya2signal_sigterm

  subroutine alya2signal_sigusr1()
    
    print *,'SIGNAL SIGUSR1: Alya will stop after current time-step, please wait'
    run_config%restart%preliminary = .true.
    kfl_stop  = 1

  end subroutine alya2signal_sigusr1
  
  subroutine alya2signal_sigusr2()
    
    print *,'SIGNAL SIGUSR2: Alya will stop after current time-step, please wait'
    run_config%restart%preliminary = .true.
    kfl_stop  = 1

  end subroutine alya2signal_sigusr2
  
end module mod_alya2signal
!> @}
