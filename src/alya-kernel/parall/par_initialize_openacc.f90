!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_initialize_openacc.f90
!> @author  Guillaume Houzeaux
!> @date    1789-05-05
!> @brief   Initialize OpenACC
!> @details Maps processes and devices
!>
!> @}
!------------------------------------------------------------------------

subroutine par_initialize_openacc()
#ifdef _OPENACC

  use def_kintyp,   only : ip
  use def_master,   only : intost, kfl_paral
  use def_parall,   only : kfl_openacc_policy
  use def_parall,   only : kfl_openacc_ranks_per_dev
  use def_parall,   only : kfl_openacc_contiguous_ranks_per_dev
  use def_parall,   only : kfl_openacc_streams_per_dev
  use def_parall,   only : kfl_openacc_multithreading
  use def_parall,   only : OPENACC_POLICY_1PPDEV
  use def_parall,   only : OPENACC_POLICY_NPPDEV
  use def_parall,   only : OPENACC_POLICY_ALL2DEV

  use mod_parall,   only : PAR_COMM_SHARED_WM
  use mod_parall,   only : PAR_SHARED_RANK_WM

  use openacc,      only : acc_device_host, acc_device_not_host, acc_get_num_devices
  use mod_communications_tools, only : PAR_COMM_RANK_AND_SIZE

  implicit none

  integer(ip)           :: ndevs, not_host_devs, host_devs, processes_to_devs, dev_num
  integer(ip)           :: rank_on_host, host_comm_size

  ndevs = 0
  not_host_devs = 0
  host_devs = 0
  processes_to_devs = 1

  if (kfl_paral < 0 .or. kfl_paral > 0) then
    rank_on_host = PAR_SHARED_RANK_WM
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_SHARED_WM,rank_on_host,host_comm_size)

    ! get the number of devices, host and not host
    not_host_devs = acc_get_num_devices(acc_device_not_host)
    ndevs = ndevs + not_host_devs

    ! The OpenACC code must be compiled for multithreding if you want to use this
    if (kfl_openacc_multithreading) then
      host_devs = acc_get_num_devices(acc_device_host)
      ndevs = ndevs + host_devs
    end if

    select case (kfl_openacc_policy)
    case (OPENACC_POLICY_NPPDEV) ! Several processes per device
      processes_to_devs = ndevs * kfl_openacc_ranks_per_dev
    case (OPENACC_POLICY_ALL2DEV) ! All the processes use a device
      processes_to_devs = host_comm_size
    case DEFAULT ! (OPENACC_POLICY_1PPDEV) ! One process per device
      processes_to_devs = ndevs
    end select

    dev_num = mod(rank_on_host/kfl_openacc_contiguous_ranks_per_dev, processes_to_devs)

    ! Set and initialize the device for openacc parallel code
    ! Fill up the GPUs first
    if (dev_num < not_host_devs) then
      call acc_init_device(dev_num, acc_device_not_host)
    else if (kfl_openacc_multithreading .and. dev_num-not_host_devs < host_devs) then
      call acc_init_device(dev_num, acc_device_host)
    end if ! The remaining processes will run serially on a cpu core

  end if

#endif
end subroutine par_initialize_openacc
