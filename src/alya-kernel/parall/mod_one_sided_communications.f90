!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @defgroup Communication_Toolbox
!> Toolbox for MPI communication, bridge to MPI
!> @{
!> @name    Parallelization toolbox
!> @file    mod_one_sided_communications.f90
!> @author  Guillaume Houzeaux
!> @date    28/09/2023
!> @brief   ToolBox for parallel one-sided communications
!> @details ToolBox for parallel one-sided communications
!------------------------------------------------------------------------

module mod_one_sided_communications

  use def_kintyp_basic,                  only : ip,rp
  use def_master,                        only : INOTMASTER
  use def_domain,                        only : ndime,npoin,nelem,nboun
  use mod_memory,                        only : memory_alloca
  use mod_communications_point_to_point, only : PAR_SEND_RECEIVE
  use def_kintyp_comm,                   only : COMM_ONE_SIDED
  use mod_parall
  use def_mpi
  use, intrinsic :: ISO_C_BINDING
#include "def_mpi.inc"

  implicit none
  private
  
  public :: par_one_sided_initialization
  public :: par_one_sided_allocate

  contains

    !-----------------------------------------------------------------------
    !> 
    !> @author  guillaume
    !> @date    2022-09-28
    !> @brief   Allocate window
    !> @details Allocate window
    !> 
    !-----------------------------------------------------------------------

    subroutine par_one_sided_initialization()
      
    end subroutine par_one_sided_initialization
    
    !-----------------------------------------------------------------------
    !> 
    !> @author  guillaume
    !> @date    2022-09-28
    !> @brief   Allocate window
    !> @details Allocate window
    !> 
    !-----------------------------------------------------------------------

    subroutine par_one_sided_allocate()

      MY_MPI_INFO                    :: info
      integer(kind=MPI_ADDRESS_KIND) :: my_size
      type(C_PTR)                    :: cptr_buf
      integer(4)                     :: disp_unit
      integer(4)                     :: ierro
      integer(ip)                    :: sizebuff,ii,kk

#ifndef MPI_OFF
      if( INOTMASTER .and. par_one_sided /= 0 ) then
         !
         ! Window
         !
         commd % type = COMM_ONE_SIDED
         sizebuff     = ndime * commd % bound_dim                ! Maximum expected dimension for buffer
         my_size      = int(kind(commd % buffe_recv)*sizebuff,8) ! Buffer size in bytes
         disp_unit    = kind(commd % buffe_recv)                 ! Units
         info         = MPI_INFO_NULL                            ! Default options

         call MPI_Win_allocate(my_size, disp_unit, info, PAR_COMM_MY_CODE_WM, cptr_buf, commd % WIN_WM , ierro)
         !
         ! Buffers
         !
         CALL C_F_POINTER(cptr_buf, commd % buffe_recv, (/sizebuff/) )

         call memory_alloca(par_memor,'displ_recv','par_one_sided_allocate',commd % displ_recv,commd % nneig)
         do ii = 1,commd % nneig
            kk = commd % bound_size(ii)-1
            call PAR_SEND_RECEIVE(kk,commd % displ_recv(ii),dom_i=commd % neights(ii))
         end do

      end if
#endif

    end subroutine par_one_sided_allocate

    
end module mod_one_sided_communications
