!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    mor_comm_basic.f90
!> @author  houzeaux
!> @date    2021-06-21
!> @brief   Simple communication
!> @details Communication using the basic communicator
!>
!-----------------------------------------------------------------------

module mod_comm_basic

  use def_kintyp_basic, only : ip,rp
  use def_kintyp_comm,  only : comm_data_par_basic
  use mod_memory_basic, only : memory_size
  use mod_memory_basic, only : memory_alloca
  use mod_memory_basic, only : memory_deallo
  use mod_parall,       only : PAR_REAL
  use def_mpi
#include "def_mpi.inc"
  implicit none
  private

  interface par_interface_exchange
     module procedure par_interface_exchange_1 
  end interface par_interface_exchange

  public :: par_interface_exchange

contains
  
  subroutine par_interface_exchange_1(xx,comm)

    real(rp),                  pointer,  intent(inout) :: xx(:)
    class(comm_data_par_basic),          intent(in)    :: comm
    integer(ip)                                        :: ii,nsize,jj
    integer(ip)                                        :: ipoin,ini,kk
    integer(4)                                         :: istat4,nsize4
    integer(4)                                         :: dom_i4
    MY_MPI_COMM                                        :: PAR_COMM_WORLD
    real(rp),                  pointer                 :: sendbuff_rp(:)
    real(rp),                  pointer                 :: recvbuff_rp(:)
    integer(8)                                         :: memor_loc(2)
    
#ifndef MPI_OFF

    PAR_COMM_WORLD = comm % PAR_COMM_WORLD
    
    if( PAR_COMM_WORLD /= PAR_COMM_NULL ) then
       !
       ! Allocate memory
       !
       nullify(sendbuff_rp)
       nullify(recvbuff_rp)

       call memory_alloca(memor_loc,'SENDBUFF_RP','par_interface_node_exchange_rp',sendbuff_rp,comm % bound_dim,'DO_NOT_INITIALIZE')
       call memory_alloca(memor_loc,'RECVBUFF_RP','par_interface_node_exchange_rp',recvbuff_rp,comm % bound_dim,'DO_NOT_INITIALIZE')
       !
       ! Save in temp_send
       !
       kk = 0
       do jj = 1,comm % bound_dim
          ipoin = comm % bound_perm(jj)
          kk              = kk + 1
          sendbuff_rp(kk) = xx(ipoin)
          recvbuff_rp(kk) = 0.0_rp
       end do
       !
       ! Send    temp_send
       ! Receive temp_recv
       !
       istat4 = 0_4
       do ii = 1,comm % nneig

          dom_i4 = int(comm % neights(ii),4)
          ini    = comm % bound_size(ii)   
          nsize  = comm % bound_size(ii+1) - ini          
          nsize4 = int(nsize,4)

          call MPI_Sendrecv(                       &
               sendbuff_rp(ini:), nsize4,          &
               PAR_REAL, dom_i4, 0_4,             &
               recvbuff_rp(ini:), nsize4,          &
               PAR_REAL, dom_i4, 0_4,             &
               PAR_COMM_WORLD, status , istat4     )
 
          if( istat4 /= MPI_SUCCESS ) stop 1

       end do

       kk = 0
       do jj = 1,comm % bound_dim
          ipoin     = comm % bound_perm(jj)
          kk        = kk + 1
          xx(ipoin) = xx(ipoin) + recvbuff_rp(kk)
       end do

       call memory_deallo(memor_loc,'SENDBUFF_RP','par_interface_node_exchange_rp',sendbuff_rp)
       call memory_deallo(memor_loc,'RECVBUFF_RP','par_interface_node_exchange_rp',recvbuff_rp)

    end if

#endif

  end subroutine par_interface_exchange_1

!!$  subroutine par_interface_exchange_2(xx,comm)
!!$
!!$    class(*),                  pointer,  intent(inout) :: xx(:,:)
!!$    type(comm_data_par_basic),           intent(in)    :: comm
!!$    integer(ip)                                        :: ii,nsize,jj,dom_i,dom_min
!!$    integer(ip)                                        :: ipoin,ini,kk,idofn
!!$    integer(4)                                         :: istat4,nsize4,count4
!!$    integer(4)                                         :: dom_i4,ndofn
!!$    integer(ip)                                        :: bound_dim
!!$    integer(4)                                         :: PAR_COMM_WORLD4
!!$    real(rp),                  pointer                 :: sendbuff_rp(:)
!!$    real(rp),                  pointer                 :: recvbuff_rp(:)
!!$    integer(ip),               pointer                 :: sendbuff_ip(:)
!!$    integer(ip),               pointer                 :: recvbuff_ip(:)
!!$    integer(4)                                         :: PAR_INTEGER4
!!$    integer(8)                                         :: memor_loc(2)
!!$    
!!$#ifndef MPI_OFF
!!$
!!$    select type ( xx )
!!$    type is ( real    (kind=rp) ) ; ndofn = memory_size(xx,1_ip)
!!$    type is ( integer (kind=ip) ) ; ndofn = memory_size(xx,1_ip)
!!$    end select
!!$    
!!$    PAR_COMM_WORLD4 = int(comm % PAR_COMM_WORLD,4)
!!$    
!!$    if( PAR_COMM_WORLD4 /= 0 ) then
!!$       !
!!$       ! Allocate memory
!!$       !
!!$       nullify(sendbuff_rp)
!!$       nullify(recvbuff_rp)
!!$       nullify(sendbuff_ip)
!!$       nullify(recvbuff_ip)
!!$       
!!$       select type ( xx )
!!$       type is ( real (kind=rp) )
!!$          call memory_alloca(memor_loc,'SENDBUFF_RP','par_interface_node_exchange_rp',sendbuff_rp,bound_dim * ndofn,'DO_NOT_INITIALIZE')
!!$          call memory_alloca(memor_loc,'RECVBUFF_RP','par_interface_node_exchange_rp',recvbuff_rp,bound_dim * ndofn,'DO_NOT_INITIALIZE')
!!$       type is ( integer (kind=ip) )
!!$          call memory_alloca(memor_loc,'SENDBUFF_IP','par_interface_node_exchange_rp',sendbuff_ip,bound_dim * ndofn,'DO_NOT_INITIALIZE')
!!$          call memory_alloca(memor_loc,'RECVBUFF_IP','par_interface_node_exchange_rp',recvbuff_ip,bound_dim * ndofn,'DO_NOT_INITIALIZE')
!!$          if( ip == 4 ) then
!!$             PAR_INTEGER4 = MPI_INTEGER
!!$          else
!!$             PAR_INTEGER4 = MPI_INTEGER8
!!$          end if          
!!$       end select
!!$       !
!!$       ! Save in temp_send
!!$       !
!!$       kk = 0
!!$       do jj = 1,comm % bound_dim
!!$          ipoin = comm % bound_perm(jj)
!!$          do idofn = 1,ndofn
!!$             kk = kk + 1
!!$             select type ( xx )
!!$             type is ( real    (kind=rp) )
!!$                sendbuff_rp(kk) = xx(idofn,ipoin)
!!$                recvbuff_rp(kk) = 0.0_rp
!!$             type is ( integer (kind=ip) )
!!$                sendbuff_ip(kk) = xx(idofn,ipoin)
!!$                recvbuff_ip(kk) = 0_ip
!!$             end select
!!$          end do
!!$       end do
!!$       !
!!$       ! Send    temp_send
!!$       ! Receive temp_recv
!!$       !
!!$       istat4 = 0_4
!!$       do ii = 1,comm % nneig
!!$          
!!$          dom_i4 = int(comm % neights(ii),4)
!!$          ini    = ndofn * ( comm % bound_size(ii)   - 1 ) + 1
!!$          nsize  = ndofn * ( comm % bound_size(ii+1) - 1 ) + 1 - ini          
!!$          nsize4 = int(nsize,4)
!!$          
!!$          select type ( xx )
!!$          type is ( real    (kind=rp) )
!!$             call MPI_Sendrecv(                       &
!!$                  sendbuff_rp(ini:), nsize4,          &
!!$                  MPI_DOUBLE_PRECISION, dom_i4, 0_4,  &
!!$                  recvbuff_rp(ini:), nsize4,          &
!!$                  MPI_DOUBLE_PRECISION, dom_i4, 0_4,  &
!!$                  PAR_COMM_WORLD4, status, istat4     )
!!$          type is ( integer (kind=ip) )
!!$             call MPI_Sendrecv(                       &
!!$                  sendbuff_ip(ini:), nsize4,          &
!!$                  PAR_INTEGER4, dom_i4, 0_4,          &
!!$                  recvbuff_ip(ini:), nsize4,          &
!!$                  PAR_INTEGER4,dom_i4, 0_4,           &
!!$                  PAR_COMM_WORLD4, status, istat4     )
!!$          end select
!!$          
!!$          if( istat4 /= MPI_SUCCESS )call PAR_MPI_RUNEND(istat4,'PAR_INTERFACE_NODE_EXCHANGE_RP')
!!$
!!$       end do
!!$
!!$       kk = 0
!!$       do jj = 1,comm % bound_dim
!!$          ipoin = comm % bound_perm(jj)
!!$          do idofn = 1,ndofn
!!$             kk = kk + 1
!!$             select type ( xx )
!!$             type is ( real    (kind=rp) ) ; xx(idofn,ipoin) = xx(idofn,ipoin) + recvbuff_rp(kk)
!!$             type is ( integer (kind=ip) ) ; xx(idofn,ipoin) = xx(idofn,ipoin) + recvbuff_ip(kk)
!!$             end select
!!$          end do
!!$       end do
!!$
!!$       call memory_deallo(memor_loc,'SENDBUFF_RP','par_interface_node_exchange_rp',sendbuff_rp)
!!$       call memory_deallo(memor_loc,'RECVBUFF_RP','par_interface_node_exchange_rp',recvbuff_rp)
!!$       call memory_deallo(memor_loc,'SENDBUFF_IP','par_interface_node_exchange_rp',sendbuff_ip)
!!$       call memory_deallo(memor_loc,'RECVBUFF_IP','par_interface_node_exchange_rp',recvbuff_ip)
!!$
!!$    end if
!!$    
!!$#endif
!!$
!!$  end subroutine par_interface_exchange_2

end module mod_comm_basic
!> @}
