!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    mod_par_global_numbering.f90
!> @author  houzeaux
!> @date    2018-06-06
!> @brief   Global numbering
!> @details Global numbering stuffs
!>
!-----------------------------------------------------------------------

module mod_par_global_numbering

  use def_kintyp,         only : ip
  use def_kintyp_comm,    only : comm_data_par
  use def_master,         only : INOTMASTER,kfl_paral,IMASTER
  use def_master,         only : ISEQUEN
  use def_master,         only : npoi1,npoi2,npoi3
  use def_master,         only : npart
  use mod_parall,         only : par_memor
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_ALLGATHER
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications, only : PAR_COMM_RANK_AND_SIZE
  use mod_communications, only : PAR_ALLGATHER
  use mod_communications, only : PAR_MAX
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_memory,         only : memory_copy
  use mod_memory,         only : memory_size
  use mod_memory,         only : memory_resize
  use def_mpi
#include "def_mpi.inc"
  implicit none
  private

  public :: par_global_numbering_nodes
  public :: par_global_numbering_elements_boundaries
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-06-06
  !> @brief   ???
  !> @details Find out a global numbering for nodes:
  !>
  !>          LNINV: Uniquely renumber the nodes
  !>
  !>          If these array are given it is assumed that they come
  !>          from a new mesh where old nodes are numbered right after
  !>          the old ones
  !> 
  !-----------------------------------------------------------------------

  subroutine par_global_numbering_nodes(npoin,lninv_loc,COMM)

    integer(ip),                            intent(in)    :: npoin
    integer(ip),                   pointer, intent(inout) :: lninv_loc(:)
    type(comm_data_par), optional,          intent(in)    :: COMM
    integer(ip)                                           :: kpoin,ipoin,ipart
    integer(ip)                                           :: npoin_total_old
    integer(ip),                   pointer                :: own_nodes(:)
    integer(ip)                                           :: npoi1_loc
    integer(ip)                                           :: npoi2_loc
    integer(ip)                                           :: npoi3_loc
    MY_MPI_COMM                                           :: PAR_COMM4
    
    nullify(own_nodes)
    if( present(COMM) ) then
       npoi1_loc = COMM % npoi1
       npoi2_loc = COMM % npoi2
       npoi3_loc = COMM % npoi3
       PAR_COMM4 = COMM % PAR_COMM_WORLD
    else
       npoi1_loc = npoi1
       npoi2_loc = npoi2
       npoi3_loc = npoi3
    end if
    !
    ! Resize array if needed
    !
    if( memory_size(lninv_loc) < npoin ) &
         call memory_resize(par_memor,'LNINV_LOC','mesh_multiplication_divide_elements',lninv_loc,npoin)

    if( INOTMASTER ) npoin_total_old = maxval( lninv_loc )
    call PAR_MAX(npoin_total_old)

    if( INOTMASTER ) then

       do ipoin = npoi1_loc+1,npoi2_loc-1
          lninv_loc(ipoin) = 0
       end do
       do ipoin = npoi3_loc+1,npoin 
          lninv_loc(ipoin) = 0
       end do
       kpoin = 0
       do ipoin = 1,npoin
          if( lninv_loc(ipoin) == 0 ) then
             if( ipoin <= npoi1_loc .or. ( ipoin >= npoi2_loc .and. ipoin <= npoi3_loc ) ) then
                kpoin = kpoin + 1
             end if
          end if
       end do
    else
       kpoin = 0
    end if 

    call memory_alloca(par_memor,'ONW_NODES','par_submsh',own_nodes,npart+1_ip,lboun=0_ip)
    if( present(COMM) ) then
       call PAR_ALLGATHER(kpoin,own_nodes,1_4,PAR_COMM_IN=PAR_COMM4)
    else
       call PAR_ALLGATHER(kpoin,own_nodes,1_4,'IN MY CODE')
    end if
    
    if( INOTMASTER ) then
       do ipart = 1,npart
          own_nodes(ipart) = own_nodes(ipart-1) + own_nodes(ipart)
       end do
       kpoin = own_nodes(kfl_paral-1) + npoin_total_old 

       do ipoin = 1,npoin
          if( lninv_loc(ipoin) == 0 ) then
             if( ipoin <= npoi1_loc .or. ( ipoin >= npoi2_loc .and. ipoin <= npoi3_loc ) ) then
                kpoin = kpoin + 1
                lninv_loc(ipoin) = kpoin
             end if
          end if
       end do
 
       if( present(COMM) ) then
          call PAR_INTERFACE_NODE_EXCHANGE(lninv_loc,'SUM',COMM=COMM)
       else
          call PAR_INTERFACE_NODE_EXCHANGE(lninv_loc,'SUM')
       end if

    end if

    call memory_deallo(par_memor,'ONW_NODES','par_submsh',own_nodes)

  end subroutine par_global_numbering_nodes

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-15
  !> @brief   Renumber elements or boundaries
  !> @details Specific treatment for numbering elements and boundaries
  !>          as they are disjoint sets
  !> 
  !-----------------------------------------------------------------------

  subroutine par_global_numbering_elements_boundaries(nxxxx,lxinv_loc)

    integer(ip),          intent(in)    :: nxxxx        !< Current number of entitites
    integer(ip), pointer, intent(inout) :: lxinv_loc(:) !< Current or old glboal numbering

    integer(ip)                         :: ipart,kxxxx,ixxxx
    integer(ip), pointer                :: nxxxx_tot(:)
    integer(ip), pointer                :: offset(:)

    nullify(nxxxx_tot)
    nullify(offset)
    !
    ! Resize array if needed
    !
    if( memory_size(lxinv_loc) < nxxxx ) &
         call memory_resize(par_memor,'LXINV_LOC','mesh_multiplication_divide_elements',lxinv_loc,nxxxx)

    call memory_alloca(par_memor,'NXXXX_TOT','par_submsh',nxxxx_tot,npart+1_ip,lboun=0_ip)
    call PAR_ALLGATHER(nxxxx,nxxxx_tot)

    if( ISEQUEN    ) then

       do ixxxx = 1,nxxxx
          lxinv_loc(ixxxx) = ixxxx
       end do
       
    else if( INOTMASTER ) then
 
       call memory_alloca(par_memor,'OFFSET','par_submsh',offset,npart+1_ip,lboun=0_ip)
       do ipart = 0,npart-1
          offset(ipart+1) = nxxxx_tot(ipart)
       end do
       do ipart = 1,npart
          offset(ipart) = offset(ipart-1) + offset(ipart) 
       end do
      
       kxxxx = offset(kfl_paral)
       do ixxxx = 1,nxxxx
          kxxxx = kxxxx + 1
          lxinv_loc(ixxxx) = kxxxx
       end do

    end if

    call memory_deallo(par_memor,'NXXXX_TOT','par_submsh',nxxxx_tot)
    call memory_deallo(par_memor,'OFFSET'   ,'par_submsh',offset)

  end subroutine par_global_numbering_elements_boundaries
  
end module mod_par_global_numbering
!> @}
 
