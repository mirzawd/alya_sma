!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_domgra.f90
!> @author  Guillaume Houzeaux
!> @brief   Zone-wise communication array
!> @details Zone-wise communication array
!> @} 
!------------------------------------------------------------------------
subroutine par_zone_communication_arrays()
  use def_kintyp_basic,   only :  ip,rp
  use def_kintyp_comm,    only :  comm_data_par
  use def_master,         only :  ISEQUEN,ISLAVE
  use def_master,         only :  current_code
  use def_master,         only :  current_zone
  use def_master,         only :  current_subd
  use def_master,         only :  gisca,npoi3
  use def_master,         only :  THIS_NODE_IS_MINE
  use def_domain,         only :  npoin,nzone
  use def_domain,         only :  nelem,lesub,nsubd,lnods
  use def_domain,         only :  lnnod,lperi,nperi
  use mod_communications, only :  PAR_DEFINE_COMMUNICATOR
  use mod_communications, only :  PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications, only :  PAR_SEND_RECEIVE
  use mod_communications, only :  PAR_MAX
  use mod_communications, only :  PAR_COMM_RANK_AND_SIZE
  use mod_parall,         only :  PAR_COMM_COLOR
  use mod_parall,         only :  PAR_COMM_COLOR_ARRAY
  use mod_parall,         only :  PAR_COMM_MY_CODE_ARRAY 
  use mod_parall,         only :  PAR_COMM_COLOR_PERM
  use mod_parall,         only :  par_world_rank_of_a_code_neighbor
  use mod_parall,         only :  PAR_MY_WORLD_RANK
  use mod_parall,         only :  par_code_zone_subd_to_color
  use mod_messages,       only :  messages_live
  use mod_communication_arrays
  use def_mpi
#include "def_mpi.inc"
  implicit none
  MY_MPI_COMM                  :: PAR_COMM_ORIGINAL
  integer(4)                   :: PAR_RANK_ZONE
  integer(ip)                  :: ierro,dom_i_world,dom_i_zone,ielem
  integer(ip)                  :: icolo,ipoin,ineig,dom_i_code
  integer(ip)                  :: inode,jpoin,iperi
  type(comm_data_par), pointer :: commu

  if( ISEQUEN ) return

  call messages_live('PARALL: CREATE INTRA-ZONE COMMUNICATION ARRAY')
  ! 
  ! Allocate memory
  !
  call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_ORIGINAL,commu)
  
  !----------------------------------------------------------------------
  !
  ! All-zone and all-subdomain communication array
  !
  !----------------------------------------------------------------------

  icolo = par_code_zone_subd_to_color(current_code,0_ip,0_ip)
  call PAR_COMM_COLOR_ARRAY(icolo) % copy(PAR_COMM_MY_CODE_ARRAY(1))

  !----------------------------------------------------------------------
  !
  ! Intra-zone communication arrays
  !
  !----------------------------------------------------------------------

  if( ISLAVE ) call memgen(1_ip,npoin,0_ip)

  do current_zone = 1,nzone

     icolo = par_code_zone_subd_to_color(current_code,current_zone,0_ip)
     PAR_RANK_ZONE = PAR_COMM_COLOR_PERM(icolo,icolo,PAR_MY_WORLD_RANK) !TODO: possible int8 to int4 conversion
     !
     ! Sub-communicator only in current zone
     !            
     if( ISLAVE ) then

        do ipoin = 1,npoin
           gisca(ipoin) = 1
        end do
        call PAR_COMM_COLOR_ARRAY(icolo) % sub (commu,gisca,COMM_NAME='PAR_COMM_COLOR_ARRAY')
        !
        ! Convert code neighbor to zone neighbor (dm_i_code is kfl_paral!!!)
        ! This is because the communicator was already defined using a split
        ! at the beginning (PAR_COMM_COLOR_ARRAY(icolo) % PAR_COMM_WORLD)
        ! and commu invovles the code partition numbering
        !
        do ineig = 1,PAR_COMM_COLOR_ARRAY(icolo) % nneig
           dom_i_code  = PAR_COMM_COLOR_ARRAY(icolo) % neights(ineig)
           dom_i_world = par_world_rank_of_a_code_neighbor(dom_i_code,current_code)
           dom_i_zone  = PAR_COMM_COLOR_PERM(icolo,icolo,dom_i_world)
           PAR_COMM_COLOR_ARRAY(icolo) % neights(ineig) = dom_i_zone
        end do
        do ipoin = 1,npoin
           gisca(ipoin) = 0
        end do
     end if
     !
     ! MPI communicator
     !
     PAR_COMM_COLOR_ARRAY(icolo) % PAR_COMM_WORLD = PAR_COMM_COLOR(icolo,icolo)

  end do

  if( ISLAVE ) call memgen(3_ip,npoin,0_ip)

  !----------------------------------------------------------------------
  !
  ! A short test to check inter-zone communication
  !
  !----------------------------------------------------------------------

  do current_zone = 1,nzone
     ierro = 0
     if( ISLAVE ) then
        call memgen(1_ip,npoin,0_ip)
        do ipoin = 1,npoi3
           gisca(ipoin) = 1
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gisca,'SUM','IN MY ZONE','SYNCHRONOUS')        
        do ipoin = 1,npoin
           if( gisca(ipoin) /= 1 ) then
              ierro = ierro + 1
           end if
        end do
        call memgen(3_ip,npoin,0_ip)
     end if
     call PAR_MAX(ierro,'IN MY CODE')
     !if( ierro /= 0 ) call runend('PROBLEM WITH ZONE COMMUNICATION ARRAY')
  end do

  !----------------------------------------------------------------------
  !
  ! Very important: CURRENT_ZONE is all zones
  !
  !----------------------------------------------------------------------

  current_zone = 0

  !---------------------------------------------------------------------------------------------------!
  !
  !  ==========================
  !  DO THE SAME FOR SUBDOMAINS
  !  ==========================
  !
  !---------------------------------------------------------------------------------------------------!

  call messages_live('PARALL: CREATE INTRA-SUBDOMAIN COMMUNICATION ARRAY')
  !
  ! Allocate memory
  !
  call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_ORIGINAL,commu)

  !----------------------------------------------------------------------
  !
  ! Inter-zone communication arrays
  !
  !----------------------------------------------------------------------

  if( ISLAVE ) call memgen(1_ip,npoin,0_ip)

  do current_subd = 1,nsubd

     icolo = par_code_zone_subd_to_color(current_code,0_ip,current_subd)
     !
     ! Sub-communicator only in current zone
     !            
     if( ISLAVE ) then

        do ielem = 1,nelem
           if( lesub(ielem) == current_subd ) then
              do inode = 1,lnnod(ielem)
                 ipoin = lnods(inode,ielem)
                 gisca(ipoin) = 1
              end do
           end if
        end do
        do iperi = 1,nperi
           ipoin = lperi(1,iperi) 
           jpoin = lperi(2,iperi)
           gisca(ipoin) = gisca(jpoin) 
        end do

        call PAR_INTERFACE_NODE_EXCHANGE(gisca,'MAX','IN MY CODE')
        call PAR_COMM_COLOR_ARRAY(icolo) % sub(commu,gisca,COMM_NAME='PAR_COMM_COLOR_ARRAY')
        !
        ! Convert code neighbor to zone neighbor (dm_i_code is kfl_paral!!!)
        ! This is because the communicator was already defined using a split
        ! at the beginning (PAR_COMM_COLOR_ARRAY(icolo) % PAR_COMM_WORLD)
        ! and commu invovles the code partition numbering
        !
        do ineig = 1,PAR_COMM_COLOR_ARRAY(icolo) % nneig
           dom_i_code  = PAR_COMM_COLOR_ARRAY(icolo) % neights(ineig)
           dom_i_world = par_world_rank_of_a_code_neighbor(dom_i_code,current_code)
           dom_i_zone  = PAR_COMM_COLOR_PERM(icolo,icolo,dom_i_world)
           PAR_COMM_COLOR_ARRAY(icolo) % neights(ineig) = dom_i_zone
        end do
        do ipoin = 1,npoin
           gisca(ipoin) = 0
        end do
     end if
     !
     ! MPI communicator
     !
     PAR_COMM_COLOR_ARRAY(icolo) % PAR_COMM_WORLD = PAR_COMM_COLOR(icolo,icolo)

  end do

  if( ISLAVE ) call memgen(3_ip,npoin,0_ip)
  
  !----------------------------------------------------------------------
  !
  ! A short test to check intra-subdomain communication
  !
  !----------------------------------------------------------------------
  
  do current_subd = 1,nsubd
     ierro = 0
     if( ISLAVE ) then
        call memgen(1_ip,npoin,0_ip)
 
        do ielem = 1,nelem
           if( lesub(ielem) == current_subd ) then
              do inode = 1,lnnod(ielem)
                 ipoin = lnods(inode,ielem)
                 if( THIS_NODE_IS_MINE(ipoin) ) then
                    gisca(ipoin) =  1
                 else
                    gisca(ipoin) = -1
                 end if
              end do
           end if
        end do

        do iperi = 1,nperi
           ipoin = lperi(1,iperi)
           jpoin = lperi(2,iperi)
           if( THIS_NODE_IS_MINE(ipoin) .and. abs(gisca(jpoin)) == 1 ) gisca(ipoin) = 1
        end do

        if( npoin > 0 ) gisca = max(0_ip,gisca)
        
        call PAR_INTERFACE_NODE_EXCHANGE(gisca,'SUM','IN MY SUBD','SYNCHRONOUS')        
        
        do ielem = 1,nelem
           if( lesub(ielem) == current_subd ) then
              do inode = 1,lnnod(ielem)
                 ipoin = lnods(inode,ielem)
                 if( gisca(ipoin) /= 1 ) then
                    ierro = ierro + 1
                 end if
              end do
           end if
        end do
        call memgen(3_ip,npoin,0_ip)
     end if
     !call PAR_MAX(ierro,'IN MY CODE')
     !if( ierro /= 0 ) call runend('PROBLEM WITH SUBDOMAIN COMMUNICATION ARRAY')
  end do

  !----------------------------------------------------------------------
  !
  ! Very important: CURRENT_SUBD is all subdomains
  !
  !----------------------------------------------------------------------

  current_subd = 0

end subroutine par_zone_communication_arrays
