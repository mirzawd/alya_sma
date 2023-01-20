!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Coupling
!> @{
!> @file    cou_communication_arrays.f90
!> @author  Guillaume Houzeaux
!> @date    12/05/2014
!> @brief   Compute some communication arrays
!> @details Compute some communication arrays need by the coupling.
!>          It uses the same MPI communicator as the intra-zone MPI 
!>          communicator
!>
!>          When using a projection method, interface nodes have
!>          non-assembled resdiual, as the contrbution come from 
!>          boundary integrals. Therefore, this one should be 
!>          assembled after the parallel exchange:
!>
!>          \verbtatim
!>
!>                        interface
!>                 CPU1     node      CPU2
!>                           /\
!>          o----x------x----o o----x------x----o
!>
!>          \endverbtatim
!>
!> @} 
!-----------------------------------------------------------------------

subroutine cou_communication_arrays()
  use def_kintyp,         only :  ip,rp
  use def_master,         only :  ISLAVE
  use def_master,         only :  current_code
  use def_master,         only :  current_zone
  use def_domain,         only :  npoin,nzone
  use mod_parall,         only :  par_code_zone_subd_to_color
  use mod_parall,         only :  PAR_COMM_COLOR_ARRAY
  use def_coupli,         only :  memor_cou,mcoup
  use def_coupli,         only :  BETWEEN_SUBDOMAINS
  use def_coupli,         only :  RESIDUAL
  use def_coupli,         only :  PROJECTION
  use def_coupli,         only :  STRESS_PROJECTION
  use def_coupli,         only :  coupling_type
  use def_coupli,         only :  COU_COMM_COUPLING_ARRAY
  use def_coupli,         only :  ncoup_implicit
  use mod_couplings,      only :  I_AM_INVOLVED_IN_A_COUPLING_TYPE
  use mod_communications, only :  PAR_INTERFACE_NODE_EXCHANGE
  use mod_memory,         only :  memory_alloca
  use mod_memory,         only :  memory_deallo
  implicit none 
  integer(ip)                  :: ipoin_wet,ipoin,icoup
  integer(ip)                  :: icolo
  integer(ip), pointer         :: gisca(:)
  !
  ! Needed for each zone
  !
  if( ISLAVE .and. ncoup_implicit > 0 ) then 

     nullify(gisca)
     call memory_alloca(memor_cou,'GISCA','par_communication_arrays',gisca,npoin) 
 
     if(  I_AM_INVOLVED_IN_A_COUPLING_TYPE(BETWEEN_SUBDOMAINS,RESIDUAL,PROJECTION) .or. &
        & I_AM_INVOLVED_IN_A_COUPLING_TYPE(BETWEEN_SUBDOMAINS,RESIDUAL,STRESS_PROJECTION)) then
        !
        ! Allocate memory for communicator and initialize it
        !
        allocate( COU_COMM_COUPLING_ARRAY(0:nzone) )
        do icoup = 0,nzone
           call COU_COMM_COUPLING_ARRAY(icoup) % init(COMM_NAME='COU_COMM_COUPLING_ARRAY')
        end do
        !
        ! Mark my wet nodes of projection type between subdomains: GISCA(IPOIN)=-1 
        !      
        do icoup = 1,mcoup
           if(       coupling_type(icoup) % kind  == BETWEEN_SUBDOMAINS .and. &
                &    coupling_type(icoup) % what  == RESIDUAL           .and. &
                &  ( coupling_type(icoup) % itype == PROJECTION         .or.  &
                &    coupling_type(icoup) % itype == STRESS_PROJECTION) ) then
              do ipoin_wet = 1,coupling_type(icoup) % wet % npoin_wet
                 ipoin = coupling_type(icoup) % wet % lpoin_wet(ipoin_wet)
                 gisca(ipoin) = -1
             end do
           end if
        end do
        !
        ! For all zones, compute the sub-communicator just for the nodes with GISCA(IPOIN)=1
        !
        do current_zone = 0,nzone
           do ipoin = 1,npoin
              gisca(ipoin) =  abs(gisca(ipoin))
           end do
           icolo = par_code_zone_subd_to_color(current_code,current_zone,0_ip)
           call COU_COMM_COUPLING_ARRAY(current_zone) % sub(PAR_COMM_COLOR_ARRAY(icolo),gisca)
           do ipoin = 1,npoin
              gisca(ipoin) = -abs(gisca(ipoin))
           end do
        end do
     end if

     call memory_deallo(memor_cou,'GISCA','par_communication_arrays',gisca)       

     current_zone = 0

  end if

end subroutine cou_communication_arrays
