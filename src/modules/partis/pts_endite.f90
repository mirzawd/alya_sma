!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Partis
!! @{
!> @name    Partis inner iteration
!! @file    pts_doiter.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!! @brief   This routine ends a time step for particles
!! @details Update residence time
!> @} 
!------------------------------------------------------------------------

subroutine pts_endite()
  use def_kintyp
  use def_master
  use def_domain
  use def_partis
  use mod_communications, only : PAR_FROM_GHOST_ELEMENT_EXCHANGE
  use mod_communications, only : PAR_FROM_GHOST_NODE_EXCHANGE
  use mod_communications, only : PAR_GHOST_NODE_EXCHANGE
  use mod_communications, only : PAR_FROM_GHOST_BOUNDARY_EXCHANGE
  use mod_pts_transport,  only : pts_transport_finalize
  implicit none
  integer(ip) :: ielem,iboun,ipoin
  !
  ! Compute some numbers and output some information
  !
  call pts_transport_finalize()
  !
  ! Residence: pass information from ghost element
  ! Pass information computed on my ghost element and send it to my neighbors
  !
  if( ISLAVE .and. kfl_resid_pts /= 0 ) then
     call PAR_FROM_GHOST_ELEMENT_EXCHANGE(resid_pts,'SUM','IN MY CODE')
     do ielem = nelem+1,nelem_2
        resid_pts(1:ntyla_pts,ielem) = 0.0_rp
     end do
  end if
  !
  ! Boundary deposition on boundaries
  !
  if( ISLAVE .and. kfl_depos_pts /= 0 ) then
     call PAR_FROM_GHOST_BOUNDARY_EXCHANGE(depob_pts,'SUM','IN MY CODE')
     do iboun = nboun+1,nboun_2
        depob_pts(1:ntyla_pts,iboun) = 0.0_rp
     end do
  end if

  !----------------------------------------------------------------------
  !
  ! Exchange source terms
  !
  !----------------------------------------------------------------------

  if( ISLAVE ) then
     
     if( kfl_momentum_sink_pts /= 0 ) then
        call PAR_FROM_GHOST_NODE_EXCHANGE(momentum_sink,'SUM','IN MY CODE')
        do ipoin = npoin+1,npoin_2
           momentum_sink(:,ipoin) = 0.0_rp
        end do
        call rhsmod(ndime,momentum_sink)
     end if 

     if( kfl_mass_sink_pts /= 0 ) then
        call PAR_FROM_GHOST_NODE_EXCHANGE(mass_sink,'SUM','IN MY CODE')
        do ipoin = npoin+1,npoin_2
           mass_sink(ipoin) = 0.0_rp
        end do
        call rhsmod(1_ip,mass_sink) 
     end if

     if( kfl_heat_sink_pts /= 0 ) then
        call PAR_FROM_GHOST_NODE_EXCHANGE(heat_sink,'SUM','IN MY CODE')
        do ipoin = npoin+1,npoin_2
           heat_sink(ipoin) = 0.0_rp
        end do
        call rhsmod(1_ip,heat_sink) 
     end if
  end if

end subroutine pts_endite
