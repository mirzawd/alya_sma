!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    pts_inibcs.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966   
!> @brief   Impose boundary conditions on Partis
!> @details Impose boundary conditions on Partis
!> @} 
!-----------------------------------------------------------------------

subroutine pts_inibcs()

  use def_parame
  use def_inpout
  use def_master
  use def_kermod
  use def_domain
  use def_partis
  use mod_memory
  use mod_opebcs
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications, only : PAR_SUM
  implicit none
  integer(ip) :: iboun,inodb,ipoin

  if( INOTMASTER ) then

     !-------------------------------------------------------------------
     !
     ! Allocate memory
     ! Boundary conditions are of size NBOUN_2 in order to treat wall 
     ! impacts for which the halo's boundary conditions are required
     !
     !-------------------------------------------------------------------

     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXBO_PTS','pts_inibcs',kfl_fixbo_pts,nboun)
     call memory_alloca(mem_modul(1:2,modul),'BVNAT_PTS'    ,'pts_inibcs',bvnat_pts,1_ip,nboun)

     !-------------------------------------------------------------------
     !
     ! Boundary codes
     !
     !-------------------------------------------------------------------

     if( kfl_icodb > 0 ) then       
        kfl_fixbo => kfl_fixbo_pts
        bvnat     => bvnat_pts
        tbcod     => tbcod_pts
        call reacod(IMPOSE_BOUNDARY_CODES)
     end if

  end if

  !------------------------------------------------------------------
  !
  ! Slip wall distance
  !
  !-------------------------------------------------------------------

  kfl_slip_wall_pts = 0
  if( nboun > 0 ) then
     kfl_slip_wall_pts = count(kfl_fixbo_pts == PTS_SLIP_CONDITION) 
  end if
  call PAR_SUM(kfl_slip_wall_pts)

  if( kfl_slip_wall_pts > 0 ) then
     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXNO_WALLD_SLIP_PTS','pts_memall',kfl_fixno_walld_slip_pts,1_ip,npoin)
     do iboun = 1,nboun
        if( kfl_fixbo_pts(iboun) == PTS_SLIP_CONDITION ) then
           do inodb = 1,lnnob(iboun)
              ipoin = lnodb(inodb,iboun)
              kfl_fixno_walld_slip_pts(1,ipoin) = 1
           end do
        end if
     end do
     call PAR_INTERFACE_NODE_EXCHANGE(kfl_fixno_walld_slip_pts,'MAX','IN MY CODE')
  end if

  !------------------------------------------------------------------
  !
  ! Bouncing wall distance
  !
  !-------------------------------------------------------------------

  kfl_bouncing_wall_pts = 0
  if( nboun > 0 ) then
     kfl_bouncing_wall_pts = count(kfl_fixbo_pts == PTS_BOUNCING_CONDITION) 
  end if
  call PAR_SUM(kfl_bouncing_wall_pts)
  if( kfl_bouncing_wall_pts > 0 ) then
     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXNO_WALLD_BOUNCING_PTS','pts_memall',kfl_fixno_walld_bouncing_pts,1_ip,npoin)
     do iboun = 1,nboun
        if( kfl_fixbo_pts(iboun) == PTS_BOUNCING_CONDITION ) then
           do inodb = 1,lnnob(iboun)
              ipoin = lnodb(inodb,iboun)
              kfl_fixno_walld_bouncing_pts(1,ipoin) = 1
           end do
        end if
     end do
     call PAR_INTERFACE_NODE_EXCHANGE(kfl_fixno_walld_bouncing_pts,'MAX','IN MY CODE')
  end if

end subroutine pts_inibcs
