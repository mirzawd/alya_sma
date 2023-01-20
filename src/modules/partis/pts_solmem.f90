!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Partis
!! @{
!> @name    Partis memory
!! @file    pts_solmem.f90
!> @author  Guillaume Houzeaux
!> @date    05/02/2020
!! @brief   This routine allocate memory for particles
!! @details Allocate memory and initialize particle type
!> @}
!------------------------------------------------------------------------

subroutine pts_solmem()
  use def_parame
  use def_domain
  use def_master
  use def_kermod
  use def_partis
  use def_solver
  use def_inpout
  use mod_memory
  use mod_messages,       only : livinf
  implicit none
  integer(ip) :: itype
  !
  ! Number of used types
  !
  ntyla_pts = 0
  do itype = 1,mtyla
     if( parttyp(itype) % kfl_exist == 1 ) ntyla_pts = max(ntyla_pts,itype)
  end do

  if( INOTMASTER ) then
     !
     ! Velocity deformation tensor
     !
     itype_loop: do itype = 1,ntyla_pts
        if( parttyp(itype) % kfl_exist /= 0 .and. parttyp(itype) % kfl_saffm /= 0 ) then
           call memory_alloca(mem_modul(1:2,modul),'DEFOR_PTS','pts_memall',defor_pts,ntens,npoin)
           exit itype_loop
        end if
     end do itype_loop
     !
     ! Wall element
     !
     call memory_alloca(mem_modul(1:2,modul),'LBOUE_PTS','pts_memall',lboue_pts,nelem)
     !
     ! Element natural length
     !
     call memory_alloca(mem_modul(1:2,modul),'HLENG_PTS','pts_memall',hleng_pts,nelem)
     !
     ! Distance to slip boundaries and friction coefficient
     !
     if( kfl_slip_wall_pts > 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'WALLD_SLIP_PTS','pts_memall',walld_slip_pts,npoin)
        call memory_alloca(mem_modul(1:2,modul),'FRICTION_PTS'  ,'pts_memall',friction_pts,npoin)
     end if
     !
     ! Distance to bouncing boundaries and friction coefficient
     !
     if( kfl_bouncing_wall_pts > 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'WALLD_BOUNCING_PTS','pts_memall',walld_bouncing_pts,npoin)
     end if

  else

     if( kfl_slip_wall_pts > 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'WALLD_SLIP_PTS','pts_memall',walld_slip_pts,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'FRICTION_PTS'  ,'pts_memall',friction_pts,1_ip)
     end if
     if( kfl_bouncing_wall_pts > 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'WALLD_BOUNCING_PTS','pts_memall',walld_bouncing_pts,1_ip)
     end if

  end if
  !
  ! Solver
  !
  solve_sol                => solve(1:1)
  solve_sol(1) % kfl_fixno => kfl_fixno_walld_slip_pts
  call soldef(4_ip)

  solve_sol                => solve(2:2)
  solve_sol(1) % kfl_fixno => kfl_fixno_walld_bouncing_pts
  call soldef(4_ip)
  
end subroutine pts_solmem
