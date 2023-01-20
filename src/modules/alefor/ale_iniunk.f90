!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Alefor
!> @{
!> @file    ale_iniunk.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966
!> @brief   This routine sets up the initial condition for the mesh velocity
!> @details This routine sets up the initial condition for the mesh velocity
!> @} 
!-----------------------------------------------------------------------
subroutine ale_iniunk()
  use def_master
  use def_domain
  use def_alefor
  use def_kermod, only : kfl_adj_prob,sens_mesh

  implicit none

  integer(ip) :: idime,ipoin
  !
  ! Reading mesh sensitivities for optimization
  !
  if (kfl_adj_prob == 1_ip) then
     do ipoin=1,npoin
        do idime = 1, ndime
           sens_mesh(idime,ipoin) = xfiel(-kfl_sensi_ale) %a(idime,ipoin,1)
        end do
     end do
  end if

  if( kfl_rstar == 0 ) then
     
     !------------------------------------------------------------------
     !
     ! Initial solution
     !
     !------------------------------------------------------------------
     !
     ! Coordinates
     !
     do ipoin = 1,npoin
        coord_ale(:,ipoin,1) = coord(:,ipoin)
        coord_ale(:,ipoin,2) = coord(:,ipoin)
        coord_ale(:,ipoin,3) = coord(:,ipoin)
        coord_ori(:,ipoin)   = coord(:,ipoin)
     end do
     
     if (moddi_ale < 0) then
        do ipoin = 1,npoin
           coord(:,ipoin) = coord(:,ipoin) + xfiel(-moddi_ale)%a(:,ipoin,1)
        end do
     end if
     !
     ! Update boundary conditions
     !
     call ale_updbcs()
     !
     ! Initial solution
     !
     call ale_updunk(ITASK_INIUNK) 
     !
     ! Body fitter
     !
     if( kfl_rigid_ale == 1 ) call ale_inirbo()

     !------------------------------------------------------------------
     !
     ! If ALE is only used for smoothing purpose
     !
     !------------------------------------------------------------------

     if( coupling('ALEFOR','SOLIDZ') == 0 ) then
        if( kfl_delay(modul) == 0 .and.  ndela(modul) >= 0 ) then
           call ale_doiter()
        end if
     end if

  end if

  call mescek(1_ip)

end subroutine ale_iniunk
