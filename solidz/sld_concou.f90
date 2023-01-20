!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_concou.f90
!> @author  Guillaume Houzeaux
!> @date    September 2006
!> @brief   Check convergence outer iterations for block coupling
!> @details Check convergence outer iterations for block coupling
!>          
!> @} 
!-----------------------------------------------------------------------

subroutine sld_concou()

  use def_kintyp, only : ip
  use def_master, only : kfl_conve, kfl_gocou
  use def_master, only : glres,coutp,routp,lun_outpu
  use def_master, only : modul
  use mod_outfor, only : outfor
  use def_solidz, only : resid_sld, cotol_sld
  
  implicit none
  !
  ! Check convergence
  !
  if( kfl_conve(modul) == 1 ) then
     if( resid_sld > cotol_sld ) kfl_gocou = 1
  end if
  glres(modul) = resid_sld
  !
  ! Output residuals
  !
  coutp(1) = 'DISPLACEMENT'
  routp(1) = resid_sld
  call outfor(9_ip,lun_outpu,' ')
    
end subroutine sld_concou
