!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    elsini.f90
!> @author  houzeaux
!> @date    2019-05-03
!> @brief   Initialize Elsest
!> @details Allocate Elsest data structure and compute preprocess
!>          for one single mesh
!> @} 
!-----------------------------------------------------------------------

subroutine elsini()

  use def_kintyp,   only : ip
  use def_master,   only : INOTMASTER
  use def_domain,   only : meshe
  use def_kermod,   only : ndivi
  use def_kermod,   only : ielse
  use def_kermod,   only : relse
  use def_kermod,   only : kfl_elses
  use def_kermod,   only : ndivi
  use mod_elsest,   only : elsest_preprocess
  use mod_elsest,   only : elsest_allocate
  use mod_messages, only : messages_live
  
  implicit none

  if( kfl_elses == 1 ) then
     call messages_live('ELSEST PREPROCESS')
     if( INOTMASTER ) then
        call elsest_allocate  (ielse,NUMBER_MESHES=1_ip)
        call elsest_preprocess(ielse,relse,meshe(ndivi),CURRENT_MESH=1_ip)
     end if
  end if

end subroutine elsini
