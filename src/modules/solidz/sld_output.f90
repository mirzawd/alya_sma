!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_output.f90
!> @author  Solidz Team
!> @date    August, 2006
!>          - Subroutine written
!> @brief   Solidz Output centralization for post-process of results
!>
!> @details
!>
!>          \verbatim
!>          Field Output:
!>          -------------
!>          Variables defined for post-processs visualization.
!>
!>          History Output:
!>          ---------------
!>          Variables defined on sets or witness points.
!>
!>          \endverbatim
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_output()

  use def_kintyp,             only : ip, rp
  use def_master,             only : ITASK_INITIA, ITASK_ENDTIM, ITASK_ENDRUN, INOTSLAVE
  use def_master,             only : ittyp
  use mod_output_postprocess, only : output_postprocess_variables
  use def_solidz,             only : kfl_foten_sld
  use mod_sld_post_reaction,  only : kfl_preac_sld, sld_post_reaction
  use mod_eccoupling,         only : kfl_exmsld_ecc, eccou_save_log

  implicit none

  external :: sld_outvar
  external :: sld_outset
  external :: sld_outwit
  external :: sld_exaerr
  external :: sld_outinf
  !
  ! Initialize flag fote (only after each time step)
  !
  if( ittyp == ITASK_INITIA ) kfl_foten_sld = 1_ip
  if( ittyp == ITASK_ENDTIM ) kfl_foten_sld = 0_ip
  if( ittyp == ITASK_ENDRUN ) kfl_foten_sld = 1_ip

  if( INOTSLAVE ) then
    if( ittyp == ITASK_INITIA ) then
        if( kfl_exmsld_ecc ) then
            call eccou_save_log()
        end if   
    end if
  end if

  !
  ! Post-process Field Output variables
  !
  call output_postprocess_variables(sld_outvar)
  !
  ! Post-process History Output variables
  !
  if( ittyp == ITASK_INITIA .or. ittyp == ITASK_ENDTIM ) then
     !
     ! Post-process on sets
     !
     call sld_outset()
     !
     ! Post-process on witness points
     !
     call sld_outwit()
     !
     ! FEM errors (End of the run)
     !
     call sld_exaerr(2_ip)
     !
     ! User-postprocess for reactions
     !
     if( kfl_preac_sld ) call sld_post_reaction()
        
  end if
  !
  ! Write information to files
  !
  if( ittyp == ITASK_ENDTIM ) then
      call sld_outinf(2_ip)
  end if
  !
  ! Write at the end of the analysis
  !
  if( ittyp == ITASK_ENDRUN ) then
      continue
  end if

end subroutine sld_output
