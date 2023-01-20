!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_doiter.f90
!> @author  Solidz Team
!> @date
!> @brief   This routine controls the internal loop of the module.
!> @details
!>          \verbatim
!>          Iterations depend on the integration scheme:
!>          Implicit scheme
!>            - Iterates until convergence is achieved
!>          Explicit scheme
!>            - Perform only one iteration (no convergence)
!>          Runge-Kutta scheme
!>            - Perform only one iteration (no convergence)
!>          \endverbatim
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_doiter

  use def_kintyp,  only : ip
  use def_master,  only : ITASK_ENDINN, ITASK_ENDITE, ITASK_BEGITE
  use def_domain,  only : kfl_elcoh
  use def_solidz,  only : kfl_rigid_sld
  use def_solidz,  only : kfl_stead_sld, kfl_goite_sld
  use mod_sld_rbo, only : sld_rbo_doiter
  use mod_timings, only : timings_ini
  use mod_timings, only : timings_end
#if defined COMMDOM && COMMDOM==2
  use def_solidz,                only : kfl_conta_sld
  use mod_sld_pdn_contact,       only : commdom_sld_nodes_release
  use mod_sld_pdn_contact_plepp, only : SLD_PDN_RIGID_TO_DEFOR
  use mod_sld_pdn_contact_plepp, only : sld_pdn_contact_releas
#endif

  implicit none

  if( kfl_stead_sld == 0_ip ) then

     if ( kfl_rigid_sld == 0_ip ) then

        !-------------------------------------------------------------------
        !
        ! Deformable body
        !
        !-------------------------------------------------------------------
        !
        ! Begin iteration
        !
        call timings_ini()
        call sld_begite()
        call timings_end(ITASK_BEGITE)
        !
        ! Iterate
        !
        do while( kfl_goite_sld == 1_ip ) ! if not converged
           call sld_solite()
           call sld_endite(ITASK_ENDINN)
#if defined COMMDOM && COMMDOM==2
           if( kfl_conta_sld == SLD_PDN_RIGID_TO_DEFOR ) then
              call sld_pdn_contact_releas()
           else
              call commdom_sld_nodes_release()
           end if
#endif
        end do
        !
        ! End iteration
        !
        call timings_ini()
        call sld_endite(ITASK_ENDITE)
        call timings_end(ITASK_ENDITE)
        !
        ! Update traction cohesive elements (XFM/CZM)
        !
        if( kfl_elcoh /= 0 ) call sld_uptcoh()

     else if ( kfl_rigid_sld == 1_ip ) then

        !-------------------------------------------------------------------
        !
        ! Rigid body
        !
        !-------------------------------------------------------------------

        call sld_rbo_doiter()

     end if

  end if

end subroutine sld_doiter
